!=======================================================================
!   SEIS_FILO: 
!   SEISmological tools for Flat Isotropic Layered structure in the Ocean
!   Copyright (C) 2019 Takeshi Akuhara
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program.  If not, see <https://www.gnu.org/licenses/>.
!
!
!   Contact information
!
!   Email  : akuhara @ eri. u-tokyo. ac. jp 
!   Address: Earthquake Research Institute, The Univesity of Tokyo
!           1-1-1, Yayoi, Bunkyo-ku, Tokyo 113-0032, Japan
!
!=======================================================================
module cls_mcmc
  use mod_random
  use cls_trans_d_model
  use cls_proposal
  use cls_hyper_model
  implicit none 
  
  type mcmc 
     private
     
     type(trans_d_model) :: tm
     type(hyper_model) :: hyp_rf, hyp_disp
     type(proposal) :: prop
     
     integer :: n_iter
     integer :: n_mod
     integer :: i_mod
     integer :: i_iter
     integer :: i_proposal_type
     integer :: n_proposal_type
     integer, allocatable :: n_propose(:)
     integer, allocatable :: n_accept(:)
     integer, allocatable :: k_saved(:)
     double precision :: log_likelihood
     double precision, allocatable :: likelihood_saved(:)
     double precision, allocatable :: temp_saved(:)
     double precision :: temp = 1.d0
     logical :: is_accepted
   contains
     procedure :: propose_model => mcmc_propose_model
     procedure :: judge_model   => mcmc_judge_model
     procedure :: one_step_summary => mcmc_one_step_summary
     procedure :: get_tm  => mcmc_get_tm
     procedure :: get_hyp_disp => mcmc_get_hyp_disp
     procedure :: get_hyp_rf => mcmc_get_hyp_rf
     procedure :: get_is_accepted => mcmc_get_is_accepted
     procedure :: get_n_iter => mcmc_get_n_iter
     procedure :: set_temp => mcmc_set_temp
     procedure :: get_temp => mcmc_get_temp
     procedure :: get_log_likelihood => mcmc_get_log_likelihood
     procedure :: get_likelihood_saved => mcmc_get_likelihood_saved
     procedure :: get_temp_saved => mcmc_get_temp_saved
     procedure :: get_n_propose => mcmc_get_n_propose
     procedure :: get_n_accept => mcmc_get_n_accept
     procedure :: finish => mcmc_finish
  end type mcmc
  
  interface mcmc
     module procedure init_mcmc
  end interface mcmc


contains

  !---------------------------------------------------------------------
  
  type(mcmc) function init_mcmc(tm, hyp_disp, hyp_rf, &
       & prop, n_iter) result(self)
    type(trans_d_model), intent(in) :: tm
    type(hyper_model), intent(in) :: hyp_rf, hyp_disp
    type(proposal), intent(in) :: prop
    integer, intent(in) :: n_iter
    
    
    self%tm = tm
    self%hyp_disp = hyp_disp
    self%hyp_rf = hyp_rf
    self%prop = prop
    self%n_iter = n_iter
    self%i_iter = 0
    self%i_mod = 0
    self%log_likelihood = -1.0d300
    self%n_proposal_type = self%prop%get_n_proposal()
    allocate(self%likelihood_saved(n_iter))
    allocate(self%temp_saved(n_iter))
    allocate(self%k_saved(n_iter))
    allocate(self%n_propose(self%prop%get_n_proposal()))
    allocate(self%n_accept(prop%get_n_proposal()))
    self%n_propose = 0
    self%n_accept = 0

    return 
  end function init_mcmc

  !---------------------------------------------------------------------
  
  subroutine mcmc_propose_model(self, tm_proposed, hyp_disp_proposed, &
       & hyp_rf_proposed, is_ok, &
       & log_prior_ratio, log_proposal_ratio)
    class(mcmc), intent(inout) :: self
    type(trans_d_model), intent(out) :: tm_proposed
    type(hyper_model), intent(out) :: hyp_disp_proposed
    type(hyper_model), intent(out) :: hyp_rf_proposed
    double precision, intent(out) :: log_prior_ratio
    double precision, intent(out) :: log_proposal_ratio
    
    logical, intent(out) :: is_ok
    integer :: itype

    itype = int(rand_u() * self%prop%get_n_proposal()) + 1
    self%i_proposal_type = itype
    tm_proposed = self%tm
    hyp_disp_proposed = self%hyp_disp
    hyp_rf_proposed = self%hyp_rf
    if (   itype == self%prop%get_i_vs() .or. &
         & itype == self%prop%get_i_vp() .or. &
         & itype == self%prop%get_i_depth()) then
       ! Perturbation of Trans-D model
       call tm_proposed%perturb(itype, is_ok, log_prior_ratio)
       log_proposal_ratio = 0.d0
    else if (itype == self%prop%get_i_birth()) then
       ! Birth
       call tm_proposed%birth(is_ok)
       ! prior & proposal ratio 
       ! (can be ignored if birth prorposal from prior distribution)
       log_prior_ratio = 0.d0
       log_proposal_ratio = 0.d0
    else if (itype == self%prop%get_i_death()) then
       ! Death
       call tm_proposed%death(is_ok)
       log_prior_ratio = 0.d0
       log_proposal_ratio = 0.0
    else if (itype == self%prop%get_i_rf_sig()) then
       call hyp_rf_proposed%perturb(is_ok)
       log_prior_ratio = 0.d0
       log_proposal_ratio = 0.0
    else if (itype == self%prop%get_i_disper_sig()) then
       call hyp_disp_proposed%perturb(is_ok)
       log_prior_ratio = 0.d0
       log_proposal_ratio = 0.0
    end if
    
    ! Count proposal
    self%n_propose(itype) = &
         & self%n_propose(itype) + 1
    
    return
  end subroutine mcmc_propose_model

  !---------------------------------------------------------------------

  subroutine mcmc_judge_model(self, tm, hyp_disp, hyp_rf, prior_ok, &
       & log_likelihood, log_prior_ratio, log_proposal_ratio)
    class(mcmc), intent(inout) :: self
    type(trans_d_model), intent(in) :: tm
    type(hyper_model), intent(in) :: hyp_disp, hyp_rf
    double precision, intent(in) :: log_likelihood, log_prior_ratio, &
         & log_proposal_ratio
    logical, intent(in) :: prior_ok
    double precision :: ratio
    double precision :: r
    double precision, parameter :: eps = epsilon(1.d0)
    
    self%is_accepted = .false.
    ratio = (log_likelihood - self%log_likelihood) / self%temp
    ratio = ratio + log_prior_ratio + log_proposal_ratio
    
    r = rand_u()
    if (r >= eps) then
       if (log(r) <= ratio) then
          self%is_accepted = .true.
       end if
    end if
       
    if (.not. prior_ok) then
       self%is_accepted = .false.
    end if

    if (self%is_accepted) then
       ! Accept model
       self%tm = tm
       self%hyp_rf = hyp_rf
       self%hyp_disp = hyp_disp
       self%log_likelihood = log_likelihood
       self%n_accept(self%i_proposal_type) = &
         & self%n_accept(self%i_proposal_type) + 1
    end if

    ! Adds iteration counter
    self%i_iter = self%i_iter + 1

    ! Save models etc.
    self%likelihood_saved(self%i_iter) = self%log_likelihood
    self%k_saved(self%i_iter) = self%tm%get_k()
    self%temp_saved(self%i_iter) = self%temp


    return 
  end subroutine mcmc_judge_model
  
  !---------------------------------------------------------------------

  subroutine mcmc_one_step_summary(self)
    class(mcmc), intent(in) :: self
    write(*,*)"Iteration    : ", self%i_iter, "/", self%n_iter
    write(*,*)"Proposal type: ", self%i_proposal_type
    write(*,*)"Accepted     : ", self%is_accepted
    write(*,*)"Likelihood   : ", self%log_likelihood
    write(*,*)
  end subroutine mcmc_one_step_summary

  !---------------------------------------------------------------------
  
  type(trans_d_model) function mcmc_get_tm(self) result(tm)
    class(mcmc), intent(in) :: self
    
    tm = self%tm

    return 
  end function mcmc_get_tm

  !---------------------------------------------------------------------

  type(hyper_model) function mcmc_get_hyp_disp(self) result(hyp_disp)
    class(mcmc), intent(in) :: self
    
    hyp_disp = self%hyp_disp

    return 
  end function mcmc_get_hyp_disp

  !---------------------------------------------------------------------

  type(hyper_model) function mcmc_get_hyp_rf(self) result(hyp_rf)
    class(mcmc), intent(in) :: self
    
    hyp_rf = self%hyp_rf
    
    return 
  end function mcmc_get_hyp_rf
  
  !---------------------------------------------------------------------
  
  logical function mcmc_get_is_accepted(self) result(is_accepted)
    class(mcmc), intent(in) :: self
    
    is_accepted = self%is_accepted

    return 
  end function mcmc_get_is_accepted

  !---------------------------------------------------------------------

  integer function mcmc_get_n_iter(self) result(n_iter)
    class(mcmc), intent(in) :: self
    
    n_iter = self%n_iter

    return 
  end function mcmc_get_n_iter

  !---------------------------------------------------------------------
  
  subroutine mcmc_set_temp(self, temp)
    class(mcmc), intent(inout) :: self
    double precision, intent(in) :: temp

    self%temp = temp

    return 
  end subroutine mcmc_set_temp

  !---------------------------------------------------------------------

  double precision function mcmc_get_temp(self) result(temp)
    class(mcmc), intent(in) :: self

    temp = self%temp
    
    return 
  end function mcmc_get_temp

  !---------------------------------------------------------------------

  double precision function mcmc_get_log_likelihood(self) result(l)
    class(mcmc), intent(in) :: self
    
    l = self%log_likelihood
    
    return 
  end function mcmc_get_log_likelihood
  
  !---------------------------------------------------------------------

  function mcmc_get_likelihood_saved(self) result(l)
    class(mcmc), intent(in) :: self
    double precision :: l(self%n_iter)
    
    l(:) = self%likelihood_saved(:)
    
    return 
  end function mcmc_get_likelihood_saved
    
  !---------------------------------------------------------------------
  
  function mcmc_get_temp_saved(self) result(t)
    class(mcmc), intent(in) :: self
    double precision :: t(self%n_iter)
    
    t(:) = self%temp_saved(:)
    
    return 
  end function mcmc_get_temp_saved
  
  !---------------------------------------------------------------------

  function mcmc_get_n_accept(self) result(n_accept)
    class(mcmc), intent(in) :: self
    integer :: n_accept(self%n_proposal_type)

    n_accept = self%n_accept

    return 
  end function mcmc_get_n_accept
  
  !---------------------------------------------------------------------

  function mcmc_get_n_propose(self) result(n_propose)
    class(mcmc), intent(in) :: self
    integer :: n_propose(self%n_proposal_type)

    n_propose = self%n_propose

    return 
  end function mcmc_get_n_propose
  
  !---------------------------------------------------------------------

  subroutine mcmc_finish(self)
    class(mcmc), intent(inout) :: self

    self%i_iter = 0
    self%i_mod  = 0
    self%temp = 1.d0
    
    deallocate(self%k_saved)
    deallocate(self%likelihood_saved)
    deallocate(self%temp_saved)
    deallocate(self%n_propose)
    deallocate(self%n_accept)

    return 
  end subroutine mcmc_finish

  !---------------------------------------------------------------------

end module cls_mcmc
