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
module mod_mcmc
  use mod_random
  use mod_trans_d_model
  implicit none 
  
  type mcmc 
     private
     
     type(trans_d_model) :: tm
     integer :: n_iter
     integer :: n_mod
     integer :: i_mod
     integer :: i_iter
     integer :: i_proposal_type
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
     procedure :: get_is_accepted => mcmc_get_is_accepted
     procedure :: get_n_iter => mcmc_get_n_iter
     procedure :: set_temp => mcmc_set_temp
     procedure :: get_temp => mcmc_get_temp
     procedure :: get_log_likelihood => mcmc_get_log_likelihood
     procedure :: get_likelihood_saved => mcmc_get_likelihood_saved
     procedure :: get_temp_saved => mcmc_get_temp_saved
     procedure :: finish => mcmc_finish
  end type mcmc
  
  interface mcmc
     module procedure init_mcmc
  end interface mcmc


contains

  !---------------------------------------------------------------------
  
  type(mcmc) function init_mcmc(tm, n_iter) result(self)
    type(trans_d_model), intent(in) :: tm
    integer, intent(in) :: n_iter
    
    
    self%tm = tm
    self%n_iter = n_iter
    self%i_iter = 0
    self%i_mod = 0
    self%log_likelihood = -1.0d300

    
    allocate(self%likelihood_saved(n_iter))
    allocate(self%temp_saved(n_iter))
    allocate(self%k_saved(n_iter))

    
    return 
  end function init_mcmc

  !---------------------------------------------------------------------
  
  subroutine mcmc_propose_model(self, tm_proposed, is_ok, &
      & log_prior_ratio, log_proposal_ratio)
    class(mcmc), intent(inout) :: self
    type(trans_d_model), intent(out) :: tm_proposed
    double precision, intent(out) :: log_prior_ratio
    double precision, intent(out) :: log_proposal_ratio
    
    logical, intent(out) :: is_ok
    integer :: nparam, iparam

    nparam = self%tm%get_n_x()
    self%i_proposal_type = int(rand_u() * (nparam + 2))
    tm_proposed = self%tm
    if (self%i_proposal_type == 0) then
       ! Birth
       call tm_proposed%birth(is_ok)
       ! prior & proposal ratio 
       ! (can be ignored if birth prorposal from prior distribution)
       log_prior_ratio = 0.d0
       log_proposal_ratio = 0.d0
    else if (self%i_proposal_type == 1) then
       ! Death
       call tm_proposed%death(is_ok)
       log_prior_ratio = 0.d0
       log_proposal_ratio = 0.0
    else
       ! Perturbation
       iparam = self%i_proposal_type - 1
       call tm_proposed%perturb(iparam, is_ok, log_prior_ratio)
       log_proposal_ratio = 0.d0
    end if
    return
    
  end subroutine mcmc_propose_model

  !---------------------------------------------------------------------

  subroutine mcmc_judge_model(self, tm, log_likelihood, &
       & log_prior_ratio, log_proposal_ratio)
    class(mcmc), intent(inout) :: self
    type(trans_d_model), intent(in) :: tm
    double precision, intent(in) :: log_likelihood, log_prior_ratio, &
         & log_proposal_ratio
    double precision :: ratio
    double precision :: r

    
    self%is_accepted = .false.

    ratio = (log_likelihood - self%log_likelihood) / self%temp
    ratio = ratio + log_prior_ratio + log_proposal_ratio
    !write(*,*)log_likelihood, self%log_likelihood
    r = log(rand_u())
    if (r <= ratio) then
       self%is_accepted = .true.
    end if

    if (self%is_accepted) then
       ! Accept model
       self%tm = tm
       self%log_likelihood = log_likelihood
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
    
    l = self%likelihood_saved
    
    return 
  end function mcmc_get_likelihood_saved
    
  !---------------------------------------------------------------------
  
  function mcmc_get_temp_saved(self) result(t)
    class(mcmc), intent(in) :: self
    double precision :: t(self%n_iter)
    
    t = self%temp_saved
    
    return 
  end function mcmc_get_temp_saved
  
  !---------------------------------------------------------------------

  subroutine mcmc_finish(self)
    class(mcmc), intent(inout) :: self

    self%i_iter = 0
    self%i_mod  = 0
    self%temp = 1.d0
    
    deallocate(self%k_saved)
    deallocate(self%likelihood_saved)

    return 
  end subroutine mcmc_finish

  !---------------------------------------------------------------------

end module mod_mcmc
