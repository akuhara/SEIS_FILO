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
     type(trans_d_model), allocatable :: tm_saved(:)
     integer :: n_iter
     integer :: n_corr = 1
     integer :: n_burn = 0
     integer :: n_mod
     integer :: i_mod
     integer :: i_iter
     integer :: i_proposal_type
     double precision :: log_likelihood
     double precision, allocatable :: likelihood_saved(:)
     logical :: is_accepted
   contains
     procedure :: propose_model => mcmc_propose_model
     procedure :: accept_model  => mcmc_accept_model
     procedure :: judge_model   => mcmc_judge_model
     procedure :: one_step_summary => mcmc_one_step_summary
  end type mcmc
  
  interface mcmc
     module procedure init_mcmc
  end interface mcmc


contains

  !---------------------------------------------------------------------
  
  type(mcmc) function init_mcmc(tm, n_iter, n_corr, n_burn)
    type(trans_d_model), intent(in) :: tm
    integer, intent(in) :: n_iter
    integer, intent(in), optional :: n_corr
    integer, intent(in), optional :: n_burn

    init_mcmc%tm = tm
    init_mcmc%n_iter = n_iter
    init_mcmc%i_iter = 0
    init_mcmc%i_mod = 0
   

    if (present(n_burn)) then
       init_mcmc%n_burn = n_burn
    end if

    if (present(n_corr)) then
       init_mcmc%n_corr = n_corr
    end if

    allocate(init_mcmc%likelihood_saved(n_iter))

    init_mcmc%n_mod = int((init_mcmc%n_iter - init_mcmc%n_burn) / &
         & init_mcmc%n_corr)
    allocate(init_mcmc%tm_saved(init_mcmc%n_mod))
    
    init_mcmc%log_likelihood = 0.d0
    
    return 
  end function init_mcmc

  !---------------------------------------------------------------------
  
  subroutine mcmc_propose_model(self, tm_proposed, is_ok)
    class(mcmc), intent(inout) :: self
    type(trans_d_model), intent(out) :: tm_proposed
    logical, intent(out) :: is_ok
    integer :: nparam, iparam

    
    nparam = self%tm%get_n_x()
    self%i_proposal_type = int(rand_u() * (nparam + 2))
    tm_proposed = self%tm
    if (self%i_proposal_type == 0) then
       ! Birth
       call tm_proposed%birth(is_ok)
    else if (self%i_proposal_type == 1) then
       ! Death
       call tm_proposed%death(is_ok)
    else
       iparam = self%i_proposal_type - 1
       call tm_proposed%perturb(iparam, is_ok)
    end if
    return
    
  end subroutine mcmc_propose_model

  !---------------------------------------------------------------------

  subroutine mcmc_accept_model(self, tm, log_likelihood)
    class(mcmc), intent(inout) :: self
    type(trans_d_model), intent(in) :: tm
    double precision, intent(in) :: log_likelihood
        
    self%tm = tm
    self%log_likelihood = log_likelihood
    
    return 
  end subroutine mcmc_accept_model

  !---------------------------------------------------------------------

  subroutine mcmc_judge_model(self, tm, log_likelihood)
    class(mcmc), intent(inout) :: self
    type(trans_d_model), intent(in) :: tm
    double precision, intent(in) :: log_likelihood
    double precision :: ratio
    double precision :: r
    
    self%is_accepted = .false.

    ratio = (log_likelihood - self%log_likelihood)
    r = log(rand_u())
    if (r <= ratio) then
       self%is_accepted = .true.
    end if

    if (self%is_accepted) then
       call self%accept_model(tm, log_likelihood)
    end if


    self%i_iter = self%i_iter + 1
    self%likelihood_saved(self%i_iter) = self%log_likelihood


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

end module mod_mcmc
