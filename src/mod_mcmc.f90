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
     
   contains
     procedure :: propose_model => mcmc_propose_model
     procedure :: accept_model  => mcmc_accept_model
  end type mcmc
  
  interface mcmc
     module procedure init_mcmc
  end interface mcmc

contains

  !---------------------------------------------------------------------
  
  type(mcmc) function init_mcmc(tm)
    type(trans_d_model), intent(in) :: tm
    
    init_mcmc%tm = tm
    
    return 
  end function init_mcmc

  !---------------------------------------------------------------------
  
  subroutine mcmc_propose_model(self, tm_proposed, is_ok)
    class(mcmc), intent(inout) :: self
    type(trans_d_model), intent(out) :: tm_proposed
    logical, intent(out) :: is_ok
    integer :: nparam, itype, iparam

    
    nparam = self%tm%get_n_x()
    itype = int(rand_u() * (nparam + 2))
    tm_proposed = self%tm
    if (itype == 0) then
       ! Birth
       call tm_proposed%birth(is_ok)
    else if (itype == 1) then
       ! Death
       call tm_proposed%death(is_ok)
    else
       iparam = itype - 1
       call tm_proposed%perturb(iparam, is_ok)
    end if
    return
    
  end subroutine mcmc_propose_model

  !---------------------------------------------------------------------

  subroutine mcmc_accept_model(self, tm)
    class(mcmc), intent(inout) :: self
    type(trans_d_model), intent(in) :: tm

    self%tm = tm

    return 
  end subroutine mcmc_accept_model

  !---------------------------------------------------------------------

  subroutine mcmc_judge_model(self, tm, rslt)
    class(mcmc), intent(inout) :: self
    logical, intent(out) :: rslt
    type(trans_d_model), intent(in) :: tm


    return 
  end subroutine mcmc_judge_model
  
  !---------------------------------------------------------------------

end module mod_mcmc
