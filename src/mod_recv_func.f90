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
module mod_recv_func
  use mod_vmodel
  use mod_signal_process
  implicit none 
  
  type recv_func
     private
     type(vmodel) :: vmodel
     double precision :: rayp
     double precision :: a_gauss
     character(len=1) :: phase
     
     logical :: deconv_flag = .true.
     
   contains
     procedure :: set_vmodel => recv_func_set_vmodel
     
  end type recv_func
  
  interface recv_func
     module procedure :: init_recv_func
  end interface recv_func
    
contains

  !---------------------------------------------------------------------

  type(recv_func) function init_recv_func(vm, rayp, a_gauss, phase, &
       & deconv_flag) result(self)
    type(vmodel), intent(in) :: vm
    double precision, intent(in) :: rayp
    double precision, intent(in) :: a_gauss
    character(*), intent(in) :: phase
    logical, intent(in), optional :: deconv_flag

    self%vmodel = vm
    self%rayp = rayp
    self%a_gauss = a_gauss
    self%phase = phase
    if (present(deconv_flag)) then
       self%deconv_flag = deconv_flag
    end if
    
    return 
  end function init_recv_func

  !---------------------------------------------------------------------

  subroutine recv_func_set_vmodel(self, vm)
    class(recv_func), intent(inout) :: self
    type(vmodel), intent(in) :: vm
    
    self%vmodel = vm
    
    return 
  end subroutine recv_func_set_vmodel

  !---------------------------------------------------------------------

end module mod_recv_func
