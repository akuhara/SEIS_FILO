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
module mod_interpreter
  use mod_trans_d_model
  use mod_vmodel
  use mod_sort
  use mod_const
  implicit none 
  
  type interpreter
     private
     integer :: nlay_max

     logical :: ocean_flag = .false.
     double precision :: ocean_thick = 0.d0
     double precision :: ocean_vp    = 1.5d0
     double precision :: ocean_rho   = 1.d0

     logical :: solve_vp = .true.

     double precision, allocatable :: wrk_vp(:), wrk_vs(:), wrk_z(:) 
   contains
     procedure :: get_vmodel => interpreter_get_vmodel
  end type interpreter
  
  interface interpreter
     module procedure init_interpreter
  end interface interpreter
  
contains

  !---------------------------------------------------------------------
  type(interpreter) function init_interpreter(nlay_max, &
       & ocean_flag, ocean_thick, ocean_vp, ocean_rho, solve_vp)
    integer, intent(in) :: nlay_max
    logical, intent(in), optional :: ocean_flag
    double precision, intent(in), optional :: ocean_thick, &
         ocean_vp, ocean_rho
    logical, intent(in), optional :: solve_vp

    init_interpreter%nlay_max = nlay_max
    allocate(init_interpreter%wrk_vp(nlay_max))
    allocate(init_interpreter%wrk_vs(nlay_max))
    allocate(init_interpreter%wrk_z(nlay_max))

    if (present(ocean_flag)) then
       init_interpreter%ocean_flag = ocean_flag
    end if

    if (present(ocean_thick)) then
       init_interpreter%ocean_thick = ocean_thick
    end if

    if (present(ocean_vp)) then
       init_interpreter%ocean_vp = ocean_vp
    end if
    
    if (present(ocean_rho)) then
       init_interpreter%ocean_rho = ocean_rho
    end if

    if (present(solve_vp)) then
       init_interpreter%solve_vp = solve_vp
    end if

    return 
  end function init_interpreter

  !---------------------------------------------------------------------
  
  type(vmodel) function interpreter_get_vmodel(self, tm) result(vm)
    class(interpreter), intent(inout) :: self
    type(trans_d_model), intent(in) :: tm
    integer :: i, i1, k

    
    k = tm%get_k()
    vm = init_vmodel()

    ! Set ocean layer
    if (self%ocean_flag) then
       call vm%set_nlay(k + 1) ! k solid layers + 1 ocean layer
       call vm%set_vp(1, self%ocean_vp)
       call vm%set_vs(1, -999.d0)
       call vm%set_rho(1, self%ocean_rho)
       call vm%set_h(1, self%ocean_thick)
       i1 = 1
    else
       call vm%set_nlay(k)
       i1 = 0
    end if
    
    ! Translate trand_d_model to vmodel
    self%wrk_z(:) = tm%get_rx(id_z)
    self%wrk_vs(:) = tm%get_rx(id_vs)
    if (self%solve_vp) then
       self%wrk_vp(:) = tm%get_rx(id_vp)
    end if
    call quick_sort(self%wrk_z, 1, k, self%wrk_vs, self%wrk_vp)
    ! Middle layers
    do i = 1, k - 1
       ! Vs
       call vm%set_vs(i+i1, self%wrk_vs(i))
       ! Vp
       if (self%solve_vp) then
          call vm%set_vp(i+i1, self%wrk_vp(i))
       else
          call vm%vs2vp_brocher(i+i1)
       end if
       ! Thickness
       call vm%set_h(i+i1,  self%wrk_z(i+1) - self%wrk_z(i))
       ! Density
       call vm%vp2rho_brocher(i+i1)
    end do
    ! Bottom layer
    ! Vs
    call vm%set_vs(k+i1, 4.6d0) ! <- Fixed
    ! Vp
    call vm%vs2vp_brocher(k+i1)
    ! Thickness
    call vm%set_h(k+i1,  -99.d0) ! <- half space
    ! Density
    call vm%vp2rho_brocher(k+i1)
    
    return 
  end function interpreter_get_vmodel
  
  !---------------------------------------------------------------------
  
end module mod_interpreter

