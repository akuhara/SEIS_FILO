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
module mod_random
  implicit none 

  integer, private :: x, y, z, w  
  double precision, parameter, private :: r31    = 2.d0 ** 31
  double precision, parameter, private :: r32    = 2.d0 ** 32 
  double precision, parameter, private :: pi2 = 2.d0 * acos(-1.d0)

contains

  !---------------------------------------------------------------------

  subroutine init_random(i1, i2, i3, i4, rank)
    integer, intent(in) :: i1, i2, i3, i4
    integer, intent(in), optional :: rank
    integer :: j

    if (present(rank)) then
       j = rank
    else
       j = 0
    end if
    x = i1 * (j + 1)**4 + 1000 * i1 * (j + 1) ** 2 + i1
    y = i2 * (j + 1)**4 + 1000 * i2 * (j + 1) ** 2 + i2
    z = i3 * (j + 1)**4 + 1000 * i3 * (j + 1) ** 2 + i3
    w = i4 * (j + 1)**4 + 1000 * i4 * (j + 1) ** 2 + i4

    return 
  end subroutine init_random

  !---------------------------------------------------------------------
  
  ! U[0, 1)
  double precision function rand_u()
    integer :: t
    
    t = ieor(x, ishft(x, 11))
    x = y
    y = z
    z = w
    w = &
         & ieor( &
         & ieor(w, ishft(w, -19)), &
         & ieor(t, ishft(t,  -8)) &
         & )
    rand_u = (dble(w) + r31) / r32
    return 
  end function rand_u

  !---------------------------------------------------------------------
  ! U(0, 1)
  double precision function rand_u2()
    integer :: t
    
    t = ieor(x, ishft(x, 11))
    x = y
    y = z
    z = w
    w = &
         & ieor( &
         & ieor(w, ishft(w, -19)), &
         & ieor(t, ishft(t,  -8)) &
         & )
    rand_u2 = (dble(w) + r31 + 0.5d0) / r32
    return 
  end function rand_u2
  !---------------------------------------------------------------------
  
  double precision function rand_g()
    double precision:: v1, v2
    
    v1 = rand_u2() 
    v2 = rand_u2()
    rand_g = sqrt(-2.d0 * log(v1)) * cos(pi2 * v2)

  end function rand_g
end module mod_random
  
