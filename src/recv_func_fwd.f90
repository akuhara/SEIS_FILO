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
program main
  use mod_param
  use mod_vmodel
  use mod_random
  use mod_recv_func
  
  implicit none   
  integer :: n_arg, n
  character(len=200) :: param_file
  type(param) :: para
  type(vmodel) :: vm
  type(recv_func) :: rf
  double precision :: rayp, a_gauss, delta, t_pre
  character(len=1) :: phase
  logical :: deconv_flag, correct_amp
  rayp = 0.05d0
  a_gauss = 8.0d0
  phase = "S"
  deconv_flag = .true.
  n = 1024
  delta = 0.05d0
  t_pre = 30.d0
  correct_amp = .true.
  call init_random(100, 20000, 3000000, 43444)
  
  ! Get parameter file name from command line argument
  n_arg = command_argument_count()
  if (n_arg /= 1) then
     write(0, *)"USAGE: rayleigh_fwd [parameter file]"
     stop
  end if
  call get_command_argument(1, param_file)
  
  ! Read parameter file
  para = init_param(param_file)

  ! Set velocity model
  call vm%read_file(para%get_vmod_in())

  ! Init RF
  rf = init_recv_func(vm, n, delta, rayp, a_gauss, &
       & phase, deconv_flag, t_pre=t_pre, correct_amp = .true.)
  
  ! Main
  call rf%compute()
  
  block 
    double precision :: rft(n), xr(n), xz(n)
    integer :: i
    ! Get result
    rft(:) = rf%get_rf_data()
    xr(:) = rf%get_t_data(1)
    xz(:) = rf%get_t_data(2)
    
    do i = 1, n
       write(111,*) (i  - 1) * delta - t_pre, rft(i)
    end do
    write(111,*)
    write(111,*)

  end block

  stop
end program main
