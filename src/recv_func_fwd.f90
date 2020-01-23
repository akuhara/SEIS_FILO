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
  use mod_signal_process
  use mod_random
  implicit none   
  integer :: n_arg
  character(len=200) :: param_file
  type(param) :: para
  type(vmodel) :: vm
  type(signal_process) :: sp
  
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


  ! Test FFT
  block 
    integer, parameter :: n = 512
    double precision, parameter :: delta = 0.05d0
    integer :: i
    double precision :: xt(n)

    do i = 1, n
       xt(i) = 2.d0 * rand_u() - 1.d0
    end do
    
    sp = init_signal_process(n, delta)
    call sp%set_t_data(xt)
    call sp%set_gaussian_filter(a_gauss=8.d0)
    
    call sp%forward_fft()
    call sp%apply_filter()
    call sp%inverse_fft()
    do i = 1, n
       write(30, *) i*delta, xt(i)
    end do
    
    
  end block
  
  stop
end program main
