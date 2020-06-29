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
  use cls_param
  use cls_vmodel
  use cls_recv_func
  use mod_random
  implicit none   
  integer :: n_arg
  character(len=200) :: param_file
  type(param) :: para
  type(vmodel) :: vm
  type(recv_func) :: rf
  logical :: is_ok
  
  ! Get parameter file name from command line argument
  n_arg = command_argument_count()
  if (n_arg /= 1) then
     write(0, *)"USAGE: rayleigh_fwd [parameter file]"
     stop
  end if
  call get_command_argument(1, param_file)
  
  ! Read parameter file
  para = param(param_file)
  call para%check_recv_func_fwd_params(is_ok)
  if (.not. is_ok) then
     write(0,*)"ERROR: while checking parameters"
     stop
  end if

  ! Set velocity model
  call vm%read_file(para%get_vmod_in())

  ! Init random generator
  call init_random(para%get_i_seed1(), &
       &           para%get_i_seed2(), &
       &           para%get_i_seed3(), &
       &           para%get_i_seed4())

  
  ! Init RF
  rf = recv_func(&
       & vm    = vm, &
       & n     = para%get_n_smp(), &
       & delta = para%get_delta(), &
       & rayp  = para%get_rayp(), &
       & a_gauss = para%get_a_gauss(), &
       & rf_phase = para%get_rf_phase(), &
       & deconv_flag = para%get_deconv_flag(), &
       & t_pre = para%get_t_pre(), &
       & correct_amp = para%get_correct_amp(), &
       & noise_added = para%get_noise_added() &
       & )
  
  ! Main
  call rf%compute()
  if (para%get_noise_added() > 0.d0) then
     call rf%add_noise()
  end if

  ! Output
  call rf%output_sac(para%get_recv_func_out())

  stop
end program main
