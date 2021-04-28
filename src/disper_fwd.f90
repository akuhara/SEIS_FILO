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
  use cls_disper
  use mod_random
  implicit none 
  include 'mpif.h'

  integer :: n_arg, ierr
  character(len=200) :: param_file
  type(param) :: para
  type(vmodel) :: vm, vm2
  type(disper) :: disp
  logical :: verb, is_ok

  verb = .true.
  
  ! MPI init
  !  just in case of erroneous exit, where mpi_finalize is called.
  call mpi_init(ierr) 
  
  ! Get parameter file name from command line argument
  n_arg = command_argument_count()
  if (n_arg /= 1) then
     error stop "USAGE: disper_fwd [parameter file]"
  end if
  call get_command_argument(1, param_file)

  ! Read parameter file
  para = init_param(param_file, verb=verb)
  
  ! Set velocity model
  vm = init_vmodel()
  call vm%read_file(para%get_vmod_in())
  if (para%get_is_sphere()) then
     call vm%sphere2flat(para%get_r_earth(), vm2)
  else 
     vm2 = vm
  end if     
  
  ! Init random generator
  call init_random(para%get_i_seed1(), &
       &           para%get_i_seed2(), &
       &           para%get_i_seed3(), &
       &           para%get_i_seed4())

    ! Calculate dispersion curve
  disp = disper(&
       & vm     = vm2, &
       & fmin   = para%get_fmin(), &
       & fmax   = para%get_fmax(), &
       & df     = para%get_df(), &
       & cmin   = para%get_cmin(), &
       & cmax   = para%get_cmax(), &
       & dc     = para%get_dc(), &
       & disper_phase = para%get_disper_phase(), &
       & n_mode = para%get_n_mode(), &
       & disper_out = para%get_disper_out(), &
       & noise_added = para%get_noise_added() &
       & )
  
  call disp%dispersion(is_ok=is_ok)
  if (.not. is_ok) then
     error stop 
  end if
    
  stop
end program main
