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
  use mod_random
  use mod_trans_d_model
  use mod_mcmc
  use mod_rayleigh
  use mod_interpreter
  use mod_const
  implicit none 
  include 'mpif.h'

  integer, parameter :: n_iter = 1000000
  integer, parameter :: k_min = 1, k_max = 21
  integer, parameter :: n_rx = 3
  double precision, parameter :: vs_min = 2.5d0, vs_max = 5.0d0
  double precision, parameter :: vp_min = 4.0d0, vp_max = 9.0d0
  double precision, parameter :: z_min = 0.d0, z_max = 50.d0
  double precision, parameter :: dev_vs = 0.05d0
  double precision, parameter :: dev_vp = 0.05d0
  double precision, parameter :: dev_z  = 0.2d0
  
  logical, parameter :: solve_vp = .true.
  logical, parameter :: ocean_flag = .true.
  double precision :: ocean_thick = 1.d0

  double precision, parameter :: fmin = 0.1d0, fmax = 1.d0, df = 0.1d0
  double precision, parameter :: cmin = 0.2d0, cmax = 1.8d0, dc = 0.01d0

  
  integer :: i
  logical :: is_ok
  type(vmodel) :: vm
  type(trans_d_model) :: tm, tm_tmp
  type(interpreter) :: intpr
  type(mcmc) :: mc
  type(rayleigh) :: ray


  
  call init_random(12345678, 23456789, 34567890, 45678901)
  
  ! Set model parameter & generate initial sample
  tm = init_trans_d_model(k_min=k_min, k_max=k_max, n_rx=n_rx)
  call tm%set_prior(id_vs, id_uni, vs_min, vs_max)
  call tm%set_prior(id_vp, id_uni, vp_min, vp_max)
  call tm%set_prior(id_z,  id_uni, z_min,  z_max )
  call tm%set_birth(id_vs, id_uni, vs_min, vs_max)
  call tm%set_birth(id_vp, id_uni, vp_min, vp_max)
  call tm%set_birth(id_z,  id_uni, z_min,  z_max )
  call tm%set_perturb(id_vs, dev_vs)
  call tm%set_perturb(id_vp, dev_vp)
  call tm%set_perturb(id_z,  dev_z)
  call tm%generate_model()

  ! Set interpreter 
  intpr = init_interpreter(nlay_max=k_max, &
       & ocean_flag =ocean_flag, ocean_thick=ocean_thick, &
       & solve_vp=solve_vp)
  vm = intpr%get_vmodel(tm)

  ! Set MCMC chain
  mc = init_mcmc(tm)
  
  ! Set forward computation
  ray = init_rayleigh(vm=vm, fmin=fmin, fmax=fmax, df=df, &
       cmin=cmin, cmax=cmax, dc=dc)
  
  ! Main
  do i = 1, n_iter
     call mc%propose_model(tm_tmp, is_ok)
     vm = intpr%get_vmodel(tm_tmp)
     call ray%set_vmodel(vm)
     
     if (is_ok) then
        call mc%accept_model(tm_tmp)
     else
        
     end if

  end do

  call vm%display()
  

  stop
end program main
