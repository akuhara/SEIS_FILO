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
  use mod_parallel
  use mod_random
  use mod_trans_d_model
  use mod_mcmc
  use mod_rayleigh
  use mod_interpreter
  use mod_const
  use mod_observation
  use mod_param
  implicit none 
  !include 'mpif.h'

  integer, parameter :: n_chain = 5
  integer, parameter :: n_cool  = 1
  double precision, parameter :: temp_high = 30.d0

  integer, parameter :: k_min = 1, k_max = 21
  integer, parameter :: n_rx = 3
  double precision, parameter :: vs_min = 2.5d0, vs_max = 5.0d0
  double precision, parameter :: vp_min = 5.0d0, vp_max = 8.5d0
  double precision, parameter :: z_min = 0.d0, z_max = 30.d0
  double precision, parameter :: dev_vs = 0.03d0
  double precision, parameter :: dev_vp = 0.03d0
  double precision, parameter :: dev_z  = 0.03d0
  integer, parameter :: nbin_z = 50, nbin_vs = 50
  

  logical, parameter :: solve_vp = .true.
  logical :: verb
  logical, parameter :: ocean_flag = .false.
  double precision :: ocean_thick = 1.d0

  double precision, parameter :: cmin = 0.2d0 * vs_min, &
       &  cmax = 4.5d0, dc = 0.005d0
  double precision :: fmin, fmax, df
  
  integer :: i, j, ierr, n_proc, rank, io_vz, io_ray, n_arg

  double precision :: log_likelihood, temp
  logical :: is_ok
  type(vmodel) :: vm
  type(trans_d_model) :: tm
  type(trans_d_model) ::  tm_tmp
  type(interpreter) :: intpr
  type(mcmc) :: mc
  type(rayleigh) :: ray, ray_tmp
  type(observation) :: obs
  type(parallel) :: pt
  type(param) :: para
  character(200) :: filename, param_file
  
  ! MPI 
  call mpi_init(ierr)
  call mpi_comm_size(MPI_COMM_WORLD, n_proc, ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)
  pt = init_parallel(n_proc=n_proc, rank=rank, n_chain = n_chain)
  if (rank == 0) then
     verb = .true.
  else
     verb = .false.
  end if
   
  ! Get parameter file name from command line argument
  n_arg = command_argument_count()
  if (n_arg /= 1) then
     write(0, *)"USAGE: rayleigh_inv [parameter file]"
     stop
  end if
  call get_command_argument(1, param_file)
 
  ! Read parameter file
  para = init_param(param_file, verb)

  ! Initialize random number sequence
  call init_random(22222322, 2246789, 123147890, 65678901, rank)
  
  ! Read observation file
  obs = init_observation("rayobs.in")
  fmin = obs%get_fmin()
  df   = obs%get_df()
  fmax = fmin + df * (obs%get_nf() - 1)

  ! Set model parameter & generate initial sample
  do i = 1, n_chain
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
     call pt%set_tm(i, tm)
     call tm%finish()
  end do

  ! Set interpreter 
  intpr = init_interpreter(nlay_max=k_max, &
       & z_min=z_min, z_max=z_max, nbin_z=nbin_z, &
       & vs_min=vs_min, vs_max=vs_max, nbin_vs=nbin_vs, &
       & ocean_flag =ocean_flag, ocean_thick=ocean_thick, &
       & solve_vp=solve_vp)
  vm = intpr%get_vmodel(pt%get_tm(1))

  ! Set forward computation
  ray = init_rayleigh(vm=vm, fmin=obs%fmin, fmax=fmax, df=df, &
       cmin=cmin, cmax=cmax, dc=dc)
  
  ! Set MCMC chain
  do i = 1, n_chain
     ! Set transdimensional model
     mc = init_mcmc(pt%get_tm(i), para%get_n_iter())
     ! Set temperatures
     if (i <= n_cool) then
        call mc%set_temp(1.d0)
     else
        temp = exp(rand_u() * log(temp_high))
        call mc%set_temp(temp)
     end if
     call pt%set_mc(i, mc)
  end do
  
  ! Output files
  write(filename,'(A10,I3.3)')"vz_models.", rank
  open(newunit=io_vz, file=filename, status="unknown", iostat=ierr)
  if (ierr /= 0) then
     write(0,*)"ERROR: cannot create ", trim(filename)
     call mpi_finalize(ierr)
     stop
  end if

  write(filename,'(A8,I3.3)')"syn_ray.", rank
  open(newunit=io_ray, file=filename, status="unknown", iostat=ierr)
  if (ierr /= 0) then
     write(0,*)"ERROR: cannot create ", trim(filename)
     call mpi_finalize(ierr)
     stop
  end if
     
  call mpi_finalize(ierr)
  stop

  ! Main
  do i = 1, para%get_n_iter()
     write(*,*)i
     ! MCMC step
     do j = 1, n_chain
        
        mc = pt%get_mc(j)

        ! Proposal
        call mc%propose_model(tm_tmp, is_ok)

        ! Forward computation
        ray_tmp = ray
        if (is_ok) then
           call forward_rayleigh(tm_tmp, intpr, obs, ray_tmp, log_likelihood)
        else
           log_likelihood = -1.d300
        end if

        ! Judege
        call mc%judge_model(tm_tmp, log_likelihood)
        if (mc%get_is_accepted()) then
           ray = ray_tmp
        end if
        call pt%set_mc(j, mc)

        ! Output
        ! ** One step summary
        if (pt%get_rank() == 0) then
           call mc%one_step_summary()
        end if
        
        ! ** Recording
        if (i > para%get_n_burn() .and. &
             & mod(i, para%get_n_corr()) == 0 .and. &
             & mc%get_temp() < 1.d0 + 1.0e-8) then
           tm = mc%get_tm()
           
           ! V-Z
           call intpr%output_vz(tm, io_vz)

           ! Synthetic data
           call ray%output(io_ray)
        end if
     end do

     ! Swap temperture
     call pt%swap_temperature(verb=.true.)
     
  end do
  close(io_vz)
  close(io_ray)


  
  ! Output (First chain only)
  if (pt%get_rank() == 0) then
     do i = 1, n_chain
        mc = pt%get_mc(i)
        
        ! K history
        write(filename, '(A10,I3.3)')'k_history.', i
        call mc%output_k_history(filename)
        
        ! Likelihood history
        write(filename, '(A19,I3.3)')'likelihood_history.', i
        call mc%output_likelihood_history(filename)
     end do
  end if
  call mpi_finalize(ierr)

  stop
end program main


!-----------------------------------------------------------------------

subroutine forward_rayleigh(tm, intpr, obs, ray, log_likelihood)
  use mod_trans_d_model
  use mod_interpreter
  use mod_observation
  use mod_rayleigh
  use mod_vmodel
  implicit none 
  type(trans_d_model), intent(in) :: tm
  type(interpreter), intent(inout) :: intpr
  type(observation), intent(in) :: obs
  type(rayleigh), intent(inout) :: ray
  double precision, intent(out) :: log_likelihood
  type(vmodel) :: vm
  integer :: i
  
  ! calculate synthetic dispersion curves
  vm = intpr%get_vmodel(tm)
  call ray%set_vmodel(vm)
  call ray%dispersion()
  
  ! calc misfit
  log_likelihood = 0.d0
  do i = 1, obs%get_nf()
     log_likelihood = &
          & log_likelihood - (ray%get_c(i) - obs%get_c(i)) ** 2 / &
          & (obs%get_sig_c(i) ** 2)
     log_likelihood = &
          & log_likelihood - (ray%get_u(i) - obs%get_u(i)) ** 2 / &
          & (obs%get_sig_u(i) ** 2)

  end do
  log_likelihood = 0.5d0 * log_likelihood

  return 
end subroutine forward_rayleigh


!-----------------------------------------------------------------------
  
