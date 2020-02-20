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
  use cls_parallel  
  use mod_random
  use cls_trans_d_model
  use cls_mcmc
  use cls_rayleigh
  use mod_const
  use cls_interpreter
  use cls_observation_disper
  use cls_param
  use mod_output
  implicit none 

  integer :: n_rx
  logical :: verb
  double precision :: fmin, fmax, df
  integer :: i, j, ierr, n_proc, rank, io_vz, io_ray, n_arg
  integer :: n_mod
  double precision :: log_likelihood, temp
  double precision :: log_prior_ratio, log_proposal_ratio
  logical :: is_ok
  type(vmodel) :: vm
  type(trans_d_model) :: tm, tm_tmp
  type(interpreter) :: intpr
  type(mcmc) :: mc
  type(rayleigh) :: ray, ray_tmp
  type(observation_disper) :: obs
  type(parallel) :: pt
  type(param) :: para
  character(200) :: filename, param_file

  
  ! Initialize MPI 
  call mpi_init(ierr)
  call mpi_comm_size(MPI_COMM_WORLD, n_proc, ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)
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
  call para%check_mcmc_params(is_ok)
  if (.not. is_ok) then
     if (verb) then
        write(0,*)"ERROR: while checking MCMC parameters"
     end if
     call mpi_finalize(ierr)
     stop
  end if
  
  ! Initialize parallel chains
  pt = init_parallel(n_proc = n_proc, rank = rank, &
       & n_chain = para%get_n_chain())
  
  ! Initialize random number sequence
  call init_random(para%get_i_seed1(), &
       &           para%get_i_seed2(), &
       &           para%get_i_seed3(), &
       &           para%get_i_seed4(), &
       &           rank)
  
  ! Read observation file
  obs = init_observation_disper(trim(para%get_disper_in()))
  fmin = obs%get_fmin()
  df   = obs%get_df()
  fmax = fmin + df * (obs%get_nf() - 1)

  ! Set interpreter 
  write(*,*)"Setting interpreter"
  intpr = init_interpreter(nlay_max= para%get_k_max(), &
       & z_min = para%get_z_min(), z_max = para%get_z_max(), &
       & n_bin_z = para%get_n_bin_z(), &
       & vs_min = para%get_vs_min(), vs_max = para%get_vs_max(), &
       & n_bin_vs = para%get_n_bin_vs(), &
       & vp_min = para%get_vp_min(), vp_max = para%get_vp_max(), &
       & n_bin_vp = para%get_n_bin_vp(), &
       & is_ocean = para%get_is_ocean(), &
       & ocean_thick = para%get_ocean_thick(), &
       & solve_vp = para%get_solve_vp())

  ! Set model parameter & generate initial sample
  write(*,*)"Setting model parameters"
  if (para%get_solve_vp()) then
     n_rx = 3
  else 
     n_rx = 2
  end if
  do i = 1, para%get_n_chain()
     tm = init_trans_d_model( &
          & k_min = para%get_k_min(), &
          & k_max = para%get_k_max(), &
          & n_rx=n_rx)
     call tm%set_prior(id_vs, id_uni, &
          & para%get_vs_min(), para%get_vs_max())
     call tm%set_prior(id_z,  id_uni, &
          & para%get_z_min(), para%get_z_max())
     call tm%set_birth(id_vs, id_uni, &
          & para%get_vs_min(), para%get_vs_max())
     call tm%set_birth(id_z,  id_uni, &
          & para%get_z_min(), para%get_z_max())
     call tm%set_perturb(id_vs, para%get_dev_vs())

     call tm%set_perturb(id_z,  para%get_dev_z())
     if (para%get_solve_vp()) then
        call tm%set_prior(id_vp, id_uni, &
             & para%get_vp_min(), para%get_vp_max()) 
        call tm%set_birth(id_vp, id_uni, &
             & para%get_vp_min(), para%get_vp_max())
        call tm%set_perturb(id_vp, para%get_dev_vp())
     end if
     call tm%generate_model()
     call pt%set_tm(i, tm)
     call tm%finish()
  end do
  ! Set forward computation
  tm = pt%get_tm(1)
  vm = intpr%get_vmodel(pt%get_tm(1))
  ray = init_rayleigh(vm=vm, fmin=obs%get_fmin(), fmax=fmax, df=df, &
         & cmin=para%get_cmin(), cmax=para%get_cmax(), &
         & dc=para%get_dc())

  ! Set MCMC chain
  do i = 1, para%get_n_chain()
     ! Set transdimensional model
     mc = init_mcmc(pt%get_tm(i), para%get_n_iter())
     ! Set temperatures
     if (i <= para%get_n_cool()) then
        call mc%set_temp(1.d0)
     else
        temp = exp((rand_u() *(1.d0 - eps) + eps) &
             & * log(para%get_temp_high()))
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
  
  ! Main
  do i = 1, para%get_n_iter()
     ! MCMC step
     do j = 1, para%get_n_chain()
        
        mc = pt%get_mc(j)

        ! Proposal
        call mc%propose_model(tm_tmp, is_ok, log_prior_ratio, &
             & log_proposal_ratio)

        ! Forward computation
        ray_tmp = ray
        if (is_ok) then
           call forward_rayleigh(tm_tmp, intpr, obs, &
                & ray_tmp, log_likelihood)
        else
           log_likelihood = minus_infty
        end if

        ! Judege
        call mc%judge_model(tm_tmp, log_likelihood, &
              & log_prior_ratio, log_proposal_ratio)
        if (mc%get_is_accepted()) then
           ray = ray_tmp
        end if
        call pt%set_mc(j, mc)

        ! One step summary
        if (pt%get_rank() == 0) then
           call mc%one_step_summary()
        end if
        
        ! Recording
        if (i > para%get_n_burn() .and. &
             & mod(i, para%get_n_corr()) == 0 .and. &
             & mc%get_temp() < 1.d0 + eps) then
           tm = mc%get_tm()
           
           ! V-Z
           call intpr%save_model(tm, io_vz)
           
           ! Synthetic data
           call ray%save_syn()
        end if
     end do
     ! Swap temperture
     call pt%swap_temperature(verb=.true.)
  end do
  close(io_vz)
  close(io_ray)

  ! Post processing following the main loop
  
  ! Total model number
  n_mod = para%get_n_cool() * n_proc * &
       & (para%get_n_iter() - para%get_n_burn()) / &
       & para%get_n_corr()

  ! K
  filename = "n_layers.ppd"
  call output_ppd_1d(filename, rank, para%get_k_max(), &
       & intpr%get_n_layers(), n_mod, &
       & dble(para%get_k_min()), 1.d0)
  
  ! marginal posterior of Vs 
  filename = "vs_z.ppd"
  call output_ppd_2d(filename, rank, para%get_n_bin_vs(), &
       & para%get_n_bin_z(), intpr%get_n_vsz(), n_mod, &
       & para%get_vs_min() + 0.5d0 * intpr%get_dvs(), &
       & intpr%get_dvs(), &
       & para%get_z_min() + 0.5d0 * intpr%get_dz(), &
       & intpr%get_dz())
  
  ! Mean Vs
  filename = "vs_z.mean"
  call output_mean_model(filename, rank, para%get_n_bin_z(), &
       & intpr%get_vsz_mean(), n_mod, &
       & para%get_z_min() + 0.5d0 * intpr%get_dz(), &
       & intpr%get_dz())

  if (para%get_solve_vp()) then
     ! marginal posterior of Vp
     filename = "vp_z.ppd"
     call output_ppd_2d(filename, rank, para%get_n_bin_vp(), &
          & para%get_n_bin_z(), intpr%get_n_vpz(), n_mod, &
          & para%get_vp_min() + 0.5d0 * intpr%get_dvp(), &
          & intpr%get_dvp(), &
          & para%get_z_min() + 0.5d0 * intpr%get_dz(), &
          & intpr%get_dz())
     ! Mean Vp
     filename = "vp_z.mean"
     call output_mean_model(filename, rank, para%get_n_bin_z(), &
          & intpr%get_vpz_mean(), n_mod, &
          & para%get_z_min() + 0.5d0 * intpr%get_dz(), &
          & intpr%get_dz())
  end if

  ! Frequency-phase velocity
  filename = "f_c.ppd"
  call output_ppd_2d(filename, rank, ray%get_nf(), &
       & ray%get_nc(), ray%get_n_fc(), &
       & n_mod, obs%get_fmin(), obs%get_df(), &
       & para%get_cmin(), para%get_dc())

  ! Frequency-group velocity
  filename = "f_u.ppd"
  call output_ppd_2d(filename, rank, ray%get_nf(), &
       & ray%get_nc(), ray%get_n_fu(), &
       & n_mod, obs%get_fmin(), obs%get_df(), &
       & para%get_cmin(), para%get_dc())

  ! Likelihood history
  filename = "likelihood.history"
  call pt%output_history(filename, 'l')
  filename = "temp.history"
  call pt%output_history(filename, 't')
  
  ! Proposal count
  filename = "proposal.count"
  call pt%output_proposal(filename)
  
  call mpi_finalize(ierr)
  
  stop
end program main

!-----------------------------------------------------------------------

subroutine forward_rayleigh(tm, intpr, obs, ray, log_likelihood)
  use cls_trans_d_model
  use cls_interpreter
  use cls_observation_disper
  use cls_rayleigh
  use cls_vmodel
  use mod_const, only: minus_infty
  
  implicit none 
  type(trans_d_model), intent(in) :: tm
  type(interpreter), intent(inout) :: intpr
  type(observation_disper), intent(in) :: obs
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
     if (ray%get_c(i) == 0.d0 .or. ray%get_u(i) == 0.d0) then
        log_likelihood = minus_infty
        exit
     end if
     if (obs%get_sig_c(i) > 0.d0) then
        log_likelihood = &
             & log_likelihood - (ray%get_c(i) - obs%get_c(i)) ** 2 / &
             & (obs%get_sig_c(i) ** 2)
     end if
     if (obs%get_sig_u(i) > 0.d0) then
        log_likelihood = &
             & log_likelihood - (ray%get_u(i) - obs%get_u(i)) ** 2 / &
             & (obs%get_sig_u(i) ** 2)
     end if

  end do
  log_likelihood = 0.5d0 * log_likelihood

  return 
end subroutine forward_rayleigh


!-----------------------------------------------------------------------
  
