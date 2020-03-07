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
  use mod_const
  use cls_interpreter
  use cls_param
  use cls_observation_recv_func
  use cls_recv_func
  use cls_covariance
  use cls_disper
  use cls_observation_disper
  use mod_output
  use mod_forward
  implicit none 

  integer :: n_rx
  logical :: verb
  integer :: i, j, k, ierr, n_proc, rank, n_arg
  integer :: n_mod
  double precision :: log_likelihood, temp, log_likelihood2
  double precision :: log_prior_ratio, log_proposal_ratio
  double precision :: del_amp
  logical :: is_ok
  type(vmodel) :: vm
  type(trans_d_model) :: tm, tm_tmp
  type(interpreter) :: intpr
  type(mcmc) :: mc
  type(parallel) :: pt
  type(param) :: para
  type(observation_recv_func) :: obs_rf
  type(recv_func), allocatable :: rf(:), rf_tmp(:)
  type(covariance), allocatable :: cov(:)
  type(disper), allocatable :: disp(:), disp_tmp(:)
  type(observation_disper) :: obs_disp
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

  
  if (verb) then
     write(*,*)"------------------------------------------------------"
     write(*,*)" SEIS_FILO:                                           "
     write(*,*)" SEISmological inversion tools for Flat and Isotropic "
     write(*,*)" Layered structure in the Ocean                       "
     write(*,*)" Copyright (C) 2019 Takeshi Akuhara                   "
     write(*,*)"------------------------------------------------------"
     write(*,*)
     write(*,*)"Start Program <joint_inv>"
     write(*,*)
  end if

  ! Get parameter file name from command line argument
  n_arg = command_argument_count()
  if (n_arg /= 1) then
     write(0, *)"USAGE: joint_inv [parameter file]"
     stop
  end if
  call get_command_argument(1, param_file)
 
  ! Read parameter file
  para = param(param_file, verb)
  call para%check_mcmc_params(is_ok)
  if (.not. is_ok) then
     if (verb) then
        write(0,*)"ERROR: while checking MCMC parameters"
     end if
     call mpi_finalize(ierr)
     stop
  end if
        
  
  ! Initialize parallel chains
  pt = parallel(                       &
       & n_proc  = n_proc,             &
       & rank    = rank,               &
       & n_chain = para%get_n_chain(), &
       & verb    = verb)
  
  call mpi_finalize(ierr)
  stop
  
  ! Initialize random number sequence
  call init_random(para%get_i_seed1(), &
       &           para%get_i_seed2(), &
       &           para%get_i_seed3(), &
       &           para%get_i_seed4(), &
       &           rank)
  
  ! Read observation file
  if (trim(para%get_recv_func_in()) /= "") then
     if (verb) write(*,*)"Reading observation file"
     obs_rf = observation_recv_func(trim(para%get_recv_func_in()))
  else
     if (verb) write(*,*)"No receiver function input"
     call obs_rf%set_n_rf(0)
  end if
  
  ! Read observation file
  if (trim(para%get_disper_in()) /= "") then
     if (verb) write(*,*)"Reading observation file"
     obs_disp = observation_disper(trim(para%get_disper_in()))
  else
     if (verb) write(*,*)"No dispersion curve input"
     call obs_disp%set_n_disp(0)
  end if
  
  ! Covariance matrix
  if (verb) write(*,*)"Constructing covariance matrix"
  allocate(cov(obs_rf%get_n_rf()))
  do i = 1, obs_rf%get_n_rf()
     cov(i) = covariance(                    &
          & n = obs_rf%get_n_smp(i),         &
          & a_gauss = obs_rf%get_a_gauss(i), &
          & delta = obs_rf%get_delta(i)      &
          & )
  end do
  

  ! Set interpreter 
  if (verb) write(*,*)"Setting interpreter"
  intpr = interpreter(nlay_max= para%get_k_max(), &
       & z_min = para%get_z_min(), z_max = para%get_z_max(), &
       & n_bin_z = para%get_n_bin_z(), &
       & vs_min = para%get_vs_min(), vs_max = para%get_vs_max(), &
       & n_bin_vs = para%get_n_bin_vs(), &
       & vp_min = para%get_vp_min(), vp_max = para%get_vp_max(), &
       & n_bin_vp = para%get_n_bin_vp(), &
       & is_ocean = para%get_is_ocean(), &
       & ocean_thick = para%get_ocean_thick(), &
       & solve_vp = para%get_solve_vp(), &
       & vp_bottom = para%get_vp_bottom(), &
       & vs_bottom = para%get_vs_bottom(), &
       & rho_bottom = para%get_rho_bottom())

  
  ! Set model parameter & generate initial sample
  if (verb) write(*,*)"Setting model parameters"
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
  

  ! Set forward computation (recv_func)
  tm = pt%get_tm(1)
  vm = intpr%get_vmodel(pt%get_tm(1))
  allocate(rf(obs_rf%get_n_rf()), rf_tmp(obs_rf%get_n_rf()))
  do i = 1, obs_rf%get_n_rf()
     rf(i) = recv_func( &
          & vm = vm, &
          & n = obs_rf%get_n_smp(i) * 2, &
          & delta = obs_rf%get_delta(i), &
          & rayp = obs_rf%get_rayp(i), &
          & a_gauss = obs_rf%get_a_gauss(i), &
          & rf_phase = obs_rf%get_rf_phase(i), &
          & deconv_flag = obs_rf%get_deconv_flag(i), &
          & t_pre = -obs_rf%get_t_start(i), &
          & correct_amp = obs_rf%get_correct_amp(i) &
          & )
     
  end do
  
  ! Set forward computation (disper)
  tm = pt%get_tm(1)
  vm = intpr%get_vmodel(pt%get_tm(1))
  allocate(disp(obs_disp%get_n_disp()), &
       & disp_tmp(obs_disp%get_n_disp()))
  do i = 1, obs_disp%get_n_disp()
     disp(i) = disper(vm=vm, &
          & fmin=obs_disp%get_fmin(i), &
          & fmax=obs_disp%get_fmax(i), &
          & df=obs_disp%get_df(i), &
          & cmin=obs_disp%get_cmin(i), &
          & cmax=obs_disp%get_cmax(i), &
          & dc=obs_disp%get_dc(i), &
          & n_mode = obs_disp%get_n_mode(i), &
          & disper_phase = obs_disp%get_disper_phase(i) &
          & )
  end do
  
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
  
  ! Main
  do i = 1, para%get_n_iter()
     ! MCMC step
     do j = 1, para%get_n_chain()
        
        mc = pt%get_mc(j)

        ! Proposal
        call mc%propose_model(tm_tmp, is_ok, log_prior_ratio, &
             & log_proposal_ratio)

        ! Forward computation
        rf_tmp = rf
        disp_tmp = disp
        if (is_ok) then
           call forward_recv_func(tm_tmp, intpr, obs_rf, &
                & rf_tmp, cov, log_likelihood)
           call forward_disper(tm_tmp, intpr, obs_disp, &
                & disp_tmp, log_likelihood2)
           log_likelihood = &
                log_likelihood + log_likelihood2
        else
           log_likelihood = minus_infty
        end if

        ! Judege
        call mc%judge_model(tm_tmp, log_likelihood, &
              & log_prior_ratio, log_proposal_ratio)
        if (mc%get_is_accepted()) then
           rf = rf_tmp
           disp = disp_tmp
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
           call intpr%save_model(tm)
           
           ! Synthetic data
           do k = 1, obs_rf%get_n_rf()
              call rf(k)%save_syn()
           end do
           do k = 1, obs_disp%get_n_disp()
              call disp(k)%save_syn()
           end do
           
        end if
     end do
     
     ! Swap temperture
     call pt%swap_temperature(verb=.true.)
  end do
  
  
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
       & 0.5d0 * intpr%get_dz(), &
       & intpr%get_dz())
  
  ! Mean Vs
  filename = "vs_z.mean"
  call output_mean_model(filename, rank, para%get_n_bin_z(), &
       & intpr%get_vsz_mean(), n_mod, &
       & 0.5d0 * intpr%get_dz(), &
       & intpr%get_dz())

  if (para%get_solve_vp()) then
     ! marginal posterior of Vp
     filename = "vp_z.ppd"
     call output_ppd_2d(filename, rank, para%get_n_bin_vp(), &
          & para%get_n_bin_z(), intpr%get_n_vpz(), n_mod, &
          & para%get_vp_min() + 0.5d0 * intpr%get_dvp(), &
          & intpr%get_dvp(), &
          & 0.5d0 * intpr%get_dz(), &
          & intpr%get_dz())
     ! Mean Vp
     filename = "vp_z.mean"
     call output_mean_model(filename, rank, para%get_n_bin_z(), &
          & intpr%get_vpz_mean(), n_mod, &
          & 0.5d0 * intpr%get_dz(), &
          & intpr%get_dz())
  end if
  
  ! Synthetic RF
  do i = 1, obs_rf%get_n_rf()
     del_amp = (rf(i)%get_amp_max() - rf(i)%get_amp_min()) / &
          & rf(i)%get_n_bin_amp()
     write(filename, '(A6,I3.3,A4)')"syn_rf", i, ".ppd"
     call output_ppd_2d(filename, rank, obs_rf%get_n_smp(i) * 2, &
          & rf(i)%get_n_bin_amp(), rf(i)%get_n_syn_rf(), &
          & n_mod, obs_rf%get_t_start(i), obs_rf%get_delta(i), &
          & rf(i)%get_amp_min(), del_amp)
  end do
  
  ! Synthetic dispersion curves
  do i = 1, obs_disp%get_n_disp()
     write(filename, '(A9,I3.3,A4)')"syn_phase", i, ".ppd"
     call output_ppd_2d(filename, rank, disp(i)%get_nf(), &
          & disp(i)%get_nc(), disp(i)%get_n_fc(), &
          & n_mod, obs_disp%get_fmin(i), obs_disp%get_df(i), &
          & obs_disp%get_cmin(i), obs_disp%get_dc(i))
     
     ! Frequency-group velocity
     write(filename, '(A9,I3.3,A4)')"syn_group", i, ".ppd"
     call output_ppd_2d(filename, rank, disp(i)%get_nf(), &
          & disp(i)%get_nc(), disp(i)%get_n_fu(), &
          & n_mod, obs_disp%get_fmin(i), obs_disp%get_df(i), &
          & obs_disp%get_cmin(i), obs_disp%get_dc(i))
  end do
  
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
