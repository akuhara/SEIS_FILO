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
  use cls_proposal
  implicit none 

  integer :: n_td, io_vmod_all
  logical :: verb
  integer :: i, j, k, ierr, n_proc, rank, n_arg
  integer :: n_mod
  double precision :: log_likelihood, temp, log_likelihood2
  double precision :: log_prior_ratio, log_proposal_ratio
  double precision :: del_amp
  logical :: is_ok, is_ok2
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
  type(proposal) :: prop
  type(hyper_model) :: hyp_rf, hyp_disp, hyp_rf_tmp, hyp_disp_tmp
  character(200) :: filename, param_file
  

  !---------------------------------------------------------------------
  ! 1. Initialize
  !---------------------------------------------------------------------

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
     if (verb) then
        write(0, *)"USAGE: joint_inv [parameter file]"
     end if
     call mpi_finalize(ierr)
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
        

  
  ! Initialize random number sequence
  call init_random(para%get_i_seed1(), &
       &           para%get_i_seed2(), &
       &           para%get_i_seed3(), &
       &           para%get_i_seed4(), &
       &           rank)
  
  ! Read observation file
  if (trim(para%get_recv_func_in()) /= "") then
     obs_rf = observation_recv_func(para%get_recv_func_in(), verb)
  else
     if (verb) write(*,*)"<< No receiver function input >>"
     call obs_rf%set_n_rf(0)
  end if
  
  ! Read observation file
  if (trim(para%get_disper_in()) /= "") then
     obs_disp = observation_disper(para%get_disper_in(), verb)
  else
     if (verb) write(*,*)"<< No dispersion curve input >>"
     call obs_disp%set_n_disp(0)
  end if
  
  ! Covariance matrix
  allocate(cov(obs_rf%get_n_rf()))
  do i = 1, obs_rf%get_n_rf()
     cov(i) = covariance(                    &
          & n       = obs_rf%get_n_smp(i),   &
          & a_gauss = obs_rf%get_a_gauss(i), &
          & delta   = obs_rf%get_delta(i),   &
          & verb  = verb                     &
          & )
  end do


  ! Set interpreter 
  intpr = interpreter(&
       & nlay_max= para%get_k_max(),                             &
       & z_min = para%get_z_min(),                               &
       & z_max = para%get_z_max(),     &
       & n_bin_z = para%get_n_bin_z(),                           &
       & vs_min = para%get_vs_min(), &
       & vs_max = para%get_vs_max(), &
       & n_disp = obs_disp%get_n_disp(),                         & 
       & n_rf = obs_rf%get_n_rf(),      &
       & n_bin_sig = para%get_n_bin_sig(),                       &
       & n_bin_vs = para%get_n_bin_vs(),                         &
       & vp_min = para%get_vp_min(), &
       & vp_max = para%get_vp_max(), &
       & n_bin_vp = para%get_n_bin_vp(),                         &
       & is_ocean = para%get_is_ocean(),                         &
       & ocean_thick = para%get_ocean_thick(),                   &
       & solve_vp = para%get_solve_vp(),                         &
       & solve_anomaly = para%get_solve_anomaly(),               &
       & dvs_sig = para%get_dvs_sig(),                           &
       & dvp_sig = para%get_dvp_sig(),                           &
       & ref_vmod_in = para%get_ref_vmod_in(),                   &
       & vp_bottom = para%get_vp_bottom(),                       &
       & vs_bottom = para%get_vs_bottom(),                       &
       & rho_bottom = para%get_rho_bottom(),                     &
       & verb = verb                                             &
       & )

  ! Set proposal
  prop = proposal( &
       & solve_vp         = para%get_solve_vp(),         &
       & solve_rf_sig     = para%get_solve_rf_sig(),     &
       & solve_disper_sig = para%get_solve_disper_sig(), &
       & n_rf             = obs_rf%get_n_rf(),           &
       & n_disp           = obs_disp%get_n_disp(),       &
       & verb             = verb                         &
       &) 

  ! Set transdimensional model prameters
  if (para%get_solve_vp()) then
     n_td = 3 ! Vs, Vp, z
  else 
     n_td = 2  ! Vs, z
  end if
  tm = trans_d_model( &
       & k_min = para%get_k_min(), &
       & k_max = para%get_k_max(), &
       & nx    = n_td,               &
       & verb  = verb              &
       & )
  call tm%set_prior(prop%get_i_vs(), id_uni, &
       & para%get_vs_min(), para%get_vs_max())
  call tm%set_prior(prop%get_i_depth(),  id_uni, &
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
  

  ! Set hyper parameters
  if (obs_rf%get_n_rf() > 0) then
     hyp_rf = hyper_model(             &
          & nx    = obs_rf%get_n_rf(), &
          & verb  = verb               &
          & )
     if (para%get_solve_rf_sig()) then
        do i = 1, obs_rf%get_n_rf()
           call hyp_rf%set_prior(i, &
                & obs_rf%get_sigma_min(i), obs_rf%get_sigma_max(i))
           call hyp_rf%set_perturb(i, obs_rf%get_dev_sigma(i))
        end do
     else 
        do i = 1, obs_rf%get_n_rf()
           call hyp_rf%set_prior(i, &
                & obs_rf%get_sigma_min(i), obs_rf%get_sigma_min(i))
        end do
     end if
  end if

  if (obs_disp%get_n_disp() > 0) then
     hyp_disp = hyper_model(                   &
          & nx    = obs_disp%get_n_disp() * 2, &
          & verb  = verb                       &
          & )
     if (para%get_solve_disper_sig()) then
        do i = 1, obs_disp%get_n_disp()
           call hyp_disp%set_prior(2*i - 1, &
                & obs_disp%get_sig_c_min(i), obs_disp%get_sig_c_max(i))
           call hyp_disp%set_perturb(2*i - 1, obs_disp%get_dev_sig_c(i))
           call hyp_disp%set_prior(2*i, &
                & obs_disp%get_sig_u_min(i), obs_disp%get_sig_u_max(i))
           call hyp_disp%set_perturb(2*i, obs_disp%get_dev_sig_u(i))
        end do
     else
        do i = 1, obs_disp%get_n_disp()
           call hyp_disp%set_prior(2*i - 1, &
                & obs_disp%get_sig_c_min(i), obs_disp%get_sig_c_min(i))
           call hyp_disp%set_prior(2*i, &
                & obs_disp%get_sig_u_min(i), obs_disp%get_sig_u_min(i))
           !call hyp_disp%set_x(2*i - 1, obs_disp%get_sig_c_min(i))
           !call hyp_disp%set_x(2*i    , obs_disp%get_sig_u_min(i))
        end do
     end if
  end if
  
  ! Set forward computation (recv_func)
  if (obs_rf%get_n_rf() > 0) then
     call tm%generate_model()
     call intpr%construct_vmodel(tm, vm, is_ok)
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
  end if
  
  ! Set forward computation (disper)
  if (obs_disp%get_n_disp() > 0) then
     call tm%generate_model()
     call intpr%construct_vmodel(tm, vm, is_ok)
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
  end if
    
  ! Initialize parallel chains
  pt = parallel(                       &
       & n_proc  = n_proc,             &
       & rank    = rank,               &
       & n_chain = para%get_n_chain(), &
       & verb    = verb)
  
  ! Set MCMC chain
  do i = 1, para%get_n_chain()
     ! Generate random model
     call tm%generate_model()
     if (obs_rf%get_n_rf() > 0) then
        call hyp_rf%generate_model()
     end if
     if (obs_disp%get_n_disp() > 0) then
        call hyp_disp%generate_model()
     end if


     ! Set transdimensional model
     mc = mcmc( &
          & tm       = tm,     &
          & hyp_disp = hyp_disp,         &
          & hyp_rf   = hyp_rf,           &
          & prop     = prop,             &
          & n_iter   = para%get_n_iter() &
          & )
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
  

  !---------------------------------------------------------------------
  ! 2. MCMC 
  !---------------------------------------------------------------------
  write(filename,'("vmod_out.",I3.3)') rank
  open(newunit = io_vmod_all, file = filename, iostat = ierr, &
       & status = "replace")
  if (ierr /=0) then
     write(0,*)"ERROR: cannot create: ", trim(filename)
     call mpi_finalize(ierr)
     stop
  end if
  do i = 1, para%get_n_iter()
     ! MCMC step
     do j = 1, para%get_n_chain()
        
        mc = pt%get_mc(j)

        ! Proposal
        call mc%propose_model(tm_tmp, hyp_disp_tmp, hyp_rf_tmp, &
             & is_ok, log_prior_ratio, &
             & log_proposal_ratio)


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !call intpr%construct_vmodel(tm_tmp, vm, is_ok2)
        !write(io_vmod_all,'("> ",E15.7)') mc%get_log_likelihood()
        !call vm%display(io_vmod_all)
        !call mpi_barrier(MPI_COMM_WORLD ,ierr)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        ! Forward computation
        if (is_ok) then
           log_likelihood = 0.d0
           log_likelihood2 = 0.d0
           if (obs_rf%get_n_rf() > 0) then
              call forward_recv_func(tm_tmp, hyp_rf_tmp, intpr, obs_rf, &
                   & rf, cov, log_likelihood, is_ok2)
           end if
           if (obs_disp%get_n_disp() > 0) then
              call forward_disper(tm_tmp, hyp_disp_tmp, intpr, obs_disp, &
                   & disp, log_likelihood2, is_ok2)
           end if
           log_likelihood = &
                log_likelihood + log_likelihood2
        else
           log_likelihood = minus_infty
        end if
        

        ! Judege
        call mc%judge_model(tm_tmp, hyp_disp_tmp, hyp_rf_tmp, is_ok, &
             & log_likelihood, log_prior_ratio, log_proposal_ratio)
        call pt%set_mc(j, mc)
        
        ! One step summary
        if (verb) then
           call mc%one_step_summary()
        end if
        
        ! Recording
        if (i > para%get_n_burn() .and. &
             & mod(i, para%get_n_corr()) == 0 .and. &
             & mc%get_temp() < 1.d0 + eps) then
           
           ! V-Z
           call intpr%save_model(mc%get_tm())

           ! Sigma
           if (para%get_solve_disper_sig()) then
              call intpr%save_disp_sigma(mc%get_hyp_disp())
           end if
           if (para%get_solve_rf_sig()) then
              call intpr%save_rf_sigma(mc%get_hyp_rf())
           end if

           ! Synthetic data 
           do k = 1, obs_rf%get_n_rf()
              call rf(k)%save_syn()
           end do
           do k = 1, obs_disp%get_n_disp()
              call disp(k)%save_syn()
           end do
           
           ! V model
           !call intpr%construct_vmodel(mc%get_tm(), vm, is_ok)
           !write(io_vmod_all,'("> ",E15.7)') mc%get_log_likelihood()
           !call vm%display(io_vmod_all)
           
        end if
     end do
     
     ! Swap temperture
     call pt%swap_temperature(verb=.true.)
  end do
  close(io_vmod_all)

  !---------------------------------------------------------------------
  ! 3. Output
  !---------------------------------------------------------------------

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
  
  ! RF
  do i = 1, obs_rf%get_n_rf()
     ! Synthetic 
     del_amp = (rf(i)%get_amp_max() - rf(i)%get_amp_min()) / &
          & rf(i)%get_n_bin_amp()
     write(filename, '(A6,I3.3,A4)')"syn_rf", i, ".ppd"
     call output_ppd_2d(filename, rank, obs_rf%get_n_smp(i) * 2, &
          & rf(i)%get_n_bin_amp(), rf(i)%get_n_syn_rf(), &
          & n_mod, obs_rf%get_t_start(i), obs_rf%get_delta(i), &
          & rf(i)%get_amp_min(), del_amp)

     ! Noise
     del_amp = (hyp_rf%get_prior_param(i, 2) - &
          & hyp_rf%get_prior_param(i, 1)) / para%get_n_bin_sig()
     write(filename, '(A,I3.3,A)')"rf_sigma", i, ".ppd"
     call output_ppd_1d(filename, rank, para%get_n_bin_sig(), &
          & intpr%get_n_rf_sig(i), n_mod, &
          & hyp_rf%get_prior_param(i, 1), del_amp)
  end do
  
  ! Dispersion
  do i = 1, obs_disp%get_n_disp()
     ! Synthetic phase velocity
     write(filename, '(A9,I3.3,A4)')"syn_phase", i, ".ppd"
     call output_ppd_2d(filename, rank, disp(i)%get_nf(), &
          & disp(i)%get_nc(), disp(i)%get_n_fc(), &
          & n_mod, obs_disp%get_fmin(i), obs_disp%get_df(i), &
          & obs_disp%get_cmin(i), obs_disp%get_dc(i))
     
     ! Synthetic group velocity
     write(filename, '(A9,I3.3,A4)')"syn_group", i, ".ppd"
     call output_ppd_2d(filename, rank, disp(i)%get_nf(), &
          & disp(i)%get_nc(), disp(i)%get_n_fu(), &
          & n_mod, obs_disp%get_fmin(i), obs_disp%get_df(i), &
          & obs_disp%get_cmin(i), obs_disp%get_dc(i))

     ! Noise phase velocity
     del_amp = (hyp_disp%get_prior_param(2*i-1, 2) - &
          & hyp_disp%get_prior_param(2*i-1, 1)) / para%get_n_bin_sig()
     write(filename, '(A,I3.3,A)')"phase_sigma", i, ".ppd"
     call output_ppd_1d(filename, rank, para%get_n_bin_sig(), &
          & intpr%get_n_disp_sig(2*i-1), n_mod, &
          & hyp_disp%get_prior_param(2*i-1, 1), del_amp)

     ! Noise group velocity
     del_amp = (hyp_disp%get_prior_param(2*i, 2) - &
          & hyp_disp%get_prior_param(2*i, 1)) / para%get_n_bin_sig()
     write(filename, '(A,I3.3,A)')"group_sigma", i, ".ppd"
     call output_ppd_1d(filename, rank, para%get_n_bin_sig(), &
          & intpr%get_n_disp_sig(2*i), n_mod, &
          & hyp_disp%get_prior_param(2*i, 1), del_amp)
  end do
  
  ! Likelihood history
  filename = "likelihood.history"
  call pt%output_history(filename, 'l')
  filename = "temp.history"
  call pt%output_history(filename, 't')
  
  ! Proposal count
  filename = "proposal.count"
  call pt%output_proposal(filename, prop%get_labels())
  
  call mpi_finalize(ierr)
  stop
end program main
