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
module cls_param
  use cls_line_text
  implicit none 
  
  type param
     private
     character(len=line_max) :: param_file
     character(len=line_max), allocatable :: given(:)
     
     ! Iteration numbers
     integer :: n_iter = -12345
     integer :: n_burn = -12345
     integer :: n_corr = -12345

     ! Parallel tempering
     double precision :: temp_high = -12345.d0
     integer          :: n_chain   = -12345
     integer          :: n_cool  = -12345
     
     !Prior probability
     integer          :: k_min = -12345
     integer          :: k_max = -12345
     double precision :: z_min = -12345.d0
     double precision :: z_max = -12345.d0
     double precision :: vs_min = -12345.d0
     double precision :: vs_max = -12345.d0
     double precision :: vp_min = -12345.d0
     double precision :: vp_max = -12345.d0 

     ! MCMC proposal
     integer :: i_seed1 = 11111111
     integer :: i_seed2 = 22222222
     integer :: i_seed3 = 33333333
     integer :: i_seed4 = 44444444

     
     
     logical :: solve_vp = .false.
     logical :: solve_anomaly = .false.
     logical :: solve_rf_sig = .true.
     logical :: solve_disper_sig = .true.
     logical :: is_sphere = .false.
     double precision :: r_earth = 6371.d0

     double precision :: dvp_sig = 0.1d0
     double precision :: dvs_sig = 0.1d0
     double precision :: rf_sig_min = 0.005d0
     double precision :: rf_sig_max = 0.1d0
     double precision :: disper_sig_min = 0.01d0
     double precision :: disper_sig_max = 0.2d0
     integer :: n_bin_z  = 100
     integer :: n_bin_vs = 100
     integer :: n_bin_vp = 100
     integer :: n_bin_sig = 100
     
     ! Parameters for random perturbation
     double precision :: dev_z  = 0.1d0
     double precision :: dev_vp = 0.04d0
     double precision :: dev_vs = 0.02d0
     double precision :: dev_rf_sig = 0.02d0
     double precision :: dev_disper_sig = 0.02d0
     
     ! Parameters for ocean layer
     logical :: is_ocean = .false.
     double precision :: ocean_thick = 0.d0
     
     ! Attenuation
     logical :: is_attenuative = .false.

     ! Parameters for receiver functions
     ! only used by recv_func_fwd
     integer          :: n_smp    = -12345
     double precision :: rayp     = -12345.d0
     double precision :: a_gauss  = -12345.d0
     double precision :: delta    = -12345.d0
     double precision :: t_pre    = -12345.d0
     character(1)     :: rf_phase = ""
     double precision :: damp     = -12345.d0 ! damping factor for water level
     
     double precision :: amp_min = -0.6d0
     double precision :: amp_max = 0.6d0

     
     logical :: deconv_flag
     logical :: correct_amp
     
     integer :: n_bin_c  = 50

     ! Parameter for dispersion
     ! only used by rayleigh_fwd 
     character(6)     :: freq_or_period = ""
     double precision :: xmin           = -12345.d0
     double precision :: xmax           = -12345.d0
     double precision :: dx             = -12345.d0
     double precision :: cmin           = -12345.d0
     double precision :: cmax           = -12345.d0
     double precision :: dc             = -12345.d0
     character(1)     :: disper_phase   = ""
     integer          :: n_mode         = -12345

     

     ! Noise level added to forward computation 
     double precision :: noise_added = 0.d0
     
     ! I/O file
     character(len=line_max) :: vmod_in = ""
     character(len=line_max) :: ref_vmod_in = ""
     character(len=line_max) :: disper_out = ""
     character(len=line_max) :: disper_in = ""
     character(len=line_max) :: recv_func_in = ""
     character(len=line_max) :: recv_func_out = ""

     ! bottom layer
     double precision :: vs_bottom = 4.6d0
     double precision :: vp_bottom = 8.1d0
     double precision :: rho_bottom = 3.3d0
     
     logical :: diagnostic_mode = .false.
     logical :: verb = .false.


   contains
     procedure :: read_file => param_read_file
     procedure :: set_value  => param_set_value
     procedure :: get_n_iter => param_get_n_iter
     procedure :: get_n_burn => param_get_n_burn
     procedure :: get_n_corr => param_get_n_corr
     procedure :: get_n_chain => param_get_n_chain
     procedure :: get_n_cool =>  param_get_n_cool
     procedure :: get_i_seed1 => param_get_i_seed1
     procedure :: get_i_seed2 => param_get_i_seed2
     procedure :: get_i_seed3 => param_get_i_seed3
     procedure :: get_i_seed4 => param_get_i_seed4
     procedure :: get_k_min => param_get_k_min
     procedure :: get_k_max => param_get_k_max
     procedure :: get_n_bin_z => param_get_n_bin_z
     procedure :: get_n_bin_vs => param_get_n_bin_vs
     procedure :: get_n_bin_vp => param_get_n_bin_vp
     procedure :: get_n_bin_sig => param_get_n_bin_sig
     procedure :: get_n_bin_c => param_get_n_bin_c
     
     procedure :: get_temp_high => param_get_temp_high
     procedure :: get_freq_or_period => param_get_freq_or_period
     procedure :: get_xmin => param_get_xmin
     procedure :: get_xmax => param_get_xmax
     procedure :: get_dx => param_get_dx
     procedure :: get_cmin => param_get_cmin
     procedure :: get_cmax => param_get_cmax
     procedure :: get_dc => param_get_dc
     procedure :: get_n_mode => param_get_n_mode
     procedure :: get_disper_phase => param_get_disper_phase
     procedure :: get_dev_vs => param_get_dev_vs
     procedure :: get_dev_vp => param_get_dev_vp
     procedure :: get_dev_z => param_get_dev_z
     procedure :: get_dev_rf_sig => param_get_dev_rf_sig
     procedure :: get_dev_disper_sig => param_get_dev_disper_sig
     procedure :: get_vs_min => param_get_vs_min
     procedure :: get_vs_max => param_get_vs_max
     procedure :: get_vp_min => param_get_vp_min
     procedure :: get_vp_max => param_get_vp_max
     procedure :: get_dvs_sig => param_get_dvs_sig
     procedure :: get_dvp_sig => param_get_dvp_sig
     procedure :: get_z_min => param_get_z_min
     procedure :: get_z_max => param_get_z_max
     procedure :: get_rf_sig_min => param_get_rf_sig_min
     procedure :: get_rf_sig_max => param_get_rf_sig_max
     procedure :: get_disper_sig_min => param_get_disper_sig_min
     procedure :: get_disper_sig_max => param_get_disper_sig_max

     procedure :: get_ocean_thick => param_get_ocean_thick

     procedure :: get_solve_vp => param_get_solve_vp
     procedure :: get_solve_anomaly => param_get_solve_anomaly
     procedure :: get_solve_rf_sig => param_get_solve_rf_sig
     procedure :: get_solve_disper_sig => param_get_solve_disper_sig
     procedure :: get_is_ocean => param_get_is_ocean

     procedure :: get_is_attenuative => param_get_is_attenuative
     
     procedure :: get_ref_vmod_in => param_get_ref_vmod_in
     procedure :: get_vmod_in => param_get_vmod_in
     procedure :: get_disper_out => param_get_disper_out
     procedure :: get_disper_in  => param_get_disper_in
     procedure :: get_recv_func_in => param_get_recv_func_in

     procedure :: get_n_smp => param_get_n_smp
     procedure :: get_rayp => param_get_rayp
     procedure :: get_a_gauss => param_get_a_gauss
     procedure :: get_t_pre => param_get_t_pre
     procedure :: get_rf_phase => param_get_rf_phase
     procedure :: get_delta => param_get_delta
     procedure :: get_amp_min => param_get_amp_min
     procedure :: get_amp_max => param_get_amp_max
     procedure :: get_damp => param_get_damp
     procedure :: get_deconv_flag => param_get_deconv_flag
     procedure :: get_correct_amp => param_get_correct_amp
     procedure :: get_recv_func_out => param_get_recv_func_out
     
     procedure :: get_vp_bottom => param_get_vp_bottom
     procedure :: get_vs_bottom => param_get_vs_bottom
     procedure :: get_rho_bottom => param_get_rho_bottom
     
     procedure :: get_noise_added => param_get_noise_added

     procedure :: get_is_sphere => param_get_is_sphere
     procedure :: get_r_earth => param_get_r_earth
     
     procedure :: get_diagnostic_mode => param_get_diagnostic_mode

     procedure :: check_mcmc_params => param_check_mcmc_params
     !procedure :: check_recv_func_fwd_params &
     !     & => param_check_recv_func_fwd_params
     !procedure :: check_disper_fwd_params &
     !     & => param_check_disper_fwd_params
     
     !procedure :: is_given_double  => param_is_given_double
     !procedure :: is_given_integer => param_is_given_integer
     !procedure :: is_given_char    => param_is_given_char
     !generic :: is_given => is_given_double, is_given_integer, &
     !     & is_given_char
     procedure :: check_given_param => param_check_given_param

  end type param
  
  interface param
     module procedure init_param
  end interface param

contains
  
  !---------------------------------------------------------------------
  
  type(param) function init_param(param_file, verb) result(self)
    character(len=*), intent(in) :: param_file
    logical, intent(in), optional :: verb

    if (present(verb)) then
       self%verb = verb
    end if
    
    if (self%verb) then
       write(*,'(3A)')"<< Reading parameters from ", &
            & trim(param_file), " >>"
    end if    
    
    self%param_file = param_file
    self%given = [ character(line_max) :: ]
    call self%read_file()

    if (self%verb) then
       write(*,*)
    end if

    return 
  end function init_param

  !---------------------------------------------------------------------
    
  subroutine param_read_file(self)
    class(param), intent(inout) :: self
    character(len=line_max) :: line, name, val
    integer :: ierr, io
    type(line_text) :: lt
    logical :: is_ok

    open(newunit = io, file = self%param_file, &
         & status = 'old', iostat = ierr)
    if (ierr /= 0) then
       if (self%verb) then
          write(0,*) "ERROR: cannot open ", trim(self%param_file)
       end if
       call mpi_finalize(ierr)
       stop
    end if
    
    do 
       read(io, '(a)', iostat=ierr) line
       if (ierr /= 0) then
          exit
       end if
       lt = init_line_text(line)
       call lt%read_value(name, val, is_ok)
       if (is_ok) then
          call self%set_value(name, val)
       end if
    end do
    close(io)
    

    

    return 
  end subroutine param_read_file

  !---------------------------------------------------------------------
  
  subroutine param_set_value(self, name, val)
    class(param), intent(inout) :: self
    character(len=*), intent(in) :: name, val
    integer :: ierr
    character(line_max), allocatable :: tmp(:)
    
    if (self%verb) then
       write(*,*)trim(name), " <- ", trim(val)
    end if
    if (name == "xmin") then
       read(val, *) self%xmin
    else if (name == "xmax") then
       read(val, *) self%xmax 
    else if (name == "dx") then
       read(val, *) self%dx 
    else if (name == "freq_or_period") then
       self%freq_or_period = val
    else if (name == "cmin") then
       read(val, *) self%cmin 
    else if (name == "cmax") then
       read(val, *) self%cmax 
    else if (name == "dc") then
       read(val, *) self%dc 
    else if (name == "n_mode") then
       read(val, *) self%n_mode
    else if (name == "disper_phase") then
       self%disper_phase = val
    else if (name == "vmod_in") then
       self%vmod_in = val
    else if (name == "disper_out") then
       self%disper_out = val
    else if (name == "n_iter") then
       read(val, *) self%n_iter
    else if (name == "n_burn") then
       read(val, *) self%n_burn 
    else if (name == "n_corr") then
       read(val, *) self%n_corr 
    else if (name == "n_chain") then
       read(val, *) self%n_chain
    else if (name == "n_cool") then
       read(val, *) self%n_cool 
    else if (name == "temp_high") then
       read(val, *) self%temp_high
    else if (name == "i_seed1") then
       read(val, *) self%i_seed1 
    else if (name == "i_seed2") then
       read(val, *) self%i_seed2 
    else if (name == "i_seed3") then
       read(val, *) self%i_seed3
    else if (name == "i_seed4") then
       read(val, *) self%i_seed4
    else if (name == "k_min") then
       read(val, *) self%k_min 
    else if (name == "k_max") then
       read(val, *) self%k_max
    else if (name == "disper_in") then
       self%disper_in = val
    else if (name == "dev_vs") then
       read(val, *) self%dev_vs
    else if (name == "dev_vp") then
       read(val, *) self%dev_vp
    else if (name == "dev_z") then
       read(val, *) self%dev_z 
    else if (name == "dev_rf_sig") then
       read(val, *) self%dev_rf_sig
    else if (name == "dev_disper_sig") then
       read(val, *) self%dev_disper_sig
    else if (name == "vs_min") then
       read(val, *) self%vs_min
    else if (name == "vs_max") then
       read(val, *) self%vs_max
    else if (name == "vp_min") then
       read(val, *) self%vp_min
    else if (name == "vp_max") then
       read(val, *) self%vp_max
    else if (name == "z_min") then
       read(val, *) self%z_min 
    else if (name == "z_max") then
       read(val, *) self%z_max 
    else if (name == "rf_sig_min") then
       read(val, *) self%rf_sig_min
    else if (name == "rf_sig_max") then
       read(val, *) self%rf_sig_max
    else if (name == "disper_sig_min") then
       read(val, *) self%disper_sig_min
    else if (name == "disper_sig_max") then
       read(val, *) self%disper_sig_max
    else if (name == "dvs_sig") then
       read(val, *) self%dvs_sig
    else if (name == "dvp_sig") then
       read(val, *) self%dvp_sig
    else if (name == "ocean_thick") then
       read(val, *) self%ocean_thick 
    else if (name == "solve_vp") then
       read(val, *) self%solve_vp
    else if (name == "solve_anomaly") then
       read(val, *) self%solve_anomaly
    else if (name == "solve_rf_sig") then
       read(val, *)self%solve_rf_sig
    else if (name == "solve_disper_sig") then
       read(val, *)self%solve_disper_sig
    else if (name == "is_sphere") then
       read(val, *) self%is_sphere
    else if (name == "r_earth") then
       read(val, *) self%r_earth
    else if (name == "is_ocean") then
       read(val, *) self%is_ocean 
    else if (name == "is_attenuative") then
       read(val, *) self%is_attenuative
    else if (name == "n_bin_z") then
       read(val, *) self%n_bin_z 
    else if (name == "n_bin_vs") then
       read(val, *) self%n_bin_vs
    else if (name == "n_bin_vp") then
       read(val, *) self%n_bin_vp
    else if (name == "n_bin_sig") then
       read(val, *) self%n_bin_sig
    else if (name == "n_bin_c") then
       read(val, *) self%n_bin_c
    else if (name == "n_smp") then
       read(val, *) self%n_smp
    else if (name == "delta") then
       read(val, *) self%delta
    else if (name == "amp_min") then
       read(val, *) self%amp_min
    else if (name == "amp_max") then
       read(val, *) self%amp_max
    else if (name == "damp") then
       read(val, *) self%damp
    else if (name == "rayp") then
       read(val, *) self%rayp
    else if (name == "a_gauss") then
       read(val, *) self%a_gauss
    else if (name == "t_pre") then
       read(val, *) self%t_pre
    else if (name == "rf_phase") then
       self%rf_phase = val
    else if (name == "deconv_flag") then
       read(val, *) self%deconv_flag
    else if (name == "correct_amp") then
       read(val, *) self%correct_amp
    else if (name == "recv_func_out") then
       self%recv_func_out = val
    else if (name == "recv_func_in") then
       self%recv_func_in = val
    else if (name == "vs_bottom") then
       read(val, *) self%vs_bottom
    else if (name == "vp_bottom") then
       read(val, *) self%vp_bottom
    else if (name == "rho_bottom") then
       read(val, *) self%rho_bottom
    else if (name == "ref_vmod_in") then
       self%ref_vmod_in = val
    else if (name == "noise_added") then
       read(val, *) self%noise_added
    else if (name == "diagnostic_mode") then
       read(val, *) self%diagnostic_mode
    else
       if (self%verb) then
          write(0,*)"ERROR: Invalid parameter name"
          write(0,*)"        : ", name, "  (?)"
       end if
       call mpi_finalize(ierr)
       stop
    end if
    
    tmp = name
    call move_alloc(self%given, tmp)
    allocate(self%given(size(tmp)+1))
    self%given(1:size(tmp)) = tmp(:)
    self%given(size(tmp)+1) = name
    
    return 
  end subroutine param_set_value

  !---------------------------------------------------------------------

  integer function param_get_n_iter(self) result(n_iter)
    class(param), intent(in) :: self
    
    n_iter = self%n_iter
    
    return 
  end function param_get_n_iter

  !---------------------------------------------------------------------
  
  integer function param_get_n_burn(self) result(n_burn)
    class(param), intent(in) :: self
    
    n_burn = self%n_burn
    
    return 
  end function param_get_n_burn

  !---------------------------------------------------------------------

  integer function param_get_n_corr(self) result(n_corr)
    class(param), intent(in) :: self
    
    n_corr = self%n_corr
    
    return 
  end function param_get_n_corr

  
  !---------------------------------------------------------------------

  integer function param_get_n_chain(self) result(n_chain)
    class(param), intent(in) :: self
    
    n_chain = self%n_chain
    
    return 
  end function param_get_n_chain
  
  !---------------------------------------------------------------------

  integer function param_get_n_cool(self) result(n_cool)
    class(param), intent(in) :: self
    
    n_cool = self%n_cool
    
    return 
  end function param_get_n_cool
  
  !---------------------------------------------------------------------

  integer function param_get_i_seed1(self) result(i_seed1)
    class(param), intent(in) :: self
    
    i_seed1 = self%i_seed1
    
    return 
  end function param_get_i_seed1
  
  !---------------------------------------------------------------------

  integer function param_get_i_seed2(self) result(i_seed2)
    class(param), intent(in) :: self
    
    i_seed2 = self%i_seed2
    
    return 
  end function param_get_i_seed2
  
  !---------------------------------------------------------------------

  integer function param_get_i_seed3(self) result(i_seed3)
    class(param), intent(in) :: self
    
    i_seed3 = self%i_seed3
    
    return 
  end function param_get_i_seed3
  
  !---------------------------------------------------------------------

  integer function param_get_i_seed4(self) result(i_seed4)
    class(param), intent(in) :: self
    
    i_seed4 = self%i_seed4
    
    return 
  end function param_get_i_seed4
  
  !---------------------------------------------------------------------

  integer function param_get_k_min(self) result(k_min)
    class(param), intent(in) :: self
    
    k_min = self%k_min
    
    return 
  end function param_get_k_min
  
  !---------------------------------------------------------------------

  integer function param_get_k_max(self) result(k_max)
    class(param), intent(in) :: self
    
    k_max = self%k_max
    
    return 
  end function param_get_k_max
  
  !---------------------------------------------------------------------

  integer function param_get_n_bin_z(self) result(n_bin_z)
    class(param), intent(in) :: self
    
    n_bin_z = self%n_bin_z
    
    return 
  end function param_get_n_bin_z
  
  
  !---------------------------------------------------------------------

  integer function param_get_n_bin_vs(self) result(n_bin_vs)
    class(param), intent(in) :: self
    
    n_bin_vs = self%n_bin_vs
    
    return 
  end function param_get_n_bin_vs
  
  !---------------------------------------------------------------------

  integer function param_get_n_bin_vp(self) result(n_bin_vp)
    class(param), intent(in) :: self
    
    n_bin_vp = self%n_bin_vp
    
    return 
  end function param_get_n_bin_vp

  !---------------------------------------------------------------------

  integer function param_get_n_bin_sig(self) result(n_bin_sig)
    class(param), intent(in) :: self
    
    n_bin_sig = self%n_bin_sig
    
    return 
  end function param_get_n_bin_sig

  !---------------------------------------------------------------------

  integer function param_get_n_bin_c(self) result(n_bin_c)
    class(param), intent(in) :: self
    
    n_bin_c = self%n_bin_c
    
    return 
  end function param_get_n_bin_c

  !---------------------------------------------------------------------

  
  double precision function param_get_temp_high(self) result(temp_high)
    class(param), intent(in) :: self

    temp_high = self%temp_high

    return
  end function param_get_temp_high
  
  !---------------------------------------------------------------------
  
  character(6) function param_get_freq_or_period(self) result(freq_or_period)
    class(param), intent(in) :: self
    
    freq_or_period = self%freq_or_period
    
    return 
  end function param_get_freq_or_period

  !---------------------------------------------------------------------

  double precision function param_get_xmin(self) result(xmin)
    class(param), intent(in) :: self

    xmin = self%xmin

    return
  end function param_get_xmin
  
  !---------------------------------------------------------------------

  double precision function param_get_xmax(self) result(xmax)
    class(param), intent(in) :: self

    xmax = self%xmax

    return
  end function param_get_xmax
  !---------------------------------------------------------------------

  double precision function param_get_dx(self) result(dx)
    class(param), intent(in) :: self

    dx = self%dx

    return
  end function param_get_dx
  
  !---------------------------------------------------------------------

  double precision function param_get_cmin(self) result(cmin)
    class(param), intent(in) :: self

    cmin = self%cmin

    return
  end function param_get_cmin
  
  !---------------------------------------------------------------------

  double precision function param_get_cmax(self) result(cmax)
    class(param), intent(in) :: self

    cmax = self%cmax

    return
  end function param_get_cmax
  !---------------------------------------------------------------------

  double precision function param_get_dc(self) result(dc)
    class(param), intent(in) :: self

    dc = self%dc

    return
  end function param_get_dc

  !---------------------------------------------------------------------
  
  integer function param_get_n_mode(self) result(n_mode)
    class(param), intent(in) :: self
    
    n_mode = self%n_mode
    
    return 
  end function param_get_n_mode

  !---------------------------------------------------------------------
  
  character(1) function param_get_disper_phase(self) &
       & result(disper_phase)
    class(param), intent(in) :: self
    
    disper_phase = self%disper_phase
    
    return 
  end function param_get_disper_phase

  !---------------------------------------------------------------------

  double precision function param_get_dev_vs(self) result(dev_vs)
    class(param), intent(in) :: self

    dev_vs = self%dev_vs

    return
  end function param_get_dev_vs

  !---------------------------------------------------------------------

  double precision function param_get_dev_vp(self) result(dev_vp)
    class(param), intent(in) :: self

    dev_vp = self%dev_vp

    return
  end function param_get_dev_vp

  !---------------------------------------------------------------------

  double precision function param_get_dev_z(self) result(dev_z)
    class(param), intent(in) :: self

    dev_z = self%dev_z

    return
  end function param_get_dev_z

  !---------------------------------------------------------------------

  double precision function param_get_dev_rf_sig(self) &
       & result(dev_rf_sig)
    class(param), intent(in) :: self

    dev_rf_sig = self%dev_rf_sig

    return
  end function param_get_dev_rf_sig

  !---------------------------------------------------------------------

  double precision function param_get_dev_disper_sig(self) &
       & result(dev_disper_sig)
    class(param), intent(in) :: self

    dev_disper_sig = self%dev_disper_sig

    return
  end function param_get_dev_disper_sig

  !---------------------------------------------------------------------
  
  double precision function param_get_vs_min(self) result(vs_min)
    class(param), intent(in) :: self

    vs_min = self%vs_min

    return
  end function param_get_vs_min

  !---------------------------------------------------------------------

  double precision function param_get_vs_max(self) result(vs_max)
    class(param), intent(in) :: self

    vs_max = self%vs_max

    return
  end function param_get_vs_max

  !---------------------------------------------------------------------
  
  double precision function param_get_vp_min(self) result(vp_min)
    class(param), intent(in) :: self

    vp_min = self%vp_min

    return
  end function param_get_vp_min

  !---------------------------------------------------------------------

  double precision function param_get_vp_max(self) result(vp_max)
    class(param), intent(in) :: self

    vp_max = self%vp_max

    return
  end function param_get_vp_max

  !---------------------------------------------------------------------

  double precision function param_get_dvs_sig(self) result(dvs_sig)
    class(param), intent(in) :: self

    dvs_sig = self%dvs_sig

    return
  end function param_get_dvs_sig

  !---------------------------------------------------------------------

  double precision function param_get_dvp_sig(self) result(dvp_sig)
    class(param), intent(in) :: self

    dvp_sig = self%dvp_sig

    return
  end function param_get_dvp_sig

  !---------------------------------------------------------------------
  
  double precision function param_get_z_min(self) result(z_min)
    class(param), intent(in) :: self

    z_min = self%z_min

    return
  end function param_get_z_min

  !---------------------------------------------------------------------

  double precision function param_get_z_max(self) result(z_max)
    class(param), intent(in) :: self

    z_max = self%z_max

    return
  end function param_get_z_max

  !---------------------------------------------------------------------

  double precision function param_get_rf_sig_min(self) result(rf_sig_min)
    class(param), intent(in) :: self

    rf_sig_min = self%rf_sig_min

    return
  end function param_get_rf_sig_min

  !---------------------------------------------------------------------

  double precision function param_get_rf_sig_max(self) result(rf_sig_max)
    class(param), intent(in) :: self

    rf_sig_max = self%rf_sig_max

    return
  end function param_get_rf_sig_max

  !---------------------------------------------------------------------

  double precision function param_get_disper_sig_min(self) &
       &result(disper_sig_min)
    class(param), intent(in) :: self

    disper_sig_min = self%disper_sig_min

    return
  end function param_get_disper_sig_min

  !---------------------------------------------------------------------

  double precision function param_get_disper_sig_max(self) &
       & result(disper_sig_max)
    class(param), intent(in) :: self

    disper_sig_max = self%disper_sig_max

    return
  end function param_get_disper_sig_max

  !---------------------------------------------------------------------

  double precision function param_get_ocean_thick(self) &
       & result(ocean_thick)
    class(param), intent(in) :: self

    ocean_thick = self%ocean_thick

    return
  end function param_get_ocean_thick

  !---------------------------------------------------------------------
  
  logical function param_get_solve_anomaly(self) result(solve_anomaly)
    class(param), intent(in) :: self

    solve_anomaly = self%solve_anomaly
    
    return 
  end function param_get_solve_anomaly

  !---------------------------------------------------------------------

  logical function param_get_solve_vp(self) result(solve_vp)
    class(param), intent(in) :: self

    solve_vp = self%solve_vp
    
    return 
  end function param_get_solve_vp
  
  !---------------------------------------------------------------------

  logical function param_get_solve_rf_sig(self) result(solve_rf_sig)
    class(param), intent(in) :: self

    solve_rf_sig = self%solve_rf_sig
    
    return 
  end function param_get_solve_rf_sig
  
  !---------------------------------------------------------------------

  logical function param_get_solve_disper_sig(self) result(solve_disper_sig)
    class(param), intent(in) :: self

    solve_disper_sig = self%solve_disper_sig
    
    return 
  end function param_get_solve_disper_sig
  
  !---------------------------------------------------------------------

  logical function param_get_is_ocean(self) result(is_ocean)
    class(param), intent(in) :: self

    is_ocean = self%is_ocean
    
    return 
  end function param_get_is_ocean

  !---------------------------------------------------------------------

  logical function param_get_is_attenuative(self) result(is_attenuative)
    class(param), intent(in) :: self

    is_attenuative = self%is_attenuative
    
    return 
  end function param_get_is_attenuative

  !---------------------------------------------------------------------
  
  character(len=line_max) function param_get_ref_vmod_in(self) &
       & result(ref_vmod_in)
    class(param), intent(in) :: self
    
    ref_vmod_in = self%ref_vmod_in
    
    return
  end function param_get_ref_vmod_in
  
  !---------------------------------------------------------------------
  
  character(len=line_max) function param_get_vmod_in(self) &
       & result(vmod_in)
    class(param), intent(in) :: self
    
    vmod_in = self%vmod_in
    
    return
  end function param_get_vmod_in
  
  !---------------------------------------------------------------------
  
  character(len=line_max) function param_get_disper_out(self) &
       & result(disper_out)
    class(param), intent(in) :: self
    
    disper_out = self%disper_out
    
    return
  end function param_get_disper_out
  
  !---------------------------------------------------------------------

  character(len=line_max) function param_get_disper_in(self) &
       & result(disper_in)
    class(param), intent(in) :: self
    
    disper_in = self%disper_in
    
    return
  end function param_get_disper_in
  
  !---------------------------------------------------------------------

  character(len=line_max) function param_get_recv_func_in(self) &
       & result(recv_func_in)
    class(param), intent(in) :: self
    
    recv_func_in = self%recv_func_in
    
    return
  end function param_get_recv_func_in
  
  !---------------------------------------------------------------------

  integer function param_get_n_smp(self) result(n_smp)
    class(param), intent(in) :: self
    
    n_smp = self%n_smp
    
    return 
  end function param_get_n_smp
  
  !---------------------------------------------------------------------

  double precision function param_get_rayp(self) result(rayp)
    class(param), intent(in) :: self

    rayp = self%rayp

    return
  end function param_get_rayp

  !---------------------------------------------------------------------

  double precision function param_get_a_gauss(self) result(a_gauss)
    class(param), intent(in) :: self

    a_gauss = self%a_gauss

    return
  end function param_get_a_gauss

  !---------------------------------------------------------------------

  double precision function param_get_delta(self) result(delta)
    class(param), intent(in) :: self

    delta = self%delta

    return
  end function param_get_delta

  !---------------------------------------------------------------------

  double precision function param_get_amp_min(self) result(amp_min)
    class(param), intent(in) :: self

    amp_min = self%amp_min

    return
  end function param_get_amp_min

  !---------------------------------------------------------------------

  double precision function param_get_amp_max(self) result(amp_max)
    class(param), intent(in) :: self

    amp_max = self%amp_max
    
    return
  end function param_get_amp_max

  !---------------------------------------------------------------------

  double precision function param_get_damp(self) result(damp)
    class(param), intent(in) :: self

    damp = self%damp
    
    return
  end function param_get_damp

  !---------------------------------------------------------------------

  double precision function param_get_t_pre(self) result(t_pre)
    class(param), intent(in) :: self

    t_pre = self%t_pre

    return
  end function param_get_t_pre

  !---------------------------------------------------------------------

  character(1) function param_get_rf_phase(self) result(rf_phase)
    class(param), intent(in) :: self

    rf_phase = self%rf_phase

    return
  end function param_get_rf_phase
  
  !---------------------------------------------------------------------
  
  logical function param_get_deconv_flag(self) result(deconv_flag)
    class(param), intent(in) :: self

    deconv_flag = self%deconv_flag

    return
  end function param_get_deconv_flag

  !---------------------------------------------------------------------

  logical function param_get_correct_amp(self) result(correct_amp)
    class(param), intent(in) :: self

    correct_amp = self%correct_amp

    return
  end function param_get_correct_amp

  !---------------------------------------------------------------------
  
  character(len=line_max) function param_get_recv_func_out(self) &
       & result(recv_func_out)
    class(param), intent(in) :: self
    
    recv_func_out = self%recv_func_out
    
    return
  end function param_get_recv_func_out

  !---------------------------------------------------------------------

  double precision function param_get_vp_bottom(self) result(vp_bottom)
    class(param), intent(in) :: self
    
    vp_bottom = self%vp_bottom 
    
    return 
  end function param_get_vp_bottom
  
  !---------------------------------------------------------------------

  double precision function param_get_vs_bottom(self) result(vs_bottom)
    class(param), intent(in) :: self
    
    vs_bottom = self%vs_bottom 
    
    return 
  end function param_get_vs_bottom
  
  !---------------------------------------------------------------------

  double precision function param_get_rho_bottom(self) result(rho_bottom)
    class(param), intent(in) :: self
    
    rho_bottom = self%rho_bottom 
    
    return 
  end function param_get_rho_bottom

  !---------------------------------------------------------------------
  
  double precision function param_get_noise_added(self) &
       & result(noise_added)
    class(param), intent(in) :: self

    noise_added = self%noise_added
    
    return 
  end function param_get_noise_added
  
  !---------------------------------------------------------------------

  logical function param_get_is_sphere(self) result(is_sphere)
    class(param), intent(in) :: self
    
    is_sphere = self%is_sphere
  
    return 
  end function param_get_is_sphere
  
  !---------------------------------------------------------------------

  double precision function param_get_r_earth(self) result(r_earth)
    class(param), intent(in) :: self
    
    r_earth = self%r_earth
    
    return 
  end function param_get_r_earth

  !---------------------------------------------------------------------

  logical function param_get_diagnostic_mode(self) result(diagnostic_mode)
    class(param), intent(in) :: self
    
    diagnostic_mode = self%diagnostic_mode
    
    return 
  end function param_get_diagnostic_mode

  !---------------------------------------------------------------------

  subroutine param_check_mcmc_params(self, is_ok)
    class(param), intent(in) :: self
    logical, intent(out) :: is_ok

    is_ok = .true.

    if (self%n_iter <= 0) then
       if (self%verb) write(0,*)"ERROR: n_iter must > 0"
       is_ok = .false.
    end if
    if (self%n_iter <= self%n_burn) then
       if (self%verb) write(0,*)"ERROR: n_iter must > n_burn"
       is_ok = .false.
    end if
    if (self%n_iter - self%n_burn <= self%n_corr) then
       if (self%verb) write(0,*)"ERROR: n_iter - n_burn must > n_corr"
       is_ok = .false.
    end if
    if (self%temp_high <= 1.d0) then
       if (self%verb) write(0,*)"ERROR: temp_high must > 1.0"
       is_ok = .false.
    end if
    if (self%n_chain <= 0) then
       if (self%verb) write(0,*)"ERROR: n_chain must > 0"
       is_ok = .false.
    end if
    if (self%n_chain < self%n_cool) then
       if (self%verb) write(0,*)"ERROR: n_chain must >= n_cool"
       is_ok = .false.
    end if
    if (self%k_min <= 0) then
       if (self%verb) write(0,*)"ERROR: k_min must > 0"
       is_ok = .false.
    end if
    if (self%k_max <= self%k_min) then
       if (self%verb) write(0,*)"ERROR: k_max must > k_min"
       is_ok = .false.
    end if
    if (self%z_min < 0.d0) then
       if (self%verb) write(0,*)"ERROR: z_min must >= 0.d0"
       is_ok = .false.
    end if
    if (self%z_max <= self%z_min) then
       if (self%verb) write(0,*)"ERROR: z_max must > z_min"
       is_ok = .false.
    end if
    if (self%vp_min < 0.d0) then
       if (self%verb) write(0,*)"ERROR: vp_min must >= 0.0"
       is_ok = .false.
    end if
    if (self%vp_max <= self%vp_min) then
       if (self%verb) write(0,*)"ERROR: vp_max must > vp_min"
       is_ok = .false.
    end if
    if (self%vs_min < 0.d0) then
       write(0,*)"ERROR: vs_min must >= 0.0"
       is_ok = .false.
    end if
    if (self%vs_max <= self%vs_min) then
       if (self%verb) write(0,*)"ERROR: vs_max must > vs_min"
       is_ok = .false.
    end if
    if (self%vs_max <= self%vs_bottom) then
       if (self%verb) write(0,*)"ERROR: vs_max must > vs_bottom"
       is_ok = .false.
    end if
    if (self%vp_max <= self%vp_bottom) then
       if (self%verb) write(0,*)"ERROR: vp_max must > vp_bottom"
       is_ok = .false.
    end if
    if (self%dev_z <= 0.d0) then
       if (self%verb) write(0,*)"ERROR: dev_z must > 0.0"
       is_ok = .false.
    end if
    if (self%dev_vp <= 0.d0) then
       if (self%verb) write(0,*)"ERROR: dev_vp must > 0.0"
       is_ok = .false.
    end if
    if (self%dev_vs <= 0.d0) then
       if (self%verb) write(0,*)"ERROR: dev_vs must > 0.0"
       is_ok = .false.
    end if
    if (self%n_bin_z <= 0) then
       if (self%verb) write(0,*)"ERROR: n_bin_z must > 0"
       is_ok = .false.
    end if
    if (self%n_bin_vp <= 0) then
       if (self%verb) write(0,*)"ERROR: n_bin_vp must > 0"
       is_ok = .false.
    end if
    if (self%n_bin_vs <= 0) then
       if (self%verb) write(0,*)"ERROR: n_bin_vs must > 0"
       is_ok = .false.
    end if

    if (self%is_ocean .and. self%ocean_thick > self%z_min) then
       if (self%verb) write(0,*)"ERROR: z_min must > " // &  
            & "ocean_thick for ocean mode"
       is_ok = .false.
    end if

    if (trim(self%recv_func_in) == "" &
         & .and. trim(self%disper_in) == "") then
       if (self%verb) write(0,*)"ERROR: either recv_func_in or " // &
            & "disper_in must be specified"
       is_ok = .false.
    end if

       

    return 
  end subroutine param_check_mcmc_params

  !---------------------------------------------------------------------

  !ogical function param_check_recv_func_fwd_params(self) result(is_ok)
  ! class(param), intent(in) :: self
  !
  ! write(*,*)self%given
  ! is_ok = .true.
  ! if (.not. self%is_given(self%recv_func_out, 'recv_func_out')) then
  !    is_ok = .false.
  ! else if (.not. self%is_given(self%vmod_in, 'vmod_in')) then
  !    is_ok = .false.
  ! else if (.not. self%is_given(self%rayp, 'rayp')) then
  !    is_ok = .false.
  ! else if (.not. self%is_given(self%rf_phase, 'rf_phase')) then
  !    is_ok = .false.
  ! else if (.not. self%is_given(self%a_gauss, 'a_gauss')) then
  !    is_ok = .false.
  ! else if (.not. self%is_given(self%t_pre, 't_pre')) then
  !    is_ok = .false.
  ! else if (.not. self%is_given(self%n_smp, 'n_smp')) then
  !    is_ok = .false.
  ! else if (.not. self%is_given(self%delta, 'delta')) then
  !    is_ok = .false.
  ! end if
  ! 
  ! return 
  !nd function param_check_recv_func_fwd_params
  !
  !---------------------------------------------------------------------
  !
  !ogical function param_check_disper_fwd_params(self) result(is_ok)
  ! class(param), intent(in) :: self
  !
  ! is_ok = .true.
  ! if (.not. self%is_given(self%disper_out, 'disper_out')) then
  !    is_ok = .false.
  ! else if (.not. self%is_given(self%vmod_in, 'vmod_in')) then
  !    is_ok = .false.
  ! else if (.not. self%is_given(self%freq_or_period, &
  !      & 'freq_or_period')) then
  !    is_ok = .false.
  ! else if (.not. self%is_given(self%xmin, 'xmin')) then
  !    is_ok = .false.
  ! else if (.not. self%is_given(self%xmax, 'xmax')) then
  !    is_ok = .false.
  ! else if (.not. self%is_given(self%dx, 'dx')) then
  !    is_ok = .false.
  ! else if (.not. self%is_given(self%cmin, 'cmin')) then
  !    is_ok = .false.
  ! else if (.not. self%is_given(self%cmax, 'cmax')) then
  !    is_ok = .false.
  ! else if (.not. self%is_given(self%dc, 'dc')) then
  !    is_ok = .false.
  ! else if (.not. self%is_given(self%disper_phase, 'disper_phase')) then
  !    is_ok = .false.
  ! else if (.not. self%is_given(self%n_mode, 'n_mode')) then
  !    is_ok = .false.
  ! end if
  !
  ! 
  ! return 
  !nd function param_check_disper_fwd_params
  !
  !---------------------------------------------------------------------
  !
  !ogical function param_is_given_double(self, x, key) result(is_given)
  ! class(param), intent(in) :: self
  ! double precision, intent(in) :: x
  ! character(*), intent(in) :: key
  !  
  ! is_given = .true.
  ! if (abs(x + 12345.d0) <= 1.0e-8) then
  !    is_given = .false.
  !    if (self%verb) then
  !       write(0,*)"ERROR: ", trim(key), " is not given"
  !     end if
  !  end if
  !  
  !  return 
  !end function param_is_given_double
  !
  !!---------------------------------------------------------------------
  !
  !logical function param_is_given_char(self, a, key) result(is_given)
  !  class(param), intent(in) :: self
  !  character,    intent(in) :: a
  !  character(*), intent(in) :: key
  !  
  !  is_given = .true.
  !  if (a ==  "") then
  !     is_given = .false.
  !     if (self%verb) then
  !        write(0,*)"ERROR: ", trim(key), " is not given"
  !     end if
  !  end if
  !
  !  return 
  !end function param_is_given_char
  !
  !!---------------------------------------------------------------------
  !
  !logical function param_is_given_integer(self, i, key) result(is_given)
  !  class(param), intent(in) :: self
  !  integer, intent(in) :: i
  !  character(*), intent(in) :: key
  !  
  !  is_given = .true.
  !  if (i == -12345) then
  !     is_given = .false.
  !     if (self%verb) then
  !        write(0,*)"ERROR: ", trim(key), " is not given"
  !     end if
  !  end if
  !  
  !  return 
  !end function param_is_given_integer
 
   !---------------------------------------------------------------------

   logical function param_check_given_param(self, required_params) &
        & result(is_ok)
     class(param), intent(in) :: self
     character(*), intent(in) :: required_params(:)
     logical :: found
     integer :: i

     is_ok = .true.
     do i = 1, size(required_params)
        found = any(self%given == required_params(i))
        if (.not. found) then
           if (self%verb) then
              write(0,*)
              write(0,*) "ERROR: ", &
                   & trim(required_params(i)), " must be given"
              write(0,*)
           end if
           is_ok = .false.
        end if
     end do
     
     return 
   end function param_check_given_param
end module cls_param
