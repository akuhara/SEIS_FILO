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
module cls_interpreter
  use cls_trans_d_model
  use cls_hyper_model
  use cls_vmodel
  use mod_sort
  use mod_const
  use cls_line_text, only: line_max
  implicit none 
  
  type interpreter
     private
     integer :: nlay_max

     logical :: is_ocean = .false.
     double precision :: ocean_thick = 0.d0
     double precision :: vp_ocean    = 1.5d0
     double precision :: rho_ocean   = 1.d0
     double precision :: vs_bottom   = 4.6d0
     double precision :: vp_bottom   = 8.1d0
     double precision :: rho_bottom  = 3.3d0

     integer :: n_bin_z
     integer :: n_bin_vs
     integer :: n_bin_vp
     integer :: n_bin_sig

     integer :: n_disp
     integer :: n_rf

     integer, allocatable :: n_vsz(:, :)
     integer, allocatable :: n_vpz(:, :)
     integer, allocatable :: n_layers(:)
     integer, allocatable :: n_rf_sig(:, :)
     integer, allocatable :: n_disp_sig(:,:)

     double precision, allocatable :: vsz_mean(:)
     double precision, allocatable :: vpz_mean(:)
     
     double precision :: z_min
     double precision :: z_max
     double precision :: dz
     double precision :: vs_min
     double precision :: vs_max
     double precision :: dvs
     double precision :: vp_min
     double precision :: vp_max
     double precision :: dvp
     
     logical :: solve_anomaly = .false.
     logical :: solve_vp = .true.
     logical :: is_sphere = .false.

     logical :: verb = .false.

     character(line_max) :: ref_vmod_in
     double precision :: dvs_sig
     double precision :: dvp_sig
     double precision :: dz_ref
     double precision, allocatable :: vs_ref(:)
     double precision, allocatable :: vp_ref(:)
     

     double precision, allocatable :: wrk_vp(:), wrk_vs(:), wrk_z(:) 
     
   contains
     procedure :: construct_vmodel => interpreter_construct_vmodel
     procedure :: save_model => interpreter_save_model
     procedure :: save_rf_sigma => interpreter_save_rf_sigma
     procedure :: save_disp_sigma => interpreter_save_disp_sigma
     procedure :: get_n_bin_z => interpreter_get_n_bin_z
     procedure :: get_n_bin_vs => interpreter_get_n_bin_vs
     procedure :: get_n_bin_vp => interpreter_get_n_bin_vp
     procedure :: get_n_vpz => interpreter_get_n_vpz
     procedure :: get_n_vsz => interpreter_get_n_vsz
     procedure :: get_vpz_mean => interpreter_get_vpz_mean
     procedure :: get_vsz_mean => interpreter_get_vsz_mean
     procedure :: get_n_layers => interpreter_get_n_layers
     procedure :: get_n_rf_sig => interpreter_get_n_rf_sig
     procedure :: get_n_disp_sig => interpreter_get_n_disp_sig
     
     procedure :: get_dvs => interpreter_get_dvs
     procedure :: get_dvp => interpreter_get_dvp
     procedure :: get_dz => interpreter_get_dz
     
     procedure :: read_ref_vmod_in => interpreter_read_ref_vmod_in
  end type interpreter
  
  interface interpreter
     module procedure init_interpreter
  end interface interpreter
  
contains

  !---------------------------------------------------------------------
  type(interpreter) function init_interpreter(nlay_max, z_min, z_max, &
       & n_bin_z, vs_min, vs_max, n_bin_vs, n_disp, n_rf, n_bin_sig, &
       & vp_min, vp_max, n_bin_vp, &
       & is_ocean, ocean_thick, vp_ocean, rho_ocean, solve_vp, &
       & solve_anomaly, is_sphere, vp_bottom, vs_bottom, rho_bottom, &
       & dvs_sig, dvp_sig, ref_vmod_in, verb) result(self)
    integer, intent(in) :: nlay_max
    integer, intent(in) :: n_bin_z, n_bin_vs
    integer, intent(in), optional :: n_bin_vp
    integer, intent(in) :: n_disp, n_rf, n_bin_sig
    double precision, intent(in) :: z_min, z_max, vs_min, vs_max
    double precision, intent(in), optional :: vp_min, vp_max
    logical, intent(in), optional :: is_ocean

    double precision, intent(in), optional :: ocean_thick, &
         & vp_ocean, rho_ocean
    double precision, intent(in), optional :: vp_bottom, &
         & vs_bottom, rho_bottom
    double precision, intent(in), optional :: dvs_sig, dvp_sig
    logical, intent(in), optional :: solve_vp
    logical, intent(in), optional :: solve_anomaly
    logical, intent(in), optional :: is_sphere

    character(line_max), intent(in), optional :: ref_vmod_in
    logical, intent(in), optional :: verb
    integer :: ierr

    if (present(verb)) then
       self%verb = verb
    end if
    
    if (self%verb) then
       write(*,'(a)')"<< Initialize model interpreter >>"
    end if
    
    self%nlay_max = nlay_max
    allocate(self%wrk_vp(nlay_max + 1))
    allocate(self%wrk_vs(nlay_max + 1))
    allocate(self%wrk_z(nlay_max + 1))

    
    self%z_min   = z_min
    self%z_max   = z_max
    self%n_bin_z  = n_bin_z
    self%vs_min  = vs_min
    self%vs_max  = vs_max
    self%n_bin_vs = n_bin_vs
    self%n_disp = n_disp
    self%n_rf   = n_rf
    self%n_bin_sig = n_bin_sig
    
    self%dz  = (z_max - 0.d0) / dble(n_bin_z)
    self%dvs = (vs_max - vs_min) / dble(n_bin_vs)
    
    allocate(self%n_vsz(n_bin_vs, n_bin_z))
    self%n_vsz = 0
    allocate(self%n_layers(self%nlay_max))
    self%n_layers = 0
    allocate(self%vsz_mean(n_bin_z))
    self%vsz_mean = 0.d0
    allocate(self%n_rf_sig(self%n_bin_sig, self%n_rf))
    self%n_rf_sig(:,:) = 0
    allocate(self%n_disp_sig(self%n_bin_sig, self%n_disp * 2))
    self%n_disp_sig(:,:) = 0

    if (present(is_ocean)) then
       self%is_ocean = is_ocean
    end if

    if (present(ocean_thick)) then
       self%ocean_thick = ocean_thick
    end if

    if (present(vp_ocean)) then
       self%vp_ocean = vp_ocean
    end if
    
    if (present(rho_ocean)) then
       self%rho_ocean = rho_ocean
    end if

    if (present(vp_bottom)) then
       self%vp_bottom = vp_bottom
    end if
    
    if (present(vs_bottom)) then
       self%vs_bottom = vs_bottom
    end if

    if (present(rho_bottom)) then
       self%rho_bottom = rho_bottom
    end if
    
    if (present(is_sphere)) then
       self%is_sphere = is_sphere
    end if

    if (present(solve_vp)) then
       if (.not. present(vp_min)) then
          write(0,*)"ERROR: vp_min is not given"
          call mpi_finalize(ierr)
          stop
       end if
       if (.not. present(vp_max)) then
          write(0,*)"ERROR: vp_max is not given"
          call mpi_finalize(ierr)
          stop
       end if
       if (.not. present(n_bin_vp)) then
          write(0,*)"ERROR: n_bin_vp is not given"
          call mpi_finalize(ierr)
          stop
       end if
       self%solve_vp = solve_vp
       if (self%solve_vp) then
          self%vp_min = vp_min
          self%vp_max = vp_max
          self%n_bin_vp = n_bin_vp
          self%dvp = (vp_max - vp_min) / dble(n_bin_vp)
          allocate(self%n_vpz(n_bin_vp, n_bin_z))
          self%n_vpz = 0
          allocate(self%vpz_mean(n_bin_z))
          self%vpz_mean = 0.d0
       end if
    end if
    
    if (present(solve_anomaly)) then
       self%solve_anomaly = solve_anomaly
    end if

    if (self%solve_anomaly) then
       if (.not. present(dvs_sig)) then
          write(0,*)"ERROR: dvs_sig is not given"
          call mpi_finalize(ierr)
          stop
       end if
       if (.not. present(ref_vmod_in)) then
          write(0,*)"ERROR: ref_vmod_in is not given"
          call mpi_finalize(ierr)
          stop
       end if
       if (self%solve_vp) then
          if (.not. present(dvp_sig)) then
             write(0,*)"ERROR: dvp_sig is not given"
             call mpi_finalize(ierr)
             stop
          end if
       end if
       self%solve_anomaly = solve_anomaly
       self%dvs_sig = dvs_sig
       self%ref_vmod_in = ref_vmod_in
       if (self%solve_vp) then
          self%dvp_sig = dvp_sig
       end if

       call self%read_ref_vmod_in()
    end if


    return 
  end function init_interpreter

  !---------------------------------------------------------------------
  
  subroutine interpreter_construct_vmodel(self, tm, vm, is_ok)
    class(interpreter), intent(inout) :: self
    type(trans_d_model), intent(in) :: tm
    type(vmodel), intent(out) :: vm
    logical, intent(out) :: is_ok
    double precision :: thick, vel
    integer :: i, i1, k, iz
    
    k = tm%get_k()
    vm = init_vmodel()
    
    is_ok = .true.
    
    ! Set ocean layer
    if (self%is_ocean) then
       call vm%set_nlay(k + 2) ! k middle layers + 1 ocean layer + 
                               ! 1 half space
       call vm%set_vp(1, self%vp_ocean)
       call vm%set_vs(1, -999.d0)
       call vm%set_rho(1, self%rho_ocean)
       call vm%set_h(1, self%ocean_thick)
       i1 = 1
    else
       call vm%set_nlay(k + 1)
       i1 = 0
    end if
    
    ! Translate trand_d_model to vmodel
    self%wrk_z(:) = tm%get_x(id_z)
    self%wrk_vs(:) = tm%get_x(id_vs)
    if (self%solve_vp) then
       self%wrk_vp(:) = tm%get_x(id_vp)
    end if
    call quick_sort(self%wrk_z, 1, k, self%wrk_vs, self%wrk_vp)
    ! Middle layers
    do i = 1, k

       ! Thickness
       if (i == 1) then
          thick = self%wrk_z(i) - self%z_min
          !call vm%set_h(i+i1, self%wrk_z(i) - self%z_min)
       else if (i == k) then
          ! deepest layer
          thick = self%z_max - self%wrk_z(i-1)
          !call vm%set_h(i+i1, self%z_max - self%wrk_z(i-1))
       else 
          ! others
          thick = self%wrk_z(i) - self%wrk_z(i-1)
          !call vm%set_h(i+i1, self%wrk_z(i) - self%wrk_z(i-1))
       end if
       if (thick <= 0.d0) then
          is_ok = .false.
          return
       end if
       call vm%set_h(i+i1, thick)
       
       ! Velocity
       if (.not. self%solve_anomaly) then
          ! Vs
          call vm%set_vs(i+i1, self%wrk_vs(i))
          ! Vp
          if (self%solve_vp) then
             call vm%set_vp(i+i1, self%wrk_vp(i))
          else
             call vm%vs2vp_brocher(i+i1)
          end if
       else
          if (i == 1) then
             iz = nint(0.5d0 * (self%wrk_z(i) + self%z_min) / &
                  & self%dz_ref) + 1
          else if (i == k) then
             iz = nint(0.5d0 * (self%wrk_z(i-1) + self%z_max) / &
                  & self%dz_ref) + 1
          else
             iz = nint(0.5d0 * (self%wrk_z(i-1) + self%wrk_z(i)) / &
                  & self%dz_ref) + 1
          end if
          ! Vs
          vel = self%vs_ref(iz) + self%wrk_vs(i)
          !vel = self%vs_ref(iz) * (1.d0 + self%wrk_vs(i) / 100.d0)
          if (vel < self%vs_min .or. vel > self%vs_max) then
             is_ok = .false.
             return
          end if
          call vm%set_vs(i+i1, vel)
          ! Vp
          if (self%solve_vp) then
             vel = self%vp_ref(iz) + self%wrk_vp(i)
             !vel = self%vp_ref(iz) * (1.d0 + self%wrk_vp(i) / 100.d0)
             if (vel < self%vp_min .or. vel > self%vp_max) then
                is_ok = .false.
                return
             end if
             call vm%set_vp(i+i1, vel)
          else
             !call vm%vs2vp_brocher(i+i1)
             vel = self%vp_ref(iz)
             call vm%set_vp(i+i1, vel)
          end if
       end if
       ! Density
       call vm%vp2rho_brocher(i+i1)
    end do
    ! Bottom layer
    ! Vs
    if (self%vs_bottom > 0.d0) then
       call vm%set_vs(k+1+i1, self%vs_bottom) ! <- Fixed
    else
       call vm%set_vs(k+1+i1, vm%get_vs(k+i1)) ! <- Same as layer above
    end if
    ! Vp
    if (self%vp_bottom > 0.d0) then
       call vm%set_vp(k+1+i1, self%vp_bottom)
    else
       if (self%solve_vp) then
          call vm%set_vp(k+1+i1, vm%get_vp(k+i1))
       else
          call vm%vs2vp_brocher(k+1+i1)
       end if
    end if
    ! Thickness
    call vm%set_h(k+1+i1,  -99.d0) ! <- half space
    ! Density
    if (self%rho_bottom > 0.d0) then
       call vm%set_rho(k+1+i1, self%rho_bottom)
    else
       call vm%vp2rho_brocher(k+1+i1)
    end if
    
    if (any(self%wrk_vs(i1+1:i1+k+1) > self%vs_max)) then
       is_ok = .false.
    end if

    !call vm%display()
    
    return 
  end subroutine interpreter_construct_vmodel

  !---------------------------------------------------------------------

  function interpreter_get_n_vsz(self) result(n_vsz)
    class(interpreter), intent(in) :: self
    integer :: n_vsz(self%n_bin_vs, self%n_bin_z)

    n_vsz = self%n_vsz

    return 
  end function interpreter_get_n_vsz
  
  !---------------------------------------------------------------------

  function interpreter_get_n_vpz(self) result(n_vpz)
    class(interpreter), intent(in) :: self
    integer :: n_vpz(self%n_bin_vp, self%n_bin_z)

    n_vpz = self%n_vpz

    return 
  end function interpreter_get_n_vpz
  
  !---------------------------------------------------------------------

  function interpreter_get_vsz_mean(self) result(vsz_mean)
    class(interpreter), intent(in) :: self
    double precision :: vsz_mean(self%n_bin_z)

    vsz_mean = self%vsz_mean

    return 
  end function interpreter_get_vsz_mean

  !---------------------------------------------------------------------

  function interpreter_get_vpz_mean(self) result(vpz_mean)
    class(interpreter), intent(in) :: self
    double precision :: vpz_mean(self%n_bin_z)

    vpz_mean = self%vpz_mean

    return 
  end function interpreter_get_vpz_mean

  !---------------------------------------------------------------------

  function interpreter_get_n_layers(self) result(n_layers)
    class(interpreter), intent(in) :: self
    integer :: n_layers(self%nlay_max)
    
    n_layers = self%n_layers

    return 
  end function interpreter_get_n_layers

  !---------------------------------------------------------------------

  function interpreter_get_n_rf_sig(self, i) result(n_rf_sig)
    class(interpreter), intent(in) :: self
    integer, intent(in) :: i
    integer :: n_rf_sig(self%n_bin_sig)
    
    n_rf_sig(:) = self%n_rf_sig(:, i)

    return 
  end function interpreter_get_n_rf_sig

  !---------------------------------------------------------------------

  function interpreter_get_n_disp_sig(self, i) result(n_disp_sig)
    class(interpreter), intent(in) :: self
    integer, intent(in) :: i
    integer :: n_disp_sig(self%n_bin_sig)
    
    n_disp_sig(:) = self%n_disp_sig(:, i)

    return 
  end function interpreter_get_n_disp_sig

  !---------------------------------------------------------------------

  subroutine interpreter_save_model(self, tm)
    class(interpreter), intent(inout) :: self
    type(trans_d_model), intent(in) :: tm
    type(vmodel) :: vm
    integer :: nlay
    integer :: ilay, iz, iz1, iz2, iv
    double precision :: tmpz, z, vp, vs
    logical :: is_ok
    
    call self%construct_vmodel(tm, vm, is_ok)
 
    nlay = vm%get_nlay()
    
    self%n_layers(tm%get_k()) = self%n_layers(tm%get_k()) + 1
    tmpz = 0.d0
    model: do ilay = 1, nlay
       iz1 = int(tmpz / self%dz) + 1
       
       if (ilay < nlay) then
          iz2 = int((tmpz + vm%get_h(ilay)) / self%dz) + 1
       else
          iz2 = self%n_bin_z + 1
       end if
       layer: do iz = iz1, iz2 - 1
          z = (iz - 1) * self%dz
          vp = vm%get_vp(ilay)
          vs = vm%get_vs(ilay)
          
          iv = int((vs - self%vs_min) / self%dvs) + 1
          
          if (iv < 1) then
             iv = 1 ! This should occur for ocean
          end if
          if (iz > self%n_bin_z) exit model
          
          self%n_vsz(iv, iz) = self%n_vsz(iv, iz) + 1
          self%vsz_mean(iz) = self%vsz_mean(iz) + vs
          if (self%solve_vp) then
             iv = int((vp - self%vp_min) / self%dvp) + 1
             if (iv < 1) then
                iv = 1 ! This should occur for ocean
             end if
             self%n_vpz(iv, iz) = self%n_vpz(iv, iz) + 1
             self%vpz_mean(iz) = self%vpz_mean(iz) + vp
          end if
          
       end do layer
       tmpz = tmpz + vm%get_h(ilay) 
    end do model
    
  end subroutine interpreter_save_model

  !---------------------------------------------------------------------

  subroutine interpreter_save_rf_sigma(self, hyp)
    class(interpreter), intent(inout) :: self
    type(hyper_model), intent(in) :: hyp
    integer :: iobs, isig
    double precision :: del_sig
    

    do iobs = 1, hyp%get_nx()
       del_sig = (hyp%get_prior_param(iobs, 2) - &
            & hyp%get_prior_param(iobs, 1)) / self%n_bin_sig
       
       isig = int((hyp%get_x(iobs) - hyp%get_prior_param(iobs, 1)) / &
            & del_sig) + 1
       self%n_rf_sig(isig, iobs) = self%n_rf_sig(isig, iobs) + 1
    end do
    
    return 
  end subroutine interpreter_save_rf_sigma

  !---------------------------------------------------------------------

  subroutine interpreter_save_disp_sigma(self, hyp)
    class(interpreter), intent(inout) :: self
    type(hyper_model), intent(in) :: hyp
    integer :: iobs, isig
    double precision :: del_sig
    

    do iobs = 1, hyp%get_nx()
       del_sig = (hyp%get_prior_param(iobs, 2) - &
            & hyp%get_prior_param(iobs, 1)) / self%n_bin_sig
       
       isig = int((hyp%get_x(iobs) - hyp%get_prior_param(iobs, 1)) / &
            & del_sig) + 1
       self%n_disp_sig(isig, iobs) = self%n_disp_sig(isig, iobs) + 1
    end do
    
    return 
  end subroutine interpreter_save_disp_sigma

  !---------------------------------------------------------------------
  
  subroutine interpreter_read_ref_vmod_in(self)
    class(interpreter), intent(inout) :: self
    integer :: io, ierr, nz, i
    double precision :: z_old, z_tmp, dz
    double precision, parameter :: eps = 1.0e-5
    
    open(newunit = io, file = self%ref_vmod_in, status = "old", &
         & iostat = ierr)
    if (ierr /= 0) then
       if (self%verb) then
          write(0,*)"ERROR: cannot open reference velocity file:", &
               & trim(self%ref_vmod_in)
          write(0,*)"       Set 'ref_vmod_in = <Filename>' or " // &
               & "'solve_anomaly = .false.' in the parameter file"
       end if
       call mpi_finalize(ierr)
       stop
    end if

    ! get nz and dz
    nz = 0
    dz = -99.d0
    z_old = -99.d0
    do 
       read(io, *, iostat= ierr) z_tmp
       if (ierr /= 0) exit
       nz = nz + 1
       if (nz > 2 .and. abs(z_tmp - z_old - dz) > eps .and. &
            & self%verb) then
          write(0,*)"ERROR: while reading reference velocity file:", &
               & trim(self%ref_vmod_in)
          write(0,*)"       Depth increment must be constant", &
               & z_tmp - z_old, dz

          call mpi_finalize(ierr)
          stop
       end if
       dz = z_tmp - z_old
       z_old = z_tmp
    end do
    self%dz_ref = dz

    allocate(self%vs_ref(nz))
    allocate(self%vp_ref(nz))
    
    rewind(io)
    do i = 1, nz
       read(io, *)z_tmp, self%vp_ref(i), self%vs_ref(i)
    end do
    close(io)

    if (z_tmp <= self%z_max) then
       if (self%verb) then
          write(0,*)"ERROR: maximum depth in reference velocity " &
               & // "model file must be deeper than z_max in " &
               & // "main parameter file." 
          write(0,*)"       Maximum depth in ", &
               & trim(self%ref_vmod_in), &
               & ": ", z_tmp 
          write(0,*)"       z_max: ", self%z_max
       end if
       call mpi_finalize(ierr)
       stop
    end if

    if (self%verb) then
       write(*,*)"----- Reference velocity -----"
       do i = 1, nz
          write(*,'(3F9.3)')(i-1) * self%dz_ref, self%vp_ref(i), &
               & self%vs_ref(i)
       end do
    end if
       
       
    
    return 
  end subroutine interpreter_read_ref_vmod_in

  !---------------------------------------------------------------------
  
  integer function interpreter_get_n_bin_z(self) result(n_bin_z)
    class(interpreter), intent(in) ::self

    n_bin_z = self%n_bin_z
    
    return 
  end function interpreter_get_n_bin_z
  
  !---------------------------------------------------------------------

  integer function interpreter_get_n_bin_vs(self) result(n_bin_vs)
    class(interpreter), intent(in) ::self

    n_bin_vs = self%n_bin_vs
    
    return 
  end function interpreter_get_n_bin_vs
  
  !---------------------------------------------------------------------

  integer function interpreter_get_n_bin_vp(self) result(n_bin_vp)
    class(interpreter), intent(in) ::self

    n_bin_vp = self%n_bin_vp
    
    return 
  end function interpreter_get_n_bin_vp
  
  !---------------------------------------------------------------------  

  double precision function interpreter_get_dvs(self) result(dvs)
    class(interpreter), intent(in) ::self

    dvs = self%dvs
    
    return 
  end function interpreter_get_dvs
  
  !---------------------------------------------------------------------  
  
  double precision function interpreter_get_dvp(self) result(dvp)
    class(interpreter), intent(in) ::self

    dvp = self%dvp
    
    return 
  end function interpreter_get_dvp
  
  !---------------------------------------------------------------------  

  double precision function interpreter_get_dz(self) result(dz)
    class(interpreter), intent(in) ::self
    
    dz = self%dz
    
    return 
  end function interpreter_get_dz
  
  !---------------------------------------------------------------------  
  
end module cls_interpreter

