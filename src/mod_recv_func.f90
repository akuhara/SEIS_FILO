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
module mod_recv_func
  use mod_vmodel
  use mod_signal_process
  implicit none 
  
  complex(kind(0d0)), private, parameter :: ei = (0.d0, 1.d0)
  double precision, private, parameter :: pi = acos(-1.d0)
  integer, private, parameter :: ir = 1, iz = 2
  
  type recv_func
     private
     type(vmodel) :: vmodel
     
     integer :: n
     double precision :: delta
     double precision :: df
     double precision :: rayp
     double precision :: a_gauss
     double precision :: damp = 0.001d0
     double precision :: t_pre = 5.d0
     character(len=1) :: phase
     
     logical :: deconv_flag = .true.
     logical :: is_ocean = .false.
     logical :: correct_amp = .false.
     complex(kind(0d0)) :: p_mat(4, 4)
     complex(kind(0d0)), allocatable :: f_data(:,:)
     double precision, allocatable :: t_data(:,:)
     double precision, allocatable :: rf_data(:)
     
   contains
     procedure :: set_vmodel => recv_func_set_vmodel
     procedure :: compute => recv_func_compute
     procedure :: do_propagation => recv_func_do_propagation
     procedure :: init_propagation => recv_func_init_propagation
     procedure :: solid_propagator => recv_func_solid_propagator
     procedure :: liquid_propagator => recv_func_liquid_propagator
     procedure :: water_level_deconv => recv_func_water_level_deconv
     procedure :: first_arrival => recv_func_first_arrival
     procedure :: shift_rf_data => recv_func_shift_rf_data
     procedure :: normalization => recv_func_normalization
     procedure :: get_t_data => recv_func_get_t_data
     procedure :: get_rf_data => recv_func_get_rf_data
     

  end type recv_func
  
  interface recv_func
     module procedure :: init_recv_func
  end interface recv_func
  
contains

  !---------------------------------------------------------------------

  type(recv_func) function init_recv_func(vm, n, delta, rayp, a_gauss, &
       & phase, deconv_flag, t_pre, correct_amp) result(self)
    type(vmodel), intent(in) :: vm
    integer, intent(in) :: n
    double precision, intent(in) :: delta
    double precision, intent(in) :: rayp
    double precision, intent(in) :: a_gauss
    character(*), intent(in) :: phase
    logical, intent(in), optional :: deconv_flag
    logical, intent(in), optional :: correct_amp
    double precision, intent(in), optional :: t_pre
    
    self%vmodel = vm
    self%n = n
    self%delta = delta
    self%df = 1.d0 / (n * delta)
    self%rayp = rayp
    self%a_gauss = a_gauss
    self%phase = phase
    if (present(deconv_flag)) then
       self%deconv_flag = deconv_flag
    end if
    if (present(correct_amp)) then
       self%correct_amp = correct_amp
    end if
    if (present(t_pre)) then
       self%t_pre = t_pre
    end if
    if (self%vmodel%get_vs(1) < 0.d0) then
       self%is_ocean = .true.
    else
       self%is_ocean = .false.
    end if
    
    

    allocate(self%rf_data(self%n))
    allocate(self%f_data(self%n, 2))
    allocate(self%t_data(self%n, 2))

    self%f_data(1:self%n,1:2) = (0.d0, 0.d0)
    self%t_data(1:self%n,1:2) = 0.d0
    self%rf_data(1:self%n) = 0.d0

    return 
  end function init_recv_func

  !---------------------------------------------------------------------

  subroutine recv_func_set_vmodel(self, vm)
    class(recv_func), intent(inout) :: self
    type(vmodel), intent(in) :: vm
    
    self%vmodel = vm
    
    return 
  end subroutine recv_func_set_vmodel

  !---------------------------------------------------------------------
  
  subroutine recv_func_compute(self)
    class(recv_func), intent(inout) :: self
    type(signal_process) :: sp
    double precision :: t_arrival, fac
    integer :: i_src, i_target
    
    sp = init_signal_process(self%n, self%delta)

    ! Calculate synthetic Green's function in the freq. domain
    call self%do_propagation()
    self%f_data(:, ir) = conjg(self%f_data(:, ir))
    self%f_data(:, iz) = -conjg(self%f_data(:, iz))
    

    ! Freq.-to-time conversion
    call sp%set_gaussian_filter(self%a_gauss)
    if (self%deconv_flag) then
       ! W/ deconvolution
       call sp%set_f_data(self%water_level_deconv()) ! deconv.
       call sp%apply_filter()
       call sp%inverse_fft()
       self%rf_data(:) = sp%get_t_data()
       call self%shift_rf_data(self%t_pre)
       if (self%correct_amp) then
          call self%normalization()
       end if
    else 
       ! W/O deconvolution 
       if (self%phase == "P") then
          i_src = iz
          i_target = ir
       else
          i_src = ir
          i_target = iz
       end if

       call sp%set_f_data(self%f_data(:, i_src))
       call sp%inverse_fft()
       fac = maxval(sp%get_t_data())
       
       if (self%phase == "P") then
          call sp%set_f_data(self%f_data(:, i_target))
       else
          write(0,*)"ERROR: Deconv_mode = .ture." // &
               & " works only for P"
          stop
       end if
       call sp%apply_filter()
       call sp%inverse_fft()
       self%rf_data(:) = sp%get_t_data() / fac
       t_arrival = self%first_arrival()
       
       call self%shift_rf_data(self%t_pre - t_arrival)
       if (self%correct_amp) then
          call self%normalization()
       end if
    end if
    
    return 
  end subroutine recv_func_compute

  !---------------------------------------------------------------------

  subroutine recv_func_do_propagation(self)
    class(recv_func), intent(inout) :: self
    integer :: iomg, ilay, ilay0, nlay
    double precision :: omega
    complex(kind(0d0)) :: denom, lq(2, 2), a, b

    
    if (self%is_ocean) then
       ilay0 = 2
    else
       ilay0 = 1
    end if
    
    nlay = self%vmodel%get_nlay()

    do iomg = 1, self%n/2 + 1
       omega = (iomg - 1) * self%df * 2.d0 * pi
       if (iomg == 1) then
          omega = 1.0e-5
       end if
       
       call self%init_propagation(omega)
       
       do ilay = nlay - 1, ilay0, -1
          self%p_mat = &
               & matmul(self%p_mat, self%solid_propagator(ilay, omega))
       end do
       
       if (.not. self%is_ocean) then
          ! On-land
          denom = self%p_mat(3,1) * self%p_mat(4,2) - &
                & self%p_mat(3,2) * self%p_mat(4,1)
          if (self%phase == "P") then
             self%f_data(iomg, ir) =   self%p_mat(4,2) / denom
             self%f_data(iomg, iz) = - self%p_mat(4,1) / denom
          else if (self%phase == "S") then
             self%f_data(iomg, ir) = - self%p_mat(3,2) / denom
             self%f_data(iomg, iz) =   self%p_mat(3,1) / denom
          else
             write(0,*)"ERROR: invalid phase type:",  self%phase
             stop
          end if
       else
          ! Ocean-bottom
          lq =  self%liquid_propagator(omega)
          a = self%p_mat(4,2)*lq(1,1) + self%p_mat(4,4)*lq(2,1)
          b = self%p_mat(3,2)*lq(1,1) + self%p_mat(3,4)*lq(2,1)
          denom = a * self%p_mat(3,1) - b * self%p_mat(4,1)
          if (self%phase == "P") then 
             self%f_data(iomg, ir) = a / denom
             self%f_data(iomg, iz) = &
                  & - lq(1,1)* self%p_mat(4,1) / denom
          else if (self%phase == "S") then
             self%f_data(iomg, ir) = - b / denom
             self%f_data(iomg, iz) = &
                  & lq(1,1) * self%p_mat(3,1) / denom
          else
             write(0,*)"ERROR: invalid phase type:",  self%phase
             stop
          end if
       end if
    end do
    
    return 
  end subroutine recv_func_do_propagation
  
  !---------------------------------------------------------------------

  subroutine recv_func_init_propagation(self, omega)
    class(recv_func), intent(inout) :: self
    double precision, intent(in) :: omega
    double precision :: alpha, beta, p, rho, eta, xi, bp
    integer :: nlay
    
    nlay  = self%vmodel%get_nlay()
    alpha = self%vmodel%get_vp(nlay)
    beta  = self%vmodel%get_vs(nlay)
    rho   = self%vmodel%get_rho(nlay)
    p     = self%rayp
    
    self%p_mat(:,:) = (0.d0, 0.d0)
    eta = sqrt(1.d0/(beta*beta) - p*p)
    xi  = sqrt(1.d0/(alpha*alpha) - p*p)
    bp = 1.d0 - 2.d0*beta*beta*p*p

    self%p_mat(1,1) = beta*beta*p/alpha
    self%p_mat(1,2) = bp/(2.d0*alpha*xi)
    self%p_mat(1,3) = -p/(2.d0*omega*rho*alpha*xi) * ei
    self%p_mat(1,4) = -1.d0/(2.d0*omega*rho*alpha) * ei
    self%p_mat(2,1) = bp / (2.d0*beta*eta)
    self%p_mat(2,2) = -beta*p
    self%p_mat(2,3) = -1.0/(2.d0*omega*rho*beta) * ei
    self%p_mat(2,4) = p/(2.d0*omega*rho*beta*eta) * ei
    self%p_mat(3,1) = self%p_mat(1,1)
    self%p_mat(3,2) = - self%p_mat(1,2)
    self%p_mat(3,3) = - self%p_mat(1,3)
    self%p_mat(3,4) = self%p_mat(1,4)
    self%p_mat(4,1) = self%p_mat(2,1)
    self%p_mat(4,2) = - self%p_mat(2,2)
    self%p_mat(4,3) = - self%p_mat(2,3)
    self%p_mat(4,4) = self%p_mat(2,4)

    return 
  end subroutine recv_func_init_propagation

  !---------------------------------------------------------------------

  function recv_func_solid_propagator(self, i, omega) result(rslt)
    class(recv_func), intent(inout) :: self
    integer, intent(in) :: i
    double precision, intent(in) :: omega
    double precision :: alpha, beta, rho, h, p
    double precision :: eta, xi, beta2, p2, bp
    double precision :: cos_xi, cos_eta, sin_xi, sin_eta
    complex(kind(0d0)) :: rslt(4, 4)

    alpha = self%vmodel%get_vp(i)
    beta  = self%vmodel%get_vs(i)
    rho   = self%vmodel%get_rho(i)
    h     = self%vmodel%get_h(i)
    p     = self%rayp
    
    beta2 = beta*beta
    p2 =p*p
    bp = 1.d0 -2.d0*beta2*p2
    eta = sqrt(1.d0/(beta2) - p2)
    xi  = sqrt(1.d0/(alpha*alpha) - p2)
    cos_xi = cos(omega*xi*h)
    cos_eta = cos(omega*eta*h)
    sin_xi = sin(omega*xi*h)
    sin_eta = sin(omega*eta*h)

    rslt(1,1) = 2.d0*beta2*p2*cos_xi + bp*cos_eta
    rslt(2,1) = p*( 2.d0*beta2*xi*sin_xi - bp/eta*sin_eta ) * ei
    rslt(3,1) = omega*rho*&
         & ( -4.d0*beta2*beta2*p2*xi*sin_xi - bp*bp/eta*sin_eta )
    rslt(4,1) = 2.d0*omega*beta2*rho*p*bp*( cos_xi - cos_eta ) * ei
    rslt(1,2) = p*( bp/xi*sin_xi - 2.d0*beta2*eta*sin_eta ) * ei
    rslt(2,2) = bp*cos_xi + 2.d0*beta2*p2*cos_eta
    rslt(3,2) = rslt(4,1)
    rslt(4,2) = -omega*rho*&
         & ( bp*bp/xi*sin_xi + 4.d0*beta2*beta2*p2*eta*sin_eta  )    
    rslt(1,3) = (p2/xi*sin_xi + eta*sin_eta)/(omega*rho)
    rslt(2,3) = p*(-cos_xi + cos_eta)/(omega*rho) * ei  
    rslt(3,3) = rslt(1,1)
    rslt(4,3) = rslt(1,2)  
    rslt(1,4) = rslt(2,3)
    rslt(2,4) = (xi*sin_xi + p2/eta*sin_eta)/(omega*rho)
    rslt(3,4) = rslt(2,1)
    rslt(4,4) = rslt(2,2)
    
    return 
  end function recv_func_solid_propagator

  !---------------------------------------------------------------------

  function recv_func_liquid_propagator(self, omega) result(rslt)
    class(recv_func), intent(inout) :: self
    double precision, intent(in) :: omega
    double precision :: alpha, h, rho, p
    double precision :: xi, cos_xi, sin_xi, g
    complex(kind(0d0)) :: rslt(2, 2)
    
    alpha = self%vmodel%get_vp(1)
    rho   = self%vmodel%get_rho(1)
    h     = self%vmodel%get_h(1)
    p     = self%rayp

    xi  = sqrt(1.d0/(alpha*alpha) - p * p)
    cos_xi = cos(omega*xi*h)
    sin_xi = sin(omega*xi*h)
    g = rho * omega / xi
    
    rslt(1,1) = cos_xi
    rslt(1,2) = sin_xi / g
    rslt(2,1) = -g * sin_xi
    rslt(2,2) = cos_xi
  
    return 
  end function recv_func_liquid_propagator

  !---------------------------------------------------------------------

  function recv_func_water_level_deconv(self) result(rslt)
    class(recv_func), intent(inout) :: self
    double precision :: zz(self%n)
    integer :: i, i_src, i_target
    double precision :: w_level
    complex(kind(0d0)) :: rslt(self%n)
    
    
    if (self%phase == "P") then
       i_src = iz
       i_target = ir
    else
       i_src = ir
       i_target = iz
    end if
    do concurrent (i=1:self%n)
       zz(i) = real(self%f_data(i, i_src)) **2 + &
            & aimag(self%f_data(i, i_src)) **2
    end do
    w_level = self%damp * maxval(zz)
    
    do concurrent (i=1:self%n)
       rslt(i) = self%f_data(i, i_target) * &
            & conjg(self%f_data(i, i_src)) / &
            & max(zz(i), w_level)
    end do
    
    if (self%phase == "S") then
       rslt(:) = -conjg(rslt)
    end if
    
    return 
  end function recv_func_water_level_deconv

  !---------------------------------------------------------------------

  double precision function recv_func_first_arrival(self) result(t)
    class(recv_func), intent(inout) :: self
    integer :: ilay, nlay, i0
    double precision :: v
    
    nlay = self%vmodel%get_nlay()

    if (self%is_ocean) then
       i0 = 2
    else
       i0 = 1
    end if
    t = 0.d0
    do ilay = i0, nlay - 1
       if (self%phase == "P") then
         v = self%vmodel%get_vp(ilay)
      else
         v = self%vmodel%get_vs(ilay)
      end if
       t = t + self%vmodel%get_h(ilay) * &
            & sqrt(1.d0 / v * v - self%rayp * self%rayp)
    end do
    
    return 
  end function recv_func_first_arrival
  
  !---------------------------------------------------------------------

  subroutine recv_func_shift_rf_data(self, t_shift)
    class(recv_func), intent(inout) :: self
    double precision, intent(in) :: t_shift ! 
    integer :: i, j, n_shift
    double precision :: tmp(self%n)
    write(*,*)self%n
    write(*,*)size(self%rf_data), size(self%f_data)


    tmp(:) = self%rf_data(:)
    
    n_shift = nint(t_shift / self%delta)
    !if (self%phase == "P") then
       do i = 1, self%n
          j = mod(self%n - n_shift + i, self%n)
          if (j == 0) then
             j = self%n
          end if
          self%rf_data(i) = tmp(j)
       end do
    !else if (self%phase == "S") then
    !   do i = 1, self%n
    !      j = mod(self%n + n_shift - i, self%n)
    !      if (j == 0) then
    !         j = self%n
    !      end if
    !      self%rf_data(i) = -tmp(j)
    !   end do
    !else
    !   write(0,*)"ERROR: invalid phase name" , self%phase 
    !   stop
    !end if

    return 
  end subroutine recv_func_shift_rf_data

  !---------------------------------------------------------------------

  subroutine recv_func_normalization(self)
    class(recv_func), intent(inout) :: self
    double precision :: fac
    if (self%a_gauss > 0.d0) then
       fac = self%a_gauss * self%delta / sqrt(pi)
    else
       fac = 1.d0
    end if
    self%rf_data(:) = self%rf_data(:) / fac

    return
  end subroutine recv_func_normalization

  !---------------------------------------------------------------------

  function recv_func_get_t_data(self, i) result(t_data)
    class(recv_func), intent(in) :: self
    integer, intent(in) :: i
    double precision :: t_data(self%n)
    
    t_data(:) = self%t_data(:,i)
    
    return 
  end function recv_func_get_t_data

  !---------------------------------------------------------------------

  function recv_func_get_rf_data(self) result(rf_data)
    class(recv_func), intent(in) :: self
    double precision :: rf_data(self%n)
    
    rf_data(1:self%n) = self%rf_data(1:self%n)
    
    return 
  end function recv_func_get_rf_data

  !---------------------------------------------------------------------
end module mod_recv_func
