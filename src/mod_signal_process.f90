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
module mod_signal_process
  implicit none
  include 'fftw3.f'

  double precision, private, parameter :: pi = acos(-1.d0)

  !---------------------------------------------------------------------
  ! signal_process
  !---------------------------------------------------------------------
  type signal_process
     private

     integer :: n
     double precision :: dt
     double precision :: df
     double precision, allocatable   :: t_data(:)
     complex(kind(0d0)), allocatable :: f_data(:)
     integer(8) :: ifft_fwd
     integer(8) :: ifft_inv
     double precision, allocatable :: filter(:)
     
   contains
     procedure :: forward_fft => signal_process_forward_fft
     procedure :: inverse_fft  => signal_process_inverse_fft
     procedure, private :: set_rx   => signal_process_set_rx
     procedure, private :: set_cx   => signal_process_set_cx
     procedure :: set_gaussian_filter => signal_process_set_gaussian_filter
     procedure :: apply_filter => signal_process_apply_filter
  end type signal_process

  interface signal_process
     module procedure init_signal_process
  end interface signal_process

contains
  !---------------------------------------------------------------------
  
  type(signal_process) function init_signal_process(n, dt) result(self)
   integer, intent(in) :: n
    double precision, intent(in) :: dt
    
    self%n = n
    self%dt = dt
    self%df = 1.d0 / (n * dt)
    allocate(self%f_data(n), self%t_data(n), self%filter(n))
    self%f_data(1:n) = (0.d0, 0.d0)
    self%t_data(1:n) = 0.d0
    self%filter(1:n) = 0.d0
    call dfftw_plan_dft_c2r_1d(self%ifft_inv, n, &
         & self%f_data, self%t_data, FFTW_ESTIMATE)
    call dfftw_plan_dft_r2c_1d(self%ifft_fwd, n, &
         & self%t_data, self%f_data, FFTW_ESTIMATE)
    
    return 
  end function init_signal_process

  !---------------------------------------------------------------------
  
  function signal_process_forward_fft(self, rx_in) result(cx_out)
    class(signal_process), intent(inout) :: self
    double precision, intent(in) :: rx_in(1:self%n)
    complex(kind(0d0)) :: cx_out(self%n)

    call self%set_rx(rx_in)
    call dfftw_execute(self%ifft_fwd)
    
    cx_out = self%f_data

    return 
  end function signal_process_forward_fft
  
  !---------------------------------------------------------------------

  function signal_process_inverse_fft(self, cx_in) result(rx_out)
    class(signal_process), intent(inout) :: self
    complex(kind(0d0)), intent(in) :: cx_in(1:self%n)
    double precision :: rx_out(self%n)
    
    call self%set_cx(cx_in)
    call dfftw_execute(self%ifft_inv)
    
    rx_out = self%t_data / self%n
    
    return 
  end function signal_process_inverse_fft

  !---------------------------------------------------------------------

  subroutine signal_process_set_rx(self, rx_in)
    class(signal_process), intent(inout) :: self
    double precision, intent(in) :: rx_in(self%n)
    
    self%t_data = rx_in

    return 
  end subroutine signal_process_set_rx
    
  !---------------------------------------------------------------------

  subroutine signal_process_set_cx(self, cx_in)
    class(signal_process), intent(inout) :: self
    complex(kind(0d0)), intent(in) :: cx_in(self%n)
    
    self%f_data = cx_in
    
    return 
  end subroutine signal_process_set_cx
  
  !---------------------------------------------------------------------
  
  function signal_process_apply_filter(self, fdata) result(out)
    class(signal_process), intent(inout) :: self
    complex(kind(0d0)), intent(inout) :: fdata(1:self%n)
    complex(kind(0d0)) :: out(1:self%n)
    
    out = fdata * self%filter
    

    return 
  end function signal_process_apply_filter

  !---------------------------------------------------------------------

  subroutine signal_process_set_gaussian_filter(self, a_gauss)
    class(signal_process), intent(inout) :: self
    double precision, intent(in) :: a_gauss
    integer :: i
    double precision :: omega, fac
    
    do i = 1, self%n / 2 + 1
       omega = (i - 1) * 2.d0 * pi * self%df
       self%filter(i) = &
            & exp(-(omega / (2.d0 * a_gauss))**2)
    end do
    

    return 
  end subroutine signal_process_set_gaussian_filter

  
end module mod_signal_process


