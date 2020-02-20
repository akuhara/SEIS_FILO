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
module cls_signal_process
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
     procedure :: set_t_data => signal_process_set_t_data
     procedure :: set_f_data => signal_process_set_f_data
     procedure :: get_t_data => signal_process_get_t_data
     procedure :: get_f_data => signal_process_get_f_data
     procedure :: forward_fft => signal_process_forward_fft
     procedure :: inverse_fft  => signal_process_inverse_fft
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
  
  subroutine signal_process_set_t_data(self, t_data)
    class(signal_process), intent(inout) :: self
    double precision, intent(in) :: t_data(self%n)
    
    self%t_data = t_data
    
    return 
  end subroutine signal_process_set_t_data

  !---------------------------------------------------------------------

  subroutine signal_process_set_f_data(self, f_data)
    class(signal_process), intent(inout) :: self
    complex(kind(0d0)), intent(in) :: f_data(self%n)

    self%f_data = f_data

    return 
  end subroutine signal_process_set_f_data

  !---------------------------------------------------------------------
  
  function signal_process_get_t_data(self) result(t_data)
    class(signal_process), intent(inout) :: self
    double precision :: t_data(self%n)
    
    t_data(:) = self%t_data(:)

    return 
  end function signal_process_get_t_data
    
  !---------------------------------------------------------------------

  function signal_process_get_f_data(self) result(f_data)
    class(signal_process), intent(inout) :: self
    complex(kind(0d0)) :: f_data(self%n)
    
    f_data(:) = self%f_data(:)
    
    return 
  end function signal_process_get_f_data
    
  !---------------------------------------------------------------------

  subroutine signal_process_forward_fft(self)
    class(signal_process), intent(inout) :: self

    call dfftw_execute(self%ifft_fwd)

    return 
  end subroutine signal_process_forward_fft

  !---------------------------------------------------------------------

  subroutine signal_process_inverse_fft(self)
    class(signal_process), intent(inout) :: self

    call dfftw_execute(self%ifft_inv)
    self%t_data(:) = self%t_data(:) / self%n

    return 
  end subroutine signal_process_inverse_fft

  !---------------------------------------------------------------------

  subroutine signal_process_apply_filter(self)
    class(signal_process), intent(inout) :: self

    self%f_data(:) = self%f_data(:) * self%filter(:)
    
    return 
  end subroutine signal_process_apply_filter

  !---------------------------------------------------------------------

  subroutine signal_process_set_gaussian_filter(self, a_gauss)
    class(signal_process), intent(inout) :: self
    double precision, intent(in) :: a_gauss
    integer :: i
    double precision :: omega
    
    if (a_gauss >= 0.d0) then
       !do i = 1, self%n / 2 + 1
       do concurrent (i=1:self%n)
          omega = (i - 1) * 2.d0 * pi * self%df
          self%filter(i) = &
               & exp(-(omega / (2.d0 * a_gauss))**2)
       end do
    else
       self%filter(:) = 1.d0
    end if
    

    return 
  end subroutine signal_process_set_gaussian_filter

  
end module cls_signal_process


