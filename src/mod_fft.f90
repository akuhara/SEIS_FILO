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
!> Wrapper of FFTW by objective-oriented programming
!> @author
!! Takeshi Akuhara
module mod_fft
  implicit none
  include 'fftw3.f'

  !---------------------------------------------------------------------
  ! fft
  !---------------------------------------------------------------------
  type fft
     private

     integer :: n
     complex(kind(0d0)), allocatable :: cx(:)
     double precision, allocatable   :: rx(:)
     integer(8) :: ifft_fwd
     integer(8) :: ifft_inv
     
   contains
     procedure :: forwawrd => fft_forward
     procedure :: inverse  => fft_inverse
     !procedure :: finish   => fft_finish
     procedure, private :: set_rx   => fft_set_rx
     procedure, private :: set_cx   => fft_set_cx
     
  end type fft

  interface fft
     module procedure init_fft
  end interface fft

contains
  
  !---------------------------------------------------------------------

  type(fft) function init_fft(n) result(self)
    integer, intent(in) :: n
    
    self%n = n

    allocate(self%cx(self%n), self%rx(self%n))
    self%cx(1:n) = (0.d0, 0.d0)
    self%rx(1:n) = 0.d0
    call dfftw_plan_dft_c2r_1d(self%ifft_inv, self%n, &
         & self%cx, self%rx, FFTW_ESTIMATE)
    call dfftw_plan_dft_r2c_1d(self%ifft_fwd, self%n, &
         & self%rx, self%cx, FFTW_ESTIMATE)
    
    
    return 
  end function init_fft

  !---------------------------------------------------------------------

  function fft_forward(self, rx_in) result(cx_out)
    class(fft), intent(inout) :: self
    double precision, intent(in) :: rx_in(1:self%n)
    complex(kind(0d0)) :: cx_out(self%n)

    call self%set_rx(rx_in)
    call dfftw_execute(self%ifft_fwd)
    
    cx_out = self%cx

    return 
  end function fft_forward
  
  !---------------------------------------------------------------------

  function fft_inverse(self, cx_in) result(rx_out)
    class(fft), intent(inout) :: self
    complex(kind(0d0)), intent(in) :: cx_in(1:self%n)
    double precision :: rx_out(self%n)

    call self%set_cx(cx_in)
    call dfftw_execute(self%ifft_inv)
    
    rx_out = self%rx
    
    return 
  end function fft_inverse

  !---------------------------------------------------------------------

  subroutine fft_set_rx(self, rx_in)
    class(fft), intent(inout) :: self
    double precision, intent(in) :: rx_in(self%n)
    
    self%rx = rx_in

    return 
  end subroutine fft_set_rx
    
  !---------------------------------------------------------------------

  subroutine fft_set_cx(self, cx_in)
    class(fft), intent(inout) :: self
    complex(kind(0d0)), intent(in) :: cx_in(self%n)
    
    self%cx = cx_in
    
    return 
  end subroutine fft_set_cx
  
  !---------------------------------------------------------------------
  
  
end module mod_fft


