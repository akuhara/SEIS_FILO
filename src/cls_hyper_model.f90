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
module cls_hyper_model
  use mod_random
  implicit none 
  
  type hyper_model
     private
     
     integer :: nx  = 0 ! Number of parameters
     
     double precision, allocatable :: x(:)
     double precision, allocatable :: prior_param(:,:)
     double precision, allocatable :: perturb_param(:) ! STDEV

     logical :: verb = .false.
     
   contains
     procedure :: set_prior => hyper_model_set_prior
     procedure :: set_perturb => hyper_model_set_perturb
     procedure :: get_prior_param => hyper_model_get_prior_param
     procedure :: get_nx => hyper_model_get_nx
     procedure :: get_x => hyper_model_get_x
     procedure :: generate_model => hyper_model_generate_model
     procedure :: perturb => hyper_model_perturb
     procedure :: display => hyper_model_display
     
  end type hyper_model
  
  interface hyper_model
     module procedure init_hyper_model
  end interface hyper_model
  
contains

  !---------------------------------------------------------------------
  
  type(hyper_model) function init_hyper_model(nx, verb) result(self)
    integer, intent(in) :: nx
    logical, intent(in), optional :: verb
    
    if (present(verb)) then
       self%verb = verb
    end if
    if (self%verb) then
       write(*,'(A)')"<< Initialize hyper model parameters >>"
    end if
    
    ! Get number of model parameters
    self%nx = nx
    allocate(self%x(nx))
    allocate(self%prior_param(nx, 2))
    allocate(self%perturb_param(nx))
    
    if (self%verb) write(*,'(A,I4)')" nx=", nx
    if (self%verb) write(*,*)
    
    return 
  end function init_hyper_model
  
  !---------------------------------------------------------------------
  
  subroutine hyper_model_set_prior(self, iparam, d1, d2)
    class(hyper_model), intent(inout) :: self
    integer, intent(in) :: iparam
    double precision, intent(in) :: d1, d2
    integer :: ierr
    
    if(iparam < 1 .or. iparam > self%nx) then
       write(0,*) "ERROR: out of range (hyper_model_set_prior)"
       write(0,*) "     : iparam=", iparam
       call mpi_finalize(ierr)
       stop
    end if

    if (d1 >= d2) then
       write(0,*)"ERROR: d2 must be > d1 (hyper_model_set_prior)"
       call mpi_finalize(ierr)
       stop
    end if

    self%prior_param(iparam, 1) = d1
    self%prior_param(iparam, 2) = d2
    
    return 
  end subroutine hyper_model_set_prior

  !---------------------------------------------------------------------
  
  subroutine hyper_model_set_perturb(self, iparam, d1)
    class(hyper_model), intent(inout) :: self
    integer, intent(in) :: iparam
    double precision, intent(in) :: d1

    self%perturb_param(iparam) = d1
    
    return 
  end subroutine hyper_model_set_perturb
  
  !---------------------------------------------------------------------

  double precision function hyper_model_get_prior_param(self, &
       & iparam, j) result(prior_param)
    class(hyper_model), intent(in) :: self
    integer, intent(in) :: iparam, j

    prior_param = self%prior_param(iparam, j)
    
    return 
  end function hyper_model_get_prior_param

  !---------------------------------------------------------------------

  integer function hyper_model_get_nx(self) result(nx)
    class(hyper_model), intent(in) :: self
    
    nx = self%nx

    return 
  end function hyper_model_get_nx

  !---------------------------------------------------------------------
  
  double precision function hyper_model_get_x(self, iparam) result(x)
    class(hyper_model), intent(in) :: self
    integer, intent(in) :: iparam
    
    x = self%x(iparam)
    
    return 
  end function hyper_model_get_x

  !---------------------------------------------------------------------

  subroutine hyper_model_generate_model(self)
    class(hyper_model), intent(inout) :: self
    integer :: i

    do i = 1, self%nx
       self%x(i) = self%prior_param(i, 1) + &
               & rand_u() * (self%prior_param(i, 2) - &
               & self%prior_param(i, 1))
    end do

    return 
  end subroutine hyper_model_generate_model
  
  !---------------------------------------------------------------------

  subroutine hyper_model_perturb(self, is_ok)
    class(hyper_model), intent(inout) :: self
    logical, intent(out) :: is_ok
    integer :: iparam

    double precision :: x_new, x_old

    iparam = 1 + int(rand_u() * self%nx)
    x_old = self%x(iparam)
    x_new = x_old + rand_g() * self%perturb_param(iparam)
    self%x(iparam) = x_new
    
    is_ok = .true.
    
    if (   self%x(iparam) < self%prior_param(iparam, 1) .or. &
         & self%x(iparam) > self%prior_param(iparam, 2)) then
       is_ok = .false.
    end if
    
    return 
  end subroutine hyper_model_perturb

  !---------------------------------------------------------------------
  
  subroutine hyper_model_display(self)
    class(hyper_model), intent(inout) :: self
    integer :: i
    
    do i = 1, self%nx
       write(*,*)"------------------------------"
       write(*,*)"Parameter", i
       write(*,*)self%x(i)
    end do

    write(*,*)
    
    return 
  end subroutine hyper_model_display

  !---------------------------------------------------------------------


end module cls_hyper_model
  
