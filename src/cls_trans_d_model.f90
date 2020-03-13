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
module cls_trans_d_model
  use mod_random
  implicit none 
  
  type trans_d_model
     private
     integer :: k
     integer :: k_min
     integer :: k_max
     
     integer :: nx  = 0 ! Number of parameters per layer
     
     integer, allocatable :: birth_type(:) ! 1: Uniform, 2: Gaussian
     integer, allocatable :: prior_type(:) ! 1: Uniform, 2: Gaussian
     
     double precision, allocatable :: x(:, :)
     
     double precision, allocatable :: birth_param(:,:)
     double precision, allocatable :: prior_param(:,:)
     double precision, allocatable :: perturb_param(:) ! STDEV

     logical :: verb = .false.
     
   contains
     procedure :: set_k => trans_d_model_set_k
     procedure :: set_birth => trans_d_model_set_birth
     procedure :: set_prior => trans_d_model_set_prior
     procedure :: set_perturb => trans_d_model_set_perturb
     procedure :: get_k => trans_d_model_get_k
     procedure :: get_k_max => trans_d_model_get_k_max
     procedure :: get_nx => trans_d_model_get_nx
     procedure :: get_x => trans_d_model_get_x
     procedure :: get_x_single => trans_d_model_get_x_single
     procedure :: generate_model => trans_d_model_generate_model
     procedure :: birth => trans_d_model_birth
     procedure :: death => trans_d_model_death
     procedure :: perturb => trans_d_model_perturb
     procedure :: display => trans_d_model_display
     procedure :: finish  => trans_d_model_finish
     
  end type trans_d_model
  
  interface trans_d_model
     module procedure init_trans_d_model
  end interface trans_d_model
  
contains

  !---------------------------------------------------------------------
  
  type(trans_d_model) function init_trans_d_model(k_min, k_max, &
       & nx, verb) result(self)
    integer, intent(in) :: k_min, k_max
    integer, intent(in) :: nx
    logical, intent(in), optional :: verb
    integer :: ierr
    
    if (present(verb)) then
       self%verb = verb
    end if
    if (self%verb) then
       write(*,'(A)')"<< Initialize trans-D model parameters >>"
    end if
    
    ! Get dimension
    if (k_max <= k_min) then
       write(0,*)"ERROR: k_max must be > k_min"
       write(0,*)"     : k_min =", k_min, "k_max=", k_max
       call mpi_finalize(ierr)
       stop
    end if
    self%k_min = k_min
    self%k_max = k_max
    if (self%verb) write(*,'(A,I4)')" k_min=", k_min
    if (self%verb) write(*,'(A,I4)')" k_max=", k_max


    ! Get number of model parameters
    self%nx = nx
    allocate(self%x(k_max, nx))
    allocate(self%birth_type(nx))
    allocate(self%prior_type(nx))
    allocate(self%birth_param(nx, 2))
    allocate(self%prior_param(nx, 2))
    allocate(self%perturb_param(nx))
    
    if (self%verb) write(*,'(A,I4)')" nx=", nx
    if (self%verb) write(*,*)
    return 
  end function init_trans_d_model
  
  !---------------------------------------------------------------------
  
  subroutine trans_d_model_set_k(self, k)
    class(trans_d_model), intent(inout) :: self
    integer, intent(in) :: k
    integer :: ierr
    
    if (k < self%k_min .or. k >= self%k_max) then
       write(0,*)"ERROR: k must be k_min <= k < k_max ", &
            & "(trans_d_model_set_k)"
       write(0,*)"     : k=", k, " k_min=", self%k_min, &
            & " k_max=", self%k_max
       call mpi_finalize(ierr)
       stop
    end if

    self%k = k
    
    return 
  end subroutine trans_d_model_set_k

  !---------------------------------------------------------------------
  
  subroutine trans_d_model_set_birth(self, iparam, itype, d1, d2)
    class(trans_d_model), intent(inout) :: self
    integer, intent(in) :: iparam, itype
    double precision, intent(in) :: d1, d2
    integer :: ierr

    if(iparam < 1 .or. iparam > self%k_max) then
       write(0,*) "ERROR: out of range (trans_d_model_set_birth)"
       call mpi_finalize(ierr)
       stop
    end if

    if (itype == 1 .and. d1 >= d2) then
       write(0,*)"ERROR: d2 must be > d1 (trans_d_model_set_birth)"
       call mpi_finalize(ierr)
       stop
    end if

    if (itype /= 1 .and. itype /= 2) then
       write(0, *)"ERROR: itype (1: uniform, 2: Gaussian)"
       write(0, *)"     : (trans_d_model_set_birth)"
       call mpi_finalize(ierr)
       stop
    end if
    
    self%birth_type(iparam) = itype
    self%birth_param(iparam, 1) = d1
    self%birth_param(iparam, 2) = d2
    
    return 
  end subroutine trans_d_model_set_birth

  !---------------------------------------------------------------------
  
  subroutine trans_d_model_set_prior(self, iparam, itype, d1, d2)
    class(trans_d_model), intent(inout) :: self
    integer, intent(in) :: iparam, itype
    double precision, intent(in) :: d1, d2
    integer :: ierr
    
    if(iparam < 1 .or. iparam > self%k_max) then
       write(0,*) "ERROR: out of range (trans_d_model_set_prior)"
       write(0,*) "     : iparam=", iparam
       call mpi_finalize(ierr)
       stop
    end if

    if (itype == 1 .and. d1 >= d2) then
       write(0,*)"ERROR: d2 must be > d1 (trans_d_model_set_prior)"
       call mpi_finalize(ierr)
       stop
    end if

    if (itype /= 1 .and. itype /= 2) then
       write(0, *)"ERROR: itype (1: uniform, 2: Gaussian)"
       write(0, *)"     : (trans_d_model_set_prior)"
       call mpi_finalize(ierr)
       stop
    end if
    
    self%prior_type(iparam) = itype
    self%prior_param(iparam, 1) = d1
    self%prior_param(iparam, 2) = d2
    
    return 
  end subroutine trans_d_model_set_prior

  !---------------------------------------------------------------------
  
  subroutine trans_d_model_set_perturb(self, iparam, d1)
    class(trans_d_model), intent(inout) :: self
    integer, intent(in) :: iparam
    double precision, intent(in) :: d1
    integer :: ierr

    if(iparam < 1 .or. iparam > self%k_max) then
       write(0,*) "ERROR: out of range (trans_d_model_set_prior)"
       write(0,*) "     : iparam=", iparam
       call mpi_finalize(ierr)
       stop
    end if
    
    self%perturb_param(iparam) = d1
    
    return 
  end subroutine trans_d_model_set_perturb
  
  !---------------------------------------------------------------------

  integer function trans_d_model_get_k(self) result(k)
    class(trans_d_model), intent(in) :: self
    
    k = self%k

    return 
  end function trans_d_model_get_k

  !---------------------------------------------------------------------
  
  integer function trans_d_model_get_k_max(self) result(k_max)
    class(trans_d_model), intent(in) :: self
    
    k_max = self%k_max

    return 
  end function trans_d_model_get_k_max

  !---------------------------------------------------------------------

  integer function trans_d_model_get_nx(self) result(nx)
    class(trans_d_model), intent(in) :: self
    
    nx = self%nx

    return 
  end function trans_d_model_get_nx

  !---------------------------------------------------------------------
  
  function trans_d_model_get_x(self, iparam) result(x)
    class(trans_d_model), intent(in) :: self
    integer, intent(in) :: iparam
    double precision :: x(self%k_max)
    
    x = self%x(1:self%k_max, iparam)
    
    return 
  end function trans_d_model_get_x

  !---------------------------------------------------------------------
  
  double precision function trans_d_model_get_x_single(self, k, &
       & iparam) result(x)
    class(trans_d_model), intent(in) :: self
    integer, intent(in) :: k, iparam
    
    x = self%x(k, iparam)

    return 
  end function trans_d_model_get_x_single

  !---------------------------------------------------------------------

  subroutine trans_d_model_generate_model(self)
    class(trans_d_model), intent(inout) :: self
    integer :: i, k
    logical :: is_ok
    
    ! select k
    k = self%k_min + int(rand_u() * (self%k_max - self%k_min))

    ! generate model parameters
    self%k = 0
    do i = 1, k
       is_ok = .false.
       do 
          call self%birth(is_ok)
          if (is_ok) exit
       end do
    end do

    return 
  end subroutine trans_d_model_generate_model
  
  !---------------------------------------------------------------------

  subroutine trans_d_model_birth(self, is_ok)
    class(trans_d_model), intent(inout) :: self
    logical, intent(out) :: is_ok
    integer :: i, k2

    if (self%k < self%k_max - 1) then
       is_ok = .true.
    else
       is_ok = .false.
       return
    end if
    k2 = self%k + 1
    do i = 1, self%nx
       if (self%birth_type(i) == 1) then
          self%x(k2, i) = &
               & self%birth_param(i, 1) + &
               & rand_u() * (self%birth_param(i, 2) - &
               & self%birth_param(i, 1))

       else if (self%birth_type(i) == 2) then
          self%x(k2, i) = &
               & self%birth_param(i, 1) + &
               & rand_g() * self%birth_param(i, 2)
       end if
       ! check prior bounds
       if (self%prior_type(i) == 1) then
          if(self%x(k2, i) < self%prior_param(i, 1) .or. &
               & self%x(k2, i) > self%prior_param(i, 2)) then
             is_ok = .false.
             return
          end if
       end if
    end do
    
    ! Final check
    if (is_ok) then
       self%k = k2
    end if
    
    return 
  end subroutine trans_d_model_birth

!---------------------------------------------------------------------

  subroutine trans_d_model_death(self, is_ok)
    class(trans_d_model), intent(inout) :: self
    logical, intent(out) :: is_ok
    integer :: i, j, k_target

    k_target = self%k_min + int(rand_u()*((self%k + 1) - self%k_min))

    
    if (self%k > self%k_min) then
       self%k = self%k - 1
       is_ok = .true.
    else
       is_ok = .false.
       return
    end if
    

    do i = 1, self%nx
       do j = k_target, self%k
          self%x(j,  i) = self%x(j + 1, i) 
       end do
    end do
    
    return 
  end subroutine trans_d_model_death

  !---------------------------------------------------------------------
  
  !subroutine trans_d_model_perturb(self, iparam, is_ok, log_prior_ratio)
  !  class(trans_d_model), intent(inout) :: self
  !  integer, intent(in) :: iparam
  !  logical, intent(out) :: is_ok
  !  double precision, intent(out) :: log_prior_ratio
  !  
  !  if (iparam <= self%n_rx) then
  !     call self%perturb_rx(iparam, &
  !          & self%rx_perturb_param(iparam), is_ok, log_prior_ratio)
  !  else
  !     call self%perturb_ix(iparam - self%n_rx, &
  !          & self%ix_perturb_param(iparam), is_ok, log_prior_ratio)
  !  end if
  !  
  !  return 
  !end subroutine trans_d_model_perturb
  
  !---------------------------------------------------------------------

  subroutine trans_d_model_perturb(self, iparam, is_ok, &
       & log_prior_ratio)
    class(trans_d_model), intent(inout) :: self
    integer, intent(in) :: iparam
    logical, intent(out) :: is_ok
    double precision, intent(out) :: log_prior_ratio
    double precision :: dev
    integer :: k_target
    double precision :: x_new, x_old, mu, sig

    dev = self%perturb_param(iparam)
    k_target = self%k_min + int(rand_u()*((self%k + 1) - self%k_min))

    x_old = self%x(k_target, iparam)
    x_new = x_old + rand_g() * dev
    self%x(k_target, iparam) = x_new
    
    is_ok = .true.
    if (self%prior_type(iparam) == 1) then
       if (self%x(k_target, iparam) < &
            & self%prior_param(iparam, 1) .or. &
            & self%x(k_target, iparam) > &
            & self%prior_param(iparam, 2)) then
          is_ok = .false.
       end if
       log_prior_ratio = 0.d0
    else if (self%prior_type(iparam) == 2) then
       mu = self%prior_param(iparam, 1)
       sig = self%prior_param(iparam, 2)
       log_prior_ratio = -((x_new - mu)**2 - (x_old - mu)**2) / &
            & (2.d0 * sig * sig)
    end if
    
    
    return 
  end subroutine trans_d_model_perturb

  !---------------------------------------------------------------------

  !subroutine trans_d_model_perturb_ix(self, iparam, dev, is_ok, &
  !     & log_prior_ratio)
  !  class(trans_d_model), intent(inout) :: self
  !  integer, intent(in) :: iparam
  !  double precision, intent(in) :: dev
  !  logical, intent(out) :: is_ok
  !  double precision, intent(out) :: log_prior_ratio
  !  integer :: k_target
  !  double precision :: mu, sig
  !  integer :: x_new, x_old
  !  
  !  k_target = self%k_min + int(rand_u()*((self%k + 1) - self%k_min))
  !  
  !  x_old = self%ix(k_target, iparam)
  !  x_new = x_old + nint(rand_g() * dev)
  !  self%ix(k_target, iparam) = x_new
  !  
  !  is_ok = .true.
  !  if (self%ix_prior_type(iparam) == 1) then
  !     if (self%ix(k_target, iparam) < &
  !          & self%ix_prior_param(iparam, 1) .or. &
  !          & self%ix(k_target, iparam) > &
  !          & self%ix_prior_param(iparam, 2)) then
  !        is_ok = .false.
  !     end if
  !     log_prior_ratio = 0.d0
  !  else if (self%ix_prior_type(iparam) == 2) then
  !     mu = self%ix_prior_param(iparam, 1)
  !     sig = self%ix_prior_param(iparam, 2)
  !     log_prior_ratio = -((x_new - mu)**2 - (x_old - mu)**2) / &
  !          & (2.d0 * sig * sig)
  !  end if
  !  
  !  return 
  !end subroutine trans_d_model_perturb_ix



  !---------------------------------------------------------------------
  
  subroutine trans_d_model_display(self)
    class(trans_d_model), intent(inout) :: self
    integer :: i, j
    write(*,*)
    write(*,*)"k = ", self%k
    
    do i = 1, self%nx
       write(*,*)"------------------------------"
       write(*,*)"Parameter", i
       do j = 1, self%k
          write(*,*)j, self%x(j, i)
       end do
    end do

    write(*,*)
    
    return 
  end subroutine trans_d_model_display

  !---------------------------------------------------------------------

  subroutine trans_d_model_finish(self)
    class(trans_d_model), intent(inout) :: self
    
    self%k = 0
    self%k_min = 0
    self%k_max = 0
    
    deallocate(self%x)
    deallocate(self%birth_type)
    deallocate(self%prior_type)
    deallocate(self%birth_param)
    deallocate(self%prior_param)
    deallocate(self%perturb_param)
    
    
    return
  end subroutine trans_d_model_finish

end module cls_trans_d_model
  
