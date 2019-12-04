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
module mod_trans_d_model
  use mod_random
  implicit none 
  
  type trans_d_model
     private
     integer :: k
     integer :: k_min
     integer :: k_max
     
     integer :: n_rx = 0! Number of real parameter 
     integer :: n_ix = 0! Number of integer parameter
     integer :: n_x  = 0 ! Number of all parameter
     
     integer, allocatable :: rx_birth_type(:) ! 1: Uniform, 2: Gaussian
     integer, allocatable :: ix_birth_type(:) ! 1: Uniform, 2: Gaussian
     integer, allocatable :: rx_prior_type(:) ! 1: Uniform, 2: Gaussian
     integer, allocatable :: ix_prior_type(:) ! 1: Uniform, 2: Gaussian
     
     double precision, allocatable :: rx(:, :) ! real Gaussian param.
     integer, allocatable          :: ix(:, :) ! int. Gaussian param.
     
     double precision, allocatable :: rx_birth_param(:,:)
     double precision, allocatable :: ix_birth_param(:,:)
     double precision, allocatable :: rx_prior_param(:,:)
     double precision, allocatable :: ix_prior_param(:,:)
     double precision, allocatable :: rx_perturb_param(:) ! STDEV
     double precision, allocatable :: ix_perturb_param(:) ! STDEV

   contains
     procedure :: set_k => trans_d_model_set_k
     procedure :: set_birth => trans_d_model_set_birth
     procedure :: set_prior => trans_d_model_set_prior
     procedure :: set_perturb => trans_d_model_set_perturb
     procedure :: get_k => trans_d_model_get_k
     procedure :: get_n_x => trans_d_model_get_n_x
     procedure :: get_rx => trans_d_model_get_rx
     procedure :: get_ix => trans_d_model_get_ix
     procedure :: generate_model => trans_d_model_generate_model
     procedure :: birth => trans_d_model_birth
     procedure :: death => trans_d_model_death
     procedure :: perturb => trans_d_model_perturb
     procedure :: perturb_rx => trans_d_model_perturb_rx
     procedure :: perturb_ix => trans_d_model_perturb_ix
     procedure :: display => trans_d_model_display
     procedure :: finish  => trans_d_model_finish
     
  end type trans_d_model
  
  interface trans_d_model
     module procedure init_trans_d_model
  end interface trans_d_model
  
contains

  !---------------------------------------------------------------------
  
  type(trans_d_model) function init_trans_d_model(k_min, k_max, &
       & n_rx, n_ix)
    integer, intent(in) :: k_min, k_max
    integer, intent(in), optional :: n_rx, n_ix
    logical :: is_given

    ! Get dimension
    if (k_max <= k_min) then
       write(0,*)"ERROR: k_max must be > k_min"
       write(0,*)"     : k_min =", k_min, "k_max=", k_max, &
            & "(init_trans_d_model)"
       stop
    end if
    init_trans_d_model%k_min = k_min
    init_trans_d_model%k_max = k_max

    ! Get number of model parameters
    init_trans_d_model%n_rx = 0
    init_trans_d_model%n_ix = 0
    is_given = .false.
    if (present(n_rx)) then
       init_trans_d_model%n_rx = n_rx
       allocate(init_trans_d_model%rx(k_max, n_rx))
       allocate(init_trans_d_model%rx_birth_type(n_rx))
       allocate(init_trans_d_model%rx_prior_type(n_rx))
       allocate(init_trans_d_model%rx_birth_param(n_rx, 2))
       allocate(init_trans_d_model%rx_prior_param(n_rx, 2))
       allocate(init_trans_d_model%rx_perturb_param(n_rx))
       is_given = .true.
    end if
    if (present(n_ix)) then
       init_trans_d_model%n_ix = n_ix
       allocate(init_trans_d_model%ix(k_max, n_ix))
       allocate(init_trans_d_model%ix_birth_type(n_ix))
       allocate(init_trans_d_model%ix_prior_type(n_ix))
       allocate(init_trans_d_model%ix_birth_param(n_ix, 2))
       allocate(init_trans_d_model%ix_prior_param(n_ix, 2))
       allocate(init_trans_d_model%ix_perturb_param(n_ix))
       is_given = .true.
    end if
    if (.not. is_given) then
       write(0,*)"ERROR: number of model parameters must be given"
       write(0,*)"     : (init_trans_d_model)"
       stop
    end if
    
    init_trans_d_model%n_x = &
         & init_trans_d_model%n_rx + &
         & init_trans_d_model%n_ix

    
    return 
  end function init_trans_d_model
  
  !---------------------------------------------------------------------
  
  subroutine trans_d_model_set_k(self, k)
    class(trans_d_model), intent(inout) :: self
    integer, intent(in) :: k
    
    if (k < self%k_min .or. k >= self%k_max) then
       write(0,*)"ERROR: k must be k_min <= k < k_max ", &
            & "(trans_d_model_set_k)"
       write(0,*)"     : k=", k, " k_min=", self%k_min, &
            & " k_max=", self%k_max
    end if

    self%k = k
    
    return 
  end subroutine trans_d_model_set_k

  !---------------------------------------------------------------------
  
  subroutine trans_d_model_set_birth(self, iparam, itype, d1, d2)
    class(trans_d_model), intent(inout) :: self
    integer, intent(in) :: iparam, itype
    double precision, intent(in) :: d1, d2

    if(iparam < 1 .or. iparam > self%k_max) then
       write(0,*) "ERROR: out of range (trans_d_model_set_birth)"
       write(0,*) "     : iparam=", iparam
       stop
    end if

    if (itype == 1 .and. d1 >= d2) then
       write(0,*)"ERROR: d2 must be > d1 (trans_d_model_set_birth)"
       stop
    end if

    if (itype /= 1 .and. itype /= 2) then
       write(0, *)"ERROR: itype (1: uniform, 2: Gaussian)"
       write(0, *)"     : (trans_d_model_set_birth)"
       stop
    end if
    
    if (iparam <= self%n_rx) then
       self%rx_birth_type(iparam) = itype
       self%rx_birth_param(iparam, 1) = d1
       self%rx_birth_param(iparam, 2) = d2
    else
       self%ix_birth_type(iparam - self%n_rx) = itype
       self%ix_birth_param(iparam - self%n_rx, 1) = d1
       self%ix_birth_param(iparam - self%n_rx, 2) = d2
    end if
    
    return 
  end subroutine trans_d_model_set_birth

  !---------------------------------------------------------------------
  
  subroutine trans_d_model_set_prior(self, iparam, itype, d1, d2)
    class(trans_d_model), intent(inout) :: self
    integer, intent(in) :: iparam, itype
    double precision, intent(in) :: d1, d2

    if(iparam < 1 .or. iparam > self%k_max) then
       write(0,*) "ERROR: out of range (trans_d_model_set_prior)"
       write(0,*) "     : iparam=", iparam
       stop
    end if

    if (itype == 1 .and. d1 >= d2) then
       write(0,*)"ERROR: d2 must be > d1 (trans_d_model_set_prior)"
       stop
    end if

    if (itype /= 1 .and. itype /= 2) then
       write(0, *)"ERROR: itype (1: uniform, 2: Gaussian)"
       write(0, *)"     : (trans_d_model_set_prior)"
       stop
    end if
    
    if (iparam <= self%n_rx) then
       self%rx_prior_type(iparam) = itype
       self%rx_prior_param(iparam, 1) = d1
       self%rx_prior_param(iparam, 2) = d2
    else
       self%ix_prior_type(iparam - self%n_rx) = itype
       self%ix_prior_param(iparam - self%n_rx, 1) = d1
       self%ix_prior_param(iparam - self%n_rx, 2) = d2
    end if
    
    return 
  end subroutine trans_d_model_set_prior

  !---------------------------------------------------------------------
  
  subroutine trans_d_model_set_perturb(self, iparam, d1)
    class(trans_d_model), intent(inout) :: self
    integer, intent(in) :: iparam
    double precision, intent(in) :: d1

    if(iparam < 1 .or. iparam > self%k_max) then
       write(0,*) "ERROR: out of range (trans_d_model_set_prior)"
       write(0,*) "     : iparam=", iparam
       stop
    end if
    
    if (iparam <= self%n_rx) then
       self%rx_perturb_param(iparam) = d1
    else
       self%ix_perturb_param(iparam - self%n_rx) = d1
    end if
    
    return 
  end subroutine trans_d_model_set_perturb
  
  !---------------------------------------------------------------------

  integer function trans_d_model_get_k(self) result(k)
    class(trans_d_model), intent(in) :: self
    
    k = self%k

    return 
  end function trans_d_model_get_k

  !---------------------------------------------------------------------

  integer function trans_d_model_get_n_x(self) result(n_x)
    class(trans_d_model), intent(in) :: self
    
    n_x = self%n_x

    return 
  end function trans_d_model_get_n_x

  !---------------------------------------------------------------------
  
  function trans_d_model_get_rx(self, iparam) result(rx)
    class(trans_d_model), intent(in) :: self
    integer, intent(in) :: iparam
    double precision :: rx(self%k_max)
    
    rx = self%rx(1:self%k_max, iparam)
    
    return 
  end function trans_d_model_get_rx

  !---------------------------------------------------------------------

  function trans_d_model_get_ix(self, iparam) result(ix)
    class(trans_d_model), intent(in) :: self
    integer, intent(in) :: iparam
    integer :: ix(self%k_max)
    
    ix = self%ix(1:self%k_max, iparam)
    
    return 
  end function trans_d_model_get_ix

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

    if (self%k < self%k_max) then
       is_ok = .true.
    else
       is_ok = .false.
       return
    end if
    k2 = self%k + 1
    do i = 1, self%n_rx
       if (self%rx_birth_type(i) == 1) then
          self%rx(k2, i) = &
               & self%rx_birth_param(i, 1) + &
               & rand_u() * (self%rx_birth_param(i, 2) - &
               & self%rx_birth_param(i, 1))

       else if (self%rx_birth_type(i) == 2) then
          self%rx(k2, i) = &
               & self%rx_birth_param(i, 1) + &
               & rand_g() * self%rx_birth_param(i, 2)
       end if
       ! check prior bounds
       if (self%rx_prior_type(i) == 1) then
          if(self%rx(k2, i) < self%rx_prior_param(i, 1) .or. &
               & self%rx(k2, i) > self%rx_prior_param(i, 2)) then
             is_ok = .false.
             return
          end if
       end if
    end do
    
    do i = 1, self%n_ix
       if (self%ix_birth_type(i) == 1) then
          self%ix(k2, i) = &
               & nint(self%ix_birth_param(i, 1)) + &
               & int(rand_u() * &
               & (nint(self%ix_birth_param(i, 2)) - &
               & nint(self%ix_birth_param(i, 1)))&
               &)
       else if (self%ix_birth_type(i) == 2) then
          self%ix(k2, i) = &
               & nint( &
               & self%ix_birth_param(i, 1) + &
               & rand_g() * self%ix_birth_param(i, 2))
       end if
       ! Check prior bounds
       if (self%ix_prior_type(i) == 1) then
          if(self%ix(k2, i) < nint(self%ix_prior_param(i, 1)) .or. &
               & self%ix(k2, i) >= nint(self%ix_prior_param(i, 2))) then
             is_ok = .false.
             return 
          end if
       end if
    end do

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
    

    do i = 1, self%n_rx
       do j = k_target, self%k
          self%rx(j,  i) = self%rx(j + 1, i) 
       end do
    end do
    do i = 1, self%n_ix
       do j = k_target, self%k
          self%ix(j,  i) = self%ix(j + 1, i) 
       end do
    end do
    
    return 
  end subroutine trans_d_model_death

  !---------------------------------------------------------------------
  
  subroutine trans_d_model_perturb(self, iparam, is_ok)
    class(trans_d_model), intent(inout) :: self
    integer, intent(in) :: iparam
    logical, intent(out) :: is_ok
    
    if (iparam <= self%n_rx) then
       call self%perturb_rx(iparam, &
            & self%rx_perturb_param(iparam), is_ok)
    else
       call self%perturb_ix(iparam - self%n_rx, &
            & self%ix_perturb_param(iparam), is_ok)
    end if
    
    return 
  end subroutine trans_d_model_perturb
  
  !---------------------------------------------------------------------

  subroutine trans_d_model_perturb_rx(self, iparam, dev, is_ok)
    class(trans_d_model), intent(inout) :: self
    integer, intent(in) :: iparam
    double precision, intent(in) :: dev
    logical, intent(out) :: is_ok
    integer :: k_target

    k_target = self%k_min + int(rand_u()*((self%k + 1) - self%k_min))
    
    self%rx(k_target, iparam) = self%rx(k_target, iparam) + &
         & rand_g() * dev

    is_ok = .true.
    if (self%rx_prior_type(iparam) == 1) then
       if (self%rx(k_target, iparam) < &
            & self%rx_prior_param(iparam, 1) .or. &
            & self%rx(k_target, iparam) > &
            & self%rx_prior_param(iparam, 2)) then
          is_ok = .false.
       end if
    end if
    
    
    return 
  end subroutine trans_d_model_perturb_rx

  !---------------------------------------------------------------------

  subroutine trans_d_model_perturb_ix(self, iparam, dev, is_ok)
    class(trans_d_model), intent(inout) :: self
    integer, intent(in) :: iparam
    double precision, intent(in) :: dev
    logical, intent(out) :: is_ok
    integer :: k_target
    
    k_target = self%k_min + int(rand_u()*((self%k + 1) - self%k_min))
    
    self%ix(k_target, iparam) = self%ix(k_target, iparam) + &
         & nint(rand_g() * dev)
    
    is_ok = .true.
    if (self%ix_prior_type(iparam) == 1) then
       if (self%ix(k_target, iparam) < &
            & self%ix_prior_param(iparam, 1) .or. &
            & self%ix(k_target, iparam) > &
            & self%ix_prior_param(iparam, 2)) then
          is_ok = .false.
       end if
    end if
    
    return 
  end subroutine trans_d_model_perturb_ix



  !---------------------------------------------------------------------
  
  subroutine trans_d_model_display(self)
    class(trans_d_model), intent(inout) :: self
    integer :: i, j
    write(*,*)
    write(*,*)"k = ", self%k
    
    do i = 1, self%n_rx
       write(*,*)"------------------------------"
       write(*,*)"Real parameter", i
       do j = 1, self%k
          write(*,*)j, self%rx(j, i)
       end do
    end do

    do i = 1, self%n_ix
       write(*,*)"------------------------------"
       write(*,*)"Integer prameter", i
       do j = 1, self%k
          write(*,*)j, self%ix(j, i)
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
    
    if (self%n_rx > 0) then
       deallocate(self%rx)
       deallocate(self%rx_birth_type)
       deallocate(self%rx_prior_type)
       deallocate(self%rx_birth_param)
       deallocate(self%rx_prior_param)
       deallocate(self%rx_perturb_param)
    end if
    if (self%n_ix > 0) then
       deallocate(self%ix)
       deallocate(self%ix_birth_type)
       deallocate(self%ix_prior_type)
       deallocate(self%ix_birth_param)
       deallocate(self%ix_prior_param)
       deallocate(self%ix_perturb_param)
    end if
    
    return
  end subroutine trans_d_model_finish

end module mod_trans_d_model
  
