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
module mod_obs
  implicit none
  
  type obs
     private
     character(len=500) :: file_name
     integer :: npts
     integer :: n
     integer :: io
     double precision, allocatable :: data(:)
     double precision :: t_start
     double precision :: t_end 
     double precision :: delta
     double precision :: b
     double precision :: e
     
   contains
     procedure :: set_file_name => obs_set_file_name
     procedure :: get_file_name => obs_get_file_name
     procedure :: get_delta     => obs_get_delta
     procedure :: get_b         => obs_get_b
     procedure :: get_t_start   => obs_get_t_start
     procedure :: get_t_end     => obs_get_t_end
     procedure :: get_n         => obs_get_n
     procedure :: get_data      => obs_get_data
     procedure, private :: read_file => obs_read_file
       
  end type obs
  
  interface obs
     module procedure init_obs
  end interface obs
  
  
contains
  
  !---------------------------------------------------------------------
  
  type(obs) function init_obs(file_name, t_start, t_end)
    
    character(len=*), intent(in) :: file_name
    double precision, intent(in), optional :: t_start, t_end
    integer :: ierr
    real(kind(0e0)) :: tmp
    
    ! check file
    call init_obs%set_file_name(file_name)
    open(newunit = init_obs%io, file = init_obs%file_name, &
         & status = "old", access = "direct", iostat = ierr, recl = 4)
    if (ierr /= 0) then
       write(0,*)"ERROR: cannot open ", trim(init_obs%file_name), "."
       stop
    end if

    ! Read some headers
    read(init_obs%io, rec = 1) tmp
    init_obs%delta = dble(tmp)
    read(init_obs%io, rec = 6) tmp
    init_obs%b = dble(tmp)
    read(init_obs%io, rec = 7) tmp
    init_obs%e = dble(tmp)
    read(init_obs%io, rec = 80) init_obs%npts

    ! Set time window
    if (present(t_start)) then
       init_obs%t_start = t_start
    else
       init_obs%t_start = init_obs%b
    end if
    if (present(t_end)) then
       init_obs%t_end = t_end
    else 
       init_obs%t_end = init_obs%e
    end if

    ! Get data
    call init_obs%read_file()
    
    close(init_obs%io)
    
  end function init_obs

  !---------------------------------------------------------------------
  ! Methods
  !---------------------------------------------------------------------
  subroutine obs_set_file_name(self, file_name)
    class(obs), intent(inout) :: self
    character(len=*), intent(in) :: file_name
    
    self%file_name = file_name
    
    return 
  end subroutine obs_set_file_name
  
  !*********************************************************************

  character(len=500) function obs_get_file_name(self) result(file_name)
    class(obs), intent(in) :: self
    
    file_name = self%file_name
    
    return 
  end function obs_get_file_name
  
  !*********************************************************************
  
  double precision function obs_get_delta(self) result(delta)
    class(obs), intent(in) :: self
    
    delta = self%delta

    return 
  end function obs_get_delta

  !*********************************************************************

  double precision function obs_get_b(self) result(b)
    class(obs), intent(in) :: self
    
    b = self%b

    return 
  end function obs_get_b

  !*********************************************************************

  double precision function obs_get_t_start(self) result(t_start)
    class(obs), intent(in) :: self
    
    t_start = self%t_start

    return 
  end function obs_get_t_start

  !*********************************************************************

  double precision function obs_get_t_end(self) result(t_end)
    class(obs), intent(in) :: self
    
    t_end = self%t_end
    
    return 
  end function obs_get_t_end

  !*********************************************************************
  
  integer function obs_get_n(self) result(n)
    class(obs), intent(in) :: self
    
    n = self%n

    return 
  end function obs_get_n

  !---------------------------------------------------------------------

  function obs_get_data(self) result (x)
    class(obs), intent(in) :: self
    double precision :: c(self%n)
    
    x(1:self%n) = self%data(1:self%n)
    
    return 
  end function obs_get_data

  !---------------------------------------------------------------------

  subroutine obs_read_file(self)
    class(obs), intent(inout) :: self
    real(kind(0e0)) :: tmp
    integer :: i_start, i_end, i
    
    i_start = nint((self%t_start - self%b) / self%delta) + 1
    i_end   = nint((self%t_end - self%b) / self%delta)
    self%n = i_end - i_start + 1
    allocate(self%data(self%n))
    do i = i_start, i_end
       read(self%io, rec = 158+i)tmp
       self%data(i - i_start + 1) = dble(tmp)
    end do
    
    return 
  end subroutine obs_read_file

  
end module mod_obs
