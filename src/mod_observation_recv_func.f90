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
module mod_observation_recv_func
  use mod_line_text, only: line_max, line_text
  implicit none
  
  type observation_recv_func
     private
     
     ! receiver function
     integer :: n_rf
     integer, allocatable :: n_smp(:)
     double precision, allocatable :: a_gauss(:)
     double precision, allocatable :: rayp(:)
     double precision, allocatable :: delta(:)
     double precision, allocatable :: sigma_min(:)
     double precision, allocatable :: sigma_max(:)
     double precision, allocatable :: rf_data(:,:)
     double precision, allocatable :: t_start(:)
     double precision, allocatable :: t_end(:)
     character(len=1), allocatable :: phase(:)
     logical, allocatable          :: deconv_flag(:)
     logical, allocatable          :: correct_amp(:)
     
     
   contains
     procedure :: read_sac => observation_recv_func_read_sac
     procedure :: get_rf_data => observation_recv_func_get_rf_data
     procedure :: get_n_rf => observation_recv_func_get_n_rf
     procedure :: get_n_smp => observation_recv_func_get_n_smp
     procedure :: get_t_start => observation_recv_func_get_t_start
     procedure :: get_delta => observation_recv_func_get_delta
     procedure :: get_rayp => observation_recv_func_get_rayp
     procedure :: get_a_gauss => observation_recv_func_get_a_gauss
     procedure :: get_phase => observation_recv_func_get_phase
     procedure :: get_deconv_flag => &
          & observation_recv_func_get_deconv_flag
     procedure :: get_correct_amp => &
          & observation_recv_func_get_correct_amp
  end type observation_recv_func
  
  interface observation_recv_func
     module procedure :: init_observation_recv_func
  end interface observation_recv_func

contains

  !--------------------------------------------------------------------- 
  
  type(observation_recv_func) function &
     & init_observation_recv_func(recv_func_in) &
     & result(self)
    character(len=*), intent(in) :: recv_func_in
    type(line_text) :: lt
    character(line_max) :: line
    character(line_max), allocatable :: filename(:)
    integer :: io, ierr, i
    

    write(*,*)"Reading observed receiver functions from ", &
         & trim(recv_func_in)

    open(newunit = io, file = recv_func_in, status = "old", &
         & iostat = ierr)
    if (ierr /= 0) then
       write(0,*)"ERROR: cannot open ", trim(recv_func_in)
       write(0,*)"     : (init_observation_recv_func)"
       stop
    end if
    
    
    ! Read observation file
    do
       read(io, '(a)')line
       lt = line_text(line, ignore_space=.false.)
       line = lt%get_line()
       if (len_trim(line) /= 0) then
          read(line, *) self%n_rf
          write(*,*)"n_rf =", self%n_rf
          exit
       end if
    end do
    
    ! allocate
    allocate(self%a_gauss(self%n_rf), self%rayp(self%n_rf), &
         & self%delta(self%n_rf))
    allocate(self%sigma_min(self%n_rf), self%sigma_max(self%n_rf))
    allocate(self%t_start(self%n_rf), self%t_end(self%n_rf))
    allocate(self%n_smp(self%n_rf))
    allocate(self%phase(self%n_rf))
    allocate(self%deconv_flag(self%n_rf))
    allocate(self%correct_amp(self%n_rf))
    allocate(filename(self%n_rf))


    do i = 1, self%n_rf
       ! get file name
       do
          read(io, '(a)')line
          lt = line_text(line, ignore_space=.false.)
          line = lt%get_line()
          if (len_trim(line) /= 0) then
             read(line, *) filename(i)
             write(*,*)"filename = ", trim(filename(i))
             exit
          end if
       end do
       ! a_gauss, rayp, delta
       do
          read(io, '(a)')line
          lt = line_text(line, ignore_space=.false.)
          line = lt%get_line()
          if (len_trim(line) /= 0) then
             read(line, *) self%a_gauss(i), &
                  & self%rayp(i), self%delta(i)
             write(*,*)self%a_gauss(i), self%rayp(i), self%delta(i)
             exit
          end if
       end do
       ! phase
       do
          read(io, '(a)')line
          lt = line_text(line, ignore_space=.false.)
          line = lt%get_line()
          if (len_trim(line) /= 0) then
             read(line, *) self%phase(i)
             write(*,*)self%phase(i)
             exit
          end if
       end do
       ! t_start, t_end
       do
          read(io, '(a)')line
          lt = line_text(line, ignore_space=.false.)
          line = lt%get_line()
          if (len_trim(line) /= 0) then
             read(line, *) self%t_start(i), self%t_end(i)
             write(*,*)self%t_start(i), self%t_end(i)
             exit
          end if
       end do
       !sigma_min, sigma_max
       do
          read(io, '(a)')line
          lt = line_text(line, ignore_space=.false.)
          line = lt%get_line()
          if (len_trim(line) /= 0) then
             read(line, *) self%sigma_min(i), self%sigma_max(i)
             write(*,*)self%sigma_min(i), self%sigma_max(i)
             exit
          end if
       end do
       ! deconv_flag, correct_amp
       do
          read(io, '(a)')line
          lt = line_text(line, ignore_space=.false.)
          line = lt%get_line()
          if (len_trim(line) /= 0) then
             read(line, *) self%deconv_flag(i), self%correct_amp(i)
             write(*,*)self%deconv_flag(i), self%correct_amp(i)
             exit
          end if
       end do
       
    end do

    ! allocate data array
    do concurrent (i=1:self%n_rf)
       self%n_smp(:) &
             & = nint((self%t_end(:) - self%t_start(:)) / self%delta(:))
    end do
    allocate(self%rf_data(maxval(self%n_smp(:)), self%n_rf))

    ! Read data
    do i = 1, self%n_rf
       call self%read_sac(filename(i), i)
    end do

    return 
  end function init_observation_recv_func
  
  !---------------------------------------------------------------------

  subroutine observation_recv_func_read_sac(self, filename, i)
    class(observation_recv_func), intent(inout) :: self
    character(*), intent(in) :: filename
    integer, intent(in) :: i
    integer :: ierr
    real :: tmp
    double precision :: b, delta
    integer :: i_start, i_end, j, io
    
    open(newunit = io, file = filename, access = 'direct', &
         & iostat = ierr, recl = 4)
    if (ierr /= 0) then
       write(0,*)"ERROR: cannot open ", trim(filename)
       stop
    end if

    
    read(io, rec = 1) tmp    
    delta = tmp
    if (abs(delta - self%delta(i)) > 1.0e-5) then
       write(0,*)"ERROR: sampling interval is not correct"
       write(0,*)" Given: ", self%delta(i)
       write(0,*)" Data : ", delta
       stop
    end if
    read(io, rec = 6) tmp
    b = tmp

    i_start = nint((self%t_start(i) - b) / self%delta(i)) + 1
    i_end   = i_start + self%n_smp(i) - 1
    do j = i_start, i_end
       read(io, rec = 158 + j)tmp
       self%rf_data(j - i_start + 1, i) = dble(tmp)
    end do
    close(io)
    
    return 
  end subroutine observation_recv_func_read_sac

  !---------------------------------------------------------------------

  function observation_recv_func_get_rf_data(self, j, i) result(rf_data)
    class(observation_recv_func), intent(in) :: self
    integer, intent(in) :: i, j
    double precision :: rf_data

    rf_data = self%rf_data(j, i)
    
    return 
  end function observation_recv_func_get_rf_data

  !---------------------------------------------------------------------

  integer function observation_recv_func_get_n_rf(self) result(n_rf)
    class(observation_recv_func), intent(in) :: self
    
    n_rf = self%n_rf
    
    return 
  end function observation_recv_func_get_n_rf
    
  !---------------------------------------------------------------------

  integer function observation_recv_func_get_n_smp(self, i) result(n_smp)
    class(observation_recv_func), intent(in) :: self
    integer, intent(in) :: i
    
    n_smp = self%n_smp(i)
    
    return 
  end function observation_recv_func_get_n_smp
  
  !---------------------------------------------------------------------

  double precision function observation_recv_func_get_t_start(self, i) &
       & result(t_start)
    class(observation_recv_func), intent(in) :: self
    integer, intent(in) :: i
    
    t_start = self%t_start(i)
    
    return 
  end function observation_recv_func_get_t_start
  
  !---------------------------------------------------------------------

  double precision function observation_recv_func_get_delta(self, i) &
       & result(delta)
    class(observation_recv_func), intent(in) :: self
    integer, intent(in) :: i
    
    delta = self%delta(i)
    
    return 
  end function observation_recv_func_get_delta
  
  !---------------------------------------------------------------------

  double precision function observation_recv_func_get_rayp(self, i) &
       & result(rayp)
    class(observation_recv_func), intent(in) :: self
    integer, intent(in) :: i
    
    rayp = self%rayp(i)
    
    return 
  end function observation_recv_func_get_rayp
  
  !---------------------------------------------------------------------

  double precision function observation_recv_func_get_a_gauss(self, i) &
       & result(a_gauss)
    class(observation_recv_func), intent(in) :: self
    integer, intent(in) :: i
    
    a_gauss = self%a_gauss(i)
    
    return 
  end function observation_recv_func_get_a_gauss
  
  !---------------------------------------------------------------------

  character(1) function observation_recv_func_get_phase(self, i) &
       & result(phase)
    class(observation_recv_func), intent(in) :: self
    integer, intent(in) :: i
    
    phase = self%phase(i)
    
    return 
  end function observation_recv_func_get_phase
  
  !---------------------------------------------------------------------

  logical function observation_recv_func_get_deconv_flag(self, i) &
       & result(deconv_flag)
    class(observation_recv_func), intent(in) :: self
    integer, intent(in) :: i
    
    deconv_flag = self%deconv_flag(i)
    
    return 
  end function observation_recv_func_get_deconv_flag
  
  !---------------------------------------------------------------------

  logical function observation_recv_func_get_correct_amp(self, i) &
       & result(correct_amp)
    class(observation_recv_func), intent(in) :: self
    integer, intent(in) :: i
    
    correct_amp = self%correct_amp(i)
    
    return 
  end function observation_recv_func_get_correct_amp
  
  !---------------------------------------------------------------------
  
end module mod_observation_recv_func
