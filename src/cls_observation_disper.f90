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
module cls_observation_disper
  use cls_line_text
  implicit none 
  
  type observation_disper
     private
     ! surface wave
     integer :: n_disp
     integer, allocatable :: nf(:)
     double precision, allocatable :: fmin(:), df(:), fmax(:)
     double precision, allocatable :: c(:,:) ! phase velocity
     double precision, allocatable :: u(:,:) ! group velocity
     double precision, allocatable :: sig_c(:) ! uncertainties in c
     double precision, allocatable :: sig_u(:) ! uncertainties in u
     double precision, allocatable :: cmin(:), cmax(:), dc(:)
     double precision, allocatable :: sig_c_min(:), sig_u_min(:)
     double precision, allocatable :: sig_c_max(:), sig_u_max(:)
     double precision, allocatable :: dev_sig_c(:), dev_sig_u(:)
     
     logical, allocatable :: c_use(:,:), u_use(:,:)
     integer, allocatable :: n_mode(:)
     character(1), allocatable :: disper_phase(:)
     logical :: verb = .false.
     

   contains
     procedure :: get_nf => observation_disper_get_nf
     procedure :: get_fmin => observation_disper_get_fmin
     procedure :: get_fmax => observation_disper_get_fmax
     procedure :: get_df => observation_disper_get_df
     procedure :: get_cmin => observation_disper_get_cmin
     procedure :: get_cmax => observation_disper_get_cmax
     procedure :: get_dc => observation_disper_get_dc
     procedure :: get_c  => observation_disper_get_c
     procedure :: get_c_array => observation_disper_get_c_array
     procedure :: get_u => observation_disper_get_u
     procedure :: get_u_array => observation_disper_get_u_array
     procedure :: get_sig_c => observation_disper_get_sig_c
     procedure :: get_sig_u => observation_disper_get_sig_u
     procedure :: get_sig_c_min => observation_disper_get_sig_c_min
     procedure :: get_sig_u_min => observation_disper_get_sig_u_min
     procedure :: get_sig_c_max => observation_disper_get_sig_c_max
     procedure :: get_sig_u_max => observation_disper_get_sig_u_max
     procedure :: get_dev_sig_c => observation_disper_get_dev_sig_c
     procedure :: get_dev_sig_u => observation_disper_get_dev_sig_u
     procedure :: get_c_use => observation_disper_get_c_use
     procedure :: get_u_use => observation_disper_get_u_use
     
     procedure :: get_n_disp => observation_disper_get_n_disp
     procedure :: set_n_disp => observation_disper_set_n_disp
     procedure :: get_disper_phase &
          & => observation_disper_get_disper_phase
     procedure :: get_n_mode => observation_disper_get_n_mode
     procedure :: read_data => observation_disper_read_data
     
  end type observation_disper
  
  interface observation_disper
     module procedure :: init_observation_disper
  end interface observation_disper
     
contains
  
  !---------------------------------------------------------------------

  type(observation_disper) &
       & function init_observation_disper(rayobs_in, verb) &
       & result(self)
    character(len=*), intent(in) :: rayobs_in
    logical, intent(in), optional :: verb
    integer :: ierr, io, i
    character(line_max) :: line
    character(line_max), allocatable :: filename(:)
    type(line_text) :: lt
    
    if (present(verb)) then
       self%verb = verb
    end if

    if (self%verb) then
       write(*,'(3A)') &
            & "<< Reading observed dispersion curve from ", &
            & trim(rayobs_in), " >>"
    end if

    open(newunit = io, file = rayobs_in, status = "old", iostat = ierr)
    if (ierr /= 0) then
       if (self%verb) then
          write(0,*)"ERROR: cannot open ", trim(rayobs_in)
       end if
       call mpi_finalize(ierr)
       stop
    end if
    
    ! Read number of dispersion curves from input file
    do
       read(io, '(a)')line
       lt = line_text(line, ignore_space=.false.)
       line = lt%get_line()
       if (len_trim(line) /= 0) then
          read(line, *) self%n_disp
          if (self%verb) then
             write(*,'(A,I5)')"# of dispersion curves (n_disp)  =", &
                  & self%n_disp
          end if
          exit
       end if
    end do
    
    ! allocate
    allocate(self%nf(self%n_disp), self%fmin(self%n_disp))
    allocate(self%fmax(self%n_disp))
    allocate(self%df(self%n_disp))
    allocate(self%disper_phase(self%n_disp))
    allocate(self%cmin(self%n_disp), self%cmax(self%n_disp))
    allocate(self%dc(self%n_disp))
    allocate(filename(self%n_disp))
    allocate(self%n_mode(self%n_disp))
    allocate(self%sig_c(self%n_disp))
    allocate(self%sig_u(self%n_disp))
    allocate(self%sig_c_min(self%n_disp))
    allocate(self%sig_u_min(self%n_disp))
    allocate(self%sig_c_max(self%n_disp))
    allocate(self%sig_u_max(self%n_disp))
    allocate(self%dev_sig_c(self%n_disp))
    allocate(self%dev_sig_u(self%n_disp))

    do i = 1, self%n_disp
       ! get file name
       if (self%verb) write(*,*)"-----------------"
       do 
          read(io, '(a)')line
          lt = line_text(line, ignore_space=.false.)
          line = lt%get_line()
          if (len_trim(line) /= 0) then
             read(line, *) filename(i)
             if (self%verb) then
                write(*,'(2A)')"filename = ", trim(filename(i))
             end if
             exit
          end if
       end do
       ! disper_phase
       do 
          read(io, '(a)')line
          lt = line_text(line, ignore_space=.false.)
          line = lt%get_line()
          if (len_trim(line) /= 0) then
             read(line, *) self%disper_phase(i), self%n_mode(i)
             if (self%verb) then
                write(*,'(2A)')"Phase type (disper_phase) = ", &
                     & self%disper_phase(i)
                write(*,'(A,I5)')"Mode (n_mode) = ", &
                     & self%n_mode(i)
             end if
             exit
          end if
       end do
       
       ! nf, fmin, df
       do 
          read(io, '(a)')line
          lt = line_text(line, ignore_space=.false.)
          line = lt%get_line()
          if (len_trim(line) == 0) cycle
          read(line, *) self%nf(i), self%fmin(i), self%df(i)
          if (self%verb) then
             write(*,'(A,I5)')"# of measurements (nf) = ", &
                  & self%nf(i)
             write(*,'(A,F12.5)')"Minimum frequency (fmin) = ", &
                  & self%fmin(i)
             write(*,'(A,F12.5)')"Frequency interval (df) = ", &
                  & self%df(i)
          end if
          self%fmax(i) = &
               & self%fmin(i) + self%df(i) * (self%nf(i) - 1)
          exit
       end do
       
       ! cmin, cmax, dc
       do 
          read(io, '(a)')line
          lt = line_text(line, ignore_space=.false.)
          line = lt%get_line()
          if (len_trim(line) == 0) cycle
          read(line, *) self%cmin(i), self%cmax(i), self%dc(i)
          if (self%verb) then
             write(*,'(A,F12.5)')"Minimum phase velocity (cmin) = ", &
                  & self%cmin(i)
             write(*,'(A,F12.5)')"Maximum phase velocity (cmax) = ", &
                  & self%cmax(i)
             write(*,'(A,F12.5)')"Phase velocity interval (dc) = ", &
                  & self%dc(i)
          end if
          exit
       end do
       
       ! sig_c_min, sig_c_max, dev_sig_c
       do 
          read(io, '(a)')line
          lt = line_text(line, ignore_space=.false.)
          line = lt%get_line()
          if (len_trim(line) == 0) cycle
          read(line, *) self%sig_c_min(i), self%sig_c_max(i), &
               & self%dev_sig_c(i)
          if (self%verb) then
             write(*,'(A,F12.5)') &
                  & "Min. sigma of phase velocity (sig_c_min) = ", &
                  & self%sig_c_min(i)
             write(*,'(A,F12.5)') &
                  & "Max. sigma of phase velocity (sig_c_max) = ", &
                  & self%sig_c_max(i)
             write(*,'(A,F12.5)')"Stdev. sigma of phase velocity " &
                  & // "(dev_sig_c) = ", &
                  & self%dev_sig_c(i)
          end if
          
          ! ERROR 
          if (self%sig_c_min(i) + self%dev_sig_c(i) > self%sig_c_max(i)) then
             write(0,*)"ERROR: following relation must be satisfied"
             write(0,*)"       sig_c_min + dev_sig_c <= sig_c_max"
             call mpi_abort(ierr)
          end if


          exit
       end do

       ! sig_u_min, sig_u_max, dev_sig_u
       do 
          read(io, '(a)')line
          lt = line_text(line, ignore_space=.false.)
          line = lt%get_line()
          if (len_trim(line) == 0) cycle
          read(line, *) self%sig_u_min(i), self%sig_u_max(i), &
               & self%dev_sig_u(i)
          if (self%verb) then
             write(*,'(A,F12.5)') &
                  & "Min. sigma of phase velocity (sig_u_min) = ", &
                  & self%sig_u_min(i)
             write(*,'(A,F12.5)') &
                  & "Max. sigma of phase velocity (sig_u_max) = ", &
                  & self%sig_u_max(i)
             write(*,'(A,F12.5)')"Stdev. sigma of phase velocity " &
                  & // "(dev_sig_u) = ", &
                  & self%dev_sig_u(i)
          end if
          
          ! ERROR 
          if (self%sig_u_min(i) + self%dev_sig_u(i) > self%sig_u_max(i)) then
             write(0,*)"ERROR: following relation must be satisfied"
             write(0,*)"       sig_u_min + dev_sig_u <= sig_u_max"
             call mpi_abort(ierr)
          end if
          
          exit
       end do

       if (self%verb) write(*,*)"-----------------"
    end do
    close(io)
    
    ! Read data files
    allocate(self%c(maxval(self%nf(:)), self%n_disp))
    allocate(self%u(maxval(self%nf(:)), self%n_disp))
    allocate(self%c_use(maxval(self%nf(:)), self%n_disp))
    allocate(self%u_use(maxval(self%nf(:)), self%n_disp))
    
    do i = 1, self%n_disp
       call self%read_data(filename(i), i)
    end do

    return 
  end function init_observation_disper
  
  !---------------------------------------------------------------------
  
  integer function observation_disper_get_nf(self, i) result(nf)
    class(observation_disper), intent(in) :: self
    integer, intent(in) :: i
    
    nf = self%nf(i)

    return 
  end function observation_disper_get_nf

  !---------------------------------------------------------------------
  
  double precision function observation_disper_get_fmin(self, i) &
       & result(fmin)
    class(observation_disper), intent(in) :: self
    integer, intent(in) :: i
    
    fmin = self%fmin(i)

    return 
  end function observation_disper_get_fmin
  
  !---------------------------------------------------------------------

  double precision function observation_disper_get_fmax(self, i) &
       & result(fmax)
    class(observation_disper), intent(in) :: self
    integer, intent(in) :: i
    
    fmax = self%fmax(i)

    return 
  end function observation_disper_get_fmax
  
  !---------------------------------------------------------------------
  
  double precision function observation_disper_get_df(self, i) &
       & result(df)
    class(observation_disper), intent(in) :: self
    integer, intent(in) :: i
    
    df = self%df(i)

    return 
  end function observation_disper_get_df

  !---------------------------------------------------------------------

  double precision function observation_disper_get_cmin(self, i) &
       & result(cmin)
    class(observation_disper), intent(in) :: self
    integer, intent(in) :: i
    
    cmin = self%cmin(i)

    return 
  end function observation_disper_get_cmin
  
  !---------------------------------------------------------------------

  double precision function observation_disper_get_cmax(self, i) &
       & result(cmax)
    class(observation_disper), intent(in) :: self
    integer, intent(in) :: i
    
    cmax = self%cmax(i)

    return 
  end function observation_disper_get_cmax
  
  !---------------------------------------------------------------------
  
  double precision function observation_disper_get_dc(self, i) &
       & result(dc)
    class(observation_disper), intent(in) :: self
    integer, intent(in) :: i
    
    dc = self%dc(i)

    return 
  end function observation_disper_get_dc

  !---------------------------------------------------------------------

  double precision function observation_disper_get_c(self, i, j) &
       & result(c)
    class(observation_disper), intent(in) :: self
    integer, intent(in) :: i, j
    
    c = self%c(i, j)
    
    return 
  end function observation_disper_get_c
  
  !---------------------------------------------------------------------
  
  function observation_disper_get_c_array(self, i) result(c)
    class(observation_disper), intent(in) :: self
    integer, intent(in) :: i
    double precision :: c(self%nf(i))
    
    c(:) = self%c(:, i)
    
    return 
  end function observation_disper_get_c_array

  !---------------------------------------------------------------------

  double precision function observation_disper_get_u(self, i, j) &
       & result(u)
    class(observation_disper), intent(in) :: self
    integer, intent(in) :: i, j
    
    u = self%u(i, j)
    
    return 
  end function observation_disper_get_u
  
  !---------------------------------------------------------------------
  
  function observation_disper_get_u_array(self, i) result(u)
    class(observation_disper), intent(in) :: self
    integer, intent(in) :: i
    double precision :: u(self%nf(i))
    
    u(:) = self%u(:, i)
    
    return 
  end function observation_disper_get_u_array
  
  !---------------------------------------------------------------------

  double precision function observation_disper_get_sig_c(self, i) &
       & result(sig_c)
    class(observation_disper), intent(in) :: self
    integer, intent(in) :: i
    
    sig_c = self%sig_c(i)
    
    return 
  end function observation_disper_get_sig_c

  !---------------------------------------------------------------------

  double precision function observation_disper_get_sig_u(self, i) &
       & result(sig_u)
    class(observation_disper), intent(in) :: self
    integer, intent(in) :: i
    
    sig_u = self%sig_u(i)
    
    return 
  end function observation_disper_get_sig_u

  !---------------------------------------------------------------------
  
  double precision function observation_disper_get_sig_c_min(self, i) &
       & result(sig_c_min)
    class(observation_disper), intent(in) :: self
    integer, intent(in) :: i
    
    sig_c_min = self%sig_c_min(i)
    
    return 
  end function observation_disper_get_sig_c_min

  !---------------------------------------------------------------------

  double precision function observation_disper_get_sig_u_min(self, i) &
       & result(sig_u_min)
    class(observation_disper), intent(in) :: self
    integer, intent(in) :: i
    
    sig_u_min = self%sig_u_min(i)
    
    return 
  end function observation_disper_get_sig_u_min

  !---------------------------------------------------------------------
  
  double precision function observation_disper_get_sig_c_max(self, i) &
       & result(sig_c_max)
    class(observation_disper), intent(in) :: self
    integer, intent(in) :: i
    
    sig_c_max = self%sig_c_max(i)
    
    return 
  end function observation_disper_get_sig_c_max

  !---------------------------------------------------------------------

  double precision function observation_disper_get_sig_u_max(self, i) &
       & result(sig_u_max)
    class(observation_disper), intent(in) :: self
    integer, intent(in) :: i
    
    sig_u_max = self%sig_u_max(i)
    
    return 
  end function observation_disper_get_sig_u_max

  !---------------------------------------------------------------------

  double precision function observation_disper_get_dev_sig_c(self, i) &
       & result(dev_sig_c)
    class(observation_disper), intent(in) :: self
    integer, intent(in) :: i
    
    dev_sig_c = self%dev_sig_c(i)
    
    return 
  end function observation_disper_get_dev_sig_c

  !---------------------------------------------------------------------

  double precision function observation_disper_get_dev_sig_u(self, i) &
       & result(dev_sig_u)
    class(observation_disper), intent(in) :: self
    integer, intent(in) :: i
    
    dev_sig_u = self%dev_sig_u(i)
    
    return 
  end function observation_disper_get_dev_sig_u

  !---------------------------------------------------------------------

  logical function observation_disper_get_c_use(self, i, j) &
       & result(c_use)
    class(observation_disper), intent(in) :: self
    integer, intent(in) :: i, j
    
    c_use = self%c_use(i, j)
    
    return 
  end function observation_disper_get_c_use

  !---------------------------------------------------------------------

  logical function observation_disper_get_u_use(self, i, j) &
       & result(u_use)
    class(observation_disper), intent(in) :: self
    integer, intent(in) :: i, j
    
    u_use = self%u_use(i, j)
    
    return 
  end function observation_disper_get_u_use

  !---------------------------------------------------------------------

  integer function observation_disper_get_n_disp(self) result(n_disp)
    class(observation_disper), intent(in) :: self
    
    n_disp = self%n_disp
    
    return 
  end function observation_disper_get_n_disp

  !---------------------------------------------------------------------

  subroutine observation_disper_set_n_disp(self, n_disp)
    class(observation_disper), intent(inout) :: self
    integer, intent(in) :: n_disp

    self%n_disp = n_disp
    
    return 
  end subroutine observation_disper_set_n_disp
  
  !---------------------------------------------------------------------
  
  character(1) function observation_disper_get_disper_phase(self, i) &
       & result(disper_phase)
    class(observation_disper), intent(in) :: self
    integer, intent(in) :: i
    
    disper_phase = self%disper_phase(i)

    return 
  end function observation_disper_get_disper_phase

  !---------------------------------------------------------------------

  integer function observation_disper_get_n_mode(self, i) result(n_mode)
    class(observation_disper), intent(in) :: self
    integer, intent(in) :: i

    n_mode = self%n_mode(i)
    
    return 
  end function observation_disper_get_n_mode

  !---------------------------------------------------------------------
  
  subroutine observation_disper_read_data(self, filename, i)
    class(observation_disper), intent(inout) :: self
    character(*), intent(in) :: filename
    integer, intent(in) :: i
    integer :: io, ierr, j
    character(line_max) :: line
    type(line_text) :: lt

    open(newunit = io, file = filename, status = 'old', &
         & iostat = ierr)
    if (ierr /= 0) then
       write(0,*)"ERROR: cannot open ", trim(filename)
       call mpi_finalize(ierr)
       stop
    end if
    do j = 1, self%nf(i)
       do
          read(io, '(a)')line
          lt = line_text(line, ignore_space = .false.)
          line = lt%get_line()
          if (len_trim(line) == 0) cycle
          read(line,*,iostat = ierr)self%c(j, i), self%c_use(j, i), &
               & self%u(j, i), self%u_use(j, i)
          if (ierr /= 0) then
             write(0,*)"ERROR: while reading ", trim(filename), & 
                  & " of Line", j, "th data"
             call mpi_finalize(ierr)
             stop
          end if
          exit
       end do
    end do
    close(io)
    
    return 
  end subroutine observation_disper_read_data
  
  !---------------------------------------------------------------------
end module cls_observation_disper
