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
     double precision, allocatable :: sig_c(:,:) ! uncertainties in c
     double precision, allocatable :: sig_u(:,:) ! uncertainties in u
     double precision, allocatable :: cmin(:), cmax(:), dc(:)
     character(1), allocatable :: phase(:)
     
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
     procedure :: get_sig_c_array => observation_disper_get_sig_c_array
     procedure :: get_sig_u => observation_disper_get_sig_u
     procedure :: get_sig_u_array => observation_disper_get_sig_u_array
     procedure :: get_n_disp => observation_disper_get_n_disp
     procedure :: read_data => observation_disper_read_data
     
  end type observation_disper
  
  interface observation_disper
     module procedure :: init_observation_disper
  end interface observation_disper
     
contains
  
  !---------------------------------------------------------------------

  type(observation_disper) function init_observation_disper(rayobs_in) &
       & result(self)
    character(len=*), intent(in) :: rayobs_in
    integer :: ierr, io, i
    character(line_max) :: line
    character(line_max), allocatable :: filename(:)
    type(line_text) :: lt
    
    write(*,*)"Reading observed Rayleigh wave dispersion curve from ", &
         & trim(rayobs_in)

    open(newunit = io, file = rayobs_in, status = "old", iostat = ierr)
    if (ierr /= 0) then
       write(0,*)"ERROR: cannot open ", trim(rayobs_in)
       stop
    end if
    
    ! Read number of dispersion curves from input file
    do
       read(io, '(a)')line
       lt = line_text(line, ignore_space=.false.)
       line = lt%get_line()
       if (len_trim(line) /= 0) then
          read(line, *) self%n_disp
          write(*,*)"n_disp =", self%n_disp
          exit
       end if
    end do

    ! allocate
    allocate(self%nf(self%n_disp), self%fmin(self%n_disp))
    allocate(self%fmax(self%n_disp))
    allocate(self%df(self%n_disp))
    allocate(self%phase(self%n_disp))
    allocate(self%cmin(self%n_disp), self%cmax(self%n_disp))
    allocate(self%dc(self%n_disp))
    allocate(filename(self%n_disp))

    do i = 1, self%n_disp
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
       ! phase
       do 
          read(io, '(a)')line
          lt = line_text(line, ignore_space=.false.)
          line = lt%get_line()
          if (len_trim(line) /= 0) then
             read(line, *) self%phase(i)
             exit
          end if
       end do
       
       ! nf, fmin, df
       do 
          read(io, '(a)')line
          lt = line_text(line, ignore_space=.false.)
          line = lt%get_line()
          if (len_trim(line) /= 0) then
             read(line, *) self%nf(i), self%fmin(i), self%df(i)
             self%fmax(i) = &
                  & self%fmin(i) + self%df(i) * (self%nf(i) - 1)
             exit
          end if
       end do
       
       ! cmin, cmax, dc
       do 
          read(io, '(a)')line
          lt = line_text(line, ignore_space=.false.)
          line = lt%get_line()
          if (len_trim(line) /= 0) then
             read(line, *) self%cmin(i), self%cmax(i), self%dc(i)
             exit
          end if
       end do

       
    end do
    close(io)
    
    ! allocate data array
    allocate(self%c(maxval(self%nf(:)), self%n_disp))
    allocate(self%u(maxval(self%nf(:)), self%n_disp))
    allocate(self%sig_c(maxval(self%nf(:)), self%n_disp))
    allocate(self%sig_u(maxval(self%nf(:)), self%n_disp))
    
    ! Read data files
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

  double precision function observation_disper_get_sig_c(self, i, j) &
       & result(sig_c)
    class(observation_disper), intent(in) :: self
    integer, intent(in) :: i, j
    
    sig_c = self%sig_c(i, j)
    
    return 
  end function observation_disper_get_sig_c

  !---------------------------------------------------------------------
  
  function observation_disper_get_sig_c_array(self, i) result(sig_c)
    class(observation_disper), intent(in) :: self
    integer, intent(in) :: i
    double precision :: sig_c(self%nf(i))

    sig_c(:) = self%sig_c(:, i)
    
    return 
  end function observation_disper_get_sig_c_array

  !---------------------------------------------------------------------

  double precision function observation_disper_get_sig_u(self, i, j) &
       & result(sig_u)
    class(observation_disper), intent(in) :: self
    integer, intent(in) :: i, j
    
    sig_u = self%sig_u(i, j)
    
    return 
  end function observation_disper_get_sig_u

  !---------------------------------------------------------------------

  function observation_disper_get_sig_u_array(self, i) result(sig_u)
    class(observation_disper), intent(in) :: self
    integer, intent(in) :: i
    double precision :: sig_u(self%nf(i))

    sig_u(:) = self%sig_u(:, i)
    
    return 
  end function observation_disper_get_sig_u_array

  !---------------------------------------------------------------------

  integer function observation_disper_get_n_disp(self) result(n_disp)
    class(observation_disper), intent(in) :: self
    
    n_disp = self%n_disp
    
    return 
  end function observation_disper_get_n_disp

  !---------------------------------------------------------------------

  subroutine observation_disper_read_data(self, filename, i)
    class(observation_disper), intent(inout) :: self
    character(*), intent(in) :: filename
    integer, intent(in) :: i
    integer :: io, ierr, j

    open(newunit = io, file = filename, status = 'old', &
         & iostat = ierr)
    if (ierr /= 0) then
       write(0,*)"ERROR: cannot open ", trim(filename)
       call mpi_finalize(ierr)
       stop
    end if
    do j = 1, self%nf(i)
       read(io,*)self%c(j, i), self%sig_c(j, i), &
            & self%u(j, i), self%sig_u(j, i)
    end do
    close(io)
    
    return 
  end subroutine observation_disper_read_data
  
  !---------------------------------------------------------------------
end module cls_observation_disper
