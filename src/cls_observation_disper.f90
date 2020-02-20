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
  
  implicit none 
  
  type observation_disper
     private
     ! surface wave
     integer :: nf
     double precision :: fmin, df
     double precision, allocatable :: c(:) ! phase velocity
     double precision, allocatable :: u(:) ! group velocity
     double precision, allocatable :: sig_c(:) ! uncertainties in c
     double precision, allocatable :: sig_u(:) ! uncertainties in u
     
   contains
     procedure :: get_nf             => observation_get_nf
     procedure :: get_fmin           => observation_get_fmin
     procedure :: get_df             => observation_get_df
     procedure :: get_c              => observation_get_c
     procedure :: get_c_array        => observation_get_c_array
     procedure :: get_u              => observation_get_u
     procedure :: get_u_array        => observation_get_u_array
     procedure :: get_sig_c          => observation_get_sig_c
     procedure :: get_sig_c_array    => observation_get_sig_c_array
     procedure :: get_sig_u          => observation_get_sig_u
     procedure :: get_sig_u_array    => observation_get_sig_u_array
     
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

    write(*,*)"Reading observed Rayleigh wave dispersion curve from ", &
         & trim(rayobs_in)

    open(newunit = io, file = rayobs_in, status = "old", iostat = ierr)
    if (ierr /= 0) then
       write(0,*)"ERROR: cannot open ", trim(rayobs_in)
       write(0,*)"     : (init_observation_disper)"
       stop
    end if
    
    ! Read file
    read(io, *)self%nf, self%fmin, &
         & self%df
    allocate(self%c(self%nf))
    allocate(self%u(self%nf))
    allocate(self%sig_c(self%nf))
    allocate(self%sig_u(self%nf))
    do i = 1, self%nf
       read(io, *)self%c(i), self%sig_c(i), &
            & self%u(i), self%sig_u(i)
    end do

    close(io)
    write(*,*)
    
    return 
  end function init_observation_disper
  
  !---------------------------------------------------------------------
  
  integer function observation_get_nf(self) result(nf)
    class(observation_disper), intent(in) :: self
    
    nf = self%nf

    return 
  end function observation_get_nf

  !---------------------------------------------------------------------
  
  double precision function observation_get_fmin(self) result(fmin)
    class(observation_disper), intent(in) :: self
    
    fmin = self%fmin

    return 
  end function observation_get_fmin
  
  !---------------------------------------------------------------------
  
  double precision function observation_get_df(self) result(df)
    class(observation_disper), intent(in) :: self
    
    df = self%df

    return 
  end function observation_get_df

  !---------------------------------------------------------------------

  double precision function observation_get_c(self, i) result(c)
    class(observation_disper), intent(in) :: self
    integer, intent(in) :: i
    
    c = self%c(i)
    
    return 
  end function observation_get_c
  
  !---------------------------------------------------------------------
  
  function observation_get_c_array(self) result(c)
    class(observation_disper), intent(in) :: self
    double precision :: c(self%nf)
    
    c(:) = self%c(:)
    
    return 
  end function observation_get_c_array

  !---------------------------------------------------------------------

  double precision function observation_get_u(self, i) result(u)
    class(observation_disper), intent(in) :: self
    integer, intent(in) :: i
    
    u = self%u(i)
    
    return 
  end function observation_get_u

  !---------------------------------------------------------------------
  
  function observation_get_u_array(self) result(u)
    class(observation_disper), intent(in) :: self
    double precision :: u(self%nf)

    u(:) = self%u(:)
    
    return 
  end function observation_get_u_array
  
  !---------------------------------------------------------------------

  double precision function observation_get_sig_c(self, i) result(sig_c)
    class(observation_disper), intent(in) :: self
    integer, intent(in) :: i
    
    sig_c = self%sig_c(i)
    
    return 
  end function observation_get_sig_c

  !---------------------------------------------------------------------
  
  function observation_get_sig_c_array(self) result(sig_c)
    class(observation_disper), intent(in) :: self
    double precision :: sig_c(self%nf)

    sig_c(:) = self%sig_c(:)
    
    return 
  end function observation_get_sig_c_array

  !---------------------------------------------------------------------

  double precision function observation_get_sig_u(self, i) result(sig_u)
    class(observation_disper), intent(in) :: self
    integer, intent(in) :: i
    
    sig_u = self%sig_u(i)
    
    return 
  end function observation_get_sig_u

  !---------------------------------------------------------------------

  function observation_get_sig_u_array(self) result(sig_u)
    class(observation_disper), intent(in) :: self
    double precision :: sig_u(self%nf)

    sig_u(:) = self%sig_u(:)
    
    return 
  end function observation_get_sig_u_array

  !---------------------------------------------------------------------

  
end module cls_observation_disper
