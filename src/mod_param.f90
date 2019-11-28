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
module mod_param
  implicit none 
  
  integer, parameter, private :: line_max = 200
  
  type param
     private
     character(len=line_max) :: param_file
     double precision :: fmin = -999.d0
     double precision :: fmax = -999.d0
     double precision :: df = -999.d0
     double precision :: cmin = -999.d0
     double precision :: cmax = -999.d0
     double precision :: dc = -999.d0
     character(len=line_max) :: vmod_in = ""
     character(len=line_max) :: ray_out = ""
   contains
     procedure :: read_file => param_read_file
     procedure :: read_line => param_read_line
     procedure :: read_value => param_read_value
     procedure :: get_fmin => param_get_fmin
     procedure :: get_fmax => param_get_fmax
     procedure :: get_df => param_get_df
     procedure :: get_cmin => param_get_cmin
     procedure :: get_cmax => param_get_cmax
     procedure :: get_dc => param_get_dc
     procedure :: get_vmod_in => param_get_vmod_in
     procedure :: get_ray_out => param_get_ray_out
  end type param
  
  interface param
     module procedure init_param
  end interface param



contains
  
  !---------------------------------------------------------------------
  
  type(param) function init_param(param_file) 
    character(len=*), intent(in) :: param_file
 
    init_param%param_file = param_file
    init_param%vmod_in = ""
    init_param%fmin = -999.d0
    init_param%fmax = -999.d0
    init_param%df = -999.d0
    init_param%cmin = -999.d0
    init_param%cmax = -999.d0
    init_param%dc = -999.d0
    call init_param%read_file()
    
    return 
  end function init_param

  !---------------------------------------------------------------------
    
  subroutine param_read_file(self)
    class(param), intent(inout) :: self
    character(len=line_max) :: line
    integer :: ierr, io

    write(*,*)"Reading parameters from ", trim(self%param_file)

    open(newunit = io, file = self%param_file, &
         & status = 'old', iostat = ierr)
    if (ierr /= 0) then
       write(0,*) "ERROR: cannot open ", trim(self%param_file)
       write(0,*) "     : (param_read_file)"
       stop
    end if
    
    do 
       read(io, '(a)', iostat=ierr) line
       if (ierr /= 0) then
          exit
       end if
       call self%read_line(line)
    end do
    close(io)
    
    write(*,*)

    return 
  end subroutine param_read_file

  !---------------------------------------------------------------------

  subroutine param_read_line(self, line)
    class(param), intent(inout) :: self
    character(len=line_max), intent(in) :: line
    integer :: i, nlen, j

    nlen = len_trim(line)
    i = 1
    j = index(line(i:nlen), "=")
    if (j == 0) then
       return
    end if
    do while (i <= nlen)
       if (line(i:i) == " ") then
          i = i + 1
          cycle
       else if (line(i:i) == "#") then
          return
       end if
       j = index(line(i:nlen), " ")
       if (j /= 0) then
          call self%read_value(line(i:i+j-2))
          i = i + j - 1
       else
          call self%read_value(line(i:nlen))
          return
       end if
    end do

    return 
  end subroutine param_read_line

  !---------------------------------------------------------------------

  subroutine param_read_value(self, str)
    class(param), intent(inout) :: self
    character(len=*), intent(in) :: str
    character(len=line_max) :: name, var
    integer :: nlen, j
    double precision :: rtmp
    
    
    nlen = len(str)
    j = index(str, "=")
    if (j == 0 .or. j == 1 .or. j == nlen) then
       write(0,*)"ERROR: invalid parameter in ", trim(self%param_file)
       write(0,*)"     : ", str, "   (?)"
       stop
    end if
    
    name = str(1:j-1)
    var = str(j+1:nlen)
    write(*,*)trim(name), " <- ", trim(var)
    if (name == "fmin") then
       read(var, *) rtmp
       self%fmin = rtmp
    else if (name == "fmax") then
       read(var, *) rtmp
       self%fmax = rtmp
    else if (name == "df") then
       read(var, *) rtmp
       self%df = rtmp
    else if (name == "cmin") then
       read(var, *) rtmp
       self%cmin = rtmp
    else if (name == "cmax") then
       read(var, *) rtmp
       self%cmax = rtmp
    else if (name == "dc") then
       read(var, *) rtmp
       self%dc = rtmp
    else if (name == "vmod_in") then
       self%vmod_in = var
    else if (name == "ray_out") then
       self%ray_out = var
    else
       write(0,*)"Warnings: Invalid parameter name"
       write(0,*)"        : ", name, "  (?)"
    end if
    return 
  end subroutine param_read_value

  !---------------------------------------------------------------------

  double precision function param_get_fmin(self) result(fmin)
    class(param), intent(inout) :: self

    fmin = self%fmin

    return
  end function param_get_fmin
  
  !---------------------------------------------------------------------

  double precision function param_get_fmax(self) result(fmax)
    class(param), intent(inout) :: self

    fmax = self%fmax

    return
  end function param_get_fmax
  !---------------------------------------------------------------------

  double precision function param_get_df(self) result(df)
    class(param), intent(inout) :: self

    df = self%df

    return
  end function param_get_df
  
  !---------------------------------------------------------------------

  double precision function param_get_cmin(self) result(cmin)
    class(param), intent(inout) :: self

    cmin = self%cmin

    return
  end function param_get_cmin
  
  !---------------------------------------------------------------------

  double precision function param_get_cmax(self) result(cmax)
    class(param), intent(inout) :: self

    cmax = self%cmax

    return
  end function param_get_cmax
  !---------------------------------------------------------------------

  double precision function param_get_dc(self) result(dc)
    class(param), intent(inout) :: self

    dc = self%dc

    return
  end function param_get_dc

  !---------------------------------------------------------------------

  character(len=line_max) function param_get_vmod_in(self) &
       & result(vmod_in)
    class(param), intent(inout) :: self
    
    vmod_in = self%vmod_in
    
    return
  end function param_get_vmod_in
  
  !---------------------------------------------------------------------
  character(len=line_max) function param_get_ray_out(self) &
       & result(ray_out)
    class(param), intent(inout) :: self
    
    ray_out = self%ray_out
    
    return
  end function param_get_ray_out
  !---------------------------------------------------------------------


end module mod_param
