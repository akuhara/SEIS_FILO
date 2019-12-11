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
     integer :: n_iter = -999
     integer :: n_burn = -999
     integer :: n_corr = -999
     integer :: n_chain = -999
     integer :: n_cool = -999
     integer :: i_seed1 = 11111111
     integer :: i_seed2 = 22222222
     integer :: i_seed3 = 33333333
     integer :: i_seed4 = 44444444
     integer :: k_min = -999
     integer :: k_max = -999
     integer :: nbin_z = -999
     integer :: nbin_vs = -999
     integer :: nbin_vp = -999
     integer :: nbin_c  = -999

     double precision :: temp_high = -999.d0
     double precision :: fmin = -999.d0
     double precision :: fmax = -999.d0
     double precision :: df = -999.d0
     double precision :: cmin = -999.d0
     double precision :: cmax = -999.d0
     double precision :: dc = -999.d0
     double precision :: dev_vs = -999.d0
     double precision :: dev_vp = -999.d0
     double precision :: dev_z  = -999.d0
     double precision :: vs_min = -999.d0
     double precision :: vs_max = -999.d0
     double precision :: vp_min = -999.d0
     double precision :: vp_max = -999.d0
     double precision :: z_min = -999.d0
     double precision :: z_max = -999.d0
     double precision :: ocean_thick = -999.d0

     character(len=line_max) :: vmod_in = ""
     character(len=line_max) :: ray_out = ""
     character(len=line_max) :: obs_in = ""

     logical :: solve_vp = .false.
     logical :: ocean_flag = .false.
     logical :: verb = .false.


   contains
     procedure :: read_file => param_read_file
     procedure :: read_line => param_read_line
     procedure :: read_value => param_read_value
     procedure :: get_n_iter => param_get_n_iter
     procedure :: get_n_burn => param_get_n_burn
     procedure :: get_n_corr => param_get_n_corr
     procedure :: get_n_chain => param_get_n_chain
     procedure :: get_n_cool =>  param_get_n_cool
     procedure :: get_i_seed1 => param_get_i_seed1
     procedure :: get_i_seed2 => param_get_i_seed2
     procedure :: get_i_seed3 => param_get_i_seed3
     procedure :: get_i_seed4 => param_get_i_seed4
     procedure :: get_k_min => param_get_k_min
     procedure :: get_k_max => param_get_k_max
     procedure :: get_nbin_z => param_get_nbin_z
     procedure :: get_nbin_vs => param_get_nbin_vs
     procedure :: get_nbin_vp => param_get_nbin_vp
     procedure :: get_nbin_c => param_get_nbin_c
     
     procedure :: get_temp_high => param_get_temp_high
     procedure :: get_fmin => param_get_fmin
     procedure :: get_fmax => param_get_fmax
     procedure :: get_df => param_get_df
     procedure :: get_cmin => param_get_cmin
     procedure :: get_cmax => param_get_cmax
     procedure :: get_dc => param_get_dc
     procedure :: get_dev_vs => param_get_dev_vs
     procedure :: get_dev_vp => param_get_dev_vp
     procedure :: get_dev_z => param_get_dev_z
     procedure :: get_vs_min => param_get_vs_min
     procedure :: get_vs_max => param_get_vs_max
     procedure :: get_vp_min => param_get_vp_min
     procedure :: get_vp_max => param_get_vp_max
     procedure :: get_z_min => param_get_z_min
     procedure :: get_z_max => param_get_z_max
     procedure :: get_ocean_thick => param_get_ocean_thick

     procedure :: get_solve_vp => param_get_solve_vp
     procedure :: get_ocean_flag => param_get_ocean_flag
     
     procedure :: get_vmod_in => param_get_vmod_in
     procedure :: get_ray_out => param_get_ray_out
     procedure :: get_obs_in  => param_get_obs_in

  end type param
  
  interface param
     module procedure init_param
  end interface param



contains
  
  !---------------------------------------------------------------------
  
  type(param) function init_param(param_file, verb) 
    character(len=*), intent(in) :: param_file
    logical, intent(in), optional :: verb

    if (present(verb)) then
       init_param%verb = verb
    end if
    
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
    integer :: nlen, j, itmp
    double precision :: rtmp
    logical :: ltmp

    
    
    nlen = len(str)
    j = index(str, "=")
    if (j == 0 .or. j == 1 .or. j == nlen) then
       write(0,*)"ERROR: invalid parameter in ", trim(self%param_file)
       write(0,*)"     : ", str, "   (?)"
       stop
    end if
    
    name = str(1:j-1)
    var = str(j+1:nlen)
    if (self%verb) then
       write(*,*)trim(name), " <- ", trim(var)
    end if
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
    else if (name == "n_iter") then
       read(var, *) itmp
       self%n_iter = itmp
    else if (name == "n_burn") then
       read(var, *) itmp
       self%n_burn = itmp
    else if (name == "n_corr") then
       read(var, *) itmp
       self%n_corr = itmp
    else if (name == "n_chain") then
       read(var, *) itmp
       self%n_chain = itmp
    else if (name == "n_cool") then
       read(var, *) itmp
       self%n_cool = itmp
    else if (name == "temp_high") then
       read(var, *) rtmp
       self%temp_high = rtmp
    else if (name == "i_seed1") then
       read(var, *) itmp
       self%i_seed1 = itmp 
    else if (name == "i_seed2") then
       read(var, *) itmp
       self%i_seed2 = itmp
    else if (name == "i_seed3") then
       read(var, *) itmp
       self%i_seed3 = itmp
    else if (name == "i_seed4") then
       read(var, *) itmp
       self%i_seed4 = itmp
    else if (name == "k_min") then
       read(var, *) itmp
       self%k_min = itmp
    else if (name == "k_max") then
       read(var, *) itmp
       self%k_max = itmp
    else if (name == "obs_in") then
       self%obs_in = var
    else if (name == "dev_vs") then
       read(var, *) rtmp
       self%dev_vs = rtmp
    else if (name == "dev_vp") then
       read(var, *) rtmp
       self%dev_vp = rtmp
    else if (name == "dev_z") then
       read(var, *) rtmp
       self%dev_z = rtmp
    else if (name == "vs_min") then
       read(var, *) rtmp
       self%vs_min = rtmp
    else if (name == "vs_max") then
       read(var, *) rtmp
       self%vs_max = rtmp
    else if (name == "vp_min") then
       read(var, *) rtmp
       self%vp_min = rtmp
    else if (name == "vp_max") then
       read(var, *) rtmp
       self%vp_max = rtmp
    else if (name == "z_min") then
       read(var, *) rtmp
       self%z_min = rtmp
    else if (name == "z_max") then
       read(var, *) rtmp
       self%z_max = rtmp
    else if (name == "ocean_thick") then
       read(var, *) rtmp
       self%ocean_thick = rtmp
    else if (name == "solve_vp") then
       read(var, *) ltmp
       self%solve_vp = ltmp
    else if (name == "ocean_flag") then
       read(var, *) ltmp
       self%ocean_flag = ltmp
    else if (name == "nbin_z") then
       read(var, *) itmp
       self%nbin_z = itmp
    else if (name == "nbin_vs") then
       read(var, *) itmp
       self%nbin_vs = itmp
    else if (name == "nbin_vp") then
       read(var, *) itmp
       self%nbin_vp = itmp
    else if (name == "nbin_c") then
       read(var, *) itmp
       self%nbin_c = itmp
    else
       write(0,*)"Warnings: Invalid parameter name"
       write(0,*)"        : ", name, "  (?)"
    end if
    return 
  end subroutine param_read_value

  !---------------------------------------------------------------------

  integer function param_get_n_iter(self) result(n_iter)
    class(param), intent(in) :: self
    
    n_iter = self%n_iter
    
    return 
  end function param_get_n_iter

  !---------------------------------------------------------------------
  
  integer function param_get_n_burn(self) result(n_burn)
    class(param), intent(in) :: self
    
    n_burn = self%n_burn
    
    return 
  end function param_get_n_burn

  !---------------------------------------------------------------------

  integer function param_get_n_corr(self) result(n_corr)
    class(param), intent(in) :: self
    
    n_corr = self%n_corr
    
    return 
  end function param_get_n_corr

  
  !---------------------------------------------------------------------

  integer function param_get_n_chain(self) result(n_chain)
    class(param), intent(in) :: self
    
    n_chain = self%n_chain
    
    return 
  end function param_get_n_chain
  
  !---------------------------------------------------------------------

  integer function param_get_n_cool(self) result(n_cool)
    class(param), intent(in) :: self
    
    n_cool = self%n_cool
    
    return 
  end function param_get_n_cool
  
  !---------------------------------------------------------------------

  integer function param_get_i_seed1(self) result(i_seed1)
    class(param), intent(in) :: self
    
    i_seed1 = self%i_seed1
    
    return 
  end function param_get_i_seed1
  
  !---------------------------------------------------------------------

  integer function param_get_i_seed2(self) result(i_seed2)
    class(param), intent(in) :: self
    
    i_seed2 = self%i_seed2
    
    return 
  end function param_get_i_seed2
  
  !---------------------------------------------------------------------

  integer function param_get_i_seed3(self) result(i_seed3)
    class(param), intent(in) :: self
    
    i_seed3 = self%i_seed3
    
    return 
  end function param_get_i_seed3
  
  !---------------------------------------------------------------------

  integer function param_get_i_seed4(self) result(i_seed4)
    class(param), intent(in) :: self
    
    i_seed4 = self%i_seed4
    
    return 
  end function param_get_i_seed4
  
  !---------------------------------------------------------------------

  integer function param_get_k_min(self) result(k_min)
    class(param), intent(in) :: self
    
    k_min = self%k_min
    
    return 
  end function param_get_k_min
  
  !---------------------------------------------------------------------

  integer function param_get_k_max(self) result(k_max)
    class(param), intent(in) :: self
    
    k_max = self%k_max
    
    return 
  end function param_get_k_max
  
  !---------------------------------------------------------------------

  integer function param_get_nbin_z(self) result(nbin_z)
    class(param), intent(in) :: self
    
    nbin_z = self%nbin_z
    
    return 
  end function param_get_nbin_z
  
  
  !---------------------------------------------------------------------

  integer function param_get_nbin_vs(self) result(nbin_vs)
    class(param), intent(in) :: self
    
    nbin_vs = self%nbin_vs
    
    return 
  end function param_get_nbin_vs
  
  !---------------------------------------------------------------------

  integer function param_get_nbin_vp(self) result(nbin_vp)
    class(param), intent(in) :: self
    
    nbin_vp = self%nbin_vp
    
    return 
  end function param_get_nbin_vp

  !---------------------------------------------------------------------

  integer function param_get_nbin_c(self) result(nbin_c)
    class(param), intent(in) :: self
    
    nbin_c = self%nbin_c
    
    return 
  end function param_get_nbin_c

  !---------------------------------------------------------------------

  
  double precision function param_get_temp_high(self) result(temp_high)
    class(param), intent(in) :: self

    temp_high = self%temp_high

    return
  end function param_get_temp_high
  
  !---------------------------------------------------------------------
  
  double precision function param_get_fmin(self) result(fmin)
    class(param), intent(in) :: self

    fmin = self%fmin

    return
  end function param_get_fmin
  
  !---------------------------------------------------------------------

  double precision function param_get_fmax(self) result(fmax)
    class(param), intent(in) :: self

    fmax = self%fmax

    return
  end function param_get_fmax
  !---------------------------------------------------------------------

  double precision function param_get_df(self) result(df)
    class(param), intent(in) :: self

    df = self%df

    return
  end function param_get_df
  
  !---------------------------------------------------------------------

  double precision function param_get_cmin(self) result(cmin)
    class(param), intent(in) :: self

    cmin = self%cmin

    return
  end function param_get_cmin
  
  !---------------------------------------------------------------------

  double precision function param_get_cmax(self) result(cmax)
    class(param), intent(in) :: self

    cmax = self%cmax

    return
  end function param_get_cmax
  !---------------------------------------------------------------------

  double precision function param_get_dc(self) result(dc)
    class(param), intent(in) :: self

    dc = self%dc

    return
  end function param_get_dc

  !---------------------------------------------------------------------

  double precision function param_get_dev_vs(self) result(dev_vs)
    class(param), intent(in) :: self

    dev_vs = self%dev_vs

    return
  end function param_get_dev_vs

  !---------------------------------------------------------------------

  double precision function param_get_dev_vp(self) result(dev_vp)
    class(param), intent(in) :: self

    dev_vp = self%dev_vp

    return
  end function param_get_dev_vp

  !---------------------------------------------------------------------

  double precision function param_get_dev_z(self) result(dev_z)
    class(param), intent(in) :: self

    dev_z = self%dev_z

    return
  end function param_get_dev_z

  !---------------------------------------------------------------------
  
  double precision function param_get_vs_min(self) result(vs_min)
    class(param), intent(in) :: self

    vs_min = self%vs_min

    return
  end function param_get_vs_min

  !---------------------------------------------------------------------

  double precision function param_get_vs_max(self) result(vs_max)
    class(param), intent(in) :: self

    vs_max = self%vs_max

    return
  end function param_get_vs_max

  !---------------------------------------------------------------------
  
  double precision function param_get_vp_min(self) result(vp_min)
    class(param), intent(in) :: self

    vp_min = self%vp_min

    return
  end function param_get_vp_min

  !---------------------------------------------------------------------

  double precision function param_get_vp_max(self) result(vp_max)
    class(param), intent(in) :: self

    vp_max = self%vp_max

    return
  end function param_get_vp_max

  !---------------------------------------------------------------------
  
  double precision function param_get_z_min(self) result(z_min)
    class(param), intent(in) :: self

    z_min = self%z_min

    return
  end function param_get_z_min

  !---------------------------------------------------------------------

  double precision function param_get_z_max(self) result(z_max)
    class(param), intent(in) :: self

    z_max = self%z_max

    return
  end function param_get_z_max

  !---------------------------------------------------------------------

  double precision function param_get_ocean_thick(self) &
       & result(ocean_thick)
    class(param), intent(in) :: self

    ocean_thick = self%ocean_thick

    return
  end function param_get_ocean_thick

  !---------------------------------------------------------------------

  logical function param_get_solve_vp(self) result(solve_vp)
    class(param), intent(in) :: self

    solve_vp = self%solve_vp
    
    return 
  end function param_get_solve_vp
  
  !---------------------------------------------------------------------

  logical function param_get_ocean_flag(self) result(ocean_flag)
    class(param), intent(in) :: self

    ocean_flag = self%ocean_flag
    
    return 
  end function param_get_ocean_flag
  
  !---------------------------------------------------------------------
  
  character(len=line_max) function param_get_vmod_in(self) &
       & result(vmod_in)
    class(param), intent(in) :: self
    
    vmod_in = self%vmod_in
    
    return
  end function param_get_vmod_in
  
  !---------------------------------------------------------------------
  
  character(len=line_max) function param_get_ray_out(self) &
       & result(ray_out)
    class(param), intent(in) :: self
    
    ray_out = self%ray_out
    
    return
  end function param_get_ray_out
  
  !---------------------------------------------------------------------

  character(len=line_max) function param_get_obs_in(self) &
       & result(obs_in)
    class(param), intent(in) :: self
    
    obs_in = self%obs_in
    
    return
  end function param_get_obs_in
  
  !---------------------------------------------------------------------


end module mod_param
