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
module cls_disper
  use cls_vmodel
  implicit none 
  
  double precision, private, parameter :: pi = acos(-1.d0)
  double precision, private, parameter :: eps = 1.0d-10
  double precision, private, parameter :: huge = 100.d0
  integer, private, parameter :: i12 = 1, i13 = 2, i14 = 3
  integer, private, parameter :: i23 = 4, i24 = 5


  type disper
     private
     type(vmodel) :: vmodel
     
     double precision :: fmin 
     double precision :: fmax
     double precision :: df  
     integer :: nf

     double precision :: cmin 
     double precision :: cmax
     double precision :: dc  
     integer :: nc
     integer :: niter = 8
     integer :: io
     integer, allocatable :: n_fc(:,:)
     integer, allocatable :: n_fu(:,:)
     character(1) :: disper_phase
     integer :: n_mode

     logical :: full_calculation = .false.

     double precision, allocatable :: c(:) ! phase velocity
     double precision, allocatable :: u(:) ! group velocity
     double precision :: y(5)
     logical :: is_ocean
     
     character(len=200) :: disper_out = ""
     logical :: out_flag = .false.

   contains
     procedure :: dispersion => disper_dispersion
     procedure :: init_propagation => disper_init_propagation
     procedure :: do_propagation => disper_do_propagation
     procedure :: solid_propagator => disper_solid_propagator
     procedure :: liquid_propagator => disper_liquid_propagator
     procedure :: find_root => disper_find_root
     procedure :: set_full_calculation => disper_set_full_calculation
     procedure :: do_full_calculation => disper_do_full_calculation
     procedure :: set_vmodel => disper_set_vmodel
     procedure :: get_nf => disper_get_nf
     procedure :: get_nc => disper_get_nc
     procedure :: get_c => disper_get_c
     procedure :: get_c_array => disper_get_c_array
     procedure :: get_u => disper_get_u
     procedure :: get_u_array => disper_get_u_array
     procedure :: save_syn => disper_save_syn
     procedure :: get_n_fc => disper_get_n_fc
     procedure :: get_n_fu => disper_get_n_fu
     
  end type disper

  interface disper
     module procedure init_disper
  end interface disper
  
contains
  
  !---------------------------------------------------------------------
  
  type(disper) function init_disper(vm, fmin, fmax, df, &
       & cmin, cmax, dc, disper_phase, n_mode, disper_out) result(self)
    type(vmodel), intent(in) :: vm
    double precision, intent(in) :: fmin, fmax, df
    double precision, intent(in) :: cmin, cmax, dc
    character(*), intent(in) :: disper_phase
    integer, intent(in) :: n_mode
    
    character(*), intent(in), optional :: disper_out
    integer :: nlay, ierr
    
    ! velocity model
    self%vmodel = vm
    nlay = self%vmodel%get_nlay()
    if (nlay < 1) then
       write(0,*)"ERROR: velocity model is not defined"
       stop
    end if
    if (self%vmodel%get_vs(1) < 0.d0) then
       self%is_ocean = .true.
    else
       self%is_ocean = .false.
    end if
    
    ! Check parameters
    self%fmin = fmin
    self%fmax = fmax
    self%df = df
    self%nf = nint((fmax - fmin) / df) + 1
    if (self%nf < 1) then
       write(0,*)"ERROR: fmax must be .gt. fmin (self)"
       stop
    end if
    allocate(self%c(self%nf), &
         &   self%u(self%nf))
    
    self%cmin = cmin
    self%cmax = cmax
    self%dc = dc
    self%nc = nint((cmax - cmin) / dc) + 1
    if (self%nc < 1) then
       write(0,*)"ERROR: cmax must be .gt. cmin (self)"
       stop
    end if
    self%n_mode = n_mode
    
    self%disper_phase = disper_phase
    if (self%disper_phase == "L") then
       write(0,*)"ERROR: Love wave it not supported yet"
       stop
    end if
    if (self%disper_phase /= "L" .and. self%disper_phase /= "R") then
       write(0,*)"ERROR: disper_phase must be L or R"
       stop
    end if

    ! output file
    if (present(disper_out)) then
       self%disper_out = disper_out
       open(newunit = self%io, status = "unknown", &
            & file = disper_out, iostat=ierr)
       if (ierr /= 0) then
          write(0, *)"ERROR: cannot creat ", trim(disper_out)
          write(0, *)"     : (init_disper)"
          stop
       end if
       self%out_flag = .true.
    end if
    
    ! for MCMC
    allocate(self%n_fc(self%nf, self%nc), self%n_fu(self%nf, self%nc))
    self%n_fc(:,:) = 0
    self%n_fu(:,:) = 0
   
    return 
  end function init_disper
    
  !---------------------------------------------------------------------

  subroutine disper_dispersion(self)
    class(disper), intent(inout) :: self
    double precision :: omega, c_start, prev_rslt, grad
    double precision :: c, u
    logical :: first_flag
    integer :: i
    
    prev_rslt = 0.d0
    c_start = self%cmin
    first_flag = .true.
    do i = 1, self%nf
       omega = 2.d0 * pi * (self%fmin  + (i - 1) * self%df)
       if (.not. self%full_calculation) then
          ! find root
          if (i /= 1) then
             first_flag = .false.
             if (self%u(i-1) /= 0.d0) then
                grad = (1.d0 - self%c(i-1) / self%u(i-1)) * &
                     & self%c(i-1) / (omega - 2.d0 * pi * self%df)
                c_start = self%c(i-1) + &
                     & 3.5d0 * grad * 2.d0 * pi * self%df
                !c_start = self%cmin
                if (c_start <= self%cmin) then
                   c_start = self%cmin
                end if
             end if
          else
             c_start = self%cmin
          end if
          call self%find_root(omega, c_start, c, u, first_flag)
          if (c == 0.d0) then
             self%c(:) = 0.d0
             self%u(:) = 0.d0
             write(*,*)"END dispersion calculation"
             return
          end if
          self%c(i) = c
          self%u(i) = u
          if (self%out_flag) then
             write(self%io, '(3F10.4)') &
                  & omega / (2.d0 * pi), self%c(i), self%u(i)
          end if
       else
          ! full calculation
          call self%do_full_calculation(omega)
       end if
    end do

    
    return 
  end subroutine disper_dispersion
  
  !---------------------------------------------------------------------
  
  subroutine disper_do_full_calculation(self, omega)
    class(disper), intent(inout) :: self
    double precision, intent(in) :: omega
    integer :: i
    double precision :: c_tmp, rslt

    do i = 1, self%nc
       c_tmp = self%cmin + (i - 1) * self%dc 
       call self%do_propagation(omega, c_tmp, rslt)
       write(self%io, *) omega / (2.d0 * pi), c_tmp, rslt
    end do
    
    return 
  end subroutine disper_do_full_calculation

  !---------------------------------------------------------------------
  
  subroutine disper_find_root(self, omega, c_start, c, u, first_flag)
    class(disper), intent(inout) :: self
    double precision, intent(in) :: omega, c_start
    double precision, intent(out) :: c, u
    logical, intent(in) :: first_flag
    integer :: i, count
    double precision :: c_tmp, rslt, prev_rslt, cmin2, cmax2
    double precision :: c1, c2, rslt1_c, rslt2_c, del_c
    double precision :: rslt_min, rslt_max
    double precision :: omega1, omega2, del_omega, rslt1_omg, rslt2_omg
    logical :: is_found

    is_found = .false.

    ! First step 
    prev_rslt = 0.d0
    count = -1
    do i = 1, self%nc
       c_tmp = c_start + (i - 1) * self%dc
       !write(*,*)"C_tmp=", c_tmp
       if (c_tmp > self%cmax) exit
       !write(*,*)"First"
       call self%do_propagation(omega, c_tmp, rslt)
       if (prev_rslt * rslt < 0) then
          count = count + 1
          if (.not. first_flag .or. self%n_mode == count) then
             is_found = .true.
             prev_rslt = rslt
             cmin2 = c_tmp - self%dc
             cmax2 = c_tmp
             exit
          end if
       end if
       prev_rslt = rslt
    end do
    if (.not. is_found) then
       write(0,*)"Warning: root is not found (disper_find_root)", &
            & " frequency = ", omega / 2.0 / pi, &
            & ", c_tmp = ", c_tmp, ", n_mode = ", &
            & self%n_mode
       c = 0.d0
       u = 0.d0
       return
    end if

    !c1 = c - 0.001d0 * self%dc 
    !call self%do_propagation(omega, c1, rslt1_c)
    !c2 = c + 0.001d0 * self%dc 
    !call self%do_propagation(omega, c2, rslt2_c)
    !del_c = (rslt2_c - rslt1_c) / (c2 - c1)
    
    ! Second step
    do i = 1, self%niter
       c_tmp = 0.5d0 * (cmin2 + cmax2)
       !write(*,*)"Second"             
       call self%do_propagation(omega, c_tmp, rslt)
       if (prev_rslt * rslt > 0.d0) then
          cmax2 = c_tmp
          rslt_max = rslt
       else
          cmin2 = c_tmp
          rslt_min = rslt
       end if
       !if (rslt / prev_rslt < abs(del_c * eps)) then
       !   write(*,*)"iter: ", i
       !   exit
       !end if
    end do
    
    c = c_tmp
    ! Third step 
    !c = (rslt_min * cmax2 + rslt_max * cmin2) / (rslt_min + rslt_max)
    !write(*,*)"Thrid", rslt_min, rslt_max, c, cmin2, cmax2
    !call self%do_propagation(omega, c, rslt)
    
    ! Group velocity
    ! derivative by c
    c1 = c - 0.001d0 * self%dc 
    call self%do_propagation(omega, c1, rslt1_c)
    c2 = c + 0.001d0 * self%dc 
    call self%do_propagation(omega, c2, rslt2_c)
    del_c = (rslt2_c - rslt1_c) / (c2 - c1)
    ! derivative by omega
    omega1 = omega - 0.001d0 * self%df * 2.d0 * pi 
    call self%do_propagation(omega1, c, rslt1_omg)
    omega2 = omega + 0.001d0 * self%df * 2.d0 * pi 
    call self%do_propagation(omega2, c, rslt2_omg)
    del_omega = (rslt2_omg - rslt1_omg) / (omega2 - omega1)
    ! group velocity from the two derivatives
    u = c / (1.d0 + omega * del_omega / (c * del_c))
    !write(*,*)"C=", c, "U=",u
    return
  end subroutine disper_find_root
  !---------------------------------------------------------------------

  subroutine disper_do_propagation(self, omega, c, rslt)
    class(disper), intent(inout) :: self
    double precision, intent(in) :: omega, c
    double precision, intent(out) :: rslt
    integer :: i, nlay, i0
    
    nlay = self%vmodel%get_nlay()

    if (self%is_ocean) then
       i0 = 2
    else
       i0 = 1
    end if
    
    call self%init_propagation(c)
    do i = nlay-1, i0, -1
       self%y = matmul(self%solid_propagator(i, omega, c), self%y)
    end do

    if (self%is_ocean) then
       self%y(1) = self%y(3)
       self%y(2) = self%y(5)
       self%y(1:2) = matmul(self%liquid_propagator(omega, c), self%y(1:2))
       rslt = self%y(2)
    else
       rslt = self%y(i24)
    end if
    
    return 
  end subroutine disper_do_propagation
  
  !---------------------------------------------------------------------
  
  subroutine disper_init_propagation(self, c)
    class(disper), intent(inout) :: self
    integer :: nlay
    double precision, intent(in) :: c

    double precision :: vp, vs, rho, ra, rb, gm

    ! get velocity of model bottom
    nlay = self%vmodel%get_nlay()
    vp = self%vmodel%get_vp(nlay)
    vs = self%vmodel%get_vs(nlay)
    rho = self%vmodel%get_rho(nlay)
    
    if (c > vp .or. c > vs) then
       write(0,*)"ERROR: phase velocity must be lower than Vp or Vs "
       write(0,*)"     : at the model bottom (disper_init_propagation)"
       write(0,*)"     : vp=", vp, "vs=", vs, "c=", c
       stop
    end if
    ra = sqrt(1.d0 - (c / vp) ** 2)
    rb = sqrt(1.d0 - (c / vs) ** 2)
    gm = 2.d0 * (vs / c) ** 2
    
    self%y(i13) = ra * rb - 1.d0 
    self%y(i14) = -rho * ra       
    self%y(i23) = -rho * rb      
    self%y(i12) = rho * (gm * self%y(i13) + 1.d0) 
    self%y(i24) = -rho * (gm * self%y(i12) + rho * (gm - 1.d0)) 
    
    return 
  end subroutine disper_init_propagation
  
  !---------------------------------------------------------------------
  
  function disper_solid_propagator(self, i, omega, c) result(p)
    class(disper), intent(inout) :: self
    integer, intent(in) :: i
    double precision, intent(in) :: omega, c
    double precision :: p(5, 5)
    double precision :: vp, vs, rho, h
    double precision :: k
    double precision :: ca, sa, cb, sb, ra2, rb2, fac_a, fac_b, fac, gm

    vp = self%vmodel%get_vp(i)
    vs = self%vmodel%get_vs(i)
    rho = self%vmodel%get_rho(i)
    h = self%vmodel%get_h(i)
    
    k = omega / c
    gm = 2.d0 * (vs / c)**2
    call calc_hyp_trig(c, vp, h, k, ca, sa, ra2, fac_a)
    call calc_hyp_trig(c, vs, h, k, cb, sb, rb2, fac_b)
    fac = fac_a * fac_b

    
    p(i12, i13) = &
         & rho* (gm * (gm - 1.d0) * (2.d0 * gm - 1.d0) * &
         & (ca * cb - fac) - ((gm - 1.d0)**3 + gm**3 * ra2 * rb2) * sa * sb)
    p(i12, i14) = &
         & (gm - 1.d0) * sa * cb - gm * rb2 * ca * sb
    p(i12, i23) = &
         & -gm * ra2 * sa * cb + (gm - 1.d0) * ca * sb
    p(i12, i24) = &
         & (-(2.d0 * gm - 1.d0) * (ca * cb - fac) + &
         & ((gm - 1.d0) + gm * ra2 * rb2) * sa* sb) / rho
    p(i13, i13) = &
         & fac + ((gm - 1.d0)**2 + gm**2) * (ca * cb - fac) - &
         & ((gm - 1.d0)**2 + gm**2 * ra2 * rb2) * sa * sb
    p(i13, i14) = &
         & (sa * cb - rb2 * ca * sb) / rho
    p(i13, i23) = &
         & (ca * sb - ra2 * sa * cb) / rho
    p(i13, i24) = &
         & (-2.d0 * (ca * cb - fac) + &
         & (1.d0 + ra2 * rb2) * sa * sb) / rho**2 
    p(i14, i13) = &
         rho* (gm**2 * ra2 * sa *cb - (gm - 1.d0)**2 * ca * sb) 
    p(i14, i14) = ca * cb 
    p(i14, i23) = ra2 * sa * sb
    p(i23, i13) = &
         & rho* (-(gm - 1.d0)**2 * sa * cb + gm**2 * rb2 * ca * sb) 
    p(i23, i14) = rb2 * sa * sb
    p(i24, i13) = &
         & rho**2 * (-2.d0 * gm**2 * (gm - 1.d0)**2 * (ca * cb - fac) + &
         & ((gm - 1.d0)**4 + gm**4 * ra2 * rb2) * sa * sb) 
    p(i12, i12) = &
         & fac - 2.d0 * (2.d0 * gm * (gm - 1.d0) * (ca * cb - fac) - &
         & ((gm - 1.d0)**2 + gm**2 * ra2 * rb2) * sa * sb)
    p(i13, i12) = 2.d0 * p(i12, i24)
    p(i14, i12) = 2.d0 * p(i12, i23)
    p(i14, i24) = p(i13, i23)
    p(i23, i12) = 2.d0 * p(i12, i14)
    p(i23, i23) = p(i14, i14)
    p(i23, i24) = p(i13, i14)
    p(i24, i12) = 2.d0 * p(i12, i13)
    p(i24, i14) = p(i23, i13)
    p(i24, i23) = p(i14, i13)
    p(i24, i24) = p(i13, i13)
    

    return
  end function disper_solid_propagator

  !---------------------------------------------------------------------
  
  function disper_liquid_propagator(self, omega, c) result(p)
    class(disper), intent(inout) :: self
    double precision, intent(in) :: omega, c
    double precision :: vp, rho, h, k
    double precision :: ca, sa, ra2, fac
    double precision :: p(2, 2)
    
    vp = self%vmodel%get_vp(1)
    rho = self%vmodel%get_rho(1)
    h = self%vmodel%get_h(1)
    
    k = omega / c
    call calc_hyp_trig(c, vp, h, k, ca, sa, ra2, fac)

    p(1, 1) = ca
    p(1, 2) = -ra2 * sa / rho
    p(2, 1) = - rho * sa
    p(2, 2) = ca
    
    return 
  end function disper_liquid_propagator

  !---------------------------------------------------------------------
  
  subroutine calc_hyp_trig(c, v, h, k, ch, sh, r2, fac)
    real(8), intent(in) :: c, v, h, k
    real(8), intent(out) :: sh, ch, r2, fac
    real(8) :: r
    
    
    if (v > c) then
       r = sqrt(1.d0 - (c/v)**2)
       r2 = r * r
       !if (r * h * k > 0.01d0) then
       ch = 1.d0
       sh = tanh(r * h * k) / r
       if (r * h * k <= huge) then
          fac = 1.d0 / cosh(r * h * k)
       else
          fac = 0.d0
       end if
       !else
       !   ch = 1.d0
       !   sh = 0.d0
       !   fac = 1.d0
       !end if
    else if (v < c) then
       ! Use cos and sin instead of cosh and sinh 
       !   to avoid complex values
       r = sqrt((c/v)**2 - 1.d0)
       r2 = - r * r
       ch = cos(r * h * k) 
       sh = sin(r * h * k) / r 
       fac = 1.d0
    else
       r2 = 0.d0
       ch = 1.d0
       sh = h * k
       fac = 1.d0
    end if
    
    return 
  end subroutine calc_hyp_trig
  
  !---------------------------------------------------------------------

  subroutine disper_set_full_calculation(self, flag)
    class(disper), intent(inout) :: self
    logical, intent(in), optional :: flag

    if (present(flag)) then
       self%full_calculation = flag
    else
       self%full_calculation = .true.
    end if
    
    return 
  end subroutine disper_set_full_calculation

  !---------------------------------------------------------------------

  subroutine disper_set_vmodel(self, vm)
    class(disper), intent(inout) :: self
    type(vmodel), intent(in) :: vm
    
    self%vmodel = vm

    return 
  end subroutine disper_set_vmodel

  !---------------------------------------------------------------------

  integer function disper_get_nf(self) result(nf)
    class(disper), intent(in) :: self

    nf = self%nf
    
    return 
  end function disper_get_nf

  !---------------------------------------------------------------------

  integer function disper_get_nc(self) result(nc)
    class(disper), intent(in) :: self

    nc = self%nc
    
    return 
  end function disper_get_nc

  !---------------------------------------------------------------------

  double precision function disper_get_c(self, i) result(c)
    class(disper), intent(in) :: self
    integer, intent(in) :: i
    
    c = self%c(i)

    return 
  end function disper_get_c

  !---------------------------------------------------------------------
  
  function disper_get_c_array(self) result(c)
    class(disper), intent(in) :: self
    double precision :: c(self%nf)
  
    c(:) = self%c(:)

    return 
  end function disper_get_c_array

  !---------------------------------------------------------------------
  
  double precision function disper_get_u(self, i) result(u)
    class(disper), intent(in) :: self
    integer, intent(in) :: i
    
    u = self%u(i)

    return 
  end function disper_get_u

  !---------------------------------------------------------------------

  function disper_get_u_array(self) result(u)
    class(disper), intent(in) :: self
    double precision :: u(self%nf)
  
    u(:) = self%u(:)

    return 
  end function disper_get_u_array

  !---------------------------------------------------------------------
  
  subroutine disper_save_syn(self)
    class(disper), intent(inout) :: self
    integer :: i, j
    
    do i = 1, self%nf
       j = int((self%c(i) - self%cmin) / self%dc) + 1
       if (j < 1 .or. j > self%nc) cycle
       self%n_fc(i, j) = self%n_fc(i, j) + 1
       
       j = int((self%u(i) - self%cmin) / self%dc) + 1
       if (j < 1 .or. j > self%nc) cycle
       self%n_fu(i, j) = self%n_fu(i, j) + 1
    end do
    
    return 
  end subroutine disper_save_syn

  !---------------------------------------------------------------------

  function disper_get_n_fc(self) result(n_fc)
    class(disper), intent(in) :: self
    integer :: n_fc(self%nf, self%nc)
    
    n_fc(:,:) = self%n_fc(:,:)

    return 
  end function disper_get_n_fc
  
  !---------------------------------------------------------------------

  function disper_get_n_fu(self) result(n_fu)
    class(disper), intent(in) :: self
    integer :: n_fu(self%nf, self%nc)
    
    n_fu(:,:) = self%n_fu(:,:)

    return 
  end function disper_get_n_fu
  
  !---------------------------------------------------------------------

end module cls_disper
