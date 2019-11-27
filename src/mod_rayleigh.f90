module mod_rayleigh
  use mod_vmodel
  implicit none 
  
  double precision, private, parameter :: pi = acos(-1.d0)
  double precision, private, parameter :: eps = 1.0d-10
  double precision, private, parameter :: huge = 100.d0
  integer, private, parameter :: i12 = 1, i13 = 2, i14 = 3
  integer, private, parameter :: i23 = 4, i24 = 5


  type rayleigh
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

     integer :: niter = 7

     integer :: io

     logical :: full_calculation = .false.

     double precision, allocatable :: c(:) ! phase velocity
     double precision, allocatable :: u(:) ! group velocity
     double precision :: y(5)
     logical :: is_ocean
   contains
     procedure :: dispersion => rayleigh_dispersion
     procedure :: init_propagation => rayleigh_init_propagation
     procedure :: do_propagation => rayleigh_do_propagation
     procedure :: solid_propagator => rayleigh_solid_propagator
     procedure :: liquid_propagator => rayleigh_liquid_propagator
     procedure :: find_root => rayleigh_find_root
     procedure :: set_full_calculation => rayleigh_set_full_calculation
     procedure :: do_full_calculation => rayleigh_do_full_calculation
     procedure :: set_vmodel => rayleigh_set_vmodel
  end type rayleigh

  interface rayleigh
     module procedure init_rayleigh
  end interface rayleigh
     
contains
  
  !---------------------------------------------------------------------
  
  type(rayleigh) function init_rayleigh(vm, fmin, fmax, df, &
       & cmin, cmax, dc)
    type(vmodel), intent(in) :: vm
    double precision, intent(in) :: fmin, fmax, df
    double precision, intent(in) :: cmin, cmax, dc
    integer :: nlay
    
    ! velocity model
    init_rayleigh%vmodel = vm
    nlay = init_rayleigh%vmodel%get_nlay()
    if (nlay < 1) then
       write(0,*)"ERROR: velocity model is not defined (mod_rayleigh)"
       stop
    end if
    if (init_rayleigh%vmodel%get_vs(1) < 0.d0) then
       init_rayleigh%is_ocean = .true.
    else
       init_rayleigh%is_ocean = .false.
    end if
    
    init_rayleigh%fmin = fmin
    init_rayleigh%fmax = fmax
    init_rayleigh%df = df
    init_rayleigh%nf = nint((fmax - fmin) / df) + 1
    if (init_rayleigh%nf < 1) then
       write(0,*)"ERROR: fmax must be .gt. fmin (init_rayleigh)"
       stop
    end if
    allocate(init_rayleigh%c(init_rayleigh%nf), &
         &   init_rayleigh%u(init_rayleigh%nf))
    
    init_rayleigh%cmin = cmin
    init_rayleigh%cmax = cmax
    init_rayleigh%dc = dc
    init_rayleigh%nc = nint((cmax - cmin) / dc) + 1
    if (init_rayleigh%nc < 1) then
      write(0,*)"ERROR: cmax must be .gt. cmin (init_rayleigh)"
      stop
   end if

   ! output file
   open(newunit = init_rayleigh%io, status = "unknown", &
        & file = "rayleigh.out")
   


    return 
  end function init_rayleigh
    
  !---------------------------------------------------------------------

  subroutine rayleigh_dispersion(self)
    class(rayleigh), intent(inout) :: self
    double precision :: omega, c_next, c_start
    logical :: is_first
    integer :: i
    
    is_first = .true.
    c_start = self%cmin
    do i = 1, self%nf
       omega = 2.d0 * pi * (self%fmin  + (i - 1) * self%df)
       if (.not. self%full_calculation) then
          ! find root
          call self%find_root(omega, c_start, is_first, &
               & self%c(i), self%u(i), c_next)
          write(self%io, *) omega / (2.d0 * pi), self%c(i), self%u(i)
          is_first = .false.
          c_start = c_next
       else
          ! full calculation
          call self%do_full_calculation(omega)
       end if

       
    end do

    
    return 
  end subroutine rayleigh_dispersion
  
  !---------------------------------------------------------------------
  
  subroutine rayleigh_do_full_calculation(self, omega)
    class(rayleigh), intent(inout) :: self
    double precision, intent(in) :: omega
    integer :: i
    double precision :: c_tmp, rslt

    do i = 1, self%nc
      c_tmp = self%cmin + (i - 1) * self%dc 
      call self%do_propagation(omega, c_tmp, rslt)
      write(self%io, *) omega / (2.d0 * pi), c_tmp, rslt
      
    end do
    
    
    return 
  end subroutine rayleigh_do_full_calculation

  !---------------------------------------------------------------------
  
  subroutine rayleigh_find_root(self, omega, c_start, is_first, &
       & c, u, c_next)
    class(rayleigh), intent(inout) :: self
    double precision, intent(in) :: omega, c_start
    double precision, intent(out) :: c, u, c_next
    logical, intent(in) :: is_first
    integer :: i
    double precision :: c_tmp, rslt, prev_rslt, cmin2, cmax2
    double precision :: c1, c2, rslt1, rslt2, del_c
    double precision :: rslt_min, rslt_max
    double precision :: omega1, omega2, del_omega
    logical :: is_found

    is_found = .false.

    ! First step
    prev_rslt = 0.d0
    do i = 1, self%nc
       if (is_first) then ! increase c (phase velocity)
          c_tmp = c_start + (i - 1) * self%dc
       else               ! decrease in c (phase velocity)
          c_tmp = c_start - (i - 1) * self%dc
       end if
       call self%do_propagation(omega, c_tmp, rslt)
       if (prev_rslt * rslt < 0) then
          is_found = .true.
          prev_rslt = rslt
          if (is_first) then
             cmin2 = c_tmp - self%dc
             cmax2 = c_tmp
          else
             cmin2 = c_tmp
             cmax2 = c_tmp + self%dc
          end if
          exit
       end if
       prev_rslt = rslt
    end do
    if (.not. is_found) then
       write(0,*)"Warning: root is not found (rayleigh_find_root)", &
            & "       : omega = ", omega
       c = -999.d0
       u = -999.d9
       return
    end if
    
       
    ! Second step
    do i = 1, self%niter
       c_tmp = 0.5d0 * (cmin2 + cmax2)
       call self%do_propagation(omega, c_tmp, rslt)
       if (.not. is_first) then
          if (prev_rslt * rslt > 0.d0) then
             cmin2 = c_tmp
             rslt_min = rslt
          else
             cmax2 = c_tmp
             rslt_max = rslt
          end if
       else
          if (prev_rslt * rslt > 0.d0) then
             cmax2 = c_tmp
             rslt_max = rslt
          else
             cmin2 = c_tmp
             rslt_min = rslt
          end if

       end if
    end do
    c = (rslt_min * cmax2 + rslt_max * cmax2) / (rslt_min + rslt_max)
    c_next = c 
    
    ! Group velocity
    c1 = c - 0.001d0 * self%dc 
    call self%do_propagation(omega, c1, rslt1)
    c2 = c + 0.001d0 * self%dc 
    call self%do_propagation(omega, c2, rslt2)
    del_c = (rslt2 - rslt1) / (c2 - c1)

    omega1 = omega - 0.001d0 * self%df * 2.d0 * pi 
    call self%do_propagation(omega1, c, rslt1)
    omega2 = omega + 0.001d0 * self%df * 2.d0 * pi 
    call self%do_propagation(omega2, c, rslt2)
    del_omega = (rslt2 - rslt1) / (omega2 - omega1)
    
    write(*,*)c, omega, del_omega, del_c, rslt2, rslt1
    u = c / (1.d0 + omega * del_omega / (c * del_c))
    

    return
  end subroutine rayleigh_find_root
  !---------------------------------------------------------------------

  subroutine rayleigh_do_propagation(self, omega, c, rslt)
    class(rayleigh), intent(inout) :: self
    double precision, intent(in) :: omega, c
    double precision, intent(out) :: rslt
    double precision :: maxamp
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
       maxamp = maxval(abs(self%y))
       self%y(1:5) = self%y(1:5) / maxamp
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
  end subroutine rayleigh_do_propagation
  
  !---------------------------------------------------------------------
  
  subroutine rayleigh_init_propagation(self, c)
    class(rayleigh), intent(inout) :: self
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
       write(0,*)"     : at the model bottom (rayleigh_init_propagation)"
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

    self%y(1:5) = self%y(1:5) * eps
    
    return 
  end subroutine rayleigh_init_propagation
  
  !---------------------------------------------------------------------
  
  function rayleigh_solid_propagator(self, i, omega, c) result(p)
    class(rayleigh), intent(inout) :: self
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
  end function rayleigh_solid_propagator

  !---------------------------------------------------------------------
  
  function rayleigh_liquid_propagator(self, omega, c) result(p)
    class(rayleigh), intent(inout) :: self
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
  end function rayleigh_liquid_propagator

  !---------------------------------------------------------------------
  
  subroutine calc_hyp_trig(c, v, h, k, ch, sh, r2, fac)
    real(8), intent(in) :: c, v, h, k
    real(8), intent(out) :: sh, ch, r2, fac
    real(8) :: r
    
    !write(*,*)c, v, h, k
    if (v > c) then
       r = sqrt(1.d0 - (c/v)**2)
       r2 = r * r
       ch = 1.d0
       sh = tanh(r * h * k) / r
       if (r * h * k <= huge) then
          fac = 1.d0 / cosh(r * h * k)
       else
          fac = 0.d0
       end if
    else 
       ! Use cos and sin instead of cosh and sinh 
       !   to avoid complex values
       r = sqrt((c/v)**2 - 1.d0)
       if (r == 0.d0) r = eps ! to avoid zero division
       r2 = - r * r
       ch = cos(r * h * k) 
       sh = sin(r * h * k) / r 
       fac = 1.d0
    end if
    
    return 
  end subroutine calc_hyp_trig
  
  !---------------------------------------------------------------------

  subroutine rayleigh_set_full_calculation(self, flag)
    class(rayleigh), intent(inout) :: self
    logical, intent(in), optional :: flag

    if (present(flag)) then
       self%full_calculation = flag
    else
       self%full_calculation = .true.
    end if
    
    return 
  end subroutine rayleigh_set_full_calculation

  !---------------------------------------------------------------------

  subroutine rayleigh_set_vmodel(self, vm)
    class(rayleigh), intent(inout) :: self
    type(vmodel), intent(in) :: vm
    
    self%vmodel = vm

    return 
  end subroutine rayleigh_set_vmodel

end module mod_rayleigh
