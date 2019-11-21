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

    return 
  end function init_rayleigh
    
  !---------------------------------------------------------------------

  subroutine rayleigh_dispersion(self)
    class(rayleigh), intent(inout) :: self
    double precision :: omega, c, rslt
    integer :: i, j
    
    do i = 1, self%nf
       omega = 2.d0 * pi * (self%fmin  + (i - 1) * self%df)
       write(*,*) omega / (2.d0 * pi)
       do j = 1, self%nc
          c = self%cmin + (j - 1) * self%dc
          call self%do_propagation(omega, c, rslt)
          write(111,*)c, omega/ (2.d0 * pi), rslt
       end do
    end do

    
    return 
  end subroutine rayleigh_dispersion
  
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
    else if (v < c) then
       ! Use cos and sin instead of cosh and sinh 
       !   to avoid complex values
       r = sqrt((c/v)**2 - 1.d0)
       r2 = - r * r
       ch = cos(r * h * k) 
       sh = sin(r * h * k) / r 
       fac = 1.d0
    else
       r = 0.d0
       r2 = 0.d0
       ch = 1.d0
       sh = 0.d0
       fac = 1.d0
    end if
    
    return 
  end subroutine calc_hyp_trig

end module mod_rayleigh
