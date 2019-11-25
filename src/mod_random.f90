module mod_random
  implicit none 
  
  double precision, parameter, private :: r_max = 2.d0 ** 32
  integer, private :: x, y, z, w
  integer, parameter, private :: k = 8
  integer, parameter, private :: n = 2**k
  double precision, parameter, private :: r = 3.6541528853610088d0
  double precision, parameter, private :: v = 0.00492867323399d0
  double precision, private :: xg(0:n)
  double precision, parameter :: pi2 = 2.d0 * acos(-1.d0)

contains

  !---------------------------------------------------------------------

  subroutine init_random(i1, i2, i3, i4, rank)
    integer, intent(in) :: i1, i2, i3, i4
    integer, intent(in), optional :: rank
    integer :: j

    if (present(rank)) then
       j = rank
    else
       j = 0
    end if
    x = i1 * (j + 1)**4 + 1000 * i1 * (j + 1) ** 2 + i1
    y = i2 * (j + 1)**4 + 1000 * i2 * (j + 1) ** 2 + i2
    z = i3 * (j + 1)**4 + 1000 * i3 * (j + 1) ** 2 + i3
    w = i4 * (j + 1)**4 + 1000 * i4 * (j + 1) ** 2 + i4

    return 
  end subroutine init_random

  !---------------------------------------------------------------------
  
  ! U[0, 1)
  double precision function rand_u()
    integer :: t
    
    

    t = ieor(x, ishft(x, 11))
    x = y
    y = z
    z = w
    w = &
         & ieor( &
         & ieor(w, ishft(w, -19)), &
         & ieor(t, ishft(t,  -8)) &
         & )
    if (w < 0) then
       rand_u = (dble(w) + r_max) / r_max
    else
       rand_u = dble(w) / r_max
    end if
    return 
  end function rand_u

  !---------------------------------------------------------------------
  
  !---------------------------------------------------------------------
  
  double precision function rand_g()
    double precision:: v1, v2
    
    v1 = rand_u() 
    v2 = rand_u()
    rand_g = sqrt(-2.d0 * log(v1)) * cos(pi2 * v2)

  end function rand_g
end module mod_random
  
