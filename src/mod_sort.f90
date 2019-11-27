module mod_sort
  implicit none
  
  public quick_sort
  private swap
  
contains
  
  !---------------------------------------------------------------------

  recursive subroutine quick_sort(a, il, ir, b, c)
    implicit none 
    integer, intent(in) :: il, ir
    real(8), intent(inout) :: a(:), b(:), c(:)
    integer :: ipiv, i, j
    real(8) :: piv
    
    if (ir - il <= 0) return
    
    ipiv = (il + ir) / 2
    piv  = a(ipiv)
    
    call swap(a, ipiv, ir)
    call swap(b, ipiv, ir)
    call swap(c, ipiv, ir)
    
    i = il
    do j = il, ir
       if (a(j) < piv) then
          call swap(a, i, j)
          call swap(b, i, j)
          call swap(c, i, j)
          i = i + 1
       end if
    end do
    
    call swap(a, i, ir)
    call swap(b, i, ir)
    call swap(c, i, ir)
    
    call quick_sort(a, il, i, b, c)
    call quick_sort(a, i + 1, ir, b, c)
    
    return 
  end subroutine quick_sort
  
  !---------------------------------------------------------------------
  
  subroutine swap(a, i1, i2)
    implicit none
    real(8), intent(inout) :: a(:)
    integer, intent(in) :: i1, i2
    real(8) :: tmp
    
    tmp = a(i1)
    a(i1) = a(i2)
    a(i2) = tmp
    
    return 
  end subroutine swap

  !---------------------------------------------------------------------
  
end module mod_sort
  
