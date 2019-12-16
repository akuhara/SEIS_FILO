module mod_line_text
  implicit none 
  
  integer, parameter :: line_max = 200, symbol_max = 10
  
  type line_text
     private
     character(len=line_max) :: line
     integer :: nlen
     character(len=1) :: separator = "="
     character(len=1) :: comment_out = "#"
     logical :: ignore_space = .true.
     
   contains
     procedure :: remove_comment => line_text_remove_comment
     procedure :: remove_space  => line_text_remove_space
     procedure :: read_value => line_text_read_value
     

  end type line_text

  interface line_text
     module procedure init_line_text
  end interface line_text
  
contains

  !---------------------------------------------------------------------

  type(line_text) function init_line_text(line, separator, &
       & comment_out, ignore_space) result(self)
    character(line_max), intent(in) :: line
    character(1), intent(in), optional :: separator, comment_out
    logical, intent(in), optional :: ignore_space
    
    self%line = line
    if (present(comment_out)) then
       self%comment_out = comment_out
    end if
    if (present(separator)) then
       self%separator = separator
    end if
    if (present(ignore_space)) then
       self%ignore_space = ignore_space
    end if

    self%nlen = len_trim(line)

    !write(*,*)self%line
    call self%remove_comment()
    !write(*,*)self%line
    call self%remove_space()
    !write(*,*)self%line
    
    return 
  end function init_line_text

  !---------------------------------------------------------------------
  
  subroutine line_text_remove_comment(self)
    class(line_text), intent(inout) :: self
    integer :: i

    do i = 1, self%nlen
       if (self%line(i:i) == self%comment_out) then
          if (i > 1) then
             self%line = self%line(1:i-1)
             return
          else
             self%line = ""
          end if
       end if
    end do
    

    return 
  end subroutine line_text_remove_comment

  !---------------------------------------------------------------------

  subroutine line_text_remove_space(self)
    class(line_text), intent(inout) :: self
    integer :: i
    
    do i = 1, self%nlen
       if (self%line(i:i) == " ") then
          self%line = self%line(1:i-1) // self%line(i+1:self%nlen)
       end if
    end do
    
    return 
  end subroutine line_text_remove_space
    
  !---------------------------------------------------------------------
  
  subroutine line_text_read_value(self, name, val, is_ok)
    class(line_text), intent(in) :: self
    character(*), intent(out) :: name, val
    logical, intent(out) :: is_ok
    integer :: j
    
    j = index(self%line, self%separator)
    if (j == 0 .or. j == 1 .or. j == self%nlen) then
       name = "undef"
       val = "undef"
       is_ok = .false.
    else
       name = self%line(1:j-1)
       val = self%line(j+1:self%nlen)
       is_ok = .true.
    end if
    
    return 
  end subroutine line_text_read_value

  !---------------------------------------------------------------------
  


end module mod_line_text
