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
module cls_line_text
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
     procedure :: get_line => line_text_get_line
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
    if (self%ignore_space) then
       call self%remove_space()
    end if
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
       val = trim(adjustl(self%line(j+1:self%nlen)))
       is_ok = .true.
    end if
    
    return 
  end subroutine line_text_read_value

  !---------------------------------------------------------------------
  
  character(line_max) function line_text_get_line(self) result(line)
    class(line_text), intent(in) :: self
    
    line = self%line
    
    return 
  end function line_text_get_line

  !-----------------------------------------------------------------------

end module cls_line_text
