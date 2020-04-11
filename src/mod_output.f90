module mod_output
  use cls_parallel ! for including mpif.h
  implicit none 
contains
  
  !---------------------------------------------------------------------

  subroutine output_ppd_1d(filename, rank, n_bin_x, n_x, n_mod, x_min, dx)
    character(*), intent(in) :: filename
    integer, intent(in) :: rank, n_bin_x, n_mod
    integer, intent(in) :: n_x(n_bin_x)
    double precision, intent(in) :: x_min, dx
    double precision :: x
    integer :: io_x, ierr, i
    integer :: n_x_all(n_bin_x)
    
    call mpi_reduce(n_x, n_x_all, n_bin_x, MPI_INTEGER4, MPI_SUM, 0, &
         & MPI_COMM_WORLD, ierr)
    
    if (rank == 0) then
       open(newunit = io_x, file = filename, status = "unknown", &
            & iostat = ierr)
       if (ierr /= 0) then
          write(0, *)"ERROR: cannot open ", trim(filename)
          call mpi_finalize(ierr)
          stop
       end if
       
       
       do i = 1, n_bin_x
          x = x_min + (i - 1) * dx
          write(io_x, '(3F13.5)') x, dble(n_x_all(i)) / dble(n_mod)
       end do
       close(io_x)
    end if
    
    return
  end subroutine output_ppd_1d
  
  !---------------------------------------------------------------------
  
  subroutine output_ppd_2d(filename, rank, n_bin_x, n_bin_y, n_xy, n_mod, &
       & x_min, dx, y_min, dy)
    include 'mpif.h'
    character(*), intent(in) :: filename
    integer, intent(in) :: rank, n_bin_x, n_bin_y, n_mod
    integer, intent(in) :: n_xy(n_bin_x, n_bin_y)
    double precision, intent(in) :: x_min, dx, y_min, dy
    double precision :: x, y
    integer :: io_xy, ierr, i, j
    integer :: n_xy_all(n_bin_x, n_bin_y)
    
    call mpi_reduce(n_xy, n_xy_all, n_bin_x * n_bin_y, MPI_INTEGER4, &
         & MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    
    if (rank == 0) then
       open(newunit = io_xy, file = filename, status = "unknown", &
            & iostat = ierr)
       if (ierr /= 0) then
          write(0, *)"ERROR: cannot open ", trim(filename)
          call mpi_finalize(ierr)
          stop
       end if
       
       do i = 1, n_bin_y
          y = y_min + (i - 1) * dy
          do j = 1, n_bin_x
             x = x_min + (j - 1) * dx
             write(io_xy, '(3F13.5)') x, y, &
                  & dble(n_xy_all(j, i)) / dble(n_mod)
          end do
       end do
       close(io_xy)
    end if
    
    return
  end subroutine output_ppd_2d
  
  !---------------------------------------------------------------------
  
  subroutine output_mean_model(filename, rank, n_bin_x, sum_x, n_mod, &
       & x_min, dx)
    character(*), intent(in) :: filename
    integer, intent(in) :: rank, n_bin_x, n_mod
    double precision, intent(in) :: sum_x(n_bin_x)
    double precision, intent(in) :: x_min, dx
    double precision :: x
    integer :: io_x, ierr, i
    double precision :: sum_x_all(n_bin_x)
    
    call mpi_reduce(sum_x, sum_x_all, n_bin_x, MPI_DOUBLE_PRECISION, &
         & MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    
    if (rank == 0) then
       open(newunit = io_x, file = filename, status = "unknown", &
            & iostat = ierr)
       if (ierr /= 0) then
          write(0, *)"ERROR: cannot open ", trim(filename)
          call mpi_finalize(ierr)
          stop
       end if
       
       
       do i = 1, n_bin_x
          x = x_min + (i - 1) * dx
          write(io_x, '(3F13.5)') x, sum_x_all(i) / dble(n_mod)
       end do
       close(io_x)
    end if
    
    return
  end subroutine output_mean_model

  !---------------------------------------------------------------------

end module mod_output
