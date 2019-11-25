program main
  use mod_random
  include 'mpif.h'
  integer :: i, n, io, ierr, rank
  integer :: iseed1, iseed2, iseed3, iseed4
  character(50) :: outfile
  
  iseed1 = 12345678
  iseed2 = 77777777
  iseed3 = 10101010
  iseed4 = 99887766

  call mpi_init(ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)

  call init_random(iseed1, iseed2, iseed3, iseed4, rank)
  n = 100000
  write(outfile, '(a,I2.2,a)') "test_random_", rank, ".out"
  open(newunit=io, file = outfile, status = "unknown")
  do i = 1, n 
     !write(*,*)i
     write(io, *)rand_g()
  end do
  close(io)
  
  call mpi_barrier(ierr)
  
  call mpi_finalize(ierr)

  stop
end program main
  
