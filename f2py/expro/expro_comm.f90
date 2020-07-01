subroutine expro_icomm(p)

  use mpi

  implicit none

  integer, intent(inout) :: p
  integer :: ierr,iproc
  logical :: flag, file_exists
  
  call MPI_INITIALIZED(flag,ierr)

  ! Check if running QLGYRO - if directly read in data
  ! Otherwise have to wait for proc 0 to reach here
  inquire(file='out.qlgyro.status', exist=file_exists)
  
  if (file_exists) flag = .false.
  
  if (flag) then
     call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
     if (iproc == 0) read(1,*) p
     call MPI_BCAST(p,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  else
     read(1,*) p
  endif

end subroutine expro_icomm

subroutine expro_rcomm(x)

  use mpi

  implicit none

  double precision, intent(inout) :: x
  integer :: ierr,iproc
  logical :: flag, file_exists
  
  call MPI_INITIALIZED(flag,ierr)
  
  inquire(file='out.qlgyro.status', exist=file_exists)
  if (file_exists) flag = .false.
  
  if (flag) then
     call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
     if (iproc == 0) read(1,10) x
     call MPI_BCAST(x,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  else
     read(1,10) x
  endif

10 format(1pe14.7)

end subroutine expro_rcomm

subroutine expro_scomm(x,n)

  use mpi

  implicit none

  integer, intent(in) :: n
  double precision, intent(inout), dimension(n) :: x
  integer :: ierr,iproc
  logical :: flag, file_exists
  
  call MPI_INITIALIZED(flag,ierr)

  inquire(file='out.qlgyro.status', exist=file_exists)
  if (file_exists) flag = .false.
  
  if (flag) then
     call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
     if (iproc == 0) read(1,*) x
     call MPI_BCAST(x,n,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  else
     read(1,*) x
  endif
  
end subroutine expro_scomm

subroutine expro_lcomm(x,n)

  use mpi

  implicit none

  integer, intent(in) :: n
  double precision, intent(inout), dimension(n) :: x
  integer :: ierr,iproc
  logical :: flag, file_exists
  
  call MPI_INITIALIZED(flag,ierr)
  
  inquire(file='out.qlgyro.status', exist=file_exists)
  if (file_exists) flag = .false.
  
  if (flag) then
     call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
     if (iproc == 0) read(1,10) x
     call MPI_BCAST(x,n,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  else
     read(1,10) x
  endif

10 format(10(1pe14.7))

end subroutine expro_lcomm

subroutine expro_tcomm(x,n)

  use mpi

  implicit none

  integer, intent(in) :: n
  character*10, intent(inout), dimension(20) :: x
  integer :: ierr,iproc
  logical :: flag, file_exists
  
  call MPI_INITIALIZED(flag,ierr)

  inquire(file='out.qlgyro.status', exist=file_exists)
  if (file_exists) flag = .false.
  
  if (flag) then
     call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
     if (iproc == 0) read(1,*) x(1:n)
     call MPI_BCAST(x(1:n),n*10,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  else
     read(1,*) x(1:n)
  endif

end subroutine expro_tcomm

subroutine expro_vcomm(x,n)

  use mpi

  implicit none

  integer :: idum,i
  integer, intent(in) :: n
  double precision, intent(inout), dimension(n) :: x
  integer :: ierr,iproc
  logical :: flag, file_exists

  call MPI_INITIALIZED(flag,ierr)
  
  inquire(file='out.qlgyro.status', exist=file_exists)
  if (file_exists) flag = .false.
  
  if (flag) then
     call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
     if (iproc == 0) then
        do i=1,n
           read(1,10) idum,x(i)
        enddo
     endif
     call MPI_BCAST(x,n,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  else
     do i=1,n
        read(1,10) idum,x(i)
     enddo
  endif

10 format(i3,1x,1pe14.7)

end subroutine expro_vcomm

subroutine expro_acomm(x,m,n)

  use mpi

  implicit none

  integer :: idum,i
  integer, intent(in) :: m,n
  double precision, intent(inout), dimension(m,n) :: x
  integer :: ierr,iproc
  logical :: flag, file_exists

  call MPI_INITIALIZED(flag,ierr)

  inquire(file='out.qlgyro.status', exist=file_exists)
  if (file_exists) flag = .false.
  
  if (flag) then
     call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
     if (iproc == 0) then
        do i=1,n
           read(1,10) idum,x(:,i)
        enddo
     endif
     call MPI_BCAST(x,m*n,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  else
     do i=1,n
        read(1,10) idum,x(:,i)
     enddo
  endif

10 format(i3,1x,10(1pe14.7,1x))

end subroutine expro_acomm

