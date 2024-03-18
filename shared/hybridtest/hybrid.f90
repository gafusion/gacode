program hybrid
  
  use mpi

  implicit none
 
  integer :: ierr,iproc
  integer :: nmpi,nomp,ntot
  integer :: i,j
  integer,external :: omp_get_max_threads
  integer,external :: omp_get_thread_num
  integer :: supported

  real :: mysum,mysums,total
  real, dimension(3,3) :: a
  real, dimension(3) :: b
  integer, dimension(3) :: work
  
  nomp = omp_get_max_threads()

  if (nomp > 1) then
     call MPI_INIT_THREAD(MPI_THREAD_FUNNELED,supported,ierr)
     if (supported < MPI_THREAD_FUNNELED) then
        print *, "ERROR: Multi-threaded MPI not supported." 
     endif
  else 
     call MPI_INIT_THREAD(MPI_THREAD_SINGLE,supported,ierr)
  endif

  call MPI_COMM_SIZE(MPI_COMM_WORLD,nmpi,ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)

  mysum = 0.0
  
!$omp parallel do private (i,a,b,work,ierr) reduction(+:mysum)
  do i=0,nomp-1
   
     a(:,:) = 0.0
     a(1,1) = 1.0
     a(2,2) = 1.0
     a(3,3) = 1.0
     
     b(:) = omp_get_thread_num()+iproc*nmpi
 
     call DGESV(3,1,a,3,work,b,3,ierr)

     mysum = mysum+sum(b)
  enddo

  call MPI_ALLREDUCE(mysum,mysums,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_FINALIZE(ierr)

  ntot = nomp*nmpi
  if (iproc == 0) then
     if (int(3*ntot*(ntot-1)/2-mysums) == 0) then
        print *,'Success'
     else
        print *,'Fail'
     endif
  endif
  
end program hybrid
