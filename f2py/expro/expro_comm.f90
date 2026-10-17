subroutine expro_icomm(p)

#ifndef EXPRO_NOMPI  
  use mpi
  use expro, only : expro_comm
#endif
  use expro, only : hasmpi

  implicit none

  integer, intent(inout) :: p
  integer :: iproc,ierr
 
  if (hasmpi) then
#ifndef EXPRO_NOMPI  
     call MPI_COMM_RANK(expro_comm,iproc,ierr)
     if (iproc == 0) read(1,*) p
     call MPI_BCAST(p,1,MPI_INTEGER,0,expro_comm,ierr)
#endif
  else
     read(1,*) p
  endif

end subroutine expro_icomm

subroutine expro_rcomm(x)

#ifndef EXPRO_NOMPI  
  use mpi
  use expro, only : expro_comm
#endif
  use expro, only : hasmpi

  implicit none

  double precision, intent(inout) :: x
  integer :: iproc,ierr
      
  if (hasmpi) then
#ifndef EXPRO_NOMPI  
     call MPI_COMM_RANK(expro_comm,iproc,ierr)
     if (iproc == 0) read(1,10) x
     call MPI_BCAST(x,1,MPI_DOUBLE_PRECISION,0,expro_comm,ierr)
#endif
  else
     read(1,10) x
  endif

10 format(1pe14.7)

end subroutine expro_rcomm

subroutine expro_scomm(x,n)

#ifndef EXPRO_NOMPI  
  use mpi
  use expro, only : expro_comm
#endif
  use expro, only : hasmpi

  implicit none

  integer, intent(in) :: n
  double precision, intent(inout), dimension(n) :: x
  integer :: iproc,ierr
    
  if (hasmpi) then
#ifndef EXPRO_NOMPI  
     call MPI_COMM_RANK(expro_comm,iproc,ierr)
     if (iproc == 0) read(1,*) x
     call MPI_BCAST(x,n,MPI_DOUBLE_PRECISION,0,expro_comm,ierr)
#endif
  else
     read(1,*) x
  endif
  
end subroutine expro_scomm

subroutine expro_lcomm(x,n)

#ifndef EXPRO_NOMPI  
  use mpi
  use expro, only : expro_comm
#endif
  use expro, only : hasmpi

  implicit none

  integer, intent(in) :: n
  double precision, intent(inout), dimension(n) :: x
  integer :: iproc,ierr
    
  if (hasmpi) then
#ifndef EXPRO_NOMPI  
     call MPI_COMM_RANK(expro_comm,iproc,ierr)
     if (iproc == 0) read(1,10) x
     call MPI_BCAST(x,n,MPI_DOUBLE_PRECISION,0,expro_comm,ierr)
#endif
  else
     read(1,10) x
  endif

10 format(10(1pe14.7))

end subroutine expro_lcomm

subroutine expro_tcomm(x,n)

#ifndef EXPRO_NOMPI  
  use mpi
  use expro, only : expro_comm
#endif
  use expro, only : expro_ion_max,hasmpi

  implicit none

  integer, intent(in) :: n
  character*10, intent(inout), dimension(expro_ion_max) :: x
  integer :: iproc,ierr
     
  if (hasmpi) then
#ifndef EXPRO_NOMPI  
     call MPI_COMM_RANK(expro_comm,iproc,ierr)
     if (iproc == 0) read(1,*) x(1:n)
     call MPI_BCAST(x(1:n),n*10,MPI_CHARACTER,0,expro_comm,ierr)
#endif
  else
     read(1,*) x(1:n)
  endif

end subroutine expro_tcomm

subroutine expro_vcomm(x,n)

#ifndef EXPRO_NOMPI
  use mpi
  use expro, only : expro_comm
#endif
  use expro, only : hasmpi

  implicit none

  integer :: idum,i
  integer, intent(in) :: n
  double precision, intent(inout), dimension(n) :: x
  integer :: iproc,ierr
    
  if (hasmpi) then
#ifndef EXPRO_NOMPI  
     call MPI_COMM_RANK(expro_comm,iproc,ierr)
     if (iproc == 0) then
        do i=1,n
           read(1,10) idum,x(i)
        enddo
     endif
     call MPI_BCAST(x,n,MPI_DOUBLE_PRECISION,0,expro_comm,ierr)
#endif
  else
     do i=1,n
        read(1,10) idum,x(i)
     enddo
  endif

10 format(i3,1x,1pe14.7)

end subroutine expro_vcomm

subroutine expro_acomm(x,m,n)

#ifndef EXPRO_NOMPI  
  use mpi
  use expro, only : expro_comm
#endif
  use expro, only : hasmpi

  implicit none

  integer :: idum,i
  integer, intent(in) :: m,n
  double precision, intent(inout), dimension(m,n) :: x
  integer :: iproc,ierr
  
  if (hasmpi) then
#ifndef EXPRO_NOMPI  
     call MPI_COMM_RANK(expro_comm,iproc,ierr)
     if (iproc == 0) then
        do i=1,n
           read(1,10) idum,x(:,i)
        enddo
     endif
     call MPI_BCAST(x,m*n,MPI_DOUBLE_PRECISION,0,expro_comm,ierr)
#endif
  else
     do i=1,n
        read(1,10) idum,x(:,i)
     enddo
  endif

10 format(i3,1x,10(1pe14.7,1x))

end subroutine expro_acomm

