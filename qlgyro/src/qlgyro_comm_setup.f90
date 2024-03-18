!-----------------------------------------------------------------
! qlgyro_comm_setup.f90
!
! PURPOSE:
!  Broadcast path and other information specific to each 
!  code instance.  Also, split QLGYRO_COMM_WORLD into 
!  required communicators:
!
!   - qlgyro_comm
!   - qlgyro_adj
!   - qlgyro_rad (for parallel block method)
!-----------------------------------------------------------------

subroutine qlgyro_comm_setup

  use mpi
  use qlgyro_globals

  implicit none

  integer :: i
  integer :: j
  integer :: i_color
  integer :: ip
  integer :: is

  integer :: low
  integer :: high
  integer :: splitkey

  integer, dimension(n_proc_global) :: colorvec
  integer, dimension(n_proc_global) :: workervec
  integer, dimension(n_proc_global) :: adjointvec
  integer, dimension(n_proc_global) :: workeradjvec


  ! Determine the number of cores for each ky 
  !-----------------------------
  ! Multi-job 
  !-----------------------------
  
  ! 1 worker; each DIR line specifies exact number of cores to GYRO
  
  !n_worker = 1

  ! Sort processors into communicators (color), workers
  ! and adjoint (to color).
  !  
  ! -   color: a specific communicator
  ! -  worker: color index at a given ky
  ! - adjoint: adjoint group to color 
  ! -   lproc: number of tasks per color
  ! -   lpath: local simulation path (fix this for GYRO) 
  
  low     = 0
  i_color = 0
  color = int(i_proc_global/procs)
  worker = 0
  adjoint = mod(i_proc_global, procs)
  workeradj = i_proc_global
  
  
  call MPI_GATHER(color,1,MPI_INTEGER,&
       colorvec,1,MPI_INTEGER,0,QLGYRO_COMM_WORLD,ierr)
  call MPI_GATHER(worker,1,MPI_INTEGER,&
       workervec,1,MPI_INTEGER,0,QLGYRO_COMM_WORLD,ierr)
  call MPI_GATHER(adjoint,1,MPI_INTEGER,&
       adjointvec,1,MPI_INTEGER,0,QLGYRO_COMM_WORLD,ierr)
  call MPI_GATHER(workeradj,1,MPI_INTEGER,&
       workeradjvec,1,MPI_INTEGER,0,QLGYRO_COMM_WORLD,ierr)

  if (i_proc_global == 0) then
     open(unit=1,file='out.qlgyro.taskmapping',status='replace')
     write(1,'(t2,a,t10,a,t18,a,t26,a,t34,a)') &
          'core','gcomm','worker','adjoint','workeradj'
     do i=1,n_proc_global
        write(1,'(5(i5,3x),a)') &
             i-1,colorvec(i),workervec(i),adjointvec(i),workeradjvec(i)
     enddo
     close(1)
  endif

  ! Split QLGYRO_COMM_WORLD into n_inst*n_worker different communicators.

  ! Choose key for task ordering
  splitkey = i_proc_global

  call MPI_COMM_SPLIT(QLGYRO_COMM_WORLD,&
       color,&
       splitkey,&
       qlgyro_comm,&
       ierr)
  if (ierr /= 0) then
     call qlgyro_catch_error('ERROR: (QLGYRO) GYRO_COMM not created') 
  endif

  call MPI_COMM_SPLIT(QLGYRO_COMM_WORLD,&
       adjoint,&
       splitkey,&
       qlgyro_adj,&
       ierr)
  if (ierr /= 0) then
     call qlgyro_catch_error('ERROR: (QLGYRO) GYRO_ADJ not created') 
  endif

  call MPI_COMM_SPLIT(QLGYRO_COMM_WORLD,&
       workeradj,&
       splitkey,&
       qlgyro_rad,&
       ierr)
  if (ierr /= 0) then
     call qlgyro_catch_error('ERROR: (QLGYRO) GYRO_RAD not created') 
  endif

  call MPI_COMM_RANK(qlgyro_comm,qlgyro_comm_rank,ierr)
  call MPI_COMM_RANK(qlgyro_adj,qlgyro_adj_rank,ierr)
  call MPI_COMM_RANK(qlgyro_rad,qlgyro_rad_rank,ierr)

  if (i_proc_global == 0) then
     open(unit=1,file=trim(runfile),position='append')
     write(1,'(t2,a)') 'INFO: (QLGYRO) MPI communicators split in QLGYRO'
     close(1)
  endif

end subroutine qlgyro_comm_setup
