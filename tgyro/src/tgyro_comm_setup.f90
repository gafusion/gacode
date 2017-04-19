!-----------------------------------------------------------------
! tgyro_comm_setup.f90
!
! PURPOSE:
!  Broadcast path and other information specific to each 
!  code instance.  Also, split MPI_COMM_WORLD into 
!  required communicators:
!
!   - gyro_comm
!   - gyro_adj
!   - gyro_rad (for parallel block method)
!-----------------------------------------------------------------

subroutine tgyro_comm_setup

  use mpi
  use tgyro_globals

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


  ! Determine the number of "workers" at each radius

  select case (tgyro_mode)

  case (1)

     !-----------------------------
     ! Local Transport
     !-----------------------------

     ! Mange which profiles to evolve (1,2,3,4,5)=(ti,te,er,ne,he)

     ip = 0
     if (loc_ti_feedback_flag == 1) then
        ip = ip+1
        evolve_indx(ip) = 1 
     endif
     if (loc_te_feedback_flag == 1) then
        ip = ip+1
        evolve_indx(ip) = 2
     endif
     if (loc_er_feedback_flag == 1) then
        ip = ip+1
        evolve_indx(ip) = 3
     endif
     do is=0,loc_n_ion
        if (evo_e(is) == 1) then
           ip = ip+1
           evolve_indx(ip) = 4+is
        endif
     enddo
     n_evolve = ip

     if (tgyro_iteration_method == 5) then

        ! Parallel Jacobian 
        n_worker = n_evolve+1

        if (n_proc_global < n_worker*n_inst) then
           call tgyro_catch_error('ERROR: (TGYRO) Bad core count')
        endif

     else

        n_worker = 1

     endif

  case (3)

     !-----------------------------
     ! Multi-job 
     !-----------------------------

     ! 1 worker; each DIR line specifies exact number of cores to GYRO

     n_worker = 1

  end select

  ! Sort processors into communicators (color), workers
  ! and adjoint (to color).
  !  
  ! -   color: a specific communicator
  ! -  worker: color index at a given radius
  ! - adjoint: adjoint group to color 
  ! -   lproc: number of tasks per color
  ! -   lpath: local simulation path (fix this for GYRO) 
  ! -     i_r: local radius index (2, ... ,n_inst+1)

  low     = 0
  i_color = 0
  do i=1,n_inst
     do j=1,n_worker
        i_color = i_color+1
        high = low+procs(i)/n_worker-1
        if (i_proc_global >= low .and. i_proc_global <= high) then
           color   = i_color-1
           worker  = j-1
           adjoint = i_proc_global-low+worker*procs(i)/n_worker
           lproc   = procs(i)/n_worker
           lpath   = paths(i)
           lcode   = code(i)
           i_r     = i+1
           workeradj = i_proc_global-low+(i-1)*procs(i)/n_worker
        endif
        low = high+1
     enddo
  enddo

  call MPI_GATHER(color,1,MPI_INTEGER,&
       colorvec,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_GATHER(worker,1,MPI_INTEGER,&
       workervec,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_GATHER(adjoint,1,MPI_INTEGER,&
       adjointvec,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_GATHER(workeradj,1,MPI_INTEGER,&
       workeradjvec,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  if (i_proc_global == 0) then
     open(unit=1,file='out.tgyro.taskmapping',status='replace')
     write(1,'(t2,a,t10,a,t18,a,t26,a,t34,a)') &
          'core','gcomm','worker','adjoint','workeradj'
     do i=1,n_proc_global
        write(1,'(5(i5,3x),a)') &
             i-1,colorvec(i),workervec(i),adjointvec(i),workeradjvec(i)
     enddo
     close(1)
  endif

  ! Split MPI_COMM_WORLD into n_inst*n_worker different communicators.

  ! Choose key for task ordering
  splitkey = i_proc_global

  call MPI_COMM_SPLIT(MPI_COMM_WORLD,&
       color,&
       splitkey,&
       gyro_comm,&
       ierr)
  if (ierr /= 0) then
     call tgyro_catch_error('ERROR: (TGYRO) GYRO_COMM not created') 
  endif

  call MPI_COMM_SPLIT(MPI_COMM_WORLD,&
       adjoint,&
       splitkey,&
       gyro_adj,&
       ierr)
  if (ierr /= 0) then
     call tgyro_catch_error('ERROR: (TGYRO) GYRO_ADJ not created') 
  endif

  call MPI_COMM_SPLIT(MPI_COMM_WORLD,&
       workeradj,&
       splitkey,&
       gyro_rad,&
       ierr)
  if (ierr /= 0) then
     call tgyro_catch_error('ERROR: (TGYRO) GYRO_RAD not created') 
  endif

  call MPI_COMM_RANK(gyro_comm,gyro_comm_rank,ierr)
  call MPI_COMM_RANK(gyro_adj,gyro_adj_rank,ierr)
  call MPI_COMM_RANK(gyro_rad,gyro_rad_rank,ierr)

  if (tgyro_iteration_method == 5) then

     ! Set Jacobian indices for each worker (0,...,n_worker-1)

     if (worker == 0) then 
        worker_index = 0
        gyro_restart_method = 1
     else
        gyro_restart_method = 2
     endif

     do ip=1,n_evolve
        if (worker == ip) worker_index = evolve_indx(ip)
     enddo

  endif

  if (i_proc_global == 0) then
     open(unit=1,file=trim(runfile),position='append')
     write(1,'(t2,a)') 'INFO: (TGYRO) MPI communicators split in TGYRO'
     close(1)
  endif

end subroutine tgyro_comm_setup
