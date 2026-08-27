program cgyro_test_alltotall
  use mpi
  implicit none

  integer :: supported, i_err
  integer :: i_proc,n_proc
  integer :: i,n,j
  integer :: cnt_err

  integer :: splitkey, i_group
  integer :: NEW_COMM_1
  integer :: req
  logical :: iflag
  integer :: istat(MPI_STATUS_SIZE)

  integer :: total_bufsize
  integer :: bufsize
  real, dimension(:,:), allocatable :: t1,t2,t3
  real :: c100,m100

  real:: m1,m2,m3

  c100 = 100.0
  m100 = 0.01

  call MPI_INIT_THREAD(MPI_THREAD_FUNNELED,supported,i_err)

  if (i_err /= 0) then
     write(*,*) 'ERROR: MPI_INIT_THREAD failed'
     call exit(1)
  endif

  call MPI_COMM_RANK(MPI_COMM_WORLD,i_proc,i_err)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,n_proc,i_err)

  if (n_proc .LT. 2) then
     write(*,*) i_proc, 'ERROR: n_proc < 2', n_proc
     call exit(1)
  endif

#ifdef DISABLE_MPI_SPLIT
  NEW_COMM_1 = MPI_COMM_WORLD
  if (i_proc==0) then
    write(*,*) "INFO: Using MPI_COMM_WORLD"
    flush(6)
  endif
#else
  splitkey =  i_proc
  i_group = 1
  call MPI_COMM_SPLIT(MPI_COMM_WORLD,&
       i_group,&
       splitkey,&
       NEW_COMM_1, &
       i_err)
  if (i_err /= 0) then
     write(*,*) i_proc, 'ERROR: MPI_COMM_SPLIT failed'
     call exit(1)
  endif

  call MPI_COMM_RANK(NEW_COMM_1,i_proc,i_err)
  call MPI_COMM_SIZE(NEW_COMM_1,n_proc,i_err)

  if (n_proc .LT. 2) then
     write(*,*) i_proc, 'ERROR: split n_proc < 2', n_proc
     call exit(1)
  endif
#endif


  if (modulo(n_proc, 2) /= 0) then
     write(*,*) i_proc, 'ERROR: split n_proc not even', n_proc
     call exit(1)
  endif

  !-------------------
  ! Initialize buffers
  !-------------------

  total_bufsize = 8*6*9*64
  total_bufsize = total_bufsize * ALLTOALLBUFMUL
  bufsize = total_bufsize/n_proc

  if (modulo(total_bufsize, n_proc*2) /=  0) then
    write(*,*) i_proc, "ERROR: total_bufsize not dividable by n_proc*2", total_bufsize, n_proc*2
    call exit(1)
  endif

  if (i_proc==0) then
    write(*,*) "n_proc=",n_proc,"bufsize=",bufsize
    flush(6)
  endif

  allocate(t1(bufsize,n_proc))
  allocate(t2(bufsize,n_proc))
  allocate(t3(bufsize,n_proc))

!$omp parallel do collapse(2)
  do n=1, n_proc
    do i=1,bufsize/2
      t1(i,n) = c100*i_proc + i + m100*n
      t2(i,n) = 0.0
      t3(i,n) = 0.0
    enddo
  enddo

#ifdef OMPGPU
!$omp target enter data map(to:t2)
!$omp target enter data map(to:t1)
!$omp target enter data map(to:t3)
#else
!$acc enter data copyin(t2)
!$acc enter data copyin(t1)
!$acc enter data copyin(t3)
#endif

#ifdef OMPGPU
!$omp target teams distribute parallel do collapse(2)
#else
#ifdef _OPENACC
!$acc parallel loop collapse(2) gang vector independent present(t1,t2,t3)
#else
!$omp parallel do collapse(2)
#endif
#endif
do n=1, n_proc
    do i=(1+bufsize/2),bufsize
      t1(i,n) = c100*i_proc + i + m100*n
      t2(i,n) = 0.0
      t3(i,n) = 0.0
    enddo
  enddo

#ifdef OMPGPU
!$omp target
#else
!$acc kernels present(t1)
#endif
  if (abs(t1(1+bufsize/2,n_proc/2) - (c100*i_proc + (1+bufsize/2) + m100*n_proc/2)) .gt. 0.00001) then
    write(*,*) i_proc, "WARNING: Initial t1(bufsize,n_proc) validation failed!", &
               t1(1+bufsize/2,n_proc/2), (c100*i_proc + (1+bufsize/2)+ m100*n_proc/2)
  endif
#ifdef OMPGPU
!$omp end target
#else
!$acc end kernels
#endif
  call flush(6)

  m1 = MPI_Wtime()

#ifdef OMPGPU
!$omp target data use_device_ptr(t1,t2)
#else
!$acc host_data use_device(t1,t2)
#endif

#ifdef NO_ASYNC_MPI
    call MPI_ALLTOALL(t1, &
         bufsize, &
         MPI_DOUBLE, &
         t2, &
         bufsize, &
         MPI_DOUBLE, &
         NEW_COMM_1, &
         i_err)
#else
    call MPI_IALLTOALL(t1, &
         bufsize, &
         MPI_DOUBLE, &
         t2, &
         bufsize, &
         MPI_DOUBLE, &
         NEW_COMM_1, &
         req, i_err)

    call MPI_REQUEST_GET_STATUS(req, iflag, istat, i_err)

    call MPI_WAIT(req, &
         istat, &
         i_err)
#endif

#ifdef OMPGPU
!$omp end target data
#else
!$acc end host_data
#endif

#ifdef OMPGPU
!$omp target
#else
!$acc kernels present(t2)
#endif
  if (abs(t2(1+bufsize/2,n_proc/2) - (c100*(n_proc/2-1) + (1+bufsize/2) + m100*(i_proc+1))) .gt. 0.00001) then
    write(*,*) i_proc, "WARNING: Initial t2(bufsize,n_proc) validation failed!", &
               t2(1+bufsize/2,n_proc/2), (c100*(n_proc/2-1) + (1+bufsize/2)+ m100*(i_proc+1))
  endif
#ifdef OMPGPU
!$omp end target
#else
!$acc end kernels
#endif

  call flush(6)

  if (i_proc==0) then
     m1 = MPI_Wtime()
     write(*,*) m1, "Initial exchange succeeded"
     call flush(6)
  endif

#ifdef OMPGPU
!$omp target data use_device_ptr(t2,t3)
#else
!$acc host_data use_device(t2,t3)
#endif

#ifdef NO_ASYNC_MPI
    call MPI_ALLTOALL(t2, &
         bufsize, &
         MPI_DOUBLE, &
         t3, &
         bufsize, &
         MPI_DOUBLE, &
         NEW_COMM_1, &
         i_err)
#else
    call MPI_IALLTOALL(t2, &
         bufsize, &
         MPI_DOUBLE, &
         t3, &
         bufsize, &
         MPI_DOUBLE, &
         NEW_COMM_1, &
         req, i_err)

    call MPI_REQUEST_GET_STATUS(req, iflag, istat, i_err)

    call MPI_WAIT(req, &
         istat, &
         i_err)
#endif

#ifdef OMPGPU
!$omp end target data
#else
!$acc end host_data
#endif

  cnt_err = 0
#ifdef OMPGPU
!$omp target teams distribute parallel do collapse(2) reduction(+:cnt_err)
#else
#ifdef _OPENACC
!$acc parallel loop collapse(2) gang vector independent present(t1,t3) copy(cnt_err) reduction(+:cnt_err)
#else
!$omp parallel do collapse(2) reduction(+:cnt_err)
#endif
#endif
do n=1, n_proc
    do i=1,bufsize
      if (t3(i,n)/=t1(i,n)) then
           write(*,*) i_proc, "WARNING: Initial t3 validation failed!",i,n,t3(i,n),t1(i,n)
           cnt_err = cnt_err + 1
       endif
    enddo
  enddo
  if (cnt_err>0) then
     write(*,*) i_proc, "WARNING: Initial t3 validation failed! CNT ",cnt_err
  endif

  call flush(6)

  if (i_proc==0) then
     m1 = MPI_Wtime()
     write(*,*) m1, "Initial backwards exchange succeeded"
     call flush(6)
  endif


  ! -------------------
  ! Now run many times for benchmark purposes
  m2 = MPI_Wtime()

  do j=1,1000
    if ((i_proc==0) .AND. (modulo(j,100) == 0)) then
      m1 = MPI_Wtime()
      write(*,*) m1, "j=", j, "dt=", m1-m2
      call flush(6)
    endif

#ifdef OMPGPU
!$omp target teams distribute parallel do private(n)
#else
#ifdef _OPENACC
!$acc parallel loop gang vector independent present(t2) private(n)
#else
!$omp parallel do private(n)
#endif
#endif
    do i=1,bufsize
!$acc loop seq
      do n=2, n_proc
        t2(i,n) = t2(i,n) +  t2(i,1)
      enddo
    enddo

#ifdef OMPGPU
!$omp target data use_device_ptr(t1,t2)
#else
!$acc host_data use_device(t1,t2)
#endif

#ifdef NO_ASYNC_MPI
      call MPI_ALLTOALL(t2, &
           bufsize, &
           MPI_DOUBLE, &
           t1, &
           bufsize, &
           MPI_DOUBLE, &
           NEW_COMM_1, &
           i_err)
#else
      call MPI_IALLTOALL(t2, &
           bufsize, &
           MPI_DOUBLE, &
           t1, &
           bufsize, &
           MPI_DOUBLE, &
           NEW_COMM_1, &
           req,i_err)

    call MPI_REQUEST_GET_STATUS(req, iflag, istat, i_err)

    call MPI_WAIT(req, &
         istat, &
         i_err)
#endif

#ifdef OMPGPU
!$omp end target data
#else
!$acc end host_data
#endif

#ifdef OMPGPU
!$omp target teams distribute parallel do private(n)
#else
#ifdef _OPENACC
!$acc parallel loop gang vector independent present(t1) private(n)
#else
!$omp parallel do private(n)
#endif
#endif
    do i=1,bufsize
!$acc loop seq
    do n=1, n_proc-1
        t1(i,n) = t1(i,n) + t1(i,n_proc)
      enddo
    enddo

#ifdef OMPGPU
!$omp target data use_device_ptr(t1,t2)
#else
!$acc host_data use_device(t1,t2)
#endif

#ifdef NO_ASYNC_MPI
      call MPI_ALLTOALL(t1, &
           bufsize, &
           MPI_DOUBLE, &
           t2, &
           bufsize, &
           MPI_DOUBLE, &
           NEW_COMM_1, &
           i_err)
#else
      call MPI_IALLTOALL(t1, &
           bufsize, &
           MPI_DOUBLE, &
           t2, &
           bufsize, &
           MPI_DOUBLE, &
           NEW_COMM_1, &
           req,i_err)

    call MPI_REQUEST_GET_STATUS(req, iflag, istat, i_err)

    call MPI_WAIT(req, &
         istat, &
         i_err)
#endif

#ifdef OMPGPU
!$omp end target data
#else
!$acc end host_data
#endif

  enddo
  m3 = MPI_Wtime()

  if (i_proc==0) then
     write(*,*) m3, "End of loop reached"
     call flush(6)
  endif

#ifdef OMPGPU
!$omp target exit data map(from:t2)
!$omp target exit data map(delete:t1)
#else
!$acc exit data delete(t1)
!$acc exit data copyout(t2)
#endif
  write(*,*) i_proc, "Success", t2(1,i_proc+1), "dt=", m3-m2

end program

