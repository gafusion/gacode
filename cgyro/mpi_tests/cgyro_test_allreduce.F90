program cgyro_test_alltotall
  use mpi
  implicit none

  integer :: supported, i_err
  integer :: i_proc,n_proc
  integer :: i_proc_1,n_proc_1
  integer :: i_proc_2,n_proc_2
  integer :: i,n,j

  integer :: splitkey_1, i_group_1
  integer :: splitkey_2, i_group_2
  integer :: NEW_COMM_1, NEW_COMM_2

  ! allreduce buffers
  integer :: redbufsize
  real, dimension(:), allocatable :: redt1,redt2

  ! alltoall buffers
  integer :: total_bufsize
  integer :: bufsize
  real, dimension(:,:), allocatable :: t1,t2

  real:: m1,m2,m3

  call MPI_INIT_THREAD(MPI_THREAD_FUNNELED,supported,i_err)

  if (i_err /= 0) then
     write(*,*) 'ERROR: MPI_INIT_THREAD failed'
     call exit(1)
  endif

  call MPI_COMM_RANK(MPI_COMM_WORLD,i_proc,i_err)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,n_proc,i_err)

  if (n_proc .LT. (2*ALLTOALLNUM)) then
     write(*,*) i_proc, 'ERROR: n_proc < 2*ALLTOALLNUM', n_proc, 2*ALLTOALLNUM
     call exit(1)
  endif

  splitkey_1 =  i_proc
  splitkey_2 =  i_proc
  i_group_1 = modulo(i_proc,ALLTOALLNUM)
  i_group_2 = i_proc/ALLTOALLNUM
  call MPI_COMM_SPLIT(MPI_COMM_WORLD,&
       i_group_1,&
       splitkey_1,&
       NEW_COMM_1, &
       i_err)
  if (i_err /= 0) then
     write(*,*) i_proc, 'ERROR: MPI_COMM_SPLIT_1 failed'
     call exit(1)
  endif

  call MPI_COMM_RANK(NEW_COMM_1,i_proc_1,i_err)
  call MPI_COMM_SIZE(NEW_COMM_1,n_proc_1,i_err)

  if (n_proc_1 .LT. 2) then
     write(*,*) i_proc, 'ERROR: split n_proc_1 < 2', n_proc_1
     call exit(1)
  endif


  if (modulo(n_proc_1, 2) /= 0) then
     write(*,*) i_proc, 'ERROR: split n_proc_1 not even', n_proc_1
     call exit(1)
  endif

  call MPI_COMM_SPLIT(MPI_COMM_WORLD,&
       i_group_2,&
       splitkey_2,&
       NEW_COMM_2, &
       i_err)
  if (i_err /= 0) then
     write(*,*) i_proc, 'ERROR: MPI_COMM_SPLIT_2 failed'
     call exit(1)
  endif

  call MPI_COMM_RANK(NEW_COMM_2,i_proc_2,i_err)
  call MPI_COMM_SIZE(NEW_COMM_2,n_proc_2,i_err)

  if (n_proc_2 .LT. 2) then
     write(*,*) i_proc, 'ERROR: split n_proc_2 < 2', n_proc_2
     call exit(1)
  endif


  if (modulo(n_proc_2, 2) /= 0) then
     write(*,*) i_proc, 'ERROR: split n_proc_2 not even', n_proc_2
     call exit(1)
  endif

  if (n_proc_2 /= ALLTOALLNUM) then
     write(*,*) i_proc, 'ERROR: split n_proc_2 /= ALLTOALLNUM', n_proc_2, ALLTOALLNUM
     call exit(1)
  endif

  !-------------------------------
  ! Initialize buffers - AllReduce
  !-------------------------------

  redbufsize = 8*6*9*64
  redbufsize = redbufsize * ALLREDUCEBUFMUL

  if (i_proc==0) then
    write(*,*) "n_proc_1=",n_proc_1,"allreduce_bufsize=",redbufsize
    flush(6)
  endif

  allocate(redt1(redbufsize))
  allocate(redt2(redbufsize))

!$omp parallel do
  do i=1,redbufsize/2
    redt1(i) = 1.0*i_proc_1 + 0.5*i + 0.01*i_proc_2
    redt2(i) = 0.0
  enddo

#if defined(OMPGPU)
!$omp target enter data map(to:redt2)
!$omp target enter data map(to:redt1)
#elif defined(_OPENACC)
!$acc enter data copyin(redt2)
!$acc enter data copyin(redt1)
#endif

#if defined(OMPGPU)
!$omp target teams distribute parallel do
#elif defined(_OPENACC)
!$acc parallel loop gang vector independent present(redt1,redt2)
#else
!$omp parallel do
#endif
  do i=(1+redbufsize/2),redbufsize
    redt1(i) = 1.0*i_proc_1 + 0.5*i + 0.01*i_proc_2
    redt2(i) = 0.0
  enddo

#if defined(OMPGPU)
!$omp target update from(redt1)
#elif defined(_OPENACC)
!$acc update host(redt1)
#endif
  if (abs(redt1(1+redbufsize/2) - (1.0*i_proc_1 + 0.5*(1+redbufsize/2) + 0.01*i_proc_2)) .gt. 0.01) then
    write(*,*) i_proc, "WARNING: Initial redt1(redbufsize) validation failed!", &
               redt1(1+redbufsize/2), (1.0*i_proc_1 + 0.5*(1+redbufsize/2)+ 0.01*i_proc_2)
    call flush(6)
  endif

  !------------------------------
  ! Initialize buffers - AllToAll
  !------------------------------

  total_bufsize = 8*6*9*64
  total_bufsize = total_bufsize * ALLTOALLBUFMUL
  bufsize = total_bufsize/n_proc_2

  if (modulo(total_bufsize, n_proc_2*2) /=  0) then
    write(*,*) i_proc, "ERROR: total_bufsize not dividable by n_proc_2*2", total_bufsize, n_proc_2*2
    call exit(1)
  endif

  if (i_proc==0) then
    write(*,*) "n_proc_2=",n_proc_2,"alltoall_bufsize=",bufsize
    flush(6)
  endif

  allocate(t1(bufsize,n_proc_2))
  allocate(t2(bufsize,n_proc_2))

!$omp parallel do collapse(2)
  do n=1, n_proc_2
    do i=1,bufsize/2
      t1(i,n) = 100.0*i_proc_2 + i + 0.01*n
      t2(i,n) = 0.0
    enddo
  enddo

#if defined(OMPGPU)
!$omp target enter data map(to:t2)
!$omp target enter data map(to:t1)
#elif defined(_OPENACC)
!$acc enter data copyin(t2)
!$acc enter data copyin(t1)
#endif

#if defined(OMPGPU)
!$omp target teams distribute parallel do
#elif defined(_OPENACC)
!$acc parallel loop collapse(2) gang vector independent present(t1,t2)
#else
!$omp parallel do collapse(2)
#endif
do n=1, n_proc_2
    do i=(1+bufsize/2),bufsize
      t1(i,n) = 100.0*i_proc_2 + i + 0.01*n
      t2(i,n) = 0.0
    enddo
  enddo

#if defined(OMPGPU)
!$omp target update from(t1)
#elif defined(_OPENACC)
!$acc update host(t1)
#endif
  if (abs(t1(1+bufsize/2,n_proc_2/2) - (100.0*i_proc_2 + 1+bufsize/2 + 0.01*n_proc_2/2)) .gt. 0.01) then
    write(*,*) i_proc, "WARNING: Initial t1(bufsize,n_proc_2) validation failed!", &
               t1(1+bufsize/2,n_proc_2/2), (100.0*i_proc_2 + 1+bufsize/2+ 0.01*n_proc_2/2)
    call flush(6)
  endif

  m1 = MPI_Wtime()

#ifdef OMPGPU
!$omp target data use_device_ptr(redt1,redt2)
#else
!$acc host_data use_device(redt1,redt2)
#endif
  call MPI_ALLREDUCE(redt1,&
       redt2,&
       redbufsize,&
       MPI_DOUBLE_PRECISION,&
       MPI_SUM,&
       NEW_COMM_1,&
       i_err)
#ifdef OMPGPU
!$omp end target data
#else
!$acc end host_data
#endif

#if defined(OMPGPU)
!$omp target update from(redt2)
#elif defined(_OPENACC)
!$acc update host(redt2)
#endif
  if (abs(redt2(1+redbufsize/2) - (0.5*1.0*n_proc_1*(n_proc_1-1) + n_proc_1*(0.5*(1+redbufsize/2) + 0.01*i_proc_2))) .gt. 0.01) then
    write(*,*) i_proc, "WARNING: Initial redt2(redbufsize) validation failed!", &
               redt2(1+redbufsize/2), (0.5*1.0*n_proc_1*(n_proc_1-1) + n_proc_1*(0.5*(1+redbufsize/2) + 0.01*i_proc_2))
    call flush(6)
  endif

#ifdef OMPGPU
!$omp target data use_device_ptr(t1,t2)
#else
!$acc host_data use_device(t1,t2)
#endif
    call MPI_ALLTOALL(t1, &
         bufsize, &
         MPI_DOUBLE, &
         t2, &
         bufsize, &
         MPI_DOUBLE, &
         NEW_COMM_2, &
         i_err)
#ifdef OMPGPU
!$omp end target data
#else
!$acc end host_data
#endif

#if defined(OMPGPU)
!$omp target update from(t2)
#elif defined(_OPENACC)
!$acc update host(t2)
#endif
  if (abs(t2(1+bufsize/2,n_proc_2/2) - (100.0*(n_proc_2/2-1) + 1+bufsize/2 + 0.01*(i_proc_2+1))) .gt. 0.01) then
    write(*,*) i_proc, "WARNING: Initial t2(bufsize,n_proc_2) validation failed!", &
               t2(1+bufsize/2,n_proc_2/2), (100.0*(n_proc_2/2-1) + 1+bufsize/2+ 0.01*(i_proc_2+1))
    call flush(6)
  endif

  if (i_proc==0) then
     m1 = MPI_Wtime()
     write(*,*) m1, "Initial exchange succeeded"
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

#if defined(OMPGPU)
!$omp target teams distribute parallel do private(n)
#elif defined(_OPENACC)
!$acc parallel loop gang vector independent present(t2) private(n)
#else
!$omp parallel do private(n)
#endif
    do i=1,bufsize
!$acc loop seq
      do n=2, n_proc_2
        t2(i,n) = t2(i,n) +  t2(i,1)
      enddo
    enddo

#ifdef OMPGPU
!$omp target data use_device_ptr(t1,t2)
#else
!$acc host_data use_device(t1,t2)
#endif
      call MPI_ALLTOALL(t2, &
           bufsize, &
           MPI_DOUBLE, &
           t1, &
           bufsize, &
           MPI_DOUBLE, &
           NEW_COMM_2, &
           i_err)
#ifdef OMPGPU
!$omp end target data
#else
!$acc end host_data
#endif

#if defined(OMPGPU)
!$omp target teams distribute parallel do private(n)
#elif defined(_OPENACC)
!$acc parallel loop gang vector independent present(redt2,redt1) private(n)
#else
!$omp parallel do private(n)
#endif
  do i=1,redbufsize
    redt1(i) = redt1(i) +  redt2(((i-1)/2)+1)/128
  enddo

#ifdef OMPGPU
!$omp target data use_device_ptr(redt1,redt2)
#else
!$acc host_data use_device(redt1,redt2)
#endif
  call MPI_ALLREDUCE(redt1,&
       redt2,&
       redbufsize,&
       MPI_DOUBLE_PRECISION,&
       MPI_SUM,&
       NEW_COMM_1,&
       i_err)
#ifdef OMPGPU
!$omp end target data
#else
!$acc end host_data
#endif

#if defined(OMPGPU)
!$omp target teams distribute parallel do private(n)
#elif defined(_OPENACC)
!$acc parallel loop gang vector independent present(t1) private(n)
#else
!$omp parallel do private(n)
#endif
    do i=1,bufsize
!$acc loop seq
    do n=1, n_proc_2-1
        t1(i,n) = t1(i,n) + t1(i,n_proc_2)
      enddo
    enddo

#ifdef OMPGPU
!$omp target data use_device_ptr(t1,t2)
#else
!$acc host_data use_device(t1,t2)
#endif
      call MPI_ALLTOALL(t1, &
           bufsize, &
           MPI_DOUBLE, &
           t2, &
           bufsize, &
           MPI_DOUBLE, &
           NEW_COMM_2, &
           i_err)
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

#if defined(OMPGPU)
!$omp target exit data map(release:t1,redt1)
!$omp target exit data map(from:t2,redt2)
#elif defined(_OPENACC)
!$acc exit data delete(t1,redt1)
!$acc exit data copyout(t2,redt2)
#endif
  write(*,*) i_proc, "Success", redt2(i_proc_1), t2(1,i_proc_2+1), "dt=", m3-m2

end program

