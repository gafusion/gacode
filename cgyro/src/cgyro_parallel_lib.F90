!-----------------------------------------------------------------
! cgyro_parallel_lib.F90
!
! PURPOSE:
!  Routines to manage the two types (lib,slib) of ALLTOALL
!  communication required by CGYRO.  This includes both GPU
!  and non-GPU options, as well as GPUDirect MPI options.
!-----------------------------------------------------------------

module parallel_lib

  implicit none

  ! flib
  ! simple Linear parallelization dimensions

  integer :: nfproc,ifproc
  integer, private :: nfi,nfj
  integer, private :: nfi_loc
  integer, private :: nfj_loc
  integer, private :: flib_comm

  ! lib
  ! aggregate Linear parallelization dimensions
  ! can be shared by many independent simulations

  integer :: nproc,iproc
  integer, private :: ni,nj
  integer, private :: ni_loc
  integer, private :: nj_loc
  integer, private :: nk_loc
  integer, private :: nk1,nk2
  integer, private :: nsm
  integer, private :: lib_comm
  integer, private :: nsend
  integer, private :: n_field
  integer, private :: nsend_real

  real, dimension(:,:,:,:,:), allocatable, private :: fsendr_real

  ! (expose these)
  complex, dimension(:,:,:,:), allocatable :: fsendf
  complex, dimension(:,:,:,:), allocatable :: fsendr

  ! clib
  ! subset of flib

  integer, private :: ncproc,icproc
  integer, private :: ns1,ns2
  integer, private :: clib_comm

  ! slib
  ! Nonlinear parallelization dimensions 

  integer :: nsproc,isproc
  integer, private :: nn
  integer, private :: nexch
  integer, private :: nkeep
  integer, private :: nsplit
  integer, private :: slib_comm

  integer, private, parameter :: default_size = 1

contains

  !=========================================================
  !  parallel_lib_f -> f(ni_loc,nj) -> g(nj_loc,ni) 
  !  parallel_lib_r -> g(nj_loc,ni) -> f(ni_loc,nj)
  !=========================================================

  subroutine parallel_lib_init(ni_in,nj_in,nj_loc_in,nk1_in,nk_loc_in,n_field_in,nsm_in,ni_loc_out,comm)

    use mpi

    implicit none

    integer, intent(in) :: ni_in,nj_in,nj_loc_in
    integer, intent(in) :: nk1_in,nk_loc_in,n_field_in,nsm_in
    integer, intent(in) :: comm
    integer, intent(inout) :: ni_loc_out
    integer, external :: parallel_dim
    integer :: ierr

    lib_comm = comm

    call MPI_COMM_RANK(lib_comm,iproc,ierr)
    call MPI_COMM_SIZE(lib_comm,nproc,ierr)

    ni = ni_in
    nj = nj_in

    ! parallel_dim(x,y) ~= x/y
    ni_loc = parallel_dim(ni,nproc)
    nj_loc = nj_loc_in
    nsm = nsm_in
    nk_loc = nk_loc_in

    ! nk1_in is typically iproc*nk_loc, but not always
    nk1 = nk1_in
    nk2 = nk1 + nk_loc -1

    ni_loc_out = ni_loc

    nsend = nj_loc*ni_loc*nk_loc

    allocate(fsendf(nj_loc,nk1:nk2,ni_loc,nproc))
    allocate(fsendr(ni_loc,nk1:nk2,nj_loc,nproc))

    n_field = n_field_in
    nsend_real = n_field*nsend
    if (.not. allocated(fsendr_real)) allocate(fsendr_real(n_field,ni_loc,nk1:nk2,nj_loc,nproc))

#if defined(OMPGPU)
!$omp target enter data map(alloc:fsendf,fsendr)
!$omp target enter data map(to:nproc,nk1,nk2,ni_loc,nsm)
#elif defined(_OPENACC)
!$acc enter data create(fsendf,fsendr)
!$acc enter data copyin(nproc,nk1,nk2,ni_loc,nsm)
#endif

  end subroutine parallel_lib_init

  !=========================================================

  subroutine parallel_lib_f_i_do(ft)

    use mpi

    implicit none

    complex, intent(inout), dimension(:,:,:) :: ft
    integer :: ierr

    call MPI_ALLTOALL(fsendf, &
         nsend, &
         MPI_DOUBLE_COMPLEX, &
         ft, &
         nsend, &
         MPI_DOUBLE_COMPLEX, &
         lib_comm, &
         ierr)

  end subroutine parallel_lib_f_i_do

  !=========================================================

  subroutine parallel_lib_f_i_do_gpu(ft)

    use mpi

    implicit none

    complex, intent(inout), dimension(:,:,:) :: ft
    integer :: ierr

#ifdef DISABLE_GPUDIRECT_MPI
#if defined(OMPGPU)
!$omp target update from(fsendf)
#elif defined(_OPENACC)
!$acc update host(fsendf)
#endif
#else
#if defined(OMPGPU)
!$omp target data use_device_addr(fsendf,ft)
#elif defined(_OPENACC)
!$acc host_data use_device(fsendf,ft)
#endif
#endif

    call MPI_ALLTOALL(fsendf, &
         nsend, &
         MPI_DOUBLE_COMPLEX, &
         ft, &
         nsend, &
         MPI_DOUBLE_COMPLEX, &
         lib_comm, &
         ierr)

#ifdef DISABLE_GPUDIRECT_MPI
#if defined(OMPGPU)
!$omp target update to(ft)
#elif defined(_OPENACC)
!$acc update device(ft)
#endif
#else
#if defined(OMPGPU)
!$omp end target data
#elif defined(_OPENACC)
!$acc end host_data
#endif
#endif

  end subroutine parallel_lib_f_i_do_gpu

  !=========================================================

  subroutine parallel_lib_r_do(f)

    use mpi

    implicit none

    complex, intent(inout), dimension(:,:,:,:) :: f
    integer :: ierr

    call MPI_ALLTOALL(fsendr, &
         nsend, &
         MPI_DOUBLE_COMPLEX, &
         f, &
         nsend, &
         MPI_DOUBLE_COMPLEX, &
         lib_comm, &
         ierr)

  end subroutine parallel_lib_r_do

  !=========================================================

  subroutine parallel_lib_r_do_gpu(f)

    use mpi

    implicit none

    complex, intent(inout), dimension(:,:,:,:) :: f

    integer :: ierr

#ifdef DISABLE_GPUDIRECT_MPI
#if defined(OMPGPU)
!$omp target update from(fsendr)
#elif defined(_OPENACC)
!$acc update host(fsendr)
#endif
#else
#if defined(OMPGPU)
!$omp target data use_device_addr(fsendr,f)
#elif defined(_OPENACC)
!$acc host_data use_device(fsendr,f)
#endif
#endif

    call MPI_ALLTOALL(fsendr, &
         nsend, &
         MPI_DOUBLE_COMPLEX, &
         f, &
         nsend, &
         MPI_DOUBLE_COMPLEX, &
         lib_comm, &
         ierr)

#ifdef DISABLE_GPUDIRECT_MPI
#if defined(OMPGPU)
!$omp target update to(f)
#elif defined(_OPENACC)
!$acc update device(f)
#endif
#else
#if defined(OMPGPU)
!$omp end target data
#elif defined(_OPENACC)
!$acc end host_data
#endif
#endif

  end subroutine parallel_lib_r_do_gpu

  !=========================================================

  subroutine parallel_lib_rtrans_pack(fin)
! -----------------------------------------
! transpose version of parallel_lib_r_pack(fin)
! -----------------------------------------
    use mpi

    implicit none

    complex, intent(in), dimension(:,:,:) :: fin
    integer :: j_loc,i,j,k,j1,j2,itor

    ! fin is assumed to be cap_h_c(nc,nv_loc,nt_loc)
    ! where nc = ni

    j1 = 1+iproc*nj_loc
    j2 = (1+iproc)*nj_loc

!$omp parallel do collapse(3) if (size(fsendr) >= default_size) default(none) &
!$omp& shared(nproc,j1,j2,ni_loc,nk1,nk2,nsm) &
!$omp& private(j,j_loc,i) &
!$omp& shared(fin,fsendr)
    do k=1,nproc
      do j=j1,j2
       do itor=nk1,nk2
          do i=1,ni_loc
             j_loc = j-j1+1 
             fsendr(i,itor,j_loc,k) = fin(i+(k-1)*ni_loc, j_loc, 1+(itor-nk1))
          enddo
       enddo
      enddo
    enddo

  end subroutine parallel_lib_rtrans_pack

  !=========================================================

  subroutine parallel_lib_rtrans_pack_gpu(fin)
! -----------------------------------------
! transpose version of parallel_lib_r_pack(fin)
! -----------------------------------------

    use mpi

    implicit none

    complex, intent(in), dimension(:,:,:) :: fin
    integer :: j_loc,i,j,k,j1,j2,itor

    ! fin is assumed to be cap_h_c(nc,nv_loc,nt_loc)
    ! where nc = ni

    j1 = 1+iproc*nj_loc
    j2 = (1+iproc)*nj_loc

#if defined(OMPGPU)
!$omp target teams distribute parallel do simd collapse(4) &
!$omp&  private(j_loc) map(to:j1,j2)
#else
!$acc parallel loop collapse(4) gang vector independent private(j_loc) &
!$acc&         present(fsendr,fin) present(nproc,nk1,nk2,ni_loc,nsm) &
!$acc&         copyin(j1,j2) default(none)
#endif
    do k=1,nproc
      do j=j1,j2
       do itor=nk1,nk2
          do i=1,ni_loc
             j_loc = j-j1+1
             fsendr(i,itor,j_loc,k) = fin(i+(k-1)*ni_loc,j_loc,1+(itor-nk1))
          enddo
       enddo
      enddo
    enddo

  end subroutine parallel_lib_rtrans_pack_gpu

  !=========================================================

  subroutine parallel_lib_rtrans(fin,f)
! -----------------------------------------
! transpose version of parallel_lib_r(fin,f)
! -----------------------------------------
    use mpi

    implicit none

    complex, intent(in), dimension(:,:,:) :: fin
    complex, intent(inout), dimension(:,:,:,:) :: f

    call parallel_lib_rtrans_pack(fin)
    call parallel_lib_r_do(f)

  end subroutine parallel_lib_rtrans

  !=========================================================

  subroutine parallel_lib_rtrans_real(fin,f)

    use mpi

    implicit none

    real, intent(in), dimension(:,:,:,:) :: fin
    real, intent(inout), dimension(:,:,:,:,:) :: f
    integer :: ierr,j_loc,i,j,k,j1,j2,itor,fi

    j1 = 1+iproc*nj_loc
    j2 = (1+iproc)*nj_loc

!$omp parallel do collapse(2) if (size(fsendr_real) >= default_size) default(none) &
!$omp& firstprivate(nproc,j1,j2,ni_loc,nk1,nk2,n_field) &
!$omp& private(j,j_loc,i,fi) &
!$omp& shared(fin,fsendr_real)
    do k=1,nproc
     do itor=nk1,nk2
       do j=j1,j2
          j_loc = j-j1+1
          do i=1,ni_loc
             do fi=1,n_field
                fsendr_real(fi,i,itor,j_loc,k) = fin(fi,i+(k-1)*ni_loc,j_loc,1+(itor-nk1)) 
             enddo
          enddo
       enddo
     enddo
    enddo

    call MPI_ALLTOALL(fsendr_real, &
         nsend_real, &
         MPI_DOUBLE_PRECISION,&
         f, &
         nsend_real, &
         MPI_DOUBLE_PRECISION, &
         lib_comm, &
         ierr)

  end subroutine parallel_lib_rtrans_real

  !=========================================================

  subroutine parallel_lib_collect_field(field_loc_v,field_v)

    use mpi

    implicit none

    complex, intent(in), dimension(:,:,:,:) :: field_loc_v
    complex, intent(inout), dimension(:,:,:,:) :: field_v
    integer :: ierr

#ifdef DISABLE_GPUDIRECT_MPI
#if defined(OMPGPU)
!$omp target update from(field_loc_v)
#elif defined(_OPENACC)
!$acc update host(field_loc_v)
#endif
#else
#if defined(OMPGPU)
!$omp target data use_device_addr(field_loc_v,field_v)
#elif defined(_OPENACC)
!$acc host_data use_device(field_loc_v,field_v)
#endif
#endif

    call MPI_ALLGATHER(field_loc_v(:,:,:,:),&
         size(field_loc_v(:,:,:,:)),&
         MPI_DOUBLE_COMPLEX,&
         field_v(:,:,:,:),&
         size(field_loc_v(:,:,:,:)),&
         MPI_DOUBLE_COMPLEX,&
         lib_comm,&
         ierr)

#ifdef DISABLE_GPUDIRECT_MPI
#if defined(OMPGPU)
!$omp target update to(field_v)
#elif defined(_OPENACC)
!$acc update device(field_v)
#endif
#else
#if defined(OMPGPU)
!$omp end target data
#elif defined(_OPENACC)
!$acc end host_data
#endif
#endif

  end subroutine parallel_lib_collect_field

  !=========================================================

  subroutine parallel_flib_init(ni_in,nj_in,ni_loc_out,nj_loc_out,comm)

    use mpi

    implicit none

    integer, intent(in) :: ni_in,nj_in
    integer, intent(in) :: comm
    integer, intent(inout) :: ni_loc_out,nj_loc_out
    integer, external :: parallel_dim
    integer :: ierr

    flib_comm = comm

    call MPI_COMM_RANK(flib_comm,ifproc,ierr)
    call MPI_COMM_SIZE(flib_comm,nfproc,ierr)

    nfi = ni_in
    nfj = nj_in

    ! parallel_dim(x,y) ~= x/y
    nfi_loc = parallel_dim(nfi,nfproc)
    nfj_loc = parallel_dim(nfj,nfproc)

    ni_loc_out = nfi_loc
    nj_loc_out = nfj_loc

  end subroutine parallel_flib_init

  !=========================================================

  ! Note: Using intent(inout) for field_loc due to possible copy from GPU to CPU memory
  ! Same for many other argumnets in other subroutines

  subroutine parallel_flib_sum_field(field_loc,field)

    use mpi

    implicit none

    complex, intent(inout), dimension(:,:,:) :: field_loc
    complex, intent(inout), dimension(:,:,:) :: field
    integer :: ierr

#ifdef DISABLE_GPUDIRECT_MPI
#if defined(OMPGPU)
!$omp target update from(field_loc)
#elif defined(_OPENACC)
!$acc update host(field_loc)
#endif
#else
#if defined(OMPGPU)
!$omp target data use_device_addr(field_loc,field)
#elif defined(_OPENACC)
!$acc host_data use_device(field_loc,field)
#endif
#endif

    call MPI_ALLREDUCE(field_loc(:,:,:),&
          field(:,:,:),&
          size(field(:,:,:)),&
          MPI_DOUBLE_COMPLEX,&
          MPI_SUM,&
          flib_comm,&
          ierr)

#ifdef DISABLE_GPUDIRECT_MPI
#if defined(OMPGPU)
!$omp target update to(field)
#elif defined(_OPENACC)
!$acc update device(field)
#endif
#else
#if defined(OMPGPU)
!$omp end target data
#elif defined(_OPENACC)
!$acc end host_data
#endif
#endif

  end subroutine parallel_flib_sum_field

  !=========================================================
  !  Species communicator
  !=========================================================

  subroutine parallel_clib_init(ns1_in,ns2_in,ns_loc_out,comm)

    use mpi

    implicit none

    integer, intent(in) :: ns1_in,ns2_in
    integer, intent(inout) :: ns_loc_out
    integer, intent(in) :: comm
    integer :: ierr

    clib_comm = comm

    call MPI_COMM_RANK(clib_comm,icproc,ierr)
    call MPI_COMM_SIZE(clib_comm,ncproc,ierr)

    ns1 = ns1_in
    ns2 = ns2_in

    ns_loc_out = ns2-ns1+1

  end subroutine parallel_clib_init

  !=========================================================

  subroutine parallel_clib_sum_upwind(upwind_loc,upwind)

    use mpi

    implicit none

    complex, intent(inout), dimension(:,:,:) :: upwind_loc
    complex, intent(inout), dimension(:,:,:) :: upwind
    integer :: ierr

#ifdef DISABLE_GPUDIRECT_MPI
#if defined(OMPGPU)
!$omp target update from(upwind_loc)
#elif defined(_OPENACC)
!$acc update host(upwind_loc)
#endif
#else
#if defined(OMPGPU)
!$omp target data use_device_addr(upwind_loc,upwind)
#elif defined(_OPENACC)
!$acc host_data use_device(upwind_loc,upwind)
#endif
#endif

    call MPI_ALLREDUCE(upwind_loc(:,:,:),&
          upwind(:,:,:),&
          size(upwind(:,:,:)),&
          MPI_DOUBLE_COMPLEX,&
          MPI_SUM,&
          clib_comm,&
          ierr)

#ifdef DISABLE_GPUDIRECT_MPI
#if defined(OMPGPU)
!$omp target update to(upwind)
#elif defined(_OPENACC)
!$acc update device(upwind)
#endif
#else
#if defined(OMPGPU)
!$omp end target data
#elif defined(_OPENACC)
!$acc end host_data
#endif
#endif

  end subroutine parallel_clib_sum_upwind

  !=========================================================
  subroutine parallel_clib_sum_upwind32(upwind_loc,upwind)

    use mpi
    use, intrinsic :: iso_fortran_env

    implicit none

    complex(KIND=REAL32), intent(inout), dimension(:,:,:) :: upwind_loc
    complex(KIND=REAL32), intent(inout), dimension(:,:,:) :: upwind
    integer :: ierr

#ifdef DISABLE_GPUDIRECT_MPI
#if defined(OMPGPU)
!$omp target update from(upwind_loc)
#elif defined(_OPENACC)
!$acc update host(upwind_loc)
#endif
#else
#if defined(OMPGPU)
!$omp target data use_device_addr(upwind_loc,upwind)
#elif defined(_OPENACC)
!$acc host_data use_device(upwind_loc,upwind)
#endif
#endif

    call MPI_ALLREDUCE(upwind_loc(:,:,:),&
          upwind(:,:,:),&
          size(upwind(:,:,:)),&
          MPI_COMPLEX,&
          MPI_SUM,&
          clib_comm,&
          ierr)

#ifdef DISABLE_GPUDIRECT_MPI
#if defined(OMPGPU)
!$omp target update to(upwind)
#elif defined(_OPENACC)
!$acc update device(upwind)
#endif
#else
#if defined(OMPGPU)
!$omp end target data
#elif defined(_OPENACC)
!$acc end host_data
#endif
#endif

  end subroutine parallel_clib_sum_upwind32

  !=========================================================

  !  parallel_slib_f -> f(ni_loc,nj) -> g(nj_loc,ni) 
  !  parallel_slib_r -> g(nj_loc,ni) -> f(ni_loc,nj)

  subroutine parallel_slib_init(nn_in,nexch_in,nkeep_in,nsplit_out,comm)

    use mpi

    !-------------------------------------------
    implicit none
    !
    integer, intent(in) :: nn_in
    integer, intent(in) :: nexch_in
    integer, intent(in) :: nkeep_in
    integer, intent(inout) :: nsplit_out
    integer, intent(in) :: comm
    integer :: ierr
    !-------------------------------------------

    slib_comm = comm

    !-------------------------------------------------
    ! Get rank and number of processors from slib_comm 
    ! (which is assumed to already exist).
    !
    call MPI_COMM_RANK(slib_comm,isproc,ierr)
    call MPI_COMM_SIZE(slib_comm,nsproc,ierr)
    !-----------------------------------------------

    ! nn_in == nsproc
    nn  = nn_in
    nexch = nexch_in
    nkeep = nkeep_in

    nsplit = 1+(nexch-1)/nn

    nsplit_out = nsplit

  end subroutine parallel_slib_init

  ! test an async req, to progress async operations
  subroutine parallel_slib_test(req)
    use mpi
    !-------------------------------------------------------
    implicit none
    !
    integer, intent(inout) :: req
    !
    logical :: iflag
    integer :: ierr
    integer :: istat(MPI_STATUS_SIZE)
    !-------------------------------------------------------

#ifndef NO_ASYNC_MPI

    call MPI_REQUEST_GET_STATUS(req, iflag, istat, ierr)
    ! we discard all the outputs... it is just a way to progress async mpi

#endif
    !else, noop

  end subroutine parallel_slib_test

!=========================================================

  subroutine parallel_slib_f_nc(nsx,x,xt)

    use mpi

    !-------------------------------------------------------
    implicit none
    !
    integer, intent(in) :: nsx
    complex, intent(inout), dimension(nkeep,nk_loc,nsx*nn) :: x
    complex, intent(inout), dimension(nkeep,nk_loc,nsx,nn) :: xt
    !
    integer :: ierr
    !-------------------------------------------------------

#ifdef DISABLE_GPUDIRECT_MPI
#if defined(OMPGPU)
!$omp target update from(x)
#elif defined(_OPENACC)
!$acc update host(x)
#endif
#else
#if defined(OMPGPU)
!$omp target data use_device_addr(x,xt)
#elif defined(_OPENACC)
!$acc host_data use_device(x,xt)
#endif
#endif

    call MPI_ALLTOALL(x, &
         nkeep*nk_loc*nsx, &
         MPI_DOUBLE_COMPLEX, &
         xt, &
         nkeep*nk_loc*nsx, &
         MPI_DOUBLE_COMPLEX, &
         slib_comm, &
         ierr)

#ifdef DISABLE_GPUDIRECT_MPI
#if defined(OMPGPU)
!$omp target update to(xt)
#elif defined(_OPENACC)
!$acc update device(xt)
#endif
#else
#if defined(OMPGPU)
!$omp end target data
#elif defined(_OPENACC)
!$acc end host_data
#endif
#endif

  end subroutine parallel_slib_f_nc

  subroutine parallel_slib_f_nc_async(nsx,x,xt,req)
    use mpi
    !-------------------------------------------------------
    implicit none
    !
    integer, intent(in) :: nsx
    complex, intent(inout), dimension(nkeep,nk_loc,nsx*nn) :: x
    complex, intent(inout), dimension(nkeep,nk_loc,nsx,nn) :: xt
    integer, intent(inout) :: req
    !
    integer :: ierr
    !-------------------------------------------------------

#ifdef NO_ASYNC_MPI
   call parallel_slib_f_nc(nsx,x,xt)

#else

#ifdef DISABLE_GPUDIRECT_MPI
#if defined(OMPGPU)
!$omp target update from(x)
#elif defined(_OPENACC)
!$acc update host(x)
#endif
#else
#if defined(OMPGPU)
!$omp target data use_device_addr(x,xt)
#elif defined(_OPENACC)
!$acc host_data use_device(x,xt)
#endif
#endif

   call MPI_IALLTOALL(x, &
         nkeep*nk_loc*nsx, &
         MPI_DOUBLE_COMPLEX, &
         xt, &
         nkeep*nk_loc*nsx, &
         MPI_DOUBLE_COMPLEX, &
         slib_comm, &
         req, &
         ierr)

#ifndef DISABLE_GPUDIRECT_MPI
#if defined(OMPGPU)
!$omp end target data
#elif defined(_OPENACC)
!$acc end host_data
#endif
#endif

#endif

  end subroutine parallel_slib_f_nc_async

  subroutine parallel_slib_f_nc32(nsx,x,xt)

    use mpi
    use, intrinsic :: iso_fortran_env

    !-------------------------------------------------------
    implicit none
    !
    integer, intent(in) :: nsx
    complex(KIND=REAL32), intent(inout), dimension(nkeep,nk_loc,nsx*nn) :: x
    complex(KIND=REAL32), intent(inout), dimension(nkeep,nk_loc,nsx,nn) :: xt
    !
    integer :: ierr
    !-------------------------------------------------------

#ifdef DISABLE_GPUDIRECT_MPI
#if defined(OMPGPU)
!$omp target update from(x)
#elif defined(_OPENACC)
!$acc update host(x)
#endif
#else
#if defined(OMPGPU)
!$omp target data use_device_addr(x,xt)
#elif defined(_OPENACC)
!$acc host_data use_device(x,xt)
#endif
#endif

    call MPI_ALLTOALL(x, &
         nkeep*nk_loc*nsx, &
         MPI_COMPLEX, &
         xt, &
         nkeep*nk_loc*nsx, &
         MPI_COMPLEX, &
         slib_comm, &
         ierr)

#ifdef DISABLE_GPUDIRECT_MPI
#if defined(OMPGPU)
!$omp target update to(xt)
#elif defined(_OPENACC)
!$acc update device(xt)
#endif
#else
#if defined(OMPGPU)
!$omp end target data
#elif defined(_OPENACC)
!$acc end host_data
#endif
#endif

  end subroutine parallel_slib_f_nc32

  subroutine parallel_slib_f_nc32_async(nsx,x,xt,req)
    use mpi
    use, intrinsic :: iso_fortran_env
    !-------------------------------------------------------
    implicit none
    !
    integer, intent(in) :: nsx
    complex(KIND=REAL32), intent(inout), dimension(nkeep,nk_loc,nsx*nn) :: x
    complex(KIND=REAL32), intent(inout), dimension(nkeep,nk_loc,nsx,nn) :: xt
    integer, intent(inout) :: req
    !
    integer :: ierr
    !-------------------------------------------------------

#ifdef NO_ASYNC_MPI
   call parallel_slib_f_nc32(nsx,x,xt)

#else

#ifdef DISABLE_GPUDIRECT_MPI
#if defined(OMPGPU)
!$omp target update from(x)
#elif defined(_OPENACC)
!$acc update host(x)
#endif
#else
#if defined(OMPGPU)
!$omp target data use_device_addr(x,xt)
#elif defined(_OPENACC)
!$acc host_data use_device(x,xt)
#endif
#endif

   call MPI_IALLTOALL(x, &
         nkeep*nk_loc*nsx, &
         MPI_COMPLEX, &
         xt, &
         nkeep*nk_loc*nsx, &
         MPI_COMPLEX, &
         slib_comm, &
         req, &
         ierr)

#ifndef DISABLE_GPUDIRECT_MPI
#if defined(OMPGPU)
!$omp end target data
#elif defined(_OPENACC)
!$acc end host_data
#endif
#endif

#endif

  end subroutine parallel_slib_f_nc32_async

  ! require x and xt to ensure they exist until this finishes
  subroutine parallel_slib_f_nc_wait(nsx,x,xt,req)
    use mpi
    !-------------------------------------------------------
    implicit none
    !
    integer, intent(in) :: nsx
    complex, intent(inout), dimension(nkeep,nk_loc,nsx*nn) :: x
    complex, intent(inout), dimension(nkeep,nk_loc,nsx,nn) :: xt
    integer, intent(inout) :: req
    !
    integer :: ierr
    integer :: istat(MPI_STATUS_SIZE)
    !-------------------------------------------------------

#ifndef NO_ASYNC_MPI

    call MPI_WAIT(req, &
         istat, &
         ierr)

#ifdef DISABLE_GPUDIRECT_MPI
#if defined(OMPGPU)
!$omp target update to(xt)
#elif defined(_OPENACC)
!$acc update device(xt)
#endif
#endif

#endif
   !else, noop

  end subroutine parallel_slib_f_nc_wait

  ! require x and xt to ensure they exist until this finishes
  subroutine parallel_slib_f_nc32_wait(nsx,x,xt,req)
    use mpi
    use, intrinsic :: iso_fortran_env
    !-------------------------------------------------------
    implicit none
    !
    integer, intent(in) :: nsx
    complex(KIND=REAL32), intent(inout), dimension(nkeep,nk_loc,nsx*nn) :: x
    complex(KIND=REAL32), intent(inout), dimension(nkeep,nk_loc,nsx,nn) :: xt
    integer, intent(inout) :: req
    !
    integer :: ierr
    integer :: istat(MPI_STATUS_SIZE)
    !-------------------------------------------------------

#ifndef NO_ASYNC_MPI

    call MPI_WAIT(req, &
         istat, &
         ierr)

#ifdef DISABLE_GPUDIRECT_MPI
#if defined(OMPGPU)
!$omp target update to(xt)
#elif defined(_OPENACC)
!$acc update device(xt)
#endif
#endif

#endif
   !else, noop

  end subroutine parallel_slib_f_nc32_wait

 !=========================================================

  subroutine parallel_slib_r_nc (nsx,xt,x)
    use mpi
    !-------------------------------------------------------
    implicit none
    !
    integer, intent(in) :: nsx
    complex, intent(inout), dimension(nkeep,nk_loc,nsx,nn) :: xt
    complex, intent(inout), dimension(nkeep,nk_loc,nsx*nn) :: x
    !
    integer :: ierr
    !-------------------------------------------------------

#ifdef DISABLE_GPUDIRECT_MPI
#if defined(OMPGPU)
!$omp target update from(xt)
#elif defined(_OPENACC)
!$acc update host(xt)
#endif
#else
#if defined(OMPGPU)
!$omp target data use_device_addr(xt,x)
#elif defined(_OPENACC)
!$acc host_data use_device(xt,x)
#endif
#endif

    call MPI_ALLTOALL(xt, &
         nkeep*nk_loc*nsx, &
         MPI_DOUBLE_COMPLEX, &
         x, &
         nkeep*nk_loc*nsx, &
         MPI_DOUBLE_COMPLEX, &
         slib_comm, &
         ierr)

#ifdef DISABLE_GPUDIRECT_MPI
#if defined(OMPGPU)
!$omp target update to(x)
#elif defined(_OPENACC)
!$acc update device(x)
#endif
#else
#if defined(OMPGPU)
!$omp end target data
#elif defined(_OPENACC)
!$acc end host_data
#endif
#endif

  end subroutine parallel_slib_r_nc

  subroutine parallel_slib_r_nc32 (nsx,xt,x)
    use mpi
    use, intrinsic :: iso_fortran_env
    !-------------------------------------------------------
    implicit none
    !
    integer, intent(in) :: nsx
    complex(KIND=REAL32), intent(inout), dimension(nkeep,nk_loc,nsx,nn) :: xt
    complex(KIND=REAL32), intent(inout), dimension(nkeep,nk_loc,nsx*nn) :: x
    !
    integer :: ierr
    !-------------------------------------------------------

#ifdef DISABLE_GPUDIRECT_MPI
#if defined(OMPGPU)
!$omp target update from(xt)
#elif defined(_OPENACC)
!$acc update host(xt)
#endif
#else
#if defined(OMPGPU)
!$omp target data use_device_addr(xt,x)
#elif defined(_OPENACC)
!$acc host_data use_device(xt,x)
#endif
#endif

    call MPI_ALLTOALL(xt, &
         nkeep*nk_loc*nsx, &
         MPI_COMPLEX, &
         x, &
         nkeep*nk_loc*nsx, &
         MPI_COMPLEX, &
         slib_comm, &
         ierr)

#ifdef DISABLE_GPUDIRECT_MPI
#if defined(OMPGPU)
!$omp target update to(x)
#elif defined(_OPENACC)
!$acc update device(x)
#endif
#else
#if defined(OMPGPU)
!$omp end target data
#elif defined(_OPENACC)
!$acc end host_data
#endif
#endif

  end subroutine parallel_slib_r_nc32

  subroutine parallel_slib_r_nc_async (nsx,xt,x,req)
    use mpi
    !-------------------------------------------------------
    implicit none
    !
    integer, intent(in) :: nsx
    complex, intent(inout), dimension(nkeep,nk_loc,nsx,nn) :: xt
    complex, intent(inout), dimension(nkeep,nk_loc,nsx*nn) :: x
    integer, intent(inout) :: req
    !
    integer :: ierr
    !-------------------------------------------------------

#ifdef NO_ASYNC_MPI
   call parallel_slib_r_nc(nsx,xt,x)

#else

#ifdef DISABLE_GPUDIRECT_MPI
#if defined(OMPGPU)
!$omp target update from(xt)
#elif defined(_OPENACC)
!$acc update host(xt)
#endif
#else
#if defined(OMPGPU)
!$omp target data use_device_addr(xt,x)
#elif defined(_OPENACC)
!$acc host_data use_device(xt,x)
#endif
#endif

    call MPI_IALLTOALL(xt, &
         nkeep*nk_loc*nsx, &
         MPI_DOUBLE_COMPLEX, &
         x, &
         nkeep*nk_loc*nsx, &
         MPI_DOUBLE_COMPLEX, &
         slib_comm, &
         req, &
         ierr)

#ifndef DISABLE_GPUDIRECT_MPI
#if defined(OMPGPU)
!$omp end target data
#elif defined(_OPENACC)
!$acc end host_data
#endif
#endif
#endif

  end subroutine parallel_slib_r_nc_async

  subroutine parallel_slib_r_nc32_async (nsx,xt,x,req)
    use mpi
    use, intrinsic :: iso_fortran_env
    !-------------------------------------------------------
    implicit none
    !
    integer, intent(in) :: nsx
    complex(KIND=REAL32), intent(inout), dimension(nkeep,nk_loc,nsx,nn) :: xt
    complex(KIND=REAL32), intent(inout), dimension(nkeep,nk_loc,nsx*nn) :: x
    integer, intent(inout) :: req
    !
    integer :: ierr
    !-------------------------------------------------------

#ifdef NO_ASYNC_MPI
   call parallel_slib_r_nc32(nsx,xt,x)

#else

#ifdef DISABLE_GPUDIRECT_MPI
#if defined(OMPGPU)
!$omp target update from(xt)
#elif defined(_OPENACC)
!$acc update host(xt)
#endif
#else
#if defined(OMPGPU)
!$omp target data use_device_addr(xt,x)
#elif defined(_OPENACC)
!$acc host_data use_device(xt,x)
#endif
#endif

    call MPI_IALLTOALL(xt, &
         nkeep*nk_loc*nsx, &
         MPI_COMPLEX, &
         x, &
         nkeep*nk_loc*nsx, &
         MPI_COMPLEX, &
         slib_comm, &
         req, &
         ierr)

#ifndef DISABLE_GPUDIRECT_MPI
#if defined(OMPGPU)
!$omp end target data
#elif defined(_OPENACC)
!$acc end host_data
#endif
#endif
#endif

  end subroutine parallel_slib_r_nc32_async

  ! require x and xt to ensure they exist until this finishes
  subroutine parallel_slib_r_nc_wait(nsx,xt,x,req)
    use mpi
    !-------------------------------------------------------
    implicit none
    !
    integer, intent(in) :: nsx
    complex, intent(inout), dimension(nkeep,nk_loc,nsx,nn) :: xt
    complex, intent(inout), dimension(nkeep,nk_loc,nsx*nn) :: x
    integer, intent(inout) :: req
    !
    integer :: ierr
    integer :: istat(MPI_STATUS_SIZE)
    !-------------------------------------------------------

#ifndef NO_ASYNC_MPI

    call MPI_WAIT(req, &
         istat, &
         ierr)

#ifdef DISABLE_GPUDIRECT_MPI
#if defined(OMPGPU)
!$omp target update to(x)
#elif defined(_OPENACC)
!$acc update device(x)
#endif
#endif

#endif
   !else, noop

  end subroutine parallel_slib_r_nc_wait

  subroutine parallel_slib_r_nc32_wait(nsx,xt,x,req)
    use mpi
    use, intrinsic :: iso_fortran_env
    !-------------------------------------------------------
    implicit none
    !
    integer, intent(in) :: nsx
    complex(KIND=REAL32), intent(inout), dimension(nkeep,nk_loc,nsx,nn) :: xt
    complex(KIND=REAL32), intent(inout), dimension(nkeep,nk_loc,nsx*nn) :: x
    integer, intent(inout) :: req
    !
    integer :: ierr
    integer :: istat(MPI_STATUS_SIZE)
    !-------------------------------------------------------

#ifndef NO_ASYNC_MPI

    call MPI_WAIT(req, &
         istat, &
         ierr)

#ifdef DISABLE_GPUDIRECT_MPI
#if defined(OMPGPU)
!$omp target update to(x)
#elif defined(_OPENACC)
!$acc update device(x)
#endif
#endif

#endif
   !else, noop

  end subroutine parallel_slib_r_nc32_wait

!=========================================================

  ! Automaticallt use CPU or GPU version, based on presence of ACC
  subroutine parallel_slib_f_fd(nels1,nels2,nels3,x,xt)

    use mpi

    !-------------------------------------------------------
    implicit none
    !
    integer, intent(in) :: nels1,nels2,nels3
    complex, intent(inout), dimension(nels1,nels2,nels3,nk_loc*nn) :: x
    complex, intent(inout), dimension(nels1,nels2,nels3,nk_loc*nn) :: xt
    !
    integer :: ierr
    !-------------------------------------------------------

#ifdef DISABLE_GPUDIRECT_MPI
#if defined(OMPGPU)
!$omp target update from(x)
#elif defined(_OPENACC)
!$acc update host(x)
#endif
#else
#if defined(OMPGPU)
!$omp target data use_device_addr(x,xt)
#elif defined(_OPENACC)
!$acc host_data use_device(x,xt)
#endif
#endif

    call MPI_ALLTOALL(x, &
         nels1*nels2*nels3*nk_loc, &
         MPI_DOUBLE_COMPLEX, &
         xt, &
         nels1*nels2*nels3*nk_loc, &
         MPI_DOUBLE_COMPLEX, &
         slib_comm, &
         ierr)

#ifdef DISABLE_GPUDIRECT_MPI
#if defined(OMPGPU)
!$omp target update to(xt)
#elif defined(_OPENACC)
!$acc update device(xt)
#endif
#else
#if defined(OMPGPU)
!$omp end target data
#elif defined(_OPENACC)
!$acc end host_data
#endif
#endif

  end subroutine parallel_slib_f_fd

  subroutine parallel_slib_f_fd_async(nels1,nels2,nels3,x,xt,req)
    use mpi
    !-------------------------------------------------------
    implicit none
    !
    integer, intent(in) :: nels1,nels2,nels3
    complex, intent(inout), dimension(nels1,nels2,nels3,nk_loc*nn) :: x
    complex, intent(inout), dimension(nels1,nels2,nels3,nk_loc*nn) :: xt
    integer, intent(inout) :: req
    !
    integer :: ierr
    !-------------------------------------------------------

#ifdef NO_ASYNC_MPI
    call parallel_slib_f_fd(nels1,nels2,nels3,x,xt)
#else

#ifdef DISABLE_GPUDIRECT_MPI
#if defined(OMPGPU)
!$omp target update from(x)
#elif defined(_OPENACC)
!$acc update host(x)
#endif
#else
#if defined(OMPGPU)
!$omp target data use_device_addr(x,xt)
#elif defined(_OPENACC)
!$acc host_data use_device(x,xt)
#endif
#endif

   call MPI_IALLTOALL(x, &
         nels1*nels2*nels3*nk_loc, &
         MPI_DOUBLE_COMPLEX, &
         xt, &
         nels1*nels2*nels3*nk_loc, &
         MPI_DOUBLE_COMPLEX, &
         slib_comm, &
         req, &
         ierr)

#ifndef DISABLE_GPUDIRECT_MPI
#if defined(OMPGPU)
!$omp end target data
#elif defined(_OPENACC)
!$acc end host_data
#endif
#endif

#endif

  end subroutine parallel_slib_f_fd_async

  subroutine parallel_slib_f_fd32(nels1,nels2,nels3,x,xt)

    use mpi
    use, intrinsic :: iso_fortran_env

    !-------------------------------------------------------
    implicit none
    !
    integer, intent(in) :: nels1,nels2,nels3
    complex(KIND=REAL32), intent(inout), dimension(nels1,nels2,nels3,nk_loc*nn) :: x
    complex(KIND=REAL32), intent(inout), dimension(nels1,nels2,nels3,nk_loc*nn) :: xt
    !
    integer :: ierr
    !-------------------------------------------------------

#ifdef DISABLE_GPUDIRECT_MPI
#if defined(OMPGPU)
!$omp target update from(x)
#elif defined(_OPENACC)
!$acc update host(x)
#endif
#else
#if defined(OMPGPU)
!$omp target data use_device_addr(x,xt)
#elif defined(_OPENACC)
!$acc host_data use_device(x,xt)
#endif
#endif

    call MPI_ALLTOALL(x, &
         nels1*nels2*nels3*nk_loc, &
         MPI_COMPLEX, &
         xt, &
         nels1*nels2*nels3*nk_loc, &
         MPI_COMPLEX, &
         slib_comm, &
         ierr)

#ifdef DISABLE_GPUDIRECT_MPI
#if defined(OMPGPU)
!$omp target update to(xt)
#elif defined(_OPENACC)
!$acc update device(xt)
#endif
#else
#if defined(OMPGPU)
!$omp end target data
#elif defined(_OPENACC)
!$acc end host_data
#endif
#endif

  end subroutine parallel_slib_f_fd32

  subroutine parallel_slib_f_fd32_async(nels1,nels2,nels3,x,xt,req)
    use mpi
    use, intrinsic :: iso_fortran_env
    !-------------------------------------------------------
    implicit none
    !
    integer, intent(in) :: nels1,nels2,nels3
    complex(KIND=REAL32), intent(inout), dimension(nels1,nels2,nels3,nk_loc*nn) :: x
    complex(KIND=REAL32), intent(inout), dimension(nels1,nels2,nels3,nk_loc*nn) :: xt
    integer, intent(inout) :: req
    !
    integer :: ierr
    !-------------------------------------------------------

#ifdef NO_ASYNC_MPI
    call parallel_slib_f_fd32(nels1,nels2,nels3,x,xt)
#else

#ifdef DISABLE_GPUDIRECT_MPI
#if defined(OMPGPU)
!$omp target update from(x)
#elif defined(_OPENACC)
!$acc update host(x)
#endif
#else
#if defined(OMPGPU)
!$omp target data use_device_addr(x,xt)
#elif defined(_OPENACC)
!$acc host_data use_device(x,xt)
#endif
#endif

   call MPI_IALLTOALL(x, &
         nels1*nels2*nels3*nk_loc, &
         MPI_COMPLEX, &
         xt, &
         nels1*nels2*nels3*nk_loc, &
         MPI_COMPLEX, &
         slib_comm, &
         req, &
         ierr)

#ifndef DISABLE_GPUDIRECT_MPI
#if defined(OMPGPU)
!$omp end target data
#elif defined(_OPENACC)
!$acc end host_data
#endif
#endif

#endif

  end subroutine parallel_slib_f_fd32_async

  ! require x and xt to ensure they exist until this finishes
  subroutine parallel_slib_f_fd_wait(nels1,nels2,nels3,x,xt,req)
    use mpi
    !-------------------------------------------------------
    implicit none
    !
    integer, intent(in) :: nels1,nels2,nels3
    complex, intent(inout), dimension(nels1,nels2,nels3,nk_loc*nn) :: x
    complex, intent(inout), dimension(nels1,nels2,nels3,nk_loc*nn) :: xt
    integer, intent(inout) :: req
    !
    integer :: ierr
    integer :: istat(MPI_STATUS_SIZE)
    !-------------------------------------------------------

#ifndef NO_ASYNC_MPI

    call MPI_WAIT(req, &
         istat, &
         ierr)

#ifdef DISABLE_GPUDIRECT_MPI
#if defined(OMPGPU)
!$omp target update to(xt)
#elif defined(_OPENACC)
!$acc update device(xt)
#endif
#endif

#endif
  ! else noop


  end subroutine parallel_slib_f_fd_wait

  ! require x and xt to ensure they exist until this finishes
  subroutine parallel_slib_f_fd32_wait(nels1,nels2,nels3,x,xt,req)
    use mpi
    use, intrinsic :: iso_fortran_env
    !-------------------------------------------------------
    implicit none
    !
    integer, intent(in) :: nels1,nels2,nels3
    complex(KIND=REAL32), intent(inout), dimension(nels1,nels2,nels3,nk_loc*nn) :: x
    complex(KIND=REAL32), intent(inout), dimension(nels1,nels2,nels3,nk_loc*nn) :: xt
    integer, intent(inout) :: req
    !
    integer :: ierr
    integer :: istat(MPI_STATUS_SIZE)
    !-------------------------------------------------------

#ifndef NO_ASYNC_MPI

    call MPI_WAIT(req, &
         istat, &
         ierr)

#ifdef DISABLE_GPUDIRECT_MPI
#if defined(OMPGPU)
!$omp target update to(xt)
#elif defined(_OPENACC)
!$acc update device(xt)
#endif
#endif

#endif
  ! else noop


  end subroutine parallel_slib_f_fd32_wait

!=========================================================

  subroutine parallel_slib_distribute_real(nels1,nels2,nels3,nels4,nels5,x)

    use mpi

    !-------------------------------------------------------
    implicit none
    !
    integer, intent(in) :: nels1,nels2,nels3,nels4,nels5
    real, intent(inout), dimension(nels1,nels2,nels3,nels4,nels5,nn) :: x
    !
    integer :: nels
    integer :: ierr
    !-------------------------------------------------------

    nels = nels1*nels2*nels3*nels4*nels5
#ifdef DISABLE_GPUDIRECT_MPI
#if defined(OMPGPU)
!$omp target update from(x)
#elif defined(_OPENACC)
!$acc update host(x)
#endif
#else
#if defined(OMPGPU)
!$omp target data use_device_addr(x)
#elif defined(_OPENACC)
!$acc host_data use_device(x)
#endif
#endif

    call MPI_ALLTOALL(MPI_IN_PLACE, &
         nels, &
         MPI_DOUBLE, &
         x, &
         nels, &
         MPI_DOUBLE, &
         slib_comm, &
         ierr)

#ifdef DISABLE_GPUDIRECT_MPI
#if defined(OMPGPU)
!$omp target update to(x)
#elif defined(_OPENACC)
!$acc update device(x)
#endif
#else
#if defined(OMPGPU)
!$omp end target data
#elif defined(_OPENACC)
!$acc end host_data
#endif
#endif

  end subroutine parallel_slib_distribute_real

  subroutine parallel_slib_distribute_real32(nels1,nels2,nels3,nels4,nels5,x)

    use mpi
    use, intrinsic :: iso_fortran_env

    !-------------------------------------------------------
    implicit none
    !
    integer, intent(in) :: nels1,nels2,nels3,nels4,nels5
    real(KIND=REAL32), intent(inout), dimension(nels1,nels2,nels3,nels4,nels5,nn) :: x
    !
    integer :: nels
    integer :: ierr
    !-------------------------------------------------------

    nels = nels1*nels2*nels3*nels4*nels5
#ifdef DISABLE_GPUDIRECT_MPI
#if defined(OMPGPU)
!$omp target update from(x)
#elif defined(_OPENACC)
!$acc update host(x)
#endif
#else
#if defined(OMPGPU)
!$omp target data use_device_addr(x)
#elif defined(_OPENACC)
!$acc host_data use_device(x)
#endif
#endif

    call MPI_ALLTOALL(MPI_IN_PLACE, &
         nels, &
         MPI_REAL, &
         x, &
         nels, &
         MPI_REAL, &
         slib_comm, &
         ierr)

#ifdef DISABLE_GPUDIRECT_MPI
#if defined(OMPGPU)
!$omp target update to(x)
#elif defined(_OPENACC)
!$acc update device(x)
#endif
#else
#if defined(OMPGPU)
!$omp end target data
#elif defined(_OPENACC)
!$acc end host_data
#endif
#endif

  end subroutine parallel_slib_distribute_real32

!=========================================================

  subroutine parallel_lib_nj_loc(nj_loc_in)

    integer, intent(inout) :: nj_loc_in

    nj_loc_in = nj_loc
    
  end subroutine parallel_lib_nj_loc

  subroutine parallel_lib_clean
    implicit none
    
    if(allocated(fsendf)) then
#if defined(OMPGPU)
!$omp target exit data map(release:fsendf)
#elif defined(_OPENACC)
!$acc exit data delete(fsendf)
#endif
       deallocate(fsendf)
    endif

    if(allocated(fsendr)) then
#if defined(OMPGPU)
!$omp target exit data map(release:fsendr)
#elif defined(_OPENACC)
!$acc exit data delete(fsendr)
#endif
       deallocate(fsendr)
    endif

    if(allocated(fsendr_real))     deallocate(fsendr_real)
    
  end subroutine parallel_lib_clean
  
end module parallel_lib
