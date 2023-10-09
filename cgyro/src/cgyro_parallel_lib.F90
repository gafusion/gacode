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

  ! lib

  integer :: nproc,iproc
  integer, private :: ni,nj
  integer, private :: ni_loc
  integer, private :: nj_loc
  integer, private :: nk_loc
  integer, private :: nk1,nk2
  integer, private :: lib_comm
  integer, private :: nsend

  real, dimension(:,:,:,:), allocatable, private :: fsendr_real

  ! (expose these)
  complex, dimension(:,:,:,:), allocatable :: fsendf
  complex, dimension(:,:,:,:), allocatable :: fsendr

  ! clib

  integer, private :: ncproc,icproc
  integer, private :: ns1,ns2
  integer, private :: clib_comm

  ! slib

  integer :: nsproc,isproc
  integer, private :: nn
  integer, private :: nexch
  integer, private :: nkeep
  integer, private :: nsplit
  integer, private :: slib_comm

  integer, private, parameter :: default_size = 1

contains

#ifdef DISABLE_GPUDIRECT_MPI

#if defined(OMPGPU)

#define cpl_use_device(fin,fout) \
!$omp target update from(fin)
#define cpl_release_device(fin,fout) \
!$omp target update to(fout)
  ! used for async, so no copy at this point
#define cpl_unbind_device(fin,fout)

#define cpl_finalize_device(fin,fout) \
!$omp target update to(fout)

#elif defined(_OPENACC)

#define cpl_use_device(fin,fout) \
!$acc update host(fin)
#define cpl_release_device(fin,fout) \
!$acc update device(fout)
  ! used for async, so no copy at this point
#define cpl_unbind_device(fin,fout)

#define cpl_finalize_device(fin,fout) \
!$acc update device(fout)

#else
  ! no devices, no-ops only
#define cpl_use_device(fin,fout) 
#define cpl_release_device(fin,fout) 
#define cpl_unbind_device(fin,fout)
#define cpl_finalize_device(fin,fout) 

#endif

#else

#if defined(OMPGPU)

#define cpl_use_device(fin,fout) \
!$omp target data use_device_ptr(fin,fout)
#define cpl_release_device(fin,fout) \
!$omp end target data
#define cpl_unbind_device(fin,fout) \
!$omp end target data
  ! no-op, as there was no copying involved
#define cpl_finalize_device(fin,fout)

#elif defined(_OPENACC)

#define cpl_use_device(fin,fout) \
!$acc host_data use_device(fin,fout)
#define cpl_release_device(fin,fout) \
!$acc end host_data
#define cpl_unbind_device(fin,fout) \
!$acc end host_data
  ! no-op, as there was no copying involved
#define cpl_finalize_device(fin,fout)

#else
  ! no devices, no-ops only
#define cpl_use_device(fin,fout) 
#define cpl_release_device(fin,fout) 
#define cpl_unbind_device(fin,fout)
#define cpl_finalize_device(fin,fout) 

#endif

#endif

  !=========================================================
  !  parallel_lib_f -> f(ni_loc,nj) -> g(nj_loc,ni) 
  !  parallel_lib_r -> g(nj_loc,ni) -> f(ni_loc,nj)
  !=========================================================

  subroutine parallel_lib_init(ni_in,nj_in,nk1_in,nk_loc_in,ni_loc_out,nj_loc_out,comm)

    use mpi

    implicit none

    integer, intent(in) :: ni_in,nj_in
    integer, intent(in) :: nk1_in,nk_loc_in
    integer, intent(in) :: comm
    integer, intent(inout) :: ni_loc_out,nj_loc_out
    integer, external :: parallel_dim
    integer :: ierr

    lib_comm = comm

    call MPI_COMM_RANK(lib_comm,iproc,ierr)
    call MPI_COMM_SIZE(lib_comm,nproc,ierr)

    ni = ni_in
    nj = nj_in

    ! parallel_dim(x,y) ~= x/y
    ni_loc = parallel_dim(ni,nproc)
    nj_loc = parallel_dim(nj,nproc)
    nk_loc = nk_loc_in

    ! nk1_in is typically iproc*nk_loc, but not always
    nk1 = nk1_in
    nk2 = nk1 + nk_loc -1

    ni_loc_out = ni_loc
    nj_loc_out = nj_loc

    nsend = nj_loc*ni_loc*nk_loc

    allocate(fsendf(nj_loc,nk1:nk2,ni_loc,nproc))
    allocate(fsendr(ni_loc,nk1:nk2,nj_loc,nproc))
    if (.not. allocated(fsendr_real)) allocate(fsendr_real(ni_loc,nk1:nk2,nj_loc,nproc))

#if defined(OMPGPU)
!$omp target enter data map(alloc:fsendf,fsendr)
!$omp target enter data map(to:nproc,nk1,nk2,ni_loc)
#elif defined(_OPENACC)
!$acc enter data create(fsendf,fsendr)
!$acc enter data copyin(nproc,nk1,nk2,ni_loc)
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

    cpl_use_device(fsendf,ft)

    call MPI_ALLTOALL(fsendf, &
         nsend, &
         MPI_DOUBLE_COMPLEX, &
         ft, &
         nsend, &
         MPI_DOUBLE_COMPLEX, &
         lib_comm, &
         ierr)

    cpl_release_device(fsendf,ft)

  end subroutine parallel_lib_f_i_do_gpu

  !=========================================================

  subroutine parallel_lib_r_do(f)

    use mpi

    implicit none

    complex, intent(inout), dimension(:,:,:) :: f
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

    complex, intent(inout), dimension(:,:,:) :: f

    integer :: ierr

    cpl_use_device(fsendr,f)

    call MPI_ALLTOALL(fsendr, &
         nsend, &
         MPI_DOUBLE_COMPLEX, &
         f, &
         nsend, &
         MPI_DOUBLE_COMPLEX, &
         lib_comm, &
         ierr)

    cpl_release_device(fsendr,f)

  end subroutine parallel_lib_r_do_gpu

  !=========================================================

  subroutine parallel_lib_r_pack(ft)

    use mpi

    implicit none

    complex, intent(in), dimension(:,:,:) :: ft
    integer :: j_loc,i,j,k,j1,j2,itor

    j1 = 1+iproc*nj_loc
    j2 = (1+iproc)*nj_loc

!$omp parallel do collapse(2) if (size(fsendr) >= default_size) default(none) &
!$omp& shared(nproc,j1,j2,ni_loc,nk1,nk2) &
!$omp& private(j,j_loc,i) &
!$omp& shared(ft,fsendr)
    do k=1,nproc
     do itor=nk1,nk2
       do j=j1,j2
          j_loc = j-j1+1 
          do i=1,ni_loc
             fsendr(i,itor,j_loc,k) = ft(j_loc,i+(k-1)*ni_loc,1+(itor-nk1))
          enddo
       enddo
     enddo
    enddo

  end subroutine parallel_lib_r_pack

  !=========================================================

  subroutine parallel_lib_rtrans_pack(fin)
! -----------------------------------------
! transpose version of parallel_lib_r_pack(fin)
! -----------------------------------------
    use mpi

    implicit none

    complex, intent(in), dimension(:,:,:) :: fin
    integer :: j_loc,i,j,k,j1,j2,itor

    j1 = 1+iproc*nj_loc
    j2 = (1+iproc)*nj_loc

!$omp parallel do collapse(2) if (size(fsendr) >= default_size) default(none) &
!$omp& shared(nproc,j1,j2,ni_loc,nk1,nk2) &
!$omp& private(j,j_loc,i) &
!$omp& shared(fin,fsendr)
    do k=1,nproc
     do itor=nk1,nk2
       do j=j1,j2
          j_loc = j-j1+1 
          do i=1,ni_loc
             fsendr(i,itor,j_loc,k) = fin(i+(k-1)*ni_loc,j_loc,1+(itor-nk1))
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

    j1 = 1+iproc*nj_loc
    j2 = (1+iproc)*nj_loc
#if defined(OMPGPU)
!$omp target teams distribute parallel do simd collapse(4) &
!$omp&  private(j_loc) map(to:j1,j2)
#else
!$acc parallel loop collapse(4) gang vector independent private(j_loc) &
!$acc&         present(fsendr,fin) present(nproc,nk1,nk2,ni_loc) &
!$acc&         copyin(j1,j2) default(none)
#endif
    do k=1,nproc
     do itor=nk1,nk2
       do j=j1,j2
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
    complex, intent(inout), dimension(:,:,:) :: f

    call parallel_lib_rtrans_pack(fin)
    call parallel_lib_r_do(f)

  end subroutine parallel_lib_rtrans

  !=========================================================

  subroutine parallel_lib_rtrans_real(fin,f)

    use mpi

    implicit none

    real, intent(in), dimension(:,:,:) :: fin
    real, intent(inout), dimension(:,:,:) :: f
    integer :: ierr,j_loc,i,j,k,j1,j2,itor

    j1 = 1+iproc*nj_loc
    j2 = (1+iproc)*nj_loc

!$omp parallel do collapse(2) if (size(fsendr_real) >= default_size) default(none) &
!$omp& shared(nproc,j1,j2,ni_loc,nk1,nk2) &
!$omp& private(j,j_loc,i) &
!$omp& shared(fin,fsendr_real)
    do k=1,nproc
     do itor=nk1,nk2
       do j=j1,j2
          j_loc = j-j1+1
          do i=1,ni_loc
             fsendr_real(i,itor,j_loc,k) = fin(i+(k-1)*ni_loc,j_loc,1+(itor-nk1)) 
          enddo
       enddo
     enddo
    enddo

    call MPI_ALLTOALL(fsendr_real, &
         nsend, &
         MPI_DOUBLE_PRECISION,&
         f, &
         nsend, &
         MPI_DOUBLE_PRECISION, &
         lib_comm, &
         ierr)

  end subroutine parallel_lib_rtrans_real

  !=========================================================

  subroutine parallel_lib_sum_field(field_loc,field)

    use mpi

    implicit none

    complex, intent(in), dimension(:,:,:) :: field_loc
    complex, intent(inout), dimension(:,:,:) :: field
    integer :: ierr


    call MPI_ALLREDUCE(field_loc(:,:,:),&
         field(:,:,:),&
         size(field(:,:,:)),&
         MPI_DOUBLE_COMPLEX,&
         MPI_SUM,&
         lib_comm,&
         ierr)

  end subroutine parallel_lib_sum_field

  !=========================================================

  subroutine parallel_lib_sum_field_gpu(field_loc,field)

    use mpi

    implicit none

    complex, intent(in), dimension(:,:,:) :: field_loc
    complex, intent(inout), dimension(:,:,:) :: field
    integer :: ierr

    cpl_use_device(field_loc,field)

    call MPI_ALLREDUCE(field_loc(:,:,:),&
          field(:,:,:),&
          size(field(:,:,:)),&
          MPI_DOUBLE_COMPLEX,&
          MPI_SUM,&
          lib_comm,&
          ierr)

    cpl_release_device(field_loc,field)

  end subroutine parallel_lib_sum_field_gpu

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

    complex, intent(in), dimension(:,:,:) :: upwind_loc
    complex, intent(inout), dimension(:,:,:) :: upwind
    integer :: ierr

    cpl_use_device(upwind_loc,upwind)

    call MPI_ALLREDUCE(upwind_loc(:,:,:),&
          upwind(:,:,:),&
          size(upwind(:,:,:)),&
          MPI_DOUBLE_COMPLEX,&
          MPI_SUM,&
          clib_comm,&
          ierr)

    cpl_release_device(upwind_loc,upwind)

  end subroutine parallel_clib_sum_upwind

  !=========================================================
  subroutine parallel_clib_sum_upwind32(upwind_loc,upwind)

    use mpi
    use, intrinsic :: iso_fortran_env

    implicit none

    complex(KIND=REAL32), intent(in), dimension(:,:,:) :: upwind_loc
    complex(KIND=REAL32), intent(inout), dimension(:,:,:) :: upwind
    integer :: ierr

    cpl_use_device(upwind_loc,upwind)

    call MPI_ALLREDUCE(upwind_loc(:,:,:),&
          upwind(:,:,:),&
          size(upwind(:,:,:)),&
          MPI_COMPLEX,&
          MPI_SUM,&
          clib_comm,&
          ierr)

    cpl_release_device(upwind_loc,upwind)

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

  subroutine parallel_slib_f_nc(x,xt)

    use mpi

    !-------------------------------------------------------
    implicit none
    !
    complex, intent(in), dimension(nkeep,nk_loc,nsplit*nn) :: x
    complex, intent(inout), dimension(nkeep,nk_loc,nsplit,nn) :: xt
    !
    integer :: ierr
    !-------------------------------------------------------

    cpl_use_device(x,xt)

    call MPI_ALLTOALL(x, &
         nkeep*nk_loc*nsplit, &
         MPI_DOUBLE_COMPLEX, &
         xt, &
         nkeep*nk_loc*nsplit, &
         MPI_DOUBLE_COMPLEX, &
         slib_comm, &
         ierr)

    cpl_release_device(x,xt)

  end subroutine parallel_slib_f_nc

  subroutine parallel_slib_f_nc_async(x,xt,req)
    use mpi
    !-------------------------------------------------------
    implicit none
    !
    complex, intent(in), dimension(nkeep,nk_loc,nsplit*nn) :: x
    complex, intent(inout), dimension(nkeep,nk_loc,nsplit,nn) :: xt
    integer, intent(inout) :: req
    !
    integer :: ierr
    !-------------------------------------------------------

#ifdef NO_ASYNC_MPI
   call parallel_slib_f_nc(x,xt)

#else

    cpl_use_device(x,xt)

   call MPI_IALLTOALL(x, &
         nkeep*nk_loc*nsplit, &
         MPI_DOUBLE_COMPLEX, &
         xt, &
         nkeep*nk_loc*nsplit, &
         MPI_DOUBLE_COMPLEX, &
         slib_comm, &
         req, &
         ierr)

    cpl_unbind_device(x,xt)

#endif

  end subroutine parallel_slib_f_nc_async

  ! require x and xt to ensure they exist until this finishes
  subroutine parallel_slib_f_nc_wait(x,xt,req)
    use mpi
    !-------------------------------------------------------
    implicit none
    !
    complex, intent(in), dimension(nkeep,nk_loc,nsplit*nn) :: x
    complex, intent(inout), dimension(nkeep,nk_loc,nsplit,nn) :: xt
    integer, intent(inout) :: req
    !
    integer :: ierr
    integer :: istat(MPI_STATUS_SIZE)
    !-------------------------------------------------------

#ifndef NO_ASYNC_MPI

    call MPI_WAIT(req, &
         istat, &
         ierr)

    cpl_finalize_device(x,xt)

#endif
   !else, noop

  end subroutine parallel_slib_f_nc_wait

 !=========================================================

  subroutine parallel_slib_r_nc (xt,x)
    use mpi
    !-------------------------------------------------------
    implicit none
    !
    complex, intent(in), dimension(nkeep,nk_loc,nsplit,nn) :: xt
    complex, intent(inout), dimension(nkeep,nk_loc,nsplit*nn) :: x
    !
    integer :: ierr
    !-------------------------------------------------------

    cpl_use_device(xt,x)

    call MPI_ALLTOALL(xt, &
         nkeep*nk_loc*nsplit, &
         MPI_DOUBLE_COMPLEX, &
         x, &
         nkeep*nk_loc*nsplit, &
         MPI_DOUBLE_COMPLEX, &
         slib_comm, &
         ierr)

    cpl_release_device(xt,x)

  end subroutine parallel_slib_r_nc

!=========================================================

  ! Automaticallt use CPU or GPU version, based on presence of ACC
  subroutine parallel_slib_f_fd(nels1,nels2,nels3,x,xt)

    use mpi

    !-------------------------------------------------------
    implicit none
    !
    integer, intent(in) :: nels1,nels2,nels3
    complex, intent(in), dimension(nels1,nels2,nels3,nk_loc*nn) :: x
    complex, intent(inout), dimension(nels1,nels2,nels3,nk_loc*nn) :: xt
    !
    integer :: ierr
    !-------------------------------------------------------

    cpl_use_device(x,xt)

    call MPI_ALLTOALL(x, &
         nels1*nels2*nels3*nk_loc, &
         MPI_DOUBLE_COMPLEX, &
         xt, &
         nels1*nels2*nels3*nk_loc, &
         MPI_DOUBLE_COMPLEX, &
         slib_comm, &
         ierr)

    cpl_release_device(x,xt)

  end subroutine parallel_slib_f_fd

  subroutine parallel_slib_f_fd_async(nels1,nels2,nels3,x,xt,req)
    use mpi
    !-------------------------------------------------------
    implicit none
    !
    integer, intent(in) :: nels1,nels2,nels3
    complex, intent(in), dimension(nels1,nels2,nels3,nk_loc*nn) :: x
    complex, intent(inout), dimension(nels1,nels2,nels3,nk_loc*nn) :: xt
    integer, intent(inout) :: req
    !
    integer :: ierr
    !-------------------------------------------------------

#ifdef NO_ASYNC_MPI
    call parallel_slib_f_fd(nels1,nels2,nels3,x,xt)
#else

    cpl_use_device(x,xt)

   call MPI_IALLTOALL(x, &
         nels1*nels2*nels3*nk_loc, &
         MPI_DOUBLE_COMPLEX, &
         xt, &
         nels1*nels2*nels3*nk_loc, &
         MPI_DOUBLE_COMPLEX, &
         slib_comm, &
         req, &
         ierr)

    cpl_unbind_device(x,xt)

#endif

  end subroutine parallel_slib_f_fd_async

  ! require x and xt to ensure they exist until this finishes
  subroutine parallel_slib_f_fd_wait(nels1,nels2,nels3,x,xt,req)
    use mpi
    !-------------------------------------------------------
    implicit none
    !
    integer, intent(in) :: nels1,nels2,nels3
    complex, intent(in), dimension(nels1,nels2,nels3,nk_loc*nn) :: x
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

    cpl_finalize_device(x,xt)

#endif
  ! else noop


  end subroutine parallel_slib_f_fd_wait

!=========================================================

  ! x is logically (nels,nn)
  subroutine parallel_slib_distribute_real(nels,x)

    use mpi

    !-------------------------------------------------------
    implicit none
    !
    integer, intent(in) :: nels
    real, intent(inout), dimension(*) :: x
    !
    integer :: ierr
    !-------------------------------------------------------

#ifdef DISABLE_GPUDIRECT_MPI

#if defined(OMPGPU)
!$omp target update from(x(1:nels))
#elif defined(_OPENACC)
!$acc update host(x(1:nels))
#endif

#else

#if defined(OMPGPU)
!$omp target data use_device_ptr(x)
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
!$omp target update to(x(1:nels))
#elif defined(_OPENACC)
!$acc update device(x(1:nels))
#endif

#else

#if defined(OMPGPU)
!$omp end target data
#elif defined(_OPENACC)
!$acc end host_data
#endif

#endif

  end subroutine parallel_slib_distribute_real

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

