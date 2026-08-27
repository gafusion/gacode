!-----------------------------------------------------------------
! cgyro_coll_data.F90
!
! PURPOSE:
!  CGYRO collision matrix storage and access helpers.
!-----------------------------------------------------------------

module cgyro_coll_data

  use, intrinsic :: iso_fortran_env

  implicit none

  !
  ! Collision operator
  integer :: n_low_energy

 ! Unlike most other arrays, last dimension is simply nt_loc
  ! due to limits of pointer handling
  real, dimension(:,:,:,:), pointer :: cmat ! only used if collision_precision_mode=0 & 64
  real(KIND=REAL32), dimension(:,:,:,:), pointer :: cmat_fp32 ! only used if collision_precision_mode=1 & 32
  !
  real(KIND=REAL32), dimension(:,:,:,:,:,:), allocatable :: cmat_stripes ! only used if collision_precision_mode=1
  real(KIND=REAL32), dimension(:,:,:,:,:,:), allocatable :: cmat_e1 ! only used if collision_precision_mode=1
  real, dimension(:,:,:,:,:,:), allocatable :: cmat_simple ! only used in collision_model=5

contains

!
! Avoid having two copies of cmat, one in CPU and one in GPU memory
! On GPU systems under OMPGPU, keep it only in GPU memory (omp_target_alloc)
! ON CPU-only systems and OpenACC, keep the old logic
!

#ifdef OMPGPU
subroutine allocate_cmat_gpu_fp32(nv,nc_loc_coll,nt1,nt2)
    use iso_c_binding
    use omp_lib
    implicit none

    ! ----------------------
    integer, intent(in)          :: nv,nc_loc_coll,nt1,nt2

    ! Pointer to array, only need raw pointer locally, can recover
    type(c_ptr) :: c_ptr_buffer
    integer :: dev_id
    ! Use c_sizeof to avoid integer wraparound
    integer(c_size_t) :: total_bytes
    integer(c_size_t) :: nv2
    integer :: nt_loc

    nv2 = nv*nv
    nt_loc = nt2-nt1+1

    total_bytes = (nv2*nc_loc_coll)*nt_loc*4 ! fp32 =  4 bytes

    ! OMP supports multi-GPU setups, we only support 1-GPU ones
    dev_id = omp_get_default_device()
    c_ptr_buffer = omp_target_alloc(total_bytes, dev_id)

    if (.not. c_associated(c_ptr_buffer)) then
         ! Catastrophic error, do not even try to process it cleanly
         write(*,*) "Error: allocate_cmat_fp32 failed to allocate memory."
        stop
    end if

    call c_f_pointer(c_ptr_buffer, cmat_fp32, [nv, nv, nc_loc_coll, nt_loc])

end subroutine allocate_cmat_gpu_fp32

subroutine allocate_cmat_gpu(nv,nc_loc_coll,nt1,nt2)
    use iso_c_binding
    use omp_lib
    implicit none

    ! ----------------------
    integer, intent(in)          :: nv,nc_loc_coll,nt1,nt2

    ! Pointer to array, only need raw pointer locally, can recover
    type(c_ptr) :: c_ptr_buffer
    integer :: dev_id
    ! Use c_sizeof to avoid integer wraparound
    integer(c_size_t) :: total_bytes
    integer(c_size_t) :: nv2
    integer :: nt_loc

    nv2 = nv*nv
    nt_loc = nt2-nt1+1

    total_bytes = (nv2*nc_loc_coll)*nt_loc*8 ! fp64 =  8 bytes

    ! OMP supports multi-GPU setups, we only support 1-GPU ones
    dev_id = omp_get_default_device()
    c_ptr_buffer = omp_target_alloc(total_bytes, dev_id)

    if (.not. c_associated(c_ptr_buffer)) then
         ! Catastrophic error, do not even try to process it cleanly
         write(*,*) "Error: allocate_cmat failed to allocate memory."
        stop
    end if

    call c_f_pointer(c_ptr_buffer, cmat, [nv, nv, nc_loc_coll, nt_loc])

end subroutine allocate_cmat_gpu

subroutine deallocate_cmat_gpu_fp32
    use iso_c_binding
    use omp_lib
    implicit none
    type(c_ptr) :: c_ptr_buffer
    integer :: dev_id

    if (associated(cmat_fp32)) then
      c_ptr_buffer = c_loc(cmat_fp32(1,1,1,1))

      ! OMP supports multi-GPU setups, we only support 1-GPU ones
      ! Note: Must be the same as the one used in alloc
      dev_id = omp_get_default_device()
      call omp_target_free(c_ptr_buffer, dev_id)
    endif

end subroutine

subroutine deallocate_cmat_gpu
    use iso_c_binding
    use omp_lib
    implicit none
    type(c_ptr) :: c_ptr_buffer
    integer :: dev_id

    if (associated(cmat)) then
      c_ptr_buffer = c_loc(cmat(1,1,1,1))

      ! OMP supports multi-GPU setups, we only support 1-GPU ones
      ! Note: Must be the same as the one used in alloc
      dev_id = omp_get_default_device()
      call omp_target_free(c_ptr_buffer, dev_id)
    endif

end subroutine

subroutine copy_into_cmat_gpu_fp32(amat, ic_loc, itor)
    use iso_c_binding
    use omp_lib
    use cgyro_globals, only : nv, nc_loc_coll, nt1
    implicit none

    real(KIND=REAL32), target, intent(in)    :: amat(nv, nv)
    integer, intent(in) :: ic_loc, itor

    integer(c_size_t)   :: bytes_to_copy, dst_offset
    integer             :: host_id, device_id
    integer             :: ierr

    ! Size of the (nv, nv) block in bytes
    bytes_to_copy = int(nv, c_size_t) * int(nv, c_size_t) * c_sizeof(amat(1,1))

    host_id   = omp_get_initial_device()
    device_id = omp_get_default_device()

    ! 3. Calculate destination offset (0-based bytes)
    ! We are targeting: cmat_fp32(1, 1, ic_loc, itor-nt1+1)
    ! Offset = [(dim4_idx - 1) * (size_dim3 * size_dim2 * size_dim1) + (dim3_idx - 1) * (size_dim2 * size_dim1)] * element_size
    dst_offset = ( int(itor - nt1, c_size_t) * int(nc_loc_coll, c_size_t) * int(nv * nv, c_size_t) + &
                   int(ic_loc - 1, c_size_t) * int(nv * nv, c_size_t) ) * c_sizeof(amat(1,1))

    ierr = omp_target_memcpy( &
        c_loc(cmat_fp32),    & ! Destination (device pointer)
        c_loc(amat),        & ! Source (host pointer)
        bytes_to_copy,      & ! Length in bytes
        dst_offset,         & ! Offset in destination
        0_c_size_t,         & ! Offset in source
        device_id,          & ! Destination device
        host_id             & ! Source device
    )

    if (ierr /= 0) then
        ! Catastrophic error, do not even try to process it cleanly
        write(*,*) "Error in omp_target_memcpy: ", ierr
        stop
    end if

end subroutine

subroutine copy_into_cmat_gpu(amat, ic_loc, itor)
    use iso_c_binding
    use omp_lib
    use cgyro_globals, only : nv, nc_loc_coll, nt1
    implicit none

    real, target, intent(in)    :: amat(nv, nv)
    integer, intent(in) :: ic_loc, itor

    integer(c_size_t)   :: bytes_to_copy, dst_offset
    integer             :: host_id, device_id
    integer             :: ierr

    ! Size of the (nv, nv) block in bytes
    bytes_to_copy = int(nv, c_size_t) * int(nv, c_size_t) * c_sizeof(amat(1,1))

    host_id   = omp_get_initial_device()
    device_id = omp_get_default_device()

    ! 3. Calculate destination offset (0-based bytes)
    ! We are targeting: cmat(1, 1, ic_loc, itor-nt1+1)
    ! Offset = [(dim4_idx - 1) * (size_dim3 * size_dim2 * size_dim1) + (dim3_idx - 1) * (size_dim2 * size_dim1)] * element_size
    dst_offset = ( int(itor - nt1, c_size_t) * int(nc_loc_coll, c_size_t) * int(nv * nv, c_size_t) + &
                   int(ic_loc - 1, c_size_t) * int(nv * nv, c_size_t) ) * c_sizeof(amat(1,1))

    ierr = omp_target_memcpy( &
        c_loc(cmat),        & ! Destination (device pointer)
        c_loc(amat),        & ! Source (host pointer)
        bytes_to_copy,      & ! Length in bytes
        dst_offset,         & ! Offset in destination
        0_c_size_t,         & ! Offset in source
        device_id,          & ! Destination device
        host_id             & ! Source device
    )

    if (ierr /= 0) then
        ! Catastrophic error, do not even try to process it cleanly
        write(*,*) "Error in omp_target_memcpy: ", ierr
        stop
    end if

end subroutine

subroutine copy_from_cmat_gpu(amat, ic_loc, itor)
    use iso_c_binding
    use omp_lib
    use cgyro_globals, only : nv, nc_loc_coll, nt1
    implicit none

    real, target, intent(inout) :: amat(nv, nv)
    integer, intent(in) :: ic_loc, itor

    integer(c_size_t)   :: bytes_to_copy, src_offset
    integer             :: host_id, device_id
    integer             :: ierr

    ! Size of the (nv, nv) block in bytes
    bytes_to_copy = int(nv, c_size_t) * int(nv, c_size_t) * c_sizeof(amat(1,1))

    host_id   = omp_get_initial_device()
    device_id = omp_get_default_device()

    ! 3. Calculate cmat offset (0-based bytes)
    ! We are targeting: cmat(1, 1, ic_loc, itor-nt1+1)
    ! Offset = [(dim4_idx - 1) * (size_dim3 * size_dim2 * size_dim1) + (dim3_idx - 1) * (size_dim2 * size_dim1)] * element_size
    src_offset = ( int(itor - nt1, c_size_t) * int(nc_loc_coll, c_size_t) * int(nv * nv, c_size_t) + &
                   int(ic_loc - 1, c_size_t) * int(nv * nv, c_size_t) ) * c_sizeof(amat(1,1))

    ierr = omp_target_memcpy( &
        c_loc(amat),        & ! Destination (host pointer)
        c_loc(cmat),        & ! Source (device pointer)
        bytes_to_copy,      & ! Length in bytes
        0_c_size_t,         & ! Offset in destination
        src_offset,         & ! Offset in source
        host_id,            & ! Destination device
        device_id           & ! Source device
    )

    if (ierr /= 0) then
        ! Catastrophic error, do not even try to process it cleanly
        write(*,*) "Error in omp_target_memcpy: ", ierr
        stop
    end if

end subroutine

subroutine copy_from_cmat_gpu_all(cmat_cpu)
    use iso_c_binding
    use omp_lib
    use cgyro_globals, only : nv, nc_loc_coll, nt_loc
    implicit none

    real, target, intent(inout) :: cmat_cpu(nv, nv, nc_loc_coll, nt_loc)

    integer(c_size_t)   :: bytes_to_copy
    integer             :: host_id, device_id
    integer             :: ierr

    ! Size of the cmat block in bytes
    bytes_to_copy = int(nv, c_size_t) * int(nv, c_size_t) * &
                    int(nc_loc_coll, c_size_t) * int(nt_loc, c_size_t) * &
                    c_sizeof(cmat_cpu(1,1,1,1))

    host_id   = omp_get_initial_device()
    device_id = omp_get_default_device()

    ierr = omp_target_memcpy( &
        c_loc(cmat_cpu),    & ! Destination (host pointer)
        c_loc(cmat),        & ! Source (device pointer)
        bytes_to_copy,      & ! Length in bytes
        0_c_size_t,         & ! Offset in destination
        0_c_size_t,         & ! Offset in source
        host_id,            & ! Destination device
        device_id           & ! Source device
    )

    if (ierr /= 0) then
        ! Catastrophic error, do not even try to process it cleanly
        write(*,*) "Error in omp_target_memcpy: ", ierr
        stop
    end if

end subroutine

#endif
! =============  end OMPGPU =============

subroutine allocate_cmat_fp32(nv,nc_loc_coll,nt1,nt2)
    use cgyro_globals, only : gpu_bigmem_flag
    implicit none

    ! ----------------------
    integer, intent(in)          :: nv,nc_loc_coll,nt1,nt2

    integer :: nt_loc
    nt_loc = nt2-nt1+1

#ifdef OMPGPU
    if (gpu_bigmem_flag > 0) then
      ! we keep only the GPU buffer with BIGMEM
      call allocate_cmat_gpu_fp32(nv,nc_loc_coll,nt1,nt2)
    else
#else
    if (.TRUE.) then
#endif
      allocate(cmat_fp32(nv,nv,nc_loc_coll,nt_loc))
#if defined(_OPENACC)
!$acc enter data create(cmat_fp32) if (gpu_bigmem_flag > 0)
#endif
    endif

end subroutine allocate_cmat_fp32

subroutine allocate_cmat(nv,nc_loc_coll,nt1,nt2)
    use cgyro_globals, only : gpu_bigmem_flag
    implicit none

    ! ----------------------
    integer, intent(in)          :: nv,nc_loc_coll,nt1,nt2

    integer :: nt_loc
    nt_loc = nt2-nt1+1

#ifdef OMPGPU
    if (gpu_bigmem_flag > 0) then
      ! we keep only the GPU buffer with BIGMEM
      call allocate_cmat_gpu(nv,nc_loc_coll,nt1,nt2)
    else
#else
    if (.TRUE.) then
#endif
      allocate(cmat(nv,nv,nc_loc_coll,nt_loc))
#if defined(_OPENACC)
!$acc enter data create(cmat) if (gpu_bigmem_flag > 0)
#endif
    endif
end subroutine allocate_cmat

subroutine deallocate_cmat_fp32
    use cgyro_globals, only : gpu_bigmem_flag
    implicit none

    if (associated(cmat_fp32)) then
#ifdef OMPGPU
      if (gpu_bigmem_flag > 0) then
        call deallocate_cmat_gpu_fp32
      else
#else
      if (.TRUE.) then
#endif
        deallocate(cmat_fp32)
      endif
    endif

end subroutine

subroutine deallocate_cmat
    use cgyro_globals, only : gpu_bigmem_flag
    implicit none

    if (associated(cmat)) then
#ifdef OMPGPU
      if (gpu_bigmem_flag > 0) then
        call deallocate_cmat_gpu
      else
#else
      if (.TRUE.) then
#endif
        deallocate(cmat)
      endif
    endif

end subroutine

subroutine copy_into_cmat_fp32(amat,ic_loc,itor)
    use cgyro_globals, only : nv, nt1, gpu_bigmem_flag
    implicit none

    ! ----------------------
    real(KIND=REAL32), target, intent(in)    :: amat(nv,nv)
    integer, intent(in) :: ic_loc,itor

#ifdef OMPGPU
    if (gpu_bigmem_flag > 0) then
      ! cmat_fp32 in GPU memory only
      call copy_into_cmat_gpu_fp32(amat,ic_loc,itor)
    else
#else
    if (.TRUE.) then
#endif
      cmat_fp32(:,:,ic_loc,itor-nt1+1) = amat(:,:)
    endif

end subroutine

subroutine copy_into_cmat(amat,ic_loc,itor)
    use cgyro_globals, only : nv, nt1, gpu_bigmem_flag
    implicit none

    ! ----------------------
    real, target, intent(in)    :: amat(nv,nv)
    integer, intent(in) :: ic_loc,itor

#ifdef OMPGPU
    if (gpu_bigmem_flag > 0) then
      ! cmat in GPU memory only
      call copy_into_cmat_gpu(amat,ic_loc,itor)
    else
#else
    if (.TRUE.) then
#endif
      cmat(:,:,ic_loc,itor-nt1+1) = amat(:,:)
    endif

end subroutine

subroutine copy_from_cmat(amat,ic_loc,itor)
    use cgyro_globals, only : nv, nt1, gpu_bigmem_flag
    implicit none

    ! ----------------------
    real, target, intent(inout)    :: amat(nv,nv)
    integer, intent(in) :: ic_loc,itor

#ifdef OMPGPU
    if (gpu_bigmem_flag > 0) then
      ! cmat in GPU memory only
      call copy_from_cmat_gpu(amat,ic_loc,itor)
    else
#else
    if (.TRUE.) then
#endif
       amat(:,:) = cmat(:,:,ic_loc,itor-nt1+1)
    endif

end subroutine

subroutine copy_from_cmat_all(cmat_dest)
    use cgyro_globals, only : nv, nc_loc_coll, nt_loc, gpu_bigmem_flag
    implicit none

    ! ----------------------
    real, target, intent(inout) :: cmat_dest(nv,nv,nc_loc_coll,nt_loc)

#ifdef OMPGPU
    if (gpu_bigmem_flag > 0) then
      ! cmat in GPU memory only
      call copy_from_cmat_gpu_all(cmat_dest)
    else
#else
    if (.TRUE.) then
#endif
      cmat_dest(:,:,:,:) = cmat(:,:,:,:)
    endif

end subroutine

! -------------------------------------
! Initialization and cleanup of cmat
! -------------------------------------

subroutine cgyro_coll_data_init

    use cgyro_globals, only : collision_model, collision_precision_mode, &
         energy, nc_loc, nc_loc_coll, n_energy, n_species, n_theta, nt1, &
         nt2, nv, n_xi
    use cgyro_io, only : cgyro_error, cgyro_info

    implicit none

    character(len=128) :: msg
    integer :: ie

    n_low_energy = 0

    if (collision_model == 5) then
       if (nc_loc_coll/=nc_loc) then
          call cgyro_error("CMAT sharing not supported for COLLISION_MODEL 5")
          return
       endif
       allocate(cmat_simple(n_xi,n_xi,n_energy,n_species,n_theta,nt1:nt2))
    else
       if (collision_precision_mode == 1) then
          if (collision_model>=6) then
            call cgyro_error("COLLISION_PRECISION_MODE 1 not supported for COLLISION_MODEL >=6")
            return
          endif
          ! the lowest energy(s) has the most spread, so treat differently
          n_low_energy = 1
          do ie=2,n_energy
            if (energy(ie)<1.0e-2) then
              n_low_energy = ie
            endif
          enddo
          call allocate_cmat_fp32(nv,nc_loc_coll,nt1,nt2)
          allocate(cmat_stripes(n_xi,n_species,(n_low_energy+1):n_energy,n_xi,nc_loc_coll,nt1:nt2))
          allocate(cmat_e1(n_xi,n_species,n_low_energy,nv,nc_loc_coll,nt1:nt2))

          write (msg, "(A,I1,A)") "Using fp32 collision precision except e<=",n_low_energy," or same e&s."
          call cgyro_info(msg)
       else if (collision_precision_mode == 32) then
          if (collision_model>=6) then
            call cgyro_error("COLLISION_PRECISION_MODE 32 not supported for COLLISION_MODEL >=6")
            return
          endif
          call allocate_cmat_fp32(nv,nc_loc_coll,nt1,nt2)
       else
          call allocate_cmat(nv,nc_loc_coll,nt1,nt2)
       endif
    endif

end subroutine cgyro_coll_data_init

subroutine cgyro_coll_data_cleanup(gpu_bigmem_flag)

    implicit none

    integer, intent(in) :: gpu_bigmem_flag

    call deallocate_cmat
    call deallocate_cmat_fp32
    if (allocated(cmat_stripes)) then
#if defined(OMPGPU)
!$omp target exit data map(release:cmat_stripes) if (gpu_bigmem_flag > 0)
#elif defined(_OPENACC)
!$acc exit data delete(cmat_stripes) if (gpu_bigmem_flag > 0)
#endif
       deallocate(cmat_stripes)
    endif
    if (allocated(cmat_simple)) then
#if defined(OMPGPU)
!$omp target exit data map(release:cmat_simple)
#elif defined(_OPENACC)
!$acc exit data delete(cmat_simple)
#endif
       deallocate(cmat_simple)
    endif

end subroutine cgyro_coll_data_cleanup

end module cgyro_coll_data
