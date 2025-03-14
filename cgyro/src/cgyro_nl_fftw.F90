!-----------------------------------------------------------------
! cgyro_nl_fftw.F90
!
! PURPOSE:
!  Evaluate nonlinear bracket with dealiased FFT.  It is natural 
!  to use the FFTW complex-to-real (c2r) transform to move to real 
!  space, compute the convolution (uv), and transform back to the 
!  spectral form using a real-to-complex transform (r2c).
!  
! NOTE: Need to be careful with (p=-nr/2,n=0) component.
!-----------------------------------------------------------------

module cgyro_nl

  implicit none

contains

#if defined(_OPENACC) || defined(OMPGPU)
#define CGYRO_GPU_FFT
#endif

!-----------------------------------------------------------------
! cgyro_nl_fftw_init
!
! Initialize all the relevant plan_r2c and plan_c2r variants
! Many different variants due to many possible FFT backends
!-----------------------------------------------------------------

#ifdef CGYRO_GPU_FFT

#if defined(MKLGPU)
include 'fftw/offload/fftw3_omp_offload.f90'
#endif

subroutine cgyro_nl_fftw_init_fp64

#if defined(HIPGPU)
  use hipfort_hipfft
#elif defined(MKLGPU)
  use fftw3_omp_offload
#else
  use cufft
#endif
  use cgyro_globals

  implicit none
#if defined(MKLGPU)
  include 'fftw/fftw3.f'
#endif

  integer :: howmany,istatus
  integer, parameter :: irank = 2
  integer, dimension(irank) :: ndim,inembed,onembed
  integer :: idist,odist
  integer, parameter :: istride = 1
  integer, parameter :: ostride = 1

#if !defined(MKLGPU)
  integer, parameter :: singlePrecision = selected_real_kind(6,30)
#endif

  !-------------------------------------------------------------------
  ! 2D
  !   input[ b*idist + (x * inembed[1] + y)*istride ]
  !  output[ b*odist + (x * onembed[1] + y)*ostride ]
  !  isign is the sign of the exponent in the formula that defines
  !  Fourier transform  -1 == FFTW_FORWARD
  !                      1 == FFTW_BACKWARD
  !-------------------------------------------------------------------

#if defined(MKLGPU)
  ! oneMKL offload uses the reverse ordering
  ndim(2) = nx
  ndim(1) = ny
#else
  ndim(1) = nx
  ndim(2) = ny
#endif
  idist = size(fxmany,1)*size(fxmany,2)
  odist = size(uxmany,1)*size(uxmany,2)
  inembed = size(fxmany,1)
  onembed = size(uxmany,1)
#if defined(MKLGPU)
  inembed(2) = size(fxmany,2)
  onembed(2) = size(uxmany,2)
#endif

#if defined(HIPGPU)
  plan_c2r_manyA = c_null_ptr
  istatus = hipfftPlanMany(&
       plan_c2r_manyA, &
       irank, &
       ndim, &
       inembed, &
       istride, &
       idist, &
       onembed, &
       ostride, &
       odist, &
       HIPFFT_Z2D, &
       nsplitA)

  if (nsplitB > 0) then ! no fft if nsplitB==0
    plan_c2r_manyB = c_null_ptr
    istatus = hipfftPlanMany(&
       plan_c2r_manyB, &
       irank, &
       ndim, &
       inembed, &
       istride, &
       idist, &
       onembed, &
       ostride, &
       odist, &
       HIPFFT_Z2D, &
       nsplitB)
  endif

  plan_c2r_manyG = c_null_ptr
  istatus = hipfftPlanMany(&
       plan_c2r_manyG, &
       irank, &
       ndim, &
       inembed, &
       istride, &
       idist, &
       onembed, &
       ostride, &
       odist, &
       HIPFFT_Z2D, &
       nsplit)
#elif defined(MKLGPU)
     plan_c2r_manyA = 0
!$omp target data map(tofrom: fymany,uymany)
     !$omp dispatch
     call dfftw_plan_many_dft_c2r(&
          plan_c2r_manyA, &
          irank, &
          ndim, &
          nsplitA, &
          fymany, &
          inembed, &
          istride, &
          idist, &
          uymany, &
          onembed, &
          ostride, &
          odist, &
          FFTW_ESTIMATE)

  if (nsplitB > 0) then ! no fft if nsplitB==0
     plan_c2r_manyB = 0
     !$omp dispatch
     call dfftw_plan_many_dft_c2r(&
          plan_c2r_manyB, &
          irank, &
          ndim, &
          nsplitB, &
          fymany, &
          inembed, &
          istride, &
          idist, &
          uymany, &
          onembed, &
          ostride, &
          odist, &
          FFTW_ESTIMATE)
  endif
!$omp end target data

     plan_c2r_manyG = 0
!$omp target data map(tofrom: gymany,vymany)
     !$omp dispatch
     call dfftw_plan_many_dft_c2r(&
          plan_c2r_manyG, &
          irank, &
          ndim, &
          nsplit, &
          gymany, &
          inembed, &
          istride, &
          idist, &
          vymany, &
          onembed, &
          ostride, &
          odist, &
          FFTW_ESTIMATE)
!$omp end target data

#else
  istatus = cufftPlanMany(&
       plan_c2r_manyA, &
       irank, &
       ndim, &
       inembed, &
       istride, &
       idist, &
       onembed, &
       ostride, &
       odist, &
       CUFFT_Z2D, &
       nsplitA)

  if (nsplitB > 0) then ! no fft if nsplitB==0
    istatus = cufftPlanMany(&
       plan_c2r_manyB, &
       irank, &
       ndim, &
       inembed, &
       istride, &
       idist, &
       onembed, &
       ostride, &
       odist, &
       CUFFT_Z2D, &
       nsplitB)
  endif

  istatus = cufftPlanMany(&
       plan_c2r_manyG, &
       irank, &
       ndim, &
       inembed, &
       istride, &
       idist, &
       onembed, &
       ostride, &
       odist, &
       CUFFT_Z2D, &
       nsplit)
#endif

  idist = size(uxmany,1)*size(uxmany,2)
  odist = size(fxmany,1)*size(fxmany,2)
  inembed = size(uxmany,1)
  onembed = size(fxmany,1) 
#if defined(MKLGPU)
  inembed(2) = size(uxmany,2)
  onembed(2) = size(fxmany,2) 
#endif

#if defined(HIPGPU)
  plan_r2c_manyA = c_null_ptr
  istatus = hipfftPlanMany(&
       plan_r2c_manyA, &
       irank, &
       ndim, &
       inembed, &
       istride, &
       idist, &
       onembed, &
       ostride, &
       odist, &
       HIPFFT_D2Z, &
       nsplitA)

  if (nsplitB > 0) then ! no fft if nsplitB==0
    plan_r2c_manyB = c_null_ptr
    istatus = hipfftPlanMany(&
       plan_r2c_manyB, &
       irank, &
       ndim, &
       inembed, &
       istride, &
       idist, &
       onembed, &
       ostride, &
       odist, &
       HIPFFT_D2Z, &
       nsplitB)
  endif
#elif defined(MKLGPU)
     plan_r2c_manyA = 0
!$omp target data map(tofrom: uvmany,fxmany)
     !$omp dispatch
     call dfftw_plan_many_dft_r2c(&
          plan_r2c_manyA, &
          irank, &
          ndim, &
          nsplitA, &
          uvmany, &
          inembed, &
          istride, &
          idist, &
          fxmany, &
          onembed, &
          ostride, &
          odist, &
          FFTW_ESTIMATE)

  if (nsplitB > 0) then ! no fft if nsplitB==0
     plan_r2c_manyB = 0
     !$omp dispatch
     call dfftw_plan_many_dft_r2c(&
          plan_r2c_manyB, &
          irank, &
          ndim, &
          nsplitB, &
          uvmany, &
          inembed, &
          istride, &
          idist, &
          fxmany, &
          onembed, &
          ostride, &
          odist, &
          FFTW_ESTIMATE)
  endif
!$omp end target data

#else
  istatus = cufftPlanMany(&
       plan_r2c_manyA, &
       irank, &
       ndim, &
       inembed, &
       istride, &
       idist, &
       onembed, &
       ostride, &
       odist, &
       CUFFT_D2Z, &
       nsplitA)

  if (nsplitB > 0) then ! no fft if nsplitB==0
    istatus = cufftPlanMany(&
       plan_r2c_manyB, &
       irank, &
       ndim, &
       inembed, &
       istride, &
       idist, &
       onembed, &
       ostride, &
       odist, &
       CUFFT_D2Z, &
       nsplitB)
  endif
#endif

end subroutine cgyro_nl_fftw_init_fp64

subroutine cgyro_nl_fftw_init_fp32

#if defined(HIPGPU)
  use hipfort_hipfft
#elif defined(MKLGPU)
  use fftw3_omp_offload
#else
  use cufft
#endif
  use cgyro_globals

  implicit none
#if defined(MKLGPU)
  include 'fftw/fftw3.f'
#endif

  integer :: howmany,istatus
  integer, parameter :: irank = 2
  integer, dimension(irank) :: ndim,inembed,onembed
  integer :: idist,odist
  integer, parameter :: istride = 1
  integer, parameter :: ostride = 1

#if !defined(MKLGPU)
  integer, parameter :: singlePrecision = selected_real_kind(6,30)
#endif

  !-------------------------------------------------------------------
  ! 2D
  !   input[ b*idist + (x * inembed[1] + y)*istride ]
  !  output[ b*odist + (x * onembed[1] + y)*ostride ]
  !  isign is the sign of the exponent in the formula that defines
  !  Fourier transform  -1 == FFTW_FORWARD
  !                      1 == FFTW_BACKWARD
  !-------------------------------------------------------------------

#if defined(MKLGPU)
  ! oneMKL offload uses the reverse ordering
  ndim(2) = nx
  ndim(1) = ny
#else
  ndim(1) = nx
  ndim(2) = ny
#endif
  idist = size(fxmany32,1)*size(fxmany32,2)
  odist = size(uxmany32,1)*size(uxmany32,2)
  inembed = size(fxmany32,1)
  onembed = size(uxmany32,1)
#if defined(MKLGPU)
  inembed(2) = size(fxmany32,2)
  onembed(2) = size(uxmany32,2)
#endif

#if defined(HIPGPU)
  plan_c2r_manyA = c_null_ptr
  istatus = hipfftPlanMany(&
       plan_c2r_manyA, &
       irank, &
       ndim, &
       inembed, &
       istride, &
       idist, &
       onembed, &
       ostride, &
       odist, &
       HIPFFT_C2R, &
       nsplitA)

  if (nsplitB > 0) then ! no fft if nsplitB==0
    plan_c2r_manyB = c_null_ptr
    istatus = hipfftPlanMany(&
       plan_c2r_manyB, &
       irank, &
       ndim, &
       inembed, &
       istride, &
       idist, &
       onembed, &
       ostride, &
       odist, &
       HIPFFT_C2R, &
       nsplitB)
  endif

  plan_c2r_manyG = c_null_ptr
  istatus = hipfftPlanMany(&
       plan_c2r_manyG, &
       irank, &
       ndim, &
       inembed, &
       istride, &
       idist, &
       onembed, &
       ostride, &
       odist, &
       HIPFFT_C2R, &
       nsplit)
#elif defined(MKLGPU)
     plan_c2r_manyA = 0
!$omp target data map(tofrom: fymany32,uymany32)
     !$omp dispatch
     call dfftw_plan_many_dft_c2r(&
          plan_c2r_manyA, &
          irank, &
          ndim, &
          nsplitA, &
          fymany32, &
          inembed, &
          istride, &
          idist, &
          uymany32, &
          onembed, &
          ostride, &
          odist, &
          FFTW_ESTIMATE)

  if (nsplitB > 0) then ! no fft if nsplitB==0
     plan_c2r_manyB = 0
     !$omp dispatch
     call dfftw_plan_many_dft_c2r(&
          plan_c2r_manyB, &
          irank, &
          ndim, &
          nsplitB, &
          fymany32, &
          inembed, &
          istride, &
          idist, &
          uymany32, &
          onembed, &
          ostride, &
          odist, &
          FFTW_ESTIMATE)
  endif
!$omp end target data

     plan_c2r_manyG = 0
!$omp target data map(tofrom: gymany32,vymany32)
     !$omp dispatch
     call dfftw_plan_many_dft_c2r(&
          plan_c2r_manyG, &
          irank, &
          ndim, &
          nsplit, &
          gymany32, &
          inembed, &
          istride, &
          idist, &
          vymany32, &
          onembed, &
          ostride, &
          odist, &
          FFTW_ESTIMATE)
!$omp end target data

#else
  istatus = cufftPlanMany(&
       plan_c2r_manyA, &
       irank, &
       ndim, &
       inembed, &
       istride, &
       idist, &
       onembed, &
       ostride, &
       odist, &
       CUFFT_C2R, &
       nsplitA)

  if (nsplitB > 0) then ! no fft if nsplitB==0
    istatus = cufftPlanMany(&
       plan_c2r_manyB, &
       irank, &
       ndim, &
       inembed, &
       istride, &
       idist, &
       onembed, &
       ostride, &
       odist, &
       CUFFT_C2R, &
       nsplitB)
  endif

  istatus = cufftPlanMany(&
       plan_c2r_manyG, &
       irank, &
       ndim, &
       inembed, &
       istride, &
       idist, &
       onembed, &
       ostride, &
       odist, &
       CUFFT_C2R, &
       nsplit)
#endif

  idist = size(uxmany32,1)*size(uxmany32,2)
  odist = size(fxmany32,1)*size(fxmany32,2)
  inembed = size(uxmany32,1)
  onembed = size(fxmany32,1) 
#if defined(MKLGPU)
  inembed(2) = size(uxmany32,2)
  onembed(2) = size(fxmany32,2) 
#endif

#if defined(HIPGPU)
  plan_r2c_manyA = c_null_ptr
  istatus = hipfftPlanMany(&
       plan_r2c_manyA, &
       irank, &
       ndim, &
       inembed, &
       istride, &
       idist, &
       onembed, &
       ostride, &
       odist, &
       HIPFFT_R2C, &
       nsplitA)

  if (nsplitB > 0) then ! no fft if nsplitB==0
    plan_r2c_manyB = c_null_ptr
    istatus = hipfftPlanMany(&
       plan_r2c_manyB, &
       irank, &
       ndim, &
       inembed, &
       istride, &
       idist, &
       onembed, &
       ostride, &
       odist, &
       HIPFFT_R2C, &
       nsplitB)
  endif
#elif defined(MKLGPU)
     plan_r2c_manyA = 0
!$omp target data map(tofrom: uvmany32,fxmany32)
     !$omp dispatch
     call dfftw_plan_many_dft_r2c(&
          plan_r2c_manyA, &
          irank, &
          ndim, &
          nsplitA, &
          uvmany32, &
          inembed, &
          istride, &
          idist, &
          fxmany32, &
          onembed, &
          ostride, &
          odist, &
          FFTW_ESTIMATE)

  if (nsplitB > 0) then ! no fft if nsplitB==0
     plan_r2c_manyB = 0
     !$omp dispatch
     call dfftw_plan_many_dft_r2c(&
          plan_r2c_manyB, &
          irank, &
          ndim, &
          nsplitB, &
          uvmany32, &
          inembed, &
          istride, &
          idist, &
          fxmany32, &
          onembed, &
          ostride, &
          odist, &
          FFTW_ESTIMATE)
  endif
!$omp end target data

#else
  istatus = cufftPlanMany(&
       plan_r2c_manyA, &
       irank, &
       ndim, &
       inembed, &
       istride, &
       idist, &
       onembed, &
       ostride, &
       odist, &
       CUFFT_R2C, &
       nsplitA)

  if (nsplitB > 0) then ! no fft if nsplitB==0
    istatus = cufftPlanMany(&
       plan_r2c_manyB, &
       irank, &
       ndim, &
       inembed, &
       istride, &
       idist, &
       onembed, &
       ostride, &
       odist, &
       CUFFT_R2C, &
       nsplitB)
  endif
#endif

end subroutine cgyro_nl_fftw_init_fp32

#else  /* if not defined CGYRO_GPU_FFT */

subroutine cgyro_nl_fftw_init_fp64
  use cgyro_globals
  implicit none
  include 'fftw3.f03'

  ! Create plans once and for all, with global arrays fx,ux
  plan_c2r = fftw_plan_dft_c2r_2d(nx,ny,fx(:,:,1),uxmany(:,:,1),FFTW_PATIENT)
  plan_r2c = fftw_plan_dft_r2c_2d(nx,ny,uv(:,:,1),fx(:,:,1),FFTW_PATIENT)

end subroutine cgyro_nl_fftw_init_fp64

subroutine cgyro_nl_fftw_init_fp32
  use cgyro_globals
  implicit none
  include 'fftw3.f03'

  ! Create plans once and for all, with global arrays fx,ux
  plan_c2r = fftwf_plan_dft_c2r_2d(nx,ny,fx32(:,:,1),uxmany32(:,:,1),FFTW_PATIENT)
  plan_r2c = fftwf_plan_dft_r2c_2d(nx,ny,uv32(:,:,1),fx32(:,:,1),FFTW_PATIENT)

end subroutine cgyro_nl_fftw_init_fp32

#endif /* CGYRO_GPU_FFT */

!-----------------------------------------------------------------
! cgyro_nl_fftw
! (and all the helper sub-functions)
!
! Do the actual FFT-base NL step
!-----------------------------------------------------------------

#ifdef CGYRO_GPU_FFT

!-----------------------------------------------------------------
!  GPU code for cgyro_nl_fftw and dependent sub-functions
!-----------------------------------------------------------------

!-----------------------------------------------------------------
! Helper vector multiply function with GPU offload
!-----------------------------------------------------------------

subroutine cgyro_nl_fftw_mul(sz,uvm,uxm,vym,uym,vxm,inv_nxny)
  implicit none

  !-----------------------------------
  integer, intent(in) :: sz
  real, dimension(*),intent(out) :: uvm
  real, dimension(*),intent(in) :: uxm,vym,uym,vxm
  real, intent(in) :: inv_nxny
  !-----------------------------------

  integer :: i

#if defined(OMPGPU)
!$omp target teams distribute parallel do simd &
!$omp&  map(to:uxm(1:sz),vym(1:sz),uym(1:sz),vxm(1:sz)) &
!$omp&  map(from:uvm(1:sz))
#else
!$acc parallel loop independent gang vector &
!$acc&         present(uvm,uxm,vym,uym,vxm) private(i)
#endif
  do i=1,sz
    uvm(i) = (uxm(i)*vym(i)-uym(i)*vxm(i))*inv_nxny
  enddo

end subroutine cgyro_nl_fftw_mul

subroutine cgyro_nl_fftw_mul32(sz,uvm,uxm,vym,uym,vxm,inv_nxny)
  use, intrinsic :: iso_fortran_env

  implicit none

  !-----------------------------------
  integer, intent(in) :: sz
  real(KIND=REAL32), dimension(*),intent(out) :: uvm
  real(KIND=REAL32), dimension(*),intent(in) :: uxm,vym,uym,vxm
  real(KIND=REAL32), intent(in) :: inv_nxny
  !-----------------------------------

  integer :: i

#if defined(OMPGPU)
!$omp target teams distribute parallel do simd &
!$omp&  map(to:uxm(1:sz),vym(1:sz),uym(1:sz),vxm(1:sz)) &
!$omp&  map(from:uvm(1:sz))
#else
!$acc parallel loop independent gang vector &
!$acc&         present(uvm,uxm,vym,uym,vxm) private(i)
#endif
  do i=1,sz
    uvm(i) = (uxm(i)*vym(i)-uym(i)*vxm(i))*inv_nxny
  enddo

end subroutine cgyro_nl_fftw_mul32

subroutine cgyro_nl_fftw_mul_sub_mean(ny,nx,nt,uvm,uxm,vym,uym,vxm,inv_nxny)
  implicit none

  integer, intent(in) :: ny,nx,nt
  real, dimension(*),intent(out) :: uvm
  real, dimension(*),intent(in) :: uxm,vym,uym,vxm
  real, intent(in) :: inv_nxny

  integer :: iy,ix,it
  integer :: sz,sz_t
  integer :: i,i_xt
  real :: y_mean_ux, y_mean_uy, y_mean_vx, y_mean_vy
  real :: r_ux, r_uy, r_vx, r_vy
  real :: inv_ny

  sz = ny*nx*nt
  sz_t = ny*nx
  inv_ny = 1.0/ny

#if defined(OMPGPU)
!$omp target teams distribute parallel do collapse(2) &
!$omp&  private(iy,i_xt,y_mean_ux,y_mean_uy,y_mean_vx,y_mean_vy) &
!$omp&  private(i,r_ux,r_uy,r_vx,r_vy) &
!$omp&  map(to:uxm(1:sz),vym(1:sz),uym(1:sz),vxm(1:sz)) &
!$omp&  map(from:uvm(1:sz))
#else
!$acc parallel loop independent gang collapse(2) &
!$acc&         private(iy,i_xt,y_mean_ux,y_mean_uy,y_mean_vx,y_mean_vy) &
!$acc&         present(uvm,uxm,vym,uym,vxm)
#endif
  do it=0,nt-1
   do ix=0,nx-1
    i_xt = (it*sz_t) + ix*ny;
    y_mean_ux = sum(uxm(i_xt+1:i_xt+ny) ) * inv_ny
    y_mean_uy = sum(uym(i_xt+1:i_xt+ny) ) * inv_ny
    y_mean_vx = sum(vxm(i_xt+1:i_xt+ny) ) * inv_ny
    y_mean_vy = sum(vym(i_xt+1:i_xt+ny) ) * inv_ny
#if (!defined(OMPGPU)) && defined(_OPENACC)
!$acc loop vector private(i,r_ux,r_uy,r_vx,r_vy)
#endif
    do iy=1,ny
      i = i_xt+iy
      ! remove ky=0
      r_ux = uxm(i) - y_mean_ux
      r_uy = uym(i) - y_mean_uy
      r_vx = vxm(i) - y_mean_vx
      r_vy = vym(i) - y_mean_vy

      ! compute and save uv
      uvm(i) = (r_ux * r_vy - r_uy * r_vx)*inv_nxny
    enddo
   enddo
  enddo

end subroutine

!-----------------------------------------------------------------
! Batch FFT invocations
! All the c2r and r2c variants, both fp64 and fp32
!
! Many different internal variants due to many possible FFT backends
!-----------------------------------------------------------------

subroutine cgyro_fft_z2d(plan, indata, outdata)

#if defined(HIPGPU)
  use hipfort
  use hipfort_hipfft
#elif defined(MKLGPU)
  use fftw3_omp_offload
#else
  use cufft
#endif
  use cgyro_nl_comm

  implicit none
#if defined(MKLGPU)
  include 'fftw/fftw3.f'
#endif
  !-----------------------------------
#if defined(HIPGPU)
  type(C_PTR), intent(inout)    :: plan
#elif defined(MKLGPU)
  INTEGER*8, intent(inout)      :: plan
#else
  integer(c_int), intent(inout) :: plan
#endif
  complex, dimension(:,:,:), intent(inout) :: indata
  real, dimension(:,:,:), intent(inout)    :: outdata
  !-----------------------------------
  integer :: rc

#if defined(OMPGPU)

#if defined(MKLGPU)
!$omp target data map(tofrom: indata, outdata)
#else
!$omp target data use_device_ptr(indata, outdata)
#endif

#else
!$acc  host_data use_device(indata, outdata)
#endif

#if defined(HIPGPU)
  rc = hipfftExecZ2D(plan,c_loc(indata),c_loc(outdata))
#elif defined(MKLGPU)
  !$omp dispatch
  call dfftw_execute_dft_c2r(plan,indata, outdata)
  rc = 0
#else
  rc = cufftExecZ2D(plan,indata, outdata)
#endif

#ifdef HIPGPU
  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()
  ! hipfftExec is asynchronous, will need the results below
  rc = hipDeviceSynchronize()
#endif

#if defined(OMPGPU)
!$omp end target data
#else
  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()
!$acc wait
!$acc end host_data
#endif

end subroutine cgyro_fft_z2d

subroutine cgyro_fft_c2r(plan, indata, outdata)

#if defined(HIPGPU)
  use hipfort
  use hipfort_hipfft
#elif defined(MKLGPU)
  use fftw3_omp_offload
#else
  use cufft
#endif
  use cgyro_nl_comm
  use, intrinsic :: iso_fortran_env

  implicit none
#if defined(MKLGPU)
  include 'fftw/fftw3.f'
#endif
  !-----------------------------------
#if defined(HIPGPU)
  type(C_PTR), intent(inout)    :: plan
#elif defined(MKLGPU)
  INTEGER*8, intent(inout)      :: plan
#else
  integer(c_int), intent(inout) :: plan
#endif
  complex(KIND=REAL32), dimension(:,:,:), intent(inout) :: indata
  real(KIND=REAL32), dimension(:,:,:), intent(inout)    :: outdata
  !-----------------------------------
  integer :: rc

#if defined(OMPGPU)

#if defined(MKLGPU)
!$omp target data map(tofrom: indata, outdata)
#else
!$omp target data use_device_ptr(indata, outdata)
#endif

#else
!$acc  host_data use_device(indata, outdata)
#endif

#if defined(HIPGPU)
  rc = hipfftExecC2R(plan,c_loc(indata),c_loc(outdata))
#elif defined(MKLGPU)
  !$omp dispatch
  call dfftw_execute_dft_c2r(plan,indata, outdata)
  rc = 0
#else
  rc = cufftExecC2R(plan,indata, outdata)
#endif

#ifdef HIPGPU
  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()
  ! hipfftExec is asynchronous, will need the results below
  rc = hipDeviceSynchronize()
#endif

#if defined(OMPGPU)
!$omp end target data
#else
  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()
!$acc wait
!$acc end host_data
#endif

end subroutine cgyro_fft_c2r

subroutine cgyro_fft_d2z(plan, indata, outdata)

#if defined(HIPGPU)
  use hipfort
  use hipfort_hipfft
#elif defined(MKLGPU)
  use fftw3_omp_offload
#else
  use cufft
#endif
  use cgyro_nl_comm

  implicit none
#if defined(MKLGPU)
  include 'fftw/fftw3.f'
#endif
  !-----------------------------------
#if defined(HIPGPU)
  type(C_PTR), intent(inout)    :: plan
#elif defined(MKLGPU)
  INTEGER*8, intent(inout)      :: plan
#else
  integer(c_int), intent(inout) :: plan
#endif
  real, dimension(:,:,:), intent(inout)    :: indata
  complex, dimension(:,:,:), intent(inout) :: outdata
  !-----------------------------------
  integer :: rc

#if defined(OMPGPU)

#if defined(MKLGPU)
!$omp target data map(tofrom: indata, outdata)
#else
!$omp target data use_device_ptr(indata, outdata)
#endif

#else
!$acc  host_data use_device(indata, outdata)
#endif

#if defined(HIPGPU)
  rc = hipfftExecD2Z(plan,c_loc(indata),c_loc(outdata))
#elif defined(MKLGPU)
  !$omp dispatch
  call dfftw_execute_dft_r2c(plan,indata, outdata)
  rc = 0
#else
  rc = cufftExecD2Z(plan,indata, outdata)
#endif

#ifdef HIPGPU
  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()
  ! hipfftExec is asynchronous, will need the results below
  rc = hipDeviceSynchronize()
#endif

#if defined(OMPGPU)
!$omp end target data
#else
  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()
!$acc wait
!$acc end host_data
#endif

end subroutine cgyro_fft_d2z

subroutine cgyro_fft_r2c(plan, indata, outdata)

#if defined(HIPGPU)
  use hipfort
  use hipfort_hipfft
#elif defined(MKLGPU)
  use fftw3_omp_offload
#else
  use cufft
#endif
  use cgyro_nl_comm
  use, intrinsic :: iso_fortran_env

  implicit none
#if defined(MKLGPU)
  include 'fftw/fftw3.f'
#endif
  !-----------------------------------
#if defined(HIPGPU)
  type(C_PTR), intent(inout)    :: plan
#elif defined(MKLGPU)
  INTEGER*8, intent(inout)      :: plan
#else
  integer(c_int), intent(inout) :: plan
#endif
  real(KIND=REAL32), dimension(:,:,:), intent(inout)    :: indata
  complex(KIND=REAL32), dimension(:,:,:), intent(inout) :: outdata
  !-----------------------------------
  integer :: rc

#if defined(OMPGPU)

#if defined(MKLGPU)
!$omp target data map(tofrom: indata, outdata)
#else
!$omp target data use_device_ptr(indata, outdata)
#endif

#else
!$acc  host_data use_device(indata, outdata)
#endif

#if defined(HIPGPU)
  rc = hipfftExecR2C(plan,c_loc(indata),c_loc(outdata))
#elif defined(MKLGPU)
  !$omp dispatch
  call dfftw_execute_dft_r2c(plan,indata, outdata)
  rc = 0
#else
  rc = cufftExecR2C(plan,indata, outdata)
#endif

#ifdef HIPGPU
  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()
  ! hipfftExec is asynchronous, will need the results below
  rc = hipDeviceSynchronize()
#endif

#if defined(OMPGPU)
!$omp end target data
#else
  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()
!$acc wait
!$acc end host_data
#endif

end subroutine cgyro_fft_r2c

!-----------------------------------------------------------------
! Helper function for zero-ing off-diagonal elements
! with proper GPU offload
! Using async compute whenever possible
! Caller responsible for synchronization
!-----------------------------------------------------------------

subroutine cgyro_zero_offdiag_async(nj,xmany,ymany)

  use cgyro_globals

  implicit none
  !-----------------------------------
  integer, intent(in) :: nj
  complex, dimension(0:ny2,0:nx-1,nj), intent(inout) :: xmany,ymany
  !-----------------------------------
  integer :: j
  integer :: ix,iy

#if defined(OMPGPU)
  !no async for OMPGPU for now
!$omp target teams distribute parallel do simd collapse(3)
#else
!$acc parallel loop gang vector independent collapse(3) async(2) &
!$acc&         private(j,ix,iy) &
!$acc&         present(xmany,ymany)
#endif
  do j=1,nj
     do ix=nx2,nx0-1
       do iy=0,ny2
         xmany(iy,ix,j) = 0
         ymany(iy,ix,j) = 0
       enddo
     enddo
  enddo

#if defined(OMPGPU)
  !no async for OMPGPU for now
!$omp target teams distribute parallel do simd collapse(3)
#else
!$acc parallel loop gang vector independent collapse(3) async(2) &
!$acc&         private(j,ix,iy) &
!$acc&         present(xmany,ymany)
#endif
  do j=1,nj
     do ix=0,nx-1
       do iy=n_toroidal,ny2
         xmany(iy,ix,j) = 0
         ymany(iy,ix,j) = 0
       enddo
     enddo
  enddo

end subroutine cgyro_zero_offdiag_async

subroutine cgyro_zero_offdiag32_async(nj,xmany,ymany)

  use cgyro_globals

  implicit none
  !-----------------------------------
  integer, intent(in) :: nj
  complex(KIND=REAL32), dimension(0:ny2,0:nx-1,nj), intent(inout) :: xmany,ymany
  !-----------------------------------
  integer :: j
  integer :: ix,iy

#if defined(OMPGPU)
  !no async for OMPGPU for now
!$omp target teams distribute parallel do simd collapse(3)
#else
!$acc parallel loop gang vector independent collapse(3) async(2) &
!$acc&         private(j,ix,iy) &
!$acc&         present(xmany,ymany)
#endif
  do j=1,nj
     do ix=nx2,nx0-1
       do iy=0,ny2
         xmany(iy,ix,j) = 0
         ymany(iy,ix,j) = 0
       enddo
     enddo
  enddo

#if defined(OMPGPU)
  !no async for OMPGPU for now
!$omp target teams distribute parallel do simd collapse(3)
#else
!$acc parallel loop gang vector independent collapse(3) async(2) &
!$acc&         private(j,ix,iy) &
!$acc&         present(xmany,ymany)
#endif
  do j=1,nj
     do ix=0,nx-1
       do iy=n_toroidal,ny2
         xmany(iy,ix,j) = 0
         ymany(iy,ix,j) = 0
       enddo
     enddo
  enddo

end subroutine cgyro_zero_offdiag32_async

!-----------------------------------------------------------------
! Helper function for filling fxmany and fymany (inputs to FFT)
! with proper GPU offload.
! Input either fA_nl or fB_nl, with nj being either nsplitA or nplitB
! Using async compute whenever possible
! Caller responsible for synchronization
!-----------------------------------------------------------------

subroutine cgyro_fmany_async(nj, f_nl)

  use cgyro_globals

  implicit none
  !-----------------------------------
  integer, intent(in) :: nj
  complex, dimension(n_radial,nt_loc,nj,n_toroidal_procs), intent(in) :: f_nl
  !-----------------------------------
  integer :: j,p
  integer :: ir,itm,itl,ix,iy
  integer :: itor, it_loc
  complex :: f0

! f_nl is (radial, nt_loc, theta, nv_loc1, toroidal_procs)
! where nv_loc1 * toroidal_procs >= nv_loc

  ! no tiling, does not seem to help
#if defined(OMPGPU)
  !no async for OMPGPU for now
!$omp target teams distribute parallel do simd collapse(4) &
!$omp&  private(j,ir,p,ix,itor,iy,f0,itm,itl)
#else
!$acc parallel loop gang vector independent collapse(4) async(2) &
!$acc&         private(j,ir,p,ix,itor,iy,f0,itm,itl) present(f_nl,fxmany,fymany)
#endif
  do j=1,nj
     do ir=1,n_radial
       do itm=1,n_toroidal_procs
         do itl=1,nt_loc
              itor=itl + (itm-1)*nt_loc
              iy = itor-1
              p  = ir-1-nx0/2
              ix = p
              if (ix < 0) ix = ix+nx

              f0 = i_c*f_nl(ir,itl,j,itm)
              fxmany(iy,ix,j) = p*f0
              fymany(iy,ix,j) = iy*f0
         enddo
       enddo
     enddo
  enddo

end subroutine cgyro_fmany_async

subroutine cgyro_fmany32_async(nj, f_nl)

  use cgyro_globals

  implicit none
  !-----------------------------------
  integer, intent(in) :: nj
  complex(KIND=REAL32), dimension(n_radial,nt_loc,nj,n_toroidal_procs), intent(in) :: f_nl
  !-----------------------------------
  integer :: j,p
  integer :: ir,itm,itl,ix,iy
  integer :: itor, it_loc
  complex(KIND=REAL32) :: f0

! f_nl is (radial, nt_loc, theta, nv_loc1, toroidal_procs)
! where nv_loc1 * toroidal_procs >= nv_loc

  ! no tiling, does not seem to help
#if defined(OMPGPU)
  !no async for OMPGPU for now
!$omp target teams distribute parallel do simd collapse(4) &
!$omp&  private(j,ir,p,ix,itor,iy,f0,itm,itl)
#else
!$acc parallel loop gang vector independent collapse(4) async(2) &
!$acc&         private(j,ir,p,ix,itor,iy,f0,itm,itl) present(f_nl,fxmany,fymany)
#endif
  do j=1,nj
     do ir=1,n_radial
       do itm=1,n_toroidal_procs
         do itl=1,nt_loc
              itor=itl + (itm-1)*nt_loc
              iy = itor-1
              p  = ir-1-nx0/2
              ix = p
              if (ix < 0) ix = ix+nx

              f0 = i_c*f_nl(ir,itl,j,itm)
              fxmany32(iy,ix,j) = p*f0
              fymany32(iy,ix,j) = iy*f0
         enddo
       enddo
     enddo
  enddo

end subroutine cgyro_fmany32_async

!-----------------------------------------------------------------
! Helper function for filling gxmany and gymany (inputs to FFT)
! with proper GPU offload.
! No arguments, uses fixed logic and cgyro_globals.
! Using async compute whenever possible
! Caller responsible for synchronization
!-----------------------------------------------------------------

subroutine cgyro_gmany_async

  use cgyro_globals

  implicit none
  !-----------------------------------
  integer :: j,p
  integer :: it, ir,itm,itl,ix,iy
  integer :: mytm, itor, it_loc
  integer :: jtheta_min
  integer :: iy0, iy1, ir0, ir1
  complex :: g0

#ifdef GACODE_GPU_AMD
  ! AMD GPU  (MI250X) optimal
  integer, parameter :: F_RADTILE = 8
  integer, parameter :: F_TORTILE = 16
#else
  ! NVIDIA GPU  (A100) optimal
  integer, parameter :: F_RADTILE = 8
  integer, parameter :: F_TORTILE = 16
#endif

! g_nl      is (n_field,n_radial,n_jtheta,nt_loc,n_toroidal_procs)
! jcev_c_nl is (n_field,n_radial,n_jtheta,nv_loc,nt_loc,n_toroidal_procs)

  ! tile for performance, since this is effectively a transpose
#if defined(OMPGPU)
  !no async for OMPGPU for now
!$omp target teams distribute parallel do simd collapse(5) &
!$omp&   private(j,p,ix,itor,mytm,iy,g0,it,iv_loc,it_loc,jtheta_min,itm,itl,ir)
#else
!$acc parallel loop gang vector independent collapse(5) async(2) &
!$acc&         private(j,p,ix,itor,mytm,iy,g0,it,iv_loc,it_loc,jtheta_min,itm,itl,ir) &
!$acc&         present(g_nl,jvec_c_nl,gxmany,gymany) &
!$acc&         present(nsplit,n_radial,n_toroidal_procs,nt_loc,nt1,n_theta,nv_loc,nx0)
#endif
  do j=1,nsplit
    do iy0=0,n_toroidal+(F_TORTILE-1)-1,F_TORTILE  ! round up
      do ir0=0,n_radial+(F_RADTILE-1)-1,F_RADTILE  ! round up
        do iy1=0,(F_TORTILE-1)   ! tile
          do ir1=0,(F_RADTILE-1)  ! tile
            iy = iy0+iy1
            ir = 1 + ir0+ir1
            if ((iy < n_toroidal) .and. (ir <= n_radial)) then
              itor = iy+1
              itm = 1 + iy/nt_loc
              itl = 1 + modulo(iy,nt_loc)
              mytm = 1 + nt1/nt_loc !my toroidal proc number
              it = 1+((mytm-1)*nsplit+j-1)/nv_loc
              iv_loc = 1+modulo((mytm-1)*nsplit+j-1,nv_loc)
              jtheta_min = 1+((mytm-1)*nsplit)/nv_loc
              it_loc = it-jtheta_min+1
              iy = itor-1

              p  = ir-1-nx0/2
              ix = p
              if (ix < 0) ix = ix+nx

              if (it > n_theta) then
                 g0 = (0.0,0.0)
              else
                 g0 = i_c*sum( jvec_c_nl(1:n_field,ir,it_loc,iv_loc,itor)*g_nl(1:n_field,ir,it_loc,itor))
              endif
              gxmany(iy,ix,j) = p*g0
              gymany(iy,ix,j) = iy*g0
            endif
          enddo
        enddo
      enddo
    enddo
  enddo
end subroutine cgyro_gmany_async

subroutine cgyro_gmany32_async

  use cgyro_globals

  implicit none
  !-----------------------------------
  integer :: j,p
  integer :: it, ir,itm,itl,ix,iy
  integer :: mytm, itor, it_loc
  integer :: jtheta_min
  integer :: iy0, iy1, ir0, ir1
  complex(KIND=REAL32) :: g0

#ifdef GACODE_GPU_AMD
  ! AMD GPU  (MI250X) optimal
  integer, parameter :: F_RADTILE = 8
  integer, parameter :: F_TORTILE = 16
#else
  ! NVIDIA GPU  (A100) optimal
  integer, parameter :: F_RADTILE = 8
  integer, parameter :: F_TORTILE = 16
#endif

! g_nl      is (n_field,n_radial,n_jtheta,nt_loc,n_toroidal_procs)
! jcev_c_nl is (n_field,n_radial,n_jtheta,nv_loc,nt_loc,n_toroidal_procs)

  ! tile for performance, since this is effectively a transpose
#if defined(OMPGPU)
  !no async for OMPGPU for now
!$omp target teams distribute parallel do simd collapse(5) &
!$omp&   private(j,p,ix,itor,mytm,iy,g0,it,iv_loc,it_loc,jtheta_min,itm,itl,ir)
#else
!$acc parallel loop gang vector independent collapse(5) async(2) &
!$acc&         private(j,p,ix,itor,mytm,iy,g0,it,iv_loc,it_loc,jtheta_min,itm,itl,ir) &
!$acc&         present(g_nl,jvec_c_nl,gxmany32,gymany32) &
!$acc&         present(nsplit,n_radial,n_toroidal_procs,nt_loc,nt1,n_theta,nv_loc,nx0)
#endif
  do j=1,nsplit
    do iy0=0,n_toroidal+(F_TORTILE-1)-1,F_TORTILE  ! round up
      do ir0=0,n_radial+(F_RADTILE-1)-1,F_RADTILE  ! round up
        do iy1=0,(F_TORTILE-1)   ! tile
          do ir1=0,(F_RADTILE-1)  ! tile
            iy = iy0+iy1
            ir = 1 + ir0+ir1
            if ((iy < n_toroidal) .and. (ir <= n_radial)) then
              itor = iy+1
              itm = 1 + iy/nt_loc
              itl = 1 + modulo(iy,nt_loc)
              mytm = 1 + nt1/nt_loc !my toroidal proc number
              it = 1+((mytm-1)*nsplit+j-1)/nv_loc
              iv_loc = 1+modulo((mytm-1)*nsplit+j-1,nv_loc)
              jtheta_min = 1+((mytm-1)*nsplit)/nv_loc
              it_loc = it-jtheta_min+1
              iy = itor-1

              p  = ir-1-nx0/2
              ix = p
              if (ix < 0) ix = ix+nx

              if (it > n_theta) then
                 g0 = (0.0,0.0)
              else
                 g0 = i_c*sum( jvec_c_nl32(1:n_field,ir,it_loc,iv_loc,itor)*g_nl32(1:n_field,ir,it_loc,itor))
              endif
              gxmany32(iy,ix,j) = p*g0
              gymany32(iy,ix,j) = iy*g0
            endif
          enddo
        enddo
      enddo
    enddo
  enddo
end subroutine cgyro_gmany32_async

!-----------------------------------------------------------------
! Helper function for averaging elements to enforce FFT preconditions
! with proper GPU offload.
! Works for both fxmany and gxmany, actual size passed as nj
! Using async compute whenever possible
! Caller responsible for synchronization
!-----------------------------------------------------------------

subroutine cgyro_sym_async(nj, many)

  use cgyro_globals

  implicit none
  !-----------------------------------
  integer, intent(in) :: nj
  complex, dimension(0:ny2,0:nx-1,nj), intent(inout) :: many
  !-----------------------------------
  integer:: j, ix
  complex :: f0

  ! Average elements so as to ensure
  !   f(kx,ky=0) = f(-kx,ky=0)^*
  ! This symmetry is required for complex input to c2r
#if defined(OMPGPU)
  !no async for OMPGPU for now
!$omp target teams distribute parallel do simd collapse(2) &
!$omp&  private(j,ix,f0) firstprivate(nj,nx) 
#else
!$acc parallel loop gang vector independent collapse(2) async(2) &
!$acc&         private(j,ix,f0) firstprivate(nj,nx) &
!$acc&         present(many)
#endif
  do j=1,nj
    do ix=1,nx/2-1
      f0 = 0.5*( many(0,ix,j)+conjg(many(0,nx-ix,j)) )
      many(0,ix   ,j) = f0
      many(0,nx-ix,j) = conjg(f0)
    enddo
  enddo

end subroutine cgyro_sym_async

subroutine cgyro_sym32_async(nj, many)

  use cgyro_globals

  implicit none
  !-----------------------------------
  integer, intent(in) :: nj
  complex(KIND=REAL32), dimension(0:ny2,0:nx-1,nj), intent(inout) :: many
  !-----------------------------------
  integer:: j, ix
  complex :: f0

  ! Average elements so as to ensure
  !   f(kx,ky=0) = f(-kx,ky=0)^*
  ! This symmetry is required for complex input to c2r
#if defined(OMPGPU)
  !no async for OMPGPU for now
!$omp target teams distribute parallel do simd collapse(2) &
!$omp&  private(j,ix,f0) firstprivate(nj,nx) 
#else
!$acc parallel loop gang vector independent collapse(2) async(2) &
!$acc&         private(j,ix,f0) firstprivate(nj,nx) &
!$acc&         present(many)
#endif
  do j=1,nj
    do ix=1,nx/2-1
      f0 = 0.5*( many(0,ix,j)+conjg(many(0,nx-ix,j)) )
      many(0,ix   ,j) = f0
      many(0,nx-ix,j) = conjg(f0)
    enddo
  enddo

end subroutine cgyro_sym32_async

!-----------------------------------------------------------------
! Helper function for transposing FFT output into MPI-friendly format
! with proper GPU offload.
! Acn be used to fill either fA_nl or fB_nl, with nl either nspitA or nsplitB
! Using async compute whenever possible
! Caller responsible for synchronization
!-----------------------------------------------------------------

subroutine cgyro_fmany_r_async(nj, fmany, f_nl)

  use cgyro_globals

  implicit none
  !-----------------------------------
  integer, intent(in) :: nj
  complex, dimension(0:ny2,0:nx-1,nj), intent(in) :: fmany
  complex, dimension(n_radial,nt_loc,nj,n_toroidal_procs), intent(inout) :: f_nl
  !-----------------------------------
  integer :: j,p
  integer :: ir,itm,itl,ix,iy
  integer :: iy0, iy1, ir0, ir1

#ifdef GACODE_GPU_AMD
  ! AMD GPU  (MI250X) optimal
  integer, parameter :: R_RADTILE = 32
  integer, parameter :: R_TORTILE = 8
#else
  ! NVIDIA GPU  (A100) optimal
  integer, parameter :: R_RADTILE = 16
  integer, parameter :: R_TORTILE = 8
#endif

#if defined(OMPGPU)
!$omp target teams distribute parallel do collapse(5) &
!$omp&   private(iy,ir,itm,itl,ix)
#else
!$acc parallel loop independent collapse(5) gang &
!$acc&         private(iy,ir,itm,itl,ix) present(f_nl,fmany)
#endif
  do j=1,nj
   do iy0=0,n_toroidal+(R_TORTILE-1)-1,R_TORTILE  ! round up
    do ir0=0,n_radial+(R_RADTILE-1)-1,R_RADTILE  ! round up
    do iy1=0,(R_TORTILE-1)   ! tile
      do ir1=0,(R_RADTILE-1)  ! tile
       iy = iy0 + iy1
       ir = 1 + ir0 + ir1
       if ((iy < n_toroidal) .and. (ir <= n_radial)) then
           ! itor = iy+1
           itm = 1 + iy/nt_loc
           itl = 1 + modulo(iy,nt_loc)
           ix = ir-1-nx0/2
           if (ix < 0) ix = ix+nx

           f_nl(ir,itl,j,itm) = fmany(iy,ix,j)
        endif
      enddo
     enddo
    enddo
   enddo
  enddo

end subroutine cgyro_fmany_r_async

 subroutine cgyro_fmany_r32_async(nj, fmany, f_nl)

  use cgyro_globals

  implicit none
  !-----------------------------------
  integer, intent(in) :: nj
  complex, dimension(0:ny2,0:nx-1,nj), intent(in) :: fmany
  complex(KIND=REAL32), dimension(n_radial,nt_loc,nj,n_toroidal_procs), intent(inout) :: f_nl
  !-----------------------------------
  integer :: j,p
  integer :: ir,itm,itl,ix,iy
  integer :: iy0, iy1, ir0, ir1

#ifdef GACODE_GPU_AMD
  ! AMD GPU  (MI250X) optimal
  integer, parameter :: R_RADTILE = 32
  integer, parameter :: R_TORTILE = 8
#else
  ! NVIDIA GPU  (A100) optimal
  integer, parameter :: R_RADTILE = 16
  integer, parameter :: R_TORTILE = 8
#endif

#if defined(OMPGPU)
!$omp target teams distribute parallel do collapse(5) &
!$omp&   private(iy,ir,itm,itl,ix)
#else
!$acc parallel loop independent collapse(5) gang &
!$acc&         private(iy,ir,itm,itl,ix) present(f_nl,fmany)
#endif
  do j=1,nj
   do iy0=0,n_toroidal+(R_TORTILE-1)-1,R_TORTILE  ! round up
    do ir0=0,n_radial+(R_RADTILE-1)-1,R_RADTILE  ! round up
    do iy1=0,(R_TORTILE-1)   ! tile
      do ir1=0,(R_RADTILE-1)  ! tile
       iy = iy0 + iy1
       ir = 1 + ir0 + ir1
       if ((iy < n_toroidal) .and. (ir <= n_radial)) then
           ! itor = iy+1
           itm = 1 + iy/nt_loc
           itl = 1 + modulo(iy,nt_loc)
           ix = ir-1-nx0/2
           if (ix < 0) ix = ix+nx

           f_nl(ir,itl,j,itm) = fmany(iy,ix,j)
        endif
      enddo
     enddo
    enddo
   enddo
  enddo

end subroutine cgyro_fmany_r32_async

subroutine cgyro_fmany32_r32_async(nj, fmany, f_nl)

  use cgyro_globals

  implicit none
  !-----------------------------------
  integer, intent(in) :: nj
  complex(KIND=REAL32), dimension(0:ny2,0:nx-1,nj), intent(in) :: fmany
  complex(KIND=REAL32), dimension(n_radial,nt_loc,nj,n_toroidal_procs), intent(inout) :: f_nl
  !-----------------------------------
  integer :: j,p
  integer :: ir,itm,itl,ix,iy
  integer :: iy0, iy1, ir0, ir1

#ifdef GACODE_GPU_AMD
  ! AMD GPU  (MI250X) optimal
  integer, parameter :: R_RADTILE = 32
  integer, parameter :: R_TORTILE = 8
#else
  ! NVIDIA GPU  (A100) optimal
  integer, parameter :: R_RADTILE = 16
  integer, parameter :: R_TORTILE = 8
#endif

#if defined(OMPGPU)
!$omp target teams distribute parallel do collapse(5) &
!$omp&   private(iy,ir,itm,itl,ix)
#else
!$acc parallel loop independent collapse(5) gang &
!$acc&         private(iy,ir,itm,itl,ix) present(f_nl,fmany)
#endif
  do j=1,nj
   do iy0=0,n_toroidal+(R_TORTILE-1)-1,R_TORTILE  ! round up
    do ir0=0,n_radial+(R_RADTILE-1)-1,R_RADTILE  ! round up
    do iy1=0,(R_TORTILE-1)   ! tile
      do ir1=0,(R_RADTILE-1)  ! tile
       iy = iy0 + iy1
       ir = 1 + ir0 + ir1
       if ((iy < n_toroidal) .and. (ir <= n_radial)) then
           ! itor = iy+1
           itm = 1 + iy/nt_loc
           itl = 1 + modulo(iy,nt_loc)
           ix = ir-1-nx0/2
           if (ix < 0) ix = ix+nx

           f_nl(ir,itl,j,itm) = fmany(iy,ix,j)
        endif
      enddo
     enddo
    enddo
   enddo
  enddo

end subroutine cgyro_fmany32_r32_async

!-----------------------------------------------------------------
! High-level cgyro_nl_fftw logic
! Do the actual FFT-base NL step
! NOTE: Expects cgyro_nl_fftw_comm1 to have been already called
!-----------------------------------------------------------------

subroutine cgyro_nl_fftw_fp64(i_triad)

  use timer_lib
  use parallel_lib
  use cgyro_nl_comm
  use cgyro_globals

  implicit none
  !-----------------------------------
  integer, intent(in) :: i_triad
  !-----------------------------------

  real :: inv_nxny
  inv_nxny = 1.0/(nx*ny)

  call timer_lib_in('nl_mem')
  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()

  ! we can zero the elements we know are zero while we wait
  call cgyro_zero_offdiag_async(nsplitA,fxmany,fymany)

  call timer_lib_out('nl_mem')

  ! time to wait for the FA_nl to become avaialble
  call timer_lib_in('nl_comm')
  call parallel_slib_f_nc_wait(nsplitA,fpackA,fA_nl,fA_req)
  fA_req_valid = .FALSE.
  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()
  call timer_lib_out('nl_comm')

  call timer_lib_in('nl')

! f_nl is (radial, nt_loc, theta, nv_loc1, toroidal_procs)
! where nv_loc1 * toroidal_procs >= nv_loc

  call cgyro_fmany_async(nsplitA, fA_nl)

#if defined(OMPGPU)
  !no async for OMPGPU for now
#else
  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()
!$acc wait(2)
#endif

  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()

  ! Average elements so as to ensure
  !   f(kx,ky=0) = f(-kx,ky=0)^*
  ! This symmetry is required for complex input to c2r
  call cgyro_sym_async(nsplitA,fxmany)

     ! --------------------------------------
     ! perform many Fourier Transforms at once
     ! --------------------------------------

  call cgyro_fft_z2d(plan_c2r_manyA,fymany,uymany)
  ! fxmany is complete now

  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()

#if !defined(OMPGPU)
!$acc wait
#endif

  call cgyro_fft_z2d(plan_c2r_manyA,fxmany,uxmany)

  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()

  ! we can zero the elements we know are zero while we wait for comm
  call cgyro_zero_offdiag_async(nsplit,gxmany,gymany)

  call timer_lib_out('nl')

  call timer_lib_in('nl_comm')
  ! time to wait for the g_nl to become avaialble
  call parallel_slib_f_fd_wait(n_field,n_radial,n_jtheta,gpack,g_nl,g_req)
  g_req_valid = .FALSE.
  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()
  call timer_lib_out('nl_comm')

  call timer_lib_in('nl')

! g_nl      is (n_field,n_radial,n_jtheta,nt_loc,n_toroidal_procs)
! jcev_c_nl is (n_field,n_radial,n_jtheta,nv_loc,nt_loc,n_toroidal_procs)

  call cgyro_gmany_async

#if defined(OMPGPU)
  !no async for OMPGPU for now
#else
  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()
!$acc wait(2)
#endif

  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()

  ! Average elements so as to ensure
  !   g(kx,ky=0) = g(-kx,ky=0)^*
  ! This symmetry is required for complex input to c2r
  call cgyro_sym_async(nsplit,gxmany)

  call cgyro_fft_z2d(plan_c2r_manyG,gymany,vymany)

  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()
  ! gxmany is complete now

  call cgyro_fft_z2d(plan_c2r_manyG,gxmany,vxmany)

  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()

  ! Poisson bracket in real space
  ! uv = (ux*vy-uy*vx)/(nx*ny)

  call cgyro_nl_fftw_mul(size(uvmany,1)*size(uvmany,2)*nsplitA, &
                         uvmany, &
                         uxmany,vymany(:,:,1:nsplitA), &
                         uymany,vxmany(:,:,1:nsplitA), &
                         inv_nxny)

  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()

  ! ------------------
  ! Transform uv to fx
  ! ------------------

#if !defined(OMPGPU)
!$acc wait
#endif

  call cgyro_fft_d2z(plan_r2c_manyA,uvmany,fxmany)

  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()

  call timer_lib_out('nl')
  call timer_lib_in('nl_mem')

  ! NOTE: The FFT will generate an unwanted n=0,p=-nr/2 component
  ! that will be filtered in the main time-stepping loop

  ! tile for performance, since this is effectively a transpose
  if (nl_single_flag .EQ. 0) then
    call cgyro_fmany_r_async(nsplitA, fxmany, fA_nl)
  else ! fp32 return
    call cgyro_fmany_r32_async(nsplitA, fxmany, fA_nl32)
  endif

  if (i_triad == 1) then
    call timer_lib_out('nl_mem')
    call timer_lib_in('nl')
    ! 2. Poisson bracket in real space with Non_Zonal pairs

    call cgyro_nl_fftw_mul_sub_mean(size(uvmany,1),size(uvmany,2),nsplitA, &
                         uvmany, &
                         uxmany,vymany(:,:,1:nsplitA), &
                         uymany,vxmany(:,:,1:nsplitA), &
                         inv_nxny)

    ! make sure reqs progress
    call cgyro_nl_fftw_comm_test()

    ! ------------------
    ! Transform uv to fx
    ! ------------------

    call cgyro_fft_d2z(plan_r2c_manyA,uvmany,fxmany)

    call timer_lib_out('nl')
    call timer_lib_in('nl_mem')

    ! tile for performance, since this is effectively a transpose
    call cgyro_fmany_r_async(nsplitA, fxmany, eA_nl)

  endif ! i_triad==1

  if (nsplitB > 0) then
    ! we can zero the elements we know are zero while we waita
    ! assuming nsplitB<=nsplitA
    call cgyro_zero_offdiag_async(nsplitB,fxmany,fymany)
  endif

  call timer_lib_out('nl_mem')

  call timer_lib_in('nl_comm')
  ! start the async reverse comm
  ! can reuse the same req, no overlap with forward fA_req
  if (nl_single_flag .EQ. 0) then
    call parallel_slib_r_nc_async(nsplitA,fA_nl,fpackA,fA_req)
  else
    call parallel_slib_r_nc32_async(nsplitA,fA_nl32,fpackA32,fA_req)
  endif
  fA_req_valid = .TRUE.

  if (i_triad == 1) then
    call parallel_slib_r_nc_async(nsplitA,eA_nl,epackA,eA_req)
  end if

  if (nsplitB > 0) then
    ! time to wait for the 2nd half of F_nl to become avaialble
    call parallel_slib_f_nc_wait(nsplitB,fpackB,fB_nl,fB_req)
    fB_req_valid = .FALSE.
  endif
  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()
  call timer_lib_out('nl_comm')

  if (nsplitB > 0) then

  call timer_lib_in('nl')

! f_nl is (radial, nt_loc, theta, nv_loc1, toroidal_procs)
! where nv_loc1 * toroidal_procs >= nv_loc

  call cgyro_fmany_async(nsplitB, fB_nl)

#if defined(OMPGPU)
  !no async for OMPGPU for now
#else
  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()
!$acc wait(2)
#endif

  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()

  ! Average elements so as to ensure
  !   f(kx,ky=0) = f(-kx,ky=0)^*
  ! This symmetry is required for complex input to c2r
  call cgyro_sym_async(nsplitB,fxmany)

     ! --------------------------------------
     ! perform many Fourier Transforms at once
     ! --------------------------------------

  call cgyro_fft_z2d(plan_c2r_manyB,fymany,uymany)

  ! fxmany is complete now

  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()

  call cgyro_fft_z2d(plan_c2r_manyB,fxmany,uxmany)

  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()

  ! Poisson bracket in real space
  ! uv = (ux*vy-uy*vx)/(nx*ny)

  call cgyro_nl_fftw_mul(size(uvmany,1)*size(uvmany,2)*nsplitB, &
                         uvmany, &
                         uxmany,vymany(:,:,(nsplitA+1):nsplit), &
                         uymany,vxmany(:,:,(nsplitA+1):nsplit), &
                         inv_nxny)

  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()

  ! ------------------
  ! Transform uv to fx
  ! ------------------

#if !defined(OMPGPU)
!$acc wait
#endif

  call cgyro_fft_d2z(plan_r2c_manyB,uvmany,fxmany)

  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()

  call timer_lib_out('nl')
  call timer_lib_in('nl_mem')

  ! NOTE: The FFT will generate an unwanted n=0,p=-nr/2 component
  ! that will be filtered in the main time-stepping loop

  ! tile for performance, since this is effectively a transpose
  if (nl_single_flag .EQ. 0) then
    call cgyro_fmany_r_async(nsplitB, fxmany, fB_nl)
  else ! fp32 return
    call cgyro_fmany_r32_async(nsplitB, fxmany, fB_nl32)
  endif

  if (i_triad == 1) then
    call timer_lib_out('nl_mem')
    call timer_lib_in('nl')
    ! 2. Poisson bracket in real space with Non_Zonal pairs
    call cgyro_nl_fftw_mul_sub_mean(size(uvmany,1),size(uvmany,2),nsplitB, &
                         uvmany, &
                         uxmany,vymany(:,:,(nsplitA+1):nsplit), &
                         uymany,vxmany(:,:,(nsplitA+1):nsplit), &
                         inv_nxny)

    ! make sure reqs progress
    call cgyro_nl_fftw_comm_test()

    ! ------------------
    ! Transform uv to fx
    ! ------------------

    call cgyro_fft_d2z(plan_r2c_manyB,uvmany,fxmany)

    call timer_lib_out('nl')
    call timer_lib_in('nl_mem')

    ! NOTE: The FFT will generate an unwanted n=0,p=-nr/2 component
    ! that will be filtered in the main time-stepping loop

    ! tile for performance, since this is effectively a transpose
    call cgyro_fmany_r_async(nsplitB, fxmany, eB_nl)

  endif ! if i_triad

  call timer_lib_out('nl_mem')

  call timer_lib_in('nl_comm')
  ! start the async reverse comm
  ! can reuse the same req, no overlap with forward fB_req
  if (nl_single_flag .EQ. 0) then
    call parallel_slib_r_nc_async(nsplitB,fB_nl,fpackB,fB_req)
  else
    call parallel_slib_r_nc32_async(nsplitB,fB_nl32,fpackB32,fB_req)
  endif
  if (i_triad == 1) then
     call parallel_slib_r_nc_async(nsplitB,eB_nl,epackB,eB_req)
  end if
  fB_req_valid = .TRUE.
  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()
  call timer_lib_out('nl_comm')

  endif ! nsplitB>0

end subroutine cgyro_nl_fftw_fp64

subroutine cgyro_nl_fftw_fp32(i_triad)

  use timer_lib
  use parallel_lib
  use cgyro_nl_comm
  use cgyro_globals

  implicit none
  !-----------------------------------
  integer, intent(in) :: i_triad
  !-----------------------------------

  real(KIND=REAL32) :: inv_nxny
  inv_nxny = 1.0/(nx*ny)

  call timer_lib_in('nl_mem')
  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()

  ! we can zero the elements we know are zero while we wait
  call cgyro_zero_offdiag32_async(nsplitA,fxmany32,fymany32)

  call timer_lib_out('nl_mem')

  ! time to wait for the FA_nl to become avaialble
  call timer_lib_in('nl_comm')
  call parallel_slib_f_nc32_wait(nsplitA,fpackA32,fA_nl32,fA_req)
  fA_req_valid = .FALSE.
  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()
  call timer_lib_out('nl_comm')

  call timer_lib_in('nl')

! f_nl is (radial, nt_loc, theta, nv_loc1, toroidal_procs)
! where nv_loc1 * toroidal_procs >= nv_loc

  call cgyro_fmany32_async(nsplitA, fA_nl32)

#if defined(OMPGPU)
  !no async for OMPGPU for now
#else
  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()
!$acc wait(2)
#endif

  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()

  ! Average elements so as to ensure
  !   f(kx,ky=0) = f(-kx,ky=0)^*
  ! This symmetry is required for complex input to c2r
  call cgyro_sym32_async(nsplitA,fxmany32)

     ! --------------------------------------
     ! perform many Fourier Transforms at once
     ! --------------------------------------

  call cgyro_fft_c2r(plan_c2r_manyA,fymany32,uymany32)
  ! fxmany is complete now

  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()

#if !defined(OMPGPU)
!$acc wait
#endif

  call cgyro_fft_c2r(plan_c2r_manyA,fxmany32,uxmany32)

  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()

  ! we can zero the elements we know are zero while we wait for comm
  call cgyro_zero_offdiag32_async(nsplit,gxmany32,gymany32)

  call timer_lib_out('nl')

  call timer_lib_in('nl_comm')
  ! time to wait for the g_nl to become avaialble
  call parallel_slib_f_fd32_wait(n_field,n_radial,n_jtheta,gpack32,g_nl32,g_req)
  g_req_valid = .FALSE.
  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()
  call timer_lib_out('nl_comm')

  call timer_lib_in('nl')

! g_nl      is (n_field,n_radial,n_jtheta,nt_loc,n_toroidal_procs)
! jcev_c_nl is (n_field,n_radial,n_jtheta,nv_loc,nt_loc,n_toroidal_procs)

  call cgyro_gmany32_async

#if defined(OMPGPU)
  !no async for OMPGPU for now
#else
  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()
!$acc wait(2)
#endif

  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()

  ! Average elements so as to ensure
  !   g(kx,ky=0) = g(-kx,ky=0)^*
  ! This symmetry is required for complex input to c2r
  call cgyro_sym32_async(nsplit,gxmany32)

  call cgyro_fft_c2r(plan_c2r_manyG,gymany32,vymany32)

  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()
  ! gxmany is complete now

  call cgyro_fft_c2r(plan_c2r_manyG,gxmany32,vxmany32)

  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()

  ! Poisson bracket in real space
  ! uv = (ux*vy-uy*vx)/(nx*ny)

  call cgyro_nl_fftw_mul32(size(uvmany32,1)*size(uvmany32,2)*nsplitA, &
                         uvmany32, &
                         uxmany32,vymany32(:,:,1:nsplitA), &
                         uymany32,vxmany32(:,:,1:nsplitA), &
                         inv_nxny)

  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()

  ! ------------------
  ! Transform uv to fx
  ! ------------------

#if !defined(OMPGPU)
!$acc wait
#endif

  call cgyro_fft_r2c(plan_r2c_manyA,uvmany32,fxmany32)

  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()

  call timer_lib_out('nl')
  call timer_lib_in('nl_mem')

  ! NOTE: The FFT will generate an unwanted n=0,p=-nr/2 component
  ! that will be filtered in the main time-stepping loop

  ! tile for performance, since this is effectively a transpose
  call cgyro_fmany32_r32_async(nsplitA, fxmany32, fA_nl32)

  if (nsplitB > 0) then
    ! we can zero the elements we know are zero while we waita
    ! assuming nsplitB<=nsplitA
    call cgyro_zero_offdiag32_async(nsplitB,fxmany32,fymany32)
  endif

  call timer_lib_out('nl_mem')

  call timer_lib_in('nl_comm')
  ! start the async reverse comm
  ! can reuse the same req, no overlap with forward fA_req
  call parallel_slib_r_nc32_async(nsplitA,fA_nl32,fpackA32,fA_req)
  fA_req_valid = .TRUE.

  if (nsplitB > 0) then
    ! time to wait for the 2nd half of F_nl to become avaialble
    call parallel_slib_f_nc32_wait(nsplitB,fpackB32,fB_nl32,fB_req)
    fB_req_valid = .FALSE.
  endif
  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()
  call timer_lib_out('nl_comm')

  if (nsplitB > 0) then

  call timer_lib_in('nl')

! f_nl is (radial, nt_loc, theta, nv_loc1, toroidal_procs)
! where nv_loc1 * toroidal_procs >= nv_loc

  call cgyro_fmany32_async(nsplitB, fB_nl32)

#if defined(OMPGPU)
  !no async for OMPGPU for now
#else
  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()
!$acc wait(2)
#endif

  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()

  ! Average elements so as to ensure
  !   f(kx,ky=0) = f(-kx,ky=0)^*
  ! This symmetry is required for complex input to c2r
  call cgyro_sym32_async(nsplitB,fxmany32)

     ! --------------------------------------
     ! perform many Fourier Transforms at once
     ! --------------------------------------

  call cgyro_fft_c2r(plan_c2r_manyB,fymany32,uymany32)

  ! fxmany is complete now

  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()

  call cgyro_fft_c2r(plan_c2r_manyB,fxmany32,uxmany32)

  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()

  ! Poisson bracket in real space
  ! uv = (ux*vy-uy*vx)/(nx*ny)

  call cgyro_nl_fftw_mul32(size(uvmany32,1)*size(uvmany32,2)*nsplitB, &
                         uvmany32, &
                         uxmany32,vymany32(:,:,(nsplitA+1):nsplit), &
                         uymany32,vxmany32(:,:,(nsplitA+1):nsplit), &
                         inv_nxny)

  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()

  ! ------------------
  ! Transform uv to fx
  ! ------------------

#if !defined(OMPGPU)
!$acc wait
#endif

  call cgyro_fft_r2c(plan_r2c_manyB,uvmany32,fxmany32)

  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()

  call timer_lib_out('nl')
  call timer_lib_in('nl_mem')

  ! NOTE: The FFT will generate an unwanted n=0,p=-nr/2 component
  ! that will be filtered in the main time-stepping loop

  ! tile for performance, since this is effectively a transpose
  call cgyro_fmany32_r32_async(nsplitB, fxmany32, fB_nl32)

  call timer_lib_out('nl_mem')

  call timer_lib_in('nl_comm')
  ! start the async reverse comm
  ! can reuse the same req, no overlap with forward fB_req
  call parallel_slib_r_nc32_async(nsplitB,fB_nl32,fpackB32,fB_req)
  fB_req_valid = .TRUE.
  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()
  call timer_lib_out('nl_comm')

  endif ! nsplitB>0

end subroutine cgyro_nl_fftw_fp32

#else  /* if not definedCGYRO_GPU_FFT */

!-----------------------------------------------------------------
!  CPU code for cgyro_nl_fftw and dependent sub-functions
!-----------------------------------------------------------------

!-----------------------------------------------------------------
! Helper function for computing the reverse FFT transform
! and then convert the output into the required MPI format
! for a single i_omp partition from a single iteration of the j loop.
! May be used for both fA_nl and fB_nl.
! Single threaded and thread safe, as it is 
! expected to be called from inside a OMP loop.
!-----------------------------------------------------------------

subroutine cgyro_nl_fftw_stepr(g_j, f_j, nl_idx, i_omp)

  use timer_lib
  use parallel_lib
  use cgyro_globals

  implicit none

  integer, intent(in) :: g_j, f_j
  integer,intent(in) :: nl_idx ! 1=>A, 2=>B
  integer,intent(in) :: i_omp
  integer :: ix,iy
  integer :: ir,itm,itl,itor

  include 'fftw3.f03'

  ! Poisson bracket in real space
  uv(:,:,i_omp) = (uxmany(:,:,f_j)*vymany(:,:,g_j)-uymany(:,:,f_j)*vxmany(:,:,g_j))/(nx*ny)

  call fftw_execute_dft_r2c(plan_r2c,uv(:,:,i_omp),fx(:,:,i_omp))

  ! NOTE: The FFT will generate an unwanted n=0,p=-nr/2 component
  ! that will be filtered in the main time-stepping loop

  ! this should really be accounted against nl_mem, but hard to do with OMP
  if (nl_single_flag .EQ. 0) then

   if (nl_idx==1) then
    do itm=1,n_toroidal_procs
     do itl=1,nt_loc
      itor=itl + (itm-1)*nt_loc
      do ir=1,n_radial
        ix = ir-1-nx0/2
        if (ix < 0) ix = ix+nx
        iy = itor-1
        fA_nl(ir,itl,f_j,itm) = fx(iy,ix,i_omp)
      enddo
     enddo
    enddo
   else ! nl_idx>1
    do itm=1,n_toroidal_procs
     do itl=1,nt_loc
      itor=itl + (itm-1)*nt_loc
      do ir=1,n_radial
        ix = ir-1-nx0/2
        if (ix < 0) ix = ix+nx
        iy = itor-1
        fB_nl(ir,itl,f_j,itm) = fx(iy,ix,i_omp)
      enddo
     enddo
    enddo
   endif

  else ! fp32 return

   if (nl_idx==1) then
    do itm=1,n_toroidal_procs
     do itl=1,nt_loc
      itor=itl + (itm-1)*nt_loc
      do ir=1,n_radial
        ix = ir-1-nx0/2
        if (ix < 0) ix = ix+nx
        iy = itor-1
        fA_nl32(ir,itl,f_j,itm) = fx(iy,ix,i_omp)
      enddo
     enddo
    enddo
   else ! nl_idx>1
    do itm=1,n_toroidal_procs
     do itl=1,nt_loc
      itor=itl + (itm-1)*nt_loc
      do ir=1,n_radial
        ix = ir-1-nx0/2
        if (ix < 0) ix = ix+nx
        iy = itor-1
        fB_nl32(ir,itl,f_j,itm) = fx(iy,ix,i_omp)
      enddo
     enddo
    enddo
   endif

  endif

end subroutine cgyro_nl_fftw_stepr

subroutine cgyro_nl_fftw_stepr32(g_j, f_j, nl_idx, i_omp)

  use timer_lib
  use parallel_lib
  use cgyro_globals

  implicit none

  integer, intent(in) :: g_j, f_j
  integer,intent(in) :: nl_idx ! 1=>A, 2=>B
  integer,intent(in) :: i_omp
  integer :: ix,iy
  integer :: ir,itm,itl,itor

  include 'fftw3.f03'

  ! Poisson bracket in real space
  uv32(:,:,i_omp) = (uxmany32(:,:,f_j)*vymany32(:,:,g_j)-uymany32(:,:,f_j)*vxmany32(:,:,g_j))/(nx*ny)

  call fftwf_execute_dft_r2c(plan_r2c,uv32(:,:,i_omp),fx32(:,:,i_omp))

  ! NOTE: The FFT will generate an unwanted n=0,p=-nr/2 component
  ! that will be filtered in the main time-stepping loop

  ! this should really be accounted against nl_mem, but hard to do with OMP
  if (nl_idx==1) then
    do itm=1,n_toroidal_procs
     do itl=1,nt_loc
      itor=itl + (itm-1)*nt_loc
      do ir=1,n_radial
        ix = ir-1-nx0/2
        if (ix < 0) ix = ix+nx
        iy = itor-1
        fA_nl32(ir,itl,f_j,itm) = fx32(iy,ix,i_omp)
      enddo
     enddo
    enddo
  else ! nl_idx>1
    do itm=1,n_toroidal_procs
     do itl=1,nt_loc
      itor=itl + (itm-1)*nt_loc
      do ir=1,n_radial
        ix = ir-1-nx0/2
        if (ix < 0) ix = ix+nx
        iy = itor-1
        fB_nl32(ir,itl,f_j,itm) = fx32(iy,ix,i_omp)
      enddo
     enddo
    enddo
  endif

end subroutine cgyro_nl_fftw_stepr32

  ! assumes cgyro_nl_fftw_stepr has already been called
subroutine cgyro_nl_fftw_stepr_triad(g_j, f_j, nl_idx, i_omp)

  use timer_lib
  use parallel_lib
  use cgyro_globals

  implicit none

  !-----------------------------------
  integer, intent(in) :: g_j, f_j
  integer,intent(in) :: nl_idx ! 1=>A, 2=>B
  integer,intent(in) :: i_omp
  !-----------------------------------
  real :: y_mean_ux, y_mean_uy, y_mean_vx, y_mean_vy
  real :: r_ux, r_uy, r_vx, r_vy
  integer :: ix,iy
  integer :: ir,itm,itl,itor
  real :: inv_nxny,inv_ny

  include 'fftw3.f03'

  inv_ny = 1.0/ny
  inv_nxny = 1.0/(nx*ny)
  ! 1. Poisson bracket in real space was done in cgyro_nl_fftw_stepr

  ! 2. Poisson bracket in real space with Non_Zonal pairs

  do ix=0,nx-1
    y_mean_ux = sum(uxmany(:,ix,f_j) ,dim=1 ) * inv_ny
    y_mean_uy = sum(uymany(:,ix,f_j) ,dim=1 ) * inv_ny
    y_mean_vx = sum(vxmany(:,ix,g_j) ,dim=1 ) * inv_ny
    y_mean_vy = sum(vymany(:,ix,g_j) ,dim=1 ) * inv_ny
    do iy=0,ny-1
      ! remove ky=0
      r_ux = uxmany(iy,ix,f_j) - y_mean_ux
      r_uy = uymany(iy,ix,f_j) - y_mean_uy
      r_vx = vxmany(iy,ix,g_j) - y_mean_vx
      r_vy = vymany(iy,ix,g_j) - y_mean_vy
      ! we could save, but we do not really need to (not used outside from this function)
      !uxmany(iy,ix,f_j) = r_ux
      !uymany(iy,ix,f_j) = r_y
      !vxmany(iy,ix,g_j) = r_vx
      !vymany(iy,ix,g_j) = r_vy

      ! compute and save uv
      uv(iy,ix,i_omp) = (r_ux*r_vy-r_uy*r_vx) * inv_nxny
    enddo
  enddo

  call fftw_execute_dft_r2c(plan_r2c,uv(:,:,i_omp),fx(:,:,i_omp))

  do itm=1,n_toroidal_procs
   do itl=1,nt_loc
    itor=itl + (itm-1)*nt_loc
    do ir=1,n_radial
     ix = ir-1-nx0/2
     if (ix < 0) ix = ix+nx
     iy = itor-1
     if (nl_idx==1) then
       eA_nl(ir,itl,f_j,itm) = fx(iy,ix,i_omp)
     else
       eB_nl(ir,itl,f_j,itm) = fx(iy,ix,i_omp)
     endif
    enddo
   enddo
  enddo

end subroutine cgyro_nl_fftw_stepr_triad

!-----------------------------------------------------------------
! Helper function for filling fx and fy (inputs to FFT)
! for a single i_omp partition from a single iteration of the j loop.
! The input matrix can be either fN_nl or fB_nl, with nj either nsplitA or nsplitB.
! Single threaded and thread safe, as it is 
! expected to be called from inside a OMP loop.
!-----------------------------------------------------------------

subroutine cgyro_f_one(i_omp, j, nj, f_nl)

  use cgyro_nl_comm
  use cgyro_globals

  implicit none
  !-----------------------------------
  integer, intent(in) :: i_omp, j, nj
  complex, dimension(n_radial,nt_loc,nj,n_toroidal_procs), intent(in) :: f_nl
  !-----------------------------------
  integer :: ix,iy
  integer :: ir,itm,itl
  integer :: itor
  integer :: p

  complex :: f0

  include 'fftw3.f03'

! f_nl is (radial, nt_loc, theta, nv_loc1, toroidal_procs)
! where nv_loc1 * toroidal_procs >= nv_loc
        ! zero elements not otherwise set below
        fx(0:ny2,nx2:nx0-1,i_omp) = 0.0
        fy(0:ny2,nx2:nx0-1,i_omp) = 0.0

        ! Array mapping
        do ir=1,n_radial
           p  = ir-1-nx0/2
           ix = p
           if (ix < 0) ix = ix+nx
           do itm=1,n_toroidal_procs
            do itl=1,nt_loc
              itor=itl + (itm-1)*nt_loc
              iy = itor-1
              f0 = i_c*f_nl(ir,itl,j,itm)
              fx(iy,ix,i_omp) = p*f0
              fy(iy,ix,i_omp) = iy*f0
            enddo
           enddo
           if ((ix/=0) .and. (ix<(nx/2))) then ! happens after ix>nx/2
             ! Average elements so as to ensure
             !   f(kx,ky=0) = f(-kx,ky=0)^*
             ! This symmetry is required for complex input to c2r
             f0 = 0.5*( fx(0,ix,i_omp)+conjg(fx(0,nx-ix,i_omp)) )
             fx(0,ix   ,i_omp) = f0
             fx(0,nx-ix,i_omp) = conjg(f0)
           endif
           fx(n_toroidal:ny2,ix,i_omp) = 0.0
           fy(n_toroidal:ny2,ix,i_omp) = 0.0
        enddo

        if (i_omp==1) then
         ! use the main thread to progress the async MPI
          call cgyro_nl_fftw_comm_test()
        endif

        call fftw_execute_dft_c2r(plan_c2r,fx(:,:,i_omp),uxmany(:,:,j))
        if (i_omp==1) then
         ! use the main thread to progress the async MPI
          call cgyro_nl_fftw_comm_test()
        endif

        call fftw_execute_dft_c2r(plan_c2r,fy(:,:,i_omp),uymany(:,:,j))

end subroutine cgyro_f_one

subroutine cgyro_f32_one(i_omp, j, nj, f_nl32)

  use cgyro_nl_comm
  use cgyro_globals

  implicit none
  !-----------------------------------
  integer, intent(in) :: i_omp, j, nj
  complex(KIND=REAL32), dimension(n_radial,nt_loc,nj,n_toroidal_procs), intent(in) :: f_nl32
  !-----------------------------------
  integer :: ix,iy
  integer :: ir,itm,itl
  integer :: itor
  integer :: p

  complex(KIND=REAL32) :: f0

  include 'fftw3.f03'

! f_nl is (radial, nt_loc, theta, nv_loc1, toroidal_procs)
! where nv_loc1 * toroidal_procs >= nv_loc
        ! zero elements not otherwise set below
        fx32(0:ny2,nx2:nx0-1,i_omp) = 0.0
        fy32(0:ny2,nx2:nx0-1,i_omp) = 0.0

        ! Array mapping
        do ir=1,n_radial
           p  = ir-1-nx0/2
           ix = p
           if (ix < 0) ix = ix+nx
           do itm=1,n_toroidal_procs
            do itl=1,nt_loc
              itor=itl + (itm-1)*nt_loc
              iy = itor-1
              f0 = i_c*f_nl32(ir,itl,j,itm)
              fx32(iy,ix,i_omp) = p*f0
              fy32(iy,ix,i_omp) = iy*f0
            enddo
           enddo
           if ((ix/=0) .and. (ix<(nx/2))) then ! happens after ix>nx/2
             ! Average elements so as to ensure
             !   f(kx,ky=0) = f(-kx,ky=0)^*
             ! This symmetry is required for complex input to c2r
             f0 = 0.5*( fx32(0,ix,i_omp)+conjg(fx32(0,nx-ix,i_omp)) )
             fx32(0,ix   ,i_omp) = f0
             fx32(0,nx-ix,i_omp) = conjg(f0)
           endif
           fx32(n_toroidal:ny2,ix,i_omp) = 0.0
           fy32(n_toroidal:ny2,ix,i_omp) = 0.0
        enddo

        if (i_omp==1) then
         ! use the main thread to progress the async MPI
          call cgyro_nl_fftw_comm_test()
        endif

        call fftwf_execute_dft_c2r(plan_c2r,fx32(:,:,i_omp),uxmany32(:,:,j))
        if (i_omp==1) then
         ! use the main thread to progress the async MPI
          call cgyro_nl_fftw_comm_test()
        endif

        call fftwf_execute_dft_c2r(plan_c2r,fy32(:,:,i_omp),uymany32(:,:,j))

end subroutine cgyro_f32_one

!-----------------------------------------------------------------
! Helper function for filling gx and gy (inputs to FFT)
! for a single i_omp partition from a single iteration of the j loop.
! Single threaded and thread safe, as it is 
! expected to be called from inside a OMP loop.
!-----------------------------------------------------------------

subroutine cgyro_g_one(i_omp, j)

  use cgyro_nl_comm
  use cgyro_globals

  implicit none
  !-----------------------------------
  integer, intent(in) :: i_omp, j
  !-----------------------------------
  integer :: ix,iy
  integer :: ir,itm,itl
  integer :: it,itor,mytm,it_loc
  integer :: p
  integer :: jtheta_min

  complex :: g0

  include 'fftw3.f03'

! g_nl      is (n_field,n_radial,n_jtheta,nt_loc,n_toroidal_procs)
! jcev_c_nl is (n_field,n_radial,n_jtheta,nv_loc,nt_loc,n_toroidal_procs)
        ! zero elements not otherwise set below
        gx(0:ny2,nx2:nx0-1,i_omp) = 0.0
        gy(0:ny2,nx2:nx0-1,i_omp) = 0.0

        ! Array mapping
        do ir=1,n_radial
           p  = ir-1-nx0/2
           ix = p
           if (ix < 0) ix = ix+nx
           do itm=1,n_toroidal_procs
            do itl=1,nt_loc
              itor = itl + (itm-1)*nt_loc
              mytm = 1 + nt1/nt_loc !my toroidal proc number
              it = 1+((mytm-1)*nsplit+j-1)/nv_loc
              iv_loc = 1+modulo((mytm-1)*nsplit+j-1,nv_loc)
              jtheta_min = 1+((mytm-1)*nsplit)/nv_loc
              it_loc = it-jtheta_min+1

              iy = itor-1
              if (it > n_theta) then
                 g0 = (0.0,0.0)
              else
                 g0 = i_c*sum( jvec_c_nl(:,ir,it_loc,iv_loc,itor)*g_nl(:,ir,it_loc,itor))
              endif
              gx(iy,ix,i_omp) = p*g0
              gy(iy,ix,i_omp) = iy*g0
            enddo
           enddo
           if ((ix/=0) .and. (ix<(nx/2))) then ! happens after ix>nx/2
              ! Average elements so as to ensure
              !   g(kx,ky=0) = g(-kx,ky=0)^*
              ! This symmetry is required for complex input to c2r
              g0 = 0.5*( gx(0,ix,i_omp)+conjg(gx(0,nx-ix,i_omp)) )
              gx(0,ix   ,i_omp) = g0
              gx(0,nx-ix,i_omp) = conjg(g0)
           endif
           gx(n_toroidal:ny2,ix,i_omp) = 0.0
           gy(n_toroidal:ny2,ix,i_omp) = 0.0
        enddo

        if (i_omp==1) then
         ! use the main thread to progress the async MPI
          call cgyro_nl_fftw_comm_test()
        endif

        call fftw_execute_dft_c2r(plan_c2r,gx(:,:,i_omp),vxmany(:,:,j))

        if (i_omp==1) then
         ! use the main thread to progress the async MPI
          call cgyro_nl_fftw_comm_test()
        endif

        call fftw_execute_dft_c2r(plan_c2r,gy(:,:,i_omp),vymany(:,:,j))
end subroutine cgyro_g_one

subroutine cgyro_g32_one(i_omp, j)

  use cgyro_nl_comm
  use cgyro_globals

  implicit none
  !-----------------------------------
  integer, intent(in) :: i_omp, j
  !-----------------------------------
  integer :: ix,iy
  integer :: ir,itm,itl
  integer :: it,itor,mytm,it_loc
  integer :: p
  integer :: jtheta_min

  complex(KIND=REAL32) :: g0

  include 'fftw3.f03'

! g_nl      is (n_field,n_radial,n_jtheta,nt_loc,n_toroidal_procs)
! jcev_c_nl is (n_field,n_radial,n_jtheta,nv_loc,nt_loc,n_toroidal_procs)
        ! zero elements not otherwise set below
        gx32(0:ny2,nx2:nx0-1,i_omp) = 0.0
        gy32(0:ny2,nx2:nx0-1,i_omp) = 0.0

        ! Array mapping
        do ir=1,n_radial
           p  = ir-1-nx0/2
           ix = p
           if (ix < 0) ix = ix+nx
           do itm=1,n_toroidal_procs
            do itl=1,nt_loc
              itor = itl + (itm-1)*nt_loc
              mytm = 1 + nt1/nt_loc !my toroidal proc number
              it = 1+((mytm-1)*nsplit+j-1)/nv_loc
              iv_loc = 1+modulo((mytm-1)*nsplit+j-1,nv_loc)
              jtheta_min = 1+((mytm-1)*nsplit)/nv_loc
              it_loc = it-jtheta_min+1

              iy = itor-1
              if (it > n_theta) then
                 g0 = (0.0,0.0)
              else
                 g0 = i_c*sum( jvec_c_nl32(:,ir,it_loc,iv_loc,itor)*g_nl32(:,ir,it_loc,itor))
              endif
              gx32(iy,ix,i_omp) = p*g0
              gy32(iy,ix,i_omp) = iy*g0
            enddo
           enddo
           if ((ix/=0) .and. (ix<(nx/2))) then ! happens after ix>nx/2
              ! Average elements so as to ensure
              !   g(kx,ky=0) = g(-kx,ky=0)^*
              ! This symmetry is required for complex input to c2r
              g0 = 0.5*( gx32(0,ix,i_omp)+conjg(gx32(0,nx-ix,i_omp)) )
              gx32(0,ix   ,i_omp) = g0
              gx32(0,nx-ix,i_omp) = conjg(g0)
           endif
           gx32(n_toroidal:ny2,ix,i_omp) = 0.0
           gy32(n_toroidal:ny2,ix,i_omp) = 0.0
        enddo

        if (i_omp==1) then
         ! use the main thread to progress the async MPI
          call cgyro_nl_fftw_comm_test()
        endif

        call fftwf_execute_dft_c2r(plan_c2r,gx32(:,:,i_omp),vxmany32(:,:,j))

        if (i_omp==1) then
         ! use the main thread to progress the async MPI
          call cgyro_nl_fftw_comm_test()
        endif

        call fftwf_execute_dft_c2r(plan_c2r,gy32(:,:,i_omp),vymany32(:,:,j))
end subroutine cgyro_g32_one

!-----------------------------------------------------------------
! High-level cgyro_nl_fftw logic
! Do the actual FFT-base NL step
! NOTE: Expects cgyro_nl_fftw_comm1 to have been already called
!-----------------------------------------------------------------

subroutine cgyro_nl_fftw_fp64(i_triad)

  use timer_lib
  use parallel_lib
  use cgyro_nl_comm
  use cgyro_globals

  implicit none

  !-----------------------------------
  integer, intent(in) :: i_triad
  !-----------------------------------
  integer :: ix,iy
  integer :: ir,it,itm,itl,it_loc
  integer :: itor,mytm
  integer :: j,p
  integer :: i_omp
  integer :: jtheta_min

  complex :: f0,g0

  integer, external :: omp_get_thread_num

  include 'fftw3.f03'
  
  call cgyro_nl_fftw_comm_test()

  ! time to wait for the FA_nl to become avaialble
  call timer_lib_in('nl_comm')
  call parallel_slib_f_nc_wait(nsplitA,fpackA,fA_nl,fA_req)
  fA_req_valid = .FALSE.
  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()
  call timer_lib_out('nl_comm')

  call timer_lib_in('nl')

! f_nl is (radial, nt_loc, theta, nv_loc1, toroidal_procs)
! where nv_loc1 * toroidal_procs >= nv_loc
! parallel do schedule(dynamic,1) private(itm,itl,itor,iy,ir,p,ix,f0,i_omp,j)
!$omp parallel do schedule(dynamic,1) private(i_omp,j)
  do j=1,nsplitA
        i_omp = omp_get_thread_num()+1
        call cgyro_f_one(i_omp, j, nsplitA, fA_nl)
  enddo ! j

  call timer_lib_out('nl')

  call timer_lib_in('nl_comm')
  ! time to wait for the g_nl to become avaialble
  call parallel_slib_f_fd_wait(n_field,n_radial,n_jtheta,gpack,g_nl,g_req)
  g_req_valid = .FALSE.
  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()
  call timer_lib_out('nl_comm')

  call timer_lib_in('nl')

! g_nl      is (n_field,n_radial,n_jtheta,nt_loc,n_toroidal_procs)
! jcev_c_nl is (n_field,n_radial,n_jtheta,nv_loc,nt_loc,n_toroidal_procs)
!$omp parallel do schedule(dynamic,1) private(i_omp,j)
  do j=1,nsplit
        i_omp = omp_get_thread_num()+1
        call cgyro_g_one(i_omp, j)
        if (i_omp==1) then
         ! use the main thread to progress the async MPI
          call cgyro_nl_fftw_comm_test()
        endif

        if (j<=nsplitA) then
           call cgyro_nl_fftw_stepr(j, j, 1, i_omp)
           if (i_omp==1) then
             ! use the main thread to progress the async MPI
             call cgyro_nl_fftw_comm_test()
           endif

           if (i_triad == 1) then
              call cgyro_nl_fftw_stepr_triad(j, j, 1, i_omp)
           endif
        endif
        ! else we will do it in the next loop
  enddo ! j

  call timer_lib_out('nl')

  call timer_lib_in('nl_comm')
  ! start the async reverse comm
  ! can reuse the same req, no overlap with forward fA_req
  if (nl_single_flag .EQ. 0) then
    call parallel_slib_r_nc_async(nsplitA,fA_nl,fpackA,fA_req)
  else
    call parallel_slib_r_nc32_async(nsplitA,fA_nl32,fpackA32,fA_req)
  endif
   if (i_triad == 1) then
     call parallel_slib_r_nc_async(nsplitA,eA_nl,epackA,eA_req)
   end if
  fA_req_valid = .TRUE.

  if (nsplitB > 0) then
    ! time to wait for the 2nd half of F_nl to become avaialble
    call parallel_slib_f_nc_wait(nsplitB,fpackB,fB_nl,fB_req)
    fB_req_valid = .FALSE.
  endif

  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()
  call timer_lib_out('nl_comm')

  if (nsplitB > 0) then
   call timer_lib_in('nl')

! f_nl is (radial, nt_loc, theta, nv_loc1, toroidal_procs)
! where nv_loc1 * toroidal_procs >= nv_loc
!$omp parallel do schedule(dynamic,1) private(i_omp,j)
   do j=1,nsplitB
        i_omp = omp_get_thread_num()+1
        call cgyro_f_one(i_omp, j, nsplitB, fB_nl)

        if (i_omp==1) then
         ! use the main thread to progress the async MPI
          call cgyro_nl_fftw_comm_test()
        endif

        call cgyro_nl_fftw_stepr(nsplitA+j, j, 2, i_omp)
        if (i_omp==1) then
         ! use the main thread to progress the async MPI
          call cgyro_nl_fftw_comm_test()
        endif

        if (i_triad == 1) then
          call cgyro_nl_fftw_stepr_triad(nsplitA+j, j, 2, i_omp)
        endif
   enddo ! j

   call timer_lib_out('nl')

   call timer_lib_in('nl_comm')
   ! start the async reverse comm
   ! can reuse the same req, no overlap with forward fB_req
   if (nl_single_flag .EQ. 0) then
     call parallel_slib_r_nc_async(nsplitB,fB_nl,fpackB,fB_req)
   else
     call parallel_slib_r_nc32_async(nsplitB,fB_nl32,fpackB32,fB_req)
   endif
   if (i_triad == 1) then
     call parallel_slib_r_nc_async(nsplitB,eB_nl,epackB,eB_req)
   end if
   fB_req_valid = .TRUE.
   ! make sure reqs progress
   call cgyro_nl_fftw_comm_test()
   call timer_lib_out('nl_comm')
  endif ! if nsplitB>0

end subroutine cgyro_nl_fftw_fp64

subroutine cgyro_nl_fftw_fp32(i_triad)

  use timer_lib
  use parallel_lib
  use cgyro_nl_comm
  use cgyro_globals

  implicit none

  !-----------------------------------
  integer, intent(in) :: i_triad
  !-----------------------------------
  integer :: ix,iy
  integer :: ir,it,itm,itl,it_loc
  integer :: itor,mytm
  integer :: j,p
  integer :: i_omp
  integer :: jtheta_min

  complex(KIND=REAL32) :: f0,g0

  integer, external :: omp_get_thread_num

  include 'fftw3.f03'
  
  call cgyro_nl_fftw_comm_test()

  ! time to wait for the FA_nl to become avaialble
  call timer_lib_in('nl_comm')
  call parallel_slib_f_nc32_wait(nsplitA,fpackA32,fA_nl32,fA_req)
  fA_req_valid = .FALSE.
  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()
  call timer_lib_out('nl_comm')

  call timer_lib_in('nl')

! f_nl is (radial, nt_loc, theta, nv_loc1, toroidal_procs)
! where nv_loc1 * toroidal_procs >= nv_loc
! parallel do schedule(dynamic,1) private(itm,itl,itor,iy,ir,p,ix,f0,i_omp,j)
!$omp parallel do schedule(dynamic,1) private(i_omp,j)
  do j=1,nsplitA
        i_omp = omp_get_thread_num()+1
        call cgyro_f32_one(i_omp, j, nsplitA, fA_nl32)
  enddo ! j

  call timer_lib_out('nl')

  call timer_lib_in('nl_comm')
  ! time to wait for the g_nl to become avaialble
  call parallel_slib_f_fd32_wait(n_field,n_radial,n_jtheta,gpack32,g_nl32,g_req)
  g_req_valid = .FALSE.
  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()
  call timer_lib_out('nl_comm')

  call timer_lib_in('nl')

! g_nl      is (n_field,n_radial,n_jtheta,nt_loc,n_toroidal_procs)
! jcev_c_nl is (n_field,n_radial,n_jtheta,nv_loc,nt_loc,n_toroidal_procs)
!$omp parallel do schedule(dynamic,1) private(i_omp,j)
  do j=1,nsplit
        i_omp = omp_get_thread_num()+1
        call cgyro_g32_one(i_omp, j)
        if (i_omp==1) then
         ! use the main thread to progress the async MPI
          call cgyro_nl_fftw_comm_test()
        endif

        if (j<=nsplitA) then
           call cgyro_nl_fftw_stepr32(j, j, 1, i_omp)
        endif
        ! else we will do it in the next loop
  enddo ! j

  call timer_lib_out('nl')

  call timer_lib_in('nl_comm')
  ! start the async reverse comm
  ! can reuse the same req, no overlap with forward fA_req
  call parallel_slib_r_nc32_async(nsplitA,fA_nl32,fpackA32,fA_req)
  fA_req_valid = .TRUE.

  if (nsplitB > 0) then
    ! time to wait for the 2nd half of F_nl to become avaialble
    call parallel_slib_f_nc32_wait(nsplitB,fpackB32,fB_nl32,fB_req)
    fB_req_valid = .FALSE.
  endif

  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()
  call timer_lib_out('nl_comm')

  if (nsplitB > 0) then
   call timer_lib_in('nl')

! f_nl is (radial, nt_loc, theta, nv_loc1, toroidal_procs)
! where nv_loc1 * toroidal_procs >= nv_loc
!$omp parallel do schedule(dynamic,1) private(i_omp,j)
   do j=1,nsplitB
        i_omp = omp_get_thread_num()+1
        call cgyro_f32_one(i_omp, j, nsplitB, fB_nl32)

        if (i_omp==1) then
         ! use the main thread to progress the async MPI
          call cgyro_nl_fftw_comm_test()
        endif

        call cgyro_nl_fftw_stepr32(nsplitA+j, j, 2, i_omp)
   enddo ! j

   call timer_lib_out('nl')

   call timer_lib_in('nl_comm')
   ! start the async reverse comm
   ! can reuse the same req, no overlap with forward fB_req
   call parallel_slib_r_nc32_async(nsplitB,fB_nl32,fpackB32,fB_req)
   fB_req_valid = .TRUE.
   ! make sure reqs progress
   call cgyro_nl_fftw_comm_test()
   call timer_lib_out('nl_comm')
  endif ! if nsplitB>0

end subroutine cgyro_nl_fftw_fp32

#endif /* CGYRO_GPU_FFT */

!-----------------------------------------------------------------
!  Shared code between CPU and GPU
!-----------------------------------------------------------------

! Initialize all the relevant plan_r2c and plan_c2r variants
subroutine cgyro_nl_fftw_init

  use cgyro_globals

  implicit none
  !-----------------------------------

  if (nl_single_flag > 1) then
    call cgyro_nl_fftw_init_fp32
  else
    call cgyro_nl_fftw_init_fp64
  endif
end subroutine cgyro_nl_fftw_init

! Do the actual FFT-base NL step
! NOTE: call cgyro_nl_fftw_comm1 before cgyro_nl_fftw
subroutine cgyro_nl_fftw(i_triad)

  use cgyro_globals

  implicit none
  !-----------------------------------
  integer, intent(in) :: i_triad
  !-----------------------------------

  if (nl_single_flag > 1) then
    if (i_triad == 1) then
       ! TODO: Implement this use case
       write(*,*) "NOT IMPLEMENTED: Cannot use triad with nl_single"
       call abort
    endif
    call cgyro_nl_fftw_fp32(i_triad)
  else
    call cgyro_nl_fftw_fp64(i_triad)
  endif
end subroutine cgyro_nl_fftw

end module cgyro_nl

