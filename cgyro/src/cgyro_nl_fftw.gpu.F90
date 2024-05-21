!-----------------------------------------------------------------
! cgyro_nl_fftw.gpu.f90 [GPU (acc-cuFFT/hipFFT) version]
!
! PURPOSE:
!  Evaluate nonlinear bracket with dealiased FFT.  It is natural 
!  to use the FFTW complex-to-real (c2r) transform to move to real 
!  space, compute the convolution (uv), and transform back to the 
!  spectral form using a real-to-complex transform (r2c).
!  
! NOTE: Need to be careful with (p=-nr/2,n=0) component.
!-----------------------------------------------------------------

#if defined(MKLGPU)
include 'fftw/offload/fftw3_omp_offload.f90'
#endif

subroutine cgyro_nl_fftw_init

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
  hip_plan_c2r_manyA = c_null_ptr
  istatus = hipfftPlanMany(&
       hip_plan_c2r_manyA, &
       irank, &
       ndim, &
       inembed, &
       istride, &
       idist, &
       onembed, &
       ostride, &
       odist, &
       merge(HIPFFT_C2R,HIPFFT_Z2D,kind(uxmany) == singlePrecision), &
       nsplitA)

  if (nsplitB > 0) then ! no fft if nsplitB==0
    hip_plan_c2r_manyB = c_null_ptr
    istatus = hipfftPlanMany(&
       hip_plan_c2r_manyB, &
       irank, &
       ndim, &
       inembed, &
       istride, &
       idist, &
       onembed, &
       ostride, &
       odist, &
       merge(HIPFFT_C2R,HIPFFT_Z2D,kind(uxmany) == singlePrecision), &
       nsplitB)
  endif

  hip_plan_c2r_manyG = c_null_ptr
  istatus = hipfftPlanMany(&
       hip_plan_c2r_manyG, &
       irank, &
       ndim, &
       inembed, &
       istride, &
       idist, &
       onembed, &
       ostride, &
       odist, &
       merge(HIPFFT_C2R,HIPFFT_Z2D,kind(uxmany) == singlePrecision), &
       nsplit)
#elif defined(MKLGPU)
     dfftw_plan_c2r_manyA = 0
!$omp target data map(tofrom: fymany,uymany)
     !$omp dispatch
     call dfftw_plan_many_dft_c2r(&
          dfftw_plan_c2r_manyA, &
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
     dfftw_plan_c2r_manyB = 0
     !$omp dispatch
     call dfftw_plan_many_dft_c2r(&
          dfftw_plan_c2r_manyB, &
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

     dfftw_plan_c2r_manyG = 0
!$omp target data map(tofrom: gymany,vymany)
     !$omp dispatch
     call dfftw_plan_many_dft_c2r(&
          dfftw_plan_c2r_manyG, &
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
       cu_plan_c2r_manyA, &
       irank, &
       ndim, &
       inembed, &
       istride, &
       idist, &
       onembed, &
       ostride, &
       odist, &
       merge(CUFFT_C2R,CUFFT_Z2D,kind(uxmany) == singlePrecision), &
       nsplitA)

  if (nsplitB > 0) then ! no fft if nsplitB==0
    istatus = cufftPlanMany(&
       cu_plan_c2r_manyB, &
       irank, &
       ndim, &
       inembed, &
       istride, &
       idist, &
       onembed, &
       ostride, &
       odist, &
       merge(CUFFT_C2R,CUFFT_Z2D,kind(uxmany) == singlePrecision), &
       nsplitB)
  endif

  istatus = cufftPlanMany(&
       cu_plan_c2r_manyG, &
       irank, &
       ndim, &
       inembed, &
       istride, &
       idist, &
       onembed, &
       ostride, &
       odist, &
       merge(CUFFT_C2R,CUFFT_Z2D,kind(uxmany) == singlePrecision), &
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
  hip_plan_r2c_manyA = c_null_ptr
  istatus = hipfftPlanMany(&
       hip_plan_r2c_manyA, &
       irank, &
       ndim, &
       inembed, &
       istride, &
       idist, &
       onembed, &
       ostride, &
       odist, &
       merge(HIPFFT_R2C,HIPFFT_D2Z,kind(uxmany) == singlePrecision), &
       nsplitA)

  if (nsplitB > 0) then ! no fft if nsplitB==0
    hip_plan_r2c_manyB = c_null_ptr
    istatus = hipfftPlanMany(&
       hip_plan_r2c_manyB, &
       irank, &
       ndim, &
       inembed, &
       istride, &
       idist, &
       onembed, &
       ostride, &
       odist, &
       merge(HIPFFT_R2C,HIPFFT_D2Z,kind(uxmany) == singlePrecision), &
       nsplitB)
  endif
#elif defined(MKLGPU)
     dfftw_plan_r2c_manyA = 0
!$omp target data map(tofrom: uvmany,fxmany)
     !$omp dispatch
     call dfftw_plan_many_dft_r2c(&
          dfftw_plan_r2c_manyA, &
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
     dfftw_plan_r2c_manyB = 0
     !$omp dispatch
     call dfftw_plan_many_dft_r2c(&
          dfftw_plan_r2c_manyB, &
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
       cu_plan_r2c_manyA, &
       irank, &
       ndim, &
       inembed, &
       istride, &
       idist, &
       onembed, &
       ostride, &
       odist, &
       merge(CUFFT_R2C,CUFFT_D2Z,kind(uxmany) == singlePrecision), &
       nsplitA)

  if (nsplitB > 0) then ! no fft if nsplitB==0
    istatus = cufftPlanMany(&
       cu_plan_r2c_manyB, &
       irank, &
       ndim, &
       inembed, &
       istride, &
       idist, &
       onembed, &
       ostride, &
       odist, &
       merge(CUFFT_R2C,CUFFT_D2Z,kind(uxmany) == singlePrecision), &
       nsplitB)
  endif
#endif

end subroutine cgyro_nl_fftw_init

subroutine cgyro_nl_fftw_mul(sz,uvm,uxm,vym,uym,vxm,inv_nxny)
  implicit none

  integer, intent(in) :: sz
  real, dimension(*),intent(out) :: uvm
  real, dimension(*),intent(in) :: uxm,vym,uym,vxm
  real, intent(in) :: inv_nxny

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

end subroutine

subroutine cgyro_nl_fftw

#if defined(HIPGPU)
  use hipfort
  use hipfort_hipfft
#elif defined(MKLGPU)
  use fftw3_omp_offload
#else
  use cufft
#endif
  use timer_lib
  use parallel_lib
  use cgyro_nl_comm
  use cgyro_globals

  implicit none
#if defined(MKLGPU)
  include 'fftw/fftw3.f'
#endif
  !-----------------------------------
  integer :: j,p,iexch
  integer :: it,ir,itm,itl,ix,iy
  integer :: itor,mytm
  integer :: i1,i2
  integer :: it_loc
  integer :: ierr
  integer :: rc
  complex :: f0,g0
  integer :: jtheta_min
  integer :: iy0, iy1, ir0, ir1

  real :: inv_nxny

#ifdef GACODE_GPU_AMD
  ! AMD GPU  (MI250X) optimal
  integer, parameter :: F_RADTILE = 8
  integer, parameter :: F_TORTILE = 16
  integer, parameter :: R_RADTILE = 32
  integer, parameter :: R_TORTILE = 8
#else
  ! NVIDIA GPU  (A100) optimal
  integer, parameter :: F_RADTILE = 8
  integer, parameter :: F_TORTILE = 16
  integer, parameter :: R_RADTILE = 16
  integer, parameter :: R_TORTILE = 8
#endif

  call timer_lib_in('nl_mem')
  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()

  ! we can zero the elements we know are zero while we wait
#if defined(OMPGPU)
  !no async for OMPGPU for now
!$omp target teams distribute parallel do simd collapse(3)
#else
!$acc parallel loop gang vector independent collapse(3) async(2) &
!$acc&         private(j,ix,iy) &
!$acc&         present(nsplitA,ny2,nx0,nx2)
#endif
  do j=1,nsplitA
     do ix=nx2,nx0-1
       do iy=0,ny2
         fxmany(iy,ix,j) = 0
         fymany(iy,ix,j) = 0
       enddo
     enddo
  enddo

#if defined(OMPGPU)
  !no async for OMPGPU for now
!$omp target teams distribute parallel do simd collapse(3)
#else
!$acc parallel loop gang vector independent collapse(3) async(2) &
!$acc&         private(j,ix,iy) &
!$acc&         present(nsplitA,ny2,n_toroidal,nx)
#endif
  do j=1,nsplitA
     do ix=0,nx-1
       do iy=n_toroidal,ny2
         fxmany(iy,ix,j) = 0
         fymany(iy,ix,j) = 0
       enddo
     enddo
  enddo
  call timer_lib_out('nl_mem')

  ! time to wait for the FA_nl to become avaialble
  call timer_lib_in('nl_comm')
  call parallel_slib_f_nc_wait(nsplitA,fpackA,fA_nl,fA_req)
  fA_req_valid = .FALSE.
  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()
  call timer_lib_out('nl_comm')

  call timer_lib_in('nl')
#if !defined(OMPGPU)
!$acc  data present(fA_nl)  &
!$acc&      present(fxmany,fymany,gxmany,gymany) &
!$acc&      present(uxmany,uymany,vxmany,vymany) &
!$acc&      present(uvmany)
#endif

! f_nl is (radial, nt_loc, theta, nv_loc1, toroidal_procs)
! where nv_loc1 * toroidal_procs >= nv_loc

  ! no tiling, does not seem to help
#if defined(OMPGPU)
  !no async for OMPGPU for now
!$omp target teams distribute parallel do simd collapse(4) &
!$omp&  private(j,ir,p,ix,itor,iy,f0,itm,itl)
#else
!$acc parallel loop gang vector independent collapse(4) async(2) &
!$acc&         private(j,ir,p,ix,itor,iy,f0,itm,itl)
#endif
  do j=1,nsplitA
     do ir=1,n_radial
       do itm=1,n_toroidal_procs
         do itl=1,nt_loc
              itor=itl + (itm-1)*nt_loc
              iy = itor-1
              p  = ir-1-nx0/2
              ix = p
              if (ix < 0) ix = ix+nx

              f0 = i_c*fA_nl(ir,itl,j,itm)
              fxmany(iy,ix,j) = p*f0
              fymany(iy,ix,j) = iy*f0
         enddo
       enddo
     enddo
  enddo

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
#if defined(OMPGPU)
  !no async for OMPGPU for now
!$omp target teams distribute parallel do simd collapse(2) &
!$omp&  private(j,ix,f0)
#else
!$acc parallel loop gang vector independent collapse(2) async(2) &
!$acc&         private(j,ix,f0) &
!$acc&         present(fxmany) &
!$acc&         present(nsplitA,nx)
#endif
  do j=1,nsplitA
    do ix=1,nx/2-1
      f0 = 0.5*( fxmany(0,ix,j)+conjg(fxmany(0,nx-ix,j)) )
      fxmany(0,ix   ,j) = f0
      fxmany(0,nx-ix,j) = conjg(f0)
    enddo
  enddo

     ! --------------------------------------
     ! perform many Fourier Transforms at once
     ! --------------------------------------

#if defined(OMPGPU)

#if defined(MKLGPU)
!$omp target data map(tofrom: fymany,uymany)
#else
!$omp target data use_device_ptr(fymany,uymany)
#endif

#else
!$acc  host_data use_device(fymany,uymany)
#endif

#if defined(HIPGPU)
  rc = hipfftExecZ2D(hip_plan_c2r_manyA,c_loc(fymany),c_loc(uymany))
#elif defined(MKLGPU)
  !$omp dispatch
  call dfftw_execute_dft_c2r(dfftw_plan_c2r_manyA,fymany,uymany)
  rc = 0
#else
  rc = cufftExecZ2D(cu_plan_c2r_manyA,fymany,uymany)
#endif

#if defined(OMPGPU)
!$omp end target data
#else
!$acc end host_data
#endif

#if defined(OMPGPU)
  !no async for OMPGPU for now
#else
  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()
!$acc wait(2)
#endif
  ! fxmany is complete now

  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()

#if defined(OMPGPU)

#if defined(MKLGPU)
!$omp target data map(tofrom: fxmany,uxmany)
#else
!$omp target data use_device_ptr(fxmany,uxmany)
#endif

#else
!$acc  host_data use_device(fxmany,uxmany)
#endif

#if defined(HIPGPU)
  rc = hipfftExecZ2D(hip_plan_c2r_manyA,c_loc(fxmany),c_loc(uxmany))
#elif defined(MKLGPU)
  !$omp dispatch
  call dfftw_execute_dft_c2r(dfftw_plan_c2r_manyA,fxmany,uxmany)
  rc = 0
#else
  rc = cufftExecZ2D(cu_plan_c2r_manyA,fxmany,uxmany)
#endif

#if defined(OMPGPU)
!$omp end target data
#else
!$acc end host_data
#endif

  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()

#if !defined(OMPGPU)
!$acc  data present(g_nl)  &
!$acc&      present(gxmany,gymany)
#endif
  ! we can zero the elements we know are zero while we wait for comm
#if defined(OMPGPU)
  !no async for OMPGPU for now
!$omp target teams distribute parallel do simd collapse(3)
#else
!$acc parallel loop gang vector independent collapse(3) async(2) &
!$acc&         private(j,ix,iy) &
!$acc&         present(nsplit,ny2,nx0,nx2)
#endif
  do j=1,nsplit
     do ix=nx2,nx0-1
       do iy=0,ny2
         gxmany(iy,ix,j) = 0
         gymany(iy,ix,j) = 0
       enddo
     enddo
  enddo

#if defined(OMPGPU)
  !no async for OMPGPU for now
!$omp target teams distribute parallel do simd collapse(3)
#else
!$acc parallel loop gang vector independent collapse(3) async(2) &
!$acc&         private(j,ix,iy) &
!$acc&         present(nsplit,ny2,n_toroidal,nx)
#endif
  do j=1,nsplit
     do ix=0,nx-1
       do iy=n_toroidal,ny2
         gxmany(iy,ix,j) = 0
         gymany(iy,ix,j) = 0
       enddo
     enddo
  enddo
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

  ! tile for performance, since this is effectively a transpose
#if defined(OMPGPU)
  !no async for OMPGPU for now
!$omp target teams distribute parallel do simd collapse(5) &
!$omp&   private(j,p,ix,itor,mytm,iy,g0,it,iv_loc,it_loc,jtheta_min,itm,itl,ir)
#else
!$acc parallel loop gang vector independent collapse(5) async(2) &
!$acc&         private(j,p,ix,itor,mytm,iy,g0,it,iv_loc,it_loc,jtheta_min,itm,itl,ir) &
!$acc&         present(g_nl,jvec_c_nl) &
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
#if defined(OMPGPU)
  !no async for OMPGPU for now
!$omp target teams distribute parallel do simd collapse(2) &
!$omp&   private(j,ix,g0)
#else
!$acc parallel loop gang vector independent collapse(2) async(2) &
!$acc&         private(j,ix,g0) &
!$acc&         present(gxmany) &
!$acc&         present(nsplit,nx)
#endif
  do j=1,nsplit
    do ix=1,nx/2-1
      g0 = 0.5*( gxmany(0,ix,j)+conjg(gxmany(0,nx-ix,j)) )
      gxmany(0,ix   ,j) = g0
      gxmany(0,nx-ix,j) = conjg(g0)
    enddo
  enddo

#if !defined(OMPGPU)
!$acc end data
#endif

#if defined(OMPGPU)

#if defined(MKLGPU)
!$omp target data map(tofrom: gymany,vymany)
#else
!$omp target data use_device_ptr(gymany,vymany)
#endif

#else
!$acc  host_data use_device(gymany,vymany)
#endif

#if defined(HIPGPU)
  rc = hipfftExecZ2D(hip_plan_c2r_manyG,c_loc(gymany),c_loc(vymany))
#elif defined(MKLGPU)
  !$omp dispatch
  call dfftw_execute_dft_c2r(dfftw_plan_c2r_manyG,gymany,vymany)
  rc = 0
#else
  rc = cufftExecZ2D(cu_plan_c2r_manyG,gymany,vymany)
#endif

#if defined(OMPGPU)
!$omp end target data
#else
!$acc end host_data
#endif

#if defined(OMPGPU)
  !no async for OMPGPU for now
#else
  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()
!$acc wait(2)
#endif

  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()
  ! gxmany is complete now

#if defined(OMPGPU)

#if defined(MKLGPU)
!$omp target data map(tofrom: gxmany,vxmany)
#else
!$omp target data use_device_ptr(gxmany,vxmany)
#endif

#else
!$acc  host_data use_device(gxmany,vxmany)
#endif

#if defined(HIPGPU)
  rc = hipfftExecZ2D(hip_plan_c2r_manyG,c_loc(gxmany),c_loc(vxmany))
#elif defined(MKLGPU)
  !$omp dispatch
  call dfftw_execute_dft_c2r(dfftw_plan_c2r_manyG,gxmany,vxmany)
  rc = 0
#else
  rc = cufftExecZ2D(cu_plan_c2r_manyG,gxmany,vxmany)
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

  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()

  ! Poisson bracket in real space
  ! uv = (ux*vy-uy*vx)/(nx*ny)

  inv_nxny = 1.0/(nx*ny)

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

#if defined(OMPGPU)

#if defined(MKLGPU)
!$omp target data map(tofrom: uvmany,fxmany)
#else
!$omp target data use_device_ptr(uvmany,fxmany)
#endif

#else
!$acc wait
!$acc  host_data use_device(uvmany,fxmany)
#endif

#if defined(HIPGPU)
  rc = hipfftExecD2Z(hip_plan_r2c_manyA,c_loc(uvmany),c_loc(fxmany))
#elif defined(MKLGPU)
  !$omp dispatch
  call dfftw_execute_dft_r2c(dfftw_plan_r2c_manyA,uvmany,fxmany)
  rc = 0
#else
  rc = cufftExecD2Z(cu_plan_r2c_manyA,uvmany,fxmany)
#endif

#ifdef HIPGPU
  ! hipfftExec is asynchronous, will need the results below
  rc = hipDeviceSynchronize()
#endif

#if defined(OMPGPU)
!$omp end target data
#else
!$acc wait
!$acc end host_data
#endif

  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()

  call timer_lib_out('nl')
  call timer_lib_in('nl_mem')

  ! NOTE: The FFT will generate an unwanted n=0,p=-nr/2 component
  ! that will be filtered in the main time-stepping loop

  ! tile for performance, since this is effectively a transpose
#if defined(OMPGPU)
!$omp target teams distribute parallel do collapse(5) &
!$omp&   private(iy,ir,itm,itl,ix)
#else
!$acc parallel loop independent collapse(5) gang &
!$acc&         private(iy,ir,itm,itl,ix) present(fA_nl,fxmany)
#endif
  do j=1,nsplitA
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

           fA_nl(ir,itl,j,itm) = fxmany(iy,ix,j)
        endif
      enddo
     enddo
    enddo
   enddo
  enddo

#if !defined(OMPGPU)
  ! end data fA_nl
!$acc end data
#endif

  if (nsplitB > 0) then
  ! we can zero the elements we know are zero while we waita
  ! assuming nsplitB<=nsplitA
#if defined(OMPGPU)
  !no async for OMPGPU for now
!$omp target teams distribute parallel do simd collapse(3)
#else
!$acc parallel loop gang vector independent collapse(3) async(2) &
!$acc&         private(j,ix,iy) &
!$acc&         present(nsplitB,ny2,nx0,nx2)
#endif
   do j=1,nsplitB
     do ix=nx2,nx0-1
       do iy=0,ny2
         fxmany(iy,ix,j) = 0
         fymany(iy,ix,j) = 0
       enddo
     enddo
   enddo

#if defined(OMPGPU)
  !no async for OMPGPU for now
!$omp target teams distribute parallel do simd collapse(3)
#else
!$acc parallel loop gang vector independent collapse(3) async(2) &
!$acc&         private(j,ix,iy) &
!$acc&         present(nsplitB,ny2,n_toroidal,nx)
#endif
   do j=1,nsplitB
     do ix=0,nx-1
       do iy=n_toroidal,ny2
         fxmany(iy,ix,j) = 0
         fymany(iy,ix,j) = 0
       enddo
     enddo
   enddo
  endif ! if nsplitB>0

  call timer_lib_out('nl_mem')

  call timer_lib_in('nl_comm')
  ! start the async reverse comm
  ! can reuse the same req, no overlap with forward fA_req
  call parallel_slib_r_nc_async(nsplitA,fA_nl,fpackA,fA_req)
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
#if !defined(OMPGPU)
!$acc  data present(fB_nl)  &
!$acc&      present(fxmany,fymany,gxmany,gymany) &
!$acc&      present(uxmany,uymany,vxmany,vymany) &
!$acc&      present(uvmany)
#endif

! f_nl is (radial, nt_loc, theta, nv_loc1, toroidal_procs)
! where nv_loc1 * toroidal_procs >= nv_loc

  ! no tiling, does not seem to help
#if defined(OMPGPU)
  !no async for OMPGPU for now
!$omp target teams distribute parallel do simd collapse(4) &
!$omp&  private(j,ir,p,ix,itor,iy,f0,g0,itm,itl)
#else
!$acc parallel loop gang vector independent collapse(4) async(2) &
!$acc&         private(j,ir,p,ix,itor,iy,f0,itm,itl)
#endif
  do j=1,nsplitB
     do ir=1,n_radial
       do itm=1,n_toroidal_procs
         do itl=1,nt_loc
              itor=itl + (itm-1)*nt_loc
              iy = itor-1
              p  = ir-1-nx0/2
              ix = p
              if (ix < 0) ix = ix+nx

              f0 = i_c*fB_nl(ir,itl,j,itm)
              fxmany(iy,ix,j) = p*f0
              fymany(iy,ix,j) = iy*f0
         enddo
       enddo
     enddo
  enddo

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
#if defined(OMPGPU)
  !no async for OMPGPU for now
!$omp target teams distribute parallel do simd collapse(2) &
!$omp&  private(j,ix,f0)
#else
!$acc parallel loop gang vector independent collapse(2) async(2) &
!$acc&         private(j,ix,f0) &
!$acc&         present(fxmany) &
!$acc&         present(nsplitB,nx)
#endif
  do j=1,nsplitB
    do ix=1,nx/2-1
      f0 = 0.5*( fxmany(0,ix,j)+conjg(fxmany(0,nx-ix,j)) )
      fxmany(0,ix   ,j) = f0
      fxmany(0,nx-ix,j) = conjg(f0)
    enddo
  enddo

     ! --------------------------------------
     ! perform many Fourier Transforms at once
     ! --------------------------------------

#if defined(OMPGPU)

#if defined(MKLGPU)
!$omp target data map(tofrom: fymany,uymany)
#else
!$omp target data use_device_ptr(fymany,uymany)
#endif

#else
!$acc  host_data use_device(fymany,uymany)
#endif

#if defined(HIPGPU)
  rc = hipfftExecZ2D(hip_plan_c2r_manyB,c_loc(fymany),c_loc(uymany))
#elif defined(MKLGPU)
  !$omp dispatch
  call dfftw_execute_dft_c2r(dfftw_plan_c2r_manyB,fymany,uymany)
  rc = 0
#else
  rc = cufftExecZ2D(cu_plan_c2r_manyB,fymany,uymany)
#endif

#if defined(OMPGPU)
!$omp end target data
#else
!$acc end host_data
#endif

#if defined(OMPGPU)
  !no async for OMPGPU for now
#else
  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()
!$acc wait(2)
#endif
  ! fxmany is complete now

  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()

#if defined(OMPGPU)

#if defined(MKLGPU)
!$omp target data map(tofrom: fxmany,uxmany)
#else
!$omp target data use_device_ptr(fxmany,uxmany)
#endif

#else
!$acc  host_data use_device(fxmany,uxmany)
#endif

#if defined(HIPGPU)
  rc = hipfftExecZ2D(hip_plan_c2r_manyB,c_loc(fxmany),c_loc(uxmany))
#elif defined(MKLGPU)
  !$omp dispatch
  call dfftw_execute_dft_c2r(dfftw_plan_c2r_manyB,fxmany,uxmany)
  rc = 0
#else
  rc = cufftExecZ2D(cu_plan_c2r_manyB,fxmany,uxmany)
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

  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()

  ! Poisson bracket in real space
  ! uv = (ux*vy-uy*vx)/(nx*ny)

  inv_nxny = 1.0/(nx*ny)

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

#if defined(OMPGPU)

#if defined(MKLGPU)
!$omp target data map(tofrom: uvmany,fxmany)
#else
!$omp target data use_device_ptr(uvmany,fxmany)
#endif

#else
!$acc wait
!$acc  host_data use_device(uvmany,fxmany)
#endif

#if defined(HIPGPU)
  rc = hipfftExecD2Z(hip_plan_r2c_manyB,c_loc(uvmany),c_loc(fxmany))
#elif defined(MKLGPU)
  !$omp dispatch
  call dfftw_execute_dft_r2c(dfftw_plan_r2c_manyB,uvmany,fxmany)
  rc = 0
#else
  rc = cufftExecD2Z(cu_plan_r2c_manyB,uvmany,fxmany)
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

  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()

  call timer_lib_out('nl')
  call timer_lib_in('nl_mem')

  ! NOTE: The FFT will generate an unwanted n=0,p=-nr/2 component
  ! that will be filtered in the main time-stepping loop

  ! tile for performance, since this is effectively a transpose
#if defined(OMPGPU)
!$omp target teams distribute parallel do collapse(5) &
!$omp&   private(iy,ir,itm,itl,ix)
#else
!$acc parallel loop independent collapse(5) gang &
!$acc&         private(iy,ir,itm,itl,ix) present(fB_nl,fxmany)
#endif
  do j=1,nsplitB
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

           fB_nl(ir,itl,j,itm) = fxmany(iy,ix,j)
        endif
      enddo
     enddo
    enddo
   enddo
  enddo

#if !defined(OMPGPU)
  ! end data fB_nl
!$acc end data
#endif

  call timer_lib_out('nl_mem')

  call timer_lib_in('nl_comm')
  ! start the async reverse comm
  ! can reuse the same req, no overlap with forward fB_req
  call parallel_slib_r_nc_async(nsplitB,fB_nl,fpackB,fB_req)
  fB_req_valid = .TRUE.
  ! make sure reqs progress
  call cgyro_nl_fftw_comm_test()
  call timer_lib_out('nl_comm')

  endif ! nsplitB>0

end subroutine cgyro_nl_fftw
