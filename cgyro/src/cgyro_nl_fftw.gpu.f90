!-----------------------------------------------------------------
! cgyro_nl_fftw.gpu.f90 [GPU (acc-cuFFT) version]
!
! PURPOSE:
!  Evaluate nonlinear bracket with dealiased FFT.  It is natural 
!  to use the FFTW complex-to-real (c2r) transform to move to real 
!  space, compute the convolution (uv), and transform back to the 
!  spectral form using a real-to-complex transform (r2c).
!  
! NOTE: Need to be careful with (p=-nr/2,n=0) component.
!-----------------------------------------------------------------

! NOTE: call cgyro_nl_fftw_comm1 before cgyro_nl_fftw
subroutine cgyro_nl_fftw_comm1
  use timer_lib
  use parallel_lib
  use cgyro_globals

  implicit none

  integer :: ir,it,iv_loc_m,ic_loc_m
  integer :: iexch

!$omp parallel do private(it,ir,iexch,ic_loc_m,iv_loc_m)
  do iv_loc_m=1,nv_loc
     do it=1,n_theta
        iexch = it + (iv_loc_m-1)*n_theta
        do ir=1,n_radial
           ic_loc_m = it + (ir-1)*n_theta
           fpack(ir,iexch) = h_x(ic_loc_m,iv_loc_m)
        enddo
     enddo
  enddo

  do iexch=nv_loc*n_theta+1,nsplit*n_toroidal
     fpack(:,iexch) = (0.0,0.0)
  enddo

  call parallel_slib_f_nc(fpack,f_nl)
end subroutine cgyro_nl_fftw_comm1

subroutine cgyro_nl_fftw_comm2
  use timer_lib
  use parallel_lib
  use cgyro_globals

  implicit none

  integer :: ir,it,iv_loc_m,ic_loc_m
  integer :: iexch

!$omp parallel do private(it,ir,iexch,ic_loc_m,iv_loc_m)
  do iv_loc_m=1,nv_loc
     do it=1,n_theta
        iexch = it + (iv_loc_m-1)*n_theta
        do ir=1,n_radial
           ic_loc_m = it + (ir-1)*n_theta
           gpack(ir,iexch) = psi(ic_loc_m,iv_loc_m)
        enddo
     enddo
  enddo

  do iexch=nv_loc*n_theta+1,nsplit*n_toroidal
     gpack(:,iexch) = (0.0,0.0)
  enddo

  call parallel_slib_f_nc(gpack,g_nl)

end subroutine cgyro_nl_fftw_comm2


subroutine cgyro_nl_fftw(ij)

  use cufft
  use timer_lib
  use parallel_lib

  use cgyro_globals
  implicit none
  include 'mpif.h'
  integer, intent(in) :: ij
  integer :: j,p,iexch
  integer :: it,ir,in,ix,iy
  integer :: i1,i2
  integer :: ierr
  integer :: rc
  complex :: f0,g0

  real :: inv_nxny

  include 'fftw3.f03'


  if (is_staggered_comm_2) then ! stagger comm2, to load ballance network traffic
     call timer_lib_in('nl_comm')
     call cgyro_nl_fftw_comm2
     call timer_lib_out('nl_comm')
  endif

  call timer_lib_in('nl')
!$acc  data pcopyin(f_nl)   &
!$acc& pcreate(fxmany,fymany,gxmany,gymany) &
!$acc& pcreate(uxmany,uymany,vxmany,vymany) &
!$acc& pcreate(uvmany)

!$acc parallel
!$acc loop gang
  do  j=lbound(fxmany,3),ubound(fxmany,3)
!$acc loop worker
     do ix=lbound(fxmany,2),ubound(fxmany,2)
!$acc loop vector
        do iy=lbound(fxmany,1),ubound(fxmany,1)
           fxmany(iy,ix,j) = 0.0
           fymany(iy,ix,j) = 0.0
           gxmany(iy,ix,j) = 0.0
           gymany(iy,ix,j) = 0.0
        enddo
     enddo
  enddo
!$acc end parallel

!$acc parallel 
!$acc loop gang
  do j=1,nsplit

   ! Array mapping
!$acc loop worker private(p,ix)
     do ir=1,n_radial

        p  = ir-1-nx0/2
        ix = p
        if (ix < 0) ix = ix+nx  
!$acc   loop vector private(iy,f0,g0)
        do in=1,n_toroidal
           iy = in-1
           f0 = i_c*f_nl(ir,j,in)
           fxmany(iy,ix,j) = p*f0
           fymany(iy,ix,j) = iy*f0
        enddo
     enddo
  enddo
!$acc end parallel

     ! --------------------------------------
     ! perform many Fourier Transforms at once
     ! --------------------------------------
!$acc wait
!$acc  host_data &
!$acc& use_device(fxmany,fymany) &
!$acc& use_device(uxmany,uymany)

  rc = cufftExecZ2D(cu_plan_c2r_many,fxmany,uxmany)
  rc = cufftExecZ2D(cu_plan_c2r_many,fymany,uymany)

!$acc wait
!$acc end host_data
  if (.not. is_staggered_comm_2) then ! stagger comm2, to load ballance network traffic
     call timer_lib_out('nl')
     call timer_lib_in('nl_comm')
     call cgyro_nl_fftw_comm2
     call timer_lib_out('nl_comm')
     call timer_lib_in('nl')
  endif

!$acc data copyin(g_nl)  

!$acc parallel 
!$acc loop gang
  do j=1,nsplit

   ! Array mapping
!$acc loop worker private(p,ix)
     do ir=1,n_radial

        p  = ir-1-nx0/2
        ix = p
        if (ix < 0) ix = ix+nx
!$acc   loop vector private(iy,f0,g0)
        do in=1,n_toroidal
           iy = in-1
           g0 = i_c*g_nl(ir,j,in)
           gxmany(iy,ix,j) = p*g0
           gymany(iy,ix,j) = iy*g0
        enddo
     enddo
  enddo
!$acc end parallel

!$acc end data

!$acc wait
!$acc  host_data &
!$acc& use_device(gxmany,gymany) &
!$acc& use_device(vxmany,vymany)

  rc = cufftExecZ2D(cu_plan_c2r_many,gxmany,vxmany)
  rc = cufftExecZ2D(cu_plan_c2r_many,gymany,vymany)

!$acc wait
!$acc end host_data

!$acc wait

  ! Poisson bracket in real space
  ! uv = (ux*vy-uy*vx)/(nx*ny)

  inv_nxny = 1.0/(nx*ny)

!$acc  parallel 
!$acc loop gang
  do j=1,nsplit
!$acc loop worker
     do ix=lbound(uvmany,2),ubound(uvmany,2)
!$acc loop vector
        do iy=lbound(uvmany,1),ubound(uvmany,1)
           uvmany(iy,ix,j) = (uxmany(iy,ix,j)*vymany(iy,ix,j)- &
                uymany(iy,ix,j)*vxmany(iy,ix,j))*inv_nxny
        enddo
     enddo
  enddo
!$acc  end parallel

  ! ------------------
  ! Transform uv to fx
  ! ------------------

!$acc wait
!$acc host_data use_device(uvmany,fxmany)
  rc = cufftExecD2Z(cu_plan_r2c_many,uvmany,fxmany)
!$acc wait
!$acc end host_data
!$acc wait

  ! NOTE: The FFT will generate an unwanted n=0,p=-nr/2 component
  ! that will be filtered in the main time-stepping loop

!$acc parallel  
!$acc loop gang
  do j=1,nsplit
!$acc loop worker private(ix)
     do ir=1,n_radial 
        ix = ir-1-nx0/2
        if (ix < 0) ix = ix+nx
!$acc loop vector private(iy)
        do in=1,n_toroidal
           iy = in-1
           g_nl(ir,j,in) = fxmany(iy,ix,j)
        enddo
     enddo
  enddo
!$acc end parallel

!$acc wait
!$acc end data
!$acc wait

  call timer_lib_out('nl')

  call timer_lib_in('nl_comm')
  call parallel_slib_r_nc(g_nl,gpack)

!$omp parallel do private(it,ir,iexch,ic_loc)
  do iv_loc=1,nv_loc
     do it=1,n_theta
        iexch = it + (iv_loc-1)*n_theta
        do ir=1,n_radial
           ic_loc = it + (ir-1)*n_theta
           psi(ic_loc,iv_loc) = gpack(ir,iexch) 
        enddo
     enddo
  enddo

  call timer_lib_out('nl_comm')

  ! RHS -> -[f,g] = [f,g]_{r,-alpha}

!$omp workshare 
  rhs(:,:,ij) = rhs(:,:,ij)+((q*rho/rmin)*(2*pi/length))*psi(:,:)
!$omp end workshare

end subroutine cgyro_nl_fftw
