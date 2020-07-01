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

  call timer_lib_in('nl_mem')

!$acc parallel loop collapse(3) independent private(iexch,ic_loc_m) &
!$acc&         present(h_x,fpack) default(none)
  do iv_loc_m=1,nv_loc
     do it=1,n_theta
        do ir=1,n_radial
           iexch = it + (iv_loc_m-1)*n_theta
           ic_loc_m = it + (ir-1)*n_theta
           fpack(ir,iexch) = h_x(ic_loc_m,iv_loc_m)
        enddo
     enddo
  enddo

!$acc parallel loop gang present(fpack) if(nv_loc*n_theta+1>=nsplit*n_toroidal)
  do iexch=nv_loc*n_theta+1,nsplit*n_toroidal
     fpack(:,iexch) = (0.0,0.0)
  enddo

  call timer_lib_out('nl_mem')
  call timer_lib_in('nl_comm')

  call parallel_slib_f_nc_gpu(fpack,f_nl)

  call timer_lib_out('nl_comm')

end subroutine cgyro_nl_fftw_comm1

subroutine cgyro_nl_fftw_comm2
  use timer_lib
  use parallel_lib
  use cgyro_globals

  implicit none

  integer :: ir,it,iv_loc_m,ic_loc_m
  integer :: iexch

  call timer_lib_in('nl_mem')

!$acc parallel loop collapse(3) independent private(iexch,ic_loc_m) &
!$acc&         present(psi,gpack) default(none)
  do iv_loc_m=1,nv_loc
     do it=1,n_theta
        do ir=1,n_radial
           iexch = it + (iv_loc_m-1)*n_theta
           ic_loc_m = it + (ir-1)*n_theta
           gpack(ir,iexch) = psi(ic_loc_m,iv_loc_m)
        enddo
     enddo
  enddo

!$acc parallel loop gang present(gpack) if(nv_loc*n_theta+1>=nsplit*n_toroidal)
  do iexch=nv_loc*n_theta+1,nsplit*n_toroidal
     gpack(:,iexch) = (0.0,0.0)
  enddo

  call timer_lib_out('nl_mem')
  call timer_lib_in('nl_comm')

  call parallel_slib_f_nc_gpu(gpack,g_nl)

  call timer_lib_out('nl_comm')

end subroutine cgyro_nl_fftw_comm2

subroutine cgyro_nl_fftw_zero4(sz,v1,v2,v3,v4)
  implicit none

  integer, intent(in) :: sz
  complex, dimension(*), intent(inout) :: v1,v2,v3,v4

  integer :: i

!$acc parallel loop independent present(v1,v2,v3,v4) private(i)
  do i=1,sz
    v1(i) = 0.0
    v2(i) = 0.0
    v3(i) = 0.0
    v4(i) = 0.0
  enddo

end subroutine

subroutine cgyro_nl_fftw_mul(sz,uvm,uxm,vym,uym,vxm,inv_nxny)
  implicit none

  integer, intent(in) :: sz
  real, dimension(*),intent(out) :: uvm
  real, dimension(*),intent(in) :: uxm,vym,uym,vxm
  real, intent(in) :: inv_nxny

  integer :: i

!$acc parallel loop independent present(uvm,uxm,vym,uym,vxm) private(i)
  do i=1,sz
    uvm(i) = (uxm(i)*vym(i)-uym(i)*vxm(i))*inv_nxny
  enddo

end subroutine

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


  if (is_staggered_comm_2) then ! stagger comm2, to load ballance network traffic
     call cgyro_nl_fftw_comm2
  endif

  call timer_lib_in('nl_mem')
!$acc  data present(f_nl)  &
!$acc&      present(fxmany,fymany,gxmany,gymany) &
!$acc&      present(uxmany,uymany,vxmany,vymany) &
!$acc&      present(uvmany)

  call timer_lib_out('nl_mem')

  call timer_lib_in('nl')

  call cgyro_nl_fftw_zero4(size(fxmany,1)*size(fxmany,2)*size(fxmany,3), &
                           fxmany,fymany,gxmany,gymany)

!$acc parallel loop independent collapse(3) private(j,ir,p,ix,in,iy,f0,g0)
  do j=1,nsplit
     do ir=1,n_radial
        do in=1,n_toroidal
           p  = ir-1-nx0/2
           ix = p
           if (ix < 0) ix = ix+nx

           iy = in-1
           f0 = i_c*f_nl(ir,j,in)
           fxmany(iy,ix,j) = p*f0
           fymany(iy,ix,j) = iy*f0
        enddo
     enddo
  enddo

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
  call timer_lib_out('nl')

  if (.not. is_staggered_comm_2) then ! stagger comm2, to load ballance network traffic
     call cgyro_nl_fftw_comm2
  endif

  call timer_lib_in('nl_mem')
!$acc data present(g_nl)  

  call timer_lib_out('nl_mem')
  call timer_lib_in('nl')


!$acc parallel loop independent collapse(3) private(j,ir,p,ix,in,iy,f0,g0)
  do j=1,nsplit
     do ir=1,n_radial
        do in=1,n_toroidal
           p  = ir-1-nx0/2
           ix = p
           if (ix < 0) ix = ix+nx

           iy = in-1
           g0 = i_c*g_nl(ir,j,in)
           gxmany(iy,ix,j) = p*g0
           gymany(iy,ix,j) = iy*g0
        enddo
     enddo
  enddo

!$acc  host_data &
!$acc& use_device(gxmany,gymany) &
!$acc& use_device(vxmany,vymany)

  rc = cufftExecZ2D(cu_plan_c2r_many,gxmany,vxmany)
  rc = cufftExecZ2D(cu_plan_c2r_many,gymany,vymany)

!$acc wait
!$acc end host_data

  ! Poisson bracket in real space
  ! uv = (ux*vy-uy*vx)/(nx*ny)

  inv_nxny = 1.0/(nx*ny)

  call cgyro_nl_fftw_mul(size(uvmany,1)*size(uvmany,2)*size(uvmany,3), &
                         uvmany,uxmany,vymany,uymany,vxmany,inv_nxny)

  ! ------------------
  ! Transform uv to fx
  ! ------------------

!$acc wait
!$acc host_data use_device(uvmany,fxmany)
  rc = cufftExecD2Z(cu_plan_r2c_many,uvmany,fxmany)
!$acc wait
!$acc end host_data
!$acc wait

  call timer_lib_out('nl')
  call timer_lib_in('nl_mem')

  ! NOTE: The FFT will generate an unwanted n=0,p=-nr/2 component
  ! that will be filtered in the main time-stepping loop

!$acc parallel loop independent collapse(3) private(j,ir,in,ix,iy)
  do j=1,nsplit
     do ir=1,n_radial 
        do in=1,n_toroidal
           ix = ir-1-nx0/2
           if (ix < 0) ix = ix+nx

           iy = in-1
           g_nl(ir,j,in) = fxmany(iy,ix,j)
        enddo
     enddo
  enddo

  ! end data g_nl
!$acc end data

  ! end data f_nl
!$acc end data

  call timer_lib_out('nl_mem')

  call timer_lib_in('nl_comm')
  call parallel_slib_r_nc_gpu(g_nl,gpack)
  call timer_lib_out('nl_comm')

  call timer_lib_in('nl_mem')

!$acc parallel loop collapse(3) independent private(iexch,ic_loc) &
!$acc&         present(psi,gpack) default(none)
  do iv_loc=1,nv_loc
     do it=1,n_theta
        do ir=1,n_radial
           iexch = it + (iv_loc-1)*n_theta
           ic_loc = it + (ir-1)*n_theta
           psi(ic_loc,iv_loc) = gpack(ir,iexch) 
        enddo
     enddo
  enddo

  call timer_lib_out('nl_mem')

  ! Filter
  if (n == 0) then
!$acc parallel loop gang vector private(ir) &
!$acc          present(ir_c,psi,px) default(none)
      do ic=1,nc
        ir = ir_c(ic)
        if (ir == 1 .or. px(ir) == 0) then
           psi(ic,:) = 0.0
        endif
     enddo
  endif

  ! RHS -> -[f,g] = [f,g]_{r,-alpha}

  call timer_lib_in('nl')

!$acc parallel loop collapse(2) independent present(rhs(:,:,ij),psi)
  do iv_loc=1,nv_loc
     do ic_loc=1,nc
        rhs(ic_loc,iv_loc,ij) = rhs(ic_loc,iv_loc,ij)+((q*rho/rmin)*(2*pi/length))*psi(ic_loc,iv_loc)
     enddo
  enddo

  call timer_lib_out('nl')

end subroutine cgyro_nl_fftw
