!-----------------------------------------------------------------
! cgyro_nl_fftw.f90
!
! PURPOSE:
!  Evaluate nonlinear bracket with dealiased FFT.  It is natural 
!  to use the FFTW complex-to-real (c2r) transform to move to real 
!  space, compute the convolution (uv), and transform back to the 
!  spectral form using a real-to-complex transform (r2c).
!  
! NOTE: Need to be careful with (p=-nr/2,n=0) component.
!-----------------------------------------------------------------

! NOTE: call cgyro_nl_fftw_comm1/2_async before cgyro_nl_fftw
subroutine cgyro_nl_fftw_comm1_async

  use timer_lib
  use parallel_lib
  use cgyro_globals

  implicit none

  integer :: ir,it,iv_loc_m
  integer :: iexch

  call timer_lib_in('nl_mem')

!$omp parallel do private(iv_loc,it,iexch,ir)
  do iv_loc_m=1,nv_loc
     iexch = (iv_loc_m-1)*n_theta
     do it=1,n_theta
        iexch = iexch+1
        do ir=1,n_radial
           fpack(ir,iexch) = h_x(ic_c(ir,it),iv_loc_m)
        enddo
     enddo
  enddo

  do iexch=nv_loc*n_theta+1,nsplit*n_toroidal
     fpack(:,iexch) = (0.0,0.0)
  enddo

  call parallel_slib_f_nc_async(fpack,f_nl,f_req)

  call timer_lib_out('nl_mem')

end subroutine cgyro_nl_fftw_comm1_async

! Note: Calling test propagates the async operations in some MPI implementations
subroutine cgyro_nl_fftw_comm1_test
  use parallel_lib
  use cgyro_globals

  implicit none

  call parallel_slib_test(f_req)

end subroutine cgyro_nl_fftw_comm1_test

subroutine cgyro_nl_fftw_comm2_async

  use timer_lib
  use parallel_lib
  use cgyro_globals

  implicit none

  integer :: ir,it,iv_loc_m
  integer :: iexch

  call timer_lib_in('nl_mem')

  if (nonlinear_field .ne. 0) then
!$omp parallel do private(iv_loc,it,iexch,ir)
     do iv_loc_m=1,nv_loc
        iexch = (iv_loc_m-1)*n_theta
        do it=1,n_theta
           iexch = iexch+1
           do ir=1,n_radial
              gpack(ir,iexch) = jvec_c(nonlinear_field,ic_c(ir,it),iv_loc_m)*field(nonlinear_field,ic_c(ir,it))
           enddo
        enddo
     enddo
  else
!$omp parallel do private(iv_loc,it,iexch,ir)
     do iv_loc_m=1,nv_loc
        iexch = (iv_loc_m-1)*n_theta
        do it=1,n_theta
           iexch = iexch+1
           do ir=1,n_radial
              gpack(ir,iexch) = psi(ic_c(ir,it),iv_loc_m)
           enddo
        enddo
     enddo
  end if

  do iexch=nv_loc*n_theta+1,nsplit*n_toroidal
     gpack(:,iexch) = (0.0,0.0)
  enddo

  call parallel_slib_f_nc_async(gpack,g_nl,g_req)

  call timer_lib_out('nl_mem')

end subroutine cgyro_nl_fftw_comm2_async

! Note: Calling test propagates the async operations in some MPI implementations
subroutine cgyro_nl_fftw_comm2_test
  use parallel_lib
  use cgyro_globals

  implicit none

  call parallel_slib_test(g_req)

end subroutine cgyro_nl_fftw_comm2_test

subroutine cgyro_nl_fftw_stepr(j, i_omp)

  use timer_lib
  use parallel_lib
  use cgyro_globals

  implicit none

  integer, intent(in) :: j, i_omp
  integer :: ix,iy
  integer :: ir,in

  include 'fftw3.f03'

  ! Poisson bracket in real space

  uv(:,:,i_omp) = (uxmany(:,:,j)*vy(:,:,i_omp)-uymany(:,:,j)*vx(:,:,i_omp))/(nx*ny)

  call fftw_execute_dft_r2c(plan_r2c,uv(:,:,i_omp),fx(:,:,i_omp))

  ! NOTE: The FFT will generate an unwanted n=0,p=-nr/2 component
  ! that will be filtered in the main time-stepping loop

  ! this should really be accounted against nl_mem, but hard to do with OMP
  do ir=1,n_radial
     ix = ir-1-nx0/2
     if (ix < 0) ix = ix+nx
     do in=1,n_toroidal
        iy = in-1
        g_nl(ir,j,in) = fx(iy,ix,i_omp)
     enddo
  enddo

end subroutine cgyro_nl_fftw_stepr

! NOTE: call cgyro_nl_fftw_comm1 before cgyro_nl_fftw
subroutine cgyro_nl_fftw(ij)

  use timer_lib
  use parallel_lib
  use cgyro_globals

  implicit none

  integer, intent(in) :: ij
  integer :: ix,iy
  integer :: ir,it,in
  integer :: j,p,iexch
  integer :: i_omp
  logical :: force_early_comm2, one_pass_fft
  integer :: o,num_one_pass

  complex :: f0,g0

  integer, external :: omp_get_thread_num

  include 'fftw3.f03'
  
  ! no staggering, will need all threads for FFTW in one pass
  force_early_comm2 = (n_omp>=(4*nsplit)) 
  one_pass_fft = force_early_comm2

  ! time to wait for the F_nl to become avaialble
  call timer_lib_in('nl_comm')
  call parallel_slib_f_nc_wait(fpack,f_nl,f_req)
  ! make sure g_req progresses
  call parallel_slib_test(g_req)
  call timer_lib_out('nl_comm')

  if (force_early_comm2) then
     ! time to wait for the g_nl to become avaialble
     call timer_lib_in('nl_comm')
     call parallel_slib_f_nc_wait(gpack,g_nl,g_req)
     call timer_lib_out('nl_comm')
  endif

  call timer_lib_in('nl')

  if (n_omp<=nsplit) then
!$omp parallel private(in,iy,ir,p,ix,f0,i_omp,j)
!$omp do schedule(dynamic,1)
     do j=1,nsplit
        i_omp = omp_get_thread_num()+1

        fx(:,:,i_omp) = 0.0
        fy(:,:,i_omp) = 0.0

        ! Array mapping
        do ir=1,n_radial
           p  = ir-1-nx0/2
           ix = p
           if (ix < 0) ix = ix+nx
           do in=1,n_toroidal
              iy = in-1
              f0 = i_c*f_nl(ir,j,in)
              fx(iy,ix,i_omp) = p*f0
              fy(iy,ix,i_omp) = iy*f0
           enddo
        enddo

        call cleanx(fx(0,:,i_omp),nx)
        call fftw_execute_dft_c2r(plan_c2r,fx(:,:,i_omp),uxmany(:,:,j))
        call cleanx(fy(0,:,i_omp),nx)
        call fftw_execute_dft_c2r(plan_c2r,fy(:,:,i_omp),uymany(:,:,j))
     enddo ! j
!$omp end do
!$omp end parallel
  else ! (n_omp>nsplit), increase parallelism
     num_one_pass = 2
     if (one_pass_fft) then
        num_one_pass = 4
     endif
!$omp parallel private(in,iy,ir,p,ix,f0,i_omp,j,o)
!$omp do schedule(dynamic,1) collapse(2)
     do j=1,nsplit
        do o=1,num_one_pass
           i_omp = j ! j<n_omp in this branch

           select case(o)
           case (1)
              fx(:,:,i_omp) = 0.0

              ! Array mapping
              do ir=1,n_radial
                 p  = ir-1-nx0/2
                 ix = p
                 if (ix < 0) ix = ix+nx
                 do in=1,n_toroidal
                    iy = in-1
                    f0 = i_c*f_nl(ir,j,in)
                    fx(iy,ix,i_omp) = p*f0
                 enddo
              enddo

              call cleanx(fx(0,:,i_omp),nx)
              call fftw_execute_dft_c2r(plan_c2r,fx(:,:,i_omp),uxmany(:,:,j))

           case (2)
              fy(:,:,i_omp) = 0.0

              ! Array mapping
              do ir=1,n_radial
                 p  = ir-1-nx0/2
                 ix = p
                 if (ix < 0) ix = ix+nx
                 do in=1,n_toroidal
                    iy = in-1
                    f0 = i_c*f_nl(ir,j,in)
                    fy(iy,ix,i_omp) = iy*f0
                 enddo
              enddo

              call cleanx(fy(0,:,i_omp),nx)
              call fftw_execute_dft_c2r(plan_c2r,fy(:,:,i_omp),uymany(:,:,j))

           case (3)
              gx(:,:,i_omp) = 0.0

              ! Array mapping
              do ir=1,n_radial
                 p  = ir-1-nx0/2
                 ix = p
                 if (ix < 0) ix = ix+nx
                 do in=1,n_toroidal
                    iy = in-1
                    g0 = i_c*g_nl(ir,j,in)
                    gx(iy,ix,i_omp) = p*g0
                 enddo
              enddo

              call cleanx(gx(0,:,i_omp),nx)
              call fftw_execute_dft_c2r(plan_c2r,gx(:,:,i_omp),vx(:,:,i_omp))

           case (4)
              gy(:,:,i_omp) = 0.0

              ! Array mapping
              do ir=1,n_radial
                 p  = ir-1-nx0/2
                 ix = p
                 if (ix < 0) ix = ix+nx
                 do in=1,n_toroidal
                    iy = in-1
                    g0 = i_c*g_nl(ir,j,in)
                    gy(iy,ix,i_omp) = iy*g0
                 enddo
              enddo

              call cleanx(gy(0,:,i_omp),nx)
              call fftw_execute_dft_c2r(plan_c2r,gy(:,:,i_omp),vy(:,:,i_omp))

           end select
        enddo ! o
     enddo ! j
!$omp end do
!$omp end parallel
  endif

  call timer_lib_out('nl')
  
  ! stagger comm2, to load ballance network traffic
  if (.not. force_early_comm2) then
     ! time to wait for the g_nl to become avaialble
     call timer_lib_in('nl_comm')
     call parallel_slib_f_nc_wait(gpack,g_nl,g_req)
     call timer_lib_out('nl_comm')
  endif

  call timer_lib_in('nl')

  if (n_omp <= nsplit) then
!$omp parallel private(in,iy,ir,p,ix,g0,i_omp,j)
!$omp do schedule(dynamic,1)
     do j=1,nsplit
        i_omp = omp_get_thread_num()+1

        gx(:,:,i_omp) = 0.0
        gy(:,:,i_omp) = 0.0

        ! Array mapping
        do ir=1,n_radial
           p  = ir-1-nx0/2
           ix = p
           if (ix < 0) ix = ix+nx  
           do in=1,n_toroidal
              iy = in-1
              g0 = i_c*g_nl(ir,j,in)
              gx(iy,ix,i_omp) = p*g0
              gy(iy,ix,i_omp) = iy*g0
           enddo
        enddo

        call cleanx(gx(0,:,i_omp),nx)
        call fftw_execute_dft_c2r(plan_c2r,gx(:,:,i_omp),vx(:,:,i_omp))
        call cleanx(gy(0,:,i_omp),nx)
        call fftw_execute_dft_c2r(plan_c2r,gy(:,:,i_omp),vy(:,:,i_omp))

        call cgyro_nl_fftw_stepr(j, i_omp)
     enddo ! j
!$omp end do
!$omp end parallel
  else ! n_omp>nsplit
     if (.not. one_pass_fft) then
!$omp parallel private(in,iy,ir,p,ix,g0,i_omp,j)
!$omp do schedule(dynamic,1) collapse(2)
        do j=1,nsplit
           do o=1,2
              i_omp = j ! j<n_omp in this branch, so we can do it

              if (o == 1) then
                 gx(:,:,i_omp) = 0.0

                 ! Array mapping
                 do ir=1,n_radial
                    p  = ir-1-nx0/2
                    ix = p
                    if (ix < 0) ix = ix+nx
                    do in=1,n_toroidal
                       iy = in-1
                       g0 = i_c*g_nl(ir,j,in)
                       gx(iy,ix,i_omp) = p*g0
                    enddo
                 enddo

                 call cleanx(gx(0,:,i_omp),nx)
                 call fftw_execute_dft_c2r(plan_c2r,gx(:,:,i_omp),vx(:,:,i_omp))
              else
                 gy(:,:,i_omp) = 0.0

                 ! Array mapping
                 do ir=1,n_radial
                    p  = ir-1-nx0/2
                    ix = p
                    if (ix < 0) ix = ix+nx
                    do in=1,n_toroidal
                       iy = in-1
                       g0 = i_c*g_nl(ir,j,in)
                       gy(iy,ix,i_omp) = iy*g0
                    enddo
                 enddo

                 call cleanx(gy(0,:,i_omp),nx)
                 call fftw_execute_dft_c2r(plan_c2r,gy(:,:,i_omp),vy(:,:,i_omp))
              endif
           enddo ! o
        enddo ! j
!$omp end do
!$omp end parallel
     endif

!$omp parallel private(j)
!$omp do schedule(dynamic,1)
     do j=1,nsplit
        call cgyro_nl_fftw_stepr(j, j) ! we used i_omp=j in the g section
     enddo ! j
!$omp end do
!$omp end parallel

  endif

  call timer_lib_out('nl')

  call timer_lib_in('nl_comm')
  call parallel_slib_r_nc(g_nl,gpack)
  call timer_lib_out('nl_comm')

  call timer_lib_in('nl')

!$omp parallel do private(iv_loc,it,iexch,ir,ic)
  do iv_loc=1,nv_loc
     iexch = (iv_loc-1)*n_theta
     do it=1,n_theta
        iexch = iexch+1
        do ir=1,n_radial
           psi(ic_c(ir,it),iv_loc) = gpack(ir,iexch)
        enddo
     enddo

     ! Filter
     if (n == 0) then
      do ic=1,nc
        ir = ir_c(ic)
        if (ir == 1 .or. px(ir) == 0) then
           psi(ic,iv_loc) = 0.0
        endif
      enddo
     endif

     ! RHS -> -[f,g] = [f,g]_{r,-alpha}

     rhs(:,iv_loc,ij) = rhs(:,iv_loc,ij)+(q*rho/rmin)*(2*pi/length)*psi(:,iv_loc)
  enddo

  call timer_lib_out('nl')

end subroutine cgyro_nl_fftw

subroutine cleanx(f,n)

  implicit none

  integer, intent(in) :: n
  complex, intent(inout) :: f(0:n-1)

  integer :: i

  ! Average elements so as to ensure
  !
  !   f(kx,ky=0) = f(-kx,ky=0)^*
  !
  ! This symmetry is required for complex input to the FFTW
  ! c2r transform 
  
  do i=1,n/2-1
     f(i)   = 0.5*( f(i)+conjg(f(n-i)) )
     f(n-i) = conjg(f(i)) 
  enddo

end subroutine cleanx
