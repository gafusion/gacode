!-----------------------------------------------------------------
! cgyro_nl_direct.f90
!
! PURPOSE:
!  Evaluate nonlinear bracket with direct dealiased convolution.
!  Need to be careful with (p=-nr/2,n=0) component.
!-----------------------------------------------------------------

subroutine cgyro_nl_direct(ij)

  use timer_lib
  use parallel_lib

  use cgyro_globals

  implicit none

  integer, intent(in) :: ij
  integer :: ix,ixp,iy,iyp,iexch
  integer :: ixpp,iypp
  integer :: ir,it,j,in
  integer, dimension(:), allocatable :: pcyc
  integer, dimension(:), allocatable :: ncyc
  complex, dimension(:,:), allocatable :: f
  complex, dimension(:,:), allocatable :: g
  complex, dimension(:,:), allocatable :: fg


  ! To keep code self-contained, just recompute these every time
  ! rather than storing in globals.  Cost is tiny.
  allocate(pcyc(-3*nx:3*nx-1))
  allocate(ncyc(-3*ny-1:3*ny+1))
  do ix=-nx,nx-1
     pcyc(ix-2*nx) = ix
     pcyc(ix) = ix
     pcyc(ix+2*nx) = ix
  enddo
  do iy=-ny,ny
     ncyc(iy-(2*ny+1)) = iy
     ncyc(iy) = iy
     ncyc(iy+(2*ny+1)) = iy
  enddo

  allocate( f(-nx:nx,-ny:ny) )
  allocate( g(-nx:nx,-ny:ny) )
  allocate(fg(-nx:nx,-ny:ny) )

  call timer_lib_in('nl_comm')
  if ((nv_loc*n_theta) /= nsplit*n_toroidal) then
!$omp parallel do private(iexch)
     do iexch=1,nv_loc*n_theta
        fpack(:,iexch) = h_x(:,iexch)
        gpack(:,iexch) = psi(:,iexch)
     enddo

     do iexch=nv_loc*n_theta+1,nsplit*n_toroidal
        fpack(:,iexch) = (0.0,0.0)
        gpack(:,iexch) = (0.0,0.0)
     enddo

     call parallel_slib_f_nc(fpack,f_nl)
     call parallel_slib_f_nc(gpack,g_nl)
  else
     call parallel_slib_f_nc(h_x,f_nl)
     call parallel_slib_f_nc(psi,g_nl)
  endif

  call timer_lib_out('nl_comm')

  call timer_lib_in('nl')
  do j=1,nsplit
     do it=1,n_theta 

        f = 0.0
        g = 0.0

        ! Array mapping
        do in=1,n_toroidal
           iy = in-1
           do ir=1,n_radial
              ic = ic_c(ir,it) 
              ix = ir-1-nx0
              f(ix,iy) = f_nl(ic,j,in)
              g(ix,iy) = g_nl(ic,j,in)
              f(-ix,-iy) = conjg(f(ix,iy))
              g(-ix,-iy) = conjg(g(ix,iy)) 
           enddo
        enddo

        ! Reality
        do ix=-nx0,nx0
           do iy=1,ny0
              f(-ix,-iy) = conjg(f(ix,iy))
              g(-ix,-iy) = conjg(g(ix,iy))
           enddo
        enddo

        fg = (0.0,0.0)
        do ix=-nx0,nx0-1
           do ixp=-nx,nx-1
              ixpp = pcyc(ix-ixp)
              do iy=0,ny0
                 do iyp=-ny+iy,ny
                    iypp = ncyc(iy-iyp)
                    fg(ix,iy) = fg(ix,iy)+f(ixpp,iypp)*g(ixp,iyp)*(iypp*ixp-iyp*ixpp)
                 enddo
              enddo
           enddo
        enddo
        ! Must annhilate n=0,p=-nr/2
        fg(-nx0,0) = 0.0

        do ir=1,n_radial
           ic = ic_c(ir,it) 
           ix = ir-1-nx0
           do in=1,n_toroidal
              iy = in-1
              g_nl(ic,j,in) = fg(ix,iy)
           enddo
        enddo

     enddo ! it
  enddo ! j
  call timer_lib_out('nl')

  call timer_lib_in('nl_comm')
  if ((nv_loc*n_theta) /= nsplit*n_toroidal) then
    call parallel_slib_r_nc(g_nl,gpack)
!$omp parallel do private(iexch)
     do iexch=1,nv_loc*n_theta
        psi(:,iexch) = gpack(:,iexch)
     enddo
  else
    call parallel_slib_r_nc(g_nl,psi)
  endif
  call timer_lib_out('nl_comm')

  deallocate( f)
  deallocate( g)
  deallocate(fg)
  
  deallocate(pcyc)
  deallocate(ncyc)

  ! RHS -> -[f,g] = (n'' p' - n' p'') f'' g' 

  rhs(:,:,ij) = rhs(:,:,ij)+(q*rho/rmin)*(2*pi/length)*psi(:,:)

end subroutine cgyro_nl_direct
