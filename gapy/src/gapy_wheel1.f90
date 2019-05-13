subroutine wheel1(nr,nth,nn,nx,nz,c,f)

  implicit none

  integer, intent(in) :: nr,nth,nn,nx,nz
  double complex, intent(in) :: c(0:nr-1,0:nth-1,0:nn-1)
  double precision, intent(inout) :: f(0:nx-1,0:nz-1)

  double complex, dimension(:,:), allocatable :: epx,c0
  double precision, dimension(:), allocatable :: eny
  double precision, dimension(:), allocatable :: z,th

  integer :: i,j,k,kc
  integer :: n,p
  double precision :: pi,fsum,xi
  double complex :: ic

  ! f2py intent(in) nr
  ! f2py intent(in) nth  
  ! f2py intent(in) nn
  ! f2py intent(in) nx
  ! f2py intent(in) nz
  ! f2py intent(in) c
  ! f2py intent(in,out) f

  allocate(epx(0:nr-1,0:nx-1))
  allocate(eny(0:nn-1))
  allocate(z(0:nz-1))
  allocate(th(0:nth))
  allocate(c0(0:nr-1,0:nn-1))

  ic = (0d0,1d0)
  pi = atan(1d0)*4d0

  do i=0,nx-1
     xi = i*2*pi/(nx-1)
     do p=0,nr-1    
        epx(p,i) = exp(ic*(p-nr/2)*xi)
     enddo
  enddo

  do kc=0,nth
     th(kc) = kc*2*pi/nth-pi
  enddo

  do k=0,nz-1
     z(k) = k*pi/(nz-1)-pi
  enddo

  if (nn > 1) then
     ! factor of 1/2 for n=0
     eny(:) = 1d0
     eny(0) = 0.5d0
  else
     eny(0) = 1d0
  endif

!$omp parallel do private(k,kc,c0,i,n,p,fsum)
  do k=0,nz-1
     do kc=0,nth-1
        if (z(k) >= th(kc) .and. z(k) <= th(kc+1)) then 
           c0(:,:) = ( c(:,kc+1,:)*(z(k)-th(kc)) &
                + c(:,kc,:)*(th(kc+1)-z(k)) )/(th(1)-th(0))
           exit
        endif
     enddo
     do i=0,nx-1
        fsum = 0d0
        if (nn > 1) then
           do n=0,nn-1
              do p=0,nr-1
                 fsum = fsum+real(c0(p,n)*epx(p,i)*eny(n))
              enddo
           enddo
        else
           n=1
           do p=0,nr-1
              fsum = fsum+real(c0(p,0)*epx(p,i)*eny(0))
           enddo
        endif
        f(i,k) = fsum
     enddo
  enddo

  deallocate(epx,eny,z,th,c0)

end subroutine wheel1
