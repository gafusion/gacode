subroutine torcut(dn,m,q,nr,nth,nn,nx,nz,g1,g2,c,f)

  implicit none

  integer, intent(in) :: dn,m
  double precision, intent(in) :: q
  integer, intent(in) :: nr,nn,nth,nx,nz
  double precision, intent(in) :: g1(0:nz-1),g2(0:nz-1)
  double complex, intent(in) :: c(0:nr-1,0:nth-1,0:nn-1)
  double precision, intent(inout) :: f(0:nx-1,0:nz-1)

  double complex, dimension(:,:), allocatable :: epx,eny,c0
  double complex, dimension(:,:,:), allocatable :: cx
  double precision, dimension(:), allocatable :: x,z,th

  integer :: i,k,kc
  integer :: n,p,pp
  double precision :: pi,fsum,x0
  double complex :: ic
  
  ! f2py intent(in) dn
  ! f2py intent(in) m
  ! f2py intent(in) q
  ! f2py intent(in) nr
  ! f2py intent(in) nth
  ! f2py intent(in) nn
  ! f2py intent(in) nx
  ! f2py intent(in) nz
  ! f2py intent(in) g1
  ! f2py intent(in) g2
  ! f2py intent(in) c
  ! f2py intent(in,out) f

  allocate(epx(0:nr-1,0:nx-1))
  allocate(eny(0:nn-1,0:nz-1))
  allocate(x(0:nx-1))
  allocate(z(0:nz-1))
  allocate(th(0:nth))
  allocate(c0(0:nr-1,0:nn-1))
  allocate(cx(0:nr-1,0:nth,0:nn-1))
  
  ic = (0d0,1d0)
  pi = atan(1d0)*4d0

  do kc=0,nth
     th(kc) = kc*2*pi/nth-pi
  enddo

  do i=0,nx-1
     x(i) = i*2*pi/(nx-1)/dn
     do p=0,nr-1    
        epx(p,i) = exp(dn*ic*(p-nr/2)*x(i))
     enddo
  enddo

  if (nn > 1) then
     do k=0,nz-1
        z(k) = k*2*pi/(nz-1)-pi
        do n=0,nn-1    
           eny(n,k) = exp(dn*ic*n*g1(k))
        enddo
     enddo
     ! factor of 1/2 for n=0
     eny(0,:) = eny(0,:)/2
  else
     do k=0,nz-1
        z(k) = k*2*pi/(nz-1)-pi
        eny(0,k) = exp(dn*ic*g1(k))
     enddo
  endif

  ! Extended coefficients
  cx(:,0:nth-1,:) = c(:,:,:) 
  if (nn > 1) then 
     do n=0,nn-1
        do p=0,nr-1
           pp = p+n*m
           if (pp > nr-1) pp = pp-nr
           cx(p,nth,n) = c(pp,0,n)*exp(-2*pi*dn*ic*n*q)
        enddo
     enddo
  else
     do p=0,nr-1
        pp = p+1*m
        if (pp > nr-1) pp = pp-nr
        cx(p,nth,0) = c(pp,0,0)*exp(-2*pi*dn*ic*q)
     enddo
  endif

!$omp parallel do private(k,kc,c0,i,fsum,x0,n,p)
  do k=0,nz-1

     do kc=0,nth-1
        if (z(k) >= th(kc) .and. z(k) <= th(kc+1)) then 
           c0(:,:) = ( cx(:,kc+1,:)*(z(k)-th(kc)) &
                + cx(:,kc,:)*(th(kc+1)-z(k)) )/(th(1)-th(0))
           exit
        endif
     enddo

     do i=0,nx-1
        fsum = 0.0
        x0 = m*x(i)*g2(k)/(2*pi)
        if (nn > 1) then
           do n=0,nn-1
              do p=0,nr-1
                 fsum = fsum+real( c0(p,n)*epx(p,i)*eny(n,k)*exp(dn*ic*n*x0) )
              enddo
           enddo
        else
           do p=0,nr-1
              fsum = fsum+real( c0(p,0)*epx(p,i)*eny(0,k)*exp(dn*ic*x0) )
           enddo
        endif
        f(i,k) = fsum
     enddo

  enddo

  deallocate(epx,eny,x,z,th,c0,cx)

end subroutine torcut
