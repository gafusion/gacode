module vis

  ! torcut
  ! torside
  ! realfluct

contains  

  subroutine torcut(dn,m,q,nr,nth,nn,nx,nz,thi,g1,g2,c,f)

    implicit none

    integer, intent(in) :: dn,m
    double precision, intent(in) :: q
    integer, intent(in) :: nr,nth,nn,nx,nz
    double precision, intent(in) :: thi(0:nth-1)
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
    ! f2py intent(in) thi
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

    th(0:nth-1) = thi
    th(nth) = -th(0)

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
          ! Interpolation of ballooning potential (on coarse, irregular th-grid)
          ! onto fine, equally-spaced z-grid
          if (z(k) >= th(kc) .and. z(k) <= th(kc+1)) then 
             c0(:,:) = ( cx(:,kc+1,:)*(z(k)-th(kc)) &
                  + cx(:,kc,:)*(th(kc+1)-z(k)) )/(th(kc+1)-th(kc))
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

  subroutine torside(dn,m,q,phimax,nr,nth,nn,nphi,nz,thi,g1,g2,c,f)

    implicit none

    integer, intent(in) :: dn,m
    double precision, intent(in) :: q
    double precision, intent(in) :: phimax
    integer, intent(in) :: nr,nth,nn,nphi,nz
    double precision, intent(in) :: thi(0:nth-1)
    double precision, intent(in) :: g1(0:nz-1),g2(0:nz-1)
    double complex, intent(in) :: c(0:nr-1,0:nth-1,0:nn-1)
    double precision, intent(inout) :: f(0:1,0:nz-1,0:nphi-1)

    double complex, dimension(:,:), allocatable :: eny,c0
    double complex, dimension(:,:,:), allocatable :: cx
    double precision, dimension(:), allocatable :: phi,z,th

    integer :: i,k,kc
    integer :: n,p,pp
    double precision :: pi,fsum1,fsum2,x1,x2,xb1,xb2
    double complex :: ic

    ! f2py intent(in) dn
    ! f2py intent(in) m
    ! f2py intent(in) q
    ! f2py intent(in) phimax
    ! f2py intent(in) nr
    ! f2py intent(in) nth
    ! f2py intent(in) nn
    ! f2py intent(in) nphi
    ! f2py intent(in) nz
    ! f2py intent(in) thi
    ! f2py intent(in) g1
    ! f2py intent(in) g2
    ! f2py intent(in) c
    ! f2py intent(in,out) f

    allocate(eny(0:nn-1,0:nz-1))
    allocate(phi(0:nphi-1))
    allocate(z(0:nz-1))
    allocate(th(0:nth))
    allocate(c0(0:nr-1,0:nn-1))
    allocate(cx(0:nr-1,0:nth,0:nn-1))

    ic = (0d0,1d0)
    pi = atan(1d0)*4d0

    th(0:nth-1) = thi
    th(nth) = -th(0)

    xb1 = 0.0
    xb2 = 2*pi/dn
    
    do i=0,nphi-1
       phi(i) = -i*2*pi/(nphi-1)*phimax/dn
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

!$omp parallel do private(k,kc,c0,i,fsum1,fsum2,x1,x2,n,p)
    do k=0,nz-1

       do kc=0,nth-1
          ! Interpolation of ballooning potential (on coarse, irregular th-grid)
          ! onto fine, equally-spaced z-grid
          if (z(k) >= th(kc) .and. z(k) <= th(kc+1)) then 
             c0(:,:) = ( cx(:,kc+1,:)*(z(k)-th(kc)) &
                  + cx(:,kc,:)*(th(kc+1)-z(k)) )/(th(kc+1)-th(kc))
             exit
          endif
       enddo

       do i=0,nphi-1
          fsum1 = 0.0
          fsum2 = 0.0
          x1 = m*g2(k)*xb1/(2*pi)-phi(i) 
          x2 = m*g2(k)*xb2/(2*pi)-phi(i) 
          if (nn > 1) then
             do n=0,nn-1
                do p=0,nr-1
                   fsum1 = fsum1+real( c0(p,n)*eny(n,k)*exp(dn*ic*n*x1) )
                   fsum2 = fsum2+real( c0(p,n)*eny(n,k)*exp(dn*ic*n*x2) )
                enddo
             enddo
          else
             do p=0,nr-1
                fsum1 = fsum1+real( c0(p,0)*eny(0,k)*exp(dn*ic*x1) )
                fsum2 = fsum2+real( c0(p,0)*eny(0,k)*exp(dn*ic*x2) )
             enddo
          endif
          f(0,k,i) = fsum1
          f(1,k,i) = fsum2
       enddo

    enddo

    deallocate(eny,phi,z,th,c0,cx)

  end subroutine torside

  subroutine realfluct(nr,nn,nx,ny,c,f)

    implicit none

    integer, intent(in) :: nr,nn,nx,ny
    double complex, intent(in) :: c(0:nr-1,0:nn-1)
    double precision, intent(inout) :: f(0:nx-1,0:ny-1)

    double complex, dimension(:,:), allocatable :: epx,eny

    integer :: i,j,n,p
    double precision :: pi,xi,yj,fsum
    double complex :: ic

    ! f2py intent(in) nr
    ! f2py intent(in) nn
    ! f2py intent(in) nx
    ! f2py intent(in) ny
    ! f2py intent(in) c
    ! f2py intent(in,out) f

    allocate(epx(0:nr-1,0:nx-1))
    allocate(eny(0:nn-1,0:ny-1))

    ic = (0d0,1d0)
    pi = atan(1d0)*4d0

    do i=0,nx-1
       xi = i*2*pi/(nx-1)
       do p=0,nr-1    
          epx(p,i) = exp(ic*(p-nr/2)*xi)
       enddo
    enddo

    if (nn > 1) then
       do j=0,ny-1
          yj = j*2*pi/(ny-1)
          do n=0,nn-1    
             eny(n,j) = exp(-ic*n*yj)
          enddo
       enddo
       ! factor of 1/2 for n=0
       eny(0,:) = eny(0,:)/2
    else
       do j=0,ny-1
          yj = j*2*pi/(ny-1)
          n=1    
          eny(0,j) = exp(-ic*n*yj)
       enddo
    endif

!$omp parallel do private(j,i,fsum,n,p)
    do j=0,ny-1
       do i=0,nx-1
          fsum = 0d0
          if (nn > 1) then
             do n=0,nn-1
                do p=0,nr-1
                   fsum = fsum+real(c(p,n)*epx(p,i)*eny(n,j))
                enddo
             enddo
          else
             n=1
             do p=0,nr-1
                fsum = fsum+real(c(p,0)*epx(p,i)*eny(0,j))
             enddo
          endif
          f(i,j) = fsum
       enddo
    enddo

    deallocate(epx,eny)

  end subroutine realfluct

 
end module vis
