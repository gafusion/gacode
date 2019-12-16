subroutine wheel2(nr,nth,nn,ny,nz,c,f)

  implicit none

  integer, intent(in) :: nr,nth,nn,ny,nz
  double complex, intent(in) :: c(0:nr-1,0:nth-1,0:nn-1)
  double precision, intent(inout) :: f(0:ny-1,0:nz-1)

  double complex, dimension(:,:), allocatable :: eny,c0
  double precision, dimension(:), allocatable :: y,z,th

  integer :: i,j,k,kc
  integer :: n,p
  double precision :: pi,fsum
  double complex :: ic

  ! f2py intent(in) nr
  ! f2py intent(in) nth  
  ! f2py intent(in) nn
  ! f2py intent(in) ny
  ! f2py intent(in) nz
  ! f2py intent(in) c
  ! f2py intent(in,out) f

  allocate(eny(0:nn-1,0:ny-1))
  allocate(c0(0:nr-1,0:nn-1))
  allocate(th(0:nth))
  allocate(y(0:ny-1))
  allocate(z(0:nz-1))

  ic = (0d0,1d0)
  pi = atan(1d0)*4d0

  do kc=0,nth
     th(kc) = kc*2*pi/nth-pi
  enddo

  if (nn > 1) then
     do j=0,ny-1
        y(j) = j*2*pi/(ny-1)
        do n=0,nn-1    
           eny(n,j) = exp(-ic*n*y(j))
        enddo
     enddo
     ! factor of 1/2 for n=0
     eny(0,:) = eny(0,:)/2
  else
     do j=0,ny-1
        y(j) = j*2*pi/(ny-1)
        n=1    
        eny(0,j) = exp(-ic*n*y(j))
     enddo
  endif

  do k=0,nz-1
     z(k) = k*pi/(nz-1)-pi
  enddo

!$omp parallel do private(k,kc,c0,j,n,p,fsum)
  do k=0,nz-1
     do kc=0,nth-1
        if (z(k) >= th(kc) .and. z(k) <= th(kc+1)) then 
           c0(:,:) = ( c(:,kc+1,:)*(z(k)-th(kc)) &
                + c(:,kc,:)*(th(kc+1)-z(k)) )/(th(1)-th(0))
           exit
        endif
     enddo
     do j=0,ny-1
        fsum = 0d0
        if (nn > 1) then
           do n=0,nn-1
              do p=0,nr-1
                 fsum = fsum+real(c0(p,n)*eny(n,j))
              enddo
           enddo
        else
           n=1
           do p=0,nr-1
              fsum = fsum+real(c0(p,0)*eny(0,j))
           enddo
        endif
        f(j,k) = fsum
     enddo
  enddo

  deallocate(eny,c0,th,y,z)

end subroutine wheel2
