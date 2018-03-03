!--------------------------------------------------
! Generate Legendre pseudospectral data on x=(-1,1) 
!
! x(n)    -> nodes 
! w(n)    -> weights with Sum(w)=2.
! d1(n,n) -> d/dx
! dl(n,n) -> d/dx (1-x^2) d/dx
!
!       1
!       /
! Note: | dx f(x) = Sum w(i) f(i)
!       /
!      -1
!---------------------------------------------------

subroutine pseudo_legendre(n,x,w,d1,dl)

  implicit none

  integer, intent(in) :: n

  real, intent(out) :: x(n)
  real, intent(out) :: w(n)
  real, intent(out) :: d1(n,n)
  real, intent(out) :: dl(n,n)

  integer :: i,j
  integer :: lwork,info
  integer, dimension(:), allocatable :: ipiv
  real, dimension(:), allocatable :: work
  real, dimension(:,:), allocatable :: c
  real, dimension(:,:), allocatable :: cp
  real, dimension(:,:), allocatable :: cl

  integer, parameter :: print_flag=0

  call gauss_legendre(-1.0,1.0,x,w,n)

  lwork = 2*n
  allocate(ipiv(n))
  allocate(work(lwork))

  allocate(c(n,n))
  allocate(cp(n,n))
  allocate(cl(n,n))

  do i=1,n
     do j=1,n
        ! d/dxi
        call pseudo_rec_legendre(j-1,x(i),c(i,j),cp(i,j))
        ! L
        cl(i,j) = -(j-1)*j*c(i,j)
     enddo
  enddo

  call DGETRF(n,n,c,n,ipiv,info)
  call DGETRI(n,c,n,ipiv,work,lwork,info)
  ! c -> C^(-1)
  call DGEMM('N','N',n,n,n,1.0,cp,n,c,n,0.0,d1,n)
  ! d1 -> (Cp C^(-1))
  call DGEMM('N','N',n,n,n,1.0,cl,n,c,n,0.0,dl,n)
  ! dl -> (L C^(-1))

  if (print_flag == 1) then
     print *
     print *,'Integration weights (x,w,e,w_e)'
     print *
     do i=1,n
        print '(i2,2x,10(1pe14.7,1x))',i,x(i),w(i)
     enddo
     print *
     print *,'Sum(w) = ',sum(w)
     print *
     print *,'Pseudospectral 1st derivative error on exp(-x)' 
     print *
     do i=1,n
        print '(i2,2x,10(1pe14.7,1x))',i,sum(d1(i,:)*exp(-x(:)))/(-exp(-x(i)))-1
     enddo
     print *
     print *,'Pseudospectral L-derivative error on exp(-x)' 
     print *
     do i=1,n
        print '(i2,2x,10(1pe14.7,1x))',i,sum(dl(i,:)*exp(-x(:)))/(exp(-x(i))*(1.0+2.0*x(i)-x(i)**2))-1
     enddo

  endif

end subroutine pseudo_legendre

!---------------------------------------------------
! Generate Maxwell pseudospectral data on e=(0,emax) 
!
! e(n)    -> nodes 
! w(n)    -> weights with Sum(w)=1.0
! d1(n,n) -> d/dx
!
! where e = emax x^2.
!
!               emax
!           2    /
! Note: -------- | de sqrt(e) exp(-e) h(x) = Sum w(i) h(i)
!       sqrt(pi) /
!                0
!----------------------------------------------------

subroutine pseudo_maxwell_new(n,emax,e,w,d1,datafile)

  use math_constants

  implicit none

  character (len=*), intent(in) :: datafile

  integer, intent(in) :: n
  real, intent(in) :: emax
  integer :: i

  real, intent(out) :: e(n)
  real, intent(out) :: w(n)
  real, intent(out) :: d1(n,n)

  open(unit=1,file=trim(datafile),status='old')
  do i=1,n
     read(1,*) e(i),w(i)
  enddo
  do i=1,n
     read(1,*) d1(i,:)
  enddo
  close(1)

end subroutine pseudo_maxwell_new

subroutine pseudo_orthog(n,nu,alpha,beta,a,b)  

  implicit none

  integer, intent(in) :: n
  real, intent(in) :: alpha(2*n-1),beta(2*n-1),nu(2*n) 
  real, intent(inout) :: a(n),b(n)  
  integer :: k,l  
  real, dimension(:,:), allocatable :: sig

  allocate(sig(2*n+1,2*n+1))  

  do l=3,2*n  
     sig(1,l) = 0.0  
  enddo
  do l=2,2*n+1  
     sig(2,l) = nu(l-1)  
  enddo
  a(1) = alpha(1)+nu(2)/nu(1)  
  b(1) = 0.0  
  do k=3,n+1  
     do l=k,2*n-k+3  
        sig(k,l) = sig(k-1,l+1)+(alpha(l-1)-a(k-2))*sig(k-1,l)  & 
             -b(k-2)*sig(k-2,l)+beta(l-1)*sig(k-1,l-1)  
     enddo
     a(k-1) = alpha(k-1)+sig(k,k+1)/sig(k,k)-sig(k-1,k)/sig(k-1,k-1)  
     b(k-1) = sig(k,k)/sig(k-1,k-1)  
  enddo

  deallocate(sig)

end subroutine pseudo_orthog

!-------------------------------------------------------
! p -> P(n,x)
! q -> P(n,x)'
!-------------------------------------------------------

subroutine pseudo_rec_legendre(n,x,p,q)

  implicit none

  integer, intent (in) :: n
  real, intent (in) :: x
  real, intent(out) :: p,q
  integer :: j
  real :: a,b
  real :: pmm,pm
  real :: qmm,qm

  if (n == 0) then
     p  = 1.0
     q  = 0.0
  else
     pmm = 0.0
     pm  = 1.0

     qmm = 0.0
     qm  = 0.0

     do j=0,n-1

        a = (2*j+1.0)/(j+1.0)
        b = j/(j+1.0)

        ! P(n,x)
        p = a*x*pm-b*pmm
        pmm = pm
        pm = p

        ! P'(n,x)
        q = a*x*qm-b*qmm + a*pmm
        qmm = qm
        qm  = q
     enddo
  endif

end subroutine pseudo_rec_legendre

real function p_legendre(n,x)

  implicit none

  integer, intent(in) :: n
  real, intent(in) :: x
  real :: p,q

  call pseudo_rec_legendre(n,x,p,q)

  p_legendre = p

end function p_legendre

!-------------------------------------------------------
! p -> Poly(n,x)
! q -> Poly(n,x)'
! r -> Poly(n,x)''
!-------------------------------------------------------
subroutine pseudo_rec_maxwell(n,x,p,q,r,am,bm,nm)

  implicit none

  integer, intent (in) :: n,nm
  real, intent (in) :: x
  real, intent (in) :: am(0:nm-1),bm(0:nm-1)
  real, intent(out) :: p,q,r
  integer :: j
  real :: pmm,pm
  real :: qmm,qm
  real :: rmm,rm

  if (n == 0) then
     p  = 1.0
     q  = 0.0
     r  = 0.0
  else
     pmm = 0.0
     pm  = 1.0
 
     qmm = 0.0
     qm  = 0.0

     rmm = 0.0
     rm  = 0.0

     do j=0,n-1

        ! P(n,x)
        p = (x-am(j))*pm-bm(j)*pmm
        pmm = pm
        pm = p

        ! P'(n,x)
        q = (x-am(j))*qm-bm(j)*qmm + pmm
        qmm = qm
        qm  = q

        ! P''(n,x)
        r = (x-am(j))*rm-bm(j)*rmm + 2*qmm
        rmm = rm
        rm  = r
     enddo
  endif

end subroutine pseudo_rec_maxwell
