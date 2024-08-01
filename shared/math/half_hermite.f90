!gfortran -march=native -O3 -fno-stack-arrays -fimplicit-none -fdefault-real-8 half_hermite.f90 -c
!gfortran -march=native -O3 -fcheck=all -ffpe-trap=invalid -fno-stack-arrays -fimplicit-none -fdefault-real-8 half_hermite.f90 -c
!intel:
!ifort -stand f15 -warn all -march=native -O3 -heap-arrays 10 -implicitnone -real-size 64 half_hermite.f90 -c
module half_hermite
  real, private, parameter :: pi1=atan(1.)*4
  integer, private, parameter :: nsafe=3,maxmem=50
  ! The code uses mainly nintervals intervals, and for these nintervals-1
  ! vertices for integration (the interior vertices). In addition, as a
  ! safety margin nsafe+1 vertices are added on each side and checked
  ! whether any such contribution is substantial, when the domain is extended
  ! by nsafe+1 intervals on the respective side.
  ! If the first interior point has index 1, the last interior vertex has
  ! index nintervals-1 and the first "safe" point has index -nsafe and the
  ! last of those is nintervals+nsafe. (That's why it's nsafe+1)
  private realloc
contains
  subroutine realloc(step,mw,mp,poly,poly1,minalloc,maxalloc,i0,nintervals,nn,vb)
    implicit none
    integer,intent(in) :: step
    integer,intent(inout) :: minalloc,maxalloc,i0,nintervals,nn
    integer :: newmin,newmax
    real, dimension(:), pointer,intent(inout) :: mw,mp,poly,poly1 !Mori weights and points
    real, dimension(:), pointer :: mw2,mp2,poly2,poly12 !Mori weights and
    !points
    integer istart,iend,isn,ien,i
    logical vb
    ! we assume nintervals has already been updated.
    ! if step=2 ==> the points are supposed to be spread with stepsize step.
    istart=i0-(2*nsafe)/step
    iend=i0+(2*(nintervals+nsafe))/step
    newmin=-10*nsafe-9
    newmax=10*nsafe+9+2*nintervals
    if (vb) print *,'Reallocating',minalloc,maxalloc,newmin,newmax,step
    if ((istart-i0)*step<newmin .or. (iend-i0)*step>newmax) then
       print *,'Serious error reallocating in half_hermite:'
       print *,'Reallocating',minalloc,i0,nintervals,maxalloc,newmin,newmax,step
       stop
    endif
    if (newmax-newmin>2048+maxmem*nn) then
       print *,'Too much memory needed in half_hermite',newmax-newmin&
            ,maxmem*nn,maxmem,nn,newmin,newmax
       stop
    endif
    allocate(mw2(newmin:newmax),mp2(newmin:newmax),poly2(newmin:newmax)&
         ,poly12(newmin:newmax))
    if (vb) then
       mw2=-500 ! for debugging purposes
       mp2=-500
       poly2=-500
       poly12=-500
    end if
    isn=(istart-i0)*step
    ien=(iend-i0)*step
    mw2(isn:ien:step)=mw(istart:iend)
    mp2(isn:ien:step)=mp(istart:iend)
    poly2(isn:ien:step)=poly(istart:iend)
    poly12(isn:ien:step)=poly1(istart:iend)

    do i=isn,ien,2*step
       if (mw2(i)==-500) then
          print *,'Very serious error i',i,isn,ien,istart,iend
          print *,mw(istart:min(istart+10,iend):2)
          stop
       endif
    end do

    deallocate(mw,mp,poly,poly1)
    mw=>mw2
    mp=>mp2
    poly=>poly2
    poly1=>poly12
    minalloc=newmin
    maxalloc=newmax
    i0=0
    if (vb) print *,'done'
  end subroutine realloc
    

  subroutine half_hermite_norm(n,xmin,xmax,sign,alpha,a1,b1,c1,a,bsq,logg,verbos&
       &e)
    implicit none
    integer, intent(in) :: n
    integer, intent(in),optional :: verbose
    real, intent(in) :: xmin,xmax,sign,alpha
    real, dimension(n), intent(out) :: a1,b1,c1,a,bsq,logg  !note: a is as
    real, dimension(:), allocatable :: gw,gp,gw2,gp2 !Gaussweights and points
    !below and bsq=sqrt(b) with b as below.
    real :: t,t0,t1,dt,dxdt,x,p0,p1,pi,s1,s2,s3,s1x,s2x,s3x,s4,s4x,s5,s5x
    real :: t00,t10
    real, parameter :: mori=.5,epsmori=1e-6,epsbase=5e-16,epsgauss=1e-18
    integer ngauss
    real w,sech2,sech,sechb,tanhp1,tanhm1,tanhb
    real eps
    integer i,j,k,ia
    integer leftnew,rightnew
    logical makenew
    logical vb

    real, dimension(:), pointer :: mw,mp,poly,poly1 !Mori weights and points
    !and polynomial recursion variables
    integer :: minalloc,maxalloc,i0,nintervals,nn ! current limits of allocation
    !i0 points to where the index 0 is.

    ! calculate recursion coefficients for orthogonal polynomials for weight
    ! function 

    !       w(x) = exp(-sign*x^2) x^alpha
    ! sign=+1 if xmax>0 and -1 otherwise.

    ! For the integration range x = 0 ... |xmax|

    ! i.e. <PQ> := integral_0...xmax w(x) P(x) Q(x) dx

    ! Each polynomial obeys the recursion relation
    ! p(i+1)=(x-a(i)) p(i) - b(i) p(i-1)
    ! p(1)=1
    ! Polynomial "i" has degree i-1 and consists of i monomials.
    ! Coefficients a(1..n) b(1..n) sufficient to generate polynomial p(n+1) are
    ! calculated. (actually sqrt(b(i)) is stored in bsq(i))
    ! In addition, the Norms <P(i)^2> are stored in g(1 ...n)
    ! 
    ! Further relations are:
    ! <P(i)^2>=<P(i)x^(i-1)>=<P(i) P(i-1) x>=b(i) <P(i-1)^2>
    ! b(i)=g(i)/g(i-1)
    ! <P(i+1) P(i)>=0 ==> <x P(i)^2>=a(i) <P(i)^2>
    ! a(i)=<x P(i)^2>/g(i)
    ! normalised polynomials pn(i)=p(i)/sqrt(g(i)) can be gained thus from the
    ! recursion
    ! a(i)=<x pn(i)^2>
    ! pn(i+1)*g(i+1)^.5=(x-a(i)) pn(i)g(i)^.5 - b(i) pn(i-1)g(i-1)^.5
    ! pn(i+1)*b(i+1)^.5=(x-a(i)) pn(i) - b(i)^.5 pn(i-1)
    ! pn(i+1)=(x-a(i))/b(i+1)^.5 pn(i) - (b(i)/b(i+1))^.5 pn(i-1)
    ! define a1(i+1)=a(i)/b(i+1)^.5 b1(i+1)=b(i)/b(i+1)^.5
    ! and c1(i+1)=1/b(i+1)^.5
    ! Then pn(i+1)=(c1(i+1)*x-a1(i+1))*pn(i)-b1(i+1) pn(i-1) (i>=1)
    ! and pn(0)=0 pn(-1)=0 (not part of the orthonormal polynomials
    ! and a1(1)=b1(1)=0 and c1(1)=1/g(1)^.5 pn(1)=c1(1) (not part of recursion)
    ! Recursive calculation:
    ! once we know pn(i) we also know a1(i) b1(i) c1(i).
    ! assume we know in addition a(i), b(i) from the previous step.
    ! b(i)^.5=1/c1(i) i>=2
    ! calculate pn(i+1)*b(i+1)^.5=(x-a(i))pn(i)-b(i)^.5 pn(i-1)
    ! can do this with preliminarily setting c1(i+1)=1 a1(i+1)=a(i)
    ! b1(i+1)=b(i)^.5
    ! a(i+1)b(i+1)= <x (pn(i+1)*b(i+1)^.5)^2 >
    ! b(i+1)= <(pn(i+1)*b(i+1)^.5)^2 >
    !
    ! Klaus Hallatschek
    !
    ! One could also normalise the polynomials with respect to their value at
    ! 0, or at xmax.
    ! Normalisation with respect to 0:
    ! pn(i+1)*n(i+1)=(x-a(i)) pn(i) n(i)-b(i) pn(i-1) n(i-1)
    ! n(i+1)=|-a(i)n(i)pn(i,x=0)-b(i)n(i-1)pn(i,x=0)|
    ! But we normalise to <P^2>


    w(x)=x**alpha*exp(-sign*x*x) ! inline function for weight
    sechb(x)=2.*x/(1.+x*x)
    sech(x)=sechb(exp(-x))
    sech2(x)=sech(abs(x))**2
    tanhb(x)=2*x/(1+x)
    tanhp1(x)=tanhb(exp(2*x))   ! tanh(x)+1 for x<0 approx. 2*exp(2*x) for x<<-1
    tanhm1(x)=-tanhb(exp(-2*x)) ! tanh(x)-1 for x>0
    leftnew=0
    rightnew=0

    if (n==0) return
    nn=n
    vb=present(verbose)

    t0=-log(-log(epsmori/2.)/mori)
    t1=-t0
    t00=t0
    t10=t1
    if (vb) then
       ngauss=n+xmax**2+xmax*sqrt(-2*log(epsgauss))-log(epsgauss)/3
       ngauss=n+ceiling(10*xmax)+ceiling(alpha/2)
       print *,'ngauss=',ngauss
       allocate(gw(ngauss),gp(ngauss),gw2(ngauss*2),gp2(ngauss*2))
       call gauss_legendre(xmin,xmax,gp,gw,ngauss)
       call gauss_legendre(xmin,xmax,gp2,gw2,2*ngauss)
    endif
    if (xmax**2 > Log(1e300)) then
       if (vb) print *,'Careful, the squared polynomials will be of order exp(xmax^2) at&
            & xmax. xmax should be smaller than',sqrt(log(1e300))
    endif
    if (vb) print *,'thp,e',tanhp1(log(epsmori/2.)/2),tanhp1(mori*sinh(t0)),epsmori

    nintervals=max(128,2*n) ! number of elementary intervals between integration points
    eps=epsbase*sqrt(1.*nintervals) !*3*n
    dt=(t1-t0)/nintervals
    !preliminary values
    a1(1)=0
    b1(1)=0
    c1(1)=1
    ! allocate enough space, so that we do not have to reallocate frequently
    minalloc=-10*nsafe-9
    maxalloc=10*nsafe+9+2*nintervals
    i0=0
    allocate(mp(minalloc:maxalloc),mw(minalloc:maxalloc),poly1(minalloc:maxalloc),poly(minalloc:maxalloc))
    if (vb) then
       mw=-500 ! for debugging purposes
       mp=-500
       poly=-500
       poly1=-500
    end if
    do i=1,n
       ! first compute g(i)=<P(i)^2> and <x P(i)^2> using double exponential
       ! integration.
       ! We transform  the interval x=0 ... xmax into
       ! x=xmax/2*(1+tanh(mori*sinh(t))) with t=-infty ... infty.
       !  =xmax/2*(1+z) 
       ! dx/dt=xmax/2*mori*cosh(t) (1+tanh^2(mori*sinh(t)))

       ! guess the minimum x as xmax*epsmori and the maximum as
       ! xmax*(1-epsmori)
       ! ==> 1st guess t integral has actually bounds t=+-log(-mori*log(2*eps))
       !   (see above definition of t0 and t1)
       ! 1st guess for number of required points is 2n+1.

       ! Another thing: we accumulate the sum from its small ends, so that
       ! numerical precision errors remain as small as possible.

       ! here we have the preliminary starting values
       !a1(i)=a(i) 
       !b1(i)=sqrt(b(i-1))
       !c1(i)=1

       s1=0
       s1x=0
       do k=1,nintervals/2
          ia=2*k+i0
          if (i==1) then
             t=t0+dt*k
             dxdt=.5*(xmax-xmin)*sech2(mori*sinh(t))*mori*cosh(t)
             ! x approx. xmax*exp(-mori*exp(-t)) for t<-3 or so
             x=.5*(xmax-xmin)*tanhp1(mori*sinh(t))+xmin
             ! calculate polynomial at x
             mp(ia)=x
             mw(ia)=dxdt*w(x)
             poly1(ia)=0
             pi=c1(1)
             if (vb) then; if (poly(ia)/=-500) then
                print *,'serious error h',ia,minalloc,maxalloc,i0,k
                stop
             end if; end if
          else
             x=mp(ia)
             pi=(c1(i)*x-a1(i))*poly(ia)-b1(i)*poly1(ia)
             poly1(ia)=poly(ia)
          endif
          poly(ia)=pi

          !if (vb) print *,'left: i,k,t,x,pi,dxdt,w',i,k,t,x,pi,dxdt,w(x)
          pi=pi*pi*mw(ia)
          s1=s1+pi
          s1x=s1x+pi*x
       enddo
       if (vb) print *,'left integral i,s,sx',i,s1,s1x
       s2=0
       s2x=0
       do k=1,(nintervals-1)/2
          ia=2*(nintervals-k)+i0
          if (i==1) then
             t=t1-dt*k
             !z=tanh(mori*sinh(t))
             dxdt=.5*(xmax-xmin)*sech2(mori*sinh(t))*mori*cosh(t)
             !x=.5*(xmax-xmin)*(1+z)+xmin
             x=.5*(xmax-xmin)*tanhm1(mori*sinh(t))+xmax
             ! calculate polynomial at x
             mp(ia)=x
             mw(ia)=dxdt*w(x)
             poly1(ia)=0
             pi=c1(1)
             if (vb) then; if (poly(ia)/=-500) then
                print *,'serious error g',ia,minalloc,maxalloc,i0,k
                stop
             end if; end if
          else
             x=mp(ia)
             pi=(c1(i)*x-a1(i))*poly(ia)-b1(i)*poly1(ia)
             poly1(ia)=poly(ia)
          end if
          poly(ia)=pi
          !if (vb) print *,'right: i,k,t,x,pi,dxdt,w',i,k,t,x,pi,dxdt,w(x)
          pi=pi*pi*mw(ia)
          s2=s2+pi
          s2x=s2x+pi*x
       enddo
       ! nintervals/2 + (nintervals-1)/2 = nintervals-1, for odd or even
       if (vb) print *,'right integral i,s,sx',i,s2,s2x

       ! now see whether we need to extend the integration domain, or whether
       ! it is accurate enough: We do nsafe extra points at both ends.
       makenew=.false.
       do
          s3=0
          s3x=0
          do k=-nsafe,0
             ia=2*k+i0
             if (makenew .or. i==1) then
                t=t0+dt*k
                dxdt=.5*(xmax-xmin)*sech2(mori*sinh(t))*mori*cosh(t)
                ! x approx. xmax*exp(-mori*exp(-t)) for t<-3 or so
                x=.5*(xmax-xmin)*tanhp1(mori*sinh(t))+xmin
                ! calculate polynomial at x
                mp(ia)=x
                mw(ia)=dxdt*w(x)
                p1=0
                pi=c1(1)
                do j=2,i
                   p0=pi
                   pi=(c1(j)*x-a1(j))*pi-b1(j)*p1
                   p1=p0
                enddo
                poly1(ia)=p1

                if (vb) then; if (poly(ia)/=-500) then
                   print *,'serious error e',ia,minalloc,maxalloc,i0,k
                   stop
                end if; end if
             else
                x=mp(ia)
                pi=(c1(i)*x-a1(i))*poly(ia)-b1(i)*poly1(ia)
                poly1(ia)=poly(ia)
             endif
             poly(ia)=pi
             !           if (vb) print *,'Extension left test',k,s1,s3,t,x,pi,w(x),dxdt
             pi=pi*pi*mw(ia)
             s3=s3+pi
             s3x=s3x+pi*x
          enddo
          !        if (vb) print *,'Extension left test',k,s1,s3,t,x,pi,w(x),dxdt
          if (abs(s3)>abs(s1+s2)*eps .or. abs(s3x)>abs(s1x+s2x)*eps) then
             ! must extend t-domain to the left
             ! this extension should never stop for alpha=-1. But it does. Why?
             s1=s1+s3
             s1x=s1x+s3x
             if (vb) print *,'Extending t-domain to left',t0,t0-(nsafe+1)&
                  *dt,s1,s3,s1x,s3x,s2,s2x
             if (i0-4*nsafe-2<minalloc) then
                if (i==1) nintervals=nintervals-nsafe-1 ! to prevent error with not yet defined
                !right boundary region.
                call realloc(1,mw,mp,poly,poly1,minalloc,maxalloc,i0,nintervals,nn,vb)
                if (i==1) nintervals=nintervals+nsafe+1
             endif
             t0=t0-(nsafe+1)*dt
             nintervals=nintervals+(nsafe+1)
             i0=i0-2*(nsafe+1)
             makenew=.true.
             leftnew=leftnew+nsafe+1
          else
             ! if (vb) print *,'NOT Extending t-domain to left',t0,t0-(nsafe+1)&
             !    *dt,s1,s3,s1x,s3x,s2,s2x,'dxdt',dxdt,poly(i0),mw(i0),x&
             !    ,w(x),alpha
             exit
          endif
       enddo
       makenew=.false.
       do
          s3=0
          s3x=0
          do k=-nsafe,0
             ia=2*(nintervals-k)+i0
             if (makenew .or. i==1) then
                t=t1-dt*k
                !z=tanh(mori*sinh(t))
                dxdt=.5*(xmax-xmin)*sech2(mori*sinh(t))*mori*cosh(t)
                !x=.5*(xmax-xmin)*(1+z)+xmin
                x=.5*(xmax-xmin)*tanhm1(mori*sinh(t))+xmax
                ! calculate polynomial at x
                mp(ia)=x
                mw(ia)=dxdt*w(x)
                p1=0
                pi=c1(1)
                do j=2,i
                   p0=pi
                   pi=(c1(j)*x-a1(j))*pi-b1(j)*p1
                   p1=p0
                enddo
                poly1(ia)=p1

                if (vb) then; if (poly(ia)/=-500) then
                   print *,'serious error f',i,ia,minalloc,maxalloc,i0,k
                   stop
                end if; end if
             else
                x=mp(ia)
                pi=(c1(i)*x-a1(i))*poly(ia)-b1(i)*poly1(ia)
                poly1(ia)=poly(ia)
             endif
             poly(ia)=pi
             pi=pi*pi*mw(ia)
             s3=s3+pi
             s3x=s3x+pi*x
          enddo
          if (abs(s3)>abs(s1+s2)*eps .or. abs(s3x)>abs(s1x+s2x)*eps) then
             ! must extend t-domain to the right
             s1=s1+s3
             s1x=s1x+s3x
             if (vb) print *,'Extending t-domain to right',t1,t1+(nsafe+1)&
                  *dt,s1,s3,s1x,s3x,s2,s2x!,'dxdt',dxdt,poly(2*nintervals&
             !+i0),mw(2*nintervals+i0),x,w(x)
             if (i0+2*(nintervals+2*nsafe+1)>maxalloc) &
                  call realloc(1,mw,mp,poly,poly1,minalloc,maxalloc,i0,nintervals,nn,vb)
             t1=t1+(nsafe+1)*dt
             nintervals=nintervals+(nsafe+1)
             makenew=.true.
             rightnew=rightnew+nsafe+1
          else
             exit
          endif
       enddo
       ! add the two halves together and apply dt:
       s3=dt*(s1+s2)
       s3x=dt*(s1x+s2x)
       if (vb) print *,'complete integrals',i,s3,s3x

       ! Now we check, whether we have to increase the number of points
       ! For this doubling of the integration points we take one more interval
       ! than above
       ! Why do we not get arbitrarily good accuracy, when we take more and
       ! more points?
       ! shouldn't also the numerical errors decrease statistically due to
       ! averaging? 
       ! We're summing only positive numbers here.
       ! agreed, there may be shot noise in the polynomial itself.
       ! in addition: the dominant effort is probably in the calculation of the
       ! vertices, not in the polynomial. The variance increases propto n.
       if (.true. .or. nintervals < 1000 .or. nintervals < maxmem*n) then
          do
             s1=0
             s1x=0
             ! The extent here goes only for nsafe points outside our
             ! interior intervals.

             do k=1-nsafe,(nintervals+1)/2
                ia=2*k-1+i0
                if (k<=leftnew .or. i==1) then
                   t=(t0-.5*dt)+dt*k
                   dxdt=.5*(xmax-xmin)*sech2(mori*sinh(t))*mori*cosh(t)
                   x=.5*(xmax-xmin)*tanhp1(mori*sinh(t))+xmin
                   ! calculate polynomial at x
                   mp(ia)=x
                   mw(ia)=dxdt*w(x)
                   p1=0
                   pi=c1(1)
                   do j=2,i
                      p0=pi
                      pi=(c1(j)*x-a1(j))*pi-b1(j)*p1
                      p1=p0
                      if (vb .and. k==1-nsafe) print *,'refiningxx',x,k,pi,p1,j
                   enddo
                   if (vb .and. k==1-nsafe) print *,'refiningyy',c1(1:i),'a',a1(2:i),'b&
                        &',b1(2:i)
                   poly1(ia)=p1

                   if (vb) then; if (poly(ia)/=-500) then
                      print *,'serious error a',ia,minalloc,maxalloc,i0,k
                      print *,'poly(',i0-nsafe*2,')',poly(i0-nsafe*2:i0)
                      stop
                   end if; end if
                else
                   if (vb) then; if (mp(ia)==-500) then
                      print *,'serious error c',ia,minalloc,maxalloc,i0,k
                      print *,'poly(',i0-nsafe*2,')',poly(i0-nsafe*2:i0)
                      stop
                   end if; end if
                   x=mp(ia)
                   pi=(c1(i)*x-a1(i))*poly(ia)-b1(i)*poly1(ia)
                   poly1(ia)=poly(ia)
                endif
                poly(ia)=pi
                !if (vb) print *,'left:
                !i,k,t,x,pi,dxdt,w',i,k,t,x,pi,dxdt,w(x)
                if (vb) then
                   if (mw(ia)==-500) then
                      print *,'error xxx'
                      stop
                   endif
                end if
                pi=pi*pi*mw(ia)
                s1=s1+pi
                s1x=s1x+pi*x
                !if (vb) print *,'refining',k,poly(ia-1:ia+1),mw(ia-1:ia+1),mp(ia-1:ia+1)
             enddo
             s2=0
             s2x=0
             do k=1-nsafe,nintervals/2
                ia=2*(nintervals-k)+1+i0
                if (k<=rightnew .or. i==1) then
                   t=(t1+.5*dt)-dt*k
                   !z=tanh(mori*sinh(t))
                   dxdt=.5*(xmax-xmin)*sech2(mori*sinh(t))*mori*cosh(t)
                   !x=.5*(xmax-xmin)*(1+z)+xmin
                   x=.5*(xmax-xmin)*tanhm1(mori*sinh(t))+xmax
                   ! calculate polynomial at x
                   mp(ia)=x
                   mw(ia)=dxdt*w(x)
                   p1=0
                   pi=c1(1)
                   do j=2,i
                      p0=pi
                      pi=(c1(j)*x-a1(j))*pi-b1(j)*p1
                      p1=p0
                   enddo
                   poly1(ia)=p1
                   if (vb) then; if (poly(ia)/=-500) then
                      print *,'serious error b',ia,minalloc,maxalloc,i0,k
                      stop
                   end if; end if
                else
                   if (vb) then; if (mp(ia)==-500) then
                      print *,'serious error d',ia,minalloc,maxalloc,i0,k
                      print *,'poly(',i0-nsafe*2,')',poly(i0-nsafe*2:i0)
                      print *,'poly(',ia-nsafe,')',poly(ia-nsafe:ia+nsafe)
                      print *,'mp(',ia-nsafe,')',mp(ia-nsafe:ia+nsafe)
                      stop
                   end if; end if
                   x=mp(ia)
                   pi=(c1(i)*x-a1(i))*poly(ia)-b1(i)*poly1(ia)
                   poly1(ia)=poly(ia)
                end if
                poly(ia)=pi
                !if (vb) print *,'right: i,k,t,x,pi,dxdt,w',i,k,t,x,pi,dxdt,w(x)
                pi=pi*pi*mw(ia)
                s2=s2+pi
                s2x=s2x+pi*x
             enddo
             ! nintervals/2 + (nintervals+1)/2 = nintervals, for odd or even
             ! so now we have in toto nintervals*2-1 points.
             s1=(s1+s2)*dt
             s1x=(s1x+s2x)*dt
             if (vb) print *,'complete intermediate integrals',i,s1,s2*dt,s1x
             s1=.5*(s1+s3)
             s1x=.5*(s1x+s3x)
             !s1 and s1x contain the up to now most accurate values.
             if (abs(s1-s3)>eps*abs(s3) .or. abs(s1x-s3x)>eps*abs(s3x)) then
                ! need to increase number of intervals
                if (vb) print *,'Increasing nintervals',dt,s1,s3,abs(s1/s3&
                     -1),s1x,s3x,nintervals,2*nintervals+1
                s3=s1
                s3x=s1x
                if (vb) print *,'old dt, refinements',dt,(t1-t0)/nintervals
                !nintervals=2*(nintervals+2*nsafe)-2*nsafe
                nintervals=2*(nintervals+nsafe)
                t0=t0-.5*dt*nsafe
                t1=t1+.5*dt*nsafe
                eps=epsbase*sqrt(1.*nintervals) !*3*n
                dt=(t1-t0)/nintervals
                i0=i0-nsafe ! new i0 index
                !print *,'prerefpoly',poly(i0-2*nsafe:i0+nsafe)
                call realloc(2,mw,mp,poly,poly1,minalloc,maxalloc,i0,nintervals,nn,vb)
                leftnew=nintervals !make all the half points new
                rightnew=nintervals
                if (vb) print *,'new dt, refinements',dt,log((t10-t00)/(dt*2*n))/log(2.)
             else
                leftnew=-nsafe !make nothing new next time.
                rightnew=-nsafe
                exit
             endif
          enddo
       else
          s1=s3
          s1x=s3x
       endif
       ! now we have the hopefully accurate integrals in s1 and s1x.
       ! that's it!
       ! compare with gauss
       if (vb) then
          s4=0
          s4x=0
          do k=ngauss,1,-1
             x=gp(k)
             ! calculate polynomial at x
             p1=0
             pi=c1(1)

             do j=2,i
                p0=pi
                pi=(c1(j)*x-a1(j))*pi-b1(j)*p1
                p1=p0
             enddo

             pi=pi*pi*w(x)

             s4=s4+pi*gw(k)
             s4x=s4x+pi*x*gw(k)
          enddo
          s5=0
          s5x=0
          do k=2*ngauss,1,-1
             x=gp2(k)
             ! calculate polynomial at x
             p1=0
             pi=c1(1)

             do j=2,i
                p0=pi
                pi=(c1(j)*x-a1(j))*pi-b1(j)*p1
                p1=p0
             enddo

             pi=pi*pi*w(x)

             s5=s5+pi*gw2(k)
             s5x=s5x+pi*x*gw2(k)
          enddo
          print *,'gauss comparison',i,s4/s3-1,s4x/s3x-1,s5/s4-1,s5x/s4x-1,s5/s3-1,s5x/s3x-1,ngauss,nintervals
       endif !vb
       logg(i)=log(s1)+logg(i)
       !a(i)=s1x/s1
       !b(i)=s1
       !final settings a1(i)=a(i-1)/b(i)^.5  b1(i)=(b(i-1)/b(i))^.5
       ! c1(i)=1/b(i)^.5
       bsq(i)=sqrt(s1)
       c1(i)=1/bsq(i)
       a1(i)=a1(i)*c1(i)
       b1(i)=b1(i)*c1(i)
       if (vb) then
          if (i0-nsafe*2<minalloc .or. i0+2*(nintervals+nsafe)>maxalloc) then
             print *,'Utmost error',i0,nsafe,nintervals,minalloc,maxalloc
             stop
          endif
          do k=i0-nsafe*2,i0+2*(nintervals+nsafe)
             if (poly(k)==-500) then
                print *,'Still error',k,i0,nsafe,nintervals,minalloc,maxalloc
                stop
             endif
          end do
       endif
       poly(i0-nsafe*2:i0+2*(nintervals+nsafe))=poly(i0-nsafe*2:i0+2*(nintervals+nsafe))*c1(i)
       ! This is only the recalculation of poly to avoid different numerical
       ! values than the old version of this routine.
!!$       do k=i0-nsafe*2,i0+2*(nintervals+nsafe)
!!$          p1=0
!!$          pi=c1(1)
!!$          x=mp(k)
!!$          do j=2,i
!!$             p0=pi
!!$             pi=(c1(j)*x-a1(j))*pi-b1(j)*p1
!!$             p1=p0
!!$          enddo
!!$          poly(k)=pi
!!$       enddo
       a(i)=s1x/s1
       if (i<n) then 
          !preliminary settings a1(i+1)=a(i)  ,  b1(i+1)=b(i)^.5  , logg.
          a1(i+1)=s1x/s1
          if (i>1) then
             b1(i+1)=sqrt(s1)
          else
             b1(i+1)=0
          end if
          c1(i+1)=1
          logg(i+1)=logg(i)
       end if
       if (vb) print *,'i,a1,b1,c1,logg',i,a1(i),b1(i),c1(i),logg(i),'nintervals=',nintervals
    enddo
    if (vb) print *,'hhn used mem',2*(nsafe+nintervals)
    deallocate(mw,mp,poly,poly1) !deallocation may be necessary because these are pointers.
  end subroutine half_hermite_norm

  subroutine gauss_nodes_functions(n,a,bsq,sp,ortho_f)
    ! given the recursion coeffs a,bsq calculates the Gauss
    ! points sp and the orthogonal functions
    ! Poly(sp)*sqrt(w(sp)) at these points  in ortho_f,
    ! where w(x) is the weight function of the set of polynomials
    use, intrinsic :: ieee_exceptions
    implicit none
    integer,intent(in) :: n
    real,intent(in) :: a(n+1),bsq(n+1)
    real,intent(out) :: sp(n)
    real,optional,target,intent(out) :: ortho_f(n,n)
    real,pointer :: projsteen(:,:)
    
    integer info
    integer,allocatable,dimension(:) :: ifail,isuppz,iwork
    real,allocatable,dimension(:) :: work
    integer liwork,lwork
    integer i,m
    logical invmode,zeromode !ieee_exception mode storage

    if (present(ortho_f)) then
       projsteen=>ortho_f
    else
       allocate(projsteen(n,n))
    end if

    !dstevr: liwork>=1,10*n
    liwork=10*n
    !dstevr: lwork >=1,20*n
    lwork=20*n
    !dstevr: isuppz: dim>=2*n
    
    !dstein: ifail: dim>=n
    !dstein: work: dim>=5*n
    !dstein: iwork: dim>=n
    
    allocate(ifail(n),isuppz(2*n),iwork(liwork),work(lwork))

    ! Set up Jacobi-Matrix. (tridiagonal, symmetric)
    ! Note: J*(pn_i-1(x))=(xp_(i-1)(x)-delta_i,n J_n,n+1 pn_n(x)) with the
    ! vector pn_(i-1)(x),i=1...n for any x.
    ! pn_i is the normalised polynomial. pn_i=p_i/sqrt(gamma_i)
    ! i.e. x pn_i-1(x)=sum_j J_ij pn_j-1(x) for i=1...n-1 and
    !      x pn_n-1(x)=sum_j J_nj pn_j-1(x) + J_n,n+1 Pn_n(x)
    
    call ieee_get_halting_mode(ieee_invalid,invmode)
    call ieee_get_halting_mode(ieee_divide_by_zero,zeromode)
    call ieee_set_halting_mode(ieee_invalid,.false.)
    call ieee_set_halting_mode(ieee_divide_by_zero,.false.)
    call dstevr('N','A',n,a(1:n),bsq(2:n),0.,0.,0,0,0.,m,sp,projsteen,n,isuppz,work,lwork,iwork,liwork,info)

    if (info/=0) then
       print *,'module half_hermite: gaus_nodes_weights: dstevr: info=',info,m,n
       stop
    end if
    
    !Using dstein for the eigenvectors instead of dstevr itself, since they are accurate
    !even in the face of very small components.
    call dstein(n,a(1:n),bsq(2:n),n,sp,(/(1,i=1,n)/),(/(n,i=1,n)/),projsteen,n,work,iwork,ifail,info)
    call ieee_set_halting_mode(ieee_invalid,invmode)
    call ieee_set_halting_mode(ieee_divide_by_zero,zeromode)
    call ieee_set_flag(ieee_all,.false.)
    ! the eigenvectors are normalised to their square sums
    ! The eigenvectors in projsteen are projsteen(i,j)=pn_(i)(x_j) sqrt(w_j) n_i
    ! with a normalisation factor n_i and the Gauss weights w_j
    ! They are normalised, such that
    ! sum_j n_i^2pn_{i}(x_j)^2w_j=1 for any i, and
    ! sum_i n_i^2pn_{i}(x_j)^2w_j=1 for any j
    ! Since pn_(1)=1/sqrt(int w(x) dx)=c11 for the weight function w(x),
    ! we have the Gauss weights as
    
    
    if (.not. present(ortho_f)) deallocate(projsteen)
    
  end subroutine gauss_nodes_functions
  
  !interface function for cgyro_init_manager: (cf also pseudo_spec_lib.f90)
  subroutine pseudo_maxwell_pliocene(n,e_max,x,w_e,d1,alpha,datafile)
    
    implicit none

    character (len=*), intent(in) :: datafile
    optional :: datafile

    integer, intent(in) :: n
    real, intent(in) :: e_max,alpha

    real, intent(out) :: x(n)
    real, intent(out) :: w_e(n)
    real, intent(out) :: d1(n,n)
    
    real, dimension(:), allocatable :: a1,b1,c1,lg,a,bsq
    real, dimension(:,:), allocatable :: projsteen,deriv(:,:)

    real xmax,v
    integer i

    xmax=sqrt(e_max)
    
    allocate(a1(n+1),b1(n+1),c1(n+1),a(n+1),bsq(n+1),lg(n+1))
    call  half_hermite_norm(n+1,0.,xmax,1.,alpha,a1,b1,c1,a,bsq,lg)
    ! calculate steen points and weights
    
    allocate(projsteen(n,n),deriv(n,n))
    call gauss_nodes_functions(n,a,bsq,x,projsteen)
    ! recalculate weights for weight function x**2 expected by the code.
    ! relative size of remainder integral to recreate balancing the
    ! integral by last element
    v=exp(-e_max)*2*xmax/sqrt(pi1)+erfc(xmax) !(1-erf(xmax))
    w_e=projsteen(1,:)**2*x**(2-alpha)
    w_e=w_e/sum(w_e)*(1-v)
    w_e(n)=w_e(n)+v

    
    ! now x contains the Gauss points as velocities
    !   w_e contains the Gauss weights.
    !   w_e=projsteen(1,:)**2/c1(1)**2 true Gauss integral weights
    ! projsteen is an orthogonal matrix:
    !       projsteen^-1=projsteen^T
    ! 
    ! Steen polynomial Pi(x_j)=projsteen(i,j)/sqrt(w_j)
    !  P1(x_j)=c1(1) (const.)
    ! With Dw=diag(sqrt(w_j)):
    ! projsteen*Dw converts polynomials given on the vertices into polynomial space
    ! Dw^-1*projsteen^-1=Dw^-1*projsteen^T=(projsteen*Dw^-1)^T converts
    ! from polynomial space into values on the vertices.
    
    deriv(1,:)=0
    if (n>=2) then
       i=2
       deriv(i,:)=c1(i)*projsteen(i-1,:)
    endif
    do i=3,n
       deriv(i,:)=c1(i)*projsteen(i-1,:)+(c1(i)*x-a1(i))*deriv(i-1,:)-b1(i)*deriv(i-2,:)
    end do
    ! deriv^T converts from polynomial space to derivatives on the vertices times sqrt(w_j)
    do i=1,n
       deriv(:,i)=deriv(:,i)/projsteen(1,i) !/sqrt(w_e(i)) /c1(1) 
    end do
    ! deriv^T now converts from polynomial space to derivatives on the vertices
    ! d1=deriv^T*projsteen*Dw
    call dgemm('t','n',n,n,n,1.,deriv,n,projsteen,n,0.,d1,n)
    do i=1,n
       d1(:,i)=d1(:,i)*projsteen(1,i) !*sqrt(w_e(i))*c1(1)
    end do
    d1=d1*xmax !code wants derivatives with respect to v multiplied with xmax
    x=x*x !code wants energy points not velocity points.
    if (present(datafile)) then
       open(unit=1,file=trim(datafile))
       do i=1,n
          write(1,2) x(i),w_e(i)
       enddo
       do i=1,n
          write(1,2) d1(i,:)
       enddo
       close(1)
    end if
2   format (*(EN25.16))
  end subroutine pseudo_maxwell_pliocene
end module half_hermite
