
      SUBROUTINE icsicu1 (x, y, nx, bpar, c, ic, ier)

! ----------------------------------------------------------------------
!
!   modified version of IMSL subroutine named ICSICU
!
!   purpose             - interpolatory approximation by cubic splines
!                           with arbitrary second derivative end
!                           conditions.
!
!   usage               - call icsicu1 (x, y, nx, bpar, c, ic, ier)
!
!   arguments    x      - vector of length nx containing the abscissae
!                           of the nx data points (x(i),y(i)) i = 1,...,
!                           nx. (input) x must be ordered so that
!                           x(i) .lt. x(i+1).
!                y      - vector of length nx containing the ordinates
!                           (or function values) of the nx data points.
!                           (input)
!                nx     - number of elements in x and y. (input) nx
!                           must be .ge. 2.
!                bpar   - vector of length 4 containing the end
!                           condition parameters. (input)
!                           2.0*spp(1)+bpar(1)*spp(2) = bpar(2),
!                           bpar(3)*spp(nx-1)+2.0*spp(nx) = bpar(4),
!                           where spp(i) = second derivative of the
!                           cubic spline function s evaluated at x(i).
!                c      - spline coefficients. (output) c is an nx-1 by
!                           3 matrix. the value of the spline
!                           approximation at t is
!                           s(t) = ((c(i,3)*d+c(i,2))*d+c(i,1))*d+y(i)
!                           where x(i) .le. t .lt. x(i+1) and
!                           d = t-x(i).
!                ic     - row dimension of matrix c exactly as
!                           specified in the dimension statement in
!                           the calling program. (input)
!                ier    - error parameter. (output)
!                         terminal error
!                           ier = 129, ic is less than nx-1
!                           ier = 130, nx is less than 2.
!                           ier = 131, input abscissa are not ordered
!                             so that x(1) .lt. x(2) ... .lt. x(nx).
!
!   copyright           - 1978 by imsl, inc. all rights reserved.
!
! ----------------------------------------------------------------------
!
      USE nrtype,                                 ONLY : DP,I4B
      IMPLICIT NONE
!     specifications for arguments
      INTEGER            nx,ic,ier
      REAL*8             x(nx),y(nx),bpar(4),c(ic,3)



!     specifications for local variables
!
      INTEGER(I4B)       i,j,nxm1
      REAL(DP)           dx,dxj,dxjp1,dxp,dyj,dyjp1,half,one,pj, &
                         six,sixi,two,yppa,yppb,zero
      EQUIVALENCE        (dxj,yppb),(pj,sixi),(dxjp1,yppa)
      DATA               zero/0.0_DP/,half/0.5_DP/,one/1.0_DP/, &
                         two/2.0_DP/,six/6.0_DP/
!
      ier = 0
!
!     check error conditions
!
      nxm1 = nx-1
      IF (ic .LT. nxm1)  go to 30
      IF (nx .LT. 2   )  go to 35
      IF (nx .EQ. 2   )  go to 10
!
!     compute coefficients and right hand side of the tridiagonal
!     system defining the second derivatives of the spline interpolant for (x,y)
!
!     c(j,1) = lambda(j)
!     c(j,2) = mu(j)
!     c(j,3) = d(j)
!
      dxj = x(2)-x(1)
      IF (dxj .LE. zero)  go to 40
      dyj = y(2)-y(1)
      DO 5 j=2,nxm1
         dxjp1 = x(j+1)-x(j)
         IF (dxjp1 .LE. zero)  go to 40
         dyjp1 = y(j+1)-y(j)
         dxp = dxj+dxjp1
         c(j,1) = dxjp1/dxp
         c(j,2) = one-c(j,1)
         c(j,3) = six*(dyjp1/dxjp1-dyj/dxj)/dxp
         dxj = dxjp1
         dyj = dyjp1
    5 CONTINUE
!
!     factor the tridiagonal matrix and solve for u
!
!     c(j,2)  = u(j)
!     c(j,1)  = q(j)
!     bpar(1) = lambda(1)
!     bpar(2) = d(1)
!     bpar(3) = mu(nx)
!     bpar(4) = d(nx)
!
   10 c(1,1) = -bpar(1)*half
      c(1,2) = bpar(2)*half
      IF (nx .EQ. 2)  go to 20
      DO 15 j=2,nxm1
         pj = c(j,2)*c(j-1,1)+two
         c(j,1) = -c(j,1)/pj
         c(j,2) = (c(j,3)-c(j,2)*c(j-1,2))/pj
   15 CONTINUE
!
!     solve for cubic coefficients of spline interpolant
!     c(j,1), c(j,2), and c(j,3)
!
   20 yppb = (bpar(4)-bpar(3)*c(nxm1,2))/(bpar(3)*c(nxm1,1)+two)
      sixi = one/six
      DO 25 i=1,nxm1
         j = nx-i
         yppa = c(j,1)*yppb+c(j,2)
         dx = x(j+1)-x(j)
         c(j,3) = sixi*(yppb-yppa)/dx
         c(j,2) = half*yppa
         c(j,1) = (y(j+1)-y(j))/dx-(c(j,2)+c(j,3)*dx)*dx
         yppb = yppa
   25 CONTINUE
      go to 9005
   30 ier = 129
      go to 9000
   35 ier = 130
      go to 9000
   40 ier = 131
!
 9000 CONTINUE
      

 9005 RETURN
!
      END SUBROUTINE icsicu1

      SUBROUTINE icsevu1 (x, y, nx, c, ic, u, s, m, ier)
! --- modified version of IMSL subroutine named ICSEVU
!
! ----------------------------------------------------------------------
!   purpose             - evaluation of a cubic spline
!
!   usage               - call icsevu1 (x, y, nx, c, ic, u, s, m, ier)
!
!   arguments    x      - vector of length nx containing the abscissae
!                           of the nx data points (x(i),y(i)) i = 1,...,
!                           nx (input). x must be ordered so that
!                           x(i) .lt. x(i+1).
!                y      - vector of length nx containing the ordinates
!                           (or function values) of the nx data points
!                           (input).
!                nx     - number of elements in x and y (input).
!                           nx must be .ge. 2.
!                c      - spline coefficients (input). c is an nx-1 by
!                           3 matrix.
!                ic     - row dimension of matrix c exactly as
!                           specified in the dimension statement
!                           in the calling program (input).
!                           ic must be .ge. nx-1
!                u      - vector of length m containing the abscissae
!                           of the m points at which the cubic spline
!                           is to be evaluated (input).
!                s      - vector of length m (output).
!                           the value of the spline approximation at
!                           u(i) is
!                           s(i) = ((c(j,3)*d+c(j,2))*d+c(j,1))*d+y(j)
!                           where x(j) .le. u(i) .lt. x(j+1) and
!                           d = u(i)-x(j).
!                m      - number of elements in u and s (input).
!                ier    - error parameter (output).
!                         warning error
!                           ier = 33, u(i) is less than x(1).
!                           ier = 34, u(i) is greater than x(nx).
!
!                           ********************************************
!                           output of warning errors has been suppressed
!                           ********************************************
!
!   notation            - information on special notation and
!                           conventions is available in the manual
!                           introduction or through IMSL routine uhelp
!
!   remarks  1.  the routine assumes that the abscissae of the nx
!                data points are ordered such that x(i) is less than
!                x(i+1) for i = 1,...,nx-1. no check of this condition
!                is made in the routine. unordered abscissae will cause
!                the algorithm to produce incorrect results.
!            2.  the routine generates two warning errors. one error
!                occurs if u(i) is less than x(1), for some i in the
!                the interval (1,m) inclusively. the other error occurs
!                if u(i) is greater than x(nx), for some i in the
!                interval (1,m) inclusively.
!            3.  the ordinate y(nx) is not used by the routine. for
!                u(k) .gt. x(nx-1), the value of the spline, s(k), is
!                given by
!                 s(k) = ((c(nx-1,3)*d+c(nx-1,2))*d+c(nx-1,1))*d+y(nx-1)
!                where d = u(k)-x(nx-1).
!
!   copyright           - 1978 by imsl, inc. all rights reserved.
!
! ----------------------------------------------------------------------
!
      USE nrtype,                                 ONLY : DP,I4B
      IMPLICIT NONE

!     specifications for arguments
      INTEGER(I4B)            nx,ic,m,ier
      REAL(DP)             x(nx),y(nx),c(ic,3),u(m),s(m)


!     specifications for local variables
      INTEGER            i,jer,ker,nxm1,k
      REAL*8             d,dd,zero
      DATA               i/1/, zero/0.0_DP/
!
!     first executable statement
!
      jer = 0
      ker = 0
      IF (m .LE. 0)  go to 9005
      nxm1 = nx-1
      IF (i .GT. nxm1)  i = 1
!
!     evaluate spline at m points
!
      DO 40 k=1,m
!
!        find the proper interval
!
         d = u(k)-x(i)
         IF (d) 5, 25, 15
    5    IF (i .EQ. 1)  go to 30
         i = i-1
         d = u(k)-x(i)
         IF (d) 5, 25, 20
   10    i = i+1
         d = dd
   15    IF (i .GE. nx)  go to 35
         dd = u(k)-x(i+1)
         IF (dd .GE. zero)  go to 10
         IF ( d .EQ. zero)  go to 25
!
!        perform evaluation
!
   20    s(k) = ((c(i,3)*d+c(i,2))*d+c(i,1))*d+y(i)
         go to 40
   25    s(k) = y(i)
         go to 40
!
!        u(k) < x(1)
!
   30    jer = 33
         go to 20
!
!        u(k) > x(nx)
!
   35    IF (dd .GT. zero)  ker = 34
         d = u(k) - x(nxm1)
         i = nxm1
         go to 20
!
   40 CONTINUE
!
      ier = MAX0 (jer, ker)

 9005 RETURN
!
      END SUBROUTINE icsevu1
