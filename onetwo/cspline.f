      subroutine cspline (x, y, n, yp1, ypn, y2)
cProlog

c
      implicit none
c
c
      integer i, n, nmax, k
      doubleprecision p, sig, qn, un
      parameter (nmax = 129)
      doubleprecision x(n), y(n), y2(n), u(nmax), yp1, ypn
c
      if (n .gt. nmax) call STOP (
     &'subroutine CSPLINE: dimensional error', 201)
c
      if (yp1 .gt. 0.99d30) then
      y2(1) = 0.d0
      u(1) = 0.d0
      else
      y2(1) = -0.5d0
      u(1) = (3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do i=2,n-1
      sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))
      p = sig*y2(i-1)+2.d0
      y2(i) = (sig-1.d0)/p
      u(i) = (6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1)) /(x(i)-x(i
     &-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      enddo
      if (ypn .gt. 0.99d30) then
      qn = 0.d0
      un = 0.d0
      else
      qn = 0.5d0
      un = (3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.d0)
      do k=n-1,1,-1
      y2(k) = y2(k)*y2(k+1)+u(k)
      enddo
      return
c
      end
