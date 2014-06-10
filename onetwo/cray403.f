      subroutine exp1 (itran)
c
      USE param
      USE soln
      USE numbrs
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      character rcs_id*63
      save      rcs_id
      data      rcs_id /
     ."$Id: cray403.f,v 1.8 2014/03/26 17:10:43 stjohn Exp $"/
c
c     source file cray403.f contains routines, based on those in the the
c     "Numerical Recipes" books, used to do integrals related to neutron
c     rates
c
c     this subroutine expands the order of the solution vector u
c     from nkt to nk
c
****  dimension  u(kk,*), usave(kk,*), itran(*)
      dimension itran(*)
c
c      include 'param.i'
c      include 'soln.i'
c      include 'numbrs.i'
c
      kt = nkt + 1
      do 40 k=nk,1,-1
        if (itran(k) .le. 0)  go to 30
        kt = kt - 1
        do j=1,nj
          u(k,j) = u(kt,j)
        end do
        go to 40
   30   do j=1,nj
          u(k,j) = usave(k,j)
          if (k .le. nprim .and. itran(k) .le. 0)  u(k,j) = en(j,k)
        end do
   40 continue
      return
c
      end

      subroutine gauleg (x1, x2, x, w, n)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      parameter (eps = 3.0e-14)
      dimension  w(n), x(n)
c
      m  = (n+1)/2
      xm = 0.5*(x2+x1)
      xl = 0.5*(x2-x1)
      do i=1,m
        z  = COS (3.141592654*(i-.25)/(n+.5))
   10   p1 = 1.0
        p2 = 0.0
        do j=1,n
          p3 = p2
          p2 = p1
          p1 = ((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j
        end do
        pp = n*(z*p1-p2)/(z*z-1.0)
        z1 = z
        z  = z1-p1/pp
        if (ABS (z-z1).gt.eps)  go to 10
        x(i)     = xm-xl*z
        x(n+1-i) = xm+xl*z
        w(i)     = 2.0*xl/((1.0-z*z)*pp*pp)
        w(n+1-i) = w(i)
      end do
      return
c
      end

      subroutine polint (xa, ya, n, x, y, dy)
c
c --- Numerical recipes
c     given arrays xa and values of the function ya  return the
c     value Y at X . Extrapolation is possible. N-1 is the degree
c     of the interpolating polynomial. Dy is an estimate of the
c     error incurred.
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      parameter (nmax = 10)
      dimension  xa(n), ya(n), c(nmax), d(nmax)
c
      if (n .gt. nmax ) then
        call STOP ('subroutine POLINT: extrapolation error', 275)
      end if
c
      ns  = 1
      dif = ABS (x-xa(1))
      do i=1,n
        dift = ABS (x-xa(i))
        if (dift .lt. dif) then
          ns  = i
          dif = dift
        end if
        c(i) = ya(i)
        d(i) = ya(i)
      end do
      y  = ya(ns)
      ns = ns - 1
      do m=1,n-1
        do i=1,n-m
          ho  = xa(i)-x
          hp  = xa(i+m)-x
          w   = c(i+1)-d(i)
          den = ho-hp
          if (den .eq. 0.0)  call STOP ('subroutine POLINT: DEN=0', 173)
          den  = w  / den
          d(i) = hp * den
          c(i) = ho * den
        end do
        if (2*ns .lt. n-m) then
          dy = c(ns+1)
        else
          dy = d(ns)
          ns = ns-1
        end if
        y = y + dy
      end do
      return
c
      end

      subroutine qgaus1 (func, a, b, ss)
c------------------------------------------------------
c -- fixed abscissa 10 point Guass quadrature
c -- is used for qguas* routines, see Piessens,et all,
c -- Quadpack, Springer,1980.
c----------------------------------------------------HSJ
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      external           func
      dimension x(5),w(5)
      data x
     .     /.1488743389,.4333953941,.6794095682,.8650633666,.9739065285/
      data w
     .     /.2955242247,.2692667193,.2190863625,.1494513491,.0666713443/
c
      xm = 0.5*(b+a)
      xr = 0.5*(b-a)
      ss = 0
      do j=1,5
        dx = xr*x(j)
        ss = ss+w(j)*(func(xm+dx)+func(xm-dx))
      end do
      ss = xr * ss
      return
c
      end

      subroutine qgaus2 (func, a, b, ss)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      external           func
      dimension x(5),w(5)
      data x
     .     /.1488743389,.4333953941,.6794095682,.8650633666,.9739065285/
      data w
     .     /.2955242247,.2692667193,.2190863625,.1494513491,.0666713443/
c
      xm = 0.5*(b+a)
      xr = 0.5*(b-a)
      ss = 0
      do j=1,5
        dx = xr*x(j)
        ss = ss+w(j)*(func(xm+dx)+func(xm-dx))
      end do
      ss = xr * ss
      return
c
      end

      subroutine qgaus3 (func, a, b, ss)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      external           func
      dimension x(5), w(5)
      data x
     .     /.1488743389,.4333953941,.6794095682,.8650633666,.9739065285/
      data w
     .     /.2955242247,.2692667193,.2190863625,.1494513491,.0666713443/
c
      xm = 0.5*(b+a)
      xr = 0.5*(b-a)
      ss = 0
      do j=1,5
        dx = xr*x(j)
        ss = ss+w(j)*(func(xm+dx)+func(xm-dx))
      end do
      ss = xr * ss
      return
c
      end

      subroutine qgaus4 (func, a, b, ss)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      external           func
      dimension x(5), w(5)
      data x
     .     /.1488743389,.4333953941,.6794095682,.8650633666,.9739065285/
      data w
     .     /.2955242247,.2692667193,.2190863625,.1494513491,.0666713443/
c
      xm = 0.5*(b+a)
      xr = 0.5*(b-a)
      ss = 0
      do j=1,5
        dx = xr*x(j)
        ss = ss+w(j)*(func(xm+dx)+func(xm-dx))
      end do
      ss = xr * ss
      return
c
      end

      subroutine qgaus5 (func, a, b, ss)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      external           func
      dimension x(5), w(5)
      data x
     .     /.1488743389,.4333953941,.6794095682,.8650633666,.9739065285/
      data w
     .     /.2955242247,.2692667193,.2190863625,.1494513491,.0666713443/
c
      xm = 0.5*(b+a)
      xr = 0.5*(b-a)
      ss = 0
      do j=1,5
        dx = xr*x(j)
        ss = ss+w(j)*(func(xm+dx)+func(xm-dx))
      end do
      ss = xr * ss
      return
c
      end

      subroutine qromb (func, a, b, ss)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
****  parameter (eps = 1.0e-6, jmax = 20, jmaxp = jmax+1, k = 5, km = 4)
      parameter (eps = 1.0e-4, jmax = 20, jmaxp = jmax+1, k = 5, km = 4)
      dimension  s(jmaxp), h(jmaxp)
      real*8     zero
      external func
c
      zero = 0.0
      h(1) = 1.0
      do j=1,jmax
        call trapzd (func, a, b, s(j), j)
        if (j .ge. k) then
          l = j - km
          call polint (h(l), s(l), k, zero, ss, dss)
          if (ABS (dss) .lt. eps * ABS (ss))  return
        end if
        s(j+1) = s(j)
        h(j+1) = 0.25 * h(j)
      end do
      call STOP ('subroutine QROMB: too many steps', 174)
c
      end

      subroutine trapzd (func, a, b, s, n)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      external           func
c
      if (n .eq. 1) then
        s  = 0.5*(b-a)*(func(a)+func(b))
        it = 1
      else
        tnm = it
        del = (b-a)/tnm
        x   = a+0.5*del
        sum = 0.
        do j=1,it
          sum = sum+func(x)
          x   = x+del
        end do
        s  = 0.5*(s+(b-a)*sum/tnm)
        it = 2*it
      end if
      return
c
      end
