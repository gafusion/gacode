
      SUBROUTINE dbcevl1 (x, nx, y, ny, c, ic, xl, yl, pds, ier, icalc)
!

!
! --- this is a special version of dbcevl, used by the ONETWO suite
! --- of subroutines. It differs from the IMSL version in the calling
! --- arguments and the calculations below are reordered for efficiency
! --- That is, many times we don't need all of the pds array. This
! --- version returns only those values of pds specifically requested
! --- using the argument icalc. This leads to about at 50% improvement
! --- in execution time. --- HSJ
!
      USE nrtype,                                    ONLY : DP,I4B

      IMPLICIT  NONE

      INTEGER(I4B)  nx,ny,ic,ier,i,j,k,km1,kp1,kp2,lxpl,lx,   &
           ly,l,lxp1,icalc

      REAL(DP)  x(*),y(*),c(2,ic,*),xl,yl,pds(6), &
               hx,hy,sux(2),suy(2),su(2),svx(2),sv(2),sxy(2), &
               u,v,spln0,spln1,spln2,s0,sh,sp0,sph,h,d
!
      spln0(s0,sh,sp0,sph,h,d) = s0+d*(h*sp0+d*(3.0*(sh-s0)- &
       (sph+2.0*sp0)*h+d*(2.0*(s0-sh)+(sph+sp0)*h)))
!
      spln1(s0,sh,sp0,sph,h,d) = sp0+d*(6.0 * (sh-s0)/h-2.0* &
       (sph+2.0*sp0)+3.0*d*(2.0*(s0-sh)/h+(sph+sp0)))
!
      spln2(s0,sh,sp0,sph,h,d) = 6.0 * (sh-s0)/h**2-2.0* &
       (sph+2.0*sp0)/h+d*(2.0*(s0-sh)/h**2+(sph+sp0)/h)*6.
!
      DATA  lx, ly /0, 0/
!
      ier = 0
!
! --- correlated table search for xl
!
      CALL tableintrp (x, nx, xl, lx)
      IF (lx  .EQ. 0 )  ier = 33
      IF (lx  .EQ. nx)  ier = 35
      IF (ier .NE. 0 )  go to 100
!
! --- correlated table search for yl
!
      CALL tableintrp (y, ny, yl, ly)
      IF ( ly .EQ. 0 )  ier = 34
      IF ( ly .EQ. ny)  ier = 36
      IF (ier .NE. 0 )  go to 100

!
      lxp1 = lx+1
      hx   = x(lxp1)-x(lx)
      hy   = y(ly+1)-y(ly)
      u    = (xl-x(lx))/hx
      v    = (yl-y(ly))/hy
      k    = 2*ly
      kp1  = k+1
      kp2  = k+2
      km1  = k-1
!
      DO 25 l=1,2
         lxpl = lx-1+l
         i    = 2*(ly-1+l)
         j    = i-1
         sv(l) = spln0(c(1,lxpl,km1),c(1,lxpl,kp1),c(1,lxpl,k), &
         c(1,lxpl,kp2),hy,v)
         svx(l) = spln0(c(2,lxpl,km1),c(2,lxpl,kp1),c(2,lxpl,k), &
         c(2,lxpl,kp2),hy,v)
      IF (icalc .LT. 3)  go to 25                    ! rearranged by HSJ
         su(l) = spln0(c(1,lx,j),c(1,lxp1,j),c(2,lx,j), &
         c(2,lxp1,j),hx,u)
         suy(l) = spln0(c(1,lx,i),c(1,lxp1,i),c(2,lx,i), &
         c(2,lxp1,i),hx,u)
      IF (icalc .LT. 4)  go to 25
         sux(l) = spln1(c(1,lx,j),c(1,lxp1,j),c(2,lx,j), &
         c(2,lxp1,j),hx,u)
         sxy(l) = spln1(c(1,lx,i),c(1,lxp1,i),c(2,lx,i), &
         c(2,lxp1,i),hx,u)
   25 CONTINUE
!
      pds(1) = spln0(sv(1),sv(2),svx(1),svx(2),hx,u)
      IF (icalc .EQ. 1)  RETURN ! special routine only calculates pds(1)
      pds(2) = spln1(sv(1),sv(2),svx(1),svx(2),hx,u)
      IF (icalc .EQ. 2)  RETURN
      pds(3) = spln1(su(1),su(2),suy(1),suy(2),hy,v)
      IF (icalc .EQ. 3)  RETURN
      pds(4) = spln1(sux(1),sux(2),sxy(1),sxy(2),hy,v)
      IF (icalc .EQ. 4)  RETURN
      pds(5) = spln2(sv(1),sv(2),svx(1),svx(2),hx,u)
      IF (icalc .EQ. 5)  RETURN
      pds(6) = spln2(su(1),su(2),suy(1),suy(2),hy,v)
  100 IF (ier .GT. 0) THEN
!***    call uertst1 (ier, 'dbcevl')
      END IF
      RETURN
!
      END
