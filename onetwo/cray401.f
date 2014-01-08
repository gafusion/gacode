      subroutine abcg (j, a, b, c, g, em)
c
      USE param
      USE ions
      USE solcon
      USE soln
      USE numbrs
      USE mesh
      USE adaptive
      USE sourc
      USE geom
      USE flags  
      USE tordlrot
      USE tcoef                          ! vpinch
      USE bd_condtn
      USE solcon_gcnmp,                 ONLY : use_stab_flux
      USE flx
      USE paleocl,                      ONLY : include_paleo   ! ,paleo_flux
      USE tmpcom
      USE fdyimp,                       ONLY : afdayi,bfdayi,cfdayi
      USE glf23 !jmp.snu
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      character rcs_id*63
      save      rcs_id
      data      rcs_id /
     ."$Id: cray401.f,v 1.140 2013/12/12 18:37:01 stjohn Exp $"/
c
c ----------------------------------------------------------------------
c     this subroutine generates the matrices a, b, and c (and em and w)
c     and the vector g (and vminus and vplus) for mesh point j.
c     the matrices and vectors (all of order nk) are calculated only
c     as needed for a particular mesh point in order to save storage.
c     Note that en,te,ti,rbp,angrot (the quantities at the
c     central,n+theta, time point) are used.
c ----------------------------------------------------------------------
c


      
c
      dimension a(kk,kk), b(kk,kk), c(kk,kk), g(kk)
      dimension em(kk,kk), w(kk,kk)
      !to save between calls:
      common /abcg_com / vminus(kk), vplus(kk)
c     note  kk = kion + 4 , used to be kion+3 before rotation was added
c
c --- fdyimp is shared only with subroutine impsrc. Impsrc moves the
c --- term (1/hp)(d/drho)(-d*bp0) to the right hand (source) side.
c
      real*8          kevperg
      real*8,save :: vpluspr
      data            kevperg /6.242e8/, xmassp /1.673e-24/
      data            five_halfs_te,five_halfs_ti /2.5,2.5/


      if(no_te_convection .eq. -1)five_halfs_te =1.5
      if(no_te_convection .eq. 1) five_halfs_te =0.0
      if(no_ti_convection .eq. -1)five_halfs_ti =1.5
      if(no_ti_convection .eq. 1)five_halfs_ti =0.0
      
c
c     calculate time weighting factorc
c     calculate inverse time step
c
      dtinv = 1.0d0/ dtt
      if(diffeq_methd .eq. 1)dtinv =0.0d0 !eliminates conribution
                                          !from em to B and g for method lines
                                          ! dtinv is local to this routine
c
      zero  = 0.0
      onemt = (1.0-theta)/theta
      thetas = theta
      if(ABS(steady_state-1.0) .gt. 1.e-5)then
         if(diffeq_methd .eq. 1)
     .    call STOP('sub abcg,diffeq_methd = 1 not allowed here',0)
         onemt = 0.0
         theta =1.0
         dtinv = 1.0   
      endif

c
c     common multipliers for toroidal rotation
c
      fmult   = angrcple*angrot(j)*kevperg*0.5*iangrot
      r2omega = angrot(j)*r2capi(j)
c
c     calculate hcap ratios
c
      if (j .eq. 1 )  go to 4
      hjm1 = 0.5*(hcap(j-1)+hcap(j))/hcap(j)
      if (j .eq. nj)  go to 6
    4 hjp1 = 0.5*(hcap(j+1)+hcap(j))/hcap(j)
c
c     calculate matrix em (i.e., matrix m in GAA report)
c
    6 do 10 k=1,nk
      do 10 l=1,nk
   10 em(k,l) = 0.0
      if(ABS(steady_state -1.0d0) .lt. 1.e-5)then !if this is time dep case 
                                                  ! then setup em,etc.
        sume = ene(j)
        sumi = 0.0
        do 20 k=1,nk-3-iangrot   !loop over primary and impurity densities
        em(k,k) = 1.0
        em(nk-2-iangrot,k) = 1.5*z(j,k)*te(j)
        em(nk-1-iangrot,k) = 1.5*ti(j)+fmult*r2omega*atw(k)*xmassp
        sume = sume+en(j,k)*te(j)*dzdte(j,k)
   20   sumi = sumi + en(j,k)
        em(nk-2-iangrot,nk-2-iangrot) = 1.5*sume
        em(nk-1-iangrot,nk-1-iangrot) = 1.5*sumi
        if (iangrot .eq. 1) then
          smassden = 0.0
          do 200 k=1,nk-4
            xmassi      = atw(k)*xmassp
            smassden    = smassden+xmassi*en(j,k)
  200     em(nk,k)      = r2omega*xmassi
          em(nk-2,nk)   = smassden*r2omega*kevperg
          em(nk,nk)     = smassden*r2capi(j)
        end if
      endif



c
c --- calculate vminus and vplus
c --- vplus is the convection matrix v evaluated at j+1/2
c --- for j = 1 to nj-1. vminus is v at mesh point j-1/2 for
c --- j = 2,nj (vplus at nj-1/2 and vminus at j=1+1/2 are
c --- not used). vpluspr and vminuspr are similar but
c --- off-diagonal elements required for toroidal rotation.
c --- Note that this relies on vplus being preserved between
c --- calls so that vminus can be set equal to it. HSJ
c
      if (j .eq. 1)  go to 22
      do 21 k=1,nk
   21 vminus(k) = vplus(k)   !vplus from previous call to this routine
      vminuspr = vpluspr
   22 if (j .eq. nj)  go to 24
      fmassden = 0.0
      do 23 k=1,nk-3-iangrot
      fmassden = fmassden+atw(k)*xmassp*flux(k,j)
   23 vplus(k) = vpinch(j)
      fmassden = 0.5*fmassden*r2omega
      vplus(nk-2-iangrot) = five_halfs_te*fluxe(j)
      if(include_paleo .eq. 1)vplus(nk-2-iangrot) = 
     .             vplus(nk-2-iangrot)  !+  paleo_flux(j) 
      vplus(nk-1-iangrot) = five_halfs_ti*fluxi(j)
      vpluspr = iangrot*angrcple*kevperg*(fmassden+fluxangv(j))
      if (iangrot .eq. 1) then
c
c --- diagonal convective term in angular momentum eq.:
c
      fmassden = 0.0
      do 26 k=1,nk-4
   26 fmassden = fmassden+atw(k)*flux(k,j)
      vplus(nk) = fmassden*xmassp*r2capi(j)
      if (ichiv_model .eq. -1) vplus(nk) = 0.0 !jmp.snu 
      end if


c
c ----------------------------------------------------------------------
c  define w matrix to treat qdimpl (implicit q-delta term)
c  implicitly in the matrix equation.
c ----------------------------------------------------------------------
c
   24 do 25 k=1,nk
      do 25 l=1,nk
   25 w(k,l) = 0.0
      w(nk-iangrot-2,nk-iangrot-2) =  qdimpl(j)
      w(nk-iangrot-2,nk-iangrot-1) = -qdimpl(j)
      w(nk-iangrot-1,nk-iangrot-2) = -qdimpl(j)
      w(nk-iangrot-1,nk-iangrot-1) =  qdimpl(j)
      if (j .eq. 1)  go to 30
      em(nk-iangrot,nk-iangrot) = 1.0 /
     .                (fcap(j) * gcap(j) * (hcap(j) * r(j))**2 * twkfar)
c       print *,'em here,nk =',em(nk-iangrot,nk-iangrot),nk
        if(ABS(steady_state -1.0d0) .gt. 1.e-5)
     .  em(nk-iangrot,nk-iangrot) =em(nk-iangrot,nk-iangrot)*
     .                                      steady_state
      go to 60
c
c  calculate matrices b and c and vector g for center mesh point (j = 1)
c  matrix a is not defined for j = 1 (since it multiplies u(K,j-1).
   30 do k=1,nk
      if (k .ne. nk-iangrot) then
c
        do 35 l=1,nk
          b(k,l) =  theta*drr(j)*rrp(j)*hjp1*d(k,l,j)
     .            + dtinv*em(k,l) + theta*w(k,l)
   35   c(k,l) = -theta*drr(j)*rrp(j)*hjp1*d(k,l,j)
        b(k,k) = b(k,k) + theta*drr(j)*rrp(j)*hjp1*dr(j)
     .                      * MAX (vplus(k), zero)
        c(k,k) = c(k,k) + theta*drr(j)*rrp(j)*hjp1*dr(j)
     .                      * MIN (vplus(k), zero)
        g(k) = s(k,j)  !here s does not include qdelt !!

c
      else
c
c --- Faraday's law:
c
        do 50 l=1,nk
          b(k,l) = 0.0
   50   c(k,l) = 0.0
          b(k,k) = 1.0
          g(k) = 0.0
      end if
      if (((angrcple .ne. 0.0) .and. (k .eq. nk-2)) .and.
     .                            (iangrot .eq. 1)) then
c
c --- add the off-diagonal ion energy term due to
c --- radial toroidal momentum transport.
c --- assume the diagonal term dominates the convective wind
c
        b(k,nk) = b(k,nk) + (theta * drr(j) * rrp(j) * hjp1
     .          * dr(j) * MAX (vpluspr, zero)) * angrcple
        c(k,nk) = c(k,nk) + (theta * drr(j) * rrp(j) * hjp1
     .          * dr(j) * MIN (vpluspr, zero)) * angrcple
      end if
      if (k .ne. nk-iangrot) then
        do l=1,nk
          g(k) = g(k)
     .          -onemt* (b(k,l)-dtinv*em(k,l))*usave(l,j)
     .          +dtinv*em(k,l)*usave(l,j)*steady_state
     .          -c(k,l)*usave(l,j+1)*onemt
        end do

      end if
      end do
      theta =thetas

      return
c
c  calculate matrices a, b, and c and vector g for interior mesh
c     points (j = 2,...,nj-1)
c
   60 if (j .eq. nj)  go to 80
c



      do   k=1,nk
        do l=1,nk
          a(k,l) = -theta*drr(j)*rrm(j)*hjm1*d(k,l,j-1)
          c(k,l) = -theta*drr(j)*rrp(j)*hjp1*d(k,l,j)
c          b(k,l) =  theta*drr(j)*rrm(j)*hjm1*d(k,l,j-1)
c     .            + theta*drr(j)*rrp(j)*hjp1*d(k,l,j)
c     .            + dtinv*em(k,l) + theta*w(k,l)
          b(k,l) = dtinv*em(k,l) + theta*w(k,l) -(c(k,l)+a(k,l))
        end do
      end do
c
c ----------------------------------------------------------------------
c  compute implicit derivative term for Faraday's law
c ----------------------------------------------------------------------
c
      if (codeid .eq. 'onedee')  go to 68
      denom = (r(j-1)-r(j))*(r(j+1)-r(j))*(r(j-1)-r(j+1))
     .                     *hcap(j)*r(j)
      cnum = (r(j-1)-r(j))**2
      anum = -(r(j+1)-r(j))**2
      bnum = -cnum-anum
      if (include_adaptive .eq. 0) then
               dscripla = -(adot+bdot*r(j+1))/
     .                        (ascrip+2.0*bscrip*r(j+1))*r(j+1)
               dscriplb = -(adot+bdot*r(j))/(ascrip+2.0*
     .                                        bscrip*r(j))*r(j)
               dscriplc = -(adot+bdot*r(j-1))/(ascrip+2.0*
     .                                           bscrip*r(j-1))
      else
         dscripla = drhodt_adaptive(j+1)
         dscriplb = drhodt_adaptive(j)
         dscriplc = drhodt_adaptive(j-1)
      end if
      cnum = -0.5*cnum/denom*dscripla
     .     /(r(j+1)*fcap(j+1)*gcap(j+1)*hcap(j+1))
c
****  dscripl = -(adot+bdot*r(j))/(ascrip+2.0*bscrip*r(j))*r(j)
      bnum = -0.5*bnum/denom*dscriplb
     .     /(r(j)*fcap(j)*gcap(j)*hcap(j))
c      dscripl = -(adot+bdot*r(j-1))/(ascrip+2.0*bscrip*r(j-1))
      anum = -0.5*anum/denom*dscriplc
     .     /(fcap(j-1)*gcap(j-1)*hcap(j-1))
      afdayi(j) = 2.0 * theta*anum
      bfdayi(j) = 2.0 * theta*bnum
      cfdayi(j) = 2.0 * theta*cnum
      a(nk-iangrot,nk-iangrot) = a(nk-iangrot,nk-iangrot)+afdayi(j)
      b(nk-iangrot,nk-iangrot) = b(nk-iangrot,nk-iangrot)+bfdayi(j)
      c(nk-iangrot,nk-iangrot) = c(nk-iangrot,nk-iangrot)+cfdayi(j)
c

c ----------------------------------------------------------------------
c

 68   continue
      do 70 k=1,nk
      a(k,k) = a(k,k) - theta*drr(j)*rrm(j)*hjm1*dr(j-1)
     .                        * MAX (vminus(k), zero)
 
      b(k,k) = b(k,k) - theta*drr(j)*rrm(j)*hjm1*dr(j-1)
     .                        * MIN (vminus(k), zero)
      b(k,k) = b(k,k) + theta*drr(j)*rrp(j)*hjp1*dr(j)
     .                        * MAX (vplus(k), zero)
 
      c(k,k) = c(k,k) + theta*drr(j)*rrp(j)*hjp1*dr(j)
     .                        * MIN (vplus(k), zero)
      xke =0.0D0 ;xki =0.0D0 ; xkw =0.0D0
      if(use_stab_flux) call source_mod_stab(-1,j,xke,xki,xkw)
   70 g(k) = s(k,j) - comp_term(k)
      if ( iangrot .eq. 0  ) go to 75
      if (angrcple .ne. 0.0) then
c
c --- add the off diagonal energy convection term due to radial
c --- toroidal momentum flux.  Assume the diagonal term dominates
c --- the convective wind:
c
        a(nk-2,nk) = a(nk-2,nk)-(theta*drr(j)*rrm(j)*hjm1
     .     *dr(j-1) * MAX (vminuspr, zero))*angrcple
        b(nk-2,nk) = b(nk-2,nk)+(-theta*drr(j)*rrm(j)*hjm1*
     .     dr(j-1) * MIN (vminuspr, zero)
     .    + theta*drr(j)*rrp(j)*hjp1*dr(j)
     .     * MAX (vpluspr, zero))*angrcple
         if(ABS(b(nk-2,nk)) .gt. 1.d100)then
             print *,'abcg 4 ',b(k,l)
             call exit(1)
          endif
        c(nk-2,nk) = c(nk-2,nk)+(theta*drr(j)*rrp(j)*hjp1*dr(j)
     .     * MIN (vpluspr, zero))*angrcple
      end if
c
   75 do 78 k=1,nk
      do 78 l=1,nk
   78 g(k) = g(k)-onemt*a(k,l)*usave(l,j-1)
     .      -onemt* (b(k,l)-dtinv*em(k,l))*usave(l,j)
     .      +dtinv*em(k,l)*usave(l,j)*steady_state
     .      -c(k,l)*usave(l,j+1)*onemt
      theta =thetas
 
      return
c
c  calculate matrices a and b and vector g for boundary mesh
c     point (j = nj)
c
   80 do k=1,nk
        do 90 l=1,nk
        a(k,l) = 0.0
   90   b(k,l) = 0.0
        b(k,k) = 1.0
        g(k) = ub(k)
      end do
      theta =thetas
      return
c
      end




      subroutine addac (b, a, nj, c)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c --- ADDAC adds the array c*b to the array a
c
      dimension  a(*), b(*)
c
      do j=1,nj
        a(j) = a(j) + c*b(j)
      end do
      return
c
      end





      subroutine adda (b, a, nj)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c --- ADDA adds the array b to the array a
c
      dimension  a(*), b(*)
c
      do j=1,nj
        a(j) = a(j) + b(j)
      end do
      return
c
      end

      real*8 function asimp (a1, b, ep, m, n, FUN)
c
      USE io,only : ncrt
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c --- author:  k. hillstrom (argonne national laboratory, chicago, illinois)
c
      external FUN
c
      real*8    a1, b, ep, FUN, a, eps, absar, est, fa, fm, fb, dx, sx,
     .          f1, f2, fbp, est2, nrtr, est1, sum, daft, esum, tsum, da
      real*8    aest2, ftst, fmax, aest1, delta, aest
      dimension f2(30), fbp(30), est2(30), nrtr(30), aest2(30), ftst(3)
c
c     the parameter setup for the initial call
c
      n7 = ncrt
      if (n .le. 0) then
        write (n7, '(a)')  ' **** ASIMP error return:  n .le. 0 ****'
        asimp = 0.0
        return
      end if
c
      if (n .gt. 3) then
        write (n7, '(a)')  ' **** ASIMP error return:  n .gt. 3 ****'
        asimp = 0.0
        return
      end if
c
      a       = a1
      eps     = ep*15.0
      esum    = 0.0
      tsum    = 0.0
      lvl     = 1
      da      = b-a
      fa      = FUN(a)
      fm      = FUN((a+b)*0.5)
      fb      = FUN(b)
      m       = 3
      fmax    = ABS (fa)
      ftst(1) = fmax
      ftst(2) = ABS (fm)
      ftst(3) = ABS (fb)
      do 10 i=2,3
        if (fmax .ge. ftst(i))  go to 10
        fmax = ftst(i)
   10 continue
      est   = (fa+4.0*fm+fb)*da/6.0
      absar = (ftst(1)+4.0*ftst(2)+ftst(3))*da/6.0
      aest  = absar
c
c 1 = recur
c
   20 dx         = da/(2.0**lvl)
      sx         = dx/6.0
      f1         = FUN(a+0.5*dx)
      f2(lvl)    = FUN(a+1.5*dx)
      est1       = sx*(fa+4.0*f1+fm)
      fbp(lvl)   = fb
      est2(lvl)  = sx*(fm+4.0*f2(lvl)+fb)
      sum        = est1+est2(lvl)
      ftst(1)    = ABS (f1)
      ftst(2)    = ABS (f2(lvl))
      ftst(3)    = ABS (fm)
      aest1      = sx*(ABS (fa)+4.0*ftst(1)+ftst(3))
      aest2(lvl) = sx*(ftst(3) +4.0*ftst(2) + ABS (fb))
      absar      = absar-aest+aest1+aest2(lvl)
      m          = m+2
      go to (60,30,70),n
   30 delta = absar
      go to 90
   60 delta = 1.0
      go to 90
   70 do 80 i=1,2
        if (fmax .ge. ftst(i))  go to 80
        fmax = ftst(i)
   80 continue
      delta = fmax
   90 diff  = ABS (est-sum)
      daft  = (est-sum)/15.0
      if (diff-eps*delta) 110,110,100
  100 if (lvl-30) 140,120,120
  110 if (lvl- 1) 120,140,120
c
c 2 = up
c
  120 a    = a+2.0*dx
  130 lvl  = lvl-1
      esum = esum+daft
      l    = nrtr(lvl)
      tsum = tsum+sum
      go to (160, 170), l
c
c 11 = r1,12=r2
c
  140 nrtr(lvl) = 1
      est = est1
      aest = aest1
      fb = fm
      fm = f1
      eps = eps/2.0
  150 lvl = lvl+1
      go to 20
  160 nrtr(lvl) = 2
      fa = fb
      fm = f2(lvl)
      fb = fbp(lvl)
      est = est2(lvl)
      aest = aest2(lvl)
      go to 150
  170 eps = 2.0 * eps
      sum = 0.0
      if (lvl-1) 180,180,130
  180 asimp = tsum
      a = ABS (esum)
      ep = diff/delta
      if (a .ge. ep)  go to 190
      asimp = asimp - esum
  190 return
c
      end

      subroutine average
c---------------------------------------------------------------
c      average the calculated solution in u over 2*ddebug(12) + 1
c      grid points, with ddebug(2) points to the left and right
c      of the current grid point. Poloidal B field is not averaged.
c      points near the magnetic axis are not averaged since zero
c      gradient is destroyed by this averaging
c      points near the boundary are not averaged due to boundary conditions
c-------------------------------------------------------------------

      USE param
      USE io
      USE soln
      USE numbrs
      USE mesh
      USE flags
      implicit  integer (i-n), real*8 (a-h, o-z)
c      include 'param.i'  
c      include 'flags.i'  !itran
c      include 'io.i'     !ddebug(12)
c      include 'mesh.i'   !te,ti,rot_index  parameters
c      include 'numbrs.i' !nk,nj,nion
c      include 'soln.i'   !u
      
      ngrid = ddebug(12)
      do k= 1, nk
          if(itran(k) .eq. 1)then
               if(k .eq. nion + 3) go to 10  !rbp
               if(k .le. nion) nedge = nj
               if(k .eq. nion+1) nedge = te_index
               if(k .eq. nion+2) nedge=  ti_index
               if(k.eq.  nion+4) nedge =  rot_index
               js = ngrid + 1 
               if(nedge .eq. nj)then
                  je = nj - ngrid
               else
                  je = nedge-1  !points j ge nedge are not aeraged
               endif
               do j= js,je   ! grid points j ge nedge are  fixed
                  sum =0.0
                  jup = MIN(j+ngrid,nj)
                  i=0
                  do l= j-ngrid,jup
                     i = i+1
                     sum=sum+u(k,l)
                  enddo
                  u(k,j) = sum/float(i)
               enddo
          endif
 10   enddo 
      return
      end

      subroutine bilin (nx, ny, x, y, u, xx, yy, uu)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c     This subroutine performs a bilinear interpolation to determine
c     uu = u(xx,yy) given u on an x-y mesh.
c ----------------------------------------------------------------------
c
      dimension u(nx,*), x(*), y(*)
c
      do i=2,nx
        if (x(i) .ge. xx)  go to 60
      end do
      i  = nx
   60 i1 = i-1
      i2 = i
      do j=2,ny
        if (y(j) .ge. yy)  go to 80
      end do
      j  = ny
   80 j1 = j-1
      j2 = j
      a1 = (x(i2)-xx)*(y(j2)-yy)
      a2 = (x(i2)-xx)*(yy-y(j1))
      a3 = (xx-x(i1))*(y(j2)-yy)
      a4 = (xx-x(i1))*(yy-y(j1))
      uu = (a1*u(i1,j1)+a2*u(i1,j2)+a3*u(i2,j1)+a4*u(i2,j2))
     .     /((x(i2)-x(i1))*(y(j2)-y(j1)))
      return
c
      end

      subroutine capcheck (fcap0, hcap0, fcap, hcap, fcapctr,
     .                     hcapctr, nj, theta, delcapc)
c
      USE soln2d
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c   subroutine checks the implicit cap parameters for convergence
c
c  input
c
c    fcap0               the cap values at the start of the predictor-
c                         corrector cycle
c    hcap0
c
c    fcap                the cap values at the end of the current
c                        predictor-corrector cycle.
c    hcap
c
c    fcapctr               the cap values at the central time point
c                          of the predictor-corrector cycle.
c    hcapctr
c
c    theta               the weight factor that determines
c                        the location of the central time point.
c    nj                  the rho grid size.
c
c  output:
c
c    delcapc             the maximum relative change in the cap
c                        parameters
c ------------------------------------------------------------------ HSJ
c
c      include 'soln2d.i'
c
      dimension   fcap0  (*), hcap0  (*),
     .            fcap   (*), hcap   (*),
     .            fcapctr(*), hcapctr(*)
      character*4 a
c
      if (.not. implicit_fh)  return
      delcapc = 0.0
      onemt   = 1.0 - theta
      df      = 0.0
      dh      = 0.0
c
      do j=1,nj
        fnewc      = onemt*fcap0(j)+theta*fcap(j)
        hnewc      = onemt*hcap0(j)+theta*hcap(j)
        df         = ABS ((fnewc-fcapctr(j))/fcapctr(j))
        dh         = ABS ((hnewc-hcapctr(j))/hcapctr(j))
        hcapctr(j) = hnewc
        fcapctr(j) = fnewc
        delcapo    = delcapc
        delcapc    = MAX (delcapc,df,dh)
        if (delcapo .ne. delcapc) then
          if (delcapc .eq. df)  a = 'fcap'
          if (delcapc .eq. dh)  a = 'hcap'
          jmaxsave = j
        end if
      end do
c
****  write  (6, 1)  a, jmaxsave, delcapc
****1 format (' a, jmax, delmax for cap parameters = ',
****.           a, 2x, i3, 2x, 1pe12.2)
      return
c
      end

      subroutine check_solution (u, usave, en, te, ti, rbp, nk, nj, kj,
     .                           kk, r, delit, theta, itran, iangrot,
     .                           angrot, relaxsol, ksave, jsave, value,
     .                           te_index,ti_index,rot_index)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c     this subroutine calculates delit, the maximum relative difference
c     between the predicted and corrected values at the half time point.
c     in addition, the predicted values are updated (by delrep)
c     points adjacent to a rbp-sign-change are excluded.
c     relaxsol is used to (under-)relax the updating of the predictor
c     (the finite difference equations are written in such a way
c     that the value at the first grid point in rho is determined
c     so that a zero gradient boundary condition is applicable.
c     Consequently the value for j = 1 is not relaxed). ------------ HSJ
c ----------------------------------------------------------------------
c
      integer te_index,ti_index,rot_index
      dimension u(kk,*), usave(kk,*), en(kj,*), te(*), ti(*), rbp(*),
     .          itran(*), r(*), ischng(3), angrot(*)
      data ncycle /0/
c
c     check for up to 3 sign changes in rbp
c
      onemt     = 1.0 - theta
      ischng(1) = 0
      ischng(2) = 0
      ischng(3) = 0
      isc       = 1
c
      do 5 j=2,nj-1
        if (usave(nk-iangrot,j)*usave(nk-iangrot,j+1) .gt. 0.0)  go to 5
        ischng(isc) = j
        isc         = isc + 1
        if (isc .gt. 3)  go to 6
    5 continue
c
    6 delit   = 0.0
      delsave = delit
      ncyle_set = 10
c
      ncycle =ncycle+1
      icycle  = MOD(ncycle,ncyle_set)

      do 4000 j=1,nj
      relaxc = relaxsol
      if(icycle .eq. 0)then 
          relaxc = 1.0
          ncycle = 0
      endif
      if (j .eq. 1)  relaxc = 1.0  ! i.e., no under relaxation

      do 2100 k=1,nk-3-iangrot     ! i.e., k=1,2..nion (nion=nprim+nimp)
        if (itran(k) .le. 0)  go to 2100
        call delrep (u(k,j), usave(k,j), en(j,k),
     .               delit, theta, relaxc, 0)
        if (delsave .ne. delit) then
          delsave = delit
          ksave   = k
          jsave   = j
          value   = en(j,k)
        end if
 2100 continue
c
      if (itran(nk-2-iangrot) .le. 0)  go to 2200    ! te
      if(j .ge. te_index) go to 2200
c     delrep doesnt change u or usave. te is updated to be the
c     newly estimated value at the central time point.
      call delrep (u(nk-2-iangrot,j), usave(nk-2-iangrot,j), te(j),
     .             delit, theta, relaxc, 1)
      if (delsave .ne. delit) then
        delsave = delit
        ksave   = nk - 2 - iangrot
        jsave   = j
        value   = te(j)
      end if
c
 2200 if (itran(nk-1-iangrot) .le. 0)  go to 2250    ! ti
      if(j .ge. ti_index) go to 2250
c     delrep doesnt change u or usave. ti is updated to be the
c     newly estimated value at the central time point.
      call delrep (u(nk-1-iangrot,j), usave(nk-1-iangrot,j), ti(j),
     .             delit, theta, relaxc, 1)
      if (delsave .ne. delit) then
        delsave = delit
        ksave   = nk-1-iangrot
        jsave   = j
        value   = ti(j)
      end if
 2250 if (iangrot .ne. 0) then
        if (itran(nk) .le. 0)  go to 2300    ! rotation
       if(j .ge. rot_index) go to 2300
c     delrep doesnt change u or usave. angrot is updated to be the
c     newly estimated value at the central time point.
        call delrep (u(nk,j), usave(nk,j), angrot(j),
     .               delit, theta, relaxc, 0)
        if (delsave .ne. delit) then
          delsave = delit
          ksave   = nk
          jsave   = j
          value   = angrot(j)
        end if
      end if
 2300 if (itran(nk-iangrot) .le. 0)  go to 4000    ! rbp
      if (j .eq. 1 .or.  j .eq. nj)  go to 3999
c
c     omit points adjacent to a sign change
c
      do is=1,3
        if (ischng(is) .eq. 0)  go to 20
        if (j .ge. (ischng(is)-1) .and.
     .      j .le. (ischng(is)+2))  go to 3999
      end do
c
c     don't update delit for points where rbp is less than 5% of
c     parabolic current profile value (but still update rbp to
c     the new predicted value at the central time point)
c
   20 rnorm2 = (r(j)/r(nj))**2
      if (ABS (rbp(j)) .lt. 0.05*rbp(nj)*rnorm2*(2.0-rnorm2)) go to 3999
c     delrep doesnt change u or usave. rbp is updated to be the
c     newly estimated value at the central time point.
      call delrep (u(nk-iangrot,j), usave(nk-iangrot,j),
     .             rbp(j), delit, theta, relaxc, 0)
      if (delsave .ne. delit) then
        delsave = delit
        ksave   = nk - iangrot
        jsave   = j
        value   = rbp(j)
      end if
      go to 4000
 3999 rbp(j) = theta*u(nk-iangrot,j)+onemt*usave(nk-iangrot,j)
 4000 continue
      return
c
      end

      subroutine chekdt (time, dt, n, timprt, mprt, timplt, mplot,
     .                   irf, dthold, dtevmin,
     .                   timpel, ipelet, npel, timrfp, iplot,
     .                   iprt, irfplt, codeid, time0, dteq,
     .                   eqtim0,taus,timmax)
c
      USE param
      USE nub
      USE numbrs
      USE verbose
      USE nub4
      USE io,             ONLY : prtlst,nprtlst_max,pltlst,npltlst_max
      USE transp,         ONLY : beam_data,use_nubeam,nubeam_dt,
     .                           nubeam0_dt,nubeam_calls,nubeam_steps


      USE solcon,only    : time_tol,dtmax,dtmin
      USE bd_condtn,only : bctime,ub,fluxb,
     .    ub_save,ub_rho_edge,bctime_zone
      USE events
      USE P_Nfreya_12_interface,        ONLY : use_P_Nfreya

      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------

c
      character*8 codeid 
      dimension   timrfp(*), irf(*)
      dimension   timpel(10),taus(*)
      integer*4 itest(5),ln,lt


c
c --- if time_tol above is changed,
c     change also in subroutines GETMHDBC and TSPLINEW
c
c assign some values to event array
c this routine sets timevent for those events where the
c time is not set apriori or the time is not given in a list
c (events that fall into these categories are the old single pulse beam
c  and rf times, both of which are set outside of this routine)
c For the other events this routine determines the next event
c time greater than the current time.
c The time step,dt,  is then selected that  is compatible with the events.
c this program checks for the following events:
c
c     1. timmax
c     2. plot  point (timplt) 3D plot time interval
c     3. print point (prtlst, timprt)
c     4. beam plot point (timbplt)
c     5. dteq
c     6. plot point (pltlst)
c     7. RF plot point (timrfp)
c     8. pellet injection time
c     9. beam on
c    10. beam off
c    11. rf(1) on
c    12. rf(1) off
c    13. rf(2) on
c    14. rf(2) off
c    15.
c      .     all these sources were put here to simulate
c      .     time-dependent sources  HSJ
c    70.  rf(30) off  (krf = 30)
c    71.  time-dependent boundary conditions (if time_dep_beam =1)
c    72.  time dependent beam pulse on  (if time_dep_beam =1)
c    73.  time dependent beam pulse off (if time_dep_beam =1)
c    74.  beam integral upper limit
c    75.  beam integral lower limit
c    76.  nubeam switching time (either pulse on or off or 
c         change in power level of any one of the beams.
c    77.  blank
c    78   P_Nf beam(1) on  ! 32 Parallel nfreya beamlets
c    79   p_Nf beam(1) off
c    .
c    .
c   141   P_Nf beam (32) off 

      time_out_of_bounds = MAX(timmax,time0) + 10.*dt
      one_third =1.d0/3.d0
      isave_event = 0

c ----------------------------------------------------------------------
c determine end of this equilibrium/transport iteration
c dthold is initialized to -1.  if the time step is adjusted dthold will
c be used to reset the time step to its value prior to this reduction
c for the next time step
c determine end of this equilibrium/transport iteration
c ----------------------------------------------------------------------
c
      timevent(5) = time_out_of_bounds
      if (codeid .eq. 'onedee')  go to 1900
      timevent(5) = eqtim0 + dteq
c
c ----------------------------------------------------------------------
c intialize event switches to zero
c ----------------------------------------------------------------------
c
 1900 do i=1,nevents
        ievent(i) = 0
      end do
c
      timnew = time + dt
c
c ----------------------------------------------------------------------
c determine next plot time
c ----------------------------------------------------------------------
c
      timevent(2) = time_out_of_bounds
      timevent(6) = 60.0
      do i=1,npltlst_max
         k=i
        if (pltlst(i) .gt. timevent(1))  go to 2030 ! timevent(1) = timmax
        if (pltlst(i) .gt. time       )  go to 2020
      end do
c
      go to 2030
 2020 timevent(6) = pltlst(k)


c
c     interval for 3D plots is timplt:
 2030 if (timplt .eq. 0.0)  go to 2040
      i    = (time   + time_tol - time0) / timplt
      inew = (timnew + time_tol - time0) / timplt
      k    = inew - i
      if (k .ne. 0)  timevent(2) = time0 + (i+1) * timplt
c
c ----------------------------------------------------------------------
c determine next print time
c ----------------------------------------------------------------------
c
 2040 timevent(3) = time_out_of_bounds
      do i=1, nprtlst_max
        if (prtlst(i) .gt. timevent(1))  go to 2070  ! timevent(1) = timmax
        if (prtlst(i) .gt. time       )  go to 2060
      end do
      go to 2070
 2060 timevent(3) = prtlst(i)
      go to 2080

c
c     timprt is time iterval between printouts
 2070 if (timprt .eq.  0.0)  go to 2080
****  i    =  time          / timprt
      i    = (time   + time_tol) / timprt     ! add-in of time_tol is experimental
****  inew =  timnew        / timprt
      inew = (timnew + time_tol) / timprt     ! add-in of time_tol is experimental
      k    = inew - i
      if (k .ne. 0)  timevent(3) = (i+1) * timprt
c
c ----------------------------------------------------------------------
c determine next plot time for beam plot file
c ----------------------------------------------------------------------
c
 2080 timevent(4) = time_out_of_bounds
      do i=1,ntimbplt
        if (timbplt(i) .gt. time)  go to 3000
      end do
      go to 3050
 3000 timevent(4) = timbplt(i)
c
c ----------------------------------------------------------------------
c determine next plot time for RF plot file
c ----------------------------------------------------------------------
c
 3050 timevent(7) = time_out_of_bounds
      do i=1,5  !again hard wired
        if (timrfp(i) .gt. time)  go to 3070
      end do
      go to 3075
 3070 timevent(7) = timrfp(i)

c----------------------------------------------------------------------
c time dependent beam. Note that this is different than the old single
c pulse beam ( which we continue to support for backward compatibility)
c Any Onetwo run will have either the old single pulse
c input or the new multi pulse input, mixing of the two is not allowed. 
c or the new nubeam input or the P_nfreya input     HSJ
c------------------------------------------------------------------------
 3075 if(time_dep_beam .eq. 1 .AND. .NOT. use_P_Nfreya) then
c        make sure old single pulse input is turned off:
         timevent(9) =  time_out_of_bounds
         timevent(10) =  time_out_of_bounds

c        now find that pulse which is the closest to 
c        in tie to switching on or off. Use that pulse to
c        determine if dt will have to be modified.
         ionoff = -1
         stime = timmax -time0 + 1.
         do l=1,nbeams
            ll=l
            do i = 1,n21s
               ii=i
               do np =1,n_pulse
                  nnpp =np
                  ontime = pbeamOn(np,i,l) - time
                  offtime = pbeamOff(np,i,l) - time
                  if(ontime .gt. 0.0e0)then
                    stime = MIN(ontime,stime)
                    if(ABS(stime - ontime) .lt. time_tol )then
                       ionoff = 1
                       go to 3071
                    endif
                  elseif(offtime .gt. 0.0e0)then
                    stime = MIN(offtime,stime)
                    if(ABS(stime- offtime) .lt. time_tol)then
                       ionoff = 0
                       go to 3071
                    endif
                  endif
               enddo ! end loop on pulse #
 3071          continue
            enddo    ! end loop on source #
          enddo      ! end loop on beam #
          if(stime .lt. dt )then
             if(ionoff .eq. 1)then
                timevent(index_pbon) = pbeamOn(nnpp,ii,ll)
                timevent(index_pboff) = time_out_of_bounds
             else
                timevent(index_pboff) = pbeamOff(nnpp,ii,ll)
                timevent(index_pbon) = time_out_of_bounds 
             endif
          else
                timevent(index_pbon) = time_out_of_bounds 
                timevent(index_pboff) = time_out_of_bounds     
          endif

c------check upper and lower limits of integration on beam integrals
c      done in subroutine beam_time_dependance. 
c      We have to monitor two additional times:
c
c      a)  the lower limit of integration reaching the cutoff value
c          of  vthermal or 0.0 .
c      b)  the upper limit of integration reaching  the
c          lower limit ( ie the integral no longer contributes) 
c          Which happens once the pulse is off and the lower limit is
c          already at the cutoff.
c      The information needed to to the calcs (eg vcrit,etc.) is not
c      available initially .

      if ( beam_pulse_control .gt. -1  .and. beam_time_init .eq. 1)then
         t_tol=time_tol
         dtnew =   dt              
         do l=1,nbeams                !loop over  all beamlines
           if( beam_thermal_cutoff .eq. -1)
     .     tau0_vcut = tau0_cut_avg(taus,vcrit,kj,kb,nj,l)
            do i=1,n21s               !loop over ns  sources per beamline
               if( timmax .lt. pbeamOn(1,i,l) )go to 10 ! this source is not used in current run
               nbe = ke
               if(neg_ion_source(l) .eq. 1)nbe =1
               do ie = 1,nbe          !loop over energies for each beam line
                 do j =1,nj           !loop over transport grid
                   vc3 = vcrit(j,l)**3
                   vthcut3 = eval_vthcut3(j)
                   if( beam_thermal_cutoff .eq. -1) then
                      tau0_vcutoff = tau0_vcut
                   else
                      tau0_vcutoff = one_third*taus(j)*
     .                            LOG((vthcut3+vc3)/vc3)
                   endif
                   do k=1,n_pulse    ! loop over input  time values in pbeamOn(k,i,l)
                     if(j*k .eq. 1 )
     .               write(6,'("j=1,k=1,i,ie,l,pssvoff =",3(i3,1x),l1)')
     .                       i,ie,l,pssvoff(k,j,i,ie,l)
                    if( .not. pssvoff(k,j,i,ie,l) )then 
                          !pssvoff = true when pulse k is contributing. 
                          !tau0_vljoff = tau0_vlj

                          dtslowboff = tau0_vljoff(k,j,i,ie,l) - 
     .                                                      tau0_vcutoff
                          dtslowbon = tau0_vlj(k,j,i,ie,l) - 
     .                                                      tau0_vcutoff
                          if( tau0_vljoff(k,j,i,ie,l) .lt. 0.0d0)then
                             !tau0_vljoff  is .lt. 0.0 if it has not been given a value in
                             !sub beam_time_depdendance (because the pulse is not
                             !yet shut off). In this case we must check to see if
                             !tau0_vlj has reached the thermal cutoff (this would happen
                             !if the pulse were on longer than the beam slowing down time)
                             dtslowb = dtslowbon  !upper limit is checked
                                iul_ll = -1
                          else
                             !here tau0_vljoff is set which means that pulse k was
                             !shut off and the lower limit of integration is at the
                             !thermal cutoff or is approaching it.
                             if(dtslowboff .le. 0.0d0)then !lower limit is at cutoff value
                                dtslowb = dtslowbon  !upper limit is checked
c                                iul_ll = -1
                                iul_ll = 1
                             else    ! lower limit not at cutoff
                                dtslowb = dtslowboff !lower limit is checked
c                                iul_ll = 1
                                iul_ll = -1
                             endif
                           endif
                          if( dtslowb  .gt. t_tol .and. dtslowb 
     .                                          .lt. dt )then
                             dtnew = MIN(dtslowb,dtnew)
                             if( ABS(dtnew -dtslowb) .lt. t_tol)then
                                itest(1) =k
                                itest(2) =j
                                itest(3) = i
                                itest(4) = ie
                                itest(5) = l
                                if (iul_ll .eq. 1)then
                                  timevent(index_beam_ll) = 
     .                                              time_out_of_bounds
                                  timevent(index_beam_ul) = time+dtnew
                                else
                                  timevent(index_beam_ul) = 
     .                                              time_out_of_bounds
                                  timevent(index_beam_ll) = time+dtnew
                                endif
                             endif
                          endif
                    endif              !end  pssvoff branch
                 enddo          !end  k loop over times
                enddo                  !end  j loop over transport grid
              enddo                    !end  ie loop on beam energies
 10         enddo                      !end  i loop on sources per beamline
          enddo                        !end  l loop on beamlines
          if(ABS( dtnew - dt ) .gt. t_tol)then
            print *,' dt greater than dtnew: dt,dtnew =',dt,dtnew
            if(iul_ll .eq. 1)then
              print *,'timevent(index_beam_ll ) = ',
     .                            timevent(index_beam_ll )
            else
              print *,'timevent(index_beam_ul ) = ',
     .                            timevent(index_beam_ul )
            endif
            print *,'^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'
          else
              timevent(index_beam_ul) = time_out_of_bounds
              timevent(index_beam_ll) = time_out_of_bounds
          endif
        endif      ! beam_pulse_control branch
      endif        ! time_dep_beam branch



c      do loop 5010 takes care of this:
c     IF(use_P_Nfreya)THEN
c       !timevent(index_p_nf_start:index_p_end ) contains 
c       ! on and off times of beamlets as set in sub init.
c       ! these times are NOT sorted so check all beamlets
c       x= HUGE(time)
c       DO i=index_p_nf_start,index_p_nf_end
c          time_pnf = timevent(i)
c          time_pnf = timevent(i) - time
c          IF(time_pnf .GE. 0.0D0)THEN
c             x = MIN(x,time_pnf)
c          ENDIF
c       ENDDO
c       
c     ENDIF


c ----------------------------------------------------------------------
c determine next boundary condition time for kinetic profiles
c (boundary conditions for MHD calculations are determined
c automatically through dteq which is event #5)
c ----------------------------------------------------------------------
c
 3080 timevent(index_bc) = time_out_of_bounds
      do i=1,nbctim
        timebcl = bctime(i)
        if (timebcl  .gt. time)  go to 4000
      end do
      go to 4010
 4000 timevent(index_bc) = timebcl
      print *,'timevent,index_bc,time =',timebcl,index_bc,time

c
c ----------------------------------------------------------------------
c determine if we will be skipping any events in this time step
c note that other events such as pellet time,pulse beam time
c and times for calling rf related calcs are included in the 5010 loop
c below in addition to the event times set above.
c ----------------------------------------------------------------------
c
 4010 dtxx   =  1.0e20   
      dthold = -1.0
      timnew = time + dt
      neve = nevents
c      if(use_nubeam)neve = neve-1
      if(use_nubeam)neve = index_p_nf_start - 2 ! skip p_nf events in this case
      do 5010 i=1,neve
        if (time        .ge. timevent(i) )  go to 5010
        if (timevent(i) .gt. timnew + time_tol)  go to 5010
        x = timevent(i) - time + time_tol
        x = timevent(i) - time - time_tol ! hsj 3/22/12
        if (x .lt. dtxx) then
          dtxx = x
          isave_event = i
c        write(969,FMT='("changed dtxx,timevent,i =",2x,1pe16.8,
c     .         2x,1pe16.8,2x,i5)')dtxx,timevent(i),i ! 8888889999
c        write(969,FMT='("time,time_tol =",2(2x,1pe16.8))')
c     .        time,time_tol   ! 8888889999
        end if
 5010 continue

    
      IF(use_P_Nfreya)THEN
c         Force a plot point next to beam switching times using this scheme:
         dtxxxx = HUGE(time)
         DO i=index_p_nf_start,index_p_nf_stop
            if(time .GE. beam_data%switch_on_off_delay(i))CYCLE
            IF(beam_data%switch_on_off_delay(i) .GT. 
     .                        timnew + time_tol)CYCLE
            x = beam_data%switch_on_off_delay(i) - time + time_tol
            IF(x .Lt. dtxxxx)THEN
               dtxxxx = x
               isavexxxx_event = i
            ENDIF
         ENDDO
         dtxx = MIN(dtxx,dtxxxx)
         IF(dtxx == dtxxxx)isave_event = isavexxxx_event
      ENDIF



c  --------------------------------------------------------------------
c     if time+dt (=timnew) .lt. the time of the next event
c     (plus time_tol) then dtxx = large number (1.e20)
c     otherwise 
c     dtxx is now the time interval defined so that the next
c     event occurs at time + dtxx
c  ---------------------------------------------------------------------




c  -------------------nubeam --------------------------------------------
c  --
c  -----------------------------------------------------------------------
      IF(use_nubeam)THEN
c         print *,'checkdt nubeam,time,timnew,dt =',time,timnew,dt
c         print *,'x,xx,dtxx =',x,xx,dtxx
cJMP         print *,'beam_data%beam_power_rise_time,nubeam_steps =',
cJMP     .       beam_data%beam_power_rise_time,nubeam_steps
c        if(allocated(nubeam_calls))print *,'nubeam_calls =',
c     .                  nubeam_calls
cJMP        if(associated(nubeam_calls))print *,'nubeam_calls =',
cJMP     .                  nubeam_calls
cJMP         print *,'beam_data%beam_times =',beam_data%beam_times
cJMP         print *,'beam_data%beam_switch_times =',
cJMP     .   beam_data%beam_switch_times
c         print *,'nubeam_dt =',nubeam_dt
c       print *,'SIZE(beam_data%beam_times) =',SIZE(beam_data%beam_times)


c        get the intervals lt,ln in which time and timnew lie:
         ln =0
         lt =0
         xxx = 1.e30
         do j=1,SIZE(beam_data%beam_times)-1
            if(beam_data%beam_times(j) .le. time .and.
     .         time .le. beam_data%beam_times(j+1))lt =j
            if(beam_data%beam_times(j) .le. timnew .and.
     .         timnew .le. beam_data%beam_times(j+1))ln =j
             if(timnew .gt.  
     .           beam_data%beam_times(SIZE(beam_data%beam_times)))
     .           ln = SIZE(beam_data%beam_times)
            IF(nubeam_steps .gt. 1)then
              x= beam_data%beam_times(j)-nubeam_calls(nubeam_steps-1)
              if(x .gt.0.0)then
                 xxx = min(x,xxx)
              endif
            ENDIF
         enddo
c 

         IF(lt*ln .eq.0 .or. ln .lt. lt)THEN
            print *,' Error in nubeam chekdt:'
             print *,'lt,ln =',lt,ln
            print *,'time,timnew =',time,timnew
            print *,'beam_data%beam_times =',beam_data%beam_times
            print *,'beam_data%beam_switch_times =',
     .                                beam_data%beam_switch_times
            print *,'nubeam_dt =',nubeam_dt
            CALL STOP('%%%%CHEKDT NUBEAM DT ERROR',1)
          ELSE   ! LN >=  LT
             x = beam_data%beam_times(lt+1) -time
          ENDIF 


c        make sure that time+dt doesnt step beyond the time
c        that nubeam has been run up to. nubeam will be called if
c        ABS(time +dt (= timnew)  -  nubeam_calls(nubeam_steps-1) ) .lt.
c        3.*time_tol    
c        change nubeam_dt from its input value only if we would otherwsie
c        cross a beam switching boundary given in beam_data%beam_times
         IF(nubeam_steps .gt. 1)THEN
               xx = nubeam_calls(nubeam_steps-1) -time
               x=MIN(x,xx)
               IF(ABS(xx) .lt.  1.5*time_tol)
     .               nubeam_dt = MIN(xxx,nubeam_dt)
         ENDIF
c        x is now the time interval such that time +x = z
c        where z is the minimum of( next beam switching time,
c        time up to which nubeam has been run.
         IF(x .lt. dtxx .or. ABS(x-dtxx) .lt. 1.5*time_tol)THEN
               dtxx =x
c               IF(ln .gt. lt) nubeam_dt = MIN(dtxx,nubeam_dt)
               isave_event = index_p_nf_start-1
         ENDIF
        
 1460    continue

      ENDIF
       nubeam_dt = MAX(nubeam_dt,beam_data%beam_power_rise_time*.5)
       print *,'end nubeam chekdt,nubeam_dt =',nubeam_dt
!--------------------end nubeam ---------------------------------------




c      if (dtxx .gt. 1.0e10 )  go to 5040 !no need to change dt
c       write(969,FMT='("dtxx -dt =",1pe16.8)')dtxx - dt
c       write(969,FMT='("dt,time =",3(2x,1pe16.8))')dt,time ! 88888999
      if (dtxx .gt. dt)  go to 5040 !  HSJ 03/04/05  no need to change dt
      if (tportvb .ge. 1)
     .    write (*, '(" dt changed due to event ", a /
     .                " new dt = ", 1pe16.8, " old dt ", 1pe16.8,
     .                " time = ", 1pe16.8)')
     .                  event_array(isave_event), dtxx, dt, time
      write(*,'("event,time,timnew,new dt ",a,3(2x,1pe16.8),5(2x,i3))')
     .event_array(isave_event),time,time+dtxx,dtxx,(itest(i),i=1,5)
c
c ----------------------------------------------------------------------
c assign new value to dt
c ----------------------------------------------------------------------
c
      dtxx = MIN(dtxx,dtmax)
c      dtxx =MAX(dtxx,dtmin)
      dthold = dt
      dt     = dtxx
      if (dt .lt. dtevmin)  dt = MIN (dthold, dtevmin)
      if(nubeam_dt .lt. dt)nubeam_dt = dt
      timnew = time + dt
      do 5030 i=1,nevents
c          if(time .gt. 1.4998 .and. time .lt. 1.5002 
c     .        .AND. i .gt. 10 .and. i .lt. 26)THEN
c             write(969,FMT='("time,timnew ,dt,i=",3(2x,1pe16.8),i5)')
c     .            time,timnew ,dt,i
c             dtp = timevent(i) -timnew
c            write(969,FMT='("timevent(i),timmnew =",3(2x,1pe16.8),i5)')
c     .        timevent(i),timnew ,dtp,i
c             dtp = time - timevent(i)
c             write(969,FMT='("timevent(i),time =",3(2x,1pe16.8),i5)')
c     .        timevent(i),time ,dtp,i
c             write(969,FMT='("ievent before =",i5,2x,i5,/)')ievent(i),i
c          ENDIF
         if (timevent(i) .gt. timnew     )  go to 5030 
         if( time - timevent(i) .GT. time_tol) go to 5030
           ievent(i) = 1
 5030   continue

c
 5040 do i=1,krf
        if (irf(i) .eq. 1)  irf(i) = 0
        if (irf(i) .eq. 2)  irf(i) = 3
        if (ievent( 9+2*i) .eq. 1 .and. irf(i) .eq. 0)  irf(i) = 1
        if (ievent(10+2*i) .eq. 1                    )  irf(i) = 2
      end do
c

      if ( ibeam .eq. 1)   ibeam = 0
      if ( ibeam .eq. 2)   ibeam = 3
      if (ievent( 9    ) .eq. 1 .and.  ibeam .eq. 0)  ibeam  = 1
      if (ievent(10    ) .eq. 1                    )  ibeam  = 2
      inubplt = ievent(4)
c
c --- the following line is no longer used.
c --- irfplt is now an input variable
c     irfplt = ievent(7)
      ipelet = ievent(8)
c
c ----------------------------------------------------------------------
c check mprt, mplot and mbug
c ----------------------------------------------------------------------
c
      nnext = n + 1
csjt**iplot = ievent(2) + 2
      iprt  = ievent(3)
      if (mprt .eq. 0)  go to 5050
      i     = (nnext/mprt)*mprt
      if (i .eq. nnext)  iprt = 1
 5050 if (mplot .eq. 0)  go to 5080
      i     = (nnext/mplot)*mplot
      if (i .eq. nnext)  iplot = 3
 5080 if (ievent(6) .eq. 1) iplot = -3
      if (ievent(1) .eq. 1) then
        iprt  =  1
        iplt  =  3
      end if



      return
c
      end

      subroutine cheku (u, nk, nj, kk, delav, delmax, kmax, jmax,
     .                  delit, iter, ihalve, iangrot)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c     check that all values of 'u' are non-negative.  if they are not,
c     set ihalve such that the time step will be halved in the calling
c     routine.  allow negative angular rotation so check only for
c     k = 1, nk-iangrot:
c     allow negative rbp, so check only for k = 1,nk-iangrot-1:
c     on option set any remaining negative values to near zero and don't
c     return an error flag.
c ----------------------------------------------------------------------
c
      dimension u(kk,*)
c
      ihalve = 0
      do   j=1,nj
        do k=1,nk-1-iangrot        ! -1 allows for negative bp
          if (u(k,j) .lt. 0) then
            ihalve = 1
            delav  = 0.0
            delmax = 1.0
            kmax   = k
            jmax   = j
            delit  = 0.0
            iter   = 0
            return
          end if
        end do
      end do
      return
c
      end

      subroutine chkcon (is, tnion, tnerr, tnrerr)
c


c
c ----------------------------------------------------------------------
c     CHKCON checks on conservation of particles.
c     The density returned by SOLVE is integrated and added to
c     the number of new particles added to the sytem.
c
c     The error (in number of particles) is computed based on:
c     ipcons = 0: assume recycling of neutrals, or
c     ipcons = 1: maintain constant number of ions in the plasma.
c
c     snadd  =       number of ions added this step
c     snaddt = total number of ions added previously
c ----------------------------------------------------------------------
c
      USE param
      USE ions
      USE neut
      USE solcon
      USE soln
      USE numbrs
      USE mesh
      USE sourc
      USE machin
      USE geom
      USE flags
      USE bd_condtn,only : bctime,ub,fluxb,
     .    ub_save,ub_rho_edge,bctime_zone
      implicit  integer (i-n), real*8 (a-h, o-z)
c      include 'param.i'
c      include 'numbrs.i'
c      include 'sourc.i'
c      include 'soln.i'
c      include 'mesh.i'
c      include 'geom.i'
c      include 'machin.i'
c      include 'solcon.i'
c      include 'bcon.i'
c      include 'neut.i'
c      include 'flags.i'
c      include 'ions.i'
c
      dimension xu(kj)
c
      do j=1,nj
        xu(j) = u(is,j)
      end do
      call trapv (r,xu,hcap,nj,tnion)
      tnion  = volfac * tnion    ! current # ions of type is in system
c
      tnerr  = 0.0
      tnrerr = 0.0
****  if (iten .eq. 0)  return
      if (itran(is) .eq. 0) return   ! HSJ 1/26/96
c
      tnionx = tnion0(is)            ! tnion0 is total # ions, type is,
c                                      in the system at time=time0 (or
c                                      at the start of any transport/mhd
c                                      if mhd cals are being done)
      if (ipcons(is) .eq. 0)         ! recycling of neutrals
     .  tnionx = tnionx + snaddt(is) + snadd(is) + sngas(is) + sn2d
c
c     tnionx= # at time0 + number added since time0 (up to last converged
c             time point) + current number being added                HSJ
c
      tnerr  = tnionx - tnion
      if (tnionx .ne. 0.0)  tnrerr = -tnerr/tnionx
      !print *,'jmp.den.tmp',fluxn(1)
      return
c
      end

      subroutine chkflx (is)
c

c
c ----------------------------------------------------------------------
c     Chkflx checks the neutral influx at the plasma boundary
c     to be sure it is above some minimum.  We don't want this
c     flux to become negative while trying to satisfy the
c     ipcons = 1 option.
c     We also don't allow the neutral flux to be zero in any case:
c     this complicates the 2 neutral species case too much.
c     Ways to get a non-zero neutral flux:
c     1) neutral gas injection (gasflx)
c     2) iten = 0: non-zero density profile, taupin<infinity
c     3) iten = 1: non-zero, non-flat density profile
c ----------------------------------------------------------------------
c
      USE param
      USE ions
      USE io
      USE neut
      USE solcon
      USE sourc
      implicit  integer (i-n), real*8 (a-h, o-z)
c      include 'param.i'
c      include 'io.i'
c      include 'ions.i'
c      include 'neut.i'
c      include 'solcon.i'
c      include 'sourc.i'
c
      if (ipcons(is) .eq. 0     )  go to 200
      if (flxadd(is) .eq. 0.0   )  go to 200
      if (fluxn (is) .ge. flxmin)  return
c
c     The neutral influx is too small for ipcons = 1
c
      if (ipfail(is) .ne. 0)  go to 100
      ipfail(is) = 1
      write (nout, 8000) is,namep(is),time,n,flxadd(is),fluxn(is),flxmin
      write (nqik, 8000) is,namep(is),time,n,flxadd(is),fluxn(is),flxmin
      write (ncrt, 8000) is,namep(is),time,n,flxadd(is),fluxn(is),flxmin
 8000 format (/' ERROR: For species',i2,' (',a2,'), the neutral influx'/
     . ' at the edge is too small to permit a constant number of ions' /
     .'  time = ',f8.3,5x,'n =',i5/
     .'  flxadd = ',e11.3,'  fluxn=',e11.3,'  flxmin=',e11.3,'/cm**2-s'/
     .' fluxn set to flxmin, and particle conservation disabled.'/)
c
  100 fluxn(is) = flxmin
      return
c
  200 if (fluxn(is) .gt. 0.0)  return
c
c     The neutral flux is zero, or worse
c
      write  (nout, 8010)  is, namep(is)
      write  (nqik, 8010)  is, namep(is)
      write  (ncrt, 8010)  is, namep(is)
 8010 format (/ ' ERROR: for species', i2, ' (', a2,
     .          '), the neutral influx' /
     .          ' at the edge is zero'  /)
      print *,'jmp.den',fluxn(is) !jmp.den
      call STOP ('subroutine CHKFLX: zero or negative flux', 81)
c
      end

      subroutine co2cal
c

c
c ----------------------------------------------------------------------
c     CO2CAL calculates the line-averaged electron density
c     along several tangential chords.  This simulates the
c     measurements of the CO2 interferometer array.
c ----------------------------------------------------------------------
c
      USE param
      USE soln
      USE numbrs
      USE mesh
      USE machin
      USE geom
      implicit  integer (i-n), real*8 (a-h, o-z)
c      include 'param.i'
      include 'co2.i'
c      include 'geom.i'
c      include 'machin.i'
c      include 'mesh.i'
c      include 'numbrs.i'
c      include 'soln.i'
c
      call zeroa (denco2,nco2)
      call zeroa (patco2,nco2)
      do 200 ico2=1,nco2
      if (rtco2(ico2) .ge. rmax)  go to 200
c
c          calculate start and end points for line integral
c
      x0 = rtco2(ico2)
      y0 = -SQRT (rmax**2 - x0**2)
      z0 = 89.0
      xf = x0
      yf = 0.0
      zf = z0
      if (rtco2(ico2) .gt. rin)  go to 100
      yf = -SQRT (rin**2 - xf**2)
c
c     initialize line integral; use 1 cm step size
c
  100 nc = yf - y0 + 1.0
      if (nc .lt. 30) nc = 30
      dy = (yf-y0)/(nc-1.0)
      y1 = y0
      call getrho(x0,y1,z0,rho1)
      call interp(rho1,r,nj,ene,ene1)
      sum = 0.0
c
c          evaluate line integral
c
      do i=2,nc
        y2   = y1 + dy
        call getrho(x0,y2,z0,rho2)
        call interp(rho2,r,nj,ene,ene2)
        sum  = sum + 0.5*(ene1+ene2)*dy
        if (ene1+ene2 .gt. 0.0)  patco2(ico2) = patco2(ico2) + dy
        y1   = y2
        ene1 = ene2
      end do
c
      denco2(ico2) = sum/patco2(ico2)
      if (rtco2(ico2) .gt. 1.0) patco2(ico2) = 2.0 * patco2(ico2)
  200 continue
      return
c
      end

      subroutine copya (old_array, new_array, number)
c
      implicit none
c
c --- subroutine COPYA copies the array OLD_ARRAY into the array NEW_ARRAY
c
      integer  index, number
      real*8   old_array(*), new_array(*)
c
      do index=1,number
        new_array(index) = old_array(index)
      end do
      return
c
      end

      subroutine copyas (old_array, stride_old,new_array,stride_new,
     .                   number)
c
      implicit none
c
c ----------------------------------------------------------------------
c     subroutine COPYAS
c     copies the array OLD_ARRAY into the array NEW_ARRAY
c     old array has stride,stride_old,new array has stride stride_new
c ----------------------------------------------------------------------
c
      integer  index,indexn, number, stride_old,stride_new, i
      real*8   old_array(*), new_array(*)
c
      index=1
      indexn=1
      do i=1,number
        new_array(indexn) = old_array(index)
        index=index+stride_old
        indexn=indexn+stride_new
      end do
      return
c
      end

      subroutine correct_u (delta_time)
c

c
c ----------------------------------------------------------------------
c     subroutine determines the density of analysis mode primary
c     ion when ifus = -1
c     only u(ii,j) and dnidt(j,ii) are used subsequently everything
c     else is scratch storage
c ----------------------------------------------------------------------
c
      USE param
      USE fusion
      USE io
      USE ions
      USE nub2
      USE solcon
      USE soln
      USE numbrs
      USE flags
      USE bd_condtn
      USE etc
      implicit  integer (i-n), real*8 (a-h, o-z)

      include 'storage.i'      ! xdum,etc.
c
      dimension    ene_timnew(kj),zeff_timnew(kj),en_timnew(kj,kion)
      equivalence (ene_timnew(1),zdum(1))
      equivalence (zeff_timnew(1),wdum(1))
      equivalence (en_timnew(1,1),tdum(1))
c
c --- load simulation density into en_timnew:
c
      do j=1,nj
         do ii=1,nprim
             if (itran(ii) .eq. 1)  en_timnew(j,ii)=u(ii,j)
         end do
      end do
c
c --- get ne and zeff at the new time:
c
      timesave=time   ! solcon.i stores these values
      time=timnew     ! tsplinew uses time
c
c --- spline time-dependent ne assumed:
c
        call tsplinew (enein,renein,ene_timnew,knotsene,kprim+kimp+1,
     .                 xdum,ydum)
c --- spline time-dependent zeff assumed:
c
            call tsplinew (zeffin, rzeffin, zeff_timnew, knotszef,
     .                     kprim+kimp+4, xdum, ydum)
c
      do j=1,nj
c
c --- given the simulation mode primary ion density (in en_timnew(j,isym))
c --- and given ene_timnew,zeff_timnew at this time,get the analysis
c --- mode primary ion density,en_timnew(j,ianal) at this time:
c
      call get_density (enbeam,enalp,ene_timnew,en_timnew,
     .                  zsq,zeff_timnew,z,itran,kj,j,nprim,nimp,
     .                  ihe,id,it,ncrt,nout, -1.0, time)
c
c --- get density of analysis mode primary ion density at the
c --- central time:
c
           do ii=1,nprim
               if (itran(ii) .eq. 0) then
                  u(ii,j) = en_timnew(j,ii)   ! density in analysis mode
                  en(j,ii)= theta*u(ii,j) + (1.0-theta)*usave(ii,j)
****              dnidt(j,ii)=(en(j,ii)-usave(ii,j))/delta_time
               end if
           end do
      end do
      time = timesave    ! restore time
      return
c
      end

      subroutine cubicextrp (f2, f3, f4, r2, r3, r4, f1, itype)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c --- subroutine gets the value of the function f at r(1) ( = 0)
c --- by extrapolation of the fitted quadratic or cubic segment:
c      f(r) = a*r**2+b              itype = 1
c or
c      f(r) = a*r**3+b*r**2+c       itype = 2
c or   the finite difference form
c      f(1) = (4/3)*f(2)-(1/3)*f(3) itype = 3
c --- note that in all cases df/dr = 0 at r = 0.0
c
c --- input
c  itype               = 1 selects quadratic
c                      = 2 selects cubic
c  f2                  value of f at r2
c  f3                                r3
c  r2                  independent variable
c  r3
c if itype = 2 then  the following also have to be input
c f4
c r4
c --- output
c f1               value of f at r1,where it is assumed that r1 = 0
c                  with df/dr = 0.0 at r = 0.0
c ------------------------------------------------------------------ HSJ
c
      if      (itype .eq. 1) then
          a    = (f3-f2)/(r3**2-r2**2)
          f1   = f2-a*r2**2
      else if (itype .eq. 2) then
          det  = (r4**3)*(r3**2-r2**2)-(r4**2)*(r3**3-r2**3)
     .         + (r3**3)*(r2**2)-(r3**2)*(r2**3)
          detc = (r4**3)*(f2*r3**2-f3*r2**2)-(r4**2)*(f2*r3**3-f3*r2**3)
     .         +  f4*((r3**3)*(r2**2)-(r3**2)*(r2**3))
          f1   = detc/det
      else             ! itype = 3
          h1   = r2
          h2   = r3 - r2
          a1   = -(2.0*h1+h2)/(h1*(h1+h2))
          b1   =  (h1+h2)/(h1*h2)
          c1   = -(h1/h2)/(h1+h2)
          f1   = -(b1*f2+c1*f3)/a1
      end if
c
      return
c
      end

      real*8 function cxra (ti, atw)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c     CXRA returns the charge exchange rate coefficient <sigma-v>
c     for any primary ion: hydrogenic and helium ions have to be
c     handled separately.
c ----------------------------------------------------------------------
c
      if (atw .gt. 3.0)  go to 10
c
c          hydrogenic
c
      cxra = cxr(ti/atw)
      return
c
c          helium
c
   10 cxra = cxrhe(ti)
      return
c
      end

      real*8 function cxr (x)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c this function calculates the charge exchange rate for hydrogen atoms
c interacting with protons in units of cm**3/s.
c x is in units of keV for 1.0e-3 .le. x .le. 100, the
c the formula is taken from the paper by r.l. freeman and e.m. jones
c clm-r 137 culham laboratory 1974.
c for x .lt. 1.0e-3, a rate coefficient derived from an analytic average
c over the approximate cross section  sigma = 0.6937e-14*(1.0-0.155*LOG10
c (e/1ev))**2 is used.  this cross section is an approximation to that
c given by riviere, nuclear fusion 11,363(1971).
c ----------------------------------------------------------------------
c
      if (x .lt.   1.0e-3)  go to 10
      if (x .gt. 100     )  go to 20
      tc  = LOG (x)+6.9077553
      dum = 0.2530205e-4-tc*0.8230751e-6
      tc  = -0.1841757e2+tc*(0.528295-tc*(0.2200477-tc*(0.9750192e-1-tc*
     .      (0.1749183e-1-tc*(0.4954298e-3+tc*(0.2174910e-3-tc*dum))))))
      tc  = EXP (tc)
      cxr = tc
      return
c
   10 tc  = 0.50654 - 6.7316e-2 * LOG  (x)
      cxr =           3.3340e-7 * SQRT (x) * (tc*tc + 7.454e-3)
      return
c
   20 cxr = 0.0
      return
c
      end

      real*8 function cxrhe (tkev)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c this function calculates the double charge exchange
c rate for helium as a function of temperature in keV
c equations were derived by doing
c a least squares curve fit on data provided by:
c     Oak Ridge national labratory
c     atomic data for controlled fusion
c          ornl-dwg 76-6797
c
c vel(thermal)  = sqr(2t/mass(he))
c 5.4676655e-11 = sqr(2 /mass(he))
c t = temperature in ergs
c mass is in grams
c ----------------------------------------------------------------------
c
      vth   = 5.4676655e-11 * SQRT (tkev/1.6e-12)
      atkev = LOG10 (tkev)
      if (tkev .gt. 1.587e3) ratehe = 0.0
      if (tkev .lt. 0.05   ) ratehe = 4.0e-16
      if (tkev .le. 90.0 .and. tkev .ge. 0.05)
     .  ratehe = vth*(10**(-15.6868+atkev*(-0.192101+atkev*2.21944e-2)))
      if (tkev .gt. 90.0 .and. tkev .le. 1.587e3)
     .  ratehe = vth*(10**(-24.9253+atkev*(9.07495-atkev*2.31366)))
      cxrhe    = ratehe
      return
c
      end

      real*8 function cxrv (x)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c  NOTE THAT FUNCTION CXR USES THE MAXWELLIAN AVERAGE SIGMA V,
c       THIS FUNCTION ( ie CXRV ) USES SIGMA*V INSTEAD.
c       HSJ 6/5/96
c this function calculates the charge exchange rate for hydrogen atoms
c interacting with protons in units of cm**3/s.
c x is in units of keV.
c for 1.0e-3 .le. x .le. 100, the
c the formula is taken from the paper by r.l. freeman and e.m. jones
c clm-r 137 culham laboratory 1974
c for x .lt. 1.0e-3, a rate coefficient derived from an analytic average
c over the approximate cross section  sigma = 0.6937e-14*(1.0-0.155*LOG10
c (e/1ev))**2 is used.  this cross section is an approximation to that
c given by riviere, nuclear fusion 11,363(1971).
c ----------------------------------------------------------------------
c
      if      (x .lt.   1.0e-3) then
        call STOP ('function CXRV: X is out of range', 225)
      else if (x .gt. 100.0   )  then
        cxrv = 0.0
      else
        ev   = 1000.0*x
        cxrv = 0.6937e-14 * (1.0 - 0.156 * LOG10 (ev))**2 /
     .        (1.0 + 0.1112e-14 * (ev**3.3))
      end if
      return
c
      end

      subroutine ddfbm (te, e0, ne, nd, nddot, rbd)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c     this subroutine calculates neutron production due to deuterium
c     neutral beams.  discussion of this computation can be found in
c     appendix b of Steve Scott's thesis entitled "an experimental
c     investigation of magnetic field ripple effects on tokamak plasmas"
c
c        input quantities:
c
c     te    - electron temperature in keV
c     e0    - injection energy of beam ion in keV
c     ne    - thermal electron density in cm**-3
c     nd    - thermal deuterium density in cm**-3
c     nddot - number of beam deuterons per sec per cm**-3
c
c        output quantities:
c
c     rbd - number of neutrons produced per sec per cm**-3
c ----------------------------------------------------------------------
c
      real*8  nddot, ne, nd, nude
c
c     constant for cm/sec to keV
c
      k = 3.1e07
c
c     sigmadd constant
c
      c = 2.18e-22
c
c     sigmadd constant
c
      b = 47.88
c
c     rate (1/sec) of energy transfer fast d to e
c
      nude = 7.59*(ne/1.0e13)/te**1.5
c
c     critical energy (keV) for deuterium
c
      ecd = 18.5 * te
      rbd =  2.0 * nddot * c * k * nd * EXP (-b / SQRT (e0)) /
     .            (nude * (1.0 + (ecd/e0)**1.5) * b)
      return
c
      end

      subroutine ddfknc (te, e0, ne, nd, nhdot, rkd)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c     this subroutine calculates neutron production due to large angle
c     collisions between fast beam hydrogen ions and thermal deuterons.
c     discussion of this computation can be found in appendix b of Steve
c     Scott's thesis entitled, "an experimental investigation of
c     magnetic field ripple effects on tokamak plasmas".
c
c        input quantities:
c
c   te-electron temperature in keV
c   e0-injection energy of beam ion in keV
c   ne-thermal electron density in cm**-3
c   nd-thermal deuterium density in cm**-3
c   nhdot-number of beam hydrogens per sec per cm**-3
c
c        output quantities
c
c   rkd-number of neutrons produced per sec per cm**-3
c
      real*8  nhdot, ne, nd, nd13, nude, nuhe
c
      nude = 7.59*(ne/1.0e13)/te**1.5
c
c  rate (1/sec) of energy transfer fast d to e
c
      nuhe = 2.0 * nude
c
c  rate (1/sec) of energy transfer fast h to e
c
      ecd = 18.5*te
c
c  critical energy (keV) for deuterium
c
      ech = 9.2*te
c
c  critical energy (keV) for hydrogen
c
      b = 47.88
c
c  constant for sigmadd
c
      em = 8.0*e0/9.0
c
c  max transfer energy (keV) h to d
c
      nd13  = nd/1.0e13
      x1    = ech/e0
      x2    = ecd/em
      term1 = nuhe*nude
      term2 = 7.29*nhdot*nd13**2 * EXP (-b / SQRT (em))
      term3 = (1.0 + x2**1.5)*(1.0 + x1**1.5)
      term4 = (1.0 + 2.0 * SQRT (em)/b)/(SQRT (em)*b**3)
      rkd   = (1.0/term1)*(term2/term3)*term4
      return
c
      end

      subroutine ddfust (tix, enx, ddfus, ddnt, time, bpol, totcur1)
c
c
c ----------------------------------------------------------------------
c     calculate the number of neutrons produced by D-D fusion if
c     deuterium is present as a primary ion.  the <sigmav> formula is
c     taken from the NRL Formulary.  the rate coefficient is taken as
c     sigmav/4 since about 50% of D-D fusions produce a neutron and both
c     reactants are D.
c
c     Updated 2/22/95 by HSJ:
c          new rate coefficient for thermal d(d,n)he3 reaction
c          is now available. The selection is made with iddfusrate
c          for the thermal rate as follows (for the beam-thermal
c          rate; see iddfusb below):
c
c             iddfusrate   =   1  new cross section (highly recommended
c                                 but not the default since we wish
c                                 to remain comensurate with other input)
c                          =   0  default, use old NRL rate coefficient
c                                 A graphical comparison of the
c                                 two rate coefficients is available in
c                                  /u/stjohn/ONETWOPAGES/neutron_rate/
c                                                thermal_rate/ddnhe.?
c                                      where ? is tex and/or ps
c
c          neutron production due to the fast deuteron distribution
c          (sourced by the beam) is controlled by the following switches
c             iddfusb      =   0   default, use the old
c                                  (i.e., previously existing)
c                                  approximation to the neutron rate for the
c                                  d(df,n)he3 reaction (df=fast deutron)
c                                  This approximation typically will
c                                  overestimate the neutron rate at low ion
c                                  temperatures and underestimate at high
c                                  (i.e., > 5 keV) temperatures. Bulk motion
c                                  of ions is negelected in this model.
c                          =   1   user selected,do a more detailed calc
c                                  of the d(df,n)he3 rate as selected by the
c                                  switches described in subroutine ddfbm1.
c                                  (the Bosch & Hale
c                                  cross sections are used throughout)
c          icalc_cxfactor          (see cray101.f)
c
c          The new RATE COEFFICIENT AND CROSS SECTIONS are from
c             Bosch & Hale, Nuc. Fus., vol32, no.4 (1992) 611
c
c     OTHER INPUT:
c     iddfus   =  0    no deuterium in system,therefore no d(d,n)he3 reactio
c              =  1    thermal species is d
c              =  2    thermal species is dt, use fraction of dt density
c                      given by fd  to specify d density
c     fd               fraction of deuterium in dt mixture
c
c     tix(j)           j=1,2..nj, ion temperature (must be in keV)
c
c     enx(j,i)         j=1,2..nj,density of ion species i,#/cm**3
c     time             current time,sec
c    bpol
c    totcur1           bpol and totcur1 are used to determine correct
c                      signs on some angles
c
c     OUTPUT:
c      thermal d,d reactants:
c             ddfus(j)   local reaction (i.e., neutron production) rate,
c                        #/cm**3/sec
c             ddnt       volume integrated rate, #/sec
c
c ----------------------------------------------------------------------
c
      USE param
      USE fusion
      USE ions
      USE neut              ! pick up neutral density here
      USE nub
      USE nub2
      USE soln
      USE numbrs
      USE mesh
      USE machin            ! btor
      USE geom
      USE colrate
      USE tordlrot
      implicit  integer (i-n), real*8 (a-h, o-z)

      dimension      tix(kj), enx(kj,kion), ddfus(kj),
     .               ddbeamint(kj), bpol(kj)
c
      data xkeverg / 6.241e+08 /
      data xmassp  / 1.673e-24 /      ! grams
c
      data isetsigvr, tionmax_table, ebeammax_table
     .    /0        , 0.0          , 85.0          / ! controls calc..
c                                                    ..sigvrddn table
c
      real*8  mrcsq            ! for new Bosch & Hale rate coefficient
c
      data    c1, c2, c3, c5, b_g, mrcsq
     .      / 5.43360e-12,  5.85778e-3,      7.68222e-3,
     .       -2.96400e-06, 31.3970    , 937814.0 /
c
c     other parameters required for neutron rate calculations
c     when iddfusb = 1
c
      ddnt = 0.0

      if (iddfus .eq. 0)  return
c
      mass_deut       = 2.0 * xmassp   ! approx mass of deuteron, grams
      mass_beam       = atw_beam * xmassp
      ibeam_chg_exchg = ibcx           ! load colrate.i
      icxcalc         = icalc_cxfactor ! move it into colrate.i
      jsavef          = 0
      ddnt            = 0.0
      vbeammin_table  = 0.0            ! range of sigvrddn table, cm/sec
c                                        may be adjusted below
      if (vbeammax_table .le. 0.0)
     .  vbeammax_table  = SQRT (2.0 * ebeammax_table
     .                   / (xkeverg * mass_beam))    ! cm/sec
      vionmin_table     = 0.0
      if (vionmax_table .le. 0.0)
     .    vionmax_table = 4.0026e-5 * SQRT (1.5 * ti(1) / mass_beam)
c
      if (iddfus .eq. 0)  return
c
      do j=1,nj
c
c        neutrons produced by bulk plasma d-d fusion
c
         if (iddfus .eq. 1)  endfus = enx(j,id)
         if (iddfus .eq. 2)  endfus = fd*enx(j,idt)
         if (iddfus .eq. 3)  endfus = enx(j,id)
         if (iddfus .eq. 4)  return ! no deuterium in system
c
         if (iddfusrate .eq. 0) then
           tmot   = 1.0 / (tix(j)**(1.0/3.0))
           sigmav = 2.33e-14 * (tmot**2) * EXP (-18.76*tmot)
           sigmav = sigmav/4.0   ! sigma counts d(d,n)he3 and d(d,p)t
         else                    ! use Bosch & Hale rate coefficient
           theta  = tix(j)*C2/(1.0+tix(j)*(C3+tix(j)*C5))
           theta  = tix(j)/(1.0-theta)
           xsi    = (B_g**2/(4.0*theta))**(0.33333333334)
           sigmav = C1 * theta * SQRT (xsi/(mrcsq*tix(j)**3))
     .                         *  EXP (-3.0*xsi)
           sigmav = sigmav/2.0   ! sigma is specific to d(d,n)he3 here
         end if
c
         ddfus(j) = (endfus**2) * sigmav   ! sigmav is in cm**3/sec
c
c        neutrons produced by beam-enhanced tail (knock-on) and beam
c        fast deuterons.  note:  the model assumes a classical, steady-
c        state distribution of fast ions.  transients are modeled by
c        using enb/taupb as the fast ion source.
c
c        The effect of fast ion charge exchange is included through
c        taupb, which is defined as
c                                    taupb = taus * encap
c        where taus is the usual Spitzer slowing down time of fast ions
c        on electrons and encap is the factor that accounts for charge
c        exchange losses.   HSJ
c
          ddknck(j)  = 0.0
          ddbeam(j)  = 0.0
          ecritlocal = ecrit(j)  ! keV, stored in colrate.i
          erot       = xkeverg*0.5*atw_beam*xmassp*vionz(j)**2
c
c         let Vf be fast ion velocity,Vr thermal ion velocity. Then
c         e0rel=0.5*mf*(Vf-Vr)dot(Vf-Vr)=ef+(mf/mi)erot
c                         -mag(Vf)*mag(Vr)*mf*cos (theta)
c         But mag(Vf)*cos (theta)*mf = momentum of fast ion in direction
c         of Vr, i.e., in toroidal direction = pb0(j,jc,jb).
c
          tauslocal = taus(j)  ! sec, stored in colrate.i, used in cxint
          do 105 jb=1,nbeams
          do 105 jc=1,3
              if (enb(j,jc,jb) .eq. 0.0)  go to 105
              sbh   = (1.0-fdbeam)*enb(j,jc,jb)/taupb(j,jc,jb)
              e0rel = ABS (ebkev(jb)/jc + erot - xkeverg*pb0(j,jc,jb)
     .                                                   *vionz(j))
              call ddfknc(te(j),e0rel,ene(j),endfus,sbh,ddkno)
              ddknck(j) = ddknck(j) + ddkno
c
c ---------------- beam-thermal neutron calculations -------------------
c
              if (iddfusb .eq. 0) then              ! old simple model..
c                                                   ..NRL cross section
                    sbd = fdbeam*enb(j,jc,jb)/taupb(j,jc,jb)
                    call ddfbm (te(j), e0rel, ene(j), endfus, sbd, ddbm)
c
              else if (iddfusb_s .eq. 0) then       ! iddfusb = 1 is..
c                                                   ..implicit here
c                   use new cross sections, simple model
c                   see description in:
c                   ONETWOPAGES/neutron_rate/beam_thermal_rate/nratef.ps
c                   fast ion source (sdot*taus) with or without
c                   charge exchange (depends on setting of ibcx
c                   in beam package and is already implicit in taupb)
c
                    sbd   = fdbeam*enb(j,jc,jb)*taus(j)/taupb(j,jc,jb)
                    e0rel = ebkev(jb)/jc   ! neglect bulk,thermal motion
                    call ddfbm1 (e0rel, endfus, sbd, ddbm)
c
              else if (IABS (iddfusb_s) .eq. 1) then ! iddfusb   = +-1..
c                                                ..and iddfusb_s = 1
c
c                   new cross section, more complex models
c                   note that we do not assume asymptotic source rates here.
c                   Consequently there will be some lack of
c                   self consistence between slowing down rates and
c                   neutron rates during transient (i.e., beam on and off)
c                   periods
c
c                   things that depend only on j but are kept inside the inner
c                   loops so that it is clear when they are used:
c
                    if (j .ne. jsavef) then
                       jsavef = j
                       z1     = 0.0
                       z2     = 0.0
                       do i=1,nion
                           z12 = enx(j,i)*zsq(j,i)
                           z1  = z1+z12/atw(i)
                           z2  = z2+z12
                       end do
                       z1 = z1*atw_beam/ene(j)
                       z2 = z2/(z1*ene(j)) ! effective charge of bg ions
                       den_neut = 0.0
                       do i=1,2        ! nneut = 2, hardcoded everywhere
                         den_neut = den_neut + enn(j,i) ! load colrate.i
                       end do
c
c                      flux surface average cos angle between B and toroidal
c                      direction is given by cosvzb. Note that <Bp**2>=
c                      Bp0**2*gcap and <Bt**2>=(btor/fcap)**2*<(R0/R)**2>
c
                       btavesq = ((btor/fcap(j))**2)*r2cap(j)
                       bpavesq = (bpol(j)*gcap(j))**2
                       cosvzb  = SQRT (btavesq/(btavesq+bpavesq))
                       cosvzb  = cosvzb * btor / ABS (btor) ! SIGN fcn..
c                                                           ..universal?
                       if (iddfusb_bulk .eq. 1) then
                         v_ion_bulk = vionz(j)*cosvzb  ! bulk v toroid..
c                                                      ..in B dir
                       else
                         v_ion_bulk = 0.0
                         e0rel      = ebkev(jb)/jc     ! beam energy
                       end if
c
c                      vthion is energy at which a fast ion is considered
c                      thermalized and therefore no longer part of the fast
c                      ion distribution function:
c
                       vthion  = 4.0026e-5
     .                         * SQRT (2.0 * ti(j) / mass_beam) ! cm/sec
                       vthelec = 1.3256e9
     .                         * SQRT (2.0 * te(j)            ) ! cm/sec
                    end if
c
c                   sbpure = true (FREYA calced/saved) source rate,#/(cm**3sec):
c                   sbsav has effects of charge exchange and time dependence
c                   already in it, so must use sbpure if icalc_cxfactor = 1
c
                    if (icalc_cxfactor .eq. 0) then
                      sbd = fdbeam*enb(j,jc,jb)*taus(j)/taupb(j,jc,jb)
                    else if (icalc_cxfactor .eq. 1) then   ! source, sbd
c                                                    ..without cx factor
                      sbd = fdbeam*sbpure(j,jc,jb)*taus(j)
                    else              ! source, sbd, has cx factor in it
                      sbd = fdbeam*enb(j,jc,jb)
                    end if
                    beam_pitch_angle = zeta(j,jc,jb) * btor * totcur1
     .                                          / ABS (btor * totcur1)
c
                    vbeam  = SQRT (2.0 * e0rel
     .                          / (xkeverg*mass_beam)) ! cm/sec
                    vcrit  = 0.09 * ((z1 / 2.0)**0.33) * vthelec
                    vcrit3 = vcrit**3
c
c                   transient beam effects are done by using the
c                   correct support set of the fast ion distribution
c                   function (which is time-dependent):
c                   beam just turned on?
c
                    tsldown = LOG ((vbeam**3+vcrit**3)/vcrit**3)
                    tsldown = tsldown*taus(j)/3.0
                    if ((time .gt. beamon(1)) .and.
     .                  (time-beamon(1)) .le. tsldown) then
c
c                        slowing down distribution is only populated
c                        down to this speed
c
                         vfast_low_lim = (vbeam**3+vcrit**3)*
     .                     EXP (-3.0*(time-beamon(1))/taus(j))-vcrit**3
                         if (vfast_low_lim .gt. 0.0) then
                           vfast_low_lim = vfast_low_lim**(1.0/3.0)
                         else
                           vfast_low_lim = vthion
                         end if
                    else   ! beam has been on long enough to establish..
c                          ..asymptotic slowing-down distribution
c
                         vfast_low_lim = vthion
                    end if
c
c                   beam just turned off?:
c
                    if (time .gt. beam_end(1)) then
c
c                      fast ions from vbeam down to this speed have
c                      all slowed down past this speed:
c
                       vfast_up_lim=(vbeam**3+vcrit**3)*
     .                   EXP (-3.0*(time-beam_end(1))/taus(j))-vcrit**3
                       if (vfast_up_lim .gt. 0.0)
     .                     vfast_up_lim = vfast_up_lim**(1.0/3.0)
                       vfast_up_lim = MAX (vfast_low_lim,vfast_up_lim)
                    else if (time .lt. beam_end(1)) then  ! beam currently on
                      vfast_up_lim = vbeam
                    end if
c
                    if (icalc_cxfactor .ne. 1) then
c
c                     the source, sbd, contains transient effects
c                     so we do not modify the integration limits
c
                      vfast_up_lim  = vbeam
                      vfast_low_lim = vthion
                    end if
c
c                   parameters required but not explicitly passed to
c                   DDFBM2 were loaded into colrate.i above
c
                    if (iddfusb_s .eq. -1) then      ! not a user option
                      call ddfbm2 (ti(j),v_ion_bulk,beam_pitch_angle,
     .                             e0rel,vbeam,endfus,sbd,ddbm)
                    else       ! iddfusb_s = 1, fast method, user option
c
c                     set up the sigvrddn array if not yet set.
c                     this table should be applicable for the entire
c                     run and thus the range of thermal velocities
c                     considered should be large enough to include all
c                     possibilities expected. (The table will be
c                     recalculated if the range is exceeded so we
c                     recover gracefully but waste some computer time)
c
                      tion      = ti(j)                ! keV
                      expd      = 0.5 * mass_beam
     .                                * xkeverg / tion ! 1/(cm/sec)**2
                      vthuplim  = 2.0 * SQRT (1.0/expd)! int.limits over
                      vthlowlim = 0.0                  ! distribution
                      if (vthuplim .gt. vionmax_table) then
                        isetsigvr     = 0
                        betamin       = 0.9 * expd
                        vionmax_table = 2.0 * SQRT (1.0/betamin)! cm/sec
                      end if
                      if (vbeam .gt. vbeammax_table) then
                        vbeammax_table = 1.01 * vbeam           ! cm/sec
                        isetsigvr      = 0
                      end if
                      if (vionmax_table .ne. v1(n1v))  isetsigvr = 0
                      if (vionmin_table .ne. v1(1)  )  isetsigvr = 0
                      if (vthuplim .gt. v1(n1v)) then
                         isetsigvr     = 0
                         vionmax_table = 1.01 * vthuplim
                      end if
                      if (isetsigvr .eq. 0) then
                        isetsigvr = 1
c
c                       load sigvrddn(n2v,n1v) (stored in colrate.i):
c                       d(d,p)t reaction is not done
c
                        call sigintddn (vionmin_table,vionmax_table,
     .                                  vbeammin_table,vbeammax_table)
                      end if
                      call ddfbm3 (vthuplim,vthlowlim,vbeam,
     .                             endfus,sbd,ddbm)
                    end if
              end if

              ddbeam(j) = ddbeam(j) + ddbm

c
  105     continue      ! end loop over beams and beam energy components
****      write (6, '(a, i2, a, f16.4)') ' ddbeam(', j, ') =', ddbeam(j)
      end do            ! end loop over spatial mesh j
c
      call trapv (r, ddfus , hcap, nj, ddnt  )
      call trapv (r, ddknck, hcap, nj, ddknct)
      call trapv (r, ddbeam, hcap, nj, ddbmt )
      ddnt   = ddnt   * volfac              ! volfac = 4 * pisq * rmajor
      ddnthm = ddnt
      ddknct = ddknct * volfac
      ddbmt  = ddbmt  * volfac
      if (iddcal .eq. 1 .or. iddcal .eq. 3)  ddnt = ddnt + ddknct
      if (iddcal .eq. 2 .or. iddcal .eq. 3)  ddnt = ddnt + ddbmt
      call trap3 (r, ddbeam, nj, hcap, volfac, ddbeamint)
      if (ddbeamint(nj) .ne. 0.0) then
        do j=1,nj
          ddbeamint(j) = ddbeamint(j) / ddbeamint(nj)
        end do
      end if
      return
c
      end

      real*8 function ddnrate (ti)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ------------------------------------------------- 11/20/95 --- HSJ ---
c             returns rate(cm**3/sec) of d(d,n)he3 reaction
c             (the division by 2 for collisions between like
c              particles is NOT accounted for here).
c             New Bosch & Hale rate coefficient:
c             Bosch & Hale, Nuc. Fus., vol32, no.4 (1992) 611
c
      real*8 mrcsq
c
c     data for d(d,n)he3
c
      data  c1, c2, c3, c5, b_gsq, mrcsq
     .    / 5.43360e-12,  5.85778e-3, 7.68222e-3,
     .     -2.96400e-06, 985.7716, 937814.0 /
c
      if (ti  .lt. 0.2 ) then
        ddnrate = 0.0
      else if (ti .le. 100.0) then
        theta   = ti*C2/(1.0+ti*(C3+ti*C5))
        theta   = ti/(1.0-theta)
        xsi     = (B_gsq/(4.0*theta))**(0.33333333334)
        ddnrate = C1 * theta * SQRT (xsi/(mrcsq*ti**3))
     .                       *  EXP (-3.0*xsi)
      else
        tidum = 100.
        theta   = tidum*C2/(1.0+tidum*(C3+tidum*C5))
        theta   = tidum/(1.0-theta)
        xsi     = (B_gsq/(4.0*theta))**(0.33333333334)
        ddnrate = C1 * theta * SQRT (xsi/(mrcsq*tidum**3))
     .                       *  EXP (-3.0*xsi)
c        call STOP ('function DDNRATE: TI out of range', 28)
      end if
      return
c
      end

      real*8 function ddprate (ti)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ------------------------------------------------- 11/20/95 --- HSJ ---
c
c             returns rate(cm**3/sec) of d(d,p)t reaction
c             (the division by 2 for collisions between like
c              particles is NOT accounted for here).
c             New Bosch & Hale rate coefficient:
c             Bosch & Hale, Nuc. Fus., vol32, no.4 (1992) 611
c
      real*8 mrcsq
c
c     data for d(d,p)t:
c
      data  C1,C2,C3,C5,B_gsq, mrcsq
     .    / 5.65718e-12,  3.41267e-3, 1.99167e-3,
     .      1.05060e-05, 985.7715, 937814.0 /
c
      if (ti .lt. 0.2) then
        ddprate = 0.0
      else if (ti .le. 100.0) then
        theta   = ti*C2/(1.0+ti*(C3+ti*C5))
        theta   = ti/(1.0-theta)
        xsi     = (B_gsq/(4.0*theta))**(0.33333333334)
        ddprate = C1 * theta * SQRT (xsi/(mrcsq*ti**3))
     .                       *  EXP (-3.0*xsi)
      else
        tidum = 100.
        theta   = tidum*C2/(1.0+tidum*(C3+tidum*C5))
        theta   = tidum/(1.0-theta)
        xsi     = (B_gsq/(4.0*theta))**(0.33333333334)
        ddprate = C1 * theta * SQRT (xsi/(mrcsq*tidum**3))
     .                       *  EXP (-3.0*xsi)
c        call STOP ('function DDPRATE: TI out of range', 29)
      end if
      return
c
      end

      subroutine delrep (u, usave, pre, delit, theta, relaxc, ier_type)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c     this subroutine, which is called by CHECK_SOLUTION, does the
c     repetitive calculations required to obtain delit and to update
c     the predictor
c ----------------------------------------------------------------------
c
      cor   = theta*u + (1.0-theta)*usave
      cor1  = cor
      if (     cor .eq. 0.0)  go to 10
      if (ier_type .gt. 0  )  cor1 = 1.0 ! use abs error if ier_type > 0
      del   = ABS (pre-cor) / cor1
      delit = MAX (delit, del)
***10 pre   = cor
   10 pre   = (1.0-relaxc)*pre + relaxc*cor    ! HSJ relax solution
      return
c
      end

      subroutine delta_sol (u, usave, nk, nj, kk, r, delav, delmax,
     .             kmax, jmax, iangrot,te_index,ti_index,
     .             rot_index,itran,analysis_check)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
      integer te_index,ti_index,rot_index,itran(*),
     .                                analysis_check
c
c ----------------------------------------------------------------------
c     this subroutine calculates delav and delmax, the average and maximum
c     relative change of u from one time step to the next.
c     however, contributions from ABS (rpb) .near. 0 are omitted.
c     (near zero crossings, or if rbp is less than 5% of the value
c     given by a parabolic current profile).
c ----------------------------------------------------------------------
c
      dimension  u(kk,*), usave(kk,*), r(*), ischng(3)
c
c     check for up to 3 sign changes in rbp
c
      ischng(1) = 0
      ischng(2) = 0
      ischng(3) = 0
      isc       = 1
      do 5 j=2,nj-1
        if (usave(nk-iangrot,j)*usave(nk-iangrot,j+1) .gt. 0.0)  go to 5
        ischng(isc) = j
        isc         = isc + 1
        if (isc .gt. 3)  go to 6
    5 continue
c
    6 delav  = 0.0
      delmax = 0.0
      do 10 j=1,nj
          rnorm2 = (r(j)/r(nj))**2
c
c ------------------------------------------------------ SCC/HSJ 2/21/92
c
****      do 10 k=1,nk-iangrot
          do 10 k=1,nk
c
             if( analysis_check .eq. 0 .and. itran(k) .eq. 0)go to 10
             if (usave(k,j) .eq. 0.0       )  go to 10
             !for fixed edge te, ti cases:
             if( j .ge. te_index .and. k .eq. nk-iangrot-2 )go to 10
             if( j .ge. ti_index  .and. k .eq. nk-iangrot-1 )go to 10
             if( j .ge. rot_index  .and. k .eq. nk .and. iangrot
     .                          .eq. 1)go to 10
             if (         k .ne. nk-iangrot)  go to 21
                do is=1,3
                  if (ischng(is) .eq. 0)  go to 20
                  if (j .ge. (ischng(is)-1) .and. j .le. 
     .                        (ischng(is)+2))  go to 10
                end do
c
   20           if (ABS (usave(k,j)) .lt. 0.05*usave(k,nj)*
     .                      rnorm2*(2.0-rnorm2)) go to 10
   21        del    = ABS ((u(k,j)-usave(k,j)) / usave(k,j))
             delav  = delav + del
             if (del .le. delmax)  go to 10
             delmax = del
             kmax   = k
             jmax   = j
   10     continue
      xkj    = nk*nj
      delav  = delav/xkj
      k      = kmax
      j      = jmax
c
c ------------------------------------------------------ SCC/HSJ 2/21/92
c
      if (delmax .eq. 0.0) then
        jmax = 1
        kmax = 1
      else
        delmax = ABS ((u(k,j)-usave(k,j)) / usave(k,j))
      end if
c
      return
c
      end


      subroutine dif2ydx(dr, a, a2prime, nj)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c --- given a profile,a(j),j = 1,2..nj,as a function of r(j),get the 
c --- second derivative on the same grid
c --- non uniform grid is allowed
c ----------------------------------------------------------------------
c
      dimension  a(*), a2prime(*), dr(*)
      do k=2,nj-1
         hf = dr(k)*0.5
         hb = dr(k-1)*0.5
         hp = hb+hf
        a2prime(k) = 2.d0*(a(k-1)/(hb*hp) -(a(k)/(hf*hb)) +
     .               a(k+1)/(hp*hf))
      end do

         !forward difference
        a2prime(1) = 2.d0*(a(1)/(hb*hp) -(a(2)/(hf*hb)) +
     .               a(3)/(hp*hf))
        !backward  difference
        a2prime(nj) = 2.d0*(a(nj-2)/(hb*hp) -(a(nj-1)/(hf*hb)) +
     .               a(nj)/(hp*hf))

      return
c
      end

      subroutine difydxhalfgrid (dr, a, aprime, nj)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c --- given a profile,a(j),j = 1,2..nj,as a function of r(j),get the derivative
c --- on the half grid points k = j+1/2,k=1,2..nj-1:
c --- non uniform grid is allowed
c ----------------------------------------------------------------------
c
      dimension  a(*), aprime(*), dr(*)
      do k=1,nj-1
        aprime(k) = (a(k+1)-a(k))/dr(k)
      end do
      return
c
      end



      subroutine divflx (flux, fluxb, hcap, r, ra, drr, k, nj,
     .                   kkloc, dflux)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c     this subroutine calculates the divergence of the flux
c     (1./H*rho)*d/drho(H* rho* gamma) for the k'th variable
c     Note that the flux is assumed to be specified on the half grid
c ----------------------------------------------------------------------
c
      dimension flux(kkloc,*), fluxb(*), hcap(*), r(*), ra(*), drr(*),
     .          dflux(*)
c
      j = 1
      dflux(j) = drr(j)*flux(k,j)*(hcap(j+1)+hcap(j))/hcap(j)
      do 10 j=2,nj-1
        havm = 0.5*(hcap(j)+hcap(j-1))
        havp = 0.5*(hcap(j)+hcap(j+1))
   10   dflux(j) = drr(j)*(ra(j)*havp*flux(k,j)-
     .           ra(j-1)*havm*flux(k,j-1))/(r(j)*hcap(j))
      j = nj
      dflux(j) = drr(j)*(r(j)*hcap(j)*fluxb(k)-
     .         ra(j-1)*havp*flux(k,j-1))/(r(j)*hcap(j))
      return
c
      end

      subroutine dndt (dtinv, ene, enesav, kj, kk, nj, nion, u, usave,
     .                 dnedt, dnidt,dtedt,dtidt,dpedt,dpidt,iangrot,
     .                 nk, dangrot)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c     This subroutine calculates dnedt, dnidt and dangrot/dt:
c ----------------------------------------------------------------------
c
      dimension ene(*), enesav(*), u(kk,*), usave(kk,*)
      dimension dnedt(*), dnidt(kj,*),dangrot(*),dpedt(*),dpidt(kj,*),
     .          dtedt(*),dtidt(*)
c
      DO j=1,nj
        dnedt(j) = (ene(j)-enesav(j))*dtinv
        dtedt(j) = (u(nion+1,j)-usave(nion+1,j))*dtinv
c        dpedt(j) = 0.75_DP*((ene(j)+enesav(j))*dtedt(j) +
c     .              (u(nion+1,j)+usave(nion+1,j))* dnedt(j))
        dpedt(j) = 1.5*(ene(j)*u(nion+1,j) - enesav(j)*
     .                     usave(nion+1,j))*dtinv
        dtidt(j) = (u(nion+2,j)-usave(nion+2,j))*dtinv
        if (iangrot .ne. 0) dangrot(j) = (u(nk,j)-usave(nk,j))*dtinv
        DO  k=1,nion 
           dnidt(j,k) = (u(k,j)-usave(k,j))*dtinv
c           dpidt(j,k) = 0.75*(dnidt(j,k)*(u(nion+2,j)+usave(nion+2,j))+
c     .                        dtidt(j)*(u(k,j)+usave(k,j)))
           dpidt(j,k)  = 1.5*(u(nion+2,j)*u(k,j) - 
     .                   usave(nion+2,j)*usave(k,j))*dtinv
        ENDDO
      ENDDO
      return
c
      end

      subroutine dtlim (iter, delit, relit, itmax, idterr)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c     DTLIM helps keep the time step to a numerically stable regime.
c     The following condition, if met, indicates dt (dtmax) is too big:
c     -limiting enabled (ilimdt = 1): assumed if dtlim is called.
c     -about to halve dt (iter = itmax)
c     -in this step, 2 attempts to achieve particle
c        conservation failed (icount .ge. 2):
c        this is symptomatic of instability
c ----------------------------------------------------------------------
c
      data icount / 0 / 
      idterr = 0
      if (iter .ne. 1)  go to 100
c
c          First iteration, reset and return
c
      icount = 0
      return
  100 if (delit .gt. relit)  go to 110
      icount = icount + 1
  110 if (iter .ne. itmax)  return
c
c          Last iteration: if not converged, check icount for
c          number of attempts at particle conservation.
c
      if (icount .lt.     2)  return
      if (delit  .le. relit)  return
      idterr = 1
      return
c
      end

      real*8 function dtrate (ti)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ------------------------------------------------- 11/20/95 --- HSJ ---
c     returns rate(cm**3/sec) of t(d,n)he4 reaction
c     New Bosch & Hale rate coefficient:
c     Bosch & Hale, Nuc. Fus., vol32, no.4 (1992) 611
c ----------------------------------------------------------------------
c
c     data for t(d,n)he4:
c
      data  C1,C2,C3,C4,C5,C6,C7, B_gsq, mrcsq
     .    / 1.17302e-9, 1.51361e-2, 7.51886e-2, 4.60643e-3,1.3500e-2,
     .     -1.06750e-4, 1.36600e-5,   1.182170e+3, 1124656.0 /
c
c     neutrons produced by bulk plasma d-t fusion:
c
      theta  = ti*(C2 +ti*(C4+ti*C6))
      theta  = theta/(1.0+ti*(C3+ti*(C5+ti*C7)))
      theta  = ti/(1.0-theta)
      xsi    = (B_gsq/(4.0*theta))**(0.33333333334)
      dtrate = C1 * theta * SQRT (xsi/(mrcsq*ti**3)) * EXP (-3.0*xsi)
      return
c
****  function dtrate (tempk)
****
****  undocumented - will not use HSJ   11/20/95
****  data  a/-0.768165/, b/4.05877/, c/-46.6242/
****
****  y = 5.5e-21
****  if (tempk .le. 1.0)  go to 2000
****  x = LOG (tempk)
****  y = (a*x+b)*x+c
****  y = EXP (y) * (tempk*tempk)
* 2000 dtrate = y
c
      end

      subroutine dumper
c
c --- Dump source and solution data to file 'solution' for debugging.
c --- This subroutine is only called by subroutine SOLVE (right now).
c --- Written by Holger St. John and Joe Freeman on 06 April 1999 to
c --- aid in debugging the port of the latest version of ONETWO from
c --- HP (HP-UX) to DEC (Digital UNIX); DEC is producing bad numbers
c --- with Masanori Murakami's test cases; HP and CRAY output are ok.
c
      USE param
      USE io
      USE soln
      USE numbrs
      USE sourc
      USE geom
      USE flags
      USE tordlrot
      USE bd_condtn,only : bctime,ub,fluxb,
     .    ub_save,ub_rho_edge,bctime_zone
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c      include 'param.i'
c      include 'bcon.i'
      include 'imsl.i'
c      include 'flags.i'
c      include 'numbrs.i'
c      include 'soln.i'
c      include 'tordlrot.i'
c      include 'io.i'
c      include 'geom.i'
c      include 'sourc.i'
c
      integer unit
      call getioun(unit,54)
c
      open  (unit = unit, file = 'solution', status = 'UNKNOWN')
c
      write (unit = unit, fmt = '(a // a / (6e13.5))')
     .      '@@@@ BEGIN @@@@',
     .      'FIRST write: (kk, kj) REAL arrays',
     .     ((u(kkv,kjv), usave(kkv,kjv), uav(kkv,kjv), uav0(kkv,kjv),
     .       dudtsv(kkv,kjv), s(kkv,kjv),
     .       kkv=1,kk), ! first  index
     .       kjv=1,kj)  ! second index
c
      write (unit = unit, fmt = '(/ a / (6e13.5))')
     .      'SECOND write: some (kj) REAL arrays',
     .      (te(kjv), ti(kjv), rbp(kjv), ene(kjv), enesav(kjv),
     .       curden(kjv), etor(kjv), eneav0(kjv),
     .       eneav1(kjv), curpar_soln(kjv), cur_tor_ps_soln(kjv),
     .       fcap0(kjv), dfdt(kjv), fcap(kjv), gcap0(kjv), dgdt(kjv),
     .       gcap(kjv), hcap0(kjv), sione(kjv), sfus(kjv), stfus(kjv),
     .       sbfus(kjv), sbeam(kjv), qdimpl(kjv), qdelt(kjv),
     .       qgam(kjv), qohm(kjv), qione(kjv), qioni(kjv),
     .       qcx(kjv), qbeame(kjv), qbeami(kjv), qrfe(kjv), qrfi(kjv),
     .       qrad(kjv), qfus(kjv), qfuse(kjv), qfusi(kjv),
     .       qbfus(kjv), qbfuse(kjv), qbfusi(kjv), qtfus(kjv),
     .       qtfuse(kjv), qtfusi(kjv), qmag(kjv), qsawe(kjv),
     .       qsawi(kjv), ssawe(kjv), qexch(kjv), sum2d(kjv), qe2d(kjv),
     .       qi2d(kjv), curohm(kjv), curbe(kjv), curbi(kjv), currf(kjv),
     .       curdri(kjv), eta(kjv), etap(kjv), esaw(kjv),
     .       wejcm3(kjv), wijcm3(kjv), dpedtc(kjv), dpidtc(kjv),
     .       wbjcm3(kjv), pohm1(kjv), pohm2(kjv), cur2(kjv),
     .       kjv=1,kj) ! only index
c
      write (unit = unit, fmt = '(/ a / (8(6e13.5/), 3e13.5/))')
     .      'THIRD write: more (kj) REAL arrays (TROUBLE)',
     .      ( conde(kjv), kjv=1,kj),
     .      ( conve(kjv), kjv=1,kj),
     .      ( condi(kjv), kjv=1,kj), ! trouble
     .      ( convi(kjv), kjv=1,kj),
     .      (qconde(kjv), kjv=1,kj),
     .      (qconve(kjv), kjv=1,kj),
     .      (qcondi(kjv), kjv=1,kj), ! trouble
     .      (qconvi(kjv), kjv=1,kj),
     .      (pconde(kjv), kjv=1,kj),
     .      (pcondi(kjv), kjv=1,kj),
     .      (pconve(kjv), kjv=1,kj),
     .      (pconvi(kjv), kjv=1,kj) ! only index
c
      write (unit = unit, fmt = '(/ a / (6e13.5))')
     .      'FOURTH write: even more (kj) REAL arrays',
     .      (curboot(kjv), qohmi(kjv), curbet(kjv),
     .       curb(kjv), scurdri(kjv), scurfast(kjv),
     .       xjbne(kjv), xjbte(kjv), xjbnf(kjv), gconde(kjv),
     .       gcondi(kjv),curdbeam(kjv), gconve(kjv), gconvi(kjv),
     .       gconvde(kjv), gconvdi(kjv), sbfsav(kjv), qbfsav(kjv),
     .       dhdt(kjv), hcap(kjv), hcapctr(kjv), fcapctr(kjv),
     .       dhdtsv(kjv), dfdtsv(kjv), dhdtcap(kjv), dfdtcap(kjv),
     .       r2cap0(kjv), dr2dt(kjv), r2cap(kjv), rdif, ddifdt,
     .       rhosp0, xi(kjv), r2capi(kjv), r2capi0(kjv), dr2idt(kjv),
     .       rcap(kjv), bsqncap(kjv), rcap0(kjv), bsq_avg_cap(kjv),
     .       drcapdt(kjv),
     .       kjv=1,kj) ! only index
c
      write (unit = unit, fmt = '(/ a / (6e13.5))')
     .      'FIFTH write: (kj, kion) REAL arrays',
     .     ((en(kjv,kionv), s2d(kjv, kion), xjbni(kjv, kion),
     .       xjbti(kjv, kion),
     .       kjv=1,kj),    ! first  index
     .       kionv=1,kion) ! second index
c
      write (unit = unit, fmt = '(/ a / (6e13.5))')
     .      'SIXTH write: (kj, 2) REAL arrays',
     .     ((siadd(kjv,itwo), sion(kjv,itwo), srecom(kjv,itwo),
     .       scx(kjv,itwo),
     .       kjv=1,kj), ! first  index
     .       itwo=1,2)  ! second index
c
      write (unit = unit, fmt = '(/ a / (6e13.5))')
     .      'SEVENTH write: (2) REAL arrays',
     .      (flxmod(itwo), tnion0(itwo), snadd(itwo), snaddt(itwo),
     .       sngas(itwo), flxadd(itwo),
     .       itwo=1,2) ! only index
c
      write (unit = unit, fmt = '(/ a / (6e13.5))')
     .      'EIGHTH write: REAL scalars',
     .       timcap, rhoa0, drhoadt_geom, rhoa, rhoa_save, sfarea,
     .       cxarea, volfac, ali, psivloop, psivlop0, dpsivlop, pvbar,
     .       bpsqsurf, betapmhd, betatmhd, totohm, totbe, totbi, totb,
     .       sn2d, totrf, totboot, xlica, totdri, totbeam, wdelt, wgam,
     .       wohm, refrad
c
      write (unit = unit, fmt = '(/ a / 5i8)')
     .      'NINTH write: INTEGER scalars',
     .       itorfluse, njcur, curtype, diffeq_methd, comp_methd_eqdsk
c
      write (unit = unit, fmt = '(/ a / 3a // a)')
     .      'TENTH write: CHARACTER scalars',
     .       codeid, machinei, eqgrdsze,
     .      '@@@@ END @@@@'
c
**** .       njqin(krf), qine(kjv, krf), qini(kjv, krf), ! overlooked
c
      call giveupus(unit)
      close (unit = unit)
      call STOP ('subroutine DUMPER: hard stop for debugging', 271)
c
      end

      real*8 function eira (te, atw)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c     EIRA computes the electron ionization rate coefficient
c     <sigma-v> for any primary ion: hydrogenic and helium ions
c     have to be handled separately.
c ----------------------------------------------------------------------
c
      if (atw .gt. 3.0)  go to 10
c
c     hydrogenic
c
      eira = eir(te)
      return
c
c     helium
c
   10 eira = eirhe(te)
      return
c
      end

      real*8 function eir (x)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c this function calculates the ionization rate of atomic hydrogen by ele
c impact in units of cm**3/s.  x is in units of keV
c the formula is taken from the paper by r.l. freeman and e.m. jones
c clm-r 137 culham laboratory 1974
c ----------------------------------------------------------------------
c
      if (x .lt. .001 .or. x .gt. 100.0)  go to 20
      ta = LOG (x)+6.9077553
      ta = -0.3173850e2+ta*(0.1143818e2-ta*(0.3833998e1-ta*(0.7046692-ta
     .    *(0.7431486e-1-ta*(0.4153749e-2-ta*0.9486967e-4)))))
      ta = EXP (ta)
      eir = ta
      return
   20 eir = 0.0
      return
c
      end

      real*8 function eirhe (x)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c this function calculates the ionization rate of atomic helium by elect
c impact in units of cm**3/s.  x is in units of keV
c the formula is taken from the paper by r.l. freeman and e.m. jones
c clm-r 137 culham laboratory 1974
c ----------------------------------------------------------------------
c
      if (x .lt.  0.025)  xx = LOG (0.025)+6.9077553
      if (x .gt. 10.0  )  xx = LOG (10.0)+6.9077553
      if (x .ge. .025 .and. x .le. 10.0) xx = LOG (x)+6.9077553
      xx = -0.4450917e2 + xx*(0.2442988e2 + xx*(-0.1025714e2
     .     + xx*(0.2470931e1 + xx*(-0.3426362 + xx*(0.2505100e-1
     .     - xx*0.7438675e-3)))))
      eirhe = EXP (xx)
      return
c
      end

      subroutine elongt (time, fkap, dfkapdt)
c
      USE param
      USE numbrs
      USE machin
      USE bd_condtn,only : bc,bctime,ub,fluxb,
     .    ub_save,ub_rho_edge,bctime_zone
      USE tordlrot,only: iangrot
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c      include 'param.i'
c      include 'bcon.i'
c      include 'numbrs.i'
c      include 'machin.i'
c
      call find (i1, i2, time, bctime, nbctim)
      fkap    = 0.0
      dfkapdt = 0.0
      if (i1 .eq. 0)
     .  call STOP ('subroutine ELONGT: i1 = 0', 82)
      if (elong(1) .lt. 0.0)  go to 50
      if (elong(2) .eq. 0.0) then
        fkap = elong(1)
        return
      end if
c
c     the following give kappa and d(kappa)/dt as a function of time,
c     as specified by the elong(m) array.
c
      fkap    = elong(i1)
      dfkapdt = 0.0
      if (i1 .eq. i2)  return
      dfkapdt = (elong(i2)-elong(i1))/(bctime(i2)-bctime(i1))
      fkap    =  elong(i1)+(time-bctime(i1))*dfkapdt
      return
c
c     The following gives kappa(t) .ge. 1. such that as the current is
c     ramped in time, the safety factor at the plasma edge is a constant
c     value equal to -elong(1).
c
c     the input total current data has been saved in bc(m,nk-iangrot) = 0.2*totcur(m)
c
   50 pcur    = 5.0*bc(i1,nk-iangrot)
      dpcurdt = 0.0
      if (i1 .eq. i2)  go to 60
      dpcurdt = 5.0*(bc(i2,nk-iangrot)-bc(i1,nk-iangrot)) /
     .                                (bctime(i2)-bctime(i1))
      pcur    = pcur+(time-bctime(i1))*dpcurdt
   60 c1      = 0.4*rmajor*(-elong(1))/(rminor*rminor*btor)
      fkap2   = c1*pcur-1.0
      fkap    = SQRT (fkap2)
      dfkapdt = c1*dpcurdt/(2.0*fkap)
      if (fkap .gt. 1.0)  return
      fkap    = 1.0
      dfkapdt = 0.0
      return
c
      end

      subroutine expand (u, usave, itran, nk, nkt, nj, kk)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c     this subroutine expands the order of the solution vector u from nkt to nk
c
      dimension  u(kk,*), usave(kk,*), itran(*)
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
        end do
   40 continue
      return
c
      end

      subroutine extrap (x1, x2, x, y1, y2, y)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c --- calculates y at x from (x1,y1) and (x2,y2) using linear extrapolation
c
      xx = (x-x1)/(x2-x1)
      y  = y1 + xx*(y2-y1)
      return
c
      end

      subroutine fhextrap (dt, fcap0, fcap, hcap0, hcap,
     .                     dfdtsv, dhdtsv, nj, implicit_fh)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c   extrapolate fcap and hcap to the central time point
c   as appropriate, using the derivative from the previous time step.
c --- input
c   dt                   time interval,sec
c   fcap0                value of fcap,gcap,hcap, at time t
c   hcap0
c   dfdtsv
c   dhdtsv               corresponding derivatives
c   nj                   r grid size
c   implicit_fh          logical,if true do the extrapolation
c                        for fcap and hcap
c
c --- output:
c    fcap                 at time t+dt
c    hcap                 at time t+dt  (if implicit_fh = true)
c
c ------------------------------------------------------------------ HSJ
c
      dimension fcap0(*),fcap(*),dfdtsv(*),
     .          hcap0(*),hcap(*),dhdtsv(*)
      logical   implicit_fh
c
      if (.not. implicit_fh)  return
      do j=1,nj
        fcap(j) = fcap0(j)+dfdtsv(j)*dt
        hcap(j) = hcap0(j)+dhdtsv(j)*dt
      end do
      return
c
      end

      subroutine fhrecall
c ----------------------------------------------------------------------
c  subroutine sets fcap = fcap0,hcap=hcap0
c  and dfdt = dfdtsv,dhdt=dhdtsv
c  input
c  INCLUDE file param:
c    kj                   dimension of cap parameters
c
c  INCLUDE file geom.i:
c    fcap0(j)
c    hcap0(j)
c    dhdtsv(j)
c    dfdtsv(j)
c
c  INCLUDE file numbrs:
c    nj                   actual size of cap parameters
c
c  output:
c  INCLUDE file geom.i:
c    fcap(j)
c    hcap(j)
c    dhdt(j)
c    dgdt(j)
c    dfdt(j)
c
c ------------------------------------------------------------------ HSJ
c

      USE param
      USE numbrs
      USE geom
      USE soln2d
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c      include 'param.i'
c      include 'geom.i'
c      include 'numbrs.i'
c      include 'soln2d.i'
c
      if (.not. implicit_fh)  return
      do j=1,nj
        fcap(j) = fcap0(j)
        dfdt(j) = dfdtsv(j)
        hcap(j) = hcap0(j)
        dhdt(j) = dhdtsv(j)
      end do
      return
c
      end
      FUNCTION fhsource(f2,dp,dbp,c,fesq)
         fhsource  = -(dp+f2*f2*c*dbp) / (fesq+c*f2)
      END FUNCTION
      FUNCTION dfhsdf(f2,dp,dbp,c,fesq)
         dfhsdf = -c*(f2*f2*c*dbp +2.0*f2*dbp*fesq-dp)/
     .                               ((f2*c+fesq)**2)
      END FUNCTION
      subroutine fhcalc (imks, iuse)
c
c
c ----------------------------------------------------------------------
c subroutine calculates the profiles fcap, hcap,using the
c the flux surface average of Grad-Shafranov equation.
c
c  input ( <a> means flux surface average of a):
c
c  argument list:
c  imks                   imks = 0 means calculate pressure in CGS unit
c                              = 1                             MKS
c  iuse                   iuse = 1 means use en,ene,te,ti for pressure
c                              = 2           u,ene,wbeam,walp
c                              = 3           usave,enesav,wbeam,walp
c
c  INCLUDE file constnts.i:
c    pisq                     pi*pi
c
c  INCLUDE file geom.i:
c    r2cap(j)                 j = 1,2..nj =<rmajor**2/r**2>
c
c  INCLUDE file machin.i:
c    rmajor                   R0,cm
c    btor                     toroidal b field,gauss
c
c  INCLUDE file mesh.i:
c    r(j)                     j = 1,2..nj current rho grid
c    dr(j)                    j = 1,2..nj-1,dr(j) = r(j+1)-r(j)
c
c  INCLUDE file soln.i:
c    u(kk,kj)                 solution vector (needed if iuse = 2)
c
c  INCLUDE file soln2d.i:
c    implicit_fh              if true calculate fh implicitly
c    rbp(j)                   j = 1,2,..nj,=fcap*gcap*hcap*rho*bp
c                             in units of gauss cm (needed if iuse = 1)
c
c  INCLUDE file numbrs:
c    nj                       transport grid size
c
c  INCLUDE file rhog.i:
c    press(j)                 j = 1,2..nj,the pressure,including beam
c                             and alpha particles due to fusion as
c                             well as impurity ions.
c
c  INCLUDE file storage.i:
c    xdum(j)                 j = 1,2..nj,xdum,ydum are temporary storage
c    ydum(j)                 vectors of min length nj.they contain no
c    zdum(j)                 useful information on output.
c    wdum(j)
c    vdum(j)
c
c  INCLUDE file tordlrot.i:
c    iangrot                   needed to get index into u
c                              iangrot = 1 if rotation is included
c                              iangrot = 0 if rotation is excluded
c
c  output:
c
c  INCLUDE file geom.i:
c    fcap
c    hcap
c    sfarea                 "surface area" (without grad rho) cm**2
c
c ------------------------------------------------------------------ HSJ
c
      USE param
      USE soln
      USE numbrs
      USE mesh
      USE machin
      USE geom
      USE tordlrot
      USE constnts
      USE soln2d
      USE rhog
      implicit  integer (i-n), real*8 (a-h, o-z)
c      include 'param.i'
c      include 'constnts.i'
c      include 'geom.i'
c      include 'machin.i'
c      include 'mesh.i'
c      include 'numbrs.i'
c      include 'rhog.i'
c      include 'soln.i'
c      include 'soln2d.i'
c      include 'storage.i'
c      include 'tordlrot.i'
c
      dimension       fsq(kj),                       ! work vector
     .             rbplcl(kj),                       ! work vector
     .               cvec(kj),                       ! work vector
     .            dbpdrho(kj),                       ! work vector
     .             dpdrho(kj),                       ! work vector
     .          fcap_save(kj),                       ! work vector
     .          hcap_save(kj)                        ! work vector
c      equivalence (    fsq(1), xdum(1))              ! temporary storage
c      equivalence ( dpdrho(1), ydum(1))              ! temporary storage
c      equivalence ( rbplcl(1), zdum(1))              ! temporary storage
c      equivalence (   cvec(1), wdum(1))              ! temporary storage
c      equivalence (dbpdrho(1), vdum(1))              ! temporary storage
c      equivalence (fcap_save(1), vdum(kj+1))         ! temporary storage
c      equivalence (hcap_save(1), vdum(3*kj+1))       ! temporary storage

c


      data         dfhmax, iterfhmax /1.0e-06, 20/
c
      if (.not. implicit_fh)  return
      call copya(fcap,fcap_save,nj)
      call copya(hcap,hcap_save,nj)    ! save in case of failure
c
c --- set up the constants that don't change during the iteration:
c --- fsq has units of (gauss-cm)**2):
c
      fedgesq = rmajor*rmajor*btor*btor
      r0sq    = rmajor*rmajor    ! cm**2
      const1  = 8.0 * pi*r0sq
c
c --- get the appropriate poloidal b field:
c
      if (iuse .eq. 1) then
          call copya(rbp,rbplcl,nj)
      else if (iuse .eq. 2) then
          do j=1,nj
             rbplcl(j) = u(nk-iangrot,j)
          end do
      else if (iuse .eq. 3) then
          do j=1,nj
             rbplcl(j) = usave(nk-iangrot,j)
          end do
      end if
c
c --- square the field and differentiate it:
c
      do j=1,nj
        rbplcl(j) = rbplcl(j)**2
      end do
      call difydx(r,rbplcl,dbpdrho,nj)
c
c --- form the required functions of bp**2:
c
      do j=2,nj
          cvec(j)    = rbplcl(j)*r0sq/(gcap(j)*r(j)*r(j))
          dbpdrho(j) = dbpdrho(j)/rbplcl(j)
      end do
      cvec(1)    = 0.0       ! (not used)
      dbpdrho(1) = 0.0       ! (not used)
c
c --- get the pressure,press, on the rho grid(units are dyne/cm**2).
c --- (results are returned by INCLUDE file rhog.i):
c
      call pressr(imks,iuse)
c
c --- next get the gradient,dpdrho:
c
      call difydx(r,press,dpdrho,nj)
      dpdrho(1) = 0.0
c
c --- form the term needed in the flux avg g.s. equation that
c --- involves the pressure gradient:
c
      do j=2,nj
        dpdrho(j) = dpdrho(j)*const1/r2cap(j)
      end do
c
c --- set the boundary condition at r(nj):
c
      fsq(nj) = 1.0
c
c ----------------------------------------------------------------------
c --- do the integration from the plasma edge inward to the axis.
c --- we solve dfsq/dr = s(r,fsq) (which is nonlinear due
c --- to the dependence of the source s on fsq) using forward
c --- (explicit)euler method on predictor step and backward
c --- (implicit) euler method on correctors steps. the corrector
c --- is converged using Newton's method.
c ----------------------------------------------------------------------
c
      jstep = nj
c
c --- evaluate the source at the plasma edge:
c
      fhsrc = fhsource(fsq(nj),dpdrho(nj),dbpdrho(nj),cvec(nj),fedgesq)
c
c --- predictor step is explicit
c
  100 continue
      iterfh  = 0    ! counts number of Newton iteration for corrector
      itry    = 0    ! allows one try to fix Newton's method if it fails
      fsqpred = fsq(jstep) - dr(jstep-1)*fhsrc     ! negative..
c                                                    ..because integrate inwards
c --- corrector is implicit and iterated
c
  120 iterfh = iterfh+1
      fhsrc  = fhsource(fsqpred,dpdrho(jstep-1),dbpdrho(jstep-1),
     .                            cvec(jstep-1),fedgesq)
      anum   = fsqpred - fsq(jstep)+dr(jstep-1)*fhsrc
      denom  = 1.0 + dr(jstep-1)*dfhsdf(fsqpred,dpdrho(jstep-1),
     .                         dbpdrho(jstep-1),cvec(jstep-1),fedgesq)
      dfsq = -anum/denom
      if (iterfh .eq. 1) then
          dfsqsv = dfsq         ! save for use if Newton fails
          anumsv = anum
      end if
      fsqcor = fsqpred+dfsq
      dfhrel = ABS (dfsq/fsqcor)
      if (dfhrel .lt. dfhmax)  go to 150    ! Newton iterates converged
      if (iterfh .gt. iterfhmax) then       ! Newton didn't converge
          if (itry .lt. 1) then
c
c             the initial Newton step was in the right direction
c             but its magnitude was so large that it carried us out
c             of the region of convergence. use this fact to bracket solution.
c
              itry  = itry + 1
              stepl = 0.0
              do while (stepl .lt. 1.0)
                  stepl = stepl + 0.03
c
c                 step in Newton direction:
c
                  fsqtry = fsq(jstep) + stepl*dfsqsv
                  fhsrc1 = fhsource(fsq(jstep),dpdrho(jstep),
     .                         dbpdrho(jstep),cvec(jstep),fedgesq)
                  anum   = fsqtry-fsq(jstep)+dr(jstep-1)*fhsrc1
                  if (anum*anumsv .le. 0.0) then    ! solution bracketed
                    iterfh  = 0
                    fsqpred = 0.5*(fsqtry+fsq(jstep)+(stepl-0.03)*dfsq)
                    go to 120
                  end if
              end do
              go to 1000   ! method failed, quit
          else
c
c            give up if the first refinement didn't work
c
             go to 1000
          end if
      end if
      fsqpred = fsqcor
      go to 120   ! iterate the corrector
c
c --- converged at rho(jstep). go to next rho, i.e., rho(jstep-1):
c
  150 jstep      = jstep-1
      fsq(jstep) = fsqcor
      if (jstep .gt. 1)  go to 100
c
c --- get fsq(1) by extrapolation. we use d(fsq)/drho = 0 at rho=0
c
      call cubicextrp(fsq(2),fsq(3),fsq(4),r(2),r(3),r(4),fsq(1),2)
c
c --- done with calculation of fsq. get new fcap and hcap
c
      do j=1,nj
        fcap(j) = 1.0 / SQRT (fsq(j))
        hcap(j) = fcap(j) / r2cap(j)
      end do
c
c --- finally get surface area of plasma
c
   10 sfarea = 4.0 * pisq * rmajor * hcap(nj) * r(nj)
c
****  write  (6, 5)                iterfh, iterfhmax, dfhrel
****5 format (' subroutine FHCALC: iterfh, iterfhmax, dfhrel = ',
**** .          2(2x,i4), 1pe12.2)
      return
c
 1000 write  (6, 6)  iterfh, dfhrel
    6 format (/ ' subroutine FHCALC reports:'                  /
     .          '   Newton''s method did not converge in ', i5,
     .            ' iterations'                                /
     .          '   error = ', 1pe12.3                         /
     .          '   ONETWO will be stopped')
c
****  call STOP ('subroutine FHCALC: failure to converge', 83)
****
****  Instead of stopping use the previous solution:
c
      call copya (fcap_save, fcap, nj)
      call copya (hcap_save, hcap, nj)
      go to 10
c
      end

      subroutine fhupdate
c ----------------------------------------------------------------------
c  subroutine sets fcap0 = fcap,hcap0=hcap
c  and dfdtsv = dfdt,dhdtsv=dhdt
c  input
c  INCLUDE file param:
c    kj                   dimension of cap parameters
c
c  INCLUDE file geom.i:
c    fcap(j)
c    hcap(j)
c    dhdt(j)
c    dfdt(j)
c
c  INCLUDE file numbrs:
c    nj                   actual size of cap parameters
c
c  INCLUDE file soln2d.i:
c    implicit_fh          if true f,h,are calculated implicitly
c
c  output:
c  INCLUDE file geom.i:
c    fcap0(j)
c    hcap0(j)
c    dhdtsv(j)
c    dfdtsv(j)
c
c ------------------------------------------------------------------ HSJ
c
c
      USE param
      USE numbrs
      USE geom
      USE soln2d
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c      include 'param.i'
c      include 'geom.i'
c      include 'numbrs.i'
c      include 'soln2d.i'
c
      if (implicit_fh) then
        do j=1,nj
          fcap0 (j) = fcap(j)
          dfdtsv(j) = dfdt(j)
          hcap0 (j) = hcap(j)
          dhdtsv(j) = dhdt(j)
        end do
      end if
      return
c
      end

      subroutine find (m, n, val, array, ia)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c     this routine finds m an n such that array(m) > val > array(n)
c     if val  =  array(i),m=n=i are returned
c     if val is not within the array limits, m = n=0 are returned
c     the routine presumes that array(i+1) .gt. array(i)
c ----------------------------------------------------------------------
c
      dimension array(*)
c
      if (val-array(1 ))  12, 15, 10
   10 if (val-array(ia))  30, 20, 12
   12 m = 0
      n = 0
      return
   15 m = 1
      n = 1
      return
   20 m = ia
      n = ia
      return
   30 m = 1
      n = ia
   35 im = (m+n)/2
      if (im .eq. m)  return
      if (val-array(im)) 40,50,45
   40 n = im
      go to 35
   45 m = im
      go to 35
   50 m = im
      n = im
      return
c
      end

      subroutine flxcal (dflux, drr, hcap, nj, r, ra, flux)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c  This subroutine calculates the flux given its divergence.
c --- given dflux(j) ,j = 1,2..nj-1 on the full grid,get flux(k),k=1,2..nj-1
c --- on the half grid by expanding the expression
c ---    (1.0/hr)d/dr(hrgamma) = dflux about full grid point j
c --- here flux(k) corresponds to full grid point j plus 1/2
c --- i.e., flux(k) is at radial position (r(j+1)+r(j))/2
c --- note that drr is inverse deltar     HSJ  (notes vol 2, pg 32)
c ----------------------------------------------------------------------
c
      dimension  dflux(*), drr(*), hcap(*), r(*), ra(*), flux(*)
c
      k       = 1
      flux(k) = hcap(k)*dflux(k)/((hcap(k)+hcap(k+1))*drr(k))
      do k=2,nj-1
        havm    = 0.5*(hcap(k-1)+hcap(k))
        havp    = 0.5*(hcap(k)+hcap(k+1))
        flux(k) = (hcap(k)*r(k)*dflux(k)+havm*ra(k-1)*drr(k)*flux(k-1))
     .          / (havp*ra(k)*drr(k))
      end do
              
      return
c
      end

      subroutine flxtitle (nout, i, nk)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      if      (nk .eq. 4) then
        write (nout, 54) (i,k , k=1,nk)
      else if (nk .eq. 5) then
        write (nout, 55) (i,k , k=1,nk)
      else if (nk .eq. 6) then
        write (nout, 56) (i,k , k=1,nk)
      else if (nk .eq. 7) then
        write (nout, 57) (i,k , k=1,nk)
      else if (nk .eq. 8) then
        write (nout, 58) (i,k , k=1,nk)
      end if
c
      return
c
   54 format (4x,'j',7x,'r',1x,4(8x,'(',i1,',',i1,')'),7x,'sum' /
     .       10x,'(cm)')
   55 format (4x,'j',7x,'r',1x,5(8x,'(',i1,',',i1,')'),7x,'sum' /
     .       10x,'(cm)')
   56 format (4x,'j',7x,'r',1x,6(8x,'(',i1,',',i1,')'),7x,'sum' /
     .       10x,'(cm)')
   57 format (4x,'j',7x,'r',1x,7(8x,'(',i1,',',i1,')'),7x,'sum' /
     .       10x,'(cm)')
   58 format (4x,'j',7x,'r',1x,8(8x,'(',i1,',',i1,')'),7x,'sum' /
     .       10x,'(cm)')
c
      end

      real*8 function fqq (rr)
c
      USE param
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c      include 'param.i'
      include 'delben.i'
c
      fqq = seval (nben, rr, rben, qben, bq, cq, dq)
      return
c
      end

      real*8 function fqqq (rr)
c
      USE param
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c      include 'param.i'
      include 'delben.i'
c
      fqqq = fqq(rr) - qs
      return
c
      end

      subroutine fred
c
      USE param
      USE io
      USE nub
      USE solcon
      USE soln
      USE numbrs
      USE mesh
      USE sourc
      USE machin
      USE geom
      USE constnts
      USE flx
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c      include 'param.i'
c      include 'nub.i'
c      include 'flx.i'
c      include 'geom.i'
c      include 'io.i'
c      include 'machin.i'
c      include 'mesh.i'
c      include 'numbrs.i'
c      include 'solcon.i'
c      include 'soln.i'
c      include 'sourc.i'
c      include 'constnts.i'
c
      dimension  power(kj), hold(kj), factor(kj), rate(kj)
c
c ----------------------------------------------------------------------
c constp = 17.6e6 * 1.6e-19 / 4.0
c ----------------------------------------------------------------------
c
      data constp /7.04e-13/
c
c ----------------------------------------------------------------------
c intialize some constants
c ----------------------------------------------------------------------
c
      constv = 2.0 * pisq*rmajor
      const0 = 2.0 * constv*1.6e-16
      tpower = 0.0
      t      = n
      call header (nout,time,t)
      write (nout, 8000)
c
      do j=1,nj-1
        rate(j)  = dtrate(ti(j))
        power(j) = constv*constp*en(j,1)*en(j,1)*rate(j)*hcap(j)*
     .                                  (r(j+1)*r(j+1)-r(j)*r(j))
        tpower   = tpower+power(j)
      end do
c
      do 21 j=1,nj
        j1prt = ((j-1)/jprt)*jprt
        if (j1prt .ne. j-1 .and. j .ne. nj)  go to 21
        write (nout,8010) j,r(j),ti(j),en(j,1),rate(j),power(j)
   21 continue
c
 8000 format (10x,'ion temperature, density, dtrate, fusion power ' //
     .         4x,'j',7x,'r',9x,'ti',8x,'en',8x,'dtrate',7x,'power' /
     .        10x,'(cm)',6x,'(keV)',5x,'(1/cm**3)',14x,'(watts)')
 8010 format (1x,i4,f10.2,1p5e12.3)
c
c ----------------------------------------------------------------------
c calculate some total powers
c ----------------------------------------------------------------------
c
      do j=1,nj
        factor(j) = hcap(j)*const0
      end do
c
      call trapv (r,qohm,factor,nj,ohmpwr)
c
      do j=1,nj
        hold(j) = qbeame(j)+qbeami(j)
      end do
c
      call trapv (r,hold,factor,nj,beamdpwr)
c
c ----------------------------------------------------------------------
c now calculate q values for fred marcus
c ----------------------------------------------------------------------
c
      pheat = 0.0
c
      do jb=1,nbeams
        pheat = pheat+pbeam(1,jb)+pbeam(2,jb)+pbeam(3,jb)+ohmpwr
      end do
c
      bigq = tpower/pheat
      q2   = tpower/(ohmpwr+beamdpwr)
      write  (nout, 8030)  tpower, pheat, bigq, q2
 8030 format (// ' fusion power = ',e12.3,3x,'beam+ohmic=',e12.3,
     .             3x,'big q = ',e12.3 // ' q2=',e12.3)
      return
c
      end




        subroutine fusion12(initbm,palpha0, vzalpha)
c---------------------------------------------------------------------------
c          FUSION
c
c calculate particle and energy sources due to thermal and fast ion fusion.
c stfus is thermal fusion rate and is assumed equal to rate at
c which alphas are produced in #/(cm**3-sec)
c qtfus energy production rate due to alpha fusion,in keV/cm**3/sec
c note that both of these quantitites are prompt (i.e., alphas born at
c 3.5 Mev)
c NOTE: ONLY THE D(T,N)HE4 REACTION IS INCLUDED IN PARTICLE BALANCE
c       AND ONLY THE 3.5 MeV ALPHA PARTICLE ENERGY IS INCLUDED IN
c       ENERGY BALANCE. OTHER REACTION RATES MAY BE CALCULATED BUT
c       ARE NOT USED IN THE TRANSPORT MODELLING AT THE PRESENT TIME!!!!
c ------------------------------------------- 11/20/95 --------- HSJ ---
c
        USE param,only:    kj
        USE extra,only:bpol
        USE fusion, only : 
     .           ifus,iddfus,fd,id,it,idt,thermal_thermal_ddntot,
     .           thermal_thermal_ddptot,thermal_thermal_dtntot,
     .           thermal_thermal_tt2ntot,thermal_thermal_he3dptot,
     .           ddnfus,ddpfus,dtnfus,ttnfus,he3dpfus,
     .           beam_thermaltth_df,idhe3fus,ihe,ihei,he3_frac,
     .           beam_thermaldth_tf,tot_th_fuse,no_beam_fusion,
     .           iddfusb_bulk,tausa,fencap,enasav,ffe,wasav,enalp,
     .           walp,tauea,wtifus,ffi,beam_thermal_long_calc,
     .           beam_beam_long_calc,ecritalpha,iaslow,taupa
        USE soln,only:     ti,en,ene,te
        USE verbose, only : fusionvb
        USE geom,only:      volfac,hcap
        USE ions,only:      nameb
        USE mesh,only:      r
        USE neut,only :     ineut
        USE numbrs,only :   nj,nion
        USE sourc,only :    sbfus,siadd,qtfus,stfus,qbfus,qtfuse,
     .                      qtfusi,qfusi,qfuse,qbfsav,sbfsav,qfus,sfus,
     .                      qbfuse,qbfusi
        USE io,only :       xdebug
        USE ions,only :     atw,zsq
        USE machin,only:    rmajor
        USE neut,only :     enn
        USE nub,only :      ibeam,fionx,nbeams
        USE nub2,only:      emzrat
        USE solcon,only:    ifreya_old,ifreya,time,timmax,time0,dt
        USE transp,only :   use_nubeam
        USE bd_condtn,only : totcur
        implicit none
        integer j,initbm
        real *8 dum(nj)
        real *8 beam_fusd,beam_fust,vdum(nj),wdum(nj),palpha0(nj),
     .          vzalpha(nj),rtstcxf,atwf0,ezero0,zsqf0,
     .          dtalpha,ffqt,ffqb

       if (ifus .eq. 0) return
c
****   do j=1,nj
****     if (ifus .eq. 2)  go to 2286
****
****     ifus=1,dt mixture
****
****     stfus(j) = fd*(1.0-fd)*en(j,idt)**2*dtrate(ti(j))
****     if (ineut(idt) .ne. 0)  siadd(j,idt) =
****  .                          siadd(j,idt)-2.0*stfus(j)-2.0*sbfus(j)
****     go to 2288
****
****     ifus=2,d,t,separate (one or the other or both)
****
* 2286   stfus(j) = en(j,id)*en(j,it)*dtrate(ti(j))
****     if (ineut(id) .ne. 0)  siadd(j,id) =
****  .                         siadd(j,id) - stfus(j) - sbfus(j)
****     if (ineut(it) .ne. 0)  siadd(j,it) =
****  .                         siadd(j,it) - stfus(j) - sbfus(j)
* 2288   qtfus(j) = 3.5e3*stfus(j)
****   end do
c
      if(fusionvb .gt. 0)print *,'in fusion12,call thermonuclear_rate'

      call thermonuclear_rate (ddnfus,ddpfus,dtnfus,ttnfus,
     .                         he3dpfus,
     .                         thermal_thermal_ddntot,
     .                         thermal_thermal_ddptot,
     .                         thermal_thermal_dtntot,
     .                         thermal_thermal_tt2ntot,
     .                         thermal_thermal_he3dptot,
     .                         iddfus,
     .                         idhe3fus,ihe,ihei,he3_frac,
     .                         ti,en,kj,nj,
     .                         id,idt,it,fd,volfac,hcap,r)
      if(fusionvb .gt. 0)print *,'in fusion12,done thermonuclear_rate'
      do j=1,nj
c           recall that only dt reactions are considered
            if      (nameb .eq. 'd' ) then           ! beam is deuterium
              beam_fusd = sbfus(j)                   ! loss of thermal t
              beam_fust = 0.0
            else if (nameb .eq. 't' ) then           ! beam is tritium
              beam_fust = sbfus(j)                   ! loss of thermal d
              beam_fusd = 0.0
            else if (nameb .eq. 'dt') then   ! beam is dt mixture
              if (ifus .gt. 0) then          ! here the thermal fluid
                    beam_fust = sbfus(j)     ! is also dt
                    beam_fusd = sbfus(j)
              else if (ifus .eq. -1) then    ! individual thermal..
c                                            ..d and t fluids
c                   loss of thermal t due to fast d:
c
                    beam_fust=beam_thermaltth_df(j,3*nbeams+1)
c
c                   loss of thermal d due to fast t:
c
                    beam_fusd=beam_thermaldth_tf(j,3*nbeams+1)
              end if
            else                        ! no beam, or beam is not d or t
              beam_fusd = 0.0
              beam_fust = 0.0
            end if
c
            if (ifus .eq. 1 .and. idt .ne. 0) then
c
c             one of the thermal species is a dt mixture the
c             rate at which d and t are lost from this mixture is twice
c             the reaction rate, since we lose 2 ions per collision.
c
              if (ineut(idt) .ne. 0 .and. iddfus .eq. 2)
     .          siadd(j,idt) = siadd(j,idt)-2.0*dtnfus(j)-2.0*beam_fusd
            end if
c
            if (idt .eq. 0 .and. ifus .gt. 1) then       ! no dt mixture
c
c             if one of the thermal species is deuterium then the
c             rate at which d is lost is (note that d-d reactions
c             are not included in the particle balance):
c
              if (id .ne. 0) then ! ifus = 2..
c                                 ..we have d or t or both d and t
****          if (ineut(id) .ne. 0 .and. iddfus .eq. 1)
c
c             d in system, no t,beam (if present at all)  must be d
c             but d(d,n)he3 is not included in particle balance
c             so nothing is added
c
              if (ineut(id) .ne. 0 .and. iddfus .eq. 3)
c
c               both d and t in system,beam can be d or t
c               burnout of thermal d is due to thermal (dtnfus)
c               and fast t (beam_fust).
c               d(d,n)he3 is neglected :
c
     .          siadd(j,id) = siadd(j,id) - dtnfus(j) - beam_fust
              end if
c
c             if one of the thermal species is tritium then the
c             rate at which t is lost is:
c
              if (it .ne. 0) then
                 if (ineut(it) .ne. 0 .and. iddfus .eq. 3)
c
c                    both d and t in system,
c                    beam, if present, can be d or t
c                    beam must be d to deplete thermal t:
c
     .               siadd(j,it) = siadd(j,it) - dtnfus(j) - beam_fust
c
c             no d in system, only count depletion due to D(T,N)HE4 reaction
c             at present so following lines for T(T,2N)HE4  inactive
c
****                  if (ineut(it) .ne. 0 .and. iddfus .eq. 4)
**** .                     siadd(j,it) =   ! t in system, no d
**** .                     siadd(j,it) - ttnfus(j) - ??
c
                 end if
            else if (idt .eq. 0 .and. ifus .eq. -1) then  ! dt beam
c
c                    both d and t ions and neutrals are present:
c
                     siadd(j,it) =
     .               siadd(j,it) - dtnfus(j) - beam_fust
                     siadd(j,id) =
     .               siadd(j,id) - dtnfus(j) - beam_fusd
            end if
c
            stfus(j) = dtnfus(j)  ! thermal dt rate,#/(cm**3-sec)
            qtfus(j) = 3.5e3*dtnfus(j)
c
c           total energy released in thermal fusion (includes
c           neutron, proton, alpha energies)
c
            tot_th_fuse(j) = dtnfus(j)*17.6e3+ddnfus(j)*4.03e3
     .                     + ttnfus(j)*11.3e3+he3dpfus(j)*18.3e3
     .                                       +ddpfus(j)*4.03e3
      end do





c
c ----------------------------------------------------------------------
c --- a quick and dirty spreading of qtfus - NOT FOR PERMANENT USE
c --- spreadqtfus returns with a box car averaged qtfus.
c --- spreadqtfus to be replaced with a proper routine in the future HSJ
c ----------------------------------------------------------------------
c
      if(fusionvb .gt. 0 .and. xdebug(1) .gt. 0.5)
     .       print *,'in fusion12,calling spreadqtfus'
      if (ABS (xdebug(1)) .gt. 0.5)
     .  call spreadqtfus (qtfus,vdum,r,hcap,rmajor,nj,xdebug(1),wdum)
      if(fusionvb .gt. 0 .and. xdebug(1) .gt. 0.5)
     .       print *,'in fusion12,done spreadqtfus'



c
c ----------------------------------------------------------------------
c  calculate particle and energy sources due to beam-thermal fusion
c  the following is an analytic fit to the results of Jassby
c  PPPL-1280 (Princeton Plasma Physics Lab. 1976).
c  We neglect thermal t beam t fusion for now.
c  note that the beam is ASSUMED to be deuterium here
c  Will update to Bosch & Hale in future.                            HSJ
c ----------------------------------------------------------------------
c
c     beam fusion is off during beam initial iterations:
c
      if(fusionvb .eq. 1)print *,'no_beam_fusion =',no_beam_fusion
      if (no_beam_fusion .eq. 1)  go to 2300
      if (         ibeam .le. 1)  go to 2300   ! beam not currently on..
c                                   ..(never was or has been turned off)
c
c     the other possibility is that the beam is turned on,
c     but FREYA has not been called yet
c
      if (initbm .eq. 1)  go to 2300
c
c     in what follows, ent is the thermal tritium density
c
c      do j=1,nj
c        if (ifus .eq. 1)  ent = (1.0 - fd) * en(j,idt)
c        if (ifus .eq. 2 .and. it .ne. 0)  ent =    en(j,it)
c        bx = ent * 0.00931 * (16.0 * te(j)        )**1.5
c     .    / (ene(j) * (1.0 + (16.0 * te(j) / 106.0)**1.5))
c        do jb=1,nbeams
c          do ic=1,3
c            by = SQRT (ebeam(ic,jb)) * (1.0+(106.0/ebeam(ic,jb))**2.78)
c            epsb         = bx / by
c            qbf(j,ic,jb) = epsb * qb(j,ic,jb)
c            qbfus(j)     = qbfus(j) + qbf(j,ic,jb) * 0.62415064e16
c          end do
c        end do
c        sbfus(j) = qbfus(j) / 3.5e3
c      end do
c
c     replace the above with what follows note that qbf(j,ic,jb) is
c     no longer used. still exists in memory however.
c
      if (ifreya_old .ne. ifreya .and.  .not. use_nubeam ) then 
                    ! do these rather time-consuming
c                   !calcs only if the beam was recalculated and
                    !we are not using nubeam
         if (iddfusb_bulk .eq. 0) then
c
c          calculate the beam-thermal rates assuming no thermal or bulk
c          motion of ions:
c
c          call beam_thermal_fus0  ! routine does not exist, no real..
c                                  ..interest in this case.
c
         else if (iddfus .gt. 0) then
c
c          calculate the beam-thermal rates assuming thermal and bulk
c          motion of ions . for qbfus and sbfus we count only the
c          fast d, thermal t and fast t, thermal d rates:
c
           if (fusionvb .eq. 1  .and. .not. use_nubeam)
     .       write (*,'(" calling BEAM_THERMAL_FUS")')
           if(.not. use_nubeam)then
               call beam_thermal_fus (time, timmax, bpol, totcur(1),
     .                            qbfus, sbfus)
               if (beam_thermal_long_calc .lt. 0)
     .         call beam_thermal_approx_fus (time, timmax, qbfus, sbfus)
           endif
         end if
c
c        the beam_beam fusion rates are calculated but not used otherwise.
c        IN PARTICULAR, FAST ION DEPLETION DUE TO BEAM FUSION IS
c        NOT ACCOUNTED FOR IN THE FAST ION DISTRIBUTION FUNCTION !!!!!!!!HSJ
c
         if (fusionvb .eq. 1 .and.  .not. use_nubeam)
     .     write (*, '(" calling BEAM_BEAM_RATE from SOURCE")')
         if( .not. use_nubeam )
     .    call beam_beam_rate (time, time0, timmax, totcur(1))
         if (beam_beam_long_calc .lt. 0 .and.  .not. use_nubeam)
     .     call beam_beam_approx_rate(time,time0,timmax)
         ifreya_old = ifreya
         call copya (qbfus, qbfsav, nj)  ! save results for case when
         call copya (sbfus, sbfsav, nj)  ! when calcs are not done
      else                               ! restore sbfus, qbfus
         call copya (qbfsav, qbfus, nj)  ! restore previous result since
         call copya (sbfsav, sbfus, nj)  ! calcs were skipped this time
      end if
c
c ----------------------------------------------------------------------
c sum thermal and beam fusion terms:
c   account for alpha particle slowing down
c   save the critical alpha energy, ecritalpha, for printout
c ----------------------------------------------------------------------
c
 2300 call zeroa (palpha0, nj)           ! ASSUME NO directed MOMENTUM..
c                                        ..FOR THE TIME BEING - HSJ
      call zeroa (vzalpha, nj)
      rtstcxf=1.0
      atwf0  =4.0
      ezero0 =3.5e+3
      zsqf0  =4.0
      if(fusionvb .gt. 0)print *,'in fusion12,calling slow1'

      call slow1 (atw, atwf0, ene, en, enn, ezero0,
     .            0, 0, 1, kj, nj,
     .            nion, palpha0, te, vzalpha, zsq, zsqf0, dum, dum,
     .            ecritalpha, emzrat, fencap, ffe, ffi, dum, tausa,
     .            fionx,rtstcxf)
      if(fusionvb .gt. 0)print *,'in fusion12,done slow1'
c

      do j=1,nj
        qfus(j) = qtfus(j) + qbfus(j) 
        sfus(j) = stfus(j) + sbfus(j)
        if (ihe .ne. 0) then
          if (ineut(ihe) .ne. 0)  siadd(j,ihe) = siadd(j,ihe) + sfus(j)
        end if
      end do
c

c
c --- taupa (in sec) is taus*N for alphas (particle slowing-down time)
c --- taupa is defined as <nfast>/ndot,<nfast> is average fast ion density
c --- tauea is 0.5 * taus * ge  (energy slowing down time)
c --- tauea is defined as <nfast*efast>/(ndot*e0) where ndot is the source
c --- rate and e0 is the source energy.
c
         dtalpha = dt 
         IF(iaslow == 0)dtalpha = 10.0         ! assume asymptotic values at all times
                                               

         if(fusionvb .gt. 0)print *,' in fusion12 calling slow2'

         call slow2 (dum, dtalpha, fencap, enasav, ffe, 0, nj,
     .            dum, qfus, sfus, dum, tausa, wasav,
     .            enalp, dum, qfus, sfus, dum, taupa, dum,
     .            tauea, walp, 0, dum, dum, dum, dum)
c

         if(fusionvb .gt. 0)print *,'fusion12, done slow2'


 2302 do j=1,nj
        if (wtifus .gt.  0.0) ffi(j) = wtifus
        ffe(j)    = 1.0 - ffi(j)
        ffqt=0.0
        ffqb=0.0
        if (qtfus(j)+qbfus(j) .gt. 0.0)
     .            ffqt      = qtfus(j)/(qtfus(j)+qbfus(j))
                  ffqb      = 1.0 - ffqt
        qtfuse(j) = ffe(j)*ffqt*qfus(j)
        qtfusi(j) = ffi(j)*ffqt*qfus(j)
        qbfuse(j) = ffe(j)*ffqb*qfus(j)
        qbfusi(j) = ffi(j)*ffqb*qfus(j)
        qfuse(j)  = qtfuse(j) + qbfuse(j)
        qfusi(j)  = qtfusi(j) + qbfusi(j)
      end do

      if(fusionvb .gt. 0)print *,'leaving fusion12'

      return
      end

      subroutine fustable
c
      USE param
      USE fusion
      USE io
      USE ions
      USE neut
      use solcon
      USE soln
      USE mhdpar
      USE nub
      USE rf
      USE extra
      USE numbrs
      USE mesh
      USE sourc
      USE machin
      USE geom
      USE flags
      USE tcoef
      USE bd_condtn,only : bctime,ub,fluxb,
     .    ub_save,ub_rho_edge,bctime_zone
      USE flx
      USE flxav
      USE pelcom
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c      include 'param.i'
c      include 'mhdpar.i'
c      include 'bcon.i'
c      include 'extra.i'
c      include 'flags.i'
c      include 'flx.i'
c      include 'flxav.i'
c      include 'fusion.i'
c      include 'geom.i'
c      include 'io.i'
c      include 'ions.i'
c      include 'machin.i'
c      include 'mesh.i'
c      include 'neut.i'
c      include 'nub.i'
c      include 'numbrs.i'
c      include 'pelcom.i'
c      include 'rf.i'
c      include 'solcon.i'
c      include 'soln.i'
c      include 'sourc.i'
c      include 'tcoef.i'
c
      write  (nout, 8100)
 8100 format ('1', 'fusion table' //
     .         4x, 'j', 7x, 'r', 7x, 'taupa', 9x, 'tauea', 7x, 'enalp',
     .         7x, 'walp', 7x, 'balpha' /
     .        10x,'(cm)')
c
      do j=1,nj,jprt
        write  (nout, 8120)  j, r(j), taupa(j), tauea(j), enalp(j),
     .                       walp(j), balpha(j)
 8120   format (i5, f9.2, 9e13.3)
      end do
c
      write  (nout, 8130) betaa
 8130 format (/ '    average', 55x, e13.3)
      return
c
      end










      subroutine get_density (enbeam, enalp, ene, en, zsq, zeff, z,
     .                        itran, kj, j, nprim, nimp, ihe, id, it,
     .                        ncrt, nout, zfrac, time)
c
      implicit none
c
c ----------------------------------------------------------------------
c   subroutine determines the densities required if ifus=-1
c   at the spatial mesh point j
c
c   input
c   enbeam(j)        j=1,2...nj fast ion density
c   enalp(j)                    fast alpha density
c   ene(j)
c   en(j,k)                     k=1,2..nion
c   zeff
c   zsq
c   kj,j
c   nprim
c   nimp
c   ihe,id,it
c   ncrt,nout
c   zfrac
c   time
c
c  OUTPUT
c
c   zeff           may be changed
c   en(j,i)        i=impurity and analysis mode hydrogenic species
c
c ------------------------------------------------------------------ HSJ
c
      logical   exists, opened
      integer   kj, nprim, nimp, j, k, ihe, id, it, itran(*),
     .          ncrt, nout, ierd, iert, itry, itrymax, i,
     .          iozeff, iostat
      real*8    enbeam(*), enalp(*), ene(*), zsq(kj,*), z(kj,*),
     .          zeff(*), en(kj,*), en_he, en_i, en_e, en_a, en_b,
     .          z_h, z_a, z_b, z_i, z_he, zsq_h, zsq_a, zsq_b, zsq_i,
     .          zsq_he, z_sv, z_eff, en_h, d0, d1, d2, dz_eff, zfrac,
     .          time
      character timestamp*20
      data      itrymax, iozeff /20, 69/
c
      ierd = 0
      iert = 0
      itry = 0
      z_sv = zeff(j)
c
c --- if he is a primary ion species, its density is known (either from
c --- the input initial profile or from the current solution for
c --- t > time0). solve for the total hydrogenic (i.e., d+t) density
c --- and the impurity ion density at each spatial grid point j:
c
c     thermal he:
c


      en_he  = 0.0
      zsq_he = 4.0
      z_he   = 2.0
      if (ihe .ne. 0) then         ! he is a primary ion species
        en_he  = en(j,ihe)
        zsq_he = zsq(j,ihe)
        z_he   = z(j,ihe)
      end if
c
c     hydrogenic species:
c
      z_h   = 1.0
      zsq_h = 1.0
c
c     fast alpha:
c
      z_a   = 2.0
      zsq_a = 4.0
      en_a  = enalp(j)
c
c     beam:
c
      z_b   = 1.0
      zsq_b = 1.0
      en_b  = enbeam(j)
c
c     impurity:
c
      k     = nprim + 1                    ! pointer to impurity
      z_i   = z(j,k)
      zsq_i = zsq(j,k)
c
c     electrons:
c
      en_e  = ene(j)
      z_eff = zeff(j)
c
c     for hydrogenic,he,beam,alpha,impurity respectively we have
c
c           n_e = z_h*n_h + z_he*n_he + z_b*n_b + z_a*n_a +z_i*n_i
c     z_eff*n_e = zsq_h*n_h + zsq_he*n_he + zsq_b*n_b + zsq_a*n_a
c                 + zsq_i*n_i
c
c     solve these equations for the hydrogenic(n_h) and
c     impurity (n_i) densities using Cramers rule:
c
  100 iert  = 0
      ierd  = 0
      d0    = zsq_h*z_i - zsq_i*z_h
      d1    =    z_eff*en_e - zsq_he*en_he - zsq_b*en_b
     .          -zsq_a*en_a
      d2    =   en_e - z_he*en_he - z_b*en_b -z_a*en_a
      en_h  = z_i*d1 - zsq_i*d2
      en_i  = zsq_h*d2 - z_h*d1
      if (d0 .ne. 0.0) then
        en_h = en_h/d0
        en_i = en_i/d0
      else
        write  (ncrt, 1)
        write  (nout, 1)
    1   format (' subroutine GET_DENSITY reports:' /
     .          '   unable to determine hydrogenic',
     .          ' and impurity densities' /
     .          '   equations are degenerate')
        call STOP ('subroutine GET_DENSITY: degenerate eqns', 213)
      end if
c
      en(j,k) = en_i                    ! set impurity density
c
      if (zfrac .gt. 0.0 .and. zfrac .lt. 1.0) then
        if (it .lt. id) then
          en(j,it)=zfrac*en_h
          en(j,id)=(1.0-zfrac)*en_h
        else
          en(j,id)=zfrac*en_h
          en(j,it)=(1.0-zfrac)*en_h
        end if
        if (en(j,id) .le. 0.0) ierd = 1
        if (en(j,it) .le. 0.0) iert = 1
      else
c
c       get the density of the hydrogenic species run in analysis mode
c
         do i=1,nprim
           if (itran(i) .eq. 1 .and. i .eq. id) then   ! t density
             en(j,it)=en_h-en(j,id)           ! is in analyis mode
             if (en(j,it) .le. 0.0)  iert = 1
            end if
           if (itran(i) .eq. 1 .and. i .eq. it) then   ! d density
             en(j,id)=en_h-en(j,it)           ! is in analyis mode
             if (en(j,id) .le. 0.0)  ierd = 1
           end if
         end do
       end if
c

       if (iert .eq. 0 .and. ierd .eq. 0) then
         if (zeff(j) .ne. z_sv) then
c
c          save the zeff corrections in a file named "zeff.mod"
c
           inquire (file = 'zeff.mod', iostat = iostat,
     .              exist = exists, opened = opened)
           if (iostat .ne. 0) then         ! problem with INQUIRE
             write (ncrt, '(/ 2a, i6)')
     .                    ' ERROR: Fatal INQUIRE failure,',
     .                           ' IOSTAT =', iostat
            call STOP ('subroutine GET_DENSITY: bad INQUIRE', 209)
           end if
   20      if (exists .and. opened) then   ! file is open
             write  (iozeff, 10) time, j, zeff(j), z_sv,
     .                           en(j,id), en(j,it), en(j,k)
   10        format ('time = ', f12.6, 4x,
     .               'j, zeff_new, zeff_old = ', i3, 2f10.3 /
     .               'en(j,id), en(j,it), en(j,k) =',
     .                3(1x, 1pe12.4) /)
           else                            ! file is not open
             if (exists) then              ! file exists; trash it
               call DESTROY ('zeff.mod')
             end if
             call getioun(iozeff,iozeff)
             open (unit = iozeff, file = 'zeff.mod',
     .             status = 'NEW', iostat = iostat)
             if (iostat .ne. 0) then
               write (ncrt, '(/ 2a, i6)')
     .                      ' ERROR: Fatal OPEN failure,',
     .                             ' IOSTAT =', iostat
               call giveupus(iozeff)
               call STOP ('subroutine GET_DENSITY: bad OPEN', 216)
             end if
             call GET_DATE_TIME (timestamp)
             write (iozeff, '(2a /)')
     .                      'FILE zeff.mod CREATED  ', timestamp
             exists = .true.
             opened = .true.
             go to 20
           end if
         end if
         return                       ! required density was found
       end if
c
c     problem with charge conservation try to correct it
c     by adjusting zeff. Assume zero impurity:
c
    
      if (itry .eq. 0) then
        iert = 0
        ierd = 0
        en_h = en_e -  z_he*en_he - z_b*en_b - z_a*en_a
        do i=1,nprim
          if (itran(i) .eq. 1 .and. i .eq. id) then   ! t density
            en(j,it) = en_h-en(j,id)        ! is in analyis mode
            if (en(j,it) .le. 0.0)  iert = 1
          end if
          if (itran(i) .eq. 1 .and. i .eq. it) then  ! d density
            en(j,id) = en_h-en(j,it)        ! is in analyis mode
            if (en(j,id) .le. 0.0)  ierd = 1
          end if
        end do
        if (ierd .ne. 0 .or. iert .ne. 0) then
          write  (ncrt, 2)
          write  (nout, 2)
    2     format (/' subroutine GET_DENSITY detected an ERROR:' /
     .             '   hydrogenic density is negative even',
     .             ' with impurity density set to zero')
          call STOP ('subroutine GET_DENSITY: zero density', 214)
        else  ! recovery is possible with some impurity density
          z_eff = (en_h+zsq_he*en_he+zsq_b*en_b+zsq_a*en_a)/en_e
          ierd  = 0
          iert  = 0
        end if
      end if
      itry  = itry + 1
      if (itry .eq. 1)  dz_eff = (z_sv-z_eff)/(itrymax+1)
      z_eff = z_sv-itry*dz_eff
c
      write (*, '(" adjusting zeff at j =",i3)')j
      write (*, '(" old zeff, new zeff = ", 2(2x, f12.6))')
     .                  zeff(j), z_eff
      write (*, '(" en(j,id) =", 1pe14.4)') en(j,id)
      write (*, '(" en(j,it) =", 1pe14.4)') en(j,it)
      write (*, '(" en_h     =", 1pe14.4)') en_h
      write (*, '(" en_i     =", 1pe14.4)') en_i
      write (*, '(" zfrac    =", 1pe14.4)') zfrac
      zeff(j) = z_eff
      if (itry .lt. itrymax)  go to 100
      call STOP ('subroutine GET_DENSITY: exceeded ITRYMAX', 215)
      return
c
      end

      subroutine get_external_rfcur
c
      USE param
      USE mhdpar
      USE io
      USE rf
      USE numbrs
      USE mesh
      USE sourc
      USE geom
      USE constnts
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c --- subroutine adds (constant in time) RF current profile read in
c     from inone. subroutine uses temporary storage wdum,vdum,sdum

      include 'storage.i'
c
      real*8       one
      dimension    rova(kj), rf_cur(kj), cspline(kpsi,3), bparrf(4)
      equivalence (cspline(1,1), vdum(1)     )
      equivalence (        rova, sdum(1)     )
      equivalence (   rf_cur(1), sdum(kj+1)  )
      equivalence (   bparrf(1), sdum(2*kj+1))
c
      logical      check_extrf_input(krf)
      data        (check_extrf_input(j), j=1,krf)
     .             / krf * .true. /           ! check only on first call
c
      one = 1.0
      do l=1,krf
      if (ABS (extcurrf(l)) .gt. 1.0e-8) then
      if (extcurrf_nj(l) .le. 2 .or. extcurrf_nj(l) .gt. kj) then
         write  (ncrt, 1) extcurrf_nj(l)
         write  (nout, 1) extcurrf_nj(l) 
    1    format (/
     .   ' subroutine GET_EXTERNAL_RFCUR 2 reports:'               /
     .   '   error in specification of number of current values' /
     .   '   extcurrf_nj(l) = ', i5)
         call STOP ('subroutine GET_EXTERNAL_RFCUR: problem #1', 185)
      else if (extcurrf_nj(l) .ne. nj) then
c
c        this branch is executed only once regardless of
c        value of check_extrf_input
c        do spline fitting first time routine is called:
c
         if (extcurrf_rho(             1,l) .ne. 0.0 .or.
     .       extcurrf_rho(extcurrf_nj(l),l) .ne. 1.0) then
           write  (ncrt, 2)
           write  (nout, 2)
    2      format (/ ' subroutine GET_EXTERNAL_RFCUR reports:' /
     .               '   start or end of extcurrf_rho',
     .                 ' not specified correctly')
           write  (ncrt, 3) (extcurrf_rho(j,l), j=1,extcurrf_nj(l))
           write  (nout, 3) (extcurrf_rho(j,l), j=1,extcurrf_nj(l))
    3      format (' extcurrf_rho =', (6(2x,1pe12.4)))
           call STOP ('subroutine GET_EXTERNAL_RFCUR: problem #2', 186)
         end if
         if (extcurrf_nj(l) .gt. kpsi) then
           write  (ncrt, 4)  kpsi, extcurrf_nj(l)
           write  (nout, 4)  kpsi, extcurrf_nj(l)
    4      format (/ ' subroutine GET_EXTERNAL_RFCUR reports:'         /
     .               '   spline representation of external RF current'
     .                 ' is limited to ', i5, ' knots.'                /
     .               '   In inone file, ', i5, ' knots are set')
           call STOP ('subroutine GET_EXTERNAL_RFCUR: problem #3', 187)
         end if
         rf_ext_curtot(l) = 0.0
         do j=1,extcurrf_nj(l)
           rf_ext_curtot(l) = rf_ext_curtot(l) + extcurrf_curr(j,l)**2
         end do
         if (rf_ext_curtot(l) .le. 0.0) then
           write  (nout, 5)
           write  (ncrt, 5)
 5         format (/ ' subroutine GET_EXTERNAL_RFCUR reports:' /
     .               '   error in input of extcurrf_curr'      /
     .               '   input current cannot be identically zero')
           call STOP ('subroutine GET_EXTERNAL_RFCUR: problem #4', 188)
         end if
         tension   = 1.0
         call copya (r, rova, nj) ! load r into LOCAL variable rova
         call multpl1 (rova, nj, 1.0/r(nj))        ! normalize rova
         rova(nj)  = 1.0                           ! roundoff correction
         call multpl1 (extcurrf_curr(l,l), extcurrf_nj(l), extcurrf(l))
         bparrf(1) = -1.0e30      ! zero gradient  at rhoa = 0
         bparrf(2) =  0.0
         bparrf(3) =  0.0         ! natural spline at rhoa = 1
         bparrf(4) =  0.0
         call tspline (extcurrf_rho(1,l),extcurrf_curr(1,l),
     .              extcurrf_nj(l),bparrf,cspline,kpsi,ier,tension,
     .              wdum(1),wdum(kpsi+1),wdum(2*kpsi+1),wdum(3*kpsi+1),
     .              wdum(4*kpsi+1),wdum(5*kpsi+1),tenmax,rova,
     .                 rf_cur,nj)
c
c        normalize if called for (note that this normalization assumes
c        that the RF current is in fact <JRF*R0/R>, which is not true but
c        is consistent with the way it's done in other sections of the code.)
c
         if (extcurrf_amps(l) .ne. 0.0) then
           call trapv (r, rf_cur, hcap, nj, rf_ext_curtot(l))
           xrfn             = extcurrf_amps(l)/(twopi*rf_ext_curtot(l))
           call multpl1 (rf_cur, nj, xrfn)
           rf_ext_curtot(l) = extcurrf_amps(l)    ! save it for printout
           extcurrf_amps(l) = 0.0                 ! for subsequent calls
         end if
         call copya (rf_cur, extcurrf_curr(1,l), nj)
         call addac (extcurrf_curr(1,l), currf, nj, one)
         extcurrf_nj(l)       =  nj               ! for subsequent calls
         check_extrf_input(l) = .false.           ! for subsequent calls
         extcurrf(l)          = 1.0               ! for subsequent calls
      else if (extcurrf_nj(l) .eq. nj) then
c
c        on first call we go through here if extcurrf_nj=nj in inone.
c        on subsequent calls to this routine  we always go through
c        here because extcurrf_nj is forced to 51 after the first call.
c
         if (check_extrf_input(l)) then
           rf_ext_curtot(l) = 0.0
           do j=1,extcurrf_nj(l)
             rf_ext_curtot(l) = rf_ext_curtot(l) + extcurrf_curr(j,l)**2
           end do
           if (rf_ext_curtot(l) .le. 0.0) then
             write  (nout, 6)
             write  (ncrt, 6)
    6        format (/ ' subroutine GET_EXTERNAL_RFCUR reports:' /
     .                 '   error in input of extcurrf_curr'      /
     .                 '   input current cannot be identically zero')
            call STOP ('subroutine GET_EXTERNAL_RFCUR: problem #5', 189)
           end if
           if (ABS (extcurrf(l)-1.0) .gt. 1.0e-5) then
             call multpl1 (extcurrf_curr(1,l), nj, extcurrf(l))
             extcurrf(l) = 1.0
           end if
           if (extcurrf_amps(l) .ne. 0.0) then
              extcurrf(l) = 1.0
              call trapv (r, extcurrf_curr(1,l), hcap, nj,
     .                    rf_ext_curtot(l))
              xrfn = extcurrf_amps(l) / (twopi*rf_ext_curtot(l))
              call multpl1 (extcurrf_curr(1,l), nj, xrfn)
              rf_ext_curtot(l) = extcurrf_amps(l) ! save it for printout
              extcurrf_amps(l) = 0.0              ! for subsequent calls
              extcurrf(l) = 1.0  ! this multiplier not used in this case
            else
              call trapv (r, extcurrf_curr(1,l), hcap, nj,
     .                    rf_ext_curtot(l))
                rf_ext_curtot(l) =  twopi*rf_ext_curtot(l)  
               extcurrf_amps(l) = 0.0              ! for subsequent calls
               extcurrf(l) = 1.0  ! this multiplier not used in this case

           end if
           check_extrf_input(l) = .false. ! no check on subsequent calls
         end if
c
c        extcurrf_curr is now properly set up so just add it in
c        (i.e., set  currf=currf+1.0*extcurrf_curr )
c
         call addac (extcurrf_curr(1,l), currf, nj, one)
      end if
      end if
      end do
c
      return
c
      end

      subroutine get_external_rfqe
c
      USE param
      USE mhdpar
      USE io
      USE rf
      USE numbrs
      USE mesh
      USE sourc
      USE machin
      USE geom
      USE constnts
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c --- subroutine adds (constant in time) RF electron heating profile
c     read in from inone. subroutine uses temporary storage wdum,vdum,sdum
c
c      include 'param.i'
c      include 'constnts.i'
c      include 'geom.i'
c      include 'numbrs.i'
c      include 'io.i'
c      include 'machin.i'
c      include 'mesh.i'
c      include 'mhdpar.i'
c      include 'rf.i'
c      include 'sourc.i'
      include 'storage.i'
c
      real*8       kev_per_joule
      dimension    rova(kj),rf_qe(kj), cspline(kpsi,3),bparrf(4)
      equivalence (cspline(1,1), vdum(1))
      equivalence (rova,sdum(1))
      equivalence (rf_qe(1),sdum(kj+1))
      equivalence (bparrf(1),sdum(2*kj+1))
      logical      check_extrf_input(krf)
      data        (check_extrf_input(j), j=1,krf)
     .             / krf * .true. /           ! check only on first call
c
      kev_per_joule = 6.242e15
c
      do l=1,krf
      if (                       ABS (extqerf(l)) .gt. 1.0e-8) then
      if (extqerf_nj(l) .le. 2 .or. extqerf_nj(l) .gt. kj    ) then
         write  (ncrt, 1) extqerf_nj(l)
         write  (nout, 1) extqerf_nj(l)
    1    format (/
     .   ' subroutine GET_EXTERNAL_RFQE reports:'                      /
     .   '   error in specification of number of power density values' /
     .   '   extqerf_nj(l) = ', i5)
         call STOP ('subroutine GET_EXTERNAL_RFQE: problem #1', 190)
      else if (extqerf_nj(l) .ne. nj) then
c
c        this branch is executed only onece regardless of
c        value of check_extrf_input because extqerf_nj will be
c        set to nj below.
c        do spline fitting first time routine is called:
c
         if (extqerf_rho(            1,l) .ne. 0.0 .or.
     .       extqerf_rho(extqerf_nj(l),l) .ne. 1.0) then
             write  (ncrt, 2)
             write  (nout, 2)
    2        format (/ ' subroutine GET_EXTERNAL_RFQE reports:' /
     .                 '   start or end of extqerf_rho',
     .                   ' not specified correctly')
             write  (ncrt, 3) (extqerf_rho(j,l),j=1,extqerf_nj(l))
             write  (nout, 3) (extqerf_rho(j,l),j=1,extqerf_nj(l))
    3        format (' extqerf_rho =',(6(2x,1pe12.4)))
             call STOP ('subroutine GET_EXTERNAL_RFQE: problem #2', 191)
         end if
         if (extqerf_nj(l) .gt. kpsi) then
             write  (ncrt, 4)  kpsi, extqerf_nj(l)
             write  (nout, 4)  kpsi, extqerf_nj(l)
    4        format (/ ' subroutine GET_EXTERNAL_RFQE reports:'       /
     .                 '   spline representation of external RF'      /
     .                 '   electron power limited to ', i5, ' knots.' /
     .                 '   In inone file ', i5, ' knots are set')
             call STOP ('subroutine GET_EXTERNAL_RFQE: problem #3', 192)
         end if
         rf_ext_qetot(l) = 0.0
         do j=1,extqerf_nj(l)
           rf_ext_qetot(l) = rf_ext_qetot(l)+extqerf_qe(j,l)**2
         end do
         if (rf_ext_qetot(l) .le. 0.0) then
            write  (nout, 5)
            write  (ncrt, 5)
 5          format (/ ' subroutine GET_EXTERNAL_RFQE reports:'       /
     .                '   error in input of extqerf_qe, input power' /
     .                '   cannot be identically equal to zero')
            call STOP ('subroutine GET_EXTERNAL_RFQE: problem #4', 193)
         end if
          tension    = 1.0
          call copya (r, rova, nj)     ! load r into LOCAL variable rova
          call multpl1 (rova, nj, 1.0/r(nj))       ! normalize rova
          rova(nj)=1.0                             ! roundoff correction
          call multpl1 (extqerf_qe(1,l), extqerf_nj(l), extqerf(l))
c
c         extqerf_qe is now in watts/cm**3 even if it was keV/(cm**3*sec)
c         on input, because it is assumed that the user has set
c         extqerf = 1.602e-16 if extqerf_qe is input in keV/(cm**3*sec)
c         in the inone file.
c
          bparrf(1) = -1.0e30     ! zero gradient  at rhoa = 0
          bparrf(2) =  0.0
          bparrf(3) =  0.0        ! natural spline at rhoa = 1
          bparrf(4) =  0.0
          call tspline (extqerf_rho(1,l),extqerf_qe(1,l),
     .             extqerf_nj(l),bparrf,cspline,kpsi,ier,tension,
     .             wdum(1),wdum(kpsi+1),wdum(2*kpsi+1),wdum(3*kpsi+1),
     .             wdum(4*kpsi+1),wdum(5*kpsi+1),tenmax,rova,
     .             rf_qe,nj)
c
c        normalize if called for(dv=4*pi**2*R0*H*r*dr)
c
         if (extqerf_watts(l) .ne. 0.0) then
            call trapv (r, rf_qe, hcap, nj,rf_ext_qetot(l))
            xrfn = extqerf_watts(l)/(twopi*twopi*rmajor*rf_ext_qetot(l))
            call multpl1 (rf_qe, nj, xrfn)
            rf_ext_qetot (l) = extqerf_watts(l)   ! save it for printout
            extqerf_watts(1) = 0.0                ! for subsequent calls
         else
            call trapv (r, rf_qe, hcap, nj,rf_ext_qetot(l))
            rf_ext_qetot(l)=twopi*twopi*rmajor*rf_ext_qetot(l)
         end if
         call copya (rf_qe, extqerf_qe(1,l), nj)
c
c        qrfe is in keV/(cm**3*sec) so convert extqerf_qe here
c        (6.242e15 converts from watts to keV/sec)
c        note that extqerf_qe remains in watts/cm**3
c
         call addac (extqerf_qe(1,l), qrfe, nj, kev_per_joule)
c
         extqerf_nj(l)        =  nj               ! for subsequent calls
         check_extrf_input(l) = .false.           ! for subsequent calls
         extqerf(l)           = 1.0               ! for subsequent calls
      else if (extqerf_nj(l) .eq. nj) then
c
c        on first call we go through here if extqerf_nj=nj in inone.
c        on subsequent calls to this routine  we always go through
c        here because extqerf_nj is forced to nj  after the first call.
c
         if (check_extrf_input(l)) then
           rf_ext_qetot(l) = 0.0
           do j=1,extqerf_nj(l)
             rf_ext_qetot(l) = rf_ext_qetot(l)+extqerf_qe(j,l)**2
           end do
           if (rf_ext_qetot(l) .le. 0.0) then
             write  (nout, 6)
             write  (ncrt, 6)
 6           format (/' subroutine GET_EXTERNAL_RFQE reports:'         /
     .                '   error in input of extqerf_qe, input current' /
     .                '   cannot be identically equal to zero')
             call STOP ('subroutine GET_EXTERNAL_RFQE: problem #5', 194)
           end if
           if (ABS (extqerf(l)-1.0) .gt. 1.0e-5) then
c
c            only on first call execute these lines
c
             call multpl1 (extqerf_qe(1,l), nj, extqerf(l))
             extqerf(l) = 1.0
           end if
           if (extqerf_watts(l) .ne. 0.0) then
             call trapv (r, extqerf_qe(1,l), hcap, nj,rf_ext_qetot(l))
             xrfn=extqerf_watts(l)/(twopi*twopi*rmajor*rf_ext_qetot(l))
             call multpl1 (extqerf_qe(1,l), nj, xrfn)
             rf_ext_qetot (l) = extqerf_watts(l)  ! save it for printout
             extqerf_watts(l) = 0.0               ! for subsequent calls
             extqerf(l) = 1.0    ! this multiplier not used in this case

            ELSE
              call trapv (r, extqerf_qe(1,l), hcap, nj,rf_ext_qetot(l))
              rf_ext_qetot (l) = twopi*twopi*rmajor*rf_ext_qetot(l) 
              extqerf_watts(l) = 0.0               ! for subsequent calls
              extqerf(l) = 1.0    ! this multiplier not used in this case
           end if
           check_extrf_input(l) = .false. ! no check on subsequent calls
         end if
c
c        extqerf_qe is now properly set up so just add it in
c        (i.e., set  qrfe = qrfe + 6.242e15 * extqerf_qe )
c
         call addac (extqerf_qe(1,l), qrfe, nj, kev_per_joule)
      end if
      end if
      end do
      return
c
      end

      subroutine get_external_rfqi
c
      USE param
      USE io
      USE mhdpar
      USE rf
      USE numbrs
      USE mesh
      USE sourc
      USE machin
      USE geom
      USE constnts
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c --- subroutine adds (constant in time) RF heating profile read in
c     from inone. subroutine uses temporary storage wdum,vdum,sdum
c
c      include 'param.i'
c      include 'constnts.i'
c      include 'geom.i'
c      include 'numbrs.i'
c      include 'io.i'
c      include 'machin.i'
c      include 'mesh.i'
c      include 'mhdpar.i'
c      include 'rf.i'
c      include 'sourc.i'
      include 'storage.i'
c
      real*8       kev_per_joule
      dimension    rova(kj),rf_qi(kj), cspline(kpsi,3),bparrf(4)
      equivalence (cspline(1,1), vdum(1))
      equivalence (rova,sdum(1))
      equivalence (rf_qi(1),sdum(kj+1))
      equivalence (bparrf(1),sdum(2*kj+1))
      logical      check_extrf_input(krf)
      data        (check_extrf_input(j), j=1,krf)
     .             / krf * .true. /           ! check only on first call
c
      kev_per_joule = 6.242e15
c
      do l=1,krf
      if (ABS (extqirf(l)) .gt. 1.0e-8) then
      if (extqirf_nj(l) .le. 2 .or. extqirf_nj(l) .gt. kj) then
         write  (ncrt, 1)
         write  (nout, 1)
    1    format (/
     .   ' subroutine GET_EXTERNAL_RFQI reports:' /
     .   '   error in specification of number of power density values' /
     .   '   extqirf_nj = ', i5)
         call STOP ('subroutine GET_EXTERNAL_RFQI: problem #6', 195)
      else if (extqirf_nj(l) .ne. nj) then
c
c        this branch is executed only onece regardless of
c        value of check_extrf_input
c        do spline fitting first time routine is called:
c
         if (extqirf_rho(         1,l) .ne. 0.0 .or.
     .       extqirf_rho(extqirf_nj(l),l) .ne. 1.0) then
             write  (ncrt, 2)
             write  (nout, 2)
    2        format (/ ' subroutine GET_EXTERNAL_RFQI reports:' /
     .                 '   start or end of extqirf_rho',
     .                   ' not specified correctly')
             write  (ncrt, 3)  (extqirf_rho(j,l), j=1,extqirf_nj(l))
             write  (nout, 3)  (extqirf_rho(j,l), j=1,extqirf_nj(l))
    3        format (' extqirf_rho =', (6(2x,1pe12.4)))
             call STOP ('subroutine GET_EXTERNAL_RFQI: problem #7', 196)
         end if
         if (extqirf_nj(l) .gt. kpsi) then
           write  (ncrt, 4)  kpsi, extqirf_nj(l)
           write  (nout, 4)  kpsi, extqirf_nj(l)
    4      format (/ ' subroutine GET_EXTERNAL_RFQI reports:'     /
     .               '   spline representation of external RF '   /
     .               '   ion power is limited to ', i5, ' knots.' /
     .               '   In inone ', i5, ' knots are set')
           call STOP ('subroutine GET_EXTERNAL_RFQI: problem #8', 197)
         end if
         rf_ext_qitot(l) = 0.0
         do j=1,extqirf_nj(l)
           rf_ext_qitot(l)=rf_ext_qitot(l)+extqirf_qi(j,l)**2
         end do
         if (rf_ext_qitot(l) .le. 0.0) then
            write  (nout, 5)
            write  (ncrt, 5)
    5       format (/ ' subroutine GET_EXTERNAL_RFQI reports:'       /
     .                '   error in input of extqirf_qi, input power' /
     .                '   cannot be identically equal to zero')
            call STOP ('subroutine GET_EXTERNAL_RFQI: problem #9', 198)
         end if
         tension  = 1.0
         call copya (r, rova, nj)     ! load r into LOCAL variable rova
         call multpl1 (rova, nj, 1.0/r(nj)) ! normalize rova
         rova(nj) = 1.0                     ! roundoff correction
         call multpl1 (extqirf_qi(1,l), extqirf_nj(l), extqirf(l))
c
c        extqirf_qi is now in watts/cm**3 even if it was keV/cm**3*sec on input
c        because it is assumed that the user has set extqirf = 1.602e-16
c        if extqirf_qi is input in keV/(cm**3*sec) in the inone file.
c
         bparrf(1) = -1.0e30      ! zero gradient  at rhoa = 0
         bparrf(2) =  0.0
         bparrf(3) =  0.0         ! natural spline at rhoa = 1
         bparrf(4) =  0.0
         call tspline (extqirf_rho(1,l),extqirf_qi(1,l),extqirf_nj(l),
     .            bparrf,cspline,kpsi,ier,tension,
     .            wdum(1),wdum(kpsi+1),wdum(2*kpsi+1),wdum(3*kpsi+1),
     .            wdum(4*kpsi+1),wdum(5*kpsi+1),tenmax,rova,
     .                 rf_qi,nj)
c
c        normalize if called for 
         if (extqirf_watts(l) .ne. 0.0) then
            call trapv (r, rf_qi, hcap, nj,rf_ext_qitot(l))
            xrfn = extqirf_watts(l)/(twopi*twopi*rmajor*rf_ext_qitot(l))
            call multpl1 (rf_qi, nj, xrfn)
            rf_ext_qitot (l) = extqirf_watts(l)   ! save it for printout
            extqirf_watts(l) = 0.0                ! for subsequent calls
         else
            call trapv (r, rf_qi, hcap, nj,rf_ext_qitot(l))
            rf_ext_qitot(l)=twopi*twopi*rmajor*rf_ext_qitot(l)
         end if
         call copya (rf_qi, extqirf_qi(1,l), nj)
c
c        qrfi is in keV/(cm**3*sec) so convert extqerf_qi here
c        (6.242e15 converts from watts to keV/sec)
c        note that extqirf_qi remains in watts/cm**3
c
         call addac (extqirf_qi(1,l), qrfi, nj, kev_per_joule)
         extqirf_nj(l)        = nj                ! for subsequent calls
         check_extrf_input(l) = .false.           ! for subsequent calls
         extqirf(l)           = 1.0               ! for subsequent calls
      else if (extqirf_nj(l) .eq. nj) then
c
c        on first call we go through here if extqirf_nj=nj in inone.
c        on subsequent calls to this routine  we always go through
c        here because extqirf_nj is forced to 51 after the first call.
c
         if (check_extrf_input(l)) then
           rf_ext_qitot(l) = 0.0
           do j=1,extqirf_nj(l)
             rf_ext_qitot(l) = rf_ext_qitot(l)+extqirf_qi(j,l)**2
           end do
           if (rf_ext_qitot(l) .le. 0.0) then
            write  (nout,6)
            write  (ncrt,6)
    6       format (/ ' subroutine GET_EXTERNAL_RFQI reports:'         /
     .                '   error in input of extqirf_qi, input current' /
     .                '   cannot be identically equal to zero')
            call STOP ('subroutine GET_EXTERNAL_RFQI: problem #10', 199)
           end if
           if (ABS (extqirf(l)-1.0) .gt. 1.0e-5) then
             call multpl1 (extqirf_qi(1,l), nj, extqirf(l))
             extqirf(l) = 1.0
           end if
           if (extqirf_watts(l) .ne. 0.0) then
              extqirf(l)       = 1.0
              call trapv (r, extqirf_qi(1,l),hcap, nj,rf_ext_qitot(l))
              xrfn             = extqirf_watts(l) / (twopi*twopi*rmajor*
     .                                               rf_ext_qitot(l))
              call multpl1 (extqirf_qi(1,l), nj, xrfn)
              rf_ext_qitot (l) = extqirf_watts(l) ! save it for printout
              extqirf_watts(l) = 0.0              ! for subsequent calls
              extqirf(l) = 1.0   ! this multiplier not used in this case
           else
              call trapv (r, extqirf_qi(1,l),hcap, nj,rf_ext_qitot(l))
              rf_ext_qitot(l) = twopi*twopi*rmajor*rf_ext_qitot(l)
              extqirf_watts(l) = 0.0              ! for subsequent calls
              extqirf(l) = 1.0   ! this multiplier not used in this case
           end if
           check_extrf_input(l) = .false. ! no check on subsequent calls
         end if
c
c        extqirf_qi is now properly set up so just add it in
c        (i.e., set  qrfi = qrfi+6.242e15*extqirf_qi )
c
         call addac (extqirf_qi(1,l), qrfi, nj, kev_per_joule)
      end if
      end if
      end do
      return
c
      end

      subroutine getrad
c
      USE param
      USE numbrs
      USE sourc
      USE bd_condtn
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c --- subroutine GETRAD interpolates for the radiation profile
c --- from the given radiation profile input data
c
c      include 'param.i'
      include 'imsl.i'
c      include 'invers.i'
c      include 'numbrs.i'
c      include 'sourc.i'
c
      imslmd = 'getrad  '
      call tsplin (qradin, qradr, nqrad, qrad)
      do j=1,nj
        qrad(j) = qrad(j) * 0.625e16
      end do
      return
c
      end

      subroutine getrho (xpos, ypos, zpos, rhopos)
c
      USE param
      USE mhdpar
      USE mhdgrid
      USE numbrs
      USE mesh
      USE machin
      USE geom
      USE rhog
      USE gpsi
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c --- subroutine GETRHO calculates the minor radius parameter
c --- at the spatial coordinate (xpos,ypos,zpos)
c
c      include 'param.i'
c      include 'mhdpar.i'
c      include 'geom.i'
c      include 'machin.i'
c      include 'mesh.i'
c      include 'numbrs.i'
c      include 'mhdgrid.i'
c      include 'rhog.i'
c      include 'small.i'
c
      rpos = SQRT (xpos**2 + ypos**2)
      if (codeid .eq. 'onedee')
     .  rhopos = SQRT (kappa*(rpos-rmajor)**2 + (zpos-zax)**2/kappa)
      if (codeid .ne. 'onedee')
     .  call bilin(nx,ny,rmhdgrid,zmhdgrid,p,rpos,zpos,ppos)
      if (codeid .ne. 'onedee')
     .  call interp(ppos,psir,nj,r,rhopos)
      return
c
      end



      integer *4 function get_time_dep_beam()
c-------------------------------------------------------------
c     fetch value of time_dep_beam from nub4.i:
      USE param
      USE nub
      USE nub4,only : time_dep_beam
c      include 'param.i'
c      include 'nub.i'
c      include 'nub4.i'


      get_time_dep_beam = time_dep_beam

      return
      end





      subroutine get_vloop_bc(t_want,v_want)
c  -------------------------------------------------------------
c     interpolate the loop voltage given in (vloop_bc_time,vloop_bc)
c     at time t_want onto v_want
c  -------------------------------------------------------------HSJ
      USE bd_condtn,only : vloop_bc,vloop_bc_time,nvloop
      IMPLICIT NONE
      real *8 t_want, v_want
          call interp (t_want,vloop_bc_time, nvloop, vloop_bc, v_want)
      return
      end







      real*8 function he3dprate (ti)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ------------------------------------------------- 11/20/95 --- HSJ ---
c             returns rate(cm**3/sec) of he3(d,p)he4 reaction
c             New Bosch & Hale rate coefficient:
c             Bosch & Hale, Nuc. Fus., vol32, no.4 (1992) 611
c
      real*8 mrcsq
c
c     data for d(d,p)t:
c
      data  C1,C2,C3,C4,C5,B_gsq, mrcsq
     .    /  5.51036e-10, 6.41918e-3, -2.02896e-3,
     .      -1.91080e-05, 1.35776e-4, 4726.672, 1124572 /
c
      if (ti .lt. 0.5) then
        he3dprate = 0.0
      else if (ti .le. 190.0) then
        theta  = ti*(C2+ti*C4)/(1.0+ti*(C3+ti*C5))
        theta  = ti/(1.0-theta)
        xsi    = (B_gsq/(4.0*theta))**(0.33333333334)
        he3dprate = C1 * theta * SQRT (xsi/(mrcsq*ti**3))
     .                       *  EXP (-3.0*xsi)
      else
        call STOP ('function He3DPRATE: TI out of range', 38)
      end if
      return
c
      end

      subroutine header (iunit, timet, t)
c
      USE param
      USE io
      USE solcon
      USE mhdpar
      USE yoka 
      USE numbrs
      USE soln2d
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c this subroutine generates the top line of each page of transport output
c ----------------------------------------------------------------------

c
      character*8 time_avg
c
      time_avg = ' '
      if (timav .ne. 0.0 .and. ilastp .eq. 1)  time_avg = 'time-avg'
      if (iyoka .ne. 0                      )  go to 20
      write  (iunit, 10)  timet, t, ieq, itre
   10 format ('1time = ', 3pf14.2, ' ms,  time point=', 0pf14.1,
     .     '   equilibrium point = ', i3,
     .     '   equilibrium iteration number = ', i2, '  transport.' /)
      return
c
   20 write  (iunit, 30)  timet, t, ishot, itime, time_avg,
     .                    versid, ieq, itre
   30 format ('1time = ', 3pf14.2, ' ms,  time point=', 0pf14.1,
     .     '   Analysis of shot', i8, ' at ', i5, ' ms  ', a8,
     .                              '  (onetwo version: ', a, ')' /
     .     '   equilibrium point = ', i3,
     .     '   equilibrium iteration number = ', i2                  /)
      return
c
      end

      subroutine icsevu1 (x, y, nx, c, ic, u, s, m, ier)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c --- modified version of IMSL subroutine named ICSEVU
c
c ----------------------------------------------------------------------
c
c   computer            - dec10/single
c
c   latest revision     - january 1, 1978
c
c   purpose             - evaluation of a cubic spline
c
c   usage               - call icsevu1 (x, y, nx, c, ic, u, s, m, ier)
c
c   arguments    x      - vector of length nx containing the abscissae
c                           of the nx data points (x(i),y(i)) i = 1,...,
c                           nx (input). x must be ordered so that
c                           x(i) .lt. x(i+1).
c                y      - vector of length nx containing the ordinates
c                           (or function values) of the nx data points
c                           (input).
c                nx     - number of elements in x and y (input).
c                           nx must be .ge. 2.
c                c      - spline coefficients (input). c is an nx-1 by
c                           3 matrix.
c                ic     - row dimension of matrix c exactly as
c                           specified in the dimension statement
c                           in the calling program (input).
c                           ic must be .ge. nx-1
c                u      - vector of length m containing the abscissae
c                           of the m points at which the cubic spline
c                           is to be evaluated (input).
c                s      - vector of length m (output).
c                           the value of the spline approximation at
c                           u(i) is
c                           s(i) = ((c(j,3)*d+c(j,2))*d+c(j,1))*d+y(j)
c                           where x(j) .le. u(i) .lt. x(j+1) and
c                           d = u(i)-x(j).
c                m      - number of elements in u and s (input).
c                ier    - error parameter (output).
c                         warning error
c                           ier = 33, u(i) is less than x(1).
c                           ier = 34, u(i) is greater than x(nx).
c
c                           ********************************************
c                           output of warning errors has been suppressed
c                           ********************************************
c
c   precision/hardware  - single and double/h32
c                       - single/h36,h48,h60
c
c   reqd. IMSL routines - uertst1,ugetio
c
c   notation            - information on special notation and
c                           conventions is available in the manual
c                           introduction or through IMSL routine uhelp
c
c   remarks  1.  the routine assumes that the abscissae of the nx
c                data points are ordered such that x(i) is less than
c                x(i+1) for i = 1,...,nx-1. no check of this condition
c                is made in the routine. unordered abscissae will cause
c                the algorithm to produce incorrect results.
c            2.  the routine generates two warning errors. one error
c                occurs if u(i) is less than x(1), for some i in the
c                the interval (1,m) inclusively. the other error occurs
c                if u(i) is greater than x(nx), for some i in the
c                interval (1,m) inclusively.
c            3.  the ordinate y(nx) is not used by the routine. for
c                u(k) .gt. x(nx-1), the value of the spline, s(k), is
c                given by
c                 s(k) = ((c(nx-1,3)*d+c(nx-1,2))*d+c(nx-1,1))*d+y(nx-1)
c                where d = u(k)-x(nx-1).
c
c   copyright           - 1978 by imsl, inc. all rights reserved.
c
c   warranty            - IMSL warrants only that IMSL testing has been
c                           applied to this code. no other warranty,
c                           expressed or implied, is applicable.
c
c ----------------------------------------------------------------------
c
c     specifications for arguments
c
      integer            nx,ic,m,ier
      real*8             x(nx),y(nx),c(ic,3),u(m),s(m)
c
c     specifications for local variables
c
      integer            i,jer,ker,nxm1,k
      real*8             d,dd,zero
      data               i/1/, zero/0.0/
c
c     first executable statement
c
      jer = 0
      ker = 0
      if (m .le. 0)  go to 9005
      nxm1 = nx-1
      if (i .gt. nxm1)  i = 1
c
c     evaluate spline at m points
c
      do 40 k=1,m
c
c        find the proper interval
c
         d = u(k)-x(i)
         if (d) 5, 25, 15
    5    if (i .eq. 1)  go to 30
         i = i-1
         d = u(k)-x(i)
         if (d) 5, 25, 20
   10    i = i+1
         d = dd
   15    if (i .ge. nx)  go to 35
         dd = u(k)-x(i+1)
         if (dd .ge. zero)  go to 10
         if ( d .eq. zero)  go to 25
c
c        perform evaluation
c
   20    s(k) = ((c(i,3)*d+c(i,2))*d+c(i,1))*d+y(i)
         go to 40
   25    s(k) = y(i)
         go to 40
c
c        u(k) < x(1)
c
   30    jer = 33
         go to 20
c
c        u(k) > x(nx)
c
   35    if (dd .gt. zero)  ker = 34
         d = u(k) - x(nxm1)
         i = nxm1
         go to 20
c
   40 continue
c
      ier = MAX0 (jer, ker)
c
****  if (jer .gt. 0)  call uertst1 (jer, 'icsevu1')
****  if (ker .gt. 0)  call uertst1 (ker, 'icsevu1')
c
 9005 return
c
      end

      subroutine icsicu1 (x, y, nx, bpar, c, ic, ier)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c --- modified version of IMSL subroutine named ICSICU
c
c ----------------------------------------------------------------------
c
c   computer            - dec10/single
c
c   latest revision     - january 1, 1978
c
c   purpose             - interpolatory approximation by cubic splines
c                           with arbitrary second derivative end
c                           conditions.
c
c   usage               - call icsicu1 (x, y, nx, bpar, c, ic, ier)
c
c   arguments    x      - vector of length nx containing the abscissae
c                           of the nx data points (x(i),y(i)) i = 1,...,
c                           nx. (input) x must be ordered so that
c                           x(i) .lt. x(i+1).
c                y      - vector of length nx containing the ordinates
c                           (or function values) of the nx data points.
c                           (input)
c                nx     - number of elements in x and y. (input) nx
c                           must be .ge. 2.
c                bpar   - vector of length 4 containing the end
c                           condition parameters. (input)
c                           2.0*spp(1)+bpar(1)*spp(2) = bpar(2),
c                           bpar(3)*spp(nx-1)+2.0*spp(nx) = bpar(4),
c                           where spp(i) = second derivative of the
c                           cubic spline function s evaluated at x(i).
c                c      - spline coefficients. (output) c is an nx-1 by
c                           3 matrix. the value of the spline
c                           approximation at t is
c                           s(t) = ((c(i,3)*d+c(i,2))*d+c(i,1))*d+y(i)
c                           where x(i) .le. t .lt. x(i+1) and
c                           d = t-x(i).
c                ic     - row dimension of matrix c exactly as
c                           specified in the dimension statement in
c                           the calling program. (input)
c                ier    - error parameter. (output)
c                         terminal error
c                           ier = 129, ic is less than nx-1
c                           ier = 130, nx is less than 2.
c                           ier = 131, input abscissa are not ordered
c                             so that x(1) .lt. x(2) ... .lt. x(nx).
c
c   precision/hardware  - single and double/h32
c                       - single/h36,h48,h60
c
c   reqd. IMSL routines - uertst1,ugetio
c
c   notation            - information on special notation and
c                           conventions is available in the manual
c                           introduction or through IMSL routine uhelp
c
c   copyright           - 1978 by imsl, inc. all rights reserved.
c
c   warranty            - IMSL warrants only that IMSL testing has been
c                           applied to this code. no other warranty,
c                           expressed or implied, is applicable.
c
c ----------------------------------------------------------------------
c
c     specifications for arguments
c
      integer            nx,ic,ier
      real*8             x(nx),y(nx),bpar(4),c(ic,3)
c
c     specifications for local variables
c
      integer            i,j,nxm1
      real*8             dx,dxj,dxjp1,dxp,dyj,dyjp1,half,one,pj,
     .                   six,sixi,two,yppa,yppb,zero
      equivalence        (dxj,yppb),(pj,sixi),(dxjp1,yppa)
      data               zero/0.0/,half/0.5/,one/1.0/,
     .                   two/2.0/,six/6.0/
c
      ier = 0
c
c     check error conditions
c
      nxm1 = nx-1
      if (ic .lt. nxm1)  go to 30
      if (nx .lt. 2   )  go to 35
      if (nx .eq. 2   )  go to 10
c
c     compute coefficients and right hand side of the tridiagonal
c     system defining the second derivatives of the spline interpolant for (x,y)
c
c     c(j,1) = lambda(j)
c     c(j,2) = mu(j)
c     c(j,3) = d(j)
c
      dxj = x(2)-x(1)
      if (dxj .le. zero)  go to 40
      dyj = y(2)-y(1)
      do 5 j=2,nxm1
         dxjp1 = x(j+1)-x(j)
         if (dxjp1 .le. zero)  go to 40
         dyjp1 = y(j+1)-y(j)
         dxp = dxj+dxjp1
         c(j,1) = dxjp1/dxp
         c(j,2) = one-c(j,1)
         c(j,3) = six*(dyjp1/dxjp1-dyj/dxj)/dxp
         dxj = dxjp1
         dyj = dyjp1
    5 continue
c
c     factor the tridiagonal matrix and solve for u
c
c     c(j,2)  = u(j)
c     c(j,1)  = q(j)
c     bpar(1) = lambda(1)
c     bpar(2) = d(1)
c     bpar(3) = mu(nx)
c     bpar(4) = d(nx)
c
   10 c(1,1) = -bpar(1)*half
      c(1,2) = bpar(2)*half
      if (nx .eq. 2)  go to 20
      do 15 j=2,nxm1
         pj = c(j,2)*c(j-1,1)+two
         c(j,1) = -c(j,1)/pj
         c(j,2) = (c(j,3)-c(j,2)*c(j-1,2))/pj
   15 continue
c
c     solve for cubic coefficients of spline interpolant
c     c(j,1), c(j,2), and c(j,3)
c
   20 yppb = (bpar(4)-bpar(3)*c(nxm1,2))/(bpar(3)*c(nxm1,1)+two)
      sixi = one/six
      do 25 i=1,nxm1
         j = nx-i
         yppa = c(j,1)*yppb+c(j,2)
         dx = x(j+1)-x(j)
         c(j,3) = sixi*(yppb-yppa)/dx
         c(j,2) = half*yppa
         c(j,1) = (y(j+1)-y(j))/dx-(c(j,2)+c(j,3)*dx)*dx
         yppb = yppa
   25 continue
      go to 9005
   30 ier = 129
      go to 9000
   35 ier = 130
      go to 9000
   40 ier = 131
c
 9000 call uertst1 (ier, 'icsicu1')
 9005 return
c
      end

      subroutine impsrc
c
      USE param
      USE solcon
      USE soln
      USE numbrs
      USE sourc
      USE geom
      USE tordlrot
      USE tmpcom
      USE fdyimp,                       ONLY : afdayi,bfdayi,cfdayi ! HSJ 12/28/08
c      USE fdyimp,                       ONLY : afdayi,bfdayi,cfdayi,
c    .                                          fday2d3,f2d3mult
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c this subroutine adds to the source those terms that
c are handled implicitly in the difference equations
c
c

c IMPSRC moves the term (1/hp)(d/drho)(-d*bp0) to the right hand (source) side
c
c
      onemt = (1.0 - theta) / theta
c
      do 10 j=1,nj
      qdelt(j) = qdimpl(j)*theta*((u(nk-iangrot-2,j)-u(nk-iangrot-1,j))
     .           + onemt*(usave(nk-iangrot-2,j)-usave(nk-iangrot-1,j)))
      s(nk-iangrot-2,j) = s(nk-iangrot-2,j) - qdelt(j)
   10 s(nk-iangrot-1,j) = s(nk-iangrot-1,j) + qdelt(j)
c
      if (codeid .eq. 'onedee')  return
c
      !fday sources  are not known  until a single time sep is taken
      IF(n > 0)THEN ! n (solcon.f90) cts transport time steps
         do j=2,nj-1
            fday2d3(j) = -(afdayi(j)*(u(nk-iangrot,j-1)+onemt*
     .                 usave(nk-iangrot,j-1))
     .      + bfdayi(j)*(u(nk-iangrot,j  )+onemt*usave(nk-iangrot,j  ))
     .      + cfdayi(j)*(u(nk-iangrot,j+1)+onemt*usave(nk-iangrot,j+1)))
            s(nk-iangrot,j) = s(nk-iangrot,j) + fday2d3(j)*f2d3mult
         end do
      ENDIF
      return
c
      end

      subroutine info (nout, nqik, ihead, n, extime, time, dt, dtt,
     .                 delav, delmax, kmax, jmax, delit, iter,
     .                 ineu, inub, irfcalc, imix, iconvg, isecrem0)
c
      USE param
      USE ions
      USE yoka
      USE solcon ,only : cpu_time
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c --- this subroutine writes out time step information
c
c
      character*8 np, namek
      dimension   irfcalc(krf)



c
      call TIMELEFT (isecrem1)
      extime = isecrem0 - isecrem1
c
      np     = '  '
      if (dtt .ne. dt)
     .np     = '.5'
c
      namek  = '  '
      if (kmax .ne. 0)
     .namek = nameu(kmax)
      ipart = iconvg - 1
c

      do nunit=nout,nqik
        if (nunit .ne. nqik .or. iyoka .eq. 0) then
          if (ihead .ne. 0)
     .    write (nunit, 50)
   50     format ('1',4x,'n',3x,'ex. time',3x,'tr. time',5x,'dt',
     .                8x,'delav',7x,'delmax',3x,'kmax', '  jmax',
     .                5x,'delit',4x,'iter',1x,'part',1x,'imix',
     .                1x,'ineu',1x,'inub', '  irf'      /
     .               12x, '(s)', 7x, '(ms)', 6x, '(ms)' /)
          write (nunit, 60) n, np(1:2), extime, time, dt, delav, delmax,
     .                      namek(1:2), jmax, delit, iter, ipart, imix,
     .                      ineu, inub, (irfcalc(k), k=1,9) ! really krf
   60     format (i6, a2, f8.2, 1x, 2f10.2, 1x, 2(1pe12.3), 3x, a2,
     .            i6, 1pe13.3, 5i5, 9i2)
        end if
      end do
c
      ihead = 0
      ineu  = 0
      inub  = 0
      do k=1,krf
        irfcalc(k) = 0
      end do
      return
c
      end

      subroutine interp1 (xval, x, n, y, yval)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c --- same as interp but:
c     a) if xval is .le. x(1) returns yval = y(1) (instead if 0.0)
c     b)            .ge. x(n)             =y(n)
c
c ----------------------------------------------------------------------
c
      dimension  x(*), y(*)
c
      if      (xval .le. x(1)) then
        yval = y(1)
      else if (xval .ge. x(n)) then
        yval = y(n)
      else
        call find (i1, i2, xval, x, n)
        if (i1 .eq. i2) then
          yval = y(i1)
        else
          yval = y(i1) + (y(i2)-y(i1)) * (xval-x(i1)) / (x(i2)-x(i1))
        end if
      end if
      return
c
      end

      subroutine interp (xval, x, n, y, yval)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      dimension  x(*), y(*)
c
      call find (i1, i2, xval, x, n)
      yval = 0.0
      if (i1 .eq. 0 )  return
      yval = y(i1)
      if (i1 .eq. i2)  return
      yval = y(i1) + (y(i2)-y(i1)) * (xval-x(i1)) / (x(i2)-x(i1))
      return
c
      end

      subroutine intrp (ilow1, ihigh1, x, y, nold, xnew, ynew, nnew)
c
c ----------------------------------------------------------------------
c  This subroutine uses the cubic spline routines to calculate
c  ynew as a function of xnew given y as a function of x.
c
c  The boundary conditions are specified as follows:
c
c       ilow: 0        set derivative at lower x value to zero
c                  -1  set derivative according to bpar
c         not 0 or -1  let derivative at lower x value be free
c
c      ihigh: 0        set derivative at upper x value to zero
c                  -1  set derivative according to bpar
c         not 0 or -1  let derivative at upper x value be free
c
c ----------------------------------------------------------------------
c
      USE param
      USE bd_condtn,only : bparenp,bpareni,bparte, bparti,
     .                     bparang, bparene,bparzeff, bparcur,
     .                     bparkpol, bparcurb_external,
     .                     iprofnbr, mtimeprf, knoterror,
     .                     profiles_bcondspl
      USE replace_imsl,    ONLY : my_icsicu,my_icsevu
      implicit  integer (i-n), real*8 (a-h, o-z)
c     include 'param.i'
c     include 'bcondspl.i'
c
      integer    ilow, ilow1, ihigh, ihigh1
      dimension  xnew(*), ynew(*), x(*), y(*)
      dimension  bpar(4)
      real *8 ,dimension(:),allocatable :: xt,yt
      real *8 ,dimension(:,:),allocatable :: c


      noldm1 = nold - 1
      allocate (c(noldm1,3),STAT = istat)
      if(istat .ne. 0)
     .          call allocate_error("c,intrp",0,istat) 
c
c     must ensure that x(1) .le. x(2) ...
c
      if (x(2) .le. x(1))  go to 10
      iord = 0
c
   60 do i=1,4
        bpar(i) = 0.0
      end do
c
      ilow  = ilow1
      ihigh = ihigh1
      if (ilow .eq. -1 .or. ihigh .eq. -1) then
        if (iprofnbr .le. kprim) then
          k = iprofnbr
          do 100 j=1,4
  100     bpar(j) = bparenp(j,k)
        else if ((iprofnbr .gt. kprim      ) .and.
     .           (iprofnbr .le. kprim+kimp)) then
          k = iprofnbr - kprim
          do 110 j=1,4
  110     bpar(j) = bpareni(j,k)
        else if (iprofnbr .eq. kprim+kimp+1) then
          do 120 j=1,4
  120     bpar(j) = bparene(j,mtimeprf)
        else if (iprofnbr .eq. kprim+kimp+2) then
          do 130 j=1,4
  130     bpar(j) = bparte(j,mtimeprf)
        else if (iprofnbr .eq. kprim+kimp+3) then
          do 140 j=1,4
  140     bpar(j) = bparti(j,mtimeprf)
        else if (iprofnbr .eq. kprim+kimp+4) then
          do 150 j=1,4
  150     bpar(j) = bparzeff(j,mtimeprf)
        else if (iprofnbr .eq. kprim+kimp+5) then
          do 160 j=1,4
  160     bpar(j) = bparcur(j,mtimeprf)
        else if (iprofnbr .eq. kprim+kimp+6) then
          do 170 j=1,4
  170     bpar(j) = bparang(j,mtimeprf)
        else if (iprofnbr .eq. kprim+kimp+7) then
          do 180 j=1,4
  180     bpar(j) = bparang(j,mtimeprf)
        else if (iprofnbr .eq. kprim+kimp+8) then
          do 190 j=1,4
  190     bpar(j) = bparcurb_external(j,mtimeprf)
        else
c          write  (n7, 45) iprofnbr
   45     format (' ERROR in spline input, iprofnbr =', i5)
          call STOP ('subroutine INTRP: spline input problem', 84)
        end if
        if (ilow .eq. -1) then
          if (bpar(1) .lt. -1.0e29) then
            bpar(1) = 1.0
            bpar(2) = 6.0*((y(2)-y(1))/(x(2)-x(1))-bpar(2))/(x(2)-x(1))
          else
            ilow    = 0
            bpar(1) = 1.0         !this sets gradient at r=0.0 to 0.0
            bpar(2) = 6.0 * (y(2)-y(1))/(x(2)-x(1))**2
          end if
        end if
        if (ihigh .eq. -1) then
          if (bpar(3) .le. -1.0e29) then
            bpar(3) = 1.0    !set gradient to value in bpar(4)
            bpar(4) = 6.0 * (bpar(4)-(y(nold)-y(nold-1)) /
     .                      (x(nold)-x(nold-1))) / (x(nold)-x(nold-1))
          else
            ihigh   = 1       !set gradient appropriate for natural cubic spline
            bpar(3) = 0.0
            bpar(4) = 0.0
          end if
        end if
        go to 8
      end if
c
      if (ilow .ne. 0)  go to 6
      bpar(1) =  1.0
      bpar(2) =  6.0*(y(2)-y(1))/(x(2)-x(1))**2
c
    6 if (ihigh .ne. 0)  go to 8
      bpar(3) =  1.0
      bpar(4) = -6.0*(y(nold)-y(nold-1))/(x(nold)-x(nold-1))**2
c
c    8 call icsicu1 (x, y, nold, bpar, c, noldm1, ier)
    8 call my_icsicu (x, y, nold, bpar, c, noldm1, ier)
      if (ier .ne. 0)  go to 20
c
c      call icsevu1 (x, y, nold, c, noldm1, xnew, ynew, nnew, ier)
      call my_icsevu(x, y, nold, c, noldm1, xnew, ynew, nnew, ier)
      deallocate (c, STAT = istat)
      if(istat .ne. 0)
     .          call deallocate_error("c,intrp",0,istat)
      if (iord .eq. 0)  return
c
      do i=1,nold
        x(i) = xt(i)
        y(i) = yt(i)
      end do
      if(allocated(xt))then
          deallocate (xt, STAT = istat)
          if(istat .ne. 0)
     .          call deallocate_error("xt,intrp",0,istat)
      endif
      if(allocated(yt))then
          deallocate (yt, STAT = istat)
          if(istat .ne. 0)
     .          call deallocate_error("yt,intrp",0,istat)
      endif
c
c --- normal return
c
      return
c
c --- handle ordering anomalies
c
   10 iord = 1
      allocate (xt(1:nold),STAT = istat)
      if(istat .ne. 0)
     .          call allocate_error("xt",0,istat) 
      allocate (yt(1:nold),STAT = istat)
      if(istat .ne. 0)
     .          call allocate_error("yt",0,istat) 
c
      do i=1,nold
        xt(i) = x(i)
        yt(i) = y(i)
      end do
c
      do i=1,nold
        x(i) = xt(nold+1-i)
        y(i) = yt(nold+1-i)
      end do
c
      go to 60
c
c --- abort, with a data dump, due to problem returning from IMSL
c
   20 write  (6, 40)  ier, iord, nold
c      write  (n7, 40)  ier, iord, nold
   40 format (/ ' INTRP   ier =', i4, '   iord =', i2, '   nold =', i3 /
     .         '    i        x           y')
      write  (6, '(i5, 2e12.4)')  (i, x(i), y(i), i=1,nold)
c      write  (n7, '(i5, 2e12.4)')  (i, x(i), y(i), i=1,nold)
      write  (6, 50)
   50 format (/ ' ERROR from subroutine INTRP:'                  /
     .          '       interpolation cannot proceed due to the' /
     .          '       non-monotonicity of the function above'  )
c      call giveupus(n7)
      call STOP ('subroutine INTRP: non-monotonic function', 85)
c
      end

      subroutine leqt1fl(a, m, n, ia, b, idgt, wkarea, ier)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c-leqt1f--------s/d-----library 2--------------------------------------
c
c   function            - linear equation solution - full storage
c                           mode - space economizer solution.
c   usage               - call leqt1f (a,m,n,ia,b,idgt,wkarea,ier)
c   parameters   a      - input matrix of dimension n by n containing
c                           the coefficient matrix of the equation
c                           ax = b.
c                         on output, a is replaced by the lu
c                           decomposition of a rowwise permutation of
c                           a.
c                m      - number of right-hand sides.(input)
c                n      - order of a and number of rows in b.(input)
c                ia     - number of rows in the dimension statement
c                           for a and b in the calling program. (input)
c                b      - input matrix of dimension n by m containing
c                           right-hand sides of the equation ax = b.
c                         on output, the n by m solution x replaces b.
c                idgt   - input option.
c                         if idgt is greater than 0, the elements of
c                           a and b are assumed to be correct to idgt
c                           decimal digits and the routine performs
c                           an accuracy test.
c                         if idgt equals zero, the accuracy test is
c                           bypassed.
c                wkarea - work area of dimension greater than or equal
c                           to n.
c                ier    - error parameter
c                         terminal error = 128+n.
c                           n = 1 indicates that a is algorithmically
c                             singular. (see the chapter l prelude).
c                         warning error = 32+n.
c                           n = 2 indicates that the accuracy test
c                                 failed.
c                                 the computed solution may be in error
c                                 by more than can be accounted for by
c                                 the uncertainty of the data.
c                                 this warning can be produced only if
c                                 idgt is greater than 0 on input.
c                                 see chapter l prelude for further
c                                 discussion.
c   precision           - single/double
c   req'd IMSL routines - ludatf1,luelmf1,uertst1
c   language            - fortran
c ----------------------------------------------------------------------
c   latest revision     - march 22,1974
c                                  dec
c
      dimension          a(ia,*),b(ia,*),wkarea(*)
****  double precision   a,b,wkarea,d1,d2,wa
c
      ier = 0
c                                  decompose a
      call ludatf1 (a,a,n,ia,idgt,d1,d2,wkarea,wkarea,wa,ier)
      if (ier .gt. 128)  go to 9000
c                                  call routine luelmf1 (forward and
c                                  backward substitutions)
      do 10 j=1,m
         call luelmf1 (a,b(1,j),wkarea,n,ia,b(1,j))
   10 continue
      if (ier .eq. 0)  go to 9005
 9000 continue
      call uertst1 (ier,'leqt1fl')
 9005 return
c
      end

      subroutine linav (ene, nj, nx, ny, p, plim, psir, rmhdgrid,
     .                  zmhdgrid, zax, enebar)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c This subroutine calculates the line-average density along a horizontal
c chord passing through the magnetic axis and the toroidal axis.
c Any contribution from the scrape-off layer is neglected.
c
      dimension  p(nx,*), rmhdgrid(*), zmhdgrid(*), ene(*), psir(*)
      data nc    / 101   /
      data plimo /-1.0e10/
c
      enebar = 0.0
      path   = 0.0
c
c determine inner and outer radii for line integral
c
      rin   = 0.0
!      if (plim .eq. plimo)  go to 80
      plimo = plim
      dx    = rmhdgrid(2)-rmhdgrid(1)
      i     = 1
      call bilin (nx, ny, rmhdgrid, zmhdgrid, p, rmhdgrid(i), zax, psi1)
      ii    = 1

      do i=2,nx
        call bilin(nx,ny,rmhdgrid,zmhdgrid,p,rmhdgrid(i),zax,psi2)
        if (  ii .eq. 2   )  go to 30
        if (psi2 .gt. plim)  go to 50
        rin = rmhdgrid(i-1) + dx*(psi1-plim)/(psi1-psi2)
        ii  = 2
        go to 50
   30   if (psi2 .lt. plim)  go to 50
        rout = rmhdgrid(i-1) + dx*(plim-psi1)/(psi2-psi1)
        go to 80
   50   psi1 = psi2
      end do
c
      rout = rmhdgrid(nx)
   80 if (rin .eq. 0.0)  return
c
c evaluate line integral
c
      dr = (rout-rin)/(nc-1.0)
c
      do i=1,nc
        wt   = 0.5
        ene1 = ene(nj)
        if (i .eq. 1 .or. i .eq. nc)  go to 120
        wt   = 1.0
        r1   = rin + (i-1) * dr
        call bilin  (nx, ny, rmhdgrid, zmhdgrid, p, r1, zax, psi1)
        call interp (psi1, psir, nj, ene, ene1)
  120   enebar = enebar + wt*ene1
        path = path + wt
      end do
c
      enebar = enebar / path

      return
c
      end

      subroutine ludatf1 (a, lu, n, ia, idgt, d1, d2, ipvt,
     .                    equil, wa, ier)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c-ludatf1--------s/d-----library 2--------------------------------------
c
c   function            - l-u decomposition by the crout algorithm
c                           with optional accuracy test.
c   usage               - call ludatf1 (a,lu,n,ia,idgt,d1,d2,ipvt,
c                                       equil,wa,ier)
c   parameters   a      - input matrix of dimension n by n containing
c                           the matrix to be decomposed
c                lu     - real output matrix of dimension n by n
c                           containing the l-u decomposition of a
c                           rowwise permutation of the input matrix.
c                           for a description of the format of lu, see
c                           example.
c                n      - input scalar containing the order of the
c                           matrix a.
c                ia     - input scalar containing the row dimension of
c                           matrices a and lu in the calling program.
c                idgt   - input option.
c                           if idgt is greater than zero, the non-zero
c                           elements of a are assumed to be correct to
c                           idgt decimal places.  ludatf1 performs an
c                           accuracy test to determine if the computed
c                           decomposition is the exact decomposition
c                           of a matrix which differs from the given on
c                           by less than its uncertainty.
c                         if idgt is equal to zero, the accuracy test i
c                           bypassed.
c                d1     - output scalar containing one of the two
c                           components of the determinant. see
c                           description of parameter d2, below.
c                d2     - output scalar containing one of the
c                           two components of the determinant. the
c                           determinant may be evaluated as (d1)(2**d2)
c                ipvt   - output vector of length n containing the
c                           permutation indices. see document
c                           (algorithm).
c                equil  - output vector of length n containing
c                           reciprocals of the absolute values of
c                           the largest (in absolute value) element
c                           in each row.
c                wa     - accuracy test parameter, output only if
c                           idgt is greater than zero.
c                           see element documentation for details.
c                ier    - error parameter
c                         terminal error = 128+n
c                           n = 1 indicates that matrix a is
c                                 algorithmically singular. (see the
c                                 chapter l prelude).
c                         warning error = 32+n
c                           n = 2 indicates that the accuracy test
c                                 failed.
c                                 the computed solution may be in error
c                                 by more than can be accounted for by
c                                 the uncertainty of the data.
c                                 this warning can be produced only if
c                                 idgt is greater than 0 on input.
c                                 see chapter l prelude for further
c                                 discussion.
c   precision           - single/double
c   req'd IMSL routines - uertst1
c   language            - fortran
c ----------------------------------------------------------------------
c
c   latest revision     - march 22,1974 (dec)
c
      dimension a(ia,*), lu(ia,*), ipvt(*), equil(*)
      real*8    lu
      data      zero, one, four, sixtn, sixth
     .         /0.0 , 1.0, 4.0 , 16.0 , 0.0625/
c
      ier  = 0
      rn   = n
      wrel = zero
      d1   = one
      d2   = zero
      biga = zero
      do 10 i=1,n
         bigg = zero
         do 5 j=1,n
            p       = a(i,j)
            lu(i,j) = p
****        p       = DABS (p)
            p       =  ABS (p)
            if (p .gt. bigg)  bigg = p
    5    continue
         if (bigg .gt.  biga) biga = bigg
         if (bigg .eq.  zero)  go to 110
         equil(i) = one/bigg
   10 continue
      do 105 j=1,n
         jm1 = j-1
         if (jm1 .lt. 1)  go to 40
c
c compute u(i,j), i = 1,...,j-1
c
         do 35 i=1,jm1
            sum = lu(i,j)
            im1 = i-1
            if (idgt .eq. 0)  go to 25
c
c with accuracy test
c
            ai = ABS (sum)
            wi = zero
            if (im1 .lt. 1)  go to 20
            do 15 k=1,im1
               t = lu(i,k)*lu(k,j)
               sum = sum-t
               wi  = wi + ABS (t)
   15       continue
            lu(i,j) = sum
   20       wi = wi + ABS (sum)
            if (ai .eq. zero) ai = biga
            test = wi/ai
            if (test .gt. wrel) wrel = test
            go to 35
c
c without accuracy
c
   25       if (im1 .lt. 1)  go to 35
            do 30 k=1,im1
               sum = sum-lu(i,k)*lu(k,j)
   30       continue
            lu(i,j) = sum
   35    continue
   40    p = zero
c
c compute u(j,j) and l(i,j), i = j+1,...
c
         do 70 i=j,n
            sum = lu(i,j)
            if (idgt .eq. 0)  go to 55
c
c with accuracy test
c
            ai = ABS (sum)
            wi = zero
            if (jm1 .lt. 1)  go to 50
            do 45 k=1,jm1
               t = lu(i,k)*lu(k,j)
               sum = sum-t
               wi = wi + ABS (t)
   45       continue
            lu(i,j) = sum
   50       wi = wi + ABS (sum)
            if (ai .eq. zero) ai = biga
            test = wi/ai
            if (test .gt. wrel) wrel = test
            go to 65
c
c without accuracy test
c
   55       if (jm1 .lt. 1)  go to 65
            do 60 k=1,jm1
               sum = sum-lu(i,k)*lu(k,j)
   60       continue
            lu(i,j) = sum
   65       q = equil(i) * ABS (sum)
            if (p .ge. q)  go to 70
            p = q
            imax = i
   70    continue
c
c test for algorithmic singularity
c
         q = rn+p
         if (q .eq. rn)  go to 110
         if (j .eq. imax)  go to 80
c
c interchange rows j and imax
c
         d1 = -d1
         do 75 k=1,n
            p = lu(imax,k)
            lu(imax,k) = lu(j,k)
            lu(j,k) = p
   75    continue
         equil(imax) = equil(j)
   80    ipvt(j) = imax
         d1 = d1*lu(j,j)
   85    if (ABS (d1) .le. one)  go to 90
         d1 = d1*sixth
         d2 = d2+four
         go to 85
   90    if (ABS (d1) .ge. sixth)  go to 95
         d1 = d1*sixtn
         d2 = d2-four
         go to 90
   95    continue
         jp1 = j+1
         if (jp1 .gt. n)  go to 105
c
c divide by pivot element u(j,j)
c
         p = lu(j,j)
         do 100 i=jp1,n
            lu(i,j) = lu(i,j)/p
  100    continue
  105 continue
c
c perform accuracy test
c
      if (idgt .eq. 0)  go to 9005
      p = 3*n+3
      wa = p*wrel
****  q = wa+10.0d0**(-idgt)
      q = wa+10.0**(-idgt)
      if (q .ne. wa)  go to 9005
      ier = 34
      go to 9000
c
c algorithmic singularity
c
  110 ier = 129
      d1 = zero
      d2 = zero
c
c print error
c
 9000 call uertst1 (ier, 'ludatf1')
 9005 return
c
      end

      subroutine luelmf1 (a, b, ipvt, n, ia, x)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c-luelmf1--------s/d-----library 2--------------------------------------
c
c   function            - elimination part of solution of ax = b -
c                           full storage mode
c   usage               - call luelmf1 (a,b,ipvt,n,ia,x)
c   parameters   a      - the result, lu, computed in the subroutine
c                           'ludatf1', where l is a lower triangular
c                           matrix with ones on the main diagonal. u is
c                           upper triangular. l and u are stored as a
c                           single matrix a, and the unit diagonal of
c                           l is not stored
c                b      - b is a vector of length n on the right hand
c                           side of the equation ax = b
c                ipvt   - the permutation matrix returned from the
c                           subroutine 'ludatf1', stored as an n length
c                           vector
c                n      - order of a and number of rows in b
c                ia     - number of rows in the dimension statement
c                           for a in the calling program.
c                x      - the result x
c   precision           - single/double
c   language            - fortran
c ----------------------------------------------------------------------
c   latest revision     - april 11,1975
c                                  dec
c
      dimension          a(ia,*),b(*),ipvt(*),x(*)
****  double precision   a,b,x,sum
c                                  solve ly = b for y
      do 5 i=1,n
    5 x(i) = b(i)
      iw = 0
      do 20 i=1,n
         ip = ipvt(i)
         sum = x(ip)
         x(ip) = x(i)
         if (iw .eq. 0)  go to 15
         im1 = i-1
         do 10 j=iw,im1
            sum = sum-a(i,j)*x(j)
   10    continue
         go to 20
   15    if (sum .ne. 0.0) iw = i
   20 x(i) = sum
c                                  solve ux = y for x
      do 30 ib=1,n
         i = n+1-ib
         ip1 = i+1
         sum = x(i)
         if (ip1 .gt. n)  go to 30
         do 25 j=ip1,n
            sum = sum-a(i,j)*x(j)
   25   continue
   30 x(i) = sum/a(i,i)
      return
c
      end

      subroutine makpro (x, y, nj, y0, ynj, yalp, ygam)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      dimension  x(*), y(*)
c
      xnji   = 1.0 / x(nj)
      do j=1,nj
        xsc  = 1.0 - (x(j)*xnji)**ygam
        if (xsc .lt. 0.0)
     .  xsc  = 0.0
        y(j) = ynj + (y0-ynj) * xsc**yalp
      end do
      return
c
      end

      subroutine maxrc (xnew, xold, nj, delmax)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c     MAXRC computes the maximum relative change in the array xold
c     as compared with the array xnew.
c     For xold(j) = 0 and xnew(j)<>0, del is set to 1
c
      dimension  xnew(*), xold(*)
c
      delmax = 0.0
      do j=1,nj
        if (xold(j) .eq. 0.0)  go to 100
        del = ABS (xnew(j)-xold(j))/xold(j)
        go to 110
  100   del = 0.0
        if (xnew(j) .ne. 0.0)  del = 1.0
  110   delmax = MAX (del, delmax)
      end do
      return
c
      end

      subroutine maxu (u, itran, nk, nj, kk, x, iangrot)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c this subroutine calculates the maximum value of each dependent
c variable for which transport is considered.
c these maxima are used as scale factors in solve.
c
      dimension  u(kk,*), itran(*), x(*)
c
      kt = 0
      do 20 k=1,nk -iangrot
      if (itran(k) .le. 0)  go to 20
      kt = kt + 1
      x(kt) = 0.0
      do 10 j=1,nj
      if (u(k,j) .lt. x(kt))  go to 10
      x(kt) = u(k,j)
   10 continue
   20 continue
c
c --- the toroidal rotation profile can be negative:
c
      if (iangrot .ne. 0 .and. itran(nk) .eq. 1) then
        kt       = kt + 1
        omegamin =  1.0d50
        omegamax = -1.0d50
        do j=1,nj
          omegamin = MIN (omegamin,u(nk,j))
          omegamax = MAX (omegamax,u(nk,j))
        end do
        ominabs = ABS (omegamin)
        omaxabs = ABS (omegamax)
        if (ominabs .gt. omaxabs) then
          x(kt) = omegamin
        else
          x(kt) = omegamax
        end if
      end if
      return
c
      end

      subroutine mescon (y, dr, nj)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c     this subroutine converts a mesh-centered array to a mesh-point
c     array. linear interpolation and extrapolation are used.
c     dr(j) = r(j+1)-r(j), j=1...nj-1
c
      dimension y(*), dr(*)
c
      x      = (dr(nj-2)+2.0*dr(nj-1)) / (dr(nj-2)+dr(nj-1))
      y(nj)  = y(nj-2) + x*(y(nj-1)-y(nj-2))
      do j=nj-1,2,-1
        x    = dr(j-1) / (dr(j-1)+dr(j))
        y(j) = y(j-1) + x*(y(j)-y(j-1))
      end do
      y(1)   = y(1) - (y(2)-y(1))
      return
c
      end

      subroutine neointrp (x0, x1, j, xav)
c
      USE param
       USE solcon
      USE mhdpar
      USE soln2d
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c      include 'param.i'
c      include 'mhdpar.i'
c      include 'solcon.i'
c      include 'soln2d.i'
c
      dimension  x0(*), x1(*)
c
      x0av = 0.5 * (x0(j) + x0(j+1))
      xav  = x0av
      if (itre .ne. 1) then
        x1av = 0.5 * (x1(j) + x1(j+1))
        xav  = x0av + (time-eqtim0) * (x1av-x0av) / dteq
      end if
      return
c
      end




      subroutine Onetwo_pedestal_driver
c  -----------------------------------------------------------------
c      subroutine sets up input for pedestal code, spwans pedestal,
c      reads the output file created by pedestal and passes the info
c      onto Onetwo proper.        HSJ 6/22/04
c  ------------------------------------------------------------------

c  Input:
c  ------
c lbound(j)  integer control array --- see more details below 
c cbound(j)  general control array --- see more details below
c
c real scalars:
c rminor     minor radius [m]
c rmajor     major radius [m]
c kappa      elongation
c delta      triangularity
c current    plasma current [MA]
c btor       toroidal field at rmajor [Tesla]
c nebar      line averaged electron density [particles/m^3]
c hydmass    average hydrogenic mass [AMU]
c zeff       effective charge at the edge of plasma
c ploss      power crossing separatrix [MW]
c
c Note: not all variable are needed for each model
c
c  Output:
c  -------
c
c real scalars:
c pth       power threshold of the transition from L-mode to H-mode [MW]
c neped     electron density at the top of the pedestal [particles/m^3]
c nGr       Greenwald density [particles/m^3]
c n_ratio   ratio of pedestal density to Greenwald density
c tped      temperature at the top of the pedestal [keV]
c dpdr_c    critical pressure gradient [Pa/m]
c width     width of the pedestal [m]
c
c integer scalars
c mode      plasma operation mode --- see more details below
c iter_t is number of iterations used in temperature calculation
c iter_s is number of iterations used in magnetic shear calculation 
c iter_q is number of iterations used in safety factor calculation 
c
c  Internal control variables:
c  ---------------------------
c lbound(j), j=1,32   integer array of control variables: 
c
c lbound(1) controls which L-H transition model is used
c                  Default lbound(1) = 1
c             = 1  L-H transtion model based on power threadhold
c lbound(2) controls which pedestal density model is used
c                  Default lbound(2) = 1
c             = 1  Use simple empirical pedestal density model
c lbound(3) controls which pedestal temperature model is used
c                  Default lbound(3) = 1
c             = 1  Use pedestal temperature model based on magnetic shear 
c                  and flow shear stabilization
c             = 2  Use pedestal temperature model based on flow shear 
c                  stabilization
c             = 3  Use pedestal temperature model based on normalized 
c                  pressure
c             = 11 Use pedestal temperature model based on the first 
c                  conduction model
c             = 12 Use pedestal temperature model based on the second 
c                  conduction model
c lbound(4) controls maximum number of iteration allowed for temperature 
c           calculation
c                  Dedault lbound(4) = 1000
c lbound(5) controls maximum number of iteration allowed for magnetic shear 
c           calculation
c                  Dedault lbound(5) = 1000
c lbound(6) controls maximum number of iteration allowed for safety factor 
c           calculation
c                  Dedault lbound(6) = 1000
c
c
c cbound(j), j=1,32   general array of control variables:
c
c cbound(1) error limit in the iterations in pedestal temperature calculation
c                  Default cbound(1) = 0.001
c cbound(2) error limit in the iterations in magnetic shear calculation
c                  Default cbound(2) = 0.001
c cbound(3) error limit in the iterations in safety factor calculation
c                  Default cbound(3) = 0.001
c
c  Output Integer:
c  --------------
c mode is a plasma operation mode
c              = 0 when plasma is in L-mode
c              = 1 when plasma is in H-mode
c
c ierr(1) is an input error flag
c             = 1  error in a minor radius input
c             = 2  error in a major radius input
c             = 3  error in an elongation input
c             = 4  error in a plasma current input
c             = 5  error in a toroidal field input
c             = 6  error in a line averaged density input
c             = 7  error in an averaged hydrogenic mass input
c             = 8  error in an effective charge input
c             = 9  error in a power cross separatrix input
c
c ierr(2) is an iteration error flag
c             = 1 eror in temperature iteration
c             = 2 error in shear iteration
c             = 3 error in safety factor iteration
c
c
      USE bd_condtn, only : lbound,cbound,totcur,fluxb,
     .                      pedestal_models,ped_nebar,ped_temp,
     .                      ped_Ngr,ped_grad,ped_mode, ped_nratio

      USE machin,    only : rminor_cm => rminor,             
     .                      kappa,rmajor_cm => rmajor,
     .                      btor_gauss=> btor
      USE psig,      only : triangnpsi
      USE extra,     only : enebar_cm3 => enebar
      USE numbrs,    ONLY : nprim,nj,nion
      USE ions,      ONLY : atw,zeff_12  => zeff
      USE soln ,     ONLY : en, ene
      USE mesh ,     ONLY : r
      USE io,        ONLY : ncrt,nout
      USE sourc,     ONLY : qrad
      USE fusion,    ONLY : thermal_thermal_dtntot,beam_thermal_dtntot,
     .                      beam_beam_dtntot
      USE constnts,  ONLY : pisq
      USE geom,      ONLY : hcap
      USE ext_prog_info, only: get_pedestal,pedestal_to_run


      IMPLICIT NONE
c
      integer j,jj,len_str,ishell,
     &  nin,
     &  mode,        iter_t,   iter_s,   iter_q,
     &  ierr(2),     i , irun_pedestal
c

c
      real*8  ele_fluxe,ion_fluxe,sfarea,pradl,qnloss,
     &  rmajor,      rminor,        btor,const,
     &  delta,       current,       nebar,
     &  hydmass,     zeff,     ploss,    pth,
     &  neped,       nGr,      n_ratio,  tped,
     &  dpdr_c,   width,tot_den, amassth

      character *72 runid_pedestal 
      logical exists
c     
       namelist /nt/ 
     &  lbound,     cbound,     rminor,     rmajor,
     &  kappa,      delta,      current,    btor,
     &  nebar,      hydmass,    zeff,       ploss 

        !enebar is calculated in fiziks. We cant do the
        !pedestal calculations until fizics has ben called:
        if(enebar_cm3 .lt. 1.e4) return

 



        lbound(1) = pedestal_models(1)
        lbound(2) = pedestal_models(2)
        lbound(3) = pedestal_models(3)
        runid_pedestal ='PEDESTAL input file created by ONETWO'  
        current = totcur(1)/1.e6
        delta = triangnpsi(1)         !1 = plasma edge
        btor =btor_gauss/1.e4
        rminor = rminor_cm/1.e2
        rmajor = rmajor_cm/1.e2
        nebar =enebar_cm3*1.e6
        amassth = 0.0d0
        tot_den = 0.0d0
        zeff =0.0
        do jj=1,nprim                     ! sum over hydrogenic species
           do j=1,nj
             zeff = zeff_12(j)*en(j,jj) + zeff
             if (atw(jj) .le. 3.1) then
                amassth = amassth+atw(jj)*en(j,jj)
                tot_den = tot_den+en(j,jj)
             end if
           enddo
        enddo
        zeff =    zeff/tot_den
        hydmass = amassth/tot_den


c       try to get an estimate of ploss. Note that none of these quantities
c       are known until we call the relevant routines. Hence some
c       startup difficulty is expected here:
c       ! electron and ion conduction and convection through the
c       !plasma boundary:
        ele_fluxe = MAX(0.0D0,fluxb(nion+1))
        ion_fluxe = MAX(0.0d0,fluxb(nion+2))
c       surface area of torus, (not the true area, but the
c       area consistent with the flux surface average calcualtions)
        const  = 4.0 * pisq*rmajor_cm*1.6e-22  
        sfarea = const* hcap(nj) * r(nj)
        ploss = sfarea*(ele_fluxe +ion_fluxe) ! MW note that this
                                              ! value has a considerable
                                              !uncertainty due to edge
                                              !effects 
c       add radiative loss:       
        call trapv (r,qrad,hcap,nj,pradl)
        pradl    = pradl*const
        ploss = ploss + pradl 
c       assumme all neutrons are lost (only dt reactions are counted)
        qnloss = (thermal_thermal_dtntot+ beam_thermal_dtntot
     .    + beam_beam_dtntot)*14100.0      ! 14.1  MeV neutrons for d(t,n)he4
        qnloss = qnloss*1.602e-22 ! MW
        ploss = ploss + qnloss  
        print *,thermal_thermal_dtntot,beam_thermal_dtntot,
     .                                     beam_beam_dtntot
        print *,'ploss ,qnloss,pradl,pcond =',
     .  ploss,qnloss,pradl,sfarea*(ele_fluxe +ion_fluxe)





c      original default values (not appropriate for Onetwo):
c      rminor      = 0.910E+00   ! Minor radius [m]
c      rmajor      = 2.920E+00   ! Major radius [m]
c      kappa       = 1.610E+00   ! Elongation
c      delta       = 0.170E+00   ! Triangluarity
c      current     = 2.550E+00   ! Plasma current [MA]
c      btor        = 2.350E+00   ! Toroidal field [T]
c      nebar       = 6.000E+19   ! Line averaged density [particles/m^-3] 
c      hydmass     = 2.000E+00   ! Averaged hydrogenic mass [AMU]
c      zeff        = 2.000E+00   ! Effective charge
c      ploss       = 1.500E+01   ! Power across separatrix [MW]


       nin = 91 
       call getioun(nin,nin)
       exists = .false.

c   destroy any existing input_1 file:
       INQUIRE(FILE = 'input_1', EXIST = exists)
       IF( exists ) CALL DESTROY ('input_1')

c   create input namelist file input_1:
       open(nin,  file='input_1',  status='new')



c      write the namelist file:
       WRITE (unit = nin, fmt = '(3x, a)') runid_pedestal
       WRITE (unit = nin, nml = nt)

       CLOSE (unit = nin)
       CALL giveupus(nin)
       


c      get pedestal code path  info for machine we are running on :
       call get_pedestal(ncrt,nout,pedestal_to_run,len_str)

c      spawn pedestal code:
       write (ncrt, '(   / a)')  ' spawning PEDESTAL CODE ...' 
       write (ncrt, '(a)')pedestal_to_run(1:len_str)
        if (ISHELL (pedestal_to_run(1:len_str)) .ne. 0)
     .  call STOP ('pedestal driver: failure of spawned PEDESTAL',1)
        write (ncrt, '(53x, a)')
     .       '...returned from PEDESTAL'
c
c
c
c      read the output of Pedestal:
       call getioun(nin,nin)
       open(nin, file='output_12', status='unknown')
       read(nin,FMT ='(4(2x,i3))')ierr(1),ierr(2),iter_t,mode
       if(ierr(1) + ierr(2) .eq. 0)then
          read(nin,FMT ='(3(2x,1pe12.6))')neped,nGr,n_ratio
          read(nin,FMT ='(3(2x,1pe12.6))')tped,width,dpdr_c
          !Set boundary condtions for Onetwo:
          ped_nebar= neped/1.e6     !1/cm^3
          ped_temp = tped
          ped_nGr  = nGr
          ped_nratio = n_ratio
          ped_grad = dpdr_c*6.24e7  ! kev/cm^4
          ped_mode = mode
       ELSE
          write(ncrt,
     . FMT='( "--------------ERROR in PEDESTAL MODULE-------------",/,
     . "   ierr(1),ierr(2) = ",/,
     . "   RESULTS ARE NOT USED")') ierr(1),ierr(2)
       ENDIF
       CLOSE (unit = nin)
       CALL giveupus(nin)
       
      return
      end
c










      SUBROUTINE piksr2(n,arr,brr)
      INTEGER n
      REAL*8 arr(n),brr(n)
      INTEGER i,j
      REAL*8  a,b
      do  j=2,n
        a=arr(j)
        b=brr(j)
        do  i=j-1,1,-1
          if(arr(i).le.a)goto 10
          arr(i+1)=arr(i)
          brr(i+1)=brr(i)
        enddo
        i=0
10      arr(i+1)=a
        brr(i+1)=b
      enddo
      return
      END

      subroutine print_matrix(a,nr,nc,kk)
      dimension a(kk,*)

      do  j=1,nr
         write(*,FMT='(2(1pe25.14,x))')(a(j,i),i=1,nc)
      enddo

      return
      end

      subroutine radfit (prad, eni, te, namei, nimp, nj, kj, ncrt, nout)
c
      USE numbrs,                  ONLY : nprim
      USE rad_loss,                ONLY : brems_nions
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c
c     this subroutine computes prad, the impurity radiation per
c     electron in watts. coeffficents for the curve fits used were
c     obtained from a Princeton PPL paper, "steady state radiative
c     cooling rates for low-density high-temperature plasmas"
c     (pppl-1352)
c
c     changed logic to use largest energy range in the table
c     instead of stopping code with error message.  HSJ 09/11/06
c
c     name(i)    = chemical abbreviation for element
c     coeff(i,1) = lower temperature range
c     coeff(i,2) = upper temperature range
c     coeff(i,3) ---> coeff(i,8) are the coefficients of a 5th degree
c
      parameter       (krows = 49)
      character*8 name,                      namei
      dimension   name(krows), coeff(krows,8)
      dimension   prad(*), eni(kj,*), te(*), namei(*)

c
c ----------------------------------------------------------------------
c
c *********  argon             **********
      data name( 1) /'ar'/, (coeff( 1,i),i = 1,8)/
     .3.000e-02,2.000e-01,-2.05304e+01,-2.83429e+00,
     . 1.50690e+01, 3.51718e+01, 2.40012e+01, 5.07272e+00/
      data name( 2) /'ar'/, (coeff( 2,i),i = 1,8)/
     .2.000e-01,2.000e+00,-1.96520e+01,-1.17276e-01,
     . 7.83322e+00,-6.35158e+00,-3.05885e+01,-1.52853e+01/
      data name( 3) /'ar'/, (coeff( 3,i),i = 1,8)/
     .2.000e+00,2.000e+01,-1.97488e+01, 2.96484e+00,
     .-8.82939e+00, 9.79100e+00,-4.96002e+00, 9.82003e-01/
      data name( 4) /'ar'/, (coeff( 4,i),i = 1,8)/
     .2.000e+01,1.000e+02,-2.11794e+01, 5.19148e+00,
     .-7.43972e+00, 4.96902e+00,-1.55318e+00, 1.87705e-01/
c *********  carbon            **********
      data name( 5) /'c '/, (coeff( 5,i),i = 1,8)/
     .2.000e-03,2.000e-02, 2.26196e+03, 5.27071e+03,
     . 4.80768e+03, 2.16700e+03, 4.83472e+02, 4.27841e+01/
      data name( 6) /'c '/, (coeff( 6,i),i = 1,8)/
     .2.000e-02,2.000e-01, 5.01266e+01, 3.32680e+02,
     . 6.00316e+02, 5.17230e+02, 2.13035e+02, 3.36380e+01/
      data name( 7) /'c '/, (coeff( 7,i),i = 1,8)/
     .2.000e-01,2.000e+00,-2.12371e+01,-3.13175e-01,
     . 7.60470e-01,-3.01657e-01, 7.63153e-02, 1.40114e-01/
      data name( 8) /'c '/, (coeff( 8,i),i = 1,8)/
     .2.000e+00,2.000e+01,-2.12367e+01,-3.44789e-01,
     . 1.03546e+00,-1.01249e+00, 5.92047e-01,-1.43559e-01/
      data name( 9) /'c '/, (coeff( 9,i),i = 1,8)/
     .  2.00000e+01, 1.00000e+04, -2.47680e+01, 9.40818e+00,
     . -9.65745e+00, 4.99916e+00, -1.23738e+00, 1.16061e-01/
c *********  chromium          **********
      data name(10) /'cr'/, (coeff(10,i),i = 1,8)/
     .2.000e-02,2.000e-01,-1.04662e+01, 4.84777e+01,
     . 1.03531e+02, 9.96556e+01, 4.53392e+01, 7.96357e+00/
      data name(11) /'cr'/, (coeff(11,i),i = 1,8)/
     .2.000e-01,2.000e+00,-1.85628e+01,-2.86143e+00,
     .-7.07898e+00, 1.04095e+01, 4.27995e+01, 2.95598e+01/
      data name(12) /'cr'/, (coeff(12,i),i = 1,8)/
     .2.000e+00,2.000e+01,-1.70389e+01,-1.81938e+01,
     . 5.07803e+01,-6.45257e+01, 3.81862e+01,-8.58752e+00/
      data name(13) /'cr'/, (coeff(13,i),i = 1,8)/
     .2.000e+01,1.000e+02,-1.41239e+01,-1.39958e+01,
     . 1.49161e+01,-8.16669e+00, 2.30136e+00,-2.62443e-01/
c *********  iron              **********
      data name(14) /'fe'/, (coeff(14,i),i = 1,8)/
     .2.000e-03,2.000e-02, 5.40647e+02, 1.31358e+03,
     . 1.22549e+03, 5.67632e+02, 1.30545e+02, 1.19286e+01/
      data name(15) /'fe'/, (coeff(15,i),i = 1,8)/
     .2.000e-02,2.000e-01,-1.19062e+01, 2.59325e+01,
     . 4.12580e+01, 3.03780e+01, 1.04477e+01, 1.39669e+00/
      data name(16) /'fe'/, (coeff(16,i),i = 1,8)/
     .2.000e-01,2.000e+00,-1.82996e+01,-1.55942e+00,
     .-6.00441e+00,-1.19285e-01, 1.70088e+01, 1.17635e+01/
      data name(17) /'fe'/, (coeff(17,i),i = 1,8)/
     .2.000e+00,2.000e+01,-1.67982e+01,-1.59773e+01,
     . 3.82500e+01,-4.25281e+01, 2.22674e+01,-4.46320e+00/
      data name(18) /'fe'/, (coeff(18,i),i = 1,8)/
     .2.000e+01,1.000e+02,-2.45396e+01, 1.79522e+01,
     .-2.35636e+01, 1.48450e+01,-4.54232e+00, 5.47746e-01/
c *********  helium            **********
      data name(19) /'he'/, (coeff(19,i),i = 1,8)/
     .2.000e-03,1.000e-02, 3.84322e+03, 8.93072e+03,
     . 8.17947e+03, 3.71287e+03, 8.35739e+02, 7.46792e+01/
      data name(20) /'he'/, (coeff(20,i),i = 1,8)/
     .1.000e-02,2.000e-01,-2.25831e+01, 1.15730e-02,
     .-8.32355e-01,-1.17916e+00,-4.74033e-01,-8.64483e-02/
      data name(21) /'he'/, (coeff(21,i),i = 1,8)/
     .2.000e-01,2.000e+00,-2.25560e+01, 3.28276e-01,
     . 1.39226e-01,-1.22085e-01,-2.76602e-01,-2.90494e-01/
      data name(22) /'he'/, (coeff(22,i),i = 1,8)/
     .2.000e+00,2.000e+01,-2.25798e+01, 4.57017e-01,
     .-1.59274e-01, 2.71952e-01,-1.71810e-01, 3.86609e-02/
      data name(23) /'he'/, (coeff(23,i),i = 1,8)/
     .2.000e+01,1.000e+02,-1.73046e+01,-1.62462e+01,
     . 2.10079e+01,-1.31207e+01, 4.06935e+00,-5.00944e-01/
c *********  krypton           **********
      data name(24) /'kr'/, (coeff(24,i),i = 1,8)/
     .5.000e-02,2.000e-01, 5.67564e+01, 4.00852e+02,
     . 8.51654e+02, 8.89400e+02, 4.54768e+02, 9.12158e+01/
      data name(25) /'kr'/, (coeff(25,i),i = 1,8)/
     .2.000e-01,2.000e+00,-1.81364e+01,-1.51898e+00,
     . 2.83187e+00, 6.12136e+00,-1.07998e+01,-1.57415e+01/
      data name(26) /'kr'/, (coeff(26,i),i = 1,8)/
     .2.000e+00,2.000e+01,-2.33830e+01, 3.90879e+01,
     .-1.05886e+02, 1.26461e+02,-7.04126e+01, 1.50163e+01/
      data name(27) /'kr'/, (coeff(27,i),i = 1,8)/
     .2.000e+01,1.000e+02,-2.63927e+01, 1.95442e+01,
     .-2.06960e+01, 1.09324e+01,-2.88191e+00, 3.06312e-01/
c *********  molybdenum        **********
      data name(28) /'mo'/, (coeff(28,i),i = 1,8)/
     .6.000e-02,2.000e-01,-1.39105e+02,-6.49334e+02,
     .-1.36584e+03,-1.40646e+03,-7.08621e+02,-1.40057e+02/
      data name(29) /'mo'/, (coeff(29,i),i = 1,8)/
     .2.000e-01,2.000e+00,-1.77259e+01,-1.05822e+00,
     .-3.58317e+00, 1.66009e+00, 8.56537e+00, 4.53291e+00/
      data name(30) /'mo'/, (coeff(30,i),i = 1,8)/
     .2.000e+00,2.000e+01,-1.38510e+01,-3.67845e+01,
     . 1.14059e+02,-1.63563e+02, 1.07626e+02,-2.64249e+01/
      data name(31) /'mo'/, (coeff(31,i),i = 1,8)/
     .2.000e+01,1.000e+02, 3.99268e+01,-1.75709e+02,
     . 2.07493e+02,-1.21459e+02, 3.53180e+01,-4.08383e+00/
c *********  nickel            **********
      data name(32) /'ni'/, (coeff(32,i),i = 1,8)/
     .3.000e-02,2.000e-01,-1.20325e+01, 3.25391e+01,
     . 6.79077e+01, 6.52992e+01, 2.97346e+01, 5.27128e+00/
      data name(33) /'ni'/, (coeff(33,i),i = 1,8)/
     .2.000e-01,2.000e+00,-1.83048e+01,-3.31924e-03,
     .-3.33231e+00,-1.11280e+01, 1.05307e-01, 9.44891e+00/
      data name(34) /'ni'/, (coeff(34,i),i = 1,8)/
     .2.000e+00,2.000e+01,-1.69768e+01,-9.49547e+00,
     . 1.10936e+01, 4.04590e-02,-6.52193e+00, 2.65491e+00/
      data name(35) /'ni'/, (coeff(35,i),i = 1,8)/
     .2.000e+01,1.000e+02,-2.86408e+01, 2.99929e+01,
     .-3.72608e+01, 2.25806e+01,-6.71660e+00, 7.91169e-01/
c *********  oxygen            **********
      data name(36) /'o '/, (coeff(36,i),i = 1,8)/
     .2.000e-03,2.000e-02,-4.52975e+02,-9.70819e+02,
     .-8.54356e+02,-3.70311e+02,-7.89495e+01,-6.59517e+00/
      data name(37) /'o '/, (coeff(37,i),i = 1,8)/
     .2.000e-02,2.000e-01,-5.85047e+01,-1.59971e+02,
     .-2.41395e+02,-1.60654e+02,-4.43559e+01,-3.33047e+00/
      data name(38) /'o '/, (coeff(38,i),i = 1,8)/
     .2.000e-01,2.000e+00,-2.07261e+01,-7.67836e-01,
     . 8.75731e-01,-1.08357e+00, 2.41126e+00, 5.75616e+00/
      data name(39) /'o '/, (coeff(39,i),i = 1,8)/
     .2.000e+00,2.000e+01,-2.06891e+01,-1.09982e+00,
     . 2.00613e+00,-1.58614e+00, 6.69091e-01,-1.11969e-01/
      data name(40) /'o '/, (coeff(40,i),i = 1,8)/
     .2.000e+01,1.000e+02,-2.78060e+01, 2.14606e+01,
     .-2.66591e+01, 1.67083e+01,-5.19194e+00, 6.41029e-01/
c *********  silicon           **********
      data name(41) /'si'/, (coeff(41,i),i = 1,8)/
     .2.000e-03,1.800e-02, 1.55215e+03, 3.32531e+03,
     . 2.75818e+03, 1.12074e+03, 2.23142e+02, 1.74141e+01/
      data name(42) /'si'/, (coeff(42,i),i = 1,8)/
     .1.800e-02,2.000e-01,-3.33775e+01,-4.15621e+01,
     .-3.35750e+01,-3.86101e-01, 9.55951e+00, 2.92688e+00/
      data name(43) /'si'/, (coeff(43,i),i = 1,8)/
     .2.000e-01,2.000e+00,-1.94215e+01,-1.56426e-01,
     .-5.29353e+00,-1.24578e+00, 2.42390e+01, 1.86013e+01/
      data name(44) /'si'/, (coeff(44,i),i = 1,8)/
     .2.000e+00,2.000e+01,-1.92475e+01,-1.98358e+00,
     . 8.42427e-01, 7.54957e-01,-6.74371e-01, 1.46365e-01/
      data name(45) /'si'/, (coeff(45,i),i = 1,8)/
     .2.000e+01,1.000e+02,-1.91107e+01,-2.66860e+00,
     . 2.43841e+00,-1.00956e+00, 2.21647e-01,-2.07805e-02/
c *********  tungsten          **********
      data name(46) /'w '/, (coeff(46,i),i = 1,8)/
     .1.000e-01,2.000e-01, 5.34083e+00, 1.56088e+02,
     . 4.17170e+02, 5.50258e+02, 3.56758e+02, 9.04279e+01/
      data name(47) /'w '/, (coeff(47,i),i = 1,8)/
     .2.000e-01,2.000e+00,-1.72389e+01, 5.42375e-02,
     .-1.22107e+00, 4.41181e-01,-4.48582e+00,-7.83614e+00/
      data name(48) /'w '/, (coeff(48,i),i = 1,8)/
     .2.000e+00,2.000e+01,-1.47488e+01,-1.43954e+01,
     . 2.10585e+01,-4.39475e+00,-1.10601e+01, 5.61699e+00/
      data name(49) /'w '/, (coeff(49,i),i = 1,8)/
     .2.000e+01,1.000e+02,-2.62426e+02, 7.12559e+02,
     .-8.25017e+02, 4.74241e+02,-1.35517e+02, 1.54189e+01/
      data nrows/49/
c
c ----------------------------------------------------------------------
c
      do 6000 j=1,nj
      tkev    = te(j)
      prad(j) = 0.0
      do 5000 imp=1,nimp
c
      do i=1,nrows
        ki = i
        if (name(i) .eq. namei(imp))  go to 2010
      end do
c
      write  (ncrt, 8000) namei(imp)
      write  (nout, 8000) namei(imp)
 8000 format (' FATAL ERROR in RADFIT'                            /
     .        ' impurity ', a8, ' has not been added to the code' /)
      call STOP ('subroutine RADFIT: problem #1', 86)
c
c     label 2010 takes care of lower limit.
c     Implicit assumption here is that table is
c     appropriately ordered(eg first entry for each impurity
c     is the lowest energy one )!!!!! HSJ 
c     Use this lower limit even if tkev is less than c(i,1)
 2010 i = ki
      if (tkev .ge. coeff(i,1))  go to 2050
      irow = i
      y    = coeff(irow,8)
      x    = LOG10 (coeff(irow,1))
      do i=8,4,-1
        y  = y*x + coeff(irow,i-1)
      end do
      rad  = 10.0**y
      rad  = tkev*(rad/coeff(irow,1))
      go to 4999
c
 2050 i0 = i
      do i=i0,nrows
        ki = i
        if (tkev .ge. coeff(i,1) .and. tkev .le. coeff(i,2)) go to 2100
c        if (name(i) .ne. namei(imp))  go to 9000      

         if (name(i) .ne. namei(imp))THEN
         ! back up to last valid range for impurity i, HSJ,09/11/06:
            ki = ki-1
            Write(nout,8020)name(ki), namei(imp), tkev, coeff(ki,1),
     .                                                coeff(ki,2)
            WRITE(ncrt,8020)name(ki), namei(imp), tkev, coeff(ki,1),
     .                                                coeff(ki,2)

            GO TO 2100
         endif
      end do
      go to 9000
c
 2100 irow = ki
      y    = coeff(irow,8)
      x    = LOG10 (tkev)
      do i=8,4,-1
        y  = y*x + coeff(irow,i-1)
      end do
      rad  = 10.0**y
c
 4999 brems_nions(j,nprim+imp) =  rad*eni(j,imp)*1.0e-7 
 5000 prad(j) = prad(j) + rad*eni(j,imp)*1.0e-7
 6000 continue
      return
c
c ----------------------------------------------------------------------
c fatal errors
c ----------------------------------------------------------------------
c
 9000 write  (nout, 8010)  name(i), namei(imp), tkev, coeff(i,1),
     .                                                coeff(i,2)
      write  (ncrt, 8010)  name(i), namei(imp), tkev, coeff(i,1),
     .                                                coeff(i,2)
 8010 format (' FATAL ERROR in RADFIT'                              /
     .        ' you have exceeded temperature range for curve fits' /
     .        ' name(i), namei(imp)', 2x, a, 2x, a                  /
     .        ' tkev, lower limit, upper limit =', 3(2x, 1pe12.4))
      call STOP ('subroutine RADFIT: problem #2', 87)
c
 8020 format ('WARNING ERROR in RADFIT'                              /
     .        ' you have exceeded temperature range for curve fits' /
     .        ' name(i), namei(imp)', 2x, a, 2x, a                  /
     .        ' tkev, lower limit, upper limit =', 3(2x, 1pe12.4),  /
     .        ' code will continue by using largest energy  value ', /
     .        ' in the tables')
      end

      real*8 function radhe (temkev)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c     this function calculates the radiation losses from helium as a
c     function of elec. temp. (keV) using a formula derived by
c     George Hopkins and John M. Rawls (ga-a14141 impurity radiations)
c
c     temkev = temperature of electrons in keV
c     value returned is in units of  (watts-cm**3)
c
      x = temkev
      if (temkev .gt. 100.00)  x = 100.00
      if (temkev .lt.   0.01)  x =   0.01
      if (     x .gt.   0.4 )  go to 10
      x = LOG10 (x)
      y = 10.0**(-35.632+x*(0.473114+x*0.521255))
      go to 20
   10 sqrtem = SQRT (x)
      y      = 2.136e-36*sqrtem+2.08e-37/sqrtem+4.8e-38/(x**1.34)
   20 radhe  = y*1.0e6
      return
c
      end

      subroutine rate (ratef)
c
      USE param
      USE ions
      USE neut
      USE soln
      USE numbrs
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c     RATE computes the neutral reaction rate coefficients <sigma-v>
c     for electron ionization and charge exchange. The units are cm**3/s
c     ratef is the correction factor for elongated plasmas.
c
c      include 'param.i'
c      include 'ions.i'
c      include 'numbrs.i'
c      include 'neut.i'
c      include 'soln.i'
c
      if (nneu .eq. 2)  go to 110
c
c          Only one neutral species - ion number in
c
      iother = 3 - in
      do j=1,njs
        eionr(j)       = ratef*eira(te(j),atw(in))
        eirate(j)      = eionr(j)*ene(j)
        cexr(j,in)     = ratef*cxra(ti(j),atw(in))
        cexr(j,iother) = 0.0
        cx12r(j)       = 0.0
      end do
      return
c
c          Two neutral species - use average atw for mixed charge
c          exchange reactions, and use only hydrogenic ionization rate.
c
  110 atwm        = (atw(1)+atw(2))/2.0
      do j=1,njs
        eionr(j)  = ratef*eira(te(j),1.0)
        eirate(j) = eionr(j)*ene(j)
        cexr(j,1) = ratef*cxra(ti(j),atw(1))
        cexr(j,2) = ratef*cxra(ti(j),atw(2))
        cx12r(j)  = ratef*cxra(ti(j),atwm)
      end do
      return
c
      end




c      subroutine read_iter_ufile_12(filename,f,iwant,tin,
c     .                              )
c---------------------------------------------------------------
c     subroutine reads ufile information for use in Onetwo
c
c-----------------------------------------------------HSJ-------
c       implicit none
c       integer task
cc       include 'param.i'
c       !include 'transp.i'
c       character *(*) filename
       
c       iounit = 95
c       call getioun(iounit,iounit)c
c       open (unit = iounit, file = filename, status = 'OLD', 
c     .                                      iostat = iostat)
c
c      if (iostat .ne. 0) then
c        write  (ncrt, 1000)  filename(1:LENGTH(filename)), nin
c 1000   format (/ ' ---- ERROR:  file "', a, '" on logical unit', i3   /
c     .                     14x, 'cannot be opened, therefore ',
cc     .                          'execution of ONETWO cannot continue;' /
c     .                     14x, 'check that this file exists ',
c     .                          'and that it is readable')
c        call giveupus(iounit)
c        call STOP ('read_iter_ufile: bad or missing file', 1)
c      end if

c      if(iwant .eq. -1)then ! first call read in time and rho arrays
c                            !and dimensions

c      else                  
c         !subsequent calls read in arrays at two time
c         !levels which bracket the time of interest
c         !start by bracketing the time tin
c         task = -1
c         call 2d_lin_interp(tin,rhon,nj,transp_times, transp_rho,
c     .                          array_in,nx,nt,array_out,i1,i2,
c     .                          task)
c         if(i1 le. 0)
c     .     call STOP ('read_iter_ufile: failed time inerpolation', 1)
c         rewind(iounit)
c         !space down to consecutive indecies  i1,i2
         
c      endif
c      !call giveupus(iounit) keep unit open

c      return
c      end

c      subroutine 2d_lin_interp(tin,rhon,nj,tarray_in, rarray_in,
c     .                          array_in,nx,nt,array_out,i1,i2,
c     .                           task)
c----------------------------------------------------------------------
c        given array_in(nx,nt) nx = space,nt= time,get array_out(nx) at
c        time tin using linear interpolation in space and time.
c INPUT:
c tin            time at which vales are to be found
c rhon(nj)       normalized [0,1] rho grid.Results will be obtained for
c                values of rho in rhon.
c nj             size of rhon
c tarray_in      input list of times at which second index of array_in is given
c                tarray_in must span the time over which the transport
c                simulation is done.Out of bounds interpolation is fatal.
c rarray_in      input list of normalized rho values at which first index of
c                array_in is given. Must cover [0,1]
c array_in(i,j)  i = 1 is assumed to corespond to rhon(1)  = magnetic axis
c                i = nx                           rhon(nj) = plasma edge
c                intermediate values, 1 < i < nj are NOT assumed
c                to match the values in rhon(i) since typically nx < nj
c nx,nt          actual space and time sizes of arrays
c     
c task           -1 only fetch i1,i2 (see output below)



c OUTPUT:
c i1,i2           time indecies only .
c                 tarray_in(i1) .le. tin .le. tarray_in(i2)
c                 failure means i1 is set to -1 (stop code)
c                 0 return ?? 
c array_out(nj) spatial values at time tin            
c------------------------------------------------------------HSJ--------
c      implicit none
c      integer nx,nt,nj,j,i1,i2,k,k1,k2
c      real *8 tin,tarray_in(nt),array_out(nx),rarray_in(nx),
c    .           rhon(nj),a,array_in(nx,nt),a1,a2,b,rhow


c      !first get the time interval that brackets tin:
c      do j=1,nt-1
c         dt = (tin-tarray_in(j))*(tin-tarray_in(j+1))
c         if(dt .le. 0.0)then
c            i1 = j
c            i2 = j+1
cc            go to 10
c         endif
c      enddo
c      i1 = -1
c      call STOP('Sub transp_2d_interp time out of bounds',0)


c      !now get the space interval and do interpolation
c 10   if(task .eq. -1)return
c      do j=2,nj-1
c         rhow = rho(j)
c         do k = 1,nx-1
c           drho  = (rhow-rarray_in(k))*(rhow-rarray_in(k+1))
c           if(drho .le. 0.0)then
c             k1 = k
c             k2 = k+1
c           endif
c           !first TIME INTERPOLATION:
c           a = (array_in(k1,i2)-array_in(k1,i1))/(tarray_in(i2)
cc    .                                          -tarray_in(i1))
c           b = array_in(k1,i1)-a*tarray_in(i1)
c           a1 = a*tin+b  ! value of array_in at rho = rarray_in(k1),time =tin
c           a = (array_in(k2,i2)-array_in(k2,i1))/(tarray_in(i2)
c    .                                          -tarray_in(i1))
c           b = array_in(k2,i1)-a*tarray_in(i1)
c           a2 = a*tin+b  ! value of array_in at rho = rarray_in(k2) ,time =tin
c           !next  SPACE INTERPOLATION:
c           a = (a2-a1)/(rarray_in(k2) -rarray_in(k1))
c           b = a1 -a*rarray_in(1)
c           array_out(j) = a*rhow + b  !interpolated value at rhow and tin
c           go to 20
c         enddo
c         call STOP('Sub transp_2d_interp rho  out of bounds',0)
c 20      continue
c      enddo

c      !value at the magnetic axis:
c      a = (array_in(1,i2)-array_in(1,i1))/(tarray_in(i2)-tarray_in(i1))
c      b = array_in(1,i1) - a*tarray_in(i1)
c      array_out(1) = a*tin + b

c      !value at the plasma edge
c      a = (array_in(nx,i2)-array_in(nx,i1))/(tarray_in(i2)-tarray_in(i1))
c      b = array_in(nx,i1) - a*tarray_in(i1)
cc      array_out(nj) = a*tin + b

c      return
c      end
      




      subroutine recall (u, usave, nk, nj, kk)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c  this subroutine sets the current solution u(k,j) equal to usave(k,j).
c
      dimension u(kk,*), usave(kk,*)
c
      do 10 j=1,nj
      do 10 k=1,nk
   10 u(k,j) = usave(k,j)
      return
c
      end

      real*8 function recrat (te)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c          RECRAT computes the local rate of radiative recombination
c          of thermal ions and electrons into neutrals (#/cm**3-s).
c          The <sigma-v> reaction rate coefficient has been fit
c          by a rational polynomial to the formula in Gordeev et. al.
c          [JETP Lett.,V.25,No.4,p.204,1977] to better than 0.1%.
c          Multiply this rate by Z of the ions in the calling routine.
c          The limits on Te are as specified in Gordeev et.al.
c          Te is in keV.
c
      t = te
      if (t .lt. 1.0e-3) t = 1.0e-3
      if (t .gt. 25.0) t = 25.
      t = LOG10 (t)
      recrat = (-15.42043+t*(-6.215405+t*(-1.559212-t*0.8086629e-1)))/
     .(1.0 + t*(0.3166686+t*(0.6894623e-1-t*0.1031934e-2)))
      recrat = 10.0**recrat
      return
c
      end

      subroutine redate (u, en, te, ti, rbp, nk, nj, kj, kk,
     .                   iangrot, angrot)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c  this subroutine redates the solution variables u
c  by setting them equal to the dependent variables (en,te,ti,rbp,angrot)
c
      dimension u(kk,*), en(kj,*), te(*), ti(*), rbp(*), angrot(*)
c
      do 20 j=1,nj
      do 10 k=1,nk-3-iangrot
   10 u(k,j) = en(j,k)
      u(nk-2-iangrot,j) = te(j)
      u(nk-1-iangrot,j) = ti(j)
      u(nk-iangrot,j)   = rbp(j)
   20 if (iangrot .eq. 1)  u(nk,j) = angrot(j)
      return
c
      end

      subroutine reduce (j, a, b, c, g, em,u, itran, nk, nj, kk,
     .                   diffeq_methd)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
      integer diffeq_methd
c
c  this subroutine reduces the order of the matrices a, b, and c
c     and the vector g from nk to nkt.
c
      dimension a(kk,*), b(kk,*), c(kk,*), g(*),em(kk,*),
     .           u(kk,*), itran(*)
c
c  calculate new elements of g corresponding to variables with
c     values that are not specified
c
      do 30 k=1,nk
      if (itran(k) .le. 0)  go to 30
      do 20 l=1,nk
      if (itran(l) .eq. 1)  go to 20
      if (j .eq. 1)  go to 10
      g(k) = g(k) - a(k,l)*u(l,j-1)
   10 g(k) = g(k) - b(k,l)*u(l,j)
      if (j .eq. nj)  go to 20
      g(k) = g(k) - c(k,l)*u(l,j+1)
   20 continue
   30 continue
c
c  reduce a, b, c, and g by eliminating all rows and columns
c     corresponding to specified variables
c
      kt = 0
      do 50 k=1,nk
      if (itran(k) .le. 0)  go to 50
      kt = kt + 1
      do 40 l=1,nk
      a(kt,l) = a(k,l)
      b(kt,l) = b(k,l)
      if(diffeq_methd .eq. 1)em(kt,l)= em(k,l)
   40 c(kt,l) = c(k,l)
      do 45 l=1,nk
      a(l,kt) = a(l,k)
      b(l,kt) = b(l,k)
      if(diffeq_methd .eq. 1)em(l,kt)= em(l,k)
   45 c(l,kt) = c(l,k)
      g(kt) = g(k)
   50 continue
      return
c
      end


       SUBROUTINE reinit_bp0
c ---------------------------------------------------------------------
c initialize rbp = F*G*H*r*Bp0 to the value that it would have
c if the driven current and profiles were frozen at the current values.
c Note that this ignores the profiles obtained from the eqdsk
c --------------------------------------------------HSJ----------------
       USE geom,   only    :  hcap,fcap
       USE numbrs, only    :  nj
       USE CONSTNTS, only  :  twopi
       USE sourc, only     :  etap,curohm,curdri,curboot   
       USE soln, only      :  rbp,curden,etor
       USE mesh, only      :  r
       USE machin, only    :  rmajor
       IMPLICIT NONE
       real *8,allocatable,dimension (:) ::  f,feta
       real *8 E_const,xint,vloop ,toti,totni

          allocate(f(nj),feta(nj))


       
c     first get the value of the steady state constant electric field,
c     E_const, in  volts/cm :

          feta(1:nj) = etap(1:nj)
          if(feta(nj) .eq.0.0)feta(nj)=feta(nj-1)
          if(feta(1) .eq.0.0)feta(1)=feta(2)
          call eff_resist(feta,r,nj,twopi,xint)
c         total non inductive current, jni:
          f(1:nj) = curboot(1:nj) + curdri(1:nj)
          call trapv (r,f , hcap, nj, totni) !assumes < jni RO/R> instead
                                             !of < Jni dot B / Bt0>
          totni = totni*twopi                !which  may be problematic
                                             !but since this is just zn
                                             !an inital gues for curden 
                                             !and rbp it may be good enough 
          toti = 5.*rbp(nj)


          E_const = (toti - totni )/xint !volts/cm     
c          vloop = twopi * rmajor * E_const
          curohm(1:nj) = E_const/(8.98755179e11*feta(1:nj)*hcap(1:nj))
          curden(1:nj) = curohm(1:nj)+ curboot(1:nj) + curdri(1:nj)
          !if curden is negative anywhere then it is an indication that
          !we will not be able to find a steady state solution.
          WHERE(curden <=0)curden = 1.0
          call trap3(r, curden, nj, hcap,twopi,f)
          print *,'totcurden,amps  in reinit_bp0 =',f(nj)
          rbp(1:nj-1) = f(1:nj-1)*0.2*fcap(1:nj-1)   !rbp in  gauss cm
          etor(1:nj) = E_const

          deallocate(f,feta)
       RETURN
       END SUBROUTINE reinit_bp0

       SUBROUTINE eff_resist(eta,r,nj,twopi,xint)
c-----------------------------------------------------------------------
       IMPLICIT none
       REAL *8 eta(nj),r(nj),twopi,xint,f(nj)
       INTEGER nj

          f(:) = r(1:nj)/eta(1:nj)
          call trapl (r, f, nj, xint)
          xint = twopi*xint/8.98755179e11         !converts eta to ohm-cm

       RETURN
       END SUBROUTINE eff_resist








      subroutine reord (i)
      USE param
      USE soln
      USE numbrs
      USE flags
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c  subroutine reorders the components of the solution vector u.
c  i = 1 ; called by reduce - reorder u so that transported
c        quantities appear first
c  i = 2 ; called by expand - reorder u to its original form
c ----------------------------------------------------------------------
c
c      include 'param.i'
c      include 'flags.i'
c      include 'numbrs.i'
c      include 'soln.i'
c
      dimension  ipos(kk), iposi(kk), ifill(kk), save1(kj), save2(kj)
c
c ----------------------------------------------------------------------
c  rearrange u so that transported quantities occur first.
c ----------------------------------------------------------------------
c
      loc = 0
c
      do 100 k=1,nk
        if (itran(k) .eq. 0)  go to 100
        loc     = loc + 1
        ipos(k) = loc
  100 continue
c
      do 110 k=1,nk
        if (itran(k) .eq. 1)  go to 110
        loc     = loc + 1
        ipos(k) = loc
  110 continue
c
      if (i .eq. 1)  go to 116
c
c ----------------------------------------------------------------------
c  determine original positions for components of u
c ----------------------------------------------------------------------
c
      do 112 k=1,nk
      ki        = ipos(k)
  112 iposi(ki) = k
      do 114 k=1,nk
  114 ipos(k) = iposi(k)
c
  116 do 120 k=1,nk
  120 ifill(k) = 0
c
c ----------------------------------------------------------------------
c  start switching components of u around
c ----------------------------------------------------------------------
c
      do 190 k=1,nk
      newk = ipos(k)
c
c ----------------------------------------------------------------------
c  does k need to be   moved ?
c  has it already been moved ?
c ----------------------------------------------------------------------
c
      if (      newk  .eq. k)  ifill(newk) = 1
      if (ifill(newk) .eq. 1)  go to 190
c
c ----------------------------------------------------------------------
c save contents of k
c ----------------------------------------------------------------------
c
      do 140 j=1,nj
  140 save2(j) = u(k,j)
c
c ----------------------------------------------------------------------
c save contents of newk
c ----------------------------------------------------------------------
c
  150 do 160 j=1,nj
  160 save1(j) = u(newk,j)
c
c ----------------------------------------------------------------------
c  store k in newk
c ----------------------------------------------------------------------
c
      do 170 j=1,nj
  170 u(newk,j)   = save2(j)
      ifill(newk) = 1
c
c ----------------------------------------------------------------------
c  determine next newk
c ----------------------------------------------------------------------
c
      newk = ipos(newk)
      if (ifill(newk) .eq. 1)  go to 190
c
c ----------------------------------------------------------------------
c  store and switch
c ----------------------------------------------------------------------
c
      do 180 j=1,nj
      save2(j)    = u(newk,j)
  180 u(newk,j)   = save1(j)
      ifill(newk) = 1
      newk        = ipos(newk)
      if (ifill(newk) .eq. 1)  go to 190
      go to 150
  190 continue
      return
c
      end



c      subroutine  reset_edge_values(te,ti,angrot,fix_edge_te,
c     .                     fix_edge_ti,fix_edge_rot,u,kk,nk,nj,iangrot)
c----------------------------------------------------------HSJ-----
c      implicit none
c      integer*4 fix_edge_te,fix_edge_ti,fix_edge_rot,
c     .                                   kk,nj,nk,iangrot,j
c      real *8 te(*),ti(*),u(kk,*),angrot(*)
c      do j = 1 ,nj
c         if(j .ge. fix_edge_te)u(nk-2-iangrot,j)= te(j)
c         if(j .ge. fix_edge_ti)u(nk-1-iangrot,j)= ti(j)
c         if(iangrot .eq. 1 .and. j .ge. fix_edge_rot)u(nk,j)= angrot(j)
c      enddo
c      return
c      end

     


c
      
      subroutine rhomsh (time)
c
      USE param
      USE soln
      USE contour
      USE limiter
      USE mhdpar
      USE tdem
      USE numbrs
      USE extra
      USE mesh
      USE adaptive
      USE machin
      USE geom
      USE flags
      USE tordlrot
      USE soln2d
      USE solcon, only : steady_state
      USE mhdcom
      USE ifs
      USE neo2d
      USE tmpcom
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c      include 'param.i'
c      include 'flags.i'    ! itran
c      include 'extra.i'    ! bpol
c      include 'geom.i'
c      include 'machin.i'   ! btor
c      include 'mesh.i'
c      include 'numbrs.i'   ! nk,
c      include 'tmpcom.i'
c      include 'soln.i'     ! rbp
c      include 'soln2d.i'
c      include 'neo2d.i'
c      include 'mhdpar.i'   ! to get mhdcom.i
c      include 'mhdcom.i'   ! mhdmethd
c      include 'tordlrot.i' ! iangrot
c      include 'contour.i'  ! for tdem.i
c      include 'limiter.i'  ! for tdem.i
c      include 'tdem.i'
c      include 'adaptive.i' ! include_adaptive
c      include 'ifs.i'
c
      data pi_squared /9.869604404/
      data ffmult,fgmult,fhmult,frhomult /1.0,1.0,1.0,1.0/
c
      if (codeid .ne. 'onedee')  go to 30
c
c ----------------------------------------------------------------------
c    this routine was really designed for the twodee
c    code to describe the time-dependence of the
c    geometric factors and the rho mesh.
c    here it allows fixed shape calculations by suitably defining gcap.
c ----------------------------------------------------------------------
c
      if (nbctim .le. 1 .or. elong(2) .eq. 0.0)  go to 19
c
c     the following gives variation of plasma elongation, at constant rminor
c
      call elongt (time, kappa, dkapdt)
      dkapdt = steady_state*dkapdt
      rhoa       = SQRT (kappa)*rminor
      rhoasv     = r(nj)
      x          = rhoa/rhoasv
      do 10 j=1,nj
      fcap(j)    = 1.0
      gcap(j)    = (1.0 + kappa*kappa)/(2.0*kappa)
      hcap(j)    = 1.0
      r2cap(j)   = 1.0
      r2capi(j)  = rmajor**2
      rcap(j)    = rmajor
      rcapi(j) = 1./rmajor
      dfdt(j)    = 0.0
      dgdt(j)    = 0.5*(1.0-1.0/(kappa*kappa))*dkapdt
      dhdt(j)    = 0.0
      dr2dt(j)   = 0.0
      dr2idt(j)  = 0.0
      drcapdt(j) = 0.0
      drcapidt(j) = 0.0
      r(j)       = x*r(j)
      drr(j)     = drr(j)/x
      rrm(j)     = rrm(j)/x
      bsqncap(j)=r2capi(j)/(rmajor+r(j))**2
      bsq_avg_cap(j)=1./bsqncap(j)
      grho1_ifs(j)=1.
      grho2_ifs(j)=1.0
      elong_r(j)=1.
      if (j .eq. nj)  go to 10
      ra(j)      = ra(j)*x
      dr(j)      = dr(j)*x
      rrp(j)     = rrp(j)/x
   10 continue
      ascrip     = 1.0/rhoa
      bscrip     = 0.0
      drhoadt_geom     = 0.5*dkapdt*rminor / SQRT (kappa)
      adot       = (-drhoadt_geom/rhoa**2) * adotmult
      bdot       = 0.0
      go to 100
c
c          (we only need to do this once for onedee cases)
c          (with constant elongation)
c
   19 rhoa       = r(nj)
      do 20 j=1,nj
      dfdt(j)    = 0.0
      dgdt(j)    = 0.0
      dhdt(j)    = 0.0
      dr2dt(j)   = 0.0
      dr2idt(j)  = 0.0
      drcapdt(j) = 0.0
      drcapidt(j) = 0.0
      fcap(j)    = 1.0
      gcap(j)    = (1.0 + kappa*kappa)/(2.0*kappa)
      hcap(j)    = 1.0
      r2cap(j)   = 1.0
      r2capi(j)  = rmajor**2
      rcap(j)    = rmajor
      rcapi(j)   = 1./rmajor
      bsqncap(j)=r2capi(j)/(rmajor+r(j))**2 !hsj 1/21/2000 added
      bsq_avg_cap(j)=1./bsqncap(j)          !hsj 1/21/2000  "
      grho1_ifs(j)=1.                       !hsj 1/21/2000  "
      grho1_ifs(j)=1.                       !hsj 1/21/2000  "
      eps(j)     = r(j)/SQRT(kappa)/rmajor  !HSJ 4/23/01
   20 continue
      drhoadt_geom     = 0.0
      ascrip     = 1.0
      go to 100
c
c this subroutine calculates the correct mesh for time = time
c the maximum rho for this time is rhoa.
c other quantities which are advanced explicitly in time
c are also calculated here:
c              r2cap,r2capi
c fcap, hcap and rho may be advanced implicitly (done in subroutine FHRCALC)
c
   30 if (mhdmethd .eq. 'tdem') then
           igetcur = 0      ! get both current and mhd parms
           call update_tdem_data(time,igetcur)
      else

              do j=1,nj
                 gcap(j) = gcap0(j)+dgdt(j)*(time-timcap)
                 if (.not. implicit_fh) then
                    fcap(j) = fcap0(j)+dfdt(j)*(time-timcap)
                    hcap(j) = hcap0(j)+dhdt(j)*(time-timcap)
                    if (xi_include) then
                       eps0(j)  = eps(j)
                       xhm20(j) = xhm2(j)
                       xi110(j) = xi11(j)
                       xi330(j) = xi33(j)
                       xips0(j) = xips(j)
                    end if
                 end if
                 r2cap(j)  = r2cap0(j)+dr2dt(j)*(time-timcap)
                 r2capi(j) = r2capi0(j)+dr2idt(j)*(time-timcap)
                 rcap(j)   = rcap0(j)+drcapdt(j)*(time-timcap)
                 rcapi(j)   = rcap0i(j)+drcapidt(j)*(time-timcap)
                 if (xi_include) then
                    eps(j)  = eps0(j)+depsdt(j)*(time-timcap)
                    xhm2(j) = xhm20(j)+dxhm2dt(j)*(time-timcap)
                    xi11(j) = xi110(j)+dxi11dt(j)*(time-timcap)
                    xi33(j) = xi330(j)+dxi33dt(j)*(time-timcap)
                    xips(j) = xips0(j)+dxipsdt(j)*(time-timcap)
                 end if
              end do
              rhoa   = rhoa0+drhoadt_geom*(time-timcap)
      end if
c
c ----------------------------------------------------------------------
c  because the original formulation transformed only time derivatives
c  at constant rho to time derivatives at constant zeta,without making
c  the corresponding change in the rho derivatives,we have to
c  define the new rho mesh as x times the old rho mesh( to
c  insure that normalized rho mantains the same value at each
c  grid point). this is required because the time derivatives are at
c  constant zeta,where zeta is an implicit variable, defined as rho/rhomax
c
c -----------------------------------------------------------HSJ--------
c
          rhoasv = r(nj)
          x      = rhoa/rhoasv
c
          if (include_adaptive .eq. 0) then ! no adaptive grid calcs
              do 60 j=1,nj
                 r(j)   = x*r(j)
                 drr(j) = drr(j)/x
                 rrm(j) = rrm(j)/x
                 if (j .eq. nj)  go to 60
                 ra(j)  = ra(j)*x
                 dr(j)  = dr(j)*x
                 rrp(j) = rrp(j)/x
   60         continue
           else
               call set_rho(time)
           end if
c
          ascrip = 1.0/rhoa
          bscrip = 0.0
          adot   = (-drhoadt_geom/rhoa**2) * adotmult
          bdot   = 0.0
c
  100 do 110 j=1,njs
  110 roa( j) = r(j)/r(nj)
      roa(nj) = 1.0
      call radius_calc    
c
       if (mhdmethd .eq. 'tdem' .and.
     .     itran(nk-iangrot) .eq. 0) then
c
c        analysis mode for rbp: get rbp from data read from netcdf file
c        done here because r(j) is now known at the desired time.
c
         do j=1,nj
            bpol(j)=bp0_tdem(j)*10000.        ! Gauss
            rbp(j)=fcap(j)*gcap(j)*hcap(j)*r(j)*bpol(j)
c
c           may want to renormalize rbp here
c
            q(j) = SIGN (q_tdem(j), btor)
            if (j.eq. 1) then      ! drbpdtdem_tdem is defined as (1/FGHr)*
               drbpdt_tdem(j)=0.0  ! d(FGHrBp0)/dt
            else
               drbpdt_tdem(j)=bpol(j)*(dfdt(j)/fcap(j)*ffmult
     .                                 +dgdt(j)/gcap(j)*fgmult
     .                                 +dhdt(j)/hcap(j)*fhmult
     .                                 -roa(j)*drhoadt_geom*frhomult)
     .                                 +dbdt_tdem(j)*10000.  ! gauss/sec
            end if
         end do
         call curcalc(rbp, fcap, hcap, gcap, r, curden,
     .                                           curpar_soln, nj)
       end if
c
      call psirho(0)     ! get psir corresponding to r mesh
c
c --- sfarea and cxarea, will be recalculated in fhcalc if appropriate
c
      volfac = 4.0 * pi_squared * rmajor
      sfarea = volfac * hcap(nj) * r(nj)
      return
c
      end

      subroutine savesol (u, usave, nk, nj, kk)
      USE solcon,only : steady_state
      USE soln,only   : enesav,ene
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c --- this subroutine saves the current solution u(k,j) in usave(k,j)
c
      dimension  u(kk,*), usave(kk,*)
c
      do   j=1,nj
         if(steady_state .lt. 1.e-5)enesav(j) = ene(j) !HSJ 8/13/2004
        do k=1,nk
          usave(k,j) = u(k,j)
        end do
      end do
      return
c
      end

      subroutine scale_length (pro, var, nj, out)
c
      implicit none
c
c ----------------------------------------------------------------------
c return out(j=1,..nj) = (-1/pro)(dpro/dvar) (i.e., the gradient scale
c                        length)
c
c input  pro(j),var(j) where pro is the profile defined
c        as a function of the independent variable variable
c -------------------------------------------------------- HSJ ---------
c
      real*8  pro(*), var(*), out(*)
      integer j, nj
c
      call difydx (var, pro, out, nj)
c
      do j=1,nj
        out(j) = -out(j) / pro(j)
      end do
      return
c
      end

      real*8 function seval (n, u, x, y, b, c, d)
c
c this subroutine evaluates the cubic spline function
c
c    seval = y(i) + b(i)*(u-x(i)) + c(i)*(u-x(i))**2 + d(i)*(u-x(i))**3
c
c    where  x(i) < u < x(i+1), using Horner's rule
c
c  if  u .lt. x(1)  then  i = 1  is used.
c  if  u .ge. x(n)  then  i = n  is used.
c
c  ---------------------------HSJ Modification 1/30/98------------------
c  The above two rules are not acceptible in the context in which
c  this routine is used in the ONETWO code. (actually they are not
c  acceptable in any context in my opinion.
c  Interpolating outside the cubic spline fit is totally bogus.)
c  I have changed the rules to read
c   if u .gt. x(n) then use Y(n) as the most reasonable
c   approximation if you dont actually want to stop the code.
c   similarly if u .lt. x(1) then use y(1) as the returned value.
c  using the first or last point may not be much better in some circumstances
c  but at least we will get physically relevant values !!!!!!!!!!!HSJ
c  Note that the switch extend_seval must be set to 1 to get this
c  behavior. Otherwise an error stop is taken when trying to interpoalte
c  outside the tables.
c -----------------------------------------------------------------------
c
c  input..
c
c    n     = the number of data points
c    u     = the abscissa at which the spline is to be evaluated
c    x,y   = the arrays of data abscissas and ordinates
c    b,c,d = arrays of spline coefficients computed by spline
c
c  if  u  is not in the same interval as the previous call, then a
c  binary search is performed to determine the proper interval.
c
c
      USE param
      USE io
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      integer n
      real*8  u, x(n), y(n), b(n), c(n), d(n)
c      include 'param.i'
c      include 'io.i'
      integer i, j, k
      real*8  dx
      data    i/1/
c
c
c  added  1/29/98  HSJ
 1    format('Subroutine Seval,a spline fit evaluater,',/,
     .       ' has detected that out of bounds interpolation',/,
     .       ' is occuring. The value at which the cubic spline',/,
     .       ' is to be evaluated is ',1pe18.8,/,
     .       ' the interpolating table extends from ',1pe18.8,
     .       ' to ',1pe18.8,/,
     .       ' You can set extend_seval=1 in the first namelist',/,
     .       ' in inone if you want to ignore this problem')
      if (u .lt. x(1)) then
         if (extend_seval .eq. 1) then
           seval = y(1)
           return
         else
           write (ncrt, 1)  u, x(1), x(n)
           write (nout, 1)  u, x(1), x(n)
           call STOP ('subroutine SEVAL: bad interpolation #1', 267)
         end if
      else if (u .gt. x(n)) then
         if (extend_seval .eq. 1) then
           seval = y(n)
           return
         else
           write (ncrt, 1)  u, x(1), x(n)
           write (nout, 1)  u, x(1), x(n)
           call STOP ('subroutine SEVAL: bad interpolation #2', 268)
         end if
      end if
c
c
      if (i .ge. n     )  i = 1
      if (u .lt. x(i  ))  go to 10
      if (u .le. x(i+1))  go to 30
c
c  binary search
c
   10 i = 1
      j = n+1
   20 k = (i+j)/2
      if (u .lt. x(k))  j = k
      if (u .ge. x(k))  i = k
      if (j .gt. i+1 )  go to 20
c
c  evaluate spline
c
   30 dx = u - x(i)
      seval = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
      return
c
      end


      subroutine solve
c---------------------------------------------------------------------------
c

c
c  this subroutine solves for the vector uu that satisfies the
c  linear matrix equation
c       aa*uu = gg
c  the matrix aa and the vectors uu and gg are of order nj*nkt.
c  moreover, the matrix aa is block tridiagonal with nj blocks
c  of order nkt.  the jth block-matrix equation can be written as
c       a(j)*u(j-1) + b(j)*u(j) + c(j)*u(j+1) = g(j) .
c
c  NOTE:
c  kk  = leading dimension of u
c  nk  = number of dependent variables in the problem (as opposed to
c        the total number of dependent variables the code allows, which
c        is set to kk)
c        The value of nk depends on how many primary ions
c        and impurities there are in the input (counting both simulation and
c        analysis mode variables)
c        nk = #prim ions  +#impurity ions   + 1 for te + 1 for ti
c                           +1 for rbp  + 1 for rotation
c        Inclusion of rotation is optional,equations for te,ti,rbp
c        and at least one primary ion are always present so min value of
c        nk is 4. Max value of nk is kk ( or kk-1 if rotation is excluded).
c
c  nkt = number of dependent variables for which a transport solution
c        is to be found (i.e., the number of variables for which
c        itran(i),i = 1,2..nk, is equal to 1).
c  the matrices a, b, c,and em and the vector g are calculated by the
c        subroutine ABCG and are of order nk.  if nkt < nk (i.e., if
c        the values of some dependent variables are specified), then
c        the arrays are reduced to order nkt before the matrix equation
c        is solved.  subsequently the solution vector u is expanded
c        back to order nk.
c
c  because aa is block tridiagonal, the solution uu can be obtained
c     efficiently by a direct method.  specifically, gaussian
c     elimination (i.e., forward elimination followed by backward
c     substitution) is used.  the kth element of the solution
c     vector for block j is stored in u(k,j).
c-------------------------------------------------------------------------
      USE param
      USE io
      USE soln
      USE numbrs
      USE mesh
      USE solcon,only : time
      USE flags
      USE tordlrot
      USE bd_condtn,only : bctime,ub,fluxb,
     .    ub_save,ub_rho_edge,bctime_zone
      USE nonlin
      implicit  integer (i-n), real*8 (a-h, o-z)
c      include 'param.i'
c      include 'bcon.i'
      include 'imsl.i'
c      include 'flags.i'
c      include 'mesh.i'  !te,ti,rot_index
c      include 'numbrs.i'
c      include 'soln.i'
c      include 'tordlrot.i'
c      include 'io.i'
c      include 'nonlin.i'

c
      parameter (kk1 = kk+1)
c

      dimension  e(kk,kk,kj)
      dimension  a(kk,kk),b(kk,kk),c(kk,kk),g(kk),bb(kk,kk),
     .           cc(kk,kk1),f(kk,kj),x(kk),wkarea(kk),em(kk,kk)
      character  ntcc_input_file*16
      data       iacc /6/
c
****  call dumper ! debugging of solution for Alpha port (6 April 1999)
c
      if (nkt .eq. 0)  return
      imslmd = 'solve'
c
      if (testing_NTCC .eq. 1) then
        ntcc_input_file = "ntcc_input_file"
c
c       output file  to be created,fresh each time:
c
        call destroy (ntcc_input_file)
        call getioun(ntcc,ntcc)
        open (unit = ntcc, file = ntcc_input_file, status = 'NEW')
      end if
c
c  calculate scale factors
c
      if(diffeq_methd .eq. 0)
     .      call maxu (u, itran, nk, nj, kk, x, iangrot)
      nkt1 = nkt + 1
      nequations  = 0             !counts number of equations if 
                                  !residual  or ml_resid is called,
                                  !saved in nonlin.i


      do 100 j=1,nj
c
c  calculate the matrices a, b, c,and em, and the vector g
c  the matrix em is required explicitely for the method of
c  lines solution. In that case (specified by diffeq_methd =1),
c  em is not included in matrix b and vector g.
      call abcg (j, a, b, c, g, em)



c     redefine a,b,c,g for those profiles which have the boundary set
c     at a value less than nj:
      if(j .ge. te_index  .and. itran(nk-2-iangrot) .eq. 1)then
         kij = kion+1    !index into matrices for te 
         kij = nion+1    !index into matrices for te 
         do ijk =1,nk   ! over row kij of a,b,c,g
            a(kij,ijk)= 0.0
            b(kij,ijk)= 0.0
            c(kij,ijk)= 0.0
            em(kij,kij) =0.0
         enddo
         b(kij,kij)= 1.
         g(kij) = bctime_zone(j,kion+1)
         if (testing_NTCC .eq. 1) then
            write(ntcc,'("j,nion+1,g(nion+1) =",
     .       i3,2x,i3,2x,1pe15.8)')j,nion+1,g(kij)
         endif
      endif

      if(j .ge. ti_index .and. itran(nk-1-iangrot) .eq. 1)then
         kij = kion+2    !index into matrices for ti 
         kij = nion+2    !index into matrices for ti 
         do ijk =1,nk   ! over row kij of a,b,c,g
            a(kij,ijk)= 0.0
            b(kij,ijk)= 0.0
            c(kij,ijk)= 0.0
            em(kij,kij) =0.0
         enddo
         b(kij,kij)= 1.
         g(kij) = bctime_zone(j,kion+2)
         if (testing_NTCC .eq. 1) then
            write(ntcc,'("j,nion+2,g(nion+2) =",
     .       i3,2x,i3,2x,1pe15.8)')j,nion+2,g(kij)
         endif
      endif

      if(j .ge. rot_index  .and.iangrot .eq. 1
     .                .and.             itran(nk) .eq. 1)then
         kij = kion+4    !index into matrices for toroidal rotation
         kij = nion+4    !index into matrices for toroidal rotation
         do ijk =1,nk   ! over row kij of a,b,c,g
            a(kij,ijk)= 0.0
            b(kij,ijk)= 0.0
            c(kij,ijk)= 0.0
         enddo
         b(kij,kij)= 1.
         g(kij) = bctime_zone(j,kion+4)
      endif

cjmp.den start
      if(j .ge. ni_index  .and. itran(1) .eq. 1)then
         kij = 1    !index into matrices for en(1)
         do ijk =1,nk   ! over row kij of a,b,c,g
            a(kij,ijk)= 0.0
            b(kij,ijk)= 0.0
            c(kij,ijk)= 0.0
            em(kij,kij) =0.0
         enddo
         b(kij,kij)= 1.
         g(kij) = u(j,1)
         if (testing_NTCC .eq. 1) then
            write(ntcc,'("j,nion+1,g(nion+1) =",
     .       i3,2x,i3,2x,1pe15.8)')j,nion+1,g(kij)
         endif
      endif
cjmp.den end

      if (testing_NTCC .eq. 1) then
         write(ntcc,'(5(2x,i5))')j,nkt,nk,kk,nj
         if (j .eq. 1) then
           write(ntcc,'(5(2x,i5))')(itran(jj),jj=1,nk)
           write(ntcc,'(5(2x,1pe16.8))')((u(k,jj),k=1,nk),jj=1,nj)
         end if
         write(ntcc,'(5(2x,1pe16.8))')time
         write(ntcc,'(5(2x,1pe16.8))')((a(jj,jjk),jjk=1,nk),jj=1,nk)
         write(ntcc,'(5(2x,1pe16.8))')((b(jj,jjk),jjk=1,nk),jj=1,nk)
         write(ntcc,'(5(2x,1pe16.8))')((c(jj,jjk),jjk=1,nk),jj=1,nk)
        itestntcc=2
        write(ntcc,'( 2x,i5)')itestntcc
         write(ntcc,'(5(2x,1pe16.8))')(g(jjk),jjk=1,nk)
         write(ntcc,'(5(2x,1pe16.8))')(x(jjk),jjk=1,nk)
      end if
c
c  reduce the arrays to order nkt if nkt < nk
c
      if (nkt .eq. nk)  go to 5
      call reduce (j, a, b, c, g,em, u, itran, nk, nj, kk,diffeq_methd)
      if (testing_NTCC .eq. 1 ) then
         write(ntcc,'(5(2x,1pe16.8))')((a(jj,jjk),jjk=1,nk),jj=1,nk)
         write(ntcc,'(5(2x,1pe16.8))')((b(jj,jjk),jjk=1,nk),jj=1,nk)
         write(ntcc,'(5(2x,1pe16.8))')((c(jj,jjk),jjk=1,nk),jj=1,nk)
        itestntcc=1
        write(ntcc,'( 2x,i5)')itestntcc
         write(ntcc,'(5(2x,1pe16.8))')(g(jjk),jjk=1,nk)
           write(ntcc,'(5(2x,1pe16.8))')((u(k,jj),k=1,nk),jj=1,nj)
      end if

 5    continue
      if(diffeq_methd .eq. 1)THEN
         !METHOD OF LINES SOLUTION
         call ml_resid(j,a,b,c,g)
         go to 100
      else if(diffeq_methd .eq. 2)then
         !NEWTON TYPE SOLVER
c         print *,'j,nkt  = ',j,nkt
c         print *,'a,b,c =',a(nkt,nkt),b(nkt,nkt),c(nkt,nkt)
         call residual(j,a,b,c,g)
         go to 100
      endif
c
c  perform forward elimination
c  PREDICTOR-CORRECTOR SOLUTION:
      do 20 k=1,nkt
      do 10 l=1,nkt
      bb(k,l) = b(k,l)
   10 cc(k,l) = c(k,l)
   20 cc(k,nkt1) = g(k)
      if (j .eq. 1)  go to 50
      do 40 k=1,nkt
      do 40 l=1,nkt
      do 30 m=1,nkt
   30 bb(k,l) = bb(k,l) - a(k,m)*e(m,l,j-1)
   40 cc(k,nkt1) = cc(k,nkt1) - a(k,l)*f(l,j-1)
   50 do 60 k=1,nkt
      do 60 l=1,nkt
   60 bb(k,l) = x(l)*bb(k,l)                ! x is set in maxu, used only if diffeq_methd = 0
          jm1=j-1
          if (jm1 .eq. 0)  jm1=1
      if (testing_NTCC .eq. 1) then
       write(ntcc,'(5(2x,1pe16.8))')((e(k,jj,jm1),k=1,nkt),jj=1,nkt)
        write(ntcc,'(5(2x,1pe16.8))')((cc(k,jj),k=1,nkt),jj=1,nkt1)
        itestntcc=3
        write(ntcc,'( 2x,i5)')itestntcc
        write(ntcc,'(5(2x,1pe16.8))')((bb(k,jj),k=1,nkt),jj=1,nkt)
        itestntcc=4
        write(ntcc,'( 2x,i5)')itestntcc
      end if

      call leqt1fl(bb, nkt1, nkt, kk, cc, iacc, wkarea, ier)
 
      if (testing_NTCC .eq. 1) then
        write(ntcc,'(5(2x,1pe16.8))')((cc(k,jj),k=1,nkt),jj=1,nkt1)
        itestntcc=5
        write(ntcc,'( 2x,i5)')itestntcc
        write(ntcc,'(5(2x,1pe16.8))')((bb(k,jj),k=1,nkt),jj=1,nkt)
      end if
      do 70 k=1,nkt
      do 70 l=1,nkt1
   70 cc(k,l) = x(k)*cc(k,l)
      do 90 k=1,nkt
      do 80 l=1,nkt
   80 e(k,l,j) = cc(k,l)
   90 f(k,j) = cc(k,nkt1)
  100 continue   !ends loop over grid points
 
      if(diffeq_methd .gt. 0 )return
c
c ----------------------------------------------------------------------
c  reorder solution vector
c ----------------------------------------------------------------------
c

      if (testing_NTCC .eq. 1) then
        itestntcc=6
        write(ntcc,'( 2x,i5)')itestntcc
               write(ntcc,'(5(2x,1pe16.8))')((u(k,j),k=1,nk),j=1,nj)
      end if
      if (nk .eq. nkt)  go to 105
      call reord(1)
  105 continue
c
c  perform backward substitution
c
      do 120 j=nj,1,-1
      do 120 k=1,nkt
      u(k,j) = f(k,j)
      if (j .eq. nj)  go to 120
      do 110 l=1,nkt
  110 u(k,j) = u(k,j) - e(k,l,j)*u(l,j+1)
  120 continue
c
c  expand the solution vector to order nk if nkt < nk 
c
      if (testing_NTCC .eq. 1) then 
        itestntcc=7
        write(ntcc,'( 2x,i5)')itestntcc
               write(ntcc,'(5(2x,1pe16.8))')((u(k,j),k=1,nk),j=1,nj)
      end if
      if (nkt .ne. nk)  then
          call expand(u, usave, itran, nk, nkt, nj, kk)
c         call exp1(itran)
      end if
      if (testing_NTCC .eq. 1) then
          write(ntcc,'(5(2x,1pe16.8))')((u(k,j),k=1,nk),j=1,nj)
          write(ntcc,'( "teindex,tiindex,rotindex =",3(1x,i3))')
     .            te_index,ti_index,rot_index
          do j = te_index-5,nj
             write(ntcc,'("j,u(nion+1,j),usave(nion+1,j) =",
     .        i3,2x,2(1pe16.8))')j,u(nion+1,j),usave(nion+1,j)
          enddo

          do j = ti_index-5,nj
             write(ntcc,'("j,u(nion+2,j),usave(nion+2,j) =",
     .        i3,2x,2(1pe16.8))')j,u(nion+2,j),usave(nion+2,j)
          enddo

          call giveupus(ntcc)
          close(unit=ntcc)
          call STOP ('subroutine SOLVE: below line-label 120', 290)
      end if


      return
c
      end


       subroutine sort_unique(a,km)
c --------------------------------------------------------------------
c sort array a in asceding order,eliminate duplicate entries
c -------------------------------------------------------------HSJ----
       USE replace_imsl,                         ONLY : my_vsrta
       IMPLICIT NONE
       REAL *8, intent(INOUT),dimension(:) :: a
       INTEGER, INTENT(OUT) :: km
       INTEGER i,j,k,kmm
    


       k=SIZE(a)
       !sort switch on times,eliminate duplicate times:
       IF(k .GT. 1)THEN
          CALL my_vsrta (a,k)   ! sort in increasing order
          km=k
          DO WHILE(1 .gt. 0 ) !lf95 doesnt like DO WHILE(1) 
             kmm=0
             DO j=2,km
                IF(a(j) .EQ. a(j-1))THEN
                   DO i=j+1,k
                      a(i-1) = a(i)
                   ENDDO
                   km = km-1
                   kmm=1
                   EXIT
                ENDIF
             ENDDO
             IF(kmm .EQ. 0) EXIT
          ENDDO
       ELSE
         km=1
       ENDIF
 
       RETURN
       END




      subroutine specify
c
c ----------------------------------------------------------------------
c calculate the specified time-dependent profiles for
c those variables which are not being transported
c (i.e. analysis mode variables, itxx = 0).
c calculate ene and zeff at t = time only if inenez .ne. 0
c
c NOTE that both u and the local working names rpb,te,..etc
c      are loaded.
c modified 9/25/96 to account for tdem analysis mode
c of Farady's law. HSJ
c 4/4/97 added interpolation of profiles for adaptive grid HSJ
c        kprim,kimp are parameters for maximum # of primary and
c        impurity ions possible
c        kprim+kimp+1 points to ene
c                  +2           te
c                  +3           ti
c                  +4           zeff
c                  +5           curden
c                  +6           angrot
c ---------------------------------------------------------------
c
      USE param
      USE ions
      USE solcon
      USE soln
      USE contour
      USE limiter
      USE mhdpar
      USE tdem
      USE extra
      USE numbrs
      USE mesh
      USE adaptive
      USE sourc
      USE machin
      USE tfact
      USE geom
      USE flags
      USE tordlrot
      USE constnts
      USE soln2d
      USE bd_condtn
      USE mhdcom
      implicit  integer (i-n), real*8 (a-h, o-z)
      INTEGER use_eqc
c

      include 'storage.i'

c
      nj_rho_edge = nj
      if (rho_edge .lt. 1.0)  nj_rho_edge = rho_edge * nj
c
      if (mhdmethd .eq. 'tdem' .and. itran(nk-iangrot) .eq. 0) then
         igetcur = 1         ! get the current related profiles only
         call update_tdem_data (time, igetcur)
         do j=1,nj
            bpol(j) = bp0_tdem(j)*10000.         ! Gauss
            rbp(j)  = fcap(j)*gcap(j)*hcap(j)*r(j)*bpol(j)
            q(j)    = SIGN (q_tdem(j), btor)
         end do
         call curcalc(rbp, fcap, hcap, gcap, r, curden, curpar_soln, nj)
         ucenter(2,nk-iangrot) = 0.0  ! in case inone input is not..
c                                     ..consistent with tdem mode
c
c        renormalize curden to total current:
c
         do  j=1,nj
           ydum(j) = 2.0 * pi * hcap(j)*r(j)*curden(j)
         end do
         call trap2 (r, ydum, rbp, nj) ! note: rbp(j) is now integrated
c
c                                              current in amps
c
         if(.not. u_vloop_bc)
     .    call interp1 (time,bctime,nbctim,bc(1,nk-iangrot),
     .                 ub(nk-iangrot))
c
c        ub(nk-iangrot) contains input current at this time /5.
c        (=fcap(nj)*gcap(nj)*hcap(nj)*r(nj)*bp0(nj) in gauss -cm )
c
         fac = ub(nk-iangrot)/rbp(nj)
         do  j=1,nj
           curden(j) = 5.0 * fac * curden(j)
           rbp(j)    = fac*fcap(j)*rbp(j)     ! get normalized rbp
           u(nk-iangrot,j) = rbp(j)
           if (nj_rho_edge .eq. j)  ub_rho_edge(nk-iangrot) = rbp(j)
         end do
         if(.not. u_vloop_bc)
     .    ub(nk-iangrot) = u(nk-iangrot,nj) ! differ only due to roundoff
      end if
c
      if (include_adaptive .gt. 0 .and. nbctim .eq. 1) go to 4000
c
      if (nbctim .eq. 1)  return
c
c     time-dependent profiles with standard (i.e., parabolic) parameters
c     angular rotation does not have parabolic input so use nk-iangrot
c
      do 2200 k=1,nk-iangrot   !note NO rho_edge option here
        if (    itran(k) .eq. 1  )  go to 2200
        if (ucenter(2,k) .eq. 0.0)  go to 2200
        call interp1(time,bctime,nbctim,ucenter(1,k),y0)
        call interp1(time,bctime,nbctim,uedge(1,k),ynj)
        call interp1(time,bctime,nbctim,ualp(1,k),yalp)
        call interp1(time,bctime,nbctim,ugam(1,k),ygam)
        call makpro(r,ydum,nj,y0,ynj,yalp,ygam)
c
c       Faraday's law: if this is a time dependent eqdsk mode run
c
        if (k .ne. nk-iangrot)  go to 2080
        do 2020 j=1,nj
          curden(j) = ydum(j)
 2020     ydum(j) = 2.0 * pi * hcap(j)*r(j)*curden(j)
        call trap2(r,ydum,rbp,nj)
        call interp1(time,bctime,nbctim,bc(1,k),ub(k))
        fac = ub(k)/rbp(nj)
        do 2040 j=1,nj
          curden(j) = 5.0 * fac*curden(j)
 2040   ydum(j) = fac*fcap(j)*rbp(j)
 2080   do j=1,nj
          u(k,j) = ydum(j)
        end do
        ub(k) = u(k,nj)
 2200 continue
c
c     time-dependent spline profiles (te and ti and rotation
c                                      if splninpt = 'old')
c     otherwise te,ti,curden,rotation
c     note that time dependent spline input of  primary and impurity
c     densities is not implemented (enp(j,ion) and eni(j,ion) are
c     not implemented as a function of time)
c
      if (itran(nk-2-iangrot) .eq. 1 .and. (te_index .lt. nj))
     .  go to 2215
c
c     te_index will be set to a number > nj if option not in effect
c
      if (itran(nk-2-iangrot) .eq. 1 .and. (nj_rho_edge .eq. nj))
     .                   go to 2220
 2215 if (njte  .eq. 0)  go to 2220
c
      if (splninpt .eq. 'old') then
        if (tein(1,2) .eq. 0.0)  go to 2220
        call tsplin   (tein, rnormin, njte, zdum)
      else
        call tsplinew (tein, rtein, zdum, knotste, kprim+kimp+2,
     .                 xdum, ydum)
      end if
c
      do   j=1,nj
        if (itran(nk-2-iangrot) .eq. 0)
     .      u(nk-2-iangrot,j) = zdum(j)
        if (j .eq. nj_rho_edge)
     .      ub_rho_edge(nk-2-iangrot) = zdum(j)
        if (j .ge. te_index)
     .      u(nk-2-iangrot,j) = zdum(j)  ! see comment above
      end do
      ub(nk-2-iangrot) = zdum(nj)
c
 2220 if (itran(nk-1-iangrot) .eq. 1 .and. (ti_index .lt. nj))
     .  go to 2225
c
c     i_index will be set to a number > nj if option not in effect
c
      if (itran(nk-1-iangrot) .eq. 1 .and. (nj_rho_edge .eq. nj))
     .                   go to 2240
 2225 if (njti  .eq. 0)  go to 2240
      if (splninpt .eq. 'old') then
        if (tiin(1,2) .eq. 0.0)  go to 2240
        call tsplin (tiin, rnormin, njti, zdum)
      else
        call tsplinew (tiin, rtiin, zdum, knotsti, kprim+kimp+3,
     .                 xdum, ydum)
      end if
c
      do  j=1,nj
        if (itran(nk-1-iangrot) .eq. 0)
     .      u(nk-1-iangrot,j) = zdum(j)
        if (j .eq. nj_rho_edge)
     .      ub_rho_edge(nk-1-iangrot) = zdum(j)
        if (j .ge. ti_index)
     .      u(nk-1-iangrot,j) = zdum(j) ! see comment above
      end do
      ub(nk-1-iangrot) = zdum(nj)
c
c --- current profile
c
 2240 if (splninpt .eq. 'new'
     .             .and. itran(nk-iangrot) .eq. 0
     .             .and. mhdmethd .ne. 'tdem') then
              !if knots are not set we take it to mean that
              !the current density from the eqdsk is to be used. HSJ,7/30/07
              use_eqc = 0
              DO j=1,nbctim
                 if(knotscur(j) .gt. 0)use_eqc =1
              ENDDO
              IF(use_eqc == 1)THEN
                call tsplinew (curdenin,rcurdein,zdum,knotscur,
     .                          kprim+kimp+5, xdum,ydum)
                do 2241 j=1,nj
                  curden(j) = zdum(j)
 2241             ydum  (j) = 2.0 * pi * hcap(j)*r(j)*curden(j)
                call trap2(r,ydum,rbp,nj) ! note redefinition of rbp here
                call interp1(time,bctime,nbctim,bc(1,nk-iangrot),
     .                                          ub(nk-iangrot))
                fac = ub(nk-iangrot)/rbp(nj)
                do 2042 j=1,nj
                  curden(j) = 5.0 * fac*curden(j)
                  rbp(j)=fac*fcap(j)*rbp(j)
                  if (j .eq. nj_rho_edge)
     .                ub_rho_edge(nk-iangrot) = rbp(j)
 2042             u(nk-iangrot,j) = rbp(j)
              ENDIF
      end if

c
c --- special input settings , set_te_to_ti  and set_ti_to_te = 0,1
c --- used if one of te,ti is evolved and the other is to be
c --- set equal to it.
c
      DO j=1,nj
          IF(set_te_to_ti ==1 .AND. itran(nk-2-iangrot) ==0)THEN
              u(nk-2-iangrot,j) =  u(nk-1-iangrot,j)
          ELSEIF(set_ti_to_te == 1 .AND. itran(nk-1-iangrot) == 0)THEN
              u(nk-1-iangrot,j)  = u(nk-2-iangrot,j)
          ENDIF
      ENDDO

c
c --- transp beam driven current profile
c
      if(external_beam_cur .gt. 0)then
          call tsplinew (curbeam_external,rcurb_external,curb_external,
     .                   knotscurb_external,kprim+kimp+8,xdum,ydum)
          call trapv(r,curb_external,hcap,nj,curbtot_external)
      endif

      if (iangrot  .eq. 0                            )  go to 2250
      if (itangrot .eq. 1 .and. (rot_index .lt.  nj))  go to 2245
      if (itangrot .eq. 1 .and. (nj_rho_edge .eq. nj))  go to 2250
c
c --- angular rotation speed profile
c
 2245 call tsplinew (angrotin,rangrot,zdum,knotsang,kprim+kimp+6,
     .               xdum,ydum)
      do j=1,nj
        if ( itran(nk) .eq. 0)u(nk,j) = zdum(j )
        if (j .eq. nj_rho_edge)ub_rho_edge(nk)=zdum(j)
        if (j .ge. rot_index)
     .      u(nk,j) = zdum(j) ! see comment above
      end do
      ub(nk) = u(nk,nj) 
c
c     now copy vector u to en,te,ti,rbp,angrot:
c
 2250 call update(u,en,te,ti,rbp,nk,nj,kj,kk,iangrot,angrot)
c
c     if inenez # 0, the ion densities are calculated from
c     ene and zeff.  if the ene and zeff profiles are
c     time-dependent they are updated for the present time.
c
      if (inenez  .eq. 0  )  go to 2500
      if (ttweak  .ne. 0.0)  go to 2500
      if (enec(2) .eq. 0.0)  go to 2300
c
c --- parabolic time-dependent ne
c
      call interp1 (time,bctime,nbctim,enec,y0)
      call interp1 (time,bctime,nbctim,eneb,ynj)
      call interp1 (time,bctime,nbctim,alpene,yalp)
      call interp1 (time,bctime,nbctim,gamene,ygam)
      call makpro (r,ene,nj,y0,ynj,yalp,ygam)
      go to 2320
 2300 if (njene .eq. 0)  go to 2320
c
c --- spline time-dependent ne
c
      if (splninpt .eq. 'old') then
        if (enein(1,2) .eq. 0.0)  go to 2320
        call tsplin(enein,rnormin,njene,ene)
      else
c
        call tsplinew (enein, renein, ene, knotsene, kprim+kimp+1,
     .                 xdum, ydum)
      end if
c
c ----------------------------------------------------------------------
c  zeff profile:
c ----------------------------------------------------------------------
c
 2320 if (njzef .eq. 0) then                    ! parabolic profile
          if (zeffc(2) .eq. 0.0)  go to 2330    ! constant in time
c
c --- parabolic time-dependent zeff
c
          call interp1(time,bctime,nbctim,zeffc,y0)
          call interp1(time,bctime,nbctim,zeffb,ynj)
          call interp1(time,bctime,nbctim,alpzef,yalp)
          call interp1(time,bctime,nbctim,gamzef,ygam)
          call makpro(r,zeff,nj,y0,ynj,yalp,ygam)
 2330     continue
      else               ! spline profile for zeff
c
c --- spline time-dependent zeff
c
          if (splninpt    .eq. 'old') then
            if (zeffin(1,2) .eq.  0.0)  go to 2500
            call tsplin (zeffin, rnormin, njzef, zeff)
          else
            call tsplinew (zeffin, rzeffin, zeff, knotszef,
     .                     kprim+kimp+4, xdum, ydum)
          end if
      end if
c
 2500 if (include_adaptive .eq. 0) return
c
c      move those profiles that are not transported to
c      the adaptive grid values. They are currently defined on
c      the mesh r_mesh (r_mesh is the radial mesh at the begining
c      of the time step) and need to be interpolated to the
c      mesh r ( the r grid is the adaptive mesh used to move the solution
c      from the begining of the time step to the end of the time step.
c      The r mesh does not change during the corrector steps)
c      Note that parabolic profiles will be spline fitted and
c      interpolated.
c
 4000 if (include_adaptive .gt. 0) then
             do i=1,nk      ! do all profiles
                  if (itran(i) .eq. 0) then
                    call copyas(u(i,1),kk,ydum,1,nj)
                    call intrp_adp (r_mesh, ydum, nj, r, zdum, nj)
                    call copyas(zdum,1,u(i,1),kk,nj)
                  end if
             end do
             if (inenez .ne. 0) then
c                 ne
                 call copya(ene,ydum,nj)
                 call intrp_adp (r_mesh, ydum, nj, r, ene, nj)
c                 zeff
                 call copya(zeff,ydum,nj)
                 call intrp_adp (r_mesh, ydum, nj, r, zeff,nj)
             end if
             call update(u,en,te,ti,rbp,nk,nj,kj,kk,iangrot,angrot)
      end if
c
      end

      subroutine spline_12 (n, x, y, b, c, d)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      integer n
      real*8  x(n), y(n), b(n), c(n), d(n)
c
c  the coefficients b(i), c(i), and d(i), i = 1,2,...,n are computed
c  for a cubic interpolating spline
c
c    s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
c
c    for  x(i) .le. x .le. x(i+1)
c
c  input..
c
c    n = the number of data points or knots (n .ge. 2)
c    x = the abscissas of the knots in strictly increasing order
c    y = the ordinates of the knots
c
c  output..
c
c    b, c, d  = arrays of spline coefficients as defined above.
c
c  using  p  to denote differentiation,
c
c    y(i) = s(x(i))
c    b(i) = sp(x(i))
c    c(i) = spp(x(i))/2
c    d(i) = sppp(x(i))/6  (derivative from the right)
c
c  the accompanying function subprogram  seval  can be used
c  to evaluate the spline.
c
      integer nm1, ib, i
      real*8  t
c
      nm1 = n - 1
      if (n .lt. 2)  return
      if (n .lt. 3)  go to 50
c
c  set up tridiagonal system
c
c  b = diagonal, d = offdiagonal, c = right hand side.
c
      d(1) = x(2) - x(1)
      c(2) = (y(2) - y(1))/d(1)
      do i=2, nm1
         d(i)   = x(i+1) - x(i)
         b(i)   = 2.0 * (d(i-1) + d(i))
         c(i+1) = (y(i+1) - y(i))/d(i)
         c(i)   = c(i+1) - c(i)
      end do
c
c  end conditions.
c  third derivatives at x(1) and x(n) obtained from divided differences
c
      b(1) = -d(1)
      b(n) = -d(n-1)
      c(1) = 0.0
      c(n) = 0.0
      if (n .eq. 3)  go to 15
      c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
      c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
      c(1) = c(1)*d(1)**2/(x(4)-x(1))
      c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
c
c  forward elimination
c
   15 do i=2, n
        t    = d(i-1)/b(i-1)
        b(i) = b(i) - t*d(i-1)
        c(i) = c(i) - t*c(i-1)
      end do
c
c  back substitution
c
      c(n) = c(n)/b(n)
      do ib=1, nm1
         i = n-ib
         c(i) = (c(i) - d(i)*c(i+1))/b(i)
      end do
c
c  c(i) is now the sigma(i) of the text
c
c  compute polynomial coefficients
c
      b(n) = (y(n) - y(nm1))/d(nm1) + d(nm1)*(c(nm1) + 2.0*c(n))
      do i=1, nm1
         b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.0*c(i))
         d(i) = (c(i+1) - c(i))/d(i)
         c(i) = 3.0*c(i)
      end do
      c(n) = 3.0*c(n)
      d(n) = d(n-1)
      return
c
   50 b(1) = (y(2)-y(1))/(x(2)-x(1))
      c(1) = 0.0
      d(1) = 0.0
      b(2) = b(1)
      c(2) = 0.0
      d(2) = 0.0
      return
c
      end

      subroutine thermonuclear_rate (ddfusn,ddpfus,dtnfus,ttnfus,
     .                         he3dpfus,ddntot,ddptot,dtntot,
     .                         ttntot,he3dptot,iddfus,
     .                         idhe3fus,ihe,ihei,he3_frac,ti,en,
     .                         kj,nj,id,
     .                         idt,it,fd,volfac,hcap,r)
c
      USE verbose
      USE numbrs,            ONLY : nprim
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ------------------------------------------------- 11/20/95 --- HSJ ---
c     BASIC THERMAL FUSION REACTION RATE CALCUALTIONS:
c
c     DDFUSN    = D(D,   N)HE3       THERMAL-THERMAL (I.E., THERMONUCLEAR)
c     DDPFUS    = D(D,   P)T           "       "             "
c     DTNFUS    = D(T,   N)HE4         "       "             "
c     TTNFUS    = T(T,  2N)HE4         "       "             "
c     HE3DPFUS  = HE3(D, P)HE4         "       "             "
c
c
c             RATE COEFFICIENT AND CROSS SECTIONS are from
c             Bosch & Hale, Nuc. Fus., vol32, no.4 (1992) 611
c             (except for t(t,2n)he4 ,see subroutine ttnrate )
c
c     INPUT
c
c     iddfus   =  0    no deuterium or tritium in system
c              =  1    thermal species is d, no tritium in system
c              =  2    thermal species is dt mixture,
c                      use fraction of dt density
c                      given by fd  to specify d density
c              =  3    thermal species are d and t separate (beam is not dt)
c              =  4    thermal species is  t (no d in system)
c              =  5    thermal species are d and t (beam is dt)
c
c     ti(j)            j=1,2..nj  ion temperature,kev
c     en(j,i)          j=1,2..nj,i=1,2..nprim, primary ion densities,#/cm**3
c
c     fd               fraction of deuterium in dt mixture
c                      (note that fd applies only for iddfus =2 )
c
c     volfac           volume factor =4.0*pi**2*rmajor
c     r(j)             j=1,2,..nj rho grid,cm
c     id               index for species deuterim
c     idt                                dt mixture
c     it                                 tritium
c     hcap
c
c     OUTPUT:
c      thermal d,d reactants (if 0 < iddfus < 3)
c             ddfusn(j)   local reaction (i.e., neutron production) rate,
c                        #/cm**3/sec
c             ddntot       volume integrated rate, #/sec
c      thermal d,d reactants (if 0 < iddfus < 3)
c             ddpfus(j)   local reaction (i.e., proton production) rate,
c                        #/cm**3/sec
c             ddptot       volume integrated rate, #/sec
c      thermal d,t reactants: (if 1 < iddfus < 4)
c             dtnfus(j)   local reaction (i.e., neutron production) rate,
c                        #/cm**3/sec
c             dtntot       volume integrated rate, #/sec
c      thermal t,t reactants: ( if 1 < iddfus <= 4 )
c             ttnfus(j)   local reaction (i.e., neutron production) rate,
c                        #/cm**3/sec
c             ttntot       volume integrated rate, #/sec
c      thermal he,d reactants: 
c             he3dpfus(j)   local reaction (i.e., proton production) rate,
c                        #/cm**3/sec
c             he3dptot       volume integrated rate, #/sec
c
c ----------------------------------------------------------------------
c
c
      dimension  ti(*),en(kj,*),ddfusn(*),dtnfus(*),ttnfus(*),
     .           he3dpfus(*),ddpfus(*),hcap(*),r(*)
c
c     zero volume integrated rates:
c
      ddntot   = 0.0
      dtntot   = 0.0
      ttntot   = 0.0
      he3dptot = 0.0
      ddptot   = 0.0
c
c     zero profiles:
c
      call zeroa (ddfusn, kj)
      call zeroa (dtnfus, kj)
      call zeroa (ttnfus, kj)
      call zeroa (ddpfus, kj)
      call zeroa (he3dpfus, kj)
c
      if (iddfus .eq. 0 )  return  ! thermal species is not d or t or dt
c                                    mixture so no fusion reactions take place
c
      iddfushold = iddfus
      if (iddfus .eq. 5)  iddfus=3
      do j=1,nj           ! begin loop over rho(j)
         if (iddfus .le. 3) then
c
c               we have thermal deuterium in system
c
c               d(d,n)he3 reaction:
c
                if (iddfus .eq. 1 .or. iddfus .eq. 3)
     .                               endfus = en(j,id)
                if (iddfus .eq. 2)  endfus =
     .                     fd*en(j,idt)  ! thermal species is dt mixture
c
                sigmav = ddnrate(ti(j)) / 2.0      ! like-like reaction,
c                                                    divide by 2
c
                ddfusn(j) = (endfus**2) * sigmav   ! sigmav in cm**3/sec
c
c               d(d,p)t reaction:
c
                sigmav=ddprate(ti(j))/2.
                ddpfus(j) = (endfus**2) * sigmav   ! sigmav in cm**3/sec
c
         end if    ! d(d,n)he3,d(d,p)t branch
c
         if (iddfus .gt. 1) then
c                  we have thermal tritium in the system
c                  Thermal tritium ions in the system  are present either
c                  as a sparate species (iddfus =3 or 4 )
c                  or as part of the dt mixture species (iddfus=2)
c                  d density is endfus,t density is endfust
c
                   if (iddfus .eq. 3 .or. iddfus .eq. 4)
     .                   entfus = en(j,it)        ! thermal species is t
                   if (iddfus .eq. 2)
     .                   entfus = (1.0-fd)*en(j,idt)      ! t density in
c                                                           dt mixture
c
             if (iddfus .ne. 4) then   ! d is present so get dt reaction
c
c              bulk plasma d-t fusion:
c
               sigmavt   = dtrate(ti(j))
               dtnfus(j) = endfus *
     .                     entfus * sigmavt    ! sigmavt is in cm**3/sec
c
             end if
c
c            bulk plasma t-t fusion, we get 2 neutrons per fusion.
c            ttnfus is 1/2 the neutron rate!!
c
             sigmavt   = ttnrate(ti(j))   ! sigmavt is in cm**3/sec
             ttnfus(j) = 0.5*entfus*entfus * sigmavt        ! like-like
c                                                             collisions
         end if

            !he3 fusion,default is idhe3fus =0,no helium in system

            IF(idhe3fus .NE. 0)call he3den(enhe3fus,j)
            IF(idhe3fus      == 1)THEN      ! prim ion d and  he3
               endfus    = en(j,id)
               !enhe3fus  = en(j,ihe)*he3_frac
            ELSEIF(idhe3fus  == 2)THEN      ! prim ion d , imp ion he3
               endfus    = en(j,id)
               !enhe3fus  = en(j,nprim+ihei)*he3_frac
            ELSEIF(idhe3fus  == 3)THEN      ! prim ion d, prim ion he3 and imp ion  he3        
               endfus     = en(j,id)
               ! assume  same he3_frac,T imp = T prim ion:
               !enhe3fus   = (en(j,ihe)+en(j,nprim+ihei))*he3_frac 

            ELSEIF(idhe3fus  == 4)THEN      ! dt mixture and prim ion he3
               endfus = fd*en(j,idt) 
               !enhe3fus  = en(j,ihe)*he3_frac
            ELSEIF(idhe3fus   == 5)THEN      ! impurity he3 and dt mixture
               endfus = fd*en(j,idt) 
               !enhe3fus  = en(j,nprim+ihei)*he3_frac
            ELSEIF(idhe3fus   == 6)THEN      ! prim he3 and imp he3 and dt mixture
               endfus = fd*en(j,idt) 
               !enhe3fus  = (en(j,ihe) + en(j,nprim+ihei))*he3_frac 
            ELSEIF(idhe3fus .NE. 0)THEN
               WRITE(ncrt,FMT='("Sub thermonuclear rate encounterd",
     .           " bad idhe3fus value = ",i6)')idhe3fus
               call STOP('thermonuclear rate,idhe3fus not set',1)
            ENDIF
            IF(idhe3fus .NE. 0)THEN ! he3(d,p)he4 ,Maxwellian at Ti
               sigmav= he3dprate(ti(j))
               he3dpfus(j) = endfus*enhe3fus*sigmav
            ENDIF
c
        end do            ! end loop over rho(j)
c
        if (iddfus .le. 3) then
           call trapv (r, ddpfus , hcap, nj, ddptot  )
           ddptot   = ddptot   * volfac     ! volfac = 4 * pisq * rmajor
           call trapv (r, ddfusn , hcap, nj, ddntot  )
           ddntot   = ddntot   * volfac
        end if
        if (iddfus .gt. 1) then
           call trapv (r, ttnfus , hcap, nj, ttntot  )
           ttntot   = ttntot   * volfac
        end if
        if (iddfus .gt. 1 .and. iddfus .ne. 4) then
           call trapv (r, dtnfus , hcap, nj, dtntot  )
           dtntot   = dtntot   * volfac
        end if

        IF(idhe3fus .ne. 0)THEN
           call trapv (r, he3dpfus , hcap, nj,he3dptot )
           he3dptot   = he3dptot   * volfac
        ENDIF
        if (fusionvb .gt. 0) then
          write (*, '(" thermal-thermal d(d,n)he3  rate per second",
     .                  1pe12.3 /
     .            " thermal-thermal d(t,n)he4  rate per second",
     .                  1pe12.3 /
     .            " thermal-thermal t(t,2n)he4 rate per second",
     .                  1pe12.3 /
     .            " thermal-thermal d(d,p)t    rate per second",
     .                  1pe12.3 /
     .            " thermal-thermal he3(d,p)he4 rate per second",
     .            1pe12.3)')  ddntot, dtntot, ttntot, ddptot,he3dptot
        end if
        iddfus = iddfushold
      return
c
      end

      subroutine timederiv (dt, n, znew, zold, dzdt, nj)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c --- subroutine calculates the derivative of cap parameters handled implicitly
c
c  input
c    dt              time,in sec. defined by the fact that znew
c                    is at time t and zold is at time t-dt
c
c    znew(j)         j = 1,2..nj, function (of rho) at time t
c    zold(j)                                              t-dt
c
c    n               for the initial time step znew(j) is not
c                    known so set dzdt(j) = 0.0 if n =0
c
c  output:
c    dzdt(j)         j = 1,2..nj =0 if n=0
c                               = (znew(j)-zold(j))/dt if n .ne. 0
c
c ------------------------------------------------------------------ HSJ
c
      dimension  znew(*), zold(*), dzdt(*)
c
      if (n .eq. 0) then
          do j=1,nj
              dzdt(j) = 0.0
          end do
      else
          do j=1,nj
             dzdt(j) = (znew(j)-zold(j))/dt
          end do
      end if
c
      return
c
      end

      subroutine trap1 (r, y, fact, nj, xint) !jmp.den
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c this subroutine integrates fact*r*y(r) with respect to r,
c from zero to rminor; the trapezoidal rule is used
c
      dimension r(*), y(*)
c
      xint = 0.0
      do 10 j=2,nj
   10 xint = xint + 0.5*(r(j-1)*y(j-1)+r(j)*y(j))*(r(j)-r(j-1))
      xint = xint * fact
      return
c
      end

      subroutine trap2 (x, y, yint, npts)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c calculates integral y*dx using trapezoidal rule
c ----------------------------------------------------------------------
c
      dimension  x(*), y(*), yint(*)
c
      yint(1) = 0.0
      do 10 i=2,npts
   10 yint(i) = (y(i)+y(i-1))*(x(i)-x(i-1))*0.5+yint(i-1)
      return
c
      end

      subroutine trap3 (r, y, nj, fact, const, yint)
c---------------------------------------------------------
c      integral of {y*r*fact}*const is returned as a function
c      of r
c--------------------------------------------------------HSJ
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      dimension  r(*), y(*), fact(*), yint(*)
c
      yint(1) = 0.0
      do j=2,nj
        xint = 0.5*(r(j-1)*fact(j-1)*y(j-1)+r(j)*y(j)
     .         *fact(j))*(r(j)-r(j-1))
        yint(j) = yint(j-1)+xint*const
      end do
      return
c
      end

      real*8 function trapf (j, r, y, fact, const)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      dimension  fact(*), r(*), y(*)
c
      xint  = 0.5*(r(j-1)*fact(j-1)*y(j-1)+r(j)*y(j)
     .                   *fact(j))*(r(j)-r(j-1))
      xint  = xint*const
      trapf = xint
      return
c
      end

      subroutine trapl (r, y, nj, xint)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c this subroutine integrates y(r) with respect to r from zero to rminor
c the trapezoidal rule is used
c
      dimension  r(*), y(*)
c
      xint = 0.0
      do 10 j=2,nj
 10   xint = xint + 0.5 * (y(j)+y(j-1)) * (r(j)-r(j-1))
      return
c
      end

      subroutine trapv (r, y, fact, nj, xint)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c     this subroutine integrates fact(r)*r*y(r) with respect to r from
c     zero to rminor.  the trapezoidal rule is used.
c
      dimension  r(*), fact(*), y(*)
c
      xint = 0.0
      do 10 j=2,nj
   10 xint = xint+0.5*(r(j-1)*fact(j-1)*y(j-1)+r(j)*y(j)
     .      *fact(j))*(r(j)-r(j-1))
      return
c
      end

      subroutine trapvp (r, y, fact, js, nj, xint) !jmp.den
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c     this subroutine integrates fact(r)*r*y(r) with respect to r from
c     zero to rminor.  the trapezoidal rule is used.
c
      dimension  r(*), fact(*), y(*)
c
      xint = 0.0
      do 10 j=js+1,nj
   10 xint = xint+0.5*(r(j-1)*fact(j-1)*y(j-1)+r(j)*y(j)
     .      *fact(j))*(r(j)-r(j-1))
      return
c
      end

      subroutine trplot (flag)

c
c ----------------------------------------------------------------------
c this routine records transport data in the file trpltfil.
c this data is later plotted by a standalone plot program named TRPLOT.
c
c  Revisions:   03/26/92 S.Thompson - Bootstrap current added to
c                                     output file
c               04/15/92 S.Thompson - Current drive, RF current, and
c                                     ohmic current added to output file
c               04/17/92 S.Thompson - XJB terms added to output file
c               04/20/92 S.Thompson - RBP added to output file
c               06/01/92 S.Thompson - GRAD_TE, CRIT_GRAD, XCHIIRL, and
c                                     XCHIERL (Rebut-Lallia) added to output
c               08/29/94 G.M.Staebler - ANGROT added to output
c               03/28/95 HSJ          DDFUSN,DDBEAM,DDKNCK (thermal, beam,and
c                                     knock on d(d,n)he3 neutron rate profiles
c                                     and DDBMT,DDNTHM,DDKNCT vol integrated
c                                     rates added.
c               01/26/96 HSJ          dt neutron rates added
c ----------------------------------------------------------------------
c
      USE param
      USE bd_condtn,                   ONLY : dpedt,dpidt
      USE fusion
      USE ions
      USE solcon
      USE glf23
      USE io
      USE neut
      USE nub
      USE nub2
      USE soln
      USE contour
      USE limiter
      USE mhdpar
      USE rf
      USE tdem
      USE extra
      USE numbrs
      USE mesh
      USE sourc
      USE machin
      USE nub4
      USE geom
      USE transp, only:beam_data,use_nubeam
      USE flags
      USE tordlrot
      USE tcoef
      USE bd_condtn,only : bctime,ub,fluxb,
     .    ub_save,ub_rho_edge,bctime_zone,totcur
      USE aid_newton, ONLY : ssqr_hist
      USE mcgo
      USE flx
      USE neo2d,  epslocal  => eps ! because  beam_plot_dat.i has eps in it
      USE tmpcom
      USE rad_loss,            ONLY : brems_nions,brems_prim,brems_imp
      USE island
      USE zeffcom
      USE Plasma_properties,        ONLY : neut_beam
      USE P_Nfreya_12_interface,    ONLY : use_P_Nfreya

      USE pelcom
      implicit  integer (i-n), real*8 (a-h, o-z)

      include 'beam_plot_dat.i'                        !nbeam_pts,timplot(),waveform()
      include 'co2.i'



      include 'ohml.i'
      include 'rebut.i'



      include 'storage.i'
      include 'sxrcom.i'



      include 'sauter.i' ! sl31,sl32,sl34
c
      dimension    holdte(10)

      REAL*8 dpidt_tot(kj)
      integer      flag, absflg, tdba,get_time_dep_beam
      real*8       kevperg
      data         kevperg /6.242e+08  /,
     .             xmassp  /1.673e-24  /,
     .             pi      /3.141592654/,
     .             tdba   /0/
      data five_halfs_te, five_halfs_ti /2.5,2.5/
      real tarray(2)
      real function etime !these are pgf90 intrinsics
      real function dtime

      if(no_te_convection .eq. -1)five_halfs_te = 1.5
      if(no_ti_convection .eq. -1)five_halfs_ti = 1.5
      if(no_te_convection .eq. 1)five_halfs_te =0.0
      if(no_ti_convection .eq. 1)five_halfs_ti =0.0
c
      zero   = 0.0
      const  = 4.0 * pi * pi*rmajor*1.6e-16
      call trapv (r,qtfusi,hcap,nj,ptfusi)
      ptfusi = ptfusi*const
      call trapv (r,qtfuse,hcap,nj,ptfuse)
      ptfuse = ptfuse*const
      call trapv (r,qbfusi,hcap,nj,pbfusi)
      pbfusi = pbfusi*const
      call trapv (r,qbfuse,hcap,nj,pbfuse)
      pbfuse = pbfuse*const
      call trapv (r,qohm,hcap,nj,pohm)
      pohm   = pohm*const
      call trapv (r,qrfe,hcap,nj,prfe)
      prfe   = prfe*const
      call trapv (r,qrfi,hcap,nj,prfi)
      prfi   = prfi*const
      pbsum  = 0.0

      call trapv(r,pwf_tot_source,hcap,nj,pwf_tot_source_intg)
      pwf_tot_source_intg = pwf_tot_source_intg *const   !watts

      call trapv(r,enbeam,hcap,nj,enbeam_intg)
      enbeam_intg = enbeam_intg *const/1.6e-16  !total # fast ions in system at this time
c
      do   jb=1,nbeams
        do ic=1,ke
          pbsum = pbsum + pbeam(ic,jb)
        end do
      end do
c
c ----------------------------------------------------------------------
c  calculate corrected fluxes at time point theta.
c  calculate energy fluxes (keV/cm**2-s)
c ----------------------------------------------------------------------
c
      onemt = 1.0 - theta
      do j=1,nj-1
        do k=1,nk
          if (itran(k) .ne. 0) then
            du = theta*(u(k,j+1)-u(k,j))+onemt*(usave(k,j+1)-usave(k,j))
            dudr(k,j) = du/dr(j)
            flux(k,j) = 0.0
            do i=1,nk
              flux(k,j) = flux(k,j) - d(k,i,j)*dudr(i,j)
            end do
          end if
        end do
      end do
c
      do j=1,nj-1
        if (iangrot .ne. 0) then
          k       = nion + 4
          angrota = (theta*(u(nk,j)+u(nk,j+1))+onemt*
     .              (usave(nk,j)+usave(nk,j+1)))*0.5
          flxmass = 0.0
          do i=1,nprim
            flxmass     = flxmass+atw(i)*flux(i,j)
          end do
          fluxangc(j) = 0.5 * (r2capi(j)
     .                      +  r2capi(j+1)) * angrota * flxmass * xmassp
          flxangce(j) = 0.5 * angrcple * angrota * kevperg * fluxangc(j)
          if (itran(k) .eq. 1) then
            fluxangv(j) = flux(k,j)
            flux  (k,j) = fluxangc(j)+fluxangv(j)
          else
            fluxangv(j) = flux(k,j)-fluxangc(j)
          end if
          omegapi(j) = angrcple * kevperg * angrota * fluxangv(j)
        end if
        k         = nion + 1
        temm      = theta*u(k,j) + onemt*usave(k,j)
        tepp      = theta*u(k,j+1) + onemt*usave(k,j+1)
        conve(j)  = five_halfs_te * (MAX (fluxe(j), zero)*temm +
     .                     MIN (fluxe(j), zero)*tepp)
        if (itran(k) .eq. 1) then
          conde(j)  = flux(k,j)                ! in this case flux = -Ddu/dr, set above
          flux(k,j) = conde(j)  + conve(j)
        else
          conde(j)  = flux(k,j) - conve(j)     ! analysis mode, assumes flux = total flux
                                               ! from power balance.

        end if
c



        k        = nion + 2
        timm     = theta*u(k,j  ) + onemt*usave(k,j  )
        tipp     = theta*u(k,j+1) + onemt*usave(k,j+1)
        convi(j) = five_halfs_ti * (MAX (fluxi(j), zero)*timm +
     .                    MIN (fluxi(j), zero)*tipp)
        if (itran(k) .eq. 1) then
          condi(j)  = flux(k,j)
          flux(k,j) =  condi(j) +  convi(j)
     .                          + (omegapi(j)+flxangce(j)) * iangrot
        else
          condi(j)  = flux(k,j) -  convi(j)
     .                          - (omegapi(j)+flxangce(j)) * iangrot
        end if
      end do
c
c extrapolate for end points
c
      do k=1,nk
        call extrap (ra(nj-2), ra(nj-1), r(nj),
     .               flux(k,nj-2), flux(k,nj-1), fluxb(k))
      end do
c
      call extrap (ra(nj-2),ra(nj-1),r(nj),conde(nj-2),
     .             conde(nj-1),condeb)
      call extrap (ra(nj-2),ra(nj-1),r(nj),conve(nj-2),
     .             conve(nj-1),conveb)
      call extrap (ra(nj-2),ra(nj-1),r(nj),condi(nj-2),
     .             condi(nj-1),condib)
      call extrap (ra(nj-2),ra(nj-1),r(nj),convi(nj-2),
     .             convi(nj-1),convib)

      call divflx (conde,condeb,hcap,r,ra,drr,1,nj,1,qconde)
      call divflx (conve,conveb,hcap,r,ra,drr,1,nj,1,qconve)
      call divflx (condi,condib,hcap,r,ra,drr,1,nj,1,qcondi)
      call divflx (convi,convib,hcap,r,ra,drr,1,nj,1,qconvi)
c
c     calculate dpedt and dpidt_tot and sum for total.
c

      pdpidt = 0.0
      pdpedt = 0.0
      pdpdt  = 0.0
      if (time .eq. time0)  go to 301
c
      do j=1,nj
        k        = nion + 1
        dpedt(j) = 1.5*(ene(j)*u(k,j)-enesav(j)*usave(k,j))/dtt
        if (itimav .eq. 1)  dpedt(j) =
     .     1.5*(eneav1(j)*uav(k,j)-eneav0(j)*uav0(k,j))*dtsumi
        k        = nion + 2
        sum      = 0.0
        do 605 i=1,nion
          dpidt(j,i)    = (u(i,j)*u(k,j)-usave(i,j)*usave(k,j))/dtt
          if (itimav .eq. 1)
     .      dpidt(j,i)  = (uav(i,j)*uav(k,j)-uav0(i,j)*uav0(k,j))*dtsumi
  605     sum    = sum + dpidt(j,i)
        dpidt_tot(j) = 1.5*sum
      end do
c
      call trapv (r,dpedt,hcap,nj,pdpedt)
      pdpedt = pdpedt*const
      call trapv (r,dpidt_tot,hcap,nj,pdpidt)
      pdpidt = pdpidt*const
      pdpdt  = pdpedt+pdpidt
c
  301 call trapv (r,qrad,hcap,nj,prad)
      prad    = prad*const


      brems_imp  = 0.0D0
      brems_prim = 0.0D0
      DO k = 1,nion
         call trapv (r,brems_nions(1,k),hcap,nj,totl)
         IF( k .LE. nprim) brems_prim = brems_prim + totl
         IF( k .GT. nprim) brems_imp  = brems_imp  + totl
      ENDDO
      brems_prim = brems_prim*const
      brems_imp  = brems_imp*const

 

      call trapv (r,qbeame,hcap,nj,pbeame)
      pbeame  = pbeame*const
      call trapv (r,qbeami,hcap,nj,pbeami)
      pbeami  = pbeami*const
      call trapv (r,qconde,hcap,nj,pcondet)
      pcondet = pcondet*const
      call trapv (r,qconve,hcap,nj,pconvet)
      pconvet = pconvet*const
      call trapv (r,qcondi,hcap,nj,pcondit)
      pcondit = pcondit*const
      call trapv (r,qconvi,hcap,nj,pconvit)
      pconvit = pconvit*const
      call trapv (r,qioni,hcap,nj,pioni)
      pioni   = pioni*const
      call trapv (r,qione,hcap,nj,pione)
      pione   = pione*const
      call trapv (r,qcx,hcap,nj,pcx)
      pcx     = pcx*const
      call trapv (r,qe2d,hcap,nj,pe2d)
      pe2d    = pe2d*const
      call trapv (r,qi2d,hcap,nj,pi2d)
      pi2d    = pi2d*const
      call trapv (r,qsawi,hcap,nj,psawi)
      psawi   = psawi*const
      call trapv (r,qsawe,hcap,nj,psawe)
      psawe   = psawe*const
      paux    = pbeame+pbeami+prfe+prfi
      ptre    = pcondet+pconvet+pione+pe2d+psawe
      ptri    = pcondit+pconvit+pioni+pi2d+psawi+pcx

      palpha  = ptfuse+pbfuse+ptfusi+pbfusi
      denom = pohm+pbsum+prfe+prfi
      if (denom .ne. 0.0)then
        qtrepl  = 5.0 * (ptfuse+ptfusi)/denom
        qbrepl  = 5.0 * (pbfuse+pbfusi)/denom
      else
        qtrepl =0.0
        qbrepl =0.0
      endif
c
      call trapv (r, curden , hcap, nj, totcurden)
      totcurden  = 2.0 * pi * totcurden
      call trapv (r, cur_tor_ps_soln , hcap, nj, totcurpar_ps_dmc)
      totcurpar_ps_dmc = 2.0 * pi * totcurpar_ps_dmc
      call trapv (r, curpar_soln , bsqncap, nj, totcurpar_dmc)
      totcurpar_dmc  = 2.0 * pi * totcurpar_dmc
      totcurpar_dmc  =  totcurpar_dmc  + totcurpar_ps_dmc
      do  j=1,nj
        if(external_beam_cur .eq. 0)then
               curb(j) = curbe(j) + curbi(j) + curbet(j)
        else
               curb(j) = curb_external(j)
        endif
      enddo
c


      call trapv (r,curohm,hcap,nj,totohm)
      call trapv (r,curboot,hcap,nj,totboot)
      call trapv (r,curb,hcap,nj,totb)
      call trapv (r,currf,hcap,nj,totrf)
      totohm  = 2.0 * pi * totohm
      totboot = 2.0 * pi * totboot
      totb    = 2.0 * pi * totb
      totrf   = 2.0 * pi * totrf

      if (curtype .ne. 0) then
c
c       _dmc total currents are diamagnetic corrected  HSJ 8/25/98
c
        call trapv (r, cur_tor_ps_soln , hcap, nj, totcurpar_ps_dmc)
        totcurpar_ps_dmc = 2.0 * pi * totcurpar_ps_dmc
        call trapv (r,curohm,bsqncap,nj,totohm_dmc)
        call trapv (r,curboot,bsqncap,nj,totboot_dmc)
        call trapv (r,curb,bsqncap,nj,totb_dmc)
        call trapv (r,currf,bsqncap,nj,totrf_dmc)
        totohm_dmc = 2.0 * pi * totohm_dmc
        totohm = totohm_dmc
        totboot_dmc = 2.0 * pi * totboot_dmc
        totboot = totboot_dmc + totcurpar_ps_dmc
        totb_dmc = 2.0 * pi * totb_dmc
        totb = totb_dmc
        totrf_dmc = 2.0 * pi * totrf_dmc
        totrf = totrf_dmc
      end if



c
c ----------------------------------------------------------------------
c  write out plot flag; flag = 2 record data for 2-d plots only;
c  otherwise record data for 2-d and 3-d plots
c ----------------------------------------------------------------------
c
      write (ntrplt, '(a)')  ' **** continue ****'

      absflg = IABS (flag)
      write (ntrplt, 8010)  absflg
****  write (6, *)  'absflg in trplot =', absflg
c
c ----------------------------------------------------------------------
c  write out current mesh values; if this is a 2-d run
c  the radial points could be changing.
c ----------------------------------------------------------------------
c
 
      write (ntrplt, 8020)  (r(j), j=1,nj)
c
c ----------------------------------------------------------------------
c ignflg = 1 indicates we have ingnition for the first time
c we will flag this for the plotting routines by setting time and
c avgtim to -time and -avgtim
c ----------------------------------------------------------------------
c
      avgtim = time - (1.0-theta)*dtt
c
c     avgtim is output to plot programs so those quantities
c     that are calculated at the theta pivot time will be
c     correctly represented as functions of time
c
      if ( time .eq. time0)  avgtim =  time
****  if (iflag .eq. -3   )  avgtim = -avgtim
****  if (iflag .eq. -3   )  time   = -time
c
c ----------------------------------------------------------------------
c output 2-d plotting information
c ----------------------------------------------------------------------
c
c ----------------------------------------------------------------------
c the next line of variable involve time derivatives
c and therefore are evaluated at the half time point
c ----------------------------------------------------------------------
c

      avalpbet = 0.0
      bpalpha  = 0.0
      write (ntrplt, 8020) avgtim, voltag, taup, taue, entaue, tauec,
     .                     qtrepl, qbrepl, tauea(1), taupa(1), walp(1),
     .                     enalp(1), avalpbet, bpalpha,
     .                     pohm, prad, paux, ptre, ptri, palpha, pdpdt,
     .                     voltag0, voltag25, voltag50, voltag75,
     .                     voltag90, voltag95,detot,deitot_extra,
     .                     entot, etot,eetot
c
c ----------------------------------------------------------------------
c variables known at the full time point
c ----------------------------------------------------------------------
c
c --- get experimental neutron rate at this time by linear interpolation
c --- from exptl_neutron_rate(1..nbctim)
c
      if      (nbctim .le. 1) then
         exptl_n_rate = exptl_neutron_rate(1)
      else if (ABS (time-bctime(1)) .lt. 5.e-8) then
         exptl_n_rate = exptl_neutron_rate(1)
      else if (ABS (time-bctime(nbctim)) .lt. 5.e-8) then
         exptl_n_rate = exptl_neutron_rate(nbctim)
      else
         call interp(time, bctime, nbctim, exptl_neutron_rate,
     .               exptl_n_rate)
      end if
      exptl_n_rate = MAX (exptl_n_rate, zero)

c
      pradn  = 0.0
      abrate = 0.0
      if(use_nubeam)pwf_tot_source_intg = beam_data%pwf_tot_intg
      IF(use_P_Nfreya)pwf_tot_source_intg = 
     .                      neut_beam%pwf_tot_source_intg
      write (ntrplt, 8021) time, enebar, ene(1), (en(1,k), k=1,nion),
     .                     te(1), ti(1), curden(1), beta, q(1),
     .                     pradn, abrate, qrfe(1), ddntot,
     .                     totcur(1), totohm, totboot, totb, totrf,
     .                     ddbmt, ddnthm, ddknct,
     .                     thermal_thermal_ddntot,
     .                     thermal_thermal_ddptot,
     .                     thermal_thermal_dtntot,
     .                     thermal_thermal_tt2ntot,
     .                     thermal_thermal_he3dptot,
     .                        beam_thermal_ddntot,
     .                        beam_thermal_dtntot,
     .                        beam_thermal_ddptot,
     .                        beam_thermal_tt2ntot,
     .                        beam_beam_ddntot,
     .                        beam_beam_dtntot,
     .                        beam_beam_ddptot,
     .                        beam_beam_tt2ntot,
     .                     dtntot, ttntot, qdd, qdt, qtt, exptl_n_rate,
     .                     ebtot,eatot,pwf_tot_source_intg,
     .                     enbeam_intg
cjmp.ibm       cpu_time = etime( tarray)
      call cpu_time_12(cpu_time) !jmp.ibm      
      write(ntrplt,8021)cpu_time     !04/08/04  HSJ
      write(ntrplt,8021)ssqr_hist    !02/17/05  HSJ
      if(use_nubeam)then
       write (ntrplt, 8021) (beam_data%pinja(k),k=1,beam_data%nbeam)
       print *,'beam_data%pinja ='
       write (6, 8021) (beam_data%pinja(k),k=1,beam_data%nbeam)
      endif
      if(tdba .le. 0) tdba= get_time_dep_beam()                      !function call
      timplta = time
      waveforma = 0.0
      if(tdba .gt. 0)then
c     get time in timplot that is  closest to time
           if(nbeam_pts .gt. 1)then
c              tdifa = 2.*time
c              do ik =1,nbeam_pts
c                 dtaaa = ABS(time - timplot(ik))
c                 tdifa = MIN(tdifa,dtaaa)
c                 if(tdifa .eq. dtaaa)il =ik
c              enddo
c              timplta = timplot(il)
c              waveforma = waveform(il)
              timplta = timplot(nbeam_pts)
              waveforma = waveform(nbeam_pts)
           else if(nbeam_pts .eq. 1)then
              timplta = timplot(1)
              waveforma = waveform(1)
           else
              timplta = time
              waveforma = 0.0
           endif
      endif
       write(ntrplt,8010)nbeam_pts,nplt_max,tdba
       write(ntrplt,8020)timplta,waveforma
       if (n_spec .gt. 0) then
         do   ik=1,nj
           do il=1,n_spec       ! user-requested values of te,ti,rotation at
             ilk = spec_profiles(il)
             if (ilk .eq. ik)   ! radial position ik as function of time
     .         write (ntrplt, 8020)  te(ik), ti(ik),angrot(ik)
           end do
         end do
       end if
c
 
      if (iangrot .eq. 0)  go to 2003
      write (ntrplt, 8020) angrot(1)
 2003 if (wisl .le. 0.0)  go to 2004
      write (ntrplt, 8020)  rsm32, rsp32, rsm21, rsp21
 2004 continue
c      time = ABS (time)
      if (nterow .eq. 0)  go to 2006
      do 2005 i=1,nterow
 2005   holdte(i) = te(jterow(i))
      write (ntrplt, 8020) (holdte(i),i=1,nterow)
 2006 if (jsxr .eq. 0)  go to 2010
      write (ntrplt, 8020) ((sxr(i,k),i=1,ndiode(k)),k=1,narray)
 2010 if (jco2 .eq. 0)  go to 2020
      write (ntrplt, 8020) (denco2(i),i=1,nco2)
 2020 if (jzeff .eq. 0)  go to 2030
      write (ntrplt, 8020) (phzeff(i),i=1,nzeff)
 2030 if (absflg .eq. 2)  return
c
c ----------------------------------------------------------------------
c output 3-d plotting information
c ----------------------------------------------------------------------
c
      write (ntrplt, 8031) (ene(j),zeff(j),te(j),ti(j),curden(j),
     .                      curboot(j),
     .                      curdri(j),currf(j),curohm(j),rbp(j),
     .                      xjbte(j),xjbnf(j), etor(j),q(j), curb(j),
     .                      eta_par_ohml(j),eta_par_tdem_ohml(j),
     .                      curboot_ohml(j),xjbne(j),
     .                      curni_ohml(j),curdrive_ohml(j),
     .                      curohmic_ohml(j),sl31(j),sl32(j), 
     .                      sl34(j),h88l31(j),h88l32(j),
     .                      curden_tdem(j),dpsidt_const_rho_tdem(j),
     .                      dpsidt_const_zeta_tdem(j),
     .                      dpsidrho_const_t_tdem(j),j = 1,nj)




c individual contributions form rf current drive models
      write(ntrplt,FMT='(5(a))')(irfmodel(i),i=1,krf)          !determines active models
      write(ntrplt,8020)(rfmodel_power_e(i),i=1,krf)           !power to electrons
      write(ntrplt,8020)(rfmodel_power_i(i),i=1,krf)           !power to ions
      write(ntrplt,8020)(rfmodel_cd(i),i=1,krf)                !total current drive by model
      do i =1,krf 
         if(irfmodel(i) .ne. no_rf)then
           write(ntrplt,FMT='(a)')irfmodel(i)
           write(ntrplt, 8020)(currfs(j,i),j=1,nj)       !current profile by model
           write(ntrplt, 8020)(jtor_eccd(j,i),j=1,nj)    !eccd toroidal current
                                                         !if model ne eccd then 
                                                         !this will be zero
           qrfes(1:nj,i) = qrfes(1:nj,i)/  0.62415064e16
           write(ntrplt, 8020)(qrfes(j,i),j=1,nj)      !elc heating,watts/cm**3  profile by model
           qrfes(1:nj,i) = qrfes(1:nj,i)* 0.62415064e16
           qrfis(1:nj,i) = qrfis(1:nj,i)/ 0.62415064e16
           write(ntrplt, 8020)(qrfis(j,i),j=1,nj)          !ion heating  profile by model
           qrfis(1:nj,i) = qrfis(1:nj,i)*0.62415064e16
         endif
      enddo
      write (ntrplt, 8020)  (cur_seed(j),j = 1,nj)


c
c  
c       glf23 diagnostic output added 4/15/02 HSJ
c       these arrays may not be allocated if the glf23 model
c       is not active. To avoid a lot of propagation
c       into the plot  codes we allocate and zerro the arrays
c       here if necerssary:
        if( .not. allocated(egamma_m))then    
             allocate (egamma_m(1:nj),STAT = istat)
             if(istat .ne. 0)
     .        call allocate_error("egamma_m,allocate_glf",0,istat)
              egamma_m(:) = 0.0
        endif
        if( .not. allocated(shat_exp))then    
             allocate (shat_exp(1:nj),STAT = istat)
             if(istat .ne. 0)
     .        call allocate_error("shat_exp,allocate_glf",0,istat)
              shat_exp(:) = 0.0
        endif
        if( .not. allocated(alpha_exp))then    
             allocate (alpha_exp(1:nj),STAT = istat)
             if(istat .ne. 0)
     .        call allocate_error("alpha_exp,allocate_glf",0,istat)
              alpha_exp(:) = 0.0
        endif
        if( .not. allocated(gamma_p_m))then    
             allocate (gamma_p_m(1:nj),STAT = istat)
             if(istat .ne. 0)
     .        call allocate_error("gamma_p_m,allocate_glf",0,istat)
              gamma_p_m(:) = 0.0
        endif
        if( .not. allocated(anrate_m))then    
             allocate (anrate_m(1:nj),STAT = istat)
             if(istat .ne. 0)
     .        call allocate_error("anrate_m,allocate_glf",0,istat)
              anrate_m(:) = 0.0
        endif
        if( .not. allocated(gamma_net_i))then    
             allocate (gamma_net_i(1:nj),STAT = istat)
             if(istat .ne. 0)
     .        call allocate_error("gamma_net_i,allocate_glf",0,istat)
              gamma_net_i(:) = 0.0
        endif
        if( .not. allocated(gamma_net_e))then    
             allocate (gamma_net_e(1:nj),STAT = istat)
             if(istat .ne. 0)
     .        call allocate_error("gamma_net_e,allocate_glf",0,istat)
              gamma_net_e(:) = 0.0
        endif
        if( .not. allocated(chie_etg_m))then    
             allocate (chie_etg_m(1:nj),STAT = istat)
             if(istat .ne. 0)
     .        call allocate_error("chie_etg_m,allocate_glf",0,istat)
              chie_etg_m(:) = 0.0
        endif
        if( .not. allocated(anfreq_m))then    
             allocate (anfreq_m(1:nj),STAT = istat)
             if(istat .ne. 0)
     .        call allocate_error("anfreq_m,allocate_glf",0,istat)
              anfreq_m(:) = 0.0
        endif
        if( .not. allocated(exch_m))then    
             allocate (exch_m(1:nj),STAT = istat)
             if(istat .ne. 0)
     .        call allocate_error("exch_m,allocate_glf",0,istat)
              exch_m(:) = 0.0
        endif

        if( .not. allocated(gamma_den))then    
             allocate (gamma_den(1:nj),STAT = istat)
             allocate (gamma_ti(1:nj),STAT = istat)
             allocate (gamma_w(1:nj),STAT = istat)
             if(istat .ne. 0)
     .        call allocate_error("gamma_*,trplot",0,istat)
              gamma_den(:) = 0.0
              gamma_ti(:) = 0.0
              gamma_w(:) = 0.0
        endif
        if( .not. allocated(vexb_den_grad))then    
             allocate (vexb_den_grad(1:nj),STAT = istat)
             allocate (vexb_ti_grad(1:nj),STAT = istat)
             allocate (vexb_w_grad(1:nj),STAT = istat)
             allocate (vexb_tot(1:nj),STAT = istat)
             if(istat .ne. 0)
     .        call allocate_error("vexb_*,trplot",0,istat)
              vexb_den_grad(:) = 0.0
              vexb_ti_grad(:)  = 0.0
              vexb_w_grad(:)   = 0.0
              vexb_tot(:)      = 0.0
        endif
 
        write (ntrplt, 8020) (xketot(i), i=1,nj)
        write (ntrplt, 8020) (xkitot(i), i=1,nj)
        write (ntrplt, 8020) (xchietot(i), i=1,nj)
        write (ntrplt, 8020) (chie_etg_m(i),i=1,nj)
        write (ntrplt, 8020) (xchiitot(i), i=1,nj)
        write (ntrplt, 8020) (xkangrot(i), i=1,nj) 
        write (ntrplt, 8020) (xdchitot(i), i=1,nj)
        write (ntrplt, 8020) (grad_te(i), i=1, nj)

        write (ntrplt, 8020) (egamma_m(i), i=1,nj)
        write (ntrplt, 8020) (gamma_den(i),i=1,nj)
        write (ntrplt, 8020) (gamma_ti(i),i=1,nj)
        write (ntrplt, 8020) (gamma_w(i),i=1,nj)
        write (ntrplt, 8020) (vexb_den_grad(i),i=1,nj)
        write (ntrplt, 8020) (vexb_ti_grad(i),i=1,nj)
        write (ntrplt, 8020) (vexb_w_grad(i),i=1,nj)
        write (ntrplt, 8020) (vexb_tot(i),i=1,nj)
        write (ntrplt, 8020) (shat_exp(i), i=1,nj)
        write (ntrplt, 8020) (alpha_exp(i), i=1,nj)
        write (ntrplt, 8020) (gamma_p_m(i), i=1,nj)
        write (ntrplt, 8020) (anrate_m(i), i=1,nj)
        write (ntrplt, 8020) (gamma_net_i(i),i=1,nj)
        write (ntrplt, 8020) (gamma_net_e(i),i=1,nj)
        write (ntrplt, 8020) (anfreq_m(i), i=1,nj)
        write (ntrplt, 8020) (exch_m(i),   i=1,nj)
c
c------------------------------



      write (ntrplt, 8020)  ((en   (j,i), j=1,nj), i=1,nion)
      write (ntrplt, 8020)  ((xjbni(j,i), j=1,nj), i=1,nion)
      write (ntrplt, 8020)  ((xjbti(j,i), j=1,nj), i=1,nion)
      write (ntrplt, 8020)  ((enn   (j,i), j=1,nj), i=1,2)
c
      if (iangrot .eq. 0)  go to 2035
      write (ntrplt, 8020) (angrot(j), j=1,nj)
c
!      storqueb(:) = -1.e11 ; enbeam(:) = -1.e13 

 2035 do k=1,krf
        if (rfon(k) .le. timmax)  go to 2036
      end do
c
      go to 2100
 2036 write (ntrplt, 8020) (qrfe(j),qrfi(j),j=1,nj)
      if (irfech .eq. 0)  go to 2100
      write (ntrplt, 8010) jresonrf,jrfmin,jrfmax
      write (ntrplt, 8020) (omoder(j),omodei(j),xmoder(j),xmodei(j),
     .                     j = 1,nrfrad)
      write (ntrplt, 8020) (rf1(j),rf2(j),rf3(j),rf4(j),
     .                     j = 1,nrfrad)
c

      
 2100 if (ifus .eq. 0)  go to 2120
      write (ntrplt, 8020) (qfuse(j),qfusi(j),j=1,nj)
      if (iaslow .eq. 0)  go to 2120
      write (ntrplt, 8020) ( tauea(j),j=1,nj)
      write (ntrplt, 8020) ( taupa(j),j=1,nj)
      write (ntrplt, 8020) ((( presprp_12_mcgo(j,ic,ib),j=1,nj),
     .                                      ic=1,ke),ib=1,nbeams)
      write (ntrplt, 8020) ((( prespar_12_mcgo(j,ic,ib),j=1,nj),
     .                                      ic=1,ke),ib=1,nbeams)
      write (ntrplt, 8020) (  walp(j),j=1,nj)
      write (ntrplt, 8020) ( enalp(j),j=1,nj)
 2120 write (ntrplt, 8020) (ddfusn(j),j=1,nj)
      write (ntrplt, 8020) (ddbeam(j),j=1,nj)
      write (ntrplt, 8020) (ddknck(j),j=1,nj)
      write (ntrplt, 8020) (ddnfus(j),j=1,nj)
      write (ntrplt, 8020) (ddpfus(j),j=1,nj)
      write (ntrplt, 8020) (dtnfus(j),j=1,nj)
      write (ntrplt, 8020) (ttnfus(j),j=1,nj)
      write (ntrplt, 8020) (he3dpfus(j),j=1,nj)
      jj = 3*nbeams+1        ! only write out the totals for the beams
      write (ntrplt, 8020) (beam_thermalddn(j,jj),j=1,nj)
      write (ntrplt, 8020) (beam_thermaltth_df(j,jj),j=1,nj)
      write (ntrplt, 8020) (beam_thermalddp(j,jj),j=1,nj)
      write (ntrplt, 8020) (beam_thermaltt2n(j,jj),j=1,nj)
      jj = ((3*nbeams+1)*3*nbeams)/2+1  ! only write out totals
 
      write (ntrplt, 8020) (beam_beamddn(j,jj),j=1,nj)
      write (ntrplt, 8020) (beam_beamdtn(j,jj),j=1,nj)
      write (ntrplt, 8020) (beam_beamddp(j,jj),j=1,nj)
      write (ntrplt, 8020) (beam_beamtt2n(j,jj),j=1,nj)
      write (ntrplt, 8020) (fday2d1(j),j=1,nj)
      write (ntrplt, 8020) (fday2d2(j),j=1,nj)
      write (ntrplt, 8020) (fday2d3(j),j=1,nj)
      write (ntrplt, 8020) (dscrip(j),j=1,nj)
      write (ntrplt, 8020) (qe2d(j),j=1,nj)
      write (ntrplt, 8020) (qi2d(j),j=1,nj)
      write (ntrplt, 8020) (  wbeam(j),j=1,nj)
      write (ntrplt, 8020) ( enbeam(j),j=1,nj)
      write (ntrplt, 8020) ( storqueb(j),j=1,nj)
      write (ntrplt, 8020) ( smagtorque(j),j=1,nj)
c
      call EMPTY (ntrplt)
      return
c
 8010 format ( 6i10)
 8020 format ( 6(1pe14.6))    !needs to match other formats for ntrplt
 8021 format ( 6(1pe16.8))
 8031 format (13(1pe12.4))
c
      end

      subroutine tsplinew (profin, rprof, prof, knots, index,
     .                     xprof, yprof)
c

c
c TSPLINEW interpolates the profile data in profin
c to construct a profile for the present time in prof.
c this routine replaces routine tsplin for the time-dependent
c profile cases when splninpt = 'new'
c
c          Inputs:
c          profin(j,m) - spline knot value j at time = bctime(m).
c                        profin must be dimensioned as below.
c          rprof(j,m)  - normalized radius of knot j,time bctime(m)
c          knots       - knots(m) gives the
c                        number of knots in rprof at time bctime(m)
c          index         pointer to the profile to be interpolated
c                        required to select the proper bpar array in intrp
c          index =1,..kprim   points to primary ion
c          index =kprim+1...kprim+kimp  points to impurity ions
c          index=kprim+kimp+1 points to ene
c                          +2           te
c                          +3           ti
c                          +4           zeff
c                          +5           curden
c                          +6           angrot
c          xprof,yprof - vectors of length nj used for storage
c
c          Output:
c          prof(j)     - profile value at the present time at
c                        r(j)  (j = 1,nj).
c
c --- INCLUDE file bcondspl.i is required here so that the time index
c --- mtimeprf and pointer iprofnbr are properly passed to intrp
c
      USE param
      USE io
      USE solcon
      USE numbrs
      USE mesh
      USE bd_condtn,only : bparenp,bpareni,bparte, bparti,
     .                     bparang, bparene,bparzeff, bparcur,
     .                     bparkpol, bparcurb_external,
     .                     iprofnbr, mtimeprf, knoterror,
     .                     profiles_bcondspl,
     .                     bctime,ub,fluxb,
     .                     ub_save,ub_rho_edge,bctime_zone

      implicit  integer (i-n), real*8 (a-h, o-z)

c
      dimension  profin(ksplin,kbctim), rprof(ksplin,kbctim),
     .           prof(kj), xprof(*), yprof(*), knots(*)
c      data       time_tol/1.0e-9/ defined in solcon
c
      iprofnbr  = index
      knoterror = 0
c


      do m=2,nbctim
        mtimeprf  = m
        if (bctime(m) .ge. time)  go to 2120
      end do
    
c
c      if (ABS (bctime(nbctim)-time) .lt. 2.0*time_tol) then
      if (ABS (bctime(nbctim)-time) .lt. dtmin) then
c
c --- time is larger than timmax, due to effect of time_tol
c
          m      = nbctim
          nprofm = knots(m)
          if (nprofm .lt. 3) then
            write (6, *)  'ERROR in TSPLINEW, number of knots =', nprofm
            write (6, *)  '      nbctim=', nbctim
            write (6, *)  '        time=', time
            call STOP ('subroutine TSPLINEW: problem #1', 88)
          end if
c
c --- fit a spline at time point m
c
          call intrp (-1, -1, rprof(1,m), profin(1,m),
     .                 nprofm, roa, prof, nj)
          return
      else
          knoterror = -1
      end if
 2120 m = mtimeprf
 
c
      if (            knoterror .ne. 0  )  go to 1510
      if (                  m-1 .lt. 1  )  go to 1510
      if (bctime(m)-bctime(m-1) .le. 0.0)  go to 1510
c
      ft = (time-bctime(m-1))/(bctime(m)-bctime(m-1))
c
c --- get number of knots
c
      nprofm = knots(m)
      if (   nprofm .lt. 3           )  knoterror = 1
      if (knoterror .ne. 0           )  go to 1510
      if (    index .eq. kprim+kimp+6)  go to 110
      do j=1,nprofm
        if (profin(j,m) .gt. 0.0)  go to 110
      end do
      knoterror = 3
      go to 1510
c
c --- fit a spline at time point m
c
  110 call intrp (-1, -1, rprof(1,m), profin(1,m), nprofm,
     .             roa, yprof, nj)
c
c --- get number of knots
c
      nprofm = knots(m-1)
      if (   nprofm .lt. 3           )  knoterror = 2
      if (knoterror .ne. 0           )  go to 1510
      if (    index .eq. kprim+kimp+6)  go to 120
      do j=1,nprofm
        if (profin(j,m-1) .gt. 0.0)  go to 120
      end do
      knoterror = 3
      go to 1510
c
c --- fit a spline at time point m-1
c
  120 mtimeprf = mtimeprf-1
      call intrp (-1,-1,rprof(1,m-1),profin(1,m-1),nprofm,
     .             roa,xprof,nj)
c
c --- interpolated the two profiles to the current time
c
      do j=1,nj
        prof(j) = xprof(j) + ft*(yprof(j)-xprof(j))
      end do
      
c      IF(iprofnbr == kprim+kimp+1 )THEN
c       print *,'===============debug=========cray401.f line 9930'
c       print *,'m-1,m,bctime(m-1),bctime(m) =',m,m-1
c       print *,'bctime(m-1),bctime(m)=',bctime(m-1),bctime(m)
c       print*,'time,ft =',time,ft
c       print*,'yprof =',yprof(1:25)
c       print*,'xprof =',xprof(1:25)
c      ENDIF

      return
c
 1510 if (nbctim .eq. 1) then
c
c       only a single boundary condition time was given.
c       this means that the profiles are constant in time
c
        nprofm = knots(nbctim)
        if (nprofm .lt. 3) then
          write (6, *)  'ERROR in TSPLINEW, number of knots =', nprofm
          write (6, *)  '      nbctim = ', nbctim
          write (6, *)  '      time   = ', time
          call STOP ('subroutine TSPLINEW: problem #2', 255)
        end if
        call intrp (-1, -1, rprof(1,1), profin(1,1),
     .                      nprofm, roa, prof, nj)
      else
        write  (ncrt, 1500)  knoterror, iprofnbr
 1500   format (/ ' subroutine TSPLINEW reports:' /
     .            '   knoterror ='       , i5     /
     .            '   for profile index ', i5)
        if (knoterror .gt. 0) then
           write (ncrt,'("fewer than 3 knots were found for" /
     .                   " profile ",a," at time ",1pe12.6)')
     .                     profiles_bcondspl(index),time
           call STOP ('subroutine TSPLINEW: problem #3', 89)
c
        else if (knoterror .lt. 0) then
          write  (ncrt, 1520)  m, time, (bctime(i),i=1,nbctim)
 1520     format ('  m = ', i5, '  time = ', 1pe16.10 /
     .            '  bctime(1 to nbctim)'             /
     .               (5(2x,1pe16.10)))
           call STOP ('subroutine TSPLINEW: problem #4', 248)
        end if
      end if
      return
c
      end

      subroutine tsplin (profin, rprof, nprof, prof)
c
      USE param
      USE solcon
      USE numbrs
      USE mesh
      USE bd_condtn,only : bctime,ub,fluxb,
     .    ub_save,ub_rho_edge,bctime_zone
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c TSPLIN interpolates the profile data in profin
c to construct a profile for the present time in prof.
c Linear interpolation is performed in time (if the
c profile is time-dependent), and cubic spline
c interpolation is performed in space.
c
c          Inputs:
c          profin(j,m) - spline knot value j at time = bctime(m).
c                        profin must be dimensioned as below.
c                        if profin(1,2) = 0, the profile is
c                        considered to be time-independent.
c          rprof(j)   -  normalized radius of knot j.
c          nprof      -  number of knots.
c
c          Output:
c          prof(j)    -  profile value at the present time at
c                        r(j)  (j = 1,nj).
c
c      include 'param.i'
c      include 'bcon.i'
c      include 'mesh.i'
c      include 'numbrs.i'
c      include 'solcon.i'
c
      dimension profin(ksplin,kbctim),rprof(nprof),prof(kj)
      dimension splint(ksplin)
c
      if (profin(1,2) .eq. 0.0)  go to 2140
c
c          The profile is time-dependent, so we do a
c          linear interpolation in time to get the profile.
c          (roa is a normalized r mesh set up by rhomsh.)
c
      do m=2,nbctim
        if (bctime(m) .gt. time)  go to 2120
      end do
 2120 ft = (time-bctime(m-1))/(bctime(m)-bctime(m-1))
      do 2130 j=1,nprof
 2130 splint(j) = profin(j,m-1) + ft*(profin(j,m)-profin(j,m-1))
      call intrp (0, 1, rprof, splint, nprof, roa, prof, nj)
      return
c
c     the profile is not time-dependent
c
 2140 call intrp (0, 1, rprof, profin, nprof, roa, prof, nj)
      return
c
      end

      real*8 function ttnrate (tikevc)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ------------------------------------------------- 11/20/95 --- HSJ ---
c function returns the t(t,2n)he4 reaction rate,cm**3/sec (averaged
c over a Maxwellian distribution at temperature tikev)
c This subroutine is a simple linear itnterpolaton of measured? and calc.?
c reaction rates "Atomic Data for Controled Fusion Research",ORNL-5207,
c E.W. Thomas et al. .NOV. 1979 for energies ge 6keV
c for the lower energies the NRL formulary table was used
c                         NOTE
c      THERE APPEARS TO BE CONSIDERABLE UNCERTAINTY IN T(T,2N)HE4 REACTION
c      RATES SO THESE VALUES MAY BE QUITE CRUDE.
c input
c   tikev       ion temperature in keV
c ----------------------------------------------------------------------
c
      parameter (ndt = 12)
c
      integer    ndtlo, ndtp1, i
      real*8     ekev_table(ndt), rttrate_table(ndt), tikev,
     .           tikevc
c
      data    ndtlo/0/ ! ndtlo should be saved between calls
c                        data statement just sets it for very first call
c
      data   (ekev_table(i), rttrate_table(i), i=1,ndt)/
     .               0.0,     0.0,
     .               1.0,     3.3e-22,
     .               2.0,     7.1e-21,
     .               5.0,     1.4e-19,
     .               6.0,     2.1e-19,
     .              10.0,     6.8e-19,
     .              20.0,     2.4e-18,
     .              40.0,     6.6e-18,
     .              70.0,     1.3e-17,
     .             100.0,     2.0e-17,
     .             200.0,     4.1e-17,
     .             400.0,     7.4e-17
     .                                 /
c

       tikev = MIN(tikevc, ekev_table(ndt))
      if (tikev .le. ekev_table(ndt)) then         ! find tikev in table
        call tableintrp (ekev_table, ndt, tikev, ndtlo)
        if (ndtlo .le. 0 .or. ndtlo .gt. ndt-1) then
          call STOP ('function TTNRATE: bugcheck', 69)
        else
          ndtp1   = ndtlo + 1
          dele    = ekev_table(ndtp1)-ekev_table(ndtlo)
          delr    = rttrate_table(ndtp1)-rttrate_table(ndtlo)
          ttnrate = rttrate_table(ndtlo)+
     .              (delr/dele)*(tikev-ekev_table(ndtlo))
        end if
      else
        print *,'tikev =',tikev
        print *,'table range =',ekev_table(ndt)
        call STOP ('function TTNRATE: TI out of range', 70)
      end if
      return
c
      end

      subroutine tweak1
c

c
c  this subroutine tweaks the following parameters:
c
c     wneo(3,3)    to get a specified fusion neutron rate (fusnin)
c     OR w3typt    OR a specified central ion temperature (ticin),
c                 (w3typt>0 takes precedence over wneo(3,3) and
c                  fusnin>0 takes precedence over ticin)
c
c     zeffc,zeffb  to get a specified one-turn voltage (voltin)
c                  at the plasma surface
c                  OR a specified average voltage (voltav)
c                  across the plasma,
c                  (voltav>0 takes precedence over voltin)
c
c     zeffc/zeffb  to get a specified q on axis (qcin).
c
      USE param
      USE fusion
      USE io
      USE ions
      USE solcon
      USE soln
      USE extra
      USE numbrs
      USE mesh
c
      USE machin
      USE tfact
      USE geom
      USE flags
      USE bd_condtn
      implicit  integer (i-n), real*8 (a-h, o-z)
c      include 'param.i'
c      include 'extra.i'
c      include 'flags.i'
c      include 'fusion.i'
c      include 'geom.i'
c      include 'invers.i'
c      include 'io.i'
c      include 'ions.i'
c      include 'machin.i'
c      include 'mesh.i'
c      include 'numbrs.i'
c      include 'solcon.i'
c      include 'soln.i'
c      include 'tfact.i'
c
      data ifirst /1/
      save ifirst
c
      one_point_01 =   1.01
      hundred      = 100.0
c
c     set w3
c
      w3 = wneo(3,3)
      if (w3typt .gt. 0.0)  w3 = w3typt
c
c     write heading for printout; save initial values
c
      if (ifirst .eq. 0)  go to 10
      ifirst = 0
      write (ntweak, 1110)
      timold = time-dt
      fusold = ddntot
      tiold  = ti(1)
      vold   = voltag
      qcold  = q(1)
      w33old = w3
      zcold  = zeff(1)
      zbold  = zeff(nj)
      zrold  = zcold/zbold
      w33zo  = w33old*zcold
c
c  write tweak status and return if time-timold .lt. ttweak
c
   10 if (time-timold .ge. ttweak-1.0e-6)  go to 20
      write (ntweak, 1120) time, ddntot, ti(1), voltag, q(1)
      return
c
   20 timold = time
c
c  tweak zeffc and zeffb to get specified average voltage: voltav
c
      if (  iten .eq. 1 .or. nimp .eq. 0)  go to 540
      if (voltav .eq. 0.0               )  go to 340
      if (voltav .eq. voltoh            )  go to 330
      itweak = 2
      zc     = zeff(1)
      zb     = zeff(nj)
      if (zc .ne. zcold .and. zb .ne. zbold)  go to 310
      zc     = (voltav/voltoh)*zc
      zb     = (voltav/voltoh)*zb
      go to 320
  310 if (voltoh .eq. vold)  go to 440
      zcn    = zcold + (zc-zcold)*(voltav-vold)/(voltoh-vold)
      if (zcn/zc .gt. 2.0)  zcn = 1.5*zc
      if (zcn/zc .lt. 0.5)  zcn = 0.7*zc
      zb     = zb*zcn/zc
      zc     = zcn
  320 zzzmax = z(1,nprim+1)-0.5
      if (zeflim .gt. 0.0)  zzzmax = zeflim
      zc     = MIN (zc, zzzmax      )
      zb     = MIN (zb, zzzmax      )
      zc     = MAX (zc, one_point_01)
      zb     = MAX (zb, one_point_01)
      zeffc(1) = zc
      zeffb(1) = zb
  330 vold     = voltoh
      go to 440
c
c  OR tweak zeffc and zeffb to get specified surface voltage: voltin
c
  340 if (voltin .eq. 0.0   )  go to 440
      if (voltag .eq. voltin)  go to 430
      itweak = 2
      zc     = zeff(1)
      zb     = zeff(nj)
      if (zc .ne. zcold .and. zb .ne. zbold)  go to 410
      zc     = (voltin/voltag)*zc
      zb     = (voltin/voltag)*zb
      go to 420
  410 if (voltag .eq. vold)  go to 440
      zcn    = zcold + (zc-zcold)*(voltin-vold)/(voltag-vold)
      if (zcn/zc .gt. 2.0)  zcn = 1.5*zc
      if (zcn/zc .lt. 0.5)  zcn = 0.7*zc
      zb     = zb*zcn/zc
      zc     = zcn
  420 zzzmax = z(1,nprim+1) - 0.5
      if (zeflim .gt. 0.0)  zzzmax = zeflim
      zc     = MIN (zc, zzzmax      )
      zb     = MIN (zb, zzzmax      )
      zc     = MAX (zc, one_point_01)
      zb     = MAX (zb, one_point_01)
      zeffc(1) = zc
      zeffb(1) = zb
  430 vold     = voltag
c
c  tweak zeffc/zeffb to get specified q on axis
c
  440 if (qcin .eq. 0.0                     )  go to 540
      if (qcin .ne. 0.0 .and. q(1) .eq. qcin)  go to 530
      itweak = 2
      zrat = zeff(1)/zeff(nj)
      zratp = zrat
      if (zrat .ne. zrold)  go to 510
      if (q(1) .gt. qcin) zrat = 0.9*zrat
      if (q(1) .lt. qcin) zrat = 1.1*zrat
      go to 520
  510 if (q(1) .eq. qcold)  go to 540
      zratn = zrold + (zrat-zrold)*(qcin-qcold)/(q(1)-qcold)
      if (zratn/zrat .gt. 2.0) zratn = 2.0 * zrat
      if (zratn/zrat .lt. 0.5) zratn = 0.5*zrat
      zrat = zratn
  520 zzzmax = z(1,nprim+1)-0.5
      if (zeflim .gt. 0.0) zzzmax = zeflim
      zrat = MIN (zrat,zzzmax)
      zrat = MAX (zrat,1.0/zzzmax)
      zb = SQRT (zratp/zrat)*zeffb(1)
      zb = MIN (zb, zzzmax      )
      zb = MAX (zb, one_point_01)
      zc = zrat * zb
      if (zc .le. zzzmax)  go to 522
      zc = zzzmax
      zb = zc/zrat
      go to 524
  522 if (zc .ge. 1.01)  go to 526
      zc = 1.01
      zb = zc/zrat
  524 zb = MIN (zb, zzzmax      )
      zb = MAX (zb, one_point_01)
  526 zeffc(1) = zc
      zeffb(1) = zb
  530 qcold    = q(1)
c
c  tweak w33, which is either wneo(3,3) or w3typt, to get specified
c     fusion neutron rate or central ion temperature.
c     Note: when itweak = 2, Zeffc has been tweaked, so we compute
c     the new product w33*zeffc, and then divide by the new zeffc.
c
  540 w33  = w3
      w33z = w3*zeff(1)
      if (fusnin .eq. 0.0)  go to 730
c
c     tweaking to obtain fusnin
c
      if (itweak .eq. 2)  go to 720
c
c     zeff is not being tweaked
c
      if (ddntot .eq. fusnin)  go to 770
      itweak = 1
      if (w33 .ne. w33old)  go to 710
      if (ddntot .gt. fusnin) w33 = 1.5*w33
      if (ddntot .lt. fusnin) w33 = w33/1.5
      go to 760
  710 if (ddntot .eq. fusold)  go to 780
      w33 = w33old + (w33-w33old)*(fusnin-fusold)/(ddntot-fusold)
      if (ddntot .lt. fusnin .and. w3 .eq. w33min) w33 = w33min
      go to 760
c
c          zeff is being tweaked
c
  720 if (ddntot .ne. fusnin)  go to 722
      w33 = w33z/zeffc(1)
      go to 760
  722 if (w33z .ne. w33zo)  go to 724
      if (ddntot .gt. fusnin) w33z = 1.5*w33z
      if (ddntot .lt. fusnin) w33z = w33z/1.5
      w33 = w33z/zeffc(1)
      go to 760
  724 if (ddntot .eq. fusold)  go to 780
      w33zn = w33zo + (w33z-w33zo)*(fusnin-fusold)/(ddntot-fusold)
      if (w33zn/w33z .gt. 2.0) w33zn = 1.5*w33z
      if (w33zn/w33z .lt. 0.5) w33zn = 0.7*w33z
      w33z = w33zn
      w33 = w33z/zeffc(1)
      if (ddntot .lt. fusnin .and. w3 .eq. w33min) w33 = w33min
      go to 760
 730  if (ticin .eq. 0.0)  go to 780
c
c             tweaking to obtain ticin
c
      if (itweak .eq. 2)  go to 750
c
c          zeff is not being tweaked
c
      if (ti(1) .eq. ticin)  go to 770
      itweak = 1
      if (w33 .ne. w33old)  go to 740
      if (ti(1) .gt. ticin) w33 = 1.5*w33
      if (ti(1) .lt. ticin) w33 = w33/1.5
      go to 760
  740 if (ti(1) .eq. tiold)  go to 780
      w33 = w33old + (w33-w33old)*(ticin-tiold)/(ti(1)-tiold)
      if (ti(1) .lt. ticin .and. w3 .eq. w33min) w33 = w33min
      go to 760
c
c          zeff is being tweaked
c
 750  if (ti(1) .ne. ticin)  go to 752
      w33 = w33z/zeffc(1)
      go to 760
 752  if (w33z .ne. w33zo)  go to 754
      if (ti(1) .gt. ticin) w33z = 1.5*w33z
      if (ti(1) .lt. ticin) w33z = w33z/1.5
      w33 = w33z/zeffc(1)
      go to 760
  754 if (ti(1) .eq. tiold)  go to 780
      w33z = w33zo + (w33z-w33zo)*(ticin-tiold)/(ti(1)-tiold)
      w33 = w33z/zeffc(1)
      if (ti(1) .lt. ticin .and. w3 .eq. w33min) w33 = w33min
 760  w33old = w3
      w33zo = w3*zeff(1)
      w33 = MAX (w33, w33min )
      w33 = MIN (w33, hundred)
      if (w3typt .eq. 0.0)  wneo(3,3) = w33
      if (w3typt .gt. 0.0)  w3typt    = w33
 770  fusold = ddntot
      tiold  = ti(1)
c
c  save zcold, zbold, and zrold; calculate new zeff profile
c     and adjust ion densities
c
  780 if (itweak .le. 1)  go to 920
      zcold = zeff(1)
      zbold = zeff(nj)
      zrold = zcold/zbold
      do j=1,nj
        xsc = 1.0 - (r(j)/r(nj))**gamzef(1)
        if (xsc .lt. 0.0) xsc = 0.0
        zeff(j) = zeffb(1) + (zeffc(1)-zeffb(1))*xsc**alpzef(1)
      end do
      call zen
c
c  write tweak status
c
  920 if (itweak .ne. 0)  go to 930
      write (ntweak,1120) time, ddntot, ti(1), voltag, q(1)
      return
  930 write (ntweak, 1120) time, ddntot, ti(1), voltag, q(1),
     .                     wneo(3,3), zeffc(1), zeffb(1)
      itweak = 0
      return
c
 1110 format (/ 6x,'time',9x,'fusn',5x,'ti(1)',4x,'voltag',6x,'q(1)',
     .          7x,'w33',5x,'zeffc',5x,'zeffb' /
     .          7x,'(s)',8x,'(s-1)',5x,'(keV)',7x,'(V)')
 1120 format (f10.4,e13.3,6f10.4)
c
      end

      subroutine tweak2 (dtecal, dtemix, dtical, dtimix, epste, epsti,
     .                   fuscal, fusmix, s3cal, s3mix, s22cal, s22mix,
     .                   s18cal, s18mix, trcal, trmix)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c  This subroutine tweaks the following parameters:
c     epste  to get a specified dtemix, s3mix, s22mix, or s18mix,
c     epsti  to get a specified dtimix, fusmix, or trmix.
c
      data ifirst /1/
c
c  save initial values
c
      if (ifirst .eq. 0)  go to 10
      ifirst = 0
      dteold = dtecal
      s3old  = s3cal
      s22old = s22cal
      s18old = s18cal
      dtiold = dtical
      fusold = fuscal
      trold  = trcal
      epsteo = epste
      epstio = epsti
c
c  tweak epste to get specified dtemix
c
   10 if (dtemix .eq. 0.0)  go to 200
      call tweakx (epste,epsteo,dtemix,dtecal,dteold,
     .             0.01,dtemix,1.5,2.0)
      go to 400
c
c  tweak epste to get specified s3mix
c
  200 if (s3mix .eq. 0.0)  go to 300
      call tweakx(epste,epsteo,s3mix,s3cal,s3old,
     .            0.01,10.0,1.5,2.0)
      go to 400
c
c  tweak epste to get specified s22mix
c
  300 if (s22mix .eq. 0.0)  go to 350
      call tweakx(epste,epsteo,s22mix,s22cal,s22old,
     .            0.01,10.0,1.5,2.0)
      go to 400
c
c  tweak epste to get specified s18mix
c
  350 if (s18mix .eq. 0.0)  go to 400
      call tweakx(epste,epsteo,s18mix,s18cal,s18old,
     .            0.01,10.0,1.5,2.0)
c
c  tweak epsti to get specified dtimix
c
  400 if (dtimix .eq. 0.0)  go to 500
      call tweakx(epsti,epstio,dtimix,dtical,dtiold,
     .            0.01,dtimix,1.5,2.0)
      return
c
c  tweak epsti to get specified fusmix
c
  500 if (fusmix .eq. 0.0)  go to 600
      call tweakx(epsti,epstio,fusmix,fuscal,fusold,
     .            0.01,10.0,1.5,2.0)
      return
c
c  tweak epsti to get specified trmix
c
  600 if (trmix .eq. 0.0)  return
      call tweakx(epsti,epstio,trmix,trcal,trold,
     .            0.01,10.0,1.5,2.0)
      return
c
      end

      subroutine tweakx (x, xold, fspec, fnow, fold,
     .                   xmin, xmax, xnudge, xrat)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c          TWEAKX tweaks x to obtain f(x) = fspec.
c          Inputs are such that  f(x) = fnow,  f(xold) = fold.
c          x is constrained to  xmin <= x <= xmax, and changes
c          in x in any one tweak are limited to a factor of xrat.
c          xnudge is the factor for the initial change in x.
c
c          The new value of x is returned, and xold and fold are
c          updated with the input x and fnow.
c
      xnow = x
      if (fnow .eq. fspec)  go to 130
      if (xnow .ne. xold )  go to 110
      if (fnow .gt. fspec) x = xnow*xnudge
      if (fnow .lt. fspec) x = xnow/xnudge
      go to 120
  110 x = xold + (xnow-xold)*(fspec-fold)/(fnow-fold)
      if (x*xnow .eq. 0.0)  go to 120
      if (x*xnow .lt. 0.0) x = xnow/xrat
      if (x/xnow .gt. xrat) x = xnow*xrat
      if (xnow/x .gt. xrat) x = xnow/xrat
  120 x = MAX (x,xmin)
      x = MIN (x,xmax)
  130 xold = xnow
      fold = fnow
      return
c
      end

      subroutine uertst1 (ier, obsolete)
c --- changed name from UERTST to UERTST1 so that only IMSL routines
c --- with sources in ONETWO will call UERTST1. Other IMSL routines used
c --- in ONETWO which are linked from the IMSL library will thus not
c --- be involved with this subroutine.
c
c --- uertst1 ---------------- library 2 -------------------------------
c
c   function            - error message generation
c   usage               - call uertst1 (ier, obsolete)
c   parameters   ier    - error parameter. type + n  where
c                           type = 128 implies terminal error
c                                   64 implies warning with fix
c                                   32 implies warning
c                              n = error code relevant to calling routine
c              obsolete - input scalar (double precision on dec)
c                         containing the name of the calling routine
c                         as a 6-character literal string. --- OBSOLETE
c   language            - Fortran
c                         DEC
c
c ----------------------------------------------------------------------
c
c   latest revision     - october 1, 1975
c
c
      USE param
      USE io
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c      include 'param.i'
      include 'imsl.i'
c      include 'io.i'
c
      integer        warn, warf, term
      dimension      ibit(4)
      character*(*)  obsolete
      equivalence   (ibit(1), warn), (ibit(2), warf), (ibit(3), term)
      data           ibit / 32, 64, 128, 0 /
c
      ier2 = ier
      if (ier2 .ge. warn)  go to 5
c
c     non-defined
c
      ier1 = 4
      go to 20
    5 if (ier2 .lt. term)  go to 10
c
c     terminal
c
      ier1 = 3
      go to 20
   10 if (ier2 .lt. warf)  go to 15
c
c     warning (with fix)
c
      ier1 = 2
      go to 20
c
c     warning
c
   15 ier1 = 1
c
c     extract 'n'
c
   20 ier2 = ier2 - ibit(ier1)
c
c     print error message
c     disable undefined error output (it's not an error)
c
      if (ier1 .eq. 4)  return
c
      write  (nout, 26)  imslmd, ier
      write  (nqik, 26)  imslmd, ier
      write  (ncrt, 26)  imslmd, ier
   26 format (/ ' ERROR message from ', a,
     .          ':  ier (error code) =', i5)
      return
c
      end

      subroutine updatb (enb, kj, ke, nj, nbeams, const, enbeam)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c --- UPDATB sums over beam and energy component at each point
c
      dimension  enb(kj,ke,*), enbeam(kj)
c
      do j=1,nj
        enbeam(j) = 0.0
        do ib=1,nbeams
          do ic=1,3
            enbeam(j) = enbeam(j) + const*enb(j,ic,ib)
          end do
        end do
      end do
c     make sure gradient is zero at the magnetic axis (the
c     following is second order accurate in grid spacing)
c     (important for p-prime and ffprime for mhd calcs) HSJ 12/4/99
c
      enbeam(1)=(4.*enbeam(2)-enbeam(3))/3.
c
      return
c
      end

      subroutine update (u, en, te, ti, rbp, nk, nj, kj, kk,
     .                   iangrot, angrot)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c this subroutine updates the dependent variables (en, te, ti, rbp and
c angrot) by setting them equal to the appropriate elements of the
c current solution vector u.c
      dimension u(kk,*), en(kj,*), te(*), ti(*), rbp(*), angrot(*)
c
      do j=1,nj
        do k=1,nk-3-iangrot
          en (j,k) = u(k,j)
        end do
        te (j  ) = u(nk-2-iangrot,j)
        ti (j  ) = u(nk-1-iangrot,j)
        rbp(j  ) = u(nk  -iangrot,j)
        if (iangrot .eq. 1)  angrot(j) = u(nk,j)
      end do
      return
c
      end

      subroutine vecprod (veca, vecb, vecout, nj, itype)
c
      implicit none
c
c ----------------------------------------------------------------------
c     if itype .eq. 0
c       return vecout(j)=(outer product) veca(j)*vecb(j)   for j=1,..nj
c       vecout and/or veca or vecb can point to same storage
c     if itype .ne. 0
c       return vecout(1)= (inner product)  veca dot vecb
c -------------------------------------------------------- HSJ ---------
c
      integer nj,j,itype
      real*8  veca(*),vecb(*),vecout(*),sum
c
      if (itype .eq. 0) then
        do j=1,nj
          vecout(j)=veca(j)*vecb(j)
        end do
      else
        sum=0.0
        do j=1,nj
          sum=sum+veca(j)*vecb(j)
        end do
        vecout(1)=sum
      end if
      return
c
      end

      subroutine wprof
c
      USE param
      USE numbrs
      USE tfact
      USE bd_condtn
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c --- WPROF interpolates for a factor profile by which
c --- to multiply the 'TYP' electron thermal conductivity
c
c      include 'param.i'
      include 'imsl.i'
c      include 'invers.i'
c      include 'numbrs.i'
c      include 'tfact.i'
c
      data ifirst /1/
c
      imslmd = 'wprof  '
      if (nw1pro .eq. 0)  go to 2
      call tsplin (w1pro, rnormin, nw1pro, w1fact)
    2 if (nw2pro .eq. 0)  go to 3
      call tsplin (w2pro, rnormin, nw2pro, w2fact)
    3 if (nw3pro .eq. 0)  go to 4
      call tsplin (w3pro, rnormin, nw3pro, w3fact)
    4 if (nw4pro .eq. 0)  go to 5
      call tsplin (w4pro, rnormin, nw4pro, w4fact)
    5 if ( nvpro .eq. 0)  return
      call tsplin ( vpro, rnormin,  nvpro, vfact )
      ifirst = 0
      return
c
      end

      subroutine zefcal
c
      USE param
      USE ions
      USE soln
      USE numbrs
      USE mesh
      USE machin
      USE geom
      USE zeffcom
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c          ZEFCAL calculates the line-integrated emission of visible
c          continuum light due to electron bremsstrahlung radiation
c          along several tangential chords.  This simulates the
c          measurements of the z-effective profile diagnostic.
c          The output is in units of photons/s incident on detector.
c

c
      data pi/3.1415926/
c
      brems(dn,zef,tek) =
     .  0.95e-13*3.0*dn**2*zef * EXP (-12.4/(zeffwl*tek))*zeffdl/
     .  (zeffwl * SQRT (tek * 1000.0))
c
      call zeroa (phzeff,nzeff)
      do 200 izef=1,nzeff
      if (rtzeff(izef) .ge. rmax)  go to 200
c
c          calculate start and end points for line integral
c
      x0 = rtzeff(izef)
      y0 =  - SQRT (rmax**2 - x0**2)
      z0 = 89.0
      xf = x0
      yf = 0.0
      zf = z0
      if (rtzeff(izef) .gt. rin)  go to 100
      yf =  - SQRT (rin**2 - xf**2)
c
c     initialize line integral; use 1 cm step size
c
 100  nc = yf - y0 + 1.0
      if (nc .lt. 30) nc = 30
      dy = (yf-y0)/(nc-1.0)
      delvol = dy * pi * zefrad**2
      y1 = y0
      call getrho(x0,y1,z0,rho1)
      call interp(rho1,r,nj,ene,ene1)
      call interp(rho1,r,nj,zeff,zeff1)
      call interp(rho1,r,nj,te,te1)
      phot1 = brems(ene1,zeff1,te1)
      sum = 0.0
c
c          evaluate line integral
c
      do 120 i=2,nc
      y2 = y1 + dy
      call getrho(x0,y2,z0,rho2)
      call interp(rho2,r,nj,ene,ene2)
      call interp(rho2,r,nj,zeff,zeff2)
      call interp(rho2,r,nj,te,te2)
      phot2 = brems(ene2,zeff2,te2)
      sum = sum + 0.5*(phot1+phot2)*delvol
      y1 = y2
      phot1 = phot2
  120 continue
      phzeff(izef) = 2.0 * sum
  200 continue
      return
c
      end

      subroutine zen
c
      USE nrtype,                   ONLY : I4B,DP
      USE param
      USE fusion
      USE io
      USE ions
      USE nub
      USE nub2
      USE solcon
      USE soln
      USE numbrs
      USE rf
      USE verbose
      USE sourc
      USE flags
      USE bd_condtn
      USE zen_state,               ONLY : zavold,set_dzdt,
     .                                        en_calc,get_charge_state
      USE mesh,                    ONLY : fix_edge_ni,roa
      implicit  integer (i-n), real*8 (a-h, o-z)
      REAL(DP) sumd(5)
      LOGICAL set_u
c
c     this subroutine calculates various densities, charge numbers,
c     and derivatives of charge numbers.
c

c


c
c ----------------------------------------------------------------------
c  calculate the following quantities for impurity ions:
c     density
c     average charge number
c     average square of charge number
c     temperature derivative of average charge number
c     time derivative of average charge number
c ----------------------------------------------------------------------
c
      one      = 1.0
      izenstop = 0
c
      if (zenvb .gt. 0) then
          write (*,'(// " on entry to ZEN we have time, nprim,inenez =",
     .                       f14.6,2x, i5,2x,i5)')time, nprim,inenez
          do i=1,nprim
             write (*, '(" i = ", i3,
     .                   " en(1:5,i) = ",5(1pe14.4))')i,(en(j,i),j=1,5)
          end do
          write (*,'(" on entry to ZEN we have nimp =", i5)') nimp
          if (nimp .gt. 0) then
            do i=1,nimp
               write (*,'(" i = ",i3,
     .                    " en(1-5,nprim+i) = ",5(1pe14.4))')i,
     .                     (en(j,nprim+i),j=1,5)
            end do
          end if
          write (*,'(" enalp (1:5) =  ", 5(1pe14.4))') (enalp (j),j=1,5)
          write (*,'(" enbeam(1:5) =  ", 5(1pe14.4))') (enbeam(j),j=1,5)
          write (*,'(" ene   (1:5) =  ", 5(1pe14.4))') (ene   (j),j=1,5)
          do j=1,5
            sumd(j)=enbeam(j)
            do i=1,nion
              sumd(j)=sumd(j)+en(j,i)
            end do
          end do
          write (*,'(" total ion density at r(1:5) =" /
     .                 5(1pe14.4))')(sumd(j),j=1,5)
      end if
c
      if (itweak .ne. 0)  go to 3010
      if (  nimp .eq. 0)  go to 4000


      CALL get_charge_state(te)


c
      timzen = timnew   ! update for next zen call (timzen local to zen)
      if (inenez .eq. 0)  go to 4000
      if (inenez .eq. -99)  go to 4001 !JMP
      if (inenez .eq. -98)  go to 4001 !JMP
c
c          inenez # 0
c          compute ion densities from ene and zeff
c
 3010  ld_en = SIZE(en,1)
       set_u = .TRUE.   ! define the u array here
       CALL en_calc(ene,zeff,en,ld_en,izenstop,set_u)

       if (izenstop .eq. 0) return
       write  (nout, 9030)
 9030  format (1x, 'ONETWO run was terminated by subroutine ZEN',
     .            ' due to negative ion density.' /
     .         3x, 'You can try to run with ADJZEFF = 1' /
     .         3x, 'or else make adjustments in ZEFF, ENE, BPTOR, etc' /
     .         3x, 'to get around this problem.')
       call STOP ('subroutine ZEN: negative ion density', 285)
c
****  inenez = 0
c     compute ene and zeff from given ion densities
c

 4000 continue

c        IF user requested initial density input multiplier
c        This is no longer effective, multipliers are done in
c        sub mod_profiles,  HSJ 4/2909
         DO k=1,nion
           IF(ABS(density_mult(k) -1.d0) .GT. 1.e-13)THEN
            DO j = 1,nj
               IF (( fix_edge_ni(1) .GT.  1.0  .AND. j .lt.
     &             fix_edge_ni(1) ) .OR. (fix_edge_ni(1) .lt. 1.0 
     &          .AND. roa(j) .lt. fix_edge_ni(1)))THEN
                  en(j,k) = en(j,k)*density_mult(k)
               ENDIF
            ENDDO
           ENDIF
         ENDDO
        density_mult(:) = 1.d0  ! do it only onece


      do j=1,nj
        ene (j) = 0.0
        zeff(j) = 0.0
        do k=1,nion
          ene(j) = ene(j)+z(j,k)*en(j,k)
          zeff(j) = zeff(j)+zsq(j,k)*en(j,k)
        end do
        ene(j)  = ene(j) + enbeam(j) + 2.0*enalp(j)
        zeff(j) = zeff(j) + enbeam(j) + 4.0 * enalp(j)
        zeff(j) = zeff(j)/ene(j)
      end do

      !if ((itran(1) .eq. 1).and.(nprim .eq .2)) then !jmp.start
      !  do j=1,nj
      !    en(j,2) = en(j,1)
      !  end do 
      !end if !jmp.end

      goto 4002 !JMP
 
 4001 do j=1,nj !JMP BLOCK START

        en(j,1) = ene(j)
        do k=nprim+1,nion
          en(j,1) = en(j,1) - z(j,k)*en(j,k)
        end do
        en(j,1) = en(j,1) - enbeam(j) - 2.0*enalp(j)

        if (nprim .eq. 2) then
          en(j,2) = (1.0 - zfrac)*en(j,1)
          en(j,1) = zfrac*en(j,1)
        end if  

      end do

      if (inenez .eq. -98) then

         do k=1,nimp
         	
         	if(namei(k).eq.'he') then
         		do j=1,nj
         		  en(j,nprim+k) = taupin * dtnfus(j)
         		end do
         	end if		  
         
         end do
      
      end if

      do j=1,nj
        zeff(j) = 0.0
        do k=1,nion
          zeff(j) = zeff(j)+zsq(j,k)*en(j,k)
        end do
        zeff(j) = zeff(j) + enbeam(j) + 4.0 * enalp(j)
        zeff(j) = zeff(j)/ene(j)
      end do

 4002 continue !JMP BLOCK END

c
      if (zenvb .gt. 0) then
        write (*, '(/ " ON EXIT FROM ZEN WE HAVE nprim = ", i5)') nprim
        write (*, '(  " ENE AND ZEFF DETERMINED FROM GIVEN",
     .               " ION DENSITIES")')
        do i=1,nprim
          write (*, '(" i = ", i3,
     .                " en(1:5,i) = ",5(1pe14.4))')i,(en(j,i),j=1,5)
        end do
        write (*, '("on exit from ZEN we have nimp = ", i5)') nimp
        if (nimp .gt. 0) then
          do i=1,nimp
            write (*,'(" i = ",i3,
     .                 " en(1-5,nprim+i) = ",5(1pe14.4))')
     .                   i, (en(j,nprim+i),j=1,5)
          end do
        end if
        write (*, '(" enalp (1:5) =  ", 5(1pe14.4))') (enalp (j),j=1,5)
        write (*, '(" enbeam(1:5) =  ", 5(1pe14.4))') (enbeam(j),j=1,5)
        write (*, '(" ene   (1:5) =  ", 5(1pe14.4))') (ene   (j),j=1,5)
        do j=1,5
          sumd(j) = enbeam(j)
          do i=1,nion
            sumd(j) = sumd(j)+en(j,i)
          end do
        end do
        write (*, '(" total ion density at r(1:5) =" /
     .                5(1pe14.4))') (sumd(j),j=1,5)
        write (*, '(///)')
        do j=1,nj
          if(ene(j) .le.0.0)then
             print *,'in zen,ene ,j =',ene(j),j
             do i=1,nion
               print *,'in zen,en(j,i) =',en(j,i)
             enddo
          endif
        enddo
      end if
      return
c
      end

      subroutine zeroa (x, n)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      dimension x(n)
c
      do i=1,n
        x(i) = 0.0
      end do
      return
c
      end

      real*8 function zeroin (ax, bx, f, tol)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      external F
c
      real*8  ax, bx, f, tol
c
c a zero of the function  f(x)  is computed in the interval ax,bx
c
c  input..
c
c  ax     left  endpoint of initial interval
c  bx     right endpoint of initial interval
c  f      function subprogram which evaluates f(x) for any x in
c         the interval  ax,bx
c  tol    desired length of the interval of uncertainty of the
c         final result ( .ge. 0.0)
c
c  output..
c
c  zeroin abscissa approximating a zero of  f  in the interval ax,bx
c
c      it is assumed  that   f(ax)   and   f(bx)   have  opposite  signs
c  without  a  check.  zeroin  returns a zero  x  in the given interval
c  ax,bx  to within a tolerance  4*macheps * ABS (x) + tol, where macheps
c  is the relative machine precision.
c  this function subprogram is a slightly-modified translation of
c  the ALGOL 60 (JESUS!) procedure  ZERO  given in  Richard Brent,
c  Algorithms for Minimization without Derivatives, Prentice-Hall, Inc. (1973).
c
      real*8  a, b, c, d, e, eps, fa, fb, fc, tol1, xm, p, q, r, s
c
c     compute eps, the relative machine precision
c
      eps  = 1.0
   10 eps  = eps / 2.0
      tol1 = 1.0 + eps
      if (tol1 .gt. 1.0)  go to 10
c
c     initialization
c
      a  = ax
      b  = bx
      fa = f(a)
      fb = f(b)
c
c     begin step
c
   20 c  = a
      fc = fa
      d  = b - a
      e  = d
   30 if (ABS (fc) .ge. ABS (fb))  go to 40
      a  = b
      b  = c
      c  = a
      fa = fb
      fb = fc
      fc = fa
c
c     convergence test
c
   40 tol1 = 2.0 * eps * ABS (b) + 0.5*tol
      xm   = 0.5*(c - b)
      if (ABS (xm) .le. tol1)  go to 90
      if (fb .eq. 0.0)  go to 90
c
c     is bisection necessary
c
      if (ABS (e) .lt. tol1)  go to 70
      if (ABS (fa) .le. ABS (fb))  go to 70
c
c     is quadratic interpolation possible
c
      if (a .ne. c)  go to 50
c
c     linear interpolation
c
      s = fb/fa
      p = 2.0 * xm*s
      q = 1.0 - s
      go to 60
c
c     inverse quadratic interpolation
c
   50 q = fa/fc
      r = fb/fc
      s = fb/fa
      p = s*(2.0*xm*q*(q - r) - (b - a)*(r - 1.0))
      q = (q - 1.0)*(r - 1.0)*(s - 1.0)
c
c     adjust signs
c
   60 if (p .gt. 0.0)  q = -q
      p = ABS (p)
c
c     is interpolation acceptable
c
      if ((2.0*p) .ge. (3.0*xm*q - ABS (tol1*q)))  go to 70
      if (p .ge. ABS (0.5*e*q))  go to 70
      e = d
      d = p/q
      go to 80
c
c     bisection
c
   70 d = xm
      e = d
c
c     complete step
c
   80 a  = b
      fa = fb
      if (ABS (d) .gt. tol1)  b = b + d
      if (ABS (d) .le. tol1)  b = b + SIGN (tol1, xm)
      fb = f(b)
      if ((fb*(fc / ABS (fc))) .gt. 0.0)  go to 20
      go to 30
c
c     done
c
   90 zeroin = b
      return
c
      end

      subroutine zero_intg (ix, n)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      dimension ix(n)
c
      do i=1,n
        ix(i) = 0.0
      end do
      return
c
      end

      subroutine zfit (zav, te, namei, nj, ncrt, nout,
     .                  te_range_check)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c this subroutine computes the average charge state for an impurity.
c coeffficents for the curve fits used were obtained from a princeton
c paper, steady state radiative cooling rates for low-density high-
c temperture plasmas (pppl-1352)
c
c name(i)    = chemical abbreviation for element
c coeff(i,1) = lower temperature range
c coeff(i,2) = upper temperature range
c coeff(i,3) ---> coeff(i,8) are the coefficients of a 5th degree
c
      parameter   (krows = 44)
      character*8  name, namei
      dimension    name(krows), coeff(krows,8)
      dimension    zav(*), te(*)
      logical te_range_check
c
c ----------------------------------------------------------------------
c
c *********  argon             **********
      data name( 1) /'ar'/, (coeff( 1,i),i = 1,8)/
     .3.000e-02,2.000e-01,-6.35188e+01,-4.14519e+02,
     .-8.50220e+02,-8.07478e+02,-3.62126e+02,-6.18783e+01/
      data name( 2) /'ar'/, (coeff( 2,i),i = 1,8)/
     .2.000e-01,2.000e+00, 1.59107e+01,-7.88677e-01,
     . 2.87454e+00, 3.36119e+01,-3.30689e+01,-7.16260e+01/
      data name( 3) /'ar'/, (coeff( 3,i),i = 1,8)/
     .2.000e+00,2.000e+01, 1.29638e+01, 1.83325e+01,
     .-2.83480e+01, 2.26716e+01,-9.21974e+00, 1.50713e+00/
      data name( 4) /'ar'/, (coeff( 4,i),i = 1,8)/
     .2.000e+01,1.000e+02,-7.89001e+01, 2.93992e+02,
     .-3.55108e+02, 2.13368e+02,-6.37657e+01, 7.58257e+00/
c *********  carbon            **********
      data name( 5) /'c '/, (coeff( 5,i),i = 1,8)/
     .2.000e-03,1.370e-02,-2.15844e+03,-4.74739e+03,
     .-4.12560e+03,-1.77619e+03,-3.79399e+02,-3.21921e+01/
      data name( 6) /'c '/, (coeff( 6,i),i = 1,8)/
     .1.370e-02,4.047e-02, 4.00000e+00, 0.         ,
     . 0.         , 0.         , 0.         , 0.         /
      data name( 7) /'c '/, (coeff( 7,i),i = 1,8)/
     .4.047e-02,1.900e-01, 9.30638e+00, 4.16782e+01,
     . 1.27151e+02, 1.57530e+02, 8.45156e+01, 1.64735e+01/
      data name( 8) /'c '/, (coeff( 8,i),i = 1,8)/
     .1.900e-01,1.000e+02, 6.00000e+00, 0.         ,
     . 0.         , 0.         , 0.         , 0.         /
c *********  chromium          **********
      data name( 9) /'cr'/, (coeff( 9,i),i = 1,8)/
     .2.000e-02,2.000e-01,-7.69779e+01,-4.03142e+02,
     .-6.77843e+02,-5.42278e+02,-2.12927e+02,-3.30324e+01/
      data name(10) /'cr'/, (coeff(10,i),i = 1,8)/
     .2.000e-01,2.000e+00, 1.99540e+01, 1.57264e+01,
     .-2.60513e+01,-5.67938e+01, 8.44905e+01, 1.33601e+02/
      data name(11) /'cr'/, (coeff(11,i),i = 1,8)/
     .2.000e+00,2.000e+01, 1.84837e+01, 3.11439e+01,
     .-1.07363e+02, 1.69120e+02,-1.18178e+02, 3.02290e+01/
      data name(12) /'cr'/, (coeff(12,i),i = 1,8)/
     .2.000e+01,1.000e+02,-3.75732e+01, 1.85446e+02,
     .-2.25017e+02, 1.36696e+02,-4.14507e+01, 5.01119e+00/
c *********  iron              **********
      data name(13) /'fe'/, (coeff(13,i),i = 1,8)/
     .2.000e-03,2.000e-02,-1.07943e+03,-2.54883e+03,
     .-2.36318e+03,-1.08483e+03,-2.47337e+02,-2.24317e+01/
      data name(14) /'fe'/, (coeff(14,i),i = 1,8)/
     .2.000e-02,2.000e-01,-3.33248e+01,-2.64547e+02,
     .-5.05432e+02,-4.27802e+02,-1.69118e+02,-2.54619e+01/
      data name(15) /'fe'/, (coeff(15,i),i = 1,8)/
     .2.000e-01,2.000e+00, 1.90230e+01, 1.79155e+01,
     . 2.39328e+01,-5.69628e+01,-1.71077e+02,-1.06708e+02/
      data name(16) /'fe'/, (coeff(16,i),i = 1,8)/
     .2.000e+00,2.000e+01, 1.94569e+01, 2.87694e+01,
     .-7.75938e+01, 1.03931e+02,-6.41141e+01, 1.47605e+01/
      data name(17) /'fe'/, (coeff(17,i),i = 1,8)/
     .2.000e+01,1.000e+02, 2.12208e+01, 6.89161e+00,
     .-4.07685e+00, 1.57717e+00,-5.13947e-01, 8.93418e-02/
c *********  helium            **********
      data name(18) /'he'/, (coeff(18,i),i = 1,8)/
     .2.000e-03,1.000e-02, 2.23785e+03, 5.36958e+03,
     . 5.11718e+03, 2.41851e+03, 5.66901e+02, 5.27445e+01/
      data name(19) /'he'/, (coeff(19,i),i = 1,8)/
     .1.000e-02,1.000e+02, 2.00000e+00, 0.         ,
     . 0.         , 0.         , 0.         , 0.         /
c *********  krypton           **********
      data name(20) /'kr'/, (coeff(20,i),i = 1,8)/
     .5.000e-02,2.000e-01, 2.73740e+02, 1.42106e+03,
     . 3.10855e+03, 3.38963e+03, 1.82089e+03, 3.83951e+02/
      data name(21) /'kr'/, (coeff(21,i),i = 1,8)/
     .2.000e-01,2.000e+00, 2.47735e+01, 1.09293e+01,
     .-3.98248e+01,-5.75164e+00, 1.58536e+02, 1.51381e+02/
      data name(22) /'kr'/, (coeff(22,i),i = 1,8)/
     .2.000e+00,2.000e+01, 4.20480e+01,-1.66164e+02,
     . 5.70516e+02,-7.92204e+02, 4.95793e+02,-1.16257e+02/
      data name(23) /'kr'/, (coeff(23,i),i = 1,8)/
     .2.000e+01,1.000e+02, 6.97899e+01,-1.22590e+02,
     . 1.55217e+02,-9.33267e+01, 2.74276e+01,-3.18901e+00/
c *********  molybdenum        **********
      data name(24) /'mo'/, (coeff(24,i),i = 1,8)/
     .6.000e-02,2.000e-01, 4.99240e+01, 8.50640e+01,
     .-7.97318e+01,-3.87039e+02,-3.71884e+02,-1.12086e+02/
      data name(25) /'mo'/, (coeff(25,i),i = 1,8)/
     .2.000e-01,2.000e+00, 2.44902e+01, 2.68426e+01,
     . 1.29935e+01,-4.14538e+01,-1.49693e+02,-1.38544e+02/
      data name(26) /'mo'/, (coeff(26,i),i = 1,8)/
     .2.000e+00,2.000e+01, 3.41244e+01, 4.46808e+00,
     .-1.30193e+02, 3.62476e+02,-3.30031e+02, 9.86647e+01/
      data name(27) /'mo'/, (coeff(27,i),i = 1,8)/
     .2.000e+01,1.000e+02, 1.36713e+02,-2.28731e+02,
     . 2.00524e+02,-7.89953e+01, 1.28869e+01,-4.47343e-01/
c *********  nickel            **********
      data name(28) /'ni'/, (coeff(28,i),i = 1,8)/
     .3.000e-02,2.000e-01, 1.03241e+01,-6.18430e+01,
     .-1.60005e+02,-1.56411e+02,-7.15398e+01,-1.26777e+01/
      data name(29) /'ni'/, (coeff(29,i),i = 1,8)/
     .2.000e-01,2.000e+00, 1.92973e+01, 1.59013e+01,
     . 4.11513e+01,-2.44004e+01,-2.10927e+02,-1.65174e+02/
      data name(30) /'ni'/, (coeff(30,i),i = 1,8)/
     .2.000e+00,2.000e+01, 1.74633e+01, 4.98534e+01,
     .-1.12517e+02, 1.16930e+02,-5.36758e+01, 8.53433e+00/
      data name(31) /'ni'/, (coeff(31,i),i = 1,8)/
     .2.000e+01,1.000e+02,-8.43480e+01, 3.26564e+02,
     .-3.83686e+02, 2.26184e+02,-6.66565e+01, 7.84120e+00/
c *********  oxygen            **********
      data name(32) /'o '/, (coeff(32,i),i = 1,8)/
     .2.000e-03,2.000e-02,-6.49058e+02,-1.65379e+03,
     .-1.64411e+03,-8.04612e+02,-1.94499e+02,-1.86055e+01/
      data name(33) /'o '/, (coeff(33,i),i = 1,8)/
     .2.000e-02,2.000e-01,-1.13871e+01,-1.00970e+02,
     .-1.97775e+02,-1.72799e+02,-6.80264e+01,-9.44028e+00/
      data name(34) /'o '/, (coeff(34,i),i = 1,8)/
     .2.000e-01,6.000e-01, 7.98714e+00, 6.12345e-02,
     . 1.95765e-01,-3.53047e-02,-3.87677e+00, 1.24360e+00/
      data name(35) /'o '/, (coeff(35,i),i = 1,8)/
     .6.000e-01,1.000e+02, 8.00000e+00, 0.         ,
     . 0.         , 0.         , 0.         , 0.         /
c *********  silicon           **********
      data name(36) /'si'/, (coeff(36,i),i = 1,8)/
     .2.000e-03,2.000e-02, 3.80969e+03, 9.13991e+03,
     . 8.70264e+03, 4.10486e+03, 9.58975e+02, 8.87920e+01/
      data name(37) /'si'/, (coeff(37,i),i = 1,8)/
     .2.000e-02,2.000e-01, 1.26254e+01,-4.21502e+01,
     .-1.32778e+02,-1.39890e+02,-6.62081e+01,-1.19168e+01/
      data name(38) /'si'/, (coeff(38,i),i = 1,8)/
     .2.000e-01,2.000e+00, 1.27831e+01, 4.00160e+00,
     . 2.22816e+00,-1.29867e+01,-1.35309e+01, 6.70687e+00/
      data name(39) /'si'/, (coeff(39,i),i = 1,8)/
     .2.000e+00,4.000e+00, 1.28282e+01, 5.46582e+00,
     .-1.08022e+01, 1.09604e+01,-5.61033e+00, 1.14926e+00/
      data name(40) /'si'/, (coeff(40,i),i = 1,8)/
     .4.000e+00,1.000e+02, 1.40000e+01, 0.         ,
     . 0.         , 0.         , 0.         , 0.         /
c *********  tungsten          **********
      data name(41) /'w '/, (coeff(41,i),i = 1,8)/
     .1.000e-01,2.000e-01, 1.55578e+02, 9.07526e+02,
     . 2.34466e+03, 3.05061e+03, 1.99510e+03, 5.25087e+02/
      data name(42) /'w '/, (coeff(42,i),i = 1,8)/
     .2.000e-01,2.000e+00, 2.56069e+01, 2.74645e+01,
     . 1.88036e+01, 2.10978e+01,-1.77969e+01,-5.67164e+01/
      data name(43) /'w '/, (coeff(43,i),i = 1,8)/
     .2.000e+00,2.000e+01,-8.08150e+01, 7.80711e+02,
     .-1.79349e+03, 1.85649e+03,-8.28859e+02, 1.20794e+02/
      data name(44) /'w '/, (coeff(44,i),i = 1,8)/
     .2.000e+01,1.000e+02, 6.27044e+03,-1.82726e+04,
     . 2.13281e+04,-1.23551e+04, 3.55778e+03,-4.07816e+02/
c
      data nrows/44/
c
      do 6000 j=1,nj
      tkev   = te(j)
c     out of range allowed:
      if(.not. te_range_check )
     .        tkev   = Min(tkev,100.D0) !jmp.ibm
      zav(j) = 0.0
c
      do i=1,nrows
        if (name(i) .eq. namei)  go to 2010
      end do
      write  (ncrt, 8000) namei
      write  (nout, 8000) namei
 8000 format (' FATAL ERROR in ZFIT'                              /
     .        ' impurity ', a8, ' has not been added to the code' /)
      call STOP ('subroutine ZFIT: problem #1', 90)
c
 2010 if (tkev .ge. coeff(i,1))  go to 2050
      irow = i
      y    = coeff(irow,8)
      x    = LOG10 (coeff(irow,1))
      do i=8,4,-1
        y    = y * x + coeff(irow,i-1)
      end do
      zav(j) = tkev*(y/coeff(irow,1))
      go to 6000
c
 2050 i0 = i
      do i=i0,nrows
        if (tkev .ge. coeff(i,1) .and. tkev .le. coeff(i,2))  go to 2100
        if (                        name(i) .ne. namei     )  go to 9000
      end do
      go to 9000
c
 2100 irow = i
      y    = coeff(irow,8)
      x    = LOG10 (tkev)
      do i=8,4,-1
        y  = y*x+coeff(irow,i-1)
      end do
      zav(j) = y
c
 6000 continue
      return
c
c ----------------------------------------------------------------------
c fatal errors
c ----------------------------------------------------------------------
c
 9000 write  (nout, 8010)
      write  (ncrt, 8010)tkev
 8010 format (' FATAL ERROR in ZFIT'                                /
     .        ' you have exceeded temperature range for curve fits' /,
     .        '  tkev = ',1pe12.4)
      call STOP ('subroutine ZFIT: problem #2', 91)
c
      end

      subroutine zsqfit (zsqav, te, namei, nj, ncrt, nout,
     .                                         te_range_check)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c this subroutine computes the mean square charge state for an impurity.
c coeffficents for the curve fits used were obtained from a princeton
c paper, steady state radiative cooling rates for low-density high-
c temperture plasmas (pppl-1352)
c
c name(i)    = chemical abbreviation for element
c coeff(i,1) = lower temperature range
c coeff(i,2) = upper temperature range
c coeff(i,3) ---> coeff(i,8) are the coefficients of a 5th degree
c
      parameter       (krows = 46)
      character*8 name, namei
      logical     te_range_check
      dimension   name(krows), coeff(krows,8)
      dimension   zsqav(*), te(*)
c
c ----------------------------------------------------------------------
c
c *********  argon             **********
c
      data name( 1) /'ar'/, (coeff( 1,i),i = 1,8)/
     .3.000e-02,2.000e-01,-9.61710e+02,-6.47898e+03,
     .-1.39467e+04,-1.37241e+04,-6.35091e+03,-1.12035e+03/
      data name( 2) /'ar'/, (coeff( 2,i),i = 1,8)/
     .2.000e-01,2.000e+00, 2.53414e+02,-2.34631e+01,
     . 7.07563e+01, 1.01411e+03,-7.34252e+02,-1.96330e+03/
      data name( 3) /'ar'/, (coeff( 3,i),i = 1,8)/
     .2.000e+00,2.000e+01, 1.58253e+02, 5.76331e+02,
     .-8.29530e+02, 5.95453e+02,-2.05795e+02, 2.60289e+01/
      data name( 4) /'ar'/, (coeff( 4,i),i = 1,8)/
     .2.000e+01,1.000e+02, 2.95994e+02, 6.59094e+01,
     .-7.35873e+01, 4.64761e+01,-1.56755e+01, 2.16710e+00/
c
c *********  carbon            **********
c
      data name( 5) /'c '/, (coeff( 5,i),i = 1,8)/
     .2.000e-03,1.300e-02,-1.19723e+04,-2.59832e+04,
     .-2.22347e+04,-9.39972e+03,-1.96582e+03,-1.62847e+02/
      data name( 6) /'c '/, (coeff( 6,i),i = 1,8)/
     .1.300e-02,4.000e-02, 1.60000e+01, 0.         ,
     . 0.         , 0.         , 0.         , 0.         /
      data name( 7) /'c '/, (coeff( 7,i),i = 1,8)/
     .4.000e-02,1.900e-01,-5.23254e+01,-1.42389e+02,
     . 2.75141e+02, 7.17200e+02, 4.86806e+02, 1.06456e+02/
      data name( 8) /'c '/, (coeff( 8,i),i = 1,8)/
     .1.900e-01,1.000e+02, 3.60000e+01, 0.         ,
     . 0.         , 0.         , 0.         , 0.         /
c
c *********  chromium          **********
c
      data name( 9) /'cr'/, (coeff( 9,i),i = 1,8)/
     .2.000e-02,2.000e-01,-2.33301e+03,-1.13137e+04,
     .-1.92235e+04,-1.55134e+04,-6.09525e+03,-9.41303e+02/
      data name(10) /'cr'/, (coeff(10,i),i = 1,8)/
     .2.000e-01,2.000e+00, 4.00504e+02, 5.98513e+02,
     .-7.84260e+02,-2.27763e+03, 2.08980e+03, 4.15338e+03/
      data name(11) /'cr'/, (coeff(11,i),i = 1,8)/
     .2.000e+00,2.000e+01, 3.29830e+02, 1.38295e+03,
     .-4.80366e+03, 7.59438e+03,-5.31360e+03, 1.35960e+03/
      data name(12) /'cr'/, (coeff(12,i),i = 1,8)/
     .2.000e+01,1.000e+02, 1.04511e+03,-1.72775e+03,
     . 2.31388e+03,-1.47903e+03, 4.59971e+02,-5.61846e+01/
c
c *********  iron              **********
c
      data name(13) /'fe'/, (coeff(13,i),i = 1,8)/
     .2.000e-03,2.000e-02,-1.40670e+04,-3.42381e+04,
     .-3.26786e+04,-1.53740e+04,-3.57550e+03,-3.29390e+02/
      data name(14) /'fe'/, (coeff(14,i),i = 1,8)/
     .2.000e-02,2.000e-01,-3.47441e+03,-1.76130e+04,
     .-3.16595e+04,-2.68416e+04,-1.09495e+04,-1.73609e+03/
      data name(15) /'fe'/, (coeff(15,i),i = 1,8)/
     .2.000e-01,2.000e+00, 3.65613e+02, 7.05549e+02,
     . 1.10748e+03,-2.11670e+03,-7.37044e+03,-5.03081e+03/
      data name(16) /'fe'/, (coeff(16,i),i = 1,8)/
     .2.000e+00,2.000e+01, 3.67376e+02, 1.32217e+03,
     .-3.57325e+03, 4.79087e+03,-2.94446e+03, 6.73278e+02/
      data name(17) /'fe'/, (coeff(17,i),i = 1,8)/
     .2.000e+01,1.000e+02,-3.24362e+02, 2.47613e+03,
     .-2.55182e+03, 1.34865e+03,-3.61953e+02, 3.92133e+01/
c
c *********  helium            **********
c
      data name(18) /'he'/, (coeff(18,i),i = 1,8)/
     .2.000e-03,2.000e-02, 6.62735e+03, 1.58070e+04,
     . 1.49612e+04, 7.01809e+03, 1.63142e+03, 1.50392e+02/
      data name(19) /'he'/, (coeff(19,i),i = 1,8)/
     .1.200e-02,1.000e+02, 4.00000e+00, 0.         ,
     . 0.         , 0.         , 0.         , 0.         /
      data name(20) /'he'/, (coeff(20,i),i = 1,8)/
     .2.000e-02,1.200e-02, 4.04302e+00, 2.10593e-01,
     . 4.07570e-01, 3.90840e-01, 1.86209e-01, 3.53502e-02/
c
c *********  krypton           **********
c
      data name(21) /'kr'/, (coeff(21,i),i = 1,8)/
     .5.000e-02,2.000e-01, 5.38182e+03, 2.74893e+04,
     . 5.85402e+04, 6.26166e+04, 3.31931e+04, 6.93547e+03/
      data name(22) /'kr'/, (coeff(22,i),i = 1,8)/
     .2.000e-01,2.000e+00, 6.14316e+02, 5.14605e+02,
     .-1.71603e+03,-6.50115e+02, 7.27807e+03, 7.34126e+03/
      data name(23) /'kr'/, (coeff(23,i),i = 1,8)/
     .2.000e+00,2.000e+01, 1.77994e+03,-1.09043e+04,
     . 3.63663e+04,-4.97846e+04, 3.08695e+04,-7.18792e+03/
      data name(24) /'kr'/, (coeff(24,i),i = 1,8)/
     .2.000e+01,1.000e+02, 6.40233e+02, 8.44951e+02,
     .-7.68434e+02, 5.66236e+02,-2.22623e+02, 3.27402e+01/
c
c *********  molybdenum        **********
c
      data name(25) /'mo'/, (coeff(25,i),i = 1,8)/
     .6.000e-02,2.000e-01, 2.01456e+03, 7.29892e+03,
     . 9.50696e+03, 2.94308e+03,-2.49193e+03,-1.36698e+03/
      data name(26) /'mo'/, (coeff(26,i),i = 1,8)/
     .2.000e-01,2.000e+00, 6.06171e+02, 1.34140e+03,
     . 1.23873e+03,-2.32759e+03,-9.06615e+03,-7.68058e+03/
      data name(27) /'mo'/, (coeff(27,i),i = 1,8)/
     .2.000e+00,2.000e+01, 1.03905e+03, 1.54394e+03,
     .-1.29818e+04, 3.07817e+04,-2.65536e+04, 7.73080e+03/
      data name(28) /'mo'/, (coeff(28,i),i = 1,8)/
     .2.000e+01,1.000e+02,-1.04401e+04, 4.14475e+04,
     .-5.56599e+04, 3.64213e+04,-1.16273e+04, 1.45479e+03/
c
c *********  nickel            **********
c
      data name(29) /'ni'/, (coeff(29,i),i = 1,8)/
     .3.000e-02,2.000e-01, 5.68427e+02, 3.33388e+02,
     .-1.01141e+03,-1.50432e+03,-7.74428e+02,-1.41858e+02/
      data name(30) /'ni'/, (coeff(30,i),i = 1,8)/
     .2.000e-01,2.000e+00, 3.75180e+02, 6.53741e+02,
     . 1.82044e+03,-7.30319e+02,-8.66849e+03,-7.22875e+03/
      data name(31) /'ni'/, (coeff(31,i),i = 1,8)/
     .2.000e+00,2.000e+01, 2.45680e+02, 2.51595e+03,
     .-5.65959e+03, 5.82043e+03,-2.60969e+03, 3.94659e+02/
      data name(32) /'ni'/, (coeff(32,i),i = 1,8)/
     .2.000e+01,1.000e+02,-5.46560e+02, 3.21744e+03,
     .-3.28362e+03, 1.73827e+03,-4.71450e+02, 5.19611e+01/
c
c *********  oxygen            **********
c
      data name(33) /'o '/, (coeff(33,i),i = 1,8)/
     .2.000e-03,2.000e-02,-1.17046e+03,-3.48550e+03,
     .-3.80639e+03,-1.98106e+03,-4.99595e+02,-4.92484e+01/
      data name(34) /'o '/, (coeff(34,i),i = 1,8)/
     .2.000e-02,2.000e-01,-1.31005e+02,-1.10200e+03,
     .-2.31452e+03,-2.14423e+03,-9.02654e+02,-1.38299e+02/
      data name(35) /'o '/, (coeff(35,i),i = 1,8)/
     .2.000e-01,2.000e+00, 6.38084e+01, 8.67344e-01,
     . 3.05377e+00, 4.64915e-01,-5.95321e+01, 1.18394e+01/
      data name(36) /'o '/, (coeff(36,i),i = 1,8)/
     .2.000e+00,2.000e+01, 6.37958e+01, 1.18364e+00,
     .-2.98014e+00, 3.80390e+00,-2.37768e+00, 5.75365e-01/
      data name(37) /'o '/, (coeff(37,i),i = 1,8)/
     .2.000e+01,1.000e+02,-3.63953e+02, 1.33378e+03,
     .-1.65417e+03, 1.02059e+03,-3.13304e+02, 3.82889e+01/
c
c *********  silicon           **********
c
      data name(38) /'si'/, (coeff(38,i),i = 1,8)/
     .2.000e-03,2.000e-02, 2.19672e+04, 5.27672e+04,
     . 5.02832e+04, 2.37329e+04, 5.54687e+03, 5.13678e+02/
      data name(39) /'si'/, (coeff(39,i),i = 1,8)/
     .2.000e-02,2.000e-01, 5.40379e+02, 7.46322e+02,
     .-1.21937e+02,-8.04741e+02,-5.27067e+02,-1.08871e+02/
      data name(40) /'si'/, (coeff(40,i),i = 1,8)/
     .2.000e-01,2.000e+00, 1.63920e+02, 1.03378e+02,
     . 6.78394e+01,-3.28525e+02,-4.22436e+02, 6.23053e+01/
      data name(41) /'si'/, (coeff(41,i),i = 1,8)/
     .2.000e+00,2.000e+01, 1.63752e+02, 1.50787e+02,
     .-2.97061e+02, 2.99360e+02,-1.51746e+02, 3.06711e+01/
      data name(42) /'si'/, (coeff(42,i),i = 1,8)/
     .2.000e+01,1.000e+02, 3.62084e+02,-5.79022e+02,
     . 7.89878e+02,-5.29103e+02, 1.74421e+02,-2.26710e+01/
c
c *********  tungsten          **********
c
      data name(43) /'w '/, (coeff(43,i),i = 1,8)/
     .1.000e-01,2.000e-01,-1.05752e+04,-6.46137e+04,
     .-1.53578e+05,-1.80606e+05,-1.05008e+05,-2.40900e+04/
      data name(44) /'w '/, (coeff(44,i),i = 1,8)/
     .2.000e-01,2.000e+00, 6.62079e+02, 1.41753e+03,
     . 1.72297e+03, 1.88480e+03, 3.56943e+02,-1.38435e+03/
      data name(45) /'w '/, (coeff(45,i),i = 1,8)/
     .2.000e+00,2.000e+01,-7.48175e+03, 5.66668e+04,
     .-1.21066e+05, 1.08444e+05,-3.34462e+04,-8.94193e+01/
      data name(46) /'w '/, (coeff(46,i),i = 1,8)/
     .2.000e+01,1.000e+02, 8.06427e+05,-2.35879e+06,
     . 2.74916e+06,-1.59020e+06, 4.57256e+05,-5.23409e+04/
c
      data nrows /46/
c
      do 6000 j=1,nj
        tkev     = te(j)
c     out of range allowed:
      if(.not. te_range_check )
     .        tkev   = Min(tkev,100.D0) !jmp.ibm
        zsqav(j) = 0.0
c
        do i=1,nrows
          if (name(i) .eq. namei)  go to 2010
        end do
        write (ncrt, 8000) namei
        write (nout, 8000) namei
 8000   format (' FATAL ERROR in ZSQFIT'                            /
     .          ' impurity ', a8, ' has not been added to the code' /)
        call STOP ('subroutine ZSQFIT: problem #1', 92)
c
 2010   if (tkev .ge. coeff(i,1))  go to 2050
        irow = i
        y    = coeff(irow,8)
        x    = LOG10 (coeff(irow,1))
        do i=8,4,-1
          y = y*x+coeff(irow,i-1)
        end do
        zsqav(j) = tkev*(y/coeff(irow,1))
        go to 6000
c
 2050   i0 = i
        do i=i0,nrows
          if (tkev .ge. coeff(i,1) .and.
     .        tkev .le. coeff(i,2))  go to 2100
          if (name(i) .ne. namei)  go to 9000
        end do
        go to 9000
c
 2100   irow = i
        y = coeff(irow,8)
        x = LOG10 (tkev)
        do i=8,4,-1
          y = y*x + coeff(irow,i-1)
        end do
        zsqav(j) = y
 6000 continue
c
      return
c
c ----------------------------------------------------------------------
c fatal errors
c ----------------------------------------------------------------------
c
 9000 write  (nout, 8010)
      write  (ncrt, 8010)
 8010 format (/ ' FATAL ERROR in ZSQFIT:'                             /
     .          ' you have exceeded temperature range for curve fits' /)
      call STOP ('subroutine ZSQFIT: problem #2', 93)
c
      end
