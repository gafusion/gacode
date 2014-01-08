c
      real*8 function alphac (shearp)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c --- calculate critical value of alpha for given shear parameter
c --- for stability to ideal ballooning - Taylor TPN/83/19
c
      dimension  alp(9), shear(9)
      data       alp/0.0,0.1,0.2,0.3,0.35,0.4,0.5,0.6,0.7/
      data       shear/0.0,0.03,0.06,0.15,0.42,0.58,0.81,1.01,1.21/
c
      do i=1,9
        if (shearp .lt. shear(i))  go to 20
      end do
c
c     use straight line
c
      alphac = 0.58*shearp
      return
c
   20 if (i .eq. 1)  go to 30
c
c     interpolate
c
      alphac = (alp(i)*(shearp-shear(i-1))+alp(i-1)*(shear(i)-shearp))/
     .         (shear(i)-shear(i-1))
      return
c
c     set alpha to zero for negative shear
c
   30 alphac = 0.0
      return
c
      end

      subroutine check_zero_bc
c
c ----------------------------------------------------------------------
c --- check to make sure boundary conditions are set properly
c ------------------------------------------------------------------ HSJ
c
      USE param
      USE io 
      USE ions
      USE solcon 
      USE numbrs
      USE soln, only :usave,u,u_continue
      USE flags
      USE bd_condtn,only : bctime,ub,fluxb,
     .    ub_save,ub_rho_edge,bctime_zone
      implicit  integer (i-n), real*8 (a-h, o-z)
c

c
      ierr = 0
      do j=1,nk
        if (ub(j) .eq. 0.0 .and. itran(j) .ne. 0) then
          write (ncrt,'("ERROR: boundary condition=0 for ",
     .                   " dependent variable: ",
     .                a," at time = ",1pe14.8)')nameu(j),time
          write (nout,'("ERROR: boundary condition=0 for ",
     .                   " dependent variable: ",
     .                a," at time = ",1pe14.8)')nameu(j),time
          write (nqik,'("ERROR: boundary condition=0 for ",
     .                   " dependent variable: ",
     .                a," at time = ",1pe14.8)')nameu(j),time
          ierr = 1
        end if

      end do
c
      if (ierr .gt. 0)
     .  call STOP ('subroutine CHECK_ZERO_BC: bad boundary cond.', 227)
      return
c
      end

      subroutine delsol (imax1, r1, q, modben, wben, wdtben, rsben,
     .                   btil, xlg, wdisl0, dpisl0)
c
      USE param
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c      include 'param.i'
      include 'delben.i'
c
      dimension  bcl(2),ro(2),dr(2),r(500,2),v(500,2),
     .           a(500,2),b(500,2),psi(500,2),
     .           e(500,2),f(500,2),delta(500),eta(500)
      dimension  r1(kj), q(kj), modben(2), xjben(kj)
      external   fqqq
c
c     set up interface
c
      nben = imax1-1
c
      do i=1,nben
        qben(i) = q(i+1)
        rben(i) = r1(i+1)/r1(imax1)
      end do
c
c     xjben(i) is diagnostic
c
      nbenm1 = nben - 1
c
      do i=1,nbenm1
        xjben(i) = qben(1)/(rben(i+1)+rben(i))*(rben(i+1)**2
     .            /qben(i+1)-rben(i)**2/qben(i))/(rben(i+1)-rben(i))
      end do
c
c     m,n mode to be calculated
c
c     calculation of the tearing mode flux function
c     outside the tearing layer and delta prime
c
c     parameters
c
      rng   = nben - 1
      drben = 1.0/rng
      alarm = 1.0e-20
      bcl(1) = 0.0
      bcl(2) = 1.0
      t0l   = 0.001
      mto   = 500-2
c
c     determination of rs
c
      rmode = modben(1)
      rnode = modben(2)
      qs    = rmode/rnode
      do 20 j=1,nben
      jkeep = j
   20 if (qben(j) .gt. qs)  go to 21
   21 if (jkeep .eq. 1 .or. jkeep .eq. nben)  go to 2000
      jk    = jkeep
      jkm   = jk - 1
      call spline_12 (nben,rben,qben,bq,cq,dq) 
      d1    = qben(jk)-qs
      d2    = qs-qben(jkm)
      ax    = rben(jkm)
      bx    = rben(jk)
      rs    = zeroin (ax, bx, fqqq, t0l)
      ro(1) = 0.0
      ro(2) = rs
      rnor  = qs/qben(jk)
      rsben = rs*r1(imax1)
c
c     setting up the grid
c
      mji    = rs*mto
      rmji   = mji
      dr(1)  = rs/rmji
      mje    = (1.0-rs)/dr(1)
      rmje   = mje
      dr(2)  = (1.0-rs)/rmje
      ndelta = mji
      if (mje .le. mji)  ndelta = mje
      mj     = mji
      if (mje .gt. mji)  mj = mje
c
      do i=1,2
        if (i .eq. 1)  mm = mji
        if (i .eq. 2)  mm = mje
        r(1,i) = ro(i)+dr(i)
        do j=2,mm
          jm     = j-1
          r(j,i) = r(jm,i)+dr(i)
        end do
      end do
c
      mjim1 = mji-1
      mjem1 = mje-1
      mjm1  = mj-1
c
c     calculation of the equilibrium
c
      do 2 i=1,2
      if (i .eq. 1) mm = mjim1
      if (i .eq. 2) mm = mjem1
      do 2 j=1,mm
      rr = r(j,i)
      jp = j+1
      jm = j-1
      rp = r(jp,i)
      if (j .eq. 1)  go to 51
      rm = r(jm,i)
      go to 52
   51 rm  = ro(i)
   52 dsq = (rr+rp)*rp/fqq(rp)+(rr+rm)*rm/fqq(rm)-4.0 * rr*rr/fqq(rr)
      dsq = 0.5 * dsq*fqq(rr)/(dr(i)*dr(i)*rr)-1.0/rr
      d   = rnode*fqq(rr)/rmode-1.0
      amd = d * d
      if (amd .le. alarm)  go to 165
      v(j,i) = dsq / d
      go to 166
c
  165 v(j,i) = 0.0
      n7 = 0 !error output to ftn0
      write  (n7, 101)
  101 format (10x, 'alarm')
      call STOP ('subroutine DELSOL: unspecified problem', 40)
  166 continue
    2 continue
c
c     calculation of the coefficients a, b, c
c
      do i=1,2
        if (i .eq. 1)  mm = mjim1
        if (i .eq. 2)  mm = mjem1
        do j=1,mm
          rr     = r(j,i)
          dd     =   2.0*rr+dr(i)
          a1     = dr(i)*dr(i)*(rmode*rmode/rr-v(j,i))+2.0*rr
          a(j,i) =   2.0*a1/dd
          b(j,i) = -(2.0*rr-dr(i))/dd
        end do
      end do
c
c     calculation of the coefficients e and f
c
      do i=1,2
        e(1,i) = 1.0/a(1,i)
        f(1,i) = -(b(1,i)*bcl(i))/a(1,i)
        if (i .eq. 1)  mm = mjim1
        if (i .eq. 2)  mm = mjem1
        do j=2,mm
          jm     = j-1
          dd     = a(j,i)+b(j,i)*e(jm,i)
          e(j,i) = 1.0/dd
          f(j,i) = -(b(j,i)*f(jm,i))/dd
        end do
      end do
c
c     calculation of the flux function psi
c
      psi(mji,1) = bcl(2)
      psi(mje,2) = 0.0
      if (xlg .eq. 1) psi(mje,2) = 1.0
      xdelt      = 1.0
      do 200 i=1,2
  201 if (i .eq. 1)  mm = mjim1
      if (i .eq. 2)  mm = mjem1
      do j=1,mm
        jj        = mm-j+1
        jjp       = jj+1
        psi(jj,i) = e(jj,i)*psi(jjp,i)+f(jj,i)
      end do
      if (xlg .eq. 0.0 .or. i .eq. 1)  go to 200
      ytest = -(psi(mje,2)-psi(mjem1,2))/(psi(mje,2)+psi(mjem1,2))*2.0/
     .          dr(2)-rmode
      xdelt = xdelt/2.0
      if (ytest .ge. 0.0)  psi(mje,2) = psi(mje,2)+xdelt
      if (ytest .lt. 0.0)  psi(mje,2) = psi(mje,2)-xdelt
      abyt  = ABS (ytest)
      ebyt  = 0.01*rmode
      if (abyt .lt. ebyt)  go to 200
      go to 201
  200 continue
c
c     calculation of delta prime
c
      do 15 l=1,ndelta
      ji   = mji-ndelta-1+l
      jip  = ji+1
      jim  = ji-1
      je   = ndelta+1-l
      jep  = je+1
      jem  = je-1
      if (l .eq. ndelta)  go to 80
      psii = (psi(jip,1)+psi(jim,1))/2.0
      psii = 1.0
      dd   = (psi(jip,1)-psi(jim,1))/dr(1)/psii
      psie = (psi(jep,2)+psi(jem,2))/2.
      psie = 1.0
      delta(l) = dd-(psi(jep,2)-psi(jem,2))/dr(2)/psie
      go to 81
   80 delta(l) = (1.0 - psi(jim,1))/dr(1)-(psi(jep,2)-1.0)/dr(2)
   81 delta(l) = -delta(l)*0.5
      eta(l)   = r(je,2)-r(ji,1)
      lm1      = l - 1
      ccs      = delta(l)*delta(lm1)
      if (ccs .le. 0.0) lkeep = l
   15 continue
c
      l      = lkeep
      lm     = l-1
      width  = delta(l)*eta(lm)-delta(lm)*eta(l)
      width  = width/(delta(l)-delta(lm))
      nn     = ndelta
      nnm1   = nn-1
      deltap = 2.0 * delta(nn)-delta(nnm1)
      if (deltap .le. 0.0) width = 0.0
      wben   = width*r1(imax1)
      wdtben = deltap
c
c     calculation of btilda/b
c
      dd     = 0.5 * (dr(1)+dr(2))
      rpdd   = rs+dd
      rmdd   = rs-dd
      qprime = 0.5 * (fqq(rpdd)-fqq(rmdd))/dd
      dpsi   = -psi(mjem1,2)/dr(2)
      btil   = -fqq(1.0)*(width**2)*rs*qprime*dpsi/64.0
c
c ----------------------------------------------------------------------
c dynamic island calculation rew 10/29/85
c ----------------------------------------------------------------------
c
      rin = rs-wdisl0/2.0/r1(imax1)
      rex = rs+wdisl0/2.0/r1(imax1)
      if (rin .lt. dr(1)) rin = dr(1)
      if (rex .gt. 1.0) rex = 1.0
      jisli = rin/dr(1)
      jisle = (rex-rs)/dr(2)
      if (jisle .eq. mje) jisle = mje-1
      jislim1 = jisli-1
      jislep1 = jisle+1
      psijisle = 1.0
      if (jisle .ne. 0)  psijisle = psi(jisle,2)
      psie = (psi(jislep1,2)+psijisle)/2.
      psii = (psi(jisli,1)+psi(jislim1,1))/2.
      dpisl0 = (psi(jislep1,2)-psijisle)/dr(2)/psie-
     .         (psi(jisli,1)-psi(jislim1,1))/dr(1)/psii
      dpisl0 = dpisl0/r1(imax1)
      return
c
 2000 wben   = 0.0
      wdtben = 0.0
      btil   = 0.0
      return
c
      end


      subroutine set_boundary_condition (timel)
c
c ---------------------------------------------------------------------
c  INPUT
c
c  timel
c  bctime(kbctim)
c  nbctim          # values in bctime
c  bc(kbctim,kk)
c  nk              nprim+nimp+1+1+1+? where ? is 1 if toroidal rotation is
c                  included and the 3 1's are for te,ti,rbp respectively.
c
c
c  OUTPUT
c    ub(k)     k=1,2...nk the value of the profile at the point
c                         r=rho_edge
c    bctime_zone( , )
c ------------------------------------------------------- HSJ-7/30/97 --
c
      USE param
      USE numbrs
      USE mesh
      USE tordlrot
      USE tfact,only : t_elms,t_elme,itot_elm,n_elms,k_elms
      USE bd_condtn
      implicit  integer (i-n), real*8 (a-h, o-z)



      real *8 ,dimension(:),allocatable :: yprof1,yprof2
c
      allocate (yprof1(1:nj),STAT = istat)
      if(istat .ne. 0)
     .          call allocate_error("yprof2 - set_bc",0,istat) 
      allocate (yprof2(1:nj),STAT = istat)
      if(istat .ne. 0)
     .          call allocate_error("yprof1 - set_bc",0,istat) 



      do k=1,nk
        ub(k) = bc(1,k)
        if (bc(2,k) .ne. 0.0 .and. k .ne. nk - iangrot)then
           call interp1 (timel, bctime, nbctim, bc(1,k), ub(k)) ! returns..
c                                                   ..ub(k) at t = timel
        
        else if (bc(2,k) .ne. 0.0 .and. k .eq. nk -iangrot )then
           if(.not. u_vloop_bc)then   !bc is total current 
              call interp1 (timel, bctime, nbctim, bc(1,k), ub(k))
           else !vloop bc 
                !this boundary condition is non linear and
                !changes as the profiles change durint the iteration
              ub(k) = 0.2*totcur(1)  !obtained from bd_condtn

           endif
        endif

c       modify total current to simulate elms if no of active elms > 0
        if(k_elms .gt. 0 .and. k .eq. nk-iangrot )then
          do j=1,k_elms
            if( t_elms(j) .le. timel .and.
     .       timel .le. t_elme(j))ub(k) = itot_elm(j)*ub(k)
          enddo
        endif


        if(use_pedestal)call Onetwo_pedestal_driver


        !determine values in the  edge zone at this time by linear 
        !interpolation from *_edge_zone
        !note that *_edge_zone(i,j) arrays are defined (in sub bctime_zone)
        ! so that i=1 at the beginning of the zone and j is the time in 
        !bctime(j)
        if(k .eq. nion +1 .and. te_var_edge .ne. 0)then  !variable edge zone for te
           if(nbctim .gt. 1)then 
              call find1 (i1, i2, timel, bctime, nbctim)
              if(i1*i2 .eq. 0)
     .        call STOP ('subroutine set_boundary_condition te',1)
              !first get the index value, te_index:
              if(i1 .ne. i2)then
                 call linear (fix_edge_te(i1),fix_edge_te(i2),
     .                   fun,bctime(i1),bctime(i2), timel)
              else
                 fun = fix_edge_te(i1)
              endif
              if(te_var_edge .eq. 1)then
                  te_index = INT(fun)
              else  ! te_var_edge .eq. -1
                 !find the index of the r value closest to the normalized flux value fun:
                 call find1 (i1r, i2r, fun, roa, nj)
                 if(i1r*i2r .eq. 0)
     .           call STOP ('subroutine set_boundary_condition te',2)
                 !pick closest value:
                 diff1 =ABS(fun -roa(i1r))
                 diff2 =ABS(fun -roa(i2r))
                 if(diff1 .ge. diff2)then
                    te_index = i2r
                 else
                    te_index = i1r
                  endif
              endif
              if(freeze_te_index)te_index =te_index_save
              !next get the profile at this time
              iprofnbr =  kprim+kimp+2  !used by intrp via bcondspl
              if(i1 .eq. i2)then
                 !get profile of te over roa at time in bctime(i2)
                 !copy this profile into bctime_zone
                 call intrp (-1,-1,rtein(1,i1),tein(1,i1),knotste(i1),
     .               roa,yprof1,nj)
                 do j=te_index,nj
                   bctime_zone(j,kion+1) = yprof1(j)
                 enddo
              else
                !get two spline profiles, one at time bctime(i1) and the oter at time
                !bctime(i2) defined over the roa grid. Linearly interpolate these 
                !two profiles in time at each grid point from te_index to nj to get the
                !required profile and save it in bctime_zone:
                call intrp (-1,-1,rtein(1,i1),tein(1,i1),knotste(i1),
     .             roa,yprof1,nj)
                call intrp (-1,-1,rtein(1,i2),tein(1,i2),knotste(i2),
     .             roa,yprof2,nj)
                do j = te_index, nj
                  call linear (yprof1(j),yprof2(j),
     .                   fun,bctime(i1),bctime(i2), timel)
                  bctime_zone(j,kion+1) = fun
                enddo
             endif
           else         !nbctim = 1
             print *,'this option requires that nbctim .gt. 1 '
             call STOP ('subroutine set_boundary_condition te nbctim',1)
           endif !     (nbctim branches )
           ub(k) = bctime_zone(nj,kion+1)  !reset from above for consistency
         endif   !(te_var_edge .ne. 0)


        if(k .eq. nion +2 .and. ti_var_edge .ne. 0)then
          if(nbctim .gt. 1)then
           call find1 (i1, i2, timel, bctime, nbctim)
           if(i1*i2 .eq. 0)
     .        call STOP ('subroutine set_boundary_condition ti',1)
           !first get the index value, ti_index:
           if(i1 .ne. i2)then
              call linear (fix_edge_ti(i1),fix_edge_ti(i2),
     .                   fun,bctime(i1),bctime(i2), timel)
           else
              fun = fix_edge_ti(i1)
           endif
           if(ti_var_edge .eq. 1)then
                  ti_index = INT(fun)
           else  ! te_var_edge .eq. -1
             !find the index of the r value closest to the normalized flux value fun:
              call find1 (i1r, i2r, fun, roa, nj)
              if(i1r*i2r .eq. 0)
     .        call STOP ('subroutine set_boundary_condition ti',2)
              !pick closest value:
              diff1 =ABS(fun -roa(i1r))
              diff2 =ABS(fun -roa(i2r))
              if(diff1 .ge. diff2)then
                 ti_index = i2r
              else
                 ti_index = i1r
              endif
           endif
           if(freeze_ti_index)ti_index =ti_index_save
           !next get the profile at this time
           iprofnbr =  kprim+kimp+3  !used by intrp via bcondspl
           if(i1 .eq. i2)then
              !get profile of ti over roa at time in bctime(i2)
              !copy this profile into bctime_zone
              call intrp (-1,-1,rtiin(1,i1),tiin(1,i1),knotsti(i1),
     .             roa,yprof1,nj)
              do j=ti_index,nj
                 bctime_zone(j,kion+2) = yprof1(j)
              enddo
           else
              !get two spline profiles, one at time bctime(i1) and the oter at time
              !bctime(i2) defined over the roa grid. Linearly interpolate these 
              !two profiles in time at each grid point from te_index to nj to get the
              !required profile and save it in bctime_zone:
              call intrp (-1,-1,rtiin(1,i1),tiin(1,i1),knotsti(i1),
     .             roa,yprof1,nj)
              call intrp (-1,-1,rtiin(1,i2),tiin(1,i2),knotsti(i2),
     .             roa,yprof2,nj)
              do j = ti_index, nj
                  call linear (yprof1(j),yprof2(j),
     .                   fun,bctime(i1),bctime(i2), timel)
                  bctime_zone(j,kion+2) = fun
              enddo
           endif
           else         !nbctim = 1
             print *,'this option requires that nbctim .gt. 1 '
             call STOP ('subroutine set_boundary_condition ti nbctim',1)
           endif
           ub(k) = bctime_zone(nj,kion+2)  !reset from above for consistency
        endif




        rot_index = nj ! needed in trplot
        if(k .eq. nk .and. iangrot .eq. 1 .and. 
     .                      rot_var_edge .ne. 0)then
          if(nbctim .gt. 1)then
           !determine values in the  edge zone at this time by linear 
           !interpolation from angrot_edge_zone
           call find1 (i1, i2, timel, bctime, nbctim)
           if(i1*i2 .eq. 0)
     .        call STOP ('subroutine set_boundary_condition rot',1)
           !first get the index value rot_index:
           if(i1 .ne. i2)then
              call linear (fix_edge_rot(i1),fix_edge_rot(i2),
     .                   fun,bctime(i1),bctime(i2), timel)
           else
              fun = fix_edge_rot(i1)
           endif
           if(rot_var_edge .eq. 1)then
                  rot_index = INT(fun)
           else  ! rot_var_edge .eq. -1
             !find the index of the r value closest to the normalized flux value fun:
              call find1 (i1r, i2r, fun, roa, nj)
              if(i1r*i2r .eq. 0)
     .        call STOP ('subroutine set_boundary_condition w',2)
              !pick closest value:
              diff1 =ABS(fun -roa(i1r))
              diff2 =ABS(fun -roa(i2r))
              if(diff1 .ge. diff2)then
                 rot_index = i2r
              else
                 rot_index = i1r
              endif
           endif
           if(freeze_rot_index)rot_index =rot_index_save
           !next get the profile at this time
           iprofnbr =  kprim+kimp+6   !used by intrp via bcondspl
           if(i1 .eq. i2)then
              !get profile of angrot over roa at time in bctime(i2)
              !copy this profile into bctime_zone
              call intrp (-1,-1,rangrot(1,i1),angrotin(1,i1),
     .                    knotsang(i1),roa,yprof1,nj)
              do j=rot_index,nj
                 bctime_zone(j,kion+4) = yprof1(j)
              enddo
           else
              !get two spline profiles, one at time bctime(i1) and the oter at time
              !bctime(i2) defined over the roa grid. Linearly interpolate these 
              !two profiles in time at each grid point from rot_index to nj to get the
              !required profile and save it in bctime_zone:
              call intrp (-1,-1,rangrot(1,i1),angrotin(1,i1),
     .                     knotsang(i1),roa,yprof1,nj)
              call intrp (-1,-1,rangrot(1,i2),angrotin(1,i2),
     .                    knotsang(i2),roa,yprof2,nj)
              do j = rot_index, nj
                  call linear (yprof1(j),yprof2(j),
     .                   fun,bctime(i1),bctime(i2), timel)
                  bctime_zone(j,kion+4) = fun
              enddo
           endif
           else         !nbctim = 1
             print *,'this option requires that nbctim .gt. 1 '
             call STOP ('subroutine set_boundary_cond. angrot nbctim',1)
           endif
           ub(k) = bctime_zone(nj,kion+4)  !reset from above for consistency
        endif
      end do      !   k=1,nk loop 


      deallocate (yprof1, STAT = istat)
      if(istat .ne. 0)
     .          call deallocate_error("yprof1,set_bc",0,istat)
      deallocate (yprof2, STAT = istat)
      if(istat .ne. 0)
     .          call deallocate_error("yprof2,set_bc",0,istat)
c
      call check_zero_bc
      return
c
      end







      subroutine set_error_field(timel)
c----------------------------------------------------------------------------
c INPUT:
c timel        desired time of interpolation

c INPUT (through tordlrot.i):
c  nt_mgbr   number of entries in berrqn_mgbr

c time_mgbr(j)  times (in sec) of berr values in berrqn_mgbr
c berrqn_mgbr(j)  in Gauss, error field value at qn surface at times in
c               time_mgbr

c  OUTPUT (into tordlrot.i)
c berrqn_val     in Gauss if time is out of range of time_mgbr then
c              berrqn_val =berrqn_val(1) is returned if time is less
c              than time_mbgr(1) 
c              berrqn_val =berrqn_val(nt_mgbr) is returned if time is greater
c              than time_mbgr(nt_mgbr) 
c 
c
c  the drag down torque term associated with this error field is determined
c  in subroutine magnetic_drag,call from subroutine source (splitting the
c  caluclations between two subroutines is done for efficiency and also
c  because q is not known until source is called
c---------------------------------------------------------------------------
c
      USE param
      USE tordlrot
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c      include 'param.i'
c      include 'tordlrot.i'
      integer *4 j
      real *8 timel

      if(nt_mgbr .eq. 0)return
      call interp1 (timel,time_mgbr, nt_mgbr,berrqn_mgbr, berrqn_val)

      return
      end







      subroutine showme(istat,iter)
c-----------------------------------------------------------------------
c     subroutine dumps predictor/corrector values for analysis
c----------------------------------------------------------------------
      USE param
      USE solcon
      USE soln
      USE numbrs
      USE mesh
      implicit  integer (i-n), real*8 (a-h, o-z)
      LOGICAL,SAVE :: idef
      DATA idef / .FALSE. /

c     save the last cycle of predictor corrector calculations
      if(istep .eq. 'pred')then
         call getioun(iounit,17)
         idef = .TRUE.
         open(unit = iounit,file = 'predcor.txt',status ='REPLACE')
         write(iounit,'(2x,i5,2x,i5,4(2x,1pe16.8),a)')iter,nj,dt,dtt,
     .                                             theta,time,istep
          do j=1,nj
            write(iounit,'(2(2x,1pe14.6))')r(j),usave(nion+1,j)
            write(iounit,'(2(2x,1pe14.6))')r(j),usave(nion+2,j)
            write(iounit,'(2(2x,1pe14.6))')r(j),u(nion+1,j)
            write(iounit,'(2(2x,1pe14.6))')r(j),u(nion+2,j)
          enddo
      elseif(istat .ne. 1)then
         if(idef)THEN
          write(iounit,'(2x,i5,1pe16.8)')iter,time
          do j=1,nj
            write(iounit,'(2(2x,1pe14.6))')r(j),u(nion+1,j)
            write(iounit,'(2(2x,1pe14.6))')r(j),u(nion+2,j)
          enddo
         ENDIF
      else
         !corrector steps are done for this delta t
         ! close up shop
         !write(iounit,'("closing unit,time = ",1pe14.8)')time
         IF(idef)THEN 
            close(unit = iounit)
         ENDIF
         idef = .FALSE.
      endif

      return
      end
c
c
c

      subroutine tport
c
      USE param
      USE fusion
      USE aid_newton  !freeze_alpha_exp
      USE io 
      USE ions
      USE neut
      USE solcon
      USE soln
      USE transp
      USE mhdpar
      USE nub
      USE nub2
      USE rf
      USE glf23, only : write_glf_namelist
      USE mhdgrid
      USE xptor_sim
      USE extra
      USE numbrs
      USE mesh
      USE verbose
      USE sourc
      USE machin
      USE nub4
      USE tfact
      USE geom
      USE flags
      USE tordlrot
      USE soln2d
      USE bd_condtn
      USE mixcom
      USE kinetic_efit,              ONLY :  wrt_kin_efit_naml,
     .                                       kine_message
      USE rhog
      USE nonlin
      USE nbi_restart,               ONLY : beam_prof_init
      USE monitr
      USE oloss_dat
      USE  gcnmp_interface,          ONLY : gcnmp_driver

      USE gpsi
      USE island
      USE pelcom
      use pellet_mod

      USE P_Nfreya_12_interface,      ONLY : use_P_Nfreya

      USE statistics,               ONLY : collect_stats,start_timer,
     .                                     stop_timer,descrip_mon,
     .                                     p_nf_mon_set,p_nf_index,
     .                                     nubeam_mon_set,ech_mon_set,
     .                                     last_mon_index,
     .                                     nubeam_index,ech_index,
     .                                     elapsed_time,o12_index


      implicit  integer (i-n), real*8 (a-h, o-z)
c
      character rcs_id*63
      save      rcs_id
      data      rcs_id /
     ."$Id: cray306.f,v 1.112 2013/09/03 18:53:21 stjohn Exp $"/
c
c ----------------------------------------------------------------------
c --- subroutine TPORT is driver for transport portion of ONETWO code.
c --- evolve the transport equations from t = time to t = time + dteq if
c --- equilibria are being calculated,
c --- otherwise evolve               from t = time to t = timmax
c ----------------------------------------------------------------------
c


      include 'imsl.i'

      include 'quit.i'

      include 'sxrcom.i'

D     include 'mpif.h'
      include 'rebut.i'



c
      external FLUSH
c
      logical    impfh
      integer    beamconvg , dont_stop,nn
      real*4     singloid ! single precision variable
      real*8     zero, one, kev_per_joule
      dimension  curp(kj), enep(kj), enpp(kj,kion), enap(kj),
     .           enbp(kj,ke,kb), etorp(kj), psisp(kj), rbpp(kj),
     .           qp(kj), wap(kj), wbp(kj,ke,kb), ppbp(kj,ke,kb)
      dimension  tnion(2),tnerr(2),tnerrl(2),tnrerr(2)
c
      integer counter,kine_stat
      character *40 profile
      character *16  kin_mesg 
      parameter (icparamx = 6)
      real *8 cparamt(icparamx)
      REAL *8 pos_def
      data    counter,nn /0,0/
      data cparamt / 0.0d0,0.01d0, 0.1d0, 0.4d0,0.7d0,1.0d0 /
      pos_def = 10._DP*epsilon(0.0_DP)
      if(steady_state .gt. 0.0 .and. diffeq_methd .eq. 2)then
         cparamt(1) = .1d0
         cparamt(2) = .4d0
         cparamt(3) = .7d0
         cparamt(4) = 1.0d0
      endif
      first_step = 1
      myid = 0
      master = myid
D     call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr) !get processor id
c
c ----------------------------------------------------------------------
c istop is set to 1 if any of the following events occur
c
c 1. dt     <  dtmin
c 2. n      >   nmax
c 3. time   > timmax
c 4. timmax <= time0
c 5. nmax   <= 0
c ----------------------------------------------------------------------
c
c
c ----------------------------------------------------------------------
c initialize solution variables
c ----------------------------------------------------------------------
c

      kin_mesg = "none    " 
      kine_stat = 1
      timabs = 0.0D0
      get_typf   = 0    !initialization for nonlinear solution method 
      istop      = 0
      dont_stop  = 0    ! 1 means don't stop (hard wired for now)
!      dont_stop  = 1
      one        = 1.0
      kev_per_joule = 0.62415064e16
      info_conv = 1   !changed only if  newton solver is called
      dtsum  = 0.0
 
      istep    = 'init'
      ilastp   = 0
      itimav   = 0
      iconvg   = 0
      imix     = 0
      idterr   = 0
      ihalve_nl = 0      !used by non linear solver
      impfh    = implicit_fh
      timsav   = time
      eqtim0   = time    ! time of last equilibrium calculation
      losscalc = 1       ! turn on orbit calc in subroutine GET_ORBLOSS
      iborb_save = iborb ! may change below due to mcgo
      dtt      = 0       ! required for mhd cases since source tests
c                          on dtt value for beam slowing down
      dtmb1 = time
      dtmb2 = time+dt ! dtmb1,2 for transp beam

c         quantities set in nubeam, but note that nubeam is not
c         called during a corrector step:
      IF(use_nubeam)then
             beam_data%pwf_tot_intg =0.0
            do j=1,beam_data%nbeam 
                beam_data%pinja(j) =0.0
            enddo
       ENDIF
c
c ----------------------------------------------------------------------
c for test case, adjust profiles so that rbp should not evolve in time.
c then evolve the rbp equation anyway, as a check
c ----------------------------------------------------------------------
c
****  if (ltest_code .eq. 1 .and. time .eq. time0)  call fixprofile
c
c ----------------------------------------------------------------------
c subroutine SAVE copies u into usave for nk dependent variables
c ----------------------------------------------------------------------
c
*%%%% write (6, *) '%%%% calling subroutine SAVE'
      call savesol (u, usave, nk, nj, kk)

c      call stop('c306,l892 ',1)
c
c ----------------------------------------------------------------------
c set up boundary conditions at t = time.
c here time is the start of the run or the start of a new
c equilibirum/transport cycle.
c Load ub(k) boundary value (at plasma edge) for dependent
c variable k at t = time
c bc(i,k) input boundary condition (i = time index ,k=dependent
c variable index )
c ----------------------------------------------------------------------
c
      call set_boundary_condition (time)
      call set_error_field(time)
c
c ----------------------------------------------------------------------
c update copies vector u into the appropriate vectors en,te,ti,etc.
c note en(j,k),j = 1..nj,k=1...nion,where nion = nprim+nimp
c ----------------------------------------------------------------------
c
      write (nmhddat, '("starting TPORT, time = ", 1pe14.8)') time
      write (nmhddat, '("rbp(2-5) = ",4(2x,1pe14.6))')   (rbp(j),j=2,5)
      write (nmhddat, '("curden(1-4) = ",4(2x,1pe14.6))')
     .                                                (curden(j),j=1,4)

c     set en(1..nj,1..?), te(1.nj),ti(1..nj),rbp(1..nj),angrot(1..nj)
c     to values saved in u(:,1..nj)
      call update (u, en, te, ti, rbp, nk, nj, kj, kk,iangrot,angrot)

      call curcalc (rbp, fcap, hcap, gcap, r, curden, curpar_soln, nj)


      write (nmhddat, '("after update+curcalc")')
      write (nmhddat, '("rbp(2-5) = ",4(2x,1pe14.6))')   (rbp(j),j=2,5)
      write (nmhddat, '("curden(1-4) = ",4(2x,1pe14.6))')
     .                                                (curden(j),j=1,4)
      write (nmhddat, '("hcap(1-4) = ",4(2x,1pe14.6))') (hcap(j),j=1,4)
      write (nmhddat, '("dhdt(1-4) = ",4(2x,1pe14.6))') (dhdt(j),j=1,4)
      write (nmhddat, '("gcap(1-4) = ",4(2x,1pe14.6))') (gcap(j),j=1,4)
      write (nmhddat, '("dgdt(1-4) = ",4(2x,1pe14.6))') (dgdt(j),j=1,4)
      write (nmhddat, '("fcap(1-4) = ",4(2x,1pe14.6))') (fcap(j),j=1,4)
      write (nmhddat, '("dfdt(1-4) = ",4(2x,1pe14.6))') (dfdt(j),j=1,4)
      write (nmhddat, '("r(2-5) = "   ,4(2x,1pe14.6))') (r   (j),j=2,5)
      write (nmhddat, '("q(1-4) = "   ,4(2x,1pe14.6))') (q   (j),j=1,4)
      write (nmhddat, '("ub(1-4) = "   ,4(2x,1pe14.6))') (ub (j),j=1,4)
      write (nmhddat, '("ub(3-5) = "   ,4(2x,1pe14.6))') (ub (j),j=3,5)
      write (nmhddat, '("psir(1-4) = "   ,4(2x,1pe14.6))') 
     .                                               (psir   (j),j=1,4)
      write (nmhddat, '("Rbp(nj-3,nj) = "   ,4(2x,1pe14.6))') 
     .                                            (rbp   (j),j=nj-3,nj)
c
c ----------------------------------------------------------------------
c RHOMESH calculates the cap quantities and the various rho grids
c that are treated explicitly in time.
c The last time these quantities were determined by an MHD calculation
c was at time timcap (the corresponding values are in the cap0 vectors,
c timcap is set in subroutine PREPAR).
c fcap and hcap can be determined at any instant without
c doing  a new MHD calculation and hence are treated
c implicitly, in subroutine FHCALC, called from RHOMESH.
c in order for fhcalc to be able to calculate gcap,hcap and fcap
c the densities, temperatures, and poloidal field must be available in
c en,te,ti,rbp
c If user selected to apply boundary condition inside the actual
c plasma edge then pick it up here for the first time by loading
c the appropriate value into ub.
c ----------------------------------------------------------------------
c
      call copya(r,r_mesh,nj) ! for adaptive grid and rho_edge < 1 calcs
      call rhomsh(time)
      if (tportvb .ge. 2)  call dump_values ('top of tport')
c
c ----------------------------------------------------------------------
c FHUPDATE sets the fcap0,and hcap0 parameters equal to the current
c values, fcap, and hcap respectively
c ----------------------------------------------------------------------
c
      if (implicit_fh) then
        call fhcalc(0,1)
        call fhupdate
      end if
c
c ----------------------------------------------------------------------
c PROPEL injects pellets if ipelet .ne. 0
c ----------------------------------------------------------------------
c
*%%%% write (6, *) '%%%% calling subroutine PROPEL'
      call propel
c
c ----------------------------------------------------------------------
c SPECIFY sets the profiles for dependent variables run in analysis
c mode if nbctim .gt. 1 (otherwise profiles are not changed)
c both u and the profiles en,te,ti,rbp and angrot are updated.
c ne and zeff are updated to t = time only if inenez .ne. 0 and ttweak = 0
c ----------------------------------------------------------------------
c
*%%%% write (6, *) '%%%% calling subroutine SPECIFY'
      call specify
c
c ----------------------------------------------------------------------
c ZEN calculates ene and zeff if no impurities (nimp = 0) are present
c otherwise zen finds dz/dte and dz/dt for impurities and
c     a) if inenez    = 0, ene and zeff are calculated
c     b) if inenez .ne. 0, en(j,k),j = 1..nj,k=1..nprim+nimp
c        is determined and u(k,j) is updated
c ----------------------------------------------------------------------
c
*%%%% write (6, *) '%%%% calling subroutine ZEN'
      call zen
c
c ----------------------------------------------------------------------
c --- save the profiles at t = time
c ----------------------------------------------------------------------
c
*%%%% write (6, *) '%%%% calling subroutine SAVE'
      call savesol  (u, usave, nk, nj, kk)  ! ZEN may have updated vector u
*%%%% write (6, *) '%%%% calling subroutine COPYA'
      call copya (ene, enesav, nj)
c
c ----------------------------------------------------------------------
c if pellet injection is active the pressure profile will be changed.
c hence recalculate the cap parameters
c ----------------------------------------------------------------------
c
      if (ipelet .ne. 0 .and. implicit_fh) then
*%%%%   write (6, *) '%%%% calling subroutine FHCALC'
        call fhcalc(0,1)
*%%%%   write (6, *) '%%%% calling subroutine FHUPDATE'
        call fhupdate
      end if
c
c ----------------------------------------------------------------------
c zero quantities which will be calculated below
c
c               ******************************
c NOTE: if tport is called multiple times due to mhd coupling
c then the zeroing of these quantities implies that certain tolerances
c such as delrf and delte will be set to 1.0 which will cause
c the RF packages to be called if relrf is small enough.
c Setting relrf to a number greater than 1 will shutoff the rf heating
c amd current drive since these arrays are zeroed here and are not
c recalculated when delrf,delte < relrf.
c               ******************************
c ----------------------------------------------------------------------
c
      if (ieq .lt. 2) then
          call zeroa (enneu(1,1),nj)
          call zeroa (enneu(1,2),nj)
          call zeroa (tineu,nj)
          call zeroa (ennub,nj)
          call zeroa (tenub,nj)
          call zeroa (enrf,kj*krf)
          call zeroa (terf,kj*krf)
          call zeroa (tirf,kj*krf)
          call zeroa (qrfes,kj*krf)
          call zeroa (qrfis,kj*krf)
          call zeroa (currfs,kj*krf)
          call zeroa (jtor_eccd,kj*krf)
          call zeroa (dnedt,nj)
!          call zeroa (dnidt,kj*nprim) 
          call zeroa (dnidt,kj*nion)
          call zeroa (dtidt,kj)
          call zeroa (dtedt,kj)
          call zeroa (dpidt,kj)
          call zeroa (dpedt,kj)
          call zeroa (ssaw,kj*kprim)
          call zeroa (ssawe,nj)
      end if 
      ipfail(1) = 0
      ipfail(2) = 0
      if (timav .eq. 0.0)  go to 2120
      timabs = ABS (timav)
      timdif = timmax - time0
      if (timabs .gt. timdif) timav = timdif*timav/timabs
      dtsum  = 0.0
      iavneu = 0
      do 2115 j=1,nj
      do 2110 k=1,nk
 2110 uav(k,j) = 0.0
      do 2112 ib=1,nbeams
      do 2112 ie=1,3
      enbav(j,ie,ib) = 0.0
      wbav(j,ie,ib)  = 0.0
 2112 ppbav(j,ie,ib) = 0.0
      enaav(j)       = 0.0
 2115 waav(j)        = 0.0
 2120 if (w2mix .eq. 0.0 .and. w3mix .eq. 0.0)  go to 2122
      if (jsxr .eq. 0 .or. jsxr .eq. 1)  write (nqik, 9015)
      if (jsxr .eq. 2                 )  write (nqik, 9017)
 9015 format ('1   time   epste  dtecal  srad3m  srad3p   s3cal',
     . '   epsti  dtical   trcal   ddntm   ddntp  fuscal')
 9017 format ('1   time   epste  dtecal    s18m    s18p  s18cal',
     . '   epsti  dtical   trcal   ddntm   ddntp  fuscal')
c
 2122 if (n .gt. 0)  go to 3000  ! n = 0 means first transport time step
c
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c take care of initial time
c if the beam is on at this time (signified by ibeam = 3) then
c we must do a special startup iteration to
c        a)converge the fast ion density (the fast ions see
c          only the thermal part of the ion density as a target.
c          but the thermal ion density must be adjusted to allow
c          for the fast ions).
c        b)bootstrap related calculations in diffuse use
c          this information,which must be made self consistent
c          before we proceed
c        c)mcgo may be involved as well
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::HSJ:::::::::::::::
c



c--------------------------------------------------------------------------
c     a test section.  load usave,u ,te,ti,etc with profiles
c     read from  (usually xptor generated results ) file.
c     That way we calculate  the sources (pbeam,qrf, etc)
c     as was done above,  with the
c     initial profiles given in inone but then switch  to
c     the  profiles  read from the file as the inital conditions.
      If(test_xptor .eq. 1)then
         call xptor_init
      endif
c--------------------------------------------------------------------------

      call copya(r,r_mesh,nj)  ! for adaptive grid calcs
*%%%% write (6, *) '%%%% calling subroutine RHOMSH'
      call rhomsh(time)
      if (implicit_fh) then
*%%%%   write (6, *) '%%%% calling subroutine FHCALC'
        call fhcalc (0, 1)
*%%%%   write (6, *) '%%%% calling subroutine FHUPDATE'
        call fhupdate
c
c       set dfdt,dhdt = 0 if corresponding parameter is implicit
c
*%%%%   write (6, *) '%%%% calling subroutine TIMEDERIV'
        call timederiv (dt, n, fcap, fcap0, dfdt, nj)
*%%%%   write (6, *) '%%%% calling subroutine TIMEDERIV'
        call timederiv (dt, n, hcap, hcap0, dhdt, nj)
      end if
c
****  imaxbeamconvg = 5      ! default set in INIT, also user-selectable
      beamconvg     = 0                               ! integer variable
      beamerror     = 1.0e6
c

      beam_iteration = .false.

c      write(888,FMT='("use_nubeam = ",l5 )')use_nubeam  ! 8888899999

      beam_startup: if(use_nubeam) then
c-----
           ! start collecting time info here
             IF(.NOT. nubeam_mon_set)THEN  
                     nubeam_index = last_mon_index+1
                     last_mon_index = nubeam_index
                     nubeam_mon_set = .TRUE.
               descrip_mon(nubeam_index) = "Nubeam elapsed  time"
c               write(888,FMT='("start nubeam_indexset in c306 ",i5)')
c     .                nubeam_index  ! 88888889999999
             ENDIF
             start_timer(nubeam_index) =.TRUE.
             stop_timer(nubeam_index)  =.FALSE.
             CALL collect_stats(nubeam_index)
             start_timer(nubeam_index) =.FALSE.        
c             write(888,FMT='("start nubeam ",1pe14.6,x,i5)')
c     .             elapsed_time(nubeam_index),nubeam_index  ! 88888889999999
c-----
           nubeam_evolve = 1
           !using nubeam, if this is not a nubeam restart case then
           !we cant get the beam consistent with the thermal ion
           !density until the beam is turned on so 
           !there is nothing to be done here:
cJMP       if(nubeam_restart == 1) call beam_prof_init(time)
           

c           write(888,FMT='("nubeam_version = ",i10)')nubeam_version
           if (nubeam_version .lt. 201107) then 
              if((nubeam_restart == 1).or.(nubeam_restart ==-1)) then
                 call beam_prof_init(time)
              endif
           else

              call nubeam_beam_startup()


           endif
              stop_timer(nubeam_index)  =.TRUE.
              CALL collect_stats(nubeam_index)
              stop_timer(nubeam_index)  =.FALSE.
c              write(888,FMT='("stop nubeam 3610 ",1pe14.6,x,i5)')
c     .                   elapsed_time(nubeam_index),nubeam_index 
      else beam_startup   ! possibly iterate freya to consistency with
                          ! thermal ion densities
c          set the beam convergence. monte carlo noise enters here.
c          the minimum achievable error is governed by the number of
c          ions run in FREYA.
c          5% is a good number in most cases (npart .ge. 50000) 
c          otherwise we will
c          exit after imaxbeamconvg passes are made
c
           if (relaxden_err .le. 0.0) then
               relaxden_err = 0.10
               if (npart .ge. 50000)  relaxden_err = 0.05
           end if
           icallnub       = 1           ! forces call to neutral beam package
           no_beam_fusion = 1           ! turn off fusion calcs during beam..
c                                  ..iteration
           beam_iteration = .true.      ! let other packages(eq rf ) know that this
                                   ! is   beam iteration phase
c
c          turn off MCGO calculations during this phase
c          (see below if this is a snapshot case):
c
           if (iborb .eq. 3)  iborb = 1
c
           do while (beamconvg .lt. imaxbeamconvg  .and.
     .          beamerror .gt. relaxden_err)
c
                counter = counter + 1 ! count how many times DIFFUS is called
c     
                cparam  =1.0d0
                icparam = 1
                if (tportvb .gt. 0)
     .              print *,'calling diffuse  during beam iterations'
                call diffus (xi_include)
                if (tportvb .gt. 0)
     .              print *,'calling source  during beam iterations'


                call source    !beam gets called from source
                if (tportvb .gt. 0)
     .              print *,'done source  during beam iterations'
c
                if (tportvb .gt. 0)
     .              print *,'calling fluxx  during beam iterations'
                call fluxx
                if (tportvb .gt. 0)
     .             print *,'done fluxx  during beam iterations'
c
c
c
c                 zen recalculates ion densities including 
c                 enbeam when inenez = 1
c
                if (tportvb .gt. 0)
     .             print *,'calling zen  during beam iterations'
                call  zen
                if (tportvb .gt. 0)
     .               print *,'done  zen  during beam iterations'
c
c               if beam is not on (ie ibeam .ne. 3), and this is an old
c               time independnet beam run (ie time_dep_beam = 0 )
c               then iterations are not necessary
c         
                if (ibeam .ne. 3 .and. time_dep_beam .eq. 0)  go to 50

c
c               if user doesn't want to do these iterations exit
c               note that iterate_beam = false is set automatically
c               if time_dep_beam mode of running is selected.
c

               if (.not. iterate_beam)  go to 50
               
c
c               don't iterate the beam after the first equilibirum cycle
c
                if (ieq .gt. 1)  go to 50
                reldifnemax = 0.0
                reldifnamax = 0.0
                if (inenez .eq. 0) then
c
c                   if inenez = 0 (meaning ene and zeff are calculated) use the
c                   change in electron density as a measure of convergence
c
                    do j=1,nj
                      if (enesav(j) .gt. 0.0)
     .                reldifne   = ABS ((ene(j)-enesav(j))/(enesav(1)))
                      reldifnemax = MAX (reldifnemax, reldifne)
                      if (reldifnemax .eq. reldifne) then
                        jrhomax = j
                        prtden  = ene(j)
                        prtbden = enbeam(j)
                      end if
                      ene(j) = relaxden*ene(j) +
     .                    (1.0-relaxden)*enesav(j)   ! relax solution
                    end do
                    call copya (ene, enesav, nj)    ! update for next pass
                end if
c
                if (inenez .ne. 0) then
c
c               if inenez .ne. 0 (meaning primary and impurity
c               densities are calculated from ene and zeff) use the change
c               in the primary ion species corresponding to the beam (which
c               is en(j,ibion),j = 1,2..nj, ibion=1 or 2)
c               as a measure of convergence
c
                  if (ibion .gt. 0)  ij = ibion
                  if ( ifus .lt. 0) then
                    ij = 0
c
c                   use the density which is run in analysis mode
c                   to determine the beam convergence.
c
                    do j=1,nprim
                      if (itran(j) .eq. 0 .and. j .eq. id)  ij = id
                      if (itran(j) .eq. 0 .and. j .eq. it)  ij = it
                    end do
                    if (ij .eq. 0)
     .                call STOP ('subroutine TPORT: IJ = 0', 212)
                 end if
                 do j=1,nj
                   if (en(j,ij) .ne. 0.0)
**** .             reldifna    = ABS ((en(j,ij)-usave(ij,j)) /
**** .                                            en(j,ij))
     .             reldifna    = ABS ((en(j,ij)-usave(ij,j)) /
     .                                            en(1,ij))
                   reldifnamax = MAX (reldifnamax, reldifna)
                   if (reldifnamax .eq. reldifna) then
                     jrhomax = j
                     prtden  = en(j,ij)
                     prtbden = enbeam(j)
                     write (ncrt, '(" en, usave =", 2(2x, 1pe12.2))')
     .                              en(j,ij), usave(ij,j)
                   end if
                   en(j,ij) = relaxden*en(j,ij)
     .                    + (1.0-relaxden)*usave(ij,j)  ! relax solution
                 end do
c
c                update for next pass
c
                 call redate (u,en,te,ti,rbp,nk,nj,kj,kk,iangrot,angrot)
                 call savesol   (u, usave, nk, nj, kk)
             end if


             beamerror = MAX (reldifnemax, reldifnamax) ! one arg will be 0
             beamconvg = beamconvg + 1             ! limit number of passes


c
c            beam pressure changes fcap and hcap
c
             if (implicit_fh) then
*%%%%           write (6, *) '%%%% calling subroutine FHCALC'
                call fhcalc (0, 1)
             end if
             istep = 'pred'        ! so beam deposition can be recalculated
             write  (ncrt, 55) beamconvg, beamerror, relaxden_err,
     .                      jrhomax, relaxden, prtden, prtbden
   55        format (' beam iteration  =', i5,
     .           '  error/specified =', 1pe12.2, ' /', 1pe12.2,
     .           '  at spatial index j =', i5 /
     .            ' relaxation fctr =', 1pe12.2,
     .           '  thermal den =', 1pe12.2, '  beam den =', 1pe12.2)
        end do  ! end do while with beamconvg

      endif beam_startup


   50 no_beam_fusion = 0      ! allow fusion calculations now 
      beam_iteration = .false.


      icallnub = 0            ! don't force a call to the beam package
      

c
c     but if MCGO is selected we have to force a call to beam package
c     now because it may not get called again (depending on how the
c     user set the inputs)
c
      if (iborb_save .eq. 3)  icallnub = 1
      iborb = iborb_save        ! enable MCGO if it was set on input
      call diffus (xi_include)  !HSJ call added  11/07/03
      print *,'calling source  beam iterations are done'

      call source

      if(bp0_ic(1:LEN_TRIM(bp0_ic)) == bp0_icv(2))then !bp0_icv(2) = 'analytic'
         call reinit_bp0
         call curcalc (rbp, fcap, hcap, gcap, r, curden,
     .                                          curpar_soln, nj)
         call redate (u,en,te,ti,rbp,nk,nj,kj,kk,iangrot,angrot)
         call savesol   (u, usave, nk, nj, kk)
         bp0_ic = 'eqdsk'                !only inital guess, dont do it again
      endif
c
c --- save the modified profiles in by calling SAVIT
c --- (written to disk)
      call savit                                          ! HSJ, 4/24/95
c
c --- check if time0-beamon is sufficiently large that
c --- asymptotic density was achieved
c

      if (ibslow .eq. 1 .and. (time0 - beamon(1)) .ge. 0.0
     .                                  .and. .not. use_nubeam ) then
        taupbmax = 0.0                      ! fast ion slowing down time
        do    jb=1,nbeams
          do  ic=1,3
            do j=1,nj
              taupbmax = MAX (taupbmax, taupb(j,ic,jb))
            end do
          end do
        end do
        irslow = (time0 - beamon(1)) / taupbmax
      else
        irslow = 10 ! taupb was not calculated so set irslow arbitrarily
      end if
c
      if (irslow .lt. 5)  write (ncrt, 60) irslow
   60 format (' WARNING: beamon(1) is not set sufficiently below'      /
     .             9x, ' time0 to yield asymptotic fast ion density' /
     .             9x, ' number of slowing down times =', i5)
 
      istep    = 'init'
      iborb    = iborb_save        ! enable MCGO if it was set on input
      icallnub =  0
      call savesol  (u, usave, nk, nj, kk)
      call copya (ene, enesav, nj)
      if (implicit_fh) then
        call fhcalc (0, 1)
        call fhupdate
        call timederiv (dt, n, fcap, fcap0, dfdt, nj)
        call timederiv (dt, n, hcap, hcap0, dhdt, nj)
      end if
****  if (twkfar .gt. 2.0)  implicit_fh = .false.

c
c ----------------------------------------------------------------------
c update fast ion quantities
c ----------------------------------------------------------------------
c
      do 2130 ib=1,nbeams
      do 2130 ie=1,3
      call copya ( enb(1,ie,ib), enbsav (1,ie,ib), nj)
      call copya (  wb(1,ie,ib), wbsav  (1,ie,ib), nj)
      call copya (  sb(1,ie,ib), sbsav  (1,ie,ib), nj)     ! HSJ 5/12/97
      call copya (pprb(1,ie,ib), pprbsav(1,ie,ib), nj)
 2130 call copya ( ppb(1,ie,ib), ppbsav (1,ie,ib), nj)
      call copya (enalp        , enasav          , nj)
      call copya ( walp        ,  wasav          , nj)
c
*%%%% write (6, *) '%%%% calling subroutine IMPSRC'
      call impsrc  ! adds qdelt and flux in Faraday's law to source term
      sn2d = 0.0
      do i=1,2
        snaddt(i) = 0.0
        flxmod(i) = 0.0
        tnerr (i) = 0.0
        tnerrl(i) = 0.0
        tnrerr(i) = 0.0
        if (ineut(i) .ne. 0) then  ! neutral species i is present
*%%%%     write (6, *) '%%%% calling subroutine CHKCON'
c
c         chkcon calculates
c         tnion0(i),(  total # of ions of type i )
c         for analysis mode variables tnerr(i) and tnrerr(i) =0 is returned
c         for simulation current # - (previous +added ions) if recycling
c         or current# -previous # for constant # ions option
c         either way tnerr is error and tnrerr is relative error
c
          call chkcon (i, tnion0(i), tnerr(i), tnrerr(i))
        end if
      end do
c
c ----------------------------------------------------------------------
c print out and plot results at initial time.
c print out only if specified in prtlst.
c summary page will be printed if ilastp = 1 and iyoka > 0
c ----------------------------------- HSJ ---- 1/8/99 ---- start
c

c
      call specify   ! set u for itran=0 variables. These variables can
c                      change due to specification of input profiles
c                      at times given in bctime(1...nbctim)
c                      copy u into  en,te,ti,rbp,angrot
c                      zeff and ene are update if inenez .ne. 0  HSJ
c
      call zen       ! get zeff and ene if inenez = 0
      call savesol  (u, usave, nk, nj, kk)  ! ZEN may have updated vector u
c
c ----------------------------------- HSJ ---- 1/8/99 ---- end
c
 
      call fiziks
c
c get the experimental profiles if this is tdem mode
c
c ----------------------------------------------------------------
c
      zero = 0.0
      call non_inductive_cd (zero)
c

c      if (time .eq. time0 .and. prtlst(1) .ne. time0)  go to 2160
      if ( ( (time .eq. time0) .and.( prtlst(1) .ne. time0) ) .AND.
     .             (timmax .GT. time0) )go to 2160
c      added above line to fix problem with timmax <= time0 issue
c      skipping inital call to out HSJ 7/27/2011 
      if (timmax .le. time0)  ilastp = 1

      if(myid .eq. 0)call out                     ! initial time point
  
 2160 ihead = 1
      if (tportvb .ge. 1)  write (*, '(" calling   TRPLOT A at time = ",
     .                                   1pe14.6)') time
      if(myid .eq. 0)call trplot (3)
      if (tportvb .ge. 1 .and. myid .eq. 0 )  
     .            write (6, '(" time after trplot " /
     .            1pe20.12," dt  = ",1pe20.12 )') time,dt
      if (tportvb .ge. 1 .and. myid .eq. 0) 
     .                         write (*, '(" done with TRPLOT")')
      if (inenez .eq. -1)  inenez = 0     ! compute ene and zeff..
c                                         ..from en in future
      if (inubplt .eq. +1 .and. myid .eq. 0) then
*%%%%   write (6, *) '%%%% calling subroutine NUBPLT'
        call nubplt
      end if
      

c
c --------------------------------------------------------------------------
c   create kinetic efit output:
        kine_stat =1
        kine_message ='after_init'
        call wrt_kin_efit_naml(1,kine_stat)



c
c ----------------------------------------------------------------------
c stop if time-dependent solution is not desired (i.e., snapshot mode)
c ----------------------------------------------------------------------
c

      if (timmax .le. time0 .or. nmax .le. 0) then
         istop = 1
         call write_restart_profs
      endif

c      if(timmax .le. time0      .or.  nmax .le. 0           )THEN
c        write(940,FMT='("leaving tport because timmax .lt. time0 ",
c     .  2(x,1pe16.8),//)')timmax,time0 ! 888888899999
c      endif

      if (timmax .le. time0                 )  return
      if (                       nmax .le. 0)  return
c
c  end of initial time calculations
c

c     a test section.  load usave,u ,te,ti,etc with profiles
c     read from  (usually xptor generated results ) file.
c     That way we calculate  the sources (pbeam,qrf, etc)
c     as was done above,  with the
c     initial profiles given in inone but then switch  to
c     the  profiles  read from the file as the inital conditions.
      If(test_xptor .eq. 1)then
         call xptor_init
      endif




c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c evolve the diffusive system of equations.
c a)predictor step:
c CHEKDT updates timnew to the value it will have after the corrector
c step has converged. (i.e., set timnew = time+dt) dt will be cut
c down if necessary to avoid crossing of source switching, print or
c plot times, switches will be set to indicate that sources are
c on or off for times greater than timenew.
c label 3000 is entry point for each new transport time step
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      ncycles_dt = 3
      idt_half = ncycles_dt !if this reaches 0 then dt will be doubled

 3000 if(n .eq. 0)ddebugsave = ddebug(3)
      iter_vloop = 0
      tipts2 = tipts2 + 1 !used in solve_newton if wrt_nwt > 1
      if(n .gt. 0) nubeam_evolve = 1       !allow calls to nubeam 
                              !(effective only if if nubeam is to be used)
                              !was set to -1 if  en_init was  called
      nubeam_dt = nubeam0_dt  !reset time interval between nubeam calls to
                              !original input value in case checkdt changed it.
      if(time + nubeam_dt .gt. timmax)  !avoid running nubeam beyond the 
     .nubeam_dt = timmax-time +time_tol !end of the run
      first_step = 1
      write (6, '(a, i4)') ' transport time step number', n
 
      if(MODULO(n,int(ddebug(11))) .eq. 0)
     .  ddebug(3) =1.0                            !No relaxation at beginning of this time step
      nback = 0 ! passed by io.i, used if ddebug(3) < 0
      ihalve_nl =0
      if (time .ge. timecrit)  dtmax = dtmaxcrit
      freeze_nl = 0
c      freeze_alpha_exp = 0 
c     if time step is becoming small freze something in the nonlinear solution model
      if(dt .lt. 2*dtmin .and. n .gt. 1)freeze_nl =  1
      if (tportvb .ge. 1)  write (6, '(" time before chekdt " /
     .            1pe20.12," dt  = ",1pe20.12 )') time,dt
   
      call chekdt (time, dt, n, timprt, mprt, timplt, mplot,
     .             irf, dthold, dtevmin, timpel,
     .             ipelet, npel, timrfp, iplot, iprt,
     .             irfplt, codeid, time0, dteq, eqtim0,taus,timmax)
      if (tportvb .ge. 1)  write (6, '(" time after chekdt " /
     .            1pe20.12," dt  = ",1pe20.12 )') time,dt

      if (tportvb .ge. 1)  write (6, '(" try to move solution " /
     .            " from time"  1pe20.12," to time = ",1pe20.12 /
     . " using pivot value ", 1pe20.12)') time,time+dt,time+theta*dt
      if (tportvb .ge. 1 .and. myid .eq. 0)
     .  call dump_values ("starting new time step in TPORT")
c
c     to check TDEM mode
c
      call save_bpo_tdem
****  call copya (r, r_mesh, nj)  ! for adaptive grid calculations
c
c --- note that in the predictor step u = usave so theta weighting
c --- of u and usave done in subroutine ABCG is the same as setting
c --- theta = 0 (which is what we want, since the predictor step is
c --- fully explicit)
c
      dtt   =  theta * dt   ! HSJ 08/14/04define central time increment for predictor
      
      istep = 'pred'
      dcparam = dcparam_in
      icut_con = 0
c
c     tnerrl is not defined until predictor step is done
c     hence flxmod =tnerrl means that we assume for the current
c     step the same rate of adding particles as in the previous step
c     this will be iterated if necessary:                   HSJ
c
      do 3010 i=1,2
 3010 if (ineut(i) .ne. 0)  flxmod(i) = tnerrl(i)
      do i=1,2
        flxmod(i) = 0.0
      end do                                         ! HSJ 2/9/96
c
      timnew = timsav + dtt    ! after predictor step
      nnew   = n
      if ( timav .eq. 0.0 .or. timnew .le. timmax-timabs)  go to 3015
      if (iavneu .eq. 1                                 )  go to 3015
      icallneu = 1
      iavneu   = 1
c
c ----------------------------------------------------------------------
c the following calls to RHOMSH, DIFFUS, SOURCE and FLUXX generate the
c information needed to define the coeffcient matrix A and source
c vector G at the initial time (i.e.,  t = time = timsav)
c at this time en,te,ti,rbp and angrot are equal to the respective
c values in the vector u (due to a previous call to update). however
c below the values of en,te,ti,rbp and angrot will be set to the
c values at the central time point. hence the routines diffus,source,
c and fluxx must use en,te,ti,rbp and angrot rather than the vector
c u when  istep = 'pred'  or  istep = 'corr'.
c fcap, and hcap may be treated implicitly in time. this can be
c done even if ifixshap = 1 (i.e., no MHD evolution) by using the
c flux surface average of the G.S. equation.
c ----------------------------------------------------------------------
c
 3015 continue
*%%%% write (6, *) '%%%% calling subroutine RHOMSH'
      call copya (r,r_mesh,nj) ! for adaptive grid calcs
      call rhomsh (time)   ! elapsed time for mhd-dependent parameters..
c                          ..is time-timcap

c
c --- get fcap, hcap at time t+dtt by extrapolation, use the
c --- derivative from the previous time step dfdtsv, etc.
c
      if (implicit_fh) then
*%%%%   write (6, *) '%%%% calling subroutine FHEXTRAP'
        call fhextrap (dtt, fcap0, fcap, hcap0, hcap,
     .                 dfdtsv, dhdtsv, nj, implicit_fh)
      end if
c
c ----------------------------------------------------------------------
c use the extrapolated cap quantities in the calculation of matrix A
c and source vector G for the predictor step
c ----------------------------------------------------------------------
c
c     set time-dependent w1-4typ (take-your-pick transport coefficient
c     multipliers)
c
      if (wtyp .ne. 0.0)  call wtypt (timnew)
c
      counter = counter + 1 ! count how many times DIFFUS is called
****  if (counter .gt. 59)  call STOP ('subroutine TPORT: DEBUG', 902)!
c
*%%%% write (6, *) '%%%% calling subroutine DIFFUS'
      if(continuation_method .eq. 1 )then
          icparam = 1
c          cparam =0.0 !only implemented for glf23 and typ at this time
          cparam = cparamt(icparam)
      else
         cparam =1.0d0
      endif
      IF(diffeq_methd .NE. 3)call diffus (xi_include)
      if(glf_debug .gt. 0)  call STOP('sub TPORT: glf_debug stop',0)

      call source

      if(write_glf_namelist .eq. 3) return ! return to get iterdb file written
c              stop would be cleaner but doesnt write iterdb file
c    &      call STOP('TPORT: write_glf_namelist specified  stop',0)
*%%%% write (6, *) '%%%% calling subroutine FLUXX'
      call fluxx
c
c --- get the boundary conditions at the central time point
c     (done using linear interpolation)
c
      call set_boundary_condition (timnew)
      call set_error_field(timnew)
c
c ----------------------------------------------------------------------
c generate the solution at the central time point,
c t = timsav + theta * dt (= timnew)
c ----------------------------------------------------------------------
c
c      write(940,FMT='("diffeq_method in tport =",i5)')diffeq_methd
      if (diffeq_methd .eq. 0 .or. diffeq_methd .eq. 2 ) then
         !use standard predictor step for both solution methods
         if (rho_edge .lt. 1.0) then
           nj_save=nj
           nj=nj_rho_edge
           call copya(ub,ub_save,kk)     ! save ub
           call copya(ub_rho_edge,ub,kk) ! set bc at rho_edge using ub
         end if
c
           call solve ! sets u = usave for analysis mode variables
c                     provided there are some. If there are no
c                     simulation variables (i.e., nkt =0) then solve
c                     returns without doing anything!!!
c
           first_step =0 
           if(tportvb .ge. 2)call showme(0,0)
      else if(diffeq_methd .eq. 1)then
           call solve_ml2
           first_step =0 
c          call solve_ml
c
c          get solution from te, etc. into u:
c
           call redate (u, en, te, ti, rbp, nk, nj, kj, kk, 
     .                                       iangrot,angrot)
           nnew = n + 1
      else if(diffeq_methd .eq. 3)then ! setup and call gcnmp
           print *,'before gcnmp call , time,timmax  =',time,timmax
           CALL gcnmp_driver
           
           IF(ABS(time-timmax) .lt. pos_def)time =timmax
           call redate (u,en,te,ti,rbp,nk,nj,kj,kk,iangrot,angrot)
           !roundoff problem encountered here:

           timnew = time
           nnew    = n + 1
           dtt = timnew - timsav
           IF(tportvb .GE. 1 )
     .      WRITE(ncrt,'(a)')'completed GCNMP interaction' 
           print *,'after gcnmp , time =',time,timmax,dtt
           IF(dtt .LT. pos_def)CALL STOP('gcnmp failure',1)
      end if
      

c     freeze_alpha_exp = 1 ! freeze alpha stab factor after initial step
ch      if(fix_edge_te .lt. nj .or. fix_edge_ti .lt. nj
ch     .     .or. fix_edge_rot  .lt. nj) 
ch     . call reset_edge_values(te,ti,angrot,fix_edge_te,
ch     .              fix_edge_ti,fix_edge_rot,u,kk,nk,nj,iangrot)

c
      if (rho_edge .lt. 1.0) then
        nj = nj_save
        call copya(ub_save,ub,kk)           ! restore ub
      end if
c
      if (tportvb .ge. 2) then
        write (6, '(a)')  ' after SOLVE, predictor:'
        call testsol (1, nj)
      end if
c
c ----------------------------------------------------------------------
c check whether any elements of solution vector
c (excluding rbp and toroidal rotation) are negative.
c for the first profile encountered that is negative cheku returns with
c
c     ihalve = 1        delav = 0.0
c     delmax = 1.0      kmax  = k of (first) profile checked that failed
c     jmax   = points to r(jmax) where problem occured
c     delit  = 0.0
c     iter   = 0
c
c ----------------------------------------------------------------------
c
*%%%% write (6, *) '%%%% calling subroutine CHEKU'
      call cheku (u, nk, nj, kk, delav, delmax, kmax, jmax,
     .            delit, iter, ihalve, iangrot)
      if (ihalve .eq. 1 .or. ihalve_nl .eq. 1) then
c        if(diffeq_methd .eq. 0)then
             if(kmax .le. nprim) then
                profile ='primary ion density '
             else if(kmax .le. nion)then
                profile ='impurity ion density'
             else if(kmax .eq. nion+1)then
                profile ='electron temperature'
             else if(kmax .eq. nion+2)then
                profile ='ion temperature'
             else if(kmax .eq. nk-iangrot)then
                profile ='poloidal B field'
             else if(kmax .eq. nk .and. iangrot .eq. 1)then
                profile ='rotation speed'
             endif
             write (6, *)
     .           'predicted solution is negative for profile k =', kmax
             write(6,*)'problem profile is :',profile
             write (6, *)  'at grid point j=', jmax
             write (6, *)  'u(kmax,jmax) = ', u(kmax,jmax)
             write (6, *)  
     .             'time step will be cut in half current dt=',dt
c         else
             !non linear solver initial guess is poor. fix latter
c             call STOP ('subroutine TPORT: inital soln guess fails', 0)
c         endif
        if (dt .gt. 2.01*dtmin) go to 5000  ! new line HSJ
****    go to 5000                            old line
        if (dont_stop .eq. 0) go to 5000
c
c            fudge the solution a bit to avoid stopping
c
             u(kmax,jmax)=usave(kmax,jmax)
             write (6, *)  'reset u(kmax,jmax) = ', u(kmax,jmax)
             ihalve=0
             ihalve_nl = 0
      end if
c
c ----------------------------------------------------------------------
c update the ...sav quantities to the values at the predictor time
c ----------------------------------------------------------------------
c
 3180 if (ibeam .eq. 3) then
        do   ib=1,nbeams
          do ie=1,3
            call copya ( enb(1,ie,ib),   enbsav(1,ie,ib), nj)
            call copya (  wb(1,ie,ib),   wbsav(1,ie,ib), nj)
            call copya (  pprb(1,ie,ib), pprbsav(1,ie,ib), nj)
            call copya ( ppb(1,ie,ib),   ppbsav(1,ie,ib), nj)
          end do
        end do
      end if
      call copya ( enalp, enasav, nj)
      call copya ( walp,  wasav, nj) ! added HSJ
      if (diffeq_methd ==  1  .OR.
     .              diffeq_methd == 3)  go to 3344   !method of lines,gcnmp
                                                     !do not use following
c
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c b)corrector steps
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      dttheta = dtt
      dtt     = dt                 ! the full time step for corrector
      iconvg  = 0
      istep   = 'corr'
      time    = timnew             ! central time (after predictor step)
      timnew  = timsav + dtt       ! after corrector step     HSJ
      nnew    = n + 1
      iter    = 0
      ddebug(3) = ddebugsave       !corrector iterations will be  underrelaxed if
                                   !ddebugsave is .lt. 1
c
c ----------------------------------------------------------------------
c update sets en,te,ti,rbp,angrot to the values in vector u,which was
c just calculated in solve,at the central time point
c solve sets u=usave for analysis mode variables  (in subroutine EXPAND)
c only if at least one variable is run in simulation . HSJ
c ----------------------------------------------------------------------
c
      call update (u, en, te, ti, rbp, nk, nj, kj, kk, iangrot, angrot)
c
****  write (6, *) '    ti(1), predictor',     ti(1)
****  write (6, *) 'qbeami(1), predictor', qbeami(1)
*%%%% write (6, *) '%%%% calling subroutine SPECIFY'
      call specify   ! set u for itran=0 variables. These variables can
c                      change due to specification of input profiles
c                      at times given in bctime(1...nbctim)
c                      copy u into  en,te,ti,rbp,angrot
c                      zeff and ene are update if inenez .ne. 0  HSJ
c
*%%%% write (6, *) '%%%% calling subroutine ZEN'
      call zen       ! get zeff and ene if inenez = 0
c
c ----------------------------------------------------------------------
c update explicit cap and mesh quantities to central time point:
c rhomesh is not included in the corrector iterations because the
c relevant parameters do not change on subsequent corrector steps.
c ----------------------------------------------------------------------
c
*%%%% write (6, *) '%%%% calling subroutine RHOMSH'
      call rhomsh (time)
c
c ----------------------------------------------------------------------
c get new estimates of the cap parameters based on the solution
c found at the central time point.
c use these new estimates to generate a new time derivative
c for the implicit cap parameters
c ----------------------------------------------------------------------
c
      if (implicit_fh) then
        call fhcalc (0, 1)    ! use en, etc. and ene
        call timederiv (dttheta, nnew, fcap, fcap0, dfdt, nj)
        call timederiv (dttheta, nnew, hcap, hcap0, dhdt, nj)
        call copya     (fcap, fcapctr, nj)
        call copya     (hcap, hcapctr, nj)
      end if
c
c ----------------------------------------------------------------------
c dndt calculates d(ene)/dt, d(angrot)/dt, ,dtedt,dtidt,dpedt,dpidt,and
c d(en)/dt. the vector u is current set for t = time and the vector
c usave is currently set for t = timsav. ene is currently set for t = time,
c enesav is currently set for t = timsav. (t = time is the central time point)
c the electron density is not changed during iterations of the corrector
c (ene will  change during the corrector iterations if inenez=0 . This
c change is neglected)        HSJ
c ----------------------------------------------------------------------
c
      dtinv = 1.0 / (time - timsav)     ! central-start
*%%%% write (6, *) '%%%% calling subroutine DNDT'
      call dndt (dtinv, ene, enesav, kj, kk, nj, nion, u, usave,
     .           dnedt, dnidt, dtedt,dtidt,dpedt,dpidt,iangrot, 
     .           nk, dangrot)
      call non_inductive_cd (time-timsav)
c
c --- get boundary conditions at time=timnew             HSJ
c
      call set_boundary_condition (timnew)
      call set_error_field(timnew)
c
      ibacks = 0  !for continuation method
 3200 iter      = iter + 1       ! start each corrector step here
      itercorct = iter           ! passed to DIFFUS by way of rebut.i
      if (iter .gt. itmax .and. continuation_method .eq. 0 ) then
        write (6, *)  'iter > itmax, so halve the value of dt'
        if(dont_stop .eq. 0) go to 5000               ! halve the time step
        !ignore the iteration limit itmax 
        go to 162
      else if(iter .gt. itmax .and. continuation_method .eq. 1)then
        !half the increment in the continuation parameter
        icparam = icparam -1
        ibacks = ibacks -10   !dont use ibacks below
        if(icparam .lt. 2) go to 5000
c        cparam =  (cparamt(icparam) + cparam)*0.5
c        dcparam = dcparam*.5
        iter = 0
        icut_con = icut_con+2
        if(icut_con .gt. 1 .and. ibacks .lt. -1)then !half cparam increment only onece
          write(6,*)  'cparam too large,cutting back timestep'
          go to 5000
        endif
         call recall (u, u_continue, nk, nj, kk)  ! copy u_continue into u
         call update (u, en, te, ti, rbp, nk, nj, kj, kk, iangrot, 
     .                                                     angrot)
         call specify
         call zen
      endif
c
c ----------------------------------------------------------------------
c --- copy the central time values into the cap vectors for use in
c --- source, etc. (if iter = 1 then fcapsv is already set to fcap, etc)
c ----------------------------------------------------------------------
c
      if (implicit_fh .and. iter .gt. 1) then
*%%%%   write (6, *) '%%%% calling subroutine COPYA'
        call copya (fcapctr,fcap,nj)
*%%%%   write (6, *) '%%%% calling subroutine COPYA'
        call copya (hcapctr,hcap,nj)
      end if

c
c ----------------------------------------------------------------------
c --- generate the information needed to get matrix A and source
c --- vector G at the central time point,using en,te,ti,rbp,angrot
c ----------------------------------------------------------------------
c
      counter = counter + 1 ! count how many times DIFFUS is called
****  if (counter .gt. 59)  call STOP ('subroutine TPORT: DEBUG', 903)!
c
*%%%% write (6, *) '%%%% calling subroutine DIFFUS'

         call diffus (xi_include) 

      call source

      write (6, *) '%%%% calling subroutine FLUXX'
      call fluxx

      call non_inductive_cd (dtt)   ! dtt is full step here
c
c ----------------------------------------------------------------------
c now advance the solution from t = timsav  to t = timnew (= timsav+dt):
c (i.e., from the previous time point to the new time point, using the
c centered coefficent matrix A and source vector G).
c ----------------------------------------------------------------------
c
      if (rho_edge .lt. 1.0) then
        nj_save=nj
        nj=nj_rho_edge
        call copya(ub,ub_save,kk)     ! save ub
        call copya(ub_rho_edge,ub,kk) ! set bc at rho_edge using ub
      end if
      IF(diffeq_methd .eq. 0)THEN
           call solve        ! sets u = usave for analysis mode variables
           info_conv = 1     ! used to distinguish diffeq_methd = 2
      ELSEIF(diffeq_methd .eq. 2)THEN
         print *,"calling solve_newton,iter,itmax =",iter,itmax
         call solve_newton(info_conv)
      ELSE
         CALL STOP('sub tport, diffeq_methd error',1)
      ENDIF
 



ch      if(fix_edge_te .lt. nj .or. fix_edge_ti .lt. nj
ch     .                               .or. fix_edge_rot .lt. nj)
ch     .call reset_edge_values(te,ti,angrot,fix_edge_te,
ch     .                 fix_edge_ti,fix_edge_rot,u,kk,nk,nj,iangrot)
      if(tportvb .ge. 2 )call showme(0,iter)
      if (rho_edge .lt. 1.0) then
         nj=nj_save
         call copya(ub_save,ub,kk)    ! restore ub
      end if
      if (tportvb .ge. 2) then
        write (6, *) ' after solve in corrector'
        call testsol(1,nj)
      end if
c
c ----------------------------------------------------------------------
c chkcon returns tnion(i),the total number of ions of primary species
c i (= 1 or 2) currently in the system,
c tnerr(i) and tnrerr(i) are returned as zero if primary ions are run
c in anlysis mode. otherwise ??
c ----------------------------------------------------------------------
c
      do i=1,2
        if (ineut(i) .ne. 0) then
*%%%%     write (6, *) '%%%% calling subroutine CHKCON'
          call chkcon (i, tnion(i), tnerr(i), tnrerr(i))
          flxmod(i) =  tnerr(i)   ! HSJ 2/9/96
        end if
      end do


c
c ----------------------------------------------------------------------
c check whether any elements of solution vector
c (excluding rbp and toroidal rotation) are negative.
c for the first profile encountered that is negative cheku returns with
c
c     ihalve = 1        delav = 0.0
c     delmax = 1.0      kmax  = k of (first) profile checked that failed
c     jmax   = points to r(jmax) where problem occured
c     delit  = 0.0
c     iter   = 0
c
c ----------------------------------------------------------------------
c

      if(steady_state  .le. 0.0  .and. diffeq_methd
     .            .eq. 2)then
        !even if we didnt converge to the required degree 
        !(indicated by info_conv .ne. 1) we return he best
        !result obtained rather than quit.
        if(info_conv .ne. 1)then
           write(6,4321)SSQR,gradmax,tot_iters
           write(nout,4321)SSQR,gradmax,tot_iters
 4321       format(2x,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',/,
     .       2x,'convergence to specified degree not achieved',/,
     .       2x,'convergence results are :',/,
     .       2x,'total sum of squares of residuals (all equations):',
     .           1pe12.2,/,
     .       2x,'maximum gradient of any equation:',1pe12.2,/,
     .       2x,'total number of iterations taken :',i5,/,
     .       2x,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
         endif
      else

         call cheku (u, nk, nj, kk, delav, delmax, kmax, jmax,
     .            delit, iter, ihalve, iangrot)
         if (ihalve .eq. 1 .or. info_conv .ne. 1) then
            if(steady_state  .le. 0.0  .and. diffeq_methd
     .            .eq. 2)then
              !steady state newton solver didn't converge
              call STOP ('subroutine TPORT: steady state soltn failed',
     .                                                               0)
            else
               if(diffeq_methd .eq. 0)then
                 write (6, *)
     .           'corrector solution is negative for profile k =', kmax
                 write (6, *)  'at grid point j = ', jmax
                 write (6, *)  '   u(k,j), time = ', u(kmax,jmax), time
                 write (6, *)  'time step will be halved; current dt =',
     .                                                              dt
               else if(diffeq_methd .eq. 2)then
                 write (6, *)  'time step will be halved; current dt =',
     .                                                              dt
               endif
           endif
           if (dt .gt. 2.01*dtmin )  go to 5000  ! new line HSJ
****       go to 5000                            ! old line
           if (dont_stop .eq. 0) go to 5000
c
c            fudge the solution a bit to avoid stopping
c
            if(kmax .ne.0)u(kmax,jmax) = usave(kmax,jmax)
            ihalve       = 0
            ihalve_nl = 0
         end if
      endif

c     steady state solution using newton method- if it converged then 
c     finish up this run:
        if(steady_state .le. 0.0  .and. diffeq_methd
     .            .eq. 2)go to 3344
c
c     get new rho mesh if adaptive grid calcs are done:
c
****  call set_rho (time)
c
c ----------------------------------------------------------------------
c recalculate the implicit cap parameters using the solution vector
c u (which corresponds to timenew):
c get new time derivatives at the central time using the new values
c fcap, gcap, hcap, and the old values fcap0, gcap0, hcap0:
c ----------------------------------------------------------------------
c
      if (implicit_fh) then
*%%%%   write (6, *) '%%%% calling subroutine FHCALC'
        call fhcalc (0, 2)                       ! use vector u, and ene
*%%%%   write (6, *) '%%%% calling subroutine TIMEDERIV'
        call timederiv (dtt, nnew, fcap, fcap0, dfdt, nj)
*%%%%   write (6, *) '%%%% calling subroutine TIMEDERIV'
        call timederiv (dtt, nnew, hcap, hcap0, dhdt, nj)
      end if
c
c ----------------------------------------------------------------------
c calculate the maximum relative change in any of the implicit
c cap parameters, delcapc.
c we have the values stored in fcap,hcap for t = timnew
c and the values stored in fcapctr,hcapctr  for t = time
c (which is the central time point). subroutine CAPCHECK updates
c fcapctr, and hcapctr to the latest estimate.
c ----------------------------------------------------------------------
c
      delcapc = 0.0
      if (implicit_fh) then
*%%%%   write (6, *) '%%%% calling subroutine CAPCHECK'
        call capcheck (fcap0, hcap0, fcap, hcap, fcapctr, hcapctr,
     .                 nj, theta, delcapc)
      end if


c
c ----------------------------------------------------------------------
c calculate delit, the maximum relative difference between predicted and
c corrected solutions at the central time point. check only those
c dependent variables run in simulation mode. also update all
c simulation mode variables to the central time point,using the newly
c found solution from the corrector, which is stored in vector u.
c (before the call to CHECK_SOLUTION the old predicted values of
c simulation mode variables are stored in en,te,ti,rbp and angrot.
c after the call to CHECK_SOLUTION en,te,ti,rbp
c and angrot are updated to form a new predictor. u and usave
c are not changed. For total analysis (i.e., all itran(i)=0) mode the
c only effect of check_solution is to return delit = 0)
c the profiles at the central time are now corrected in CHECK_SOLUTION
c and the largest change, delit, is obtained.
c The call to GET_DENSITY is required for particle transport runs in
c which not all primary ions are run in simulation mode. To calculate
c correct flux for the simulation mode ion we need to have the correct
c density gradient for all primary ions, including those run in analysis
c mode ifus = -1 is only case of this type for now                   HSJ
c ----------------------------------------------------------------------
c
*%%%% write (6, *) '%%%% calling subroutine CHECK_SOLUTION'

      call check_solution (u, usave, en, te, ti, rbp, nk, nj, kj, kk,
     .                     r, delit, theta, itran, iangrot, angrot,
     .                     relaxsol, kprof, jloc, value_sol,
     .                     te_index,ti_index,rot_index)
c
      if (ifus .eq. -1) then ! only implemented for this case at present
        call correct_u (dtt)
        if (tportvb .ge. 2) then
          write (6, *) 'after CORRECT_U'
          call testsol (1, nj)
        end if
      end if
c
c ----------------------------------------------------------------------
c check for convergence of corrector (delit .le. relit).
c The corrector must also satisfy particle conservation
c (tnrerr .le. erneumax).  flxmod is the number of particles in error,
c which adjusts the recycling rate in SOURCE.
c DTLIM checks for signs of numerical instability,if the user
c requested it (by setting ilimdt = 1).
c dtlim returns idterr = 1 if
c               a) iter   = itmax, and
c               b) delit  > relit, and
c               c) icount = 2,     where icount is internal to dtlim
c                  icount is zeroed whenever iter = 1 and
c                  is incremented by 1 whenever delit .le. relit
c otherwise dtlim returns idterr = 0
c if idterr = 1 is set then dtmax will be decreased below.
c iconvg counts number of corrector iterations required to achieve
c particle conservation.
c ----------------------------------------------------------------------
c
      if (ilimdt .eq. 1) then
*%%%%   write (6, *) '%%%% calling subroutine DTLIM' 
        call dtlim (iter, delit, relit, itmax, idterr)
      end if
c

      igoto3200 = 0            ! assume corrector step not required  HSJ
      if (delit .gt. relit) then              ! iterate on the corrector
c         if(info_conv .ne. 1 .or. iteration_method .eq. 1)then
         if(diffeq_methd .eq. 0)then        !predictor- corrector  solution method
             if (tportvb .ge. 1) then
                write  (6, '(a, i8, 2e16.8)') ' iter, delit, relit =',
     .                                              iter, delit, relit
                write  (6, 160)  nameu(kprof), r(jloc),jloc,nj,value_sol
  160           format (' for profile ', a2, ' @ r =', f10.3,
     .            ' j/nj = ', i3, '/', i3,'   value = ', 1pe14.6)
                write (6,'("time =,cparam =",2(2x,1pe16.8))')time,cparam
                write  (6, '(a)')  ' corrector will be iterated'
             end if
****         go to 3200
             igoto3200         = 1 ! corrector step is required
             rmonu(jloc,kprof) = rmonu(jloc,kprof)+1.
             imontf            = imontf + 1
          else !newton solution method overide the check if newton solution  is converged
             if(info_conv .eq. 1) delit = .9*relit      
          endif
      end if
c
c     implicit cap parameters not converged
c

      if (implicit_fh .and. delcapc .gt. delcap) then
        if (tportvb .ge. 1)
     .    write (6, '(" delcap > delcapc ", 2(f10.2, 2x))')
     .                  delcapc, delcap
****    go to 3200
        igoto3200 = 1 ! corrector step is required
      end if
c
c     made it to here if delit < relit and cap parms are converged
c     hence lack of convergence (determined below) must be due to
c     (lack of) particle conservation. Note that during corrector
c     steps neucg2 is not called only renormalization with subroutine
c     NEUDEN is done. flxmod/dt is added to the incoming neutral flux
c     to yield total particle balance.                            HSJ
c
      iconvg = iconvg + 1
      do 3240 i=1,2
        if (ineut(i) .eq. 0)  go to 3240
****    flxmod(i) = flxmod(i) + tnerr(i)   ! HSJ 2/8/96
        flxmod(i) =             tnerr(i)
        tnerrl(i) =             tnerr(i)
 3240 continue
c
      do 3250 i=1,2
        if (      ineut(i)  .eq. 0       )  go to 3250
        if (     ipfail(i)  .ne. 0       )  go to 3250
        if (ABS (tnrerr(i)) .gt. erneumax) then
           if (tportvb .ge. 1)
     .        write (6, '(" tnrerr(i), erneumax = ", 2(1pe12.4))')
     .                      tnrerr(i), erneumax
c          go to 3200    ! iterate on corrector
           igoto3200 = 1 ! corrector step is required
           imonn       = imonn  + 1
           imontf      = imontf + 1
        end if
 3250 continue
      print *,'delit,relit,igoto3200 ',delit,relit,igoto3200
      if (igoto3200 .gt. 0)  go to 3200 ! perform corrector step



c-------for continuation method we must continue the correctors until
c-------cparam .eq. 1
      if(continuation_method .eq. 1 .and. cparam .lt. 1.0d0)then 
c            cparam_previous = cparam
c            cparam = cparam + dcparam
            icparam =icparam +1
            icparam = Min(icparam,icparamx)
            cparam = cparamt(icparam)
c            if(cparam .gt. 0.0d0)then
c               cparam = cparam *2.5
c            else
c               cparam =0.005
c            endif
            cparam = MIN(cparam,1.0d0)
            iter = 0
            call recall(u_continue, u, nk, nj, kk)  !u_continue, is set to u
            icut_con = 0
 
            go to 3200                        ! continue with correctors
      endif

      icparam = 1
                                           

c
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c corrector converged and conserved particles
c print out time step information
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
*%%%% write (6, *) '%%%% calling subroutine DELTA'
c     delta returns delav and delmax for both simulation and analysis
c     mode variables.
c
 3344 call delta_sol (u, usave, nk, nj, kk,  r, delav, delmax,
     .      kmax, jmax, iangrot,te_index,ti_index,rot_index,itran,
     .       analysis_check)
*%%%% write (6, *) '%%%% calling subroutine INFO'
****  if (counter .gt. 59) call STOP ('subroutine TPORT', 987)! MATCH
      call info  (nout, nqik, ihead, nnew, extime, timnew, dt, dtt,
     .            delav, delmax, kmax, jmax, delit, iter, ineu, inub,
     .            irfcalc, imix, iconvg, isecrem0)
****  if (counter .gt. 58) call STOP ('subroutine TPORT', 987)! NOMATCH
      imix = 0
c
c ----------------------------------------------------------------------
c check change in solution: 
c   even though we converged,
c   the time step will be halved if things change too rapidly
c the check is done even in total analysis mode since time derivatives
c of these quantities must be calculated (see fiziks for example).
c ----------------------------------------------------------------------
c
      if(diffeq_methd .ge. 2  .and. steady_state .le. 0.0 )go to 162

      if (ABS (delmax) .gt. relmax) then
        write  (6,   *)  'converged at time ', time, ' but:'
        write  (6,   *)  'delmax > relmax, so cut dt in half'
        write  (6,   *)  'delmax, relmax =', delmax, relmax
        write  (6, 161)   nameu(kmax), r(jmax)/r(nj)
  161   format (' for dependent variable ', a2, ' at r/a = ', f10.3)
        write  (6,   *)  'usave(k,j),u(k,j) =',usave(kmax,jmax),
     .                                         u(kmax,jmax)
c        write(940,FMT='("go to  5000 half the  time step")')
        go to 5000                                 ! halve the time step
      end if
c
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c made it through another time step;
c now for some physics calculations and general bookkeeping
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
 162  time = timnew
      if(steady_state .le. 0.0 .and. diffeq_methd .gt. 0)then
         call savesol  (u, usave, nk, nj, kk)
         time = timmax
         if(info_conv .eq. 1)then
         write(6,'("converged to steady state solution")')
         else
         write(6,'("not converged, best result is returned")')
         endif
         write(6,'(" there is no guarantee that " ,
     .             "this solution is unique")')
      else
         write (6, '(a, e14.4,2x,e14.4)')  
     .   ' converged at time,dt = ', time,dt
      endif
      if(tportvb .ge. 2 )call showme(1,iter)
c
c --- flush standard output channel ncrt for interactive viewing of results
c

      if (MOD (n, 5) .eq. 0)  call FLUSH (ncrt, iflush)

c
      n  = nnew
      nn = nn + 1
 
c
c ----------------------------------------------------------------------
c en, te, ti, rbp, and angrot were evaluated at the central time point
c during the corrector steps. we now update these profiles to
c the new time using subroutine UPDATE
c subroutine average averages the solution stored in u(k,j)
c ----------------------------------------------------------------------
c
      if(ddebug(12) .gt. 0) call average
      call update (u, en, te, ti, rbp, nk, nj, kj, kk, iangrot, angrot)
      if (tportvb .ge. 1)  write (*, '(" calling   SPECIFY")')
      call specify
      if (tportvb .ge. 1)  write (*, '(" done with SPECIFY")')
      if (tportvb .ge. 1)  write (*, '(" calling   ZEN")')
      if (tportvb .ge. 1)  call dump_values ('done with this time step')
      call zen           ! note that if inenez=0 ZEN has updated ene HSJ
      if(steady_state .lt. 1.e-5)call copya (ene, enesav, nj)      
      if (tportvb .ge. 1)  write (*, '(" done with ZEN")')
      dtinv = 1.0 / dt
      call dndt (dtinv, ene, enesav, kj, kk, nj, nion, u, usave,
     .           dnedt, dnidt,dtedt,dtidt,dpedt,dpidt, iangrot,
     .           nk, dangrot)
      if (inenez .eq. 0)                      ! HSJ 2/14/96
     .  call reset_fluxe (dnedt,sione,ssawe,r,ra,drr,hcap,nj,itimav)
      if (tportvb .ge. 2) then
        write (6, *)  ' after RESET_FLUXE'
        call testsol (1, nj)
      end if


c
c ----------------------------------------------------------------------
c --- set the cap0 parameters equal to the cap parameters
c --- and update the derivatives dfdtsv, etc as appropriate
c ----------------------------------------------------------------------
c
      if (implicit_fh)  call fhupdate
c
c ----------------------------------------------------------------------
c add in source terms that are handled implicitly
c in the difference equations (electron ion energy exchange term
c and part of flux in Faradays law)
c ----------------------------------------------------------------------
c
      if (tportvb .ge. 1)  write (*, '(" calling   IMPSRC")')
      call impsrc
      if (tportvb .ge. 1)  write (*, '(" done with IMPSRC")')
c
c ----------------------------------------------------------------------
c update number of particles added to system and
c calculate some output quantities
c ----------------------------------------------------------------------
c
      do 3526 i=1,2
        if (ineut(i) .eq. 0)  go to 3526
        snaddt(i) = snaddt(i) + snadd(i) + sngas(i) + sn2d
 3526 continue
c
      if (tportvb .ge. 1)  write (*, '(" calling   FIZIKS")')

      call fiziks
      if (tportvb .ge. 1)  write (*, '(" done with FIZIKS 1 ")')
      timsav = time
 
c
c get the experimental profiles if this is tdem mode
c ----------------------------------------------------------------
c
****  call non_inductive_cd (dtt)
c
c ----------------------------------------------------------------------
c calculate time-average of solution if time > timmax - ABS (timav)
c (time-1.0e-8 corrects for 'tol' used in chekdt)
c ----------------------------------------------------------------------
c

      if (time - 1.0e-8  .le.  timmax - timabs)  go to 3620
      if (dtsum .ne. 0.0)  go to 3605
      call savesol  (usave, uav0, nk, nj, kk)
      call copya (enesav, eneav0, nj)
      call copya (enbeam, enbav0, nj)
      call copya ( enalp, enaav0, nj)
 3605 dtsum = dtsum + dt
      do 3615 j=1,nj
      do 3610 k=1,nk
 3610 uav(k,j) = uav(k,j) + 0.5 * dt * (usave(k,j)+u(k,j))
      do 3612 ib=1,nbeams
      do 3612 ie=1,3
       enbav(j,ie,ib) =  enbav(j,ie,ib) + 0.5 * dt * ( enbsav(j,ie,ib)
     .                                                 +  enb(j,ie,ib))
        wbav(j,ie,ib) =   wbav(j,ie,ib) + 0.5 * dt * (  wbsav(j,ie,ib)
     .                                                 +   wb(j,ie,ib))
      pprbav(j,ie,ib) = pprbav(j,ie,ib) + 0.5 * dt * (pprbsav(j,ie,ib)
     .                                                 + pprb(j,ie,ib))
 3612  ppbav(j,ie,ib) =  ppbav(j,ie,ib) + 0.5 * dt * ( ppbsav(j,ie,ib)
     .                                                 +  ppb(j,ie,ib))
      enaav(j)        = enaav(j) +        0.5 * dt *(enasav(j)+enalp(j))
 3615  waav(j)        =  waav(j) +        0.5 * dt * (wasav(j)+walp(j))
c
c ----------------------------------------------------------------------
c carry out magnetic reconnection and plasma mixing calculation if wmix .ne. 0
c ----------------------------------------------------------------------
c
 3620 imix   = 0
      if (wmix .eq. 0.0)  go to 3650
      timmxo = timmix
      imslmd = 'mix'
c
      if (tportvb .ge. 1)  write (*, '(" calling   MIX")')
      call mix (atw, btor, curden, dt, ene, en, enalp, enb,
     .          etor, fcap, gcap, hcap,
     .          nbeams, nion, nj, nmix,  ppb, psis, q,
     .          r, rbp, rmajor, te, ti, time,
     .          walp, wb, zeff,
     .          curp, enep, enpp, enap, enbp, esaw, etorp, ppbp,
     .          psisp, qp, qmag, qsawe, qsawi, rbpp, rmixxx, rsxxx,
     .          wap, wbp)
      if (tportvb .ge. 1)  write (*, '(" done with MIX")')
c
      if (time .ne. timmix)  go to 3650
      imix   = 1
      rmixx  = rmixxx
      rsx    = rsxxx
      if (w2mix .eq. 0.0 .and. w3mix .eq. 0.0)  go to 3628
      dtecal = tem(1)-tep(1)
      dtical = tim(1)-tip(1)
      trcal  = 0.0
      if (dtecal .ne. 0.0)  trcal = dtical/dtecal
      if (  jsxr .eq. 0.0)  go to 3624
      call sxrcal (codeid, kappa, ene, jsxr, nj, nw, nh,
     .             p, pmax, psir, r, rmajor, rin, rmax, tem,
     .             rmhdgrid, zmhdgrid, zax, zmin, zmax,
     .             idiode, narray, namar, ndiode, roamin, sxrm)
      if (tportvb .ge. 1)  write (*, '(" calling   DDFUST")')
      call ddfust (tim, en, ddfusm, ddntm, time, bpol, totcur(1))
      if (tportvb .ge. 1)  write (*, '(" done with DDFUST")')
      call sxrcal (codeid, kappa, enep, jsxr, nj, nw, nh,
     .             p, pmax, psir, r, rmajor, rin, rmax, tep,
     .             rmhdgrid, zmhdgrid, zax, zmin, zmax,
     .             idiode, narray, namar, ndiode, roamin, sxrp)
      if (tportvb .ge. 1)  write (*, '(" calling   DDFUST")')
      call ddfust (tip, enpp, ddfusp, ddntp, time, bpol, totcur(1))
      if (tportvb .ge. 1)  write (*, '(" done with DDFUST")')
      do 3622 k=1,kar
      do 3622 i=1,kdi
        sxcal(i,k) = 0.5 * (sxrm(i,k)+sxrp(i,k))
 3622 if (sxcal(i,k) .ne. 0.0)
     .  sxcal(i,k) = (sxrm(i,k)-sxrp(i,k))/sxcal(i,k)
       s3cal = sxcal(3,1)
      s71cal = sxcal(2,3)
      s18cal = sxcal(2,2)
      fuscal = 0.5 * (ddntm+ddntp)
      if (fuscal .ne. 0.0)  fuscal = (ddntm-ddntp) / fuscal
c
 3624 if (w2mix .ge. 0. .and. w3mix .ge. 0.0)  go to 3625
      call tweak2 (dtecal, dtemix, dtical, dtimix, epste, epsti,
     .             fuscal, fusmix, s3cal, s3mix, s71cal, s71mix,
     .             s18cal, s18mix, trcal, trmix)
c

 3625 ddntm = ddntm * 1.0e-10
      ddntp = ddntp * 1.0e-10
      if (jsxr .ne. 2)  write (nqik, 9020)
     .  time, epste, dtecal, sxrm(3,1), sxrp(3,1),
     .  s3cal, epsti, dtical, trcal, ddntm, ddntp, fuscal
      if (jsxr .eq. 2)  write (nqik, 9020)
     .  time, epste, dtecal, sxrm(2,2), sxrp(2,2),
     .  s18cal, epsti, dtical, trcal, ddntm, ddntp, fuscal
 9020 format (13f8.3)
      ddntm = ddntm * 1.0e10
      ddntp = ddntp * 1.0e10
c
 3628 call copya (ene, ssawe, nj)
      dtmixi = 0.0
      if (time .ne. timmxo)  dtmixi = 1.0 / (time-timmxo)
      timmxo = time
      call copya (curp , curden, nj)
      call copya (enep , ene   , nj)
      call copya (enap , enalp , nj)
      call copya (etorp, etor  , nj)
      call copya (psisp, psis  , nj)
      call copya (rbpp , rbp   , nj)
      call copya (qp   , q     , nj)
      cconst = -1.0
      if (btor .lt. 0.0)  call multpl1 (q, nj, cconst)
      call copya (wap  , walp  , nj)
      if (w2mix .gt. 0.0)  call copya (tep, te, nj)
      if (w3mix .gt. 0.0)  call copya (tip, ti, nj)
      do 3630 i=1,nion
 3630 if (itran(i) .eq. 0)   call copya (en(1,i),enpp(1,i),nj)
      do 3632 i=1,nprim
      do 3632 j=1,nj
 3632 ssaw(j,i) = (enpp(j,i)-en(j,i))*dtmixi
      do 3635 i=1,nion
 3635 call copya (enpp(1,i),en(1,i),nj)
      if (w5mix .eq. 0.0)  go to 3641
      do 3640 ib=1,nbeams
      do 3640 ie=1,3
      call copya (enbp(1,ie,ib), enb(1,ie,ib), nj)
      call copya ( wbp(1,ie,ib),  wb(1,ie,ib), nj)
 3640 call copya (ppbp(1,ie,ib), ppb(1,ie,ib), nj)
c
 3641 call updatb (enb, kj, ke, nj, nbeams, one       , enbeam)
      call updatb ( wb, kj, ke, nj, nbeams, kev_per_joule,  wbeam)
      call zen
      do 3645 j=1,nj
 3645 ssawe(j) = (ene(j)-ssawe(j))*dtmixi
      call redate (u, en, te, ti, rbp, nk, nj, kj, kk, iangrot, angrot)
      if (tportvb .ge. 1)  write (*, '(" calling   FIZIKS")')
      call fiziks
      if (tportvb .ge. 1)  write (*, '(" done with FIZIKS 2 ")')
c
c ----------------------------------------------------------------
c get the experimental profiles if this is tdem mode
c ----------------------------------------------------------------
c
****  call non_inductive_cd (dtt)
c
c ----------------------------------------------------------------------
c write out magnetic island information
c ----------------------------------------------------------------------
c
 3650 if (wisl .le. 0.0)  go to 3700
      rsm21 = risl(1)-wdisl(1)/2.0
      rsp21 = risl(1)+wdisl(1)/2.0
      rsm31 = risl(2)-wdisl(2)/2.0
      rsp32 = risl(2)+wdisl(2)/2.0
c
      write  (nitre , '(1p4e15.4)')  rsm21, rsp21, rsm32, rsp32
      write  (nitrex,   8141)       (risl (i), i=1, 5)
      write  (nitrex,   8141)       (wdisl(i), i=1, 5)
      write  (nitrex,   8141)       (dpisl(i), i=1, 5)
      write  (nitrex,   8141)       (risl (i), i=6,10)
      write  (nitrex,   8141)       (wdisl(i), i=6,10)
      write  (nitrex,   8141)       (dpisl(i), i=6,10)
 8141 format (1p5e13.2)


 
c
c ----------------------------------------------------------------------
c impose constraint for stability to ideal ballooning modes
c if ibaloo .ne. 0 and plasma has a circular cross section
c ----------------------------------------------------------------------
c
 3700 if (ibaloo .eq.  0      )  go to 3800
      if (codeid .ne. 'onedee')  go to 3800
      call baloo
c
c ----------------------------------------------------------------------
c tweak various parameters if ttweak .ne. 0
c ----------------------------------------------------------------------
c
 3800 if (ttweak .eq. 0.0)  go to 3900
      call tweak1
c
c ----------------------------------------------------------------------
c simulate pellet injection at desired times
c ----------------------------------------------------------------------
c
 3900 call propel
      if (ipelet .eq. 1)  call redate (u,en,te,ti,rbp,nk,nj,kj,kk,
     .                                 iangrot,angrot)
c
c ----------------------------------------------------------------------
c dthold > 0.0 implies time step was reduced for
c informational purposes only; if this is the case
c then we will restore the time step to its
c original value
c ----------------------------------------------------------------------
c

      if (dthold .le. 0.0)  go to 4000
      dt = dthold
      go to 4020
c
c ----------------------------------------------------------------------
c if change solution < relmin double time step
c ----------------------------------------------------------------------
c
 4000 if (delmax .ge. relmin)  go to 4020
      if(idt_half .eq. 0)then
         dt = 2.0 * dt
         idt_half = ncycles_dt
      else
         idt_half = idt_half-1
      endif
      dt = MIN (dt, dtmax)
c
c ----------------------------------------------------------------------
c check timmax, nmax and dteq
c ----------------------------------------------------------------------
c
 4020 if (   n .ge. nmax .and. steady_state .gt. 0.0 )  go to 6000
      if (time .ge. timmax)  go to 6000
      if (codeid .eq. 'onedee')  go to 4040
      if (ifixshap .eq.1) go to 4040
           eqtime = time - eqtim0
           if (ABS (eqtime-dteq) .le. dtevmin+dtmin)  go to 6000
           if (             time .gt. eqtim0+dteq  )  go to 6000
c
c ----------------------------------------------------------------------
c if inubplt = 1 write out data for FREYA-type plots
c ----------------------------------------------------------------------
c

 4040 if (tportvb .ge. 1 .and. inubplt .eq. 1)
     .                     write (*, '(" calling   NUBPLT")')
      if (inubplt .eq. 1)  call nubplt
      if (tportvb .ge. 1 .and. inubplt .eq. 1)
     .                     write (*, '(" done with NUBPLT")')


c
c ----------------------------------------------------------------------
c call output routines if iprt = 1
c ----------------------------------------------------------------------
c
      if ( ignflg .eq. 1)  iprt = 1
      if (tportvb .ge. 1 .and. iprt .eq. 1)
     .                     write (*, '(" calling   OUT")')
 
      if (   iprt .eq. 1)THEN 
c      write(940,FMT='("c306,line 2677 ,totrf =",1pe12.4)')totrf ! 888888889999
c      write(940,FMT ='("c306 currf(nj/2) =",1pe12.4)')currf(nj/2) ! 888889999
         call out
c      write(940,FMT='("c306,line 2679 ,totrf =",1pe12.4)')totrf ! 888888889999
         call Statefile_proc
      endif

      if (tportvb .ge. 1 .and. iprt .eq. 1)
     .                     write (*, '(" done with OUT")')
 
c
c ----------------------------------------------------------------------
c record plot information (control generated file size  with steps_per_plot)
c ----------------------------------------------------------------------
c

      if (MOD (n, steps_per_plot) .eq. 0) then
        if (tportvb .ge. 1)  write (*, '(" calling   TRPLOTB")')
        iplot = 3  ! normally set in CHEKDT

        call trplot (iplot)
        if (tportvb .ge. 1)  write (*, '(" done with TRPLOT")')
 
      end if 
c
c ----------------------------------------------------------------------
c ibeam = 1 turn beam on; ibeam = 2 turn off beam
c irf   = 1 turn RF on;   irf   = 2 turn off rf
c ----------------------------------------------------------------------
c

      DO  k=1,krf
          if (irf(k) .eq. 1)  dt = 0.01*dtmax
          if (irf(k) .eq. 2)  dt = 0.01*dtmax
      ENDDO
      if (ipelet .eq. 1)  dt = 0.01*dtmax

      DO k=1,krf
        if (irf(k) .eq. 1)  irf(k) =  3
        if (irf(k) .eq. 2)  irf(k) = -1
      ENDDO


      IF(.NOT. use_P_Nfreya)THEN
!     Equivalent of following is done in source
!     for P_Nfreya

         IF(.not. use_nubeam )then
            if ( ibeam .eq. 1)  dt = 0.01*dtmax
            if ( ibeam .eq. 2)  dt = 0.01*dtmax
         ENDIF

  
         if ( ibeam .ne. 2)  go to 4210
c
c     The beam is now off so we zero the saved sources here. Just to
c     to make it a bit harder to understand we assign values to the
c     saved sources in subroutine SOURCE if the beam is still on. HSJ
c
         DO jb=1,nbeams
            DO ic=1,3
               call zeroa (  sbsav(1,ic,jb), nj)
               call zeroa (  qbsav(1,ic,jb), nj)
               call zeroa (spbrsav(1,ic,jb), nj)
               call zeroa ( spbsav(1,ic,jb), nj)
            ENDDO
         ENDDO
c
 4210          if (ibeam .eq. 1)  ibeam =  3
               if (ibeam .eq. 2)  ibeam = -1
c

      ENDIF

c
      if (ignflg .eq. 1) ignflg = -1
****  prado    = pradn
****  dprado   = dpradn
      icallneu = 0        ! do not force calls to NEUCG, FREYA and RF
      icallnub = 0
      icallrf  = 0
c
c ----------------------------------------------------------------------
c save current solution
c ----------------------------------------------------------------------
c
      call savesol  (u, usave, nk, nj, kk)
      call copya (ene,enesav, nj)
      do 4525 ib=1,nbeams
      do 4525 ie=1,3
      call copya ( enb(1,ie,ib),  enbsav(1,ie,ib), nj)
      call copya (  wb(1,ie,ib),   wbsav(1,ie,ib), nj)
      call copya (pprb(1,ie,ib), pprbsav(1,ie,ib), nj)
 4525 call copya ( ppb(1,ie,ib),  ppbsav(1,ie,ib), nj)
      call copya (        enalp,           enasav, nj)
      call copya (         walp,            wasav, nj)
c      write(940,FMT='("back to 3000 next imte step")')
      go to 3000                     ! go back and do the next time step
c
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c time step must be halved
c iconvg =0 if no corrector iterations were required to maintain
c particle conservation.
c We end up here if we didnt converge OR if the change in the converged
c solution was too large compared to the previous value. HSJ
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
 5000 if (iconvg .eq. 0)  icallneu = 1           ! forces call to neucg
      if (iconvg .eq. 0 .and. inub .ne. 0)
     .                    icallnub = 1           ! forces call to FREYA
      if (iconvg .eq. 0)  icallrf  = 1           ! forces call to RF pkg
      call showme(1,iter)
      time   = timsav
      call copya(r_mesh,r,nj)
      timnew = time
      call recall (u, usave, nk, nj, kk)
      call update (u, en, te, ti, rbp, nk, nj, kj, kk, iangrot, angrot)
      if (tportvb .gt. 0) then
         if (iter .gt. itmax) then
            jmax=1
            kmax=1
         end if
         write (6, '(" restored u, updated en, te, ti, rbp" /
     .               " u(kmax,jmax), time = ", 2(1x,1pe24.12))')
     .                 u(kmax,jmax), time
      end if
      call specify
      call zen
c
c --- set the cap parameters (and derivatives) back to the
c --- values they had before the time step was taken
c
      if (implicit_fh)  call fhrecall
c      if(ihalve_nl .gt. 0 .and. diffeq_methd .eq. 2)then
c         dt = dtmin  !go immediately to smallest time 
c      else
         dt = 0.5 * dt
         idt_half = idt_half+1
c      endif
      if (ilimdt .ne. 1)  go to 5010  !ilimdt is user controlled by namelist input
      if (idterr .eq. 0)  go to 5010
c
c --- try for integer number of milliseconds, and usually round down
c
      if (dt .gt. 0.005) then
        singloid =              dt * 1000.0 + 0.2
        dt       = DBLE (INT (singloid         )) / 1000.0
****    dt       = FLOAT (INT (dt * 1000.0 + 0.2)) / 1000.0
      end if
      dtmax  = dt
      idterr = 0
c
      write  (nqik, 9010) dtmax
      write  (nout, 9010) dtmax
      write  (ncrt, 9010) dtmax
 9010 format (' dtmax limited to ', 1pe12.5, ' s')
c
c ----------------------------------------------------------------------
c time step has been reduced below dtmin
c ----------------------------------------------------------------------
c
 5010 if (dt .gt. dtmin)  go to 3000  ! try again with smaller time step
      write  (ncrt, 8200)
      write  (nout, 8200)
      write  (nqik, 8200)
 8200 format (' time step has been reduced below dtmin' /
     .        ' ONETWO execution will be halted')
      istop = 1
      time_step_too_small = .true.
c
      write (nmhddat, '("leaving TPORT,dt too small, time = ",
     .                                                  1pe14.8)') time
      write (nmhddat, '("rbp(2-5) = "   ,4(2x,1pe14.6))')(rbp (j),j=2,5)
      write (nmhddat, '("curden(1-4) = ",4(2x,1pe14.6))')
     .                                                 (curden(j),j=1,4)
      write (nmhddat, '("hcap(1-4) = "  ,4(2x,1pe14.6))')(hcap(j),j=1,4)
      write (nmhddat, '("dhdt(1-4) = "  ,4(2x,1pe14.6))')(dhdt(j),j=1,4)
      write (nmhddat, '("gcap(1-4) = "  ,4(2x,1pe14.6))')(gcap(j),j=1,4)
      write (nmhddat, '("dgdt(1-4) = "  ,4(2x,1pe14.6))')(dgdt(j),j=1,4)
      write (nmhddat, '("fcap(1-4) = "  ,4(2x,1pe14.6))')(fcap(j),j=1,4)
      write (nmhddat, '("dfdt(1-4) = "  ,4(2x,1pe14.6))')(dfdt(j),j=1,4)
      write (nmhddat, '("r(2-5) = "     ,4(2x,1pe14.6))')(r   (j),j=2,5)
      write (nmhddat, '("q(1-4) = "   ,4(2x,1pe14.6))')  (q   (j),j=1,4)
      write (nmhddat, '("ub(1-4) = "   ,4(2x,1pe14.6))') (ub  (j),j=1,4)
      write (nmhddat, '("ub(3-5) = "   ,4(2x,1pe14.6))') (ub (j),j=3,5)
      write (nmhddat, '("psir(1-4) = "   ,4(2x,1pe14.6))') 
     .                                                (psir   (j),j=1,4)
      write (nmhddat, '("Rbp(nj-3,nj) = "   ,4(2x,1pe14.6))') 
     .                                             (rbp   (j),j=nj-3,nj)
c
      call write_restart_profs

      return
c
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c take care of final time
c either we are finished or else it's time for a new equilibrium calculation
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c ----------------------------------------------------------------------
c generate time-average information if timav .ne. 0
c ----------------------------------------------------------------------
c
 6000 continue

      if (timav .eq. 0.0)  go to 6100
      if (dtsum .eq. 0.0)  go to 6100
c      write(940,FMT='("c306,line 2883 ,totrf =",1pe12.4)')totrf ! 888888889999
      if (timav .lt. 0.0)  call out

      dtsumi = 1.0 / dtsum
      do 6020 j=1,nj
      do 6010 k=1,nk
 6010 uav(k,j) = uav(k,j)*dtsumi
      do 6015 ib=1,nbeams
      do 6015 ie=1,3
      enbav(j,ie,ib)  = enbav(j,ie,ib)*dtsumi
      pprbav(j,ie,ib) = pprbav(j,ie,ib)*dtsumi
      wbav(j,ie,ib)   = wbav(j,ie,ib)*dtsumi
 6015 ppbav(j,ie,ib)  = ppbav(j,ie,ib)*dtsumi
      enaav(j)        = enaav(j)*dtsumi
 6020 waav(j)         = waav(j)*dtsumi
      call copya  (ene,eneav1,nj)
      call copya  (enbeam,enbav1,nj)
      call update (uav,en,te,ti,rbp,nk,nj,kj,kk,iangrot,angrot)
      call updatb (enbav, kj, ke, nj, nbeams, one, enbeam)
      call zen
c
c     copy uav into u and usave, and save present u in uav
c
      call savesol  (uav,usave,nk,nj,kk)
      call savesol  (u,uav,nk,nj,kk)
      call savesol  (usave,u,nk,nj,kk)
      call copya (ene,enesav,nj)
      do 6030 k=1,nimp
 6030 call zeroa (dzdtim(1,k),nj)




      call dndt  (dtsumi,eneav1,eneav0,kj,kk,nj,nion,uav,uav0,
     .            dnedt,dnidt,dtedt,dtidt,dpedt,dpidt,iangrot,
     .            nk,dangrot)
      do 6040 ib=1,nbeams
      do 6040 ie=1,3
      call copya (enbav(1,ie,ib),enbsav(1,ie,ib),nj)
      call copya (wbav(1,ie,ib),wbsav(1,ie,ib),nj)
      call copya (pprbav(1,ie,ib),pprbsav(1,ie,ib),nj)
 6040 call copya (ppbav(1,ie,ib),ppbsav(1,ie,ib),nj)
      call copya (enaav,enasav,nj)
      call copya (enalp,enaav,nj)
      call copya (waav,wasav,nj)
c
      counter = counter + 1 ! count how many times DIFFUS is called
****  if (counter .gt. 59)  call STOP ('subroutine TPORT: DEBUG', 904)!
      cparam =1.0d0
      IF(diffeq_methd .ne. 3)call diffus (xi_include)
      dt     = 0.0
      dtt    = 0.0
      itimav = 1
      istep  = 'pred'
      print *,' calling source before returning to runtwo'


      call source
      call fluxx
      dt     = 1.0e-10
      dtt    = 1.0e-10
      delav  = 0.0
      delmax = 0.0
      kmax   = 0
      jmax   = 0
      delit  = 0.0
      iter   = 0
      imix   = 0
      iconvg = 1
      call info (nout, nqik, ihead, nnew, extime, timnew, dt, dtt,
     .           delav, delmax, kmax, jmax, delit, iter, ineu, inub,
     .           irfcalc, imix, iconvg, isecrem0)
      call impsrc
      nsave = n
      n     = 0
      call fiziks
c
c get the experimental profiles if this is tdem mode
c
c ----------------------------------------------------------------
c
c
      n = nsave
c
c ----------------------------------------------------------------------
c call output routines at final time point
c ----------------------------------------------------------------------
c
 6100  continue
       if (impfh) then
        delcapc     =  0.0
        implicit_fh = .true.
        call fhcalc (0, 1)
        do j=1,nj
          dhcap   = ABS ((hcap(j)-hcap0(j))/hcap(j))
          dfcap   = ABS ((fcap(j)-fcap0(j))/fcap(j))
          delcapc = MAX (delcapc, dfcap, dhcap)
        end do
        delcapc   = 100.0 * delcapc
        write  (ncrt, 8099) delcapc
 8099   format (' maximum percent change in fcap, hcap = ', f10.2)
      end if
      ilastp = 1


      call out

      kine_message  = 'final'
      kine_stat = 1
      call wrt_kin_efit_naml(1,kine_stat)

      if (tportvb .ge. 1)  write (*, '(" calling   TRPLOTC at time = ",
     .                                   1pe14.6)') time
      call trplot (3)
      if (tportvb .ge. 1)  write (*, '(" done with TRPLOT")')

      if (inubplt .eq. 1)  call nubplt

      if (n .ge. nmax) then
        istop = 1
        write  (ncrt, 8201)  n, nmax
        write  (nout, 8201)  n, nmax
        write  (nqik, 8201)  n, nmax
 8201   format (' *************  n >= nmax, RUN WILL BE STOPPED',
     .         '  *************' /
     .          ' n, nmax =', 2(2x, i5))
      end if
      if (time .ge. timmax)  istop = 1
c
      write (nmhddat, '("leaving TPORT,time .ge. timmax, time = ",
     .                                                   1pe14.8)') time
      write (nmhddat, '("rbp(2-5) = "   ,4(2x,1pe14.6))')(rbp (j),j=2,5)
      write (nmhddat, '("curden(1-4) = ",4(2x,1pe14.6))')
     .                                                 (curden(j),j=1,4)
      write (nmhddat, '("hcap(1-4) = "  ,4(2x,1pe14.6))')(hcap(j),j=1,4)
      write (nmhddat, '("dhdt(1-4) = "  ,4(2x,1pe14.6))')(dhdt(j),j=1,4)
      write (nmhddat, '("gcap(1-4) = "  ,4(2x,1pe14.6))')(gcap(j),j=1,4)
      write (nmhddat, '("dgdt(1-4) = "  ,4(2x,1pe14.6))')(dgdt(j),j=1,4)
      write (nmhddat, '("fcap(1-4) = "  ,4(2x,1pe14.6))')(fcap(j),j=1,4)
      write (nmhddat, '("dfdt(1-4) = "  ,4(2x,1pe14.6))')(dfdt(j),j=1,4)
      write (nmhddat, '("r(2-5) = "     ,4(2x,1pe14.6))')(r   (j),j=2,5)
      write (nmhddat, '("q(1-4) = "   ,4(2x,1pe14.6))')  (q   (j),j=1,4)
      write (nmhddat, '("ub(1-4) = "   ,4(2x,1pe14.6))') (ub  (j),j=1,4)
      write (nmhddat, '("ub(3-5) = "   ,4(2x,1pe14.6))') (ub (j),j=3,5)
      write (nmhddat, '("Rbp(nj-3-nj) = "   ,4(2x,1pe14.6))') 
     .                                             (rbp   (j),j=nj-3,nj)

 
      call write_restart_profs


c

 
      return
c
      end





      subroutine xptor_init
c----------------------------------------------------------------------------
c     a test section.  load usave,u ,te,ti,etc with profiles
c     read from  (usually xptor generated results ) file.
       USE param
       USE soln
       USE xptor_sim
       USE numbrs
      USE mesh
      USE tordlrot
       implicit  integer (i-n), real*8 (a-h, o-z)
c       include 'numbrs.i'  !nj
c       include  'mesh.i'
c       include 'tordlrot.i'  !angrot,iangrot
       logical no_eof,exists
       character *128 input_string

       iout = 91
       allocate (rxpt(kj),STAT = istat)
       allocate (texpt(kj),STAT = istat)
       allocate (tixpt(kj),STAT = istat)
       allocate (qbeame_fixed(kj),STAT = istat)
       allocate (qbeami_fixed(kj),STAT = istat)
       allocate (qrfe_fixed(kj),STAT = istat)
       allocate (qrfi_fixed(kj),STAT = istat)



       read_qbeami = 0
       read_qbeami = 0
       read_qrfe = 0
       read_qrfi = 0




       call getioun(iout,iout)
       exists = .false.
       INQUIRE(FILE = 'profile.dat', EXIST = exists)
!       if(profile.dat doesnt exist in current working directory then
!      assume the user does not want to initalize te,ti using this file.
!      Instead use the initalizations in inone as usuall:
       if(exists)then
           open (unit = iout, file = 'profile.dat', status = 'OLD', 
     .                                 iostat = iostat)
           do j =1,nj
             read(iout,*,err = 2,end= 3)rxpt(j),texpt(j),tixpt(j)
           enddo
           close(unit = iout)
           call giveupus(iout)
           go to 10

 2          call stop('xptor_init error reading profile.dat',1)
 3          call stop('xptor_init error reading profile.dat',2)
 10         continue

            do j=1,nj
               diff = abs(r(j)/r(nj) -rxpt(j))
               if (diff .gt. 1.e-4)
     .          call stop('xptor_init, r grid not comensurate',1)
           
            enddo
            !copy texpt,etc into u:
            call redate (u,en,texpt,tixpt,rbp,nk,nj,kj,kk,
     .                                    iangrot,angrot)
            !copy  u into usave:

            call savesol (u, usave, nk, nj, kk)
        
            !copy texpt to te,etc:
            call copya (texpt, te, nj)
            call copya (tixpt, ti, nj)
        endif
c


c      get qbeame,qbeami,qrfe,qrfi form iterdb type file.
c      dont use the iterdb read subroutine  because we dont
c      want to store directly into the common blocks:
c      iterdb_xptor should be an iterdb file  created by Onetwo
c      and used by xptor in its calculations:
       exists = .false.
       INQUIRE(FILE = 'iterdb_xptor',EXIST  = exists)
c       if(exists)then
           call getioun(iout,iout)
           open (unit = iout, file = 'iterdb_xptor', status = 'OLD', 
     .                                 iostat = iostat)

       
           !loop over lines in file:
           no_eof = .true.
           do while(no_eof)
             read(iout,FMT ='(a)')input_string
             if(input_string   .eq.
     .             "*  power to elec. from beam, watts/meter**3")
     .        then
                read(iout, 15,end = 20) (qbeame_fixed(j),j=1,nj)
                read_qbeame = 1
             endif
             if(input_string .eq.
     .             '*  power to ions from beam, watts/meter**3' )
     .        then
                read(iout, 15,end = 20) (qbeami_fixed(j), j=1,nj)
               read_qbeami = 1
             endif




!             RF power:

             if(input_string   .eq.
     .             '*  qrfe, RF electron heating, watts/meter**3')
     .        then
                read(iout, 15,end = 40) (qrfe_fixed(j),j=1,nj)
                read_qrfe = 1
             endif

             if(input_string   .eq.
     .            '*  qrfi, -RF ion heating, watts/meter**3')
     .        then
                read(iout, 15,end = 40) (qrfi_fixed(j),j=1,nj)
                read_qrfi = 1
             endif

             if(read_qbeami + read_qbeame 
     .          + read_qrfe + read_qrfi .eq. 4) no_eof = .false.


           enddo  !do while no_eof


c           convert from w/m^3  to (kev/(cm^3 sec) for use in sub source:
            qbeame_fixed(:) = qbeame_fixed(:)*6.2415097e+09
            qbeami_fixed(:) = qbeami_fixed(:)*6.2415097e+09
            qrfe_fixed(:) = qrfe_fixed(:)*6.2415097e+09
            qrfi_fixed(:) = qrfi_fixed(:)*6.2415097e+09


            close(unit = iout)
            call giveupus(iout)

 20         if(read_qbeami .eq. 0 .or. read_qbeame .eq. 0)
     .      call STOP('xptor_init file iterdb_xptor,no beam dazta',1)
 40         if(read_qrfe  .eq. 0 .or. read_qrfi .eq. 0)
     .      call STOP('xptor_init file iterdb_xptor,no rf data',1)
c        endif  ! iterdb_xptor exists 
        deallocate(rxpt)
        deallocate(texpt)
        deallocate(tixpt)
        
        return



 15     format (5(2x,1pe14.4))  ! common floating  write/read format
        end
