
      subroutine central_value (valprev, valcur, theta, timecur, dt,
     .                          valthet)
c
c ----------------------------------------------HSJ--9/11/95------------
c
c     this subroutine calculates the value at the "n+theta" time
c     given the value at n-1 and at n
c
c     INPUT
c        valprev           The value at time point n-1
c        valcur            The value at time point n
c        timecur           The time at time point n
c        dt                dt =time(n)-time(n-1)
c        theta             The pivot point
c
c     OUTPUT
c        valthet           value at time(n-1)+theta*dt
c
c ----------------------------------------------------------------------
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      valthet = theta*valcur+(1.-theta)*valprev
      return
c
      end



      subroutine diffus (xi_include)
      USE param
      USE constnts,ONLY : pi,rootpi,root2,cee,xmasse,xmassp,charge,
     .                    charg4,rot2pi, rotpi2
      USE fusion
      USE aid_newton
      USE io 
      USE ions
      USE nub
      USE nub2
      USE solcon
      USE solcon_gcnmp,              ONLY : use_stab_flux
      USE soln
      USE mhdpar
      USE extra
      USE numbrs
      USE mesh
      USE paleocl ,only : include_paleo,paleocl_tr
c
      USE sourc
      USE machin
      USE tfact
      USE geom
      USE flags
      USE tordlrot
      USE tcoef
      USE bd_condtn
      USE rhog
      USE ifs
      USE staebler
      USE flx
      USE neo2d
      USE kapisn_module,               ONLY : kapisn_driver
      USE mmm95_mod !JMP
      USE brainfuse_mod !OM
      USE glf23 !jmp.snu
#ifdef PARALLEL_12
      USE mpi12 !jmp.snu.par
#endif
      USE weiland
      USE island
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      character rcs_id*63
      save      rcs_id
      data      rcs_id /
     ."$Id: cray307.F,v 1.37 2014/05/08 17:48:11 stjohn Exp $"/
c
c     This subroutine calculates transport (or diffusion) coefficients.
c     Note that the logic in TPORT is such that
c     DIFFUS, FLUXX and SOURCE are always called at the time stored
c     in variable time (in INCLUDE file solcon.i)
c



      include 'mod_gbohm.i'

      include 'rebut.i'

      include 'shay.i'
      include 'storage.i'

     

c
      dimension n_mode(10),m_mode(10)
      dimension rdels(kj),qdels(kj),tieff(kj)
      dimension dtest(kj),dtestt(kj),xktest(kj),colfac(kj,2)
      dimension xkbtst(kj),alpha_exp_save(kj)
      dimension dd(kj),dsav(kk,kk,kjm1),qexsav(kj)
****  dimension modben(2), slen(kion)
      dimension modben(2),dwifs(kj)
      dimension rzsqa(kion),den(kion)
      dimension ena(kion), dzdtea(kion), za(kion), zsqa(kion),
     .          di(kion), rho(kion), xnu(kion), vth(kion)
      dimension xk0(3,3), xk01(3,3), xk0i(3,3),
     .          a(3,3), a1(3,3), a4(3,3), b(3,3), b1(3,3), b4(3,3),
     .          c(3,3), c1(3,3), c4(3,3)
      dimension xk(3,3), xia(3,3)
      dimension tot_primary_ion_den(kj),tot_thermal_ion_den(kj)
c
c --- the following is only for use in debugging subroutine DIFFUS
c
      dimension  xiaj(kj), xkj(kj), xk0j(kj), aj(kj), bj(kj), cj(kj),
     .           xnuej(kj), rdeltaj(kj), xhm2aj(kj), xlamj(kj),
     .           safacj(kj), vthej(kj), vthj(kj), rhoj(kj), xnuj(kj),
     .           bpaj(kj), zinj(kj), rhoej(kj), xnusij(kj), del32j(kj),
     .           eneaj(kj), teaj(kj)
c
c --- the following for Hirshman bootstrap calculation
c
      dimension zstar(kion+1),xra(kion+1),xnush(kion+1),xlogl(kion+1)
      dimension xmassa(kion+1),xta(kion+1),xna(kion+1),zeda(kion+1)
      dimension xu1(kion+1),xu2(kion+1),xu1s(kion+1),xu2s(kion+1)
      dimension bp(kj),aabsq(kj),glf_work(kj)
      real *8  xketyp_prev(kj),xketotm,xkitotm,xkangtotm
      equivalence (tieff(1),rdum(1))
****  equivalence (aabsq(1),udum(1))
c
      logical xi_include
c
      integer diffus_counter    , ip_chi2_counter
      data    diffus_counter /0/, ip_chi2_counter /0/
c
      data  puny   / 1.0e-6      /
      data xk01/1.04,0.  ,0.  ,1.20,2.55,0.  ,2.30,4.19,1.83/,
     .     xk0i/0.73,0.  ,0.  ,0.73,1.46,0.  ,1.46,3.65,1.46/,
     .     a1  /2.01,0.  ,0.  ,0.76,0.45,0.  ,1.02,0.57,0.68/,
     .     a4  /2.30,0.  ,0.  ,0.80,0.46,0.  ,0.79,0.48,0.47/,
     .     b1  /1.53,0.  ,0.  ,0.67,0.34,0.  ,0.75,0.38,0.32/,
     .     b4  /0.98,0.  ,0.  ,0.42,0.22,0.  ,0.56,0.33,0.20/,
     .     c1  /0.89,0.  ,0.  ,0.56,0.43,0.  ,1.07,0.61,0.66/,
     .     c4  /0.74,0.  ,0.  ,0.48,0.30,0.  ,0.51,0.28,0.51/
      data xk20,a2,b2,c2/0.66,1.03,0.29,0.74/
c
      f1i(zin, x1, xidum) = xidum + zin * (x1-xidum)
      f14(zin, x1, x4) = 1.33333*(1.0-zin)*x4 +0.333333*(4.0*zin-1.0)*x1
c

c
cjek - glf23 variables
c
      integer  jshoot, jmm, jmaxm, igrad,
     .        idengrad_glf
c
      real*8 rmajor_glf, rminor_glf, btor_exp, 
     .   zimp_glf, amassimp_glf, arho_iglf,
     .   alphad_glf, alphae_glf, 
     .   etaphim, etaparm, etaperm, exchm,
     .   fc, akappa1, akappa2, alpha_neo, vstar_sign, corot,
     .   pgeo_local, rdrho_local, rdrho_local_p1,zval

c
cJMP  mmm variables
c
      integer matdim_mm, npoints_mm, nprout_mm 
      parameter (matdim_mm=12)
      integer lprint_mm, nerr_mm , lsuper_mm, lreset_mm, lswitch_mm(5)
      real*8  cswitch_mm(23), fig_mm(4), frb_mm(4), fkb_mm(4)

      real*8 rminor_mm, rmajor_mm, elong_mm
      real*8 dense_mm, densh_mm, densimp_mm, densfe_mm, denstot_mm
      real*8 densh_m2, densimp_m2, denstot_m2
      real*8 xzeff_mm, tekev_mm, tikev_mm, q_mm, btor_mm

      real*8 avezimp_mm, amassimp_mm, amasshyd_mm, aimass_mm
      real*8 wexbs_mm

      real*8 grdne_mm, grdni_mm, grdnh_mm, grdnz_mm
      real*8 grdte_mm, grdti_mm, grdq_mm

      real*8 thiig_mm, thdig_mm, theig_mm, ythzig_mm
      real*8 thirb_mm, thdrb_mm, therb_mm, ythzrb_mm
      real*8 thikb_mm, thdkb_mm, thekb_mm, ythzkb_mm

      real*8 gamma_mm(matdim_mm), omega_mm(matdim_mm)
      real*8 difthi_mm(matdim_mm,matdim_mm)
      real*8 velthi_mm(matdim_mm), vflux_mm(matdim_mm)

      real*8 rda_mm,dens_min_mm

      real*8 chi_e_ifs_save(kj)
      real*8 chi_i_ifs_save(kj)
      real*8 xkang_ifs_save(kj)
      real*8 d_ifs_save(kj)
      real*8 diff_ms(kj) !jmp.den
      real*8 swork(kj),pdum(kj) !jmp.den
      integer icalled
      save icalled
      data icalled /0/

      !jmp.ibm.par start
      integer jin,jout  
      integer IPARAM(32)
      real*8  RPARAM(32)
      integer iflag_mpistop(1)
      !jmp.ibm.par end


c
c ----------------------------------------------------------------------
c --- debug print parameters. turn on and off with ddt
c
      data idbgprt, jdbgprt, iodbgprt, ibrkpt, ldbglnth
     .    /      0,       0,       87,      0,   200000/
c
****  if (ip_chi2_counter .ge. 2950)  call STOP ('DIFFUS: debug', 991)!
      diffus_counter = diffus_counter + 1
c
      zero          = 0.0D0  !jmp.ibm
      one           = 1.0D0  !jmp.ibm
      one_hundredth = 0.01D0 !jmp.ibm





c
c --- create a breakpoint label
c
      if (ibrkpt .ne. 0)  ibrkpt = ibrkpt + 1
c
      if (idbgprt .eq. 1 .and. jdbgprt .eq. 0)then
        call getioun(iodbgprt,iodbgprt)
        open (unit = iodbgprt, file = 'dbgprint', status = 'UNKNOWN')
      endif
c
c ----------------------------------------------------------------------
c calculate safety factor ; use quadratic extrapolation on axis
c bp(j) is local to this routine
c ----------------------------------------------------------------------
c
      j_clamp = -1
      do 10 j=2,nj
         bp(j) = rbp(j)/(fcap(j) * gcap(j) * hcap(j) * r(j))
         if(r(j)/r(nj) .le. r_clamp)j_clamp = j
   10    q(j) = r(j) * btor / (rmajor * bp(j))
      bp(1) = 0.0
       q(1) = (q(2)*r(3)**2-q(3)*r(2)**2)/(r(3)**2-r(2)**2)
c
c     keep  ABS (q(j)) < 1000
c
      do j=1,nj
        absq = ABS (q(j))
        if (absq .ge. 1.0e3)  q(j) = 1.0e3*q(j) / absq
        aabsq(j) = ABS (q(j))
      end do
c
c ----------------------------------------------------------------------
c                        island-induced transport
c ----------------------------------------------------------------------
c set up and call 2/1 and 3/2 island width calculation;
c     model is valid only for circular geometry
c ----------------------------------------------------------------------
c
      nnisl   = nisl
      if (wisl .le. 0.0D0)  go to 50
c
c    A vacuum boundary with njb grid points can be added to DELSOL.
c    It corresponds to putting conducting wall farther away.
c    Note nj+njb must be less than kj
c
      njb   = ddebug(16)
      ndels = nj + njb
      dq0   = 0.1
      dlam  = 0.25
      nnq0  = 1
      nnlam = 1
c
      if (nisl .lt. -10) then
        nnq0  =  6
        nnlam = 17
      end if
c
      do 51 nq0=1,nnq0
      anq0 = nq0
      do 51 nlam=1,nnlam
      anlam = nlam
c
      do jj=1,nj
        rdels(jj) =      r(jj)
        qdels(jj) =      q(jj)
        xlg       = ddebug(24)
c
c       test ORNL q-profile
c
        if (ddebug(21) .ne. 0.0) then
          qaornl    = ddebug(21)
          q0ornl    = ddebug(22) +  dq0 * ( anq0 - 1.0)
          xlornl    = ddebug(23) + dlam * (anlam - 1.0)
          r0ornl    = r(nj)
     .              * ((qaornl/q0ornl)**xlornl-1.0)**(-1.0/2.0/xlornl)
          qdels(jj) = q0ornl
     .              * (1.0+(r(jj)/r0ornl)**(2.0*xlornl))**(1.0/xlornl)
        end if
      end do
c
      delsdr = r(nj) - r(nj-1)
      if (njb .le. 0)  go to 49
      do jb=1,njb
        jj        = nj + jb
        ajb       = jb
        rdels(jj) = r(nj)+delsdr*ajb
        qdels(jj) = qdels(nj)*rdels(jj)**2/r(nj)**2
      end do
c
   49 n_mode( 1) = 1
      n_mode( 2) = 1
      n_mode( 3) = 1
      n_mode( 4) = 1
      n_mode( 5) = 1
      n_mode( 6) = 2
      n_mode( 7) = 2
      n_mode( 8) = 2
      n_mode( 9) = 2
      n_mode(10) = 2
c
      m_mode( 1) = 2
      m_mode( 2) = 3
      m_mode( 3) = 4
      m_mode( 4) = 5
      m_mode( 5) = 6
      m_mode( 6) = 3
      m_mode( 7) = 5
      m_mode( 8) = 7
      m_mode( 9) = 9
      m_mode(10) =11
c
      if (time .lt. timisl)  go to 45
      dtisl   = time-timisl
      timeisl = time
      nnisl   = nisl
      if (nisl .le. 0)  nnisl = -nisl
      if (nisl .eq. 0)  nnisl =  2
      do jisl=1,nnisl
        modben(1) = m_mode(jisl)
        modben(2) = n_mode(jisl)
        if (jisl .eq. 1 .and. nisl .eq. 0)  modben(1) = m_mode(1)
        if (jisl .eq. 1 .and. nisl .eq. 0)  modben(2) = n_mode(2)
        if (jisl .eq. 2 .and. nisl .eq. 0)  modben(1) = m_mode(6)
        if (jisl .eq. 2 .and. nisl .eq. 0)  modben(2) = n_mode(6)
        wdisl0 = wdisl(jisl)
        call delsol (ndels,rdels,qdels,modben,wmn,delpmn,rsmn,btlmn,xlg,
     .               wdisl0,dpisl0)
        js     = rsmn/(r(2)-r(1)) + 1
        ratisl = 1.66 * cee**2 * 0.51 * 1.15e-14 * 15.0 /
     .                       (te(js) * 1000.0)**1.5 / 4.0 / pi
        wdisl(jisl) = 0.3*(wdisl0+ratisl*dpisl0*dtisl)+(1.0-0.3)*wdisl0
        if (wdisl(jisl) .lt. 0.0)  wdisl(jisl) = 0.0
         risl(jisl) = rsmn
        dpisl(jisl) = dpisl0
        if (nisl .lt. 0)  wdisl(jisl) = wmn
        if (nisl .lt. 0)  dpisl(jisl) = delpmn
      end do
c
   45 if (nisl .ge. -10)  go to 51
      overlap = (risl(1)-wdisl(1)/2.0)-(risl(2)+wdisl(2)/2.0)
      write  (nout, 9000)  q0ornl, xlornl, overlap, risl(1), wdisl(1),
     .                                              risl(2), wdisl(2)
      write  (ncrt, 9000)  q0ornl, xlornl, overlap, risl(1), wdisl(1),
     .                                              risl(2), wdisl(2)
 9000 format (1p7e10.3)
   51 continue
c
      if (nisl .lt. -10)  nisl = -10
c
c  ion cyclotron frequency
c
   50 xmasip = atw(1)*xmassp
      omci   = charge*btor/(xmasip*cee)
c
c ----------------------------------------------------------------------
c initialize transport coefficients
c ----------------------------------------------------------------------
c
      wneo1 = wneo(1,1)+wneo(1,2)+wneo(1,3)+wneo(1,4)+wneo(1,5)
      wneo2 = wneo(2,1)+wneo(2,2)+wneo(2,3)+wneo(2,4)+wneo(2,5)
      wneo3 = wneo(3,1)+wneo(3,2)+wneo(3,3)+wneo(3,4)+wneo(3,5)
      wneo5 = wneo(5,1)+wneo(5,2)+wneo(5,3)+wneo(5,4)+wneo(5,5)
      const = 4.0 * pi * pi*rmajor*1.60217733e-16
c
c --- integrate the power if it is used in take your pick transport
c
        if ((typa(16) .ne. 0.0 .and.     wtyp .ne. 0.0)  .and.
     .      (typa( 9) .ne. 0.0 .or.  typa(10) .ne. 0.0)) then
          call trapv (r, qbeame, hcap, nj, pbeame)
          call trapv (r, qbeami, hcap, nj, pbeami)
          call trapv (r, qrfe  , hcap, nj, prfe  )
          call trapv (r, qrfi  , hcap, nj, prfi  )
c
c --- convert total power to watts
c
          pbeames = pbeame * const
          pbeamis = pbeami * const
          prfes   =   prfe * const
          prfis   =   prfi * const
c
c --- pohm is calculated in subroutine SOURCE (in units of watts)
c
          pohms   = pohm
        end if
      qexsav(nj) = qexch(nj)
c
c --- check for flow shear turblence suppression
c --- if it is turned on, load the multipliers into common in staebler.i
c --- use these multipliers below to modify diffusivities appropriately
c
      ifsflag = INT (fs(1))  ! =1 if flow shear turblence suppress is on
      do j=1,nj-1
        coefa(j) = 1.0
        coefb(j) = 1.0
        coefc(j) = 1.0
        coefd(j) = 1.0
      end do
c
c     staebler/hinton flow shear suppression factor
c     (modfies coefa, coefb, coefc, coefd):
c
      if (ifsflag .gt. 0)  call cstaebler
c
c     DeBoo (phenomenological) flow shear suppression factor (modifies coefa,
c     coefb, acts in place of or in addition to the staebler/hinton model,
c     depending on which models are selected):
c
      if (irgc .gt. 0)  call cdeboo (coefa, coefb)
      tstaebler = time
c
c ----------------------------------------------------------------------
c       start loop over mesh points on half grid
c ----------------------------------------------------------------------
c
      do 2000 j=1,nj-1
           do 60 k=1,nk
           do 60 i=1,nk
                dsav(i,k,j) = d(i,k,j)
   60           d(i,k,j)  = 0.0
      qexsav(j) = qexch(j)
      qieneo(j) = 0.0
      dneo(j)   = 0.0
      xkeneo(j) = 0.0
      xkineo(j) = 0.0
      chiineo(j)= 0.0
      dtyp(j)   = 0.0
      xketyp(j) = 0.0
      xkityp(j) = 0.0
      dsaw(j)   = 0.0
      xkesaw(j) = 0.0
      xkisaw(j) = 0.0
      xkeisl(j) = 0.0
      xkiisl(j) = 0.0
      disl(j)   = 0.0
      xkepp(j)  = 0.0
      xkeohk(j) = 0.0
      xkicha(j) = 0.0
      xkecar(j) = 0.0
      xkedom(j) = 0.0
c
c     calculate dependent variables at mesh centers
c
      enea    = 0.5 * ( ene(j)+ ene(j+1))
      tea     = 0.5 * (  te(j)+  te(j+1))
      tia     = 0.5 * (  ti(j)+  ti(j+1))
      curdna  = 0.5 * (curden(j)+curden(j+1))
      dte     = (te(j+1)- te(j))/dr(j)
      dti     = (ti(j+1)- ti(j))/dr(j)
      zeffa   = 0.5 * (zeff(j)+zeff(j+1))
      zin     = 1.0 / zeffa
      angrota = 0.5 * (angrot(j)+angrot(j+1))
c
c calculate some more quantities at mesh centers
c
      enasum = 0.0
      do k=1,nion
        za    (k)  = 0.5 * ( z(j,k)+ z(j+1,k))
        zsqa  (k)  = 0.5 * (zsq(j,k)+zsq(j+1,k))
        rzsqa (k)  = SQRT (zsqa(k))
        dzdtea(k)  = 0.5 * (dzdte(j,k)+dzdte(j+1,k))
        ena(k)     = 0.5 * (en(j,k)+en(j+1,k))
        enasum = enasum+ena(k) ! thermal dif should include all ions HSJ
        den(k)     = (en(j+1,k)-en(j,k))/dr(j)
        dudr(k,j)  = den(k)
        etaim(j,k) = 0.0
        dentia     = tia*(ene(j+1)-ene(j))/dr(j)
        if (dentia .ne. 0.0 .and. k .le. nprim)
     .    etaim(j,k) = enea * dti / dentia
      end do
c
      k           = nion + 1
      dudr(k,j)   = dte
      k           = nion + 2
      dudr(k,j)   = dti
      k           = nion + 3
      drbp        = (rbp(j+1)-rbp(j))/dr(j)
      dudr(k,j)   = drbp
      vionzgrd(j) = 0.0
      if (iangrot .eq. 0)  go to 85
      k           = nion+4
      drangrot    = (angrot(j+1)-angrot(j))/dr(j)
      dudr(k,j)   = drangrot
      vionzgrd(j) = (angrot(j+1)*r2capi(j+1)/rcap(j+1)
     .              -angrot(j)*r2capi(j)/rcap(j))/dr(j)
c
c  obtain delta and flux surface integrals
c
   85 if (codeid .ne. 'onedee')  go to 90
      delta    = ra(j) / (SQRT (kappa) * rmajor)
      rdelta   = SQRT (delta)
      delta2   = delta**2
      xhm2a    = 1.0
      xia(1,1) = 1.95*rdelta
      xia(3,3) = 1.95*rdelta
      xipsa    = 2.00*delta2
      go to 95
c
c --- eps,xhm2,xi11,xi33,xips are updated in rhomesh if xi_include = true
c
   90 if (.not. xi_include) then
        call neointrp (eps0, eps, j, epsa)
      else
        epsa = (eps(j+1)+eps(j)) * 0.5
      end if
      delta  = epsa
      rdelta = SQRT (delta)
      delta2 = delta**2
      if (.not. xi_include) then
        call neointrp (xhm20, xhm2, j, xhm2a   )
        call neointrp (xi110, xi11, j, xia(1,1))
        call neointrp (xi330, xi33, j, xia(3,3))
        call neointrp (xips0, xips, j, xipsa   )
      else
        xhm2a    = 0.5 * (xhm2(j+1)+xhm2(j))
        xia(1,1) = 0.5 * (xi11(j+1)+xi11(j))
        xia(3,3) = 0.5 * (xi33(j+1)+xi33(j))
        xipsa    = 0.5 * (xips(j+1)+xips(j))
      end if
   95 xia(1,2) = xia(1,1)
      xia(1,3) = xia(3,3) / xhm2a
      xia(2,2) = xia(1,1)
      xia(2,3) = xia(1,3)
c
c --- calculate trapped/circulating fraction HSJ
c
****  if (jhirsh .eq. 1 .or. jhirsh .eq. 2 ) then
      if (ftcalc .eq. 'analytic') then
c
c --- see for example Hirshman et al N.F. 17,3 (1977) pg 611, eq. 8, for
c --- the following approximate form of xft
c
          xft = 1.0 - (1.0 - delta)**2 / SQRT (1.0 - delta2)/
     .                                   (1.0 + 1.46*rdelta)
      else    ! exact calculation (when jhirsh .ne. 1 or 2)
          if (codeid .eq. 'onedee') then
            xft = 0.75*xia(3,3)/xhm2a
          else
            xft = 0.75*(xi33(j+1)+xi33(j))/(xhm2(j+1)+xhm2(j))
          end if
      end if
      xfc     = 1.0 - xft
      ftfc(j) = 0.0
      if (xfc .gt. 0.0)  ftfc(j) = xft/xfc
c
c  Save trapped fraction for Houlberg Boostrap Model - DFF
c
      ftncl(j) = xft
c
c  obtain geometric factors
c
      del32  = delta*rdelta
      del32i = 1.0 / del32
      del52  = delta*del32
      gcapa  = 0.5 * (gcap(j)+gcap(j+1))
      fcapa  = 0.5 * (fcap(j)+fcap(j+1))
      hcapa  = 0.5 * (hcap(j)+hcap(j+1))
      r2capa = 0.5 * (r2cap(j)+r2cap(j+1))
      grsq   = gcapa/r2capa
c
c  calculate bp, safac and shear parameter
c
      bp1       = rbp(j)/((r(j)+puny)*fcap(j)*gcap(j)*hcap(j))
      bp2       = rbp(j+1)/(r(j+1)*fcap(j+1)*gcap(j+1)*hcap(j+1))
      bpa       = 0.5 * (bp1+bp2)
      rbpa      = ra(j)*bpa
      safac     = ra(j)*btor/(rmajor*bpa)
      safac     = ABS (safac)
      shearp(j) = ra(j) * ABS (q(j+1)-q(j))/(safac*dr(j))
      shearp(j) = MAX (shearp(j), one_hundredth)
c
c     calculate Larmor radii and various collision frequencies
c
      teaerg = 1.6e-9 * tea
      vthe   = SQRT (2.0*teaerg/xmasse)
      rhoe   = xmasse*cee*vthe/(charge*bpa)
      xlam   = 24.0 - LOG (SQRT (enea)/(1.0e3*tea))
      xnue   = 1.33333*rot2pi*enea*zeffa*charg4*xlam
     .                   / SQRT (xmasse*teaerg**3)
      xnuse(j) = root2*rmajor*safac*xnue/(del32*vthe)
      tiaerg   = 1.6e-9*tia
      xnubyv   = 1.33333*rotpi2*enea*zeffa**3*charg4*xlam/tiaerg**2
      xnusi(j) = root2*rmajor*safac*xnubyv/del32
      do 100 i=1,nion
        xmassi = atw(i)*xmassp
        vth(i) = SQRT (2.0 * tiaerg / xmassi)
        rho(i) = xmassi*cee*vth(i)/(rzsqa(i)*charge*bpa)
        xnu(i) = 1.33333*rootpi*ena(i)*zsqa(i)**2*charg4*xlam
     .         / SQRT (xmassi*tiaerg**3)
  100   xnus(i,j) = root2*rmajor*safac*xnu(i)/(del32*vth(i))
c
c  calculate fit parameters of Rawls, Chu, and Hinton for given zeff
c
      do 106 m=1,3
      do 106 k=m,3
      xk0(m,k) = f1i(zin,xk01(m,k),xk0i(m,k))
      a(m,k)   = f14(zin,a1(m,k),a4(m,k))
      b(m,k)   = f14(zin,b1(m,k),b4(m,k))
      c(m,k)   = f14(zin,c1(m,k),c4(m,k))
      if (k .eq. 3)  go to 104
      xk(m,k) = xk0(m,k)
     . * ((1.0/(1.0 + a(m,k) * SQRT (xnuse(j))+b(m,k)*xnuse(j)))
     . * xia(m,k)/(1.95 * rdelta) + 0.5 * c(m,k)**2*xnuse(j)*delta*xipsa
     . / (b(m,k)*(1.0 + c(m,k)*xnuse(j)*del32)))
      go to 106
  104 xk(m,k) =
     .  xk0(m,k)*(1.0/(1.0 + a(m,k) * SQRT (xnuse(j))+b(m,k)*xnuse(j)))
     . *(1.0 / (1.0 + c(m,k)*xnuse(j)*del32))*xia(m,k) / (1.95 * rdelta)
  106 continue
c
c  calculate resistivity
c
      if (resistive .eq. 'hinton') then       ! old, small inv. aspect..
c                                             ..ratio approximate
c        calculate trapped electron fraction
c        (see eq. 4.2-47, gaa)
c
          ftrap(j) =
     .      wneot*xk(3,3)*(1.0-(1.0-1.0/xk0(3,3))*rdelta)*rdelta/xhm2a
          xkap11   = 1.0/(0.29+0.46/(1.08+zeffa))
          ftrap(j) = MIN(ftrap(j),0.995D0) !jmp.ibm
          eta(j)   = wneo(4,4)*xmasse*xnue
     .             /(xkap11*enea*charge**2*(1.0-ftrap(j)))
      else    ! new, valid for inv aspect ratio at least up to 0.4
c
c               note that for zeffa >3.0 the conductivity reduction
c               due to e-e collisons is assumed negligible.
c               neoresist returns effective trapped particle fraction
c               ftrap,etaspitzer,eta(j),and xkap11. (xkap11 is numerically
c               essentially identical to the xkap11 defined above.)
c
          call neoresist (zeffa,wneot,xft,
     .                    xnue,xnuse(j),xmasse,enea,charge,
     .                    ftrap(j),etaspitzer,eta(j),xkap11)
          eta(j) = wneo(4,4)*eta(j)
      end if
c      if( roa(j) .gt. r_elm )then
c         call baloo        ! determines  alphar(1,..nj) 
c         alphar(j) = 1.    ! use time control instead 
c         if(alphar(j) .ge. 0.99 .and. 
c     .    (t_elms .le. time .and. time .le. t_elme))
c     .                          eta(j) = etam_elm*eta(j)
c      endif 


c
c --- debug print
c
      if (idbgprt .eq. 1 .and. jdbgprt .eq. 0)  go to 444
      go to 445
  444 xiaj(j)    = xia(3,3)
      xkj(j)     = xk(3,3)
      xk0j(j)    = xk0(3,3)
      aj(j)      = a(3,3)
      bj(j)      = b(3,3)
      cj(j)      = c(3,3)
      rdeltaj(j) = rdelta
      xhm2aj(j)  = xhm2a
      xnuej(j)   = xnue
      xlamj(j)   = xlam
      safacj(j)  = safac
      vthej(j)   = vthe
      vthj(j)    = vth(1)
      rhoj(j)    = rho(1)
      xnuj(j)    = xnu(1)
      bpaj(j)    = bpa
      del32j(j)  = del32
      eneaj(j)   = enea
      teaj(j)    = tea
      zinj(j)    = zin
      rhoej(j)   = rhoe
c
c  calculate scale lengths of densities, temperatures, and pressure.
c     smooth the derivatives over 2*jstep+1 intervals.
c
  445 jstep     = 0
      jl        = MAX0 (j-jstep  , 1 )
      ju        = MIN0 (j+jstep+1, nj)
      drs       = 0.0
      do 110 js=jl,ju-1
  110 drs       = drs + dr(js)
      deriv     = ABS (ene(ju)-ene(jl))/drs+puny
      slene(j)  = ABS (enea/deriv)
      do 112 i=1,nion
      deriv     = ABS (en(ju,i)-en(jl,i))/drs + puny
  112 slen(j,i) = ena(i)/deriv
      deriv     = ABS (te(ju)-te(jl))/drs + puny
      slte(j)   = tea/deriv
      deriv     = ABS (ti(ju)-ti(jl))/drs + puny
      slti(j)   = tia/deriv
      prese1    = ene(jl)*te(jl)
      prese2    = ene(ju)*te(ju)
      presi1    = 0.666667*(wbeam(jl)+walp(jl))
      presi2    = 0.666667*(wbeam(ju)+walp(ju))
      rhom0     = xmassp*(atw_beam*enbeam(1)+4.0 * enalp(1))
      do 120 i=1,nion
      presi1    = presi1+en(jl,i)*ti(jl)
      presi2    = presi2+en(ju,i)*ti(ju)
  120 rhom0     = rhom0 + xmassp*atw(i)*en(1,i)
      presr1    = prese1+presi1
      presr2    = prese2+presi2
      presra    = 0.5 * (presr1+presr2)
      deriv     = ABS (presr2-presr1)/drs + puny
      slpres(j) = presra/deriv
      slpe      = 1.0/(1.0/slene(j)+1.0/slte(j))
      slpi      = 1.0/(1.0/slene(j)+1.0/slti(j))
c
c  calculate the bohm diffusion coefficient
c
      dbohm = cee*teaerg/(16.0 * charge*btor)
c
c          calculate diagnostic transport coefficients:
c          anomalous electron and ion thermal conductivities
c          xkepp:          Parail and Pogutse model
c          xkeohk:         Ohkawa model
c          xkecha,xkicha:  Chang model
c          xkecar:         Carreras-Diamond model
c          xkedom:         Dominguez model
c          xkwmatdm:       Mattor-Diamond model (phys fluids 31,1180(1988))
c          xkematdm:       Mattor-Diamond model (phys fluids 31,1180(1988))
c
      wpe2      = 4.0 * pi * enea*charge**2/xmasse
      d0        = enea*(cee**2*vthe/(wpe2*rmajor*safac))
      xkepp (j) = d0 * shearp(j)
      xkeohk(j) = d0 * SQRT (1.0 + 0.5 * (xnuse(j)*del32)**2) / pi
      xkparb    = shearp(j)/(rmajor*safac)
      if (shearp(j) .eq. 0.0)  go to 160
c
c     note: evaluation of Z is left for later
c
      zetaim(j) = xnue/(xkparb*vthe)
      do i=1,nion
        rhot      = rho(i)*bpa/btor
        xkicha(j) = xkicha(j) + ena(i)*(rhot**2*xkparb*vth(i)*
     .                         (xkparb*slen(j,i))**(-0.666667))
      end do
      if (shearp(j)*slpres(j) .le. 0.0)  go to 160
      taua0     = rmajor * SQRT (4.0 * pi*rhom0)/btor
      tres      = 4.0 * pi * r(nj)**2/(eta(j)*cee**2)
      betax     = 8.0 * pi * 1.6e-9*presra/btor**2
      betac     = shearp(j)*slpres(j)/(rmajor*safac**2)
      xkecar(j) = 1.5*enea*vthe*r(nj)*safac*(taua0/tres)
     .            *(betax/betac)**1.5
c
c++ form Dominguez electron conductivity
c
      cs        = SQRT (teaerg/xmasip)
      rhos      = cs/omci
      sls       = rmajor*safac/shearp(j)
      xkedom(j) = 9.0*enea*(cs*rhos**2/slpres(j))
     .                * SQRT (tia*slpres(j)/(tea*sls))
  160 continue
c
c --- evaluate the Mattor-Diamond shear viscosity and electron and
c --- ion thermal diffusivites. The average of ion thermal diffusivites
c --- is stored in xkwmatdm(j,nprim+1)
c
      sls     = rmajor*safac/shearp(j)
      tautite = tia/tea
      viona   = 0.0
      if (iangrot .ne. 0)
     .viona = 0.5 * (angrot(j+1)*r2capi(j+1)/rcap(j+1)
     .        + angrot(j)*r2capi(j)/rcap(j))
      kp1             = nprim+1
      xkwmatdm(j,kp1) = 0.0
      xkematdm(j)     = 0.0
      sumden          = 0.0
      do 162 k=1,nprim
      xkwmatdm(j,k) = 0.0
      etaiabs       = ABS (etaim(j,k))
      if (etaiabs .lt. etaioff)  go to 162
      xmassion = atw(k)*xmassp
      csion    = SQRT (teaerg/xmassion)
      if (ABS (viona) .gt. csion)  go to 162
      etaf1 = (1.0 + etaiabs)*tautite*fimpurty
      if (etaf1 .le. 0.0)  go to 162
      etaf2 = csion*etaf1**1.5/(tautite*slen(j,k))
      if (vionzgrd(j) .gt. etaf2)  go to 162
      omcion  = charge*btor/(xmassion*cee)
      rhosion = csion/omcion
      etaf3   = (slen(j,k)*vionzgrd(j)/(2.0*csion))**2
      xkwmatdm(j,k) = wetai(k)*1.3*((etaf1+etaf3**2)**2
     .                *rhosion**2*csion/sls )
      xkwmatdm(j,kp1) = xkwmatdm(j,kp1)+xkwmatdm(j,k)*xmassion*ena(k)
      xkematdm(j) = (wetaie*10.0*del32*(etaf1+etaf3)**3*(rhosion*csion
     .               /sls)**2/xnue)*xmassion*ena(k)+xkematdm(j)
      sumden = xmassion*ena(k)+sumden
  162 continue
      if (sumden .ne. 0.0)  xkematdm(j    ) = xkematdm(j    ) / sumden
      if (sumden .ne. 0.0)  xkwmatdm(j,kp1) = xkwmatdm(j,kp1) / sumden
c
c ----------------------------------------------------------------------
c                        neoclassical transport
c ----------------------------------------------------------------------
c                                            ! jneo = 0 for rawls but jneo > 0
       neorawls :      if (jneo .le. 7) then ! needed to supply missing values 
                                             ! for other models of neoclassical
                                             ! transport such as kapisn.
                                             
c              if (j .eq. 1 ) write(*,170)
c  170         format(/,' Calculating neoclassical transport with Rawls model')
c
c              calculate miscellaneous quantities
c
               tt    = 1.0 + tia/(zeffa*tea)
               sumdz = 0.0
               do 204 k=1,nion
  204          sumdz = sumdz+dzdtea(k)*ena(k)
               sumdz = sumdz*tea/enea
               de    = rdelta*rhoe**2*xnue
               do 206 i=1,nion
  206          di(i) = rdelta*rho(i)**2*xnu(i)
               x     = SQRT (xnusi(j))
               y     = xnusi(j)**2*delta**3
               betag = ((1.17-0.35*x)/(1.0 + 0.7*x)-2.1*y)/(1.0 + y)
               xnu23 = 1.0 + xnuse(j)**2*delta**3
               betagx = betag/xnu23
               zetax  = 1.0/(xkap11*(1.0-ftrap(j)))
               betape = 8.0 * pi*enea*teaerg/bpa**2
c
c              calculate dmn factors
c
               d11 = xk(1,1) + 0.5 * zetax*rdelta*xk(1,3)**2/xhm2a
               d12 = xk(1,2) + 0.5 * zetax*rdelta*xk(1,3)*xk(2,3)/xhm2a
               d22 = xk(2,2) + 0.5 * zetax*rdelta*xk(2,3)**2/xhm2a
               d13 = zetax*xk(1,3)/(xhm2a*betape)
               d23 = zetax*xk(2,3)/(xhm2a*betape)
c
c              calculate the transport coefficients for ion particle transport.
c
c              obtain terms arising from electron-ion collisions.
c              use expressions due to Rawls, Chu, and Hinton.
c
               if (wneo1 .eq. 0.0)  go to 240
               do 220 i=1,nion
                  x        = ena(i)*de*za(i)/zeffa
                  d(i,i,j) = d(i,i,j) + wneo(1,1)*x*tia*d11/
     .                                          (za(i)*tea*ena(i))
                  do k=1,nion
                     xmult = 1.0
                     if (ifus .eq. -1 .and. k .ne. i)  xmult = 0.0
                     d(i,k,j) = d(i,k,j) + xmult * wneo(1,1) 
     .                           * x * za(k) * d11 / enea
                  end do
c
                  k        = nion + 1
                  d(i,k,j) = d(i,k,j) + wneo(1,2)*x*((-1.5+sumdz)*
     .                                                    d11+d12)/tea
                  k        = nion + 2
                  d(i,k,j) = d(i,k,j) + wneo(1,3)*x*(1.0-betagx)*d11/
     .                                                     (za(i)*tea)
                  k        = nion + 3
 220              d(i,k,j) = d(i,k,j) + wneo(1,4)*x*d13/
     .                        (fcapa**2*hcapa*rbpa)
c
c                 save diffusion coefficient of first ion species
c
                  dneo(j) = d(1,1,j)
 240              continue
c
c               calculate the transport coefficients for electron energy transport.
c               obtain terms arising from electron-ion collisions.
c               use expressions due to Rawls, Chu, and Hinton, except that
c               convection is subtracted out.
c
                  if (wneo2 .eq. 0.0)  go to 440
                  i = nion+1
                  x = enea*tea*de
                  do 420 k=1,nion
 420              d(i,k,j) = d(i,k,j)+ wneo(2,1)*x*tt*za(k)*
     .                                        (d12-2.5*d11)/enea
                  k        = nion+1
                  d(i,k,j) = d(i,k,j)+ wneo(2,2)*x*((-1.5+tt*sumdz)*
     .                       (d12-2.5*d11)  +(d22-2.5*d12))/tea
                  k        = nion+2
                  d(i,k,j) = d(i,k,j)+ wneo(2,3)*x*(1.0-betagx)
     .                               *(d12-2.5*d11) / (zeffa*tea)
                  k        = nion+3
                  d(i,k,j) = d(i,k,j)+ wneo(2,4)*x*(d23-2.5*d13)
     .                                /(fcapa**2*hcapa*rbpa)
c
c                  save electron thermal conductivity
c
                  dwifs(j)  = d(i,i,j)
                  xkeneo(j) = d(i,i,j)
c
c                calculate the transport coefficients for ion energy transport.
c                obtain terms arising from electron-ion collisions.
c                use expressions due to Rawls, Chu, and Hinton.
c
 440              i = nion + 2
                  x = -enea*tia*de*betagx/zeffa
                  if (wneo3 .eq. 0.0)  go to 640
                  sumen = 0.0
                  do 612 k=1,nion
                     sumen = sumen + ena(k)
 612                 d(i,k,j) = d(i,k,j) +x*wneo(3,1)*(za(k)+tia/tea)*
     .                                        d11/enea
                  k        = nion + 1
                  d(i,k,j) = d(i,k,j) +x*wneo(3,2)*((-1.5+sumdz)*
     .                       d11+d12)/tea
                  k        = nion + 2
                  d(i,k,j) = d(i,k,j) +x*wneo(3,3)*(1.0-betagx)*sumen*
     .                                       d11/(enea*tea)
                  k        = nion + 3
                  d(i,k,j) = d(i,k,j) +x*wneo(3,4)*d13/
     .                       (fcapa**2*hcapa*rbpa)
c
c                 calculate partial ion heat flux for use in source
c
                  do 615 k=1,nk
 615                 qieneo(j) = qieneo(j) - d(i,k,j)*dudr(k,j)
c
c                  obtain terms arising from like ion-ion collisions.
c                  use expression due to Chang and Hinton.
c
                  zefcor = zeffa*enea/(zsqa(1)**2*ena(1))
                  chfact = 1.0 + 2.85*rdelta - 2.33*delta
                  xk2    = xk20*((chfact/(1.0 + a2 * SQRT (xnus(1,j)
     .                                *zefcor)
     .                 +b2*xnus(1,j)*zefcor))
     .                 *xia(1,1)/(1.95*rdelta)
     .                 +0.5 * c2**2*xnus(1,j)*zefcor*delta*xipsa
     .                 /(b2*(1.0 + c2*xnus(1,j)*zefcor*del32)))
                  do 620 l=1,nion
 620                 xkineo(j) = xkineo(j) +
     .            wneo(3,3)*xk2*zeffa*enea*di(l)/zsqa(l)**2
                  chiineo(j)=xkineo(j)/enasum ! HSJ 9/12/96
                  if(j .eq. j_clamp)chi_clampr = chiineo(j)
                  if ((fusnin .ne. 0.0 .or. ticin .ne. 0.0) .and. 
     .                                w3typ(1) .ne. 0)
     .            go to 640
                  d(i,i,j) = d(i,i,j) + xkineo(j)
 640              continue
c
                  if (ABS (0.5*(q(j+1)+q(j))) .lt. 1.0    )
     .                 dwifs(j) = xkineo(j) ! set dwifs to ion neoclassical value
                  if (ABS (wneo_elct_ifs)     .gt. 1.0e-10)
     .                 dwifs(j) = wneo_elct_ifs * xkineo(j)
c
c                 add classical ion conduction;
c                 helium plasmas may go wrong here (among other places)
c
                  if (w3cla .eq. 0.0)  go to 660
                  xkicla(j) = 0.0
                  do 655 ll=1,nion
                     if (ll .le. nprim)  zefcor = zeffa*enea/
     .                               (zsqa(ll)**2*ena(ll))
                     if (ll .gt. nprim)  zefcor = 1.0
                     xkicla(j) = xkicla(j)+w3cla*gcapa*ena(ll)*xnu(ll)/
     .                    (9.58e3*btor/atw(ll))**2
     .                    /zsqa(ll)*2.0*tia*1.0e3*(9.76e5)**2/atw(ll)*
     .                     zefcor
 655              continue
                  d(i,i,j) = d(i,i,j)+xkicla(j)
c
c                 calculate the transport coefficients corresponding to Ohm's law.
c                 use expressions due to Rawls, Chu, and Hinton.
c
 660              i = nion + 3
                  x = zetax*de*bpa/ra(j)
c                if( roa(j) .gt. r_elm )then
c                      call baloo        ! determines  alphar(1,..nj) 
c                      alphar(j) = 1.    ! use time control instead 
c                      if(alphar(j) .ge. 0.99 .and. 
c     .                (t_elms .le. time .and. time .le. t_elme))
c     .                  x  = etam_elm*x
c                endif 


                  if (jhirsh .ne. 0)  go to 899 ! jhirsh = 0 is Hinton+Hazel. bootstrap
                  do 820 k=1,nion
 820                 d(i,k,j) = d(i,k,j) + x*wneo(4,1)*tt*za(k)*
     .                           xk(1,3)/(2.0*enea)
                  k        = nion+1
                  d(i,k,j) = d(i,k,j) + x*wneo(4,2)*((-1.5+tt*sumdz)
     .                       *xk(1,3) +xk(2,3))/(2.0*tea)
                  k        = nion+2
                  d(i,k,j) = d(i,k,j) + x*wneo(4,3)*(1.0-betagx)*xk(1,3)
     .                       /(2.0*zeffa*tea)
 899              k        = nion + 3 !pick up only diagonal term if jhirsh > 0
                  d(i,k,j) = d(i,k,j) + x*wneo(4,4)
     .                       /(rdelta*betape*fcapa**2*hcapa*rbpa)
c                 print *,'d(kfar,kfar,j),j =', d(i,k,j),j
c
            end if  neorawls




c----------------------------------------------------------------------------
c ----  KAPISN neoclassical at zone center  grid point j : HSj 1/22/09
c ----  determines only values for chiine,chieneo. We assume chiwneo
c ----  (toroidal rotation chi) is equal to chiineo as determined here.
c ---- all other (particle and poloidal b field) required chis are from above
c ----  rawlse all model.
c ---------------------------------------------------------------------------
        IF( jneo .gt. 2 .and. jneo .le. 7) call kapisn_driver(j,jneo)




c -------------------------------------------------------------------------
c ---------------------neoclassical nclass goes here only if single rid point
c nclass is deeloped---------------------
c
c      if (jneo .eq. 1 .and. j .eq. 1 ) write(*,180)
c  180 format (/ ' Calculating neoclassical transport with NCLASS model')
c
C-------------------not wired up -------------------------------
c          call nclass_dr(j)





        bootstrap :  if (0 .lt. jhirsh .and. jhirsh .lt. 30) then
c ----------------------------------------------------------------------
c     --- JHIRSH .ne. 0, HIRSHMAN AAND SAUTER BOOTSTRAP CURRENT MODELS HSJ
c     ****** HIRSHMAN 1978 BOOTSTRAP CURRENT CALCULATION. ******
c     ****** by  C D Challis                              ******
c     Alternative bootstrap current calculation following
c     Hirshman (Phys. Fluids 21 (1978) 1295)
c     (eqn. numbers are given from the above paper)
c
c     Note that this is still a large aspect ratio approx.
c     but it deals with multi-species plasmas.
c
c     This section has been checked by hand at a specific
c     pulse number/time slice/radial location, and against
c     the old version at Zeff = 1 where they should be similar.
c     Beyond this no tests have been performed!!
c     Detailed outputs are installed to help you check that
c     your results are reasonable.
c
c     The electrons are divided between the ion species as
c     above.
c     Here nion+1 represents the electrons.
c     Assume no contribution from the fast ions themselves if
c     jhirsh = 1.0
c     If jhirsh = 2 add fast ions to thermal ions.
c     A future addition to this model should be a proper
c     treatment of the fast ions.
c     Assume beam Z = 1.0
c     Assume F = Bt0*R0/Fcap
c     Assume d(PSI)/dr = Bp0*R0
c     Soften current close to the axis
c ----------------------------------------------------------------------
c
c     Set up initial quantities.
c     XUNITS changes units from  keV.cm-4.gauss-1  to  amp.cm-2,
c     then cancels  statamp  to  amp  correction in SOURCE.
c     XFT is the trapped fraction of thermal particles neglecting
c     collisions which are taken into account elsewhere.
c
      xunits = 1.602e-8*2.998e+9
      dpsidr = -bpa*rmajor
      xres   = cee*eta(j)/ra(j)
c
c  Thermal ions and impuritiy quantities.
c  Coulomb logarithm for self collisions taken from Spitzer.
c
      do 890 i=1,nion
        xmassa(i) = atw(i)*xmassp
        xna(i)    = ena(i)   ! thermal ion density, species k, half grid
        xnfa      = (enbeam(j+1) + enbeam(j))/2.0 +
     .               enalp(j+1) + enalp(j)
        if (i .eq. 1) then
          if (jhirsh .eq. 2 .or. jhirsh .eq. 22)  xna(1) = xna(1) + xnfa
        end if
        xta(i)  = tia
        zeda(i) = za(i)
        xlamii  = 9.745e+9
     .           * SQRT ((1.0e3*xta(i))**3/xna(i)/2.0)/zeda(i)**2
        if (xlamii .gt. 1.0) then
          xlogl(i) = LOG (xlamii)
        else
          xlogl(i) = 0.0
        end if
  890 continue
c
c  Electron quantities.
c
      xmassa(nion+1) = xmasse
      xna(nion+1)    = enea
      xta(nion+1)    = tea
      zeda(nion+1)   = -1.0
      xlogl(nion+1)  = xlam
c
c  Calculate the effective charge scattering species, zstar,
c  and light particle collision parameter, xra, as per p1299,
c  and collisionality parameter from eqn. 37
c
      do 892 i=1,nion+1
        zstar(i) = 0.0
        rasum1   = 0.0
        rasum2   = 0.0
        do 891 k=1,nion+1
          if (xmassa(k) .gt. xmassa(i)) then
            zstar(i) = zstar(i) + xna(k)*zeda(k)*zeda(k)
          end if
          if (xmassa(k) .ge. xmassa(i)) then
            rasum1 = rasum1 + xna(k)*zeda(k)
            rasum2 = rasum2 + xna(k)*zeda(k)*zeda(k)
          end if
 891    continue
        zstar(i) = zstar(i)/xna(i)/zeda(i)/zeda(i)
        if (zstar(i) .lt. 1.0e-10) zstar(i) = 1.0e-10
        xra(i)   = 1.0 - zeda(i)*rasum1/rasum2
        xnush(i) = 6.93e-20*rmajor*safac*xna(i)*zeda(i)**4*
     .             xlogl(i)/(del32*xta(i)*xta(i))
 892  continue
c
c  calculate coefficients from table 1
c
      do 893 i=1,nion+1
        xk0h11 = 0.53 + zstar(i)
        xk0h12 = 0.71 + zstar(i)
        xk0h13 = 1.0 + (0.53 + zstar(i))/zstar(i)/(1.0 + 1.35*zstar(i))
        xk0h23 = 2.5 + (0.72 - .127*zstar(i))/zstar(i)/
     .           (1.0 + 0.28*zstar(i))
c
        xah11  = (0.52 - .42*zstar(i))/(1.0 + 1.35*zstar(i)) +
     .           3.44*zstar(i)
        xah12  = (0.1 + 0.084*zstar(i))/(1.0 + 1.3*zstar(i)) +
     .           0.31*zstar(i)
        xah13  = 0.52 + zstar(i)*(0.13 + 1.53/(1.0 + 2.1*zstar(i))**2)
        xah23  = 0.12 + 0.041*zstar(i)
c
        xbh11  = 0.56 + 0.96*zstar(i)
        xbh12  = 0.26 + 0.35*zstar(i)
        xbh13  = 0.56 + 0.25*zstar(i)
        xbh23  = 0.26 + 0.141*zstar(i)
c
        xch11  = (0.14 + 0.55*zstar(i))/(1.0 + 5.0 * zstar(i)) +
     .           .25*zstar(i)
        xch12  = (0.072 + 0.15*zstar(i))/(1.0 + 3.0*zstar(i)) +
     .           .081*zstar(i)
        xch13  = 0.25 + 0.051*zstar(i) -
     .   (0.11 + 0.74*zstar(i))/(1.0 + 2.29*zstar(i) + 0.85*zstar(i)**2)
        xch23  = 0.21 + 0.022*zstar(i) -
     .   (0.14 + 0.84*zstar(i))/(1.0 + 3.28*zstar(i) + 1.22*zstar(i)**2)
c
        xdh11  = (0.7  + 0.78*zstar(i))/(1.0 + zstar(i)) + 0.51*zstar(i)
        xdh12  = (0.42 + 0.62*zstar(i))/(1.0 + 1.42*zstar(i)) +
     .            0.29*zstar(i)
        xdh13  = (0.7  + 2.99*zstar(i))/(1.0 + 5.0 * zstar(i)) +
     .            0.22*zstar(i)
        xdh23  = (0.43 + 1.95*zstar(i))/(1.0 + 5.0 * zstar(i)) +
     .            0.15*zstar(i)
c
c  calculate analytic functions of eqn 40
c
        xk(1,1) = xk0h11/(1.0 + SQRT (xah11*xnush(i)) + xbh11*xnush(i))/
     .       (1.0 + SQRT (xch11*xnush(i)*del32) + xdh11*xnush(i)*del32)
        xk(1,2) = xk0h12/(1.0 + SQRT (xah12*xnush(i)) + xbh12*xnush(i))/
     .       (1.0 + SQRT (xch12*xnush(i)*del32) + xdh12*xnush(i)*del32)
        xk(1,3) = xk0h13/(1.0 + SQRT (xah13*xnush(i)) + xbh13*xnush(i))/
     .       (1.0 + SQRT (xch13*xnush(i)*del32) + xdh13*xnush(i)*del32)
        xk(2,3) = xk0h23/(1.0 + SQRT (xah23*xnush(i)) + xbh23*xnush(i))/
     .       (1.0 + SQRT (xch23*xnush(i)*del32) + xdh23*xnush(i)*del32)
c
c  calculate coefficients from eqn 39
c
        xu1s(i) = xft*xra(i)*xk(1,3)
        xu2s(i) = xft*xra(i)*(xk(2,3) - 2.5*xk(1,3))
        if (i .eq. nion+1) then
          xu1(i) = 1.38*xna(i)*xta(i)*rmajor*safac*xnush(i)*
     .             xk(1,1)/vthe
          xu2(i) = 1.38*xna(i)*xta(i)*rmajor*safac*xnush(i)*
     .             (xk(1,2) - 2.5*xk(1,1))/vthe
        else
          xu1(i) = 1.38*xna(i)*xta(i)*rmajor*safac*xnush(i)*
     .             xk(1,1)/vth(i)
          xu2(i) = 1.38*xna(i)*xta(i)*rmajor*safac*xnush(i)*
     .             (xk(1,2) - 2.5*xk(1,1))/vth(i)
        end if
 893  continue
c
c  Calculate transport coefficients from eqn 19 multiplied by R0/F(psilim)
c  Change units from  keV.cm-4.gauss-1  to  amp.cm-2, then cancel
c  statamp to amp correction in SOURCE.
c
      do i=1,nion+1
        usum0 = 0.0
        do k=1,nion+1
          usum0 = usum0 + xu1(k)
        end do
        usum1 = 0.0
        usum2 = 0.0
        do k=1,nion+1
          usum1 = usum1 + zeda(k)*xna(k)*(xu1(i)*xu1s(k)/usum0)
          usum2 = usum2 + zeda(k)*xna(k)*(xu2(i)*xu1s(k)/usum0)
        end do
        usum1 = usum1 - zeda(i)*xna(i)*xu1s(i)
        usum2 = usum2 - zeda(i)*xna(i)*xu2s(i)
        xl31  = rmajor*xta(i)*usum1*xunits/zeda(i)/fcapa
        xl32  = rmajor*xta(i)*usum2*xunits/zeda(i)/fcapa
c
c       As the bootstrap calculation is invalid near the axis,
c       we add a little realism by multiplying the coefficients by the
c       ratio r/(poloidal field gyro-radius) if it is less than unity.
c
        rlar   = 2.5e+15 * SQRT (2.0*xmassa(i)*xta(i))/zeda(i)/bpa
        rratio = ra(j)/rlar
c**HSJ
        if (jhirsh .eq. 11 .or. jhirsh .eq. 22)  rratio = 1.0
        axcoef = MIN (one, ABS (rratio))
        xl31   = xl31 * axcoef
        xl32   = xl32 * axcoef
c
c       For diagnostic purposes only normalised L31 and L32 are
c       calculated so that they are multiplied by
c           RTn*dlog( )/d(psi)
c       to get <j.B/Bto> for each term.
c
        if (i .eq. nion+1) then
          xl31ne(j)  = usum1*axcoef/xna(i)/zeda(i)/fcapa
          xl32ne(j)  = usum2*axcoef/xna(i)/zeda(i)/fcapa
        else
          xl31n(i,j) = usum1*axcoef/xna(i)/zeda(i)/fcapa
          xl32n(i,j) = usum2*axcoef/xna(i)/zeda(i)/fcapa
        end if
c
c       End of literal coding of Hirshman.  Hereafter the coefficients
c       are converted to a form suitable for SOURCE.
c
c       Calculate coefficients for bootstrap calculation converting
c       to d( )/dr from dlog( )/d(psi) and normalising to the
c       conductivity as per the form in SOURCE.
c
        if (i .eq. nion+1) then
          xden(j)   = xres* xl31 / xna(i) / dpsidr
c         xdet(j)   = xres*(xl31 + xl32)/xta(i)/dpsidr
          xdet(j)   = xres*(xl31 *(1.0 + sumdz)+ xl32)/xta(i)/dpsidr
        else
          xdin(i,j) = xres* xl31 / xna(i) / dpsidr
          xdit(i,j) = xres*(xl31 + xl32)/xta(i)/dpsidr
        end if
c
      end do
c
c  Rearrange the coefficients to divide the electron term
c  between the ions and leave the fast ion electrons to the
c  end.  A new coefficient 'dfast' is added to account for
c  them.  Add the wneo switches as in the old calculation.
c
c  k = nion+1: electron temperature gradient
c
      d(nion+3,nion+1,j) = wneo(4,2)*xdet(j)
      d(nion+3,nion+2,j) = 0.0
c
      do i=1,nion
c
c       k = nion+2: ion temperature gradient
c
        d(nion+3,nion+2,j) = d(nion+3,nion+2,j) + wneo(4,3)*xdit(i,j)
c
c       k = 1 to nion: ion density gradients with electron density
c                      gradients included
c
        d(nion+3,i,j) = wneo(4,1)*(xdin(i,j) + xden(j)*zeda(i))
c
      end do
c
c     fast ion electron density gradient
c
      dfast(j) = wneo(4,1)*xden(j)
c
c     fast ion density gradient
c
      dfion(j) = 0.0
      if (jhirsh .eq. 2 .or. jhirsh .eq. 22)
     .     dfion(j) = wneo(4,1) * xdin(1,j)
c
c ****** END OF BOOTSTRAP CURRENT MODIFICATION USING HIRSHMAN ******
c
        else if (jhirsh .ge. 80 .and. jhirsh .lt. 95) then  bootstrap 
c
c ----------------------------------------------------------------------
c --- Hirshman arbitrary aspect ratio, small collisonality limit,
c --- electron, ion plasma model
c --- reference: Hirshman, Phys. Fluids, 31(10), Oct 88
c --- the calculations in subroutine HIRSH88 are patterned after the hirshman 78
c --- model given above so that no changes have to be made in subroutine SOURCE
c --- HSJ
c ----------------------------------------------------------------------
c
              call hirsh88 (jhirsh,j,nion,fcapa,enea,tea,
     .                      ena,tia,zeffa,bpa,cee,eta,ra,rmajor,
     .                    enbeam,enalp,sumdz,za,wneo,xia,xhm2a,codeid)
c                           d,dfast,dfion,xden,xdit,xdet,xdin)
c
       else if (jhirsh .gt. 109 .and. jhirsh .le. 112) then bootstrap 
c            bootstrap from Sauter 99 model:
             rrmajor = 0.5*(ravg_r(j+1)+ravg_r(j))
             epsa = (eps(j+1)+eps(j)) * 0.5
             call sauter  (jhirsh, j, nion, fcapa, enea, tea, ena, tia,
     .                     zeffa,q,xft,epsa,rrmajor,
     .                     cee, eta, ra, rmajor, enbeam,
     .                     enalp, sumdz, za, wneo,btor,
     .                     d, dfast, dfion, xden, xdit, xdet, xdin)

      end if   bootstrap



c
c ----------------------------------------------------------------------
c --- calculate transport coefficients for toroidal rotation.
c --- assume diagonal model for now with neoclassical expression given by
c --- Hinton and Wong and Wong (Phys. Fluids 28,3082(1985) and 30,818(1988)
c --- drift wave model of mattor-diamond (phys fluids 31,1180(1988)
c --- or use empirical take your pick type input model:
c --- NOTE THAT IF THE STAEBLER-HINTON MODIFICATION FOR SHEAR
c --- TURBULENCE SUPPRESSION IS ON (I.E., IFSFLAG > 0) THEN
c --- THE MOMENTUM DIFFUSIVITY IS ADDED IN THE RLW OR HSIEH SECTIONS BELOW !!!!
c --- a take-your-pick momentum diffusivity may also be added
c ----------------------------------------------------------------------
c
      tor_rot :    if (iangrot .ne. 0) then
          i        = nion + 4
          dmassden(j) = 0.0
c
c ---     determine momentum diffusivity from neoclassical theory
c
          xmuib      = 0.1 *  delta2
          xmuip      = 0.6 * (bpa/btor)**2
          chiwneo(j) = 0.0
c
         do k=1,nprim
c
c           banana regime
c
            if      (xnus(k,j) .le. 1.0) then
              xmui = xmuib
            else if (1.0 .lt. xnus(k,j) .and. 
     .                             xnus(k,j) .lt. del32i) then
c
c             plateau regime: assume linear connection (banana to Pfirsch-Schluter)
c
              xmui = xmuib + (xnus(k,j)-1.0)*(xmuip-xmuib)/(del32i-1.0)
            else             ! Pfirsch-Schluter regime
              xmui = xmuip
            end if
            xmui       = xmui * rho(k)**2 * xnu(k)
            chiwneo(j) = chiwneo(j)+xmui*atw(k)*xmassp*ena(k)*rmajor**2 !(g cm)/sec
            dmassden(j)   = dmassden(j)+atw(k)*xmassp*ena(k)
          end do
c
          chiwneo(j)   = chiwneo(j)/(rmajor**2*dmassden(j)) !cm**2/sec
          if (xkwmatdm(j,nprim+1) .eq. 0.0) then
            chiwmatm(j) = chiwneo(j)
          else
            chiwmatm(j) = xkwmatdm(j,nprim+1)
          end if
c
c --- get the diffusion coefficient and the momentum diffusivity
c --- based on neoclassical and/or Mattor-Diamond models
c
          if (iwangrot .eq. 0) then
c
c ---         neoclassical only
c
c ---         momentum diffusivity may be added later in TYP, RLW or Hsieh models
c
              xkangrot(j) = chiwneo(j)*wneo(5,5)
              d(i,i,j)    = d(i,i,j)+xkangrot(j)*dmassden(j)*rmajor**2
          else if (iwangrot .eq. -2) then
c
c ---         add the Mattor-Diamond model to the neoclassical, note that
c ---         wneo(5,5) can be zero due to a)wneo(5,5) = 0 on input
c ---         or b)xmtmdifs(nprim+3) = 0  (recall that
c ---         wneo(5,5) = wneo(5,5)*xmtmdifs(nprim+3)  if iwangrot = -2)
c
              xkangrot(j) = wneo(5,5)*chiwneo(j) + xkwmatdm(j,nprim+1)
              d(i,i,j)    = d(i,i,j)+(dmassden(j)*rmajor**2)*xkangrot(j)
          end if
c
c --- xkangrot is set by empirical input model
c
          if (iwangrot .ne. 0 .and. iwangrot .ne. -2) then
               xkangavg    = (xkangrot(j)+xkangrot(j+1))*0.5
              d(i,i,j)    = dmassden(j)*xkangavg*rmajor**2+d(i,i,j)
              dkapomeg(j) = d(i,i,j) / (2.0*xmassp)          ! save for output
          end if

          if(iwangrot .eq. -3)then
              chiwneo(j)  = chiineo(j)        !cm**2/sec
              xkangrot(j) = chiwneo(j)*wneo(5,5)
              d(i,i,j)    = dmassden(j)*xkangrot(j)*rmajor**2 ! (g cm)/sec
          endif
c
c --- explicitly zero off-diagonal elements
c
          do k=1,i-1
            d(i,k,j) = 0.0
            d(k,i,j) = 0.0
          end do
c
c --- extrapolate to boundary (j = nj, used for printout only)
c
          if (j .eq. nj-1) then
            call extrap (ra(nj-2), ra(nj-1), r(nj),
     .                   xkangrot(nj-2), xkangrot(nj-1), xkangrob)
            xkangrot(nj) = xkangrob
          end if



      end if tor_rot







c
c ----------------------------------------------------------------------
c                       take-your-pick transport
c ----------------------------------------------------------------------
c
c  calculate the transport coefficients
c
      if (wtyp .eq. 0.0)  go to 1130
c
c  map splines for w1-4typ,vtyp profiles to mesh
c
      call wprof
      xkx = 5.0e17 * (enea*1.0e-13)**typa(1)
     .             * tea**typa(2)
     .             * tia**typa(3)
     .             * ABS (curdna)**typa(4)
     .             * safac**typa(5)
     .             *(ra(j)/r(nj))**typa(6)
     .  *((1.0-typa(8))*(1.0-ra(j)**2/r(nj)**2)**typa(7)+typa(8))
      pheatmul = 1.0
      if (typa(9) .ne. 0.0 .or. typa(10) .ne. 0.0) then
        if (typa(16) .ne. 0.0) then
c
c --- set up the total heating powers (units are watts)
c
          pohmsa   = pohms
          pbeamesa = pbeames
          pbeamisa = pbeamis
          prfesa   = prfes
          prfisa   = prfi
        else
c
c --- set up the local heating powers (units are watts)
c
          pohmsa   = (qohm  (j) + qohm  (j+1)) * 0.5
          pbeamesa = (qbeame(j) + qbeame(j+1)) * 0.5
          pbeamisa = (qbeami(j) + qbeami(j+1)) * 0.5
          prfesa   = (qrfe  (j) + qrfe  (j+1)) * 0.5
          prfisa   = (qrfi  (j) + qrfi  (j+1)) * 0.5
c
c --- convert from keV/cm**3*sec to watts, dv = 4 * pi * pi*R0*h*r*dr
c
          const1   = const    * hcapa * ra(j) * (r(j+1)-r(j))
          pohmsa   = pohmsa   * const1
          pbeamesa = pbeamesa * const1
          pbeamisa = pbeamisa * const1
          prfesa   = prfesa   * const1
          prfisa   = prfisa   * const1
        end if
c
c --- get the local or global power
c
        pheattot = pohmsa+pbeamesa+pbeamisa+prfesa+prfisa
        if (pheattot .le. pohmsa)  pwrexp = typa(9)
        if (pheattot .gt. pohmsa)  pwrexp = typa(10)
c
c --- get the scaling factor in units of Mw
c
        pheattot = (typa(11)*pohmsa+typa(12)*pbeamesa
     .             +typa(13)*pbeamisa+typa(14)*prfesa
     .             +typa(15)*prfisa)/1.0e6
        if (pheattot .gt. 0.0 .and. pwrexp .ne. 0.0)
     .  pheatmul = pheattot**pwrexp
      end if
      xkx = xkx * pheatmul
c
      do jisl=1,nnisl
****    psmmn = risl(jisl)-wdisl(jisl)/2.0
        rsmmn = risl(jisl)-wdisl(jisl)/2.0
        rspmn = risl(jisl)+wdisl(jisl)/2.0
        rsmn  = risl(jisl)
        wmn   = wdisl(jisl)
        if (r(j) .gt. rsmmn .and. r(j) .lt. rspmn) xkx = xkx
     .      *(1.0 + (wisl) * EXP (-3.0*((r(j)-rsmn)/0.5/(wmn+0.1))**2))
        if (xkx .gt. 1.0)  go to 851
      end do
c
  851 dx          = xkx/enea
c     to avoid negative diffusion coefficients:
      xkdum       = coefc(j)*w1typt*0.5*(w1fact(j)+w1fact(j+1))*dx
      dtyp(j)     = MAX(xkdum,w1typmin)
      vpinch(j)   = vtypt*0.5*(vfact(j)+vfact(j+1))*ra(j)/r(nj)
      xkdum       = coefa(j)*w2typt*0.5*(w2fact(j)+w2fact(j+1))*xkx
      xketyp(j)   = MAX(xkdum,w2typmin)
      xkdum       = coefb(j)*w3typt*0.5*(w3fact(j)+w3fact(j+1))*xkx
      xkityp(j)   = MAX(xkdum,w3typmin)
      xkdum       = coefd(j)*w4typt*0.5*(w4fact(j)+w4fact(j+1))*
     .              xkx/5.0e13
      xkangtyp(j) = MAX(xkdum,w4typmin)
c     add critical gradient option:
      if(freeze_xte .eq.0)then
      if(w2cg .gt. 0.0)then
         tegrad = (te(j+1)-te(j))/(dr(j))
         if(ABS(tegrad) .gt. w2cg)xketyp(j) = w2cgm*cparam*xketyp(j)
         xketyp_prev(j)= xketyp(j)
      else if(w2cg .lt. -100.)then
c        set single step change in chi at fixed rho:
         if(r(j)/r(nj) .gt.-100. -w2cg)xketyp(j) = 
     .                                  w2cgm*cparam*xketyp(j)
         xketyp_prev(j)= xketyp(j)
      endif
      else
         xketyp(j)= xketyp_prev(j)
      endif
      if(w3cg .gt. 0.0)then
         tigrad = (ti(j+1)-ti(j))/(dr(j))
         if(ABS(tigrad) .gt. w3cg)xkityp(j) = w3cgm*cparam*xkityp(j)
      endif
c
c  relaxation step
c
      if (itercorct .gt. 1) then
        xketyp(j)   = (1.0 - relaxtyp)*xketypv(j) +
     .                       relaxtyp*xketyp(j)
        xkityp(j)   = (1.0 - relaxtyp)*xkitypv(j) +
     .                       relaxtyp*xkityp(j)
        xkangtyp(j) = (1.0 - relaxtyp)*xkangtypv(j) +
     .                       relaxtyp*xkangtyp(j)
        dtyp(j)     = (1.0 - relaxtyp)*dtypv(j) +
     .                       relaxtyp*dtyp(j)
      end if
      xketypv(j)   = xketyp(j)
      xkitypv(j)   = xkityp(j)
      xkangtypv(j) = xkangtyp(j)
      dtypv(j)     = dtyp(j)
c
      do 1120 i=1,nion
      d(i,i,j)  = d(i,i,j) + dtyp(j)
      d(i,nion+1,j) = d(i,nion+1,j) + w12typ*dx*ena(i)/tea
 1120 d(i,nion+2,j) = d(i,nion+2,j) + w13typ*dx*ena(i)/tia
      i        = nion + 1
      d(i,i,j) = d(i,i,j) + xketyp(j)
      i        = nion + 2
      d(i,i,j) = d(i,i,j) + xkityp(j)
      i        = nion + 4
      d(i,i,j) = d(i,i,j) + xkangtyp(j)*dmassden(j)*rmajor**2
c
c ----------------------------------------------------------------------
c     standard drift-ballooning model
c ----------------------------------------------------------------------
c
 1130 if (wstd .eq. 0)  go to 1150
      xkp    = 1.0/3.0/rhos
      omgstr = xkp*rhos*cs/slene(j)
      omgoms = 1.0/(1.0 + xkp**2*rhos**2)
      omgt   = vthe/rmajor/safac
      etae   = (ene(j)+ene(j+1))/(ene(j)-ene(j+1)+puny)
     .        *(te(j)-te(j+1))/(te(j)+te(j+1))
      etai   = (ene(j)+ene(j+1))/(ene(j)-ene(j+1)+puny)
     .        *(ti(j)-ti(j+1))/(ti(j)+ti(j+1))
c
      facbeta =  1.0
      if (alphar(j) .lt. 1.0)
     .facbeta = (1.0 / (1.0-alphar(j)))**ddebug(43)
      if (alphar(j) .ge. 1.0 .and. ddebug(43) .gt. 0.0)
     .facbeta =  100.0
c
      wetainw = ddebug(8)
      etaths  = 1.0 + ddebug(6)
      xths    = etai - etaths
      if (xths .gt.  10.0)  xths =  10.0
      if (xths .lt. -10.0)  xths = -10.0
      fths   = 1.0/(1.0 + EXP (-6.0/(1.9+ddebug(28))*xths))
      xketai = fths*root2*(slene(j)/rmajor)**0.5 * (ABS (etai))**0.5
     .        *(tia/tea)**0.5/xkp/rhos/(1.0+(xkp*rhos)**2)**0.5
     .        *cs/slene(j)*rhos**2
      xketai = enea*xketai*facbeta*safac**ddebug(44)
c
      colfac(j,1) = xnue/omgt
      colmul      = colfac(j,1)
      xone        = 1.0 - ddebug(27)
      if (colmul .lt. xone)
     .colmul      = 1.0 - ddebug(34)
      colfac(j,2) = root2/xnuse(j)/omgt*omgstr/rdelta
      tpfac       = colfac(j,2)
      yone        = 1.0 - ddebug(46)
      if (tpfac .gt. yone)
     .tpfac       = yone
      tufac       = 1.0
      zoone       = (rmajor*slpi/slene(j)**2)**(0.5)
      if (tpfac .gt. zoone) tufac = zoone
      dstdtp      = (1.0 + ddebug(1))*tpfac*rdelta*omgstr/xkp**2
      dstdcp      = (1.0 + ddebug(7))*colmul*omgstr/omgt*omgstr/xkp**2
      xths        = ddebug(41)-etae
      if (xths .gt.  10.0) xths =  10.0
      if (xths .lt. -10.0) xths = -10.0
      feths = 1.0/(1.0 + EXP (-6.0/(1.0 + ddebug(42))*xths))
      if (ddebug(41) .eq. 0.0) feths = 1.0
      dstdcp = dstdcp*feths
      xetaif = (1.0 - ddebug(29) * fths
     .       + fths * tufac * ddebug(20) * slene(j)**2 / rmajor/slpi)
      dstdtp = dstdtp*xetaif*facbeta*safac**ddebug(44)
      dstdcp = dstdcp*(xetaif-ddebug(45))*facbeta*safac**ddebug(44)
c
      tpinch = 1.0
      rpinch = ddebug(5)*r(nj)
      if (ra(j) .gt. rpinch) tpinch = 0.0
      dstdn  = dstdtp*(2.0*(1.0-omgoms))*2.0/rootpi
     .       + dstdcp*((1.0-omgoms)) * pi * 2.0/rootpi
      dstdt  = dstdtp * 3.0 * 2.0 / rootpi
     .       - dstdcp*1.0/2.0*tpinch * pi * 2.0/rootpi
      dstd   = dstdn+etae*dstdt
      xkstdt = dstdtp*63.0/2.0*2.0/rootpi+dstdcp*9.0/4.0 * pi*2.0/rootpi
      xkstdt = xkstdt*enea
      xkstdn = dstdtp*7.0*(1.0-omgoms)*2.0/rootpi
     .       - dstdcp*1.0/2.0*(1.0-omgoms) * pi * 2.0/rootpi
      xkstdn = xkstdn*enea
      dtest(j)  = dstd
      dtestt(j) = dstdt*etae
      xktest(j) = xkstdt+xkstdn/etae
c
c --- the following is a simplified model for electron drift waves
c
      if (ddebug(25) .eq. 0.0)  go to 1131
      dstdt  = 0.0
      xkstdn = 0.0
      dstdcp = dstdcp*delta**ddebug(36)
      dstd   = dstdtp+dstdcp
      if (ddebug(35) .eq. 1.0)  dstd = MAX (dstdtp, dstdcp)
      dstdn  = dstd
      xkstdt = enea*5.0/2.0*dstd
c
 1131 omg    = omgoms * omgstr
      delomg = ABS (dstd)*xkp**2
      dstdi  = delomg/(omg**2+delomg**2)*omgstr**2/xkp**2
      xkistd = ABS (dstd)*(5.0/2.0*etai)*(1.0-ddebug(26))
      xkistd = xkistd*enea/etai
      exchg  = enea*tea * ABS (dstd)/slene(j)**2*omgoms*(1.0-ddebug(26))
      qexch(j)  = wstd*exchg
      qexch(nj) = qexch(nj-1)
      betcrt    = 2.0 * slpe / (rmajor * safac**2)
      betae     = 8.0 * pi * 1.6e-9 * enea * tea / btor**2
      xxkp      = xkp
      xkrmhd    = 1.0/rhos*(slpe/rmajor/safac)**0.5 * (xnue*slpe/cs
     . *xmasse/xmasip)**(-0.5)*(betae/betcrt)**(-0.25)*(betae/2.0)**0.25
      if (xkrmhd .lt. xkp)  xxkp = xkrmhd
      xkbb      = rmajor*safac*vthe*(betae/betcrt/rmajor/safac/xxkp)**2
      xkbb      = enea*xkbb
      xkbtst(j) = xkbb
      wbb       = ddebug(9)
      wmrnv     = ddebug(4)
      wmrnvi    = ddebug(4)*ddebug(19)
      fact2     = 0.0
      immax     = 10-ddebug(33)
c
      do im=2,immax
        xm1 = im
        xm1 = xm1-ddebug(32)
        xm2 = xm1+1.0-ddebug(18)
        if (safac .ge. xm1 .and. safac .le. xm2)
     .    fact2 = 4.0/(xm1-xm2)**2*(safac-xm1)*(xm2-safac)
      end do
c
      if (ddebug(18) .eq. 0.0)  fact2 = MAX (safac-2.0, zero)
      xkmrnv = fact2*1.0e19
      wpse   = ddebug(17)
      dpse   = safac**2*xnue*(xmasse/xmasip)**0.5 * rhos**2
      xkpsee = (5.0/2.0)*enea*dpse
      xkpsei = (5.0/2.0)*enea*dpse
c
      ydebug(j, 1) = colfac(j,1)
      ydebug(j, 2) = colfac(j,2)
      ydebug(j, 3) = betae/betcrt
      ydebug(j, 4) = xkp/xkrmhd
      ydebug(j, 5) = etai
      ydebug(j, 9) = xkstdt*wstd
      ydebug(j,10) = enea*5.0/2.0*dstdtp*wstd
      ydebug(j,11) = enea*5.0/2.0*dstdcp*wstd
      ydebug(j,12) = xketai*wetainw
      ydebug(j,13) = xkmrnv*wmrnv
      ydebug(j,14) = xkpsee*wpse
c
      do i=1,nion
        d(i,i,j)      = d(i,i,j)+wstd*dstdn+wpse*dpse
        d(i,nion+1,j) = d(i,nion+1,j)+wstd*ena(i)/tea*dstdt
        d(nion+1,i,j) = d(nion+1,i,j)+wstd*xkstdn/enea*tea
      end do
c
      i        = nion+1
      d(i,i,j) = d(i,i,j)+wstd*xkstdt+wbb*xkbb
     .         + wmrnv*xkmrnv+wpse*xkpsee
      i        = nion + 2
      d(i,i,j) = d(i,i,j)+wstd*xkistd+wetainw*xketai+wpse*xkpsei
     .         + wmrnvi*xkmrnv
c
 1150 if (ddebug(38) .eq. 0.0)  go to 1151
c
c ----------------------------------------------------------------------
c recalculate Ohkawa model 9/24/85 rew
c ----------------------------------------------------------------------
c
      alf2   = ddebug(38)**2
      gamohk = slene(j)/slte(j)*alf2*vth(1)**2/rmajor**2/safac**2
     .       * rdelta*delta/xnue*delta**ddebug(37)
      xnuohk = xnue/delta
      pfohk  = (qohm(j)+qbeame(j))/enea/tea
      xkeohk(j) =
     .  enea*16.0*dbohm * SQRT (rdelta*pfohk*(1.0/gamohk+1.0/xnuohk))
c
      i        = nion+1
      d(i,i,j) = d(i,i,j)+ddebug(39)*xkeohk(j)
      i        = nion+2
      d(i,i,j) = d(i,i,j)+ddebug(40)*xkeohk(j)
c
c ----------------------------------------------------------------------
c                sawtooth-induced transport (simple model)
c ----------------------------------------------------------------------
c
 1151 if (wsaw .eq. 0.0)  go to 1240
      fact1     = MAX (1.0 - safac, zero) * grsq
      if (fact1 .ne. 0.0) jsaw = j
      dsaw(j)   = w1saw*fact1*dbohm
      xkesaw(j) = w2saw*fact1*enea*dbohm
      enasum    = 0.0D0
      do 1210 i=1,nion
 1210 enasum    = enasum + ena(i)
      xkisaw(j) = w3saw*fact1*enasum*dbohm
      do 1220 i=1,nion
 1220 d(i,i,j)  = d(i,i,j) + dsaw(j)
      i         = nion+1
      d(i,i,j)  = d(i,i,j) + xkesaw(j)
      i         = nion+2
      d(i,i,j)  = d(i,i,j) + xkisaw(j)
c
c ----------------------------------------------------------------------
c                ballooning mode-induced transport
c ----------------------------------------------------------------------
c
 1240 if (ibaloo .ne. 0) then
        i        = nion+1
        d(i,i,j) = d(i,i,j) + xkebal(j)
        i        = nion+2
        d(i,i,j) = d(i,i,j) + xkibal(j)
      end if








c
c ----------------------------------------------------------------------
c end loop over mesh points label 2000
c ----------------------------------------------------------------------
c
 2000 continue





      if(j_clamp .gt. 0)then
         slope_clamp = (chi_clampr - chi_clamp)/r(j_clamp)
      else
         slope_clamp =0.0
      endif
      
      if (idbgprt .eq. 0 .or. jdbgprt .eq. 1)  go to 446
  447 format (/// (2x,i3,7(2x,1pe15.7)))
      jdbgprt = 1
      write (iodbgprt, 447) (j,r(j),xiaj(j),xkj(j),xk0j(j),aj(j),
     .  bj(j),cj(j),j = 1,nj)
      write (iodbgprt, 447) (j,r(j),ra(j),rdeltaj(j),xhm2aj(j),xnuse(j),
     .  xnuej(j),xlamj(j),j = 1,nj)
      write (iodbgprt, 447) (j,r(j),safacj(j),te(j),ti(j),vthej(j),
     .  vthj(j),rhoj(j),j = 1,nj)
      write (iodbgprt, 447) (j,r(j),xnuj(j),bpaj(j),zinj(j),ftrap(j),
     .  rhoej(j),xnusij(j),j = 1,nj)
      write (iodbgprt, 447) (j,r(j),eta(j),curden(j),time,
     .  del32j(j),eneaj(j),teaj(j),j = 1,nj)
      call giveupus(iodbgprt)
      close (unit = iodbgprt)
c
  446 if (wsaw .eq. 0.0 .or. jsaw .le. 1)  go to 2100
c
c ----------------------------------------------------------------------
c     adjust ftrap and eta to give a constant resistivity inside the
c     q=1 - surface.
c     the constant value of eta is chosen such that the Ohmic power
c     is unchanged by the adjustment.
c     ftrap is changed accordingly.
c ----------------------------------------------------------------------
c
      do 2090 j=1,jsaw
      cur2(j) = curohm(j)*curohm(j)
c
c     at the initial step, curohm has not yet been set.
c
      if (time .eq. time0)  cur2(j) = curden(j)*curden(j)
 2090 pohm1(j) = eta(j)*cur2(j)
c
      call trapv (r,cur2,hcap,jsaw,cur2t)
      call trapv (r,pohm1,hcap,jsaw,pohm1t)
      etabar = pohm1t/cur2t
      i      = nion+3
      do 2095 j=1,jsaw
      ftrapn =  1.0-(1.0-ftrap(j))*eta(j)/etabar
      frat   = (1.0-ftrap(j))/(1.0-ftrapn)
c
c     reset transport coefficient corresponding to Ohm's law resistivity.
c     account for the adjustment in d(4,1 = 1,..nion+2,j) needed
c     because resistivity changed for jhirsh > 0 models. HSJ
c
      d(i,i,j) = d(i,i,j)*frat
          if (jhirsh .gt. 0) then
              dumyf = etabar/eta(j)
              do in=1,nion+2
                d(i,in,j) = d(i,in,j)*dumyf
              end do
              dfast(j) = dfast(j)*dumyf
              dfion(j) = dfion(j)*dumyf
          end if
      ftrap(j) = ftrapn
      eta  (j) = etabar
 2095 continue
c
c ----------------------------------------------------------------------
c --- add Rebut-Lallia-Watkins model:
c --- eta, q and ftrap are on the half grid.
c --- curden,fcap,te,ti,zeff,ene are on the full grid
c --- rebutlal1 returns xchierl,xchiirl ,drebut on the half grid,
c --- in cm**2/sec,through include file rebut.i
c --- d(i,i,j) is in 1.0/(cm sec) for i = nion+1,nion+2 (i.e., te,ti)
c --- d(i,i,j) is in cm**2/sec for i = 1,2..nion
c --- these results will affect the calculations only if wrebut .ne. 0.0
c --- and the appropriate simulation mode is set. we calculate them even
c --- if wrebut = 0.0 however so that an appropriate printout and plot
c --- can be made. timerebut is used to indicate the time when
c --- subroutine was last called.
c --- there are two single fluid chi effectives. the first, xeffrl,
c --- calculated below, is based on the RLW models of chie and chii.
c --- the second, xeffsf, cannot be calculated here since we may want
c --- to get xeffsf from experimental rather than theoretical data.
c --- thus we have to calculate xeffsf in subroutine PSOURC
c ----------------------------------------------------------------------
c
 2100       if (tirlw .eq. 0) then
              call copya (ti, tieff, nj)      ! use real ti in RLW model
            else                              ! use eff. ti in RLW model
              do j=1,nj
                ensum = 0.0
                  do i=1,nprim      ! have to use nprim to be consistent
                    ensum = en(j,i) + ensum
                  end do
                  tieff(j) = (ti(j)*ensum+0.666667*(wbeam(j)+walp(j)))
     .                     / (ensum+enbeam(j)+enalp(j))
              end do
            end if
c
            timerebut = time
            call rebutlal (eta,curden,fcap,te,tieff,q,zeff,ftrap,
     .                     ene,r,dr,atw(1),nj,nout,ncrt)
            do j=1,nj-1
                enea = (ene(j+1)+ene(j))*0.5
                teprime = (te(j+1)-te(j))/dr(j)
                tiprime = (tieff(j+1)-tieff(j))/dr(j)
                ensum = 0.0
                do i=1,nprim    ! have to use nprim to be consistent
****            do i=1,nion     ! with calculation of xkiinv (in PSOURC)
                     enav = 0.5 * (en(j,i)+en(j+1,i))
                     ensum = enav +ensum
                     d(i,i,j) = d(i,i,j)+wrebut*drebut(j)
                end do
c
c               modify the RLW xchie, xchii and di for the
c               flow shear suppression effect if ifsflag > 0
c
                if (wrebut .gt. 0.0) then
                    xchierl(j) = xchierl(j)*coefa(j)
                    drebut(j)  = drebut(j)*coefc(j)
c
c --- the momentum diffusivity is assumed proportional to the ion thermal
c --- conductivity
c
                    xkangrl(j) = xchiirl(j) * coefd(j)
                    xchiirl(j) = xchiirl(j) * coefb(j)
                end if
c
c --- relax just the anomalous part of the diffusivity (i.e., without the
c --- neoclassical term) note that at this point xchierl and xchiirl
c --- are just the anomalous part returned by subroutine rebutlal.
c --- similarly, if itercorct > 1, then xchierlsv and xchiirlsv are
c --- just the anomalous parts from the previous iteration.
c
                if (itercorct .gt. 1) then
                  xchierl(j) = (1.0 - relaxrebut)*xchierlsv(j) +
     .                                relaxrebut *xchierl(j)
                  xchiirl(j) = (1.0 - relaxrebut)*xchiirlsv(j) +
     .                                relaxrebut *xchiirl(j)
                  xkangrl(j) = (1.0 - relaxrebut)*xkangrlsv(j) +
     .                                relaxrebut *xkangrl(j)
                end if
c
                xchierlsv(j)       = xchierl(j)
                xchiirlsv(j)       = xchiirl(j)
                xkerebut(j)        = xchierl(j)*enea
                xkangrlsv(j)       = xkangrl(j)
                d(nion+1,nion+1,j) = d(nion+1,nion+1,j)
     .                             + wrebut*xkerebut(j)
                xkirebut(j)        = xchiirl(j)*ensum
                d(nion+2,nion+2,j) = d(nion+2,nion+2,j)
     .                             + wrebut*xkirebut(j)
                if (iangrot .ne. 0)
     .          d(nion+4,nion+4,j) = d(nion+4,nion+4,j)
     .                          +wrebut*xkangrl(j)*dmassden(j)*rmajor**2
                xkerebut(j)        = xkerebut(j) + xkeneo(j)
                xkirebut(j)        = xkirebut(j) + xkineo(j)
                xchierl(j)         = xchierl(j)  + xkeneo(j)/enea
                xchiirl(j)         = xchiirl(j)  + xkineo(j)/ensum
                if (teprime .eq. 0.0 .and. tiprime .eq. 0.0)  go to 101
                xeffrl(j) = (enea*xchierl(j)*teprime+ensum*xchiirl(j)
     .                      *tiprime)/(enea*teprime+ensum*tiprime)
                go to 102
  101           if (j .eq. 1) then
                    xeffrl(j) = 0.0
               else
                    xeffrl(j) = xeffrl(j-1)
               end if
  102          continue
            end do
c
c --------------- HSIEH MODEL OF KE ------------------------------------
c --- NOTE THAT THE PARTICLE DIFFUSIVITY DUE TO THE HSIEH MODEL,
c --- DSHAY(J) IS CURRENTLY ZEROED, SINCE WE DON'T KNOW WHAT IT
c --- SHOULD BE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c ----------------------------------------------------------------------
c
      if (wshay .ne. 0.0) then
          timeshay = time
          call shay_chie (te, ene, bp, r, ra, dr, rmajor, rmbpav, nj,
     .                    smult, snexp, sbpexp, srexp, sbigrexp, stexp,
     .                    sdtdrexp, srin, srout, suserho, schie, schii,
     .                    skimult, irinshay, iroutshay, nout, ncrt,
     .                    sdenscale, ishayform, skimass)
          wmults = 1.0
c
c         if scsmult is true then we just calculatd the shape
c         factor in subroutine SHAY_CHIE (because smult is set
c         to 1.0 in this case). since the schie and schii are not
c         normalized we don't add them to d
c
          if (scsmult)  wmults = 0.0
            do j=1,nj-1               ! loop over the "half" grid
                dshay(j) = 0.0
                enea     = (ene(j+1)+ene(j))*0.5
                teprime  = (te(j+1)-te(j))/dr(j)
                tiprime  = (ti(j+1)-ti(j))/dr(j)
                ensum    = 0.0
c
c --- particle transport
c
                do i=1,nprim    ! have to use nprim to be consistent
c               do i=1,nion     ! with calculation of xkiinv (in PSOURC)
                  enav     = 0.5 * (en(j,i) + en(j+1,i))
                  ensum    = enav + ensum
                  d(i,i,j) = d(i,i,j)+wshay*dshay(j)
c                  dmassden = dmassden+atw(i)*xmassp*enav
                end do
c
c               modify the Hsieh schie, schii, and di, for the
c               flow shear suppression effect if ifsflag > 0
c
                schie(j) = schie(j) * coefa(j)
c
c --- the momentum diffusivity is assumed proportional to the ion thermal
c --- conductivity
c
                xkangsch(j) = coefd(j)*schii(j)
                schii(j)    = schii(j)*coefb(j)
                dshay(j)    = dshay(j)*coefc(j)
c
c --- relax just the anomalous part of the diffusivity (i.e., without
c --- the neoclassical term) note that at this point schie and schii
c --- are just the anomalous part returned by subroutine SHAY_CHIE.
c --- similarly, if itercorct > 1, then schiesv and schiisv are
c --- just the anomalous parts from the previous iteration
c
                if (itercorct .gt. 1) then
                  schie(j)    = (1.0 - relaxshay)*schiesv(j) +
     .                                 relaxshay *schie(j)
                  schii(j)    = (1.0 - relaxshay)*schiisv(j) +
     .                                 relaxshay *schii(j)
                  xkangsch(j) = (1.0 - relaxshay)*xkangschsv(j) +
     .                                 relaxshay *xkangsch(j)
                end if
c
                schiesv(j)         = schie(j)
                schiisv(j)         = schii(j)
                xkangschsv(j)      = xkangsch(j)
                xkeshay(j)         = schie(j) * enea
                d(nion+1,nion+1,j) = d(nion+1,nion+1,j)+wmults*wshay*
     .                                                       xkeshay(j)
                xkishay(j)         = schii(j)*ensum
                d(nion+2,nion+2,j) = d(nion+2,nion+2,j)+wmults*wshay*
     .                                                       xkishay(j)
                if (iangrot .ne. 0)
     .          d(nion+4,nion+4,j) = d(nion+4,nion+4,j)
     .            + wmults*wshay*coefd(j)*schii(j)*dmassden(j)*rmajor**2
                xkeshay(j)         = xkeshay(j) + wmults*xkeneo(j)
                xkishay(j)         = xkishay(j) + wmults*xkineo(j)
                schie(j)           = schie(j) + wmults*xkeneo(j)/enea
                schii(j)           = schii(j)+wmults*xkineo(j)/ensum
                if (teprime .eq. 0.0 .and. tiprime .eq. 0.0)  go to 103
                scheff(j) = wmults*(enea*schie(j)*teprime+ensum*schii(j)
     .                            *tiprime)/(enea*teprime+ensum*tiprime)
                go to 105
  103           if (j .eq. 1) then
                    scheff(j) = 0.0
                else
                    scheff(j) = scheff(j-1)
                end if
  105           continue
            end do
      end if
c
c ----------------------------------------------------------------------
c --- Weiland,Nordman model:
c --- weiland_setup returns xchie_weiland,xchii_weiland ,d_weiland on
c --- the half grid, in cm**2/sec,through include file weiland.i
c --- d(i,i,j) is in 1.0/(cm sec) for i = nion+1,nion+2 (i.e., te,ti)
c --- d(i,i,j) is in cm**2/sec for i = 1,2..nion
c --- these results will affect the calculations only if wweiland .ne. 0.0
c --- and the appropriate simulation mode is set. we calculate them even
c --- if wweiland = 0.0 however so that an appropriate printout and plot
c --- can be made. time_weiland is used to indicate the time when
c --- subroutine was last called.
c --- there are two single fluid chi effectives. the first, xeffwl,
c --- calculated below, is based on the Weiland  models of chie and chii.
c --- the second, xeffsf, cannot be calculated here since we may want
c --- to get xeffsf from experimental rather than theoretical data.
c --- thus we have to calculate xeffsf in subroutine PSOURC
c ----------------------------------------------------------------------
c
      if (include_weiland .eq. 1) then
c
c --- constants:
c
            time_weiland=time
            k1=1               ! NOTE may need to change HSJ
            azin_exp   = atw(nprim+1) / atw(k1)
            arho_exp   =  r(nj) / 100.0
            rmajor_exp = rmajor / 100.0
            do j=1,nj-1
c
                call weiland_setup (j, azin_exp, arho_exp, rmajor_exp)
c
                enea = (ene(j+1)+ene(j))*0.5
                teprime = (te(j+1)-te(j))/dr(j)
                tiprime = (ti(j+1)-ti(j))/dr(j)
                ensum = 0.0
c                dmassden = 0.0
                do i=1,nprim    ! have to use nprim to be consistent
****            do i=1,nion     ! with calculation of xkiinv (in PSOURC)
                  enav = 0.5 * (en(j,i)+en(j+1,i))
                  ensum = enav +ensum
                  d(i,i,j) = d(i,i,j)+wweiland*d_weiland(j)
c                  dmassden = dmassden+atw(i)*xmassp*enav
                end do
c
c               the momentum diffusivity is assumed proportional to the
c               ion thermal conductivity
c
                xkangwl(j) = xchii_weiland(j)
c
c               modify the Weiland xchie, xchii and di for the
c               flow shear suppression effect if ifsflag > 0
c               (if ifsflag =0 the coefficients are all unity)
c
                if (wweiland .gt. 0.0) then
                  xchie_weiland(j) = xchie_weiland(j)*coefa(j)
                  xkangwl(j)       = xkangwl(j) * coefd(j)
                  xchii_weiland(j) = xchii_weiland(j) * coefb(j)
                end if
c
c --- relax just the anomalous part of the diffusivity (i.e., without the
c --- neoclassical term) note that at this point xchie_weiland and xchii_weiland
c --- are just the anomalous part returned by subroutine weiland_setup.
c --- similarly, if itercorct > 1, then xchie_weilandsv and xchii_weilandsv are
c --- just the anomalous parts from the previous iteration.
c
              if (itercorct .gt. 1) then
                xchie_weiland(j) =(1.0 - relaxrebut)*xchie_weilandsv(j)+
     .                                   relaxrebut *xchie_weiland(j)
                xchii_weiland(j) =(1.0 - relaxrebut)*xchii_weilandsv(j)+
     .                                   relaxrebut *xchii_weiland(j)
                xkangwl(j)       =(1.0 - relaxrebut)*xkangwlsv(j) +
     .                                   relaxrebut *xkangwl(j)
              end if
c
                xchie_weilandsv(j)   = xchie_weiland(j)
                xchii_weilandsv(j)   = xchii_weiland(j)
                xke_weiland(j)       = xchie_weiland(j)*enea
                xki_weiland(j)       = xchii_weiland(j)*ensum
                xkangwlsv(j)         = xkangwl(j)
c
                d(nion+1,nion+1,j)   = d(nion+1,nion+1,j)           ! te
     .                                   + wweiland*xke_weiland(j)
                xketot(j)            = d(nion+1,nion+1,j)
                xchietot(j)          = xketot(j)/enea
c
                d(nion+2,nion+2,j)   = d(nion+2,nion+2,j)           ! ti
     .                                   + wweiland*xki_weiland(j)
                xkitot(j)            = d(nion+2,nion+2,j)
                xchiitot(j)          = xkitot(j)/ensum
c
                if (iangrot .ne. 0)
     .          d(nion+4,nion+4,j) = d(nion+4,nion+4,j) ! tord. rot.
     .                      +wweiland*xkangwl(j)*dmassden(j)*rmajor**2
                xkangtot(j)        = d(nion+4,nion+4,j)
                xdchitot(j)        = d(1,1,j)                    ! ni
c
                if (teprime .eq. 0.0 .and. tiprime .eq. 0.0) then
                    if (j .eq. 1) then
                        xeffwl(j) = 0.0
                    else
                        xeffwl(j) = xeffwl(j-1)
                    end if
                else
                    xeffwl(j) = (xketot(j)*teprime+xkitot(j)*tiprime)/
     .                                    (enea*teprime+ensum*tiprime)
                end if
            end do                   ! end loop over r(j), Weiland model
        end if                       ! end include_weiland branch
c 
c
c --------------------------------------------------------------------
        paleo:  if(include_paleo .gt. 0)then
                    call paleocl_tr(0,nion)
                endif paleo

c  --------------------------------------------------------------------
        if (include_ifs .eq. 1) then ! ======== BEGIN include_ifs branch
c
c           this model is active when Lti < Lticrit
c                           (i.e., dTi/dr > dTi/drcrit).
c           get the required scale lengths on the full grid
c           (input to IP_CHI2 is calculated on the full grid and
c           converted to half grid for use. Output from IP_CHI2 is
c           calculated on the half grid and converted to the full
c           grid for storage in ifs.i below)
c
            time_ifs = time
c
c           rhod_ifs is normalized, rho defined as d/D
c
            call scale_length (ti, rhod_ifs, nj, rlt_ifs)
c
c           to get nondimensional quantities,
c           multiply by rmajor/rhod_max_ifs
c
            rmajor_n = rmajor/rhod_max_ifs
            call multpl1 (rlt_ifs,nj,rmajor_n)
            call scale_length (ene,rhod_ifs,nj,rlne_ifs)
            call multpl1 (rlne_ifs,nj,rmajor_n)
            call scale_length (aabsq,rhod_ifs,nj,shat_ifs)
            call vecprod (shat_ifs,rhod_ifs,shat_ifs,nj,0)
            call multpl1 (shat_ifs, nj, -one) ! used to be -1.0
            do j=1,nj
c
c              thermal zeff as required by IFS model:
c              assume impurity is carbon
c
               zth_ifs(j) = (en(j,1)+36.*en(j,nprim+1))/
     .                                        (ene(j)-enbeam(j))
               zth_ifs(j) = zeff(j) ! to match Waltz
c
c              collisionality parameter:
c
               gnu_ifs(j) = 2.5e-15*ene(j)*rmajor/
     .                      ((te(j)**1.5)*(ti(j)**.5))
c
c              mass average thermal ion speed:
c
               amassth = 0.0
               tot_den = 0.0
               do jj=1,nprim       ! sum over hydrogenic species
                 if (atw(jj) .le. 3.1) then
                      amassth = amassth+atw(jj)*en(j,jj)*xmassp
                      tot_den = tot_den+en(j,jj)
                 end if
               end do
               amassth = amassth/tot_den   ! in grams
               vthi_ifs(j) = SQRT (ti(j) * 1.602e-9 / amassth)  ! CM/SEC
****           b_tot       = SQRT (btor**2+bp(j)**2)
               b_tot       = ABS (btor)
               omci_ifs(j) = charge*b_tot/(amassth*cee)
               tot_primary_ion_den(j) = tot_den
               tot_thermal_ion_den(j) = tot_den
               do jj=nprim+1,nion
                 tot_thermal_ion_den(j) =tot_thermal_ion_den(j)+en(j,jj)
               end do
            end do
c
c           scale length for ions,uses total primary ion density:
c
            call scale_length (tot_primary_ion_den,rhod_ifs,nj,rln_ifs)
            call multpl1 (rln_ifs,nj,rmajor_n)
            if (exbmult_ifs .ne. 0.0)
     .         call csperp  ! new from gary 12/15/98
            do j=1,nj-1     ! loop over the half grid
               rhodl_ifs = 0.5*(rhod_ifs(j+1)+rhod_ifs(j))
               rltl_ifs  = 0.5*(rlt_ifs(j)+rlt_ifs(j+1))
               rlnl_ifs  = 0.5*(rln_ifs(j)+rln_ifs(j+1))
               rlnel_ifs = 0.5*(rlne_ifs(j)+rlne_ifs(j+1))
c
c              CURRENT IFS NOT VALID FOR NEGATIVE (OR WEAK) SHEAR:
c              message from Dorland says use sthat instead of shat
c              see IP_CHI2 routine:
c
               shatl_ifs  = 0.5*(shat_ifs(j+1) + shat_ifs(j))
               taul_ifs   = (ti(j+1)+ti(j))/(te(j+1)+te(j))
               epsl_ifs   = rhodl_ifs*rhod_max_ifs/rmajor
               enbl_ifs   = (enbeam(j+1)+enbeam(j))/(ene(j+1)+ene(j))
               elongl_ifs = 0.5*(elong_r(j+1)+elong_r(j))
               zthl_ifs   = 0.5*(zth_ifs(j+1)+zth_ifs(j))
               gnul_ifs   = 0.5*(gnu_ifs(j+1)+gnu_ifs(j))
****           gperpl_ifs = 0.0   ! take from Waltz paper ?
               gperpl_ifs = sperp(j)*exbmult_ifs  ! Gary S. 12/15/98
               vthil_ifs  = 0.5*(vthi_ifs(j+1)+vthi_ifs(j))
               rho_il_ifs = 2.0 * vthil_ifs ! avg thermal ion gyroradius
     .                        / (omci_ifs(j+1)+omci_ifs(j))
               ql_ifs     = 0.5*(aabsq(j+1)+aabsq(j))
               grho1l_ifs =
     .            0.5*(grho1_ifs(j+1)+grho1_ifs(j)) !  <grad rho>
               grho2l_ifs =
     .            0.5*(grho2_ifs(j+1)+grho2_ifs(j)) ! <(grad rho)**2>
               drho_drhod =
     .            (r(j+1)-r(j))/(rhod_ifs(j+1)-rhod_ifs(j)) ! drho/drhod
               drho_drhod = drho_drhod/rhod_max_ifs
               itest_chi_ifs = 0
c
c              Call the IFS model.
c              Returned values are the last 7 arguments in CGS units.
c
               call ip_chi2 (itest_chi_ifs, rltl_ifs, rlnl_ifs,
     .                       rlnel_ifs, ql_ifs, elongl_ifs,
     .                       shatl_ifs, zthl_ifs, enbl_ifs,
     .                       taul_ifs, epsl_ifs, gnul_ifs,
     .                       gperpl_ifs, rmajor, rho_il_ifs,
     .                       vthil_ifs,                   ! end of input
     .                       RLTcrit_ifs(j), RLTcritz_ifs(j),
     .                       chi_0_ifs, g_ifs(j),
     .                       ip_chi2_counter,
     .                       gamma_ifs, chi_i_ifs(j), chi_e_ifs(j))
c
****  write (ncrt, '(a, i5)') ' IP_CHI2_COUNTER =', ip_chi2_counter
c
               chi_i_ifs(j) = chi_i_ifs(j) * grho1l_ifs * drho_drhod
               chi_e_ifs(j) = chi_e_ifs(j) * grho1l_ifs * drho_drhod
c
c              get results into d for simulation purposes:
c
c               NOTE      MY BIG ASSUMPTIONS ... HSJ !!!!!!!
c???????????????????????????????????????????????????????????????????????
c
                d_ifs(j)     = chi_i_ifs(j)
                xkang_ifs(j) = chi_i_ifs(j)
c
c               modify the chie, chii, and di, for the
c               flow shear suppression effect if ifsflag > 0
c
                chi_e_ifs(j) = chi_e_ifs(j) * coefa(j) ! = 1.0 if..
c                                              ..cstaebler is not called
                xkang_ifs(j) = xkang_ifs(j) * coefd(j)
                chi_i_ifs(j) = chi_i_ifs(j) * coefb(j)
                d_ifs(j)     = d_ifs(j)     * coefc(j)
c
c???????????????????????????????????????????????????????????????????????
c               END OF MY ASSUMPTIONS   HSJ   !!!!!!!
c
c --- relax just the anomalous part of the diffusivity (i.e., without the
c --- neoclassical term) note that at this point chi_e_ifs and chi_i_ifs
c --- are just the anomalous part returned by subroutine IP_CHI2.
c --- we use the Weiland storage vectors for the IFS model here (since
c --- there is no conflict,only one of the models is active in a given run)
c --- If itercorct > 1, then xchie_weilandsv and xchii_weilandsv are
c --- just the anomalous parts from the previous iteration.
c
               if (itercorct .gt. 1) then
                 d_ifs(j)     = (1.0 - relaxrebut) *xchii_weilandsv(j) +
     .                                 relaxrebut *d_ifs(j)
                 chi_e_ifs(j) = (1.0 - relaxrebut) *xchie_weilandsv(j) +
     .                                 relaxrebut *chi_e_ifs(j)
                 chi_i_ifs(j) = (1.0 - relaxrebut) *xchii_weilandsv(j) +
     .                                 relaxrebut *chi_i_ifs(j)
                 xkang_ifs(j) = (1.0 - relaxrebut) *xkangwlsv(j)       +
     .                                 relaxrebut *xkang_ifs(j)
               end if
c
               xchie_weilandsv(j) = chi_e_ifs(j)
               xchii_weilandsv(j) = chi_i_ifs(j)
               xkangwlsv(j)       = xkang_ifs(j)
c
               enea     = 0.5 * (ene(j+1)+ene(j))
               ensum    = 0.0
c               dmassden = 0.0
               do i=1,nprim    ! have to use nprim to be consistent
****           do i=1,nion     ! with calculation of xkiinv (in PSOURC)
                 enav     = 0.5 * (en(j,i)+en(j+1,i))
                 ensum    = enav +ensum
                 d(i,i,j) = d(i,i,j)+dorl_kotch*d_ifs(j)
c                 dmassden = dmassden+atw(i)*xmassp*enav
               end do
c
               xke_ifs(j)           = chi_e_ifs(j)*0.5*(ene(j+1)+
     .                                                  ene(j))
               d(nion+1,nion+1,j)   = d(nion+1,nion+1,j)           ! te
     .                                   + dorl_kotche*xke_ifs(j)
c
c              to check with ron's code:
c
****           d(nion+1,nion+1,j)   = dwifs(j)           ! te
**** .                                   + dorl_kotche*xke_ifs(j)
c
****           xki_ifs(j)           = chi_i_ifs(j)*0.5*(en(j+1,1)+
**** .                                                   en(j,1))
c
               xki_ifs(j)           = chi_i_ifs(j) * 0.5 *
     .           (tot_primary_ion_den(j) + tot_primary_ion_den(j+1))
               d(nion+2,nion+2,j)   = d(nion+2,nion+2,j)           ! ti
     .                                   + dorl_kotchi*xki_ifs(j)
               if (iangrot .ne. 0)
     .           d(nion+4,nion+4,j) = d(nion+4,nion+4,j) ! tord. rot.
     .                  +dorl_kotchi*xkang_ifs(j)*dmassden(j)*rmajor**2
c
               xdchitot(j)          = d(1,1,j)
               xketot(j)            = d(nion+1,nion+1,j)
               xchietot(j)          = xketot(j)/enea
               xkitot(j)            = d(nion+2,nion+2,j)
               xchiitot(j)          = xkitot(j)/ensum
               xkangtot(j)          = d(nion+4,nion+4,j)
c
            end do ! end loop over r(j) half grid
c
        end if ! ================================ END include_ifs branch

c---------------------------------------------------------------------------
c     GLF23model JEK 7/18/00
c     renormalized version from NTCC web site 05/01/02  HSJ
c     GLF23, Feb 19,2003 version has temporary corrections
c     for negative shear ( uses alpha_e =1.35, ifs =1)
c     Ref. is Waltz et. al. , Phys Plas, Vol 4. No. 7,July `997,2484
c---------------------------------------------------------------------------
c

      glf_sect:    if(include_glf .eq. 1 )then
c          write(6,*) 'Calculating transport using GLF23 model ...'
c          call radius_calc done in tport HSJ
          irotstab_save = irotstab
          jshoot=0
          jmm=0          !use diff_m quantities,full grid calc
          jmaxm=nj-1     ! grid size (0 to jmaxm in glf)
          igrad=0        !set to 1 for dv method below
c ----    set transport flags consistent with onetwo HSJ 3/15/01:
          itport_pt(1)= itran(1)
          itport_pt(2)= itran(nion+1)
          itport_pt(3)= itran(nion+2)
          itport_pt(4)= itran(nion+4)  
          itport_pt(5)= 0 !vtheta transport
          if (iglf_eq .eq. 1) then !jmp.den end
            itport_pt(1)= 1 !itran(1)
            itport_pt(2)= 1 !itran(nion+1)
            itport_pt(3)= 1 !itran(nion+2)
            itport_pt(4)= 1 !itran(nion+4)
          end if !jmp.den end    
          if (ivphi_glf.gt.0) then
c             write(6,*) 'Transporting vphi using GLF23 model ...'
             itport_pt(4)=1
          endif
          alphae_glf = exbmult_glf   !note exbmult_glf = alpha_e in glf model
          alphad_glf=1.          ! include diamagnetic terms in gamma_E
          rmajor_glf=rmajor/100. ! major radius, R (m)
          rminor_glf=rminor/100. ! minor radius, a (m)
c          btor_exp=-btor/1.e4   ! toroidal field, Bt (T)
          btor_exp= btor/1.e4    ! toroidal field, Bt (T)
          arho_iglf=r(nj)/100.   ! toroidal flux at LCFS, rho (m)
          k = 1                  ! use first primary ion ??
          amassgas_glf=2.0       ! mass of main ions
          amassimp_glf=12.0      ! mass of impurity
          zimp_glf=6.0           ! charge of impurity
c ----    set above consistent with onetwo HSJ 3/15/01:
c          zimp_glf = z(1,nprim+1) ! cant pass in profile
          amassgas_glf = atw(k)
c ----    amassgas_glf should be constant over grid points in glf23 module
c         if it is not constant, the average value is used. 
          if(nprim .gt. 1) then !JMP BLOCK START
             if ((nprim .eq. 2).and.
     .           ((inenez .eq. 1).or.(inenez .lt. -1))) then
                amassgas_glf = atw(1)*zfrac + atw(2)*(1.0 - zfrac)
             else
                amassgas_glf = 0.0
                do j = 2,nj
                   amassgas_tmp = 0.0
                   tot_den = 0.0
                   do jj=1,nprim
                      if (atw(jj) .le. 3.1) then
                         amassgas_tmp = amassgas_tmp + atw(jj)*en(j,jj)
                         tot_den=tot_den+en(j,jj)
                      end if
                   end do
                   amassgas_tmp = amassgas_tmp/tot_den
                   amassgas_glf = amassgas_glf + amassgas_tmp
                end do   
                amassgas_glf = amassgas_glf/(nj-1.0)      
c               print *,nprim,inenez,amassgas_glf
c               pause
             end if   
c            print *,nprim,inenez,amassgas_glf
c            pause
          end if !JMP BLOCK END

c         check to make sure density gradient scale length is not zero:
c          (a zero value will crash glf2d)
c             !$OMP PARALLEL DO PRIVATE(diff)
            te_glf(1)=te(1)
            ti_glf(1)=ti(1)
          do j = 2,nj
             te_glf(j)=te(j)
             ti_glf(j)=ti(j)
             diff = en(j,k)-en(j-1,k)
             if(abs(diff) .lt. 1.e-8)
     .                en(j,k)= 0.99999*en(j,k)
             diff = ene(j)-ene(j-1)
             if(abs(diff) .lt. 1.e-8)
     .                ene(j)= 0.99999*ene(j)
             diff = te_glf(j)-te_glf(j-1)
             if(abs(diff) .lt. 1.e-8)
     .                te_glf(j)= 0.99999*te_glf(j)

             diff = ti_glf(j)-ti_glf(j-1)
             if(abs(diff) .lt. 1.e-8)
     .                ti_glf(j)= 0.99999*ti_glf(j)

           enddo
c
c..compute geometric quantities, densities, temperatures, beta, q
c  note: ni is primary ion density
c
c         !$OMP PARALLEL DO PRIVATE(amassth,tot_den)
          do j=1,nj
            rho_glf(j)=max( r(j)/(100.*arho_iglf), 1.D-6 ) ! rho_glf =[0,1] !jmp.ibm
            tau_glf(j)=ti_glf(j)/te_glf(j)
            ne_glf(j)=1.e-13*ene(j)
            !JMP ni_glf(j)=1.e-13*en(j,k)
            ni_glf(j)=0.0 !JMP BLOCK START
            do jj=1,nprim
              if (atw(jj) .le. 3.1) then
                 ni_glf(j)=ni_glf(j)+en(j,jj)
              end if
            end do
            ni_glf(j)    = 1.e-13*ni_glf(j) !JMP BLOCK END
            ns_glf(j)    = 1.e-13*enbeam(j) !beam density
            zeff_exp(j)  = zeff(j)
            betae_glf(j) = 400.*ne_glf(j)*te_glf(j)/(1.e5*btor_exp**2)
            betai_glf(j) = 400.*ni_glf(j)*ti_glf(j)/(1.e5*btor_exp**2)
            beta_glf(j)  = betae_glf(j)+betai_glf(j)
c            q_glf(j)=-q(j)
            q_glf(j)=ABS(q(j))  !HSJ sign of q bas
            rmin_glf(j) = max( rminor_r(j)/100., 1.D-6) !jmp.ibm
            rmaj_glf(j) = ravg_r(j)/100.
c            grho1_glf(j) = grho1_ifs(j)*100./r(nj) !  <grad rho>/a in m^-1
c            grho2_glf(j) = grho2_ifs(j)*(100./r(nj))**2 !<(grad rho)**2>/a**2 in m^-2
            grho1_glf(j) = grho1_ifs(j)
            grho2_glf(j) = grho2_ifs(j)
            elong_glf(j) = kappa_nj(j)
            amassth=0.0 
            tot_den=0.0
            do jj=1,nprim       ! sum over hydrogenic species
              if (atw(jj) .le. 3.1) then
                 amassth=amassth+atw(jj)*en(j,jj)*xmassp
                 tot_den=tot_den+en(j,jj)
              end if
            end do
            amassth=amassth/tot_den   ! in grams
            tot_primary_ion_den(j)=tot_den
            tot_thermal_ion_den(j)=tot_den
            do jj=nprim+1,nion
              tot_thermal_ion_den(j)=tot_thermal_ion_den(j)+en(j,jj)
            end do
          enddo
c         !$OMP PARALLEL DO
          do j=2,nj
            drhodr(j)=(rho_glf(j)-rho_glf(j-1))*arho_iglf/
     .         (rmin_glf(j)-rmin_glf(j-1)+1.e-6)
            drhodrrrho(j)=drhodr(j)*rmin_glf(j)/arho_iglf/rho_glf(j)
          enddo
          drhodr(1)=drhodr(2)
          drhodrrrho(1)=drhodrrrho(2)
c
c..compute csda and rhosda using effective B-field
c  if ibtflag.gt.0 (default)
c
c         !$OMP PARALLEL DO
          do j=1,nj
            csda_glf(j)=9.79e5*(te_glf(j)*1.e3)**.5/
     .         (arho_iglf*100.)/amassgas_glf**.5             ! 1/sec
            rhosda_glf(j)=((1.02e2*(te_glf(j)*1.e3)**.5)/
     .         btor_exp/1.e4)*(amassgas_glf**.5)/(arho_iglf*100.)
          enddo
          if(ibtflag_glf.gt.0) then
c          !$OMP PARALLEL DO
            do j=1,nj
              bteff_glf(j)=btor_exp*rho_glf(j)*arho_iglf/
     .                   rmin_glf(j)*drhodr(j)
              rhosda_glf(j)=((1.02e2*(te_glf(j)*1.e3)**.5)/
     .         bteff_glf(j)/1.e4)*(amassgas_glf**.5)/(arho_iglf*100.)
            enddo
          endif
c
c..compute normalized gradients, shear, alpha
c
c         !$OMP PARALLEL DO PRIVATE(dqdrho,drhodq)
          do j=2,nj
            drhogf(j)=rho_glf(j)-rho_glf(j-1)  !rho_glf is [0,1]
            !forward difference on derivatives to get match with xptor: 
c            zpte_glf(j-1)=-(log(te_glf(j))-log(te_glf(j-1)))/drhogf(j) ! a/LTe
c            zpti_glf(j-1)=-(log(ti_glf(j))-log(ti_glf(j-1)))/drhogf(j) ! a/LTi
c            zpne_glf(j-1)=-(log(ne_glf(j))-log(ne_glf(j-1)))/drhogf(j) ! a/Lne
c            zpni_glf(j-1)=-(log(ni_glf(j))-log(ni_glf(j-1)))/drhogf(j) ! a/Lni

            tea = 0.5*(te(j)+te(j-1))
            tia = 0.5*(ti(j)+ti(j-1))
            ane = 0.5*(ne_glf(j)+ne_glf(j-1)) 
            ani = 0.5*(ni_glf(j)+ni_glf(j-1)) 
            zpte_glf(j-1)=(-te_glf(j)+te_glf(j-1))/(tea*drhogf(j)) ! 1/LTe
            zpti_glf(j-1)=(-ti_glf(j)+ti_glf(j-1))/(tia*drhogf(j)) ! 1/LTi
            zpne_glf(j-1)=(-ne_glf(j)+ne_glf(j-1))/(ane*drhogf(j)) ! 1/Lne
            zpni_glf(j-1)=(-ni_glf(j)+ni_glf(j-1))/(ani*drhogf(j)) ! 1/Lni
            dqdrho=-(q(j)-q(j-1))/drhogf(j)  
            drhodq=-(rho_glf(j)+rho_glf(j-1))/(q(j)+q(j-1)) 
            shat_exp(j)=drhodq*dqdrho      !(r/q)dq/dr
            if (abs(shat_exp(j)).lt.1.e-6) shat_exp(j)=1.e-6
            alpha_exp(j)=-drhodr(j)*q_glf(j)**2*rmajor_glf/arho_iglf*
     .                   (beta_glf(j)-beta_glf(j-1))/drhogf(j)
            if(freeze_alpha_exp .eq.0)then
                alpha_exp_save(j)=alpha_exp(j)
            else if(freeze_alpha_exp .eq. 1)then
                alpha_exp(j) = alpha_exp_save(j)
            else 
               if(shat_exp(j) .lt. 0.0)then
                 alpha_exp(j) = MIN(0.999D0,alpha_exp(j)) !jmp.ibm
                 if(freeze_alpha_exp .gt. 2)then
                 print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                 print *,'alpha_exp =1.0 s,r =',shat_exp(j),r(j)
                 endif
               endif
            endif
          enddo


          if(ddebug(14) .ge. 1.)then
             call copya(alpha_exp,glf_work,nj)  !to avoid compound smoothing
             jstart = ddebug(14) + 1
             jend = nj - ddebug(14)  !OK to do average past fix_edgete,etc
             jsteps = 2*ddebug(14) + 1
             do j=jstart,jend
                sum = 0.0d0
                jl = j - ddebug(14)
                jr = j + ddebug(14)
                do jj = jl,jr
                   sum = sum + glf_work(jj)
                enddo         
                alpha_exp(j) = sum/jsteps
             enddo
          endif

          if(ddebug(15) .ge. 1.)then
             call copya(shat_exp,glf_work,nj)  !to avoid compound smoothing
             jstart = ddebug(15) + 1
             jend = nj - ddebug(15)  !OK to do average past fix_edgete,etc
             jsteps = 2*ddebug(15) + 1
             do j=jstart,jend
                sum = 0.0d0
                jl = j - ddebug(15)
                jr = j + ddebug(15)
                do jj = jl,jr
                   sum = sum + glf_work(jj)
                enddo         
                shat_exp(j) = sum/jsteps
             enddo
          endif
 

          drhogf(1)=1.e-6
          if(freeze_alpha_exp .eq.0)then
              alpha_exp(1)=1.e-6
              alpha_exp_save(1)=alpha_exp(1)
           else if(freeze_alpha_exp .eq. 1)then
                alpha_exp(1) = alpha_exp_save(1)
           endif
 999      format(i2,2x,0p1f4.2,0p6f14.5)
c
c..compute plasma toroidal angular velocity (angrotp) in 1/s
c  and parallel (vpar), perpendicular(vper), and toroidal (vphi) 
c  velocities in m/s
c
c         !$OMP PARALLEL DO PRIVATE(fc,akappa1,akappa2,alpha_neo)
          do j=1,nj        ! loop over the full grid
            fc=1-1.46*(rmin_glf(j)/rmaj_glf(j))**0.5+
     .         0.46*(rmin_glf(j)/rmaj_glf(j))**1.5
            akappa1=0.8839*fc/(0.3477+0.4058*fc)
            akappa2=(1.-fc)/(1.+1.1671*fc)
            alpha_neo=-akappa1+1.   !alpha1neo in Waltz ref.
c
            fcglf(j)= fc
            ak1glf(j)=akappa1
            ak2glf(j)=akappa2
            aneoglf(j)=alpha_neo
c
            vstar_sign=-1.
            corot=1.
            pgeo_local=drhodr(j)
            rdrho_local=rmin_glf(j)/arho_iglf/rho_glf(j)
            rdrho_local_p1=rdrho_local
            !angular rotaton version of eq 16c in phys plas, 4,7,97,2486. HSJ:
            angrotp_exp(j)=angrot(j)+
     .         akappa2*3./2.*csda_glf(j)*zpti_glf(j)*
     .         tau_glf(j)/rho_glf(j)*q_glf(j)*
     .         rhosda_glf(j)*pgeo_local/rdrho_local
 
            if (angrot(j).eq.0.) angrotp_exp(j)=0.
            angrotp_exp(j)=corot*angrotp_exp(j)
            !determine Vexb (= ve)
            !orginal form (by J.Kinsey)
c            ve_glf(j)=-tau_glf(j)*csda_glf(j)*arho_iglf*
c     .         rhosda_glf(j)*(zpni_glf(j)+zpti_glf(j)+(1.-
c     .         rdrho_local*rho_glf(j)*arho_iglf/rmajor_glf/q_glf(j))*
c     .         (alpha_neo-1.)*zpti_glf(j))*vstar_sign*pgeo_local

            !set ve to value used in callglf2d  HSJ 05/0103:
            ve_glf(j)=-tau_glf(j)*csda_glf(j)*arho_iglf*
     .         rhosda_glf(j)*(zpni_glf(j) +
     .         alpha_neo*zpti_glf(j))*vstar_sign*pgeo_local
            vexb_den_grad(j) = -tau_glf(j)*csda_glf(j)*arho_iglf*
     .         rhosda_glf(j)*zpni_glf(j)*vstar_sign*pgeo_local
            vexb_ti_grad(j) =-tau_glf(j)*csda_glf(j)*arho_iglf*
     .         rhosda_glf(j)*zpti_glf(j)*pgeo_local*alpha_neo *
     .         vstar_sign -rdrho_local*rho_glf(j)*
     .         arho_iglf/rmajor_glf/q_glf(j)*
     .         rmajor_glf*(angrotp_exp(j) - angrot(j))
            vexb_w_grad(j) = -rdrho_local*rho_glf(j)*arho_iglf/
     .                  rmajor_glf/q_glf(j)*rmajor_glf*angrot(j) 
              !for either form of the above we now 
              !add sigma*(rminor/(Rmajor*q))Vtoroidal to Vexb.
              !Note that sigma is assumed  +1 here because the
              !toroidal rotation velocity is a signed quantitiy.
              !This term is identical to eq. 16a in Waltz,
              !phys plas, 4,7,97,2486. HSJ:
            ve_glf(j)=ve_glf(j)-rdrho_local*rho_glf(j)*
     .         arho_iglf/rmajor_glf/q_glf(j)*
     .         rmajor_glf*angrotp_exp(j) 
            vexb_tot(j) =  vexb_ti_grad(j) +vexb_den_grad(j) 
     .             + vexb_w_grad(j)
c            test = ve_glf(j) - vexb_tot(j)
c            print *,'test in diffuse =',test
            vparm_glf(j)=rmajor_glf*angrotp_exp(j)+
     .         vstar_sign*((1.-alpha_neo)*zpti_glf(j))*
     .         tau_glf(j)*csda_glf(j)*arho_iglf*
     .         rhosda_glf(j)*pgeo_local*rho_glf(j)*arho_iglf
     .         /rmajor_glf/q_glf(j)*rdrho_local
            vperm_glf(j)=-tau_glf(j)*csda_glf(j)*arho_iglf*
     .         rhosda_glf(j)*((1.-rdrho_local*rho_glf(j)*arho_iglf/
     .         rmajor_glf/q_glf(j))*(alpha_neo-1.)*zpti_glf(j))*
     .         vstar_sign*pgeo_local
            vperm_glf(j)=vperm_glf(j)-rdrho_local*rho_glf(j)*
     .         arho_iglf/rmajor_glf/q_glf(j)*
     .         rmajor_glf*angrotp_exp(j)
            vphim_glf(j)=rmajor_glf*angrotp_exp(j)
          enddo
c
c..compute ExB and parallel velocity shear rates
c
c         !$OMP PARALLEL DO PRIVATE(rdrho_local,rdrho_local_p1)
          do j=1,nj-1        
            rdrho_local=rmin_glf(j)/arho_iglf/rho_glf(j)
c            rdrho_local_p1=rdrho_local  
            rdrho_local_p1=rmin_glf(j+1)/arho_iglf/rho_glf(j+1)   !HSJ 4/30/03
            egamma_exp(j)=drhodrrrho(j)*(rho_glf(j)+rho_glf(j+1))/
     .         (q_glf(j)+q_glf(j+1))*(ve_glf(j+1)*q_glf(j+1)/
     .         rho_glf(j+1)/rdrho_local_p1-ve_glf(j)*q_glf(j)/
     .         rho_glf(j)/rdrho_local)/(rho_glf(j+1)-rho_glf(j))/
     .         arho_iglf/csda_glf(j)
c           individual componenets of egamma_exp:
            gamma_p_exp(j)=-drhodr(j)*(vparm_glf(j+1)-vparm_glf(j))/
     .         (rho_glf(j+1)-rho_glf(j))/arho_iglf/csda_glf(j)
            gamma_den(j)= drhodrrrho(j)*(rho_glf(j)+rho_glf(j+1))/
     .         (q_glf(j)+q_glf(j+1))*(vexb_den_grad(j+1)*q_glf(j+1)/
     .         rho_glf(j+1)/rdrho_local_p1-vexb_den_grad(j)*q_glf(j)/
     .         rho_glf(j)/rdrho_local)/(rho_glf(j+1)-rho_glf(j))/
     .         arho_iglf/csda_glf(j)
            gamma_ti(j)= drhodrrrho(j)*(rho_glf(j)+rho_glf(j+1))/
     .         (q_glf(j)+q_glf(j+1))*(vexb_ti_grad(j+1)*q_glf(j+1)/
     .         rho_glf(j+1)/rdrho_local_p1-vexb_ti_grad(j)*q_glf(j)/
     .         rho_glf(j)/rdrho_local)/(rho_glf(j+1)-rho_glf(j))/
     .         arho_iglf/csda_glf(j)
            gamma_w(j)= drhodrrrho(j)*(rho_glf(j)+rho_glf(j+1))/
     .         (q_glf(j)+q_glf(j+1))*(vexb_w_grad(j+1)*q_glf(j+1)/
     .         rho_glf(j+1)/rdrho_local_p1-vexb_w_grad(j)*q_glf(j)/
     .         rho_glf(j)/rdrho_local)/(rho_glf(j+1)-rho_glf(j))/
     .         arho_iglf/csda_glf(j)
c            test2 = egamma_exp(j) - gamma_den(j) - gamma_ti(j) 
c     .                                            - gamma_w(j)
c            print *,'test2 in diffuse =',test2
          enddo
          if(irotstab_save .eq. -1 .and. init_egamma_exp .eq. 0 )then
             egamma_exp_initial(:) = egamma_exp(:)
             init_egamma_exp = 1
          endif
          if(irotstab_save .eq. -1)
     .    egamma_exp(:) = egamma_exp_initial(:)

          if(ddebug(13) .ge. 1.)then
             call copya(egamma_exp,glf_work,nj)  !to avoid compound smoothing
             jstart = ddebug(13) + 1
             jend = nj-1 - ddebug(13)  !OK to do average past fix_edgete,etc
             jsteps = 2*ddebug(13) + 1
             do j=jstart,jend
                sum = 0.0d0
                jl = j - ddebug(13)
                jr = j + ddebug(13)
                do jj = jl,jr
                   sum = sum + glf_work(jj)
                enddo         
                egamma_exp(j) = sum/jsteps
             enddo
          endif
          if(irotstab_save .eq. -1)irotstab = 0 !passed to glf23 as 0
c
c..for impurity eqns (limp_glf > 0)
c
          idengrad_glf=2  ! simple dilution
          if(limp_glf.gt.0) then
            write(6,*) 'Turning on impurity dynamics in GLF23 model ...'
            jroot_iglf=12
            idengrad_glf=3
          endif
c

c
       
                 njm1 =nj-1
           if(write_glf_namelist .gt. 0)call writeglf2d(             
     &                               jeigen_iglf,jroot_iglf, glf23_iglf,       
     &                                  jshoot, jmm, jmaxm, itport_pt  , 
     &                                 irotstab, te_glf, ti_glf, ne_glf, 
     &                     ni_glf, ns_glf,igrad, idengrad_glf, zpte_glf,
     &                                     zpti_glf, zpne_glf, zpni_glf,
     &                             angrotp_exp, egamma_exp, gamma_p_exp,
     &                                  vphim_glf, vparm_glf, vperm_glf,
     &                        zeff_exp, btor_exp, ibtflag_glf, rho_glf,
     &                                  arho_iglf, grho1_glf, grho2_glf,
     &                         rmin_glf, rmaj_glf, rmajor_glf, zimp_glf, 
     &                         amassimp_glf, q_glf, shat_exp, alpha_exp,
     &                               elong_glf, amassgas_glf,alphae_glf,
     &                  x_alpha_glf, i_delay,lprint_glf,itport_glf,njm1)
            if(write_glf_namelist .eq. 2)
     &      call STOP('DIFFUS: write_glf_namelist specified  stop',0)
c
       if(freeze_xte .eq. 0 .or. freeze_xti .eq. 0 
     .                                      .or. freeze_xwr .eq. 0) then
             i_delay = 0                    ! no delay in ExB egamma_d

             if(glf23_ncpus  .eq. 1 )then     !single procesor version
               If(.not. dv_method)then

#ifndef PARALLEL_12
c single call to glf does entire interior grid from 2  to nj-1 ( which is 1 to nj-2
c in callglf2d where the grid runs from 0 to nj-1)
                          call callglf2d(           !inputs
     &                               jin,jout, !jmp.ibm.par
     &                               jeigen_iglf,jroot_iglf, glf23_iglf,         
     &                                  jshoot, jmm, jmaxm, itport_pt  , 
     &                                 irotstab, te_glf, ti_glf, ne_glf, 
     &                     ni_glf, ns_glf,igrad, idengrad_glf, zpte_glf,
     &                                     zpti_glf, zpne_glf, zpni_glf,
     &                             angrotp_exp, egamma_exp, gamma_p_exp,
     &                                  vphim_glf, vparm_glf, vperm_glf,
     &                        zeff_exp, btor_exp, ibtflag_glf, rho_glf,
     &                                  arho_iglf, grho1_glf, grho2_glf,
     &                         rmin_glf, rmaj_glf, rmajor_glf, zimp_glf, 
     &                         amassimp_glf, q_glf, shat_exp, alpha_exp,
     &                               elong_glf, amassgas_glf,alphae_glf,
     &                                             x_alpha_glf, i_delay,      !outputs:
     &                      diffnem, chietem,chie_etg_tem,chiitim, 
     &                      etaphim, etaparm, 
     &                      etaperm,exchm, diff_m, chie_m,chie_etg_m,
     &                      chii_m,etaphi_m, etapar_m, etaper_m,
     &                      exch_m, egamma_m, egamma_d, gamma_p_m,
     &                      anrate_m, anrate2_m, anfreq_m, anfreq2_m,
     &                      gamma_net_i,gamma_net_e  )
#else
!------------------------------------------------------------------------------
!
!
      !print *,'*** glf23 section entered'
      start_time = mpi_wtime()

      iflag_mpistop(1) = 0
      CALL mpi12_send_int(iflag_mpistop,1,99)
      !print *,'*** after first send'  

      iparam( 1) = jeigen_iglf    
      iparam( 2) = jroot_iglf     
      iparam( 3) = glf23_iglf    
      iparam( 4) = jshoot        
      iparam( 5) = jmm           
      iparam( 6) = jmaxm         
      iparam( 7) = itport_pt(1)  
      iparam( 8) = itport_pt(2)  
      iparam( 9) = itport_pt(3)  
      iparam(10) = itport_pt(4)  
      iparam(11) = itport_pt(5)  
      iparam(12) = irotstab      
      iparam(13) = igrad         
      iparam(14) = idengrad_glf  
      iparam(15) = ibtflag_glf   
      iparam(16) = i_delay       

      CALL mpi12_send_int(iparam,32,1)

      rparam( 1) = zpte_glf(1)  
      rparam( 2) = zpti_glf(1)  
      rparam( 3) = zpne_glf(1)  
      rparam( 4) = zpni_glf(1)  
      rparam( 5) = btor_exp     
      rparam( 7) = arho_iglf    
      rparam( 8) = rmajor_glf   
      rparam( 9) = zimp_glf     
      rparam(10) = amassimp_glf 
      rparam(11) = amassgas_glf 
      rparam(12) = alphae_glf !exbmult_glf  
      rparam(13) = x_alpha_glf  

      CALL mpi12_send_real(rparam,32,2)

      end_time = mpi_wtime()
      mpi_etimes(1) = mpi_etimes(1) + end_time-start_time

      start_time = mpi_wtime()

      CALL mpi12_d_send(te_glf     ,3)
      CALL mpi12_d_send(ti_glf     ,3)
      CALL mpi12_d_send(ne_glf     ,3)
      CALL mpi12_d_send(ni_glf     ,3)
      CALL mpi12_d_send(ns_glf     ,3)
      CALL mpi12_d_send(angrotp_exp,3)
      CALL mpi12_d_send(egamma_exp ,3)
      CALL mpi12_d_send(gamma_p_exp,3)
      CALL mpi12_d_send(vphim_glf  ,3)
      CALL mpi12_d_send(vparm_glf  ,3)
      CALL mpi12_d_send(vperm_glf  ,3)
      CALL mpi12_d_send(zeff_exp   ,3)
      CALL mpi12_d_send(rho_glf    ,3)  
      CALL mpi12_d_send(grho1_glf  ,3)
      CALL mpi12_d_send(grho2_glf  ,3)
      CALL mpi12_d_send(rmin_glf   ,3)
      CALL mpi12_d_send(rmaj_glf   ,3)
      CALL mpi12_d_send(q_glf      ,3)
      CALL mpi12_d_send(shat_exp   ,3)
      CALL mpi12_d_send(alpha_exp  ,3)
      CALL mpi12_d_send(elong_glf  ,3)

      !CALL MPI_BCAST(te_glf     ,51, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
      !CALL MPI_BCAST(ti_glf     ,51, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
      !CALL MPI_BCAST(ne_glf     ,51, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
      !CALL MPI_BCAST(ni_glf     ,51, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
      !CALL MPI_BCAST(ns_glf     ,51, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
      !CALL MPI_BCAST(angrotp_exp,51, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
      !CALL MPI_BCAST(egamma_exp ,51, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
      !CALL MPI_BCAST(gamma_p_exp,51, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
      !CALL MPI_BCAST(vphim_glf  ,51, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
      !CALL MPI_BCAST(vparm_glf  ,51, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
      !CALL MPI_BCAST(vperm_glf  ,51, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
      !CALL MPI_BCAST(zeff_exp   ,51, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
      !CALL MPI_BCAST(rho_glf    ,51, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
      !CALL MPI_BCAST(grho1_glf  ,51, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
      !CALL MPI_BCAST(grho2_glf  ,51, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
      !CALL MPI_BCAST(rmin_glf   ,51, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
      !CALL MPI_BCAST(rmaj_glf   ,51, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
      !CALL MPI_BCAST(q_glf      ,51, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
      !CALL MPI_BCAST(shat_exp   ,51, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
      !CALL MPI_BCAST(alpha_exp  ,51, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
      !CALL MPI_BCAST(elong_glf  ,51, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

      end_time = mpi_wtime()
      mpi_etimes(2) = mpi_etimes(2) + end_time-start_time

      !print *,'*** before callglf2d'  
      start_time = mpi_wtime()
      jin = 1
      jout = js_mpi12(myid+1)-2
      mpi_counts(1) = jout-jin+1
      jmm = -1
      CALL callglf2d(    ! INPUT------------------------------------
     .                   !
     .     jin,jout,
     .     jeigen_iglf,  ! 0
     .     jroot_iglf,   ! 8
     .     glf23_iglf,   ! 1
     .     jshoot,       ! 0
     .     jmm,          ! 0
     .     jmaxm,        ! nj-1
     .     itport_pt,    ! dimension 5
     .     irotstab,     ! 1
     .     te_glf,       
     .     ti_glf,       
     .     ne_glf,       
     .     ni_glf,       
     .     ns_glf,       
     .     igrad,        ! 0 
     .     idengrad_glf, ! 2
     .     zpte_glf(1),  ! not used since igrad = 0
     .     zpti_glf(1),  ! not used since igrad = 0
     .     zpne_glf(1),  ! not used since igrad = 0
     .     zpni_glf(1),  ! not used since igrad = 0
     .     angrotp_exp,  
     .     egamma_exp,   
     .     gamma_p_exp,  
     .     vphim_glf,    
     .     vparm_glf,    
     .     vperm_glf,    
     .     zeff_exp,     
     .     btor_exp,     ! constant
     .     ibtflag_glf,  ! 1
     .     rho_glf,      ! normalized
     .     arho_iglf,    ! r(nj)
     .     grho1_glf,    
     .     grho2_glf,    
     .     rmin_glf,     
     .     rmaj_glf,     
     .     rmajor_glf,   ! constant
     .     zimp_glf,     ! constant
     .     amassimp_glf, ! constant
     .     q_glf,        
     .     shat_exp,     
     .     alpha_exp,    
     .     elong_glf,    
     .     amassgas_glf, ! constant
     .     alphae_glf,   ! constant  exbmult_glf
     .     x_alpha_glf,  ! constant 
     .     i_delay,      ! 0
     .                   !
     .                   ! OUTPUT-----------------------------------
     .                   !                   
     .     diffnem,      ! constat not used
     .     chietem,      ! constat not used
     .     chie_etg_tem, ! constat not used
     .     chiitim,      ! constat not used
     .     etaphim,      ! constat not used
     .     etaparm,      ! constat not used
     .     etaperm,      ! constat not used
     .     exchm,        ! constat not used  
     .     diff_m,       ! not used ******************
     .     chie_m,       ! ***************************
     .     chie_etg_m,   ! ***************************
     .     chii_m,       ! ***************************
     .     etaphi_m,     ! ***************************
     .     etapar_m,     ! not used
     .     etaper_m,     ! not used
     .     exch_m,       ! not used 
     .     egamma_m,     ! not used
     .     egamma_d,     ! not used
     .     gamma_p_m,    ! not used
     .     anrate_m,     ! not used
     .     anrate2_m,    ! not used
     .     anfreq_m,     ! not used
     .     anfreq2_m,    ! not used
     .     gamma_net_i,  ! not used  
     .     gamma_net_e,  ! not used  
     .     )
      end_time = mpi_wtime()
      mpi_etimes(4) = mpi_etimes(4) + end_time-start_time
      glf_load_buf = end_time-start_time
      !print *,'*** after callglf2d'  

      if(iloadbal.eq.1) then
         CALL dynamic_balance(glf_load_buf)
         !print *,'*** after dynamic balance'
      endif

      start_time = mpi_wtime()
      CALL mpi12_g_recv(diff_m     ,4)
      CALL mpi12_g_recv(chie_m     ,4)
      CALL mpi12_g_recv(chii_m     ,4)
      CALL mpi12_g_recv(etaphi_m   ,4)
      end_time = mpi_wtime()
      mpi_etimes(3) = mpi_etimes(3) + end_time-start_time
      !print *,'*** after g_recv'        

      if(icalled .eq. 0) then
      !print *,'before print'
        open(71,file='glf_start.dat',status='unknown',iostat=ierr)
        write(71,*) js_mpi12(0:MAXPROC)
        do j=1,nj-1
          write(71,*) j,chii_m(j)
        end do 
        icalled = 1
        close(71)
      !print *,'end print'        
      end if  
!
!
!------------------------------------------------------------------------------
#endif

                    do j= 1,nj-1 
                       if(ichiv_chii .eq. 1) etaphi_m(j) = chii_m(j) !jmp.snu

!                       if(diff_m(j) .lt. 0.0)then !jmp.snu
!                            vpinch(j) = diff_m(j)*zpni_glf(j)/
!     .                                          arho_iglf             !m/sec
!     .                                          *100. + vpinch(j)
!                            diff_m(j) = 0.0
!                       endif

                       chie_etg_m(j) = 1.e4*chie_etg_m(j)              ! convert  to cm^2/sec
                       if(chie_m(j) .lt. 0.0)then                      ! includes chie_etg_m HSJ
                          if (no_te_convection .eq. 1)then
                            chie_m(j) = 0.0
                            ve_eff_m(j) =0.0
                          else
                            ve_eff_m(j) = chie_m(j)*zpte_glf(j)/
     .                                           arho_iglf             !m/sec
                            ve_eff_m(j)  = ve_eff_m(j)*100.
                            chie_m(j) = 0.0
                          endif
                       else
                           ve_eff_m(j) = 0.0
                       endif
                       if(chii_m(j) .lt. 0.0)then
                          if (no_ti_convection .eq. 1)then
                            chii_m(j) = 0.0
                            vi_eff_m(j) =0.0
                          else
                            vi_eff_m(j) = chii_m(j)*zpti_glf(j)/
     .                                           arho_iglf             !m/sec
                            vi_eff_m(j)  = vi_eff_m(j)*100.
                            chii_m(j) = 0.0
                          endif
                       else
                           vi_eff_m(j) = 0.0
                       endif
                       if(etaphi_m(j) .lt. 0.0)then
                          if (no_w_convection .eq. 1)then
                            etaphi_m(j) = 0.0
                            vrot_eff_m(j) =0.0
                          else
                            vrot_eff_m(j) = etaphi_m(j)*zpti_glf(j)/  
     .                             arho_iglf                   !m/sec
                            etaphi_m(j) = 0.0
                          endif
                       else
                            vrot_eff_m(j) = 0.0  !flux for vrot not implemented
                       endif

                    enddo  ! j loop
               else  !dv_method

#ifdef PARALLEL_12 
               call STOP('Parallel ONETWO does not support dv_method',1)
#endif
                  igrad = 1
                  iptrpa = 2
                  if(iptrpa .eq. 2)then
                     jstart = 0
                     jend = 1
                  else if(iptrpa .eq. 3)then
                     jstart = -1
                     jend = 1
                     call STOP('diffuse,iptrpa = 3 not implemented',1)
                  else
                     call STOP('diffuse,iptrpa not set ',1)
                  endif
                  vi_eff_m(1) = 0.0
                  ve_eff_m(1) = 0.0
                  vrot_eff_m(1) = 0.0
                  if(first_step .ne. 1  .and. ddebug(50) .gt. 0.99 )
     .                           jstart = jend !do one pass only
                  do j = 2,jmaxm-1
                     jmm = j
                     do iptrb =jstart,jend
                        if(itte_dv .eq. 1)then
                           zptej_glf = (1.+iptrb*dv_delt) *zpte_glf(j)
                        else
                            zptej_glf = zpte_glf(j)
                        endif
                        if(itti_dv .eq. 1)then
                            zptij_glf = (1.+iptrb*dv_delt) *zpti_glf(j)
                        else
                            zptij_glf = zpti_glf(j)
                        endif
                        zptij_p(iptrb) = zptij_glf
                        if(itenp_dv .eq. 1)then
                            zpnij_glf =  (1.+iptrb*dv_delt) *zpni_glf(j)
                        else
                            zpnij_glf = zpni_glf(j)
                        endif
                        if(itene_dv .eq. 1)then
                            zpnej_glf =  (1.+iptrb*dv_delt) *zpne_glf(j)
                        else
                            zpnej_glf = zpne_glf(j)
                        endif
                        if(itangrot_dv .eq. 1)then
                            zpnij_glf = (1.+iptrb*dv_delt) *zpni_glf(j)
                        else
                            zpnij_glf = zpni_glf(j)
                        endif

                          call callglf2d(           !inputs
     &                               jin,jout, !jmp.ibm.par
     &                               jeigen_iglf,jroot_iglf, glf23_iglf,         
     &                                  jshoot, jmm, jmaxm, itport_pt  , 
     &                                 irotstab, te_glf, ti_glf, ne_glf, 
     &                    ni_glf, ns_glf,igrad, idengrad_glf, zptej_glf,
     &                                  zptij_glf, zpnej_glf, zpnij_glf,
     &                             angrotp_exp, egamma_exp, gamma_p_exp,
     &                                  vphim_glf, vparm_glf, vperm_glf,
     &                        zeff_exp, btor_exp, ibtflag_glf, rho_glf,
     &                                  arho_iglf, grho1_glf, grho2_glf,
     &                         rmin_glf, rmaj_glf, rmajor_glf, zimp_glf, 
     &                         amassimp_glf, q_glf, shat_exp, alpha_exp,
     &                               elong_glf, amassgas_glf,alphae_glf,
     &                                             x_alpha_glf, i_delay,      !outputs:
     &                                  diffnemdv, chietemdv,
     &                         chie_etg_tem,chiitimdv, 
     &                          etaphimdv, etaparmdv, etapermdv,exchmdv, 
     &                                   diff_mdv, chie_mdv,chie_etg_m,
     &                    chii_mdv,etaphi_mdv, etapar_mdv, etaper_mdv,
     &                    exch_mdv, egamma_mdv, egamma_d, gamma_p_mdv,
     &               anrate_mdv, anrate2_mdv, anfreq_mdv, anfreq2_mdv,
     &               gamma_net_i,gamma_net_e)
c                           write(6,*) j,egamma_m(j+1), 'egamma_m-j+1'
                            chie_dv(iptrb) = chietemdv
                            chii_dv(iptrb) = chiitimdv                       
                            chitorot_dv(iptrb)  = etaphimdv
                            if(ddebug(50) .gt. 0.99)then !for itti_dv only
                               if(first_step .eq. 1 .and. 
     .                                                iptrb .eq. 0)then
                                  chii_save(j) = chii_dv(0)
                                  zpti_save(j) = zptij_p(0)
                               else if(first_step .eq. 0)then
                                   !only iprtrb =1 pass was made, define
                                   !index 0 arrays in this case:
                                  chii_dv(0) = chii_save(j)
                                  zptij_p(0)  = zpti_save(j)
                                  if(itte_dv .eq. 0)chie_dv(0) = 
     .                                                     chie_dv(1)
                                  if(itangrot_dv .eq. 0)
     .                                  chitorot_dv(0) = chitorot_dv(1)
                                endif
                            endif
                        enddo     !iptrb loop
                        !output arrays ending in dv are zeroed each time
                        !callglf2d is called. this necessitates the following:
                        if(irotstab .eq. 0)then
                           egamma_m(j) = egamma_exp(j)
                        else
                           egamma_m(j) = egamma_mdv(j+1)
                        endif
                        anfreq_m(j) = anfreq_mdv(j+1)
                        anrate_m(j) = anrate_mdv(j+1)
                        ve_eff_m(j) =0.0
                        if(itte_dv .eq. 1)then
                           !assumes iptrpa = 2
                           chie_m(j) = (chie_dv(1)*(1.+dv_delt)
     .                                   - chie_dv(0))/dv_delt    !m**2/sec
                            if(chie_m(j) .gt. 0.0)then
                              ve_eff_m(j)  = (chie_dv(0) -chie_m(j))*
     .                                         zpte_glf(j)/arho_iglf   
                              ve_eff_m(j)  = ve_eff_m(j)*100. !cm/sec
                           else   !
                              chie_m(j) = 0.0 
                              ve_eff_m(j) = chie_dv(0)*zpte_glf(j)/
     .                                           arho_iglf             !m/sec
                              ve_eff_m(j)  = ve_eff_m(j)*100.          !cm/sec
                           endif
                        else  !without dv method
                           if(chie_dv(0) .ge. 0.0)then
                              chie_m(j) = chie_dv(0)
                              ve_eff_m(j)  = 0.0
                           else
                              chie_m(j) = 0.0
                              ve_eff_m(j) = chie_dv(0)*zpte_glf(j)/
     .                                           arho_iglf             !m/sec
                              ve_eff_m(j)  = ve_eff_m(j)*100.
                           endif
                        endif
c                        print *,'j,chie = ',j,chie_m(j),ve_eff_m(j)

                        vi_eff_m(j) =0.0
                        if(itti_dv .eq. 1)then
                           !assumes iptrpa = 2
                           chii_m(j) = (chii_dv(1)*zptij_p(1)
     .                                      - chii_dv(0)*zptij_p(0))
     .                                      /(zptij_p(1)-zptij_p(0))
                           if(chii_m(j) .gt. 0.0)then
                              vi_eff_m(j)  = (chii_dv(0) -chii_m(j))*
     .                                         zpti_glf(j)/arho_iglf   !m/sec
                              vi_eff_m(j)  = vi_eff_m(j)*100.          !cm/sec
                           else   !
                              chii_m(j) = 0.0 
                              vi_eff_m(j) = chii_dv(0)*zpti_glf(j)/
     .                                           arho_iglf             !m/sec
                              vi_eff_m(j)  = vi_eff_m(j)*100.          !cm/sec
                           endif
                        else   !without dv method
                           if(chii_dv(0) .ge. 0.0)then
                               chii_m(j) = chii_dv(0)                  !m**2/sec
                           else
                               chii_m(j) = 0.0
                               vi_eff_m(j) = chii_dv(0)*zpti_glf(j)/
     .                                           arho_iglf 
                               vi_eff_m(j)  = vi_eff_m(j)*100.          !cm/sec
                           endif
                        endif

                        if(no_te_convection .eq. 1)ve_eff_m(j) = 0.0
                        if(no_ti_convection .eq. 1)vi_eff_m(j) = 0.0
                        vrot_eff_m(j) =0.0
                        if(itangrot_dv .eq. 1)then
                  call STOP('dv method not implemented for rotation',0)
                           !assumes iptrpa = 2
                           etaphi_m(j) = (chitorot_dv(1)*(1.+dv_delt) 
     .                                   - chitorot_dv(0))/dv_delt
                           if(etaphi_m(j) .gt. 0.0)then
                              vrot_eff_m(j)  = (chitorot_dv(0) 
     .                                          -etaphi_m(j))*
     .                                         zpti_glf(j)/arho_iglf   !m/sec
                              vrot_eff_m(j)  = vrot_eff_m(j)*100.      !cm/sec
                           else   !
                              etaphi_m(j) = 0.0 
                              vrot_eff_m(j) = chitorot_dv(0)*
     .                                        zpti_glf(j)/arho_iglf    !m/sec
                              vrot_eff_m(j)  = vrot_eff_m(j)*100.      !cm/sec
                           endif
                        else
                           etaphi_m(j) = chitorot_dv(0)                 !m**2/sec
                        endif
                     enddo  ! j loop


              endif    !dv_method

           else                           !multiple processor (mpi version)
                 !in order to repeatedly start up and shutdown MPI
                 !we spawn an external process here that calls callglf2d
                 !under control of MPI:
                 !on taurus /usr/local/mpich-1.2.4/bin/mpirun -np glf23_ncpus
!                all of the arguments in writeglf2d are inputs to 
!                callglf2d execept litport,njm1 

           endif
       endif



          if(glf_debug .eq. 1)then
c                debug glf23 section
                 write(ncrt,'("glf_debug selected, file "
     .                               "glfdebug created")')
                 iishot = 0
                 ioglf = 44

                 call getioun(ioglf,ioglf)
                 open(unit = ioglf, file = 'glfdebug', 
     .                status = 'UNKNOWN', iostat = iostat)
 

                 do j=1,jmaxm
                    if(j .eq. 1)then
                       write(ioglf,'(a12,i9)')'shot      = ',iishot
                       write(ioglf,'(a12,i3)')'leigen    = ',jeigen_iglf            !new,HSJ
                       write(ioglf,'(a12,i3)')'jmaxm     = ',jmaxm                 !new,HSJ
                       write(ioglf,'(a12,i3)')'i_delay   = ',i_delay             !new,HSJ
                       write(ioglf,'(a12,5(x,i3))')'itport_pt = ',
     &                                             (itport_pt(i),i=1,5)    !new,HSJ
                       write(ioglf,'(a12,i3)')'nroot      = ',jroot_iglf
                       write(ioglf,'(a12,i3)')'jshoot     = ',jshoot
                       write(ioglf,'(a12,i3)')'jmm        = ',jmm
                       write(ioglf,'(a12,i3)')'irotstab   = ',irotstab
                       write(ioglf,'(a12,i3)')'igrad      = ',igrad
                       write(ioglf,'(a12,i3)')'idengrad = ',idengrad_glf
                       write(ioglf,'(a12,i3)')'bt_flag   = ',ibtflag_glf
                       write(ioglf,'(a12,f10.6)')'bt_exp    = ',btor_exp
                       write(ioglf,'(a12,f10.6)')'arho_exp = ',arho_iglf
                      write(ioglf,'(a12,f10.6)')'amassgas=',amassgas_glf
                   write(ioglf,'(a12,f10.6)')'amassimp  = ',amassimp_glf
                       write(ioglf,'(a12,f10.6)')'zimp_exp  = ',zimp_glf
                     write(ioglf,'(a12,f10.6)')'alpha_e   = ',alphae_glf
                       write(ioglf,'(a12,f10.6)')'x_alpha= ',x_alpha_glf
                    endif
                    rova = r(j)/r(nj)          !this value is not passed to glf23 (rho_glf is) HSJ
                    write(ioglf,4321)j,rova,te_glf(j),ti_glf(j),
     .                     ne_glf(j),ni_glf(j),zeff_exp(j),ns_glf(j)        !modified HSJ
                    write(ioglf,4321)j,rova,rmin_glf(j),rmaj_glf(j),
     .                           grho1_glf(j),grho2_glf(j),
     .                           elong_glf(j)
                    write(ioglf,4321)j,rova,elong_glf(j),q_glf(j),
     .                     shat_exp(j),
     .                     alpha_exp(j),angrotp_exp(j),egamma_exp(j)
                    write(ioglf,4321)j,rova,zpne_glf(j),zpni_glf(j),
     .              zpte_glf(j),zpti_glf(j),vphim_glf(j),egamma_exp(j)

c                   if jmm  .ne. 0 then the following contains output for the
c                   single grid point given in jmm
                    write(ioglf,4321)j,rova,diffnem,chietem,chiitim,
     .                                                        exchm
                    write(ioglf,4321)j,rova,etaphim,etaparm,etaperm    !new HSJ

                    write(ioglf,4321)j,rova,diff_m(j),chie_m(j),
     .                                                   chii_m(j)
                    write(ioglf,4321)j,rova,exch_m(j),egamma_m(j),
     .                                                egamma_d(j,1)
                    write(ioglf,4321)j,rova,gamma_p_m(j),anrate_m(j),
     .                                                    anrate2_m(j)
                    write(ioglf,4321)j,rova, anfreq_m(j), anfreq2_m(j)
                    write(ioglf,4321)j,rova, gamma_p_exp(j), 
     .                            vparm_glf(j),vperm_glf(j),rho_glf(j)  !new HSJ
                    write(ioglf,4321)j,rova, etaphi_m(j),etapar_m(j),
     .                                        etaper_m(j)               !new HSJ
                 enddo
 4321            format(2x,i3,2x,e8.3,6(2x,e12.5))
c                call exit                      ! stop after single call
                 close(unit = ioglf)
          endif                            ! if(glf_debug .eq. 1)





 4322     continue
          navg_chi = navg_chi+1
          if(freeze_xte .eq. 0)then
            chie_ms(1:nj) = chie_m(:)  !save for possible reuse next iteratin
            chie_avg(1:nj) = (chie_avg(1:nj)*(navg_chi-1)
     .                                        +chie_m(1:nj))/navg_chi 
          else
             if(use_avg_chi .eq. 1)then
                     chie_m(1:nj) = chie_avg(1:nj)
             else
                     chie_m(:) = chie_ms(1:nj)  !set to saved value if 
                                                ! chie_m if tte is frozen
             endif
          endif
          if(freeze_xti .eq. 0)then
            chii_ms(1:nj) = chii_m(:)  !save copies for possible reuse next iteratin
            chii_avg(1:nj) = (chii_avg(1:nj)*(navg_chi-1)
     .                                        +chii_m(1:nj))/navg_chi 
          else 
             if(use_avg_chi .eq. 1)then
                    chii_m(1:nj) = chii_avg(1:nj)
             else
                    chii_m(1:nj) = chii_ms(1:nj)  !set to saved value if 
                                            ! chii_m if tti is frozen
             endif
          endif

          if(freeze_xwr .eq. 0)then
            etaphi_ms(1:nj) = etaphi_m(:)
          else
            etaphi_m(:) = etaphi_ms(1:nj)
          endif

          !-jmp.den start-------------
         
          if(freeze_xni .eq. 0)then  
            diff_ms(1:nj) = diff_m(:)
          else
            diff_m(:) = diff_ms(1:nj)
          endif

          if(time .lt. time0+iglf_idt) diff_m(:) = 0.0
          if(time .lt. time0+iglf_idt) chii_m(:) = 0.0
          if(time .lt. time0+iglf_idt) chie_m(:) = 0.0
          if(time .lt. time0+iglf_idt) etaphi_m(:) = 0.0

          diff_m(1:nj) = diff_m(1:nj) + chid_glf_d         !*****(1)
          diff_m(1:nj) = diff_m(1:nj) + chid_glf_e*chie_m(1:nj) !*****(2)

          if(ichid_smooth .eq. 1) then !*****(3)
             do j=2,nj-2
             	 chid_avg(j) = diff_m(j-1)+diff_m(j)+diff_m(j+1)
             	 chid_avg(j) = chid_avg(j)/3.0
                !if(chid_avg(j) .le. 0.0) chid_avg(j) = 0.0
             end do
             diff_m(2:nj-2) = chid_avg(2:nj-2)
          endif   

          !-jmp.den start-------------

c
c --- cparam is continuation parameter, normally =1.0 HSJ 3/27/01 ----------` 
c --- Reuse IFS variables names:
          cparam1 =cparam 
          do j=1,nj-1
               if(chii_m(j) .gt. 0.0)then
                    chi_i_ifs(j) = cparam1*chii_m(j) * 1.e4
                    chi_i_ifs(j) = MIN(chi_i_ifs(j),1.D6) !jmp.ibm
               else
                    chi_i_ifs(j) = 0.0
               endif


               if(chie_m(j) .gt. 0.0)then
                    chi_e_ifs(j) =  cparam1*chie_m(j) * 1.e4
                    chi_e_ifs(j) = MIN(chi_e_ifs(j),1.D6) !jmp.ibm
               else
                    chi_e_ifs(j) = 0.0
               endif




c
c              get results into d and xkang for simulation purposes:
c

cjmp.den    d_ifs(j)     = chi_i_ifs(j)
            d_ifs(j)     = diff_m(j) * grho2_glf(j) * 1.e4 !jmp.den



            xkang_ifs(j) = etaphi_m(j) * 1.e4

c               modify the chie, chii, and di, for the
c               flow shear suppression effect if ifsflag > 0
c
c                chi_e_ifs(j) = chi_e_ifs(j) * coefa(j)
c                xkang_ifs(j) = xkang_ifs(j) * coefd(j)
c                chi_i_ifs(j) = chi_i_ifs(j) * coefb(j)
c                d_ifs(j)     = d_ifs(j)     * coefc(j)
c
c???????????????????????????????????????????????????????????????????????
c
c --- relax just the anomalous part of the diffusivity (i.e., without the
c --- neoclassical term) note that at this point chi_e_ifs and chi_i_ifs
c --- are just the anomalous part returned by subroutine callglf2d.
c --- we use the Weiland storage vectors for the GLF23 model here (since
c --- there is no conflict,only one of the models is active in a given run)
c --- If itercorct > 1, then xchie_weilandsv and xchii_weilandsv are
c --- just the anomalous parts from the previous iteration.
c
              if (itercorct .gt. 1) then
                d_ifs(j)     =(1.0 - relaxrebut) *xchii_weilandsv(j)+
     .                                relaxrebut *d_ifs(j)
                chi_e_ifs(j) =(1.0 - relaxrebut) *xchie_weilandsv(j)+
     .                                relaxrebut *chi_e_ifs(j)
                chi_i_ifs(j) =(1.0 - relaxrebut) *xchii_weilandsv(j)+  
     .                                relaxrebut *chi_i_ifs(j)
                xkang_ifs(j) =(1.0 - relaxrebut) *xkangwlsv(j) +
     .                                relaxrebut *xkang_ifs(j)
              end if
c
                xchie_weilandsv(j)   = chi_e_ifs(j)
                xchii_weilandsv(j)   = chi_i_ifs(j)
                xkangwlsv(j)         = xkang_ifs(j)
c
                enea=0.5*(ene(j+1)+ene(j))
                ensum = 0.0
c                dmassden = 0.0
                do i=1,nprim    ! have to use nprim to be consistent
****            do i=1,nion     ! with calculation of xkiinv (in PSOURC)
                  enav = 0.5 * (en(j,i)+en(j+1,i))
                  ensum = enav +ensum
                  d(i,i,j) = d(i,i,j) + d_ifs(j)
c                  dmassden = dmassden+atw(i)*xmassp*enav
                end do
c
               xke_ifs(j)           = chi_e_ifs(j)*0.5*(ene(j+1)+
     .                                                  ene(j))
c             add < (grad rho)**2 >  non circular correction factor
c             Note that it is added to the d's but is not otherwise
c             added to the chis that are printed out,in orer to facilitate
c             comparison with other codes. But this factor is included in
c             the xketot,xkitot,xchietot,xchiitot,xkangtot values below HSJ
c             also extrapolate chiineo and  chiwneo to near zero
              if(j .le. j_clamp)then
                 if(d(nion+2,nion+2,j) .gt. 0.0)then
                  clamp_multti = (slope_clamp*r(j)+chi_clamp)*ensum
     .                                             /d(nion+2,nion+2,j)
                 else
                   clamp_multti = 1.
                 endif
                 if(d(nion+4,nion+4,j) .gt. 0.0)then
                  clamp_multw  = (slope_clamp*r(j)+
     .                                chi_clamp)*dmassden(j)*rmajor**2
     .                                             /d(nion+4,nion+4,j)
                 else
                    clamp_multw = 1.0
                 endif
              else
                  clamp_multti = 1.0
                  clamp_multw  = 1.0
              endif
               d(nion+1,nion+1,j)   = d(nion+1,nion+1,j)           ! te
     .                               + iglf_chie*grho2_glf(j)*xke_ifs(j)

               xki_ifs(j)           = chi_i_ifs(j)*0.5*(
     .           tot_primary_ion_den(j) + tot_primary_ion_den(j+1))
               d(nion+2,nion+2,j)   = d(nion+2,nion+2,j)*clamp_multti ! ti
     .                               + iglf_chii*grho2_glf(j)*xki_ifs(j)
               if (iangrot .ne. 0)
     .         d(nion+4,nion+4,j)   = d(nion+4,nion+4,j)*clamp_multw ! tord. rot.
     .        +iglf_chiv*grho2_glf(j)*xkang_ifs(j)*dmassden(j)*rmajor**2
c

               xdchitot(j)          = MAX(d(1,1,j),0.0_DP)
               xketot(j)            = d(nion+1,nion+1,j)
               xchietot(j)          = xketot(j)/enea
               xkitot(j)            = d(nion+2,nion+2,j)
               xchiitot(j)          = xkitot(j)/ensum
               xkangtot(j)          = d(nion+4,nion+4,j)

          enddo

          do j =1, jelc_clamp
             d(nion+1,nion+1,j) =  d(nion+1,nion+1,jelc_clamp)
             xketot(j)          =  d(nion+1,nion+1,j)
             xchietot(j)        =  xchietot(jelc_clamp)
          enddo
          do j =1, jion_clamp
             d(nion+2,nion+2,j) =  d(nion+2,nion+2,jion_clamp)
             xkitot(j)          =  d(nion+2,nion+2,j)
             xchiitot(j)        =  xchiitot(jion_clamp)
          enddo
c
c 
c ----------------------------------------------------------------------
c numerical stabilization 
c ----------------------------------------------------------------------

               IF(use_stab_flux)THEN
                  IF(itncount ==1)THEN               ! set new level of dbar  on first iteration 
                                                     ! of each time step
                                                     ! note: intcount is defined in solve_newton
                     xketotm   = MAXVAL(xketot)
                     xkitotm   = MAXVAL(xkitot)   
                     xkangtotm = MAXVAL(xkangtot)
                     call source_mod_stab(1,j,xketotm,
     .                                    xkitotm,xkangtotm)
                  ENDIF
                  DO  j=1,nj-1    ! adjust d on every iteration
                     i = nion +1 
                     d(i,i,njin)  = dbar(i,njin) + d(i,i,njin)
                     i = nion + 2
                     d(i,i,njin)  = dbar(i,njin) + d(i,i,njin)
                     i= nion+4
                     d(i,i,njin)  = dbar(i,njin) + d(i,i,njin)
                  ENDDO

               ENDIF



c
c..diagnostic printout
c 
         if ( lprint_glf .gt. 0 .and. time .ge. timmax-dt ) then
          open(unit=12,file='glf.out',status='unknown')
          write(12,*) ' '
          write(12,*) ' GLF23 diagnostic printout'
          write(12,*) ' time = ',time,' time-max = ',timmax
          write(12,*) ' R = ',rmajor_glf
          write(12,*) ' a = ',rminor_glf
          write(12,*) ' B = ',btor_exp
          write(12,*) ' arho = ',arho_iglf
          write(12,*) ' amass = ',amassgas_glf
          write(12,*) ' jmaxm = ',jmaxm
c
          write(12,190)
          do j=1,nj
            write(12,199) j, rho_glf(j), te_glf(j), ti_glf(j),
     .        ne_glf(j), ni_glf(j), enbeam(j), zeff_exp(j)
          enddo
          write(12,191)
          do j=1,nj
            write(12,199) j, rho_glf(j), rmin_glf(j), rmaj_glf(j),
     .        grho1_glf(j), grho2_glf(j), elong_glf(j)
          enddo
          write(12,192)
          do j=1,nj
            write(12,231) j,rho_glf(j),csda_glf(j),rhosda_glf(j),
     .        bteff_glf(j), drhodr(j), drhodrrrho(j)
          enddo
          write(12,193)
          do j=1,nj
            write(12,231) j, rho_glf(j), zpte_glf(j), zpti_glf(j),
     .        zpne_glf(j), zpni_glf(j)
          enddo
          write(12,194)
          do j=1,nj
            write(12,231) j, rho_glf(j), q_glf(j), shat_exp(j),
     .        beta_glf(j), alpha_exp(j)
          enddo
          write(12,195)
          do j=1,nj
            write(12,231) j, rho_glf(j), vparm_glf(j),
     .         vperm_glf(j), vphim_glf(j), ve_glf(j)
          enddo
          write(12,196)
          do j=1,nj
            write(12,231) j, rho_glf(j), angrot(j), angrotp_exp(j),
     .        egamma_exp(j), gamma_p_exp(j)
          enddo
          write(12,208)
          do j=1,nj
            write(12,230) j, rho_glf(j), chie_m(j), chii_m(j),
     .        etaphi_m(j), etapar_m(j), etaper_m(j)
          enddo
c
c
 190      format(/,1x,'j',5x,'rho',7x,'te',10x,'ti',10x,'ne',
     .           10x,'ni',10x,'ns'10x,'zeff')
 191      format(/,1x,'j',5x,'rho',8x,'rminor',6x,'rmajor',6x,
     .           'grho1',7x,'grho2',7x,'kappa')
 192      format(/,1x,'j',5x,'rho',9x,'csda',7x,'rhosda',6x,
     .           'bteff',6x,'drhodr',4x,'drhodrrrho')
 193      format(/,1x,'j',5x,'rho',10x,'zpte',9x,'zpti',9x,'zpne',
     .           9x,'zpni')
 194      format(/,1x,'j',5x,'rho',11x,'q',11x,'shat',9x,'beta',
     .           9x,'alpha')
 195      format(/,1x,'j',5x,'rho',10x,'vpar',9x,'vper',9x,
     .           'vphi',10x,'ve')
 196      format(/,1x,'j',5x,'rho',10x,'angrot',7x,'anrotp',
     .           7x,'egamma',7x,'gamma_p')
 199      format(i2,2x,f8.4,2x,6(1pe10.4,2x))
 200      format(/,1x,'j',5x,'rho',9x,'te',10x,'ti',10x,'ne',
     .           10x,'ni',10x,'zeff')
 201      format(/,1x,'j',5x,'rho',8x,'rminor',6x,'rmajor',7x,
     .           'grho1',7x,'grho2')
 202      format(/,1x,'j',5x,'rho',9x,'zpte',8x,'zpti',8x,'zpne',
     .           8x,'zpni')
 203      format(/,1x,'j',5x,'rho',11x,'q',11x,'shat',9x,'beta',
     .           9x,'alpha')
 208      format(/,1x,'j',5x,'rho',9x,'chie',8x,'chii',
     .           7x,'etaphi',6x,'etapar',6x,'etaper')
 209      format(/,1x,'j',5x,'rho',9x,'csda',7x,'rhosda',6x,
     .           'angrot',6x,'anrotp')
 210      format(/,1x,'j',5x,'rho',10x,'vpar',9x,'vper',9x,
     .           'vphi',10x,'ve')
 211      format(/,1x,'j',5x,'rho',11x,'fc',9x,'akappa1',6x,
     .           'akappa2',6x,'alpha_neo')
 230      format(i2,2x,6(1pe10.4,2x))
 231      format(i2,2x,6(1pe11.4,2x))
         endif  !end print block
c

          irotstab = irotstab_save
        endif  glf_sect  !end glf23 block

c---------------------------------------------------------------------------
c     BRAINFUSE OM 24/01/14
c---------------------------------------------------------------------------
        brainfuse_sect:    IF(include_brainfuse .EQ. 1 )THEN
!     
!     Deuterium/Tritium densities

           do jj=1,nprim        ! loop to find index of deuterium
              if(atw(jj) .eq. 2.)  exit
           enddo

           CALL BRAINFUSE(nj,   ! number of points
     .          r(nj)*0.01,     ! a : minor radius [m]
     .          btor*1.0e-4,    ! Bt: Toroidal magnetic field on axis [T]
     .          r*0.01,         ! r[j] : radius [m]
     .          ravg_r*0.01,    ! R[j] : major radius [m] = rmajor+shif(j)
     .          kappa_nj,       ! kappa[j] : elongation
     .          ene*1E6,        ! ne[j] : electron density [m-3]
     .          en(:,jj)*1E6,   ! ni[j] : ion density [m-3]
     .          te,             ! Te[j] : electron temperature [keV]
     .          ti,             ! Ti[j]: ion temperature [keV]
     .          q               ! q[j] : safety factor with sign
     .          )
        ENDIF  brainfuse_sect

c---------------------------------------------------------------------------
c     MMM95 JMP 3/23/06
c---------------------------------------------------------------------------

      mmm_sect:    if(include_mmm .eq. 1 )then

ccc---ExB shear calculation ( see glf_sect )--------------------------------

          rmajor_glf=rmajor/100. ! major radius, R (m)
          btor_exp= btor/1.e4    ! toroidal field, Bt (T)
          arho_iglf=r(nj)/100.   ! toroidal flux at LCFS, rho (m)

c..compute geometric quantities, densities, temperatures, q
c  note: ni is primary ion density

          do j=1,nj
            rho_glf(j)=max( r(j)/(100.*arho_iglf), 1.D-6 ) !jmp.ibm
            tau_glf(j)=ti(j)/te(j)
            ni_glf(j)=0.0
            do jj=1,nprim
              if (atw(jj) .le. 3.1) then
                 ni_glf(j)=ni_glf(j)+en(j,jj)
              end if
            end do
            ni_glf(j)=1.e-13*ni_glf(j)
            q_glf(j)=ABS(q(j))
            rmin_glf(j) = max( rminor_r(j)/100., 1.D-6) !jmp.ibm
            rmaj_glf(j) = ravg_r(j)/100.
          enddo

          do j=2,nj
            drhodr(j)=(rho_glf(j)-rho_glf(j-1))*arho_iglf/
     .         (rmin_glf(j)-rmin_glf(j-1)+1.e-6)
            drhodrrrho(j)=drhodr(j)*rmin_glf(j)/arho_iglf/rho_glf(j)
          enddo
          drhodr(1)=drhodr(2)
          drhodrrrho(1)=drhodrrrho(2)

c..compute csda and rhosda using effective B-field
c  if ibtflag.gt.0 (default)

          do j=1,nj

            amassgas_glf = 0.0
            tot_den = 0.0
            do jj=1,nprim
               if (atw(jj) .le. 3.1) then
                  amassgas_glf = amassgas_glf + atw(jj)*en(j,jj)
                  tot_den=tot_den+en(j,jj)
               end if
            end do
            amassgas_glf = amassgas_glf/tot_den

            csda_glf(j)=9.79e5*(te(j)*1.e3)**.5/
     .         (arho_iglf*100.)/amassgas_glf**.5

c           if(ibtflag_glf.gt.0) then
              bteff=btor_exp*rho_glf(j)*arho_iglf/
     .                   rmin_glf(j)*drhodr(j)
c           else 
c             bteff=botr_exp
c           endif
            rhosda_glf(j)=((1.02e2*(te(j)*1.e3)**.5)/
     .         bteff/1.e4)*(amassgas_glf**.5)/(arho_iglf*100.)
          enddo
c
c..compute normalized gradients
c
          do j=2,nj
            drho_glf=rho_glf(j)-rho_glf(j-1)  !rho_glf is [0,1]
            tia = 0.5*(ti(j)+ti(j-1))
            ani = 0.5*(ni_glf(j)+ni_glf(j-1)) 
            zpti_glf(j-1)=(-ti(j)+ti(j-1))/(tia*drho_glf) ! a/LTi
            zpni_glf(j-1)=(-ni_glf(j)+ni_glf(j-1))/(ani*drho_glf) ! a/Lni
          enddo

c..compute plasma toroidal angular velocity (angrotp) in 1/s
c  and Vexb in m/s

          do j=1,nj

            fc=1-1.46*(rmin_glf(j)/rmaj_glf(j))**0.5+
     .         0.46*(rmin_glf(j)/rmaj_glf(j))**1.5
            akappa1=0.8839*fc/(0.3477+0.4058*fc)
            alpha_neo=-akappa1+1.
c
            vstar_sign=-1.
            corot=1.
            pgeo_local=drhodr(j)
            rdrho_local=rmin_glf(j)/arho_iglf/rho_glf(j)
            angrotp_exp(j)=corot*angrot(j)

            !set ve to value used in callglf2d  HSJ 05/0103:
            ve_glf(j)=-tau_glf(j)*csda_glf(j)*arho_iglf*
     .         rhosda_glf(j)*(zpni_glf(j) +
     .         alpha_neo*zpti_glf(j))*vstar_sign*pgeo_local

            !add sigma*(rminor/(Rmajor*q))Vtoroidal to Vexb.
            !Note that sigma is assumed  +1 here because the
            !toroidal rotation velocity is a signed quantitiy.
            !This term is identical to eq. 16a in Waltz,
            !phys plas, 4,7,97,2486. HSJ:
            ve_glf(j)=ve_glf(j)-rdrho_local*rho_glf(j)*
     .         arho_iglf/rmajor_glf/q_glf(j)*
     .         rmajor_glf*angrotp_exp(j)
          enddo
c
c..compute ExB shear rates
c
          do j=1,nj-1        
            rdrho_local=rmin_glf(j)/arho_iglf/rho_glf(j)
            rdrho_local_p1=rmin_glf(j+1)/arho_iglf/rho_glf(j+1)   !HSJ 4/30/03
            egamma_exp(j)=drhodrrrho(j)*(rho_glf(j)+rho_glf(j+1))/
     .         (q_glf(j)+q_glf(j+1))*(ve_glf(j+1)*q_glf(j+1)/
     .         rho_glf(j+1)/rdrho_local_p1-ve_glf(j)*q_glf(j)/
     .         rho_glf(j)/rdrho_local)/(rho_glf(j+1)-rho_glf(j))/
     .         arho_iglf
c    .         arho_iglf/csda_glf(j)
          enddo

ccc----------------------------------------------------------------------------

c     print *,'mmm95 section JMP'
      nprout_mm = 6 
      lprint_mm = 0
      lsuper_mm = 0
      lreset_mm = 0

      npoints_mm = 1

      do j=1,nj-1

      grho2_mm = grho2_ifs(j)

      rminor_mm = rminor*0.01
      rmajor_mm = ravg_r(j)*0.01 !rmajor+shif(j)
      elong_mm = kappa_nj(j)
      dense_mm = 1.e6*ene(j)
      xzeff_mm  = zeff(j)
      tekev_mm = te(j)
      tikev_mm = ti(j)
      q_mm = ABS(q(j))
      btor_mm = btor*1.0e-4

      dens_min_mm = 1.e-5

c-----hydrogen density

      densh_mm = 0.0
      densh_m2 = 0.0
      do jj=1,nprim
         if(atw(jj) .le. 3.1) then
            densh_mm = densh_mm+en(j  ,jj)
            densh_m2 = densh_m2+en(j+1,jj)
         end if
      end do

      if ( densh_mm .lt. dens_min_mm ) densh_mm = dens_min_mm
      if ( densh_m2 .lt. dens_min_mm ) densh_m2 = dens_min_mm

c-----impurity density

      densimp_mm = 0.0
      densimp_m2 = 0.0
      zdensimp_mm = 0.0
      zdensimp_m2 = 0.0
      do jj=nprim+1,nion
         densimp_mm = densimp_mm+en(j  ,jj)
         densimp_m2 = densimp_m2+en(j+1,jj)
         zdensimp_mm = zdensimp_mm+z(j,jj  )*en(j  ,jj)
         zdensimp_m2 = zdensimp_m2+z(j+1,jj)*en(j+1,jj)
      end do

      if ( densimp_mm .lt. dens_min_mm ) densimp_mm = dens_min_mm
      if ( densimp_m2 .lt. dens_min_mm ) densimp_m2 = dens_min_mm

c-----fast ion density

      densfe_mm= enbeam(j)

c-----total density

      denstot_mm = 0.0
      denstot_m2 = 0.0
      do jj=1,nion
         denstot_mm = denstot_mm + en(j  ,jj)
         denstot_m2 = denstot_m2 + en(j+1,jj)
      end do

c-----average density weighted charge & atomic mass of impurities
      
      avezimp_mm = 0.0
      amassimp_mm = 0.0 
      do jj=nprim+1,nion
         avezimp_mm = avezimp_mm + z(j,jj)*en(j,jj)
         amassimp_mm = amassimp_mm + atw(jj)*en(j,jj)
      end do
      avezimp_mm = avezimp_mm/densimp_mm
      amassimp_mm = amassimp_mm/densimp_mm
      
c-----average density weighted charge & atomic mass of hydragen

      do jj=1,nprim
         if(atw(jj) .le. 3.1) then
            amasshyd_mm = amasshyd_mm+atw(jj)*en(j,jj)
         end if
      end do
      amasshyd_mm = amasshyd_mm/densh_mm

c-----density weighted atomic mass of total ions

      aimass_mm = 0.0
      do jj=1,nion
         aimass_mm =  aimass_mm + atw(jj)*en(j,jj)
      end do
      aimass_mm =  aimass_mm/denstot_mm
      
c-----normalized gradients

c     rda_mm = -(2.0*rmajor+shif(j+1)+shif(j))/(r(j+1)-r(j))
      rda_mm = -(ravg_r(j+1)+ravg_r(j))/(r(j+1)-r(j))

      grdne_mm = rda_mm*(ene(j+1)-ene(j))/(ene(j+1)+ene(j))
      grdni_mm = rda_mm*(denstot_m2-denstot_mm)/(denstot_m2+denstot_mm)
      grdnh_mm = rda_mm*(densh_m2-densh_mm)/(densh_m2+densh_mm)
      grdnz_mm = rda_mm*(zdensimp_m2-zdensimp_mm)
     .                 /(zdensimp_m2+zdensimp_mm)
      grdte_mm = rda_mm*(te(j+1)-te(j))/(te(j+1)+te(j))
      grdti_mm = rda_mm*(ti(j+1)-ti(j))/(ti(j+1)+ti(j))
      grdq_mm  = rda_mm*(ABS(q(j+1))-ABS(q(j)))
     .                 /(ABS(q(j+1))+ABS(q(j)))

      densh_mm = 1.e6*densh_mm
      densh_m2 = 1.e6*densh_m2
      densimp_mm = 1.e6*densimp_mm
      densfe_mm = 1.e6*densfe_mm 

c-----ExB shear

      wexbs_mm = exbmult_mmm*abs(egamma_exp(j))

c----call mmm95

c     print *,'before mmm95JMP'
      call mmm95 (
     &   rminor_mm, rmajor_mm, elong_mm
     & , dense_mm, densh_mm, densimp_mm, densfe_mm
     & , xzeff_mm, tekev_mm, tikev_mm, q_mm, btor_mm
     & , avezimp_mm, amassimp_mm, amasshyd_mm, aimass_mm, wexbs_mm
     & , grdne_mm, grdni_mm, grdnh_mm, grdnz_mm
     & , grdte_mm, grdti_mm, grdq_mm
     & , thiig_mm, thdig_mm, theig_mm, thzig_mm
     & , thirb_mm, thdrb_mm, therb_mm, thzrb_mm
     & , thikb_mm, thdkb_mm, thekb_mm, thzkb_mm
     & , gamma_mm, omega_mm, difthi_mm, velthi_mm, vflux_mm
     & , matdim_mm, npoints_mm, nprout_mm, lprint_mm,  nerr_mm
     & , lsuper_mm, lreset_mm, lswitch_mm, cswitch_mm
     & , fig_mm, frb_mm, fkb_mm)
c     print *,thiig_mm,thirb_mm,thikb_mm,'end mmm95JMP'

      chi_i_ifs(j) = (thiig_mm + thirb_mm + thikb_mm)*1.e4
      chi_e_ifs(j) = (theig_mm + therb_mm + thekb_mm)*1.e4
      d_ifs(j)     = chi_i_ifs(j)
      xkang_ifs(j) = chi_i_ifs(j)

c     print *,j,thiig_mm,thirb_mm,thikb_mm,'JMP'

c      if (j.eq.48) then
c      print *,
c     &   rminor_mm, rmajor_mm, elong_mm
c     & , dense_mm, densh_mm, densimp_mm, densfe_mm
c     & , xzeff_mm, tekev_mm, tikev_mm, q_mm, btor_mm
c     & , avezimp_mm, amassimp_mm, amasshyd_mm, aimass_mm, wexbs_mm
c      pause
c      end if

c-----freeze jmp.den

      if(freeze_xte .eq. 0) then
         chi_e_ifs_save(j) = chi_e_ifs(j)
      else
         chi_e_ifs(j) = chi_e_ifs_save(j)
      endif

      if(freeze_xti .eq. 0) then
      	 chi_i_ifs_save(j) = chi_i_ifs(j)
      	 d_ifs_save(j) = d_ifs(j)
      else
         chi_i_ifs(j) = chi_i_ifs_save(j)
      	 d_ifs(j) = d_ifs_save(j)
      end if 	   
         
      if(freeze_xwr .eq. 0)then
         xkang_ifs_save(j) = xkang_ifs(j)
      else
         xkang_ifs(j) = xkang_ifs_save(j)
      endif

c-----common post process

              if (itercorct .gt. 1) then
                d_ifs(j)     =(1.0 - relaxrebut) *xchii_weilandsv(j)+
     .                                relaxrebut *d_ifs(j)
                chi_e_ifs(j) =(1.0 - relaxrebut) *xchie_weilandsv(j)+
     .                                relaxrebut *chi_e_ifs(j)
                chi_i_ifs(j) =(1.0 - relaxrebut) *xchii_weilandsv(j)+  
     .                                relaxrebut *chi_i_ifs(j)
                xkang_ifs(j) =(1.0 - relaxrebut) *xkangwlsv(j) +
     .                                relaxrebut *xkang_ifs(j)
              end if
c
                xchie_weilandsv(j)   = chi_e_ifs(j)
                xchii_weilandsv(j)   = chi_i_ifs(j)
                xkangwlsv(j)         = xkang_ifs(j)
c
                enea=0.5*(ene(j+1)+ene(j))
                ensum = 0.0
c                dmassden = 0.0
                do i=1,nprim    ! have to use nprim to be consistent
****            do i=1,nion     ! with calculation of xkiinv (in PSOURC)
                  enav = 0.5 * (en(j,i)+en(j+1,i))
                  ensum = enav +ensum
                  d(i,i,j) = d(i,i,j) + d_ifs(j)
c                  dmassden = dmassden+atw(i)*xmassp*enav
                end do
c
               xke_ifs(j)           = chi_e_ifs(j)*0.5*(ene(j+1)+
     .                                                  ene(j))
c             add < (grad rho)**2 >  non circular correction factor
c             Note that it is added to the d's but is not otherwise
c             added to the chis that are printed out,in orer to facilitate
c             comparison with other codes. But this factor is included in
c             the xketot,xkitot,xchietot,xchiitot,xkangtot values below HSJ
c             also extrapolate chiineo and  chiwneo to near zero
              if(j .le. j_clamp)then
                 if(d(nion+2,nion+2,j) .gt. 0.0)then
                  clamp_multti = (slope_clamp*r(j)+chi_clamp)*ensum
     .                                             /d(nion+2,nion+2,j)
                 else
                   clamp_multti = 1.
                 endif
                 if(d(nion+4,nion+4,j) .gt. 0.0)then
                  clamp_multw  = (slope_clamp*r(j)+
     .                                chi_clamp)*dmassden(j)*rmajor**2
     .                                             /d(nion+4,nion+4,j)
                 else
                    clamp_multw = 1.0
                 endif
              else
                  clamp_multti = 1.0
                  clamp_multw  = 1.0
              endif

               d(nion+1,nion+1,j)   = d(nion+1,nion+1,j)           ! te
     .                               + immm_chie*grho2_mm*xke_ifs(j)

               xki_ifs(j)           = chi_i_ifs(j)*0.5*(
     .            densh_mm + densh_m2)*1.0e-6
               d(nion+2,nion+2,j)   = d(nion+2,nion+2,j)*clamp_multti ! ti
     .                               + immm_chii*grho2_mm*xki_ifs(j)
               if (iangrot .ne. 0)
     .         d(nion+4,nion+4,j)   = d(nion+4,nion+4,j)*clamp_multw ! tord. rot.
     .        +immm_chiv*grho2_mm*xkang_ifs(j)*dmassden(j)*rmajor**2
c

               xdchitot(j)          = d(1,1,j)
               xketot(j)            = d(nion+1,nion+1,j)
               xchietot(j)          = xketot(j)/enea
               xkitot(j)            = d(nion+2,nion+2,j)
               xchiitot(j)          = xkitot(j)/ensum
               xkangtot(j)          = d(nion+4,nion+4,j)

c      write(*,2002),j,xchietot(j),xchiitot(j) 
c2002  format (2x,i3,2(2x,1pe15.7))
c-----common post process

      end do

      end if mmm_sect


c
c ----------------------------------------------------------------------
c     modified gyro Bohm and other ITB models HSJ 12/10/98
c ----------------------------------------------------------------------
c
      if (include_itb .eq. 1) then
         write(ncrt,'(" model itb model is not completely coded",/,
     .                2x, "onetwo must stop")')
        call STOP ('subroutine DIFFUS: ITB MODEL', 0)
        do j=1,nj-1
c
c         modified gyro Bohm:
c
c         calculate the gyro Bohm diffusion coefficient
c         (== dbohm*rhostar)
c
          dgbohm = dbohm*rhos
          qa     = 0.5*(q(j+1)+q(j))
          dqdrho = (q(j+1)-q(j))/dr(j) ! shear defined above = abs value
          shear  = (ra(j)/qa)*dqdrho   ! here we need it with sign
          f_mgb  = 1.0 / (1.0 + 2.25*(shear-0.6666)**2)
          tea    = 0.5*(te(j+1)+te(j))
          tia    = 0.5*(ti(j+1)-ti(j))
          chie_mgb(j) = ce0_mgb*f_mgb*((tea/tia)**alpha_mgb)
     .                 *tea*sqrt(tea)*qa*qa/(btor*btor) ! check units
****      chii_mgb(j) = ci1_mgb*chie_neo(j)
          if (dqdrho .gt. 0)
     .        chii_mgb(j) = chii_mgb(j) + ci2_mgb * chie_mgb(j)
c
c         Bohm gyro Bohm ITB model:
c
****      chie_bohm(j) = ???? not yet done
****      chie_gb(j)   = ????
          gammacr = c1_g + c2_g * omse   ! C1_G and OMSE not yet defined
          hfact   = ce_bgb * omecb / gammacr     ! OMECB not yet defined
          if (hfact .lt. 1.0) then
****        chii_bohm(j) = ???
            chii_bitb(j) = chii_bohm(j)*(1.0-hfact)
          else
            chii_bitb(j) = 0
          end if
          chie_bgb(j)= chie_bohm(j)+chie_gb(j)
****      chii_bgbitb(j)=chii_bitb(j)+chii_gb(j)+chi_neo(j)
          chief_mgb(j)=cfe_mgb*chie_mgb(j)+cfe_bgb*chie_bgb(j)
          xke_mgb(j)= chief_mgb(j)*0.5*(ene(j+1)+ ene(j))
          d(nion+1,nion+1,j)   = d(nion+1,nion+1,j) + xke_mgb(j)
****      chiif_mgb(j)=cfi_mgb*chii_mgb(j)+cfi_bgb*chii_bgitb(j)
          xki_mgb(j)           = chiif_mgb(j)*0.5*
     .           (tot_primary_ion_den(j) + tot_primary_ion_den(j+1))
          d(nion+2,nion+2,j)   = d(nion+2,nion+2,j)+xki_mgb(j)
        end do ! end loop over half grid
      end if   ! end include_itb model

c---------------------to match some ptor results------------------------
      if(neomult .ne. 0)then
         do j =1,nj-1
            qa = .5*(q(j+1)+q(j))
            if (ABS(qa) .lt. qneomult)then
                xkeneo(j) = neomult1*xkineo(j)
            else
                xkeneo(j) = neomult2*xkineo(j)
            endif
            d(nion+1,nion+1,j) = d(nion+1,nion+1,j)+xkeneo(j)
            xketot(j)          = d(nion+1,nion+1,j)
            xchietot(j)        = 2.0 * xketot(j)/(ene(j+1)+ene(j))
         enddo
      endif

c
c -----------------set_chie_chii----------------------------------------
c This model adds a term set_chie_chii*chii to the electron chi.
c Normally you do not want anything additional in chie. But for generallity
c we allow the total chie to include additional terms. Thus we have
c total chie = whateever chie was set above + set_chie_chii*chii
c by doing this here we pick up the total chii from above. HSJ 6/25/98
c ----------------------------------------------------------------------
c
        if (ABS (set_chie_chii) .gt. 1.0e-5) then   ! model is active
          do j=1,nj-1               ! loop over the half grid
            d(nion+1,nion+1,j) = d(nion+1,nion+1,j) +
     .                           set_chie_chii * d(nion+2,nion+2,j)
            xketot(j)          = d(nion+1,nion+1,j)
            xchietot(j)        = 2.0 * xketot(j)/(ene(j+1)+ene(j))
          end do
        end if


      nup  = nion + 2
      if(ddebug(3) .ge. 0.0)then
c
c         average transport coefficients over 2*jsteps+1 intervals
c
          rlxd = 1.0 - ddebug(3)         ! ddebug vector is input from inone

          if (iangrot .eq. 1)  nup = nion + 4
c
          do 3000 i=1,nup
             if (i .eq. nion+3)  go to 3000 !skip rbp 
c
             do 3001 k=1,nup
                if (k .eq. nion+3)  go to 3001
c
                do j=1,nj-1
c
                   jstep =  ABS (ddebug(2))     ! ddebug vector is input from inone
                   jl    = MAX0 (j-jstep  ,  1)
                   ju    = MIN0 (j+jstep+1, nj)         ! symmetric average about j
                   if (ddebug(10) .gt. 0.0)  jl = j     ! forward   average only
                   if (ddebug(10) .lt. 0.0)  ju = j     ! backward  average only
                   dd(j) = 0.0
                   dave  = 0.0
c
                   do js=jl,ju-1
                       if (ddebug(2) .ge. 0.0) then
                          dd(j) = dd(j) + d(i,k,js)
                       else                                 ! smooth logarithm
c
c                         1.0d-100 will be lost due to roundoff for nonzero d
c
                          argument = d(i,k,js) + 1.0d-100
                          dd(j)    = dd(j) + LOG10 (argument)
                       end if
                       dave = dave + 1.0
                   end do
c
                   if      (dave      .gt. 0.0) then
                        dd(j) = dd(j)/dave
                   else if (ddebug(2) .gt. 0.0) then
                        dd(j) = d(i,k,j)
                   else if (ddebug(2) .lt. 0.0) then
                        argument = d(i,k,j) + 1.0d-100
                        dd(j)    = LOG10 (argument)
                   end if
c
                   if (ddebug(2) .lt. 0.0) then
                        dd(j) = 10.0**dd(j)
                        if (dd(j) .lt. 1.0d-80)  dd(j) = 0.0     ! restore 0 value
                   end if
c
                end do
c
                do j=1,nj-1
                   if (dsav(i,k,j) .eq. 0.0) dsav(i,k,j) = dd(j)
                   d(i,k,j) = rlxd*dsav(i,k,j)+(1.0-rlxd)*dd(j)
                end do
 3001        continue
 3000     continue
          do j=1,nj-1
               jstep =  ABS (ddebug(2))
               jl    = MAX0 (j-jstep,1)
               ju    = MIN0 (j+jstep+1,nj)
               dd(j) = 0.0
               dave  = 0.0
               do js=jl,ju-1
                  dd(j) = dd(j) + qexch(js)
                  dave  = dave  + 1.0
               end do
                  dd(j)   = dd(j) / dave
          end do
          do j=1,nj-1
               if (qexsav(j) .eq. 0.0)  qexsav(j) = dd(j)
              qexch(j) = rlxd*qexsav(j)+(1.0-rlxd)*dd(j)
          end do


      else !ddebug(3) .lt. 0.0 branch
           nback = nback+1
           ! back average all d's :
           do  i=1,nup

           do  k=1,nup

c
           do j=1,nj-1
              !this is the average d after nback iterations:
              !where dsav is the average d after nback-1 iterations
              d(i,k,j) = (d(i,k,j) + (nback-1)*dsav(i,k,j))/nback
              !dsav(i,k,j) =d(i,k,j) !for next iteration set on entry above
           enddo !j
           enddo ! i
           enddo !k
      endif !ddebug(3) .lt. 0.0 branch







      do j =1,nj-1
        xketot(j)          = d(nion+1,nion+1,j)
        enea =0.5*(ene(j+1)+ene(j))
        xchietot(j)          = xketot(j)/enea
        xkitot(j)          = d(nion+2,nion+2,j)
c        dmassden = 0.0
        ensum=0.0
        do i=1,nprim    ! have to use nprim to be consistent
****    do i=1,nion     ! with calculation of xkiinv (in PSOURC)
           enav = 0.5 * (en(j,i)+en(j+1,i))
           ensum = enav +ensum
c           dmassden = dmassden+atw(i)*xmassp*enav
        end do
        xchiitot(j)          = xkitot(j)/ensum 
        xkangtot(j)        = d(nion+4,nion+4,j)
      enddo

      if (include_mmm .eq. 1) then !JMP START
        do j =1,nj-1
          if(xchietot(j) .gt. 1.0e5 ) xchietot(j) = 1.0e5
          if(xchiitot(j) .gt. 1.0e5 ) xchiitot(j) = 1.0e5
        end do
      end if !JMP END
c
c --- set up xkangrot for printout 
c

      do j=1,nj-1
        if(iangrot .ne. 0)then
          xkangrot(j) = d(nion+4,nion+4,j)/(dmassden(j)*rmajor**2)
        else
          xkangrot(j) = 0.0 !toroidal rotation excluded by user
        endif
      end do

c
****  write (ncrt, '(a, i7)') ' DIFFUS_COUNTER =', diffus_counter
****  if (ip_chi2_counter .ge. 2951)  call STOP ('DIFFUS: debug', 992)!
      return
c
      end








      subroutine exptlprofiles (timet)
c
      USE param
      USE solcon
      USE numbrs
      USE tordlrot
      USE bd_condtn
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c the purpose of this subroutine is to calculate experimental values
c of the input profiles at time t = time. the experimental values are
c obtained by linear interpolation in time from the input values of the
c profile given at times bctime(j),j = 1,2..nbctim.
c for simulation mode runs this information may not be available since
c it is not required input (in inone) for simulation. In this case
c we simply return a vector (of length nj) with each element equal to 0.
c ----------------------------------------------------------------------
c
c      include 'param.i'
      include 'exptlprof.i'       ! exptlxx profiles will be stored here
      include 'storage.i'
c      include 'bcondspr.i'
c      include 'invers.i'
c      include 'numbrs.i'
c      include 'solcon.i'
c      include 'tordlrot.i'
c
c --- first zero all the experimental profiles
c
      do j=1,nj
        exptlne    (j) = 0.0
        exptlte    (j) = 0.0
        exptlti    (j) = 0.0
        exptlangrot(j) = 0.0
        exptlcur   (j) = 0.0
        excurb_external(j) = 0.0 
      end do
c
c --- save the time
c
      svtim = time
      time  = timet            ! so tsplinew will use timet
c
c --- electron density:
c
      indexx = kprim + kimp + 1           ! set pointer so bparne will..
c                                         ..be used in INTRP

      if (njene .ne. 0)
     .call tsplinew (enein, renein, exptlne, knotsene, indexx,xdum,ydum)
c
c electron temp
c
      indexx = kprim + kimp + 2           ! set pointer so bparte will..
c                                         ..be used in INTRP
      if (njte .ne. 0)
     .call tsplinew (tein, rtein, exptlte, knotste, indexx, xdum,ydum)
c
c ion temp
c
      indexx = kprim + kimp + 3           ! set pointer so bparti will..
c                                         ..be used in INTRP
      if (njti .ne. 0)
     .call tsplinew (tiin, rtiin, exptlti, knotsti, indexx, xdum,ydum)
c
c toroidal angular rotation
c
      if (iangrot .eq. 1) then
        indexx = kprim + kimp + 6        ! set pointer so bparang will..
c                                        ..be used in INTRP
        if (njinang .ne. 0) then
          call tsplinew (angrotin, rangrot, exptlangrot, knotsang,
     .                   indexx, xdum, ydum)
        end if
      end if
c
c current density (experimental current density is not yet available.)
c
****  indexx = kprim + kimp+5 ! set pointer so bpar?? will be used in INTRP
****  if (njcur .ne. 0)
**** .call tsplinew (curdenin,rcurdein,exptlcur,knotscur,indexx,xdum,ydum)

c
c     transp beam driven current:
      if (external_beam_cur  .gt. 0) then
        indexx = kprim + kimp + 8        ! set pointer so bparcurb_external will..
c                                        ..be used in INTRP
        if (njcurb_external  .ne. 0) then
          call tsplinew (curbeam_external, rcurb_external, 
     .                   excurb_external, 
     .                   knotscurb_external,indexx, xdum, ydum)
        end if
      end if
       
c
c --- restore time
c
      time = svtim
      return
c
      end





      subroutine fiziks
c
      USE param
      USE fusion
      USE io 
      USE ions
      USE neut
      USE solcon
      USE soln
      USE mhdpar
      USE mhdgrid
      USE nub
      USE nub2
      USE extra
      USE numbrs
      USE mesh
      USE verbose
      USE sourc
      USE transp,only : use_nubeam
      USE machin
      USE geom
      USE flags
      USE tordlrot
      USE soln2d
      USE bd_condtn,           ONLY  : bctime,ub,fluxb,
     .                                 ub_save,ub_rho_edge,
     .                                 bctime_zone,totcur
      USE rhog
      USE flxav
      USE gpsi
      USE zeffcom
      USE nrtype ,            ONLY :        DP,I4B


      implicit  integer (i-n), real*8 (a-h, o-z)
    
      REAL(DP) taueprev,entaueprev,qtotprev, 
     .         detotprev,etotprev 

      data taueprev,entaueprev,qtotprev, 
     .         detotprev,etotprev /5*0.0_DP/

c
c this subroutine calculates some important physical quantities
c that are needed by subroutines OUT and PLOT
c

      include 'co2.i'
      include 'sxrcom.i'

c
      dimension qeloc(kj), qloc(kj),eneprev(kj),factprev(kj)
      dimension dum(kjp)
      data     pi /3.141592654  /
      data xmassp /1.6726231e-24/



      ksymp2  = 3*nbeams+1                      !for beam_thermal
c
c ----------------------------------------------------------------------
c               quantities at time point n-1/2
c ----------------------------------------------------------------------
c entot : total number of electrons in system
c dentot: time derivative (1/s) of entot
c stot  : total electron particle source (1/s)
c taup  : electron particle confinement time (s)
c etot  : total energy (j) in system
c detot : time derivative (W) of etot
c qtot  : gross energy source (W)
c taue  : energy confinement time (s)
c entaue: confinement parameter (s/cm**2)
c ----------------------------------------------------------------------
c
c set the time
c
      timet = time
      if (n .ne. 0)  timet = time - (1.0 - theta) * dtt
c
c calculate local values for the electron particle source (sione),
c the electron thermal energy (eeloc),
c the plasma thermal energy (eloc), and
c the gross energy source (qloc)
c local (angmtm(j)) and total (angmtot) angular momentum
c moment of inertia density amtinrta:
c
      do 2120 j=1,nj
      eeloc(j)  = 1.5*ene(j)*te(j)
      eloc(j)   = eeloc(j)
      angmtm(j) = 0.0
      do 2100 k=1,nion
      angmtm(j)   = angmtm(j)+atw(k)*xmassp*en(j,k)
 2100 eloc(j)     = eloc(j)+1.5*en(j,k)*ti(j)
      qeloc(j)    = qohm(j) + qbeame(j) + qrfe(j) + qfuse(j)
      amtinrta(j) = angmtm(j)*r2capi(j)*iangrot
      angmtm(j)   = angmtm(j)*iangrot*r2capi(j)*angrot(j)
 2120 qloc(j)     = qohm(j)+qbeame(j)+qbeami(j)+qrfe(j)+qrfi(j)+qfuse(j)
     .             + qfusi(j)
c
c integrate local quantities using trapezoidal rule
c
c ----------------------------------------------------------------------
c source terms are at half time step
c densities and energies are at full time step
c ----------------------------------------------------------------------
c

      call rhomsh(time)

c
      m1 = 1
      if (nneu .eq. 1) m1 = in
c
      do j=1,nj
        fact(j) = volfac * hcap(j)
      end do
c
      call trapv (r,ene,fact,nj,entotn)

      if (fiziksvb .eq. 1) then
           jemax = 0
           enediffmax =0.0D0
           jfmax=0 
           factdiffmax =0.0D0
          if (n .gt. 0) then
             do j=1,nj
                enediffmax1 = ABS (ene(j) - eneprev(j))
                enediffmax  = MAX (enediffmax, enediffmax1)
                if (enediffmax .eq. enediffmax1)  jemax = j
                factdiffmax1 = ABS (fact(j) - factprev(j))
                factdiffmax  = MAX (factdiffmax, factdiffmax1)
                if (factdiffmax .eq. factdiffmax1)  jfmax = j
             end do
          end if
          write (*,'("FIZIKS: time =",1pe14.6," entotn =",1pe14.6  /
     .               "max diff in ene at j =",i5,"value =",1pe14.6 /
     .               "max diff in fact at j=",i5,"value =",1pe14.6)')
     .                time,entotn,jemax,enediffmax,jfmax,factdiffmax
          call copya (ene,eneprev,nj)
          call copya (fact,factprev,nj)
      end if
      call trapv (r,eeloc,fact,nj,eetotn)
      call trapv (r,eloc,fact,nj,etotn)
      call trapv (r,en(1,m1),fact,nj,enitn)
c
      if (iangrot .ne. 0) then
        call trapv (r,angmtm,fact,nj,angmtotn)
        angmtotn = angmtotn*1.0e-7
        call trapv (r,amtinrta,fact,nj,totinrta)
        totinrta = totinrta*1.0e-7
      end if 
c
c --- moment of inertia (totinrta) in kg*m**2
c --- angular momentum             in kg*m**2/sec
c
      cconst = 1.0
      call trap1 (r, fact, cconst, nj, volume)
      call rhomsh (timet)
c
      do j=1,nj
        fact(j) = volfac * hcap(j)
      end do
c
      call trapv (r,sione,fact,nj,stot)
      call trapv (r,qeloc,fact,nj,qetot)
      call trapv (r,qloc,fact,nj,qtot)
****  if (iangrot .ne. 0) then
        call trapv (r, storque , fact, nj, storquet)
        storquet = storquet * 1.0e-7
        call trapv (r, storqueb, fact, nj, beamtorq)
        beamtorq = beamtorq * 1.0e-7
c
c torque in nt-m
c
****  end if
c
      cconst = 1.0
      call trap1 (r, fact, cconst, nj, volumh)
c
c convert energies to MKS units
c
      eetotn = 1.60217733e-16*eetotn
      qetot  = 1.60217733e-16*qetot
      etotn  = 1.60217733e-16*etotn
      qtot   = 1.60217733e-16*qtot
      ecn    = 1.60217733e-16*eloc(1)
      qcen   = 1.60217733e-16*qloc(1)
      eecn   = 1.60217733e-16*eeloc(1)
      qecen  = 1.60217733e-16*qeloc(1)
c
c  obtain values and derivatives of entot, eetot, etot, and ec
c
      if (n .gt. 0)  go to 2200
      entot    = entotn
      eetot    = eetotn
      angmtot  = angmtotn
      etot     = etotn
      ec       = ecn
      eec      = eecn
      dentot   = 0.0
      deetot   = 0.0
      detot    = 0.0
      dec      = 0.0
      deec     = 0.0
      dangmtot = 0.0
      go to 2220
c
 2200 entot    = 0.5 * (entoto   + entotn) ! entoto etc. were set in
      eetot    = 0.5 * (eetoto + eetotn) ! previous call to this routine
      etot     = 0.5 * (etoto    + etotn ) ! so they refer to time step
      ec       = 0.5 * (eco      + ecn   ) ! n-1   note the assumption
      eec      = 0.5 * (eeco     + eecn  ) ! of theta =0.5 here HSJ
      angmtot  = 0.5 * (angmtotn + angmtoto)
c
c     get the central (i.e., pivot) values HSJ 9/10/96
c
      call central_value(entoto,entotn,theta,time,dtt,entot)
      call central_value(eetoto,eetotn,theta,time,dtt,eetot)
      call central_value(etoto,etotn,theta,time,dtt,etot)
      call central_value(eco,ecn,theta,time,dtt,ec)
      call central_value(eeco,eecn,theta,time,dtt,eec)
      call central_value(angmtoto,angmtotn,theta,time,dtt,angmtot)
c
      dentot   = (entotn   - entoto  ) / dtt
      deetot   = (eetotn   - eetoto  ) / dtt
      detot    = (etotn    - etoto   ) / dtt
      deitot_extra   =    detot - deetot
      dec      = (ecn      - eco     ) / dtt
      deec     = (eecn     - eeco    ) / dtt
      dangmtot = (angmtotn - angmtoto) / dtt
c      print *,'angmtotn,angmtoto ,dangmtot=',angmtotn,angmtoto,dangmtot
c      print *,'n,dtt =',n,dtt
c
 2220 entoto   = entotn       ! entoto,..etc. when used above
      eetoto   = eetotn       ! are from the n-1 time step.
      etoto    = etotn        ! now update these quantities to
      eco      = ecn          ! time step n
      eeco     = eecn
      angmtoto = angmtotn
c      print *,'angmtotn,angmtoto ,dangmtot=',angmtotn,angmtoto,dangmtot
c      print *,'n,dtt =',n,dtt
c
c calculate one-turn surface voltage
c
      etorb    = etor(nj)
      voltag   = 2.0 * pi * rmajor*etorb
      voltag0  = 2.0 * pi * rmajor*etor(1)
      voltag25 = 2.0 * pi * rmajor*etor(nj/4)
      voltag50 = 2.0 * pi * rmajor*etor(nj/2)
      voltag75 = 2.0 * pi * rmajor*etor(3*nj/4)
      voltag90 = 2.0 * pi * rmajor*etor(90*nj/100)
      voltag95 = 2.0 * pi * rmajor*etor(95*nj/100)
c
c ----------------------------------------------------------------------
c                quantities at full time point
c ----------------------------------------------------------------------
c
c  calculate volume-average and line-average electron densities
c
      eneav = entotn/volume
      if (codeid .ne. 'onedee')  go to 2310
      call trapl(r,ene,nj,enebar)
      enebar = enebar/r(nj)
      do 2301 k=1,kion
      call trapl(r,en(1,k),nj,enbar(k))
 2301 enbar(k) = enbar(k)/r(nj)
      go to 2320
 2310 call linav(ene,nj,nw,nh,p,pmax,psir,rmhdgrid,zmhdgrid,
     .          ymagn1,enebar)
      do k=1,kion
        call linav (en(1,k),nj,nw,nh,p,pmax,psir,rmhdgrid,zmhdgrid,
     .              ymagn1,enbar(k))
      end do

c
c  calculate various beta values
c
 2320 call rhomsh(time)
      do 2410 j=1,nj
      fact(j)   = 4.0 * pi**2 * rmajor*hcap(j)
 2410 balpha(j) = 8.0 * pi*1.6e-9*walp(j)/(1.5*btor**2)
      call trapv (r,wbeam,fact,nj,ebtot)
      call trapv (r,walp,fact,nj,eatot)
      const  = 8.0 * pi*1.6e-9/(1.5*volume*btgeom**2)
      betae  = const * 0.62415064e16 *  eetotn
      betai  = const * 0.62415064e16 * (etotn-eetotn)
      betab  = const*ebtot
      betaa  = const*eatot
      beta   = betae + betai + betab + betaa
      betan  = ABS(beta*rminor*0.01*btgeom*1.e-4/(totcur(1)*1.e-6))
      betan  = betan*100.
      const  = const*volume
      betae0 = const*eeloc(1)
      betai0 = const*(eloc(1)-eeloc(1))
      betab0 = const*wbeam(1)
      betaa0 = const*walp(1)
      beta0  = betae0 + betai0 + betab0 + betaa0
c
c  calculate betap; use the flux surface average of bp**2
c
c++  use elliptical circumference (cm) if running 'onedee':
c
      if (circum .le. 0.0)
     .    circum = 2.0 * pi * rminor * SQRT ((1.0 + kappa**2) / 2.0)
c
c++  average B poloidal (G):
c
      bpbar = (4.0*pi / 1000.0)*totcur(1)/(circum*0.01)
c
c++  average beta poloidal
c
      betap = beta*(btor/bpbar)**2
c
c  add fast ion energies to total energy
c
      eatot = eatot*1.60217733e-16
      ebtot = ebtot*1.60217733e-16
      etot  = etot + eatot + ebtot
      ec    = ec + 1.60217733e-16*(walp(1)+wbeam(1))
c
c calculate confinement times
c
      taup = 0.0
      xx   = stot-dentot
      if (xx .ne. 0.0)  taup = entot/xx
c
      tauee = 0.0
      xx    = qetot-deetot
      if (xx .ne. 0.0)  tauee = eetot/xx
c
      taue = 0.0
      xx   = qtot-detot
      if (xx .ne. 0.0)  taue = etot/xx
      entaue = entot*taue/volumh
c ----
      if (fiziksvb .eq. 1) then
          if (n .le. 1) then
               taueprev   = taue
               entaueprev = entaue
               qtotprev   = qtot
               detotprev  = detot
               etotprev   = etot
          end if
          write (*,'("FIZIKS: time,timet  =",1pe14.6,2x,1pe14.6   /
     .               " taue =",1pe14.6," taueprev =",1pe14.6      /
     .               " entaue =",1pe14.6," entaueprev = ",1pe14.6 /
     .               " qtot =",1pe14.6," qtotprev = ",1pe14.6     /
     .               " detot =",1pe14.6," detotprev =",1pe14.6    /
     .               " etot =",1pe14.6," etotprev =",1pe14.6)')
     .            time,timet,taue,taueprev,entaue,entaueprev,qtot,
     .            qtotprev,detot,detotprev,etot,etotprev
          entaueprev=entaue
          taueprev=taue
          qtotprev=qtot
          detotprev=detot
          etotprev=etot
      end if
c
      tauec = 0.0
      xx    = qcen-dec
      if (xx .ne. 0.0)  tauec = ec/xx
c
      taueec = 0.0
      xx     = qecen-deec
      if (xx .ne. 0.0)  taueec = eec/xx
c
      tauangt = 0.0
      xx      = beamtorq-dangmtot
      if (xx .ne. 0) tauangt = angmtot*iangrot/xx
c      print *,'tauangt,beamtorq,dangmtot,angmtot =',
c     . tauangt,beamtorq,dangmtot,angmtot
c       print *,'angmtoto,n,dtt 2  =',angmtoto,n,dtt
c      print *,'angmtotn,angmtoto =',angmtotn,angmtoto 
c
c ----------------------------------------------------------------------
c  calculate bpoloidal
c ----------------------------------------------------------------------
c
      bpol(1) = 0.0
      do j=2,nj
        bpol(j) = rbp(j)/(r(j)*fcap(j)*gcap(j)*hcap(j))
      end do
c
c  calculate total current
c
      totcur(1) = 5.0*rbp(nj)/fcap(nj)
c
c safety factor (PAP 10/89); quadratic extrapolation to axis
c     limit maximum q to 1000
c
      do 3600 j=2,nj
 3600    q(j) = r(j)*btor/(rmajor*bpol(j))
      q(1) = (q(2)*r(3)**2-q(3)*r(2)**2)/(r(3)**2-r(2)**2)
c
      do j=1,nj
        absq = ABS (q(j))
        if (absq .ge. 1.0e3) q(j) = 1.0e3*q(j)/absq
      end do
c
c current density (toroidal) (PAP 10/89)
c
      do 3001 j=1,nj
 3001    dum(j) = rbp(j)/fcap(j)
      do 3002 j=2,nj-1
         curden(j) = (1.0/(0.4 * pi * hcap(j)*r(j)))
     .     *(-dum(j-1)*(r(j+1)-r(j))**2+dum(j)*((r(j+1)-r(j))**2
     .      -(r(j)-r(j-1))**2)+dum(j+1)*(r(j)-r(j-1))**2)
     .      /((r(j+1)-r(j))*(r(j)-r(j-1))*(r(j+1)-r(j-1)))
 3002 continue
c      curden(nj) = (1.0/(0.4 * pi * hcap(nj)*r(nj)))
c     .      * (dum(nj-2)*(r(nj)-r(nj-1))**2-dum(nj-1)*(r(nj)-r(nj-2))**2
c     .      + dum(nj)*(r(nj-1)-r(nj-2))*(2.0*r(nj)-r(nj-1)-r(nj-2)))
c     .      / ((r(nj)-r(nj-1))*(r(nj-1)-r(nj-2))*(r(nj)-r(nj-2)))
c      curden(1) = (1.0 / (0.4 * pi * hcap(1))) * 2.0
c     .          * (dum(3) * r(2) - dum(2) * r(3))
c     .          / ((r(3) - r(2)) * r(3) * r(2))
      call curcalc(rbp,fcap,hcap,gcap,r,curden,curpar_soln,nj)
c
c  calculate li from Bpoloidal - unchecked
c
      xsum1 = 0.0
      xsum2 = 0.0
c
      do j=1,nj-1
        zra = (r(j)+r(j+1))/2.0
        zdra = r(j+1)-r(j)
        zbpa = (bpol(j)+bpol(j+1))/2.0
        zhcapa = (hcap(j)+hcap(j+1))/2.0
        zgcapa = (gcap(j)+gcap(j+1))/2.0
        zdvoln = zra*zdra*zhcapa
        xsum1 = xsum1+zdvoln*zgcapa*zbpa*zbpa
        xsum2 = xsum2+zdvoln
      end do
c
      xlica = xsum1/bpbar/bpbar/xsum2
c
c     calculate the number of neutrons produced by D-D and D-T fusion
c
****  call ddfust (ti, en, ddfusn, ddntot, time, bpol, totcur(1))
c
c --- New fusion calcs, 1/25/96 HSJ 
c
      if (fusionvb .gt. 0)  write (*, '(" calling THERMONUCLEAR_RATE")')
 
   
      call thermonuclear_rate (ddfusn,ddpfus,dtnfus,ttnfus,
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
 

       if (ibeam .gt. 1 .and. .not. use_nubeam) then


       if (fusionvb .gt. 0)  write (*, '(" calling BEAM_THERMAL_FUS")')
       call beam_thermal_fus (time,timmax, bpol, totcur(1),qbfus,sbfus)
       if (beam_thermal_long_calc .lt. 0)
     .   call beam_thermal_approx_fus (time, timmax, qbfus, sbfus)
       call copya(beam_thermalddn(1,ksymp2),ddbeam,nj)
       if (fusionvb .gt. 0)
     .   write (*, '(" calling BEAM_BEAM_RATE from FIZIKS")')
       call beam_beam_rate (time, time0, timmax, totcur(1))
       if (beam_beam_long_calc .lt. 0)
     .   call beam_beam_approx_rate(time,time0,timmax)
c       if (fusionvb .gt. 0)  write (*, '(" beam-thermal d(d,n)he3",
c     .                                " from DDFUST =", 1pe14.3)') ddbmt
       if (fusionvb .gt. 0)  write (*, '(" d(d,n)he3 from ",
     .                "BEAM_THERMAL_FUS", 1pe14.3)') beam_thermal_ddntot

       ELSE if(use_nubeam)then
           call copya(beam_thermalddn(1,ksymp2),ddbeam,nj)
       end if
c
c --- total d(d,n)he3 neutron rate:
c
       ddntot = thermal_thermal_ddntot + beam_thermal_ddntot
     .          + beam_beam_ddntot
c       print *,'ddntot,thermal_thermal_ddntot ='
c       print *,ddntot,thermal_thermal_ddntot
c       print *,beam_thermal_ddntot ,beam_beam_ddntot
c       call stop('fiziks',1)
       ddnthm = thermal_thermal_ddntot
c
c --- total d(t,n)he4 rate
c
       dtntot = thermal_thermal_dtntot + beam_thermal_dtntot
     .          + beam_beam_dtntot
c
c --- total t(t,2n)he4 rate:
c
       ttntot = thermal_thermal_tt2ntot + beam_thermal_tt2ntot
     .          + beam_beam_tt2ntot
c
c --- total d(d,p)t rate
c
       ddptot = thermal_thermal_ddptot + beam_thermal_ddptot
     .          + beam_beam_ddptot
c
c
c --- total He3(d,p)He4 rate
c     NOTE: beam_beam_he3 is not programmed since we assume single species beam
      beam_beam_he3dptot =0.0D0
      he3dptot = thermal_thermal_he3dptot + beam_thermal_he3dptot
     .           + beam_beam_he3dptot

       qdd = ddntot*3270.0               ! 3.27 MeV for d(d,n)he3
     .     + ddptot*4030.0               ! 4.03 MeV for d(d,p)t
c
       qdt = dtntot*17600.0              ! 17.6 MeV for d(t,n)he4
       qtt = ttntot*11300.0              ! 11.3 MeV for d(d,n)he3
       qdd = qdd*1.602e-16 ! watts, divide by input power elsewhere
       qdt = qdt*1.602e-16
       qtt = qtt*1.602e-16
       qhe3d = he3dptot *1.602e-16
       qhe3d = qhe3d*18300.0             ! 18.3 mev for  he3(d,p)He4

       call getpow (time, pwrinpt)
       if (pwrinpt .gt. 1.0) then
         qdd    = qdd/pwrinpt
         qtt    = qtt/pwrinpt
         qdt    = qdt/pwrinpt
         qhe3d  = qhe3d/pwrinpt
       else
          qdd   = 0.0
          qtt   = 0.0
          qdt   = 0.0
          qhe3d = 0.0
       end if

c
c ----------------------------------------------------------------------
c     calculate sxr, co2 and zeff diagnostic information
c ----------------------------------------------------------------------
c

      if (jsxr .ne. 0)
     .  call sxrcal (codeid, kappa, ene, jsxr, nj, nw, nh,
     .               p, pmax, psir, r, rmajor, rin, rmax, te, rmhdgrid,
     .               zmhdgrid, zax, zmin, zmax,
     .               idiode, narray, namar, ndiode, roamin, sxr)
      if ( jco2 .ne. 0)  call co2cal
      if (jzeff .ne. 0)  call zefcal

      return

      end

      subroutine fluxx
c
      USE param
      USE io 
      USE ions
      USE soln
      USE numbrs
      USE mesh
      USE sourc
      USE machin
      USE geom
      USE flags
      USE tordlrot
      USE tcoef
      USE bd_condtn
      USE common_constants,                ONLY : zeroc

      USE paleocl,                         ONLY : include_paleo,
     .                                            paleo_mult,
     .                                            paleo_pinch_mult

      USE flx
      USE tordlrot,                        ONLY : aniwdwdt,wdnidt,
     .                                            qomegapi,qangce  ! HSJ 4/7/14
      USE glf23!jmp.snu
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c     This subroutine calculates the particle, energy and angular
c     momentum fluxes, as well as the flux in Faraday's Law related to
c     the electric field.
c

c      include 'flx.i'

      include 'storage.i'

c
      dimension              dflux(kj),fluxp(kjm1),dnidtm(kj),
     .                       qdelt_l(kj)
      equivalence  (xdum(1), dflux(1))
      real*8        kevperg
      data          kevperg /6.2415064e+08/
      data          xmassp, five_halfs /1.6727231e-24, 2.5/



c
c  initialize fluxes
c

      call zeroa (flux,kk*kjm1)
      call zeroa (fluxi,nj-1)
      call zeroa (dnidtm,nj-1)
      five_halfs_te = five_halfs
      five_halfs_ti = five_halfs
      if (no_te_convection .eq. 1)  five_halfs_te = 0.0
      if (no_ti_convection .eq. 1)  five_halfs_ti = 0.0
      if (no_te_convection .eq. -1)  five_halfs_te = 1.5
      if (no_ti_convection .eq. -1)  five_halfs_ti = 1.5
c
      if (iangrot .ne. 0) then
         do 45 j=1,nj-1
            do 45 k=1,nprim
               dnidtm(j) = dnidtm(j)+dnidt(j,k)*xmassp*atw(k)
     .                   * 0.5 * (angrot(j)+angrot(j+1))
     .                   * 0.5 * (r2capi(j)+r2capi(j+1))
   45   continue
      end if
c
c  calculate particle fluxes
c

      do 300 i=1,nprim
      if (ichiv_model .eq. 1) go to 998 !jmp.snu
      if (itran(i) .eq. 0)  go to 200 !jmp.den
      !jmp.snu do 100 j=1,nj-1
  998 do 100 j=1,nj-1 !jmp.snu
      ena = 0.5 * (en(j,i)+en(j+1,i))
      flux(i,j) = flux(i,j) + ena*vpinch(j)
      do 100 k=1,nk
  100 flux(i,j) = flux(i,j) - d(i,k,j)*dudr(k,j)
      go to 300
c
c --- for analysis mode variables
c
  200 do 220 j=1,nj
      dflux(j) = s(i,j) - dnidt(j,i)

  220 if (itimav .eq. 1)  dflux(j) = dflux(j) + ssaw(j,i)
      call flxcal (dflux,drr,hcap,nj,r,ra,fluxp)
      do 260 j=1,nj-1
cjmp.den      do 260 j=45,nj-1	
      flux(i,j) = fluxp(j)
      do 240 k=1,nk
  240 d(i,k,j) = 0.0
  260 if (dudr(i,j) .ne. 0.0) d(i,i,j) = -flux(i,j)/dudr(i,j)
c
  300 continue    ! end primary ion do loop
c
      IF  (no_ti_convection .eq. 1) THEN
         call zeroa (fluxi,nj-1)
      ELSE
         DO j=1,nj-1 
            DO i =1,nprim
              fluxi(j) = fluxi(j) + flux(i,j) !convective total, due to all primary ions
            ENDDO
         ENDDO
      ENDIF
       
c     
c
c --- get the electron particle flux on the half grid:
c

      IF (no_te_convection .eq. 1)THEN
         call zeroa (fluxe,nj-1)
      ELSE
         do j=1,nj
            dflux(j) = sione(j) - dnedt(j)
            if (itimav .eq. 1)
     .           dflux(j) = dflux(j) + ssawe(j)
         end do
         call flxcal (dflux,drr,hcap,nj,r,ra,fluxe)
      ENDIF




c
c     calculate energy fluxes and flux in Faraday's law
c     (note that there should be an additional term in the Faraday's law
c     flux proportional to the fast ion density gradient, but it is
c     included as a source term instead (PAP 10/89))
c
      

      DO  j=1,nj  ! may 2d terms to qdelt_l, see impsrc
         qdelt_l(j) = qdimpl(j)*(u(nk-iangrot-2,j)-u(nk-iangrot-1,j))
      ENDDO

      i = nion+1
      IF(itran(i) == 0)THEN                 ! analysis mode for te 
         dflux(:) = s(i,:) - dpedt(:) - qdelt_l(:)
         CALL  flxcal (dflux,drr,hcap,nj,r,ra,fluxp)
         flux(i,:)  = fluxp(:)              ! total flux in analysis mode
      ELSE                                  ! simulation  mode for te  
         DO  j=1,nj-1
            tea = 0.5 * (te(j)+te(j+1))
            enea =0.5 * (ene(j)+ene(j+1))
            if(include_glf .eq. 1)then
               flux(nion+1,j) = flux(nion+1,j) + ve_eff_m(j)*enea*tea
            elseif(include_paleo .eq. 1 )then
               flux(nion+1,j) = flux(nion+1,j) + tea*( five_halfs_te *
     .              fluxe(j))   ! +    paleo_mult*paleo_pinch_mult * paleo_flux(j))    
            else
               flux(nion+1,j) = flux(nion+1,j) + five_halfs_te * 
     .                                               tea * fluxe(j)
            endif
            DO k=1,nk-iangrot
               flux(i,j) = flux(i,j) - d(i,k,j)*dudr(k,j) 
            ENDDO
          ENDDO
       ENDIF

      i = nion+2
      IF(itran(i) == 0)THEN                   ! analysis mode for ti
         dflux(:) = zeroc
         DO l=1,nion
            dflux(:) = dflux(:)+ dpidt(:,l)
         ENDDO

         call rot_terms

         dflux(:) = 1.5*dflux(:) + wdnidt(:)+aniwdwdt(:)  ! hsj 4/7/14

 
         dflux(:) = s(i,:) - dflux(:) + qdelt_l(:) 


         CALL  flxcal (dflux,drr,hcap,nj,r,ra,fluxp)
c        note that flux(nion+2,:) = condi+convi+ omegapi+flxangce

         flux(i,:) = fluxp(:)

      ELSE
         DO  j=1,nj-1
            tia = 0.5 * (ti(j)+ti(j+1))
            if(include_glf .eq. 1 )then

                enpa =0.0
             do i =1 ,nprim
                enpa = enpa + (en(j+1,i)+en(j,i))*0.5
             enddo
             flux(nion+2,j) = flux(nion+2,j) + vi_eff_m(j)*enpa*tia

            else
             flux(nion+2,j) = flux(nion+2,j) + five_halfs_ti * 
     .                                                tia * fluxi(j)
            endif
            DO k=1,nk-iangrot
               flux(i,j) = flux(i,j) - d(i,k,j)*dudr(k,j) 
            ENDDO
         ENDDO
      ENDIF




      i = nion+3                                 ! Faraday's law
      DO  j=1,nj-1
         IF(itran(i) == 0)THEN    ! analysis mode
            flux(i,j) = zeroc     ! not implemented
         ELSE
            DO k=1,nk-iangrot
               flux(i,j) = flux(i,j) - d(i,k,j)*dudr(k,j) 
            ENDDO
         ENDIF
      ENDDO


c
c calculate fluxes at the boundary using linear extrapolation
c

      do k=1,nk-1-iangrot
        call extrap (ra(nj-2), ra(nj-1), r(nj),
     .               flux(k,nj-2), flux(k,nj-1), fluxb(k))
      end do
      if (no_te_convection .eq. 1)  fluxb(nion+1) = 0.0
      if (no_ti_convection .eq. 1)  fluxb(nion+2) = 0.0




c
c --- angular momentum flux calculations
c
      do j=1,nj
        flxangce(j) = 0.0
        omegapi(j) = 0.0
      end do
c
      momtfl:   if (iangrot .ne. 0) then
        i = nk
c
c ---   convective part of angular momentum flux
c
        do 650 j=1,nj-1
          angrota = 0.5 * (angrot(j)+angrot(j+1))
          r2capia = 0.5 * (r2capi(j)+r2capi(j+1))
          fluxangc(j) = 0.0
          do 660 k=1,nprim
  660       fluxangc(j) = fluxangc(j)+atw(k)*flux(k,j)
          fluxangc(j) = fluxangc(j)*xmassp*r2capia*angrota
  650   flux(i,j) = fluxangc(j)
        if (itran(i) .eq. 1) then
c
c           simulation mode, get the flux
c
            if(include_glf .eq. 1 )then
              do j=1,nj-1
                 do k=1,nk
                    flux(i,j) = flux(i,j)-d(i,k,j)*dudr(k,j)
                 end do
                 fluxangv(j) = flux(i,j)-fluxangc(j)+ vrot_eff_m(j)*.0
              end do
            else
              do j=1,nj-1
                 do k=1,nk
                    flux(i,j) = flux(i,j)-d(i,k,j)*dudr(k,j)
                 end do
                 fluxangv(j) = flux(i,j)-fluxangc(j)
              end do
            endif
c
        else
c
c         analysis mode first get the flux (given the divergence)
c

          do 1200 j=1,nj-1
 1200     dflux(j) = s(i,j)-dnidtm(j)
          call flxcal (dflux,drr,hcap,nj,r,ra,fluxp)
          do 1300 j=1,nj-1
            flux(i,j) = fluxp(j)
            fluxangv(j) = fluxp(j)-fluxangc(j)
c
c ---       now get the diagonal diffusion coefficient
c
            do 1400 k=1,nk
 1400         d(i,k,j) = 0.0
 1300       if (dudr(i,j) .ne. 0.0)  d(i,i,j) = -fluxangv(j)/dudr(i,j)
        end if
c
c ---   correction to the ion energy flux
c
        if (angrcple .ne. 0.0 .AND. itran(nion+2) == 1) then
            do j=1,nj-1
              angrota     = (angrot(j)+angrot(j+1))*0.5
              omegapi(j)  = angrcple*kevperg*angrota*fluxangv(j)
              flxangce(j) = 0.5 * angrcple*angrota*kevperg*fluxangc(j)
              flux(nion+2,j) = flux(nion+2,j)+omegapi(j)+flxangce(j)
            end do
c
c ---       get new extrapolated edge value
c
            call extrap (ra(nj-2),ra(nj-1),r(nj),
     .      flux(nion+2,nj-2),flux(nion+2,nj-1),fluxb(nion+2))
c
c ---       extrapolate omegapi and flxangce
c
            call extrap (ra(nj-2),ra(nj-1),r(nj),
     .                   omegapi(nj-2),omegapi(nj-1),omegapi(nj))
            call extrap (ra(nj-2),ra(nj-1),r(nj),
     .                   flxangce(nj-2),flxangce(nj-1),flxangce(nj))
        end if
c
c ---   extrapolation of momentum flux
c
        call extrap (ra(nj-2),ra(nj-1),r(nj),
     .                          flux(i,nj-2),flux(i,nj-1),fluxb(i))
      end if    momtfl

c
c
c
c --- calculate divergence of electron and ion energy flux to get
c --- the source that would be required in order for the given
c --- te and ti profiles to be in steady state. these calculations
c --- are not used elsewhere in the code
c
c      call divflx(flux,fluxb,hcap,r,ra,drr,nion+1,nj,kk,ydum)   ! for te
c      call divflx(flux,fluxb,hcap,r,ra,drr,nion+2,nj,kk,zdum)   ! for ti
 
      return
c
      end





      subroutine hirsh88 (jhirsh, j, nion, fcapa, enea, tea, ena, tia,
     .                    zeffa, bpa, cee, eta, ra, rmajor, enbeam,
     .                    enalp, sumdz, za, wneo, xia, xhm2a, codeid)
c                         d,dfast, dfion, xden, xdit, xdet, xdin)
c
c ----------------------------------------------------------------------
c  this subroutine sets up the diffusion coefficients so that the
c  bootstrap current model given by hirshman,phys fluids,vol31(10),1988
c  will be used.
c  eq 10 of hirshman is transformed into the following form
c  for use in ONETWO:
c    <jboot dot b/bto> = a*{b*dni/drho+c*dte/drho
c                           +d*dti/drho +e*dnf/drho}
c  where
c      <jboot dot b/bto> has units of [keV/ cm**4*gauss]
c        a = -ne*tea*x/(fcapa*bp0)      {keV/(cm**3*gauss)}
c        x = fraction trapped/fraction circulating
c        b = l31*tia/(zi*tea*ni) + l31*zi/ne
c        c = l31/tea + l32/tea + l31*zi/ne
c        d = l31/(z*tea)3+/l32*alphai/(zi*tea)
c        e = l31*zf/nea+f
c        f = l31*tia/(zi*tea*ni)
c    l31 and l32 as defined here are dimensionless numbers
c    note that l31 = L31/(J0*X) and l32 =L32/(J0*X)
c    where  L31 and L32 are given by eq 11a,11b in hirshman
c    and J0*X   =  -f(psi)*ne*te*x , (a=J0*X/(R0*Bt0*Bp0, f(psilim) = R0*Bt0,
c    and fcapa  =   f(psilim)/f(psi)  )
c    note that electron density gradients are written in terms
c    of equivalent ion density gradients (since ONETWO does not
c    have a separate electron density equation).  Term e is the
c    electron density contribution from fast ions and (term f) the
c    fast ion density contribution.
c    Term f results from assuming that the fast ion density
c    can be treated as a thermal species.  Jhirsh = 88 sets term f
c    to zero, jhirsh = 89 includes term f.
c    This model assumes a pure electron,single ion plasma.
c    Thus we interpret ni as the sum over all ion species,zi
c    as zeff.  We break the individual
c    ion species density gradients out to be comensurate with the
c    rest of the code.  However it is only the sum that has meaning
c    in this model.
c
c    in ONETWO the bootstrap is calculated as
c    <jboot dot b/bto>= -s*(sum over ions {d(kfar,i,j)*dni/drho )
c                           +d(kfar,nion+1,j)*dte/drho
c                           +d(kfar,nion+2,j)*dti/drho
c                           +dcoef(j)*dnf/drho  )
c    where
c    s = rho/(cee*eta)  (a pure number)
c    sum over ions is sum from i = 1 to i=nion,and includes
c    primary as well as impurity species.
c    kfar = nion+3, indexes Faraday's law.
c    dcoef is a special coefficient used to model electron
c    density due to fast ions and fast ion density (on option)
c    it is calculated in subroutine SOURCE as dcoef = dfast+dfion
c    dfast is always used whereas dfion,the contribution from
c    fast ions,will be zero if jhirsh = 88.
c
c --- input through argument list:
c
c  jhirsh                selects bootstrap current model.
c                        jhirsh=88 means do not include fast ion
c                        density gradient in calculation.
c                        jhirsh = 89 means inlcude this term (it is
c                        included in an ad hock manner,which may
c                        not yield the correct fast ion term ).
c
c  j                   is the grid index. note that we are working
c                      on the half grid here so j ranges from 1 to nj-1
c
c  nion                number of ions (primary plus impurity)
c
c  fcapa               f(psilim)/f(psi) on the (half) grid point j
c
c  enea                electron density, #/cm**3, at j
c
c  tea                 electron temp, keV, at j
c
c  ena(i)              i = 1,2..nion,density,#/cm**3 at j
c
c  tia                 common ion temp, keV, at j
c
c  zeffa               zeff at j
c
c  bpa                 poloidal b field,gauss, at j
c
c  cee                 speed of light,cm/sec
c
c  eta(j)              resistivity,sec,on half grid
c
c  ra(j)               rho,cm,on half grid
c
c  enbeam(j)           beam density,on full grid,#/cm**3
c
c  enalp(j)            alpha density,#/cm**3,on full grid
c
c  sumdz               sum of (ni*d<z>/dte)*te/ne at j
c
c  za(i)               i = 1,2,..nion z of ion i at j
c
c  wneo                neoclassical multiplier array
c
c  xia(3,3)            takes the place of xi33 for trapped
c                      particle calculations if codeid = 'onedee'
c
c  xhm2a               takes the place of xhm2 for trapped
c                      particle calculations if codeid = 'onedee'
c
c  codeid              distinguishes between onedee and twodee runs
c
c --- input from include files
c    param.i:
c       kj               size of rho grid
c
c    neo2d.i:
c       xi33(j)          j = 1,2,...nj, see eq. 4.2-58 of gaa16178
c       xhm2(j)          the exact trapped particle fraction is
c                        given by ft = 0.75*(xi33(j)/xhm2(j))
c                        (if codeid = 'twodee')
c
c --- output:
c  through argument list:
c
c diffusion coefficients are defined so that each one makes
c a contribution of the form -s*d(kfar,i,j) to <jboot dot b/bto>
c see above for definition of s
c  d(kfar,i,j)          i = 1,2,..nion ,statamp.cm**2
c                       gives density contribution at half grid pt j
c  d(kfar,i,j)          i = nion+1 diffusion coefficient,statamp/(cm keV)
c                       gives te contribution at half grid pt j
c  d(kfar,i,j)          i = nion+2 diffusion coefficient,statamp/(cm keV)
c                       gives ti contribution at half grid point j
c
c  dfast(j)             electron "diffusion" coefficient for electrons
c                       sourced by the fast ions. It is assumed that zbeam = 1
c                       dfast has units of statamp cm**2
c                       the contribution to <jboot dot b/bto> is
c                                -dfast*(dnf/drho)/xres
c                       A SOURCE TERM
c                            (1.0/(HCAPA*RHO)*(D/DRHO)FA
c                       (WHERE
c                            FA = (HCAPA*RHO*DFAST*ZBEAM*(DNF/DRHO)/XRES
c                        AND DNF/DRHO IS THE FAST ION DENSITY GRADIENT )
c                        IS ADDED TO THE RHS OF FARADAY'S LAW
c
c  dfion(j)              fast ion "diffusion" coefficient.
c                        dfion has units of statamps cm**2
c                        dfion is set to zero if jhirsh = 88
c                        otherwise,if jhirsh = 89, then
c                        the contribution to <jboot dot b/bto> is
c                            -dfion*(dnf/drho)/xres
c                       A SOURCE TERM
c                            (1.0/(HCAPA*RHO)*(D/DRHO)FA
c                       (WHERE
c                            FA = (HCAPA*RHO*DFION(DNF/DRHO)/XRES
c                        AND DNF/DRHO IS THE FAST ION DENSITY GRADIENT )
c                        IS ADDED TO THE RHS OF FARADAY'S LAW
c                        also if dfion is used (because jhirsh = 89)
c                        then the total ion density includes the fast ions
c
c  ancillary arrays:
c
c  xden(j)                 is defined so that ion species i contributes
c                               -xden*zi*(dni/drho)/xres
c                          to <jboot dot b/bto>
c                          (accounts for electron density gradients)
c
c  xdin(i,j)               is defined so that ion species i contributes
c                               -xdin*(dni/drho)/xres
c                          to jboot dot b/bto>
c                          (accounts for ion density gradients)
c
c  xdet(j)                 is the same as d(kfar,nion+1,j)
c
c  xdit(i,j)               is defined so that ion species i contributes
c                               -xdit*(dti/drho)/xres
c                          to <jboot dot b/bto>
c
c --- output to include files:
c
c  neo2d.i:
c      ftfc(j)              j = 1,2..nj-1 the trapped/circulating
c                           particle ratio
c      h88l31(j)                          l31
c      h88l32(j)                          l32 saved for comparison with
c                                         other models
c
c ------------------------------------------------------------------ HSJ
c
      USE param
      USE neo2d
      USE tcoef, ONLY :                          dfast, dfion, xden,
     .                                           xdit, xdet, xdin,d
      implicit none
c
      integer j,kfar,l,nion,jhirsh
c
      real*8  a31,a32,alphai,df,l31,l32,bpa,zeffa,ena,xunits,
     .        cee,eta,ra,nia,wneo,enf,rmajor,
     .        za,enalp,enbeam,enea,fc,ft,fcapa,sumdz,tea,
     .        tia,effz,
     .        xres,x,z,a,const,denom,xia,xhm2a
c
      character*8 codeid
c
c      include 'param.i'
c      include 'neo2d.i'
c
      dimension eta(*),ra(*),ena(*),enbeam(*),enalp(*),
     .          xia(3,*),wneo(5,*),za(*)
c
c --- statement functions (x = trapped particle fraction, z = charge #,
c --- see equations 11a,11b,11c,and equation for D(X) ):
c
      a31(x,z)  = ((z+2.21)*z+0.754)+x*((z+1.243)*z+0.348)
      a32(z)    =  0.884+2.074*z
      alphai(x) = -1.172/(1.0 + 0.462*x)
****  df(x,z)   = (1.414+z)*z+(x*((z+1.243)*z+3.48) ! fixed 8/25/99 HSJ
      df(x,z)   = (1.414+z)*z+(x*((z+1.243)*z+0.348)
     .          + ((2.0*z+2.657)*z+0.754))*x
c
c  we have 1.602e-8 amps/cm**2 = 1 keV/(cm**4 gauss)
c          1 amp = 2.998e9 statamps
c  multiplication by xunits changes keV/(cm**4 gauss) to statamps/cm**2
c
      xunits = 1.602e-8 * 2.998e9
      kfar   = nion + 3
      if (codeid .eq. 'onedee') then
        ft = 0.75 * xia(3,3) / xhm2a
      else
        ft = 0.75*(xi33(j+1)+xi33(j))/(xhm2(j+1)+xhm2(j))
      end if
c
      fc      = 1.0 - ft             ! fraction of particles circulating
      x       = ft / fc              ! trapped/circulating ratio
      if (jhirsh .eq. 87)  x = SQRT (ra(j)/rmajor) * 1.46
      ftfc(j) = x                    ! save for printout
      a = xunits * enea * tea * x /
     .           (fcapa * bpa)       ! see def of A above, statamps/cm
      xres = cee * eta(j) / ra(j)    ! pure number
c
c --- get the total thermal density
c
      nia = 0.0
      do l=1,nion
        nia = nia + ena(l)
      end do
c
c --- add fast ion density if called for
c
      IF (jhirsh .eq. 89) then
        enf = (enbeam(j+1) + enbeam(j)) * 0.5
        enf = enf + enalp(j+1) + enalp(j)    ! z = 2 for alpha particles
        nia = nia + enf
      ELSE
         enf = 0.0
      ENDIF
      denom =   df(x,zeffa)
      l31   =  a31(x,zeffa) / denom
      l32   = -a32(  zeffa) / denom
      h88l31(j) = l31*x
      h88l32(j) = l32*x
     
****  effz  = zeffa       ! for strict Hirshman 88
      effz  = enea / nia  ! to account for multiple ion species (approx)
c
c --- electron temp gradient term(units on xdet are statamp/(cm*keV))
c --- the electron temp grad contribution to jboot is -xdet/xres
c --- see subroutine SOURCE
c
      if (jhirsh .ne. 87) then
        xdet(j) = a*xres*(l31*(1.0 + sumdz)+l32)/tea
        xdet(j) = xdet(j)*wneo(4,2)
      else
        xdet(j) = a*xres*0.608/(1.46*tea)
        xdet(j) = xdet(j)*wneo(4,2)
      end if
      d(kfar,nion+1,j) = xdet(j)
c
c --- ion temperature gradient term
c --- the ion temp grad contribution to jboot is -sum over i of
c                                              {xdit(i,j)/xres}
c --- see subroutine SOURCE. xdin has units of statamp/(cm keV)
c
      if (jhirsh .ne. 87) then
        xdit(1,j) = (a*xres*l31/(effz*tea))*(1.0 + alphai(x))
        xdit(1,j) = xdit(1,j)*wneo(4,3)
      else
        xdit(1,j) = -a*xres*0.412/(1.46*tea)
        xdit(1,j) = xdit(1,j)*wneo(4,3)
      end if
      d(kfar,nion+2,j) = xdit(1,j)
c
c --- ratio out the total so that each ion species receives
c --- the fraction of the total density as its bootstrap contribution.
c --- note that this is arbitrary and is done only for printout purposes
c
      do l=1,nion
        xdit(l,j) = d(kfar,nion+2,j)*ena(l)/nia
      end do
c
c --- add fast ion contribution to primary ion #1
c
      xdit(1,j) = xdit(1,j) + d(kfar,nion+2,j) * enf / nia
c
c --- the electron density gradient term. It is used only to split out
c --- the electron density contribution in the printout.  Otherwise
c --- it is included in terms of the equivalent ion density terms.
c --- xden has units of statamps cm**2
c
      xden(j) = a * xres * l31 / enea
      xden(j) = xden(j) * wneo(4,1)
      if (jhirsh .eq. 87)
     .xden(j) = a * xres * 2.397 / (1.46*enea) * wneo(4,1)
c
c --- the (thermal) ion density gradient terms. only the total,single ion
c --- fluid, model is physical. the split into the various ion
c --- species terms is not significant in this model and is done
c --- only to be comensurate with the remainder of the code.
c --- the ion den grad contribution to jboot is -sum over i of {
c                                              xdin(i,j)/xres  }
c --- see subroutine SOURCE. xdin has units of statamps cm**2
c
      const = a *  l31 * xres *   tia / (effz * tea * nia)
      if (jhirsh .eq. 87)
     .const = a * xres *  tia * 2.397 / (1.46 * tea * nia)
      do l=1,nion
        xdin(  l,j) = const
        xdin(  l,j) = wneo(4,1)*xdin(l,j)
        d(kfar,l,j) = xdin(l,j)+xden(j)*za(l)
      end do
c
c --- density gradient term of electrons due to fast ions
c --- the bootstrap current contribution is -(xden/xres)*zbeam
c --- in units of statamps/cm**2 . It is assumed that zbeam = 1
c
      dfast(j) = wneo(4,1) * xden(j)
c
c --- finally, the fast ion gradient term.  It is assumed that the
c --- fast ion term can be treated as part of the thermal species.
c --- the bootstrap current contribution is -(dfion/xres)
c --- in units of statamps/cm**2
c
      dfion(j) = 0.0
      if (jhirsh .eq. 89) then
        dfion(j) = wneo(4,1)*a*xres*l31*tia/(effz*tea*nia)
      end if
      return
c
      end

      subroutine neoresist (zeffa, wneot, xft, xnue, xnuse, xmasse,
     .                      enea, charge, fftrap, etaspitzer,
     .                      etaparallel, xkap11)
c
c ----------------------------------------------------------------------
c subroutine calculates the neoclassical resistivity with finite
c inverse aspect ratio and collisional effects.
c see Hirshman et al.,N.F. 17,3,(1977)pg. 611
c                     N.F. 21,9,(1981)pg. 1158,eq. 7.41
c  Note that this formulation is valid only for zeff< 3.0
c  for zeff > = 3.0 we set cr (the conductivity reduction due to
c  e-e collisions) equal to 0.0
c
c --- input
c
c  zeffa        zeff on half grid
c  xft          trapped particle fraction (collisionless)
c  wneot        multiplier,normally set to 1.0,for effective trapped
c               electron fraction.
c  xnue         e-e collision frequency (includes zeff term)
c  xnuse        collision/bounce frequency (includes zeff)
c  xmasse       electron mass,grams
c  enea         electron density, cm-3, on half grid
c  charge       electron charge, esu
c
c --- output
c  fftrap        effective trapped electron fraction
c
c  etaparallel   parallel resistivity (in sec)
c
c  etaspitzer    Spitzer resistivity
c  xkap11        zeff dependent factor of etaspitzer (note that zeff
c                also appears in defn of xnue,xnuse)
c ------------------------------------------------------------------ HSJ
c
      implicit none
c
      real*8 xft,zeffa,fftrap,etaspitzer,etaparallel,xnue,xnuse,
     .       xnusem,enea,charge,xmasse,z,cr,zeta,zr,ftad,zta,cra,
     .       wneot,xkap11
c
      cr(z)   = 0.56*(3.0-z)/(3.0+z)/z
      zeta(z) = 0.58+0.20*z
c
c --- eq. 7.36, hirshman N.F. 1981, pg 1157
c
      zr(z) = ((0.222*z+1.198)*z+1.0)/((0.753*z+2.966)*z+1.0)
c
c     note: tauee = zeffa/xnue
c
      xnusem     = xnuse/zeffa
      xkap11     = 1.0 / zr(zeffa)
      etaspitzer = xmasse*xnue/(enea*charge**2)         ! in sec
      etaspitzer = etaspitzer/xkap11 ! eq 7.36, pg 1157, NF 1981
      zta        = zeta(zeffa)
      if (zeffa .lt. 3.0) then
        cra = cr(zeffa)
      else
        cra = 0.0
      end if
      ftad        = xft/(1.0 + zta*xnusem)
      fftrap      = wneot*ftad*(1.0 + cra*(1.0-ftad))
      etaparallel = etaspitzer/(1.0-fftrap)
      if(etaparallel .lt. 0.0)then
         print *,'sub neoresist, fftrap , etaspitzer =',
     .         fftrap,etaspitzer
         call STOP ('subroutine neoresist: failed resistivity', 0)
      endif
      return
c
      end

      subroutine radius_calc
c
c ----------------------------------------------------------------------
c     subroutine calculates major and minor radii as described below
c ----------------------------------------------------------------------
c
c  INPUT  ( all input is from common blocks):
c  geom.i
c  codeid               "onedee" use simplified 1d geometry
c                       otherwise use 2d mhd calculated  geometry
c  nj                # of grid points in r grid
c
c  from flxav.i:
c      npsi              # of values in psival,etc.
c
c  from IO.I:
c      ncrt              fortran unit # associated with terminal output
c
c  from mhdpar.i:
c      kpsi              #declared size of arrays used below
c
c  from numbrs.i:
c      nj                # grid point in rho array (ie r)
c  from PSIG.I:
c      psival(k)         k=1,2..npsi the mhd psi grid
c      rminavnpsi(k)     k=1,2..npsi, the average minor radius
c      rmajavnpsi(k)     k=1,2..npsi, the major radius of each flux surface
c
c  from RHOG.I:
c      rmajorvec(j)      j=1,..nj
c  from storage.i:
c      work vectors used below in equivalence statements
c
c  OUTPUT ( to rhog.i))
c     rmajor_r(j)    j=1..nj  rmajor in cm, where rmajor(j) corresponds tp
c                    r(j) (r(j) is the  rho grid,in cm,  used in the
c                    transport calcs) rmajor_r is the same as
c                    rmajorvec(j) calculated in getrmaj1.
c     rminor_r(j)    j=1,..nj is the average minor radius of the flux
c                    surface that corresponds to r(j). It is obtained by
c                    interpolation from rminavnpsi
c     ravg_r(j)      j=1,..nj is the major radius of each flux surface
c                    in cm. It is obtained by interpolation from
c                    rmajavnpsi
c
c
c----------------------------------------------------------HSJ-3-3-99
c
      USE param
      USE glf23
      USE io 
      USE mhdpar
      USE numbrs
      USE mesh
      USE machin
      USE geom
      USE psig
      USE rhog
      USE flxav
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c      include 'param.i'
c      include 'flxav.i'
c      include 'geom.i'
c      include 'machin.i'
c      include 'mesh.i'
c      include 'io.i'
c      include 'mhdpar.i'
c      include 'numbrs.i'
c      include 'psig.i'
c      include 'rhog.i'
      include 'storage.i'
c      include 'glf23.i'
c

      real*8
     .           aspline(kpsi), bspline(kpsi), cspline(kpsi),
     .           cs2spline(kpsi,3), dspline(kpsi), espline(kpsi),
     .           fspline(kpsi), bpar(4),psi_local(kpsi), psi_jk(kpsi)
      equivalence (aspline(1)    , xdum(1))
      equivalence (bspline(1)    , xdum(kpsi+1))
      equivalence (cspline(1)    , wdum(1))
      equivalence (dspline(1)    , vdum(1))
      equivalence (espline(1)    , ydum(1))
      equivalence (fspline(1)    , ydum(kpsi+1))
      equivalence (cs2spline(1,1), sdum(1))
c


      if (codeid .eq. 'onedee')then
         do j=1,nj
            rmajor_r(j) = rmajor
            rminor_r(j) = r(j)
            ravg_r(j)   = rmajor
         enddo
         
      else  !2d mhd geometry   
      tension = 0.0               ! don't use tension option of tspline
      tmax    = 0.0               ! max allowed tension
      bpar(1) = 0.0               ! set boundary conditions on spline
      bpar(2) = 0.0
      bpar(3) = 0.0
      bpar(4) = 0.0
      ier     = 0.0
c
          dpsi_local=(psir(nj)-psir(1))/(npsi-1)
          psi_local(1)=psir(1)
          do i=2,npsi-1
               psi_local(i)=psi_local(i-1)+dpsi_local
          enddo
          psi_local(npsi)=psir(nj)
cjek
          psi_jk(1)=psir(nj)
          do i=2,npsi-1
               psi_jk(i)=psi_jk(i-1)-dpsi_local
          enddo
          psi_jk(npsi)=psir(1)
c
          call tspline (psi_jk,rminavnpsi,npsi,bpar,cs2spline,kpsi,
     .                  ier,tension,aspline,bspline,cspline,dspline,
     .                  espline,fspline,tmax,psir,rminor_r,nj)
      if(ier .ne. 0 .and. ier .ne. 131)
     .  call STOP ('subroutine radius_calc: failed rminor_r', 0)

          call tspline (psi_jk,rmajavnpsi,npsi,bpar,cs2spline,kpsi,
     .                  ier,tension,aspline,bspline,cspline,dspline,
     .                  espline,fspline,tmax,psir,ravg_r,nj)
                call tspline (psi_jk,elongx,npsi,bpar,cs2spline,kpsi,
     .                  ier,tension,aspline,bspline,cspline,dspline,
     .                  espline,fspline,tmax,psir,kappa_nj,nj)
              
      if(ier .ne. 0 .and. ier .ne. 131)
     .  call STOP ('subroutine radius_calc: failed ravg_r', 0)
c      load rmajorvec:
c      call getrmaj1 (ncrt, nj)
c     copy rmajorvec to laocl array:
      call copya(rmajorvec,rmajor_r,nj)
c
c      open(unit=13,file='rcalc.out',status='unknown')
c      do i=1,npsi-1
c        write(13,*) i,psi_jk(i), psi_local(i)
c      enddo


      endif  !1d,2d geometry branches
      return
      end

      subroutine sauter99(nein,niin,tein,tiin,zz,zpi,ft,qq,eps,rr,
     .                    ell31,ell32,ell34,aone,sigma,jhirsh)
c-----------------------------------------------------------------------
c     copied from nvloop.pro (P.Politzers code)    HSJ 1/9/01
c
c  sauter99.pro
c
c  4/22/99 p a politzer
c  5/11/99 repaired error in jboot formula
c  modified for use in ONETWO HSJ
c
c  bootstrap current, collisionality, conductivity
c      according to Sauter, et al.
c      O. Sauter, C. Angioni, and Y. R. Lin-Liu
c      CRPP/EPFL report LRP 630/99, February 1999
c  jboot = <Jbs.B>/B_T0
c  this paper includes collisions, geometry
c  this paper applies to a single ion speciesc
c     multiple ions are incorporated here by replacing z with zeff
c
c  inputs:
c    nein         = electron density ( cm^-3)
c    niin         = total ion density
c    tein         = electron temperature (keV)
c    tiin         = ion temperature (keV)
c    zz           = average charge
c    zpi          = prim ion charge
c    ft           = trapped particle fraction
c    qq           = safety factor
c    eps          = inverse aspect ratio
c    rr           = major radius to center of surface
c    jhirsh       =112 is used to set z to main ion charge in nu*i and lnlni
c                  if jhirsh ne 112 then zeff is used in these two formulas
c
c----------------------------------------------------------------------
       implicit none

       real *8
     . te,tein,ti,tiin,nee,nein,nii,niin,zz,
     . nustare,nustari,qq,rr,lnlame,lnlamii,eps,
     . aqnue,sqnui,sigsp,azero,aone,ft,
     . f31,f32ee,f32ei,f33,f34,zp1,a1,a2,a3,
     . b1,b2,b3,ff32ee,ff32ei,ell31,ell34,ell32,sigma,
     . sqnue,zpi,zzm
       integer jhirsh
c      renormalize temperature and density to eV and #/m^3
c      (these units are used in the original report)
       te = tein*1000.
       ti = tiin*1000.
       nee = nein*1.e6
       nii = niin*1.e6

c      ion density (use ni = n_fuel + n_impurity)
c      zimpurity = 6.0
c      nii = nee*(1.-(zz-1.)/zimpurity)
c      pratio = nii*ti/nee/te


       lnlame = 31.3-log(sqrt(nee)/te)

c       eps = eps > 1.e-10
       nustare = 6.921e-18*abs(qq)*rr*nee*zz*lnlame/te**2/eps**1.5


       IF(jhirsh  == 112)THEN
          zzm = zpi
       ELSE
          zzm = zz
       ENDIF
       lnlamii = 30.-log(zzm**3*sqrt(nii)/ti**1.5)
       nustari = 4.90e-18*abs(qq)*rr*nii*zzm**4*lnlamii/ti**2/eps**1.5



       sqnue = sqrt(nustare)
       sqnui = sqrt(nustari)

       sigsp = 1.9012e4*te**1.5/zz/lnlame/(0.58+0.74/(0.76+zz))

        azero = -1.17*(1.-ft)/(1.-0.22*ft-0.19*ft**2)
c       aone = ((azero+0.25*(1.-ft**2)*sqnui)/(1.+0.5*sqnui) 
c     .     - 0.315*nustari**2*ft**6)/(1.+0.15*nustari**2*ft**6)
c      changed above to below due to errata,dec. 2002, Physic of Plasmas HSJ
       aone = ((azero+0.25*(1.-ft**2)*sqnui)/(1.+0.5*sqnui) 
     .     + 0.315*nustari**2*ft**6)/(1.+0.15*nustari**2*ft**6)

       f31 = ft/(1.+(1.-0.1*ft)*sqnue+0.5*(1.-ft)*nustare/zz)
       f32ee = ft/(1.+0.26*(1.-ft)*sqnue+0.18*(1.-0.37*ft)
     .                                        *nustare/sqrt(zz))
       f32ei = ft/(1.+(1.+0.6*ft)*sqnue+0.85*(1.-0.37*ft)*
     .                                           nustare*(1.+zz))
       f33 = ft/(1.+(0.55-0.1*ft)*sqnue+0.45*(1.-ft)*nustare/zz**1.5)
       f34 = ft/(1.+(1.-0.1*ft)*sqnue+0.5*(1.-0.5*ft)*nustare/zz)

       zp1 = zz+1.
       ell31 = (1.+1.4/zp1)*f31-1.9/zp1*f31**2+0.3/zp1*f31**3
     .                                       +0.2/zp1*f31**4
       a1 = (0.05+0.62*zz)/zz/(1.+0.44*zz)
       a2 = 1./(1.+0.22*zz)
       a3 = 1.2/(1.+0.5*zz)
       b1 = (0.56+1.93*zz)/zz/(1.+0.44*zz)
       b2 = 4.95/(1.+2.48*zz)
       b3 = 1.2/(1.+0.5*zz)
       ff32ee = a1*f32ee*(1.-f32ee**3)+a2*f32ee**2*(1.-1.2*f32ee+
     .               0.2*f32ee**2) +a3*f32ee**4
       ff32ei = -b1*f32ei*(1.-f32ei**3)+b2*f32ei**2*(1.-0.55*f32ei
     .            -0.45*f32ei**2) -b3*f32ei**4
       ell32 = ff32ee+ff32ei
c       ell34 = (1.+1.4/zp1)*f34-1.9/zp1*f34**2+0.3/zp1*f34**3+
c     .                                             0.2/zp1*f34**4
      ell34 = (((0.2/zp1*f34+0.3/zp1)*f34-1.9/zp1)*f34+(1.+1.4/zp1))*f34
c       sigma = sigsp*(1.-(1.+0.36/zz)*f33+0.59/zz*f33**2-0.23/zz*f33**3)
      sigma = sigsp*(1.+((-0.23/zz*f33 +0.59/zz)*f33-(1.+0.36/zz))*f33) 
c      ignore gradient of zeff
c      (else add -ell31*ti/zimpurity/te*d(zz)/d(rhohat))
c      jb0 = 1.6022e-19*nee*te*qq*factor
c      jbootc = -jb0*(ell31*(1.+pratio)*neph 
c    .                 +(ell31+ell32)*teph 
c     .                 +pratio*(ell31+aone*ell34)*tiph)

c  return,{nustare:nustare, sigma:sigma, jbootc:jbootc}
      return
      end


        SUBROUTINE  rot_terms
c  -------------------------------------------------------------------
c       Determine wdnidt and aniwdwdt
c  ---------------------------------------------------------HSJ-4/7/14
        USE param,               ONLY : kj
        USE solcon,              ONLY : dtsumi,dtt
        USE io,                  ONLY : itimav
        USE numbrs,              ONLY : nprim,nion,nk,nj
        USE soln,                ONLY : u,usave,uav,uav0
        USE ions,                ONLY : atw
        USE constnts,            ONLY : xmassp,kevperg
        USE geom,                ONLY : r2capi
        USE tordlrot,            ONLY : angrcple,wdnidt,aniwdwdt
        IMPLICIT none
        REAL*8 dtavg,sum,avg1,avg0,dndt,angterm,dnt
        REAL*8 amassden(kj) ! move this to ions module
        INTEGER jin,k,i

        DO jin=1,nj
           amassden(jin) = 0.0D0
           DO k =1,nion
              amassden(jin) = amassden(jin) +atw(k)*u(k,jin)
           ENDDO
           amassden(jin) = amassden(jin) * xmassp 

           dtavg = 0.0D0
           if (   dtt .ne. 0.0D0)  dtavg = 1.0D0 / dtt
           if (itimav .eq. 1  )  dtavg = dtsumi
           sum = 0.0D0

           DO i=1,nion
              avg1 = u(i,jin)
              avg0 = usave(i,jin)
              if (itimav .eq. 1)  avg1 = uav(i,jin)
              if (itimav .eq. 1)  avg0 = uav0(i,jin)
              IF (i .LE. nprim)THEN
                 dndt = (avg1-avg0)*dtavg*atw(i)
                 sum  = sum + dndt
              ENDIF
           ENDDO

           avg1 = u(nk,jin)
           avg0 = usave(nk,jin)
           if (itimav .eq. 1)  avg1 = uav (nk,jin)
           if (itimav .eq. 1)  avg0 = uav0(nk,jin)
           dnt         = (avg1-avg0)*dtavg
           sum         = sum * 0.5D0 * r2capi(jin)*avg1**2*xmassp
           angterm     = amassden(jin)*r2capi(jin)*avg1*dnt
           wdnidt(jin)   = sum*angrcple*kevperg
           aniwdwdt(jin) = angterm*angrcple*kevperg

        ENDDO

        RETURN
        END


        SUBROUTINE visc_flux
c  -------------------------------------------------------------------
c       Determine energy fluxes associated with rotation
c  ---------------------------------------------------------HSJ-5/5/14
        USE bd_condtn,           ONLY : fluxb
        USE constnts,            ONLY : kevperg,jouperg
        USE flx,                 ONLY : flux
        USE geom,                ONLY : hcap
        USE mesh,                ONLY : r,ra,drr
        USE numbrs,              ONLY : nprim,nion,nk,nj
        USE solcon,              ONLY : theta
        USE soln,                ONLY : u,usave,uav,uav0
        USE tcoef,               ONLY : dudr
        USE tordlrot,            ONLY : angrcple,omegapi,fluxangv,
     .                                  vischeat,flxangce,qangce,
     .                                  qomegapi,fluxangc

        IMPLICIT none
        REAL*8 angrota,omegapib,onemt
        INTEGER j

          onemt = 1.D0 - theta
          if (angrcple .ne. 0.0) then
            DO j=1,nj-1
                angrota    = (theta*(u(nk,j)+u(nk,j+1))+onemt *
     .                     (usave(nk,j) + usave(nk,j+1))) * 0.5
                omegapi(j) = angrcple*kevperg*angrota*fluxangv(j)
c --- viscous heating term in w/cm**3:
                vischeat(j) = -fluxangv(j)*jouperg*angrcple*dudr(nk,j)
                flxangce(j) = 0.5*angrcple*angrota*kevperg*fluxangc(j)
c                flux(nion+2,j) = flux(nion+2,j)+omegapi(j)+flxangce(j)

            ENDDO 

              call extrap (ra(nj-2),ra(nj-1),r(nj),flux(nion+2,nj-2),
     .                     flux(nion+2,nj-1),fluxb(nion+2))
              call extrap (ra(nj-2),ra(nj-1),r(nj),
     .                     flxangce(j-2),flxangce(j-1),flxangce(nj))
              call extrap (ra(nj-2),ra(nj-1),r(nj),omegapi(nj-2),
     .                     omegapi(nj-1),omegapib)
              omegapi(nj) = omegapib

              call divflx (omegapi,omegapib,hcap,r,ra,drr,1,nj,1,
     .                     qomegapi)
              call divflx (flxangce,flxangce(nj),hcap,r,ra,drr,1,nj,1,
     .                     qangce)
           end if

        RETURN
        END SUBROUTINE visc_flux


      subroutine sauter  (jhirsh, j, nion, fcapa, enea, tea, ena, tia,
     .                    zeffa,q,xft,epsa,rrmajor,
     .                    cee, eta, ra, rmajor, enbeam,
     .                    enalp, sumdz, za, wneo,btor,
     .                    d, dfast, dfion, xden, xdit, xdet, xdin)
c
c ---------------------------------------------------------HSJ----------
c  this subroutine sets up the diffusion coefficients so that the
c  bootstrap current model given by Sauter (CRPP/EPFL report LRP 630/99,
c   February 1999 ) will be used.

c
c    in ONETWO the bootstrap is calculated as
c    <jboot dot b/bto>= -s*(sum over ions {d(kfar,i,j)*dni/drho )
c                           +d(kfar,nion+1,j)*dte/drho
c                           +d(kfar,nion+2,j)*dti/drho
c                           +dcoef(j)*dnf/drho  )
c    where
c    s = rho/(cee*eta)  (a pure number)
c    sum over ions is sum from i = 1 to i=nion,and includes
c    primary as well as impurity species.
c    kfar = nion+3, indexes Faraday's law.
c    dcoef is a special coefficient used to model electron
c    density due to fast ions and fast ion density (on option)
c    it is calculated in subroutine SOURCE as dcoef = dfast+dfion
c    dfion,the contribution from fast ions, is not implemented .
c
c --- input through argument list:
c
c  jhirsh                selects bootstrap current model.
c                        jhirsh= 110  use Sauter bootstrap but
c                        dont use resistivity predicted by this model
c                        jhirsh = 111  means return the resistivity
c                        eta(j) as determined by the Sauter model.
c                        jhirsh = 112 passes z of primary ion into
c                        sauter99 routine directly ( 111 option is not available
c                        in this case)
c
c  j                   is the grid index. note that we are working
c                      on the half grid here so j ranges from 1 to nj-1
c
c  nion                number of ions (primary plus impurity)
c
c  fcapa               f(psilim)/f(psi) on the (half) grid point j
c
c  enea                electron density, #/cm**3, at j
c
c  tea                 electron temp, keV, at j
c
c  ena(i)              i = 1,2..nion,density,#/cm**3 at j
c
c  tia                 common ion temp, keV, at j
c
c  zeffa               zeff at j
c
c  q                   ( signed) safety factor at j
c
c  xft                 trapped particle fraction
c
c  epsa                inverse aspect ratio at j (ie the half grid)
c
c  rrmajor             major radius of flux surface j,cm
c
c  cee                 speed of light,cm/sec
c
c  eta(j)              resistivity,sec,on half grid
c
c  ra(j)               rho,cm,on half grid
c
c  enbeam(j)           beam density,on full grid,#/cm**3
c
c  enalp(j)            alpha density,#/cm**3,on full grid
c
c  sumdz               sum of (ni*d<z>/dte)*te/ne at j
c
c  za(i)               i = 1,2,..nion z of ion i at j
c                      
c  wneo                neoclassical multiplier array
c  btor                toroidal b field at rmajor, in vaccuum, gauss
c  rmajor              major radius where btor is quoted,cm
c
c --- input from include files
c    param.i:
c       kj               size of rho grid
c

c --- output:
c  through argument list:
c
c diffusion coefficients are defined so that each one makes
c a contribution of the form -s*d(kfar,i,j) to <jboot dot b/bto>
c see above for definition of s
c  d(kfar,i,j)          i = 1,2,..nion ,statamp.cm**2
c                       gives density contribution at half grid pt j
c  d(kfar,i,j)          i = nion+1 diffusion coefficient,statamp/(cm keV)
c                       gives the contribution at half grid pt j
c  d(kfar,i,j)          i = nion+2 diffusion coefficient,statamp/(cm keV)
c                       gives ti contribution at half grid point j
c
c  dfast(j)             electron "diffusion" coefficient for electrons
c                       sourced by the fast ions. It is assumed that zbeam = 1
c                       dfast has units of statamp cm**2
c                       the contribution to <jboot dot b/bto> is
c                                -dfast*(dnf/drho)/xres
c                       A SOURCE TERM
c                            (1.0/(HCAPA*RHO)*(D/DRHO)FA
c                       (WHERE
c                            FA = (HCAPA*RHO*DFAST*ZBEAM*(DNF/DRHO)/XRES
c                        AND DNF/DRHO IS THE FAST ION DENSITY GRADIENT )
c                        IS ADDED TO THE RHS OF FARADAY'S LAW
c
c  dfion(j)              fast ion "diffusion" coefficient.
c                        dfion has units of statamps cm**2
c                        dfion is set to zero for now.
c                        (the bootstrap from the fast ions is not
c                         represented correctly by the current theory
c                         and it is best to just ignore it.
c                        to include it:
c                        the contribution to <jboot dot b/bto> is
c                            -dfion*(dnf/drho)/xres
c                       A SOURCE TERM
c                            (1.0/(HCAPA*RHO)*(D/DRHO)FA
c                       (WHERE
c                            FA = (HCAPA*RHO*DFION(DNF/DRHO)/XRES
c                        AND DNF/DRHO IS THE FAST ION DENSITY GRADIENT )
c                        IS ADDED TO THE RHS OF FARADAY'S LAW
c                        also if dfion is used (because jhirsh = 89)
c                        then the total ion density includes the fast ions
c
c  ancillary arrays:
c
c  xden(j)                 is defined so that ion species i contributes
c                               -xden*zi*(dni/drho)/xres
c                          to <jboot dot b/bto>
c                          (accounts for electron density gradients)
c
c  xdin(i,j)               is defined so that ion species i contributes
c                               -xdin*(dni/drho)/xres
c                          to jboot dot b/bto>
c                          (accounts for ion density gradients)
c
c  xdet(j)                 is the same as d(kfar,nion+1,j)
c
c  xdit(i,j)               is defined so that ion species i contributes
c                               -xdit*(dti/drho)/xres
c                          to <jboot dot b/bto>
c

c ------------------------------------------------------------------ HSJ
c
      USE param
      implicit none
c
      integer j,kfar,l,nion,jhirsh

c
c      include 'param.i'
      include 'sauter.i' !sl31,sl32

c
      integer k
      real *8    
     .          cee,nia,enf,enea,fcapa,sumdz,tea,
     .          tia,effz,xres,
     .          q(*),xft,epsa,rrmajor,rmajor,zeffa,
     .          const31,const32,const34,temp,
     .          eta(*),ra(*),ena(*),enbeam(*),enalp(*),
     .          d(kk,kk,*),xden(*),xdet(*),xdin(kion,*),
     .          xdit(kion,*),wneo(5,*),za(*),dfast(*),dfion(*),
     .          xunits,pea,beta,ell31,ell32,ell34,aone,sigma,
     .          btor,electron_part,ion_part,safac,rrmeters,zpi

      kfar =nion+3   !index to poloidal B field evolution equation

      zpi = za(1)    ! set z of primary ion for use in sauter99 if jhirsh = 112

c
c --- get the total thermal density
c
      nia = 0.0
      do l=1,nion
        nia = nia + ena(l)
      end do
c
c --- add fast ion density if called for
c
c      if (jhirsh .eq. 111) then
c        enf = (enbeam(j+1) + enbeam(j)) * 0.5
c        enf = enf + enalp(j+1) + enalp(j)    ! z = 2 for alpha particles
c        nia = nia + enf
c      end if
      if(jhirsh .eq. 111)then
c         effz  = enea / nia  ! to account for multiple ion species (approx)
          effz  = zeffa       !changed 7/24/01 HSJ
      else
         effz  = zeffa
      endif
c
c
c  we have 1.602e-8 amps/cm**2 = 1 keV/(cm**4 gauss)
c          1 amp = 2.998e9 statamps
c  multiplication by xunits changes keV/(cm**4 gauss) to statamps/cm**2
c
      xunits = 1.602e-8 * 2.998e9
      pea =enea * tea    ! kev/cm**3
c     note q and btor are signed quantities. 
      safac = 0.5*(q(j+1)+q(j))
      beta = -(safac/fcapa)*(rmajor/ra(j))*xunits*pea/btor   ! statamps/cm 
      rrmeters = rrmajor*0.01
 
      !jhirsh = 112 ,uses z main ion, internal to sauter99
      call sauter99(enea,nia,tea,tia,effz,zpi,xft,safac,epsa,rrmeters,
     .              ell31,ell32,ell34,aone,sigma,jhirsh)

       sl31(j) = ell31
       sl32(j) = ell32
       sl34(j) = ell34
c     ellij output by sauter99 are dimensionless quantities
c     can use Sauter sigma in following (sigma is in ohm^-1*m^-1
      if(jhirsh .eq. 111)eta(j) = 1./sigma/8.98755179e9  !sec
      xres = cee * eta(j) / ra(j)    ! pure number
      const31 = beta*ell31 ! statamps/cm 
      const32 = beta*ell32 ! statamps/cm 
      const34 = beta*ell34*aone ! statamps/cm 
      temp = -wneo(4,1)*const31*xres/pea ! statamps *cm^2/kev
      xden(j) = temp*tea                 ! statamps *cm^2
      !note za(k) is added to xden  in bootstrap contribution in source 
      dfion(j) =0.0                        !no contribution from fast ions
c     get contribution for fast electrons ( the contribution for electrons
c     grouped with thermal ions was done above):
c      dfast(j) = -wneo(4,1)*const31*xres*tea/pea      !statamps cm**2 
      dfast(j) = xden(j)                               !statamps cm**2

      do k = 1,nion
c        these d have units of statamp.cm**2
         electron_part = tea*za(k)
         ion_part = tia
         d(kfar,k,j) =  temp*(electron_part+ion_part)
         xdin(k,j) = temp*ion_part
      enddo

c     d has units of statamp/(cm keV) for te,ti terms:
c        Te term:
      d(kfar,nion+1,j) = -wneo(4,2) *xres *(const31 + const32)/tea
         xdet(j) = d(kfar,nion+1,j)
c        T term:
      d(kfar,nion+2,j) = -wneo(4,3)* xres *nia*(const31 +const34 )/pea
c --- ratio out the total so that each ion species receives
c --- the fraction of the total density as its bootstrap contribution.
c --- note that this is arbitrary and is done only for printout purposes
      do l=1,nion
        xdit(l,j) = d(kfar,nion+2,j)*ena(l)/nia
      end do
c     contribution due to dz/drho gradient is treated as a source term:
c     dz/drho = dz/dte * dte/drho 
c      xdiz(j) =  0.0 !not calculated for now (use sumdz to do it)


      return
c
      end

