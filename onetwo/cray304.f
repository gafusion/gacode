      subroutine nclass_dr (wneo, jneo, jhirsh,
     .                      curbs, curbs_bt, curbs_bp, eta)
c
c ----------------------------------------------------------------------
c  NCLASS neoclassical transport model
c  Implemented by Jon Kinsey, Oct. 1999
c  Minor format changes by Joe Freeman, 16 February 2000
c
c  This subroutine computes the bootstrap current, parallel
c  resistivity and neoclassical fluxes for a multi-species plasma
c  of arbitrary aspect ratio, geometry, and collisionality using
c  the NCLASS model (W. Houlberg, K. Shaing, S. Hirshman, M. Zarnstorff,
c  Phys. Plasmas 4, 3230 (1997).
c
c  The total bootstrap current density <J.B>_bs is given as the
c  flux surface average, p_bsjb (A-T/m2). It is transformed into
c  the form <J_bs dot B/Bt0> (A/cm2) and stored as curbs(j).
c  The bootstrap current is also computed as a set of coefficients
c  of the pressure and temperature gradients of each species,
c  bsjbt_s and bsjbp_s. This is then transformed into the following
c  form for use in ONETWO:
c
c    <jboot dot b/bto>= -s*(sum over ions {d(kfar,i,j)*dni/drho )
c                           +d(kfar,nion+1,j)*dte/drho
c                           +d(kfar,nion+2,j)*dti/drho
c                           +dcoef(j)*dnf/drho  )
c    where s = rho/(cee*eta)  (a pure number).
c    sum over ions is sum from i = 1 to i=nion,and includes
c    primary as well as impurity species.
c
c  Note: Currently, potato orbit effects have been disabled
c  (k_potato=0).
c
c --- input through argument list:
c
c  wneo                  transport coefficients
c
c  jneo                  switch for neoclassical transport
c                             0 = do not account for transport
c                             1 = add into total transport
c
c  jhirsh                switch for bootstrap current model
c                        elsewhere in Onetwo. Here, it is a
c                        switch for fast ion contributions:
c                             99 = no fast ions
c                            100 = include fast ions
c
c --- output:
c
c  curbs(j)              total bootstrap current density
c                        <J_bs dot B/Bt0> at zone j (A-T/cm2)
c
c  curbs_bt(j)           bootstrap current density driven
c                        by temperature gradients
c                        <J_bs dot B/Bt0> at zone j (A-T/cm2)
c
c  curbs_bp(j)           bootstrap current density driven
c                        by pressure gradients
c                        <J_bs dot B/Bt0> at zone j (A-T/cm2)
c
c  eta(j)                parallel resistivity (s)
c
c ----------------------------------------------------------------------
c
c ONETWO include files used to pass variables through common blocks
c
c
      USE param
      USE fusion
      USE ions
      USE soln
      USE mhdpar
      USE nub2
      USE extra
      USE numbrs
      USE mesh
      USE machin
      USE geom
      USE tordlrot
      USE tcoef
      USE metrics
      USE neo2d
      USE cer
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c
c Common blocks for transport, d(i,j,k), etc.
c
c      include 'tcoef.i'
c
c More terms for the beam and alpha components, Added 8-21-95 by DFF
c
c      include 'nub2.i'
c      include 'fusion.i'
c
c  CER velocities and magnetic components
c
c      include 'cer.i'
c      include 'tordlrot.i'
c
c  Metrics calculated in ONETWO to be used in calculation
c
c      include 'metrics.i'
c
c  variables to be passed to NCLASS
c
      include 'houlberg.i'
c     include 'pamx_mi.inc'
c     include 'pamx_mz.inc'
c     include 'pamx_ms.inc'
c
c  Variables used locally
c
      integer ima, k, j, itersqz
      integer k_order, k_potato, m_i, m_z, im, iza, iz, idum(8)
      integer iflag, lprintcl_test, kfar
      real*8  c_den, c_potb, c_potl, p_b2, p_bm2, p_eb
      real*8  p_grbm2, p_ngrth, p_fhat, p_ft, p_grphi, p_gr2phi
      real*8  p_fm(3)
      real*8  curbs(kj), curbs_bt(kj), curbs_bp(kj), xmks, dvedr
      real*8  fcapa, rbp1, rbp2, rbpa, dpsidr
      real*8  rbta, udia_da(kj), Kpol_da(kj)
      real*8  udia_ca(kj), Kpol_ca(kj)
      real*8  sqzmax, sqzmax_old
      real*8  Epsi_expa(kj), Epsia(kj)
      real*8  angrot_ca(kj), angrot_da(kj)
      real*8  sqzmin, sqz_da(kj), ugrta(kj)
      real*8  xjdotb_bt, xjdotb_bp, sum(8)
      real*8  chiieff_ncl(kj), vieff_ncl(kj),
     .   chieeff_ncl(kj), veeff_ncl(kj), ni_ncl(kj)
      real*8  dq_s(mx_ms), vq_s(mx_ms)
      real*8  wneo(5,*), xres, sumdz, enea, tea, tia, nia, effz
      real*8  eta(kj), xdne(kj), xdni(kj), xdte(kj), xdti(kj)
      real*8  ena(kion), dzdtea(kion)
c
c  Variables output from NCLASS
c
      data cee    / 2.99792458e+10 /    ! speed of light  (cm/s)
      integer m_s
      integer jm_s(mx_ms), jz_s(mx_ms)
      real*8 p_bsjb, p_etap, p_exjb
      real*8 calm_i(3,3,mx_mi),
     .   caln_ii(3,3,mx_mi,mx_mi), capm_ii(3,3,mx_mi,mx_mi),
     .   capn_ii(3,3,mx_mi,mx_mi),
     .   bsjbp_s(mx_ms), bsjbt_s(mx_ms),
     .   dn_s(mx_ms), gfl_s(5,mx_ms),
     .   qfl_s(5,mx_ms), sqz_s(mx_ms),
     .   upar_s(3,3,mx_ms), utheta_s(3,3,mx_ms),
     .   vn_s(mx_ms), veb_s(mx_ms),
     .   qeb_s(mx_ms), xi_s(mx_ms), ymu_s(3,3,mx_ms),
     .   chip_ss(mx_ms,mx_ms), chit_ss(mx_ms,mx_ms),
     .   dp_ss(mx_ms,mx_ms), dt_ss(mx_ms,mx_ms)
c
      save    init
      data    init /0/
c
      iflag = 0
      lprintcl_test = 0
      mimz4     = 4 * mxmi * mxmz
      xmks      = 1.0e6
      inclprint = 0
      kncord    = 2
      sqzmin    = 1.0e-6
      zero      = 0.D0
c
c  Initialize output variables
c
      m_s = 0
      do j = 1, mx_ms
        jm_s(j) = zero
        jz_s(j) = zero
        bsjbp_s(j) = zero
        bsjbt_s(j) = zero
        dn_s(j) = zero
        sqz_s(j) = zero
        vn_s(j) = zero
        veb_s(j) = zero
        qeb_s(j) = zero
        xi_s(j) = zero
c
        do k = 1, 3
          do l = 1, 3
            upar_s(k,l,mx_ms) = zero
            utheta_s(k,l,mx_ms) = zero
            ymu_s(k,l,mx_ms) = zero
          enddo
        enddo
c
        do k = 1, 5
          gfl_s(k,j) = zero
          qfl_s(k,j) = zero
        enddo
c
        do k = 1, mx_ms
          chip_ss(j,k) = zero
          chit_ss(j,k) = zero
          dp_ss(j,k) = zero
          dt_ss(j,k) = zero
        enddo
      enddo
c
      do j = 1, mx_mi
        do k = 1, mx_mi
          do l= 1, 3
            do m = 1, 3
               calm_i(m,l,k) = zero
               caln_ii(m,l,k,j) = zero
               capm_ii(m,l,k,j) = zero
               capn_ii(m,l,k,j) = zero
            enddo
          enddo
        enddo
      enddo
c
      p_bsjb = zero
      p_etap = zero
      p_exjb = zero
c
      if (init .eq. 0) then
c
c ***** Initialize fixed parameters.
c
c  Set constants
c
        call inicon
c
c  Set cutoff density
c
        dencut = 1.0e10
c
c  Set maximum charge state mz = 10 (up to Neon)
c
        mz = 10 ! should be set to mxmz HSJ
        mz = mxmz ! HSj 12/16/98
c
c  Set masses and charges. mith is number of THERMAL species!
c
        mith     = 1 + nion
        amuai(1) = elemas / promas
        do k = 2, mith
          amuai(k) = atw(k-1)
        end do
        init = 1
c
      end if
c
      if (inclprint .eq. 1) then
        write (6, *) '================================================'
        write (6, '(1x, 14(a10, 2x))')
     .               'Rho', 'xngrth', 'xfmsum', 'xb2', 'xbm2', 'xfhat',
     .               'xft', 'xgrbm2', 'denai(2,1)', 'tai(2)',
     .               'xgrp(2,1)','xgrt(2)', 'curbs'
      end if
c
c Determine if fast ions are to be included.
c Set IFASTION = 1 if JHIRSH = 100
c
      if (jhirsh .eq. 100)  ifastion = 1
c
c ----------------------------------------------------------------------
c Start big loop over centered-mesh
c ----------------------------------------------------------------------
c
      itersqz    = 0
      sqzmax     = 1.0
      sqzmax_old = 1.0
      call sscal (mi  , zero, tai      , 1)
      mimz       = mxmi * mz
      call sscal (mimz, zero, denai    , 1)
      call sscal (mimz, zero, fex      , 1)
      call sscal (mimz, zero, xzi      , 1)
      call sscal (mi  , zero, xgrt     , 1)
      call sscal (mimz, zero, xgrp     , 1)
      call sscal (kj  , zero, Epsi     , 1)
      call sscal (kj  , zero, Epsia    , 1)
      call sscal (kj  , zero, Epsi_expa, 1)
      call sscal (kj  , zero, Kpol_da  , 1)
      call sscal (kj  , zero, Kpol_ca  , 1)
      call sscal (kj  , zero, angrot_da, 1)
      call sscal (kj  , zero, angrot_ca, 1)
      call sscal (kj  , zero, udia_ca  , 1)
      call sscal (kj  , zero, udia_da  , 1)
      call sscal (kj  , zero, ugrta    , 1)
c
  100 do j=1,nj-1
        drho = dr(j) / 100.0
c
c       Set metrics:
c
        xft    = 0.5 * (ftncl(j+1)+ftncl(j) )
        xb2    = 0.5 * (bsq(j+1) +bsq(j) )
        xbm2   = 0.5 * (bmsq(j+1) +bmsq(j) )
        xgrbm2 = 0.5 * (grbmsq(j+1) +grbmsq(j) )
        xngrth = 0.5 * (grth(j+1) + grth(j) )
        do m=1,3
          p_fm(m) = gfm(m,j)
        end do
c
c       Calculate xfhat - NOTE that units cancel:
c
        fcapa  = 0.5 * (fcap(j+1) +fcap(j))
        rbp1   = rbp(j)/(fcap(j)*gcap(j)*hcap(j))
        rbp2   = rbp(j+1)/(fcap(j+1)*gcap(j+1)*hcap(j+1))
        rbpa   = 0.5 * (rbp1 + rbp2)
        dpsidr = (rbpa / ra(j)) * rmajor
        xfhat  = btor * rmajor / (fcapa * dpsidr)
c
c      Calculate normalized Electric Field derivative at mesh center
c
        rbta  = btor*rmajor/(fcapa*xmks)
        dvedr = (promas*rbta*xfhat/(coulom*(btor*1.0e-4)**2))*
     .          (Epsi(j+1) - Epsi(j))/drho
c
c **** Calculate some dependent variables at mesh centers
c
c   First Electrons:
c      (Determine from ion densities)
c
        ene0       = zero
        ene1       = zero
        do k=1,nion
          ene0     = ene0 + en(j,k)  * z(j,k)
          ene1     = ene1 + en(j+1,k)* z(j+1,k)
        end do
        ene0       = ene0 + 1.0 * enbeam(j) + 2.0 * enalp(j)
        ene1       = ene1 + 1.0 * enbeam(j) + 2.0 * enalp(j)
        denai(1,1) = xmks*(ene1 + ene0) / 2.0
        fex(1,1,1) = zero
        fex(2,1,1) = zero
        fex(3,1,1) = zero
        tai  (1  ) = 0.5 * (te(j+1) + te (j))
        xgrp (1,1) = xmks*(ene1 * te(j+1) - ene0 * te(j))/drho
        xgrt (1  ) = (te(j+1) - te(j)) / drho
        xzi  (1,1) = 1.0
        sqz  (1,1) = 1.0 + amuai(1)*dvedr
        if (sqz(1,1) .lt. sqzmin)
     .  sqz  (1,1) = sqzmin
        sqz  (1,1) = (sqz(1,1))**1.5
        sqzmax     =  MAX (sqzmax, sqz(1,1))
c
c   Ions and Impurities (k=1,..NION):
c     For now we assume that only a single (fully stripped)
c     charge state exists for each isotope, xzi=1.0
c
        do k=1,nion
          ima           = k+1
          ia            = ANINT (z(j,k))
          denai(ima,ia) = 0.5 * xmks*(en(j+1,k)+en(j,k))
          fex(1,ima,ia) = zero
          fex(2,ima,ia) = zero
          fex(3,ima,ia) = zero
          tai(ima)      = 0.5 * (ti(j+1) + ti(j))
          xgrt(ima)     = (ti(j+1)-ti(j))/drho
          xgrp(ima,ia)  = xmks*(en(j+1,k)*ti(j+1)-en(j,k)*ti(j))/drho
          xzi(ima,ia)   = 1.0
          sqz(ima,ia)   = 1.0 - amuai(ima)*dvedr/z(j,k)
          if (sqz(ima,ia) .lt. sqzmin)
     .    sqz(ima,ia)   = sqzmin
          if (k .eq. 1)  sqz_da(j) = sqz(ima,ia)
          sqz(ima,ia)   = (sqz(ima,ia))**1.5
          sqzmax        =  MAX (sqzmax, sqz(ima,ia))
        end do
c
c   Set the number of ion species to the number of thermal species
c       (May be incremented if beam and alpha ions are present)
c
        mi = mith
c
c ------------ FAST ION CONTRIBUTIONS ----------------------------------
c   If fastion = 1 then Add Beam and Alpha particles
c   For now assume an effective temperature derived from the  averaged
c   stored energy density of the beam and alphas, ie:
c                    press = nkT = 0.666667*W
c   The (2/3) factor arises from Kinetic Theory.
c   We have to be careful since enbeam or enalp may be zero at some points
c   on the grid. Therefore, only include these ion densities in the mi-long
c   arrays if greater than dencut (Houlberg's routines will crash for
c   zero denisity species).
c
        if (ifastion .eq. 1) then
c
c         Beam Particles:
c
          enbeama = 0.5 * xmks*(enbeam(j+1)+enbeam(j))
          if (enbeama .ge. dencut) then
            if (enbeam(j) .ne. 0) then
              tbeamj0     = 0.666667 * wbeam(j)  /enbeam(j)
              tbsave      = MIN (tbeamj0, tbsave)
            else
              tbeamj0     = tbsave
            end if
            if (enbeam(j+1) .ne. 0) then
              tbeamj1     = 0.666667 * wbeam(j+1)/enbeam(j+1)
            else
              tbeamj1     = tbsave
            end if
            mi  = mi + 1
            ima = mi
            ia  = 1
            amuai(ima)    = 2.0
            denai(ima,ia) = enbeama
            fex(1,ima,ia) = zero
            fex(2,ima,ia) = zero
            fex(3,ima,ia) = zero
            tai  (ima)    = 0.5 * (tbeamj0 + tbeamj1)
            xgrt (ima)    = (tbeamj1 - tbeamj0)/drho
            xgrp (ima,ia) = xmks*0.666667*(wbeam(j+1)-wbeam(j))/drho
            xzi  (ima,ia) = 1.0
            sqz  (ima,ia) = 1.0 - amuai(ima)*dvedr
            if (sqz(ima,ia) .lt. sqzmin)
     .      sqz  (ima,ia) =  sqzmin
            sqz  (ima,ia) = (sqz(ima,ia))**1.5
            sqzmax        = MAX (sqzmax, sqz(ima,ia))
          end if
c
c         Alpha Particles:
c
          enalpa = 0.5 * xmks*(enalp(j+1)+enalp(j))
          if (enalpa .ge. dencut) then
            if (enalp(j) .ne. 0) then
              talpj0      = 0.666667 * walp(j)  /enalp(j)
              talpsave    = MIN (talpj0, talpsave)
            else
              talpj0      = talpsave
            end if
            if (enalp(j+1) .ne. 0) then
              talpj1      = 0.666667 * walp(j+1)/enalp(j+1)
            else
              talpj1      = talpsave
            end if
            mi            = mi + 1
            ima           = mi
            ia            = 2
            amuai(ima)    = 4.0
            denai(ima,ia) = enalpa
            fex(1,ima,ia) = zero
            fex(2,ima,ia) = zero
            fex(3,ima,ia) = zero
            tai  (ima)    = 0.5 * (talpj0 + talpj1)
            xgrt (ima)    = (talpj1 - talpj0)/drho
            xgrp (ima,ia) = xmks*0.666667*(walp(j+1)-walp(j))/drho
            xzi  (ima,ia) = 1.0
            sqz  (ima,ia) = 1.0 - amuai(ima)*dvedr/2.0
            if (sqz(ima,ia) .lt. sqzmin)
     .      sqz  (ima,ia) =  sqzmin
            sqz  (ima,ia) = (sqz(ima,ia))**1.5
            sqzmax        = MAX (sqzmax, sqz(ima,ia))
          end if
c
        end if
c
c -------------- end of fast ion contributions -------------------------
c
c       Set <E.B>:
c
        xedotb = 1.e-10
c
        k_order = 2     ! order of v-moments: u, p_q only
        k_potato = 0    ! do not include potato orbits
        m_i = mi        ! number of isotopes
        m_z = mz        ! highest charge state
        c_den = 1.e10   ! density cutoff
        c_potb = elong_r(1)*(btor/1.e4)/(2*q(1)**2)
        c_potl = q(1)*rmajor/100.D0
        p_b2 = xb2      ! <B**2>
        p_bm2 = xbm2    ! <B**-2>
        p_eb = xedotb   ! <E.B>
        p_fhat = xfhat  ! mu_0*F/dpsidr
        p_ft = xft      ! trapped particle fraction
        p_grbm2 = xgrbm2  ! <|grad psi|**2/B**2>
        p_ngrth = -xngrth ! <n.grad(theta)>
        p_grphi = 0.
        p_gr2phi = 0.
c       amu_i = amuai ! amuai(mxmi) in houlberg.i
c       grt_i = xgrt  ! xgrt(mxmi) in houlberg.i
c       temp_i = tai  ! tai(mxmi) in houlberg.i
c       den_iz = denai ! denai(mxmi,mxmz) in houlberg.i
c       fex_iz = fex  ! fex(3,mxmi,mxmz) in houlberg.i
c       grp_iz = xgrp ! xgrp(mxmi,mxmz) in houlberg.i
c
c  Thermal velocities-m/s:
c
        do ima=1,mi
          vtai(ima) = SQRT (2.0*xj7kv*tai(ima)/(amuai(ima)*promas))
        end do
c
c  Compute sum of (ni*d<z>/dte)*te/ne at j
c  dzdte(j,k) is computed in cray401 and in common block ions.i
c  Currently, dzdte set to zero in cray102 -> sumdz = 0
c
        enea = 0.5 * ( ene(j)+ ene(j+1))
        tea  = 0.5 * ( te(j)+ te(j+1))
        tia  = 0.5 * ( ti(j)+ ti(j+1))
c
        nia = 0.0   ! total thermal ion density
        sumdz = 0.0
        do k=1,nion
          dzdtea(k)  = 0.5 * (dzdte(j,k)+dzdte(j+1,k))
          ena(k)     = 0.5 * (en(j,k)+en(j+1,k))
          nia = nia + ena(k)
        enddo
        effz = enea / nia
        do k=1,nion
          sumdz = sumdz+dzdtea(k)*ena(k)
        enddo
        sumdz = sumdz*tea/enea
c
        if (lprintcl_test .gt. 0) then
          write(6,*) ' m_i = ',m_i
          write(6,*) ' m_z = ',m_z
          write(6,*) ' mx_mi = ',mx_mi
          write(6,*) ' mx_mz = ',mx_mz
          write(6,*) ' nion = ',nion
          write(6,*) ' q = ',q(j)
          write(6,*) ' btor = ',btor
          write(6,*) ' c_den = ',c_den
          write(6,*) ' c_potb = ',c_potb
          write(6,*) ' c_potl = ',c_potl
          write(6,*) ' p_b2 = ',p_b2
          write(6,*) ' p_bm2 = ',p_bm2
          write(6,*) ' p_eb = ',p_eb
          write(6,*) ' p_fhat = ',p_fhat
          write(6,*) ' p_fm = ',(p_fm(m), m=1,3)
          write(6,*) ' p_ft = ',p_ft
          write(6,*) ' p_grbm2 = ',p_grbm2
          write(6,*) ' p_ngrth = ',p_ngrth
          write(6,*) ' amu_i = ',(amuai(m), m=1,m_i)
          write(6,*) ' grt_i = ',(xgrt(m), m=1,m_i)
          write(6,*) ' temp_i = ',(tai(m), m=1,m_i)
          write(6,*) ' drho, r(j) = ',drho, r(j)
          do n=1,m_z
            write(6,*) ' den_iz = ',(denai(m,n), m=1,m_i)
          enddo
          do n=1,m_z
            write(6,*) ' fex_iz-1 = ',(fex(1,m,n), m=1,m_i)
            write(6,*) ' fex_iz-2 = ',(fex(2,m,n), m=1,m_i)
            write(6,*) ' fex_iz-3 = ',(fex(3,m,n), m=1,m_i)
          enddo
          do n=1,m_z
            write(6,*) ' grp_iz = ',(xgrp(m,n), m=1,m_i)
          enddo
        endif
c
c -------------------------------------------------------------------
c Now call NCLASS to get chi's, J_bs, etc. on the half-grid
c -------------------------------------------------------------------
c
        call nclass(k_order,k_potato,m_i,m_z,c_den,c_potb,c_potl,
     .              p_b2,p_bm2,p_eb,p_fhat,p_fm,p_ft,p_grbm2,
     .              p_grphi,p_gr2phi,p_ngrth,amuai,xgrt,tai,
     .              denai,fex,xgrp,m_s,jm_s,jz_s,p_bsjb,
     .              p_etap,p_exjb,calm_i,caln_ii,capm_ii,capn_ii,
     .              bsjbp_s,bsjbt_s,dn_s,gfl_s,qfl_s,sqz_s,upar_s,
     .              utheta_s,vn_s,veb_s,qeb_s,xi_s,ymu_s,chip_ss,
     .              chit_ss,dp_ss,dt_ss,iflag)
c
c -------------------------------------------------------------------
c  Find diffusivities and convective velocities.
c  Here, qfl_s is summed to get the total conduction and
c  the diagonal diffusivities and convective velocities are dq_s, vq_s
       sum(1)=zero
       do n=1,m_s
         im=jm_s(n)
         iza=iabs(jz_s(n))
         idum(n)=n
         sum(1)=rarray_sum(5,qfl_s(1,n),1)
         dq_s(n)=chit_ss(n,n)+chip_ss(n,n)
         vq_s(n)=sum(1)/(denai(im,iza)*1.6022e-16*tai(im))
     .           +dq_s(n)*xgrt(im)/tai(im)
       enddo
c
c  Compute effective diffusivities
c
       call rarray_zero(4,sum)
       do n=1,m_s
         im=jm_s(n)
         iz=jz_s(n)
         iza=iabs(iz)
         if (iz.gt.0) then
           sum(1)=sum(1)+(chit_ss(n,n)+chip_ss(n,n))*denai(im,iza)
           sum(2)=sum(2)+rarray_sum(5,qfl_s(1,n),1)/tai(im)/1.6022e-16
           sum(3)=sum(3)+denai(im,iza)
           sum(4)=sum(4)+(chit_ss(n,n)+chip_ss(n,n))*denai(im,iza)*
     .          xgrt(im)/tai(im)
         endif
       enddo
       ni_ncl(j)=sum(3)
       chieeff_ncl(j)=dq_s(1)*1.e4
       chiieff_ncl(j)=(sum(1)/sum(3))*1.e4
       veeff_ncl(j)=vq_s(1)
       vieff_ncl(j)=(sum(2)+sum(4))/sum(3)
c
c -------------------------------------------------------------------
c  Compute conductivity and diffusivity terms
c
       if (jneo .eq. 1) then
c
c..Coefficients for electron thermal transport
c  Note: n=1 is for electrons, n=2->m_s are for ions
c
c   chie: grad-n terms
c   d(nion+1,1,j) and d(nion+1,2,j) in 1/(cm-s) from
c   chip_ss(n,1) where n=2 to m_s for nion=2 case
c
       sum(1)=zero
       do n=2,m_s
         im=jm_s(n)
         iz=jz_s(n)
         iza=iabs(iz)
         sum(1) = chip_ss(n,1)*1.e4*denai(im,iza)*1.e-6
c        write(*,201) j,' xke-gradni = ',
c    .                sum(1), d(nion+1,n-1,j),
c    .                'for species i=',1,
c    .                ' driven by grad of j=',n
         d(nion+1,n-1,j)=d(nion+1,n-1,j)+sum(1)
       enddo
c
c   chie: grad-Te term
c   d(nion+1,nion+1,j)
c
c      write(*,200) j,' xke-gradTe = ', (dq_s(1)*1.e4)*enea,
c    .            d(nion+1,nion+1,j)
       d(nion+1,nion+1,j) = d(nion+1,nion+1,j) +
     .      dq_s(1)*1.e4*enea
c
c   chie: grad-Ti term
c   d(nion+1,nion+2,j)
c
       sum(1)=zero
       do n=2,m_s
         sum(1) = sum(1) +
     .      (chit_ss(n,1)+chip_ss(n,1))*1.e4*enea
       enddo
c
c      write(*,200) j,' xke-gradTi = ', sum(1),
c    .            d(nion+1,nion+2,j)
       d(nion+1,1,j) = d(nion+1,nion+2,j) + sum(1)
c
c..Coefficients for ion thermal transport
c
c  chii: grad-n terms
c  d(nion+2,1,j) and d(nion+2,2,j) from
c  sum over n of chip_ss(k,n) where k=2 to m_s for nion=2 case
c
       do k=2,m_s
         sum(1)=zero
         do n=2,m_s
           im=jm_s(n)
           iz=jz_s(n)
           iza=iabs(iz)
           sum(1) = sum(1) + chip_ss(k,n)*1.e4*denai(im,iza)*1.e-6
         enddo
c        if (k.eq.2 ) write(*,200) j,' xki-gradnh = ',
c    .                sum(1), d(nion+1,k-1,j)
         d(nion+2,k-1,j) = d(nion+1,k-1,j) + sum(1)
       enddo
c
c  chii: grad-Te term
c  d(nion+2,nion+1,j)
c
       sum(1)=zero
       do k=2,m_s
         im=jm_s(k)
         iz=jz_s(k)
         iza=iabs(iz)
         sum(1) = sum(1) +
     .        (chit_ss(1,k)+chip_ss(1,k))*1.e4*
     .        denai(im,iza)*1.e-6
       enddo
c
c      write(*,200) j,' xki-gradTe = ', sum(1),
c    .            d(nion+2,nion+1,j)
       d(nion+2,nion+1,j) = d(nion+2,nion+1,j) + sum(1)
c
c  chii: grad-Ti term
c  d(nion+2,nion+2,j)
c
       sum(2)=zero
       do k=2,m_s
         sum(1)=zero
         do n=2,m_s
           im=jm_s(n)
           iz=jz_s(n)
           iza=iabs(iz)
           sum(1) = sum(1) +
     .        (chit_ss(n,k)+chip_ss(n,k))*1.e4
         enddo
         sum(2)=sum(2)+sum(1)*denai(im,iza)*1.e-6
       enddo
c
c      write(*,200) j,' xki-gradTi = ', sum(2),
c    .            d(nion+2,nion+2,j)
       d(nion+2,nion+2,j) = d(nion+2,nion+2,j) + sum(2)
c
c..Coefficients for particle transport
c
c  di: grad-ni terms for hydrogenic and impurity particle transport
c  d(1,1,j) and d(1,2,j), d(2,1,j) and d(2,2,j) for nion=2 (m_s=3)
c
       do k=2,m_s
         do n=2,m_s
           sum(1) = dp_ss(n,k)*1.e4
c          write(*,201) j,' xdi-gradni = ',
c    .                sum(1), d(k-1,n-1,j),
c    .                'for species i=',k,
c    .                ' driven by grad of j=',n
           d(k-1,n-1,j)=d(k-1,n-1,j)+sum(1)
         enddo
       enddo
c
c  di: grad-Te term for hydrogenic and impurity particle transport
c  d(1,nion+1,j) and d(2,nion+1,j) for nion=2 (m_s=3)
c
       do k=2,m_s
         sum(1) = (dt_ss(1,k)+dp_ss(1,k))*1.e4
c        write(*,200) j,' xdi-gradTe = ',
c    .                sum(1), d(k-1,nion+1,j)
         d(k-1,nion+1,j)=d(k-1,nion+1,j)+sum(1)
       enddo
c
c  di: grad-Ti term for hydrogenic and impurity particle transport
c  d(1,nion+2,j) and d(2,nion+2,j) for nion=2 (m_s=3)
c
       do k=2,m_s
         sum(1)=zero
         do n=2,m_s
           sum(1) = sum(1) +
     .             (dt_ss(n,k)+dp_ss(n,k))*1.e4
         enddo
c        write(*,200) j,' xdi-gradTi = ',
c    .                sum(1), d(k-1,nion+2,j)
         d(k-1,nion+2,j)=d(k-1,nion+2,j)+sum(1)
       enddo
c
       endif
c
c      write(*,*) j,' chiieff = ',chiieff_ncl(j)
c      write(*,*) j,' chieeff = ',chieeff_ncl(j)
c      write(*,*) j,' dq_s = ', (dq_s(n),n=1,m_s)
c      write(*,*) j,' vq_s = ', (vq_s(n),n=1,m_s)
c
c -------------------------------------------------------------------
c    Transform the NCLASS bootstrap current <J_bs dot B> (A-T/m2)
c    into the form used by ONETWO (A-T/cm2) as well as convert
c    parallel resisitivity from (ohm*m) to (s) ...
c
c      convert units  A-T/m2 -> A-T/cm2:
C
       xjdotb = p_bsjb * 1.e-4
c
c      Divide by Bt0 to get bootstrap current in form used by ONETWO:
c      <J_bs dot B/Bt0> (A/cm2). Note that Btor is in Gauss.
c
       curbs(j) = xjdotb / (btor*1.0e-4)
c
       xjdotb_bt = 0.0
       xjdotb_bp = 0.0
       do n=1,m_s
         im=jm_s(n)
         iza=iabs(jz_s(n))
         xjdotb_bt=xjdotb_bt-bsjbt_s(n)*xgrt(im)/tai(im)
         xjdotb_bp=xjdotb_bp-bsjbp_s(n)*xgrp(im,iza)
     .             /denai(im,iza)/tai(im)
       enddo
       xjdotb_bt=xjdotb_bt * 1.e-4
       xjdotb_bp=xjdotb_bp * 1.e-4
c
       curbs_bt(j) = xjdotb_bt / (btor*1.0e-4)
       curbs_bp(j) = xjdotb_bp / (btor*1.0e-4)
c
c   Convert parallel resistivity from ohm*m to secs
c
       eta(j) = p_etap / 9.e9
c
c      if (j.eq.1) write(*,*) 'NCLASS boostrap current ...'
c      write(*,*) '  iflag = ', iflag
c      write(*,*) j,'  <J_bs dot B> (A-T/cm2) = ', xjdotb
c      write(*,*) j,'  <J_bs dot B> from grad-T (A-T/cm2) = ', xjdotb_bt
c      write(*,*) j,'  <J_bs dot B> from grad-P (A-T/cm2) = ', xjdotb_bp
c      write(*,*) j,'  eta-parallel (s) = ', eta(j)
c
c -------------------------------------------------------------------
c  Contributions to Faraday's law
c
c    In ONETWO, the bootstrap is calculated as
c    <jboot dot b/bto>= -s*(sum over ions {d(kfar,i,j)*dni/drho )
c                       +d(kfar,nion+1,j)*dte/drho
c                       +d(kfar,nion+2,j)*dti/drho
c                       +dcoef(j)*dnf/drho  )
c    where the sum over ions is a sum from i = 1 to i=nion and
c    includes primary as well as impurity species.
c    dcoef is a special coefficient used to model electron
c    density due to fast ions and fast ion density (optional)
c    and is calculated in subroutine SOURCE as dcoef = dfast+dfion
c    dfast is always used whereas dfion is the contribution from
c    fast ions.
c
c     Unit conversions -
c        1 A = 2.998e9 statamps
c        10**4 G = T, T = N/(C*m/s) = N/(A*m)
c        1.602e-8 A/cm**2 = 1 keV/(cm**4*gauss)
c        Multiply by xunits to change from keV/(cm**4 gauss)
c        to statamps/cm**2 : xunits = 1.602e-8 * 2.998e9
c        eta(s) = (1/9.e9)*eta(ohm*m)
c
c     kfar   = nion + 3, indexes Faraday's law.
c     s = 1/xres ... xres is a number: xres = cee * eta(j) / ra(j)
c     where cee=2.99792458e10 is the speed of light (cm/s), eta is
c     the resistivity (s), and ra is the rho value (cm).
c
       kfar=nion+3
       xres = cee * eta(j) / ra(j)    ! pure number
c
c --- Electron density gradient term.
c     It is used only to split out the electron density contribution in
c     the printout.  Otherwise, it is included in terms of the equivalent
c     ion density terms. Here, first species, s, are electrons.
c
       xdne(j) = xres*bsjbp_s(1)/denai(1,1)
c
c --- Ion (thermal) density gradient terms:
c     d(kfar,i,j), i=1,nion, units = statamps cm**2
c     The ion density gradient contribution to jboot is -sum over i of
c     {xres*xdin(j)*(dni/drho)}.
c
       xdni(j)=0.0
       do n=1,m_s
        im=jm_s(n)
        iza=iabs(jz_s(n))
        xdni(j)=xdni(j)+bsjbp_s(n)*(1.0+iza*denai(im,iza)/denai(1,1))
        d(kfar,n,j)=wneo(4,1)*2.998e11*xres*xdni(j)/(btor*nia)
c       write(*,*) j,' denai(im,iza) = ', denai(im,iza),
c    .             ' im = ',im, ' iza = ',iza, ' effz = ',effz
c       write(*,*) j,' ne = ',ene(j)*1.e6, ' ni = ',en(j,1)*1.e6,
c    .             ' nimp = ',en(j,2)*1.e6, ' nia = ', nia*1.e6
       end do
c
c --- Electron temp gradient term:
c     d(kfar,nion+1,j), units = statamp/(cm*keV)
c     The electron temp grad contribution to jboot is -sum over i of
c     {xres*xdte*(dte/drho)} where xdte(s)=bsjbt_s+bsjbp_s
c     Note: Assume first species, s, are electrons
c
       xdte(j)=bsjbt_s(1)+bsjbp_s(1)*(1.0+sumdz)
       d(kfar,nion+1,j) = wneo(4,2)*2.998e11*xres*xdte(j)/(btor*tai(1))
c
c --- Ion temperature gradient term:
c     d(kfar,nion+2,j), units = statamp/(cm*keV)
c     The ion temp gradient contribution to jboot is -sum over i of
c     {xres*xdti*(dti/drho)} where xdti(s)=bsjbt_s+bsjbp_s
c     Note: Assume equal Ti for all ions, tai(2)
c
       xdti(j)=0.0
       do n=2,m_s
         xdti(j)=xdti(j) + (bsjbt_s(n)+bsjbp_s(n))
       enddo
       d(kfar,nion+2,j) = wneo(4,3)*2.998e11*xres*xdti(j)/(btor*tai(2))
c
c      write(*,*) j,' d(kfar,1,j) = ',d(kfar,1,j),
c    .            ' wneo41 = ',wneo(4,1)
c      write(*,*) j,' d(kfar,2,j) = ',d(kfar,2,j),
c    .            ' wneo41 = ',wneo(4,1)
c      write(*,*) j,' d(kfar,nion+1,j) = ',d(kfar,nion+1,j),
c    .            ' wneo42 = ',wneo(4,2)
c      write(*,*) j,' d(kfar,nion+2,j) = ',d(kfar,nion+2,j),
c    .            ' wneo43 = ',wneo(4,3)
c
c --- Density gradient term of electrons due to fast ions:
c     dfast(j), units = statamps/cm**2
c     The ne gradient due to fast ions contributes to jboot as
c     -(xres*xdne(j)*(dnf/drho))*zbeam where it is
c     assumed that zbeam=1.
c
       dfast(j) = wneo(4,1) * xdne(j)
c
c --- Fast ion density gradient term:
c     dfion(j), units = statamps/cm**2
c     It is assumed that the fast ion term can be treated as part of the
c     thermal species where the bootstrap current contribution is
c     -(xres*dfion(j)*(dnf/drho))
c     Note: disabled for now
c
       dfion(j) = 0.0
c
c
      end do
c
c *** end of loop over mesh ***
c
 200   format(i2,a15,2x,2(1pe12.5,2x))
 201   format(i2,a15,2x,2(1pe12.5,2x),a14,i1,a21,i1)
c
      return
c
      end
