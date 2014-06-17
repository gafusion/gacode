      subroutine azimuth_integral (v, nv2, norder, sigvrd, itype)
c
      USE colrate
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c --- subroutine AZIMUTH_INTEGRAL constructs a table of
c ---         Integral from -1 to 1 {sigma(Erel)*vrel*P_n(Zeta)}
c --- P_n is the Legendre Polynomial of degree n
c --- the integration is over zeta = COS(theta), where theta is the angle
c --- between the two velocity vectors of the two distributions
c --- The result depends only on vrel which is symmetric in v1 and v2
c --- consequently we use SYMMETRIC STORAGE MODE in sigvr. The index
c --- ksym, corresponding to v1=v1_list(i) and v2=v2_list(j) is given by
c --- ksym=(i(i-1)/2)+j for i .ge. j
c
c --- input
c     v_list   in cm/sec of length nv1 and nv2 respectively.
c     norder   order of Legendre polynomial
c     itype    type of reaction
c
c --- output
c     sigvrd(i,j)      i=1,2..n1,j=0,2..norder the value of the reaction
c                      rate integral[cm/sec]*[units of sigma]
c
c REF. pg 117, Vol 2, Transport Notes.
c ------------------------------------------------- 11/30/95 --- HSJ ---
c
c      include 'colrate.i'    ! pick up xzeta,wzeta,nzeta,nzeta_set
c
      dimension v(*),sigvrd(1:n2vtble,0:nlegendre)
c
c     first get the weights for the zeta quadrature rule
c     if they do not yet exist:
c
      if (nzeta_set .eq. 0) then
        xx1 = -1.
        xx2 =  1.
        call gauleg (xx1, xx2, xzeta, wzeta, nzeta)
        nzeta_set = 1
      end if
c
c     now do the integral over zeta for speed v(i) in distribution #1
c     and speed v(j) in distribution #2:
c
      ksym=0
      do i=1,n2v
        do j=1,i
          ksym=ksym+1               ! symmetric storage mode
          sigvrd(ksym,norder)=0.0
          do k=1,nzeta
            if      (itype .eq. 1) then
              factr = sgv_ddnhe3(xzeta(k),v(j),v(i))
            else if (itype .eq. 2) then
              factr = sgv_dtnhe4(xzeta(k),v(j),v(i))
            else if (itype .eq. 3) then
              factr = sgv_tt2nhe4(xzeta(k),v(j),v(i))
            else if (itype .eq. 4) then
              factr = sgv_ddpt(xzeta(k),v(j),v(i))
            else
              print *,'itype =',itype
              call STOP ('subroutine AZIMUTH_INTEGRAL: bad ITYPE', 182)
            end if
            factr = factr * 78.95683 / (2*norder+1)  ! 8pisq/(2L+1)
            sigvrd(ksym,norder)=sigvrd(ksym,norder)+(v(i)*v(j))**2
     .           *(p_legendre(norder,xzeta(k))*wzeta(k) * factr)
          end do
        end do
      end do
      return
c
      end

      subroutine beam_beam_approx_rate (time, time0, timmax)
c
c
c ------------------------------------------------------------------
c     this subroutine calculates beam-beam fusion rates for the
c     following reactions:
c            D(D, n)He3
c            D(T, n)He4
c            T(T,2n)He4
c            D(D, p)T
c     subroutine  assumes that the fast ion distribution is given by
c     (for example) Gaffey, J. Plasma Physics, vol 16, 1976, pg 149
c     This subroutine uses the approximation that the fast ion
c     distribution doesnt change (as a function of te for example)
c     See sub beam_beam rate
c
c    INPUT
c    beam_beamddn_scale
c    beam_beamdtn_scale
c    beam_beamddp_scale
c    beam_beamtt2n_scale      these are set up in sub beam_beam_rate
c
c    OUTPUT
c
c         beam_beamddn(j,ksym)      the rates in #/cm**3sec
c         beam_beamddp(j,ksym)
c         beam_beamdtn(j,ksym)      for r(j=1) to r(j=nj)
c         beam_beamtt2n(j,ksym)     ksym is the symmetric storage
c                                   mode applied to 3*nbeams
c                                   fast ion distributions.
c                                   The total,summed over all distributions,
c                                   at each spatial mesh point j is in
c                                   element(j,ksym+1)
c
c ------------------------------------------------- 11/23/95 --- HSJ ---
c
      USE param
      USE fusion
      USE ions
      USE nub
      USE nub2
      USE soln
      USE neut
      USE numbrs
      USE mesh
      USE verbose
      USE machin
      USE geom
      USE colrate
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c      include 'param.i'
c      include 'fusion.i'
c      include 'geom.i'
c      include 'ions.i'
c      include 'machin.i'     ! btor
c      include 'mesh.i'
c      include 'neut.i'       ! pick up neutral density here
c      include 'nub.i'
c      include 'nub2.i'
c      include 'numbrs.i'
c      include 'soln.i'
c      include 'colrate.i'
c      include 'verbose.i'
c
      data  idfdf, idftf, itftf, ixvfast_set /0, 0, 0, 0/
      data  xkeverg, xmassp /6.241e+08, 1.6726e-24/
c
c     load these next two into /colrate/ (in file "colrate.i")
c
      twopi = 6.283185308 
      umdd  = 5.2175e-16       ! (keV/(cm/sec)**2), for d-d
c
      icxcalc   = icalc_cxfactor     ! load colrate.i
      mass_deut = 3.3435e-24         ! collect into a common area?
      mass_trit = 5.007289e-24
      mass_beam = xmassp*atw_beam
      nterms    = nlegendre ! terms retained in Legendre expansion of
c                             reaction rate integrals. This is not the
c                             same as # of terms retained in expansion
c                             of the distribution functions !!!
c
      call zeroa (beam_beamddn , kj*ksymbp)
      call zeroa (beam_beamddp , kj*ksymbp)
      call zeroa (beam_beamdtn , kj*ksymbp)
      call zeroa (beam_beamtt2n, kj*ksymbp)
      beam_beam_ddntot  = 0.0
      beam_beam_ddptot  = 0.0
      beam_beam_dtntot  = 0.0
      beam_beam_tt2ntot = 0.0
c
      if (beam_beam_fusion .eq.  0 )  return
      if (beamon(1) .ge. timmax)  return  ! beam not on in this run
c
      if (fusionvb .eq. 1)
     .  write (*, '(" calling BEAM_BEAM_APPROX_RATE")')
c
      do j=1,nj    ! loop over spatial mesh
c
         ksym = 0
c
c       loop over 3 energies for each beam
c
         do ie1=1,3*nbeams        ! first distribution
                 je = MOD (ie1, 3)
                 if (je .eq. 0) then
                    je=3
                    je1=ie1/3
                 else                       ! je1 is beam number,je is
                    je1=ie1/3+1             ! full half or third energy
                 end if                     ! of this beam
                 if (neg_ion_source(je1) .eq. 0 .or.
     .              (neg_ion_source(je1) .gt. 0 .and. je .eq. 1)) then
                 euplim1 = ebkev(je1)/je    ! no upscatter for now
c
c                sbpure = true (FREYA calced/saved) source rate,#/(cm**3sec):
c                sbsav has effects of charge exchange and time dependence
c                already in it, so must use sbpure if icalc_cxfactor = 1
c
                 if (icalc_cxfactor .eq. 0) then   ! source, sbd,..
c                                                  ..has cx factor in it
                     sbd1 = enb(j,je,je1)*taus(j)/taupb(j,je,je1)
                 else if (icalc_cxfactor .eq. 1) then   ! source, sbd,..
c                                                    ..without cx factor
                     sbd1 = sbpure(j,je,je1)*taus(j)
                 else
                     sbd1 = enb(j,je,je1)
                 end if
                 sbd1=sbd1/twopi  ! due to normalization of distribution
c
            do ie2=ie1,3*nbeams   ! second distribution
                 ksym=ksym+1
                 je = MOD (ie2, 3)
                 if (je .eq. 0) then
                    je2=ie2/3
                    je=3
                 else
                    je2=ie2/3+1
                 end if
                 if (neg_ion_source(je2) .eq. 0 .or.
     .              (neg_ion_source(je2) .gt. 0 .and. je .eq. 1)) then
c
c                sbpure = true (FREYA calc/saved) source rate,#/(cm**3sec):
c                sbsav has effects of charge exchange and time dependence
c                already in it, so must use sbpure if icalc_cxfactor = 1
c
                 if (icalc_cxfactor .eq. 0) then   ! source, sbd,..
c                                                  ..has cx factor in it
                     sbd2 = enb(j,je,je2)*taus(j)/taupb(j,je,je2)
                 else if (icalc_cxfactor .eq. 1) then   ! source, sbd,..
c                                                    ..without cx factor
                     sbd2 = sbpure(j,je,je2)*taus(j)
                 else
                     sbd2 = enb(j,je,je2)
                 end if
c
                 sbd2 = sbd2 / twopi      ! due to norm. of distribution
c
****           if (sbd1 * sbd2 .gt. 0.0) then     ! do only if necessary
c
                 factr = sbd1*sbd2*1.0e-27 ! 1.0e-27 sigma in millibarns
                 if (ie1 .eq. ie2)  factr = 0.5 * factr ! like-like..
c                                                       .. collisions
                 if (fusionvb .gt. 0 .and. j .eq. 1) then
                   write (*, '(a / 4(1pe14.3, 1x), 2(i3, 1x))')
     .                       ' sbd1, sbd2, factr, time, ie1, ie2',
     .                         sbd1, sbd2, factr, time, ie1, ie2
                 end if
                 if (fusionvb .gt. 0 .and. ksym .eq. 1) then
                   write (*, '(" tvalddn,factr =", 2(1pe14.3,2x))')
     .                           beam_beamddn_scale(j,ksym),factr
                 end if
                 beam_beamddn(j,ksym)  = beam_beamddn_scale(j,ksym)
     .                                 * factr
                 beam_beamddp(j,ksym)  = beam_beamddp_scale(j,ksym)
     .                                 * factr
                 beam_beamdtn(j,ksym)  = beam_beamdtn_scale(j,ksym)
     .                                 * factr
                 beam_beamtt2n(j,ksym) = beam_beamtt2n_scale(j,ksym)
     .                                 * factr
****           end if
               end if                  ! end neg_ion_source case
             end do
             end if                    ! end neg_ion_source case
          end do
c
c         store the sum of all beam progenated neutrons
c         in the last element, ksymp1 of each array
c
          ksymp1                  = ksym + 1
          beam_beamddn (j,ksymp1) = 0.0
          beam_beamddp (j,ksymp1) = 0.0
          beam_beamdtn (j,ksymp1) = 0.0
          beam_beamtt2n(j,ksymp1) = 0.0
c
          do jj=1,ksymp1-1
             beam_beamddn(j,ksymp1)  = beam_beamddn(j,ksymp1)+
     .                                 beam_beamddn(j,jj)
             beam_beamddp(j,ksymp1)  = beam_beamddp(j,ksymp1)+
     .                                 beam_beamddp(j,jj)
             beam_beamdtn(j,ksymp1)  = beam_beamdtn(j,ksymp1)+
     .                                 beam_beamdtn(j,jj)
             beam_beamtt2n(j,ksymp1) = beam_beamtt2n(j,ksymp1)+
     .                                 beam_beamtt2n(j,jj)
          end do
      end do
c
c     integrate the totals
c
      call trapv (r, beam_beamddn(1,ksymp1) , hcap, nj,
     .               beam_beam_ddntot)
      call trapv (r, beam_beamddp(1,ksymp1) , hcap, nj,
     .               beam_beam_ddptot)
      call trapv (r, beam_beamdtn(1,ksymp1) , hcap, nj,
     .               beam_beam_dtntot)
      call trapv (r, beam_beamtt2n(1,ksymp1), hcap, nj,
     .               beam_beam_tt2ntot)
c
      beam_beam_ddntot  = beam_beam_ddntot  * volfac ! volfac =
c                                                      4 * pisq * rmajor
      beam_beam_ddptot  = beam_beam_ddptot  * volfac
      beam_beam_dtntot  = beam_beam_dtntot  * volfac
      beam_beam_tt2ntot = beam_beam_tt2ntot * volfac
c
c     guard against single precision underflow (in subsequent programs)
c
      if (beam_beam_ddptot  .lt. 1.0e-30)  beam_beam_ddptot  = 0.0
      if (beam_beam_dtntot  .lt. 1.0e-30)  beam_beam_dtntot  = 0.0
      if (beam_beam_tt2ntot .lt. 1.0e-30)  beam_beam_tt2ntot = 0.0
      if (fusionvb .gt. 0) then
        write (*, '(" beam_beam_ddntot  =", 1pe14.3)') beam_beam_ddntot
        write (*, '(" beam_beam_ddptot  =", 1pe14.3)') beam_beam_ddptot
        write (*, '(" beam_beam_dtntot  =", 1pe14.3)') beam_beam_dtntot
        write (*, '(" beam_beam_tt2ntot =", 1pe14.3)') beam_beam_tt2ntot
      end if
      return
c
      end

      subroutine beam_beam_intg1 (vlow1, vup1, vfast_low_lim2,
     .                            vfast_up_lim2, vbeam13, vbeam23,
     .                            beam_pitch_angle1, beam_pitch_angle2,
     .                            tvalddn, tvaldtn, tvaltt2n,tvalddp)
c
      USE colrate
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ------------------------------------------------- 11/23/95 --- HSJ ---
c     returns the reaction rates
c         valddn    d(d,n)he3 if ib_d  =1 in colrate.i
c         valdtn    d(t,n)he4    ib_dt =1
c         valtt2n   t(t,2n)he4   ib_t  =1
c
c ----------------------------------------------------------------------
c
c      include 'colrate.i'
c
      tvalddn=0.0
      tvalddp=0.0
      tvaldtn=0.0
      tvaltt2n=0.0
      dv=vup1-vlow1
      do l=0,nterms            ! loop over Legendre Polynomial expansion
         valddn=0.0
         valddp=0.0
         valdtn=0.0
         valtt2n=0.0
         do j=1,nvfast     ! evaluate first integral over beam ion speed
                v1scaled=xvfast(j)*dv+vlow1
                v1scaled3=v1scaled*v1scaled*v1scaled
                if (icxcalc .eq. 1) then
                  vb = vup1       ! set upper integration limit in cxint
                  factr=cxint(v1scaled)/(v1scaled3+vcrit3) ! source does
c                                                    not include cx part
                 else
                   factr = 1.0/(v1scaled3+vcrit3)          ! source does
c                                                        include cx part
                end if
                factr = factr * (2.*l+1.)*0.5
                factr = factr * ((v1scaled3/vbeam13)*(vbeam13+vcrit3)
     .                        /(v1scaled3+vcrit3)) **(l*(l+1.)*z2/6.)
                call beam_beam_intg2(v1scaled,l,vfast_low_lim2,
     .                               vfast_up_lim2, vbeam23,
     .                               beam_pitch_angle2,valddnl,
     .                               valdtnl,valtt2nl,valddpl)
                if (ib_d .eq.1 ) then
                    valddn=valddn+valddnl*wxvfast(j)*factr
                    valddp=valddp+valddpl*wxvfast(j)*factr
                end if
                if (ib_dt .eq.1 ) then
                     valdtn=valdtn+valdtnl*wxvfast(j)*factr
                end if
                if (ib_t .eq.1 ) then
                    valtt2n=valtt2n+valtt2nl*wxvfast(j)*factr
                end if
          end do
c
          pleg     = p_legendre(l,beam_pitch_angle1)
          tvalddn  = valddn*pleg +tvalddn
          tvalddp  = valddp*pleg +tvalddp
          tvaldtn  = valdtn*pleg +tvaldtn
          tvaltt2n = valtt2n*pleg +tvaltt2n
      end do
      tvalddn  = tvalddn*dv
      tvalddp  = tvalddp*dv
      tvaldtn  = tvaldtn*dv
      tvaltt2n = tvaltt2n*dv
      return
c
      end

      subroutine beam_beam_intg2 (v1scaled, norder, vlow2, vup2,
     .                            vbeam23, beam_pitch_angle2, valddnl,
     .                            valdtnl, valtt2nl,valddpl)
c
      USE colrate
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ------------------------------------------------- 11/23/95 --- HSJ ---
c
c      include 'colrate.i'
c
      dv       = vup2-vlow2
      valddnl  = 0.0
      valdtnl  = 0.0
      valtt2nl = 0.0
      valddpl  = 0.0
      do j=1,nvfast  ! evaluate second integral over beam ion speed
         v2scaled = xvfast(j)*dv+vlow2
c
c        get ddn,dtn,and tt2n rates simultaneoulsy:
c
         call beam_beam_zeta_intg (v1scaled,v2scaled,norder,valddnll,
     .                            valdtnll,valtt2nll,valddpll)
c
c        the charge exchange integral, cxint,
c        is not precomputed because it depends on te (through vcrit).
c        Thus we would need a separate table in each region
c        and would have to recompute anyway each time te changes
c
         v2scaled3 = v2scaled*v2scaled*v2scaled
         if (icxcalc .eq. 1) then
           vb    = vup2           ! set upper integration limit in cxint
           factr = cxint(v2scaled)/(v2scaled3+vcrit3) ! source does not
c                                                       include cx part
         else
           factr = 1.0/(v2scaled3+vcrit3)              ! source does
c                                                       includes cx part
         end if
****     factr = factr*p_legendre(norder,beam_pitch_angle2)
         factr = factr*(2.*norder+1.)*0.5
         factr = factr*((v2scaled3/vbeam23)*(vbeam23+vcrit3)
     .                 /(v2scaled3+vcrit3))**(norder*(norder+1)*z2/6.0)
         if (ib_d .eq. 1) then
             valddnll  = valddnll*wxvfast(j)*factr
             valddnl   = valddnl+valddnll
             valddpll  = valddpll*wxvfast(j)*factr
             valddpl   = valddpl+valddpll
         end if
         if (ib_dt .eq. 1) then
             valdtnll  = valdtnll*wxvfast(j)*factr
             valdtnl   = valdtnl+valdtnll
         end if
         if (ib_t .eq. 1) then
             valtt2nll = valtt2nll*wxvfast(j)*factr
             valtt2nl  = valtt2nl+valtt2nll
         end if
      end do
!      pleg    = p_legendre(norder,beam_pitch_angle2)
      pleg     = 1.0d0 ! Only norder =0 contributes
      valddnl  = dv*pleg*valddnl          ! due to scaling of weights
      valddpl  = dv*pleg*valddpl
      valdtnl  = dv*pleg*valdtnl
      valtt2nl = dv*pleg*valtt2nl
      return
c
      end

      subroutine beam_beam_rate (time, time0, timmax, totcur1)
c

c
c ----------------------------------------------------------------------
c     this subroutine calculates beam-beam fusion rates for the
c     following reactions:
c            D(D, n)He3
c            D(T, n)He4
c            T(T,2n)He4
c            D(D, p)T
c     subroutine assumes that the fast ion distribution is given by
c     (for example) Gaffey, J.Plasma Physics, vol 16, 1976, pg 149
c
c    INPUT
c    time
c    time0
c    timmax        times used to establish when these calcs are to be done
c    totcur(1)     used only for sign info on beam pitch angle
c
c    OUTPUT
c
c         beam_beamddn(j,ksym)      the rates in #/cm**3
c         beam_beamddp(j,ksym)
c         beam_beamdtn(j,ksym)      for r(j = 1) to r(j = nj)
c         beam_beamtt2n(j,ksym)     ksym is the symmetric storage
c                                   mode applied to 3*nbeams
c                                   fast ion distributions.
c                                   The total,summed over all distributions,
c                                   at each spatial mesh point j is in
c                                   element(j,ksym+1)
c
c ------------------------------------------------- 11/23/95 --- HSJ ---
c
      USE param
      USE fusion
      USE ions
      USE neut
      USE nub
      USE nub2
      USE soln
      USE numbrs
      USE mesh
      USE verbose
      USE machin
      USE geom
      USE colrate
      implicit  integer (i-n), real*8 (a-h, o-z)
c      include 'param.i'
c      include 'fusion.i'
c      include 'geom.i'
c      include 'ions.i'
c      include 'machin.i'     ! btor
c      include 'mesh.i'
c      include 'neut.i'       ! pick up neutral density here
c      include 'nub.i'
c      include 'nub2.i'
c      include 'numbrs.i'
c      include 'soln.i'
c      include 'colrate.i'
c      include 'verbose.i'
c
      data  idfdf, idftf, itftf, ixvfast_set, imod_calc
     .     /0    , 0    , 0    , 0          , 0/
      data  xkeverg  , xmassp
     .     /6.241e+08, 1.6726e-24/
c
      icxcalc   = icalc_cxfactor     ! load colrate.i
      mass_deut = 3.3435e-24         ! collect into a common area?
      mass_trit = 5.007289e-24
      mass_beam = xmassp*atw_beam
      nterms    = nlegendre ! terms retained in Legendre expansion of
c                             reaction rate integrals. This is not the
c                             same as # of terms retained in expansion
c                             of the distribution functions !!!
c
      call zeroa (beam_beamddn , kj*ksymbp)
      call zeroa (beam_beamddp , kj*ksymbp)
      call zeroa (beam_beamdtn , kj*ksymbp)
      call zeroa (beam_beamtt2n, kj*ksymbp)
      beam_beam_ddntot    =  0.0
      beam_beam_ddptot    =  0.0
      beam_beam_dtntot    =  0.0
      beam_beam_tt2ntot   =  0.0
      beam_beam_long_calc = -1
c
      if (beam_beam_fusion .eq.  0                       )  return
c
****  if (beam_beam_fusion .eq. -1 .and. time .gt. time0 )  return
****  if (beam_beam_fusion .eq. -2 .and. time .lt. timmax)  return
****  if (beam_beam_fusion .eq. -3 ) then
****    if (.not.(time .eq. time0 .or. time .ge. timmax))   return
****  end if
c
      if (beamon(1) .ge. timmax)  return  ! beam not on in this run
      if (beam_beam_fusion .gt. 0) then
         imod_calc=imod_calc+1
         if (imod_calc .eq. 1) go to 10 ! always do it on the first call
         if (MOD (imod_calc, beam_beam_fusion) .eq. 0)  go to 10
         if (time .ge. timmax) go to 10 ! do it at last time point
         return
      end if
c
   10 if (fusionvb .eq. 1)  write (*, '(" IN  BEAM_BEAM_RATE")')
c
c     set up the distribution independent tables,beam_beam_*,
c     the first time this routine is called:
c     even though there are nbeams,the species in each is the same:
c
      ib_d  = 0
      ib_dt = 0
      ib_t  = 0
      if (nameb .eq. 'd' ) ib_d = 1
      if (nameb .eq. 't' ) ib_t = 1
      if (nameb .eq. 'dt') then
        ib_dt = 1
        ib_d  = 1
        ib_t  = 1
      end if
      ebeammax_table = 0.0
      do jb=1,nbeams
        ebeammax_table = MAX (ebeammax_table, ebkev(jb))
      end do
      ebeammax_table = 1.05*ebeammax_table
      ebeammin_table = 0.01 ! keV
      ebeammin_table = 0.001 ! keV  HSJ 5/1/01
      vbeammax_table = SQRT (2.0 * ebeammax_table
     .               / (xkeverg * mass_beam)) ! cm/sec, upper limit
      vbeammin_table = SQRT (2.0 * ebeammin_table
     .               / (xkeverg * mass_beam)) ! cm/sec, lower limit
      dv             = (vbeammax_table-vbeammin_table)/(n2v-1)
      do j=1,n2v
        vbeam_beam(j) = vbeammin_table+(j-1)*dv
      end do
      vbeam_beam(n2v) = vbeammax_table
c
      do j=0,nterms            ! loop over Legendre polynomials
c
c       if the beam has deuterium in it, get ddnhe3 and ddpt tables:
c
        if (ib_d .eq. 1 .and. idfdf .eq. 0) then
          call azimuth_integral (vbeam_beam,n2v, j, beam_beam_ddnhe3, 1)
          call azimuth_integral (vbeam_beam,n2v, j, beam_beam_ddpt  , 4)
        end if
c
c       if the beam has deuterium and tritium in it, get dtnhe4 tables:
c
        if (ib_dt .eq. 1 .and. idftf .eq. 0)
     .    call azimuth_integral (vbeam_beam,n2v,j,beam_beam_dtnhe4,2)
c
c       if the beam has tritium in it,get tt2nhe4 tables:
c
        if (ib_t .eq. 1 .and. itftf .eq. 0)
     .    call azimuth_integral (vbeam_beam,n2v,j,beam_beam_tt2nhe4,3)
      end do
c
      idfdf = 1
      itftf = 1
      idftf = 1
c
c     get the Gauss-Legendre weights and abscissa values for the
c     integration over the beam ion speed. Since the lower and
c     upper limits of integration can vary (due to beam turn on
c     and off and beam energy components we get weights
c     scaled to [0,1] and do the transformation
c     when we do the integration. xvfast and wxvfast are passed in
c     colrate.i:
c
      if (ixvfast_set .eq. 0) then
         ixvfast_set=1
         xx1 = 0.0
         xx2 = 1.0
         call gauleg (xx1, xx2, xvfast, wxvfast, nvfast)
      end if
      if (fusionvb .eq. 1) 
     .       write (*, '(" IN  BEAM_BEAM_RATE,loop over j")')
c
      do j=1,nj                    ! loop over spatial mesh
            z1     = 0.0
            z2     = 0.0
            do i=1,nion
                 z12 = en(j,i)*zsq(j,i)
                 z1  = z1+z12/atw(i)
                 z2  = z2+z12
            end do
            z1 = z1*atw_beam/ene(j)
            z2 = z2/(z1*ene(j))    ! effective charge of bg ions
            den_neut = 0.0
            do i=1,2               ! nneut = 2, hardcoded everywhere
              den_neut = den_neut + enn(j,i)        ! load colrate.i
            end do
c
c           vthion is energy at which a fast ion is considered
c           thermalized and therefore no longer part of the fast
c           ion distribution function:
c
            vthion  = 4.0026e-5 * SQRT (2.0*ti(j)/mass_beam) ! cm/sec
            vthion  = vbeammin_table                         ! CHANGE???
            vthelec = 1.3256e9 * SQRT (2.0*te(j))            ! cm/sec
            vcrit   = 0.09*((z1/atw(ibion))**0.33)*vthelec
            vcrit3  = vcrit**3
            tauslocal = taus(j)            ! sec, stored in colrate.i,..
c                                          ..used in cxint
            ksym      = 0
c
c       loop over 3 energies for each beam
c
         do ie1=1,3*nbeams                 ! first distribution
                 je = MOD (ie1, 3)
                 if (je .eq. 0) then
                    je=3
                    je1=ie1/3
                 else                      ! je1 is beam number,je is
                    je1=ie1/3+1            ! full half or third energy
                 end if                    ! of this beam
                 euplim1  = ebkev(je1)/je  ! no upscatter for now
c
                 if (neg_ion_source(je1) .eq. 0 .or.
     .              (neg_ion_source(je1) .gt. 0  .and. je .eq. 1)) then
c
c                sbpure = true (FREYA calced/saved) source rate,#/(cm**3sec):
c                sbsav has effects of charge exchange and time dependence
c                already in it, so must use sbpure if icalc_cxfactor = 1
c
                 if (icalc_cxfactor .eq. 0) then   ! source, sbd,..
c                                                  ..has cx factor in it
                     sbd1 = enb(j,je,je1)*taus(j)/taupb(j,je,je1)
                 else if (icalc_cxfactor .eq. 1) then   ! source, sbd,..
c                                                    ..without cx factor
                     sbd1 = sbpure(j,je,je1)*taus(j)
                 else
                     sbd1 = enb(j,je,je1)
                 end if
                 sbd1=sbd1/twopi ! due to normalization of distribution
                 beam_pitch_angle1 = zeta(j,je,je1) * btor * totcur1
     .                                          / ABS (btor * totcur1)
c
                 vbeam1 = SQRT (2.*euplim1/(xkeverg*mass_beam)) ! cm/sec
                 call set_beam_limits (vbeam1,vcrit,taus(j),beamon(1),
     .             beam_end(1),time,vthion,vfast_low_lim1,vfast_up_lim1,
     .                 icalc_cxfactor)
c
            do ie2=ie1,3*nbeams            ! second distribution
                 ksym = ksym+1
                 je   = MOD (ie2, 3)
                 if (je .eq. 0) then
                    je2 = ie2/3
                    je  = 3
                 else
                    je2 = ie2/3+1
                 end if
c
                 if (neg_ion_source(je2) .eq. 0 .or.
     .              (neg_ion_source(je2) .gt. 0 .and. je .eq. 1)) then
c
                 euplim2 = ebkev(je2)/je  ! no upscatter for now
c
c                sbpure = true (FREYA calc/saved) source rate,#/(cm**3sec):
c                sbsav has effects of charge exchange and time dependence
c                already in it, so must use sbpure if icalc_cxfactor = 1
c
                 if (icalc_cxfactor .eq. 0) then   ! source, sbd,..
c                                                  ..has cx factor in it
                     sbd2 = enb(j,je,je2)*taus(j)/taupb(j,je,je2)
                 else if (icalc_cxfactor .eq. 1) then   ! source, sbd,..
c                                                    ..without cx factor
                     sbd2 = sbpure(j,je,je2)*taus(j)
                 else
                     sbd2 = enb(j,je,je2)
                 end if
c
                sbd2 = sbd2/twopi ! due to normalization of distribution
****            if (sbd1*sbd2 .gt. 0.0) then  ! do only if necessary
                beam_pitch_angle2 = zeta(j,je,je2) * btor * totcur1
     .                                        / ABS (btor * totcur1)
c
                vbeam2 = SQRT (2.*euplim2/(xkeverg*mass_beam)) ! cm/sec
                call set_beam_limits (vbeam2,vcrit,taus(j),beamon(1),
     .             beam_end(1),time,vthion,vfast_low_lim2,vfast_up_lim2,
     .                             icalc_cxfactor)
c
c               do all the beam-beam reaction rate calculations simultaneously
c
                vbeam13 = vbeam1*vbeam1*vbeam1
                vbeam23 = vbeam2*vbeam2*vbeam2
                call beam_beam_intg1 (vfast_low_lim1,vfast_up_lim1,
     .                                vfast_low_lim2,vfast_up_lim2,
     .                                vbeam13,vbeam23,
     .                              beam_pitch_angle1,beam_pitch_angle2,
     .                              tvalddn,tvaldtn,tvaltt2n,tvalddp)
                factr = sbd1*sbd2*1.0e-27  ! 1.0e-27 sigma in millibarns
                if (ie1 .eq. ie2)
     .            factr = 0.5 * factr             ! like-like collisions
                if (fusionvb .gt. 0 .and. j .eq. 1) then
                   write (*, '(a / 6(1pe14.3, 1x), 2(i3, 1x))')
     .         ' sbd1, sbd2, factr, time, euplim1, euplim2, ie1, ie2 =',
     .           sbd1, sbd2, factr, time, euplim1, euplim2, ie1, ie2
                end if
c                if (fusionvb .gt. 0 .and. ksym .eq. 1) then
                if (fusionvb .gt. 0 ) then
                  write (*, '(a / 4(1pe14.3, 2x),i3)')
     .            ' tvalddn, factr, den_neut, tauslocal,j  =',
     .              tvalddn, factr, den_neut, tauslocal,j
                end if
                beam_beamddn (j,ksym) = tvalddn*factr*fdbeam*fdbeam
                beam_beamddp (j,ksym) = tvalddp*factr*fdbeam*fdbeam
                beam_beamdtn (j,ksym) = tvaldtn*factr*fdbeam
     .                                                *(1.0-fdbeam)
                beam_beamtt2n(j,ksym) = tvaltt2n*factr*(1.0-fdbeam)
     .                                                *(1.0-fdbeam)
c
c --- the beam beam rate is usually a small contribution to the total.
c --- consequently it may not be necessary to calculate it very
c --- accurately. If we neglect changes in fast ion charge exchange
c --- and assume that the initial critical velocity doesnt change too
c --- much then the integrals performed above are in fact constant,
c --- we need only scale by the parameter factr defined above. So let's
c --- save the relevant info:
c
                beam_beamddn_scale(j,ksym)=tvalddn*fdbeam*fdbeam
                beam_beamddp_scale(j,ksym)=tvalddp*fdbeam*fdbeam
                beam_beamdtn_scale(j,ksym)=tvaldtn*fdbeam*
     .                                                   (1.0-fdbeam)
                beam_beamtt2n_scale(j,ksym)=tvaltt2n*(1.-fdbeam)*
     .                                                   (1.0-fdbeam)
                beam_beam_long_calc=1
c
c ---  to use the above arrays rather than do the integrals
c ---  just multiply by factr as defined above. HSJ
c
****            end if          ! sbd1*sbd2 = 0.0 clause
             end if             ! neg_ion_source, second beam
             end do
             end if             ! neg_ion_source, second beam
          end do
c
c         store the sum of all beam progenated fast ion distributions
c         in the last element, ksymp1 of each array
c
          ksymp1                  = ksym+1
          beam_beamddn(j,ksymp1)  = 0.0
          beam_beamddp(j,ksymp1)  = 0.0
          beam_beamdtn(j,ksymp1)  = 0.0
          beam_beamtt2n(j,ksymp1) = 0.0
          do jj=1,ksymp1-1
             beam_beamddn(j,ksymp1)  = beam_beamddn(j,ksymp1)+
     .                                 beam_beamddn(j,jj)
             beam_beamddp(j,ksymp1)  = beam_beamddp(j,ksymp1)+
     .                                 beam_beamddp(j,jj)
             beam_beamdtn(j,ksymp1)  = beam_beamdtn(j,ksymp1)+
     .                                 beam_beamdtn(j,jj)
             beam_beamtt2n(j,ksymp1) = beam_beamtt2n(j,ksymp1)+
     .                                 beam_beamtt2n(j,jj)
          end do
      end do
c
c     integrate the totals
c
      call trapv (r, beam_beamddn(1,ksymp1) , hcap, nj,
     .                               beam_beam_ddntot)
      call trapv (r, beam_beamddp(1,ksymp1) , hcap, nj,
     .                               beam_beam_ddptot)
      call trapv (r, beam_beamdtn(1,ksymp1), hcap, nj,
     .                               beam_beam_dtntot)
      call trapv (r, beam_beamtt2n(1,ksymp1), hcap, nj,
     .                               beam_beam_tt2ntot)
c
      beam_beam_ddntot  = beam_beam_ddntot  * volfac ! volfac =
c                                                      4 * pisq * rmajor
      beam_beam_ddptot  = beam_beam_ddptot  * volfac
      beam_beam_dtntot  = beam_beam_dtntot  * volfac
      beam_beam_tt2ntot = beam_beam_tt2ntot * volfac
c
c     guard against single precision underflow (for subsequent programs):
c
      if (beam_beam_ddptot  .lt. 1.0e-30)  beam_beam_ddptot  = 0.0
      if (beam_beam_dtntot  .lt. 1.0e-30)  beam_beam_dtntot  = 0.0
      if (beam_beam_tt2ntot .lt. 1.0e-30)  beam_beam_tt2ntot = 0.0
      if (fusionvb .gt. 0) then
        write (*, '(" beam_beam_ddntot  =", 1pe14.3)') beam_beam_ddntot
        write (*, '(" beam_beam_ddptot  =", 1pe14.3)') beam_beam_ddptot
        write (*, '(" beam_beam_dtntot  =", 1pe14.3)') beam_beam_dtntot
        write (*, '(" beam_beam_tt2ntot =", 1pe14.3)') beam_beam_tt2ntot
      end if
      return
c
      end

      subroutine beam_beam_zeta_intg (v1s, v2s, norder, valddnll,
     .                                valdtnll, valtt2nll,valddpll)
c
      USE colrate
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c --- evaluate the beam ion dist.
c --- times the integral defined by azimuth_integral
c ------------------------------------------------------------------ HSJ
c
c      include 'colrate.i'
c
c --- get the value of the integral over angle between the two velocity
c --- vectors by interpolation in the tables sig.
c
      call tableintrp (vbeam_beam, n2v, v1s, j1)
      if (j1 .le. 0 .or. j1 .ge. n2v) then
        call STOP ('subroutine BEAM_BEAM_ZETA_INTG: problem #1', 200)
      end if
      call tableintrp (vbeam_beam, n2v, v2s, j2)
      if (j2 .le. 0 .or. j2 .ge. n2v) then
        call STOP ('subroutine BEAM_BEAM_ZETA_INTG: problem #2', 201)
      end if
c
      dv1  = v1s-vbeam_beam(j1)
      dv2  = v2s-vbeam_beam(j2)
      a1   = dv1*dv2
      a2   = (vbeam_beam(j1+1)-v1s)*dv2
      a3   = (vbeam_beam(j1+1)-v1s)*(vbeam_beam(j2+1)-v2s)
      a4   = dv1*(vbeam_beam(j2+1)-v2s)
      area = a1+a2+a3+a4
c
c     use quasilinear interpolation in the precomputed beam_beam_*
c     arrays
c     get the index of the four points that form the square that
c     embraces the point (v1,v2). Recall that we are using symmetric
c     storage here:
c
      js    = MIN (j1, j2)
      jh    = MAX (j1, j2)
      k1sym = (jh * (jh-1)) / 2 + js
      j1    = j1 + 1
      js    = MIN (j1, j2)
      jh    = MAX (j1, j2)
      k2sym = (jh * (jh-1)) / 2 + js
      j2    = j2 + 1
      js    = MIN (j1, j2)
      jh    = MAX (j1, j2)
      k3sym = (jh * (jh-1)) / 2 + js
      j1    = j1 - 1
      js    = MIN (j1, j2)
      jh    = MAX (j1, j2)
      k4sym = (jh * (jh-1)) / 2 + js
c
      if (ib_d .eq. 1) then    ! if beam has d get ddn reaction
          val1=beam_beam_ddnhe3(k1sym,norder)
          val2=beam_beam_ddnhe3(k2sym,norder)
          val3=beam_beam_ddnhe3(k3sym,norder)
          val4=beam_beam_ddnhe3(k4sym,norder)
          valddnll  = (a1*val3+a2*val4+a3*val1+a4*val2)/area
          val1=beam_beam_ddpt(k1sym,norder)
          val2=beam_beam_ddpt(k2sym,norder)
          val3=beam_beam_ddpt(k3sym,norder)
          val4=beam_beam_ddpt(k4sym,norder)
          valddpll  = (a1*val3+a2*val4+a3*val1+a4*val2)/area
      else
          valddnll = 0.0
          valddpll = 0.0
      end if
c
      if (ib_dt .eq. 1) then    ! if beam has d and t get dtn reaction
          val1=beam_beam_dtnhe4(k1sym,norder)
          val2=beam_beam_dtnhe4(k2sym,norder)
          val3=beam_beam_dtnhe4(k3sym,norder)
          val4=beam_beam_dtnhe4(k4sym,norder)
          valdtnll  = (a1*val3+a2*val4+a3*val1+a4*val2)/area
      else
          valdtnll = 0.0
      end if
c
      if (ib_t .eq. 1) then   ! if beam has t get tt2n reaction
          val1=beam_beam_tt2nhe4(k1sym,norder)
          val2=beam_beam_tt2nhe4(k2sym,norder)
          val3=beam_beam_tt2nhe4(k3sym,norder)
          val4=beam_beam_tt2nhe4(k4sym,norder)
          valtt2nll  = (a1*val3+a2*val4+a3*val1+a4*val2)/area
      else
          valtt2nll = 0.0
      end if
c
      return
c
      end

      subroutine beam_thermal_approx_fus (time, timmax, qbfus, sbfus)
c

c
c -------------------------------------------------- 12/5/95 --- HSJ ---
c     this subroutine calculates beam-thermal fusion rates for the
c     following reactions:             output vector:
c     --------------------            ---------------------
c          Df(Dth,n)He3                  beam_thermalddn
c          Df(Dth,p)T                    beam_thermalddp
c          Df(Tth,n)He4                  beam_thermaltth_df
c          Tf(Dth,n)He4                  beam_thermaldth_tf
c          Tf(Tth,2n)He4                 beam_thermaltt2n
c
c     CROSS SECTIONS are from
c     Bosch & Hale, Nuc. Fus., vol32, no.4 (1992) 611
c
c     INPUT:
c     iddfus   =  0    no thermal d,t or dt mixture
c              =  1    thermal species is d
c              =  2    thermal species is dt, use fraction of dt density
c                      given by fd  to specify d density
c              =  3    both d and t are present as separate thermal species
c              =  4    t is a thermal species,d is not present
c              =  5    allows two fluid d,t thermal and single
c                      fluid  'dt' beam
c
c     beam_thermalddn_scale(j,ksym)
c     beam_thermalddp_scale(j,ksym)
c     beam_thermaltth_df_scale(j,ksym)
c     beam_thermaltt2n_scale(j,ksym)
c     beam_thermaldth_tf_scale(j,ksym)
c
c     OUTPUT:
c     sbfus(j)  beam fusion rate density counts total of
c               all d(t,n)he4 reactions
c     qbfus(j)  alpha power density from sbfus
c
c         if beam d and thermal d are present:
c         beam_thermalddn(j,k)     #/cm**3/sec
c                                  j=1,2,..nj,k=1,2,..3*nbeams
c             beam_thermal_ddntot       volume integrated rate, #/sec
c         if beam d and thermal t are prsent:
c         beam_thermaltth_df(j,k)
c             beam_thermal_dtntot       volume integrated rate, #/sec
c         if beam t and therma t are present:
c         beam_thermaltt2n(j,k)
c             beam_thermal_tt2ntot       volume integrated rate, #/sec
c         if beam t and thermal d are present
c         beam_thermaldth_tf(j,k)
c             beam_thermaldth_tftot      volume integrated rate, #/sec
c
c ----------------------------------------------------------------------
c
      USE param
      USE fusion
      USE ions
      USE neut
      USE nub
      USE nub2
      USE soln
      USE numbrs
      USE mesh
      USE verbose
      USE machin
      USE geom
      USE colrate
      USE tordlrot
      implicit  integer (i-n), real*8 (a-h, o-z)
c      include 'param.i'
c      include 'fusion.i'
c      include 'geom.i'
c      include 'ions.i'
c      include 'machin.i'     ! btor
c      include 'mesh.i'
c      include 'neut.i'       ! pick up neutral density here
c      include 'nub.i'
c      include 'nub2.i'
c      include 'numbrs.i'
c      include 'soln.i'
c      include 'tordlrot.i'
c      include 'colrate.i'
c      include 'verbose.i'
c
      real*8  qbfus(*),sbfus(*)
c
      if (beam_thermal_fusion .eq. 0     )  return
      if (             iddfus .eq. 0     )  return
      if (             beamon(1) .ge. timmax)  return ! beam is not on..
c                                                  ..in this run
      if (           fusionvb .gt. 0     )
     .  write (*, '(" in BEAM_THERMAL_APPROX_FUS")')
c
      do j=1,nj
          k         = 0
          do jb=1,nbeams
             iec=3
             if (neg_ion_source(jb) .gt. 0)  iec = 1
            do jc=1,iec
              k = k + 1
c
c                  sbpure = true (FREYA calced/saved) source rate,#/(cm**3sec):
c                  sbsav has effects of charge exchange and time dependence
c                  already in it, so must use sbpure if icalc_cxfactor = 1
c
                   if (icalc_cxfactor .eq. 0) then
                     sbd = enb(j,jc,jb)*taus(j)/taupb(j,jc,jb)
                   else if (icalc_cxfactor .eq. 1) then ! source, sbd,..
c                                                    ..without cx factor
                     sbd = sbpure(j,jc,jb)*taus(j)
                   else               ! source, sbd, has cx factor in it
                     sbd = enb(j,jc,jb)
                   end if
c
                    if (iddfus .eq. 1) then
                      endth = en(j,id)
                      entth = 0.0
                    end if
c
                    if (iddfus .eq. 2) then
                      endth = fd*en(j,idt)
                      entth = (1.0-fd)*en(j,idt)
                    end if
c
                    if (iddfus .eq. 3) then
                      endth = en(j,id)
                      entth = en(j,it)
                    end if
c
                    if (iddfus .eq. 4) then
                      entth = en(j,it)
                      endth = 0.0
                    end if
c
                    if (iddfus .eq. 5) then
                      endth = en(j,id)
                      entth = en(j,it)
                    end if
c
                      beam_thermalddn(j,k)    = endth*sbd*
     .                                        beam_thermalddn_scale(j,k)
                      beam_thermaltth_df(j,k)    = entth*sbd*
     .                                     beam_thermaltth_df_scale(j,k)
                      beam_thermaldth_tf(j,k) = endth*sbd*
     .                                     beam_thermaldth_tf_scale(j,k)
                      beam_thermaltt2n(j,k)   = entth*sbd*
     .                                       beam_thermaltt2n_scale(j,k)
                      beam_thermalddp(j,k)    = endth*sbd*
     .                                        beam_thermalddp_scale(j,k)
            end do
          end do
c
          k                       = k + 1
          beam_thermalddn(j,k)    = 0.0
          beam_thermaltth_df(j,k)    = 0.0
          beam_thermaltt2n(j,k)   = 0.0
          beam_thermaldth_tf(j,k) = 0.0
          beam_thermalddp(j,k)    = 0.0
c
c         sum over beams and energy components and save in last column:
c
          do jj=1,k-1
              beam_thermalddn(j,k)=beam_thermalddn(j,k)
     .                   +beam_thermalddn(j,jj)
              beam_thermalddp(j,k)=beam_thermalddp(j,k)
     .                   +beam_thermalddp(j,jj)
              beam_thermaltth_df(j,k)=beam_thermaltth_df(j,k)
     .                   +beam_thermaltth_df(j,jj)
              beam_thermaldth_tf(j,k)=beam_thermaldth_tf(j,k)
     .                   +beam_thermaldth_tf(j,jj)
              beam_thermaltt2n(j,k)=beam_thermaltt2n(j,k)
     .                   +beam_thermaltt2n(j,jj)
          end do
c         include fast d,thermal t and fast t,thermal d in sbfus:
          sbfus(j)=beam_thermaltth_df(j,k)+beam_thermaldth_tf(j,k)
          qbfus(j) = sbfus(j)*3.5e3                      ! keV/cm**3/sec
      end do                              ! end loop over spatial mesh j
c
      call trapv (r, beam_thermalddn(1,k) , hcap, nj,
     .                               beam_thermal_ddntot  )
      call trapv (r, beam_thermalddp(1,k) , hcap, nj,
     .                               beam_thermal_ddptot  )
      call trapv (r, beam_thermaltth_df(1,k), hcap, nj,
     .                               beam_thermal_dtntot  )
      call trapv (r, beam_thermaltt2n(1,k), hcap, nj,
     .                               beam_thermal_tt2ntot  )
      call trapv (r, beam_thermaldth_tf(1,k), hcap, nj,
     .                               beam_thermaldth_tftot  )
      beam_thermal_ddntot   = beam_thermal_ddntot*volfac ! volfac = 4 *
c                                                          pisq * rmajor
      beam_thermal_ddptot   = beam_thermal_ddptot*volfac
      beam_thermal_dtntot   = beam_thermal_dtntot*volfac
      beam_thermal_tt2ntot  = beam_thermal_tt2ntot*volfac
      beam_thermaldth_tftot = beam_thermaldth_tftot*volfac
c
c     guard against underflow for subsequent 32-bit programs:
c
      if (beam_thermal_ddptot   .lt. 1.0e-30)  beam_thermal_ddptot  =0.0
      if (beam_thermal_dtntot   .lt. 1.0e-30)  beam_thermal_dtntot  =0.0
      if (beam_thermal_tt2ntot  .lt. 1.0e-30)  beam_thermal_tt2ntot =0.0
      if (beam_thermaldth_tftot .lt. 1.0e-30)  beam_thermaldth_tftot=0.0
c
      if (fusionvb .gt. 0) then
         write (*,'(" beam_thermal_ddntot   =", 1pe14.3)')
     .                beam_thermal_ddntot
         write (*,'(" beam_thermal_ddptot   =", 1pe14.3)')
     .                beam_thermal_ddptot
         write (*,'(" beam_thermal_dtntot   =", 1pe14.3)')
     .                beam_thermal_dtntot
         write (*,'(" beam_thermal_tt2ntot  =", 1pe14.3)')
     .                beam_thermal_tt2ntot
         write (*,'(" beam_thermaldth_tftot =", 1pe14.3)')
     .                beam_thermaldth_tftot
      end if
      return
c
      end

      subroutine beam_thermal_fus (time, timmax, bpol, totcur1,
     .                             qbfus, sbfus)
c

c
c -------------------------------------------------- 12/5/95 --- HSJ ---
c     this subroutine calculates beam-thermal fusion rates for the
c     following reactions:             output vector:
c     --------------------            ---------------------
c          Df(Dth,n)He3                  beam_thermalddn
c          Df(Dth,p)T                    beam_thermalddp
c          Df(Tth,n)He4                  beam_thermaltth_df
c          Tf(Dth,n)He4                  beam_thermaldth_tf
c          Tf(Tth,2n)He4                 beam_thermaltt2n
c          Df(he3,p)he4                  beam_thermalhe3dp
c
c     CROSS SECTIONS are from
c     Bosch & Hale, Nuc. Fus., vol32, no.4 (1992) 611
c
c     INPUT:
c     iddfus   =  0    no thermal d,t or dt mixture
c              =  1    thermal species is d
c              =  2    thermal species is dt, use fraction of dt density
c                      given by fd  to specify d density
c              =  3    both d and t are present as separate thermal species
c              =  4    t is a thermal species,d is not present
c              =  5    allows two fluid d,t thermal and single
c                      fluid  'dt' beam
c     fd               fraction of deuterium in dt mixture
c
c     ti(j)           j=1,2..nj, ion temperature (must be in keV)
c     te(j)           j=1,2..nj, electron temperature (must be in keV)
c     en(j,i)         j=1,2..nj,density of ion species i,#/cm**3
c     ene(j)                    density of electrons,#/cm**3
c     zsq(j,i)
c     atw(i)
c     nameb
c     time             current time,sec
c     bpol
c     totcur1           bpol and totcur1 are used to determine correct
c                      signs on some angles
c
c     OUTPUT:
c         if beam d and thermal d are present:
c         beam_thermalddn(j,k)     #/cm**3/sec
c                                  j=1,2,..nj,k=1,2,..3*nbeams
c             beam_thermal_ddntot       volume integrated rate, #/sec
c         if beam d and thermal t are prsent:
c         beam_thermaltth_df(j,k)
c             beam_thermal_dtntot       volume integrated rate, #/sec
c         if beam t and thermal t are present:
c         beam_thermaltt2n(j,k)
c             beam_thermal_tt2ntot       volume integrated rate, #/sec
c         if beam t and thermal d are present
c         beam_thermaldth_tf(j,k)
c             beam_thermaldth_tftot      volume integrated rate, #/sec
c         if beam d and thermal he3 are present
c         beam_thermalhe3th_df(j,k)
c             beam_thermalhe3th_dftot      volume integrated rate, #/sec
c
c     sbfus(j)  beam fusion rate density counts total of
c               all d(t,n)he4 reactions
c     qbfus(j)  alpha power density from sbfus
c
c     beam_thermalddn_scale(j,ksym)
c     beam_thermalddp_scale(j,ksym)
c     beam_thermaltth_df_scale(j,ksym)
c     beam_thermaltt2n_scale(j,ksym)
c     beam_thermaldth_tf_scale(j,ksym)
c     beam_thermalhe3th_df_scale(j,ksym)
c ----------------------------------------------------------------------
c
      USE param
      USE fusion
      USE ions
      USE neut
      USE nub
      USE nub2
      USE soln
      USE numbrs
      USE mesh
      USE verbose
      USE machin
      USE geom
      USE colrate
      USE tordlrot
      USE common_constants,             ONLY: AMU_Value ! in kg
      implicit  integer (i-n), real*8 (a-h, o-z)

c
      real*8  bpol(*), qbfus(*), sbfus(*)
      data    idfdf, idftf, itftf, ixvfast_set /0, 0, 0, 0/
      data    xkeverg, xmassp /6.241e+08, 1.6726e-24/
      data    imod_calc /0/
c
      icxcalc   = icalc_cxfactor     ! load colrate.i
      mass_deut = 3.3435e-24         ! collect into a common area?
      mass_trit = 5.007289e-24
      mass_beam = xmassp*atw_beam
      mass_he3  = AMU_Value*1000.*3.0160293  ! grams
      nterms    = nlegendre  ! terms retained in Legendre expansion of
c                              reaction rate integrals. This is not the
c                              same as # of terms retained in expansion
c                              of the distribution functions !!!
c
      call zeroa (beam_thermalddn,kj*3*kb)
      call zeroa (beam_thermalddp,kj*3*kb)
      call zeroa (beam_thermaltth_df,kj*3*kb)
      call zeroa (beam_thermaltt2n,kj*3*kb)
      call zeroa (beam_thermaldth_tf,kj*3*kb)
      call zeroa (beam_thermalhe3th_df,kj*3*kb)
c
      beam_thermal_long_calc = -1    ! determines if full or approximate
c                                      calcs will be done
      beam_thermal_ddptot     = 0.0
      beam_thermal_dtntot     = 0.0
      beam_thermal_tt2ntot    = 0.0
      beam_thermaldth_tftot   = 0.0
      beam_thermalhe3th_dftot = 0.0
      if (beam_thermal_fusion .le. 0)  return
      IF (iddfus .eq. 0 .AND. idhe3fus .EQ. 0    )  return
      if (beamon(1) .ge. timmax)  return       ! beam is not on in this run
c
      if (beam_thermal_fusion .gt. 0) then ! determine frequncey of calcs
        imod_calc = imod_calc + 1
        if (imod_calc .eq. 1) go to 10 ! always do it on the first call
        if (MOD (imod_calc, beam_thermal_fusion) .eq. 0) go to 10
        if (time .ge. timmax) go to 10 ! always do it at last time point
        print *,'returning from beam_thermal-fusion'
        return     ! do approximate calcs in sub BEAM_THERMAL_APPROX_FUS
      end if
c
c     set up the distribution independent tables,beam_thermal_*,
c     the first time this routine is called:
c     even though there are nbeams, the species in each is the same
c
   10 ebeammax_table = 0.0         ! we will do th efull calcs this time
      do jb=1,nbeams
        ebeammax_table = MAX (ebeammax_table, ebkev(jb))
      end do
      ebeammax_table = 1.05 * ebeammax_table
      ebeammin_table = 0.01        ! keV
      ebeammin_table = 0.001        ! keV HSJ 5/1/01
      ib_d  = 0
      ib_dt = 0
      ib_t  = 0
      if (nameb .eq. 'd' ) ib_d = 1
      if (nameb .eq. 't' ) ib_t = 1
      if (nameb .eq. 'dt') then
        ib_dt = 1
        ib_d  = 1
        ib_t  = 1
      end if
c
      ibeam_chg_exchg = ibcx            ! load colrate.i
      icxcalc         = icalc_cxfactor  ! move it into colrate.i
      ddnt            = 0.0
      vbeammax_table  = SQRT (2.0 * ebeammax_table
     .                / (xkeverg * mass_beam)) ! cm/sec, upper limit
      vbeammin_table =  SQRT (2.0 * ebeammin_table
     .                / (xkeverg * mass_beam)) ! cm/sec, lower limit
      vionmin_table     = 0.0
      vionmax_table = 4.0026e-5 * SQRT (1.5 * ti(1) / xmassp)
c
      do j=1,nj
           z1 = 0.0
           z2 = 0.0
           do i=1,nion
             z12 = en(j,i)*zsq(j,i)
             z1  = z1+z12/atw(i)
             z2  = z2+z12
           end do
           z1 = z1*atw_beam/ene(j)
           z2 = z2/(z1*ene(j)) ! effective charge of bg ions
           den_neut = 0.0
           do i=1,2            ! nneut = 2, hardcoded everywhere
             den_neut = den_neut + enn(j,i)     ! load colrate.i
           end do
c
c          flux surface average cos angle between B and toroidal
c          direction is given by cosvzb. Note that <Bp**2>=
c          Bp0**2*gcap and <Bt**2>=(btor/fcap)**2*<(R0/R)**2>
c
           btavesq = ((btor/fcap(j))**2)*r2cap(j)
           bpavesq = (bpol(j)*gcap(j))**2
           cosvzb  = SQRT (btavesq/(btavesq+bpavesq))
           cosvzb  = cosvzb * btor / ABS (btor)
c
c         let Vf be fast ion velocity,Vr thermal ion velocity. Then
c         e0rel=0.5*mf*(Vf-Vr)dot(Vf-Vr)=ef+(mf/mi)erot
c                         -mag(Vf)*mag(Vr)*mf*cos (theta)
c         But mag(Vf)*cos (theta)*mf = momentum of fast ion in direction
c         of Vr, i.e., in toroidal direction = pb0(j,jc,jb).
c
          tauslocal = taus(j)  ! sec, stored in colrate.i, used in cxint
          erot      = xkeverg*0.5*atw_beam*xmassp*vionz(j)**2
          k         = 0
          do jb=1,nbeams
            iec=3
            if (neg_ion_source(jb) .gt. 0)  iec = 1
            do jc=1,iec
              k     = k+1
c
                e0rel = ABS (ebkev(jb)/jc + erot - xkeverg*pb0(j,jc,jb)
     .                                                   *vionz(j))
                       if (iddfusb_bulk .eq. 1) then
                         v_ion_bulk = vionz(j)*cosvzb  ! bulk v toroid..
c                                                      ..in B direction
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
     .                         * SQRT (2.0*ti(j)/mass_beam)    ! cm/sec
                       vthion  = vbeammin_table                ! CHANGE?
                       vthelec = 1.3256e9
     .                         * SQRT (2.0*te(j))              ! cm/sec
c
c                   sbpure = true (FREYA calc/saved) source rate,#/(cm**3sec):
c                   sbsav has effects of charge exchange and time dependence
c                   already in it, so must use sbpure if icalc_cxfactor = 1
c
                   if (icalc_cxfactor .eq. 0 .and. taupb(j,jc,jb) 
     .                                                 .ne. 0.0) then
                     sbd = enb(j,jc,jb)*taus(j)/taupb(j,jc,jb)
                   else if (icalc_cxfactor .eq. 1) then ! source, sbd,..
c                                                    ..without cx factor
                     sbd = sbpure(j,jc,jb)*taus(j)
                   else              ! source, sbd, has cx factor in it
                     sbd = enb(j,jc,jb)
                   end if
c
                   beam_pitch_angle = zeta(j,jc,jb) * btor * totcur1
     .                                         / ABS (btor * totcur1)
                   vbeam  = SQRT (2.*e0rel/(xkeverg*mass_beam)) ! cm/sec
                   vcrit  = 0.09*((z1/2.)**0.33)*vthelec
                   vcrit3 = vcrit**3
c
                   call set_beam_limits (vbeam,vcrit,taus(j), beamon(1),
     .                                   beam_end(1), time, vthion,
     .                                   vfast_low_lim, vfast_up_lim,
     .                                   icalc_cxfactor)
c
c                set up the sigvrddn and sigvrdtn arrays if not yet set.
c                this table should be applicable for the entire
c                run and thus the range of thermal velocities
c                considered should be large enough to include all
c                possibilities expected. (The table will be
c                recalculated if the range is exceeded so we
c                recover gracefully but waste some computer time)
c
                 tion       = ti(j)                      ! keV
                 expd       = 0.5*mass_beam*xkeverg/tion ! 1/(cm/sec)**2
****             expt       = 0.5*mass_trit*xkeverg/tion ! 1/(cm/sec)**2
                 expt       = expd
                 expdt      = expd
                 exphe3     = 0.5*mass_he3*xkeverg/tion  ! 1/(cm/sec)**2
                 vthuplim   = 2.0 * SQRT (1.0/expd) ! int. limits over..
                 vthlowlim  = 0.0                   ! ..distribution
                 vthuplim   = 2.0 *vthuplim
c
c                if the integration limits are not consistent with
c                the limits in the tables, adjust them here:
c
                 if (vthuplim .gt. vionmax_table) then
                   isetsigvr     = 0
                   betamin       = 0.9 * expd
                   vionmax_table = 2.0 * SQRT (1.0/betamin) ! cm/sec
                   vionmax_table = 1.01 * vthuplim
                 end if
                 if (vbeam .gt. vbeammax_table) then
                   vbeammax_table = 1.01*vbeam              ! cm/sec
                   isetsigvr      = 0
                 end if
c
                 IF (isetsigvr .eq. 0) then
                   isetsigvr = 1
c
c                  load sigvrddn(n2v,n1v)and sigvrddp (stored in colrate.i)
c
                   if (ib_d .ne. 0) then   ! beam has deuterium
c
c                    create table for d(d,n)he3 reaction
c
                     call sigintddn (vionmin_table,vionmax_table,
     .                               vbeammin_table,vbeammax_table)
c
c                    create table for d(d,p)t reaction
c
                     call sigintddp (vionmin_table,vionmax_table,
     .                               vbeammin_table,vbeammax_table)
c                    create table for df(he3,p)he4 reaction
                     IF(idhe3fus .NE. 0)THEN
                        call sigintdhe3 (vionmin_table,vionmax_table,
     .                               vbeammin_table,vbeammax_table)
                     ENDIF
                   ENDIF
c
c                  load sigvrdtn(n2v,n1v) (stored in colrate.i):
c
****               if (ib_t .ne. 0 .or. ib_dt .ne. 0) then  ! beam has..
c                                                           ..tritium
                     if (iddfus .ne. 4)  then      ! must have thermal d
c
c                      create table for d(t,n)he4 reaction:
c
                       call sigintdtn (vionmin_table,vionmax_table,
     .                                 vbeammin_table,vbeammax_table)
                     end if
                     if (iddfus .gt. 1) then      ! must have thermal t
c
c                      create table for t(t,2n)he4 reaction:
c
                       call siginttt2n (vionmin_table,vionmax_table,
     .                                  vbeammin_table,vbeammax_table)
                     end if
****               end if
                 end if
c
c                     do the integrals for all reactions simultaneously:
c
                      call beam_thermal_int1 (iddfus,vthuplim,vthlowlim,
     .                     vbeam,rddn,rdtn,rtt2n,rdth_tf,rddp,rhe3th_df)
c
                      if (iddfus .eq. 1) then
                        endth = en(j,id)
                        entth = 0.0
                      end if
                      if (iddfus .eq. 2) then
                        endth = fd*en(j,idt)
                        entth = (1.0-fd)*en(j,idt)
                      end if
                      if (iddfus .eq. 3) then
                        endth = en(j,id)
                        entth = en(j,it)
                      end if
                      if (iddfus .eq. 4) then
                        entth = en(j,it)
                        endth = 0.0
                      end if
                      if (iddfus .eq. 5) then
                        endth = en(j,id)
                        entth = en(j,it)
                      end if
                      if(idhe3fus .ne. 0)call he3den(enhe3th,j)
                      beam_thermalddn(j,k)    = rddn*endth*sbd*fdbeam
                      beam_thermaltth_df(j,k) = rdtn*entth*sbd*fdbeam
                      beam_thermalhe3th_df(j,k) = rhe3th_df*enhe3th*
     .                                    sbd*fdbeam *d_beam_spin_pol
                      beam_thermaldth_tf(j,k) = rdth_tf*endth
     .                                        * sbd*(1.0-fdbeam)
                      beam_thermaltt2n(j,k)   = rtt2n*entth
     .                                        * sbd*(1.0-fdbeam)
                      beam_thermalddp(j,k)    = rddp*endth*sbd*fdbeam 
            rhe3max = max(rhe3max,rhe3th_df)
            rddnmax = max(rddnmax,rddn)
c
c ---- If we neglect changes in fast and thermal distributions
c ---- and assume that the initial critical velocity doesnt change too
c ---- much then the integrals performed above me be treated as constant,
c ---- So let's save the relevant info:
c
                      beam_thermalddn_scale(j,k)=rddn*fdbeam
                      beam_thermalhe3th_df_scale(j,k)=rhe3th_df*fdbeam
                      beam_thermalddp_scale(j,k)=rddp*fdbeam
                      beam_thermaltth_df_scale(j,k)=rdtn*fdbeam
                      beam_thermaldth_tf_scale(j,k)=rdth_tf*(1.0-fdbeam)
                      beam_thermaltt2n_scale(j,k)=rtt2n*(1.0-fdbeam)
                      beam_thermal_long_calc=1
            end do ! loop over energy components
          end do   ! loop over beams
c
          k                            = k + 1
          beam_thermalddn(j,k)         = 0.0
          beam_thermaltth_df(j,k)      = 0.0
          beam_thermaltt2n(j,k)        = 0.0
          beam_thermaldth_tf(j,k)      = 0.0
          beam_thermalddp(j,k)         = 0.0
          beam_thermalhe3th_df(j,k)    = 0.0
c
c         sum over beams and energy components and save in last column:
c
          do jj=1,k-1
              beam_thermalddn(j,k)=beam_thermalddn(j,k)
     .                   +beam_thermalddn(j,jj)
              beam_thermalddp(j,k)=beam_thermalddp(j,k)
     .                   +beam_thermalddp(j,jj)
              beam_thermaltth_df(j,k)=beam_thermaltth_df(j,k)
     .                   +beam_thermaltth_df(j,jj)
              beam_thermaldth_tf(j,k)=beam_thermaldth_tf(j,k)
     .                   +beam_thermaldth_tf(j,jj)
              beam_thermaltt2n(j,k)=beam_thermaltt2n(j,k)
     .                   +beam_thermaltt2n(j,jj)
              beam_thermalhe3th_df(j,k)=beam_thermalhe3th_df(j,k)
     .                   +beam_thermalhe3th_df(j,jj)
          end do
c         include fast d,thermal t and fast t,thermal d in sbfus,
c         do not include fastd thermal he3
          sbfus(j)=beam_thermaltth_df(j,k)+beam_thermaldth_tf(j,k)
          qbfus(j) = sbfus(j)*3.5e3                      ! keV/cm**3/sec
      end do                              ! end loop over spatial mesh j
c
      call trapv (r, beam_thermalddn(1,k)   , hcap, nj,
     .               beam_thermal_ddntot  )
      call trapv (r, beam_thermalddp(1,k)   , hcap, nj,
     .               beam_thermal_ddptot  )
      call trapv (r, beam_thermaltth_df(1,k)   , hcap, nj,
     .               beam_thermal_dtntot  )
      call trapv (r, beam_thermaltt2n(1,k)  , hcap, nj,
     .               beam_thermal_tt2ntot )
      call trapv (r, beam_thermaldth_tf(1,k), hcap, nj,
     .               beam_thermaldth_tftot)
      call trapv (r, beam_thermalhe3th_df(1,k), hcap, nj,
     .               beam_thermalhe3th_dftot)


      beam_thermal_ddntot     = beam_thermal_ddntot*volfac ! volfac = 4 *
c                                                          pisq * rmajor
      beam_thermal_ddptot     = beam_thermal_ddptot*volfac
      beam_thermal_dtntot     = beam_thermal_dtntot*volfac
      beam_thermal_tt2ntot    = beam_thermal_tt2ntot*volfac
      beam_thermaldth_tftot   = beam_thermaldth_tftot*volfac
      beam_thermalhe3th_dftot = beam_thermalhe3th_dftot*volfac
c
c     guard against underflow for subsequent 32-bit programs:
c
      if (beam_thermal_ddptot   .lt. 1.0e-30)beam_thermal_ddptot  =0.0
      if (beam_thermal_dtntot   .lt. 1.0e-30)beam_thermal_dtntot  =0.0
      if (beam_thermal_tt2ntot  .lt. 1.0e-30)beam_thermal_tt2ntot =0.0
      if (beam_thermaldth_tftot .lt. 1.0e-30)beam_thermaldth_tftot=0.0
      if (beam_thermalhe3th_dftot .lt. 1.0e-30)  
     .                                     beam_thermalhe3th_dftot=0.0
c
      if (fusionvb .gt. 0) then
         write (*, '(" beam_thermal_ddntot   =", 1pe14.3)')
     .                 beam_thermal_ddntot
         write (*, '(" beam_thermal_ddptot   =", 1pe14.3)')
     .                 beam_thermal_ddptot
         write (*, '(" beam_thermal_dtntot   =", 1pe14.3)')
     .                 beam_thermal_dtntot
         write (*, '(" beam_thermal_tt2ntot  =", 1pe14.3)')
     .                 beam_thermal_tt2ntot
         write (*, '(" beam_thermaldth_tftot =", 1pe14.3)')
     .                 beam_thermaldth_tftot
         write (*, '(" beam_thermalhe3th_dftot =", 1pe14.3)')
     .                 beam_thermalhe3th_dftot
 
      
      end if
      return
c
      end

      subroutine beam_thermal_int1 (iddfus, vthuplim, vthlowlim, vbeam,
     .                       rddn, rdtn, rtt2n, rdth_tf,rddp,rhe3th_df)
c
c
c ----------------------------------------------------------------------
c    integrate the thermdist function
c    the weight function,thweight,and evaluation points,thval,were
c    determined on a [0,1] integration interval. Hence they should
c    be used to evaluate integrals like
c        integral from 0 to 1 of { f(x)dx} = sum from 1 to n of {
c                                              thweight(i)*f(thval(i)) }
c
c    Here we must scale
c    the results to the actual interval of integration [0,vthmax]
c
c --- input:
c
c --- output:
c
c        rddn            reactivities rate (cm**3/sec)
c        rdtn
c        rtt2n
c        rdth_tf
c        rddp
c        rhe3th_df
c ------------------------------------------------------------------ HSJ
      USE fusion,                         ONLY : idhe3fus
      USE colrate
c
      implicit  integer (i-n), real*8 (a-h, o-z)
      data idid /0/
c
c     generate the weights,wxvtherml, and evaluation points, xvtherml,
c     for quadrature done below. we use the same set of integration
c     points for both d and t maxwellian distributions (which is obviously
c     not optimal but I have to stop somwehere).
c     do this only once per run:
c
      if (idid .eq. 0) then
         idid = 1
         xx1  = 0.0
         xx2  = 1.0
         call gauleg(xx1,xx2,xvtherml,wxvtherml,nvtherml)
      end if
c
c     evaluate the integral by summing with the appropriate weights:
c
      valddn   = 0.0
      valdtn   = 0.0
      valtt2n  = 0.0
      valddp   = 0.0
      valhe3p  = 0.0
      dv       = vthuplim-vthlowlim
      do j=1,nvtherml
        vth     = xvtherml(j) * dv + vthlowlim ! scale to interval..
c                                              ..pass in colrate.i
        call thermal_distb (iddfus, valddnl, valdtnl, valtt2nl,
     .                      valddpl,valhe3pl)
        valddn   = valddn  + wxvtherml(j) * valddnl
        valdtn   = valdtn  + wxvtherml(j) * valdtnl
        valtt2n  = valtt2n + wxvtherml(j) * valtt2nl
        valddp   = valddp  + wxvtherml(j) * valddpl
        valhe3p  = valhe3p + wxvtherml(j) * valhe3pl
      end do

c
c     done with integrations do some normalization of results:
c
      if (iddfus .eq. 1) then
c
c         for thermal d (no t ,beam must be d):
c
          alphad = (expd/(twopi/2.0))**1.5  ! norm factor for Maxwellian
          rddn   = valddn * dv * twopi * alphad * 1.0e-27
          rddp   = valddp * dv * twopi * alphad * 1.0e-27
      else if (iddfus .eq. 2) then
c
c         for thermal dt mixture (beam must also be dt mixture):
c
          alphadt = (expdt/(twopi/2.))**1.5 ! norm factor for Maxwellian
          rddn    = valddn  * dv * twopi * alphadt * 1.0e-27   ! fast d,
c                                                             thermal d
          rddp    = valddp  * dv * twopi * alphadt * 1.0e-27   ! fast d,
c                                                             thermal d
          rdtn    = valdtn  * dv * twopi * alphadt * 1.0e-27   ! fast d,
c                                                             thermal t
          rtt2n   = valtt2n * dv * twopi * alphadt * 1.0e-27   ! fast t,
c                                                             thermal t
          rdth_tf = rdtn                                       ! fast t,
c                                                             thermal d
      else if (iddfus .eq. 3) then
c
c        for individual thermal d and t species,beam is
c        either d or t (but not a combination of both d and t):
c
         alphad  = (expd/(twopi/2.0))**1.5  ! norm factor for Maxwellian
         alphat  = (expt/(twopi/2.0))**1.5  ! norm factor for Maxwellian
         rddn    = valddn  * dv * twopi * alphad * 1.0e-27
         rddp    = valddp  * dv * twopi * alphad * 1.0e-27
         rdtn    = valdtn  * dv * twopi * alphat * 1.0e-27
         rtt2n   = valtt2n * dv * twopi * alphat * 1.0e-27
         rdth_tf = valdtn  * dv * twopi * alphad * 1.0e-27
      else if (iddfus .eq. 4) then
c
c        for thermal t (no thermal d,beam must be t):
c
         alphat = (expt/(twopi/2.0))**1.5   ! norm factor for Maxwellian
         rtt2n  = valtt2n * dv * twopi * alphat * 1.0e-27
      else if (iddfus .eq. 5) then
c
c        for individual thermal d and t species,beam is dt mixture:
c
         alphad  = (expd/(twopi/2.0))**1.5  ! norm factor for Maxwellian
         alphat  = (expt/(twopi/2.0))**1.5  ! norm factor for Maxwellian
         rddn    = valddn  * dv * twopi * alphad * 1.0e-27
         rddp    = valddp  * dv * twopi * alphad * 1.0e-27
         rdtn    = valdtn  * dv * twopi * alphat * 1.0e-27
         rtt2n   = valtt2n * dv * twopi * alphat * 1.0e-27
         rdth_tf = valdtn  * dv * twopi * alphad * 1.0e-27
      end if
     

      IF(idhe3fus .NE. 0)THEN  
         alphahe3  = (exphe3/(twopi/2.0))**1.5  ! norm factor for Maxwellian
         rhe3th_df = valhe3p * dv * twopi * alphahe3* 1.0e-27
      ENDIF

      return
c
      end

      subroutine beam_thermal_int2 (vthh, valddnl, valdtnl, valtt2nl,
     .                                             valddpl,valhe3pl) 
c
      USE colrate
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c     integrate over fast ion distribution
c ------------------------------------------------------------------ HSJ
c
c      include 'colrate.i'
c
      data idid /0/
c
c     generate the weights, wxvfast, and evaluation points, xvfast,
c     for quadrature done below. do this only once per run:
c
      if (idid .eq. 0) then
         idid = 1
         if (nvfast .eq. nvtherml) then  ! weights = thermal case
            do i=1,nvfast
               xvfast(i) =  xvtherml(i)
               wxvfast(i) = wxvtherml(i)
            end do
         else                            ! get different weights
            xx1 = 0.0
            xx2 = 1.0
            call gauleg (xx1, xx2, xvfast, wxvfast, nvfast)
         end if
      end if
c
c     integrate fast ion distribution function over fast ion speed
c
      valddnl  = 0.0
      valdtnl  = 0.0
      valtt2nl = 0.0
      valddpl  = 0.0
      valhe3pl = 0.0
      dv = vfast_up_lim - vfast_low_lim
      vb = vfast_up_lim          ! set upper integration limit for cxint
      do j=1,nvfast
        vf       = xvfast(j) * dv
     .           + vfast_low_lim ! scale to actual integration limits
        call fast_ion_distb(vf,valddnll,valtt2nll,valdtnll,valddpll)
        valddnl   = valddnl  + wxvfast(j)*valddnll   ! gaussian integration
        valdtnl   = valdtnl  + wxvfast(j)*valdtnll   ! gaussian integration
        valtt2nl  = valtt2nl + wxvfast(j)*valtt2nll  ! gaussian integration
        valddpl   = valddpl + wxvfast(j)*valddpll    ! gaussian integration
        valhe3pl  = valddpl + wxvfast(j)*valdhe3pll  ! gaussian integration
      end do
      valddnl   = valddnl  * dv
      valdtnl   = valdtnl  * dv
      valtt2nl  = valtt2nl * dv
      valddpl   = valddpl  * dv
      valhe3pl  = valhe3pl * dv
      return
c
      end

      real*8 function colf1 (phithh)
c
      USE colrate
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c --- colf1 is the second function to be integrated (over the thermal
c --- ion azimuthal angle phif)
c ----------------------------------------------------------------------
c
c      include 'colrate.i'
c
      phith   = phithh              ! set phith in common block
      a       = colint1 (phithh)
      if (wx .eq. 0.0 .and. wy .eq. 0.0) then
        colf1 = a
      else
c
c       note that betasq = 2.0 * expd * SQRT (1-zetath**2)*vth
c
        colf1 = a * EXP (betasq*(wx * COS (phithh)+wy * SIN (phithh))
     .                 - expd*(wx**2+wy**2))
      end if
      return
c
      end

      real*8 function colf2 (zetaff)
c
      USE colrate
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c --- function is the integrand for colint3
c ------------------------------------------------------------------ HSJ
c
c      include 'colrate.i'
c
      ztf   = zetath*zetaff
      zetaf = zetaff               ! load common block variable
      sqpd  = SQRT ((1.0-zetath**2)*(1.0-zetaf**2))
      answ  = colint2(zetaff)
      colf2 = answ*fastiondist(vf,zetaf)
      return
c
      end

      real*8 function colf3 (vff)
c
      USE colrate
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c --- colf3 is the integrand for colint4
c ----------------------------------------------------------------------
c
c      include 'colrate.i'
c
c     some fast ion distb parms (used in colf2)
c
      vf   = vff
      vf3  = vf**3
      vc3  = vf3 + vcrit3
      vc3i = 1.0 / vc3
      vc4  = ((vb3 + vcrit3) / vc3) * vf3 / vb3
c
      answ = colint3(vf)
      if (ibeam_chg_exchg .gt. 0)
     .  answ = answ * cxint(vf)    ! charge exchange factor cxint
      colf3 = answ*vf*vf
      return
c
      end

      real*8 function colf4 (zetathh)
c
      USE colrate
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c --- function is the integrand for colint5
c ----------------------------------------------------------------------
c
c      include 'colrate.i'
c
      zetath = zetathh
      betasq = 2.0 * expd * SQRT (1.0-zetath*zetath) * vth
      colf4  = colint4(zetathh)
      colf4  = colf4 * EXP (expd*wz*(2.0*vth*zetath -wz))
      return
c
      end

      real*8 function colf5 (vthh)
c
      USE colrate
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c --- colf5 is the integrand for integration over thermal ion speed
c --- called by dfbm2 (through qromb integration routine)
c ------------------------------------------------------------------ HSJ
c
c      include 'colrate.i'
c
      vth   = vthh                       ! load common block
      vthsq = vthh*vthh
      colf5 = vthsq * EXP (-expd*vthsq)
      colf5 = colf5*colint5(vthh)
      return
c
      end

      real*8 function colint1 (phithh)
c
      USE colrate
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c --- subroutine evaluates sigma*v
c ------------------------------------------------------------------ HSJ
c
c      include 'colrate.i'
c
****  sqpd    = SQRT ((1-zetaf*2)*(1-zetath**2))
****  ztf     = zetath*zetaf
****  svthf   = vf**2+vth**2
****  vthf    = 2*vf*vth
c
      cgam    = ztf + sqpd * COS (phith)
      vrel    = SQRT (svthf-vthf*cgam)  ! relative speed
      e       = umdd * vrel * vrel      ! relative energy in c.o.m., keV
      colf00  = ddnhe3(e)
      colf00  = colf00*vrel
      colint1 = colf00
      return
c
      end

      real*8 function colint2 (zetaff)
c
      USE colrate
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c --- subroutine does the integration
c ------------------------------------------------------------------ HSJ
c
c      include 'colrate.i'
c
      real*8   zeroc
      external colf1
c
      zeroc = 0.0
      call qgaus1 (colf1, zeroc, twopi/2.0, answ)
      colint2 = 2.0 * answ
      return
c
      end

      real*8 function colint3 (vff)
c
      USE colrate
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c --- function integrates colf2 over the fast ion pitch angle,zetaf
c ------------------------------------------------------------------ HSJ
c
c      include 'colrate.i'
c
      real*8   one
      external colf2
c
      one = 1.0
      if (vf .gt. nvfact*vcutnterm)
     .  taufctr = acore * LOG ((1.0+vcrit3/vf**3)/(1.0+vcrit3/vb3))
      vthf      = 2.0*vff*vth
      svthf     = vff*vff+vth*vth
      if (vf .lt. vbfrac*vb) then
        call qgaus3 (colf2, -one, one, answ)
        colint3 = answ
      else                  ! add upscattered case here
        zetaf   = zetab
        colint3 = colf2 (zetaf)
      end if
      return
c
      end

      real*8 function colint4 (zetathh)
c
      USE colrate
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c --- function integrates colf3 over fast ion speed:
c ----------------------------------------------------------------------
c
c      include 'colrate.i'
c
      external colf3
c
      call qgaus4 (colf3, vfast_low_lim, vfast_up_lim, answ)
      colint4 = answ
      return
c
      end

      real*8 function colint5 (vthh)
c
      USE colrate
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c --- colint5 integrates colf4 over the  thermal pitch  angle
c ----------------------------------------------------------------------
c
c      include 'colrate.i'
c
      real*8   one
      external colf4
c
      one = 1.0
      call qgaus5 (colf4, -one, one, answ)
      colint5 = answ
      return
c
      end

      real*8 function cxeval (vff)
c
      USE colrate
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c --- evaluates the function which is integrated by cxint
c     cxeval = (v**2/(v**3+vc**3))(sigcx*v)
c --- where sigcx*v is the charge exchange rate in cm**3/sec
c --- sigcx*v is evaluated by function cxr.
c ----------------------------------------------------------------------
c
c      include 'colrate.i'
c
      erel   = 0.5*mass_beam*vff*vff*6.242e8 ! in keV note APPROXIMATION
      sigcx  = cxr(erel)                     ! so debug can look at it
      cxeval = (vff**2/(vff**3+vcrit3))*sigcx
      return
c
      end

      real*8 function cxint (vff)
c
      USE colrate
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c --- function calculates the exponential integrating factor due to
c --- loss of fast ions by charge exchange.
c ----------------------------------------------------------------------
c
c      include 'colrate.i'
c
      external cxeval
c
      call qgaus2 (cxeval, vff, vb, answ)
      answ  = answ * den_neut     ! den_neut is neutral density, #/cm**3
      answ  = tauslocal * answ
      cxint = EXP (-answ)
      return
c
      end

      subroutine ddfbm1 (e0, nd, nddot, rbd)
c
c ----------------------------------------------------------------------
c     this subroutine calculates neutron production due to deuterium
c     neutral beams.  It is a newer (feb 22,'95) version of ddfbm
c     that features more accurate cross sections (taken from
c             Bosch & Hale, Nuc.Fus., vol32, no.4(1992)611 )
c     and allows for inclusion of the high energy fast ion tail,
c             where  v  > vbeam injected
c     This subroutine is called only if the user has selected
c     iddfusb = 1, iddfusb_s = 0
c
c    INPUT PARAMETERS:
c
c        e0    - injection energy of beam ion in keV
c        nd    - thermal deuterium density in cm**-3
c        nddot - number of beam deuterons per sec per cm**-3
c
c        model selection parameters
c
c        iddfusb       = 1, selects this subroutine
c        the following are required if iddfusb = 1:
c           iddfusb_s  = 0, use model which neglects bulk and thermal
c                         motion of ions and charge exchange
c                         effects. This model is similar to the
c                         old model (called when iddfusb = 0, see sub
c                         DDFBM) but uses Bosch & Hale cross
c                         sections and does the integral over
c                         the energy of the fast ion slowing down
c                         distribution numerically (i.e., essentially
c                         exactly,without the approximations involved
c                         in order to get an analytically integrable
c                         result). The result obtained with this
c                         method is correct only if the thermal ions
c                         have negligible motion and hence the
c                         actual neutron rate will tend to be
c                         underestimated at all temperatures.
c
c    OUTPUT PARAMETERS:
c
c     rbd - number of neutrons produced per sec per cm**-3
c
c ------------------------------------------------------ 2/22/95 --- HSJ
c
      USE param
      USE fusion
      USE colrate
      implicit  integer (i-n), real*8 (a-h, o-z)
      real*8  nddot, nd
c
c      include 'param.i'
c      include 'fusion.i'
c      include 'colrate.i'
c
      elowlim = 0.01        ! keV, cross section ~zero below here anyway
      euplim  = e0*ddfusb_t
      if (iddfusb_s .eq. 0) then    ! neglect all motion of thermal ions
            rbd = ddfbm1_s0(elowlim,euplim)
            rbd = rbd*nddot         ! fast ion density factor
            rbd = rbd*nd            ! thermal ion density factor
            if (ddfusb_t .gt. 1.0)  ! include upscaterred fast ions
     .        rbd = rbd+ddfbm1_s0(e0,ddfusb_t*e0)     ! NOT COMPLETE
      else                          ! include motion of thermal ions
              rbd = 0.0             ! just temporarily
      end if
      return
c
      end

      real*8 function ddfbm1_s0a (ekev)
c
      USE colrate
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c --- function is the integrand for ddfbm1_s0
c ------------------------------------------------------------------ HSJ
c
c      include 'colrate.i'               ! pick up ecritlocal
c
      denom      = SQRT (ekev) * (1.0 + (ecritlocal/ekev)**1.5)
      ecom       = ekev*0.5             ! center of mass energy
      a          = ddnhe3(ecom)/denom
      ddfbm1_s0a = a
      return
c
      end

      real*8 function ddfbm1_s0 (elowlim, ebkev)
c
      USE colrate
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c --- function integrates over the fast ion energy.
c --- elowlim,ebkev are limits of integration in keV
c --- this subroutine and ddfbm1_s0a are used when the simple model
c --- of neutron rate calculations that neglects both thermal and
c --- bulk motion of the thermal ions is selected
c --- (iddfusb = 1, iddfusb_s = 0)
c --------------------------------------------------- 2/22/95 ------ HSJ
c
c      include 'colrate.i'
c
      external ddfbm1_s0a
c
c     do the integral over the fast ion energy distribution
c
      call qromb(ddfbm1_s0a,elowlim,ebkev,answ)
      answ      = answ * 1.0e-27   ! cross section in DDNEHE 3 in mbarns
      answ      = answ / (SQRT (2.0 * mass_beam))       ! mass_beam in g
      answ      = answ *  SQRT (1.602e-9)               ! convert to cgs
      ddfbm1_s0 = answ
      return
c
      end

      subroutine ddfbm2 (ti, v_ion_bulk, beam_pitch_angle,
     .                   e0, vbeam, nd, nddot, rbd)
c
      USE colrate
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c     this subroutine calculates neutron production due to deuterium
c     neutral beams.  It is a newer (feb 22,'95) version of ddfbm
c     that features more accurate cross sections (taken from
c             Bosch & Hale, Nuc. Fus., vol 32, no. 4 (1992) 611 )
c     and allows for
c                  a) bulk and thermal motion of thermal ions
c                  b) depletion of the stored fast ion density
c                     by charge exchange with neutrals
c                  c) inclusion of the high energy fast ion tail,
c                     where  v  > vbeam injected
c     This subroutine is called only if the user has selected
c     iddfusb = 1, iddfusb_s = 1
c
c    INPUT PARAMETERS:
c
c        ti            - ion      temperature in keV
c  v_ion_bulk          - bulk ion speed,cm/sec assumed to be in
c                        (+/-) toroidal direction.
c  beam_pitch_angle    - local COS(theta) of beam relative to b field
c        e0            - injection energy of beam ion in keV
c        vbeam         - speed of beam ion in cm/sec
c        nd            - thermal deuterium density in cm**-3
c        nddot         - number of beam deuterons per sec per cm**-3
c        den_neut      - neutral density (for charge exchange losses
c                        of beam)
c        vfast_low_lim -
c        vfast_up_lim  - limits on integral over fast ion distribution
c                        function.
c                        (will vary due to beam turn on and off)
c
c        model selection parameters
c
c        iddfusb         = 1,  together with iddfusb_s=1 selects this
c                              subroutine
c        the following are required to specify the model:
c
c               ddfusb_t = le 1.0, (default = 1.0)
c                         neglect neutrons produced by part of fast
c                         ion distribution function that is above the
c                         energy ddfusb_t*e0,where e0 is the
c                         injected energy normally ddfusb_t = 1.0
c                         should be used. Values less than 1 are not
c                         useful generally.
c                         (NOTE: the model assumes this
c                         is due to collisions with electrons only,the
c                         collision of fast ions amongst themselves,
c                         indpendent of the three injection energies,
c                         is neglected (i.e., the density of fast ions
c                         is supposed to be small for the analytic
c                         distributions used here to make sense).
c               ddfusb_t > 1.0 (user specified)
c                         account for neutrons produced by
c                         part of distribution function between
c                         the injected energy e0 and iddfusb_t*e0.
c                         (normally,due to the rapid decay of the
c                         fast ion distribution function above the
c                         injected energy,
c                         ddfusb_t = 1.25 should suffice)
c
c    OUTPUT PARAMETERS:
c
c     rbd - number of neutrons produced per sec per cm**-3
c
c ------------------------------------------------------ 2/22/95 --- HSJ
c
      real*8   nddot, nd
c
c      include 'colrate.i'
c
      external colf5
c
c     most of following goes into colrate.i
c
      zeroc      = 0.0
      tion      = ti                            ! keV
      expd      = 0.5*mass_beam*6.241e8/tion    ! [1/(cm/sec)**2]
      wflowv    = ABS (v_ion_bulk)              ! cm/sec
      wz        = v_ion_bulk         ! wz is projection onto B direction
      wx        = 0.0
      wy        = 0.0                ! flow perpendicular to B
      vthuplim  = wflowv
     .          + 5.0 * SQRT (1.0/expd)! integration limits over thermal
      vthlowlim = wflowv
     .          - 5.0 * SQRT (1.0/expd)! ion distribution speed
      vthlowlim = MAX (zeroc, vthlowlim)
      alpha     = (expd/(twopi/2.0))**1.5
      zetab     = beam_pitch_angle       ! in [-1,1]
      vb        = vbeam
      vb3       = vb**3
      nterms    = 15                     ! terms to retain in series
      nvfact    =  1                     ! ??
      vbfrac    =  0.99                  ! ??
      acore     = z2*(1.0-zetab**2)/6.0  ! Cores paper, beta=z2/6
****  umdd      = 5.2175e-16             ! (keV/(cm/sec)**2) set in DATA
      abszb     = (1.0 - ABS (zetab)) / (1.0 + ABS (zetab))
c
c     max value of vcutnterm is vb:
c
      vcutnterm = ((1+vcrit3/vb3)* EXP (0.375*z2*abszb)-1.0)**(1.0/3.0)
      vcutnterm = vcrit/vcutnterm  ! determines # terms in Legendre poly
c
c     we limit vcutnterm to vbfrac*vb. Above this value a delta
c     function dependence is assumed.
c     for  vcutnterm .le. vfast .le. vbfrac*vb  a Gaussian dependence
c     is assumed for the fast ion distribution function. This gaussian
c     region will shrink to zero as vcutnterm approaches vbfrac*vb
c
      vcutnterm = MIN (vbfrac*vb, vcutnterm)
c
c     integrate over thermal ion speed (the remainining integrals
c     appear as colf5 calls them)
c
      call qromb (colf5, vthlowlim, vthuplim, answ)
c
      rbd = answ*alpha ! alpha is normalization factor for thermal maxwl
      rbd = rbd*1.0e-27! cross section in sub DDNEHE 3 is in millibarns
      rbd = rbd*nddot  ! fast ion source rate normalization factr
      rbd = rbd*nd     ! thermal deuteron density factor
      return
c
      end

      subroutine ddfbm3 (vthuplim, vthlowlim, vbeam, nd, nddot, rbd)
c
      USE colrate
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c    integrate the thermdist function
c    the weight function,thweight,and evaluation points,thval,were
c    determined on a [0,1] integration interval. Hence they should
c    be used to evaluate integrals like
c        integral from 0 to 1 of { f(x)dx} = sum from 1 to n of {
c                                              thweight(i)*f(thval(i)) }
c
c    Here we must scale
c    the results to the actual interval of integration [0,vthmax]
c
c --- input:
c        nd            - thermal deuterium density in cm**-3
c        nddot         - number of beam deuterons per sec per cm**-3
c
c --- output:
c
c        rbd            neutron rate (#/(cm**3*sec))
c ------------------------------------------------------------------ HSJ
c
c      include 'colrate.i'
c
      real*8  nd, nddot
      data  idid /0/
c
c     generate the weights, wxvtherml, and evaluation points, xvtherml,
c     for quadrature done below. do this only once per run:
c
      if (idid .eq. 0) then
         idid = 1
         xx1  = 0.0
         xx2  = 1.0
         call gauleg(xx1,xx2,xvtherml,wxvtherml,nvtherml)
      end if
c
      alpha = (expd/(twopi/2.0))**1.5       ! norm factor for Maxwellian
      vb    = vbeam
c
c     evaluate the integral by summing with the appropriate weights:
c
      answ = 0.0
      do j=1,nvtherml
        vth  = xvtherml(j) * (vthuplim - vthlowlim)
     .                                 + vthlowlim   ! scale to interval
        answ = answ + wxvtherml(j) * thermdist(vth)
      end do
      answ = answ * (vthuplim - vthlowlim)
      rbd  = answ * twopi * alpha * nddot * nd * 1.0e-27
      return
c
      end

      real*8 function ddnhe3 (e)
c
      implicit none
c
c ----------------------------------------------------------------------
c --- function returns cross section for d(d,n)he3 reaction.
c --- parameterization is taken from
c --- Bosch & Hale, Nuc. Fus., Vol.32, No.4 (1992)
c --- input argument e is energy in keV (in com system)
c --- output is cross section in millibarns (i.e., 10^-27 cm^2)
c --- Results are valid in [0.5,4900] keV. We arbitrarily assume
c --- that extrapolation down to el (=0.1 keV) is valid!
c ------------------------------------------------------------------ HSJ
c
      real*8 a(5), bg, el, eh, s, arg, e, estar, answ
      data   a   , bg, el, eh
     .      /5.3701e+4,  3.3027e2, -1.2706e-1,    2.9327e-5,
     .      -2.5151e-9, 31.3970  ,  0.1      , 4900.0/
c
      estar = MAX (el   , e )
      estar = MIN (estar, eh)
c
****  e must be within range of parmeterization
****  if (e .lt. el .or. e .gt. eh)
**** .  call STOP ('subroutine DDNHE3: unspecified problem', 972)
c
      s      = (((estar*a(5)+a(4))*estar+a(3))*estar+a(2))*estar+a(1)
      arg    = bg / SQRT (estar)
      answ   = s * EXP (-arg) / estar
      ddnhe3 = answ
      return
c
      end

      real*8 function ddpt (e)
c
      implicit none
c
c ----------------------------------------------------------------------
c --- function returns cross section for d(d,p)t reaction.
c --- parameterization is taken from
c --- Bosch & Hale, Nuc. Fus., Vol.32, No.4 (1992)
c --- input argument e is energy in keV (in com system)
c --- output is cross section in millibarns (i.e., 10^-27 cm^2)
c --- Results are valid in [0.5,5000] keV. We arbitrarily assume
c --- that extrapolation down to el (=0.1 keV) is valid!
c ------------------------------------------------------------------ HSJ
c
      real*8 a(5), bg, el, eh, s, arg, e, estar, answ
      data   a   , bg, el, eh
     .      /5.5576e+4,  2.1054e2, -3.2638e-2,    1.4987e-6,
     .       1.8181e-10, 31.3970  ,  0.1      , 5000.0/
c
      estar = MAX (el   , e )
      estar = MIN (estar, eh)
c
****  e must be within range of parmeterization
****  if (e .lt. el .or. e .gt. eh)
**** .  call STOP ('subroutine DDNHE3: unspecified problem', 973)
c
      s    = (((estar*a(5)+a(4))*estar+a(3))*estar+a(2))*estar+a(1)
      arg  = bg / SQRT (estar)
      answ = s * EXP (-arg) / estar
      ddpt = answ
      return
c
      end


      real*8 function dhe3p (e)
c
      implicit none
c
c ----------------------------------------------------------------------
c --- function returns cross section for d(he3,p)he4 reaction.
c --- parameterization is taken from
c --- Bosch & Hale, Nuc. Fus., Vol.32, No.4 (1992)
c --- input argument e is energy in keV (in com system)
c --- output is cross section in millibarns (i.e., 10^-27 cm^2)
c --- Results are valid in [el,eh] keV. We arbitrarily assume
c --- that if e < el then dhe3p =0
c ---------------------------------------------------------HSJ--9/17/13--
c
      real*8 a(5), b(4), bg,el, eh, s, arg, e, estar, answ
      data   a 
     .      /5.7501e6,  2.5226e3,4.5566e1, 0.0,0.0/
      data   b
     .      / -3.1995e-3,-8.5530e-6,5.9014e-8,0.0/
 
      data  bg,el, eh
     .       / 68.7508, 0.3, 900.0/
c
      estar = MAX (el   , e )
      estar = MIN (estar, eh)
c
****  e must be within range of parmeterization
        IF ( e .gt. eh)THEN
           write(6,FMT='("rel energy in dhe3,emin,emax =",
     .       3(x,1pe12.2))') e,el,eh
           call STOP('subroutine Dhe3p: erel out of range(0.3,900kev',
     .             973)
        ENDIF

        IF(e .lt. el)THEN
          dhe3p = 0.0
        ELSE
          s = (((estar*a(5)+a(4))*estar+a(3))*estar+a(2))*estar+a(1)
          s = s/(1. +(((b(3)+estar*b(4))*estar+b(2))*estar+b(1))*estar)
          arg  = bg / SQRT (estar)
          dhe3p = s * EXP (-arg) / estar
        ENDIF
      return
c
      end

      real*8 function dtnhe4 (e)
c
      implicit none
c
c ----------------------------------------------------------------------
c --- function returns cross section for d(t,n)he4 reaction.
c --- parameterization is taken from
c --- Bosch & Hale, Nuc. Fus., Vol.32, No.4 (1992)
c --- input argument e is energy in keV (in com system)
c --- output is cross section in millibarns (i.e., 10^-27 cm^2)
c --- Results are valid in [0.5,550] keV. We arbitrarily assume
c --- that extrapolation down to el (=0.1 keV) is valid!
c ------------------------------------------------------------------ HSJ
c
      real*8 a(5), b(4), bg, el, eh, s, arg, e, estar, answ
      data   a   , bg, el, eh
     .      /6.6927e+4,  7.454e8, 2.050e+6,    5.2002e+4,
     .       0.0, 34.3827  ,  0.1      , 550.0/
      data   b
     .     / 63.8, -.995, 6.981e-5, 1.728e-4 /
c
      estar = MAX (el   , e )
      estar = MIN (estar, eh)
c
****  e must be within range of parmeterization
****  if (e .lt. el .or. e .gt. eh)
**** .  call STOP ('subroutine DTNHE4: unspecified problem', 974)
c
      s      = (((estar*a(5)+a(4))*estar+a(3))*estar+a(2))*estar+a(1)
      s      = s/((((estar*b(4)+b(3))*estar+b(2))*estar+b(1))*estar+1.)
      arg    = bg / SQRT (estar)
      answ   = s * EXP (-arg) / estar
      dtnhe4 = answ
      return
c
      end

      real*8 function fastdist (vff)
c
      USE colrate
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c --- evaluate the fast ion dist. times the integral defined by sigintddn
c ------------------------------------------------------------------ HSJ
c
c      include 'colrate.i'
c
c --- get the value of the integral over angle between the two velocity
c --- vectors by interpolation in the table sigvrddn(vf,vth).
c
      call tableintrp (v1, n1v, vth, jth)
      if (jth .le. 0 .or. jth .ge. n1v) then
        call STOP ('function FASTDIST: unspecified problem #1', 172)
      end if
      call tableintrp (v2, n2v, vf, jf)
      if (jf .le. 0 .or. jf .ge. n2v) then
        call STOP ('function FASTDIST: unspecified problem #2', 176)
      end if
c
c     use quasilinear interpolation
c
      dv1  = vth-v1(jth)
      dv2  = vf-v2(jf)
      a1   = dv1*dv2
      a2   = (v1(jth+1)-vth)*dv2
      a3   = (v1(jth+1)-vth)*(v2(jf+1)-vf)
      a4   = dv1*(v2(jf+1)-vf)
      area = a1+a2+a3+a4
      val  = (a1*sigvrddn(jf+1,jth+1)+a2*sigvrddn(jf+1,jth)
     .     + a3*sigvrddn(jf,jth)+a4*sigvrddn(jf,jth+1))/area
c
c     the charge exchange integral, cxint,
c     is not precomputed because it depends on te (through vcrit).
c     Thus we would need a separate table in each region
c     and would have to recompute anyway each time te changes
c
      vfsq = vf*vf
      val  = val/(vfsq*vf+vcrit3)
      if (icxcalc .eq. 1) then         ! assume cx function of vf
        fastdist = vfsq*val*cxint(vf)  ! source does not include cx part
      else
        fastdist = vfsq*val            ! assume cx independent of vf
      end if                           ! source includes the cx part
      return
c
      end

      real*8 function fastint (vthh)
c
      USE colrate
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c     integrate function fastdist
c ------------------------------------------------------------------ HSJ
c
c      include 'colrate.i'
c
      data idid /0/
c
c     generate the weights, wxvfast, and evaluation points, xvfast,
c     for quadrature done below. do this only once per run:
c
      if (idid .eq. 0) then
         idid = 1
         if (nvfast .eq. nvtherml) then  ! weights = thermal case
            do i=1,nvfast
               xvfast(i) =  xvtherml(i)
              wxvfast(i) = wxvtherml(i)
            end do
         else                            ! get different weights
            xx1 = 0.0
            xx2 = 1.0
            call gauleg (xx1, xx2, xvfast, wxvfast, nvfast)
         end if
      end if
c
c     integrate fast ion distribution function over fast ion speed
c
      answ = 0.0
      do j=1,nvfast
        vf   = xvfast(j) * (vfast_up_lim - vfast_low_lim)
     .        + vfast_low_lim       ! scale to actual integration limits
        answ = answ+wxvfast(j)*fastdist(vf)       ! gaussian integration
      end do
      fastint = answ * (vfast_up_lim - vfast_low_lim)
      return
c
      end

      subroutine fast_ion_distb (vff, valddnll, valtt2nll, valdtnll,
     .                                                     valddpll)
c
      USE colrate
      implicit  integer (i-n), real*8 (a-h, o-z)
      jf=0
c
c ----------------------------------------------------------------------
c --- evaluate the fast ion distribution times the integral defined by
c --- vrel*sigma
c ------------------------------------------------------------------ HSJ
c
c --- get the value of the integral over angle between the two velocity
c --- vectors by interpolation in the table sigvrddn(vf,vth)
c
      call tableintrp (v1, n1v, vth, jth)
      if (jth .le. 0 .or. jth .ge. n1v) then
        print *,'v1(1),v1(n1v) =',v1(1),v1(n1v)
        print *,'vth =',vth
        print *,'vth must be between v1(1) and v1(n1v)'
        call STOP ('subroutine FAST_ION_DISTB: problem #1', 175)
      end if
      if(jf .gt. n2v .or. jf .lt. 1)jf = n2v/2   !HSJ 05/02502
      call tableintrp (v2, n2v, vf, jf)
      if (jf .le. 0 .or. jf .ge. n2v) then
        print *,'v2(1),v2(n2v) =',v2(1),v2(n2v)
        print *,'vf =',vf
        print *,'vf must be between v2(1) and v2(n2v)'
        call STOP ('subroutine FAST_ION_DISTB: problem #2', 202)
      end if
c
c     use quasilinear interpolation
c
      dv1  = vth-v1(jth)
      dv2  = vf-v2(jf)
      a1   = dv1*dv2
      a2   = (v1(jth+1)-vth)*dv2
      a3   = (v1(jth+1)-vth)*(v2(jf+1)-vf)
      a4   = dv1*(v2(jf+1)-vf)
      area = a1+a2+a3+a4
c     if beam has deuterium (then thermal species must also have deuterium)
      if (ib_d .ne. 0) then
       valddnll  = (a1*sigvrddn(jf+1,jth+1)+a2*sigvrddn(jf+1,jth)
     .     + a3*sigvrddn(jf,jth)+a4*sigvrddn(jf,jth+1))/area
       valddpll  = (a1*sigvrddp(jf+1,jth+1)+a2*sigvrddp(jf+1,jth)
     .     + a3*sigvrddp(jf,jth)+a4*sigvrddp(jf,jth+1))/area
      end if
c
c     if beam has tritium (then thermal species must also have tritium)
c
      if (ib_dt .ne. 0)
     . valtt2nll  = (a1*sigvrtt2n(jf+1,jth+1)+a2*sigvrtt2n(jf+1,jth)
     .     + a3*sigvrtt2n(jf,jth)+a4*sigvrtt2n(jf,jth+1))/area
c     if beam has d and or t and thermal has d and or t then need:
      valdtnll  = (a1*sigvrdtn(jf+1,jth+1)+a2*sigvrdtn(jf+1,jth)
     .     + a3*sigvrdtn(jf,jth)+a4*sigvrdtn(jf,jth+1))/area
c
c     the charge exchange integral, cxint,
c     is not precomputed because it depends on te (through vcrit).
c     Thus we would need a separate table in each region
c     and would have to recompute anyway each time te changes
c
      vfsq = vf*vf
      factr=1./(vfsq*vf+vcrit3)
      valddnll=factr*valddnll
      valddpll=factr*valddpll
      valdtnll=factr*valdtnll
      valtt2nll=factr*valtt2nll
      if (icxcalc .eq. 1) then         ! assume cx function of vf
        factr=vfsq*cxint(vf)           ! source does not include cx part
        valddnll=valddnll*factr
        valddpll=valddpll*factr
        valdtnll=valdtnll*factr
        valtt2nll=valtt2nll*factr
      else                             ! source includes the cx part
        valddnll=valddnll*vfsq
        valddpll=valddpll*vfsq
        valdtnll=valdtnll*vfsq
        valtt2nll=valtt2nll*vfsq       ! assume cx independent of vf
      end if
      return
c
      end

      real*8 function fastiondist (vff, zetaff)
c
      USE colrate
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c --- evaluate the fast ion distribution function
c --- zetaff is the cos pitch angle relative to beam direction at
c --- the local (R,phi,Z) point.
c --- vff is the fast ion speed.
c ----------------------------------------------------------------------
c
c      include 'colrate.i'
c
      if (vff .le. nvfact*vcutnterm) then    ! nvfact = 1, .le. required
c
c       evaluate the Legendre Polynomials
c
        sumz = 0.5
        if (nterms .eq. 1)  go to 10
        j    = 1
        gz   = 1.
        hz   = 1.
        gb   = 1.
        hb   = 1.
        fz   = zetaf
        fb   = zetab       ! pitch angle of beam (rel to mag axis)
        sumz = sumz+fz*fb*1.5*vc4**(z2/3.)
        if (nterms .eq. 2)  go to 10
        gz   = fz
        gb   = fb
    5   j    = j+1
        fz   = ((2*j-1)*zetaf*gz-(j-1)*hz)/j
        fb   = ((2*j-1)*zetab*gb-(j-1)*hb)/j
        addterm = (fz*fb*(2*j+1)/2)*(vc4**((j*(j+1))*z2/6.))
        sumz    = sumz+addterm
        if (j .eq. nterms)  go to 10
        hz = gz
        gz = fz
        hb = gb
        gb = fb
        go to 5
   10   fastiondist = vc3i*sumz
c
      else if (vff .le. vbfrac*vb) then
c
c       vff is large enough so that ion distribution is peaked in
c       direction of zetab
c
        fastiondist = vc3i *  EXP (-(zetaf-zetab)**2/(4.0*taufctr))
     .              / (2.0 * SQRT (twopi*taufctr/2.0))
c
      else if (vff .lt. (2.0-vbfrac)*vb) then
c
c       near injected velocity, use delta functions
c
        fastiondist = vc3i
      else                  ! above injected velocity
        continue            ! okay, what goes here?
      end if
      return
c
      end
     

      SUBROUTINE he3den(enhe3th,gp)
!-------------------------------------------------------------------------
! -- return thermal he3 density at grid point gp
!-------------------------------------------------------------------------

      USE fusion,                         ONLY : idhe3fus,ihe,ihei,
     .                                           he3_frac,   
     .                                           he3_thermal_spin_pol
      USE soln ,                          ONLY : en      !en(kj,kion)

      USE numbrs,                         ONLY : nprim,nj


      IMPLICIT NONE
      INTEGER gp
      REAl*8 enhe3th

      IF(idhe3fus ==1)THEN         ! prim d and prim he
            enhe3th  =  en(gp,ihe) 
      ELSEIF(idhe3fus ==2)THEN     ! prim d impurity he
            enhe3th  =  en(gp,nprim+ihei)
      ELSEIF(idhe3fus ==3)THEN     ! prim d,he and impurity he
            enhe3th  =  (en(gp,nprim+ihei) + en(gp,nprim+ihe))
      ELSEIF(idhe3fus ==4)THEN     ! prim dt and prim he
            enhe3th  =  en(gp,ihe)
      ELSEIF(idhe3fus ==5)THEN     ! prim dt and imp he
            enhe3th  =  en(gp,nprim+ihei)
      ELSEIF(idhe3fus ==6)THEN     ! prim dt and he and imp he
            enhe3th  =  (en(gp,nprim+ihei) + en(gp,nprim+ihe))
      ENDIF
      enhe3th  = enhe3th*he3_thermal_spin_pol*he3_frac
      RETURN
      END SUBROUTINE he3den

      subroutine non_inductive_cd (dt_tdem)
c
c ----------------------------------------------------------------------
c    This subroutine determines the resistivity and ohmic
c    current that satisfies Faraday's law.
c    (Subroutine returns with all output zeroed  if tdem mode
c    is not set)
c    The diffeq solved in this subroutine is documented in the
c    section "Faraday's Law In Tdem Mode " in the
c    ....ONETWOPAGES/sources_o12.ps document.   HSJ
c
c  INPUT:
c
c          etap(j)           neoclassical resistivity in units of sec
c
c                   The following input current densities are used
c                   To enable the separation of the product of
c                   resistivity and ohmic current:
c          curpar_soln(j)     <J dot B /Bt0>,parallel current density,
c                             amps/cm**2,on full grid
c          curdri(j)          driven (rf and beam) parallel current
c                             density,amps/cm**2
c          curboot(j)         bootstrap current density
c
c  OUTPUT (to ohml.i):
c
c          eta_par_ohml(j) the experimental resistivity ,ohm-cm
c                          Given the beam,rf,and bootstrap current
c                          profiles (as determined by the input
c                          models selected by the user) get the
c                          resistivity that satisfies Faraday's Law.
c
c               parallel current density profiles in  amps/cm**2
c                    based on various assumptions:
c         ----------------------------------------------------
c          curohmic_ohml(j)                  ohmic current
c                                   (uses whatever model of resistivity
c                                    is currently stored in etap to
c                                    separate the product of
c                                    eta*curden_ohmic)
c          curni_ohml(j)                  noninductive current
c                                   (uses whatever model of resistivity
c                                    is currently stored in etap to
c                                    separate the product of
c                                    eta*curden_ohmic)
c
c          curboot_ohml(j)                bootstrap current
c                                    (uses eta as above and additionally
c                                     uses models of beam and rf driven
c                                     current as stored in curdri to separate
c                                     out the bootstrap current)
c
c          curdrive_ohml(j)               driven current
c                                    (uses eta as above and the model of
c                                     the bootstrap current as given in
c                                     curboot to determine the driven current)
c
c          see ONETWOPAGES reference cited above for exact details
c ----------------------------------------------------- 8/5/96 -- HSJ --
c
      USE param
      USE soln
      USE contour
      USE limiter
      USE mhdpar
      USE tdem
      USE numbrs
      USE mesh
      USE sourc
      USE machin
      USE geom
      USE flags
      USE tordlrot
      USE constnts
      USE mhdcom
      USE tmpcom

      implicit  integer (i-n), real*8 (a-h, o-z)
c

      include 'ohml.i'     ! eta_par_ohml,curni_ohml,..etc.

      include 'storage.i'  ! rdum(kstore),zdum(kstore)
c      include 'tdem.i'     ! drbpdt_tdem (set in rhomsh = (1/FGHr)*dFGHrBp0/dt
c      include 'tmpcom.i'   ! dscrip(kj)
c      include 'tordlrot.i' ! iangrot
c
      dimension rhs_ohml(kjm1),soln_ohml(kj),dbpnodta(kj)
      dimension soln_ohml1(kjm1),soln_ohml2(kj),bpnohma(kj)
      dimension bpoohma(kj)
      data itest    /  0 /    ! light speed
c
      equivalence (rhs_ohml(1),rdum(1))
      equivalence (soln_ohml(1),rdum(kj+1))
      equivalence (soln_ohml1(1),zdum(1))
      equivalence (soln_ohml2(1),zdum(kj+1))
      equivalence (dbpnodta(1),zdum(2*kj+2))
c
      ifd_ohml=itran(nk-iangrot)
      zeroc=0.0
      call multpl1(eta_par_ohml, kj, zeroc)
      call multpl1(curni_ohml,   kj, zeroc)
      call multpl1(curboot_ohml, kj, zeroc)
      call multpl1(curdrive_ohml,kj, zeroc)
      call multpl1(soln_ohml,kj, zeroc)



****  if (ifd_ohml .gt.  0    ) return       ! not analysis mode
      if (mhdmethd .ne. 'tdem') return       ! info not available



      if (intg_tdem(7) .lt. 0 .or.
     .    intg_tdem(8) .lt. 0) return        ! smooth fits not available
      vloop_ohml=dpsidt_tdem*twopi           ! volts
      vloop_ohml=vloop_ohml/300.             ! statvolts
      xnormal=vloop_ohml*cee/(twopi*rmajor)
      soln_ohml(nj)=xnormal                  ! statvolt/sec (peuck)
c     etor is in v/cm,mult by 0.0033333 to convert to statvolts/cm
c     another bc: assume all non inductive current is zero
c     at the plasma edge. then ohmic=total parallel there
c     so soln_ohml(nj)=cee*H*etapar*curpar_soln
c      soln_ohml(nj)=cee*hcap(nj)*etap(nj)*curpar_soln(nj)*2.9979e9
c ------------------try setting central value---------------------------
c
       soln_ohml(1)=cee*hcap(1)*etap(1)*curohm(1)*2.9979e9
       do j=1,nj-1  ! half grid
         if (j .gt. 1) then
         bpnohma(j)  =  (u(nk-iangrot,j)/(fcap(j)*gcap(j)*hcap(j)*r(j))
     .    +u(nk-iangrot,j+1)/(fcap(j+1)*gcap(j+1)*hcap(j+1)*r(j+1)))*0.5
         else
            bpnohma(j) =(u(nk-iangrot,j+1)/
     .                  (fcap(j+1)*gcap(j+1)*hcap(j+1)*r(j+1)))*0.5
         end if
         bpoohma(j)      = (bpo_save_tdem(j) + bpo_save_tdem(j+1))*0.5
         if (dt_tdem .ne. 0.0)
     .    dbpnodta(j)       = (bpnohma(j) - bpoohma(j))/dt_tdem
          dbpdra = (bpo_save_tdem(j+1)-bpo_save_tdem(j))/dr(j)
          dbdta = dbpnodta(j)
         dscripa     = (dscrip(j) + dscrip(j+1))*0.5  ! cm/sec
         rhs_ohml(j) =  dbdta - dscripa * dbpdra      ! gauss/sec
       end do
       curohmic_ohml(1)=curohm(1)
       do j=2,nj
         soln_ohml(j)=soln_ohml(j-1)+dr(j-1)*rhs_ohml(j-1)
         curohmic_ohml(j)=
     .     soln_ohml(j)/(cee*hcap(j)*etap(j))/2.9979e9 ! amps/cm**2
       end do
c
c ------------------ end central value calcs ---------------------------
c
      soln_ohml(nj)=xnormal            ! statvolt/sec (peuck)
      do j=nj-1,1,-1                   ! get rhs_ohml on the half grid
         hcapa          = (hcap(j)+hcap(j+1))*0.5
         fday2d1a       = (fday2d1(j)+fday2d1(j+1))*0.5*f2d1mult
         fday2d2a       = (fday2d2(j)+fday2d2(j+1))*0.5*f2d2mult
         fday2d3a       = (fday2d3(j)+fday2d3(j+1))*0.5*f2d3mult
         drbpdta        = (drbpdt_tdem(j)+drbpdt_tdem(j+1))*0.5
         if (j .gt. 1) then
         bpnohma(j)  =  (u(nk-iangrot,j)/(fcap(j)*gcap(j)*hcap(j)*r(j))
     .    +u(nk-iangrot,j+1)/(fcap(j+1)*gcap(j+1)*hcap(j+1)*r(j+1)))*0.5
         else
            bpnohma(j) =(u(nk-iangrot,j+1)/
     .                  (fcap(j+1)*gcap(j+1)*hcap(j+1)*r(j+1)))*0.5
         end if
         bpoohma(j)    = (bpo_save_tdem(j) + bpo_save_tdem(j+1))*0.5
c
         if (dt_tdem .ne. 0.0)
     .    dbpnodta(j) = (bpnohma(j) - bpoohma(j))/dt_tdem
c
c        here is one way to define the rhs:
         rhs_ohml(j)  =  drbpdta   ! construct source term,gauss/sec
     .                   - (fday2d1a + fday2d2a + fday2d3a )*hcapa*ra(j)
c                            solution on full grid ,(gauss cm)/sec
         dbdta  = (dbdt_tdem(j)+dbdt_tdem(j+1))*5000.0       ! gauss/sec
         dbpdra = (bp0_tdem(j+1)-bp0_tdem(j))*1.e4/dr(j)     ! gauss/cm
         if (itest .eq. 1 .and. dt_tdem .ne. 0.0) then
              dbdta  = dbpnodta(j)
              dbpdra = (bpo_save_tdem(j+1)-bpo_save_tdem(j))/dr(j)
         end if
         dscripa        = (dscrip(j) + dscrip(j+1))*0.5      ! cm/sec
c        here is another way to define the rhs. This way eliminates some
c        of the noise that might be associated with derivatives of FGHr,etc.
         rhs_ohml(j)    =  dbdta - dscripa * dbpdra          ! gauss/sec
         soln_ohml(j)   = -rhs_ohml(j)*dr(j)+soln_ohml(j+1)  ! gauss/sec
         etaa           = (eta(j) + eta(j+1))*0.5
         curohma        = (curohm(j) + curohm(j+1))*0.5
         soln_ohml1(j)  = cee*hcapa*etaa*curohma*2.9975e9
c
      end do
c
c     get c*(d/drho)(H*eta*Johmic):
c
      soln_ohml2(1)    =  0.0
      soln_ohml2(nj)   =  0.0
      do j=2,nj-1
         soln_ohml2(j)=(soln_ohml1(j)-soln_ohml1(j-1))/dr(j-1)
      end do
c
c     soln_ohml is now c*H*eta_par*<JohmicB/Bt0>
c     with the boundary condition
c                    vloop/(2PIR0)=E0=eta_par*H*<JohmicB/Bt0> at rho edge
c     where vloop is the loop voltage on the plasma SURFACE !!!
c
c     Now come various ways to use this experimental knowledge of
c     the Ohmic current:
c
      do j=1,nj
c
c        given the neoclassical resistivity get the ohmic and
c        non inductive current:
c
         if (etap(j) .ne. 0.0) then
           curohmic_ohml(j) =
     .       soln_ohml(j)/(cee*hcap(j)*etap(j))/2.9979e9  ! amps/cm**2
           curni_ohml(j)    = curpar_soln(j) - curohmic_ohml(j)
         else
           curni_ohml(j)    = 0.0
         end if
c
c        given eta and the rf and beam driven current get the
c        bootstrap current that satisfies Faraday's law:
c
         curboot_ohml(j)=curni_ohml(j)-curdri(j)
c
c        given the bootstrap current and eta get the driven current:
         curdrive_ohml(j)=curni_ohml(j)-curboot(j)
c
c        finally given the currents from models,get the resistivity'
c        that satisfies Faraday's Law:
         if (curohm(j) .ne. 0.0) then
           eta_par_ohml(j) = soln_ohml(j)/
     .                      (cee*hcap(j)*curohm(j)*2.9979e9)
         else
           eta_par_ohml(j) = 0.0
         end if
c
      end do
c
c     use dpsidt at constant rho directly:
c
      do j=1,nj
          if (j .eq. 1) then
             etaa=eta(j)    ! eta(j) is on half grid from 1 to nj-1
          else if ( j .lt. nj) then
             etaa=0.5*(eta(j+1)+eta(j))
          else
             etaa=eta(nj-1)
          end if
          etaa = etap(j)    ! full grid value of eta
          eta_par_tdem_ohml(j)=etap(j)*9.0e11
          curohmic_ohml(j)= dpsidt_const_rho_tdem(j)/
     .                     (hcap(j)*rmajor*etaa*9.0e11)
c
c         given the estimated ohmic current the non inductive becomes:
c
         curni_ohml(j)= curpar_soln(j) - curohmic_ohml(j)
c
c        given eta and the rf and beam driven current get the
c        bootstrap current that satisfies Faraday's law:
c
         curboot_ohml(j)=curni_ohml(j)-curdri(j)
c
c        given the bootstrap current and eta get the driven current:
c
         curdrive_ohml(j)=curni_ohml(j)-curboot(j)
c
c        finally given the currents from models,get the resistivity'
c        that satisfies Faraday's Law:
c
         if (curohm(j) .ne. 0.0) then
           eta_par_ohml(j) = dpsidt_const_rho_tdem(j)/
     .                       (curohm(j)*hcap(j)*rmajor)
         else
           eta_par_ohml(j) = 0.0
         end if
      end do
      return
c
      end

      real*8 function p_legendre (l, zeta)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c --- evaluate the Legendre polynomial of (integer) order l at
c --- cos(theta)=zeta
c --- copyleft StJohn numerical factory
c ------------------------------------------------------------------ HSJ
c
c       evaluate the Legendre Polynomial of integer order l:
c
        f = zeta
        g = f
        if (l .eq. 1)  go to 10
        j    = 1
        f    = 1.
        h    = 1.
        if (l .eq. 0)  go to 10
    5   j    = j+1
        f   = ((2*j-1)*zeta*g-(j-1)*h)/j
        if (j .eq. l)  go to 10
        h = g
        g = f
        go to 5
c
c     viola 
c
   10 p_legendre = f
      return
c
      end

      subroutine save_bpo_tdem
c
c ----------------------------------------------------------------------
c
      USE param
      USE solcon
      USE soln
      USE contour
      USE limiter
      USE mhdpar
      USE tdem
      USE numbrs
      USE mesh
      USE geom
      USE tordlrot
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c      include 'param.i'
c      include 'solcon.i'
c      include 'soln.i'
c      include 'geom.i'
c      include 'limiter.i'
c      include 'mesh.i'
c      include 'numbrs.i'
c      include 'mhdpar.i'
c      include 'contour.i'
c      include 'tordlrot.i'
c      include 'tdem.i'  ! most include files above are required to get
c                         the params for tdem.i
c
      bpo_save_tdem(1)=0.0
      do j=2,nj
        bpo_save_tdem(j) = usave(nk-iangrot,j) /
     .                     (fcap(j)*gcap(j)*hcap(j)*r(j))
      end do
      return
c
      end

      subroutine set_beam_limits (vbeam, vcrit, tausj, beamon, beamof,
     .                            time, vthion, vfast_low_lim,
     .                            vfast_up_lim, icalc_cxfactor)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c
c     transient beam effects are done by using the
c     correct support set of the fast ion distribution
c     function (which is time-dependent):
c     beam just turned on?
c
                    tsldown = LOG ((vbeam**3+vcrit**3)/vcrit**3)
                    tsldown = tsldown*tausj/3.0
                    if ((time .ge. beamon) .and.
     .                  (time-beamon) .le. tsldown) then
c
c                        slowing down distribution is only populated
c                        down to this speed
c
                         vfast_low_lim = (vbeam**3+vcrit**3)*
     .                     EXP (-3.*(time-beamon)/tausj)-vcrit**3
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
                    if (time .ge. beamof) then
c
c                      fast ions from vbeam down to this speed have
c                      all slowed down past this speed:
c
                       vfast_up_lim=(vbeam**3+vcrit**3)*
     .                   EXP (-3.*(time-beamof)/tausj)-vcrit**3
                       if (vfast_up_lim .gt. 0.0)
     .                     vfast_up_lim = vfast_up_lim**(1.0/3.0)
                       vfast_up_lim = MAX (vfast_low_lim,vfast_up_lim)
                    else if (time .lt. beamof) then ! beam currently on
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
      return
c
      end

      real*8 function sgv_ddnhe3 (zeta, vb1, vb2)
c
      USE colrate
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c --- return sigma*vrel for ddnhe3 reaction
c --- zeta is angle between vf (fast ion) and vth (thermal ion) vectors
c ----------------------------------------------------------------------
c
c      include 'colrate.i'
c
      vrelsq     = vb1**2 + vb2**2 - 2.0*vb1*vb2*zeta
      ecom       = umdd * vrelsq              ! keV
      sgv_ddnhe3 = ddnhe3 (ecom) * SQRT (vrelsq)
****  sgv_ddnhe3 = 1.0
      return
c
      end

      real*8 function sgv_ddpt (zeta, vb1, vb2)
c
      USE colrate
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c --- return sigma*vrel for ddpt reaction
c --- zeta is angle between vf (fast ion) and vth (thermal ion) vectors
c ----------------------------------------------------------------------
c
c      include 'colrate.i'
c
****  umdd      = 5.2175e-16             ! (keV/(cm/sec)**2) set in DATA
      vrelsq    = vb1**2 + vb2**2 - 2.0*vb1*vb2*zeta
      ecom      = umdd * vrelsq          ! keV
      sgv_ddpt  = ddpt (ecom) * SQRT (vrelsq)
      return
c
      end



      real*8 function sgv_dhe3p (zeta, vb1, vb2)
c
      USE colrate
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c --- return sigma*vrel for ddpt reaction
c --- zeta is angle between vf (fast ion) and vth (thermal ion) vectors
c ----------------------------------------------------------------------
c
c      include 'colrate.i'
c
****  umdd      = 5.2175e-16             ! (keV/(cm/sec)**2) set in DATA
      vrelsq    = vb1**2 + vb2**2 - 2.0*vb1*vb2*zeta
      ecom      = umdd * vrelsq          ! keV
      sgv_dhe3p  = dhe3p (ecom) * SQRT (vrelsq)
      return
c
      end


      real*8 function sgv_dtnhe4 (zeta, vb1, vb2)
c
      USE colrate
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c --- return sigma*vrel for dtnhe4 reaction
c --- zeta is angle between vf (fast ion) and vth (thermal ion) vectors
c ----------------------------------------------------------------------
c
c      include 'colrate.i'                   ! umdt
c
      umdt       = 6.2571e-16               ! (keV/(cm/sec)**2), for d-t
      vrelsq     = vb1**2 + vb2**2 - 2.0*vb1*vb2*zeta
      ecom       = umdt*vrelsq              ! com energy, keV
      sgv_dtnhe4 = dtnhe4 (ecom) * SQRT (vrelsq)
      return
c
      end

      real*8 function sgv_tt2nhe4 (zeta, vb1, vb2)
c
      USE colrate
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c --- return sigma*vrel for tt2nhe4 reaction
c --- zeta is angle between vf (fast ion) and vth (thermal ion) vectors
c ----------------------------------------------------------------------
c
c      include 'colrate.i'                   ! umtt
c
      umtt        = 7.8138e-16              ! (keV/(cm/sec)**2), for t-t
      vrelsq      = vb1**2 + vb2**2 - 2.0*vb1*vb2*zeta
      ecom        = umtt * vrelsq           ! keV
      sgv_tt2nhe4 = tt2nhe4 (ecom) * SQRT (vrelsq)
      return
c
      end

      subroutine sigintddn (v1min, v1max, v2min, v2max)
c
      USE colrate
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ---------------------------------------------------------------------
c --- SIGINTDDN constructs a table of
c ---          Integral from -1 to 1 {sigma(Erel)*vrel}
c ---  the integration is over COS(theta), where theta is the angle
c ---  between the two velocity vectors of the two distributions
c
c --- input
c     v1max
c     v1min
c     v2min      the max and min speeds to be considered for the
c     v2max      two distributions. [ in cm/sec]
c
c --- Parameters defined in colrate.i
c     n1v        #of intervals for grid of first speed
c     n2v        #                         second
c
c --- output
c     sigvrddn(i,j)      i=1,2..n1,j=1,2..n2 the value of the reaction
c                      rate integral[cm/sec]*[units of sigma]
c
c ---------------------------------------------------------------------
c
c      include 'colrate.i'
c
c     first get the weights for the zeta quadrature rule
c     we use nzeta points
c
      xx1 = -1.
      xx2 =  1.
      call gauleg(xx1,xx2,xzeta,wzeta,nzeta)
      dv1=(v1max-v1min)/(n1v-1)
      dv2=(v2max-v2min)/(n2v-1)
      do j=1,n2v
         vth=v1min+(j-1)*dv1
         v1(j)=vth
         do i=1,n1v
             vf=v2min+(i-1)*dv2      ! pass v1=vth, v2=vf with colrate.i
             v2(i)=vf
             sigvrddn(i,j)=0.0
             do k=1,nzeta
                 sigvrddn(i,j)=sigvrddn(i,j)+wzeta(k)
     .                            *sgv_ddnhe3(xzeta(k),vth,vf)
             end do
         end do
      end do
      return
c
      end

      subroutine sigintddp (v1min, v1max, v2min, v2max)
c
      USE colrate
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ---------------------------------------------------------------------
c --- SIGINTDDP constructs a table of
c ---          Integral from -1 to 1 {sigma(Erel)*vrel}
c ---  the integration is over COS(theta), where theta is the angle
c ---  between the two velocity vectors of the two distributions
c
c --- input
c     v1max
c     v1min
c     v2min      the max and min speeds to be considered for the
c     v2max      two distributions. [ in cm/sec]
c
c --- Parameters defined in colrate.i
c     n1v        #of intervals for grid of first speed
c     n2v        #                         second
c
c --- output
c     sigvrddp(i,j)      i=1,2..n1,j=1,2..n2 the value of the reaction
c                      rate integral[cm/sec]*[units of sigma]
c
c ---------------------------------------------------------------------
c
c      include 'colrate.i'
c
c     first get the weights for the zeta quadrature rule
c     we use nzeta points
c
      xx1 = -1.
      xx2 =  1.
      call gauleg(xx1,xx2,xzeta,wzeta,nzeta)
      dv1=(v1max-v1min)/(n1v-1)
      dv2=(v2max-v2min)/(n2v-1)
      do j=1,n2v
         vth=v1min+(j-1)*dv1
         v1(j)=vth
         do i=1,n1v
             vf=v2min+(i-1)*dv2      ! pass v1=vth, v2=vf with colrate.i
             v2(i)=vf
             sigvrddp(i,j)=0.0
             do k=1,nzeta
                 sigvrddp(i,j)=sigvrddp(i,j)+wzeta(k)
     .                            *sgv_ddpt(xzeta(k),vth,vf)
             end do
         end do
      end do
      return
c
      end

      subroutine sigintdtn (v1min, v1max, v2min, v2max)
c
      USE colrate
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ---------------------------------------------------------------------
c --- SIGINTDTN constructs a table of
c ---          Integral from -1 to 1 {sigma(Erel)*vrel}
c ---  the integration is over COS(theta), where theta is the angle
c ---  between the two velocity vectors of the two distributions
c
c --- input
c     v1max
c     v1min
c     v2min      the max and min speeds to be considered for the
c     v2max      two distributions. [ in cm/sec]
c
c --- Parameters defined in colrate.i
c     n1v        #of intervals for grid of first speed
c     n2v        #                         second
c
c --- output
c     sigvrdtn(i,j)      i=1,2..n1,j=1,2..n2 the value of the reaction
c                      rate integral[cm/sec]*[units of sigma]
c
c ---------------------------------------------------------------------
c
c
c     first get the weights for the zeta quadrature rule
c     we use nzeta points
c
      xx1 = -1.
      xx2 =  1.
      call gauleg(xx1,xx2,xzeta,wzeta,nzeta)
      dv1=(v1max-v1min)/(n1v-1)
      dv2=(v2max-v2min)/(n2v-1)
      do j=1,n2v
         vth=v1min+(j-1)*dv1
         v1(j)=vth
         do i=1,n1v
             vf=v2min+(i-1)*dv2      ! pass v1=vth, v2=vf with colrate.i
             v2(i)=vf
             sigvrdtn(i,j)=0.0
             do k=1,nzeta
                 sigvrdtn(i,j)=sigvrdtn(i,j)+wzeta(k)
     .                            *sgv_dtnhe4(xzeta(k),vth,vf)
             end do
         end do
      end do
      return
c
      end

      subroutine siginttt2n (v1min, v1max, v2min, v2max)
c
      USE colrate
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ---------------------------------------------------------------------
c --- SIGINTTT2N constructs a table of
c ---          Integral from -1 to 1 {sigma(Erel)*vrel}
c ---  the integration is over COS(theta), where theta is the angle
c ---  between the two velocity vectors of the two distributions
c
c --- input
c     v1max
c     v1min
c     v2min      the max and min speeds to be considered for the
c     v2max      two distributions. [ in cm/sec]
c
c --- Parameters defined in colrate.i
c     n1v        #of intervals for grid of first speed
c     n2v        #                         second
c
c --- output
c     sigvrtt2n(i,j)      i=1,2..n1,j=1,2..n2 the value of the reaction
c                      rate integral[cm/sec]*[units of sigma]
c
c ---------------------------------------------------------------------
c
c
c     first get the weights for the zeta quadrature rule
c     we use nzeta points
c
      xx1 = -1.
      xx2 =  1.
      call gauleg(xx1,xx2,xzeta,wzeta,nzeta)
      dv1=(v1max-v1min)/(n1v-1)
      dv2=(v2max-v2min)/(n2v-1)
      do j=1,n2v
         vth=v1min+(j-1)*dv1
         v1(j)=vth
         do i=1,n1v
             vf=v2min+(i-1)*dv2      ! pass v1=vth, v2=vf with colrate.i
             v2(i)=vf
             sigvrtt2n(i,j)=0.0
             do k=1,nzeta
                 sigvrtt2n(i,j)=sigvrtt2n(i,j)+wzeta(k)
     .                            *sgv_tt2nhe4(xzeta(k),vth,vf)
             end do
         end do
      end do
      return
c
      end



      subroutine sigintdhe3 (v1min, v1max, v2min, v2max)
c
      USE colrate
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ---------------------------------------------------------------------
c --- SIGINTDHE3 constructs a table of
c ---          Integral from -1 to 1 {sigma(Erel)*vrel}
c ---  the integration is over COS(theta), where theta is the angle
c ---  between the two velocity vectors of the two distributions
c
c --- input
c     v1max
c     v1min
c     v2min      the max and min speeds to be considered for the
c     v2max      two distributions. [ in cm/sec]
c
c --- Parameters defined in colrate.i
c     n1v        #of intervals for grid of first speed
c     n2v        #                         second
c
c --- output
c     sigvrdhe3p(i,j)      i=1,2..n1,j=1,2..n2 the value of the reaction
c                      rate integral[cm/sec]*[units of sigma]
c
c ---------------------------------------------------------------------
c
c     first get the weights for the zeta quadrature rule
c     we use nzeta points
c
      xx1 = -1.
      xx2 =  1.
      call gauleg(xx1,xx2,xzeta,wzeta,nzeta)
      dv1=(v1max-v1min)/(n1v-1)
      dv2=(v2max-v2min)/(n2v-1)
      do j=1,n2v
         vth=v1min+(j-1)*dv1
         v1(j)=vth
         do i=1,n1v
             vf=v2min+(i-1)*dv2      ! pass v1=vth, v2=vf with colrate.i
             v2(i)=vf
             sigvrdhe3p(i,j)=0.0
             do k=1,nzeta
                 sigvrdhe3p(i,j)=sigvrdhe3p(i,j)+wzeta(k)
     .                            *sgv_dhe3p(xzeta(k),vth,vf)
             end do
         end do
      end do
      return
c
      end

      subroutine thermal_distb (iddfus,valddnl, valdtnl, valtt2nl,
     .                                          valddpl,valhe3pl)
c
      USE colrate
      USE fusion ,                                ONLY : idhe3fus
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c     thermal distribution times fastint
c ------------------------------------------------------------------ HSJ
c
c
      vthsq = vth * vth
      call beam_thermal_int2(vth,valddnl,valdtnl,valtt2nl,valddpl,
     .                       valhe3pl) 

c
      if (iddfus .eq. 1) then
c
c         for thermal d (no t ,beam must be d):
c
          thermdist = vthsq * EXP (-expd*vthsq)
          valddnl   = valddnl * thermdist
          valddpl   = valddpl * thermdist
      else if (iddfus .eq. 2) then
c
c         for thermal dt mixture (beam must also be dt mixture):
c
          thermdistdt = vthsq * EXP (-expdt*vthsq)
          valddnl     = valddnl* thermdistdt
          valddpl     = valddpl* thermdistdt
          valtt2nl    = valtt2nl*thermdistdt
          valdtnl     = valdtnl* thermdistdt
      else if (iddfus .eq. 3) then
c
c        for individual thermal d and t species,beam is
c        either d or t (but not a combination of both d and t):
c
         thermdistd  = vthsq * EXP (-expd*vthsq)
         valddnl     = valddnl*thermdistd
         valddpl     = valddpl*thermdistd
         thermdistt  = vthsq * EXP (-expt*vthsq)
         valdtnl     = valdtnl*thermdistt
      else if (iddfus .eq. 4) then
c
c        for thermal t (no thermal d,beam must be t):
c
         thermdistt = vthsq * EXP (-expt*vthsq)
         valtt2nl   = valtt2nl*thermdistt
      else if (iddfus .eq. 5) then
c
c        for individual thermal d and t species,beam is dt mixture:
c
         thermdistd  = vthsq * EXP (-expd*vthsq)
         valddnl     = valddnl*thermdistd
         valddpl     = valddpl*thermdistd
         thermdistt  = vthsq * EXP (-expt*vthsq)
         valdtnl     = valdtnl*thermdistt
         valtt2nl    = valtt2nl*thermdistt
      else
         call STOP ('subroutine THERMAL_DISTB: IDDFUS not set', 205)
      end if
      
      IF(idhe3fus .ne. 0)THEN
         thermdisthe3  = vthsq * EXP (-exphe3*vthsq)
         valhe3pl      = thermdisthe3 * valhe3pl
      ENDIF
      return
c
      end

      real*8 function thermdist (vthh)
c
      USE colrate
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c     thermal distribution times fastint
c ------------------------------------------------------------------ HSJ
c
c
      vthsq     =   vth * vth
      thermdist = vthsq * EXP (-expd*vthsq)*fastint(vth)
      return
c
      end

      real*8 function tt2nhe4 (e)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ------------------------------------------------------------------ HSJ
c
c function returns the t(t,2n)he4 cross section in milli-barns
c This subroutine is a simple linear interpolation of measured? and calc.?
c values from  "Atomic Data for Controled Fusion Research",ORNL-5207,
c E.W. Thomas et al. .NOV. 1979 for energies ge 6kev
c for the lower energies the NRL formulary table was used
c                         NOTE
c      THERE APPEARS TO BE CONSIDERABLE UNCERATAINTY IN T(T,2N)HE4
c      CROSS SECTIONS  SO THESE VALUES MAY BE QUITE CRUDE.
c input
c   e       com energy in keV
c ----------------------------------------------------------------------
c
      parameter (ndt = 16)
      real*8     ekev_table(ndt), rttxsect_table(ndt), e
      integer    ndtlo, ndtp1, i
      data       ndtlo /0/  ! ndtlo should be saved between calls
c                             DATA statement sets it for very first call
c
      data (ekev_table(i),rttxsect_table(i), i=1,ndt)/
     .               0.0,     0.0,         ! 0.0 means don't have values
     .               1.0,     0.0,
     .               2.0,     0.0,
     .               5.0,     0.0,
     .               6.0,     2.1e-3,
     .               8.0,     1.4e-2,
     .              10.0,     4.5e-2,
     .              20.0,     8.0e-1,
     .              40.0,     5.0,
     .              70.0,     13.0,
     .             100.0,     19.0,
     .             200.0,     32.0,
     .             400.0,     50.0,
     .             700.0,     66.0,
     .             1000.0,     97.0,
     .             3000.0,      97.0 /   !fudged 3000,97 is made up HSJ (to accomodate iter)
c
      if (e .le. ekev_table(ndt)) then
c
c       find tikev in table
c
        call tableintrp(ekev_table,ndt,e,ndtlo)
        if (ndtlo .le. 0 .or. ndtlo .gt. ndt-1) then
          print *,'unknown bug  in tt2nhe4'
          call STOP ('function TT2NHE4: bugcheck', 71)
        else
          ndtp1   = ndtlo + 1
          dele    = ekev_table(ndtp1)-ekev_table(ndtlo)
          delr    = rttxsect_table(ndtp1)-rttxsect_table(ndtlo)
          tt2nhe4 = rttxsect_table(ndtlo)
     .            + (delr/dele)*(e-ekev_table(ndtlo))
        end if
      else
         print *,'energy out of range in tt2nhe4'
         print *,'max table value =',ekev_table(ndt)
         print *,'requested value =',e
        call STOP ('function TT2NHE4: energy out of range', 203)
      end if
      return
c
      end
