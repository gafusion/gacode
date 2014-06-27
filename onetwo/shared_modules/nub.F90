      MODULE nub
#ifdef ONETWO
      USE param, only : kb,nap,kj,kz,ke,kion,kcm,kf

      logical      iterate_beam,beam_iteration
      integer      fast_ion_target, neg_ion_source(kb)
      integer, parameter :: ntimbplt = 15
      integer                                                         &
             imaxbeamconvg, naptr, ne_tk, ibeam, ibion,nbion,         &
             inubplt, mf, mfm1, nbeams, npart, npskip, nsourc
      character*8  nbshape(kb), nashape(nap,kb)
!
      real *8,public,save ::                                          &
             anglev(kb), angleh(kb),                                  &
             aheigh(nap, kb), awidth(nap, kb), bcur(kb),              &
             alen(nap, kb), bleni(kb), blenp(kb), bhofset(kb),        &
             bvofset(kb), bptor(kb), bheigh(kb),                      &
             bwidth(kb), bhfoc(kb), bvfoc(kb), bhdiv(kb),             &
             bvdiv(kb), beamon(kb), beam_end(kb), btime(kb),          &
             ebeam(ke,kb), beamoff(kb),beam_cycles(kb),               &
             ebkev(kb), ennub(kj), fbcur(ke, kb), fap(ke,kb),         &
             fwall(ke,kb), hibr(kj,ke,kb), hibrz(kz,ke,kb),           &
             hicm(kj,ke,kb,kcm), hicmz(kz,ke,kb,kcm), ds_tk,          &
             fe_tk, zeta(kj,ke,kb), csgn, hibrs(kj,kb),               &
             ftrapfi(kz,ke,kb), ftrapfit(kj,ke,kb),                   &
             angmpf(kj,ke,kb), angmpz(kz,ke,kb),                      &
             psif(kf), rowpsi(kj), ranseed, relnub, sfrac1(kb),       &
             rpivot(kb), zpivot(kb), pbeam(ke, kb), tenub(kj),        &
             timbplt(ntimbplt), sb(kj,ke,kb), sbcx(kj,2), sbion(kj),  &
             spb(kj,ke,kb), qb(kj,ke,kb), qbbe(kj,ke,kb),             &
             qbbi(kj,ke,kb), qbf(kj,ke,kb), zzi(kz,kion+1),           &
             zne(kz), zni(kz,kion), zte(kz), zangrot(kz),             &
             qbeami_rot(kj), qbeame_rot(kj), psivol(kz),              &
             freyr(kf), fionx, hdepsmth, sbpure(kj,ke,kb),            &
             atw_beam, rtstcx, relaxden, relaxden_err
!
!
!   iterate_beam = .true.  on startup, iterate the beam to
!                          consistency with thermal ion density
!   imaxbeamconvg          max number of iterations of beam allowed
!                          (used only if iterate_beam =.true.)
!
#elif defined NFREYA
        USE nrtype,                 ONLY : DP,I4B
        REAL(DP)  hdepsmth
#endif
        LOGICAL fidiff_on,& ! set fast ion diffusion smoothng in P_Nfreya(only)
            bfr_neutrlz      ! fbcur(ffulla,fhalfa in nubeam) specified 
                            ! before/after neutralizer
       REAL*8 d_beam_spin_pol

       END MODULE nub
 
