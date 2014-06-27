!
!
     MODULE fusion
!--- INCLUDE file fusion.i     
      USE param, only : kj,kb,ksymbp,kbctim
!
      integer beam_beam_fusion, beam_thermal_fusion, no_beam_fusion,&
              beam_beam_long_calc, beam_thermal_long_calc,          &
              iddfusrate, iddfusb, iddfusb_s, iddfusb_bulk,         &
              icalc_cxfactor,                                       &
              iaslow, ifus, id, it, idt, ihe, itfus, ignflg, iddfus,&
              iddcal,idhe3fus,ihei
      real *8 ,public,save ::                                       &
          fd, fencap(kj), ffe(kj), ffi(kj),                         &
          tfusbb, wtifus,pfusdt,pfusdd,pfustt,pfushe3d,             &
          enalp(kj), enasav(kj), enaav(kj), enaav0(kj),             &
          taupa(kj), tauea(kj), walp(kj), wasav(kj),tausa(kj),      &
          waav(kj), balpha(kj), ddfusn(kj), ddntot,he3_frac,        &
          ddfusm(kj), ddntm, ddfusp(kj), ddntp,ddptot,he3dptot,     &
          fdbeam, ddknct, ddknck(kj), ddbmt, ddbeam(kj),            &
          beam_beamddn(kj,ksymbp),                                  &
          beam_beamdtn(kj,ksymbp),                                  &
          beam_beamddp(kj,ksymbp),                                  &
          beam_beam_ddptot,                                         &
          beam_beamtt2n(kj,ksymbp),                                 &
          beam_beamntot(kj),                                        &
          beam_beamnstot,                                           &
          ddpfus(kj), tot_th_fuse(kj), press_alpha(kj),             &
          ddnfus(kj), dtnfus(kj),tt2nfus(kj),                       &
          ttnfus(kj), he3dpfus(kj), beam_thermal_dtntot,            &
          beam_thermal_tt2ntot, beam_thermaldth_tftot,              &
          beam_thermalhe3th_dftot,                                  &
          beam_thermal_ddntot, beam_beam_ddntot, beam_beam_dtntot,  &
          beam_beam_tt2ntot, ddnthm, ecritalpha(kj), ddfusb_t,      &
          thermal_thermal_ddntot, dtntot, ttntot,                   &
          thermal_thermal_dtntot,                                   &
          thermal_thermal_tt2ntot,                                  &
          thermal_thermal_he3dptot,                                 &
          thermal_thermal_ddptot,                                   &
          exptl_neutron_rate(kbctim),                               &
          qdd, qdt, qtt,qhe3d,nalp_thresh,                          &
          beam_beamddn_scale(kj,ksymbp),                            &
          beam_beamdtn_scale(kj,ksymbp),                            &
          beam_beamddp_scale(kj,ksymbp),                            &
          beam_beamtt2n_scale(kj,ksymbp),                           &
          beam_thermalddn_scale(kj,3*kb),                           &
          beam_thermalddp_scale(kj,3*kb),                           &
          beam_thermaltth_df_scale(kj,3*kb),                        &
          beam_thermaltt2n_scale(kj,3*kb),                          &
          beam_thermaldth_tf_scale(kj,3*kb),                        &
          beam_thermalhe3th_df_scale(kj,3*kb+1),                    &
          beam_thermalddn (kj,3*kb+1),                              &
          beam_thermaltth_df(kj,3*kb+1),                            &
          beam_thermalddp (kj,3*kb+1),                              &
          beam_thermalhe3th_df(kj,3*kb+1),                          &
          beam_thermal_ddptot,                                      &
          beam_thermaltt2n(kj,3*kb+1),                              &
          beam_thermaldth_tf(kj,3*kb+1),                            &
          beam_thermalntot(kj),                                     &
!         NOTE: beam_beam_he3dptot is currently a scalar == 0
!         (we assume that beam d and beam he3 do not exist together) HSJ
          beam_beam_he3dptot,he3_thermal_spin_pol
!
      data walp /kj*0.0 /
      data enalp /kj*0.0 /
      data beam_beam_dtntot,beam_beam_ddptot,beam_beamnstot,        &
           beam_thermal_dtntot, beam_thermal_tt2ntot,               &
           beam_thermaldth_tftot,beam_thermal_ddntot,               &
           beam_beam_tt2ntot,                                       &
           thermal_thermal_ddntot,thermal_thermal_dtntot,           &
           thermal_thermal_tt2ntot,thermal_thermal_he3dptot,        &
           thermal_thermal_ddptot                                   &
           /13*0.0d0/
!
! --- ksymbp is ibe(ibe+1)/2 where ibe is the number of beamlines
! --- (given by parameter kb)
! --- times the number of energies per beamline (which is fixed at 3)
! --- so  ksymb = 3*kb*(3*kb+1)/2 = 21  for the current value of kb
!
       END MODULE fusion
 
