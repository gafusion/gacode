
    MODULE set_12_gcnmp_vars

!    ----------------------------------------------------------------------
!    modules:
!    ------------------------------------------------------------HSJ-------

     USE nrtype,                     ONLY : I4B,DP
     USE fdyimp,                     ONLY : dfdt_ncd => dfdt,                  &
                                            dgdt_ncd => dgdt,                  &
                                            dhdt_ncd => dhdt

     USE tension_spline,             ONLY : tspline90,t716_TSPSI,              &
                                            t716_TSVAL1,t716_TSPSS

     USE common_constants,           ONLY : M2cm,T2gauss,im32icm3,cm2m,        &
                                            im22icm2,Tm2gcm,M22cm2,M32cm3,     &
                                            convert,zeroc,g2kg,joupkev,        &
                                            ppmks,ffpmks,pconvert,kgauss2t,    &
                                            kevpjou,M32CM3,izero

     USE transp,                     ONLY : use_nubeam,beam_data,nameb_nubeam, &
                                            enbeam_species_c,nubeam_version

     USE bd_condtn,                  ONLY : en_bc_inpt,ren_bc_inpt,ken_bc_inpt,&
                                            enein,tein,tiin,curdenin,zeffin,   &
                                            rtein,rtiin,renein,rzeffin,        &
                                            rcurdein,fix_edge_ni_bc_inpt,      &
                                            totcur,bctime,dnedt,dnidt,vloop_bc

     USE param,                      ONLY : ksplin
   
     USE mhdcom,                     ONLY : psi
 
     USE mhdpar,                     ONLY : nw,nh,kpsi

     USE tcoef,                      ONLY : vpinch

!    ----------------------------------------------------------------------
!    shared  modules:
!    ----------------------------------------------------------------------
     USE vector_class,               ONLY : delete_vector,new_vector,       &
                                            Vector_mult_real,assign_Vector, &
                                            assignrt_Vector,zero_Vector,    &
                                            copy_Vector, delete_vector_nf,  &
                                            Vector_mult

     USE grid_class,                 ONLY : nj_ncd => nj,psir_grid, &
                                            rho_grid,rho_gridn,     &
                                            xhm2_ncd => xhm2,       &
                                            xi11_ncd => xi11,       &
                                            xi33_ncd => xi33,       &
                                            xips_ncd => xips,       &
                                            xhm20,xi110,            &
                                            xi330,xips0,            &
                                            eps_ncd     => eps,     &
                                            r2cap_ncd   => r2cap,   &
                                            r2capi_ncd  => r2capi,  &
                                            rcap_ncd    => rcap,    &
                                            rho_mhd_gridnpsi

     USE ions_gcnmp,                 ONLY : namep_ncd => namep,     &
                                            namei_ncd => namei,     &
                                            nion_ncd  => nion,      &
                                            namen_ncd => namen,     &
                                            nprim_ncd => nprim,     &
                                            nimp_ncd  => nimp,      &  
                                            nneu_ncd  => nneu,      &
                                            z_ncd     => z,         &
                                            zsq_ncd   => zsq,       &
                                            zeff_ncd => zeff,       &
                                            name_size,              &
                                            atw_ncd => atw,         &
                                            fd_thermal_ncd => fd_thermal, &
                                            nprimp1

     USE neutral_beams,             ONLY  : nameb_ncd => nameb,       &
                                            nbion_ncd => nbion ,      &
                                            fd_beam_ncd =>fd_beam ,   &
                                            wbeam_ncd    => wbeam,    & 
                                            enbeam_ncd => enbeam,     &
                                             enbeam_tot,              &
                                            storqueb_ncd => storqueb, &
                                            bptor,neut_beam_allocate

     USE Plasma_properties ,        ONLY  : dischg,profile,mhd_dat, &
                                            get_values,diffuse,     &
                                            pwrden, wpdot,          &
                                            pellet,prtcl_src,       &
                                            fus_prod,neut_beam

      USE pl_freq,                  ONLY :  common_frequencies

      USE dep_var,                  ONLY :                             &
                                           dnedt_ncd => dnedt,         &
                                           dnidt_ncd => dnidt,         &
                                           dangrotdt_ncd => dangrotdt, &
                                           dp4
                                           ! (enesave is here and in soln)

     USE solcon_gcnmp,              ONLY : time_ncd     => time,    &
                                           eqtime_ncd   => eqtime,  &
                                           tGCNMf,tGCNMs,time_max

     USE source_terms_gcnmp,        ONLY : stsource,                &
                                           scx_ncd    => scx,       &
                                           sion_ncd   => sion,      &
                                         !  srecom_ncd => srecom,    &
                                           sbcx_ncd   => sbcx,      &
                                           sbeame,                  &
                                           dudtsv_ncd => dudtsv,    &
                                           sbeam_ncd  => sbeam
 
     USE bc_values_gcnmp,           ONLY :  totcur_bc,time_bc,      &
                                            fix_edge_te_bc,         &
                                            fix_edge_ti_bc,         &
                                            fix_edge_rot_bc,        &
                                            fix_edge_ni_bc,         &
                                            zeff_bc,ene_bc,         &
                                            te_bc,ti_bc,            &
                                            angrot_bc,en_bc,        &
                                            vloop_bc_ncd => vloop_bc, &
                                            flux_bc,       &
                                            nbctim_ncd => nbctim

     USE fast_ion_data_gcnmp,       ONLY : walp_ncd     => walp,    &
                                           enalp_ncd    => enalp

     USE neutral_data,              ONLY : enn_ncd      => enn,     &
                                           ennw_ncd     => ennw,    &
                                           volsn_ncd    => volsn,   &
                                           ennv_ncd     => ennv

     USE curden_terms,              ONLY : currf_ncd => currf,      &
                                           curbeam,ibcur,irfc

     USE fusion_gcnmp,              ONLY : neutr_ddn_th,neutr_ddn_beam_beam,   & ! neutron rates
                                           neutr_ddn_beam_thermal,             &
                                           neutr_ddn_knock,neutr_ddn_tot,      &
                                           he3dp_fus,ddp_fus,ddn_fus,          &
                                           dtn_fus,tt2n_fus,                   &
                                           he3_frac_ncd   => he3_frac,         &
                                           he3_thermal_spin_pol_ncd  =>        &
                                           he3_thermal_spin_pol

     USE tglfin,                    ONLY : tglf_p_output,tglf_e_output,        &
                                           tglf_m_output

     USE tcoef,                     ONLY : chie_paleo

     USE shot_info,                 ONLY : shot_id

     USE glf23_gcnmp,               ONLY : glf_p_output,glf_e_output,glf_m_output,    &
                                           glf_gamma_net_e_output,glf_anfreq_output,  &
                                           glf_anfreq2_output,glf_anrate_output,      &
                                           glf_anrate2_output,glf_gamma_net_i_output, &
                                           glf_etg_output,glf_eff_chi_output


!    ----------------------------------------------------------------------
     ! ONETWO modules:
!    ----------------------------------------------------------------------
     USE numbrs,                    ONLY :  nprim, nimp, nion,       &
                                            nneu, nk, nkt, nj,nbctim

     USE contour,                   ONLY :  nplasbdry,rcontr,zcontr

     USE nub,                       ONLY :  nbion,sbcx,sbion,nbeams, &
                                            d_beam_spin_pol

     USE nub2,                      ONLY :  enbeam,wbeam

     USE ions,                      ONLY :  fd_thermal,namep,namen,  &
                                            namei,nameb,zeff,z,zsq,  &
                                            atw

     USE fusion,                    ONLY :  enalp,walp,fdbeam,fd,          &
                                            ddfusn,ddbeam,beam_thermalddn, &
                                            ddnthm,ddknct,ddbmt,ddntot,    &
                                            ddpfus,he3dpfus,dtnfus,        &
                                            beam_thermal_ddntot,           &
                                            beam_beamddn,ddknck,ddnfus,    &
                                            beam_beam_ddntot,tt2nfus,      &
                                            thermal_thermal_ddntot,        &
                                            thermal_thermal_ddptot,        &
                                            thermal_thermal_he3dptot,      &
                                            beam_thermalhe3th_df,          &
                                            beam_thermalhe3th_dftot,       &
                                            he3_frac,he3_thermal_spin_pol

     USE psig,                      ONLY :  pindentnpsi,psivolp,     &
                                            rmajavnpsi,              &
                                            triangnpsi,triangnpsi_l, &
                                            grho1npsi,grho2npsi,     &
                                            cxareanpsi,sfareanpsi,   &
                                            rminavnpsi,psival,rho,   &
                                            elongx,ravgi,ravg,fpsi,  &
                                            qpsi,ffppsi,presspsi,    &
                                            pppsi
     USE flxav,                     ONLY :  npsi,xmagn1,ymagn1,rhomax

     USE gpsi,                      ONLY :  xax,yax,p

     USE machin,                    ONLY :  rmajor,kappa,btor,volume,&
                                            flim

     USE mhdbcdtn,                  ONLY :  mhdmultp

     USE constnts,                  ONLY :  psikgaus,psimks  

     USE etc,                       ONLY :  tocur

     USE sourc,                     ONLY :  totohm,totbeam,totrf,    &
                                            totboot,sion,srecom,scx, &
                                            s,dudtsv,curohm,curboot, &
                                            currf,curdbeam,sbeam,    &
                                            dpedtc,dpidtc,qconde,    &
                                            qcondi,qconve,qconvi,    &
                                            qdelt,qbeame,qbeami,qrfe,&
                                            qrfi,qione,qioni,qsawe,  &
                                            qsawi,qcx,qe2d,qi2d,     &
                                            qfuse,qfusi,qbfuse,      &
                                            qbfusi,qmag,qrad,qohm,   &
                                            stfus,sbfus,qexch,sfus,  &
                                            spellet,eta,qbth_intg

     USE extra,                     ONLY :  betap,beta,q,rgeom,      &
                                            btgeom,circum

     USE geom,                      ONLY :  ali,fcap,gcap,hcap,rcap, &
                                            r2cap,r2capi,dfdt,dgdt,  &
                                            dhdt,drhoadt_geom,rcap0, &
                                            drcapdt,drcapidt,r2capi0,&
                                            dr2idt,rcapi,rcap0i,     &
                                            rhoa0,rhoa

     USE rhog,                      ONLY :  psir,pprim,ffprim,press, &
                                            pressb,btotrmaj,bprmaj,  &
                                            curold,ffpold,ppold,rold,&
                                            qpsir

     USE mesh,                      ONLY :  r,roa,fix_edge_te,       &
                                            fix_edge_ti,fix_edge_rot,&
                                            fix_edge_ni,ni_index

     USE soln,                      ONLY :  te,ti,ene,en,curden,etor,&
                                            rbp,curpar_soln

     USE flx,                       ONLY :  flux,anal_eng_flux_e,          &
                                            anal_eng_flux_i,fluxe,fluxi

     USE neut,                      ONLY :  enn,ennw,volsn,ennv

     USE tordlrot,                  ONLY :  angrot,storqueb,angrotin,      &
                                            rangrot,omegale,qomegapi,      &
                                            qangce,sprcxre,spreimpe,       &
                                            sprcxree,fluxangc

     USE tcoef,                     ONLY :  chieinv,chiinv,xkeneo,         &
                                            xkineo,d,xdchitot,ftrap,       &
                                            xnuse,xnus

     USE neo2d,                     ONLY :  xips,xhm2,xi11,xi33,eps
     USE contour,                   ONLY :  nplasbdry,rplasbdry,           &
                                            zplasbdry,rsep,zsep, rplasmin, &
                                            rplasmax, zplasmin, zplasmax
 
     USE io,                        ONLY :  nout,ncrt,nitre

     USE solcon,                    ONLY :  time,timmax
     USE flags,                     ONLY :  inenez,itran
 
     USE mhdcom,                    ONLY :  psiaxis,psibdry,rma,zma
         
     USE pelcom,                    ONLY : nampel

     USE etc,                       ONLY : bp

     USE limiter,                   ONLY : nlimiter,xlimiter,ylimiter

     USE iterdbmd,                  ONLY : statefile_name

     USE mhdgrid,                   ONLY : rmhdgrid,zmhdgrid

     USE P_Nfreya_12_interface,     ONLY : P_Nfreya_read

     USE zen_state,                 ONLY : get_charge_state,en_calc,set_dzdt

     USE rad_loss,                  ONLY : brems_nions

     USE yoka,                      ONLY : ishot

     IMPLICIT NONE
     INTEGER(I4B)  j,jj,ntot,ncd,iendc,iflag,ibctim,knots,ok
     REAL(DP),DIMENSION(:),ALLOCATABLE :: work,rev_psival,rev_generic
     REAL(DP) bpar(4),tmax,cfa,sm,smtol,smlevel
     LOGICAL periodic,uniform,per


    CONTAINS

    SUBROUTINE set_onetwo_vars
!-------------------------------------------------------------------------
! --- Subroutine takes data from statefile( netcdf or text)
! --- and maps it into onetwo variables with appropriate units.
! --- It is assumed that we are reading a restart file and hence the onetwo
! --- thermal species,etc need to be set.
! --- If instead we are reading a statefile generated by P_Nfreya then
! --- we want to update only the Monte Carlo derived beam sources .
!
!
!  not sure how to best use these:
!        nkt       =                              !cts itran = 1 
!        name_size
!        title of dataset (defined locallly in iter_dbase_nc)
!        shot ditto
!        time_ncd
!        eqtime_ncd 
!        vloop_bc ==0 option not implemented
!        fix_edge_te_bc,fix_edge_ti_bc,fix_edge_rot_bc,fix_edge_ni_bc
!        s(j,..) need to extend to nion+4
!        dangrotdt not used ?
!        enpb 

!        xii*0   not used ??
!  since onetwo is driving GCNMP boundary condition values from
!  gcnmp are not passed back to onetwo variables:
!        tGCNMF
!---------------------------------------------------------------

        USE gpsi,                           ONLY : pppsi_eqdsk,presspsi_eqdsk, &
                                                   fpsi_eqdsk,ffppsi_eqdsk,    &
                                                   qpsi_eqdsk,nxeqd_eqdsk,     &
                                                   psival_eqdsk
        REAL(DP) sum,sum1

         IF(nj .NE. nj_ncd)THEN  ! nj is from inone (or default value),nj_ncd is from state file
            WRITE(nout,8)statefile_name(1:LEN_TRIM(statefile_name)),nj_ncd,nj
            WRITE(ncrt,8)statefile_name(1:LEN_TRIM(statefile_name)),nj_ncd,nj
8           FORMAT(2x," ERROR in using statefile ",a,/,  &
                   2X," size of vectors in statefile = ",i5,/, &
                   2X," size of vectors in inone     = ",i5,/, &
                   2x," must be the same ")
            CALL EXIT
         ENDIF



         ene_bc(:) = ene_bc(:)*im32icm3 

 
         IF(P_Nfreya_read)THEN ! anything readin in here must be changed back to units
                               ! units used in statefile output.
 
 
            ! lhs cgs  <== rhs mks
            neut_beam%sb(:,:,:)    = neut_beam%sb(:,:,:) / M32CM3 ! #/(m^3 sec) ==> #/(cm^3 sec)
            neut_beam%qb(:,:,:)    = neut_beam%qb(:,:,:)/M32CM3   ! Watts/m^3 ==> watts/cm^3 
            neut_beam%spb(:,:,:)   = neut_beam%spb(:,:,:)/(M22cm2*g2kg)   !Kg/(m2-s2) ==>g/(cm^2 sec^2)
            neut_beam%spbr(:,:,:)  = neut_beam%spbr(:,:,:)/(g2kg*M2cm)    !Kg/(m-s2)  ==> g/(cm sec2)
            neut_beam%pb0(:,:,:)   = neut_beam%pb0(:,:,:)*M2cm/g2kg    !Kg m/s  ==> g cm/sec

            RETURN  ! 8888999
         ENDIF



        time       = time_ncd                      ! will be reset if statefile is used for startup
        nj         = nj_ncd
        nion       = nion_ncd
        nprim      = nprim_ncd
        nimp       = nimp_ncd
        nneu       = nneu_ncd
        nk         = nion + dp4                    !cts total vars 
        nbion      = nbion_ncd
        nplasbdry  = dischg%nplasbdry 
        fd_thermal = fd_thermal_ncd
        fdbeam     = fd_beam_ncd
        ntot       = nion + dp4
        !mhd_dat%npsi = npsi
        npsi       = mhd_dat%npsi
        namep(1:nprim)  = namep_ncd(1:nprim_ncd)
        namei(1:nimp)   = namei_ncd(1:nimp_ncd)
        namen(1:nneu)   = namen_ncd(1:nneu_ncd)
        nameb_nubeam(1:nbion)  = nameb_ncd(1:nbion_ncd)
        he3_frac =  he3_frac_ncd
        he3_thermal_spin_pol =  he3_thermal_spin_pol_ncd
        IF(nbion == 1)nameb = nameb_ncd(1)

        ishot            = shot_id%shot_nmbr
        rmajor           = dischg%rmajor     * m2cm
        kappa            = dischg%kappa

        pindentnpsi(1)   = dischg%pindento
        psivolp(1)       = dischg%volo/im32icm3 
        volume           = psivolp(1)
        cxareanpsi(1)    = dischg%areao/im22icm2
        triangnpsi(1)    = dischg%deltao
        circum           = dischg%circum  * m2cm



        !plasma_prop.f90 defines npsi 
        !flxav.f90 defines npsi for onetwo (which most likely 
        !is not the same as in gcnmp. hence we need to interpolate):
        psir             = get_values(psir_grid)   
        psir(:)          = psir(:) * psikgaus  ! from mks to cgs
 
        IF(ALLOCATED(work))DEALLOCATE(work)
        ALLOCATE(work(nj_ncd))

 
        IF( npsi == 100000*mhd_dat%npsi)THEN  !disabled
             work           = get_values(dischg%rmajavnpsi )* m2cm
             periodic       = .FALSE. !interpolant is not periodic
             uniform        = .TRUE.  ! tension is uniform
             ncd            = 2       ! # continuous derivatives
             iendc          = 0       ! natural spline bc
            !iendc          = 1       ! give deriv bc at end points
             bpar(1)        = 0.0     ! used only if iendc =1,2 
             bpar(2)        = 0.0     ! used only if iendc =1,2 
            ! NOTE -- TENSION PARAMETER IS OPTIONAL ARGUMENT IN 
            ! t716_TSVAL1 :
             CALL  t716_TSPSI (mhd_dat%npsi,psir,work,ncd,iendc, &
                               periodic,uniform,bpar)
            ! IFLAG = 0 if values of H are to be computed.
            ! IFLAG = 1 if first derivative values are to be computed.
            ! IFLAG = 2 if second derivative values are to be computed.

              iflag = 0
             CALL  t716_TSVAL1 (mhd_dat%npsi,psir,work,IFLAG,npsi,    &
                                psival,rmajavnpsi)

             work           = get_values(dischg%triangnpsi_u )
             CALL  t716_TSPSI (mhd_dat%npsi,psir,work,ncd,iendc, &
                               periodic,uniform,bpar)
             CALL  t716_TSVAL1 (mhd_dat%npsi,psir,work,IFLAG,npsi,    &
                                psival,triangnpsi)

             
             work           = get_values(dischg%triangnpsi_l )
             CALL  t716_TSPSI (mhd_dat%npsi,psir,work,ncd,iendc, &
                               periodic,uniform,bpar)
             CALL  t716_TSVAL1 (mhd_dat%npsi,psir,work,IFLAG,npsi,    &
                                psival,triangnpsi_l)

             
             work           = get_values(dischg%rminavnpsi )* m2cm
             CALL  t716_TSPSI (mhd_dat%npsi,psir,work,ncd,iendc, &
                               periodic,uniform,bpar)
             CALL  t716_TSVAL1 (mhd_dat%npsi,psir,work,IFLAG,npsi,    &
                                psival,rminavnpsi)
             work           = get_values(dischg%psivolpnpsi )/im32icm3
             CALL  t716_TSPSI (mhd_dat%npsi,psir,work,ncd,iendc, &
                               periodic,uniform,bpar)
             CALL  t716_TSVAL1 (mhd_dat%npsi,psir,work,IFLAG,npsi,    &
                                psival,psivolp)
             
             work           = get_values(dischg%elongxnpsi )
             CALL  t716_TSPSI (mhd_dat%npsi,psir,work,ncd,iendc, &
                               periodic,uniform,bpar)
             CALL  t716_TSVAL1 (mhd_dat%npsi,psir,work,IFLAG,npsi,    &
                                psival,elongx)
             
             work           = get_values(dischg%pindentnpsi )
             CALL  t716_TSPSI (mhd_dat%npsi,psir,work,ncd,iendc, &
                               periodic,uniform,bpar)
             CALL  t716_TSVAL1 (mhd_dat%npsi,psir,work,IFLAG,npsi,    &
                                psival,pindentnpsi)
             
             work           = get_values(dischg%sfareanpsi ) *m22cm2
             CALL  t716_TSPSI (mhd_dat%npsi,psir,work,ncd,iendc, &
                               periodic,uniform,bpar)
             CALL  t716_TSVAL1 (mhd_dat%npsi,psir,work,IFLAG,npsi,    &
                                psival,sfareanpsi)
             
             work           = get_values(dischg%cxareanpsi ) *m22cm2
             CALL  t716_TSPSI (mhd_dat%npsi,psir,work,ncd,iendc, &
                               periodic,uniform,bpar)
             CALL  t716_TSVAL1 (mhd_dat%npsi,psir,work,IFLAG,npsi,    &
                                psival,cxareanpsi)  
           
             work           = get_values(dischg%grho1npsi )
             CALL  t716_TSPSI (mhd_dat%npsi,psir,work,ncd,iendc, &
                               periodic,uniform,bpar)
             CALL  t716_TSVAL1 (mhd_dat%npsi,psir,work,IFLAG,npsi,    &
                                psival,grho1npsi)
           
             work           = get_values(dischg%grho2npsi )
             CALL  t716_TSPSI (mhd_dat%npsi,psir,work,ncd,iendc, &
                               periodic,uniform,bpar)
             CALL  t716_TSVAL1 (mhd_dat%npsi,psir,work,IFLAG,npsi,    &
                                psival,grho2npsi)

        ENDIF

        IF(nh .NE. dischg%nz_mhd .OR. nw  .NE. dischg%nr_mhd)THEN
           WRITE(nout,100)nw,nh,dischg%nr_mhd,dischg%nz_mhd
           WRITE(ncrt,100)nw,nh,dischg%nr_mhd,dischg%nz_mhd
           WRITE(nitre,100)nw,nh,dischg%nr_mhd,dischg%nz_mhd
           CALL STOP('set_onetwo_vars',1)
100        FORMAT(2x,'ERROR: Inconsisent,dimensions: onetwo:',2(2x,i5),/,&
           2x,'statefile:',2(2x,i5))
        ENDIF

        cxareanpsi(1:npsi)  =  get_VAlues(dischg%cxareanpsi)
        cxareanpsi(1:npsi)  =  cxareanpsi(1:npsi)*M22cm2
        sfareanpsi(1:npsi)  =  get_VAlues(dischg%sfareanpsi)
        sfareanpsi(1:npsi)  =  sfareanpsi(1:npsi)*M22cm2
        rmajavnpsi(1:npsi)  =  get_VAlues(dischg%rmajavnpsi) *M2cm
        rminavnpsi(1:npsi)  =  get_VAlues(dischg%rminavnpsi) *M2cm
        triangnpsi(1:npsi)  =  get_VAlues(dischg%triangnpsi_u)
        triangnpsi_l(1:npsi)=  get_VAlues(dischg%triangnpsi_l)
        psivolp(1:npsi)     =  get_VAlues(dischg%psivolpnpsi)/im32icm3
        elongx(1:npsi)      =  get_VAlues(dischg%elongxnpsi)
        pindentnpsi(1:npsi) =  get_VAlues(dischg%pindentnpsi)
        grho1npsi(1:npsi)   =  get_VAlues(dischg%grho1npsi)
        grho2npsi(1:npsi)   =  get_VAlues(dischg%grho2npsi)




        rplasbdry              = get_values(dischg%rplasbdry)
        zplasbdry              = get_values(dischg%zplasbdry)

        nplasbdry              = dischg%nplasbdry
        rplasbdry(1:nplasbdry) = rplasbdry(1:nplasbdry)    *M2cm
        zplasbdry(1:nplasbdry) = zplasbdry(1:nplasbdry)    *M2cm
        rcontr(1:nplasbdry)    = rplasbdry(1:nplasbdry)
        zcontr(1:nplasbdry)    = zplasbdry(1:nplasbdry)

        rplasmin               = dischg%rplasmin*M2cm
        rplasmax               = dischg%rplasmax*M2cm
        zplasmin               = dischg%zplasmin*M2cm
        zplasmax               = dischg%zplasmax*M2cm
        rma                    = dischg%rma*M2cm
        zma                    = dischg%zma*M2cm
        xax(1)                 = rma
        yax(1)                 = zma
        xmagn1                 = rma
        ymagn1                 = zma
        xax(2)                 = zeroc
        yax(2)                 = zeroc
        rsep                   = dischg%rsep*M2cm
        zsep                   = dischg%zsep*M2cm

        ! nh = dischg%nr_mhd,   nz = dischg%nz_mhd is automatic
        psi(:,:)         = mhd_dat%psi(:,:)/psimks ! psi is in kguass cm^2
        p(:,:)           = psi(:,:)

        btor             = mhd_dat%btor * T2gauss ! btor is in guass
        flim             = btor*rmajor
        tocur            = mhd_dat%tot_cur
        mhdmultp = 1
        IF (tocur .LT. 0.0) THEN
          tocur = -tocur
          mhdmultp = -1
        END IF
        totcur(1)        = tocur
        totohm           = mhd_dat%totohm_cur
        totboot          = mhd_dat%totboot_cur
        totbeam          = mhd_dat%totbeam_cur
        totrf            = mhd_dat%totrf_cur
        betap            = mhd_dat%betap
        beta             = mhd_dat%beta
        ali              = mhd_dat%ali
        rcap             = get_values(mhd_dat%rcap)  !rcap = < R>,m 8888888888
        rcap(1:nj)       = rcap(1:nj) * M2cm
        rcap0(:)         = rcap(:)
        drcapdt(:)       = zeroc
      
        nxeqd_eqdsk = izero
        IF(ALLOCATED(presspsi_eqdsk))nxeqd_eqdsk = SIZE(presspsi_eqdsk)
        IF(nxeqd_eqdsk .GT. izero .AND. nxeqd_eqdsk  .NE. mhd_dat%npsi) &
           DEALLOCATE(presspsi_eqdsk,qpsi_eqdsk,pppsi_eqdsk,fpsi_eqdsk, &
                ffppsi_eqdsk,psival_eqdsk)

        IF(.NOT. ALLOCATED(presspsi_eqdsk))THEN
           nxeqd_eqdsk = mhd_dat%npsi
           ALLOCATE(presspsi_eqdsk(nxeqd_eqdsk))           
           ALLOCATE(qpsi_eqdsk(nxeqd_eqdsk))
           ALLOCATE(pppsi_eqdsk(nxeqd_eqdsk))
           ALLOCATE(fpsi_eqdsk(nxeqd_eqdsk))
           ALLOCATE(ffppsi_eqdsk(nxeqd_eqdsk))
           ALLOCATE(psival_eqdsk(nxeqd_eqdsk))
        ENDIF

!        CALL reverse(nxeqd_eqdsk,mhd_dat%qpsinpsi%data,qpsi_eqdsk)
         qpsi_eqdsk(:) = mhd_dat%qpsinpsi%data(:)

!        CALL reverse(nxeqd_eqdsk,mhd_dat%pressnpsi%data,presspsi_eqdsk)

        presspsi_eqdsk(:) = mhd_dat%pressnpsi%data(:)
        presspsi_eqdsk(:) = presspsi_eqdsk(:)*10._DP  ! ergs/cm^3 <== nt/m^2
 
!        CALL reverse(nxeqd_eqdsk,mhd_dat%ffprimnpsi%data,ffppsi_eqdsk)
         ffppsi_eqdsk(:) = mhd_dat%ffprimnpsi%data(:)
         ffppsi_eqdsk(:) = ffppsi_eqdsk(:) *1.0E+04_DP ! gauss <==kg/(A sec^2)
 
!        CALL reverse(nxeqd_eqdsk,mhd_dat%pprimnpsi%data,pppsi_eqdsk)
        pppsi_eqdsk(:) =  mhd_dat%pprimnpsi%data(:)
        pppsi_eqdsk = pppsi_eqdsk *1.0E+07_DP !gram/gauss*cm**3*sec**2 <== A/m^3


!        CALL reverse(nxeqd_eqdsk,mhd_dat%fpsinpsi%data,fpsi_eqdsk)
        fpsi_eqdsk(:)          = mhd_dat%fpsinpsi%data(:)
        fpsi_eqdsk(:)          = fpsi_eqdsk(:)/Tm2gcm

!        CALL reverse(nxeqd_eqdsk,mhd_dat%psivalnpsi%data,psival_eqdsk)
        psival_eqdsk(:)          = mhd_dat%psivalnpsi%data(:)
        psival_eqdsk(:)          = psival_eqdsk(:)/psimks
        npsi = mhd_dat%psivalnpsi%size
        IF(npsi .GT. kpsi)CALL STOP("npsi,kpsi error in sub set_onetwo_vars",1)
        psival(1:npsi) = psival_eqdsk(:)
        presspsi(1:npsi) = presspsi_eqdsk(:)
        pppsi(1:npsi) = pppsi_eqdsk(:)
        ffppsi(1:npsi)= ffppsi_eqdsk(:)
        fpsi(1:npsi)  = fpsi_eqdsk(:)
        qpsi(1:npsi)  = qpsi_eqdsk(:)



        r2cap            = get_values(mhd_dat%r2cap) !r2cap  = <R0**2/R**2>
        r2capi           = get_values(mhd_dat%r2capi) 
        r2capi0          = r2capi
        dr2idt(:)        = zeroc

        rcapi(:)         = get_values(mhd_dat%rcapi) ! rcapi = <1/R>
        rcapi(:)         = rcapi(:)/M2cm
        rcap0i           = rcapi
        drcapidt(:)      = zeroc

        curden           = get_values(mhd_dat%curden)
        curden(1:nj)     = curden(1:nj)            *im22icm2
        curold(1:nj)     = curden(1:nj)
        curohm           = get_values(mhd_dat%curohm)
        curohm(1:nj)     = curohm(1:nj)             *im22icm2
        curboot          = get_values(mhd_dat%curboot)
        curboot(1:nj)    = curboot(1:nj)             *im22icm2
        currf(1:nj)      = currf_ncd(1:nj)           *im22icm2
        curpar_soln(1:nj)= get_values(mhd_dat%curpar)*im22icm2

        r                = get_values(rho_grid)
        r(:)             = r(:)               *    M2cm
        rold(:)          = r(:)
        rho              = get_values(rho_mhd_gridnpsi)
        rho(:)           = rho(:)             *M2cm
        rhomax           = rho(1)
        rhoa             = rhomax
        rhoa0            = rhoa
        drhoadt_geom     = zeroc
        roa              = get_values(rho_gridn)   
        fcap             = get_values(mhd_dat%fcap)
        gcap             = get_values(mhd_dat%gcap)
        hcap             = get_values(mhd_dat%hcap)
        eps(1:nj)        = eps_ncd(1:nj)
        dfdt(1:nj)       = dfdt_ncd(1:nj)
        dgdt(1:nj)       = dgdt_ncd(1:nj)
        dhdt(1:nj)       = dhdt_ncd(1:nj)
        q                = get_values(mhd_dat%q_value)  ! may need ABS(q)
        qpsir(:)            = q(:)
        rbp              = get_values(mhd_dat%rbp)      !mhd_dat%rbp in tesla meters
        rbp(1:nj)        = rbp(1:nj)/ Tm2gcm            !rbp(=G*H*r*Bp) in gauss cm
        btotrmaj         = get_values(mhd_dat%btotrmaj)
        btotrmaj(:)      = btotrmaj(:)/kgauss2t
        bprmaj           = get_values(mhd_dat%bprmaj)
        bprmaj(:)        = bprmaj(:)/kgauss2t
        bp               = get_values(mhd_dat%bp)
        bp(:)            = bp (:)/kgauss2t
        ffprim           = get_values(mhd_dat%ffprim)
        ffprim(:)        = ffprim(:)/ffpmks
        ffpold(:)        = ffprim(:)
        pprim            = get_values(mhd_dat%pprim)
        pprim(:)         = pprim(:)/ppmks
        ppold(:)         = pprim(:)

  
        psival           = get_values(mhd_dat%psivalnpsi)
        psival(:)        = psival(:)/psimks
        psiaxis          = mhd_dat%psiaxis/psimks
        psibdry          = mhd_dat%psibdry/psimks
        ravgi            = get_values(mhd_dat%ravginpsi)
        ravgi(:)         = ravgi(:)/M2cm
        ravg             = get_values(mhd_dat%ravgnpsi)
        ravg(:)          = ravg(:)*M2cm
 
        press            = get_values(profile%press)
        walp(:)          = get_values(profile%walp) ! kev/M^3
        walp(:)          = walp(:)*im32icm3         ! kev/cm^3
        press(:)         = press(:)/pconvert
        pressb           = get_values(profile%pressb)
        pressb(:)        = pressb(:)/pconvert
        te               = get_values(profile%te)
        ti               = get_values(profile%ti)
        ene              = get_values(profile%ene)
        etor             = get_values(profile%etor) / M2cm   ! etor in volts/cm
        angrot           = get_values(profile%angrot)
        angrot(1:nj)     = angrot(1:nj)        ! in rad/sec
        ene(:)           = ene(:)              *   im32icm3
             PER         = .FALSE. 
             UNIFORM     = .TRUE.              ! set tension tmax  in true case
             SM          =  0.005
             SMTOL       =  0.01_DP
             smlevel     =  1./(ene(1)**2)     ! weights may want inverse xchi here 
             work(1:nj)  = ene(1:nj)
             tmax        = 25._DP               !tension 0.<= tmax <= 85.
           !  CALL t716_TSPSS(nj,roa,work,PER,UNIFORM,SM,SMTOL,ene,smlevel,tmax)
        zeff(1:nj)       = zeff_ncd(1:nj)
!        dnedt(1:nj)      = dnedt_ncd(1:nj)     *   im32icm3
        dnedt            = get_values(wpdot%dnedt)*im32icm3
        DO jj =1,nion
            dnidt(1:nj,jj) = dudtsv_ncd(jj,1:nj)
!           dnidt(1:nj,jj) = dnidt_ncd(1:nj,jj)     *   im32icm3
        ENDDO
!       dangrotdt(1:nj)  = dangrotdt_ncd(1:nj) ! 888888888888 not defined

        nlimiter         = dischg%nlimiter
        xlimiter         = get_VAlues(dischg%rlimiter)
        ylimiter         = get_VAlues(dischg%zlimiter)
!        xlimiter(:)      = xlimiter(:)*M2cm
!        ylimiter(:)      = ylimiter(:)*M2cm  ! limiter coords are in meters in onetwo ?
        rgeom            = dischg%rgeom*M2cm
        btgeom           = dischg%btgeom*T2gauss
        rmhdgrid         = get_Values(dischg%rmhdgrid)
        zmhdgrid         = get_Values(dischg%zmhdgrid)
        rmhdgrid(:)      = rmhdgrid(:)*M2cm
        zmhdgrid(:)      = zmhdgrid(:)*M2cm


        xi11(1:nj)       = xi11_ncd(1:nj)
        xi33(1:nj)       = xi33_ncd(1:nj)
        xhm2(1:nj)       = xhm2_ncd(1:nj)
        xips(1:nj)       = xips_ncd(1:nj)


        DO j=1,ntot
           DO jj=1,ntot
              IF(j .LE. nion)THEN  ! nprim+nimp densities
                 IF( jj .LE. nion)THEN
                    d(j,jj,1:nj-1)=  diffuse%dcoef(j,jj,1:nj-1)*M22CM2                   ! d is in cm^2/sec
                 ELSE IF( jj == nion+1 .OR. jj == nion+2)THEN
                    d(j,jj,1:nj-1) =  diffuse%dcoef(j,jj,1:nj-1) /M2cm                   ! d is in (1/(kev cm sec))
                 ELSE IF( jj == nion+3)THEN
                   d(j,jj,1:nj-1) =  diffuse%dcoef(j,jj,1:nj-1)/(T2gauss*M22cm2)        ! d is in (1/(gauss cm^2 sec)
                 ELSE IF( jj == nion+4)THEN
                    d(j,jj,1:nj-1) =  diffuse%dcoef(j,jj,1:nj-1)/M2cm                    ! d is in (1/cm)
                 ENDIF
              ELSEIF( j == nion+1 .OR. j == nion+2)THEN ! Te,Ti
                     IF( jj .LE. nion)THEN
                        d(j,jj,1:nj-1) =  diffuse%dcoef(j,jj,1:nj-1)/im22icm2   ! d is in (cm^2 kev/sec)
                     ELSE IF(  jj == nion+1 .OR. jj == nion+2)THEN
                        d(j,jj,1:nj-1) =  diffuse%dcoef(j,jj,1:nj-1) /M2cm  ! d is in (1/(cm sec))
                     ELSE IF( jj == nion+3)THEN
                        d(j,jj,1:nj-1) =  diffuse%dcoef(j,jj,1:nj-1)/(T2gauss*M22cm2) ! d is in (kev gauss/(cm^2 sec))
                     ELSE IF(  jj == nion+4)THEN
                        d(j,jj,1:nj-1) =  diffuse%dcoef(j,jj,1:nj-1)/M2cm ! d is in (kev/cm)
                     ENDIF
              ELSEIF(j == nion+3)THEN ! rbp
                     IF( jj .LE. nion)THEN 
                        d(j,jj,1:nj-1)=  diffuse%dcoef(j,jj,1:nj-1)*T2gauss/(M22cm2*M22cm2)  ! d is in (gauss cm^4/sec)
                     ELSE IF(  jj == nion+1 .OR. jj == nion+2)THEN
                        d(j,jj,1:nj-1)=  diffuse%dcoef(j,jj,1:nj-1)/Tm2gcm ! d is in (gauss cm /(kev sec))
                     ELSE IF( jj == nion+3)THEN
                        d(j,jj,1:nj-1) =  diffuse%dcoef(j,jj,1:nj-1) ! d is in (1/sec)
                     ELSE IF(  jj == nion+4)THEN
                        d(j,jj,1:nj-1) =  diffuse%dcoef(j,jj,1:nj-1)/Tm2gcm ! d is in (gauss cm)
                     ENDIF
              ELSEIF(j == nion+4)THEN ! toroidal rotation
                     IF( jj .LE. nion)THEN
                        d(j,jj,1:nj-1) =  diffuse%dcoef(j,jj,1:nj-1)/(im22icm2*im22icm2)  ! d is in (g cm^4/sec^2)
                     ELSE IF(  jj == nion+1 .OR. jj == nion+2)THEN
                        d(j,jj,1:nj-1) =  diffuse%dcoef(j,jj,1:nj-1)* M2cm/g2kg           ! d is in (g cm/(kev sec^2))
                     ELSE IF( jj == nion+3)THEN
                        d(j,jj,1:nj-1) =  diffuse%dcoef(j,jj,1:nj-1)* T2gauss/g2kg         ! d is in ( g /(gauss sec^2))
                     ELSE IF(  jj == nion+4)THEN
                        d(j,jj,1:nj-1) =  diffuse%dcoef(j,jj,1:nj-1) /(g2kg *M2cm)        ! d is in ( g cm/sec)
                     ENDIF
              ENDIF
           ENDDO
        ENDDO



        xdchitot(1:nj-1) = d(1,1,1:nj-1) ! particle diffusivity for output
                                         ! only first ion is saved 

        chie_paleo(1:nj) = get_values(diffuse%chie_paleo)*M22cm2 !cm**2/sec
        chieinv          = get_values(diffuse%chieinv)
        chieinv(1:nj)    = chieinv(1:nj) *       M22cm2            !cm**2/sec
        chiinv           = get_values(diffuse%chiinv)
        chiinv(1:nj)     = chiinv(1:nj) *       M22cm2             !cm**2/sec
        xkeneo           = get_values(diffuse%xkeneo)
        xkeneo(1:nj)     = xkeneo(1:nj)        / M2cm
        xkineo           = get_values(diffuse%xkineo)
        xkineo(1:nj)     = xkineo(1:nj)        /M2cm               !1/(cm sec)


        dpidtc(:)        = zeroc
        DO j=1,nion
           dpidtc(1:nj)     = dpidtc(1:nj)+get_values(wpdot%dpidt(j))  !wdot, ions, watts/meter**3
        ENDDO
        dpidtc(1:nj)     = dpidtc(1:nj)/(1.5*convert)      !keV/sec/cm*3


        dpedtc           = get_values(wpdot%dpedt)  !wdot, elec, watts/meter**3
        dpedtc(1:nj)     = dpedtc(1:nj)/(1.5*convert)


        qconde           = get_values(pwrden%qconde)
        qconde(1:nj)     = qconde(1:nj)/convert
        qcondi           = get_values(pwrden%qcondi)
        qcondi(1:nj)     = qcondi(1:nj)/convert
        qcondi           = get_values(pwrden%qcondi)
        qcondi(1:nj)     = qcondi(1:nj)/convert
        qconve           = get_values(pwrden%qconve)
        qconve(1:nj)     = qconve(1:nj)/convert
        qconvi           = get_values(pwrden%qconvi)
        qconvi(1:nj)     = qconvi(1:nj)/convert
        qbeame           = get_values(pwrden%qbeame)
        qbeame(1:nj)     = qbeame(1:nj)/convert
        qbeami           = get_values(pwrden%qbeami)
        qbeami(1:nj)     = qbeami(1:nj)/convert
        qdelt(1:nj)      = get_values(pwrden%qdelt)
        qdelt(1:nj)      = qdelt(1:nj)/convert
        qrfe             = get_values(pwrden%qrfe)   
        qrfe(1:nj)       = qrfe(1:nj) /convert
        qrfi             = get_values(pwrden%qrfi)
        qrfi(1:nj)       = qrfi(1:nj)/convert
        qione            = get_values(pwrden%qione) 
        qione(1:nj)      = qione(1:nj)/convert
        qioni            = get_values(pwrden%qioni)
        qioni(1:nj)      = qioni(1:nj)/convert
        qcx              = get_values(pwrden%qcx) 
        qcx(1:nj)        = qcx(1:nj)/convert
        qe2d             = get_values(pwrden%qe2d)
        qe2d(1:nj)       = qe2d(1:nj)/convert
        qi2d             = get_values(pwrden%qi2d)
        qi2d(1:nj)       = qi2d(1:nj)/convert
        qfuse            = get_values(pwrden%qfuse)
        qfuse(1:nj)      = qfuse(1:nj)/convert
        qfusi            = get_values(pwrden%qfusi)
        qfusi(1:nj)      = qfusi(1:nj)/convert
        qbfuse           = get_values(pwrden%qbfuse)
        qbfuse(1:nj)     = qbfuse(1:nj)/convert
        qbfusi           = get_values(pwrden%qbfusi)
        qbfusi(1:nj)     = qbfusi(1:nj)/convert
        qsawe            = get_values(pwrden%qsawe)
        qsawe(1:nj)      = qsawe(1:nj)/convert
        qsawi            = get_values(pwrden%qsawi)
        qsawi(1:nj)      = qsawi(1:nj)/convert
        qmag             = get_values(pwrden%qmag)
        qmag(1:nj)       = qmag(1:nj)/convert
        qohm             = get_values(pwrden%qohm)
        qohm(1:nj)       = qohm(1:nj)/convert
        qrad             = get_values(pwrden%qrad)
        qrad(1:nj)       = qrad(1:nj)/convert

        ddnfus(1:nj)     = ddn_fus(1:nj)*im32icm3   ! convert to 1/cm^3 sec
        ddpfus(1:nj)     = ddp_fus(1:nj)*im32icm3 
        dtnfus(1:nj)     = dtn_fus(1:nj)*im32icm3 
        tt2nfus(1:nj)    = tt2n_fus(1:nj)*im32icm3 
        he3dpfus(1:nj)   = he3dp_fus(1:nj)*im32icm3 



        storqueb(1:nj)  = storqueb_ncd(1:nj)/pconvert
        ! storqueb in cgs, storqueb_ncd in mks
        stfus(1:nj)     = get_values(prtcl_src%stfuse)/im32icm3
        sbfus(1:nj)     = get_values(prtcl_src%sbfuse)/im32icm3
        sfus(1:nj)      = stfus(1:nj) + sbfus(1:nj)
        spellet(1:nj)   = get_values(prtcl_src%spellet)/im32icm3
        qexch           = get_values(pwrden%qexch)
        qexch(:)        = qexch(1:)*convert
  
        omegale(1:nj)   = get_values(pwrden%omegale)   * kevpjou*im32icm3
        qangce(1:nj)    = get_values(pwrden%qangce)    * kevpjou*im32icm3
        qomegapi(1:nj)  = get_values(pwrden%qomegapi)  * kevpjou*im32icm3
        sprcxre(1:nj)   = get_values(pwrden%sprcxre)   * kevpjou*im32icm3
        spreimpe(1:nj)  = get_values(pwrden%spreimpe)  * kevpjou*im32icm3
        sprcxree(1:nj)  = get_values(pwrden%sprcxree)  * kevpjou*im32icm3


        DO jj =1,nj
           IF(roa(jj)-fix_edge_ni(1) .GE. 0.0_DP)THEN
              j=jj -1
              EXIT
           ENDIF
        ENDDO
 

        DO jj = 1,nion
            en(1:nj,jj)     = get_values(profile%en(jj))
            en(1:nj,jj)     = en(:,jj)           *   im32icm3
!             PER         = .FALSE.
!             UNIFORM     = .TRUE.              ! set tension tmax  in true case
!             SM          =  0.001
!             SMTOL       =  0.01_DP
!             smlevel     =  1./(en(1,jj)**2)                ! weights may want inverse xchi here 
!             work(1:nj)  = en(1:nj,jj)
!             tmax        = 0.0_DP               !tension 0.<= tmax <= 85.
!             CALL t716_TSPSS(j,roa,work,PER,UNIFORM,SM,SMTOL,en(:,jj),smlevel,tmax)
!             CALL  t716_TSVAL1 (nj,roa,en(:,jj),IFLAG,nj,roa,work)
        ENDDO

        DO jj = 1 ,ntot
            DO j=1,nj-1
               flux(jj,j)       = profile%flux(jj)%data(j)
               flux(jj,j)       = flux(jj,j)*im22icm2 
               dudtsv(jj,j)     = dudtsv_ncd(jj,j)
            ENDDO
            dudtsv(jj,nj)       = dudtsv_ncd(jj,nj)
        ENDDO


        work(1:nj)      = get_values(profile%fluxe)
        DO j=1,nj-1
           fluxe(j)     = 0.5_DP*(work(j)+work(j+1))/M22cm2
        ENDDO
        work(1:nj)      = get_values(profile%fluxi)
        DO j=1,nj-1
           fluxi(j)     = 0.5_DP*(work(j)+work(j+1))/M22cm2
        ENDDO

       DO  j=1,nj-1
           sum = zeroc ; sum1 = zeroc
           DO jj=1,nion_ncd
              sum  = sum  + profile%flux_conv(jj)%data(j)/profile%en(jj)%data(j)
              sum1 = sum1 + profile%flux_conv(jj)%data(j+1)/profile%en(jj)%data(j+1)
           ENDDO
           vpinch(j) = 0.5_DP*(sum+sum1)*M2cm/nion_ncd   ! average of all species pinches
           fluxangc(j) = profile%flux_conv(nion+4)%data(j)
       ENDDO
 


        DO jj= 1,MIN(2,nion) ! in onetwo only 2
           sion(1:nj,jj)    = sion_ncd(1:nj,jj)   *   im32icm3
           !srecom(1:nj,jj) = srecom_ncd(1:nj,jj) *   im32icm3
           srecom(1:nj,jj)  = prtcl_src%srecom(jj)%data(1:nj) *im32icm3
           scx(1:nj,jj)     = scx_ncd(1:nj,jj)    *   im32icm3
           sbcx(1:nj,jj)    = sbcx_ncd(1:nj,jj)   *   im32icm3
        ENDDO

        DO jj=1,nion 
           s(jj,1:nj)       = stsource(jj,1:nj)   
           z(1:nj,jj)       = z_ncd(1:nj,jj)
           zsq(1:nj,jj)     = zsq_ncd(1:nj,jj)
           ! Note in onetwo sion(1:kk=9,1;kj)  88888888888888
        ENDDO
        !Note: enbeam_species_c might not be allocated if nbeam is not selected
        IF(.NOT. ALLOCATED(enbeam_species_c))ALLOCATE(enbeam_species_c(nj,nbion))
        DO jj=1,nbion 
          enbeam_species_c(1:nj,jj) = enbeam_ncd(1:nj,jj)  * im32icm3 !one species case
          sbeam(1:nj)        = sbeam_ncd(1:nj,jj)* im32icm3   !beam ion source multipl species 
        ENDDO
 
        enbeam(1:nj)        = enbeam_ncd(1:nj,1)  * im32icm3 !one species case
        sbion(1:nj)         = sbeame(1:nj)        * im32icm3 !beam electron source
        curdbeam(1:nj)      = curbeam(1:nj)       * im22icm2
        wbeam(1:nj)         = wbeam_ncd(1:nj)     * im32icm3
        walp(1:nj)          = walp_ncd(1:nj)      * im32icm3
        enalp(1:nj)         = enalp_ncd(1:nj)     * im32icm3

       DO jj=1,nneu
          enn(1:nj,jj)         = enn_ncd(:,jj)       * im32icm3
          ennv(1:nj,jj)        = ennv_ncd(:,jj)      * im32icm3
          ennw(1:nj,jj)        = ennw_ncd(:,jj)      * im32icm3
          volsn(1:nj,jj)       = volsn_ncd(:,jj)     * im32icm3
       ENDDO

        CALL  check_statefile_equilibrium(xax(1),yax(1),psiaxis)
        rma = xax(1) ; zma = yax(1)
        xmagn1 = rma ; ymagn1 = zma
 
!       release storage used for statefile reads:
        CALL delete_Vector (psir_grid)
        CALL delete_Vector (rho_grid)
        CALL delete_Vector (rho_gridn)
        CALL delete_Vector (mhd_dat%hcap)
        CALL delete_Vector (mhd_dat%gcap)
        CALL delete_Vector (mhd_dat%fcap)
        CALL delete_Vector (mhd_dat%q_value)
        CALL delete_Vector (mhd_dat%curden)
        CALL delete_Vector (mhd_dat%curpar)
        CALL delete_Vector (mhd_dat%curohm)
        CALL delete_Vector (mhd_dat%curboot)
        CALL delete_Vector (dischg%rmajavnpsi)
        CALL delete_Vector (dischg%rminavnpsi)
        CALL delete_Vector (profile%te)
        CALL delete_Vector (profile%ti)
        CALL delete_Vector (profile%ene)
        CALL delete_Vector (profile%etor)
        CALL delete_Vector (profile%angrot)
        CALL delete_Vector (diffuse%chieinv)
        CALL delete_Vector (diffuse%chiinv)
        CALL delete_Vector (diffuse%xkeneo)
        CALL delete_Vector (diffuse%xkineo)
        DO j=1,nion
           CALL delete_Vector (wpdot%dpidt(j))
        ENDDO
        CALL delete_Vector (wpdot%dpedt)
        CALL delete_Vector (wpdot%dnedt)
        CALL delete_Vector (pwrden%qconde)
        CALL delete_Vector (pwrden%qcondi)
        CALL delete_Vector (pwrden%qconve)
        CALL delete_Vector (pwrden%qconvi)
        CALL delete_Vector (pwrden%qbeame)
        CALL delete_Vector (pwrden%qbeami)
        CALL delete_Vector (pwrden%qdelt)
        CALL delete_Vector (pwrden%qrfe)
        CALL delete_Vector (pwrden%qrfi)
        CALL delete_Vector (pwrden%qione)
        CALL delete_Vector (pwrden%qioni)
        CALL delete_Vector (pwrden%qe2d)
        CALL delete_Vector (pwrden%qi2d)
        CALL delete_Vector (pwrden%qcx)
        CALL delete_Vector (pwrden%qfuse)
        CALL delete_Vector (pwrden%qfusi)
        CALL delete_Vector (pwrden%qbfuse)
        CALL delete_Vector (pwrden%qbfusi)
        CALL delete_Vector (pwrden%qmag)
        CALL delete_Vector (pwrden%qsawe)
        CALL delete_Vector (pwrden%qsawi)
        CALL delete_Vector (pwrden%qrad)
        CALL delete_Vector (pwrden%qohm)

        CALL delete_Vector (dischg%psivolpnpsi)
        CALL delete_Vector (dischg%elongxnpsi)
        CALL delete_Vector (dischg%triangnpsi_u)
        CALL delete_Vector (dischg%triangnpsi_l)
        CALL delete_Vector (dischg%pindentnpsi)
        CALL delete_Vector (dischg%sfareanpsi)
        CALL delete_Vector (dischg%cxareanpsi)
        CALL delete_Vector (dischg%grho1npsi)
        CALL delete_Vector (dischg%grho2npsi)
        CALL delete_Vector (dischg%rplasbdry)
        CALL delete_Vector (dischg%zplasbdry)
        CALL delete_Vector (rho_mhd_gridnpsi)
        CALL delete_Vector (pwrden%qexch)
        CALL delete_vector (mhd_dat%btotrmaj)
        CALL delete_vector (mhd_dat%bprmaj)
        CALL delete_vector (mhd_dat%bp)
        CALL delete_vector (mhd_dat%ffprim)
        CALL delete_vector (mhd_dat%pprim)
        CALL delete_vector (mhd_dat%ravgnpsi)
        CALL delete_vector (mhd_dat%ravginpsi)
        CALL delete_vector (mhd_dat%rcap)
        CALL delete_vector (mhd_dat%rcapi)
        CALL delete_vector (profile%press)
        CALL delete_vector (profile%pressb)
        CALL delete_vector (dischg%rlimiter)
        CALL delete_vector (dischg%zlimiter)
        CALL delete_vector (dischg%rmhdgrid)
        CALL delete_vector (dischg%zmhdgrid)





!        DEALLOCATE(scx_ncd,srecom_ncd,sion_ncd,sbcx_ncd ,      &
!                   stsource,dudtsv_ncd,enbeam_ncd,enn_ncd)    
        DEALLOCATE(scx_ncd,sion_ncd,sbcx_ncd ,                  &
                   stsource,dudtsv_ncd,enbeam_ncd,enn_ncd)    
        DEALLOCATE (ennv_ncd,ennw_ncd,volsn_ncd,sbeame,        &
                   sbeam_ncd,currf_ncd,zeff_ncd,storqueb_ncd,  &
                   z_ncd,zsq_ncd,wbeam_ncd,walp_ncd)          
        DEALLOCATE (eps_ncd,xhm2_ncd,                          &
                   xi33_ncd,xi11_ncd,work,mhd_dat%psi)
        DO jj = 1,nion
           CALL delete_Vector (profile%en(jj))
        ENDDO
        DO jj =1 ,ntot
            CALL delete_vector(profile%flux(jj))
            CALL delete_vector(profile%flux_conv(jj))
        ENDDO
        
        DEALLOCATE(profile%flux)
        DEALLOCATE(profile%flux_conv)
        DEALLOCATE(profile%en)


     RETURN
     END SUBROUTINE set_onetwo_vars

     




     SUBROUTINE set_gcnmp_vars
!-----------------------------------------------------------------------------
!--- we have to map between gcnmp and onetwo variables due to
!    1) DIFFERENT UNITS,GRIDS
!    2) DIFFERENT ASSUMPTIONS MADE BY NETCDF AND TEXT FILE ITERDB OUTPUT
!       ROUTINES. THE TEXTFILE READ AND WRITE ROUTINE, iter_dbase_txt,
!       USES ONETWO VARIABLES FOR BOTH INPUT AND OUTPUT.
!       THE NETCDF READ/WRITE ROUTINE ,iter_dbase_nc, USES GCNMP VARIABLES
!       WHICH ARE TRANSLATED TO ONETWO VARIABLES BY SET_ONETWO_VARS
!       AND FOR WRITING ONETWO VARS WHICH ARE TRANSLATED TO GCNMP VARS BY THIS ROUTINE
!------------------------------------------------------------------------------
 

     USE PARAM,                                  ONLY : kj,kion,ke

     USE mhdcom,                                 ONLY : mhdmethd
     USE gpsi,                                   ONLY : p
     USE yoka,                                   ONLY : pbe,pbi
     USE flags,                                  ONLY : no_te_convection, &
                                                        no_ti_convection
     USE transp,                                 ONLY : nubeam_version


     IMPLICIT NONE 
     INTEGER(I4B) izenstop,lden_bc,isz,nb
     LOGICAL set_u
     REAL(DP) sumpow,five_halfs_ti,five_halfs_te

        shot_id%shot_nmbr         = ishot ! from yoka
        time_ncd                  = time
        tGCNMF                    = timmax
        nj_ncd                    = nj
        nion_ncd                  = nion
        nprim_ncd                 = nprim
        nimp_ncd                  = nimp
        nneu_ncd                  = nneu
        he3_frac_ncd              = he3_frac
        he3_thermal_spin_pol_ncd  = he3_thermal_spin_pol
        nk                        = nion + dp4  
        neut_beam%nj_beam         = nj
        neut_beam%d_beam_spin_pol = d_beam_spin_pol
      IF(use_nubeam)THEN
         isz = MIN(SIZE(nameb_nubeam),SIZE(beam_data%label_beamions))
         nbion_ncd              = beam_data%nbeam_species
         fd_beam_ncd            = 1.0_DP    ! no 'dt' beams for nubeam
         nameb_nubeam(1:isz)    = beam_data%label_beamions(1:isz)
         neut_beam%nbeams       = beam_data%nbeam ! for statefile write

         CALL  neut_beam_allocate

     
         sumpow = zeroc
         DO nb = 1,neut_beam%nbeams

            neut_beam%fbcur(1,nb)   = beam_data%ffulla(nb)
            neut_beam%fbcur(2,nb)   = beam_data%fhalfa(nb)
            neut_beam%fbcur(3,nb)   = MAX(0.0_DP,1._DP - beam_data%ffulla(nb) - beam_data%fhalfa(nb))

          neut_beam%prompt_pwr_in_plasma(1,nb) = beam_data%pinja(nb)* neut_beam%fbcur(1,nb) 
          neut_beam%prompt_pwr_in_plasma(2,nb) = beam_data%pinja(nb)* neut_beam%fbcur(2,nb) 
          neut_beam%prompt_pwr_in_plasma(3,nb) = beam_data%pinja(nb)* neut_beam%fbcur(3,nb) 
          neut_beam%ebeam(1,nb) =  beam_data%einja(nb)/1000._DP
          neut_beam%ebeam(2,nb) =  beam_data%einja(nb)/2000._DP
          neut_beam%ebeam(3,nb) =  beam_data%einja(nb)/3000._DP
          sumpow = sumpow + beam_data%pinja(nb)
         ENDDO
         !Note pbeam is supposed to be power to apertures
         ! use this until we get info from nubeam
         neut_beam%pbeam(:,:)   =  neut_beam%prompt_pwr_in_plasma(:,:)


          ! approximate beam power going into torus by assuming
          ! that it is the same as thermalized power:
          ! user injected power fractions to assign to each beamlet
          ! bptor is supposed to be power past apertures  
          ! note that if thebeam was just turned on then]
          ! the thermalized power is small and hence a veruy poor approximation.
          DO nb =1,neut_beam%nbeams
             bptor(nb)             = 1.e6_DP*(pbe+pbi)*beam_data%pinja(nb)/sumpow 
          ENDDO                                          



      ELSE
         nbion_ncd              = 1 
         fd_beam_ncd            = fdbeam
         nameb_nubeam(1)        = nameb
      ENDIF


 
        dischg%nr_mhd           =  nw               !required below, dont move
        dischg%nz_mhd           =  nh
        dischg%nplasbdry        =  nplasbdry
        fd_thermal_ncd          =  fd_thermal
        ntot                    =  nion + dp4   


! ----ALLOCATE based on above settings of **_ncd -------------------
           IF(.NOT. ALLOCATED(work))ALLOCATE (work(nj)) 
           IF( .NOT. ASSOCIATED(namep_ncd))   &
              ALLOCATE ( namep_ncd(nprim),namei_ncd(nimp), &
                         namen_ncd(nneu),nameb_ncd(nbion_ncd))
           IF( .NOT. ASSOCIATED(profile%en))THEN
              ALLOCATE (profile%en(nion))                                               !2d array 
              DO j=1,nion
                 profile%en(j) = zero_Vector(nj)
              ENDDO
           ENDIF
           IF( .NOT. ASSOCIATED(profile%flux))THEN
              ALLOCATE (profile%flux(ntot))                                            !2d array
              DO j=1,ntot
                 profile%flux(j) = zero_Vector(nj)
              ENDDO
           ENDIF



           IF( .NOT. ASSOCIATED(profile%flux_conv))THEN
              ALLOCATE (profile%flux_conv(ntot))                                       
              DO j=1,ntot
                 profile%flux_conv(j) = zero_Vector(nj)
              ENDDO
           ENDIF




!--- tglf fluxes are not know when we are going to write a statefile
!--- because Onetwo does not calculate them
!--- (tglf fluxes are  calcualted in gcnmp)
!--- it is assumed that the second dimension is 3 always
           IF(.NOT. ALLOCATED(tglf_p_output))THEN
              ALLOCATE(tglf_p_output(nj,3))
              tglf_p_output(:,:) = zeroc
          ELSEIF(SIZE(tglf_p_output,1) .NE. nj)THEN
              DEALLOCATE(tglf_p_output)
              ALLOCATE(tglf_p_output(nj,3))
              tglf_P_output(:,:) = zeroc   !! should interpolate old to new 
           ENDIF
           IF(.NOT. ALLOCATED(tglf_e_output))THEN
              ALLOCATE(tglf_e_output(nj,3))
              tglf_e_output(:,:) = zeroc
           ELSEIF(SIZE(tglf_e_output,1) .NE. nj)THEN
              DEALLOCATE(tglf_e_output)
              ALLOCATE(tglf_e_output(nj,3))
              tglf_e_output(:,:) = zeroc    !! should interpolate old to new
           ENDIF
           IF(.NOT. ALLOCATED(tglf_m_output))THEN
              ALLOCATE(tglf_m_output(nj,3))
              tglf_m_output(:,:) = zeroc
           ELSEIF(SIZE(tglf_m_output,1) .NE. nj)THEN
              DEALLOCATE(tglf_m_output)
              ALLOCATE(tglf_m_output(nj,3))
              tglf_m_output(:,:) = zeroc    !! should interpolate old to new
           ENDIF


!---- glf fluxes assume they are calcualted in gcnmp so just zero them here:
          IF(.NOT. ALLOCATED(glf_p_output))THEN
                   ALLOCATE(glf_p_output(nj,3))
                   glf_p_output(:,:) = zeroc
          ELSEIF(SIZE(glf_p_output,1) .NE. nj)THEN
              DEALLOCATE(glf_p_output)
              ALLOCATE(glf_p_output(nj,3))
           ENDIF

          IF(.NOT. ALLOCATED(glf_e_output))THEN
                   ALLOCATE(glf_e_output(nj,3))
          ELSEIF(SIZE(glf_e_output,1) .NE. nj)THEN
              DEALLOCATE(glf_e_output)
              ALLOCATE(glf_e_output(nj,3))
           ENDIF

          IF(.NOT. ALLOCATED(glf_m_output))THEN
                   ALLOCATE(glf_m_output(nj,3))
          ELSEIF(SIZE(glf_m_output,1) .NE. nj)THEN
              DEALLOCATE(glf_m_output)
              ALLOCATE(glf_m_output(nj,3))
           ENDIF

           glf_p_output(:,:) = zeroc
           glf_e_output(:,:) = zeroc 
           glf_m_output(:,:) = zeroc
 
          IF(.NOT. ALLOCATED(glf_etg_output))THEN
                   ALLOCATE(glf_etg_output(nj))
          ELSEIF(SIZE(glf_etg_output,1) .NE. nj)THEN
              DEALLOCATE(glf_etg_output)
              ALLOCATE(glf_etg_output(nj))
           ENDIF

          IF(.NOT. ALLOCATED(glf_anfreq_output))THEN
                   ALLOCATE(glf_anfreq_output(nj))
          ELSEIF(SIZE(glf_anfreq_output,1) .NE. nj)THEN
              DEALLOCATE(glf_anfreq_output)
              ALLOCATE(glf_anfreq_output(nj))
           ENDIF

          IF(.NOT. ALLOCATED(glf_anfreq2_output))THEN
                   ALLOCATE(glf_anfreq2_output(nj))
          ELSEIF(SIZE(glf_anfreq2_output,1) .NE. nj)THEN
              DEALLOCATE(glf_anfreq2_output)
              ALLOCATE(glf_anfreq2_output(nj))
           ENDIF


          IF(.NOT. ALLOCATED(glf_anrate_output))THEN
                   ALLOCATE(glf_anrate_output(nj))
          ELSEIF(SIZE(glf_anrate_output,1) .NE. nj)THEN
              DEALLOCATE(glf_anrate_output)
              ALLOCATE(glf_anrate_output(nj))
           ENDIF

          IF(.NOT. ALLOCATED(glf_anrate2_output))THEN
                   ALLOCATE(glf_anrate2_output(nj))
          ELSEIF(SIZE(glf_anrate2_output,1) .NE. nj)THEN
              DEALLOCATE(glf_anrate2_output)
              ALLOCATE(glf_anrate2_output(nj))
           ENDIF
           glf_etg_output(:)     = zeroc
           glf_anrate_output(:)  = zeroc
           glf_anrate2_output(:) = zeroc
           glf_anfreq_output(:)  = zeroc
           glf_anfreq2_output(:) = zeroc

          IF(.NOT. ALLOCATED(glf_gamma_net_i_output))THEN
                   ALLOCATE(glf_gamma_net_i_output(nj))
          ELSEIF(SIZE(glf_gamma_net_i_output,1) .NE. nj)THEN
              DEALLOCATE(glf_gamma_net_i_output)
              ALLOCATE(glf_gamma_net_i_output(nj))
          ENDIF

          IF(.NOT. ALLOCATED(glf_gamma_net_e_output))THEN
                   ALLOCATE(glf_gamma_net_e_output(nj))
          ELSEIF(SIZE(glf_gamma_net_e_output,1) .NE. nj)THEN
              DEALLOCATE(glf_gamma_net_e_output)
              ALLOCATE(glf_gamma_net_e_output(nj))
          ENDIF
          glf_gamma_net_i_output(:) = zeroc
          glf_gamma_net_e_output(:) = zeroc





           IF(.NOT. ASSOCIATED(prtcl_src%srecom))ALLOCATE(prtcl_src%srecom(nion))
           IF(.NOT. ASSOCIATED(wpdot%dpidt))ALLOCATE(wpdot%dpidt(nion))
           IF(.NOT. ALLOCATED(te_bc))ALLOCATE(te_bc(nj))
           IF(.NOT. ALLOCATED(ti_bc))ALLOCATE(ti_bc(nj))
           IF(.NOT. ALLOCATED(ene_bc))ALLOCATE(ene_bc(nj))
           IF(.NOT. ALLOCATED(zeff_bc))ALLOCATE(zeff_bc(nj))
           IF(.NOT. ALLOCATED(angrot_bc))ALLOCATE(angrot_bc(nj))
           IF(.NOT. ALLOCATED(eps_ncd))ALLOCATE(eps_ncd(nj))
           IF(.NOT. ALLOCATED(storqueb_ncd))ALLOCATE(storqueb_ncd(nj))
           IF(.NOT. ALLOCATED(wbeam_ncd))ALLOCATE(wbeam_ncd(nj))
           IF(.NOT. ALLOCATED(walp_ncd))ALLOCATE(walp_ncd(nj))
           IF(.NOT. ALLOCATED(enalp_ncd))ALLOCATE(enalp_ncd(nj))
           IF(.NOT. ALLOCATED(dnedt_ncd))ALLOCATE(dnedt_ncd(nj))
           IF(.NOT. ALLOCATED(dangrotdt_ncd))ALLOCATE(dangrotdt_ncd(nj))
           IF(.NOT. ASSOCIATED(zeff_ncd))ALLOCATE(zeff_ncd(nj))
           IF(.NOT. ASSOCIATED(z_ncd))ALLOCATE(z_ncd(nj,nion))
           IF(.NOT. ALLOCATED(en_bc))ALLOCATE(en_bc(nj,nion))
           IF(.NOT. ALLOCATED(flux_bc))ALLOCATE(flux_bc(nj,ntot))
           IF(.NOT. ALLOCATED(dnidt_ncd))ALLOCATE(dnidt_ncd(nj,nion))
           IF(.NOT. ASSOCIATED(zsq_ncd))ALLOCATE(zsq_ncd(nj,nion))
           IF(ALLOCATED(rev_psival))DEALLOCATE(rev_psival)
           ALLOCATE(rev_psival(npsi))
           IF(ALLOCATED(rev_generic))DEALLOCATE(rev_generic)
           ALLOCATE(rev_generic(npsi))
           IF(.NOT. ALLOCATED(xi11_ncd))ALLOCATE(xi11_ncd(nj_ncd))
           IF(.NOT. ALLOCATED(xi33_ncd))ALLOCATE(xi33_ncd(nj_ncd))
           IF(.NOT. ALLOCATED(xhm2_ncd))ALLOCATE(xhm2_ncd(nj_ncd))
           IF(.NOT. ALLOCATED(xips_ncd))ALLOCATE(xips_ncd(nj_ncd))
           IF(.NOT. ALLOCATED(currf_ncd))ALLOCATE(currf_ncd(nj_ncd))
           IF(.NOT. ALLOCATED(curbeam))ALLOCATE(curbeam(nj_ncd))
           IF(.NOT. ALLOCATED(enn_ncd))ALLOCATE(enn_ncd(nj,nneu))
           IF(.NOT. ALLOCATED(ennv_ncd))ALLOCATE(ennv_ncd(nj,nneu))
           IF(.NOT. ALLOCATED(ennw_ncd))ALLOCATE(ennw_ncd(nj,nneu))
           IF(.NOT. ALLOCATED(volsn_ncd))ALLOCATE(volsn_ncd(nj,nneu))
           IF(.NOT. ALLOCATED(sbeame))ALLOCATE(sbeame(nj))
           IF(.NOT. ALLOCATED(enbeam_ncd))ALLOCATE(enbeam_ncd(nj,nbion_ncd))
           IF(.NOT. ALLOCATED(sbeam_ncd))ALLOCATE(sbeam_ncd(nj,nbion_ncd))
           IF(.NOT. ALLOCATED(stsource))ALLOCATE(stsource(nion,nj))
           IF(.NOT. ALLOCATED(sion_ncd))ALLOCATE(sion_ncd(nj,nion))
!           IF(.NOT. ALLOCATED(srecom_ncd))ALLOCATE(srecom_ncd(nj,nion))
           IF(.NOT. ALLOCATED(sbcx_ncd))ALLOCATE(sbcx_ncd(nj,nion))
           IF(.NOT. ALLOCATED(scx_ncd))ALLOCATE(scx_ncd(nj,nion))
           IF(.NOT. ALLOCATED(dudtsv_ncd))ALLOCATE(dudtsv_ncd(ntot,nj))
           IF(.NOT. ALLOCATED(fix_edge_ni_bc))ALLOCATE(fix_edge_ni_bc(nion_ncd))
           IF(.NOT. ALLOCATED(rcap_ncd))ALLOCATE(rcap_ncd(nj_ncd))
           IF(.NOT. ALLOCATED(r2cap_ncd))ALLOCATE(r2cap_ncd(nj_ncd))
           IF(.NOT. ALLOCATED(r2capi_ncd))ALLOCATE(r2capi_ncd(nj_ncd))
           IF(.NOT. ALLOCATED(dfdt_ncd))ALLOCATE(dfdt_ncd(nj_ncd))
           IF(.NOT. ALLOCATED(dgdt_ncd))ALLOCATE(dgdt_ncd(nj_ncd))
           IF(.NOT. ALLOCATED(dhdt_ncd))ALLOCATE(dhdt_ncd(nj_ncd))
           IF(.NOT. ASSOCIATED( mhd_dat%psi)) ALLOCATE(mhd_dat%psi(dischg%nr_mhd,dischg%nz_mhd))
           IF(.NOT. ASSOCIATED(diffuse%dcoef))THEN 
              ALLOCATE(diffuse%dcoef(ntot,ntot,nj_ncd)) ! only nj-1 elements Are used
                                  ! set to  size nj for output purposes
              diffuse%dcoef(:,:,:) = 0.0_DP
           ENDIF
           IF(.NOT. ASSOCIATED(diffuse%xnus))THEN 
              ALLOCATE(diffuse%xnus(nion_ncd))
              diffuse%xnus(:) = zero_vector(nj_ncd)
           ENDIF
          DO jj =1,nion_ncd
             diffuse%xnus(jj)%data(1:nj_ncd-1) = xnus(jj,:)
          ENDDO
              
          diffuse%ftrap               =  new_Vector(nj,ftrap)
          
          !diffuse%chie_paleo         =  assignrt_Vector(chie_paleo,im22icm2) ! diffuse%chie_paleo  in m^2/sec
          !diffuse%eta                =  assignrt_Vector(eta,8.98755e9) ! ohm m
          diffuse%eta                 =  zero_Vector(nj_ncd)
          diffuse%chie_paleo          =  zero_Vector(nj_ncd)
          !eta,chie_paleo are  dimmed to kj (not nj):
          DO  j=1,nj_ncd
             diffuse%eta%data(j)           =  eta(j)*8.98755e9 ! ohm m, 
             diffuse%chie_paleo%data(j)    =  chie_paleo(j)*im22icm2 ! diffuse%chie_paleo  in m^2/sec
          ENDDO

          diffuse%xnuse                  =  zero_Vector(nj_ncd)
          diffuse%xnuse%data(1:nj_ncd-1) =  xnuse(1:nj_ncd-1)

         DO jj =1,nbion_ncd                       ! nbion_ncd =1 if nubeam_version .ge. 201107
            IF(use_nubeam) THEN                   ! because nubeam only puts out single density
               IF (nubeam_version .lt. 201107) then 
                  enbeam_ncd(1:nj,jj) = enbeam_species_c(1:nj,jj)/im32icm3 
               ELSE
                   enbeam_ncd(1:nj,jj) = enbeam(1:nj)/im32icm3 
               ENDIF
            ELSE
               enbeam_ncd(1:nj,jj)    = enbeam(1:nj)/im32icm3
            ENDIF
            sbeam_ncd(1:nj,jj)     = sbeam(1:nj)/im32icm3 ! no multiple beam species
         ENDDO

         sbeame(:)                 = sbion(:)/im32icm3
         IF(use_nubeam .and. nubeam_version .ge. 201107)THEN
             sbeame(:)             = sbeam(1:nj)/im32icm3
         ENDIF


         storqueb_ncd(1:nj)        = storqueb(1:nj)*pconvert
         !storqueb_ncd in mks (nt m/m^3) , storqueb in dyne cm /cm**3
         enalp_ncd(1:nj_ncd)       = enalp(1:nj)/im32icm3
         wbeam_ncd(1:nj_ncd)       = wbeam(1:nj)/im32icm3   ! kev/m**3
         walp_ncd(1:nj_ncd)        = walp(1:nj)/im32icm3    ! kev/m**3
         profile%walp              = new_Vector(nj,walp_ncd)

       DO jj=1,nneu
          enn_ncd(1:nj,jj)   = enn(1:nj,jj)/im32icm3
          ennv_ncd(1:nj,jj)  = ennv(1:nj,jj)/im32icm3
          ennw_ncd(1:,jj)    = ennw(1:,jj)/im32icm3
          volsn_ncd(1:nj,jj) = volsn(1:nj,jj)/im32icm3
       ENDDO

       DO jj=1,nion
          zsq_ncd(1:nj,jj)= zsq(1:nj,jj)
          z_ncd(1:nj,jj)   = z(1:nj,jj)
       ENDDO

           !the species are ordered as 1..nprim,1..nimp
           !where nion = nprim+nimp
           DO jj=1,nion_ncd
              profile%en(jj)    = new_Vector(nj,en(1,jj))
              profile%en(jj)    = Vector_mult_real(profile%en(jj),1./im32icm3)
              en_bc(1:nj,jj)       = en(1:nj,jj)/im32icm3 ! there currently is no time
                                                    ! dependent input in inone for 
                                                    ! densities so just pass  
           ENDDO                                    ! on the time independent profiles

      
           ! thermal-thermal ddn rate:
           ! fus_prod%neutr_ddn_th = assignrt_Vector(ddfusn,1.e6_DP) #/(m^3 sec) wont work 
           ! if kj ne nj  so do the following instead:
           fus_prod%neutr_ddn_th                = zero_Vector(nj_ncd)
          DO  j=1,nj_ncd
             fus_prod%neutr_ddn_th%data(j)      = ddfusn(j)*1.e6_DP ! #/(m^3 sec) dimmed to kj (not nj)
          ENDDO

 
           ! beam-beam ddn rate summed over all beamlines:
           work(1:nj)                          = beam_beamddn(1:nj,((3*nbeams+1)*3*nbeams)/2+1)
           fus_prod%neutr_ddn_beam_beam        = assignrt_Vector(work,1.e6_DP)   ! #/(m^3 sec)

           ! beam thermal ddn  rate  summed over all beam lines:
           work(1:nj)                          = beam_thermalddn(1:nj,3*nbeams+1 )
           fus_prod%neutr_ddn_beam_thermal     = assignrt_Vector(work,1.e6_DP)   ! #/(m^3 sec)
 
           ! knock on ddn reactions:


           !fus_prod%neutr_ddn_knock            = assignrt_Vector(ddknck,1.e6_DP) ! #/(m^3 sec)
          fus_prod%neutr_ddn_knock              =  zero_Vector(nj_ncd)
          DO  j=1,nj_ncd
             fus_prod%neutr_ddn_knock%data(j)   = ddknck(j)*1.e6_DP ! #/(m^3 sec) dimmed to kj (not nj)
          ENDDO

           ! total neutron rate:
           fus_prod%neutr_ddn_tot              = copy_Vector(fus_prod%neutr_ddn_beam_beam)
           fus_prod%neutr_ddn_tot%data(:)      = fus_prod%neutr_ddn_th%data(:)           +    &
                                                 fus_prod%neutr_ddn_beam_thermal%data(:) +    &
                                                 fus_prod%neutr_ddn_tot%data(:)         ! #/(m^3 sec)



      !     fus_prod%alpha_dthe4_th               = new_Vector(nj,ddfusn?)
      !     fus_prod%alpha_dthe4_beam             = new_Vector(nj,ddfusn?)
      !     fus_prod%alpha_dthe4_beam_thermal     = new_Vector(nj,ddfusn?)
      !     fus_prod%alpha_dthe4_knock            = new_Vector(nj,ddfusn?)
      !     fus_prod%alpha_dthe4_tot              = new_Vector(nj,ddfusn?)




           fus_prod%p_he3d_th                = new_Vector(nj,he3dpfus)
           CALL Vector_mult(fus_prod%p_he3d_th,M32cm3) !1/M^3 sec
           fus_prod%p_ddt                 = new_Vector(nj,ddpfus)
           CALL Vector_mult(fus_prod%p_ddt,M32cm3) !1/M^3 sec
           fus_prod%total_p_ddt           = thermal_thermal_ddptot ! note only thermal at preent
           fus_prod%total_p_he3d              = thermal_thermal_he3dptot + beam_thermalhe3th_dftot       
           fus_prod%total_p_he3d_th           = thermal_thermal_he3dptot
           fus_prod%total_p_he3d_beam_thermal = beam_thermalhe3th_dftot

           fus_prod%p_he3d_beam_thermal   = new_Vector(nj,beam_thermalhe3th_df(1,nbeams+1))
           CALL Vector_mult(fus_prod%p_he3d_beam_thermal,M32cm3) !1/M^3 sec

           fus_prod%total_neutr_ddn_th            = thermal_thermal_ddntot 
           fus_prod%total_neutr_ddn_beam_beam     = beam_beam_ddntot 
           fus_prod%total_neutr_ddn_beam_thermal  = beam_thermal_ddntot 
           fus_prod%total_neutr_ddn_knock         = ddknct  
           fus_prod%total_neutr_ddn               = ddntot ! see sub fiziks 307.F 




           flux_bc(1:nj,1:ntot) = zeroc            


           !the species are ordered as 1..nprim,1..nimp
           !where nion = nprim+nimp , ntot = nion + 4
           work(nj) = zeroc
           DO jj=1,ntot
              work(1:nj-1) = flux(jj,1:nj-1)
              profile%flux(jj)  = new_Vector(nj,work)
              IF(jj == nion_ncd +1 .AND. itran(jj) == 0)THEN
                   profile%flux(jj)       = new_Vector(nj,anal_eng_flux_e)  
              ENDIF
              IF(jj == nion_ncd + 2 .AND. itran(jj) == 0)THEN
                   profile%flux(jj)       = new_Vector(nj,anal_eng_flux_i)
              ENDIF
           ENDDO
           DO jj =1 ,nion_ncd
              profile%flux(jj)   = Vector_mult_real(profile%flux(jj),M22cm2)
           ENDDO
           profile%flux(nion_ncd+1)  = Vector_mult_real(profile%flux(nion_ncd+1),joupkev*M22cm2)
           profile%flux(nion_ncd+2)  = Vector_mult_real(profile%flux(nion_ncd+2),joupkev*M22cm2)
           profile%flux(nion_ncd+3)  = Vector_mult_real(profile%flux(nion_ncd+3),T2gauss)
           profile%flux(nion_ncd+4)  = Vector_mult_real(profile%flux(nion_ncd+4),g2kg)



           DO jj=1,nion
              profile%flux_conv(jj)%data(1) = zeroc
              DO j = 2,nj-1
                 profile%flux_conv(jj)%data(j) = 0.5_DP*(vpinch(j)+ vpinch(j-1))*en(j,jj)*M22cm2
              ENDDO
              profile%flux_conv(jj)%data(nj) = profile%flux_conv(jj)%data(nj-1) 
            ENDDO
            five_halfs_te = 2.5_DP
            five_halfs_ti = 2.5_DP
            IF(no_te_convection .EQ. -1)five_halfs_te = 1.5_DP
            IF(no_te_convection .EQ. 1) five_halfs_te = zeroc
            IF(no_ti_convection .EQ. -1)five_halfs_ti = 1.5_DP
            IF(no_ti_convection .EQ. 1) five_halfs_ti  = zeroc

           profile%fluxe                         = zero_Vector(nj)
           profile%flux_conv(nion_ncd+1)         = zero_Vector(nj)
           profile%fluxe%data(1)                 = zeroc
           profile%flux_conv(nion_ncd+1)%data(1) = zeroc
           DO j = 2,nj-1
              profile%fluxe%data(j) = (fluxe(j-1)+fluxe(j))*0.5_DP
              profile%flux_conv(nion_ncd+1)%data(j)= five_halfs_te * profile%fluxe%data(j)*te(j)
           ENDDO
           profile%fluxe%data(nj)                = fluxe(nj-1) ! approximate
           profile%flux_conv(nion_ncd+1)%data(nj)= five_halfs_te * profile%fluxe%data(nj)
           profile%fluxe                         = Vector_mult_real(profile%fluxe,M22cm2)
           profile%flux_conv(nion_ncd+1)         = Vector_mult_real(profile%flux_conv(nion_ncd+1),M22cm2*joupkev)



           profile%fluxi                         = zero_Vector(nj)
           profile%flux_conv(nion_ncd+2)         = zero_Vector(nj)
           profile%fluxi%data(1)                 = zeroc
           profile%flux_conv(nion_ncd+2)%data(1) = zeroc
           DO j=2,nj-1
              profile%fluxi%data(j) = (fluxi(j-1)+fluxi(j))*0.5_DP
              profile%flux_conv(nion+2)%data(j)= five_halfs_ti * profile%fluxi%data(j)*ti(j)
           ENDDO
           profile%fluxi%data(nj)   = fluxi(nj-1) ! approximate
           profile%flux_conv(nion+2)%data(nj)= five_halfs_ti * profile%fluxi%data(nj)
           profile%fluxi            = Vector_mult_real(profile%fluxi,M22cm2)
           profile%flux_conv(nion+2)= Vector_mult_real(profile%flux_conv(nion_ncd+2),M22cm2*joupkev)
 

           profile%flux_conv(nion_ncd+3)         = zero_Vector(nj) ! no converctiver term for Faraday's law

          
           ! convective term in toroidal rotation:
           profile%flux_conv(nion_ncd+4)         = zero_Vector(nj)
           profile%flux_conv(nion_ncd+4)%data(1) = zeroc
           DO jj=2,nj-1
              profile%flux_conv(nion_ncd+4)%data(j)= (fluxangc(j-1)+fluxangc(j))*0.5_DP
           ENDDO
           profile%flux_conv(nion_ncd+4)%data(nj)= profile%flux_conv(nion_ncd+4)%data(nj-1)
           profile%flux_conv(nion_ncd+4)= Vector_mult_real(profile%flux_conv(nion_ncd+4),0.001_DP) ! kg/sec^2
 



          sion_ncd(:,:)   = zeroc
!          srecom_ncd(:,:) = zeroc
          sbcx_ncd(:,:)   = zeroc
          scx_ncd(:,:)    = zeroc
          stsource(:,:)   = zeroc
          DO jj=1,MIN(2,nion)     ! onetwo only has 2 species here
             sion_ncd(1:nj,jj)    = sion(1:nj,jj)/im32icm3
!             srecom_ncd(1:nj,jj) = srecom(1:nj,jj)/im32icm3
             sbcx_ncd(1:nj,jj)    = sbcx(1:nj,jj)/im32icm3
             scx_ncd(1:nj,jj)     = scx(1:nj,jj)/im32icm3
             stsource(jj,1:nj)    = s(jj,1:nj)/im32icm3
          ENDDO

          DO jj=1,ntot
             IF(jj .LE. nion)cfa    = 1._DP/im32icm3
             IF((jj == nion+1)  .OR.  ( jj == nion+2))cfa = joupkev
             IF( jj == nion+3)  cfa = Tm2gcm
             IF( jj == nion+4)  cfa = 1._DP
             dudtsv_ncd(jj,1:nj)    = dudtsv(jj,1:nj)*cfa
          ENDDO

        namep_ncd(1:nprim_ncd) =  namep(1:nprim)
        namei_ncd(1:nimp_ncd)  =  namei(1:nimp) 
        namen_ncd(1:nneu_ncd)  =  namen(1:nneu) 
        nameb_ncd(1:nbion_ncd) =  nameb_nubeam(1:nbion)

        profile%te             =  new_vector(nj,te)
        profile%ti             =  new_vector(nj,ti)
        profile%ene            =  new_vector(nj,ene)
        profile%ene            =  Vector_mult_real (profile%ene,1._DP/im32icm3)
        profile%etor           =  new_vector(nj,etor)
        profile%etor           =  Vector_mult_real(profile%etor,M2cm)
        profile%angrot         =  new_vector(nj,angrot)
        profile%press          =  new_Vector(nj,press)
        profile%press          =  Vector_mult_real (profile%press, pconvert)

        profile%pressb         = new_Vector(nj,pressb)
        profile%pressb         = Vector_mult_real (profile%pressb, pconvert)



        
        zeff_ncd(1:nj)         = zeff(1:nj)
!        dnedt_ncd(1:nj)       = dnedt(1:nj)/im32icm3
        wpdot%dnedt            = new_Vector(nj,dnedt)
        wpdot%dnedt            = Vector_mult_real (wpdot%dnedt, 1._DP/im32icm3)
        DO jj =1,nion
           dnidt_ncd(1:nj,jj)  = dnidt(1:nj,jj)/ im32icm3
        ENDDO
        xi11_ncd(1:nj) = xi11(1:nj)
        xi33_ncd(1:nj) = xi33(1:nj) 
        xhm2_ncd(1:nj) = xhm2(1:nj)
        xips_ncd(1:nj) = xips(1:nj)
        

          DO j=1,ntot
           DO jj=1,ntot
              IF(j .LE. nion)THEN     !density equations
                 IF( jj .LE. nion)THEN
                    diffuse%dcoef(j,jj,1:nj-1) = d(j,jj,1:nj-1)/M22CM2    ! dcoef is in M^2/sec
                 ELSE IF( jj == nion+1 .OR. jj == nion+2)THEN
                    diffuse%dcoef(j,jj,1:nj-1) = d(j,jj,1:nj-1)*M2cm         ! work is in 1/(kev m sec)
                 ELSE IF( jj == nion+3)THEN
                    diffuse%dcoef(j,jj,1:nj-1) = d(j,jj,1:nj-1)*T2gauss*M22cm2  ! work is in 1/(tesla  m^2 sec)
                 ELSE IF( jj == nion+4)THEN
                    diffuse%dcoef(j,jj,1:nj-1) = d(j,jj,1:nj-1)*M2cm ! work is in 1/m
                 ENDIF
              ELSEIF( j == nion+1 .OR. j == nion+2)THEN
                     IF( jj .LE. nion)THEN
                        diffuse%dcoef(j,jj,1:nj-1) = d(j,jj,1:nj-1) *im22icm2 ! work is in (M^2 kev/sec)
                     ELSE IF(  jj == nion+1 .OR. jj == nion+2)THEN
                        diffuse%dcoef(j,jj,1:nj-1) = d(j,jj,1:nj-1)*M2cm ! work is in  1/(m sec)
                     ELSE IF( jj == nion+3)THEN
                        diffuse%dcoef(j,jj,1:nj-1) = d(j,jj,1:nj-1)* T2gauss*M22cm2 ! work is in (kev tesla)/(m^2 sec)
                     ELSE IF(  jj == nion+4)THEN
                        diffuse%dcoef(j,jj,1:nj-1) = d(j,jj,1:nj-1)* M2cm ! work is in (kev/m)
                     ENDIF
              ELSEIF(j == nion+3)THEN
                     IF( jj .LE. nion)THEN
                        diffuse%dcoef(j,jj,1:nj-1) = d(j,jj,1:nj-1)* M22cm2*M22cm2/T2gauss ! work is in (tesla M^4/sec)
                     ELSE IF(  jj == nion+1 .OR. jj == nion+2)THEN
                        diffuse%dcoef(j,jj,1:nj-1) = d(j,jj,1:nj-1)*Tm2gcm ! work is in (Tesla m)/(kev sec)
                     ELSE IF( jj == nion+3)THEN
                        diffuse%dcoef(j,jj,1:nj-1) = d(j,jj,1:nj-1)  ! work is in 1/sec
                     ELSE IF(  jj == nion+4)THEN
                        diffuse%dcoef(j,jj,1:nj-1) = d(j,jj,1:nj-1) * Tm2gcm ! work is in (Tesla m)
                     ENDIF
              ELSEIF(j == nion+4)THEN
                     IF( jj .LE. nion)THEN
                        diffuse%dcoef(j,jj,1:nj-1) = d(j,jj,1:nj-1) * im22icm2*im22icm2 ! work is in (kg m^4/sec^2)
                     ELSE IF(  jj == nion+1 .OR. jj == nion+2)THEN
                        diffuse%dcoef(j,jj,1:nj-1) = d(j,jj,1:nj-1)* g2kg/M2cm         ! work is in (kg m)/(kev sec^2)
                     ELSE IF( jj == nion+3)THEN
                        diffuse%dcoef(j,jj,1:nj-1) = d(j,jj,1:nj-1)* g2kg*T2gauss      ! work is in (kg/(tesla sec^2))
                     ELSE IF(  jj == nion+4)THEN
                        diffuse%dcoef(j,jj,1:nj-1) = d(j,jj,1:nj-1)* g2kg *M2cm        ! work is in (kg m/sec)
                     ENDIF
              ENDIF
           ENDDO
         ENDDO



        diffuse%chieinv  = new_vector(nj,chieinv)
        diffuse%chieinv  = Vector_mult_real(diffuse%chieinv,1./m22cm2)
        diffuse%chiinv   = new_vector(nj,chiinv)
        diffuse%chiinv   = Vector_mult_real(diffuse%chiinv,1./m22cm2)
        diffuse%xkineo   = new_vector(nj,xkineo)
        diffuse%xkineo   = Vector_mult_real(diffuse%xkineo,M2cm)
        diffuse%xkeneo   = new_vector(nj,xkeneo)
 
        diffuse%xkeneo   = Vector_mult_real(diffuse%xkeneo,M2cm)
 

        IF(.NOT. ALLOCATED(glf_eff_chi_output))THEN
           ALLOCATE(glf_eff_chi_output(nj,ntot))
           glf_eff_chi_output(:,:) = zeroc
        ENDIF


        ! total ion pressure, load only first comonent:
        DO j=1,nion
           work(:)         = dnidt(:,j)*ti(:)
           wpdot%dpidt(j)  = new_Vector(nj,work)
           wpdot%dpidt(j)  = Vector_mult_real(wpdot%dpidt(j),1.5*convert)
        ENDDO



        wpdot%dpedt     = new_Vector(nj,dpedtc)
        wpdot%dpedt     = Vector_mult_real(wpdot%dpedt,1.5*convert)
        work(:)         = dpedtc(:)/(te(:)*im32icm3)  ! #/m^3
        wpdot%dnedt     = new_Vector(nj,work)


        
        dischg%rmajor          =  rmajor / m2cm 
        dischg%kappa           =  kappa 
        dischg%pindento        =  pindentnpsi(1)
        dischg%volo            =  psivolp(1) * im32icm3  
        dischg%areao           =  cxareanpsi(1)* im22icm2
        dischg%deltao          =  triangnpsi(1)
        dischg%circum          =  circum/ m2cm

!        dischg%rgeom          =  rmajavnpsi(1)/m2cm 
!        dischg%rmag           =  rmajavnpsi(npsi)/m2cm  
        dischg%rmag            =  rma/m2cm  
        dischg%rplasbdry       =  new_Vector(dischg%nplasbdry,rplasbdry)
        dischg%rplasbdry       =  Vector_mult_real (dischg%rplasbdry,1./M2cm)
        dischg%zplasbdry       =  new_Vector(dischg%nplasbdry,zplasbdry)
        dischg%zplasbdry       =  Vector_mult_real (dischg%zplasbdry,1./M2cm)
        !find boundary condition time closest to time up to
        !which gcnmp is to be run
        !gcnmp will be run from time0 to time_max (as set in 
        !namelist write routine write_gcnmp_namelist)

        IF(time_max == 0.0_DP)time_max = tGCNMF 
        ibctim =0
        IF(nbctim == 1)THEN
           time_bc             =  bctime(nbctim)
           ibctim              =  1
        ELSE
           DO jj = 1,nbctim-1
              IF(bctime(jj) .LE. time_max .AND.  time_max   &
                   .LE. bctime(jj+1))ibctim = jj+1
           ENDDO
           
           time_bc             =  bctime(ibctim) ! pass onetwo bc vars at this time
        ENDIF
        nbctim_ncd             = nbctim
        IF(ibctim == 0)ibctim  = nbctim ! HSJ 2/15/2012
!
! --- take care of a roundoff problem
! --- changed to relative measure HSJ 8/25/98
!
      IF  (psir(nj) .NE. 0.0) THEN
        IF (ABS((psival(1)-psir(nj))/(psir(nj))) .LT. 1.0e-6) &
                                           psival(1) = psir(nj)
      ELSE
        IF (ABS(psival(1)-psir(nj)) .LT. 1.0e-10) &
                                           psival(1) = psir(nj)
      END IF
      psir_grid              =  new_Vector(nj,psir)
      psir_grid              =  Vector_mult_real (psir_grid,psimks)


      rho_grid               =  new_Vector(nj,r)
      rho_grid               =  Vector_mult_real (rho_grid,1./M2cm)


!
! --- convert vectors defined on the npsi (i.e., MHD) grid to corresponding
! --- quantities on the rho (i.e., transport) grid
!

      bpar(1) = 0.0               ! set boundary conditions on spline
      bpar(2) = 0.0
      bpar(3) = 0.0
      bpar(4) = 0.0



!------------------------------------------------------------------------
       go to 10000
      IF(ALLOCATED(work))DEALLOCATE(work)
      ALLOCATE(work(nj_ncd))

             periodic       = .FALSE. !interpolant is not periodic
             uniform        = .TRUE.  ! tension is uniform
             ncd            = 2       ! # continuous derivatives
             iendc          = 0       ! natural spline bc
            !iendc          = 1       ! give deriv bc at end points
             bpar(1)        = 0.0     ! used only if iendc =1,2 
             bpar(2)        = 0.0     ! used only if iendc =1,2 
            ! NOTE -- TENSION PARAMETER IS OPTIONAL ARGUMENT IN 
            ! t716_TSVAL1 :
             CALL reverse(npsi,psival,rev_psival)
             CALL reverse(npsi,rmajavnpsi,rev_generic)
             CALL  t716_TSPSI (npsi,rev_psival,rev_generic,ncd,iendc, &
                               periodic,uniform,bpar)
            ! IFLAG = 0 if values of function are to be computed.
            ! IFLAG = 1 if first derivative values are to be computed.
            ! IFLAG = 2 if second derivative values are to be computed.

             iflag = 0
             CALL  t716_TSVAL1 (npsi,rev_psival,rev_generic,IFLAG,nj_ncd,    &
                                psir,work)
             dischg%rmajavnpsi        =  new_Vector(nj_ncd,work)
             dischg%rmajavnpsi        =  Vector_mult_real (dischg%rmajavnpsi,1._dp/m2cm)

             CALL reverse(npsi,triangnpsi,rev_generic)
             CALL  t716_TSPSI (npsi,rev_psival,rev_generic,ncd,iendc, &
                               periodic,uniform,bpar)
             CALL  t716_TSVAL1 (npsi,rev_psival,rev_generic,IFLAG,nj_ncd,    &
                                psir,work)
             dischg%triangnpsi_u   = new_Vector(nj_ncd,work)

             CALL reverse(npsi,triangnpsi_l,rev_generic)
             CALL  t716_TSPSI (npsi,rev_psival,rev_generic,ncd,iendc,       &
                               periodic,uniform,bpar)
             CALL  t716_TSVAL1 (npsi,rev_psival,rev_generic,IFLAG,nj_ncd,    &
                                psir,work)
             dischg%triangnpsi_l  = new_Vector(nj_ncd,work)


 

             CALL reverse(npsi,rminavnpsi,rev_generic)
             CALL  t716_TSPSI (npsi,rev_psival,rev_generic,ncd,iendc, &
                               periodic,uniform,bpar)
             CALL  t716_TSVAL1 (npsi,rev_psival,rev_generic,IFLAG,nj_ncd,    &
                                psir,work)
             dischg%rminavnpsi =   new_Vector(nj_ncd,work)
             dischg%rminavnpsi = Vector_mult_real (dischg%rminavnpsi,1._dp/m2cm)   
             CALL reverse(npsi,psivolp,rev_generic)
             CALL  t716_TSPSI (npsi,rev_psival,rev_generic,ncd,iendc, &
                               periodic,uniform,bpar)
             CALL  t716_TSVAL1 (npsi,rev_psival,rev_generic,IFLAG,nj_ncd,    &
                                psir,work)             
             dischg%psivolpnpsi        =   new_Vector(nj_ncd,work)
             dischg%psivolpnpsi        =   Vector_mult_real (dischg%psivolpnpsi, im32icm3)


             CALL reverse(npsi,elongx,rev_generic)
             CALL  t716_TSPSI (npsi,rev_psival,rev_generic,ncd,iendc, &
                               periodic,uniform,bpar)
             CALL  t716_TSVAL1 (npsi,rev_psival,rev_generic,IFLAG,nj_ncd,    &
                                psir,work)             
             dischg%elongxnpsi        =   new_Vector(nj_ncd,work)

                                      

             CALL reverse(npsi,pindentnpsi,rev_generic)
             CALL  t716_TSPSI (npsi,rev_psival,rev_generic,ncd,iendc, &
                               periodic,uniform,bpar)
             CALL  t716_TSVAL1 (npsi,rev_psival,rev_generic,IFLAG,nj_ncd,    &
                                psir,work)             
             dischg%pindentnpsi       =   new_Vector(nj_ncd,work)

                                                  
             CALL reverse(npsi,sfareanpsi,rev_generic) 
             CALL  t716_TSPSI (npsi,rev_psival,rev_generic,ncd,iendc, &
                               periodic,uniform,bpar)
             CALL  t716_TSVAL1 (npsi,rev_psival,rev_generic,IFLAG,nj_ncd,    &
                                psir,work)             
             dischg%sfareanpsi       =   new_Vector(nj_ncd,work)
             dischg%sfareanpsi       = Vector_mult_real (dischg%sfareanpsi, im22icm2)


             CALL reverse(npsi,cxareanpsi,rev_generic)
             CALL  t716_TSPSI (npsi,rev_psival,rev_generic,ncd,iendc, &
                               periodic,uniform,bpar)
             CALL  t716_TSVAL1 (npsi,rev_psival,rev_generic,IFLAG,nj_ncd,    &
                                psir,work)             
             dischg%cxareanpsi       =   new_Vector(nj_ncd,work)
             dischg%cxareanpsi       =   Vector_mult_real (dischg%cxareanpsi, im22icm2)
           

             CALL reverse(npsi,grho1npsi,rev_generic)
             CALL  t716_TSPSI (npsi,rev_psival,rev_generic ,ncd,iendc, &
                               periodic,uniform,bpar)
             CALL  t716_TSVAL1 (npsi,rev_psival,rev_generic,IFLAG,nj_ncd,    &
                                psir,work)             
             dischg%grho1npsi      =   new_Vector(nj_ncd,work)

                      
             CALL reverse(npsi,grho2npsi,rev_generic)
             CALL  t716_TSPSI (npsi,rev_psival,rev_generic,ncd,iendc, &
                               periodic,uniform,bpar)
             CALL  t716_TSVAL1 (npsi,rev_psival,rev_generic,IFLAG,nj_ncd,    &
                                psir,work)             
             dischg%grho2npsi       =   new_Vector(nj_ncd,work)

!-----------------------------------------------------------------------------
! --- output full npsi grid profiles in same revere order as they exists in Onetwo:
!-----------------------------------------------------------------------------
10000        CONTINUE

             dischg%rmajavnpsi     = new_Vector(npsi,rmajavnpsi)
             dischg%rmajavnpsi     = Vector_mult_real (dischg%rmajavnpsi,1._dp/m2cm)
             dischg%triangnpsi_u   = new_Vector(npsi,triangnpsi)
             dischg%triangnpsi_l   = new_Vector(npsi,triangnpsi_l)
             dischg%rminavnpsi     = new_Vector(npsi,rminavnpsi)
             dischg%rminavnpsi     = Vector_mult_real (dischg%rminavnpsi,1._dp/m2cm)
             dischg%psivolpnpsi    = new_Vector(npsi,psivolp)
             dischg%psivolpnpsi    = Vector_mult_real (dischg%psivolpnpsi, im32icm3)
             dischg%elongxnpsi     = new_Vector(npsi,elongx)
             dischg%pindentnpsi    = new_Vector(npsi,pindentnpsi)
             dischg%sfareanpsi     = new_Vector(npsi,sfareanpsi)
             dischg%sfareanpsi     = Vector_mult_real (dischg%sfareanpsi, im22icm2)
             dischg%cxareanpsi     = new_Vector(npsi,cxareanpsi)
             dischg%cxareanpsi     = Vector_mult_real (dischg%cxareanpsi, im22icm2)
             dischg%grho1npsi      = new_Vector(npsi,grho1npsi)
             dischg%grho2npsi      = new_Vector(npsi,grho2npsi)

             dischg%rma            = rma *cm2M
             dischg%zma            = zma *cm2M
             dischg%rsep           = rsep *cm2M
             dischg%zsep           = zsep *cm2M
             dischg%rplasmax       = rplasmax *cm2M
             dischg%zplasmax       = zplasmax *cm2M
             dischg%rplasmin       = rplasmin *cm2M
             dischg%zplasmin       = zplasmin *cm2M
             dischg%rgeom          = rgeom*cm2m
             dischg%btgeom         = btgeom/T2gauss
             dischg%rmhdgrid       = new_Vector(dischg%nr_mhd,rmhdgrid)
             dischg%rmhdgrid       = Vector_mult_real(dischg%rmhdgrid,cm2m)
             dischg%zmhdgrid       = new_Vector(dischg%nz_mhd,zmhdgrid)
             dischg%zmhdgrid       = Vector_mult_real(dischg%zmhdgrid,cm2m)
             
             dischg%nlimiter       = nlimiter
             dischg%rlimiter       = new_Vector(nlimiter,xlimiter)
!             dischg%rlimiter       = Vector_mult_real(dischg%rlimiter,cm2m)
             dischg%zlimiter       = new_Vector(nlimiter,ylimiter)
!             dischg%zlimiter       = Vector_mult_real(dischg%zlimiter,cm2m)
!           xlimiter ,ylimiter are in meters here  HSJ , 1/20/2011


!------------------------------------------------------------------------------
! --- load mhd_dat
!------------------------------------------------------------------------------

             mhd_dat%npsi           =  npsi
             mhd_dat%btor           =  btor/T2GAUSS
             mhd_dat%rbp            =  new_Vector(nj,rbp) 
             mhd_dat%rbp            =  Vector_mult_real (mhd_dat%rbp,Tm2gcm)
             mhd_dat%betan          =  zero_vector(nj) ! betan calcs are all in gcnmp

             mhd_dat%tot_cur        =  tocur
             mhd_dat%totohm_cur     =  totohm
             mhd_dat%totboot_cur    =  totboot
             mhd_dat%totbeam_cur    =  totbeam
             mhd_dat%totrf_cur      =  totrf
             mhd_dat%betap          =  betap
             mhd_dat%beta           =  beta
             mhd_dat%ali            =  ali

             mhd_dat%rcap           =  new_Vector(nj,rcap)
             mhd_dat%rcap           =  Vector_mult_real (mhd_dat%rcap, cm2m)
             mhd_dat%rcapi          =  new_Vector(nj,rcapi)
             mhd_dat%rcapi          =  Vector_mult_real (mhd_dat%rcapi, M2cm)
             rcap_ncd               =  get_values(mhd_dat%rcap)
             mhd_dat%r2cap          =  new_Vector(nj,r2cap)
             r2cap_ncd              =  get_values(mhd_dat%r2cap)
             mhd_dat%r2capi         =  new_Vector(nj,r2capi)
             mhd_dat%r2capi         =  Vector_mult_real(mhd_dat%r2capi,im22icm2)
             r2capi_ncd             =  get_values(mhd_dat%r2capi)

             eps_ncd(1:nj_ncd)      =  eps(1:nj_ncd)
             dfdt_ncd(1:nj_ncd)     =  dfdt(1:nj_ncd)
             dgdt_ncd(1:nj_ncd)     =  dgdt(1:nj_ncd)
             dhdt_ncd(1:nj_ncd)     =  dhdt(1:nj_ncd)
             mhd_dat%fcap           =  new_Vector(nj,fcap)
             mhd_dat%gcap           =  new_Vector(nj,gcap)
             mhd_dat%hcap           =  new_Vector(nj,hcap)
             mhd_dat%q_value        =  new_Vector(nj,q)



             mhd_dat%curden         =  new_Vector(nj,curden)  
             mhd_dat%curden         =  Vector_mult_real (mhd_dat%curden, 1./im22icm2)

             mhd_dat%curpar         = new_Vector(nj,curpar_soln)
             mhd_dat%curpar         = Vector_mult_real (mhd_dat%curpar, 1./im22icm2)

             mhd_dat%curohm         =  new_Vector(nj,curohm)
             mhd_dat%curohm         =  Vector_mult_real (mhd_dat%curohm, 1./im22icm2)

             mhd_dat%curboot        =  new_Vector(nj,curboot)
             mhd_dat%curboot        =  Vector_mult_real (mhd_dat%curboot, 1./im22icm2)

             mhd_dat%pprim     = new_Vector(nj,pprim)
             mhd_dat%pprim     = Vector_mult_real(mhd_dat%pprim,ppmks)

             mhd_dat%ffprim    = new_Vector(nj,ffprim)
             mhd_dat%ffprim    = Vector_mult_real(mhd_dat%ffprim,ffpmks)




             mhd_dat%bp        = new_Vector(nj,bp)
             mhd_dat%bp        = Vector_mult_real(mhd_dat%bp,kgauss2t)


             mhd_dat%bprmaj    = new_Vector(nj,bprmaj)
             mhd_dat%bprmaj    = Vector_mult_real(mhd_dat%bprmaj, kgauss2t)
            

             mhd_dat%btotrmaj  = new_Vector(nj,btotrmaj)
             mhd_dat%btotrmaj  = Vector_mult_real(mhd_dat%btotrmaj,kgauss2t)

             mhd_dat%ravgnpsi      = new_Vector(npsi,ravg)
             mhd_dat%ravgnpsi      = Vector_mult_real (mhd_dat%ravgnpsi,cm2M)

             mhd_dat%ravginpsi     = new_Vector(npsi,ravgi)
             mhd_dat%ravginpsi     = Vector_mult_real (mhd_dat%ravginpsi,M2cm)

 
             IF(mhdmethd .NE. 'tdem')THEN
                  mhd_dat%psi(:,:)  = psi(:,:) * psimks
             ELSE
                mhd_dat%psi(:,:)  = p(:,:)* psimks
             ENDIF

             mhd_dat%psiaxis   = psiaxis  * psimks
             mhd_dat%psibdry   = psibdry  * psimks

             ok = delete_vector_nf(mhd_dat%psivalnpsi)
             ok = delete_vector_nf(mhd_dat%qpsinpsi)
             ok = delete_vector_nf(mhd_dat%pressnpsi)
             ok = delete_vector_nf(mhd_dat%ffprimnpsi)
             ok = delete_vector_nf(mhd_dat%pprimnpsi)
             ok = delete_Vector_nf (mhd_dat%fpsinpsi)
             mhd_dat%qpsinpsi    = new_Vector(npsi,qpsi)
             mhd_dat%pressnpsi   = new_Vector(npsi,presspsi)    ! ergs/cm^3 ==> nt/m^2
             mhd_dat%pressnpsi   = Vector_mult_real(mhd_dat%pressnpsi,0.1_DP)
             mhd_dat%ffprimnpsi  = new_Vector(npsi,ffppsi)      ! gauss ==> kg/(A sec^2)
             mhd_dat%ffprimnpsi  = Vector_mult_real ( mhd_dat%ffprimnpsi,1.e-04_DP)
             mhd_dat%pprimnpsi   = new_Vector(npsi,pppsi)       ! gram/gauss*cm**3*sec**2 => A/m^3
             mhd_dat%pprimnpsi   = Vector_mult_real ( mhd_dat%pprimnpsi,1.e7_DP) 
             mhd_dat%psivalnpsi  = new_Vector(npsi,psival)
             mhd_dat%psivalnpsi  = Vector_mult_real (mhd_dat%psivalnpsi, psimks)
             mhd_dat%fpsinpsi    = new_Vector(npsi,fpsi)
             mhd_dat%fpsinpsi    = Vector_mult_real(mhd_dat%fpsinpsi, Tm2gcm) !tesal m
!              call set_eqdsk_statefile_values
!             mhd_dat%fpsinpsi      = new_Vector(npsi,fpsi)
!             mhd_dat%fpsinpsi      = Vector_mult_real(mhd_dat%fpsinpsi, Tm2gcm) !tesal m


             currf_ncd(1:nj_ncd)    =  currf(1:nj)/im22icm2
             curbeam(1:nj_ncd)      =  curdbeam(1:nj)/im22icm2
             rho_mhd_gridnpsi       =  new_Vector(npsi,rho)
             rho_mhd_gridnpsi       =  Vector_mult_real (rho_mhd_gridnpsi,cm2M)

             rho_gridn              =  new_Vector(nj,roa)


             IF(ALLOCATED(work))DEALLOCATE(work)
             ALLOCATE(work(nj))
!             DO jj=1,nion        ! HSJ 7/13/11 arrays like srecom(kj,2) ,see sourc.f90
             DO jj=1,nprim        ! hence cant use nion here
                   work(1:nj)            =  srecom(1:nj,jj)/im32icm3
                   prtcl_src%srecom(jj)  =  assign_Vector(work)
             ENDDO
             work(1:nj)            =  zeroc
             DO jj=nprim+1,nion
                   prtcl_src%srecom(jj)  =  assign_Vector(work) ! no sources for impurities
             ENDDO






             prtcl_src%stfuse       = new_Vector(nj,stfus)
             prtcl_src%sbfuse       = new_Vector(nj,sbfus)
             prtcl_src%spellet      = new_Vector(nj,spellet)
             prtcl_src%stfuse       = Vector_mult_real(prtcl_src%stfuse,1._DP/im32icm3 )
             prtcl_src%sbfuse       = Vector_mult_real(prtcl_src%sbfuse,1._DP/im32icm3 )
             prtcl_src%spellet      = Vector_mult_real(prtcl_src%spellet,1._DP/im32icm3 )

             pellet%name            = nampel


!------------------------------------------------------------------------------
! --- load pwrden and related
!------------------------------------------------------------------------------

             IF(.NOT. ASSOCIATED(pwrden%brems_nions))THEN
                ALLOCATE(pwrden%brems_nions(nion))

                DO j=1,nion
                   pwrden%brems_nions(j) = assignrt_Vector (brems_nions(:,j),joupkev)   ! watts/m^3
                ENDDO
             ENDIF

             pwrden%qconde          = new_Vector(nj,qconde)
             pwrden%qconde          = Vector_mult_real(pwrden%qconde,convert)

             pwrden%qcondi          = new_Vector(nj,qcondi)
 
             pwrden%qcondi          = Vector_mult_real(pwrden%qcondi,convert)

             pwrden%qconve          = new_Vector(nj,qconve)
             pwrden%qconve          = Vector_mult_real(pwrden%qconve,convert)
             pwrden%qconvi          = new_Vector(nj,qconvi)
             pwrden%qconvi          = Vector_mult_real(pwrden%qconvi,convert)

             pwrden%qbeame          = new_Vector(nj,qbeame)
             pwrden%qbeame          = Vector_mult_real(pwrden%qbeame,convert)
             pwrden%qbeami          = new_Vector(nj,qbeami)
             pwrden%qbeami          = Vector_mult_real(pwrden%qbeami,convert)

             pwrden%qrfe          = new_Vector(nj,qrfe)
             pwrden%qrfe          = Vector_mult_real(pwrden%qrfe,convert)
             pwrden%qrfi          = new_Vector(nj,qrfi)
             pwrden%qrfi          = Vector_mult_real(pwrden%qrfi,convert)

             pwrden%qione         = new_Vector(nj,qione)
             pwrden%qione         = Vector_mult_real(pwrden%qione,convert)
             pwrden%qioni         = new_Vector(nj,qioni)
             pwrden%qioni         = Vector_mult_real(pwrden%qioni,convert)

             pwrden%qe2d         = new_Vector(nj,qe2d)
             pwrden%qe2d         = Vector_mult_real(pwrden%qe2d,convert)
             pwrden%qi2d         = new_Vector(nj,qi2d)
             pwrden%qi2d         = Vector_mult_real(pwrden%qi2d,convert)

             pwrden%qcx         = new_Vector(nj,qcx)
             pwrden%qcx         = Vector_mult_real(pwrden%qcx,convert) 

             pwrden%qfuse        = new_Vector(nj,qfuse)
             pwrden%qfuse        = Vector_mult_real(pwrden%qfuse,convert)
             pwrden%qfusi        = new_Vector(nj,qfusi)
             pwrden%qfusi        = Vector_mult_real(pwrden%qfusi,convert)

             pwrden%qbfuse       = new_Vector(nj,qbfuse)
             pwrden%qbfuse       = Vector_mult_real(pwrden%qbfuse,convert)
             pwrden%qbfusi       = new_Vector(nj,qbfusi)
             pwrden%qbfusi       = Vector_mult_real(pwrden%qbfusi,convert)

             pwrden%qmag         = new_Vector(nj,qmag)
             pwrden%qmag         = Vector_mult_real(pwrden%qmag,convert)

             pwrden%qsawe        = new_Vector(nj,qsawe)
             pwrden%qsawe        = Vector_mult_real(pwrden%qsawe,convert)
             pwrden%qsawi        = new_Vector(nj,qsawi)
             pwrden%qsawi        = Vector_mult_real(pwrden%qsawi,convert)


             pwrden%qrad         = new_Vector(nj,qrad)
             pwrden%qrad         = Vector_mult_real(pwrden%qrad,convert)
             pwrden%qohm         = new_Vector(nj,qohm)
             pwrden%qohm         = Vector_mult_real(pwrden%qohm,convert)
             pwrden%qdelt        = new_Vector(nj,qdelt)
             pwrden%qdelt        = Vector_mult_real(pwrden%qdelt,convert)
             pwrden%qexch        = new_Vector(nj,qexch)
             pwrden%qexch        = Vector_mult_real(pwrden%qexch,convert)



             pwrden%omegale      = new_Vector(nj,omegale)
             pwrden%omegale      = Vector_mult_real(pwrden%omegale,convert)
             pwrden%qangce       = new_Vector(nj,qangce)
             pwrden%qangce       = Vector_mult_real(pwrden%qangce,convert)      
             pwrden%qomegapi     = new_Vector(nj,qomegapi)
             pwrden%qomegapi     = Vector_mult_real(pwrden%qomegapi,convert)            
             pwrden%sprcxre      = new_Vector(nj,sprcxre)
             pwrden%sprcxre      = Vector_mult_real(pwrden%sprcxre,convert)
             pwrden%spreimpe     = new_Vector(nj,spreimpe)
             pwrden%spreimpe     = Vector_mult_real(pwrden%spreimpe,convert)
             pwrden%sprcxree     = new_Vector(nj,sprcxree)
             pwrden%sprcxree     = Vector_mult_real(pwrden%sprcxree,convert)


        !----allocate and load plasma frequencies
         nprimp1 = nprim+1 ! ==>ions_gcnmp
         IF(.NOT. ASSOCIATED(atw_ncd))ALLOCATE(atw_ncd(nion))
         atw_ncd(:) = atw(1:nion)
         CALL common_frequencies



             IF(ABS(totcur(ibctim)) .GT. 1._DP)THEN
                   totcur_bc = totcur(ibctim)
             ELSE
                   totcur_bc = totcur(1)
             ENDIF

             vloop_bc_ncd        = vloop_bc(ibctim)     
             fix_edge_te_bc      = fix_edge_te(ibctim)
             fix_edge_ti_bc      = fix_edge_ti(ibctim)
             fix_edge_rot_bc     = fix_edge_rot(ibctim)

             !assumes all nion_ncd values are equal to single 
             !input value in onetwo:
             fix_edge_ni_bc(1:nion_ncd) = fix_edge_ni_bc_inpt(1:nion_ncd,ibctim) 

            ! use convention that any values > 1.0 are outp


            ! onetwo input vars (tein,tiin,etc.)
            ! are defined on spline knot grid. we do not pass splines to
            ! iterdb files. Hence need to account for this here:
            ! tspline settings:
              periodic       = .FALSE. !interpolant is not periodic
              uniform        = .TRUE.  ! tension is uniform
              ncd            = 2       ! # continuous derivatives
              iendc          = 0       ! natural spline bc
            ! iendc         = 1        ! give deriv bc at end points
              bpar(1)        = 0.0     ! used only if iendc =1,2 
              bpar(2)        = 0.0     ! used only if iendc =1,2 
            ! NOTE -- TENSION PARAMETER IS OPTIONAL ARGUMENT IN 
            ! t716_TSVAL1 :
             iflag = 0                 !means get function values (instead of derivs)

         ! LOAD te_bc:
              CALL  get_knot_num(rtein(:,ibctim),knots,ksplin)

              CALL  t716_TSPSI (knots,rtein(:,ibctim),tein(:,ibctim),&
                               ncd,iendc, periodic,uniform,bpar)

              CALL  t716_TSVAL1 (knots,rtein(:,ibctim),tein(:,ibctim),            &
                                IFLAG,nj_ncd,roa,te_bc)

         ! LOAD ti_bc:
              CALL  get_knot_num(rtiin(:,ibctim),knots,ksplin)
              CALL  t716_TSPSI (knots,rtiin(:,ibctim),tiin(:,ibctim),             &
                               ncd,iendc, periodic,uniform,bpar)
              CALL  t716_TSVAL1 (knots,rtiin(:,ibctim),tiin(:,ibctim),            &
                                 IFLAG,nj_ncd,roa,ti_bc)

         ! LOAD angrot_bc:
               CALL  get_knot_num(rangrot(:,ibctim),knots,ksplin)
               CALL  t716_TSPSI (knots,rangrot(:,ibctim),angrotin(:,ibctim),      &
                               ncd,iendc, periodic,uniform,bpar)
               CALL  t716_TSVAL1 (knots,rangrot(:,ibctim),angrotin(:,ibctim),     &
                                IFLAG,nj_ncd,roa,angrot_bc)


      IF(inenez == 1 .OR. inenez == -99)THEN  !enein and zeff were input:
         ! LOAD ene_bc:
              CALL  get_knot_num(renein(:,ibctim),knots,ksplin)
              CALL  t716_TSPSI (knots,renein(:,ibctim),enein(:,ibctim),           &
                               ncd,iendc, periodic,uniform,bpar)
              CALL  t716_TSVAL1 (knots,renein(:,ibctim),enein(:,ibctim),          &
                                IFLAG,nj_ncd,roa,ene_bc)
               ene_bc(:) = ene_bc(:)/im32icm3   ! ene_bc in 1/m**3
         ! LOAD zeff_bc from input values:
              IF(inenez == 1)THEN
                 CALL  get_knot_num(rzeffin(:,ibctim),knots,ksplin)
                 CALL  t716_TSPSI (knots,rzeffin(:,ibctim),zeffin(:,ibctim),      &
                      ncd,iendc, periodic,uniform,bpar)
                 CALL  t716_TSVAL1 (knots,rzeffin(:,ibctim),zeffin(:,ibctim),     &
                      IFLAG,nj_ncd,roa,zeff_bc)
         ! LOAD ion densities,en_bc, at bctime(ibctim) from ene,zeff input:
         ! given ene_bc,zeff_bc,te_bc  from above get z,zsq and then en_bc:
         ! We do not know enbeam and enalpha at time bctime(ibctim) SO
         ! the best estimate is their current value.

 
                 set_dzdt =.FALSE.
                 CALL get_charge_state(te_bc)
                 lden_bc = SIZE(en_bc,1)

                 set_u = .FALSE. ! Do not reset u on this call
                 CALL en_calc(ene_bc,zeff_bc,en_bc,lden_bc,izenstop,set_u) ! en_bc in #/m^3 
 
                 set_dzdt = .TRUE.

              ENDIF


      ENDIF

      IF(inenez == 0)THEN    ! primary and impurity densities were input
         ! LOAD ion densities:
              DO j=1,nion
                 CALL  t716_TSPSI (ken_bc_inpt(j,ibctim),ren_bc_inpt(:,j,ibctim), &
                                  en_bc_inpt(:,j,ibctim),ncd,iendc, periodic,     &
                                  uniform,bpar)
                 CALL  t716_TSVAL1 (ken_bc_inpt(j,ibctim),ren_bc_inpt(:,j,ibctim),&
                                    en_bc_inpt(:,j,ibctim),IFLAG,nj_ncd,roa,      &
                                    en_bc(:,j))
                 en_bc(:,j) = en_bc(:,j)/im32icm3
              ENDDO

 
          ! LOAD ene_bc,zeff_bc at time bctime(ibctim) from 
          ! en_bc,enalp_ncd,enbeam_ncd (all are in 1/m^3) here
              CALL get_ene_bc
      
       ENDIF

       IF (inenez == -99)THEN        !ene and impurities were input,
                                     !get primary and zeff. Picked up ene
                                     !under inenez = 1. Now get impurities
                                     ! from input:
          DO j=nprim+1,nion
                 CALL  t716_TSPSI (ken_bc_inpt(j,ibctim),ren_bc_inpt(:,j,ibctim), &
                                  en_bc_inpt(:,j,ibctim),ncd,iendc, periodic,     &
                                  uniform,bpar)
                 CALL  t716_TSVAL1 (ken_bc_inpt(j,ibctim),ren_bc_inpt(:,j,ibctim),&
                                    en_bc_inpt(:,j,ibctim),IFLAG,nj_ncd,roa,      &
                                    en_bc(:,j))
                 en_bc(:,j) = en_bc(:,j)/im32icm3
          ENDDO

         
       ENDIF
      IF(.NOT.(inenez .EQ. 0 .OR. inenez .EQ. 1 .OR. inenez .EQ. -99))THEN
         WRITE(nout,1000)inenez
         WRITE(ncrt,1000)inenez
1000     FORMAT(2x,'ERROR : inenez =',i5,' not curently implemented')
         CALL STOP('set_gcnmp_vars',1)
         
      ENDIF


      IF(ALLOCATED(rev_psival))DEALLOCATE(rev_psival)
      IF(ALLOCATED(rev_generic))DEALLOCATE(rev_generic)

   RETURN

   END SUBROUTINE set_gcnmp_vars


   SUBROUTINE get_ene_bc
! -----------------------------------------------------------------------
! --- given thermal and fast ion densities get electron density:
! --- Note that this routine relies on modules included in common 
! --- section  at beginning of this module.
! --- beam chg state == 1 !!
! --- units for en_bc,enbeam_bc,enalp_bc will determine units of ene_bc
! -----------------------------------------------------------------------
   IMPLICIT NONE
   INTEGER(I4B) j, k
      ene_bc(:)  = zeroc
      zeff_bc(:) = zeroc
      DO j=1,nj
        DO k=1,nion
          ene_bc(j) = ene_bc(j)+z(j,k)*en_bc(j,k)
          zeff_bc(j) = zeff_bc(j)+zsq(j,k)*en_bc(j,k)
        END DO
        DO k = 1,nbion_ncd
           ene_bc(j) = ene_bc(j) + enbeam_ncd(j,k)
           zeff_bc(j) = zeff_bc(j) + enbeam_ncd(j,k)
        ENDDO
        ene_bc(j)  = ene_bc(j)  + 2.0 * enalp_ncd(j)
        zeff_bc(j) = zeff_bc(j) + 4.0 * enalp_ncd(j)
        zeff_bc(j) = zeff_bc(j)/ene_bc(j)
      ENDDO

   RETURN
   END SUBROUTINE get_ene_bc


   SUBROUTINE  check_statefile_equilibrium(rma_in,zma_in,psiax_in)
!-----------------------------------------------------------------------------
! -- consitency check of equilibirum in  statfile
!-----------------------------------------------------------------------------
      USE bicube,                        ONLY : cspln,wnoperm,pds,n2cspln,nh2
      USE replace_imsl,                  ONLY : my_ibcccu,my_dbcevl1
      USE mhdpar,                        ONLY : nh,nwh,nw
      use mhdcom,                        ONLY : zero
      USE limiter,                       ONLY : nlimiter,xlimiter,ylimiter

      IMPLICIT NONE
      INTEGER(i4B) ier ,icalc,iknowax,ispln,isignpsi, &
                   ioout,nxs,nys
      REAL(DP) psiax_in,psiaxisl,elongaxl,rma_in,zma_in,rmal,zmal,psimal,psi_in, &
           rdum,zdum
      REAL(DP),DIMENSION(:),ALLOCATABLE ::  xliml,yliml


      ALLOCATE(xliml(nlimiter),yliml(nlimiter))
      xliml = xlimiter*100._DP ; yliml = ylimiter*100._DP


      nxs = size(rmhdgrid)   ; nys =size(zmhdgrid)
      ioout    = 777
!      WRITE(ioout,FMT='("nw,nh,nxs,nys =",4(x,i3))')nw,nh,nxs,nys

      iknowax  = 0  ! mag axis not known
      ispln    = 0  ! cspln not known
      isignpsi = -1  ! search for min in psi
      rmal     = 0.0_DP  ; zmal =0.0_DP

      ! --- set up array zero ( = 1 inside limiter, =0 outside limiter)
      iflag = 1
      call zlim (zero,nw,nh,nlimiter,xliml,yliml,rmhdgrid, &
                      zmhdgrid,iflag)
!      write(ioout,FMT='("nlimiter =",i5)')nlimiter
!      DO j=1,nlimiter
!         WRITE(ioout,FMT='("j,xlim,ylim =",i5,x,2(x,1pe12.4))')j,xliml(j),yliml(j)
!      ENDDO
 
      call magax (psi,rmhdgrid,zmhdgrid,nxs,nys,isignpsi,iknowax,  &
                 ioout,ispln,zero,rmal,zmal,cspln,n2cspln,nh2,     &
                 psiaxisl,elongaxl)
!      WRITE(ioout,FMT='("magax reports xaxis,yaxis,psiaxis,elong =",4(x,1pe16.8))')rmal,zmal,psiaxisl,elongaxl
      if (iknowax .eq. 0) &
        call STOP ('sub check_statefile_equilibrium: problem 1', 1)

      call my_ibcccu (psi,rmhdgrid,nxs,zmhdgrid,nys,cspln,nw,wnoperm,ier)
      icalc = 1
      call my_dbcevl1 (rmhdgrid,nxs,zmhdgrid,nys,cspln,nw,rmal,zmal,pds, &
                       ier,icalc)
      psimal = pds(1)
      call my_dbcevl1 (rmhdgrid,nxs,zmhdgrid,nys,cspln,nw,rma_in,zma_in,pds, &
                       ier,icalc)
      psi_in = pds(1)
!     WRITE(ioout,FMT='("rma_in,zma_in,psiax calc =",3(x,1pe16.8))')rma_in,zma_in,pds(1)
     DEALLOCATE(xliml,yliml)
!!     WRITE(ioout,FMT='("rma_in,rmal =",2(x,1pe16.8))')rma_in,rmal    
!     WRITE(ioout,FMT='("zma_in,zmal =",2(x,1pe16.8))')zma_in,zmal
!     WRITE(ioout,FMT='("psiax_in, psiax(rmal,zmal) =",2(x,1pe16.8))')psiax_in,psimal
!     WRITE(ioout,FMT='("psiax_in, psiax(rma_in,zma_in) =",2(x,1pe16.8))')psiax_in,psi_in
!     WRITE(ioout,FMT='("psiax(rma_in,zma_in)  - psi(rmal,zmal) =",1pe16.8)')psi_in - psimal

     rma_in = rmal
     zma_in = zmal

   RETURN
   END SUBROUTINE  check_statefile_equilibrium


   END MODULE set_12_gcnmp_vars
