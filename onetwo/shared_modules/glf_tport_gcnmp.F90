 MODULE GLF_tport_gcnmp
     USE nrtype,                 ONLY : Dp,I4B

 
     USE dep_var,                ONLY : te,ti,angrot,ene,en,dp4

     USE error_handler,          ONLY : iomaxerr, lerrno,terminate

     USE glf23_gcnmp,            ONLY : te_glf,ti_glf,shat_exp,alpha_exp,        &
                                        rho_glf,drhogf,ni_glf,ne_glf,rmin_glf,   &
                                        rmaj_glf,drhodr,chie_m,q_glf,vp_eff_m,   &
                                        ve_eff_m,vi_eff_m,vrot_eff_m,chii_m,     &
                                        zpti_glf,zpte_glf,zpne_glf,zpni_glf,     &
                                        etaphi_m,beta_glf,zpti_save,zptij_p,     &
                                        zptej_p,chie_dv,chii_dv,chii_save,       &
                                        grho1_glf,chitorot_dv,betae_glf,zeff_exp,&
                                        elong_glf,vphim_glf,betai_glf,           &
                                        gamma_p_exp,anfreq_m,anrate_m,glf_grho2, &
                                        egamma_m,egamma_exp,angrotp_exp,         &
                                        grho2_glf,ns_glf,itport_pt,vparm_glf,    &
                                        csda_glf,vperm_glf,rhosda_glf,vpolm_glf, &
                                        drhodrrrho,bteff_glf,ve_glf,tau_glf,     &
                                        ak1glf,ak2glf,aneoglf,fcglf,             &
                                        egamma_mdv,anfreq_mdv,anrate_mdv,        &
                                        diff_m,exch_m,gamma_p_m,anfreq2_m,       &
                                        etapar_m,egamma_d,anrate2_m,etaper_m,    &
                                        iglf_chie,iglf_chii,iglf_chiv,iglf_d,    &
                                        anfreq2_mdv,anrate2_mdv,                 &
                                        gamma_p_mdv,exch_mdv,etaper_mdv,         &
                                        etapar_mdv,etaphi_mdv,chii_mdv,          &
                                        chie_mdv,diff_mdv,itport_glf,            &
                                        i_delay,itene_dv,itti_dv,itangrot_dv,    &
                                        itenp_dv, itte_dv,dcoefp_glf23,          &
                                        jroot_glf,x_alpha_glf,init_egamma_exp,   &
                                        egamma_exp_initial,limp_glf,glf23_iglf,  &
                                        jeigen_glf,ibtflag_glf,irotstab,         &
                                        exbmult_glf,pdelay,etaphi_avg,           &
                                        chi_e_avg, chi_i_avg,exchmdv,navg_chi,   &
                                        chie_ms,chii_ms,first_step,              &
                                        freeze_alpha_exp,chie_avg,chii_avg,      &
                                        etaphi_ms,d_avg,xkang_avg,xchie_sv,      &
                                        xchii_sv,xkangwlsv,xke_glf23,xki_glf23,  &
                                        xkw_glf23,dv_method,spline_gsc,dv_delt,  &
                                        reset_mask,use_mask,alpha_exp_save,      &
                                        t_delay_glf23,diff_ms,diff_avg,          &
                                        niflux_gf, nimpflux_gf,tiflux_gf,        &
                                        teflux_gf,vparflux_gf,vperflux_gf,       &
                                        vphiflux_gf,etgflux_gf,                  &
                                        niflux_gf_m, nimpflux_gf_m,tiflux_gf_m,  &
                                        teflux_gf_m,vparflux_gf_m,vperflux_gf_m, &
                                        vphiflux_gf_m,etgflux_gf_m,              &
                                        niflux_gf_m_loc, nimpflux_gf_m_loc,      &
                                        tiflux_gf_m_loc,teflux_gf_m_loc,         &
                                        vparflux_gf_m_loc,vperflux_gf_m_loc,     &
                                        vphiflux_gf_m_loc,etgflux_gf_m_loc,d_sv, &
                                        gamma_net_i,gamma_net_e,gamma_net_i_loc, &
                                        gamma_net_e_loc,glf_p_flux,glf_e_flux,   &
                                        glf_m_flux,glf_send_packed,              &
                                        nspecies_glf,glf_val_pert,glf_val_base,  &
                                        glf_typ_val,glf_max_pack,vexb_base,      &
                                        glf_receive_packed,dqrm_drho,            &
                                        dvexb_base_drho,vexb_term1,ddiamag_drho, &
                                        dvexb_term1_drho,dvpar_term1_drho,       &
                                        diamag_term,vparm_term1, vperm_term1,    &
                                        dvphim_drho,egamma_exp_zc,gamma_p_exp_zc,&
                                        glf_corot,glf_gamma_net_i_output,        &
                                        glf_gamma_net_e_output,                  & 
                                        glf_anrate_output,                       &
                                        glf_anrate2_output,                      &
                                        glf_anfreq_output,                       &
                                        glf_anfreq2_output,                      &
                                        glf_anrate,glf_anrate2,                  &
                                        glf_anfreq,glf_anfreq2,                  &
                                        glf_etg_output,glf_eff_chi_output,         &
                                        rdrho


     USE tension_spline,         ONLY :  s_coef,dsp,d2sp,sp,s_bc,bc_0,bc_1,rho_w,&
                                         get_scoef,spline_intrp
     USE common_constants,       ONLY :  Permeability,joupkev,Electron_charge,   &
                                         Proton_Mass,zeroc,izero

     USE neutral_beams,          ONLY : enbeam_tot,nbion

     USE grid_class,             ONLY : nj,rminor_r,ravg_r,r,dr

     USE io_gcnmp,               ONLY : nlog

     USE ions_gcnmp,             ONLY : nprim,nion,z,zeff,atw,dmassden,          &
                                        tot_primary_ion_den,tot_thermal_ion_den, &
                                        fi_index,nimp
                                       
     USE plasma_properties,      ONLY : dischg,mhd_dat,dcoefp,xdchitot,xketot,   &
                                        xchietot,xchiitot,xkangrot,xkitot
                                      
     USE solcon_gcnmp,                 ONLY : freeze_xwr,freeze_xti,freeze_xte,  &
                                        use_avg_glf_chi,itran,icalln,use_glf23,  &
                                        cm_xke,cm_xki,cm_xkw,freeze_xni,         &
                                        use_glf23_flux,ntot,select_solver,dudr

     USE curden_terms,           ONLY : q
                                       
     USE vector_class,           ONLY : get_values

! ----------------------The following are local to sub glf_driver_single -----
! ----------------------declared here so info gets saved between calls -------
      
      IMPLICIT NONE

      REAL(DP) etapermdv,etaparmdv,etaphimdv,chiitimdv,diffnemdv,                &
               zptej_glf,zptij_glf, rmajor_glf,                                  &
               btor_exp,arho_iglf,etaphim, etaparm,                              &
               etaperm, exchm,zpnej_glf,zpnij_glf,chietemdv,cparam,cparam1,      &
               rdrho_local_p1,chiitim,chietem, vstar_sign,          &
               corot,pgeo_local,fc,akappa1, akappa2, alpha_neo,sum,ane,ani,      &
               dqdrho,drhodq,omci,tea,tia,btsq,u0,tmu0,tiny,tot_den,             &
               diffnem,amassimp_glf,zimp_glf,amassgas_glf,xmassp,diff,charge,    &
               relaxd,enea,ensum,enav,dbetadrho



      INTEGER(i4B) iptrb,idengrad_glf,                                           &
                   njm1,jr,jstart,jend,jsteps,jl,jj,j,k,irotstab_save,           &
                   jshoot,jmm,jmaxm,igrad, no_w_convection,                      &
                   iptrpa,i,jelc_clamp,jion_clamp

! ----------------------------------------------------------------------------




     CONTAINS

      SUBROUTINE  create_glf_output
! ----------------------------------------------------------------------------
! --  load some glf output profiles to be written to state file
! ----------------------------------------------------------------------------



 
        glf_gamma_net_e_output(:)     = gamma_net_e(:)
        glf_gamma_net_i_output(:)     = gamma_net_i(:)
 ! following are not currently set when using diffusive model
 ! so compensate here (need to fix this) HSJ 11/5/12
        IF(ALLOCATED(glf_anrate))THEN
           glf_anrate_output(:)          = glf_anrate(:)
           glf_anrate2_output(:)         = glf_anrate2(:)
           glf_anfreq_output(:)          = glf_anfreq(:)
           glf_anfreq2_output(:)         = glf_anfreq2(:)
           glf_etg_output(:)             = etgflux_gf_m(:)
           CALL set_glf_eff_chi
        ELSE
           glf_anrate_output(:)          = zeroc
           glf_anrate2_output(:)         = zeroc
           glf_anfreq_output(:)          = zeroc
           glf_anfreq2_output(:)         = zeroc
           glf_etg_output(:)             = zeroc
           glf_eff_chi_output(:,:)         = zeroc
        ENDIF
       END SUBROUTINE  create_glf_output



       SUBROUTINE gamma_ptrb(pertrb,jmm,stepsize,omega_fg,gamma_e_zc,gamma_p_zc)
! ----------------------------------------------------------------------------
! -- glf_egamma_exp depends on toroidal rotation and its gradient.
! -- (glf_gamma_p_exp depends only on the gradient )
! -- Here we perturb these field  values  and recalculate  glf_egamma_exp,glf_gamma_p_exp
! -- INPUT
! --     jmm
! --     stepsize
! --     omega_fg(0:nj-1), (gamma_ptrb assumes diamagnetic correction 
! --                        is included in omega_fg)
! --
! -- OUTPUT
! --    gamma_e(0:nj-1)
! --    gamma_p(0:nj-1)
! ----------------------------------------------------------------------------
           USE glf23_gcnmp,                    ONLY : glf_vparm, vexb_base,vexb_term1
           
           IMPLICIT NONE
           INTEGER(I4b)  jmm,pertrb
           REAL(DP) omega_fg(0:nj-1),gamma_e_zc(0:nj-2),gamma_p_zc(0:nj-2)

           REAL(DP) vexb_term1a,vphima,stepsize,qa,rmina,vexba,rhoa,drhodrrra,  &
                    vphim_prtrb,delta1, grad_omega,csdaa,omega,dvphim_prtrb,    &
                    vpar_base,vparmod,gpe,gpp , diamag_terma  
           REAL(DP) vexb(0:nj-1) ! locall array

           !---------------------------------------------------------------------------
           ! note:  gcnmp grid defined quantites:
           ! rho_glf(1:nj),rmin_glf(1:nj),q_glf(1:nj),csda_glf(1:nj),vexb_base(1:nj),
           ! ve_glf(1:nj),drhodrrrho(1:nj-1)
           !---------------------------------------------------------------------------



           vexb_term1a      =  0.5_DP*(vexb_term1(jmm+1)  + vexb_term1(jmm+2))
           vexb(0:nj-1)     =  ve_glf(1:nj) 
           vexba            =  0.5_DP*(vexb(jmm+1)      +  vexb(jmm)) 
           csdaa            =  0.5_DP*(csda_glf(jmm+2)  +  csda_glf(jmm+1))
           omega            =  0.5_DP*(omega_fg(jmm+1)  +  omega_fg(jmm))
           qa               =  0.5_DP*(q_glf(jmm+1)     +  q_glf(jmm+2))
           rmina            =  0.5_DP*(rmin_glf(jmm+1)  +  rmin_glf(jmm+2))
           vphima           =  0.5_DP*(vphim_glf(jmm+1) + vphim_glf(jmm+2))
           rhoa             =  0.5_DP*(rho_glf(jmm+2)   + rho_glf(jmm+1))
           diamag_terma     =  0.5_DP*(diamag_term(jmm+1)+diamag_term(jmm+2))
           drhodrrra        =  drhodrrrho(jmm+1)
           delta1           =  rhoa*drhodrrra/rmina/csdaa
           vpar_base        =  rmajor_glf*omega_fg(jmm)


           IF(pertrb == 5)THEN  
              !-----------------------------------------------------------------------  
              ! perturb the field value omega , get perturbed gamma_e,
              ! gamma_p depends only on grad of 
              ! omega so no perturbation for it here  but vpar is perturbed
              !-----------------------------------------------------------------------
              omega_fg(jmm)    = omega_fg(jmm)+ stepsize 
              vparmod          = rmajor_glf*0.5_DP*(omega_fg(jmm)+omega_fg(jmm+1))
              glf_vparm(jmm)   = glf_vparm(jmm) - vpar_base  ! remove base value
              glf_vparm(jmm)   = glf_vparm(jmm) + vparmod    ! add perturbed value
  
              vphim_prtrb      = vparmod + rmajor_glf*diamag_terma
              gamma_e_zc(jmm)  = delta1 *(dvexb_base_drho(jmm)                         &
                 + vexb_term1a*dvphim_drho(jmm) + vphim_prtrb*dvexb_term1_drho(jmm)    &
                 + (vexba*rmina/qa)*dqrm_drho(jmm))




           ELSEIF(pertrb == -5)THEN                       ! perturb gradient  omega


              grad_omega       = (omega_fg(jmm+1)-omega_fg(jmm))/(rho_glf(jmm+2) - rho_glf(jmm+1))
              dvphim_prtrb     = rmajor_glf*(grad_omega + stepsize)


              gamma_e_zc(jmm)  = delta1 *(dvexb_base_drho(jmm)                         &
                 + vexb_term1a*(dvphim_prtrb + rmajor_glf*ddiamag_drho(jmm)) + vphima*dvexb_term1_drho(jmm)    &
                 + (vexba*rmina/qa)*dqrm_drho(jmm))
              gamma_p_zc(jmm)  = -delta1 * (dvpar_term1_drho(jmm)   &
                 +  dvphim_prtrb + ddiamag_drho(jmm))

           ENDIF

  
 

         RETURN

       END SUBROUTINE  gamma_ptrb




       SUBROUTINE glf_driver_single(grid_pt,glf_setup)
!-------------------------------------------------------------------------------
! -- returns both d's and chi's as well as fluxes directly
! -- actual choice as to which will be used in the calculations
! -- is done in module fluxx
!
! -- Fluxes are returned in:
! --                niflux_gf_m
! --                nimpflux_gf_m
! --                tiflux_gf_m
! --                teflux_gf_m
! --                vparflux_gf_m
! --                vperflux_gf_m
! --                vphiflux_gf_m
! --                etgflux_gf_m
! -- diffusion coefficients  are returned in:
! --                dcoefp_glf23  (glf23 part only)
! --                dcoefp  (total, including dcoefp_glf23)
! --                xdchitot(j)          = dcoefp(1,1,j)
! --                xketot(j)            = dcoefp(nion+1,nion+1,j)
! --                xchietot(j)          = xketot(j)/enea
! --                xkitot(j)            = dcoefp(nion+2,nion+2,j)
! --                xchiitot(j)          = xkitot(j)/ensum
! --                xkangrot(j)          = dcoefp(nion+dp4,nion+dp4,j)
!-------------------------------------------------------------------------------

         INTEGER(I4B) grid_pt , glf_setup,k


 
         IF(glf_setup ==izero) CALL glf_init(glf_setup)

!         jmm   = grid_pt-1
         jmm   = grid_pt         
         jmaxm = nj-1               !grid size (0 to jmaxm in glf)


         !
         NOFREEZE: IF(freeze_xte .EQ. 0 .OR. freeze_xti .EQ. 0 &
              .OR. freeze_xwr .EQ. 0 .OR. freeze_xni .EQ. 0 ) THEN
            i_delay = 0                    ! no delay in ExB egamma_d
   
               !arrays dimensioned (0:jmaxp) in callglf2d are passed
               !by address and hence can be numbered (1:jmaxp+1) = (1,nj)
               !outside the glf routines. This fact is used explicitely below

               CALL callglf2d(                                                & !inputs
                    jeigen_glf,jroot_glf, glf23_iglf,                         &
                    jshoot, jmm, jmaxm, itport_pt  ,                          &
                    irotstab, te_glf, ti_glf, ne_glf,                         &
                    ni_glf, ns_glf,igrad, idengrad_glf,zpte_glf(jmm),         &
                    zpti_glf(jmm), zpne_glf(jmm), zpni_glf(jmm),              &
                    angrotp_exp, egamma_exp, gamma_p_exp,                     &
                    vphim_glf, vparm_glf, vperm_glf,                          &
                    zeff_exp, btor_exp, ibtflag_glf, rho_glf,                 &
                    arho_iglf, grho1_glf, grho2_glf,                          &
                    rmin_glf, rmaj_glf, rmajor_glf, zimp_glf,                 &
                    amassimp_glf, q_glf, shat_exp, alpha_exp,                 &
                    elong_glf, amassgas_glf,exbmult_glf,                      &
                    x_alpha_glf, i_delay,                                     &
                    diffnem, chietem,chiitim, etaphim, etaparm,               & ! outputs
                    etaperm,exchm, diff_m, chie_m, chii_m,                    &
                    etaphi_m, etapar_m, etaper_m,                             &
                    exch_m, egamma_m, egamma_d, gamma_p_m,                    &
                    anrate_m, anrate2_m, anfreq_m, anfreq2_m,                 &
                    niflux_gf, niflux_gf_m, nimpflux_gf,nimpflux_gf_m,        &
                    tiflux_gf,tiflux_gf_m,teflux_gf,teflux_gf_m,              &     
                    vparflux_gf,vparflux_gf_m,vperflux_gf,vperflux_gf_m,      &
                    vphiflux_gf,vphiflux_gf_m,etgflux_gf,etgflux_gf_m,        &
                    gamma_net_i,gamma_net_e  )

               j= grid_pt

  
 


               vphiflux_gf_m(j) = vphiflux_gf_m(j) * dmassden(j+1)*dischg%rmag &   ! dmassden(1:nj-1)
                                           *0.5_DP*(grho2_glf(j+1)+ grho2_glf(j))  ! kg/sec^2
               vparflux_gf_m(j) = vparflux_gf_m(j) * dmassden(j+1)*dischg%rmag &
                                           *0.5_DP*(grho2_glf(j+1)+ grho2_glf(j))  ! kg/sec^2
               vperflux_gf_m(j) = vperflux_gf_m(j) * dmassden(j+1)*dischg%rmag &
                                           *0.5_DP*(grho2_glf(j+1)+ grho2_glf(j))  ! kg/sec^2

               IF(chie_m(j) .LT. 0.0_DP)THEN
                  ve_eff_m(j)  = chie_m(j)*zpte_glf(j)/arho_iglf             
                  ve_eff_m(j)  = ve_eff_m(j)*100._DP
                  chie_m(j)    = 0.0_DP
               ELSE
                  ve_eff_m(j) = 0.0_DP
               ENDIF
               IF(chii_m(j) .LT. 0.0_DP)THEN
                  vi_eff_m(j) = chii_m(j)*zpti_glf(j)/arho_iglf 
                  vi_eff_m(j)  = vi_eff_m(j)*100._DP
                  chii_m(j) = 0.0_DP
               ELSE
                  vi_eff_m(j) = 0.0_DP
               ENDIF
               IF(etaphi_m(j) .LT. 0.0_DP)THEN
                  IF (no_w_convection .EQ. 1)THEN
                     etaphi_m(j) = 0.0_DP
                     vrot_eff_m(j) =0.0_DP
                  ELSE
                     vrot_eff_m(j) = etaphi_m(j)*zpti_glf(j)/arho_iglf !m/sec
                     etaphi_m(j) = 0.0_DP
                  ENDIF
               ELSE
                  vrot_eff_m(j) = 0.0_DP  !flux for vrot not impl
               ENDIF

               vp_eff_m(j) = 0.0_DP
               IF(use_glf23_flux(1) > 0)diff_m(j) = zeroc
               IF(diff_m(j) .LE. zeroc)THEN
                  vp_eff_m(j)   = diff_m(j)*zpni_glf(j)/arho_iglf  ! m/sec
                  diff_m(j)     = 0.0_DP
               ENDIF


         ENDIF NOFREEZE


!--------------------------------------------------------------------------
! -- load some quantities that may (or may not) be used by all solvers:
!--------------------------------------------------------------------------

            navg_chi = navg_chi+1
            IF(freeze_xte .EQ. 0)THEN
               chie_ms(grid_pt) = chie_m(grid_pt)  !save for possible reuse next ite
               chie_avg(grid_pt) = (chie_avg(grid_pt)*(navg_chi-1) &
                    +chie_m(grid_pt))/navg_chi
            ELSE
               IF(use_avg_glf_chi .EQ. 1)THEN
                  chie_m(grid_pt) = chie_avg(grid_pt)  
               ELSE
                  chie_m(grid_pt) = chie_ms(grid_pt)              !set to saved value
               ENDIF
            ENDIF
            IF(freeze_xti .EQ. 0)THEN
               chii_ms(grid_pt) = chii_m(grid_pt)  !save copies for possible reuse 
               chii_avg(grid_pt) = (chii_avg(grid_pt)*(navg_chi-1) &
                    +chii_m(grid_pt))/navg_chi
            ELSE
               IF(use_avg_glf_chi .EQ. 1)THEN
                  chii_m(grid_pt) = chii_avg(grid_pt)
               ELSE
                  chii_m(grid_pt) = chii_ms(grid_pt)      !set to saved value
               ENDIF
            ENDIF

            IF(freeze_xwr .EQ. 0)THEN
               etaphi_ms(grid_pt) = etaphi_m(grid_pt)
               etaphi_avg(grid_pt) = (etaphi_avg(grid_pt)*(navg_chi-1) &
                    +etaphi_m(grid_pt))/navg_chi
            ELSE
               IF(use_avg_glf_chi .EQ. 1)THEN
                  etaphi_m(grid_pt) = etaphi_avg(grid_pt)
               ELSE
                  etaphi_m(grid_pt) = etaphi_ms(grid_pt)        !set to saved value
               ENDIF
            ENDIF


            IF(freeze_xni .EQ. 0)THEN
               diff_ms(grid_pt) = diff_m(grid_pt)
               diff_avg(grid_pt) = (diff_avg(grid_pt)*(navg_chi-1) &
                    +diff_m(grid_pt))/navg_chi
            ELSE
               IF(use_avg_glf_chi .EQ. 1)THEN
                  diff_m(grid_pt) = diff_avg(grid_pt)
               ELSE
                  diff_m(grid_pt) = diff_ms(grid_pt)        !set to saved value
               ENDIF
            ENDIF
 

            !
            ! --- cparam is continuation PARAMETER, normally =1.0 HSJ  ----------`
            ! 
            cparam1 =cparam
            j = grid_pt
            IF(chii_m(j) .GT. 0.0)THEN
               chi_i_avg(j) = cparam1*chii_m(j)
               !chi_i_avg(j) = MIN(chi_i_avg(j),1.d6) !limit max val here
            ELSE
               chi_i_avg(j) = 0.0
            ENDIF


            IF(chie_m(j) .GT. 0.0)THEN
               chi_e_avg(j) =  cparam1*chie_m(j)
               !chi_e_avg(j) = MIN(chi_e_avg(j),1.d6) !limit max val here
            ELSE
               chi_e_avg(j) = 0.0
            ENDIF




            !
            !              get results into d and xkang for simulation purposes:
            !
            d_avg(j)     = diff_m(j)

            xkang_avg(j) = etaphi_m(j)




            !
            ! --- relax just the anomalous part of the diffusivity (i.e., without the
            ! --- neoclassical term) note that at this point chi_e_avg and chi_i_avg
            ! --- are just the anomalous part returned by SUBROUTINE callglf2d.
            ! --- IF icalln > 1, THEN xchie_sv and xchii_sv are
            ! --- just the anomalous parts from the previous iteration.
            !
            IF (icalln  .GT. 1) THEN
               d_avg(j)     =(1.0_DP - relaxd) *d_sv(j)+ &
                    relaxd *d_avg(j)
               chi_e_avg(j) =(1.0_DP - relaxd) *xchie_sv(j)+ &
                    relaxd *chi_e_avg(j)
               chi_i_avg(j) =(1.0_DP - relaxd) *xchii_sv(j)+ &
                    relaxd *chi_i_avg(j)
               xkang_avg(j) =(1.0_DP - relaxd) *xkangwlsv(j) + &
                    relaxd *xkang_avg(j)
            END IF
            !
            d_sv(j)       = diff_m(j)
            xchie_sv(j)   = chi_e_avg(j)
            xchii_sv(j)   = chi_i_avg(j)
            xkangwlsv(j)  = xkang_avg(j)
            !
            enea=0.5*(ene(j+1)+ene(j))
            ensum = 0.0
            DO i=1,nion             
               enav = 0.5 * (en(j,i)+en(j+1,i))     ! note that this must use previously set
                                                    ! en because of domain decomp in effect if
                                                    ! more than 1 cpu is used
               ensum = enav +ensum
            END DO
            !

            !The following cm_**      are active if the user 
            !has set use_mask > 0 
            !If a transport coefficient is 0.0  on the
            !initial iteration of a time step then
            !do not let it become active on subsequqnt 
            !iterations of the same time step:
            IF(use_mask > 0 .AND. reset_mask == 1) THEN
               IF(xke_glf23(j) .LE. 1.e-8_DP) cm_xke(j) = zeroc
               IF(xki_glf23(j) .LE. 1.e-8_DP) cm_xki(j) = zeroc
               IF(xkw_glf23(j) .LE. 1.e-8_DP) cm_xkw(j) = zeroc
            ENDIF

            xke_glf23(j)           = chi_e_avg(j)*0.5*(ene(j+1) + ene(j))
            !              add < (grad rho)**2 >  non circular correction factor

            xke_glf23(j) = iglf_chie*grho2_glf(j)*xke_glf23(j)*cm_xke(j)

            dcoefp(nion+1,nion+1,j)   = dcoefp(nion+1,nion+1,j) + xke_glf23(j) !TE
            dcoefp_glf23(nion+1,j)    = xke_glf23(j)

            xki_glf23(j)   = chi_i_avg(j)*0.5*( &
                 tot_primary_ion_den(j) + tot_primary_ion_den(j+1))
            xki_glf23(j)              = iglf_chii*grho2_glf(j)*xki_glf23(j)*cm_xki(j)
 
            dcoefp(nion+2,nion+2,j)   = dcoefp(nion+2,nion+2,j) + xki_glf23(j)  !TI
            dcoefp_glf23(nion+2,j)    = xki_glf23(j)

            xkw_glf23(j)  =  xkang_avg(j)*dmassden(j)*dischg%rmag**2 !(m^2/sec)*(kg/m^3)*m^2
            xkw_glf23(j)  =  iglf_chiv*grho2_glf(j)*xkw_glf23(j)*cm_xkw(j)
            dcoefp(nion+dp4,nion+dp4,j)   = dcoefp(nion+dp4,nion+dp4,j)+xkw_glf23(j) !ROTATION
            dcoefp_glf23(nion+dp4,j)      = xkw_glf23(j)                 ! kg m/sec

            ! ------------------------------------------------------------------------------------
            ! -- particle D's  here
            ! -- density, diff_m from glf is same for all ions
            ! -- if use_glf23_flux(1) > 0 then contribution from  diff_m is not to be used
            ! -- but neoclassical part allready present in dcoeff is retained.
            ! ------------------------------------------------------------------------------------
            IF(use_glf23_flux(1) == 0)THEN
               DO i=1,nion
                  dcoefp(i,i,j)     = dcoefp(i,i,j)+iglf_d*grho2_glf(j)*d_avg(j) 
                  dcoefp_glf23(i,j) = iglf_d*grho2_glf(j)*d_avg(j) 
               ENDDO
            ELSE
               DO i=1,nion
                  IF(dudr(i,j) .NE. zeroc)THEN
                     dcoefp_glf23(i,j)  = -niflux_gf_m(j)/dudr(i,j) 
                  ELSE
                     dcoefp_glf23(i,j)  = zeroc
                  ENDIF
                  dcoefp(i,i,j)         =  dcoefp(i,i,j) + dcoefp_glf23(i,j)
               ENDDO
            ENDIF

 
            xdchitot(j)          = dcoefp(1,1,j)
            xketot(j)            = dcoefp(nion+1,nion+1,j)
            xchietot(j)          = xketot(j)/enea
            xkitot(j)            = dcoefp(nion+2,nion+2,j)
            xchiitot(j)          = xkitot(j)/ensum
            xkangrot(j)          = dcoefp(nion+dp4,nion+dp4,j)
            !



            IF(grid_pt .LE. jelc_clamp)THEN
               j = grid_pt
               dcoefp(nion+1,nion+1,j)  =  dcoefp(nion+1,nion+1,jelc_clamp)
               dcoefp_glf23(nion+1,j)   =  dcoefp_glf23(nion+1,jelc_clamp)
               xketot(j)                =  dcoefp(nion+1,nion+1,j)
               xchietot(j)              =  xchietot(jelc_clamp)
            ENDIF



            !
            IF( grid_pt .LE. jion_clamp)THEN
               j = grid_pt
               dcoefp(nion+2,nion+2,j) =  dcoefp(nion+2,nion+2,jion_clamp)
               dcoefp_glf23(nion+2,j)  =  dcoefp_glf23(nion+2,jion_clamp)
               xkitot(j)               =  dcoefp(nion+2,nion+2,j)
               xchiitot(j)             =  xchiitot(jion_clamp)
            ENDIF



            irotstab = irotstab_save
 
   END SUBROUTINE glf_driver_single

#ifdef GCNMP

       SUBROUTINE glf_driver_single_new(grid_pt,glf_setup)
!-------------------------------------------------------------------------------
! -- returns both d's and chi's as well as fluxes directly
! -- actual choice as to which will be used in the calculations
! -- is done in module fluxx
!
! -- Fluxes are returned in (thorugh glf_driver_pert):
! --                niflux_gf_m
! --                nimpflux_gf_m
! --                tiflux_gf_m
! --                teflux_gf_m
! --                vparflux_gf_m
! --                vperflux_gf_m
! --                vphiflux_gf_m
! --                etgflux_gf_m
! --                glf_p_flux(:,1:2,0) = niflux,nimpflux
! --                glf_e_flux(:,1:2,0) = teflux_gf_m,tiflux_gf_m
! --                glf_m_flux(:,2:2,0) = vphiflux_gf_m
! -- diffusion coefficients  are returned in:
! --                dcoefp_glf23  (glf23 part only)
! --                dcoefp  (total, including dcoefp_glf23)
! --                xdchitot(j)          = dcoefp(1,1,j)
! --                xketot(j)            = dcoefp(nion+1,nion+1,j)
! --                xchietot(j)          = xketot(j)/enea
! --                xkitot(j)            = dcoefp(nion+2,nion+2,j)
! --                xchiitot(j)          = xkitot(j)/ensum
! --                xkangrot(j)          = dcoefp(nion+dp4,nion+dp4,j)
!-------------------------------------------------------------------------------

         INTEGER(I4B) grid_pt , glf_setup,pertrb


!         write(173,FMT='("in glf_driver_single_new,grid_point,glf_setup,pertrb =",x,i5,x,i5,x,i5)')grid_pt,glf_setup,pertrb ! 888888899999
         IF(glf_setup ==izero) CALL glf_init(glf_setup) ! set glf params at start of each time step

         jmm      = grid_pt-1          ! 0 to nj-1        
         jmaxm    = nj-1               ! grid size (0 to jmaxm in glf)
         pertrb   = -100               ! turn off perturbation calcualtions
         j= grid_pt
         !
         NOFREEZE: IF(freeze_xte .EQ. 0 .OR. freeze_xti .EQ. 0 &
              .OR. freeze_xwr .EQ. 0 .OR. freeze_xni .EQ. 0 ) THEN
            i_delay = 0                    ! no delay in ExB egamma_d


               CALL glf_driver_pert(jmm,pertrb)



         ENDIF NOFREEZE


!--------------------------------------------------------------------------
! -- load some quantities that may (or may not) be used by all solvers:
!--------------------------------------------------------------------------

            navg_chi = navg_chi+1
            IF(freeze_xte .EQ. 0)THEN
               chie_ms(grid_pt) = chie_m(grid_pt)  !save for possible reuse next ite
               chie_avg(grid_pt) = (chie_avg(grid_pt)*(navg_chi-1) &
                    +chie_m(grid_pt))/navg_chi
            ELSE
               IF(use_avg_glf_chi .EQ. 1)THEN
                  chie_m(grid_pt) = chie_avg(grid_pt)  
               ELSE
                  chie_m(grid_pt) = chie_ms(grid_pt)              !set to saved value
               ENDIF
            ENDIF
            IF(freeze_xti .EQ. 0)THEN
               chii_ms(grid_pt) = chii_m(grid_pt)  !save copies for possible reuse 
               chii_avg(grid_pt) = (chii_avg(grid_pt)*(navg_chi-1) &
                    +chii_m(grid_pt))/navg_chi
            ELSE
               IF(use_avg_glf_chi .EQ. 1)THEN
                  chii_m(grid_pt) = chii_avg(grid_pt)
               ELSE
                  chii_m(grid_pt) = chii_ms(grid_pt)      !set to saved value
               ENDIF
            ENDIF

            IF(freeze_xwr .EQ. 0)THEN
               etaphi_ms(grid_pt) = etaphi_m(grid_pt)
               etaphi_avg(grid_pt) = (etaphi_avg(grid_pt)*(navg_chi-1) &
                    +etaphi_m(grid_pt))/navg_chi
            ELSE
               IF(use_avg_glf_chi .EQ. 1)THEN
                  etaphi_m(grid_pt) = etaphi_avg(grid_pt)
               ELSE
                  etaphi_m(grid_pt) = etaphi_ms(grid_pt)        !set to saved value
               ENDIF
            ENDIF


            IF(freeze_xni .EQ. 0)THEN
               diff_ms(grid_pt) = diff_m(grid_pt)
               diff_avg(grid_pt) = (diff_avg(grid_pt)*(navg_chi-1) &
                    +diff_m(grid_pt))/navg_chi
            ELSE
               IF(use_avg_glf_chi .EQ. 1)THEN
                  diff_m(grid_pt) = diff_avg(grid_pt)
               ELSE
                  diff_m(grid_pt) = diff_ms(grid_pt)        !set to saved value
               ENDIF
            ENDIF


            !
            ! --- cparam is continuation PARAMETER, normally =1.0 HSJ  ----------`
            ! 
            cparam1 =cparam
            j = grid_pt
            IF(chii_m(j) .GT. 0.0)THEN
               chi_i_avg(j) = cparam1*chii_m(j)
               !chi_i_avg(j) = MIN(chi_i_avg(j),1.d6) !limit max val here
            ELSE
               chi_i_avg(j) = 0.0_DP
            ENDIF


            IF(chie_m(j) .GT. 0.0)THEN
               chi_e_avg(j) =  cparam1*chie_m(j)
               !chi_e_avg(j) = MIN(chi_e_avg(j),1.d6) !limit max val here
            ELSE
               chi_e_avg(j) = 0.0_DP
            ENDIF




            !
            !              get results into d and xkang for simulation purposes:
            !
            d_avg(j)     = diff_m(j)

            xkang_avg(j) = etaphi_m(j)




            !
            ! --- relax just the anomalous part of the diffusivity (i.e., without the
            ! --- neoclassical term) note that at this point chi_e_avg and chi_i_avg
            ! --- are just the anomalous part returned by SUBROUTINE callglf2d.
            ! --- IF icalln > 1, THEN xchie_sv and xchii_sv are
            ! --- just the anomalous parts from the previous iteration.
            !
            IF (icalln  .GT. 1) THEN
               d_avg(j)     =(1.0_DP - relaxd) *d_sv(j)+ &
                    relaxd *d_avg(j)
               chi_e_avg(j) =(1.0_DP - relaxd) *xchie_sv(j)+ &
                    relaxd *chi_e_avg(j)
               chi_i_avg(j) =(1.0_DP - relaxd) *xchii_sv(j)+ &
                    relaxd *chi_i_avg(j)
               xkang_avg(j) =(1.0_DP - relaxd) *xkangwlsv(j) + &
                    relaxd *xkang_avg(j)
            END IF
            !
            d_sv(j)       = diff_m(j)
            xchie_sv(j)   = chi_e_avg(j)
            xchii_sv(j)   = chi_i_avg(j)
            xkangwlsv(j)  = xkang_avg(j)
            !
            enea=0.5*(ene(j+1)+ene(j))
            ensum = 0.0
            DO i=1,nion             
               enav = 0.5 * (en(j,i)+en(j+1,i))     ! note that this must use previously set
                                                    ! en because of domain decomp in effect if
                                                    ! more than 1 cpu is used
               ensum = enav +ensum
            END DO
            !

            !The following cm_**      are active if the user 
            !has set use_mask > 0 
            !If a transport coefficient is 0.0  on the
            !initial iteration of a time step then
            !do not let it become active on subsequqnt 
            !iterations of the same time step:
            IF(use_mask > 0 .AND. reset_mask == 1) THEN
               IF(xke_glf23(j) .LE. 1.e-8_DP) cm_xke(j) = zeroc
               IF(xki_glf23(j) .LE. 1.e-8_DP) cm_xki(j) = zeroc
               IF(xkw_glf23(j) .LE. 1.e-8_DP) cm_xkw(j) = zeroc
            ENDIF

            xke_glf23(j)           = chi_e_avg(j)*0.5*(ene(j+1) + ene(j))
            !              add < (grad rho)**2 >  non circular correction factor

            xke_glf23(j) = iglf_chie*grho2_glf(j)*xke_glf23(j)*cm_xke(j)

            dcoefp(nion+1,nion+1,j)   = dcoefp(nion+1,nion+1,j) + xke_glf23(j) !TE
            dcoefp_glf23(nion+1,j)    = xke_glf23(j)

            xki_glf23(j)   = chi_i_avg(j)*0.5*( &
                 tot_primary_ion_den(j) + tot_primary_ion_den(j+1))
            xki_glf23(j)              = iglf_chii*grho2_glf(j)*xki_glf23(j)*cm_xki(j)
            dcoefp(nion+2,nion+2,j)   = dcoefp(nion+2,nion+2,j) + xki_glf23(j)  !TI
            dcoefp_glf23(nion+2,j)    = xki_glf23(j)

            xkw_glf23(j)  =  xkang_avg(j)*dmassden(j)*dischg%rmag**2 !(m^2/sec)*(kg/m^3)*m^2
            xkw_glf23(j)  =  iglf_chiv*grho2_glf(j)*xkw_glf23(j)*cm_xkw(j)
            dcoefp(nion+dp4,nion+dp4,j)   = dcoefp(nion+dp4,nion+dp4,j)+xkw_glf23(j) !ROTATION
            dcoefp_glf23(nion+dp4,j)      = xkw_glf23(j)                 ! kg m/sec

            ! ------------------------------------------------------------------------------------
            ! -- particle D's  here
            ! -- density, diff_m from glf is same for all ions
            ! -- if use_glf23_flux(1) > 0 then contribution from  diff_m is not to be used
            ! -- but neoclassical part allready present in dcoeff is retained.
            ! ------------------------------------------------------------------------------------
            IF(use_glf23_flux(1) == 0)THEN
               DO i=1,nion
                  dcoefp(i,i,j)     = dcoefp(i,i,j)+iglf_d*grho2_glf(j)*d_avg(j) 
                  dcoefp_glf23(i,j) = iglf_d*grho2_glf(j)*d_avg(j) 
               ENDDO
            ELSE
               DO i=1,nion
                  IF(dudr(i,j) .NE. zeroc)THEN
                     dcoefp_glf23(i,j)  = -niflux_gf_m(j)/dudr(i,j) 
                  ELSE
                     dcoefp_glf23(i,j)  = zeroc
                  ENDIF
                  dcoefp(i,i,j)         =  dcoefp(i,i,j) + dcoefp_glf23(i,j)
               ENDDO
            ENDIF
            ! lhs [0,nj-1] == rhs [1:nj]
            glf_p_flux(j-1,2,0)    = niflux_gf_m(j) 
            glf_p_flux(j-1,3,0)    = nimpflux_gf_m(j)
            glf_e_flux(j-1,1,0)    = teflux_gf_m(j)
            glf_e_flux(j-1,2,0)    = tiflux_gf_m(j)
            glf_m_flux(j-1,2,0)    = vphiflux_gf_m(j)
           
 
            xdchitot(j)          = dcoefp(1,1,j)
            xketot(j)            = dcoefp(nion+1,nion+1,j)
            xchietot(j)          = xketot(j)/enea
            xkitot(j)            = dcoefp(nion+2,nion+2,j)
            xchiitot(j)          = xkitot(j)/ensum
            xkangrot(j)          = dcoefp(nion+dp4,nion+dp4,j)
            !



            IF(grid_pt .LE. jelc_clamp)THEN
               j = grid_pt
               dcoefp(nion+1,nion+1,j)  =  dcoefp(nion+1,nion+1,jelc_clamp)
               dcoefp_glf23(nion+1,j)   =  dcoefp_glf23(nion+1,jelc_clamp)
               xketot(j)                =  dcoefp(nion+1,nion+1,j)
               xchietot(j)              =  xchietot(jelc_clamp)
            ENDIF



            !
            IF( grid_pt .LE. jion_clamp)THEN
               j = grid_pt
               dcoefp(nion+2,nion+2,j) =  dcoefp(nion+2,nion+2,jion_clamp)
               dcoefp_glf23(nion+2,j)  =  dcoefp_glf23(nion+2,jion_clamp)
               xkitot(j)               =  dcoefp(nion+2,nion+2,j)
               xchiitot(j)             =  xchiitot(jion_clamp)
            ENDIF




            irotstab = irotstab_save

       RETURN
       
   END SUBROUTINE glf_driver_single_new







       SUBROUTINE glf_driver_pert(grid_pt,pertrb)
!-----------------------------------------------------------------------------------------------------
! -- use this driver for the perturbation approach. THIS is a flux based method
! -- so we enforce use_glf23_flux(i) = 1 for all i (which are in simulation mode)
! -- returns both d's and chi's as well as fluxes directly
! -- actual choice as to which will be used in the calculations
! -- is done in module fluxx. 
! -- Note transport flags for glf23:
! --        itport_pt(1)= itran(1)           ! density
! --        itport_pt(2)= itran(nion+1)      ! Te
! --        itport_pt(3)= itran(nion+2)      ! Ti
! --        itport_pt(4)= itran(nion+dp4)    ! tor rot 
! --        itport_pt(5)=0                   ! vtheta transport  ( not  currently done) 
! -- This suborutine is set up to run with itport_pt(4) =0,itport_pt(5) =0 OR
! -- itport_pt(4)= 1 and itport_pt(5) = 0  All other combinations of itport_pt(4,5) 
! -- are not useable at this time
! --
! -- Fluxes are returned in:
! --                niflux_gf_m
! --                nimpflux_gf_m
! --                tiflux_gf_m
! --                teflux_gf_m
! --                vparflux_gf_m
! --                vperflux_gf_m
! --                vphiflux_gf_m
! --                etgflux_gf_m
! -- diffusion coefficients  are returned in:
! --                dcoefp_glf23  (glf23 part only)
! --                dcoefp  (total, including dcoefp_glf23)
! --                xdchitot(j)          = dcoefp(1,1,j)
! --                xketot(j)            = dcoefp(nion+1,nion+1,j)
! --                xchietot(j)          = xketot(j)/enea
! --                xkitot(j)            = dcoefp(nion+2,nion+2,j)
! --                xchiitot(j)          = xkitot(j)/ensum
! --                xkangrot(j)          = dcoefp(nion+dp4,nion+dp4,j)
!----------------------------------------------------------------------------------------HSJ-----------

        USE glf23_gcnmp,                    ONLY :                                                  &
                                                   glf_te_fg,glf_ti_fg,glf_ne_fg,glf_ni_fg,         &
                                                   glf_ns_fg,glf_angrot_fg,glf_exch_m,              &
                                                   glf_egamma_exp,glf_gamma_p_exp,glf_vphim,        &
                                                   glf_vparm,glf_vperm,glf_zeff,glf_rho,glf_grho1,  &
                                                   glf_grho2,glf_rmin,glf_rmaj,glf_q,glf_shat,      &
                                                   glf_alpha,glf_elong,glf_diff_m,glf_chie_m,       &
                                                   glf_chii_m,glf_etaphi_m,glf_etapar_m,            &
                                                   glf_etaper_m,glf_egamma_m,glf_egamma_d,          &
                                                   glf_gamma_p_m,glf_anrate,glf_anrate2,glf_anfreq, &
                                                   glf_anfreq2,glf_ni_flux,glf_nimp_flux,           &
                                                   glf_ti_flux,glf_te_flux,glf_vparflux_gf_m,       &
                                                   glf_vperflux_gf_m,glf_vphiflux_gf_m,             &
                                                   glf_etgflux_gf_m,glf_gamma_net_i,glf_gamma_net_e,&
                                                   nspecies_glf,glf_egamma_zc,glf_gamma_p_zc

         USE error_handler,                 ONLY : iomaxerr, lerrno,terminate

         USE io_gcnmp,                      ONLY : nlog

         USE common_constants,              ONLY : izero,zeroc

         USE solcon_gcnmp,                  ONLY : itran_max,dudr,stepsize_factor

        USE plasma_properties,              ONLY : dischg

        USE ions_gcnmp,                     ONLY : dmassden
        
        USE grid_class,                     ONLY : r2capi

        !USE MPI_data,                      ONLY :  myid,master,mpiierr 

         IMPLICIT NONE

         REAL(DP) ani,tea,tia,ane,zpte_in,zpti_in,zpne_in,zpni_in,drho,sqrteta 
         REAL(DP) denfactor,ni_mod,ne_mod,vphi_mod,te_mod,ti_mod,stepsize
         REAL(DP) grad_angrot,gradvphim,gradnem,gradnim,gradtem,gradtim
   
         INTEGER(I4B) grid_pt , pertrb,jmm,sj,j
 
  


          ! for nwt_pred_core solver use irotstab = 0 to input egamma_exp and gamma_p_exp
          ! broken down according to field and field gradient dependences
          ! for maccormack solver this can be 0 or one as input in namelist:
          ! irtostab =0    changed 9/18/2012 HSJ
          ! IF(pertrb .GT. -50)irotstab = 0   


         jmaxm       = nj-1                   ! grid size (0 to jmaxm in glf)
         jmm         = grid_pt                ! in callglf2dcrct jmm ranges from 0 to jmaxm-2    
         i_delay     = izero                  ! no delay in ExB egamma_d

         denfactor   = 1.e-19_DP


         IF(itport_pt(4) .LT. izero .OR.  itport_pt(5) .NE. izero  )THEN
            lerrno   =  iomaxerr + 155                     
            CALL terminate(lerrno,nlog)
         ENDIF



 
 
         IF(ALLOCATED(glf_te_fg))THEN
            sj = SIZE(glf_te_fg)
            IF(sj .NE. nj)CALL glf_deallocate_lcl
         ENDIF
 
         IF(.NOT. ALLOCATED(glf_te_fg))CALL glf_allocate_lcl

 
         ! each time routine is called set current base values  (determined in glf_init)
         ! because pertrubations done below may have changed  some of these
         glf_te_fg(0:nj-1)         = te_glf(1:nj)
         glf_ti_fg(0:nj-1)         = ti_glf(1:nj)
         glf_ne_fg(0:nj-1)         = ne_glf(1:nj)          ! in units of 10e19
         glf_ni_fg(0:nj-1)         = ni_glf(1:nj)          ! in units of 10e19
         glf_ns_fg(0:nj-1)         = ns_glf(1:nj)          ! in units of 10e19
         glf_angrot_fg(0:nj-1)     = angrot(1:nj)          ! no diamagnetic corrected
         glf_egamma_exp(0:nj-1)    = egamma_exp(1:nj)  
         glf_egamma_zc(0:nj-2)     = egamma_exp_zc(0:nj-2) 
         glf_egamma_d(0:nj-1,1:10) = egamma_d(1:nj,1:10)
         glf_gamma_p_exp(0:nj-1)   = gamma_p_exp(1:nj)
         glf_gamma_p_zc(0:nj-2)    = gamma_p_exp_zc(0:nj-2) 
         glf_exch_m(0:nj-1)        = exch_m(1:nj)
         glf_zeff(0:nj-1)          = zeff_exp(1:nj)
         glf_rho(0:nj-1)           = rho_glf(1:nj)
         glf_grho1(0:nj-1)         = grho1_glf(1:nj)
         glf_grho2(0:nj-1)         = grho2_glf(1:nj)
         glf_rmin(0:nj-1)          = rmin_glf(1:nj)
         glf_rmaj(0:nj-1)          = rmaj_glf(1:nj)
         glf_q(0:nj-1)             = q_glf(1:nj)
         glf_shat(0:nj-1)          = shat_exp(1:nj)
         glf_alpha(0:nj-1)         = alpha_exp(1:nj)
         glf_elong(0:nj-1)         = elong_glf(1:nj) 

        !---------------------------------------------------------------------------- 
        ! the following 3 speeds are reset in callglf2d if itport_pt(4)=itport_pt(5) =0
        ! callglf2d output momentum chi and fluxes are zero in this case.
             glf_vphim(0:nj-1)         = vphim_glf(1:nj)
             glf_vparm(0:nj-1)         = vparm_glf(1:nj)
             glf_vperm(0:nj-1)         = vperm_glf(1:nj)
        ! if itport(4) = 1 and itport(5) = 0 then 
        ! callglf2d internally calculates vexb (=ve) and vpar using input  glf_vphim
        !                                 vpar using glf_vphim 
        !                                 glf_vparm is reset to vpar
        !                                 glf_vperm is reset to vexb + corrections
        !                                 note that irotstab =0 prevents recalclation
        !                                 of glf_egamma_exp and glf_gamma_p_exp however
        !---------------------------------------------------------------------------- 



!CHSJ            the CHSJ flag maches the one in callgfl2dcrct
!CHSJ    drho                  = glf_rho(jmm-1)-glf_rho(jmm)
!CHSJ    zpni_in               = -(dlog(glf_ni_fg(jmm-1))-dlog(glf_ni_fg(jmm)))/drho
!CHSJ    zpne_in               = -(dlog(glf_ne_fg(jmm-1))-dlog(glf_ne_fg(jmm)))/drho
!CHSJ    zpti_in               = -(dlog(glf_ti_fg(jmm-1))-dlog(glf_ti_fg(jmm)))/drho
!CHSJ    zpte_in               = -(dlog(glf_te_fg(jmm-1))-dlog(glf_te_fg(jmm)))/drho

         igrad = 1             ! always input the gradients to callglf2d  ,zpte_in,zpti_in
                               ! zpne_in,zpni_in             
         drho                  = glf_rho(jmm)-glf_rho(jmm+1)         ! rho_glf in (0.0,1.0]
         zpni_in               = -(dlog(glf_ni_fg(jmm))-dlog(glf_ni_fg(jmm+1)))/drho
         zpne_in               = -(dlog(glf_ne_fg(jmm))-dlog(glf_ne_fg(jmm+1)))/drho
         zpti_in               = -(dlog(glf_ti_fg(jmm))-dlog(glf_ti_fg(jmm+1)))/drho
         zpte_in               = -(dlog(glf_te_fg(jmm))-dlog(glf_te_fg(jmm+1)))/drho

    IF(pertrb .GT. -50)THEN        ! prtrb  .GE. -5 if perturbation solver in use
                                   ! otherwise skip over this section

         ! set the perturbation variables:
         sqrteta = 1.e-2_DP*stepsize_factor ! for 1.e-n perturbation  will be in digit no n
 
         IF(pertrb == 0)THEN
            glf_val_base(1,jmm)   = glf_ni_fg(jmm)/denfactor     ! undo glf scaling
            glf_val_base(2,jmm)   = glf_te_fg(jmm)     ! Kev
            glf_val_base(3,jmm)   = glf_ti_fg(jmm)     ! Kev
            glf_val_base(4,jmm)   = zeroc
            glf_val_base(5,jmm)   = glf_angrot_fg(jmm)      ! no diamagnetic contorbution
 
            gradnem               = (glf_ne_fg(jmm+1) - glf_ne_fg(jmm))/dr(jmm+1) ! 1/m, ( dr is in m) 
                                                                              ! densities in units of  1.e19
            gradnim               = (glf_ni_fg(jmm+1) - glf_ni_fg(jmm ))/dr(jmm+1)
            gradtem               = (glf_te_fg(jmm+1) - glf_te_fg(jmm))/dr(jmm+1)
            gradtim               = (glf_ti_fg(jmm+1) - glf_ti_fg(jmm))/dr(jmm+1)
            grad_angrot           = (glf_angrot_fg(jmm+1)- glf_angrot_fg(jmm))/dr(jmm+1) ! no diamagnetic contorbution/dr(jmm+1)
            glf_val_base(-1,jmm)  = gradnim/denfactor                         ! undo tglf scaling
            glf_val_base(-2,jmm)  = gradtem
            glf_val_base(-3,jmm)  = gradtim
            glf_val_base(-4,jmm)  = zeroc
            glf_val_base(-5,jmm)  = grad_angrot
         ENDIF

 
         stepsize              = sqrteta*MAX(ABS(glf_val_base(ABS(pertrb),jmm)),glf_typ_val(ABS(pertrb))) &
                                                    *SIGN(1.D0,glf_val_base(ABS(pertrb),jmm))

         glf_val_pert(pertrb,jmm) = stepsize ! if itran(i) =0 then glf_val_pert(pertrb,jmm)=0
                                             ! for all pertrb and jmm is set in pturb_vars


         IF( pertrb          == 1 )THEN
             ni_mod           = glf_ni_fg(jmm) + stepsize*denfactor
             glf_ni_fg(jmm)   = ni_mod

            
         ELSE IF( pertrb     == 2 )THEN
             te_mod           = glf_te_fg(jmm) + stepsize
             glf_te_fg(jmm)   = te_mod

         ELSE IF( pertrb     == 3 )THEN
             ti_mod           = glf_ti_fg(jmm) + stepsize
             glf_te_fg(jmm)   = ti_mod

         ELSE IF( pertrb     == 4 )THEN
             glf_val_pert(pertrb,jmm) = zeroc
         ELSE IF( pertrb     == 5 )THEN 

             CALL gamma_ptrb(pertrb,jmm,stepsize,glf_angrot_fg,                     &
                                             glf_egamma_zc,glf_gamma_p_zc)

            ! here we do not change glf_egamma_exp,glf_gamma_p_exp
            ! because these quantities depend only on the gradient
            ! of vphi ( = Ra*angrot, irotstab = 0 so callglf2d does not redefine these)
         ENDIF
    


 
         IF( jmm .LT. nj-1)THEN    
               IF( pertrb      == -1 )THEN
                  ni_mod        = glf_ni_fg(jmm) + stepsize*denfactor
                  !CHSJ          zpni_in       = -(dlog(glf_ni_fg(jmm-1))-dlog(ni_mod))/drho
                  zpni_in       = -(dlog(glf_ni_fg(jmm))-dlog(ni_mod))/drho

               ELSE IF( pertrb == -2 )THEN

                  !CHSJ                te_mod        = glf_te_fg(jmm) + stepsize
                  !CHSJ               zpte_in       = -(dlog(glf_te_fg(jmm-1))-dlog(te_mod))/drho
                  te_mod        = glf_te_fg(jmm+1) + stepsize
                  zpte_in       = -(dlog(glf_te_fg(jmm))-dlog(te_mod))/drho

                                                
               ELSE IF( pertrb == -3 )THEN
                  !CHSJ                ti_mod        = glf_ti_fg(jmm) + stepsize
                  !CHSJ                zpti_in       = -(dlog(glf_ti_fg(jmm-1))-dlog(ti_mod))/drho
                  ti_mod        = glf_ti_fg(jmm+1) + stepsize
                  zpti_in       = -(dlog(glf_ti_fg(jmm))-dlog(ti_mod))/drho
               ELSE IF( pertrb == -4 )THEN
                  glf_val_pert(pertrb,jmm)  = zeroc

               ELSE IF( pertrb == -5 )THEN             ! toroidal rotation gradient perturbation
                  ! must use stepsize ( so glf_val_per is set)
                  ! here we perturb egamma_exp,gammap_exp 
                  ! irotstab = 0 so callglf2d does not redefine these
                  CALL gamma_ptrb(pertrb,jmm,stepsize,glf_angrot_fg,                     &
                       glf_egamma_zc,glf_gamma_p_zc)
               ENDIF
          ENDIF         ! jmm .LT. nj-1
       ENDIF            ! pertrb > -50



          validjmm : IF( jmm .LT. nj-1)THEN  

!CHSJ          this is the corrected modified routine  for callglf2d: 


               CALL callglf2dcrct(                                                & !inputs
                    jeigen_glf,jroot_glf, glf23_iglf,                             &
                    jshoot, jmm, jmaxm, itport_pt  ,                              &
                    irotstab, glf_te_fg, glf_ti_fg, glf_ne_fg,                    &
                    glf_ni_fg,glf_ns_fg,igrad, idengrad_glf,zpte_in,              &
                    zpti_in, zpne_in, zpni_in,                                    &
                    glf_angrot_fg,glf_egamma_zc, glf_gamma_p_zc,                  &
                    glf_vphim, glf_vparm,glf_vperm,                               &
                    glf_zeff, btor_exp, ibtflag_glf, glf_rho,                     &
                    arho_iglf, glf_grho1, glf_grho2,                              &
                    glf_rmin, glf_rmaj, rmajor_glf, zimp_glf,                     &
                    amassimp_glf, glf_q, glf_shat, glf_alpha,                     &
                    glf_elong, amassgas_glf,exbmult_glf,                          &
                    x_alpha_glf, i_delay,                                         &
                    diffnem, chietem,chiitim, etaphim, etaparm,                   & ! outputs
                    etaperm,exchm, glf_diff_m, glf_chie_m, glf_chii_m,            &
                    glf_etaphi_m, glf_etapar_m, glf_etaper_m,                     &
                    glf_exch_m, glf_egamma_m, glf_egamma_d, glf_gamma_p_m,        &
                    glf_anrate, glf_anrate2, glf_anfreq, glf_anfreq2,             &
                    niflux_gf, glf_ni_flux, nimpflux_gf,glf_nimp_flux,            &
                    tiflux_gf,glf_ti_flux,teflux_gf,glf_te_flux,                  &    
                    vparflux_gf,glf_vparflux_gf_m,vperflux_gf,glf_vperflux_gf_m,  &
                    vphiflux_gf,glf_vphiflux_gf_m,etgflux_gf,                     &
                    glf_etgflux_gf_m,glf_gamma_net_i,glf_gamma_net_e  )

              CALL glf_limit_values (                                             &
                    jmm,jmaxm,                                                    &
                    diffnem, chietem,chiitim, etaphim, etaparm,                   &
                    etaperm,exchm, glf_diff_m, glf_chie_m, glf_chii_m,            &
                    glf_etaphi_m, glf_etapar_m, glf_etaper_m,                     &
                    glf_exch_m, glf_egamma_m, glf_egamma_d, glf_gamma_p_m,        &
                    glf_anrate, glf_anrate2, glf_anfreq, glf_anfreq2,             &
                    niflux_gf, glf_ni_flux, nimpflux_gf,glf_nimp_flux,            &
                    tiflux_gf,glf_ti_flux,teflux_gf,glf_te_flux,                  &    
                    vparflux_gf,glf_vparflux_gf_m,vperflux_gf,glf_vperflux_gf_m,  &
                    vphiflux_gf,glf_vphiflux_gf_m,etgflux_gf,                     &
                    glf_etgflux_gf_m,glf_gamma_net_i,glf_gamma_net_e  )

 
     !------------------------------------------------------------------------------------------
     ! handle conversion to toroidal momentum diffusivity
     ! for use in  nwt_pred_core we do not add diffusivities into the
     ! total d because the turbulent part is handled separately in terms of fluxes:
     !------------------------------------------------------------------------------------------
                   j= grid_pt           ! dmassden(1:nj-1),glf_grho2(0:nj-1), here j=0 to nj-2
                   glf_vphiflux_gf_m(j) = glf_vphiflux_gf_m(j) * dmassden(j+1)*dischg%rmag &
                                           *0.5_DP*(glf_grho2(j+1)+ glf_grho2(j))  ! kg/sec^2

                   glf_vparflux_gf_m(j) = glf_vparflux_gf_m(j) * dmassden(j+1)*dischg%rmag &
                                           *0.5_DP*(glf_grho2(j+1)+ glf_grho2(j))  ! kg/sec^2

                   glf_vperflux_gf_m(j) = glf_vperflux_gf_m(j) * dmassden(j+1)*dischg%rmag &
                                           *0.5_DP*(glf_grho2(j+1)+ glf_grho2(j))  ! kg/sec^2
                   vphiflux_gf_m(j+1) = glf_vphiflux_gf_m(j) ! kg/sec^2
                   vparflux_gf_m(j+1) = glf_vparflux_gf_m(j)
                   vperflux_gf_m(j+1) = glf_vperflux_gf_m(j)
                   xkw_glf23(j+1)     =  etaphim*dmassden(j+1)*dischg%rmag**2 ! (m^2/sec)*(kg/m^3)*m^2
                   xkw_glf23(j+1)     =  0.5_DP*(glf_grho2(j+1)+ glf_grho2(j))*xkw_glf23(j+1)
                   
      !            dcoefp(nion+dp4,nion+dp4,j) = dcoefp(nion+dp4,nion+dp4,j)+xkw_glf23(j+1) !ROTATION
                   dcoefp_glf23(nion+dp4,j+1)  = xkw_glf23(j+1)                  ! kg m/sec


 
     pertrb0_cond :            IF( pertrb == izero .OR. pertrb .LT. -50)THEN
                    j= grid_pt
                    diff_m(j+1)              = glf_diff_m(j)
                    chie_m(j+1)              = glf_chie_m(j)
                    chii_m(j+1)              = glf_chii_m(j)
                    etaphi_m(j+1)            = glf_etaphi_m(j)
                    etapar_m(j+1)            = glf_etapar_m(j)
                    etaper_m(j+1)            = glf_etaper_m(j)
                    egamma_m(j+1)            = glf_egamma_m(j)
                    egamma_d(j+1,1:10)       = glf_egamma_d(j,1:10)
                    gamma_p_m(j+1)           = glf_gamma_p_m(j)
                    anrate_m(j+1)            = glf_anrate(j)
                    anrate2_m(j+1)           = glf_anrate2(j)
                    anfreq_m(j+1)            = glf_anfreq(j)
                    anfreq2_m(j+1)           = glf_anfreq2(j)
                    niflux_gf_m(j+1)         = glf_ni_flux(j)
                    nimpflux_gf_m(j+1)       = glf_nimp_flux(j)
                    tiflux_gf_m(j+1)         = glf_ti_flux(j)
                    teflux_gf_m(j+1)         = glf_te_flux(j)
                    vparflux_gf_m(j+1)       = glf_vparflux_gf_m(j)
                    vperflux_gf_m(j+1)       = glf_vperflux_gf_m(j) 
                    vphiflux_gf_m(j+1)       = glf_vphiflux_gf_m(j)  
                    etgflux_gf_m(j+1)        = glf_etgflux_gf_m(j)
                    gamma_net_i(j+1)         = glf_gamma_net_i(j)
                    gamma_net_e(j+1)         = glf_gamma_net_e(j)

                    IF(chie_m(j+1) .LT. 0.0_DP)THEN
                       ve_eff_m(j+1)  = chie_m(j+1)*zpte_glf(j+1)/arho_iglf             
                       ve_eff_m(j+1)  = ve_eff_m(j+1)*100._DP
                       chie_m(j+1)    = 0.0_DP
                    ELSE
                       ve_eff_m(j+1) = 0.0_DP
                    ENDIF
                    IF(chii_m(j+1) .LT. 0.0_DP)THEN
                       vi_eff_m(j+1) = chii_m(j+1)*zpti_glf(j+1)/arho_iglf 
                       vi_eff_m(j+1)  = vi_eff_m(j+1)*100._DP
                       chii_m(j+1) = 0.0_DP
                    ELSE
                       vi_eff_m(j+1) = 0.0_DP
                    ENDIF
                    IF(etaphi_m(j+1) .LT. 0.0_DP)THEN
                       IF (no_w_convection .EQ. 1)THEN
                          etaphi_m(j+1) = 0.0_DP
                          vrot_eff_m(j+1) =0.0_DP
                       ELSE
                          vrot_eff_m(j+1) = etaphi_m(j+1)*zpti_glf(j+1)/arho_iglf !m/sec
                          etaphi_m(j+1) = 0.0_DP
                       ENDIF
                    ELSE
                       vrot_eff_m(j+1) = 0.0_DP  !flux for vrot not impl
                    ENDIF

                    vp_eff_m(j+1) = 0.0_DP
                    IF(use_glf23_flux(1) > 0)diff_m(j+1) = zeroc
                    IF(diff_m(j+1) .LE. zeroc)THEN
                       vp_eff_m(j+1)   = diff_m(j+1)*zpni_glf(j+1)/arho_iglf  ! m/sec
                       diff_m(j+1)     = 0.0_DP
                    ENDIF

                ENDIF pertrb0_cond

        ELSE           validjmm                   ! jmm >= nj-1 ,error condition
            lerrno = iomaxerr + 45                    
            CALL terminate(lerrno,nlog)
        ENDIF          validjmm

 
     RETURN
     END     SUBROUTINE glf_driver_pert

#endif


      SUBROUTINE glf_allocate_wrap(njnew)
!------------------------------------------------------------------------------------
! -- Argument version of glf_allocate:
!------------------------------------------------------------------------------------
         USE nrtype,                                 ONLY :  DP,I4B

         USE grid_class,                             ONLY :  nj

         IMPLICIT NONE
         INTEGER(I4B) njnew,njsave
         

         njsave = nj
         nj     = njnew
         CALL glf_allocate
         nj     = njsave


         RETURN
      END   SUBROUTINE glf_allocate_wrap



      SUBROUTINE allocate_glf_output(grid_size)
        USE nrtype,                                 ONLY :  DP,I4B

        USE glf23_gcnmp,                            ONLY :  glf_p_output,           &
                                                            glf_e_output,           &
                                                            glf_m_output,           &
                                                            glf_etg_output,         &
                                                            glf_anrate_output,      &
                                                            glf_anrate2_output,     &
                                                            glf_anfreq_output,      &
                                                            glf_anfreq2_output,     &
                                                            glf_gamma_net_i_output, &
                                                            glf_gamma_net_E_output, &
                                                            glf_eff_chi_output

        USE solcon_gcnmp,                           ONLY :  ntot

        IMPLICIT NONE
        INTEGER(I4B) grid_size
         




          IF(ALLOCATED(glf_p_output))DEALLOCATE(glf_p_output)
              ALLOCATE(glf_p_output(grid_size,3))
          IF(ALLOCATED(glf_e_output))DEALLOCATE(glf_e_output)
              ALLOCATE(glf_e_output(grid_size,3))
          IF(ALLOCATED(glf_m_output))DEALLOCATE(glf_m_output)
              ALLOCATE(glf_m_output(grid_size,3))
          IF(ALLOCATED(glf_etg_output))DEALLOCATE(glf_etg_output)
              ALLOCATE(glf_etg_output(grid_size))
          IF(ALLOCATED(glf_anrate_output))DEALLOCATE(glf_anrate_output)
              ALLOCATE(glf_anrate_output(grid_size))
          IF(ALLOCATED(glf_anrate2_output))DEALLOCATE(glf_anrate2_output)
              ALLOCATE(glf_anrate2_output(grid_size))
          IF(ALLOCATED(glf_anfreq_output))DEALLOCATE(glf_anfreq_output)
              ALLOCATE(glf_anfreq_output(grid_size))
          IF(ALLOCATED(glf_anfreq2_output))DEALLOCATE(glf_anfreq2_output)
              ALLOCATE(glf_anfreq2_output(grid_size))
          IF(ALLOCATED(glf_gamma_net_i_output))DEALLOCATE(glf_gamma_net_i_output)
              ALLOCATE(glf_gamma_net_i_output(grid_size))
          IF(ALLOCATED(glf_gamma_net_e_output))DEALLOCATE(glf_gamma_net_e_output)
              ALLOCATE(glf_gamma_net_e_output(grid_size))
          IF(ALLOCATED(glf_eff_chi_output))DEALLOCATE(glf_eff_chi_output)
              ALLOCATE(glf_eff_chi_output(grid_size,ntot))


               glf_p_output(:,:)         = zeroc             ! these arrays are sent to the state file
               glf_e_output(:,:)         = zeroc             ! as zeros if glf is not active
               glf_m_output(:,:)         = zeroc  
               glf_etg_output(:)         = zeroc 
               glf_anrate_output(:)      = zeroc
               glf_anrate2_output(:)     = zeroc
               glf_anfreq_output(:)      = zeroc
               glf_anfreq2_output(:)     = zeroc
               glf_gamma_net_e_output(:) = zeroc
               glf_gamma_net_i_output(:) = zeroc
               glf_eff_chi_output(:,:)     = zeroc

         RETURN
      END SUBROUTINE allocate_glf_output



      SUBROUTINE glf_allocate
!-------------------------------------------------------------------------------------------------
! --
!-------------------------------------------------------------------------------------------------
          USE solcon_gcnmp,                                  ONLY : itran_max,n_slave_grid_tot, &
                                                                    t_slave_compt_tot
          USE MPI_data,                                      ONLY :  myid,master,numprocs


          IF( .NOT. ALLOCATED(itport_pt))ALLOCATE(itport_pt(itran_max))
          itport_pt(:) = izero
          IF( .NOT. ALLOCATED(rho_glf))ALLOCATE(rho_glf(nj))
          rho_glf(:) = zeroc
          IF( .NOT. ALLOCATED(tau_glf))ALLOCATE(tau_glf(nj))
          tau_glf(:) = zeroc
          IF( .NOT. ALLOCATED(ne_glf))ALLOCATE(ne_glf(nj))
          ne_glf(:) = zeroc
          IF( .NOT. ALLOCATED(ni_glf))ALLOCATE(ni_glf(nj))
          ni_glf(:) = zeroc
          IF( .NOT. ALLOCATED(ns_glf))ALLOCATE(ns_glf(nj))
          ns_glf(:) = zeroc
          IF( .NOT. ALLOCATED(zeff_exp))ALLOCATE(zeff_exp(nj))
          zeff_exp(:) = zeroc
          IF( .NOT. ALLOCATED(betae_glf))ALLOCATE(betae_glf(nj))
          betae_glf(:) = zeroc
          IF( .NOT. ALLOCATED(betai_glf))ALLOCATE(betai_glf(nj))
          betai_glf(:) = zeroc
          IF( .NOT. ALLOCATED(beta_glf))ALLOCATE(beta_glf(nj))
          beta_glf(:) = zeroc
          IF( .NOT. ALLOCATED(q_glf))ALLOCATE(q_glf(nj))
          q_glf(:) = zeroc
          IF( .NOT. ALLOCATED(rmin_glf))ALLOCATE(rmin_glf(nj))
          rmin_glf(:) = zeroc
          IF( .NOT. ALLOCATED(rmaj_glf))ALLOCATE(rmaj_glf(nj))
          rmaj_glf(:) = zeroc

          IF( .NOT. ALLOCATED(grho1_glf))ALLOCATE(grho1_glf(nj))
            grho1_glf(:) = dischg%grho1nj%data(:)      ! <abs(grad rho_glf) > 1/m
          IF( .NOT. ALLOCATED(grho2_glf))ALLOCATE(grho2_glf(nj))
            grho2_glf(:) = dischg%grho2nj%data(:)      ! <(grad rho_glf)**2 > 1/m**2

          IF( .NOT. ALLOCATED(elong_glf))ALLOCATE(elong_glf(nj))
          elong_glf(:) = zeroc
          IF( .NOT. ALLOCATED(drhodr))ALLOCATE(drhodr(nj))
          drhodr(:) = zeroc
          IF( .NOT. ALLOCATED(drhodrrrho))ALLOCATE(drhodrrrho(nj))
          drhodrrrho(:) = zeroc
          IF( .NOT. ALLOCATED(csda_glf))ALLOCATE(csda_glf(nj))
          csda_glf(:) = zeroc 
          IF( .NOT. ALLOCATED(rhosda_glf))ALLOCATE(rhosda_glf(nj))
          rhosda_glf(:) = zeroc 
          IF( .NOT. ALLOCATED(bteff_glf))ALLOCATE(bteff_glf(nj))
          bteff_glf(:) = zeroc 
          IF( .NOT. ALLOCATED(drhogf))ALLOCATE(drhogf(nj))
          drhogf(:) = zeroc 
          IF( .NOT. ALLOCATED(zpte_glf))ALLOCATE(zpte_glf(nj))
          zpte_glf(:) = zeroc 
          IF( .NOT. ALLOCATED(zpti_glf))ALLOCATE(zpti_glf(nj))
          zpti_glf(:) = zeroc 
          IF( .NOT. ALLOCATED(zpne_glf))ALLOCATE(zpne_glf(nj))
          zpne_glf(:) = zeroc 
          IF( .NOT. ALLOCATED(zpni_glf))ALLOCATE(zpni_glf(nj))
          zpni_glf(:) = zeroc 
          IF( .NOT. ALLOCATED(shat_exp))ALLOCATE(shat_exp(nj))
          shat_exp(:) = zeroc 
          IF( .NOT. ALLOCATED(alpha_exp))ALLOCATE(alpha_exp(nj))
          alpha_exp(:) = zeroc 
          IF( .NOT. ALLOCATED(alpha_exp_save))ALLOCATE(alpha_exp_save(nj))
          alpha_exp_save(:) = zeroc 
          IF( .NOT. ALLOCATED(fcglf))ALLOCATE(fcglf(nj))
          fcglf(:) = zeroc 
          IF( .NOT. ALLOCATED(ak1glf))ALLOCATE(ak1glf(nj))
          ak1glf(:) = zeroc 
          IF( .NOT. ALLOCATED(ak2glf))ALLOCATE(ak2glf(nj))
          ak2glf(:) = zeroc 
          IF( .NOT. ALLOCATED(aneoglf))ALLOCATE(aneoglf(nj))
          aneoglf(:) = zeroc
          IF( .NOT. ALLOCATED(angrotp_exp))ALLOCATE(angrotp_exp(nj))
          angrotp_exp(:) = zeroc
          IF( .NOT. ALLOCATED(ve_glf))ALLOCATE(ve_glf(nj))
          ve_glf(:) = zeroc
          IF( .NOT. ALLOCATED(vparm_glf))ALLOCATE(vparm_glf(nj))
          vparm_glf(:) = zeroc
          IF( .NOT. ALLOCATED(vparm_term1))ALLOCATE(vparm_term1(nj))
          vparm_term1(:) = zeroc
          IF( .NOT. ALLOCATED(vperm_glf))ALLOCATE(vperm_glf(nj))
          vperm_glf(:) = zeroc
          IF( .NOT. ALLOCATED(vpolm_glf))ALLOCATE(vpolm_glf(nj))
          vpolm_glf(:) = zeroc
          IF( .NOT. ALLOCATED(vperm_term1))ALLOCATE(vperm_term1(nj))
          vperm_term1(:) = zeroc
          IF( .NOT. ALLOCATED(vphim_glf))ALLOCATE(vphim_glf(nj))
          vphim_glf(:) = zeroc
          IF( .NOT. ALLOCATED(egamma_exp))ALLOCATE(egamma_exp(nj))
          egamma_exp(:) = zeroc
          IF( .NOT. ALLOCATED(egamma_exp_zc))ALLOCATE(egamma_exp_zc(0:nj-2))
          egamma_exp_zc(:) = zeroc
          IF( .NOT. ALLOCATED(gamma_p_exp))ALLOCATE(gamma_p_exp(nj))
          gamma_p_exp(:) = zeroc
          IF( .NOT. ALLOCATED(gamma_p_exp_zc))ALLOCATE(gamma_p_exp_zc(0:nj-2))
          gamma_p_exp_zc(:) = zeroc
          IF( .NOT. ALLOCATED(diff_m))ALLOCATE(diff_m(nj))
          diff_m(:) = zeroc
          IF( .NOT. ALLOCATED(diff_mdv))ALLOCATE(diff_mdv(nj))
          diff_mdv(:) = zeroc
          IF( .NOT. ALLOCATED(chie_m))ALLOCATE(chie_m(nj))
          chie_m(:) = zeroc
          IF( .NOT. ALLOCATED(chie_mdv))ALLOCATE(chie_mdv(nj))
          chie_mdv(:) = zeroc
          IF( .NOT. ALLOCATED(chii_m))ALLOCATE(chii_m(nj))
          chii_m(:) = zeroc
          IF( .NOT. ALLOCATED(chii_mdv))ALLOCATE(chii_mdv(nj))
          chii_mdv(:) = zeroc
          IF( .NOT. ALLOCATED(etaphi_m))ALLOCATE(etaphi_m(nj))
          etaphi_m(:) = zeroc
          IF( .NOT. ALLOCATED(etaphi_mdv))ALLOCATE(etaphi_mdv(nj))
          etaphi_mdv(:) = zeroc
          IF( .NOT. ALLOCATED(etapar_m))ALLOCATE(etapar_m(nj))
          etapar_m(:) = zeroc
          IF( .NOT. ALLOCATED(etapar_mdv))ALLOCATE(etapar_mdv(nj))
          etapar_mdv(:) = zeroc
          IF( .NOT. ALLOCATED(etaper_m))ALLOCATE(etaper_m(nj))
          etaper_m(:) = zeroc
          IF( .NOT. ALLOCATED(etaper_mdv))ALLOCATE(etaper_mdv(nj))
          etaper_mdv(:) = zeroc
          IF( .NOT. ALLOCATED(exch_m))ALLOCATE(exch_m(nj))
          exch_m(:) = zeroc
          IF( .NOT. ALLOCATED(egamma_m))ALLOCATE(egamma_m(nj))
          egamma_m(:) = zeroc
          IF( .NOT. ALLOCATED(egamma_mdv))ALLOCATE(egamma_mdv(nj))
          egamma_mdv(:) = zeroc
          IF( .NOT. ALLOCATED(egamma_d))ALLOCATE(egamma_d(nj,pdelay))
          egamma_d(:,:) = zeroc
          IF( .NOT. ALLOCATED(gamma_p_m))ALLOCATE(gamma_p_m(nj))
          gamma_p_m(:) = zeroc
          IF( .NOT. ALLOCATED(gamma_p_mdv))ALLOCATE(gamma_p_mdv(nj))
          gamma_p_mdv(:) = zeroc
          IF( .NOT. ALLOCATED(anrate_m))THEN 
             ALLOCATE(anrate_m(nj))
             anrate_m(:) = zeroc
          ENDIF
          IF( .NOT. ALLOCATED(anrate2_m))THEN
            ALLOCATE(anrate2_m(nj))
            anrate2_m(:)  = zeroc
          ENDIF

          IF( .NOT. ALLOCATED(anrate_mdv))ALLOCATE(anrate_mdv(nj))
           anrate_mdv(:)  = zeroc
          IF( .NOT. ALLOCATED(anrate2_mdv))ALLOCATE(anrate2_mdv(nj))
          anrate2_mdv(:) = zeroc

          IF( .NOT. ALLOCATED(anfreq_m))THEN
             ALLOCATE(anfreq_m(nj))
             anfreq_m(:) = zeroc
          ENDIF
          IF( .NOT. ALLOCATED(anfreq2_m))THEN
             ALLOCATE(anfreq2_m(nj))
             anfreq2_m(:) = zeroc
          ENDIF


          IF( .NOT. ALLOCATED(anfreq_mdv))ALLOCATE(anfreq_mdv(nj))
           anfreq_mdv(:)  = zeroc
          IF( .NOT. ALLOCATED(anfreq2_mdv))ALLOCATE(anfreq2_mdv(nj))
           anfreq2_mdv(:) = zeroc
          IF( .NOT. ALLOCATED(vp_eff_m))ALLOCATE(vp_eff_m(nj))
          vp_eff_m(:) = zeroc
          IF( .NOT. ALLOCATED(ve_eff_m))ALLOCATE(ve_eff_m(nj))
          ve_eff_m(:) = zeroc
          IF( .NOT. ALLOCATED(vi_eff_m))ALLOCATE(vi_eff_m(nj))
           vi_eff_m(:)    = zeroc 
          IF( .NOT. ALLOCATED(vrot_eff_m))ALLOCATE(vrot_eff_m(nj))
           vrot_eff_m(:)  = zeroc
          IF( .NOT. ALLOCATED(chi_e_avg))ALLOCATE(chi_e_avg(nj))
          chi_e_avg(:) = zeroc
          IF( .NOT. ALLOCATED(chi_i_avg))ALLOCATE(chi_i_avg(nj))
          chi_i_avg(:) = zeroc
          IF( .NOT. ALLOCATED(exchmdv))ALLOCATE(exchmdv(nj))
           exchmdv(:)     = zeroc
          IF( .NOT. ALLOCATED(exch_mdv))ALLOCATE(exch_mdv(nj))
           exch_mdv(:)    = zeroc
          IF( .NOT. ALLOCATED(chie_ms))ALLOCATE(chie_ms(nj))
          chie_ms(:) = zeroc
          IF( .NOT. ALLOCATED(chii_ms))ALLOCATE(chii_ms(nj))
          chii_ms(:) = zeroc
          IF( .NOT. ALLOCATED(etaphi_ms))ALLOCATE(etaphi_ms(nj))
          etaphi_ms(:) = zeroc
          IF( .NOT. ALLOCATED(diff_ms))ALLOCATE(diff_ms(nj))
          diff_ms(:) = zeroc
          IF( .NOT. ALLOCATED(diff_avg))ALLOCATE(diff_avg(nj))
           diff_avg(:)       = zeroc
          IF( .NOT. ALLOCATED(d_avg))ALLOCATE(d_avg(nj))
           d_avg(:)       = zeroc
          IF( .NOT. ALLOCATED(xkang_avg))ALLOCATE(xkang_avg(nj))
          xkang_avg(:)       = zeroc
          IF( .NOT. ALLOCATED(xchii_sv))ALLOCATE(xchii_sv(nj))
          xchii_sv(:)       = zeroc
          IF( .NOT. ALLOCATED(d_sv))ALLOCATE(d_sv(nj))
          d_sv(:)       = zeroc
          IF( .NOT. ALLOCATED(xchie_sv))ALLOCATE(xchie_sv(nj))
          xchie_sv(:) = zeroc
          IF( .NOT. ALLOCATED(xkangwlsv))ALLOCATE(xkangwlsv(nj))
          xkangwlsv(:) = zeroc
          IF( .NOT. ALLOCATED(xke_glf23))ALLOCATE(xke_glf23(nj))
          xke_glf23(:) = zeroc
          IF( .NOT. ALLOCATED(xki_glf23))ALLOCATE(xki_glf23(nj))
          xki_glf23(:) = zeroc
          IF( .NOT. ALLOCATED(xkw_glf23))ALLOCATE(xkw_glf23(nj))
          xkw_glf23(:) = zeroc
          IF( .NOT. ALLOCATED(te_glf))ALLOCATE(te_glf(nj))
          te_glf(:) = zeroc
          IF( .NOT. ALLOCATED(ti_glf))ALLOCATE(ti_glf(nj))
          ti_glf=zeroc

          IF( .NOT. ALLOCATED(dqrm_drho))ALLOCATE(dqrm_drho(0:nj-2))
          IF( .NOT. ALLOCATED(vexb_base))ALLOCATE(vexb_base(nj))
          IF( .NOT. ALLOCATED(vexb_term1))ALLOCATE(vexb_term1(nj))
          IF( .NOT. ALLOCATED(diamag_term))ALLOCATE(diamag_term(nj))
          IF( .NOT. ALLOCATED(rdrho))ALLOCATE(rdrho(nj))
          IF( .NOT. ALLOCATED(dvexb_base_drho))ALLOCATE(dvexb_base_drho(0:nj-2))
          IF( .NOT. ALLOCATED(ddiamag_drho))ALLOCATE(ddiamag_drho(0:nj-2))
          IF( .NOT. ALLOCATED(dvexb_term1_drho))ALLOCATE(dvexb_term1_drho(0:nj-2))
          IF( .NOT. ALLOCATED(dvpar_term1_drho))ALLOCATE(dvpar_term1_drho(0:nj-2))
          IF( .NOT. ALLOCATED(dvphim_drho))ALLOCATE(dvphim_drho(0:nj-2))
          IF(.NOT. ALLOCATED(dcoefp_glf23))ALLOCATE(dcoefp_glf23(ntot,nj))
          dcoefp_glf23(:,:)= zeroc
          IF( .NOT. ALLOCATED(gamma_net_i_loc))ALLOCATE(gamma_net_i_loc(nj))
          IF( .NOT. ALLOCATED(gamma_net_e_loc))ALLOCATE(gamma_net_e_loc(nj))
             gamma_net_i_loc(:) = zeroc ; gamma_net_e_loc(:) =zeroc
       
            
             IF( .NOT. ALLOCATED(niflux_gf_m))THEN
                ALLOCATE(niflux_gf_m(nj))
                niflux_gf_m(:) = zeroc
             ENDIF
             IF( .NOT. ALLOCATED(nimpflux_gf_m))THEN
                ALLOCATE(nimpflux_gf_m(nj))
                nimpflux_gf_m = zeroc
             ENDIF

             IF( .NOT. ALLOCATED(etgflux_gf_m))THEN
                ALLOCATE(etgflux_gf_m(nj))
                etgflux_gf_m(:) = zeroc
             ENDIF

             IF( .NOT. ALLOCATED(tiflux_gf_m)) THEN
                      ALLOCATE(tiflux_gf_m(nj))
                      tiflux_gf_m(:) = zeroc
             ENDIF
             IF( .NOT. ALLOCATED(teflux_gf_m)) THEN
                      ALLOCATE(teflux_gf_m(nj))
                      teflux_gf_m(:) = zeroc
             ENDIF
             IF( .NOT. ALLOCATED(vparflux_gf_m))THEN
                ALLOCATE(vparflux_gf_m(nj))
                vparflux_gf_m(:) = zeroc
             ENDIF
             IF( .NOT. ALLOCATED(vperflux_gf_m))THEN
                 ALLOCATE(vperflux_gf_m(nj))
                 vperflux_gf_m(:) = zeroc
             ENDIF
             IF( .NOT. ALLOCATED(vphiflux_gf_m))THEN
                ALLOCATE(vphiflux_gf_m(nj))
                vphiflux_gf_m(:) = zeroc
             ENDIF
             IF( .NOT. ALLOCATED(gamma_net_i))THEN
                ALLOCATE(gamma_net_i(nj))
                    gamma_net_i(:) =zeroc
             ENDIF
             IF( .NOT. ALLOCATED(gamma_net_e))THEN 
                    ALLOCATE(gamma_net_e(nj))
                    gamma_net_e(:) =zeroc
             ENDIF

            !Arrays that will be passed around to processors
            glf_max_pack      = nspecies_glf*11*3+2*11 ! cummulative size of following arrays:
            IF( .NOT. ALLOCATED(glf_p_flux))   ALLOCATE(glf_p_flux(0:nj-1,nspecies_glf,-itran_max:itran_max))
            IF( .NOT. ALLOCATED(glf_e_flux))   ALLOCATE(glf_e_flux(0:nj-1,nspecies_glf,-itran_max:itran_max))
            IF( .NOT. ALLOCATED(glf_m_flux))   ALLOCATE(glf_m_flux(0:nj-1,nspecies_glf,-itran_max:itran_max))
            IF( .NOT. ALLOCATED(glf_val_pert)) ALLOCATE(glf_val_pert(-itran_max:itran_max,0:nj-1))
            IF( .NOT. ALLOCATED(glf_val_base)) ALLOCATE(glf_val_base(-itran_max:itran_max,0:nj-1))
             glf_val_base(:,:) = zeroc
             glf_val_pert(:,:) = zeroc
            glf_max_pack       = glf_max_pack + 7
            ! 7 => jm + anrate +anrate2 + anfreq +anfreq2 + slave time + slave grid points
            ! anrate,anrate2,anfreq,anfreq2 count as single element scalar quantites 
            ! here becasue only one lement of each is packed and sent out at a time.
            IF( ALLOCATED(glf_receive_packed))DEALLOCATE(glf_receive_packed)
            IF( ALLOCATED(glf_send_packed)) DEALLOCATE(glf_send_packed)
            IF( ALLOCATED(n_slave_grid_tot))DEALLOCATE(n_slave_grid_tot)
            IF( ALLOCATED(t_slave_compt_tot))DEALLOCATE(t_slave_compt_tot)
            glf_max_pack       = glf_max_pack +1 ! niflux_gf_m(j) 
            glf_max_pack       = glf_max_pack +1 ! nimpflux_gf_m(j) 
            glf_max_pack       = glf_max_pack +1 ! tiflux_gf_m(j) 
            glf_max_pack       = glf_max_pack +1 ! teflux_gf_m(j) 
            glf_max_pack       = glf_max_pack +1 ! vparflux_gf_m(j) 
            glf_max_pack       = glf_max_pack +1 ! vperflux_gf_m(j) 
            glf_max_pack       = glf_max_pack +1 ! vphiflux_gf_m(j) 
            glf_max_pack       = glf_max_pack +1 ! etg_flux_gf_m(j) 
            glf_max_pack       = glf_max_pack +1 ! gamma_net_i(j) 
            glf_max_pack       = glf_max_pack +1 ! gamma_net_e(j) 
          IF( myid .NE. master)THEN
             ALLOCATE(glf_send_packed(glf_max_pack))
          ELSE
             ALLOCATE(glf_receive_packed(glf_max_pack))
             ALLOCATE(n_slave_grid_tot(0:numprocs-1),t_slave_compt_tot(0:numprocs-1))
             t_slave_compt_tot(:)    = zeroc 
             n_slave_grid_tot(:)     = izero
          ENDIF

             etaphi_mdv(:)  = zeroc


           vrot_eff_m(:)  = zeroc
           chie_m(:)      = zeroc
           chii_m(:)      = zeroc
           etaphi_m(:)    = zeroc
           xkang_avg(:)   = zeroc

!           vphiflux_gf_m(:) = zeroc
!           vperflux_gf_m(:) = zeroc
!           vparflux_gf_m(:) = zeroc
!           tiflux_gf_m(:)   = zeroc
!           teflux_gf_m(:)   = zeroc
!           etgflux_gf_m(:)  = zeroc
           
          IF( .NOT. ALLOCATED(chie_avg))THEN 
            ALLOCATE(chie_avg(nj))
            chie_avg(:) =zeroc
          ENDIF
          IF(.NOT. ALLOCATED(chii_avg))THEN 
           ALLOCATE(chii_avg(nj))
           chii_avg(:) =zeroc
          ENDIF
          IF(.NOT. ALLOCATED(etaphi_avg))THEN 
           ALLOCATE(etaphi_avg(nj))
           etaphi_avg(:) =zeroc
          ENDIF
       END SUBROUTINE glf_allocate


      SUBROUTINE glf_allocate_lcl
!----------------------------------------------------------------------------------------------------
! -- Allocate arrays used only locally in glf
!----------------------------------------------------------------------------------------------------

     USE glf23_gcnmp,                        ONLY :                                                 &
                                                   glf_te_fg,glf_ti_fg,glf_ne_fg,glf_ni_fg,         &
                                                   glf_ns_fg,glf_angrot_fg,glf_exch_m,              &
                                                   glf_egamma_exp,glf_gamma_p_exp,glf_vphim,        &
                                                   glf_vparm,glf_vperm,glf_zeff,glf_rho,glf_grho1,  &
                                                   glf_grho2,glf_rmin,glf_rmaj,glf_q,glf_shat,      &
                                                   glf_alpha,glf_elong,glf_diff_m,glf_chie_m,       &
                                                   glf_chii_m,glf_etaphi_m,glf_etapar_m,            &
                                                   glf_etaper_m,glf_egamma_m,glf_egamma_d,          &
                                                   glf_gamma_p_m,glf_anrate,glf_anrate2,glf_anfreq, &
                                                   glf_anfreq2,glf_ni_flux,glf_nimp_flux,           &
                                                   glf_ti_flux,glf_te_flux,glf_vparflux_gf_m,       &
                                                   glf_vperflux_gf_m,glf_vphiflux_gf_m,             &
                                                   glf_etgflux_gf_m,glf_gamma_net_i,glf_gamma_net_e,&
                                                   glf_egamma_zc,glf_gamma_p_zc
                                               

     USE grid_class,                         ONLY : nj

          IF(.NOT. ALLOCATED(glf_te_fg))           ALLOCATE(glf_te_fg(0:nj-1))
          IF(.NOT. ALLOCATED(glf_ti_fg))           ALLOCATE(glf_ti_fg(0:nj-1))
          IF(.NOT. ALLOCATED(glf_ne_fg))           ALLOCATE(glf_ne_fg(0:nj-1))
          IF(.NOT. ALLOCATED(glf_ni_fg))           ALLOCATE(glf_ni_fg(0:nj-1))
          IF(.NOT. ALLOCATED(glf_ns_fg))           ALLOCATE(glf_ns_fg(0:nj-1))
          IF(.NOT. ALLOCATED(glf_angrot_fg))       ALLOCATE(glf_angrot_fg(0:nj-1))
          IF(.NOT. ALLOCATED(glf_egamma_exp))      ALLOCATE(glf_egamma_exp(0:nj-1))
          IF(.NOT. ALLOCATED(glf_gamma_p_exp))     ALLOCATE(glf_gamma_p_exp(0:nj-1))
          IF(.NOT. ALLOCATED(glf_vphim))           ALLOCATE(glf_vphim(0:nj-1))
          IF(.NOT. ALLOCATED(glf_vparm))           ALLOCATE(glf_vparm(0:nj-1))
          IF(.NOT. ALLOCATED(glf_vperm))           ALLOCATE(glf_vperm(0:nj-1))
          IF(.NOT. ALLOCATED(glf_zeff))            ALLOCATE(glf_zeff(0:nj-1))
          IF(.NOT. ALLOCATED(glf_rho))             ALLOCATE(glf_rho(0:nj-1))
          IF(.NOT. ALLOCATED(glf_grho1))           ALLOCATE(glf_grho1(0:nj-1))
          IF(.NOT. ALLOCATED(glf_grho2))           ALLOCATE(glf_grho2(0:nj-1))
          IF(.NOT. ALLOCATED(glf_rmin))            ALLOCATE(glf_rmin(0:nj-1))
          IF(.NOT. ALLOCATED(glf_rmaj))            ALLOCATE(glf_rmaj(0:nj-1))
          IF(.NOT. ALLOCATED(glf_q))               ALLOCATE(glf_q(0:nj-1))
          IF(.NOT. ALLOCATED(glf_shat))            ALLOCATE(glf_shat(0:nj-1))
          IF(.NOT. ALLOCATED(glf_alpha))           ALLOCATE(glf_alpha(0:nj-1))
          IF(.NOT. ALLOCATED(glf_elong))           ALLOCATE(glf_elong(0:nj-1))

          IF(.NOT. ALLOCATED(glf_diff_m))          ALLOCATE(glf_diff_m(0:nj-1))
          IF(.NOT. ALLOCATED(glf_chie_m))          ALLOCATE(glf_chie_m(0:nj-1))
          IF(.NOT. ALLOCATED(glf_chii_m))          ALLOCATE(glf_chii_m(0:nj-1))
          IF(.NOT. ALLOCATED(glf_etaphi_m))        ALLOCATE(glf_etaphi_m(0:nj-1))
          IF(.NOT. ALLOCATED(glf_etapar_m))        ALLOCATE(glf_etapar_m(0:nj-1))
          IF(.NOT. ALLOCATED(glf_etaper_m))        ALLOCATE(glf_etaper_m(0:nj-1))

          IF(.NOT. ALLOCATED(glf_egamma_m))        ALLOCATE(glf_egamma_m(0:nj-1))
          IF(.NOT. ALLOCATED(glf_egamma_d))        ALLOCATE(glf_egamma_d(0:nj-1,1:10))
          IF(.NOT. ALLOCATED(glf_exch_m))          ALLOCATE(glf_exch_m(0:nj-1))

          IF(.NOT. ALLOCATED(glf_gamma_p_m))       ALLOCATE(glf_gamma_p_m(0:nj-1))
          IF(.NOT. ALLOCATED(glf_anrate))          ALLOCATE(glf_anrate(0:nj-1))
          IF(.NOT. ALLOCATED(glf_anrate2))         ALLOCATE(glf_anrate2(0:nj-1))
          IF(.NOT. ALLOCATED(glf_anfreq))          ALLOCATE(glf_anfreq(0:nj-1))
          IF(.NOT. ALLOCATED(glf_anfreq2))         ALLOCATE(glf_anfreq2(0:nj-1))
          IF(.NOT. ALLOCATED(glf_ni_flux))         ALLOCATE(glf_ni_flux(0:nj-1))
          IF(.NOT. ALLOCATED(glf_nimp_flux))       ALLOCATE(glf_nimp_flux(0:nj-1))
          IF(.NOT. ALLOCATED(glf_ti_flux))         ALLOCATE(glf_ti_flux(0:nj-1))
          IF(.NOT. ALLOCATED(glf_te_flux))         ALLOCATE(glf_te_flux(0:nj-1))
          IF(.NOT. ALLOCATED(glf_vparflux_gf_m))   ALLOCATE(glf_vparflux_gf_m(0:nj-1))
          IF(.NOT. ALLOCATED(glf_vperflux_gf_m ))  ALLOCATE(glf_vperflux_gf_m(0:nj-1))
          IF(.NOT. ALLOCATED(glf_vphiflux_gf_m))   ALLOCATE(glf_vphiflux_gf_m(0:nj-1))
          IF(.NOT. ALLOCATED(glf_etgflux_gf_m))    ALLOCATE(glf_etgflux_gf_m(0:nj-1))
          IF(.NOT. ALLOCATED(glf_gamma_net_i))     ALLOCATE(glf_gamma_net_i(0:nj-1))
          IF(.NOT. ALLOCATED(glf_gamma_net_e))     ALLOCATE(glf_gamma_net_e(0:nj-1))
          IF(.NOT. ALLOCATED(glf_egamma_zc))       ALLOCATE(glf_egamma_zc(0:nj-2))
          IF(.NOT. ALLOCATED(glf_gamma_p_zc))      ALLOCATE(glf_gamma_p_zc(0:nj-2))
          glf_vperflux_gf_m(:) = zeroc         ;   glf_vphiflux_gf_m(:) = zeroc
          glf_vparflux_gf_m(:) = zeroc         ;   glf_etgflux_gf_m(:)   = zeroc
          glf_gamma_net_i(:)   = zeroc         ;   glf_gamma_net_e(:)    = zeroc
          glf_te_flux(:)       = zeroc         ;   glf_ti_flux(:)        = zeroc
          glf_nimp_flux(:)     = zeroc         ;   glf_ni_flux(:)        = zeroc

 
          RETURN

       END SUBROUTINE glf_allocate_lcl


      SUBROUTINE glf_deallocate_lcl
!----------------------------------------------------------------------------------------------------
! -- Deallocate arrays used only locally in glf
!----------------------------------------------------------------------------------------------------

     USE glf23_gcnmp,                        ONLY :                                                 &
                                                   glf_te_fg,glf_ti_fg,glf_ne_fg,glf_ni_fg,         &
                                                   glf_ns_fg,glf_angrot_fg,glf_exch_m,              &
                                                   glf_egamma_exp,glf_gamma_p_exp,glf_vphim,        &
                                                   glf_vparm,glf_vperm,glf_zeff,glf_rho,glf_grho1,  &
                                                   glf_grho2,glf_rmin,glf_rmaj,glf_q,glf_shat,      &
                                                   glf_alpha,glf_elong,glf_diff_m,glf_chie_m,       &
                                                   glf_chii_m,glf_etaphi_m,glf_etapar_m,            &
                                                   glf_etaper_m,glf_egamma_m,glf_egamma_d,          &
                                                   glf_gamma_p_m,glf_anrate,glf_anrate2,glf_anfreq, &
                                                   glf_anfreq2,glf_ni_flux,glf_nimp_flux,           &
                                                   glf_ti_flux,glf_te_flux,glf_vparflux_gf_m,       &
                                                   glf_vperflux_gf_m,glf_vphiflux_gf_m,             &
                                                   glf_etgflux_gf_m,glf_gamma_net_i,                &
                                                   glf_gamma_net_e,glf_egamma_zc,glf_gamma_p_zc



          IF(ALLOCATED(glf_te_fg))           DEALLOCATE(glf_te_fg )
          IF(ALLOCATED(glf_ti_fg))           DEALLOCATE(glf_ti_fg )
          IF(ALLOCATED(glf_ne_fg))           DEALLOCATE(glf_ne_fg )
          IF(ALLOCATED(glf_ni_fg))           DEALLOCATE(glf_ni_fg )
          IF(ALLOCATED(glf_ns_fg))           DEALLOCATE(glf_ns_fg )
          IF(ALLOCATED(glf_angrot_fg))       DEALLOCATE(glf_angrot_fg )
          IF(ALLOCATED(glf_egamma_exp))      DEALLOCATE(glf_egamma_exp )
          IF(ALLOCATED(glf_gamma_p_exp))     DEALLOCATE(glf_gamma_p_exp )
          IF(ALLOCATED(glf_vphim))           DEALLOCATE(glf_vphim )
          IF(ALLOCATED(glf_vparm))           DEALLOCATE(glf_vparm )
          IF(ALLOCATED(glf_vperm))           DEALLOCATE(glf_vperm )
          IF(ALLOCATED(glf_zeff))            DEALLOCATE(glf_zeff )
          IF(ALLOCATED(glf_rho))             DEALLOCATE(glf_rho )
          IF(ALLOCATED(glf_grho1))           DEALLOCATE(glf_grho1 )
          IF(ALLOCATED(glf_grho2))           DEALLOCATE(glf_grho2 )
          IF(ALLOCATED(glf_rmin))            DEALLOCATE(glf_rmin )
          IF(ALLOCATED(glf_rmaj))            DEALLOCATE(glf_rmaj )
          IF(ALLOCATED(glf_q))               DEALLOCATE(glf_q )
          IF(ALLOCATED(glf_shat))            DEALLOCATE(glf_shat )
          IF(ALLOCATED(glf_alpha))           DEALLOCATE(glf_alpha )
          IF(ALLOCATED(glf_elong))           DEALLOCATE(glf_elong )

          IF(ALLOCATED(glf_diff_m))          DEALLOCATE(glf_diff_m )
          IF(ALLOCATED(glf_chie_m))          DEALLOCATE(glf_chie_m )
          IF(ALLOCATED(glf_chii_m))          DEALLOCATE(glf_chii_m )
          IF(ALLOCATED(glf_etaphi_m))        DEALLOCATE(glf_etaphi_m )
          IF(ALLOCATED(glf_etapar_m))        DEALLOCATE(glf_etapar_m )
          IF(ALLOCATED(glf_etaper_m))        DEALLOCATE(glf_etaper_m )

          IF(ALLOCATED(glf_egamma_m))        DEALLOCATE(glf_egamma_m )
          IF(ALLOCATED(glf_egamma_d))        DEALLOCATE(glf_egamma_d )
          IF(ALLOCATED(glf_exch_m))          DEALLOCATE(glf_exch_m)
          IF(ALLOCATED(glf_gamma_p_m))       DEALLOCATE(glf_gamma_p_m )
         !   IF(ALLOCATED(glf_anrate))          DEALLOCATE(glf_anrate )
         !   IF(ALLOCATED(glf_anrate2))         DEALLOCATE(glf_anrate2 )  
         !   IF(ALLOCATED(glf_anfreq))          DEALLOCATE(glf_anfreq )
         !   IF(ALLOCATED(glf_anfreq2))         DEALLOCATE(glf_anfreq2 )
          IF(ALLOCATED(glf_ni_flux))         DEALLOCATE(glf_ni_flux )
          IF(ALLOCATED(glf_nimp_flux))       DEALLOCATE(glf_nimp_flux )
          IF(ALLOCATED(glf_ti_flux))         DEALLOCATE(glf_ti_flux )
          IF(ALLOCATED(glf_te_flux))         DEALLOCATE(glf_te_flux )
          IF(ALLOCATED(glf_vparflux_gf_m))   DEALLOCATE(glf_vparflux_gf_m )
          IF(ALLOCATED(glf_vperflux_gf_m ))  DEALLOCATE(glf_vperflux_gf_m )
          IF(ALLOCATED(glf_vphiflux_gf_m))   DEALLOCATE(glf_vphiflux_gf_m )
          IF(ALLOCATED(glf_etgflux_gf_m))    DEALLOCATE(glf_etgflux_gf_m )
          IF(ALLOCATED(glf_gamma_net_i))     DEALLOCATE(glf_gamma_net_i )
          IF(ALLOCATED(glf_gamma_net_e))     DEALLOCATE(glf_gamma_net_e )
          IF(ALLOCATED(glf_egamma_zc))       DEALLOCATE(glf_egamma_zc)
          IF(ALLOCATED(glf_gamma_p_zc))      DEALLOCATE(glf_gamma_p_zc)  


          RETURN

       END SUBROUTINE glf_deallocate_lcl



       SUBROUTINE glf_input_defaults
!------------------------------------------------------------------------------
! -- 
!------------------------------------------------------------------------------
          USE nrtype,                 ONLY : I4B,DP
          USE glf23_gcnmp,                  ONLY : exbmult_glf,x_alpha_glf,                    &
                                             jroot_glf,limp_glf,ibtflag_glf,irotstab,    &
                                              dv_method,spline_gsc,dv_delt,              &
                                             itte_dv,itti_dv,itangrot_dv,itenp_dv,       &
                                             itene_dv,jeigen_glf,i_delay,use_mask,       &
                                             t_delay_glf23
          USE common_constants,       ONLY : izero,zeroc
          IMPLICIT NONE
          
             exbmult_glf = 1.35_DP    ! alphae in glf23 ( E x B shear stabilization)
             x_alpha_glf = -1.0_DP    ! alpha_stabilization
             jroot_glf   = 8          ! 8 equation model (no impurity dynamics )  
             jeigen_glf  = izero      ! selects eigenvalue solver
             limp_glf    = izero      ! no impurity dynamics
             ibtflag_glf = 1          ! switch for effective B field 
             irotstab    = 1          ! rotational stabilization
             dv_method   = izero      ! dv_method is not implemented
             spline_gsc  = izero      ! dont spline profiles to get derivatives
             dv_delt     = 0.01
             itte_dv     = izero 
             itti_dv     = izero  
             itangrot_dv = izero 
             itenp_dv    = izero 
             itene_dv    = izero 
             i_delay     = izero 
             use_mask    = izero
             t_delay_glf23 = zeroc    ! delay time relative to time0 at which point
                                      ! glf23 transport will be turned on (effective
                                      ! only if glf23 transport is selected. This
                                      ! switch does not determie whether or not glf23 is
                                      ! is used). 
           RETURN
       END SUBROUTINE glf_input_defaults

 




       SUBROUTINE glf_init(glf_setup)      

         USE solcon_gcnmp,                                  ONLY : cm_xke,cm_xki,cm_xkw,itran_max 
         USE dep_var,                                       ONLY : dangrot_drho
  USE MPI_data,                      ONLY :  myid,master,mpiierr  ! 888889999 debug only
         IMPLICIT NONE
         INTEGER(I4B) glf_setup
         REAL(DP) denfactor,dra,rhoa,rmina,qa,drhodrrra,csdaa,delta1,vexb_term1a, &
                  vphima,vexba
              



! -------------------------------------------------------------------------------------
! --  IF confinement model returns flux directly then 
! --  there is a contribution from conduction and convection only
! --  due to the neoclassical background. Iglf_* quantites
! --  are set accordingly here: 
! --------------------------------------------------------------------------------------


          denfactor               = 1.e-19_DP   
          IF(.NOT. ALLOCATED(cm_xke))THEN
             !these arrays are used only if masking of points is in effect
             !here we unmask every point for use in the maccormack solver
                   ALLOCATE(cm_xke(nj))
                   ALLOCATE(cm_xki(nj))
                   ALLOCATE(cm_xkw(nj))
                   cm_xke(:) = 1._DP
                   cm_xki(:) = 1._DP
                   cm_xkw(:) = 1._DP
          ENDIF
         !first defaults:
         iglf_chie = zeroc  ; iglf_chii = zeroc
         iglf_chiv = zeroc  ; iglf_d    = zeroc
         !next turn on if glf is selected
         IF(use_glf23(1) > 0) iglf_d    = 1.0_DP 
         IF(use_glf23(2) > 0) iglf_chie = 1.0_DP
         IF(use_glf23(3) > 0) iglf_chii = 1.0_DP
         IF(use_glf23(5) > 0) iglf_chiv = 1.0_DP
         !next explicitely turn off if glf flux model is selected instead of diffusion model
         ! this model still requires the background flux to be added in since
         ! the glf flux will be zero if there are no unstable modes:
         IF(use_glf23_flux(1) > 0) iglf_d       = zeroc
         IF(use_glf23_flux(2) > 0) iglf_chie    = zeroc
         IF(use_glf23_flux(3) > 0) iglf_chii    = zeroc  !(4) is rbp which is not used  in glf
         IF(use_glf23_flux(5) > 0) iglf_chiv    = zeroc





         freeze_alpha_exp = 0 
         jelc_clamp = 0 ; jion_clamp = 0
         no_w_convection = 0_I4B
         tiny = 1.e-8
         relaxd = 1.0_DP            !disables under relaxation method
         ibtflag_glf     = 1        !switch for effective B field 
         glf23_iglf      = 1        !newest version of glf23 selected
         cparam          = 1.0_DP       ; cparam1 = cparam
         u0 = Permeability ; tmu0 = 2._DP*u0                   !H/m
         charge= -Electron_Charge               ! ABS(qe) coul
         njm1            = nj -1_I4B
         xmassp = Proton_Mass                   ! kg
         ! proton  gyro freq, 1/sec:
         omci            = ABS(charge*mhd_dat%btor/xmassp )      
         irotstab_save = irotstab
         jshoot=0
         igrad=0              !set to 1 for dv method below
         ! ----    set transport flags consistent with itran input:
         ! itport_pt(1)=1 plasma transport on; itport_pt(2:3)=1; electron and ion energy transport
         ! itport_pt(4)=1 ;itport_pt(5)=0  transport vphi with neoclassical determining vtheta
         ! itport_pt(4)=1 ;itport_pt(5)=1  transport vphi and vtheta with fast time scale
         !       neoclassical drag built into vphi and vtheta transport equations...
         !       consult G.M. Staebler
         itport_pt(1)= itran(1)           ! density
         itport_pt(2)= itran(nion+1)      ! Te
         itport_pt(3)= itran(nion+2)      ! Ti
         itport_pt(4)= itran(nion+dp4)    ! tor rot
         itport_pt(5)=0                   ! neoclassical vtheta transport

         rmajor_glf=dischg%rmag           ! major radius at mag axis, R (m)
         btor_exp  = mhd_dat%btor         ! toroidal field, Bt (T)
         arho_iglf = r(nj)                ! toroidal flux at LCFS, rho, (m)


         amassimp_glf = zeroc
         zimp_glf     = zeroc
         IF(nimp .GT. 0)THEN
            DO j = 1,nimp
               amassimp_glf = amassimp_glf + atw(nprim+j)
               zimp_glf    = zimp_glf     + z(1,nprim+j)
            ENDDO
            amassimp_glf = amassimp_glf/nprim    ! mass of impurity
            zimp_glf     = zimp_glf/nprim
         ENDIF



         !         check to make sure gradient scale lengths  are  not zero:
         !          (a zero value will crash glf2d)


         te_glf(:)      = te(:)
         ti_glf(:)      = ti(:)


         !         glf23 single species thermal ion :
         ni_glf(:) = zeroc
         amassgas_glf = zeroc
         DO j = 1,nprim                        ! effective primary ion
            amassgas_glf = amassgas_glf + atw(j)
            ni_glf(:) = ni_glf(:) + en(:,j)
         ENDDO
         amassgas_glf = amassgas_glf/nprim                  ! mass of effective ion species
                                                            ! with density ni_glf
         ne_glf(:) = ene(:)  


         !glf balks if gradient is zero:                
         DO j = 2,nj
            diff = ni_glf(j)- ni_glf(j-1)
            IF(ABS(diff) .LT. 1.e-8) &
                 ni_glf(j)= 0.99999*ni_glf(j-1)
            diff = ne_glf(j)-ne_glf(j-1)
            IF(ABS(diff) .LT. 1.e-8) &
                 ne_glf(j)= 0.99999*ne_glf(j-1)
            diff = te_glf(j)-te_glf(j-1)
            IF(ABS(diff) .LT. 1.e-8) &
                 te_glf(j)= 0.99999*te_glf(j-1)

            diff = ti_glf(j)-ti_glf(j-1)
            IF(ABS(diff) .LT. 1.e-8) &
                 ti_glf(j)= 0.99999*ti_glf(j-1)


            ! rotation is done through angrotp_exp
            
         ENDDO



         !
         !        compute geometric quantities, densities, temperatures, beta, q
         !        note: ni is primary ion density
         !

         btsq = btor_exp*btor_exp
         elong_glf(:) = get_values(dischg%elongxnj)

         DO j=1,nj
            betae_glf(j) = ne_glf(j)*te_glf(j)*joupkev * tmu0 /btsq
            betai_glf(j) = ni_glf(j)*ti_glf(j)*joupkev * tmu0 /btsq
            beta_glf(j)  = betae_glf(j)+betai_glf(j)
            rho_glf(j)   = MAX( r(j)/arho_iglf, tiny )   ! rho_glf in (0.0,1.0]
            tau_glf(j)   = ti_glf(j)/te_glf(j)
            ne_glf(j)    = denfactor*ne_glf(j)              !for input into callglf
            ni_glf(j)    = denfactor*ni_glf(j)              !required in units of 1e19
            ns_glf(j)    = denfactor*enbeam_tot(j)          !total beam density
            zeff_exp(j)  = zeff(j)
            q_glf(j)     = ABS(q(j))               !HSJ sign of q has to be positive
            rmin_glf(j)  = MAX( rminor_r(j), tiny)
            rmaj_glf(j)  = ravg_r(j)
            tot_den=zeroc
            DO jj=1,nprim       ! sum over hydrogenic species
               IF (atw(jj) .LE. 3.1) THEN
                  tot_den=tot_den+en(j,jj)
               END IF
            END DO
            tot_primary_ion_den(j)=tot_den
            tot_thermal_ion_den(j)=tot_den
            DO jj=nprim+1,nion
               tot_thermal_ion_den(j)=tot_thermal_ion_den(j)+en(j,jj)
            END DO
         ENDDO
  
         DO j=2,nj
            drhodr(j)=(rho_glf(j)-rho_glf(j-1))*arho_iglf/ &
                 (rmin_glf(j)-rmin_glf(j-1)+1.e-6)
            drhodrrrho(j)=drhodr(j)*rmin_glf(j)/arho_iglf/rho_glf(j)
            !get d(q/rminor)/drho on zone center grid labelled (0,nj-2)
            dqrm_drho(j-2)= (q_glf(j)/rmin_glf(j)-q_glf(j-1)/rmin_glf(j-1))/(rho_glf(j)-rho_glf(j-1))
         ENDDO
         drhodr(1)=drhodr(2)
         drhodrrrho(2)=drhodrrrho(1)


         DO j=1,nj-1
            drhodr(j)     =(rho_glf(j+1)-rho_glf(j))*arho_iglf/ &
                             (rmin_glf(j+1)-rmin_glf(j)+1.e-6)
            drhodrrrho(j) = drhodr(j)*(rmin_glf(j+1)+rmin_glf(j))/arho_iglf/(rho_glf(j+1)+rho_glf(j))
            !get d(q/rminor)/drho on zone center grid labelled (0,nj-2)
            dqrm_drho(j-1)= (q_glf(j+1)/rmin_glf(j+1)-q_glf(j)/rmin_glf(j))/(rho_glf(j+1)-rho_glf(j))
         ENDDO

         drhodr(nj)=drhodr(nj-1)
         drhodrrrho(nj)=drhodrrrho(nj-1)
         !
         !         compute csda and rhosda : 
         !
         DO j=1,nj
            !            csda_glf(j)=9.79e5*(te_glf(j)*1.e3)**.5/ &
            !               (arho_iglf*100.)/amassgas_glf**.5
            csda_glf(j) = SQRT(te_glf(j)*joupkev/(xmassp*amassgas_glf))/arho_iglf !1/sec
            !            rhosda_glf(j)=((1.02e2*(te_glf(j)*1.e3)**.5)/ &
            !               btor_exp/1.e4)*(amassgas_glf**.5)/(arho_iglf*100.)
            rhosda_glf(j) = csda_glf(j)/omci
         ENDDO
         IF(ibtflag_glf .GT.0 )THEN
            DO j=1,nj
               bteff_glf(j)=btor_exp*rho_glf(j)*arho_iglf/rmin_glf(j)*drhodr(j)
               !              rhosda_glf(j)=((1.02e2*(te_glf(j)*1.e3)**.5)/ &
               !               bteff_glf(j)/1.e4)*(amassgas_glf**.5)/(arho_iglf*100.)

               rhosda_glf(j)= SQRT(te_glf(j)*joupkev/(xmassp*amassgas_glf))/ &
                    (arho_iglf*omci*bteff_glf(j)/mhd_dat%btor)
            ENDDO
         ENDIF
         !
         !         compute normalized gradients, shear, alpha
         !
     sdvm: IF( spline_gsc  == 1 )THEN
            !spline smoothing of profiles input for
            !better behavior in glf23:
            IF(.NOT. ALLOCATED(s_coef)) ALLOCATE(s_coef(nj))
            s_bc = 1        ! spline bc, = 1 means first derivative given at ends
            bc_0 = 0.0      ! first deriv at rhjo =0
            !first deriv at rho =1 :   
            bc_1 = (te_glf(nj)-te_glf(nj-1))/(rho_glf(nj)-rho_glf(nj-1))
            get_scoef = .TRUE. 
            CALL spline_intrp (nj,rho_glf,te_glf,s_coef,get_scoef,s_bc,  &
                 bc_0,bc_1,rho_w,sp,dsp,d2sp)
            get_scoef = .FALSE.
            DO j=2,nj
               rho_w = (rho_glf(j)+rho_glf(j-1))*0.5_DP
               CALL spline_intrp (nj,rho_glf,te_glf,s_coef,get_scoef,s_bc,  &
                    bc_0,bc_1,rho_w,sp,dsp,d2sp)
               zpte_glf(j-1)= -dsp/sp
            ENDDO
            bc_1 = (ti_glf(nj)-ti_glf(nj-1))/(rho_glf(nj)-rho_glf(nj-1))
            get_scoef = .TRUE. 
            CALL spline_intrp (nj,rho_glf,ti_glf,s_coef,get_scoef,s_bc,  &
                 bc_0,bc_1,rho_w,sp,dsp,d2sp)
            get_scoef = .FALSE.
            DO j=2,nj
               rho_w = (rho_glf(j)+rho_glf(j-1))*0.5_DP
               CALL spline_intrp (nj,rho_glf,ti_glf,s_coef,get_scoef,s_bc,  &
                    bc_0,bc_1,rho_w,sp,dsp,d2sp)
               zpti_glf(j-1)= -dsp/sp
            ENDDO

            bc_1 = (ne_glf(nj)-ne_glf(nj-1))/(rho_glf(nj)-rho_glf(nj-1))
            get_scoef = .TRUE. 
            CALL spline_intrp (nj,rho_glf,ne_glf,s_coef,get_scoef,s_bc,  &
                 bc_0,bc_1,rho_w,sp,dsp,d2sp)
            get_scoef = .FALSE.
            DO j=2,nj
               rho_w = (rho_glf(j)+rho_glf(j-1))*0.5_DP
               CALL spline_intrp (nj,rho_glf,ne_glf,s_coef,get_scoef,s_bc,  &
                    bc_0,bc_1,rho_w,sp,dsp,d2sp)
               zpne_glf(j-1)= -dsp/sp
            ENDDO

            bc_1 = (ni_glf(nj)-ni_glf(nj-1))/(rho_glf(nj)-rho_glf(nj-1))
            get_scoef = .TRUE. 
            CALL spline_intrp (nj,rho_glf,ni_glf,s_coef,get_scoef,s_bc,  &
                 bc_0,bc_1,rho_w,sp,dsp,d2sp)
            get_scoef = .FALSE.
            DO j=2,nj
               rho_w = (rho_glf(j)+rho_glf(j-1))*0.5_DP
               CALL spline_intrp (nj,rho_glf,ni_glf,s_coef,get_scoef,s_bc,  &
                    bc_0,bc_1,rho_w,sp,dsp,d2sp)
               zpni_glf(j-1)= -dsp/sp
            ENDDO

            bc_1 = (q(nj)-q(nj-1))/(rho_glf(nj)-rho_glf(nj-1))
            get_scoef = .TRUE. 
            CALL spline_intrp (nj,rho_glf,q,s_coef,get_scoef,s_bc,  &
                 bc_0,bc_1,rho_w,sp,dsp,d2sp)
            get_scoef = .FALSE.
            DO j=2,nj
               rho_w = (rho_glf(j)+rho_glf(j-1))*0.5_DP
               CALL spline_intrp (nj,rho_glf,q,s_coef,get_scoef,s_bc,  &
                    bc_0,bc_1,rho_w,sp,dsp,d2sp)
               dqdrho = -dsp/sp
               drhodq=-(rho_glf(j)+rho_glf(j-1))/(q(j)+q(j-1))
               shat_exp(j)=drhodq*dqdrho
               IF (ABS(shat_exp(j)).LT.1.e-6) shat_exp(j)=1.e-6
            ENDDO

            bc_1 = (beta_glf(nj)-beta_glf(nj-1))/(rho_glf(nj)-rho_glf(nj-1))
            get_scoef = .TRUE. 
            CALL spline_intrp (nj,rho_glf,beta_glf,s_coef,get_scoef,s_bc,  &
                 bc_0,bc_1,rho_w,sp,dsp,d2sp)
            get_scoef = .FALSE.
            DO j=2,nj
               rho_w = (rho_glf(j)+rho_glf(j-1))*0.5_DP
               CALL spline_intrp (nj,rho_glf,beta_glf,s_coef,get_scoef,s_bc,  &
                    bc_0,bc_1,rho_w,sp,dsp,d2sp)
               dbetadrho = dsp/sp
               alpha_exp(j)=-drhodr(j)*q_glf(j)**2*rmajor_glf/arho_iglf* &
                    dbetadrho
               IF(freeze_alpha_exp .EQ.0)THEN
                  alpha_exp_save(j)= alpha_exp(j)
               ELSE IF(freeze_alpha_exp .EQ. 1)THEN
                  alpha_exp(j) = alpha_exp_save(j)
               ELSE
                  IF(shat_exp(j) .LT. 0.0)THEN
                     alpha_exp(j) = MIN(0.999_DP,alpha_exp(j))
                  ENDIF
               ENDIF
            ENDDO
            zpte_glf(nj)=zpte_glf(nj-1)
            zpti_glf(nj)=zpti_glf(nj-1)
            zpne_glf(nj)=zpne_glf(nj-1)
            zpni_glf(nj)=zpni_glf(nj-1)
            shat_exp(1)  = shat_exp(2)
     ELSE  sdvm
            DO j=2,nj
               drhogf(j)=rho_glf(j)-rho_glf(j-1)  !rho_glf is [0,1]
               !forward difference on derivatives to get match with xptor:
               tea = 0.5*(te(j)+te(j-1))
               tia = 0.5*(ti(j)+ti(j-1))
               ane = 0.5*(ne_glf(j)+ne_glf(j-1))
               ani = 0.5*(ni_glf(j)+ni_glf(j-1))
               zpte_glf(j-1)=(-te_glf(j)+te_glf(j-1))/(tea*drhogf(j)) ! 1/LTe
               zpti_glf(j-1)=(-ti_glf(j)+ti_glf(j-1))/(tia*drhogf(j)) ! 1/LTi
               zpne_glf(j-1)=(-ne_glf(j)+ne_glf(j-1))/(ane*drhogf(j)) ! 1/Lne
               zpni_glf(j-1)=(-ni_glf(j)+ni_glf(j-1))/(ani*drhogf(j)) ! 1/Lni

!               zpte_glf(j-1)=ABS((-te_glf(j)+te_glf(j-1))/(tea*drhogf(j))) ! 1/LTe ! 88889999 added ABS()
!               zpti_glf(j-1)=ABS((-ti_glf(j)+ti_glf(j-1))/(tia*drhogf(j))) ! 1/LTi
!               zpne_glf(j-1)=ABS((-ne_glf(j)+ne_glf(j-1))/(ane*drhogf(j))) ! 1/Lne
!               zpni_glf(j-1)=ABS((-ni_glf(j)+ni_glf(j-1))/(ani*drhogf(j))) ! 1/Lni
              
               dqdrho=-(q(j)-q(j-1))/drhogf(j)
               drhodq=-(rho_glf(j)+rho_glf(j-1))/(q(j)+q(j-1))
               shat_exp(j)=drhodq*dqdrho                                   !(r/q)dq/dr
               IF (ABS(shat_exp(j)).LT.1.e-6) shat_exp(j)=1.e-6
               alpha_exp(j)=-drhodr(j)*q_glf(j)**2*rmajor_glf/arho_iglf* &
                    (beta_glf(j)-beta_glf(j-1))/drhogf(j)
               IF(freeze_alpha_exp .EQ.0)THEN
                  alpha_exp_save(j)= alpha_exp(j)
               ELSE IF(freeze_alpha_exp .EQ. 1)THEN
                  alpha_exp(j) = alpha_exp_save(j)
               ELSE
                  IF(shat_exp(j) .LT. 0.0)THEN
                     alpha_exp(j) = MIN(0.999_DP,alpha_exp(j))
                  ENDIF
               ENDIF
            ENDDO

            zpte_glf(nj)=zpte_glf(nj-1)
            zpti_glf(nj)=zpti_glf(nj-1)
            zpne_glf(nj)=zpne_glf(nj-1)
            zpni_glf(nj)=zpni_glf(nj-1)
            shat_exp(1)  = shat_exp(2)

    ENDIF sdvm



         drhogf(1)=1.e-6
         IF(freeze_alpha_exp .EQ.0)THEN
            alpha_exp(1)=1.e-6
            alpha_exp_save(1)=alpha_exp(1)
         ELSE IF(freeze_alpha_exp .EQ. 1)THEN
            alpha_exp(1) = alpha_exp_save(1)
         ENDIF
         !
         !    compute plasma toroidal angular velocity (angrotp) in 1/s
         !    and parallel (vpar), perpendicular(vper), and toroidal (vphi)
         !    velocities in m/s
            vstar_sign=-1.
            corot=1._DP ! this is sign for co/counter injection.
            !we assume that value of central rotation gives proper direction.
            !because we assume that angrt is a signed quantity we do not
            ! mutiply it by corot below
            IF(angrot(1) .LT. zeroc)corot = -1._DP
            glf_corot = corot
         DO j=1,nj        ! loop over the full grid
            fc=1-1.46*(rmin_glf(j)/rmaj_glf(j))**0.5+ &
                 0.46*(rmin_glf(j)/rmaj_glf(j))**1.5
            akappa1=0.8839*fc/(0.3477+0.4058*fc)
            akappa2=(1.-fc)/(1.+1.1671*fc)
            alpha_neo=-akappa1+1.   !alpha1neo in Waltz ref.
            !
            fcglf(j)=fc
            ak1glf(j)=akappa1
            ak2glf(j)=akappa2
            aneoglf(j)=alpha_neo
            !

            pgeo_local=drhodr(j)
            rdrho(j)=rmin_glf(j)/arho_iglf/rho_glf(j)
            ! correct for diamagnetic terms, see Waltz,Phys Plas,4(7)1997,pg2482, eq.(16c)
            ! corrected to divide by (r/Rq) instead of multiply .
!            angrotp_exp(j)=angrot(j)+ &
!                 akappa2*3./2.*csda_glf(j)*zpti_glf(j)* &
!                 tau_glf(j)/rho_glf(j)*q_glf(j)* &
!                 rhosda_glf(j)*pgeo_local/rdrho
!            angrotp_exp(j)=corot*angrotp_exp(j)
             diamag_term(j)  = corot*akappa2*3./2.*csda_glf(j)*zpti_glf(j)* &
                 tau_glf(j)/rho_glf(j)*q_glf(j)* &
                 rhosda_glf(j)*pgeo_local/rdrho(j)

            angrotp_exp(j) =  angrot(j)+ diamag_term(j)  ! no corot factor in angrot (see above )
            IF (angrot(j).EQ. zeroc )THEN
               angrotp_exp(j)  = zeroc
               diamag_term(j)  = zeroc
            ENDIF
            vphim_glf(j)      = rmajor_glf*angrotp_exp(j)  ! this is eq 16c corrected for
                                                           ! division by r/q (16c has this factor 
                                                           ! as a multipier
            !as below note
            !determine Vexb (= ve)
            !orginal form (by J.Kinsey)
            !            ve_glf(j)=-tau_glf(j)*csda_glf(j)*arho_iglf*
            !     .         rhosda_glf(j)*(zpni_glf(j)+zpti_glf(j)+(1.-
            !     .         rdrho*rho_glf(j)*arho_iglf/rmajor_glf/q_glf(j))*
            !     .         (alpha_neo-1.)*zpti_glf(j))*vstar_sign*pgeo_local

            !set ve to value used in callglf2d  HSJ 05/0103:
            ve_glf(j)=-tau_glf(j)*csda_glf(j)*arho_iglf* &
                 rhosda_glf(j)*(zpni_glf(j) + &
                 alpha_neo*zpti_glf(j))*vstar_sign*pgeo_local
            vexb_base(j) = ve_glf(j)

            !for either form of ve above we now
            !add sigma*(rminor/(Dischg%rmajor*q))Vtoroidal to Vexb.
            !Note that sigma is assumed  +1 here because the
            !toroidal rotation velocity is a signed quantitiy.
            !This term is identical to eq. 16a in Waltz,
            !phys plas, 4,7,97,2486. HSJ:
!            ve_glf(j)=ve_glf(j)-rdrho*rho_glf(j)* &
!                 arho_iglf/rmajor_glf/q_glf(j)*vphim_glf(j)
            vexb_term1(j) = -corot*rdrho(j)*rho_glf(j)* &
                 arho_iglf/rmajor_glf/q_glf(j)
            ve_glf(j) = vexb_base(j) + vexb_term1(j)*vphim_glf(j)

            vparm_term1(j) = vstar_sign*((1.-alpha_neo)*zpti_glf(j))* &
                 tau_glf(j)*csda_glf(j)*arho_iglf* &
                 rhosda_glf(j)*pgeo_local*rho_glf(j)*arho_iglf &
                 /rmajor_glf/q_glf(j)*rdrho(j)
             vparm_glf(j) = vphim_glf(j) + vparm_term1(j)

!            vparm_glf(j)= vphim_glf(j) + &
!                 vstar_sign*((1.-alpha_neo)*zpti_glf(j))* &
!                 tau_glf(j)*csda_glf(j)*arho_iglf* &
!                 rhosda_glf(j)*pgeo_local*rho_glf(j)*arho_iglf &
!                 /rmajor_glf/q_glf(j)*rdrho


            vperm_term1(j)=-tau_glf(j)*csda_glf(j)*arho_iglf* &
                 rhosda_glf(j)*((1.-rdrho(j)*rho_glf(j)*arho_iglf/ &
                 rmajor_glf/q_glf(j))*(alpha_neo-1.)*zpti_glf(j))* &
                 vstar_sign*pgeo_local
!            vperm_glf(j) =vperm_term1(j) - rdrho*rho_glf(j)* &
!                 arho_iglf/rmajor_glf/q_glf(j)* &
!                 rmajor_glf*angrotp_exp(j)
            vperm_glf(j) = vperm_term1(j) + vexb_term1(j)*vphim_glf(j)
            ! note vpolm_glf Is not used in calculartions but is set for output here
            ! approximate Bp/B with rho/qR0:
            vpolm_glf(j) = vperm_glf(j) -  r(j)/(q_glf(j)*rmajor_glf)*vparm_glf(j)
         ENDDO

            IF(myid == master) THEN
!              WRITE(173,FMT='("j,vparm,vphim,vpol,vperm,vperm_term1,vexb_term1,diamag_term= ")')
              DO j=1,nj
!                 WRITE(173,FMT='(x,i3,x,7(x,1pe12.4))')j,vparm_glf(j),vphim_glf(j),vpolm_glf(j),vperm_glf(j),vperm_term1(j),vexb_term1(j),diamag_term(j)   ! 8888889999
              ENDDO
!              WRITE(173,FMT='("j,rho,akappa2, zpti_glf,ti,rdrho =")')
              DO j=1,nj
 !               WRITE(173,FMT='(x,i3,5(x,1pe12.4))')j,rho_glf(j),ak2glf(j), ti_glf(j),zpti_glf(j),rdrho(j)
              ENDDO
!              WRITE(173,FMT='("rho,tau,q,rhosda,pgeo =")')
              DO j=1,nj
 !                WRITE(173,FMT='(x,i3,5(x,1pe12.4))')j,rho_glf(j),tau_glf(j),q_glf(j),rhosda_glf(j),drhodr(j)  ! 8888899999
              ENDDO
           ENDIF
 !----------------------------------------------------------------
 ! compute ExB and parallel velocity shear rates
 ! Use perturbative formulas for egamma_exp and gamma_p_exp
 ! construct rho derivatives first:
 !---------------------------------------------------------------
         DO j=1,nj-1
            dra                    =  rho_glf(j+1)             - rho_glf(j)
            rhoa                   =  0.5_DP*(rho_glf(j+1)     + rho_glf(j))
            rmina                  =  0.5_DP*(rmin_glf(j+1)    + rmin_glf(j))
            qa                     =  0.5_DP*(q_glf(j+1)       + q_glf(j))
            drhodrrra              =  0.5_DP*(drhodrrrho(j+1)  + drhodrrrho(j))
            drhodrrra              =  drhodrrrho(j)
            csdaa                  =  0.5_DP*(csda_glf(j+1)    + csda_glf(j))
            vexb_term1a            =  0.5_DP*(vexb_term1(j+1)  + vexb_term1(j))
            vphima                 =  0.5_DP*(vphim_glf(j+1)   + vphim_glf(j))
            vexba                  =  0.5_DP*(ve_glf(j+1)      + ve_glf(j))
            delta1                 =  rhoa*drhodrrra/rmina/csdaa

            dvexb_base_drho(j-1)   =  (vexb_base(j+1)          - vexb_base(j))/dra
            dvphim_drho(j-1)       =  (vphim_glf(j+1)          - vphim_glf(j))/dra
            ddiamag_drho(j-1)      =  (diamag_term(j+1)        - diamag_term(j))/dra
            dvexb_term1_drho(j-1)  =  (vexb_term1(j+1)         - vexb_term1(j))/dra
            dvpar_term1_drho(j-1)  =  (vparm_term1(j+1)        - vparm_term1(j))/dra

            dangrot_drho(j)        =  (angrot(j+1)             - angrot(j))/dra !dangrot_drho(1..nj-1)
            rdrho(j)=rmin_glf(j)/arho_iglf/rho_glf(j)
            rdrho_local_p1=rmin_glf(j+1)/arho_iglf/rho_glf(j+1) 
            ! originals:
            egamma_exp(j)=drhodrrrho(j)*(rho_glf(j)+rho_glf(j+1))/ &
                 (q_glf(j)+q_glf(j+1))*(ve_glf(j+1)*q_glf(j+1)/ &
                 rho_glf(j+1)/rdrho_local_p1-ve_glf(j)*q_glf(j)/ &
                 rho_glf(j)/rdrho(j))/(rho_glf(j+1)-rho_glf(j))/ &
                 arho_iglf/csda_glf(j)
 
            gamma_p_exp(j)=-drhodr(j)*(vparm_glf(j+1)-vparm_glf(j))/ &
                 (rho_glf(j+1)-rho_glf(j))/arho_iglf/csda_glf(j)

            ! zone center versions of the above:
            egamma_exp_zc(j-1) = delta1 *(dvexb_base_drho(j-1) &
                 + vexb_term1a*dvphim_drho(j-1) + vphima*dvexb_term1_drho(j-1)  &
                 + (vexba*rmina/qa)*dqrm_drho(j-1))
            gamma_p_exp_zc(j-1)= -delta1 * (dvpar_term1_drho(j-1)   &
                 +  rmajor_glf * dangrot_drho(j) + corot*ddiamag_drho(j-1))

         ENDDO

         IF(irotstab_save .EQ. -1 .AND. init_egamma_exp .EQ. 0 )THEN
            egamma_exp_initial(:) = egamma_exp(:)
            init_egamma_exp = 1
         ENDIF
         IF(irotstab_save .EQ. -1) &
              egamma_exp(:) = egamma_exp_initial(:)


         IF(irotstab_save .EQ. -1)irotstab = 0 !passed to glf23 as 0
         !
         !     for impurity eqns (limp_glf > 0)
         !
         idengrad_glf=2  ! simple dilution
         IF(limp_glf.GT.0) THEN
            WRITE(6,*) 'Turning on impurity dynamics in GLF23 model ...'
            jroot_glf=12
            idengrad_glf=3
         ENDIF
         !

         !

         njm1      = nj-1
         glf_setup = 1           ! prevents further calls to this routine

        glf_e_flux(:,:,:) = zeroc                     ! fluxes to be loaded by glf
        glf_p_flux(:,:,:) = zeroc
        glf_m_flux(:,:,:) = zeroc


        IF( .NOT. ALLOCATED(glf_typ_val))THEN
            ALLOCATE(glf_typ_val(-itran_max:itran_max))
            glf_typ_val(:)   = zeroc
            glf_typ_val(1)   = MAXVAL(ABS(ni_glf))/denfactor  ! densities are in units of denfactor
            glf_typ_val(2)   = MAXVAL(ABS(te_glf))
            glf_typ_val(3)   = MAXVAL(ABS(ti_glf))
            glf_typ_val(4)   = zeroc !  RBp
            glf_typ_val(5)   = MAXVAL(ABS(angrot))

            glf_typ_val(-1)  = ABS((ni_glf(2) - ni_glf(1))/dr(1))/denfactor
            glf_typ_val(-2)  = ABS((te_glf(2) - te_glf(1))/dr(1))
            glf_typ_val(-3)  = ABS((ti_glf(2) - ti_glf(1))/dr(1))
            glf_typ_val(-4)  = zeroc ! RBp
            glf_typ_val(-5)  = ABS((angrot(2) - angrot(1))/dr(1)) 
            DO j=1 ,nj-2 
                  glf_typ_val(-1)  = MAX(glf_typ_val(-1),ABS((ni_glf(j+1) - ni_glf(j))/dr(j+1)/denfactor))
                  glf_typ_val(-2)  = MAX(glf_typ_val(-2),ABS((te_glf(j+1) - te_glf(j))/dr(j+1)))
                  glf_typ_val(-3)  = MAX(glf_typ_val(-3),ABS((ti_glf(j+1) - ti_glf(j))/dr(j+1)))
                  glf_typ_val(-5)  = MAX(glf_typ_val(-5),ABS((angrot(j+1) - angrot(j))/dr(j+1)))
            ENDDO

        ENDIF

       RETURN

       END SUBROUTINE glf_init



       SUBROUTINE glf_flux_equiv
!---------------------------------------------------------------------------------
! -- glf23 can be run in flux mode where fluxes are returned directly 
! -- or it can be run in the chi mode where chi's are returned and fluxes
! -- are calculated from that using Fick's law. (Note that the grad rho
! -- factor is included in the chis returned from glf_driver_single).
! -- In the latter case the fluxes  output directly by glf23:
! --                niflux_gf_m
! --                nimpflux_gf_m
! --                tiflux_gf_m
! --                teflux_gf_m
! --                vparflux_gf_m
! --                vperflux_gf_m
! --                vphiflux_gf_m
! --                etgflux_gf_m
! -- are set but are not used. This subroutine  overwrites these values 
! -- with equivalent fluxes based on Ficks law assumption with the
! -- glf23 determiend d's (in docefp_glf23)
!---------------------------------------------------------------------------------
      USE nrtype,                                             ONLY : DP,I4B
      USE plasma_properties,                                  ONLY : diffuse
      USE solcon_gcnmp,                                       ONLY : ntot,use_glf23_flux,dudr

      USE glf23_gcnmp,                                        ONLY : niflux_gf_m, nimpflux_gf_m,    &
                                                                     tiflux_gf_m,teflux_gf_m,       &
                                                                     vperflux_gf_m,vparflux_gf_m,   &
                                                                     vphiflux_gf_m,etgflux_gf_m,    &
                                                                     anrate_m,anrate2_m,            &
                                                                     anfreq_m,anfreq2_m,            &
                                                                     glf_anrate,glf_anrate2,        &
                                                                     glf_anfreq,glf_anfreq2

 
      IMPLICIT NONE
      INTEGER(i4B) j,i,k
      REAL(DP) epsilon,enpa

      
      
          epsilon = 1.e-10_DP
          DO j = 1,nj-1
             k = 1
             IF(use_glf23_flux(k) == 1)THEN
                !fluxes were loaded by glf23, get an equivalent d
               DO k = 1,nion
                  IF(ABS(dudr(k,j)) .GT. epsilon )THEN
                     dcoefp_glf23(k,j)  = -niflux_gf_m(j)/dudr(k,j) 
                  ELSE
                     dcoefp_glf23(k,j)  = zeroc
                  ENDIF
               ENDDO
            ELSE
               !d's were loaded get an equivalent flux
                enpa = zeroc
                DO i=1,nion
                   enpa = enpa + 0.5_DP*(en(j+1,i)+en(j,i))
                ENDDO
                  IF(vp_eff_m(j) .NE.  zeroc)THEN
                     niflux_gf_m(j) =   enpa*vp_eff_m(j)
                  ELSE
                     niflux_gf_m(j) = - dcoefp_glf23(1,j)*dudr(1,j) 
                  ENDIF
            ENDIF


            k=nion+1
            IF(use_glf23_flux(k-nion+1) == 1)THEN
                  IF(ABS(dudr(k,j)) .GT. epsilon)THEN 
                     dcoefp_glf23(k,j)   = -teflux_gf_m(j)/dudr(k,j)
                  ELSE
                     dcoefp_glf23(k,j)    = zeroc
                  ENDIF
            ELSE
                  IF(ve_eff_m(j) .NE. zeroc) THEN
                     tea = 0.5*(te(j+1)+te(j))
                     enea=0.5*(ene(j+1)+ene(j))
                     teflux_gf_m(j) = ve_eff_m(j)*enea*tea
                  ELSE
                     teflux_gf_m(j) =  - dcoefp_glf23(k,j)*dudr(k,j)
                  ENDIF
            ENDIF


            k=nion+2
            IF(use_glf23_flux(k-nion+1) == 1)THEN
               IF(ABS(dudr(k,j)) .GT. epsilon)THEN 
                   dcoefp_glf23(k,j)   = -tiflux_gf_m(j)/dudr(k,j)
               ELSE
                   dcoefp_glf23(k,j)    = zeroc
               ENDIF

            ELSE
               IF(vi_eff_m(j) .NE. zeroc) THEN
                     tia = 0.5*(ti(j+1)+ti(j))
                     tiflux_gf_m(j) =  vi_eff_m(j)*enpa*tia
               ELSE
                     tiflux_gf_m(j) =  - dcoefp_glf23(k,j)*dudr(k,j)
               ENDIF
            ENDIF

            k = nion + dp4
            IF(use_glf23_flux(k-nion+1) == 1)THEN
               IF(ABS(dudr(k,j)) .GT. epsilon)THEN 
                   dcoefp_glf23(k,j)   = -vphiflux_gf_m(j)/dudr(k,j)   ! vphiflux_gf in kg/sec^2 
               ELSE
                   dcoefp_glf23(k,j)    = zeroc
               ENDIF
            ELSE
                 vphiflux_gf_m(j)    = - dcoefp_glf23(k,j)*dudr(k,j)   ! (kg m/sec)*(1/(m sec) = kg/sec^2 momentum flux 
                 !leave the following at the output values from glf23
                 ! since we dont have corresponding d's:
!                 vperflux_gf_m(j)    = zeroc
!                 vparflux_gf_m(j)    = zeroc
!                 etgflux_gf_m(j)     = zeroc
            ENDIF
         ENDDO

 
         ! glf_driver_single loads dcoefp_glf23:
           diffuse%dcoef_glf23(:,:) = dcoefp_glf23(:,:)



        glf_anrate(:) = anrate_m(:)
        glf_anrate2(:) = anrate2_m(:)
        glf_anfreq(:)  = anfreq_m(:)
        glf_anfreq2(:) = anfreq2_m(:)

            ! lhs [0,nj-1] == rhs [1:nj]   
            glf_p_flux(0:nj-1,2,0)   = niflux_gf_m(1:nj) 
            glf_p_flux(0:j-1,3,0)    = nimpflux_gf_m(1:nj)
            glf_e_flux(0:j-1,1,0)    = teflux_gf_m(1:nj)
            glf_e_flux(0:j-1,2,0)    = tiflux_gf_m(1:nj)
            glf_m_flux(0:j-1,2,0)    = vphiflux_gf_m(1:nj)
           
!         CALL dump_test(0) 

   END SUBROUTINE glf_flux_equiv




  SUBROUTINE dump_test(flag)
!---------------------------------------------------------------------------------------------
! --
!--------------------------------------------------------------------------------------------
      USE MPI_data,                                  ONLY : myid,master
      USE glf23_gcnmp,                               ONLY : niflux_gf_m,nimpflux_gf_m,teflux_gf_m, &
                                                            tiflux_gf_m,etgflux_gf_m,anrate_m,     &
                                                            glf_anrate,glf_anrate2,anrate2_m,      &
                                                            anfreq_m,glf_anfreq,anfreq2_m,         &
                                                            glf_anfreq2,gamma_net_i,gamma_net_e,   &
                                                            glf_e_flux,glf_m_flux
                    
      INTEGER(I4b) flag
      IF(myid == master)THEN
         IF(flag == 0)WRITE(888,FMT='("===================called from glf_flux_equiv======================")')
         IF(flag == 100)WRITE(888,FMT='("===================called from set_profiless before refactor======================")')
         IF(flag == 200)WRITE(888,FMT='("===================called from set_profiles after refactor======================")')
         WRITE(888,FMT='("niflux_gf_m =")')
         WRITE(888,1)niflux_gf_m
         WRITE(888,FMT='("nimpflux_gf_m =")')
         WRITE(888,1)nimpflux_gf_m
         WRITE(888,FMT='("teflux_gf_m =")')
         WRITE(888,1)teflux_gf_m
         WRITE(888,FMT='("tiflux_gf_m =")')
         WRITE(888,1)tiflux_gf_m
         WRITE(888,FMT='("etglux_gf_m =")')
         WRITE(888,1)etgflux_gf_m
         WRITE(888,FMT='("anrate_m =")')
         WRITE(888,1)anrate_m
         WRITE(888,FMT='("glf_anrate =")')
         WRITE(888,FMT='("size glf_anrate in dump_test",i5)'),SIZE(glf_anrate)
         WRITE(888,1)glf_anrate
         WRITE(888,FMT='("anrate2_m =")')
         WRITE(888,1)anrate2_m
         WRITE(888,FMT='("glf_anrate2 =")')
         WRITE(888,FMT='("size glf_anrate2 in dump_test",i5)'),SIZE(glf_anrate2) 
         WRITE(888,1)glf_anrate2
         WRITE(888,FMT='("anfreq_m =")')
         WRITE(888,1)anfreq_m
         WRITE(888,FMT='("glf_anfreq =")')
         WRITE(888,FMT='("size glf_anrfreq in dump_test",i5)'),SIZE(glf_anfreq) 
         WRITE(888,1)glf_anfreq
         WRITE(888,FMT='("anfreq2_m =")')
         WRITE(888,1)anfreq2_m
         WRITE(888,FMT='("glf_anfreq2 =")')
         WRITE(888,FMT='("size glf_anfreq2 in dump_test",i5)'),SIZE(glf_anfreq2) 
         WRITE(888,1)glf_anfreq2
         WRITE(888,FMT='("gamma_net_i =")')
         WRITE(888,1)gamma_net_i
         WRITE(888,FMT='("gamma_net_e =")')
         WRITE(888,1)gamma_net_e
         WRITE(888,FMT='("glf_e_flux =")')
         WRITE(888,1)glf_e_flux(:,1,0)
         WRITE(888,FMT='("glf_m_flux  =")')
         WRITE(888,1)glf_m_flux(:,1,0)
         WRITE(888,FMT='("glf_m_flux prim ions  =")')
         WRITE(888,1)glf_m_flux(:,2,0)
         WRITE(888,FMT='("glf_m_flux imp ions  =")')
         WRITE(888,1)glf_m_flux(:,3,0)


1        FORMAT(5(1pe12.4,2x))


     ENDIF

     RETURN

  END  SUBROUTINE dump_test

  SUBROUTINE set_glf_eff_chi
  !----------------------------------------------------------------------
  ! define effective conductivity for glf23 model
  ! used for output to statefile
  !----------------------------------------------------------------------
    USE MPI_data,                      ONLY :  myid,master,mpiierr  ! 888889999 debug only
    REAL(DP) grade,gradi,enea
    REAL(DP) enia
    INTEGER(I4B) j,k

         glf_eff_chi_output(:,:) = zeroc

         DO i=1,ntot
            DO j=2,SIZE(dudr,2)
                  IF(i == nion+1)THEN
                     enea = (ene(j)+ ene(j+1))*0.5_DP
                     grade = dudr(i,j)*enea
                     IF(ABS(grade) .GT. 0.00001_DP)&
                     glf_eff_chi_output(j,i)  = -teflux_gf_m(j)/grade
                  ELSEIF(i == nion+2)THEN
                     enia = zeroc
                     DO k=1,nion
                        enia = enia+en(j+1,k)+en(j,k)
                     ENDDO
                     enia = enia*0.5_DP
                     gradi = dudr(i,j)*enia
                     IF(ABS(gradi) .GT. 0.00001_DP)&
                     glf_eff_chi_output(j,i)  = -tiflux_gf_m(j)/gradi
                  ELSE
                     glf_eff_chi_output(j,i)  = zeroc ! not yet programmed
                  ENDIF
            ENDDO
         ENDDO

     RETURN

  END  SUBROUTINE set_glf_eff_chi

 


 
  

      SUBROUTINE  glf_limit_values (   &
        jmm,jmaxm,                     &
        diffnem,                       & ! ion plasma diffusivity in m**2/sec
        chietem,                       & ! electron ENERGY diffuivity in m**2/sec
        chiitim,                       & ! ion      ENERGY diffuivity in m**2/sec
        etaphim,                       & ! toroidal velocity diffusivity in m**2/sec
        etaparm,                       & ! parallel velocity diffusivity in m**2/sec
        etaperm,                       & ! perpendicular velocity diffusivity in m**2/sec
        exchm,                         & ! turbulent electron to ion ENERGY exchange in MW/m**3
        diff_m,                        & ! 0:jmaxm values (diffnem)
        chie_m,                        & ! 0:jmaxm values (chietem)
        chii_m,                        & ! 0:jmaxm values (chiitim)
        etaphi_m,                      & ! 0:jmaxm values (etaphim)
        etapar_m,                      & ! 0:jmaxm values (etaparm)
        etaper_m,                      & ! 0:jmaxm values (etaperm)
        exch_m,                        & ! 0:jmaxm values (exchm)
        egamma_m,                      & !0:jmaxm exb shear rate in units of local csda_m
        egamma_d,                      & !0:jmaxm exb shear rate delayed by i_delay steps
        gamma_p_m,                     & !0:jmaxm par. vel. shear rate in units of local csda_m
        anrate_m,                      & !0:jmaxm leading mode rate in units of local csda_m
        anrate2_m,                     & !0:jmaxm 2nd mode rate in units of local csda_m
        anfreq_m,                      & !0:jmaxm leading mode frequency
        anfreq2_m,                     & !0:jmaxm 2nd mode frequency
        niflux_gfl,                    &  ! ion particle flux
        niflux_gf_m,                   &  ! 0:jmaxm
        nimpflux_gfl,                  &  ! impurity particle flux
        nimpflux_gf_m,                 &
        tiflux_gfl,                    & ! ion  energy flux
        tiflux_gf_m,                   &
        teflux_gfl,                    & ! electron  energy flux
        teflux_gf_m,                   &
        vparflux_gfl,                  & ! parallel viscous stress
        vparflux_gf_m,                 &
        vperflux_gfl,                  & ! pependicularr viscous sress
        vperflux_gf_m,                 &
        vphiflux_gfl,                  & ! toroidal viscous stress
        vphiflux_gf_m,                 &
        etgflux_gfl,                   & ! contribution from high k etg modes to electron energy flux
        etgflux_gf_m,                  & ! electron temperature gradient mode flux
        gamma_net_i,                   & ! net growth rate ion modes only
        gamma_net_e                    & ! net growth rates electron modes only
        )
!-----------------------------------------------------------------------------------------
! apply lower and upper limits to glf output quantites
!------------------------------------------------------------------------------------------
      IMPLICIT NONE 
      INTEGER(I4B) jmaxm,jmm
      REAL(DP)                                  &
            diffnem,chietem, chiitim, etaphim,  &
            etaparm, etaperm, exchm,            &
            tiflux_gfl,teflux_gfl,              &
            vparflux_gfl,                       &
            vperflux_gfl,vphiflux_gfl,          &
            etgflux_gfl,niflux_gfl,nimflux_gfl

      REAL(DP)                                  & ! zone edge grid [0:nj-1]
            diff_m(0:jmaxm),chie_m(0:jmaxm),    &
            chii_m(0:jmaxm),etaphi_m(0:jmaxm),  &
            etapar_m(0:jmaxm),etaper_m(0:jmaxm),&
            exch_m(0:jmaxm),gamma_m(0:jmaxm),   &
            egamma_m(0:jmaxm),                  &
            egamma_d(0:jmaxm),                  &
            gamma_p_m(0:jmaxm),                 &
            anrate_m(0:jmaxm),                  &
            anrate2_m(0:jmaxm),                 &
            anfreq_m(0:jmaxm),                  &
            anfreq2_m(0:jmaxm),                 &
            niflux_gf_m(0:jmaxm),               &
            nimpflux_gf_m(0:jmaxm),             &
            tiflux_gf_m(0:jmaxm),               &
            teflux_gf_m(0:jmaxm),               &
            vparflux_gf_m(0:jmaxm),             &
            vperflux_gf_m(0:jmaxm),             &
            vphiflux_gf_m(0:jmaxm),             &
            etgflux_gf_m(0:jmaxm),              &
            gamma_net_i(0:jmaxm),               &
            gamma_net_e(0:jmaxm)


      REAL(DP)                                  &
            diffnem_glf_min,diffnem_glf_max,    &
            chie_glf_max,chie_glf_min,          &
            chii_glf_max,chii_glf_min,          &
            etaphim_glf_max,etaphim_glf_min,    &
            etaparm_glf_max,etaparm_glf_min,    &
            etaperm_glf_max,etaperm_glf_min,    &
            niflux_glf_max,niflux_glf_min,      &
            efluxe_glf_max,efluxe_glf_min,      &
            efluxi_glf_max,efluxi_glf_min,      &
            nimpflux_gfl

            diffnem_glf_max  = 25._DP  ;   diffnem_glf_min   = -diffnem_glf_max   ! m^2/sec
            chie_glf_max     = 25._DP  ;   chie_glf_min      = -chie_glf_max      ! m^2/sec
            chii_glf_max     = 25._DP  ;   chii_glf_min      = -chii_glf_max      ! m^2/sec
            etaphim_glf_max  = 25._DP  ;   etaphim_glf_min   = -etaphim_glf_max   ! m^2/sec
            etaparm_glf_max  = 25._DP  ;   etaparm_glf_min   = -etaparm_glf_max   ! m^2/sec
            etaperm_glf_max  = 25._DP  ;   etaperm_glf_min   = -etaperm_glf_max   ! m^2/sec
            niflux_glf_max   =  1.e24_DP   ! fix this
            niflux_glf_min   = -1.e24_DP
            efluxe_glf_max   =  ABS(chie_glf_max*ene(jmm+1)*dudr(nion+1,jmm+1))
            efluxe_glf_min   = -ABS(chie_glf_min*ene(jmm+1)*dudr(nion+1,jmm+1))
 !           efluxe_glf_min   = 0.0_DP
            efluxi_glf_max   =  ABS(chii_glf_max*en(jmm+1,1)*dudr(nion+2,jmm+1))
            efluxi_glf_min   = -ABS(chii_glf_min*en(jmm+1,1)*dudr(nion+2,jmm+1))
 !           efluxi_glf_min   = 0.0_DP

            diffnem        = MIN(diffnem,diffnem_glf_max)
            diffnem        = MAX(diffnem,diffnem_glf_min)

            chietem        = MIN(chietem,chie_glf_max)
            chietem        = MAX(chietem,chie_glf_min)

            chiitim        = MIN(chiitim,chii_glf_max)
            chiitim        = MAX(chiitim,chii_glf_min)

            etaphim        = MIN(etaphim,etaphim_glf_max) 
            etaphim        = MAX(etaphim,etaphim_glf_min) 
                     
            etaparm        = MIN(etaparm,etaparm_glf_max) 
            etaparm        = MAX(etaparm,etaparm_glf_min)   
                     
            etaperm        = MIN(etaperm,etaperm_glf_max) 
            etaperm        = MAX(etaperm,etaperm_glf_min)   

            niflux_gfl     = MIN(niflux_gfl,niflux_glf_max)
            niflux_gfl     = MAX(niflux_gfl,niflux_glf_min)

            nimpflux_gfl   = MIN(nimpflux_gfl,niflux_glf_max)
            nimpflux_gfl   = MAX(nimpflux_gfl,niflux_glf_min)
            tiflux_gfl     = MIN(tiflux_gfl,efluxi_glf_max)
            tiflux_gfl     = MAX(tiflux_gfl,efluxi_glf_min)
            teflux_gfl     = MIN(teflux_gfl,efluxe_glf_max)
            teflux_gfl     = MAX(teflux_gfl,efluxe_glf_min)
            etgflux_gfl    = MIN(etgflux_gfl,efluxe_glf_max)
            etgflux_gfl    = MAX(etgflux_gfl,efluxe_glf_min)
 RETURN
            DO j=0,jmaxm
               diff_m(j)           = MIN(diff_m(j),diffnem_glf_max)
               diff_m(j)           = MAX(diff_m(j),diffnem_glf_min)
               chie_m(j)           = MIN(chie_m(j),chie_glf_max)
               chie_m(j)           = MAX(chie_m(j),chie_glf_min)
               chii_m(j)           = MIN(chii_m(j),chii_glf_max)
               chii_m(j)           = MAX(chii_m(j),chii_glf_min)
               etaphi_m(j)         = MIN(etaphi_m(j),etaphim_glf_max)
               etaphi_m(j)         = MAX(etaphi_m(j),etaphim_glf_min)
               etapar_m(j)         = MIN(etapar_m(j),etaparm_glf_max)
               etapar_m(j)         = MAX(etapar_m(j),etaparm_glf_min)
               etaper_m(j)         = MIN(etaper_m(j),etaperm_glf_max)
               etaper_m(j)         = MAX(etaper_m(j),etaperm_glf_min)
               niflux_gf_m(j)      = MIN(niflux_gf_m(j),niflux_glf_max)
               niflux_gf_m(j)      = MAX(niflux_gf_m(j),niflux_glf_min)
               nimpflux_gf_m(j)    = MIN(nimpflux_gf_m(j),niflux_glf_max)
               nimpflux_gf_m(j)    = MAX(nimpflux_gf_m(j),niflux_glf_min)
               tiflux_gf_m(j)      = MIN(tiflux_gf_m(j),efluxi_glf_max)
               tiflux_gf_m(j)      = MAX(tiflux_gf_m(j),efluxi_glf_min)
               teflux_gf_m(j)      = MIN(teflux_gf_m(j),efluxe_glf_max)
               teflux_gf_m(j)      = MAX(teflux_gf_m(j),efluxe_glf_min)
               etgflux_gf_m(j)     = MIN(etgflux_gf_m(j),efluxe_glf_max)
               etgflux_gf_m(j)     = MAX(etgflux_gf_m(j),efluxe_glf_min) 
            ENDDo

  END  SUBROUTINE glf_limit_values 


   END MODULE glf_tport_gcnmp



