
      subroutine allocate_glf
c-----------------------------------------------------------------------
c -- create arrays needed for glf23 
c ----------------------------------------------------------------------
c
      USE param 
      USE allocate_err
      USE glf23
      USE numbrs
      implicit  integer (i-n), real*8 (a-h, o-z)
c      include 'param.i'
c      include 'numbrs.i'
      integer pdelay
      parameter (pdelay =10)
c     pdelay is sized according to i_delay
c     the value 10 is hardwired in glf23 for dimension of egamma_d
      if(pdelay .lt. i_delay)
     . call STOP('subroutine ALLOCATE_GLF: i_dealy problem',0)

      
        allocate (te_glf(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("te_glf,allocate_glf",0,istat)

        allocate (ti_glf(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("ti_glf,allocate_glf",0,istat)

        allocate (ne_glf(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("ne_glf,allocate_glf",0,istat)

        allocate (ni_glf(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("ni_glf,allocate_glf",0,istat)

        allocate (q_glf(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("q_glf,allocate_glf",0,istat)

        allocate (shat_exp(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("shat_exp,allocate_glf",0,istat)

        allocate (elong_glf(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("elong_glf,allocate_glf",0,istat)

        allocate (zeff_exp(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("zeff_exp,allocate_glf",0,istat)

        allocate (rho_glf(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("rh_glf,allocate_glf",0,istat)

        allocate (grho1_glf(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("grho1_glf,allocate_glf",0,istat)

        allocate (grho2_glf(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("grho2_glf,allocate_glf",0,istat)

        allocate (gradte_glf(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("gradte_glf,allocate_glf",0,istat)

        allocate (zpte_glf(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("zpte_glf,allocate_glf",0,istat)

        allocate (gradti_glf(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("gradti_glf,allocate_glf",0,istat)

        allocate (zpti_glf(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("zpti_glf,allocate_glf",0,istat)

        allocate (gradne_glf(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("gradne_glf,allocate_glf",0,istat)

        allocate (zpne_glf(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("zpne_glf,allocate_glf",0,istat)

        allocate (gradni_glf(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("gradni_glf,allocate_glf",0,istat)

        allocate (tau_glf(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("tau_glf,allocate_glf",0,istat)

        allocate (zpni_glf(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("zpni_glf,allocate_glf",0,istat)

        allocate (betae_glf(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("betae_glf,allocate_glf",0,istat)

        allocate (betai_glf(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("betai_glf,allocate_glf",0,istat)

        allocate (beta_glf(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("beta_glf,allocate_glf",0,istat)

        allocate (alpha_exp(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("alpha_exp,allocate_glf",0,istat)

        allocate (rmin_glf(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("rmin_glf,allocate_glf",0,istat)

        allocate (rmaj_glf(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("rmaj_glf,allocate_glf",0,istat)

        allocate (drhogf(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("drhogf,allocate_glf",0,istat)

        allocate (ns_glf(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("ns_glf,allocate_glf",0,istat)

        allocate (bteff_glf(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("bteff_glf,allocate_glf",0,istat)

        allocate (csda_glf(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("csda_glf,allocate_glf",0,istat)

        allocate (rhosda_glf(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("rhosda_glf,allocate_glf",0,istat)

        allocate (drhodr(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("drhodr,allocate_glf",0,istat)

        allocate (drhodrrrho(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("drhodrrrho,allocate_glf",0,istat)

        allocate (angrotp_exp(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("angrotp_exp,allocate_glf",0,istat)

        allocate (egamma_exp(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("egamma_exp,allocate_glf",0,istat)

        if(irotstab .le. 0)then
           allocate (egamma_exp_initial(1:nj),STAT = istat)
           if(istat .ne. 0)
     .        call allocate_error("egamma_exp_initial,allocate_glf",
     .                                                      0,istat)
        endif

        allocate (gamma_p_exp(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("gamma_p_exp,allocate_glf",0,istat)

        allocate (vphim_glf(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("vphim_glf,allocate_glf",0,istat)

        allocate (vparm_glf(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("vparm_glf,allocate_glf",0,istat)

      allocate (vperm_glf(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("vperm_glf,allocate_glf",0,istat)

        allocate (ve_glf(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("ve_glf,allocate_glf",0,istat)

        allocate (fcglf(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("fcglf,allocate_glf",0,istat)

        allocate (ak1glf(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("ak1glf,allocate_glf",0,istat)

        allocate (ak2glf(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("ak2glf,allocate_glf",0,istat)

        allocate (aneoglf(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("aneoglf,allocate_glf",0,istat)

        allocate (diff_m(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("diff_m,allocate_glf",0,istat)

        allocate (chie_m(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("chie_m,allocate_glf",0,istat)

        allocate (chie_etg_m(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("chie_etg_m,allocate_glf",0,istat)

        allocate (chii_m(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("chii_m,allocate_glf",0,istat)

        allocate (etaphi_m(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("etaphi_m,allocate_glf",0,istat)

        allocate (etapar_m(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("etapar_m,allocate_glf",0,istat)

        allocate (etaper_m(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("etaper_m,allocate_glf",0,istat)


        allocate (exch_m(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("exch_m,allocate_glf",0,istat)

        allocate (egamma_m(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("egamma_m,allocate_glf",0,istat)

        allocate (gamma_den(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("gamma_den,allocate_glf",0,istat)

        allocate (gamma_ti(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("gamma_ti,allocate_glf",0,istat)

        allocate (gamma_w(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("gamma_w,allocate_glf",0,istat)

        allocate (vexb_den_grad(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("vexb_den_grad,allocate_glf",0,istat)

        allocate (vexb_ti_grad(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("vexb_ti_grad,allocate_glf",0,istat)

        allocate (vexb_w_grad(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("vexb_w_grad,allocate_glf",0,istat)

        allocate (vexb_tot(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("vexb_tot,allocate_glf",0,istat)

        allocate (gamma_p_m(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("gamma_p__m,allocate_glf",0,istat)

        allocate (anrate_m(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("anrate_m,allocate_glf",0,istat)

        allocate (gamma_net_i(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("gamma_net_i,allocate_glf",0,istat)

        allocate (gamma_net_e(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("gamma_net_e,allocate_glf",0,istat)

        allocate (anrate2_m(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("anrate2_m,allocate_glf",0,istat)

        allocate (anfreq_m(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("anfreq_m,allocate_glf",0,istat)

        allocate (anfreq2_m(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("anfreq2_m,allocate_glf",0,istat)

        allocate (egamma_d(1:nj,pdelay),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("egamma_d,allocate_glf",0,istat)

        allocate (itport_pt(1:itport_glf),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("itport_pt,allocate_glf",0,istat)

        allocate (ve_eff_m(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("ve_eff_mdv ,allocate_glf",
     .                                                      0,istat)

        allocate (vi_eff_m(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("vi_eff_mdv,allocate_glf",
     .                                                    0,istat)

        allocate (vrot_eff_m(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("vrot_eff_mdv,allocate_glf",
     .                                                    0,istat)



        if(dv_method )then

        allocate (diff_mdv(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("diff_mdv,allocate_glf",0,istat)

        allocate (chie_mdv(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("chie_mdv,allocate_glf",0,istat)

        allocate (chii_mdv(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("chii_mdv,allocate_glf",0,istat)

        allocate (etaphi_mdv(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("etaphi_mdv,allocate_glf",0,istat)

        allocate (etapar_mdv(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("etapar_mdv,allocate_glf",0,istat)

        allocate (etaper_mdv(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("etaper_mdv,allocate_glf",0,istat)


        allocate (exch_mdv(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("exch_m,allocate_glf",0,istat)

        allocate (egamma_mdv(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("egamma_mdv,allocate_glf",0,istat)

        allocate (gamma_p_mdv(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("gamma_p__mdv,allocate_glf",0,istat)

        allocate (anrate_mdv(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("anrate_mdv,allocate_glf",0,istat)

        allocate (anrate2_mdv(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("anrate2_mdv,allocate_glf",0,istat)

        allocate (anfreq_mdv(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("anfreq_mdv,allocate_glf",0,istat)

        allocate (anfreq2_mdv(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("anfreq2_mdv,allocate_glf",0,istat)

        allocate (zpti_save(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("zpti_save,allocate_glf",0,istat)

        allocate (chii_save(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("chii_save,allocate_glf",0,istat)

        endif


      return
      end




      subroutine init (tohmwant_dummy, run_preplt_dummy)
c 
c ----------------------------------------------------------------------      @
c --- SUMMARY OF * MAJOR * NEW ONETWO VERSION 3.0 RELEASED FEB 2003            @
c -----------------------------------------------------------HSJ-----------    @
c   !) method of lines
c   2) globally convergent trust region solutions
c   3) glf23
c   4) Sauter bootstrap
c   5) Nclass bootstrap
c   6) tdem for rf codes toray,curray,fastwave
c   7) time dept beam
c   8) nubeam
c   9) partial, ongoing comversion to f90 with modules
c  10) MPI,OMP optimization (doesnt work very well yet)
c  11) dozens of minor  updates/bug fixes, etc
c
c

c ----------------------------------------------------------------------       @
c --- SUMMARY OF * MAJOR * NEW ONETWO VERSION RELEASED LATE JUNE 1997          @
c ----------------------------------------------------------------------       @
c  1) TDEM mode of operations (i.e., time-dependent eqdsk mode)                @
c  2) IFS (Dorland-Kotschenreuther) confinement model                          @
c     (full 2d, 1d model implemented for circular plasmas only !!!!!!!)        @
c  3) New Weiland model (11 equations instead of 6)                            @
c ----------------------------------------------------------------------       @
c  ---- SUMMARY OF MAJOR NEW ONETWO VERSION RELEASED MID-MARCH 1996              @
c ----------------------------------------------------------------------       @
c                                                                              @
c  1) SFAREA in summary page was changed to be the physical area.              @
c  2) All new input variables are described in file "cray102.f" as usual.      @
c  3) An example            case is in file "inone.march"  for shot 87977.     @
c  4) An example d-t fusion case is in file "inone.fusion" for shot 87977.     @
c  5) Use the ONETWO only with the new version of PREPLT and TRPLOT.           @
c  6) The 65x65 version of ONETWO will now process 33x65 EQDSKs as well.       @
c  7) The Weiland confinement model is installed but is not yet tested!!!      @
c  8) Beam-beam, beam-thermal, and thermonuclear d-d and d-t fusion rates      @
c     are calculated with Bosch & Hale cross sections.                         @
c                                                                              @
c ----------------------------------------------------------------------       @
c                                                                              @
c  1) 8/4/99 HSJ fixed w1-4typ input so that time-dependent aspects work       @
c                                                                              @
c ----------------------------------------------------------------------       @
c         
      USE param
      USE aid_newton  
      USE fusion
      USE gcnmp_input,             ONLY : itran_gcnmp,eq_split_gcnmp,    
     .                                    save_incr_restart_file_gcnmp,
     .                                    gcnmp_macro_dt,dtmin_gcnmp ,
     .                                    gcnmp_host,gcnmp_nml_filename,
     .                                    write_iterdb_txt,gcnmp_nprocs,
     .                                    gcnmp_remote_dir,
     .                                    switch_iterdb_output,
     .                                    gcnmp_iterdb_filename
   

      USE fdyimp,                 ONLY :  allocate_fdy_arrays
      USE glf23
      USE io
      USE ions
      USE neut
      USE transp
cJMP  USE nbi_dimensions 
      USE nub  
      USE nub2
      USE ext_prog_info
      USE solcon
      USE solcon_gcnmp,              ONLY : ts_smfactor,
     .                                      single_density_simulation,
     .                                      use_stab_flux,dbar_stab
      USE soln
      USE contour
      USE limiter
      USE mhdpar
      USE mhdgrid
      USE rf
      USE tdem
      USE ename
      USE extra
      USE xptor_sim ,only : test_xptor
      USE yoka
      USE numbrs
      USE toq_12_interface
      USE mesh
      USE verbose
      USE adaptive
      USE sourc
      USE machin
      USE nub4
      USE tfact
      USE geom
      USE events

      USE flags
      USE tordlrot
      USE constnts
      USE nub3
      USE soln2d
      USE tcoef
      USE bd_condtn
      USE mixcom
      USE kinetic_efit,only : wrt_kinetic_efit
      USE rhog
      USE mhdcom
      USE paleocl,            only : include_paleo,paleo_mult,
     .                               paleo_pinch_mult,r_paleo_pinch_cut,
     .                               paleo_pinch_mult2

      USE grid_class,         ONLY : nj_ncd => nj,njm1_ncd=>njm1,
     .                               npsi_ncd => npsi

      USE echdat_module,      ONLY : echdat_read,nsys,set_irf
      USE nonlin
      USE ifs
      USE mcgo
      USE staebler
      USE shapctr
      USE flxav
      USE iterdbmd
      USE cer
      USE tmpcom
      USE mmm95_mod !JMP
      USE brainfuse_mod         !OM
      !USE platform !jmp.ibm.par
      USE rad_loss,                    ONLY : brems_nions          
      USE gpsi
      USE weiland
      USE island
      USE zeffcom
      USE zen_state,                    ONLY : set_dzdt
      USE mhdbcdtn
      USE pelcom

      USE P_Nfreya_12_interface,        ONLY : use_P_Nfreya,          
     .                                         P_Nfreya_read,
     .                                         P_Nf_events,sent_ufile
      USE P_Nfreya_rpc_interface,       ONLY : setup_run_P_Nfreya
      USE Nfreya_namelist,              ONLY : 
     .                                   write_P_Nfreya_run_directives,
     .                                   P_Nfreya_dt

       USe replace_imsl,                ONLY : my_usmnmx,my_vsrta

       USE file_proc,                   ONLY : delete_file

      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c 1) initializes input variables to their default values and reads
c    the overriding input file normally named "inone"
c 2) edits the namelist onto all required output files
c 3) does initialization for job set-up.
c ----------------------------------------------------------------------
c

D      include 'mpif.h'               !required for MPI

      include 'co2.i'
      include 'gauss_info.i'
      include 'glf.m'

      include 'imsl.i'


      include 'mod_gbohm.i'
      include 'quit.i'
      include 'sxrcom.i'

      include 'storage.i'
      include 'rebut.i'    ! to get wrebut into common block

      include 'pckcom.i'


      include 'spare.i'

      include 'fixbdry.i'
      include 'shay.i'       ! 11/11/93
      include 'fitparms.i'

      include 'shape.i'


      include 'wedgerf.i' 
c
      integer  sizel,userid(2),beam_restart_file_length,
     .         taskl,
     .         GETUID
      external RANDOM12,             ! portable random number generator
     .         GETUID                ! get user's identification number
      logical fcd_path
                                             
      integer  strleng         ! length of onetwo (below)
c
      logical      run_preplt, run_preplt_dummy,file_exists
      dimension    itenp(kprim),iteni(kimp)
      dimension    profin(kj)
      dimension    rfzone(10),nprf(10)
c      dimension    eni(kj,kimp),enp(kj,kprim)
      character*7  cursrce, btorsrce

c
      character*8  iflag, title(10)
      data         iflag /'--------'/

      character *16 task
    

c ----------------------------------------------------------------------       @
c --- GEOMETRY FLAG AND MACHINE NAME                                           @
c ----------------------------------------------------------------------       @
c  The following flags appear in the FIRST LINE of the input file "inone",     @
c  in the order given. Free-field format decoding is done to get the names.    @
c                                                                              @
c% codeid        Geometry identification flag                                  @
c    DEFAULT -> 'onedee'   : 1-D     run; circular or elliptical flux surfaces @
c                  'dee'   : 1-1/2-D run;             dee-shaped flux surfaces @
c                                                                              @
c     runid      characer string (129 or less) identifying the run
c         DEFAULT -> 'Onetwo_run'
c% machinei      Machine identification flag                                   @
c               'doub-iii' : old Doublet III machine                           @
c    DEFAULT -> 'diii-d'   : current  DIII-D machine                           @
c
c    functions qq and qqq (in cray401.f) were renamed  fqq and fqqq            @
c    respectively. This was done to eliminate a compiler warning regarding     @
c    the use of qq as both a functio and as an array.                          @
c    Appropriate changes were made in cray306.f (subroutine DELSOL)            @
c                                                                              @
c ----------------------------------------------------------------------       @
c   


c ----------------------------------------------------------------------       @
c --- DEFINE THE THREE NAMELISTS (NAMELIS1, NAMELIS2, NAMELIS3)                @
c ----------------------------------------------------------------------       @
c                                                                              @
      NAMELIST /namelis1/                                                      %
     .  rmajor, rminor, elong, btor, nprim, nimp, namep, namei, fd,            %
     .  itenp, iteni, itte, itti, itxj, taupin, icenez, inenez, zfrac,         %
     .  adjzeff, w1saw, w2saw, w3saw, wneo, wneot, w3cla, w1typ, w12typ,       %
     .  w13typ, w2typ, w3typ, w4typ,w1typmin,w2typmin,w3typmin,w4typmin,       % 
     .  typa, wstd, w1isl, w2isl, w3isl,w2cgm,w3cgm,w2cg,w3cg,                 %
     .  nisl, w1mix, w2mix, bparzeff, w3mix, w4mix, w5mix, qmix, rsmixx,       %
     .  tsmix, tdmix, ipmix, w0mix, dtemix, dtimix, fusmix, mixpro,            %
     .  s3mix, s71mix, s18mix, trmix, jsxr, jco2, jzeff, jterow, nterow,       %
     .  imesh, nj, rnormin, njenp, njeni, njte, njti, njene, njzef,            %
     .  njcur, r, dr, rin, rout, zax, ene, enec, eneb, alpene, gamene,         %
     .  enp, enpb, enpc, alpenp, gamenp, eni, enib, enic, alpeni,              %
     .  gameni, te, teb, tec, alpte, gamte, ti, tib, tic, alpti, gamti,        %
     .  curden, xjb, xjc, alpxj, gamxj, zeff, zeffb, zeffc, alpzef,            %
     .  gamzef, xmtmdifs, rmtmdifs, bparenp, bparte, nbctim, bctime,           %
     .  totcur, iffar, enein, tein, tiin, zeffin, bpareni, bparti,             %
     .  time0, timmax, nmax, dt, dtmin, dtevmin, dtmax, relmin, relmax,        %
     .  relit, itmax, theta, timav, timprt, timplt, mprt, mplot, prtlst,       %
     .  pltlst, curdenin, jflux, jcoef, jsourc, jbal, jtfus, nneu,             %
     .  namen, iref, ipcons, gasflx, recyc, twall, wion, wrad, raneut,         %
     .  relneu, erneumax, idiagn, nengn, englstn, jprt, rtandn, rhdn,          %
     .  widths, nps, enes, tes, ttweak, relaxtyp,diffeq_methd,                 %
     .  nw1pro, w1pro, nw2pro, w2pro, nw3pro, w3pro, nw4pro, w4pro,            %
     .  bl2in, fusnin, bparcur, ticin, voltin, qcin,r_clamp,chi_clamp,         %
     .  zeflim, bparang, w33min, xdebug, ifred, ilimdt, vtyp, nvpro,           %
     .  vpro, ibaloo, voltav, betlim, betlim0, iangrot, itangrot,              %
     .  angrotin, rangrot, iwangrot, bparene, etaioff, renpin, reniin,         %
     .  rtein, rtiin,renein, rzeffin, rcurdein,splninpt, jgboot,               %
     .  jhirsh, jneo, implicit_fh, wrebut, relaxsol, relaxrebut, ddebug,       %
     .  resistive, ftcalc, qrebsmth, wshay, smult, scsmult, snexp,             %
     .  sbpexp, srexp, sbigrexp, stexp, sdtdrexp, srin, srout, suserho,        %
     .  skimult, srincev, sroutcev, relaxshay, tohmwant, ifsflag, aeh,         %
     .  ael, aih, ail, bh, bl, ch, cl, alfae, alfai, betah, sigma,             %
     .  gammah, lsfctr, tirlw, rlw_model, run_preplt, xrot, xeden,             %
     .  xsecder, ishayform, skimass, times_rgc, timee_rgc, irgc, rgc,          %
     .  rgca, rgcb, rgcc, rgcd, rgce_mult, rgci_mult, rgc_string,              %
     .  fusionvb, tportvb, namelistvb, zenvb, exptl_neutron_rate,              %
     .  neucgvb, wweiland, include_weiland, cer_ion, tdemvb, fiziksvb,         %
     .  squeeze, ikpol, kpolin, rkpol, bparkpol,wneo_elct_ifs,                 %
     .  include_ifs, dorl_kotch, eps_adaptive, gridgenvb,freeze_type,          %
     .  speed_adaptive, freeze_adaptive, curve_eps_adaptive,freeze_xsn,        %
     .  spec_profiles, no_te_convection, no_ti_convection, rho_edge,           %
     .  timecrit, dtmaxcrit, fix_edge_te, fix_edge_ti, extend_seval,           %
     .  set_chie_chii,testing_NTCC,ce0_mgb,alpha_mgb,ci1_mgb,neomult,          %
     .  ci2_mgb,ce_bgb,cfe_mgb,cfe_bgb,cfi_mgb,cfi_bgbc1_g,c2_g,c_theta,       %
     .  include_itb, include_glf, dorl_kotche,dorl_kotchi,exbmult_ifs,         %
     .  exbmult_glf,limp_glf,ivphi_glf,iglf_chie,iglf_chii,iglf_chiv,          %
     .  ibtflag_glf,lprint_glf,steps_per_plot, momtm_file, do_eqplot,          %
     .  neomult1,neomult2,qneomult,dcparam_in,continuation_method,             %
     .  steady_state,iteration_method,rbsaxis,non_lin_method,ssqrmin,          %
     .  steptol,gradtol,fvectol,bandwidth,maxfev,fdigits,switch_method,        %
     .  tot_iters_max,jmm_glf23,zeff_mult,jacobian_type , fix_edge_rot,
     .  cb_mgbr,qn_mgbr,lqn_mgbr,time_mgbr,berrqn_mgbr,iters_freeze,
     .  glf23_ncpus, glf23_iglf,jroot_iglf, jeigen_iglf,cximult,
     .  irotstab,i_delay,x_alpha_glf,glf_debug,jac_skip,jelc_clamp,
     .  random_pert,freeze_alpha_exp,conv_skip,te_range_check,
     .  q0_max,q0_radius,q0_mult, runid, write_profiles,analysis_check,
     .  write_glf_namelist,newtonvb,use_avg_chi,test_xptor,
     .  itte_dv,itti_dv,itenp_dv,itene_dv,itangrot_dv,dv_delt,
     .  jion_clamp,
     .  nalp_thresh,nb_thresh, dn0out,nlfbmflr,bp0_ic,r_elm,t_elms,
     .  t_elme,etam_elm,itot_elm,vloop_bc,vloop_bc_time,vloopvb,
     .  use_pedestal,pedestal_path,pedestal_models,include_paleo, 
     .  paleo_mult,paleo_pinch_mult,wrt_nwt,r_paleo_pinch_cut,dbar_stab,
     .  paleo_pinch_mult2,save_scratch1,mult_curboot,use_stab_flux,
     .  include_mmm,immm_chie,immm_chii,immm_chiv,exbmult_mmm,  !JMP
     .  fix_edge_ni,iperp,iglf_eq,iglf_idt, !jmp.den
     .  chid_glf_d,chid_glf_e,ichid_smooth,ichiv_model,ichiv_chii, !jmp.den jmp.snu
     .  host_name_ext                                              !jmp.ibm.par
     .  itran_gcnmp,save_incr_restart_file_gcnmp,gcnmp_macro_dt,
     .  gcnmp_host,gcnmp_nml_filename,write_iterdb_txt,gcnmp_nprocs,
     .  switch_iterdb_output,gcnmp_remote_dir,gcnmp_iterdb_filename,
     .  en_bc_inpt,ren_bc_inpt,iglf_d,fix_edge_ni_bc_inpt,ts_smfactor, 
     .  single_density_simulation,density_mult,set_te_to_ti,
     .  set_ti_to_te,te_mult,ti_mult,ang_mult,ene_mult,
     .  include_ntv,c_p,delta_b_sqr,mp_polnum,mp_tornum , !C.K. Pan
     .  he3_frac,he3_thermal_spin_pol,
     .  include_brainfuse, brainfuse_path
     
c                                                                          
     %
      NAMELIST /namelis2/  beam_thermal_fusion,                                            
     .  timbplt, beamon, btime, nameb, relnub, tfusbb, anglev, angleh,         %
     .  nashape, aheigh, awidth, bcur, bptor, blenp, nbshape, bleni,           %
     .  bheigh, bwidth, bhfoc, bvfoc, bhdiv, bvdiv, ebkev, fbcur,              %
     .  nbeams, naptr, alen, bvofset, bhofset, nsourc, sfrac1, mf,             %
     .  npart, npskip, rpivot, zpivot, ranseed, fionx, iddcal, fdbeam ,        %
     .  nbinject, rfmode, rfon, rftime, turnonp, turnonc, rfpow, idamp,        %
     .  xec, zec, irfplt, gafsep, freq, rfrad1, rfrad2, wrfo, wrfx,            %
     .  rnormin, njqin, rfrad1ic, rfrad2ic, necsmth, qine, qini, xdebug,       %
     .  a1rf, a2rf, wrfe, wrfi, nrfzon, rfzone, ichmod, betalm, relrf,         %
     .  nprf, iside, xkpar, nhigh, ykperp, navg, thetec, phaiec, hlwec,        %
     .  nampel, pelrad, vpel, nbgpel, timpel, pelmod,                          %
     .  inubpat, npat, ratwec,ech_input,relrf_pow,                             %
     .  nray, ifus, iaslow, wtifus, ibcur, ibcx, ibslow, iborb, iyoka,         %
     .  ishot, itime, itrapech, itrapfi, wdelt, wgam, wohm, nqrad,             %
     .  qradr, qradin, refrad, ds_tk, fe_tk, ne_tk, iexcit, ilorent,           %
     .  mstate, ncont, izstrp, kdene, nbeamtcx, kdeni, kdenz, ksvi,            %
     .  ksvz, ksve, krad, ngh, ngl, nouthx, rnp, irfcur, ifb, rfcur,           %
     .  lmode, ifbprof, alphaf, hdepsmth, zrffw, nzrffw, pzrffw, htsfw,        %
     .  iswchfw, lifw, nihfw, timrfp, freqfw, rnpfw, impath, angrm2d,          %
     .  angrcple, nrayptrt, powersrt, nnkparrt, heightrt, maxrefrt,            %
     .  islofart, anzinfrt, anzsuprt, nthinrt, nfwsmth,rframp_timedown,        %
     .  nkfcd, xntor, rpant, gamloss, iterate_beam, rframp_timeup,             %
     .  iddfusrate, iddfusb, iddfusb_s, ddfusb_t, iddfusb_bulk,                %
     .  icalc_cxfactor, imaxbeamconvg, beam_beam_fusion,                       %
     .  extcurrf, extqerf, extqirf, extcurrf_id, extqerf_id, extqirf_id,       %
     .  extcurrf_amps, extqerf_watts, extqirf_watts, extcurrf_rho,             %
     .  extqerf_rho, extqirf_rho, extcurrf_curr, extqerf_qe, extqirf_qi,       %
     .  extqerf_nj, extqirf_nj, extcurrf_nj, freyavb,nicsmth,                  %
     .   fast_ion_target, rtstcx, relaxden,curray_fi,                          %
     .  relaxden_err,npart_mcgo,mcgo_output_file, use_Callen,                  %
     .  neg_ion_source, time_dep_beam,beam_thermal_cutoff,                     %
     .  beam_mode,beam_restart_file,source2_phase, external_beam_cur,          %
     .  beamoff,ngauss,beam_init_restart_file,beam_end,                        %
     .  beam_thermal_speed,beam_pulse_control,beam_cycles,                     %
     .  knotscurb_external,rcurb_external,beam_data_namelist,                  %
     .  bparcurb_external,curbeam_external,beam_data_ufile,                    %
     .  use_nubeam, toray_version,toray_path,                                  %
     .  echin_save ,fastcd_path,nnkpolrt,anpinfrt,anpsuprt,                    %
     .  save_curray_input,psistep,pkexpnt,epserr,epser1,igraph,                %
     .  iprint,idcur,indvar,ichois,ichoisrt,igrill,modcd,idmpsw,nminor,        %
     .  kalfa,idmpswrt,use_P_Nfreya,                                            %
     .  nspect,bmaxrt,curray_path,incrt,icurdrrt,thgrilrt,psi_startrt,         %
     .  irayiort,beam_spec, dtn_orbit,nptcls,nptclf,nubeam_restart,            %
     .  wrt_restart_file_time,nubeam_state_path,nubeam_xplasma_path,           %
     .  nubeam_path,save_nubeam_input,nubeam_setup_ext,tmin_curray,            %
     .  enbmin_curray,wrt_kinetic_efit,genray_path,save_genray_io,             %
     .  genray_fi,genraydat,avg_nubeam_torque,fidiff_on,bfr_neutrlz,
     .  nubeam_back_delt,nubeam_fix_t,ifix_nubeam_dt,use_ufile, !JMP
     .  nubeam_back_average,P_Nfreya_dt,nubeam_version,    !JMP
     .  d_beam_spin_pol,he3_thermal_spin_pol


c                                                                      
      NAMELIST /namelis3/                                                      %
     .  ifixshap, mhdmode, xdim, ydim, redge, nlimiter, xlimiter,              %
     .  ylimiter, ifill, greentab, isym, mhdmethd, fixfcoil, ifitpsi,          %
     .  ifitprob, iecurr, ivessel, timeqbcd, flxeqbcd, curmax, toleq,          %
     .  tolcur, iteq, itcur, omeq, omcur, npsi, ieqmax, deqlst, delcap,        %
     .  delrho, minchisq, maxitr, ibypas, dtmine, ecurrt, fcoilcur,            %
     .  vesscur, irguess, eqdskin, pcurmhd, btormhd, vloopmhd, expmp2,         %
     .  volaray, voladj, volnudge, itvol, tolvol, limpos, dvadjmax,            %
     .  ieqprt, j2prt, iprtit, ieqdsk, rf_output, ispare, aspare,              %
     .  xdebug, ddebug, rmagax, zmagax, rscale, zscale, mhdonly,               %
     .  fitfcur, errpsilp, errmprbe, errfcoil, errecoil, errvescr,             %
     .  ivertsbl, rvloop, zvloop, ltest_code, optwi, optwe, optomegai,         %
     .  optomegae, tstark, sigstark, fwtstark, rstark, zstark, a1stark,        %
     .  a2stark, a3stark, a4stark, use_stark, timestark, xi_include,           %
     .  psifctr, tensionspl, mresidual, nfitpoints, splinefitbd,               %
     .  fitboundary, derwght, adotmult, f2d1mult, f2d2mult, f2d3mult,          %
     .  use_cnt1, use_cnt2, iterdb, iterdsc, itorfluse, irwflag,               %
     .  comp_methd_eqdsk, use_efit_cntr, curtype, pol_flux_lim,use_Bp0,        %
     .  cap_mult, deltat_fixed_boundary, qe2dmult, qi2dmult,dlnhdtmult,        %
     .  dlnhdrmult, dlnrdtmult, renormalize_rbp,toq_path,ieqdtoq,              %         
     .  npsi_toq,nthet_toq,fixcur, baxis0,imislo,alpsi,modelbnd,ishape,        %
     .  eshape,xshape,equiltype,modelf,nffp,ffpcoef,modelp,modelq,             %
     .  nppcoef,betafix,betiso,ieqdsk_toq, prtbal,prtboot,prtgato,             %
     .  prteqdata, prtfixb,prtflux,bavmg,minmg, maxmg,tolmg,loopinmg,          %
     .  nhior, nbndry, nbug, updownsym,psifactr_toq,dsrat, igimblett,          %
     .  ppcoef,zeff_toq,qedge_toq, ieqdsk_toqp,iyoka,ishot,itime,              %
     .  wrt_kinetic_efit ,create_GCNMP_input,smth_mhd_parms,                   %
     .  initialize_from_statefile,statefile_name,write_iterdb_txt,             %
     .  create_XPTOR_input ,switch_iterdb_output
*                                                                              @
c **********************************************************************       @
c *** FIRST NAMELIST (NAMELIS1) ****************************************       @
c **********************************************************************       @
c                                                                              @

        
c ----------------------------------------------------------------------       @
c --- MACHINE PARAMETERS                                                       @
c ----------------------------------------------------------------------       @
c% rmajor    Major radius (cm) to magnetic axis in 1-D runs.                   @
c            In 1-1/2-D runs this radius is used only to specify btor.         @
c            Moreover, if irguess < 0 and rmajor = 0, then rmajor is           @
c            read from the eqdsk file.                                         @
c% rminor    Horizontal minor radius (cm) in 1-D runs (recalculated in         @
c            1-1/2-D runs)                                                     @
c% elong(m)  Elongation (height-to-width ratio) of elliptical flux             @
c            surfaces in 1-D runs. If nbctim > 1, elong(2) .ne. 0              @
c            specifies that this quantity is time-dependent (rminor            @
c            is held fixed).  See below on nbctim and bctime(m).               @
c% rin       Inside major radius (cm) of vacuum vessel used to determine       @
c            whether neutral beams or SXR chords are reentrant . Note that     @
c            this value is also used in subroutine inject to determine the     @
c            injection geometry. For 2D it must be set to the inside radius    @
c            of the mhd grid.                                                  @
c% rout      Outside major radius (cm) of vacuum vessel used only by           @
c            neutral beam plotting code NUBPLT                                 @
c% zax       Vertical position (cm) of magnetic axis in 1-D runs               @
c            (recalculated in 1-1/2-D runs)                                     @
c% btor      Toroidal magnetic field (Gauss) at rmajor.  In 1-1/2-D            @
c            runs if irguess < 0 and btor = 0.0, then btor is read from        @
c            the eqdsk file. if btor = -1.0e30 on input then btor is           @
c            obtained from the third namelist of inone from the vector         @
c            btormhd, in a fashion similar to the treatment of totcur          @
c            (see totcur description). (btor is nominally time-independent     @
c            for normal discharges. We allow time-dependent input              @
c            for generality using btormhd).                                    @
c% tportvb   integer: if set to 1, gives some detailed output to screen        @
c            (primarily for diagnostic use)                                    @
c% namelistvb (an integer):                                                    @
c%            If set to 1, it will print a message to the                      @
c             screen as each namelist is read in. Use this to isolate          @
c             undefined variables in the namelists.                            @
c             It also causes the three namelists to be dumped consecutively    @
c             to a file named "namelists".                                     @
c% zenvb      an integer if set to 1 will print first few elements of          @
c             subroutine ZEN (i.e., density) calculations to screen            @
c% tdemvb     an integer if set to 1 will cause dump_data to be called         @
c             (used for debugging)                                             @
c% fiziksvb   an integer if set to 1 will cause FIZIKS to dump some info       @
c             to be called (used for debugging)                                @


      host_name_ext = 'undefined' !jmp.ibm.par

c ----------------------------------------------------------------------       
c --- get  location of cross section data files
c ----------------------------------------------------------------------
c
cjmp.ibm.par       call get_xsect_path(ncrt,nout) ! will be called 
                                                  ! after reading namelist

      point99   = 0.99
      point999999 = 0.999999
c


c
c
c ======================================================================
c
      if (machinei .eq. 'doub-iii') then
        rmajor   = 143.0
        rminor   =  45.0
        elong(1) =   1.0
        rin      =  96.0
        rout     = 192.0
        zax      =  89.0
        btor     =  24.0e3
      else
c
c       DIII-D parameters
c
        rmajor   = 169.5
        rminor   =  66.5
        elong(1) =   1.0
        rin      = 101.39      ! old value
        rin      = 105.313135  ! this must be .ge. inside edge of eqdsk
c                                for 2D cases. It is used in inject
c                                to set up critical injection parameters ... HSJ
        rout     = 237.56   ! note that it is rmax, not rout, that goes with rin
        zax      =   0.0
        btor     =   0.0    !changed from 24e3 06/11/02 HSJ
      end if
*                                                                              @
c ----------------------------------------------------------------------       @
c --- ION PARAMETERS                                                           @
c ----------------------------------------------------------------------       @
c% nprim       Number of      primary ion species (1 to 3)                     @
c% nimp        Number of     impurity ion species (0 to 2)                     @
c% namep(i)    Name   of ith  primary ion species                              @
c              'h' , protons                                                   @
c              'd' , deuterons                                                 @
c              't' , tritons                                                   @
c              'dt', mixture of d and t                                        @
c              'he', thermal alphas                                            @
c% namei(i)    Name   of ith impurity ion species                              @
c              'he', helium                                                    @
c              'c' , carbon                                                    @
c              'o' , oxygen                                                    @
c              'si', silicon                                                   @
c              'ar', argon                                                     @
c              'ti', titanium                                                  @
c              'cr', chromium                                                  @
c              'fe', iron                                                      @
c              'ni', nickel                                                    @
c              'kr', krypton                                                   @
c              'mo', molybdenum                                                @
c              'w' , tungsten                                                  @
c% fd          Number fraction of deuterons in d-t mixture                     @
c ----------------------------------------------------------------------       @
c                                                                              @
      do i=1,kprim                                                             %
        namep(i) = ' '                                                         %
      end do                                                                   %
      do i=1,kimp                                                              %
        namei(i) = ' '                                                         %
      end do                                                                   %
      nprim      = 1                                                           %
      nimp       = 0                                                           %
      namep(1)   = 'h'                                                         %
      fd         = 0.5                                                         %
      tportvb    = 0    ! controls diagnostic output                           %
      fiziksvb   = 0                                                           %
      namelistvb = 0                                                           %
      zenvb      = 0                                                           %
      neucgvb    = 0                                                           %
      tdemvb     = 0                                                           %
*                                                                              @
c ----------------------------------------------------------------------       @
c --- TRANSPORT FLAGS                                                          @
c ----------------------------------------------------------------------       @
c  analysis_check if = 1   (which is the default) the
c                 profiles run in analysis mode are included in the
c                 check for the relative maximum change, see relmax.
c                 This means that the time step might be cut back
c                 even if we are not solving the corresponding
c                 diffusion equation. This is useful for checking the
c                 consistency of the time dependent profiles that are 
c                 specified in inone. If this effect is not wanted
c                 then set analysis_check = 0
c                 
c% testing_NTCC  special flag used to compare solutios with NTCC code          @
c                (intended for use by developers only)                         @
c% itenp(i)   1, transport of primary ion species i is considered,             @
c                i.e., profile is calculated                                   @
c             0, transport of primary ion species i is not considered,         @
c                i.e., profile is specified                                    @
c% itte       Same as above for electron temperature                           @
c% no_te_convection  Set to 1 to turn OFF convective ELECTRON energy flux      @
c                    set to -1 to use 3/2 instead of 5/2 for convection        @
c                                                                              @
c% itti       Same as above for ion temperature                                @
c 
c set_te_to_ti
c set_ti_to_te   integer switches 0 (default) means no effect.
c                =1 means one of te,ti is set to the other.
c                (only one should be set to 1 and the would presumably
c                 be evolved,but can also be used to match boundary condition
c                 input profiles )
c                useful when evolving one profile and want te other to follow.
c% no_ti_convection  Set to 1 to turn OFF convective      ION energy flux      @
c                    set to -1 to use 3/2 instead of 5/2 for convection        @
c% itxj       Same as above for current density                                @
c% iangrot    switch which includes ( = 1) or excludes (=0) the toroidal       @
c             angular momentum.  If iangrot = 0 (default), the code runs       @
c             in its old mode of operation (i.e., toroidal angular             @
c             momentum effects are totally absent).  if iangrot = 1, then      @
c             toroidal angular momentum effects are to be included,            @
c             depending on the setting of switch itangrot.  (the purpose       @
c             of switch iangrot is to allow automatic setting of               @
c             compatibility with older versions of the code which do not       @
c             include toroidal rotation.  As a consequence, two switches       @
c             are required to include toroidal rotation).                      @
c  no_w_convection        Set to 1 to turn OFF convective  momentum  flux      @
c                                                                              @
c momtm_file  if iangrot =1 then a detailed description of the toroidal        @
c             rotation is given in file momtmout provided that momtm_file =1   @
c             default is 0                                                     @
c itangrot   switch which is used only if iangrot = 1. If iangrot=1 and       @
c             itangrot = 1 the toroidal angular momentum equation is           @
c             included in the set of transport equations to be solved.         @
c             (i.e., the simulation mode).                                     @
c             if iangrot = 1 and itangrot=0 the toroidal angular momentum      @
c             is not transported (the initial angular rotation speed           @
c             profile is independent of time,or input as a function of         @
c             time, (i.e., the analysis mode).                                 @
c             The source terms due to the presence of a nonzero rotation       @
c             speed are incorporated in the remaining equations (the           @
c             sources can be selectively turned off in both analysis and       @
c             simulation runs-see below).Only first order corrections          @
c             in rotation speed / ion thermal spee are included.               @
c             default is itangrot = 1.0                                        @
c% ikpol      If ikpol = 1, spline data for kpolin is expected                 @
c             (poloidal velocity of CER ion/Bp)  (default = 0)                 @
c% nbeamtcx   switch related to torque calculation.  should be set if          @
c             angular rotation is used.  see description under beam            @
c             input (note: nbeamtcx is input in the second namelist).          @
c ----------------------------------------------------------------------       @
c                                                                              @
      do index=1,kprim                                                         %
        itenp(index) = 1                                                       %
      end do                                                                   %
      analysis_check = 1
      ikpol =0
c                                                                              %
      testing_NTCC     = 0 ! no tests are performed                            %
      itte             = 1                                                     %
      itti             = 1                                                     %
      itxj             = 1                                                     %
      iangrot          = 0                                                     %
      itangrot         = 1                                                     %
      implicit_fh      = .false.                                               %
      momtm_file       = 0                                                     %
      no_te_convection = 0                                                     %
      no_ti_convection = 0                                                     %
      no_w_convection  = 0
      set_te_to_ti     = 0
      set_ti_to_te     = 0
*                                                                              @
c ----------------------------------------------------------------------       @
c --- ANALYSIS MODE PARAMETERS                                                 @
c ----------------------------------------------------------------------       @
c                                                                              @
c% taupin     Particle confinement time (s) for each primary species           @
c             with itenp(i) = 0 (same as global electron confinement time)     @
c% ttweak     Time step (s) between adjustments (tweaks) to obtain             @
c             a specified value for one of several parameters                  @
c% fusnin     Specified fusion neutron rate (s-1) obtained by                  @
c             adjusting wneo(3,3) or w3typ (w3typ>0 takes precedence           @
c             over wneo(3,3); fusnin>0 takes precedence over ticin);           @
c             see also the neutral beam parameters iddcal and fdbeam in        @
c             the second NAMELIST                                              @
c% ticin      Specified central ion temperature (keV) obtained by              @
c             adjusting wneo(3,3) or w3typ (w3typ>0 takes precedence           @
c             over wneo(3,3))                                                  @
c% voltin     Specified one-turn voltage (V) obtained by adjusting             @
c             zeffc and zeffb                                                  @
c% voltav     Specified average voltage (V) across plasma obtained by          @
c             adjusting zeffc and zeffb (voltav>0 takes precedence             @
c             over voltin)                                                     @
c% qcin       Specified q on axis obtained by adjusting zeffc/zeffb            @
c% w33min     Lower limit allowed for wneo(3,3) or w3typ                       @
c% zeflim     Upper limit allowed for zeff (if non-zero)                       @
c% tohmwant   Total desired ohmic current in steady state (can be negative)    @
c               tohmwant is used only in subroutine SSCURDRV to determine      @
c               required driven current in steady state                        @
c ----------------------------------------------------------------------       @
c                                                                              @
c      tohmwant = -1.0e+100    ! turns off call to SSCURDRV                     %
      tohmwant = -1.0e+30    ! turns off call to SSCURDRV  for lf95            %
      ttweak   =  0.05                                                         %
      taupin   =  1.0                                                          %
      fusnin   =  0.0                                                          %
      ticin    =  0.0                                                          %
      voltav   =  0.0                                                          %
      voltin   =  0.0                                                          %
      qcin     =  0.0                                                          %
      w33min   =  1.0                                                          %
      zeflim   =  0.0                                                          %
*                                                                              @
c ----------------------------------------------------------------------       @
c --- FLAG FOR TOGGLING THE RUNNING OF PREPLT                                  @
c ----------------------------------------------------------------------       @
c                                                                              @
c%    run_preplt = .true.  => DO     call PREPLT postprocessor right before    @
c                             the end of ONETWO's execution, so that TRPLOT    @
c                             can be run subsequently. This is the DEFAULT.    @
c                                                                              @
c     run_preplt = .false. => DO NOT call PREPLT postprocessor right before    @
c                             the end of ONETWO's execution.                   @
c                             If you later decide that you want to run TRPLOT, @
c                             remember to first run (standalone) PREPLT, since @
c                             it would not have already been called by ONETWO. @
c                                                                              @
      run_preplt = .true.   ! PREPLT will be run at the end of the ONETWO run  %
*                                                                              @
c ----------------------------------------------------------------------       @
c --- PARAMETERS FOR BASIC TRANSPORT MODELS                                    @
c ----------------------------------------------------------------------       @
c% wneo                                                                        @
c  ((wneo(i,k), i=1,5),     Weights for neoclassical transport                 @
c               k=1,5)       if iangrot = 0 (i.e., toroidal momentum effects   @
c                              are absent), wneo is input exactly as           @
c                              before.  the increase in dimension from         @
c                              4 to 5 is accounted for in the code.            @
c                            if iangrot = 1 wneo must be input as a 5x5        @
c                              array!  This means adding a trailing            @
c                              number to each row of the old 4 row wneo        @
c                              array and adding a 5th row of 5 elements.       @
c                              Note that this scheme will allow input of       @
c                              indicies to the wneo array by assigning         @
c                              each input value only if iangrot = 1.  if       @
c                              iangrot = 0, statements of the type             @
c                              wneo(3,2) = 5, for example, will not work!      @
c                              for iangrot = 0, the entire 4x4 wneo array      @
c                              must be input.  indicies assigned to the        @
c                              elements of the array must be obtained          @
c                              through the standard Fortran default.           @
c                                                                              @
c  jneo        Selection Switch for neoclassical transport                     @
c              = 0  (Default) Rawls, Chu, and Hinton model                     @
c              = 1  NCLASS model   (not operational)                           @
c              = 2  NOT USED
c                         
!             KAPISN (neoclassical)  Models:
c  jneo        = 3  ==>   nkimod = 1 Rutherford neoclassical ion conductivity - Julich
c              = 4        nkimod = 2 Modified Hazeltine-Hinton as used in BALDUR code
c              = 5        nkimod = 3 Bolton model implemented by McCune 7/81
c              = 6        nkimod = 4 Chang-Hinton model, Phys Fluids 25, 1982
c              = 7        nkimod = 5 Modified Chang-Hinton model for Zeff > 1, ppf 1986
c                                                     @
c% jhirsh      Selection Switch for Bootstrap Current Model                    @
c  jboot       jboot is a pseudonym for jhirsh
c              = 0  (Default) The calculation is based on Hinton & Hazeltine   @
c                      (Rev. Mod. Phys. 48 (1976) 239).                        @
c              = 1  The bootstrap current is modeled using the formulae of     @
c                      Hirshman (Phys. Fluids 21 (1978) 1295).                 @
c              = 2  Similar to jhirsh = 1 but includes the fast ions in        @
c                      with the thermal ions.                                  @
c              NOTE:   jhirsh = 1,2 may modify the bootstrap                   @
c                      current calculations near the magnetic                  @
c                      axis by the poloidal field gyroradius.                  @
c                      if this modification is not wanted use                  @
c                      jhirsh = 11 instead of jhirsh = 1 and                   @
c                      jhirsh = 22 instead of jhirsh = 2 ... HSJ               @
c              = 88 The small collisionality, arbitrary aspect ratio,          @
c                      single ion fluid model of Hirshman 1988.                @
c                      (Phys. Fluids, vol 31 (10) 1988)                        @
c                      No corrections for the behavior near the magnetic       @
c                      axis due to finite poloidal gyroradius are done.        @
c              = 89 Similar to jhirsh = 88 but includes the fast-ion density   @
c                      gradient, assumed to be part of the single-ion fluid.   @
c              = 95 The arbitrary collisionality, arbitrary aspect ratio,      @
c                      and multiple ion species model of Houlberg, 1995        @
c                      (Nuclear Fusion, to be published).                      @
c                      Considered the most comprehensive model.                @
c                      No provision has been made to account for the fast ion  @
c                      bootstrap (Electrons from fast ions and alphas are      @
c                      included.)                                              @
c              NOTE: Houlberg bootstrap implementation has  max charge of      @
c                    impurity = 8 so limited to oxygen
c                    (see houlberg.i,parameter mxmz)
c                    attempt to run with higher z is not checked for ! HSJ 1/16/13
c              = 96 Houlberg model, including a fast ion pressure term for     @
c                      the beam and alpha particles. The beam and alpha        @
c                      pressures are derived from the averaged stored energy   @
c                      density (keV/cm3) of the beam and alphas, i.e.:         @
c                                 press = nkT = 0.666667*W                     @
c                      See subroutine NCLBOOT for more information.            @
c             NOTE: Both jhirsh = 95 and 96 use the formula of Y. R. Lin-Liu   @
c                      R. L. Miller for the trapped fraction calculation       @
c                      unless  ftcalc = 'exact'  is set.                       @
c             99,100 are not debugged . There are problems with the 
c                      implementation. Do not use these !!
c              = 99 NCLASS neoclassical transport model including              @
c                      bootstrap current and resistivity                       @
c              = 100 NCLASS neoclassical transport model including             @
c                      fast ion contributions                                  @
c              = 110  Sauter bootstrap model 
c              = 111  Sauter bootstrap and Sauter form of resistivity
c                     is used in code
c                     (no option for includeing the ion part of the fast
c                      ion bootstrap is implemented).
c              = 112  Sauter same as 110 but uses z of first primary ion
c                     rather than zeffa . Note that zeff is used in nvloop
c                     (Option 111 is not available if jboot = 112 is used)
c
c  mult_curboot         This is a bootstrap current multiplier
c                       used only for NCLASS models 95,96.
c                       For other models  you input the multipliers for
c                       the individual components of the bootstrap
c                       through the wneo(4,j) array. The
c                       NCLASSS 95,96  model does not break down the 
c                       bootstrap in this way so we use the averall
c                       mult_curboot factor for this case.

 
c% resistive               character variable selects the resistivity          @
c                          to use according to:                                @
c                          = 'hinton' old (small inverse aspect ratio)         @
c                            'hinton' is the DEFAULT, valid for eps < 0.05     @
c                          = 'hirshman'   (arbitrary     aspect ratio)         @
c                          = 'kim'        (arbitrary     aspect ratio,         @
c                                          valid in banana regime only)        @
c                            'kim' is NOT CURRENTLY IMPLEMENTED                @
c                          = 'nclass'
c                                                                              @
c% ftcalc                   character variable,                                @
c                           ftcalc = 'exact'  means  evaluate the integrals    @
c                           along the psi contours to determine the            @
c                           non collisional trapped particle fraction.         @
c                           ftcalc = 'analytic' means use the (circular)       @
c                           inverse aspect ratio expansion:                    @
c                           ft = 1.0-(1-delta)**2/(SQRT (1.0-delta**2)/        @
c                               (1.0 + 1.46 * SQRT (delta))                    @
c                           where delta is horizontal inverse aspect ratio.    @
c                           The formula of Lin-Lui and Miller is used if       @
c                           jhirsh = 95 or 96                                  @
c                           ftcalc = 'analytic' is the DEFAULT                 @
c                     NOTE: this affects bootstrap current calculations        @
c                           for models with jhirsh > 0 and also                @
c                           the plasma resistivity if resistive = 'hirshman'   @
c                                                                              @
c% cer_ion            identifies the ion species which the charge              @
c                     exchange recombination (CER) diagnostic measured.        @
c                     The input angrot is then assumed to be Vtor/R along      @
c                     the Z=0 cord for the cer_ion species. This is used       @
c                     to compute the radial electric field if jhirsh > 95      @
c                     The poloidal and toroidal rotation of the main and       @
c                     cer ions is output to the netCDF file iterdb.nc .        @
c                     cer_ion can take the same values as namep, namei.        @
c                     The default = ' ' . 
c
c
c                     Note that the angrot profile is defined
c                     as the toroidal bulk  rotation. Hence the above
c                     interpretation as the cer_ion rotation is inconsistent 
c                     with the Onetwo usage of this quantitiy. 
c                     The cer_ion only comes into play when the ifs-ppl
c                     confinement model is selected so dont select this model
c                     to avoid the inconsistency HSJ 06/23/04  
c                                                                              @
c% squeeze            If .true. and cer_ion is set then ion orbit              @
c                     squeezing is applied to the NCLBOOT calculations.        @
c                     default is .false.                                       @
c                                                                              @
c% wneot                  Weight for neoclassical enhancement of               @
c                           resistivity due to trapped electrons               @
c% w1typ(nbctim),         Weights for empirical ("take-your-pick")             @
c% w2typ(nbctim),         transport. The diagonal elements may be              @
c% w3typ(nbctim),         time-dependent if w*typ(2) .ne. 0 and                @
c% w4typ(nbctim)          nbctim > 1 (see nbctim below)                        @
c% w1typmin - w4typmin    time independent quantities. These give the minimum  @
c                         value that each of the w1-4 typ profiles is allowed  @
c                         to have over the rho grid. Each local value of these @
c                         profiles will be replaced by this minimum value if   @
c                         the calcualted value is less than this minimum.      @
c                         (typically this could happen due to wiggly splines)  @
c                         Defaults are -1.e100 (so effectively off)            @
c  w2cg ,w3cg             critical gradient values, kev/cm, for including a    @
c                         critical gradient dependence in te and ti. If the    @
c                         gradient of  te or ti profiles is > w2,3cg in        @
c                         absolute value then xketyp and or xkitype is enchaced  @
c                         by the multiplicative factor w2cgm(for te) or        @
c                         w3cgm (for ti).                                      @
c% w12typ, w13typ         Off-diagonal empirical weights                       @
c% vtyp(nbctim)           Edge velocity (cm/s) for empirical pinch             @
c% typa(i), i=1,16           Exponents for empirical transport                  @
c                           typa(i), i=9,16:                                   @
c                             typa(9) and typa(10) are exponents for the       @
c                             total or local integrated heating power,         @
c                             p, in Mw.  Whether or not p is the local         @
c                             heating (in an annular volume region of          @
c                             width drho) or the global heating for the        @
c                             entire torus is determined by typa(16).          @
c                             if typa(16) .ne. 0, the global value is          @
c                             selected.  if typa(16) = 0, the local value      @
c                             is selected.  p is defined as,                   @
c                               p = typa(11)*pohm   +                          @
c                                   typa(12)*pbeame +                          @
c                                   typa(13)*pbeami +                          @
c                                   typa(14)*prfe   +                          @
c                                   typa(15)*prfi                              @
c                             the multiplier of the transport                  @
c                             coefficient is p**typa(9) or p**typa(10),        @
c                             depending on the value of,                       @
c                                pheat = pohm+pbeame+pbeami+prfe+prfi          @
c                             if pheat = pohm, typa( 9) is used.               @
c                             if pheat > pohm, typa(10) is used.               @
c                             no allowance for fusion currently exists.        @
c                             note that the factors pohm,pbeame,...,           @
c                             appearing above are local (typa(16) = 0.0)       @
c                             or global (typa(16) > 0.0).  In effect,          @
c                             typa(9) and typa(10) are time-dependent          @
c                             exponents.  This may cause some problems         @
c                             (oscillations?) though none are expected.        @
c                             since typa(9),typa(10) do not depend on          @
c                             the dependent variables (ni,te,..etc.),          @
c                             there will be a one time step lag in             @
c                             switching from typa(9) to typa(10) since         @
c                             the transport coefficients are calculated        @
c                             before the powers are known.  Again, this        @
c                             is not expected to cause difficulties.           @
c                             This argument also applies at the initial        @
c                             time step (pohm and/or paux is zero at the       @
c                             initial time when the transport                  @
c                             coefficients are set up).                        @
c% nw1pro, w1pro          Number of points (nw1pro) and scale factor           @
c                           profile (w1pro) for 'typ' particle diffusion,      @
c                           which may be time-dependent.                       @
c% nvpro , vpro           Number of points (nvpro) and scale factor            @
c                           profile (vpro) for 'typ' pinch velocity,           @
c                           which may be time-dependent.                       @
c% nw2pro, w2pro          Number of points (nw2pro) and scale factor           @
c                           profile (w2pro) for 'typ' electron thermal         @
c                           conductivity, which may be time-dependent,         @
c                           analogous to nqrad and qradin.                     @
c% nw3pro, w3pro          Number of points (nw3pro) and scale factor           @
c                           profile (w3pro) for 'typ' ion thermal              @
c                           conductivity, which may be time-dependent.         @
c% nw4pro, w4pro          Number of points (nw4pro) and scale factor           @
c                           profile (w4pro) for 'typ' momentum diffusivity     @
c                           in units of 10**4 cm**2/sec                        @
c% relaxtyp             under-relaxation parameter for take-your-pick model    @
c                       during the iteration of the corrector steps the Rebut  @
c                       diffusivity will be under-relaxed if this parameter    @
c                       is < 1.  set to 1 for no relaxation.  under-relaxation @
c                       may be necessary to get smooth solution.               @
c                       see also paramater relaxsol and ddebug(3).             @
c% wstd                 Weight for standard drift-ballooning model.            @
c% betlim, betlim0      An input limit on beta.  When the total plasma         @
c                       beta passes beyond betlim, the 'typ'                   @
c                       transport is increased by the factor                   @
c                       (1.0 + (beta-betlim)/betlim*betlim0).                  @
c                       betlim0 give the strength of the beta-limit            @
c                       turn-on.                                               @
c                                                                              @
c% set_chie_chii        Switch added 6/25/98.                                  @
c                       This is a more general and convenient  way                @
c                       of settig chie = const*xchi than using the wneo(2,2)   @
c                       value. If set to a non zero value then                 @
c                       the electron chi will have added to it the term        @
c                       set_chie_chii*chi. The total electron diffusivity      @
c                       will be the sum of all active electron chi models as   @
c                       usual.                                                 @
c  use_stab_flux        Introduce a stabilizing flux for turbulent transport
c                       models based on dbar_stab. 
c                       values of dbar (diffusion doefficients)
c                       are determined internally to yield stability.
c  dbar_stab            if set then  use it instead of internal values.
c
c% neomult              Integer, set to 1 to use neomult1 and neomult2         @
c                       multipliers (default = 0)                              @
c% neomult1             Real quantity. To set electron diffusivity to          @
c                       neomult1*(ion neoclassical diffusivity)                @
c                       for q < qneomult.                                      @
c                       (no default)                                           @
c                       (No additive models are incorporated)                  @
c% neomult2             use this multiplier of q > qneomult                    @
c                       (no default value)                                     @
c% qneomult             use neomult1 for q < qneomult,use neomult2 for         @
c                       q > qneomult  (no default)                             @


c                                                                              @
c% iwangrot             angular momentum diffusivity model switch              @
c                       (ignored in analysis mode)                             @
c                           = -2 (Mattor-Diamond model):                       @
c                             Note this model predicts momentum                @
c                             diffusivity = ion thermal diffusivity            @
c                             and also gives a value for the electron          @
c                             thermal diffusivity.  At present the ion         @
c                             and electron thermal diffusivities are           @
c                             calculated and printed out, but not used         @
c                             otherwise.  Ref: Phys Fluids 31(1988)1180.       @
c                             The Mattor-Diamond model,as implemented here     @
c                             requires the specification of some multipliers   @
c                             as follows:                                      @
c% xmtmdifs(i)                   i = 1,nprim (nprim is number of primary ions) @
c                                xmtmdifs(i) is an arbitrary multiplier of     @
c                                chi mometum (and chi ion thermal diffus.)     @
c                                for primary ion species i. That is,each       @
c                                primary ion species yields a chi given        @
c                                by eq. 43 (above ref.) A weighted sum of      @
c                                these individual diffusivites,with weight     @
c                                xmtmdifs(i),is used as the Mattor-            @
c                                Diamond momentum diffusivity.                 @
c  xmtmdifs(nprim+1)                An impurity enhancement factor.            @
c                                   According to ref                           @
c                                   This factor should multiply the            @
c                                   expression                                 @
c                                   to account for impurities. Normally        @
c                                   this factor would be set to 1.0            @
c  xmtmdifs(nprim+2)                A threshold value for the etai modes.      @
c                                   If etai is less than xmtmdifs(nprim+2)     @
c                                   then the model gives zero diffusivity.     @
c                                   Note: at present there is no soft turn-    @
c                                   on of the etai modes in this model. This   @
c                                   can cause oscillations in the angular      @
c                                   momentum but does not affect the thermal   @
c                                   ion and electron energy equations directly @
c                                   (there is some coupling through the ion    @
c                                   energy equation but this should not        @
c                                   cause any problems. If oscillations do     @
c                                   show up a soft turn on model will have     @
c                                   to be developed in subroutine DIFFUS.)     @
c  xmtmdifs(nprim+3)                  A multiplier for the neoclassical        @
c                                     momentum diffusivity which will be       @
c                                     added to the Mattor-Diamond diffusivity. @
c                                     Thus for example we can set xmtmdifs(    @
c                                     nprim+3) = 1.0 and xmtmdifs(i) = 0.0,i=  @
c                                     1,..nprim, to get the neoclassical       @
c                                     diffusivity with the switch iwangrot     @
c                                     set equal to 2. (of course we would      @
c                                     normally achieve this more directly      @
c                                     by setting iwangrot = 0 to begin with)   @
c                                     Note that the neoclassical diffusivity   @
c                                     is always multiplied by wneo(5,5) so     @
c                                     that the actuall multiplier becomes      @
c                                     wneo(5,5)*xmtmdifs(nprim+3).             @
c  xmtmdifs(nprim+4)                  Multiplier for electron thermal          @
c                                     diffusivity calculated from eq. (51)     @
c                                     (above reference). Not significant at    @
c                                     this time since the electron (and ion)   @
c                                     thermal diffusivites are not used in     @
c                                     any calculations.(printed out only)      @
c                                                                              @
c            THE FOLLOWING MODELS APPLY ONLY TO MOMENTUM DIFFUSIVITY           @
c                          parabolic diffusivity:                              @
c             iwangrot = -3 set non turbulent part to ion neoclassical chi
c                           (additional contributions from other
c                            models will be added in if the model
c                            provides such a chi)
c             iwangrot = -1                                                    @
c                          the model is defined by                             @
c                          xkangrot(j) = (xmtmdifs(1)-xmtmdifs(2))*(1.0-       @
c                          roa(j)**xmtmdifs(3))**xmtmdifs(4)-xmtmdifs(2)       @
c                                                                              @
c             iwangrot = 0 (default)                                           @
C                          the model is neoclassical (see subroutine DIFFUS)   @
c                                                                              @
c             iwangrot = 1                                                     @
c                          the angular momemtum diffusivity is constant        @
c                          xkangrot(j) = xmtmdifs(1) for all j                 @
c                                                                              @
c             iwangrot = 2                                                     @
c                          the model is linear                                 @
c                          xkangrot(j) = xmtmdifs(1)*roa(j)+xmtmdifs(2)        @
c                                                                              @
c             iwangrot> = 3                                                    @
c                          the model is a spline                               @
c                          in this case a minimum of three knots must be       @
c                          specified in the array rmtmdifs. The first knot     @
c                          must be at 0.0 and the last knot must be at 1.0     @
c                                                                              @
c  NOTE                    the units of xkangrot are assumed to be             @
c                          cm**2/sec so set xmtmdifs accordingly.              @
c% xmtmdifs                                                                    @
c% rmtmdifs                 input arrays that define the angular momentum      @
c                           diffusivity model (see iwangrot above).            @
c                                                                              @
c





c-------------------------------Magnetic Breaking Model -------------------
c empirical model for magnetic drag on toroidal rotation speed of plasma
c the model is 
c     cb_mgbr*(berrqn_mgbr(t)/wqn_mgbr) *exp(((rho-rhoqn)**2/lqn_mgbr**2)
c
c
c     the code determines rhoqn (the rho value that corresponds to qn_mgbr)
c     
c input parameters in the model are user selected:
c cb_mgbr   this parameter should be choosen so that the units on
c           the above term represent a torque density ( g/(cm*sec**2) )
c qn_mgbr   q value at breaking surface
c lqn_mgbr  scale parameter (~width of exponential)
c 
c time_mgbr(n_mgbr)    times in secs n_mgbr is defined in tordlrot.i
c berrqn_mgbr(n_mgbr)  error field (GAUSS !!!!) at times in time_mgbr

      toff_mgbr = -1.e35
      do j=1,n_mgbr
       time_mgbr(j) = toff_mgbr
      enddo

      use_stab_flux = .FALSE.
      dbar_stab = 0.0D0



c --------------------- Rebut-Lallia-Watkins Model ---------------------       @
c                                                                              @
c% wrebut               weight for Rebut-Lallia-Watkins diffusivities          @
c                       and diffusion coeffients. set wrebut = 0 (default)     @
c                       if this option is not wanted. set wrebut = 1.0 to      @
c                       get full Rebut model. wrebut multiplies the model      @
c                       so an arbitrary part can be added.                     @
c% relaxrebut           under-relaxation parameter for Rebut-Lallia-Watkins    @
c                       model.  during the iteration of the                    @
c                       corrector steps the Rebut diffusivity will be          @
c                       under relaxed if this parameter is less than 1.        @
c                       set to 1 for no relaxation. under relaxation           @
c                       may be necessary to get smooth solution.               @
c                       see also paramater relaxsol and ddebug(3).             @
c                       note that relaxrebut refers only to the                @
c                       anomalous part of the RLW diffusivity.                 @
c                       if it is desired to relax (in addition to,             @
c                       or in place of, just the anomalous RLW value),         @
c                       then use ddebug(3).                                    @
c% qrebsmth       = 0 no effect,                                               @
c                 = 1 smooth the q profile used in Rebut-Lallia model          @
c                    (note that q is not smoothed elsewhere in the code)       @
c                                                                              @
c% tirlw          = 0 (default) means use the real ion temperature,            @
c                      ti, for the evaluation of the RLW model.                @
c                 = 1  use an effective ti, defined as                         @
c                                                                              @
c                            ti*nth + 0.66*(wbeam+walp)                        @
c                      ti =  -----------------------------                     @
c                                (nth+nbeam+nalp)                              @
c                                                                              @
c                      where nth is the total thermal ion density              @
c                      (includes impurities)                                   @
c                      nbeam is the fast ion density due to beams              @
c                      nalp is the fast alpha particle density                 @
c                      wbeam is the fast ion stored energy density             @
c                      due to beams                                            @
c                      walp is the stored fast alpha particle energy density   @
c                                                                              @
c% rlw_model = 'old'  ==>  original 1988 RLW model                             @
c            = 'new'  ==>  with corrections for ion rho-star scaling (DEFAULT) @
c                                                                              @
c               REFERENCES:                                                    @
c                                                                              @
c               Rebut, P.H., Lallia, P.P., Watkins, M.L., Plasma Physics       @
c               and Controlled Nuclear Fusion Research, (Proc. 12th Int.       @
c               Conf., Nice, 1988), Vol. 2, IAEA, Vienna (1989) 191            @
c                                                                              @
c               Rebut, P.H., Watkins, M.L, Gambier, D.J., Boucher, D.,         @
c               Physics of Fluids b 3 (1991) 2209                              @
c                                                                              @
c               Boucher, D., personal communication                            @
c                                                                              @
c               For the HSJ 8/22/95 update (rlw_model = 'new'),                @
c               the reference is Seville, IAEA, 1994, Rosenbluth et al.        @
c                                                                              @
c -------------------- HSIEH ("SHAY") MODEL OF CONDUCTIVITY ------------       @
c                                                                              @
c% wshay               multiplier for Hsieh model of ke. set wshay = 0.0       @
c                      to not use this model ,set wshay = 1.0 to use full      @
c                      Hsieh model (in addition to any other models that       @
c                      are turned on).                                         @
c% smult               constant in Hsieh model (units of this constant         @
c                      depend on exponents choosen for Hsieh model.            @
c% scsmult             logical,if true constant multiplier for Hsieh           @
c                      model, smult, is to be found from power balance         @
c                      analysis (using least squares method)                   @
c                      NOTE that if scsmult is set then                        @
c                                   itte=0                                     @
c                                   itti=0                                     @
c                                   iterate_beam = .true.                      @
c                      is internally set by the code.                          @
c% snexp               exponent of electron density in Hsieh model             @
c% sbpexp                            bp                                        @
c% srexp                             r                                         @
c% sbigrexp                          R                                         @
c% stexp                             Te                                        @
c% sdtdrexp                          d(Te)/dr                                  @
c% srin                starting value of NORMALIZED rho at which Hsieh model   @
c                      becomes active                                          @
c% srout               ending value of NORMALIZED rho at which Hsieh is active @
c% srincev             these parameters are similar to srin and srout          @
c% sroutcev            they are used only if scsmult = true. in this           @
c                      case srincev is the starting value of NORMALIZED        @
c                      rho at which data for the least squares evaluation      @
c                      will be used. Similarly srout is the ending value       @
c                      of NORMALIZED rho at which data for the least           @
c                      squares evaluation will be used. It is required         @
c                      that srincev .ge. srin and sroutcev .le. srout.         @
c                      Also we must have srincev .le. sroutcev. If             @
c                      srincev = sroutcev then the single nearest point        @
c                      in NORMALIZED rho will be used. In this case            @
c                      a least squares anaylis is not possible because         @
c                      there are not enough degrees of freedom.                @
c% suserho             logical, if true, then use flux surface average form    @
c                      of Hsieh model, including (grad-rho)**2 term.           @
c% skimult             set  ion k =skimult * shay (ke)                         @
c                      (on top of any other ion ki that is active).            @
c                      of Hsieh model, including (grad-rho)**2 term.           @
c                                                                              @
c ================================                                             @
c 1995 Upgrade of the Hsieh Model:                                             @
c ================================                                             @
c                                                                              @
c% ishayform    = 1    Use the dimensionally correct version of the model      @
c               = 0    Use the old version with funky units                    @
c               NOTE:  The dimensionally correct form is valid only for the    @
c                      default exponents                                       @
c% skimass             Since the dimensionally correct version of the Hsieh    @
c   (amu)              model for electron thermal conductivity now has a       @
c                      1/SQRT(elec mass) dependence, the ion multiplier        @
c                      SKIMULT should have a similar dependence on the ion     @
c                      mass. This is accomplished by using an effective        @
c                      atomic mass for the ions, in a.m.u. (The ratio          @
c                      sqrt(me/mp) is assumed to be included in the value      @
c                      of SKIMULT)                                             @
c                      For a deuterium plasma, use skimass= 2.0 (default),     @
c                      a hydrogen plasma use 1.0, etc. Its up to you to        @
c                      decide what to use in a multiple-species plasma!        @
c                      NOTE: skimass is auotmatically reset to 1.0 for         @
c                            ISHAYFORM = 0 to keep old form of model           @
c                                                         -Daniel Finkenthal   @
c                                                                              @
c   ------------------------------------------------------------------------   @
c   ---------------------------- WEILAND MODEL -----------------------------   @
c   This model is described in:                                                @
c                                                                              @
c ...    J. Weiland and H. Nordman, "Drift wave model for inward               @
c ...    energy transport in tokamak plasmas," Institute for                   @
c ...    Electromagnetic Field Theory and Plasma Physics,                      @
c ...    Gothenburg, Sweden, (1992) CTH-IEFT/PP-1992-13 ISSN.                  @
c                                                                              @
c ...    J. Weiland, A.B. Jarm\'{e}n, and H. Nordman, "Diffusive               @
c ...    particle and heat pinch effects in toroidal plasmas,"                 @
c ...    Nucl. Fusion {\bf 29} (1989) 1810--1814.                              @
c                                                                              @
c ...    and many more (see notes of Batemen)                                  @
c                                                                              @
c                 ONETWO INPUT SWITCHES FOR THE WEILAND MODEL                  @
c    ------------------------------------------------------------------------  @
c% include_weiland = 1 must be set (in BOTH  analysis and simulation           @
c                      mode) in order to activate this model.                  @
c                                                                              @
c               the Weiland particle and energy diffusivities will be used     @
c               in simulation mode for any dependent variable for which the    @
c               simulation mode flag is set (i.e., itenp, itte, itti) if       @
c% wweiland .ne. 0.0                                                           @
c               in this case the term                                          @
c                    wweiland*(d and/or chie and/or chi weiland)               @
c               WILL BE ADDED TO any other models that are turned on.          @
c                           (So make sure other models are turned off )        @
c                                                                              @
c               the switch relaxrebut is used for the Weiland model as well    @
c               see comments under the RLW model for meaning of relaxrebut.    @
c                                                                              @
c                                                                              @
c   In ANALYSIS MODE,or if wweiland =0.0, the d,chie,and chii for the          @
c   Weiland model are calculated,printed and plotted for inspection            @
c   (provided that include_weiland=1 of course)  BUT ARE NOT USED              @
c   IN ANY OTHER WAY.                                                          @
c                                                                              @
c ----------------------------------------------------------HSJ/2/16/96-----   @
c                                                                              @
c  OCTOBER 1995 MODIFICATION OF ELECTRON AND ION THERMAL CONDUCTIVITIES        @
c  DUE TO A PHENOMENOLOGICAL CRITICAL GRADIENT IN THE TOROIDAL ROTATION        @
c  SPEED PROFILE                                                               @
c                    RGC (Rotational-Gradient-Critical) Model                  @
c             also known as the "DeBoo Phenomenological Model"                 @
c                                                                              @
c     This model multiplies the anomalous part of the electron and/or ion      @
c     thermal conductivity by the factor                                       @
c                                                    rgca                      @
c                       rgc_mult(rho)= ----------------------------------      @
c                                      rgcb+rgcc*[abs(rg)/abs(rgc)]**rgcd      @
c                                                                              @
c     Here rg and rgc are gradients (derivatives wrt rho, not normalized rho)  @
c     and rgca, rgcb, rgcc, rgcd are input constants.                          @
c     abs(rg(rho)) is taken as the maximum of abs(rg(rho)) and abs(rgc(rho))   @
c     so that their ratio is never less than 1. To reduce rgc_mult to 1        @
c     you will typically want to take rgca=rgcb+rgcc, but that's up to you.    @
c                                                                              @
c     The idea behind using this model is that rgc(rho) is determined a priori @
c     over some time interval (i.e., H mode).  Then the model is applied over  @
c     a different time interval (i.e., VH mode), or a different case with the  @
c     now known  rgc(rho) factor. Because we cannot switch from analysis to    @
c     simulation mode during a onetwo run this typicaly means that you must do @
c     two runs, one to determine rgc, which would typically be done in         @
c     analysis mode, and another run to use rgc in simulation mode.            @
c     There is an exception to this rule because rgc is determined a priori    @
c     (i.e., before we take any time steps away from time0).                   @
c     Hence you can set up your input                                          @
c     file with kinetic profiles that are given for times before time0,        @
c     and of course at least all the way  up to timmax (use the  bctimes       @
c     vector to specify the times involved).                                   @
c     Using times_rgc and timee_rgc you then select over what time interval    @
c     you want rgc(rho) averaged. Typically you would not want the interval    @
c     (times_rgc, timee_rgc) to overlap with (time0, timmax) but again that's  @
c     up to you.                                                               @
c                                                                              @
c%   irgc          turns on/off  the model as follows:                         @
c          irgc=0  (default) model not used                                    @
c          irgc=-1 don't apply the model,just determine                        @
c                  the critical gradient factor rgc(i),i=1,2..nj               @
c                  and print it out in outone {times_rgc and                   @
c                  timee_rgc must be set within bctime(1) to bctime(nbctim)}   @
c          irgc=+1 Apply the  critical gradient factor                         @
c                  to  the electron and ion  thermal conductivity              @
c                  selected. If one or both times_rgc and timee_rgc are set to @
c                  times outside the interval bctime(1) to bctime(nbctim)      @
c                  it is taken to mean that rgc(rho) is supplied in the inone  @
c                  file and should not be calculated.                          @
c                  If times_rgc and timee_rgc are such                         @
c                  that they form a subset of bctime(1),bctime(nbctim)         @
c                  then rgc will be calculated first, before the simulation    @
c                  is performed.                                               @
c                  (irgc is independent of ifsflag, so you could have          @
c                  both flow shear suppression effects on)                     @
c                                                                              @
c%   rgce_mult         take electron conductivity to be                        @
c                  xketotal=xkeano*rgce_mult*rgc_mult+xkeneo                   @
c%   rgci_mult         take ion conductivity to be                             @
c                  xkitotal=xkiano*rgci_mult*rgc_mult +xkineo                  @
c                  here xkeano and xkiano are the anomalous part               @
c                  of the thermal conductivity (i.e., HSIEH, RLW, etc.)        @
c                                                                              @
c                               AT PRESENT RGC_MULT                            @
c                         IS  WIRED UP TO THE (ANOMALOUS PART)                 @
c                     OF THE HSIEH, RLW, AND TAKE-YOUR-PICK MODELS!!           @
c%   rgca                                                                      @
c%   rgcb                                                                      @
c%   rgcc                                                                      @
c%   rgcd          input factors see above                                     @
c%   rgc(i)        i=1,2..nj, see explanation above. NOTE if rgc is input      @
c                  in inone teh gradient should be wrt normalized rho.         @
c%   times_rgc     start time for determination of critical gradient           @
c                  factor rgc (sec)                                            @
c%   timee_rgc     end time for determination of critical gradient             @
c                  factor rgc (sec)                                            @
c                  YOU MUST SUPPLY APPROPRIATE EMPIRICAL TOROIDAL ROTATION     @
c                  SPEED PROFILES (SEE DESCRIPTION OF ANGROTIN...)             @
c                  THAT COVER THE RANGE OF TIMES SELECTED BY times_rgc         @
c                  and timee_rgc in order to determine rgc!                    @
c                                                                              @
c%   rgc_string    character variable (up to 128 characters in length)         @
c                  will be echoed on output listing exactly as input           @
c                                                                              @
c   ------------------------------------------------------------------------   @
c   ------------------------------------------------------------------------   @
c   ---------------------------- IFS MODEL ---------------------------------   @
c   This model is described in:                                                @
c                                                                              @
c ...    Dorland et. al. IFSR#677 "Comparison of Nonlinear                     @
c        Toroidal Turbulence Simulations with Experiments"                     @
c                                                                              @
c        Kotchenreuther et al. "Quantitative Predictions of Tokamak            @
c        Energy Confinment From First Principles Simulations with              @
c        Kinetic Effects." Phys. Plasmas 2,(6),june 1995                       @
c                                                                              @
c         ONETWO INPUT SWITCHES FOR THE IFS MODEL                               @
c                                                                              @
c         (THIS MODEL IS INCOMPLETE AT THIS TIME, FOR EXAMPLE IT WILL          @
c         NOT!!!!  HANDLE NEGATIVE SHEAR . WE WILL HAVE TO AWAIT THE           @
c         FUTURE DEVELOPMENTS FROM DORLAND AND KOTCHENREUTHER)                 @
c         Full 2d model is implemented. For 1d however we are currently        @
c         limited to cirular plasmas ONLY!!                                    @
c                                                                              @
c    ------------------------------------------------------------------------  @
c% include_ifs = 1  must be set (in BOTH  analysis and simulation modes)       @
c                   in order to activate this model (default = 0)              @
c               the IFS energy diffusivities will be used                      @
c               in simulation mode for any dependent variable for which the    @
c               simulation mode flag is set (i.e., itte, itti) if              @
c     flow shear suppression with the Staebler-Hinton model is possible in     @
c     the IFS model, see switch ifsflg                                         @
c                                                                              @
c% dorl_kotch .ne. 0.0 (default = 0.0)                                         @
c              in this case the term                                           @
c                    dorl_kotch*(  chie and/or chi IFS )                       @
c               WILL BE ADDED TO any other models that are turned on.          @
c                           (So make sure other models are turned off )        @
c                                                                              @
c               the switch relaxrebut is used for the IFS model as well        @
c               see comments under the RLW model for meaning of relaxrebut.    @
c   separate multipliers for the electron and ion channels are now (12/16/98)  @
c   available :                                                                @
c% dorl_kotche .ne. 0.0       (default = 0.0) electrons                        @
c% dorl_kotchi .ne. 0.0       (default = 0.0) ions                             @
c   if dorl_kotch is set to a non zero value the code will set                 @
c               dorl_kotche = dorl_kotch                                       @
c               dorl_kotchi = dorl_kotch                                       @
c   hence to use dorl_kotche and dorl_kotchi separately make sure              @
c   dorl_kotch = 0 ( which is the default)                                     @
c   note that dorlkotchi is used as a multiplier for ion and momentum          @
c   diffusivities as well.                                                     @
c                                                                              @
c% exbmult_ifs     E X B flow shear suppression multiplier(exclusive to        @
c                  the IFS model). This value is                               @
c                  defaulted to 0.0 meaning flow shear sup. is turned off.     @
c                  set to 1.0 to get the full effect, etc.                     @
c                  note that cer_ion must be set also for this to work         @
c                  use jhirsh =95,96 with this option                          @
c                                                                              @
c                                                                              @
c% wneo_elct_ifs   normally the IFS electron diffusivity is assumed            @
c                  to be sitting on top of wneo(2,2)*xchie_neo                 @
c                  set wneo_elct_ifs to make the ELECTRON diffusivity          @
c                  sit on top of wneo_elct_ifs*xchii_neo                       @
c   In ANALYSIS MODE, the chie, and chii for the                               @
c   IFS model are calculated,printed and plotted for inspection                @
c   (provided that include_ifs =1 of course)  BUT ARE NOT USED                 @
c   IN ANY OTHER WAY.                                                          @
c                                                                              @
c ------------------------------------------------- HSJ/2/10/97 --------       @

c NTV effect
c -------------------------------------------C.K. Pan /01/03/10 -------
c include_ntv =1 to turn on this effect
c c_p             neoclassical coefficient
c delta_b_sqr     the volume-averaged squared magnetic perturbation strength
c mp_polnum,
c mp_tornum       the poloidal and toroidal number of the magnetic perturbation
      include_ntv = 0
      c_p = 0.0
      delta_b_sqr = 0.0
      mp_potnum = 0
      mp_tornum = 0
c---------------------------------------------------------------------------      
c ---------------------------------------------------------------------------  @
c Paleoclassical transport
c see "Paleoclassical transport in low collisionality toroidal palsmas"
c  to be published. J.D> Callen Unv Wisc., Oct. 25, 2004
c
c  include_paleo   =1 to turn on this model
c
c

c  paleo_mult      =  arbitrary multiplier for xchie_paleo

c  paleo_pinch_mult      =  arbitrary multiplier for paleo pinch term
c                           The total multuplier for the pinch term is
c                           taken as paleo_mult*paleo_pinch_mult.
c                           Hence you can eliminate the pinch term by setting
c                           paleo_pinch_mult = 0.0 
c
c  paleo_pinch_mult2     = an additional multiplier for the pinch term
c                          for rho < r_paleo_pinch_mult the multiplier on
c                          the picnh term is 
c                          paleo_mult*paleo_pinch_mult*paleo_pinch_mult2
c
c  r_paleo_pinch_cut       value of normalized rho for limiting the paleo 
c                          flux and chi. For rho less than this value a 
c                          limit is palced on the transport. Supplied by
c                          Callen email,02/022/05.  -- HSJ
c ------------------------------------------------------------------------------
c
       include_paleo = 0

       paleo_mult = 1.0
       paleo_pinch_mult = 1.0     
       paleo_pinch_mult2 =1.0
       r_paleo_pinch_cut = -1.


c                                                                              @
c% ddebug(i)               i = 1,2,..50                                        @
c                          these are special switches used to turn on          @
c                          various models. We will eventually document         @
c                          all these switches as the need arises.              @
c                          The original programmer did not supply any          @
c                          info on these switches so the following is          @
c                          my interpretation of the ones I have used --- HSJ   @
c                                                                              @
c  ddebug(1)               unknown (not used?)                                 @
c  ddebug(2)               sets the value of jsteps (note floating             @
c                          to integer conversion). The diffusion               @
c                          coefficients d(i,k,j), i,k = 1,2..nion+2,           @
c                          i,k = nion+4,j=1,2,..nj,are averaged over           @
c                          2*jsteps+1 spatial intervals. the rbp equation,     @
c                          nion+3, is not averaged. (the average can be        @
c                          just jsteps forwards or backwards, see ddebug(10))  @
c  ddebug(3)               a relaxation parameter for d(i,j,k). Normally       @
c                          you would underrelax d so 0 .lt. ddebug(3)          @
c                          .le. 1 makes sense here. default is 1.0             @
c                          (1.0 means no relaxation, 0.0 means take            @
c                          all of the old d and none of the new d).            @
c    NOTE:                 ddebug(2) and ddebug(3) are also in effect          @
c                          for qexch, the anomalous energy exchange term       @
c                          between electrons and ions. qexch is zero           @
c                          unless wstd, the drift ballooning model, is         @
c                          turned on.                                          @
c  ddebug(10)              has an effect only if ddebug(2) is set.             @
c                          in this case the following is done:                 @
c                          a) if ddebug(10) = 0.0   (default)                  @
c                               then the diffusion coefficient is averaged     @
c                               symmetrically about each grid point            @
c                          b) if ddebug(10) .gt. 0.0                           @
c                               then the diffusion coeffient is averaged       @
c                               only in the direction of increasing rho        @
c                          c) if ddebug(10) .lt. 0.0                           @
c                               then the diffusion coeffient is averaged       @
c                               only in the direction of decreasing rho        @
c  ddebug(11)              integer, used if ddebug(3) .ne. 1.0.
c                          if Mod(n,ddebug(11)) = 0 then the realaxation
c                          of the diffusion coefficient is restarted.
c                          (n is the time step number). Normally ddebug(11)
c                          =1 would be used but more than one time step can
c                          be used in the relaxation by setting ddebug(11)
c                          .gt. 1.
c  ddebug(12)              the number of grid points over which
c                          the solution (te,ti rotation) will be symmetrically 
c                          averaged is 2*ddebug(12)+1
c                          Normally no averaging is to be done so ddebug(12)
c                          should be set to 0.
c                          Small values like 1,2,3 may be required to match some codes.
c
c  ddebug(13)             egamma_exp ,E x B shear rate, averaging in glf23.
c                         If ddebug(13) is
c                         non zero then 2*ddebug(13) + 1 spatial points
c                         are averaged to get gamme_exp(j)
c  ddebug(14)             alpha_exp MHD alpha averaging in glf23. If 
c                         ddebug(14) is
c                         non zero then 2*ddebug(14) + 1 spatial points
c                         are averaged to get alpha_exp(j)
c
c  ddebug(15)             shat _exp magnetic shear averaging in glf23. If 
c                         ddebug(15) is
c                         non zero then 2*ddebug(15) + 1 spatial points
c                         are averaged to get shat _exp(j)
c  ddebug(16-48)          used in cray307.f
c
c
c  ddebug(48)            set equal to 1 to get diagnostic paleoclassical
c                        transport output in file outone . 
c                        (include_paleo must be =1)
c ----------------------------------------------------------------------       @
c ------------------------------------------------------------------------     @
c ------------------------------------------------------------------------     @
c --------- modified Bohm-gyro Bohm internal transport barrier model           @
c                                                                              @
c% include_itb   set to 1 to turn this model on                                @
c                multipliers see the documentation for definitions.            @
c                       ce0_mgb,alpha_mgb,ci1_mgb,                             @
c                       ci2_mgb,ce_bgb,cfe_mgb,cfe_bgb,                        @
c                       cfi_mgb,cfi_bgbc1_g,c2_g,c_theta,                      @
c                                                                              @
c                                                                              @
c   ---------------------------- GLF23 MODEL --------------------------------- @
c   J. Kinsey, 7/19/00                                                         @
c   GLF23 transport model for energy and toroidal momentum transport.          @
c   This version does not include density transport or poloidal                @
c   momentum transport.                                                        @
c                                                                              @
c   The model is described in:                                                 @
c       R.E. Waltz, G.M. Staebler, W. Dorland, G.W. Hammett,                   @
c       M. Kotschenreuther, and J.A. Konings,"A Gyro-Landau-fluid              @
c       Transport Model," Phys. Plasmas 4, 2482 (1997).                        @
c                                                                              @
c         ONETWO INPUT SWITCHES FOR THE GLF23 MODEL                            @
c    ------------------------------------------------------------------------  @
c% include_glf = 1 must be set (in BOTH analysis and simulation modes)         @
c                in order to activate this model.                              @
c                (default = 0)                                                 @
c  glf_debug     0 (default) no effect,integer quantitiy                       @
c                1 run in debug mode (for testing only) makes a single call    @
c                  to glf23( after the beam is converged),
c                  writes results to file glfdebug and exits.                  @
c  glf23_ncpus   = number of cpus to use in glf23 clacs
c                   Valid only on multiple cpu machines with mpi
c  write_glf_namelist integer,0,1, or 2,( default 0)
c                   0 means do not write a namelist (in file
c                     glf23_namelist) 
c                   1 means do write the namelist.
c                    This namelist is a convenient way of getting
c                     a single time slice of information from onetwo
c                    into the glf23 stand alone code (as released to
c                    the NTCC). In that stand alone code the
c                    executable testglf can read this namelist (after changeing
c                    the name from "glf_23_namelist" to "in" 
c                    and setting lprint = 1). In this mode only the
c                    results fromn the most current time are saved.
c                 2  write the namelist and terminate Onetwo
c                   right after the very first time that glf is called.
c                3  write the namelist and terminate Onetwo
c                   at the start of the inital time step.
c                   This differs from the  write_glf_namelist =2 
c                   option in that all sources are (rf,beam)
c                   are called at the initial time first, which changes
c                   some input into glf such as the electron density
c                   due to beam efects.
c  glf23_iglf      =0 run old, prior to Feb. 14,2003 version of
c                     glf23
c                  =1 run new version, post Feb. 14,2003. This version
c                     has some fixes for negative shear cases
c  jroot_iglf      = 6,8,12, number of roots to be found
c
c  jeigen_iglf = 1   you have your choice of
c                    cgg or toms eigenvalue solvers, 0 =cgg,1= toms,
c                    undocumented value 2 = zgeev
c
c
c% iglf_chie     multiplier on model chi-e                                     @
c% iglf_chii     multiplier on model chi-i                                     @
c% iglf_chiv     multiplier on model chi-phi                                   @
c% exbmult_glf   E X B flow shear suppression multiplier (default = 0.0)       @
c                for the GLF23 model). This version uses the real geometry     @
c                extension for gamma_E and the nominal multiplier is           @
c                should be 0.60 (for the non renomrmed model                   @
c                this parameter is the same as alpha_e in glf23                @
c  x_alpha_glf       real varible, switch for alpha stabilization,             @
c                =0.0 (default ) stabilization is off                            @
c                =1.0                             on                             @
c                = -1.0  self consistent alpha stabilization                     @
c% limp_glf      integer switch to turn on impurity dynamics (default = 0)     @
c                Set > 0 for impurity effects and real dilution.               @
c                Note: transport coefficients for transport of                 @
c                is currently not available. Also, including impurity          @ 
c                dynamics significantly slows down the code and is             @
c                usually a small effect for small to moderate Zeff.            @
c% ivphi_glf     integer switch to turn on toroidal momentum transport.        @
c                (default = 0). Set > 0 to transport vphi.                     @
c% lprint_glf    integer switch for diagnostic printout (default = 0).         @
c                Set to 1 for output to file 'glf.out' (unit=12).              @
c  ibtflag_glf    switch for effective B field. if 0 use Bt otherwise use Bteff@
c                 (see glf23 documentation for further details) Default = 1
c                 if retuned glf23 is used then this is forced to 1 inside glf23

c  jmm_glf23     set to 1 to generate glf model on a grid point by gridpoint   @
c                case. 
c  i_delay       delay mechanism of some sort? for glf23, see glf23
c                subroutine. normally use 0 
c
c  irotstab      rotational stabilization method, 0 and - 1  means 
c                shear is determined outside of glf23.
c                1 means calculate it internally to glf23.
c                if irotstab = -1 the value determined
c                by the initial conditions is used and 
c                egamma_exp is never changed once it is determined initially!!
c                if irotstab = 0 then the shear is calculated at each time
c                step but smoothing can be applied before it is passed
c                to glf23. Smoothing is done with ddebug(13) and is applied 
c                only when irotstab =-1,0 .
c
c
c following are for testing soltuion only do not use:
c continuation_method   set to 1 to use continuation parameter for resolution  @
c                of the nonlinear xchi model. The increment in chi during each @
c                iteration is then governed by dcparam_in                      @
c iteration_method  used if continuation_methd =0


c random_pert       integer,maximum number of times a random perurbation
c                   of the solution vector is tried before stopping.
c                   (such random perturbations may help free the solver from
c                   a trapped local minimum )
c jelc_clamp        for all grid point 1 .le. j .le. jelc_clamp clamp
c                   the elctron diffusivity at the value give at grid
c                   point jelc_clamp. Not used if set to 0 (default)
c jion_clamp        for all grid point 1 .le. j .le. jion_clamp clamp
c                   the ion  diffusivity at the value give at grid
c                   point jion_clamp. Not used if set to 0 (default)
c CONVERGENCE CONTROL if 

c                set to 1 to use predictor/corrector scheme with non lin solver@
c dcparam_in        increment by which chi is increased each iteration         @
c                (default = 0.1) obsolete!!!!!!
c
c  use_avg_chi    default 0, if set to 1 then glf  chi is set to
c                 fixed avg of chie,chii  over this time step
c                 if iters_freeze is exceeded and freeze_xsn =1 is set
c test_xptor      default 0, set to 1 to read file profile.dat
c                 see sub  xptor_init . (This is a debugging aid only)
c newtonvb    set to 1 for verbose output to screen
c                   (default 0, minimal output to screen)
c itenp_dv
c itene_dv
c itte_dv
c itti_dv
c itangrot_dv     set to 1 to use dv method forni, te,ti,tord. rot. equations
c itangrot_dv     
c                 NOTE dv switches should only be set so that iti_dv =1
c                 all others are untested and are not used in XPTOR
c                 and hence probably shouldn't be used in onetwo either.
c
c
c
c dv_delt         The factor by which the gradient scale length will
c                 be perturbed using the dv method (dfault =0.03)
c ---------------default values added  3/15/01 HSJ ----------------------
      dv_delt =0.03
      init_egamma_exp = 0   !monitors if egamma_exp_initial is set
      jelc_clamp =0
      jion_clamp =0
      random_pert = 10
      newtonvb = 0
      itenp_dv = 0
      itene_dv = 0
      itte_dv =0
      itti_dv =0
      itangrot_dv =0
      dv_method = .false.
      jneo = 0
      include_glf  = 0 
      lprint_glf   = 0
      ivphi_glf    = 0
      limp_glf     = 0
      exbmult_glf  = 0.0
      x_alpha_glf  = 0.0
      iglf_chiv    = 1.0
      iglf_chii    = 1.0   
      iglf_chie    = 1.0 
      jroot_iglf = 8
      glf23_ncpus = 1 
      write_glf_namelist =0
      glf23_iglf  = 1
      jeigen_iglf = 0 
      ibtflag_glf  = 1 
      irotstab = 1 
      i_delay = 0 
      glf_debug = 0
      conv_skip = 0 
      use_avg_chi = 0
      test_xptor = 0                      
      
      chid_glf_d = 0.0 !jmp.den 
      chid_glf_e = 0.0 !jmp.den
      chiv_glf_d = 0.0 !jmp.den 
      chiv_glf_e = 0.0 !jmp.den
      ichid_smooth = 0 !jmp.den
      ichiv_model = 0  !jmp.snu
      ichiv_chii = 0   !jmp.snu
      iglf_eq = 0      !jmp.den
      iglf_idt = 0.0   !jmp.den

c   ---------------------------- MMM95 MODEL --------------------------------- @
c   OM, 24/01/14                                                               @
c                                                                              @
c                                                                              @
c                 ONETWO INPUT SWITCHES FOR THE BRAINFUSE MODEL                @
c    ------------------------------------------------------------------------  @
      include_brainfuse = 1
      brainfuse_path = '/u/meneghini/onetwo/brainfuse.net'

c   ---------------------------- MMM95 MODEL --------------------------------- @
c   JMP, 3/23/06                                                               @
c                                                                              @
c   The model is described in:                                                 @
c   Glenn Bateman et al "Predicting temperature and density profiles           @
c   in tokamaks," Physics of Plasmas, 5 (1998) 1793-1799.                      @
c                                                                              @
c                 ONETWO INPUT SWITCHES FOR THE MMM95 MODEL                    @
c    ------------------------------------------------------------------------  @
c% include_mmm = 1 must be set (in BOTH analysis and simulation modes)         @
c                in order to activate this model.                              @
c                (default = 0)                                                 @
c% iglf_chie     multiplier on model chi-e (defalut = 1.0)                     @
c% iglf_chii     multiplier on model chi-i (defalut = 1.0)                     @                                     @
c% iglf_chiv     multiplier on model chi-phi (defalut = 1.0)                   @
c% exbmult_glf   E X B flow shear suppression multiplier (default = 1.0)       @
c                for the MMM95 model). 

      include_mmm = 0
      immm_chie = 1.0
      immm_chii = 1.0
      immm_chiv = 1.0
      exbmult_mmm = 1.0
  
c ----------------------------------------------------------------------       @
c 
      continuation_method = 0
      non_lin_method   = 2
      iteration_method = 0
      dcparam_in = 0.1
      neomult = 0                                                              %
      jmm_glf23 = 0
      set_chie_chii = 0.0                                                      %
      irgc      =  0                                                           %
      rgca      =  2.0                                                         %
      rgcb      =  1.0                                                         %
      rgcc      =  1.0                                                         %
      rgcd      =  2.0                                                         %
      times_rgc = -1.0e35          ! checked in subroutine CDEBOO              %
      timee_rgc = -1.0e35                                                      %
      do j=1,kj                                                                %
        rgc(j)  = -1.0e35          ! inititalize to non-physical number        %
      end do                                                                   %
c                                                                              %
      do j=1,kddebug               ! zero debug switches                       %
        ddebug(j) = 0.0                                                        %
      end do                                                                   %
      ddebug(3) = 1.0     ! means all of new, none of old value                %
      ddebug(11) = 1.0
      ddebug(12) = 0.0
      do   i=1,4                                                               %
        do j=1,4                                                               %
          wneo(i,j) = 1.0                                                      %
        end do                                                                 %
      end do                                                                   %
      do i=1,4                                                                 %
        wneo(i,5) = 0.0                                                        %
        wneo(5,i) = 0.0                                                        %
      end do                                                                   %
      wneo(5,5)   = 1.0                                                        %
      wneot       = 1.0                                                        %
      w12typ      = 0.0                                                        %
      w13typ      = 0.0                                                        %
      wneo_elct_ifs = 0.0                                                      %
      do i=1,2                                                                 %
        vtyp (i) = 0.0                                                         %
        w1typ(i) = 0.0                                                         %
        w2typ(i) = 0.0                                                         %
        w3typ(i) = 0.0                                                         %
        w4typ(i) = 0.0                                                         %
      end do                                                                   %
c      w1typmin = -1.e100                                                      %
c      w2typmin = -1.e100                                                      %
c      w3typmin = -1.e100                                                      %
c      w4typmin = -1.e100  
      w1typmin = -1.e38                                                        %
      w2typmin = -1.e38                                                        %
      w3typmin = -1.e38                                                        %
      w4typmin = -1.e38                                                        %
      wrebut     = 0.0    !               Rebut-Lallia-Watkins model off       %
      relaxrebut = 1.0    ! relaxation of Rebut-Lallia diffusivity   off       %
      qrebsmth   = 0      ! no smoothing of q in RLW model                     %
      tirlw      = 0          ! use real ti in RLW model                       %
      rlw_model  = 'new'      ! use corrections for ion rho-star scaling       %
      resistive  = 'hinton'   ! in tfact.i                                     %
      ftcalc     = 'analytic' ! in tfact.i                                     %
      wweiland   = 1.0                                                         %
      include_weiland = 0     ! exclude Weiland model calculations             %
      include_ifs     = 0     ! exclude IFS model                              %
      include_glf   = 0        ! exclude glf23 model                           %
      dorl_kotch      = 0.0                                                    %
      dorl_kotche      = 0.0                                                   %
      dorl_kotchi      = 0.0                                                   %
c                                                                              %
c --- modified gyro Bohm defaults                                              %
c                                                                              %
      ce0_mgb   = 1.0                                                          %
      alpha_mgb = 1.0                                                          %
      ci1_mgb   = 1.0                                                          %
      ci2_mgb   = 1.0                                                          %
      ce_bgb    = 1.0                                                          %
      cfe_mgb   = 0.0 ! cfe_mgb or cfe_bgb must be non zero if include_itb=1   %
      cfe_bgb   = 0.0                                                          %
      cfi_mgb   = 0.0 ! cfi_mgb or cfi_bgb must be non zero if include_itb=1   %
      cfi_bgb   = 0.0                                                          %
      c1_g      = 100.0 ! units are kHz                                        %
      c2_g      = 0.0                                                          %
      c_theta   = 0.0                                                          %
      include_itb = 0  ! do not include model by default                       %
c                                                                              %
c --- Hsieh (aka "shay") model defaults                                        %
c                                                                              %
      wshay     =  0.0     ! turn off Hsieh model                              %
      scsmult   = .false.  ! don't calculate multiplier for Hsieh model        %
      suserho   = .false.  ! use simple (rather than flux avg) form            %
      smult     =  0.0     ! constant in Hsieh model                           %
      snexp     =  2.0     ! density exponent                                  %
      sbpexp    =  2.0     ! bp exponent                                       %
      srexp     =  3.0     ! rho exponent                                      %
      sbigrexp  =  1.0     ! R exponent                                        %
      stexp     =  1.0     ! te exponent                                       %
      sdtdrexp  =  2.0     ! dte/dr exponent                                   %
      srin      = -1.0     !  inside radius (normalized rho)                   %
      srout     = -1.0     ! outside radius (normalized rho)                   %
      skimult   =  1.0     ! ki = 1.0 * ke                                     %
      sdenscale =  1.0e13  ! Hsieh model elec dens in units of sdenscale       %
c                                                                              %
c --- New 1995 upgrade of the Hsieh model                                      %
c                                                                              %
      ishayform =  1       ! Use the dimensionally correct form of the model   %
      skimass   =  2.0     ! Effective ion mass (amu) for the ion model        %
c                                                                              %
c ----------------------------------------------------------------------       %
c                                                                              %
      do i=1,16                                                                %
        typa(i) = 0.0                                                          %
      end do                                                                   %
      nw1pro    = 0                                                            %
      nvpro     = 0                                                            %
      nw2pro    = 0                                                            %
      nw3pro    = 0                                                            %
      nw4pro    = 0                                                            %
      do   j=1,ksplin                                                          %
        do m=1,kbctim                                                          %
          w1pro(j,m) = 0.0                                                     %
           vpro(j,m) = 0.0                                                     %
          w2pro(j,m) = 0.0                                                     %
          w3pro(j,m) = 0.0                                                     %
          w4pro(j,m) = 0.0                                                     %
        end do                                                                 %
      end do                                                                   %
      relaxtyp   = 1.0                                                         %
      wstd       =   0.0                                                       %
      betlim     = 100.0                                                       %
      betlim0    =  20.0                                                       %
      jgboot     =   0                                                         %
      jhirsh     =   0                                                         %
      jboot      =   -1
      cer_ion    =  ' '                                                        %
      mult_curboot = 1.d0
      squeeze    =  .false.                                                    %
      iwangrot   =   0                                                         %
c                                                                              %
c --- multipliers for Mattor-Diamond diffusivities (etai modes)                %
c                                                                              %
      fimpurty = 1.0                                                           %
      wetaie   = 1.0                                                           %
      etaioff  = 1.0                                                           %
      do j=1,kion                                                              %
        wetai(j) = 1.0                                                         %
      end do                                                                   %
*                                                                              @
c ----------------------------------------------------------------------       @
c --- PARAMETER FOR CLASSICAL ION THERMAL CONDUCTIVITY                         @
c ----------------------------------------------------------------------       @
c% w3cla            Weight for classical ion thermal conductivity              @
c ----------------------------------------------------------------------       @
c                                                                              @
      w3cla = 0.0                                                              %
c                                                                              @
c ----------------------------------------------------------------------       @
c       PARAMETERS FOR TRANSPORT MODELS ASSOCIATED WITH MHD ACTIVITY           @
c ----------------------------------------------------------------------       @
c% w1isl,      Weights for   island-induced transport                          @
c% w2isl,                                                                      @
c% w3isl                                                                       @
c                                                                              @
c% w1saw,      Weights for sawtooth-induced transport                          @
c% w2saw,                                                                      @
c% w3saw            (simple model)                                             @
c                   the simple model modifies                                  @
c                   a) particle diffusion coeffient by w1saw                   @
c                   b) electron diffusivity by w2saw                           @
c                   c) ion diffusivity by w3saw                                @
c                   d) resistivity                                             @
c                   if the profiles are given (analysis mode) then             @
c                   the  only effect is the change in resistivity which        @
c                   modifies  the current diffusion (provided itxj = 1)         @
c                                                                              @
c                   Weights for sawtooth-induced transport                     @
c                   (mixing model) values are defaulted to 0 (meaning
c                    corresponding action is not taken)                        @
c% w1mix,     =1   mix the thermal particle densities while conserving particles
c% w2mix,     =1         mix the electron temperature while conserving electron energy
c% w3mix,     =1         mix the ion temperature while conserving ion energy         @
c% w4mix,     =1         not used                                                    @
c% w5mix      =1           mix the fast ion density, energy, and parallel momentum   @
c%                  for w(1-5)mix .ne. 0 mixpro determines how the
c%                  associated profile is mixed. The choices are:
c% mixpro        =  2, Generate mixed profiles consistent with helical         @
c                      flux conservation as suggested by Kadomtsev for         @
c                      a single-valued q profile and since extended by         @
c                      Parail and Pereverzev to the case of a double-          @
c                      valued q profile.  Mix temperatures instead of          @
c                      energy densities; then renormalize to conserve          @
c                      energy.                                                 @
c                =  1, Generate mixed profiles consistent with helical         @
c                      flux conservation as above, but mix energy densi-       @
c                      ties.  This can lead to unphysical temperature          @
c                      profiles if the density profiles are peaked and         @
c                      not mixed.                                              @
c                =   0, Generate mixed profiles for densities and temper-       @
c                      atures that are flat.                                   @
c                =  -1, Generate mixed profiles for densities and energies      @
c                      that are flat.                                          @
c% dtemix           Change in electron temperature (keV) on axis during        @
c                   sawtooth oscillation (used when w2mix<0)                   @
c% dtimix           Change in ion temperature (keV) on axis during             @
c                   sawtooth oscillation (used when w3mix<0)                   @
c% fusmix           Relative change (dNdot/Ndot) in neutron rate during        @
c                   sawtooth oscillation (used when w3mix<0 and                @
c                   dtimix = 0)                                                @
c% qmix             Critical value of the safety factor for                    @
c                   triggering a sawtooth disruption                           
c  q0_max           
c  q0_radius        try to limit max value of q on axis to
c                   approximately q0_max. Do his by adding sufficient
c                   seed current to make etor flat out to a normalized
c                   radius of q0_radius. if q0_max <= 0.0 then no
c                   seed current is added.
c  q0_mult          multiplier for seed current for q0_max activation

c% rsmixx           Critical value of the normalized q = 1 radius (rs/a)       @
c                   for triggering a sawtooth disruption                       @
c% s3mix            Relative change (dS/S) in signal from SXR radial           @
c                   diode 3 in old Doublet III diode array during sawtooth     @
c                   oscillation (used when w2mix<0, dtemix = 0, and jsxr=1)    @
c% s71mix           Relative change (dS/S) in signal from SXR vertical         @
c                   diode 71 in old Doublet III diode array during sawtooth    @
c                   oscillation (used when w2mix<0, dtemix = 0, s3mix=0, and   @
c                   jsxr = 1)                                                  @
c% s18mix           Relative change (dS/S) in signal from SXR side diode       @
c                   18 in new Doublet III diode array during sawtooth          @
c                   oscillation (used when w2mix<0, dtemix = 0, and jsxr=2)    @
c% trmix            Ratio of ion temperature change to electron tempera-       @
c                   ture change on axis during sawtooth oscillation            @
c                   (used when w3mix<0, dtimix = 0, and fusmix=0)              @
c% tsmix            Critical value of the sawtooth period (s)                  @
c                   for triggering a single sawtooth disruption                @
c% tdmix(1),(2)     Critical values of the sawtooth period (s)                 @
c                   for triggering double sawtooth disruptions                 @
c% ipmix            Flag for triggering a sawtooth disruption using            @
c                   the prescription of Parail and Pereverzev                  @
c% w0mix            Initial value for the width (cm) of the m = 1 island       @
c                   following a sawtooth disruption                            @
c                   sawtooth disruption                                        @
c ----------------------------------------------------------------------       @
c                                                                              @
      w1isl    = 0.0                                                           %
      w2isl    = 0.0                                                           %
      w3isl    = 0.0 
      nisl     = 0                                                          %
      w1saw    = 0.0                                                           %
      w2saw    = 0.0                                                           %
      w3saw    = 0.0                                                           %
      w1mix    = 0.0                                                           %
      w2mix    = 0.0                                                           %
      w3mix    = 0.0                                                           %
      w4mix    = 0.0                                                           %
      w5mix    = 0.0                                                           %
      dtemix   = 0.0                                                           %
      dtimix   = 0.0                                                           %
      fusmix   = 0.0                                                           %
      mixpro   = 2                                                             %
      qmix     = 0.0                                                           %
      q0_max    = 0.0
      q0_radius = 0.0
      q0_mult = 1.0
      rsmixx   = 0.0                                                           %
      s3mix    = 0.0                                                           %
      s71mix   = 0.0                                                           %
      s18mix   = 0.0                                                           %
      trmix    = 0.0                                                           %
      tsmix    = 0.0                                                           %
      tdmix(1) = 0.0                                                           %
      tdmix(2) = 0.0                                                           %
      ipmix    = 0                                                             %
      w0mix    = 0.0                                                           %
c ----------------------------------------------------------------------       @
c --- STAEBLER-HINTON MODEL                                                    @
c ----------------------------------------------------------------------       @
c                                                                              @
c  Input variables aeh, ael, aih, ail, bh, bl, ch, cl, alfae, alfai,           @
c  betah,sigma and gammah are defined by the following model                   @
c  developed by  G.M. Staebler and F.L. Hinton*:                               @
c                                                                              @
c  (xe is electron energy difusivity)                                          @
c  Xe = A * Xe  where  A = aeh +            ael                                @
c                                    ---------------------------               @
c                                    1 + (alfae * sperp)**gammah               @
c                                                                              @
c  (xi is ion energy diffusivity)                                              @
c  Xi = B * Xi  where  B = aih +            ail                                @
c                                    ---------------------------               @
c                                    1 + (alfai * sperp)**gammah               @
c                                                                              @
c  (Di is particle diffusivity)                                                @
c  Di = C * Di  where  C = bh +             bl                                 @
c                                    ---------------------------               @
c                                    1 + (betah * sperp)**gammah               @
c                                                                              @
c  ( ui is momentum diffusivity)                                               @
c  ui = D * Xi  where  D = ch +             cl                                 @
c                                    ---------------------------               @
c                                    1 + (sigma * sperp)**gammah               @
c                                                                              @
c  If the flow shear suppression flag, ifsflag, is set to 0 the                @
c  model will not be used; if it is set to 1, the coefficients                 @
c  A,B,C,D, will be applied . At the present time                              @
c  this model acts only on the RLW and HSIU models of diffusivity.             @
c  Neoclassical values of xe,xi,di,ui are added to the above turbulent         @
c  values after the above modification is made. Thus for example the           @
c  total momentum diffusivity is                                               @
c                   ui = D*xi +uneo                                            @
c  where uneo is the neoclassical momentum diffusivity ,D is the               @
c  multiplier defined above and xi is either the RLW or HSIU model             @
c  of ion energy diffusivity.                                                  @
c                                                                              @
c  The constant variables are saved in array fs in the following               @
c  order:                                                                      @
c                                                                              @
c%    fs( 1) =  ifsflag                                                        @
c     fs( 2) =  aeh                                                            @
c     fs( 3) =  ael                                                            @
c     fs( 4) =  aih                                                            @
c     fs( 5) =  ail                                                            @
c     fs( 6) =  bh                                                             @
c     fs( 7) =  bl                                                             @
c     fs( 8) =  ch                                                             @
c     fs( 9) =  cl                                                             @
c     fs(10) =  alfae                                                          @
c     fs(11) =  alfai                                                          @
c     fs(12) =  betah                                                          @
c     fs(13) =  sigma                                                          @
c     fs(14) =  gammah                                                         @
c     fs(15) =  lsfctr                                                         @
c%    lsfctr    is a switch used in determination of sperp                     @
c               (see subroutine csperp for actual defn of sperp)               @
c                if lsfctr = 0 then the magnetic shear scale length            @
c                used in sperp is given by                                     @
c                       ls = (R*q**2)/( r*dq/dr)                               @
c                if lsfctr = 1 then we take the constant length                @
c                       ls = Rmajor                                            @
c                                                                              @
c     fs(16) =  xrot       * Multiplier of toroidal frequency derivative       @
c                            component in Sperp term. Set to 1.0 by default.   @
c     fs(17) =  xeden      * Multiplier of electron density   derivative       @
c                            component in Sperp term. Set to 1.0 by default.   @
c     fs(18) =  xsecder    * Multiplier of             second derivative       @
c                            component in Sperp term. Set to 1.0 by default.   @
c                                                                              @
c  *Reference: G.M. Staebler and F.L. Hinton, Particle and Energy              @
c              Confinement Bifurcation in Tokamaks, Phys. Fluids B 5           @
c              (4), p. 1281-1288, April 1993.                                  @
c ----------------------------------------------------------------------       @
c                                                                              @
      ifsflag = 0                                                              %
      aeh     = 0.0                                                            %
      ael     = 1.0                                                            %
      aih     = 0.0                                                            %
      ail     = 1.0                                                            %
      bh      = 0.0                                                            %
      bl      = 1.0                                                            %
      ch      = 0.0                                                            %
      cl      = 1.0                                                            %
      alfae   = 1.0                                                            %
      alfai   = 1.0                                                            %
      betah   = 1.0                                                            %
      sigma   = 1.0                                                            %
      gammah  = 2.0                                                            %
      lsfctr  = 1                                                              %
      xrot    = 1.0                                                            %
      xeden   = 1.0                                                            %
      xsecder = 1.0                                                            %
*                                                                              @
*                                                                              @
c ----------------------------------------------------------------------       @
c --- FLAG FOR IDEAL BALLOONING MODE STABILITY CONSTRAINT                      @
c ----------------------------------------------------------------------       @
c% ibaloo   1, calculate the maximum plasma pressure gradient allowed          @
c              for stability to ideal ballooning modes, and increase           @
c              energy transport, if necessary, to prevent this pressure        @
c              gradient from being exceeded                                    @
c           0, do not perform any calculations related to ideal                @
c              ballooning modes                                                @
c          -1, calculate the maximum plasma pressure gradient allowed          @
c              for stability to ideal ballooning modes, and print this         @
c              out, but do not increase energy transport                       @
c  r_elm       phenomenological elm modeling. set r_elm = to normalized
c              rho value such that elms extens in from rho=1 to rho =r_elm
c  t_elms(n_elms) ,t_elme(n_elms)   start and end time of elms (sec)
c              up to n_elms cycles are allowed
c              The times must be ordered:
c              t_elms(1) < t_elme(1) < t_elms(2) < t_elme(2) ... etc
c               
c  etam_elm(n_elms)    multipier for resistivity during the period t_elms to t_elme
c  itot_elm(n_elms)   multiplier for total current
c ----------------------------------------------------------------------       @
c                                                                              @
      ibaloo = 0                                                               %
      do j=1,n_elms
        etam_elm(j) = 1.0
        itot_elm(j) = 1.0
        r_elm(j) = 2.  !defaulted to off(> 1.)
        t_elms(j) = 1.e30
        t_elme(j) = 1.e30
      enddo                                        
c ----------------------------------------------------------------------       @
c --- MESH PARAMETERS                                                          @
c ----------------------------------------------------------------------       @
c% imesh    0, uniform mesh                                                    @
c              USE ONLY UNIFORM MESH IN TIME-DEPENDENT EQDSK MODE, HSJ         @
c           1, nonuniform mesh with r(j) input                                 @
c           2, nonuniform mesh with dr(j) input                                @
c           3,    uniform mesh in r(j)**2                                      @
c% nj        Number of mesh points (3 to kj)                                   @
c% r(j)     Position (cm) of jth mesh point                                    @
c% dr(j)    Width (cm) of jth mesh interval                                    @
c% rho_edge   A normalized value of rho that defines the plasma boundary       @
c             FOR PURPOSES OF SOLVING THE DIFFUSION EQUATIONS (ONLY)!          @
c             This option is valid only in combination with the tdem mode.     @ 
c             (The required input for TDEM mode is given in NAMELIST 3).       @
c             Normally rho_edge is left at its default value of 1.0            @
c             which means the equations are solved on a grid which is          @
c             determined by the last closed flux surface in the eqdsk.         @
c             You may elect (for the IFS model for example) to solve           @
c             over a smaller interval, expressed in normalized rho space,      @
c                             rho=[0., rho_edge]                               @
c             instead by specifying 0< rho_edge <= 1.                          @
c             The boundary condition at rho_edge for each of the profiles       @
c             run in simulation mode is obtained by interpolation from         @
c             the profiles given in the inone file. Note that this             @
c             means that the input profiles given in inone  MUST!!!!!          @
c             be defined over rho=[0.,1.] even if rho_edge < 1.0 is input.     @
c                                                                              @
c             IT IS ASSUMED THAT PROFILES OF THE QUANTITIES TO                 @
c             BE TRANSPORTED ARE INPUT AS A FUNCTION OF SPACE AND TIME         @
c             JUST LIKE THEY WOULD BE IN ANALYSIS MODE.!!!!!!                  @
c                                                                              @
c             These profiles will                                              @
c             be used only to generate the required boundary condition.        @
c             (For Faraday's law we get the rbp profile as a function of       @
c             space and time from the tdem mode)                               @
c             if you set rho_edge < 0.0 or rho_edge > 1.0 the code will        @
c             ignore your request and use the default of rho_edge=1 !          @
c   fix_edge_ni(m)
c   fix_edge_te(m)                                                              
c   fix_edge_ti(m)                                                              
c   fix_edge_rot(m)
c             set these values to a grid point less than nj to fix the te,ti   @
c             temperatures,and toroidal rotation
c             at input values for j >= fix_edge_te,ti,fix_edge_rot             @
c             the values will be interpolated from the spline input of te,ti   @
c             using this option disables the rho_edge option !!!!!!!!!         @
c             values for te,ti,rotation do not have to be the same             @
c             and can vary over bctime
c             values can also be input in normalized flux as of version        @
c             3/12/02 . For a given profile inputs must all be in grid points  @
c             or all in normalized flux. For the densities,fix_edge_ni(1..kbcim)
c             is the same for all ion species!!!!!
c
c   r_clamp   a normalized rho value . If this is set to a number in (0.0,1.0)
c             then for 0< rho < r_clamp the NEOCLASSICAL  values of the ion
c             temperature diffusivity and the toroidal momentum diffusivity 
c             are linearly extrapolated to chi_clamp at rho =0 
c              implemented only for the glf23 model at this time
c   chi_clamp  cm**2/sec see r_clamp
c% eps_adaptive(i) ,i=1...nion+4                                               @
c            selects the adaptive grid calculations.                           @
c            These calculations dynamically adjust the rho                     @
c            grid spacing as the solution evolves,placing more grid points     @
c            near regions of large variation of the profiles selected by       @
c            eps_adaptive(i).                                                  @
c            the index i corresponds to the                                    @
c            profiles and hence depends on how the input is set up             @
c            i=1 to nprim corresponds to primary ion densities                 @
c            i=nprim+1 to nion (where nion = nprim+nimp)                       @
c            corresponds to the impurity densities                             @
c              i=nion+1 corresponds to Te                                      @
c              i=nion+2 corresponds to Ti                                      @
c              i=nion+3 corresponds to RBP                                     @
c              i=nion+4 corresponds to Toroidal rotation                       @
c             eps_adaptive(i) specifies the weight of each GRADIENT            @
c             to be included in determining the adaptive grid.                 @
c             Note that you can include  analysis mode profiles in this        @
c             scheme.(I dont know why you would want to however).              @
c                                  AN EXAMPLE:                                 @
c                 suppose we have an inone file with nprim=2,nimp=1            @
c                 eps_adaptive (1) =0.0   ! ignore 1'st primary ion density    @
c                                         ! gradient in grid determination     @
c                              (2) =1.0   ! include 2'nd primary ion density   @
c                                         ! gradient in grid determination     @
c                              (3) =0.0   ! ignore impurity ion density        @
c                                         ! gradient in grid determination     @
c                              (4) =2.0   ! include te gradient with twice the @
c                                         ! density gradient importance        @
c                              (5) =1.0   ! include ti gradient                @
c                              (6) =0.0   ! exclude rbp gradient               @
c                              (7) =0.5   ! include rotation gradient at half  @
c                                         ! the importance of density gradient @
c                                         ! this value would be included only  @
c                                         ! if iangrot=1 in the input file     @
c                 Typically only the te gradient would be selected with all    @
c                 other values set to 0.0 But the above generality is coded in @
c                 if eps_adaptive(i)=0 for all i (default) then the            @
c                 adaptive grid calculations are not used.                     @
c                 in addition to (or in place of) the gradient weighting       @
c                 given by eps_adaptive above there is identical weighting     @
c                 for curvature given by curve_eps_adaptive (1 to 7)           @
c                                                                              @
c%    speed_adaptive  a number between 0.0 and 1.0 that controls the           @
c                  rate at which the r grid moves. numbers near zero           @
c                  cause slow evolvement of the grid. speed_adaptive=1.0       @
c                  means switch to the newly calculated grid immediately.      @
c                  This generally is to be avoided since the term              @
c                  drho/dt at constant zeta would then introduce a large       @
c                  time dependent source perturbation.                         @
c                  Default is 0.01. (we want the grid to evolve but            @
c                  generally it should evolve slowly enough that the           @
c                  code doesnt start changing time steps,etc due to            @
c                  sources caused by the changing grid).                       @
c%    gridgenvb   set to 1 to get some output to the screen                    @
c                 regarding the grid generation. (primarily for                @
c                 developers)                                                  @
c%    freeze_adaptive  adaptive grid generation stops at time                  @
c                 t .gt. freeze_adaptive (within a time step,no effort         @
c                 is made to stop exactly at the specified time)               @
c ----------------------------------------------------------------------       @
c                                                                              @
      imesh = 0                                                                %
      nj    = kj                                                               %
      r_clamp = -1.  !no clamping in effect
      chi_clamp = 0.0
      do i=1,kk                                                                %
              eps_adaptive(i) = 0.0                                            %
        curve_eps_adaptive(i) = 0.0                                            %
      end do                                                                   %
       speed_adaptive =   0.05                                                 %
      freeze_adaptive =  -1.0e30   ! disable                                   %
      gridgenvb       =   0                                                    %
      rho_edge        = 1.0                                                    %
      do j = 1,kbctim 
         fix_edge_ni(j)     = DFLOAT(kj + 1)
         fix_edge_te(j)     = DFLOAT(kj + 1) !values are checked in bc_zone    %
         fix_edge_ti(j)     = DFLOAT(kj + 1 )                                  %
         fix_edge_rot(j)    = DFLOAT(kj + 1 )                                  %
      enddo
      te_index = kj+1       !set in set_boundary_conditions to proper values
      ti_index = kj+1
      rot_index = kj+1
      te_index_save = te_index
      ti_index_save = ti_index
      rot_index_save= rot_index
      freeze_te_index =.False.
      freeze_ti_index =.False.
      freeze_rot_index=.False.
      freeze_xte = 0
      freeze_xti = 0
      freeze_xrbp = 0
      freeze_xwr = 0
      freeze_xni = 0
      freeze_xnue = 0
      freeze_alpha_exp = 0
      freeze_xsn = 0

*                                                                              @
c ----------------------------------------------------------------------       @
c --- INITIAL CONDITIONS (PROFILES) AND BOUNDARY CONDITIONS                    @
c ----------------------------------------------------------------------       @
c  The "Old standard"  profiles for the particle densities, temperatures,      @
c     current density, and Zeff have the following form:                       @
c        u(k,j) = ub(m,k) + (uc(m,k)-ub(m,k))*x(j)**alpu(m,k) ,                @
c     with                                                                     @
c        x(j) = 1.0 - (r(j)/r(nj))**gamu(m,k) .                                @
c                                                                              @
c  The central values, boundary values, and exponents at time point            @
c     m (to be described shortly) are as follows:                              @
c                                                                              @
c% enec(m)     , eneb(m)     , alpene(m)    , gamene(m)                        @
c% enpc(m,i)   , enpb(m,i)   , alpenp(m,i)  , gamenp(m,i)                      @
c% enic(m,i)   , enib(m,i)   , alpeni(m,i)  , gameni(m,i)                      @
c% tec(m)      , teb(m)      , alpte(m)     , gamte(m)                         @
c% tic(m)      , tib(m)      , alpti(m)     , gamti(m)                         @
c% xjc(m)      , xjb(m)      , alpxj(m)     , gamxj(m)                         @
c% zeffc(m)    , zeffb(m)    , alpzef(m)    , gamzef(m)                        @
c                                                                              @
c  These produce profiles for the following variables:                         @
c                                                                              @
c% ene(j)      Density of electrons (1/cm**3)                                  @
c% enp(j,i)    Density of ith primary species (1/cm**3)                        @
c% eni(j,i)    Density of ith impurity species (1/cm**3)                       @
c% te(j)       Electron temperature (keV)                                      @
c% ti(j)       Ion temperature (keV)                                           @
c% curden(j)   Current density (see below)                                     @
c% zeff(j)     Effective charge number of ions                                 @
c                                                                              @
c  The units for xjb and xjc are arbitrary, but curden is in A/cm**2,          @
c     which is obtained by normalization to                                    @
c                                                                              @
c% totcur(m)   Total plasma current (A)                                        @
c                                                                              @
c      In 1-1/2-D runs if irguess < 0  and totcur(1) = 0, then                 @
c      totcur(1) is read from the eqdsk file.                                  @
c      if totcur(1) = -1.0e30 then the current is obtained from the third      @
c      namelist of inone, from vector pcurmhd, and if appropriate, pcurmhd     @
c      is interpolated onto the bctime array to generate time-dependent        @
c      totcur vector. bctime and timeqbcd must be commensurate for the         @
c      interpolation to work. If not an error message will be generated.       @
c      totcur has precedence over pcurmhd so pcurmhd will not be used if       @
c      totcur(1) is set to something other than -1.0e30                        @
c      Note that he number of values of totcur must be consistent with nbctim  @
c
c vloop_bc(m),volts  :
c     Instead of totcur, vloop in volts, can be specified
c     at the times given in bctime. If only one value is given then
c     it is assumed that vloop is constant. Vloop_bc takes precedence
c     over totcur. Note that as the resistivity changes vloop will
c     drive different amounts of ohmic current and hence the total current will
c     float. To use this option specify at least one value of vloop_bc
c     in inone.  ( HSJ 6/07/04)
c vloop_bc_time(1:kbctim) ,sec:
c     vloop_bc can be given in bctime or it can be given on a separate
c     time base.
c     vloop_bc_time has the same maximum length as bctime (ie kbctim)
c     the code will detect how many valid entries there are in 
c     vloop_bc_time. The time array must be in 1:1 correspondence
c     with the vloop_bc array but the two arrays do not have to be in
c     any special order. (Onetwo will time sort  the arrays)
c     Obviously some vloop_bc_time(?) .le. time0 .le. timmax .le. 
c     vloop_bc_time(??) must be observed fortwo of the elements.
c     If vloop_bc_time is not used then vloop_bc must correspond
c     to the elements in bctime and bctime itself must be ordered
c     monotonically increasing in time (bctime is not time sorted).
c
c
c vloopvb   set to 1 to get special monitoring output to file
c           vloop_monitor.txt. This is used primarily as a debugging
c           aid.
c
c                                                                              @
c  The particle densities may be specified in different ways, depending        @
c     upon the value of the following parameter:                               @
c                                                                              @
c  bp0_ic   'analytic' or 'eqdsk'
c           see notes on how to use these
c          
c
c% inenez   1, input electron density and zeff, and calculate the              @
c              hydrogen and impurity densities.  up to two                     @
c              hydrogen species and a single impurity species are              @
c              allowed.  for two hydrogen species the fraction                 @
c              of the first is specified by zfrac.  (icenez = 0)               @
c                  Note: In simulation mode the electron density is held       @
c                  fixed in time. Hence as the primary ions diffuse            @
c                      the impurity density changes to keep zeff from          @
c                      changing.                                               @
c                  NOTE: if ifus = -1 ,up to three primary ions may be input   @
c                        see ifus=-1 for details                               @
c           0, calculate electron density and zeff from primary and            @
c              impurity ion densities.  (icenez = 1)                           @
c          -1, input electron density and zeff, and calculate the              @
c              hydrogen and impurity ion densities for initial time            @
c              only.  after initial time the electron density and              @
c              zeff are calculated from hydrogen and impurity                  @
c              densities.  (icenez = -1)                                       @
c         -99, input electron and impurity densities  JMP START   
c              -> calculate the hydrogen density
c         -98, same as inezez = -99, but calculate He density JMP END
c              at each time step by nHe/taup* = int nD nT <sigv> dV.
c% adjzeff         if inenez = 1 or -1 zeff can lead to negative               @
c                  ion densitites. by setting adjzeff (which is an             @
c                  integer variable) = 1,the code will decrease zeff           @
c                  locally until a value for zeff is found that will           @
c                  lead to positive primary ion and impurity densitites.       @
c                  If zeff attains its minimum value of 1.0 no further         @
c                  adjustment is made.                                         @
c                  adjzeff = 0 (default ) does not make this adjustment.       @
c                  DO NOT USE ADJZEFF IN COMBINATION WITH TWEAKING OF ZEFF!    @
c
c        -- modify input profiles using *_mult settings if
c        -- called for. This is not intended to be a rewrite of the inone
c        -- file! Use this to make simple modifications of the input
c        -- profiles. These profiles changes will not be  visible in your
c        -- inone input file.  But they will affect the output.
c        -- IT IS THE USERS RESPONSIBILITY to use this option wisely 
c        -- It is suggested that appropriate comments be put in the inone file
c        -- for later reference when seemingly  "results not consistent with the profiles
c        -- in inone are encountered"
c        -- guessing at these multiplier values will most often result in
c        -- the "garbage in - garbage out syndrome" . Use these switches responsibly
c        -- by cautiously setting them in inone  and then running a SINGLE STEP
c        -- ONLY to check the outputs.
c        -- the density multipliers ene_mult and density_mult are
c        -- NOT IMPLEMENTED for inenez =-99,-98 at this time.
c        -- None of these multipliers are applied if the input profile is parabolic
c        -- instead of a spline.

c density_mult(1..nprim)       multiplier  for input primary ion  densities (used only if inenez =0)
c density_mult(nprim+1,nion)   multiplier  for input impurity ion  densities (used only if inenez =0)
c                              note that nion = nprim+nimp
c ene_mult                     same for ene (used only if inenez >0)
c te_mult                      same for te   always applied 
c ti_mult                      same for ti   always applied 
c ang_mult                     same for ang rotation speed profiles always applied 
c
c                                                                              
c
c
c
c -------------------------PEDESTAL CODE ---------------------------------------
c use_pedestal   Logical variable, set to .true.  to use the NTCC Pedestal 
c                code ( The Pedestal code uses semiempirical methods,
c                calibrated against Type I Elmy H mode Plasmas).  
c                Models for the L-H transition,
c                pedestal density and pedestal tempeatures are available.
c                A complete description of the models is given in the
c                writeup  of the NTCC module. It is also included in the
c                Onetwo documentation as file pedestal_models.ps
c                The namelist file that the pedestal
c                code uses in stand alone mode,input_1, will be written out
c                by Onetwo based on current parameters.
c                The file will be overwritten each time a new
c                call to Pedestal is made. A user friendy output file
c                output_1 is created together with a machine readable file
c                output_12 tha Onetwo reads.
c
c
c pedestal_path  The default path is set internally to
c         point at curret  version of pedestal. If you  want to run a
c         specific different verion  then you can set the path to that version
c         here. 
c         pedestal_path is a 256 (or less) character variable.
c
c
c pedestal_models(1) :
c          controls which L-H transition model is used
c                  Default pedestal_models(1) = 1
c             = 1  L-H transtion model based on power threadhold
c                  (Other models not available at this time)
c
c pedestal_models(2) :
c          controls which pedestal density model is used
c                  Default pedestal_models(2) = 1
c             = 1  Use simple empirical pedestal density model
c                  (Other models not available at this time)
c
c pedestal_models(3) 
c            controls which pedestal temperature model is used
c                  Default pedestal_models(3) = 1
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
c
c
c The entry into the pedestal code is through 
c subroutine Onetwo_pedestal_driver in cray401.f. The remaining
c pedestal code subroutines are stored in pedestal.f
c
c
c The input to the Pedestal code is listed here for reference only. 
c All quantities are known
c from other input in inone: Onetwo will use the internal values of 
c these parameters at the
c appropriate times to create an input file for the Pedestal code:
c                      rminor     minor radius [m]          from inone
c                      rmajor     major radius [m]          from inone
c                      kappa      elongation                calculated in Onetwo
c                      delta      triangularity             calculated in Onetwo
c                      current    plasma current [MA]       from inone or eqdsk
c                      btor       toroidal field at rmajor [Tesla] eqdsk or inone
c                      nebar      line averaged electron density [particles/m^3] calculated in Onetwo
c                      hydmass    average hydrogenic mass [AMU]            calculated in Onetwo
c                      zeff       effective charge at the edge of plasma   calculated in Onetwo
c                     ploss      power crossing separatrix [MW]            calculated in Onetwo
c                     Note: not all variable are needed for each model
c
c  Output:
c  -------
c The output is used in Onetwo as appropriate. It is listed here for reference only
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
c
c --------------------END  PEDESTAL CODE DESCRIPTION  --------------------------
c
c
c
c
c
c

c  Instead of the standard profiles, spline profiles may be specified.         @
c     These profiles are controlled by the following parameters:               @
c                                                                              @
c% rnormin(j)  normalized radii (from 0 to 1) at which input data              @
c              are specified; these locations are called knots                 @
c% enp(j,i)    density of primary ion species i at knot j                      @
c% njenp(i)    >0, number of knots in spline profile                           @
c              =0, spline profile is not generated.                            @

c% eni(j,i),njeni(i) - same for impurity ion species i.                        @
c% te(j),njte        - same for electron temperature.                          @
c% ti(j),njti        - same for ion temperature.                               @
c% ene(j),njene      - same for electron density.                              @
c% zeff(j),njzef     - same for zeff.                                          @
c% curden(j),njcur   - same for current density.                               @
c  Spline profiles specified by the preceding controls are                     @
c     time-independent.  Time-dependent spline profiles may be                 @
c     specified for the following quantities.  (The spline knot                @
c     values are interpolated linearly in time.):                              @
c           
c       
c  en_bc_inpt(j,i,m) New time dependent boundary condition density profile input.   @
C               This array will currently only be used if gcnmp solver is 
c               selected (diffeq_methd =3).
c               en_bc_inpt(1:nj,1:nion,1:nbctim) is dynamically re allocated after
c               the namelist has ben read. SO in use ,nj,nion,nbctim are the
c               upper limits on this array. Note that no distinction is made
c               on the basis of primary or impurity ions. The convention is
c               that i=1 to nprim for primaries, i=nprim+1 to nion for 
c               impurities. index m is the time index (see below m ==>bctime(m))
c ren_bc_inpt(j,i,m) normalized rho locations
c               The number of knots is automatically determined internally.


c               INPUT RESTRICTIONS  are present for  en_bc_inpt and ren_bc_inpt.
c               That is, if these arrays are used at all then you must
c               specify a complete set. NO DEFAULTS  are allowed!!!
c               This means that nbctim*nion arrays are supplied!!!
c               Each of these arrays contains arbitrary number of spline 
c               points. Each ren_bc_inpt array must start at 0.0 and end at 1.0.
c               Thes earrays are input in rho space so we used zero gradient
c               conditons at  rho=0 and natural spline at rho =1.
c               Setting the number of knots equal to the number of gris points
c               (nj) means that spline interpoaltion will be ineffective but
c                is allowed.
c               
c fix_edge_ni_bc_inpt(1:nion,1:kbctim) fix the edge values,(see fix_edge_te)

c% tein(j,m)  - electron temperature at knot j at time point m.                @
c% njte       - same as defined above.                                         @
c% tiin(j,m),njti   - same for ion temperature.                                @
c% enein(j,m),njene - same for electron density.                               @
c% zeffin(j,m),njzef- same for zeff  
c  zeff_mult is a (constant in time) multiplier for zeffin
c  CAUTION: all non-zero njxx MUST have the same value if splninpt = 'old'     @
c  if splninpt = 'new' this restriction is removed,and curden is also          @
c  allowed to be time-dependent. splninpt = 'old' is retained to allow full    @
c  backwards compatibility. However it is recommended that splninpt = 'new'    @
c  is used for time-dependent spline profiles from here on out.                @
c  default is splninpt = 'old' to be consistent with old input files           @
c **********************************************************************       @
c *** SPLNINPT = OLD IS NO LONGER ACCEPTED *****************************       @
c **********************************************************************       @
c  specifically the following input can be set using different knot sets:      @
c% enp(j,i)  renpin(j,i)  njenp(i)  bparenp(4,i)  (no time dep.)               @
c% eni(j,i)  reniin(j,i)  njeni(i)  bpareni(4,i)  (no time dep.)               @
c% curdenin(j,m) rcurdein(j,m)  njcur       bparcur(4,m)                       @
c% tein(j,m)     rtein(j,m)     njte        bparte(4,m)                        @
c% tiin(j,m)     rtiin(j,m)     njti        bparti(4,m)                        @
c% enein(j,m)    renein(j,m)    njene       bparene(4,m)                       @
c% zeffin(j,m)   rzeffin(j,m)   njzef       bparzeff(4,m)                      @
c% angrotin(j,m) rangrot(j,m)               bparang(4,m)   (see below)         @
c% kpolin(j,m)   rkpol(j,m)                 bparkpol(4,m)  (set ikpol=1)       @
c  tein, tiin, enein, zeffin and curdenin may be input as                      @
c  te  , ti  , ene  , zeff   and curden   respectively,                        @
c  if there is no time dependence.                                             @
c  here index j is the knot and i is the species as usual                      @
c  index m is the time index (see bctime).                                     @
c  rxx specifies the knots (0.0 to 1.0) as rnormin but with                    @
c  possibly different knots for each input profile.                            @
c  If njxx .ne. 0 but the corresponding rjxx is not set,                       @
c  then rnormin will be used for that profile.                                 @
c  the njxx do not have a time index m because they are not used to            @
c  determine the number of knots in rxx(j,m) at time point m. Rather,since     @
c  the first knot must be at zero and the last knot must be at 1.0,the code    @
c  discovers for itself what the number of knots is at any time m > 1.         @
c  however,in order to let the code know that spline (rather than parabolic)   @
c  input profiles are given,each of the njcur,njte,njti,njene,and njzef must   @
c  be set to the initial number of knots at time m = 1.0                       @
c  for the rotation speed profile parabolic input is not allowed so no njxx    @
c  switch is present in the input.                                             @
c                                                                              @
c                                                                              @
c  The bparxx arrays are used to set boundary conditions.                      @
c  (on the spline input,-has nothing to do with boundary conditions on         @
c   transport equations)                                                       @
c  normally the b.c. at rho = 0 is zero gradient so bparxx(1) and bparxx(2)    @
c  do not have to be set since zero gradient at rho=0 is the default.          @
c  however to set the gradient at rho = 0 to the value x,specify               @
c         bparxx(1) = -1.0e30,bparxx(2) = x                                    @
c  at rho = 1.0 the default boundary condition is the natural spline           @
c  to change it to a gradient condition,say gradient = y, proceed as above:    @
c         bparxx(3) = -1.0e30,bparxx(4) = y                                    @
c  NOTE: at present only a gradient b.c. is allowed as input.                  @
c  when using the splninpt='new' option,you may,but do not have to             @
c  specify profiles for all times nbctim. Any profile and knot set             @
c  which is required but is not given in the input will be copied              @
c  from the last previous time point at which it is available. Hence           @
c  profiles which do not change (and the corresponding knot set ) need         @
c  only be given onece. extensive error checking is done                       @
c  and a summary print of the input spline profile parameters is given         @
c  in the outone file.                                                         @
c  because there is no easy way to default or check the bparxx arrays,         @
c  these arrays must be given explicitly at all times m! The only              @
c  exception is if bparxx(j,m) is zero for all m and j = 1...4  In this        @
c  case the bparxx array need not be specified at all.                         @
c                                                                              @
c% angrotin(j,m)    angular rotation speed,rad/sec,at knot rangrot(j,m)        @
c                   at time bctime(m) of the bulk ions. If angrotin is time-independent but     @
c                   nbctim .gt. 1 then the values of angrotin for other times  @
c                   do not have to be supplied(see rangrot for further info)   @
c% rangrot(j,m)      similar to rnormin but for angular rotation speed         @
c                    profile. rangrot(j,m) gives the j'th knot                 @
c                    at time bctime(m). Since the first knot must              @
c                    always be at zero and the last knot must always           @
c                    be at one the code automatically determines the number    @
c                    of knots input so no variable specifying the number       @
c                    of knots at time bctime(m) is required.                   @
c                    the range on j is [3,ksplin]                              @
c                    if the angular rotation speed is to be time-independent   @
c                    then simply specify angrotin(j,m) with m = 1. for m = 2.. @
c                    nbctim leave angrotin(j,m) unspecified (i.e., 0.0).       @
c                    the knot locations and number of knots at each            @
c                    time point m are not required to be the same. if they     @
c                    are the same however the user need only give the values   @
c                    for m = 1 in rangrot(j,m). More generally if angrotin     @
c                    and/or rangrot is the same for all times greater than     @
c                    the last specified time the values do not have to be      @
c                    repeated.                                                 @
c                    there is no plan at present to allow input of rotation    @
c                    speed profiles in terms of the parabolic parameter set.   @
c% bparang(4,kbctim)     boundary condition array for angrotin.                @
c  note that angrotin, bparang and rangrot are independent of the setting      @
c  of splninpt (i.e., works with either splninpt = 'old' or 'new')             @
c  (there is no backward compatibility problem here).                          @
c                                                                              @
c% kpolin(j,m)        poloidal velocity of cer ion/Bp at knots. This is only   @
c                     used to compute the radial electric field from the CER   @
c                     data if cer_ion is set. Output of Epsi_exp and Kpol_exp  @
c                     is to the netCDF file iterdb.nc (must set iterdb = -1)   @
c                                                                              @
c% rkpol(j,m)         knot location for kpolin similar to rangrot above        @
c                                                                              @
c% bparkpol(4,kbctim) boundary conditions for splines similar to bparang above @
c                                                                              @
c  Time-dependent profiles and boundary conditions are controlled by           @
c     two parameters:                                                          @
c                                                                              @
c% nbctim      Number of times at which boundary conditions                    @
c              are input (up to kbctim).  See note below.                      @
c% bctime(m)   Times (sec) at which boundary conditions are input              @
c                                                                              @
c  All preceding quantities with subscript m,                                  @
c  plus elong(m), gasflx(m,i), w*typ(m), w*pro(j,m)  and qradin(j,m),          @
c  are given at the time values in bctime.                                     @
c  If nbctim = 1, all quantities are time-independent.                         @
c  If nbctim>1, the second value of each quantity specifies whether            @
c  it is time-dependent.  If the second value is 0 (default), the              @
c  quantity is time-independent, and the first value is used                   @
c  at all times.  If the second value is not zero, it is assumed               @
c  that nbctim values have been given for the quantity.                        @
c  (Use a very small number if the second time-dependent quantity is 0)        @
c                                                                              @
      adjzeff = 0                                                              %
      zeff_mult =1.0
      pfact   = 1.0                                                            %
      vloop_thresh= HUGE(1.d0)
      vloopvb = 0
      density_mult(:) = 1.0D0
      te_mult         = 1.0D0
      ti_mult         = 1.0D0
      ene_mult        = 1.0D0
      ang_mult        = 1.0D0
      zeff_mult       = 1.0D0
      he3_frac        = 0.0D0
c                             
      do i=1,kbctim                                                            %
        if (i .gt. 1)  pfact = 0.0    ! else logic in SPECIFY won't work       %
        do k=1,kprim                                                           %
          enpc(i,k)   = pfact * 5.0e13                                         %
          enpb(i,k)   = pfact * 1.0e13                                         %
          alpenp(i,k) = pfact * 1.0                                            %
          gamenp(i,k) = pfact * 2.0                                            %
        end do                                                                 %
        do k=1,kimp                                                            %
          enic(i,k)   = pfact *  5.0e10                                        %
          enib(i,k)   = 0.0                                                    %
          alpeni(i,k) = pfact * 1.0                                            %
          gameni(i,k) = pfact * 2.0                                            %
        end do                                                                 %
        tec(i)    = pfact * 1.0                                                %
        teb(i)    = pfact * 0.2                                                %
        alpte(i)  = pfact * 1.0                                                %
        gamte(i)  = pfact * 2.0                                                %
        tic(i)    = pfact * 0.4                                                %
        tib(i)    = pfact * 0.2                                                %
        alpti(i)  = pfact * 1.0                                                %
        gamti(i)  = pfact * 2.0                                                %
        xjc(i)    = pfact * 1.0                                                %
        xjb(i)    = pfact * 0.1                                                %
        alpxj(i)  = pfact * 1.0                                                %
        gamxj(i)  = pfact * 2.0                                                %
        enec(i)   = pfact * 5.0e13                                             %
        eneb(i)   = pfact * 1.0e13                                             %
        alpene(i) = pfact * 1.0                                                %
        gamene(i) = pfact * 2.0                                                %
        zeffc(i)  = pfact * 1.0                                                %
        zeffb(i)  = pfact * 1.0                                                %
        alpzef(i) = pfact * 1.0                                                %
        gamzef(i) = pfact * 2.0                                                %
        totcur(i) = pfact * 7.0e5                                              %
        vloop_bc(i) = vloop_thresh   ! any value less  than this will kick in
                                     ! vloop boundary condition
        vloop_bc_time(i) = vloop_thresh
      end do    
      use_pedestal = .false.
      pedestal_path =''              !default must be zero length string
      fix_edge_ni_bc_inpt(:,:) = 1.0  !default all values to plasma edge

!      IF(.NOT. ALLOCATED(en_bc_inpt)) 
!     .                  ALLOCATE(en_bc_inpt(kj,kion,kbctim))
!      IF(.NOT. ALLOCATED(ren_bc_inpt)) 
!     .                  ALLOCATE(ren_bc_inpt(kj,kion,kbctim))
!      IF(.NOT. ALLOCATED(ken_bc_inpt)) 
!     .                  ALLOCATE(ken_bc_inpt(kion,kbctim))
      ren_bc_inpt(:,:,:) = -1. ; en_bc_inpt(:,:,:) = -1. 
      ken_bc_inpt(:,:) = -1.
c                  
      bp0_ic ='eqdsk'                                      
      nbctim   = 1                                                             %
      inenez   = 1                                                             %
      do 1830 i=1,kprim                                                        %
 1830 njenp(i) = 0                                                             %
      do 1840 i=1,kimp                                                         %
 1840 njeni(i) = 0                                                             %
      njte     = 0                                                             %
      njti     = 0                                                             %
      njene    = 0                                                             %
      njzef    = 0                                                             %
      njcur    = 0                                                             %
      splninpt = 'old'   ! don't set to new (will defeat error checking)       %
      call zeroa (tein,ksplin*kbctim)                                          %
      call zeroa (curbeam_external,ksplin*kbctim)                              %
      call zeroa (tiin,ksplin*kbctim)                                          %
      call zeroa (enein,ksplin*kbctim)                                         %
      call zeroa (zeffin,ksplin*kbctim)                                        %
      call zeroa (angrotin,ksplin*kbctim)                                      %
      call zeroa (kpolin,ksplin*kbctim)                                        %
      call zeroa (bparkpol,4*kbctim)                                           %
      call zeroa (bparang,4*kbctim)                                            %
      call zeroa (bparenp,4*kprim)                                             %
      call zeroa (bpareni,4*kimp)                                              %
      call zeroa (bparte,4*kbctim)                                             %
      call zeroa (bparti,4*kbctim)                                             %
      call zeroa (bparcur,4*kbctim)                                            %
      call zeroa (bparzeff,4*kbctim)                                           %
      call zeroa (bparene,4*kbctim)                                            %
      call zeroa (bparcurb_external,4*kbctim)                                  %

c ----------------------------------------------------------------------       @
c --- SOLUTION CONTROL PARAMETERS                                              @
c ----------------------------------------------------------------------       @
c% time0     Initial time (sec)                                                @
c% timmax    Maximum time (sec)  set timmax .le. time0 for snapshot            @
c              mode. note that in snapshot mode all dependent variables        @
c              must be run in analysis mode and it is implicitly               @
c              assumed that all time derivatives are zero. At present          @
c              specifying profiles as a function of time will not yield        @
c              the desired time derivatives in the snapshot mode.              @
c              also note that energy source terms due to the rate of change    @
c              of beam or fusion ion densities must be neglected since         @
c              these rates are not available at the initial time.              @
c              Due to roundoff it is best if timmax is set slightly less       @
c              (say one millisec ) than bctime(nbctim) if multiple             @
c              boundary  values are specified on input.                        @
c  timecrit    The intent of this value is to allow larger time steps after
c              the plasma settles down. This is accomplished by setting
c              dtmax = dtmaxcrit for all times > timecrit.
c              Note that timecrit must be .gt. time0 
c  dtmaxcrit   see timecrit
c% nmax      Maximum number of time steps                                      @
c% dt        Initial time step (sec)                                           @
c% dtmin     Minimum time step (sec)  (suggest 1.0e-6)                         @
c% dtevmin   Minimum time step for plotting, printing, beam on (sec)           @
c% dtmax     Maximum time step (sec) value depends on confinement models       @
c            used in the simulations. For stiff models a small value << 1 msec @
c            may be necessary. Also may be changed using dtmaxcrit             @
c% relmin    Minimum delmax, where delmax is the maximum relative              @
c            change from one time point to the next in any dependent           @
c            variable.  If delmax < relmin, the time step is doubled.          @
c            (suggest 0.10)                                                    @
c% relmax    Maximum delmax.  If delmax > relmax, the time step is             @
c            halved.  (suggest 0.30)                                           @
c            Note: both simulation and analysis mode profiles are              @
c            checked. As a consequence the time step will be cut               @
c            if the change in a profile is greater than relmax,even if         @
c            that profile is not run in simulation mode. In particular if      @
c            you run in total analysis mode (i.e., all itran(i) =0) then       @
c            no profiles are transported but the change in                     @
c            input profiles and beams may still cut down the time step.       @
c            The neutral density equation has to be solved even in this        @
c            case and may also result in a decreased time step.                @
c            Set relmax to a large number to avoid this as appropriate.        @
c% relit     Maximum delit, where delit is the maximum relative                @
c            difference between the predictor and corrector in any             @
c            dependent variable.  If delit > relit, the corrector              @
c            calculation is iterated.  (suggest 0.01)                          @
c            set relit = a large number and dtmax a small time step            @
c            to mimick Euler type solution method HSJ                          @
c% itmax     Maximum number of corrector iterations  (suggest 10)              @
c            (Note: includes particle conservation iterations.)                @
c% theta     Time weighting of diffusion term in difference                    @
c            equations:                                                        @
c            0.0, fully explicit (unstable)                                    @
c            0.5, Crank-Nicolson (marginally stable)                           @
c            1.0, fully implicit (stable)                                      @
c            (suggest 0.8)                                                     @
c% ilimdt    1, allow restriction of dtmax to numerically stable               @
c            regime, as determined by the routine DTLIM (default).             @
c            0, disable this: dtmax is not touched.                            @
c             The standard Euler approach can be mimicked by setting           @
c             ilimdt=1, relit = a large number,dtmax a small time step         @
c             (small means relative to the rate at which the sources are       @
c              forcing the solution to change,typically .1 to 10  ms)  HSJ     @
c% timav     >0, time average solution over the interval timav;                @
c                output time-average results at last time point                @
c             0, no effect  (default)                                          @
c            <0, time average solution over the interval timav;                @
c                output last time results in addition to time-average          @
c                results                                                       @
c% iffar     1, Fast Faraday's law:  speeds up the current diffusion in        @
c               Faraday's law.  Useful when steady state is desired.           @
c               This option is invoked automatically when 'tweaking'           @
c               with voltin or qcin non-zero.                                  @
c            0  no effect  (default)                                           @
c                                                                              @
c% relaxsol        relaxation parameter for the solution during                @
c                  the corrector iterations. set relaxsol = 1.0 to             @
c                  disable the relaxation. otherwise relaxsol should           @
c                  be set to a number between 0 and 1. to underelax            @
c                  the solution,if it oscillates due to a nonlinearity         @
c                  in the transport coefficients. relaxsol only                @
c                  needs to be used for models such as Rebut-Lallia            @
c                  when running te in simulation mode to prevent               @
c                  extremely short time steps.                                 @
c% steady_state   Multiplier on time derivative terms for all simulation       @
c%                eqautions. Normally equal to 1.0 but can be set to a small   @
c%                number to try and approach steady state rapidly. Used only   @
c                 if diffeq_method = 2    !!!!!!!!!!!!!!!!!!!!!                @
c                 Use nmax to control the number of                            @
c                 iterations allowed . If steady_state =0.0 then the time
c                 derivative in all simulation equations is set to 0.0.
c
c
c                  The choices available for solving the coupled set of
c                  parabolic equations are:
c     diffeq_methd  =0 default,predictor-corrector ,Crank-Nichelson solution   @
c                      method. 
c                      (does not work well with stiff models like GLF23)       @
c                   =1 method of lines for solving systems of ODE's            @
c                      There are four  ODE integrators implemented:
c                        1) Radau5  (works fine, select with non_lin_method  
c                        2) Dlsoibt (works fine, select with non_lin_method
c                        3) Dlsodi (needs some additional coding )
c                        4) Dgear  (only for te,ti simulation) select with
c                            non_lin_method 
c                         Dlsodi and Dlsoibt are from odepack (Public Domain)  @
c                         (they are modern versions of Dgear)
c                         Radau5 is from Hairer and Wanner "Solving ODEs II
c                         Stiff and Differential Algebraic Problems", Springer,1996
c                         Available at 
c                         http://www.unige.ch/math/folks/hairer/software.html
c                   =2 non linear solver (newton type approach)                @ 
c                      a sequeqnce of algorithms (developed by HSJ for Cerfit originally), 
c                      based on work by Dennis, Schnabel and others.
c                      This implementation seems to work much better than any
c                      I tried from minpack and various other sources.                       
c                      Some problems work best with this solver, others        @
c                      are more efficiently solved with one of the ODE solvers @
c                      No general guide lines are availble at this time.       @
c                      But:                                                    @
c                      The ODE solvers pick their own time step based          @
c                      on error minimization methods and this may be           @
c                      advantageous on occasion. However it also makes these   @
c                      methods quite sensitive to the input tolerances!        @
c                      Newton solver can give the the best accuracy            @
c                      but may take excessive time to run.                     @
c                      The predictor-corrector method normally requires additional @
c                      under ralaxation and/or averaging to work wiht stiff models.@
c                      (the ODe solves also use some form of newton solver     @
c                        internally)                                           @
c                      For simple problems use predictor-corrector,            @
c                      for other cases try the newton and method of lines solvers.
c                      For any given case one method is sure to work best.
c                  =3  Use gcnmp for solving equations, same as 2 but uses 
c                      GCNMP instead. Reading and writing statefiles is implied here.
c                      The intialization can be through standard inone files
c                      or trhough a previosuly created statefile, see
c                      initialize_from_statefile. 
c                      NOTE: There are some important input variables 
c                            introduced in gcnmp that have no counterpart 
c                            in Onetwo: 
c    itenp,iteni,itte,itti,itxj : when run with diffeq_methd  =3  these 
c                                 variables can be -1,0,1
c                                 (recall that when run in Onetwo the 
c                                  value -1 is NOT allowed.)
c                                 The +-1 settings of these variables in gcnmp
c                                 allows you to subdivide equations run in
c                                 simulation mode into two disparate sets.
c                                 Each set will be solved simultaneously but
c                                 on interleaving time scales. For example,
c                                 all equations with +1 settings will be run from
c                                 t to t+dt, while all equations with -1 will
c                                 be run from t+dt/2 to t+3*dt/2 using fixed values
c                                 for the +1 set evaluated at t+dt.
c    save_incr_restart_file_gcnmp,INTEGER, if diffeq_method = 3 then this value can be
c                                used to save intermediate results from gcnmp 
c                                 curently this option is not active.
c    gcnmp_macro_dt              Lenght of time,sec for which gcnmp is run on
c                                each call from Onetwo
c    gcnmp_host                  (character <= 128) name of host 
c                                 on which  gcnmp will be run. Valid hosts are ones 
c                                 where  gcnmp_server is running. [If in doubt use
c                                 the command
c                                 'python gcnmp_client.py gcnmp_nprocs inputstatefile gcnmp_nml_filename  remote_dir local_dir gcnmp_dir lohan5 '
c                                 to check if server is running on lohan5 for example
c                                 before you run Onetwo with gcnmp. if you get 
c                                   error,connection on  lohan5  refused then the server
c                                   is not running]
c                                  If server is not running you may be able
c                                 to start it using 'python gcnmp_server.py'
c                                 (on the server machine,lohan5 in the above example).
c                                 Normally Valid names are
c                                 lohan4-6, benten,  because we use 64 bit gcnmp.
c                                 Also see gcnmp_client.py for current list of valid hosts.
c    gcnmp_nprocs                  # processors to use (GCNMP is mpi based, 1 or more is valid)
c    gcnmp_nml_filename          (character <= 128) name of gcnmp namelist input file
c
c    write_iterdb_txt            LOGICAL, set to true to write text file gcnmp iterdb
c                                FALSE means write netcdf file
c    switch_iterdb_output        Set to 0 or 1. 0 means gcnmp output will be in same form
c                                (eg either netcdf or text) as input . 1 means switch 
c                                  output to text if input was in netcdf form or switch 
c                                 output of gcnmp to netcdf if input was in text form.
c                                 (This option is primarily to aid in development
c                                  and testing) - Can  be used to translate iterdb files
c    gcnmp_iterdb_filename       (CHARACTER) Name of statefile to be created by onetwo 
c                                and read by gcnmp
c
c
c non_lin_method       solution method used if diffeq_method = 1 or 2          @
c                       set as follows:
c                        if diffeq_methd =1 then
c                            non_lin_method     =1 selects Radau5
c                            non_lin_method     =2 selects Dlsoibt
c                            non_lin_method     =3 selects Dgear
c                            non_lin_method     =4 selects Dlsodi (dont use for now)
c                        if diffeq_methd = 2 then
c                            non_lin_method     =1 selects linesearch
c                            non_lin_method     =2 selects hookstep
c                            non_lin_method     =3 selects dogleg
c                            Normally for diffeq_method =2 a round robin method
c                            that cycles through the above three methods is used and
c                            non_lin_method only determines the initial
c                            method used for the first maxfev iterations.
c                            The method will be cycled  if maxfev iterations
c                            are done before convergence is achieved.
c                            To stick with a single mehod, as selected by non_lin_method,
c                            set maxfev .eq. tot_iters_max
c                  
c     
c wrt_nwt                 set to 1 to write the file nwt.test
c                         implemented only for diffeq_methd = 2 
c
c                                                                              @
c             IF diffeq_methd = 2 the follwoing input is required
c  ssqrmin        minimum sum of squares of residuals of all equations
c                 convergence is assumed if the value is less than ssqrmin
c  steptol     convergence is assumed if step size calcualted is less than step tol
c  gradtol     convergence is assumed if gradient of sum of squares function is
c              less than this value
c  fvectol     ssq of any single equation 
c  maxfev      the max number of iterations allowed for a single solution
c              method as selected by non_lin_method. If convergence is
c               not achieved before maxev iteratiosn are performed the
c              soltution method is cyclically changed to the the
c              next method in a round robin fashion
c  bandwith    this is in addition to the bandwith 2*nkt -1 used internally
c              it is experimental and tests the sensititvity to splines,e etc.
c               that couple grid points further away than the finite difference 
c              approximations.
c 
c jac_skip     the jacobian will be calculated fresh each time
c              Mod(iters,jac_skip) = 0. 
c              for time dependent problems if the time step is small
c              enough it may be possible to save execution time by 
c              setting jac_skip to a number greater than 1.
c              most likely jac_skip =1,2,3, will be most useful but
c              you can try larger numbers. note that jac_skip =1 means 
c              calculate the jacobian on every iteration.
c fdigits      number of good digits in fucntion evaluations.
c jacobian_type   0 full jacobian, full calcs
c                  1 full jacobian, factored secant, Broyden update
c                 2 sparse jacobian, full calcs
c                 3 sparse jacobian , factored secant, Broyden update
c               set to 1 to use full, factored secant jacobian
c switch_method integer, gives number of times we can switch methods  in round
c               robin style from dogleg to hookstep to linesearch. set
c               to 1 to prevent switching methods.  Each method will be run
c               for a maximum of maxfev times (unless a converged result
c               terminates the process).
c tot_iters_max maximum number of passes through the nonlinear solver allowed
c iters_freeze  for time dependent run freeze the models according to 
c               freeze_type after iters_freeze iterations have been taken
c               without satisfactory convergence.
c freeze_type   integer allows freezing of profiles during non linear
c               iterations as follows. The profiles or anomalous part of the
c               diffusion coefficient are frozen at their current values for the
c               remainder of the current time step.
c               =- 5 freeze te,ti,wr
c                 -4        te,ti
c                 -3        te
c                 -2        ti
c                 -1        wr
c                  0        dont freeze anything
c                  1  freeze (anomalous part of) wr chi
c                  2                             ti
c                  3                             te
c                  4                             te and ti
c                  5                             te and ti and wr
c (wr is toroidal rotation)
c  freeze_xsn     set equal to 1 to enable freezing of xnue(in qdelt) and
c                 explicit source terms in sub source. This option
c                 is implemented only for time dependent runs.
c                 At the start of each time step freeze_xnue and
c                 freeze_source are unfrozen. After (5?) iterations
c                 freeeze_source and freeze_xnue may be turned on if
c                 freeze_xsn = 1. 
c  conv_skip      A CYCLE  consists of maxfev iterations ,unless
c                 convergence is achieved before maxfev iterations
c                 are done. In that case conv_skip is ignored for this CYCLE.
c                 Otherwise if conv_skip = 1 then 
c                 skip convergence test at the end of each CYCLE.
c                 That is, after performing a maximum of tot_iters_max
c                 iterations go on to the next time step even if
c                 convergence conditions were not meet.
c                 ( done in solve_newton, and applies only if diffeq_methd = 2)

c                 
c 
c  te_range_check logical,default true, set to false to disable range
c                 checking in sub zfit and asqfit. (if te goes out of
c                 range in  then the largest possible  value will be used )
c                 NOTE; tis is not for radfit !
c
c
c  write_profiles logical,creates output of densitites,temperatures,
c                 current density and rotation in form suitable for use in
c                 inone at the final time.
c ----------------------------------------------------------------------       @
c                                                                                 
      write_profiles =.false.
      wrt_nwt = 0
      steptol =0.0d0
      gradtol = 0.0d0
      fvectol =0.0d0
      tot_iters_max = 0
      iters_freeze = 20         !set to > tot_iters_max in inone  to turn off
      maxfev = 0
      bandwidth  = 0
      fdigits =0 
      jacobian_type = 0 
      jac_skip = 1
      switch_method = 0
      steady_state = 1.0   
      freeze_type = 0                                                          %
      time0      =   0.0                                                       %
      timmax     = 200.0e-3                                                    %
      timecrit   =   1.0d+39                                                   %
      nmax       = 200                                                         %
      dt         =   1.0e-3                                                    %
      dtmin      =   1.0e-6                                                    %
      dtevmin    = dtmin                                                       %
      dtmax      = 100.0e-3                                                    %
      dtmaxcrit  = dtmax                                                       %
      relmin     =   0.10                                                      %
      relmax     =   0.30                                                      %
      relit      =   0.01                                                      %
      theta      =   0.8                                                       %
      itmax      =  10                                                         %
      ilimdt     =   1                                                         %
      timav      =   0.0                                                       %
      iffar      =   0                                                         %
      idterr     =   0                                                         %
      ineucg     =   0                                                         %
      ifreya     =   0                                                         %
      ifreya_old =   0                                                         %
      do k=1,krf                                                               %
        irfcnt(k) = 0                                                          %
      end do                                                                   %
      diffeq_methd = 0                                                         %
      relaxsol = 1.0                                                           %
      te_range_check = .true.
      gcnmp_host = 'localhost'
      gcnmp_nml_filename ='gcnmp_nml_input'
      write_iterdb_txt = .FALSE.
      gcnmp_nprocs = 1
      switch_iterdb_output = 0
      gcnmp_macro_dt = 0.005
      save_incr_restart_file_gcnmp = 0
      gcnmp_remote_dir = '/usr/tmp'
      statefile_name = 'iterdb_12_statefile'
      statefile_type = -1     ! -1 ==> type unknown,0 ==> text ,==>1 netcdf
      gcnmp_iterdb_filename = 'statefile'
c ----------------------------------------------------------------------       @
c --- OUTPUT PARAMETERS                                                        @
c ----------------------------------------------------------------------       @
c% timprt   Time interval (sec) between printout                               @
c% timplt   Time interval (sec) between 3-D plots                              @
c% mprt     Number of time steps between printout                              @
c% mplot    Number of time steps between 3-D plots                             @
c% jterow   List of grid points for te plot -                                  @
c           plot te(jterow(i)) versus time                                     @
c           for i = 1 to nterow                                                @
c% nterow   See above comment                                                  @
c% prtlst   List of times (sec) for printout  (nprtlst_max values)             @
c           prtlst also controls times at which the "analysis mode" variable   @
c           parameters are plotted                                             @
c           These are the sequence of plots thermal conductivites,             @
c           temperatures,densities,confinement times ...etc that follow the    @
c           3D plots in the trplot.cgm file output. THE ACTUAL TIME            @
C           PLOTTED IS THE CORRECTOR  PIVOT POINT TIME                         @
c           GIVEN BY PRTLST(I)-(1-THETA)*(CURRENT TIME STEP)                   @
c           (see desription of theta ) HENCE THE TIMES IN THE PLOTS WON'T      @
c           MATCH THE PRTLST TIMES EXACTLY. (LIMIT THE MAXIMUM CURRENT TIME    @
c           STEP WITH DTMAX IF THIS IS A PROBLEM FOR YOU)                      @
c                                                                              @
c% pltlst   List of times (sec) for 3-D plots (pltlst_max values)                       @
c            Also see above prtlst list!                                       @
c            ONE and only ONE  3d plot will be made for each                   @
c            of the parameters which is designated as a 3D quantity.           @
c            pltlst guarantees that the profile as a function of rho at the    @
c            times pltlst(i) will be part of this single 3d plot, however      @
c            the dimensions of the arrays may conspire to defeat this command) @
c            Since info for the 3d plots is generated at the end of each       @
c            time step this option really only makes sense if there are        @
c            too few time steps in your analysis. However the number of time   @
c            steps is best controlled using dtmax so, in my (HSJ) opinion,     @
c            this switch is rather useless. We may remove it in the future.    @
c            it is here only for backward compatibility. -- HSJ)               @
c                                                                              @
c% jprt     Number of mesh intervals between printout                          @
c% jflux    0, do not print out components of fluxes                           @
c           1, do     print out components of fluxes                           @
c% jcoef    0, do not print out transport coefficients                         @
c           1, do     print out transport coefficients                         @
c% jsourc   0, do not print out comps. of particle and energy sources          @
c           1, do     print out comps. of particle and energy sources          @
c% jbal     0, do not print out balance tables                                 @
c           1, do     print out balance tables                                 @
c% jtfus    0, do not print out detailed fusion information                    @
c           1, do     print out detailed fusion information                    @
c                     if ifus=1,jtfus=1 is forced.                             @
c           JSXR IS OUT OF DATE (OBSOLETE, USES DOUBLET III STUFF              @
c                               don't use this option blindly) HSJ             @
c% jsxr     0, do not calculate and print SXR info                             @
c           1, do     calculate and print SXR info for old Doublet III         @
c                     diode arrays                                             @
c           2, do     calculate and print SXR info for new Doublet III         @
c                     diode arrays                                             @
c          LEAVE JSXR SET AT 0 UNTIL FURTHER NOTICE.  HSJ                      @
c% jco2     0, do not calculate and print CO2 info                             @
c           1, do     calculate and print CO2 info for Doublet III             @
c                     (line-averaged electron density for various              @
c                     tangency radii)                                          @
c% jzeff    0, do not calculate and print ZEFF info                            @
c           1, do     calculate and print ZEFF info for Doublet III            @
c                     (visible continuum light emission for various            @
c                     tangency radii)                                          @
c          LEAVE JZEFF SET AT 0 UNTIL FURTHER NOTICE.  HSJ                     @
c% spec_profiles(j)   j=1,...up to nj                                          @
c                     spec_profiles(j)=k is used to print out the              @
c                     values te(k) and ti(k) at radial grid point k            @
c                     as a function of time. (note that this                   @
c                     is at constant zeta, not necessarily at constant rho)    @
c% steps_per_plot     integer switch,default value of 1,that determines the    @
c                     number of time steps between calls to the plotting       @
c                     routine trplot. By setting steps_per_plot to an          @
c                     integer greater than 1 the output plot file can be       @
c                     reduced in size if coarser plots are acceptable.         @
c% do_eqplot          the equilibrium plot file,eqpltfil, can get very large.  @
c                     set do_eqplot = 0 to truncate this file.                 @
c                     default =1 , all the equilibrium info is written out.    @
c                     eqpltfil is plotted with program eqplot. Because         @
c                     eqplot uses display it is available only on hydra.       @
c                     /u/stjohn/onetwo_v3/eqplot/hp/eqplot_129_129             @
c                     note that eqplot is unfinished , plots are not labeled,  @
c                     etc. (it is intended to move eqplot to an open source    @
c                     plotting package - hence I will not update it)           @
c
c
c   save_scratch1     if zero delete scratch1 file after run is completed
c                     if = 1 then save this file on disk.
c                     file scratch 1 contains some info about flux contours
c                     originally required by some rf codes.
c                     It is now used to process contours outside of Onetwo.
c                     Program read_scratch can be used to read this information.
c ----------------------------------------------------------------------       @
c                                                                              @
      call zero_intg (spec_profiles, kj)                                       %
      timprt = 1.0e20                                                           %
      mprt   = 0                                                               %
      jprt   = 5                                                               %
      timplt = 0.050                                                           %
      mplot  = 0                                                               %
      nterow = 0                                                               %
      jflux  = 0                                                               %
      jcoef  = 0                                                               %
      jsourc = 1                                                               %
      jbal   = 0                                                               %
      jtfus  = 0                                                               %
      jsxr   = 0                                                               %
      jco2   = 0                                                               %
      jzeff  = 0                                                               %
      steps_per_plot = 1                                                       %
      do_eqplot = 1                                                            
      prtlst(:) = HUGE(1.D0)
      pltlst(:) = HUGE(1.D0)
      save_scratch1 = 0
c ----------------------------------------------------------------------       @
c --- NEUTRAL TRANSPORT PARAMETERS                                             @
c ----------------------------------------------------------------------       @
c% nneu        Number of neutral species (0, 1 or 2)                           @
c% namen       Names  of neutral species:  For these primary species,          @
c              neutral transport and radiative recombination of ions           @
c              will be modeled in the code.  (Note that He is treated          @
c              as if there were no singly-ionized state.)                      @
c              Permitted for primary species #1 and #2, only.                  @
c              FOR TWO NEUTRAL SPECIES INPUT THE NEUTRALS IN THE SAME ORDER    @
c              AS THE PRIMARY IONS                  HSJ                        @
c     NOTE     The following parameters that control neutral profiles          @
c              take effect only if the primary ions are run in                 @
c              simulation mode. For analysis mode use taupin to                @
c              get some control.                                               @
c% gasflx(m,i) Flux of injected neutral gas for primary species i.             @
c              (1/cm**2-sec).  If time-dependent, specify nbctim and           @
c              bctime(m).                                                      @
c% ipcons(i)   0, recycling of species i (default) with coeff recyc (iref = 1) @
c              1, adjust recycling of species i to maintain a constant         @
c                 number of ions in the plasma.                                @
c% recyc(i)    recycling coefficient of species i. (default = 1.0)             @
c% twall       Temperature (keV) of neutrals emitted from wall                 @
c% wion        Electron energy loss (keV) per ionization due to                @
c              ionization                                                      @
c% wrad        Electron energy loss (keV) per ionization due to                @
c              radiation                                                       @
c% raneut      Effective plasma radius (cm) seen by neutrals;                  @
c              if input as zero, raneut is calculated internally as            @
c              SQRT (elong(1))*rminor   (elong may be time-dependent).         @
c% relneu      Minimum relative change in ion density or 1/4 of ion            @
c              temperature between neutral transport calculations              @
c              (suggest 0.10)                                                  @
c% erneumax    Maximum epsneu, where epsneu is the relative error in           @
c              the total number of primary ions in the plasma.  If             @
c              epsneu > erneumax, the neutral density profile is               @
c              adjusted and the corrector calculation is iterated.             @
c              (suggest 0.005)                                                 @
c                                                                              @
c% widths      Width of plasma scrape-off layer surrounding the                @
c              main plasma (cm).  Default = 0.0, no scrape-off layer.          @
c              Transport is not modeled in this layer: it is used              @
c              to obtain more realistic neutral densities within the           @
c              main plasma.                                                    @
c% nps         Number of points in the scrape-off layer (suggest 4).           @
c              These are added to the nj points in the main plasma.            @
c% enes        Electron density in the scrape-off layer (1/cm**3).             @
c              Default = 0.0 causes eneb to be used.  The charge               @
c              state of impurities and proportions of primary ion              @
c              species used in the layer are the same as those at              @
c              the edge of the main plasma.                                    @
c% tes         Electron and ion temperature in the scrape-off                  @
c              layer (keV).  Default = 0.0 causes teb and tib                  @
c              to be used.                                                     @
c                                                                              @
c% idiagn      Neutral diagnostic flag                                         @
c              0: do not calculate neutral diagnostic quantities               @
c              1: calculate specific flux of outgoing neutrals -               @
c                 rtandn,rhdn define the detector line-of-sight chord          @
c% rtandn      Tangency chord of detector chord in major plane (cm).           @
c              Default = 0.0, chord is radial.                                 @
c% rhdn        Height above midplane of detector chord (cm).                   @
c              Default = 0.0, chord is in midplane.                            @
c% nengn       Number of energies at which spectrum is calculated              @
c              (max = 50)                                                      @
c% englstn     List of energies (keV) for neutral spectrum                     @
c% neucgvb     verbose switch for neucg. set to 1 to get messages to screen    @
c ----------------------------------------------------------------------       @
c                                                                              @
      nneu     = 0                                                             %
      namen(1) = 'h'                                                           %
      namen(2) = ' '                                                           %

      erneumax = 0.005                                                         %
      do i=1,2                                                                 %
        ipcons(i) = 0                                                          %
        recyc (i) = 1.0                                                        %
        do m=1,kbctim                                                          %
          gasflx(m,i) = 0.0                                                    %
        end do                                                                 %
      end do                                                                   %
      twall  = 0.010                                                           %
      wion   = 0.014                                                           %
      wrad   = 0.006                                                           %
      raneut = 0.0                                                             %
      relneu = 0.10                                                            %
      widths = 0.0                                                             %
      nps    = 4                                                               %
      enes   = 0.0                                                             %
      tes    = 0.0                                                             %
      idiagn = 0                                                               %
      nengn  = 0                                                               %
      rtandn = 0.0                                                             %
      rhdn   = 0.0   
      iperp  = .false.                                                         %
c                                                                              %
c ----------------------------------------------------------------------       %
c *** debug information ***                                                    %
c ----------------------------------------------------------------------       %
c                                                                              %
      do i=1,20                                                                %
        xdebug(i) = 0.0                                                        %
      end do                                                                   %
c                                                                              @
c ----------------------------------------------------------------------       @
c% ifred:  1 turn on Fred Marcus output;  0 turn off Fred Marcus output        @
c ----------------------------------------------------------------------       @
c                                                                              @
      ifred = 0                                                                %
*                                                                              @
c **********************************************************************       @
c *** SECOND NAMELIST (NAMELIS2) ***************************************       @
c **********************************************************************       @
c                                                                              @
c ----------------------------------------------------------------------       @
c                     YOKA ANALYSIS RUN PARAMETERS                             @
c    (input in namelist 2 or 3)
c ----------------------------------------------------------------------       @
c% iyoka       0, Not a Yoka analysis run (default)                            @
c              1, Yoka analysis run: compute and print out quantities          @
c                 to be stored in the Yoka Neutral Beam data bank.             @
c                 (Sets jsourc = 1 to compute needed quantities.)              @
c              2, Same as 1, but print an extra copy of the transport          @
c                 summary page at the end of the qikone file.                  @
c% ishot       Doublet III or DIII-D shot number for this analysis run.        @
c                      NOTE: ishot is used as part of the identifier           @
c                            for the eqdsk name for eqdsks created by          @
c                            ONETWO (this has nothing to do with input         @
c                            eqdsks which are required if irguess .ne. 0)      @
c                            Hence you may wish to assign some sort of         @
c                            identifying number to ishot even if this is       @
c                            not an official Yoka-type run.                    @
c                            In tdem mode ishot should be set to the           @
c                            shot number in the netCDF file.                   @
c% itime       Time at which the shot is being analyzed (msec).                @
c  
c  wrt_kinetic_efit   logical variable ,set to true to write file
c              12_Kshot.time which contains two namelists that EFIT
c              uses for kinetic efit analysis. Note that ishot
c              must be given in order to construct the file name!!!!
c              (time is used instead of itime since more tha one file
c               can be created).
c smth_mhd_parms This is a smoothing parameter for the triangularity
c                and elongation put out in the iterdb file. 
c                smth_mhd_parms =0.0 means no smoothing (default) and
c                smth_mhd_parms > 0.0 means some smoothing. A good value
c                is smth_mhd_parms = 0.0001
c ts_smfactor   smoothing for gcnmp code profiles (used only if gcnmp is involed)
c               0.0 no smoothing (default)
c               ts_smfactor =0.001 significant smoothing,etc
c single_density_simulation :
c               used only in gcnmp. if a single ion species is run
c               in simulation mode then setting this value to true will
c               cause all other ion densities to be ratioed according to their
c               initial ratios before any transport.
c ----------------------------------------------------------------------       @
c                                                                              @
      iyoka = 0                                                                %
      ishot = 0                                                                %
      itime = 0                                                                %
      wrt_kinetic_efit = .false.  !dont wrt the files
      smth_mhd_parms   = 0.0
      ts_smfactor      = 0.0
      single_density_simulation = .FALSE.
*                                                                              @
c ----------------------------------------------------------------------       @
c --- SOURCE MULTIPLIERS                                                       @
c ----------------------------------------------------------------------       @
c% wdelt: multiply qdelt by wdelt                                              @
c% wgam : multiply qgam by wgam                                                @
c% wohm : multiply qohm by wohm                                                @
c% angrm2d(i)      i = 1,4 multipliers for 2d angular rotation source terms    @
c                  see documentation for explanation.                          @
c% angrcple         multiplier for energy flux due to angular rotation         @
c                  in ion energy equation.                                     @
c  cximult         multiplier for qcx - charge exchange source term in         @
c                   ion energy equation .                                      @
c ----------------------------------------------------------------------       @
c                                                                              @
      wdelt    = 1.0                                                           %
      wgam       = 0.0  !used to be 1.0,but was not used since qgam was        %
                        !not coded up. Added qgam 11/22/01 so now wgam should
                        !be defaulted to 0.0 to be consisten with old cases.
      wohm       = 1.0                                                         %
      cximult    = 1.0  !multiplier for qcx (charge exchange source term in 
                        !ion energy equation
      angrm2d(1) = 1.0                                                         %
      angrm2d(2) = 1.0                                                         %
      angrm2d(3) = 1.0                                                         %
      angrm2d(4) = 1.0                                                         %
      angrcple   = 1.0                                                         %
*                                                                              @
c ----------------------------------------------------------------------       @
c --- RADIATION PARAMETERS                                                     @
c ----------------------------------------------------------------------       @
c% nqrad       = 0, calculate radiation power density                          @
c              > 0, number of values in an input table of radiation            @
c                     power density (watts/cm**3) as a function of             @
c                     normalized radius                                        @
c% qradr         Normalized r values for table.                                @
c% qradin        Power density values (watts/cm**3).                           @
c                   These may be time-dependent:                               @
c                   qradin(1 - nqrad,m) are for bctime(m).                     @
c                   A linear interpolation in time is used.                    @
c% refrad        Wall reflection coefficient for synchrotron radiation.        @
c                This parameter lies between 0 and 1, with the latter          @
c                value turning off the radiation.                              @
c ----------------------------------------------------------------------       @
c                                                                              @
      nqrad = 0                                                                %
      do   j=1,ksplin                                                          %
        do m=1,kbctim                                                          %
          qradin(j,m) = 0.0                                                    %
        end do                                                                 %
      end do                                                                   %
      refrad = 1.0                                                             %
*                                                                              @
c ----------------------------------------------------------------------       @
c --- NEUTRAL BEAM HEATING PARAMETERS                                          @
c ----------------------------------------------------------------------       @
c     The default is for co-injection. To get ctr injection                    @
c     set angleh to a negative value (usually -19.5 degrees)                   @
c     and switch the sources (sfrac(1) becomes sfrac(2), etc.)                 @
c     Be sure to check the graphic output from program NUBPLT if you           @
c     change the default parameters!                                           @
c                                                                              @
c% iterate_beam  logical, if .true., allow up to imaxbeamconvg iterations      @
c                default: iterate_beam  = .false.                              @
c              Iterate_beam = false is set internally (and user input value is @
c              ignored) if time_dep_beam = 1 is set or if nubeam is selected.  @
c              IF time_dep_beam = 1 and P_Nfreya is selected however then
c              beam iterations are valid.
c% imaxbeamconvg iterations are done to get consistency                        @
c                in the thermal and fast ion densities (ONLY).                 @
c                default: imaxbeamconvg =  5                                   @
c% relaxden      relaxation parameter for beam iteration,must be in (0,1]      @
c                with 1 meaning no relaxation.                                 @
c% relaxden_err  relative error for beam convergence (The relative error is    @
c                measured in the change in electron density if icenez=0        @
c                and as a relative change in the thermal ion species           @
c                corresponding to the beam if icenez=1)                        @
c% timbplt       Times (up to ntimbplt) to produce data for Freya-like plots          @
c                of beam deposition. output is processed by the NUBPLT code.   @
c                Defaulted to OFF.                                             @
c                timbplt(1) .le. time0 .and. beamon .lt. time0 gives o/p       @
c                for initial time. 
c
c                If  time_dep_beam = 0 only index i=1 will be used in the      @
c                following!!
c% beamon(i)     Time (sec) at which beam heating is turned on for all sources
c                of beam line i                                                @
c% btime(i)      Time interval (sec) during which beam heating is on
c                this is the pulse length for each pulse from all sources      @
c                associated with beam line i
c

c  external_beam_cur = 0,1, use this switch to input an external beam driven
c                current as a function of space and time.
c                default =0 , which  calculates the beam driven current internally
c                if external_beam_cur = 1 then the electron drag,the
c                trapped electron correction to the electron drag and the
c                pure ion current  (curbe,curbet and curbi) are not calculated
c                individually and hence will appear as 0.0 in the
c                output.
c                The time and space dependent input quantities for the transp
c                beam drive current are given by
c                         knotscurb_external(1...kbctim)
c                         bparcurb_external(kbctim)
c                         rcurb_external(ksplin,kbctim)
c                         curbeam_external(ksplin,kbctim) 
c                where the first index,ksplin, ranges over knots and
c                second index, kbctime ranges over times in bctime (this
c                structure is identical to the other profiles inputs such as
c                tein, described above)
c                
c  d_beam_spin_pol  The "spin polarization"  of the deuterium beam.
c                 This is used to determine the Dbeam(he3,p)He4
c                 reaction rate. At present this is just a multiplier
c                 on that reaction rate. The default is 1.0
c                 Spin polarized cross sections are not in the code yet. 
c                 NOTE that since this is a simple multiplier a value of zero
c                 will turn off the fusion rate!!! So 1 really means
c                 randomized spin using the spin state averaged cross
c                 section (Bosch and Hale)      HSJ
c  he3_thermal_spin_pol same as d_beam_spin_pol but for thermal he3 density

c   -----            NEW INPUTS USED ONLY IF time_dep_beam =1 --------
c% beamoff(i)     the length of time                                           @
c                all the sources of beam i are off before they turn on again   @
c                Note that all sources of a given beam line are controlled by  @
c                beamon(i),beamoff(i) ,btime(i) ,and beam_end(i)               @
c                Obviously single pulse operation is achieved by setting either@
c                btime or beamoff to a value larger than timmax                @
c                (HSJ fix sub getpow for pulsed beams )                        @###
c  beam_end(i)   the time beyond which beam i is shutdown (e.g. never comes on again)
c  beam_cycles(i)An alternative to inputting beam_end is to give the number of
c                cycles that the beam goes through. This is not an integer,
c                fractionla parts of a cycle are allowed. A cycle is defined as the
c                sum of the on and off time (eg. btime(i)+beamoff(i) )
c                if both beam_end and beam_cycles is set then beam_end is used
c% source2_phase(i)  This is the time lead (if source2_phase(i) is input as    @ 
c                a negative number) or time lag ((if source2_phase(i) is input @
c                as a positive number). Of source 2, beamline i, relative to   @
c                source 1,beamline i in seconds.                               @
c                (if source2_phase(1) = -0.050 for example it means that       @
c                 source 2 of beam line 1 is started 50 msec before source 1   @
c                 of this beam line starts. The start time of source 1 is      @
c                 given by beamon(1)                                           @
c                 if the input nsourc is set to 1 instead of 2 then this has   @
c                 no effect. By making source2_phase(i) a number larger than   @
c                 timmax you can effectively eliminate the second source of    @
c                 any beam line.
c  beam_thermal_cutoff flag for determination of lower limit on slowing down   @
c                 distribution of fast ions, set to -1,0 or 1, 1 is default    @
c                 If   beam_thermal_cutoff = 1 then the thermal speed          @
c                of the ions is used as the lower limit of integration         @
c                (see beam_pulse_control below).                               @
c                if beam_thermal_cutoff =0 then the lower limit of integration @
c                over the fast ion    @
c                distribution function is taken as  therm_frac* vthermal ion   @
c                therm_frac is DEFAULTED TO  1.0                               @
c                Both beam_thermal_cutoff=0 and beam_thermal_cutoff=1 will     @
c                result in a lower limit that depends on Ti and hence rho.     @
c                NOTE that in either case the Fokker-Planck equation used to   @
c                determine the slowing down distribution is not valid    at    @
c                these energies and hence this is an approximation
c                To get a cutoff independent of rho set beam_thermal_cutoff =  @
c                -1 and input the value of the thermal ion speed to be used    @
c                in  beam_thermal_speed ( see below)                           @
c                The other factors that go into the integration limit (ne,te)  @
c                are averaged over the grid, using a straight average at       @
c                present. see tau0_avg subroutine                              @
c  beam_mode     character variable set to "test" or "run"                     @
c                ( "run" is default) in test mode Freya is not called. this    @
c                 mode is intended primarily to build up waveforms of beam     @
c                input before doing a transport run for developers only.       @
c beam_restart_file     A file, if present in the working directory, 
c                       will be used to set beam initial conditions, 
c                       provided time0 matches the  input values in inone 
c                       and beam_init_restart_file = 0.
c                       ( beam_restart_file can hold up to 256 characters 
c                         - no check or error exit at this time )
c                       Note that this is NOT the restart file for Nubeam !!!!!
c beam_init_restart_file   
c                Set to 1 to create file beam_restart_file in this run.
c                (default 0) if file exists it will be replaced by a new one.
c
c
c ngauss         the number of gaussian quadrature points to be used in        @
c                carrying out the integrations. This is a accuracy versus      @
c                execution speed issue. This value probably does not have to   @
c                be changed by the causual user. The maximum value is          @
c                ngauss_max (see  gauss_info.i)                                @
c
c beam_pulse_control  This integer valued switch determines how the upper and  @
c                lower limits of integration on the fast ion distribution      @
c                function (and all required moments thereof) are handled in    @
c                determination of the fundamental time step dt.
c                if beam_thermal_cutoff = -1 then beam_pulse_control is forced @
c                to be 0 ( no matter what the user set it to)

c                Besides the    @
c                obvious fact that all pulse on and off times must be          @
c                reflected in the time                                         @
c                steps taken by the code we also have to consider those times  @
c                at which the upper and lower limits of integration reach      @
c                the saturated value of the thermal cutoff.
c                But this thermal cutoff normally depends on the local         @
c                temperature and hence is different at every grid point.       @
c                Thus we see that dt will have to be modified many times to    @
c                account for each beam,injection energy, source,grid point and @
c                pulse number. This means that many time steps will be         @
c                generated internally. This is what will happen if             @
c                beam_pulse_control is left at its default value of 0          @
c                and beam_thermal_cutoff = 0 or 1.

c                Set  beam_pulse_control = -1 if you      want to skip the     @
c                modification of the time step due to this effect altogether   @
c                In this case there will be some error in the calculated       @
c                beam integrals (because we wont hit the thermalization time   @
c                exactly)The amount of this error is then controlled           @
c                by the maximum allowed time step, dtmax. If you select        @
c                beam_pulse_control = -1 then you should cut dtmax  down from  @
c                its default value appropriately.                              @

c                beam_pulse_control = n, where n >0 and n < nj, ( nj           @
c                is the number of rho grid points)  means that the plasma is   @
c                broken up into n regions, with the beam thermalization        @
c                time calulated at the                                         @WARNING not yet done
c                mid point of each region. This means that instead of nj mods  @
c                of the time step for each beam,energy,source,pulse, there will@
c                only be n such interruptions. What we have in mind here is that@
c                n=1,2,3,4,5, etc but any n up to nj-1 can be used.            @
c                For n=1 for example the slowing down parmeters will be        @
c                evaluated at nj/2 and held constant at that value at all      @
c                grid points assuming beam_thermal_cutoff =0 or 1              @ 
c
c                Note that  beam_thermal_cutoff = -1 is the same as setting    @
c                beam_pulse-control = 1 but  using beam_thermal_speed for the  @
c                single value to be used as the lower limit of integration     @
c                at all the rho grid points.
c
c
c    
c beam_thermal_speed          in cm/sec, used when beam_pulse_control =0       @
c                the default value corresponds to an ion temperature of 50ev   @
c                

c           
c -----            END of       NEW INPUTS USED ONLY IF time_dep_beam =1 ----      
c           
c                
c% ibcur         Flag for neutral beam driven current                          @
c                  1, include beam driven current                              @
c                  0, neglect beam driven current                              @
c%  use_Callen    integer,has meaning only if mcgo run is done.                @
c                 Intended for checking out some model details.                @
c                 Users should leave it set at the default value of 0.         @
c% ibcx          Flag for secondary charge exchange between fast ions          @
c                  and thermal neutrals                                        @
c                  1, include secondary charge exchange (default)              @
c                  0, neglect secondary charge exchange                        @
c                  NOTE: thermal ion particle source due to secondary charge   @
c                        exchange is not included in the ion density equation, @
c                        independent of the setting of ibcx.                   @
c                                                                              @
c% nbeamtcx    switch used to determine the form of the beam torque.           @
c              if nbeamtcx = 0 (default) then the loss of beam torque due      @
c              to secondary charge exchange of fast ions with thermal          @
c              neutrals is neglected. if nbeamtcx = 1 the transfer of angular  @
c              momentum from fast ions to thermal ions and electrons is        @
c              modified to account for the cx losses (done with array ssprcxl) @
c% ibslow      Flag for neutral beam slowing down                              @
c                1, include neutral beam slowing down                          @
c                0, neglect neutral beam slowing down                          @
c                If the new multiple beam pulse model is selected              @
c                (time_dep_beam =1) then ibslow =1 is enforced by the code     @
c% fionx       Allows testing for non-classical slowing down                   @
c                  (see subroutine SLOW1).                                     @
c                                                                              @
c% fast_ion_target  integer,if set to 1 will let the incoming neutral beam     @
c                   see the stored fast ion density                            @
c                   no attempt to modify the stopping rates due to the         @
c                   fact that the fast ion distribution is not Maxwellian      @
c                   has been made!!!!!!!!!!!!!!!!!!!!!!!                       @
c                   Note that if the stored fast ion density is large enough   @
c                   that it is of concern then it is also large enough so that @
c                   a simple one pass linearization doesnt make sense. So you  @
c                   should use this option in conjunction with the iterate     @
c                   beam option.                                               @
c                   For Mcgo coupled runs the birth points of the fast ions    @
c                   are determined by the Freya in Onetwo. Hence this switch   @
c                   will affect the Mcgo results as well. Note that Mcgo       @
c                   itself does not let the Monte Carlo particles see the      @
c                   stored fast ion density. As an approximate remedy for      @
c                   this situation you can set fast_ion_target = -1 when       @
c                   running a Mcgo coupled case. This will add the stored fast @
c                   ion density determined in Freya to the thermal ion density @
c                   used in Mcgo, so that the Mcgo particles see a correct     @
c                   total density.                                             @
c                   (If Mcgo is run with fast_ion_target = 1 then              @
c                   the birth points feed to Mcgo will account for the stored  @
c                   fast ions but subsequent slowing down of the Monte Carlo   @
c                   particles  will see only the thermal ions)                 @
c% rtstcx       A factor between zero and 1 used to fudge secondary charge     @
c               of fast ions. ie the original fast ion charge exchanges        @
c               with a thermal neutral. The resulting fast neutral             @
c               will however most likely be reionized. To model this           @
c               reionization in Onetwo we don't do the secondary charge        @
c               exchange calculation, which requires a Monte Carlo approach.   @
c               Instead we assume that the probability against charge exchange @
c               of a fast ion, normally taken as                               @
c               {(vbeam**3+vcrit**3)/(v**3+vcrit**3)}**A                       @
c               where   A=(-taus/(3.*taucx)) is the ratio of slowing down      @
c               to charge exchange lifetime is given instead by                @
c               {(vbeam**3+vcrit**3)/(v**3+vcrit**3)}**B                       @
c               where B=rtstcx*A. To decrease charge exchange losses           @
c               set rtstcx less than 1. Note that this can be interpreted      @
c               as assuming that the neutral density, enn, is replaced by      @
c               rtstcx*enn for the purposes of charge exchange probability     @
c               calculations ONLY. The actual neutral density is not changed!  @
c                                                                              @
c               Set rtstcx to a negative number to use sigma*v, rather         @
c               than <sigma*v> (i.e., Maxwellian average), to determine        @
c               the mean lifetime of a fast ion against charge exchange.       @
c               The energy and velocity used will be that of the beam          @
c               corrected for bulk rotation of ions (note the assumption       @
c               that neutrals rotate as ions do).                              @
c               If rtstcx is set to a number less than -5. then                @
c               the average fast ion energy (again corrected for rotation)     @
c               for each beam slowing down distribution is used.               @
c                                                                              @
c% nameb        Name of neutral species in beam                                @
c              'h', 'd', 't', 'dt'                                             @
c               MUST BE PRIMARY ION SPECIES  #1 or #2.                         @
c               if nameb = 'h' then fdbeam defaults to natural isotopic        @
c                  content of d in 'h'.                                        @
c               if nameb =  't' fdbeam is explicitly set to 0                  @
c               if nameb = 'dt' fdbeam is set to fd, (the fraction             @
c               of d in thermal dt mixture)  USER HAS NO CHOICE                @
c                                            UNLESS IFUS = -1 is selected      @
c               NOTE that nameb=dt selects a SINGLE fast ion fluid             @
c               with a fictitious mass !!!!                                    @
c% fdbeam      Fraction of deuterium atoms in neutral beam.                    @
c% relnub      Minimum relative change in ion density or electron              @
c              temperature between neutral beam calculations                   @
c              (0.10 is suggested)                                             @
c% tfusbb      'thermal fusion balance for beams', fraction by which           @
c              the net energy gain from thermal fusion must exceed             @
c              the net energy loss for automatic beam turnoff.                 @
c              If tfusbb = 0, automatic beam turnoff is not done.              @
c% iddcal     Flag controlling treatment of beam effects on fuscal,            @
c             the calculated fusion neutron rate:                              @
c             0 = do not include knock-on or beam-d neutrons in fuscal         @
c             1 = include only knock-on neutrons in fuscal                     @
c             2 = include only beam-d neutrons in fuscal                       @
c             3 = include both knock-on and beam-d neutrons in fuscal          @
c% ranseed    Starting seed for random12 number generator used in the            @
c               FREYA determination of the beam deposition                     @
c% npart      Number of particles followed into plasma (suggest 10000)         @
c% npskip     Ratio of number of particles followed into plasma                @
c               to number of source particles (suggest 1)                      @
c             npskip =1 is set in Nfreya if this is a parallel run with more   @
c             than one processor                                               @
c% iborb      Flag for modeling orbit effects on beam-generated fast ions      @
c             3 = use the Monte Carlo Code MCGO to model the fast              @
c                 ion slowing down. (The ionization of the injected neutrals   @
c                 is still done by FREYA as usuall) This is a very expensive   @
c                 option and shouldn't be used indiscriminately. It will take  @
c                 about 30 mins of CPU time (ON WHAT MACHINE?) to follow 2500  @
c                 ions and generate fast ion density and deposition profiles.  @
c                 If iborb=3 is used the following input is required.          @
c                      npart_mcgo     the number of fast ions to be followed   @
c                                     during the slowing down process. MCGO    @
c                                     will follow the guiding center of each   @
c                                     ion until it is thermalized or lost      @
c                                     (by bad orbit or charge exchange)        @
c                                      The maximum you can set here is 10000   @
c                                      Typically 2000-5000 ions are sufficient @
c                                      to yield a reasonable sampling.         @
c                                      DO NOT CONFUSE NPART_MCGO with the      @
c                                      FREYA NPART SETTING. Npart and          @
c                                      npart_mcgo do not have to be the same   @
c                                      but npart MUST be >= npart_mcgo.        @
c                                      MCGO will take a random sample of       @
c                                      npart_mcgo initial fast ion launch      @
c                                      parameters from the list of length      @
c                                      npart generated by FREYA. Proper        @
c                                      consideration is given to the full      @
c                                      half and third energy components of the @
c                                      beam. (Npart_mcgo is the total for all  @
c                                      three components)                       @
c                                                                              @
c             -3      iborb= -3 means that all of the mcgo preparatory         @
c                     calcs will be done and the input files required to run   @
c                     mcgo will be generated but mcgo itself will not be run   @
c                     Using this option you can run mcgo manually,             @
c                     (at a later time) rather than                            @
c                     having onetwo control the execution. Note that iborb =-3 @
c                     disables most of the freya calculations in onetwo so it  @
c                     should not be used for normal transport runs.            @
c                                                                              @
c                     In order to use this option the files                    @
c                        mcgo_output_file(1,2)                                 @
c                     mentioned below must not be specified in inone           @
c                     (or specify them as zero-length files:                   @
c                        mcgo_output_file(1) = '' )                            @
c                                                                              @
c                     Another way to run onetwo is to use pre-existing mcgo    @
c                     results rather than have onetwo spawn mcgo during a      @
c                     transport run. To use this option set iborb =-3 and      @
c                     specify the name of the file(s) that onetwo is to read   @
c                     instead of running mcgo. One such file will be read for  @
c                     each beam so the number of files you specify has to      @
c                     equal the value of nbeams in inone. The files are        @
c                     specified in namelist 2 according to                     @
c                          mcgo_output_file(1) = 'XXXXX'                       @
c                          mcgo_output_file(2) = 'XXXXX'                       @
c                     (substitute actual file names for the XXXXX )            @
c                     A time-dependent version of this has NOT been programed! @
c                     So the same files will be read over and over again in a  @
c                     time-dependent ONETWO run. Obviously the propensity to   @
c                     get out of sync with the inone file is great here so use @
c                     this option cautiously !!!!!!                            @
c             mcgo_12_output       MCGO will create an output file that        @
c                                      can be read by Onetwo.                  @
c                                      You can use                             @
c                                      this file on subsequent Onetwo runs     @
c                                      as well, without running MCGO,          @
c                                      by specifying the file name of the      @
c                                      MCGO data in the variable               @
c                                       mcgo_results_12                        @
c                                      if mcgo_results_12.nc is set to         @
c                                      a valid file                            @
c                                      then this file will be used even if     @
c                                      iborb=3 (and MCGO will not be executed) @
c                                      If mcgo_results_12    is not a valid    @
c                                      file name                               @
c                                      then MCGO will be executed only if      @
c                                      iborb=3. In this case a file with       @
c                                      the name mcgo_onetwo.dat will be        @
c                                      created.                                @
c                       We would like to do this using PVM (parallel           @
c                       virtual machine) instead of a file                     @
c                       interface eventually so the above may change.          @
c                       Other information that MCGO needs is part of the       @
c                       standard Onetwo input so no further user intervention  @
c                       is required.                                           @
c                                                                              @
c --- MCGO FILE USAGE SUMMARY --- begin -------------------------------------- @
c                                                                              @
c FORTRAN VARIABLE      Function and Value of Variable                         @
c ----------------      ------------------------------                         @
c                                                                              @
c% mcgo_input_file(ib)      is set equal to the name 'mcgo_input_12_beam1' 
c                       and 'mcgo_input_12_beam2' in cray001.f  @
c                       It is an input file for MCGO created by ONETWO         @
c                       (the namelist input)                                   @
c                                                                              @
c% mcgo_input_file2(ib) these files are set equal to 'mcgo_input_freya_beam1'  @
c                       and (if there are two beams) 'mcgo_input_freya_beam2'  @
c                       in cray001.f                                           @
c                       These are additional input files for MCGO created by   @
c                       ONETWO. (They contain the fast ion deposition          @
c                       determined by FREYA. These files are unformatted for   @
c                       speed. There is one such file for each beam.)          @
c                                                                              @
c% mcgo_output_file(ib) The name used in inone to specify files created        @
c                       by MCGO. The data in these files is to be read         @
c                       into ONETWO (instead of spawning MCGO from ONETWO).    @
c                       Unless the user changes these names on disk, the       @
c                       values will be 'mcgo_12_output_beam1' and,             @
c                       for two beams, 'mcgo_12_output_beam2'                  @
c                                                                              @
c --- MCGO FILE USAGE SUMMARY --- end ---------------------------------------- @
c                                                                              @
c% iborb      2 = new prompt loss model of Stambaugh (implemented only for     @
c                 codeid = dee i.e., two d case)                               @
c              (REF: "Calculating Orbit Loss in a Separatrix Bounded Tokamak", @
c                     DIII-D Physics Memo No. 9303, 16 March 1993)             @
c             1 => do     model orbit effects (this is the original method)    @
c             0 => do not model orbit effects                                  @
c% itrapfi    Flag for trapped fast ions.  If iborb = 1 then                   @
c             setting itrapfi = 1 will modify the beam driven                  @
c             current by the initial trapped ion fraction.                     @
c             (Subsequent pitch angle diffusion is NOT taken into account).    @
c             itrapfi = 0 neglects this effect.  If iborb = 0,                 @
c             itrapfi has no effect.  itrapfi = 0 is default.                  @
c% ds_tk     maximum trajectory increment (cm) used in subroutine INJECT to    @
c            calculate psi(s), where s is the neutral trajectory               @
c            pathlength from the first closed flux surface it                  @
c            encounters.  used for non-zero toroidal rotation, where           @
c            mean free path as a function of path length is required.          @
c% fe_tk     factor (>1) to set upper limit of energy in n*sigma array.        @
c            max(ebins) = max(ebkev)*fe_tk.  required for nonzero              @
c            rotation cases. fe_tk=1.1 is default                              @
c% ne_tk     number of equi-width energy bins used in forming n*sigma          @
c            array.  required for nonzero rotation cases.  is internally       @
c            reset to zero (used as flag to turn off rotational effects        @
c            on neutral stopping) if angular rotation is not present           @
c            (iangrot = 0).                                                    @
c% iexcit    Selection switch for atomic cross section data and model:         @
c            0,  Use the fundamental atomic data of Freeman & Jones (1972).    @
c                  This option is not advised since this atomic data is        @
c                  considered outdated and excited state effects are           @
c                  not considered.                                             @
c            1,  Use hexnb package but do not include excitations in           @
c                  its calculation of cross sections.                          @
c            2,  Use hexnb package, include excitations in calculation         @
c                  of cross sections.                                          @
c                  NOTE: Options 1 and 2 are not advised because the hexnb     @
c                  model has been found to be flawed and is based on           @
c                  atomic data that is considered outdated.                    @
c                  (Boley et al, Nuclear Fusion 29, 1984)                      @
c            5,  Use the JET-ADAS effective stopping cross sections.           @
c                  This model is preferred since it is based on the most       @
c                  recent atomic data available and includes multi-step        @
c                  ionization processes due to excited states. This is an      @
c                  important effect for considering neutral beam penetration   @
c                  into fusion-grade plasmas. For a detailed description       @
c                  of the ADAS Atomic Data and Analysis Structure (ADAS)       @
c                  developed by JET, see Finkenthal, 1994.                     @
c                  Note: All ions are considered fully stripped.               @
c                  (Daniel Finkenthal Ph.D. Thesis, UC Berkeley, 1994)         @

c            6     Boley's parameterization including MSI effects:
c                                                                              @
c% neg_ion_source(i)  i=1,...,nbeams                                           @
c                     an integer array, set to 1 to indicate                   @
c                     that a NEGATIVE ion source is used for                   @
c                     neutral beam line i. At this time, the only              @
c                     effect of this switch is to set the neutralization       @
c                     efficiency, arbitrarily, to 98%, independent of the      @
c                     negative ions energy and to eliminate the second and     @
c                     third energy components from the beam.                   @
c% time_dep_beam      set to 1 to indicate that new time dependent beam input  @
c                     will be used. (default = 0).(NOTE: if this is a snapshot @
c                     run then the code will set time_dep_beam =0 even if user @
c                     selected time_dep_beam =1)
c                                                                              @
c% izstrp(i) i = 1,2..nimp set to 0 for coronal equilibrium                    @
c            values of <z> and <zsq>. set to 1 for fully stripped impurities   @
c% mstate    principal quantum number n above which excitations                @
c            count as ionizations.                                             @
c% inubpat   two-dimensional beam deposition option                            @
c            0,  do not calculate beam deposition in (r,z) coordinates         @
c            1,  calculate beam deposition on (r,z) grid (default = eqdsk      @
c                grid).  write deposition array and n = 3 excited state        @
c                fraction to file 'beamdep' for standalone analysis.           @
c% npat      neutral beam deposition (r,z) grid dimensions, used if            @
c            inubpat = 1.0                                                     @
c            npat(1)  =  number of elements in 'r' (<=2*nw)                    @
c            npat(2)  =  number of elements in 'z' (<=2*nh)                    @
c            defaults, npat(1) = nw, npat(2) = nh                              @
c% mf        Number of flux zones plus 1 (max = 81)                            @
c                                                                              @
c  In the following list the index ib designates the beam injector,            @
c  while ie refers to one of the three energy components.                      @
c  iap refers to one of the apertures.                                         @
c                                                                              @
c% nbeams          Number of neutral beam injectors (1 to kb)
c% nsourc          Number of sources per beamline.                             @
c                    If 1, source is centered on beamline axis.                @
c                    If nsourc = 2, distinguish between the beamline           @
c                    axis and the source centerline (optical axis).            @
c                    The two sources are assumed to be mirror images           @
c                    through the beamline axis.                                @
c                    In either case, the exit grid plane is perpendicular      @
c                    to the beamline axis, and contains the source             @
c                    exit grid center(s).                                      @
c                    If nsourc = 2, the alignment of the sources w.r.t.        @
c                    the beamline axis is specified through bhofset,           @
c                    bvofset, and bleni (described further below).             @
c% naptr           Total number of apertures encountered by a particle         @
c                    as is moves from the source into the plasma chamber.      @
c                    Maximum is specified by parameter nap ( = 4) .            @
c                    First set of apertures encountered by the particles       @
c                    are assumed centered on the source axis, and subsequent   @
c                    apertures are centered on the beamline axis;              @
c                    the distinction is made through ashape.                   @
c% anglev(ib)      Vertical angle (degrees) between optical axis               @
c                    and horizontal plane; a positive value indicates          @
c                    particles move upward                                     @
c% angleh(ib)      Horizontal angle (degrees) between optical axis and         @
c                    vertical plane passing through pivot point and            @
c                    toroidal axis; a zero value denotes perpendicular         @
c                    injection, while a positive value indicates par-          @
c                    ticles move in the co-current direction                   @
c% bvofset(ib)     Vertical offset from beamline axis to center                @
c                    of each source (cm; used only for nsourc = 2)             @
c% bhofset(ib)     Horizontal offset from beamline axis to center              @
c                    of each source (cm; used only for nsourc = 2)             @
c% bleni(ib)       Length along source centerline (source optical axis) from   @
c                    source to intersection point with the beamline axis.      @
c% sfrac1(ib)      Fraction of source current per beamline coming              @
c                    from upper source (used only for nsourc = 2)              @
c                    (the upper source is the more perpendicular,              @
c                      or right source normally)                               @
c% bcur(ib)        Total current (a) in ion beam (used only if bptor           @
c                    is zero)                                                  @
c% bptor(ib)       Total power (w) through aperture into torus; when           @
c                    nonzero, bptor takes precedence over bcur                 @
c% nbshape(ib)      Beam shape                                                 @
c                    'circ':  circular                                         @
c                    'rect':  rectangular                                      @
c                'rect-lps':  rect. long pulse source (DIII-D only)            @
c                             a choice of short or long pulse sources is       @
c                             available by injector (not by source).  one      @
c                             or both injectors may be long pulse by           @
c                             setting one or both to 'rect-lps'                @
c                  CAUTION:  DIII-D sources are defaulted to lps specs.        @
c                  It is the user's responsibility to overide these for        @
c                  sps configuration(s).    
c                  Default is nbshape(1:nbeams) = 'rect-lps'              
c% bheigh(ib)      Height of source (cm)                                       @
c% bwidth(ib)      Width of source (cm); diameter for circular source.         @
c% bhfoc (ib)      Horizontal focal length of source (cm)                      @
c% bvfoc (ib)      Vertical focal length of source (cm)                        @
c% bhdiv (ib)      Horizontal divergence of source (degrees)                   @
c% bvdiv (ib)      Vertical divergence of source (degrees)                     @
c% ebkev (ib)      Maximum particle energy in source (keV)                     @
c% fbcur (ie,ib)   Fraction of current at energy ebkev/ie                      @
c                  Note that this is the current fraction at the source,       @
c                  before it enters the neutralizer.                           @
c% ashape(iap,ib)  Aperture shape.                                             @
c                   Prefix 's-' indicates source axis centered.                @
c                   Prefix 'b-' indicates beamline axis centered.              @
c                     's-circ'          'b-circ'                               @
c                     's-rect'          'b-rect'                               @
c                     's-vert'          'b-vert'                               @
c                     's-horiz'         'b-horiz'                              @
c                                       'b-d3d'                                @
c                    circ = circular aperture,                                 @
c                    rect = rectangular,                                       @
c                    vert = limits vertical height of source particles,        @
c                   horiz = limits horizontal height of source particles,      @
c                     d3d = special DIII-D polygonal aperture                  @
c% aheigh(iap,ib)  Height of aperture (cm)                                     @
c% awidth(iap,ib)  Width  of aperture (cm); diameter for circular              @
c                    aperture                                                  @
c% alen(iap,ib)    Length from source to aperture for 's-type' apertures,      @
c                    and from exit grid plane along beamline axis for          @
c                    'b-type' apertures.                                       @
c% blenp(ib)       Distance along beamline axis from source exit               @
c                    plane to the fiducial "pivot" point.                      @
c% rpivot(ib)      Radial position of pivot (cm)                               @
c% zpivot(ib)      Axial position of pivot (cm)                                @
c                                                                              @
c% hdepsmth    set this param to a positive value (gt.0.0 and .le.  10) to turn off         @
c              the smoothing of hibrz and hdepz in subroutine POSTNUB.         @
c              if this option is used then enough zones must be specified      @
c              for adequate resolution (zones = number of radial grid points)  @
c              and enough injected neutrals must be followed to minimize       @
c              the statistical noise enough so that no greatly uneven          @
c              profiles result! this option was added because the smoothing    @
c              of the profiles by subroutine SMOOTH can lead to unphysical     @
c              peaking of the birth and deposition profiles.                   @
c              Matching of Monte Carlo results for fast ion distributions,     @
c              especially in the presence of mhd activity, indicates that      @
c              smoothing with hdepsmth may be required. Hence                  @
c              if hdepsmth .gt. 10 then ipass = hdepsmth -10 passes are made   @
c              over the data to get a smooth profile. Each pass averages       @
c              nhdep = 6 (not adjustable) grid points together. Note in        @
c              particular that if the grid is coarse then 6 points cover       @
c              a wider range in rho than if the grid were finer. hence varying @
c              the grid size  probably also means adjusting hdepsmth to        @
c              maintain a more or less                                         @
c              constant deposition profile. Note that the smoothed profiles    @
c              are renormalized to the plasma volume.                          @
c              It is possible to account for the fixed nhdep by adjusting      @
c              hdepsmth. As an example it was observed that with 51 grid       @
c              points and hdepsmth =45 a fast ion deposition profile resulted  @
c              which  was subsequently matched using 201 grid points by        @
c              changing hdepsmth from 45 to 600. (mf = 12 and typical DIII-D   @
c              densities for h mode shots)                                     @
c  hdepsmth in P_Nfreya.  P_Nfreya has an option to spread (smooth) hibrz 
c  fidiff_on               and hdepz  using hdepsmth and/or using fast ion diffusion
c              THESE RULES APPLY TO P_NFREYA not to Nfreya or nubeam !!!
c              The rules are as follows and apply only to P_Nfreya:
c                     1) IF hdepsmth is set to -1 in inone no smoothing        
c                        of hibrz and hdepz will be done at all, 
c                        regradles of setting of fidiff_on 
c                     2) if fidiff_on = FALSE then smoothing is through
c                        hdepsmth. In this case valid values of hdepsmth
c                        range from 11. to 999. ( A good choice is 12.)
c                        The value of hdepsmth determines the number of
c                        smoothing passes that are made over the hibrz
c                        and hdepz profiles: no passes = hdepsmth -10.
c                     3) if fidiff_on = TRUE and hdepsmth .GT. 10.0 
c                        then instead of making smoothing
c                        passes over the data as above, fast ion diffusion
c                        is used to smooth the profiles. The number of diffusive
c                        steps is as in 2 above (eg hdepsmth -10)
c                     4) if fidiff_on = TRUE and  hdepsmth = 0.0  then
c                        a full steady state fast ion diffusion model is used.
c                        This option is  under development so dont use
c                        it until further notice. 
c                     5) fidiff_on = TRUE  requires that a space,energy and
c                        beamline dependent diffusion coefficient is given.
c                        For compatibility this is done through the nubeam scheme
c                        whereby the diffusion coefficent is given by
c                        D   = ADIFF_A + (ADIFF_0 - ADIFF_A) * 
c                                  (1- rho ^ ADIFF_XPIN) ^ ADIFF_XPOUT (in cm^2/s)
c                     NOTE  that the parameters adiff_a,adiff_0,adiff_xpin,
c                         adiff_xpout are, in P_Nfreya only, assumed to be 
c                         function of beam energy component and beamlet number
c                         Hence the input could be, for example, adiff_a(1:3,1:14)
c                         for the three  beam energy componenets and the (currently 14)
c                         defined beamlets. 
c                         adiff_a,adiff_0,adiff_xpin,and adiff_xpout are input
c                         in the nubeam namelist !!! NOT in INONE
c                         The nubeam namelist is in a file identified in inone
c                         by the variable BEAM_DATA_NAMELIST.  The input for  
c                         addif_a and adiff_0 is in cm^2/sec.
c                    6) fidiff_on = .FALSE. means that even if addif_a and
c                        adiff_0 are set (in BEAM_DATA_NAMELIST file)
c                        to give non zero d, no diffusive smoothing will be done.
c                        IN other words the hierarchy is
c                            a) hdepsmth = -1 NO smoothing of any kind
c                            b) fidiff_on = FALSE no smoohting with d
c                     NOTE: The default (as set here) is to TURN OFF
c                        anomalous beam ion diffusion (eg fidiff_on = .FALSE.)   
c% freyavb     integer switch used for debug purposes                          @
c              if  freyavb .ne. 0  then some FREYA-related diagnostic          @
c              output will be written to the screen                            @
c  bfr_neutrlz   USED only in P_Nfreya, bfr_neutrlz = TRUE  means that the current
c               fractions,fbcur, to be used in P_Nfreya are given before the
c               neutralizer. IF bfr_neutrlz = FALSE then it is assumed that
c               the input fbcur is given after the neutralizer. Note that when
c               running P_Nfreya it is assumed that fbcur is specified in the
c               nubeam namelist input files where fbcur is identified with 
c               ffulla,fhalfa
c  ----------------------------------------------------------------------------
c  P_Nfreya neutral beam interface:
c  ----------------------------------------------------------------------------

       use_P_Nfreya  = .FALSE.
       P_Nfreya_read = .FALSE.   ! not an input, toggled in appropriate palces
       P_Nfreya_dt   = 0.015     ! time interval between P_Nfreya calls
                                 ! also conrolled by relnub
       sent_ufile    = .FALSE.   ! send ufile for P_Nfreya only once

c -------- NUBEAM related inputs for beams and related -------------------------
c ---------active only if use_nubeam = true -----------------------------
c ----- defaults are set here ( may be reset  by  namelist reads)
    
        use_nubeam = .false.
                    !set to true to use NTCC  
                    !nubeam package (default .FALSE>)


        nubeam_path =''
c                     specifies the fully qualified executeable name of the nubeam
c                   driver that Onetwo will call. Only users that are capable of
c                   switching out nubeam for a different version will be
c                   able to use this option. For all other users the default
c                   should not be changed. (Needs to be initialized because Len_trim
c                   function wont return correct result otherwise!)

       nubeam_setup_ext ='' 
c                   the name of an executable file that will be sourced before
c                   nubeam is executed for the first time. It is intended that
c                   special environmental variables required by nubeam 
c                   are set in this file.
c                   You need to set this value only if you want to bypass the
c                   standard setup.

       beam_data_namelist = 'set_this_in_inone'
                     !this is the file  written by nblist.for (a transp routine)
                     !it has in it the beam geometry and the order of the
                     !beam data in the ufile, the beam on and off times, etc.
                     !This file MUST  be  available if use_nubeam = .true.
                     !if use_nubeam = .false. and this file is present then
                     !it will be used to generate beam input for the
                     !standard Onetwo nfreya package! <== NOT implemented
                     !The parameters nbeam, tbona(1:nbeam),tboffa(1:nbeam),
                     !pinja(1:nbeam),einja(1:nbeam), ffulla(1:nbeam) and
                     !fhalfa(1:nbeam) are  in this file and
                     !the values are picked up as a function of time from file
                     !beam_data_ufile instead. However if these parameters
                     !have acceptable values for all nbeam then it is assumed
                     !that the user  wants to use these (constant) values
                     !instead of the time dependent ones given in 
                     !beam_data_ufile (see below) . In this case 
                     !eam_data_ufile will not
                     !be used. and need not be present.
                     !IF use_nubeam = F but file beam_data_namelist is 
                     !present and has valid input then the beam input in file
                     !inone will be overwritten with the values in 
                     !beam_data_namelist . If the input in beam_data_namelist
                     !is not valid then an error exit is taken.
 


       beam_data_ufile = beam_data_namelist
                     !beam_data_ufile supplies the beam data in 2d ufile
                     !format. This file should normally  be present 
                     ! if use_nubeam = T. However see the description  
                     ! above for file beam_data_namelist for the case 
                     ! where beam_data_ufile will not be required.
                    
                     !IF use_nubeam = F then this file is optional.
                     !if it is present then it will be used to construct
                     ! the necessary beam data (and any specification of
                     ! beam data given in inone will be ingnored!!!!!)
                     ! For each beam,  beam_data_ufile contains the four items
                     ! pinja,einja,ffulla,fhalfa as a function of time.
                     ! This information is linearly interpolated in time to get
                     ! the values required by Onetwo at any specific time.
                     ! (the file name  given by beam_data_ufile and 
                     ! beam_data_namelist may  include a path specfication)

        
       nubeam_restart  =0 ! normal onetwo/nubeam run. The restart files 
                          ! (for nubeam), if they exist, are not used. Instead
                          ! the restart files are  created as part of this run.
                          ! Any existing files with the names 
                          !     nubeam_12_nubeam_state.cdf
                          !     nubeam_12_xplasma_state.cdf
                          ! are OVERWRITTEN  by the new ones.

                          ! If nubeam_restart=1 
                          !  Use the previously created restart files
                          ! as startup for nubeam (note that Onetwo
                          ! will start with the inone specified profiles.
                          ! The restart only applies to nubeam)
                          ! If the times in inone do not match the time
                          ! in the restart files a warning message will
                          ! be printed by nubeam. However it is save to
                          ! ignore that message.

                          ! if nubeam_restart =1 is selected then a fully
                          ! qualified path name can be given that specifies
                          ! the location of these files. 
                          ! For example:

                          !  nubeam_state_path  = '/u/username/mynubeam_12_nubeam_state.cdf'
                          !  nubeam_xplasma_path ='/u/username/mynubeam_12/mynubeam_12_xplasma_state.cd
                          !   nubeam_profile_path = '/u/username/mynubeam_12/mynubeam_12_profiles_path
                          ! The default is the current working directory.

                          ! If  nubeam_state_path and or nubeam_xplasma_path
                          ! and/or nubeam_profile_path 
                          ! is not specified then files with the names 
                          !    nubeam_12_nubeam_state.cdf
                          !    nubeam_12_xplasma_state.cdf
                          !    nubeam_12_restart_profs.txt
                          ! are looked for in the current working directory.
                          ! If these files are found they are copied to
                          ! nubeam_12_nubeam_state.cdf_orig 
                          !               and 
                          ! nubeam_xplasma_path.cdf_orig
                          !               and
                          ! nubeam_12_restart_profs.txt_orig
                          ! so  that Onetwo wont overwrite them.

                          ! If nubeam_state_path , nubeam_xplasma_path ,and nubeam_profile_path
                          ! point to valid files they will be copied to the
                          ! current working directory and renamed to
                          !    nubeam_12_nubeam_state.cdf
                          !    nubeam_12_xplasma_state.cdf
                          !    nubeam_12_restart_profs.txt


                          !  Note that the final vesion of these
                          ! files that reside in the curent working directory
                          ! are at the final  time that the nubeam 
                          ! package is called. Additionally if the variable
                          ! wrt_restart_file_time is set then the restart files 
                          ! at that time will be saved.

                          !  Use the previously created restart files
                          ! as startup for nubeam (note that Onetwo
                          ! will start with the inone specified profiles.
                          ! The restart only applies to nubeam)
                          ! If the times in inone do not match the time
                          ! in the restart files a warning message will
                          ! be printed by nubeam. However it is save to
                          ! ignore that message.

                          ! JMP START
                          ! If nubeam_restart=-1
                          ! A stand-alone Nubeam run 
                          ! form time0 - nubeam_back_delt to time0 using time0 
                          ! information to provide an initial guess on fast-ion distribution 
                          ! at time0.

                          ! All plasma profiles and mhd equilibrium are fixed with time0 values 
                          ! during this stand-alone Nubeam run.If use_ufile = .true., 
                          ! however, the beam power is time-dependent 
                          ! according to beam_data_ufile. 
                          ! If use_ufile = .false., tbona should be time0.  
 
                          ! nubeam_restart = -1 can be used in TDEM mode 
                          ! with use_ufile = .true. or .false.
                          ! nubeam_restart = -1 is a default setting for kinetic EFIT.

                          ! nubeam_restart = -2
                          ! same as nubeam_restart = 0, but bypass a test on tbona > time0  
                          ! JMP END

      nubeam_back_delt = 0.0 
                     ! JMP
                     ! used when nubeam_restart = -1

      use_ufile = .false. 
                     ! JMP
                     ! If use_ufile = T, the beam power (pinja) is 
                     ! read from the UFILE as a function of time along 
                     ! with the beam energy (einja) and energy distribution (ffulla and fhalfa). 
                     ! The UFILE name is specified by the switch beam_data_ufile. 
                     ! In this case, the relevant variables in the beam_data_namelist file 
                     ! (tbona, tboffa, pinja, einja, ffulla and fhalfa) are neglected. 

      ifix_nubeam_dt = 0 
                     ! JMP
                     ! If ifix_nubeam_dt = 1, 
                     ! the nubeam time step is set by nubeam0_dt in the beam_data_namelist 
                     ! instead of the time points of beam_data_ufile. 
      
      nubeam_fix_t = 1.0e10  
                     ! JMP
                     ! If time > nubeam_fix_t, nubeam will not be called anymore. 
                     ! All neutral beam related variables (fast ion density, 
                     ! heating profile, driven current, …) are fixed 
                     ! with the values at nubeam_fix_t. 
                     
      nubeam_back_average = 0.0

      nubeam_state_path   = './nubeam_12_nubeam_state.cdf'
      nubeam_xplasma_path = './nubeam_12_xplasma_state.cdf'
      nubeam_profile_path ='./nubeam_12_restart_profs.txt'
      wrt_restart_file_time = -1.e30 
                            ! As explaind above,restart files are overwritten 
                            ! as they are generated.
                            ! Hence at the end of a run the 3 restart files
                            ! (nubeam_state_path,nubeam_xplasma_path,and
                            ! nubeam_profile_path)
                            ! will reflect the situation at the end of the run.
                            ! wrt_restart_file_time is intended to additionally
                            ! save the restart files at a time close  to the time
                            ! requested. (the actual time will be the next 
                            ! call to nubeam made after time .ge. wrt_restart_file_time





      save_nubeam_input =0 ! If set to 1 all input files passed to nubeam will
                           ! be saved. CAUTION, one such file is created each
                           ! time nubeam is called. This could potentially
                           ! lead to a large number of files. These files,
                           ! either with or without the nubeam restart files
                           ! can be used to run nubeam in stand alone mode.
                           ! At this time this option is not available. 

      avg_nubeam_torque = 0 ! This switch is for nubeam only. It is however set in 
                           !inone so that the nubeam namelist can remain undisturbed.
                           !avg_nubeam_torque is an integer that specifies the number
                           !of nubeam calls that are used to form a running average
                           !of the beam torque density returned by nubeam. This is
                           !sometimes required due to the noisy nubeam results.
                           !During steady state operation this averaging is justified.
                           !During beam buildup it is not and you must view the
                           !results with this in mind. Each nubeam call advances the
                           !beam progenated quantites by nubeam_dt sec. (nubeam_dt
                           !is set in the nubeam namelist input). Hence the averaing is over 
                           !avg_nubeam_torque*nubeam_dt sec. Using this option it 
                           !makes sense to have nubeam_dt shorther than usual.

      nubeam_version = 0 ! for example, 201201 (means 2012 Jan version)
                         ! latest version is 201201

c             

      external_beam_cur = 0
      beam_init_restart_file = 0  
      nbi_init = .false.           !flag for trasnp beam initialization         %
      beam_pulse_control = 0                                                   %
      beam_thermal_speed = 0.0d0      
      beam_mode = "run"
      ngauss = 50
      do i=1,ntimbplt                                                                %
        timbplt(i)  = 1.0e6                                                    %
      end do                                                                   %
      do i=1,kb                                                                %
        neg_ion_source(i) = 0        ! default all beam lines to std sources   %
        beamon(i)        = 1.0e18                                              %
        btime(i)         = 0.0                                                  %
        beamoff(i) =0.0
        beam_cycles(i) = 0.0
        beam_end(i) = -1.e18
      end do                                                                   %
      iterate_beam  = .false.                                                  %
      beam_thermal_cutoff =1                                                   %
      therm_frac = 1.0                                                         %
      rtstcx        = 1.0                                                      %
      relaxden      = 1.0            ! no relaxation in beam iterations        %
      relaxden_err  = 0.1            ! default is large,depends on npart       %
      imaxbeamconvg = 5                                                        %
      fast_ion_target = 0            ! neglect fast ions as a beam target      %
      ibcur         = 1                                                        %
      itrapfi       = 0                                                        %
      use_Callen    = 0                                                        %
      ibcx          = 1                                                        %
      ibslow        = 1                                                        %
      nbeamtcx      = 0                                                        %
      fionx         = 0.0                                                      %
      nameb         = 'h'                                                      %
      relnub        = 0.1                                                      %
      tfusbb        = 0.0                                                      %
      iddcal        = 3                                                        %
      fdbeam        = 0.150e-3    ! isotopic content of d in h                 %
c                                   this default is for h beams                %
      d_beam_spin_pol      = 1.0
      he3_thermal_spin_pol = 1.0

      ranseed       = 7**7                                                     %
      npart         = 10000                                                    %
      npart_mcgo    =  3000                                                    %
      hdepsmth      = -1.0   
      fidiff_on     = .FALSE.
      bfr_neutrlz   = .TRUE.      ! fbcur input is before enutralizer          
      freyavb       = 0                                                        %
      time_dep_beam = 0           !default to old single pulse beam model      %
c                                                                              %
c --- smoothing is normally on by above line                                   %
c                                                                              %
      npskip  =  1                                                             %
      iborb   =  1                                                             %
      inubpat =  0                                                             %
      npat(1) = nw                                                             %
      npat(2) = nh                                                             %
      mf      = 41                                                             %
      nbeams  =  1                                                             %
      nsourc  =  2                                                             %
c                                                                              %
c     DIII beam input                                                          %
c                                                                              %
      if (machinei .ne. 'doub-iii')  go to 9020                                %
      naptr = 2                                                                %
      do i=1,kb                                                                %
        anglev(i)    =   0.0                                                   %
        angleh(i)    =  13.5                                                   %
        bvofset(i)   =  39.75                                                  %
        bhofset(i)   =   0.0                                                   %
        bleni(i)     = 553.88                                                  %
        bcur(i)      = 110.0                                                   %
        bptor(i)     =   3.5e6                                                 %
        nbshape(i)   = 'rect'                                                  %
        bheigh(i)    =  10.0                                                   %
        bwidth(i)    =  40.0                                                   %
        bhfoc(i)     = 480.0                                                   %
        bvfoc(i)     = 550.0                                                   %
        bhdiv(i)     =   1.4                                                   %
        bvdiv(i)     =   0.45                                                  %
        ebkev(i)     =  80.0                                                   %
        fbcur(1,i)   =   0.6                                                   %
        fbcur(2,i)   =   0.3                                                   %
        fbcur(3,i)   =   0.1                                                   %
        sfrac1(i)    =   0.5                                                   %
        blenp(i)     = 486.13                                                  %
        rpivot(i)    = 270.0                                                   %
        zpivot(i)    =  89.0                                                   %
        nashape(1,i) = 's-vert'                                                %
        nashape(2,i) = 'b-horiz'                                               %
        aheigh(1,i)  =  10.1                                                   %
        awidth(1,i)  =   0.0                                                   %
        aheigh(2,1)  =   0.0                                                   %
        awidth(2,i)  =  32.0                                                   %
        alen(1,i)    = 442.0 * 1.0026                                          %
        alen(2,i)    = 456.0 * 1.0026                                          %
      end do                                                                   %
      go to 9030                                                               %
c                                                                              %
c     DIII-D beam input.  DEFAULT IS LONG-PULSE-SOURCE SPECIFICATIONS.         %
c                                                                              %
 9020 naptr = 4                                                                %
c                                                                              %
      do i=1,kb                                                                %
        anglev(i)    =   0.0                                                   %
        angleh(i)    =  19.5                                                   %
        bvofset(i)   =   0.0                                                   %
        bhofset(i)   =  42.074                                                 %
        bleni(i)     = 556.808                                                 %
        bcur(i)      = 110.0   ! dont change,see logic in freya

        bptor(i)     =  0.001                                                  
        nbshape(i)   = 'rect-lps'                                              %
        bheigh(i)    =  48.0                                                   %
        bwidth(i)    =  12.0                                                   %
        bhdiv(i)     =   0.50       !degrees                                   %
        bvdiv(i)     =   1.3        !degrees                                   %
        fbcur(1,i)   =   0.7                                                   %
        fbcur(2,i)   =   0.2                                                   %
        fbcur(3,i)   =   0.1                                                   %
        bhfoc(i)     =   1.0d100                                               %
        bvfoc(i)     =   1.0e3                                                 %
        ebkev(i)     =  75.0                                                   %
        ebkev(i)     =   0.0001         ! changed HSJ                          %
        sfrac1(i)    =   0.5                                                   %
        nashape(1,i) = 's-rect'                                                %
        nashape(2,i) = 's-rect'                                                %
        nashape(3,i) = 'b-d3d'                                                 %
        nashape(4,i) = 'b-circ'                                                %
        aheigh(1,i)  =  47.8  
        aheigh(2,i)  =  48.0                                                  
        aheigh(3,i)  =  0.0   ! not set
        aheigh(4,i)  =  0.0   ! not set 
        awidth(1,i)  =  13.8 
        awidth(2,i)  =  17.7 
        awidth(3,i)  =  0.0   ! not set 
        awidth(4,i)  =  50.9
        alen(1,i)    = 186.1                             
        alen(2,i)    = 346.0                                                   %
        alen(3,i)    = 449.0 
        alen(4,i)    = 500.0                                                   %   
        blenp(i)     = 539.0                                                   %
        rpivot(i)    = 286.6                                                   %
        zpivot(i)    =   0.0                                                   %
      end do                                                                   %
c                                                                              %
c  parameters for including rotation in neutral stopping                       %
c                                                                              %
 9030 ds_tk = 5.0                                                              %
      fe_tk = 1.1  !passed to nbsgxn where it is used as ebfac                 %
      ne_tk = 20   !passed to nbsgxn where it is used as nebin                 %
c                                                                              %
c --- following parameters are used in subroutine HEXNB                        %
c                                                                              %
      kdene   =  1                                                             %
      kdeni   =  1                                                             %
      kdenz   =  1                                                             %
      ksvi    =  0                                                             %
      ksvz    =  0                                                             %
      ksve    =  0                                                             %
      krad    =  1                                                             %
      ngh     = 10                                                             %
      ngl     = 10                                                             %
      iexcit  =  5                                                             %
      ilorent =  0                                                             %
      mstate  =  4                                                             %
      ncont   = 30                                                             %
      do 1710 j=1,kprim                                                        %
      znipm(j) = 0.0                                                           %
 1710 atwpm(j) = 0.0                                                           %
      do 1720 j=1,kimp                                                         %
      iz(j)     = 0                                                            %
      izstrp(j) = 0                                                            %
c                                                                              %
c      note izstrp = 0 implies coronal equilibrium for impurity j in hexnb     %
c                                                                              %
      atwim(j) = 0.0                                                           %
 1720 zniim(j) = 0.0                                                           %
*                                                                              @
c ----------------------------------------------------------------------       @
c --- RF HEATING PARAMETERS                                                    @
c ----------------------------------------------------------------------       @
c  PARAMETERS USED BY ALL OR SEVERAL RF MODELS:                                @
c  Parameter krf (=30) gives possible number of RF models to be used.          @
c  nth elements of quantities dimensioned by krf are to be associated          @
c    with each other.                                                          @
c% rfmode(krf)    Flag identifying RF heating model                            @
c    'ech'      : electron cyclotron heating - ray tracing code                @
c    'ich'      : ion cyclotron heating - T.K. Mau model                       @
c    'input'    : input source                                                 @
c    'wedge'    : input wedge-shaped source                                    @
c    'fb'       : give profiles of RF current, or RF power, as input.          @
c                 current drive efficiency according to Fisch-Boozer.          @
c                 (Harvey-Freije)                                              @
c    'fastwave' : 1D or 2D fastwave heating and current drive.                 @
c                 (S.C. Chiu and Bob Harvey)                                   @
c    'raytrace' : LH and fast wave raytracing model                            @
c                 (adapted from Brambilla, S.C. Chiu, T.K. Mau, Bob Harvey)    @
c    'fastcd'   : Kupfer's Ergodic Fast Wave Heating and Current Drive model   @
c    'genray'   : Smirnov-Harvey all-frequency ray tracing.
c
c       NOTE:     Turning models on and off instantaneously may cause
c                 problems in the code since there is no way to uniquely 
c                 determine the on or off state at the transition time.
c                 An easy solution is to
c                 make sure that different models do not start/stop instanaeously
c                 exactly at the same time. Stagger them by .01 msec or something
c                 like that to avoid these problems. 
c% rfon(krf)      Time          (sec) at which RF heating is turned on         @
c% rftime(krf)    Time interval (sec) during which RF heating is on            @
c% turnonp(krf)   Time during which RF power is ramped up to full power,       @
c                 i.e., power is multiplied by                                 @
c                 0.5*(1.0-COS ((time-rfon)/turnonp)).   (default = 0.0)       @
c% turnonc(krf)   Same as turnonp, except for RF current.                      @
c
c
c                 To specify a linear ramp up or ramp down in power rather     @
c                 than the cosine  form use rframp_timeup and rframp_timedown  @
c                 instead of turnonp,turnonc. There is no provision for
c                  rf current turnon in this case. If rframp_timeup and or rframp_timedown
c                 is set then turnonp and turnonc will be ignored by the code.
c                 The initial power at time rfon will always be zero           @
c                 (you can add another model with a square wave power to achieve
c                  any desired offset at ramp up or ramp down time so this
c                  ramp model doesnt allow for offsets)
c  rframp_timeup(krf) the duration of the ramp up
c  rframp_timedown(krf) the duration of the ramp down (sec)
c                 the effect of rframp_timeup,rframp-timedown, is to multiply
c                 rfpow by the time dependent factor tpowrf,defined as         @
c                 tpowrf  = 0.0      for time < rfon
c           for rfon <= time < = rfon+rframp_timeup the ramp up multiplier is
c                 tpowrf   = (time - rfon)/rframp_timeup
c           for rfon + rframp_timeup <= time < = rfon + rftime - rframp_timedown
c                   (ie,the flat top region) the multiplier is
c                  tpowrf = 1.0
c           the ramp down region ( time > rfon + rftime - rframp_timedown <=  rfon+rftime )
c                  tpowrf = (rfon +rfime -time)/rframp_timedown
c           and after the pulse (for time > rfon+rftime)
c                 tpowrf  = 0.0 
c           the flat top section and rampdown can be eliminated by setting
c           rfon = rframp_timeup, rframpdown = 0.0 to get a traingular ramp up that
c           drops back to zero immediately after the ramp is done 
c           and similarly setting rfon = rframp_timedown,rframptimeup =0
c           you get a trianglar ramp down.     
c  rframp_powf(krf)  The final power after ramp down will always be zero       @
c  rup_powi rdown_powf
c  f
c% rfpow(krf)     Total RF power (W)                                           @
c% rnp(krf)       Parallel refractive index                                    @
c% freq(krf)      Wave frequency  (Hz)                                         @
c% irfcur(krf)    current drive switch , changed to floating point number      @
c                 2/12/03 HSJ
c% relrf          Minimum relative change in electron density and              @
c                 temperature profiles between ech calculations                @
c% relrf_pow      FOR ech time dependent power inputs only:
c%                (eg ech_input points to a valid netdcf file)
c%                Minimum relative change in input power which will
c%                force a call to fetch new powers for all gyrotrons.      
c                                                                              @
c  PARAMETERS USED BY INPUT MODEL                                              @
c% rnormin   Normalized radii at which power densities are input;              @
c            may be input in the first namelist instead                        @
c  Set following array elements corresponding to rfmode(k) = 'input'           @
c% njqin(krf)   Number of radial points in input profiles                      @
c% qine(kj,krf) Power density profile (W/cm**3) input to electrons             @
c% qini(kj,krf) Power density profile (W/cm**3) input to ions                  @
c            Note: If rfpow > 0, then qine and qini are normalized such        @
c            that the total source power equals rfpow.                         @
c                                                                              @
c  PARAMETERS USED BY WEDGE MODEL                                              @
c% rfrad1(krf)    Inner minor radius of wedge (cm).                            @
c% rfrad2(krf)    Outer minor radius of wedge (cm).                            @
c            rfrad1,2 are used to set up qini and qine in this subroutine.     @
c            Unfortunately the point at which this is done uses the 1D         @
c            rminor grid to establish the support set for qine,qine.           @
c            The eqdsk r grid can be approximated by setting                   @
c                                rminor = eqdsk rminor value                   @
c                (this value isnt know until you run the code at least once).  @
c% a1rf(krf)      Source magnitude at rfrad1 (arb. units)                      @
c% a2rf(krf)      Source magnitude at rfrad2 (arb. units)                      @
c% wrfe(krf)      Fraction of power deposited in electrons.                    @
c% wrfi(krf)      Fraction of power deposited in ions.                         @
c            Note: the source profile is normalized such that the total        @
c            input power equals rfpow.                                         @
c            Note:  The wedge model is an alternate way of setting the         @
c                   qine(i,k) and qini(i,k), variables, i = 1,nj.              @
c                   Set rfmode(k) = 'wedge'.                                   @
c            Note:  If irfcur(model) .ne. 0 then subroutine RFCUR_MODEL        @
c                      will be called to calculate an approximate RF           @
c                      driven current. The current drive model in this         @
c                      subroutine is only a crude approximation so             @
c                      do not use it indiscriminately HSJ!!!!!!!!!!!!!!!!!!    @
c            Note:  more than one wedge model can be active simultaneously     @
c                     thus a more realistic energy deposition can be achieved  @
c                     by combining multiple wedges HSJ                         @
c                                                                              @
c  PARAMETERS USED BY FISCH-BOOZER MODULE.                                     @
c% rfon,rftime,turnonp,turnonc,rfpow,rfmode,relrf, 
c%  rframp_timeup,rframp_timedown ,  basically as above.                       @
c% rfmode = 'fb' for fisch-boozer model, with cordey et al efficiencies.       @
c% lmode = 0, for lower hybrid efficiency; harmonic number (1 to 3) for ech.   @
c% irfcur = 1.0/0.0/-1.0 for co/no/counter                                           @
c% itrapech     switch for modifying ech current drive by parabola             @
c               suggested by v. chan.  if itrapech = 1 and RF cur drive is     @
c               on then the calculated RF current is multiplied by             @
c               (1.0-eps**0.7)**0.73 where eps is inverse aspect ratio.        @
c                the powers .7 and .73 are a crude fit to the curve            @
c                fig.2 in cohen 'effects of trapped electrons on               @
c                current drive',to be published in phys of fluids.             @
c% rnp(krf) =  positive:                                                       @
c      parallel refractive index.  used with te and ene to find efficiency     @
c      according to the Fisch-Boozer v-parallel current drive formula.         @
c         negative:                                                            @
c        absolute value is a multiplier of SQRT (2*te/me) giving parallel      @
c        phase velocity of wave in terms of the local thermal velocity.        @
c% ifb= 1, specify total RF power and profile factors, and calculate           @
c          RF power density and current profiles and total RF current.         @
c       2, specify total RF current and profile factors, and calculate         @
c          RF current and power density profiles, and total RF power.          @
c% ifbprof= 1,parabolic profile.                                               @
c           (1.0-(r**2/a**2)**alphaf(2))**alphaf(1)                            @
c         2,guassian profile.                                                  @
c           EXP (-(r/(alphaf(1)*a))**2  - boundary value.                      @
c         3,shifted guassian profile.                                          @
c           EXP (-((r-alphaf(2)*a)**2/(alphaf(1)*a)**2) - boundary value,      @
c           except, constant inside r = alphaf(3)*a at value of above          @
c                   function.                                                  @
c% alphaf= as above.                                                           @
c% rfcur= total current, in case current to be specified.                      @
c% rfpow= total power, in case power profile to be specified.                  @
c                                                                              @
c  PARAMETERS USED BY ICH MODELS:                                              @
c% freq(krf) Frequency (1/sec) of applied wave                                 @
c% iside     +1 for high field launch; -1 for low field launch                 @
c% navg      apply moving average to output radial power absorbtion            @
c            profiles by averaging over navg adjacent points                   @
c            (total of 2*navg+1 points in average).                            @
c                                                                              @
c  PARAMETERS USED BY ECH (RAY-TRACING) MODEL ONLY:                            @
c  ech_input This switch gives a netcdf file name that is to be used to 
c     read in ech data. Any ech data in inone WILL BE IGNORED IF THIS SWITCH
c             IS SET !


c% irfplt    If irfplt .ne. 0, then plots are made each time TORAY is called.  @
c            Otherwise, the output files from the last call to TORAY can       @
c            be processed by manually by running XPLOT after the completion    @
c            of ONETWO.  Note:  if multiple TORAYs are run only the last       @
c            TORAY run remains, not all of the TORAYs at the last time.        @
c  The logic implementing the following definitions of IRFPLT and TIMRFP       @
c  was found to have a bug (TCL 6/17/91).  It has been disabled (CRAY40)       @
c  and the above definition is implemented in versions later than 6/10/91.     @
c                                                                              @
c% irfplt    If .ne. 0, plotted o/p from the ray tracing code is enabled.      @
c            Plots are made at the first call to ech, and thereafter           @
c              whenever the central temperature changes by more that           @
c              20%, unless timrfp(1) is < 1.0e-6.                              @
c            Plots are only made at times timrfp(i),i = 1,5, for each          @
c              time which is less than 1.0e6.                                  @
c% timrfp(5) plotted output from ray tracing is given at these times,          @
c            when set to less than the default value (1.0e6).                  @
c% gafsep    For codeid .ne. 'onedee', then this is the fraction               @
c            of poloidal flux smoothed over at the plasma                      @
c            magnetic axis, and sets the mesh spacing for passage              @
c            of output profiles from TORAY to ONETWO, i.e. 1/gafsep points.    @
c            If  = 0.0 and/or codeid = 'onedee', then the anayltic,            @
c            circular, concentric flux surface model is used for the           @
c            ray tracing.                                                      @
c% freq(krf) Frequency (1/sec) of applied wave                                 @
c% idamp = 0 no damping calculation                                            @
c  idamp = 1 near first  harmonics, using Matsuda & Hsu's routine (dampga)     @
c  idamp = 2 near second harmonics, using dampga                               @
c  idamp = 3 near third  harmonics, using damprm                               @
c  idamp = 7 good for between 2nd and 3rd, using damprm                        @
c  idamp = 8 near first harmonics, using Mazzucato's routine                   @
c           for fully relativistic damping.                                    @
c           Fidone, and Granata code.(Physics of Fluids, p3745, 1988)).        @
c           (Added by Bob H. 7-26-88; see helpme in l??ech??                   @
c  idamp = negative of above numbers, then damping is as for the               @
c          corresponding positive number, but the EQDSK is not                 @
c          reprocessed with GAFIT to form the PSICOF coefficient file.         @
c          This saves computer time and is valid if the equilibrium has        @
c          not changed. (the file involved is "mhddat".)                       @
c            within the ONETWO libfile).                                       @
c% xec(krf)   Major radius of source.                                          @
c% zec(krf)   Height above midplane (cm) of source.                            @
c% wrfo(krf)  Fraction of RF power launched as O-mode.                         @
c             Note that if wrfo > 0.95, X-mode power will be neglected,        @
c                   and if wrfo < 0.05, O-mode power will be neglected.        @
c                                                                              @
c% necsmth    smoothing control parameter:                                     @
c               necsmth = 0 (the default)                                      @
c                  = => no smoothing                                           @
c               necsmth = j  (a positive integer)                              @
c                  = => apply j passes of a boxcar smoother                    @
c                     to the eccd current and ech heating profile.             @
c                     Note that a large j will flatten out the profiles.       @
c                     We probably don't want to use j > 3                      @
c                                                                              @
c% hlwec     angular dispersion of rays ??? degrees ??                         @
c% ratwec    specifies antenna shape:                                          @
c              ratwec = 1  =>  circular antenna shape                          @
c              ratwec > 1  =>  more elongated   shape ?                        @
c                                 
c                                                                              @
c              New Input For Toray, 01/17/03  -- HSJ
c
c         There are three ways to spawn a version of Toray.


c         The first is to just rely on the default settings in
c         Onetwo (see default setting of the switches). This way
c         should almost always be what the user wants. It runs the 
c         latest version of toray available on the machine that 
c         Onetwo is executing on. To see what toray versions are known
c         to Onetwo you can browse the file ext_prog_info.f90

c         The second is to specify that a particular version
c         of Toray is to be run,using switch toray_version. Onetwo will
c         search its Toray paths on the architecture it is running on 
c         and if an appropriate version with the right mhd and transport grid
c         sizes is found, it will be used. If this fails the user
c         will be informed and the code will quit.( It is becoming quite
c         challenging to keep up with all Toray versions, grid sizes, 
c         and architectures. Hence users that require some special version/
c         architecture may specifically have to request it).

c         The third way is to specify a path to a custom built Toray.
c         (This allows the latest  versions of Toray to be used before
c         they become public for example). Set toray_path to the complete
c         executable name. For example( on Hydra):
c         /u/stjohn/toray/toray/129_51/hp/toray
c         This will force onetwo to use the executable called toray in
c         directory /u/stjohn/toray/toray/129_51/hp
c         For your custom toray you must be sure that the Onetwo and
c         toray grid sizes are commensurate. If they are not Toray (or
c         Gafit) will complain and return an error code to Onetwo.
c
c         Toray_path,(see below) if set, takes precedence. If not set 
c         (in inone)
c         then toray_version is checked. If this is not set then
c         the default mode is used. Issues involved with tdem mode
c         of operation of Onetwo are handled automatically if the user
c         has choosen to run tdem mode.
c

c toray_version  default is version 1.41 . Specify a
c         number . ge. 1.4 for new Toray f90 version with an internal
c         gafit . Specify a 
c         number less than 1.4 to get the old version of toray with
c         gafit run as a separate program.

c toray_path  The default path is set internally to
c         point at the selected version of toray. If you  want to run a
c         specific version then you can set the path to that version
c         here. 
c         toray_path is a 256 (or less) character variable.
c         toray_path takes precendence over toray_version. But
c         if toray_path is set then you must also specify a version
c         number of toray that toray_path is pointing to. In this case
c         the significance of the version number lies in the fact
c         that for version numbers > 1.3 gafit is assumed internal to
c         toray (and hence requires a different interface in onetwo)
c         For version < 1.3 gafit is external to Toray and the old
c         interface musr be used. 

c echin_save  integer, either 0 or 1, default 0. 
c         Toray can be run in
c         standalone mode using just the files echin and mhddat. Hence
c         if the files echin and mhddat  are saved as Onetwo executes it is
c         possible to go back and rerun toray (perhaps changing 
c         something in the echin file). echin\_save  =0 retains only the
c         last version of echin,mhddat created by \ot. echin\_save =1 saves
c         all versions of these files indexed by a time stamp. Notice
c         that if you have a time sequence of eqdsks processed using
c         TDEM  mode then this option is a convenient way to get time
c         interpolated mhddat files. Note however that mhddat is a
c         binary file which means its non portable across
c         machines.( A netcdf file would make more sense here I think
c         but thats a Toray issue). Finally recall that
c         Toray can be run using only the input files echin and mhddat
c         only if toray.in has igafit =0 (which is the default in the
c         Toray v1.4  code).



c  PARAMETERS USED BY ICH MODEL ONLY:                                          @
c% xkpar     k-parallel (1/cm)                                                 @
c% nhigh     number of terms to retain in K                                    @
c  Following two varialbes dimensioned, in view of WEDGE model, above.         @
c  Set values appropriate to rfmode(k).                                        @
c% rfrad1ic    Inner major radius (cm) for RF grid points                      @
c% rfrad2ic    Outer major radius (cm) for RF grid points                      @
c            If rfrad1ic,rfrad2ic = 0 (default) they are calculated as the     @
c            intersection of the horizontal chord at ylaunch with:             @
c            for onedee case- the outer flux surface;                          @
c            for other cases- the limiter surface.                             @
c% nrfzon    Number of zones in the RF grid.  (Each zone may have              @
c            a different density of grid points.)                              @
c% nprf(1)   Should always = 1                                                 @
c  nprf(i)   (i = 2,nrfzon+1)  Point number corresponding to the               @
c            location of rfzone(i).  Thus nprf(nrfzon+1) = total number        @
c            of RF grid points, and there are nprf(i+1)-nprf(i) steps          @
c            in the zone between rfzone(i) and rfzone(i+1).                    @
c% rfzone(i) Location of zone boundaries for RF grid, normalized so            @
c            that -1 corresponds to rfrad1ic, and +1 corresponds to rfrad2ic.  @
c            Note that the RF grid does not have to stretch from rfrad1ic      @
c            to rfrad2ic  (e.g. rfzone(1) may be >-1.0).                       @
c% relrf     Minimum relative change in electron density and                   @
c            temperature and ion temperature profiles between                  @
c            ich calculations.                                                 @
c% ichmod    1 to use t.k.mau code only                                        @
c            2 to use s.c.chiu code only                                       @
c            3 to use t.k.mau code when beta-ion(0) < betalm,                  @
c            and use s.c.chiu code when beta-ion(0) >= betalm                  @
c% betalm    Limit for beta-ion(0) when ichmod = 3                             @
c% ykperp    k-poloidal (not implemented)                                      @
c% ylaunch   Vertical launch point for wave (cm above midplane)                @
c            (Not implemented.)                                                @
c            (Note that midplane differs for onedee and 2d cases.)             @
c                                                                              @
c  PARAMETERS USED BY FASTWAVE MODULE                                          @
c% codeid, rfmode, rfon, rftime, turnonp, turnonc,
c%  rframp_timeup,rframp_timedown , rfpow, irfcur, relrf        @
c% nzrffw     number of fw energy channels ( .le. kzrf = 11)                   @
c% zrffw(kzrf)   centerline height of the fw energy channel (cm)               @
c% freqfw(kzrf)  frequency of each energy channel                              @
c% rnpfw(kzrf)   parallel refractive index of each energy channel              @
c% pzrffw(kzrf)  fraction of rfpow into fw energy channel                      @
c% htsfw(kzrf)   vertical extent of RF channel at outboard, rmajor, and        @
c             inboard positions (cm).  (Same values used for each energy       @
c             channel)                                                         @
c% iswchfw    1, ttmp+landau damping (use for ICRF).                           @
c             2, landau damping (use for LHRF).                                @
c% nspfw      number of primary ion species on which there is ion damping.     @
c             (nspfw was eliminated by s.c.chiu on 21 feb 92)                  @
c% lifw(kprim)  starting harmonic number for damping calculations.             @
c% nihfw(kprim) number of harmonics calculated.                                @
c% impath      multiple absorption switch, =1 for single path absorption.      @
c                                          =2 for multiple path                @
c                                                                              @
c
c
c
c
c
c  --------------------------------------------------------------------------  @
c  PARAMETERS USED BY LH AND FASTWAVE, RAYTRACE MODULE (CURRAY):               @
c  (rfmode = 'raytrace')                                                       @  
c  --------------------------------------------------------------------------  @
c
c

c  codeid, rfmode, rfon, rftime, freq, irfcur, relrf,rfpow                     @
c  note that rfmode,rfon,rftime,freq,irfcur,rfpow are common to all rf models  @
c  and hence are indexed by their position in rfmode. However, quantities     @
c  which are specific to curray are indexed by counting only the occurrences   @
c  of 'raytrace' in rfmodel. For example if rfmodel(3),(7),and (11) are set
c  to 'raytrace' (even if they are turned off by rftime) then
c  powersrt(1:irt,1) =====> rfmodel(3)   ======>  1st curray model specified
c  powersrt(1:irt,2) =====> rfmodel(7)   ======>  2cd curray model specified
c  powersrt(1:irt,3) =====> rfmodel(11)  ======>  3rd curray model specified
c  freq(3) ===============> rfmodel(3)   ======>  1st curray model specified
c  freq(7) ===============> rfmodel(7)   ======>  2cd curray model specified
c  freq(11) ==============> rfmodel(11)  ======>  3rd curray model specified
c             etc.
c
c  
c
c  curray_path  character variable (of length 256) normally not used. But if   @
c               you have a special version of curray you want to run instead   @
c               the default version that Onetwo knows about then you can set   @
c               the fully qualified executeable name here.                     @relaxden
c
c save_curray_input  = 1 if the files curray_in,trxpl.out and the eqdsk        @
c               are to be saved for later use in running curray in stand       @
c               alone mode. Default is  save_curray_input  = 0, which          @
c               overwrites these files each time curray is called. Note that   @
c               a potentially large number of files could be created using     @
c               this switch so beware.                                         @
c  curray_fi    used to renomralize the fraction of heating power given to     @
c               electrons and ions. In any bin the total power is              @
c               pwt =pwe +pwi . Onetwo will take pwi = curray_fi *pwt and      @
c               pwe = pwt - pwi if curray_fi is not zero.                      @
c  psistep  stepsize in rho for ray advancing                                  @
c                                            written  directly to curray_in    @
c  pkexpnt  parameter for adjusting psistep                                    @
c                                            written  directly to curray_in    @
c  epserr   first error limit for evaluating dipsersion relation               @
c                                            written  directly to curray_in    @
c  epser1   second error limit ( = 0.2 * epserr) default = 0.004               @       
c                                            written  directly to curray_in    @
c  igraph         =0    suppress graphical output                              @
c                 =1    normal graphical output                                @
c                                            written  directly to curray_in    @
c  iprint         =0    output suppressed                                      @
c                 =1    normal output                                          @
c                 =2    more detailed output                                   @
c                 =3    determine output at each call of OUTPUT1(x)            @
c                 =-1   condensed information on ray tracing                   @
c                                            written  directly to curray_in    @
c  idcur          =1 :    use analytical formula for j/p                       @
c                 =2 :    not working                                          @
c                 =3 :    use 'ADJ' for j/p  (Use for high-beta NSTX)          @
c  indvar         =0 :  poloidal phase is used as independent variable in      @
c                       ray tracing                                            @
c                 =1 :  total phase used as independent variable (recommended) @
c  ichoisrt(i)         =1    high frequency dispersion used in wave propagation     @
c                    =2    no frequency limitation  (default)                     @
c  ichois               for backward compatibility ichois can also
c                       be used (instead of ichoisrt) to input just a single
c                       case. If it is set then ichoisrt(1) = ichois
c                       is forced, even if ichoisrt(1) is also set !!!!!!!!!
c                       can have fast and slow wave calls to curray
c  igrill         =-1   analytically fitted spectrum                           @
c                 =-3   read in toroidal and poloidal spectrum                 @
c                       (nnkpar,nnkpol etc. -3 recommended)                    @
c  modcd          =0    use Ehst-Karney model for analytic j/p calculation     @
c                 =1    use Chiu-Karney-Mau model                              @
c                       (for analytic j/p, set icurdr=0, idcur=1)              @
c  idmpswrt(i)          specify idmpsw for raytrace model i                    @
c                       if idmpsw is set then idmpswrt(1) = idmpsw is forced   @
c  idmpsw         =0    ion damping not calculated                             @
c                 =1    magnetized ion damping calculated (recommended)        @
c                 =2    unmagnetized ion damping (use with caution)            @
c  nminor               number of minority species                             @
c                       (note: nspec=nprim+nimp=nsp1+nminor, the use           @
c                       of nminor is mainly in ion-ion-hybrid                  @
c                       resonance when the cold dielectrics have               @
c                       singularities due to minority fundamental              @
c                       resonance. To avoid that in ray-marching,              @
c                       set nminor=1 and then the ray-tracing                  @
c                       calculation ignores the last ion-species;              @
c                       the absorption and polarization still include          @
c                       this species.)                                         @
c  kalfa          =0    ignore fast ion damping from slowing down distribution @
c                        (forced default)                                      @
c                 =1     calculate unmagnetized fast ion damping               @
c                        kalfa =1  is meant for a slowing down                 @
c                        distribution and should not be used at this time.     @
c                        Use only kalfa =0. The fast ions are treated as a     @
c                        a separate species with maxwellian temperature when   @
c                        kalfa =0.                                             @
c  beam_spec(4)         Integer input determines fast ions contribution when   @
c                       kalfa = 0.
c                    set beam_spec(1)  = 1 to use full beam energy componenet  @
c                    set beam_spec(2)  = 1 to use 1/2  beam energy componenet  @
c                    set beam_spec(3)  = 1 to use 1/3  beam energy componenet  @
c                    set beam_spec(4)  = 1 to include fusion alpha componenet  @
c                    Components that are included are treated as local         @
c                    Maxwellian distributions in Curray. The effective temp    @
c                    of each maxwellian is obtained from the stored energy     @
c                    density:
c                              Ti = (2/3)*wbeam(r)                             @
c                                 or (2/3)*walpha(r)                           @
c                    The curray input file trxpl.out will show these           @
c                    contributions as individual species, for each beam line   @
c                    and each beam energy component selected, as well as the   @
c                    fusion alpha componenet (if slected)                      @
c                    If treatment of individual components of the beam is not  @
c                    desired then set beam_spec(1) = -1. this will cause       @
c                    the partial pressures of all beam components and injectors@
c                    to be combined into a single effective beam species.      @
c                    Note that beam_spec(1) = -1 sums the components for each  @
c                    beam but separate beams are not combined.
c                    Fast alphas are still governed by the setting of          @
c                    beam-spec(4) in this case.                                @
c    tmin_curray    KEV, the minimum temperatrue for fast ions
c    enbmin_curray  #/cm**3 the minimum density for fast ions in curray
c                   tmin_curray and enbmin_curray are FUDGE FACTORS
c                   required to make  curray run in some cases where the
c                   fast ion density doesnt cover the whole rangfe of rho
c                   from [0,1]. In this case these miniml values are
c                   assigned.
c
c% nrayptrt      maximum number of steps along each ray (step size             @
c                set elsewhere (in a raytracing input file).                   @
c                (nraypts is the curray name for nrayptrt)                     @
c  nspect        number of toroidal wave number bins                           @
c                in onetwo this number is fixed at parameter krf . The actual  @
c                number of bins in use is controlled by specifiying the        @
c                power in each bin with powersrt . A zero value in powersrt    @
c                effectively eliminates that bin.                              @
c                The bins must be contiguous !!!                               @
c
c  kcrm          kcrm is a parameter set in param.f90 it is the maximum        @
c                number of curray models that can be declared in               @
c                rfmode. (even if the model is turned off). Note that kcrm is  @
c                not an input quantitity. The number of curray models is       @
c                obtained internally by counting the number of cases that have @
c                rfmode(i) = 'raytrace'                                        @
c                
c% powersrt(krt,kcrm) power (Mwatt) of each n parallel interval.               @
c                The sum over all bins  is renormalized to rfpow(model)        @
c                The curray name is powers(nspect)                             @
c% nnkparrt(krt,kcrm) number of rays representing an n_parallel interval.      @
c% anzinfrt(krt,kcrm) lower (in absolute value) n_parallel of interval.        @
c% anzsuprt(krt,kcrm) upper (in absolute value) n_parallel of interval.        @
c% nthinrt(kcrm)       number of poloidal locations along antenna,             @
c                      represented by rays.                                    @
c% nnkpolrt(krt,kcrm) number of n-poloidal's in krt-th bin                     @
c% anpinfrt(krt,kcrm) lower limit of poloidal refractive index                 @
c% anpsuprt(krt,kcrm) upper limit of poloidal refractive                       @

c               [note: The total number of rays is
c                sum(1 to nspect) of nnkpar(j)*nnkpol(j)*nthin ]
c% nfwsmth       smoothing parameter for fast waves; meaning the same          @
c                as necsmth for ech (0 means     no smoothing,                 @
c                             n .ne. 0 means n-path smoothing)                 @
c  bmaxrt         integer,maximum Bessel function order                        @
c% heightrt(1...kcrm)      full height of the antennas (cms).                  @
c% maxrefrt(1...kcrm)      Maximum number of reflections from plasma edge.     @
c%                (1 for single pass for LH, 2 for single pass for FW).        @
c% islofart(1...kcrm) =1, for launch in LH mode, -1 for lauch in FW mode.    @
c% irayiort     =0    maximum data output to 'rayop' file, with absorption     @
c%                    and current profiles plus ray data.                      @
c%                =1    only ray data with no absorption and current profile   @
c%                    to 'rayop'.                                              @
c% incrt          number of harmonics in ion-damping is 2*incrt +1             @
c% icurdrrt     =0      adjoint Green's function not calculated                @
c%               =1      adjoint Green's function calculated                   @

c%  nalp_thresh                                                                @
c%  nb_thresh         real numbers, the density threshold for                  @
c%                    fusion alphas and beam densitites below which            @
c%                    the contribution to the ion damping is neglected.        @
c%                    The density must be less than this value across the      @
c%                    entire grid before the contribution is ignored           @
c%

c% nicsmth    smoothing control parameter( same as nfwsmth):                                     @
c               nicsmth = 0 (the default)                                      @
c                  = => no smoothing                                           @
c               nicsmth = j  (a positive integer)                              @
c                  = => apply j passes of a boxcar smoother                    @
c                     to the iccd current and ich heating profile.             @
c                     Note that a large j will flatten out the profiles.       @
c                     We probably don't want to use j > 3                      @



c  .......       There are other seldom changed input variables in             @
c                the input (command) file to the raytracing code.              @
c                      
c ----  additional input for curray_in -------------------------------         @
c ----  curray_in is one of  the input files  for program curray.              @
c ----  Curray_in is read by curray in stand alone mode as well as when it     @
c ----  is called from Onetwo. When called from Onetwo, subroutine             @
c ----  wrt_curray_in is used to write the file curray_in from the             @
c ----  information supplied in the second namelist of inone.                  @
c ----  The useer contolled infromation for curray_in is described above.
c ----  Because curray has no internal defaults the curray_in file must
c ----  set all input parameters. In addition to the above user supplied
c ----  inputs the follwoing are set in the curray_in file written by 
c ----  Onetwo.
c ----  The following parameters are written to curray_in without user         @
c ----  intervention, see wrt_curray_in.                                       @
c
c ----  iprof =1    means profiles for te,etc, are to be read from file        @
c                   trxpl.out                                                  @
c ----  idrive =0,1 current drive switch, idrive is set to irfcur(model)       @
c ----              where model points to a curray case                        @
c ----              Note that irfcur(model) alows fractional current drive     @
c ----              so the convention is that whenever irfcur(model) .ne. 0.0  @
c ----              idrive is set to 1, meaning that curray calculates the     @
c ----              current drive. In Onetwo the cd profile returned by        @
c ----              curray will be multiplied by irfcur(model) to allow        @
c ----              more or less current to be included in the Onetwo          @
c ----              calculattions.                                             @
c ---- ioread = 0   set to not  read file raytrin. Raytrin is no longer        @
c ----              written by onetwo.                                         @
c ---- nprim        number of primary ion species, same as Onetwo input        @
c ----              nprim.                                                     @
c ---- atm(i)       ion species mass numbers obtained from                     @
c ----              Onetwo variable atw(i)                                     @
c ---- azi(i)       ion species charge numbers obtained from Onetwo variable   @
c ----              z(rho=0,ion)
c ---- kboot        bootstrap current, set to zero to turn off curray bootstrap@
c ----              calculations. (bootstrap current is calcuated in Onetwo    @
c ----               instead)                                                  @
c ---- atmf         fast ion mass number (only 1fast ion specieas is allowed   @
c ---- azf          fast ion charge number                                     @
c ---- ebkev        fast ion species energy,kev same as Onetwo input ebkev     @ 
c
c                   Only one fast species is allowed at present                @
c
c ---- nalfa   =1, number fast ion species fast ion species                    @
c ---- freqcy        wave frequency (Hz) set from Onetwo variable freq(model)  @
c ---- nraypts      set to onetwo variable nrayptrt                            @
c
c

c ----       The following variables can be present in curray_in. They describe@
c ----       analytic profiles input. In Onetwo such analyic profiels are not  @
c ----       passed to curray. hence these parameters are simply skipped:      @
c ----       See curray_input file in curray repository for description if you @
c ----       want to use thes in a stand alone curray run:                     @

c ---- denc
c ---- tempec
c ---- aconc
c ---- tite
c ---- denedg
c ---- temedg
c ---- frcne
c ---- gamte
c ---- alpte
c ---- gamte1
c ---- alpte1
c ---- gamene
c ---- alpene
c ---- gamene1
c ---- alpene1
c
c ----  end curray input ------------------------------HSJ-05/22/03-----------

c
c
c  PARAMETERS USED BY FASTCD MODULE                                            @
c  (rfmode = 'fastcd')                                                         @
c  codeid, rfmode, rfon, rftime, rfpow, turnonp, turnonc,                      @
c%  rframp_timeup,rframp_timedown , irfcur, relrf,                             @
c% nkfcd(krf)       number of toroidal modes used in the calculation           @
c% xntor(kzrf,krf)  toroidal mode numbers of the antenna spectrum              @
c                   xntor = (n_parallel * omega * rmajor/c)                    @
c                   c = 3.0e8 m/sec, omega = freq/2pi                          @
c% rpant(kzrf,krf)  relative power on each toroidal mode number,               @
c                   in arbitrary units                                         @
c% gamloss(krf)     phenomenological parasitic damping rate                    @
c  fastcd_path     character*256  give full executable path to
c                   user preferred fastcd code. Normally you do not set
c                   his value since Onetwo knows the path (it is 
c                   architecture dependent and is defined in ext_prog_info.f90)
c                   If you want to test a special version then you can do 
c                   so with this switch setting.

c
c
c  INONE NAMELIST PARAMETERS USED BY THE GENRAY RF MODULE                      @
c  (rfmode() = 'genray')                                                       @
c  genray_path()    a 256 character (or less) path to the genray executable    @
c
c  [BH, 050330, save_genray_io not fully implemented:]
c  save_genray_io  = 1 if the files genray_profs_in.txt and the eqdsk          @
c               are to be saved for later use in running genray in stand       @
c               alone mode. Default is save_genray_io  = 0, which              @
c               overwrites these files each time genray is called. Note that   @
c               a potentially large number of files could be created using     @
c               this switch so beware.                                         @

c  rfmode, rfon, rftime, rfpow, turnonp, turnonc, rframp_timeup,               @
c    rframp_timedown, irfcur, relrf, genray_fi (defns as above, except         @
c    genray_fi has same meaning for genray as curray_fi has for curray).       @
c    default genray_fi=0.0.                                                    @
c    MAX(nfwsmth,nicsmth) smoothing applied to rf profiles, as above for       @
c    curray. If above MAX = 0, the necsmth applied as above.  Defaults are     @
c    no smothing, nfwsmth=nicsmth=necsmth=0.                                   @
c    
c  All other variables required for specification of the rays (frequency,      @
c    mode, ray configuration, launcher locations, etc.) are input              @
c    through genray.dat or genray.in.   The namelist file genray.dat is        @
c    the old (Aug '07) input format (mixed cgs and mks units) whereas          @
c    genray.in uses mksa/keV units. If genray.dat is not given, then           @
c    genray.in is expected.                                                    @
c    The rfpower setting in the genray.dat/genray.in file is scaled to the     @
c    inone namelist value rfpow.                                               @
c
c  genraydat(1:krf) is array of character*128 giving path genray.dat/genray.in @
c    input filenames, where the index k corresponds rfmode(k)="genray".        @
c    Different paths can be used for each rfmode(k)="genray", enabling         @
c    application of genray for different launch configurations.                @
c    The default value for all genraydat='genray.dat'.                         @
c
c  Annotated versions ('template' files in the genray distribution)            @
c    of the genray.dat/genray.in file are available, and also                  @
c    there are detailed explanations in the Genray Manual                      @
c    (see CompX report  Genray_manual_CompX-2001-1-V2.pdf, and updates).       @
c
c
c                                                                              @
c  PARAMETERS USED TO SUPPLY EXTERNAL RF CURRENT FROM INONE FILE               @
c% extcurrf(l)   master switch, if = 0.0 then model l is turned off and        @
c                following parameters are not used. If .ne. 0.0 then           @
c                extcurrf*extcurrf_curr is added to the currf array.           @
c                (The currf array contains the sum total of all RF current     @
c                sources. No time dependence of extcurrf_curr                  @
c                current is available. It is always on or always off.)         @
c% extcurrf_amps(l) if equal to zero then you get a total RF current           @
c                given by the integral of extcurrf*extcurrf_curr               @
c                if extcurrf_amps is .ne. 0.0 then the product extcurrf*       @
c                extcurrf_curr is renormalized to give extcurrf_amps total     @
c                current. (in this case the value of extcurrf is irrelevant    @
c                so long as it is not 0.0)                                     @
c% extcurrf_id(l)   a character variable (must be less than 80 characters)     @
c                used to identify extcurrf_curr in output                      @
c% extcurrf_nj(l)   an integer which gives number of values in extcurrf_curr.  @
c                if extcurrf_nj equals nj then extcurrf_rho is assumed equal   @
c                to the internal rho (so extcurrf_rho need not be given).      @
c                If extcurrf_nj is less than nj then it is assumed that        @
c                extcurrf_curr                                                 @
c                is a spline representation and exactly extcurrf_nj values     @
c                of NORMALIZED extcurrf_rho must be given in inone,with the    @
c                first value exactly 0.0 and the last value exactly 1.0        @
c                The number of values in the spline representation must be     @
c                .le. kpsi (values of extcurrf_nj less than 33 are always ok,  @
c                see include file mhdpar.i for exact value of kpsi )           @
c% extcurrf_rho(j,l)   j=1,2..extcurrf_nj NORMALIZED rho values for spline     @
c                representation of extcurrf_curr.(zero gradient at rho=0       @
c                and natural spline at plasma edge are used)                   @
c% extcurrf_curr(j,l) j=1,2...extcurrf_nj  The RF current in amps/cm**2        @
c                that is to be added to currf. (this will be the total         @
c                RF current if no RF models are turned on, or will be an       @
c                additional RF current if other RF models are in use as well)  @
c                See subroutine GET_EXTERNAL_RFCUR for further details.        @
c  NOTE          the subscript l ranges from 1 to krf and corresponds to       @
c                the rf models given above. Thus an external current and       @
c                electron, ion heating components can be input for each model @
c                                                                              @
c  AN IDENTICAL SET OF INPUTS FOR (RF) HEATING OF ELECTRONS AND IONS IS        @
c  ALSO AVAILABLE. THE RELEVANT INPUT PARAMETERS ARE                           @
c  for electrons:                                                              @
c  variable            analogous meaning as variable above for currf           @
c  ------------        ---------------------------------------------           @
c% extqerf                     extcurrf                                        @
c% extqerf_rho                 extcurrf_rho                                    @
c% extqerf_qe                  extcurrf_curr                                   @
c% extqerf_nj                  extcurrf_nj                                     @
c% extqerf_id                  extcurrf_id                                     @
c% extqerf_watts               extcurrf_amps                                   @
c                                                                              @
c  extqerf_qe is input in watts/(cm**3 ) (you can input this in keV/cm**3 sec  @
c  by setting the multiplier extqerf to 1.602e-16 joules/keV)                  @
c  extqerf_watts is input in watts                                             @
c                                                                              @
c  for ions:                                                                   @
c  variable            analogous meaning as variable above for currf           @
c  ------------        ---------------------------------------------           @
c  extqirf                     extcurrf                                        @
c  extqirf_rho                 extcurrf_rho                                    @
c  extqirf_qi                  extcurrf_curr                                   @
c  extqirf_nj                  extcurrf_nj                                     @
c  extqirf_id                  extcurrf_id                                     @
c  extqirf_watts               extcurrf_amps                                   @
c                                                                              @
c  extqirf_qe is input in watts/(cm**3 ) (you can input this in keV/cm**3 sec  @
c  by setting the multiplier extqirf to 1.602e-16 joules/keV)                  @
c  extqirf_watts is input in watts                                             @
c                                                                              @
c ----------------------------------------------------------------------       @
      toray_path = ''          !default must be zero length string
      echin_save = 0
      toray_version = 0.0      !code will use latest version if not set
      call set_rf_data
      do k=1,krf                                                               %
        rfmode(k)  =  no_rf                                                    %
        rfon(k)    =  0.5_DP*HUGE(1._DP)                                                 %
        rftime(k)  =  0.0                                                      %
        turnonp(k) =  0.0                                                      %
        turnonc(k) =  0.0                                                      %
        rframp_timeup(k) = -1.0                                                %
        rframp_timedown(k) = -1.0                                              %
        rfpow(k)   =  0.0                                                      %
        rnp(k)     =  1.0                                                      %
        freq(k)    = 60.0e9                                                    %
        irfcur(k)  =  0.0                                                      %
        xec(k)     = rmajor - rminor                                           %
        zec(k)     =  0.0                                                      %
        wrfo(k)    =  0.0                                                      %
        njqin(k)   =  0                                                        %
        a1rf(k)    =  0.0                                                      %
        a2rf(k)    =  0.0                                                      %
        wrfe(k)    =  0.0                                                      %
        wrfi(k)    =  0.0                                                      %
        rfrad1(k)  =  0.0                                                      %
        rfrad2(k)  =  0.0                                                      %
      end do   



c    curray defaults:
         save_curray_input  =  0                                               %
         nalp_thresh = 1.e9        ! #/cm**3                                   %
         nb_thresh =   1.e9        ! #/cm**3                                   %
         curray_path = ''          !default must be zero length string         %
         psistep = 0.025                                                       %
         pkexpnt = 0.001                                                       %
         epserr = 0.02                                                         %
         epser1 = 0.004                                                        %
         igraph = 0                                                            %
         iprint = 0                                                            %
         idcur  = 1                                                            %
         incrt  = 1                                                            %
         icurdrrt = 0                                                          %
         indvar = 1                                                            %
         ichois = -1  ! not valid input, forces users to supply it             %
         ichoisrt(:) = -1                                                      %
         igrill = -3                                                           %
         modcd  = 0                                                            %
         idmpsw = -1  ! not valid input, forces users to supply it             %
         idmpswrt(:) = -1                                                      %
         nminor = 0                                                            %
         kalfa  = 0                                                            %
         bmaxrt = 30                                                           %
         beam_spec(1) = 1                                                      %
         beam_spec(2) = 1                                                      %
         beam_spec(3) = 1                                                      %
         beam_spec(4) = 1                                                      %
         curray_fi = 0.0 
      do l=1,kcrm                                                              %
        heightrt(l)  =  20.0                                                   %
        nthinrt(l)   =   1                                                     %
        islofart(l)  =  -1                                                     %
        thgrilrt(l)   =  0.0                                                   %
        maxrefrt(l) =   2                                                      %
        do k=1,krt                                                             %
          powersrt(k,l) = 0.0                                                  %
          nnkparrt(k,l) = 0                                                    %
          nnkpolrt(k,l) = 0                                                    %
          anzinfrt(k,l) = 0.0                                                  %
          anzsuprt(k,l) = 0.0                                                  %
          anpinfrt(k,l) = 0.0                                                  %
          anpsuprt(k,l) = 0.0                                                  %
        end do                                                                 %
      enddo                                                                    %
c    end curray defaults                                                      
     
     
     
c    genray defaults:
      save_genray_io  =  0
      genray_path = ''          !default must be zero length string
      genray_fi = 0.0
      do k=1,krf
         genraydat(k)='genray.dat'
      enddo
c    end genray defaults



c         
      necsmth   = 0             ! 0 means no smoothing of qrfe and currf       %
      nicsmth   = 0             ! no smoothing same function as nfwsmth
      nfwsmth   = 0             ! no smoothing                                 %
      irfplt    = 0                                                            %
      iside     = 1                                                            %
      navg      = 0                                                            %
      gafsep    = 0.02                                                         %
      gafsep    = 1.e-6         !HSJ 3/06/03 set same as gasep in cray331.f
      ylaunch   = 0.0                                                          %
      xkpar     = 0.0                                                          %
      ykperp    = 0.0                                                          %
      ichmod    = 1                                                            %
      betalm    = 0.0                                                          %
      nhigh     = 4                                                            %
      rfrad1ic  = 0.0                                                          %
      rfrad2ic  = 0.0                                                          %
      relrf     = 0.1                                                          %
      relrf_pow = 0.05
      nrfzon    = 1                                                            %
      nprf(1)   = 1                                                            %
      nprf(2)   = 101                                                          %
      nprf(2)   = kj       ! I think this should be kj ... HSJ                 %
      rfzone(1) = -1.0                                                         %
      rfzone(2) = +1.0                                                         %
      do 1750 i=1,5                                                            %
 1750 timrfp(i) = 1.0e6                                                        %
      itrapech  = 0                                                            %
      ifb       = 1                                                            %
      rfcur     = 0.0                                                          %
      lmode     = 0                                                            %
      ifbprof   = 2                                                            %
      alphaf(1) = 0.2                                                          %
      alphaf(2) = 0.3                                                          %
      alphaf(3) = 0.2                                                          %
c                                                                              %
      nzrffw = 1                                                               %
      do k=1,kzrf                                                              %
        zrffw(k)  =   0.0                                                      %
        pzrffw(k) =   1.0                                                      %
        freqfw(k) = 900.0e6                                                    %
        rnpfw(k)  =   2.0                                                      %
        htsfw(k)  =  20.0                                                      %
      end do                                                                   %
      iswchfw   =   1                                                          %
****  nspfw     =   1                                                          %
      lifw(1)   =   1                                                          %
      nihfw(1)  =  20                                                          %
      impath    =   1                                                          %
      nrayptrt  = 600                                                          %
c              

c                                                                              %

      nrfrad   = 101                                                           %
c                                                                              %
c     initialize input parameters for FASTCD model (Y. R. Lin-liu)             %
c                                                                              %
      do l=1,krf                                                               %
        nkfcd(l) = 1                                                           %
        do k=1,kzrf                                                            %
          xntor(k,l) = 0.0                                                     %
          rpant(k,l) = 0.0                                                     %
        end do                                                                 %
        gamloss(l) = 0.0                                                       %
      end do                                                                   %
c                                                                              %
      do l=1,krf                                                               %
        extcurrf(l)      = 0.0    ! no external RF current from inone          %
        extqerf(l)       = 0.0    ! no external RF heating of electrons        %
        extqirf(l)       = 0.0    ! no external RF heating of ions             %
        extcurrf_id(l)   = 'NONE' ! external RF current                        %
        extqerf_id(l)    = 'NONE' ! external RF electron heating               %
        extqirf_id(l)    = 'NONE' ! external RF ion      heating               %
        extcurrf_amps(l) = 0.0                                                 %
        extqerf_watts(l) = 0.0                                                 %
        extqirf_watts(l) = 0.0                                                 %
        rf_ext_curtot(l) = 0.0                                                 %
        rf_ext_qetot(l)  = 0.0                                                 %
        rf_ext_qitot(l)  = 0.0                                                 %
        cconst = 0.0                                                           %
        call multpl1 (extcurrf_curr(1,l),kj,cconst)  ! zero for starters       %
        call multpl1 (extqerf_qe(1,l)   ,kj,cconst)  ! zero for starters       %
        call multpl1 (extqirf_qi(1,l)   ,kj,cconst)  ! zero for starters       %
      end do                                                                   %
*                                                                              @
c ----------------------------------------------------------------------       @
c --- PELLET FUELING PARAMETERS                                                @
c ----------------------------------------------------------------------       @
c                                                                              @
c  PARAMETERS USED BY OAK RIDGE PELLET MODEL (HOULBERG ET AL.)                 @
c                                                                              @
c  Note that this model is discontinuous with respect to plasma                @
c  transport time: the plasma profiles are changed between time steps.         @
c  The pellet is injected perpendicularly (along a major radius).              @
c                                                                              @
c% nampel    Name of pellet species                                            @
c            'h', 'd', 't', 'dt'                                               @
c            Must be primary species #1 or #2.                                 @
c% pelrad    Initial pellet radius (cm)                                        @
c% vpel      Pellet velocity (cm/sec)                                          @
c% timpel(i) Pellet fueling times (sec, up to 10 by increasing time)           @
c% nbgpel    Number of fast beam ion energy groups (suggest 18)                @
c% pelmod    Model to be used to predict  pellet ablation rate                 @
c            default = 3, Macaulay hydrogenic model                            @
c                                                                              @
c ----------------------------------------------------------------------       @
c                                                                              @
      nampel = 'none'                                                          %
      pelrad =   0.07                                                          %
      vpel   = 800.0e2                                                         %
      nbgpel =  18                                                             %
      do 1760 i=1,10                                                           %
 1760 timpel(i) = HUGE(1.e0)                                                    %
      pelmod = 3                                                               %
*                                                                              @
c ----------------------------------------------------------------------       @
c --- FUSION FLAGS                                                             @
c ----------------------------------------------------------------------       @
c% ifus      1, calculate sources and sinks associated with fusion             @
c            0, do not calculate sources and sinks associated with fusion      @
c               NOTE THAT ifus=1 IS DEFAULT so use this switch to turn         @
c               fusion calcs off.                                              @
c            -1 allows input of 'dt' beam and 'd','t','he' thermal species.    @
c               NOTE: USING HE AS THE THIRD ION IS NOT YET POSSIBLE.           @
c                               SELECT ONLY D AND T AS PRIMARY IONS FOR NOW:   @
c                  ifus =-1 requires that you input the following              @
c                  1) primary ions are individual d and t fluids               @
c                     he is an optional third primary ion fluid                @
c                  2) the two neutral species are d and t,no exceptions        @
c                     (it is not possible to model neutral he at present,      @
c                      probably never will be)                                 @
c                  3) a single impurity species is present (could be he        @
c                     if he is not selected as a primary species)              @
c                  4) electron density and zeff profiles (time dependent)      @
c                     MUST be given with  inenez set to 1                      @
c                  5) if he is not a primary ion you must evolve either        @
c                     the d or t density in simulation mode.The primary        @
c                     species which is run in analysis mode is then obtained   @
c                     by charge balance                                        @
c                  6) if he is a primary ion you must evolve he and one        @
c                     of the two hydrogenic species d or t.                    @
c                     The remaining one will again be abtained                 @
c                     by charge balance.                                       @
c                  7) the code will automatically adjust zeff if necessary     @
c                     to get charge balance (if possible). It is however not   @
c                     reasonable to account for the implied dynamics that  is  @
c                     associated with such a change in zeff simultaneously.    @
c                     Instead a file "zeff.mod" is created which details the   @
c                     changes that were made. You can incorporate these        @
c                     changes as appropriate in the inone file                 @
c                     and do the run again.                                    @
c                     (non-sequential time values can occur in this file due   @
c                     to predictor/corrector iterations. Simply use the        @
c                     sequential values)                                       @
c                     You are responsible for managing this file between runs. @
c                     The code will destroy any existing "zeff.mod" file and   @
c                     create a new one for the current run only if necessary.  @
c                     Hence an old "zeff.mod" file might erroneously be        @
c                     interpreted as being associated with the current run     @
c                     Note that the file begins with a date/time stamp.        @
c                  8) the beam species must be (single fluid) dt               @
c                     with fraction of d given by fdbeam                       @
c                     (if the beam is not on then ifus = -1 reproduces one     @
c                     of the other modes of operation)                         @
c                  9) at the initital time you must give the initial condition @
c                     (i.e., initital density) of he if it is a primary ion    @
c                 10) you must  either give the initital condition for         @
c                     the hydrogenic species which will be run in simulation   @
c                     or you can specify that the ratio of the two hydrogenic  @
c                     species is zfrac (AT THE INITIAL TIME ONLY).             @
c                     If zfrac is set between 0.0 and 1.0 then it is assumed   @
c                     that you want to use zfrac.(set inenez=+1 NOT -1)        @
c                     The boundary conditions  for the hydrogenic              @
c                     densities at r/a=1.0 will be determined by zfrac if      @
c                     selected or by the value of the initial profile at       @
c                     r/a=1.0 if zfrac is not used. The boundary condition     @
c                     for he is given by the initial he density profile at     @
c                     r/a=1. These boundary conditions will be assumed time    @
c                     independent unless you give non zero boundary values in  @
c                     enpb(n,k) where n corresponds to bctime(n) and k is the  @
c                     index for the hydrogenic and he densities to be          @
c                     run in simulation mode.                                  @
c                     (at r/a=0.0 zero gradient in rho is enforced as usual)   @
c                     taupin controls the edge neutral flux for the            @
c                     hydrogenic species run in analysis mode.                 @
c                                                                              @
c% exptl_neutron_rate(j)    j=1,2..nbctim,measured neutron rate, #/sec,        @
c             as a function of bctime(j). This rate will be transferred to     @
c             the output routines for comparison purposes (only).              @
c                                                                              @
c% beam_beam_fusion (integer)                                                  @
c               if ifus ne 0, then setting the integer beam_beam_fusion = 1    @
c               will calculate the beam-beam interactions. These calculations  @
c               take some time (about 1 min per beamline, each time NFREYA     @
c               is called).                                                    @
c               Typically the beam-beam contribution will be less than a 10    @
c               percent of the thermonuclear contribution. For this reason     @
c               you will most likely want an indication of what the beam-beam  @
c               rate is only at the initital and final time points. set:       @
c                                                                              @
c                   beam_beam_fusion = i , where i=1,2,...                     @
c                                         to  calculate the beam beam          @
c                                         reaction every i'th time that the    @
c                                         beam package is called. In between   @
c                                         these time the approximation         @
c                                         described below is used. If you do   @
c                                         not want to use the approximation    @
c                                         simply set beam_beam_fusion=1 which  @
c                                         will force a call to the package     @
c                                         each time the beam is called.        @
c                                         The exact calculations will always   @
c                                         be done at t=time0.                  @
c                                         PROVIDED that the run makes it       @
c                                         through  to t=timmax the exact calc. @
c                                         will also always be done at that     @
c                                         time.                                @
c                                         If the run terminates prematurely    @
c                                         the last beam-beam (and also         @
c                                         beam-thermal,see below) may be an    @
c                                         exact calculation or an approximate  @
c                                         result,there is no way of knowing    @
c                                         a priori.                            @
c                                                                              @
c               The beam_beam rates are                                        @
c               APPROXIMATED (at times that they are not actually computed)    @
c               by assuming that the dependence of the beam integrated         @
c               reaction rate has not changed (due to te evolvement for        @
c               example- which affects the critical fast ion speed)            @
c               from the last time the calculations were done. (Scaling with   @
c               the appropriate fast ion densities will be done however)       @
c                                                                              @
c                   beam_beam_fusion =  0 will skip the calculation altogether @
c               (the results are stored in beam_beamddn, beam_beamdtn          @
c                                                    and beam_beam_tt2n)       @
c% beam_thermal_fusion (integer)                                               @
c               same as beam_beam_fusion except for beam-thermal interactions  @
c                 iddfusb_bulk = 1 MUST BE SET (see description of             @
c                                               iddfusb_bulk elsewhere)        @
c% iaslow    1, account for alpha particle slowing down                        @
c            0, do not account    for alpha particle slowing down              @
c% wtifus    if wtifus > 0.0 then wtifus fraction of the alpha energy          @
c            will go to the ions                                               @
c% fusionvb  integer, if set to 1 gives some detailed output to screen         @
c            (primarly for diagnostic use)                                     @
c                                                                              @
c% iddfusrate 1  use new (thermal) d(d,n)he3 rate coefficient                  @
c                (highly recommended but not the default since                 @
c                 we wish to remain comensurate with other input)              @
c                (see Bosch & Hale, Nuc. Fus. 32, No.4 (1992) 611              @
c                iddfusrate applies only to fusion rate calculations           @
c                done as a result of mhd mixing. The thermal rate and          @
c                beam thermal rates are done with the Bosch+Hale cross         @
c                 sections, independently of this switch.                      @
c                SOME (SEE DDNRATE,DDpRATE FOR EXAMPLE) RATES HAVE LIMITED     @
C                TEMPERATURE RANGES WHICH COULD BE EXCEEDED(EQ 100KEV).        @
c                RATHER THAN STOP                                              @
C                THE CODE WE NOW RETURN VALUES OF THE RATES AT THE LARGEST     @
C                TEMPERATURE ALLOWED BY THE MODEL OF THE DATA.

c                                                                              @
c             0  (default) use old NRL Formulary rate coefficient              @
c                       A graphical comparison of the                          @
c                       two rate coefficients is available in                  @
c                       /u/stjohn/ONETWOPAGES/neutron_rate/                    @
c                                                thermal_rate/ddnhe.?          @
c                       where ? is tex and/or  ps                              @
c                                                                              @
c             BEAM-THERMAL ION FUSION IS MODELED BY ASSUMING THAT THE FAST ION @
c             DISTRIBUTION FUNCTION CAN BE REPRESENTED BY THE ANALYTICAL       @
c             FORM, GIVEN FOR EXAMPLE IN GAFFEY, J., JOUR. PLASMA PHYS. (1976) @
c             VOL 16, PART 2, PG. 149                                          @
c                                                                              @
c             iddfusb      =   0  default, use the old (as presented in        @
c                                 S. SCOTT's thesis),                          @
c                                 approximation to the neutron rate for the    @
c                                 d(df,n)he3 reaction (df = fast deutron)      @
c                                 Thermal motion of ions is neglected          @
c                                 in this model, the effect of bulk rotation   @
c                                 of thermal ions is approximated by shifting  @
c                                 the injected beam energy.                    @
c                                                                              @
c                          =   1  user selected, do a more detailed calc       @
c                                 of the d(df,n)he3 rate as selected by the    @
c                                 switches described below.                    @
c                                 (the Bosch & Hale                            @
c                                 cross sections are used throughout)          @
c              ibslow    NOTE     if ibslow=0 is selected( see beam input      @
c                                 section) the code forces iddfusb=0           @
c                                                                              @
c        the following are required if iddfusb = 1:                            @
c                                                                              @
c  +-------  iddfusb_s = 0,                                                    @
c  |                       use model which neglects bulk and thermal           @
c  |                       motion of ions.                                     @
c  |                       This model is similar to the S. SCOTT model         @
c  |                       (called when iddfusb = 0, see above                 @
c  |                       but uses Bosch & Hale cross                         @
c  |                       sections and does the integral over                 @
c  |                       the energy of the fast ion slowing down             @
c  |                       distribution numerically (i.e., essentially         @
c  |                       exactly,without the approximations involved         @
c  |                       in order to get an analytically integrable          @
c  |                       result). The result obtained with this              @
c  |                       method is correct only if the thermal ions          @
c  |                       have negligible motion and hence the                @
c  |                       actual neutron rate will be                         @
c  |                       underestimated at all temperatures. This model      @
c  |                       does not attempt to account for bulk rotation of    @
c  |                       thermal ions since such a result is misleading      @
c  |                       unless the thermal motion is also included. The     @
c  |                       intended primary use of this option is to           @
c  |                       compare ONETWO results with other codes.            @
c  |                                                                           @
c  |             ibcx      ibcx is an input in the beam section used           @
c  |                       to control inclusion/exclusion of fast ion          @
c  |                       charge exchange. This switch affects the neutron    @
c  |                       rate.                                               @
c  |                       Note that if ibcx=0 is selected the code will       @
c  |                       force iddfusb_s=0 for consistency                   @
c  |                                                                           @
c  |--------- iddfusb_s    = 1 use more refined model as determined by         @
c  |                           following set of switches:                      @
c  |                                                                           @
c  |          iddfusb_bulk = 0 neglect bulk motion (NEUTRON RATE CALCS ONLY)   @
c  |                       = 1 include bulk motion (default)                   @
c  |                           Note that this only refers to the model         @
c  |                           obtained when iddfusb_s=1.                      @
c  |                           bulk motion is always negelected in the model   @
c  |                           used when iddfusb_s=0                           @
c  |                            THIS SWITCH ONLY CONTROLS THE BEAM ENERGY      @
c  |                            USED IN THE NEUTRON RATE CALCS. THE            @
c  |                            MODIFICATION OF THE FAST ION DISTRIBUTION      @
c  |                            FUNCTION {DUE TO AN ALTERED DEPOSITION         @
c  |                            PROFILE OF FAST IONS CAUSED BY CHANGES IN      @
c  |                            CHARGE EXCHANGE CROSS SECTIONS} IS  NOT        @
c  |                            AFFECTED BY THIS SWITCH. TO GET TOTAL          @
c  |                            SELF CONSISTENCEY WITH NO BULK MOTION YOU MUST @
c  |                            RUN WITH NO TOROIDAL ROTATION !!!!!!!!!        @
c  |                                                                           @
c  |   icalc_cxfactor = 0     Use this switch to control charge exchange       @
c  |                          treatment. icalc_cxfactor=0 (default ) means     @
c  |                          that the charge exchange probability is taken    @
c  |                          as independent of the fast ion speed.            @
c  |                          This assumption is made in the ONETWO fast ion   @
c  |                          model and thus should also be used when the      @
c  |                          neutron rate is calculated.                      @
c  |                                                                           @
c  |                  = 1     Use icalc_cxfactor=1 if you want to do the       @
c  |                          charge exchange probability integral instead     @
c  |                          of assuming it to be constant.  This method      @
c  |                          is of course slower than the above method        @
c  |                          and is not fully consistent with other fast ion  @
c  |                          "stuff" in the code . It will give you an        @
c  |                          indication how important the assumption of       @
c  |                          constant cx probability is. See the ONETWOPAGES  @
c  |                          docs for an exact description of this switch.    @
c  |                                                                           @
c  |  THE FOLLOWING IS CURRENTLY NOT AVAILABLE                                 @
c  |                                                                           @
c  |             ddfusb_t = le 1.0, (default = 1.0)                            @
c  |                       neglect neutrons produced by part of fast           @
c  |                       ion distribution function that is above the         @
c  |                       energy ddfusb_t*e0,where e0 is the                  @
c  |                       injected energy normally ddfusb_t =1.0              @
c  |                       should be used. Values less than 1 are not          @
c  |                       useful generally.                                   @
c  |                       (note : the model assumes the fast ion high-energy  @
c  |                       tail is due to collisions with electrons only.      @
c  |                       The collision of fast ions amongst themselves,      @
c  |                       indpendent of the three injection energies,         @
c  |                       is neglected (i.e., the density of fast ions        @
c  |                       is supposed to be small for the analytic            @
c  |                       distributions used here to make sense).             @
c  |             ddfusb_t > 1.0 (user-specified)                               @
c  |                       account for neutrons produced by                    @
c  |                       part of distribution function between               @
c  |                       the injected energy e0 and ddfusb_t*e0.             @
c  |                       (normally,due to the rapid decay of the             @
c  |                       fast ion distribution function above the injected   @
c  |                       energy, ddfusb_t=1.25 should suffice)               @
c  |                                                                           @
c             some info and graphs are available in ONETWOPAGES                @
c                                                                              @
c ----------------------------------------------------------------------       @
c                                                                              @
      ifus                = 1                                                  %
      iaslow              = 1                                                  %
      wtifus              = 0.0                                                %
      iddfusrate          = 0    ! governs thermal d(d,n)he3 only              %
      iddfusb             = 0                                                  %
      iddfusb_s           = 0                                                  %
      ddfusb_t            = 1.0  ! neglect fast ions above inject energy       %
      iddfusb_bulk        = 1    ! include bulk motion of ions                 %
      icalc_cxfactor      = 0    ! use simple cx model                         %
      beam_beam_fusion    = 0    ! no beam_beam fusion                         %
      beam_thermal_fusion = 1    ! include beam_thermal fusion                 %
      fusionvb            = 0                                                  %
      do j=1,kbctim                                                            %
        exptl_neutron_rate(j) = 0.0                                            %
      end do                                                                   %
*                                                                              @
c **********************************************************************       @
c *** THIRD NAMELIST (NAMELIS3) ****************************************       @
c **********************************************************************       @
c                                                                              @
c ----------------------------------------------------------------------       @
c                          RUN TYPE FLAGS                                      @
c ----------------------------------------------------------------------       @
c% ltest_code    set to 1 to do some mhd/transport coupled testing.            @
c                the user not involved in development of the code              @
c                should leave ltest_code = 0 (its default value).              @
c% use_Bp0       differentiation of the pressure can be done with              @
c                finite differences by using the pressure defined over the     @
c                psi grid. This is the default and is used if use_Bp0 =0       @
c                Alternatively,if use_Bp0 >=1, then chain rule :               @
c                   d Press /d psi = d Press /d rho * d rho/d psi              @
c                is used. In this case d Press / d rho is still determined     @
c                using finite differences but                                  @
c                       dhro/d psi = 1/(R0 Bp0)   is used                      @
c                To handle the situation near the magnetic axis                @
c                three  options are available                                  @
c                   1) use_bp0 =1:                                             @
c                      evaluate B0 one half grid space away from magnetic      @
c                      axis and use that value at the magnetic axis to         @
c                      determine the derivative                                @
c                   2) use_bp0 = 2:                                            @
c                      a quadratic  fit to d Press/d psi for grid points 2,3,4 @
c                      is done. This polynomial is extrapolated to             @
c                      psi  = axis   to  get d press /d psi at the mag axis    @ 
c                   3) use_Bp0 > 2                                             @
c                      set d press /dpsi at axis  = dpress / dpsi one grid     @
c                      point away                                              @
c% curtype       set to 0 (default) to use <J*R0/R> in determination           @
c                of ohmic current. ( this is the old onetwo way)               @
c                set curtype to 1 to use < J dot Btotal/Bt0>.                  @
c                setting curtype=1 also modifies the calculation               @
c                of the total current contributions by accounting for          @
c                the diamagnetic effect. But there will be some                @
c                error in the ohmic heating which is given by                  @
c                Qohm = Etor * <J*R0/R> 
c% ifixshap    1: Fix flux surfaces in time.  Perform only a single,           @
c                 initial equilibrium calculation or read the equilibrium      @
c                 from an eqdsk,see irguess below.                             @
c              0: Evolve flux surfaces dynamically.  Perform several           @
c                 equilibrium calculations.                                    @
c                 use ifixshap = 0 if tdem mode is selected, see below.        @
c                      ( Tdem mode is the time dependent eqdisk mode)          @
c% mhdmode     'no coils': The flux (in volt-sec/rad) is specified             @
c                  on the boundary of the MHD grid (using array flxeqbcd).     @
c                  the coils are not explicitly modeled in this case.          @
c                  You must simulate effects of coils outside the              @
c                  MHD grid indirectly by using appropriate boundary           @
c                  conditions. Toroidal currents (other than that due to       @
c                  the plasma) are not allowed inside the mhdgrid.             @
c                  Consequently it is usually necessary to modify              @
c                  the mhdgrid. (You MUST set xdim,ydim and redge!)            @
c                  Use 'no coils' option if tdem mode is selected (see below)  @
c              'coils': Coil option -- model field-shaping coils,using         @
c                  psi loop and/or magnetic probe values,see below.            @
c                  this option requires a full coil set (through the           @
c                  associated Green's table) and is implemented only           @
c                  for  DIII-D type equilibria. non DIII-D runs                @
c                  must use the 'no coils' option.                             @
c% mhdmethd         set mhdmethd = 'green' to get solution of Grad-Shafranov   @
c                   equation using Green's function for entire                 @
c                   MHD grid (this is the slowww method.It should not          @
c                   be used routinely. It was included here primarily          @
c                   as an aid to verify solutions obtained with the            @
c                   new cyclic reduction solver used in subroutine CMPLTCYR.)  @
c                   This method is much faster and should be used routinely.   @
c                      mhdmode  = 'no coils' for fixed boundary cases          @
c                      mhdmethd = 'cycred' yields the cyclic reduction         @
c                                  solution of the fixed boundary problem.     @
c                      mhdmethd = 'sorpicrd' yields the successive             @
c                                  overrelaxation method.                      @
c                      mhdmethd = 'toq' to get inverse equilibrium solver TOQ  @
c                                 TOQ will be spawned from Onetwo. The         @
c                                 necessary input file for TOQ (intoq) will    @
c                                 be generated by Onetwo from information     @
c                                 supplied below(see TOQ input)                @
c                      mhdmethd ='tdem' for time-dependent eqdsk mode          @
c                if mhdmethd ='tdem' (i.e., the multiple eqdsk case)           @
c                then this list of times is not used                           @
c                (the netCDF file, which has a name that eqdskin points to)    @
c                  supplies the data instead)                                  @
c                The tdem mode also requires that you set the following        @
c                switches (these are set automatically for you even if         @
c                you set them differently or didn't set them at all)           @
c                if mhdmethd = 'tdem' then the following settings are FORCED   @
c                   ieqdsk   =  0                                              @
c                   mhdmode  = 'no coils'                                      @
c                   ifixshap =  0                                              @
c                  also eqdskin must be set to an appropriate netCDF           @
c                  file  which we CANNOT set for the user                      @
c% ishot          ishot should be set to shot number in netCDF file            @
c                The tdem mode NEVER solves the MHD (i.e., Grad-Shafranov)     @
c                equation. Instead it reads a precomputed file of MHD          @
c                data (in netCDF format) as a function of SPACE AND TIME       @
c                and interpolates this data for the values required at each    @
c                time step.                                                    @
c                                                                              @
c              if mhdmode = 'no coils' and mhdmethd = 'sorpicrd' then          @
c% optwe        optwe = .true. is used to tell the code that optimal           @
c              external relaxation parameter is to be found. if                @
c              if optwe = .false. then it is assumed that an optimal           @
c% optomegae    value is input in optomegae. failing this the code will        @
c              use a crude estimate to get optomege.                           @
c% optwi        optwi = .true. is used to tell the code that optimal           @
c              internal relaxation parameter is to be found. if                @
c              optwi = .false. then it is assumed that an optimal              @
c              optomegai value is input in optomegai. failing this the code    @
c              will use a crude estimate to get optomegai                      @
c                                                                              @
c --- THE FOLLOWING IS IMPLEMENTED FOR SYMMETRIC SOR ONLY:                     @
c                                                                              @
c mresidual                                                                    @
c nfitpoints                                                                   @
c splinefitbd                                                                  @
c fitboundary                                                                  @
c                                                                              @
c  NOTE             irguess specifies the source of the initial                @
c                   equilibrium. Thereafter ifixshap determines if             @
c                   further equilibria are to be calculated. Both              @
c                   the initial and subsequent equilibria depend               @
c                   on how the boundary condition for the poloidal             @
c                   flux in the Grad-Shafranov equation is specified           @
c                   using the switch mhdmode.                                  @
c                   if irguess .lt. 0 and ifixshap = 1 then the Grad-Shafranov @
c                   equation will never have to be solved and boundary         @
c                   conditions are irrelevant.                                 @
c                   however if ifixshap = 0 AND/OR irguess .ge. 0              @
c                   then the Grad-Shafranov equation will have to be solved    @
c                   at least once and consequently you must supply             @
c                   all information required to arrive at a solution           @
c                   (in particular since this is a free boundary code          @
c                   the limiter points as well as the boundary conditions      @
c                   determine the solution)                                    @
c                                                                              @
c% fixfcoil          INTEGER input number, pertains to treatment of            @
c                    fcoil currents in the mhdmode = 'coils' option.           @
c                    not used if mhdmode = 'no coils'.                         @
c                    (also not used if ifixshap = 1 and irguess .ne. 0)        @
c                    set fixfcoil = 0 if the fcoil currents                    @
c                    are to be calculated as part of the MHD solution          @
c                    procedure. This is done by simultaneously adjusting       @
c                    all fcoil currents in such a way that the value of        @
c                    psi on the psi loops, magnetic probes                     @
c                    and experimental fcoil currents                           @
c                    values are reproduced in a least squares sense.           @
c                    the input                                                 @
c                    fcoil currents will be compared to the calculated         @
c                    fcoil currents in the output file and plot and            @
c                    a figure of merit describing the discrepancy              @
c                    between the calculated and input values will be           @
c                    generated if the f coil current values are input          @
c                    (in the third namelist of inone, as fcoilcur).            @
c                                                                              @
c                    set fixfcoil = 1 if the f coil currents are to be         @
c                    taken from the input file and therefore NOT               @
c                    CALCULATED during the MHD solution process.               @
c                    in this case the calculated psi loop and probe            @
c                    values will generally not match the input values          @
c                    The calculated and input psi loop and probe values        @
c                    will be compared in the output file and plots.            @
c                    A figure of merit describing the discrepancy              @
c                    between the calculated and input values                   @
c                    will be generated for this case if the input              @
c                    values are present in the third namelist of inone.        @
c                                                                              @
c                    No effort other than the adjustment of fcoil              @
c                    currents described above is made to minimize              @
c                    the discrepcancy between the calculated and               @
c                    measured f coil currents or between the calculated        @
c                    and measured psi loop and probe values. We do not         @
c                    have a feedback mechanism in place to adjust              @
c                    resistivity or other parameters to force                  @
c                    agreement.                                                @
c                    YOU SHOULD NORMALLY USE FIXFCOIL = 0                      @
c                                                                              @
c                        if fixfcoil = 0 then the user may select the          @
c                        experimental data used to determine the fcoil         @
c                        currents as follows:                                  @
c% ifitpsi           =1  means use experimental psi loop values                @
c                        (obtained by interpolation in time from the input     @
c                         values stored in flxeqbcd) in determination          @
c                         of the fcoil currents.                               @
c  ifitpsi           =0  means ignore psi loop values in determination         @
c                        of fcoil currents.                                    @
c% ifitprob          =1  means use experimental magnetic probe values          @
c                        (obtained by interpolation in time from the input     @
c                         values stored in expmp2) in determination            @
c                         of the fcoil currents.                               @
c  ifitprob          =0  means ignore magnetic probe values in determination   @
c                        of fcoil currents.                                    @
c% fitfcur           =1  means use experimental fcoil currents                 @
c                        (obtained by interpolation in time from the input     @
c                         values stored in fcoilcur) in determination          @
c                         of the fcoil currents.                               @
c  fitfcur          =0  means ignore experimental fcoil currents               @
c                                                                              @
c                     NOTE: the difference between fixfcoil and fitfcur        @
c                       is as follows. Fitfcur only has meaning if             @
c                       fixfcoil=0. For only then do we calculate fcoil        @
c                       currents. There are three sources of information       @
c                       available for determination of f coil currents:        @
c                            a) the psi loops                                  @
c                            b) the magnetic probes                            @
c                            c) the experimentally measured f coil currents    @
c                            d) Rogowski coils (not implemented yet)           @
c                       It is generally not a good idea to set the calculated  @
c                       fcoil currents equal to the measured ones (by setting  @
c                       fixfcoil = 0). Rather the best results are abtained    @
c                       by determining the fcoil currents through a least      @
c                       squares determination. fitfcur is the switch that      @
c                       determines wheter or not the experimental f coil       @
c                       currents will be included as part of the data in       @
c                       the least squares fit.                                 @
c                      YOU SHOULD NORMALLY USE FITFCUR = 0 (the experimental   @
c                      fcoil currents apparently are responding to eddy        @
c                      effects not modeled in the MHD calculations so          @
c                      it is best not to use them. Quite often using           @
c                      the experimental f coil currents causses the            @
c                      MHD equilibrium not to converge).                       @
c                                                                              @
c                          at least one of ifitpsi,ifitprob must be =1         @
c                          otherwise there is no information with which        @
c                          the fcoil currents can be determined.               @
c                          feedback from plasma resistivity or (are we         @
c                          crazy enough?) transport parameters is miles        @
c                          or years down the road.                             @
c% iecurr            =1  means include ecoil  currents in MHD solution.        @
c  iecurr            =0  means exclude ecoil  currents in MHD solution.        @
c% ivessel           =1  means include vessel currents in MHD solution.        @
c  ivessel           =0  means exclude vessel currents in MHD solution.        @
c                  if the ecoil and or vessel currents are included            @
c                  in the MHD model they must be given explicitly as           @
c                  input functions of time (in arrays ecurrt and               @
c                  vesscur respectively). No models exist in ONETWO            @
c                  to allow calculation of these quantities at present.        @
c                  furthermore the vessel and ecoil currents are treated       @
c                  as exact (i.e., not measured) values.                       @
c  if the fcoil currents are fit (i.e., fixfcoil = 0) then the following       @
c  least squares weight factors must also be supplied. These are               @
c  necessary so that psiloop,probe and fcoil current data can be               @
c  handeled in a way that gives reasonable weight to each. Without             @
c  weighting the experimental fcoil currents would force the solution          @
c  to ignore the psiloop and probe values.                                     @
c                                                                              @
c% errpsilp            weight for psi loops                                    @
c% errmprge                       magnetic probes                              @
c% errfcoil                       experimental f coil currents                 @
c% errvescr                 experimental vessel currents                       @
c% errecoil                              ecoil                                 @
c                      each of the above five weights is to be input           @
c                      as a fraction of the measured value. for example        @
c                      if errpsilp = 0.01 then it will be assumed that all     @
c                      psiloop values have a standard deviation of 1% of       @
c                      their input values.                                     @
c% minchisq            set minchisq = 1 to roughly minimize                    @
c                      chisq each time the MHD package is called.              @
c                      set minchisq = 0 to turn off this calculation.          @
c                      NOTE: the only degree of freedom we have to             @
c                            minimize chisq is the total toroidal current,     @
c                            (which enters the transport calcs only as a       @
c                            boundary condition). Hence using minchisq = 1     @
c                            has the effect of changing the boundary           @
c                            condition (on Faraday's law). Using this          @
c                            option is expensive,but may lead to better        @
c                            statistical agreement with the measured           @
c                            magnetics.  However if chisq is too               @
c                            large then just changing the current              @
c                            normalization may not reduce chisq enough.        @
c                            The only recourse we have then is to modify       @
c                            the transport model (particularly Faraday's law)  @
c                                                                              @



c  START TOQ input ---http://fusion.gat.com/toq/-----------                    @
c  Toq has some settings that are indepedent of Onetwo. These must be supplied @
c  here (or the defaults given here will be used):                             @
c  
c     COMMON INPUT (REPEATED FROM ABOVE) REQUIRED TO RUN TOQ EQUILIBRIA:       @
c     toq_path =' ' set this to a fully qualified execcutable toq if you do not@
c                     want the default TOQ that Onetwo will find automatically @
c                   (file rfmod.f90 has the paths that Onetwo searches) 
c     ifixshap =0: Evolve flux surfaces dynamically.  Perform several          @
c                  equilibrium calculations.                                   @
c     mhdmode  = 'no coils' for fixed boundary cases                           @
c     mhdmethd = 'toq' to get inverse equilibrium solver TOQ                   @
c     mhdonly        =0 no effect                                              @
c                    =1 run only a single equilibrium calculation, do not      @
c                       do any transport!                                      @
c     isym         updown symmetric=1, updown asymmetric=0                     @
c                  see updownsym below                                         @
c     eqdsk        name of gfile to be read in (fneqdsk in TOQ)                @
c     comp_methd_eqdsk GOVERNS HOW FFP,PPRIME,F,P are iterpolated onto         @
c                 an EFIT type eqdsk. see description below.                    @
c                 NOTE that you can also rewrite any EFIT type eqdsk           @
c                 produced by Onetwo by using the program ploteq.              @
c                 Ploteq plots up to 5 eqdsk together for easy comparison.     @
c                 One option in ploteq is to rewrite the eqdsk with the        @
c                 on axis values smoothly extrapolated.                        @
c                 /u/stjohn/ploteq/ploteq_129_129                              @
c                 Ploteq also writes a file qplotdata.eqdsk that has in it     @
c                 tables of vnormal,rhonormal,psi,curden,etc.                  @
c     use_bp0      see above                                                   @
c                 
c     do_eqplot    see above                                                   @
c
c     TOQ SPECIFIC VALUES (USED ONLY IF TOQ IS SELECTED WITH ABOVE SWITCHES)   @
c     ieqdtoq=1          1-> read in gfile, 0->don't read in gfile
c                        if reading in gfile, namelist "input" 
c                        must appear in this file
c                        WARNING: ieqdtoq=1 will overwrite dskeqdata!

c     npsi_toq=19   nthet_toq=17      mesh size-- psi & theta points
c     npsi_toq=35   nthet_toq=33       acceptable values are listed here
c     npsi_toq=67   nthet_toq=65
c     npsi_toq=131  nthet_toq=129
c     npsi_toq=259  nthet_toq=257
c     THE FOLLOWING TOQ INPUTS ARE SUPPLIED BY ONETWO :
c     rmax            minor radius (  Not an input value will be supplied      @
c                                     by Onetwo)
c     rzero           major radius input this as rmajor (see above)
c     baxis0          vacuum b field in gauss
c     fixcur=.true.   will adjust magnitude of ffprime to fix total current
c     totcur          totcur is set above
c                     if fixcur=.true.


c     imislo              sometimes data near axis is garbage.
c                         in that case set imislo=3 to use 
c                         extrapolation to axis which should 
c                         clean it up,DEFAULT 0
c     alpsi               determines the relationship between psic--psi 
c                         coordinate and psiv--psi value=poloidal flux(*2pi)
c                         alpsi=0 makes psic = sqrt(psiv) normalized to 0-1
c                         for alpsi<0
c                         psi_norm(i)=(sin(psic(i)*alpsi*pi*0.5))**2/
c                                      sin(-alpsi*pi*0.5)**2
c                         idea for alpsi~-1 is to put lots of surfaces 
c                         near axis and near edge, DEFAULT 0.0
C     modelbnd              modelbnd=2 get boundary shape from dskeqdata file
C                          modelbnd=0 use shape formulas below to get boundary
C     ishape=4             for dee shape with elongation=eshape and 
c                         triangularity=sin(xshape)
c                         set ishape=0 to get a circle
c     eshape=2.0000          these two parameters not used for ishape=0
c     xshape=0.4             eshape=2.0 xshape=0.5236 to get triangularity=0.5
c     equiltype='ffprime'    ffprime: use ff' and p' to determine equilibrium
c                            other equiltype values are 'qsolver' and 'jdotb'
c  modelf=2               polynomial model for ffprime--see fsetup.f
c  nffp=2                 ffprime polynomial degree
c  ffpcoef= 1.0 -1.0      ffprime polynomial coefficients
c  modelp=2               polynomial model for pprime--see psetup.f
c  modelq=2               q from eqdata file
c  nppcoef=2              pprime polynomial degree
c  ppcoef=1.0 -1.0        pprime polynomial coefficients
c  betafix='no'           will adjust magnitude of pressure to achieve
c                         beta on axis value of betiso.
c                         if betafix='none' pressure not modified This is no correct HSJ
c                         betafix is 2 character variable !!!!!
c   betiso=0.001          desired beta on axis
c  ieqdsk_toq=0               if =2 use dskeqi to get initial guess
c                         =0 don't read dskeqi
c                         dskeqi is dskeqo generated from a previous run
c  prtbal=.true.          print out file dskbalnew for baloo
c  prtboot=.true.         print out file dskboot for bootstrap current analysis
c  prtgato=.true.         print out file dskgato for gato
c
c  prteqdata=.true.       print out file dskeqdatao for iput to another toq run
c                         perhaps after modifying p,p',f,f',q or boundary
c                         to use in a later run, change dskeqdatao to dskeqdata
c                         and set desired model parameter to 2,
c                         modelx=2 means get arrays from dskeqdata file
c                         modelp=2  get pressure profile from dskeqdata
c                         for other modelp options see routine psetup.f
c                         modelf=2  get f and ff'  profile from dskeqdata
c                         for other modelf options see routine sigset.f
c                         modelbnd=2 get boundary shape from dskeqdata file
c                         modelq=2 get q profile from dskeqdata file
c  prtfixb=.true.         print out dskfixb file which can be used with
c                         fixbnd.x to generate an eqdsk
c  prtflux=.true.         write out flux at each iteration in fort.10
c  bavmg=0.3              multigrid solver back-averaging
c  minmg=3                minimum number of multigrid cycles
c  maxmg=20               maximum number of multigrid cycles
c  tolmg=1.e-4            multigrid solver tolerance
c  loopinmg=2     
c  nhior=15               don't ask you don't want to know
c  nbndry=199             don't ask you don't want to know
c  nbug(1)=0              controls some multigrid printout--set to zero to get
c                         less output to stdout
c  nbug(2)=-1
c  nbug(3)=0
c  nbug(4)=0
c  toleq_toq =1.e-5       equilibrium convergence criterion. Value for this is 
c                         taken from Onetwo variable with same name. 
c  iteqmx=50              maximum number of equilibrium iterations. Value for this 
c                         is taken from Onetwo variable iteq.     
c
c   rmax_toq             toq variable rmax ... not defined in toq ??
c   
c   zeff_toq=1.5          single zeff value  passed to toq as zeff
c                        (used only when toq is run from Onetwo
c                         without onetwo specific input)

c             INPUTS FOR SECOND NAMELIST IN INTOQ FILE:
c
c  fneqdsk                value will be taken from eqdskin input name or
c                         generated internally by Onetwo
c  psifactr_toq =0.99999
c  dsrat=0.005
c  updownsym='a'          'u' = use upper boundary, 'd' = use lower boundary
c                         'a' = updown asymmetric
c                          Set internally depending on value of isim
c                          if isym =1 then updownsym ='a' is set internally.
c                          if isym =0 then updownsym must be set to 'u' or 'd'
c                           in inone.
c  npts=129               the number of points on the boundary specification
c                         stored in dskeqdata
c                         this should be >= nthet, but not a requirement
c  qedge_toq=0.
c  igimblett=.false.
c 

c  END TOQ input -------------------------                                     @





c  MOTIONAL STARK EFFECT:                                                      @
c --- motional stark effect vectors are as follows:                            @
c ---  nstark   max number of channels (parameter,defined in mhdpar.i)         @
c ---  Tstark(nstark,mxtbcmhd) tangent of mag. pitch angle                     @
c ---  sigstark(nstark,mxtbcmhd) error in tstark                               @
c ---  fwtstark(nstark,mxtbcmhd) fit weight of each channel                    @
c ---  rstark(nstark,mxtbcmhd) major rad. of measurement                       @
c ---  zstark(nstark,mxtbcmhd) elevation of measurement                        @
c --- coefficient vectors for calculation of tstark.                           @
c --- a1-4stark(nstark,mxtbcmhd)                                               @
c% timestark(i) i = 1,2..mx of mxtbcmhd times at which MSE data are given      @
c% use_stark    logical, if true stark effect measurements                     @
c               will be used, otherwise they will be ignored.                  @
c               use_stark = .false.  is the default                            @
c ----------------------------------------------------------------------       @
c                                                                              @
      optwe       = .false.                                                    %
      optwi       = .false.                                                    %
      optomegai   = 0.0                                                        %
      optomegae   = 0.0                                                        %
      splinefitbd = .false.                                                    %
      fitboundary = .false.                                                    %
      ltest_code  = 0                                                          %
      ifixshap    = 1                                                          %
      curtype     = 0                                                          %
      fixfcoil    = 0                                                          %
      fitfcur     = 0                                                          %
      errpsilp    = 0.03                                                       %
      errmprbe    = 0.05                                                       %
      errfcoil    = 0.10                                                       %
      errvescr    = 0.10                                                       %
      errecoil    = 0.10                                                       %
      ifitpsi     = 1                                                          %
      ifitprob    = 1                                                          %
      iecurr      = 1                                                          %
      ivessel     = 0                                                          %
      use_Bp0     = 0                                                          %
      call zeroa (fwtstark,mxtbcmhd*nstark)                                    %
      call zeroa (rstark,mxtbcmhd*nstark)                                      %
      call zeroa (zstark,mxtbcmhd*nstark)                                      %
      call zeroa (tstark,mxtbcmhd*nstark)                                      %
      call zeroa (timestark,nstark)                                            %
      call zeroa (sigstark,mxtbcmhd*nstark)                                    %
      call zeroa (a1stark,mxtbcmhd*nstark)                                     %
      call zeroa (a2stark,mxtbcmhd*nstark)                                     %
      call zeroa (a3stark,mxtbcmhd*nstark)                                     %
      call zeroa (a4stark,mxtbcmhd*nstark)                                     %
      use_stark = .false.                                                      %
      mhdmode   = 'coils'                                                      %
      mhdmethd  = 'cycred'                                                     %
*                                                                              @
c ----------------------------------------------------------------------       @
c --- GEOMETRIC PARAMETERS                                                     @
c ----------------------------------------------------------------------       @
c% xdim          Width  of box (cm) when MHDMODE = 'no coils'                  @
c% ydim          Height of box (cm) when MHDMODE = 'no coils'                  @
c% redge  Inside radius of box (cm) when MHDMODE = 'no coils'                  @
c                   Note: when using the 'no coils' option, xdim, ydim         @
c                         and redge define the computational MHD grid          @
c                         vectors rmhdgrid(i), zmhdgrid(j), i = 1,..nw,        @
c                         j = 1,..nh. Since there are no coil currents         @
c                         in this option it usually does not make sense        @
c                         to define xdim,ydim,redge in the same way as         @
c                         is done with the coils option. Instead you will      @
c                         probably find that the boundary of the grid          @
c                         must be very close to the plasma,in order that       @
c                         the interpolation between boundary values makes      @
c                         sense and that shape control is possible.            @
c                         don't forget that the limiter still limits the       @
c                         radial and vertical extent of the plasma. In         @
c                         particular setting up an mhdgrid with xdim,ydim,     @
c                         redge means that you will probably also have to      @
c                         define a consistent limiter point set as well.       @
c                         if mhdmode = 'no coils' xdim,ydim,redge MUST be      @
c                         input  (there is no default here).                   @
c                         if mhdmode = 'coils' xdim,ydim,redge will not        @
c                         be  used(the grid will be generated from the         @
c                         eqdsk or the Green's table as appropriate).          @
c                                                                              @
c% nlimiter     Number of points defining limiter (max = maxlimpt)             @
c                   the limiter can be specified in one of three ways:         @
c                   a)in the third namelist of file inone. to use              @
c                     this option set nlimiter equal to the number             @
c                     of (xlimiter,ylimiter) points given in the               @
c                     third namelist of inone (and specify nlimiter            @
c                     values of xlimiter,ylimiter).                            @
c                   b)the default values given below may be used.              @
c                     leave nlimiter UNSPECIFIED and DO NOT include            @
c                     the (xlimter,ylimiter) points  in file inone to get      @
c                     this option.                                             @
c                   c)the limiter points are read from the eqdsk.              @
c                     set nlimiter equal to a negative number                  @
c                     in the third namelist of inone to get this option.       @
c                     you must of course supply an eqdsk for this option.      @
c                     (most eqdsks have limiter points in them. If it          @
c                     does not a warning will be printed and the code          @
c                     will stop). Use this option for 'tdem' mode also         @
c                     (see RESTART PARAMETERS below for eqdsk specification)   @
c   Note: When running with an eqdsk (i.e., irguess .ne. 0) it is recommended  @
c                   that option c be used since this will presumably be the    @
c                   correct limiter position for the equilibrium to be         @
c                   analyzed. If MHD calculations are done in ONETWO           @
c                   the results WILL DEPEND on the limiter position            @
c                   if the plasma touches the limiter at any time during       @
c                   the analysis. It is up to you to make sure that            @
c                   the limiter points used are correct! Consequently          @
c                   option b,above, should be used only if a fixed             @
c                   equilibrium is used,since in that case,the limiter         @
c                   position is irrelevant.                                    @
c                  It is now possible to generate the initial                  @
c                  equilibrium in ONETWO so no eqdsk exists a priori.          @
c                  (this is done using irguess = 0,see below).                 @
c                  in this case option a,or c can be used.                     @
c                  if  option c is used then the eqdsk can be any              @
c                  eqdsk  which has the appropriate limiter points.            @
c                  all other information in the eqdsk will be ignored          @
c                  provided that irguess = 0.                                  @
c                  [IRGUESS = 0 IS CURRENTLY NOT IMPLEMENTED]                  @
c                                                                              @
c% xlimiter(i)  List of x coordinates for limiter (m)                          @
c% ylimiter(i)  List of y coordinates for limiter (m)                          @
c% ifill        1, Fill in additional limiter points                           @
c               0, Do not fill in additional limiter points                    @
c  ifill is not used                                                           @
c                  ifill is ignored in input file, no filling of               @
c                  limiter will be done                                        @
c% greentab       Name of Green's function table                               @
c                      unless the Green's table name was changed from          @
c                      the standard defaults, the code will determine          @
c                      the proper name if greentab is not specified            @
c                      in the input. The default name will change with         @
c                      the equilibrium grid size. It consists of 8             @
c                      characters, starting with a t, followed by number       @
c                      of grid points in r direction, nw, followed by number   @
c                      of grid points in z direction, nh, followed by          @
c                      however many characters of 'd3d' required to            @
c                      complete an 8-character name. example:                  @
c                      't65129d3' would be the default name for                @
c                      for a 65 by 129 Green's table for d3d. No option        @
c                      for machines other than d3d exists.                     @
c               NOTE   YOU CAN ONLY USE A GREEN'S TABLE WHICH                  @
c                      MATCHES IN SIZE (NW,NH,NCOILS,ETC)                      @
c                      THE PARAMETERS SET IN ONETWO. A CHECK FOR THIS          @
c                      IS MADE IN SUBROUTINE READGREN AND ERROR MESSAGE        @
c                      WILL BE ISSUED IF INCONSISTENCIES ARE DETECTED.         @
c                      A GREEN'S TABLE IS REQUIRED WHEN                        @
c                            mhdmode ="coils" and irguess .ge. 0               @
c                            (independent of ifixshap)                         @
c if irguess < 0 and ifixshap = 1, then a Green's table is not required        @
c                      NO GREEN'S TABLE IS REQUIRED IF MHDMODE ="NO COILS".    @
c                      HOWEVER YOU MUST SUPPLY BOUNDARY CONDITIONS             @
c                      FOR PSI IN THIS CASE (SEE MHDMODE)                      @
c ----------------------------------------------------------------------       @
c                                                                              @
      if (machinei .ne. 'doub-iii')  go to 9080                                %
      xdim     = 0.0                                                           %
      ydim     = 0.0                                                           %
      redge    = 0.0                                                           %
      greentab = 'default'                                                     %
      ifill    = 1                                                             %
      nlimiter = 69                                                            %
c                                                                              %
      xlimiter( 1) =  0.961                                                    %
      xlimiter( 2) =  0.961                                                    %
      xlimiter( 3) =  0.961                                                    %
      xlimiter( 4) =  0.961                                                    %
      xlimiter( 5) =  0.961                                                    %
      xlimiter( 6) =  0.961                                                    %
      xlimiter( 7) =  0.961                                                    %
      xlimiter( 8) =  0.961                                                    %
      xlimiter( 9) =  0.961                                                    %
      xlimiter(10) =  0.961                                                    %
      xlimiter(11) =  0.961                                                    %
      xlimiter(12) =  0.961                                                    %
      xlimiter(13) =  0.961                                                    %
      xlimiter(14) =  0.961                                                    %
      xlimiter(15) =  0.961                                                    %
      xlimiter(16) =  1.134                                                    %
      xlimiter(17) =  1.300                                                    %
      xlimiter(18) =  1.400                                                    %
      xlimiter(19) =  1.500                                                    %
      xlimiter(20) =  1.600                                                    %
      xlimiter(21) =  1.680                                                    %
      xlimiter(22) =  1.908                                                    %
      xlimiter(23) =  1.850                                                    %
      xlimiter(24) =  1.855                                                    %
      xlimiter(25) =  1.859                                                    %
      xlimiter(26) =  1.861                                                    %
      xlimiter(27) =  1.861                                                    %
      xlimiter(28) =  1.860                                                    %
      xlimiter(29) =  1.858                                                    %
      xlimiter(30) =  1.854                                                    %
      xlimiter(31) =  1.849                                                    %
      xlimiter(32) =  1.842                                                    %
      xlimiter(33) =  1.833                                                    %
      xlimiter(34) =  1.823                                                    %
      xlimiter(35) =  1.812                                                    %
      xlimiter(36) =  1.799                                                    %
      xlimiter(37) =  1.784                                                    %
      xlimiter(38) =  1.767                                                    %
      xlimiter(39) =  1.749                                                    %
      xlimiter(40) =  1.729                                                    %
      xlimiter(41) =  1.707                                                    %
      xlimiter(42) =  1.724                                                    %
      xlimiter(43) =  1.666                                                    %
      xlimiter(44) =  1.624                                                    %
      xlimiter(45) =  1.624                                                    %
      xlimiter(46) =  1.666                                                    %
      xlimiter(47) =  1.724                                                    %
      xlimiter(48) =  1.689                                                    %
      xlimiter(49) =  1.712                                                    %
      xlimiter(50) =  1.733                                                    %
      xlimiter(51) =  1.752                                                    %
      xlimiter(52) =  1.770                                                    %
      xlimiter(53) =  1.785                                                    %
      xlimiter(54) =  1.799                                                    %
      xlimiter(55) =  1.812                                                    %
      xlimiter(56) =  1.822                                                    %
      xlimiter(57) =  1.832                                                    %
      xlimiter(58) =  1.839                                                    %
      xlimiter(59) =  1.845                                                    %
      xlimiter(60) =  1.850                                                    %
      xlimiter(61) =  1.853                                                    %
      xlimiter(62) =  1.855                                                    %
      xlimiter(63) =  1.855                                                    %
      xlimiter(64) =  1.854                                                    %
      xlimiter(65) =  1.851                                                    %
      xlimiter(66) =  1.847                                                    %
      xlimiter(67) =  1.908                                                    %
      xlimiter(68) =  1.680                                                    %
      xlimiter(69) =  1.134                                                    %
c                                                                              %
      ylimiter( 1) = -1.300                                                    %
      ylimiter( 2) = -1.000                                                    %
      ylimiter( 3) = -0.900                                                    %
      ylimiter( 4) = -0.800                                                    %
      ylimiter( 5) = -0.600                                                    %
      ylimiter( 6) = -0.400                                                    %
      ylimiter( 7) = -0.200                                                    %
      ylimiter( 8) =  0.000                                                    %
      ylimiter( 9) =  0.200                                                    %
      ylimiter(10) =  0.400                                                    %
      ylimiter(11) =  0.600                                                    %
      ylimiter(12) =  0.800                                                    %
      ylimiter(13) =  0.900                                                    %
      ylimiter(14) =  1.000                                                    %
      ylimiter(15) =  1.300                                                    %
      ylimiter(16) =  1.455                                                    %
      ylimiter(17) =  1.455                                                    %
      ylimiter(18) =  1.455                                                    %
      ylimiter(19) =  1.455                                                    %
      ylimiter(20) =  1.455                                                    %
      ylimiter(21) =  1.455                                                    %
      ylimiter(22) =  1.220                                                    %
      ylimiter(23) =  1.220                                                    %
      ylimiter(24) =  1.170                                                    %
      ylimiter(25) =  1.120                                                    %
      ylimiter(26) =  1.070                                                    %
      ylimiter(27) =  1.020                                                    %
      ylimiter(28) =  0.970                                                    %
      ylimiter(29) =  0.920                                                    %
      ylimiter(30) =  0.870                                                    %
      ylimiter(31) =  0.820                                                    %
      ylimiter(32) =  0.770                                                    %
      ylimiter(33) =  0.720                                                    %
      ylimiter(34) =  0.670                                                    %
      ylimiter(35) =  0.620                                                    %
      ylimiter(36) =  0.570                                                    %
      ylimiter(37) =  0.520                                                    %
      ylimiter(38) =  0.470                                                    %
      ylimiter(39) =  0.420                                                    %
      ylimiter(40) =  0.370                                                    %
      ylimiter(41) =  0.320                                                    %
      ylimiter(42) =  0.320                                                    %
      ylimiter(43) =  0.230                                                    %
      ylimiter(44) =  0.064                                                    %
      ylimiter(45) = -0.064                                                    %
      ylimiter(46) = -0.230                                                    %
      ylimiter(47) = -0.320                                                    %
      ylimiter(48) = -0.320                                                    %
      ylimiter(49) = -0.370                                                    %
      ylimiter(50) = -0.420                                                    %
      ylimiter(51) = -0.470                                                    %
      ylimiter(52) = -0.520                                                    %
      ylimiter(53) = -0.570                                                    %
      ylimiter(54) = -0.620                                                    %
      ylimiter(55) = -0.670                                                    %
      ylimiter(56) = -0.720                                                    %
      ylimiter(57) = -0.770                                                    %
      ylimiter(58) = -0.820                                                    %
      ylimiter(59) = -0.870                                                    %
      ylimiter(60) = -0.920                                                    %
      ylimiter(61) = -0.970                                                    %
      ylimiter(62) = -1.020                                                    %
      ylimiter(63) = -1.070                                                    %
      ylimiter(64) = -1.120                                                    %
      ylimiter(65) = -1.170                                                    %
      ylimiter(66) = -1.220                                                    %
      ylimiter(67) = -1.220                                                    %
      ylimiter(68) = -1.455                                                    %
      ylimiter(69) = -1.455                                                    %
      go to 9090                                                               %
c                                                                              %
c --- d3d limiter points                                                       %
c                                                                              %
 9080 xdim     =  1.70                                                         %
      ydim     =  3.20                                                         %
      redge    =  0.84                                                         %
      greentab = 'default'                                                     %
      ifill    =  1                                                            %
      nlimiter = 33                                                            %
c                                                                              %
      xlimiter( 1) =  1.0139                                                   %
      xlimiter( 2) =  1.0139                                                   %
      xlimiter( 3) =  0.9989                                                   %
      xlimiter( 4) =  0.9989                                                   %
      xlimiter( 5) =  1.1378                                                   %
      xlimiter( 6) =  1.6703                                                   %
      xlimiter( 7) =  1.9055                                                   %
      xlimiter( 8) =  2.1406                                                   %
      xlimiter( 9) =  2.2035                                                   %
      xlimiter(10) =  2.3293                                                   %
      xlimiter(11) =  2.4001                                                   %
      xlimiter(12) =  2.3480                                                   %
      xlimiter(13) =  2.3590                                                   %
      xlimiter(14) =  2.3800                                                   %
      xlimiter(15) =  2.3820                                                   %
      xlimiter(16) =  2.3800                                                   %
      xlimiter(17) =  2.3590                                                   %
      xlimiter(18) =  2.3480                                                   %
      xlimiter(19) =  2.4001                                                   %
      xlimiter(20) =  2.3289                                                   %
      xlimiter(21) =  2.2024                                                   %
      xlimiter(22) =  2.1391                                                   %
      xlimiter(23) =  1.8221                                                   %
      xlimiter(24) =  1.8221                                                   %
      xlimiter(25) =  1.5630                                                   %
      xlimiter(26) =  1.5630                                                   %
      xlimiter(27) =  1.2750                                                   %
      xlimiter(28) =  1.2750                                                   %
      xlimiter(29) =  1.1440                                                   %
      xlimiter(30) =  1.0201                                                   %
      xlimiter(31) =  0.9989                                                   %
      xlimiter(32) =  0.9989                                                   %
      xlimiter(33) =  1.0139                                                   %
c                                                                              %
      ylimiter( 1) =  0.0000                                                   %
      ylimiter( 2) =  0.3998                                                   %
      ylimiter( 3) =  0.4148                                                   %
      ylimiter( 4) =  1.2443                                                   %
      ylimiter( 5) =  1.3832                                                   %
      ylimiter( 6) =  1.3832                                                   %
      ylimiter( 7) =  1.1930                                                   %
      ylimiter( 8) =  1.0027                                                   %
      ylimiter( 9) =  0.8501                                                   %
      ylimiter(10) =  0.5449                                                   %
      ylimiter(11) =  0.3922                                                   %
      ylimiter(12) =  0.3260                                                   %
      ylimiter(13) =  0.2220                                                   %
      ylimiter(14) =  0.1120                                                   %
      ylimiter(15) =  0.0000                                                   %
      ylimiter(16) = -0.1120                                                   %
      ylimiter(17) = -0.2220                                                   %
      ylimiter(18) = -0.3260                                                   %
      ylimiter(19) = -0.3922                                                   %
      ylimiter(20) = -0.5458                                                   %
      ylimiter(21) = -0.8529                                                   %
      ylimiter(22) = -1.0064                                                   %
      ylimiter(23) = -1.3832                                                   %
      ylimiter(24) = -1.3682                                                   %
      ylimiter(25) = -1.3682                                                   %
      ylimiter(26) = -1.3832                                                   %
      ylimiter(27) = -1.3832                                                   %
      ylimiter(28) = -1.3682                                                   %
      ylimiter(29) = -1.3682                                                   %
      ylimiter(30) = -1.2443                                                   %
      ylimiter(31) = -1.2443                                                   %
      ylimiter(32) = -0.4148                                                   %
      ylimiter(33) = -0.3998                                                   %
*                                                                              @
c ----------------------------------------------------------------------       @
c --- BOUNDARY CONDITIONS                                                      @
c ----------------------------------------------------------------------       @
c% deltat_fixed_boundary(mxtbcmhd)                                             @
c                time intervals (in sec) between fixed boundary                @
c                equilibrium calculations. Note that this input is valid       @
c                only for the fixed boundary case( free boundary and           @
c                tdem cases must input the times explicitely in timeqbcd)      @
c                the first equilibrium is at time0 (but see irguess which      @
c                could eliminate the first equilibrium time) the last          @
c                equilibrium is at time timmax PROVIDED that mxtbcmhd          @
c                times deltat_fixed_boundary plus time0 is large enough        @
c                to span this time interval. (if it is not an error exit is    @
c                taken) . If some values are set in timeqbcd then              @
c                deltat_fixed_boundary(1) becomes the first interval           @
c                after the last time given in timeqbcd.                        @
c                default deltat_fixed_boundary(j)=0.0,j=1...mxtbcmhd           @
c                                                                              @
c% timeqbcd(m)   the times for equilibrium                                     @
c                calculations can be (must be if deltat_fixed_boundary         @
c                is not set ,see above)  input in this vector.                 @
c                timeqbcd is a                                                 @
c                list of times (sec) for time-dependent boundary               @
c                conditions on the magnetic flux.  The maximum                 @
c                number of times allowed is mxtbcmhd.                          @
c                if mhdmode ='no coils' then this is a list of times           @
c                at which equilibrium calculations are to be done              @
c                (and possibly eqdsks written). the code may do                @
c                additional MHD calculations, if convergence is                @
c                not achieved.                                                 @
c          NOTE: if mhdmethd ='tdem' (i.e., the multiple eqdsk case)           @
c                then this list of times is not used                           @
c                (the netCDF file  supplies the data instead)                  @
c                                                                              @
c% flxeqbcd(i,m) Flux boundary conditions for MHD equilibrium.                 @
c                index m corresponds to time timeqbcd(m) and index i           @
c                corresponds to the psi loop values if mhdmode = coils         @
c                For mhdmode = 'no coils', flxeqbcd is not used.               @
c                                                                              @
c            The no coils option solves a fixed boundary equilibrium           @
c            problem. the plasma boundary is read from the eqdsk,the           @
c            mhdgrid (rmhdgrid,zmhdgrid) is calculated from the eqdsk          @
c            and initially psincbcd is set to the boundary values of           @
c            psi on the MHD grid as follows:                                   @
c                                                                              @
c                    nh > > > > > nh+nw-1                                      @
c                     ^           v                                            @
c                     ^           v                                            @
c                     ^           v                                            @
c                     ^           v                                            @
c                     1 < < < < < nh+nh+nw-1                                   @
c                                                                              @
c        That is, psincbcd has the flux (in v-sec/rad) stored                  @
c        (for time point 1) starting at the lower left corner                  @
c        and progressing clockwise. there are                                  @
c        2*(nh+nw-2) values of psi stored in psincbcd in all.                  @
c         after the initial time, psincbcd is internally set to                @
c         whatever values are required in order to force the given             @
c         plasma boundary to be a flux surface. thus no boundary               @
c         condition input is required from the user.                           @
c             Point 1 above corresponds to (rmhdgrid(1),zmhdgrid(1)),          @
c            point nh+nw-1 corresponds to (rmhdgrid(nw),zmhdgrid(nh)) etc.     @
c            (if you want to change the number of grid points from nw,nh       @
c            to something else you will have to recompile the code. BUT        @
c            nh must still be 2**m+1 for some integer m, nw is arbitrary).     @
c% curmax(i) maximum current (amps) that coil number i can live with.          @
c            (ignored in this version of the code.)                            @
c ----------------------------------------------------------------------       @
c                                                                              @
 9090 icoil     = 0                                                            %
      nfcoilmax = nfcoil                                                       %
      if (nflxeqbcd .lt. nfcoil)  nfcoilmax = nflxeqbcd                        %
      do 2610 i=1,nfcoilmax                                                    %
      curmax(i) = 1.0e30                                                       %
      do 2610 j=1,mxtbcmhd                                                     %
 2610 flxeqbcd(i,j) = 0.0                                                      %
      do 2620 i=1,mxtbcmhd                                                     %
         deltat_fixed_boundary(i)=0.0                                          %
 2620 timeqbcd(i)   = -1.e30                                                   %
*                                                                              @
c ----------------------------------------------------------------------       @
c --- SOLUTION CONTROL PARAMETERS                                              @
c ----------------------------------------------------------------------       @
c% iteq       Maximum number of inner iterations to obtain a free              @
c             boundary equilibrium for a given current density                 @
c             (suggest 55)                                                     @
c% itcur      Maximum number of outer iterations to obtain a self-             @
c             consistent current density, i.e., to find an ff' consistent      @
c             with the p' and <j*Ro/R> obtained from transport                 @
c             (suggest 3)                                                      @
c% toleq      Tolerance for inner iteration loop (suggest 1.0e-5)              @
c% tolcur     Tolerance for outer, current iteration loop (suggest 1.0e-3)     @
c% omeq       Relaxation parameter for inner iterations (suggest 0.8)          @
c             set omeq = 0.0 to get dynamic relaxation (recommended)           @
c% omcur      Relaxation parameter for outer, current iterations               @
c             omcur is used in subroutine PFPRIM (cray209.f) to relax          @
c             pprim. pprim is defined over the psir(1...nj) grid               @
c             (suggest 0.5)                                                    @
c% isym       1: Symmetrize solution                                           @
c             0: Do not symmetrize solution                                    @
c             note: at the present time the code makes                         @
c                      the solution up/down symmetric about the grid point     @
c                      nh/2+1 by first calculating the total solution and      @
c                      then overwriting the lower half of the solution with    @
c                      a copy of the upper half.                               @
c% npsi       Number of flux surfaces on which the averages are done           @
c             (suggest 51, the same as standard transport grid, max=kpsi       @
c% ivertsbl         switch for vertical stability                              @
c                   set ivertsbl = 0 to disable the internal coding            @
c                   that tries to stop vertical drift of the plasma            @
c                   if it is sensed to occur. set ivertsbl = 1 if vertical     @
c                   drift is a problem and stabilization may help.             @
c                      (this is a real problem for low current highly          @
c                       elongated plasmas. EFIT also sees this drift)          @
c                   The method used can be found in subroutine FREEBDRY.       @
c                   IVERTSBL = 1 HAS AN EFFECT ONLY IF FIXFCOIL=0              @
c                     USE ONLY IVERTSBL = 0 UNTIL FURTHER NOTICE. HSJ          @
c% tensionspl      tension spline parameter. if not equal to 0.0               @
c                  then a tension spline will be used to interpolate           @
c                  q from the MHD to the transport grid. to get                @
c                  infinite tension set tensionspl = 1.0e30. otherwise         @
c                  reasonable values range from 1 to 2000.                     @
c                  Set tensionspl to a negative number, to let the code        @
c                  decide what tension to use.                                 @
c                  NOTE: normally you want to use tensionspl = 0               @
c                        tensionspl is intended only for problem cases where   @
c                        q should be monotonic but isnt due to rapid variaton  @
c                        of the spline representation of q near the boundary.  @
c% extend_seval    integer,set to 1 to allow spline interpolation outside      @
c                  the interpolating table in subroutine seval (not a good     @
c                  idea in general) At this time seval is used only to         @
c                  interpolate beam stopping cross sections for the ADAS       @
c                  cross section set !!!! Any attempt to go outside the        @
c                  range of the table will cause the limiting value in         @
c                  the table to be returned. THERE IS NO EXTRAPOLATION OF THE  @
c                  TABLE !!!!!                                                 @
c                  Default is 0 which stops the code when an out of            @
c                  bounds interpolation is performed.                          @
c                                                                              @
c ----------------------------------------------------------------------       @
c                                                                              @



      extend_seval = 0                                                         %
      iteq       = 60                                                          %
      itcur      =  3                                                          %
      ivertsbl   =  0                                                          %
      toleq      =  1.0e-6                                                     %
      tolcur     =  1.0e-4                                                     %
      omeq       =  0.75                                                       %
      omcur      =  0.5                                                        %
      isym       =  0                                                          %
      npsi       = kpsi    ! parameter kpsi is set in MHDPAR                   %
      tensionspl =  0.0    ! don't use tension spline interpolation            %
*                                                                              @
c ----------------------------------------------------------------------       @
c --- DYNAMIC CONTROL PARAMETERS                                               @
c ----------------------------------------------------------------------       @
c% ieqmax       Number of additional equilibria to be calculated, at times     @
c               given by deqlst(j),j = 1,2..ieqmax . for the fixed             @
c               boundary case (mhdmode = 'no coils') ieqmax should             @
c               be set to zero. for the free boundary case                     @
c               (i.e., mhdmode ='coils') ieqmax may be zero or not             @
c               (see below).                                                   @
c% deqlst(i)    i = 1,2...ieqmax. list of times(sec) for additional equilibria @
c               calculations.                                                  @
c               If no additional equilibria are required then                  @
c               deqlst need not be set and ieqmax should be set to 0.          @
c                                           NOTE                               @
c                    equilibria are always calculated at times                 @
c                    given by timeqbcd(i), i = 1,2..itbcmhd.                   @
c                    deqlst provides a way to get additional equilibria        @
c                    at times other than those at which the MHD                @
c                    boundary conditions are specified. If ibypass             @
c                    is set to 0 then the code may perform additional          @
c                    equilibrium calculations because the time step            @
c                    between equilibria may be reduced due to lack             @
c                    of convergence. an absolute maximum of keqmax             @
c                    equilibria may be calculated. if keqmax is exceeded       @
c                    the code will stop.                                       @
c                    (keqmax is not an input, it is a fixed parameter).        @
c                    for the fixed boundary case (i.e., mhdmode = 'no coils')  @
c                    deqlst is redundant and hence is not used! instead        @
c                    timeqbcd is used to determine when equilibria are         @
c                    to be calculated (again the code may calculate            @
c                    equilibria at additional times if necessary).             @
c                    Note that in the fixed boundary case timeqbcd is          @
c                    the time at which equilibria are calculated rather        @
c                    than the time at which boundary conditions are given      @
c                    and (consequently) equilibria are calculated.             @
c                   An eqdsk with name g.time may be written for each          @
c                   of the times given in timeqbcd and in deqlst.              @
c                   These eqdsks may be used                                   @
c                   to restart the transport code, see restart                 @
c                   parameters, and setting of ieqdsk.                         @
c                   The eqdsks generated are the same as EFIT eqdsks           @
c                   but contain additional information appended to the         @
c                   the normal EFIT eqdsk.                                     @
c% delcap       Maximum relative change allowed in geometric factors           @
c               before and after an equilibrium calculation                    @
c               (suggest 0.10)                                                 @
c% renormalize_rbp  an integer, either 0 (the default) or 1.                   @
c                   If set to 1 then the rbp profile is not changed between    @
c                      successive transport/mhd calls.                         @
c                   If set to 0 then the rbp profile is recalculated from the  @
c                      information just generated by the mhd equilibrium       @
c                      calculation.                                            @
c               This new rbp profile then serves as the initial condition for  @
c               Faraday's law in the transport section of the code. This may   @
c               cause discontinuous jumps in the q profile and current density @
c               due to the fact that the mhd/transport coupling is imperfect.  @
c               If renormalize_rbp = 1 then  rbp is not changed but the various @
c               metrics and flux surface averages are still determined and     @
c               used on the subsequent transport step.                         @
c rbsaxis       a normalized rho value for extrapolating the bootstrap current.@
c               If this value is gt 0.0 then for 0<= rhonormalized             @
c               <= rbsaxis the bootstrap current is replaced by its parabolic  @
c               extrapolation from rbsaxis to the magnetic axis.               @
c% xi_include   logical, set to .true. if eps, xhm2, xi11, xi33, and xips      @
c               are to be included under delcap convergence criterion.         @
c% delrho       Maximum relative change allowed in the value of rho            @
c               at the plasma boundary  before and after an equilibrium        @
c               calculation (suggest 0.005)                                    @
c% maxitr       Maximum number of passes made over a single transport/mhd      @
c               time step, to converge rho on plasma                           @
c               boundary and geometric parameters used in transport, suggest 3 @
c% ibypas       Bypass flag                                                    @
c               0: equilibrium time step is halved down if convergence         @
c                  criteria on delcap and delrho are not met after maxitr      @
c                  passes                                                      @
c               1: convergence test is bypassed after maxitr passes            @
c               NOTE: RUNNING WITH IBYPAS = 1 SHOULD BE DONE ONLY              @
c                  WITH THE INTENTION OF "LETS SEE WHAT HAPPENS".              @
c                  IT DOES NOT IMPLY THAT THE RESULTS WILL BE CORRECT !        @
c% dtmine       Minimum time step between equilibrium calculations             @
c% psifctr      controls the outermost psi surface value at which              @
c               FLUXAV does the contour integral calculations.                 @
c               set psifctr = 0.0 to use the boundary value of psi.            @
c                     {Note that this boundary value of psi                    @
c                      could be different than                                 @
c                      the original value of psilim                            @
c                      on the eqdsk because pol_flux_lim,                      @
c                      see below, was used.}                                    @
c               psifctr = 0.01 means boundary will be pulled in                @
c               0.01*(delta psi). where delta psi is (psimag-psilim)/(npsi-1)  @
c               since npsi is an input param you can use it in conjunction     @
c               with psifctr to determine the outermost flux surface to use.   @
c               reasonable values for psifctr are 0.0 to about 0.9 (if         @
c               you set psifctr = 1 then the last two values of psi on         @
c               the grid become identical so the code won't let you            @
c               set psifctr greater than 0.99)                                 @
c                                                                              @
c% pol_flux_lim A switch introduced so that onetwo could match other codes     @
c               in the treatment of the plasma boundary. The value of psi      @
c               at the plasma boundary is taken as                             @
c               psilim = psiaxis + pol_flux_lim*(psilim_eqdsk-psiaxis)         @
c               where psiaxis is the mag axis psi value, psilim_eqdsk is       @
c               the value of the plasma boundary as specified in the eqdsk.    @
c               (This switch is not availble in tdem mode.)                    @
c               Default value is pol_flux_lim = 1.0                            @
c% derwght      a relaxation parameter for the cap (i.e., fcap, gcap,          @
c               hcap) parameters. derwght is defaulted to -0.5.                @
c               which means 1/2 old plus 1/2 new estimate is                   @
c               taken. to get mostly old estimate use                          @
c               derwght = -0.95, to get mostly new use derwght = -0.05         @
c               any value between -1.0 and 0.0 is valid.                       @
c% cap_mult     a multiplier for the time derivatives                          @
c               dgdt,dhdt,dfdt,dr2dt,drcapdt,dr2idt,depsdt, rcapi,                @
c                  dxhm2dt,dxi11dt,dxi33dt,dxipsdt,drhoa_dt_geom               @
c               most users should leave cap_mult set at its default            @
c               value of 1.0 (it is used for code verification purposes only)  @
c                                                                              @
c% adotmult     multiplier for adot. defaulted to 1.0                          @
c                                                                              @
c% f2d1mult                                                                    @
c% f2d2mult                                                                    @
c% f2d3mult     multipliers for 2d sources in Faraday's Law                    @
c                                                                              @
c% qe2dmult                                                                    @
c% qi2dmult     the electron and ion 2d heating source multipliers             @
c% dlnhdtmult   multiplier for d(ln H)/dt                                      @
c% dlnhdrmult                  d(ln H)/dr                                      @
c% dlnrdtmult                  d(ln r)/dt                                      @
c                                                                              @
c% use_cnt1     logical                                                        @
c% use_cnt2     logical: set only one of use_cnt1 or use_cnt2 to true.         @
c               set the other to false.                                        @
c               if use_cnt1 is true then the modified CNTOUR routine,          @
c               which finds boundary points using steps perpendicular          @
c               to the gradient of psi, is used.                               @
c               if use_cnt2 is true then the ray search method is used.        @
c ----------------------------------------------------------------------       @
c                                                                              @
      use_cnt1     = .true.                                                    %
      use_cnt2     = .false.                                                   %
      timlong      =  1.0d100                                                  %
      ieqmax       =  1                                                        %
      do 2520 i=1,keqmax                                                       %
 2520 deqlst(i)    = timlong                                                   %
      maxitr       =  2                                                        %
      delcap       =  0.01                                                     %
      renormalize_rbp = 0                                                      %
      delrho       =  0.001                                                    %
      derwght      = -0.5                                                      %
      cap_mult     =  1.0                                                      %
      adotmult     =  1.0                                                      %
      f2d1mult     =  1.0                                                      %
      f2d2mult     =  1.0                                                      %
      f2d3mult     =  1.0                                                      %
      qe2dmult     =  1.0                                                      %
      qi2dmult     =  1.0                                                      %
      dlnhdtmult   =  1.0                                                      %
      dlnhdrmult   =  1.0                                                      %
      dlnrdtmult   =  1.0                                                      %
      ibypas       =  0                                                        %
      xi_include   = .false.                                                   %
      psifctr      =  0.0                                                      %
      pol_flux_lim =  1.0                                                      %
*                                                                              @
c ----------------------------------------------------------------------       @
c --- STARTUP/RESTART PARAMETERS                                               @
c ----------------------------------------------------------------------       @
c% irguess    specifies the source of the initial equilibrium                  @
c             AND SOURCE OF INITIAL CONDITIONS FOR TRANSPORT EQUATIONS.        @
c         0:  initial equilibrium is not available. Therefore                  @
c             construct this equilibrium from scratch,                         @
c             using current and pressure profiles in the                       @
c             first namelist of inone.(An internal guess routine               @
c              is used to generate the initial equilibrium. The                @
c              starting shape for psi is controlled with rmagax,               @
c              zmagax,rcsale,zscale,see below. The initial conditions          @
c              for the transport equations are obtained from file inone        @
c              and are used to generate the current profile required           @
c              to solve the Grad-Shafranov equation. An iterative solution     @
c              of the Grad-Shafranov equation is required because the initial  @
c              transport profiles are given as a function of rho               @
c              which is not known until psi(r,z) is established. Psi           @
c              also depends on boundary conditions supplied for the            @
c              Grad-Shafranov equation (in the third namelist of inone) so in  @
c              effect we are trying to force the psi to rho mapping            @
c              to be consistent with the current profile and the MHD           @
c              boundary conditions (by forcing the poloidal current            @
c              function f(psi) to the required values). This may or            @
c              may not lead to a converged solution. If convergence            @
c              is not achieved you may have to alter the current profile       @
c              in inone and/or the boundary conditions for the Grad-Shafranov  @
c              equation. It is also possible (but not likely) that             @
c              a different initial guess for psi will lead to convergence,     @
c              (see rmagax,zmagax,rscale,zscale below).                        @
c              DO NOT USE IRGUESS = 0 UNTIL FURTHER NOTICE! HSJ                @
c        = 1: read initial equilibrium data from file eqdskin                  @
c                  and do initial equilibrium calculation.                     @
c                    (this should reproduce the equilibrium read in,           @
c                     and serves as a check on the input equilibrium)          @
c             irguess =1 will cause initialization from statefile              @
c             instead of an eqdsk if initialize_from_statefile = .TRUE.        @
c        = 2: read p prime                                                     @
c                     and f-fprime from eqdsk. construct curden from           @
c                     this pprime and ffprime. then use kinetic data           @
c                     from first namelist of inone to get new p-prime.         @
c                     use new p-prime together with curden to generate         @
c                     an f-fprime that is consistent with curden               @
c                     calculated from eqdsk and kinetic pressure profiles.     @
c        = -1: read initial equilibrium data from file eqdskin;                @
c                  use input equilibrium for initial time point,               @
c                  and read the parameters rmajor, btor, and totcur from       @
c                  eqdskin if their values are zero in the first               @
c                  namelist . If the values are not zero in the                @
c                  first namelist then rescale the initial equilibrium         @
c                  accordingly. Use this initial equilibrium to generate       @
c                  the MHD dependent geometric coefficients for the            @
c                  transport equations directly.                               @
c                  The  initial conditions required for transport equations    @
c                  are obtained from the first namelist of inone.              @
c                  (this is the method which was used in the                   @
c                  past to make time-independent 1.5 d runs,WITH THE           @
c                    FOLLOWING IMPORTANT EXCEPTION:                            @
c                    The current profile,curden,is no longer an input          @
c                    quantity. The current profile is determined from          @
c                    the eqdsk so that it is consistent with the geometric     @
c                    parameters (mostly gcap). In the old version the eqdsk    @
c                    current profile was used only if you specified the        @
c                    input current profile as a parabola. If you used          @
c                    spline input for curden,the eqdsk current profile         @
c                    given by pprim and ffprim was not used.                   @
c                                                                              @
c                  If irguess = 0 you have some limited control over the       @
c                  intial guess used for psi,by specifying the following       @
c                  parameters:                                                 @
c% rmagax          rmagax and zmagax are guesses for the location of           @
c% zmagax          the magnetic axis (in cm,must be consistent with            @
c                  definition of mhdgrid.)                                     @
c% rscale          rscale and zscale are scaling factors for the               @
c% zscale          intial psi guess.                                           @
c                  The initial shape for psi is                                @
c                                                                              @
c              psi(R,z) = a * EXP (-(R-RMAGAX)**2/RSCALE-(Z-ZMAGAX)**2/ZSCALE) @
c                                                                              @
c                  the constant a is calculated by the code (so that           @
c                  it is consistent with totcur)                               @
c                if rmagax,zmagax,rscale,zscale are not set but are            @
c                required because irguess = 0 the code will pick               @
c                appropriate? values.                                          @
c                                                                              @
c% mhdonly           =0 no effect                                              @
c                    =1 run only a single equilibrium calculation, do not      @
c                       do any transport! This option was created primarily    @
c                       to check out MHD calculations in the context of the    @
c                       transport code. It may be useful in certain            @
c                       circumstances since it allows more general             @
c                       specification of current profiles than does EFIT.      @
c                           if a large number of cases with mhdonly = 1        @
c                           are to be run you may wish to make a               @
c                           special version of ONETWO that does not            @
c                           include the transport section. to do this          @
c                           simply compile and link the code using only        @
c                           cray10,cray20 and cray40. the resulting            @
c                           executable will have some unsatisfied external     @
c                           references (which are not used) but will           @
c                           nevertheless run correctly.                        @
c                                                                              @
c% eqdskin    Name of input eqdsk file (default = 'eqdskin')                   @
c             NOTE: If this is a multiple eqdsk run (mhdmethd = 'tdem') then   @
c                   eqdskin is the name of the data file to be read, NOT the   @
c                   name of any particluar eqdsk. The data file is created     @
c                   by the standalone code MEPC. the file name is usually      @
c                  "shot_12.cdf"                                               @
c% use_efit_cntr  set to 0 to find the plasma boundary from the information    @
c                 on the eqdsk.                                                @
c                 set to 1 to use the cntour given in the eqdsk rather than    @
c                 finding it again                                             @
c ----------------------------------------------------------------------       @
c                                                                              @
      irguess = 0                                                              %
      mhdonly = 0                                                              %
      eqdskin = 'eqdskin'                                                      %
      rmagax  = 0.0                                                            %
      zmagax  = 0.0                                                            %
      rscale  = 0.0                                                            %
      zscale  = 0.0                                                            %
      use_efit_cntr =0                                                         %
*                                                                              @

c ----------------------------------------------------------------------       @
c --- VOLUME CONTROL PARAMETERS                                                @
c ----------------------------------------------------------------------       @
c                                                                              @
c  Turn volume control OFF. (LEAVE IT OFF UNTIL FURTHER NOTICE - HSJ)          @
c% volaray(m)    Volume (cm**3) desired at time timeqbcd(m)                    @
c% voladj        Volume adjustment parameter.  A value of 0. is                @
c                suggested initially.  After a restart the last                @
c                value calculated for voladj should be input.                  @
c% volnudge      Change in volume adjustment parameter (suggest 0.01)          @
c% dvadjmax      ???                                                           @
c% tolvol        Maximum relative difference allowed between actual            @
c                and desired volumes.  A value of 0.005 is suggested           @
c                if volume control is desired.  A large value will             @
c% itvol         Maximum number of attempts made to obtain desired             @
c                volume (suggest 6)                                            @
c% limpos        Side on which plasma touches limiter.  Specify 'left'         @
c                or 'right' if volume control is desired.  A value of          @
c                'none' will turn volume control off.                          @
c% itorfluse     set to 1 to calculate the toroidal flux inside                @
c                the outermost flux surface using rectangular area             @
c                elements. This method avoids problems in integration          @
c                near the x point but may introduce some inaccuracy            @
c                due to the coarse [ 33 by 65] grid normally used              @
c                for DIII-D. The usual flux surface average method             @
c                of calculating the total toroidal flux is (essentially)       @
c                independent of the grid size but may be somewhat inaccurate   @
c                due to integration near the x point (the actual contour       @
c                used is pulled away from the x point some small amount)       @
c                For plasmas without an x point itorfluse = 0 (which is        @
c                the default) is recommended. If there is/are x point(s)       @
c                then itorfluse = 1 may be necessary. if in doubt try both!    @
c                (the max value of rho is what is affected)                    @
c                                                                              @
c ----------------------------------------------------------------------       @
c                                                                              @
      limpos    = 'none'                                                       %
      itorfluse = 0                                                            %
      do i=1,10                                                                %
        volaray(i) = 1.0e30                                                    %
      end do                                                                   %
      tolvol   = 1.0e30                                                        %
      voladj   = 0.0                                                           %
      volnudge = 0.01                                                          %
      itvol    = 6                                                             %
      dvadjmax = 0.0001                                                        %
*                                                                              @
c ----------------------------------------------------------------------       @
c --- OUTPUT PARAMETERS                                                        @
c ----------------------------------------------------------------------       @
c% ieqprt    1: produce copious output concerning equilibrium and              @
c               flux surface averages                                          @
c            0: (default) produce limited output                               @
c% iterdb    1: means create a file of ITER database quantities                @
c               A file named "iterdb" will be created with the results at      @
c               all times in it. Note that any previously existing "iterdb"    @
c               file will be destroyed when the new one is created. The file   @
c               will be written at the end of the run if MHD calcs are not     @
c               done, so that only a single time slice is created. If MHD      @
c               calcs are done, then the file will be written at the end of    @
c               each transport/equilibrium cycle.                              @
c               iterdb = 1  is implemented only for  codeid = dee              @
c               NOTE: iterdb =1 WAS CHANGED ON 4/03/07 TO POINT BACK TO        @
C                     VERSION 1.29 OF CRAY204.F WHICH IS NOW CALLED            @
C                     iterdb_v129.f . TO GET THE VERSION WITH                  @
C                     SMOOTHED TRIANGULARITY,ETC. YOU MUST NOW USE ITERDB = 2  @
C
C                     
c            2: USE THIS VALUE TO GET UPPER/LOWER SMOOTHED TRIANGULARITY       @
C               VERSION OF ITERDB STORED AS CRAY204.F                          @
C
c            0: (default) do not create the database file                      @
c           -1: same as 1, but the file is instead named "iterdb.nc", and is   @
c               in portable, self-describing "netCDF" format instead of ASCII  @

c          ==================================================================  @
c          as off onetwo version 4.24 we have the following:
c           iterdb  
c             -1     calls ier_dbase1,cay205.f , old netcdf version
c              0     no iterdb output
c              1     iter_dbase129,iter_dbase_129.f old version 1.29
c              2     iter_dbase ,cray204.f, upper/lower smoothed 
c                    triangularity, etc. OBSOLETE !!!
c              1 AND create_gcnmp_input =.TRUE., iterdb_gcnmp.f90, the statefile for
c                    Onetwo, netcdf and text version of iterdb.
c                    (see write_iterdb_txt). Use this instead of old iterdb = 2
c                    setting.
c
c          ==================================================================  @
c
c create_GCNMP_input true or false, creates Onewo state file for use in
c               stand alone NTCC - GCNMP solver .
c create_XPTOR_input creates special iterdb text file for use in XPTOR.
c              this must be set to .TRUE. if xptor file is wanted.
c              Default is false.
c
c  initialize_from_statefile  true or false,read iterdb file to 
c               initialize onetwo. This file supplies initial and possibly
c               bounday conditions for Onetwo. The normal usage would be to
c               just supply initial conditions.
c               Note that only the iterdb file created with
c               create_GCNMP_input = .true. and iterdb =1 can be used as
c               such a statefile. Other iterdb files do not have sufficient
c               information. This option is currently not available
c               if tdem mode is in effect. 
c               If   initialize_from_statefile = .TRUE. then the following
c               additional input is required:
c statefile_name = character variable gives name
c                   of statefile to read. If the statefile is in
c                   netcdf format it must have the ".nc" extension


                       
c% iterdsc   1: print out a short description of each ITER database            @
c               variable before printing out the values                        @
c               (iterdb must be set of course)                                 @
c            0: (default) don't print the description                          @
c               iterdsc is currently not an option.                            @
c               the code sets it to 1 (to be compatible with                   @
c               read operation; see irwflag)                                   @
c% irwflag      flag used to indicate if read (irwflag=0)                      @
c               or write (irwflag=1) operation is to be performed              @
c               on the file "iterdb" (this option is                           @
c               included so other codes can easily read the file               @
c               using the same subroutine that originally wrote                @
c               the file. At present there is no way to use                    @
c               the information read from file "iterdb"                        @
c               in ONETWO so this flag should be left at 0                     @
c                         -which is the default value-   ).                    @
c% j2prt     Number of flux surface intervals between output                   @
c% iprtit    1: print output at every iteration in the                         @
c               equilibrium/transport cycle                                    @
c            0: print output only for the last iteration in each               @
c               equilibrium/transport cycle                                    @
c% ieqdsk    1: write eqdsk files at times given by timeqbcd and deqlst        @
c               (this could lead to a large number of files)                   @
c            0: do not write eqdsk file                                        @
c           -1: generate only one eqdsk file at the final time                 @
c               (or at the last time a successful step was taken)              @
c% comp_methd_eqdsk 0: use standard method to generate p and f (for print out  @
c                      to eqdsk only)                                          @
c                      integrate xpp and xffp to get pressure and f            @
c                      unfortunately the uniform psi grid spacing              @
c                      required by the eqdsk                                   @
c                      is typically too coarse (if nw=33 or even 65) to        @
c                      get the leveling off of the pressure profile typically @
c                      found near the magnetic axis                            @
c                      when the pressure profile is determined as a            @
c                      function of RHO . typically the first psi point away    @
c                      from the magnetic axis accounts for as many as 10 rho   @
c                      grid points in the same interval so the structure       @
c                      in press(rho) is lost in press(psi). The following      @
c                      calcuations either get the pressure and f (i.e., sp     @
c                      and sf) self-consistently by integrating the            @
c                      interpolated  pprim and ffprim (i.e., xpp and xffp)     @
c                      to form sp and sf if comp_methd_eqdsk =0                @
c                                OR                                            @
c                      sp and sf are found by interpolation of press           @
c                      and fpsi from the rho to the psi grid and then xpp      @
c                      and xffp are determined by differentiating              @
c                      sp and sf, if comp_methd_eqdsk=1                        @
c                      for fine enough grids the two methods should            @
c                      agree quite well but for the grids typically            @
c                      used in mhd calculations there is some difference       @
c                      between the methods.                                    @
c



c                   1: use modified method SEE SUBROUTINE WRTEQDSK FOR         @
c                      DEFINITIONS OF THE TWO THEORETICALLY EQUIVALENT,        @
c                      BUT COMPUTATIONALLY DIFFERENT, METHODS                  @
c% rf_output 0:  discard standard output of TORAY and GAFIT                    @
c            1:    allow standard output of TORAY and GAFIT to flow <- DEFAULT @
c            2: redirect standard output of TORAY and GAFIT to                 @
c               "toray.log"    and "gafit.log"    respectively (  overwriting) @
c            3: redirect standard output of TORAY and GAFIT to                 @
c               "toray.log"    and "gafit.log"    respectively (concatenating) @
c            4: redirect standard output of TORAY and GAFIT to                 @
c               "toray_NN.log" and "gafit_NN.log" respectively, NN = 01,02,... @
c                                                                              @
c ----------------------------------------------------------------------       @
c                                                                              @
      ieqprt    = 0                                                            %
      j2prt     = 5                                                            %
      iprtit    = 0                                                            %
      ieqdsk    = 0        ! changed from 1 to 0 by HSJ on 11/29/94            %
      rf_output = 1                                                            %
      iterdb    = 0                                                            %
      iterdsc   = 0                                                            %
      irwflag   = 0        ! WRITE (rather than READ) ITER database file       %
      read_iterdb_input = ' '
      comp_methd_eqdsk = 0      
      create_GCNMP_input = .TRUE.  ! both old and new iterdb files are written
      create_XPTOR_input = .TRUE.  ! by defualt
      initialize_from_statefile = .FALSE.

c                                                                              %
c ----------------------------------------------------------------------
c user-defined parameters
c ----------------------------------------------------------------------
c
      do j=1,10
        ispare(j) = 0
        aspare(j) = 0.0
      end do
      myid =0
      numprocs=1

c
c ----------------------------------------------------------------------
c read inone (main input control file)
c ----------------------------------------------------------------------
c6

      rewind (unit = nin)
 2650 read   (       nin, 8020    )  (title(i), i=1,10)
 8020 format (10a8)
      if (title(1) .ne. iflag)  go to 2650
c
c --- read the first namelist, named NAMELIS1
c
      ierr =0              !HSJ 05/022/02
      if(myid .eq. 0)
     . read (unit = nin, nml = namelis1)

      if(jboot .eq. -1)then
         jboot = jhirsh         ! user input jhirsh,not jboot 
      else
         jhirsh = jboot         ! user input jboot
      endif
      if(jhirsh .eq. 99  .or. jhirsh .eq. 100)then
         write(nout,'("jhirsh =99,100 not allowed at this time",/,
     .         "please use 95,96,110,111, or 112 instead")')
         write(ncrt,'("jhirsh =99,100 not allowed at this time",/,
     .         "please use 95,96,110,111, or 112 instead")')
         call STOP('subroutine INIT: invalid specification', 279)
      endif
   


c
c  check for elms
c  k_elms is the number of active elms we will encounter
c
       ierr =0 
       k_elms =0 
       do j=1,n_elms
         if(r_elm(j) .lt. 1.0)ibaloo = -1
         if(t_elms(j) .le. timmax)k_elms = j
       enddo 
       if(k_elms .gt. 0)then
         DO j=2,k_elms 
            IF(t_elms(j-1) .lt. t_elme(j-1) .and. 
     .             t_elme(j-1) .lt. t_elms(j))cycle
            ierr = 1
         ENDDO
         IF(.not. (t_elme(k_elms) .gt. t_elms(k_elms) .and.
     .                 ierr .eq. 0 ))THEN
            DO  j=1,k_elms
                 print *,'t_elms,t_elme =',t_elms(j),t_elme(j)
          ENDDO
            print *,'time0,timmax = ',time0,timmax
            print *,'elm times must be ordered consistently'
            CALL STOP('Elm input error ',1)
         ENDIF
       ENDIF


c
c     multipliers for e, ion channels in IFS model
c
      if (dorl_kotch .ne. 0.0) then
        dorl_kotche = dorl_kotch
        dorl_kotchi = dorl_kotch
      end if
      if (include_itb .eq. 1) then
         ierr=0
         if (cfe_mgb .eq. 0 .and. cfe_bgb .eq. 0) then
            write  (nout, 46)
            write  (*   , 46)
   46       format (' ERROR: cfe_mgb and cfe_bgb cant',
     .                     ' simultaneously = 0')
            ierr = 1
         end if
         if (cfi_mgb .eq. 0 .and. cfi_bgb .eq. 0) then
            write  (nout, 47)
            write  (*   , 47)
   47       format (' ERROR: cfi_mgb and cfi_bgb cant',
     .                     ' simultaneously = 0')
            ierr=1
         end if
         if (ierr .eq. 1)
     .     call STOP ('subroutine INIT: invalid specification', 279)
      end if
      n_spec = 0                      ! in numbrs.i
      do j=1,nj
        if (spec_profiles(j) .gt. 0)  n_spec = n_spec + 1
      end do
c
c --- ensure that user has input only 'old' or 'new' for rlw_model
c
      if (rlw_model .ne. 'old' .and. rlw_model .ne. 'new') then
        write  (ncrt, 43) rlw_model
   43   format (/ ' ERROR: The input value of RLW_MODEL is "', a, '".' /
     .            '        This is not a valid selection.'             /
     .            '        Only "old" and "new" are valid, where the'  /
     .            '        default is the latter.'                     )
        call STOP ('subroutine INIT: invalid RLW_MODEL input', 53)
      end if

c
c--check magnetic braking input
c
      nt_mgbr =0   !no of input error field values
      do j = 1,n_mgbr
         if(time_mgbr(j) .gt. toff_mgbr)nt_mgbr = nt_mgbr+1
      enddo

c
c --- print some diagnostic info on time steps:                      HSJ
c
      dtmin_gcnmp = dtmin
      if (dt .lt. dtmin) then
         write  (ncrt, 51)  dt, dtmin
         write  (nout, 51)  dt, dtmin
   51    format (' WARNING: initial time step is less than dtmin ' /
     .           ' time steps after the first one will be at least ' /
     .           ' dtmin in size . dt,dtmin = ',2(1pe14.8))
      end if
      if (dt .gt. dtmax) then
         write  (ncrt, 52)  dt, dtmax
         write  (nout, 52)  dt, dtmax
   52    format (' WARNING: initial time step is greater than dtmax ' /
     .           ' time steps after the first one will be at most ' /
     .           ' dtmax in size . dt,dtmax = ',2(1pe14.8))
      end if
c
      include_adaptive = 0
      do i=1,nion+4
         if (      eps_adaptive(i) .ne. 0.0) include_adaptive = 1
         if (curve_eps_adaptive(i) .ne. 0.0) then
           include_adaptive  = 1
           include_curvature = 1
         end if
      end do
      if (include_adaptive .eq. 1)  imesh = 0
c
     
      if (namelistvb .gt. 0) then
        write (ncrt, '(" read first  namelist,ierr =",i5)')ierr
c        open  (unit = n42, file = 'namelists', status = 'UNKNOWN')
c        write (unit = n42,  nml =  namelis1)
      end if

c
      if (freeze_adaptive .lt. -1.0e29)
     .    freeze_adaptive = timmax + 100.0 ! turn it off
       jac_skip=MAX(jac_skip,1)


      IF ( inenez .NE. 0 .OR. inenez .NE. 1)THEN
         DO j=1,nion
            IF(ABS(density_mult(j)-1.0D0) .GT. 1.e-13)THEN
               write(nout,7601)
               write(ncrt,7601)
               write(nitre,7601)
 7601          Format(2x,'ERROR, density_mult (1..nion) must = 1.0',/,
     .                2x,'if inenez is not equal to 0 or 1')
               CALL STOP ('subroutine INIT: density multiplier error '
     .                                                          , 255)
            ENDIF
         ENDDO
      ENDIF


cjmp.ibm.par_start 
c ----------------------------------------------------------------------       
c --- get  location of cross section data files
c ----------------------------------------------------------------------
c
         call get_xsect_path(ncrt,nout)
c
cjmp.ibm.par_end

c -----------------------------------------------------------------------
c --- read the second namelist, named NAMELIS2
c -----------------------------------------------------------------------
      if(myid .eq. 0)
     . read (unit = nin, nml = namelis2)
      if (namelistvb .gt. 0) then
        write (ncrt, '(" read second namelist")')
c        write (unit = n42, nml = namelis2)
      end if

c
c --- force namelist 2 parameters to be compatible
c
      if (ifus .ne. 0) jtfus = 1 ! turn on print for fusion if necessary
c
c --- get the value of the MCGO_PATH environment variable
c
      mcgo_path = ' '
      if (IABS (iborb) .eq. 3) then
        ibslow = 1                                     ! needed for MCGO
        call GETENV (mcgo_env_path_name, mcgo_path)
        if (LEN_TRIM (mcgo_path) .eq. 0) then
          call STOP ('subroutine INIT: MCGO_PATH not set', 280)
        end if
        print *, ' MCGO_PATH = ', mcgo_path
      end if
c

      if(ichois .gt. -1 )ichoisrt(1) = ichois
      if(idmpsw .gt. -1 )idmpswrt(1) = idmpsw


      if (ibslow .eq. 0)  iddfusb = 0 ! no fancy neutron rate..
c                                     ..if no slowing down
      mcgo_fast_ion_target = 0
      if (fast_ion_target .eq.  -1) then
          fast_ion_target      = 1
          mcgo_fast_ion_target = 1
      end if

      IF(use_P_Nfreya)THEN
          time_dep_beam =1
          !use_ufile = .TRUE.     ! ufile read governed by nlbdat
                                  ! in nubeam namelist input hsj 8/1/2011
      ENDIF 

c     if this is a snapshot run then we use the old
c     time independent beam method, even if user selected time
c     dependent beam method:
      if(timmax .le. time0)then 
           fix_edge_te(:) = kj
           fix_edge_ti(:) = kj
           fix_edge_rot(:) = kj   !  < nj not allowed in snapshot mode
           fix_edge_ni(:)  = kj
          if(time_dep_beam .gt. 0)then
              print *,'time dependent beam is not allowed in SNAPSHOT'
              print *,'MODE. time_dep_beam option will be ignored'
              time_dep_beam = 0
          endif
      endif



c ----------------------------------------------------------------------
      if( nbeams .gt.  kb )then
               ierr = 1
               print *,'input error, nbeams <= ',kb,' required'
               print *,'you have nbeams =',nbeams
      endif
c     check input if time dependent beam option is selected:
      
      time_dep_beam_sec : if ( time_dep_beam .eq. 1)then
            !CALL setup_run_P_Nfreya(time0,0)
            beam_time_init = 0
            n21s =1 ! see comments in nub4.i
c           overide some user input if not set correctly:
            no_fast_losses = .false.
            IF(.NOT. use_P_Nfreya)iterate_beam = .false.
            if(ibcx .eq. 0) no_fast_losses = .true.
            ibslow = 1
            if(beam_thermal_cutoff .eq. -1)then
               beam_pulse_control =0
               if(beam_thermal_speed .le. 0.0d0)then
                  ierr = 1
                  print *,'beam_thermal_speed, cm/sec, must be set'
               endif
            endif
            if( nbeams .gt.  kb )then
               ierr = 1
               print *,'input error, nbeams <= ',kb,' required'
            endif
            if(nsourc > kbs)then
               ierr = 1
               print *,'input error, nsourc <= ',kbs,' required'
            endif
            if(ngauss > ngauss_max)then
               ierr = 1
               print *,'input error, ngauss  <= ',ngauss_max,' required'
               print *,'ngauss_max = ',ngauss_max
               print *,'ngauss = ',ngauss
            endif
            ierr1=0
            do j=1,nbeams
               do i=1,nsourc
                 offset =0.0d0
                 if(i .ge. 2)offset =  source2_phase(j)
                 bstart = MIN (beamon(j),beamon(j)+offset)
                 if( bstart   .lt. time0)ierr1 = 1
               enddo
            enddo
            task = 'ignore file'
            if(ierr1 .eq. 0)go to 7 !no restart file requiredc
c                                    !ignore it if present . Only this option
c                                   !is valid in this version of the code.
              call STOP('sub INIT,beamon before time0 not allowed',0)



c             starting the beam before time0 is not allowed in this version
c             hence the follwoing code is currently not used.  HSJ
            if(beam_init_restart_file == 0 )then
               if(beam_restart_file(1:4) == 'null')then
                  !no file name and beam_init_restart_file .eq. 0
                  !this implies that file is to be ignored.
                  !we must have beamon(i) .ge. time0 for all i
                  !otherwise this is an error
                    task = 'ignore file'
                  if(ierr1 .gt. 0)then
                       print *,' error, beamon > time0 required '
                       print *,'     when beam_init_restart_file = 0'
                       print *,'     and beam_restart_file name is not',
     .                              ' specified'
                       ierr =ierr1
                  endif
               else
                  !file name is present and beam_init_restart_file .eq. 0
                  !read the file,check that time0 matches, if not either
                  !flag as error or do interpolation??
                  !remove leading blanks:
                  task = 'use file'
                  beam_restart_file = adjustl(beam_restart_file)
                  !get actual (as opposed to declared) length
                  beam_restart_file_length = 
     .                        len_trim(beam_restart_file)
                  open (unit = nb_strt, file = beam_restart_file(1:
     .                  beam_restart_file_length),
     .                  status ='OLD',ERR = 10,FORM = 'UNFORMATTED' )
                  go to 11
 10               ierr = 1
                  print *,'could not open file: ',beam_restart_file(1:
     .                                      beam_restart_file_length)
 11               continue
               endif
            else            !beam_init_restart_file .ne. 0
                            !create the beam restart file
               if(beam_restart_file(1:4) == 'null')then
                  print *,'beam_init_restart_file =1 requires that'
                  print *,'beam_restart_file has a valid filename value'
                  ierr =1
               else
c                   !f90 checks for file,deletes old one and/or creats new one:
                    task = 'create file'
                    open (unit = nb_strt, file = beam_restart_file,
     .                  status ='REPLACE',ERR = 12,FORM = 'UNFORMATTED')
                    go to 13
 12                 print *,'error file ',beam_restart_file
                    print *,' was not created'
                    ierr =1
 13                 continue
               endif
            endif
c
c      end of process namelist input file
 7     taskl  = len_trim(task)
       if(ierr /= 0)
     .    call STOP ('subroutine INIT: time dep beam input error', 0)
       call beam_init(task,taskl,beam_restart_file_length,
     .                beamon,beamoff,btime,nbeams,nsourc,time0,
     .                timmax,bctime,beam_end,beam_cycles) ! cray321.f
       
       endif time_dep_beam_sec





       do k = 1,krf                                 
        if(rframp_timeup(k) .ge. 0.0  .or. 
     .                  rframp_timedown(k) .ge. 0.0) then
           turnonp(k) =  0.0                                                     
           turnonc(k) =  0.0 
           rft = rftime(k) - rframp_timeup(k) - rframp_timedown(k)
           if(rft .lt. 0.0)then 
              write(6,'("error in rf ramp input")')
              write(6,
     .       '("rftime must be .ge. rframp_timeup + rframp_timedown")')
              write(6,'("problem with model number ",i5)')k
              ierr = 1
           endif
        endif
       enddo


       if(ierr .eq. 1)
     . call STOP('subroutine INIT: rframp specification  error', 226)






c --------------------------------------------------------------------
c --- set defaults for Toq and read the third namelist, 
c --- named NAMELIS3, then rewind the file
c --------------------------------------------------------------------
c
       call set_toq_default_input
      if(myid .eq. 0)
     .read   (unit = nin, nml = namelis3)
      if (namelistvb .gt. 0) then
        write (ncrt, '(" read third  namelist")')
c        write (unit = n42, nml = namelis3)
c        call giveupus(n42)
c        close (unit = n42)
      end if
c
      rewind (unit = nin)
c
c     set some switches for the user in case they were set incorrectly
c     if the tdem mode is selected
c


c---------------------------------------mhdmethod = tdem start ----------------
      if (mhdmethd .eq. 'tdem') then
        IF(initialize_from_statefile)THEN
           CALL STOP('subroutine INIT: tdem and  initialize_from_statefi 
     .le not yet compatible',228)

        ENDIF
        if (codeid .ne. 'dee') then
c
c         stop code because too many errors otherwise
c
          call STOP ('subroutine INIT: CODEID = dee required for tdem',
     .                                                             226)
        end if
c
        if ((nstark .ne. nstark1) .or.
     .      (magpri .ne. magpri1) .or.
     .      (nflxeqbcd .ne.  kside)) then
          write (nout,'(" ERROR: nstark,magpri,nflxeqbcd " /
     .                   "need to be reset and code recompiled")')
          call STOP ('subroutine INIT: parameter error', 281)
        end if

c
        ieqdsk      =  0
        ifixshap    =  0
        irguess     = -1
        mhdmode     = 'no coils'
        implicit_fh = .false.
        itdem_flag  = 1    ! set for initial reads of eqdskin file
c       get rime,limiter,volume,circum,mag axis, psilim,psisep,psiaxis,
c       totcur,baxis beq, rhoa,and misc fit info:
        call read_tdem (itdem_flag, eqdskin, ishot)
      end if 
c---------------------------------------mhdmethod = tdem end ----------------



      toq_set: if (mhdmethd .eq. 'toq') then
         toleq_toq = toleq
         equiltype = '"'//equiltype(1:LEN_TRIM(equiltype))//'"'
         betafix = '"'//betafix(1:LEN_TRIM(betafix))//'"'
         jbtype  = '"'//jbtype(1:LEN_TRIM(jbtype))//'"'
         fneqdsk = eqdskin(1:Len(fneqdsk))
         fneqdsk = '"'//fneqdsk(1:LEN_TRIM(fneqdsk))//'"'
         isym_toq = isym
         iteqmx = iteq
         if(isym .eq. 1)then
            !symmetric solution specified. For this case updownsym must
            !be set to 'd' or 'u' for Toq (the default setting is 'a')
            if(updownsym .eq. 'a')then
               write(ncrt,1111)
               write(nout,1111)
 1111          format("error, isym and updownsym not",
     .         " consistent. ",/, " must have updownsym ='d' or 'u'",/,
     .          " when isym = 1")
               write(nout,1112)isym,updownsym
               write(ncrt,1112)isym,updownsym
 1112          format(2x,'isym =',i5,2x,'updownsym =',a)
               call STOP('init, updownsym error',1)
            endif
          else
             !asymmetric case
             updownsym ='a'
c             updownsym ='u'   !force toq to run in symmetric mode
c             isym_toq = 1     !force toq to run in symmetric mode
             updownsym = '"'//updownsym(1:LEN_TRIM(updownsym))//'"'
          endif
      end if toq_set

c---------------------------------mhdmethd .eq. 'toq' end --------------


c     make sure eqdskfilename gets set even if eqdsks are not written:
      if(ieqdsk .eq. 0)eqdskfilename = eqdskin



      !if steady state case and mhd is on set ibypass to make sure
      !that no time stepping can occur.
      !Transient slowing down of alphas from fusion should be in
      !asymptotic form as well so force  iaslow =0
      IF(steady_state .lt. 1.e-5)THEN 
         ibypass =1   
         iaslow  =0
      ENDIF
        if (wrebut .gt. 0.0 .and. wshay .gt. 0.0) then
          write (nout,'(" ERROR: wshay and wrebut" /
     .                   "cannot be set simultaneously")')
          call STOP ('subroutine INIT: confinement model error #1', 249)
        end if
c
        if (wrebut .gt. 0.0 .and. include_weiland .gt. 0.0) then
          write (nout,'(" ERROR: wrebut and include_weiland" /
     .                   "cannot be set simultaneously")')
          call STOP ('subroutine INIT: confinement model error #2', 250)
        end if
c
        if (wrebut .gt. 0.0 .and. include_ifs .gt. 0.0) then
          write (nout,'(" ERROR: wrebut and include_ifs" /
     .                   "cannot be set simultaneously")')
          call STOP ('subroutine INIT: confinement model error #3', 251)
        end if
c
        if (wshay .gt. 0.0 .and. include_weiland .gt. 0.0) then
          write (nout,'(" ERROR: wshay and include_weiland" /
     .                   "cannot be set simultaneously")')
          call STOP ('subroutine INIT: confinement model error #4', 252)
        end if
c
        if (wshay .gt. 0.0 .and. include_ifs .gt. 0.0) then
          write (nout,'(" ERROR: wshay and include_ifs" /
     .                   "cannot be set simultaneously")')
          call STOP ('subroutine INIT: confinement model error #5', 253)
        end if
c
        if (include_weiland .gt. 0 .and. include_ifs .gt. 0) then
          write (nout,'(" ERROR: include_ifs and include_weiland" /
     .                   "cannot be set simultaneously")')
          call STOP ('subroutine INIT: confinement model error #6', 254)
        end if

        if (include_weiland .gt. 0 .and. include_glf .gt. 0) then
          write (nout,'(" ERROR: include_glf and include_weiland" /
     .                   "cannot be set simultaneously")')
          call STOP ('subroutine INIT: confinement model error #6', 254)
        end if
        if (include_ifs .gt. 0 .and. include_glf .gt. 0) then
          write (nout,'(" ERROR: include_glf and include_ifs" /
     .                   "cannot be set simultaneously")')
          call STOP ('subroutine INIT: confinement model error #6', 254)
        end if
cJMP START
        if (include_weiland .gt. 0 .and. include_mmm .gt. 0) then
          write (nout,'(" ERROR: include_mmm and include_weiland" /
     .                   "cannot be set simultaneously")')
          call STOP ('subroutine INIT: confinement model error #7', 254)
        end if
        if (include_ifs .gt. 0 .and. include_mmm .gt. 0) then
          write (nout,'(" ERROR: include_mmm and include_ifs" /
     .                   "cannot be set simultaneously")')
          call STOP ('subroutine INIT: confinement model error #7', 254)
        end if
        if (include_glf .gt. 0 .and. include_mmm .gt. 0) then
          write (nout,'(" ERROR: include_mmm and include_glf" /
     .                   "cannot be set simultaneously")')
          call STOP ('subroutine INIT: confinement model error #7', 254)
        end if
cJMP END
        if (   rho_edge .lt. 0.0 .or.    rho_edge .gt. 1.0)
     .    rho_edge = 1.0
        do j = 1,kbctim
           if (fix_edge_te(j) .lt. kj  .or.
     .              fix_edge_ti(j) .lt. kj .or.
     .              fix_edge_rot(j) .lt. kj)then
              rho_edge = 1.0
           if (nbctim .lt. 2 .and. timmax .gt. time0)      !snapshot mode
                                                           !requires nbctim =1 
     .     call STOP ('subroutine INIT: nbctim < 2 not good here', 282)
           endif
        enddo


      if(include_glf .eq. 1)then
           if(itte_dv .eq. 1)dv_method = .true.
           if(itti_dv .eq. 1)dv_method = .true.
           if(itangrot_dv .eq. 1)then
             dv_method = .true.
             call STOP('dv method not implemented for toroidal rotation'
     .                                                              ,1)
            endif
      else  !make sure these are off even if user has them on inadvertenly
           itte_dv =0
           itti_dv =0
           itangrot_dv =0
      endif

c

c
      if (IABS (iborb) .eq. 3) then
           spawn_mcgo=1
           if (iborb .eq. -3) then
             spawn_mcgo = 0  !used to sense mcgo  file writes 
             do ib = 1,nbeams
               if (LEN_TRIM (mcgo_output_file(ib)) .gt. 0) then
                 read_mcgo_file(ib) = 1
               end if
             end do
           end if
           npart_mcgo = MIN (npart_mcgo, npart, npart_mcgo_max)
           iborb = 3
      end if
c
c --- in case user set print options wrong
c
      if (j2prt .eq. 0)  j2prt = 5
      if ( jprt .eq. 0)   jprt = 5
c
      write  (nout  , 8000)
      write  (nqik  , 8000)
      write  (nitre , 8000)
      write  (nmix  , 8000)
      write  (ntweak, 8000)
 8000 format (/ 1x, 11('-') / 1x, 'input file:' / 1x, 11('-'))
c
 2670 read   (nin   , 8020, end = 2680) (title(i), i=1,10)
      write  (nout  , 8020)             (title(i), i=1,10)
      write  (nqik  , 8020)             (title(i), i=1,10)
      write  (nitre , 8020)             (title(i), i=1,10)
      write  (nmix  , 8020)             (title(i), i=1,10)
      write  (ntweak, 8020)             (title(i), i=1,10)
      go to 2670
c
c ----------------------------------------------------------------------
c --- set error flag to zero and start input checks
c ----------------------------------------------------------------------
c
 2680 psifctr          = MIN (psifctr, point99)
      run_preplt_dummy = run_preplt ! dummy version for argument passing
      tohmwant_dummy   = tohmwant   ! dummy version for argument passing
      ierr = 0
      if (codeid .eq. 'dee' .and. irguess .eq. 0) then
        ierr = 1
        write  (nout, 2687)  irguess
        write  (ncrt, 2687)  irguess
 2687   format (' ERROR: irguess = ', i3, ' is not implemented')
      end if

      nj_ncd = nj
      njm1_ncd = nj-1
      npsi_ncd = npsi
c
c ----------------------------------------------------------------------
c --- check for Hsieh model:
c     wshay has to be set so that subroutine SHAY_CHIE will be called
c ----------------------------------------------------------------------
c
      if (scsmult       ) wshay = 1.0
      if (wshay .ne. 0.0) then
          jsourc = 1             ! make sure subroutine PSOURC is called
          ktimes = 0
          do j=1,maxtimes
            smulta(j)       = smult
            skimulta(j)     = skimult
            smultstder(j)   = 0.0
            skimultstder(j) = 0.0
            slim95(j)       = 0.0
            skilim95(j)     = 0.0
          end do
          if (scsmult) then            ! set some switches for obtaining
            smult        =  1.0      ! the multiplier in Hsieh model
            skimult      =  1.0      ! ion multiplier in Hsieh model
            iterate_beam = .true.
            itte         =  0
            itti         =  0
             srincev  = MAX (srincev , srin )
            sroutcev  = MIN (sroutcev, srout)
            relaxshay = 1.0           ! no under-relaxation in this case
          end if
          if (ishayform .eq.  0  ) skimass = 1.0 ! use old form of model
          if (    smult .eq.  0.0) ierr    = 1
          if (   sbpexp .eq.  0.0) ierr    = 1
          if (     srin .eq. -1.0) ierr    = 1
          if (    srout .eq. -1.0) ierr    = 1
          if (     ierr .eq.  1  ) then
            write  (nout, 2345)
            write  (ncrt, 2345)
 2345       format (' ERROR in Hsieh model specifications')
          end if
      end if
 
c
c ----------------------------------------------------------------------
c  Houlberg bootstrap current model requires that the exact form of the
c  trapped ion fraction be used (ftcalc = 'exact')
c ----------------------------------------------------------------------
c
c gms annotation below is Gary M. Staebler (He commented out the lines)
c
* gms if (jhirsh .ge. 95 .and. ftcalc .ne. 'exact') then
* gms   write  (nout, 2347)
* gms   write  (ncrt, 2347)
*2347   format (
* gms.  ' NOTE: The Houlberg Bootstrap Current model requires that' /
* gms.  '       FTCALC be set to use the exact form of the trapped'
* gms.        ' ion fraction.                                     ' /
* gms.  '       FTCALC is being reset to "exact" for this run!'     )
* gms   ftcalc = 'exact'
* gms end if
c
c ----------------------------------------------------------------------
c  Load Staebler-Hinton model constants into common block variables
c  for later calculations (if ifsflag = 1)
c ----------------------------------------------------------------------
c
      fs( 1) = ifsflag
      fs( 2) = aeh
      fs( 3) = ael
      fs( 4) = aih
      fs( 5) = ail
      fs( 6) = bh
      fs( 7) = bl
      fs( 8) = ch
      fs( 9) = cl
      fs(10) = alfae
      fs(11) = alfai
      fs(12) = betah
      fs(13) = sigma
      fs(14) = gammah
      fs(15) = lsfctr
      fs(16) = xrot
      fs(17) = xeden
      fs(18) = xsecder
c
      if (lsfctr .gt. 0)  staeblrmodl(15) = 'used R0'
c
      if     (resistive .ne. 'hinton'  ) then
        if   (resistive .ne. 'hirshman') then
          if   (resistive .ne. 'nclass') then
c         if (resistive .ne. 'kim'     ) then
            ierr = 1
            write  (nout, 2688) resistive
            write  (ncrt, 2688) resistive
 2688       format (' ERROR: resistive =', a, ' is not allowed')
          end if
        end if
      end if
c
      if (ftcalc .ne. 'exact') then
          if (ftcalc .ne. 'analytic') then
            ierr = 1
            write  (nout, 2689) ftcalc
            write  (ncrt, 2689) ftcalc
 2689       format (' ERROR: ftcalc =', a, ' is not allowed')
          end if
      end if
c
c ---------------------------------------------------------------------
c --- allocate space for glf23 model if necessary
c ---------------------------------------------------------------------
      if((include_glf .eq. 1).or.(include_mmm .eq. 1)) !JMP
     .  call allocate_glf




c----------------------------------------------------------------------
c --- allocate Faraday law related routines;
c----------------------------------------------------------------------
      call allocate_fdy_arrays


c
c ----------------------------------------------------------------------
c --- check for MHD-only run
c ----------------------------------------------------------------------
c
      if (mhdonly .eq. 1 .and. codeid .ne. 'dee') then
        ierr = 1
        write  (nout, 2681)  codeid
        write  (ncrt, 2681)  codeid
 2681   format (' ERROR: codeid =', a8 /
     .          ' codeid must be set to "dee" if mhdonly = 1 is set')
      end if
c
      if (codeid .ne. 'dee')  iterdb = 0  ! turn off iterdb if necessary
c
c --- implicit_fh can at present not be used with onedee
c --- (subroutine RHOMSH must be modified for onedee case)
c
      if (codeid  .eq. 'onedee')  implicit_fh = .false.
      if (mhdonly .eq.   1)  iyoka   =  0  ! prevent main calling PSUMRY
      if (irguess .eq. -99)  irguess = -1  ! if user set old definition
      if (ltest_code .eq. 1) then
          if(u_vloop_bc)CALL STOP(' u_vloop not with ltest_code',1)
          irguess   =  1              ! test using eqdsk equilibrium
          beamon(1)    = timmax + 50.0   ! no beam for test case
          totcur(1) = -1.0e30         ! current to be taken from pcurmhd
          do j=1,mxtbcmhd
            if (pcurmhd(j) .ne. 0.0)  go to 2682
          end do
          ierr = 1
          write  (nout, 2683)
          write  (ncrt, 2683)
 2683     format (' ERROR: pcurmhd must be set for ltest_code option')
 2682     if (mhdonly .eq. 1) then
            ierr = 1
            write  (nout, 2684)
            write  (ncrt, 2684)
 2684       format
     .      (' ERROR: mhdonly must not be set for ltest_code option')
          end if
          zeffin(1,1) = 0.0
          do j=1,nj
            zeff(j)   = 1.0
          end do
          nimp   = 0      ! no impurities for test case
          nbctim = 1      ! kinetic profiles run in analysis mode..
c                         ..are fixed in time
c                           zeff set here, te set after initial MHD calculations
      end if
c
c ----------------------------------------------------------------------
c check times for time-dependent quantities
c ----------------------------------------------------------------------
c
      if (nbctim .le. kbctim)  go to 1954
      write  (nout, 7950) kbctim
      write  (ncrt, 7950) kbctim
 7950 format (' ERROR: nbctim must be less than or equal to ', i8)
      ierr   = 1
      nbctim = 1
 1954 nbctim = MAX0 (nbctim, 1)
      if (nbctim .eq. 1) then ! only one bctime given, must be equal to time0
        bctime(1) = time0
        IF(ABS(bctime(1)-time0) .GT. 1.e-10)THEN ! not active
          write  (nout, 7959)
          write  (ncrt, 7959)
 7959     format (' *** ERROR: bctime(1) must be equal to time0 here')
          ierr = 1
        ENDIF
      else
        if (bctime(1) .le. time0)  go to 1956
        write  (nout, 7956)
        write  (ncrt, 7956)
 7956   format (' *** ERROR: bctime(1) must be <= time0')
        ierr = 1
 1956   if (bctime(nbctim) .lt. timmax) then
          write  (nout, 7958)
          write  (ncrt, 7958)
 7958     format (' *** ERROR: bctime(nbctim) must be .ge. timmax')
          ierr = 1
          bctime(nbctim) = timmax
        end if
      end if
      IF(ierr .ne. 0)CALL STOP('nbctim,bctime input error',1)



c
c ---------------------------------------------------------------------
c check vloop_bc input.
c setup_vloop sets up  vloop_bc and vloop_bc_time arrays,nvloop,
c and u_vloop_bc = .true.  u_vloop_bc_time =.true. if appropriate:
      call setup_vloop(time0,timmax,ncrt,nbctim) !in bd_condtn.f90
c     tdem mode not compatible at this time:
      if( mhdmethd .eq. 'tdem' .and. u_vloop_bc) then
         write(ncrt,177)
 177      format(2x,'Error:  time dependent eqdsk mode and vloop',/,
     .    8x,'boundary condition input are not compatible')
         write(nout,177)
         CALL STOP('vloop bc error',1)
      endif



c
c ----------------------------------------------------------------------
c use eqdsk current by setting totcur(1) = 0.0
c
c ----------------------------------------------------------------------
c
C      if (  time0 .gt. bctime(1) .and.
C     .     codeid .eq. 'dee'     .and.
C    .  totcur(1) .gt. -1.0e30)  totcur(1) = 0.0
C      totcur1 = totcur(1)
C      print *,'totcur(1) =',totcur(1)  replaced with following, 5/16/05 : HSJ

       if (  time0 .gt. bctime(1) .and.
     .           codeid .eq. 'dee'     .and.
     .                totcur(1) .gt. -1.0e30)then
          if(nbctim .le. 1 )then   !assume user wants eqdsk current
                totcur(1) = 0.0
          else !time interpolate the current      
              if(bctime(nbctim) .lt. time0)then
               write(nout,867)time0,bctime
               write(nqik,867)time0,bctime
               write(ncrt,867)time0,bctime
               write(nitre,867)time0,bctime
 867           FORMAT(2x,'ERROR in input specification of ',/,
     .                2x,'time dependent total current ',/,
     .                2x, 'initial time is :',2x,1pe12.6,/,
     .                2x, 'bctime does not embrace this time ',/,
     .                2x, 'bctime = ',(5(2x,1pe12.5)),
     .                2x, 'and hence the total current cannot',/,   
     .                2x, ' cannot be determined')
                      CALL STOP('init, error with totcur and bctime',1)
              else
                   do m=2,nbctim
                    if (bctime(m) .ge. time0)then
                      totcur1 = (totcur(m)-totcur(m-1))*
     .                     (time0-bctime(m-1))/(bctime(m)-bctime(m-1))
     .                         + totcur(m-1)
                          EXIT
                     endif
                   end do
               endif
              totcur(1) = totcur1
          endif
        endif
       totcur1 = totcur(1)
   



c
c ----------------------------------------------------------------------
c --- fully implicit solution is not consistent with event flag logic.
c --- therefore avoid theta = 1.0:
c
      theta = MIN (theta, point999999)
c
c ----------------------------------------------------------------------
c  set up time-dependent boundary conditions for equilibrium calculations
c  (some results obtained here are used below,
c  so don't move this section from here)
c ----------------------------------------------------------------------
c

      cursrce  = 'inone'
      btorsrce = 'inone'         ! specifies source of current and btor
      if (totcur(1) .eq.  0.0   )  cursrce = 'eqdsk'
      if (totcur(1) .lt. -0.9e+30)  cursrce = 'pcurmhd'
      if (mhdmethd  .eq. 'tdem' )  cursrce = 'eqdsk'
      if (u_vloop_bc) cursrce = 'internl'
      if (ifixshap .eq. 1 .and. cursrce .eq. 'pcurmhd') then
        write  (ncrt, 6273)
        write  (nout, 6273)
 6273   format (' ERROR: for ifixshap = 1, totcur must be' /
     .          ' specified for nbctim values of bctime'   /
     .          ' in the first namelist of inone')
        ierr = 1
      end if
      if (    btor .eq.  0.0   )  btorsrce = 'eqdsk'
      if (    btor .lt. -0.9e+30)  btorsrce = 'btormhd'
      if (mhdmethd .eq. 'tdem' )  btorsrce = 'eqdsk'
c
c     start codeid branch
c
      if (codeid .eq. 'dee') then
c
c     start ifixshap branch
c
      mhd_methd: if (ifixshap .ne. 1 .and. mhdmode .eq. 'coils') then
c
        if ((nstark .ne. nstark1) .or.
     .      (magpri .ne. magpri1) .or.
     .      (nflxeqbcd .ne.  kside)) then
          write (nout, '("ERROR: nstark,magpri,nflxeqbcd " /
     .                   "need to be reset and code recompiled")')
          call STOP ('subroutine INIT: parameter error', 283)
        end if
c
c --- time-dependent equilibrium case:
c --- get the number of times at which boundary conditions are specified:
c
        do j=mxtbcmhd,1,-1
          i = j
          if (timeqbcd(j) .ne. 0.0)  go to 6255
        end do
c
 6255   itbcmhd = i
        if (itbcmhd .eq. 1) then
             write (nout, 6278)
             write (ncrt, 6278)
 6278        format
     .       (' WARNING THIS CASE REQUIRES THAT MHD BOUNDARY'         /
     .        ' CONDITIONS BE GIVEN AT MORE THAN ONE TIME. CODE WILL' /
     .        ' ASSUME BOUNDARY CONDITIONS ARE CONSTANT IN TIME.')
          timeqbcd(itbcmhd) = time0
          itbcmhd           = 2
          timeqbcd(itbcmhd) = timmax
        else
          if (timeqbcd(itbcmhd) .lt. timmax) then
            itbcmhd = MIN0 (itbcmhd+1, mxtbcmhd)
            write  (nout, 6277) itbcmhd
            write  (ncrt, 6277) itbcmhd
 6277       format (' WARNING TIMEQBCD(',I3,') RESET TO TIMMAX')
            timeqbcd(itbcmhd) = timmax
          end if
          if (timeqbcd(1) .gt. time0) then
             write  (nout, 6276)
             write  (ncrt, 6276)
 6276        format (' WARNING TIMEQBCD(1) RESET TO TIME0')
            timeqbcd(1) = time0
          end if
        end if
c
c --- check and set up the deqlst vector. merge the values in timeqbcd
c --- with the values in deqlst in time-sequenced order:
c
          ij = 0
          do j=1,keqmax
              if (deqlst(j) .ne. timlong)  ij = ij + 1
          end do
          if (ij .ne. ieqmax) then
              write  (ncrt, 4500)  ieqmax, ij
              write  (nout, 4500)  ieqmax, ij
 4500         format (' ERROR: ieqmax = '  , i5 /
     .                ' times in deqlst = ', i5 /
     .                ' run will be stopped')
              ierr = 1
          end if
          if (ieqmax .eq. 0) then  ! no additional times given in deqlst
            do j=1,itbcmhd
              deqlst(j) = timeqbcd(j)
            end do
            ieqmax = itbcmhd
          else                ! additional times specified, process them
c
c             first check for duplicate times (i.e., within dtmine)
c
              do j=1,itbcmhd
                 do i=1,ieqmax
                     if (ABS (timeqbcd(j)-deqlst(i)) .lt.  dtmine
     .               .or. deqlst(i) .lt. timeqbcd(1)
     .               .or. deqlst(i) .gt. timeqbcd(itbcmhd)) then
                         ieqmax = ieqmax-1
                         do ij=i,ieqmax
                             deqlst(ij) = deqlst(ij+1)
                         end do
                         go to 4510
                     end if
                 end do
 4510            continue
              end do
c
c           times which were the same within dtmine in deqlst and
c           timeqbcd have now been eliminated from deqlst.
c           next merge deqlst with timeqbcd and sort in ascending
c           (i.e., chronological) order:
c
               ij = ieqmax+itbcmhd
               if (ij .gt. keqmax) then
                  write  (ncrt, 4530)  ieqmax, keqmax
                  write  (nout, 4530)  ieqmax, keqmax
 4530             format (' ERROR: ieqmax =',   i5 /
     .                    ' keqmax = '        , i5 /
     .                    ' must have ieqmax .le. keqmax')
                  ierr = 1
               end if
               do j=ieqmax+1,ij
                 deqlst(j) = timeqbcd(j-ieqmax)
               end do
               ieqmax = ij
               call my_vsrta (deqlst, ieqmax)     ! sort by algebraic value
c
c            deqlst(i),i = 1,ieqmax now represents a chronological
c            list of times at which equilibria will be calculated
c            (deqlst(1) = timeqbcd(1),deqlst(ieqmax) = timeqbcd(itbcmhd)
c            and timeqbcd is a subset of deqlst.
c            if the time steps between equilibria are too large the
c            code may cut these times steps down,resulting in more
c            equilibrium calculations )
c
              write (nout, 4540)  (deqlst(j), j=1,ieqmax)
              write (ncrt, 4540)  (deqlst(j), j=1,ieqmax)
 4540 format ('  initial list of times for equilibrium calculations' /
     .          (5(2x, 1f12.6)))
          end if
          if (mhdmethd .eq. 'cycred' .or.
     .        mhdmethd .eq. 'green')  go to 6337
            ierr = 1
            write  (nout, 6338) mhdmethd
            write  (ncrt, 6338) mhdmethd
 6338       format (' ERROR: mhdmethd =', a, ' is not valid here')
 6337     continue

      else if ((ifixshap .ne. 1 .and.   mhdmode .eq. 'no coils')
     .                 .and. (mhdmethd .ne. 'tdem')) then  mhd_methd 

c
c --- fixed boundary MHD calculations (mhdmode .eq. 'no coils'
c --- and ifixshap = 0,mhdmethd ne 'tdem' ).  
c --- check the times at which 
c     equilibria are to be
c --- generated. if times are not given then create them now.
c
c --- for this case it is possible to just specify the time interval between
c --- between equilibrium calculations, deltat_fixed_boundary. Times may also
c --- be explicitly given in timeqbcd if desired.
c
          if (deltat_fixed_boundary(1) .gt. 0.0) then ! user didn't..
c                                                     ..specify timeqbcd
c            expand deltat_fixed_boundary
c
             do j=mxtbcmhd,1,-1 !find the lasst input value
               jj =j
               if (deltat_fixed_boundary(j) .gt. 0.0)  go to 6133
             end do
c
             call STOP ('subroutine INIT: deltat_fixed_boundary', 284)
c
 6133        do j= jj+1,mxtbcmhd !fill remainder with last input value
                deltat_fixed_boundary(j)=deltat_fixed_boundary(jj)
             end do
c
c            check if any times were specified:
c
             do j=mxtbcmhd,1,-1 !find the lasst input value
                jj =j
                if (timeqbcd(j) .gt. -1.e30 )go to 6137
             end do
             jj=1  !fell through because no times were input in timeqbcd
             timeqbcd(jj)=time0
 6137        jj=jj+1  !start index for extending timeqbcd
             do j=jj,mxtbcmhd !now setup timeqbcd
                timeqbcd(j)=timeqbcd(j-1)+ deltat_fixed_boundary(j-jj+1)
                if (timeqbcd(j) .ge. timmax)  go to 6134
             end do
             if (timmax .gt. timeqbcd(mxtbcmhd-1)) then
                write (*,'("ERROR:  using deltat_fixed_boundary = ",
     .                 1pe12.8 /
     .                 "will not span the interval [time0,timemax]" /
     .             i5, " maximum intervals are allowed")')
     .                 deltat_fixed_boundary(mxtbcmhd-1),mxtbcmhd-1
               call STOP ('subroutine INIT: ONETWO init errors', 83)
             end if
 6134        continue
          end if
          do j=mxtbcmhd,1,-1
            itbcmhd = j
            if (timeqbcd(j) .gt. -1.e30)  go to 6251
          end do
          itbcmhd = 0
c
 6251     if (itbcmhd .le. 1) then      ! user didn't set, so force them
              write  (nout, 6252) mxtbcmhd
              write  (ncrt, 6252) mxtbcmhd
 6252         format (' ERROR: timeqbcd(i), i=1,2..'                /
     .                '   must be set for this case.'               /
     .                '   A maximum of', i3, ' times may be given.' /
     .                '   The code will calculate equilibria',
     .                  ' at these times.')
              ierr = 1
          else
              call my_vsrta (timeqbcd, itbcmhd)   ! sort in increasing time
              if (timeqbcd(1) .lt. time0   .or.
     .            timeqbcd(1) .gt. timmax) then
                ierr = 1
                write  (nout, 6253)  time0, timeqbcd(1), timmax
                write  (ncrt, 6253)  time0, timeqbcd(1), timmax
 6253           format (' ERROR: timeqbcd(1) must be .ge. time0' /
     .                  ' and less than timmax'                  /
     .                  ' time0 = '      , 1pe12.6               /
     .                  ' timeqbcd(1) = ', 1pe12.6               /
     .                  ' timmax = '     , 1pe12.6)
              end if
              if (timeqbcd(itbcmhd) .lt. time0) then
                ierr = 1
                write  (nout,6254)itbcmhd,time0,timeqbcd(itbcmhd),timmax
                write  (ncrt,6254)itbcmhd,time0,timeqbcd(itbcmhd),timmax
 6254           format (' ERROR: timeqbcd(', i5,
     .                  ') must be .ge. time0'  /
     .                  ' and less than timmax' /
     .                  '    time0 = ', 1pe12.6 /
     .                  ' timeqbcd = ', 1pe12.6 /
     .                  '   timmax = ', 1pe12.6)
              end if
          end if
          ieqmax = itbcmhd
          do j=1,itbcmhd
            deqlst(j) = timeqbcd(j)
          end do
            if (mhdmethd .eq. 'cycred' .or.
     .          mhdmethd .eq. 'sorpicrd' .or.
     .          mhdmethd .eq.'toq' )  go to 6335
            ierr = 1
            write  (nout, 6336) mhdmethd
            write  (ncrt, 6336) mhdmethd
 6336       format (' ERROR: mhdmethd =', a, ' is not valid here')
 6335       continue
      else if (ifixshap .ne. 1 .and. mhdmethd .eq. 'tdem') then 
     .                                                 mhd_methd 
!           time-dependent eqdsk mode
            write (6, '("tdem mode selected,eqdskin=", a)') eqdskin
      else  mhd_methd  ! ifixshap = 1
c
c --- fixed equilibrium case
c     requires boundary conditions if irguess .ge. 0
c
        itbcmhd           = 1
        timeqbcd(itbcmhd) = time0
      end if  mhd_methd
c
c end of ifixshap branch; continue with codeid = dee branch
c
 4543    if ( mhdmode .eq. 'coils' .and. ((ifixshap .eq. 0)  .or.
     .      (ifixshap .eq.  1      .and. irguess .gt. 0))) then
c
c --- adjust data matrices for zero data:
c --- a)the psiloop values:
c

      call expchk (flxeqbcd, kside, mxtbcmhd, nsilop, itbcmhd, xdum)
c
c --- b)the probe values:
c
      call expchk (expmp2, magpri, mxtbcmhd, magpr2, itbcmhd, xdum)
c
c --- c)the fcoil current values:
c
      call expchk (fcoilcur, nfcoil, mxtbcmhd, nfcoil, itbcmhd, xdum)
c
c --- d)the ecoil current values:
c
      call expchk (ecurrt, nesum, mxtbcmhd, nesum, itbcmhd, xdum)
c
c --- e)the vessel current values:
c
      call expchk (vescurrt, nvessel, mxtbcmhd, nvessel, itbcmhd, xdum)
c
c --- further check on the MHD boundary conditions
c
      do j=1,itbcmhd
c
c --- a)the psi-loop values:
c
          sum = 0.0
          if (ifitpsi .eq. 1) then
            do 6275 k=1,nsilop
 6275       sum = sum + flxeqbcd(k,j)
            if (sum .eq. 0.0) then
              if (j .eq. 1)  go to 6285
              do 6295 k=1,nsilop
 6295           flxeqbcd(k,j) = flxeqbcd(k,j-1)
            end if
          end if
c
c --- b)the magnetic probe values:
c
          if (ifitprob .eq. 1) then
            sum = 0.0
            do 6176 k=1,magpr2
 6176       sum = sum + expmp2(k,j)
              if (sum .eq. 0.0) then
                if (j .eq. 1)  go to 6285
                do 6177 k=1,magpr2
 6177             expmp2(k,j) = expmp2(k,j-1)
              end if
          end if
c
c --- c)the f-coil current  values:
c
          if (fixfcoil .eq. 1 .or. fitfcur .eq. 1) then
            sum = 0.0
            do 6180 k=1,nfcoil
 6180       sum = sum + ABS (fcoilcur(k,j))
              if (sum .eq. 0.0) then
                if (j .eq. 1)  go to 6285
                do 6181 k=1,nfcoil
 6181           fcoilcur(k,j) = fcoilcur(k,j-1)
              end if
          end if
c
c --- d)the ecoil current  values:
c
          if (iecurr .eq. 1) then
            sum = 0.0
            do 6186 k=1,nesum
 6186       sum = sum + ABS (ecurrt(k,j)  )
              if (sum .eq. 0.0) then
                if (j .eq. 1)  go to 6285
                do 6187 k=1,nesum
 6187           ecurrt(k,j) = ecurrt(k,j-1)
              end if
          end if
c
c --- e)the vessel current values
c
          if (ivessel .eq. 1) then
            sum = 0.0
            do 6188 k=1,nvessel
 6188       sum = sum + ABS (vescurrt(k,j))
              if (sum .eq. 0.0) then
                if (j .eq. 1)  go to 6285
                do 6189 k=1,nvessel
 6189             vescurrt(k,j) = vescurrt(k,j-1)
              end if
          end if
c
c --- f)set up pcurmhd(the toroidal plasma  current  values).
c ---   these will be copied into totcur vector below if appropriate
c
          if (totcur(1) .lt. -0.9e+30) then
            if (pcurmhd(j) .eq. 0.0) then
              if (j .eq. 1)  go to 6285
                pcurmhd(j) = pcurmhd(j-1)
            end if
          end if
c
c --- g)the toroidal bfield
c
          if (btor .lt. -0.9e+30) then
            if (btormhd(j) .eq. 0.0) then
              if (j .eq. 1)  go to 6285
                btormhd(j) = btormhd(j-1)
            end if
          end if
      end do
c
      go to 6305
c
 6285 ierr = 1
      write  (nout, 6315)
      write  (ncrt, 6315)
 6315 format (' ERROR in specification of MHD boundary conditions')
c
 6305 continue
      end if          ! mhdmode = 'coils' branch
c
c --- check the MSE data. since we interpolate the values for each
c --- MSE channel in time we must guard against channels that may
c --- not be valid during some times. if such a condition is found
c --- then that channel is not used at all.
c
      if (use_stark) then
          do j=mxtbcmhd,1,-1
            istark = j
            if (timestark(j) .ne. 0.0)  go to 6333
          end do
          ierr = 1
          write  (nout, 6339) mxtbcmhd
          write  (ncrt, 6339) mxtbcmhd
 6339     format (' ERROR:  use_stark = .true. but timestark(i) = 0',
     .            ' for i = 1 to ', i5)
          istark = 1    ! so we can do rest of error checks
c
 6333     do j=1,istark-1
            if (timestark(j+1) .le. timestark(j)) then
              ierr = 1
              write  (nout, 6345)  j, timestark(j), timestark(j+1)
              write  (ncrt, 6345)  j, timestark(j), timestark(j+1)
 6345         format (' ERROR: timestark vector set incorrectly' /
     .                '     j, timestark(j), timestark(j+1) =',
     .                      i5, 2x, 2(1pe12.6))
            end if
          end do
          call expchk (tstark, nstark, mxtbcmhd, nstark, istark, xdum)
          do i=1,istark
              xdum(i) = 0.0
              do j=1,nstark
                  if (tstark(j,i) .eq. 0.0) then
                    fwtstark(j,i) = 0.0
                    sigstark(j,i) = 0.0
                  else
                    xdum(i)       = xdum(i) + 1.0
                  end if
              end do
          end do
          do i=1,istark-1
            if (ABS (xdum(i)-xdum(i+1)) .gt.  1.0e-5) then
              ierr = 1
              write  (nout, 6351)
              write  (ncrt, 6351)
 6351         format (' ERROR in MSE input data:' /
     .                '       number of channels varies with time')
            end if
          end do
          if (istark .eq. 1) then
            ierr = 1
            write (nout, 6352)
            write (ncrt, 6352)
 6352       format
     .      (' ERROR: must have at least two timepoints for MSE data')
          end if
          if (timestark(1) .gt. time0) then
            ierr = 1
            write  (nout, 6353)
            write  (ncrt, 6353)
 6353       format (' ERROR: timestark(1) must be .le. time0')
          end if
          if (timestark(istark) .lt. timmax .and. istark .gt. 1) then
            ierr = 1
            write  (nout, 6354)  istark, timmax
            write  (ncrt, 6354)  istark, timmax
 6354       format (' ERROR: timestark(', i5, ')',
     .              ' must be .ge. timmax =', 1pe12.4)
          end if
      end if
c
c --- set up the totcur vector if cursrce .eq. 'pcurmhd'
c --- if cursrce = 'eqdsk' we must obtain totcur(1) later
c

      if (cursrce .eq. 'pcurmhd') then
        call zeroa (totcur, nbctim)
c
c --- for nbctim values of time given in bctime,get totcur by interpolation
c --- from the itbcmhd values of pcurmhd:
c
      if (itbcmhd .eq. 1 .and. nbctim .gt. 1) then
        ierr = 1
        write  (nout, 6191)
        write  (ncrt, 6191)
 6191   format (' ERROR: totcur is to be taken from MHD section, but' /
     .          ' MHD section does not contain necessary information')
      else if (itbcmhd .eq. 1 .and. nbctim .eq. 1) then
          totcur(1) = pcurmhd(1)
****      totcur(1) = ABS (totcur(1))
      else
        do   j=1,nbctim
          do i=1,itbcmhd-1
            dtestl = bctime(j)-timeqbcd(i)
            dtestu = bctime(j)-timeqbcd(i+1)
            if (dtestl*dtestu .le. 0.0) then
              slope = (pcurmhd(i+1)-pcurmhd(i))/(timeqbcd(i+1)-
     .                                           timeqbcd(i))
              totcur(j) = pcurmhd(i)+slope*dtestl
****          totcur(j) = ABS (totcur(j))
            end if
          end do
        end do
        do j=2,nbctim
          totcur(j) = ABS (totcur(j))
        end do
      end if
c
      if (time0 .ne. bctime(1)) then
        totcur(1) = 0.0
c
c --- set totcur(1) to the value at time0
c
        do i=1,itbcmhd-1
           dtestl = time0-timeqbcd(i)
           dtestu = time0-timeqbcd(i+1)
           if (dtestl*dtestu .le. 0.0) then
             slope = (pcurmhd(i+1)-pcurmhd(i))/(timeqbcd(i+1)-
     .                                          timeqbcd(i  ))
             totcur(1) = pcurmhd(i)+slope*dtestl
****         totcur(1) = ABS (totcur(1))
           end if
         end do
         if (ltest_code .eq. 1) then
           do j=1,mxtbcmhd
             pcurmhd(j) = totcur(1)
           end do
         end if
         if (totcur(1) .eq. 0.0) then
           ierr = 1
           write  (nout, 6221)
           write  (ncrt, 6221)
 6221      format (' ERROR: totcur(1) not set; check bctime, pcurmhd',
     .             ' time0, timeqbcd')
         end if
        end if
        if (nbctim .gt. 1) then
          nbctims = nbctim
          do j=2,nbctims
            if (totcur(j) .eq. 0.0) then
              ierr   = 1    ! code will stop after all input is checked,
              nbctim = 1    ! so that we can check the rest of the input
              write  (nout, 6229) bctime(j)
              write  (ncrt, 6229) bctime(j)
 6229         format (' ERROR: totcur  = 0.0 at time ', f16.5)
            end if
          end do
        end if
      end if
      if (cursrce .ne. 'pcurmhd' .and. nbctim .gt. 1  .and.
     .                                 time .ne. bctime(1)) then
c
c --- get totcur(1) at t = time by interpolation from totcur vector
c
        do i=1,nbctim-1
          dtestl = time-bctime(i)
          dtestu = time-bctime(i+1)
          if (dtestl*dtestu .le. 0.0) then
            slope = (totcur(i+1)-totcur(i))/(bctime(i+1)-bctime(i))
            totcur(1) = totcur(i)+slope*dtestl
          end if
        end do
      end if
c
      if (cursrce .eq. 'inone') then
          ierrc = 0
          do j=1,nbctim
             if (totcur(j) .le. 0.0)  ierrc = 1
          end do
          if (ierrc .gt. 0) then
            write  (nout, 6240)
            write  (ncrt, 6240)
 6240       format (' ERROR: totcur is not set for all times')
            ierr = ierrc
          end if
      end if
c
c --- next get the toroidal b field if necessary. We take the average
c --- of the MHD values since btor is assumed time-independent:
c
      if (btorsrce .eq. 'btormhd') then
         sum = 0.0
         do 6193 j=1,itbcmhd
 6193      sum = sum+btormhd(j)
         btor = 1.0e4*sum/itbcmhd     ! in gauss
         if (ltest_code .eq. 1) then
           do j=1,mxtbcmhd
             btormhd(j) = btor / 1.0e4
           end do
         end if
      end if
      end if   ! end codeid = 'dee' case
c
c ----------------------------------------------------------------------
c set miscellaneous parameters
c ----------------------------------------------------------------------
c
      deltar       = rminor/(nj-1)
      rmajorvec(1) = rmajor
      do j=2,nj
        rmajorvec(j) = rmajorvec(j-1) + deltar   ! stored in file rhog.i
      end do
      time   = time0
      flim   = rmajor*btor
      btgeom = btor
      rgeom  = rmajor
      xmagn1 = rmajor
      ymagn1 = zax
c
c     including the elongation kappa
c
      kappa = elong(1)
      if (codeid .eq. 'onedee') then
        call elongt(time,kappa,dkapdt)
        do j=1,kj
          grho1_mesh(j)=1.0
          grho2_mesh(j) = 1.0
        enddo
      end if
c
c  set minimum on itmax to include particle conservation iterations
c          minimum reduced to 3 from 10--smw 1/6/89
c
c     for some non lin solve problems itmax =1 is appropriate
      if (itmax .lt. 3 .and. diffeq_methd .eq. 0)  itmax = 3
      
      inrad=0                                     !hsj 8/10/98
      if (nqrad .gt. 0)  inrad = 1
c
c  ensure jsourc = 1 for Yoka analysis runs
c
      if (iyoka .ne. 0)  jsourc = 1
c
c ----------------------------------------------------------------------
c set tweak flags
c ----------------------------------------------------------------------
c
      itweak = 0
      if (fusnin .eq. 0.0 .and.  ticin .eq. 0.0 .and.
     .    voltin .eq. 0.0 .and.   qcin .eq. 0.0 .and.
     .                          voltav .eq. 0.0)  ttweak = 0.0
      twkfar = 1.0
      if (voltin .ne. 0.0 .or. qcin .ne. 0.0) twkfar = 100.0
      if (iffar .eq.  1)  twkfar = 100.0  ! evolves bp 100 times as fast
      if (iffar .eq. -1)  twkfar =   1.0             ! what is this for?
      if (iffar .le. -2)  twkfar = IABS (iffar)      ! for MHD/transport
c
c     above, set twkfar smaller than 100 because 100 changes bp too fast
c
c ----------------------------------------------------------------------
c compute mesh quantities
c ----------------------------------------------------------------------
c
      call meshgeninitial (imesh, reqsp, kappa, rminor, reff,
     .                     r, dr, widths, delr, dzeta_adaptive,
     .                     include_adaptive, nout, ncrt, kj,
     .                     nj, njs, nps, ierr)
c
c ----------------------------------------------------------------------
c calculate some mesh-related quantities.
c ----------------------------------------------------------------------
c

      call meshgen (r, roa, ra, dr, rrp, drr, rrm, nj, codeid, time,
     .               rhod_ifs,rhod_max_ifs)
c
c ----------------------------------------------------------------------
c calculate some geometric quantities for 1-d case
c ----------------------------------------------------------------------
c
      rmin   = rmajor - rminor
      rmax   = rmajor + rminor
      zmin   = zax - kappa*rminor
      zmax   = zax + kappa*rminor
      volume = 2.0 * pisq * rmajor * reff**2
c
c ----------------------------------------------------------------------
c determine atomic weights of primary ions
c if dt is one species other species cannot be d or t
c ----------------------------------------------------------------------
c
      fd_thermal =0.0
      do 3410 i=1,nprim
        atw(i)    = 0.0
        atomno(i) = 1.0
        if (namep(i) .eq. 'he')  atomno(i) = 2.0
        if (namep(i) .eq. 'h' )  atw(i) = 1.0
        IF (namep(i) .eq. 'd' ) then
           atw(i) = 2.0
           fd_thermal = 1.0
        ENDIF
        if (namep(i) .eq. 't' )  atw(i) = 3.0
      
        IF (namep(i) .eq. 'dt') then
           atw(i) = fd * 2.0 + (1.0-fd) * 3.0
           fd_thermal =fd
        ENDIF
c
c       note that above sets effective mass for thermal fluid
c       beam dt composition is forced to be the same
c
        if (namep(i) .eq. 'he')  atw(i) = 4.0
        if (  atw(i) .ne.  0.0)  go to 3410
        ierr = 1
        write  (nout, 8420)  i, namep(i)
        write  (ncrt, 8420)  i, namep(i)
 8420   format (/ ' namep(', i2, ') = ', a8, ' is unrecognizable')
 3410 continue
c
c ----------------------------------------------------------------------
c determine atomic weights and charge  of impurity ions
c ----------------------------------------------------------------------
c
      if (nimp .eq. 0)  go to 3430
      do 3420 i=1,nimp
      k      = nprim + i
      atw(k) = 0.0
      atomno(k) = 0.0
      if (namei(i) .eq. 'he') then 
                  atw(k) =   4.0
                  atomno(k) = 2.0
      elseif (namei(i) .eq. 'c' ) then 
                  atw(k) =  12.0
                  atomno(k) = 6.0
      elseif (namei(i) .eq. 'o' )  then
                   atw(k) =  16.0
                   atomno(k) = 8.0
      elseif (namei(i) .eq. 'si')then
                   atw(k) =  28.0
                   atomno(k) = 14.0
      elseif (namei(i) .eq. 'ar')then
                   atw(k) =  40.0
                   atomno(k) = 18.
      elseif (namei(i) .eq. 'cr')then
                    atw(k) =  52.0
                   atomno(k) = 24.
      elseif (namei(i) .eq. 'fe')then
                    atw(k) =  56.0
                   atomno(k) = 26.0
      elseif (namei(i) .eq. 'ni')then
                    atw(k) =  59.0
                   atomno(k) = 28.0
      elseif (namei(i) .eq. 'kr')then
                    atw(k) =  84.0
                   atomno(k) = 36.0
      elseif (namei(i) .eq. 'mo')then
                    atw(k) =  96.0
                   atomno(k) = 42.0
      elseif (namei(i) .eq. 'w' )then
                    atw(k) = 184.0
                   atomno(k) = 74.0
      endif
      if (  atw(k) .ne.  0.0)  go to 3420
      ierr = 1
      write  (nout, 8430)  i, namei(i)
      write  (ncrt, 8430)  i, namei(i)
 8430 format (/ ' namei(',i2,') = ', a8, ' is unrecognizable')
 3420 continue
c
c ----------------------------------------------------------------------
c set charge parameters for primary ions
c ----------------------------------------------------------------------
c
 3430 do j=1,nj
        do k=1,nprim
          z(j,k) = 1.0
          if (atw(k) .ge. 4.0) z(j,k) = 2.0
          zsq(j,k) = z(j,k)**2
          dzdte(j,k) = 0.0      ! assumed fully ionized so dzdte never..
c                               ..changes for primary ions
        end do
c
c ----------------------------------------------------------------------
c initialize charge parameter for impurity ions
c ----------------------------------------------------------------------
c
        if (nimp .ne. 0) then
          do i=1,nimp
            dzdtim(j,i) = 0.0
          end do
        end if
      end do
c
c ------------------------------------------------------------------------
c determine the ion index for the CER ion
c ------------------------------------------------------------------------
c
      j_cer = 0
      if (cer_ion .ne. ' ') then
        do i=1,nprim
          if (cer_ion .eq. namep(i))  j_cer = i
        end do
        if (nimp .ne. 0) then
          do i=1,nimp
            if (cer_ion .eq. namei(i))  j_cer = i + nprim
          end do
        end if
      end if





c
c ----------------------------------------------------------------------
c now that species info is known we can check some beam items:
c check if transp style beam data input is requested
c Results in fatal error if use_nubeam is set and data is not found:
c ----------------------------------------------------------------------
      if(hdepsmth .LT. 0)fidiff_on = .FALSE.    ! This overides all fast ion smooting
      beam_data%data_allocated = .FALSE.
      call get_beam_data(beam_data)             ! ufiles_12.f90
 
      !make  sure nubeam0_dt is compatible with dtmax, dtmaxcrit:
      !IF(dtmax  .gt. nubeam0_dt )nubeam0_dt = 1.5*dtmax
      !IF(timmax .gt. timecrit .and. nubeam0_dt .lt. dtmaxcrit)
      !.                               nubeam0_dt = 1.5*dtmaxcrit
      
      IF(use_P_Nfreya)CALL P_Nf_events

c
c ----------------------------------------------------------------------
c  calculate number of dependent variables
c  note:iangrot = 1 if angular momentum effects are     considered,
c               = 0                             are not considered!
c               any other value of iangrot will not work!
c  if angular momentum coupling to energy equations is turned off
c  (angrcple = 0) then the coupling to the 2d source terms is also turned off.
c ----------------------------------------------------------------------
c
      if (iangrot .eq. 0)  go to 3074
      if (iangrot .eq. 1)  go to 3074
      write  (nout, 8432) iangrot
      write  (ncrt, 8432) iangrot
 8432 format (' ERROR: iangrot = ', i5, ' must be zero or one.')
      ierr    = 1
      iangrot = 0
 3074 nion    = nprim + nimp
      nk      = nion + 3
      nk      = nk + iangrot
      if (angrcple .eq. 0.0) then
        do j=1,3
          angrm2d(j) = 0.0
        end do
      end if

      IF(.NOT. ALLOCATED(brems_nions))ALLOCATE(brems_nions(nj,nion))
c
c ----------------------------------------------------------------------
c calculate flag for primary ion particle transport
c ----------------------------------------------------------------------
c
      idtrp =0
      iten = 0
      do i=1,nprim
        iten = MAX0 (itenp(i), iten)
        IF(itenp(i) .NE.  0)idtrp = idtrp+1
      end do
      IF(single_density_simulation .AND. idtrp .GT. 1)THEN
         ierr =1
         write(ncrt,2234)
         write(nout,2234)
 2234    Format(2x,'ERROR only 1 ion species transport is ',/,
     .    2x,'allowed if single_density_simulation =.true.')
      ENDIF
c
c ----------------------------------------------------------------------
c prevent impurity ion particle transport
c ----------------------------------------------------------------------
c
      do i=1,nimp
        iteni(i) = 0
      end do
c 
c ----------------------------------------------------------------------
c set up itran, a vector of transport flags
c allow -1 settings for gcnmp:
c ----------------------------------------------------------------------
c
      do i=1,nprim
        itran_gcnmp(i) = itenp(i)
        itenp(i)       = ABS(itenp(i))
        itran(i)       = itenp(i)
      end do

      DO  i=1,nimp
         k                   = nprim + i
         itran_gcnmp(k)      = iteni(i)
         iteni(i)            = ABS(iteni(i))
         itran(k)            = iteni(i)
      ENDDO
      itran_gcnmp(nk-2-iangrot) = itte
      itran_gcnmp(nk-1-iangrot) = itti
      itran_gcnmp(nk-iangrot) = itxj
      if (iangrot .eq. 1)  itran_gcnmp(nk) = itangrot
      itte = ABS(itte)
      itti = ABS(itti)
      itxj = ABS(itxj)
      itangrot = ABS(itangrot) 
      itran(nk-2-iangrot) = itte
      itran(nk-1-iangrot) = itti
      itran(nk  -iangrot) = itxj
      if (iangrot .eq. 1)  itran(nk) = itangrot
      IF(diffeq_methd == 3)THEN
         if(iangrot == 0)
     .   CALL STOP(' diffeq_methd =3 requires iangrot=1',1)
         DO i = 1,nk
            IF(itran_gcnmp(i) .LT. 0)eq_split_gcnmp = 1
         ENDDO
      ENDIF




c
c ----------------------------------------------------------------------
c calculate number of dependent variables for which a
c transport solution is desired
c ----------------------------------------------------------------------
c
      nkt = 0
      do k=1,nk
        if (itran(k) .eq. 1) nkt = nkt+1
      end do
      do j=nk+1,kk
        itranflag(j) = '     '
      end do
      do j=1,nk
        itranflag(j) = '(anl)'
        if (itran(j) .eq. 1)
     .  itranflag(j) = '(sim)'
      end do
c If faraday's law is the only simulation being done make sure
c glf23 is not called. Arrays allocated above if include_glf =1
c in inone should be, but are not removed at this time.
      if(nkt .eq. 1 .and. itxj .eq. 1)include_glf =0

c
c ----------------------------------------------------------------------
c save the name of each dependent variable
c ----------------------------------------------------------------------
c

      do i=1,kprim+kimp+6
         profiles_bcondspl(i) = ' '
      end do
      do i=1,nprim
         nameu(i)             = namep(i)
         profiles_bcondspl(i) = namep(i)
      end do
      k = nprim
      if (nimp .eq. 0)  go to 4044
      do i=1,nimp
        k                    = nprim + i
        profiles_bcondspl(k) = namei(i)
        nameu(k)             = namei(i)
      end do
 4044 nameu(k+1) = 'te'
      nameu(k+2) = 'ti'
      nameu(k+3) = 'bp'
      profiles_bcondspl(kprim+kimp+1) = 'ene'
      profiles_bcondspl(kprim+kimp+2) = 'te'
      profiles_bcondspl(kprim+kimp+3) = 'ti'
      profiles_bcondspl(kprim+kimp+4) = 'zeff'
      profiles_bcondspl(kprim+kimp+5) = 'curden'
      profiles_bcondspl(kprim+kimp+6) = 'angrot'
      if (iangrot .eq. 1)  nameu(k+4) = 'wr'
c
c ----------------------------------------------------------------------
c --- adjust the wneo array to account for the new 5 by 5 dimension if
c --- angular rotation is neglected.
c ----------------------------------------------------------------------
c
      if (iangrot .eq. 0) then
        ij = 0
        do 4046 j=1,5
        do 4046 i=1,5
          ij       = ij + 1
 4046     xdum(ij) = wneo(i,j)
      wneo(5,1)    = 0.0
      ij           = 4
      do 4047 j=2,5
        do 4048 i=1,4
          ij = ij + 1
 4048   wneo(i,j) = xdum(ij)
 4047 wneo(5,j) = 0.0
      end if
c
c ----------------------------------------------------------------------
c set bootstrap current flag
c turn bootstrap off for test case run
c ----------------------------------------------------------------------
c
      iboot = 0
      do k=1,3
        if (ltest_code .eq. 1  )  wneo(4,k) = 0.0
        if ( wneo(4,k) .ne. 0.0)  iboot     = 1
      end do
c
c ----------------------------------------------------------------------
c set anomalous transport flags
c if mhdonly = 1, transport will not be done at all.
c hence do not call wtypt in this case:
c ----------------------------------------------------------------------
c     interpolate for time-dependent w1typt ....
c
      w1typt = w1typ(1)
      w2typt = w2typ(1)
      w3typt = w3typ(1)
      w4typt = w4typ(1)
       vtypt =  vtyp(1)
      wtyp1  = ABS (w1typ(1)) + ABS (w2typ(1))
     .       + ABS (w3typ(1)) + ABS (w4typ(1))
      wtyp2  = ABS (w1typ(2)) + ABS (w2typ(2))
     .       + ABS (w3typ(2)) + ABS (w4typ(2))
      if (nbctim .gt. 1 .and. wtyp2 .gt. 0.0 .and. mhdonly .eq. 0)
     .  call wtypt(time)
c
      wtyp = wtyp1 + ABS (w12typ) + ABS (w13typ) + ABS (vtyp(1))
      wsaw = w1saw + w2saw + w3saw   ! simple mixing model
      wisl = w1isl + w2isl + w3isl
      wmix = w1mix + ABS (w2mix) + ABS (w3mix) + w4mix + w5mix
c
c ----------------------------------------------------------------------
c check inenez logic
c ----------------------------------------------------------------------
c
      set_dzdt = .TRUE.
      if (icenez .eq.  1)  inenez =  0
      if (icenez .eq. -1)  inenez = -1
      if (inenez .eq.  0)  go to 1950
      if (inenez .eq. -99)  go to 1950 !JMP
      if (inenez .eq. -98)  go to 1950 !JMP      
      if ( nprim .gt.  2 .and.  ifus .ge. 0)  go to 1940 
      if (  nimp .ne.  1)  go to 1940
      if (atw(nprim+1) .lt. 4.0)  go to 1940
      do 1930 i=1,nprim
 1930 if (atw(i) .ge. 4.0)  go to 1940
      if (nprim .eq. 1)  go to 1950
      if (zfrac .gt. 0.0 .and. zfrac .lt. 1.0)  go to 1950
      if (ifus  .lt. 0)   go to 1950
 1940 write (nout, 7940)  inenez,nprim,nimp,zfrac,(atw(i),i=1,3)
      write (ncrt, 7940)  inenez,nprim,nimp,zfrac,(atw(i),i=1,3)
 7940 format (// ' ERROR with inenez = ',i3,':'                       /
     .   ' Use nprim = 1,2; nimp=1; primary ions must be hydrogenic;' /
     .   ' if nprim = 2, use 0.<zfrac<1.'                             /
     .   ' nprim = ',i2,' nimp=',i2,' zfrac=',f6.3,' atw=',3f5.1)
      ierr = 1
c
c ---------------------------------------------------------------------
c CHECK FOR SINGLE FLUID BEAM,MULTIPLE FLUID THERMAL MODE  HSJ 1/5/96
c ---------------------------------------------------------------------
 1950 if (ifus .lt. 0) then
         ierrf=0
         iddd=0
         ittt=0
         ihehe=0
         itrtr=0
         if (nameb .ne. 'dt')  ierrf = 1
         if ( nimp .ne.  1  )  ierrf = 1
         if (nprim .eq.  0 .or. nprim .gt. 3)  ierrf = 1
            do j=1,nprim
               if (itran(j) .eq. 1   ) itrtr=itrtr+1
               if (namep(j) .eq. 'd' ) iddd=iddd+1
               if (namep(j) .eq. 't ') ittt=ittt+1
               if (namep(j) .eq. 'he') ihehe=ihehe+1
            end do
          if (itrtr .ne. nprim-1)  ierrf = 1
          if ( iddd .ne. 1      )  ierr  = 1
          if ( ittt .ne. 1      )  ierr  = 1
          if (nprim .eq. 3 .and. ihehe .ne. 1)  ierrf = 1
          if (nprim .eq. 2 .and. ihehe .eq. 1)  ierrf = 1
            do j=1,nneu
               if (namen(j) .eq. 'd')  iddd = iddd + 1
               if (namen(j) .eq. 't')  ittt = ittt + 1
            end do
           if (nneu .eq. 2 .and. iddd .ne. 2)  ierrf = 1
           if (nneu .eq. 2 .and. ittt .ne. 2)  ierrf = 1
           if (nneu .eq. 1 .and. (iddd .ne. 2 .or. ittt .ne. 2)) ierrf=1
           if (ierrf .ne. 0) then
              write  (ncrt, 1500)
              write  (nout, 1500)
              write  (nqik, 1500)
 1500         format (' ERROR in input with IFUS =-1 option')
           end if
           ierr = MAX (ierr, ierrf)
      end if

c
c ----------------------------------------------------------------------
c  CHECK FOR "SNAPSHOT MODE":
c --- For timmax .le. time0 we want to just evaluate the source terms
c --- and various balance tables  and quit.
c --- this mode makes sense only if running in analysis mode since
c --- the time derivatives are not known a priori at time zero.
c ----------------------------------------------------------------------
c
      if (timmax .le. time0) then
        prtlst(1) = time0           ! make sure we get printout at time0
        timmax    = time0 - 0.0001  ! otherwise getpsh may cause trouble
        if (nkt .gt. 0) then
          ierr = 1
          write (ncrt, 1977)
          write (nout, 1977)
 1977     format (
     .    ' subroutine INIT reports:'                              /
     .    '   SNAPSHOT MODE SELECTED (because timmax .le. time0)'  /
     .    '   this mode requires that all dependent variables'     /
     .    '   be set to run in analysis mode. (in simulation mode' /
     .    '   we do not have an estimate of the time derivatives'  /
     .    '   at time zero so the balance tables become arbitrary)')
          do j=1,nprim
             if (itenp(j) .ne. 0)  write (ncrt, 1978) j
             if (itenp(j) .ne. 0)  write (nout, 1978) j
 1978        format ('  itenp(',i1,') must be set to 0')
          end do
          if (itte .ne. 0)  write (ncrt, 1979)
          if (itte .ne. 0)  write (nout, 1979)
 1979     format ('  itte must be set to 0')
          if (itti .ne. 0)  write (ncrt, 1981)
          if (itti .ne. 0)  write (nout, 1981)
 1981     format ('  itti must be set to 0')
          if (itxj .ne. 0)  write (ncrt, 1982)
          if (itxj .ne. 0)  write (nout, 1982)
 1982     format ('  itxj must be set to 0')
          if (iangrot .eq. 1 .and. itangrot .ne. 0)  write (ncrt, 1983)
          if (iangrot .eq. 1 .and. itangrot .ne. 0)  write (nout, 1983)
 1983     format ('  itangrot must be set to 0')
        end if
        if (nbctim .gt. 1) then
            ierr = 1
            write (nout, 1991) nbctim
            write (ncrt, 1991) nbctim
 1991       format (' ERROR in snapshot mode, nbctim =', i5 /
     .              ' nbctim must be set to 1 for this case')
        end if
      end if



      rmhdgrid(nw) = -1.0e30          ! used as a flag until reset below
      zmhdgrid(nh) = -1.0e30

c--------------------------------------------------------------------
c --- read statefile and  and set profiles (if initialize_from_statefile
c --- is true)
      CALL set_initial_profiles 

c--------------------------------------------------------------------



c  neutral densities are calculated later in neucg: 
 4365 do 4370 i=1,2
      do 4370 j=1,nj
        enn (j,i) = 0.0
        ennw(j,i) = 0.0
        ennv(j,i) = 0.0
 4370   tn  (j,i) = 0.0



      if (wmix .ne. 0.0) then
        write  (ncrt, 1962)
        write  (nout, 1962)
 1962   format (' WARNING: mixing of angular momemtum not implemented')
      end if








c
c ----------------------------------------------------------------------
c define the angular momentum diffusivity model if necessary
c ----------------------------------------------------------------------
c
      if (iwangrot .ne. 0 .and. itangrot .eq. 1) then
        if (iwangrot .lt. 3 .and. iwangrot .gt. -2) then
          do 1963 j=1,nj
c
c --- parabolic model
c
            if (iwangrot .eq. -1)  xkangrot(j) =
     .        (xmtmdifs(1)-xmtmdifs(2))*(1.0-roa(j)**xmtmdifs(3))
     .         **xmtmdifs(4)+xmtmdifs(2)
c
c --- constant model
c
            if (iwangrot .eq. 1)  xkangrot(j) = xmtmdifs(1)
c
c --- linear model
c
            if (iwangrot .eq. 2)
     .      xkangrot(j) = xmtmdifs(1)*roa(j)+xmtmdifs(2)
 1963     continue
        end if
        if ((iwangrot .ge. 3) .and. (iwangrot .le. ksplin)) then
c
c --- spline model:
c
          if (rmtmdifs(1) .ne. 0.0)  go to 1966
          do j=3,ksplin
            k = j
            if (rmtmdifs(j-1) .ge. rmtmdifs(j))  go to 1966
            if (rmtmdifs(j) .eq. 1.0)  go to 1968
          end do
          go to 1966
 1968     call intrp (0, 1, rmtmdifs, xmtmdifs, k, roa, xkangrot, nj)
          go to 1969
 1966     ierr = 1
          write (nout, 1971)
          write (ncrt, 1971)
 1971     format (' ERROR in specification of rmtmdifs array')
 1969     continue
        end if
      if (iwangrot .eq. -2) then
        totweght = xmtmdifs(nprim+3)
        do 1973 j=1,nprim
        wetai(j) = xmtmdifs(j)
 1973   totweght = totweght+xmtmdifs(j)
      if (totweght .eq. 0.0) then
        ierr = 1
        write  (nout, 1972)
        write  (ncrt, 1972)
 1972   format (' ERROR in specification of Mattor-Diamond model' /
     .          '       total weight = 0.0')
      end if
      wneo(5,5) = wneo(5,5)*xmtmdifs(nprim+3)
      wetaie    = xmtmdifs(nprim+4)
      fimpurty  = xmtmdifs(nprim+1)
      if (fimpurty .eq. 0.0)  fimpurty = 1.0
      etaioff   = xmtmdifs(nprim+2)
      if ( etaioff .eq. 0.0)  etaioff = 1.0
      end if
        if ((iwangrot .lt. -3) .or. (iwangrot .gt. ksplin)) then
          ierr = 1
          write (nout, 1964) ksplin
          write (ncrt, 1964) ksplin
 1964     format (' ERROR in specification of iwangrot' /
     .            ' must be .ge. -2 and .lt. ', i5)
        end if
      end if

c
c ----------------------------------------------------------------------
c initial Kpol profile (m/sec/Tesla)
c ----------------------------------------------------------------------
c
 1952 if (ikpol .ne. 0) then
        m = 1
        if (rkpol(1,m) .ne. 0.0) then
          write  (nout, 8948)  m, rkpol(1,m)
          write  (ncrt, 8948)  m, rkpol(1,m)
 8948     format (' ERROR: rkpol(1,', i2, ') = ', e14.8, ' must be 0.0')
          ierr = 1
          rkpol(1,m) = 0.0
        end if
        njinkpol = 0
        do j=3,ksplin
          if (rkpol(j,m) .eq. 1.0)  njinkpol = j
        end do
        if (njinkpol .eq. 0) then
c
c --- rkpol(j,m) is not set correctly, try to copy rnormin
c
          if (jnormin .ge. 3) then
            call copya (rnormin, rkpol(1,m), jnormin)
            njinkpol = jnormin
          else
            write  (nout, 8949)  ksplin, m
            write  (ncrt, 8949)  ksplin, m
 8949       format (' ERROR: array rkpol(1..',i2,',',i2,') is not'  /
     .              '        set correctly. Last value at any time' /
     .              '        point must equal 1.0 exactly.')
            ierr     = 1
            njinkpol = 3
          end if
        end if
        iprofnbr = kprim + kimp + 7
        call intrp (-1, -1, rkpol(1,m), kpolin(1,m), njinkpol, roa,
     .               Kpol_exp, nj)
      end if
c
c --- end of initial Kpol profile
c
c ----------------------------------------------------------------------
c set up time-dependent profile parameters and boundary conditions,
c the latter in array bc
c ----------------------------------------------------------------------
c
      if (splninpt .eq. 'new') then
c
c --- set up and check time-dependent profiles:
c
      if (ikpol   .ne. 0)  call chksplnt (nbctim, kpolin, ksplin,
     .             kbctim, njinkpol, rkpol, nout, ncrt,
     .             ierr, knotskpol, 'Kpol')
      if (njcurb_external  .ne. 0)  call chksplnt(nbctim,
     .             curbeam_external,ksplin,kbctim,njcurb_external,
     .             rcurb_external,nout,ncrt,
     .             ierr,knotscurb_external,'trnspb')
c
      if (time .gt. bctime(1) .and. nbctim .GT. 1) then
c
        if (time .le. bctime(nbctim)) then
c
c --- get profiles at initial time:
c
            if (ikpol   .ne. 0)  call tsplinew (kpolin, rkpol, Kpol_exp,
     .                                knotskpol, kprim+kimp+7,
     .                                                       xdum, ydum)
            if (njcurb_external    .ne. 0)  
     .             call tsplinew (curbeam_external, rcurb_external, 
     .                   curb_external, knotscurb_external, 
     .                   kprim+kimp+6, xdum, ydum)
          else
            ierr = 1
            write  (nout, 1961)  time, bctime(nbctim)
            write  (ncrt, 1961)  time, bctime(nbctim)
 1961       format (' ERROR: time, bctime(nbctim) =', 2(2x, f12.5))
          end if
        else if (time .lt. bctime(1)) then
          ierr = 1
          write  (nout, 1959)  time, bctime(1)
          write  (ncrt, 1959)  time, bctime(1)
 1959     format (' ERROR: time, bctime(1) =', 2(2x, f12.5))
        end if
      end if








c
c     was rgc read in ? if so nrgc must equal nj
c
      nrgc=0
      do j=1,nj
        if (rgc(j) .gt. -1.0e30)  nrgc = nrgc + 1
      end do
c
c     calculate rgc if it is used (irgc .ne. 0)  and it is not
c     read in (nrgc .ne. nj)
      if (nrgc .ne. nj .and. irgc .ne. 0)
     .                            call time_grad_avg(nbctim,bctime,
     .                               angrotin,rangrot,kprim+kimp+6,
     .                               angrot,ksplin,kbctim,nj,nout,
     .                               knotsang,times_rgc,timee_rgc,
     .                               time0,timmax,roa,rgc)
c

c
c  -------------------------------------------------------------------
c  set boundary conditions 
c  --------------------------------------------------------------------
      do 4720 m=1,nbctim
         do k=1,nprim
            ucenter(m,k) = enpc(m,k)
            uedge(m,k) = enpb(m,k)
            ualp(m,k) = alpenp(m,k)
            ugam(m,k) = gamenp(m,k)
            bc(m,k) = enpb(m,k)
         end do
         if (nimp .eq. 0)  go to 4680
         do 4670 i=1,nimp
            k = nprim + i
            ucenter(m,k) = enic(m,i)
            uedge(m,k) = enib(m,i)
            ualp(m,k) = alpeni(m,i)
            ugam(m,k) = gameni(m,i)
            bc(m,k) = enib(m,i)
 4670    continue
 4680    k = nk - 2 - iangrot
         ucenter(m,k) = tec(m)
         uedge(m,k) = teb(m)
         ualp(m,k) = alpte(m)
         ugam(m,k) = gamte(m)
         bc(m,k) = teb(m)
         if (njte .gt. 0 .and.
     .        tein(1,2) .ne. 0.0 .and. splninpt .eq. 'old') then
            bc(m,k) = tein(njte,m)
         else if (njte .gt. 0   .and. splninpt .eq. 'new') then
            bc(m,k) = tein(knotste(m),m)
         end if
 4690    k = nk - 1 - iangrot
         ucenter(m,k) = tic(m)
         uedge(m,k) = tib(m)
         ualp(m,k) = alpti(m)
         ugam(m,k) = gamti(m)
         bc(m,k) = tib(m)
         if (njti .gt. 0 .and.
     .    tiin(1,2) .ne. 0.0 .and. splninpt .eq. 'old') then
            bc(m,k) = tiin(njti,m)
         else if (njti .gt. 0   .and. splninpt .eq. 'new') then
            bc(m,k) = tiin(knotsti(m),m)
         end if
 4700    k = nk-iangrot
         ucenter(m,k) = xjc(m)
         uedge(m,k) = xjb(m)
         ualp(m,k) = alpxj(m)
         ugam(m,k) = gamxj(m)
         if(u_vloop_bc)then
            call get_vloop_bc(bctime(m),bc(m,k)) ! returns bc(m,k) in volts
         else
            bc(m,k) = 0.2*totcur(m) ! gauss cm
            if (m .eq. 1)  bc(m,k) = 0.2*totcur1
         endif

         if (iangrot .eq. 0)  go to 4720
         k = nk
c
c --- no parabolic profiles for angular rotation:
c

c
         ucenter(m,k) = 0.0
         uedge  (m,k) = 0.0
         ualp   (m,k) = 0.0
         ugam   (m,k) = 0.0
         bc     (m,k) = angrotin(knotsang(m),m)
 4720 continue





c
c ----------------------------------------------------------------------
c integrate current density profile to obtain rbp profile, where
c rbp is r times the poloidal magnetic field.
c normalize curden (A/cm**2) and rbp (G-cm) to totcur (A).
c this calculation assumes fcap(j) = hcap(j) = 1 and so is correct
c only for 1-D runs.  for 1-1/2-D runs, rbp and curden are
c recalculated later.
C IF VLOOP BC IS USED THEN THESE CALCS ARE SKIPPED
C SO RBP == 0, CURDEN == 0
c ----------------------------------------------------------------------
c
      IF( .not. u_vloop_bc .AND. initialize_from_statefile 
     .                                        == .FALSE.)THEN
         do j=1,nj
            profin(j) = 2.0 * pi * r(j) * curden(j)
         end do
         call trap2 (r, profin, rbp, nj)
         fac = 0.0
         if (rbp(nj) .ne. 0.0) then
            fac = bc(1,nk-iangrot) / rbp(nj)
            if (ub(nk-iangrot) .eq. 0.0 .and. codeid .eq. 'onedee')
     .           ub(nk-iangrot) = bc(1,nk-iangrot)
            if (time0 .ne. bctime(1))  
     .                    call interp1 (time, bctime, nbctim,
     .                       bc(1,nk-iangrot), ub(nk-iangrot))
            fac = ub(nk-iangrot) / rbp(nj)
         end if
         if (fac .eq. 0.0)  go to 4770
         do 4760 j=1,nj
            curden(j) = 5.0 * fac * curden(j) ! amp/cm**2
 4760    rbp(j)    =       fac *    rbp(j) ! gauss cm
c

 4770    continue
      ENDIF




      if (splninpt .eq. 'new') then
c ---------------------------------------------------------------------
c  load  the edge zone for use as boundary conditions
c  the zone consists of all point j = fix_edge_te(m) to j = nj .
c  te_edge_zone
c  ti-edge_zone
c  ni-edge-zone
c  angrot_edge_zone
c  are set up to times in bctime(1..nbctime)
c  zone  values at arbitrary times  ( which may be > bctime(1) and must be <= bcime(nbctime)))
c  are determined by subroutine set_boundary_condition (the zone values on the radial
c  grid at the appopriate times are not
c  needed for initial problem setup like the profiles are so we dont dont
c  determine them here. 
c ---------------------------------------------------------------------
         j = kprim+kimp+2
       if(njte .ne. 0)
     .    call bc_zone(tein,knotste,
     .      rtein, bparte,j,ksplin, kbctim,bctime,
     .      nbctim,fix_edge_te,kj,nj,nout,
     .      profiles_bcondspl(j),roa,te_var_edge) 
         j = j+1
      if(njti .ne. 0)
     .    call bc_zone(tiin,knotsti,
     .      rtiin, bparti,j,ksplin, kbctim,bctime,
     .      nbctim,fix_edge_ti,kj,nj,nout, 
     .      profiles_bcondspl(j),roa,ti_var_edge)
         j = j+3
      if( iangrot .ne. 0)
     .    call bc_zone(angrotin,knotsang,
     .      rangrot, bparang,j,ksplin, kbctim,bctime,
     .      nbctim,fix_edge_rot,kj,nj,nout,
     .      profiles_bcondspl(j),roa,rot_var_edge)
cjmp.den start
c        j = kprim+kimp+1
c     if (njene .ne. 0)
c    .    call bc_zone(enein,knotsene,
c    .      renein, bparene,j,ksplin, kbctim,bctime,
c    .      nbctim,fix_edge_ni,kj,nj,nout,
c    .      profiles_bcondspl(j),roa,ni_var_edge) 
c      IF(fix_edge_ni(1) .GE. 1.0)THEN
       IF(fix_edge_ni(1) .GT. 2.0)THEN ! assumes grid pt number is input
         ni_index = fix_edge_ni(1)
      ELSE
         ni_index =0
         DO j=nj,1,-1
            if(roa(j) .le. fix_edge_ni(1)
     .      .AND.  ni_index == 0 )ni_index = j+1 ! assumes rhoa val is input
         ENDDO
      ENDIF

cjmp.den end
      !call stop("test",1)

      end if

c
c ----------------------------------------------------------------------
c  initialize solution vector
c  subroutine REDATE copies en,te,...etc into vector u
c  note that rbp may be modified so u must be updated later
c  to account for new rbp
c ----------------------------------------------------------------------
c
      id     = 0
      it     = 0
      idt    = 0
      ihe    = 0
      ihei   = 0    ! He can be primary or imourity
c
      do i=1,nprim
        if (namep(i) .eq. 'd' )  id  = i
        if (namep(i) .eq. 't' )  it  = i
        if (namep(i) .eq. 'dt')  idt = i
        if (namep(i) .eq. 'he')  ihe = i
      end do

 
      do i=1,nimp
         if(namei(i) .eq. 'he') ihei = i
      enddo
c

      call zen
      call redate (u, en, te, ti, rbp, nk, nj, kj, kk, iangrot, angrot)
      if (inenez .eq. 0)  go to 4830
      do k=1,nk-3-iangrot
        bc(1,k) = en(nj,k)  ! if inenez .eq. 1 zen has changed en
      end do                ! but if beam is on it will change again
c
c ----------------------------------------------------------------------
c initialize w1fact, vfact, w2fact, w3fact and w4fact,
c the profile factors multiplying w1typ, vtyp, w2typ, w3typ and w4typ
c ----------------------------------------------------------------------
c
 4830 do j=1,nj
        w1fact(j) = 1.0
         vfact(j) = 1.0
        w2fact(j) = 1.0
        w3fact(j) = 1.0
        w4fact(j) = 1.0
      end do
c
c ----------------------------------------------------------------------
c initialize quantities related to sawtooth oscillations
c ----------------------------------------------------------------------
c
      do j=1,nj
        esaw (j) = 0.0
        qmag (j) = 0.0
        qsawe(j) = 0.0
        qsawi(j) = 0.0
      end do
c
c ----------------------------------------------------------------------
c set up tweak2 initial values
c ----------------------------------------------------------------------
c
      if ( jsxr .ne. 1  )   s3mix = 0.0
      if ( jsxr .ne. 1  )  s71mix = 0.0
      if ( jsxr .ne. 2  )  s18mix = 0.0
      if (w2mix .ge. 0.0 .and.
     .    w3mix .ge. 0.0)   go to 4410
      epste = 0.3 * dtemix
      if (epste .eq. 0.0)   epste = 0.4 *  s3mix * te(1)
      if (epste .eq. 0.0)   epste = 0.6 * s71mix * te(1)
      if (epste .eq. 0.0)   epste = 0.6 * s18mix * te(1)
      epsti = 0.3 * dtimix
      if (epsti .eq. 0.0)   epsti = 0.1 * fusmix * ti(1)
      if (epsti .eq. 0.0)   epsti = trmix * epste
c
c ----------------------------------------------------------------------
c initialize sxr diode names
c do not call SXRCAL if only MHD is to be run (i.e., mhdonly=1)
c ----------------------------------------------------------------------
c
 4410 if (mhdonly .eq. 0)
     .call sxrcal (codeid,kappa,ene,-jsxr,nj,nw,nh,p,pmax,psir,r,rmajor,
     .        rin,rmax,te,rmhdgrid,zmhdgrid,zax,zmin,zmax,idiode,narray,
     .        namar,ndiode,roamin,sxr)
c
c ----------------------------------------------------------------------
c initialize CO2 interferometer array geometry
c ----------------------------------------------------------------------
c
      nco2     =   5
      rtco2(1) = 102.0
      rtco2(2) = 129.0
      rtco2(3) = 168.0
      rtco2(4) = 179.0
      rtco2(5) =   0.0
      call zeroa (denco2, nco2)
c
c ----------------------------------------------------------------------
c initialize z-effective diagnostic array geometry
c ----------------------------------------------------------------------
c
      nzeff      =   16
      zefrad     =    0.00163
      zeffwl     = 5230.0
      zeffdl     =  100.0
      rtzeff( 1) =  186.5
      rtzeff( 2) =  183.6
      rtzeff( 3) =  180.9
      rtzeff( 4) =  178.1
      rtzeff( 5) =  175.2
      rtzeff( 6) =  172.2
      rtzeff( 7) =  169.3
      rtzeff( 8) =  166.0
      rtzeff( 9) =  162.9
      rtzeff(10) =  159.8
      rtzeff(11) =  156.6
      rtzeff(12) =  153.5
      rtzeff(13) =  150.3
      rtzeff(14) =  146.8
      rtzeff(15) =  143.5
      rtzeff(16) =  139.6
      call zeroa (phzeff, nzeff)
c
c ----------------------------------------------------------------------
c check nterow and jterow
c ----------------------------------------------------------------------
c
      if (nterow .lt. 1 )  nterow =  0
      if (nterow .gt. 10)  nterow = 10
      if (nterow .eq. 0 )  go to 1996
      do i=1,nterow
        if (jterow(i) .lt. 1 )  jterow(i) = 1
        if (jterow(i) .gt. nj)  jterow(i) = nj
      end do

c
c ----------------------------------------------------------------------
c determine indices of the primary ion species corresponding to the
c neutral species and beam species - must not be species #3.
c ineut(i) indicates whether species i has neutrals present.
c When nneu = 1, 'in' has the index of the species with neutrals.
c ----------------------------------------------------------------------
c
 1996 if (nprim    .lt.       3 )  go to 4039
      if (namen(1) .eq. namep(3))  go to 4038
      if (namen(2) .eq. namep(3))  go to 4038
      if (nameb    .eq. namep(3))  go to 4038
      if (nampel   .eq. namep(3))  go to 4038
      go to 4039
c
 4038 write  (nout, 8038)
      write  (ncrt, 8038)
 8038 format (/ ' Species #3 must not have neutrals, beam or pellet.')
      ierr = 1
c
c --- if the beam is specified as deuterium set fdbeam to 1.
c
 4039 if (nameb .eq. 'd' )  fdbeam = 1.0
      if (nameb .eq. 't' )  fdbeam = 0.0                  ! 11/20/95 HSJ
      if (nameb .eq. 'dt' .and. ifus .ge. 0)  fdbeam = fd ! 11/20/95 HSJ
      in    = 1
      ibion = 1
      ipel  = 1
c
      do i=1,nprim
        ineut(i) = 0
        if (namen(1) .eq. namep(i))  in    = i
        if (nameb    .eq. namep(i))  ibion = i
        if (nampel   .eq. namep(i))  ipel  = i
      end do
      atw_beam=atw(ibion)
      if (nameb .eq. 'dt' .and. ifus .lt. 0) then
        ibion=-1
        atw_beam=fdbeam*2.0+(1.-fdbeam)*3.0
      end if
c
      if (btime(1) .gt. 0.0 .and. ifus .ge. 0) then
        if ((nprim .eq. 1  .and.  nameb .ne. namep(1)) .or.
     .      (nprim .eq. 2  .and. (nameb .ne. namep(1) .and.
     .                            nameb .ne. namep(2)      ))) then
          ierr = 1
          if (nprim .ge. 2) then
            write  (*   , 4037)
            write  (nout, 4037)
            write  (nqik, 4037)
 4037       format (' ERROR: beam must be primary ion #1 or #2')
          else
            write  (*   , 4035)
            write  (nout, 4035)
            write  (nqik, 4035)
 4035       format (' ERROR: beam must be same as primary ion species')
          end if
        end if
      else if (btime(1) .gt. 0.0 .and. ifus .lt. 0.0) then
c
c         beam must be dt
c
          if (nameb .ne. 'dt') then
             ierr=1
             write  (nout, 4041)
             write  (nqik, 4041)
             write  (*   , 4041)
 4041        format (' ERROR: ifus =-1 requires nameb = dt')
           end if
c
      end if
c
      if (nneu .eq. 1)  ineut(in) = 1
      if (nneu .eq. 2)  ineut(1 ) = 1
      if (nneu .eq. 2)  ineut(2 ) = 1
c
c ----------------------------------------------------------------------
c check on gas flux logic
c ----------------------------------------------------------------------
c

      do 1911 i=1,2
        if (ineut (i) .eq. 0)  go to 1911
        if (ipcons(i) .eq. 0)  go to 1911
        do m=1,nbctim
          if (gasflx(m,i) .ne. 0.0)  go to 1915
        end do
 1911 continue
      go to 1920
 1915 write  (nout, 7910)
      write  (nqik, 7910)
 7910 format (/ ' Note: injected gas flux is adjusted when ipcons=1.' /)
c
c ----------------------------------------------------------------------
c initialize some quantities for the neutral transport calculation
c ----------------------------------------------------------------------
c
 1920 ineu = 0
      if (codeid .eq. 'onedee' .and. raneut .eq. 0.0)
     .raneut = SQRT (kappa) * rminor
      reflec = 1.0
      flxmin = 1.0
      if (nengn .gt. 50)  nengn = 50
c
c ----------------------------------------------------------------------
c initialize some quantities for the neutral beam calculation
c ----------------------------------------------------------------------
c
      inub = 0
c
      do j=1,nj
        enbeam(j) = 0.0
        enbs  (j) = 0.0
        wbeam (j) = 0.0
        do jb=1,nbeams
          do ie=1,3
            enbsav(j,ie,jb) = 0.0
            wbsav (j,ie,jb) = 0.0
          end do
        end do
      end do
c
c     make sure the half and third energy components of the beam line
c     are turned off for negative ion sources
c
      do i=1,kb
        if (neg_ion_source(i) .gt. 0) then
          do j=2,ke
            fbcur(j,i) = 0.0
          end do
        end if
      end do
c
      if (iangrot .eq. 0)  ne_tk = 0

c
c ----------------------------------------------------------------------
c calculate some neutral beam parameters
c ----------------------------------------------------------------------
c
      mfm1   = mf - 1
      ibeam  = 0
      rmax   = rmajor + rminor
****  call ranset (ranseed)
      dummy  = RANDOM12 (ranseed)     ! initialize random number generator

D     call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr) !get processor id
D     call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr) !get num processors
C     instead of using one process to generate all random numbers,use
C     a separately initialized rand number generator for each process   HSJ
D      if(numprocs .gt. 1)then
D        do j =0, numprocs-1
D              ranseed = (j+1)*ranseed
D              if(myid .eq. j)then
D                  dummy = RANDOM12(ranseed)
D                  print *,'j,myid = ',j,myid 
D              endif
D        enddo
D      endif

      if (btime(1) .le. 0.0)  beamon(1) = timmax + 1.e6 * dtmax ! beam never on
      beam_end(1) = btime(1) + beamon(1)
      call my_vsrta (timbplt, ntimbplt)       ! make sure times are ascending
      kbm      = ntimbplt
      do j=1,ntimbplt
        if (ABS (timbplt(j) - beamon(1)) .lt. 1.0e-7)
     .           timbplt(j) = beamon(1) + 0.0001      ! required for chekdt
        if (timbplt(j) .lt. beamon(1))  kbm = kbm - 1 ! kbm gives # of valid times
      end do
      if (kbm .ne. ntimbplt) then
        do j=1,kbm
          timbplt(j) = timbplt(ntimbplt-kbm+j)
        end do
        do j=kbm+1,ntimbplt
          timbplt(j) = timmax + 2.0 * dtmax
        end do
      end if
      if (timbplt(1) .le. time0)  inubplt = 1
      if (    beamon(1) .le. time0)  ibeam   = 3



c
c ----------------------------------------------------------------------
c calculate some RF heating parameters
c ----------------------------------------------------------------------
c

      do k=1,krf
        irf(k)   = 0
        rfoff(k) = rfon(k) + rftime(k)
        if (rfon(k) .le. time0)  irf(k) = 3
        if(rfoff(k) .le. time0)  irf(k) = 0
      end do

c-----------------------------------------------------------------------
c -- read ech netcdf input if set
c     see ech_input gives a valid file name
c     (eq it does not include the string 'none':
C----------------------------------------------------------------------- 
      IF(ech_input .NE. 'none')THEN

             call echdat_read(ech_input)
             call set_irf
             ! this option creates a toray.in for each gyrotron
             ! at each time. So make sure no toray.in exists at the start: HSJ
             CALL delete_file('toray.in') ! non fatal delete

      ENDIF

c
      irfech = 0
      nrfrad = nprf(nrfzon+1)
      jrfmin = 1 
      jrfmax = nrfrad
c
c ----------------------------------------------------------------------
c check number of points in RF grid
c ----------------------------------------------------------------------
c
      active_rf_models =0
      fcd_path = .false. 
      icrm = 0 
      irferr = 0
      do k=1,krf
        indx_curray(k) = 0 
        irfmodel(k) = rfmode(k)
          if(irfmodel(k) .ne. no_rf)active_rf_models=
     .                                active_rf_models+1
          if(irfmodel(k) .eq. 'raytrace')then !set indecies into curray models
             icrm = icrm+1
             if(icrm .gt. kcrm)
     .         call STOP('no of curray models exceeds kcrm',1)
             indx_curray(k) = icrm
          endif

          if(irfmodel(k) .eq. 'fastcd' .and. .not. fcd_path)
     .    call get_fastcd(ncrt,nout,fcd_path)
          !add similar logic for toray and curray
          !that way we find out up front if things are set correctly

        rfmodel_power_e(k) =0.0 
        rfmodel_cd(k) =0.0       
        rfmodel_power_i(k) =0.0  
        if (irfmodel(k) .ne. 'ich')  go to 4085
        irferr = 0
        if (nrfrad .le. krfrad)  go to 1960
        write  (nout, 7960)  krfrad, nrfrad, imesh
        write  (ncrt, 7960)  krfrad, nrfrad, imesh
 7960   format (' too many RF grid points:'     /
     .          ' krfrad, nrfrad, imesh =', 3i6 /)
        ierr   = 1
        irferr = 1

c
c ----------------------------------------------------------------------
c       construct RF grid for ICH model:
c       if rfrad1ic and rfrad2ic are not specified, compute them here
c       for the 1-D case.  (the 1-1/2-D case is handled later.)
c       construct rfrad and rfrow.
c ----------------------------------------------------------------------
c

 1960   if (codeid .ne. 'onedee')  go to 4085
        if (irferr .ne.  0      )  go to 4085
        if (rfrad1ic .ne. 0.0 .and. rfrad2ic .ne. 0.0)  go to 4070
        dum = SQRT (reff**2-ylaunch**2)
        if (ylaunch .eq. 0.0)  dum = reff
        rfrad1ic = rmajor - dum
        rfrad2ic = rmajor + dum
c
 4070   do i=1,nrfzon
          delrad = 0.5 * (rfrad2ic-rfrad1ic)*(rfzone(i+1)-rfzone(i)) /
     .                   (nprf(i+1)-nprf(i))
          r0 = rfrad1ic + 0.5 * (rfrad2ic-rfrad1ic)*(1.0 + rfzone(i))
          do j=nprf(i),nprf(i+1)
            rfrad(j) = r0 + (j-nprf(i)) * delrad
            rfrow(j) = SQRT ((rfrad(j)-rmajor)**2+ylaunch**2)
          end do
        end do
c
        do j=1,nrfrad
          if (rfrow(jrfmin) .le. reff)  go to 4082
          jrfmin = jrfmin+1
        end do
c
 4082   do j=1,nrfrad
          if (rfrow(jrfmax) .le. reff)  go to 4085
          jrfmax = jrfmax-1
        end do
c
 4085   continue
      end do

c
c ----------------------------------------------------------------------
c set source profiles for 'wedge' RF model
c ----------------------------------------------------------------------
c
      call wedge_model (qine, qini, nj, r, irfmodel, -1)
c




c create merged print list (adds times in rfon and tdem mode):

      call assign_events


c ----------------------------------------------------------------------
c assign some values to event array
c this program checks for the following events
c     1. timmax
c     2. plot point (timplt) #D plto interval
c     3. print point (prtlst,timprt)
c     4. beam plot point (timbplt)
c     5. dteq
c     6. plot point (pltlst)
c     7. RF plot point (timrfp)
c     8. pellet injection time
c     9. beam on
c    10. beam of
c    11. rf(1) on
c    12. rf(1) off
c    13. rf(2) on
c    14. rf(2) off
c    15.
c      .
c      .
c      .
c    70.  rf(30) off  (krf = 30)
c    71.  time-dependent boundary conditions (if time_dep_beam =1)
c    72.  time dependent beam pulse on  (if time_dep_beam =1)
c    73.  time dependent beam pulse off (if time_dep_beam =1)
c    74.  beam integral upper limit
c    75.  beam integral lower limit
c    76.  nubeam switching  time
c    77.  blank
c    78   P_Nf beam(1) on  ! 32 Parallel nfreya beamlets
c    79   p_Nf beam(1) off
c    .
c    .
c   141   P_Nf beam (32) off 
c
c    values not set here are set in subroutine CHEKDT
c    see subroutine CHEKDT for details
c ----------------------------------------------------------------------
c
      nevents      = 10 + 2*krf + 1 ! old, with single pulse beam

      nevents      = nevents + 4  ! with multiple pulsed beam
      nevents      = nevents + 2*kb ! with P_Nfreya
      if(nevents .gt. kevents) 
     .call STOP ('subroutine INIT: kevents too small', 0)
      timevent( 1) = timmax
      timevent( 8) = timpel(1)
      timevent( 9) = beamon(1)
      timevent(10) = beam_end(1)  ! not used for P_Nfreya beams
      IF(use_P_Nfreya)THEN        
         timevent( 9) = HUGE(time0)
         timevent(10) = HUGE(time0)
         k= SIZE(beam_data%tbonac)
         ALLOCATE(beam_data%switch_on_off_delay
     .                   (index_p_nf_start:index_p_nf_stop))
         k =0

         DO j = index_p_nf_start,index_p_nf_stop,2
            k=k+1
            IF(k .GT. SIZE(beam_data%tbonac)) EXIT

            timevent(j)   = beam_data%tbonac(k)
            timevent(j+1) = beam_data%tboffac(k)
            beam_data%switch_on_off_delay(j)   =
     .               beam_data%tbonac(k) +10.*time_tol
            beam_data%switch_on_off_delay(j+1) = 
     .               beam_data%tboffac(k) +10.* time_tol
         ENDDO
      ELSE  ! P_Nfreya not used put time values beyond timmax:
         DO j = index_p_nf_start,index_p_nf_stop
            timevent( j) = HUGE(time0)
         ENDDO
      ENDIF


      do k=1,krf
        timevent( 9+2*k) = rfon (k)
        timevent(10+2*k) = rfoff(k)
      end do
      ipelet = 0
      npel   = 0
      if (timpel(1) .le. time0)  ipelet = 1
      index_bc = 10+2*krf + 1
      index_pbon = index_bc +1
      index_pboff = index_pbon + 1
      index_beam_ul = index_pboff + 1
      index_beam_ll = index_beam_ul + 1
c ---------------------------------------------------------------
c -- check for too close spacing in rf times
c -- this causes problems in events 
c -------------------------------------------------HSJ-3/20/12---
      ! check on time reosolution:
      DO k=1,krf
         if(rfon(k) .LT. 0.1*HUGE(time0))THEN
            DO j = k,krf
                 resol = ABS(rfon(k)-rfon(j))
                 IF(resol .lt. 2.*dt  
     .                    .AND. resol .GT. 10.*TINY(time0))THEN
                           dt = MIN(dt,resol)
                 ENDIF
            ENDDO
          ENDIF
       ENDDO
      ! check of time resolution:
      DO k=1,krf
         if(rfon(k) .LT. 0.1*HUGE(time0))THEN
            DO j = k,krf
                 resol = ABS(rfon(k)-rfoff(j))
                 IF(resol .lt. 2.*dt  
     .                     .AND. resol .GT. 10.*TINY(time0))THEN
                            dt = MIN(dt,resol)
                 ENDIF
            ENDDO
          ENDIF
       ENDDO
c       write(969,FMT='("dt,dtmin beforee = ",2(2x,1pe14.8))')dt,dtmin

       dtmin = MIN(dtmin,0.5*dt)
c       write(969,FMT='("dt,dtmin,after  = ",2(2x,1pe14.8))')dt,dtmin

c
c ----------------------------------------------------------------------
c set some fusion parameters
c ignflg: flag is set to 1 the first time we have ignition
c         before we have ignition ignflg = 0 , afterward -1
c ----------------------------------------------------------------------
c
 4110 ignflg = 0
c      id     = 0
c      it     = 0
c      idt    = 0
c      ihe    = 0
c
c      do i=1,nprim
c        if (namep(i) .eq. 'd' )  id  = i
c        if (namep(i) .eq. 't' )  it  = i
c        if (namep(i) .eq. 'dt')  idt = i
c        if (namep(i) .eq. 'he')  ihe = i
c      end do
c
c     check that dt and d or t do not coexist
c
      if (idt .ne. 0 .and. (id .ne. 0 .or. it .ne. 0)) then
        ierr = 1
        write  (*   , 4111)
        write  (nout, 4111)
        write  (nqik, 4111)
 4111   format (' ERROR: dt MIXTURE cannot be specified'     /
     .          ' in combination with individual d and/or t' /
     .          ' thermal ion species')
      end if
c
      if (ifus .eq. 0)  go to 4130     ! note that default is ifus = 1
      if (ifus .gt. 0) then
        ifus = 0
        if (idt .ne. 0                )  ifus = 1
****    if (id  .ne. 0 .and. it .ne. 0)  ifus = 2
        if (id  .ne. 0 .or.  it .ne. 0)  ifus = 2
      end if
c
 4130 iddfus  = 0
      idhe3fus = 0
      if (id  .ne. 0 .and. it .eq. 0)   iddfus   = 1
      if (idt .ne. 0                )   iddfus   = 2
      if (id  .ne. 0 .and. it .ne. 0)   iddfus   = 3
      if (id  .eq. 0 .and. it .ne. 0)   iddfus   = 4     ! for completeness
      IF(id .ne.   0 .AND. ihe .ne. 0)  idhe3fus = 1     ! prim d and he
      IF(id .ne.   0 .AND. ihei .ne. 0) idhe3fus = 2     ! prim d impurity he
      IF(id .ne.   0 .AND. ihe .ne. 0 
     .                .and. ihei .ne. 0)idhe3fus = 3     ! prim d,he and impurity he
      IF(idt .ne.  0 .AND. ihe  .ne. 0) idhe3fus = 4     ! prim dt and prim he
      IF(idt .ne.  0 .AND. ihei .ne. 0) idhe3fus = 5     ! prim dt and imp he
      IF(idt .ne.  0 .AND. ihe  .ne. 0 
     .               .AND. ihei .ne. 0)idhe3fus = 6      ! prim dt and he and imp he


      if (ifus .lt. 0 ) then
        iddfus  = 5
        adjzeff = 1
      end if
      do j=1,kj
        enalp (j) = 0.0
        enasav(j) = 0.0
        walp  (j) = 0.0
        wasav (j) = 0.0
      end do
c      print *,'iddfus,id,idt,it',iddfus,id,idt,it
c
c ----------------------------------------------------------------------
c  perform some calculations necessary for 1-1/2-d version of code
c ----------------------------------------------------------------------
c


      if (codeid .eq. 'onedee')  go to 5030
!      rmhdgrid(nw) = -1.0e30          ! used as a flag until reset below
!      zmhdgrid(nh) = -1.0e30
      if (irguess .lt. 0 .and. ifixshap .eq. 1)  go to 4980
c
c     read the Green's table (also returns mhdgrid in cm):
c
      if (mhdmode .eq. 'coils') then
        call readgren (ierrr)
        if (ierrr .eq. 0)  go to 4962
        write  (ncrt, 7000)  greentab
        write  (nout, 7000)  greentab
 7000   format (' FATAL ERROR: cannot read Green''s function table ',a8)
        ierr = 1
        go to 5030
c
 4962   ncoil2 = nfcoil / 2
        ncoilq = nfcoil * nfcoil
        xdim   = hotw
        ydim   = hoth
        redge  = hotrad    ! in cm
      end if
c
c force limiter to be read from eqdsk for "no coils" and "tdem" option
c
 4980 if (mhdmode .eq. 'no coils')  nlimiter = -1
c
c if mhdmode .eq. 'no coils' grid will be set in subroutine REQDSK
c
****  mhdmethd = 'cycred'
****  hotw = xdim
****  hoth = ydim
****  hotrad = redge        ! in cm
****  drmhd = hotw/(nw-1)
****  dzmhd = hoth/(nh-1)
****  do 4985 i=1,nw
*4985 rmhdgrid(i) = hotrad+(i-1)*drmhd
****  do 4990 j=1,nh
*4990 zmhdgrid(i) = -hoth*0.5+(j-1)*dzmhd
****  end if
c
c        eqdsk_tdem is stored in ename.mod to pass info about netCDF file
c        to ech:
c
         eqdsk_tdem = mhdmethd(1:LEN_TRIM(mhdmethd))
c
c --- if not defined above, mhdgrid will be set in subroutine REQDSK
c
c ----------------------------------------------------------------------
c --- limiter related calculations:
c --- get limiter points from eqdsk if necessary
c --- fill in additional limiter points if requested
c --- close the limiter (start and end on same (x,y) point )
c --- get min and max radial and vertical extension of limiter
c --- store these at the end of the limiter vectors
c --- set up the current scrape off array zero
c --- statefile(active if initialize+from statefile = true)
c --- does not have limiter points so we must still pick them
c --- up here
c ----------------------------------------------------------------------
c
      if (nlimiter .le. 0 .and. mhdmethd .ne. 'tdem' 
     .    .AND. .NOT. initialize_from_statefile )  call reqdsk
      if (nlimiter .le. 0 .and. mhdmethd .eq. 'tdem'
     .      .AND. .NOT. initialize_from_statefile )call get_cdf_data

c-----------------------------------------------------------------------
c     make sure xdim,ydim,redge set in inone are overwritten with
c     values consistent with statefile:
      IF(initialize_from_statefile)THEN
            xdim  = rmhdgrid(nw)-rmhdgrid(1)
            ydim  = zmhdgrid(nh)-zmhdgrid(1)
            redge = rmhdgrid(1)
      ENDIF
c------------------------------------------------------------------------
c
c --- if mhdgrid is still not set try to set it now, in cm:
c
      if (rmhdgrid(nw) .lt. -0.9e+30) then
        if (xdim .eq. 0.0) then
          write  (nout, 7412)
          write  (ncrt, 7412)
 7412     format (' ERROR: xdim, ydim or redge is not set correctly')
          ierr = 1
        end if
        if ( xdim .gt. 10.0)   xdim = xdim *0.01   ! convert to meters
        if (redge .gt. 10.0)  redge = redge*0.01
        drmhdgrd = xdim / (nw-1)
        do 7413 j=1,nw
 7413   rmhdgrid(j) = 100.0 * (redge+(j-1)*drmhdgrd)
      end if
      if (zmhdgrid(nh) .lt. -0.9e+30) then
        if (ydim .eq. 0.0) then
          write (nout, 7412)
          write (ncrt, 7412)
          ierr = 1
        end if
        if (ydim .gt. 50.0)  ydim = ydim * 0.01
        dzmhdgrd = ydim/(nh-1)
        do 7414 j=1,nh
 7414   zmhdgrid(j) = 100.0 * (-0.5*ydim+(j-1)*dzmhdgrd)
      end if
c
c --- if limiter points were not found and are not necessary break out
c
      if (nlimiter .le. 0) then
        if (ifixshap .eq. 1 .and. irguess .lt. 0)  go to 4997
        ierr = 1
        write  (nout, 7105)  TRIM(eqdskin)
        write  (ncrt, 7105)  TRIM(eqdskin)
 7105   format (/ ' ERROR: limiter points were not found in input' /
     .                8x, 'eqdsk file "', a, '"'                   /)
        go to 4997
      end if
c
      j     = nlimiter
      del   = 0.02*xdim/nw
      ifill = 0
      if (ifill .eq. 0)  go to 4988
      call fillp(xlimiter,ylimiter,j,nlimiter,del,xdum,ydum,maxlimpt)
c
 4988 dtoclose = (xlimiter(1)-xlimiter(nlimiter))**2+(ylimiter(1)-
     .            ylimiter(nlimiter))**2
      if (dtoclose .lt. 1.0e-10)  go to 4989
      nlimiter = MIN0 (nlimiter+1, maxlimpt)
      ylimiter(nlimiter) = ylimiter(1)
      xlimiter(nlimiter) = xlimiter(1)
c
c --- usmnmx (IMSL routine) gets min and max of a vector:
c

 4989 call my_usmnmx (xlimiter,nlimiter,1,xlmin,xlmax)
      call my_usmnmx (ylimiter,nlimiter,1,ylmin,ylmax)



      if (nlimiter .lt. maxlimpt-2) then
        xlimiter(nlimiter+1) = xlmin
        xlimiter(nlimiter+2) = xlmax
        ylimiter(nlimiter+1) = ylmin
        ylimiter(nlimiter+2) = ylmax
      else
        ierr = 1
        write  (nout, 8507)  nlimiter, maxlimpt
 8507   format (' ERROR: too many limiter points to fit'     /
     .              8x, 'must have nlimiter+2 .le. maxlimpt' /
     .              8x, 'nlimiter =', i5, '  maxlimpt = ', i5)
      end if
      xlimpos(3) = xlmax
      xlimpos(1) = xlmin
      xlimpos(2) = 0.5 * (xlimpos(1) + xlimpos(3))
      do i=1,3
        xlimpos(i) = xlimpos(i) * 100.0
      end do
      if (limpos .eq. 'left' )  go to 4993
      if (limpos .eq. 'right')  go to 4993
      tolvol = 1.0e30
c
c --- convert mhdgrid to meters, so it can be used in the following
c
 4993 cconst = 0.01
      call multpl1 (zmhdgrid, nh, cconst)
      call multpl1 (rmhdgrid, nw, cconst)
c
c --- check if limiter is consistent with MHD grid
c
      arbox = (xlmax - xlmin) * (ylmax - ylmin)
      ierrr = 0
      if (rmhdgrid(1)  .lt. xlmin)  go to 6110
      ierrr = 1
 6110 if (rmhdgrid(nw) .gt. xlmax)  go to 6120
      ierrr = 1
 6120 if (zmhdgrid(1)  .lt. ylmin)  go to 6130
      ierrr = 1
 6130 if (zmhdgrid(nh) .gt. ylmax)  go to 6140
      ierrr = 1
 6140 if (ierrr .eq. 1) then
        ierr = 1
        write  (nout, 8510)  xlmin, rmhdgrid(1), xlmax, rmhdgrid(nw)
        write  (nout, 8500)  ylmin, zmhdgrid(1), ylmax, zmhdgrid(nh)
        write  (ncrt, 8510)  xlmin, rmhdgrid(1), xlmax, rmhdgrid(nw)
        write  (ncrt, 8500)  ylmin, zmhdgrid(1), ylmax, zmhdgrid(nh)
 8500   format (' ERROR: limiter must be within MHD grid'          /
     .          '        ylim min, zgridmin, ylim max, zgrid max:' /
     .                   4x, 4(2x, 1pe12.3))
 8510   format (' ERROR: limiter must be within MHD grid'          /
     .          '        xlim min, rgridmin, xlim max, rgrid max:' /
     .                   4x, 4(2x, 1pe12.3))
      end if
c
c     set limcheck flag for use in subroutine bound
c
      limcheck = 0
      if (kappa .ne. 1.0)  limcheck = 1
c
c     construct zero array to scrape off current outside the limiter
c

      if (rmhdgrid(nw) .ne. -1.0e30) then
        call zlim (zero, nw, nh, nlimiter, xlimiter, ylimiter,
     .             rmhdgrid, zmhdgrid, 2)
      else
        ierr = 1
        write  (nout, 5150)
 5150   format (' ERROR: the MHD grid is not defined')
      end if
 
c
c --- set mhdgrid back to cm
c
      cconst = 100.0
      call multpl1 (zmhdgrid, nh, cconst)
      call multpl1 (rmhdgrid, nw, cconst)




c
c ----------------------------------------------------------------------
c --- some RF-related initialization
c ----------------------------------------------------------------------
c
 4997 do k=1,krf
        if (irfmodel(k) .ne. 'ich')  go to 4996
      end do
      go to 5030
 4996 if (irferr .ne. 0)  go to 5030
      if (rfrad1ic .ne. 0.0 .and. rfrad2ic .ne. 0.0)  go to 5010
      ychord  = ylaunch    * 0.01
      xlimcen = xlimpos(2) * 0.01
      call chord (nlimiter, xlimiter, ylimiter, xlimcen, ychord,
     .            rfrad1ic, rfrad2ic)
c
 5010 do i=1,nrfzon
        delrad = 0.5 * (rfrad2ic-rfrad1ic)*(rfzone(i+1)-rfzone(i))
     .               / (nprf(i+1)-nprf(i))
        r0 = rfrad1ic + 0.5*(rfrad2ic-rfrad1ic)*(1.0 + rfzone(i))
        do j=nprf(i),nprf(i+1)
          rfrad(j) = r0 + (j-nprf(i))*delrad
          rfrow(j) = SQRT ((rfrad(j)-rmajor)**2+ylaunch**2)
        end do
      end do




c ----------------------------------------------------------------------
c ---  NUBEAM -related initialization
c     beam on off times are given by tbona(:) and tboffa(:)
c     nubeam will be called when nubeam_on is set 
c ----------------------------------------------------------------------
c

      if(use_nubeam)then   ! move this
         nevents = nevents + 1
         if(nevents .gt. kevents) 
     .        call STOP ('subroutine INIT: kevents too small 2', 0)
         time_dep_beam =0
         nubeam_evolve = 0        !will be set to 1 in tport as necessary
         nubeam_on = MINVAL(beam_data%tbona(1:beam_data%nbeam))
         if(nubeam_on .gt. timmax)then             !no beams
            use_nubeam = .false.
            go to 5030
!         else if(time0 .gt. nubeam_on -.010)then
         else if(time0 .gt. nubeam_on .and. nubeam_restart.eq. 0)then
            print *,'beam_data%tbona =',beam_data%tbona
            print *,'Error, time0  must be leq  than'
            print *,'First beam on time'
            print *,'time0,beam on time =',time0,nubeam_on
            call STOP('init time problem',1)
         endif
         If(timmax .le. time0)then
            print *,' Snapshot mode (selected by timmax .le. time0)'
            print *,' is not possible when use_nubeam = .TRUE.'
            call STOP('nubeam init problem',1)
         ENDIF
         iterate_beam = .false.
         nbi_init = .false.  
         nzones = beam_data%nzone_nb       
         beam_data%nrhix    = nimp     

         nerngfi  = 5
         nbbox    = 0
         nsbeam = nbeams          !transp uses both nsbeams and nsbeam ?
         mig      = 7
         mis = 3*mig
         xp_nj = nj
         xp_npsi = npsi            !nubeam fails if this is larger than
                                   !eqdsk dimension ??

         miz = 2
!        mibs = 7                  ! max total number of fast ion species.
                                   !including beam ions,fusion product ions, 
                                   !and RF minority ions. defined in nbi_dimensions
         ng =0                     !current number of gases


         xp_eqdsk_name = eqdskin(1:len_trim(eqdskin))
         nseed    =  ranseed

         mrstrt   = 0             !checkpoint restart ability turned off

         ebdmax=1010*MAXVAL(ebkev(1:nbeams))              ! r4 variable                 
               ! maximum energy in distribution 
               ! function for beam ion species
               ! the maximum for fusion product species is set by the fusion
               ! product source energy (e.g. 3.5e6 eV for DT alphas-- adjusted
               ! upward to 5.0e6 allow for possible energy diffusion).
               ! **energy in eV**

         if(isym == 1) nlsym2b = .true.



         !thermal primary ion species:
c         ngmax = nprim +1     !no. primary thermal  species (add 1 in
                              !case primary species is dt mixture, 
                              !see sub set-thermal_density
         nprim_tr = nprim







             ! aperture half-width (or radius if circular), cm
             ! value of zero means: no 2nd aperture




!neutrals toroidal rotation speed:
         !note that #neutrals in nubeam is assumed equal to # primary ions
         allocate (omega_neut(1:xp_nj,nprim),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("omega_neut,sub init",0,istat)



         allocate (erngfi(1:nerngfi),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("erngfi,sub init",0,istat)

         allocate (nbebox(1:nbbox),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("nbebox,sub init",0,istat)

         allocate (nbsbox(1:nbbox),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("nbsbox,sub init",0,istat)


         allocate (xzmini(1:nmini),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("xzmini,sub init",0,istat)

         allocate (amini(1:nmini),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("amini,sub init",0,istat)


!!HSJ     all references to mj are taken out here 
!!        (because JMP removed nbi_dimensions module so mj and related
!!         parameers are undefined) HSJ 12/7/2011

!!HSJ         allocate (rhmin(1:mj,2),STAT = istat)
!!HSJ        if(istat .ne. 0)
!!HSJ     .        call allocate_error("rhmin,sub init",0,istat)

!!HSJ         allocate (den0mn(1:mig,1:mj,1:mis),STAT = istat)
!!HSJ        if(istat .ne. 0)
!!HSJ     .        call allocate_error("denomn,sub init",0,istat)

!!HSJ         allocate (omg0mn(1:mig,1:mj,1:mis),STAT = istat)
!!HSJ        if(istat .ne. 0)
!!HSJ     .        call allocate_error("omg0nm,sub init",0,istat)

!!HSJ         allocate (en0mn(1:mig,1:mj,1:mis),STAT = istat)
!!HSJ        if(istat .ne. 0)
!!HSJ     .        call allocate_error("en0mn,sub init",0,istat)

!!HSJ         if(nlntmj)then
!!HSJ            allocate (tmjsm(1:mj),STAT = istat)
!!HSJ            if(istat .ne. 0)
!!HSJ    .        call allocate_error("tmjsm,sub init",0,istat)
!!HSJ         endif


!!HSJ         allocate (xzimpj(1:mj),STAT = istat)
!!HSJ        if(istat .ne. 0)
!!HSJ     .        call allocate_error("xzimpj,sub init",0,istat)

!!HSJ         allocate (rhi(1:mj,1:miz),STAT = istat)
!!HSJ        if(istat .ne. 0)
!!HSJ     .        call allocate_error("rhi,sub init",0,istat)

         allocate (aimpj(1:xp_nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("aimpj,sub init",0,istat)

!!HSJ         allocate (nlfbmfpp(1:mibs),STAT = istat)
!!HSJ        if(istat .ne. 0)
!!HSJ     .        call allocate_error("nlfbmfpp,sub init",0,istat)






!!HSJ         allocate (rhbs(1:mj,1:miz,1:mibs),STAT = istat)
!!HSJ        if(istat .ne. 0)
!!HSJ     .        call allocate_error("rhbs,sub init",0,istat)

         allocate (omegag(1:xp_nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("rhi,sub init",0,istat)


         allocate (phiprg(1:xp_nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("phiprg,sub init",0,istat)
         phiprg(:) = 0.0d0

         allocate (curt(1:xp_nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("curt,sub init",0,istat)

         allocate (vpoh(1:xp_nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("vpoh,sub init",0,istat)

!         allocate (rhob(1:mj,1:mig,1:miz),STAT = istat)
!        if(istat .ne. 0)
!     .        call allocate_error("rhob,sub init",0,istat)

!!HSJ          allocate (rhix(1:mj,1:miz,1:mimpt),STAT = istat)
!!HSJ         if(istat .ne. 0)
!!HSJ      .        call allocate_error("rhix,sub init",0,istat)

         allocate (fracmini(1:nmini),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("fracmini,sub init",0,istat)
 

         allocate (curb_nub(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("curb_nub,sub init",0,istat)
 

         allocate (curb_nubp(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("curb_nubp,sub init",0,istat)

         allocate (qbeame_nub(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("qbeame_nub,sub init",0,istat)

         allocate (qbeame_nubp(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("qbeamep_nub,sub init",0,istat)

         allocate (qbeami_nub(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("qbeami_nub,sub init",0,istat)

         allocate (qbeami_nubp(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("qbeami_nubp,sub init",0,istat)

        if(istat .ne. 0)
     .        call allocate_error("qbeami_nub,sub init",0,istat)


         allocate (qbth_nub(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("qbth_nub,sub init",0,istat)

         allocate (qbth_nubp(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("qbth_nubp,sub init",0,istat)

         allocate (storqueb_nub(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("storqueb_nub,sub init",0,istat)

         allocate (storqueb_nubp(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("storoqueb_nubp,sub init",0,istat)


         allocate (sprbeame_nub(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("sprbeame_nub,sub init",0,istat)

         allocate (sprbeame_nubp(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("sprbeame_nub,sub init",0,istat)

         allocate (sprbeami_nub(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("sprbeami_nubp,sub init",0,istat)

         allocate (sprbeami_nubp(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("sprbeami_nubp,sub init",0,istat)


         allocate (wbeam_nub(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("wbeam_nub,sub init",0,istat)

         allocate (wbeam_nubp(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("wbeam_nubp,sub init",0,istat)

         allocate (enbeam_nub(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("enbeam_nub,sub init",0,istat)

         allocate (enbeam_species(1:nj,mfi),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("enbeam_species,sub init",0,istat)

         allocate (enbeam_species_p(1:nj,mfi),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("enbeam_species_p,sub init",0,istat)

         allocate (enbeam_species_c(1:nj,mfi),STAT = istat)
         enbeam_species(:,:)   = 0.0
         enbeam_species_p(:,:) = 0.0
         enbeam_species_c(:,:) = 0.0
        if(istat .ne. 0)
     .        call allocate_error("enbeam_species_c,sub init",0,istat)

         allocate (enbeam_nubp(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("enbeam_nubp,sub init",0,istat)





         allocate (beam_beamdtn_nub(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("beam_beamdtn__nub,sub init",0,istat)
         allocate (beam_beamdtn_nubp(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("beam_beamdtn__nubp,sub init",0,istat)


         allocate (beam_beamddp_nub(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("beam_beamddp__nub,sub init",0,istat)
         allocate (beam_beamddp_nubp(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("beam_beamddp__nubp,sub init",0,istat)




         allocate (beam_beamddn_nub(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("beam_beamddn__nub,sub init",0,istat)
         allocate (beam_beamddn_nubp(1:nj),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("beam_beamddn__nubp,sub init",0,istat)


         allocate (beam_beamtt2n_nub(1:nj),STAT = istat)
        if(istat .ne. 0)
     .      call allocate_error("beam_beamtt2n__nub,sub init",0,istat)
         allocate (beam_beamtt2n_nubp(1:nj),STAT = istat)
        if(istat .ne. 0)
     .      call allocate_error("beam_beamtt2n__nubp,sub init",0,istat)


         allocate (beam_thermaltth_df_nub(1:nj),STAT = istat)
        if(istat .ne. 0)
     .    call allocate_error("beam_thermaltth_df_nub,sub init",0,istat)
         allocate (beam_thermaltth_df_nubp(1:nj),STAT = istat)
        if(istat .ne. 0)
     .    call allocate_error("beam_thermaltth_df_nubp,sub init",
     .                                                          0,istat)

         allocate (beam_thermalddp_nub(1:nj),STAT = istat)
        if(istat .ne. 0)
     .    call allocate_error("beam_thermalddp_nub,sub init",0,istat)
         allocate (beam_thermalddp_nubp(1:nj),STAT = istat)
        if(istat .ne. 0)
     .    call allocate_error("beam_thermalddp_nubp,sub init",0,istat)

         allocate (beam_thermalddn_nub(1:nj),STAT = istat)
        if(istat .ne. 0)
     .    call allocate_error("beam_thermalddn_nub,sub init",0,istat)
         allocate (beam_thermalddn_nubp(1:nj),STAT = istat)
        if(istat .ne. 0)
     .    call allocate_error("beam_thermalddn_nubp,sub init",0,istat)

         allocate (beam_thermaltt2n_nub(1:nj),STAT = istat)
        if(istat .ne. 0)
     .    call allocate_error("beam_thermaltt2n_nub,sub init",0,istat)
         allocate (beam_thermaltt2n_nubp(1:nj),STAT = istat)
        if(istat .ne. 0)
     .    call allocate_error("beam_thermaltt2n_nubp,sub init",0,istat)

         allocate (beam_thermaldth_tf_nub(1:nj),STAT = istat)
        if(istat .ne. 0)
     .    call allocate_error("beam_thermaldth_tf_nub,sub init",0,istat)
         allocate (beam_thermaldth_tf_nubp(1:nj),STAT = istat)
        if(istat .ne. 0)
     .    call allocate_error("beam_thermaldth_tf_nubp,sub init",
     .                                                          0,istat)

cjmp.den start
         allocate (sorbn0_nub(1:nj),STAT = istat)
        if(istat .ne. 0)
     .    call allocate_error("sorbn0_nub,sub init",0,istat)
         allocate (sorbn0_nubp(1:nj),STAT = istat)
        if(istat .ne. 0)
     .    call allocate_error("sorbn0_nubp,sub init",
     .                                                          0,istat)

         allocate (sorbh_nub(1:nj),STAT = istat)
        if(istat .ne. 0)
     .    call allocate_error("sorbh_nub,sub init",0,istat)
         allocate (sorbh_nubp(1:nj),STAT = istat)
        if(istat .ne. 0)
     .    call allocate_error("sorbh_nubp,sub init",
     .                                                          0,istat)
cjmp.den end

      endif


 5030 continue
      if(use_nubeam)beam_beam_fusion = 1 



c
c ----------------------------------------------------------------------
c STOP  if error flag is not zero
c ----------------------------------------------------------------------
c

      if (ierr .eq. 0)  return
      call STOP ('subroutine INIT: ONETWO initialization errors', 2)
c
      end











      subroutine wedge_model (qine, qini, nj, r, irfmodel, modelrf)
c
      USE param
      implicit none
c
      integer nj, k, j, modelrf
c
c      include 'param.i'
      include 'wedgerf.i'
c
      integer       iwrfe_error, iwrfi_error
      real*8        qine(kj,*), qini(kj,krf), r(kj), frac, qfrac,
     .              sume, sumi
      character*(*) irfmodel(*)
c
      if (modelrf .lt. 0) then
          do k=1,krf
              if (irfmodel(k) .eq. 'wedge') then
                call zeroa (qine(1,k), nj)
                call zeroa (qini(1,k), nj)
                sume = 0.0
                sumi = 0.0
                do j=1,nj
                  if (r(j) .lt. rfrad1(k))  go to 10
                  if (r(j) .gt. rfrad2(k))  go to 10
                  frac      = (r(j)-rfrad1(k))/(rfrad2(k)-rfrad1(k))
                  qfrac     = a1rf(k) + frac*(a2rf(k)-a1rf(k))
                  qine(j,k) = qfrac * wrfe(  k)
                  qini(j,k) = qfrac * wrfi(  k)
                  sume      = sume  + qine(j,k)
                  qini(j,k) = qfrac * wrfi(  k)
                  sumi      = sumi  + qine(j,k)
   10             continue
                end do
                irfmodel(k) = 'input'
c
c             if wedge model is active, check that we actually
c             got something into qine and/or qini
c
              iwrfe_error = 0
              iwrfi_error = 0
              if (wrfe(k) .gt. 0.0 .and. sume .le. 0.0)  iwrfe_error = 1
              if (wrfi(k) .gt. 0.0 .and. sumi .le. 0.0)  iwrfi_error = 1
              if (iwrfe_error .eq. 1 .or. iwrfi_error .eq. 1)
     .        call STOP ('subroutine WEDGE_MODEL: input error', 278)
              end if
          end do
      else
        k = modelrf
        call zeroa (qine(1,k), nj)
        call zeroa (qini(1,k), nj)
        do j=1,nj
          if (r(j) .lt. rfrad1(k))  go to 20
          if (r(j) .gt. rfrad2(k))  go to 20
           frac     = (r(j)-rfrad1(k))/(rfrad2(k)-rfrad1(k))
          qfrac     = a1rf(k) + frac*(a2rf(k)-a1rf(k))
          qine(j,k) = qfrac * wrfe(k)
          qini(j,k) = qfrac * wrfi(k)
   20     continue
        end do
      end if
      return
c
      end
