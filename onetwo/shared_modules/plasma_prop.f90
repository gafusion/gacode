  MODULE  Plasma_properties
  USE nrtype,                                 ONLY : DP,I4B
  USE Vector_class
  USE ions_gcnmp,                             ONLY : name_size

  IMPLICIT NONE

  SAVE
  INTEGER,PARAMETER                  :: equil_type_size = 24
  REAL(DP),ALLOCATABLE,DIMENSION(:,:,:)    :: dcoefp,dcoefp_neo
  REAL(DP),ALLOCATABLE,DIMENSION(:,:)      :: ategrads,dtegrads
  REAL(DP),ALLOCATABLE,DIMENSION(:)        :: xketot,xdchitot,xchietot,xchiitot, &
                                              xkangrot,chiwneo,xkitot,dte_frozen
  TYPE frequency
     TYPE(Vector) :: omega_pe
     TYPE(Vector) :: omega_ce
     TYPE(vector),DIMENSION(:) ,POINTER   ::   omega_pi,omega_ci,omega_lh,omega_uh 
                                             ! omega_pi(i)%data(1:nj),etc

  END TYPE frequency


  TYPE geometry
     INTEGER(I4B)    nplasbdry,nr_mhd,nz_mhd,nlimiter
     REAL(DP)                                                              &
                     rgeom,rmag,deltao,areao,pindento,volo,rmajor,         &
                     kappa,rminor,rma,zma,rsep,zsep,rplasmax,zplasmax,     &
                     rplasmin,zplasmin,btgeom,circum

     TYPE(Vector) :: psivolpnpsi 
     TYPE(Vector) :: pindentnpsi
     TYPE(Vector) :: cxareanpsi
     TYPE(Vector) :: triangnpsi_u
     TYPE(Vector) :: triangnpsi_l
     TYPE(Vector) :: rmajavnpsi 
     TYPE(Vector) :: rminavnpsi
     TYPE(Vector) :: elongxnpsi 
     TYPE(Vector) :: sfareanpsi
     TYPE(Vector) :: grho1npsi
     TYPE(Vector) :: grho2npsi


     TYPE(Vector) :: mag_bp_contr
     TYPE(Vector) :: psivolpnj
     TYPE(Vector) :: pindentnj
     TYPE(Vector) :: cxareanj
     TYPE(Vector) :: triangnj_u
     TYPE(Vector) :: triangnj_l
     TYPE(Vector) :: rmajavnj 
     TYPE(Vector) :: rminavnj
     TYPE(Vector) :: rmin_half_width
     TYPE(Vector) :: rmaj_geom
     TYPE(Vector) :: elongxnj 
     TYPE(Vector) :: sfareanj
     TYPE(Vector) :: grho1nj
     TYPE(Vector) :: grho2nj

     TYPE(Vector) :: rplasbdry
     TYPE(Vector) :: zplasbdry
     TYPE(Vector) :: rmhdgrid
     TYPE(Vector) :: zmhdgrid     
     TYPE(Vector) :: rlimiter
     TYPE(Vector) :: zlimiter
  END TYPE geometry

  TYPE kinetic
     REAL(DP)  te0,ti0,p_intg,ene_bar_inbd,ene_bar_outbd
     TYPE(vector),DIMENSION(:) ,POINTER   :: en,flux,zsq,z,flux_conv  !en(1)%data(1:nj),etc
     TYPE(Vector)  ::  te
     TYPE(Vector)  ::  ti
     TYPE(Vector)  ::  ene
     TYPE(vector)  ::  angrot
     TYPE(vector)  ::  etor,press,pressb,walp,wbeam
     TYPE(vector)  ::  tot_press,ion_press,elect_press
     TYPE(vector)  ::  zeff
     TYPE(vector)  ::  vpol_nclass
     TYPE(vector)  ::  vpol
     TYPE(vector)  ::  vpar_nclass
     TYPE(vector)  ::  vpar
     TYPE(vector)  ::  er_tot_nclass
     TYPE(vector)  ::  er
     TYPE(vector)  ::  fluxe ! electron particle flux 
     TYPE(vector)  ::  fluxi ! ion  
     ! the following are for forcebal
     ! these quantities may not be available in the statefile:
     TYPE(vector)  ::  omega_imp        
     TYPE(vector)  ::  diagnostic_imp
     TYPE(vector)  ::  tor_imp_rot_outbd
     TYPE(vector)  ::  pol_imp_rot_outbd
     TYPE(vector)  ::  pol_imp_rot_avg


  END TYPE kinetic
      

  TYPE magnetic
     REAL(DP)  btor ,ali,tot_cur,totohm_cur,totboot_cur,totpar_cur,        &
          totbeam_cur,totrf_cur,betap,beta,vloop,psiaxis,psibdry,          &
          betan_ped,rhon_betan_ped,betan_global,R0
     INTEGER(I4B) set_cap1 ,npsi,ni_spline,mhd_ctr,run_mhd
     TYPE(Vector) :: fcap,gcap,hcap,rcap,r2capi,r2cap,ravgnpsi,ravginpsi,rcapi
     TYPE(Vector) :: q_value,rbp,curden,curohm,curboot,curpar,psivalnpsi,fpsinpsi
     TYPE(Vector) :: pprim,ffprim,bp,bprmaj,btotrmaj,betan,pressnpsi,qpsinpsi
     TYPE(Vector) :: ffprimnpsi,pprimnpsi
     TYPE(Vector) :: rhon_spline,ni_knot
     REAL(DP),DIMENSION(:,:),  POINTER  ::  psi
     CHARACTER(len = equil_type_size)   ::  equil_type
     CHARACTER(len = equil_type_size),DIMENSION(2)  ::  equil_avail  ! 'fixbd_var_grid'  or 'fixbd_unf_grid'
     LOGICAL mhd_info_avail
  END TYPE magnetic

  TYPE diffusivity                
     TYPE(Vector) :: chiinv
     TYPE(Vector) :: chieinv
     TYPE(Vector) :: xkineo
     TYPE(Vector) :: xkeneo
     TYPE(Vector) :: ftrap,eta,xnuse,chie_paleo
     TYPE(Vector) :: fixed_elct_eng_d_zct,fixed_ion_eng_d_zct
     TYPE(Vector) :: fixed_rbp_d_zct,fixed_torrot_d_zct
     TYPE(vector),DIMENSION(:) ,POINTER   :: xnus        !xnus(1)%data(1:nj),etc
     REAL(DP),DIMENSION(:),POINTER :: xkangrot,chiwneo,xchiitot,           &
                                      xdchitot, xchietot,xketot,           &
                                      xkitot
     REAL(DP),DIMENSION(:,:,:),POINTER :: dcoef
     REAL(DP),DIMENSION(:,:,:),POINTER :: dcoef_sav
     REAL(DP),DIMENSION(:,:,:),POINTER :: dcoef_neo
     REAL(DP),DIMENSION(:,:),  POINTER :: dcoef_glf23
     REAL(DP),DIMENSION(:,:),  POINTER :: dcoef_anal

     TYPE(vector),DIMENSION(:) ,POINTER   :: vpinch_nclass
                                !***%vpinch_nclass(1:nion)%data(1:nj)
     TYPE(vector),DIMENSION(:) ,POINTER   :: mmm_gammaW20
                                !***%mmm_gammaW20(1:nion)%data(1:nj)
     TYPE(vector),DIMENSION(:) ,POINTER   :: mmm_omegaW20
                               
     TYPE(vector),DIMENSION(:) ,POINTER   :: mmm_vconv
                                
     TYPE(vector),DIMENSION(:) ,POINTER   :: mmm_vflux

     TYPE(vector),DIMENSION(:) ,POINTER   :: fixed_ion_den_d_zct
 
     

     TYPE(Vector) :: mmm_gammaDBM,mmm_omegaDBM,mmm_xdi,mmm_xti,   &
                     mmm_xte,mmm_xdz,mmm_xvt,mmm_xvp,mmm_xtiW20,  &
                     mmm_xteW20,mmm_xdiW20,mmm_xtiDBM,mmm_xteDBM, &
                     mmm_xdiDBM,mmm_xteETG


   END TYPE diffusivity

  TYPE fusion_prod
     TYPE(Vector) :: neutr_ddn_th                     ! neutron rate profiles
     TYPE(Vector) :: neutr_ddn_beam_beam
     TYPE(Vector) :: neutr_ddn_beam_thermal
     TYPE(Vector) :: neutr_ddn_knock
     TYPE(Vector) :: neutr_ddn_tot
     TYPE(Vector) :: alpha_dthe4_th                   ! fusion alpha rates
     TYPE(Vector) :: alpha_dthe4_beam
     TYPE(Vector) :: alpha_dthe4_beam_thermal
     TYPE(Vector) :: alpha_dthe4_knock
     TYPE(Vector) :: alpa_dthe4_tot
     TYPE(vector) :: p_he3d_th                        ! proton production rates
     TYPE(vector) :: p_ddt
     TYPE(vector) :: p_he3d_beam_thermal
     REAL(DP) total_neutr_ddn_th,total_neutr_ddn_beam_beam,                 & ! integrated values
              total_neutr_ddn_beam_thermal,total_neutr_ddn_knock,           &
              total_neutr_ddn,                                              &
              total_alpha_dthe4_th,total_alpha_dthe4_beam,                  &  
              total_alpha_dthe4_beam_th,total_alpha_dthe4_knock,            &
              total_alpha_dthe4_tot,total_p_ddt,total_p_he3d,               &
              total_p_he3d_th,total_p_he3d_beam_thermal
  END TYPE

  TYPE power 
     TYPE(Vector)                         :: dpedt 
     TYPE(Vector)                         :: dnedt
     TYPE(vector),DIMENSION(:) ,POINTER   :: dpidt        !dpidt(1)%data(1:nj),etc
  END TYPE power

  TYPE power_den
     TYPE(Vector) :: qconde
     TYPE(Vector) :: qcondi
     TYPE(Vector) :: qconve
     TYPE(Vector) :: qconvi
     TYPE(Vector) :: qbeame
     TYPE(Vector) :: qbeami
     TYPE(Vector) :: qdelt
     TYPE(Vector) :: qione
     TYPE(Vector) :: qioni
     TYPE(Vector) :: qcx
     TYPE(Vector) :: qe2d
     TYPE(Vector) :: qi2d
     TYPE(Vector) :: qfuse
     TYPE(Vector) :: qfusi
     TYPE(Vector) :: qbfuse
     TYPE(Vector) :: qbfusi
     TYPE(Vector) :: qmag
     TYPE(Vector) :: qsawe
     TYPE(Vector) :: qsawi
     TYPE(Vector) :: qrad
     TYPE(vector),DIMENSION(:) ,POINTER   :: brems_nions 
                                !power_den%brems_nions(1:nion)%data(1:nj)
     TYPE(Vector) :: qohm
     TYPE(Vector) :: qrfi
     TYPE(Vector) :: qrfe
     TYPE(Vector) :: qexch


     TYPE(Vector) :: omegale
     TYPE(Vector) :: qomegapi
     TYPE(Vector) :: qangce
     TYPE(Vector) :: sprcxre
     TYPE(Vector) :: spreimpe
     TYPE(Vector) :: sprcxree



   END TYPE power_den

  TYPE particle_source
     TYPE(Vector) :: stfuse
     TYPE(Vector) :: sbfuse
     TYPE(Vector) :: sfus
     TYPE(Vector) :: spellet	
     TYPE(vector),DIMENSION(:) ,POINTER   :: srecom  ! particle_source%srecom(1)%data(1:nj)		
 END TYPE particle_source

 TYPE pellet_params
    LOGICAL       :: inject
    CHARACTER(len = name_size) name
    REAL(DP)      :: np
    REAL(DP)      :: dep_width
    REAL(DP)      :: rmin
    REAL(DP)      :: freq
 END TYPE pellet_params

 TYPE neutral_beam
    ! stores items appearing in statfile
    ! otherwise beam related items are in
    ! module neutral_beams.f90
     INTEGER(I4B)                              :: nbeams
     INTEGER(I4B),DIMENSION(:,:),      POINTER :: nmbrz
     INTEGER(I4B),DIMENSION(:,:),      POINTER :: nsample_izpt
     REAL(DP)     pwf_tot_source_intg
     REAL(DP),    DIMENSION(:,:),      POINTER :: prompt_pwr_in_plasma 
     REAL(DP),    DIMENSION(:,:),      POINTER :: stdev_nmbrz
     REAL(DP),    DIMENSION(:,:),      POINTER :: fap
     REAL(DP),    DIMENSION(:,:),      POINTER :: fbcur
     REAL(DP),    DIMENSION(:,:),      POINTER :: stdev_fap
     REAL(DP),    DIMENSION(:,:),      POINTER :: forb
     REAL(DP),    DIMENSION(:,:),      POINTER :: stdev_forb
     REAL(DP),    DIMENSION(:,:),      POINTER :: fwall
     REAL(DP),    DIMENSION(:,:),      POINTER :: stdev_fwall
     REAL(DP),    DIMENSION(:,:),      POINTER :: ebeam
     REAL(DP),    DIMENSION(:,:),      POINTER :: pbeam
     REAL(DP),    DIMENSION(:,:),      POINTER :: vbeam
     REAL(DP),    DIMENSION(:,:),      POINTER :: bneut
     REAL(DP),    DIMENSION(:,:),      POINTER :: stdev_bneut
     REAL(DP),    DIMENSION(:,:),      POINTER :: bion
     REAL(DP),    DIMENSION(:,:),      POINTER :: stdev_bion
     REAL(DP),    DIMENSION(:,:),      POINTER :: fb00
     REAL(DP),    DIMENSION(:,:),      POINTER :: stdev_fb00
     REAL(DP),    DIMENSION(:,:),      POINTER :: fb01
     REAL(DP),    DIMENSION(:,:),      POINTER :: stdev_fb01
     REAL(DP),    DIMENSION(:,:),      POINTER :: fb10
     REAL(DP),    DIMENSION(:,:),      POINTER :: stdev_fb10
     REAL(DP),    DIMENSION(:,:),      POINTER :: fb11
     REAL(DP),    DIMENSION(:,:),      POINTER :: stdev_fb11
     REAL(DP),    DIMENSION(:,:),      POINTER :: fber
     REAL(DP),    DIMENSION(:,:),      POINTER :: stdev_fber
     REAL(DP),    DIMENSION(:,:),      POINTER :: wb00
     REAL(DP),    DIMENSION(:,:),      POINTER :: stdev_wb00
     REAL(DP),    DIMENSION(:,:),      POINTER :: wb01
     REAL(DP),    DIMENSION(:,:),      POINTER :: stdev_wb01
     REAL(DP),    DIMENSION(:,:),      POINTER :: wb10
     REAL(DP),    DIMENSION(:,:),      POINTER :: stdev_wb10
     REAL(DP),    DIMENSION(:,:),      POINTER :: wb11
     REAL(DP),    DIMENSION(:,:),      POINTER :: stdev_wb11

     REAL(DP),    DIMENSION(:,:,:),    POINTER :: ftrapfi
     REAL(DP),    DIMENSION(:,:,:),    POINTER :: stdev_ftrapfi
     REAL(DP),    DIMENSION(:,:,:),    POINTER :: ftrapfit

     REAL(DP),    DIMENSION(:,:,:),    POINTER :: hibrz
     REAL(DP),    DIMENSION(:,:,:),    POINTER :: stdev_hibrz
     REAL(DP),    DIMENSION(:,:,:),    POINTER :: hibr

     REAL(DP),    DIMENSION(:,:,:),    POINTER :: hdepz
     REAL(DP),    DIMENSION(:,:,:),    POINTER :: stdev_hdepz
     REAL(DP),    DIMENSION(:,:,:),    POINTER :: hdep

     REAL(DP),    DIMENSION(:,:,:),    POINTER :: angmpz
     REAL(DP),    DIMENSION(:,:,:),    POINTER :: stdev_angmpz
     REAL(DP),    DIMENSION(:,:,:),    POINTER :: angmpf

     REAL(DP),    DIMENSION(:,:,:),    POINTER :: zetaz
     REAL(DP),    DIMENSION(:,:,:),    POINTER :: stdev_zetaz
     REAL(DP),    DIMENSION(:,:,:),    POINTER :: zeta

     REAL(DP),    DIMENSION(:,:,:,:),  POINTER :: hicmz
     REAL(DP),    DIMENSION(:,:,:,:),  POINTER :: stdev_hicmz
     REAL(DP),    DIMENSION(:,:,:,:),  POINTER :: hicm

     REAL(DP),    DIMENSION(:,:,:),    POINTER :: qb
     REAL(DP),    DIMENSION(:,:,:),    POINTER :: sb
     REAL(DP),    DIMENSION(:,:,:),    POINTER :: prompt_sb
     REAL(DP),    DIMENSION(:,:,:),    POINTER :: sbpure
     REAL(DP),    DIMENSION(:,:,:),    POINTER :: spb
     REAL(DP),    DIMENSION(:,:,:),    POINTER :: spbr
     REAL(DP),    DIMENSION(:,:,:),    POINTER :: pb0
     REAL(DP),    DIMENSION(:),        POINTER :: rhog_beam
     INTEGER(I4B) nj_beam  ! size of rhog_beam  on which above arrays are defined
     REAL(DP) d_beam_spin_pol

 END TYPE neutral_beam

  TYPE (geometry) dischg ; TYPE(kinetic) profile ; TYPE(magnetic) :: mhd_dat
  TYPE (diffusivity) diffuse ; TYPE (power) wpdot ; TYPE (power_den) pwrden
  TYPE (particle_source) prtcl_src
  TYPE (pellet_params)   pellet 
  TYPE (fusion_prod)     fus_prod
  TYPE (neutral_beam)    neut_beam
  TYPE (frequency)       plasma_frequencies

! vectors are nullified in sub null_pointers !!!!!!!
  
  END MODULE  Plasma_properties


