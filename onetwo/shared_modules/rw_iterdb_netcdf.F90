SUBROUTINE iter_dbase_nc
  ! ----------------------------------------------------------------------
  ! --- collects and writes out (or reads in) data for ITER database
  ! --- netcdf version. The netcdf input has different units,names, etc.
  ! --- from the Onetwo variables. Hence we use the modules from  gcnmp, 
  ! --- read in and isolate these quantities and then move them
  ! --- appropriately to Onetwo variables (in set_onetwo_vars)
  ! --- This routine is considerably extended over the older iterdb files
  ! --- and consequently may not work with codes that expect the older versions.
  ! --- On the other hand if the netcdf implementation for reading the 
  ! --- iterdb file is done  correctly in these other codes  then the
  ! ---  additions presented here will 
  ! --- not cause problems since netcdf will ignore variables that are not
  ! --- defined. 
  ! -------------------------------------------------------HSJ-8/15/07------
  ! 
  ! --- INPUT/OUTPUT
  !

  !     irwflag          read/write flag, included to make the reading
  !                      of the file created by this subroutine
  !                      consistent with the writing of the file.
  !                      irwflag = 0  ==>  WRITE the data
  !                      irwflag = 1  ==>  READ  the data
  !
  !     iterdb_file_name   name of file which will be read 
  !                      (if called in read mode)
  !     iterdb_outpt     name of output file (created in this routine,
  !                      iterdb_outpt = iterdb_file_name (could add //time)
  !     iterdsc          switch used to output description of parameter
  !
  !
  !     nplasbdry
  !     rplasbdry(j)
  !     zplasbdry(j)     j=1,2..nplasbdry describes the plasma boundary
  !
  !     q(j)             j=1,2,...nj, the safety factor on rho grid
  !     betap
  !     beta
  !
  !     tocur            total plasma current, amps
  !

  !     rho(j)           j=1,2,..npsi the rho values corresponding
  !                      to the psi grid in meter
  !     qpsi(j)          the corresponding q values
  !     psivolpnpsi(j)       volume, meter^3, enclosed by flux surface
  !     grho1(j)         < ABS (grad rho)>
  !     grho2(j)         <(grad rho)**2>
  !     rmajavnpsi(j)    j=1,2,..npsi average major radius
  !                      at elevation of magnetic axis
  !     rminavnpsi(j)    same for minor radius
  !     psivalnpsi(j)        j=1,2..npsi
  !
  !     triangnpsi_u(j)    (upper)triangularity
  !     triangnpsi_l(j)  (lower)  triangularity
  !     pindentnpsi(j)   indentation
  !     sfareanpsi(j)    surface area of flux surfaces
  !     cxareanpsi(j)    cross-sectional area
  !
  !     btor             mag field at R0, in vacuum, tesla
  !     elong(j)         elongation
  !
  !
  !     fcap(j)
  !     gcap(j)
  !     hcap(j)          the geometric factors
  !     ali              plasma inductance
  !

  !     namep(i)         i=1,..nprim
  !     fd_thermal       if primary species is 'dt' then fd_thermal gives
  !                      fraction of d in 'dt' mixture. If primary species 
  !                      is/is not 'd'  then set fd_thermal 1/0
  !     namei(i)         i=1,..nimp
  !     namen(i)         i=1,..nneu 
  !     zeff(j)          j=1,2...nj
  !     nameb(i)         i = 1,nbion name of beam ions ('dt', 'd','t','h','he')
  !     nbion            number of fast beam ion  species. For a 'dt'
  !                      beam mixture use nbion =1. The fraction of
  !     fd_beam          d in the 'dt' mixture is fd_beam. For single species
  !                      beam fd_beam set fdbeam = 1.0
  !     sbeame            electron source from beam
  !     sbeam(j,i)       thermal ion source from beam, grid point j, species nameb(i)
  !
  !     enbeam(j,i)      j=1,2..nj,i=1..nbion  fast ion density due to beam,
  !                      grid point j, species nameb(i).
  !     wbeam(j)         beam stored energy density
  !     enalp(J)         alpha particle  density
  !     walp(j)          alpha particle stored energy density
  !                      (boundary condition for w_alpha, if used)
  !
  !     enn(j,k)         j=1,2..nj(njs),k=1,2 neutral density
  !     ennw             neutral density due to wall source
  !     ennv                                    volume
  !     volsn            volume source of neutrals
  ! 
  !     rmajor 
  !     kappa
  !     btor
  !
  !     r(j)             j=1,2...nj,rho grid
  !
  !     time            current time, sec
  !     tGCNMf        the transport equations are evolved from time
  !                     to tGCNMf (sec).
  !                     The interval (time,tGCNMf] will be broken into
  !                     smaller intervals by the solver as appropriate.
  !                     results will be returned in the form of a new
  !                     iterdb file at tGCNMf. It is  assumed that
  !                     all particle and  energy sources are constant
  !                     in this time interval. Metrics introduced by
  !                     flux surface averaging (eq *cap quantites)
  !                     are also assumed constant in time over 
  !                     this interval
  ! 
  !     time_bc         boundary condition time. This time must be 
  !                     >= tGCNMf!!
  !                     It is assumed that the
  !                     interval (time,time_bc) is small enough that
  !                     linear inerpolation of boundary conditions is
  !                     reasonble. 
  !                     Initial conditions required
  !                     at t = time are taken from en(j,k),ene,zeff,te,ti,
  !                     curden,angrot (see below)  as appropriate. 
  !                     Boundary conditions
  !                     required to evolve from time to time_bc are
  !                     obtained by linear interpolation  between these 
  !                     profiles and the values given in en_bc(j,k)
  !                     ene_bc,zeff_bc,te_bc,ti_bc,angrot_bc
  !                     (see below).  Although only
  !                     information for the plasma edge is required for
  !                     the boundary conditions it is assumed that
  !                     the *_bc profiles are given over the entire rho grid
  !                     This allows us to apply boundary conditions 
  !                     for each dependent variable at a different rho
  !                     value, rather than limiting the boundary to just 
  !                     rho= 1.0. The switches  
  !                         fix-edge_te_bc  (= any value in (0.0,1.0] )     
  !                         fix_edge_ti_bc  (= any value in (0.0,1.0] )
  !                         fix_edge_rot_bc (= any value in (0.0,1.0] )
  !                         fix_edge_ni_bc  (= any value in (0.0,1.0] )
  !                         interact with the *_bc profiles in that
  !                         for rho > fix_edge_** the value of the *_bc
  !                         profile is used. To be specifie consider
  !                         for example the te profile: 
  !                         Suppose  fix_eddge_te =0.91 . Lets assume 
  !                         that the closest
  !                         normalized rho grid point is at rho =0.925.
  !                         then the boundary condtion for te at t= time is
  !                         takens as te(rho=0.925) . At t= time_bc the
  !                         bc is taken as te_bc(rho = 0.925). For arbitrary
  !                         times in the interval [time,time_bc] the bc for
  !                         te is obtained by lineat interpolation in time 
  !                         at the grid point rho = 0.925. (Because the rho grid
  !                         may change in time normalized rho is used.)
  !                         For all rho grid points > 0.925 the te profile
  !                         is also obtained by this linear inerpolation process
  !                         But note that the temperature at these points
  !                         is given a priori by the profiles te and te_bc.
  !                         The te profile for rho < 0.925 is determined 
  !                         by solving the PDE's. Note that we always solve
  !                         Faraday's law out to the plasma edge( rho =1.0) and
  !                         rhis requires te, etc all the way out to the edge.
  !                         Hence zpecification of te_bc in the edge zone will 
  !                         still influence the solution. Finally if te is run
  !                         in analysis mode then the linear inerpolation in 
  !                         time discussed above is done for all rho values from
  !                         0.0 to 1.0  inclusive. (The above holds for all
  !                         profiles, not just te).
  !                      



  !                      
  !     eqtime          time of last equilibirum 
  !     nj              transport grid size
  !     nprim
  !     nneu
  !     nimp
  !     nion
  !
  !     psir_grid(j)         j=1,2...nj psi on r (i.e., rho) grid
  !
  !     profiles at t = time (used as inital conditions)
  !     te(j)           electron and ion temperatures,keV
  !     ti(j)
  !     ene(j)          electron density
  !     en(j,k)         j=1,2..nj,k=1..nion,where nion=nprim+nimp is the
  !                     number of primary plus impurity ion species
  !     rbp(j)          rho*bp field
  !     profiles at t=time_bc (used to determine boundary condtions)
  !         te_bc
  !         ti_bc
  !         ene_bc
  !         en(j,k)_bc
  !         zeff_bc!
  !         angrot_bc

  !     
  !     fixe_edge_**  switches allow  solving the problem from
  !     rho =0 to rho = fixedge_**  .le. 1.0
  !     for fix_edge < rho <=1.0 the code takes the input profile values.  
  !     these values are constant during a call to gcnm but could
  !     change from one call to the next,as could rho,etcc so this
  !     is a fully dynamicallay controlled model.     
  !     fix-edge_te_bc      	    
  !     fix_edge_ti_bc           
  !     fix_edge_rot_bc 
  !     fix_edge_ni_bc        





  !     currf(j)
  !     curbe(i)
  !     curbi(j)
  !     qrad(j)         radiated power
  !     totohm
  !     totbeam
  !     totboot
  !     totrf           total currents, amps
  !     qconde,qcondi
  !     qconve,qconvi
  !     qdelt
  !     qexch    not currently used
  !     qione,qioni
  !     qbeame,qbeami
  !     qrfe, qrfi
  !     qe2d, qi2d
  !     qtfuse,qtfusi
  !     qmag
  !     qsawe,qsawi
  !     qbfusi
  !     sbeam
  !     dnedt(j)         change in electron density at this time 1/(M**3sec)
  !     dangrotdt(J)     change in toroidal rotation
  !     cheinv(j)
  !     chiinv(j)        electron and ion thermal cond.
  !     xkineo(j)        ion neoclassical thermal conductivity
  !     xkeneo(j)        electron neoclassical thermal conductivity
  !
  !     angrot(j)       angular rotation speed, radian/second
  !     storqueb(j)     beam torque density nt-m/m**3 
  !
  !     eps(j)          horizontal inverse aspect ratio = (rmax-rmin)/(rmax+rmin)
  !                     on each flux surface. Add to iterdb file
  !     xhm2
  !     xi11
  !     xi33
  !     xips 
  ! -------------------------------------------------------- HSJ
  !
  USE typeSizes                  !part of netcdf 
  USE netcdf                     !ditto

  !
  USE nrtype,                    ONLY : DP,I4B,I2B


  USE Plasma_properties ,        ONLY : dischg,profile,mhd_dat,                &
                                        diffuse,pwrden, wpdot,prtcl_src,       &
                                        pellet,fus_prod,neut_beam,             &
                                        plasma_frequencies

  USE shot_info,                 ONLY : shot_id
  USE vector_class,              ONLY : new_Vector,get_element,                &
                                        real_mult_Vector,list,zero_Vector,     &
                                        length_vector,get_values,              &
                                        delete_Vector_nf,delete_Vector

  USE grid_class,                ONLY : nj,psir_grid,rho_grid,rho_gridn,       &
                                        eps,xhm2,xi11,xi33,xips,xhm20,xi110,   &
                                        xi330,xips0,rho_mhd_gridnpsi

  USE error_handler,             ONLY : lerrno,terminate,iomaxerr

  USE nf_param,                  ONLY : kcm

  USE pl_freq,                   ONLY : allocate_plasma_freq,id_plas_freq,      &
                                        id_ci_freq,id_lh_freq,id_uh_freq

#ifdef GCNMP
#define USESUBS
#endif
#ifdef NFREYA
#define USESUBS
  USE nfreya_version,            ONLY : p_nfreya_ver
#endif

#ifdef USESUBS
  USE io_gcnmp,                  ONLY : ncrt,nlog
  USE gcnmp_version,             ONLY : gcnmp_ver
#else
  USE io,                        ONLY : ncrt,versid,nlog => nitre
#endif


  USE source_terms_gcnmp,        ONLY : stsource,scx,sion,srecom,sbcx,sbeame,        &
                                        dudtsv, sbeam, qrad_tot, brems_tot

  USE iterdbmd_gcnmp,            ONLY : iterdsc,irwflag,iterdb_file_name,            &
                                        iterdb_outpt,ncid

  USE ions_gcnmp,                ONLY : namep,namei,nion,        &
                                        namen,nprim,nimp,nneu,z,zsq,zeff,            &
                                        name_size,fd_thermal,nprimp1

  USE MPI_data,                  ONLY : mpiierr,myid,master



  USE solcon_gcnmp,              ONLY : time,eqtime,tGCNMf,tGCNMs

  USE fusion_gcnmp,              ONLY : pfuse_tot,he3_frac,he3_thermal_spin_pol

  USE neutral_data,              ONLY : enn,ennw,volsn,ennv

  USE fast_ion_data_gcnmp,       ONLY : walp,enalp,w_alpha

  USE neutral_beams,             ONLY : nbion,fd_beam,nameb,enbeam,storqueb, &
                                        wbeam, bptor,enbeam_tot,             &
                                        neut_beam_allocate

  USE curden_terms,              ONLY : currf,curbeam,ibcur,irfc

  USE common_constants,          ONLY : convert,zeroc,izero

  USE bc_values_gcnmp,           ONLY : totcur_bc,time_bc,                 &
                                        fix_edge_te_bc,fix_edge_ti_bc,     &
                                        fix_edge_rot_bc,fix_edge_ni_bc,    &
                                        zeff_bc,ene_bc,te_bc,ti_bc,        &
                                        angrot_bc,en_bc,vloop_bc,flux_bc,  &
                                        convert_edge_indecies

  USE dep_var,                   ONLY : dp4

  USE fdyimp,                    ONLY : dfdt,dgdt,dhdt

  USE tglfin,                    ONLY : tglf_p_output,tglf_e_output,tglf_m_output

 USE glf23_gcnmp,                ONLY : glf_etg_output,glf_gamma_net_i_output,               &
                                        glf_gamma_net_e_output,glf_anfreq_output,            &
                                        glf_anfreq2_output,glf_anrate_output,                &
                                        glf_anrate2_output,                                  &
                                        glf_p_output,glf_e_output,glf_m_output,              &
                                        glf_eff_chi_output

  
  !USE tcoef,                    ONLY : chie_paleo





  IMPLICIT  NONE


  REAL(DP), ALLOCATABLE,DIMENSION(:)   :: work_nj,               &
       work_npsi,             &
       work_nplasbdry,        &
       work_nlimiter,         &
       work_nr,               &
       work_nz,               &
       work_bptor
  REAL(DP), ALLOCATABLE,DIMENSION(:,:)   :: work_nj_nion,work_nj_nionp1
  REAL(DP), ALLOCATABLE,DIMENSION(:,:,:) :: work3d
  INTEGER(I4B), ALLOCATABLE,DIMENSION(:) :: id_glf_eff_chi,id_glf_eff_chi_err
  
  REAL(DP) elmt
  INTEGER(I4B)  j,jp,jj,ji,jn,ntot,titlen,flag,id_brems_flag,         &
                id_qrad_tot_flag,id_brems_tot_flag,id_pfuse_tot_flag, &
                id_curpar_flag,id_qpsinpsi_flag,id_ffprimnpsi_flag,   &
                id_pressnpsi_flag,id_pprimnpsi_flag,id_omegale_flag,  &
                id_qomegapi_flag,id_qangce_flag,id_sprcxre_flag,      &
                id_sprcxree_flag,id_spreimpe_flag,ngW20,oknf,nionp1,ierr
  INTEGER(I2B) i,errflg,l1,l2,l3
  INTEGER(I4B), PARAMETER :: one = 1_I4B
  INTEGER(I4B), PARAMETER :: nf_max_name = 128 ! should be nf90_max_name but cant find
  !any definition of this in netcdf modules
  !(v3.6)
  CHARACTER(len = nf_max_name) :: gen_name=' '
  CHARACTER(len = name_size + 12) ::  omega_pi_name,omega_ci_name
  CHARACTER(len = name_size + 12) ::  omega_lh_name,omega_uh_name
  CHARACTER(len = name_size + 12) ::  omega_ce_name,omega_pe_name

  INTEGER(I4B) k
  INTEGER(I4B)              &         ! netcdf required values
       id0,                 &         ! used for intent inout calls below
       rcode,               &
       id_shot,             &
       id_nj,               &
       id_time,             &
       id_rgeom,            &
       id_btgeom,           &
!       id_rmag,            &
       id_rma,              &
       id_zma,              &
       id_rsep,             &
       id_zsep,             &
       id_rmajor,           &
       id_rplasmin,         &
       id_rplasmax,         &
       id_zplasmin,         &
       id_zplasmax,         &
       id_kappa,            &
       id_psiaxis,          &
       id_psibdry,          &
       id_deltao,           &
       id_pindento,         &
       id_volume,           &
       id_circum,           &
       id_areao,            &
       id_nion,             &
       id_nprim,            &
       id_nneu,             &
       id_enn,              &
       id_ennw,             &
       id_ennv,             &
       id_nbion ,           &
       id_namep,            &
       id_namen,            &
       id_namepel,          &
       id_tgcnmf,           &
       id_time_bc,          &
       id_fd_thermal,       &
       id_nimp,             &
       id_fd_beam,          &
       id_namei,            &
       id_nameb,            &
       id_fpsinpsi,         &
       id_psivalnpsi,       &
       id_ravgnpsi,         &
       id_ravginpsi,        &
       id_ffprim,           &
       id_pprim,            &
       id_btor,             &
       id_fbcur,            &
       id_fbcur_er,         &
       id_prompt_nb_pwr,    &
       id_prompt_nb_pwr_er, &
       id_bp,               &
       id_bprmaj,           &
       id_btotrmaj,         &
       id_tot_cur,          &
       id_totohm_cur,       &
       id_totboot_cur,      &
       id_totbeam_cur,      &
       id_totrf_cur,        &
       id_betap,            &
       id_beta,             &
       id_ali,              &
       id_te0,              &
       id_ti0,              &
       id_psi,              &
       id_qpsinpsi,         &
       id_pressnpsi,        &
       id_ffprimnpsi,       &
       id_pprimnpsi,        &
       id_psir_grid,        &
       id_rho_grid,         &
       id_rho_mhd_gridnpsi, &
       id_rmhdgrid,         &
       id_zmhdgrid,         &
       id_fcap,             &
       id_gcap,             &
       id_hcap,             &
       id_betan,            &
       id_betan_er,         &
       id_te,               &
       id_ti,               &
       id_q_value,          & 
       id_ene,              &
       id_enion,            &
       id_flux_elct,        &
       id_flux_elct_er,     &
       id_flux_ion,         &
       id_flux_ion_er,      &
       id_pflux,            &
       id_e_fluxe,           &
       id_e_fluxi,           &
       id_fdyflux,          &
       id_rotflux,          &
       id_pflux_conv,       &
       id_pflux_conv_er,    &
       id_e_fluxe_conv,      &
       id_e_fluxe_conv_er,   &
       id_e_fluxi_conv,      &
       id_e_fluxi_conv_er,   &
       id_fdyflux_conv,     &
       id_fdyflux_conv_er,  &
       id_rotflux_conv,     &
       id_rotflux_conv_er,  &
       id_tglf_p_fluxe,     &
       id_tglf_e_fluxe,     &
       id_tglf_m_fluxe,     & !toroidal  momentum flux electron ion
       id_tglf_m_fluxe_er,  & ! error ahndler
       id_tglf_p_fluxp,     &
       id_tglf_e_fluxp,     &
       id_tglf_m_fluxp,     & !toroidal  momentum flux primary ion
       id_tglf_m_fluxp_er,  & 
       id_tglf_p_fluxi,     &
       id_tglf_e_fluxi,     &
       id_tglf_m_fluxi,     & !toroidal  momentum flux impurity ion
       id_tglf_m_fluxi_er,  &
       id_sion,             &
       id_srecom,           &
       id_scx,              &
       id_sbcx,             &
       id_stsource,         &
       id_dudtsv,           &
       id_volsn,            &
       id_stfuse,           &
       id_sbfuse,           &
       id_sbeame,           &
       id_sbeam,            &
       id_spellet,          &
       id_curden,           &
       id_curpar,           &
       id_curohm,           &
       id_curboot,          &
       id_curbeam,          &
       id_currf,            &
       id_etor,             &
       id_rbp,              &
       id_zeff,             &
       id_vpol,id_vpol_er,  &
       id_vpol_nclass,      &
       id_vpol_nclass_er,   &
       id_vpar,id_vpar_er,  &
       id_vpar_nclass,      &
       id_vpar_nclass_er,   &
       id_er_tot_nclass,    &
       id_er_tot_nclass_er, &
       id_angrot,           &
       id_d,                &
       id_chieinv,          &
       id_chiinve,          &
       id_xkeneo,           &
       id_dpedt,            &
       id_dpidt,            &
       id_qconde,           &
       id_qcondi,           &
       id_qconve,           &
       id_qconvi,           &
       id_qbeame,           &
       id_qbeami,           &
       id_qdelt,            &
       id_qexch,            &
       id_qrfe,             &
       id_qrfi,             &
       id_qione,            &
       id_qioni,            &
       id_qcx,              &
       id_qe2d,             &
       id_qi2d,             &
       id_qfuse,            &
       id_qfusi,            &
       id_qbfue,            &
       id_qbfusi,           &
       id_qmag,             &
       id_qsawe,            &
       id_qsawi,            &
       id_qrad,             &
       id_omegale,          &
       id_qomegapi,         &
       id_qangce,           &
       id_sprcxre,          &
       id_spreimpe,         &
       id_sprcxree,         &
       id_brems_nions,      &
       id_vpinch_nions,     &
       id_vpinch_flag,      &
       id_qrad_tot,         &
       id_brems_tot,        &
       id_pfuse_tot,        &
       id_qohm,             &
       id_rmajavnpsi,       &
       id_chiinv,           &
       id_xkineo,           &
       id_qbfuse,           &
       id_qbfusei,          &
       id_rminavnpsi,       &
       id_psivolpnpsi,      &
       id_elongxnpsi,       &
       id_triangnpsi_u,     &
       id_triangnpsi_l,     &
       id_pindentnpsi,      &
       id_sfareanpsi,       &
       id_cxareanpsi,       &
       id_grho1npsi,        &
       id_grho2npsi,        &
       !id_nplasbdry,       &
       id_rplasbdry,        &
       id_zplasbdry,        &
       id_nlimiter,         &
       id_rlimiter,         &
       id_zlimiter,         &
       id_storqueb,         &
       id_totcur_bc,        &
       id_vloop_bc,         &
       id_fix_edge_te_bc,   &
       id_fix_edge_ti_bc,   &
       id_fix_edge_rot_bc,  &
       id_fix_edge_ni_bc,   &
       id_te_bc,            &
       id_ti_bc,            &
       id_press,            &
       id_pressb,           &
       id_ene_bc,           &
       id_zeff_bc,          &
       id_angrot_bc,        &
       id_en_bc,            &
       id_flux_bc,          &
       id_z,                &
       id_zsq,              &
       id_wbeam,            &
       id_walp,             &
       id_enalp,            &
       id_dnedt,            &
       id_eps,              &
       id_rcap,             &
       id_rcapi,            &
       id_r2capi,           &
       id_r2cap,            &
       id_xi11,             &
       id_xi33,             &
       id_xips,             &
       id_xhm2,             &
       id_enbeam,           &
       id_dfdt,             &
       id_dgdt,             &
       id_dhdt,             &
       idim_char,           &
       idim_rho,            &
       idim_beam_rho,       &
       idim_nr_mhd,         &
       idim_nz_mhd,         &
       idim_nion,           &
       idim_nionp1,         &
       idim_nprim,          &
       idim_nimp,           &
       idim_nneu,           &
       idim_nbion,          &
       idim_nplasbdry,      &
       idim_ntot,           &
       idim_npsi,           &
       idim_nlimiter,       &
       nDimensions,         &
       nVariables,          &
       nAttributes,         &
       unlimitedDimId,      &
       formatNum,           &
       idlen,               &
       index,               &
       id_ftrap,            &
       id_xnuse,            &
       id_xnus,             &
       id_eta,              &
       id_neutr_ddn_th,                     &
       id_neutr_ddn_beam_beam,              &
       id_neutr_ddn_beam_thermal,           &
       id_neutr_ddn_knock,                  &
       id_neutr_ddn_tot,                    &
       id_p_ddt,                            &
       id_p_he3d_th,                        &
       id_p_he3d_beam_th,                   &
       id_total_p_he3d_beam_th,             &
       id_total_p_ddt,                      &
       id_total_p_he3d_th,                  &
       id_total_p_he3d,                     &
       id_he3_frac,                         &
       id_he3_thermal_spin_pol,             &
       id_d_beam_spin_pol,                  &
       id_total_neutr_ddn_th,               &
       id_total_neutr_ddn_beam_beam,        &
       id_total_neutr_ddn_beam_thermal,     &
       id_total_neutr_ddn_knock,            &
       id_total_neutr_ddn,                  &
       id_nbeams,                           &
       id_nbeams_er,                        &
       id_ebeam,                            &
       id_pbeam,                            &
       id_bptor,                            &
       id_bneut,                            &
       id_bion,                             &
       id_fap,                              &
       id_forb,                             &
       id_fber,                             &
       id_fb00,                             &
       id_fb01,                             &
       id_fb10,                             &
       id_fb11,                             &
       id_wb00,                             &
       id_wb01,                             &
       id_wb10,                             &
       id_wb11,                             &
       id_fwall,                            &
       id_sb,                               &
       id_qb,                               &
       id_angmpf,                           &
       id_pb0,                              &
       id_spb,                              &
       id_spbr,                             &
       id_hibr,                             &
       id_hdep,                             &
       id_zeta,                             &
       id_hicme,                            &
       id_hicmp1,                           &
       id_hicmp2,                           &
       idim_nbeams,                         &
       idim_nbeams_er,                      &  ! _er used for error handling
       id_rhog_beam,id_rhog_beam_er,        &
       idim_ke_bm,                          &
       idim_ke_bm_er,                       &
       id_chiepc,                           &
       id_glf_etg_flux,id_glf_etg_flux_er,  &
       id_glf_gam_net_i,id_glf_gam_net_i_er,&
       id_glf_gam_net_e, id_glf_gam_net_e_er,                     &
       id_glf_anfreq,id_glf_anfreq_er ,                           &
       id_glf_anfreq2, id_glf_anfreq2_er ,                        &
       id_glf_anrate,   id_glf_anrate_er ,                        &
       id_glf_anrate2,  id_glf_anrate2_er ,                       &
       id_glf_primion_eng_flux,id_glf_primion_eng_flux_er ,       &
       id_glf_impion_eng_flux,id_glf_impion_eng_flux_er ,         &
       id_glf_elct_eng_flux,id_glf_elct_eng_flux_er ,             &
       id_glf_elct_momtm_flux,id_glf_elct_momtm_flux_er ,         &
       id_glf_primion_momtm_flux,id_glf_primion_momtm_flux_er ,   &
       id_glf_impion_momtm_flux,  id_glf_impion_momtm_flux_er ,   &
       id_glf_elct_partcl_flux, id_glf_elct_partcl_flux_er ,      &
       id_glf_primion_partcl_flux,id_glf_primion_partcl_flux_er , &
       id_glf_impion_partcl_flux,id_glf_impion_partcl_flux_er,    &               
       id_mmm_gammaDBM,id_mmm_gammaDBM_er,                        &
       id_mmm_omegaDBM,id_mmm_omegaDBM_er,                        &
       id_mmm_xdi,id_mmm_xdi_er,                                  &
       id_mmm_xti, id_mmm_xti_er,                                 &
       id_mmm_xte,id_mmm_xte_er,                                  &
       id_mmm_xdz,id_mmm_xdz_er,                                  &
       id_mmm_xvt,id_mmm_xvt_er,                                  &
       id_mmm_xvp,id_mmm_xvp_er,                                  &
       id_mmm_xtiW20,id_mmm_xtiW20_er,                            &
       id_mmm_xdiW20,id_mmm_xdiW20_er,                            &
       id_mmm_xteW20,id_mmm_xteW20_er,                            &
       id_mmm_xtiDBM,id_mmm_xtiDBM_er,                            &
       id_mmm_xdiDBM,id_mmm_xdiDBM_er,                            &
       id_mmm_xteDBM,id_mmm_xteDBM_er,                            &
       id_mmm_xteETG,id_mmm_xteETG_er,                            &
       id_mmm_gamma_i1_W20,id_mmm_gamma_i1_W20_er,                &
       id_mmm_gamma_e1_W20,id_mmm_gamma_e1_W20_er,                &
       id_mmm_gamma_i2_W20,id_mmm_gamma_i2_W20_er,                &
       id_mmm_gamma_e2_W20,id_mmm_gamma_e2_W20_er,                &
       id_mmm_omega_i1_W20,id_mmm_omega_i1_W20_er,                &
       id_mmm_omega_e1_W20,id_mmm_omega_e1_W20_er,                &
       id_mmm_omega_i2_W20,id_mmm_omega_i2_W20_er,                &
       id_mmm_omega_e2_W20,id_mmm_omega_e2_W20_er,                &
       id_mmm_flux_ith,id_mmm_flux_ith_er,                        &
       id_mmm_flux_ip,id_mmm_flux_ip_er,                          &
       id_mmm_flux_eth,id_mmm_flux_eth_er,                        &
       id_mmm_flux_imp,id_mmm_flux_imp_er,                        &
       id_mmm_vconv_ith,id_mmm_vconv_ith_er,                      &
       id_mmm_vconv_ip,id_mmm_vconv_ip_er,                        &
       id_mmm_vconv_eth,id_mmm_vconv_eth_er,                      &
       id_mmm_vconv_imp,id_mmm_vconv_imp_er,                      &
       id_mmm_vmtmt,id_mmm_vmtmt_er,                              &
       id_mmm_vmtmp,id_mmm_vmtmp_er,                              &
       id_ce_freq,id_pe_freq

 




  INTEGER(I4B),PARAMETER :: max_dims = 17
  INTEGER(I4B),PARAMETER :: dim_size = 16
  CHARACTER(LEN=dim_size) dim_names(max_dims)
  LOGICAL    monotonic
  INTEGER(I4B),PARAMETER :: strlen = 132
  INTEGER(I4B),PARAMETER :: ke_bm     = 3         ! # beam energies
  INTEGER(I4B) ke ! used to check ke_bm
  INTEGER(I4B) ok
  CHARACTER  starflag*2, headerline*132,line*132,time_str*36
  CHARACTER  (LEN= strlen) label ,tlabel,base_label
  CHARACTER   bc_asc_time*24,st_asc_time*24 
  CHARACTER(len = name_size) :: tname
  CHARACTER(LEN = strlen) :: title
  CHARACTER(LEN = *), PARAMETER :: dimensionless = 'dimensionless'

  DATA dim_names /"dim_rho","dim_ion","dim_nprim","dim_char",                &
       "dim_nimp", "dim_neu","dim_fi","dim_nplasbdry",            &
       "dim_ntot","dim_nr_mhd","dim_nz_mhd","dim_npsi",           &
       "dim_nlimiter","dim_nbeams","dim_ke", "dim_nionp1","dim_beam_rho"/
  INTERFACE 
     SUBROUTINE check_monotonic (array, n, monotonic, incr) 
       USE nrtype,   ONLY : I4B,DP
       INTEGER(I4B), INTENT(IN) :: n,incr
       REAL(DP),INTENT(IN),DIMENSION(:) :: array(n)
       LOGICAL, INTENT(OUT) ::    monotonic
     END SUBROUTINE check_monotonic

     SUBROUTINE netcdf_err(status,flag)
       USE nrtype,   ONLY : I4B,DP
       IMPLICIT NONE
       INTEGER, INTENT ( in) :: status
       INTEGER(I4B), OPTIONAL, INTENT(inout)     :: flag
       
     END SUBROUTINE netcdf_err


  END INTERFACE

  iterdsc = 1                   ! always write descriptor,reads will fail otherwise
  
7 FORMAT (         a   )
8 FORMAT (5(2x,    a  ))        ! common character write/read format
9 FORMAT (5(2x,   i6  ))        ! common integer   write/read format
10 FORMAT (5(2x,1pe14.4))        ! common floating  write/read format
11 FORMAT (2x,i6,2x,1pe14.4)

  !     IF(ALLOCATED(work))DEALLOCATE(work)
  !     ALLOCATE (work(MAX(nj,mhd_dat%npsi,dischg%nplasbdry,dischg%nlimiter))) 
  !     universal array wont work, netcdf checks sizes internally:
  IF(ALLOCATED(work_nj))DEALLOCATE(work_nj)
  IF(ALLOCATED(work_npsi))DEALLOCATE(work_npsi)
  IF(ALLOCATED(work_nplasbdry))DEALLOCATE(work_nplasbdry)
  IF(ALLOCATED(work_nlimiter))DEALLOCATE(work_nlimiter)
  IF(ALLOCATED(work_nr))DEALLOCATE(work_nr)
  IF(ALLOCATED(work_nz))DEALLOCATE(work_nz)
  IF(ALLOCATED(work_nj_nion))DEALLOCATE(work_nj_nion)
  IF(ALLOCATED(work_nj_nionp1))DEALLOCATE(work_nj_nionp1)
  IF(ALLOCATED(work_bptor))DEALLOCATE(work_bptor)
 

  !If input statefile was in text form these arrays were not created.
  !nprimp1 is known in this case and we are can allocate the id_** arrays
  !so that the netcdf output section will work. 
  !If input statefile is in netcdf form then nprimp1 is not yet known so
  ! we have to allocate these arrays below
  IF(.NOT. ALLOCATED(id_plas_freq) .AND. nprimp1 .GT. 1) &
       CALL allocate_plasma_freq



  !
  IF (irwflag .EQ. 0) THEN      ! WRITE to netcdf file 
     iterdb_outpt = iterdb_file_name(1:LEN_TRIM(iterdb_file_name))
     rcode  = NF90_CREATE ( iterdb_outpt,NF90_CLOBBER, ncid) !enters define mode
     IF (rcode .NE.  nf90_noerr )THEN
        lerrno = 3_I4B
        CALL  terminate(lerrno,nlog)
     ENDIF
     WRITE(label,FMT='(1pe16.8)')time_bc          !get  bc time in string form for output
     READ(label,FMT='(a)')bc_asc_time
     bc_asc_time = ADJUSTL(bc_asc_time)
  ELSE                     ! READ from existing file 

     rcode =  NF90_OPEN (iterdb_file_name,NF90_NOWRITE, ncid)
     IF (rcode == nf90_noerr ) go to 3
2    WRITE  (ncrt, 4)  iterdb_file_name(1:LEN_TRIM(iterdb_file_name))
     WRITE  (nlog, 4)  iterdb_file_name(1:LEN_TRIM(iterdb_file_name))
4    FORMAT (                                                     / &
          ' ERROR: subroutine ITER_DBASE_NC has encountered an error' / &
          '        the ITER database file "', a, '" cannot be opened')
     lerrno = 3_I4B
     CALL  terminate(lerrno,nlog)
3    CONTINUE

  END IF




  !   ------------------------------------------------------------------

  !     Netcdf file is ready for reading or writing:

  !   ------------------------------------------------------------------

  rwncd:  IF (irwflag .EQ. 0) THEN       ! write the data

     ALLOCATE(work_nj(nj))
     ALLOCATE(work_npsi(mhd_dat%npsi))
     ALLOCATE(work_nplasbdry(dischg%nplasbdry))
     ALLOCATE(work_nlimiter(dischg%nlimiter))
     ALLOCATE(work_nr(dischg%nr_mhd))
     ALLOCATE(work_nz(dischg%nz_mhd))
     ALLOCATE(work_nj_nion(nj,nion))
     ALLOCATE(work_nj_nionp1(nj,nion+1))
     ALLOCATE(work_bptor(1:neut_beam%nbeams))
     !     in the future we may introduce time dependent metrics.
     !     then dfdt (= d/dt FCAP ,etc) will be calculated internally in this code
     !     by reading in fcap_bc, etc.
     !     in demo code just assume time independent for each
     !     slice that the solver is called.
     ntot = nion + dp4
     IF(.NOT. ALLOCATED(id_glf_eff_chi))THEN
         ALLOCATE(id_glf_eff_chi(ntot),id_glf_eff_chi_err(ntot))
     ELSE
         IF(SIZE(id_glf_eff_chi) .NE. ntot)THEN
            DEALLOCATE(id_glf_eff_chi,id_glf_eff_chi_err)
            ALLOCATE(id_glf_eff_chi(ntot),id_glf_eff_chi_err(ntot))
         ENDIF
     ENDIF
     IF(.NOT. ALLOCATED(dfdt))THEN 
        ALLOCATE(dfdt(nj))
        dfdt(:) = zeroc
     ENDIF
     IF(.NOT. ALLOCATED(dgdt))THEN 
        ALLOCATE(dgdt(nj))
        dgdt(:) = zeroc
     ENDIF
     IF(.NOT. ALLOCATED(dhdt))THEN
        ALLOCATE(dhdt(nj))
        dhdt(:) = zeroc
     ENDIF
     ! title of dataset:
#ifdef USESUBS
#ifdef GCNMP
     title = 'GCNMP  netCDF output file '//gcnmp_ver(1:LEN_TRIM(gcnmp_ver))
#endif
#ifdef NFREYA
     title = 'P_Nfreya created statefile'//p_nfreya_ver(1:LEN_TRIM(p_nfreya_ver))
#endif
#else
     title = 'ONETWO  netCDF state file '//versid
#endif

     CALL netcdf_err( nf90_put_att(ncid,  NF90_GLOBAL,"title",title))

     ! define dimensions:

     CALL netcdf_err( nf90_def_dim(ncid, dim_names(1), nj, idim_rho),idim_rho)                      !"dim_rho" 

     CALL netcdf_err( nf90_def_dim(ncid, dim_names(2), nion, idim_nion),idim_nion)                   !"dim_ion"

     CALL netcdf_err( nf90_def_dim(ncid, dim_names(3), nprim, idim_nprim),idim_nprim)                 !"dim_nprim"

     CALL netcdf_err( nf90_def_dim(ncid, dim_names(4), name_size, idim_char),idim_char)              !"dim_char"

     CALL netcdf_err( nf90_def_dim(ncid, dim_names(5), nimp, idim_nimp),idim_nimp)                   !"dim_nimp"

     IF ( nneu > 0 ) THEN 
        
        CALL netcdf_err( nf90_def_dim(ncid, dim_names(6), nneu, idim_nneu),idim_nneu)                   !"dim_neu"
        
     ENDIF

     CALL netcdf_err( nf90_def_dim(ncid, dim_names(7), nbion, idim_nbion),idim_nbion)                 !"dim_fi"

     CALL netcdf_err( nf90_def_dim(ncid, dim_names(8), dischg%nplasbdry, &
          idim_nplasbdry),idim_nplasbdry)                   !"dim_nplasbdry"
     ntot = nion + dp4
     CALL netcdf_err( nf90_def_dim(ncid, dim_names(9), ntot, idim_ntot),idim_ntot) !"dim_ntot"

     CALL netcdf_err( nf90_def_dim(ncid, dim_names(10),dischg%nr_mhd, idim_nr_mhd),idim_nr_mhd) !"dim_nr_mhd"
     CALL netcdf_err( nf90_def_dim(ncid, dim_names(11),dischg%nz_mhd, idim_nz_mhd),idim_nz_mhd) !"dim_nz_mhd"

     CALL netcdf_err( nf90_def_dim(ncid, dim_names(12),mhd_dat%npsi, idim_npsi),idim_npsi) !"dim_npsi"
     CALL netcdf_err( nf90_def_dim(ncid, dim_names(13),dischg%nlimiter, idim_nlimiter),idim_nlimiter) !"dim_nlimiter"


!-----------------------------------------------------------------------------------------------------------
!beams may or may not be present
!-----------------------------------------------------------------------------------------------------------

     idim_nbeams_er = izero  ! ==> zero action in sub netcdf_err
     ! if the function nf90_def_dim returns an error idim_nbeams_er = -1 after call to netcdf_err 

     CALL netcdf_err( nf90_def_dim(ncid, dim_names(14),neut_beam%nbeams, idim_nbeams),idim_nbeams_er) !"dim_nbeams"

      IF(idim_nbeams_er == -1)neut_beam%nbeams =0 ! use this to discover that neut_beam%nbeams is not set

     idim_ke_bm_er = izero      ! ==> zero action in sub netcdf_err
     CALL netcdf_err( nf90_def_dim(ncid, dim_names(15),ke_bm, idim_ke_bm),idim_ke_bm_er) !"dim_ke"

     nionp1 = nion +1
     CALL netcdf_err( nf90_def_dim(ncid, dim_names(16), nionp1, idim_nionp1),idim_nionp1) !"dim_nionp1"
 
     ! rho grid dimension of beam related items  - these values are not regridded so need to keep
     ! separate rho grid and dim for them.  nj_beam is defined by codes that call the beam
     ! and passed to here:
     IF(neut_beam%nj_beam  .LE. 0 .AND. neut_beam%nbeams .GT. 0)neut_beam%nj_beam = nj

     IF(neut_beam%nbeams .NE. 0)& ! if neut_beam%nbeams = 0 then  dim_names(14) is defined as unlimited dimension
                                  ! hence cant define dim_names(17) also with 0 dim
                                  ! this happens even if dim_names .ne. 'unlim'
     CALL netcdf_err( nf90_def_dim(ncid, dim_names(17), neut_beam%nj_beam, idim_beam_rho),idim_beam_rho) 


!   print *,'idim_beam_rho =',neut_beam%nj_beam ! 17
!   print *,'idim_nionp1 =',nionp1     ! 16
!   print *,'idim_ke_bm =',ke_bm       ! 15
!   print *,'idim_nbeams=',neut_beam%nbeams     ! 14
!   print *,'idim_nlimiter=',dischg%nlimiter  ! 13
!   print *,'idim_npsi=',mhd_dat%npsi         ! 12
!   print *,'idim_nz_mhd=', dischg%nz_mhd   !11
!   print *,'idim_nr_mhd=',dischg%nr_mhd     !10
!   print *,'idim_ntot=',ntot         !9
!   print *,'idim_nplasbdry=',dischg%nplasbdry ! 8
!   print *,'idim_nbion=', nbion        ! 7 
!   print *,'idim__nneu=',nneu         ! 6
!   print *,'idim_nimp=',nimp           !5
!   print *,'idim_char,=',name_size         ! 4
!   print *,'idim_nprim=', nprim         ! 3
!   print *,'idim_nion =',nion         ! 2
!   print *,'idim_rho=', nj          ! 1


     ! define variables:

     CALL netcdf_err( nf90_def_var(ncid, "shot", nf90_int,id_shot),id_shot)
     CALL netcdf_err(  nf90_put_att(ncid,id_shot,'long_name','shot number'),id_shot )
     CALL netcdf_err(  nf90_put_att(ncid,id_shot,'units', dimensionless  ))

     CALL netcdf_err( nf90_def_var(ncid, "nj", nf90_int, id_nj))
     CALL netcdf_err(  nf90_put_att(ncid,id_nj,'long_name',           &
          'the size of quantities defined on rho grid ') )

     CALL netcdf_err( nf90_def_var(ncid, "time", nf90_double,id_time))
     CALL netcdf_err(  nf90_put_att(ncid,id_time,'long_name',           &
          'time at which data is printed ') )

     CALL netcdf_err( nf90_def_var(ncid, "tGCNMf", nf90_double,id_tGCNMf))
     CALL netcdf_err(  nf90_put_att(ncid,id_tGCNMf,'long_name',           &
          'GCNMP will evolve solution from time to tGCNMf (unless changed in namelist) ') )

     CALL netcdf_err( nf90_def_var(ncid, "time_bc", nf90_double,id_time_bc))
     CALL netcdf_err(  nf90_put_att(ncid,id_time_bc,'long_name',           &
          'Boundary condition time') )

 

     CALL netcdf_err( nf90_def_var(ncid, "psiaxis", nf90_double,id_psiaxis))
     CALL netcdf_err(  nf90_put_att(ncid,id_psiaxis,'long_name',           &
          'psi value on magnetic axis, Volt sec/rad') )
     CALL netcdf_err(  nf90_put_att(ncid,id_psiaxis,'units', 'volt second/radian'  ))


     CALL netcdf_err( nf90_def_var(ncid, "psibdry", nf90_double,id_psibdry))
     CALL netcdf_err(  nf90_put_att(ncid,id_psibdry,'long_name',           &
          'psi value on plasma edge (separatrix), volt sec/rad') )
     CALL netcdf_err(  nf90_put_att(ncid,id_psibdry,'units', 'volt second/radian'  ))


     CALL netcdf_err( nf90_def_var(ncid, "rgeom", nf90_double,id_rgeom))
     CALL netcdf_err(  nf90_put_att(ncid,id_rgeom,'long_name',           &
          'major radius of geometric center at elevation of magnetic axis') )
     CALL netcdf_err(  nf90_put_att(ncid,id_rgeom,'units', 'meter'  ))


     CALL netcdf_err( nf90_def_var(ncid, "btgeom", nf90_double,id_btgeom))
     CALL netcdf_err(  nf90_put_att(ncid,id_btgeom,'long_name',           &
          'toroidal b field at geometric center rgeom, Tesla') )
     CALL netcdf_err(  nf90_put_att(ncid,id_btgeom,'units', 'tesla'  ))


     !          CALL netcdf_err( nf90_def_var(ncid, "rmag", nf90_double,id_rmag))
     !          CALL netcdf_err(  nf90_put_att(ncid,id_rmag,'long_name',           &
     !           'major radius of magnetic axis') )
     !          CALL netcdf_err(  nf90_put_att(ncid,id_rmag,'units', 'meter'  ))

     CALL netcdf_err( nf90_def_var(ncid, "rma", nf90_double,id_rma))
     CALL netcdf_err(  nf90_put_att(ncid,id_rma,'long_name',           &
          'major radius of magnetic axis') )
     CALL netcdf_err(  nf90_put_att(ncid,id_rma,'units', 'meter'  ))

     CALL netcdf_err( nf90_def_var(ncid, "zma", nf90_double,id_zma))
     CALL netcdf_err(  nf90_put_att(ncid,id_zma,'long_name',           &
          'z of magnetic axis') )
     CALL netcdf_err(  nf90_put_att(ncid,id_zma,'units', 'meter'  ))

     CALL netcdf_err( nf90_def_var(ncid, "rsep", nf90_double,id_rsep))
     CALL netcdf_err(  nf90_put_att(ncid,id_rsep,'long_name',           &
          'r of separatrix x point, m') )
     CALL netcdf_err(  nf90_put_att(ncid,id_rsep,'units', 'meter'  ))


     CALL netcdf_err( nf90_def_var(ncid, "zsep", nf90_double,id_zsep))
     CALL netcdf_err(  nf90_put_att(ncid,id_zsep,'long_name',           &
          'z of separatrix x point, m') )
     CALL netcdf_err(  nf90_put_att(ncid,id_zsep,'units', 'meter'  ))


     CALL netcdf_err( nf90_def_var(ncid, "rmajor", nf90_double,id_rmajor))
     CALL netcdf_err(  nf90_put_att(ncid,id_rmajor,'long_name',           &
          'major radius of vacuum BT0 reference == R0') )
     CALL netcdf_err(  nf90_put_att(ncid,id_rmajor,'units', 'meter'  ))

     CALL netcdf_err( nf90_def_var(ncid, "rplasmin", nf90_double,id_rplasmin))
     CALL netcdf_err(  nf90_put_att(ncid,id_rplasmin,'long_name',        &
          'min R of plasma,m') )
     CALL netcdf_err(  nf90_put_att(ncid,id_rplasmin,'units', 'meter'  ))

     CALL netcdf_err( nf90_def_var(ncid, "rplasmax", nf90_double,id_rplasmax))
     CALL netcdf_err(  nf90_put_att(ncid,id_rplasmax,'long_name',        &
          'max R of plasma,m') )
     CALL netcdf_err(  nf90_put_att(ncid,id_rplasmax,'units', 'meter'  ))

     CALL netcdf_err( nf90_def_var(ncid, "zplasmin", nf90_double,id_zplasmin))
     CALL netcdf_err(  nf90_put_att(ncid,id_zplasmin,'long_name',        &
          'min Z of plasma,m') )
     CALL netcdf_err(  nf90_put_att(ncid,id_zplasmin,'units', 'meter'  ))

     CALL netcdf_err( nf90_def_var(ncid, "zplasmax", nf90_double,id_zplasmax))
     CALL netcdf_err(  nf90_put_att(ncid,id_zplasmax,'long_name',        &
          'max Z of plasma,m') )
     CALL netcdf_err(  nf90_put_att(ncid,id_zplasmax,'units', 'meter'  ))



     CALL netcdf_err( nf90_def_var(ncid, "kappa", nf90_double,id_kappa))
     CALL netcdf_err(  nf90_put_att(ncid,id_kappa,'long_name',           &
          'plasma elongation') )
     CALL netcdf_err(  nf90_put_att(ncid,id_kappa,'units', dimensionless  ))

     CALL netcdf_err( nf90_def_var(ncid, "deltao", nf90_double,id_deltao))
     CALL netcdf_err(  nf90_put_att(ncid,id_deltao,'long_name',           &
          '*  deltao  : plasma(upper ) triangularity on axis ') )
     CALL netcdf_err(  nf90_put_att(ncid,id_deltao,'units', dimensionless  ))


     CALL netcdf_err( nf90_def_var(ncid, "pindento", nf90_double,id_pindento))
     CALL netcdf_err(  nf90_put_att(ncid,id_pindento,'long_name',           &
          '*  pindento  : on axis plasma indentation  ') )
     CALL netcdf_err(  nf90_put_att(ncid,id_pindento,'units', dimensionless  ))


     CALL netcdf_err( nf90_def_var(ncid, "volume", nf90_double,id_volume))
     CALL netcdf_err(  nf90_put_att(ncid,id_volume,'long_name',           &
          'plasma volume') )
     CALL netcdf_err(  nf90_put_att(ncid,id_volume,'units', 'meter^3'  ))

     CALL netcdf_err( nf90_def_var(ncid, "circum", nf90_double,id_circum))
     CALL netcdf_err(  nf90_put_att(ncid,id_circum,'long_name',           &
          'plasma crcumference, meter') )
     CALL netcdf_err(  nf90_put_att(ncid,id_circum,'units', 'meter'  ))


     CALL netcdf_err( nf90_def_var(ncid, "areao", nf90_double,id_areao))
     CALL netcdf_err(  nf90_put_att(ncid,id_areao,'long_name',           &
          'plasma cross sectional area') )
     CALL netcdf_err(  nf90_put_att(ncid,id_areao,'units', 'meter^2'  ))

     CALL netcdf_err( nf90_def_var(ncid, "nion", nf90_int, id_nion))
     CALL netcdf_err(  nf90_put_att(ncid,id_nion,'long_name',      &
          'the number of ion  species') )
     CALL netcdf_err(  nf90_put_att(ncid,id_nion,'units', dimensionless  ))

     CALL netcdf_err( nf90_def_var(ncid, "nprim", nf90_int, id_nprim))
     CALL netcdf_err(  nf90_put_att(ncid,id_nprim,'long_name',      &
          'the number of primary ion species') )
     CALL netcdf_err(  nf90_put_att(ncid,id_nprim,'units', dimensionless  ))

     CALL netcdf_err( nf90_def_var(ncid, "fd_thermal", nf90_double, id_fd_thermal))
     CALL netcdf_err(  nf90_put_att(ncid,id_fd_thermal,'long_name',      &
          'fraction of d in thermal dt mixture if input as one species') )
     CALL netcdf_err(  nf90_put_att(ncid,id_fd_thermal,'units', dimensionless  ))

     CALL netcdf_err( nf90_def_var(ncid, "nimp", nf90_int, id_nimp))
     CALL netcdf_err(  nf90_put_att(ncid,id_nimp,'long_name',      &
          'the number of impurity ion species') )
     CALL netcdf_err(  nf90_put_att(ncid,id_nimp,'units', dimensionless  ))

     CALL netcdf_err( nf90_def_var(ncid, "nneu", nf90_int, id_nneu))
     CALL netcdf_err(  nf90_put_att(ncid,id_nneu,'long_name',      &
          'the number of neutral species') )
     CALL netcdf_err(  nf90_put_att(ncid,id_nneu,'units', dimensionless  ))

     CALL netcdf_err( nf90_def_var(ncid, "nbion", nf90_int, id_nbion))
     CALL netcdf_err(  nf90_put_att(ncid,id_nbion,'long_name',      &
          'number of fast beam ion  species ') )
     CALL netcdf_err(  nf90_put_att(ncid,id_nbion,'units', dimensionless  ))

     CALL netcdf_err( nf90_def_var(ncid, "fd_beam", nf90_double, id_fd_beam))
     CALL netcdf_err(  nf90_put_att(ncid,id_fd_beam,'long_name',      &
          'fraction of d in dt fast ion mixture') )
     CALL netcdf_err(  nf90_put_att(ncid,id_fd_beam,'units', dimensionless  ))

     CALL netcdf_err( nf90_def_var(ncid, "nbeams", nf90_int, id_nbeams))
     CALL netcdf_err(  nf90_put_att(ncid,id_nbeams,'long_name',      &
          'Number of neutral beam injectors in this dataset') )
     CALL netcdf_err(  nf90_put_att(ncid,id_nbeams,'units', dimensionless  ))


 

     CALL netcdf_err( nf90_def_var(ncid, "namep", nf90_char,        &
          DIMIDS = (/idim_char,idim_nprim/),VARID= id_namep))
     CALL netcdf_err(  nf90_put_att(ncid,id_namep,'long_name',      &
          'name of primary ion species') )
     CALL netcdf_err(  nf90_put_att(ncid,id_namep,'units', dimensionless  ))

     CALL netcdf_err( nf90_def_var(ncid, "namei", nf90_char,        &
          DIMIDS = (/idim_char,idim_nimp/),VARID= id_namei))
     CALL netcdf_err(  nf90_put_att(ncid,id_namei,'long_name',      &
          'name of impurity ion species') )
     CALL netcdf_err(  nf90_put_att(ncid,id_namei,'units', dimensionless  ))
     
     IF ( nneu > 0 ) THEN
        
        CALL netcdf_err( nf90_def_var(ncid, "namen", nf90_char,        &
             DIMIDS = (/idim_char,idim_nneu/),VARID=id_namen))
        CALL netcdf_err(  nf90_put_att(ncid,id_namen,'long_name',      &
             'name of neutral species') )
        CALL netcdf_err(  nf90_put_att(ncid,id_namen,'units', dimensionless  ))
     
     ENDIF

     CALL netcdf_err( nf90_def_var(ncid, "nameb", nf90_char,        &
          DIMIDS = (/idim_char,idim_nbion/),VARID=   id_nameb))
     CALL netcdf_err(  nf90_put_att(ncid,id_nameb,'long_name',      &
          'name of beam  species') )
     CALL netcdf_err(  nf90_put_att(ncid,id_nameb,'units', dimensionless  ))


     CALL netcdf_err( nf90_def_var(ncid, "namepel", nf90_char,        &
          DIMIDS = (/idim_char/),VARID=   id_namepel))
     CALL netcdf_err(  nf90_put_att(ncid,id_namepel,'long_name',      &
          'name of pellet  species') )
     CALL netcdf_err(  nf90_put_att(ncid,id_namepel,'units', dimensionless  ))



     CALL netcdf_err( nf90_def_var(ncid, "btor", nf90_double,id_btor))
     CALL netcdf_err(  nf90_put_att(ncid,id_btor,'long_name',           &
          '*  Btor : vacuum toroidal field at R0, tesla') )
     CALL netcdf_err(  nf90_put_att(ncid,id_btor,'units', 'tesla'  ))


     CALL netcdf_err( nf90_def_var(ncid, "tot_cur", nf90_double,id_tot_cur))
     CALL netcdf_err(  nf90_put_att(ncid,id_tot_cur,'long_name',           &
          '*  total plasma current,amps') )
     CALL netcdf_err(  nf90_put_att(ncid,id_tot_cur,'units', 'amp'  ))

     CALL netcdf_err( nf90_def_var(ncid, "totohm_cur", nf90_double,id_totohm_cur))
     CALL netcdf_err(  nf90_put_att(ncid,id_totohm_cur,'long_name',           &
          '* total ohmic plasma current, amps') )
     CALL netcdf_err(  nf90_put_att(ncid,id_totohm_cur,'units', 'amp'  ))


     CALL netcdf_err( nf90_def_var(ncid, "totboot_cur", nf90_double,id_totboot_cur))
     CALL netcdf_err(  nf90_put_att(ncid,id_totboot_cur,'long_name',           &
          '* total bootstrap  current, amps') )
     CALL netcdf_err(  nf90_put_att(ncid,id_totboot_cur,'units', 'amp'  ))

     CALL netcdf_err( nf90_def_var(ncid, "totbeam_cur", nf90_double,id_totbeam_cur))
     CALL netcdf_err(  nf90_put_att(ncid,id_totbeam_cur,'long_name',           &
          '* beam driven  current, amps') )
     CALL netcdf_err(  nf90_put_att(ncid,id_totbeam_cur,'units', 'amp'  ))

     CALL netcdf_err( nf90_def_var(ncid, "totrf_cur", nf90_double,id_totrf_cur))
     CALL netcdf_err(  nf90_put_att(ncid,id_totrf_cur,'long_name',           &
          '* rf  driven  current, amps') )
     CALL netcdf_err(  nf90_put_att(ncid,id_totrf_cur,'units', 'amp'  ))

     CALL netcdf_err( nf90_def_var(ncid, "betap", nf90_double,id_betap))
     CALL netcdf_err(  nf90_put_att(ncid,id_betap,'long_name',           &
          '*  betap : poloidal beta') )
     CALL netcdf_err(  nf90_put_att(ncid,id_betap,'units', dimensionless  ))


     CALL netcdf_err( nf90_def_var(ncid, "beta", nf90_double,id_beta))
     CALL netcdf_err(  nf90_put_att(ncid,id_beta,'long_name',           &
          '*  beta : toroidal beta') )
     CALL netcdf_err(  nf90_put_att(ncid,id_beta,'units', dimensionless  ))


     CALL netcdf_err( nf90_def_var(ncid, "ali", nf90_double,id_ali))
     CALL netcdf_err(  nf90_put_att(ncid,id_ali,'long_name',           &
          '*  ali : plasma inductance') )
     CALL netcdf_err(  nf90_put_att(ncid,id_ali,'units', dimensionless  ))


     CALL netcdf_err( nf90_def_var(ncid, "te0", nf90_double,id_te0))
     CALL netcdf_err(  nf90_put_att(ncid,id_te0,'long_name',           &
          '*  te0 : central electron temperature') )
     CALL netcdf_err(  nf90_put_att(ncid,id_te0,'units', 'keV'  ))


     CALL netcdf_err( nf90_def_var(ncid, "ti0", nf90_double,id_ti0))
     CALL netcdf_err(  nf90_put_att(ncid,id_ti0,'long_name',           &
          '*  ti0 : central ion  temperature') )
     CALL netcdf_err(  nf90_put_att(ncid,id_ti0,'units', 'keV'  ))


     CALL netcdf_err( nf90_def_var(ncid, "psi", nf90_double,    &
          DIMIDS = (/idim_nr_mhd,idim_nz_mhd/),VARID=id_psi))
     CALL netcdf_err( nf90_put_att(ncid,id_psi,'long_name',    &
          '*  psi on (R,Z) grid, volt*second/radian'))
     CALL netcdf_err(  nf90_put_att(ncid,id_psi,'units',       &
          'volt*second/radian'))


     CALL netcdf_err( nf90_def_var(ncid, "psir_grid", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_psir_grid))
     CALL netcdf_err( nf90_put_att(ncid,id_psir_grid,'long_name',    &
          '*  psi on rho grid, volt*second/radian'))
     CALL netcdf_err(  nf90_put_att(ncid,id_psir_grid,'units',       &
          'volt*second/radian'))

     CALL netcdf_err( nf90_def_var(ncid, "rho_mhd_gridnpsi", nf90_double,    &
          DIMIDS = (/idim_npsi/),VARID=id_rho_mhd_gridnpsi))
     CALL netcdf_err( nf90_put_att(ncid,id_rho_mhd_gridnpsi,'long_name',    &
          '* rho grid corresponding to mhd psival grid, meter '))
     CALL netcdf_err(  nf90_put_att(ncid,id_rho_mhd_gridnpsi,'units',       &
          'meter'))


     CALL netcdf_err( nf90_def_var(ncid, "rmhdgrid", nf90_double,    &
          DIMIDS = (/idim_nr_mhd/),VARID=id_rmhdgrid))
     CALL netcdf_err( nf90_put_att(ncid,id_rmhdgrid,'long_name',    &
          '* R grid, meter '))
     CALL netcdf_err(  nf90_put_att(ncid,id_rmhdgrid,'units',       &
          'meter'))

     CALL netcdf_err( nf90_def_var(ncid, "zmhdgrid", nf90_double,    &
          DIMIDS = (/idim_nz_mhd/),VARID=id_zmhdgrid))
     CALL netcdf_err( nf90_put_att(ncid,id_zmhdgrid,'long_name',    &
          '* Z grid, meter '))
     CALL netcdf_err(  nf90_put_att(ncid,id_zmhdgrid,'units',       &
          'meter'))

     CALL netcdf_err( nf90_def_var(ncid, "rho_grid", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_rho_grid))
     CALL netcdf_err( nf90_put_att(ncid,id_rho_grid,'long_name',    &
          '* rho grid, meter '))
     CALL netcdf_err(  nf90_put_att(ncid,id_rho_grid,'units',       &
          'meter'))


     CALL netcdf_err( nf90_def_var(ncid, "fcap", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_fcap))
     CALL netcdf_err( nf90_put_att(ncid,id_fcap,'long_name',    &
          '*  fcap, (i.e., f(psilim)/f(psi))'))
     CALL netcdf_err(  nf90_put_att(ncid,id_fcap,'units', dimensionless  ))

     CALL netcdf_err( nf90_def_var(ncid, "gcap", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_gcap))
     CALL netcdf_err( nf90_put_att(ncid,id_gcap,'long_name',    &
          '*  gcap, (i.e., <(grad rho)**2*(R0/R)**2>)'))
     CALL netcdf_err(  nf90_put_att(ncid,id_gcap,'units', dimensionless  ))

     CALL netcdf_err( nf90_def_var(ncid, "hcap", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_hcap))
     CALL netcdf_err( nf90_put_att(ncid,id_hcap,'long_name',    &
          '*  hcap, (i.e., (dvolume/drho)/(4*pi*pi*R0*rho))'))
     CALL netcdf_err(  nf90_put_att(ncid,id_hcap,'units', dimensionless  ))


     CALL netcdf_err( nf90_def_var(ncid, "betan", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_betan))
     CALL netcdf_err( nf90_put_att(ncid,id_betan,'long_name',    &
          '*  Normalized beta, P/BTGEOM^2/2u0/I/a/BTGEOM'))
     CALL netcdf_err(  nf90_put_att(ncid,id_fcap,'units', dimensionless  ))


     CALL netcdf_err( nf90_def_var(ncid, "Te", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_te))
     CALL netcdf_err( nf90_put_att(ncid,id_te,'long_name',    &
          '*  electron temperature, keV'  ))
     CALL netcdf_err(  nf90_put_att(ncid,id_te,'units',       &
          'keV'))

     CALL netcdf_err( nf90_def_var(ncid, "Ti", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_ti))
     CALL netcdf_err( nf90_put_att(ncid,id_ti,'long_name',    &
          '*  ion temperature, keV'  ))
     CALL netcdf_err(  nf90_put_att(ncid,id_ti,'units',       &
          'keV'))

     CALL netcdf_err( nf90_def_var(ncid, "press", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_press))
     CALL netcdf_err( nf90_put_att(ncid,id_press,'long_name',    &
          '*  total pressure on transport rho grid, nt/m^2'  ))
     CALL netcdf_err(  nf90_put_att(ncid,id_press,'units',       &
          'newton/meter^2'))

     CALL netcdf_err( nf90_def_var(ncid, "pressb", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_pressb))
     CALL netcdf_err( nf90_put_att(ncid,id_pressb,'long_name',    &
          '* beam  pressure on transport rho grid nt/m^2'  ))
     CALL netcdf_err(  nf90_put_att(ncid,id_pressb,'units',       &
          'newton/meter^2'))

     CALL netcdf_err( nf90_def_var(ncid, "q_value", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_q_value))
     CALL netcdf_err( nf90_put_att(ncid,id_q_value,'long_name',    &
          '*  safety factor '  ))
     CALL netcdf_err(  nf90_put_att(ncid,id_q_value,'units', dimensionless  ))

     CALL netcdf_err( nf90_def_var(ncid, "ene", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_ene))
     CALL netcdf_err( nf90_put_att(ncid,id_ene,'long_name',    &
          '*  electron density, #/meter^3'  ))
     CALL netcdf_err(  nf90_put_att(ncid,id_ene,'units',       &
          '1/meter^3'))

     CALL netcdf_err( nf90_def_var(ncid, "p_flux_elct", nf90_double, &
          DIMIDS = (/idim_rho/),VARID=id_flux_elct))
     CALL netcdf_err( nf90_put_att(ncid,id_flux_elct,'long_name',    &
          '*  electron  particle flux, #/meter^2 sec'  ))
     CALL netcdf_err(  nf90_put_att(ncid,id_flux_elct,'units',       &
          '1/(meter^s sec)'))


     CALL netcdf_err( nf90_def_var(ncid, "p_flux_ion", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_flux_ion))
     CALL netcdf_err( nf90_put_att(ncid,id_flux_ion,'long_name',    &
          '*  ion total particle flux, #/(meter^2 sec)'  ))
     CALL netcdf_err(  nf90_put_att(ncid,id_flux_ion,'units',       &
          '1/(meter^2 sec)'))

     !thermal ion densities:
     base_label = '* thermal ion densities, species: '
     CALL set_label(label,base_label,strlen,'nion')
     CALL netcdf_err( nf90_def_var(ncid, "enion", nf90_double,    &
          DIMIDS = (/idim_rho,idim_nion/),VARID=id_enion))
     CALL netcdf_err( nf90_put_att(ncid,id_enion,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_enion,'units',       &
          '1/meter^3'))

     !thermal ion particle flux:
     base_label = '* thermal ion particle flux, species: '
     CALL set_label(label,base_label,strlen,'nion')
     CALL netcdf_err( nf90_def_var(ncid, "p_flux", nf90_double,    &
          DIMIDS = (/idim_rho,idim_nion/),VARID=id_pflux))
     CALL netcdf_err( nf90_put_att(ncid,id_pflux,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_pflux,'units',       &
          '1/(meter^2 sec)'))

     !thermal ion convective particle flux:
     base_label = '* thermal ion convective flux, species: '
     CALL set_label(label,base_label,strlen,'nion')
     CALL netcdf_err( nf90_def_var(ncid, "p_flux_conv", nf90_double,    &
          DIMIDS = (/idim_rho,idim_nion/),VARID=id_pflux_conv))
     CALL netcdf_err( nf90_put_att(ncid,id_pflux_conv,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_pflux,'units',       &
          '1/(meter^2 sec)'))

     !electron energy  flux:
     base_label = '* electron energy flux :  '
     CALL set_label(label,base_label,strlen,'electron')
     CALL netcdf_err( nf90_def_var(ncid, "e_fluxe", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_e_fluxe))
     CALL netcdf_err( nf90_put_att(ncid,id_e_fluxe,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_e_fluxe,'units',       &
          'joule/(meter^2 sec)'))

    !electron convective energy  flux:
     base_label = '* electron convective energy flux :  '
     CALL set_label(label,base_label,strlen,'electron')
     CALL netcdf_err( nf90_def_var(ncid, "e_fluxe_conv", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_e_fluxe_conv))
     CALL netcdf_err( nf90_put_att(ncid,id_e_fluxe_conv,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_e_fluxe_conv,'units',       &
          'joule/(meter^2 sec)'))



     !glf electron energy  flux all modes, excludes neoclassical:
     base_label = '* total glf electron energy flux :  '
     CALL set_label(label,base_label,strlen,'electron')
     CALL netcdf_err( nf90_def_var(ncid, "glf_elct_eng_flux", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_glf_elct_eng_flux))
     CALL netcdf_err( nf90_put_att(ncid,id_glf_elct_eng_flux,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_glf_elct_eng_flux,'units',       &
          'joule/(meter^2 sec)'))

     !glf primary ion  energy  flux all modes, excludes neoclassical:
     base_label = '* total glf primary ion  energy flux :  '
     CALL set_label(label,base_label,strlen,'effective prim ion')
     CALL netcdf_err( nf90_def_var(ncid, "glf_primion_eng_flux", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_glf_primion_eng_flux))
     CALL netcdf_err( nf90_put_att(ncid,id_glf_primion_eng_flux,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_glf_primion_eng_flux,'units',       &
          'joule/(meter^2 sec)'))

     !glf impurity ion  energy  flux all modes, excludes neoclassical:
     base_label = '* total glf impurity ion  energy flux :  '
     CALL set_label(label,base_label,strlen,'effective impurity')
     CALL netcdf_err( nf90_def_var(ncid, "glf_impion_eng_flux", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_glf_impion_eng_flux))
     CALL netcdf_err( nf90_put_att(ncid,id_glf_impion_eng_flux,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_glf_impion_eng_flux,'units',       &
          'joule/(meter^2 sec)'))

     !glf electron momentum   flux all modes, excludes neoclassical:
     base_label = '* total glf elect momtm flux :  '
     CALL set_label(label,base_label,strlen,'electron')
     CALL netcdf_err( nf90_def_var(ncid, "glf_elct_momtm_flux", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_glf_elct_momtm_flux))
     CALL netcdf_err( nf90_put_att(ncid,id_glf_elct_momtm_flux,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_glf_elct_momtm_flux,'units',       &
          'kg/(sec^2)'))

     !glf primary ion  momentum   flux all modes, excludes neoclassical:
     base_label = '* total glf primary ion  momtm flux :  '
     CALL set_label(label,base_label,strlen,'effective prim ion')
     CALL netcdf_err( nf90_def_var(ncid, "glf_primion_momtm_flux", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_glf_primion_momtm_flux))
     CALL netcdf_err( nf90_put_att(ncid,id_glf_primion_momtm_flux,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_glf_primion_momtm_flux,'units',       &
          'kg/(sec^2)'))
 
     !glf imp ion  momentum   flux all modes, excludes neoclassical:
     base_label = '* total glf imp ion  momtm flux :  '
     CALL set_label(label,base_label,strlen,'effective impurity')
     CALL netcdf_err( nf90_def_var(ncid, "glf_impion_momtm_flux", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_glf_impion_momtm_flux))
     CALL netcdf_err( nf90_put_att(ncid,id_glf_impion_momtm_flux,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_glf_impion_momtm_flux,'units',       &
          'kg/(sec^2)'))

     !glf etg contribution to electron energy  flux :
     base_label = '* glf etg(high k) part of electron energy flux :  '
     CALL set_label(label,base_label,strlen,'glf-etg')
     CALL netcdf_err( nf90_def_var(ncid, "glf_etg_eng_flux", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_glf_etg_flux))
     CALL netcdf_err( nf90_put_att(ncid,id_glf_etg_flux,'long_name',label))
     CALL netcdf_err(  nf90_put_att(ncid,id_glf_etg_flux,'units',       &
          'joule/(meter^2 sec)'))

     !glf net growth rates :
     base_label = '* glf net growth rate ions:  '
     CALL set_label(label,base_label,strlen,'glf_gam_net_i')
     CALL netcdf_err( nf90_def_var(ncid, "glf_gamma_net_i", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_glf_gam_net_i))
     CALL netcdf_err( nf90_put_att(ncid,id_glf_gam_net_i,'long_name',label))
     CALL netcdf_err(  nf90_put_att(ncid,id_glf_gam_net_i,'units',       &
          '1/sec'))
     base_label = '* glf net growth rate electrons:  '
     CALL set_label(label,base_label,strlen,'glf_gam_net_e')
     CALL netcdf_err( nf90_def_var(ncid, "glf_gamma_net_e", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_glf_gam_net_e))
     CALL netcdf_err( nf90_put_att(ncid,id_glf_gam_net_e,'long_name',label))
     CALL netcdf_err(  nf90_put_att(ncid,id_glf_gam_net_e,'units',       &
          '1/sec'))


     base_label = '* glf leading mode growth rate in units of local csda_m:  '
     CALL set_label(label,base_label,strlen,'glf_anrate')
     CALL netcdf_err( nf90_def_var(ncid, "glf_anrate", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_glf_anrate))
     CALL netcdf_err( nf90_put_att(ncid,id_glf_anrate,'long_name',label))
     CALL netcdf_err(  nf90_put_att(ncid,id_glf_anrate,'units',       &
          'csda'))
     base_label = '* glf second  mode growth rate in units of local csda_m:  '
     CALL set_label(label,base_label,strlen,'glf_anrate2')
     CALL netcdf_err( nf90_def_var(ncid, "glf_anrate2", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_glf_anrate2))
     CALL netcdf_err( nf90_put_att(ncid,id_glf_anrate2,'long_name',label))
     CALL netcdf_err(  nf90_put_att(ncid,id_glf_anrate2,'units',       &
          'csda'))

     base_label = '* glf leading mode frequency:  '
     CALL set_label(label,base_label,strlen,'glf_anfreq')
     CALL netcdf_err( nf90_def_var(ncid, "glf_anfreq", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_glf_anfreq))
     CALL netcdf_err( nf90_put_att(ncid,id_glf_anfreq,'long_name',label))
     CALL netcdf_err(  nf90_put_att(ncid,id_glf_anfreq,'units',       &
          '1/sec'))
     base_label = '* glf second  mode frequency:  '
     CALL set_label(label,base_label,strlen,'glf_anfreq2')
     CALL netcdf_err( nf90_def_var(ncid, "glf_anfreq2", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_glf_anfreq2))
     CALL netcdf_err( nf90_put_att(ncid,id_glf_anfreq2,'long_name',label))
     CALL netcdf_err(  nf90_put_att(ncid,id_glf_anfreq2,'units',       &
          '1/sec'))

     !glf turbulent particle  flux:
     base_label = '* glf turbulent particle flux, species: '
     CALL set_label(label,base_label,strlen,'electron')
     CALL netcdf_err( nf90_def_var(ncid,"glf_elct_partcl_flux", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_glf_elct_partcl_flux))
     CALL netcdf_err( nf90_put_att(ncid,id_glf_elct_partcl_flux,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_glf_elct_partcl_flux,'units',       &
          '1/(meter^2 second)'))

    !glf turbulent particle  flux:
     base_label = '* glf turbulent particle flux, species: '
     CALL set_label(label,base_label,strlen,'effective prim ion')
     CALL netcdf_err( nf90_def_var(ncid, "glf_primion_partcl_flux", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_glf_primion_partcl_flux))
     CALL netcdf_err( nf90_put_att(ncid,id_glf_primion_partcl_flux,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_glf_primion_partcl_flux,'units',       &
          '1/(meter^2 second)'))

    !glf turbulent particle  flux:
     base_label = '* glf turbulent particle flux, species: '
     CALL set_label(label,base_label,strlen,'effective impurity')
     CALL netcdf_err( nf90_def_var(ncid, "glf_impion_partcl_flux", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_glf_impion_partcl_flux))
     CALL netcdf_err( nf90_put_att(ncid,id_glf_impion_partcl_flux,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_glf_impion_partcl_flux,'units',       &
          '1/(meter^2 second)'))


     ! glf eff primary ion  particle  diffusivites:

     DO i=1,nprim
        label = '* glf effective particle diffusivities,species :'//namep(i)
        base_label(1:LEN(base_label))=''
        base_label = "glf_eff_particle_chi_"//namep(i)
        CALL netcdf_err( nf90_def_var(ncid, base_label, nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_glf_eff_chi(i)))
        CALL netcdf_err( nf90_put_att(ncid,id_glf_eff_chi(i),'long_name',label ))
        CALL netcdf_err(  nf90_put_att(ncid,id_glf_eff_chi(i),'units',       &
          '1/(meter^2 second)'))
     ENDDO

     DO i=nprim+1,nion
        label = '* glf effective particle diffusivities,species :'//namei(i-nprim)
        base_label(1:LEN(base_label))=''
        base_label = "glf_eff_particle_chi_"//namei(i-nprim)(1:len_trim(namei(i-nprim)))//"_imp"

        CALL netcdf_err( nf90_def_var(ncid, base_label, nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_glf_eff_chi(i)))
        CALL netcdf_err( nf90_put_att(ncid,id_glf_eff_chi(i),'long_name',label ))
        CALL netcdf_err(  nf90_put_att(ncid,id_glf_eff_chi(i),'units',       &
          '1/(meter^2 second)'))
     ENDDO
 
     i=nion+1
     label = '* glf effective electron energy diffusivity'
     CALL netcdf_err( nf90_def_var(ncid, "glf_eff_elect_eng_chi", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_glf_eff_chi(i)))
     CALL netcdf_err( nf90_put_att(ncid,id_glf_eff_chi(i),'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_glf_eff_chi(i),'units',       &
          '1/(meter^2 second)'))

     i=nion+2
     label = '* glf effective ion energy diffusivity'
     CALL netcdf_err( nf90_def_var(ncid, "glf_eff_ion_eng_chi", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_glf_eff_chi(i)))
     CALL netcdf_err( nf90_put_att(ncid,id_glf_eff_chi(i),'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_glf_eff_chi(i),'units',       &
          '1/(meter^2 second)'))

     i =nion+3  ! no glf output for Faradays law, put something useful here
     label = '* glf mag field diffusivity, not set'
     CALL netcdf_err( nf90_def_var(ncid, "glf_eff_BP_chi", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_glf_eff_chi(i)))
     CALL netcdf_err( nf90_put_att(ncid,id_glf_eff_chi(i),'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_glf_eff_chi(i),'units',       &
          'not set '))

     i = nion+4
     label = '* glf effective momentum  diffusivity'
     CALL netcdf_err( nf90_def_var(ncid, "glf_eff_momtm_chi", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_glf_eff_chi(i)))
     CALL netcdf_err( nf90_put_att(ncid,id_glf_eff_chi(i),'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_glf_eff_chi(i),'units',       &
          '1/(meter^2 second)'))


     !total thermal ion  energy  flux:
     base_label = '* total thermal ion energy flux :  '
     CALL set_label(label,base_label,strlen,'nion')
     CALL netcdf_err( nf90_def_var(ncid, "e_fluxi", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_e_fluxi))
     CALL netcdf_err( nf90_put_att(ncid,id_e_fluxi,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_e_fluxi,'units',       &
          'joule/(meter^2 sec)'))

     !convective  thermal ion  energy  flux:
     base_label = '*convective thermal ion energy flux :  '
     CALL set_label(label,base_label,strlen,'nion')
     CALL netcdf_err( nf90_def_var(ncid, "e_fluxi_conv", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_e_fluxi_conv))
     CALL netcdf_err( nf90_put_att(ncid,id_e_fluxi_conv,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_e_fluxi_conv,'units',       &
          'joule/(meter^2 sec)'))



     !Flux associated with fdays law:
     base_label = '* Faradays law associated flux'
     CALL netcdf_err( nf90_def_var(ncid, "fday_flux", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_fdyflux))
     CALL netcdf_err( nf90_put_att(ncid,id_fdyflux,'long_name',base_label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_fdyflux,'units',       &
          'tesla/second'))

     !convective Flux associated with fdays law:
     base_label = '* Faradays law associated convective flux'
     CALL netcdf_err( nf90_def_var(ncid, "fday_flux_conv", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_fdyflux_conv))
     CALL netcdf_err( nf90_put_att(ncid,id_fdyflux_conv,'long_name',base_label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_fdyflux_conv,'units',       &
          'tesla/second'))


     !toroidal rotation flux:
     base_label = '* flux associated with toroidal rotation '
     CALL netcdf_err( nf90_def_var(ncid, "rot_flux", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_rotflux))
     CALL netcdf_err( nf90_put_att(ncid,id_rotflux,'long_name',base_label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_rotflux,'units',       &
          'kg/(second^2)'))

     !toroidal rotation convective flux:
     base_label = '* convective flux associated with toroidal rotation '
     CALL netcdf_err( nf90_def_var(ncid, "rot_flux_conv", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_rotflux_conv))
     CALL netcdf_err( nf90_put_att(ncid,id_rotflux_conv,'long_name',base_label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_rotflux_conv,'units',       &
          'kg/(second^2)'))

     !tglf turbulent particle  flux:
     base_label = '* tglf turbulent particle flux, species: '
     CALL set_label(label,base_label,strlen,'electron')
     CALL netcdf_err( nf90_def_var(ncid,"tglf_elct_p_flux", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_tglf_p_fluxe))
     CALL netcdf_err( nf90_put_att(ncid,id_tglf_p_fluxe,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_tglf_p_fluxe,'units',       &
          '1/(meter^2 second)'))

    !tglf turbulent particle  flux:
     base_label = '* tglf turbulent particle flux, species: '
     CALL set_label(label,base_label,strlen,'effective prim ion')
     CALL netcdf_err( nf90_def_var(ncid, "tglf_ion_p_flux", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_tglf_p_fluxp))
     CALL netcdf_err( nf90_put_att(ncid,id_tglf_p_fluxp,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_tglf_p_fluxp,'units',       &
          '1/(meter^2 second)'))

    !tglf turbulent particle  flux:
     base_label = '* tglf turbulent  particle flux, species: '
     CALL set_label(label,base_label,strlen,'effective impurity')
     CALL netcdf_err( nf90_def_var(ncid,"tglf_imp_p_flux", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_tglf_p_fluxi))
     CALL netcdf_err( nf90_put_att(ncid,id_tglf_p_fluxi,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_tglf_p_fluxi,'units',       &
          '1/(meter^2 second)'))

    !tglf turbulent energy  flux:
     base_label = '* tglf turbulent energy  flux, species: '
     CALL set_label(label,base_label,strlen,'electron')
     CALL netcdf_err( nf90_def_var(ncid, "tglf_elc_e_flux", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_tglf_e_fluxe))
     CALL netcdf_err( nf90_put_att(ncid,id_tglf_e_fluxe,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_tglf_e_fluxe,'units',       &
          'Jou/(meter^2 second)'))

    !tglf turbulent energy  flux:
     base_label = '* tglf turbulent energy  flux, species: '
     CALL set_label(label,base_label,strlen,'effective prim ion')
     CALL netcdf_err( nf90_def_var(ncid, "tglf_ion_e_flux", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_tglf_e_fluxp))
     CALL netcdf_err( nf90_put_att(ncid,id_tglf_e_fluxp,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_tglf_e_fluxp,'units',       &
          'Jou/(meter^2 second)'))

    !tglf turbulent energy  flux:
     base_label = '* tglf turbulent  energy  flux, species: '
     CALL set_label(label,base_label,strlen,'effective impurity')
     CALL netcdf_err( nf90_def_var(ncid,"tglf_imp_e_flux", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_tglf_e_fluxi))
     CALL netcdf_err( nf90_put_att(ncid,id_tglf_e_fluxi,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_tglf_e_fluxi,'units',       &
          'Jou/(meter^2 second)'))

   !tglf turbulent toroidal momentum  flux:
     base_label = '* tglf turbulent toroidal momentum  flux, species: '
     CALL set_label(label,base_label,strlen,'electron')
     CALL netcdf_err( nf90_def_var(ncid, "tglf_elc_m_flux", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_tglf_m_fluxe))
     CALL netcdf_err( nf90_put_att(ncid,id_tglf_m_fluxe,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_tglf_m_fluxe,'units',       &
          'kg/sec^2'))

    !tglf turbulent toroidal momentum flux:
     base_label = '* tglf turbulent toroidal momentum flux, species: '
     CALL set_label(label,base_label,strlen,'effective prim ion')
     CALL netcdf_err( nf90_def_var(ncid, "tglf_ion_m_flux", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_tglf_m_fluxp))
     CALL netcdf_err( nf90_put_att(ncid,id_tglf_m_fluxp,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_tglf_m_fluxp,'units',       &
          'kg/sec^2'))

    !tglf turbulent toroidal momentum   flux:
     base_label = '* tglf turbulent toroidal momentum  flux, species: '
     CALL set_label(label,base_label,strlen,'effective impurity')
     CALL netcdf_err( nf90_def_var(ncid,"tglf_imp_m_flux", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_tglf_m_fluxi))
     CALL netcdf_err( nf90_put_att(ncid,id_tglf_m_fluxi,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_tglf_m_fluxi,'units',       &
          'kg/sec^2)'))


     !electron density time deriv:
     CALL netcdf_err( nf90_def_var(ncid,"dnedt", nf90_double,          &
          DIMIDS = (/idim_rho/),VARID=id_dnedt))
     WRITE(label,FMT='(1pe16.8)')tGCNMs
     READ(label,FMT='(a)')st_asc_time
     st_asc_time = ADJUSTL(st_asc_time)
     label = '* rate of change of ene  1/(m**3 sec) at time ' &
          //st_asc_time(1:LEN_TRIM(st_asc_time))
     CALL netcdf_err(  nf90_put_att(ncid,id_dnedt,'long_name',label))
     CALL netcdf_err(  nf90_put_att(ncid,id_dnedt ,'units',        &
          '1/meter^3 second'))

!!$     !thermal ion densities time derivatives:
!!$     base_label = '* thermal ion dni/dt by  species: '
!!$     CALL set_label(label,base_label,strlen,'nion')
!!$     CALL netcdf_err( nf90_def_var(ncid, "dnidt", nf90_double,    &
!!$          DIMIDS = (/idim_rho,idim_nion/),VARID=id_dnidt))
!!$     CALL netcdf_err( nf90_put_att(ncid,id_dnidt,'long_name',label ))
!!$     CALL netcdf_err(  nf90_put_att(ncid,id_dnidt,'units',       &
!!$          '1./(meter^3 sec)'))

     !thermal ion densities boundary condititons:
     base_label = " bc profles at time:"//bc_asc_time(1:LEN_TRIM(bc_asc_time)) &
          //',species: '
     CALL set_label(label,base_label,strlen,'nion')
     CALL netcdf_err( nf90_def_var(ncid, "en_bc", nf90_double,    &
          DIMIDS = (/idim_rho,idim_nion/),VARID=id_en_bc))

     CALL netcdf_err( nf90_put_att(ncid,id_en_bc,'long_name',label))
     CALL netcdf_err(  nf90_put_att(ncid,id_en_bc,'units','1/meter^3'))

     !thermal ion flux(or mixed)  boundary condititons:
     base_label = " bc profiles at time:"//bc_asc_time(1:LEN_TRIM(bc_asc_time))&
          //',flux for species: '
     CALL set_label(label,base_label,strlen,'nion')
     CALL netcdf_err( nf90_def_var(ncid, "flux_bc", nf90_double,    &
          DIMIDS = (/idim_rho,idim_ntot/),VARID=id_flux_bc))

     CALL netcdf_err( nf90_put_att(ncid,id_flux_bc,'long_name',label))
     CALL netcdf_err( nf90_put_att(ncid,id_flux_bc,'units','keV/(meter^2 sec)'))

     !thermal ion source:
     base_label = '  sion : source due to ionization,'
     CALL set_label(label,base_label,strlen,'nion')
     CALL netcdf_err( nf90_def_var(ncid, "sion", nf90_double,    &
          DIMIDS = (/idim_rho,idim_nion/),VARID=id_sion))
     CALL netcdf_err( nf90_put_att(ncid,id_sion,'long_name',label))
     CALL netcdf_err(  nf90_put_att(ncid,id_sion,'units','1/(meter^3 sec)'))

     !thermal ion recombination sink:
     base_label = '  srecom : (-)source due to recombinatioin,'
     CALL set_label(label,base_label,strlen,'nion')
     CALL netcdf_err( nf90_def_var(ncid, "srecom", nf90_double,    &
          DIMIDS = (/idim_rho,idim_nion/),VARID=id_srecom))
     CALL netcdf_err( nf90_put_att(ncid,id_srecom,'long_name',label))
     CALL netcdf_err(  nf90_put_att(ncid,id_srecom,'units','1/(meter^3 sec)'))

     !source due to cx thermal neut
     base_label = '*  scx : source due to cx thermal neut.,'
     CALL set_label(label,base_label,strlen,'nion')
     CALL netcdf_err( nf90_def_var(ncid, "scx", nf90_double,    &
          DIMIDS = (/idim_rho,idim_nion/),VARID=id_scx))
     CALL netcdf_err( nf90_put_att(ncid,id_scx,'long_name',label))
     CALL netcdf_err( nf90_put_att(ncid,id_scx,'units','1/(meter^3 sec)'))

     base_label = '*  sbcx : sink due to cx with beam neut.,'
     CALL set_label(label,base_label,strlen,'nion')
     CALL netcdf_err( nf90_def_var(ncid, "sbcx", nf90_double,    &
          DIMIDS = (/idim_rho,idim_nion/),VARID=id_sbcx))
     CALL netcdf_err( nf90_put_att(ncid,id_sbcx,'long_name',label))
     CALL netcdf_err( nf90_put_att(ncid,id_sbcx,'units','1/(meter^3 sec)'))


     base_label = '*  stsource  : total source rate,'
     CALL set_label(label,base_label,strlen,'nprim')
     CALL netcdf_err( nf90_def_var(ncid, "stsource", nf90_double,    &
          DIMIDS = (/idim_nion,idim_rho/),VARID=id_stsource))
     CALL netcdf_err( nf90_put_att(ncid,id_stsource,'long_name',label))
     CALL netcdf_err( nf90_put_att(ncid,id_stsource,'units','1/(meter^3 sec)'))


     base_label = '*  dudt : s dot, (# or  energy) /(meter^3*second)'
     CALL set_label(label,base_label,strlen,'ntot')
     CALL netcdf_err( nf90_def_var(ncid, "dudtsv", nf90_double,    &
          DIMIDS = (/idim_ntot,idim_rho/),VARID=id_dudtsv))
     CALL netcdf_err( nf90_put_att(ncid,id_dudtsv,'long_name',label))
     CALL netcdf_err( nf90_put_att(ncid,id_dudtsv,'units','(1 or keV)/(meter^3 sec)'))

     !neutral  densities:
     IF ( nneu > 0 ) THEN 
         
        base_label = '*  neutral density, #/meter^3, species: '
        CALL set_label(label,base_label,strlen,'nneu')
        CALL netcdf_err( nf90_def_var(ncid, "enn", nf90_double,    &
             DIMIDS = (/idim_rho,idim_nneu/),VARID=id_enn))
        CALL netcdf_err( nf90_put_att(ncid,id_enn,'long_name',label))
        CALL netcdf_err( nf90_put_att(ncid,id_enn,'units',       &
             '#/meter^3'))

        base_label = '*  neutral density,due to wall source , species: '
        CALL set_label(label,base_label,strlen,'nneu')
        CALL netcdf_err( nf90_def_var(ncid, "ennw", nf90_double,    &
             DIMIDS = (/idim_rho,idim_nneu/),VARID=id_ennw))
        CALL netcdf_err( nf90_put_att(ncid,id_ennw,'long_name',label))
        CALL netcdf_err( nf90_put_att(ncid,id_ennw,'units',       &
             '#/meter^3'))

        base_label = '*  neutral density,due to volume  source , species: '
        CALL set_label(label,base_label,strlen,'nneu')
        CALL netcdf_err( nf90_def_var(ncid, "ennv", nf90_double,    &
             DIMIDS = (/idim_rho,idim_nneu/),VARID=id_ennv))
        CALL netcdf_err( nf90_put_att(ncid,id_ennv,'long_name',label))
        CALL netcdf_err( nf90_put_att(ncid,id_ennv,'units',       &
             '#/meter^3'))

        base_label = '* volume  source of neutrals , species: '
        CALL set_label(label,base_label,strlen,'nneu')
        CALL netcdf_err( nf90_def_var(ncid, "volsn", nf90_double,    &
             DIMIDS = (/idim_rho,idim_nneu/),VARID=id_volsn))
        CALL netcdf_err( nf90_put_att(ncid,id_volsn,'long_name',label  ))
        CALL netcdf_err( nf90_put_att(ncid,id_volsn,'units',       &
             '#/(meter^3*second)'))
     
     ENDIF !Neutrals
     
     !z and zsq:
     base_label = '* charge  z ,species : '
     CALL set_label(label,base_label,strlen,'nion')
     CALL netcdf_err( nf90_def_var(ncid, "z", nf90_double,    &
          DIMIDS = (/idim_rho,idim_nion/),VARID=id_z))
     CALL netcdf_err( nf90_put_att(ncid,id_z,'long_name',label))
     CALL netcdf_err( nf90_def_var(ncid, "zsq", nf90_double,    &
          DIMIDS = (/idim_rho,idim_nion/),VARID=id_zsq))
     CALL netcdf_err( nf90_put_att(ncid,id_zsq,'long_name',label))


     CALL netcdf_err( nf90_def_var(ncid, "stfuse", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_stfuse))
     CALL netcdf_err( nf90_put_att(ncid,id_stfuse,'long_name',    &
          '*  stfuse : thermal fusion rate, #/(meter^3*second)'  ))
     CALL netcdf_err( nf90_put_att(ncid,id_stfuse,'units',       &
          '#/(meter^3*second)'))


     CALL netcdf_err( nf90_def_var(ncid, "sbfuse", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_sbfuse))
     CALL netcdf_err( nf90_put_att(ncid,id_sbfuse,'long_name',    &
          '*  sbfuse : beam-thermal fusion rate, #/(meter^3*second)'  ))
     CALL netcdf_err( nf90_put_att(ncid,id_sbfuse,'units',       &
          '#/(meter^3*second)'))


     CALL netcdf_err( nf90_def_var(ncid, "spellet", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_spellet))
     CALL netcdf_err( nf90_put_att(ncid,id_spellet,'long_name',    &
          '*  spellet : THERMAL ion source due to pellets, #/(meter^3*second)'  ))
     CALL netcdf_err( nf90_put_att(ncid,id_spellet,'units',       &
          '#/(meter^3*second)'))




     CALL netcdf_err( nf90_def_var(ncid, "sbeame", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_sbeame))
     CALL netcdf_err( nf90_put_att(ncid,id_sbeame,'long_name',    &
          '*  sbeame : beam electron source,#/(meter^3*second)'  ))
     CALL netcdf_err( nf90_put_att(ncid,id_sbeame,'units',       &
          '#/(meter^3*second)'))

     base_label = '*  sbeam : beam ion  source,#/(meter^3*second), species: '
     CALL set_label(label,base_label,strlen,'nf')
     CALL netcdf_err( nf90_def_var(ncid, "sbeam", nf90_double,    &
          DIMIDS = (/idim_rho,idim_nbion/),VARID=id_sbeam))
     CALL netcdf_err( nf90_put_att(ncid,id_sbeam,'long_name',label))
     CALL netcdf_err( nf90_put_att(ncid,id_sbeam,'units',       &
          '#/(meter^3*second)'))

     base_label = '*  fast ion density, #/meter^3, species: '
     CALL set_label(label,base_label,strlen,'nf')
     CALL netcdf_err( nf90_def_var(ncid, "enbeam", nf90_double,    &
          DIMIDS = (/idim_rho,idim_nbion/),VARID=id_enbeam))
     CALL netcdf_err( nf90_put_att(ncid,id_enbeam,'long_name',label  ))
     CALL netcdf_err( nf90_put_att(ncid,id_enbeam,'units',       &
          '#/(meter^3)'))



     CALL netcdf_err( nf90_def_var(ncid, "curden", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_curden))
     CALL netcdf_err(  nf90_put_att(ncid,id_curden,'long_name',    &
          '*  total toroidal current density, <Jtor R0/R>, amp/meter**2'))
     CALL netcdf_err(  nf90_put_att(ncid,id_curden,'units','amp/meter^2)'))

    CALL netcdf_err( nf90_def_var(ncid, "curpar", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_curpar))
     CALL netcdf_err(  nf90_put_att(ncid,id_curpar,'long_name',    &
          '*  total parallel current density, <J.B/Bt0>, amp/meter**2'))
     CALL netcdf_err(  nf90_put_att(ncid,id_curpar,'units','amp/meter^2)'))

     CALL netcdf_err( nf90_def_var(ncid, "curohm", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_curohm))
     CALL netcdf_err(  nf90_put_att(ncid,id_curohm,'long_name',    &
          '* ohmic current density, amp/meter**2'))
     CALL netcdf_err(  nf90_put_att(ncid,id_curohm,'units',        &
          'amp/meter^2)'))

     CALL netcdf_err( nf90_def_var(ncid, "curboot", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_curboot))
     CALL netcdf_err(  nf90_put_att(ncid,id_curboot,'long_name',    &
          '* bootstrap current density, amp/meter**2'))
     CALL netcdf_err(  nf90_put_att(ncid,id_curboot,'units',        &
          'amp/meter^2)'))
 
     CALL netcdf_err( nf90_def_var(ncid, "curbeam", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_curbeam))
     CALL netcdf_err(  nf90_put_att(ncid,id_curbeam,'long_name',    &
          '* beam driven  current density, amp/meter**2'))
     CALL netcdf_err(  nf90_put_att(ncid,id_curbeam,'units',        &
          'amp/meter^2)'))


     CALL netcdf_err( nf90_def_var(ncid, "currf", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_currf))
     CALL netcdf_err(  nf90_put_att(ncid,id_currf,'long_name',    &
          '* rf driven  current density, amp/meter**2'))
     CALL netcdf_err(  nf90_put_att(ncid,id_currf,'units',        &
          'amp/meter^2)'))


     CALL netcdf_err( nf90_def_var(ncid, "etor", nf90_double,     &
          DIMIDS = (/idim_rho/),VARID=id_etor))
     CALL netcdf_err(  nf90_put_att(ncid,id_etor,'long_name',    &
          '*  toroidal electric field profile, V/m'))
     CALL netcdf_err(  nf90_put_att(ncid,id_etor,'units',        &
          'volt/meter)'))



     CALL netcdf_err( nf90_def_var(ncid, "rbp", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_rbp))
     CALL netcdf_err(  nf90_put_att(ncid,id_rbp,'long_name',    &
          '*  rho*bp0*fcap*gcap*hcap, tesla*meter'))
     CALL netcdf_err(  nf90_put_att(ncid,id_rbp,'units',        &
          'tesla*meter)'))

     CALL netcdf_err( nf90_def_var(ncid, "psivalnpsi", nf90_double,    &
          DIMIDS = (/idim_npsi/),VARID=id_psivalnpsi))
     CALL netcdf_err(  nf90_put_att(ncid,id_psivalnpsi,'long_name',    &
          '* psivalnpsi(npsi) grid, edge to magnetic axis, volt sec/rad'))
     CALL netcdf_err(  nf90_put_att(ncid,id_psivalnpsi,'units', &
          'volt second/radian'))


     CALL netcdf_err( nf90_def_var(ncid, "ravgnpsi", nf90_double,    &
          DIMIDS = (/idim_npsi/),VARID=id_ravgnpsi))
     CALL netcdf_err(  nf90_put_att(ncid,id_ravgnpsi,'long_name',    &
          '* ravg: <R> avg radius on mhd grid, m'))
     CALL netcdf_err(  nf90_put_att(ncid,id_ravgnpsi,'units','meter'))

     CALL netcdf_err( nf90_def_var(ncid, "ravginpsi", nf90_double,    &
          DIMIDS = (/idim_npsi/),VARID=id_ravginpsi))
     CALL netcdf_err(  nf90_put_att(ncid,id_ravginpsi,'long_name',    &
          '* ravgi:<1/R>  on mhd grid, 1/m'))
     CALL netcdf_err(  nf90_put_att(ncid,id_ravginpsi,'units','meter'))


     CALL netcdf_err( nf90_def_var(ncid, "fpsinpsi", nf90_double,    &
          DIMIDS = (/idim_npsi/),VARID=id_fpsinpsi))
     CALL netcdf_err(  nf90_put_att(ncid,id_fpsinpsi,'long_name',    &
          '* fpsi: f of psi (= R*Bt) on psivalnpsi(npsi) grid, tesla meter'))
     CALL netcdf_err(  nf90_put_att(ncid,id_fpsinpsi,'units','tesla*meter'))
 
     CALL netcdf_err( nf90_def_var(ncid, "pprim", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_pprim))
     CALL netcdf_err(  nf90_put_att(ncid,id_pprim,'long_name',    &
          '* pprim: dp/dpsi on transport grid, nt/(m**2-volt-sec) = amp/m**3'))
     CALL netcdf_err(  nf90_put_att(ncid,id_pprim,'units',&
          'newton/(meter^2*volt*sec'))


     CALL netcdf_err( nf90_def_var(ncid, "ffprim", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_ffprim))
     CALL netcdf_err(  nf90_put_att(ncid,id_ffprim,'long_name',    &
          '* f*df/dpsi: on transport grid, kg/(A sec^2)'))
     CALL netcdf_err(  nf90_put_att(ncid,id_ffprim,'units', 'kg/(amp sec^2)'))


     CALL netcdf_err( nf90_def_var(ncid, "bp", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_bp))
     CALL netcdf_err(  nf90_put_att(ncid,id_bp,'long_name',    &
          '*<Bp> : flux avg B poloidal on transport grid, Tesla'))
     CALL netcdf_err(  nf90_put_att(ncid,id_bp,'units',        &
          'tesla'))

     CALL netcdf_err( nf90_def_var(ncid, "bprmaj", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_bprmaj))
     CALL netcdf_err(  nf90_put_att(ncid,id_bprmaj,'long_name',    &
          '* B poloidal on rmaj (and transport rho)  grid, tesla'))
     CALL netcdf_err(  nf90_put_att(ncid,id_bprmaj,'units',        &
          'tesla)'))


     CALL netcdf_err( nf90_def_var(ncid, "btotrmaj", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_btotrmaj))
     CALL netcdf_err(  nf90_put_att(ncid,id_btotrmaj,'long_name',    &
          '* Btotal on transport rho grid, tesla'))
     CALL netcdf_err(  nf90_put_att(ncid,id_btotrmaj,'units',        &
          'tesla)'))



     CALL netcdf_err( nf90_def_var(ncid, "zeff", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_zeff))
     CALL netcdf_err(  nf90_put_att(ncid,id_zeff,'long_name',    &
          '*  zeff profile'))
     CALL netcdf_err(  nf90_put_att(ncid,id_zeff,'units', dimensionless  ))



     CALL netcdf_err( nf90_def_var(ncid, "vpol", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_vpol))
     CALL netcdf_err(  nf90_put_att(ncid,id_vpol,'long_name',    &
          '*  poloidal velocity profile'))
     CALL netcdf_err(  nf90_put_att(ncid,id_vpol,'units', 'm/sec'  ))

     CALL netcdf_err( nf90_def_var(ncid, "vpol_nclass", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_vpol_nclass))
     CALL netcdf_err(  nf90_put_att(ncid,id_vpol_nclass,'long_name',    &
          '*  poloidal velocity profile, forcebal Nclass model'))
     CALL netcdf_err(  nf90_put_att(ncid,id_vpol_nclass,'units', 'm/sec'  ))




     CALL netcdf_err( nf90_def_var(ncid, "vpar", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_vpar))
     CALL netcdf_err(  nf90_put_att(ncid,id_vpar,'long_name',    &
          '*  parallel velocity profile'))
     CALL netcdf_err(  nf90_put_att(ncid,id_vpar,'units', 'm/sec'  ))

     CALL netcdf_err( nf90_def_var(ncid, "vpar_nclass", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_vpar_nclass))
     CALL netcdf_err(  nf90_put_att(ncid,id_vpar_nclass,'long_name',    &
          '*  parallel velocity profile, forcebal Nclass model'))
     CALL netcdf_err(  nf90_put_att(ncid,id_vpar_nclass,'units', 'm/sec'  ))

     CALL netcdf_err( nf90_def_var(ncid, "er_tot_nclass", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_er_tot_nclass))
     CALL netcdf_err(  nf90_put_att(ncid,id_er_tot_nclass,'long_name',    &
          '* Nclass total radial electric field '))
     CALL netcdf_err(  nf90_put_att(ncid,id_er_tot_nclass,'units', 'v/m'  ))




     CALL netcdf_err( nf90_def_var(ncid, "angrot", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_angrot))
     CALL netcdf_err(  nf90_put_att(ncid,id_angrot,'long_name',    &
          '*  angular rotation speed profile, rad/sec'))
     CALL netcdf_err(  nf90_put_att(ncid,id_angrot,'units',        &
          'radian/second'))

     !new 8/06/07:
     CALL netcdf_err( nf90_def_var(ncid, "d", nf90_double,    &
          DIMIDS = (/idim_ntot,idim_ntot,idim_rho/),VARID=id_d))
     CALL netcdf_err(  nf90_put_att(ncid,id_d,'long_name',    &
          '*  diffusion matrix(ntot,ntot,nj) , on half grid'))
     CALL netcdf_err(  nf90_put_att(ncid,id_d,'units',        &
          'mixed'))


     CALL netcdf_err( nf90_def_var(ncid, "chieinv", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_chieinv))
     CALL netcdf_err(  nf90_put_att(ncid,id_chieinv,'long_name',    &
          '*  electron thermal diffusivity, meter**2/sec, on half grid'))
     CALL netcdf_err(  nf90_put_att(ncid,id_chieinv,'units',        &
          'meter^2/sec'))

     CALL netcdf_err( nf90_def_var(ncid, "chiinv", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_chiinv))
     CALL netcdf_err(  nf90_put_att(ncid,id_chiinv,'long_name',    &
          '*  ion thermal diffusivity, meter**2/sec, on half grid'))
     CALL netcdf_err(  nf90_put_att(ncid,id_chiinv,'units',        &
          'meter^2/sec'))

     CALL netcdf_err( nf90_def_var(ncid, "xkineo", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_xkineo))
     CALL netcdf_err(  nf90_put_att(ncid,id_xkineo,'long_name',    &
          '*  ion neoclassical thermal conductivity'))
     CALL netcdf_err(  nf90_put_att(ncid,id_xkineo,'units',        &
          '1/(meter*second'))

     CALL netcdf_err( nf90_def_var(ncid, "xkeneo", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_xkeneo))
     CALL netcdf_err(  nf90_put_att(ncid,id_xkeneo,'long_name',    &
          '*  electron neoclassical thermal conductivity'))
     CALL netcdf_err(  nf90_put_att(ncid,id_xkeneo,'units',        &
          '1/(meter*second'))

     CALL netcdf_err( nf90_def_var(ncid, "dpedt", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_dpedt))
     CALL netcdf_err(  nf90_put_att(ncid,id_dpedt,'long_name',    &
          '*  1.5*dpedt, power density due to electron pressure, watts/meter^3'))
     CALL netcdf_err(  nf90_put_att(ncid,id_dpedt,'units',        &
          'watts/meter^3'))


     CALL netcdf_err( nf90_def_var(ncid, "dpidt", nf90_double,    &
          DIMIDS = (/idim_rho,idim_nion/),VARID=id_dpidt))
     CALL netcdf_err(  nf90_put_att(ncid,id_dpidt,'long_name',    &
          '* 1.5*dpidt,power density due to ion pressure, watts/meter^3'))
     CALL netcdf_err(  nf90_put_att(ncid,id_dpidt,'units',        &
          'watts/meter^3'))

     CALL netcdf_err( nf90_def_var(ncid, "qconde", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_qconde))
     CALL netcdf_err(  nf90_put_att(ncid,id_qconde,'long_name',    &
          '*  electron conduction, watts/meter^3'))
     CALL netcdf_err(  nf90_put_att(ncid,id_qconde,'units',        &
          'watts/meter^3'))

     CALL netcdf_err( nf90_def_var(ncid, "qcondi", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_qcondi))
     CALL netcdf_err(  nf90_put_att(ncid,id_qcondi,'long_name',    &
          '*  iom conduction, watts/meter^3'))
     CALL netcdf_err(  nf90_put_att(ncid,id_qcondi,'units',        &
          'watts/meter^3'))

     CALL netcdf_err( nf90_def_var(ncid, "qconve", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_qconve))
     CALL netcdf_err(  nf90_put_att(ncid,id_qconve,'long_name',    &
          '*  electron convection, watts/meter^3'))
     CALL netcdf_err(  nf90_put_att(ncid,id_qconve,'units',        &
          'watts/meter^3'))

     CALL netcdf_err( nf90_def_var(ncid, "qconvi", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_qconvi))
     CALL netcdf_err(  nf90_put_att(ncid,id_qconvi,'long_name',    &
          '*  ion convection, watts/meter^3'))
     CALL netcdf_err(  nf90_put_att(ncid,id_qconvi,'units',        &
          'watts/meter^3'))


     CALL netcdf_err( nf90_def_var(ncid, "qbeame", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_qbeame))
     CALL netcdf_err(  nf90_put_att(ncid,id_qbeame,'long_name',    &
          '*  power to elec. from beam, watts/meter^3'))
     CALL netcdf_err(  nf90_put_att(ncid,id_qbeame,'units',        &
          'watts/meter^3'))

     CALL netcdf_err( nf90_def_var(ncid, "qbeami", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_qbeami))
     CALL netcdf_err(  nf90_put_att(ncid,id_qbeami,'long_name',    &
          '*  power to ions  from beam, watts/meter^3'))
     CALL netcdf_err(  nf90_put_att(ncid,id_qbeami,'units',        &
          'watts/meter^3'))

     CALL netcdf_err( nf90_def_var(ncid, "qdelt", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_qdelt))
     CALL netcdf_err(  nf90_put_att(ncid,id_qdelt,'long_name',    &
          '*  electron ion energy exchange term, watts/meter^3'))
     CALL netcdf_err(  nf90_put_att(ncid,id_qdelt,'units',        &
          'watts/meter^3'))

     CALL netcdf_err( nf90_def_var(ncid, "qexch", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_qexch))
     CALL netcdf_err(  nf90_put_att(ncid,id_qexch,'long_name',    &
          '* qexch, anomalous electron-ion energy exchange term, watts/meter^3'))
     CALL netcdf_err(  nf90_put_att(ncid,id_qexch,'units',        &
          'watts/meter^3'))


     CALL netcdf_err( nf90_def_var(ncid, "qrfe", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_qrfe))
     CALL netcdf_err(  nf90_put_att(ncid,id_qrfe,'long_name',    &
          '*  qrfe, RF electron heating, watts/meter^3'))
     CALL netcdf_err(  nf90_put_att(ncid,id_qrfe,'units',        &
          'watts/meter^3'))
 
     CALL netcdf_err( nf90_def_var(ncid, "qrfi", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_qrfi))
     CALL netcdf_err(  nf90_put_att(ncid,id_qrfi,'long_name',    &
          '*  qrfi, RF electron heating, watts/meter^3'))
     CALL netcdf_err(  nf90_put_att(ncid,id_qrfi,'units',        &
          'watts/meter^3'))


     CALL netcdf_err( nf90_def_var(ncid, "qione", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_qione))
     CALL netcdf_err(  nf90_put_att(ncid,id_qione,'long_name',    &
          '*  qione, electron power density due to  recombination and impact ionization, watts/meter^3'))
     CALL netcdf_err(  nf90_put_att(ncid,id_qione,'units',        &
          'watts/meter^3'))


     CALL netcdf_err( nf90_def_var(ncid, "qioni", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_qioni))
     CALL netcdf_err(  nf90_put_att(ncid,id_qioni,'long_name',    &
          '*  qioni, ion power density due to  recombination and impact ionization, watts/meter^3'))
     CALL netcdf_err(  nf90_put_att(ncid,id_qioni,'units',        &
          'watts/meter^3'))


     CALL netcdf_err( nf90_def_var(ncid, "qcx", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_qcx))
     CALL netcdf_err(  nf90_put_att(ncid,id_qcx,'long_name',    &
          '*  qcx, ion power density due to neutral-ion charge exchange, watts/meter^3'))
     CALL netcdf_err(  nf90_put_att(ncid,id_qcx,'units',        &
          'watts/meter^3'))



     CALL netcdf_err( nf90_def_var(ncid, "qe2d", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_qe2d))
     CALL netcdf_err(  nf90_put_att(ncid,id_qe2d,'long_name',    &
          '*  2d electron heating, watts/meter^3'))
     CALL netcdf_err(  nf90_put_att(ncid,id_qe2d,'units',        &
          'watts/meter^3'))



     CALL netcdf_err( nf90_def_var(ncid, "qi2d", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_qi2d))
     CALL netcdf_err(  nf90_put_att(ncid,id_qi2d,'long_name',    &
          '*  2d ion heating, watts/meter^3'))
     CALL netcdf_err(  nf90_put_att(ncid,id_qi2d,'units',        &
          'watts/meter^3'))



     CALL netcdf_err( nf90_def_var(ncid, "qfuse", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_qfuse))
     CALL netcdf_err(  nf90_put_att(ncid,id_qfuse,'long_name',    &
          '* total fusion electron heating, watts/meter^3'))
     CALL netcdf_err(  nf90_put_att(ncid,id_qfuse,'units',        &
          'watts/meter^3'))


     CALL netcdf_err( nf90_def_var(ncid, "qfusi", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_qfusi))
     CALL netcdf_err(  nf90_put_att(ncid,id_qfusi,'long_name',    &
          '* total fusion ion heating, watts/meter^3'))
     CALL netcdf_err(  nf90_put_att(ncid,id_qfusi,'units',        &
          'watts/meter^3'))



     CALL netcdf_err( nf90_def_var(ncid, "qbfuse", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_qbfuse))
     CALL netcdf_err(  nf90_put_att(ncid,id_qbfuse,'long_name',    &
          '*  beam fusion electron heating, watts/meter^3'))
     CALL netcdf_err(  nf90_put_att(ncid,id_qbfuse,'units',        &
          'watts/meter^3'))



     CALL netcdf_err( nf90_def_var(ncid, "qbfusi", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_qbfusi))
     CALL netcdf_err(  nf90_put_att(ncid,id_qbfusi,'long_name',    &
          '* beam fusion ion heating, watts/meter^3'))
     CALL netcdf_err(  nf90_put_att(ncid,id_qbfusi,'units',        &
          'watts/meter^3'))



     CALL netcdf_err( nf90_def_var(ncid,"qmag", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_qmag))
     CALL netcdf_err(  nf90_put_att(ncid,id_qmag,'long_name',    &
          '*  qmag electron heating, watts/meter^3'))
     CALL netcdf_err(  nf90_put_att(ncid,id_qmag,'units',        &
          'watts/meter^3'))

     CALL netcdf_err( nf90_def_var(ncid,"qsawe", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_qsawe))
     CALL netcdf_err(  nf90_put_att(ncid,id_qsawe,'long_name',    &
          '*  sawtooth electron heating, watts/meter^3'))
     CALL netcdf_err(  nf90_put_att(ncid,id_qsawe,'units',        &
          'watts/meter^3'))

     CALL netcdf_err( nf90_def_var(ncid,"qsawi", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_qsawi))
     CALL netcdf_err(  nf90_put_att(ncid,id_qsawi,'long_name',    &
          '*  sawtooth ion heating, watts/meter^3'))
     CALL netcdf_err(  nf90_put_att(ncid,id_qsawi,'units',        &
          'watts/meter^3'))


     CALL netcdf_err( nf90_def_var(ncid,"qrad", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_qrad))
     CALL netcdf_err(  nf90_put_att(ncid,id_qrad,'long_name',    &
          '*  radiated power density, watts/meter^3'))
     CALL netcdf_err(  nf90_put_att(ncid,id_qrad,'units',        &
          'watts/meter^3'))

    ! Nclass derived pinch velocities,ions and electrons
     base_label = '* nclass derived pinch velocity: '
     CALL set_label(label,base_label,strlen,'vpinch Nclass')
     CALL netcdf_err( nf90_def_var(ncid, "vpinch_nclass", nf90_double,    &
          DIMIDS = (/idim_rho,idim_nionp1/),VARID=id_vpinch_nions))
     CALL netcdf_err( nf90_put_att(ncid,id_vpinch_nions,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_vpinch_nions,'units',       &
          'm/s'))

    !-----------------------------------------------------------------------
    ! multimode related inputs/outputs:
    !-----------------------------------------------------------------------
     base_label = '* multimode - growth rate most unstable DRIB mode [1/sec]: '
     CALL set_label(label,base_label,strlen,'mmm_gammaDBM')
     CALL netcdf_err( nf90_def_var(ncid, "mmm_gammaDBM", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_mmm_gammaDBM))
     CALL netcdf_err( nf90_put_att(ncid,id_mmm_gammaDBM,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_mmm_gammaDBM,'units',       &
          '1/s'))

     base_label = '* multimode - freq  most unstable DRIB mode [rad/sec]: '
     CALL set_label(label,base_label,strlen,'mmm_omegaDBM')
     CALL netcdf_err( nf90_def_var(ncid, "mmm_omegaDBM", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_mmm_omegaDBM))
     CALL netcdf_err( nf90_put_att(ncid,id_mmm_omegaDBM,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_mmm_omegaDBM,'units',       &
          'rad/s'))

     base_label = '* multimode - eff ion hydrogenic  diffusivity,M^2/sec: '
     CALL set_label(label,base_label,strlen,'mmm_xdi')
     CALL netcdf_err( nf90_def_var(ncid, "mmm_xdi", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_mmm_xdi))
     CALL netcdf_err( nf90_put_att(ncid,id_mmm_xdi,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_mmm_xdi,'units',       &
          'M^2/sec'))

     base_label = '* multimode - eff ion thermal diffusivity,M^2/sec: '
     CALL set_label(label,base_label,strlen,'mmm_xti')
     CALL netcdf_err( nf90_def_var(ncid, "mmm_xti", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_mmm_xti))
     CALL netcdf_err( nf90_put_att(ncid,id_mmm_xti,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_mmm_xti,'units',       &
          'M^2/sec'))

     base_label = '* multimode - eff electron thermal diffusivity,M^2/sec '
     CALL set_label(label,base_label,strlen,'mmm_xte')
     CALL netcdf_err( nf90_def_var(ncid, "mmm_xte", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_mmm_xte))
     CALL netcdf_err( nf90_put_att(ncid,id_mmm_xte,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_mmm_xte,'units',       &
          'M^2/sec'))

     base_label = '* multimode - impurity ion  diffusivity Weiland model,M^2/sec '
     CALL set_label(label,base_label,strlen,'mmm_xdz')
     CALL netcdf_err( nf90_def_var(ncid, "mmm_xdz", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_mmm_xdz))
     CALL netcdf_err( nf90_put_att(ncid,id_mmm_xdz,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_mmm_xdz,'units',       &
          'M^2/sec'))
 
     base_label = '* multimode - toroidal momentum transport, Weiland model  diffusivity,M^2/sec'
     CALL set_label(label,base_label,strlen,'mmm_xvt')
     CALL netcdf_err( nf90_def_var(ncid, "mmm_xvt", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_mmm_xvt))
     CALL netcdf_err( nf90_put_att(ncid,id_mmm_xvt,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_mmm_xvt,'units',       &
          'M^2/sec'))

     base_label = '* multimode - poloidal  momentum transport, Weiland model  diffusivity,M^2/sec'
     CALL set_label(label,base_label,strlen,'mmm_xvp')
     CALL netcdf_err( nf90_def_var(ncid, "mmm_xvp", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_mmm_xvp))
     CALL netcdf_err( nf90_put_att(ncid,id_mmm_xvp,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_mmm_xvp,'units',       &
          'M^2/sec'))

     base_label = '* multimode - ion thermal diffusivity,M^2/sec, Weiland model part'
     CALL set_label(label,base_label,strlen,'mmm_xtiW20')
     CALL netcdf_err( nf90_def_var(ncid, "mmm_xtiW20", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_mmm_xtiW20))
     CALL netcdf_err( nf90_put_att(ncid,id_mmm_xtiW20,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_mmm_xtiW20,'units',       &
          'M^2/sec'))

     base_label = '* multimode - particle  diffusivity,M^2/sec, Weiland model part'
     CALL set_label(label,base_label,strlen,'mmm_xdiW20')
     CALL netcdf_err( nf90_def_var(ncid, "mmm_xdiW20", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_mmm_xdiW20))
     CALL netcdf_err( nf90_put_att(ncid,id_mmm_xdiW20,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_mmm_xdiW20,'units',       &
          'M^2/sec'))

     base_label = '* multimode - electron  thermal diffusivity,M^2/sec, Weiland model part'
     CALL set_label(label,base_label,strlen,'mmm_xteW20')
     CALL netcdf_err( nf90_def_var(ncid, "mmm_xteW20", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_mmm_xteW20))
     CALL netcdf_err( nf90_put_att(ncid,id_mmm_xteW20,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_mmm_xteW20,'units',       &
          'M^2/sec'))

     base_label = '* multimode - ion thermal diffusivity,M^2/sec, DRB  model part'
     CALL set_label(label,base_label,strlen,'mmm_xtiDBM')
     CALL netcdf_err( nf90_def_var(ncid, "mmm_xtiDBM", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_mmm_xtiDBM))
     CALL netcdf_err( nf90_put_att(ncid,id_mmm_xtiDBM,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_mmm_xtiDBM,'units',       &
          'M^2/sec'))

     base_label = '* multimode - particle  diffusivity,M^2/sec, DRB model part'
     CALL set_label(label,base_label,strlen,'mmm_xdiDBM')
     CALL netcdf_err( nf90_def_var(ncid, "mmm_xdiDBM", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_mmm_xdiDBM))
     CALL netcdf_err( nf90_put_att(ncid,id_mmm_xdiDBM,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_mmm_xdiDBM,'units',       &
          'M^2/sec'))

     base_label = '* multimode - electron  thermal diffusivity,M^2/sec, DRBmodel part'
     CALL set_label(label,base_label,strlen,'mmm_xteDBM')
     CALL netcdf_err( nf90_def_var(ncid, "mmm_xteDBM", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_mmm_xteDBM))
     CALL netcdf_err( nf90_put_att(ncid,id_mmm_xteDBM,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_mmm_xteDBM,'units',       &
          'M^2/sec'))

     base_label = '* multimode - electron  thermal diffusivity,M^2/sec, ETG model part'
     CALL set_label(label,base_label,strlen,'mmm_xteETG')
     CALL netcdf_err( nf90_def_var(ncid, "mmm_xteETG", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_mmm_xteETG))
     CALL netcdf_err( nf90_put_att(ncid,id_mmm_xteETG,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_mmm_xteETG,'units',       &
          'M^2/sec'))


     base_label = '* multimode - growth rate most unstable ion mode Weiland positive freq direction [1/sec]'
     CALL set_label(label,base_label,strlen,'mmm_gamma_i1_W20')
     CALL netcdf_err( nf90_def_var(ncid, "mmm_gamma_i1_W20", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_mmm_gamma_i1_W20))
     CALL netcdf_err( nf90_put_att(ncid,id_mmm_gamma_i1_W20,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_mmm_gamma_i1_W20,'units',       &
          '1/sec'))

     base_label = '* multimode - growth rate most unstable elec mode Weiland positive freq direction [1/sec]'
     CALL set_label(label,base_label,strlen,'mmm_gamma_e1_W20')
     CALL netcdf_err( nf90_def_var(ncid, "mmm_gamma_e1_W20", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_mmm_gamma_e1_W20))
     CALL netcdf_err( nf90_put_att(ncid,id_mmm_gamma_e1_W20,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_mmm_gamma_e1_W20,'units',       &
          '1/sec'))

     base_label = '* multimode - growth rate most unstable ion mode Weiland negative freq direction [1/sec]'
     CALL set_label(label,base_label,strlen,'mmm_gamma_i2_W20')
     CALL netcdf_err( nf90_def_var(ncid, "mmm_gamma_i2_W20", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_mmm_gamma_i2_W20))
     CALL netcdf_err( nf90_put_att(ncid,id_mmm_gamma_i2_W20,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_mmm_gamma_i2_W20,'units',       &
          '1/sec'))

     base_label = '* multimode - growth rate most unstable elec mode Weiland negative freq direction [1/sec]'
     CALL set_label(label,base_label,strlen,'mmm_gamma_e2_W20')
     CALL netcdf_err( nf90_def_var(ncid, "mmm_gamma_e2_W20", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_mmm_gamma_e2_W20))
     CALL netcdf_err( nf90_put_att(ncid,id_mmm_gamma_e2_W20,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_mmm_gamma_e2_W20,'units',       &
          '1/sec'))

     base_label = '* multimode - freq,rad/sec, most unstable ion mode Weiland positive freq direction'
     CALL set_label(label,base_label,strlen,'mmm_omega_i1_W20')
     CALL netcdf_err( nf90_def_var(ncid, "mmm_omega_i1_W20", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_mmm_omega_i1_W20))
     CALL netcdf_err( nf90_put_att(ncid,id_mmm_omega_i1_W20,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_mmm_omega_i1_W20,'units',       &
          'rad/sec'))

     base_label = '* multimode - freq,rad/sec, most unstable elec mode Weiland positive freq direction'
     CALL set_label(label,base_label,strlen,'mmm_omega_e1_W20')
     CALL netcdf_err( nf90_def_var(ncid, "mmm_omega_e1_W20", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_mmm_omega_e1_W20))
     CALL netcdf_err( nf90_put_att(ncid,id_mmm_omega_e1_W20,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_mmm_omega_e1_W20,'units',       &
          'rad/sec'))

     base_label = '* multimode - freq,rad/sec, most unstable ion mode Weiland negative freq direction'
     CALL set_label(label,base_label,strlen,'mmm_omega_i2_W20')
     CALL netcdf_err( nf90_def_var(ncid, "mmm_omega_i2_W20", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_mmm_omega_i2_W20))
     CALL netcdf_err( nf90_put_att(ncid,id_mmm_omega_i2_W20,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_mmm_omega_i2_W20,'units',       &
          'rad/sec'))

     base_label = '* multimode - freq,rad/sec,most unstable elec mode Weiland negative freq direction'
     CALL set_label(label,base_label,strlen,'mmm_omega_e2_W20')
     CALL netcdf_err( nf90_def_var(ncid, "mmm_omega_e2_W20", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_mmm_omega_e2_W20))
     CALL netcdf_err( nf90_put_att(ncid,id_mmm_omega_e2_W20,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_mmm_omega_e2_W20,'units',       &
          'rad/sec'))

     base_label = '* multimode - total ion thermal flux,W/m^2'
     CALL set_label(label,base_label,strlen,'mmm_flux_ith')
     CALL netcdf_err( nf90_def_var(ncid, "mmm_flux_ith", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_mmm_flux_ith))
     CALL netcdf_err( nf90_put_att(ncid,id_mmm_flux_ith,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_mmm_flux_ith,'units',       &
          'W/m^2'))

     base_label = '* multimode - total hydrogenic ion flux flux,1/(m^2 sec)'
     CALL set_label(label,base_label,strlen,'mmm_flux_ip')
     CALL netcdf_err( nf90_def_var(ncid, "mmm_flux_ip", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_mmm_flux_ip))
     CALL netcdf_err( nf90_put_att(ncid,id_mmm_flux_ip,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_mmm_flux_ip,'units',       &
          '1/(m^2 sec'))

     base_label = '* multimode - total electron thermal  flux flux,W/m^2'
     CALL set_label(label,base_label,strlen,'mmm_flux_eth')
     CALL netcdf_err( nf90_def_var(ncid, "mmm_flux_eth", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_mmm_flux_eth))
     CALL netcdf_err( nf90_put_att(ncid,id_mmm_flux_eth,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_mmm_flux_eth,'units',       &
          'W/m^2'))

    base_label = '* multimode - total impurity ion flux 1/(m^2 sec)'
     CALL set_label(label,base_label,strlen,'mmm_flux_imp')
     CALL netcdf_err( nf90_def_var(ncid, "mmm_flux_imp", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_mmm_flux_imp))
     CALL netcdf_err( nf90_put_att(ncid,id_mmm_flux_imp,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_mmm_flux_imp,'units',       &
          '1/(m^2 sec'))

     base_label = '* multimode - ion thermal convective velocity,m/sec'
     CALL set_label(label,base_label,strlen,'mmm_vconv_ith')
     CALL netcdf_err( nf90_def_var(ncid, "mmm_vconv_ith", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_mmm_vconv_ith))
     CALL netcdf_err( nf90_put_att(ncid,id_mmm_vconv_ith,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_mmm_vconv_ith,'units',       &
          'm/sec'))

     base_label = '* multimode - hydrogenic ion particle  convective velocity,m/sec'
     CALL set_label(label,base_label,strlen,'mmm_vconv_ip')
     CALL netcdf_err( nf90_def_var(ncid, "mmm_vconv_ip", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_mmm_vconv_ip))
     CALL netcdf_err( nf90_put_att(ncid,id_mmm_vconv_ip,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_mmm_vconv_ip,'units',       &
          'm/sec'))

     base_label = '* multimode - electron thermal convective velocity,m/sec'
     CALL set_label(label,base_label,strlen,'mmm_vconv_eth')
     CALL netcdf_err( nf90_def_var(ncid, "mmm_vconv_eth", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_mmm_vconv_eth))
     CALL netcdf_err( nf90_put_att(ncid,id_mmm_vconv_eth,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_mmm_vconv_eth,'units',       &
          'm/sec'))

     base_label = '* multimode - impurity ion particle  convective velocity,m/sec'
     CALL set_label(label,base_label,strlen,'mmm_vconv_imp')
     CALL netcdf_err( nf90_def_var(ncid, "mmm_vconv_imp", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_mmm_vconv_imp))
     CALL netcdf_err( nf90_put_att(ncid,id_mmm_vconv_imp,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_mmm_vconv_imp,'units',       &
          'm/sec'))



     base_label = '* multimode - toroidal momentum pinch ,m/sec'
     CALL set_label(label,base_label,strlen,'mmm_vmtmt')
     CALL netcdf_err( nf90_def_var(ncid, "mmm_vmtmt", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_mmm_vmtmt))
     CALL netcdf_err( nf90_put_att(ncid,id_mmm_vmtmt,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_mmm_vmtmt,'units',       &
          'm/sec'))

     base_label = '* multimode - poloidal momentum pinch ,m/sec'
     CALL set_label(label,base_label,strlen,'mmm_vmtmp')
     CALL netcdf_err( nf90_def_var(ncid, "mmm_vmtmp", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_mmm_vmtmp))
     CALL netcdf_err( nf90_put_att(ncid,id_mmm_vmtmp,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_mmm_vmtmp,'units',       &
          'm/sec'))


    ! end multimode declarations

 
     ! electron radiative losses due to all ions:
     base_label = '* radiative loss  species: '
     CALL set_label(label,base_label,strlen,'bremssthralung')
     CALL netcdf_err( nf90_def_var(ncid, "brems_nions", nf90_double,    &
          DIMIDS = (/idim_rho,idim_nion/),VARID=id_brems_nions))
     CALL netcdf_err( nf90_put_att(ncid,id_brems_nions,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_brems_nions,'units',       &
          'watts/meter^3'))

     CALL netcdf_err( nf90_def_var(ncid,"omegale", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_omegale))
     CALL netcdf_err(  nf90_put_att(ncid,id_omegale,'long_name',    &
          '* omegale, beam electron energy correction due to rotation, watts/meter^3'))
     CALL netcdf_err(  nf90_put_att(ncid,id_qrad,'units',        &
          'watts/meter^3'))


     CALL netcdf_err( nf90_def_var(ncid,"qomegapi", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_qomegapi))
     CALL netcdf_err(  nf90_put_att(ncid,id_qomegapi,'long_name',    &
          '* qomegapi, beam ion energy correction due to rotation, watts/meter^3'))
     CALL netcdf_err(  nf90_put_att(ncid,id_qomegapi,'units',        &
          'watts/meter^3'))

     CALL netcdf_err( nf90_def_var(ncid,"qangce", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_qangce))
     CALL netcdf_err(  nf90_put_att(ncid,id_qangce,'long_name',    &
          '* qangce, beam ion energy correction due to rotation, watts/meter^3'))
     CALL netcdf_err(  nf90_put_att(ncid,id_qangce,'units',        &
          'watts/meter^3'))

     CALL netcdf_err( nf90_def_var(ncid,"sprcxre", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_sprcxre))
     CALL netcdf_err(  nf90_put_att(ncid,id_sprcxre,'long_name',    &
          '* sprcxre, beam ion energy correction due to rotation, watts/meter^3'))
     CALL netcdf_err(  nf90_put_att(ncid,id_sprcxre,'units',        &
          'watts/meter^3'))

     CALL netcdf_err( nf90_def_var(ncid,"sprcxree", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_sprcxree))
     CALL netcdf_err(  nf90_put_att(ncid,id_sprcxree,'long_name',    &
          '* sprcxree, beam ion energy correction due to rotation, watts/meter^3'))
     CALL netcdf_err(  nf90_put_att(ncid,id_sprcxree,'units',        &
          'watts/meter^3'))

     CALL netcdf_err( nf90_def_var(ncid,"spreimpe", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_spreimpe))
     CALL netcdf_err(  nf90_put_att(ncid,id_spreimpe,'long_name',    &
          '* spreimpe, beam ion energy correction due to rotation, watts/meter^3'))
     CALL netcdf_err(  nf90_put_att(ncid,id_spreimpe,'units',        &
          'watts/meter^3'))


     ! total thermal d(t,n)he fusion heating power (electrons plus ions) 
     CALL netcdf_err( nf90_def_var(ncid, "pfuse_tot", nf90_double, id_pfuse_tot))
     CALL netcdf_err(  nf90_put_att(ncid,id_pfuse_tot,'long_name',      &
          'total electron plus ion thermal fusion heating power') )
     CALL netcdf_err(  nf90_put_att(ncid,id_pfuse_tot,'units', 'watts' ))

     ! totals for radiative losses
     CALL netcdf_err( nf90_def_var(ncid, "qrad_tot", nf90_double, id_qrad_tot))
     CALL netcdf_err(  nf90_put_att(ncid,id_qrad_tot,'long_name',      &
          'total electron power radiated from plasma') )
     CALL netcdf_err(  nf90_put_att(ncid,id_qrad_tot,'units', 'watts' ))

     label ='* total electron radiated power due to ion species  '
     CALL netcdf_err( nf90_def_var(ncid,"brems_tot", nf90_double, &
          DIMIDS = (/idim_nion /),VARID=id_brems_tot ))
     CALL netcdf_err( nf90_put_att(ncid,id_brems_tot,'long_name',label))
     CALL netcdf_err( nf90_put_att(ncid,id_brems_tot,'units',        &
          'watts'))

 
     CALL netcdf_err( nf90_def_var(ncid,"qohm", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_qohm))
     CALL netcdf_err(  nf90_put_att(ncid,id_qohm,'long_name',    &
          '*  (electron) ohmic power density, watts/meter^3'))
     CALL netcdf_err(  nf90_put_att(ncid,id_qohm,'units',        &
          'watts/meter^3'))

     CALL netcdf_err( nf90_def_var(ncid,"rmajavnpsi", nf90_double,    &
          DIMIDS = (/idim_npsi/),VARID=id_rmajavnpsi))
     CALL netcdf_err(  nf90_put_att(ncid,id_rmajavnpsi,'long_name',    &
          '*  average major radius of each flux surface, meter, evaluated at elevation of magnetic axis'))
     CALL netcdf_err(  nf90_put_att(ncid,id_rmajavnpsi,'units',        &
          'meter'))


     CALL netcdf_err( nf90_def_var(ncid,"rminavnpsi", nf90_double,    &
          DIMIDS = (/idim_npsi/),VARID=id_rminavnpsi))
     CALL netcdf_err(  nf90_put_att(ncid,id_rminavnpsi,'long_name',    &
          '*  average minor radius of each flux surface, meter, evaluated at elevation of magnetic axis'))
     CALL netcdf_err(  nf90_put_att(ncid,id_rminavnpsi,'units',        &
          ' meter'))



     CALL netcdf_err( nf90_def_var(ncid,"psivolpnpsi", nf90_double,    &
          DIMIDS = (/idim_npsi/),VARID=id_psivolpnpsi))
     CALL netcdf_err(  nf90_put_att(ncid,id_psivolpnpsi,'long_name',    &
          '*  volume of each flux surface, meter**3'))
     CALL netcdf_err(  nf90_put_att(ncid,id_psivolpnpsi,'units',        &
          ' meter^3'))

     CALL netcdf_err( nf90_def_var(ncid,"elongxnpsi", nf90_double,    &
          DIMIDS = (/idim_npsi/),VARID=id_elongxnpsi))
     CALL netcdf_err(  nf90_put_att(ncid,id_elongxnpsi,'long_name',    &
          '*  elongation of each flux surface'))
     CALL netcdf_err(  nf90_put_att(ncid,id_elongxnpsi,'units',        &
          dimensionless))

     CALL netcdf_err( nf90_def_var(ncid,"triangnpsi_u", nf90_double,    &
          DIMIDS = (/idim_npsi/),VARID=id_triangnpsi_u))
     CALL netcdf_err(  nf90_put_att(ncid,id_triangnpsi_u,'long_name',    &
          '*  upper triangularity of each flux surface'))
     CALL netcdf_err(  nf90_put_att(ncid,id_triangnpsi_u,'units',        &
          dimensionless))

     CALL netcdf_err( nf90_def_var(ncid,"triangnpsi_l", nf90_double,    &
          DIMIDS = (/idim_npsi/),VARID=id_triangnpsi_l))
     CALL netcdf_err(  nf90_put_att(ncid,id_triangnpsi_l,'long_name',    &
          '*  lower triangularity of each flux surface'))
     CALL netcdf_err(  nf90_put_att(ncid,id_triangnpsi_l,'units',        &
          dimensionless))

     CALL netcdf_err( nf90_def_var(ncid,"pindentnpsi", nf90_double,    &
          DIMIDS = (/idim_npsi/),VARID=id_pindentnpsi))
     CALL netcdf_err(  nf90_put_att(ncid,id_pindentnpsi,'long_name',    &
          '*  indentation of each flux surface'))
     CALL netcdf_err(  nf90_put_att(ncid,id_pindentnpsi,'units',        &
          dimensionless))

     CALL netcdf_err( nf90_def_var(ncid,"sfareanpsi", nf90_double,    &
          DIMIDS = (/idim_npsi/),VARID=id_sfareanpsi))
     CALL netcdf_err(  nf90_put_att(ncid,id_sfareanpsi,'long_name',    &
          '*  surface area of each flux surface, this is 4*pi*pi*R0*hcap*rho*<ABS(grad rho)>'))
     CALL netcdf_err(  nf90_put_att(ncid,id_sfareanpsi,'units',        &
          ' meter^2'))

     CALL netcdf_err( nf90_def_var(ncid,"cxareanpsi", nf90_double,    &
          DIMIDS = (/idim_npsi/),VARID=id_cxareanpsi))
     CALL netcdf_err(  nf90_put_att(ncid,id_cxareanpsi,'long_name',    &
          '*  cross-sectional area of each flux'))
     CALL netcdf_err(  nf90_put_att(ncid,id_cxareanpsi,'units',        &
          ' meter^2'))

 
     CALL netcdf_err( nf90_def_var(ncid,"grho1npsi", nf90_double,    &
          DIMIDS = (/idim_npsi/),VARID=id_grho1npsi))
     CALL netcdf_err(  nf90_put_att(ncid,id_grho1npsi,'long_name',    &
          '*  flux surface average absolute grad rho'))
     CALL netcdf_err(  nf90_put_att(ncid,id_grho1npsi,'units',        &
          dimensionless))

     CALL netcdf_err( nf90_def_var(ncid,"grho2npsi", nf90_double,    &
          DIMIDS = (/idim_npsi/),VARID=id_grho2npsi))
     CALL netcdf_err(  nf90_put_att(ncid,id_grho2npsi,'long_name',    &
          '*  flux surface average (grad rho)**2'))
     CALL netcdf_err(  nf90_put_att(ncid,id_grho2npsi,'units',        &
          dimensionless))

     !         CALL netcdf_err( nf90_def_var(ncid, "nplasbdry", nf90_int, id_nplasbdry))
     !         CALL netcdf_err(  nf90_put_att(ncid,id_nplasbdry,'long_name',      &
     !               '*  nplasdry : number of points on plasma boundary') )

    CALL netcdf_err( nf90_def_var(ncid,"qpsinpsi", nf90_double,         &
          DIMIDS = (/idim_npsi/),VARID=id_qpsinpsi))
    CALL netcdf_err(  nf90_put_att(ncid,id_qpsinpsi,'long_name',    &
          '*  q on eqdsk psigrid'))
    CALL netcdf_err(  nf90_put_att(ncid,id_qpsinpsi,'units',        &
          dimensionless))

    CALL netcdf_err( nf90_def_var(ncid,"pressnpsi", nf90_double,         &
          DIMIDS = (/idim_npsi/),VARID=id_pressnpsi))
    CALL netcdf_err(  nf90_put_att(ncid,id_pressnpsi,'long_name',    &
          '* pressure on eqdsk psigrid'))
    CALL netcdf_err(  nf90_put_att(ncid,id_pressnpsi,'units',        &
          'nt/m^2'))

    CALL netcdf_err( nf90_def_var(ncid,"ffprimnpsi", nf90_double,         &
          DIMIDS = (/idim_npsi/),VARID=id_ffprimnpsi))
    CALL netcdf_err(  nf90_put_att(ncid,id_ffprimnpsi,'long_name',    &
          '* ffprime  on eqdsk psigrid'))
    CALL netcdf_err(  nf90_put_att(ncid,id_ffprimnpsi,'units',        &
          'kg /( A s^2)'))
 


    CALL netcdf_err( nf90_def_var(ncid,"pprimnpsi", nf90_double,         &
          DIMIDS = (/idim_npsi/),VARID=id_pprimnpsi))
    CALL netcdf_err(  nf90_put_att(ncid,id_pprimnpsi,'long_name',    &
          '* pprime  on eqdsk psigrid'))
    CALL netcdf_err(  nf90_put_att(ncid,id_pprimnpsi,'units',        &
          ' A/m^3'))


     CALL netcdf_err( nf90_def_var(ncid,"rplasbdry", nf90_double,         &
          DIMIDS = (/idim_nplasbdry/),VARID=id_rplasbdry))
     CALL netcdf_err(  nf90_put_att(ncid,id_rplasbdry,'long_name',    &
          '*  r points for plasma boundary, meter'))
     CALL netcdf_err(  nf90_put_att(ncid,id_rplasbdry,'units',        &
          'meter'))
     CALL netcdf_err( nf90_def_var(ncid,"zplasbdry", nf90_double,         &
          DIMIDS = (/idim_nplasbdry/),VARID=id_zplasbdry))
     CALL netcdf_err(  nf90_put_att(ncid,id_zplasbdry,'long_name',    &
          '*  z points for plasma boundary, meter'))
     CALL netcdf_err(  nf90_put_att(ncid,id_zplasbdry,'units',        &
          'meter'))

     CALL netcdf_err( nf90_def_var(ncid,"rlimiter", nf90_double,         &
          DIMIDS = (/idim_nlimiter/),VARID=id_rlimiter))
     CALL netcdf_err(  nf90_put_att(ncid,id_rlimiter,'long_name',    &
          '*  R points for limiter, meter'))
     CALL netcdf_err(  nf90_put_att(ncid,id_rlimiter,'units',        &
          'meter'))

     CALL netcdf_err( nf90_def_var(ncid,"zlimiter", nf90_double,         &
          DIMIDS = (/idim_nlimiter/),VARID=id_zlimiter))
     CALL netcdf_err(  nf90_put_att(ncid,id_zlimiter,'long_name',    &
          '*  Z points for limiter, meter'))
     CALL netcdf_err(  nf90_put_att(ncid,id_zlimiter,'units',        &
          'meter'))

     CALL netcdf_err( nf90_def_var(ncid,"storqueb", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_storqueb))
     CALL netcdf_err(  nf90_put_att(ncid,id_storqueb,'long_name',    &
          '* beam torque density, newton*m/m**3'))
     CALL netcdf_err(  nf90_put_att(ncid,id_storqueb,'units',        &
          ' newton*meter/meter^3'))

     CALL netcdf_err( nf90_def_var(ncid,"totcur_bc ", nf90_double,          &
          VARID=id_totcur_bc ))
     label = '* boundary condition total current, amp, at time '     &
          //bc_asc_time(1:LEN_TRIM(bc_asc_time))
     CALL netcdf_err(  nf90_put_att(ncid,id_totcur_bc ,'long_name',label))
     CALL netcdf_err(  nf90_put_att(ncid,id_totcur_bc ,'units',        &
          'amp'))


     CALL netcdf_err( nf90_def_var(ncid,"vloop_bc ", nf90_double,          &
          VARID=id_vloop_bc ))
     label = '* boundary condition loop voltage, Volts, at time '//bc_asc_time(1:LEN_TRIM(bc_asc_time))

     CALL netcdf_err(  nf90_put_att(ncid,id_vloop_bc ,'long_name',label))
     CALL netcdf_err(  nf90_put_att(ncid,id_vloop_bc ,'units',        &
          'volt'))
 

     CALL netcdf_err( nf90_def_var(ncid,"fix_edge_te_bc", nf90_double,          &
          VARID=id_fix_edge_te_bc))
     label ='* te boundary condition rho flags at time ' &
          //bc_asc_time(1:LEN_TRIM(bc_asc_time))
     CALL netcdf_err( nf90_put_att(ncid,id_fix_edge_te_bc,'long_name',label))
     CALL netcdf_err( nf90_put_att(ncid,id_fix_edge_te_bc,'units',        &
          dimensionless))

     label ='* ti boundary condition rho flags at time ' &
          //bc_asc_time(1:LEN_TRIM(bc_asc_time))
     CALL netcdf_err( nf90_def_var(ncid,"fix_edge_ti_bc", nf90_double,    &
          VARID=id_fix_edge_ti_bc ))
     CALL netcdf_err( nf90_put_att(ncid,id_fix_edge_ti_bc,'long_name',label))
     CALL netcdf_err( nf90_put_att(ncid,id_fix_edge_ti_bc,'units',        &
          dimensionless))

     label ='* rot boundary condition rho flags at time ' &
          //bc_asc_time(1:LEN_TRIM(bc_asc_time))
     CALL netcdf_err( nf90_def_var(ncid,"fix_edge_rot_bc ", nf90_double,   &
          VARID=id_fix_edge_rot_bc ))
     CALL netcdf_err( nf90_put_att(ncid,id_fix_edge_rot_bc,'long_name',label))
     CALL netcdf_err( nf90_put_att(ncid,id_fix_edge_rot_bc,'units',        &
          dimensionless))

     label ='* ni boundary condition rho flags at time ' &
          //bc_asc_time(1:LEN_TRIM(bc_asc_time))
     CALL netcdf_err( nf90_def_var(ncid,"fix_edge_ni_bc ", nf90_double, &
          DIMIDS = (/idim_nion /),VARID=id_fix_edge_ni_bc ))
     CALL netcdf_err( nf90_put_att(ncid,id_fix_edge_ni_bc,'long_name',label))
     CALL netcdf_err( nf90_put_att(ncid,id_fix_edge_ni_bc,'units',        &
          dimensionless))


     CALL netcdf_err( nf90_def_var(ncid,"te_bc", nf90_double,          &
          DIMIDS = (/idim_rho/),VARID=id_te_bc ))
     label = '* boundary condition Te, keV, at time '//bc_asc_time(1:LEN_TRIM(bc_asc_time))
     CALL netcdf_err(  nf90_put_att(ncid,id_te_bc,'long_name',label))
     CALL netcdf_err(  nf90_put_att(ncid,id_te_bc ,'units',        &
          'keV'))

     CALL netcdf_err( nf90_def_var(ncid,"ti_bc", nf90_double,          &
          DIMIDS = (/idim_rho/),VARID=id_ti_bc ))
     label = '* boundary condition Ti, keV, at time '//bc_asc_time(1:LEN_TRIM(bc_asc_time))
     CALL netcdf_err(  nf90_put_att(ncid,id_ti_bc,'long_name',label))
     CALL netcdf_err(  nf90_put_att(ncid,id_ti_bc ,'units',        &
          'keV'))

     CALL netcdf_err( nf90_def_var(ncid,"ene_bc", nf90_double,          &
          DIMIDS = (/idim_rho/),VARID=id_ene_bc ))
     label ='* bc profile: ene, 1/m**3 at time  '//bc_asc_time(1:LEN_TRIM(bc_asc_time))

     CALL netcdf_err(  nf90_put_att(ncid,id_ene_bc,'long_name',label))
     CALL netcdf_err(  nf90_put_att(ncid,id_ene_bc ,'units',        &
          '1/meter^3'))

     CALL netcdf_err( nf90_def_var(ncid,"zeff_bc", nf90_double,          &
          DIMIDS = (/idim_rho/),VARID=id_zeff_bc ))
     label = '* bc profile: zeff, at time  '//bc_asc_time(1:LEN_TRIM(bc_asc_time))
     CALL netcdf_err(  nf90_put_att(ncid,id_zeff_bc,'long_name',label))
     CALL netcdf_err(  nf90_put_att(ncid,id_zeff_bc ,'units',        &
          dimensionless))
 
     CALL netcdf_err( nf90_def_var(ncid,"angrot_bc", nf90_double,          &
          DIMIDS = (/idim_rho/),VARID=id_angrot_bc ))
     label = '*  angular rotation speed profile, rad/sec, at time '        &
          //bc_asc_time(1:LEN_TRIM(bc_asc_time))
     CALL netcdf_err( nf90_put_att(ncid,id_angrot_bc,'long_name',label))
     CALL netcdf_err( nf90_put_att(ncid,id_angrot_bc,'units',        &
          'radian/second'))

     CALL netcdf_err( nf90_def_var(ncid,"wbeam", nf90_double,          &
          DIMIDS = (/idim_rho/),VARID=id_wbeam ))
     label = '* fast ion stored energy density KEV/m**3 '
     CALL netcdf_err(  nf90_put_att(ncid,id_wbeam,'long_name',label))
     CALL netcdf_err(  nf90_put_att(ncid,id_wbeam ,'units',        &
          'keV/meter^3'))

     CALL netcdf_err( nf90_def_var(ncid,"walp", nf90_double,          &
          DIMIDS = (/idim_rho/),VARID=id_walp))
     label = '* fast alpha stored energy density KEV/m**3 ' 
     CALL netcdf_err(  nf90_put_att(ncid,id_walp,'long_name',label))
     CALL netcdf_err(  nf90_put_att(ncid,id_walp ,'units',        &
          'keV/meter^3'))

     CALL netcdf_err( nf90_def_var(ncid,"enalp", nf90_double,          &
          DIMIDS = (/idim_rho/),VARID=id_enalp))
     label = '* fast alpha density 1/m**3'
     CALL netcdf_err(  nf90_put_att(ncid,id_enalp,'long_name',label))
     CALL netcdf_err(  nf90_put_att(ncid,id_enalp ,'units',        &
          '1/meter^3'))



     CALL netcdf_err( nf90_def_var(ncid,"eps", nf90_double,          &
          DIMIDS = (/idim_rho/),VARID=id_eps))
     label = '* horizontal inverse aspect ratio = (rmax-rmin)/(rmax+rmin)'
     CALL netcdf_err( nf90_put_att(ncid,id_eps,'long_name',label))
     CALL netcdf_err( nf90_put_att(ncid,id_eps,'units',        &
          dimensionless))





     CALL netcdf_err( nf90_def_var(ncid,"rcap", nf90_double,            &
          DIMIDS = (/idim_rho/),VARID=id_rcap))
     label = '* rcap = < R>, m'
     CALL netcdf_err(  nf90_put_att(ncid,id_rcap,'long_name',label))
     CALL netcdf_err(  nf90_put_att(ncid,id_rcap ,'units',        &
          'meter'))

 
     CALL netcdf_err( nf90_def_var(ncid,"rcapi", nf90_double,            &
          DIMIDS = (/idim_rho/),VARID=id_rcapi))
     label = '* rcap i= <1/ R>, 1/m'
     CALL netcdf_err(  nf90_put_att(ncid,id_rcapi,'long_name',label))
     CALL netcdf_err(  nf90_put_att(ncid,id_rcapi ,'units',        &
          '1/meter'))




     CALL netcdf_err( nf90_def_var(ncid,"r2cap", nf90_double,           &
          DIMIDS = (/idim_rho/),VARID=id_r2cap))
     label = '*r2cap = <R0**2/R**2>'
     CALL netcdf_err( nf90_put_att(ncid,id_r2cap,'long_name',label))
     CALL netcdf_err( nf90_put_att(ncid,id_r2cap,'units',        &
          dimensionless))

     CALL netcdf_err( nf90_def_var(ncid,"r2capi", nf90_double,          &
          DIMIDS = (/idim_rho/),VARID=id_r2capi))
     label = '*r2capi = <R**2>, m**2'
     CALL netcdf_err(  nf90_put_att(ncid,id_r2capi,'long_name',label))
     CALL netcdf_err(  nf90_put_att(ncid,id_r2capi ,'units',            &
          'meter^2'))



     CALL netcdf_err( nf90_def_var(ncid,"xhm2", nf90_double,          &
          DIMIDS = (/idim_rho/),VARID=id_xhm2))
     label = '*xhm2 = < (B total/ B axis)**2 > (=1 for circular plasmas)'
     CALL netcdf_err( nf90_put_att(ncid,id_xhm2,'long_name',label))
     CALL netcdf_err( nf90_put_att(ncid,id_xhm2,'units',        &
          dimensionless))


     CALL netcdf_err( nf90_def_var(ncid,"xi11", nf90_double,          &
          DIMIDS = (/idim_rho/),VARID=id_xi11))
     label ='*xi11 ( = 1.95 sqrt(eps)for circular plasmas)'
     CALL netcdf_err( nf90_put_att(ncid,id_xi11,'long_name',label))
     CALL netcdf_err( nf90_put_att(ncid,id_xi11,'units',        &
          dimensionless))



     CALL netcdf_err( nf90_def_var(ncid,"xi33", nf90_double,          &
          DIMIDS = (/idim_rho/),VARID=id_xi33))
     label ='*xi33 ( = 1.95 sqrt(eps)for circular plasmas)'
     CALL netcdf_err( nf90_put_att(ncid,id_xi33,'long_name',label))
     CALL netcdf_err( nf90_put_att(ncid,id_xi33,'units',        &
          dimensionless))

     CALL netcdf_err( nf90_def_var(ncid,"xips", nf90_double,          &
          DIMIDS = (/idim_rho/),VARID=id_xips))
     label = '*xips = <(Baxis/B)**2)> - 1./(<(B/Baxis)**2> )( = 2 eps**2 for circular plasmas)'
     CALL netcdf_err( nf90_put_att(ncid,id_xips,'long_name',label))
     CALL netcdf_err( nf90_put_att(ncid,id_xips,'units',        &
          dimensionless))




     CALL netcdf_err( nf90_def_var(ncid,"dfdt", nf90_double,          &
          DIMIDS = (/idim_rho/),VARID=id_dfdt))
     label = '(d/dt)Fcap '
     CALL netcdf_err(  nf90_put_att(ncid,id_dfdt,'long_name',label))
     CALL netcdf_err(  nf90_put_att(ncid,id_dfdt,'units',        &
          '1/sec'))



     CALL netcdf_err( nf90_def_var(ncid,"dgdt", nf90_double,          &
          DIMIDS = (/idim_rho/),VARID=id_dgdt))
     label = '(d/dt)Gcap '
     CALL netcdf_err(  nf90_put_att(ncid,id_dgdt,'long_name',label))
     CALL netcdf_err(  nf90_put_att(ncid,id_dgdt ,'units',        &
          '1/sec'))



     CALL netcdf_err( nf90_def_var(ncid,"dhdt", nf90_double,          &
          DIMIDS = (/idim_rho/),VARID=id_dhdt))
     label = '(d/dt)Hcap '
     CALL netcdf_err(  nf90_put_att(ncid,id_dhdt,'long_name',label))
     CALL netcdf_err(  nf90_put_att(ncid,id_dhdt,'units',        &
          '1/sec'))
 

     !nu*  ion collison/bounce frequencies:
     base_label = 'nu* ion collison/bounce freq, species: '
     CALL set_label(label,base_label,strlen,'nion')
     CALL netcdf_err( nf90_def_var(ncid, "xnus", nf90_double,    &
          DIMIDS = (/idim_rho,idim_nion/),VARID=id_xnus))
     CALL netcdf_err( nf90_put_att(ncid,id_xnus,'long_name',label ))
     CALL netcdf_err( nf90_put_att(ncid,id_xnus,'units',        &
          dimensionless))


     !nue*  electron collison/bounce frequencies:
     label = 'nu*e  electron collison/bounce freq '
     CALL netcdf_err( nf90_def_var(ncid, "xnuse", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_xnuse))
     CALL netcdf_err( nf90_put_att(ncid,id_xnuse,'long_name',label ))
     CALL netcdf_err( nf90_put_att(ncid,id_xnuse,'units',        &
          dimensionless))


     !ftrap electron trapped particle fraction:
     label = 'electron trapped particle fraction '
     CALL netcdf_err( nf90_def_var(ncid, "ftrap", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_ftrap))
     CALL netcdf_err( nf90_put_att(ncid,id_ftrap,'long_name',label ))
     CALL netcdf_err( nf90_put_att(ncid,id_ftrap,'units',        &
          dimensionless))


     !eta, resistivity:
     label = 'eta  resistivity ohm m'
     CALL netcdf_err( nf90_def_var(ncid, "eta", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_eta))
     CALL netcdf_err( nf90_put_att(ncid,id_eta,'long_name',label ))
     CALL netcdf_err( nf90_put_att(ncid,id_eta ,'units',        &
          'ohm*meter'))
     
     !chiepc, paleoclassical diffusivity
     label = 'chiepc Electron Paleoclassical Diffusivity m^2/s'
     CALL netcdf_err( nf90_def_var(ncid,"chiepc", nf90_double,  &
          DIMIDS = (/idim_rho/),VARID=id_chiepc))
     CALL netcdf_err( nf90_put_att(ncid,id_chiepc,'long_name',label ))
     CALL netcdf_err( nf90_put_att(ncid,id_chiepc,'units',      &
          'meter^2/sec'))


     !------------------------------------------------------------------------
     ! neutron rates, profiles and totals:
     !------------------------------------------------------------------------
     label = 'thermal -thermal neutron rate #/m^3 sec'
     CALL netcdf_err( nf90_def_var(ncid, "neutr_ddn_th", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_neutr_ddn_th))
     CALL netcdf_err( nf90_put_att(ncid,id_neutr_ddn_th,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_neutr_ddn_th,'units',        &
          '#/(m^3 sec)'))

     label = 'beam -thermal neutron rate #/m^3 sec'
     CALL netcdf_err( nf90_def_var(ncid, "neutr_ddn_beam_thermal", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_neutr_ddn_beam_thermal))
     CALL netcdf_err( nf90_put_att(ncid,id_neutr_ddn_beam_thermal,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_neutr_ddn_beam_thermal,'units',        &
          '#/(m^3 sec)'))
 

     label = 'beam - beam neutron rate #/m^3 sec'
     CALL netcdf_err( nf90_def_var(ncid, "neutr_ddn_beam_beam", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_neutr_ddn_beam_beam))
     CALL netcdf_err( nf90_put_att(ncid,id_neutr_ddn_beam_beam,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_neutr_ddn_beam_beam,'units',        &
          '#/(m^3 sec)'))

     label = 'knock on neutron rate #/m^3 sec'
     CALL netcdf_err( nf90_def_var(ncid, "neutr_ddn_knock", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_neutr_ddn_knock))
     CALL netcdf_err( nf90_put_att(ncid,id_neutr_ddn_knock,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_neutr_ddn_knock,'units',        &
          '#/(m^3 sec)')) ! ==> unterminated string message here in pgf95

     label = 'total neutron rate #/m^3 sec'
     CALL netcdf_err( nf90_def_var(ncid, "neutr_ddn_tot", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_neutr_ddn_tot))
     CALL netcdf_err( nf90_put_att(ncid,id_neutr_ddn_tot,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_neutr_ddn_tot,'units',        &
          '#/(m^3 sec)'))

     CALL netcdf_err( nf90_def_var(ncid, "total_neutr_ddn_th", nf90_double,id_total_neutr_ddn_th))
     CALL netcdf_err(  nf90_put_att(ncid,id_total_neutr_ddn_th,'long_name',           &
          '*  total_neutr_ddn_th: total thermal-thermal neutron rate 1/sec') )
     CALL netcdf_err(  nf90_put_att(ncid,id_total_neutr_ddn_th,'units', '1/sec'  ))

     CALL netcdf_err( nf90_def_var(ncid, "total_neutr_ddn_beam_beam", nf90_double,id_total_neutr_ddn_beam_beam))
     CALL netcdf_err(  nf90_put_att(ncid,id_total_neutr_ddn_beam_beam,'long_name',           &
          '*  total_neutr_ddn_beam_beam: total beam - beam neutron rate 1/sec') )
     CALL netcdf_err(  nf90_put_att(ncid,id_total_neutr_ddn_beam_beam,'units', '1/sec'  ))

     CALL netcdf_err( nf90_def_var(ncid, "total_neutr_ddn_beam_thermal", nf90_double,id_total_neutr_ddn_beam_thermal))
     CALL netcdf_err(  nf90_put_att(ncid,id_total_neutr_ddn_beam_thermal,'long_name',           &
          '*  total_neutr_ddn_beam_thermal: total beam - thermal neutron rate 1/sec') )
     CALL netcdf_err(  nf90_put_att(ncid,id_total_neutr_ddn_beam_thermal,'units', '1/sec'  ))

     CALL netcdf_err( nf90_def_var(ncid, "total_neutr_ddn_knock", nf90_double,id_total_neutr_ddn_knock))
     CALL netcdf_err(  nf90_put_att(ncid,id_total_neutr_ddn_knock,'long_name',           &
          '*  total_neutr_ddn_knock: total knock on  neutron rate, 1/sec') )
     CALL netcdf_err(  nf90_put_att(ncid,id_total_neutr_ddn_knock,'units', '1/sec'  ))

     CALL netcdf_err( nf90_def_var(ncid, "total_neutr_ddn", nf90_double,id_total_neutr_ddn))
     CALL netcdf_err(  nf90_put_att(ncid,id_total_neutr_ddn,'long_name',           &
          '*  total_neutr_ddn : total  neutron rate, 1/sec') )
     CALL netcdf_err(  nf90_put_att(ncid,id_total_neutr_ddn,'units', '1/sec'  ))


     !------------------------------------------------------------------------
     ! aneutronic fusion, d(d,p)t,he3(d,p)he4 profiles and totals:
     !------------------------------------------------------------------------
 
     label = 'd(d,p)t rate #/(m^3 sec)'
     CALL netcdf_err( nf90_def_var(ncid, "ddpt", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_p_ddt))
     CALL netcdf_err( nf90_put_att(ncid,id_p_ddt,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_p_ddt,'units',        &
          '#/(m^3 sec)')) 

     CALL netcdf_err( nf90_def_var(ncid, "ddpt_tot", nf90_double,id_total_p_ddt))
     CALL netcdf_err(  nf90_put_att(ncid,id_total_p_ddt,'long_name',           &
          '*  total d(d,p)t reaction rate, 1/sec') )

! --- he3(dp,he4) thermal and beam thermal:
     label = 'he3(d,p)he4 thermal rate #/(m^3 sec)'
     CALL netcdf_err( nf90_def_var(ncid, "he3dp_th", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_p_he3d_th))
     CALL netcdf_err( nf90_put_att(ncid,id_p_he3d_th,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_p_he3d_th,'units',        &
          '#/(m^3 sec)')) 
     CALL netcdf_err( nf90_def_var(ncid, "he3dp_th_tot", nf90_double,id_total_p_he3d_th))
     CALL netcdf_err(  nf90_put_att(ncid,id_total_p_he3d_th,'long_name',           &
          '*  total thermal He3(d,p)He4 reaction rate, 1/sec') )

     label = 'he3(d,p)he4 beam-thermal rate #/(m^3 sec)'
     CALL netcdf_err( nf90_def_var(ncid, "he3dp_beam_th", nf90_double,    &
          DIMIDS = (/idim_rho/),VARID=id_p_he3d_beam_th))
     CALL netcdf_err( nf90_put_att(ncid,id_p_he3d_beam_th,'long_name',label ))
     CALL netcdf_err(  nf90_put_att(ncid,id_p_he3d_beam_th,'units',        &
          '#/(m^3 sec)')) 

     CALL netcdf_err( nf90_def_var(ncid, "he3dp_beam_th_tot", nf90_double,id_total_p_he3d_beam_th))
     CALL netcdf_err(  nf90_put_att(ncid,id_total_p_he3d_beam_th,'long_name',           &
          '*  total beam thermal He3(d,p)He4 reaction rate, 1/sec') )

     CALL netcdf_err( nf90_def_var(ncid, "he3dp_tot", nf90_double,id_total_p_he3d))
     CALL netcdf_err(  nf90_put_att(ncid,id_total_p_he3d,'long_name',   &
          '*  total beam plus thermal He3(d,p)He4 reaction rate, 1/sec') )

     CALL netcdf_err( nf90_def_var(ncid, "he3_frac", nf90_double,id_he3_frac))
     CALL netcdf_err(  nf90_put_att(ncid,id_he3_frac,'long_name',   &
          '*  fraction of He that is He3') )

     CALL netcdf_err( nf90_def_var(ncid, "he3_thermal_spin_pol", nf90_double,id_he3_thermal_spin_pol))
     CALL netcdf_err(  nf90_put_att(ncid,id_he3_thermal_spin_pol,'long_name',   &
          '* thermal He3 spin polarization relative to magnetic field direction') )

     CALL netcdf_err( nf90_def_var(ncid, "d_beam_spin_pol", nf90_double,id_d_beam_spin_pol))
     CALL netcdf_err(  nf90_put_att(ncid,id_d_beam_spin_pol,'long_name',   &
          '*  spin polarization of injected d beam') )

     !------------------------------------------------------------------------------
     !neutral beam quantities:
     !NOTE: neut_beam%nbeams > 0 if nubeam OR P_Nfreya is used
     !      But nubeam does not supply all the necessary info
     !      ALSO rho grid dimension is fixed at input value
     !      and identified as nj_beam with id_nbeam_rho
     !-----------------------------------------------------------------------------------
 
     nb_loop: IF(neut_beam%nbeams .GT. 0)THEN
        CALL  neut_beam_allocate

        base_label = '* neutral beam power before aperture: '
        CALL set_label(label,base_label,strlen,'pbeam')
        CALL netcdf_err( nf90_def_var(ncid, "pbeam", nf90_double,    &
             DIMIDS = (/idim_ke_bm,idim_nbeams/),VARID=id_pbeam))
        CALL netcdf_err( nf90_put_att(ncid,id_pbeam,'long_name',label ))
        CALL netcdf_err(  nf90_put_att(ncid,id_pbeam,'units',       &
                          'Watts'))

        base_label = '* neutral beam power into torus or  thermalized power if nubeam used: '
        CALL set_label(label,base_label,strlen,'bptor')
        CALL netcdf_err( nf90_def_var(ncid, "bptor", nf90_double,    &
             DIMIDS = (/idim_nbeams/),VARID=id_bptor))
        CALL netcdf_err( nf90_put_att(ncid,id_bptor,'long_name',label ))
        CALL netcdf_err(  nf90_put_att(ncid,id_bptor,'units',       &
                          'Watts'))

        base_label = '* neutral beam  current fractions: '
        CALL set_label(label,base_label,strlen,'fbcur')
        CALL netcdf_err( nf90_def_var(ncid, "fbcur", nf90_double,    &
             DIMIDS = (/idim_ke_bm,idim_nbeams/),VARID=id_fbcur))
        CALL netcdf_err( nf90_put_att(ncid,id_fbcur,'long_name',label ))
        CALL netcdf_err(  nf90_put_att(ncid,id_fbcur,'units',       &
                          'None'))

        base_label = '* prompt neutral beam  power in plasma: '
        CALL set_label(label,base_label,strlen,'prompt_nb_pwr')
        CALL netcdf_err( nf90_def_var(ncid, "prompt_nb_pwr", nf90_double,    &
             DIMIDS = (/idim_ke_bm,idim_nbeams/),VARID=id_prompt_nb_pwr))
        CALL netcdf_err( nf90_put_att(ncid,id_prompt_nb_pwr,'long_name',label ))
        CALL netcdf_err(  nf90_put_att(ncid,id_prompt_nb_pwr,'units',       &
                          'Watts'))


        base_label = '* neutral beam energies: '
        CALL set_label(label,base_label,strlen,'ebeam')
        CALL netcdf_err( nf90_def_var(ncid, "ebeam", nf90_double,    &
             DIMIDS = (/idim_ke_bm,idim_nbeams/),VARID=id_ebeam))
        CALL netcdf_err( nf90_put_att(ncid,id_ebeam,'long_name',label ))
        CALL netcdf_err(  nf90_put_att(ncid,id_ebeam,'units',       &
                          'Kev'))

        base_label = '* neutral beam neutral intensity: '
        CALL set_label(label,base_label,strlen,'bneut')
        CALL netcdf_err( nf90_def_var(ncid, "bneut", nf90_double,    &
             DIMIDS = (/idim_ke_bm,idim_nbeams/),VARID=id_bneut))
        CALL netcdf_err( nf90_put_att(ncid,id_bneut,'long_name',label ))
!        CALL netcdf_err(  nf90_put_att(ncid,id_bneut,'units',       &
!                          '#/sec'))

        base_label = '* neutral beam ion intensity: '
        CALL set_label(label,base_label,strlen,'bion')
        CALL netcdf_err( nf90_def_var(ncid, "bion", nf90_double,    &
             DIMIDS = (/idim_ke_bm,idim_nbeams/),VARID=id_bion))
        CALL netcdf_err( nf90_put_att(ncid,id_bion,'long_name',label ))
!        CALL netcdf_err(  nf90_put_att(ncid,id_bion,'units',       &
!                          '#/sec'))

        base_label = '* neutral beam aperture loss: '
        CALL set_label(label,base_label,strlen,'fap')
        CALL netcdf_err( nf90_def_var(ncid, "fap", nf90_double,    &
             DIMIDS = (/idim_ke_bm,idim_nbeams/),VARID=id_fap))
        CALL netcdf_err( nf90_put_att(ncid,id_fap,'long_name',label ))
!        CALL netcdf_err(  nf90_put_att(ncid,id_fap,'units',       &
!                          'none'))

        base_label = '* neutral beam wall loss: '
        CALL set_label(label,base_label,strlen,'fwall')
        CALL netcdf_err( nf90_def_var(ncid, "fwall", nf90_double,    &
             DIMIDS = (/idim_ke_bm,idim_nbeams/),VARID=id_fwall))
        CALL netcdf_err( nf90_put_att(ncid,id_fwall,'long_name',label ))
!        CALL netcdf_err(  nf90_put_att(ncid,id_fwall,'units',       &
!                          'none'))

        base_label = '* neutral beam orbit loss: '
        CALL set_label(label,base_label,strlen,'forb')
        CALL netcdf_err( nf90_def_var(ncid, "forb", nf90_double,    &
             DIMIDS = (/idim_ke_bm,idim_nbeams/),VARID=id_forb))
        CALL netcdf_err( nf90_put_att(ncid,id_forb,'long_name',label ))
!        CALL netcdf_err(  nf90_put_att(ncid,id_forb,'units',       &
!                          'none'))

        base_label = '* neutral beam fraction of ions trapped for which error was detected: '
        CALL set_label(label,base_label,strlen,'fber')
        CALL netcdf_err( nf90_def_var(ncid, "fber", nf90_double,    &
             DIMIDS = (/idim_ke_bm,idim_nbeams/),VARID=id_fber))
        CALL netcdf_err( nf90_put_att(ncid,id_fber,'long_name',label ))
        CALL netcdf_err(  nf90_put_att(ncid,id_fber,'units',       &
                          dimensionless))

        base_label = '*  neutral beam fraction of ions trapped and not encircling : '
        CALL set_label(label,base_label,strlen,'fb00')
        CALL netcdf_err( nf90_def_var(ncid, "fb00", nf90_double,    &
             DIMIDS = (/idim_ke_bm,idim_nbeams/),VARID=id_fb00))
        CALL netcdf_err( nf90_put_att(ncid,id_fb00,'long_name',label ))
        CALL netcdf_err(  nf90_put_att(ncid,id_fb00,'units',       &
                          dimensionless))
 
        base_label = '* neutral beam fraction of ions trapped and axis-encircling: '
        CALL set_label(label,base_label,strlen,'fb01')
        CALL netcdf_err( nf90_def_var(ncid, "fb01", nf90_double,    &
             DIMIDS = (/idim_ke_bm,idim_nbeams/),VARID=id_fb01))
        CALL netcdf_err( nf90_put_att(ncid,id_fb01,'long_name',label ))
        CALL netcdf_err(  nf90_put_att(ncid,id_fb01,'units',       &
                          dimensionless))

        base_label = '* Neutral beam fraction of ions passing and not encircling : ' 
        CALL set_label(label,base_label,strlen,'fb10')
        CALL netcdf_err( nf90_def_var(ncid, "fb10", nf90_double,    &
             DIMIDS = (/idim_ke_bm,idim_nbeams/),VARID=id_fb10))
        CALL netcdf_err( nf90_put_att(ncid,id_fb10,'long_name',label ))
        CALL netcdf_err(  nf90_put_att(ncid,id_fb10,'units',       &
                          dimensionless))

        base_label = '* Neutral beam fraction of ions passing and axis-encircling: ' 
        CALL set_label(label,base_label,strlen,'fb11')
        CALL netcdf_err( nf90_def_var(ncid, "fb11", nf90_double,    &
             DIMIDS = (/idim_ke_bm,idim_nbeams/),VARID=id_fb11))
        CALL netcdf_err( nf90_put_att(ncid,id_fb11,'long_name',label ))
        CALL netcdf_err(  nf90_put_att(ncid,id_fb11,'units',       &
                          dimensionless))


        base_label = '* Neutral beam orbit width of (m) of trapped and not encircling ions : '
        CALL set_label(label,base_label,strlen,'wb00')
        CALL netcdf_err( nf90_def_var(ncid, "wb00", nf90_double,    &
             DIMIDS = (/idim_ke_bm,idim_nbeams/),VARID=id_wb00))
        CALL netcdf_err( nf90_put_att(ncid,id_wb00,'long_name',label ))
        CALL netcdf_err(  nf90_put_att(ncid,id_wb00,'units','m'))

        base_label = '* Neutral beam orbit width (m) of ions trapped and axis-encircling:'
        CALL set_label(label,base_label,strlen,'wb01')
        CALL netcdf_err( nf90_def_var(ncid, "wb01", nf90_double,    &
             DIMIDS = (/idim_ke_bm,idim_nbeams/),VARID=id_wb01))
        CALL netcdf_err( nf90_put_att(ncid,id_wb01,'long_name',label ))
        CALL netcdf_err(  nf90_put_att(ncid,id_wb01,'units', 'm'))

        base_label = '* Neutral beam orbit width (m) of ions passing and not encircling:'
        CALL set_label(label,base_label,strlen,'wb10')
        CALL netcdf_err( nf90_def_var(ncid, "wb10", nf90_double,    &
             DIMIDS = (/idim_ke_bm,idim_nbeams/),VARID=id_wb10))
        CALL netcdf_err( nf90_put_att(ncid,id_wb10,'long_name',label ))
        CALL netcdf_err(  nf90_put_att(ncid,id_wb10,'units', 'm'))


        base_label = '* Neutral beam orbit width of ions passing and axis-encircling:'
        CALL set_label(label,base_label,strlen,'wb11')
        CALL netcdf_err( nf90_def_var(ncid, "wb11", nf90_double,    &
             DIMIDS = (/idim_ke_bm,idim_nbeams/),VARID=id_wb11))
        CALL netcdf_err( nf90_put_att(ncid,id_wb11,'long_name',label ))
        CALL netcdf_err(  nf90_put_att(ncid,id_wb11,'units', 'm'))
 
        base_label = '* neutral beam fast ion source: '
        CALL set_label(label,base_label,strlen,'sb')
        CALL netcdf_err( nf90_def_var(ncid, "sb", nf90_double,    &
             DIMIDS = (/idim_beam_rho,idim_ke_bm,idim_nbeams/),VARID=id_sb))
 
        CALL netcdf_err( nf90_put_att(ncid,id_sb,'long_name',label ))
        CALL netcdf_err(  nf90_put_att(ncid,id_sb,'units',       &
                          '1/meter^3 sec'))

        base_label = '* neutral beam fast ion power source: '
        CALL set_label(label,base_label,strlen,'qb')
        CALL netcdf_err( nf90_def_var(ncid, "qb", nf90_double,    &
             DIMIDS = (/idim_beam_rho,idim_ke_bm,idim_nbeams/),VARID=id_qb))
        CALL netcdf_err( nf90_put_att(ncid,id_qb,'long_name',label ))
        CALL netcdf_err(  nf90_put_att(ncid,id_qb,'units',       &
                          'watts/meter^3 '))
 
        base_label = '* neutral beam prompt parallel momentum rate density,  kg/(m2-s2) '
        CALL set_label(label,base_label,strlen,'spb')
        CALL netcdf_err( nf90_def_var(ncid, "spb", nf90_double,    &
             DIMIDS = (/idim_beam_rho,idim_ke_bm,idim_nbeams/),VARID=id_spb))
        CALL netcdf_err( nf90_put_att(ncid,id_spb,'long_name',label ))
        CALL netcdf_err(  nf90_put_att(ncid,id_spb,'units',       &
                          'kg/(meter^2 sec^2) '))
 
        base_label = '* neutral beam prompt angular momentum rate density,  kg/(m s^2)'
        CALL set_label(label,base_label,strlen,'spbr')
        CALL netcdf_err( nf90_def_var(ncid, "spbr", nf90_double,    &
             DIMIDS = (/idim_beam_rho,idim_ke_bm,idim_nbeams/),VARID=id_spbr))
        CALL netcdf_err( nf90_put_att(ncid,id_spbr,'long_name',label ))
        CALL netcdf_err(  nf90_put_att(ncid,id_spbr,'units',       &
                          'kg/(meter^2 sec^2) '))

        base_label = '* neutral beam total ang momt in zone/deposition kg m^2/s '
        CALL set_label(label,base_label,strlen,'angmpf')
        CALL netcdf_err( nf90_def_var(ncid, "angmpf", nf90_double,    &
             DIMIDS = (/idim_beam_rho,idim_ke_bm,idim_nbeams/),VARID=id_angmpf))
        CALL netcdf_err( nf90_put_att(ncid,id_angmpf,'long_name',label ))
        CALL netcdf_err(  nf90_put_att(ncid,id_angmpf,'units',       &
                          'kg meter^2/second^2) '))
 

        base_label = '* neutral beam flux avg parallel beam mometum kg m/s'
        CALL set_label(label,base_label,strlen,'pb0')
        CALL netcdf_err( nf90_def_var(ncid, "pb0", nf90_double,    &
             DIMIDS = (/idim_beam_rho,idim_ke_bm,idim_nbeams/),VARID=id_pb0))
        CALL netcdf_err( nf90_put_att(ncid,id_pb0,'long_name',label ))
        CALL netcdf_err(  nf90_put_att(ncid,id_pb0,'units',       &
                          'kg*meter/second '))


        base_label = '* neutral beam normalized hot ion birth rate'
        CALL set_label(label,base_label,strlen,'hibr')
        CALL netcdf_err( nf90_def_var(ncid, "hibr", nf90_double,    &
             DIMIDS = (/idim_beam_rho,idim_ke_bm,idim_nbeams/),VARID=id_hibr))
        CALL netcdf_err( nf90_put_att(ncid,id_hibr,'long_name',label ))
        CALL netcdf_err(  nf90_put_att(ncid,id_hibr,'units',       &
                          dimensionless))


        base_label = '* neutral beam normalized hot ion depostion rate'
        CALL set_label(label,base_label,strlen,'hdep')
        CALL netcdf_err( nf90_def_var(ncid, "hdep", nf90_double,    &
             DIMIDS = (/idim_beam_rho,idim_ke_bm,idim_nbeams/),VARID=id_hdep))
        CALL netcdf_err( nf90_put_att(ncid,id_hdep,'long_name',label ))
        CALL netcdf_err(  nf90_put_att(ncid,id_hdep,'units',       &
                          dimensionless))

        base_label = '* neutral beam fast ion pitch angle cosine'
        CALL set_label(label,base_label,strlen,'zeta')
        CALL netcdf_err( nf90_def_var(ncid, "zeta", nf90_double,    &
             DIMIDS = (/idim_beam_rho,idim_ke_bm,idim_nbeams/),VARID=id_zeta))
        CALL netcdf_err( nf90_put_att(ncid,id_zeta,'long_name',label ))
        CALL netcdf_err(  nf90_put_att(ncid,id_zeta,'units',       &
                          dimensionless))
 
        base_label = '* hot ion creation mode,fraction of reactions producing electrons'
        CALL set_label(label,base_label,strlen,'hicme')
        CALL netcdf_err( nf90_def_var(ncid, "hicme", nf90_double,    &
             DIMIDS = (/idim_beam_rho,idim_ke_bm,idim_nbeams/),VARID=id_hicme))
        CALL netcdf_err( nf90_put_att(ncid,id_hicme,'long_name',label ))
!        CALL netcdf_err(  nf90_put_att(ncid,id_hicme,'units',       &
!                          'none '))

        base_label = '* hot ion creation mode,fraction of reactions producing ions species 1'
        CALL set_label(label,base_label,strlen,'hicmp1')
        CALL netcdf_err( nf90_def_var(ncid, "hicmp1", nf90_double,    &
             DIMIDS = (/idim_beam_rho,idim_ke_bm,idim_nbeams/),VARID=id_hicmp1))
        CALL netcdf_err( nf90_put_att(ncid,id_hicmp1,'long_name',label ))
!        CALL netcdf_err(  nf90_put_att(ncid,id_hicmp1,'units',       &
!                          'none '))


 
        base_label = '* hot ion creation mode,fraction of reactions producing ions species 2'
        CALL set_label(label,base_label,strlen,'hicmp2')
        CALL netcdf_err( nf90_def_var(ncid, "hicmp2", nf90_double,    &
             DIMIDS = (/idim_beam_rho,idim_ke_bm,idim_nbeams/),VARID=id_hicmp2))
        CALL netcdf_err( nf90_put_att(ncid,id_hicmp2,'long_name',label ))
!        CALL netcdf_err(  nf90_put_att(ncid,id_hicmp2,'units',       &
!                          'none '))

        base_label = '* rho grid for NB deposition arrays, m'
        CALL set_label(label,base_label,strlen,'rhog_beam')
        CALL netcdf_err( nf90_def_var(ncid, "rhog_beam", nf90_double,    &
             DIMIDS = (/idim_beam_rho/),VARID=id_rhog_beam))
        CALL netcdf_err( nf90_put_att(ncid,id_rhog_beam,'long_name',label ))
        CALL netcdf_err(  nf90_put_att(ncid,id_rhog_beam,'units',       &
                          'm'))
        

     ENDIF nb_loop

     !------------------------------------------------------------
     !define output for common frequencies:
     !------------------------------------------------------------
        k = 0
        DO i = 1,nprim
           k=k+1
           IF(namep(i) == 'dt')THEN
              omega_pi_name = 'omega_pi_D'
              base_label = 'plasma frequency species: D '
              !CALL set_label(label,base_label,strlen,omega_pi_name)
              CALL netcdf_err( nf90_def_var(ncid, omega_pi_name, nf90_double,    &
                   DIMIDS = (/idim_rho/),VARID=id_plas_freq(k)))
              CALL netcdf_err( nf90_put_att(ncid,id_plas_freq(k),'long_name',base_label ))
              CALL netcdf_err(  nf90_put_att(ncid,id_plas_freq(k),'units',       &
                          'rad/sec'))

              omega_ci_name = 'omega_ci_D'
              base_label = 'ion cyclotron frequency species: D '
              CALL netcdf_err( nf90_def_var(ncid, omega_ci_name, nf90_double,    &
                   DIMIDS = (/idim_rho/),VARID=id_ci_freq(k)))
              CALL netcdf_err( nf90_put_att(ncid,id_ci_freq(k),'long_name',base_label ))
              CALL netcdf_err(  nf90_put_att(ncid,id_ci_freq(k),'units',       &
                          'rad/sec'))

              omega_lh_name = 'omega_lh_D'
              base_label = 'lower hybrid  frequency species: D '
              CALL netcdf_err( nf90_def_var(ncid, omega_lh_name, nf90_double,    &
                   DIMIDS = (/idim_rho/),VARID=id_lh_freq(k)))
              CALL netcdf_err( nf90_put_att(ncid,id_lh_freq(k),'long_name',base_label ))
              CALL netcdf_err(  nf90_put_att(ncid,id_lh_freq(k),'units',       &
                          'rad/sec'))

              omega_uh_name = 'omega_uh_D'
              base_label = 'upper hybrid  frequency species: D '
              CALL netcdf_err( nf90_def_var(ncid, omega_uh_name, nf90_double,    &
                   DIMIDS = (/idim_rho/),VARID=id_uh_freq(k)))
              CALL netcdf_err( nf90_put_att(ncid,id_uh_freq(k),'long_name',base_label ))
              CALL netcdf_err(  nf90_put_att(ncid,id_uh_freq(k),'units',       &
                          'rad/sec'))



              k = k+1
              omega_pi_name = 'omega_pi_T'
              base_label = '* plasma frequency species : T'
              CALL set_label(label,base_label,strlen,omega_pi_name)
              CALL netcdf_err( nf90_def_var(ncid, omega_pi_name, nf90_double,    &
                   DIMIDS = (/idim_rho/),VARID=id_plas_freq(k)))
              CALL netcdf_err( nf90_put_att(ncid,id_plas_freq(k),'long_name',base_label ))
              CALL netcdf_err(  nf90_put_att(ncid,id_plas_freq(k),'units',       &
                          'rad/sec'))

              omega_ci_name = 'omega_ci_T'
              base_label = 'ion cyclotron frequency species: T '
              CALL netcdf_err( nf90_def_var(ncid, omega_ci_name, nf90_double,    &
                   DIMIDS = (/idim_rho/),VARID=id_ci_freq(k)))
              CALL netcdf_err( nf90_put_att(ncid,id_ci_freq(k),'long_name',base_label ))
              CALL netcdf_err(  nf90_put_att(ncid,id_ci_freq(k),'units',       &
                          'rad/sec'))

              omega_lh_name = 'omega_lh_T'
              base_label = 'lower hybrid frequency species: T '
              CALL netcdf_err( nf90_def_var(ncid, omega_lh_name, nf90_double,    &
                   DIMIDS = (/idim_rho/),VARID=id_lh_freq(k)))
              CALL netcdf_err( nf90_put_att(ncid,id_lh_freq(k),'long_name',base_label ))
              CALL netcdf_err(  nf90_put_att(ncid,id_lh_freq(k),'units',       &
                          'rad/sec'))

              omega_uh_name = 'omega_uh_T'
              base_label = 'upper hybrid frequency species: T '
              CALL netcdf_err( nf90_def_var(ncid, omega_uh_name, nf90_double,    &
                   DIMIDS = (/idim_rho/),VARID=id_uh_freq(k)))
              CALL netcdf_err( nf90_put_att(ncid,id_uh_freq(k),'long_name',base_label ))
              CALL netcdf_err(  nf90_put_att(ncid,id_uh_freq(k),'units',       &
                          'rad/sec'))

           ELSE

              omega_pi_name = 'omega_pi_'//namep(i)
              base_label = '* plasma frequency species :'//namep(i)
              !CALL set_label(label,base_label,strlen,omega_pi_name)
              CALL netcdf_err( nf90_def_var(ncid,omega_pi_name, nf90_double,    &
                   DIMIDS = (/idim_rho/),VARID=id_plas_freq(k)))
              CALL netcdf_err( nf90_put_att(ncid,id_plas_freq(k),'long_name',base_label ))
              CALL netcdf_err(  nf90_put_att(ncid,id_plas_freq(k),'units',       &
                          'rad/sec'))

              omega_ci_name = 'omega_ci_'//namep(i)
              base_label = '* ion cyclotron  frequency species :'//namep(i)
              !CALL set_label(label,base_label,strlen,omega_pi_name)
              CALL netcdf_err( nf90_def_var(ncid,omega_ci_name, nf90_double,    &
                   DIMIDS = (/idim_rho/),VARID=id_ci_freq(k)))
              CALL netcdf_err( nf90_put_att(ncid,id_ci_freq(k),'long_name',base_label ))
              CALL netcdf_err(  nf90_put_att(ncid,id_ci_freq(k),'units',       &
                          'rad/sec'))

              omega_lh_name = 'omega_lh_'//namep(i)
              base_label = 'lower hybrid  frequency species :'//namep(i)
              !CALL set_label(label,base_label,strlen,omega_pi_name)
              CALL netcdf_err( nf90_def_var(ncid,omega_lh_name, nf90_double,    &
                   DIMIDS = (/idim_rho/),VARID=id_lh_freq(k)))
              CALL netcdf_err( nf90_put_att(ncid,id_lh_freq(k),'long_name',base_label ))
              CALL netcdf_err(  nf90_put_att(ncid,id_lh_freq(k),'units',       &
                          'rad/sec'))

              omega_uh_name = 'omega_uh_'//namep(i)
              base_label = 'upper hybrid  frequency species :'//namep(i)
              !CALL set_label(label,base_label,strlen,omega_uh_name)
              CALL netcdf_err( nf90_def_var(ncid,omega_uh_name, nf90_double,    &
                   DIMIDS = (/idim_rho/),VARID=id_uh_freq(k)))
              CALL netcdf_err( nf90_put_att(ncid,id_uh_freq(k),'long_name',base_label ))
              CALL netcdf_err(  nf90_put_att(ncid,id_uh_freq(k),'units',       &
                          'rad/sec'))
           ENDIF
        ENDDO

        omega_ce_name = 'omega_ce'
        base_label = 'electron cyclotron  frequency'
        !CALL set_label(label,base_label,strlen,omega_ce_name)
        CALL netcdf_err( nf90_def_var(ncid,omega_ce_name, nf90_double,    &
                   DIMIDS = (/idim_rho/),VARID=id_ce_freq))
        CALL netcdf_err( nf90_put_att(ncid,id_ce_freq,'long_name',base_label ))
        CALL netcdf_err(  nf90_put_att(ncid,id_ce_freq,'units',       &
                          'rad/sec'))

        omega_pe_name = 'omega_pe'
        base_label = 'electron plasma frequency'
        !CALL set_label(label,base_label,strlen,omega_pe_name)
        CALL netcdf_err( nf90_def_var(ncid,omega_pe_name, nf90_double,    &
                   DIMIDS = (/idim_rho/),VARID=id_pe_freq))
        CALL netcdf_err( nf90_put_att(ncid,id_pe_freq,'long_name',base_label ))
        CALL netcdf_err(  nf90_put_att(ncid,id_pe_freq,'units',       &
                          'rad/sec'))





    !---------------------------------------------------------------
     CALL netcdf_err(nf90_enddef(ncid))   ! leave netcdf define mode
    !---------------------------------------------------------------


     !         put values of variables

     profile%te0     = get_element(profile%te,1)
     profile%ti0     = get_element(profile%ti,1)

     CALL netcdf_err( nf90_put_var(ncid,id_shot,shot_id%shot_nmbr),id_shot)

     CALL netcdf_err( nf90_put_var(ncid,id_time,time),id_time)

     CALL netcdf_err( nf90_put_var(ncid,id_tgcnmf,tGCNMF),id_tgcnmf)

     CALL netcdf_err( nf90_put_var(ncid,id_time_bc,time_bc),id_time_bc)

     CALL netcdf_err( nf90_put_var(ncid,id_nj,nj),id_nj)

     CALL netcdf_err(nf90_put_var(ncid,id_psiaxis,mhd_dat%psiaxis),id_psiaxis)

     CALL netcdf_err(nf90_put_var(ncid,id_psibdry,mhd_dat%psibdry),id_psibdry)

     CALL netcdf_err(nf90_put_var(ncid,id_rgeom, dischg%rgeom),id_rgeom)

     CALL netcdf_err(nf90_put_var(ncid,id_btgeom, dischg%btgeom),id_btgeom)

     CALL netcdf_err(nf90_put_var(ncid,id_rma, dischg%rma),id_rma)

     CALL netcdf_err(nf90_put_var(ncid,id_zma, dischg%zma),id_zma)

     !         CALL netcdf_err(nf90_put_var(ncid,id_rmag, dischg%rmag),id_rmag)

     CALL netcdf_err(nf90_put_var(ncid,id_rmajor, dischg%rmajor),id_rmajor)

     CALL netcdf_err(nf90_put_var(ncid,id_rsep, dischg%rsep),id_rsep)
     CALL netcdf_err(nf90_put_var(ncid,id_zsep, dischg%zsep),id_zsep)

     CALL netcdf_err(nf90_put_var(ncid,id_rplasmin, dischg%rplasmin),id_rplasmin)
     CALL netcdf_err(nf90_put_var(ncid,id_rplasmax, dischg%rplasmax),id_rplasmax)
     CALL netcdf_err(nf90_put_var(ncid,id_zplasmin, dischg%zplasmin),id_zplasmin)
     CALL netcdf_err(nf90_put_var(ncid,id_zplasmax, dischg%zplasmax),id_zplasmax)


     CALL netcdf_err(nf90_put_var(ncid,id_kappa, dischg%kappa),id_kappa)

     CALL netcdf_err(nf90_put_var(ncid,id_deltao, dischg%deltao),id_deltao)

     CALL netcdf_err(nf90_put_var(ncid,id_pindento, dischg%pindento),id_pindento)

     CALL netcdf_err(nf90_put_var(ncid,id_volume, dischg%volo),id_volume)

     CALL netcdf_err(nf90_put_var(ncid,id_circum, dischg%circum),id_circum)

     CALL netcdf_err(nf90_put_var(ncid,id_areao, dischg%areao),id_areao)





     CALL netcdf_err( nf90_put_var(ncid,id_nion,nion),id_nion)

     CALL netcdf_err( nf90_put_var(ncid,id_nprim,nprim),id_nprim)

     CALL netcdf_err( nf90_put_var(ncid,id_fd_thermal,fd_thermal),id_fd_thermal)

     CALL netcdf_err( nf90_put_var(ncid,id_nimp,nimp),id_nimp)

     CALL netcdf_err( nf90_put_var(ncid,id_nneu,nneu),id_nneu)

     CALL netcdf_err( nf90_put_var(ncid,id_nbion,nbion),id_nbion)

     CALL netcdf_err( nf90_put_var(ncid,id_fd_beam,fd_beam),id_fd_beam)

     !CONVENTION name in netcdf data set is nbeams, value of nbeams is  neut_beam%nbeams

 
     CALL netcdf_err( nf90_put_var(ncid,id_nbeams,neut_beam%nbeams),id_nbeams)

     CALL netcdf_err( nf90_put_var(ncid,id_namep,namep),id_namep)

     CALL netcdf_err( nf90_put_var(ncid,id_namei,namei),id_namei)

     IF ( nneu > 0 ) THEN
     
        CALL netcdf_err( nf90_put_var(ncid,id_namen,namen ),id_namen)
     
     ENDIF

     CALL netcdf_err( nf90_put_var(ncid,id_nameb,nameb ),id_nameb)

     CALL netcdf_err( nf90_put_var(ncid,id_namepel,pellet%name),id_namepel)

     CALL netcdf_err( nf90_put_var(ncid,id_btor,   mhd_dat%btor ),id_btor)
     CALL netcdf_err( nf90_put_var(ncid,id_tot_cur,mhd_dat%tot_cur ),id_tot_cur)
     CALL netcdf_err( nf90_put_var(ncid,id_totohm_cur,mhd_dat%totohm_cur ),id_totohm_cur)
     CALL netcdf_err( nf90_put_var(ncid,id_totboot_cur,mhd_dat%totboot_cur ),id_totboot_cur)
     CALL netcdf_err( nf90_put_var(ncid,id_totbeam_cur,mhd_dat%totbeam_cur ),id_totbeam_cur)
     CALL netcdf_err( nf90_put_var(ncid,id_totrf_cur,  mhd_dat%totrf_cur ),id_totrf_cur)
     CALL netcdf_err( nf90_put_var(ncid,id_betap,  mhd_dat%betap ),id_betap)
     CALL netcdf_err( nf90_put_var(ncid,id_beta,  mhd_dat%beta ),id_beta)
     CALL netcdf_err( nf90_put_var(ncid,id_ali,    mhd_dat%ali ),id_ali)

     CALL netcdf_err( nf90_put_var(ncid,id_te0,profile%te0 ),id_te0)
     CALL netcdf_err( nf90_put_var(ncid,id_ti0,profile%ti0 ),id_ti0)

     CALL netcdf_err( nf90_put_var(ncid,id_pfuse_tot,pfuse_tot ),id_pfuse_tot)

     CALL netcdf_err( nf90_put_var(ncid,id_qrad_tot,qrad_tot ),id_qrad_tot)

     CALL netcdf_err( nf90_put_var(ncid,id_brems_tot,brems_tot ),id_brems_tot)





     ! ---     check that psir_grid is monotonic; it may not be in certain cases where
     ! ---     the current profile was evolved
     !

     work_nj(1:nj) = get_values(psir_grid)
     CALL check_monotonic (work_nj, nj, monotonic, 1)
     IF (.NOT. monotonic) THEN    ! psir_grid is not monotonic
        WRITE  (nlog, 5) time
        WRITE  (ncrt, 5) time
5       FORMAT (' the psi grid, calculated from the poloidal' / &
             ' B field evolution, is not monotonic'        / &
             ' therefore the data at this time (', 1pe12.6, &
             ') was not calculated')
        go to 2000
     END IF

     CALL netcdf_err( nf90_put_var(ncid,id_psir_grid,work_nj ),id_psir_grid)

     CALL netcdf_err( nf90_put_var(ncid,id_psi,mhd_dat%psi ),id_psi)

     work_npsi = get_values(rho_mhd_gridnpsi)
     CALL netcdf_err( nf90_put_var(ncid,id_rho_mhd_gridnpsi,work_npsi ),id_rho_mhd_gridnpsi)

     work_nr = get_values(dischg%rmhdgrid)
     CALL netcdf_err( nf90_put_var(ncid,id_rmhdgrid,work_nr ),id_rmhdgrid)

     work_nz = get_values(dischg%zmhdgrid)
     CALL netcdf_err( nf90_put_var(ncid,id_zmhdgrid,work_nz ),id_zmhdgrid)

     work_nj = get_values(rho_grid)
     CALL netcdf_err( nf90_put_var(ncid,id_rho_grid,work_nj ),id_rho_grid)
 

     work_nj(1:nj) = get_values(mhd_dat%fcap)
     CALL netcdf_err( nf90_put_var(ncid,id_fcap,work_nj ),id_fcap)

     work_nj(1:nj) = get_values(mhd_dat%gcap)
     CALL netcdf_err( nf90_put_var(ncid,id_gcap,work_nj ),id_gcap)

     work_nj(1:nj) = get_values(mhd_dat%hcap)
     CALL netcdf_err( nf90_put_var(ncid,id_hcap,work_nj ),id_hcap)

     work_nj(1:nj) = get_values(mhd_dat%betan)
     CALL netcdf_err( nf90_put_var(ncid,id_betan,work_nj ),id_betan)

     work_nj(1:nj) = get_values(profile%te)
     CALL netcdf_err( nf90_put_var(ncid,id_te,work_nj ),id_te)

     work_nj(1:nj) = get_values(profile%ti)
     CALL netcdf_err( nf90_put_var(ncid,id_ti,work_nj ),id_ti)

     work_nj(1:nj) = get_values(profile%press)
     CALL netcdf_err( nf90_put_var(ncid,id_press,work_nj ),id_press)



     work_nj(1:nj) = get_values(profile%pressb)
     CALL netcdf_err( nf90_put_var(ncid,id_pressb,work_nj ),id_pressb)

     work_nj = get_values(mhd_dat%q_value)
     CALL netcdf_err( nf90_put_var(ncid,id_q_value,work_nj ),id_q_value)

     work_nj(1:nj) = get_values(profile%ene)
     CALL netcdf_err( nf90_put_var(ncid,id_ene,work_nj ),id_ene)

     work_nj(1:nj) = get_values(profile%fluxe)
     CALL netcdf_err( nf90_put_var(ncid,id_flux_elct,work_nj ),id_flux_elct)

     work_nj(1:nj) = get_values(profile%fluxi)
     CALL netcdf_err( nf90_put_var(ncid,id_flux_ion,work_nj ),id_flux_ion)

     DO j=1,nion
        work_nj_nion(1:nj,j) = get_values(profile%en(j))
     ENDDO
     !the species are ordered as 1..nprim,1..nimp
     !where nion = nprim+nimp
     CALL netcdf_err( nf90_put_var(ncid,id_enion,work_nj_nion ),id_enion)


     DO j=1,nion
        work_nj_nion(1:nj,j) = get_values(diffuse%xnus(j))
     ENDDO
     !the species are ordered as 1..nprim,1..nimp
     !where nion = nprim+nimp
     CALL netcdf_err( nf90_put_var(ncid,id_xnus,work_nj_nion ),id_xnus)

     work_nj(:) = get_values(diffuse%xnuse)
     CALL netcdf_err( nf90_put_var(ncid,id_xnuse,work_nj),id_xnuse)

     work_nj(:) = get_values(diffuse%ftrap)
     CALL netcdf_err( nf90_put_var(ncid,id_ftrap,work_nj),id_ftrap)



    work_nj(:) = get_values(diffuse%eta) 

     CALL netcdf_err( nf90_put_var(ncid,id_eta,work_nj),id_eta)

     work_nj(:)  = get_values(diffuse%chie_paleo)

     CALL netcdf_err( nf90_put_var(ncid,id_chiepc,work_nj),id_chiepc)



!    work_nj(:) = get_values(fus_prod%neutr_ddn_th) this assignment causes seg fault
     do j=1,nj
       work_nj(j) = fus_prod%neutr_ddn_th%data(j)
     enddo
     CALL netcdf_err( nf90_put_var(ncid,id_neutr_ddn_th,work_nj),id_neutr_ddn_th)
     work_nj(:) = get_values(fus_prod%neutr_ddn_beam_beam)
     CALL netcdf_err( nf90_put_var(ncid,id_neutr_ddn_beam_beam,work_nj),id_neutr_ddn_beam_beam)

     work_nj(:) = get_values(fus_prod%neutr_ddn_beam_thermal)
     CALL netcdf_err( nf90_put_var(ncid,id_neutr_ddn_beam_thermal,work_nj),id_neutr_ddn_beam_thermal)

     work_nj(:) = get_values(fus_prod%neutr_ddn_knock) 
     CALL netcdf_err( nf90_put_var(ncid,id_neutr_ddn_knock,work_nj),id_neutr_ddn_knock)


     work_nj(:) = get_values(fus_prod%neutr_ddn_tot)
     CALL netcdf_err( nf90_put_var(ncid,id_neutr_ddn_tot,work_nj),id_neutr_ddn_tot)
     CALL netcdf_err( nf90_put_var(ncid,id_total_neutr_ddn_th,fus_prod%total_neutr_ddn_th),id_total_neutr_ddn_th)
     CALL netcdf_err( nf90_put_var(ncid,id_total_neutr_ddn_beam_beam,fus_prod%total_neutr_ddn_beam_beam), &
                                        id_total_neutr_ddn_beam_beam)

     CALL netcdf_err( nf90_put_var(ncid,id_total_neutr_ddn_beam_thermal,fus_prod%total_neutr_ddn_beam_thermal), &
                                        id_total_neutr_ddn_beam_thermal)
     CALL netcdf_err( nf90_put_var(ncid,id_total_neutr_ddn_knock,fus_prod%total_neutr_ddn_knock),id_total_neutr_ddn_knock)
     CALL netcdf_err( nf90_put_var(ncid,id_total_neutr_ddn,fus_prod%total_neutr_ddn),id_total_neutr_ddn)

! some codes may not have these allocated:
        IF(.NOT.  ASSOCIATED(fus_prod%p_ddt%data))THEN 
           fus_prod%p_ddt = zero_vector(nj)
           fus_prod%total_p_ddt = zeroc

        ENDIF
        IF(.NOT. ASSOCIATED(fus_prod%p_he3d_th%data))THEN 
           fus_prod%p_he3d_th = zero_vector(nj)
           fus_prod%total_p_he3d_th = zeroc
           fus_prod%total_p_he3d = zeroc
        ENDIF
        IF(.NOT. ASSOCIATED(fus_prod%p_he3d_beam_thermal%data))THEN 
           fus_prod%p_he3d_beam_thermal = zero_vector(nj)
           fus_prod%total_p_he3d_beam_thermal = zeroc
           fus_prod%total_p_he3d = zeroc
        ENDIF
 

     work_nj(:) = get_values(fus_prod%p_ddt)

     CALL netcdf_err( nf90_put_var(ncid,id_p_ddt,work_nj),id_p_ddt)

     CALL netcdf_err( nf90_put_var(ncid,id_total_p_ddt,fus_prod%total_p_ddt),id_total_p_ddt)

     work_nj(:) = get_values(fus_prod%p_he3d_th)
     CALL netcdf_err( nf90_put_var(ncid,id_p_he3d_th,work_nj),id_p_he3d_th)
     CALL netcdf_err( nf90_put_var(ncid,id_total_p_he3d_th,fus_prod%total_p_he3d_th),id_total_p_he3d_th)

     work_nj(:) = get_values(fus_prod%p_he3d_beam_thermal)
     CALL netcdf_err( nf90_put_var(ncid,id_p_he3d_beam_th,work_nj),id_p_he3d_beam_th)
     CALL netcdf_err( nf90_put_var(ncid,id_total_p_he3d_beam_th,fus_prod%total_p_he3d_beam_thermal),id_total_p_he3d_beam_th)
     CALL netcdf_err( nf90_put_var(ncid,id_total_p_he3d,fus_prod%total_p_he3d),id_total_p_he3d)
     CALL netcdf_err( nf90_put_var(ncid,id_he3_frac,he3_frac),id_he3_frac)
     CALL netcdf_err( nf90_put_var(ncid,id_he3_thermal_spin_pol,he3_thermal_spin_pol),id_he3_thermal_spin_pol)
     CALL netcdf_err( nf90_put_var(ncid,id_d_beam_spin_pol,neut_beam%d_beam_spin_pol),id_d_beam_spin_pol)



     work_nj(:) = glf_p_output(:,1)
     CALL netcdf_err( nf90_put_var(ncid,id_glf_elct_partcl_flux,work_nj),id_glf_elct_partcl_flux)
     work_nj(:) = glf_p_output(:,2)
     CALL netcdf_err( nf90_put_var(ncid,id_glf_primion_partcl_flux,work_nj),id_glf_primion_partcl_flux)
     work_nj(:) = glf_p_output(:,3)
     CALL netcdf_err( nf90_put_var(ncid,id_glf_impion_partcl_flux,work_nj),id_glf_impion_partcl_flux)


 
 

    !---------------------------------------------------------------------
    ! multimode output (put zeros if multimode was not used)
    !---------------------------------------------------------------------
     IF(.NOT. ASSOCIATED(diffuse%mmm_gammaDBM%data) &
          .OR. diffuse%mmm_gammaDBM%size .NE. nj  )THEN
        ! when running glf or tglf this array may not be the right size
        ! if grid size changes are involved. 
        diffuse%mmm_gammaDBM = zero_vector(nj)
     ENDIF
      work_nj(:) = diffuse%mmm_gammaDBM%data(:)
      CALL netcdf_err(nf90_put_var(ncid,id_mmm_gammaDBM,work_nj),id_mmm_gammaDBM)

     IF(.NOT. ASSOCIATED(diffuse%mmm_omegaDBM%data) &
         .OR. diffuse%mmm_omegaDBM%size .NE. nj  )THEN
        ! when running glf or tglf this array may not be the right size
        ! if grid size changes are involved. 
        diffuse%mmm_omegaDBM = zero_vector(nj)
     ENDIF
      work_nj(:) = diffuse%mmm_omegaDBM%data(:)
      CALL netcdf_err(nf90_put_var(ncid,id_mmm_omegaDBM,work_nj),id_mmm_gammaDBM)
     IF(.NOT. ASSOCIATED(diffuse%mmm_xdi%data) &
          .OR. diffuse%mmm_xdi%size .ne. nj)THEN
        diffuse%mmm_xdi = zero_vector(nj)
     ENDIF
      work_nj(:) = diffuse%mmm_xdi%data(:)
      CALL netcdf_err(nf90_put_var(ncid,id_mmm_xdi,work_nj),id_mmm_xdi)

     IF(.NOT. ASSOCIATED(diffuse%mmm_xti%data) &
          .OR. diffuse%mmm_xti%size .ne. nj)THEN
        diffuse%mmm_xti = zero_vector(nj)
     ENDIF
      work_nj(:) = diffuse%mmm_xti%data(:)
      CALL netcdf_err(nf90_put_var(ncid,id_mmm_xti,work_nj),id_mmm_xti)

     IF(.NOT. ASSOCIATED(diffuse%mmm_xte%data) &
          .OR. diffuse%mmm_xte%size .ne. nj)THEN
        diffuse%mmm_xte = zero_vector(nj)
     ENDIF
      work_nj(:) = diffuse%mmm_xte%data(:)
      CALL netcdf_err(nf90_put_var(ncid,id_mmm_xte,work_nj),id_mmm_xte)

     IF(.NOT. ASSOCIATED(diffuse%mmm_xdz%data) &
          .OR. diffuse%mmm_xdz%size .ne. nj)THEN
        diffuse%mmm_xdz = zero_vector(nj)
     ENDIF
      work_nj(:) = diffuse%mmm_xdz%data(:)
      CALL netcdf_err(nf90_put_var(ncid,id_mmm_xdz,work_nj),id_mmm_xdz)

     IF(.NOT. ASSOCIATED(diffuse%mmm_xvt%data) &
        .OR. diffuse%mmm_xvt%size .ne. nj)THEN
        diffuse%mmm_xvt = zero_vector(nj)
     ENDIF
      work_nj(:) = diffuse%mmm_xvt%data(:)
      CALL netcdf_err(nf90_put_var(ncid,id_mmm_xvt,work_nj),id_mmm_xvt)

     IF(.NOT. ASSOCIATED(diffuse%mmm_xvp%data) &
          .OR. diffuse%mmm_xvp%size .ne. nj)THEN
        diffuse%mmm_xvp = zero_vector(nj)
     ENDIF
      work_nj(:) = diffuse%mmm_xvp%data(:)
      CALL netcdf_err(nf90_put_var(ncid,id_mmm_xvp,work_nj),id_mmm_xvp)

     IF(.NOT. ASSOCIATED(diffuse%mmm_xtiW20%data) &
        .OR. diffuse%mmm_xtiw20%size .ne. nj)THEN
        diffuse%mmm_xtiW20 = zero_vector(nj)
     ENDIF
      work_nj(:) = diffuse%mmm_xtiW20%data(:)
      CALL netcdf_err(nf90_put_var(ncid,id_mmm_xtiW20,work_nj),id_mmm_xtiW20)

     IF(.NOT. ASSOCIATED(diffuse%mmm_xdiW20%data) &
          .OR. diffuse%mmm_xdiw20%size .ne. nj)THEN
        diffuse%mmm_xdiW20 = zero_vector(nj)
     ENDIF
      work_nj(:) = diffuse%mmm_xdiW20%data(:)
      CALL netcdf_err(nf90_put_var(ncid,id_mmm_xdiW20,work_nj),id_mmm_xdiW20)

     IF(.NOT. ASSOCIATED(diffuse%mmm_xteW20%data) &
          .OR. diffuse%mmm_xteW20%size .ne. nj)THEN
        diffuse%mmm_xteW20 = zero_vector(nj)
     ENDIF
      work_nj(:) = diffuse%mmm_xteW20%data(:)
      CALL netcdf_err(nf90_put_var(ncid,id_mmm_xteW20,work_nj),id_mmm_xteW20)

     IF(.NOT. ASSOCIATED(diffuse%mmm_xtiDBM%data) &
          .OR. diffuse%mmm_xtiDBM%size .ne. nj)THEN
        diffuse%mmm_xtiDBM = zero_vector(nj)
     ENDIF
      work_nj(:) = diffuse%mmm_xtiDBM%data(:)
      CALL netcdf_err(nf90_put_var(ncid,id_mmm_xtiDBM,work_nj),id_mmm_xtiDBM)

     IF(.NOT. ASSOCIATED(diffuse%mmm_xdiDBM%data) &
          .OR. diffuse%mmm_xdiDBM%size .ne. nj)THEN
        diffuse%mmm_xdiDBM = zero_vector(nj)
     ENDIF
      work_nj(:) = diffuse%mmm_xdiDBM%data(:)
      CALL netcdf_err(nf90_put_var(ncid,id_mmm_xdiDBM,work_nj),id_mmm_xdiDBM)

     IF(.NOT. ASSOCIATED(diffuse%mmm_xteDBM%data) &
          .OR. diffuse%mmm_xteDBM%size .ne. nj)THEN
        diffuse%mmm_xteDBM = zero_vector(nj)
     ENDIF
      work_nj(:) = diffuse%mmm_xteDBM%data(:)
      CALL netcdf_err(nf90_put_var(ncid,id_mmm_xteDBM,work_nj),id_mmm_xteDBM)

     IF(.NOT. ASSOCIATED(diffuse%mmm_xteETG%data) &
          .OR. diffuse%mmm_xteETG%size .ne. nj)THEN
        diffuse%mmm_xteETG = zero_vector(nj)
     ENDIF
      work_nj(:) = diffuse%mmm_xteETG%data(:)
      CALL netcdf_err(nf90_put_var(ncid,id_mmm_xteETG,work_nj),id_mmm_xteETG)


      ngW20 = 4
      IF( .NOT. ASSOCIATED(diffuse%mmm_gammaW20))THEN
          ALLOCATE(diffuse%mmm_gammaW20(ngW20))   
       ELSEIF(diffuse%mmm_gammaW20(1)%size .ne. nj )THEN
          DO jj =1,ngw20
             CALL delete_Vector(diffuse%mmm_gammaW20(jj))
          ENDDO
       ENDIF
       IF(diffuse%mmm_gammaW20(1)%size .ne. nj )THEN
          DO jj =1,ngW20
             diffuse%mmm_gammaW20(jj)  = zero_vector(nj)
          ENDDO
       ENDIF
      DO jj = 1,ngw20
         work_nj(:) = diffuse%mmm_gammaW20(jj)%data(:)
         IF(jj == 1) &
         CALL netcdf_err(nf90_put_var(ncid,id_mmm_gamma_i1_W20,work_nj),id_mmm_gamma_i1_W20)
         IF(jj == 2) &
         CALL netcdf_err(nf90_put_var(ncid,id_mmm_gamma_e1_W20,work_nj),id_mmm_gamma_e1_W20)
         IF(jj == 3) &
         CALL netcdf_err(nf90_put_var(ncid,id_mmm_gamma_i2_W20,work_nj),id_mmm_gamma_i2_W20)
         IF(jj == 4) &
         CALL netcdf_err(nf90_put_var(ncid,id_mmm_gamma_e2_W20,work_nj),id_mmm_gamma_e2_W20)
      ENDDO

      ngW20 = 4
      IF( .NOT. ASSOCIATED(diffuse%mmm_omegaW20))THEN
         ALLOCATE(diffuse%mmm_omegaW20(ngW20))  
      ELSEIF(diffuse%mmm_omegaW20(1)%size .ne. nj )THEN
         DO jj =1,ngw20
            CALL delete_Vector(diffuse%mmm_omegaW20(jj))
         ENDDO
      ENDIF
      IF(diffuse%mmm_omegaW20(1)%size .ne. nj )THEN
         DO jj =1,ngW20
            diffuse%mmm_omegaW20(jj)  = zero_vector(nj)
         ENDDO
      ENDIF
      DO jj = 1,ngw20
         work_nj(:) = diffuse%mmm_omegaW20(jj)%data(:)
         IF(jj == 1) &
         CALL netcdf_err(nf90_put_var(ncid,id_mmm_omega_i1_W20,work_nj),id_mmm_omega_i1_W20)
         IF(jj == 2) &
         CALL netcdf_err(nf90_put_var(ncid,id_mmm_omega_e1_W20,work_nj),id_mmm_omega_e1_W20)
         IF(jj == 3) &
         CALL netcdf_err(nf90_put_var(ncid,id_mmm_omega_i2_W20,work_nj),id_mmm_omega_i2_W20)
         IF(jj == 4) &
         CALL netcdf_err(nf90_put_var(ncid,id_mmm_omega_e2_W20,work_nj),id_mmm_omega_e2_W20)
      ENDDO



      ngW20 = 4
      IF( .NOT. ASSOCIATED(diffuse%mmm_vflux))THEN
         ALLOCATE(diffuse%mmm_vflux(ngW20))  
      ELSEIF(diffuse%mmm_vflux(1)%size .ne. nj )THEN
         DO jj =1,ngW20
            CALL delete_Vector(diffuse%mmm_vflux(jj))
         ENDDO
      ENDIF
      IF(diffuse%mmm_vflux(1)%size .ne. nj )THEN
         DO jj =1,ngW20
            diffuse%mmm_vflux(jj)  = zero_vector(nj)
         ENDDO
      ENDIF
      DO jj = 1,ngw20
         work_nj(:) = diffuse%mmm_vflux(jj)%data(:)
         IF(jj == 1) &
         CALL netcdf_err(nf90_put_var(ncid,id_mmm_flux_ith,work_nj),id_mmm_flux_ith)
         IF(jj == 2) &
         CALL netcdf_err(nf90_put_var(ncid,id_mmm_flux_ip,work_nj),id_mmm_flux_ip)
         IF(jj == 3) &
         CALL netcdf_err(nf90_put_var(ncid,id_mmm_flux_eth,work_nj),id_mmm_flux_eth)
         IF(jj == 4) &
         CALL netcdf_err(nf90_put_var(ncid,id_mmm_flux_imp,work_nj),id_mmm_flux_imp)
      ENDDO


      ngW20 = 6
      IF( .NOT. ASSOCIATED(diffuse%mmm_vconv))THEN
         ALLOCATE(diffuse%mmm_vconv(ngW20))   
      ELSEIF(diffuse%mmm_vconv(1)%size .ne. nj)THEN 
         DO jj =1,ngw20
            CALL delete_vector(diffuse%mmm_vconv(jj))
         ENDDO
      ENDIF
      IF(diffuse%mmm_vconv(1)%size .ne. nj)THEN
         DO jj =1,ngW20
            diffuse%mmm_vconv(jj)  = zero_vector(nj)
         ENDDO
      ENDIF

      DO jj = 1,ngw20
  !NOTE jj=1,2,3 info is put out as zero by current version of multimode
         work_nj(:) = diffuse%mmm_vconv(jj)%data(:)
         IF(jj == 1) &
         CALL netcdf_err(nf90_put_var(ncid,id_mmm_vconv_ith,work_nj),id_mmm_vconv_ith)
         IF(jj == 2) &
         CALL netcdf_err(nf90_put_var(ncid,id_mmm_vconv_ip,work_nj),id_mmm_vconv_ip)
         IF(jj == 3) &
         CALL netcdf_err(nf90_put_var(ncid,id_mmm_vconv_eth,work_nj),id_mmm_vconv_eth)
         IF(jj == 4) &
         CALL netcdf_err(nf90_put_var(ncid,id_mmm_vconv_imp,work_nj),id_mmm_vconv_imp)
         IF(jj == 5) &
         CALL netcdf_err(nf90_put_var(ncid,id_mmm_vmtmt,work_nj),id_mmm_vmtmt)
         IF(jj == 6) &
         CALL netcdf_err(nf90_put_var(ncid,id_mmm_vmtmp,work_nj),id_mmm_vmtmp)

      ENDDO


    ! end multimode output (put zeros if multimode was not used)


     DO j=1,nion
        work_nj_nion(1:nj,j) = get_values(profile%flux(j))      ! corrected 7/29/13 HSJ
     ENDDO
        CALL netcdf_err( nf90_put_var(ncid,id_pflux,work_nj_nion),id_pflux)

     DO j=1,nion 
        work_nj_nion(1:nj,j) = get_values(profile%flux_conv(j)) ! corrected 7/29/13 HSJ
     ENDDO
        CALL netcdf_err( nf90_put_var(ncid,id_pflux_conv,work_nj_nion),id_pflux_conv)

 

     work_nj(:) = get_values(profile%flux(nion+1)) 
     CALL netcdf_err( nf90_put_var(ncid,id_e_fluxe,work_nj),id_e_fluxe)

     work_nj(:) = get_values(profile%flux_conv(nion+1)) 
     CALL netcdf_err( nf90_put_var(ncid,id_e_fluxe_conv,work_nj),id_e_fluxe_conv)

     work_nj(:) = glf_e_output(:,1)
     CALL netcdf_err( nf90_put_var(ncid,id_glf_elct_eng_flux,work_nj),id_glf_elct_eng_flux)
     work_nj(:) = glf_e_output(:,2)
     CALL netcdf_err( nf90_put_var(ncid,id_glf_primion_eng_flux,work_nj),id_glf_primion_eng_flux)
     work_nj(:) = glf_e_output(:,3)
     CALL netcdf_err( nf90_put_var(ncid,id_glf_impion_eng_flux,work_nj),id_glf_impion_eng_flux)
     work_nj(:) = glf_m_output(:,1)
     CALL netcdf_err( nf90_put_var(ncid,id_glf_elct_momtm_flux,work_nj),id_glf_elct_momtm_flux)

     work_nj(:) = glf_m_output(:,2)
     CALL netcdf_err( nf90_put_var(ncid,id_glf_primion_momtm_flux,work_nj),id_glf_primion_momtm_flux)

     work_nj(:) = glf_m_output(:,3)
     CALL netcdf_err( nf90_put_var(ncid,id_glf_impion_momtm_flux,work_nj),id_glf_impion_momtm_flux)

     DO i=1,ntot
        work_nj(:) = glf_eff_chi_output(:,i)
        CALL netcdf_err( nf90_put_var(ncid,id_glf_eff_chi(i),work_nj),id_glf_eff_chi(i))
     ENDDO

     CALL netcdf_err( nf90_put_var(ncid,id_glf_etg_flux,glf_etg_output),id_glf_etg_flux)
     CALL netcdf_err( nf90_put_var(ncid,id_glf_gam_net_i,glf_gamma_net_i_output),id_glf_gam_net_i)
     CALL netcdf_err( nf90_put_var(ncid,id_glf_gam_net_e,glf_gamma_net_e_output),id_glf_gam_net_e)
     CALL netcdf_err( nf90_put_var(ncid,id_glf_anfreq,glf_anfreq_output),id_glf_anfreq)
     CALL netcdf_err( nf90_put_var(ncid,id_glf_anfreq2,glf_anfreq2_output),id_glf_anfreq2)
     CALL netcdf_err( nf90_put_var(ncid,id_glf_anrate,glf_anrate_output),id_glf_anrate)
     CALL netcdf_err( nf90_put_var(ncid,id_glf_anrate2,glf_anrate2_output),id_glf_anrate2)

     work_nj(:) = get_values(profile%flux(nion+2)) 
     CALL netcdf_err( nf90_put_var(ncid,id_e_fluxi,work_nj),id_e_fluxi)

     work_nj(:) = get_values(profile%flux_conv(nion+2)) 
     CALL netcdf_err( nf90_put_var(ncid,id_e_fluxi_conv,work_nj),id_e_fluxi_conv)

     work_nj(:) = get_values(profile%flux(nion+3)) 
     CALL netcdf_err( nf90_put_var(ncid,id_fdyflux,work_nj),id_fdyflux)

     work_nj(:) = get_values(profile%flux_conv(nion+3)) 
     CALL netcdf_err( nf90_put_var(ncid,id_fdyflux_conv,work_nj),id_fdyflux_conv)

     work_nj(:) = get_values(profile%flux(nion+4)) 
     CALL netcdf_err( nf90_put_var(ncid,id_rotflux,work_nj),id_rotflux)

     work_nj(:) = get_values(profile%flux_conv(nion+4)) 
     CALL netcdf_err( nf90_put_var(ncid,id_rotflux_conv,work_nj),id_rotflux_conv)



     work_nj(:) = tglf_p_output(:,1)
     CALL netcdf_err( nf90_put_var(ncid,id_tglf_p_fluxe,work_nj ),id_tglf_p_fluxe)
     work_nj(:) = tglf_p_output(:,2)
     CALL netcdf_err( nf90_put_var(ncid,id_tglf_p_fluxp,work_nj ),id_tglf_p_fluxp)
     work_nj(:) = tglf_p_output(:,3)
     CALL netcdf_err( nf90_put_var(ncid,id_tglf_p_fluxi,work_nj ),id_tglf_p_fluxi)


     work_nj(:) = tglf_e_output(:,1) ! tglf electron energy flux
     CALL netcdf_err( nf90_put_var(ncid,id_tglf_e_fluxe,work_nj ),id_tglf_e_fluxe)
     work_nj(:) = tglf_e_output(:,2) ! tglf effective  primary ion energy  flux
     CALL netcdf_err( nf90_put_var(ncid,id_tglf_e_fluxp,work_nj ),id_tglf_e_fluxp)
     work_nj(:) = tglf_e_output(:,3) ! tglf effective impurity energy flux
     CALL netcdf_err( nf90_put_var(ncid,id_tglf_e_fluxi,work_nj ),id_tglf_e_fluxi)

     work_nj(:) = tglf_m_output(:,1) ! tglf electron momentum flux
     CALL netcdf_err( nf90_put_var(ncid,id_tglf_m_fluxe,work_nj ),id_tglf_m_fluxe)
     work_nj(:) = tglf_m_output(:,2) ! tglf effective  primary ion momentum  flux
     CALL netcdf_err( nf90_put_var(ncid,id_tglf_m_fluxp,work_nj ),id_tglf_m_fluxp)
     work_nj(:) = tglf_m_output(:,3) ! tglf effective impurity momentum flux
     CALL netcdf_err( nf90_put_var(ncid,id_tglf_m_fluxi,work_nj ),id_tglf_m_fluxi)
     
     ! Neutral densities/sources
     IF ( nneu > 0 ) THEN 
        
        CALL netcdf_err( nf90_put_var(ncid,id_enn,enn ),id_enn)

        CALL netcdf_err( nf90_put_var(ncid,id_ennw,ennw ),id_ennw)

        CALL netcdf_err( nf90_put_var(ncid,id_ennv,ennv ),id_ennv)

        CALL netcdf_err( nf90_put_var(ncid,id_volsn,volsn ),id_volsn)

     ENDIF !Neutrals

     work_nj = get_values(prtcl_src%stfuse) 
     CALL netcdf_err( nf90_put_var(ncid,id_stfuse,work_nj),id_stfuse)

     work_nj = get_values(prtcl_src%sbfuse) 
     CALL netcdf_err( nf90_put_var(ncid,id_sbfuse,work_nj),id_sbfuse)

     work_nj = get_values(prtcl_src%spellet) 
     CALL netcdf_err( nf90_put_var(ncid,id_spellet,work_nj),id_spellet)

     CALL netcdf_err( nf90_put_var(ncid,id_sbeame,sbeame),id_sbeame)

     CALL netcdf_err( nf90_put_var(ncid,id_sbeam,sbeam),id_sbeam)

     CALL netcdf_err( nf90_put_var(ncid,id_sion,sion),id_sion)

!     CALL netcdf_err( nf90_put_var(ncid,id_srecom,srecom),id_srecom)
     DO jj=1,nion
        work_nj_nion(:,jj) = get_values(prtcl_src%srecom(jj))
     ENDDO
     CALL netcdf_err( nf90_put_var(ncid,id_srecom,work_nj_nion),id_srecom)


     CALL netcdf_err( nf90_put_var(ncid,id_scx,scx),id_scx)

     CALL netcdf_err( nf90_put_var(ncid,id_sbcx,sbcx),id_sbcx)

     CALL netcdf_err( nf90_put_var(ncid,id_stsource,stsource),id_stsource)

     CALL netcdf_err( nf90_put_var(ncid,id_dudtsv,dudtsv),id_dudtsv)

     CALL netcdf_err( nf90_put_var(ncid,id_enbeam,enbeam),id_enbeam)

     work_nj(1:nj) = get_values(mhd_dat%curden)
     CALL netcdf_err( nf90_put_var(ncid,id_curden,work_nj),id_curden)

     work_nj(1:nj) = get_values(mhd_dat%curpar)
     CALL netcdf_err( nf90_put_var(ncid,id_curpar,work_nj),id_curpar)

     work_nj(1:nj) = get_values(mhd_dat%curohm)
     CALL netcdf_err( nf90_put_var(ncid,id_curohm,work_nj),id_curohm)

     work_nj(1:nj) = get_values(mhd_dat%curboot)
     CALL netcdf_err( nf90_put_var(ncid,id_curboot,work_nj),id_curboot)

     CALL netcdf_err( nf90_put_var(ncid,id_curbeam,curbeam),id_curbeam)

     CALL netcdf_err( nf90_put_var(ncid,id_currf,currf),id_currf)

     work_nj(1:nj) = get_values(profile%etor)
     CALL netcdf_err( nf90_put_var(ncid,id_etor,work_nj),id_etor)

     work_nj(1:nj) = get_values(mhd_dat%rbp)
     CALL netcdf_err( nf90_put_var(ncid,id_rbp,work_nj),id_rbp)

     work_npsi(1:mhd_dat%npsi) = get_values(mhd_dat%ravgnpsi)
     CALL netcdf_err( nf90_put_var(ncid,id_ravgnpsi,work_npsi),id_ravgnpsi)

     work_npsi(1:mhd_dat%npsi) = get_values(mhd_dat%ravginpsi)
     CALL netcdf_err( nf90_put_var(ncid,id_ravginpsi,work_npsi),id_ravginpsi)

     work_npsi(1:mhd_dat%npsi) = get_values(mhd_dat%psivalnpsi)
     CALL netcdf_err( nf90_put_var(ncid,id_psivalnpsi,work_npsi),id_psivalnpsi)

     work_npsi(1:mhd_dat%npsi) = get_values(mhd_dat%fpsinpsi)
     CALL netcdf_err( nf90_put_var(ncid,id_fpsinpsi,work_npsi),id_fpsinpsi)

     work_nj(1:nj) = get_values(mhd_dat%pprim)
     CALL netcdf_err( nf90_put_var(ncid,id_pprim,work_nj),id_pprim)

     work_nj(1:nj) = get_values(mhd_dat%ffprim)
     CALL netcdf_err( nf90_put_var(ncid,id_ffprim,work_nj),id_ffprim)

     work_nj(1:nj) = get_values(mhd_dat%bp)
     CALL netcdf_err( nf90_put_var(ncid,id_bp,work_nj),id_bp)

     work_nj(1:nj) = get_values(mhd_dat%bprmaj)
     CALL netcdf_err( nf90_put_var(ncid,id_bprmaj,work_nj),id_bprmaj)

     work_nj(1:nj) = get_values(mhd_dat%btotrmaj)
     CALL netcdf_err( nf90_put_var(ncid,id_btotrmaj,work_nj),id_btotrmaj)

     work_nj(1:nj) = get_values(profile%angrot)
     CALL netcdf_err( nf90_put_var(ncid,id_angrot,work_nj),id_angrot)

     CALL netcdf_err( nf90_put_var(ncid,id_zeff,zeff),id_zeff)

     IF(.NOT. ASSOCIATED(profile%vpol%data))THEN
        profile%vpol = zero_vector(nj)
     ENDIF
     work_nj(1:nj) = get_values(profile%vpol)
     CALL netcdf_err( nf90_put_var(ncid,id_vpol,work_nj),id_vpol)

     IF(.NOT. ASSOCIATED(profile%vpol_nclass%data))THEN
        profile%vpol_nclass = zero_vector(nj)
     ENDIF
     work_nj(1:nj) = get_values(profile%vpol_nclass)
     CALL netcdf_err( nf90_put_var(ncid,id_vpol_nclass,work_nj),id_vpol_nclass)


     IF(.NOT. ASSOCIATED(profile%vpar%data))THEN
        profile%vpar = zero_vector(nj)
     ENDIF
     work_nj(1:nj) = get_values(profile%vpar)
     CALL netcdf_err( nf90_put_var(ncid,id_vpar,work_nj),id_vpar)

     IF(.NOT. ASSOCIATED(profile%vpar_nclass%data))THEN
        profile%vpar_nclass = zero_vector(nj)
     ENDIF
     work_nj(1:nj) = get_values(profile%vpar_nclass)
     CALL netcdf_err( nf90_put_var(ncid,id_vpar_nclass,work_nj),id_vpar_nclass)


     IF(.NOT. ASSOCIATED(profile%er_tot_nclass%data))THEN
        profile%er_tot_nclass = zero_vector(nj)
     ENDIF
     work_nj(1:nj) = get_values(profile%er_tot_nclass)
     CALL netcdf_err( nf90_put_var(ncid,id_er_tot_nclass,work_nj),id_er_tot_nclass)

     CALL netcdf_err( nf90_put_var(ncid,id_d,diffuse%dcoef),id_d)

     work_nj(1:nj) = get_values(diffuse%chieinv)
     CALL netcdf_err( nf90_put_var(ncid,id_chieinv,work_nj),id_chieinv)


     work_nj(1:nj) = get_values(diffuse%chiinv)
     CALL netcdf_err( nf90_put_var(ncid,id_chiinv,work_nj),id_chiinv)

     work_nj(1:nj) = get_values(diffuse%xkineo)
     CALL netcdf_err( nf90_put_var(ncid,id_xkineo,work_nj),id_xkineo)

     work_nj(1:nj) = get_values(diffuse%xkeneo)
     CALL netcdf_err( nf90_put_var(ncid,id_xkeneo,work_nj),id_xkeneo)


     work_nj(1:nj) = get_values(wpdot%dpedt)
     CALL netcdf_err( nf90_put_var(ncid,id_dpedt,work_nj),id_dpedt)

     DO jj = 1,nion
        work_nj_nion(1:nj,jj) = get_values(wpdot%dpidt(jj))
     ENDDO
     CALL netcdf_err( nf90_put_var(ncid,id_dpidt,work_nj_nion),id_dpidt)

     work_nj(1:nj) = get_values(pwrden%qconde)
     CALL netcdf_err( nf90_put_var(ncid,id_qconde,work_nj),id_qconde)


     work_nj(1:nj) = get_values(pwrden%qcondi)
     CALL netcdf_err( nf90_put_var(ncid,id_qcondi,work_nj),id_qcondi)

     work_nj(1:nj) = get_values(pwrden%qconve)
     CALL netcdf_err( nf90_put_var(ncid,id_qconve,work_nj),id_qconve)


     work_nj(1:nj) = get_values(pwrden%qconvi)
     CALL netcdf_err( nf90_put_var(ncid,id_qconvi,work_nj),id_qconvi)


     work_nj(1:nj) = get_values(pwrden%qbeame)
     CALL netcdf_err( nf90_put_var(ncid,id_qbeame,work_nj),id_qbeame)

     work_nj(1:nj) = get_values(pwrden%qbeami)
     CALL netcdf_err( nf90_put_var(ncid,id_qbeami,work_nj),id_qbeami)


     work_nj(1:nj) = get_values(pwrden%qdelt)
     CALL netcdf_err( nf90_put_var(ncid,id_qdelt,work_nj),id_qdelt)

     work_nj(1:nj) = get_values(pwrden%qexch)
     CALL netcdf_err( nf90_put_var(ncid,id_qexch,work_nj),id_qexch)

     work_nj(1:nj) = get_values(pwrden%qrfe)
     CALL netcdf_err( nf90_put_var(ncid,id_qrfe,work_nj),id_qrfe)

     work_nj(1:nj) = get_values(pwrden%qrfi)
     CALL netcdf_err( nf90_put_var(ncid,id_qrfi,work_nj),id_qrfi)

     work_nj(1:nj) = get_values(pwrden%qione)
     CALL netcdf_err( nf90_put_var(ncid,id_qione,work_nj),id_qione)

     work_nj(1:nj) = get_values(pwrden%qioni)
     CALL netcdf_err( nf90_put_var(ncid,id_qioni,work_nj),id_qioni)

     work_nj(1:nj) = get_values(pwrden%qcx)
     CALL netcdf_err( nf90_put_var(ncid,id_qcx,work_nj),id_qcx)

     work_nj(1:nj) = get_values(pwrden%qe2d)
     CALL netcdf_err( nf90_put_var(ncid,id_qe2d,work_nj),id_qe2d)
     work_nj(1:nj) = get_values(pwrden%qi2d)
     CALL netcdf_err( nf90_put_var(ncid,id_qi2d,work_nj),id_qi2d)


     work_nj(1:nj) = get_values(pwrden%qfuse)
     CALL netcdf_err( nf90_put_var(ncid,id_qfuse,work_nj),id_qfuse)
     work_nj(1:nj) = get_values(pwrden%qfusi)
     CALL netcdf_err( nf90_put_var(ncid,id_qfusi,work_nj),id_qfusi)

     !
     ! --- beam fusion electron ,ion heating profile
     ! --- fraction of beam fusion energy deposited on thermal electron and distributions:
     work_nj(1:nj) = get_values(pwrden%qbfuse)
     CALL netcdf_err( nf90_put_var(ncid,id_qbfuse,work_nj),id_qbfuse)
     work_nj(1:nj) = get_values(pwrden%qbfusi)
     CALL netcdf_err( nf90_put_var(ncid,id_qbfusi,work_nj),id_qbfusi)



     work_nj(1:nj) = get_values(pwrden%qmag)
     CALL netcdf_err( nf90_put_var(ncid,id_qmag,work_nj),id_qmag)
     work_nj(1:nj) = get_values(pwrden%qsawe)
     CALL netcdf_err( nf90_put_var(ncid,id_qsawe,work_nj),id_qsawe)
     work_nj(1:nj) = get_values(pwrden%qsawi)
     CALL netcdf_err( nf90_put_var(ncid,id_qsawi,work_nj),id_qsawi)
     work_nj(1:nj) = get_values(pwrden%qrad)
     CALL netcdf_err( nf90_put_var(ncid,id_qrad,work_nj),id_qrad)


     IF(ASSOCIATED(pwrden%omegale%data))THEN
        work_nj(1:nj) = get_values(pwrden%omegale)
     ELSE
        work_nj(1:nj) =zeroc
     ENDIF
     CALL netcdf_err( nf90_put_var(ncid,id_omegale,work_nj),id_omegale)

     IF(ASSOCIATED(pwrden%qomegapi%data))THEN
        work_nj(1:nj) = get_values(pwrden%qomegapi)
     ELSE
        work_nj(1:nj) =zeroc
     ENDIF
     CALL netcdf_err( nf90_put_var(ncid,id_qomegapi,work_nj),id_qomegapi)


     IF(ASSOCIATED(pwrden%qangce%data))THEN
        work_nj(1:nj) = get_values(pwrden%qangce)
     ELSE
        work_nj(1:nj) =zeroc
     ENDIF
     CALL netcdf_err( nf90_put_var(ncid,id_qangce,work_nj),id_qangce)


     IF(ASSOCIATED(pwrden%sprcxre%data))THEN
        work_nj(1:nj) = get_values(pwrden%sprcxre)
     ELSE
        work_nj(1:nj) =zeroc
     ENDIF
     CALL netcdf_err( nf90_put_var(ncid,id_sprcxre,work_nj),id_sprcxre)

     IF(ASSOCIATED(pwrden%sprcxree%data))THEN
        work_nj(1:nj) = get_values(pwrden%sprcxree)
     ELSE
        work_nj(1:nj) =zeroc
     ENDIF
     CALL netcdf_err( nf90_put_var(ncid,id_sprcxree,work_nj),id_sprcxree)

 
     IF(ASSOCIATED(pwrden%spreimpe%data))THEN
        work_nj(1:nj) = get_values(pwrden%spreimpe)
     ELSE
        work_nj(1:nj) =zeroc
     ENDIF
     CALL netcdf_err( nf90_put_var(ncid,id_spreimpe,work_nj),id_spreimpe)


    DO j=1,nion
        work_nj_nion(1:nj,j) = get_values(pwrden%brems_nions(j))
     ENDDO
 
     !the species are ordered as 1..nprim,1..nimp
     !where nion = nprim+nimp
     CALL netcdf_err( nf90_put_var(ncid,id_brems_nions,work_nj_nion ),id_brems_nions)
  

! add nclass vpinch  output 5/3/2012 HSJ
     IF( .NOT. ASSOCIATED(diffuse%vpinch_nclass))THEN
          ALLOCATE(diffuse%vpinch_nclass(nion+1))         ! + 1 for electron species
          DO jj =1,nion+1
             diffuse%vpinch_nclass(jj)  = zero_vector(nj)
          ENDDO
     ENDIF


    DO j=1,nion+1
        work_nj_nionp1(1:nj,j) = get_values(diffuse%vpinch_nclass(j))
     ENDDO
 
     !the species are ordered as 1..nprim,1..nimp,el
     !where nion = nprim+nimp
     CALL netcdf_err( nf90_put_var(ncid,id_vpinch_nions,work_nj_nionp1 ),id_vpinch_nions)
  

     work_nj(1:nj) = get_values(pwrden%qohm)
     CALL netcdf_err( nf90_put_var(ncid,id_qohm,work_nj),id_qohm)

     work_npsi(1:mhd_dat%npsi) = get_values(dischg%rmajavnpsi)
     CALL netcdf_err( nf90_put_var(ncid,id_rmajavnpsi,work_npsi),id_rmajavnpsi)

     work_npsi(1:mhd_dat%npsi) = get_values(dischg%rminavnpsi)
     CALL netcdf_err( nf90_put_var(ncid,id_rminavnpsi,work_npsi),id_rminavnpsi)

     work_npsi(1:mhd_dat%npsi) = get_values(dischg%psivolpnpsi)
     CALL netcdf_err( nf90_put_var(ncid,id_psivolpnpsi,work_npsi),id_psivolpnpsi)

     work_npsi(1:mhd_dat%npsi) = get_values(dischg%elongxnpsi)
     CALL netcdf_err( nf90_put_var(ncid,id_elongxnpsi,work_npsi),id_elongxnpsi)

     work_npsi(1:mhd_dat%npsi) = get_values(dischg%triangnpsi_u)
     CALL netcdf_err( nf90_put_var(ncid,id_triangnpsi_u,work_npsi),id_triangnpsi_u)

     work_npsi(1:mhd_dat%npsi) = get_values(dischg%triangnpsi_l)
     CALL netcdf_err( nf90_put_var(ncid,id_triangnpsi_l,work_npsi),id_triangnpsi_l)

     work_npsi(1:mhd_dat%npsi) = get_values(dischg%pindentnpsi)
     CALL netcdf_err( nf90_put_var(ncid,id_pindentnpsi,work_npsi),id_pindentnpsi)

     work_npsi(1:mhd_dat%npsi) = get_values(dischg%sfareanpsi)
     CALL netcdf_err( nf90_put_var(ncid,id_sfareanpsi,work_npsi),id_sfareanpsi)

     work_npsi(1:mhd_dat%npsi) = get_values(dischg%cxareanpsi)
     CALL netcdf_err( nf90_put_var(ncid,id_cxareanpsi,work_npsi),id_cxareanpsi)

     work_npsi(1:mhd_dat%npsi) = get_values(dischg%grho1npsi)
     CALL netcdf_err( nf90_put_var(ncid,id_grho1npsi,work_npsi),id_grho1npsi)
     work_npsi(1:mhd_dat%npsi) = get_values(dischg%grho2npsi)
     CALL netcdf_err( nf90_put_var(ncid,id_grho2npsi,work_npsi),id_grho2npsi)

     IF(ASSOCIATED(mhd_dat%qpsinpsi%data))THEN
        work_npsi(1:mhd_dat%npsi) = get_values(mhd_dat%qpsinpsi)
     ELSE
        work_npsi(1:mhd_dat%npsi) =zeroc
     ENDIF
     CALL netcdf_err( nf90_put_var(ncid,id_qpsinpsi,work_npsi),id_qpsinpsi)

     IF(ASSOCIATED(mhd_dat%pressnpsi%data))THEN
        work_npsi(1:mhd_dat%npsi) = get_values(mhd_dat%pressnpsi)
     ELSE
        work_npsi(1:mhd_dat%npsi) =zeroc
     ENDIF
     CALL netcdf_err( nf90_put_var(ncid,id_pressnpsi,work_npsi),id_pressnpsi)

     IF(ASSOCIATED(mhd_dat%ffprimnpsi%data))THEN
        work_npsi(1:mhd_dat%npsi) = get_values(mhd_dat%ffprimnpsi)
     ELSE
        work_npsi(1:mhd_dat%npsi) =zeroc
     ENDIF
     CALL netcdf_err( nf90_put_var(ncid,id_ffprimnpsi,work_npsi),id_ffprimnpsi)

     IF(ASSOCIATED(mhd_dat%pprimnpsi%data))THEN
        work_npsi(1:mhd_dat%npsi) = get_values(mhd_dat%pprimnpsi)
     ELSE
        work_npsi(1:mhd_dat%npsi) =zeroc
     ENDIF
     CALL netcdf_err( nf90_put_var(ncid,id_pprimnpsi,work_npsi),id_pprimnpsi)



     CALL netcdf_err( nf90_put_var(ncid,id_storqueb,storqueb),id_storqueb)

     CALL netcdf_err( nf90_put_var(ncid,id_totcur_bc,totcur_bc),id_totcur_bc)

     CALL netcdf_err( nf90_put_var(ncid,id_vloop_bc,vloop_bc),id_vloop_bc)

     CALL netcdf_err( nf90_put_var(ncid,id_fix_edge_te_bc,fix_edge_te_bc),id_fix_edge_te_bc)
     CALL netcdf_err( nf90_put_var(ncid,id_fix_edge_ti_bc,fix_edge_ti_bc),id_fix_edge_ti_bc)
     CALL netcdf_err( nf90_put_var(ncid,id_fix_edge_rot_bc,fix_edge_rot_bc),id_fix_edge_rot_bc)
     CALL netcdf_err( nf90_put_var(ncid,id_fix_edge_ni_bc,fix_edge_ni_bc),id_fix_edge_ni_bc)



     CALL netcdf_err( nf90_put_var(ncid,id_te_bc,te_bc),id_te_bc)

     CALL netcdf_err( nf90_put_var(ncid,id_ti_bc,ti_bc),id_ti_bc)




     CALL netcdf_err( nf90_put_var(ncid,id_ene_bc,ene_bc),id_ene_bc)

     CALL netcdf_err( nf90_put_var(ncid,id_zeff_bc,zeff_bc),id_zeff_bc)

     CALL netcdf_err( nf90_put_var(ncid,id_angrot_bc,angrot_bc),id_angrot_bc)

     CALL netcdf_err( nf90_put_var(ncid,id_en_bc,en_bc),id_en_bc)


     CALL netcdf_err( nf90_put_var(ncid,id_flux_bc,flux_bc),id_flux_bc)

     CALL netcdf_err( nf90_put_var(ncid,id_z,z),id_z)
     CALL netcdf_err( nf90_put_var(ncid,id_zsq,zsq),id_zsq)

     CALL netcdf_err( nf90_put_var(ncid,id_wbeam,wbeam),id_wbeam)

     work_nj(1:nj) = get_values(profile%walp)
     WHERE(work_nj < zeroc)work_nj= zeroc     ! interpolation error
     CALL netcdf_err( nf90_put_var(ncid,id_walp,work_nj),id_walp)

     CALL netcdf_err( nf90_put_var(ncid,id_enalp,enalp),id_enalp)

     work_nj(1:nj) = get_values(wpdot%dnedt)
     CALL netcdf_err( nf90_put_var(ncid,id_dnedt,work_nj),id_dnedt) 

     CALL netcdf_err( nf90_put_var(ncid,id_eps,eps),id_eps)
     work_nj(1:nj) = get_values(mhd_dat%rcap)
     CALL netcdf_err( nf90_put_var(ncid,id_rcap,work_nj),id_rcap)
     work_nj(1:nj) = get_values(mhd_dat%rcapi)
     CALL netcdf_err( nf90_put_var(ncid,id_rcapi,work_nj),id_rcapi)
     work_nj(1:nj) = get_values(mhd_dat%r2capi)
     CALL netcdf_err( nf90_put_var(ncid,id_r2capi,work_nj),id_r2capi)
     work_nj(1:nj) = get_values(mhd_dat%r2cap)
     CALL netcdf_err( nf90_put_var(ncid,id_r2cap,work_nj),id_r2cap)

     CALL netcdf_err( nf90_put_var(ncid,id_xi11,xi11),id_xi11)
     CALL netcdf_err( nf90_put_var(ncid,id_xi33,xi33),id_xi33)
     CALL netcdf_err( nf90_put_var(ncid,id_xips,xips),id_xips)
     CALL netcdf_err( nf90_put_var(ncid,id_xhm2,xhm2),id_xhm2)

     CALL netcdf_err( nf90_put_var(ncid,id_dfdt,dfdt),id_dfdt)
     CALL netcdf_err( nf90_put_var(ncid,id_dgdt,dgdt),id_dgdt)       
     CALL netcdf_err( nf90_put_var(ncid,id_dhdt,dhdt),id_dhdt)

     !       CALL netcdf_err( nf90_put_var(ncid,id_nplasbdry,dischg%nplasbdry),id_nplasbdry)
     work_nplasbdry(1:dischg%nplasbdry) = get_values(dischg%rplasbdry)
     CALL netcdf_err( nf90_put_var(ncid,id_rplasbdry,work_nplasbdry),id_rplasbdry)
     work_nplasbdry(1:dischg%nplasbdry) = get_values(dischg%zplasbdry)
     CALL netcdf_err( nf90_put_var(ncid,id_zplasbdry,work_nplasbdry),id_zplasbdry)



     work_nlimiter(1:dischg%nlimiter) = get_values(dischg%rlimiter)
     CALL netcdf_err( nf90_put_var(ncid,id_rlimiter,work_nlimiter),id_rlimiter)
     work_nlimiter(1:dischg%nlimiter) = get_values(dischg%zlimiter)
     CALL netcdf_err( nf90_put_var(ncid,id_zlimiter,work_nlimiter),id_zlimiter)



     IF(neut_beam%nbeams .GT. izero )THEN
        IF(SIZE(bptor) .LT. neut_beam%nbeams)THEN
           ! bptor is dimmed to  kb in nf_param
           WRITE(ncrt,FMT ='( "Error,inconistency in bptor dimension")')
           lerrno = iomaxerr + 358_I4B
           CALL  terminate(lerrno,nlog)
        ENDIF

        CALL netcdf_err( nf90_put_var(ncid,id_pbeam,neut_beam%pbeam),id_pbeam)
        CALL netcdf_err( nf90_put_var(ncid,id_ebeam,neut_beam%ebeam),id_ebeam)
        work_bptor(1:neut_beam%nbeams) = bptor(1:neut_beam%nbeams)
        CALL netcdf_err( nf90_put_var(ncid,id_bptor,work_bptor),id_bptor)


        IF(.NOT. ASSOCIATED(neut_beam%fbcur))THEN
           ke = ke_bm
           ALLOCATE(neut_beam%fbcur(ke,neut_beam%nbeams))
           neut_beam%fbcur(:,:) = zeroc
        ENDIF
        CALL netcdf_err( nf90_put_var(ncid,id_fbcur,neut_beam%fbcur),id_fbcur)

        IF(.NOT. ASSOCIATED(neut_beam%prompt_pwr_in_plasma))THEN
           ALLOCATE(neut_beam%prompt_pwr_in_plasma(ke,neut_beam%nbeams))
           neut_beam%prompt_pwr_in_plasma(:,:) = zeroc
        ENDIF
        CALL netcdf_err( nf90_put_var(ncid,id_prompt_nb_pwr,neut_beam%prompt_pwr_in_plasma),id_prompt_nb_pwr)

        CALL netcdf_err( nf90_put_var(ncid,id_bneut,neut_beam%bneut),id_bneut)
        CALL netcdf_err( nf90_put_var(ncid,id_bion,neut_beam%bion),id_bion)


        CALL netcdf_err( nf90_put_var(ncid,id_fap,neut_beam%fap),id_fap)
        CALL netcdf_err( nf90_put_var(ncid,id_fwall,neut_beam%fwall),id_fwall)
        CALL netcdf_err( nf90_put_var(ncid,id_forb,neut_beam%forb),id_forb)

        CALL netcdf_err( nf90_put_var(ncid,id_fber,neut_beam%fber),id_fber)
        CALL netcdf_err( nf90_put_var(ncid,id_fb00,neut_beam%fb00),id_fb00)
        CALL netcdf_err( nf90_put_var(ncid,id_fb01,neut_beam%fb01),id_fb01)
        CALL netcdf_err( nf90_put_var(ncid,id_fb10,neut_beam%fb10),id_fb10)
        CALL netcdf_err( nf90_put_var(ncid,id_fb11,neut_beam%fb11),id_fb11)
        CALL netcdf_err( nf90_put_var(ncid,id_wb00,neut_beam%wb00),id_wb00)
        CALL netcdf_err( nf90_put_var(ncid,id_wb01,neut_beam%wb01),id_wb01)
        CALL netcdf_err( nf90_put_var(ncid,id_wb10,neut_beam%wb10),id_wb10)
        CALL netcdf_err( nf90_put_var(ncid,id_wb11,neut_beam%wb11),id_wb11)
        CALL netcdf_err( nf90_put_var(ncid,id_sb,neut_beam%sb),id_sb)
        CALL netcdf_err( nf90_put_var(ncid,id_qb,neut_beam%qb),id_qb) 
        CALL netcdf_err( nf90_put_var(ncid,id_spb,neut_beam%spb),id_spb)
        CALL netcdf_err( nf90_put_var(ncid,id_spbr,neut_beam%spbr),id_spbr) 
        CALL netcdf_err( nf90_put_var(ncid,id_pb0,neut_beam%pb0),id_pb0)
        CALL netcdf_err( nf90_put_var(ncid,id_angmpf,neut_beam%angmpf),id_angmpf)
 
        CALL netcdf_err( nf90_put_var(ncid,id_hibr,neut_beam%hibr),id_hibr)
        CALL netcdf_err( nf90_put_var(ncid,id_hdep,neut_beam%hdep),id_hdep)
        CALL netcdf_err( nf90_put_var(ncid,id_zeta,neut_beam%zeta),id_zeta)
        CALL netcdf_err( nf90_put_var(ncid,id_rhog_beam,neut_beam%rhog_beam),id_rhog_beam)
        IF(ALLOCATED(work3d))DEALLOCATE(work3d)
        ALLOCATE(work3d(nj,ke_bm,neut_beam%nbeams))
 
        work3d(:,:,:) = neut_beam%hicm(:,:,:,1)
        CALL netcdf_err( nf90_put_var(ncid,id_hicme,work3d),id_hicme)
        work3d(:,:,:) = neut_beam%hicm(:,:,:,2)
        CALL netcdf_err( nf90_put_var(ncid,id_hicmp1,work3d),id_hicmp1)
        work3d(:,:,:) = neut_beam%hicm(:,:,:,3)
        CALL netcdf_err( nf90_put_var(ncid,id_hicmp2,work3d),id_hicmp2)
        DEALLOCATE(work3d)
     ENDIF



     !------------------------------------------------------------
     ! put out common frequencies:
     ! If these are not defined then define them and set to zero
     ! frequencies that are not flux functions are on the rho grid at zma=0
     ! on the outboard side
     !-----------------------------------------------------------------
 
  IF(.NOT. ALLOCATED(id_plas_freq) .AND. nprimp1 .GT. 1) &
      CALL allocate_plasma_freq
        
     k=0
     DO i=1,nprim
        k = k+1
        IF(namep(i) == 'dt')THEN

           work_nj(:) = get_values(plasma_frequencies%omega_pi(k))
           CALL netcdf_err( nf90_put_var(ncid,id_plas_freq(k),work_nj),id_plas_freq(k))

           work_nj(:) = get_values(plasma_frequencies%omega_ci(k))
           CALL netcdf_err( nf90_put_var(ncid,id_ci_freq(k),work_nj),id_ci_freq(k))

           work_nj(:) = get_values(plasma_frequencies%omega_lh(k))
           CALL netcdf_err( nf90_put_var(ncid,id_lh_freq(k),work_nj),id_lh_freq(k))

           work_nj(:) = get_values(plasma_frequencies%omega_uh(k))
           CALL netcdf_err( nf90_put_var(ncid,id_uh_freq(k),work_nj),id_uh_freq(k))

           k =k+1
           work_nj(:) = get_values(plasma_frequencies%omega_pi(k))
           CALL netcdf_err( nf90_put_var(ncid,id_plas_freq(k),work_nj),id_plas_freq(k))

           work_nj(:) = get_values(plasma_frequencies%omega_ci(k))
           CALL netcdf_err( nf90_put_var(ncid,id_ci_freq(k),work_nj),id_ci_freq(k))

           work_nj(:) = get_values(plasma_frequencies%omega_lh(k))
           CALL netcdf_err( nf90_put_var(ncid,id_lh_freq(k),work_nj),id_lh_freq(k))

           work_nj(:) = get_values(plasma_frequencies%omega_uh(k))
           CALL netcdf_err( nf90_put_var(ncid,id_uh_freq(k),work_nj),id_uh_freq(k))

        ELSE
           work_nj(:) = get_values(plasma_frequencies%omega_pi(k))
           CALL netcdf_err( nf90_put_var(ncid,id_plas_freq(k),work_nj),id_plas_freq(k))

           work_nj(:) = get_values(plasma_frequencies%omega_ci(k))
           CALL netcdf_err( nf90_put_var(ncid,id_ci_freq(k),work_nj),id_ci_freq(k))

           work_nj(:) = get_values(plasma_frequencies%omega_lh(k))
           CALL netcdf_err( nf90_put_var(ncid,id_lh_freq(k),work_nj),id_lh_freq(k))

           work_nj(:) = get_values(plasma_frequencies%omega_uh(k))
           CALL netcdf_err( nf90_put_var(ncid,id_uh_freq(k),work_nj),id_uh_freq(k))

        ENDIF
     ENDDO

     work_nj(:) = get_values(plasma_frequencies%omega_ce)
     CALL netcdf_err( nf90_put_var(ncid,id_ce_freq,work_nj),id_ce_freq)

     work_nj(:) = get_values(plasma_frequencies%omega_pe)
     CALL netcdf_err( nf90_put_var(ncid,id_pe_freq,work_nj),id_pe_freq)







     !-------------------------------------------------------------------------
     !
     !     netcdf read file 
     ! ------------------------------------------------------------------------

  ELSE  rwncd
     !          READ netcdf file: Cant assume that dim ids are as set in this
     !           subroutine because netcdf file most likely is  foreign.


     CALL netcdf_err( nf90_inquire(ncid, nDimensions, nVariables, nAttributes, &
          unlimitedDimId))


     DO j=1,nDimensions  
        CALL netcdf_err( NF90_INQUIRE_DIMENSION(ncid, j, gen_name, idlen)) 
        CALL search_dims(gen_name,dim_names,nf_max_name,dim_size ,max_dims,index )
        label = "error in determination of required dimensions"
        Assign_dim : SELECT CASE(index)
        CASE(0) ! netcdf file does not have this dimension defined.
           WRITE(ncrt,FMT='(a)')label
           WRITE(nlog,FMT='(a)')label
           IF( j == 17)THEN   ! new beam dimension not present in old files
               neut_beam%nj_beam = -1 ! set to nj below 
           ELSE
              lerrno = 33
              CALL  terminate(lerrno,nlog)
           ENDIF
        CASE(1) 
           nj     = idlen
        CASE(2) 
           nion   = idlen
        CASE(3) 
           nprim  = idlen
        CASE(4) 
           IF(name_size .NE. idlen)THEN
              label ='error in character array declaration,name_size ,idlen ='
              WRITE(ncrt,FMT='(a,5x,i5,x,i5)')label,name_size,idlen
              WRITE(nlog,FMT='(a,5x,i5,x,i5)')label,name_size,idlen
              lerrno = 34
              CALL  terminate(lerrno,nlog)
           ENDIF
        CASE(5) 
           nimp  = idlen
        CASE(6) 
           nneu  = idlen
        CASE(7) 
           nbion = idlen
        CASE(8) 
           dischg%nplasbdry = idlen
        CASE(9) 
           ntot = idlen
        CASE(10)
           dischg%nr_mhd = idlen
        CASE(11)
           dischg%nz_mhd = idlen
        CASE(12)
           mhd_dat%npsi =  idlen
        CASE(13)
           dischg%nlimiter = idlen
        CASE(14)
           neut_beam%nbeams = idlen
        CASE(15)
           ke  = ke_bm
           IF(ke .ne. ke_bm)THEN
              WRITE(ncrt,FMT='(a)')label
              WRITE(nlog,FMT='(a)')label
              lerrno = 33
              CALL  terminate(lerrno,nlog)
           ENDIF
        CASE(16)
           nionp1 = idlen
        CASE(17)
           neut_beam%nj_beam = idlen   ! may  or may not equal nj
        CASE DEFAULT
           WRITE(ncrt,FMT='(a)')label
           WRITE(nlog,FMT='(a)')label
           lerrno = 33
           CALL  terminate(lerrno,nlog)
        END SELECT Assign_dim
     ENDDO




     IF(ALLOCATED(work_nj))DEALLOCATE(work_nj)
     IF(ALLOCATED(work_npsi))DEALLOCATE(work_npsi)
     IF(ALLOCATED(work_nplasbdry))DEALLOCATE(work_nplasbdry)
     IF(ALLOCATED(work_nlimiter))DEALLOCATE(work_nlimiter)
     IF(ALLOCATED(work_nr))DEALLOCATE(work_nr)
     IF(ALLOCATED(work_nz))DEALLOCATE(work_nz)


     IF(ASSOCIATED(mhd_dat%psi))DEALLOCATE(mhd_dat%psi)

     ALLOCATE(work_nj(nj))
     ALLOCATE(work_npsi(mhd_dat%npsi))
     ALLOCATE(work_nplasbdry(dischg%nplasbdry))

     ALLOCATE(work_nlimiter(dischg%nlimiter))
     ALLOCATE(work_nr(dischg%nr_mhd))
     ALLOCATE(work_nz(dischg%nz_mhd))


     ALLOCATE(mhd_dat%psi(dischg%nr_mhd,dischg%nz_mhd))

     !     in the future we may introduce time dependent metrics.
     !     then dfdt (= d/dt FCAP ,etc) will be calculated internally in this code
     !     by reading in fcap_bc, etc.
     !     in demo code just assume time independent for each
     !     slice that the solver is called.
     IF(.NOT. ALLOCATED(dfdt))THEN 
        ALLOCATE(dfdt(nj))
        dfdt(:) = zeroc
     ENDIF
     IF(.NOT. ALLOCATED(dgdt))THEN 
        ALLOCATE(dgdt(nj))
        dgdt(:) = zeroc
     ENDIF
     IF(.NOT. ALLOCATED(dhdt))THEN
        ALLOCATE(dhdt(nj))
        dhdt(:) = zeroc
     ENDIF


     label(:) = " "
     CALL netcdf_err(  &
          nf90_inquire_attribute(ncid,NF90_GLOBAL,"title", len = titlen))
     IF( titlen .LE. strlen) &
          CALL netcdf_err( nf90_get_att(ncid,NF90_GLOBAL,"title",title))

     CALL netcdf_err( nf90_inq_varid(ncid,"shot", id_shot),id_shot) !given name get id
     CALL netcdf_err( nf90_get_var(ncid,id_shot,shot_id%shot_nmbr),id_shot) !given id, get value(s) 
     CALL netcdf_err( nf90_inq_varid(ncid,"nion", id_nion),id_nion) 
     CALL netcdf_err( nf90_get_var(ncid,id_nion,nion),id_nion) 
     CALL netcdf_err( nf90_inq_varid(ncid,"nprim",id_nprim),id_nprim) 




     CALL netcdf_err( nf90_get_var(ncid,id_nprim,nprim),id_nprim)
     nprimp1 = nprim+1 ! code always starts by reading a statefile so this is OK

 
     CALL netcdf_err( nf90_inq_varid(ncid,"fd_thermal",id_fd_thermal),id_fd_thermal) 
     CALL netcdf_err( nf90_get_var(ncid,id_fd_thermal,fd_thermal),id_fd_thermal)
     CALL netcdf_err( nf90_inq_varid(ncid,"nimp",id_nimp),id_nimp)
     CALL netcdf_err( nf90_get_var(ncid,id_nimp,nimp),id_nimp)
     CALL netcdf_err( nf90_inq_varid(ncid,"nneu",id_nneu),id_nneu)
     CALL netcdf_err( nf90_get_var(ncid,id_nneu,nneu),id_nneu)
     CALL netcdf_err( nf90_inq_varid(ncid,"nbion",id_nbion),id_nbion)
     CALL netcdf_err( nf90_get_var(ncid,id_nbion,nbion),id_nbion)
     CALL netcdf_err( nf90_inq_varid(ncid,"fd_beam",id_fd_beam),id_fd_beam)
     CALL netcdf_err( nf90_get_var(ncid,id_fd_beam,fd_beam),id_fd_beam)
     id_nbeams_er = izero
     CALL netcdf_err( nf90_inq_varid(ncid,"nbeams",id_nbeams),id_nbeams_er)
     IF(id_nbeams_er == -1)THEN
        neut_beam%nbeams =0
     ELSE
        id_nbeams_er = izero
        CALL netcdf_err( nf90_get_var(ncid,id_nbeams,neut_beam%nbeams),id_nbeams_er)
     ENDIF
    

     IF( .NOT. ASSOCIATED(namep))&
        ALLOCATE ( namep(nprim),namei(nimp),namen(nneu),nameb(nbion))
     IF( .NOT. ASSOCIATED(profile%en))&
        ALLOCATE (profile%en(nion))   !2d array 
     IF( .NOT. ASSOCIATED(profile%flux))&
        ALLOCATE (profile%flux(ntot)) !2d array
     IF( .NOT. ASSOCIATED(profile%zsq))&
        ALLOCATE (profile%zsq(ntot))  !2d array 
     IF( .NOT. ASSOCIATED(profile%z))&
        ALLOCATE (profile%z(ntot))    !2d array 
    


     ! nj,nion are now known  so allocate some arrays: 
     IF( ALLOCATED(work_nj_nion))DEALLOCATE(work_nj_nion)
     ALLOCATE(work_nj_nion(nj,nion))
     IF( ALLOCATED(work_nj_nionp1))DEALLOCATE(work_nj_nionp1)
     ALLOCATE(work_nj_nionp1(nj,nion+1))
     IF(.NOT. ALLOCATED(te_bc))ALLOCATE(te_bc(nj))
     IF(.NOT. ALLOCATED(ti_bc))ALLOCATE(ti_bc(nj))
     IF(.NOT. ALLOCATED(ene_bc))ALLOCATE(ene_bc(nj))
     IF(.NOT. ALLOCATED(zeff_bc))ALLOCATE(zeff_bc(nj))
     IF(.NOT. ALLOCATED(angrot_bc))ALLOCATE(angrot_bc(nj))
     IF(.NOT. ALLOCATED(eps))ALLOCATE(eps(nj))
     IF(.NOT. ALLOCATED(storqueb))ALLOCATE(storqueb(nj))
     IF(.NOT. ALLOCATED(wbeam))ALLOCATE(wbeam(nj))
     IF(.NOT. ALLOCATED(walp))ALLOCATE(walp(nj))
     IF(.NOT. ALLOCATED(w_alpha))ALLOCATE(w_alpha(nj))
     IF(.NOT. ALLOCATED(enalp))ALLOCATE(enalp(nj))
!!$     IF(.NOT. ALLOCATED(dnedt))ALLOCATE(dnedt(nj))
!!$     IF(.NOT. ALLOCATED(dangrotdt))ALLOCATE(dangrotdt(nj))
     IF(.NOT. ASSOCIATED(zeff))ALLOCATE(zeff(nj))
     IF(.NOT. ASSOCIATED(z))ALLOCATE(z(nj,nion))
     IF(.NOT. ALLOCATED(en_bc))ALLOCATE(en_bc(nj,nion))
     IF(.NOT. ALLOCATED(flux_bc))ALLOCATE(flux_bc(nj,ntot))
!!$     IF(.NOT. ALLOCATED(dnidt))ALLOCATE(dnidt(nj,nion))
     IF(.NOT. ASSOCIATED(zsq))ALLOCATE(zsq(nj,nion))


     CALL netcdf_err( nf90_inq_varid(ncid,"namep",id_namep),id_namep)
     CALL netcdf_err( nf90_get_var(ncid,id_namep,namep),id_namep)
     CALL netcdf_err( nf90_inq_varid(ncid,"namei",id_namei))
     CALL netcdf_err( nf90_get_var(ncid,id_namei,namei),id_namei)
     CALL netcdf_err( nf90_inq_varid(ncid,"namen",id_namen))
     CALL netcdf_err( nf90_get_var(ncid,id_namen,namen),id_namen)
     CALL netcdf_err( nf90_inq_varid(ncid,"nameb",id_nameb))
     CALL netcdf_err( nf90_get_var(ncid,id_nameb,nameb),id_nameb)
     CALL netcdf_err( nf90_inq_varid(ncid,"namepel",id_namepel))
     CALL netcdf_err( nf90_get_var(ncid,id_namepel,pellet%name),id_namepel)

     CALL netcdf_err( nf90_inq_varid(ncid,"time", id_time),id_time) 
     CALL netcdf_err( nf90_get_var(ncid,id_time,time),id_time) 

     eqtime = time       !equilibirum time. It is assumed that the quatities
     !fcap,gcap,hcap,eps,xhm2,xi11,xi33,xips
     !were calculated at this time.

     CALL netcdf_err( nf90_inq_varid(ncid,"psiaxis", id_psiaxis),id_psiaxis) 
     CALL netcdf_err( nf90_get_var(ncid,id_psiaxis,mhd_dat%psiaxis),id_psiaxis) 

     CALL netcdf_err( nf90_inq_varid(ncid,"psibdry", id_psibdry),id_psibdry) 
     CALL netcdf_err( nf90_get_var(ncid,id_psibdry,mhd_dat%psibdry),id_psibdry)

     CALL netcdf_err( nf90_inq_varid(ncid,"rgeom", id_rgeom),id_rgeom) 
     CALL netcdf_err( nf90_get_var(ncid,id_rgeom,dischg%rgeom),id_rgeom) 

     CALL netcdf_err( nf90_inq_varid(ncid,"btgeom", id_btgeom),id_btgeom) 
     CALL netcdf_err( nf90_get_var(ncid,id_btgeom,dischg%btgeom),id_btgeom) 

     !           CALL netcdf_err( nf90_inq_varid(ncid,"rmag", id_rmag),id_rmag))
     !           CALL netcdf_err( nf90_get_var(ncid,id_rmag,dischg%rmag),id_rmag) 

     CALL netcdf_err( nf90_inq_varid(ncid,"rma", id_rma),id_rma) 
     CALL netcdf_err( nf90_get_var(ncid,id_rma,dischg%rma),id_rma)     
     dischg%rmag = dischg%rma
     CALL netcdf_err( nf90_inq_varid(ncid,"zma", id_zma),id_zma) 
     CALL netcdf_err( nf90_get_var(ncid,id_zma,dischg%zma),id_zma) 


     CALL netcdf_err( nf90_inq_varid(ncid,"rsep", id_rsep),id_rsep) 
     CALL netcdf_err( nf90_get_var(ncid,id_rsep,dischg%rsep),id_rsep) 

     CALL netcdf_err( nf90_inq_varid(ncid,"zsep", id_zsep),id_zsep) 
     CALL netcdf_err( nf90_get_var(ncid,id_zsep,dischg%zsep),id_zsep) 

     CALL netcdf_err( nf90_inq_varid(ncid,"rmajor", id_rmajor),id_rmajor) 
     CALL netcdf_err( nf90_get_var(ncid,id_rmajor,dischg%rmajor),id_rmajor) 
     mhd_dat%R0 = dischg%rmajor

     CALL netcdf_err( nf90_inq_varid(ncid,"rplasmin", id_rplasmin),id_rplasmin) 
     CALL netcdf_err( nf90_get_var(ncid,id_rplasmin,dischg%rplasmin),id_rplasmin) 
     CALL netcdf_err( nf90_inq_varid(ncid,"rplasmax", id_rplasmax),id_rplasmax) 
     CALL netcdf_err( nf90_get_var(ncid,id_rplasmax,dischg%rplasmax),id_rplasmax)

     CALL netcdf_err( nf90_inq_varid(ncid,"zplasmin", id_zplasmin),id_zplasmin) 
     CALL netcdf_err( nf90_get_var(ncid,id_zplasmin,dischg%zplasmin),id_zplasmin) 
     CALL netcdf_err( nf90_inq_varid(ncid,"zplasmax", id_zplasmax),id_zplasmax) 
     CALL netcdf_err( nf90_get_var(ncid,id_zplasmax,dischg%zplasmax),id_zplasmax)

     CALL netcdf_err( nf90_inq_varid(ncid,"kappa", id_kappa),id_kappa) 
     CALL netcdf_err( nf90_get_var(ncid,id_kappa,dischg%kappa),id_kappa) 
     CALL netcdf_err( nf90_inq_varid(ncid,"deltao", id_deltao),id_deltao) 
     CALL netcdf_err( nf90_get_var(ncid,id_deltao,dischg%deltao),id_deltao) 
     CALL netcdf_err( nf90_inq_varid(ncid,"pindento", id_pindento),id_pindento) 
     CALL netcdf_err( nf90_get_var(ncid,id_pindento,dischg%pindento),id_pindento) 
     CALL netcdf_err( nf90_inq_varid(ncid,"volume", id_volume),id_volume) 
     CALL netcdf_err( nf90_get_var(ncid,id_volume,dischg%volo),id_volume) 
     CALL netcdf_err( nf90_inq_varid(ncid,"circum", id_circum),id_circum) 
     CALL netcdf_err( nf90_get_var(ncid,id_circum,dischg%circum),id_circum) 
     CALL netcdf_err( nf90_inq_varid(ncid,"areao", id_areao),id_areao) 
     CALL netcdf_err( nf90_get_var(ncid,id_areao,dischg%areao),id_areao)
     CALL netcdf_err( nf90_inq_varid(ncid,"btor", id_btor),id_btor) 
     CALL netcdf_err( nf90_get_var(ncid,id_btor,mhd_dat%btor),id_btor)
     CALL netcdf_err( nf90_inq_varid(ncid,"tot_cur", id_tot_cur),id_tot_cur) 
     CALL netcdf_err( nf90_get_var(ncid,id_tot_cur,mhd_dat%tot_cur),id_tot_cur)
     CALL netcdf_err( nf90_inq_varid(ncid,"totohm_cur", id_totohm_cur),id_totohm_cur) 
     CALL netcdf_err( nf90_get_var(ncid,id_totohm_cur,mhd_dat%totohm_cur),id_totohm_cur)
     CALL netcdf_err( nf90_inq_varid(ncid,"totboot_cur", id_totboot_cur),id_totboot_cur) 
     CALL netcdf_err( nf90_get_var(ncid,id_totboot_cur,mhd_dat%totboot_cur),id_totboot_cur)
     CALL netcdf_err( nf90_inq_varid(ncid,"totbeam_cur", id_totbeam_cur),id_totbeam_cur) 
     CALL netcdf_err( nf90_get_var(ncid,id_totbeam_cur,mhd_dat%totbeam_cur),id_totbeam_cur)
     CALL netcdf_err( nf90_inq_varid(ncid,"totrf_cur", id_totrf_cur),id_totrf_cur) 
     CALL netcdf_err( nf90_get_var(ncid,id_totrf_cur,mhd_dat%totrf_cur),id_totrf_cur)
     ibcur =1_I4B ;irfc = 1_I4B
     IF(ABS(mhd_dat%totbeam_cur) .LT. 1.e-5)ibcur = 0
     IF(ABS(mhd_dat%totrf_cur)   .LT. 1.e-5)irfc  = 0

     CALL netcdf_err( nf90_inq_varid(ncid,"betap", id_betap),id_betap) 
     CALL netcdf_err( nf90_get_var(ncid,id_betap,mhd_dat%betap),id_betap)
     CALL netcdf_err( nf90_inq_varid(ncid,"beta", id_beta),id_beta) 
     CALL netcdf_err( nf90_get_var(ncid,id_beta,mhd_dat%beta),id_beta)
     CALL netcdf_err( nf90_inq_varid(ncid,"ali", id_ali),id_ali) 
     CALL netcdf_err( nf90_get_var(ncid,id_ali,mhd_dat%ali),id_ali)
     CALL netcdf_err( nf90_inq_varid(ncid,"te0", id_te0),id_te0) 
     CALL netcdf_err( nf90_get_var(ncid,id_te0,profile%te0),id_te0)
     CALL netcdf_err( nf90_inq_varid(ncid,"ti0", id_ti0),id_ti0) 
     CALL netcdf_err( nf90_get_var(ncid,id_ti0,profile%ti0),id_ti0)


     id_pfuse_tot_flag  = izero
     CALL netcdf_err( nf90_inq_varid(ncid,"pfuse_tot", id_pfuse_tot),id_pfuse_tot_flag)
     IF(id_pfuse_tot_flag == -1)THEN
        pfuse_tot = zeroc
     ELSE
        CALL netcdf_err( nf90_get_var(ncid,id_pfuse_tot,pfuse_tot),id_pfuse_tot)
     ENDIF


     id_qrad_tot_flag  = izero
     id_brems_tot_flag = izero
     CALL netcdf_err( nf90_inq_varid(ncid,"qrad_tot", id_qrad_tot),id_qrad_tot_flag) 
     IF(id_qrad_tot_flag == -1)THEN
        qrad_tot = zeroc
     ELSE
        CALL netcdf_err( nf90_get_var(ncid,id_qrad_tot,qrad_tot),id_qrad_tot)
     ENDIF

     CALL netcdf_err( nf90_inq_varid(ncid,"brems_tot", id_brems_tot),id_brems_tot_flag)
     IF(.NOT. ALLOCATED(brems_tot))ALLOCATE(brems_tot(nion))
     IF(id_brems_tot_flag == -1)THEN 
          brems_tot(:) = zeroc
     ELSE
        CALL netcdf_err( nf90_get_var(ncid,id_brems_tot,brems_tot),id_brems_tot)
     ENDIF


     

     CALL netcdf_err( nf90_inq_varid(ncid,"psir_grid", id_psir_grid),id_psir_grid) 
     CALL netcdf_err( nf90_get_var(ncid,id_psir_grid,work_nj),id_psir_grid)
     psir_grid = new_Vector(nj,work_nj) ! [volt*sec/rad]

     CALL check_monotonic( work_nj, nj, monotonic,one)


     IF (.NOT. monotonic) THEN    ! psir_grid is not monotonic
        WRITE  (ncrt, 17) time
17       FORMAT (' the psi grid,psir_grid, read in  from iterdb statefile' / &
             ' is not monotonic'        /                                    &
             ' therefore  we do not continue with this run',/,               &
             ' the time value is ', 1pe14.8)
        WRITE(nlog,17)time
        lerrno = 13_I4B
        CALL  terminate(lerrno,nlog)
     END IF

     id_qpsinpsi_flag = izero
     CALL netcdf_err( nf90_inq_varid(ncid,"qpsinpsi", id_qpsinpsi),id_qpsinpsi_flag) 
     IF(id_qpsinpsi_flag .NE. -1)THEN 
        CALL netcdf_err( nf90_get_var(ncid,id_qpsinpsi,work_npsi),id_qpsinpsi)
        mhd_dat%mhd_info_avail = .TRUE.
     ELSE
        mhd_dat%mhd_info_avail = .FALSE.
        work_npsi(:) = zeroc
     ENDIF
     mhd_dat%qpsinpsi  = new_Vector(mhd_dat%npsi,work_npsi)

     id_pressnpsi_flag = izero
     CALL netcdf_err( nf90_inq_varid(ncid,"pressnpsi", id_pressnpsi),id_pressnpsi_flag) 
     IF(id_pressnpsi_flag .NE. -1)THEN 
        CALL netcdf_err( nf90_get_var(ncid,id_pressnpsi,work_npsi),id_pressnpsi)
     ELSE
        work_npsi(:) = zeroc
     ENDIF
     mhd_dat%pressnpsi  = new_Vector(mhd_dat%npsi,work_npsi)
 

     id_ffprimnpsi_flag = izero
     CALL netcdf_err( nf90_inq_varid(ncid,"ffprimnpsi", id_ffprimnpsi),id_ffprimnpsi_flag) 
     IF(id_ffprimnpsi_flag .NE. -1)THEN 
        CALL netcdf_err( nf90_get_var(ncid,id_ffprimnpsi,work_npsi),id_ffprimnpsi)
     ELSE
        work_npsi(:) = zeroc
     ENDIF
     mhd_dat%ffprimnpsi  = new_Vector(mhd_dat%npsi,work_npsi)

     id_pprimnpsi_flag = izero
     CALL netcdf_err( nf90_inq_varid(ncid,"pprimnpsi", id_pprimnpsi),id_pprimnpsi_flag) 
     IF(id_pprimnpsi_flag .NE. -1)THEN 
        CALL netcdf_err( nf90_get_var(ncid,id_pprimnpsi,work_npsi),id_pprimnpsi)
     ELSE
        work_npsi(:) = zeroc
     ENDIF
     mhd_dat%pprimnpsi  = new_Vector(mhd_dat%npsi,work_npsi)


 
     CALL netcdf_err( nf90_inq_varid(ncid,"psi", id_psi),id_psi) 
     CALL netcdf_err( nf90_get_var(ncid,id_psi,mhd_dat%psi),id_psi) ! [volt*sec/rad]

     CALL netcdf_err( nf90_inq_varid(ncid,"rho_mhd_gridnpsi", id_rho_mhd_gridnpsi),id_rho_mhd_gridnpsi) 
     CALL netcdf_err( nf90_get_var(ncid,id_rho_mhd_gridnpsi,work_npsi),id_rho_mhd_gridnpsi)
     rho_mhd_gridnpsi = new_Vector(mhd_dat%npsi,work_npsi)


     CALL netcdf_err( nf90_inq_varid(ncid,"rmhdgrid", id_rmhdgrid),id_rmhdgrid) 
     CALL netcdf_err( nf90_get_var(ncid,id_rmhdgrid,work_nr),id_rmhdgrid)
     dischg%rmhdgrid = new_Vector(dischg%nr_mhd,work_nr)

     CALL netcdf_err( nf90_inq_varid(ncid,"zmhdgrid", id_zmhdgrid),id_zmhdgrid) 
     CALL netcdf_err( nf90_get_var(ncid,id_zmhdgrid,work_nz),id_zmhdgrid)
     dischg%zmhdgrid = new_Vector(dischg%nz_mhd,work_nz)

     CALL netcdf_err( nf90_inq_varid(ncid,"rho_grid", id_rho_grid),id_rho_grid) 
     CALL netcdf_err( nf90_get_var(ncid,id_rho_grid,work_nj),id_rho_grid)
     rho_grid = new_Vector(nj,work_nj)
     elmt = get_element(rho_grid,nj)
     elmt = 1.0_DP/elmt
     rho_gridn = real_mult_Vector (elmt,rho_grid)


     CALL netcdf_err( nf90_inq_varid(ncid,"fcap", id_fcap),id_fcap) 
     CALL netcdf_err( nf90_get_var(ncid,id_fcap,work_nj), id_fcap)
     mhd_dat%fcap = new_Vector(nj,work_nj)

     CALL netcdf_err( nf90_inq_varid(ncid,"gcap", id_gcap),id_gcap) 
     CALL netcdf_err( nf90_get_var(ncid,id_gcap,work_nj),id_gcap)
     mhd_dat%gcap = new_Vector(nj,work_nj)

     CALL netcdf_err( nf90_inq_varid(ncid,"hcap", id_hcap),id_hcap) 
     CALL netcdf_err( nf90_get_var(ncid,id_hcap,work_nj),id_hcap)
     mhd_dat%hcap = new_Vector(nj,work_nj)

     CALL netcdf_err( nf90_inq_varid(ncid,"q_value", id_q_value),id_q_value) 
     CALL netcdf_err( nf90_get_var(ncid,id_q_value,work_nj),id_q_value)
     mhd_dat%q_value = new_Vector(nj,work_nj)

     id_betan_er = izero
     CALL netcdf_err( nf90_inq_varid(ncid,"betan", id_betan),id_betan_er) 
     IF(id_betan_er .NE. -1)THEN
       CALL netcdf_err( nf90_get_var(ncid,id_betan,work_nj),id_betan)
       mhd_dat%betan = new_Vector(nj,work_nj)
     ELSE
       mhd_dat%betan = zero_Vector(nj)
     ENDIF



     CALL netcdf_err( nf90_inq_varid(ncid,"Te", id_te),id_te) 
     CALL netcdf_err( nf90_get_var(ncid,id_te,work_nj),id_te)
     profile%te = new_Vector(nj,work_nj)

     CALL netcdf_err( nf90_inq_varid(ncid,"Ti", id_ti),id_ti) 
     CALL netcdf_err( nf90_get_var(ncid,id_ti,work_nj),id_ti)
     profile%ti = new_Vector(nj,work_nj)

     CALL netcdf_err( nf90_inq_varid(ncid,"press", id_press),id_press) 
     CALL netcdf_err( nf90_get_var(ncid,id_press,work_nj),id_press)
     profile%press = new_Vector(nj,work_nj)

     CALL netcdf_err( nf90_inq_varid(ncid,"pressb", id_pressb),id_pressb) 
     CALL netcdf_err( nf90_get_var(ncid,id_pressb,work_nj),id_pressb)
     profile%pressb = new_Vector(nj,work_nj)





     CALL netcdf_err( nf90_inq_varid(ncid,"ene", id_ene),id_ene) 
     CALL netcdf_err( nf90_get_var(ncid,id_ene,work_nj),id_ene)
     profile%ene = new_Vector(nj,work_nj)

     id_flux_elct_er = izero
     CALL netcdf_err( nf90_inq_varid(ncid,"p_flux_elct", id_flux_elct),id_flux_elct_er) 
     IF(id_flux_elct_er .NE. -1)THEN
        CALL netcdf_err( nf90_get_var(ncid,id_flux_elct,work_nj),id_flux_elct)
        profile%fluxe = new_Vector(nj,work_nj)
     ELSE
        profile%fluxe = zero_Vector(nj)
     ENDIF

     id_flux_ion_er = izero
     CALL netcdf_err( nf90_inq_varid(ncid,"p_flux_ion", id_flux_ion),id_flux_ion_er) 
     IF(id_flux_ion_er .NE. -1)THEN
        CALL netcdf_err( nf90_get_var(ncid,id_flux_ion,work_nj),id_flux_ion)
        profile%fluxi = new_Vector(nj,work_nj)
     ELSE
        profile%fluxi = zero_Vector(nj)
     ENDIF

    
     CALL netcdf_err( nf90_inq_varid(ncid,"enion", id_enion),id_enion) 
     CALL netcdf_err( nf90_get_var(ncid,id_enion,work_nj_nion),id_enion)
     !the species are ordered as 1..npriim,1..nimp
     !where nion = nprim+nimp
     DO jj=1,nion
        profile%en(jj)  = new_Vector(nj,work_nj_nion(1,jj))
     ENDDO




     flag = 0
     CALL netcdf_err( nf90_inq_varid(ncid,"xnus", id_xnus),flag) 
     IF(.NOT. ASSOCIATED(diffuse%xnus))THEN 
         ALLOCATE(diffuse%xnus(nion))
     ENDIF
     IF(flag == -1)THEN ! file does not contain xnus 
         DO jj =1,nion
             diffuse%xnus(jj) = zero_Vector(nj)
         ENDDO
     ELSE
         CALL netcdf_err( nf90_get_var(ncid,id_xnus,work_nj_nion),id_xnus)
         !the species are ordered as 1..nprim,1..nimp
         !where nion = nprim+nimp
         DO jj=1,nion
            diffuse%xnus(jj)  = new_Vector(nj,work_nj_nion(1,jj))
         ENDDO
     ENDIF
 
 
     flag = 0
     CALL netcdf_err( nf90_inq_varid(ncid,"xnuse", id_xnuse),flag) 
     IF(flag == -1)THEN
        diffuse%xnuse = zero_vector(nj) 
     ELSE
        CALL netcdf_err( nf90_get_var(ncid,id_xnuse,work_nj),id_xnuse)
        diffuse%xnuse = new_Vector(nj,work_nj)
     ENDIF

     flag = 0
     CALL netcdf_err( nf90_inq_varid(ncid,"ftrap", id_ftrap),flag) 
     IF(flag == -1)THEN
        diffuse%ftrap = zero_vector(nj) 
     ELSE
        CALL netcdf_err( nf90_get_var(ncid,id_ftrap,work_nj),id_ftrap)
        diffuse%ftrap = new_Vector(nj,work_nj)
     ENDIF


     flag =0
     CALL netcdf_err( nf90_inq_varid(ncid,"eta", id_eta),flag) 
     IF(flag == -1)THEN
        diffuse%eta = zero_vector(nj) 
     ELSE
        CALL netcdf_err( nf90_get_var(ncid,id_eta,work_nj),id_eta)
        diffuse%eta = new_Vector(nj,work_nj)
     ENDIF

     flag =0
     CALL netcdf_err( nf90_inq_varid(ncid,"chiepc", id_chiepc),flag) 
     IF(flag == -1)THEN
        diffuse%chie_paleo = zero_Vector(nj)
     ELSE
        CALL netcdf_err( nf90_get_var(ncid,id_chiepc,work_nj),id_chiepc)
        diffuse%chie_paleo = new_Vector(nj,work_nj)
        !chie_paleo(:) = chie_paleo(:)*1.e4_DP !Convert to cm^2/sec
     ENDIF


     IF(.NOT. ASSOCIATED(profile%flux)) &
                        ALLOCATE (profile%flux(ntot)) ! 2d array 
     CALL netcdf_err( nf90_inq_varid(ncid,"p_flux", id_pflux),id_pflux) 
     CALL netcdf_err( nf90_get_var(ncid,id_pflux,work_nj_nion),id_pflux)
     DO jj=1,nion
        profile%flux(jj)  = new_Vector(nj,work_nj_nion(1,jj))
     ENDDO


     id_pflux_conv_er = izero
     IF(.NOT. ASSOCIATED(profile%flux_conv)) &
                        ALLOCATE (profile%flux_conv(ntot)) ! 2d array 
     CALL netcdf_err( nf90_inq_varid(ncid,"p_flux_conv", id_pflux_conv),id_pflux_conv_er) 
     IF(id_pflux_conv_er .NE. -1)THEN
        CALL netcdf_err( nf90_get_var(ncid,id_pflux_conv,work_nj_nion),id_pflux_conv)
        DO jj=1,nion
           profile%flux_conv(jj)  = new_Vector(nj,work_nj_nion(1,jj))
        ENDDO
     ELSE
        DO jj=1,nion
           profile%flux_conv(jj)  = zero_Vector(nj)
        ENDDO
     ENDIF


     CALL netcdf_err( nf90_inq_varid(ncid,"e_fluxe", id_e_fluxe),id_e_fluxe) 
     CALL netcdf_err( nf90_get_var(ncid,id_e_fluxe,work_nj),id_e_fluxe)
     profile%flux(nion+1)  = new_Vector(nj,work_nj)


     id_e_fluxe_conv_er = izero
     CALL netcdf_err( nf90_inq_varid(ncid,"e_fluxe_conv", id_e_fluxe_conv),id_e_fluxe_conv_er) 
     IF(id_e_fluxe_conv_er .NE. -1)THEN
        CALL netcdf_err( nf90_get_var(ncid,id_e_fluxe_conv,work_nj),id_e_fluxe_conv)
        profile%flux_conv(nion+1)  = new_Vector(nj,work_nj)

     ELSE
        profile%flux_conv(nion+1)  = zero_Vector(nj)

     ENDIF


     !---------------------------------------------------------------------------
     ! multimode read section
     !---------------------------------------------------------------------------
    
     oknf = delete_Vector_nf(diffuse%mmm_gammaDBM) ! _nf ==> non fatal if vector is not ASSOCIATED
     id_mmm_gammaDBM_er = izero
     CALL netcdf_err( nf90_inq_varid(ncid,"mmm_gammaDBM", id_mmm_gammaDBM),id_mmm_gammaDBM_er) 
     IF(id_mmm_gammaDBM_er .ne. -1)THEN
        CALL netcdf_err( nf90_get_var(ncid,id_mmm_gammaDBM,work_nj),id_mmm_gammaDBM)
        diffuse%mmm_gammaDBM  = new_Vector(nj,work_nj)
     ELSE
        diffuse%mmm_gammaDBM = zero_vector(nj)
     ENDIF

     oknf = delete_Vector_nf(diffuse%mmm_omegaDBM) ! _nf ==> non fatal if vector is not ASSOCIATED
     id_mmm_omegaDBM_er = izero
     CALL netcdf_err( nf90_inq_varid(ncid,"mmm_omegaDBM", id_mmm_omegaDBM),id_mmm_omegaDBM_er) 
     IF(id_mmm_omegaDBM_er .NE. -1)THEN
        CALL netcdf_err( nf90_get_var(ncid,id_mmm_omegaDBM,work_nj),id_mmm_omegaDBM)
        diffuse%mmm_omegaDBM  = new_Vector(nj,work_nj)
     ELSE
        diffuse%mmm_omegaDBM = zero_vector(nj)
     ENDIF


     oknf = delete_Vector_nf(diffuse%mmm_xdi) ! _nf ==> non fatal if vector is not ASSOCIATED
     id_mmm_xdi_er = izero
     CALL netcdf_err( nf90_inq_varid(ncid,"mmm_xdi", id_mmm_xdi),id_mmm_xdi_er) 
     IF(id_mmm_xdi_er .NE. -1)THEN
        CALL netcdf_err( nf90_get_var(ncid,id_mmm_xdi,work_nj),id_mmm_xdi)
        diffuse%mmm_xdi  = new_Vector(nj,work_nj)
     ELSE
        diffuse%mmm_xdi  = zero_Vector(nj)
     ENDIF

     oknf = delete_Vector_nf(diffuse%mmm_xti) ! _nf ==> non fatal if vector is not ASSOCIATED
     id_mmm_xti_er = izero
     CALL netcdf_err( nf90_inq_varid(ncid,"mmm_xti", id_mmm_xti),id_mmm_xti_er)
     IF(id_mmm_xti_er .NE. -1)THEN
        CALL netcdf_err( nf90_get_var(ncid,id_mmm_xti,work_nj),id_mmm_xti)
        diffuse%mmm_xti  = new_Vector(nj,work_nj)
     ELSE
        diffuse%mmm_xti = zero_vector(nj)
     ENDIF

     oknf = delete_Vector_nf(diffuse%mmm_xte) ! _nf ==> non fatal if vector is not ASSOCIATED
     id_mmm_xte_er  = izero
     CALL netcdf_err( nf90_inq_varid(ncid,"mmm_xte", id_mmm_xte),id_mmm_xte_er) 
     IF(id_mmm_xte_er .NE. -1)THEN
        CALL netcdf_err( nf90_get_var(ncid,id_mmm_xte,work_nj),id_mmm_xte)
        diffuse%mmm_xte  = new_Vector(nj,work_nj)
     ELSE
        diffuse%mmm_xte = zero_vector(nj)
     ENDIF

     oknf = delete_Vector_nf(diffuse%mmm_xdz) ! _nf ==> non fatal if vector is not ASSOCIATED
     id_mmm_xdz_er  = izero
     CALL netcdf_err( nf90_inq_varid(ncid,"mmm_xdz", id_mmm_xdz),id_mmm_xdz_er) 
     IF(id_mmm_xdz_er .NE. -1)THEN
        CALL netcdf_err( nf90_get_var(ncid,id_mmm_xdz,work_nj),id_mmm_xdz)
        diffuse%mmm_xdz  = new_Vector(nj,work_nj)
     ELSE
        diffuse%mmm_xdz = zero_vector(nj)
     ENDIF

     oknf = delete_Vector_nf(diffuse%mmm_xvt) ! _nf ==> non fatal if vector is not ASSOCIATED
     id_mmm_xvt_er  = izero
     CALL netcdf_err( nf90_inq_varid(ncid,"mmm_xvt", id_mmm_xvt),id_mmm_xvt_er) 
     IF(id_mmm_xvt_er .NE. -1)THEN
        CALL netcdf_err( nf90_get_var(ncid,id_mmm_xvt,work_nj),id_mmm_xvt)
        diffuse%mmm_xvt  = new_Vector(nj,work_nj)
     ELSE
        diffuse%mmm_xvt = zero_vector(nj)
     ENDIF

     oknf = delete_Vector_nf(diffuse%mmm_xvp) ! _nf ==> non fatal if vector is not ASSOCIATED
     id_mmm_xvp_er  = izero
     CALL netcdf_err( nf90_inq_varid(ncid,"mmm_xvp", id_mmm_xvp),id_mmm_xvp_er) 
     IF(id_mmm_xvp_er .NE. -1)THEN
        CALL netcdf_err( nf90_get_var(ncid,id_mmm_xvp,work_nj),id_mmm_xvp)
        diffuse%mmm_xvp  = new_Vector(nj,work_nj)
     ELSE
        diffuse%mmm_xvp = zero_vector(nj)
     ENDIF


     oknf = delete_Vector_nf(diffuse%mmm_xtiW20) ! _nf ==> non fatal if vector is not ASSOCIATED
     id_mmm_xtiW20_er = izero
     CALL netcdf_err( nf90_inq_varid(ncid,"mmm_xtiW20", id_mmm_xtiW20),id_mmm_xtiW20_er) 
     IF(id_mmm_xtiW20_er .NE. -1)THEN
        CALL netcdf_err( nf90_get_var(ncid,id_mmm_xtiW20,work_nj),id_mmm_xtiW20)
        diffuse%mmm_xtiW20  = new_Vector(nj,work_nj)
     ELSE
        diffuse%mmm_xtiW20 = zero_vector(nj)
     ENDIF

     oknf = delete_Vector_nf(diffuse%mmm_xdiW20) ! _nf ==> non fatal if vector is not ASSOCIATED
     id_mmm_xdiW20_er  = izero
     CALL netcdf_err( nf90_inq_varid(ncid,"mmm_xdiW20", id_mmm_xdiW20),id_mmm_xdiW20_er) 
     IF(id_mmm_xdiW20_er .NE. -1)THEN
        CALL netcdf_err( nf90_get_var(ncid,id_mmm_xdiW20,work_nj),id_mmm_xdiW20)
        diffuse%mmm_xdiW20  = new_Vector(nj,work_nj)
     ELSE
        diffuse%mmm_xdiW20 = zero_vector(nj)
     ENDIF

     oknf = delete_Vector_nf(diffuse%mmm_xteW20) ! _nf ==> non fatal if vector is not ASSOCIATED
     id_mmm_xteW20_er  = izero
     CALL netcdf_err( nf90_inq_varid(ncid,"mmm_xteW20", id_mmm_xteW20),id_mmm_xteW20_er) 
     IF(id_mmm_xteW20_er .NE. -1)THEN
        CALL netcdf_err( nf90_get_var(ncid,id_mmm_xteW20,work_nj),id_mmm_xteW20)
        diffuse%mmm_xteW20  = new_Vector(nj,work_nj)
     ELSE
         diffuse%mmm_xteW20 = zero_vector(nj)
     ENDIF

     oknf = delete_Vector_nf(diffuse%mmm_xtiDBM) ! _nf ==> non fatal if vector is not ASSOCIATED
     id_mmm_xtiDBM_er  = izero
     CALL netcdf_err( nf90_inq_varid(ncid,"mmm_xtiDBM", id_mmm_xtiDBM),id_mmm_xtiDBM_er)
     IF(id_mmm_xtiDBM_er .NE. -1)THEN
        CALL netcdf_err( nf90_get_var(ncid,id_mmm_xtiDBM,work_nj),id_mmm_xtiDBM)
        diffuse%mmm_xtiDBM  = new_Vector(nj,work_nj)
     ELSE
        diffuse%mmm_xtiDBM = zero_vector(nj)
     ENDIF

     oknf = delete_Vector_nf(diffuse%mmm_xdiDBM) ! _nf ==> non fatal if vector is not ASSOCIATED
     id_mmm_xdiDBM_er  = izero
     CALL netcdf_err( nf90_inq_varid(ncid,"mmm_xdiDBM", id_mmm_xdiDBM),id_mmm_xdiDBM_er)
     IF(id_mmm_xdiDBM_er .NE. -1)THEN
        CALL netcdf_err( nf90_get_var(ncid,id_mmm_xdiDBM,work_nj),id_mmm_xdiDBM)
        diffuse%mmm_xdiDBM  = new_Vector(nj,work_nj)
     ELSE
        diffuse%mmm_xdiDBM = zero_vector(nj)
     ENDIF

     oknf = delete_Vector_nf(diffuse%mmm_xteDBM) ! _nf ==> non fatal if vector is not ASSOCIATED
     id_mmm_xteDBM_er = izero
     CALL netcdf_err( nf90_inq_varid(ncid,"mmm_xteDBM", id_mmm_xteDBM),id_mmm_xteDBM_er) 
     IF(id_mmm_xteDBM_er .NE. -1)THEN
        CALL netcdf_err( nf90_get_var(ncid,id_mmm_xteDBM,work_nj),id_mmm_xteDBM)
        diffuse%mmm_xteDBM  = new_Vector(nj,work_nj)
     ELSE
        diffuse%mmm_xteDBM = zero_vector(nj)
     ENDIF

     oknf = delete_Vector_nf(diffuse%mmm_xteETG) ! _nf ==> non fatal if vector is not ASSOCIATED
     id_mmm_xteETG_er  = izero
     CALL netcdf_err( nf90_inq_varid(ncid,"mmm_xteETG", id_mmm_xteETG),id_mmm_xteETG_er) 
     IF(id_mmm_xteETG_er .NE. -1)THEN
        CALL netcdf_err( nf90_get_var(ncid,id_mmm_xteETG,work_nj),id_mmm_xteETG)
        diffuse%mmm_xteETG  = new_Vector(nj,work_nj)
     ELSE
        diffuse%mmm_xteETG = zero_vector(nj)
     ENDIF


      ngW20 = 4
      IF(ASSOCIATED(diffuse%mmm_gammaW20))DEALLOCATE(diffuse%mmm_gammaW20)
      ALLOCATE(diffuse%mmm_gammaW20(ngW20))         
        DO jj =1,ngW20
         IF(jj == 1) THEN
           id_mmm_gamma_i1_W20_er = izero
           CALL netcdf_err( nf90_inq_varid(ncid,"mmm_gamma_i1_W20", id_mmm_gamma_i1_W20),id_mmm_gamma_i1_W20_er) 
           IF(id_mmm_gamma_i1_W20_er .NE. -1)THEN
              CALL netcdf_err( nf90_get_var(ncid,id_mmm_gamma_i1_W20,work_nj),id_mmm_gamma_i1_W20)
              diffuse%mmm_gammaW20(jj)  = new_Vector(nj,work_nj)
           ELSE
              diffuse%mmm_gammaW20(jj)  = zero_vector(nj)
           ENDIF
         ELSEIF(jj == 2)THEN
           id_mmm_gamma_e1_W20_er = izero
           CALL netcdf_err( nf90_inq_varid(ncid,"mmm_gamma_e1_W20", id_mmm_gamma_e1_W20),id_mmm_gamma_e1_W20_er) 
           IF(id_mmm_gamma_e1_W20_er .NE. -1)THEN
              CALL netcdf_err( nf90_get_var(ncid,id_mmm_gamma_e1_W20,work_nj),id_mmm_gamma_e1_W20)
              diffuse%mmm_gammaW20(jj)  = new_Vector(nj,work_nj)
           ELSE
              diffuse%mmm_gammaW20(jj)  = zero_vector(nj)
           ENDIF
         ELSEIF(jj == 3)THEN
           id_mmm_gamma_i2_W20_er = izero
           CALL netcdf_err( nf90_inq_varid(ncid,"mmm_gamma_i2_W20", id_mmm_gamma_i2_W20),id_mmm_gamma_i2_W20_er) 
           IF(id_mmm_gamma_i2_W20_er .NE. -1)THEN
              CALL netcdf_err( nf90_get_var(ncid,id_mmm_gamma_i2_W20,work_nj),id_mmm_gamma_i2_W20)
              diffuse%mmm_gammaW20(jj)  = new_Vector(nj,work_nj)
           ELSE
              diffuse%mmm_gammaW20(jj)  = zero_vector(nj)
           ENDIF
         ELSEIF(jj == 4)THEN
           id_mmm_gamma_e2_W20_er = izero
           CALL netcdf_err( nf90_inq_varid(ncid,"mmm_gamma_e2_W20", id_mmm_gamma_e2_W20),id_mmm_gamma_e2_W20_er) 
           IF(id_mmm_gamma_e2_W20_er .NE. -1)THEN
              CALL netcdf_err( nf90_get_var(ncid,id_mmm_gamma_e2_W20,work_nj),id_mmm_gamma_e2_W20)
              diffuse%mmm_gammaW20(jj)  = new_Vector(nj,work_nj)
           ELSE
              diffuse%mmm_gammaW20(jj)  = zero_vector(nj)
           ENDIF
         ENDIF
        ENDDO



      ngW20 = 4
      IF(ASSOCIATED(diffuse%mmm_omegaW20))DEALLOCATE(diffuse%mmm_omegaW20)
      ALLOCATE(diffuse%mmm_omegaW20(ngW20))         
        DO jj =1,ngW20
         IF(jj == 1) THEN
           id_mmm_omega_i1_W20_er = izero
           CALL netcdf_err( nf90_inq_varid(ncid,"mmm_omega_i1_W20", id_mmm_omega_i1_W20),id_mmm_omega_i1_W20_er) 
           IF(id_mmm_omega_i1_W20_er .NE. -1)THEN
              CALL netcdf_err( nf90_get_var(ncid,id_mmm_omega_i1_W20,work_nj),id_mmm_omega_i1_W20)
              diffuse%mmm_omegaW20(jj)  = new_Vector(nj,work_nj)
           ELSE
              diffuse%mmm_omegaW20(jj)  = zero_vector(nj)
           ENDIF
         ELSEIF(jj == 2)THEN
           id_mmm_omega_e1_W20_er = izero
           CALL netcdf_err( nf90_inq_varid(ncid,"mmm_omega_e1_W20", id_mmm_omega_e1_W20),id_mmm_omega_e1_W20_er) 
           IF(id_mmm_omega_e1_W20_er .NE. -1)THEN
              CALL netcdf_err( nf90_get_var(ncid,id_mmm_omega_e1_W20,work_nj),id_mmm_omega_e1_W20)
              diffuse%mmm_omegaW20(jj)  = new_Vector(nj,work_nj)
           ELSE
              diffuse%mmm_omegaW20(jj)  = zero_vector(nj)
           ENDIF
         ELSEIF(jj == 3)THEN
           id_mmm_omega_i2_W20_er = izero
           CALL netcdf_err( nf90_inq_varid(ncid,"mmm_omega_i2_W20", id_mmm_omega_i2_W20),id_mmm_omega_i2_W20_er) 
           IF(id_mmm_omega_i2_W20_er .NE. -1)THEN
              CALL netcdf_err( nf90_get_var(ncid,id_mmm_omega_i2_W20,work_nj),id_mmm_omega_i2_W20)
              diffuse%mmm_omegaW20(jj)  = new_Vector(nj,work_nj)
           ELSE
              diffuse%mmm_omegaW20(jj)  = zero_vector(nj)
           ENDIF
         ELSEIF(jj == 4)THEN
           id_mmm_omega_e2_W20_er = izero
           CALL netcdf_err( nf90_inq_varid(ncid,"mmm_omega_e2_W20", id_mmm_omega_e2_W20),id_mmm_omega_e2_W20_er) 
           IF(id_mmm_omega_e2_W20_er .NE. -1)THEN
              CALL netcdf_err( nf90_get_var(ncid,id_mmm_omega_e2_W20,work_nj),id_mmm_omega_e2_W20)
              diffuse%mmm_omegaW20(jj)  = new_Vector(nj,work_nj)
           ELSE
              diffuse%mmm_omegaW20(jj)  = zero_vector(nj)
           ENDIF
         ENDIF
        ENDDO


      

      ngW20 = 4
      IF(ASSOCIATED(diffuse%mmm_vflux))DEALLOCATE(diffuse%mmm_vflux)
      ALLOCATE(diffuse%mmm_vflux(ngW20))         
        DO jj =1,ngW20
         IF(jj == 1) THEN
           id_mmm_flux_ith_er = izero
           CALL netcdf_err( nf90_inq_varid(ncid,"mmm_flux_ith",id_mmm_flux_ith),id_mmm_flux_ith_er) 
           IF(id_mmm_flux_ith_er .NE. -1)THEN
              CALL netcdf_err( nf90_get_var(ncid,id_mmm_flux_ith,work_nj),id_mmm_flux_ith)
              diffuse%mmm_vflux(jj)  = new_Vector(nj,work_nj)
           ELSE
              diffuse%mmm_vflux(jj)  = zero_vector(nj)
           ENDIF
         ELSEIF(jj == 2)THEN
           id_mmm_flux_ip_er = izero
           CALL netcdf_err( nf90_inq_varid(ncid,"mmm_flux_ip", id_mmm_flux_ip),id_mmm_flux_ip_er) 
           IF(id_mmm_flux_ip_er .NE. -1)THEN
              CALL netcdf_err( nf90_get_var(ncid,id_mmm_flux_ip,work_nj),id_mmm_flux_ip)
              diffuse%mmm_vflux(jj)  = new_Vector(nj,work_nj)
           ELSE
              diffuse%mmm_vflux(jj)  = zero_vector(nj)
           ENDIF
         ELSEIF(jj == 3)THEN
           id_mmm_flux_eth_er = izero
           CALL netcdf_err( nf90_inq_varid(ncid,"mmm_flux_eth", id_mmm_flux_eth),id_mmm_flux_eth_er) 
           IF(id_mmm_flux_eth_er .NE. -1)THEN
              CALL netcdf_err( nf90_get_var(ncid,id_mmm_flux_eth,work_nj),id_mmm_flux_eth)
              diffuse%mmm_vflux(jj)  = new_Vector(nj,work_nj)
           ELSE
              diffuse%mmm_vflux(jj)  = zero_vector(nj)
           ENDIF
         ELSEIF(jj == 4)THEN
           id_mmm_flux_imp_er = izero
           CALL netcdf_err( nf90_inq_varid(ncid,"mmm_flux_imp", id_mmm_flux_imp),id_mmm_flux_imp_er) 
           IF(id_mmm_flux_imp_er .NE. -1)THEN
              CALL netcdf_err( nf90_get_var(ncid,id_mmm_flux_imp,work_nj),id_mmm_flux_imp)
              diffuse%mmm_vflux(jj)  = new_Vector(nj,work_nj)
           ELSE
              diffuse%mmm_vflux(jj)  = zero_vector(nj)
           ENDIF

         ENDIF
        ENDDO

 

      ngW20 = 6
      IF(ASSOCIATED(diffuse%mmm_vconv))DEALLOCATE(diffuse%mmm_vconv)
      ALLOCATE(diffuse%mmm_vconv(ngW20))

      DO jj =1,ngW20
         IF(jj == 1) THEN
           id_mmm_vconv_ith_er = izero
           CALL netcdf_err( nf90_inq_varid(ncid,"mmm_vconv_ith",id_mmm_vconv_ith),id_mmm_vconv_ith_er) 
           IF(id_mmm_vconv_ith_er .NE. -1)THEN
              CALL netcdf_err( nf90_get_var(ncid,id_mmm_vconv_ith,work_nj),id_mmm_vconv_ith)
              diffuse%mmm_vconv(jj)  = new_Vector(nj,work_nj)
           ELSE
              diffuse%mmm_vconv(jj)  = zero_vector(nj)
           ENDIF


         ELSEIF(jj == 2)THEN
           id_mmm_vconv_ip_er = izero
           CALL netcdf_err( nf90_inq_varid(ncid,"mmm_vconv_ip", id_mmm_vconv_ip),id_mmm_vconv_ip_er) 
           IF(id_mmm_vconv_ip_er .NE. -1)THEN
              CALL netcdf_err( nf90_get_var(ncid,id_mmm_vconv_ip,work_nj),id_mmm_vconv_ip)
              diffuse%mmm_vconv(jj)  = new_Vector(nj,work_nj)
           ELSE
              diffuse%mmm_vconv(jj)  = zero_vector(nj)
           ENDIF

         ELSEIF(jj == 3)THEN
           id_mmm_vconv_eth_er = izero
           CALL netcdf_err( nf90_inq_varid(ncid,"mmm_vconv_eth", id_mmm_vconv_eth),id_mmm_vconv_eth_er) 
           IF(id_mmm_vconv_eth_er .NE. -1)THEN
              CALL netcdf_err( nf90_get_var(ncid,id_mmm_vconv_eth,work_nj),id_mmm_vconv_eth)
              diffuse%mmm_vconv(jj)  = new_Vector(nj,work_nj)
           ELSE
              diffuse%mmm_vconv(jj)  = zero_vector(nj)
           ENDIF

         ELSEIF(jj == 4)THEN
           id_mmm_vconv_imp_er = izero
           CALL netcdf_err( nf90_inq_varid(ncid,"mmm_vconv_imp", id_mmm_vconv_imp),id_mmm_vconv_imp_er) 
           IF(id_mmm_vconv_imp_er .NE. -1)THEN
              CALL netcdf_err( nf90_get_var(ncid,id_mmm_vconv_imp,work_nj),id_mmm_vconv_imp)
              diffuse%mmm_vconv(jj)  = new_Vector(nj,work_nj)
           ELSE
              diffuse%mmm_vconv(jj)  = zero_vector(nj)
           ENDIF

         ELSEIF(jj == 5)THEN
           id_mmm_vmtmt_er = izero
           CALL netcdf_err( nf90_inq_varid(ncid,"mmm_vmtmt", id_mmm_vmtmt),id_mmm_vmtmt_er) 
           IF(id_mmm_vmtmt_er .NE. -1)THEN
              CALL netcdf_err( nf90_get_var(ncid,id_mmm_vmtmt,work_nj),id_mmm_vmtmt)
              diffuse%mmm_vconv(jj)  = new_Vector(nj,work_nj)
           ELSE
              diffuse%mmm_vconv(jj)  = zero_vector(nj)
           ENDIF
         ELSEIF(jj == 6)THEN
           id_mmm_vmtmp_er = izero
           CALL netcdf_err( nf90_inq_varid(ncid,"mmm_vmtmp", id_mmm_vmtmp),id_mmm_vmtmp_er) 
           IF(id_mmm_vmtmp_er .NE. -1)THEN
              CALL netcdf_err( nf90_get_var(ncid,id_mmm_vmtmp,work_nj),id_mmm_vmtmp)
              diffuse%mmm_vconv(jj)  = new_Vector(nj,work_nj)
           ELSE
              diffuse%mmm_vconv(jj)  = zero_vector(nj)
           ENDIF

         ENDIF



        ENDDO



     ! end multimode read section


     IF(ALLOCATED(glf_e_output))DEALLOCATE(glf_e_output)
     ALLOCATE(glf_e_output(nj,3))
     glf_e_output(:,:) = zeroc
     id_glf_elct_eng_flux_er = izero
     CALL netcdf_err( nf90_inq_varid(ncid,"glf_elct_eng_flux", id_glf_elct_eng_flux),id_glf_elct_eng_flux_er) 
     IF(id_glf_elct_eng_flux_er .NE. -1)THEN
         CALL netcdf_err( nf90_get_var(ncid,id_glf_elct_eng_flux,work_nj),id_glf_elct_eng_flux)
         glf_e_output(:,1) = work_nj(:)
     ENDIF

     id_glf_primion_eng_flux_er = izero
     CALL netcdf_err( nf90_inq_varid(ncid,"glf_primion_eng_flux",id_glf_primion_eng_flux),id_glf_primion_eng_flux_er)
     IF(id_glf_primion_eng_flux_er .NE. -1)THEN
        CALL netcdf_err( nf90_get_var(ncid,id_glf_primion_eng_flux,work_nj),id_glf_primion_eng_flux)
        glf_e_output(:,2) = work_nj(:)
     ENDIF
     id_glf_impion_eng_flux_er = izero
     CALL netcdf_err( nf90_inq_varid(ncid,"glf_impion_eng_flux", id_glf_impion_eng_flux),id_glf_impion_eng_flux_er) 
     IF(id_glf_impion_eng_flux_er .NE. -1)THEN
        CALL netcdf_err( nf90_get_var(ncid,id_glf_impion_eng_flux,work_nj),id_glf_impion_eng_flux)
        glf_e_output(:,3) = work_nj(:)
     ENDIF

    IF(ALLOCATED(glf_p_output))DEALLOCATE(glf_p_output)
    ALLOCATE(glf_p_output(nj,3))
    glf_p_output(:,:) = zeroc
    id_glf_elct_partcl_flux_er = izero
    CALL netcdf_err( nf90_inq_varid(ncid,"glf_elct_partcl_flux", id_glf_elct_partcl_flux),id_glf_elct_partcl_flux_er) 
    IF(id_glf_elct_partcl_flux_er .NE. -1)THEN
        CALL netcdf_err( nf90_get_var(ncid,id_glf_elct_partcl_flux,work_nj),id_glf_elct_partcl_flux)
        glf_p_output(:,1) = work_nj(:)
    ENDIF

    id_glf_primion_partcl_flux_er = izero
    CALL netcdf_err( nf90_inq_varid(ncid,"glf_primion_partcl_flux", id_glf_primion_partcl_flux),id_glf_primion_partcl_flux_er) 
    IF(id_glf_primion_partcl_flux_er .NE. -1)THEN
        CALL netcdf_err( nf90_get_var(ncid,id_glf_primion_partcl_flux,work_nj),id_glf_primion_partcl_flux)
        glf_p_output(:,2) = work_nj(:)
    ENDIF
    id_glf_impion_partcl_flux_er = izero
    CALL netcdf_err( nf90_inq_varid(ncid,"glf_impion_partcl_flux", id_glf_impion_partcl_flux),id_glf_impion_partcl_flux_er)
    IF(id_glf_impion_partcl_flux_er .NE. -1)THEN 
        CALL netcdf_err( nf90_get_var(ncid,id_glf_impion_partcl_flux,work_nj),id_glf_impion_partcl_flux)
        glf_p_output(:,3) = work_nj(:)
    ENDIF

 

    IF(ALLOCATED(glf_m_output))DEALLOCATE(glf_m_output)
    ALLOCATE(glf_m_output(nj,3))
    glf_m_output(:,:) = zeroc
    id_glf_elct_momtm_flux_er = izero
    CALL netcdf_err( nf90_inq_varid(ncid,"glf_elct_momtm_flux", id_glf_elct_momtm_flux),id_glf_elct_momtm_flux_er) 
    IF(id_glf_elct_momtm_flux_er .NE. -1)THEN
       CALL netcdf_err( nf90_get_var(ncid,id_glf_elct_momtm_flux,work_nj),id_glf_elct_momtm_flux)
       glf_m_output(:,1) = work_nj(:)
    ENDIF

    id_glf_primion_momtm_flux_er = izero
    CALL netcdf_err( nf90_inq_varid(ncid,"glf_primion_momtm_flux", id_glf_primion_momtm_flux),id_glf_primion_momtm_flux) 
    IF(id_glf_primion_momtm_flux_er .NE. -1)THEN
        CALL netcdf_err( nf90_get_var(ncid,id_glf_primion_momtm_flux,work_nj),id_glf_primion_momtm_flux)
        glf_m_output(:,2) = work_nj(:)
    ENDIF

    id_glf_impion_momtm_flux_er = izero
    CALL netcdf_err( nf90_inq_varid(ncid,"glf_impion_momtm_flux", id_glf_impion_momtm_flux),id_glf_impion_momtm_flux_er) 
    IF(id_glf_impion_momtm_flux_er .NE. -1)THEN
        CALL netcdf_err( nf90_get_var(ncid,id_glf_impion_momtm_flux,work_nj),id_glf_impion_momtm_flux)
        glf_m_output(:,3) = work_nj(:)
    ENDIF

    ! read glf eff diffusivites 
    IF(ALLOCATED(glf_eff_chi_output))DEALLOCATE(glf_eff_chi_output)
    ALLOCATE(glf_eff_chi_output(ntot,nj))
    IF(ALLOCATED(id_glf_eff_chi_err))DEALLOCATE(id_glf_eff_chi_err)
    ALLOCATE(id_glf_eff_chi_err(ntot))
    IF(ALLOCATED(id_glf_eff_chi))DEALLOCATE(id_glf_eff_chi)
    ALLOCATE(id_glf_eff_chi(ntot))
    DO i =1,nprim
       id_glf_eff_chi_err(i) = izero 
        base_label(1:LEN(base_label))=''
        base_label = "glf_eff_particle_chi_"//namep(i)
       CALL netcdf_err( nf90_inq_varid(ncid,base_label, id_glf_eff_chi(i)),id_glf_eff_chi_err(i)) 
       IF(id_glf_eff_chi_err(i) .NE. -1)THEN
           CALL netcdf_err( nf90_get_var(ncid,id_glf_eff_chi(i),work_nj),id_glf_eff_chi(i))
           glf_eff_chi_output(:,i) = work_nj(:)
       ELSE
          glf_eff_chi_output(:,i) = zeroc
       ENDIF
     ENDDO
     DO i = nprim+1,nion
       id_glf_eff_chi_err(i) = izero 
        base_label(1:LEN(base_label))=''
        base_label = "glf_eff_particle_chi_"//namei(i-nprim)
       CALL netcdf_err( nf90_inq_varid(ncid,base_label, id_glf_eff_chi(i)),id_glf_eff_chi_err(i)) 
       IF(id_glf_eff_chi_err(i) .NE. -1)THEN
           CALL netcdf_err( nf90_get_var(ncid,id_glf_eff_chi(i),work_nj),id_glf_eff_chi(i))
           glf_eff_chi_output(:,i) = work_nj(:)
       ELSE
          glf_eff_chi_output(:,i) = zeroc
       ENDIF
     ENDDO
     i =nion+1
     id_glf_eff_chi_err(i) = izero
     CALL netcdf_err( nf90_inq_varid(ncid,"glf_eff_elect_eng_chi", id_glf_eff_chi(i)),id_glf_eff_chi_err(i)) 
     IF(id_glf_eff_chi_err(i) .NE. -1)THEN
           CALL netcdf_err( nf90_get_var(ncid,id_glf_eff_chi(i),work_nj),id_glf_eff_chi(i))
           glf_eff_chi_output(:,i) = work_nj(:)
     ELSE
          glf_eff_chi_output(:,i) = zeroc
     ENDIF

     i =nion+2
     id_glf_eff_chi_err(i) = izero
     CALL netcdf_err( nf90_inq_varid(ncid,"glf_eff_ion_eng_chi", id_glf_eff_chi(i)),id_glf_eff_chi_err(i)) 
     IF(id_glf_eff_chi_err(i) .NE. -1)THEN
           CALL netcdf_err( nf90_get_var(ncid,id_glf_eff_chi(i),work_nj),id_glf_eff_chi(i))
           glf_eff_chi_output(:,i) = work_nj(:)
     ELSE
          glf_eff_chi_output(:,i) = zeroc
     ENDIF

     i =nion+3  ! glf Bp not set , reads zeros
     id_glf_eff_chi_err(i) = izero
     CALL netcdf_err( nf90_inq_varid(ncid,"glf_eff_ion_eng_chi", id_glf_eff_chi(i)),id_glf_eff_chi_err(i)) 
     IF(id_glf_eff_chi_err(i) .NE. -1)THEN
           CALL netcdf_err( nf90_get_var(ncid,id_glf_eff_chi(i),work_nj),id_glf_eff_chi(i))
           glf_eff_chi_output(:,i) = work_nj(:)
     ELSE
          glf_eff_chi_output(:,i) = zeroc
     ENDIF



     i =nion+4
     id_glf_eff_chi_err(i) = izero
     CALL netcdf_err( nf90_inq_varid(ncid,"glf_eff_momtm_chi", id_glf_eff_chi(i)),id_glf_eff_chi_err(i)) 
     IF(id_glf_eff_chi_err(i) .NE. -1)THEN
           CALL netcdf_err( nf90_get_var(ncid,id_glf_eff_chi(i),work_nj),id_glf_eff_chi(i))
           glf_eff_chi_output(:,i) = work_nj(:)
     ELSE
          glf_eff_chi_output(:,i) = zeroc
     ENDIF


     IF(ALLOCATED(glf_etg_output))DEALLOCATE(glf_etg_output)
     ALLOCATE(glf_etg_output(nj))
     id_glf_etg_flux_er = izero ! ==> take no action in netcdf_err
     CALL netcdf_err( nf90_inq_varid(ncid,"glf_etg_eng_flux", id_glf_etg_flux),id_glf_etg_flux_er)
     IF(id_glf_etg_flux_er .NE. -1)THEN 
        CALL netcdf_err( nf90_get_var(ncid,id_glf_etg_flux,glf_etg_output),id_glf_etg_flux)
     ELSE
        glf_etg_output(:) = zeroc ! netcdf file does not contain variable glf_etg_eng_flux
     ENDIF

     IF(ALLOCATED(glf_gamma_net_i_output))DEALLOCATE(glf_gamma_net_i_output)
     ALLOCATE(glf_gamma_net_i_output(nj))
     id_glf_gam_net_i_er = zeroc
     CALL netcdf_err( nf90_inq_varid(ncid,"glf_gamma_net_i", id_glf_gam_net_i),id_glf_gam_net_i_er) 
     IF(id_glf_gam_net_i_er .NE. -1)THEN
         CALL netcdf_err( nf90_get_var(ncid,id_glf_gam_net_i,glf_gamma_net_i_output),id_glf_gam_net_i)
     ELSE
         glf_gamma_net_i_output(:) = zeroc
     ENDIF

     IF(ALLOCATED(glf_gamma_net_e_output))DEALLOCATE(glf_gamma_net_e_output)
     ALLOCATE(glf_gamma_net_e_output(nj))
     id_glf_gam_net_e_er = zeroc
     CALL netcdf_err( nf90_inq_varid(ncid,"glf_gamma_net_e", id_glf_gam_net_e),id_glf_gam_net_e_er)
     IF(id_glf_gam_net_e_er .NE. -1)THEN
         CALL netcdf_err( nf90_get_var(ncid,id_glf_gam_net_e,glf_gamma_net_e_output),id_glf_gam_net_e)
     ELSE
         glf_gamma_net_e_output(:) = zeroc
     ENDIF

     IF(ALLOCATED(glf_anfreq_output))DEALLOCATE(glf_anfreq_output)
     ALLOCATE(glf_anfreq_output(nj))
     id_glf_anfreq_er = zeroc
     CALL netcdf_err( nf90_inq_varid(ncid,"glf_anfreq", id_glf_anfreq),id_glf_anfreq_er)
     IF(id_glf_anfreq_er .NE. -1)Then
        CALL netcdf_err( nf90_get_var(ncid,id_glf_anfreq,glf_anfreq_output),id_glf_anfreq)
     ELSE
        glf_anfreq_output(:) = zeroc
     ENDIF


     IF(ALLOCATED(glf_anfreq2_output))DEALLOCATE(glf_anfreq2_output)
     ALLOCATE(glf_anfreq2_output(nj))
     id_glf_anfreq2_er = zeroc
     CALL netcdf_err( nf90_inq_varid(ncid,"glf_anfreq2", id_glf_anfreq2),id_glf_anfreq2_er)

     IF(id_glf_anfreq2_er .NE. -1)Then

        CALL netcdf_err( nf90_get_var(ncid,id_glf_anfreq2,glf_anfreq2_output),id_glf_anfreq2)

     ELSE
        glf_anfreq2_output(:) = zeroc
     ENDIF


     IF(ALLOCATED(glf_anrate_output))DEALLOCATE(glf_anrate_output)
     ALLOCATE(glf_anrate_output(nj))
     id_glf_anrate_er = izero
     CALL netcdf_err( nf90_inq_varid(ncid,"glf_anrate", id_glf_anrate),id_glf_anrate_er)
     IF(id_glf_anrate_er .NE. -1)THEN
        CALL netcdf_err( nf90_get_var(ncid,id_glf_anrate,glf_anrate_output),id_glf_anrate)
     ELSE
        glf_anrate_output(:) = zeroc
     ENDIF

     IF(ALLOCATED(glf_anrate2_output))DEALLOCATE(glf_anrate2_output)
     ALLOCATE(glf_anrate2_output(nj))
     id_glf_anrate2_er = izero
     CALL netcdf_err( nf90_inq_varid(ncid,"glf_anrate2", id_glf_anrate2),id_glf_anrate2_er)
     IF(id_glf_anrate2_er .NE. -1)THEN
        CALL netcdf_err( nf90_get_var(ncid,id_glf_anrate2,glf_anrate2_output),id_glf_anrate2)
     ELSE
        glf_anrate2_output(:) = zeroc
     ENDIF

     CALL netcdf_err( nf90_inq_varid(ncid,"e_fluxi", id_e_fluxi),id_e_fluxi) 
     CALL netcdf_err( nf90_get_var(ncid,id_e_fluxi,work_nj),id_e_fluxi)
     profile%flux(nion+2)  = new_Vector(nj,work_nj)

     id_e_fluxi_conv_er = izero 
     CALL netcdf_err( nf90_inq_varid(ncid,"e_fluxi_conv", id_e_fluxi_conv),id_e_fluxi_conv_er) 
     IF(id_e_fluxi_conv_er .NE. -1)THEN
        CALL netcdf_err( nf90_get_var(ncid,id_e_fluxi_conv,work_nj),id_e_fluxi_conv)
        profile%flux_conv(nion+2)  = new_Vector(nj,work_nj)

     ELSE
        profile%flux_conv(nion+2)  = zero_Vector(nj)
     ENDIF

     CALL netcdf_err( nf90_inq_varid(ncid,"fday_flux", id_fdyflux),id_fdyflux) 
     CALL netcdf_err( nf90_get_var(ncid,id_fdyflux,work_nj),id_fdyflux)
     profile%flux(nion+3)  = new_Vector(nj,work_nj)

     id_fdyflux_conv_er = izero
     CALL netcdf_err( nf90_inq_varid(ncid,"fday_flux_conv", id_fdyflux_conv),id_fdyflux_conv_er) 
     IF(id_fdyflux_conv_er .NE. -1)THEN
        CALL netcdf_err( nf90_get_var(ncid,id_fdyflux_conv,work_nj),id_fdyflux_conv)
        profile%flux_conv(nion+3)  = new_Vector(nj,work_nj)
     ELSE
        profile%flux_conv(nion+3)  = zero_Vector(nj)
     ENDIF


     CALL netcdf_err( nf90_inq_varid(ncid,"rot_flux", id_rotflux),id_rotflux) 
     CALL netcdf_err( nf90_get_var(ncid,id_rotflux,work_nj),id_rotflux)
     profile%flux(nion+4)  = new_Vector(nj,work_nj)

     id_rotflux_conv_er = izero
     CALL netcdf_err( nf90_inq_varid(ncid,"rot_flux_conv", id_rotflux_conv),id_rotflux_conv_er) 
     IF(id_rotflux_conv_er .NE. -1)THEN
        CALL netcdf_err( nf90_get_var(ncid,id_rotflux_conv,work_nj),id_rotflux_conv)
        profile%flux_conv(nion+4)  = new_Vector(nj,work_nj)
     ELSE
        profile%flux_conv(nion+4)  = zero_Vector(nj)
     ENDIF

 
     IF(ALLOCATED(tglf_p_output))DEALLOCATE(tglf_p_output)
     ALLOCATE(tglf_p_output(nj,3))
     CALL netcdf_err( nf90_inq_varid(ncid,"tglf_elct_p_flux", id_tglf_p_fluxe),id_tglf_p_fluxe) 
     CALL netcdf_err( nf90_get_var(ncid,id_tglf_p_fluxe,work_nj),id_tglf_p_fluxe)
     tglf_p_output(:,1) = work_nj(:)  ! NOTE : to load tglf_p,e_flux* need to convert to zc grid
     CALL netcdf_err( nf90_inq_varid(ncid,"tglf_ion_p_flux",id_tglf_p_fluxp), id_tglf_p_fluxp) 
     CALL netcdf_err( nf90_get_var(ncid,id_tglf_p_fluxp,work_nj),id_tglf_p_fluxp)
     tglf_p_output(:,2) = work_nj(:)
     CALL netcdf_err( nf90_inq_varid(ncid,"tglf_imp_p_flux", id_tglf_p_fluxi),id_tglf_p_fluxi) 
     CALL netcdf_err( nf90_get_var(ncid,id_tglf_p_fluxi,work_nj),id_tglf_p_fluxi)
     tglf_p_output(:,3) = work_nj(:)

     IF(ALLOCATED(tglf_e_output))DEALLOCATE(tglf_e_output)
     ALLOCATE(tglf_e_output(nj,3))
     CALL netcdf_err( nf90_inq_varid(ncid,"tglf_elc_e_flux", id_tglf_e_fluxe),id_tglf_e_fluxe) 
     CALL netcdf_err( nf90_get_var(ncid,id_tglf_e_fluxe,work_nj),id_tglf_e_fluxe)
     tglf_e_output(:,1) = work_nj(:)
     CALL netcdf_err( nf90_inq_varid(ncid,"tglf_ion_e_flux", id_tglf_e_fluxp),id_tglf_e_fluxp) 
     CALL netcdf_err( nf90_get_var(ncid,id_tglf_e_fluxp,work_nj),id_tglf_e_fluxp)
     tglf_e_output(:,2) = work_nj(:)
     CALL netcdf_err( nf90_inq_varid(ncid,"tglf_imp_e_flux", id_tglf_e_fluxi),id_tglf_e_fluxi) 
     CALL netcdf_err( nf90_get_var(ncid,id_tglf_e_fluxi,work_nj),id_tglf_e_fluxi)
     tglf_e_output(:,3) = work_nj(:)


     IF(ALLOCATED(tglf_m_output))DEALLOCATE(tglf_m_output)
     ALLOCATE(tglf_m_output(nj,3))
!     id_tglf_m_fluxe_er = -id_tglf_m_fluxe
      id_tglf_m_fluxe_er= izero

     CALL netcdf_err( nf90_inq_varid(ncid,"tglf_elc_m_flux", id_tglf_m_fluxe),id_tglf_m_fluxe_er) 
     IF(id_tglf_m_fluxe_er .NE. -1)THEN

        CALL netcdf_err( nf90_get_var(ncid,id_tglf_m_fluxe,work_nj),id_tglf_m_fluxe_er)
        tglf_m_output(:,1) = work_nj(:)
     ELSE
        tglf_m_output(:,1) = zeroc
     ENDIF

!     id_tglf_m_fluxp_er = - id_tglf_m_fluxp
     id_tglf_m_fluxp_er = izero
     CALL netcdf_err( nf90_inq_varid(ncid,"tglf_ion_m_flux", id_tglf_m_fluxp),id_tglf_m_fluxp_er) 
     IF(id_tglf_m_fluxp_er .NE. -1)THEN
        CALL netcdf_err( nf90_get_var(ncid,id_tglf_m_fluxp,work_nj),id_tglf_m_fluxp_er)
        tglf_m_output(:,2) = work_nj(:)
     ELSE
        tglf_m_output(:,2) = ZEROC
     endif

   !  id_tglf_m_fluxi_er = - id_tglf_m_fluxi
     id_tglf_m_fluxi_er = izero
     CALL netcdf_err( nf90_inq_varid(ncid,"tglf_imp_m_flux", id_tglf_m_fluxi),id_tglf_m_fluxi_er) 

     IF(id_tglf_m_fluxi_er .NE.  -1)THEN
        CALL netcdf_err( nf90_get_var(ncid,id_tglf_m_fluxi,work_nj),id_tglf_m_fluxi_er)
        tglf_m_output(:,3) = work_nj(:)
     ELSE
        tglf_m_output(:,3) = zeroc
     ENDIF

     IF( .NOT. ALLOCATED(stsource))THEN
        ntot = nion+dp4
        ALLOCATE(stsource(nion,nj),sion(nj,nion),       &
             sbcx(nj,nion),scx(nj,nion),dudtsv(ntot,nj))
             sion(:,:) = zeroc ;scx(:,:) = zeroc
             sbcx(:,:) = zeroc  ; dudtsv(:,:) = zeroc
     ENDIF
     IF(.NOT. ASSOCIATED(prtcl_src%srecom))ALLOCATE(prtcl_src%srecom(nion))

     CALL netcdf_err( nf90_inq_varid(ncid,"sion", id_sion),id_sion) 
     CALL netcdf_err( nf90_get_var(ncid,id_sion,sion),id_sion)
     CALL netcdf_err( nf90_inq_varid(ncid,"srecom", id_srecom),id_srecom) 
!     CALL netcdf_err( nf90_get_var(ncid,id_srecom,srecom),id_srecom)
     CALL netcdf_err( nf90_get_var(ncid,id_srecom,work_nj_nion),id_srecom)
     DO jj=1,nion
       prtcl_src%srecom(jj) = new_Vector(nj,work_nj_nion(1,jj))
     ENDDO
     CALL netcdf_err( nf90_inq_varid(ncid,"scx", id_scx),id_scx) 
     CALL netcdf_err( nf90_get_var(ncid,id_scx,scx),id_scx)
     CALL netcdf_err( nf90_inq_varid(ncid,"sbcx", id_sbcx),id_sbcx) 
     CALL netcdf_err( nf90_get_var(ncid,id_sbcx,sbcx),id_sbcx)
     CALL netcdf_err( nf90_inq_varid(ncid,"stsource", id_stsource),id_stsource) 
     CALL netcdf_err( nf90_get_var(ncid,id_stsource,stsource),id_stsource)
     CALL netcdf_err( nf90_inq_varid(ncid,"dudtsv", id_dudtsv),id_dudtsv) 
     CALL netcdf_err( nf90_get_var(ncid,id_dudtsv,dudtsv),id_dudtsv)
     CALL netcdf_err( nf90_inq_varid(ncid,"tGCNMf", id_tGCNMf),id_tGCNMf) 
     CALL netcdf_err( nf90_get_var(ncid,id_tGCNMf,tGCNMf),id_tGCNMf) 

     CALL netcdf_err( nf90_inq_varid(ncid,"time_bc", id_time_bc),id_time_bc) 
     CALL netcdf_err( nf90_get_var(ncid,id_time_bc,time_bc),id_time_bc) 


     CALL netcdf_err( nf90_inq_varid(ncid,"rmajavnpsi",id_rmajavnpsi),id_rmajavnpsi)
     CALL netcdf_err( nf90_get_var(ncid,id_rmajavnpsi,work_npsi),id_rmajavnpsi)
     dischg%rmajavnpsi = new_Vector(mhd_dat%npsi,work_npsi)

     !
     ! --- fast ion density
     ! 
     IF( .NOT. ALLOCATED(enbeam)) ALLOCATE(enbeam(nj,nbion))
     IF( .NOT. ALLOCATED(enbeam_tot)) ALLOCATE(enbeam_tot(nj))
     CALL netcdf_err( nf90_inq_varid(ncid,"enbeam",id_enbeam),id_enbeam)
     CALL netcdf_err( nf90_get_var(ncid,id_enbeam,enbeam),id_enbeam)
     enbeam_tot(:) = zeroc

     DO jn =1,nbion
        enbeam_tot(:) = enbeam_tot(:) + enbeam(:,jn)
     ENDDO
     IF(.NOT. ALLOCATED(enn))ALLOCATE(enn(nj,nneu))
     CALL netcdf_err( nf90_inq_varid(ncid,"enn",id_enn),id_enn)
     CALL netcdf_err( nf90_get_var(ncid,id_enn,enn),id0)
     IF(.NOT. ALLOCATED(ennw))ALLOCATE(ennw(nj,nneu))
     CALL netcdf_err( nf90_inq_varid(ncid,"ennw",id_ennw))
     CALL netcdf_err( nf90_get_var(ncid,id_ennw,ennw),id0)
     IF(.NOT. ALLOCATED(ennv))ALLOCATE(ennv(nj,nneu))
     CALL netcdf_err( nf90_inq_varid(ncid,"ennv",id_ennv))
     CALL netcdf_err( nf90_get_var(ncid,id_ennv,ennv),id0)
     IF(.NOT. ALLOCATED(volsn))ALLOCATE(volsn(nj,nneu))
     volsn(:,:) = zeroc
     CALL netcdf_err( nf90_inq_varid(ncid,"volsn",id_volsn))
     CALL netcdf_err( nf90_get_var(ncid,id_volsn,volsn),id0)

     CALL netcdf_err( nf90_inq_varid(ncid,"stfuse",id_stfuse),id_stfuse)
     CALL netcdf_err( nf90_get_var(ncid,id_stfuse,work_nj),id_stfuse)
     prtcl_src%stfuse  = new_vector(nj,work_nj)

     CALL netcdf_err( nf90_inq_varid(ncid,"sbfuse",id_sbfuse),id_sbfuse)
     CALL netcdf_err( nf90_get_var(ncid,id_sbfuse,work_nj),id_sbfuse)
     prtcl_src%sbfuse  = new_vector(nj,work_nj)

     CALL netcdf_err( nf90_inq_varid(ncid,"spellet",id_spellet),id_spellet)
     CALL netcdf_err( nf90_get_var(ncid,id_spellet,work_nj),id_spellet)
     prtcl_src%spellet  = new_vector(nj,work_nj)


     IF(.NOT. ALLOCATED(sbeame))ALLOCATE(sbeame(nj))
     sbeame = zeroc
     CALL netcdf_err( nf90_inq_varid(ncid,"sbeame",id_sbeame),id_sbeame)
     CALL netcdf_err( nf90_get_var(ncid,id_sbeame,sbeame),id_sbeame)
     !
     ! --- thermal ion source due to beams
     !

     IF(.NOT. ALLOCATED(sbeam))ALLOCATE(sbeam(nj,nbion))
 
     sbeam(:,:) = zeroc
     CALL netcdf_err( nf90_inq_varid(ncid,"sbeam",id_sbeam),id_sbeam)
     CALL netcdf_err( nf90_get_var(ncid,id_sbeam,sbeam),id_sbeam)


     IF(neut_beam%nbeams .GT. 0)THEN
        IF(neut_beam%nj_beam .LE. izero ) neut_beam%nj_beam = nj

        IF(ASSOCIATED(neut_beam%pbeam))DEALLOCATE(neut_beam%pbeam)
        ALLOCATE(neut_beam%pbeam(ke,neut_beam%nbeams))
        neut_beam%pbeam(:,:) = zeroc
        CALL netcdf_err( nf90_inq_varid(ncid,"pbeam",id_pbeam),id_pbeam)
        CALL netcdf_err( nf90_get_var(ncid,id_pbeam,neut_beam%pbeam),id_pbeam)


        CALL netcdf_err( nf90_inq_varid(ncid,"bptor",id_bptor),id_bptor)
        IF(.NOT. ALLOCATED(work_bptor)) &
                    ALLOCATE(work_bptor(1:neut_beam%nbeams))
        bptor(:) = 1.E-03_DP
        CALL netcdf_err( nf90_get_var(ncid,id_bptor,work_bptor),id_bptor)
        bptor(1:neut_beam%nbeams) = work_bptor(1:neut_beam%nbeams)
        DEALLOCATE(work_bptor)

        IF(ASSOCIATED(neut_beam%bneut))DEALLOCATE(neut_beam%bneut)
        ALLOCATE(neut_beam%bneut(ke,neut_beam%nbeams))
        neut_beam%bneut(:,:) = zeroc
        CALL netcdf_err( nf90_inq_varid(ncid,"bneut",id_bneut),id_bneut)
        CALL netcdf_err( nf90_get_var(ncid,id_bneut,neut_beam%bneut),id_bneut)

        IF(ASSOCIATED(neut_beam%bion))DEALLOCATE(neut_beam%bion)
        ALLOCATE(neut_beam%bion(ke,neut_beam%nbeams))
        neut_beam%bion(:,:) = zeroc
        CALL netcdf_err( nf90_inq_varid(ncid,"bion",id_bion),id_bion)
        CALL netcdf_err( nf90_get_var(ncid,id_bneut,neut_beam%bion),id_bion)

        IF(ASSOCIATED(neut_beam%fap))DEALLOCATE(neut_beam%fap)
        ALLOCATE(neut_beam%fap(ke,neut_beam%nbeams))
        neut_beam%fap(:,:) = zeroc
        CALL netcdf_err( nf90_inq_varid(ncid,"fap",id_fap),id_fap)
        CALL netcdf_err( nf90_get_var(ncid,id_fap,neut_beam%fap),id_fap)

        IF(ASSOCIATED(neut_beam%fwall))DEALLOCATE(neut_beam%fwall)
        ALLOCATE(neut_beam%fwall(ke,neut_beam%nbeams))
        neut_beam%fwall(:,:) = zeroc
        CALL netcdf_err( nf90_inq_varid(ncid,"fwall",id_fwall),id_fwall)
        CALL netcdf_err( nf90_get_var(ncid,id_fwall,neut_beam%fwall),id_fwall)


        IF(ASSOCIATED(neut_beam%forb))DEALLOCATE(neut_beam%forb)
        ALLOCATE(neut_beam%forb(ke,neut_beam%nbeams))
        neut_beam%forb(:,:) = zeroc
        CALL netcdf_err( nf90_inq_varid(ncid,"forb",id_forb),id_forb)
        CALL netcdf_err( nf90_get_var(ncid,id_forb,neut_beam%forb),id_forb)


        IF(ASSOCIATED(neut_beam%fber))DEALLOCATE(neut_beam%fber)
        ALLOCATE(neut_beam%fber(ke,neut_beam%nbeams))
        neut_beam%fber(:,:) = zeroc
        CALL netcdf_err( nf90_inq_varid(ncid,"fber",id_fber),id_fber)
        CALL netcdf_err( nf90_get_var(ncid,id_fber,neut_beam%fber),id_fber)


        IF(ASSOCIATED(neut_beam%fb00))DEALLOCATE(neut_beam%fb00)
        ALLOCATE(neut_beam%fb00(ke,neut_beam%nbeams))
        neut_beam%fb00(:,:) = zeroc
        CALL netcdf_err( nf90_inq_varid(ncid,"fb00",id_fb00),id_fb00)
        CALL netcdf_err( nf90_get_var(ncid,id_fb00,neut_beam%fb00),id_fb00)

        IF(ASSOCIATED(neut_beam%fb01))DEALLOCATE(neut_beam%fb01)
        ALLOCATE(neut_beam%fb01(ke,neut_beam%nbeams))
        neut_beam%fb01(:,:) = zeroc
        CALL netcdf_err( nf90_inq_varid(ncid,"fb01",id_fb01),id_fb01)
        CALL netcdf_err( nf90_get_var(ncid,id_fb01,neut_beam%fb01),id_fb01)


        IF(ASSOCIATED(neut_beam%fb10))DEALLOCATE(neut_beam%fb10)
        ALLOCATE(neut_beam%fb10(ke,neut_beam%nbeams))
        neut_beam%fb10(:,:) = zeroc
        CALL netcdf_err( nf90_inq_varid(ncid,"fb10",id_fb10),id_fb10)
        CALL netcdf_err( nf90_get_var(ncid,id_fb10,neut_beam%fb10),id_fb10)


        IF(ASSOCIATED(neut_beam%fb11))DEALLOCATE(neut_beam%fb11)
        ALLOCATE(neut_beam%fb11(ke,neut_beam%nbeams))
        neut_beam%fb11(:,:) = zeroc
        CALL netcdf_err( nf90_inq_varid(ncid,"fb11",id_fb11),id_fb11)
        CALL netcdf_err( nf90_get_var(ncid,id_fb11,neut_beam%fb11),id_fb11)


        IF(ASSOCIATED(neut_beam%wb00))DEALLOCATE(neut_beam%wb00)
        ALLOCATE(neut_beam%wb00(ke,neut_beam%nbeams))
        neut_beam%wb00(:,:) = zeroc
        CALL netcdf_err( nf90_inq_varid(ncid,"wb00",id_wb00),id_wb00)
        CALL netcdf_err( nf90_get_var(ncid,id_wb00,neut_beam%wb00),id_wb00)

        IF(ASSOCIATED(neut_beam%wb01))DEALLOCATE(neut_beam%wb01)
        ALLOCATE(neut_beam%wb01(ke,neut_beam%nbeams))
        neut_beam%wb01(:,:) = zeroc
        CALL netcdf_err( nf90_inq_varid(ncid,"wb01",id_wb01),id_wb01)
        CALL netcdf_err( nf90_get_var(ncid,id_wb01,neut_beam%wb01),id_wb01)


        IF(ASSOCIATED(neut_beam%wb10))DEALLOCATE(neut_beam%wb10)
        ALLOCATE(neut_beam%wb10(ke,neut_beam%nbeams))
        neut_beam%wb10(:,:) = zeroc
        CALL netcdf_err( nf90_inq_varid(ncid,"wb10",id_wb10),id_wb10)
        CALL netcdf_err( nf90_get_var(ncid,id_wb10,neut_beam%wb10),id_wb10)


        IF(ASSOCIATED(neut_beam%wb11))DEALLOCATE(neut_beam%wb11)
        ALLOCATE(neut_beam%wb11(ke,neut_beam%nbeams))
        neut_beam%wb11(:,:) = zeroc
        CALL netcdf_err( nf90_inq_varid(ncid,"wb11",id_wb11),id_wb11)
        CALL netcdf_err( nf90_get_var(ncid,id_wb11,neut_beam%wb11),id_wb11)


        IF(ASSOCIATED(neut_beam%ebeam))DEALLOCATE(neut_beam%ebeam)
        ALLOCATE(neut_beam%ebeam(ke,neut_beam%nbeams))
        neut_beam%ebeam(:,:) = zeroc
        CALL netcdf_err( nf90_inq_varid(ncid,"ebeam",id_ebeam),id_ebeam)
        CALL netcdf_err( nf90_get_var(ncid,id_ebeam,neut_beam%ebeam),id_ebeam)

        IF(ASSOCIATED(neut_beam%sb))DEALLOCATE(neut_beam%sb)
        ALLOCATE(neut_beam%sb(neut_beam%nj_beam,ke,neut_beam%nbeams))
        neut_beam%sb(:,:,:) = zeroc
        CALL netcdf_err( nf90_inq_varid(ncid,"sb",id_sb),id_sb)
        CALL netcdf_err( nf90_get_var(ncid,id_sb,neut_beam%sb),id_sb)

        IF(ASSOCIATED(neut_beam%qb))DEALLOCATE(neut_beam%qb)
        ALLOCATE(neut_beam%qb(neut_beam%nj_beam,ke,neut_beam%nbeams))
        neut_beam%qb(:,:,:) = zeroc
        CALL netcdf_err( nf90_inq_varid(ncid,"qb",id_qb),id_qb)
        CALL netcdf_err( nf90_get_var(ncid,id_qb,neut_beam%qb),id_qb)

        IF(ASSOCIATED(neut_beam%spb))DEALLOCATE(neut_beam%spb)
        ALLOCATE(neut_beam%spb(neut_beam%nj_beam,ke,neut_beam%nbeams))
        neut_beam%spb(:,:,:) = zeroc
        CALL netcdf_err( nf90_inq_varid(ncid,"spb",id_spb),id_spb)
        CALL netcdf_err( nf90_get_var(ncid,id_spb,neut_beam%spb),id_spb)

        IF(ASSOCIATED(neut_beam%spbr))DEALLOCATE(neut_beam%spbr)
        ALLOCATE(neut_beam%spbr(neut_beam%nj_beam,ke,neut_beam%nbeams))
        neut_beam%spbr(:,:,:) = zeroc
        CALL netcdf_err( nf90_inq_varid(ncid,"spbr",id_spbr),id_spbr)
        CALL netcdf_err( nf90_get_var(ncid,id_spbr,neut_beam%spbr),id_spbr)

        IF(ASSOCIATED(neut_beam%pb0))DEALLOCATE(neut_beam%pb0)
        ALLOCATE(neut_beam%pb0(neut_beam%nj_beam,ke,neut_beam%nbeams))
        neut_beam%pb0(:,:,:) = zeroc
        CALL netcdf_err( nf90_inq_varid(ncid,"pb0",id_pb0),id_pb0)
        CALL netcdf_err( nf90_get_var(ncid,id_pb0,neut_beam%pb0),id_pb0)

        IF(ASSOCIATED(neut_beam%angmpf))DEALLOCATE(neut_beam%angmpf)
        ALLOCATE(neut_beam%angmpf(neut_beam%nj_beam,ke,neut_beam%nbeams))
        neut_beam%angmpf(:,:,:) = zeroc
        CALL netcdf_err( nf90_inq_varid(ncid,"angmpf",id_angmpf),id_angmpf)
        CALL netcdf_err( nf90_get_var(ncid,id_angmpf,neut_beam%angmpf),id_angmpf)


       IF(ASSOCIATED(neut_beam%hibr))DEALLOCATE(neut_beam%hibr)
        ALLOCATE(neut_beam%hibr(neut_beam%nj_beam,ke,neut_beam%nbeams))
        neut_beam%hibr(:,:,:) = zeroc
        CALL netcdf_err( nf90_inq_varid(ncid,"hibr",id_hibr),id_hibr)
        CALL netcdf_err( nf90_get_var(ncid,id_hibr,neut_beam%hibr),id_hibr)

       IF(ASSOCIATED(neut_beam%hdep))DEALLOCATE(neut_beam%hdep)
        ALLOCATE(neut_beam%hdep(neut_beam%nj_beam,ke,neut_beam%nbeams))
        neut_beam%hdep(:,:,:) = zeroc
        CALL netcdf_err( nf90_inq_varid(ncid,"hdep",id_hdep),id_hdep)
        CALL netcdf_err( nf90_get_var(ncid,id_hdep,neut_beam%hdep),id_hdep)


       IF(ASSOCIATED(neut_beam%zeta))DEALLOCATE(neut_beam%zeta)
        ALLOCATE(neut_beam%zeta(neut_beam%nj_beam,ke,neut_beam%nbeams))
        neut_beam%zeta(:,:,:) = zeroc
        CALL netcdf_err( nf90_inq_varid(ncid,"zeta",id_zeta),id_zeta)
        CALL netcdf_err( nf90_get_var(ncid,id_zeta,neut_beam%zeta),id_zeta)

       IF(ASSOCIATED(neut_beam%fbcur))DEALLOCATE(neut_beam%fbcur)
        ALLOCATE(neut_beam%fbcur(ke,neut_beam%nbeams))
        neut_beam%fbcur(:,:) = zeroc
        id_fbcur_er = izero
        CALL netcdf_err( nf90_inq_varid(ncid,"fbcur",id_fbcur),id_fbcur_er)
        IF(id_fbcur_er .NE. -1)THEN
           CALL netcdf_err( nf90_get_var(ncid,id_fbcur,neut_beam%fbcur),id_fbcur)
        ENDIF

        IF(ASSOCIATED(neut_beam%prompt_pwr_in_plasma))     &
                                            DEALLOCATE(neut_beam%prompt_pwr_in_plasma)
        ALLOCATE(neut_beam%prompt_pwr_in_plasma(ke,neut_beam%nbeams))
        neut_beam%prompt_pwr_in_plasma(:,:) = zeroc
        id_prompt_nb_pwr_er = izero
        CALL netcdf_err( nf90_inq_varid(ncid,"prompt_nb_pwr",id_prompt_nb_pwr),id_prompt_nb_pwr_er)
        IF(id_prompt_nb_pwr_er .NE. -1)THEN
           CALL netcdf_err( nf90_get_var(ncid,id_prompt_nb_pwr,neut_beam%prompt_pwr_in_plasma),id_prompt_nb_pwr)


        ENDIF
         IF(neut_beam%nj_beam .LE. izero)neut_beam%nj_beam =nj

        IF(ASSOCIATED(neut_beam%rhog_beam))DEALLOCATE(neut_beam%rhog_beam)
        ALLOCATE(neut_beam%rhog_beam(neut_beam%nj_beam))
        neut_beam%rhog_beam  = get_values(rho_grid)
        id_rhog_beam_er = izero
        CALL netcdf_err( nf90_inq_varid(ncid,"rhog_beam",id_rhog_beam),id_rhog_beam_er)
        IF(id_rhog_beam_er .NE. -1)THEN

            CALL netcdf_err( nf90_get_var(ncid,id_rhog_beam,neut_beam%rhog_beam),id_rhog_beam)
        ENDIF


        IF(ASSOCIATED(neut_beam%hicm))DEALLOCATE(neut_beam%hicm)  
        ALLOCATE(neut_beam%hicm(neut_beam%nj_beam,ke_bm,neut_beam%nbeams,kcm))
        neut_beam%hicm(:,:,:,:) = zeroc
        IF(ALLOCATED(work3d))DEALLOCATE(work3d)
        ALLOCATE(work3d(neut_beam%nj_beam,ke_bm,neut_beam%nbeams))

        CALL netcdf_err( nf90_inq_varid(ncid,"hicme",id_hicme),id_hicme)
        CALL netcdf_err( nf90_get_var(ncid,id_hicme,work3d),id_hicme)
        neut_beam%hicm(:,:,:,1) = work3d(:,:,:)

        CALL netcdf_err( nf90_inq_varid(ncid,"hicmp1",id_hicmp1),id_hicmp1)
        CALL netcdf_err( nf90_get_var(ncid,id_hicmp1,work3d),id_hicmp1)
        neut_beam%hicm(:,:,:,2) = work3d(:,:,:)

        CALL netcdf_err( nf90_inq_varid(ncid,"hicmp2",id_hicmp2),id_hicmp2)
        CALL netcdf_err( nf90_get_var(ncid,id_hicmp2,work3d),id_hicmp2)
        neut_beam%hicm(:,:,:,3) = work3d(:,:,:)

        DEALLOCATE(work3d)


     ENDIF
 
     ! 
     ! ---  current profiles  :
     ! 

     CALL netcdf_err( nf90_inq_varid(ncid,"curden",id_curden),id_curden)

     CALL netcdf_err( nf90_get_var(ncid,id_curden,work_nj),id_curden)
     mhd_dat%curden = new_Vector(nj,work_nj)

     ! non fatal read of curpar (for older versions of code)
     id_curpar_flag = izero
     CALL netcdf_err( nf90_inq_varid(ncid,"curpar",id_curpar),id_curpar_flag)

     IF(id_curpar_flag .NE. -1)THEN 
        CALL netcdf_err( nf90_get_var(ncid,id_curpar,work_nj),id_curpar)
     ELSE
        work_nj(:) = zeroc
     ENDIF
     mhd_dat%curpar  = new_Vector(nj,work_nj)


     CALL netcdf_err( nf90_inq_varid(ncid,"curohm",id_curohm),id_curohm)
     CALL netcdf_err( nf90_get_var(ncid,id_curohm,work_nj),id_curohm)
     mhd_dat%curohm = new_Vector(nj,work_nj)

     CALL netcdf_err( nf90_inq_varid(ncid,"curboot",id_curboot),id_curboot)
     CALL netcdf_err( nf90_get_var(ncid,id_curboot,work_nj),id_curboot)
     mhd_dat%curboot = new_Vector(nj,work_nj)

     !beam current density
     !no mhd_dat%curbeam because beam current is fixed in gcnm code
     IF(.NOT. ALLOCATED(curbeam))THEN
        ALLOCATE(curbeam(nj))
        curbeam(:) = zeroc
     ENDIF
     CALL netcdf_err( nf90_inq_varid(ncid,"curbeam",id_curbeam),id_curbeam)
     CALL netcdf_err( nf90_get_var(ncid,id_curbeam,curbeam),id_curbeam)
 
     ! RF current density
     ! no mhd_dat%currf  because rf  current is fixed in gcnm code
     IF(.NOT. ALLOCATED(currf))THEN
        ALLOCATE(currf(nj))
        currf(:) = zeroc
     ENDIF
     ! < Jrf dot B/Bt0>        
     CALL netcdf_err( nf90_inq_varid(ncid,"currf",id_currf),id_currf)
     CALL netcdf_err( nf90_get_var(ncid,id_currf,currf),id_currf)


     ! toroidal electric field < E dot B /Bt0 >
     CALL netcdf_err( nf90_inq_varid(ncid,"etor",id_etor),id_etor)
     CALL netcdf_err( nf90_get_var(ncid,id_etor,work_nj),id_etor)
     profile%etor = new_Vector(nj,work_nj)

     ! rho*bp0*fcap*gcap*hcap
     CALL netcdf_err( nf90_inq_varid(ncid,"rbp",id_rbp),id_rbp)
     CALL netcdf_err( nf90_get_var(ncid,id_rbp,work_nj),id_rbp)
     mhd_dat%rbp = new_Vector(nj,work_nj)

     !ravgnpsi (1d from edge to mag axis)
     CALL netcdf_err( nf90_inq_varid(ncid,"ravgnpsi",id_ravgnpsi),id_ravgnpsi)
     CALL netcdf_err( nf90_get_var(ncid,id_ravgnpsi,work_npsi),id_ravgnpsi)
     mhd_dat%ravgnpsi = new_Vector(mhd_dat%npsi,work_npsi)

     !ravginpsi (1d from edge to mag axis)
     CALL netcdf_err( nf90_inq_varid(ncid,"ravginpsi",id_ravginpsi),id_ravginpsi)
     CALL netcdf_err( nf90_get_var(ncid,id_ravginpsi,work_npsi),id_ravginpsi)
     mhd_dat%ravginpsi = new_Vector(mhd_dat%npsi,work_npsi)

 

     !psi grid (1d from edge to mag axis)
     CALL netcdf_err( nf90_inq_varid(ncid,"psivalnpsi",id_psivalnpsi),id_psivalnpsi)
     CALL netcdf_err( nf90_get_var(ncid,id_psivalnpsi,work_npsi),id_psivalnpsi)
     mhd_dat%psivalnpsi = new_Vector(mhd_dat%npsi,work_npsi)

     !fpsi   on  psivalnpsi
     CALL netcdf_err( nf90_inq_varid(ncid,"fpsinpsi",id_fpsinpsi),id_fpsinpsi)
     CALL netcdf_err( nf90_get_var(ncid,id_fpsinpsi,work_npsi),id_fpsinpsi)
     mhd_dat%fpsinpsi = new_Vector(mhd_dat%npsi,work_npsi)

     !pprim   on  rho grid
     CALL netcdf_err( nf90_inq_varid(ncid,"pprim",id_pprim),id_pprim)
     CALL netcdf_err( nf90_get_var(ncid,id_pprim,work_nj),id_pprim)
     mhd_dat%pprim = new_Vector(nj,work_nj)

     ! ffprim   on  rho grid
     CALL netcdf_err( nf90_inq_varid(ncid,"ffprim",id_ffprim),id_ffprim)
     CALL netcdf_err( nf90_get_var(ncid,id_ffprim,work_nj),id_ffprim)
     mhd_dat%ffprim = new_Vector(nj,work_nj)

     !<Bp>   on  rho grid
     CALL netcdf_err( nf90_inq_varid(ncid,"bp",id_bp),id_bp)
     CALL netcdf_err( nf90_get_var(ncid,id_bp,work_nj),id_bp)
     mhd_dat%bp = new_Vector(nj,work_nj)


     !bprmaj on rmaj corresponding to rho grid
     CALL netcdf_err( nf90_inq_varid(ncid,"bprmaj",id_bprmaj),id_bprmaj)
     CALL netcdf_err( nf90_get_var(ncid,id_bprmaj,work_nj),id_bprmaj)
     mhd_dat%bprmaj = new_Vector(nj,work_nj)

     !btotrmaj on rmaj corresponding to rho grid
     CALL netcdf_err( nf90_inq_varid(ncid,"btotrmaj",id_btotrmaj),id_btotrmaj)
     CALL netcdf_err( nf90_get_var(ncid,id_btotrmaj,work_nj),id_btotrmaj)
     mhd_dat%btotrmaj = new_Vector(nj,work_nj)

     !zeff profiles
     CALL netcdf_err( nf90_inq_varid(ncid,"zeff",id_zeff),id_zeff)
     CALL netcdf_err( nf90_get_var(ncid,id_zeff,zeff),id_zeff)
     profile%zeff   = new_Vector(nj,zeff)

     ! angular rotation speed profile
     CALL netcdf_err( nf90_inq_varid(ncid,"angrot",id_angrot),id_angrot)
     CALL netcdf_err( nf90_get_var(ncid,id_angrot,work_nj),id_angrot)
     profile%angrot = new_Vector(nj,work_nj)

     ! poloidal rotation speed profile
     id_vpol_er = izero
     CALL netcdf_err( nf90_inq_varid(ncid,"vpol",id_vpol),id_vpol_er)
     IF(id_vpol_er .NE.  -1)THEN
        CALL netcdf_err( nf90_get_var(ncid,id_vpol,work_nj),id_vpol)
        profile%vpol = new_Vector(nj,work_nj)
     ELSE
        profile%vpol = zero_Vector(nj)
     ENDIF


     ! Nclass poloidal rotation speed profile
     id_vpol_nclass_er = izero
     CALL netcdf_err( nf90_inq_varid(ncid,"vpol_nclass",id_vpol_nclass),id_vpol_nclass_er)
     IF(id_vpol_nclass_er .NE.  -1)THEN
       CALL netcdf_err( nf90_get_var(ncid,id_vpol_nclass,work_nj),id_vpol_nclass)
       profile%vpol_nclass = new_Vector(nj,work_nj)
     ELSE
       profile%vpol_nclass = zero_Vector(nj)
     ENDIF

     ! parallel rotation speed profile
     id_vpar_er = izero
     CALL netcdf_err( nf90_inq_varid(ncid,"vpar",id_vpar),id_vpar_er)
     IF(id_vpar_er .NE.  -1)THEN
        CALL netcdf_err( nf90_get_var(ncid,id_vpar,work_nj),id_vpar)
        profile%vpar = new_Vector(nj,work_nj)
     ELSE
        profile%vpar = zero_Vector(nj)
     ENDIF


     ! Nclass poloidal rotation speed profile
     id_vpar_nclass_er = izero
     CALL netcdf_err( nf90_inq_varid(ncid,"vpar_nclass",id_vpar_nclass),id_vpar_nclass_er)
     IF(id_vpar_nclass_er .NE.  -1)THEN
       CALL netcdf_err( nf90_get_var(ncid,id_vpar_nclass,work_nj),id_vpar_nclass)
       profile%vpar_nclass = new_Vector(nj,work_nj)
     ELSE
       profile%vpar_nclass = zero_Vector(nj)
     ENDIF

     ! Nclass total radial electric field 
     id_er_tot_nclass_er = izero
     CALL netcdf_err( nf90_inq_varid(ncid,"er_tot_nclass",id_er_tot_nclass),id_er_tot_nclass_er)
     IF(id_er_tot_nclass_er .NE.  -1)THEN
       CALL netcdf_err( nf90_get_var(ncid,id_er_tot_nclass,work_nj),id_er_tot_nclass)
       profile%er_tot_nclass = new_Vector(nj,work_nj)
     ELSE
       profile%er_tot_nclass = zero_Vector(nj)
     ENDIF


     ! diffusivity matrix 
     CALL netcdf_err( nf90_inq_varid(ncid,"d",id_d),id_d)
     IF( ASSOCIATED (diffuse%dcoef))DEALLOCATE(diffuse%dcoef)
     ALLOCATE(diffuse%dcoef(ntot,ntot,nj))
     CALL netcdf_err( nf90_get_var(ncid,id_d,diffuse%dcoef),id_d)


     ! thermal diff. profiles, electron and ion
     CALL netcdf_err( nf90_inq_varid(ncid,"chieinv",id_chieinv),id_chieinv)
     CALL netcdf_err( nf90_get_var(ncid,id_chieinv,work_nj),id_chieinv)
     diffuse%chieinv  = new_Vector(nj,work_nj)

     CALL netcdf_err( nf90_inq_varid(ncid,"chiinv",id_chiinv),id_chiinv)
     CALL netcdf_err( nf90_get_var(ncid,id_chiinv,work_nj),id_chiinv)
     diffuse%chiinv  = new_Vector(nj,work_nj)

     ! --- ion neoclassical thermal conductivity
     CALL netcdf_err( nf90_inq_varid(ncid,"xkineo",id_xkineo),id_xkineo)
     CALL netcdf_err( nf90_get_var(ncid,id_xkineo,work_nj),id_xkineo)
     diffuse%xkineo  = new_Vector(nj,work_nj)

     ! --- electron  neoclassical thermal conductivity
     CALL netcdf_err( nf90_inq_varid(ncid,"xkeneo",id_xkeneo),id_xkeneo)
     CALL netcdf_err( nf90_get_var(ncid,id_xkeneo,work_nj),id_xkeneo)
     diffuse%xkeneo  = new_Vector(nj,work_nj)

     !---not  read in(these are pointers NOT vectors):
     !---diffuse%xkitot ,diffuse%xketot plus others

     !d(electron energy)/dt profile
     CALL netcdf_err( nf90_inq_varid(ncid,"dpedt",id_dpedt),id_dpedt)
     CALL netcdf_err( nf90_get_var(ncid,id_dpedt,work_nj),id_dpedt)
     wpdot%dpedt = new_Vector(nj,work_nj)

     ! --- d(ion energy)/dt profile
     IF(.NOT. ASSOCIATED(wpdot%dpidt))ALLOCATE(wpdot%dpidt(nion))
     CALL netcdf_err( nf90_inq_varid(ncid,"dpidt",id_dpidt),id_dpidt)
     CALL netcdf_err( nf90_get_var(ncid,id_dpidt,work_nj_nion),id_dpidt)
     DO jj = 1,nion
        wpdot%dpidt(jj)  = new_Vector(nj,work_nj_nion(:,jj))
     ENDDO

     ! --- electron conduction profile
     CALL netcdf_err( nf90_inq_varid(ncid,"qconde",id_qconde),id_qconde)
     CALL netcdf_err( nf90_get_var(ncid,id_qconde,work_nj),id_qconde)
     pwrden%qconde = new_Vector(nj,work_nj)

     ! --- ion conduction profile
     CALL netcdf_err( nf90_inq_varid(ncid,"qcondi",id_qcondi),id_qcondi)
     CALL netcdf_err( nf90_get_var(ncid,id_qcondi,work_nj),id_qcondi)
     pwrden%qcondi = new_Vector(nj,work_nj)

     ! --- electron convection profile
     CALL netcdf_err( nf90_inq_varid(ncid,"qconve",id_qconve),id_qconve)
     CALL netcdf_err( nf90_get_var(ncid,id_qconve,work_nj),id_qconve)
     pwrden%qconve = new_Vector(nj,work_nj)

     ! --- ion convection profile
     CALL netcdf_err( nf90_inq_varid(ncid,"qconvi",id_qconvi),id_qconvi)
     CALL netcdf_err( nf90_get_var(ncid,id_qconvi,work_nj),id_qconvi)
     pwrden%qconvi = new_Vector(nj,work_nj)

     ! --- beam electron profile
     CALL netcdf_err( nf90_inq_varid(ncid,"qbeame",id_qbeame),id_qbeame)
     CALL netcdf_err( nf90_get_var(ncid,id_qbeame,work_nj),id_qbeame)
     pwrden%qbeame = new_Vector(nj,work_nj)

     ! --- electron ion equilibration profile
     CALL netcdf_err( nf90_inq_varid(ncid,"qdelt",id_qdelt),id_qdelt)
     CALL netcdf_err( nf90_get_var(ncid,id_qdelt,work_nj),id_qdelt)
     pwrden%qdelt = new_Vector(nj,work_nj)

     ! --- beam ion profile
     CALL netcdf_err( nf90_inq_varid(ncid,"qbeami",id_qbeami),id_qbeami)
     CALL netcdf_err( nf90_get_var(ncid,id_qbeami,work_nj),id_qbeami)
     pwrden%qbeami = new_Vector(nj,work_nj)

     ! --- anomalous electron ion energy exchange term,
     ! --- due for example to drift ballooning mode
     CALL netcdf_err( nf90_inq_varid(ncid,"qexch",id_qexch),id_qexch)
     CALL netcdf_err( nf90_get_var(ncid,id_qexch,work_nj),id_qexch)
     pwrden%qexch = new_Vector(nj,work_nj)


     ! --- RF electron heating profile
     CALL netcdf_err( nf90_inq_varid(ncid,"qrfe",id_qrfe),id_qrfe)
     CALL netcdf_err( nf90_get_var(ncid,id_qrfe,work_nj),id_qrfe)
     pwrden%qrfe = new_Vector(nj,work_nj)

     ! --- RF ion heating profile
     CALL netcdf_err( nf90_inq_varid(ncid,"qrfi",id_qrfi),id_qrfi)
     CALL netcdf_err( nf90_get_var(ncid,id_qrfi,work_nj),id_qrfi)
     pwrden%qrfi = new_Vector(nj,work_nj)

     ! --- qione heating profile
     CALL netcdf_err( nf90_inq_varid(ncid,"qione",id_qione),id_qione)
     CALL netcdf_err( nf90_get_var(ncid,id_qione,work_nj),id_qione)
     pwrden%qione = new_Vector(nj,work_nj)

     ! --- qioni heating profile
     CALL netcdf_err( nf90_inq_varid(ncid,"qioni",id_qioni),id_qioni)
     CALL netcdf_err( nf90_get_var(ncid,id_qioni,work_nj),id_qioni)
     pwrden%qioni = new_Vector(nj,work_nj)



     ! --- qxc, ion heating profile
     CALL netcdf_err( nf90_inq_varid(ncid,"qcx",id_qcx),id_qcx)
     CALL netcdf_err( nf90_get_var(ncid,id_qcx,work_nj),id_qcx)
     pwrden%qcx   = new_Vector(nj,work_nj)

     ! --- 2d electron heating profile
     CALL netcdf_err( nf90_inq_varid(ncid,"qe2d",id_qe2d),id_qe2d)
     CALL netcdf_err( nf90_get_var(ncid,id_qe2d,work_nj),id_qe2d)
     pwrden%qe2d   = new_Vector(nj,work_nj)

     CALL netcdf_err( nf90_inq_varid(ncid,"qi2d",id_qi2d),id_qi2d)
     CALL netcdf_err( nf90_get_var(ncid,id_qi2d,work_nj),id_qi2d)
     pwrden%qi2d   = new_Vector(nj,work_nj)

     ! --- fusion electron heating  profile
     !      qfuse = qtfuse + qbfuse 
     !      (sum of thermal and fast ion progenated fusion heating to electrons
     !       fast part includes beam and alphas)
     !      (since qbfuse and qfuse are input, qtfuse is implicitely determined)
     CALL netcdf_err( nf90_inq_varid(ncid,"qfuse",id_qfuse),id_qfuse)
     CALL netcdf_err( nf90_get_var(ncid,id_qfuse,work_nj),id_qfuse)
     pwrden%qfuse  = new_Vector(nj,work_nj)

     ! --- fusion ion heating profile
     !      qfusi = qtfusi + qbfusi 
     !      (sum of thermal and fast ion progenated fusion heating to ions
     !       fast part includes beam and alphas)
     !      (since qbfusi and qfusi are input, qtfusi is implicitely determined)
     CALL netcdf_err( nf90_inq_varid(ncid,"qfusi",id_qfusi),id_qfusi)
     CALL netcdf_err( nf90_get_var(ncid,id_qfusi,work_nj),id_qfusi)
     pwrden%qfusi  = new_Vector(nj,work_nj)

     ! --- beam fusion electron heating profile
     ! --- (fraction of beam fusion energy deposited on 
     ! --- thermal electron distribution
     CALL netcdf_err( nf90_inq_varid(ncid,"qbfuse",id_qbfuse),id_qbfuse)
     CALL netcdf_err( nf90_get_var(ncid,id_qbfuse,work_nj),id_qbfuse)
     pwrden%qbfuse  = new_Vector(nj,work_nj)

     ! --- beam fusion ion heating profile
     ! --- (fraction of beam fusion energy deposited on thermal ion distribution)
     CALL netcdf_err( nf90_inq_varid(ncid,"qbfusi",id_qbfusi),id_qbfusi)
     CALL netcdf_err( nf90_get_var(ncid,id_qbfusi,work_nj),id_qbfusi)
     pwrden%qbfusi  = new_Vector(nj,work_nj)

     ! --- mag electron heating profile
     CALL netcdf_err( nf90_inq_varid(ncid,"qmag",id_qmag),id_qmag)
     CALL netcdf_err( nf90_get_var(ncid,id_qmag,work_nj),id_qmag)
     pwrden%qmag  = new_Vector(nj,work_nj)

     ! --- sawtooth electron heating profile
     CALL netcdf_err( nf90_inq_varid(ncid,"qsawe",id_qsawe),id_qsawe)
     CALL netcdf_err( nf90_get_var(ncid,id_qsawe,work_nj),id_qsawe)
     pwrden%qsawe  = new_Vector(nj,work_nj)
 
     CALL netcdf_err( nf90_inq_varid(ncid,"qsawi",id_qsawi),id_qsawi)
     CALL netcdf_err( nf90_get_var(ncid,id_qsawi,work_nj),id_qsawi)
     pwrden%qsawi  = new_Vector(nj,work_nj)

     ! --- radiated power density
     CALL netcdf_err( nf90_inq_varid(ncid,"qrad",id_qrad),id_qrad)
     CALL netcdf_err( nf90_get_var(ncid,id_qrad,work_nj),id_qrad)
     pwrden%qrad  = new_Vector(nj,work_nj)

     ! components of qrad:
     IF(.NOT. ASSOCIATED(pwrden%brems_nions)) &
                               ALLOCATE(pwrden%brems_nions(nion))
     id_brems_flag =0                ! prevent error exit if not present:
     CALL netcdf_err( nf90_inq_varid(ncid,"brems_nions", id_brems_nions),id_brems_flag) 
     IF(id_brems_flag == -1)THEN     ! variable not present in data set,not fatal
        work_nj_nion(1:nj,1:nion) = zeroc
     ELSE                            ! variable is in data set
        CALL netcdf_err( nf90_get_var(ncid,id_brems_nions,work_nj_nion),id_brems_nions)
     ENDIF
     DO jj =1,nion
         pwrden%brems_nions(jj)  = new_Vector(nj,work_nj_nion(1,jj))
     ENDDO



     ! nclass pinch velocity
     IF(.NOT. ASSOCIATED(diffuse%vpinch_nclass)) &
                               ALLOCATE(diffuse%vpinch_nclass(nion+1))
     id_vpinch_flag =0                ! prevent error exit if not present:
     CALL netcdf_err( nf90_inq_varid(ncid,"vpinch_nclass", id_vpinch_nions),id_vpinch_flag) 
     IF(id_vpinch_flag == -1)THEN     ! variable not present in data set,not fatal
        work_nj_nionp1(1:nj,1:nion+1) = zeroc
     ELSE                            ! variable is in data set
        CALL netcdf_err( nf90_get_var(ncid,id_vpinch_nions,work_nj_nionp1),id_vpinch_nions)
     ENDIF
     DO jj =1,nion+1
         diffuse%vpinch_nclass(jj)  = new_Vector(nj,work_nj_nionp1(1,jj))
     ENDDO

 
    









     ! --- qohm,ohmic heating profile
     CALL netcdf_err( nf90_inq_varid(ncid,"qohm",id_qohm),id_qohm)
     CALL netcdf_err( nf90_get_var(ncid,id_qohm,work_nj),id_qohm)
     pwrden%qohm  = new_Vector(nj,work_nj)


     id_omegale_flag = izero
     CALL netcdf_err( nf90_inq_varid(ncid,"omegale", id_omegale),id_omegale_flag) 
     IF(id_omegale_flag .NE. -1)THEN 
        CALL netcdf_err( nf90_get_var(ncid,id_omegale,work_nj),id_omegale)
     ELSE
        work_nj(:) = zeroc
     ENDIF
     pwrden%omegale  = new_Vector(nj,work_nj)


     id_qangce_flag = izero
     CALL netcdf_err( nf90_inq_varid(ncid,"qangce", id_qangce),id_qangce_flag) 
     IF(id_qangce_flag .NE. -1)THEN 
        CALL netcdf_err( nf90_get_var(ncid,id_qangce,work_nj),id_qangce)
     ELSE
        work_nj(:) = zeroc
     ENDIF
     pwrden%qangce  = new_Vector(nj,work_nj)

     id_qomegapi_flag = izero
     CALL netcdf_err( nf90_inq_varid(ncid,"qomegapi", id_qomegapi),id_qomegapi_flag) 
     IF(id_qomegapi_flag .NE. -1)THEN 
        CALL netcdf_err( nf90_get_var(ncid,id_qomegapi,work_nj),id_qomegapi)
     ELSE
        work_nj(:) = zeroc
     ENDIF
     pwrden%qomegapi  = new_Vector(nj,work_nj)

     id_sprcxre_flag = izero
     CALL netcdf_err( nf90_inq_varid(ncid,"sprcxre", id_sprcxre),id_sprcxre_flag) 
     IF(id_sprcxre_flag .NE. -1)THEN 
        CALL netcdf_err( nf90_get_var(ncid,id_sprcxre,work_nj),id_sprcxre)
     ELSE
        work_nj(:) = zeroc
     ENDIF
     pwrden%sprcxre  = new_Vector(nj,work_nj)


     id_sprcxree_flag = izero
     CALL netcdf_err( nf90_inq_varid(ncid,"sprcxree", id_sprcxree),id_sprcxree_flag) 
     IF(id_sprcxree_flag .NE. -1)THEN 
        CALL netcdf_err( nf90_get_var(ncid,id_sprcxree,work_nj),id_sprcxree)
     ELSE
        work_nj(:) = zeroc
     ENDIF
     pwrden%sprcxree  = new_Vector(nj,work_nj)



     id_spreimpe_flag = izero
     CALL netcdf_err( nf90_inq_varid(ncid,"spreimpe", id_spreimpe),id_spreimpe_flag) 
     IF(id_spreimpe_flag .NE. -1)THEN 
        CALL netcdf_err( nf90_get_var(ncid,id_spreimpe,work_nj),id_spreimpe)
     ELSE
        work_nj(:) = zeroc
     ENDIF
     pwrden%spreimpe  = new_Vector(nj,work_nj)


     ! --- average minor radius
     CALL netcdf_err( nf90_inq_varid(ncid,"rminavnpsi",id_rminavnpsi),id_rminavnpsi)
     CALL netcdf_err( nf90_get_var(ncid,id_rminavnpsi,work_npsi),id_rminavnpsi)
     dischg%rminavnpsi = new_Vector(mhd_dat%npsi,work_npsi)

     CALL netcdf_err( nf90_inq_varid(ncid,"psivolpnpsi",id_psivolpnpsi),id_psivolpnpsi)
     CALL netcdf_err( nf90_get_var(ncid,id_psivolpnpsi,work_npsi),id_psivolpnpsi)
     dischg%psivolpnpsi = new_Vector(mhd_dat%npsi,work_npsi)

     ! ---  elongation
     CALL netcdf_err( nf90_inq_varid(ncid,"elongxnpsi",id_elongxnpsi),id_elongxnpsi)
     CALL netcdf_err( nf90_get_var(ncid,id_elongxnpsi,work_npsi),id_elongxnpsi)
     dischg%elongxnpsi = new_Vector(mhd_dat%npsi,work_npsi)

     ! --- triangularity
     CALL netcdf_err( nf90_inq_varid(ncid,"triangnpsi_u",id_triangnpsi_u),id_triangnpsi_u)
     CALL netcdf_err( nf90_get_var(ncid,id_triangnpsi_u,work_npsi),id_triangnpsi_u)
     dischg%triangnpsi_u = new_Vector(mhd_dat%npsi,work_npsi)
     CALL netcdf_err( nf90_inq_varid(ncid,"triangnpsi_l",id_triangnpsi_l),id_triangnpsi_l)
     CALL netcdf_err( nf90_get_var(ncid,id_triangnpsi_l,work_npsi),id_triangnpsi_l)
     dischg%triangnpsi_l = new_Vector(mhd_dat%npsi,work_npsi)
     CALL netcdf_err( nf90_inq_varid(ncid,"pindentnpsi",id_pindentnpsi),id_pindentnpsi)
     CALL netcdf_err( nf90_get_var(ncid,id_pindentnpsi,work_npsi),id_pindentnpsi)
     dischg%pindentnpsi  = new_Vector(mhd_dat%npsi,work_npsi)

     ! --- surface area
     CALL netcdf_err( nf90_inq_varid(ncid,"sfareanpsi",id_sfareanpsi),id_sfareanpsi)
     CALL netcdf_err( nf90_get_var(ncid,id_sfareanpsi,work_npsi),id_sfareanpsi)
     dischg%sfareanpsi = new_Vector(mhd_dat%npsi,work_npsi)

     ! --- cross-sectional area
     CALL netcdf_err( nf90_inq_varid(ncid,"cxareanpsi",id_cxareanpsi),id_cxareanpsi)
     CALL netcdf_err( nf90_get_var(ncid,id_cxareanpsi,work_npsi),id_cxareanpsi)
     dischg%cxareanpsi = new_Vector(mhd_dat%npsi,work_npsi)

     ! --- flux surface average grad rho
     CALL netcdf_err( nf90_inq_varid(ncid,"grho1npsi",id_grho1npsi),id_grho1npsi)
     CALL netcdf_err( nf90_get_var(ncid,id_grho1npsi,work_npsi),id_grho1npsi)
     dischg%grho1npsi = new_Vector(mhd_dat%npsi,work_npsi)


     ! --- flux surface average (grad rho)**2
     CALL netcdf_err( nf90_inq_varid(ncid,"grho2npsi",id_grho2npsi),id_grho2npsi)
     CALL netcdf_err( nf90_get_var(ncid,id_grho2npsi,work_npsi),id_grho2npsi)
     dischg%grho2npsi = new_Vector(mhd_dat%npsi,work_npsi)

     ! --- beam torque density - new input added here so that backward
     !     compatibility is maintained
     CALL netcdf_err( nf90_inq_varid(ncid,"storqueb",id_storqueb),id_storqueb)
     CALL netcdf_err( nf90_get_var(ncid,id_storqueb ,storqueb),id_storqueb)



     !     in the future we may introduce time dependent metrics.
     !     then dfdt (= d/dt FCAP ,etc) will be calculated internally in this code
     !     by reading in fcap_bc, etc.
     !     for now just assume time independent for each
     !     slice that the solver is called.
     IF(.NOT. ALLOCATED(dfdt))ALLOCATE(dfdt(nj))
     IF(.NOT. ALLOCATED(dgdt))ALLOCATE(dgdt(nj))
     IF(.NOT. ALLOCATED(dhdt))ALLOCATE(dhdt(nj))

     !
     CALL netcdf_err( nf90_inq_varid(ncid,"dfdt",id_dfdt),id_dfdt)
     rcode = nf90_get_var(ncid,id_dfdt ,dfdt )
     IF(rcode /= nf90_noerr)dfdt(:) = zeroc

     CALL netcdf_err( nf90_inq_varid(ncid,"dgdt",id_dgdt),id_dgdt)
     rcode =  nf90_get_var(ncid,id_dgdt ,dgdt )
     IF(rcode /= nf90_noerr)dgdt(:) = zeroc

     CALL netcdf_err( nf90_inq_varid(ncid,"dhdt",id_dhdt),id_dhdt)
     rcode =   nf90_get_var(ncid,id_dhdt ,dhdt )
     IF(rcode /= nf90_noerr)dhdt(:) = zeroc


     !**_bc means ** is for  boundary condtions 
     !we give the profile over the entire rho grid but most often
     !only the edge values are actually used to get the 
     !boundary conditions. SOme model such as glf23 may have
     !boundaries that are set as far in as rho 0.6 so the **_bc
     !profiles must have appropriate values at least in that range.

     !boundary condition total current :
     CALL netcdf_err( nf90_inq_varid(ncid,"totcur_bc",id_totcur_bc),id_totcur_bc)
     CALL netcdf_err( nf90_get_var(ncid,id_totcur_bc ,totcur_bc ),id_totcur_bc)
     !boundary condition loop voltage :
     CALL netcdf_err( nf90_inq_varid(ncid,"vloop_bc",id_vloop_bc),id_vloop_bc)
     CALL netcdf_err( nf90_get_var(ncid,id_vloop_bc ,vloop_bc ),id_vloop_bc)


     !boundary indicators (determine rho point at which bc is applied)
     CALL netcdf_err( nf90_inq_varid(ncid,"fix_edge_te_bc",id_fix_edge_te_bc),id_fix_edge_te_bc)
     CALL netcdf_err( nf90_get_var(ncid,id_fix_edge_te_bc ,fix_edge_te_bc ),id_fix_edge_te_bc)

     CALL netcdf_err( nf90_inq_varid(ncid,"fix_edge_ti_bc",id_fix_edge_ti_bc),id_fix_edge_ti_bc)
     CALL netcdf_err( nf90_get_var(ncid,id_fix_edge_ti_bc ,fix_edge_ti_bc ),id_fix_edge_ti_bc)

     CALL netcdf_err( nf90_inq_varid(ncid,"fix_edge_rot_bc",id_fix_edge_rot_bc),id_fix_edge_rot_bc)
     CALL netcdf_err( nf90_get_var(ncid,id_fix_edge_rot_bc ,fix_edge_rot_bc ),id_fix_edge_rot_bc)

     IF(.NOT. ALLOCATED(fix_edge_ni_bc))ALLOCATE(fix_edge_ni_bc(nion))
     CALL netcdf_err( nf90_inq_varid(ncid,"fix_edge_ni_bc",id_fix_edge_ni_bc),id_fix_edge_ni_bc)
     CALL netcdf_err( nf90_get_var(ncid,id_fix_edge_ni_bc ,fix_edge_ni_bc ),id_fix_edge_ni_bc)

     ! convert indecies input as grid points equivalent  normalized rho values
     CALL convert_edge_indecies

     ! boundary condition TE,TI kev,ene 1/m**3
     CALL netcdf_err( nf90_inq_varid(ncid,"te_bc",id_te_bc),id_te_bc)
     CALL netcdf_err( nf90_get_var(ncid,id_te_bc,te_bc ),id_te_bc)

     CALL netcdf_err( nf90_inq_varid(ncid,"ti_bc",id_ti_bc),id_ti_bc)
     CALL netcdf_err( nf90_get_var(ncid,id_ti_bc,ti_bc ),id_ti_bc)


     CALL netcdf_err( nf90_inq_varid(ncid,"ene_bc",id_ene_bc),id_ene_bc)
     CALL netcdf_err( nf90_get_var(ncid,id_ene_bc,ene_bc ),id_ene_bc)


     CALL netcdf_err( nf90_inq_varid(ncid,"zeff_bc",id_zeff_bc),id_zeff_bc)
     CALL netcdf_err( nf90_get_var(ncid,id_zeff_bc,zeff_bc ),id_zeff_bc)

     ! boundary condition,toroidal rotation rad/sec:
     CALL netcdf_err( nf90_inq_varid(ncid,"angrot_bc",id_angrot_bc),id_angrot_bc)
     CALL netcdf_err( nf90_get_var(ncid,id_angrot_bc,angrot_bc ),id_angrot_bc)

     ! boundary condition,ion densities,1/m**3:
     CALL netcdf_err( nf90_inq_varid(ncid,"en_bc",id_en_bc),id_en_bc)
     CALL netcdf_err( nf90_get_var(ncid,id_en_bc,en_bc ),id_en_bc)

     !ion flux  boundary condition,1/(m^2 sec):
     CALL netcdf_err( nf90_inq_varid(ncid,"flux_bc",id_flux_bc),id_flux_bc)
     CALL netcdf_err( nf90_get_var(ncid,id_flux_bc,flux_bc ),id_flux_bc)

     ! --- primary, impurity charge,charge square
     CALL netcdf_err( nf90_inq_varid(ncid,"z",id_z),id_z)
     CALL netcdf_err( nf90_get_var(ncid,id_z,z ),id_z)
     CALL netcdf_err( nf90_inq_varid(ncid,"zsq",id_zsq),id_zsq)
     CALL netcdf_err( nf90_get_var(ncid,id_zsq,zsq ),id_zsq)

     DO jj=1,nion
      profile%zsq(jj)          = zero_Vector(nj)  ! ifort wants it this way
      profile%z(jj)            = zero_Vector(nj)
      profile%zsq(jj)%data(:)  = zsq(:,jj)
      profile%z(jj)%data(:)    = z(:,jj)
     ENDDO

     !fast ion stored energy density KEV/m**3 ' 
     CALL netcdf_err( nf90_inq_varid(ncid,"wbeam",id_wbeam),id_wbeam)
     CALL netcdf_err( nf90_get_var(ncid,id_wbeam,wbeam ),id_wbeam)
     profile%wbeam = new_Vector(nj,wbeam)
     !fast alpha particle stored energy density KEV/m**3 ' 
     CALL netcdf_err( nf90_inq_varid(ncid,"walp",id_walp),id_walp)
     CALL netcdf_err( nf90_get_var(ncid,id_walp,walp ),id_walp)
     WHERE(walp < zeroc)walp= zeroc  
     profile%walp = new_Vector(nj,walp) 

     !fast alpha density 1/m**3 ' 
     CALL netcdf_err( nf90_inq_varid(ncid,"enalp",id_enalp),id_enalp)
     CALL netcdf_err( nf90_get_var(ncid,id_enalp,enalp ),id_enalp)


     ! rate of change of ene  1/(m**3 sec)
     CALL netcdf_err( nf90_inq_varid(ncid,"dnedt",id_dnedt),id_dnedt)
     CALL netcdf_err( nf90_get_var(ncid,id_dnedt,work_nj),id_dnedt)
     wpdot%dnedt  =  new_Vector(nj,work_nj)


     !horizontal inverse aspect ratio = (rmax-rmin)/(rmax+rmin)
     CALL netcdf_err( nf90_inq_varid(ncid,"eps",id_eps),id_eps)
     CALL netcdf_err( nf90_get_var(ncid,id_eps,eps),id_eps)

     !rcap = < R>,m
     CALL netcdf_err( nf90_inq_varid(ncid,"rcap",id_rcap),id_rcap)
     CALL netcdf_err( nf90_get_var(ncid,id_rcap,work_nj),id_rcap)
     mhd_dat%rcap = new_Vector(nj,work_nj)


     !rcapi = < 1/R>,m
     CALL netcdf_err( nf90_inq_varid(ncid,"rcapi",id_rcapi),id_rcapi)
     CALL netcdf_err( nf90_get_var(ncid,id_rcapi,work_nj),id_rcapi)
     mhd_dat%rcapi = new_Vector(nj,work_nj)


     !r2cap  = <R0**2/R**2>
     CALL netcdf_err( nf90_inq_varid(ncid,"r2cap",id_r2cap),id_r2cap)
     CALL netcdf_err( nf90_get_var(ncid,id_r2cap,work_nj),id_r2cap)
     mhd_dat%r2cap = new_Vector(nj,work_nj)


     !r2capi = <R**2>       M**2'
     CALL netcdf_err( nf90_inq_varid(ncid,"r2capi",id_r2capi),id_r2capi)
     CALL netcdf_err( nf90_get_var(ncid,id_r2capi,work_nj),id_r2capi)
     mhd_dat%r2capi = new_Vector(nj,work_nj)

     !     some quantities required for trapped particle fraction, and
     !     neocl;assical transport

     IF(.NOT. ALLOCATED(xips))ALLOCATE(xips(nj))
     IF(.NOT. ALLOCATED(xips0))ALLOCATE(xips0(nj))
     IF(.NOT. ALLOCATED(xhm2))ALLOCATE(xhm2(nj))
     IF(.NOT. ALLOCATED(xhm20))ALLOCATE(xhm20(nj))
     IF(.NOT. ALLOCATED(xi11))ALLOCATE(xi11(nj))
     IF(.NOT. ALLOCATED(xi110))ALLOCATE(xi110(nj))
     IF(.NOT. ALLOCATED(xi33))ALLOCATE(xi33(nj))
     IF(.NOT. ALLOCATED(xi330))ALLOCATE(xi330(nj))


     !xhm2 = < (B total/ B axis)**2 > (=1 for circular plasmas)
     CALL netcdf_err( nf90_inq_varid(ncid,"xhm2",id_xhm2),id_xhm2)
     CALL netcdf_err( nf90_get_var(ncid,id_xhm2,xhm2),id_xhm2)

     !xi11  = 1.95 sqrt(eps)for circular plasmas (see doc for complete desc.)
     CALL netcdf_err( nf90_inq_varid(ncid,"xi11",id_xi11),id_xi11)
     CALL netcdf_err( nf90_get_var(ncid,id_xi11,xi11),id_xi11)

     !xi33 ( = 1.95 sqrt(eps)
     CALL netcdf_err( nf90_inq_varid(ncid,"xi33",id_xi33),id_xi33)
     CALL netcdf_err( nf90_get_var(ncid,id_xi33,xi33),id_xi33)

     !xips = <(Baxis/B)**2)> - 1./(<(B/Baxis)**2> )
     CALL netcdf_err( nf90_inq_varid(ncid,"xips",id_xips),id_xips)
     CALL netcdf_err( nf90_get_var(ncid,id_xips,xips),id_xips)

     xhm20(:) = xhm2(:)
     xi110(:) = xi11(:)
     xi330(:) = xi33(:)
     xips0(:) = xips(:)



     ! --- plasma boundary

     !          CALL netcdf_err( nf90_inq_varid(ncid,"nplasbdry",id_nplasbdry),id_nplasbdry)
     !          CALL netcdf_err( nf90_get_var(ncid,id_nplasbdry,dischg%nplasbdry),id_nplasbdry)
     CALL netcdf_err( nf90_inq_varid(ncid,"rplasbdry",id_rplasbdry),id_rplasbdry)
     CALL netcdf_err( nf90_get_var(ncid,id_rplasbdry,work_nplasbdry),id_rplasbdry)
     dischg%rplasbdry  = new_Vector(dischg%nplasbdry,work_nplasbdry)

     CALL netcdf_err( nf90_inq_varid(ncid,"zplasbdry",id_zplasbdry),id_zplasbdry)
     CALL netcdf_err( nf90_get_var(ncid,id_zplasbdry,work_nplasbdry),id_zplasbdry)
     dischg%zplasbdry  = new_Vector(dischg%nplasbdry,work_nplasbdry)


     ! ---  limiter for plasma:

        ! CALL netcdf_err( nf90_inq_varid(ncid,"nlimiter",id_nlimiter),id_nlimiter)
        ! CALL netcdf_err( nf90_get_var(ncid,id_nlimiter,dischg%nlimiter),id_nlimiter)
 
     CALL netcdf_err( nf90_inq_varid(ncid,"rlimiter",id_rlimiter),id_rlimiter)
     CALL netcdf_err( nf90_get_var(ncid,id_rlimiter,work_nlimiter),id_rlimiter)
 
     dischg%rlimiter = new_Vector(dischg%nlimiter,work_nlimiter)

     CALL netcdf_err( nf90_inq_varid(ncid,"zlimiter",id_zlimiter),id_zlimiter)
     CALL netcdf_err( nf90_get_var(ncid,id_zlimiter,work_nlimiter),id_zlimiter)
     dischg%zlimiter  = new_Vector(dischg%nlimiter,work_nlimiter)


     flag = 0
     CALL netcdf_err( nf90_inq_varid(ncid,"neutr_ddn_th", id_neutr_ddn_th),flag) 
     IF(flag == -1)THEN
        fus_prod%neutr_ddn_th = zero_Vector(nj)
     ELSE
        CALL netcdf_err( nf90_get_var(ncid,id_neutr_ddn_th,work_nj),id_neutr_ddn_th)
        fus_prod%neutr_ddn_th = new_Vector(nj,work_nj)
     ENDIF

     flag  = 0
     CALL netcdf_err( nf90_inq_varid(ncid,"neutr_ddn_beam_beam", id_neutr_ddn_beam_beam),flag) 
     IF(flag == -1)THEN
       fus_prod%neutr_ddn_beam_beam = zero_Vector(nj)
     ELSE
        CALL netcdf_err( nf90_get_var(ncid,id_neutr_ddn_beam_beam,work_nj),id_neutr_ddn_beam_beam)
        fus_prod%neutr_ddn_beam_beam = new_Vector(nj,work_nj)
     ENDIF

     flag = 0
     CALL netcdf_err( nf90_inq_varid(ncid,"neutr_ddn_beam_thermal", id_neutr_ddn_beam_thermal),flag) 
     IF(flag == -1)THEN
        fus_prod%neutr_ddn_beam_thermal = zero_Vector(nj)
     ELSE
        CALL netcdf_err( nf90_get_var(ncid,id_neutr_ddn_beam_thermal,work_nj),id_neutr_ddn_beam_thermal)
        fus_prod%neutr_ddn_beam_thermal = new_Vector(nj,work_nj)
     ENDIF

     flag = 0
     CALL netcdf_err( nf90_inq_varid(ncid,"neutr_ddn_knock", id_neutr_ddn_knock),flag) 
     IF(flag == -1)THEN
        fus_prod%neutr_ddn_knock = zero_Vector(nj)
     ELSE
        CALL netcdf_err( nf90_get_var(ncid,id_neutr_ddn_knock,work_nj),id_neutr_ddn_knock)
        fus_prod%neutr_ddn_knock = new_Vector(nj,work_nj)
     ENDIF

     flag = 0
     CALL netcdf_err( nf90_inq_varid(ncid,"neutr_ddn_tot", id_neutr_ddn_tot),flag) 
     IF(flag == -1)THEN
        fus_prod%neutr_ddn_tot = zero_Vector(nj)
     ELSE
        CALL netcdf_err( nf90_get_var(ncid,id_neutr_ddn_tot,work_nj),id_neutr_ddn_tot)
        fus_prod%neutr_ddn_tot = new_Vector(nj,work_nj)
     ENDIF

     flag = 0
     CALL netcdf_err( nf90_inq_varid(ncid,"total_neutr_ddn_th", id_total_neutr_ddn_th),flag)
     IF(flag  == -1)THEN
        fus_prod%total_neutr_ddn_th = zeroc
     ELSE
        CALL netcdf_err( nf90_get_var(ncid,id_total_neutr_ddn_th,fus_prod%total_neutr_ddn_th),id_total_neutr_ddn_th)
     ENDIF

     flag = 0
     CALL netcdf_err( nf90_inq_varid(ncid,"total_neutr_ddn_beam_beam", id_total_neutr_ddn_beam_beam),flag) 
     IF(flag  == -1)THEN
         fus_prod%total_neutr_ddn_beam_beam = zeroc
     ELSE
         CALL netcdf_err( nf90_get_var(ncid,id_total_neutr_ddn_beam_beam,fus_prod%total_neutr_ddn_beam_beam), &
                                            id_total_neutr_ddn_beam_beam)
     ENDIF

     flag = 0
     CALL netcdf_err( nf90_inq_varid(ncid,"total_neutr_ddn_beam_thermal", id_total_neutr_ddn_beam_thermal),flag)
     IF(flag  == -1)THEN
          fus_prod%total_neutr_ddn_beam_thermal = zeroc
     ELSE
          CALL netcdf_err( nf90_get_var(ncid,id_total_neutr_ddn_beam_thermal,fus_prod%total_neutr_ddn_beam_thermal), &
                                             id_total_neutr_ddn_beam_thermal)
     ENDIF

     flag = 0
     CALL netcdf_err( nf90_inq_varid(ncid,"total_neutr_ddn_knock", id_total_neutr_ddn_knock),flag) 
     IF(flag  == -1)THEN
          fus_prod%total_neutr_ddn_knock = zeroc
     ELSE
          CALL netcdf_err( nf90_get_var(ncid,id_total_neutr_ddn_knock,fus_prod%total_neutr_ddn_knock),id_total_neutr_ddn_knock)
     ENDIF

     flag = 0
     CALL netcdf_err( nf90_inq_varid(ncid,"total_neutr_ddn", id_total_neutr_ddn),flag) 
     IF(flag  == -1)THEN
          fus_prod%total_neutr_ddn = zeroc
     ELSE
          CALL netcdf_err( nf90_get_var(ncid,id_total_neutr_ddn,fus_prod%total_neutr_ddn),id_total_neutr_ddn)
     ENDIF

 
! --- clear out old values if present
! --- replace with new ones read in (even if new ones are zero)
     ok = delete_vector_nf(fus_prod%p_ddt)
     ok = delete_vector_nf(fus_prod%p_he3d_th)
     ok = delete_vector_nf(fus_prod%p_he3d_beam_thermal)
     flag = 0
     CALL netcdf_err( nf90_inq_varid(ncid,"ddpt", id_p_ddt),flag) 
 
     IF(flag == -1)THEN
        fus_prod%p_ddt = zero_Vector(nj)
     ELSE
        CALL netcdf_err( nf90_get_var(ncid,id_p_ddt,work_nj),id_p_ddt)
        fus_prod%p_ddt = new_Vector(nj,work_nj)
     ENDIF
 
     flag = 0
     CALL netcdf_err( nf90_inq_varid(ncid,"ddpt_tot", id_total_p_ddt),flag) 
     IF(flag  == -1)THEN
          fus_prod%total_p_ddt = zeroc
     ELSE
          CALL netcdf_err( nf90_get_var(ncid,id_total_p_ddt,fus_prod%total_p_ddt),id_total_p_ddt)
     ENDIF



     flag = 0
     CALL netcdf_err( nf90_inq_varid(ncid,"he3dp_th_tot", id_total_p_he3d_th),flag) 
     IF(flag  == -1)THEN
          fus_prod%total_p_he3d_th = zeroc
     ELSE
 
          CALL netcdf_err( nf90_get_var(ncid,id_total_p_he3d_th,fus_prod%total_p_he3d_th),id_total_p_he3d_th)
     ENDIF
 

     flag = 0
     CALL netcdf_err( nf90_inq_varid(ncid,"he3dp_beam_th", id_p_he3d_beam_th),flag) 
 
     IF(flag == -1)THEN

        fus_prod%p_he3d_beam_thermal = zero_Vector(nj)
     ELSE

        CALL netcdf_err( nf90_get_var(ncid,id_p_he3d_beam_th,work_nj),id_p_he3d_beam_th)
        fus_prod%p_he3d_beam_thermal = new_Vector(nj,work_nj)
     ENDIF
     flag = 0

     CALL netcdf_err( nf90_inq_varid(ncid,"he3dp_beam_th_tot", id_total_p_he3d_beam_th),flag) 
     IF(flag  == -1)THEN
          fus_prod%total_p_he3d_beam_thermal = zeroc
     ELSE

          CALL netcdf_err( nf90_get_var(ncid,id_total_p_he3d_beam_th,fus_prod%total_p_he3d_beam_thermal),id_total_p_he3d_beam_th)
     ENDIF


     flag = 0
     CALL netcdf_err( nf90_inq_varid(ncid,"he3dp_tot", id_total_p_he3d),flag) 


     IF(flag  == -1)THEN
          fus_prod%total_p_he3d = zeroc
     ELSE
          CALL netcdf_err( nf90_get_var(ncid,id_total_p_he3d,fus_prod%total_p_he3d),id_total_p_he3d)
     ENDIF


     flag = 0
     CALL netcdf_err( nf90_inq_varid(ncid,"he3_frac", id_he3_frac),flag) 


     IF(flag  == -1)THEN
          he3_frac = zeroc
     ELSE
          CALL netcdf_err( nf90_get_var(ncid,id_he3_frac,he3_frac),id_he3_frac)
     ENDIF

     flag = 0
     CALL netcdf_err( nf90_inq_varid(ncid,"he3_thermal_spin_pol", id_he3_thermal_spin_pol),flag) 


     IF(flag  == -1)THEN
          he3_thermal_spin_pol = zeroc
     ELSE
          CALL netcdf_err( nf90_get_var(ncid,id_he3_thermal_spin_pol,he3_thermal_spin_pol),id_he3_thermal_spin_pol)
     ENDIF

     flag = 0
     CALL netcdf_err( nf90_inq_varid(ncid,"d_beam_spin_pol", id_d_beam_spin_pol),flag) 
     IF(flag  == -1)THEN
          neut_beam%d_beam_spin_pol = 1.0
     ELSE
          CALL netcdf_err( nf90_get_var(ncid,id_d_beam_spin_pol,neut_beam%d_beam_spin_pol),id_d_beam_spin_pol)
     ENDIF






     IF(.NOT. ALLOCATED(id_plas_freq)) CALL allocate_plasma_freq  

     k=0
     DO i=1,nprim
        k = k+1

        IF(namep(i) == 'dt')THEN
           omega_pi_name = 'omega_pi_D'
           errflg = izero
           CALL netcdf_err( nf90_inq_varid(ncid,omega_pi_name,id_plas_freq(k)),errflg)
           IF(errflg == -1)THEN
              work_nj(:) = zeroc
           ELSE
              CALL netcdf_err( nf90_get_var(ncid,id_plas_freq(k),work_nj),id_plas_freq(k))
           ENDIF
           plasma_frequencies%omega_pi(k) = new_vector(nj,work_nj)

           omega_ci_name = 'omega_ci_D'
           errflg = izero
           CALL netcdf_err( nf90_inq_varid(ncid,omega_ci_name,id_ci_freq(k)),errflg)
           IF(errflg == -1)THEN
              work_nj(:) = zeroc
           ELSE
              CALL netcdf_err( nf90_get_var(ncid,id_ci_freq(k),work_nj),id_ci_freq(k))
           ENDIF
           plasma_frequencies%omega_ci(k) = new_vector(nj,work_nj)

           omega_lh_name = 'omega_lh_D'
           errflg = izero
           CALL netcdf_err( nf90_inq_varid(ncid,omega_lh_name,id_lh_freq(k)),errflg)
           IF(errflg == -1)THEN
              work_nj(:) = zeroc
           ELSE
              CALL netcdf_err( nf90_get_var(ncid,id_lh_freq(k),work_nj),id_lh_freq(k))
           ENDIF
           plasma_frequencies%omega_lh(k) = new_vector(nj,work_nj)

           omega_uh_name = 'omega_uh_D'
           errflg = izero
           CALL netcdf_err( nf90_inq_varid(ncid,omega_uh_name,id_uh_freq(k)),errflg)
           IF(errflg == -1)THEN
              work_nj(:) = zeroc
           ELSE
              CALL netcdf_err( nf90_get_var(ncid,id_uh_freq(k),work_nj),id_uh_freq(k))
           ENDIF
           plasma_frequencies%omega_uh(k) = new_vector(nj,work_nj)

           k =k+1
           omega_pi_name = 'omega_pi_T'
           errflg = izero
           CALL netcdf_err( nf90_inq_varid(ncid,omega_pi_name,id_plas_freq(k)),errflg)
           IF(errflg == -1)THEN
              work_nj(:) = zeroc
           ELSE
              CALL netcdf_err( nf90_get_var(ncid,id_plas_freq(k),work_nj),id_plas_freq(k))
           ENDIF
           plasma_frequencies%omega_pi(k) = new_vector(nj,work_nj)

           omega_ci_name = 'omega_ci_T'
           errflg = izero
           CALL netcdf_err( nf90_inq_varid(ncid,omega_ci_name,id_ci_freq(k)),errflg)
           IF(errflg == -1)THEN
              work_nj(:) = zeroc
           ELSE
              CALL netcdf_err( nf90_get_var(ncid,id_ci_freq(k),work_nj),id_ci_freq(k))
           ENDIF
           plasma_frequencies%omega_ci(k) = new_vector(nj,work_nj)
          
           omega_lh_name = 'omega_lh_T'
           errflg = izero
           CALL netcdf_err( nf90_inq_varid(ncid,omega_lh_name,id_lh_freq(k)),errflg)
           IF(errflg == -1)THEN
              work_nj(:) = zeroc
           ELSE
              CALL netcdf_err( nf90_get_var(ncid,id_lh_freq(k),work_nj),id_lh_freq(k))
           ENDIF
           plasma_frequencies%omega_lh(k) = new_vector(nj,work_nj)
         

           omega_uh_name = 'omega_uh_T'
           errflg = izero
           CALL netcdf_err( nf90_inq_varid(ncid,omega_uh_name,id_uh_freq(k)),id_uh_freq(k))
           IF(errflg == -1)THEN
              work_nj(:) = zeroc
           ELSE
              CALL netcdf_err( nf90_get_var(ncid,id_uh_freq(k),work_nj),id_uh_freq(k))
           ENDIF
           plasma_frequencies%omega_uh(k) = new_vector(nj,work_nj)
          
        ELSE

           omega_pi_name = 'omega_pi_'//namep(i)
           errflg = izero
           CALL netcdf_err( nf90_inq_varid(ncid,omega_pi_name,id_plas_freq(k)),errflg)
           IF(errflg == -1)THEN
              work_nj(:) = zeroc
           ELSE
              CALL netcdf_err( nf90_get_var(ncid,id_plas_freq(k),work_nj),id_plas_freq(k))
           ENDIF
           plasma_frequencies%omega_pi(k) = new_vector(nj,work_nj)

           omega_ci_name = 'omega_ci_'//namep(i)
           errflg = izero
           CALL netcdf_err( nf90_inq_varid(ncid,omega_ci_name,id_ci_freq(k)),errflg)
           IF(errflg == -1)THEN
              work_nj(:) = zeroc
           ELSE
              CALL netcdf_err( nf90_get_var(ncid,id_ci_freq(k),work_nj),id_ci_freq(k))
           ENDIF
           plasma_frequencies%omega_ci(k) = new_vector(nj,work_nj)

           omega_lh_name = 'omega_lh_'//namep(i)
           errflg = izero
           CALL netcdf_err( nf90_inq_varid(ncid,omega_lh_name,id_lh_freq(k)),errflg)
           IF(errflg == -1)THEN
              work_nj(:) = zeroc
           ELSE
              CALL netcdf_err( nf90_get_var(ncid,id_lh_freq(k),work_nj),id_lh_freq(k))
           ENDIF
           plasma_frequencies%omega_lh(k) = new_vector(nj,work_nj)

           omega_uh_name = 'omega_uh_'//namep(i)
           errflg = izero
           CALL netcdf_err( nf90_inq_varid(ncid,omega_uh_name,id_uh_freq(k)),errflg)
           IF(errflg == -1)THEN
              work_nj(:) = zeroc
           ELSE
              CALL netcdf_err( nf90_get_var(ncid,id_uh_freq(k),work_nj),id_uh_freq(k))
           ENDIF
           plasma_frequencies%omega_uh(k) = new_vector(nj,work_nj)

        ENDIF
     ENDDO

     omega_ce_name = 'omega_ce'
     errflg = izero
     CALL netcdf_err( nf90_inq_varid(ncid,omega_ce_name,id_ce_freq),errflg)
     IF(errflg == -1)THEN
        work_nj(:) = zeroc
     ELSE
        CALL netcdf_err( nf90_get_var(ncid,id_ce_freq,work_nj),id_ce_freq)
     ENDIF
     plasma_frequencies%omega_ce = new_vector(nj,work_nj)

     omega_pe_name = 'omega_pe'
     errflg = izero
     CALL netcdf_err( nf90_inq_varid(ncid,omega_pe_name,id_pe_freq),errflg)
     IF(errflg == -1)THEN
        work_nj(:) = zeroc
     ELSE
        CALL netcdf_err( nf90_get_var(ncid,id_pe_freq,work_nj),id_pe_freq)

     ENDIF

     plasma_frequencies%omega_pe = new_vector(nj,work_nj)
     IF(k .ne. nprimp1)THEN
        k = k+1
        plasma_frequencies%omega_pi(k) = zero_vector(nj)
        plasma_frequencies%omega_ci(k) = zero_vector(nj)
        plasma_frequencies%omega_lh(k) = zero_vector(nj)
        plasma_frequencies%omega_uh(k) = zero_vector(nj)
     ENDIF




     ! -------------------------------------------------------------------------
     ! ---  END NETCDF file read/ write processing 
     ! -------------------------------------------------------------------------
  END IF  rwncd






  eqtime = time            !equilibirum time. It is assumed that the quatities
  !fcap,gcap,hcap,eps,xhm2,xi11,xi33,xips
  !were calcualted at this time. Hence subrotuine 
  !neointrp does nothing








2000 CONTINUE





  CALL netcdf_err(nf90_close(ncid) )         !close netcdf file
  
  IF(ALLOCATED(work_nj))DEALLOCATE(work_nj)
  IF(ALLOCATED(work_npsi))DEALLOCATE(work_npsi)
  IF(ALLOCATED(work_nplasbdry))DEALLOCATE(work_nplasbdry)
  IF(ALLOCATED(work_nlimiter))DEALLOCATE(work_nlimiter)
  IF(ALLOCATED(work_nr))DEALLOCATE(work_nr)
  IF(ALLOCATED(work_nz))DEALLOCATE(work_nz)
  IF(ALLOCATED(work_nj_nion))DEALLOCATE(work_nj_nion)
  IF(ALLOCATED(work_nj_nionP1))DEALLOCATE(work_nj_nionp1)

  WRITE(nlog,33)time
  WRITE(ncrt,33)time
33 FORMAT("Done with netcdf read/write at time ",1pe14.6)

  RETURN

END SUBROUTINE iter_dbase_nc


SUBROUTINE netcdf_err(status,flag)

  USE nrtype,                       ONLY : I4B,DP

  USE typeSizes                     !part of netcdf 
  USE netcdf                        !part of netcdf 

  USE iterdbmd_gcnmp,               ONLY : ncid
  USE error_handler,                ONLY : lerrno,terminate

  USE io_gcnmp,                      ONLY : ncrt,nlog

  IMPLICIT NONE
  INTEGER, INTENT ( in) :: status
  INTEGER(I4B)  action,vartyp,nvdims,nvatts,rcode
  INTEGER,PARAMETER                       :: MAXNCNAM = 128 !these are supposed to be
  INTEGER,PARAMETER                       :: MAXVDIMS = 32  !defined in netcdf  but are not ??? 
  INTEGER(I4B) ,DIMENSION(MAXVDIMS)       :: vdims
  INTEGER(I4B), OPTIONAL, INTENT(inout)   :: flag           ! dont call by value (eq negative number)
  CHARACTER(LEN=MAXNCNAM)                 :: varname
  action = 10000
  IF(PRESENT(flag))THEN
        action = flag
  ENDIF
  !print *,'action =',action ! 88888899999
  !print *,'status, nf90_noerr  =',status, nf90_noerr ! 88888899999
  IF(status /= nf90_noerr .AND. action .NE.  0 ) THEN
     WRITE(ncrt,FMT='("  netdcf erro = ",i5)')status
     WRITE(nlog,FMT='("  netdcf erro = ",i5)')status
     WRITE(ncrt,1)TRIM(nf90_strerror(status))
     WRITE(ncrt,FMT='(a)')trim(nf90_strerror(status))
     WRITE(nlog,1)TRIM(nf90_strerror(status)) 
     WRITE(nlog,FMT='(a)')trim(nf90_strerror(status))
!     IF(action > 0)THEN
        CALL ncvinq(ncid,action,varname,vartyp,nvdims,vdims,nvatts,rcode)
        WRITE(ncrt,2)action,varname(1:LEN_TRIM(varname))
        WRITE(nlog,2)action,varname(1:LEN_TRIM(varname))
!     ELSE IF(action < 0)THEN
!        action = -action
!        CALL ncvinq(ncid,action,varname,vartyp,nvdims,vdims,nvatts,rcode)
!        WRITE(ncrt,3)action,varname(1:LEN_TRIM(varname))
!        WRITE(ncrt,FMT='(a)')trim(nf90_strerror(status))
!        WRITE(nlog,3)action,varname(1:LEN_TRIM(varname))
!        WRITE(nlog,FMT='(a)')trim(nf90_strerror(status))
!     ENDIF
     lerrno = 32
     CALL  terminate(lerrno,nlog)
  ELSE IF(status /= nf90_noerr .AND. action ==  0 ) THEN
      flag =  -1
  END IF
  RETURN

!3    FORMAT(2x,'READ  error caused by netcdf variable no :',i5,' name: ',a)
2    FORMAT(2x,'WRITE or READ  error caused by netcdf variable no :',i5,' name: ',a)
1    FORMAT(2x,a)

END SUBROUTINE netcdf_err



SUBROUTINE set_label(label,base_label,strlen,indicator)
  !  --------------------------------------------------------------------
  USE nrtype,              ONLY : I4B
  USE ions_gcnmp,          ONLY : nprim,nimp,nion,namei,  &
                                  name_size,namep,namen,nneu
  USE neutral_beams,       ONLY:  nameb,nbion
  USE io_gcnmp,            ONLY : ncrt,nlog

  USE error_handler,       ONLY : lerrno,terminate

  IMPLICIT NONE
  INTEGER(I4B) jp,ji,rlab,blen,jj,strlen
  CHARACTER(LEN = strlen )label,base_label
  CHARACTER(LEN=*) indicator


  jp = 0
  ji = 0
  blen = LEN_TRIM(base_label)
  label(1:blen) = base_label(1:blen)
  label(blen+1:strlen) = ' '
  rlab = blen
  IF(indicator  == 'nion' .OR. indicator == 'ntot' .OR. indicator == 'bremssthralung') THEN
     DO jj=1,nion
        IF (jj .LE. nprim)  jp = jp + 1
        IF (jj .GT. nprim)  ji = ji + 1
        IF ( jj .LE. nprim .AND. rlab .LT. strlen - name_size- 1 )THEN
           label = label(1:rlab)//" "//namep(jp)
           rlab = LEN_TRIM(label)
        ELSE IF ( jj .LE. nprim .AND. rlab .GE. strlen - name_size -1 )THEN
           label = "Warning:  truncated label nion due to overflow"
           WRITE(ncrt,FMT='(a)')label
           WRITE(nlog,FMT='(a)')label
        ELSE IF (jj .GT. nprim .AND. rlab .LT. strlen - name_size -1 )THEN
           label = label(1:rlab)//" "//namei(ji)
           rlab = LEN_TRIM(label)
        ELSE IF (jj .GT. nprim .AND. rlab .GE. strlen - name_size -1 )THEN
           label = "Warning:  truncated label nion due to overflow"
           WRITE(ncrt,FMT='(a)')label
           WRITE(nlog,FMT='(a)')label
        ENDIF
     ENDDO
  ELSE IF (indicator == 'nf' ) THEN
     DO jj=1,nbion
        IF ( rlab .LT. strlen - name_size- 1 )THEN
           label = label(1:rlab)//" "//nameb(jj)
           rlab = LEN_TRIM(label)
        ELSE IF ( rlab .GE. strlen - name_size -1 )THEN
           label = "Warning:  truncated label nbion  due to overflow"
           WRITE(ncrt,FMT='(a)')label
           WRITE(nlog,FMT='(a)')label
        ENDIF
     ENDDO

  ELSE IF (indicator == 'nprim' ) THEN
     DO jj=1,nprim
        IF ( rlab .LT. strlen - name_size- 1 )THEN
           label = label(1:rlab)//" "//namep(jj)
           rlab = LEN_TRIM(label)
        ELSE IF ( rlab .GE. strlen - name_size -1 )THEN
           label = "Warning:  truncated label  namep  due to overflow"
           WRITE(ncrt,FMT='(a)')label
           WRITE(nlog,FMT='(a)')label
        ENDIF
     ENDDO
  ELSE IF (indicator == 'nneu' ) THEN
     DO jj=1,nneu
        IF ( rlab .LT. strlen - name_size- 1 )THEN
           label = label(1:rlab)//" "//namen(jj)
           rlab = LEN_TRIM(label)
        ELSE IF ( rlab .GE. strlen - name_size -1 )THEN
           label = "Warning:  truncated label  namen  due to overflow"
           WRITE(ncrt,FMT='(a)')label
           WRITE(nlog,FMT='(a)')label
        ENDIF
     ENDDO
  ELSE IF (indicator == 'electron' ) THEN
        IF ( rlab .LT. strlen - name_size- 1 )THEN
           label = label(1:rlab)//" "//'elct'
           rlab = LEN_TRIM(label)
        ELSE IF ( rlab .GE. strlen - name_size -1 )THEN
           label = "Warning:  truncated label electron  due to overflow"
           WRITE(ncrt,FMT='(a)')label
           WRITE(nlog,FMT='(a)')label
        ENDIF
  ELSE IF (indicator == 'glf-etg' )THEN
        IF ( rlab .LT. strlen - name_size- 1 )THEN
           rlab = LEN_TRIM(label)
        ELSE IF ( rlab .GE. strlen - name_size -1 )THEN
           label = "Warning:  truncated label glf etg  due to overflow"
           WRITE(ncrt,FMT='(a)')label
           WRITE(nlog,FMT='(a)')label
        ENDIF
  ELSE IF (indicator == 'glf_gam_net_i' .OR. indicator == 'glf_gam_net_e' )THEN
        IF ( rlab .LT. strlen - name_size- 1 )THEN
           rlab = LEN_TRIM(label)
        ELSE IF ( rlab .GE. strlen - name_size -1 )THEN
           label = "Warning:  truncated label glf net growth rate   due to overflow"
           WRITE(ncrt,FMT='(a)')label
           WRITE(nlog,FMT='(a)')label
        ENDIF

  ELSE IF (indicator == 'glf_anrate' .OR. indicator == 'glf_anrate2' )THEN
        IF ( rlab .LT. strlen - name_size- 1 )THEN
           rlab = LEN_TRIM(label)
        ELSE IF ( rlab .GE. strlen - name_size -1 )THEN
           label = "Warning:  truncated label glf leading mode growth rate   due to overflow"
           WRITE(ncrt,FMT='(a)')label
           WRITE(nlog,FMT='(a)')label
        ENDIF
  ELSE IF (indicator == 'glf_anfreq' .OR. indicator == 'glf_anfreq2' )THEN
        IF ( rlab .LT. strlen - name_size- 1 )THEN
           rlab = LEN_TRIM(label)
        ELSE IF ( rlab .GE. strlen - name_size -1 )THEN
           label = "Warning:  truncated label glf leading frequency due to overflow"
           WRITE(ncrt,FMT='(a)')label
           WRITE(nlog,FMT='(a)')label
        ENDIF
  ELSE IF (indicator == 'effective prim ion' ) THEN
        IF ( rlab .LT. strlen - name_size- 1 )THEN
           label = label(1:rlab)//" "//'eff prim ion'
           rlab = LEN_TRIM(label)
        ELSE IF ( rlab .GE. strlen - name_size -1 )THEN
           label = "Warning:  truncated label effective prim ion  due to overflow"
           WRITE(ncrt,FMT='(a)')label
           WRITE(nlog,FMT='(a)')label
        ENDIF
  ELSE IF (indicator == 'effective impurity' ) THEN
        IF ( rlab .LT. strlen - name_size- 1 )THEN
           label = label(1:rlab)//" "//'eff impurity'
           rlab = LEN_TRIM(label)
        ELSE IF ( rlab .GE. strlen - name_size -1 )THEN
           label = "Warning:  truncated label effective impurity due to overflow"
           WRITE(ncrt,FMT='(a)')label
           WRITE(nlog,FMT='(a)')label
        ENDIF
  ELSE IF (indicator == 'sb')THEN
       IF ( rlab .LT. strlen - name_size- 1 )THEN
           label = label(1:rlab)//" "//'sb #/m3 sec'
           rlab = LEN_TRIM(label)
        ELSE IF ( rlab .GE. strlen - name_size -1 )THEN
           label = "Warning:  truncated label for neutral beam sb"
           WRITE(ncrt,FMT='(a)')label
           WRITE(nlog,FMT='(a)')label
        ENDIF
  ELSE IF (indicator == 'qb')THEN
       IF ( rlab .LT. strlen - name_size- 1 )THEN
           label = label(1:rlab)//" "//'qb watts/m3'
           rlab = LEN_TRIM(label)
        ELSE IF ( rlab .GE. strlen - name_size -1 )THEN
           label = "Warning:  truncated label for neutral beam qb"
           WRITE(ncrt,FMT='(a)')label
           WRITE(nlog,FMT='(a)')label
        ENDIF
  ELSE IF (indicator == 'spb')THEN
       IF ( rlab .LT. strlen - name_size- 1 )THEN
           label = label(1:rlab)//" "//'spb kg/(m2-s2)'
           rlab = LEN_TRIM(label)
        ELSE IF ( rlab .GE. strlen - name_size -1 )THEN
           label = "Warning:  truncated label for neutral beam spb"
           WRITE(ncrt,FMT='(a)')label
           WRITE(nlog,FMT='(a)')label
        ENDIF
  ELSE IF (indicator == 'spbr')THEN
       IF ( rlab .LT. strlen - name_size- 1 )THEN
           label = label(1:rlab)//" "//'spbr kg/(m-s2)'
           rlab = LEN_TRIM(label)
        ELSE IF ( rlab .GE. strlen - name_size -1 )THEN
           label = "Warning:  truncated label for neutral beam spbr"
           WRITE(ncrt,FMT='(a)')label
           WRITE(nlog,FMT='(a)')label
        ENDIF
  ELSE IF (indicator == 'angmpf')THEN
       IF ( rlab .LT. strlen - name_size- 1 )THEN
           label = label(1:rlab)//" "//'angmpf kg m2/s)'
           rlab = LEN_TRIM(label)
        ELSE IF ( rlab .GE. strlen - name_size -1 )THEN
           label = "Warning:  truncated label for neutral beam angmpf"
           WRITE(ncrt,FMT='(a)')label
           WRITE(nlog,FMT='(a)')label
        ENDIF
  ELSE IF (indicator == 'pb0')THEN
       IF ( rlab .LT. strlen - name_size- 1 )THEN
           label = label(1:rlab)//" "//'pb0 kg m/s)'
           rlab = LEN_TRIM(label)
        ELSE IF ( rlab .GE. strlen - name_size -1 )THEN
           label = "Warning:  truncated label for neutral beam pb0"
           WRITE(ncrt,FMT='(a)')label
           WRITE(nlog,FMT='(a)')label
        ENDIF
  ELSE IF (indicator == 'hibr')THEN
       IF ( rlab .LT. strlen - name_size- 1 )THEN
           label = label(1:rlab)//" "//'hibr)'
           rlab = LEN_TRIM(label)
        ELSE IF ( rlab .GE. strlen - name_size -1 )THEN
           label = "Warning:  truncated label for neutral beam hibr"
           WRITE(ncrt,FMT='(a)')label
           WRITE(nlog,FMT='(a)')label
        ENDIF
  ELSE IF (indicator == 'hdep')THEN
       IF ( rlab .LT. strlen - name_size- 1 )THEN
           label = label(1:rlab)//" "//'hdep)'
           rlab = LEN_TRIM(label)
        ELSE IF ( rlab .GE. strlen - name_size -1 )THEN
           label = "Warning:  truncated label for neutral beam hdep"
           WRITE(ncrt,FMT='(a)')label
           WRITE(nlog,FMT='(a)')label
        ENDIF
 ELSE IF (indicator == 'zeta')THEN
       IF ( rlab .LT. strlen - name_size- 1 )THEN
           label = label(1:rlab)//" "//'zeta)'
           rlab = LEN_TRIM(label)
        ELSE IF ( rlab .GE. strlen - name_size -1 )THEN
           label = "Warning:  truncated label for neutral beam zeta"
           WRITE(ncrt,FMT='(a)')label
           WRITE(nlog,FMT='(a)')label
        ENDIF
  ELSE IF (indicator == 'ebeam')THEN
        IF ( rlab .LT. strlen - name_size- 1 )THEN
           label = label(1:rlab)//" "//'ebeam Kev)'
           rlab = LEN_TRIM(label)
        ELSE IF ( rlab .GE. strlen - name_size -1 )THEN
           label = "Warning:  truncated label for neutral beam ebeam"
           WRITE(ncrt,FMT='(a)')label
           WRITE(nlog,FMT='(a)')label
        ENDIF
  ELSE IF (indicator == 'pbeam')THEN
        IF ( rlab .LT. strlen - name_size- 1 )THEN
           label = label(1:rlab)//" "//'pbeam Watts)'
           rlab = LEN_TRIM(label)
        ELSE IF ( rlab .GE. strlen - name_size -1 )THEN
           label = "Warning:  truncated label for neutral beam pbeam"
           WRITE(ncrt,FMT='(a)')label
           WRITE(nlog,FMT='(a)')label
        ENDIF
  ELSE IF (indicator == 'bptor')THEN
        IF ( rlab .LT. strlen - name_size- 1 )THEN
           label = label(1:rlab)//" "//'bptor Watts)'
           rlab = LEN_TRIM(label)
        ELSE IF ( rlab .GE. strlen - name_size -1 )THEN
           label = "Warning:  truncated label for neutral beam bptor"
           WRITE(ncrt,FMT='(a)')label
           WRITE(nlog,FMT='(a)')label
        ENDIF
 ELSE IF (indicator == 'fbcur')THEN
        IF ( rlab .LT. strlen - name_size- 1 )THEN
           label = label(1:rlab)//" "//'fbcur)'
           rlab = LEN_TRIM(label)
        ELSE IF ( rlab .GE. strlen - name_size -1 )THEN
           label = "Warning:  truncated label for neutral beam fbcur"
           WRITE(ncrt,FMT='(a)')label
           WRITE(nlog,FMT='(a)')label
        ENDIF
 ELSE IF (indicator == 'prompt_nb_pwr')THEN 
        IF ( rlab .LT. strlen - name_size- 1 )THEN
           label = label(1:rlab)//" "//'prompt_nb_pwr)'
           rlab = LEN_TRIM(label)
        ELSE IF ( rlab .GE. strlen - name_size -1 )THEN
           label = "Warning:  truncated label for neutral beam prompt_pwr_in_plasma"
           WRITE(ncrt,FMT='(a)')label
           WRITE(nlog,FMT='(a)')label
        ENDIF
 ELSE IF (indicator == 'bneut')THEN
        IF ( rlab .LT. strlen - name_size- 1 )THEN
           label = label(1:rlab)//" "//'neutral intensity)'
           rlab = LEN_TRIM(label)
        ELSE IF ( rlab .GE. strlen - name_size -1 )THEN
           label = "Warning:  truncated label for neutral beam bneut"
           WRITE(ncrt,FMT='(a)')label
           WRITE(nlog,FMT='(a)')label
        ENDIF
 ELSE IF (indicator == 'bion')THEN
        IF ( rlab .LT. strlen - name_size- 1 )THEN
           label = label(1:rlab)//" "//'neutral beam ion intensity)'
           rlab = LEN_TRIM(label)
        ELSE IF ( rlab .GE. strlen - name_size -1 )THEN
           label = "Warning:  truncated label for neutral beam bion"
           WRITE(ncrt,FMT='(a)')label
           WRITE(nlog,FMT='(a)')label
        ENDIF
  ELSE IF (indicator == 'fap')THEN
        IF ( rlab .LT. strlen - name_size- 1 )THEN
           label = label(1:rlab)//" "//'aperture loss)'
           rlab = LEN_TRIM(label)
        ELSE IF ( rlab .GE. strlen - name_size -1 )THEN
           label = "Warning:  truncated label for neutral beam fap"
           WRITE(ncrt,FMT='(a)')label
           WRITE(nlog,FMT='(a)')label
        ENDIF
  ELSE IF (indicator == 'fwall')THEN
        IF ( rlab .LT. strlen - name_size- 1 )THEN
           label = label(1:rlab)//" "//'wall loss)'
           rlab = LEN_TRIM(label)
        ELSE IF ( rlab .GE. strlen - name_size -1 )THEN
           label = "Warning:  truncated label for neutral beam fwall"
           WRITE(ncrt,FMT='(a)')label
           WRITE(nlog,FMT='(a)')label
        ENDIF

  ELSE IF (indicator == 'forb')THEN
        IF ( rlab .LT. strlen - name_size- 1 )THEN
           label = label(1:rlab)//" "//'orbit loss)'
           rlab = LEN_TRIM(label)
        ELSE IF ( rlab .GE. strlen - name_size -1 )THEN
           label = "Warning:  truncated label for neutral beam forb"
           WRITE(ncrt,FMT='(a)')label
           WRITE(nlog,FMT='(a)')label
        ENDIF
 ELSE IF (indicator == 'fber')THEN
        IF ( rlab .LT. strlen - name_size- 1 )THEN
           label = label(1:rlab)//" "//'error in trapped ion calcs)'
           rlab = LEN_TRIM(label)
        ELSE IF ( rlab .GE. strlen - name_size -1 )THEN
           label = "Warning:  truncated label for neutral beam fber"
           WRITE(ncrt,FMT='(a)')label
           WRITE(nlog,FMT='(a)')label
        ENDIF
 ELSE IF (indicator == 'fb00')THEN
        IF ( rlab .LT. strlen - name_size- 1 )THEN
           label = label(1:rlab)//" "//')'
           rlab = LEN_TRIM(label)
        ELSE IF ( rlab .GE. strlen - name_size -1 )THEN
           label = "Warning:  truncated label for neutral beam fb00"
           WRITE(ncrt,FMT='(a)')label
           WRITE(nlog,FMT='(a)')label
        ENDIF
 ELSE IF (indicator == 'fb01')THEN
        IF ( rlab .LT. strlen - name_size- 1 )THEN
           label = label(1:rlab)//" "//')'
           rlab = LEN_TRIM(label)
        ELSE IF ( rlab .GE. strlen - name_size -1 )THEN
           label = "Warning:  truncated label for neutral beam fb01"
           WRITE(ncrt,FMT='(a)')label
           WRITE(nlog,FMT='(a)')label
        ENDIF
 ELSE IF (indicator == 'fb10')THEN
        IF ( rlab .LT. strlen - name_size- 1 )THEN
           label = label(1:rlab)//" "//')'
           rlab = LEN_TRIM(label)
        ELSE IF ( rlab .GE. strlen - name_size -1 )THEN
           label = "Warning:  truncated label for neutral beam fb10"
           WRITE(ncrt,FMT='(a)')label
           WRITE(nlog,FMT='(a)')label
        ENDIF
 ELSE IF (indicator == 'fb11')THEN
        IF ( rlab .LT. strlen - name_size- 1 )THEN
           label = label(1:rlab)//" "//')'
           rlab = LEN_TRIM(label)
        ELSE IF ( rlab .GE. strlen - name_size -1 )THEN
           label = "Warning:  truncated label for neutral beam fb11"
           WRITE(ncrt,FMT='(a)')label
           WRITE(nlog,FMT='(a)')label
        ENDIF
 ELSE IF (indicator == 'wb00')THEN
        IF ( rlab .LT. strlen - name_size- 1 )THEN
           label = label(1:rlab)//" "//')'
           rlab = LEN_TRIM(label)
        ELSE IF ( rlab .GE. strlen - name_size -1 )THEN
           label = "Warning:  truncated label for neutral beam wb00"
           WRITE(ncrt,FMT='(a)')label
           WRITE(nlog,FMT='(a)')label
        ENDIF
 ELSE IF (indicator == 'wb01')THEN
        IF ( rlab .LT. strlen - name_size- 1 )THEN
           label = label(1:rlab)//" "//')'
           rlab = LEN_TRIM(label)
        ELSE IF ( rlab .GE. strlen - name_size -1 )THEN
           label = "Warning:  truncated label for neutral beam wb01"
           WRITE(ncrt,FMT='(a)')label
           WRITE(nlog,FMT='(a)')label
        ENDIF
 ELSE IF (indicator == 'wb10')THEN
        IF ( rlab .LT. strlen - name_size- 1 )THEN
           label = label(1:rlab)//" "//')'
           rlab = LEN_TRIM(label)
        ELSE IF ( rlab .GE. strlen - name_size -1 )THEN
           label = "Warning:  truncated label for neutral beam wb10"
           WRITE(ncrt,FMT='(a)')label
           WRITE(nlog,FMT='(a)')label
        ENDIF
 ELSE IF (indicator == 'wb11')THEN
        IF ( rlab .LT. strlen - name_size- 1 )THEN
           label = label(1:rlab)//" "//')'
           rlab = LEN_TRIM(label)
        ELSE IF ( rlab .GE. strlen - name_size -1 )THEN
           label = "Warning:  truncated label for neutral beam wb11"
           WRITE(ncrt,FMT='(a)')label
           WRITE(nlog,FMT='(a)')label
        ENDIF

  ELSE IF (indicator == 'hicme')THEN
        IF ( rlab .LT. strlen - name_size- 1 )THEN
           label = label(1:rlab)//" "//'hicm electron creation mode)'
           rlab = LEN_TRIM(label)
        ELSE IF ( rlab .GE. strlen - name_size -1 )THEN
           label = "Warning:  truncated label for hicme"
           WRITE(ncrt,FMT='(a)')label
           WRITE(nlog,FMT='(a)')label
        ENDIF
  ELSE IF (indicator == 'hicmp1')THEN
        IF ( rlab .LT. strlen - name_size- 1 )THEN
           label = label(1:rlab)//" "//'hicm primary ion 1  creation mode)'
           rlab = LEN_TRIM(label)
        ELSE IF ( rlab .GE. strlen - name_size -1 )THEN
           label = "Warning:  truncated label for hicmp1"
           WRITE(ncrt,FMT='(a)')label
           WRITE(nlog,FMT='(a)')label
        ENDIF
  ELSE IF (indicator == 'hicmp2')THEN
        IF ( rlab .LT. strlen - name_size- 1 )THEN
           label = label(1:rlab)//" "//'hicm primary ion 2  creation mode)'
           rlab = LEN_TRIM(label)
        ELSE IF ( rlab .GE. strlen - name_size -1 )THEN
           label = "Warning:  truncated label for hicmp2"
           WRITE(ncrt,FMT='(a)')label
           WRITE(nlog,FMT='(a)')label
        ENDIF
  ELSE IF (indicator == 'vpinch Nclass')THEN
        IF ( rlab .LT. strlen - name_size- 1 )THEN
           label = label(1:rlab)//" "//'e,plus ion species'
           rlab = LEN_TRIM(label)
        ELSE IF ( rlab .GE. strlen - name_size -1 )THEN
           label = "Warning:  truncated label for vpinch Nclass"
           WRITE(ncrt,FMT='(a)')label
           WRITE(nlog,FMT='(a)')label
        ENDIF
  ELSE IF (indicator == 'mmm_gammaDBM')THEN
    RETURN
  ELSE IF (indicator == 'mmm_omegaDBM')THEN
    RETURN
  ELSE IF (indicator == 'mmm_xdi')THEN
    RETURN
  ELSE IF (indicator == 'mmm_xte')THEN  
    RETURN
  ELSE IF (indicator == 'mmm_xti')THEN 
    RETURN
  ELSE IF (indicator == 'mmm_xdz')THEN 
    RETURN
  ELSE IF (indicator == 'mmm_xvt')THEN 
    RETURN
  ELSE IF (indicator == 'mmm_xvp')THEN 
    RETURN
  ELSE IF (indicator == 'mmm_xtiW20')THEN 
    RETURN
  ELSE IF (indicator == 'mmm_xteW20')THEN 
    RETURN
  ELSE IF (indicator == 'mmm_xdiW20')THEN 
    RETURN
  ELSE IF (indicator == 'mmm_xtiDBM')THEN 
    RETURN
  ELSE IF (indicator == 'mmm_xteDBM')THEN  
    RETURN 
  ELSE IF (indicator == 'mmm_xdiDBM')THEN  
    RETURN  
  ELSE IF (indicator == 'mmm_xteETG')THEN  
    RETURN  
  ELSE IF (indicator == 'mmm_gamma_i1_W20')THEN 
    RETURN
  ELSE IF (indicator == 'mmm_gamma_i2_W20')THEN 
    RETURN
  ELSE IF (indicator == 'mmm_gamma_e1_W20')THEN 
    RETURN
  ELSE IF (indicator == 'mmm_gamma_e2_W20')THEN    
    RETURN
  ELSE IF (indicator == 'mmm_omega_i1_W20')THEN 
    RETURN
  ELSE IF (indicator == 'mmm_omega_i2_W20')THEN
    RETURN 
  ELSE IF (indicator == 'mmm_omega_e1_W20')THEN 
    RETURN
  ELSE IF (indicator == 'mmm_omega_e2_W20')THEN
    RETURN    
  ELSE IF (indicator == 'mmm_flux_ith')THEN
    RETURN
  ELSE IF (indicator == 'mmm_flux_eth')THEN
    RETURN
  ELSE IF (indicator == 'mmm_flux_ip')THEN
    RETURN
  ELSE IF (indicator == 'mmm_flux_imp')THEN
    RETURN
  ELSE IF (indicator == 'mmm_vconv_ith')THEN
    RETURN
  ELSE IF (indicator == 'mmm_vconv_eth')THEN
    RETURN
  ELSE IF (indicator == 'mmm_vconv_ip')THEN
    RETURN
  ELSE IF (indicator == 'mmm_vconv_imp')THEN
    RETURN
  ELSE IF (indicator == 'mmm_vmtmt')THEN
    RETURN
  ELSE IF (indicator == 'mmm_vmtmp')THEN
    RETURN
  ELSE IF (indicator == 'rhog_beam')THEN
    RETURN

  ELSE IF (indicator == 'omega_pi_D' .OR. indicator == 'omega_pi_T' )THEN
        IF ( rlab .LT. strlen - name_size - 13 )THEN
           label = label(1:rlab)//" "//indicator
           rlab = LEN_TRIM(label)
        ELSE IF ( rlab .GE. strlen - name_size -1 )THEN
           label = "Warning:  truncated label for vpinch Nclass"
           WRITE(ncrt,FMT='(a)')label
           WRITE(nlog,FMT='(a)')label
        ENDIF

  ELSE
     label = "ERROR, indicator = "//indicator
     WRITE(ncrt,FMT='(a)')label
     WRITE(nlog,FMT='(a)')label
     lerrno =35
     CALL  terminate(lerrno,nlog)
  ENDIF
  RETURN
END SUBROUTINE set_label




SUBROUTINE  search_dims(gen_name,dim_names,nf_max_name,dim_size ,max_dims,index)
  USE nrtype,               ONLY : I4B
#ifdef GCNMP
#define USESUBS
#endif
#ifdef NFREYA
#define USESUBS
#endif

#ifdef USESUBS
  USE io_gcnmp,             ONLY : ncrt,nlog
#else
  USE io,                   ONLY : ncrt,nlog => nitre
#endif
  IMPLICIT NONE
  INTEGER  dim_size,max_dims,nf_max_name
  INTEGER  index
  CHARACTER(LEN=dim_size) dim_names(max_dims)
  CHARACTER(LEN=nf_max_name) gen_name
  CHARACTER(LEN = 72) message
  INTEGER(I4B) j,lent
  index = 0
  DO j=1, max_dims
     dim_names(j) = ADJUSTL(dim_names(j))
     gen_name  = ADJUSTL(gen_name)
     lent = LEN_TRIM(dim_names(j))
     IF(gen_name(1:lent) == dim_names(j)(1:lent)) THEN
        index =j
        RETURN  
     ENDIF
  ENDDO
  message ='error did not find dimension :'//gen_name(1:lent)
  WRITE(ncrt,FMT='(a)')message
  WRITE(nlog,FMT='(a)')message
  RETURN
END SUBROUTINE search_dims
    
