   MODULE fusion_gcnmp
      USE nrtype,                                  ONLY : DP,I4B
      USE common_constants,                        ONLY : zeroc,joupkev

      USE grid_class,                              ONLY : nj,r,hcap,volfac

      IMPLICIT NONE

      REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:)    ::         &
           ddn_fus,dtn_fus,tt2n_fus,ddp_fus,dhe_fus,         &
           he3dp_fus

      REAL(DP),PUBLIC                             ::         &
           ddn_tot,dtn_tot,tt2n_tot,dhe_tot,ddp_tot,         &
           ddnp_branch_ratio,pfuse_tot,he3_frac,             &
           he3_thermal_spin_pol
      REAL(DP),PUBLIC,ALLOCATABLE,DIMENSION(:)    ::         &
           neutr_ddn_th,neutr_ddn_beam_beam,                 &  ! neutron rates
           neutr_ddn_beam_thermal,                           &
           neutr_ddn_knock,neutr_ddn_tot,qfuse_thermal

      LOGICAL internal_thermal_fusion,differential_burnup

      INTEGER,PARAMETER   ::  alpha_slow_size = 10
      CHARACTER(len = alpha_slow_size) :: alpha_slow

      CONTAINS

      SUBROUTINE allocate_fusion

        IF(.NOT.  ALLOCATED(ddn_fus))THEN
           ALLOCATE(ddn_fus(nj))
           ALLOCATE(dtn_fus(nj))
           ALLOCATE(tt2n_fus(nj))
           ALLOCATE(ddp_fus(nj))
           ALLOCATE(dhe_fus(nj))
           ALLOCATE(he3dp_fus(nj))
           ALLOCATE(neutr_ddn_th(nj))
           ALLOCATE(neutr_ddn_beam_beam(nj))
           ALLOCATE(neutr_ddn_beam_thermal(nj))
           ALLOCATE(neutr_ddn_knock(nj))
           ALLOCATE(neutr_ddn_tot(nj))
           ddnp_branch_ratio = 0.5_DP
           ddn_fus = zeroc ; dtn_fus = zeroc ; tt2n_fus = zeroc
           ddp_fus = zeroc ; dhe_fus = zeroc ; neutr_ddn_th = zeroc
           neutr_ddn_beam_beam = zeroc ;  neutr_ddn_knock= zeroc 
           neutr_ddn_tot = zeroc 
        ENDIF

        RETURN

      END SUBROUTINE allocate_fusion

           

      SUBROUTINE thermonuclear_rate
! ------------------------------------------------- 10/28/07--- HSJ ---
!     BASIC THERMAL FUSION REACTION RATE CALCULATIONS:
!
!     DDFUSN  = D(D,   N)HE3       THERMAL-THERMAL (I.E., THERMONUCLEAR)
!     DDPFUS  = D(D,   P)T           "       "             "
!     DTNFUS  = D(T,   N)HE4         "       "             "
!     TTNFUS  = T(T,  2N)HE4         "       "             "
!     HDPFUS  = HE3(D, P)HE4         "       "             "
!
!             RATE COEFFICIENT AND CROSS SECTIONS are from
!             Bosch & Hale, Nuc. Fus., vol32, no.4 (1992) 611
!             (except for t(t,2n)he4 ,see subroutine tt2nrate )
!
!     INPUT
!
!     ti(j)            j=1,2..nj  ion temperature,kev
!     en(j,i)          j=1,2..nj,i=1,2..nprim, primary ion densities,#/m**3
!
!     fd_thermal       fraction of deuterium in dt mixture
!
!
!     volfac           volume factor =4.0*pi^2*rmajor
!     r(j)             j=1,2,..nj rho grid,m
!     d_index                  index for species deuterim
!     dt_index                                   dt mixture
!     t_index                                 tritium
!     he_index
!     hcap
!
!     OUTPUT:
!      thermal d,d reactants 
!             ddn_fus(j)   local reaction (i.e., neutron production) rate,
!                        #/m**3/sec
!             ddn_tot       volume integrated rate, #/sec
!      thermal d,d reactants 
!             ddp_fus(j)   local reaction (i.e., proton production) rate,
!                        #/m**3/sec
!             ddp_tot       volume integrated rate, #/sec
!      thermal d,t reactants: 
!             dtn_fus(j)   local reaction (i.e., neutron production) rate,
!                        #/m**3/sec
!             dtn_tot       volume integrated rate, #/sec
!      thermal t,t reactants: 
!             ttn_fus(j)   local reaction (i.e., neutron production) rate,
!                        #/m**3/sec
!             ttn_tot       volume integrated rate, #/sec
!      thermal he,d reactants: 
!             dhe_fus(j)   local reaction (i.e., neutron production) rate,
!                        #/m**3/sec
!             dhe_tot       volume integrated rate, #/sec
!
! ----------------------------------------------------------------------
!
!


      USE ions_gcnmp,                              ONLY : d_index,t_index,he_index,    &
                                                          dt_index,fd_thermal



      USE dep_var,                                 ONLY : te,ti,en

      IMPLICIT  NONE

      REAL(DP)   sigmav_ddn,sigmav_ddp,sigmav_dtn,sigmavd_total,sigmavt_total,          &
                 sigmav_tt2n,sigmav_dhe
      REAL(DP)   nd,nt,nhe

      INTEGER(I4B) j 

!
!     zero volume integrated rates:
!
      nhe      = zeroc
      nd       = zeroc
      nt       = zeroc
      ddn_tot  = zeroc
      dtn_tot  = zeroc
      tt2n_tot = zeroc
      dhe_tot  = zeroc
      ddp_tot  = zeroc


!     get fusion rates as function of rho:
      DO j = 1,nj
         sigmav_ddn    = zeroc 
         sigmav_ddp    = zeroc 
         sigmav_dtn    = zeroc 
         sigmav_dhe    = zeroc
         sigmav_tt2n   = zeroc                              ! sigmav in m**3/sec
         IF(d_index .GT. 0)THEN                             ! get d(d,n)he3 and d(d,p)t rate

                !d(d,n)he3 reaction: 
                sigmav_ddn = ddnrate(ti(j)) / 2.0_DP         ! like-like reaction  divide by 2

                !d(d,p)t reaction:
                sigmav_ddp=ddprate(ti(j))/2.0_DP             ! like-like reaction  divide by 2
                     
                nd       = en(j,d_index)
         ENDIF

         IF(t_index .GT. 0)THEN                             ! get t(t,2n) rate
                sigmav_tt2n  = tt2nrate(ti(j))/2.0_DP       ! like-like reaction  divide by 2
                nt          = en(j,t_index)
         ENDIF

         IF(t_index .GT. 0 .AND. d_index .GT. 0)THEN        ! get d(t,n)he4 rate
            sigmav_dtn = dtnrate(ti(j))
         ENDIF



         IF(dt_index .GT. 0)THEN                            ! get all of above from thermal mixture
            sigmav_dtn      = dtnrate(ti(j))
            sigmav_ddn      = ddnrate(ti(j))  / 2.0_DP
            sigmav_tt2n     = tt2nrate(ti(j)) / 2.0_DP       ! like-like reaction  divide by 2
            nt              = (1._DP - fd_thermal)*en(j,dt_index)
            nd              = fd_thermal*en(j,dt_index)
         ENDIF

         IF(he_index .GT. 0 .AND. d_index .GT. 0)THEN        ! get he3(d,p)he4 rate
            sigmav_dhe     = hdprate(ti(j))       
            nhe            = en(j,he_index)                  !assumes he is primary ion
                                                             ! see sub set_ion_prop
         ENDIF

        sigmavd_total  = sigmav_ddn  + sigmav_ddp + sigmav_dtn + sigmav_dhe  ! valid because ti same for all
        IF(sigmavd_total > zeroc)THEN
           ddn_fus(j)     =  nd*nd  * sigmav_ddn    ! sigmav in m**3/sec
           ddp_fus(j)     =  nd*nd  * sigmav_ddp
           dtn_fus(j)     =  nd*nt  * sigmav_dtn
           dhe_fus(j)     =  nd*nhe * sigmav_dhe
        ELSE
           ddn_fus(j)     = zeroc
           ddp_fus(j)     = zeroc
           dtn_fus(j)     = zeroc
           dhe_fus(j)     = zeroc
        ENDIF

        sigmavt_total     = sigmav_tt2n + sigmav_dtn                 ! valid because ti same for all
        IF(sigmavt_total > zeroc)THEN
           tt2n_fus(j)    = nt*nt   * sigmav_tt2n
        ELSE
           tt2n_fus(j) = zeroc
        ENDIF

      ENDDO ! over grid points

      RETURN
      END   SUBROUTINE thermonuclear_rate



      SUBROUTINE integrated_fusion_rates
! --------------------------------------------------------------
! -- do volume integrals of fusion rates
! --------------------------------------------------------------

        CALL trapv(r,ddp_fus,hcap,nj,ddp_tot)
        ddp_tot = ddp_tot*volfac

        CALL trapv(r,ddn_fus,hcap,nj,ddn_tot)
        ddn_tot = ddn_tot*volfac
        
        CALL trapv(r,dtn_fus,hcap,nj,dtn_tot)
        dtn_tot = dtn_tot*volfac
                
        CALL trapv(r,tt2n_fus,hcap,nj,tt2n_tot)
        tt2n_tot = tt2n_tot*volfac
        
        CALL trapv(r,dhe_fus,hcap,nj,dhe_tot)
        dhe_tot = dhe_tot*volfac
        
        RETURN
      END SUBROUTINE integrated_fusion_rates



      REAL(DP) FUNCTION ddnrate (ti)
!
!
! ------------------------------------------------- 10/29/07 --- HSJ ---
!             returns rate(m**3/sec) of d(d,n)he3 reaction
!             (the division by 2 for collisions between like
!              particles is NOT accounted for here).
!             NEW Bosch & Hale rate coefficient:
!             Bosch & Hale, Nuc. Fus., vol32, no.4 (1992) 611
! NOTE: IF input ti > 100 kev then 100 kev is used !!!!!!!!!!!!!!!!!!!!!!
! ----------------------------------------------------------------------
      IMPLICIT  NONE 
      REAL(DP) mrcsq,c1, c2, c3, c5, b_gsq,ti,theta,xsi,tidum
!
!     DATA for d(d,n)he3
!
      DATA  c1, c2, c3, c5, b_gsq, mrcsq              &
         / 5.43360e-12,  5.85778e-3, 7.68222e-3,      &
           -2.96400e-06, 985.7716, 937814.0 /
!
      IF (ti  .LT. 0.2_DP ) THEN
        ddnrate = zeroc
      ELSE IF (ti .LE. 100.0_DP) THEN
        theta   = ti*C2/(1.0_DP+ti*(C3+ti*C5))
        theta   = ti/(1.0-theta)
        xsi     = (B_gsq/(4.0_DP*theta))**(0.33333333334_DP)
        ddnrate = C1 * theta * SQRT (xsi/(mrcsq*ti**3))  *  EXP (-3.0_DP*xsi)
      ELSE
        tidum = 100.
        theta   = tidum*C2/(1.0_DP+tidum*(C3+tidum*C5))
        theta   = tidum/(1.0_DP-theta)
        xsi     = (B_gsq/(4.0_DP*theta))**(0.33333333334)
        ddnrate = C1 * theta * SQRT (xsi/(mrcsq*tidum**3))  *  EXP (-3.0_DP*xsi)
      END IF
      ddnrate = 1.e-6*ddnrate                              ! convert to m^3/sec
      RETURN
!
      END  FUNCTION ddnrate




      REAL(DP) FUNCTION ddprate (ti)
!  ------------------------------------------------- 10/29/07 --- HSJ ---
!             returns rate(m**3/sec) of d(d,p)t reaction
!             (the division by 2 for collisions between like
!              particles is NOT accounted for here).
!             NEW Bosch & Hale rate coefficient:
!             Bosch & Hale, Nuc. Fus., vol32, no.4 (1992) 611
! NOTE: IF input ti > 100 kev then 100 kev is used !!!!!!!!!!!!!!!!!!!!!!
! ---------------------------------------------------------------------------
      IMPLICIT  NONE
      REAL(DP) mrcsq,C1,C2,C3,C5,B_gsq,theta,xsi,ti,tidum

!     DATA for d(d,p)t:

      DATA  C1,C2,C3,C5,B_gsq, mrcsq                        &
         / 5.65718e-12,  3.41267e-3, 1.99167e-3,            &
           1.05060e-05, 985.7715, 937814.0 /       

      IF (ti .LT. 0.2_DP) THEN
        ddprate = zeroc
      ELSE IF (ti .LE. 100.0_DP) THEN
        theta   = ti*C2/(1.0_DP+ti*(C3+ti*C5))
        theta   = ti/(1.0_DP-theta)
        xsi     = (B_gsq/(4.0_DP*theta))**(0.33333333334_DP)
        ddprate = C1 * theta * SQRT (xsi/(mrcsq*ti**3)) *  EXP (-3.0_DP*xsi)
      ELSE
        tidum = 100._DP
        theta   = tidum*C2/(1.0_DP+tidum*(C3+tidum*C5))
        theta   = tidum/(1.0_DP-theta)
        xsi     = (B_gsq/(4.0_DP*theta))**(0.33333333334_DP)
        ddprate = C1 * theta * SQRT (xsi/(mrcsq*tidum**3)) *  EXP (-3.0_DP*xsi)
      END IF
      ddprate = 1.e-6*ddprate                    ! convert to m^3/sec
      RETURN

      END FUNCTION ddprate


      REAL(DP) FUNCTION tt2nrate (tikevc)
!
! ------------------------------------------------- 10/29/07 --- HSJ ---
! FUNCTION returns the t(t,2n)he4 reaction rate,m**3/sec (averaged
! over a Maxwellian distribution at temperature tikev)
! This SUBROUTINE is a simple linear itnterpolaton of measured? and calc.?
! reaction rates "Atomic Data for Controled Fusion Research",ORNL-5207,
! E.W. Thomas et al. .NOV. 1979 for energies ge 6keV
! for the lower energies the NRL formulary table was used
!                         NOTE
!      THERE APPEARS TO BE CONSIDERABLE UNCERTAINTY IN T(T,2N)HE4 REACTION
!      RATES SO THESE VALUES MAY BE QUITE CRUDE.
! input
!   tikev       ion temperature in keV
!  NOTE 400 kev maximum is used !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ----------------------------------------------------------------------

      IMPLICIT NONE
      INTEGER(I4B),PARAMETER ::         ndt = 12 

      INTEGER(I4B)    ndtlo, ndtp1, i
      REAL(DP)     ekev_table(ndt), rttrate_table(ndt), tikev,    &
                   tikevc,dele,delr

      DATA    ndtlo/0/ ! ndtlo should be saved between calls
!                        DATA statement just sets it for very first CALL

      DATA   (ekev_table(i), rttrate_table(i), i=1,ndt)/         &
                    0.0,     0.0,                                &
                    1.0,     3.3e-22,                            &
                    2.0,     7.1e-21,                            &
                    5.0,     1.4e-19,                            &
                    6.0,     2.1e-19,                            &
                   10.0,     6.8e-19,                            &
                   20.0,     2.4e-18,                            &
                   40.0,     6.6e-18,                            &
                   70.0,     1.3e-17,                            &
                  100.0,     2.0e-17,                            &
                  200.0,     4.1e-17,                            &
                  400.0,     7.4e-17                             &
                                      /


       tikev = MIN(tikevc, ekev_table(ndt))
       ! find tikev in table
        CALL tableintrp (ekev_table, ndt, tikev, ndtlo)
          ndtp1   = ndtlo + 1
          dele    = ekev_table(ndtp1)-ekev_table(ndtlo)
          delr    = rttrate_table(ndtp1)-rttrate_table(ndtlo)
          tt2nrate = rttrate_table(ndtlo)+  (delr/dele)*(tikev-ekev_table(ndtlo))
          tt2nrate = 1.d-6*tt2nrate    ! convert to m^3/sec
      RETURN
      END  FUNCTION tt2nrate



      REAL(DP) FUNCTION dtnrate (ti)
! ------------------------------------------------- 10/29/07 --- HSJ ---
!     returns rate(m**3/sec) of t(d,n)he4 reaction
!     NEW Bosch & Hale rate coefficient:
!     Bosch & Hale, Nuc. Fus., vol32, no.4 (1992) 611
! ----------------------------------------------------------------------
      USE NRTYPE,                           ONLY : DP,I4B
      IMPLICIT  NONE
      REAL(DP) theta,xsi,C1,C2,C3,C4,C5,C6,C7, B_gsq, mrcsq ,ti  
!     DATA for t(d,n)he4:

      DATA  C1,C2,C3,C4,C5,C6,C7, B_gsq, mrcsq                        &
         / 1.17302e-9, 1.51361e-2, 7.51886e-2, 4.60643e-3,1.3500e-2,  &
          -1.06750e-4, 1.36600e-5,   1.182170e+3, 1124656.0 /

!     neutrons produced by bulk plasma d-t fusion:
!
      theta  = ti*(C2 +ti*(C4+ti*C6))
      theta  = theta/(1.0+ti*(C3+ti*(C5+ti*C7)))
      theta  = ti/(1.0-theta)
      xsi    = (B_gsq/(4.0*theta))**(0.33333333334)
      dtnrate = C1 * theta * SQRT (xsi/(mrcsq*ti**3)) * EXP (-3.0*xsi)
      dtnrate = dtnrate*1.D-6                            ! convert to m^3/sec
      RETURN

      END  FUNCTION dtnrate



      REAL*8 FUNCTION hdprate (ti)
! ------------------------------------------------- --- HSJ ---
!             returns rate(m**3/sec) of he3(d,p)he4 reaction
!             New Bosch & Hale rate coefficient:
!             Bosch & Hale, Nuc. Fus., vol32, no.4 (1992) 611
!  NOTE returns value for timax if ti > timax !!!!!!!!!!!!!!!!!!!!
! --------------------------------------------------------------


      IMPLICIT  NONE

      REAL(DP) mrcsq,ti,C1,C2,C3,C4,C5,B_gsq
      REAL(DP) theta,xsi,tidum,timax
!
!     data for he3(d,p)he4:

      DATA  C1,C2,C3,C4,C5,B_gsq, mrcsq                       &
        /  5.51036e-10, 6.41918e-3, -2.02896e-3,              &
           -1.91080e-05, 1.35776e-4, 4726.672, 1124572 /
      DATA timax /190.0/

      tidum = MIN(ti,timax)
      IF (tidum .LT. 0.5_DP) THEN
        hdprate = zeroc
      ELSE 
        theta  = tidum*(C2+ti*C4)/(1.0_DP+ti*(C3+ti*C5))
        theta  = tidum/(1.0_DP-theta)
        xsi    = (B_gsq/(4.0_DP*theta))**(0.33333333334_DP)
        hdprate = C1 * theta * SQRT (xsi/(mrcsq*tidum**3)) *  EXP (-3.0_DP*xsi)
      END IF
      hdprate = hdprate *1.e-6
      RETURN

      END   FUNCTION hdprate



   SUBROUTINE LOAD_FUSION_SOURCE
!-------------------------------------------------------------------------------
! -- load array stfuse if thernmal fusion calcs internal to gcnmp are requested
! -- (eg internal_thermal_fusion = true)
! -- sdburn is burnup rate of deuterium
! -- stburn is burnup rate of tritium
! -- he3,he4 and t produced in fusion are fast ions which are currently
! -- not treated consistently (eg no slowing down sources are included for these
! -- fusion progenated ions).
!-------------------------------------------------------------------------------
      USE source_terms_gcnmp,  ONLY : stfuse ,sfus,sdfuse,sdtfuse


         IF( .NOT. differential_burnup )THEN  ! assume burnup due to dt only
            ! subtract thermal part out from sfus (leaving beam-beam and beam thermal)
            ! on first pass this will be the iterdb rate. This is a workaround for he
            ! fact hat beam-beam and beam-themal is not separately avaialble:
            sfus(:)   = sfus(:) - stfuse(:)
            sdfuse(:) = dtn_fus(:)
            stfuse(:) = sdfuse(:)
            ! now add in thermal fusion:
            sfus(:) = sfus(:) + stfuse(:)
            sdtfuse(:) = sdfuse(:)
         ELSE
            sdfuse(:) = ddn_fus(:)  + dtn_fus(:) + dhe_fus(:) + ddp_fus(:)
            stfuse(:) = tt2n_fus(:) + dtn_fus(:)
            sdtfuse(:) = zeroc
         ENDIF


   END SUBROUTINE LOAD_FUSION_SOURCE


   SUBROUTINE load_thermal_fusion_energy_source
!----------------------------------------------------------------------------------------------
! -- load arrays  qtfuse,qtfusi  if thernmal fusion calcs internal to gcnmp
! -- are requested (eg internal_thermal_fusion = true - the default)
! -- he3,he4 and t produced in fusion are fast ions which are currently
! -- not treated consistently (eg no slowing down particle or energy sources are 
! -- included for these fusion  progenated ions).
! -- To get the frictional heating that occurs during the alpha particle
! -- slowing down we need to have a value for the average stored energy density,
! -- walp,kev/m^3. The initial value of walp comes from iterdb file each time
! -- gcnmp is called. Thereafter we evovle the equation
! --      dw_alpha/dt + w_alpha/taus_eng = sfe
! -- with w_alpha equal to the fast alpha energy density at time t. (The inital
! -- value of w_alpha is walp), sfe is the prompt energy source 
! -- (= 3.5Mev* thermal D-T fusion rate density)
! -- and taus_eng is the characteristic energy transfer time.
! -- Hence the buildup and decay of the alpha stored energy density as function
! -- of time at each radial position is determiend without other losses being
! -- taken into account(pitch angle scattering,etc).
! -- The thermal fusion heating rates are then
! -- qtfuse = fe*w_alpha/taus_eng
! -- qtfusi = fi*w_alpha/taus_eng
! -- We note  that for large t we have w_alpha = sfe*taus_eng 
! -- where fi and fe the fractions of the energy delivered to ions and electrons respectivly.
!-----------------------------------------------------------------------------------------------
      USE source_terms_gcnmp,                        ONLY : qtfuse,qtfusi
      USE fast_ion_data_gcnmp,                       ONLY : walp,w_alpha
      USE grid_class,                                ONLY : nj
      USE solcon_gcnmp,                              ONLY : time,time0
      USE solcon_gcnmp,                              ONLY : steady_state
         IMPLICIT NONE
         REAL(DP) sfe,t_el,ge,gi,taus_eng,dfctr,xpont
         INTEGER(I4B) j

                                                      
         IF(steady_state ==  1)THEN
            t_el = time -time0            ! elapsed time since walp was set in
                                          ! input iterdb file. Note that time0
                                          ! is the start time of gcnmp as set
                                          ! in namelist input and time is
                                          ! current time so this assumes that
                                          ! alpha slowing down rate is not
                                          ! asymptotic
           IF(alpha_slow == 'asymptotic')   &
           t_el = 1.e6  ! 8888999999      !assumes alpha slowing down rate is asymptotic.
         ELSE
             t_el = 1000000._DP
         ENDIF

         DO j =1,nj
            CALL slowd_parms(taus_eng,ge,gi,j)
            sfe        = 3.5e3*dtn_fus(j)    ! prompt alpha energy density,Kev/(m^3 sec)
            xpont      = t_el/taus_eng
            dfctr      =  0.0_DP
            IF(xpont .LT. 50._DP) dfctr      = EXP(-xpont)
            w_alpha(j) = walp(j)*dfctr + sfe*taus_eng*(1._Dp - dfctr)
            qtfuse(j)  = ge*w_alpha(j)/taus_eng
            qtfusi(j)  = gi*w_alpha(j)/taus_eng

  IF( ABS(qtfuse(j)) > 1.e25)THEN 
     print *,'qtfuse(j),j=',qtfuse(j),j
     print *,'ge,w_alpha(j),taus_eng =',ge,w_alpha(j),taus_eng
     print *,'walp(j),dfctr,sfe =',walp(j),dfctr,sfe
  ENDIF
  IF( ABS(qtfusi(j)) > 1.e25)THEN 
     print *,'qtfusi(j),j=',qtfusi(j),j
     print *,'gi,w_alpha(j),taus_eng =',gi,w_alpha(j),taus_eng
  ENDIF
         ENDDO
   END SUBROUTINE load_thermal_fusion_energy_source


   SUBROUTINE slowd_parms(taus_eng,ge,gi,j)
!----------------------------------------------------------------------------------------------
! -- Appoxiamte fast ion slowing down parameters, see for example
! -- Callen,Controlled Nucl.Fusion Research Vol I,645(1974)
! -- Note that this is specialized to alpha paricles 
! -- (eq no charge exchange  losses and a= 4)
!----------------------------------------------------------------------------------------------
         USE common_constants,                         ONLY : pi,Proton_mass,       &
                                                              Electron_Rest_Mass,   &
                                                              zeroc,sqrt2,          &
                                                              Permittivity,rootpi,  &
                                                              Electron_Charge
         USE ions_gcnmp,                               ONLY : nion,zsq,atw,zeff

         USE dep_var,                                  ONLY : en,ene,te

         IMPLICIT NONE
         REAL(DP) ecrit,const,t_el,ge,gi,taus_eng,one_thrd,two_thirds,   &
                  sum2,sum1,atwa,erel0,ezero,erot,tconst,zsqa,tej,tegt0,      &
                  taue,xlam,one_third,taus,vrat,taurat
         INTEGER(I4B) j,k

         one_third = 1._DP/3._DP ; two_thirds = 2._DP*one_third
         const     = 0.5_DP * (4.5_DP * pi * Proton_mass/Electron_Rest_Mass) ** one_third  ! ~ 14.8
         atwa      = 4._DP   ; zsqa = 4._DP
         ezero     = 3.5e3_DP                                      ! kev,  prompt alpha energy



         !---------------------------------------------------------------------------------------
         ! -- taus is  the Spitzer momentum exchange time for electron-ion collisions
         ! -- taue is  electron effective collison time with all ions
         ! -- ecrit is   the "critical energy" at which fast ions
         ! -- transfer equal energy to electrons and thermal ions;
         ! -- emzrat = mi*<Z>/(mf*[Z]) 
         ! -- xlam .Coulomb logarithm, from NRL Plas Form. valid for Ti/1836 < 10 Z^2 (in ev) < Te 
         !---------------------------------------------------------------------------------------
             tegt0   = MAX (te(j), 0.0001_DP)                 ! keV, avoid probs near plasma edge
             tej     = 1.6021765e-16*tegt0                    ! jolues
             xlam    = 24.0_DP  - LOG (SQRT (1.e-6_DP*ene(j))/(1.0e3_DP*tegt0)) 
             taue    = 12._DP*pi*rootpi*Permittivity*Permittivity*           &
                            SQRT(Electron_Rest_Mass)*tej**1.5_DP/            &
                            (sqrt2*ene(j)*zeff(j)*Electron_Charge**4 *xlam)
             tconst  = atwa * Proton_mass/(zsqa*Electron_Rest_Mass)
             taus    = tconst*zeff(j)*taue

         sum1      = zeroc
         sum2      = zeroc
         DO  k=1,nion
           sum1    = sum1 + en(j,k)*zsq(j,k)
           sum2    = sum2 + en(j,k)*zsq(j,k)/atw(k)
         ENDDO
!         emzrat(j) = sum1 / (atwa*sum2)
         ecrit     = const*te(j)*atwa*(sum2/ene(j))**two_thirds       !kev
         erot      = zeroc                                            !neglect bulk plasma rotation
         erel0     = ABS (ezero - erot)
         vrat      = SQRT (ecrit/erel0)
         !taurat is zero if we neglect charge exchange during slowing down
         taurat    = zeroc  ! for alpha particles the ratio  taus/taupcx =0
         ge        = gef_alpha(vrat,taurat)*erel0/ezero
         gi        = erel0/ezero - (1.0_DP + 0.5_DP * taurat)*ge 
         IF (gi  .LT. zeroc)  gi  = zeroc
         taus_eng  = taus*ge*0.5_DP
   END SUBROUTINE slowd_parms


   REAL(DP)  FUNCTION gef_alpha (vpar, tpar)

      IMPLICIT  NONE
      REAL(DP) root13,rooty,ATAN,y,term2,term3,term4,fact4,vpar,tpar

! ----------------------------------------------------------------------
!  This routine evaluates the function Ge, which is related to the
!     transfer of energy from fast ions slowing down on electrons.  The
!     function is defined by Callen et al., IAEA Tokyo, Vol. I, 645
!     (1974). It is specialized to alpha particles.
! ----------------------------------------------------------------------
 

       root13 = SQRT (1.0_DP/3.0_DP)
       fact4  = root13 * ATAN (root13)


        rooty = 1.0_DP / vpar
        y     = rooty**2
        term2 = LOG ((1.0_DP-rooty+y)/(1.0_DP + rooty)**2)/(6.0_DP * y)
        term3 = root13* ATAN (root13*(2.0_DP*rooty-1.0_DP))/y
        term4 = fact4/y
        gef_alpha   = 2.0_DP * (0.5_DP-term2-term3-term4)

      END FUNCTION gef_alpha

      SUBROUTINE tot_thermal_fusion_power
!-------------------------------------------------------------------------------
!-- Get volume intetal of thermal fusion hetaing rate
!-------------------------------------------------------------------------------
      USE source_terms_gcnmp,                        ONLY : qtfuse,qtfusi

      IF(.NOT. ALLOCATED(qfuse_thermal))THEN
         ALLOCATE(qfuse_thermal(SIZE(qtfuse)))
      ELSE
         IF(SIZE(qfuse_thermal) .NE. SIZE(qtfuse))THEN
            DEALLOCATE(qfuse_thermal)
            ALLOCATE(qfuse_thermal(SIZE(qtfuse)))
         ENDIF
      ENDIF

       qfuse_thermal(:) = qtfuse(:) + qtfusi(:)
       CALL trapv(r,qfuse_thermal,hcap,nj,pfuse_tot)
       pfuse_tot  = pfuse_tot*volfac*joupkev
      
      END       SUBROUTINE tot_thermal_fusion_power


   END MODULE fusion_gcnmp



