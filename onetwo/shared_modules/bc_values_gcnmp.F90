  MODULE BC_values_gcnmp

    USE nrtype,             ONLY : DP,I4B
    USE solcon_gcnmp,       ONLY : itran_max
    IMPLICIT NONE
    REAL(DP),SAVE,ALLOCATABLE, DIMENSION(:) :: ub,ub_lcl,zeff_bc,ene_bc,    &
                  te_bc,ti_bc,angrot_bc,bctime,fix_edge_te,fix_edge_ti,     &
                  fix_edge_rot, mult_flux,mult_den,gammat_bc,dudt0,dudt0_save
    REAL(DP),SAVE,ALLOCATABLE, DIMENSION(:,:) :: fix_edge_ni,flux_bc
    REAL(DP),SAVE,ALLOCATABLE, DIMENSION(:,:) :: en_bc,bc,bctime_zone
    REAL(DP),SAVE      :: totcur_bc,fix_edge_te_bc,fix_edge_ti_bc,    &
                          fix_edge_rot_bc,time_bc, vloop_bc         
 
    INTEGER(I4B),SAVE  :: nbctim,te_var_edge,ti_var_edge,rot_var_edge, u_vloop_bc, &
                          nb_robin,nb_robin_nj
    INTEGER(I4B),SAVE  :: te_index,ti_index,rot_index,iter_vloop,axis_bc_type
    INTEGER(I4B),SAVE,ALLOCATABLE,DIMENSION(:)   :: ni_index, ni_var_edge,bc_type 
    REAL(DP),SAVE,ALLOCATABLE,DIMENSION(:)       :: fix_edge_ni_bc
    REAL(DP),SAVE,ALLOCATABLE,DIMENSION(:)       :: flux_ub_conserve
    REAL(DP),SAVE,DIMENSION(itran_max)           :: flux_mult_bc 
    LOGICAL, SAVE      :: set_bc, robin_bc,conserve_bc_flux


  CONTAINS
!no or clause in preprocessor?

#ifdef GCNMP
#define USESUBS
#endif

#ifdef NFREYA
#define USESUBS
#endif
#ifdef ONETWO
#define NOUSESUBS
#endif
#if defined NOUSESUBS
  SUBROUTINE bc_conditions_init
    !Onetwo does not call this
    ! but it is needed for linking
  END SUBROUTINE bc_conditions_init
     
  SUBROUTINE set_bc_type
   !Onetwo does not call this
   ! but it is needed for linking
  END SUBROUTINE set_bc_type
#endif

#if defined USESUBS

  SUBROUTINE bc_conditions_init
! ---------------------------------------------------------------------
! input
!  ntot
!  set_bc            logical, set to false on first call
!                    returned as true. initializes the boundary arrays
! OUTPUT: 
!  ub,ub_lcl
!  bc
!  bctime
!  fix_edge_**(nbctim) these are currently set to constant in time values
!  ni_var_edge(1:nion)
!  te_var_edge
!  ti_var_edge
!  rot_var_edge   = 0  no zonal boundary condition for this profile 
!                      was specified
!                 = 1  zonal boundary conditions given in terms 
!                      of grid points
!                 = -1 zonal boundary conditions given in terms         
!                      normalized flux
! ---------------------------------------------------------------------
      USE solcon_gcnmp,                               ONLY : ntot,tGCNMS,tGCNMf,nameu
      USE ions_gcnmp,                                 ONLY : nion,nprim,nimp
      USE plasma_properties,                          ONLY : profile, mhd_dat,get_element
      USE grid_class,                                 ONLY : nj,roa,njm1
      USe dep_var,                                    ONLY : dp4
      USE common_constants,                           ONLY : zeroc



      IMPLICIT NONE

      INTEGER j
      IF(.NOT. set_bc)THEN
         nbctim = 2_I4B
         IF(.NOT. ALLOCATED(ub))ALLOCATE(ub(ntot))
             ub(:) = zeroc
         IF(.NOT. ALLOCATED(ub_lcl))ALLOCATE(ub_lcl(ntot))
             ub_lcl(:) = zeroc
         IF(.NOT. ALLOCATED(dudt0))ALLOCATE(dudt0(ntot))  ! dudt at rho =0
             dudt0(:) = zeroc
         IF(.NOT. ALLOCATED(dudt0_save))ALLOCATE(dudt0_save(ntot))  ! dudt at rho =0
             dudt0_save =  zeroc
         IF(.NOT. ALLOCATED(bc))ALLOCATE(bc(nbctim,ntot))
              bc(nbctim,ntot) = zeroc
         IF(.NOT. ALLOCATED(bctime))ALLOCATE(bctime(nbctim))
             bctime(:) = zeroc
         IF(.NOT. ALLOCATED(fix_edge_te))ALLOCATE(fix_edge_te(nbctim))
            fix_edge_te(:) = zeroc
         IF(.NOT. ALLOCATED(fix_edge_ti))ALLOCATE(fix_edge_ti(nbctim))
            fix_edge_ti(:) = zeroc
         IF(.NOT. ALLOCATED(fix_edge_rot))ALLOCATE(fix_edge_rot(nbctim))
            fix_edge_rot(:) = zeroc
         IF(.NOT. ALLOCATED(fix_edge_ni))ALLOCATE(fix_edge_ni(nbctim,nion))
            fix_edge_ni(:,:) = zeroc
      ENDIF

      bctime(1)      = tGCNMs
      bctime(nbctim) = time_bc
      DO j=1,nion
         IF(bc_type(j) == 0) THEN                     ! bc(1:nbctim,1:ntot)
           bc(1,j) = mult_den(j)*profile%en(j)%data(nj)
           bc(2,j) = mult_den(j)*en_bc(nj,j)
         ELSE IF(bc_type(j) == -1)THEN
           bc(1,j) =  mult_flux(j)*profile%flux(j)%data(nj)
           bc(2,j) =  mult_flux(j)*flux_bc(nj,j)
         ELSE
           bc(1,j) = gammat_bc(j)/mult_flux(j)
           bc(2,j) = bc(1,j)                          ! constant in time
         ENDIF
      ENDDO

      bc(1,nion+1) = get_element(profile%te,nj)
      bc(2,nion+1) = te_bc(nj)
      bc(1,nion+2) = get_element(profile%ti,nj)
      bc(2,nion+2) = ti_bc(nj)
      IF(u_vloop_bc == 0_I4B)THEN
         bc(1,nion+3) = mhd_dat%tot_cur
         bc(2,nion+3) = totcur_bc
      ELSE ! vloop BC 
         bc(1,nion+3) = mhd_dat%vloop
         bc(2,nion+3) = vloop_bc
      ENDIF
      bc(1,nion+dp4) = get_element(profile%angrot,nj)
      bc(2,nion+dp4) = angrot_bc(nj)

      !here is where the  assumption of time independent
      !fix_edge parameters is used. We may allow for more
      !flexibility in the future by having the fix_edge_**_bc
      !values input as a function of time:
      fix_edge_te(1:nbctim)   = fix_edge_te_bc
      fix_edge_ti(1:nbctim)   = fix_edge_ti_bc
      fix_edge_rot(1:nbctim)  = fix_edge_rot_bc
      DO j = 1, nion
         fix_edge_ni(1:nbctim,j)   = fix_edge_ni_bc(j) ! may be different for each ion
                                                       ! but is constant in time for
                                                       ! all ions.
      ENDDO



      CALL bc_zone(bctime,nbctim,fix_edge_te,nj,               &
                         nameu(nion+1),roa,te_var_edge)

      CALL bc_zone(bctime,nbctim,fix_edge_ti,nj,               &
                         nameu(nion+2),roa,ti_var_edge)

      CALL bc_zone(bctime,nbctim,fix_edge_rot,nj,               &
                         nameu(nion+dp4),roa,rot_var_edge)

      IF(.NOT. ALLOCATED(ni_var_edge))ALLOCATE(ni_var_edge(nion))
      !For each ion species, 1..nprim+nimp(=nion) set ni_var_edge
      !to -1,0,1 according to how fix_edge_ni(j) is specifying boundary condition. 
      DO j = 1,nion
         CALL bc_zone(bctime,nbctim,fix_edge_ni(1,j),nj,nameu(j),roa,ni_var_edge(j))
      ENDDO
      

      njm1 = nj-1
      set_bc = .TRUE.

      RETURN

  END SUBROUTINE bc_conditions_init



  SUBROUTINE bc_conditions(t_bc)
!
! ---------------------------------------------------------------------
!  INPUT
!
!  t_bc              time at which bc is required
!  bctime(nbctim)    times at which bc is available
!  nbctim            number of  values in bctime (GCNMP assumes this is always 2)
!  bc(nbctim,k)      k=1,..ntot. The actual bc values at time bctime(nbctim)
!                    These are combinations of densities, particle fluxes,temperatures,
!                    rotation speed, total current, depending on index k. They are set 
!                    in sub bc_conditions_init to the values at j = nj
!  ntot              nprim+nimp+1+1+1+1 (ie te,ti,rbp,w)
!                    included and the 4 1's are for te,ti,rbp,rotation respectively.
!
!  te_var_edge          These are bc indicators  for each dependent variable
!  ti_var_edge          They are set in sub bc_zone (called by bc_conditions_init)
!  rot_var_edge         and have values of -1,0, or 1. If the value is zero then
!  ni_var_edge(1:nion)  we simply time interpolate bc(:,k) to get ub(k). Otherwise
!                       the boundary is at a rho location less than rho =a (grid point nj)
!                       and we have to make some adjustments.
!
!
!  OUTPUT
!    ub(k)             k=1,2...ntot the value of the boundary condition  at the point
!                      r=rho_edge (not necessarily at 1.0)
!    bctime_zone(j,k)  k=1,2...nion+dp4 , j = bc value to nj
!
! ------------------------------------------------------- HSJ-7/02/05--
!
      USE error_handler
      USE grid_class,                           ONLY : nj
      USE solcon_gcnmp,                         ONLY : usave,dt,ntot,iangrot,itran,dudtb
      USE ions_gcnmp,                           ONLY : nion
      USE grid_class,                           ONLY : roa
      USE plasma_properties,                    ONLY : profile
      USE Vector_class,                         ONLY : get_values
      USE io_gcnmp,                             ONLY : nlog
      USE curden_terms,                         ONLY : tot_cur
      USE dep_var,                              ONLY : dp4
      USE common_constants,                     ONLY : zeroc


      IMPLICIT  NONE 

      REAL(DP) t_bc,fun,diff1,diff2,IP_tor
      INTEGER(I4B) k,i1,i2,i2r,i1r,j,nivaredge

      REAL(DP) work(nj)                             !temporary arrays
      
      IF(.NOT. ALLOCATED(bctime_zone ))ALLOCATE(bctime_zone(nj,ntot))
      bctime_zone(:,:) = zeroc

  DEPVAR_LOOP: DO k=1,ntot
        ub(k) = bc(1,k)                    ! see bc_conditions init for defn of bc
        ub_lcl(k) = ub(k)                  ! reset below as necessary
        dudtb(k)  = zeroc                  ! reset below as necessary
        !for all profiles except rbp (pointer to rbp is ntot-iangrot):
        IF (  k .NE. ntot - iangrot)THEN
           CALL interp_bc (t_bc, bctime, nbctim, bc(1,k), ub(k)) ! return ub(k) at t = t_bc
           dudtb(k) = (bc(2,k)-bc(1,k))/(bctime(2)-bctime(1))    ! valid at j = nj
        ELSE 
           IF( u_vloop_bc  == 0_I4B )THEN   !rbp,bc is total current
              CALL interp_bc (t_bc, bctime, nbctim, bc(1,k), IP_tor)
              dudtb(k) = (bc(2,k)-bc(1,k))/(bctime(2)-bctime(1)) ! valid at j =nj
           ELSE !vloop bc
                !this boundary condition is non linear and
                !changes as the profiles change during the iteration
               iter_vloop = 0_I4B
              !get IP_tor from vloop_c:
              CALL get_IP_from_vloop(iter_vloop,t_bc)
              IP_tor = tot_cur
           ENDIF
           !IP_tor is boundary condition current at time t_bc.
           !convert this to boundary condition for dependent variable
           !rbp(nj) = G*H*r*Bp0:

           ub(k) = IP_tor*.2e-6    !.2e-6 = u0/2pi,ub = GHrBp0 in tesla m

        ENDIF





!--------------------------------------------------------------------------------
! -- All  profiles with bc at r(nj) are done, Now handle those that have
! -- bc at  r < r(nj):
!--------------------------------------------------------------------------------


        IF(nbctim .GT. 1)THEN           !NOTE we must have nbctim >= 2
             CALL find1 (i1, i2, t_bc, bctime, nbctim) ! get indecies that bride t_bc in bctime
              IF(i1*i2 .EQ. 0)THEN
                PRINT *,'ERROR, t_bc =',t_bc
                PRINT *,'bctime =',bctime(1:nbctim)
                lerrno = 25
                CALL terminate(lerrno,nlog)
              ENDIF

         ENDIF


!--densities:
        nivaredge=0
        IF(k .le. nion)nivaredge  =  ni_var_edge(k) 
        IF(k .LE. nion  .AND. nivaredge  .NE. 0 .AND. itran(k) &
                                                .EQ. 1)THEN  !variable edge location
           IF(nbctim .GT. 1)THEN                             !for ion species k

              !first get the index value, ni_index, which is the grid point
              !number at which  (and beyond which) ni is given by boundary conditions 
              IF(i1 .NE. i2)THEN
                 CALL linear (fix_edge_ni(i1,k),fix_edge_ni(i2,k), &
                         fun,bctime(i1),bctime(i2), t_bc)
              ELSE
                 fun = fix_edge_ni(i1,k)
              ENDIF
              IF(ni_var_edge(k) .EQ. 1)THEN
                  ni_index(k) = NINT(fun)
              ELSE  ! ni_var_edge(k) .eq. -1
                 !find the index of the r value closest to the normalized roa
                 CALL find1 (i1r, i2r, fun, roa, nj)
                 IF(i1r*i2r .EQ. 0)THEN
                   lerrno = 26
                   CALL terminate(lerrno,nlog)
                 ENDIF
                 !pick closest value:
                 diff1 =ABS(fun -roa(i1r))
                 diff2 =ABS(fun -roa(i2r))
                 IF(diff1 .GE. diff2)THEN
                    ni_index(k) = i2r
                 ELSE
                    ni_index(k) = i1r
                  ENDIF
              ENDIF
                !Linearly interpolate the two profiles
                ! (one at the start time and the other  at  bctime)
                ! at each grid point from rho(ni_index)  to rho(nj)
                ! and save it in bctime_zone:
                ! There are 3 type of bc for densities:
                ! (a) given flux at ni_index(k)
                ! (b) given density at ni_index(k)
                ! (c) given alpha*density + beta*flux (eq mixed bc)
                ! for cases (a) and (c) it is assumed that the input bc
                ! values are stored in flux_bc(1:nj,1:nion)
                IF(bc_type(k) == -1_I4B  )THEN
                   !flux bc
                    work(1:nj) = get_values(profile%flux(k))  
                    !time interpolate profile to get values at time = t_bc.
                    DO j = ni_index(k), nj
                         CALL linear (work(j),flux_bc(j,k), &
                         fun,bctime(i1),bctime(i2), t_bc)
                         bctime_zone(j,k) = fun*mult_flux(k)
                    ENDDO
                ELSE IF(bc_type(k) == 0)THEN
                    !density bc . Load density at time bctime(i1):
                    work(1:nj) = get_values(profile%en(k))  
                    !time interpolate profile to get values at time = t_bc.
 
                    DO j = ni_index(k), nj
                         CALL linear (work(j),en_bc(j,k), &
                         fun,bctime(i1),bctime(i2), t_bc)
                         bctime_zone(j,k) = fun*mult_den(k)
                    ENDDO
                ELSE IF(bc_type(k) == 1)THEN
                   ! mixed flux and density
                   DO j = ni_index(k), nj ! currently no time dependance
                      bctime_zone(j,k) = gammat_bc(k)/mult_flux(k)
                   ENDDO
                ENDIF

           ENDIF !     (nbctim branches )
           !bctime_zone(ni_index:nj,k) now holds  the boundary values of ni profile
           !(pointed to by k)

           dudtb(k) = (bctime_zone(ni_index(k),k) - usave(k,ni_index(k)))/dt
           ub_lcl(k) = bctime_zone(ni_index(k),k)
           ub(k) = bctime_zone(nj,k)  
         ENDIF   !(ni_var_edge .ne. 0)

!--Te:
        !determine values in the  edge zone at this time by linear
        !interpolation  **_var_edge = -1,0,1 are valid inputs here
        IF(k .EQ. nion +1 .AND. te_var_edge .NE. 0 .AND. itran(k) .EQ. 1)THEN  !variable edge

           IF(nbctim .GT. 1)THEN               

              !first get the index value, te_index, which is the grid point
              !number at which  (and beyond which) te is given by boundary conditions 

              IF(i1 .NE. i2)THEN
                 CALL linear (fix_edge_te(i1),fix_edge_te(i2), &
                         fun,bctime(i1),bctime(i2), t_bc)

              ELSE
                 fun = fix_edge_te(i1)
              ENDIF


              IF(te_var_edge .EQ. 1)THEN
                  te_index = NINT(fun)
              ELSE  ! te_var_edge .eq. -1
                 !find the index of the r value closest to the normalized fl
                 CALL find1 (i1r, i2r, fun, roa, nj)
                 IF(i1r*i2r .EQ. 0)THEN
                   lerrno = 26
                   CALL terminate(lerrno,nlog)
                 ENDIF
                 !pick closest value:
                 diff1 =ABS(fun -roa(i1r))
                 diff2 =ABS(fun -roa(i2r))
                 IF(diff1 .GE. diff2)THEN
                    te_index = i2r
                 ELSE
                    te_index = i1r
                  ENDIF
              ENDIF

                !Linearly interpolate the two profiles in 
                !time at each grid point from te_index to nj
                !required profile and save it in bctime_zone:
                work(1:nj) = get_values(profile%te)

                DO j = te_index, nj
                  CALL linear (work(j),te_bc(j), &
                         fun,bctime(i1),bctime(i2), t_bc)
                  bctime_zone(j,nion+1) = fun

                ENDDO

           ENDIF !     (nbctim branches )
           !bctime_zone(te_index:nj,nion+1) now holds  the boundary values of te profile
           !(pointed to by nion+1)
           dudtb(k) = (bctime_zone(te_index,k) - usave(k,te_index))/dt
           ub_lcl(k) = bctime_zone(te_index,k)
           ub(k) = bctime_zone(nj,nion+1)  
 
         ENDIF   !(te_var_edge .ne. 0)
 
!--Ti:
        IF(k .EQ. nion +2 .AND. ti_var_edge .NE. 0 .AND. itran(k) .EQ. 1)THEN
          IF(nbctim .GT. 1)THEN
             !first get the index value, ti_index:
             IF(i1 .NE. i2)THEN
                CALL linear (fix_edge_ti(i1),fix_edge_ti(i2), &
                     fun,bctime(i1),bctime(i2), t_bc)
             ELSE
                fun = fix_edge_ti(i1)
             ENDIF
             IF(ti_var_edge .EQ. 1)THEN
                ti_index = NINT(fun)
             ELSE  ! ti_var_edge .eq. -1
                !find the index of the r value closest to the normalized fl
                CALL find1 (i1r, i2r, fun, roa, nj)
                IF(i1r*i2r .EQ. 0)THEN
                   lerrno = 28
                   CALL terminate(lerrno,nlog)
                ENDIF
                !pick closest value:
                diff1 =ABS(fun -roa(i1r))
                diff2 =ABS(fun -roa(i2r))
                IF(diff1 .GE. diff2)THEN
                   ti_index = i2r
                ELSE
                   ti_index = i1r
                ENDIF
             ENDIF
             !next get the profile at this time
             work(1:nj) = get_values(profile%ti)
             DO j = ti_index, nj
                CALL linear (work(j),ti_bc(j), &
                     fun,bctime(i1),bctime(i2), t_bc)
                bctime_zone(j,nion+2) = fun
             ENDDO
             dudtb(k) = (bctime_zone(ti_index,k) - usave(k,ti_index))/dt
             ub_lcl(k) = bctime_zone(ti_index,k)
           ENDIF
           ub(k) = bctime_zone(nj,nion+2)
 
        ENDIF !ti_var_edge .ne. 0

        IF(k == ntot-1 .AND. itran(k) ==1)THEN
           ! MACCORMACK REQUIRES COMMON EDGE FOR ALL PROFILES
           ! THIS EDGE IS DEFINED BY THE TURBULENT TRANSPORT MODELS
           ! load bctime_zone for maccormack solver (not used for newton method)
           ! Hence we may have boundary for rbp at rho < 1.0 for this case
           ! codeing not in place yet ! 88888999999
        ENDIF

!--Tor rot:
        IF(k .EQ. ntot .AND. iangrot .EQ. 1 .AND. &
                            rot_var_edge .NE. 0 .AND. itran(k) .EQ. 1)THEN
          IF(nbctim .GT. 1)THEN
           !determine values in the  edge zone at this time by linear
           !interpolation from angrot_edge_zone

           !first get the index value rot_index:
           IF(i1 .NE. i2)THEN
              CALL linear (fix_edge_rot(i1),fix_edge_rot(i2), &
                         fun,bctime(i1),bctime(i2), t_bc)
           ELSE
              fun = fix_edge_rot(i1)
           ENDIF
           IF(rot_var_edge .EQ. 1)THEN
                  rot_index = NINT(fun)
           ELSE  ! rot_var_edge .eq. -1
             !find the index of the r value closest to the normalized fl
              CALL find1 (i1r, i2r, fun, roa, nj)
              IF(i1r*i2r .EQ. 0)THEN
                lerrno = 30
                CALL terminate(lerrno,nlog)
              ENDIF

              !pick closest value:
              diff1 =ABS(fun -roa(i1r))
              diff2 =ABS(fun -roa(i2r))
              IF(diff1 .GE. diff2)THEN
                 rot_index = i2r
              ELSE
                 rot_index = i1r
              ENDIF
           ENDIF
           !next get the profile at this time
              work(1:nj) = get_values(profile%angrot)
              DO j = rot_index, nj
                  CALL linear (work(j),angrot_bc(j), &
                         fun,bctime(i1),bctime(i2), t_bc)
                  bctime_zone(j,nion+dp4) = fun
              ENDDO
           ENDIF
            dudtb(k) = (bctime_zone(rot_index,nion+dp4) - usave(k,rot_index))/dt
            ub_lcl(k) = bctime_zone(rot_index,nion+dp4)
            ub(k) = bctime_zone(nj,nion+dp4)  !reset from above for consist
        ENDIF
  END DO DEPVAR_LOOP      !  loop over all dependent variables, k=1,ntot


  CALL check_zero_bc(t_bc)



  RETURN
  END   SUBROUTINE bc_conditions

      SUBROUTINE boundary_flux(flux_ub,ntot,nb,ndep_var)
!-------------------------------------------------------------------------------
!   Get the steady state flux at the  boundary point nb for variable
!   ndep_var. Note that ndep_var = 1.. ntot
!   flux = [Integral from o to rho(nb) of S*H*rho*drho]/(h(nb)*rho(nb))
!   NOTE that this routine asumes that  the source term,s, is  loaded
!   Hence call this only affter call to source is made.
!------------------------------------------------------------HSJ-------------
      USE nrtype,                                    ONLY : DP,I4B

      USE common_constants,                          ONLY : zeroc,kevpjou,Proton_Mass

      USE dep_var,                                   ONLY : ene,te,ti,en,angrot

      USE MPI_data,                                  ONLY : myid,master

      USE grid_class,                                ONLY : r,fcap,gcap,hcap,nj,r2capi

      USE source_terms_gcnmp,                        ONLY : s,imp_called,qdelt

      USE error_handler,                             ONLY : lerrno,terminate
                                                   
      USE ions_gcnmp,                                ONLY : nimp,nion,z,atw,dzdte

      USE io_gcnmp,                                  ONLY : nlog,ncrt

      USE solcon_gcnmp,                              ONLY : dudt,time_deriv_in_flux,   &
                                                            angrcple,iangrot

      IMPLICIT NONE
      INTEGER(I4B) nb,ndep_var,ntot,j,k,l
      REAL(DP) flux_ub(ntot),stot(nj),td(nj)
      REAL(DP) fmult,sume,sumi,r2omega,smassden,xmassi


          td(:) = zeroc
          IF(1 .LE. ndep_var .AND. ndep_var .LE. ntot)THEN 

              IF(time_deriv_in_flux)THEN
                 DO j=1,nj
                    td(j) = zeroc
                    IF(ndep_var .LE. nion)THEN     ! density equations
                       td(:) = dudt(ndep_var,:)    ! em*dudt row ndep_var, diagonal for densities
                    ELSEIF(ndep_var == nion+1)THEN ! Te equation
                       DO k=1,ntot                 ! em*dudt row nion+1,col k
                          IF(k .LE. nion)THEN
                             td(j) = td(j) + 1.5_DP*te(j)*z(j,k)*dudt(k,j)
                          ELSEIF( k == nion+1)THEN
                             sume = ene(j)
                             DO l=1,nion
                                sume = sume + en(j,l)*te(j)*dzdte(j,l)
                             ENDDO
                             td(j) = td(j) + 1.5_DP*sume*dudt(k,j)
                          ENDIF ! no contribution if k > nion +1
                       ENDDO
                    ELSEIF(ndep_var == nion+2)THEN   ! TI equation
                       fmult   = angrcple*angrot(j)*kevpjou*0.5_DP*iangrot
                       r2omega = angrot(j)*r2capi(j)
                       DO k=1,ntot                   ! em*dudt row nion+1,col k
                          IF(k .LE. nion)THEN
                             td(j) = td(j) + (1.5_DP*ti(j)+fmult*r2omega*atw(k)*Proton_Mass)*dudt(k,j)
                          ELSEIF(k == nion+1)THEN
                             td(j) = zeroc
                          ELSEIF( k == nion+2)THEN
                             sumi = zeroc
                             DO l=1,nion
                                sumi = sumi + en(j,l)
                             ENDDO
                             td(j) = td(j) + 1.5_DP*sumi*dudt(k,j)
                          ELSEIF( k == nion +3)THEN
                             td(j) = zeroc
                          ELSEIF( k == nion +4)THEN
                             smassden = 0.0_DP
                             DO l=1,ntot-4
                                xmassi   = atw(l)*Proton_mass     
                                smassden = smassden+xmassi*en(j,l)
                                td(j)    = td(j) + r2omega*xmassi*dudt(l,j)
                             ENDDO
                             td(j)    = td(j) + smassden*r2omega*kevpjou*r2capi(j) *dudt(k,j)
                          ENDIF
                       ENDDO

                    ELSEIF(ndep_var == nion+3 .AND. j .GT. 1 )THEN  ! rbp equation
                          td(j) = (1.0_DP / (fcap(j) * gcap(j) * (hcap(j) * r(j))**2 ))*dudt(ndep_var,j)
                    ELSEIF(ndep_var == nion+4)THEN
                          smassden = 0.0_DP
                          td(j) = zeroc
                          DO k=1,ntot-4
                             xmassi   = atw(k)*Proton_mass     
                             smassden = smassden+xmassi*en(j,k)
                             td(j)    = td(j) + r2omega*xmassi*dudt(k,j)
                          ENDDO
                             td(j)    = td(j) + smassden*r2capi(j) *dudt(ndep_var,j)
                   ENDIF
                 ENDDO  ! loop over grid points 
              ENDIF     ! time_deriv_in_flux

              CALL impsrc   ! load qdelt if not loaded 
              stot(:) = s(ndep_var,:) - td(:)  ! td == 0 if time_deriv_in_flux = false
              IF(ndep_var == nion+1)stot(:) = stot(:) -qdelt(:) 
              IF(ndep_var == nion+2)stot(:) = stot(:) +qdelt(:)              
              CALL trapv (r,stot,hcap, nb, flux_ub(ndep_var))   ! utils .F90
              flux_ub(ndep_var) = flux_ub(ndep_var)/(hcap(nb)*r(nb))
write(173,FMT='("flux_ub =",x,1pe12.2)')flux_ub(ndep_var)
           ELSE
              IF(myid == master)WRITE(ncrt,FMT='("ERROR in sub ss_boundary_flux")')
              lerrno = 69
              CALL terminate(lerrno,nlog)
           ENDIF

        RETURN

      END SUBROUTINE boundary_flux





      SUBROUTINE interp_bc (xval, x, n, y, yval)
      USE nrtype,           ONLY : DP,I4B
      IMPLICIT NONE
      INTEGER(I4B) n,i1,i2
      REAL(DP) xval,yval,x(n),y(n)


!  ----------------------------------------------------------------------
! --- same as interp but:
!     a) if xval is .le. x(1) returns yval = y(1) (instead if 0.0)
!     b)            .ge. x(n)              = y(n)
!
! ----------------------------------------------------------------------
!
 
!
      IF      (xval .LE. x(1)) THEN
        yval = y(1)
      ELSE IF (xval .GE. x(n)) THEN
        yval = y(n)
      ELSE
        CALL find (i1, i2, xval, x, n)
        IF (i1 .EQ. i2) THEN
          yval = y(i1)
        ELSE
          yval = y(i1) + (y(i2)-y(i1)) * (xval-x(i1)) / (x(i2)-x(i1))
        END IF
      END IF
      RETURN
!
      END SUBROUTINE interp_bc


      SUBROUTINE bc_zone(bctime,nbctim, nz,nj, &
                         name_profile,roa,variable_edge)
!---------------------------------------------------------------------------
!     subroutine determines the edge zone values of profiles 
!     The zone extends from j =nz to j=nj
!     Note that different profiles can have different values of nz
!     and for a given profile nz may  vary with time.
!     for a given profile nz is in normalized flux (0,1.]
!     or in grid points
!     Spline the results onto a [0,1] r grid
! INPUT
! nj               rho grid size (1 ==> nj)
! nz(1:nbctim)     grid point or normalized rho values at which
!                  boundary condition is to be applied, as a function
!                  of time. 
! OUTPUT:
! variable_edge   = 0  no zonal boundary condition for this profile 
!                      was specified
!                 = 1  zonal boundary conditions given in terms 
!                      of grid points
!                 = -1 zonal boundary conditions given in terms         
!                      normalized flux
! nz(1:nbctim)     may be reset if out of range
!---------------------------------------------------------------HSJ-------
  USE nrtype,                                  ONLY : DP,I4B
  USE error_handler 
  USE io_gcnmp,                                 ONLY : nlog
      IMPLICIT NONE
      INTEGER(I4B)                                                  &
                  nj,j,nbctim,i1,i2,variable_edge,index

      REAL(DP)    roa(nj),nz(*),fun,bctime(*),fnj,diff1,diff2
      CHARACTER *(*) name_profile
      fnj = REAL(nj,DP)
      variable_edge = 0
      DO j=1,nbctim ! check if user set any of the nz values
         IF( nz(j) .LT. fnj )THEN   ! works for (0,1] or (1,nj] input
              variable_edge = 1
         ENDIF
      ENDDO
 

      IF(variable_edge .EQ. 0)RETURN  ! no boundary zone for this case


!     Here at least one nz(j) was set by the user so we have to process 
!     the data at all times.
  DO j=1,nbctim
         IF(j .EQ. 1)THEN  !first value MUST be set,we will use it to determine
                           !if input is in normalized flux or in grid point
            IF(2. .LE. nz(j).AND. nz(j) .LT. fnj)THEN
                  variable_edge = 1            !input is in grid points
            ELSE IF(0.0 .LT. nz(j).AND. nz(j)  .LT. 2.0) THEN 
                  variable_edge = -1           !input is in normalized rho
            ELSE
               WRITE(*,'(" Error in zone edge boundary input ",/, &
                     &    " for profile ",a)')name_profile
               lerrno = iomaxerr + 40
               CALL terminate (lerrno,nlog)
            ENDIF
         ELSEIF( j .GT. 1)THEN  ! input is set by user for j > 1 (ie time gt than initial time)
            IF(variable_edge .EQ. 1)THEN  ! set above Note that variable_edge is constant in time
               IF(nz(j) .LT. 2.D0 .OR. nz(j) .GT. fnj)THEN
               WRITE(*,'(" Input error for profile ",a)')name_profile
               WRITE(*,'(" At bctime point # ",i5)')j
               lerrno = iomaxerr + 40
               CALL terminate (lerrno,nlog)
               ENDIF
            ELSE  ! normalized rho
               IF(nz(j) .LE. 0.0D0 .OR. nz(j) .GT. 1.D0)THEN
                WRITE(*,'(" Input error for profile ",a)')name_profile
                WRITE(*,'(" At bctime point # ",i5)')j
               lerrno = iomaxerr + 40
               CALL terminate (lerrno,nlog)
               ENDIF
            ENDIF
         ENDIF
         IF(variable_edge .EQ. 1) THEN
               nz(j) = MIN(fnj,nz(j))      !gridpoint input
         ELSE  !variable_edge .eq. -1
               nz(j)  = MIN(1.0D0,nz(j))   !normalized flux input
         ENDIF

                 
  ENDDO


      RETURN
      END SUBROUTINE bc_zone



      SUBROUTINE get_vloop_bc(t_want,v_want)

!  -------------------------------------------------------------
!     return  the loop voltage given at time t_want in  v_want
!  -------------------------------------------------------------HSJ
      USE nrtype,                ONLY : DP,I4B
      USE ions_gcnmp,                  ONLY : nprim,nimp

      IMPLICIT NONE
      REAL(DP)  t_want, v_want
      INTEGER(I4B) k


          k = nprim+nimp+3_I4B
          CALL interp_bc (t_want, bctime, nbctim, bc(1,k), v_want)

      RETURN
      END       SUBROUTINE get_vloop_bc




      SUBROUTINE check_zero_bc(t_bc)
!
! ----------------------------------------------------------------------
! --- check to make sure boundary conditions are set properly
! --- IF there is a file associated with nlog write the error
! --- to the file, otherwise write it to ncrt
! ------------------------------------------------------------------ HSJ
!
      USE nrtype,                                 ONLY : DP,I4B
      USE solcon_gcnmp,                           ONLY : ntot,itran,time,nameu
      USE error_handler,                          ONLY : iomaxerr, lerrno,terminate
      USE io_gcnmp,                               ONLY : nlog
     
      IMPLICIT  NONE 
      INTEGER(I4B) j
      REAL(DP) t_bc
      lerrno  = 0


      DO j=1,ntot -1       ! zero bc for rotation is acceptable
        IF (ub(j) .EQ. 0.0 .AND. itran(j) .NE. 0) THEN
          WRITE (nlog,'("ERROR: boundary condition=0 for ", &
                &         " dependent variable: ", &
                &      a," at time = ",1pe14.8)')nameu(j),time
          lerrno = iomaxerr + 80_I4B
          WRITE(nlog,'("t_bc,bctime =",3(2x,1pe14.8))')t_bc,bctime
          WRITE(nlog,'("bc values are :",2(2x,1pe14.8))')bc(1,j),bc(2,j)
        END IF

      END DO
!
      IF (lerrno .GT. 0) &
        CALL terminate(lerrno,nlog)
      RETURN
!
      END SUBROUTINE check_zero_bc



      SUBROUTINE find1 (m, n, val, array, ia)
!
      USE nrtype,     ONLY:DP,I4B
      IMPLICIT  NONE
      INTEGER(i4B) n,m,ia,im,in
      REAL(DP)  is,val
!
      REAL(DP)  array (ia)
!
!  this routine finds m an n such that array(m) < val < array(n)
!  if val = array(i), m=n=i are returned
!  if val is not within the array limits, m = n=0 are returned
!  the routine works whether array is ordered largest to least
!  or least to largest
!
      is = array(ia)-array(1)
      IF (is .EQ. 0.0_I4B)  go to 12
      is = is / ABS (is)
      IF (is*(val-array(1 )))  12, 15, 10
   10 IF (is*(val-array(ia)))  30, 20, 14
   12 CONTINUE
      IF(ABS(val - array(1))  .LT. 10.*SPACING(val))THEN
         m = 1 ; n =1
      ELSE
         m = 0 ; n = 0
      ENDIF
      RETURN
   14 CONTINUE
      IF(ABS(val - array(ia))  .LT. 10.*SPACING(val))THEN
         m =ia ; n= ia
      ELSE
         m = 0 ; n = 0
      ENDIF
      RETURN
   15 m = 1
      n = 1
      RETURN
   20 m = ia
      n = ia
      RETURN
   30 IF (is .LT. 0)  go to 60
      m = 1
      n = ia
   35 im = (m+n)/2
      IF (im .EQ. m)  RETURN
      IF (val-array(im))  40, 50, 45
   40 n = im
      go to 35
   45 m = im
      go to 35
   50 m = im
      n = im
      RETURN
   60 m = ia
      n = 1
   65 in = (m+n)/2
      IF (in .EQ. n)  RETURN
      IF (val-array(in))  70, 80, 75
   70 n = in
      go to 65
   75 m = in
      go to 65
   80 m = in
      n = in
      RETURN
!
      END SUBROUTINE find1

      SUBROUTINE find (m, n, val, array, ia)
      USE nrtype,      ONLY : DP,I4B
!
      IMPLICIT  NONE
      REAL(DP) val, array(*)
      INTEGER(I4B) m, n,ia,im
!
! ----------------------------------------------------------------------
!     this routine finds m an n such that array(m) > val > array(n)
!     if val  =  array(i),m=n=i are returned
!     if val is not within the array limits, m = n=0 are returned
!     the routine presumes that array(i+1) .gt. array(i)
! ----------------------------------------------------------------------
!
      IF (val-array(1 ))  12, 15, 10
   10 IF (val-array(ia))  30, 20, 12
   12 m = 0
      n = 0
      RETURN
   15 m = 1
      n = 1
      RETURN
   20 m = ia
      n = ia
      RETURN
   30 m = 1
      n = ia
   35 im = (m+n)/2
      IF (im .EQ. m)  RETURN
      IF (val-array(im)) 40,50,45
   40 n = im
      go to 35
   45 m = im
      go to 35
   50 m = im
      n = im
      RETURN

      END SUBROUTINE find




      SUBROUTINE get_IP_from_vloop(iter_vloop,time_c)
!  ----------------------------------------------------------------------------
!     determine tot_cur from vloop value
!  -------------------------------------------------------------------HSJ------
         USE nrtype,                      ONLY : DP,I4B,I2B


         USE solcon_gcnmp,                      ONLY : time,ntot,iangrot,u,dudr, &
                                                 vloopvb

         USE dep_var,                     ONLY : etor,rbp

         USE grid_class,                  ONLY : nj,hcap,dr,r

         USE plasma_properties,           ONLY : dischg,diffuse

         USE common_constants,            ONLY : pi,twopi

         USE source_terms_gcnmp,                ONLY : etap

         USE error_handler,               ONLY : lerrno,iomaxerr,terminate

         USE io_gcnmp,                    ONLY : nlog
       
         USE curden_terms,                ONLY : tot_cur,curohm,curden

         IMPLICIT NONE
         INTEGER(I4B) iter_vloop,kfar,iodump,iostat,UNLINK
         INTEGER(I2b) get_next_io_unit
         REAL(DP) vloop_current,curohm_nj,tot_curden,vloop_obtained,    &
                  etor_njm1,ub_p,curden_nj,rbp_nj,time_c
         LOGICAL   exists, opened


         rbp_nj = rbp(nj)
         kfar = ntot - iangrot
         iter_vloop = iter_vloop + 1_I4B
         CALL get_vloop_bc(time_c,vloop_current)
         etor(nj) = vloop_current/(twopi*dischg%rmajor)       ! by definition
         etor_njm1 = etor(nj-1)
         vloop_obtained = vloop_current
         !use vloop_current/(2pi*dischg%rmajor) = Etor(nj) = eta(nj)*Hcap(nj)*<Johmic B/BT0>(j)
         ! to get curohm(nj) consistent with current temperatures,etc:
         curohm_nj = curohm(nj)
         curohm(nj) =  vloop_current/(2.*pi*dischg%rmajor*etap(nj)*hcap(nj))
         !total edge current consistent with vloop_current
         curden_nj = curden(nj)
         curden(nj) = curden(nj) -curohm_nj + curohm(nj)  !amps/m**2
         !integrate to get  tot_curden:
         CALL trapv (r, curden , hcap, nj, tot_curden)
         tot_curden = tot_curden*twopi                     !Amps
         rbp(nj) = 0.2*tot_curden 
         rbp(nj) = (hcap(nj)*r(nj) + hcap(nj-1)*r(nj-1))*0.5*diffuse%dcoef(kfar,kfar,nj-1)
         rbp(nj) = dr(nj-1)*0.5*(vloop_current/(twopi*dischg%rmajor)    &
                   + etor(nj-1) +0.0*etor_njm1 )/rbp(nj)
         rbp(nj) = -rbp(nj)/300. +  rbp(nj-1)
         rbp(nj) = 0.2*tot_curden
         rbp(nj) =0.5*(rbp(nj)  +rbp_nj)
         u(ntot-iangrot,nj) = rbp(nj)
         dudr(ntot-iangrot,nj-1) =  (u(ntot-iangrot,nj) -  u(ntot-iangrot,nj-1))/dr(nj-1)
         ub(ntot-iangrot) = rbp(nj)
         ub_p = ub(ntot-iangrot)
!
       IF( vloopvb .GT. 0)THEN
            iodump = get_next_io_unit ()
            INQUIRE (file = 'vloop_monitor.txt', iostat = iostat,            &
                                 exist = exists, opened = opened)
            IF (iostat .NE. 0) THEN         ! problem with INQUIRE
               lerrno = iomaxerr + 85
               CALL terminate(lerrno, nlog)
            END IF
!
            IF (exists .AND. .NOT. opened) THEN              ! file exists from previous case, trash it
                 iostat = UNLINK ('vloop_monitor.txt')
                 exists = .FALSE.
            END IF
            IF(.NOT. exists .AND. .NOT. opened)THEN          !create the file
                 iodump = get_next_io_unit ()
                 OPEN (unit = iodump, file = 'vloop_monitor.txt',                          &
                             status = 'NEW', iostat = iostat)
                 IF (iostat .NE. 0) THEN
                    lerrno = iomaxerr + 86
                    CALL terminate(lerrno, nlog)
                 END IF
                 CALL timestamp(iodump)
                 exists = .TRUE.
                 opened = .TRUE.
            END IF
!
            IF(opened)THEN
                   WRITE (iodump,'(" time = ", 1pe14.6)')time
                   WRITE (iodump,'(" iteration # :",i5," vloop_obtained :",1pe14.6,        &
                     &   " vloop bc :",1pe14.6)')iter_vloop,vloop_obtained,vloop_current
                   WRITE (iodump,'(" total current before correction: ",1pe14.6)') tot_cur
                   WRITE (iodump,'(" total current after correction: ",1pe14.6)') tot_curden
                   WRITE (iodump,'(" old current density (a/cm**2) at ra  :",1pe14.6 )')curden_nj
                   WRITE (iodump,'(" new current density (a/cm**2) at ra  :",1pe14.6 )')curden(nj)
                   WRITE (iodump,'(" old ohmic current density (a/cm**2) at ra  :",1pe14.6 )')curohm_nj
                   WRITE (iodump,'(" new ohmic current density (a/cm**2) at ra  :",1pe14.6 )')curohm(nj)
                   WRITE (iodump,'(" etor (V/cm) at ra-1 :",1pe14.6)')etor_njm1
                   WRITE (iodump,'(" etor (V/cm) at ra :",1pe14.6)')etor(nj)
                   WRITE (iodump,'(" old  rbp (gauss cm)  at ra :",1pe14.6)')rbp_nj
                   WRITE (iodump,'(" new  rbp (gauss cm)  at ra :",1pe14.6)')rbp(nj)
                   WRITE (iodump,'(" old  ub gauss cm)  at ra :",1pe14.6)')ub_p
                   WRITE (iodump,'(" new  ub gauss cm)  at ra :",1pe14.6)')ub(ntot-iangrot)
            ELSE
                    lerrno = iomaxerr + 86
                    CALL terminate(lerrno, nlog)
            ENDIF
         ENDIF
         tot_cur = tot_curden  

      END SUBROUTINE get_IP_from_vloop



      SUBROUTINE set_bc_values
!  -----------------------------------------------------------------------
!  -- load  usave and u at plasma edge
!  -----------------------------------------------------------------------

        USE solcon_gcnmp,              ONLY : ntot,u,usave,edge_flux,edge_flux_save
        USE grid_class,                ONLY : nj
        IMPLICIT NONE


        usave(1:ntot,1:nj) = u(1:ntot,1:nj)
        u(1:ntot,nj)       = ub(1:ntot)
        dudt0_save(:)      = dudt0(:)
        IF(.NOT. ALLOCATED(edge_flux))ALLOCATE(edge_flux(ntot),edge_flux_save(ntot))
        RETURN
      END SUBROUTINE set_bc_values




      SUBROUTINE set_bc_type
!  -----------------------------------------------------------------------
!  -- load  bc_type , assumed time independent
!  --------------------------------------------------------------HSJ------
      USE nrtype,                               ONLY : DP,I4B
      USE ions_gcnmp,                           ONLY : nion
      USE common_constants,                     ONLY : zeroc
        IMPLICIT NONE
        INTEGER(I4B) j
        bc_type(:) = -5 
        DO j = 1_I4B,nion
           IF( ABS(mult_den(j) - zeroc ) .lt. 10._DP * SPACING(0.0))THEN
              !we get to here only if mult_flux .ne. 0 . Hence only choice is
              bc_type(j) = -1  ! flux boundary condition
           ELSE
              bc_type(j) = 1 ! mixed boundary condition
              IF( ABS(mult_flux(j) - zeroc ) .lt. 10._DP * SPACING(0.0))      &
                   bc_type(j) = 0_I4B ! density boundary condition
           ENDIF
           !if bc_type(j) == 1 then a sensible value of gamma_bc(j) is
           !is expected . But we dont check for this since  a sensible
           !check for a sensible value is not well defined.
        ENDDO


        RETURN
      END SUBROUTINE set_bc_type
#endif


    SUBROUTINE convert_edge_indecies
!-----------------------------------------------------------------------
!     edge boundary values input as grid point numbers are converted
!     to equivalent normalized rho grid values:
!-----------------------------------------------------------------------
      USE nrtype,                               ONLY : DP,I4B
#ifdef GCNMP
 
      USE grid_class,                           ONLY : rho_gridn,nj
      USE ions_gcnmp,                           ONLY : nion
      IMPLICIT  NONE
      INTEGER(I4B) j 
      REAL(DP),DIMENSION(nj) :: roalcl
      INTEGER(I4b) gpt
      roalcl(:) = rho_gridn%data(:)

#elif defined NFREYA
      USE grid_class,                           ONLY : rho_gridn,nj
      USE ions_gcnmp,                           ONLY : nion
      IMPLICIT  NONE
      INTEGER(I4B) j 
      REAL(DP),DIMENSION(nj) :: roalcl
      INTEGER(I4b) gpt
      roalcl(:) = rho_gridn%data(:)
#else

      USE numbrs,                               ONLY : nion,nj
      USE mesh,                                 ONLY : roa
      IMPLICIT  NONE
      INTEGER(I4B) j 
      REAL(DP),DIMENSION(nj) :: roalcl
      INTEGER(I4b) gpt
      roalcl(:) = roa(:)

#endif


      IF(fix_edge_te_bc > 1._DP)THEN 
                    gpt = MIN(INT(fix_edge_te_bc),nj)                   
                    fix_edge_te_bc = roalcl(gpt)
      ENDIF
      IF(fix_edge_ti_bc > 1._DP)THEN 
                    gpt = MIN(INT(fix_edge_ti_bc),nj)                   
                    fix_edge_ti_bc = roalcl(gpt)
      ENDIF
      IF(fix_edge_rot_bc > 1._DP)THEN
                    gpt = MIN(INT(fix_edge_rot_bc),nj)                   
                    fix_edge_rot_bc = roalcl(gpt)
      ENDIF
      DO j=1,nion
         IF(fix_edge_ni_bc (j) > 1._DP)THEN
                    gpt = MIN(INT(fix_edge_ni_bc(j)),nj) 
                    fix_edge_ni_bc(j) = roalcl(gpt)
         ENDIF
      ENDDO

      RETURN
   END SUBROUTINE convert_edge_indecies





   INTEGER(I4B) FUNCTION common_bdry(njin)
  
      USE grid_class,                                        ONLY : nj
      USE common_constants,                                  ONLY : izero
      USE ions_gcnmp,                                        ONLY : nion
      USE solcon_gcnmp,                                      ONLY : ntot,itran
      USE error_handler,                                     ONLY : lerrno,iomaxerr,terminate
      USE io_gcnmp,                                          ONLY : nlog
      USE MPI_data,                                          ONLY : myid,master

      IMPLICIT NONE
      INTEGER(i4B) njin,njbl,jmi,jmx,k,ierr,itr,jni
      njbl = njin
    go to 10 ! 8888889999
      
      jmi =  njbl ; jmx =izero ; ierr =0
      Do k=1,nion
         if(itran(k) .ne. 0)THEN
            itr = 1
            jmi = MIN(jmi,ni_index(k))
            jmx = MAX(jmx,ni_index(k))
         ENDIF
      ENDDO
      IF(jmi .NE. jmx .AND. itr == 1)ierr       = 1
      IF(jmi .NE. 0 .AND. itr == 1)njbl = jmi
 
      IF(itran(nion+1) .NE.0)THEN
         IF(jmi .gt. 0 .AND. te_index .ne. jmi)ierr=2
         IF(itran(nion+2) .NE. 0 .AND. te_index .ne. ti_index)ierr= 2
         IF(itran(nion+4) .NE. 0 .AND. te_index .ne. rot_index)ierr= 2
         njbl = MIN(njbl,te_index)
      Endif
      IF(itran(nion+2) .NE.0)THEN
         IF(jmi .gt. 0 .AND. ti_index .ne. jmi)ierr=3
         IF(itran(nion+4) .NE. 0 .AND. ti_index .ne. rot_index)ierr= 3
         njbl = MIN(njbl,ti_index)
      Endif
      IF(itran(nion+4) .NE.0)THEN
         IF(jmi .gt. 0 .AND. rot_index .ne. jmi)ierr=4
         njbl = MIN(njbl,rot_index)
      Endif

 10 CONTINUE !888889999

      jni = njin
      Do k=1,nion
         if(itran(k) .ne. 0)THEN
            jni = MIN(jni,ni_index(k))
         ENDIF
      ENDDO
 
      IF(itran(nion+1) .NE. 0)njbl = MIN(jni,te_index)

      IF(itran(nion+2) .NE. 0)njbl = MIN(njbl,ti_index)
      IF(itran(nion+4) .NE. 0)njbl = MIN(njbl,rot_index)

      ! special case for itxj =1 allother it* =0:
      jni = izero
      DO k=1,nion+4
         IF(k == nion+3)CYCLE
         IF( itran(k) .NE. 0)jni = jni+1
      ENDDO
      IF(jni == izero)njbl = nj   ! only itxj (rbp equation) is active



      Do k=1,nion
         if(itran(k) .ne. 0)THEN
            ni_index(k) = njbl
         ENDIF
      ENDDO


      IF(itran(nion+1) .NE. 0)te_index = njbl
      IF(itran(nion+2) .NE. 0)ti_index = njbl
      IF(itran(nion+4) .NE. 0)rot_index = njbl
      ierr = 0
      IF(ierr .NE. 0)THEN
           PRINT *,'te_index,ti_index,rot_index =',te_index,ti_index,rot_index
           PRINT *,'ni_index(1:nion) = ',ni_index(:)
           lerrno = iomaxerr + 139_I4B
           CALL terminate(lerrno,nlog) 
      ENDIF

      common_bdry = njbl
      IF(myid == master) &
          WRITE(nlog,FMT='("selected common boundary is set to njb = ",i5)')common_bdry 

   END FUNCTION common_bdry




   SUBROUTINE get_nb(k,nb)
! ----------------------------------------------------------
! -- get outside rho boundary point nb  for variable type k:
! -- INPUT 
! --    k    defines the vairable type , ni(1:nion)
!            k = nion + 1  ==> Te
!            k = nion + 2  ==> Ti
!            k = nion + 3  ==> rbp
!            k = nion + 4  ==> angrot  (toroidal rotation)
! ----------------------------------------------------------

     USE grid_class,                         ONLY : nj
     USE ions_gcnmp,                         ONLY : nion

     IMPLICIT NONE
     INTEGER(I4b) k,nb

             IF(k .LE. nion)THEN    
                nb = ni_index(k)
             ELSEIF(k == nion+1)THEN
                nb = te_index
             ELSEIF(k == nion+2)THEN
                nb = ti_index
             ELSEIF( k == nion+3)THEN
                nb = nj
             ELSEIF( k == nion+4)THEN
                nb = rot_index
             ENDIF
     RETURN
   END SUBROUTINE get_nb

  END MODULE BC_values_gcnmp
