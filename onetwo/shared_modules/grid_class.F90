  MODULE grid_class
    USE Vector_class 
    USE nrtype,              ONLY : I4B,DP
    USE plasma_properties,   ONLY : mhd_dat,dischg
    USE common_constants,    ONLY : pisq,zeroc
    USE error_handler,       ONLY : lerrno,terminate,iomaxerr,dbg_print
    USE MPI_data,            ONLY : myid,master 
!!    USE BC_values_gcnmp,     ONLY : robin_bc,nb_robin,nb_robin_nj cant use because circular ref
#ifdef GCNMP
    USE io_gcnmp,            ONLY : ncrt,nlog 
#elif defined  NFREYA
    USE io_gcnmp,            ONLY : ncrt,nlog 
#else
    USE  io,                 ONLY : ncrt,nlog => nlog_gcnmp
#endif
    IMPLICIT NONE
    !localize definition of grids to this routine:
    TYPE(Vector) rho_gridn
    TYPE(Vector) psi_gridn
    TYPE(Vector) rho_grid
    TYPE(Vector) psir_grid
    TYPE(Vector) rho_mhd_gridnpsi
    TYPE(Vector) R_chord,psi_chord,rho_chord,ene_chord,te_chord,ti_chord

    INTEGER(I4B) nj,nj_out,nj_save,nj_max,nj_start,njm1,njb,set_cap1,npsi, &
                 use_compact_schemes,j_ped_index
    LOGICAL allow_regrid,uniform_mesh
!    r,rold,r_input are rho grids in m. r_input is the statefile input rho grid
!    it will be used to write the new statfile on output.
!    initially r = rold = r_input. If the grid is changed then rold will hold the
!    previous version of the grid and r will hold the current grid. r_input is
!    never changed.
    REAL(DP),ALLOCATABLE, DIMENSION(:) :: ra,r,drr,dr,rrp,rrm,roa,rold,r_input
    REAL(DP),ALLOCATABLE, DIMENSION(:) :: fcap,gcap,hcap,eps,rcap,r2capi,r2cap
    REAL(DP),ALLOCATABLE, DIMENSION(:) :: rcapi,hcapra
    REAL(DP),ALLOCATABLE, DIMENSION(:) :: xhm2,xhm20,xi11,xi110,xi33,xi330
    REAL(DP),ALLOCATABLE, DIMENSION(:) :: xips,xips0
    REAL(DP),ALLOCATABLE, DIMENSION(:) :: rminor_r, ravg_r
    REAL(DP),ALLOCATABLE, DIMENSION(:) :: rold_zc, rnew_zc 
    REAL(DP) volfac                     ! 4 pisq*R0
    REAL(DP) frac_core_points
    LOGICAL  fit_grid
!   eps(i) horizontal inverse aspect ratio = (rmax-rmin)/(rmax+rmin)

  CONTAINS
    SUBROUTINE meshgen
      INTEGER(I4B) j
      REAL(DP)  csum

      CALL mesh_allocate

        fcap(1:nj)     =  get_values(mhd_dat%fcap)
        gcap(1:nj)     =  get_values(mhd_dat%gcap)
        hcap(1:nj)     =  get_values(mhd_dat%hcap)
        r(1:nj)        =  get_values(rho_grid)
        r_input(:)     =  r(:)
        rcap(1:nj)     =  get_values(mhd_dat%rcap)
        rcapi(1:nj)    =  get_values(mhd_dat%rcapi)
        r2capi(1:nj)   =  get_values(mhd_dat%r2capi)
        r2cap(1:nj)    =  get_values(mhd_dat%r2cap)
        rminor_r(1:nj) =  get_values(dischg%rminavnj)
        ravg_r(1:nj)   =  get_values(dischg%rmajavnj)
        roa(1:nj)      =  r(1:nj)/r(nj)
        roa(1)         =  0.0_DP  
        roa(nj)        =  1.0_DP   
        csum = zeroc
        DO j= 1,nj-1
           ra(j) = 0.5*(r(j+1)+r(j))
           dr(j) = r(j+1) - r(j)
           hcapra(j) = 0.5_DP*(hcap(j+1)+hcap(j))*ra(j)
           IF(j .GT. 1)csum = csum +ABS(dr(j)-dr(j-1))
        ENDDO
 
        uniform_mesh = .TRUE.
        IF(csum > dr(1))uniform_mesh = .FALSE.


        drr(1) = 2.0_DP / dr(1)
        rrp(1) = 2.0_DP / dr(1)
        rrm(1) = zeroc
        DO  j=2,nj
          rrm(j)  = (ra(j-1)/r(j))/dr(j-1)
          IF (j .NE. nj) THEN
            drr(j)  = 2.0_DP/(dr(j-1)+dr(j))
            rrp(j)  = (ra(j)/r(j))/dr(j)
          END IF
        END DO
        drr(nj) = 2.0_DP/ dr(nj-1)
        volfac  = 4._DP*pisq*dischg%rmajor
        rold(:)  = r(:)
   END SUBROUTINE meshgen
   
   SUBROUTINE psi_contours
!---------------------------------------------------------------
! -- generate psi(r,z) contours
!--------------------------------------------------------HSJ----
     

       RETURN
   END    SUBROUTINE psi_contours

   SUBROUTINE reset_cap_parms(set_cap1)
        INTEGER(I4B) set_cap1
        mhd_dat%set_cap1 = set_cap1 
        IF(mhd_dat%set_cap1 == 1 ) THEN
             hcap(1:nj) = 1.0_DP
             fcap(1:nj) = 1.0_DP   
             gcap(1:nj) = 1.0_DP
             mhd_dat%fcap%data(:) = fcap(:)
             mhd_dat%hcap%data(:) = hcap(:)
             mhd_dat%gcap%data(:) = gcap(:)
        ENDIF


     END SUBROUTINE reset_cap_parms

     SUBROUTINE regrid(njn,njprev,factor)
!   -----------------------------------------------------------------------
!   regrid all quantities that are defined on grid njprev onto a new grid
!   of size njn.which is either half(factor =0.5)  or double (factor =2)
!   the current grid size.
!   Afterwards set njprev = njn. Note that rho grid may be non uniform
!   which dictates the approach taken below.
!   _______________________________________________________________________


        USE common_constants,     ONLY : zeroc

          INTEGER(i4B) njn,j,njprev         
          REAL(DP)factor,dru
  
          njb = 0 ! common boundary location index for flux solver
                  ! 0  calculation of njb in other routines
 
          IF(SIZE(rold) .NE.  SIZE(r))THEN
             DEALLOCATE(rold)
             ALLOCATE(rold(SIZE(r)))
             IF(SIZE(r) .NE. njprev) THEN
                WRITE(ncrt,FMT = '("single_density_simulation error")')
                lerrno = iomaxerr+108
                CALL terminate(lerrno,nlog)
             ENDIF
          ENDIF
          rold(1:njprev) = r(1:njprev) 

          DEALLOCATE(r)
          ALLOCATE(r(njn))  
          r(1) = rold(1)
          r(njn) = rold(njprev)
          IF(zeroc .LT. factor .AND. factor .LT. 1.0_DP)THEN
             DO j=2,njn-1
                r(j) = rold(2*j-1)
             ENDDO
          ELSEIF(factor .GT. 1._DP)THEN
             DO j=2,njn-1
                r(j) =  (rold(j/2) +rold(j/2+1))*0.5_DP
             ENDDO
          ELSE          ! factor < 0.0 (EG njprev not of 2^n+1 form)
                        ! set up a uniform grid even if input
                        ! grid is not uniform
             dru = (r(njn)-r(1))/(njn-1)
             DO j = 2,njn-1
                r(j) = r(j-1)+dru
             ENDDO
           ENDIF

          DEALLOCATE(roa)
          ALLOCATE(roa(njn))
          roa(1:njn-1)=r(1:njn-1)/r(njn)
          roa(njn) = 1.0_DP

          DEALLOCATE(ra,dr)
          ALLOCATE(ra(njn-1),dr(njn-1))
        DO j= 1,njn-1
           ra(j) = 0.5*(r(j+1)+r(j))
           dr(j) = r(j+1) - r(j)
        ENDDO

        DEALLOCATE(drr,rrp,rrm)
        ALLOCATE(drr(njn),rrm(njn),rrp(njn))
        drr(1)   = 2.0_DP / dr(1)
        rrp(1)   = 2.0_DP / dr(1)
        rrm(1)   = zeroc   ! not used
        rrm(njn) = zeroc ! not used
        rrp(njn) = zeroc ! not used

        DO  j=2,njn
          rrm(j)  = (ra(j-1)/r(j))/dr(j-1)
          IF (j .NE. njn) THEN
            drr(j)  = 2.0_DP/(dr(j-1)+dr(j))
            rrp(j)  = (ra(j)/r(j))/dr(j)
          END IF
        END DO
        drr(njn) = 2.0_DP / dr(njn-1)
 
   END SUBROUTINE regrid


   SUBROUTINE mhd2tport
! ---------------------------------------------------------------------------
! --- quantities of size npsi (mhd_dat%npsi) are on the
! --- mhd grid and are numbered from the edge to the magnetic axis
! --- the transport equations  require these profiles on the
! --- transport grid, of size nj, numbered form the axis to the
! --- plasma edge. This module accomplishes that task.
! --- Smoothing of some mhd input quantities is done if smooth_mhd_input= .TRUE.
! ---------------------------------------------------------------------------
     USE tension_spline,             ONLY : tspline90,t716_TSPSI,     &
                                            t716_TSVAL1,t716_TSPSS
     USE solcon_gcnmp,               ONLY : smooth_mhd_input

     LOGICAL periodic,uniform
     INTEGER(I4B) ncd,iendc,iflag
     REAL(DP) bpar(2)
     REAL(DP), ALLOCATABLE,DIMENSION(:)   :: psival,rev_psival,work,  &
                                             rev_generic,psir_locl

     npsi = mhd_dat%npsi
 
     ALLOCATE( psival(npsi),rev_psival(npsi),work(nj),                &
               rev_generic(npsi),psir_locl(nj))

             periodic       = .FALSE. !interpolant is not periodic
             uniform        = .TRUE.  ! tension is uniform
             ncd            = 2       ! # continuous derivatives
             iendc          = 0       ! natural spline bc
            !iendc          = 1       ! give deriv bc at end points
             bpar(1)        = 0.0     ! used only if iendc =1,2 
             bpar(2)        = 0.0     ! used only if iendc =1,2 
            ! NOTE -- TENSION PARAMETER IS OPTIONAL ARGUMENT IN 
            ! t716_TSVAL1 :
             rev_psival(:)  = get_values(mhd_dat%psivalnpsi)
             psir_locl(:)   = get_values(psir_grid)
             CALL reverse(npsi,rev_psival)
             rev_generic(:)    = get_values(dischg%rmajavnpsi)
             CALL reverse(npsi,rev_generic)
             CALL  t716_TSPSI (npsi,rev_psival,rev_generic,ncd,iendc, &
                               periodic,uniform,bpar) 
            ! IFLAG = 0 if values of function are to be computed.
            ! IFLAG = 1 if first derivative values are to be computed.
            ! IFLAG = 2 if second derivative values are to be computed.

             iflag = 0
             CALL  t716_TSVAL1 (npsi,rev_psival,rev_generic,IFLAG,nj,    &
                                psir_locl,work) 

              IF(smooth_mhd_input)      CALL average7_1d(work,nj)
              dischg%rmajavnj        =  new_Vector(nj,work)


             rev_generic(:)    = get_values(dischg%triangnpsi_u)
             CALL reverse(npsi,rev_generic)
             CALL  t716_TSPSI (npsi,rev_psival,rev_generic,ncd,iendc, &
                               periodic,uniform,bpar)
             CALL  t716_TSVAL1 (npsi,rev_psival,rev_generic,IFLAG,nj,    &
                                psir_locl,work)
              IF(smooth_mhd_input)      CALL average7_1d(work,nj)
             dischg%triangnj_u   = new_Vector(nj,work)



             rev_generic(:)    = get_values( dischg%triangnpsi_l)
             CALL reverse(npsi,rev_generic)
             CALL  t716_TSPSI (npsi,rev_psival,rev_generic,ncd,iendc,       &
                               periodic,uniform,bpar)
             CALL  t716_TSVAL1 (npsi,rev_psival,rev_generic,IFLAG,nj,    &
                                psir_locl,work)
              IF(smooth_mhd_input)      CALL average7_1d(work,nj)
             dischg%triangnj_l  = new_Vector(nj,work)


 
             rev_generic(:)    = get_values(dischg%rminavnpsi)
             CALL reverse(npsi,rev_generic)
             CALL  t716_TSPSI (npsi,rev_psival,rev_generic,ncd,iendc, &
                               periodic,uniform,bpar)
             CALL  t716_TSVAL1 (npsi,rev_psival,rev_generic,IFLAG,nj,    &
                                psir_locl,work)
              IF(smooth_mhd_input)      CALL average7_1d(work,nj)
             dischg%rminavnj =   new_Vector(nj,work)


             rev_generic(:)    = get_values(dischg%psivolpnpsi)
             CALL reverse(npsi,rev_generic)
             CALL  t716_TSPSI (npsi,rev_psival,rev_generic,ncd,iendc, &
                               periodic,uniform,bpar)
             CALL  t716_TSVAL1 (npsi,rev_psival,rev_generic,IFLAG,nj,    &
                                psir_locl,work)             
             dischg%psivolpnj      =   new_Vector(nj,work)


             rev_generic(:)    = get_values(dischg%elongxnpsi)
             CALL reverse(npsi,rev_generic)
             CALL  t716_TSPSI (npsi,rev_psival,rev_generic,ncd,iendc, &
                               periodic,uniform,bpar)
             CALL  t716_TSVAL1 (npsi,rev_psival,rev_generic,IFLAG,nj,    &
                                psir_locl,work)   
             IF(smooth_mhd_input)      CALL average7_1d(work,nj)
             dischg%elongxnj        =   new_Vector(nj,work)

                                      
             rev_generic(:)    = get_values(dischg%pindentnpsi)
             CALL reverse(npsi,rev_generic)
             CALL  t716_TSPSI (npsi,rev_psival,rev_generic,ncd,iendc, &
                               periodic,uniform,bpar)
             CALL  t716_TSVAL1 (npsi,rev_psival,rev_generic,IFLAG,nj,    &
                                psir_locl,work)      
              IF(smooth_mhd_input)      CALL average7_1d(work,nj)       
             dischg%pindentnj       =   new_Vector(nj,work)



             rev_generic(:)    = get_values(dischg%sfareanpsi)                                                  
             CALL reverse(npsi,rev_generic) 
             CALL  t716_TSPSI (npsi,rev_psival,rev_generic,ncd,iendc, &
                               periodic,uniform,bpar)
             CALL  t716_TSVAL1 (npsi,rev_psival,rev_generic,IFLAG,nj,    &
                                psir_locl,work)             
             dischg%sfareanj       =   new_Vector(nj,work)




             rev_generic(:)    = get_values(dischg%cxareanpsi)
             CALL reverse(npsi,rev_generic)
             CALL  t716_TSPSI (npsi,rev_psival,rev_generic,ncd,iendc, &
                               periodic,uniform,bpar)
             CALL  t716_TSVAL1 (npsi,rev_psival,rev_generic,IFLAG,nj,    &
                                psir_locl,work)             
             dischg%cxareanj       =   new_Vector(nj,work)
 
           

             rev_generic(:)    = get_values(dischg%grho1npsi)
             CALL reverse(npsi,rev_generic)
             CALL  t716_TSPSI (npsi,rev_psival,rev_generic ,ncd,iendc, &
                               periodic,uniform,bpar)
             CALL  t716_TSVAL1 (npsi,rev_psival,rev_generic,IFLAG,nj,    &
                                psir_locl,work)             
             dischg%grho1nj      =   new_Vector(nj,work)


             rev_generic(:)    = get_values(dischg%grho2npsi)                     
             CALL reverse(npsi,rev_generic)
             CALL  t716_TSPSI (npsi,rev_psival,rev_generic,ncd,iendc, &
                               periodic,uniform,bpar)
             CALL  t716_TSVAL1 (npsi,rev_psival,rev_generic,IFLAG,nj,    &
                                psir_locl,work)             
             dischg%grho2nj       =   new_Vector(nj,work)


             DEALLOCATE(psival)

             DEALLOCATE(rev_psival)

             DEALLOCATE(rev_generic)
 
             DEALLOCATE(psir_locl)

             DEALLOCATE(  work)

             RETURN

   END  SUBROUTINE mhd2tport


   SUBROUTINE grid_interface(initial,grid_dim)
     !---------------------------------------------------------------------------
     ! -- Entry point for change grid size and or spacing calculations.
     !---------------------------------------------------------------------------
     REAL(DP) factor
     INTEGER(I4B),INTENT(IN),OPTIONAL :: grid_dim ! pgf90 proglem with optional ??
     INTEGER(I4B) k,i,n
     LOGICAL initial,skip_vectors

  
     skip_vectors = .FALSE.
     nj_save = nj


     IF(initial)factor = -1._DP

     IF(PRESENT(grid_dim)  .AND.  grid_dim .NE. 0)THEN ! pgf90 forces this
        nj = grid_dim
        IF(nj == 2**(nj_save+1)+1)THEN            ! output grid is double current grid
         do i =1,nj_save+1
 
         enddo

           factor = 1.5_DP
        ELSEIF(nj == 2**(nj_save-1)+1)THEN        ! output grid is half  current grid
           factor = 0.5_DP 
        ELSE
           factor = -1.0_DP                       ! output grid not simply related to current grid
        ENDIF
        skip_vectors = .TRUE.                     ! user defined typedef vectors are  skipped in
     ELSE                                         ! refactor_profiles
        nj = rgrid_size(factor)
     ENDIF
!     factor = -1._DP ! disable the above
!     nj = rgrid_size(factor)
     IF(initial .AND. allow_regrid)THEN
         nj = nj_start
     ELSEIF(allow_regrid)THEN
        nj = grid_dim
     ENDIF
     njm1 = nj - 1

        IF(dbg_print .AND. myid == master)  &
             WRITE(nlog,FMT='("CALLING: regrid from grid_interface",i5,x,i5)')nj,nj_save
        IF(dbg_print .AND. myid == master)  &
             WRITE(ncrt,FMT='("CALLING: regrid from grid_interface",i5,x,i5)')nj,nj_save

     IF(nj .NE. nj_save)CALL regrid(nj,nj_save,factor)

        IF(dbg_print .AND. myid == master)  &
             WRITE(nlog,FMT='("DONE: regrid from grid_interface",i5,x,i5)')nj,nj_save
        IF(dbg_print .AND. myid == master)  &
             WRITE(ncrt,FMT='("DONE: regrid from grid_interface",i5,x,i5)')nj,nj_save


#ifdef GCNMP
     IF(nj .NE. nj_save)THEN 
        IF(dbg_print .AND. myid == master)  &
             WRITE(nlog,FMT='("CALLING: refactor_profiles from grid_interface",i5,x,i5)')nj,nj_save
        IF(dbg_print .AND. myid == master)  &
             WRITE(ncrt,FMT='("CALLING: refactor_profiles from grid_interface",i5,x,i5)')nj,nj_save

        CALL refactor_profiles(nj,nj_save,factor,skip_vectors)

!!        IF(robin_bc) CALL change_nb_robin(nj,nj_save) ! temp disable 12/12/13 HSJ 8888888899999
        IF(dbg_print .AND. myid == master)  &
             WRITE(nlog,FMT='("DONE: refactor_profiles from grid_interface",i5,x,i5)')nj,nj_save
        IF(dbg_print .AND. myid == master)  &
             WRITE(ncrt,FMT='("DONE: refactor_profiles from grid_interface",i5,x,i5)')nj,nj_save
     ENDIF

#endif
     RETURN

   END  SUBROUTINE grid_interface




   SUBROUTINE reverse (n,x)
!
            REAL(DP) work(n), x(n)
            INTEGER(I4B) k,i,n
!
            work(:) = x(:)

            k = 0
            DO i=n,1,-1
               k    = k + 1
               x(k) = work(i)
            END DO

           RETURN
!
      END  SUBROUTINE reverse

     INTEGER(I4B) FUNCTION  rgrid_size(factor)
!   -----------------------------------------------------------------------
!   Define a new grid size, 3 options:
!   a) either half (0 < factor < 1.0)  
!      or 
!   b) double (factor > 1.0)
!       or
!   c) the nearest 2^n+1 value to the given input nj_start (factor < 0.0)
!   options a and b are valid only if option c was called first
!   Note that rho grid may be non uniform
!   which dictates the approach taken below.
!   _______________________________________________________________________
          USE common_constants,                 ONLY : zeroc
          INTEGER(i4B) njn,j,njm
          REAL(DP)factor
 
          IF(zeroc  .LT. factor .AND. factor .LT. 1.0_DP)THEN
             njn = nj/2+1
          ELSE IF(factor > 1.0)THEN
             njn = 2*nj-1
          ELSE
             ! find n nearest to input nj such that nj <= 2^n + 1
             DO j =1,10000 !HUGE(1) does not work consisently ??
                njn = j
                IF(2**njn +1  .GE. nj_start)EXIT
             ENDDO
             DO j =1,10000
                njm = j
                IF(2**njm +1  .GE. nj_max)EXIT
             ENDDO
          ENDIF
          njn = MIN(njn,njm )
          rgrid_size = 2**njn + 1 
          rgrid_size = nj_start             ! disable above logic 
 
   END FUNCTION rgrid_size

   SUBROUTINE mesh_allocate
!----------------------------------------------------------------------
! --
!----------------------------------------------------------------------
      IF(.NOT. ALLOCATED(r)) ALLOCATE(r(nj))
      IF(.NOT. ALLOCATED(rold)) ALLOCATE(rold(nj))
      IF(.NOT. ALLOCATED(dr)) ALLOCATE(dr(nj-1))
      IF(.NOT. ALLOCATED(roa)) ALLOCATE(roa(nj))
      IF(.NOT. ALLOCATED(ra))ALLOCATE(ra(nj-1))
      IF(.NOT. ALLOCATED(drr))ALLOCATE(drr(nj))
      IF(.NOT. ALLOCATED(rrp))ALLOCATE(rrp(nj))
      IF(.NOT. ALLOCATED(rrm))ALLOCATE(rrm(nj))
      IF(.NOT. ALLOCATED(fcap))ALLOCATE(fcap(nj))
      IF(.NOT. ALLOCATED(gcap))ALLOCATE(gcap(nj))
      IF(.NOT. ALLOCATED(hcap))ALLOCATE(hcap(nj))
      IF(.NOT. ALLOCATED(hcapra))ALLOCATE(hcapra(nj-1))
      IF(.NOT. ALLOCATED(rcap))ALLOCATE(rcap(nj))
      IF(.NOT. ALLOCATED(rcapi))ALLOCATE(rcapi(nj))
      IF(.NOT. ALLOCATED(r2capi))ALLOCATE(r2capi(nj))      
      IF(.NOT. ALLOCATED(r2cap))ALLOCATE(r2cap(nj)) 
      IF(.NOT. ALLOCATED(ravg_r))ALLOCATE(ravg_r(nj)) 
      IF(.NOT. ALLOCATED(rminor_r))ALLOCATE(rminor_r(nj)) 
      IF(.NOT. ALLOCATED(r_input))ALLOCATE(r_input(nj)) 

   END SUBROUTINE mesh_allocate


   SUBROUTINE mescon (y, dr, nj)
!----------------------------------------------------------------------
!     this SUBROUTINE converts a mesh-centered array to a mesh-point
!     array. linear interpolation and extrapolation are used.
!     dr(j) = r(j+1)-r(j), j=1...nj-1
!----------------------------------------------------------------------
      USE nrtype,          ONLY : I4B,DP
      IMPLICIT  NONE

      INTEGER j,nj 
      REAL(DP) y(nj),dr(nj),x
!
      x      = (dr(nj-2)+2.0*dr(nj-1)) / (dr(nj-2)+dr(nj-1))
      y(nj)  = y(nj-2) + x*(y(nj-1)-y(nj-2))
      DO j=nj-1,2,-1
        x    = dr(j-1) / (dr(j-1)+dr(j))
        y(j) = y(j-1) + x*(y(j)-y(j-1))
      END DO
      y(1)   = y(1) - (y(2)-y(1))

      RETURN
!
   END SUBROUTINE mescon

      SUBROUTINE average7_1d(f,nk)
!*****************************************************
!
! performs a two pass seven point average
!
!*****************************************************
      USE nrtype,                     ONLY : I4b,dp
      IMPLICIT NONE
!
      INTEGER nk,k,i,m
      REAL(DP) f(nk),g(nk)
!
! check if grid is too small
      IF(nk.LT.7)RETURN
!
        g(1)=f(1)
        g(nk)=f(nk)
        k=2
        g(k)=(f(k)+f(k+1)+f(k-1))/3.D0
        k=3
        g(k)=(f(k)+f(k+1)+f(k-1)+f(k+2)+f(k-2))/5.D0
        k=nk-1
        g(k)=(f(k)+f(k+1)+f(k-1))/3.D0
        k=nk-2
        g(k)=(f(k)+f(k+1)+f(k-1)+f(k+2)+f(k-2))/5.D0
        DO k=4,nk-3
          g(k)=(f(k)+f(k-1)+f(k+1)+f(k-2)+f(k+2) &
                + f(k+3)+f(k-3))/7.D0
        ENDDO
        k=2
        f(k)=(g(k)+g(k+1)+g(k-1))/3.D0
        k=3
        f(k)= (g(k)+g(k+1)+g(k-1)+g(k+2)+g(k-2))/5.D0
        k=nk-1
        f(k)=(g(k)+g(k+1)+g(k-1))/3.D0
        k=nk-2
        f(k)= (g(k)+g(k+1)+g(k-1)+g(k+2)+g(k-2))/5.D0
        DO k=4,nk-3
          f(k)= (g(k)+g(k+1)+g(k-1)+g(k+2)+g(k-2) &
                  + g(k+3)+g(k-3))/7.D0
        ENDDO
      RETURN
      END  SUBROUTINE average7_1d

!      SUBROUTINE change_nb_robin(nj,nj_save)
      !-----------------------------------------------------------------------
      !  if robin_bc is in effect the currently nb_robin
      ! is based on grid size nj_save. But we now need it on grid_size nj
      ! so change it here
      !-----------------------------------------------------------------------
       
!        IF(size(rold) == nj_save)THEN
!        rho_nb_save = rold(nb_robin)
!        DO j=1,nj
!           IF(rho(j) .lt. ?rhold(nb_robin)
!               nb_robin_nj = j ??
!        ENDO
!        ELSE

!        ENDIF
        
!      END  SUBROUTINE change_nb_robin

  END MODULE grid_class
