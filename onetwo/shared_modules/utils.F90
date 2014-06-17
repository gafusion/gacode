
      SUBROUTINE  col_member(x,j,jm)
!  -----------------------------------------------------------------------
!  --  x(i)= k means column  k  of the Jacobian is to be processed
!  --  determine if grid point j is associated with a column k.
!  -----------------------------------------------------HSJ---1/31/06-----
      USE nrtype,                     ONLY : I4b,dp
      USE MPI_data,                   ONLY : mpiierr
      USE solcon_gcnmp,               ONLY : grid_point
      IMPLICIT NONE

      INTEGER(I4B),INTENT(IN)::  j
      INTEGER(I4B),INTENT(INOUT) :: jm
      INTEGER(I4B),INTENT(IN),DIMENSION(:)  :: x 
      INTEGER(I4B) i ,ii

      jm = 0 
      DO  i =1,SIZE(x)
        ii = x(i)
        IF (ii .GT. 0 .AND. grid_point(ii)  == j) THEN
           jm =1
           EXIT
        ENDIF
      ENDDO

      RETURN
 
      END       SUBROUTINE col_member



      SUBROUTINE cubicextrp (f2, f3, f4, r2, r3, r4, f1, itype)
!
      USE nrtype,                                  ONLY : Dp,I4B
      IMPLICIT NONE 
!
! ----------------------------------------------------------------------
! --- subroutine gets the value of the function f at r(1) ( = 0)
! --- by extrapolation of the fitted quadratic or cubic segment:
!      f(r) = a*r**2+b              itype = 1
! or
!      f(r) = a*r**3+b*r**2+c       itype = 2
! or   the finite difference form
!      f(1) = (4/3)*f(2)-(1/3)*f(3) itype = 3
! note that in all cases df/dr = 0 at r = 0.0
! 
! --- input
!  itype               = 1 selects quadratic
!                      = 2 selects cubic
!  f2                  value of f at r2
!  f3                                r3
!  r2                  independent variable
!  r3
! if itype = 2 then  the following also have to be input
!  f4
!  r4
! --- output
! f1               value of f at r1,where it is assumed that r1 = 0
!                  with df/dr = 0.0 at r = 0.0
! ------------------------------------------------------------------ HSJ
!

      REAL(DP) f1,f2,f3,f4,r2,r3,r4,h1,h2,a1,b1,c1,detc,a,det
      INTEGER(I4B) itype

      IF      (itype .EQ. 1) THEN
          a    = (f3-f2)/(r3**2-r2**2)
          f1   = f2-a*r2**2
      ELSE IF (itype .EQ. 2) THEN
          det  = (r4**3)*(r3**2-r2**2)-(r4**2)*(r3**3-r2**3)              &
              + (r3**3)*(r2**2)-(r3**2)*(r2**3)
          detc = (r4**3)*(f2*r3**2-f3*r2**2)-(r4**2)*(f2*r3**3-f3*r2**3)  &
              +  f4*((r3**3)*(r2**2)-(r3**2)*(r2**3))
          f1   = detc/det
      ELSE             ! itype = 3
          h1   = r2
          h2   = r3 - r2
          a1   = -(2.0*h1+h2)/(h1*(h1+h2))
          b1   =  (h1+h2)/(h1*h2)
          c1   = -(h1/h2)/(h1+h2)
          f1   = -(b1*f2+c1*f3)/a1
      END IF
!
      RETURN
!
      END



      SUBROUTINE delete_file(filename)
! ----------------------------------------------------------------------
! unlink not standard ?? So use this instead - HSJ
! _____________________________________________________________________

        USE nrtype,                      ONLY : I4B
        USE io_gcnmp,                     ONLY : io_repeat
        IMPLICIT NONE
        CHARACTER(LEN = *) filename
        LOGICAL exists,opnd
        INTEGER(I4B) iounit                         


        INQUIRE(FILE=filename,EXIST = exists, OPENED = opnd, NUMBER = iounit)

        IF(exists)THEN
           IF( .NOT. opnd)THEN
               iounit = INT(io_repeat,I4B)
               OPEN(UNIT=iounit,FILE=filename)
           ENDIF
           CLOSE(UNIT=iounit,STATUS = "DELETE")
        ENDIF
        RETURN
      END SUBROUTINE delete_file





      SUBROUTINE difydx (x, y, yprime, nj)
!______________________________________________________________________

      USE nrtype,              ONLY : DP,I4B
      IMPLICIT  NONE

! ----------------------------------------------------------------------
! --- subroutine differentiates y w.r.t. x, for nj values of x
! --- and returns result in yprime.
! --- we use the basic quadratic langrangian interpolation formula
! --- for non-uniform mesh spacing in x to model y and then
! --- differentiate this model analytically. This yields central
! --- differences, accurate to order of the mesh spacing squared,
! --- appropriate forward (at j = 1) and backward (at j = nj) formulae
! --- are used.
! --- input
!  x(1..nj)                   independent variable
!  y(1..nj)                   dependent variable
!  nj                         size of x,y
! --- output
!  yprime(1..nj)              dy/dx
! ------------------------------------------------------------------ HSJ
!
      INTEGER(I4B) j,nj
      REAL(DP)  x(nj), y(nj), yprime(nj)
      REAL(DP)  dxm,dxp,dxr

!
      dxm       = x(2) - x(1)
      dxp       = x(3) - x(2)
      yprime(1) = -y(1)*(2.0*dxm+dxp)/(dxm*(dxm+dxp)) &
                 + y(2)*(dxm+dxp)/(dxm*dxp) &
                 - y(3)*(dxm/dxp)/(dxp+dxm)
!
      DO j=2,nj-1
        dxr = dxp/dxm
        yprime(j) = -y(j-1)*dxr/(dxp+dxm) &
          +y(j)*(1.0/dxm-1.0/dxp) & 
          +y(j+1)/(dxr*(dxp+dxm))
        IF (j .NE. nj-1) THEN
          dxm = dxp
          dxp = x(j+2)-x(j+1)
        END IF
      END DO
!
      yprime(nj) = y(nj-2)/(dxm+dxp) &
                 - y(nj-1)*(dxp+dxm)/(dxp*dxm) &
                 + y(nj)*(2.0*dxp+dxm)/(dxp*(dxp+dxm))
!
!     yprime(1 ) = (y(2)-y(1))/(x(2)-x(1))
!     yprime(nj) = (y(nj)-y(nj-1))/(x(nj)-x(nj-1))
      RETURN
!
      END



  SUBROUTINE extract_name(line,name)
! -----------------------------------------------------------
     USE nrtype,              ONLY : I4B
     USE ions_gcnmp,                ONLY : name_size

     IMPLICIT NONE
     CHARACTER(len=*) line,name
     INTEGER(I4B) i,j


          i  = INDEX(line,': ', BACK=.TRUE.)
          i = i+1
          j = i  + name_size 
          name  = ADJUSTL(line(i:j))

      RETURN
  END SUBROUTINE extract_name



      SUBROUTINE mescon (y, dr, nj)
        USE nrtype,                       ONLY : DP, I4B
!
      IMPLICIT  NONE
!
!     this subroutine converts a mesh-centered array to a mesh-point
!     array. linear interpolation and extrapolation are used.
!     dr(j) = r(j+1)-r(j), j=1...nj-1
!
      REAL(DP)  y(*), dr(*),x
      INTEGER(I4B) j,nj
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
      END




  SUBROUTINE   MYMPE_DECOMP1D(istart,n,nump,rank,ns,ne)
! ---------------------------------------------------'------------------------
! takes the place of mpe version on machines that dont have mpe
! Slightly more general, see below  to mimick mpe version
! istart is starting index (eq for arrays this would normally be 1 or 0)
! other input same as mpe version
! --------------------------------------------------------------------------
    USE nrtype,                       ONLY : DP,I4B
    IMPLICIT NONE
    INTEGER(I4B) istart,n,nump,rank,ns,ne,k,l,j
    INTEGER(i4b) , DIMENSION(:) :: start(nump),END(nump)
    

    k= n/nump   ! truncated # grid points per processor
    l = n-nump*k
    start(1) = istart
    DO j=2,nump
     IF(l .GT. 0)THEN
       start(j) = start(j-1) + k +1
       END(j-1) = start(j)   - 1
       l= l-1
     ELSE
       start(j) = start(j-1) + k
       END(j-1) = start(j)   - 1
     ENDIF
    ENDDO

    END(nump)   = n -1 + istart


    ns = start(rank+1)
    ne = END(rank+1)

  RETURN
  END   SUBROUTINE MYMPE_DECOMP1D




#if defined (USEMPI)
      SUBROUTINE proc_checkin(var_mon)
!     -----------------------------------------------------------
!     For debuging, each processor checks in here
!     and sends its value of var_mon to the master process.
!     It is assumed that var_mon should be the same on all
!     processors and hence the max minus the min( which should be zero)
!     of the values is  print out. 
!     -----------------------------------------------------------
     

          USE nrtype,                       ONLY : DP, I4B,I2B
          USE common_constants,             ONLY : zeroc
          USE solcon_gcnmp,                 ONLY : msg_monitr
!

          USE mpi


          USE     MPI_data,                 ONLY : mpiierr,myid,master,  &
                                                   numprocs,sizeof_dp,   &
                                                   mstr_window,mstr_active
          IMPLICIT  NONE
          REAL(DP) ,INTENT(IN)  :: var_mon  ! variable to monitor
          INTEGER(I4B) size
          REAL(DP) maxdif,mindif
          IF(numprocs ==1)RETURN




#ifdef MPI2_EXT
          !Initialization , executed only once during run:

          IF(.NOT. ALLOCATED(msg_monitr))THEN 
                ALLOCATE(msg_monitr(0:numprocs-1))
                msg_monitr(:) = zeroc
                CALL MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, sizeof_dp, mpiierr) 
             IF(myid == master)THEN
                size = numprocs*sizeof_dp
                CALL MPI_WIN_CREATE(msg_monitr,size,sizeof_dp,MPI_INFO_NULL,  &
                     MPI_COMM_WORLD,mstr_window,mpiierr)
             ELSE 
                !Slaves don't have corresponding window:
                size = 0 ! means no window is specified on slaves
                CALL MPI_WIN_CREATE(msg_monitr,size,sizeof_dp,MPI_INFO_NULL,   &
                     MPI_COMM_WORLD,mstr_window,mpiierr)
             ENDIF
             mstr_active = .TRUE.
          ENDIF

         
           CALL MPI_WIN_FENCE(0, mstr_window, mpiierr) 
          IF(myid .NE.  master )THEN
             !transfer var_mon  to  master process in slot
             !msg_monitr(myid) (msg_monitr is implicit through mstr_window) :
             CALL MPI_PUT(var_mon,1,MPI_DOUBLE_PRECISION,master,myid, &
                          1,MPI_DOUBLE_PRECISION,mstr_window,mpiierr)
          ELSE
             msg_monitr(myid)= var_mon
          ENDIF
          CALL MPI_WIN_FENCE(0, mstr_window, mpiierr) 

          IF(myid == master)THEN
              maxdif = MAXVAL(msg_monitr)
              mindif = MINVAL(msg_monitr)
              PRINT*,'msg_monitr, max-min =',maxdif-mindif
          ENDIF
#endif



         RETURN

      END SUBROUTINE   proc_checkin
#endif



      SUBROUTINE print_matrix(mat,nr,nc)
      USE nrtype,                       ONLY : DP, I4B,I2B
!
      IMPLICIT  NONE
!
      CHARACTER(LEN = 32) filename,fmtstr
      INTEGER(I2B) iounit,npl,get_next_io_unit 
      INTEGER(I4B) nr,nc,i,j,k,kmax,jr
      REAL(DP) mat(nr,nc),sum,sum_min
      
      DATA filename /'matrix_printout'/


      npl = 10
      sum_min = 1.e38
      WRITE(fmtstr,FMT='(i5,"(1pe10.2,2x)")')npl


      CALL delete_file(filename)                              !overwrite previous file
      iounit = get_next_io_unit ()
      OPEN (unit = iounit, file = filename, status = 'NEW') 

      WRITE(iounit,FMT='(i5,3x,i5)')nr,nc    
      DO j=1,nr
         WRITE(iounit,FMT='(" row number ;",i5)')j
         sum =0.0_DP
         DO k=1,nc
            IF(j.NE. k)THEN
               SUM=SUM+ABS(mat(j,k))
            ELSE
               sum =sum -ABS(mat(j,j))
            ENDIF
         ENDDO
         sum = -sum
         sum_min =MIN(sum_min,sum)
         IF( sum == sum_min) jr =j
         DO k = 1,nc,npl
            kmax = MIN(nc,k+npl-1)
            WRITE(iounit,FMT=fmtstr)(mat(j,i),i=k,kmax)  
         ENDDO
      ENDDO
      WRITE(iounit,FMT='(" diag dom ,row = ",1pe12.2,2x,i5)')sum_min,jr
      CLOSE(unit=iounit)

      RETURN
      END


      function reldif (x1, x2)
!
      USE nrtype,        ONLY : DP,I4B
      IMPLICIT NONE
      REAL(DP) del,xabs,reldif,x1,x2
      del    = ABS (x1-x2)
      xabs   = ABS (x1   )
      if (xabs .gt. 1.0e-20_DP)  del = del / xabs
      reldif = del
      return
!
      end



      SUBROUTINE trap2 (x, y, yint, npts)
!
      USE nrtype,           ONLY : DP,I4B

      IMPLICIT  NONE
!
! ----------------------------------------------------------------------
! calculates integral y*dx using trapezoidal rule
! ----------------------------------------------------------------------
!
      REAL(DP)  x(*), y(*), yint(*)
      INTEGER(I4B) i,npts
!
      yint(1) = 0.0
      DO 10 i=2,npts
   10 yint(i) = (y(i)+y(i-1))*(x(i)-x(i-1))*0.5_DP+yint(i-1)
      RETURN
!
      END

      SUBROUTINE trap3 (r, y, nj, fact, const, yint)
!---------------------------------------------------------
!      integral of {y*r*fact}*const is returned as a function
!      of r
!--------------------------------------------------------HSJ
!
      USE nrtype,             ONLY : DP,I4B
      IMPLICIT  NONE
!
      REAL(DP)  r(*), y(*), fact(*), yint(*),xint,          const
      INTEGER(I4B) j,nj
!
      yint(1) = 0.0
      DO j=2,nj
        xint = 0.5*(r(j-1)*fact(j-1)*y(j-1)+r(j)*y(j) &
               *fact(j))*(r(j)-r(j-1))
        yint(j) = yint(j-1)+xint*const
      END DO
      RETURN
!
      END

      SUBROUTINE trapv (r, y, fact, nj, xint)
        USE nrtype,                         ONLY : DP,I4B
!
      IMPLICIT  NONE
!
!     this subroutine integrates fact(r)*r*y(r) with respect to r from
!     zero to rminor.  the trapezoidal rule is used.
!
      REAL(DP)   r(*), fact(*), y(*),xint
      INTEGER(I4b) nj,j
!
      xint = 0.0_DP
      DO 10 j=2,nj
   10 xint = xint+0.5_DP*(r(j-1)*fact(j-1)*y(j-1)+r(j)*y(j)       &
           *fact(j))*(r(j)-r(j-1))
      RETURN
!
      END


      SUBROUTINE scale_length (pro, var, nj, out)
!
      USE nrtype,                ONLY : DP,I4B

      IMPLICIT NONE
!
! ----------------------------------------------------------------------
! RETURN out(j=1,..nj) = (-1/pro)(dpro/dvar) (i.e., the gradient scale
!                        length)
!
! input  pro(j),var(j) WHERE pro is the profile defined
!        as a FUNCTION of the INDEPENDENT variable variable
! -------------------------------------------------------- HSJ ---------
!
      REAL(DP)  pro(*), var(*), out(*)
      INTEGER(I4B) j, nj
!
      CALL difydx (var, pro, out, nj)
!

      DO j=1,nj
        out(j) = -out(j) / pro(j)
      END DO
      RETURN
!
      END




  FUNCTION get_next_io_unit () RESULT (next)
! -------------------------------------------------------------
! find a unit number available for i/o action
! --------------------------------------------------------------
  USE nrtype,         ONLY : I2B

  USE io_gcnmp,       ONLY : niterdb,nout,nlog,ncrt,ioplot

  USE error_handler,  ONLY : lerrno, terminate

  USE MPI_data,       ONLY : myid,mpiierr,master

#if defined (USEMPI)
      USE mpi
#endif

  IMPLICIT NONE
  INTEGER(I2B) :: next   ! the next available unit number 

  INTEGER(I2B), PARAMETER :: min_unit = 100_I2B, max_unit = 999_I2B
  INTEGER(I2B), SAVE      :: last_unit = ioplot  ! start  with unit ioplot+1
  ! (unit 6 is ncrt above)
  INTEGER(I2B)           :: count                ! number of failures
  LOGICAL            :: OPEN                     ! file status

#if defined USEMPI  
  IF(myid .NE. master)THEN
    PRINT *,'myid =',myid , ' is calling get_next_io_unit'
    PRINT *,'only master(process # ',master,' )' 
    PRINT *,' is allowed to call get_next_io_unit'
    CALL MPI_ABORT(MPI_COMM_WORLD,lerrno,mpiierr)
    CALL EXIT(1)
  ENDIF
#endif

  count = 0_I2B ; next = min_unit - 1_I2B
  IF ( last_unit > 0_I2B ) THEN ! check next in line
     next = last_unit + 1_I2B
     INQUIRE (unit=next, opened=OPEN)
     IF ( .NOT. OPEN ) last_unit = next ! found it
     RETURN
  ELSE ! loop through allowed units
     DO ! forever
        next = next + 1_I2B
        INQUIRE (unit=next, opened=OPEN)
        IF ( .NOT. OPEN ) THEN 
           last_unit = next     ! found it
           EXIT ! the unit loop
        END IF
        IF ( next == max_unit ) THEN ! attempt reset 3 times
           last_unit = 0
           count     = count + 1_I2B
           IF ( count <= 3 ) next = min_unit - 1_I2B
        END IF ! reset try
        IF ( next > max_unit ) THEN ! abort
           lerrno =6
           CALL terminate(lerrno,nlog)
        END IF ! abort
     END DO ! over unit numbers
  END IF ! last_unit
END FUNCTION get_next_io_unit


    INTEGER FUNCTION STR_LENGTH (string)
!  -------------------------------------------------------------------
!  return trimmed (no trailing spaces or tabs) length of a string ---
!  needed because LEN_TRIM f90 intrinsic function doesnt seem to work
!  on some platforms
! -------------------------------------------------------------------
    IMPLICIT NONE

    CHARACTER  last*1, string*(*)

    STR_LENGTH = LEN (string)
    IF (STR_LENGTH .NE. 0) THEN
      last = string(STR_LENGTH:STR_LENGTH)
      DO WHILE (last .EQ. ' ' .OR. last .EQ. '        ')
        STR_LENGTH = STR_LENGTH - 1
        IF (STR_LENGTH .NE. 0) THEN
          last = string(STR_LENGTH:STR_LENGTH)
        ELSE
          last = '#'
        END IF
      END DO
    END IF
    RETURN

    END





   SUBROUTINE null_pointers
!  -------------------------------------------------------------------
!  change pointer state from undefined to null
!  -------------------------------------------------------------------

  USE ions_gcnmp ,                  ONLY : namep,namei,namen,            &
                                           zeff,dzdte,zsq,z,dzdtim,      &
                                           dmassden,atomno,atw,          &
                                           tot_primary_ion_den,          &
                                           tot_thermal_ion_den,          &
                                           fi_index,z_back,ni_sc

  USE neutral_beams,                ONLY : nameb

  USE plasma_properties,            ONLY : profile,diffuse,dischg,       &
                                           mhd_dat,wpdot,pwrden,         &
                                           prtcl_src

  USE solcon_gcnmp,                 ONLY : time_hist,ssqr_time_hist,     &
                                           te0_hist,ti0_hist,rbp0_hist,  &
                                           angrot0_hist
#ifdef P_Nfreya
  USE P_nfreya_interface,           ONLY : d_fast_ion
    NULLIFY(d_fast_ion%fi_d)
#endif

    NULLIFY(namep,namei,namen,nameb,fi_index)
    NULLIFY(zeff,dzdte,zsq,z,z_back,dzdtim,ni_sc)
    NULLIFY(dmassden,atomno,atw)
    NULLIFY(tot_primary_ion_den,tot_thermal_ion_den)


    NULLIFY(profile%en,profile%flux)
    NULLIFY(profile%angrot%data,profile%etor%data)
    NULLIFY(profile%te%data,profile%ti%data,profile%ene%data)
    NULLIFY(profile%etor%data,profile%press%data,profile%pressb%data)





    NULLIFY(diffuse%xkangrot,diffuse%chiwneo)
    NULLIFY(diffuse%xchiitot)
    NULLIFY(diffuse%xdchitot,diffuse%xchietot)
    NULLIFY(diffuse%dcoef,diffuse%dcoef_glf23)
    NULLIFY(diffuse%chiinv%data,diffuse%chieinv%data,diffuse%xkineo%data)
    NULLIFY(diffuse%xkeneo%data,diffuse%xkangrot,diffuse%chiwneo)
    NULLIFY(diffuse%xchiitot,diffuse%xdchitot,diffuse%xchietot)
    NULLIFY(diffuse%xketot,diffuse%xkitot)
    NULLIFY(diffuse%dcoef,diffuse%dcoef_sav)
    NULLIFY(diffuse%dcoef,diffuse%dcoef_neo)
    NULLIFY(diffuse%dcoef,diffuse%dcoef_anal)
    NULLIFY(diffuse%vpinch_nclass)

    NULLIFY(dischg%psivolpnpsi%data)
    NULLIFY(dischg%pindentnpsi%data)
    NULLIFY(dischg%cxareanpsi%data)
    NULLIFY(dischg%triangnpsi_u%data)
    NULLIFY(dischg%triangnpsi_l%data)
    NULLIFY(dischg%rmajavnpsi%data)
    NULLIFY(dischg%rminavnj%data)
    NULLIFY(dischg%elongxnpsi%data)
    NULLIFY(dischg%sfareanpsi%data)
    NULLIFY(dischg%grho1npsi%data)
    NULLIFY(dischg%grho2npsi%data)

    NULLIFY(dischg%psivolpnj%data)
    NULLIFY(dischg%pindentnj%data)
    NULLIFY(dischg%cxareanj%data)
    NULLIFY(dischg%triangnj_u%data)
    NULLIFY(dischg%triangnj_l%data)
    NULLIFY(dischg%rmajavnj%data)
    NULLIFY(dischg%rminavnj%data)
    NULLIFY(dischg%elongxnj%data)
    NULLIFY(dischg%rmaj_geom%data)
    NULLIFY(dischg%rmin_half_width%data)
    NULLIFY(dischg%mag_bp_contr%data)
    NULLIFY(dischg%sfareanj%data)
    NULLIFY(dischg%grho1nj%data)
    NULLIFY(dischg%grho2nj%data)

    NULLIFY(dischg%rplasbdry%data)
    NULLIFY(dischg%zplasbdry%data)
    NULLIFY(dischg%rmhdgrid%data)
    NULLIFY(dischg%zmhdgrid%data)
    NULLIFY(dischg%rlimiter%data)
    NULLIFY(dischg%zlimiter%data)



    NULLIFY(mhd_dat%fcap%data,mhd_dat%gcap%data,mhd_dat%hcap%data)
    NULLIFY(mhd_dat%rcap%data,mhd_dat%r2capi%data,mhd_dat%r2cap%data)
    NULLIFY(mhd_dat%rcapi%data,mhd_dat%ravgnpsi%data,mhd_dat%ravginpsi%data)
    NULLIFY(mhd_dat%q_value%data,mhd_dat%rbp%data,mhd_dat%curden%data)
    NULLIFY(mhd_dat%curohm%data,mhd_dat%curboot%data,mhd_dat%psivalnpsi%data)
    NULLIFY(mhd_dat%fpsinpsi%data,mhd_dat%pprim%data,mhd_dat%ffprim%data)
    NULLIFY(mhd_dat%bp%data,mhd_dat%bprmaj%data,mhd_dat%btotrmaj%data)
    NULLIFY(mhd_dat%psi)


    NULLIFY(pwrden%qconde%data,pwrden%qcondi%data,pwrden%qconve%data)
    NULLIFY(pwrden%qconvi%data,pwrden%qbeame%data,pwrden%qbeami%data)
    NULLIFY(pwrden%qdelt%data,pwrden%qione%data,pwrden%qioni%data)
    NULLIFY(pwrden%qcx%data,pwrden%qe2d%data,pwrden%qi2d%data)
    NULLIFY(pwrden%qfuse%data,pwrden%qfusi%data,pwrden%qbfuse%data)
    NULLIFY(pwrden%qbfusi%data,pwrden%qmag%data,pwrden%qsawe%data)
    NULLIFY(pwrden%qsawi%data,pwrden%qrad%data,pwrden%qohm%data)
    NULLIFY(pwrden%qrfe%data,pwrden%qrfi%data,pwrden%qexch%data)
    NULLIFY(pwrden%brems_nions)



    NULLIFY(wpdot%dpedt%data)
    NULLIFY(wpdot%dpidt)
    NULLIFY(wpdot%dnedt%data)

    NULLIFY(prtcl_src%srecom)
    NULLIFY(prtcl_src%stfuse%data,prtcl_src%sbfuse%data, &
         prtcl_src%sfus%data,prtcl_src%spellet%data)

    NULLIFY(time_hist,ssqr_time_hist,te0_hist)
    NULLIFY(ti0_hist,rbp0_hist,angrot0_hist)


  RETURN
 
END SUBROUTINE null_pointers



  SUBROUTINE extrap (x1, x2, x, y1, y2, y)
!
    USE NRTYPE ,         ONLY : I4B,DP
    IMPLICIT  NONE
    REAL(DP) x1, x2, x, y1, y2, y,XX
!
! --- calculates y at x from (x1,y1) and (x2,y2) using linear extrapolation
!
        xx = (x-x1)/(x2-x1)
        y  = y1 + xx*(y2-y1)
    RETURN
!
  END SUBROUTINE extrap

      SUBROUTINE cor_lin_interp(fin,rin,nin,fout,rout,nout)
!
! ----------------------------------------------------------------------
! --- linear iterpolation with rin,rout monotonic
! --- rout(1)= rin(1),rout(nout) = rin(nin) is assumed
! ------------------------------------------------------------------ HSJ
      USE nrtype,                           ONLY   : I4b,dp

      IMPLICIT NONE 
      INTEGER(I4b)  n, jlo, j, jmid, inc
      REAL(DP),DIMENSION(nin),INTENT(IN)    :: fin,rin
      REAL(DP),DIMENSION(nout),INTENT(OUT)  :: fout,rout
      INTEGER(I4B),INTENT(IN)               :: nin,nout
      
       jlo = 0
       DO j = 2,nout-1
          CALL tableintrp (rin, nin, rout(j), jlo)
          CALL linear (fin(jlo), fin(jlo+1), fout(j),rin(jlo),  &
               rin(jlo+1),rout(j))
       ENDDO
       fout(1)    = fin(1)
       fout(nout) = fin(nin)
    
       RETURN
      END SUBROUTINE cor_lin_interp

      SUBROUTINE linear (fun1, fun2, fun, x1, x2, x)
!
      USE nrtype,           ONLY : DP,I4B
      IMPLICIT  NONE 
      REAL(DP) fun1, fun2, fun, x1, x2, x,dx
!
! ----------------------------------------------------------------------
! linear interpolation to determine fun (x)
! ----------------------------------------------------------------------
!
      dx = x2-x1
      IF(dx .GT. 0.0_DP)THEN 
         fun = ((x-x1)*fun2+(x2-x)*fun1)/dx
      ELSE
         fun = fun1
      ENDIF
      RETURN
!
      END SUBROUTINE linear

      SUBROUTINE tableintrp (xx, n, x, jlo)
!
      USE nrtype,                          ONLY : DP,I4B
      IMPLICIT NONE
!
! ----------------------------------------------------------------------
! --- correlated table search routine. use jlo from previous call to
! --- get jlo for current value of x. if jlo from previous call is
! --- no good(jlo=0 or jlo=n) then do a binary search.
! --- this routine returns jlo so that
!   xx(jlo) .le. x .le. xx(jlo+1) if jlo = 1
!   xx(jlo) .lt. x .le. xx(jlo+1) if jlo =2,3...n-1
!   it is assumed that xx(j),j = 1,2..n is monotonically
!   increasing; this is NOT checked for.
!   this is a modified version of the Numerical Recipes subroutine HUNT
! ------------------------------------------------------------------ HSJ
!
      INTEGER(I4b)  n, jlo, jhi, jmid, inc
      REAL(DP)   x, xx(n)
!
      IF      (x .LT. xx(1)) THEN
        jlo = 0                     ! indicates out of range below xx(1)
      ELSE IF (x .LE. xx(2)) THEN
        jlo = 1                     ! xx(1) .le. x .le. xx(2)
      ELSE IF (x .LE. xx(n)) THEN   ! xx(2) .lt. x .le. xx(n)
!
!       check if jlo from previous call is usable:
!
        IF (jlo .LE. 0 .OR. jlo .GE. n) THEN
          jlo = 2
          jhi = n
          go to 15           ! no correlation go directly to bisection
        ELSE                 ! 1 .le. jlo .lt. n
!
!         bracket x, then use bisection:
!         start with  jlo  from previous call:
!
          inc = 1
          IF (x .GT. xx(jlo)) THEN    ! search up
    4       jhi = jlo + inc
            IF (jhi .GE. n) THEN
              jhi = n
            ELSE IF (x .GT. xx(jhi)) THEN
              inc = inc + inc
              jlo = jhi
              go to 4
            END IF
          ELSE
    5       jhi = jlo
            jlo = jlo-inc
            IF (jlo .LE. 1) THEN
              jlo = 1
            ELSE IF (x .LE. xx(jlo)) THEN
              inc = inc + inc
              go to 5
            END IF
          END IF
        END IF
!
!         bisection:
!
   10     IF (jhi-jlo .EQ. 1)  RETURN
   15     jmid = (jhi + jlo) / 2
          IF (x .GT. xx(jmid)) THEN
            jlo = jmid
          ELSE           ! x .le. xx(jmid)
            jhi = jmid
          END IF
          go to 10
!
      ELSE               ! x .gt. xx(n)
        jlo = n
      END IF
!
      RETURN
!
      END



      SUBROUTINE to_upper_case(string)
! -------------------------------------------------------------------
!
      USE nrtype,            ONLY : I4B
      IMPLICIT NONE
      INTEGER(I4B) l
      CHARACTER*(*), INTENT (INOUT) :: string
      DO  l =1,LEN_TRIM(string),1
         IF(LGE (string(l:l),'a') .AND. &
              LLE(string(l:l),'z'))     &
              string(l:l) = ACHAR(IACHAR(string(l:l)) - 32)
      ENDDO
      RETURN
      END SUBROUTINE to_upper_case


      SUBROUTINE to_lower_case(string)
! -------------------------------------------------------------------
!
      USE nrtype,            ONLY : I4B
      IMPLICIT NONE
      INTEGER(I4B) l
      CHARACTER*(*), INTENT (INOUT) :: string
      DO  l =1,LEN_TRIM(string),1
         IF(LGE (string(l:l),'A') .AND. &
              LLE(string(l:l),'Z'))     &
              string(l:l) = ACHAR(IACHAR(string(l:l)) + 32)
      ENDDO
      RETURN
      END SUBROUTINE to_lower_case


integer function getfreeunit(unit)
  ! This is a simple function to search for an available unit.
  ! LUN_MIN and LUN_MAX define the range of possible LUNs to check.
  ! The UNIT value is returned by the function, and also by the optional
  ! argument. This allows the function to be used directly in an OPEN
  ! statement, and optionally save the result in a local variable.
  ! If no units are available, -1 is returned.
  integer, intent(out), optional :: unit
  ! local
  integer, parameter :: LUN_MIN=10, LUN_MAX=1000
  logical :: opened
  integer :: lun
  ! begin
  getfreeunit=-1
  do lun=LUN_MIN,LUN_MAX
    inquire(unit=lun,opened=opened)
    if (.not. opened) then
      getfreeunit=lun
      exit
    end if
  end do
  if (present(unit)) unit=getfreeunit
end function getfreeunit
