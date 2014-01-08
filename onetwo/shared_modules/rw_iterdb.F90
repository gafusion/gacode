 SUBROUTINE iter_dbase_txt 

  !
  ! ----------------------------------------------------------------------
  ! ---  Extended iterdb file (statefile for Onetwo) collects and writes out
  ! ---  or reads in  data for ITER database. This module is used by Onetwo
  ! ---  and GCNMP at present. This routine reads/writes text files.
  ! ---  See iter_dbase_nc for netcdf files
  ! ----------------------------------------------------------------------
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
  !                      iterdb_outpt = iterdb_file_name // time 
  !     iterdsc          switch used to output description of parameters
  !     niterdb          Fortran I/O unit number for input/output
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
  !                      to the psi grid in meters
  !     qpsi(j)          the corresponding q values
  !     psivolpnpsi(j)       volume, meter**3, enclosed by flux surface
  !     grho1(j)         < ABS (grad rho)>
  !     grho2(j)         <(grad rho)**2>
  !     rmajavnpsi(j)    j=1,2,..npsi average major radius
  !                      at elevation of magnetic axis
  !     rminavnpsi(j)    same for minor radius
  !     psival(j)        j=1,2..npsi
  !
  !     triangnpsi(j)    (upper)triangularity
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
  !                      beam fd_beam set fd_beam = 1.0
  !     sbeame            electron source form beam
  !     sbeam(j,i)       thermal ion source from beam, grid point j, species nameb(i)
  !
  !     enbeam(j,i)      j=1,2..nj,i=1..nbion  fast ion density due to beam,
  !                      grid point j, species nameb(i).
  !     wbeam(j)         beam stored energy density
  !     enalp(J)         alpha particel den sity
  !     walp(j)          walp particle stored energy density
  !                      (walp is initial cond for w_alpha)
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
  !                         fix_edge_te_bc  (= any value in (0.0,1.0] )     
  !                         fix_edge_ti_bc  (= any value in (0.0,1.0] )
  !                         fix_edge_rot_bc (= any value in (0.0,1.0] )
  !                         fix_edge_ni_bc(nion)
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
  !     these values are cosntant during a call to gcnm but could
  !     change from one call to the next,as could rho,etcc so this
  !     is a fully dynamicallay controlled model.     
  !     fix_edge_te_bc      	    
  !     fix_edge_ti_bc           
  !     fix_edge_rot_bc 
  !     fix_edge_ni_bc(nion)


  !     currf(j)
  !     curbe(i)
  !     curbi(j)
  !     qrad(j)         radiated power
  !     totohm
  !     totbeam
  !     totboot
  !     totrf           total currents, amps
  !     dpedt(j), dpidt(j)
  !     qconde,qcondi
  !     qconve,qconvi
  !     qdelt
  !     qexch   
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
  !     cheinv(j)
  !     chiinv(j)        electron and ion thermal cond.
  !     xkineo(j)        ion neoclassical thermal conductivity
  !     chiepc(j)        electron paleoclassical thermal diffusivity
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
  ! ----------------------------------------------------- 11/22/94 --- HSJ
  !


  !
  USE nrtype,                    ONLY : DP,I4B,I2B

  USE error_handler,             ONLY : lerrno,terminate,iomaxerr

  USE iterdbmd_gcnmp,            ONLY : iterdsc,irwflag,iterdb_file_name,     &
                                        iterdb_outpt


#ifdef GCNMP
#define USESUBS
#endif
#ifdef NFREYA
#define USESUBS
  USE nfreya_version,            ONLY : p_nfreya_ver
#endif

#ifdef USESUBS
  USE io_gcnmp,                  ONLY : ncrt,nlog,niterdb
  USE gcnmp_version,             ONLY : gcnmp_ver
#else
  USE io,                        ONLY : ncrt,versid,nlog => nitre
  USE iterdbmd,                  ONLY : niterdb
#endif



  USE grid_class,                ONLY : nj,psir_grid,rho_grid,rho_gridn,     &
                                        eps,xhm2,xi11,xi33,xips,xhm20,xi110, &
                                        xi330,xips0,rho_mhd_gridnpsi



  USE vector_class,              ONLY : new_Vector,get_element,  &
                                        real_mult_Vector,list,   &
                                        zero_Vector,load_Vector, &
                                        delete_vector,get_values,&
                                        delete_vector_nf

  USE ions_gcnmp,                ONLY : namep,namei,nion,                     &
                                        namen,nprim,nimp,nneu,z,zsq,zeff,     &
                                        name_size,fd_thermal,nprimp1

  USE shot_info,                 ONLY : shot_id

  USE MPI_data,                  ONLY : mpiierr,myid

  USE Plasma_properties ,        ONLY : dischg,profile,mhd_dat,              &
                                        diffuse,pwrden, wpdot,prtcl_src,     &
                                        fus_prod,pellet,neut_beam,           &
                                        plasma_frequencies

  USE solcon_gcnmp,              ONLY : time,eqtime,tGCNMf,tGCNMs

  USE source_terms_gcnmp,        ONLY : stsource,scx,sion,srecom,sbcx,sbeame, &
                                        dudtsv,sbeam,stfuse,sbfuse,brems_tot, &
                                        qrad_tot

  USE neutral_data,              ONLY : enn,ennw,volsn,ennv

  USE fast_ion_data_gcnmp,       ONLY : walp,enalp,w_alpha

  USE  neutral_beams,            ONLY : nbion,fd_beam,nameb,enbeam,storqueb,  &
                                        enbeam_tot,wbeam,neut_beam_allocate, &
                                        bptor 

!  USE nub,                       ONLY : bptor  for old Nfreya

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

  USE nf_param,                  ONLY : kcm

  USE tglfin,                    ONLY : tglf_p_output,tglf_e_output,tglf_m_output
  
  USE fusion_gcnmp,              ONLY : pfuse_tot,ddn_fus,dtn_fus,tt2n_fus,                  &
                                        ddp_fus,dhe_fus,he3_frac,he3_thermal_spin_pol

  USE glf23_gcnmp,               ONLY : glf_etg_output,glf_gamma_net_i_output,               &
                                        glf_gamma_net_e_output,glf_anfreq_output,            &
                                        glf_anfreq2_output,glf_anrate_output,                &
                                        glf_anrate2_output,                                  &
                                        glf_p_output,glf_e_output,glf_m_output,              &
                                        glf_eff_chi_output

  USE pl_freq,                   ONLY : allocate_plasma_freq

  IMPLICIT  NONE
  





  REAL(DP),ALLOCATABLE,DIMENSION(:) ::   work,worknpsi



  REAL(DP) elmt
  INTEGER(I4B) ie, j,jb,jp,jj,ji,jn,k,ntot,deptot,nmax,ngW20,ok
  INTEGER(I4B), PARAMETER :: one = 1_I4B
  INTEGER(I4B), PARAMETER :: ke_bm  = 3_I4B ! 3 beam energy componenets
  INTEGER(I2B) i
  INTEGER(I2B) get_next_io_unit
  LOGICAL    monotonic
  CHARACTER  starflag*2, headerline*132,line*132,time_str*36
  CHARACTER  label*132,bc_asc_time*24,st_asc_time*24,tlabel*132
  CHARACTER(len = name_size) :: tname
  CHARACTER(LEN = *), PARAMETER :: dimensionless = 'dimensionless'
  INTERFACE      
     SUBROUTINE check_monotonic (array, n, monotonic, incr)   
       USE nrtype,   ONLY : I4B,DP
       INTEGER(I4B), INTENT(IN) :: n,incr
       REAL(DP),INTENT(IN),DIMENSION(:) :: array(n)
       LOGICAL, INTENT(OUT) ::    monotonic
     END SUBROUTINE check_monotonic
  END INTERFACE

  iterdsc = 1                   ! always write descriptor,reads will fail otherwise
  !
7 FORMAT (         a   )
8 FORMAT (5(2x,    a  ))        ! common character write/read format
9 FORMAT (5(2x,   i6  ))        ! common integer   write/read format
!10 FORMAT (5(2x,1pe14.4))        ! common floating  write/read format
10 FORMAT (5(2x,1pe14.6))        ! common floating  write/read format hsj 10/1/10
11 FORMAT (2x,i6,2x,1pe14.4)
  !
  ! OPENing method depends on whether read or write, and whether or not first time

  niterdb = get_next_io_unit ()


  IF (irwflag .EQ. 0) THEN      ! WRITE to file
     iterdb_outpt = iterdb_file_name(1:LEN_TRIM(iterdb_file_name))
     OPEN (unit= niterdb, file = iterdb_outpt, status = 'UNKNOWN',ERR=88)
     go to 89
88   lerrno = 5_I4B
     CALL  terminate(lerrno,nlog)
89   CONTINUE
  ELSE                          ! READ from existing file
     OPEN (unit = niterdb, file = iterdb_file_name, status = 'UNKNOWN', &
          err = 2)
     go to 3
2    WRITE  (ncrt, 4)  iterdb_file_name(1:LEN_TRIM(iterdb_file_name))
     WRITE  (nlog, 4)  iterdb_file_name(1:LEN_TRIM(iterdb_file_name))
4    FORMAT (                                                     / &
          ' ERROR: subroutine ITER_DBASE has encountered an error' / &
          '        the ITER database file "', a, '" cannot be opened')
     lerrno = 3_I4B
     CALL  terminate(lerrno,nlog)
3    CONTINUE
  END IF



  !
  ! write a header line each time routine is called:
  !   header lines are identified by     **                  in first two columns
  !   other comment lines (if any) have  *b  (where b=blank) in first two columns
  !


  IF (irwflag .EQ. 0) THEN
#ifdef GCNMP
     WRITE  (niterdb, 1) time,tGCNMf,time_bc, & ! for gcnmp  statefile write
          TRIM(gcnmp_ver)
1    FORMAT ('**   time = ', 1pe14.6, &
          ' tGCNMf = ',1pe14.6, ' time_bc = ',1pe14.6,5X,' GCNMP :',a)
#endif
#ifdef NFREYA
     WRITE  (niterdb, 21) time,tGCNMf,time_bc, & ! for P_Nfreya statefile write
          TRIM(p_nfreya_ver)
21    FORMAT ('**   time = ', 1pe14.6, &
          ' tGCNMf = ',1pe14.6, ' time_bc = ',1pe14.6,5X,'P_Nfreya :',a)
#endif
#ifdef ONETWO
     WRITE  (niterdb, 14) time,tGCNMf,time_bc, &  ! for onetwo statefile write
          TRIM(versid)
14    FORMAT ('**   time = ', 1pe14.6, &
          ' tGCNMf = ',1pe14.6, ' time_bc = ',1pe14.6,5X,'ONETWO :',a)
#endif

     !
     ! --- check that psir_grid is monotonic; it may not be in certain cases where
     ! --- the current profile was evolved
     !

     IF(.NOT. ALLOCATED(work))ALLOCATE (work(nj))  ! for writes nj is known 
     work(1:nj) = get_values(psir_grid)
     CALL check_monotonic (work, nj, monotonic, one)

     IF (.NOT. monotonic) THEN    ! psir_grid is not monotonic
        WRITE  (niterdb, 5) time
        WRITE  (ncrt,5)time
5       FORMAT (' the psi grid, calculated from the poloidal' / &
             ' B field evolution, is not monotonic'        / &
             ' therefore the data at this time (', 1pe12.6, &
             ') was not calculated')
        lerrno = 41_I4B
        CALL  terminate(lerrno,nlog)
     END IF
  ELSE
     READ (niterdb, 7) headerline
  END IF






  !
  ! --- some scalar quantities:
  ! --- integer and character parameters
  !

  rwblock_A  :  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 250)
250  FORMAT ('*  ishot  : shot number')
     WRITE  (niterdb, 9) shot_id%shot_nmbr
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 300)
300  FORMAT ('*  nj : the size of the vectors printed', &
          ' in this file')
     WRITE  (niterdb, 9) nj
     WRITE  (niterdb, 305)
305  FORMAT ('*  npsi : the size of mhd  vectors printed', &
          ' in this file')
     WRITE  (niterdb, 9) mhd_dat%npsi
     WRITE  (niterdb, 307)
307  FORMAT ('*nr,nz : the size of (R,Z) grid for mhd calculations')
     WRITE  (niterdb, 9) dischg%nr_mhd,dischg%nz_mhd
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 308)
308  FORMAT ('* nlimiter : the size of (R,Z) limiter contour')
     WRITE  (niterdb, 9)dischg%nlimiter
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 310)
310  FORMAT ('*  nion : the number of ion species')
     WRITE  (niterdb, 9) nion
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 320)
320  FORMAT ('*  nprim : the number of primary ion species,frac d in dt mix')
     WRITE  (niterdb, 11) nprim, fd_thermal
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 330)
330  FORMAT ('*  nimp : the number of impurity ion species')
     WRITE  (niterdb, 9) nimp
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 340)
340  FORMAT ('*  nneu : the number of neutral ion species')
     WRITE  (niterdb, 9) nneu
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 350)
350  FORMAT ('*  nbion :  number of beam ion species and d fraction ')
     WRITE  (niterdb, 11) nbion , fd_beam
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 360)
360  FORMAT ('*  namep : name(s) of primary ion species')
     WRITE  (niterdb, 8) (namep(i),i=1,nprim)
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 370)
370  FORMAT ('*  namei : name(s) of impurity ion species')
     WRITE  (niterdb, 8) (namei(i),i=1,nimp)
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 380)
380  FORMAT ('*  namen : name(s) of neutral ion species')
     WRITE  (niterdb, 8) (namen(i), i=1,nneu)

     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 385)
385  FORMAT ('*  nameb : name(s) of beam ion species')
     WRITE  (niterdb, 8) (nameb(i), i=1,nbion)

     IF (iterdsc .NE. 0) &
          WRITE  (niterdb,386)
386  FORMAT('* name pellet species: ' )
     WRITE  (niterdb,8)pellet%name

     WRITE  (niterdb, 388)
388  FORMAT ('* total   number of dep variables ')
     ! tot dep vars in Gcnmp
     ntot = nion+dp4      
     WRITE  (niterdb, 11) ntot
  ELSE rwblock_A
     READ   (niterdb, 7) starflag
     READ   (niterdb, 9) shot_id%shot_nmbr
     READ   (niterdb, 7) starflag
     READ   (niterdb, 9) nj 
     READ   (niterdb, 7) starflag
     READ   (niterdb, 9) mhd_dat%npsi
     READ   (niterdb, 7) starflag
     READ   (niterdb, 9) dischg%nr_mhd,dischg%nz_mhd
     READ   (niterdb, 7) starflag
     READ   (niterdb, 9) dischg%nlimiter
     READ   (niterdb, 7) starflag
     READ   (niterdb, 9) nion
     ntot = nion+dp4 
     READ   (niterdb, 7) starflag
     READ   (niterdb, 11) nprim,fd_thermal
     READ   (niterdb, 7) starflag
     READ   (niterdb, 9) nimp
     READ   (niterdb, 7) starflag
     READ   (niterdb, 9) nneu
     READ   (niterdb, 7) starflag
     READ   (niterdb, 11) nbion,fd_beam
     READ   (niterdb, 7) starflag
     nprimp1 = nprim+1

     IF( .NOT. ASSOCIATED(namep))THEN
        ALLOCATE ( namep(nprim),namei(nimp),namen(nneu),nameb(nbion))
        ALLOCATE (profile%en(nion))        ! 2d array 
        ALLOCATE (profile%flux(ntot))      ! 2d array 
        ALLOCATE (profile%flux_conv(ntot)) ! 2d array 
        DO jj=1,ntot
           profile%flux_conv(jj) = zero_vector(nj)
        ENDDO
        ALLOCATE (profile%zsq(ntot))       ! 2d array 
        ALLOCATE (profile%z(ntot))         ! 2d array 
     ENDIF

     READ   (niterdb, 8) (namep(i),i=1,nprim)
     READ   (niterdb, 7) starflag
     READ   (niterdb, 8) (namei(i),i=1,nimp)
     READ   (niterdb, 7) starflag
     READ   (niterdb, 8) (namen(i),i=1,nneu)
     READ   (niterdb, 7) starflag
     READ   (niterdb, 8) (nameb(i),i=1,nbion)
     READ   (niterdb, 7) starflag
     READ   (niterdb, 8) pellet%name
     READ   (niterdb, 7) starflag
     READ   (niterdb, 9) deptot

     IF(deptot .NE. ntot)THEN
        WRITE(ncrt,495)deptot,ntot
        WRITE(nlog,495)ntot,deptot
495     FORMAT (' ERROR: subroutine ITER_DBASE has encountered an error',/, &
             ' number of dependent variables in GCNMp = ',i5, /, &
             ' Number given in iterdb files is',i5)
        lerrno = 7_I4B
        CALL  terminate(lerrno,nlog)
     ENDIF
  END IF rwblock_a




  nmax = MAX(nj, mhd_dat%npsi,dischg%nr_mhd,dischg%nz_mhd,dischg%nlimiter)
  IF(ALLOCATED(work))DEALLOCATE (work)  
  ALLOCATE(work(nmax))                                     




  !
  ! --- real (i.e., floating point) parameters
  !
  rwblock_B: IF (irwflag .EQ. 0) THEN
     !  
     profile%te0     = profile%te%data(1)

     profile%ti0     = profile%ti%data(1)

     !
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 490)
490  FORMAT ('*  time : initial conditions and time of GCNM start')
     WRITE  (niterdb, 10) time

     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 499)
499  FORMAT ('*  Rgeom : major radius of geometric', &
          ' center at elevation of magnetic axis, meters')
     WRITE  (niterdb, 10) dischg%rgeom

     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 500)
500  FORMAT ('*  Btgeom : BT at rgeom, tesla')
     WRITE  (niterdb, 10) dischg%btgeom


!     IF (iterdsc .NE. 0) &
!          WRITE  (niterdb, 501)
!501  FORMAT ('*  Rmag : major radius of magnetic axis, meters')
!     WRITE  (niterdb, 10) dischg%rmag


     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 502)
502  FORMAT ('*  R0 : major radius of vacuum btor ref', &
          ' location, meters')
     WRITE  (niterdb, 10) dischg%rmajor
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 503)
503  FORMAT ('*  kappa : plasma elongation')
     WRITE  (niterdb, 10) dischg%kappa
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 504)
504  FORMAT ('*  delta : plasma triangularity')
     WRITE  (niterdb, 10) dischg%deltao
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 505)
505  FORMAT ('*  pindent : plasma indentation')
     WRITE  (niterdb, 10) dischg%pindento
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 506)
506  FORMAT ('*  volo : plasma volume, meters**3')
     WRITE  (niterdb, 10) dischg%volo

     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 507)
507  FORMAT ('*  cxareao : plasma cross-sectional area, meters**2')
     WRITE  (niterdb, 10) dischg%areao
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 508)
508  FORMAT ('*  Btor : vacuum toroidal field at R0, tesla')
     WRITE  (niterdb, 10) mhd_dat%btor
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 509)
509  FORMAT ('*  total, ohmic, bootstrap, beam and RF', &
          ' currents, amps')
     WRITE  (niterdb, 10) mhd_dat%tot_cur, mhd_dat%totohm_cur, &
          mhd_dat%totboot_cur, mhd_dat%totbeam_cur, mhd_dat%totrf_cur
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 510)
510  FORMAT ('*  betap : poloidal beta, ')
     WRITE  (niterdb, 10) mhd_dat%betap
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 511)
511  FORMAT ('*  beta : toroidal beta, ')
     WRITE  (niterdb, 10)  mhd_dat%beta
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 512)
512  FORMAT ('*  ali : plasma inductance, ')
     WRITE  (niterdb, 10) mhd_dat%ali
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 513)
513  FORMAT ('*  te0 : central electron temperature, keV')
     WRITE  (niterdb, 10) profile%te0
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 514)
514  FORMAT ('*  ti0 : central ion temperature, keV')
     WRITE  (niterdb, 10) profile%ti0
     !8/16/07
     WRITE  (niterdb, 515)
515  FORMAT ('* rma,zma : magnetic axis r,z coords, meters')
     WRITE  (niterdb, 10) dischg%rma,dischg%zma
     WRITE  (niterdb, 516)
516  FORMAT ('* rsep,zsep : x point r,z,coords, meters')
     WRITE  (niterdb, 10) dischg%rsep,dischg%zsep
     WRITE  (niterdb, 517)
517  FORMAT ('* rplasmax,zplasmax: max R,Z of plasma, meters')
     WRITE  (niterdb, 10) dischg%rplasmax,dischg%zplasmax
     WRITE  (niterdb, 518)
518  FORMAT ('* rplasmin,zplasmin: min R,Z,of plasma, meters')
     WRITE  (niterdb, 10) dischg%rplasmin,dischg%zplasmin
     WRITE  (niterdb, 519)
519  FORMAT ('* psiaxis ,psibdry : psi axis and boundary values, volt sec/rad')
     WRITE  (niterdb, 10) mhd_dat%psiaxis,mhd_dat%psibdry


  ELSE rwblock_B
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) time ;  shot_id%shot_time = time ; tGCNMs = time
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) dischg%rgeom
     READ   (niterdb,  7)starflag
     READ   (niterdb, 10) dischg%btgeom
!     READ   (niterdb,  7) starflag
!     READ   (niterdb, 10) dischg%rmag
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) dischg%rmajor ; mhd_dat%R0 = dischg%rmajor
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) dischg%kappa
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) dischg%deltao
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) dischg%pindento
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) dischg%volo

     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) dischg%areao
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) mhd_dat%btor
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) mhd_dat%tot_cur, mhd_dat%totohm_cur, &
          mhd_dat%totboot_cur, mhd_dat%totbeam_cur, mhd_dat%totrf_cur
     ibcur =1_I4B ;irfc = 1_I4B
     IF(ABS(mhd_dat%totbeam_cur) .LT. 1.e-5)ibcur = 0
     IF(ABS(mhd_dat%totrf_cur)   .LT. 1.e-5)irfc  = 0
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) mhd_dat%betap
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) mhd_dat%beta
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) mhd_dat%ali
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) profile%te0
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) profile%ti0
     !8/16/07
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) dischg%rma,dischg%zma
     dischg%rmag = dischg%rma
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) dischg%rsep,dischg%zsep
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) dischg%rplasmax,dischg%zplasmax
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) dischg%rplasmin,dischg%zplasmin
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) mhd_dat%psiaxis, mhd_dat%psibdry
  END IF rwblock_B




  eqtime = time            !equilibirum time. It is assumed that the quatities
  !fcap,gcap,hcap,eps,xhm2,xi11,xi33,xips
  !were calcualted at this time. Hence subrotuine 
  !neointrp does nothing




  !
  ! --- psir grid
  !
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 1010)
1010 FORMAT ('*  psi on rho grid, volt*second/radian')
     work(1:nj) = get_values(psir_grid)
     WRITE  (niterdb, 10) (work(j), j=1,nj)
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (work(j), j=1,nj)
     psir_grid = new_Vector(nj,work)
     CALL check_monotonic (work, nj, monotonic, one)
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
  END IF


  !
  ! --- rho grid
  !
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 1020)
1020 FORMAT ('*  rho grid, meters')
     work(1:nj) = get_values(rho_grid)
     WRITE (niterdb, 10) (work(j), j=1,nj)
  ELSE
     READ  (niterdb,  7) starflag
     READ  (niterdb, 10) (work(j)   , j=1,nj)
     rho_grid = new_Vector(nj,work)
     elmt = get_element(rho_grid,nj)
     elmt = 1.0_DP/elmt
     rho_gridn = real_mult_Vector (elmt,rho_grid)
  END IF


  !
  ! --- r,z,mhdgrid
  !
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 1021)
1021 FORMAT ('* R(1:nr_mhd):  mhd grid, meters')
     WRITE  (niterdb, 10) (dischg%rmhdgrid%data(j), j=1,dischg%nr_mhd)
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (work(j),j=1,dischg%nr_mhd)
     dischg%rmhdgrid = new_Vector(dischg%nr_mhd,work)
  END IF
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) & 
          WRITE  (niterdb, 1022)
1022 FORMAT ('* Z(1:nz_mhd):  mhd grid, meters')
     WRITE  (niterdb, 10) (dischg%zmhdgrid%data(j), j=1,dischg%nz_mhd)
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (work(j),j=1,dischg%nz_mhd)
     dischg%zmhdgrid = new_Vector(dischg%nz_mhd,work)
  END IF


  !
  ! --- limiter points
  !
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 1023)
1023 FORMAT ('* Rlimiter(1:nlimiter) : r coords of limiter, meters')
     WRITE  (niterdb, 10) (dischg%rlimiter%data(j), j=1,dischg%nlimiter)
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (work(j), j=1,dischg%nlimiter)
     dischg%rlimiter = new_Vector(dischg%nlimiter,work)
  END IF
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) & 
          WRITE  (niterdb, 1024)
1024 FORMAT ('* Zlimiter(1:nlimiter) : z coords of limiter, meters')
     WRITE  (niterdb, 10) (dischg%zlimiter%data(j), j=1,dischg%nlimiter)
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (work(j), j=1,dischg%nlimiter)
     dischg%zlimiter = new_Vector(dischg%nlimiter,work)
  END IF

  !
  ! --- rho on mhd grid:
  !
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 1025)
1025 FORMAT ('* rho(1:npsi) : rho on psival mhd grid, meters')
     WRITE  (niterdb, 10) (rho_mhd_gridnpsi%data(j), j=1,mhd_dat%npsi)
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (work(j), j=1,mhd_dat%npsi)
     rho_mhd_gridnpsi = new_Vector(mhd_dat%npsi,work)
  END IF



  !
  ! psival, mhd psi grid
  !

  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 1026)
1026 FORMAT ('* psival(1:npsi): psi on mhd grid, volt sec/rad')
     WRITE (niterdb, 10) (mhd_dat%psivalnpsi%data(j), j=1,mhd_dat%npsi)
  ELSE
     READ  (niterdb,  7) starflag
     READ  (niterdb, 10) (work(j), j=1,mhd_dat%npsi) 
     mhd_dat%psivalnpsi = new_Vector(mhd_dat%npsi,work)
  END IF


  !
  ! ravg,on  mhd psi grid
  !

  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 1027)
1027 FORMAT ('* ravg(1:npsi): <R> avg radius on mhd grid, meters')
     WRITE (niterdb, 10) (mhd_dat%ravgnpsi%data(j), j=1,mhd_dat%npsi)
  ELSE
     READ  (niterdb,  7) starflag
     READ  (niterdb, 10) (work(j), j=1,mhd_dat%npsi)
     mhd_dat%ravgnpsi = new_Vector(mhd_dat%npsi,work)
  END IF
  !
  ! ravgi,on  mhd psi grid
  !

  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 1028)
1028 FORMAT ('* ravgi(1:npsi):<1/R>  on mhd grid, 1/meters')
     WRITE (niterdb, 10) (mhd_dat%ravginpsi%data(j), j=1,mhd_dat%npsi)
  ELSE
     READ  (niterdb,  7) starflag
     READ  (niterdb, 10) (work(j), j=1,mhd_dat%npsi)
     mhd_dat%ravginpsi  = new_Vector(mhd_dat%npsi,work)
  END IF


  !
  ! psi values over mhdgrid:
  !
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 1030)
1030 FORMAT ('* psi(1:nr,1:nz) over (R,Z) grid, volt sec/rad')
     WRITE (niterdb, 10) ((mhd_dat%psi(i,j),i=1,dischg%nr_mhd), j=1,dischg%nz_mhd)
  ELSE
     READ  (niterdb,  7) starflag
     IF(ASSOCIATED(mhd_dat%psi))DEALLOCATE(mhd_dat%psi)
     ALLOCATE(mhd_dat%psi(dischg%nr_mhd,dischg%nz_mhd))
     READ  (niterdb, 10) ((mhd_dat%psi(i,j),i=1,dischg%nr_mhd), j=1,dischg%nz_mhd) 
  END IF

  !
  ! fpsi
  !
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 1031)
1031 FORMAT ('* fpsi(1:npsi): f of psi (= R*Bt) on mhd grid, tesla meters')
     WRITE (niterdb, 10) (mhd_dat%fpsinpsi%data(j),j=1,mhd_dat%npsi)
  ELSE
     READ  (niterdb,  7) starflag
     READ  (niterdb, 10) (work(j),j=1,mhd_dat%npsi)
     mhd_dat%fpsinpsi   = new_Vector(mhd_dat%npsi,work)
  END IF

  !
  ! pprim
  !
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 1033)
1033 FORMAT ('* pprim(1:nj): dp/dpsi on transport grid, newton/(m**2-volt-sec) = amp/m**3')
     WRITE (niterdb, 10) (mhd_dat%pprim%data(j),j=1,nj)
  ELSE
     READ  (niterdb,  7) starflag
     READ  (niterdb, 10) (work(j),j=1,nj)
     mhd_dat%pprim = new_Vector(nj,work)
  END IF

  !
  ! ffprim
  !
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 1034)
1034 FORMAT ('* f*df/dpsi(1:nj): on transport grid, kg/(A sec^2) ') 
     WRITE (niterdb, 10) (mhd_dat%ffprim%data(j),j=1,nj)
  ELSE
     READ  (niterdb,  7) starflag
     READ  (niterdb, 10) (work(j),j=1,nj)
     mhd_dat%ffprim  = new_Vector(nj,work)
  END IF
  !
  ! flux surface average bpoloidal
  ! 

  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 1035)
1035 FORMAT ('* <Bp>(1:nj): flux avg B poloidal on transport grid, Tesla') 
     WRITE (niterdb, 10) (mhd_dat%bp%data(j),j=1,nj)
  ELSE
     READ  (niterdb,  7) starflag
     READ  (niterdb, 10) (work(j),j=1,nj)
     mhd_dat%bp = new_Vector(nj,work)
  END IF

  !
  ! bpoloidal on rmajor (NOT flux surface average)
  ! 

  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 1036)
1036 FORMAT ('*  Bprmaj(1:nj): B poloidal on transport grid (not a flux surface average), Tesla') 
     WRITE (niterdb, 10) (mhd_dat%bprmaj%data(j),j=1,nj)
  ELSE
     READ  (niterdb,  7) starflag
     READ  (niterdb, 10) (work(j),j=1,nj)
     mhd_dat%bprmaj = new_Vector(nj,work)
  END IF

  !
  ! B  total on rmajor (NOT flux surface average)
  ! 

  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 1037)
1037 FORMAT ('*  Btotrmaj(1:nj): B total on transport grid (not a flux surface average), Tesla') 
     WRITE (niterdb, 10) (mhd_dat%btotrmaj%data(j),j=1,nj)
  ELSE
     READ  (niterdb,  7) starflag
     READ  (niterdb, 10) (work(j),j=1,nj)
     mhd_dat%btotrmaj = new_Vector(nj,work)
  END IF



  !
  ! mhd  pressure
  !
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 1038)
1038 FORMAT ('* press(1:nj): pressure on transport grid, newton/m^2 ') 
     WRITE (niterdb, 10) (profile%press%data(j),j=1,nj)
  ELSE
     READ  (niterdb,  7) starflag
     READ  (niterdb, 10) (work(j),j=1,nj)
     profile%press = new_Vector(nj,work)
  END IF

  !
  ! beam  pressure
  !
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 1040)
1040 FORMAT ('* pressb(1:nj):  beam pressure on transport grid, newton/m^2 ') 
     WRITE (niterdb, 10) (profile%pressb%data(j),j=1,nj)
  ELSE
     READ  (niterdb,  7) starflag
     READ  (niterdb, 10) (work(j),j=1,nj)
     profile%pressb = new_Vector(nj,work)
  END IF



  !
  ! --- fcap
  !
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 1051)
1051 FORMAT ('*  fcap(1:nj), (i.e., f(psilim)/f(psi)), ')
     work(1:nj) = get_values(mhd_dat%fcap)
     WRITE  (niterdb, 10) (work(j), j=1,nj)
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (work(j), j=1,nj)
     mhd_dat%fcap = new_Vector(nj,work)
  END IF


  !
  ! --- gcap
  !
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 1052)
1052 FORMAT ('*  gcap(1:nj), (i.e., <(grad rho)**2*(R0/R)**2>), ')
     work(1:nj) = get_values(mhd_dat%gcap)
     WRITE  (niterdb, 10) (work(j), j=1,nj)
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (work(j), j=1,nj)
     mhd_dat%gcap = new_Vector(nj,work)
  END IF


  !
  ! --- hcap
  !
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 1053)
1053 FORMAT ('*  hcap(1:nj), (i.e., (dvolume/drho)/(4*pi*pi*R0*rho)), ')
     work(1:nj) = get_values(mhd_dat%hcap)
     WRITE  (niterdb, 10) (work(j), j=1,nj)
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (work(j), j=1,nj)
     mhd_dat%hcap = new_Vector(nj,work)
  END IF


  !
  ! --- te (in keV)
  !
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 1045)
1045 FORMAT ('* Te(1:nj): electron temperature, keV')
     work(1:nj) = get_values(profile%te)
     WRITE  (niterdb, 10) (work(j), j=1,nj)
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (work(j), j=1,nj)
     profile%te = new_Vector(nj,work)
  END IF

  !
  ! --- ti (in keV)
  !
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 1047)
1047 FORMAT ('* Ti(1:nj)  ion temperature, keV')
     work(1:nj) = get_values(profile%ti)
     WRITE  (niterdb, 10) (work(j), j=1,nj)
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (work(j), j=1,nj)
     profile%ti = new_Vector(nj,work)
  END IF


  !
  ! --- safety factor
  !
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 1050)
1050 FORMAT ('*  q(1:nj)  (i.e., safety factor) profile, ')
     work(1:nj) = get_values(mhd_dat%q_value)
     WRITE  (niterdb, 10) (ABS (work(j)), j=1,nj)
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (work(j) , j=1,nj)
     mhd_dat%q_value = new_Vector(nj,work)
  END IF


  !
  ! --- electron density
  !
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 3054)
3054 FORMAT ('* ene(1:nj):  electron density, #/meter**3')
     work(1:nj) = get_values(profile%ene)
     WRITE  (niterdb, 10) (work(j), j=1,nj)
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (work(j)       , j=1,nj)
     profile%ene = new_Vector(nj,work)
  END IF



  !
  ! --- primary, impurity densities
  !
  jp = 0
  ji = 0
  DO jj=1,nion
     IF (jj .LE. nprim)  jp = jp + 1
     IF (jj .GT. nprim)  ji = ji + 1
     IF (irwflag .EQ. 0) THEN
        IF (iterdsc .NE. 0 .AND. jj .LE. nprim) &
             WRITE  (niterdb, 3055) namep(jp)
3055    FORMAT ('* en(1:nj)  primary ion density,', &
             ' #/meter**3, species: ', a)
        IF (iterdsc .NE. 0 .AND. jj .GT. nprim) &
             WRITE  (niterdb, 3056) namei(ji)
3056    FORMAT ('* eni(1:nj)  impurity ion density,', &
             ' #/meter**3, species: ', a)
        work(1:nj) = get_values(profile%en(jj))
        WRITE  (niterdb, 10) (work(j), j=1,nj)
     ELSE
        READ   (niterdb,  7) line
        !               CALL   extract_name(line,tname)
        !               IF(jj .LE. nprim)THEN 
        !                  namep(jj) = tname
        !               ELSE
        !                  namei(ji) = tname
        !               ENDIF
        READ   (niterdb, 10) (work(j)       , j=1,nj)
        profile%en(jj) = new_Vector(nj,work)
     END IF
  END DO


  !
  ! fusion rate density
  ! 

  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 3500)
3500 FORMAT ('* stfuse(1:nj): thermal fusion rate #/(m^3 sec)') 
     WRITE (niterdb, 10) (prtcl_src%stfuse%data(j),j=1,nj)
  ELSE
     READ  (niterdb,  7) starflag
     READ  (niterdb, 10) (work(j),j=1,nj)
     prtcl_src%stfuse = new_Vector(nj,work)
  END IF

  !
  ! beam fusion rate density
  !
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 4010)
4010 FORMAT ('* sbfuse(1:nj): beam-thermal  fusion rate #/(m^3 sec)') 
     WRITE (niterdb, 10) (prtcl_src%sbfuse%data(j),j=1,nj)
  ELSE
     READ  (niterdb,  7) starflag
     READ  (niterdb, 10) (work(j),j=1,nj)
     prtcl_src%sbfuse = new_Vector(nj,work)
  END IF

!
! pellet source density
! 
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 3037)pellet%name
3037 FORMAT ('* spellet(1:nj): pellet THERMAL  ion source  #/(m^3 sec),species: ',a) 
     WRITE (niterdb, 10) (prtcl_src%spellet%data(j),j=1,nj)
  ELSE
     READ  (niterdb,  7) starflag
     READ  (niterdb, 10) (work(j),j=1,nj)
     prtcl_src%spellet = new_Vector(nj,work)
  END IF



  !
  ! --- particle and energy  flux
  !

  jp = 0
  ji = 0
  DO jj=1,ntot
     IF (jj .LE. nprim)  jp = jp + 1
     IF (jj .GT. nprim)  ji = ji + 1
     IF (irwflag .EQ. 0) THEN
        IF (iterdsc .NE. 0 .AND. jj .LE. nprim) &
             WRITE  (niterdb, 4057) namep(jp)
4057    FORMAT ('* Flx(1:nj):  primary ion flux,', &
             ' #/(meter2 sec), species: ', a)
        IF (iterdsc .NE. 0 .AND. jj .GT. nprim .AND. jj .LE. nion) &
             WRITE  (niterdb, 4058) namei(ji)
4058    FORMAT ('* Flxi(1:nj):  impurity ion flux,', &
             ' #/(meter**2sec), species: ', a)
        IF (iterdsc .NE. 0 .AND. jj .EQ.  nion+1 ) &
             WRITE  (niterdb, 4059)
4059    FORMAT('* Flxee(1:nj):electron energy flux, Joules/(m**2 sec)')
        IF (iterdsc .NE. 0 .AND. jj .EQ.  nion+2 ) &
             WRITE  (niterdb, 4060)
4060    FORMAT('* Flxie(1:nj): ion energy flux, joules/(m**2 sec)')
        IF (iterdsc .NE. 0 .AND. jj .EQ.  nion+3 ) &
             WRITE  (niterdb, 4061)
4061    FORMAT('* Flxef(1:nj):Effective flux in Faradys law Tesla/sec')
        IF (iterdsc .NE. 0 .AND. jj .EQ.  nion+4 ) &
             WRITE  (niterdb, 4062)
4062    FORMAT('* Flxa(1:nj): toroidal momentum flux Kg/sec^2 ')
        work(1:nj) = get_values(profile%flux(jj))
        WRITE  (niterdb, 10) (work(j), j=1,nj)
     ELSE
        READ   (niterdb,  7) line
        READ   (niterdb, 10) (work(j)       , j=1,nj)
        profile%flux(jj) = new_Vector(nj,work)
     END IF
  END DO



  !
  ! --- tglf turbulent particle and energy  flux
  ! --- It is assumed that tglf always has 3 species;
  ! --- 1) electron
  ! --- 2) effective ion
  ! ----3) effective impurity, (see sub set_tglf_vars)
     IF (irwflag .EQ. 0) THEN
        IF (iterdsc .NE.0 )THEN
             WRITE  (niterdb, 9057)
9057         FORMAT ('*  tglf electron flux,', &
             ' #/(meter2 sec)' )
             WRITE  (niterdb, 10) (tglf_p_output(j,1), j=1,nj)
        ENDIF
        IF(iterdsc .NE. 0 )THEN
             WRITE  (niterdb, 9058)
9058         FORMAT ('*  tglf (effective) primary ion  flux,', &
             ' #/(meter2 sec)' )
             WRITE  (niterdb, 10) (tglf_p_output(j,2), j=1,nj)
        ENDIF
        IF(iterdsc .NE. 0 )THEN
             WRITE  (niterdb, 9059)
9059         FORMAT ('*  tglf (effective) impurity  ion  flux,', &
             ' #/(meter2 sec)' )
             WRITE  (niterdb, 10) (tglf_p_output(j,3), j=1,nj)
        ENDIF

        IF (iterdsc .NE. 0)THEN
             WRITE  (niterdb, 9061)
9061         FORMAT ('* tglf  electron energy flux', &
                     ' Joules/(meter**2sec)')
             WRITE  (niterdb, 10) (tglf_e_output(j,1), j=1,nj)
        ENDIF
        IF (iterdsc .NE. 0)THEN
             WRITE  (niterdb, 9062)
9062         FORMAT ('* tglf  effective primary ion  energy flux', &
                     ' Joules/(meter**2sec)')
             WRITE  (niterdb, 10) (tglf_e_output(j,2), j=1,nj)
        ENDIF
        IF (iterdsc .NE. 0)THEN
             WRITE  (niterdb, 9063)
9063         FORMAT ('* tglf effective impurity ion  energy flux', &
                     ' Joules/(meter**2sec)')
             WRITE  (niterdb, 10) (tglf_e_output(j,3), j=1,nj)
        ENDIF
     ELSE
        IF(ALLOCATED(tglf_p_output))DEALLOCATE(tglf_p_output)
        ALLOCATE(tglf_p_output(nj,3))
        READ   (niterdb,  7) line
        READ   (niterdb, 10) (tglf_p_output(j,1)       , j=1,nj)
        READ   (niterdb,  7) line
        READ   (niterdb, 10) (tglf_p_output(j,2)       , j=1,nj)
        READ   (niterdb,  7) line
        READ   (niterdb, 10) (tglf_p_output(j,3)       , j=1,nj)
        IF(ALLOCATED(tglf_e_output))DEALLOCATE(tglf_e_output)
        ALLOCATE(tglf_e_output(nj,3))
        READ   (niterdb,  7) line
        READ   (niterdb, 10) (tglf_e_output(j,1)       , j=1,nj)
        READ   (niterdb,  7) line
        READ   (niterdb, 10) (tglf_e_output(j,2)       , j=1,nj)
        READ   (niterdb,  7) line
        READ   (niterdb, 10) (tglf_e_output(j,3)       , j=1,nj)
     END IF

  !
  ! --- ion species sources
  ! 

  IF( .NOT. ALLOCATED(stsource))THEN 
     ALLOCATE(stsource(nion,nj))
     stsource(:,:)  = zeroc
  ENDIF
  IF( .NOT. ALLOCATED(sion))THEN
     ALLOCATE(sion(nj,nion))
     sion(:,:)      = zeroc
  ENDIF

!  IF( .NOT. ALLOCATED(srecom))THEN
!      ALLOCATE(srecom(nj,nion))
!      srecom(:,:) = zeroc
!  ENDIF
  IF(.NOT. ASSOCIATED(prtcl_src%srecom))ALLOCATE(prtcl_src%srecom(nion))

  IF( .NOT. ALLOCATED(sbcx))THEN
      ALLOCATE(sbcx(nj,nion))
      sbcx(:,:) = zeroc
  ENDIF
  IF( .NOT. ALLOCATED(scx))THEN
      ALLOCATE(scx(nj,nion))
      scx(:,:)   =  zeroc
  ENDIF
  IF( .NOT. ALLOCATED(dudtsv))THEN
      ALLOCATE(dudtsv(ntot,nj))
      dudtsv(:,:) = zeroc
  ENDIF



  ji = 0
  DO jj=1,nion
     IF(nprim .LT. jj .AND. jj .LE. nion)ji=ji+1
     IF (irwflag .EQ. 0) THEN
        label = '*  sion : source due to ionization,'      &
             //' #/(meter**3*second), species: '
        IF(jj .LE. nprim)label = label(1:LEN_TRIM(label))//namep(jj)
        IF(jj .GT. nprim)label = label(1:LEN_TRIM(label))//namei(ji)
        IF (iterdsc .NE. 0) &
             WRITE  (niterdb, FMT='(a)')label(1:LEN_TRIM(label))
        WRITE  (niterdb, 10) (sion(j,jj), j=1,nj)
        !
        label = '*  srecom : source due to recombination,'  &
             //' #/(meter**3*second), species: '
        IF(jj .LE. nprim)label = label(1:LEN_TRIM(label))//namep(jj)
        IF(jj .GT. nprim)label = label(1:LEN_TRIM(label))//namei(ji)
        IF (iterdsc .NE. 0) &
             WRITE  (niterdb, FMT='(a)')label(1:LEN_TRIM(label))
        work(1:nj) = get_values(prtcl_src%srecom(jj))
!        WRITE  (niterdb, 10) (srecom(j,jj), j=1,nj)
        WRITE  (niterdb, 10) (work(j), j=1,nj)

        label = '*  scx : source due to cx thermal neut.,'  &
             // ' #/(meter**3*second), species: '
        IF(jj .LE. nprim)label = label(1:LEN_TRIM(label))//namep(jj)
        IF(jj .GT. nprim)label = label(1:LEN_TRIM(label))//namei(ji)
        IF (iterdsc .NE. 0) &
             WRITE  (niterdb, FMT='(a)')label(1:LEN_TRIM(label))
        WRITE  (niterdb, 10) (scx(j,jj), j=1,nj)


        label = '*  sbcx : sink due to cx with beam neut.,'  &
             //' #/(meter**3*second), species: '
        IF(jj .LE. nprim)label = label(1:LEN_TRIM(label))//namep(jj)
        IF(jj .GT. nprim)label = label(1:LEN_TRIM(label))//namei(ji)
        IF (iterdsc .NE. 0) &
             WRITE  (niterdb, FMT='(a)')label(1:LEN_TRIM(label))
        WRITE  (niterdb, 10) (sbcx(j,jj), j=1,nj)

        label =  '*  s : total source rate,' &
             // ' #/(meter**3*second), species: '
        IF(jj .LE. nprim)label = label(1:LEN_TRIM(label))//namep(jj)
        IF(jj .GT. nprim)label = label(1:LEN_TRIM(label))//namei(ji)
        IF (iterdsc .NE. 0) &
             WRITE  (niterdb, FMT='(a)')label(1:LEN_TRIM(label))
        WRITE  (niterdb, 10) (stsource(jj,j), j=1,nj)
        !
        ! 
     ELSE
        READ   (niterdb,  7) starflag
        READ   (niterdb, 10) (sion(j,jj)  , j=1,nj)
        READ   (niterdb,  7) starflag
!        READ   (niterdb, 10) (srecom(j,jj), j=1,nj)
        READ   (niterdb, 10)  (work(j), j=1,nj)
        prtcl_src%srecom(jj) = new_Vector(nj,work)
        READ   (niterdb,  7) starflag
        READ   (niterdb, 10) (scx(j,jj)   , j=1,nj)
        READ   (niterdb,  7) starflag
        READ   (niterdb, 10) (sbcx(j,jj)  , j=1,nj)
        READ   (niterdb,  7) starflag
        READ   (niterdb, 10) (stsource(jj,j)     , j=1,nj)
     END IF
  END DO


  ji= 0
  DO jj=1,ntot ! ntot = nion +dp4 . ntot is not read in, must agree with caller
     IF(nprim .LT. jj .AND. jj .LE. nion)ji=ji+1
     IF (irwflag .EQ. 0) THEN
        label = '*  dudt : s dot #/(meter**3*second), species: '
        IF (iterdsc .NE. 0 .AND. jj .LE. nprim)THEN
           label = label(1:LEN_TRIM(label))//namep(jj)
        ELSEIF (iterdsc .NE. 0 .AND. jj .LE. nion)THEN
           label = label(1:LEN_TRIM(label))//namei(ji)
        ELSE 
           IF( jj == nion +1)THEN
              label = '*  dtedt : s dot,Te Kev/sec:'
           ELSEIF( jj == nion+2)THEN
              label = '*  dtidt : s dot,Ti Kev/sec:'
           ELSEIF( jj == nion+3)THEN
              label = '*  drbpdt : s dot,rbp tesla m/sec:'
           ELSEIF( jj == nion+dp4)THEN
              label = '*  dwdt : s dot,angrot,1/sec^2:'
           ENDIF
        ENDIF
        WRITE (niterdb,FMT='(a)')label(1:LEN_TRIM(label))
        WRITE  (niterdb, 10) (dudtsv(jj,j), j=1,nj)
     ELSE
        READ   (niterdb,  7) starflag
        READ   (niterdb, 10) (dudtsv(jj,j), j=1,nj)
     ENDIF
  ENDDO


  !
  ! --- fast ion density
  !
  IF( .NOT. ALLOCATED(enbeam)) ALLOCATE(enbeam(nj,nbion))
  IF( .NOT. ALLOCATED(enbeam_tot)) ALLOCATE(enbeam_tot(nj))
  enbeam_tot(:) = zeroc
  DO jn =1,nbion
     IF (irwflag .EQ. 0) THEN
        IF (iterdsc .NE. 0) &
             WRITE  (niterdb, 4000) nameb(jn)
4000    FORMAT ('*  fast ion density, #/meter**3, species: ', a)
        WRITE  (niterdb, 10) (enbeam(j,jn), j=1,nj)
     ELSE
        READ   (niterdb,FMT='(a)')line
        !              CALL   extract_name(line,nameb(jn))
        READ   (niterdb, 10) (enbeam(j,jn),j=1,nj)
        enbeam_tot(:) = enbeam_tot(:) + enbeam(:,jn)
     END IF
  ENDDO




  !
  ! --- neutral densities
  !
  IF(.NOT. ALLOCATED(enn))ALLOCATE(enn(nj,nneu))
  DO jn=1,nneu
     IF (irwflag .EQ. 0) THEN
        WRITE  (niterdb, 3057) namen(jn)
3057    FORMAT ('*  neutral density, #/meter**3, species: ', a)
        WRITE  (niterdb, 10) (enn(j,jn), j=1,nj)
     ELSE  
        READ   (niterdb,  7) line
        !               CALL   extract_name(line,namen(jn))
        READ   (niterdb, 10) (enn(j,jn)       , j=1,nj)
     END IF
  END DO
  !
  ! --- neutral density
  !
  IF(.NOT. ALLOCATED(ennw))ALLOCATE(ennw(nj,nneu))
  DO jn=1,nneu
     IF (irwflag .EQ. 0) THEN
        WRITE  (niterdb, 3200) namen(jn)
3200    FORMAT ('*  neutral density from wall source,', &
             ' #/meter**3, species: ', a)
        WRITE  (niterdb, 10) (ennw(j,jn), j=1,nj)
     ELSE 
        READ   (niterdb,  7) starflag
        READ   (niterdb, 10) (ennw(j,jn)       , j=1,nj)
     END IF
  END DO
  !
  ! --- neutral density
  !
  IF(.NOT. ALLOCATED(ennv))ALLOCATE(ennv(nj,nneu))
  DO jn=1,nneu
     IF (irwflag .EQ. 0) THEN
        WRITE  (niterdb, 3210) namen(jn)
3210    FORMAT ('*  neutral density from volume source,', &
             ' #/meter**3, species: ', a)
        WRITE  (niterdb, 10) (ennv(j,jn), j=1,nj)
     ELSE 
        READ   (niterdb,  7) line
        READ   (niterdb, 10) (ennv(j,jn)       , j=1,nj)
     END IF
  END DO
  !
  ! --- neutral source
  !
  IF(.NOT. ALLOCATED(volsn))ALLOCATE(volsn(nj,nneu))
  DO jn=1,nneu
     IF (irwflag .EQ. 0) THEN
        WRITE  (niterdb, 3420) namen(jn)
3420    FORMAT ('*  volume source of neutrals,', &
             ' #/(meter**3*second), species: ', a)
        WRITE  (niterdb, 10) (volsn(j,jn), j=1,nj)
     ELSE 
        READ   (niterdb,  7) line
        READ   (niterdb, 10) (volsn(j,jn)       , j=1,nj)
     END IF
  END DO



  !
  ! --- electron source due to beams
  !
  IF(.NOT. ALLOCATED(sbeame))ALLOCATE(sbeame(nj))
  IF (irwflag .EQ. 0) THEN
     WRITE  (niterdb, 3230)
3230 FORMAT ('*  sbeame : beam electron source,', &
          ' #/(meter**3*second)')
     WRITE  (niterdb, 10) (sbeame(j), j=1,nj)
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (sbeame(j)       , j=1,nj)
  END IF


  !
  ! --- thermal ion source due to beams
  !
  IF(nbion >  0)THEN
     IF(.NOT. ALLOCATED(sbeam))ALLOCATE(sbeam(nj,nbion))
     DO jn = 1,nbion
        label = '* sbeam : beam thermal ion source, #/(meter**3*second),species: '
        label = label(1:LEN_TRIM(label))//nameb(jn)
        IF (irwflag .EQ. 0) THEN 
           WRITE  (niterdb, FMT='(a)')label(1:LEN_TRIM(label))
           WRITE  (niterdb, 10) (sbeam(j,jn), j=1,nj)
        ELSE
           READ   (niterdb,  7) starflag 
           READ   (niterdb, 10) (sbeam(j,jn)       , j=1,nj)
        END IF
     ENDDO
  ENDIF



  ! 
  ! --- current density
  !
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 3058)
3058 FORMAT ('*  total current density, amps/meter**2')
     work(1:nj) = get_values(mhd_dat%curden)
     WRITE  (niterdb, 10) (work(j), j=1,nj)
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (work(j)      , j=1,nj)
     mhd_dat%curden = new_Vector(nj,work)
  END IF


  ! 
  ! --- ohmic current density
  !
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 3059)
3059 FORMAT ('*  ohmic current density, amps/meter**2')
     work(1:nj) = get_values(mhd_dat%curohm)
     WRITE  (niterdb, 10) (work(j), j=1,nj)
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (work(j)      , j=1,nj)
     mhd_dat%curohm = new_Vector(nj,work)
  END IF



  !
  ! --- bootstrap current density
  !
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 3060)
3060 FORMAT ('*  bootstrap current density, amps/meter**2')
     work(1:nj) = get_values(mhd_dat%curboot)
     WRITE  (niterdb, 10) (work(j), j=1,nj)
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (work(j)      , j=1,nj)
     mhd_dat%curboot = new_Vector(nj,work)
  END IF


  !
  ! --- beam current density
  ! --- no mhd_dat%curbeam because beam current is fixed in gcnm code
  !
  IF(.NOT. ALLOCATED(curbeam))ALLOCATE(curbeam(nj))
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 3061)
3061 FORMAT ('*  beam-driven current density,<Jbeam dot B/Bt0> amps/meter**2')
     WRITE  (niterdb, 10) (curbeam(j), j=1,nj)
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (curbeam(j)      , j=1,nj)
  END IF



  !
  ! --- RF current density
  !
  IF(.NOT. ALLOCATED(currf))ALLOCATE(currf(nj))
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 3070)
3070 FORMAT ('*  RF current density, amps/meter**2')
     WRITE  (niterdb, 10) (currf(j), j=1,nj)
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (currf(j)      , j=1,nj)
  END IF





  !
  ! --- toroidal electric field, profile
  !
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 4120)
4120 FORMAT ('*  toroidal electric field profile, V/m')
     work(1:nj) = get_values(profile%etor)
     WRITE  (niterdb, 10) (work(j), j=1,nj)
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (work(j), j=1,nj)
     profile%etor = new_Vector(nj,work)
  END IF


  !
  ! --- rho*bp0*fcap*gcap*hcap
  !
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 3080)
3080 FORMAT ('*  rho*bp0*fcap*gcap*hcap, tesla*meters') 
     work(1:nj) = get_values(mhd_dat%rbp)
     WRITE  (niterdb, 10) (work(j), j=1,nj)
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (work(j), j=1,nj)
     mhd_dat%rbp  = new_Vector(nj,work)
  END IF



  !
  ! --- zeff profile
  !
  IF(.NOT. ASSOCIATED(zeff))ALLOCATE(zeff(nj))
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 3110)
3110 FORMAT ('*  zeff profile')
     WRITE  (niterdb, 10) (zeff(j), j=1,nj)
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (zeff(j), j=1,nj)
     profile%zeff =  new_Vector(nj,zeff)
  END IF



  !
  ! --- angular rotation speed profile
  !
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 3120)
3120 FORMAT ('*  angular rotation speed profile, rad/sec')
     work(1:nj) = get_values(profile%angrot)
     WRITE  (niterdb, 10) (work(j), j=1,nj)
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (work(j), j=1,nj)
     profile%angrot = new_Vector(nj,work)
  END IF


  !
  ! --- total (nprim+nimp+te+ti+rbp+w) diffusivity matrix:
  ! --- dependent variable, UNITS :
  ! --- density,flux [1/(m^2 sec)]
  ! --- d(1:nion,1:nion,1:nj-1)   [m^2/sec] 
  ! --- d(1:nion,1:nion,nion+1:nion+2,1:nj-1) [1/kev m sec] 
  ! --- d(1:nion,1:nion,nion+3,1:nj-1) [1/tesla  m ^2sec] 
  ! --- d(1:nion,1:nion,nion+4,1:nj-1)  [1/m]

  ! --- electron and ion temperatures,flux [kev/(m^2 sec)]:
  ! --- d(nion+1:nion+2,1:nion,1:nj-1) [m^2 kev/sec] 
  ! --- d(nion+1:nion+2,nion+1:nion+2,1:nj-1) [1/(m sec)] 
  ! --- d(nion+1:nion+2,nion+3,1:nj-1) [(kev tesla)/(m^2 sec)] 
  ! --- d(nion+1:nion+2,nion+4,1:nj-1) [kev/m]

  ! --- r*poloidal b field, flux [Tesla/sec  ]
  ! --- d(nion+3,1:nion,1:nj-1) [Tesla m^4/sec]
  ! --- d(nion+3,nion+1:nion+2,1:nj-1) [Tesla m^4/sec] 
  ! --- d(nion+3,nion+3,1:nj-1)[1/sec] [1/sec]
  ! --- d(nion+3,nion+4,1:nj-1)[Tesla m ]


  ! --- toroidal rotation:
  ! --- d(nion+4,1:nion,1:nj-1) [kg m^4/sec^2]
  ! --- d(nion+4,nion+1:nion+2,1:nj-1) [(kg m)/(kev sec^2)] 
  ! --- d(nion+4,nion+3,1:nj-1) [(kg/(Tesla sec^2)] 
  ! --- d(nion+4,nion+4,1:nj)[(kg m)/sec] 

  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 3122)
3122 FORMAT ('* total diffusivity matrix , mixed units', &
          ' on half grid')
     DO j=1,ntot
        DO jj=1,ntot
           IF(j .LE. nion)THEN
              IF(  jj .LE. nion)THEN
                 WRITE(niterdb, FMT= '("*  d(1:nion,1:nion,1:nj) m^2/sec")' )
              ELSE IF( jj == nion+1 .OR. jj == nion+2)THEN
                 WRITE(niterdb, FMT= '("*  d(1:nion,nion+1:nion+2,1:nj) 1/kev m sec")' )
              ELSEIF(jj == nion+3)THEN
                 WRITE(niterdb, FMT= '("*  d(1:nion,nion+3,1:nj) (kev tesla)/(m^2 sec)")' )
              ELSE
                 WRITE(niterdb, FMT= '("*  d(1:nion,nion+4,1:nj) 1/m )")')
              ENDIF
           ELSEIF( j == nion+1 .OR.  j== nion+2)THEN
              IF(  jj .LE. nion)THEN
                 WRITE(niterdb, FMT= '("*  d(nion+1:nion+2,1:nion,1:nj) m^2 kev/sec )")')
              ELSE IF( jj == nion+1 .OR. jj == nion+2)THEN
                 WRITE(niterdb, FMT= '("*  d(nion+1:nion+2,nion+1:nion+2,1:nj)  1/(m sec) )")')
              ELSEIF(jj == nion+3)THEN
                 WRITE(niterdb, FMT= '("*  d(nion+1:nion+2,nion+3,1:nj) (kev tesla)/(m^2 sec) ")')
              ELSE
                 WRITE(niterdb, FMT= '("*  d(nion+1:nion+2,nion+4,1:nj) kev/m ")')
              ENDIF
           ELSEIF( j == nion+3)THEN
              IF(  jj .LE. nion)THEN
                 WRITE(niterdb, FMT= '("*  d(nion+3,1:nion,1:nj) Tesla m^4/sec )")')
              ELSE IF( jj == nion+1 .OR. jj == nion+2)THEN
                 WRITE(niterdb, FMT= '("*  d(nion+3,nion+1:nion+2,1:nj) Tesla m^4/sec )")')
              ELSEIF(jj == nion+3)THEN
                 WRITE(niterdb, FMT= '("*  d(nion+3,nion+3,1:nj) 1/sec")')
              ELSE
                 WRITE(niterdb, FMT= '("*  d(nion+3,nion+4,1:nj) Tesla m")')
              ENDIF
           ELSE
              IF(  jj .LE. nion)THEN
                 WRITE(niterdb, FMT= '("*  d(nion+4,1:nion,1:nj) kg m^4/sec^2 )")')
              ELSE IF( jj == nion+1 .OR. jj == nion+2)THEN
                 WRITE(niterdb, FMT= '("*  d(nion+4,nion+1:nion+2,1:nj) (kg m)/(kev sec^2) )")')
              ELSEIF(jj == nion+3)THEN
                 WRITE(niterdb, FMT= '("*  d(nion+4,nion+3,1:nj)(kg/(Tesla sec^2)  )")')
              ELSE
                 WRITE(niterdb, FMT= '("*  d(nion+4,nion+4,1:nj) (kg m)/sec )")')
              ENDIF
           ENDIF
           WRITE  (niterdb, 10) (diffuse%dcoef(j,jj,k), k=1,nj)
        ENDDO
     ENDDO
  ELSE
     IF(ASSOCIATED(diffuse%dcoef))DEALLOCATE(diffuse%dcoef)
     ALLOCATE(diffuse%dcoef(ntot,ntot,nj))
     READ(niterdb,7)starflag
     DO j=1,ntot
        DO jj=1,ntot
           READ   (niterdb,7) starflag
           READ   (niterdb, 10) (diffuse%dcoef(j,jj,k), k=1,nj)
        ENDDO
     ENDDO
  END IF

  !
  ! --- thermal diff. profiles, electron and ion
  !
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 3130)
3130 FORMAT ('*  electron thermal diffusivity, meters**2/sec', &
          ' on half grid')
     work(1:nj) = get_values(diffuse%chieinv)
     WRITE  (niterdb, 10) (work(j), j=1,nj)
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (work(j)       , j=1,nj)
     diffuse%chieinv = new_Vector(nj,work)
  END IF




  !
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 3140)
3140 FORMAT ('*  ion thermal diffusivity, meters**2/second', &
          ' on half grid')
     work(1:nj) = get_values(diffuse%chiinv)
     WRITE  (niterdb, 10) (work(j), j=1,nj)
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (work(j)       , j=1,nj)
     diffuse%chiinv = new_Vector(nj,work)
  END IF



  !
  ! --- ion neoclassical thermal conductivity
  !
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) WRITE (niterdb, 3145)
3145 FORMAT ('*  ion neoclassical thermal conductivity,', &
          ' 1/(meter*second), on half grid')
     work(1:nj) = get_values(diffuse%xkineo)
     WRITE  (niterdb, 10) (work(j), j=1,nj)
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (work(j)      , j=1,nj)
     diffuse%xkineo = new_Vector(nj,work)
  END IF

  !
  ! --- electron  neoclassical thermal conductivity
  !
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) WRITE (niterdb, 3146)
3146 FORMAT ('*  electron neoclassical thermal conductivity,', &
          ' 1/(meter*second), on half grid')
     work(1:nj) = get_values(diffuse%xkeneo)
     WRITE  (niterdb, 10) (work(j), j=1,nj)
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (work(j)      , j=1,nj)
     diffuse%xkeneo = new_Vector(nj,work)
  END IF

  !NOTE: total electron and ion diffusivity not read in. defined elsewhere
  ! diffuse%xkitot, diffuse%xketot ,... see plasma_properties, type diffusivity
  !(these are pointers NOT vectors):
  

  !
  ! --- d(electron pressure)/dt profile
  !
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 3150)
3150 FORMAT ('*  1.5*dpedt,power density due to electron pressure watts/(meter**3)')
     work(1:nj) = get_values(wpdot%dpedt)
     WRITE  (niterdb, 10) (work(j), j=1,nj)
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (work(j), j=1,nj)
     wpdot%dpedt = new_Vector(nj,work)
  END IF




  !
  ! --- electron conduction profile
  !
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 3170)
3170 FORMAT ('*  electron conduction, watts/meter**3')
     work(1:nj) = get_values(pwrden%qconde)
     WRITE  (niterdb, 10) (work(j), j=1,nj)
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (work(j)        , j=1,nj)
     pwrden%qconde = new_Vector(nj,work)
  END IF


  !
  ! --- ion conduction profile
  !
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 3180)
3180 FORMAT ('*  ion conduction, watts/meter**3')
     work(1:nj) = get_values(pwrden%qcondi)
     WRITE  (niterdb, 10) (work(j), j=1,nj)
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (work(j)        , j=1,nj)
     pwrden%qcondi = new_Vector(nj,work)
  END IF



  !
  ! --- electron convection profile
  !
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 3190)
3190 FORMAT ('*  electron convection, watts/meter**3')
     work(1:nj) = get_values(pwrden%qconve)
     WRITE  (niterdb, 10) (work(j), j=1,nj)
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (work(j)        , j=1,nj)
     pwrden%qconve = new_Vector(nj,work)
  END IF


  !
  ! --- ion convection profile
  !
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 3410)
3410 FORMAT ('*  ion convection, watts/meter**3')
     work(1:nj) = get_values(pwrden%qconvi)
     WRITE  (niterdb, 10) (work(j), j=1,nj)
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (work(j)        , j=1,nj)
     pwrden%qconvi = new_Vector(nj,work)
  END IF


  !
  ! --- beam electron heating profile
  !
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 3215)
3215 FORMAT ('*  power to elec. from beam, watts/meter**3')
     work(1:nj) = get_values(pwrden%qbeame)
     WRITE  (niterdb, 10) (work(j), j=1,nj)
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (work(j)        , j=1,nj)
     pwrden%qbeame = new_Vector(nj,work)
  END IF


  !
  ! --- electron ion energgy exchange  profile
  !
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 3220)
3220 FORMAT ('*  qdelt : electron-ion enrgy exchange,', &
          ' watts/meter**3')
     work(1:nj) = get_values(pwrden%qdelt)
     WRITE  (niterdb, 10) (work(j), j=1,nj)
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) ( work(j)        , j=1,nj)
     pwrden%qdelt = new_Vector(nj,work)
  END IF


  !
  ! --- beam ion heating profile
  !
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 3450)
3450 FORMAT ('*  power to ions from beam, watts/meter**3')
     work(1:nj) = get_values(pwrden%qbeami)
     WRITE  (niterdb, 10) (work(j), j=1,nj)
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (work(j)        , j=1,nj)
     pwrden%qbeami = new_Vector(nj,work)
  END IF


  !
  ! --- RF electron heating profile
  !
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 3240)
3240 FORMAT ('*  qrfe, RF electron heating, watts/meter**3')
     work(1:nj) = get_values(pwrden%qrfe)
     WRITE  (niterdb, 10) (work(j), j=1,nj)
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (work(j)        , j=1,nj)
     pwrden%qrfe = new_Vector(nj,work)
  END IF



  !
  ! --- RF ion heating profile
  !
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 3250)
3250 FORMAT ('*  qrfi,  RF ion heating, watts/meter**3')
     work(1:nj) = get_values(pwrden%qrfi)
     WRITE  (niterdb, 10) (work(j), j=1,nj)
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (work(j)        , j=1,nj)
     pwrden%qrfi = new_Vector(nj,work)
  END IF



  !
  ! --- qione heating profile
  !
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 3260)
3260 FORMAT ('*  qione, electron power density due to', &
          ' recombination and impact ionization,', &
          ' watts/meter**3')
     work(1:nj) = get_values(pwrden%qione)
     WRITE  (niterdb, 10) (work(j), j=1,nj)
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (work(j)        , j=1,nj)
     pwrden%qione = new_Vector(nj,work)
  END IF


  !
  ! --- qioni heating profile
  !
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 3270)
3270 FORMAT ('*  qioni, ion power density due to', &
          ' recombination and impact ionization,', &
          ' watts/meter**3')
     work(1:nj) = get_values(pwrden%qioni)
     WRITE  (niterdb, 10) (work(j), j=1,nj)
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (work(j)        , j=1,nj)
     pwrden%qioni = new_Vector(nj,work)
  END IF



  !
  ! --- qxc, ion heating profile
  !
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 3275)
3275 FORMAT ('*  qcx, ion power density due to', &
          ' neutral-thermal ion charge exchange, watts/meter**3')
     work(1:nj) = get_values(pwrden%qcx)
     WRITE  (niterdb, 10) (work(j), j=1,nj)
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (work(j)        , j=1,nj)
     pwrden%qcx   = new_Vector(nj,work)
  END IF


  !
  ! --- 2d electron heating profile
  !
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 3280)
3280 FORMAT ('*  2d electron heating, watts/meter**3')
     work(1:nj) = get_values(pwrden%qe2d)
     WRITE  (niterdb, 10) (work(j), j=1,nj)
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (work(j)        , j=1,nj)
     pwrden%qe2d   = new_Vector(nj,work)
  END IF


  !
  ! --- 2d ion heating profile
  !
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 3290)
3290 FORMAT ('*  2d ion heating, watts/meter**3')
     work(1:nj) = get_values(pwrden%qi2d)
     WRITE  (niterdb, 10) (work(j), j=1,nj)
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (work(j)        , j=1,nj)
     pwrden%qi2d   = new_Vector(nj,work)
  END IF


  !
  ! --- fusion electron heating  profile
  !      qfuse = qtfuse + qbfuse 
  !      (sum of thermal and fast ion progenated fusion heating to electrons
  !       fast part includes beam and alphas)
  !      (since qbfuse and qfuse are input, qtfuse is implicitely determined)
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 3300)
3300 FORMAT ('* total  fusion electron heating, watts/meter**3')
     work(1:nj) = get_values(pwrden%qfuse)
     WRITE  (niterdb, 10) (work(j), j=1,nj)
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (work(j)        , j=1,nj)
     pwrden%qfuse  = new_Vector(nj,work)
  END IF


  !
  ! --- fusion ion heating profile
  !      qfusi = qtfusi + qbfusi 
  !      (sum of thermal and fast ion progenated fusion heating to ions
  !       fast part includes beam and alphas)
  !      (since qbfusi and qfusi are input, qtfusi is implicitely determined)
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 3310)
3310 FORMAT ('* total  fusion ion heating, watts/meter**3')
     work(1:nj) = get_values(pwrden%qfusi)
     WRITE  (niterdb, 10) (work(j), j=1,nj)
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (work(j)        , j=1,nj)
     pwrden%qfusi  = new_Vector(nj,work)
  END IF


  !
  ! --- beam fusion electron heating profile
  ! --- (fraction of beam fusion energy deposited on thermal electron distribution
  !
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 3312)
3312 FORMAT ('*  beam fusion electron heating, watts/meter**3')
     work(1:nj) = get_values(pwrden%qbfuse)
     WRITE  (niterdb, 10) (work(j), j=1,nj)
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (work(j)        , j=1,nj)
     pwrden%qbfuse  = new_Vector(nj,work)
  END IF


  !
  ! --- beam fusion ion heating profile
  ! --- (fraction of beam fusion energy deposited on thermal ion distribution)
  !
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 3315)
3315 FORMAT ('*  beam fusion ion heating, watts/meter**3')
     work(1:nj) = get_values(pwrden%qbfusi)
     WRITE  (niterdb, 10) (work(j), j=1,nj)
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (work(j)        , j=1,nj)
     pwrden%qbfusi  = new_Vector(nj,work)
  END IF


  !
  ! --- mag electron heating profile,change in magnetic energy associated with 
  !     mhd mixing (qmag is normally added to the electron heating term)
  !
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 3320)
3320 FORMAT ('*  qmag electron heating, watts/meter**3')
     work(1:nj) = get_values(pwrden%qmag)
     WRITE  (niterdb, 10) (work(j), j=1,nj)
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (work(j)        , j=1,nj)
     pwrden%qmag  = new_Vector(nj,work)
  END IF

  !
  ! --- qexch anomalous electron ion energy exchange
  !
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 3325)
3325 FORMAT ('*  qexch anomlaous electron-ion energy exchange, watts/meter**3')
     WRITE  (niterdb, 10) (pwrden%qexch%data(j), j=1,nj)
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (work(j), j=1,nj)
     pwrden%qexch   = new_Vector(nj,work)
  END IF

  !
  ! --- sawtooth electron heating profile
  !
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 3330)
3330 FORMAT ('*  sawtooth electron heating, watts/meter**3')
     work(1:nj) = get_values(pwrden%qsawe)
     WRITE  (niterdb, 10) (work(j), j=1,nj)
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (work(j)        , j=1,nj)
     pwrden%qsawe  = new_Vector(nj,work)
  END IF


  !
  ! --- sawtooth ion heating profile
  !
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 3340)
3340 FORMAT ('*  sawtooth ion  heating, watts/meter**3')
     work(1:nj) = get_values(pwrden%qsawi)
     WRITE  (niterdb, 10) (work(j),j=1,nj)
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (work(j)        , j=1,nj)
     pwrden%qsawi  = new_Vector(nj,work)
  END IF


  !
  ! --- radiated power density
  !

  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 3345)
3345 FORMAT ('*  radiated power density, watts/meter**3')
     work(1:nj) = get_values(pwrden%qrad)
     WRITE  (niterdb, 10) (-ABS(work(j)), j=1,nj)
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (work(j)         , j=1,nj)
     pwrden%qrad  = new_Vector(nj,work)
  END IF


  !
  ! --- qohm,ohmic heating profile
  !
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 3350)
3350 FORMAT ('*  (electron) ohmic power density, watts/meter**3')
     work(1:nj) = get_values(pwrden%qohm)
     WRITE  (niterdb, 10) (work(j), j=1,nj)
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (work(j)        , j=1,nj)
     pwrden%qohm  = new_Vector(nj,work)
  END IF


  !
  ! --- convert vectors defined on the npsi (i.e., MHD) grid to corresponding
  ! --- quantities on the rho (i.e., transport) grid
  !
  !      tension = 0.0               ! don't use tension option of tspline
  !      tmax    = 0.0               ! max allowed tension
  !      bpar(1) = 0.0               ! set boundary conditions on spline
  !      bpar(2) = 0.0
  !      bpar(3) = 0.0
  !      bpar(4) = 0.0
  !
  ! --- take care of a roundoff problem
  ! --- changed to relative measure HSJ 8/25/98
  !
  !      IF  (psir(nj) .NE. 0.0) THEN
  !        IF (ABS((psival(1)-psir(nj))/(psir(nj))) .LT. 1.0e-10) &
  !                                           psival(1) = psir(nj)
  !      ELSE
  !        IF (ABS(psival(1)-psir(nj)) .LT. 1.0e-10) &
  !                                           psival(1) = psir(nj)
  !      END IF








  !
  ! --- average major radius
  !
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 1100)
1100 FORMAT ('*  average major radius of each flux surface,', &
          ' meters, evaluated at elevation of magnetic axis')
     WRITE  (niterdb, 10) (dischg%rmajavnpsi%data(j),j=1,mhd_dat%npsi)
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (work(j), j=1,mhd_dat%npsi)
     dischg%rmajavnpsi = new_Vector(mhd_dat%npsi,work)
  END IF




  !
  ! --- average minor radius
  !
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 1110)
1110 FORMAT ('*  average minor radius of each flux surface, ', &
          'meters, evaluated at elevation of magnetic axis')
     WRITE (niterdb, 10) (dischg%rminavnpsi%data(j), j=1,mhd_dat%npsi)
  ELSE
     READ  (niterdb,  7) starflag
     READ  (niterdb, 10) (work(j), j=1,mhd_dat%npsi)
     dischg%rminavnpsi = new_Vector(mhd_dat%npsi,work)
  END IF


  !
  ! --- volume
  !
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 1120)
1120 FORMAT ('*  volume of each flux surface, meters**3')
     WRITE (niterdb, 10) (dischg%psivolpnpsi%data(j), j=1,mhd_dat%npsi)
  ELSE
     READ  (niterdb,  7) starflag
     READ   (niterdb, 10) (work(j), j=1,mhd_dat%npsi)
     dischg%psivolpnpsi = new_Vector(mhd_dat%npsi,work)
  END IF


  !
  ! ---  elongation
  !
  IF (irwflag .EQ. 0) THEN
     !      dissable the following for NTCC veersion
     !             CALL tspline (psival,elongx,npsi,bpar,cs2spline,kpsi, &
     !                           ier,tension,aspline,bspline,cspline,dspline, &
     !                           espline,fspline,tmax,psir,work,nj)
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 1130)
1130 FORMAT ('*  elongation of each flux surface')
     WRITE  (niterdb, 10) (dischg%elongxnpsi%data(j)  , j=1,mhd_dat%npsi)
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (work(j), j=1,mhd_dat%npsi)
     dischg%elongxnpsi = new_Vector(mhd_dat%npsi,work)
  END IF






  !
  ! --- triangularity
  !
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 1140)
1140 FORMAT ('* upper   triangularity of each flux surface')
     WRITE  (niterdb, 10) (dischg%triangnpsi_u%data(j), j=1,mhd_dat%npsi)

     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 1141)
1141 FORMAT ('* lower triangularity of each flux surface')
     WRITE  (niterdb, 10) (dischg%triangnpsi_l%data(j), j=1,mhd_dat%npsi)

  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (work(j), j=1,mhd_dat%npsi)
     dischg%triangnpsi_u = new_Vector(mhd_dat%npsi,work)
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (work(j), j=1,mhd_dat%npsi)
     dischg%triangnpsi_l = new_Vector(mhd_dat%npsi,work)
  END IF


  !
  ! --- indentation
  !
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 1150)
1150 FORMAT ('*  indentation of each flux surface')
     WRITE  (niterdb, 10) (dischg%pindentnpsi%data(j), j=1,mhd_dat%npsi)
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (work(j), j=1,mhd_dat%npsi)
     dischg%pindentnpsi  = new_Vector(mhd_dat%npsi,work)
  END IF


  !
  ! --- surface area
  !
  IF (irwflag .EQ. 0) THEN

     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 1160)
1160 FORMAT ('*  surface area of each flux surface, meter**2' / &
          '*  this is 4*pi*pi*R0*hcap*rho*<ABS(grad rho)>;', &
          ' on mhd grid')

     WRITE (niterdb, 10) (dischg%sfareanpsi%data(j), j=1,mhd_dat%npsi)
  ELSE
     READ  (niterdb,  7)  starflag
     READ  (niterdb,  7)  starflag
     READ  (niterdb, 10) (work(j), j=1,mhd_dat%npsi)
     dischg%sfareanpsi = new_Vector(mhd_dat%npsi,work)
  END IF


  !
  ! --- cross-sectional area
  !
  IF (irwflag .EQ. 0) THEN

     IF (iterdsc .NE. 0) &
          
          WRITE  (niterdb, 1165)
1165 FORMAT ('*  cross-sectional area of each flux', &
          ' on mhd grid surface, meters**2')

     WRITE (niterdb, 10) (dischg%cxareanpsi%data(j), j=1,mhd_dat%npsi)
  ELSE
     READ  (niterdb,  7) starflag
     READ  (niterdb, 10) (work(j), j=1,mhd_dat%npsi)
     dischg%cxareanpsi = new_Vector(mhd_dat%npsi,work)
  END IF



  !
  ! --- flux surface average grad rho
  !
  IF (irwflag .EQ. 0) THEN

     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 1170) dimensionless
1170 FORMAT ('*  flux surface average absolute grad rho on mhd grid, ',a)

     WRITE  (niterdb, 10) (dischg%grho1npsi%data(j), j=1,mhd_dat%npsi)
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (work(j), j=1,mhd_dat%npsi)
     dischg%grho1npsi  = new_Vector(mhd_dat%npsi,work)
  END IF





  !
  ! --- flux surface average (grad rho)**2
  !
  IF (irwflag .EQ. 0) THEN

     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 1180)
1180 FORMAT ('*  flux surface average (grad rho)**2 on mhd grid ')
     WRITE  (niterdb, 10) (dischg%grho2npsi%data(j), j=1,mhd_dat%npsi)
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (work(j), j=1,mhd_dat%npsi)
     dischg%grho2npsi = new_Vector(mhd_dat%npsi,work)
  END IF



  !
  ! --- plasma boundary
  !

  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 1200)
1200 FORMAT ('*  nplasdry : number of points on plasma boundary')
     WRITE  (niterdb, 9) dischg%nplasbdry
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 1210)
1210 FORMAT ('*  r points for plasma boundary, meters')
     IF(SIZE(work) .LT. dischg%nplasbdry)THEN
        IF(ALLOCATED(work))DEALLOCATE(work)
        ALLOCATE(work(dischg%nplasbdry))
     ENDIF
     work(1:dischg%nplasbdry) = get_values(dischg%rplasbdry)
     WRITE  (niterdb, 10) (work(j), j=1,dischg%nplasbdry)
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 1215)
1215 FORMAT ('*  z points for plasma boundary, meters')
     work(1:dischg%nplasbdry) = get_values(dischg%zplasbdry)
     WRITE  (niterdb, 10) (work(j), j=1,dischg%nplasbdry)
  ELSE
     READ   (niterdb,  7)  starflag
     READ   (niterdb,  9)  dischg%nplasbdry

     DEALLOCATE(work)
     ALLOCATE(work(dischg%nplasbdry))

     READ   (niterdb,  7)  starflag
     READ   (niterdb, 10) (work(j), j=1,dischg%nplasbdry)
     dischg%rplasbdry  = new_Vector(dischg%nplasbdry,work)
     READ   (niterdb,  7)  starflag
     READ   (niterdb, 10) (work(j), j=1,dischg%nplasbdry)
     dischg%zplasbdry  = new_Vector(dischg%nplasbdry,work)
  END IF



   IF(ALLOCATED(work))THEN
     DEALLOCATE(work)
     ALLOCATE(work(nj)) ! assumes all usage of work below is of size nj
   ENDIF


  !
  ! --- beam torque density - new input added here so that backward
  !     compatibility is maintained
  !
  IF(.NOT. ALLOCATED(storqueb))ALLOCATE(storqueb(nj))
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 3125)
3125 FORMAT ('* beam   torque density, nt-m/m**3')
     WRITE  (niterdb, 10) (storqueb(j), j=1,nj) 

  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (storqueb(j), j=1,nj)

  END IF



  IF(.NOT. ALLOCATED(te_bc))ALLOCATE(te_bc(nj))
  IF(.NOT. ALLOCATED(ti_bc))ALLOCATE(ti_bc(nj))
  IF(.NOT. ALLOCATED(ene_bc))ALLOCATE(ene_bc(nj))
  IF(.NOT. ALLOCATED(zeff_bc))ALLOCATE(zeff_bc(nj))
  IF(.NOT. ALLOCATED(angrot_bc))ALLOCATE(angrot_bc(nj))
  IF(.NOT. ALLOCATED(en_bc))ALLOCATE(en_bc(nj,nion))
  IF(.NOT. ALLOCATED(flux_bc))ALLOCATE(flux_bc(nj,ntot))
  IF(.NOT. ALLOCATED(wbeam))ALLOCATE(wbeam(nj))
  IF(.NOT. ALLOCATED(walp))ALLOCATE(walp(nj))
  IF(.NOT. ALLOCATED(w_alpha))ALLOCATE(w_alpha(nj))
  IF(.NOT. ALLOCATED(enalp))ALLOCATE(enalp(nj))
  IF(.NOT. ASSOCIATED(zsq))ALLOCATE(zsq(nj,nion))
  IF(.NOT. ASSOCIATED(z))ALLOCATE(z(nj,nion))
  IF(.NOT. ASSOCIATED(wpdot%dpidt))ALLOCATE(wpdot%dpidt(nion))


  !     in the future we may introduce time dependent metrics.
  !     then dfdt (= d/dt FCAP ,etc) will be calculated internally in this code
  !     by reading in fcap_bc, etc.
  !     in demo code just assume time independent for each
  !     slice that the solver is called.
  IF(.NOT. ALLOCATED(dfdt))THEN 
     ALLOCATE(dfdt(nj))
     dfdt(:) = 0.0_DP
  ENDIF
  IF(.NOT. ALLOCATED(dgdt))THEN 
     ALLOCATE(dgdt(nj))
     dgdt(:) = 0.0_DP
  ENDIF
  IF(.NOT. ALLOCATED(dhdt))THEN
     ALLOCATE(dhdt(nj))
     dhdt(:) = 0.0_DP
  ENDIF


  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 3640)
3640 FORMAT ('* suggested  end of simulation time for next cycle,sec')
     WRITE  (niterdb, 10)tGCNMf
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) tGCNMf 
  END IF



  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 3642)
3642 FORMAT ('* boundary condition time ,sec')
     WRITE  (niterdb, 10)time_bc
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) time_bc 
  END IF


  WRITE(label,FMT='(1pe16.8)')time_bc
  READ(label,FMT='(a)')bc_asc_time
  bc_asc_time = ADJUSTL(bc_asc_time)


  IF (irwflag .EQ. 0) THEN
     label = '* totcur_bc: boundary condition total current,amps,at time '//bc_asc_time(1:LEN_TRIM(bc_asc_time))
     IF (iterdsc .NE. 0) WRITE(niterdb,FMT='(A)')label(1:LEN_TRIM(label))
     WRITE  (niterdb, 10)totcur_bc
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) totcur_bc
  END IF



  IF (irwflag .EQ. 0) THEN
     label = '* boundary condition loop voltage,Volts,at time '//bc_asc_time(1:LEN_TRIM(bc_asc_time))
     IF (iterdsc .NE. 0) WRITE(niterdb,FMT='(A)')label(1:LEN_TRIM(label))
     WRITE  (niterdb, 10)vloop_bc
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) vloop_bc
  END IF

 
  IF(.NOT. ALLOCATED(fix_edge_ni_bc))ALLOCATE(fix_edge_ni_bc(nion))
 
  IF (irwflag .EQ. 0) THEN
     label ='* rho boundary(te,ti,rot,ni1,ni2..) flags at time '//bc_asc_time(1:LEN_TRIM(bc_asc_time))
     IF (iterdsc .NE. 0) WRITE(niterdb,FMT='(A)')label(1:LEN_TRIM(label))
     WRITE  (niterdb, 10)fix_edge_te_bc,fix_edge_ti_bc,             &
          fix_edge_rot_bc,fix_edge_ni_bc(1:nion) 
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) fix_edge_te_bc,fix_edge_ti_bc,             &
          fix_edge_rot_bc,fix_edge_ni_bc(1:nion) 
  ENDIF
 !convert grid point based fix_edge values to normalized rho point values:

  CALL convert_edge_indecies

  !      te_bc(:) = 0.0
  !      te_bc(fix_edge_te_bc:nj) = bctime_zone(fix_edge_te_bc:nj,kion+1)
  IF (irwflag .EQ. 0) THEN
     label = '* boundary condition TE,kev,at time '//bc_asc_time(1:LEN_TRIM(bc_asc_time))
     IF (iterdsc .NE. 0) WRITE(niterdb,FMT='(A)')label(1:LEN_TRIM(label))
     WRITE  (niterdb, 10) (te_bc(j), j=1,nj) 
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (te_bc(j), j=1,nj)
  END IF




  !      ti_bc(:) = 0.0
  !      ti_bc(fix_edge_ti_bc:nj) = bctime_zone(fix_edge_ti_bc:nj,kion+2)
  IF (irwflag .EQ. 0) THEN
     label = '* boundary condition TI,kev,at time '//bc_asc_time(1:LEN_TRIM(bc_asc_time))
     IF (iterdsc .NE. 0) WRITE(niterdb,FMT='(A)')label(1:LEN_TRIM(label))
     WRITE  (niterdb, 10) (ti_bc(j), j=1,nj) 
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (ti_bc(j), j=1,nj)
  END IF




  IF (irwflag .EQ. 0) THEN
     label ='* bc profile: ene, 1/m**3 at time  '//bc_asc_time(1:LEN_TRIM(bc_asc_time))
     IF (iterdsc .NE. 0) WRITE(niterdb,FMT='(A)')label(1:LEN_TRIM(label))
     WRITE  (niterdb, 10) (ene_bc(j), j=1,nj) 
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (ene_bc(j), j=1,nj)
  END IF





  IF (irwflag .EQ. 0) THEN
     label = '* bc profile: zeff, at time  '//bc_asc_time(1:LEN_TRIM(bc_asc_time))
     IF (iterdsc .NE. 0) WRITE(niterdb,FMT='(A)')label(1:LEN_TRIM(label))
     WRITE  (niterdb, 10) (zeff_bc(j), j=1,nj) 
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (zeff_bc(j), j=1,nj)
  END IF





  !      angrot_bc(:) = 0.0
  !      angrot_bc(fix_edge_rot_bc:nj) = bctime_zone(fix_edge_rot_bc:nj,kion+dp4)
  IF (irwflag .EQ. 0) THEN
     label = '*  angular rotation speed profile, rad/sec at time '//bc_asc_time(1:LEN_TRIM(bc_asc_time))
     IF (iterdsc .NE. 0) WRITE(niterdb,FMT='(A)')label(1:LEN_TRIM(label))
     WRITE  (niterdb, 10) (angrot_bc(j), j=1,nj) 
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (angrot_bc(j), j=1,nj)
  END IF




  !
  ! --- primary, impurity densities boundary condition values
  ! --- the whole profile is specified but only values from rho=1.0
  ! --- down to some lower value of rho, determined by fix boundary parameters, 
  ! --- is used.
  jp = 0
  ji = 0
  DO jj=1,nion
     IF (jj .LE. nprim)  jp = jp + 1
     IF (jj .GT. nprim)  ji = ji + 1
     IF (irwflag .EQ. 0) THEN

        IF (iterdsc .NE. 0 .AND. jj .LE. nprim) THEN
           label = '* ion density, #/m^3 at time '//bc_asc_time(1:LEN_TRIM(bc_asc_time))//'species : '//namep(jp)
           WRITE(niterdb,FMT='(A)')label(1:LEN_TRIM(label))
        ENDIF
        IF (iterdsc .NE. 0 .AND. jj .GT. nprim)THEN
           label = '* ion density, #/m^3 at time '//bc_asc_time(1:LEN_TRIM(bc_asc_time))//'species : '//namei(ji)
           WRITE(niterdb,FMT='(A)')label(1:LEN_TRIM(label))
        ENDIF
        WRITE  (niterdb, 10) (en_bc(j,jj), j=1,nj)
     ELSE 
        READ   (niterdb,  7) line
        READ   (niterdb, 10) (en_bc(j,jj), j=1,nj)
     END IF
  END DO



  !      ditto for flux_bc 
  jp = 0
  ji = 0
  DO jj=1,ntot
     IF (jj .LE. nprim)  jp = jp + 1
     IF (jj .GT. nprim)  ji = ji + 1
     IF (irwflag .EQ. 0) THEN

        IF (iterdsc .NE. 0 .AND. jj .LE. nprim) &
             label = '* ion flux, #/(m^2 sec) at time '//bc_asc_time(1:LEN_TRIM(bc_asc_time))//'species : '//namep(jp)
        IF (iterdsc .NE. 0 .AND. jj .GT. nprim .AND. jj .LE. nion) &
             label = '* ion flux, #/(m^2 sec) at time '//bc_asc_time(1:LEN_TRIM(bc_asc_time))//'species : '//namei(ji)
        IF (iterdsc .NE. 0 .AND. jj ==  nion+1 )    &
             label = '* electron energy  flux, kev/(m^2 sec) at time '//bc_asc_time(1:LEN_TRIM(bc_asc_time))
        IF (iterdsc .NE. 0 .AND. jj ==  nion+2 )    &
             label = '* total ion  energy  flux, kev/(m^2 sec) at time '//bc_asc_time(1:LEN_TRIM(bc_asc_time))
        IF (iterdsc .NE. 0 .AND. jj ==  nion+3 )    &
             label = '* effective flux in Faradys law Tesla/(sec) at time '//bc_asc_time(1:LEN_TRIM(bc_asc_time))
        IF (iterdsc .NE. 0 .AND. jj ==  nion+4 )    &
             label = '* momentum flux  kg/(sec^2) at time '//bc_asc_time(1:LEN_TRIM(bc_asc_time))
        WRITE(niterdb,FMT='(A)')label(1:LEN_TRIM(label))

        WRITE  (niterdb, 10) (flux_bc(j,jj), j=1,nj)
     ELSE 
        READ   (niterdb,  7) line
        READ   (niterdb, 10) (flux_bc(j,jj), j=1,nj)
     END IF
  END DO
 




  !
  ! --- primary, impurity charge,charge square
  !
  jp = 0
  ji = 0
  DO jj=1,nion
     IF (jj .LE. nprim)  jp = jp + 1
     IF (jj .GT. nprim)  ji = ji + 1
     IF (irwflag .EQ. 0) THEN
        label ='* charge  z '
        IF (iterdsc .NE. 0 .AND. jj .LE. nprim)THEN
           label = '* charge  z,species : '//namep(jp)
        ELSE IF (iterdsc .NE. 0 .AND. jj .GT. nprim)THEN
           label = '* charge  z ,species : '//namei(ji)
        ENDIF
        WRITE  (niterdb, FMT='(A)')label(1:LEN_TRIM(label))
        WRITE  (niterdb, 10) (z(j,jj), j=1,nj)
        label ='* charge square,  zsq '
        IF (iterdsc .NE. 0 .AND. jj .LE. nprim)THEN
           label = '* charge square, zsq,species :  '//namep(jp)
        ELSE IF (iterdsc .NE. 0 .AND. jj .GT. nprim)THEN
           label = '* charge square,zsq ,species :   '//namei(ji)
        ENDIF
        WRITE  (niterdb, FMT='(A)')label(1:LEN_TRIM(label))
        WRITE  (niterdb, 10) (zsq(j,jj), j=1,nj)
     ELSE 
        READ   (niterdb,  7) starflag
        READ   (niterdb, 10) (z(j,jj), j=1,nj)
        profile%z(jj)    = zero_Vector(nj)         ! ifort is fussy about this
        profile%z(jj)%data(:) = z(:,jj)
        READ   (niterdb,  7) starflag
        READ   (niterdb, 10) (zsq(j,jj), j=1,nj)
        profile%zsq(jj)  = zero_Vector(nj)
        profile%zsq(jj)%data(:)  = zsq(:,jj)
     END IF
  END DO







  !
  ! --- fast ion stored energy density:
  !
  IF (irwflag .EQ. 0) THEN
     label = '* fast ion stored energy density KEV/m**3 ' 
     IF (iterdsc .NE. 0) WRITE(niterdb,FMT='(A)')label(1:LEN_TRIM(label))
     WRITE  (niterdb, 10) (wbeam(j), j=1,nj) 
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (wbeam(j), j=1,nj)
     profile%wbeam = new_Vector(nj,wbeam)
  END IF





  !
  ! --- fast alpha stored energy density:
  !
  IF (irwflag .EQ. 0) THEN
     label = '* fast alpha stored energy density KEV/m**3 ' 
     IF (iterdsc .NE. 0) WRITE(niterdb,FMT='(A)')label(1:LEN_TRIM(label))
     work(1:nj) = get_values(profile%walp)
     WHERE(work < zeroc)work = zeroc
     WRITE  (niterdb, 10) (work(j), j=1,nj) 
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (walp(j), j=1,nj)  ! walp used internally as initial condition for w_alpha
     WHERE(walp < zeroc) walp=zeroc          ! interpolation errors
     profile%walp   = new_Vector(nj,walp)    ! profile%walp used for io only
     
  END IF

  IF (irwflag .EQ. 0) THEN
     label = '* fast alpha density 1/m**3 ' 
     IF (iterdsc .NE. 0) WRITE(niterdb,FMT='(A)')label(1:LEN_TRIM(label))
     WRITE  (niterdb, 10) (enalp(j), j=1,nj) 
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (enalp(j), j=1,nj)
  END IF




  IF (irwflag .EQ. 0) THEN
     WRITE(label,FMT='(1pe16.8)')tGCNMs
     READ(label,FMT='(a)')st_asc_time
     st_asc_time = ADJUSTL(st_asc_time)
     label = '* rate of change of ene  1/(m**3 sec) at time ' &
          //st_asc_time(1:LEN_TRIM(st_asc_time))
     IF (iterdsc .NE. 0) WRITE(niterdb,FMT='(A)')label(1:LEN_TRIM(label))
     work(1:nj) = get_values(wpdot%dnedt)
     WRITE  (niterdb, 10) (work(j), j=1,nj) 
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (work(j), j=1,nj)
     wpdot%dnedt = new_Vector(nj,work)
  END IF


  !
  ! --- primary, impurity density derivatives
  !

!!$   WRITE(label,FMT='(1pe16.8)')tGCNMs
!!$   READ(label,FMT='(a)')st_asc_time
!!$   st_asc_time = ADJUSTL(st_asc_time)
!!$   label = '* rate of change of ion density   1/(m**3 sec) at time ' &
!!$       //st_asc_time(1:LEN_TRIM(st_asc_time))//' for species '
!!$   label = ADJUSTL(label)
!!$   jp = 0
!!$   ji = 0
!!$   DO jj=1,nion
!!$     IF (jj .LE. nprim)  jp = jp + 1
!!$     IF (jj .GT. nprim)  ji = ji + 1
!!$     IF (irwflag .EQ. 0) THEN
!!$        IF(jj .LE. nprim)THEN
!!$           tlabel = label(1:LEN_TRIM(label))//namep(jp)
!!$        ELSE
!!$           tlabel = label(1:LEN_TRIM(label))//namei(ji)
!!$        ENDIF
!!$        tlabel = ADJUSTL(tlabel)
!!$        IF (iterdsc .NE. 0 )&
!!$             WRITE  (niterdb, FMT='(A)') tlabel(1:LEN_TRIM(tlabel))
!!$        WRITE  (niterdb, 10) (dnidt(j,jj), j=1,nj)
!!$     ELSE
!!$        READ   (niterdb,  7) starflag
!!$        READ   (niterdb, 10) (dnidt(j,jj), j=1,nj)
!!$     END IF
!!$  END DO


  !
  ! --- primary, impurity pressure  derivatives
  !

  WRITE(label,FMT='(1pe16.8)')tGCNMs
  READ(label,FMT='(a)')st_asc_time
  st_asc_time = ADJUSTL(st_asc_time)
  label = '* 1.5*dpidt,power density due to  ion pressure  watts/(m**3) at time ' &
       //st_asc_time(1:LEN_TRIM(st_asc_time))//' for species : '
  label = ADJUSTL(label)
  jp = 0
  ji = 0
  DO jj=1,nion
     IF (jj .LE. nprim)  jp = jp + 1
     IF (jj .GT. nprim)  ji = ji + 1
     IF (irwflag .EQ. 0) THEN
        IF(jj .LE. nprim)THEN
           tlabel = label(1:LEN_TRIM(label))//namep(jp)
        ELSE
           tlabel = label(1:LEN_TRIM(label))//namei(ji)
        ENDIF
        tlabel = ADJUSTL(tlabel)
        IF (iterdsc .NE. 0 )&
             WRITE  (niterdb, FMT='(A)') tlabel(1:LEN_TRIM(tlabel))
        work(1:nj) = wpdot%dpidt(jj)%data(1:nj)
        WRITE  (niterdb, 10) (work(j), j=1,nj)
     ELSE
        READ   (niterdb,  7) starflag
        READ   (niterdb, 10) (work(j), j=1,nj)
        wpdot%dpidt(jj) = new_Vector(nj,work)
     END IF
  END DO






  IF(.NOT. ALLOCATED(eps))ALLOCATE(eps(nj))
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 3126)
3126 FORMAT ('*eps  horizontal inverse aspect ratio = (rmax-rmin)/(rmax+rmin)')
     WRITE  (niterdb, 10) (eps(j), j=1,nj) 
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (eps(j), j=1,nj)
  END IF



  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 3660)
3660 FORMAT ('*rcap = < R>,m')
     WRITE  (niterdb, 10) (mhd_dat%rcap%data(j), j=1,nj) 
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (work(j), j=1,nj)
     mhd_dat%rcap = new_Vector(nj,work)
  END IF

 IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 3662)
3662 FORMAT ('*rcapi = < 1/R>,1./m')
     WRITE  (niterdb, 10) (mhd_dat%rcapi%data(j), j=1,nj) 
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (work(j), j=1,nj)
     mhd_dat%rcapi = new_Vector(nj,work)
  END IF


  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 3670)
3670 FORMAT ('*r2capi = <R**2>       M**2')
     WRITE  (niterdb, 10) (mhd_dat%r2capi%data(j), j=1,nj) 
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (work(j), j=1,nj)
     mhd_dat%r2capi = new_Vector(nj,work)
  END IF




  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 3672)
3672 FORMAT ('*r2cap = <R0**2/R**2> ')
     WRITE  (niterdb, 10) (mhd_dat%r2cap%data(j), j=1,nj) 
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (work(j), j=1,nj)
     mhd_dat%r2cap = new_Vector(nj,work)
  END IF
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




  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 3680)
3680 FORMAT ('*xhm2 = < (B total/ B axis)**2 > (=1 for circular plasmas)')
     WRITE  (niterdb, 10) (xhm2(j), j=1,nj) 
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (xhm2(j), j=1,nj)
  END IF




  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 3690)
3690 FORMAT ('*xi11 ( = 1.95 sqrt(eps)for circular plasmas)')
     WRITE  (niterdb, 10) (xi11(j), j=1,nj) 
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (xi11(j), j=1,nj)
  END IF


  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 3700)
3700 FORMAT ('*xi33 ( = 1.95 sqrt(eps)for circular plasmas)')
     WRITE  (niterdb, 10) (xi33(j), j=1,nj) 
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (xi33(j), j=1,nj)
  END IF


  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 3710)
3710 FORMAT ('*xips = <(Baxis/B)**2)> - 1./(<(B/Baxis)**2> )( = 2 eps**2 for circular plasmas)')
     WRITE  (niterdb, 10) (xips(j), j=1,nj) 
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (xips(j), j=1,nj)
  END IF

  xhm20(:) = xhm2(:)
  xi110(:) = xi11(:)
  xi330(:) = xi33(:)
  xips0(:) = xips(:)


  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 3730)
3730 FORMAT ('*dfdt rate of change of fcap,1/sec')
     WRITE  (niterdb, 10) (dfdt(j), j=1,nj) 
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (dfdt(j), j=1,nj)
  END IF


  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 3740)
3740 FORMAT ('*dgdt rate of change of gcap,1/sec')
     WRITE  (niterdb, 10) (dgdt(j), j=1,nj) 
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (dgdt(j), j=1,nj)
  END IF

  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 3750)
3750 FORMAT ('*dhdt rate of change of hcap,1/sec')
     WRITE  (niterdb, 10) (dhdt(j), j=1,nj) 
  ELSE
     READ   (niterdb,  7) starflag
     READ   (niterdb, 10) (dhdt(j), j=1,nj)
  END IF


! added 9/10/10 
 IF(.NOT. ASSOCIATED(diffuse%eta%data))diffuse%eta     = zero_vector(nj) 
 IF(.NOT. ASSOCIATED(diffuse%ftrap%data))diffuse%ftrap = zero_vector(nj)
 IF(.NOT. ASSOCIATED(diffuse%xnuse%data))diffuse%xnuse = zero_vector(nj)
 IF(.NOT. ASSOCIATED(diffuse%xnus))THEN
        ALLOCATE(diffuse%xnus(nion))
        DO jj =1,nion
           diffuse%xnus(jj) = zero_Vector(nj)
        ENDDO
  ENDIF
  IF(.NOT. ASSOCIATED(fus_prod%neutr_ddn_th%data))fus_prod%neutr_ddn_th = zero_vector(nj)
  IF(.NOT. ASSOCIATED(fus_prod%neutr_ddn_beam_beam%data))fus_prod%neutr_ddn_beam_beam = zero_vector(nj)
  IF(.NOT. ASSOCIATED(fus_prod%neutr_ddn_beam_thermal%data))fus_prod%neutr_ddn_beam_thermal = zero_vector(nj)
  IF(.NOT. ASSOCIATED(fus_prod%neutr_ddn_knock%data))fus_prod%neutr_ddn_knock = zero_vector(nj)
  IF(.NOT. ASSOCIATED(fus_prod%neutr_ddn_tot%data))fus_prod%neutr_ddn_tot = zero_vector(nj)

  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 3760)
3760  FORMAT ('*  circum : plasma circumference, meters')
     WRITE  (niterdb, 10) dischg%circum 
  ELSE
     !Protected read for older statefiles
     dischg%circum = zeroc
     READ   (niterdb,  7,END = 2000) starflag
     READ   (niterdb, 10,END = 2000) dischg%circum
  END IF
! added 10/1/10
 
 
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 3770) 
3770  FORMAT ('* eta : resistivity, ohm m')
     WRITE  (niterdb, 10) (diffuse%eta%data(j), j=1,nj) 
  ELSE
     !Protected read for older statefiles
     READ   (niterdb,  7,END = 2000) starflag
     READ   (niterdb, 10,END = 2000) (work(j),j=1,nj)
     !diffuse%eta%data = load_Vector(work) ?? doesnt work
     CALL delete_vector(diffuse%eta)
     diffuse%eta = new_Vector(nj,work)
  END IF

  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 3780)
3780  FORMAT ('*  ftrap : trapped electron fraction')
     WRITE  (niterdb, 10) (diffuse%ftrap%data(j), j=1,nj) 
  ELSE
     !Protected read for older statefiles
     READ   (niterdb,  7,END = 2000) starflag
     READ   (niterdb, 10,END = 2000) (work(j),j=1,nj)
     !diffuse%ftrap = load_Vector(work)
     CALL delete_vector(diffuse%ftrap)
     diffuse%ftrap = new_Vector(nj,work)
  END IF

  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 3790)
3790  FORMAT ('*  xnues (eg nu*e): elect  coll freq/bounce freq')
     WRITE  (niterdb, 10) (diffuse%xnuse%data(j), j=1,nj) 
  ELSE
     !Protected read for older statefiles
     READ   (niterdb,  7,END = 2000) starflag
     READ   (niterdb, 10,END = 2000) (work(j),j=1,nj)
     !diffuse%xnuse           = load_Vector(work)
     CALL delete_vector(diffuse%xnuse)
     diffuse%xnuse = new_Vector(nj,work)
  END IF




  WRITE(label,FMT='(1pe16.8)')tGCNMs
  READ(label,FMT='(a)')st_asc_time
  st_asc_time = ADJUSTL(st_asc_time)
  label = '* xnus (eg nu*i): ion  coll freq/bounce freq ' &
       //st_asc_time(1:LEN_TRIM(st_asc_time))//' for species : '
  label = ADJUSTL(label)
  jp = 0
  ji = 0
  DO jj=1,nion
     IF (jj .LE. nprim)  jp = jp + 1
     IF (jj .GT. nprim)  ji = ji + 1
     IF (irwflag .EQ. 0) THEN
        IF(jj .LE. nprim)THEN
           tlabel = label(1:LEN_TRIM(label))//namep(jp)
        ELSE
           tlabel = label(1:LEN_TRIM(label))//namei(ji)
        ENDIF
        tlabel = ADJUSTL(tlabel)
        IF (iterdsc .NE. 0 )&
             WRITE  (niterdb, FMT='(A)') tlabel(1:LEN_TRIM(tlabel))
        work(1:nj) = diffuse%xnus(jj)%data(1:nj)
        WRITE  (niterdb, 10) (work(j), j=1,nj)
     ELSE
     !Protected read for older statefiles
        READ   (niterdb,  7,END = 2000) starflag
        READ   (niterdb, 10,END = 2000) (work(j), j=1,nj)
        !diffuse%xnus(jj) = load_Vector(work)
        CALL delete_vector(diffuse%xnus(jj))
        diffuse%xnus(jj) = new_Vector(nj,work)
     END IF
  END DO

! added 10/6/10

  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 3800)
3800  FORMAT ('* thermal-thermal ddn neutron production  rate')
     WRITE  (niterdb, 10) (fus_prod%neutr_ddn_th%data(j), j=1,nj) 
  ELSE
     !Protected read for older statefiles
     READ   (niterdb,  7,END = 2000) starflag
     READ   (niterdb, 10,END = 2000) (work(j),j=1,nj)
     !fus_prod%neutr_ddn_th   = load_Vector(work)
     CALL delete_vector(fus_prod%neutr_ddn_th  )
     fus_prod%neutr_ddn_th   = new_Vector(nj,work)
  END IF

 
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 3810)
3810  FORMAT ('* beam - beam  ddn neutron production  rate')
     WRITE  (niterdb, 10) (fus_prod%neutr_ddn_beam_beam%data(j), j=1,nj) 
  ELSE
     !Protected read for older statefiles
     READ   (niterdb,  7,END = 2000) starflag
     READ   (niterdb, 10,END = 2000) (work(j),j=1,nj)
     !fus_prod%neutr_ddn_beam_beam  = load_Vector(work)
     CALL delete_vector(fus_prod%neutr_ddn_beam_beam  )
     fus_prod%neutr_ddn_beam_beam   = new_Vector(nj,work)
  END IF


  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 3820)
3820  FORMAT ('* beam - thermal  ddn neutron production  rate')
     WRITE  (niterdb, 10) (fus_prod%neutr_ddn_beam_thermal%data(j), j=1,nj) 
  ELSE
     !Protected read for older statefiles
     READ   (niterdb,  7,END = 2000) starflag
     READ   (niterdb, 10,END = 2000) (work(j),j=1,nj)
     !fus_prod%neutr_ddn_beam_thermal  = load_Vector(work)
     CALL delete_vector(fus_prod%neutr_ddn_beam_thermal  )
     fus_prod%neutr_ddn_beam_thermal  = new_Vector(nj,work)
  END IF


  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 3830)
3830  FORMAT ('* knock on   ddn neutron production  rate')
     WRITE  (niterdb, 10) (fus_prod%neutr_ddn_knock%data(j), j=1,nj) 
  ELSE
     !Protected read for older statefiles
     READ   (niterdb,  7,END = 2000) starflag
     READ   (niterdb, 10,END = 2000) (work(j),j=1,nj)
     !fus_prod%neutr_ddn_knock  = load_Vector(work)
     CALL delete_vector(fus_prod%neutr_ddn_knock )
     fus_prod%neutr_ddn_knock   = new_Vector(nj,work)
  END IF


 IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 3840)
3840  FORMAT ('* ddn total  production  rate')
     WRITE  (niterdb, 10) (fus_prod%neutr_ddn_tot%data(j), j=1,nj) 
  ELSE
     !Protected read for older statefiles
     READ   (niterdb,  7,END = 2000) starflag
     READ   (niterdb, 10,END = 2000) (work(j),j=1,nj)
     !fus_prod%neutr_ddn_tot  = load_Vector(work)
     CALL delete_vector(fus_prod%neutr_ddn_tot )
     fus_prod%neutr_ddn_tot   = new_Vector(nj,work)
  END IF
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 610)
610  FORMAT ('*  total neutron rates,#/sec: thermal,beam-beam,beam-thermal, all')
     WRITE  (niterdb, 10)           fus_prod%total_neutr_ddn_th,            &
                                    fus_prod%total_neutr_ddn_beam_beam,    &
                                    fus_prod%total_neutr_ddn_beam_thermal, &
                                    fus_prod%total_neutr_ddn 
  ELSE
     READ   (niterdb,  7,END = 2000) starflag
     READ   (niterdb, 10,END= 2000) fus_prod%total_neutr_ddn_th,            &
                                    fus_prod%total_neutr_ddn_beam_beam,    &
                                    fus_prod%total_neutr_ddn_beam_thermal, &
                                    fus_prod%total_neutr_ddn
  ENDIF


  !-----------------------------------------------------------
  !new neutral beam items P_nfreya, added 3/18/11 HSJ
  ! modified 4/11/12
  !-----------------------------------------------------------


 IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 3850)
3850  FORMAT ('*  number of neutral beam injectors,grid size for beam')
     WRITE  (niterdb, 9)neut_beam%nbeams,neut_beam%nj_beam ! value of  zero is valid
  ELSE
     neut_beam%nj_beam = -1
     READ   (niterdb,  7,END = 2000) starflag
     READ   (niterdb, 9,END= 2000,ERR = 1999) neut_beam%nbeams,neut_beam%nj_beam
     go to 1998
1999 CONTINUE  ! this error path is not taken even if the field read is blank ???
               ! apparently ifort thinks blanck  is same as zero ????
     neut_beam%nj_beam = nj
1998 CONTINUE
     IF(neut_beam%nj_beam == 0)neut_beam%nj_beam = nj
  ENDIF



 !--------------------------------------------------------------------------------------------
 ! if we reading in the data then we need to make sure space is allocated for the following
 ! if we are writing out the data then these arrays allready exist,unless beams are
 ! not present. Handle this situation by setting neut_beam%nbeams =0
 ! Note that neut_beam%nbeams > 0 if either nubeam or P_Nfreya is called.
 ! But nubeam does not supply all the output in which case values are set to zero.
 !--------------------------------------------------------------------------------------------


   CALL neut_beam_allocate

nbgt0:  IF(neut_beam%nbeams .GT.0)THEN
            IF(SIZE(bptor) .LT. neut_beam%nbeams)THEN
              ! bptor is dimmed to  kb in nf_param
               WRITE(ncrt,FMT ='( "Error,inconistency in bptor dimension")')
              lerrno = iomaxerr + 358_I4B
              CALL  terminate(lerrno,nlog)
            ELSEIF(SIZE(bptor) .GT. neut_beam%nbeams)THEN
              bptor(neut_beam%nbeams+1:SIZE(bptor)) = 1.0e-03_DP
            ENDIF


    DO j=1,neut_beam%nbeams
       
       IF (irwflag .EQ. 0) THEN
          IF (iterdsc .NE. 0) &
               WRITE  (niterdb, 3920)
3920      FORMAT ('* Neutral beam power before aperture loss, watts:')
          WRITE  (niterdb, 10)(neut_beam%pbeam(i,j),i=1,ke_bm)
       ELSE
          READ   (niterdb,  7,END = 2000) starflag
          READ   (niterdb, 10,END = 2000) (neut_beam%pbeam(i,j),i=1,ke_bm)
       ENDIF

       IF (irwflag .EQ. 0) THEN
          IF (iterdsc .NE. 0) &
               WRITE  (niterdb, 3921)
3921      FORMAT ('* Neutral beam power into torus  watts or  thermalized power if nubeam used:')
          WRITE  (niterdb, 10)bptor(j)
       ELSE
          READ   (niterdb,  7,END = 2000) starflag
          READ   (niterdb, 10,END = 2000) bptor(j)
        ENDIF
    ENDDO


    DO j=1,neut_beam%nbeams
       IF (irwflag .EQ. 0) THEN
          IF (iterdsc .NE. 0) &
               WRITE  (niterdb, 3951)
3951      FORMAT ('* Neutral beam neutral intensity  bneut, #/sec:')
          WRITE  (niterdb, 10)(neut_beam%bneut(i,j),i=1,ke_bm)
       ELSE
          READ   (niterdb,  7,END = 2000) starflag
          READ   (niterdb, 10,END = 2000) (neut_beam%bneut(i,j),i=1,ke_bm)
       ENDIF
    ENDDO
    DO j=1,neut_beam%nbeams
       IF (irwflag .EQ. 0) THEN
          IF (iterdsc .NE. 0) &
               WRITE  (niterdb, 3950)
3950      FORMAT ('* Neutral beam ion intensity  bion, #/sec:')
          WRITE  (niterdb, 10)(neut_beam%bion(i,j),i=1,ke_bm)
       ELSE
          READ   (niterdb,  7,END = 2000) starflag
          READ   (niterdb, 10,END = 2000) (neut_beam%bion(i,j),i=1,ke_bm)
       ENDIF
    ENDDO
    DO j=1,neut_beam%nbeams
       IF (irwflag .EQ. 0) THEN
          IF (iterdsc .NE. 0) &
               WRITE  (niterdb, 3930)
3930      FORMAT ('* Neutral beam aperture loss:')
          WRITE  (niterdb, 10)(neut_beam%fap(i,j),i=1,ke_bm)
       ELSE
          READ   (niterdb,  7,END = 2000) starflag
          READ   (niterdb, 10,END = 2000) (neut_beam%fap(i,j),i=1,ke_bm)
       ENDIF
    ENDDO

    DO j=1,neut_beam%nbeams
       IF (irwflag .EQ. 0) THEN
          IF (iterdsc .NE. 0) &
               WRITE  (niterdb, 3932)
3932      FORMAT ('* Neutral beam wall loss:')
          WRITE  (niterdb, 10)(neut_beam%fwall(i,j),i=1,ke_bm)
       ELSE
          READ   (niterdb,  7,END = 2000) starflag
          READ   (niterdb, 10,END = 2000) (neut_beam%fwall(i,j),i=1,ke_bm)
       ENDIF
    ENDDO

    DO j=1,neut_beam%nbeams
       IF (irwflag .EQ. 0) THEN
          IF (iterdsc .NE. 0) &
               WRITE  (niterdb, 3934)
3934      FORMAT ('* Neutral beam orbit loss:')
          WRITE  (niterdb, 10)(neut_beam%forb(i,j),i=1,ke_bm)
       ELSE
          READ   (niterdb,  7,END = 2000) starflag
          READ   (niterdb, 10,END = 2000) (neut_beam%forb(i,j),i=1,ke_bm)
       ENDIF
    ENDDO




    DO j=1,neut_beam%nbeams
       IF (irwflag .EQ. 0) THEN
          IF (iterdsc .NE. 0) &
               WRITE  (niterdb, 3990)
3990      FORMAT ('* Neutral beam fraction of ions trapped for which error was detected :')
          WRITE  (niterdb, 10)(neut_beam%fber(i,j),i=1,ke_bm)
       ELSE
          READ   (niterdb,  7,END = 2000) starflag
          READ   (niterdb, 10,END = 2000) (neut_beam%fber(i,j),i=1,ke_bm)
       ENDIF
    ENDDO

    DO j=1,neut_beam%nbeams
       IF (irwflag .EQ. 0) THEN
          IF (iterdsc .NE. 0) &
               WRITE  (niterdb, 3991)
3991      FORMAT ('* Neutral beam fraction of ions trapped and not encircling :')
          WRITE  (niterdb, 10)(neut_beam%fb00(i,j),i=1,ke_bm)
       ELSE
          READ   (niterdb,  7,END = 2000) starflag
          READ   (niterdb, 10,END = 2000) (neut_beam%fb00(i,j),i=1,ke_bm)
       ENDIF
    ENDDO

    DO j=1,neut_beam%nbeams
       IF (irwflag .EQ. 0) THEN
          IF (iterdsc .NE. 0) &
               WRITE  (niterdb, 3992)
3992      FORMAT ('* Neutral beam fraction of ions trapped and axis-encircling:')
          WRITE  (niterdb, 10)(neut_beam%fb01(i,j),i=1,ke_bm)
       ELSE
          READ   (niterdb,  7,END = 2000) starflag
          READ   (niterdb, 10,END = 2000) (neut_beam%fb01(i,j),i=1,ke_bm)
       ENDIF
    ENDDO

    DO j=1,neut_beam%nbeams
       IF (irwflag .EQ. 0) THEN
          IF (iterdsc .NE. 0) &
               WRITE  (niterdb, 3993)
3993      FORMAT ('* Neutral beam fraction of ions passing and not encircling:')
          WRITE  (niterdb, 10)(neut_beam%fb10(i,j),i=1,ke_bm)
       ELSE
          READ   (niterdb,  7,END = 2000) starflag
          READ   (niterdb, 10,END = 2000) (neut_beam%fb10(i,j),i=1,ke_bm)
       ENDIF
    ENDDO

    DO j=1,neut_beam%nbeams
       IF (irwflag .EQ. 0) THEN
          IF (iterdsc .NE. 0) &
               WRITE  (niterdb, 3994)
3994      FORMAT ('* Neutral beam fraction of ions passing and axis-encircling:')
          WRITE  (niterdb, 10)(neut_beam%fb11(i,j),i=1,ke_bm)
       ELSE
          READ   (niterdb,  7,END = 2000) starflag
          READ   (niterdb, 10,END = 2000) (neut_beam%fb11(i,j),i=1,ke_bm)
       ENDIF
    ENDDO



    DO j=1,neut_beam%nbeams
       IF (irwflag .EQ. 0) THEN
          IF (iterdsc .NE. 0) &
               WRITE  (niterdb, 3995)
3995      FORMAT ('* Neutral beam orbit width of (m) of trapped and not encircling ions :')
          WRITE  (niterdb, 10)(neut_beam%wb00(i,j),i=1,ke_bm)
       ELSE
          READ   (niterdb,  7,END = 2000) starflag
          READ   (niterdb, 10,END = 2000) (neut_beam%wb00(i,j),i=1,ke_bm)
       ENDIF
    ENDDO

    DO j=1,neut_beam%nbeams
       IF (irwflag .EQ. 0) THEN
          IF (iterdsc .NE. 0) &
               WRITE  (niterdb, 3996)
3996      FORMAT ('* Neutral beam orbit width (m) of ions trapped and axis-encircling:')
          WRITE  (niterdb, 10)(neut_beam%wb01(i,j),i=1,ke_bm)
       ELSE
          READ   (niterdb,  7,END = 2000) starflag
          READ   (niterdb, 10,END = 2000) (neut_beam%wb01(i,j),i=1,ke_bm)
       ENDIF
    ENDDO

    DO j=1,neut_beam%nbeams
       IF (irwflag .EQ. 0) THEN
          IF (iterdsc .NE. 0) &
               WRITE  (niterdb, 3997)
3997      FORMAT ('* Neutral beam orbit width (m) of ions passing and not encircling:')
          WRITE  (niterdb, 10)(neut_beam%wb10(i,j),i=1,ke_bm)
       ELSE
          READ   (niterdb,  7,END = 2000) starflag
          READ   (niterdb, 10,END = 2000) (neut_beam%wb10(i,j),i=1,ke_bm)
       ENDIF
    ENDDO

    DO j=1,neut_beam%nbeams
       IF (irwflag .EQ. 0) THEN
          IF (iterdsc .NE. 0) &
               WRITE  (niterdb, 3998)
3998      FORMAT ('* Neutral beam orbit width of ions passing and axis-encircling:')
          WRITE  (niterdb, 10)(neut_beam%wb11(i,j),i=1,ke_bm)
       ELSE
          READ   (niterdb,  7,END = 2000) starflag
          READ   (niterdb, 10,END = 2000) (neut_beam%wb11(i,j),i=1,ke_bm)
       ENDIF
    ENDDO




    DO jb = 1,neut_beam%nbeams
       DO ie = 1,ke_bm
          IF (irwflag .EQ. 0) THEN
             IF (iterdsc .NE. 0) &
                  WRITE  (niterdb, 3860)jb,ie
3860         FORMAT ('* Neut beam prompt fast ion source #/(m^3 sec), injector,energy:',i5,x,i5)
             WRITE  (niterdb, 10)(neut_beam%sb(j,ie,jb), j=1,neut_beam%nj_beam) 
          ELSE
             READ   (niterdb,  7,END = 2000) starflag
             READ   (niterdb, 10,END= 2000) (neut_beam%sb(j,ie,jb), j=1,neut_beam%nj_beam)
          ENDIF
       ENDDO
    ENDDO


    DO jb = 1,neut_beam%nbeams
       DO ie = 1,ke_bm
          IF (irwflag .EQ. 0) THEN
             IF (iterdsc .NE. 0) &
                  WRITE  (niterdb, 3960)jb,ie
3960         FORMAT ('* Neut beam fast ion normalized birth rate (unitless), injector,energy:',i5,x,i5)
             WRITE  (niterdb, 10)(neut_beam%hibr(j,ie,jb), j=1,neut_beam%nj_beam) 
          ELSE
             READ   (niterdb,  7,END = 2000) starflag
             READ   (niterdb, 10,END= 2000) (neut_beam%hibr(j,ie,jb), j=1,neut_beam%nj_beam)
          ENDIF
       ENDDO
    ENDDO

    DO jb = 1,neut_beam%nbeams
       DO ie = 1,ke_bm
          IF (irwflag .EQ. 0) THEN
             IF (iterdsc .NE. 0) &
                  WRITE  (niterdb, 3972)jb,ie
3972         FORMAT ('* Neut beam normalized fast ion deposition rate (unitless), injector,energy:',i5,x,i5)
             WRITE  (niterdb, 10)(neut_beam%hdep(j,ie,jb), j=1,neut_beam%nj_beam) 
          ELSE
             READ   (niterdb,  7,END = 2000) starflag
             READ   (niterdb, 10,END= 2000) (neut_beam%hdep(j,ie,jb), j=1,neut_beam%nj_beam)
          ENDIF
       ENDDO
    ENDDO

    DO jb = 1,neut_beam%nbeams
       DO ie = 1,ke_bm
          IF (irwflag .EQ. 0) THEN
             IF (iterdsc .NE. 0) &
                  WRITE  (niterdb, 3970)jb,ie
3970         FORMAT ('* Neut beam pitch angle cosine , injector,energy:',i5,x,i5)
             WRITE  (niterdb, 10)(neut_beam%zeta(j,ie,jb), j=1,neut_beam%nj_beam) 
          ELSE
             READ   (niterdb,  7,END = 2000) starflag
             READ   (niterdb, 10,END= 2000) (neut_beam%zeta(j,ie,jb), j=1,neut_beam%nj_beam)
          ENDIF
       ENDDO
    ENDDO

    DO jb = 1,neut_beam%nbeams
       DO ie = 1,ke_bm
          IF (irwflag .EQ. 0) THEN
             IF (iterdsc .NE. 0) &
                  WRITE  (niterdb, 3870)jb,ie
3870         FORMAT ('* Neut beam prompt fast energy source Watts/m^3, injector,energy:',i5,x,i5)
             WRITE  (niterdb, 10)(neut_beam%qb(j,ie,jb), j=1,neut_beam%nj_beam) 
          ELSE
             READ   (niterdb,  7,END = 2000) starflag
             READ   (niterdb, 10,END= 2000) (neut_beam%qb(j,ie,jb), j=1,neut_beam%nj_beam)
          ENDIF
       ENDDO
    ENDDO

    DO jb = 1,neut_beam%nbeams
       DO ie = 1,ke_bm
          IF (irwflag .EQ. 0) THEN
             IF (iterdsc .NE. 0) &
                  WRITE  (niterdb, 3880)jb,ie
3880         FORMAT ('* Neut beam prompt parallel momentum rate density  Kg/(m2-s2), injector,energy:',i5,x,i5)
             WRITE  (niterdb, 10)(neut_beam%spb(j,ie,jb), j=1,neut_beam%nj_beam)
          ELSE
             READ   (niterdb,  7,END = 2000) starflag
             READ   (niterdb, 10,END= 2000) (neut_beam%spb(j,ie,jb), j=1,neut_beam%nj_beam)
          ENDIF
       ENDDO
    ENDDO

    DO jb = 1,neut_beam%nbeams
       DO ie = 1,ke_bm
          IF (irwflag .EQ. 0) THEN
             IF (iterdsc .NE. 0) &
                  WRITE  (niterdb, 3890)jb,ie
3890         FORMAT ('* Neut beam prompt angular momentum rate density  Kg/(m-s2), injector,energy:',i5,x,i5)
             WRITE  (niterdb, 10)(neut_beam%spbr(j,ie,jb), j=1,neut_beam%nj_beam) 
          ELSE
             READ   (niterdb,  7,END = 2000) starflag
             READ   (niterdb, 10,END= 2000) (neut_beam%spbr(j,ie,jb), j=1,neut_beam%nj_beam)
          ENDIF
       ENDDO
    ENDDO


    DO jb = 1,neut_beam%nbeams
       DO ie = 1,ke_bm
          IF (irwflag .EQ. 0) THEN
             IF (iterdsc .NE. 0) &
                  WRITE  (niterdb, 3900)jb,ie
3900         FORMAT ('* total ang momt in zone/deposition Kg m2/s, injector,energy:',i5,x,i5)
             WRITE  (niterdb, 10)(neut_beam%angmpf(j,ie,jb), j=1,neut_beam%nj_beam) 
          ELSE
             READ   (niterdb,  7,END = 2000) starflag
             READ   (niterdb, 10,END= 2000) (neut_beam%angmpf(j,ie,jb), j=1,neut_beam%nj_beam)
          ENDIF
       ENDDO
    ENDDO


    DO jb = 1,neut_beam%nbeams
       DO ie = 1,ke_bm
          IF (irwflag .EQ. 0) THEN
             IF (iterdsc .NE. 0) &
                  WRITE  (niterdb, 3910)jb,ie
3910         FORMAT ('* flux avg parallel beam mometum Kg m/s, injector,energy:',i5,x,i5)
             WRITE  (niterdb, 10)(neut_beam%pb0(j,ie,jb), j=1,neut_beam%nj_beam) 
          ELSE
             READ   (niterdb,  7,END = 2000) starflag
             READ   (niterdb, 10,END= 2000) (neut_beam%pb0(j,ie,jb), j=1,neut_beam%nj_beam)
          ENDIF
       ENDDO
    ENDDO
    !-----------------------------------------------------------------------
    ! Hot ion creation mode
    ! = 1, fraction of reactions producing electrons
    !                   = 2, fraction of reactions producing species 1 ion
    !                   = 3, fraction of reactions producing species 2 ion
    !-----------------------------------------------------------------------
    DO jb = 1,neut_beam%nbeams
       DO ie = 1,ke_bm
          DO k=1,kcm
             IF (irwflag .EQ. 0) THEN
                IF (iterdsc .NE. 0 .AND.  k == 1) &
                     WRITE  (niterdb, 3911)jb,ie,k
3911            FORMAT ('* hot ion creation mode: injector,energy,electrons:',i5,x,i5,x,i5)
                IF (iterdsc .NE. 0 .AND.  k == 2) &
                     WRITE  (niterdb, 3912)jb,ie,k
3912            FORMAT ('* hot ion creation mode: injector,energy,prim ion 1:',i5,x,i5,x,i5)
                IF (iterdsc .NE. 0 .AND.  k == 3) &
                     WRITE  (niterdb, 3913)jb,ie,k
3913            FORMAT ('* hot ion creation mode: injector,energy,prim ion 2:',i5,x,i5,x,i5)
                WRITE  (niterdb, 10)(neut_beam%hicm(j,ie,jb,k), j=1,neut_beam%nj_beam) 
             ELSE
                READ   (niterdb,  7,END = 2000) starflag
                READ   (niterdb, 10,END= 2000) (neut_beam%hicm(j,ie,jb,k), j=1,neut_beam%nj_beam)
             ENDIF
          ENDDO
       ENDDO
    ENDDO
 ENDIF  nbgt0 ! end writing neut_beam data if nbeams > 0




 

  !
  ! Electron paleoclassical thermal diffusivity
  !
  IF (irwflag .EQ. 0) THEN
     IF (iterdsc .NE. 0) WRITE (niterdb, 3147)
3147 FORMAT ('* electron paleoclassical thermal diffusivity,', &
            ' meter^2/second, on half grid')
     WRITE (niterdb, 10) (diffuse%chie_paleo%data(j)   , j=1,nj)
  ELSE
     READ  (niterdb,  7,END = 2042) starflag

     READ  (niterdb, 10,END = 2042) (work(j)   , j=1,nj)
     diffuse%chie_paleo  = new_Vector(nj,work)
  END IF
2042 continue
! add radiative loss output 8/31/2011 HSJ
     IF( .NOT. ASSOCIATED(pwrden%brems_nions))THEN
          ALLOCATE(pwrden%brems_nions(nion))
          DO jj =1,nion
             pwrden%brems_nions(jj)  = zero_vector(nj)
          ENDDO

     ENDIF
     IF(.NOT. ALLOCATED(brems_tot))THEN
         ALLOCATE(brems_tot(nion))
         brems_tot(:) = zeroc
     ENDIF

  label = '* bremssthralung and impurity radiation watts/(m^3 sec) for species : '
  label = ADJUSTL(label)
  jp = 0
  ji = 0
  DO jj=1,nion
     IF (jj .LE. nprim)  jp = jp + 1
     IF (jj .GT. nprim)  ji = ji + 1
     IF (irwflag .EQ. 0) THEN
        IF(jj .LE. nprim)THEN
           tlabel = label(1:LEN_TRIM(label))//namep(jp)
        ELSE
           tlabel = label(1:LEN_TRIM(label))//namei(ji)
        ENDIF
        tlabel = ADJUSTL(tlabel)
        IF (iterdsc .NE. 0 )&
             WRITE  (niterdb, FMT='(A)') tlabel(1:LEN_TRIM(tlabel))
        IF(jj == 1) &
        WRITE  (niterdb, 10)qrad_tot,pfuse_tot
        work(1:nj) = get_values(pwrden%brems_nions(jj))
        WRITE  (niterdb, 10) (work(j), j=1,nj)
        WRITE  (niterdb, 10)brems_tot(jj)
     ELSE
     !Protected read for older statefiles
        READ   (niterdb,  7,END = 2000) starflag
        IF(jj == 1) &
        READ   (niterdb, 10,END = 2043) qrad_tot,pfuse_tot
        READ   (niterdb, 10,END = 2043) (work(j), j=1,nj)
        READ   (niterdb, 10,END = 2043) brems_tot(jj)
        CALL delete_vector(pwrden%brems_nions(jj))
        pwrden%brems_nions(jj) = new_Vector(nj,work)
     END IF
  END DO
2043 CONTINUE
  !----------------------------------------------------------------------
  ! --- glf turbulent particle,energy and momentum   flux
  ! --- It is assumed that glf always has 3 species;
  ! --- 1) electron
  ! --- 2) effective ion
  ! ----3) effective impurity
  !---------------------------------------------------------------------
     IF (irwflag .EQ. 0) THEN
        IF (iterdsc .NE.0 )THEN
             WRITE  (niterdb, 1157)
1157         FORMAT ('*  glf electron particle flux,', &
             ' #/(meter2 sec)' )
             WRITE  (niterdb, 10) (glf_p_output(j,1), j=1,nj)
        ENDIF
        IF(iterdsc .NE. 0 )THEN
             WRITE  (niterdb, 1158)
1158         FORMAT ('*  glf (effective) primary ion particle  flux,', &
             ' #/(meter2 sec)' )
             WRITE  (niterdb, 10) (glf_p_output(j,2), j=1,nj)
        ENDIF
        IF(iterdsc .NE. 0 )THEN
             WRITE  (niterdb, 1159)
1159         FORMAT ('*  glf (effective) impurity  ion particle  flux,', &
             ' #/(meter2 sec)' )
             WRITE  (niterdb, 10) (glf_p_output(j,3), j=1,nj)
        ENDIF

        IF (iterdsc .NE. 0)THEN
             WRITE  (niterdb, 1161)
1161         FORMAT ('* glf  electron energy flux', &
                     ' Joules/(meter**2sec)')
             WRITE  (niterdb, 10) (glf_e_output(j,1), j=1,nj)
        ENDIF
        IF (iterdsc .NE. 0)THEN
             WRITE  (niterdb, 1164)
1164         FORMAT ('* glf  etg contribution to energy flux', &
                     ' Joules/(meter**2sec)')
             WRITE  (niterdb, 10) (glf_etg_output(j), j=1,nj)
        ENDIF

        IF (iterdsc .NE. 0)THEN
             WRITE  (niterdb, 1162)
1162         FORMAT ('* glf  effective primary ion  energy flux', &
                     ' Joules/(meter**2sec)')
             WRITE  (niterdb, 10) (glf_e_output(j,2), j=1,nj)
        ENDIF
        IF (iterdsc .NE. 0)THEN
             WRITE  (niterdb, 1163)
1163         FORMAT ('* glf effective impurity ion  energy flux', &
                     ' Joules/(meter**2sec)')
             WRITE  (niterdb, 10) (glf_e_output(j,3), j=1,nj)
        ENDIF

        IF (iterdsc .NE. 0)THEN
             WRITE  (niterdb, 1261)
1261         FORMAT ('* glf  elct momentum flux', &
                     ' kg/sec**2')
             WRITE  (niterdb, 10) (glf_m_output(j,1), j=1,nj)
        ENDIF
        IF (iterdsc .NE. 0)THEN
             WRITE  (niterdb, 1262)
1262         FORMAT ('* glf  effective primary ion momentum flux', &
                     ' kg/(sec**2)')
             WRITE  (niterdb, 10) (glf_m_output(j,2), j=1,nj)
        ENDIF
        IF (iterdsc .NE. 0)THEN
             WRITE  (niterdb, 1263)
1263         FORMAT ('* glf effective impurity ion momentum flux', &
                     ' kg/(sec**2)')
             WRITE  (niterdb, 10) (glf_m_output(j,3), j=1,nj)
        ENDIF

        IF (iterdsc .NE. 0)THEN
             WRITE  (niterdb, 1265)
1265         FORMAT ('* glf highest frequency', &
                     ' 1/sec')
             WRITE  (niterdb, 10) (glf_anfreq_output(j), j=1,nj)
        ENDIF
        IF (iterdsc .NE. 0)THEN
             WRITE  (niterdb, 1266)
1266         FORMAT ('* glf second highest frequency', &
                     ' 1/(sec')
             WRITE  (niterdb, 10) (glf_anfreq2_output(j), j=1,nj)
        ENDIF
       IF (iterdsc .NE. 0)THEN
             WRITE  (niterdb, 1267)
1267         FORMAT ('* glf growth rate leading mode', &
                     ' 1/sec')
             WRITE  (niterdb, 10) (glf_anrate_output(j), j=1,nj)
        ENDIF
        IF (iterdsc .NE. 0)THEN
             WRITE  (niterdb, 1268)
1268         FORMAT ('* glf growth rate second mode', &
                     ' 1/sec')
             WRITE  (niterdb, 10) (glf_anrate2_output(j), j=1,nj)
        ENDIF
        IF (iterdsc .NE. 0)THEN
             WRITE  (niterdb, 1269)
1269         FORMAT ('* glf net growth rate ions ', &
                     '1/sec')
             WRITE  (niterdb, 10) (glf_gamma_net_i_output(j), j=1,nj)
        ENDIF
        IF (iterdsc .NE. 0)THEN
             WRITE  (niterdb, 1270)
1270         FORMAT ('* glf net growth rate electrons ', &
                     ' 1/sec')
             WRITE  (niterdb, 10) (glf_gamma_net_e_output(j), j=1,nj)
        ENDIF
     ELSE
        IF(ALLOCATED(glf_p_output))DEALLOCATE(glf_p_output)
        ALLOCATE(glf_p_output(nj,3))
        glf_p_output  = zeroc
        READ   (niterdb,  7,END = 2000) line
        READ   (niterdb, 10,END = 2000) (glf_p_output(j,1)       , j=1,nj)
        READ   (niterdb,  7,END = 2000) line
        READ   (niterdb, 10,END = 2000) (glf_p_output(j,2)       , j=1,nj)
        READ   (niterdb,  7,END = 2000) line
        READ   (niterdb, 10,END = 2000) (glf_p_output(j,3)       , j=1,nj)


        IF(ALLOCATED(glf_e_output))DEALLOCATE(glf_e_output)
        ALLOCATE(glf_e_output(nj,3))
        glf_e_output = zeroc
        READ   (niterdb,  7,END = 2000) line
        READ   (niterdb, 10,END = 2000) (glf_e_output(j,1)       , j=1,nj)

        IF(ALLOCATED(glf_etg_output))DEALLOCATE(glf_etg_output)
        ALLOCATE(glf_etg_output(nj))
        glf_etg_output = zeroc
        READ   (niterdb,  7,END = 2000) line
        READ   (niterdb, 10,END = 2000) (glf_etg_output(j)       , j=1,nj)

        READ   (niterdb,  7,END = 2000) line
        READ   (niterdb, 10,END = 2000) (glf_e_output(j,2)       , j=1,nj)
        READ   (niterdb,  7,END = 2000) line
        READ   (niterdb, 10,END = 2000) (glf_e_output(j,3)       , j=1,nj)


        IF(ALLOCATED(glf_m_output))DEALLOCATE(glf_m_output)
        ALLOCATE(glf_m_output(nj,3))
        glf_m_output = zeroc
        READ   (niterdb,  7,END = 2000) line
        READ   (niterdb, 10,END = 2000) (glf_m_output(j,1)       , j=1,nj)
        READ   (niterdb,  7,END = 2000) line
        READ   (niterdb, 10,END = 2000) (glf_m_output(j,2)       , j=1,nj)
        READ   (niterdb,  7,END = 2000) line
        READ   (niterdb, 10,END = 2000) (glf_m_output(j,3)       , j=1,nj)

        IF(ALLOCATED(glf_anfreq_output))DEALLOCATE(glf_anfreq_output)
        ALLOCATE(glf_anfreq_output(nj))
        glf_anfreq_output(:) = zeroc
        READ   (niterdb,  7,END = 2000) line
        READ   (niterdb, 10,END = 2000) (glf_anfreq_output(j)       , j=1,nj)


        IF(ALLOCATED(glf_anfreq2_output))DEALLOCATE(glf_anfreq2_output)
        ALLOCATE(glf_anfreq2_output(nj))
        glf_anfreq2_output(:) = zeroc
        READ   (niterdb,  7,END = 2000) line
        READ   (niterdb, 10,END = 2000) (glf_anfreq2_output(j)      , j=1,nj)

        IF(ALLOCATED(glf_anrate_output))DEALLOCATE(glf_anrate_output)
        ALLOCATE(glf_anrate_output(nj))
        glf_anrate_output(:) = zeroc
        READ   (niterdb,  7,END = 2000) line
        READ   (niterdb, 10,END = 2000) (glf_anrate_output(j)       , j=1,nj)


        IF(ALLOCATED(glf_anrate2_output))DEALLOCATE(glf_anrate2_output)
        ALLOCATE(glf_anrate2_output(nj))
        glf_anrate2_output(:) = zeroc
        READ   (niterdb,  7,END = 2000) line
        READ   (niterdb, 10,END = 2000) (glf_anrate2_output(j)      , j=1,nj)

        IF(ALLOCATED(glf_gamma_net_i_output))DEALLOCATE(glf_gamma_net_i_output)
        ALLOCATE(glf_gamma_net_i_output(nj))
        glf_gamma_net_i_output(:) = zeroc
        READ   (niterdb,  7,END = 2000) line
        READ   (niterdb, 10,END = 2000) (glf_gamma_net_i_output(j)      , j=1,nj)


        IF(ALLOCATED(glf_gamma_net_e_output))DEALLOCATE(glf_gamma_net_e_output)
        ALLOCATE(glf_gamma_net_e_output(nj))
        glf_gamma_net_e_output  = zeroc
        READ   (niterdb,  7,END = 2000) line
        READ   (niterdb, 10,END = 2000) (glf_gamma_net_e_output(j)      , j=1,nj)

     END IF


2000 CONTINUE

  IF (irwflag .EQ. 0) THEN
     IF(mhd_dat%betan%size .LE. 0)THEN
        mhd_dat%betan = zero_vector(nj)
     ENDIF
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 6000)
6000 FORMAT ('* local betan (pressure divided by btgeom**2/2u0/I/aB ')
     WRITE  (niterdb, 10) (mhd_dat%betan%data(j), j=1,nj) 
  ELSE
     mhd_dat%betan%size = nj              ! if not present in state file set zero
     mhd_dat%betan = zero_Vector(nj) ! needed because we dstribute the values

     READ   (niterdb,  7,END = 2041) starflag
     READ   (niterdb, 10,END = 2041) (work(j), j=1,nj)
     mhd_dat%betan = new_Vector(nj,work)
  END IF
2041  CONTINUE


  IF (irwflag .EQ. 0) THEN
     IF(.NOT. ASSOCIATED(mhd_dat%curpar%data))THEN
        mhd_dat%curpar = zero_vector(nj)
     ENDIF
     IF (iterdsc .NE. 0) &
     WRITE  (niterdb, 6005)
6005 FORMAT ('* parallel current density <J.B/Bt0>')
     WRITE  (niterdb, 10) (mhd_dat%curpar%data(j), j=1,nj) 
  ELSE
     mhd_dat%curpar%size = nj              ! if not present in state file set zero
     mhd_dat%curpar = zero_Vector(nj)      ! needed because we dstribute the values
     READ   (niterdb,  7,END = 2001) starflag
     READ   (niterdb, 10,END = 2001) (work(j), j=1,nj)
     mhd_dat%curpar = new_Vector(nj,work)
  END IF
  
     IF(ALLOCATED(worknpsi))THEN
       DEALLOCATE(worknpsi)
     ENDIF
     ALLOCATE(worknpsi(mhd_dat%npsi)) ! assumes all usage of work below is of size nj

2001 CONTINUE
  IF (irwflag .EQ. 0) THEN
     IF(.NOT. ASSOCIATED(mhd_dat%qpsinpsi%data))THEN
        mhd_dat%qpsinpsi = zero_vector(mhd_dat%npsi)
     ENDIF
     IF (iterdsc .NE. 0) &
     WRITE  (niterdb, 6006)
6006 FORMAT ('* q value on mhd grid')
     WRITE  (niterdb, 10) (mhd_dat%qpsinpsi%data(j), j=1,mhd_dat%npsi) 
  ELSE
     !mhd_dat%qpsinpsi%size = mhd_dat%npsi              ! if not present in state file set zero
     mhd_dat%qpsinpsi = zero_Vector(mhd_dat%npsi)      ! needed because we dstribute the values
     mhd_dat%mhd_info_avail = .FALSE.
     READ   (niterdb,  7,END = 2002) starflag
     READ   (niterdb, 10,END = 2002) (worknpsi(j), j=1,mhd_dat%npsi)
     mhd_dat%qpsinpsi = new_Vector(mhd_dat%npsi,worknpsi)
     mhd_dat%mhd_info_avail = .TRUE. ! used by forcebal
  END IF

2002 CONTINUE
  IF (irwflag .EQ. 0) THEN
     IF(.NOT. ASSOCIATED(mhd_dat%ffprimnpsi%data))THEN
        mhd_dat%ffprimnpsi = zero_vector(mhd_dat%npsi)
     ENDIF
     IF (iterdsc .NE. 0) &
     WRITE  (niterdb, 6007)
6007 FORMAT ('* ffprim  on mhd grid,kg /( A s^2)')
     WRITE  (niterdb, 10) (mhd_dat%ffprimnpsi%data(j), j=1,mhd_dat%npsi) 
  ELSE
                                                     ! if not present in state file set zero
     mhd_dat%ffprimnpsi = zero_Vector(mhd_dat%npsi)      ! needed because we dstribute the values
     READ   (niterdb,  7,END = 2003) starflag
     READ   (niterdb, 10,END = 2003) (worknpsi(j), j=1,mhd_dat%npsi)
     mhd_dat%ffprimnpsi = new_Vector(mhd_dat%npsi,worknpsi)
  END IF


2003 CONTINUE
  IF (irwflag .EQ. 0) THEN
     IF(.NOT. ASSOCIATED(mhd_dat%pressnpsi%data))THEN
        mhd_dat%pressnpsi = zero_vector(mhd_dat%npsi)
     ENDIF
     IF (iterdsc .NE. 0) &
     WRITE  (niterdb, 6008)
6008 FORMAT ('* pressure   on mhd grid nt/m^2')
     WRITE  (niterdb, 10) (mhd_dat%pressnpsi%data(j), j=1,mhd_dat%npsi) 
  ELSE
     mhd_dat%pressnpsi  = zero_Vector(mhd_dat%npsi)
     READ   (niterdb,  7,END = 2004) starflag
     READ   (niterdb, 10,END = 2004) (worknpsi(j), j=1,mhd_dat%npsi)
     mhd_dat%pressnpsi = new_Vector(mhd_dat%npsi,worknpsi)
  END IF

2004 CONTINUE
  IF (irwflag .EQ. 0) THEN
     IF(.NOT. ASSOCIATED(mhd_dat%pprimnpsi%data))THEN
        mhd_dat%pprimnpsi = zero_vector(mhd_dat%npsi)
     ENDIF
     IF (iterdsc .NE. 0) &
     WRITE  (niterdb, 6009)
6009 FORMAT ('* p-prime (eg dpress/dpsi)   on mhd grid A/m^3')
     WRITE  (niterdb, 10) (mhd_dat%pprimnpsi%data(j), j=1,mhd_dat%npsi) 
  ELSE
     mhd_dat%pprimnpsi = zero_Vector(mhd_dat%npsi)
     READ   (niterdb,  7,END = 2005) starflag
     READ   (niterdb, 10,END = 2005) (worknpsi(j), j=1,mhd_dat%npsi)
     mhd_dat%pprimnpsi = new_Vector(mhd_dat%npsi,worknpsi)
  END IF

 


2005 CONTINUE
     
  IF (irwflag .EQ. 0) THEN
     IF(.NOT. ASSOCIATED(pwrden%omegale%data))THEN
        pwrden%omegale = zero_vector(nj)
     ENDIF
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 6010)
6010 FORMAT ('* omegale, beam electron energy correction due to rotation,w/m^3 ')
     WRITE  (niterdb, 10) (pwrden%omegale%data(j), j=1,nj) 
  ELSE
     pwrden%omegale%size = nj              ! if not present in state file set zero
     pwrden%omegale = zero_Vector(nj)      ! needed because we dstribute the values

     READ   (niterdb,  7,END = 2006) starflag
     READ   (niterdb, 10,END = 2006) (work(j), j=1,nj)
     pwrden%omegale%data(1:nj)  = work(1:nj)
  END IF


2006 CONTINUE


  IF (irwflag .EQ. 0) THEN
     IF(.NOT. ASSOCIATED(pwrden%qomegapi%data))THEN
        pwrden%qomegapi = zero_vector(nj)
     ENDIF
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 6011)
6011 FORMAT ('* qomegapi , beam ion  energy correction due to rotation ,w/m^3 ')
     WRITE  (niterdb, 10) (pwrden%qomegapi%data(j), j=1,nj) 
  ELSE
     pwrden%qomegapi%size = nj              ! if not present in state file set zero
     pwrden%qomegapi = zero_Vector(nj)      ! needed because we dstribute the values

     READ   (niterdb,  7,END = 2007) starflag
     READ   (niterdb, 10,END = 2007) (work(j), j=1,nj)
     pwrden%qomegapi%data(1:nj)  = work(1:nj)
  END IF


2007 CONTINUE


  IF (irwflag .EQ. 0) THEN
     IF(.NOT. ASSOCIATED(pwrden%qangce%data))THEN
        pwrden%qangce = zero_vector(nj)
     ENDIF
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 6012)
6012 FORMAT ('* qangce , beam ion  energy correction due to rotation ,w/m^3 ')
     WRITE  (niterdb, 10) (pwrden%qangce%data(j), j=1,nj) 
  ELSE
     pwrden%qangce%size = nj              ! if not present in state file set zero
     pwrden%qangce = zero_Vector(nj)      ! needed because we dstribute the values

     READ   (niterdb,  7,END = 2008) starflag
     READ   (niterdb, 10,END = 2008) (work(j), j=1,nj)
     pwrden%qangce%data(1:nj)  = work(1:nj)
  END IF


2008 CONTINUE


  IF (irwflag .EQ. 0) THEN
     IF(.NOT. ASSOCIATED(pwrden%sprcxre%data))THEN
        pwrden%sprcxre = zero_vector(nj)
     ENDIF
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 6013)
6013 FORMAT ('* sprcxre , beam ion  energy correction due to rotation ,w/m^3 ')
     WRITE  (niterdb, 10) (pwrden%sprcxre%data(j), j=1,nj) 
  ELSE
     pwrden%sprcxre%size = nj              ! if not present in state file set zero
     pwrden%sprcxre = zero_Vector(nj)      ! needed because we dstribute the values

     READ   (niterdb,  7,END = 2009) starflag
     READ   (niterdb, 10,END = 2009) (work(j), j=1,nj)
     pwrden%sprcxre%data(:)  = work(1:nj)
  END IF


2009 CONTINUE

  IF (irwflag .EQ. 0) THEN
     IF(.NOT. ASSOCIATED(pwrden%spreimpe%data))THEN
        pwrden%spreimpe = zero_vector(nj)
     ENDIF
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 6014)
6014 FORMAT ('* spreimpe , beam ion  energy correction due to rotation ,w/m^3 ')
     WRITE  (niterdb, 10) (pwrden%spreimpe%data(j), j=1,nj) 
  ELSE
     pwrden%spreimpe%size = nj              ! if not present in state file set zero
     pwrden%spreimpe = zero_Vector(nj)      ! needed because we dstribute the values

     READ   (niterdb,  7,END = 2010) starflag
     READ   (niterdb, 10,END = 2010) (work(j), j=1,nj)
     pwrden%spreimpe%data(1:nj)  = work(1:nj)
  END IF


2010 CONTINUE

  IF (irwflag .EQ. 0) THEN
     IF(.NOT. ASSOCIATED(pwrden%sprcxree%data))THEN
        pwrden%sprcxree = zero_vector(nj)
     ENDIF
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 6015)
6015 FORMAT ('* sprcxree , beam ion  energy correction due to rotation ,w/m^3 ')
     WRITE  (niterdb, 10) (pwrden%sprcxree%data(j), j=1,nj) 
  ELSE
     pwrden%sprcxree%size = nj              ! if not present in state file set zero
     pwrden%sprcxree = zero_Vector(nj)      ! needed because we dstribute the values
                                            ! (guard against end statment in read)

     READ   (niterdb,  7,END = 2011) starflag
     READ   (niterdb, 10,END = 2011) (work(j), j=1,nj)
     pwrden%sprcxree%data(:)  = work(1:nj)
  END IF


2011 CONTINUE


  If(neut_beam%nbeams .GT. 0)THEN
     IF(.NOT. ASSOCIATED(neut_beam%fbcur))THEN
        ALLOCATE(neut_beam%fbcur(ke_bm,neut_beam%nbeams))
        neut_beam%fbcur(:,:) = zeroc
     ENDIF
     IF (irwflag .EQ. 0) THEN
        IF (iterdsc .NE. 0) &
        WRITE  (niterdb, 6016)
6016    FORMAT ('* full,half,third neutral beam current fractions ')
        DO j=1,neut_beam%nbeams
           WRITE  (niterdb, 10) (neut_beam%fbcur(ie,j),ie=1,ke_bm) 
        ENDDO
     ELSE      
        READ   (niterdb,  7,END = 2012) starflag
        DO j=1,neut_beam%nbeams
           READ(niterdb, 10,END = 2012) (neut_beam%fbcur(ie,j),ie=1,ke_bm) 
        ENDDO
     END IF



    IF(.NOT. ASSOCIATED(neut_beam%prompt_pwr_in_plasma))THEN
        ALLOCATE(neut_beam%prompt_pwr_in_plasma(ke_bm,neut_beam%nbeams))
        neut_beam%prompt_pwr_in_plasma(:,:) = zeroc
    ENDIF
     IF (irwflag .EQ. 0) THEN
        IF (iterdsc .NE. 0) &
        WRITE  (niterdb, 6017)
6017    FORMAT ('* prompt power in plasma,by energy and beamlet W/m^3 ')
        DO j=1,neut_beam%nbeams
           WRITE  (niterdb, 10) (neut_beam%prompt_pwr_in_plasma(ie,j),ie=1,ke_bm) 
        ENDDO
     ELSE      
        READ   (niterdb,  7,END = 2012) starflag
        DO j=1,neut_beam%nbeams
           READ(niterdb, 10,END = 2012) (neut_beam%prompt_pwr_in_plasma(ie,j),ie=1,ke_bm) 
        ENDDO
     END IF
     IF(.NOT. ASSOCIATED(neut_beam%ebeam))THEN
        ALLOCATE(neut_beam%ebeam(ke_bm,neut_beam%nbeams))
        neut_beam%ebeam(:,:) = zeroc
     ENDIF
     IF (irwflag .EQ. 0) THEN
        IF (iterdsc .NE. 0) &
        WRITE  (niterdb, 6018)
6018    FORMAT ('* beam energy,full half third,kev ')
        DO j=1,neut_beam%nbeams
           WRITE  (niterdb, 10) (neut_beam%ebeam(ie,j),ie=1,ke_bm) 
        ENDDO
     ELSE      
        READ   (niterdb,  7,END = 2012) starflag
        DO j=1,neut_beam%nbeams
           READ(niterdb, 10,END = 2012) (neut_beam%ebeam(ie,j),ie=1,ke_bm) 
        ENDDO
     END IF

  ENDIF


2012 CONTINUE



  IF (irwflag .EQ. 0) THEN
     IF(.NOT. ASSOCIATED(profile%vpol%data))THEN
        profile%vpol = zero_vector(nj)
     ENDIF
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 6019)
6019 FORMAT ('* poloidal rotation velocity, m/sec ')
     WRITE  (niterdb, 10) (profile%vpol%data(j), j=1,nj) 
  ELSE
     profile%vpol%size = nj                 ! if not present in state file set zero
     profile%vpol = zero_Vector(nj)         ! needed because we distribute the values
                                            ! (may encounter end statement in read)

     READ   (niterdb,  7,END = 2013) starflag
     READ   (niterdb, 10,END = 2013) (work(j), j=1,nj)
     profile%vpol%data(1:nj)  = work(1:nj)
  END IF


2013 CONTINUE




  IF (irwflag .EQ. 0) THEN
     IF(.NOT. ASSOCIATED(profile%vpol_nclass%data))THEN
        profile%vpol_nclass = zero_vector(nj)
     ENDIF
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 6020)
6020 FORMAT ('* poloidal rotation velocity,forcebal/Nclass model, m/sec ')
     WRITE  (niterdb, 10) (profile%vpol_nclass%data(j), j=1,nj) 
  ELSE
     profile%vpol_nclass%size = nj                 ! if not present in state file set zero
     profile%vpol_nclass = zero_Vector(nj)         ! needed because we distribute the values
                                                   ! (may encounter end statement in read)
     READ   (niterdb,  7,END = 2014) starflag
     READ   (niterdb, 10,END = 2014) (work(j), j=1,nj)
     profile%vpol_nclass%data(1:nj)  = work(1:nj)
  END IF


2014 CONTINUE


  IF (irwflag .EQ. 0) THEN
     IF(.NOT. ASSOCIATED(profile%vpar%data))THEN
        profile%vpar = zero_vector(nj)
     ENDIF
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 6021)
6021 FORMAT ('* parallel  rotation velocity, m/sec ')
     WRITE  (niterdb, 10) (profile%vpar%data(j), j=1,nj) 
  ELSE
     profile%vpar%size = nj                 ! if not present in state file set zero
     profile%vpar = zero_Vector(nj)         ! needed because we distribute the values
                                            ! (may encounter end statement in read)

     READ   (niterdb,  7,END = 2015) starflag
     READ   (niterdb, 10,END = 2015) (work(j), j=1,nj)
     profile%vpar%data(1:nj)  = work(1:nj)
  END IF


2015 CONTINUE




  IF (irwflag .EQ. 0) THEN
     IF(.NOT. ASSOCIATED(profile%vpar_nclass%data))THEN
        profile%vpar_nclass = zero_vector(nj)
     ENDIF
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, 6022)
6022 FORMAT ('* parallel rotation velocity,forcebal/Nclass model, m/sec ')
     WRITE  (niterdb, 10) (profile%vpar_nclass%data(j), j=1,nj) 
  ELSE
     profile%vpar_nclass%size = nj                 ! if not present in state file set zero
     profile%vpar_nclass = zero_Vector(nj)         ! needed because we distribute the values
                                                   ! (may encounter end statement in read)
     READ   (niterdb,  7,END = 2016) starflag
     READ   (niterdb, 10,END = 2016) (work(j), j=1,nj)
     profile%vpar_nclass%data(1:nj)  = work(1:nj)
  END IF


2016 CONTINUE


! add nclass vpinch  output 5/3/2012 HSJ
     IF( .NOT. ASSOCIATED(diffuse%vpinch_nclass))THEN
          ALLOCATE(diffuse%vpinch_nclass(nion+1))         ! + 1 for electron species
          DO jj =1,nion+1
             diffuse%vpinch_nclass(jj)  = zero_vector(nj)
          ENDDO
     ENDIF

  label = '* Nclass derived total pinch velocity,m/sec, for species : '
  label = ADJUSTL(label)
  jp = 0
  ji = 0
  DO jj=1,nion+1
     IF (jj .LE. nprim)  jp = jp + 1
     IF (jj .GT. nprim)  ji = ji + 1

     IF (irwflag .EQ. 0) THEN
        IF(jj .LE. nprim)THEN
           tlabel = label(1:LEN_TRIM(label))//namep(jp)
        ELSEIf(jj .le. nion)THEN
           tlabel = label(1:LEN_TRIM(label))//namei(ji)
        ELSE
           tlabel = label(1:LEN_TRIM(label))//'electrons'
        ENDIF

        tlabel = ADJUSTL(tlabel)
        IF (iterdsc .NE. 0 )&
             WRITE  (niterdb, FMT='(A)') tlabel(1:LEN_TRIM(tlabel))
        work(1:nj) = get_values(diffuse%vpinch_nclass(jj))
        WRITE  (niterdb, 10) (work(j), j=1,nj)
     ELSE
     !Protected read for older statefiles
        READ   (niterdb,  7,END = 2017) starflag
        READ   (niterdb, 10,END = 2017) (work(j), j=1,nj)
        CALL delete_vector(diffuse%vpinch_nclass(jj))
        diffuse%vpinch_nclass(jj) = new_Vector(nj,work)
     END IF
  END DO


2017 CONTINUE

  !---------------------------------------------------------------------
  ! multimode related output:
  !---------------------------------------------------------------------
  IF (irwflag .EQ. 0) THEN

     IF(.NOT. ASSOCIATED(diffuse%mmm_gammaDBM%data) &
          .OR. diffuse%mmm_gammaDBM%size .NE. nj  )THEN
        ! when running glf or tglf this array may not be the right size
        ! if grid size changes are involved. 
        diffuse%mmm_gammaDBM = zero_vector(nj)
     ENDIF
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb,FMT='("* multimode - growth rate most unstable DRIB mode [1/sec]:")')
     WRITE  (niterdb, 10) (diffuse%mmm_gammaDBM%data(j), j=1,nj) 
  ELSE
     diffuse%mmm_gammaDBM%size = nj                 ! if not present in state file set to zero
     diffuse%mmm_gammaDBM = zero_Vector(nj)         ! needed because we distribute the values
                                                    ! (may encounter end statement in read)
     READ   (niterdb,  7,END = 2018) starflag
     READ   (niterdb, 10,END = 2018) (work(j), j=1,nj)
     diffuse%mmm_gammaDBM%data(1:nj)  = work(1:nj)
  END IF

2018 CONTINUE

  IF (irwflag .EQ. 0) THEN
     IF(.NOT. ASSOCIATED(diffuse%mmm_omegaDBM%data) &
         .OR. diffuse%mmm_omegaDBM%size .NE. nj  )THEN
        ! when running glf or tglf this array may not be the right size
        ! if grid size changes are involved. 
        diffuse%mmm_omegaDBM = zero_vector(nj)
     ENDIF
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb,FMT='("* multimode - freq,rad/sec, most unstable DRIB mode")')

     WRITE  (niterdb, 10) (diffuse%mmm_omegaDBM%data(j), j=1,nj) 
  ELSE
     diffuse%mmm_omegaDBM%size = nj                 ! if not present in state file set to zero
     diffuse%mmm_omegaDBM = zero_Vector(nj)         ! needed because we distribute the values
                                                   ! (may encounter end statement in read)
     READ   (niterdb,  7,END = 2019) starflag
     READ   (niterdb, 10,END = 2019) (work(j), j=1,nj)
     diffuse%mmm_omegaDBM%data(1:nj)  = work(1:nj)
  END IF
2019 CONTINUE

  IF (irwflag .EQ. 0) THEN
     IF(.NOT. ASSOCIATED(diffuse%mmm_xdi%data) &
          .OR. diffuse%mmm_xdi%size .ne. nj)THEN
        diffuse%mmm_xdi = zero_vector(nj)
     ENDIF
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, FMT='("* multimode - eff ion hydrogenic  diffusivity,M^2/sec")')
     WRITE  (niterdb, 10) (diffuse%mmm_xdi%data(j), j=1,nj) 
  ELSE
     diffuse%mmm_xdi%size = nj                 ! if not present in state file set to zero
     diffuse%mmm_xdi = zero_Vector(nj)         ! needed because we distribute the values
                                                   ! (may encounter end statement in read)
     READ   (niterdb,  7,END = 2020) starflag
     READ   (niterdb, 10,END = 2020) (work(j), j=1,nj)
     diffuse%mmm_xdi%data(1:nj)  = work(1:nj)
  END IF
2020 CONTINUE


  IF (irwflag .EQ. 0) THEN
     IF(.NOT. ASSOCIATED(diffuse%mmm_xti%data) &
          .OR. diffuse%mmm_xti%size .ne. nj)THEN
        diffuse%mmm_xti = zero_vector(nj)
     ENDIF
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, FMT='("* multimode - eff ion thermal diffusivity,M^2/sec")')
     WRITE  (niterdb, 10) (diffuse%mmm_xti%data(j), j=1,nj) 
  ELSE
     diffuse%mmm_xti%size = nj                 ! if not present in state file set to zero
     diffuse%mmm_xti = zero_Vector(nj)         ! needed because we distribute the values
                                                   ! (may encounter end statement in read)
     READ   (niterdb,  7,END = 2021) starflag
     READ   (niterdb, 10,END = 2021) (work(j), j=1,nj)
     diffuse%mmm_xti%data(1:nj)  = work(1:nj)
  END IF
2021 CONTINUE

  IF (irwflag .EQ. 0) THEN
     IF(.NOT. ASSOCIATED(diffuse%mmm_xte%data) &
          .OR. diffuse%mmm_xte%size .ne. nj)THEN
        diffuse%mmm_xte = zero_vector(nj)
     ENDIF
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, FMT='("* multimode - eff electron thermal diffusivity,M^2/sec")')
     WRITE  (niterdb, 10) (diffuse%mmm_xte%data(j), j=1,nj) 
  ELSE
     diffuse%mmm_xte%size = nj                 ! if not present in state file set to zero
     diffuse%mmm_xte = zero_Vector(nj)         ! needed because we distribute the values
                                                   ! (may encounter end statement in read)
     READ   (niterdb,  7,END = 2022) starflag
     READ   (niterdb, 10,END = 2022) (work(j), j=1,nj)
     diffuse%mmm_xte%data(1:nj)  = work(1:nj)
  END IF
2022 CONTINUE

  IF (irwflag .EQ. 0) THEN
     IF(.NOT. ASSOCIATED(diffuse%mmm_xdz%data) &
          .OR. diffuse%mmm_xdz%size .ne. nj)THEN
        diffuse%mmm_xdz = zero_vector(nj)
     ENDIF
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, FMT='("* multimode - impurity ion  diffusivity Weiland model,M^2/sec")')
     WRITE  (niterdb, 10) (diffuse%mmm_xdz%data(j), j=1,nj) 
  ELSE
     diffuse%mmm_xdz%size = nj                 ! if not present in state file set to zero
     diffuse%mmm_xdz = zero_Vector(nj)         ! needed because we distribute the values
                                                   ! (may encounter end statement in read)
     READ   (niterdb,  7,END = 2023) starflag
     READ   (niterdb, 10,END = 2023) (work(j), j=1,nj)
     diffuse%mmm_xdz%data(1:nj)  = work(1:nj)
  END IF
2023 CONTINUE

  IF (irwflag .EQ. 0) THEN
     IF(.NOT. ASSOCIATED(diffuse%mmm_xvt%data) &
        .OR. diffuse%mmm_xvt%size .ne. nj)THEN
        diffuse%mmm_xvt = zero_vector(nj)
     ENDIF
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, FMT='("* multimode - toroidal momentum transport, Weiland model  diffusivity,M^2/sec")')
     WRITE  (niterdb, 10) (diffuse%mmm_xvt%data(j), j=1,nj) 
  ELSE
     diffuse%mmm_xvt%size = nj                 ! if not present in state file set to zero
     diffuse%mmm_xvt = zero_Vector(nj)         ! needed because we distribute the values
                                                   ! (may encounter end statement in read)
     READ   (niterdb,  7,END = 2024) starflag
     READ   (niterdb, 10,END = 2024) (work(j), j=1,nj)
     diffuse%mmm_xvt%data(1:nj)  = work(1:nj)
  END IF
2024 CONTINUE

  IF (irwflag .EQ. 0) THEN
     IF(.NOT. ASSOCIATED(diffuse%mmm_xvp%data) &
          .OR. diffuse%mmm_xvp%size .ne. nj)THEN
        diffuse%mmm_xvp = zero_vector(nj)
     ENDIF
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, FMT='("* multimode - poloidal  momentum transport, Weiland model  diffusivity,M^2/sec")')
     WRITE  (niterdb, 10) (diffuse%mmm_xvp%data(j), j=1,nj) 
  ELSE
     diffuse%mmm_xvp%size = nj                 ! if not present in state file set to zero
     diffuse%mmm_xvp = zero_Vector(nj)         ! needed because we distribute the values
                                                   ! (may encounter end statement in read)
     READ   (niterdb,  7,END = 2025) starflag
     READ   (niterdb, 10,END = 2025) (work(j), j=1,nj)
     diffuse%mmm_xvp%data(1:nj)  = work(1:nj)
  END IF
2025 CONTINUE

  IF (irwflag .EQ. 0) THEN
     IF(.NOT. ASSOCIATED(diffuse%mmm_xtiW20%data) &
        .OR. diffuse%mmm_xtiw20%size .ne. nj)THEN
        diffuse%mmm_xtiW20 = zero_vector(nj)
     ENDIF
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, FMT='("* multimode - ion thermal diffusivity,M^2/sec, Weiland model part")')
     WRITE  (niterdb, 10) (diffuse%mmm_xtiW20%data(j), j=1,nj) 
  ELSE
     diffuse%mmm_xtiW20%size = nj                 ! if not present in state file set to zero
     diffuse%mmm_xtiW20 = zero_Vector(nj)         ! needed because we distribute the values
                                                   ! (may encounter end statement in read)
     READ   (niterdb,  7,END = 2026) starflag
     READ   (niterdb, 10,END = 2026) (work(j), j=1,nj)
     diffuse%mmm_xtiW20%data(1:nj)  = work(1:nj)
  END IF
2026 CONTINUE

  IF (irwflag .EQ. 0) THEN
     IF(.NOT. ASSOCIATED(diffuse%mmm_xdiW20%data) &
          .OR. diffuse%mmm_xdiw20%size .ne. nj)THEN
        diffuse%mmm_xdiW20 = zero_vector(nj)
     ENDIF
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, FMT='("* multimode - particle  diffusivity,M^2/sec, Weiland model part")')
     WRITE  (niterdb, 10) (diffuse%mmm_xdiW20%data(j), j=1,nj) 
  ELSE
     diffuse%mmm_xdiW20%size = nj                 ! if not present in state file set to zero
     diffuse%mmm_xdiW20 = zero_Vector(nj)         ! needed because we distribute the values
                                                   ! (may encounter end statement in read)
     READ   (niterdb,  7,END = 2027) starflag
     READ   (niterdb, 10,END = 2027) (work(j), j=1,nj)
     diffuse%mmm_xdiW20%data(1:nj)  = work(1:nj)
  END IF
2027 CONTINUE

  IF (irwflag .EQ. 0) THEN
     IF(.NOT. ASSOCIATED(diffuse%mmm_xteW20%data) &
          .OR. diffuse%mmm_xteW20%size .ne. nj)THEN
        diffuse%mmm_xteW20 = zero_vector(nj)
     ENDIF
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, FMT='("* multimode - electron  thermal diffusivity,M^2/sec, Weiland model part")')
     WRITE  (niterdb, 10) (diffuse%mmm_xteW20%data(j), j=1,nj) 
  ELSE
     diffuse%mmm_xteW20%size = nj                 ! if not present in state file set to zero
     diffuse%mmm_xteW20 = zero_Vector(nj)         ! needed because we distribute the values
                                                   ! (may encounter end statement in read)
     READ   (niterdb,  7,END = 2028) starflag
     READ   (niterdb, 10,END = 2028) (work(j), j=1,nj)
     diffuse%mmm_xteW20%data(1:nj)  = work(1:nj)
  END IF
2028 CONTINUE

  IF (irwflag .EQ. 0) THEN
     IF(.NOT. ASSOCIATED(diffuse%mmm_xtiDBM%data) &
          .OR. diffuse%mmm_xtiDBM%size .ne. nj)THEN
        diffuse%mmm_xtiDBM = zero_vector(nj)
     ENDIF
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, FMT='("* multimode - ion thermal diffusivity,M^2/sec, DRB  model part")')
     WRITE  (niterdb, 10) (diffuse%mmm_xtiDBM%data(j), j=1,nj) 
  ELSE
     diffuse%mmm_xtiDBM%size = nj                 ! if not present in state file set to zero
     diffuse%mmm_xtiDBM = zero_Vector(nj)         ! needed because we distribute the values
                                                  ! (may encounter end statement in read)
     READ   (niterdb,  7,END = 2029) starflag
     READ   (niterdb, 10,END = 2029) (work(j), j=1,nj)
     diffuse%mmm_xtiDBM%data(1:nj)  = work(1:nj)
  END IF
2029 CONTINUE

  IF (irwflag .EQ. 0) THEN
     IF(.NOT. ASSOCIATED(diffuse%mmm_xdiDBM%data) &
          .OR. diffuse%mmm_xdiDBM%size .ne. nj)THEN
        diffuse%mmm_xdiDBM = zero_vector(nj)
     ENDIF
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, FMT='("* multimode - particle  diffusivity,M^2/sec, DRB model part")')
     WRITE  (niterdb, 10) (diffuse%mmm_xdiDBM%data(j), j=1,nj) 
  ELSE
     diffuse%mmm_xdiDBM%size = nj                 ! if not present in state file set to zero
     diffuse%mmm_xdiDBM = zero_Vector(nj)         ! needed because we distribute the values
                                                   ! (may encounter end statement in read)
     READ   (niterdb,  7,END = 2030) starflag
     READ   (niterdb, 10,END = 2030) (work(j), j=1,nj)
     diffuse%mmm_xdiDBM%data(1:nj)  = work(1:nj)
  END IF
2030 CONTINUE


  IF (irwflag .EQ. 0) THEN
     IF(.NOT. ASSOCIATED(diffuse%mmm_xteDBM%data) &
          .OR. diffuse%mmm_xteDBM%size .ne. nj)THEN
        diffuse%mmm_xteDBM = zero_vector(nj)
     ENDIF
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, FMT='("* multimode - electron  thermal diffusivity,M^2/sec, DRBmodel part")')
     WRITE  (niterdb, 10) (diffuse%mmm_xteDBM%data(j), j=1,nj) 
  ELSE
     diffuse%mmm_xteDBM%size = nj                 ! if not present in state file set to zero
     diffuse%mmm_xteDBM = zero_Vector(nj)         ! needed because we distribute the values
                                                   ! (may encounter end statement in read)
     READ   (niterdb,  7,END = 2031) starflag
     READ   (niterdb, 10,END = 2031) (work(j), j=1,nj)
     diffuse%mmm_xteDBM%data(1:nj)  = work(1:nj)
  END IF
2031 CONTINUE

  IF (irwflag .EQ. 0) THEN
     IF(.NOT. ASSOCIATED(diffuse%mmm_xteETG%data) &
          .OR. diffuse%mmm_xteETG%size .ne. nj)THEN
        diffuse%mmm_xteETG = zero_vector(nj)
     ENDIF
     IF (iterdsc .NE. 0) &
          WRITE  (niterdb, FMT='("* multimode - electron  thermal diffusivity,M^2/sec, ETG model part")')
     WRITE  (niterdb, 10) (diffuse%mmm_xteETG%data(j), j=1,nj) 
  ELSE
     diffuse%mmm_xteETG%size = nj                 ! if not present in state file set to zero
     diffuse%mmm_xteETG = zero_Vector(nj)         ! needed because we distribute the values
                                                   ! (may encounter end statement in read)
     READ   (niterdb,  7,END = 2032) starflag
     READ   (niterdb, 10,END = 2032) (work(j), j=1,nj)
     diffuse%mmm_xteETG%data(1:nj)  = work(1:nj)
  END IF
2032 CONTINUE

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

  DO jj=1,ngW20
     IF (irwflag .EQ. 0) THEN
        IF( jj == 1)tlabel = '* multimode - growth rate most unstable ion mode Weiland positive freq direction [1/sec] '
        IF( jj == 2)tlabel = '* multimode - growth rate most unstable elec mode Weiland positive freq direction [1/sec]'
        IF( jj == 3)tlabel = '* multimode - growth rate most unstable ion mode Weiland negative freq direction [1/sec]'
        IF( jj == 4)tlabel = '* multimode - growth rate most unstable elec mode Weiland negative freq direction [1/sec]'
        tlabel = ADJUSTL(tlabel)
        IF (iterdsc .NE. 0 )&
             WRITE  (niterdb, FMT='(A)') tlabel(1:LEN_TRIM(tlabel))
        work(1:nj) = get_values(diffuse%mmm_gammaW20(jj))
        WRITE  (niterdb, 10) (work(j), j=1,nj)
     ELSE
     !Protected read for older statefiles
        READ   (niterdb,  7,END = 2033) starflag
        READ   (niterdb, 10,END = 2033) (work(j), j=1,nj)
        CALL delete_vector(diffuse%mmm_gammaW20(jj))
        diffuse%mmm_gammaW20(jj) = new_Vector(nj,work)
     END IF
2033 CONTINUE

  END DO

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

  DO jj=1,ngW20
     IF (irwflag .EQ. 0) THEN
        IF( jj == 1)tlabel = '* multimode - freq,rad/sec, most unstable ion mode Weiland positive freq direction'
        IF( jj == 2)tlabel = '* multimode - freq,rad/sec, most unstable elec mode Weiland positive freq direction'
        IF( jj == 3)tlabel = '* multimode - freq,rad/sec, most unstable ion mode Weiland negative freq direction'
        IF( jj == 4)tlabel = '* multimode - freq,rad/sec,most unstable elec mode Weiland negative freq direction '
        tlabel = ADJUSTL(tlabel)
        IF (iterdsc .NE. 0 )&
             WRITE  (niterdb, FMT='(A)') tlabel(1:LEN_TRIM(tlabel))
        work(1:nj) = get_values(diffuse%mmm_omegaW20(jj))
        WRITE  (niterdb, 10) (work(j), j=1,nj)
     ELSE
     !Protected read for older statefiles
        READ   (niterdb,  7,END = 2034) starflag
        READ   (niterdb, 10,END = 2034) (work(j), j=1,nj)
        CALL delete_vector(diffuse%mmm_omegaW20(jj))
        diffuse%mmm_omegaW20(jj) = new_Vector(nj,work)
     END IF
2034 CONTINUE

  END DO

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
  DO jj=1,ngW20
     IF (irwflag .EQ. 0) THEN
        IF( jj == 1)tlabel = '* multimode - total ion thermal flux,W/m^2:'
        IF( jj == 2)tlabel = '* multimode - total hydrogenic ion flux flux,1/(m^2 sec)'
        IF( jj == 3)tlabel = '* multimode - total electron thermal  flux flux,W/m^2'
        IF( jj == 4)tlabel = '* multimode - total impurity ion flux 1/(m^2 sec)'
        tlabel = ADJUSTL(tlabel)
        IF (iterdsc .NE. 0 )&
             WRITE  (niterdb, FMT='(A)') tlabel(1:LEN_TRIM(tlabel))
        work(1:nj) = get_values(diffuse%mmm_vflux(jj))
        WRITE  (niterdb, 10) (work(j), j=1,nj)
     ELSE
     !Protected read for older statefiles
        READ   (niterdb,  7,END = 2035) starflag
        READ   (niterdb, 10,END = 2035) (work(j), j=1,nj)
        CALL delete_vector(diffuse%mmm_vflux(jj))
        diffuse%mmm_vflux(jj) = new_Vector(nj,work)
     END IF
2035 CONTINUE

  END DO


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

  DO jj=1,ngW20
     IF (irwflag .EQ. 0) THEN
        IF( jj == 1)tlabel = '* multimode - ion thermal convective velocity,m/sec'
        IF( jj == 2)tlabel = '* multimode - hydrogenic ion particle  convective velocity,m/sec'
        IF( jj == 3)tlabel = '* multimode - electron thermal convective velocity,m/sec'
        IF( jj == 4)tlabel = '* multimode - impurity ion particle  convective velocity,m/sec'
        IF( jj == 5)tlabel = '* multimode - toroidal momentum pinch ,m/sec'
        IF( jj == 6)tlabel = '* multimode - poloidal momentum pinch ,m/sec'
  !NOTE jj=1,2,3 info is put out as zero by current version of multimode
        tlabel = ADJUSTL(tlabel)
        IF (iterdsc .NE. 0 )&
             WRITE  (niterdb, FMT='(A)') tlabel(1:LEN_TRIM(tlabel))
        work(1:nj) = get_values(diffuse%mmm_vconv(jj))
        WRITE  (niterdb, 10) (work(j), j=1,nj)
     ELSE
     !Protected read for older statefiles
        READ   (niterdb,  7,END = 2036) starflag
        READ   (niterdb, 10,END = 2036) (work(j), j=1,nj)
        CALL delete_vector(diffuse%mmm_vconv(jj))
        diffuse%mmm_vconv(jj) = new_Vector(nj,work)
     END IF
2036 CONTINUE

  END DO


  !rho grid for neutral beam (normally same as rho grid used above
  ! but some codes change grids on the fly so need this)
  IF(neut_beam%nbeams .GT. 0 ) THEN
     IF (irwflag .EQ. 0) THEN
        IF (iterdsc .NE. 0) &
             WRITE  (niterdb, FMT='("* rho grid for neutral beam")')
        WRITE  (niterdb, 10) (neut_beam%rhog_beam(j), j=1,nj) 
     ELSE
        neut_beam%rhog_beam(:)  = get_values(rho_grid)
        READ   (niterdb,  7,END = 2037) starflag
        READ   (niterdb, 10,END = 2037) (neut_beam%rhog_beam(j), j=1,nj)
     END IF
 ENDIF

2037 CONTINUE


  jp = 0
  ji = 0
  DO jj=1,ntot
     IF (jj .LE. nprim)  jp = jp + 1
     IF (jj .GT. nprim)  ji = ji + 1
     IF (irwflag .EQ. 0) THEN
        IF (iterdsc .NE. 0 .AND. jj .LE. nprim) &
             WRITE  (niterdb, 9017) namep(jp)
9017    FORMAT ('* convective  flux:  primary ion,', &
             ' #/(meter2 sec), species: ', a)
        IF (iterdsc .NE. 0 .AND. jj .GT. nprim .AND. jj .LE. nion) &
             WRITE  (niterdb, 9018) namei(ji)
9018    FORMAT ('* convective flux  impurity ion', &
             ' #/(meter**2sec), species: ', a)
        IF (iterdsc .NE. 0 .AND. jj .EQ.  nion+1 ) &
             WRITE  (niterdb, 9019)
9019    FORMAT('* convective flux,electron energy, Joules/m**2 sec')
        IF (iterdsc .NE. 0 .AND. jj .EQ.  nion+2 ) &
             WRITE  (niterdb, 9020)
9020    FORMAT('* convective flux ion energy flux, joules/m**2 sec')
        IF (iterdsc .NE. 0 .AND. jj .EQ.  nion+3 ) &
             WRITE  (niterdb, 9021)
9021    FORMAT('* convective flux in Faradys law Tesla/sec')
        IF (iterdsc .NE. 0 .AND. jj .EQ.  nion+4 ) &
             WRITE  (niterdb, 9022)
9022    FORMAT('* convective toroidal momentum flux Kg/sec^2 ')

        work(1:nj) = get_values(profile%flux_conv(jj))
        WRITE  (niterdb, 10) (work(j), j=1,nj)
     ELSE
        READ   (niterdb,  7,END = 2038) line
        READ   (niterdb, 10,END = 2038) (work(j)       , j=1,nj)
        CALL delete_Vector(profile%flux_conv(jj))  ! made it to here so delete and reassign
        profile%flux_conv(jj) = new_Vector(nj,work)
     END IF
  END DO

2038 CONTINUE

     IF (irwflag .EQ. 0) THEN
        IF (iterdsc .NE. 0) &
             WRITE  (niterdb, FMT='("* electron particle flux zone edge grid #/(M^2 sec)")')
        WRITE  (niterdb, 10) (profile%fluxe%data(j), j=1,nj) 
     ELSE
        profile%fluxe   = zero_Vector(nj)
        READ   (niterdb,  7,END = 2039) starflag
        READ   (niterdb, 10,END = 2039) (profile%fluxe%data(j), j=1,nj)
     END IF

2039 CONTINUE

     IF (irwflag .EQ. 0) THEN
        IF (iterdsc .NE. 0) &
             WRITE  (niterdb, FMT='("* total ion particle flux,zone edge grid #/(M^2 sec)")')
        WRITE  (niterdb, 10) (profile%fluxi%data(j), j=1,nj) 
     ELSE
        profile%fluxi   = zero_Vector(nj)
        READ   (niterdb,  7,END = 2040) starflag
        READ   (niterdb, 10,END = 2040) (profile%fluxi%data(j), j=1,nj)
     END IF

2040 CONTINUE

     IF (irwflag .EQ. 0) THEN
        IF (iterdsc .NE. 0) &
             WRITE  (niterdb, FMT='("* tglf toroidal momentum flux,electrons kg/sec^2")')
        WRITE  (niterdb, 10) (tglf_m_output(j,1), j=1,nj) 
        IF (iterdsc .NE. 0) &
             WRITE  (niterdb, FMT='("* tglf toroidal momentum flux,prim ions kg/sec^2")')
        WRITE  (niterdb, 10) (tglf_m_output(j,2), j=1,nj) 
        IF (iterdsc .NE. 0) &
             WRITE  (niterdb, FMT='("* tglf toroidal momentum flux,imp ions kg/sec^2")')
        WRITE  (niterdb, 10) (tglf_m_output(j,3), j=1,nj) 
     ELSE
        IF(ALLOCATED(tglf_m_output))DEALLOCATE(tglf_m_output)
        ALLOCATE(tglf_m_output(nj,3))
        tglf_m_output(:,:)  = zeroc
        READ   (niterdb,  7,END = 2050) starflag
        READ   (niterdb, 10,END = 2050) (tglf_m_output(j,1),j=1,nj)
        READ   (niterdb,  7,END = 2050) starflag
        READ   (niterdb, 10,END = 2050) (tglf_m_output(j,2),j=1,nj)
        READ   (niterdb,  7,END = 2050) starflag
        READ   (niterdb, 10,END = 2050) (tglf_m_output(j,3),j=1,nj)
     END IF

2050 CONTINUE

     IF (irwflag .EQ. 0) THEN
     IF(.NOT. ASSOCIATED(profile%er_tot_nclass%data))THEN
        profile%er_tot_nclass = zero_vector(nj)
     ENDIF
        IF (iterdsc .NE. 0) &
             WRITE  (niterdb, FMT='("* Nclass total radial electric field,v/m")')
        WRITE  (niterdb, 10) (profile%er_tot_nclass%data(j), j=1,nj) 
     ELSE
        profile%er_tot_nclass = zero_Vector(nj)
        READ   (niterdb,  7,END = 2051) starflag
        READ   (niterdb, 10,END = 2051) (profile%er_tot_nclass%data(j),j=1,nj)
     ENDIF

2051 CONTINUE
     !------------------------------------------------------------
     ! put out common frequencies:
     ! If these are not defined then define them and set to zero
     ! frequencies that are not flux functions are on the rho grid at zma=0
     ! on the outboard side
     !-----------------------------------------------------------------
 
     IF(.NOT. ASSOCIATED(plasma_frequencies%omega_pi)) &
          CALL allocate_plasma_freq

     IF (irwflag .EQ. 0) THEN
        IF (iterdsc .NE. 0) &
             WRITE  (niterdb, FMT='("* frequencies  outboard at zma,not flux function")')
        k=0
        DO i=1,nprim
           k = k+1
           IF(namep(i) == 'dt')THEN  ! first d :
              WRITE  (niterdb,FMT='("* ion plasma freq,rad/sec species ",a)')'D'
              WRITE  (niterdb, 10) (plasma_frequencies%omega_pi(k)%data(j), j=1,nj) 

              WRITE  (niterdb,FMT='("* ion cyclotron freq,rad/sec species ",a)')'D'
              WRITE  (niterdb, 10) (plasma_frequencies%omega_ci(k)%data(j), j=1,nj)

              WRITE  (niterdb,FMT='("* ion lower hybrid freq,rad/sec species ",a)')'D'
              WRITE  (niterdb, 10) (plasma_frequencies%omega_lh(k)%data(j), j=1,nj)

              WRITE  (niterdb,FMT='("* ion upper hybrid freq,rad/sec species ",a)')'D'
              WRITE  (niterdb, 10) (plasma_frequencies%omega_uh(k)%data(j), j=1,nj)


              k=k+1                  ! then t:
              WRITE  (niterdb,FMT='("* ion plasma freq,rad/sec species ",a)')'T'
              WRITE  (niterdb, 10) (plasma_frequencies%omega_pi(k)%data(j), j=1,nj) 

              WRITE  (niterdb,FMT='("* ion cyclotron freq,rad/sec species ",a)')'T'
              WRITE  (niterdb, 10) (plasma_frequencies%omega_ci(k)%data(j), j=1,nj) 

              WRITE  (niterdb,FMT='("* ion lower hybrid freq,rad/sec species ",a)')'T'
              WRITE  (niterdb, 10) (plasma_frequencies%omega_lh(k)%data(j), j=1,nj) 

              WRITE  (niterdb,FMT='("* ion upper hybrid freq,rad/sec species ",a)')'T'
              WRITE  (niterdb, 10) (plasma_frequencies%omega_uh(k)%data(j), j=1,nj)
           ELSE
              WRITE  (niterdb,FMT='("* ion plasma freq,rad/sec species ",a)')namep(i)
              WRITE  (niterdb, 10) (plasma_frequencies%omega_pi(k)%data(j), j=1,nj)

              WRITE  (niterdb,FMT='("* ion cyclotron freq,rad/sec species ",a)')namep(i) 
              WRITE  (niterdb, 10) (plasma_frequencies%omega_ci(k)%data(j), j=1,nj)
 
              WRITE  (niterdb,FMT='("* lower hybrid freq,rad/sec species ",a)')namep(i) 
              WRITE  (niterdb, 10) (plasma_frequencies%omega_lh(k)%data(j), j=1,nj) 

              WRITE  (niterdb,FMT='("* upper hybrid freq,rad/sec species ",a)')namep(i) 
              WRITE  (niterdb, 10) (plasma_frequencies%omega_uh(k)%data(j), j=1,nj)

           ENDIF
        ENDDO

        WRITE  (niterdb,FMT='("* electron cyclotron freq,rad/sec ")')
        WRITE  (niterdb, 10) (plasma_frequencies%omega_ce%data(j), j=1,nj)

        WRITE  (niterdb,FMT='("* electron plasma  freq,rad/sec")')
        WRITE  (niterdb, 10) (plasma_frequencies%omega_pe%data(j), j=1,nj) 
     ELSE
        READ   (niterdb,  7,END = 2052) starflag
        k=0
        DO i=1,nprim
           k = k+1
           IF(namep(i) == 'dt')THEN  ! first d :
              READ   (niterdb,  7,END  = 2052) starflag
              READ   (niterdb, 10,END  = 2052) (plasma_frequencies%omega_pi(k)%data(j), j=1,nj) 

              READ   (niterdb,  7,END  = 2052) starflag
              READ   (niterdb, 10,END=2052) (plasma_frequencies%omega_ci(k)%data(j), j=1,nj)

              READ   (niterdb,  7,END  = 2052) starflag
              READ   (niterdb, 10,END  = 2052) (plasma_frequencies%omega_lh(k)%data(j), j=1,nj)

              READ   (niterdb,  7,END  = 2052) starflag
              READ   (niterdb, 10,END  = 2052) (plasma_frequencies%omega_uh(k)%data(j), j=1,nj)


              k=k+1                  ! then t:
              READ   (niterdb,  7,END  = 2052) starflag
              READ   (niterdb, 10,END  = 2052) (plasma_frequencies%omega_pi(k)%data(j), j=1,nj) 

              READ   (niterdb,  7,END  = 2052) starflag
              READ   (niterdb, 10,END  = 2052) (plasma_frequencies%omega_ci(k)%data(j), j=1,nj) 

              READ   (niterdb,  7,END  = 2052) starflag
              READ   (niterdb, 10,END  = 2052) (plasma_frequencies%omega_lh(k)%data(j), j=1,nj) 

              READ   (niterdb,  7,END  = 2052) starflag
              READ   (niterdb, 10,END  = 2052) (plasma_frequencies%omega_uh(k)%data(j), j=1,nj)
           ELSE
              READ   (niterdb,  7,END  = 2052) starflag
              READ   (niterdb, 10,END  = 2052) (plasma_frequencies%omega_pi(k)%data(j), j=1,nj)

              READ   (niterdb,  7,END  = 2052) starflag
              READ   (niterdb, 10,END  = 2052) (plasma_frequencies%omega_ci(k)%data(j), j=1,nj)
 
              READ   (niterdb,  7,END  = 2052) starflag
              READ   (niterdb, 10,END  = 2052) (plasma_frequencies%omega_lh(k)%data(j), j=1,nj) 

              READ   (niterdb,  7,END  = 2052) starflag
              READ   (niterdb, 10,END  = 2052) (plasma_frequencies%omega_uh(k)%data(j), j=1,nj)

           ENDIF
        ENDDO

        READ   (niterdb,  7,END  = 2052) starflag
        READ   (niterdb, 10,END  = 2052) (plasma_frequencies%omega_ce%data(j), j=1,nj)

        READ   (niterdb,  7,END  = 2052) starflag
        READ   (niterdb, 10,END  = 2052) (plasma_frequencies%omega_pe%data(j), j=1,nj) 
        IF(k .ne. nprimp1)THEN
           k = k+1
           plasma_frequencies%omega_pi(k) = zero_vector(nj)
           plasma_frequencies%omega_ci(k) = zero_vector(nj)
           plasma_frequencies%omega_lh(k) = zero_vector(nj)
           plasma_frequencies%omega_uh(k) = zero_vector(nj)
        ENDIF
     ENDIF
2052 CONTINUE

!-------------------------------------------------------------------------
! -- put out glf eff chi values these are defined as the fluxes
! -- divided by the gradients
!-------------------------------------------------------------------------
    IF (irwflag .EQ. 0) THEN
        IF(.NOT. ALLOCATED(glf_eff_chi_output))THEN
           ALLOCATE(glf_eff_chi_output(ntot,nj)) 
           glf_eff_chi_output(:,:) = zeroc ! print out zeroes if not set
        ENDIF



        do j=1,nprim
           WRITE(niterdb,FMT='("* Effective glf prim ion diffusivity M^2/sec species: ",a)')namep(j)
           WRITE(niterdb,10)(glf_eff_chi_output(i,j),i=1,nj)
        ENDDO

        DO j=nprim+1,nion
           WRITE(niterdb,FMT='("* Effective glf imp ion diffusivity M^2/sec species: ",a)')namei(j-nprim)
           WRITE(niterdb,10)(glf_eff_chi_output(i,j),i=1,nj)
        ENDDO

        j= nion+1
        WRITE(niterdb,FMT='("* Effective glf elect energy  diffusivity M^2/sec")')
        WRITE(niterdb,10)(glf_eff_chi_output(i,j),i=1,nj)
        j=nion+2
        WRITE(niterdb,FMT='("* Effective glf ion energy  diffusivity M^2/sec")')
        WRITE(niterdb,10)(glf_eff_chi_output(i,j),i=1,nj)
        j=nion+3
        WRITE(niterdb,FMT='("* Effective glf empty slot for Bp)")')
        WRITE(niterdb,10)(glf_eff_chi_output(i,j),i=1,nj)
        j=nion+4
        WRITE(niterdb,FMT='("* Effective glf momentum  diffusivity M^2/sec")')
        WRITE(niterdb,10)(glf_eff_chi_output(i,j),i=1,nj)
     ELSE
        IF(ALLOCATED(glf_eff_chi_output))DEALLOCATE(glf_eff_chi_output)
        ALLOCATE(glf_eff_chi_output(ntot,nj))
        DO j=1,nion+4
           READ   (niterdb,  7,END = 2053) starflag
           READ   (niterdb,  10,END = 2053)(glf_eff_chi_output(i,j),i=1,nj)
        ENDDO
     ENDIF

2053 CONTINUE

!------------------------------------------------------------------------
! -- Aneutronic fusion reaction rates d(d,p)t,he3(d,p)he4
! -- see neutron rates for ddn reactions (fus_prod%..) above
! -- NOTE: gcnmp currently does not have these  rates
! -- NOTE d(d,p)t is  thermal-thermal rates only at this time
!------------------------------------------------------------------------

      IF (irwflag .EQ. 0) THEN 
        IF(.NOT. ASSOCIATED(fus_prod%p_ddt%data))THEN 
           fus_prod%p_ddt = zero_vector(nj)
           fus_prod%total_p_ddt = zeroc
        ENDIF
        IF(.NOT. ASSOCIATED(fus_prod%p_he3d_th%data))THEN 
           fus_prod%p_he3d_th = zero_vector(nj)
           fus_prod%total_p_he3d_th = zeroc
        ENDIF
        IF(.NOT. ASSOCIATED(fus_prod%p_he3d_beam_thermal%data))THEN 
           fus_prod%p_he3d_beam_thermal = zero_vector(nj)
           fus_prod%total_p_he3d_beam_thermal = zeroc
        ENDIF

        WRITE(niterdb,FMT='("* d(d,p)t reaction rate, #/(M^3 sec)",a)')
        WRITE(niterdb,10)(fus_prod%p_ddt%data(i),i=1,nj)
        WRITE(niterdb,FMT='("* d(d,p)t total reaction rate, #/sec",a)')
        WRITE(niterdb,10)fus_prod%total_p_ddt


        WRITE(niterdb,FMT='("* he3(d,p)he4 total,thermal,beam_thermal, reaction rates, #/sec",a)')
        WRITE(niterdb,10)fus_prod%total_p_he3d, fus_prod%total_p_he3d_th,fus_prod%total_p_he3d_beam_thermal  
        WRITE(niterdb,FMT='("* he3(d,p)he4 thermal-thermal reaction rate, #/(M^3 sec",a)')
        WRITE(niterdb,10)(fus_prod%p_he3d_th%data(i),i=1,nj)

        WRITE(niterdb,FMT='("* d(he3,p)he4 beam-thermal  summed over all beams  reaction rate, #/M^3 sec",a)')
        WRITE(niterdb,10)(fus_prod%p_he3d_beam_thermal%data(i),i=1,nj) ! beam-thermal
        WRITE(niterdb,FMT='("* d beam spin polarization, he3_frac for he3(d,p)he4 reactions")')
        WRITE(niterdb,10)neut_beam%d_beam_spin_pol,he3_frac,he3_thermal_spin_pol

     ELSE
        ok = delete_vector_nf(fus_prod%p_ddt)
        ok = delete_vector_nf(fus_prod%p_he3d_th)
        ok = delete_vector_nf(fus_prod%p_he3d_beam_thermal)
        fus_prod%p_ddt                     =  zero_Vector(nj)
        fus_prod%p_he3d_th                 =  zero_Vector(nj)
        fus_prod%p_he3d_beam_thermal       =  zero_Vector(nj)
        fus_prod%total_p_he3d              = zeroc
        fus_prod%total_p_he3d_th           = zeroc
        fus_prod%total_p_he3d_beam_thermal = zeroc

        READ   (niterdb,  7,END = 2054) starflag
        READ   (niterdb,  10,END = 2054)(fus_prod%p_ddt%data(j),j=1,nj)
        READ   (niterdb,  7,END = 2054) starflag
        READ   (niterdb,  10,END = 2054)fus_prod%total_p_ddt

        READ   (niterdb,  7,END = 2054) starflag
        READ   (niterdb,  10,END = 2054)fus_prod%total_p_he3d,fus_prod%total_p_he3d_th,fus_prod%total_p_he3d_beam_thermal
        READ   (niterdb,  7,END = 2054) starflag
        READ   (niterdb,  10,END = 2054)(fus_prod%p_he3d_th%data(j),j=1,nj)

        READ   (niterdb,  7,END = 2054) starflag
        READ   (niterdb,  10,END = 2054)(fus_prod%p_he3d_beam_thermal%data(j),j=1,nj)
        READ   (niterdb,  7,END = 2054) starflag
        READ   (niterdb,  10,END = 2054)neut_beam%d_beam_spin_pol,he3_frac,he3_thermal_spin_pol

     ENDIF

2054 CONTINUE


  CLOSE (unit = niterdb, status = 'KEEP')
  IF(ALLOCATED(work))DEALLOCATE(work)
  IF(ALLOCATED(worknpsi))DEALLOCATE(worknpsi)
  
  RETURN
 END SUBROUTINE iter_dbase_txt




SUBROUTINE check_monotonic (array, n, monotonic, incr)
  !
  USE nrtype,   ONLY : I4B,DP
  IMPLICIT  NONE
  INTEGER(I4B), INTENT(IN) :: n,incr
  INTEGER(I4B) j
  REAL(DP),INTENT(IN),DIMENSION(:) :: array(n)
  LOGICAL, INTENT(OUT) ::    monotonic


  !
  ! ----------------------------------------------------------------------
  ! incr = +1 check for monotonic increasing
  ! incr = -1                     decreasing
  ! return monotonic = .true.  if array is monotonic
  ! return             .false. otherwise
  ! ----------------------------------------------------------------------
  !
  !

  monotonic = .TRUE.
  IF (incr .GT. 0) THEN    ! check for increasing array
     DO j=1,n-1
        IF (array(j) .GE. array(j+1))  monotonic = .FALSE.
     END DO
  ELSE                     ! check for decreasing array
     DO j=1,n-1
        IF (array(j) .LE. array(j+1))  monotonic = .FALSE.
     END DO
  END IF

  RETURN
 
END SUBROUTINE check_monotonic
