      subroutine iter_dbase1
c

c
c ----------------------------------------------------------------------
c --- collects and writes out (or reads in) data for ITER database
c ---
c --- similar to routine ITER_DBASE except that output file pointed
c --- by iterdbfilename is in netCDF format and file is "iterdb.nc"
c ---
c ----------------------------------------------------------------------
c
c --- INPUT
c
c INCLUDE file constnts.i
c     psimks           conversion factor convers psi from kg/cm**2
c                      to volt*second/radian
c
c INCLUDE file iterdb.i:
c     irwflag  read/write flag, included to make the reading
c              of the file created by this subroutine
c              consistent with the writing of the file.
c              irwflag = 0  ==>  WRITE the data in netcdf form
c              irwflag = 1  ==>  READ  the data from text file ??? - do not use
c              if this subroutine is used in read mode note that
c              the relevant arrays that are read in are defined
c              in INCLUDE files as described below.
c              NOTE IF IRWFLAG=1, THE I/O UNIT NUMBER ASSIGNED
c              TO THE FILE FOR READING PURPOSES, NITERDB, IS
c              CURRENTLY HARD CODED AT 11 !!!!!!!!!!!!!!!!!!
c
c     iterdb   if set to -1 then this subroutine will execute.
c              Otherwise this subroutine returns immediately without
c              performing any calculations. (by doing it this way
c              the necessary INCLUDE files can all be brought in
c              at this subroutine, rather than at the calling level)
c
c     iterdbfilename  name of file which will be created or appended to
c     iterdsc         switch used to output description of parameters
c     niterdb         Fortran I/O unit number for output OR INPUT
c
c INCLUDE file contour.i:
c     nplasbdry
c     rplasbdry(j)
c     zplasbdry(j)     j=1,2..nplasbdry describes the plasma boundary
c
c INCLUDE file extra.i:
c     q(j)             j=1,2,...nj, the safety factor on rho grid
c     betap
c     beta
c
c INCLUDE file etc.i:
c     totcur            total plasma current, amps
c
c INCLUDE file io.i
c     ncrt
c     nout             unit numbers for diagnostic output
c INCLUDE file psig.i:
c     rho(j)           j=1,2,..npsi the rho values corresponding
c                      to the psi grid in meters
c     qpsi(j)          the corresponding q values
c     psivolp(j)       volume, meter**3, enclosed by flux surface
c     grho1(j)         < ABS (grad rho)>
c     grho2(j)         <(grad rho)**2>
c     rmajavnpsi(j)    j=1,2,..npsi average major radius
c                      at elevation of magnetic axis
c     rminavnpsi(j)    same for minor radius
c     psivolp(j)       plasma volume
c
c     triangnpsi(j)    upper triangularity
c     triangnpsi_l(j)    lower triangularity
c     pindentnpsi(j)   indentation
c     sfareanpsi(j)    surface area of flux surfaces
c     cxareanpsi(j)    cross-sectional area
c
c     btor             mag field at R0, in vacuum, tesla
c     elong(j)         elongation
c
c INCLUDE file geom.i:
c     fcap(j)
c     gcap(j)
c     hcap(j)          the geometric factors
c     ali              plasma inductance
c
c INCLUDE file ions.i:
c     namep(i)         i=1,..nprim
c     namei(i)         i=1,..nimp
c     namen(i)         i=1,..nneu
c     zeff(j)          j=1,2...nj
c
c INCLUDE file nub.i:
c     ibion            index (into namep) of beam species
c     sbion
c
c INCLUDE file nub2.i:
c     enbeam(j)        j=1,2..nj  fast ion density due to beam
c
c INCLUDE file neut.i:
c     enn(j,k)         j=1,2..nj(njs),k=1,2 neutral density
c     ennw             neutral density due to wall source
c     ennv                                    volume
c     volsn            volume source of neutrals
c
c INCLUDE file machin.i:
c     rmajor
c     kappa
c     btor
c
c INCLUDE file mesh.i
c     r(j)             j=1,2...nj,rho grid
c
c INCLUDE file solcon.i
c     time            current time, sec
c
c INCLUDE file numbrs.i
c     nj              transport grid size
c     nprim
c     nneu
c     nimp
c     nion
c
c INCLUDE file rhog.i
c     psir(j)         j=1,2...nj psi on r (i.e., rho) grid
c
c INCLUDE file soln.i
c     te(j)           electron and ion temperatures,keV
c     ti(j)
c     ene(j)          electron density
c     en(j,k)         j=1,2..nj,k=1..nion,where nion=nprim+nimp is the
c                     number of primary plus impurity ion species
c     curden(j)       current  density
c     rbp(j)          rho*bp field
c
c INCLUDE file sourc.i:
c     curboot(j)
c     curohm(j)
c     currf(j)
c     curbe(i)
c     curbi(j)
c     qrad(j)         radiated power
c     totohm
c     totbeam
c     totboot
c     totrf           total currents, amps
c     dpedtc(j), dpidtc(j)
c     qconde,qcondi
c     qconve,qconvi
c     qdelt
c     qexch
c     qione,qioni
c     qbeame,qbeami
c     qrfe, qrfi
c     qe2d, qi2d
c     qtfuse,qtfusi
c     qmag
c     qsawe,qsawi
c     qbfusi
c     sbeam
c
c INCLUDE file tcoef.i:
c     cheinv(j)
c     chiinv(j)        electron and ion thermal cond.
c     xkineo(j)        ion neoclassical thermal conductivity
c
c INCLUDE file tordlrot.i:
c     angrot(j)       angular rotation speed, radian/second
c     storqueb(j)     beam torque density nt-m/m**3
c
c INCLUDE file cer.i
c     Kpol_c(j)       cer impurity poloidal velocity/Bpol
c     Kpol_d(j)       main ion     poloidal velocity/Bpol
c     angrot_c(j)     cer impurity toroidal rotation radians/second
c     angrot_d(j)     main ion     toroidal rotation radians/second
c     ave_vpar_d(j)   main ion flux surface average parallel velocity
c                     meters/second
c     udia_d(j)       main ion diamagnetic velocity/Bpol
c     udia_c(j)       cer  ion diamagnetic velocity/Bpol
c     ugrt(j)         ion temperature gradient/(e*R*Bpol)
c     Epsi(j)         Er/RBp poloidal flux derivative of potential
c                     from NCLASS
c     Epsi_exp(j)     Er/RBp from cer data
c     sqz_d(j)        main ion orbit squeeze factor
c     cer_bp(j)       Bp   along Z=0 cord Tesla
c     cer_btdr(j)     Bt/R along Z=0 cord Tesla/meter
c
c INCLUDE file  yoka.i
c     ishot           shot number
c
c --- OUTPUT
c
c     ALL OUTPUT IS TO FILE POINTED TO BY ITERDBFILENAME
c     THIS ROUTINE DOES NOT (AND SHOULD NOT IN THE FUTURE)
c     PERFORM ANY CALCULATIONS WHICH DEFINE NEW PARAMETERS
c     BECAUSE THIS IS ONLY AN OUTPUT ROUTINE!
c
c ----------------------------------------------------- 11/22/94 --- HSJ
c
      USE typeSizes                  !part of netcdf 
      USE netcdf                     !ditto

      USE param
      USE ions
      USE neut
      USE solcon
      USE soln
      USE io
      USE contour
      USE mhdpar
      USE nub  
      USE nub2
      USE extra
      USE yoka
      USE numbrs
      USE mesh
      USE sourc
      USE machin
      USE geom
      USE tordlrot
      USE constnts
      USE tcoef
      USE bd_condtn, only : totcur
      USE psig
      USE rhog
      USE flxav
      USE iterdbmd
      USE etc
      USE cer
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      character rcs_id*63
      save      rcs_id
      data      rcs_id /
     ."$Id: cray205.f,v 1.33 2010/07/30 16:40:07 stjohn Exp $"/

c      include 'netcdf.inc'   ! from the netCDF package..
c                            ..this must be symlinked to local directory
c
c     F90 automatic (temporary) arrays:
      dimension  aspline(kpsi), bspline(kpsi), cspline(kpsi),
     .           cs2spline(kpsi,3), dspline(kpsi), espline(kpsi),
     .           fspline(kpsi), bpar(4), work(kj),eni(kpsi,kimp),
     .           tmpdata(kpsi)
      logical    monotonic
      character  starflag*2, title*39, headerline*132
c
      integer    rcode                                ! error code
      integer    dim_njtime(2), dim_plasbdrytime(2)   ! variable shape
      integer    dim_njprimtime(3), dim_njimptime(3), dim_njneutime(3)
c
c     matrix corner, initialize right before the first usage
c
      integer    c1t(2), cnt(2), cbt(2)
      integer    c11t(3), cnpt(3), cnit(3), cnnt(3)
      integer    cmap(3)                              ! mapping vector
c
c      equivalence (aspline(1)    , xdum(1))
c      equivalence (bspline(1)    , xdum(kpsi+1))
c      equivalence (cspline(1)    , wdum(1))
c      equivalence (dspline(1)    , vdum(1))
c      equivalence (espline(1)    , ydum(1))
c      equivalence (fspline(1)    , ydum(kpsi+1))
c      equivalence (work(1)       , zdum(1))
c      equivalence (cs2spline(1,1), sdum(1))
c      equivalence (bpar(1)       , rdum(1))
c
      data icall/0/, convert/1.60217733e-10/ ! keV/sec/cm*3 > watts/m**3
      save icall   , convert
c
      if (iterdb .ne. -1)  return
c     iterdsc = 1                   ! always write descriptor
      call NCPOPT (NCVERBOS)        ! if an error occurs let netCDF
c                                     write an error message but tell
c                                     netCDF not to consider it fatal
    7 format (         a   )
    8 format (5(2x,    a  ))        ! common character write/read format
    9 format (5(2x,   i6  ))        ! common integer   write/read format
   10 format (5(2x,1pe14.4))        ! common floating  write/read format
c
c     OPENing method depends on whether read or write,
c     and whether or not first time
c
      if (irwflag .eq. 0) then    ! WRITE to file, might be existing
        if (icall .eq. 0) then    ! first time, so open even file exists
          niterdb = NCCRE ('iterdb.nc', NCCLOB, rcode)    ! create,..
c                                                         .. define mode
        else                      ! open to append to existing file
          niterdb = NCOPN (iterdbfilename, NCWRITE, rcode)
        end if
      else                        ! READ from existing file
        niterdb = NCOPN (iterdbfilename, NCNOWRIT, rcode)
        if (rcode .ne. 0) go to 3 ! no error, continue
    2   write  (ncrt, 4)  TRIM(iterdbfilename)
        write  (nout, 4)  TRIM(iterdbfilename)
    4   format (                                                      /
     .      ' ERROR: subroutine ITER_DBASE1 has encountered an error' /
     .      '        the ITER database file "', a, '" cannot be opened')
        call STOP ('subroutine ITER_DBASE1: cannot open database', 220)
    3   continue
      end if
c
c     write a header line each time routine is called:
c     header lines are identified by    **            in first 2 columns
c     other comment lines (if any) have *b  (b=blank) in first 2 columns
c
      if (irwflag .eq. 0) then
c
c --- check that psir is monotonic; it may not be in certain cases where
c --- the current profile was evolved
c
        call check_monotonic (psir, nj, monotonic, 1)
        if (.not. monotonic) then    ! psir is not monotonic
          write  (niterdb, 5) time
    5     format (' the psi grid, calculated from the poloidal' /
     .            ' B field evolution, is not monotonic'        /
     .            ' therefore the data at this time (', 1pe12.6,
     .            ') was not calculated')
          go to 100
        end if
      else
        read (niterdb, 7) headerline
      end if
c
c --- some scalar quantities:
c --- integer and character parameters
c
      if (irwflag .eq. 0) then
        if (icall .eq. 0) then          ! first time, so define dim, var
          idim_1    = NCDDEF (niterdb, 'dim_scalar', 1, rcode)     ! dim
          idim_time = NCDDEF (niterdb, 'dim_time', nbctim, rcode)
c
          title = 'ONETWO netCDF output file iterdb.nc    '
          call NCAPTC (niterdb, NC_GLOBAL, 'title', NCCHAR, 39,
     .                 title, rcode)
c
          id_ishot = NCVDEF (niterdb, 'shot', NCLONG,
     .                       1, idim_1, rcode)
          call NCAPTC (niterdb, id_ishot, 'long_name', NCCHAR, 11,
     .                'shot number', rcode)
c
          id_nj = NCVDEF (niterdb, 'nj', NCLONG, 1, idim_1, rcode)
          call NCAPTC (niterdb, id_nj, 'long_name', NCCHAR, 44,
     .                'the size of the vectors printed in this file',
     .                 rcode)
          idim_nj = NCDDEF (niterdb, 'dim_nj', nj, rcode)
c
          id_nion = NCVDEF (niterdb, 'nion', NCLONG, 1, idim_1, rcode)
          call NCAPTC (niterdb, id_nion,  'long_name', NCCHAR, 25,
     .                'the number of ion species', rcode)
          idim_ion = NCDDEF (niterdb, 'dim_ion', nion, rcode)
c
          id_nprim = NCVDEF (niterdb, 'nprim', NCLONG, 1, idim_1,rcode)
          call NCAPTC (niterdb, id_nprim, 'long_name', NCCHAR, 33,
     .                'the number of primary ion species', rcode)
          idim_prim= NCDDEF (niterdb, 'dim_prim', nprim, rcode)
c
          id_nimp = NCVDEF (niterdb, 'nimp', NCLONG, 1, idim_1, rcode)
          call NCAPTC (niterdb, id_nimp, 'long_name', NCCHAR, 34,
     .                'the number of impurity ion species', rcode)
          idim_imp = NCDDEF (niterdb, 'dim_imp', nimp, rcode)
c
          id_nneu = NCVDEF (niterdb, 'nneu', NCLONG, 1, idim_1, rcode)
          call NCAPTC (niterdb, id_nneu, 'long_name', NCCHAR, 33,
     .                'the number of neutral ion species', rcode)
          idim_neu = NCDDEF (niterdb, 'dim_neu', nneu, rcode)
c
          id_ibion = NCVDEF (niterdb, 'ibion', NCLONG, 1, idim_1,rcode)
          call NCAPTC (niterdb, id_ibion, 'long_name', NCCHAR, 21,
     .                'index of beam species', rcode)
c
          id_namep = NCVDEF (niterdb, 'namep', NCCHAR,1,idim_prim,rcode)
          call NCAPTC (niterdb, id_namep, 'long_name', NCCHAR, 30,
     .                'name(s) of primary ion species', rcode)
c
          id_namei = NCVDEF (niterdb, 'namei', NCCHAR,1,idim_imp, rcode)
          call NCAPTC (niterdb, id_namei, 'long_name', NCCHAR, 31,
     .                'name(s) of impurity ion species', rcode)
c
          id_namen = NCVDEF (niterdb, 'namen', NCCHAR,1,idim_neu, rcode)
          call NCAPTC (niterdb, id_namen, 'long_name', NCCHAR, 30,
     .                'name(s) of neutral ion species', rcode)
c
          call NCENDF (niterdb, rcode)  ! leave define mode
c
c put values of variables
c
          call NCVPT1 (niterdb, id_ishot, 1, ishot, rcode)
          call NCVPT1 (niterdb, id_nj   , 1, nj   , rcode)
          call NCVPT1 (niterdb, id_nion , 1, nion , rcode)
          call NCVPT1 (niterdb, id_nprim, 1, nprim, rcode)
          call NCVPT1 (niterdb, id_nimp , 1, nimp , rcode)
          call NCVPT1 (niterdb, id_nneu , 1, nneu , rcode)
          call NCVPT1 (niterdb, id_ibion, 1, ibion, rcode)
          call NCVPTC (niterdb, id_namep, 1, nprim, namep,
     .                 8*nprim, rcode)
          call NCVPTC (niterdb, id_namei, 1, nimp , namei,
     .                 8*nimp , rcode)
          call NCVPTC (niterdb, id_namen, 1, nneu , namen,
     .                 8*nneu , rcode)
        end if
c
       else
          read (niterdb, 7)  starflag
          read (niterdb, 9)  ishot
          read (niterdb, 7)  starflag
          read (niterdb, 9)  nj
          read (niterdb, 7)  starflag
          read (niterdb, 9)  nion
          read (niterdb, 7)  starflag
          read (niterdb, 9)  nprim
          read (niterdb, 7)  starflag
          read (niterdb, 9)  nimp
          read (niterdb, 7)  starflag
          read (niterdb, 9)  nneu
          read (niterdb, 7)  starflag
          read (niterdb, 9)  ibion
          read (niterdb, 7)  starflag
          read (niterdb, 8) (namep(i), i=1,nprim)
          read (niterdb, 7)  starflag
          read (niterdb, 8) (namei(i), i=1,nimp)
          read (niterdb, 7)  starflag
          read (niterdb, 8) (namen(i), i=1,nneu)
       end if
c
c --- real (i.e., floating point) parameters
c
      rgeom    = rmajavnpsi(1)
      rmag     = rmajavnpsi(npsi)
      deltao   = triangnpsi(1)
      areao    = cxareanpsi(1)
      pindento = pindentnpsi(1)
      volo     = psivolp(1)
      te0      = te(1)
      ti0      = ti(1)
c
      if (irwflag .eq. 0) then
c
        if (icall .eq. 0) then          ! first time, so define dim, var
         call NCREDF (niterdb, rcode)   ! enter define mode
         id_time = NCVDEF (niterdb, 'time', NCDOUBLE, 1,
     .                     idim_time, rcode)
         call NCAPTC (niterdb, id_time , 'long_name', NCCHAR, 29,
     .               'time at which data is printed', rcode)
c
         id_rgeom = NCVDEF (niterdb, 'rgeom', NCDOUBLE, 1,
     .                      idim_time, rcode)
         call NCAPTC (niterdb, id_rgeom, 'long_name', NCCHAR, 62,
     . 'major radius of geometric center at elevation of magnetic axis',
     .                rcode)
         call NCAPTC (niterdb, id_rgeom, 'units', NCCHAR, 6,
     .               'meters', rcode)
c
         id_rmag = NCVDEF (niterdb, 'rmag' , NCDOUBLE, 1,
     .                     idim_time, rcode)
         call NCAPTC (niterdb, id_rmag , 'long_name', NCCHAR, 29,
     .               'major radius of magnetic axis', rcode)
         call NCAPTC (niterdb, id_rmag , 'units', NCCHAR, 6,
     .               'meters', rcode)
c
         id_rmajor = NCVDEF (niterdb, 'rmajor', NCDOUBLE, 1,
     .                       idim_time, rcode)
         call NCAPTC (niterdb, id_rmajor, 'long_name', NCCHAR, 40,
     .               'major radius of vacuum btor ref location', rcode)
         call NCAPTC (niterdb, id_rmajor, 'units', NCCHAR, 6,
     .               'meters', rcode)
c
         id_kappa = NCVDEF (niterdb, 'kappa', NCDOUBLE, 1,
     .                      idim_time, rcode)
         call NCAPTC (niterdb, id_kappa, 'long_name', NCCHAR, 17,
     .               'plasma elongation', rcode)
c
         id_deltao = NCVDEF (niterdb, 'deltao', NCDOUBLE, 1,
     .                       idim_time, rcode)
         call NCAPTC (niterdb, id_deltao, 'long_name', NCCHAR, 20,
     .               'plasma triangularity', rcode)
c
         id_pindento = NCVDEF (niterdb, 'pindento', NCDOUBLE, 1,
     .                         idim_time, rcode)
         call NCAPTC (niterdb, id_pindento, 'long_name', NCCHAR, 18,
     .               'plasma indentation', rcode)
c
         id_volo = NCVDEF (niterdb, 'volo', NCDOUBLE, 1,
     .                     idim_time, rcode)
         call NCAPTC (niterdb, id_volo, 'long_name', NCCHAR, 13,
     .               'plasma volume', rcode)
         call NCAPTC (niterdb, id_volo, 'units', NCCHAR, 8,
     .               'meters^3', rcode)
c
         id_areao = NCVDEF (niterdb, 'areao', NCDOUBLE, 1,
     .                      idim_time, rcode)
         call NCAPTC (niterdb, id_areao, 'long_name', NCCHAR, 27,
     .               'plasma cross-sectional area', rcode)
         call NCAPTC (niterdb, id_areao, 'units', NCCHAR, 8,
     .               'meters^2', rcode)
c
         id_btor = NCVDEF (niterdb, 'btor', NCDOUBLE, 1,
     .                     idim_time, rcode)
         call NCAPTC (niterdb, id_btor, 'long_name', NCCHAR, 31,
     .               'vacuum toroidal field at rmajor', rcode)
         call NCAPTC (niterdb, id_btor, 'units', NCCHAR, 5,
     .               'tesla', rcode)
c
         id_totcur = NCVDEF (niterdb, 'totcur', NCDOUBLE, 1,
     .                       idim_time, rcode)
         call NCAPTC (niterdb, id_totcur, 'long_name', NCCHAR, 14,
     .               'total currents', rcode)
         call NCAPTC (niterdb, id_totcur, 'units', NCCHAR, 4,
     .               'amps', rcode)
c
         id_totohm = NCVDEF (niterdb, 'totohm', NCDOUBLE, 1,
     .                       idim_time, rcode)
         call NCAPTC (niterdb, id_totohm, 'long_name', NCCHAR, 14,
     .               'ohmic currents', rcode)
         call NCAPTC (niterdb, id_totohm, 'units', NCCHAR, 4,
     .               'amps', rcode)
c
         id_totboot = NCVDEF (niterdb, 'totboot', NCDOUBLE, 1,
     .                        idim_time, rcode)
         call NCAPTC (niterdb, id_totboot, 'long_name', NCCHAR, 18,
     .               'bootstrap currents', rcode)
         call NCAPTC (niterdb, id_totboot, 'units', NCCHAR, 4,
     .               'amps', rcode)
c
         id_totbeam = NCVDEF (niterdb, 'totbeam', NCDOUBLE, 1,
     .                        idim_time, rcode)
         call NCAPTC (niterdb, id_totbeam, 'long_name', NCCHAR, 13,
     .               'beam currents', rcode)
         call NCAPTC (niterdb, id_totbeam, 'units', NCCHAR, 4,
     .               'amps', rcode)
c
         id_totrf = NCVDEF (niterdb, 'totrf', NCDOUBLE, 1,
     .                      idim_time, rcode)
         call NCAPTC (niterdb, id_totrf, 'long_name', NCCHAR, 11,
     .               'RF currents', rcode)
         call NCAPTC (niterdb, id_totrf, 'units', NCCHAR, 4,
     .               'amps', rcode)
c
         id_betap = NCVDEF (niterdb, 'betap', NCDOUBLE, 1,
     .                      idim_time, rcode)
         call NCAPTC (niterdb, id_betap, 'long_name', NCCHAR, 13,
     .               'poloidal beta', rcode)
c
         id_beta = NCVDEF (niterdb, 'beta', NCDOUBLE, 1,
     .                     idim_time, rcode)
         call NCAPTC (niterdb, id_beta, 'long_name', NCCHAR, 13,
     .               'toroidal beta', rcode)
c
         id_ali = NCVDEF (niterdb, 'ali', NCDOUBLE, 1,
     .                    idim_time, rcode)
         call NCAPTC (niterdb, id_ali, 'long_name', NCCHAR, 17,
     .               'plasma inductance', rcode)
c
         id_te0 = NCVDEF (niterdb, 'te0', NCDOUBLE, 1,
     .                    idim_time, rcode)
         call NCAPTC (niterdb, id_te0, 'long_name', NCCHAR, 28,
     .               'central electron temperature', rcode)
c
         id_ti0 = NCVDEF (niterdb, 'ti0', NCDOUBLE, 1,
     .                    idim_time, rcode)
         call NCAPTC (niterdb, id_ti0, 'long_name', NCCHAR, 23,
     .               'central ion temperature', rcode)
c
         call NCENDF (niterdb, rcode)  ! leave define mode
       end if
c
c      put values of variables defined when icall = 0
c
         call NCVPT1 (niterdb, id_time    , icall+1, time    , rcode)
         call NCVPT1 (niterdb, id_rgeom   , icall+1,
     .                rgeom   *0.01                          , rcode)
         call NCVPT1 (niterdb, id_rmag    , icall+1,
     .                rmag    *0.01                          , rcode)
         call NCVPT1 (niterdb, id_rmajor  , icall+1,
     .                rmajor  *0.01                          , rcode)
         call NCVPT1 (niterdb, id_kappa   , icall+1, kappa   , rcode)
         call NCVPT1 (niterdb, id_deltao  , icall+1, deltao  , rcode)
         call NCVPT1 (niterdb, id_pindento, icall+1, pindento, rcode)
         call NCVPT1 (niterdb, id_volo    , icall+1,
     .                volo    *1.0e-6                        , rcode)
         call NCVPT1 (niterdb, id_areao   , icall+1,
     .                areao   *1.0e-4                        , rcode)
         call NCVPT1 (niterdb, id_btor    , icall+1,
     .                btor    *1.0e-4                        , rcode)
         call NCVPT1 (niterdb, id_totcur  , icall+1, totcur  , rcode)
         call NCVPT1 (niterdb, id_totohm  , icall+1, totohm  , rcode)
         call NCVPT1 (niterdb, id_totboot , icall+1, totboot , rcode)
         call NCVPT1 (niterdb, id_totbeam , icall+1, totbeam , rcode)
         call NCVPT1 (niterdb, id_totrf   , icall+1, totrf   , rcode)
         call NCVPT1 (niterdb, id_betap   , icall+1, betap   , rcode)
         call NCVPT1 (niterdb, id_beta    , icall+1, beta    , rcode)
         call NCVPT1 (niterdb, id_ali     , icall+1, ali     , rcode)
         call NCVPT1 (niterdb, id_te0     , icall+1, te0     , rcode)
         call NCVPT1 (niterdb, id_ti0     , icall+1, ti0     , rcode)
      else
         read (niterdb,  7) starflag
         read (niterdb, 10) time
         read (niterdb,  7) starflag
         read (niterdb, 10) rgeom
         read (niterdb,  7) starflag
         read (niterdb, 10) rmag
         read (niterdb,  7) starflag
         read (niterdb, 10) rmajor
         read (niterdb,  7) starflag
         read (niterdb, 10) kappa
         read (niterdb,  7) starflag
         read (niterdb, 10) deltao
         read (niterdb,  7) starflag
         read (niterdb, 10) pindento
         read (niterdb,  7) starflag
         read (niterdb, 10) volo
         read (niterdb,  7) starflag
         read (niterdb, 10) areao
         read (niterdb,  7) starflag
         read (niterdb, 10) btor
         read (niterdb,  7) starflag
         read (niterdb, 10) totcur, totohm, totboot, totbeam, totrf
         read (niterdb,  7) starflag
         read (niterdb, 10) betap
         read (niterdb,  7) starflag
         read (niterdb, 10) beta
         read (niterdb,  7) starflag
         read (niterdb, 10) ali
         read (niterdb,  7) starflag
         read (niterdb, 10) te0
         read (niterdb,  7) starflag
         read (niterdb, 10) ti0
      end if
c
      if (irwflag .eq. 0) then
        if (icall .eq. 0) then          ! first time, so define dim, var
          call NCREDF (niterdb, rcode)  ! enter define mode
c
c --- psir grid
c
          dim_njtime(1) = idim_nj
          dim_njtime(2) = idim_time
          id_psir = NCVDEF (niterdb, 'psir', NCDOUBLE, 2,
     .                      dim_njtime, rcode)
          call NCAPTC (niterdb, id_psir, 'long_name', NCCHAR, 15,
     .                'psi on rho grid', rcode)
          call NCAPTC (niterdb, id_psir, 'units', NCCHAR, 18,
     .                'volt*second/radian', rcode)
c
c --- rho grid
c
          id_r = NCVDEF (niterdb, 'r', NCDOUBLE, 2,
     .                   dim_njtime, rcode)
          call NCAPTC (niterdb, id_r, 'long_name', NCCHAR, 8,
     .                'rho grid', rcode)
          call NCAPTC (niterdb, id_r, 'units', NCCHAR, 6,
     .                'meters', rcode)
c
c --- fcap
c
          id_fcap = NCVDEF (niterdb, 'fcap', NCDOUBLE, 2,
     .                      dim_njtime, rcode)
          call NCAPTC (niterdb, id_fcap, 'long_name', NCCHAR, 30,
     .                'fcap, (i.e., f(psilim)/f(psi))', rcode)
c
c --- gcap
c
          id_gcap = NCVDEF (niterdb, 'gcap', NCDOUBLE, 2,
     .                      dim_njtime, rcode)
          call NCAPTC (niterdb, id_gcap, 'long_name', NCCHAR, 39,
     .                'gcap, (i.e., <(grad rho)**2*(R0/R)**2>)', rcode)
c
c --- hcap
c
          id_hcap = NCVDEF (niterdb, 'hcap', NCDOUBLE, 2,
     .                      dim_njtime, rcode)
          call NCAPTC (niterdb, id_hcap, 'long_name', NCCHAR, 39,
     .                'hcap, (i.e., (dvolume/drho)/(4*pi*pi*R0*rho))',
     .                 rcode)
c
c --- te (in keV)
c
          id_te = NCVDEF (niterdb, 'te', NCDOUBLE, 2,
     .                    dim_njtime, rcode)
          call NCAPTC (niterdb, id_te, 'long_name', NCCHAR, 20,
     .                'electron temperature', rcode)
          call NCAPTC (niterdb, id_te, 'units', NCCHAR, 3,
     .                'keV', rcode)
c
c --- ti (in keV)
c
          id_ti = NCVDEF (niterdb, 'ti', NCDOUBLE, 2,
     .                    dim_njtime, rcode)
          call NCAPTC (niterdb, id_ti, 'long_name', NCCHAR, 15,
     .                'ion temperature', rcode)
          call NCAPTC (niterdb, id_ti, 'units', NCCHAR, 3,
     .                'keV', rcode)
c
c --- safety factor
c
          id_q = NCVDEF (niterdb, 'q', NCDOUBLE, 2,
     .                   dim_njtime, rcode)
          call NCAPTC (niterdb, id_q, 'long_name', NCCHAR, 21,
     .                'safety factor profile', rcode)
c
c --- electron density
c
          id_ene = NCVDEF (niterdb, 'ene', NCDOUBLE, 2,
     .                     dim_njtime, rcode)
          call NCAPTC (niterdb, id_ene, 'long_name', NCCHAR, 16,
     .                'electron density', rcode)
          call NCAPTC (niterdb, id_ene, 'units', NCCHAR, 9,
     .                '#/meter^3', rcode)
c
          call NCENDF (niterdb, rcode)  ! leave define mode
        end if
c
          c1t(1) = 1
          c1t(2) = icall + 1
          cnt(1) = nj
          cnt(2) = icall + 1
          call multpl2 (psir, tmpdata, nj, psimks)
          call NCVPT (niterdb, id_psir, c1t, cnt, tmpdata, rcode)
          const  = 1.0e-02
          call multpl2 (r, tmpdata, nj, const)
          call NCVPT (niterdb, id_r   , c1t, cnt, tmpdata, rcode)
          call NCVPT (niterdb, id_fcap, c1t, cnt, fcap   , rcode)
          call NCVPT (niterdb, id_gcap, c1t, cnt, gcap   , rcode)
          call NCVPT (niterdb, id_hcap, c1t, cnt, hcap   , rcode)
          call NCVPT (niterdb, id_te  , c1t, cnt, te     , rcode)
          call NCVPT (niterdb, id_ti  , c1t, cnt, ti     , rcode)
          do i=1,nj
             tmpdata(i) = ABS (q(i))
          end do
          call NCVPT (niterdb, id_q   , c1t, cnt, tmpdata, rcode)
          const = 1.0e+6
          call multpl2 (ene, tmpdata, nj, const)
          call NCVPT (niterdb, id_ene , c1t, cnt, tmpdata, rcode)
      else
          read (niterdb,  7)  starflag
          read (niterdb, 10) (psir(j), j=1,nj)
          read (niterdb,  7)  starflag
          read (niterdb, 10) (r(j)   , j=1,nj)
          read (niterdb,  7)  starflag
          read (niterdb, 10) (fcap(j), j=1,nj)
          read (niterdb,  7)  starflag
          read (niterdb, 10) (gcap(j), j=1,nj)
          read (niterdb,  7)  starflag
          read (niterdb, 10) (hcap(j), j=1,nj)
          read (niterdb,  7)  starflag
          read (niterdb, 10) (te(j), j=1,nj)
          read (niterdb,  7)  starflag
          read (niterdb, 10) (ti(j), j=1,nj)
          read (niterdb,  7)  starflag
          read (niterdb, 10) (q(j) , j=1,nj)
          read (niterdb,  7)  starflag
          read (niterdb, 10) (ene(j), j=1,nj)
      end if
c
c --- primary, impurity densities
c
      if (irwflag .eq. 0) then
        if (icall .eq. 0) then          ! first time, so define dim, var
          call NCREDF (niterdb, rcode)  ! enter define mode
c
c --- primary
c
          dim_njprimtime(1) = idim_nj
          dim_njprimtime(2) = idim_prim
          dim_njprimtime(3) = idim_time
          id_enp = NCVDEF (niterdb, 'enp', NCDOUBLE, 3,
     .                     dim_njprimtime, rcode)
          call NCAPTC (niterdb, id_enp, 'long_name', NCCHAR, 19,
     .                'primary ion density', rcode)
          call NCAPTC (niterdb, id_enp, 'units', NCCHAR, 9,
     .                '#/meter^3', rcode)
          call NCAPT  (niterdb, id_enp, 'species', NCCHAR,
     .                 nprim, namep, rcode)
c
c --- impurity
c
          dim_njimptime(1) = idim_nj
          dim_njimptime(2) = idim_imp
          dim_njimptime(3) = idim_time
          id_eni = NCVDEF (niterdb, 'eni', NCDOUBLE, 3,
     .                     dim_njimptime, rcode)
          call NCAPTC (niterdb, id_eni, 'long_name', NCCHAR, 20,
     .                'impurity ion density', rcode)
          call NCAPTC (niterdb, id_eni, 'units', NCCHAR, 9,
     .                '#/meter^3', rcode)
          call NCAPT  (niterdb, id_eni, 'species', NCCHAR,
     .                 nimp, namei, rcode)
c
          call NCENDF (niterdb, rcode)  ! leave define mode
        end if
c
c --- put values in the variables defined when icall .eq. 0
c
          c11t(1) = 1
          c11t(2) = 1
          c11t(3) = icall + 1
          cnpt(1) = nj
          cnpt(2) = nprim
          cnpt(3) = icall + 1
          const   = 1.0e+6
          call multpl2 (en, tmpdata, nj, const)
          call NCVPT (niterdb, id_enp, c11t, cnpt,tmpdata, rcode)
          cnit(1) = nj
          cnit(2) = nimp
          cnit(3) = icall + 1
          do i=1,nj
             do j=1,nimp
                eni(i,j) = en(i,j+nprim) *1.0e+6
             end do
          end do
          call NCVPT (niterdb, id_eni, c11t, cnit, eni, rcode)
      else
          do jj=1,nion
               read (niterdb,  7)  starflag
               read (niterdb, 10) (en(j,jj)       , j=1,nj)
          end do
      end if
c
c --- primary ion species sources
c
      if (irwflag .eq. 0) then
        if (icall .eq. 0) then          ! first time, so define dim, var
          call NCREDF (niterdb, rcode)  ! enter define mode
c
          id_sion = NCVDEF (niterdb, 'sion', NCDOUBLE, 3,
     .                      dim_njprimtime, rcode)
          call NCAPTC (niterdb, id_sion, 'long_name', NCCHAR, 24,
     .                'source due to ionization', rcode)
          call NCAPTC (niterdb, id_sion, 'units', NCCHAR, 18,
     .                '#/(meter^3*second)', rcode)
          call NCAPT  (niterdb, id_sion, 'species', NCCHAR,
     .                 nprim, namep, rcode)
c
          id_srecom = NCVDEF (niterdb, 'srecom', NCDOUBLE, 3,
     .                        dim_njprimtime, rcode)
          call NCAPTC (niterdb, id_srecom, 'long_name', NCCHAR, 27,
     .                'source due to recombination', rcode)
          call NCAPTC (niterdb, id_srecom, 'units', NCCHAR, 18,
     .                '#/(meter^3*second)', rcode)
          call NCAPT  (niterdb, id_srecom, 'species', NCCHAR,
     .                 nprim, namep, rcode)
c
          id_scx = NCVDEF (niterdb, 'scx', NCDOUBLE, 3,
     .                     dim_njprimtime, rcode)
          call NCAPTC (niterdb, id_scx, 'long_name', NCCHAR, 30,
     .                'source due to cx thermal neut.', rcode)
          call NCAPTC (niterdb, id_scx, 'units', NCCHAR, 18,
     .                '#/(meter^3*second)', rcode)
          call NCAPT  (niterdb, id_scx, 'species', NCCHAR,
     .                 nprim, namep, rcode)
c
          id_sbcx = NCVDEF (niterdb, 'sbcx', NCDOUBLE, 3,
     .                      dim_njprimtime, rcode)
          call NCAPTC (niterdb, id_sbcx, 'long_name', NCCHAR, 30,
     .                'sink due to cx with beam neut.', rcode)
          call NCAPTC (niterdb, id_sbcx, 'units', NCCHAR, 18,
     .                '#/(meter^3*second)', rcode)
          call NCAPT  (niterdb, id_sbcx, 'species', NCCHAR,
     .                 nprim, namep, rcode)
c
          id_s = NCVDEF (niterdb, 's', NCDOUBLE, 3,
     .                   dim_njprimtime, rcode)
          call NCAPTC (niterdb, id_s, 'long_name', NCCHAR, 17,
     .                'total source rate', rcode)
          call NCAPTC (niterdb, id_s, 'units', NCCHAR, 18,
     .                '#/(meter^3*second)', rcode)
          call NCAPT  (niterdb, id_s, 'species', NCCHAR,
     .                 nprim, namep, rcode)
c
          id_dudtsv = NCVDEF (niterdb, 'dudtsv', NCDOUBLE, 3,
     .                        dim_njprimtime, rcode)
          call NCAPTC (niterdb, id_dudtsv, 'long_name', NCCHAR, 5,
     .                's dot', rcode)
          call NCAPTC (niterdb, id_dudtsv, 'units', NCCHAR, 18,
     .                '#/(meter^3*second)', rcode)
          call NCAPT  (niterdb, id_dudtsv, 'species', NCCHAR,
     .                 nprim, namep, rcode)
          call NCENDF (niterdb, rcode)  ! leave define mode
        end if
c
          const = 1.0e+6
          call multpl2 (sion, tmpdata, nj, const)
          call NCVPT (niterdb, id_sion  , c11t, cnpt, tmpdata, rcode)
          call multpl2 (srecom, tmpdata, nj, const)
          call NCVPT (niterdb, id_srecom, c11t, cnpt, tmpdata, rcode)
          call multpl2 (scx, tmpdata, nj, const)
          call NCVPT (niterdb, id_scx   , c11t, cnpt, tmpdata, rcode)
          call multpl2 (sbcx, tmpdata, nj, const)
          call NCVPT (niterdb, id_sbcx  , c11t, cnpt, tmpdata, rcode)
c
c ---     dimension of s and dudtsv (kk,kj) are the reverse of their
c ---     netCDF variables (nj,nprim); set the mapping vector
c
          cmap(3) = 0
          cmap(2) = 8                   ! size of real*8
          cmap(1) = cmap(2) * kk        ! kk - first dimension
          const   = 1.0e+6
          call multpl2 (s, tmpdata, nj, const)
          call NCVPTG (niterdb, id_s     , c11t, cnpt, 0, cmap,
     .                 tmpdata                       , rcode)
          call multpl2 (dudtsv, tmpdata, nj, const)
          call NCVPTG (niterdb, id_dudtsv, c11t, cnpt, 0, cmap,
     .                 tmpdata                       , rcode)
c
      else
        do jj=1,nprim
          read (niterdb,  7)  starflag
          read (niterdb, 10) (sion(j,jj)  , j=1,nj)
          read (niterdb,  7)  starflag
          read (niterdb, 10) (srecom(j,jj), j=1,nj)
          read (niterdb,  7)  starflag
          read (niterdb, 10) (scx(j,jj)   , j=1,nj)
          read (niterdb,  7)  starflag
          read (niterdb, 10) (sbcx(j,jj)  , j=1,nj)
          read (niterdb,  7)  starflag
          read (niterdb, 10) (s(jj,j)     , j=1,nj)
          read (niterdb,  7)  starflag
          read (niterdb, 10) (dudtsv(jj,j), j=1,nj)
        end do
      end if
c
       if (irwflag .eq. 0) then
        if (icall .eq. 0) then          ! first time, so define dim, var
          call NCREDF (niterdb, rcode)  ! enter define mode
c
c --- fast ion density
c
          id_enbeam = NCVDEF (niterdb, 'enbeam', NCDOUBLE, 2,
     .                        dim_njtime, rcode)
          call NCAPTC (niterdb, id_enbeam, 'long_name', NCCHAR, 16,
     .                'fast ion density', rcode)
          call NCAPTC (niterdb, id_enbeam, 'units', NCCHAR, 9,
     .                '#/meter^3', rcode)
          call NCAPT  (niterdb, id_enbeam, 'species', NCCHAR,
     .                 1, namep(ibion), rcode)
c
c --- neutral densities
c
          dim_njneutime(1) = idim_nj
          dim_njneutime(2) = idim_neu
          dim_njneutime(3) = idim_time
          id_enn = NCVDEF (niterdb, 'enn', NCDOUBLE, 3,
     .                     dim_njneutime, rcode)
          call NCAPTC (niterdb, id_enn, 'long_name', NCCHAR, 15,
     .                'neutral density', rcode)
          call NCAPTC (niterdb, id_enn, 'units', NCCHAR, 9,
     .                '#/meter^3', rcode)
          call NCAPT  (niterdb, id_enn, 'species', NCCHAR,
     .                 nneu, namen, rcode)
c
c --- neutral densities
c
          id_ennw = NCVDEF (niterdb, 'ennw', NCDOUBLE, 3,
     .                      dim_njneutime, rcode)
          call NCAPTC (niterdb, id_ennw, 'long_name', NCCHAR, 32,
     .                'neutral density from wall source', rcode)
          call NCAPTC (niterdb, id_ennw, 'units', NCCHAR, 9,
     .                '#/meter^3', rcode)
          call NCAPT  (niterdb, id_ennw, 'species', NCCHAR,
     .                 nneu, namen, rcode)
c
c --- neutral densities
c
          id_ennv = NCVDEF (niterdb, 'ennv', NCDOUBLE, 3,
     .                      dim_njneutime, rcode)
          call NCAPTC (niterdb, id_ennv, 'long_name', NCCHAR, 34,
     .                'neutral density from volume source', rcode)
          call NCAPTC (niterdb, id_ennv, 'units', NCCHAR, 9,
     .                '#/meter^3', rcode)
          call NCAPT  (niterdb, id_ennv, 'species', NCCHAR,
     .                 nneu, namen, rcode)
c
c --- neutral source
c
          id_volsn = NCVDEF (niterdb, 'volsn', NCDOUBLE, 3,
     .                       dim_njneutime, rcode)
          call NCAPTC (niterdb, id_volsn, 'long_name', NCCHAR, 25,
     .                'volume source of neutrals', rcode)
          call NCAPTC (niterdb, id_volsn, 'units', NCCHAR, 18,
     .                '#/(meter^3*second)', rcode)
          call NCAPT  (niterdb, id_volsn, 'species', NCCHAR,
     .                 nneu, namen, rcode)
c
c --- electron source due to beams
c
          id_sbion = NCVDEF (niterdb, 'sbion', NCDOUBLE, 2,
     .                       dim_njtime, rcode)
          call NCAPTC (niterdb, id_sbion, 'long_name', NCCHAR, 20,
     .                'beam electron source', rcode)
          call NCAPTC (niterdb, id_sbion, 'units', NCCHAR, 18,
     .                '#/(meter^3*second)', rcode)
c
c --- thermal ion source due to beams
c
          id_sbeam = NCVDEF (niterdb, 'sbeam', NCDOUBLE, 2,
     .                       dim_njtime, rcode)
          call NCAPTC (niterdb, id_sbeam, 'long_name', NCCHAR,
     .                 23, 'beam thermal ion source', rcode)
          call NCAPTC (niterdb, id_sbeam, 'units', NCCHAR, 18,
     .                '#/(meter^3*second)', rcode)
c
c --- current density
c
          id_curden = NCVDEF (niterdb, 'curden', NCDOUBLE, 2,
     .                        dim_njtime, rcode)
          call NCAPTC (niterdb, id_curden, 'long_name', NCCHAR, 21,
     .                'total current density', rcode)
          call NCAPTC (niterdb, id_curden, 'units', NCCHAR, 12,
     .                'amps/meter^2', rcode)
c
c --- ohmic current density
c
          id_curohm = NCVDEF (niterdb, 'curohm', NCDOUBLE, 2,
     .                        dim_njtime, rcode)
          call NCAPTC (niterdb, id_curohm, 'long_name', NCCHAR, 21,
     .                'ohmic current density', rcode)
          call NCAPTC (niterdb, id_curohm, 'units', NCCHAR, 12,
     .                'amps/meter^2', rcode)
c
c --- bootstrap current density
c
          id_curboot = NCVDEF (niterdb, 'curboot', NCDOUBLE, 2,
     .                         dim_njtime, rcode)
          call NCAPTC (niterdb, id_curboot, 'long_name', NCCHAR, 25,
     .                'bootstrap current density', rcode)
          call NCAPTC (niterdb, id_curboot, 'units', NCCHAR, 12,
     .                'amps/meter^2', rcode)
c
c --- beam current density
c
          id_curdbeam = NCVDEF (niterdb, 'curdbeam', NCDOUBLE, 2,
     .                          dim_njtime, rcode)
          call NCAPTC (niterdb, id_curdbeam, 'long_name', NCCHAR, 27,
     .                'beam-driven current density', rcode)
          call NCAPTC (niterdb, id_curdbeam, 'units', NCCHAR, 12,
     .                'amps/meter^2', rcode)
c
c --- RF current density
c
          id_currf = NCVDEF (niterdb, 'currf', NCDOUBLE, 2,
     .                       dim_njtime, rcode)
          call NCAPTC (niterdb, id_currf, 'long_name', NCCHAR, 18,
     .                'RF current density', rcode)
          call NCAPTC (niterdb, id_currf, 'units', NCCHAR, 12,
     .                'amps/meter^2', rcode)
c
c --- rho*bp0*fcap*gcap*hcap
c
          id_rbp = NCVDEF (niterdb, 'rbp', NCDOUBLE, 2,
     .                     dim_njtime, rcode)
          call NCAPTC (niterdb, id_rbp, 'long_name', NCCHAR, 22,
     .                'rho*bp0*fcap*gcap*hcap', rcode)
          call NCAPTC (niterdb, id_rbp, 'units', NCCHAR, 12,
     .                'tesla*meters', rcode)
c
c --- zeff profile
c
          id_zeff = NCVDEF (niterdb, 'zeff', NCDOUBLE, 2,
     .                      dim_njtime, rcode)
          call NCAPTC (niterdb, id_zeff, 'long_name', NCCHAR, 12,
     .                'zeff profile', rcode)
c
c --- angular rotation speed profile
c
          id_angrot = NCVDEF (niterdb, 'angrot', NCDOUBLE, 2,
     .                        dim_njtime, rcode)
          call NCAPTC (niterdb, id_angrot, 'long_name', NCCHAR, 30,
     .                'angular rotation speed profile', rcode)
          call NCAPTC (niterdb, id_angrot, 'units', NCCHAR, 7,
     .                'rad/sec', rcode)
c
c --- fast ion torque density
c
          id_storqueb = NCVDEF (niterdb, 'storqueb', NCDOUBLE, 2,
     .                        dim_njtime, rcode)
          call NCAPTC (niterdb, id_storqueb, 'long_name', NCCHAR, 19,
     .                'beam torque density', rcode)
          call NCAPTC (niterdb, id_storqueb, 'units', NCCHAR, 7,
     .                'nt-m/m**3', rcode)
c
c --- cer impurity ion poloidal velocity/Bpol
c
          id_Kpol_c = NCVDEF (niterdb, 'Kpol_c', NCDOUBLE,
     .            2, dim_njtime, rcode)
          call NCAPTC (niterdb, id_Kpol_c, 'long_name', NCCHAR,
     .            35, 'CER impurity poloidal velocity/Bpol', rcode)
          call NCAPTC (niterdb, id_Kpol_c, 'units', NCCHAR,
     .            19, 'meters/second/Tesla', rcode)
c
c --- Experimental cer impurity ion poloidal velocity/Bpol
c
          id_Kpol_exp = NCVDEF (niterdb, 'Kpol_exp', NCDOUBLE,
     .            2, dim_njtime, rcode)
          call NCAPTC (niterdb, id_Kpol_exp, 'long_name', NCCHAR,
     .     44, 'Measured CER impurity poloidal velocity/Bpol', rcode)
          call NCAPTC (niterdb, id_Kpol_exp, 'units', NCCHAR,
     .            19, 'meters/second/Tesla', rcode)
c
c --- main ion poloidal velocity/Bpol
c
          id_Kpol_d = NCVDEF (niterdb, 'Kpol_d', NCDOUBLE,
     .            2, dim_njtime, rcode)
          call NCAPTC (niterdb, id_Kpol_d, 'long_name', NCCHAR,
     .            31, 'main ion poloidal velocity/Bpol', rcode)
          call NCAPTC (niterdb, id_Kpol_d, 'units', NCCHAR,
     .            19, 'meters/second/Tesla', rcode)
c
c --- cer impurity ion toroidal angular velocity
c
          id_angrot_c = NCVDEF (niterdb, 'angrot_c', NCDOUBLE,
     .            2, dim_njtime, rcode)
          call NCAPTC (niterdb, id_angrot_c, 'long_name', NCCHAR,
     .            38, 'CER impurity toroidal angular velocity', rcode)
          call NCAPTC (niterdb, id_angrot_c, 'units', NCCHAR,
     .            14, 'radians/second', rcode)
c
c --- main ion toroidal angular velocity
c
          id_angrot_d = NCVDEF (niterdb, 'angrot_d', NCDOUBLE,
     .            2, dim_njtime, rcode)
          call NCAPTC (niterdb, id_angrot_d, 'long_name', NCCHAR,
     .            34, 'main ion toroidal angular velocity', rcode)
          call NCAPTC (niterdb, id_angrot_d, 'units', NCCHAR,
     .            14, 'radians/second', rcode)
c
c --- main ion flux average parallel velocity
c
          id_ave_vpar_d = NCVDEF (niterdb, 'ave_vpar_d', NCDOUBLE,
     .             2, dim_njtime, rcode)
          call NCAPTC (niterdb, id_ave_vpar_d, 'long_name', NCCHAR,
     .            34, 'main ion average parallel velocity', rcode)
          call NCAPTC (niterdb, id_ave_vpar_d, 'units', NCCHAR,
     .            13, 'meters/second', rcode)
c
c --- main ion diamagnetic velocity/RBpol
c
          id_udia_d = NCVDEF (niterdb, 'udia_d', NCDOUBLE,
     .            2, dim_njtime, rcode)
          call NCAPTC (niterdb, id_udia_d, 'long_name', NCCHAR,
     .            35, 'main ion diamagnetic velocity/RBpol', rcode)
          call NCAPTC (niterdb, id_udia_d, 'units', NCCHAR,
     .            8, '1/second', rcode)
c
c --- cer ion diamagnetic velocity/RBpol
c
          id_udia_c = NCVDEF (niterdb, 'udia_c', NCDOUBLE,
     .            2, dim_njtime, rcode)
          call NCAPTC (niterdb, id_udia_c, 'long_name', NCCHAR,
     .            33, 'CER ion diamagnetic velocity/RBpol', rcode)
          call NCAPTC (niterdb, id_udia_c, 'units', NCCHAR,
     .            8, '1/second', rcode)
c
c --- ion temperature gradient/(eRBpol)
c
          id_ugrt = NCVDEF (niterdb, 'ugrt', NCDOUBLE,
     .            2, dim_njtime, rcode)
          call NCAPTC (niterdb, id_ugrt, 'long_name', NCCHAR,
     .            33, 'ion temperature gradient/(eRBpol)', rcode)
          call NCAPTC (niterdb, id_ugrt, 'units', NCCHAR,
     .            8, '1/second', rcode)
c
c --- main ion orbit squeeze factor
c
          id_sqz_d = NCVDEF (niterdb, 'sqz_d', NCDOUBLE,
     .             2, dim_njtime, rcode)
          call NCAPTC (niterdb, id_sqz_d, 'long_name', NCCHAR,
     .            31, 'main ion orbit squeezing factor', rcode)
          call NCAPTC (niterdb, id_sqz_d, 'units', NCCHAR,
     .            13, 'dimensionless', rcode)
c
c --- Er/RBp neoclassical
c
          id_Epsi = NCVDEF (niterdb, 'Epsi', NCDOUBLE,
     .            2, dim_njtime, rcode)
          call NCAPTC (niterdb, id_Epsi, 'long_name', NCCHAR,
     .            19, 'Er/RBp neoclassical', rcode)
          call NCAPTC (niterdb, id_Epsi, 'units', NCCHAR,
     .            8, '1/second', rcode)
c
c --- Er/RBp experimental
c
          id_Epsi_exp = NCVDEF (niterdb, 'Epsi_exp', NCDOUBLE,
     .            2, dim_njtime, rcode)
          call NCAPTC (niterdb, id_Epsi_exp, 'long_name', NCCHAR,
     .            19, 'Er/RBp experimental', rcode)
          call NCAPTC (niterdb, id_Epsi_exp, 'units', NCCHAR,
     .            8, '1/second', rcode)
c
c --- Bp along Z=0 cord
c
          id_cer_bp = NCVDEF (niterdb, 'cer_bp', NCDOUBLE,
     .            2, dim_njtime, rcode)
          call NCAPTC (niterdb, id_cer_bp, 'long_name', NCCHAR,
     .            17, 'Bp along Z=0 cord', rcode)
          call NCAPTC (niterdb, id_cer_bp, 'units', NCCHAR,
     .            5, 'Tesla', rcode)
c
c --- Bt/R along Z=0 cord
c
          id_cer_btdr = NCVDEF (niterdb, 'cer_btdr', NCDOUBLE,
     .            2, dim_njtime, rcode)
          call NCAPTC (niterdb, id_cer_btdr, 'long_name', NCCHAR,
     .            19, 'Bt/R along Z=0 cord', rcode)
          call NCAPTC (niterdb, id_cer_btdr, 'units', NCCHAR,
     .            11, 'Tesla/meter', rcode)
c
c --- thermal diff. profiles, electron and ion
c
          id_chieinv = NCVDEF (niterdb, 'chieinv', NCDOUBLE, 2,
     .                         dim_njtime, rcode)
          call NCAPTC (niterdb, id_chieinv, 'long_name', NCCHAR, 28,
     .                'electron thermal diffusivity', rcode)
          call NCAPTC (niterdb, id_chieinv, 'units', NCCHAR, 25,
     .                'meters^2/sec on half grid', rcode)
c
          id_chiinv = NCVDEF (niterdb, 'chiinv', NCDOUBLE, 2,
     .                        dim_njtime, rcode)
          call NCAPTC (niterdb, id_chiinv, 'long_name', NCCHAR, 23,
     .                'ion thermal diffusivity', rcode)
          call NCAPTC (niterdb, id_chiinv, 'units', NCCHAR, 25,
     .                'meters^2/sec on half grid', rcode)
c
c --- ion neoclassical thermal conductivity
c
          id_xkineo = NCVDEF (niterdb, 'xkineo', NCDOUBLE, 2,
     .                        dim_njtime, rcode)
          call NCAPTC (niterdb, id_xkineo, 'long_name', NCCHAR, 37,
     .                'ion neoclassical thermal conductivity', rcode)
          call NCAPTC (niterdb, id_xkineo, 'units', NCCHAR, 29,
     .                '1/(meter*second) on half grid', rcode)
c
c --- d(electron energy)/dt profile
c
          id_dpedtc = NCVDEF (niterdb, 'dpedtc', NCDOUBLE, 2,
     .                        dim_njtime, rcode)
          call NCAPTC (niterdb, id_dpedtc, 'long_name', NCCHAR, 15,
     .                'wdot, electrons', rcode)
          call NCAPTC (niterdb, id_dpedtc, 'units', NCCHAR, 13,
     .                'watts/meter^3', rcode)
c
c --- d(ion energy)/dt profile
c
          id_dpidtc = NCVDEF (niterdb, 'dpidtc', NCDOUBLE, 2,
     .                        dim_njtime, rcode)
          call NCAPTC (niterdb, id_dpidtc, 'long_name', NCCHAR, 10,
     .                'wdot, ions', rcode)
          call NCAPTC (niterdb, id_dpidtc, 'units', NCCHAR, 13,
     .                'watts/meter^3', rcode)
c
c --- electron conduction profile
c
          id_qconde = NCVDEF (niterdb, 'qconde', NCDOUBLE, 2,
     .                        dim_njtime, rcode)
          call NCAPTC (niterdb, id_qconde, 'long_name', NCCHAR, 19,
     .                'electron conduction', rcode)
          call NCAPTC (niterdb, id_qconde, 'units', NCCHAR, 13,
     .                'watts/meter^3', rcode)
c
c --- ion conduction profile
c
          id_qcondi = NCVDEF (niterdb, 'qcondi', NCDOUBLE, 2,
     .                        dim_njtime, rcode)
          call NCAPTC (niterdb, id_qcondi, 'long_name', NCCHAR, 14,
     .                'ion conduction', rcode)
          call NCAPTC (niterdb, id_qcondi, 'units', NCCHAR, 13,
     .                'watts/meter^3', rcode)
c
c --- electron convection profile
c
          id_qconve = NCVDEF (niterdb, 'qconve', NCDOUBLE, 2,
     .                        dim_njtime, rcode)
          call NCAPTC (niterdb, id_qconve, 'long_name', NCCHAR, 19,
     .                'electron convection', rcode)
          call NCAPTC (niterdb, id_qconve, 'units', NCCHAR, 13,
     .                'watts/meter^3', rcode)
c
c --- ion convection profile
c
          id_qconvi = NCVDEF (niterdb, 'qconvi', NCDOUBLE, 2,
     .                        dim_njtime, rcode)
          call NCAPTC (niterdb, id_qconvi, 'long_name', NCCHAR, 14,
     .                'ion convection', rcode)
          call NCAPTC (niterdb, id_qconvi, 'units', NCCHAR, 13,
     .                'watts/meter^3', rcode)
c
c --- beam electron profile
c
          id_qbeame = NCVDEF (niterdb, 'qbeame', NCDOUBLE, 2,
     .                        dim_njtime, rcode)
          call NCAPTC (niterdb, id_qbeame, 'long_name', NCCHAR, 28,
     .                'power to electrons from beam', rcode)
          call NCAPTC (niterdb, id_qbeame, 'units', NCCHAR, 13,
     .                'watts/meter^3', rcode)
c
c --- electron ion equilibration profile
c
          id_qdelt = NCVDEF (niterdb, 'qdelt', NCDOUBLE, 2,
     .                       dim_njtime, rcode)
          call NCAPTC (niterdb, id_qdelt, 'long_name', NCCHAR, 26,
     .                'electron-ion equilibration', rcode)
          call NCAPTC (niterdb, id_qdelt, 'units', NCCHAR, 13,
     .                'watts/meter^3', rcode)
c
c --- beam ion profile
c
          id_qbeami = NCVDEF (niterdb, 'qbeami', NCDOUBLE, 2,
     .                        dim_njtime, rcode)
          call NCAPTC (niterdb, id_qbeami, 'long_name', NCCHAR, 23,
     .                'power to ions from beam', rcode)
          call NCAPTC (niterdb, id_qbeami, 'units', NCCHAR, 13,
     .                'watts/meter^3', rcode)
c
c --- RF electron heating profile
c
          id_qrfe = NCVDEF (niterdb, 'qrfe', NCDOUBLE, 2,
     .                      dim_njtime, rcode)
          call NCAPTC (niterdb, id_qrfe, 'long_name', NCCHAR, 19,
     .                'RF electron heating', rcode)
          call NCAPTC (niterdb, id_qrfe, 'units', NCCHAR, 13,
     .                'watts/meter^3', rcode)
c
c --- RF ion heating profile
c
          id_qrfi = NCVDEF (niterdb, 'qrfi', NCDOUBLE, 2,
     .                      dim_njtime, rcode)
          call NCAPTC (niterdb, id_qrfi, 'long_name', NCCHAR, 15,
     .                '-RF ion heating', rcode)
          call NCAPTC (niterdb, id_qrfi, 'units', NCCHAR, 13,
     .                'watts/meter^3', rcode)
c
c --- qione heating profile
c
          id_qione = NCVDEF (niterdb, 'qione', NCDOUBLE, 2,
     .                       dim_njtime, rcode)
          call NCAPTC (niterdb, id_qione, 'long_name', NCCHAR, 65,
     .    'electron power density due to recombination and impact ioniza
     .tion', rcode)
          call NCAPTC (niterdb, id_qione, 'units', NCCHAR, 13,
     .                'watts/meter^3', rcode)
c
c --- qioni heating profile
c
          id_qioni = NCVDEF (niterdb, 'qioni', NCDOUBLE, 2,
     .                       dim_njtime, rcode)
          call NCAPTC (niterdb, id_qioni, 'long_name', NCCHAR, 60,
     .   'ion power density due to recombination and impact ionization',
     .                 rcode)
          call NCAPTC (niterdb, id_qioni, 'units', NCCHAR, 13,
     .                'watts/meter^3', rcode)
c
c --- qxc, ion heating profile
c
          id_qcx = NCVDEF (niterdb, 'qcx', NCDOUBLE, 2,
     .                     dim_njtime, rcode)
          call NCAPTC (niterdb, id_qcx, 'long_name', NCCHAR, 52,
     .        'ion power density due to neutral-ion charge exchange',
     .         rcode)
          call NCAPTC (niterdb, id_qcx, 'units', NCCHAR, 13,
     .                'watts/meter^3', rcode)
c
c --- 2d electron heating profile
c
          id_qe2d = NCVDEF (niterdb, 'qe2d', NCDOUBLE, 2,
     .                      dim_njtime, rcode)
          call NCAPTC (niterdb, id_qe2d, 'long_name', NCCHAR, 19,
     .                '2d electron heating', rcode)
          call NCAPTC (niterdb, id_qe2d, 'units', NCCHAR, 13,
     .                'watts/meter^3', rcode)
c
c --- 2d ion heating profile
c
          id_qi2d = NCVDEF (niterdb, 'qi2d', NCDOUBLE, 2,
     .                      dim_njtime, rcode)
          call NCAPTC (niterdb, id_qi2d, 'long_name', NCCHAR, 14,
     .                '2d ion heating', rcode)
          call NCAPTC (niterdb, id_qi2d, 'units', NCCHAR, 13,
     .                'watts/meter^3', rcode)
c
c --- fusion electron heating  profile
c
          id_qfuse = NCVDEF (niterdb, 'qfuse', NCDOUBLE, 2,
     .                       dim_njtime, rcode)
          call NCAPTC (niterdb, id_qfuse, 'long_name', NCCHAR, 23,
     .                'fusion electron heating', rcode)
          call NCAPTC (niterdb, id_qfuse, 'units', NCCHAR, 13,
     .                'watts/meter^3', rcode)
c
c --- fusion ion heating profile
c
          id_qfusi = NCVDEF (niterdb, 'qfusi', NCDOUBLE, 2,
     .                       dim_njtime, rcode)
          call NCAPTC (niterdb, id_qfusi, 'long_name', NCCHAR, 18,
     .                'fusion ion heating', rcode)
          call NCAPTC (niterdb, id_qfusi, 'units', NCCHAR, 13,
     .                'watts/meter^3', rcode)
c
c --- beam fusion electron heating profile (fraction of beam fusion
c     energy deposited on thermal electron distribution
c
          id_qbfuse = NCVDEF (niterdb, 'qbfuse', NCDOUBLE, 2,
     .                        dim_njtime, rcode)
          call NCAPTC (niterdb, id_qbfuse, 'long_name', NCCHAR, 28,
     .                'beam fusion electron heating', rcode)
          call NCAPTC (niterdb, id_qbfuse, 'units', NCCHAR, 13,
     .                'watts/meter^3', rcode)
c
c --- beam fusion ion heating profile (fraction of beam fusion
c     energy deposited on thermal ion distribution)
c
          id_qbfusi = NCVDEF (niterdb, 'qbfusi', NCDOUBLE, 2,
     .                        dim_njtime, rcode)
          call NCAPTC (niterdb, id_qbfusi, 'long_name', NCCHAR, 23,
     .                'beam fusion ion heating', rcode)
          call NCAPTC (niterdb, id_qbfusi, 'units', NCCHAR, 13,
     .                'watts/meter^3', rcode)
c
c --- mag electron heating profile
c
          id_qmag = NCVDEF (niterdb, 'qmag', NCDOUBLE, 2,
     .                      dim_njtime, rcode)
          call NCAPTC (niterdb, id_qmag, 'long_name', NCCHAR, 21,
     .                'qmag electron heating', rcode)
          call NCAPTC (niterdb, id_qmag, 'units', NCCHAR, 13,
     .                'watts/meter^3', rcode)
c
c --- sawtooth electron heating profile
c
          id_qsawe = NCVDEF (niterdb, 'qsawe', NCDOUBLE, 2,
     .                       dim_njtime, rcode)
          call NCAPTC (niterdb, id_qsawe, 'long_name', NCCHAR, 25,
     .                'sawtooth electron heating', rcode)
          call NCAPTC (niterdb, id_qsawe, 'units', NCCHAR, 13,
     .                'watts/meter^3', rcode)
c
c --- sawtooth ion heating profile
c
          id_qsawi = NCVDEF (niterdb, 'qsawi', NCDOUBLE, 2,
     .                       dim_njtime, rcode)
          call NCAPTC (niterdb, id_qsawi, 'long_name', NCCHAR, 20,
     .                'sawtooth ion heating', rcode)
          call NCAPTC (niterdb, id_qsawi, 'units', NCCHAR, 13,
     .                'watts/meter^3', rcode)
c
c --- radiated power density
c
          id_qrad = NCVDEF (niterdb, 'qrad', NCDOUBLE, 2,
     .                      dim_njtime, rcode)
          call NCAPTC (niterdb, id_qrad, 'long_name', NCCHAR, 22,
     .                'radiated power density', rcode)
          call NCAPTC (niterdb, id_qrad, 'units', NCCHAR, 13,
     .                'watts/meter^3', rcode)
c
c --- qohm,ohmic heating profile
c
          id_qohm = NCVDEF (niterdb, 'qohm', NCDOUBLE, 2,
     .                      dim_njtime, rcode)
          call NCAPTC (niterdb, id_qohm, 'long_name', NCCHAR, 30,
     .                '(electron) ohmic power density', rcode)
          call NCAPTC (niterdb, id_qohm, 'units', NCCHAR, 13,
     .                'watts/meter^3', rcode)
c
          call NCENDF (niterdb, rcode)  ! leave define mode
        end if
c
c --- put values of variables defined above
c
          cnnt(1) = nj
          cnnt(2) = nneu
          cnnt(3) = icall + 1
          const   = 1.0e+6
          call multpl2 (enbeam, tmpdata, nj, const)
          call NCVPT (niterdb, id_enbeam  , c1t, cnt, tmpdata , rcode)
          call multpl2 (enn      ,tmpdata,nj,const)
          call NCVPT (niterdb, id_enn     ,c11t,cnnt, tmpdata , rcode)
          call multpl2 (ennw     ,tmpdata,nj,const)
          call NCVPT (niterdb, id_ennw    ,c11t,cnnt, tmpdata , rcode)
          call multpl2 (ennv     ,tmpdata,nj,const)
          call NCVPT (niterdb, id_ennv    ,c11t,cnnt, tmpdata , rcode)
          call multpl2 (volsn    ,tmpdata,nj,const)
          call NCVPT (niterdb, id_volsn   ,c11t,cnnt, tmpdata , rcode)
          call multpl2 (sbion    ,tmpdata,nj,const)
          call NCVPT (niterdb, id_sbion   , c1t, cnt, tmpdata , rcode)
          call multpl2 (sbeam    ,tmpdata,nj,const)
          call NCVPT (niterdb, id_sbeam   , c1t, cnt, tmpdata , rcode)
          const = 1.0e+4
          call multpl2 (curden, tmpdata, nj, const)
          call NCVPT (niterdb, id_curden  , c1t, cnt, tmpdata , rcode)
          call multpl2 (curohm   ,tmpdata,nj,const)
          call NCVPT (niterdb, id_curohm  , c1t, cnt, tmpdata , rcode)
          call multpl2 (curboot  ,tmpdata,nj,const)
          call NCVPT (niterdb, id_curboot , c1t, cnt, tmpdata , rcode)
          call multpl2 (curdbeam ,tmpdata,nj,const)
          call NCVPT (niterdb, id_curdbeam, c1t, cnt, tmpdata , rcode)
          call multpl2 (currf    ,tmpdata,nj,const)
          call NCVPT (niterdb, id_currf   , c1t, cnt, tmpdata , rcode)
          const = 1.0e-6
          call multpl2 (rbp, tmpdata, nj, const)
          call NCVPT (niterdb, id_rbp     , c1t, cnt, tmpdata , rcode)
          call NCVPT (niterdb, id_zeff    , c1t, cnt, zeff    , rcode)
          call NCVPT (niterdb, id_angrot  , c1t, cnt, angrot  , rcode)
          tmpdata(1:nj) = 0.1*storqueb(1:nj)  ! corrected HSJ, 9/7/07
          call NCVPT (niterdb, id_storqueb , c1t, cnt,tmpdata, rcode)
          call NCVPT (niterdb, id_Kpol_c, c1t, cnt,
     .                  Kpol_c, rcode)
          call NCVPT (niterdb, id_Kpol_exp, c1t, cnt,
     .                  Kpol_exp, rcode)
          call NCVPT (niterdb, id_Kpol_d, c1t, cnt,
     .                  Kpol_d, rcode)
          call NCVPT (niterdb, id_angrot_c, c1t, cnt,
     .                  angrot_c, rcode)
          call NCVPT (niterdb, id_angrot_d, c1t, cnt,
     .                  angrot_d, rcode)
          call NCVPT (niterdb, id_ave_vpar_d, c1t, cnt,
     .                  ave_vpar_d, rcode)
          call NCVPT (niterdb, id_sqz_d, c1t, cnt,
     .                  sqz_d, rcode)
          call NCVPT (niterdb, id_Epsi, c1t, cnt, Epsi, rcode)
          call NCVPT (niterdb, id_Epsi_exp, c1t, cnt, Epsi_exp, rcode)
          call NCVPT (niterdb, id_udia_d, c1t, cnt, udia_d, rcode)
          call NCVPT (niterdb, id_udia_c, c1t, cnt, udia_c, rcode)
          call NCVPT (niterdb, id_ugrt  , c1t, cnt, ugrt, rcode)
          call NCVPT (niterdb, id_cer_bp, c1t, cnt, cer_bp, rcode)
          call NCVPT (niterdb, id_cer_btdr, c1t, cnt, cer_btdr, rcode)
          const = 1.0e-4
          call multpl2 (chieinv, tmpdata, nj, const)
          call NCVPT (niterdb, id_chieinv , c1t, cnt, tmpdata , rcode)
          call multpl2 (chiinv   ,tmpdata,nj,const)
          call NCVPT (niterdb, id_chiinv  , c1t, cnt, tmpdata , rcode)
          const = 100.0
          call multpl2 (xkineo, tmpdata, nj, const)
          call NCVPT (niterdb, id_xkineo  , c1t, cnt, tmpdata , rcode)
          call multpl2 (dpedtc   ,tmpdata,nj,convert)
          call NCVPT (niterdb, id_dpedtc  , c1t, cnt, tmpdata , rcode)
          call multpl2 (dpidtc   ,tmpdata,nj,convert)
          call NCVPT (niterdb, id_dpidtc  , c1t, cnt, tmpdata , rcode)
          call multpl2 (qconde   ,tmpdata,nj,convert)
          call NCVPT (niterdb, id_qconde  , c1t, cnt, tmpdata , rcode)
          call multpl2 (qcondi   ,tmpdata,nj,convert)
          call NCVPT (niterdb, id_qcondi  , c1t, cnt, tmpdata , rcode)
          call multpl2 (qconve   ,tmpdata,nj,convert)
          call NCVPT (niterdb, id_qconve  , c1t, cnt, tmpdata , rcode)
          call multpl2 (qconvi   ,tmpdata,nj,convert)
          call NCVPT (niterdb, id_qconvi  , c1t, cnt, tmpdata , rcode)
          call multpl2 (qbeame   ,tmpdata,nj,convert)
          call NCVPT (niterdb, id_qbeame  , c1t, cnt, tmpdata , rcode)
          call multpl2 (qdelt    ,tmpdata,nj,-convert)
          call NCVPT (niterdb, id_qdelt   , c1t, cnt, tmpdata , rcode)
          call multpl2 (qbeami   ,tmpdata,nj,convert)
          call NCVPT (niterdb, id_qbeami  , c1t, cnt, tmpdata , rcode)
          call multpl2 (qrfe     ,tmpdata,nj,convert)
          call NCVPT (niterdb, id_qrfe    , c1t, cnt, tmpdata , rcode)
          call multpl2 (qrfi     ,tmpdata,nj,convert)
          call NCVPT (niterdb, id_qrfi    , c1t, cnt, tmpdata , rcode)
          call multpl2 (qione    ,tmpdata,nj,-convert)
          call NCVPT (niterdb, id_qione   , c1t, cnt, tmpdata , rcode)
          call multpl2 (qioni    ,tmpdata,nj,convert)
          call NCVPT (niterdb, id_qioni   , c1t, cnt, tmpdata , rcode)
          call multpl2 (qcx      ,tmpdata,nj,-convert)
          call NCVPT (niterdb, id_qcx     , c1t, cnt, tmpdata , rcode)
          call multpl2 (qe2d     ,tmpdata,nj,convert)
          call NCVPT (niterdb, id_qe2d    , c1t, cnt, tmpdata , rcode)
          call multpl2 (qi2d     ,tmpdata,nj,convert)
          call NCVPT (niterdb, id_qi2d    , c1t, cnt, tmpdata , rcode)
          call multpl2 (qfuse    ,tmpdata,nj,convert)
          call NCVPT (niterdb, id_qfuse   , c1t, cnt, tmpdata , rcode)
          call multpl2 (qfusi    ,tmpdata,nj,convert)
          call NCVPT (niterdb, id_qfusi   , c1t, cnt, tmpdata , rcode)
          call multpl2 (qbfuse   ,tmpdata,nj,convert)
          call NCVPT (niterdb, id_qbfuse  , c1t, cnt, tmpdata , rcode)
          call multpl2 (qbfusi   ,tmpdata,nj,convert)
          call NCVPT (niterdb, id_qbfusi  , c1t, cnt, tmpdata , rcode)
          call multpl2 (qmag     ,tmpdata,nj,convert)
          call NCVPT (niterdb, id_qmag    , c1t, cnt, tmpdata , rcode)
          call multpl2 (qsawe    ,tmpdata,nj,convert)
          call NCVPT (niterdb, id_qsawe   , c1t, cnt, tmpdata , rcode)
          call multpl2 (qsawi   ,tmpdata,nj,convert)
          call NCVPT (niterdb, id_qsawi   , c1t, cnt, tmpdata , rcode)
          do i=1,nj
             tmpdata(i) = -ABS (qrad(i)) * convert
          end do
          call NCVPT (niterdb, id_qrad    , c1t, cnt, tmpdata , rcode)
          call multpl2 (qohm    ,tmpdata,nj,convert)
          call NCVPT (niterdb, id_qohm    , c1t, cnt, tmpdata , rcode)
       else
          read (niterdb,  7)  starflag
          read (niterdb, 10) (enbeam(j)       , j=1,nj)
          do jn=1,nneu
               read (niterdb,  7)  starflag
               read (niterdb, 10) (enn(j,jn)       , j=1,nj)
               read (niterdb,  7)  starflag
               read (niterdb, 10) (ennw(j,jn)       , j=1,nj)
               read (niterdb,  7)  starflag
               read (niterdb, 10) (ennv(j,jn)       , j=1,nj)
               read (niterdb,  7)  starflag
               read (niterdb, 10) (volsn(j,jn)       , j=1,nj)
          end do
          read (niterdb,  7)  starflag
          read (niterdb, 10) (sbion(j)       , j=1,nj)
          read (niterdb,  7)  starflag
          read (niterdb, 10) (sbeam(j)       , j=1,nj)
          read (niterdb,  7)  starflag
          read (niterdb, 10) (curden(j)      , j=1,nj)
          read (niterdb,  7)  starflag
          read (niterdb, 10) (curohm(j)      , j=1,nj)
          read (niterdb,  7)  starflag
          read (niterdb, 10) (curboot(j)      , j=1,nj)
          read (niterdb,  7)  starflag
          read (niterdb, 10) (curdbeam(j)      , j=1,nj)
          read (niterdb,  7)  starflag
          read (niterdb, 10) (currf(j)      , j=1,nj)
          read (niterdb,  7)  starflag
          read (niterdb, 10) (rbp(j)       , j=1,nj)
          read (niterdb,  7)  starflag
          read (niterdb, 10) (zeff(j), j=1,nj)
          read (niterdb,  7)  starflag
          read (niterdb, 10) (angrot(j), j=1,nj)
          read (niterdb,  7)  starflag
          read (niterdb, 10) (storqueb(j), j=1,nj)
          read (niterdb,  7)  starflag
          read (niterdb, 10) (Kpol_c(j), j=1,nj)
          read (niterdb,  7)  starflag
          read (niterdb, 10) (Kpol_exp(j), j=1,nj)
          read (niterdb,  7)  starflag
          read (niterdb, 10) (Kpol_d(j), j=1,nj)
          read (niterdb,  7)  starflag
          read (niterdb, 10) (angrot_c(j), j=1,nj)
          read (niterdb,  7)  starflag
          read (niterdb, 10) (angrot_d(j), j=1,nj)
          read (niterdb, 10) (ave_vpar_d(j), j=1,nj)
          read (niterdb,  7)  starflag
          read (niterdb, 10) (sqz_d(j), j=1,nj)
          read (niterdb,  7)  starflag
          read (niterdb, 10) (Epsi(j), j=1,nj)
          read (niterdb,  7)  starflag
          read (niterdb, 10) (Epsi_exp(j), j=1,nj)
          read (niterdb,  7)  starflag
          read (niterdb, 10) (udia_d(j), j=1,nj)
          read (niterdb,  7)  starflag
          read (niterdb, 10) (udia_c(j), j=1,nj)
          read (niterdb,  7)  starflag
          read (niterdb, 10) (ugrt(j), j=1,nj)
          read (niterdb,  7)  starflag
          read (niterdb, 10) (cer_bp(j), j=1,nj)
          read (niterdb,  7)  starflag
          read (niterdb, 10) (cer_btdr(j), j=1,nj)
          read (niterdb,  7)  starflag
          read (niterdb, 10) (chieinv(j)       , j=1,nj)
          read (niterdb,  7)  starflag
          read (niterdb, 10) (chiinv(j)       , j=1,nj)
          read (niterdb,  7)  starflag
          read (niterdb, 10) (xkineo(j)      , j=1,nj)
          read (niterdb,  7)  starflag
          read (niterdb, 10) (dpedtc(j)        , j=1,nj)
          read (niterdb,  7)  starflag
          read (niterdb, 10) (dpidtc(j)        , j=1,nj)
          read (niterdb,  7)  starflag
          read (niterdb, 10) (qconde(j)        , j=1,nj)
          read (niterdb,  7)  starflag
          read (niterdb, 10) (qcondi(j)        , j=1,nj)
          read (niterdb,  7)  starflag
          read (niterdb, 10) (qconve(j)        , j=1,nj)
          read (niterdb,  7)  starflag
          read (niterdb, 10) (qconvi(j)        , j=1,nj)
          read (niterdb,  7)  starflag
          read (niterdb, 10) (qbeame(j)        , j=1,nj)
          read (niterdb,  7)  starflag
          read (niterdb, 10) ( qdelt(j)        , j=1,nj)
          read (niterdb,  7)  starflag
          read (niterdb, 10) (qbeami(j)        , j=1,nj)
          read (niterdb,  7)  starflag
          read (niterdb, 10) (qrfe(j)        , j=1,nj)
          read (niterdb,  7)  starflag
          read (niterdb, 10) (qrfi(j)        , j=1,nj)
          read (niterdb,  7)  starflag
          read (niterdb, 10) ( qione(j)        , j=1,nj)
          read (niterdb,  7)  starflag
          read (niterdb, 10) (qioni(j)        , j=1,nj)
          read (niterdb,  7)  starflag
          read (niterdb, 10) ( qcx(j)        , j=1,nj)
          read (niterdb,  7)  starflag
          read (niterdb, 10) (qe2d(j)        , j=1,nj)
          read (niterdb,  7)  starflag
          read (niterdb, 10) (qi2d(j)        , j=1,nj)
          read (niterdb,  7)  starflag
          read (niterdb, 10) (qfuse(j)        , j=1,nj)
          read (niterdb,  7)  starflag
          read (niterdb, 10) (qfusi(j)        , j=1,nj)
          read (niterdb,  7)  starflag
          read (niterdb, 10) (qbfuse(j)        , j=1,nj)
          read (niterdb,  7)  starflag
          read (niterdb, 10) (qbfusi(j)        , j=1,nj)
          read (niterdb,  7)  starflag
          read (niterdb, 10) (qmag(j)        , j=1,nj)
          read (niterdb,  7)  starflag
          read (niterdb, 10) (qsawe(j)        , j=1,nj)
          read (niterdb,  7)  starflag
          read (niterdb, 10) (qsawi(j)        , j=1,nj)
          read (niterdb,  7)  starflag
          read (niterdb, 10) (qrad(j)         , j=1,nj)
          read (niterdb,  7)  starflag
          read (niterdb, 10) (qohm(j)        , j=1,nj)
       end if
c
c --- convert vectors defined on the npsi (i.e., MHD      ) grid to
c     corresponding quantities on the rho (i.e., transport) grid
c
      tension = 0.0               ! don't use tension option of tspline
      tmax    = 0.0               ! max allowed tension
      bpar(1) = 0.0               ! set boundary conditions on spline
      bpar(2) = 0.0
      bpar(3) = 0.0
      bpar(4) = 0.0
c
c --- take care of a roundoff problem
c
      if (ABS (psival(1)-psir(nj)) .lt. 1.0e-10)  psival(1) = psir(nj)
c
c --- average major radius
c
      if (irwflag .eq. 0) then
        if (icall .eq. 0) then          ! first time, so define dim, var
          call NCREDF (niterdb, rcode)  ! enter define mode
          id_rmajavnpsi = NCVDEF (niterdb, 'rmajavnpsi', NCDOUBLE,
     .                            2, dim_njtime, rcode)
          call NCAPTC (niterdb, id_rmajavnpsi, 'long_name', NCCHAR, 81,
     .        'average major radius of each flux surface evaluated at el
     .evation of magnetic axis', rcode)
          call NCAPTC (niterdb, id_rmajavnpsi, 'units', NCCHAR, 6,
     .                'meters', rcode)
          call NCENDF (niterdb, rcode)  ! leave define mode
        end if
          call tspline (psival,rmajavnpsi,npsi,bpar,cs2spline,kpsi,
     .                  ier,tension,aspline,bspline,cspline,dspline,
     .                  espline,fspline,tmax,psir,work,nj)
          const = 1.0e-2
          call multpl1 (work, nj, const)   ! convert to meters
          call NCVPT (niterdb, id_rmajavnpsi, c1t, cnt, work, rcode)
      else
        read (niterdb,  7)  starflag
        read (niterdb, 10) (rmajavnpsi(j), j=1,nj)
      end if
c
c --- average minor radius
c
      if (irwflag .eq. 0) then
        if (icall .eq. 0) then          ! first time, so define dim, var
          call NCREDF (niterdb, rcode)  ! enter define mode
          id_rminavnpsi = NCVDEF (niterdb, 'rminavnpsi', NCDOUBLE,
     .                            2, dim_njtime, rcode)
          call NCAPTC (niterdb, id_rminavnpsi, 'long_name', NCCHAR,
     .    81, 'average minor radius of each flux surface evaluated at el
     .evation of magnetic axis', rcode)
          call NCAPTC (niterdb, id_rminavnpsi, 'units', NCCHAR,
     .    6, 'meters', rcode)
          call NCENDF (niterdb, rcode)  ! leave define mode
        end if
          call tspline (psival,rminavnpsi,npsi,bpar,cs2spline,kpsi,
     .                  ier,tension,aspline,bspline,cspline,dspline,
     .                  espline,fspline,tmax,psir,work,nj)
          const = 1.0e-2
          call multpl1 (work, nj, const)   ! convert to meters
          call NCVPT (niterdb, id_rminavnpsi, c1t, cnt, work, rcode)
      else
        read (niterdb,  7)  starflag
        read (niterdb, 10) (rminavnpsi(j), j=1,nj)
      end if
c
c --- volume
c
      if (irwflag .eq. 0) then
        if (icall .eq. 0) then          ! first time, so define dim, var
          call NCREDF (niterdb, rcode)  ! enter define mode
          id_psivolp = NCVDEF (niterdb, 'psivolp', NCDOUBLE,
     .                         2, dim_njtime, rcode)
          call NCAPTC (niterdb, id_psivolp, 'long_name', NCCHAR,
     .    27, 'volume of each flux surface', rcode)
          call NCAPTC (niterdb, id_psivolp, 'units', NCCHAR,
     .    8, 'meters^3', rcode)
          call NCENDF (niterdb, rcode)  ! leave define mode
        end if
          call tspline (psival,psivolp,npsi,bpar,cs2spline,kpsi,
     .                  ier,tension,aspline,bspline,cspline,dspline,
     .                  espline,fspline,tmax,psir,work,nj)
          const = 1.0e-6
          call multpl1 (work, nj, const)   ! convert to meters^3
          call NCVPT (niterdb, id_psivolp, c1t, cnt, work, rcode)
      else
        read (niterdb,  7)  starflag
        read (niterdb, 10) (psivolp(j), j=1,nj)
      end if
c
c --- elongation
c
      if (irwflag .eq. 0) then
        if (icall .eq. 0) then          ! first time, so define dim, var
          call NCREDF (niterdb, rcode)  ! enter define mode
          id_elongx = NCVDEF (niterdb, 'elongx', NCDOUBLE,
     .                        2, dim_njtime, rcode)
          call NCAPTC (niterdb, id_elongx, 'long_name', NCCHAR,
     .    31, 'elongation of each flux surface', rcode)
          call NCENDF (niterdb, rcode)  ! leave define mode
        end if
          call tspline (psival,elongx,npsi,bpar,cs2spline,kpsi,
     .                  ier,tension,aspline,bspline,cspline,dspline,
     .                  espline,fspline,tmax,psir,work,nj)
          call NCVPT (niterdb, id_elongx, c1t, cnt, work, rcode)
      else
        read (niterdb,  7)  starflag
        read (niterdb, 10) (elongx(j), j=1,nj)
      end if
c
c --- triangularity
c
      if (irwflag .eq. 0) then
        if (icall .eq. 0) then          ! first time, so define dim, var
          call NCREDF (niterdb, rcode)  ! enter define mode
          id_triangnpsi = NCVDEF (niterdb, 'triangnpsi', NCDOUBLE,
     .                            2, dim_njtime, rcode)
          call NCAPTC (niterdb, id_triangnpsi, 'long_name', NCCHAR,
     .    34, 'triangularity of each flux surface', rcode)
!          call NCENDF (niterdb, rcode)  ! leave define mode
          id_triangnpsi_l = NCVDEF (niterdb, 'triangnpsi_l', NCDOUBLE,
     .                            2, dim_njtime, rcode)
          call NCAPTC (niterdb, id_triangnpsi_l, 'long_name', NCCHAR,
     .    40, 'lower triangularity of each flux surface', rcode)
          call NCENDF (niterdb, rcode)  ! leave define mode
        end if
        call tspline (psival,triangnpsi,npsi,bpar,cs2spline,kpsi,
     .                ier,tension,aspline,bspline,cspline,dspline,
     .                espline,fspline,tmax,psir,work,nj)
        call NCVPT (niterdb, id_triangnpsi, c1t, cnt, work, rcode)
        call tspline (psival,triangnpsi_l,npsi,bpar,cs2spline,kpsi,
     .                ier,tension,aspline,bspline,cspline,dspline,
     .                espline,fspline,tmax,psir,work,nj)
        call NCVPT (niterdb, id_triangnpsi_l, c1t, cnt, work, rcode)
      else
        read (niterdb,  7)  starflag
        read (niterdb, 10) (triangnpsi(j), j=1,nj)
        read (niterdb,  7)  starflag
        read (niterdb, 10) (triangnpsi_l(j), j=1,nj)
      end if
c
c --- indentation
c
      if (irwflag .eq. 0) then
        if (icall .eq. 0) then          ! first time, so define dim, var
          call NCREDF (niterdb, rcode)  ! enter define mode
          id_pindentnpsi = NCVDEF (niterdb, 'pindentnpsi', NCDOUBLE, 2,
     .                             dim_njtime, rcode)
          call NCAPTC (niterdb, id_pindentnpsi, 'long_name', NCCHAR, 32,
     .                'indentation of each flux surface', rcode)
          call NCENDF (niterdb, rcode)  ! leave define mode
        end if
        call tspline (psival, pindentnpsi, npsi, bpar, cs2spline, kpsi,
     .                ier, tension, aspline, bspline, cspline, dspline,
     .                espline, fspline, tmax, psir, work, nj)
        call NCVPT (niterdb, id_pindentnpsi, c1t, cnt, work, rcode)
      else
        read (niterdb,  7)  starflag
        read (niterdb, 10) (pindentnpsi(j), j=1,nj)
      end if
c
c --- surface area
c
      if (irwflag .eq. 0) then
c
****        call tspline (psival,sfareanpsi,npsi,bpar,cs2spline,kpsi,
****  .                   ier,tension,aspline,bspline,cspline,dspline,
****  .                   espline,fspline,tmax,psir,work,nj)
c
c           change definition of surface area calc
c           from using sfareanpsi to using
c
c               4 * pi * pi * R0 * h * r * ABS (grad rho)
c
c           to be consistent with other laboratories:
c
        if (icall .eq. 0) then          ! first time, so define dim, var
          call NCREDF (niterdb, rcode)  ! enter define mode
          id_sfareanpsi = NCVDEF (niterdb, 'sfareanpsi', NCDOUBLE,
     .                            2, dim_njtime, rcode)
          call NCAPTC (niterdb, id_sfareanpsi, 'long_name', NCCHAR,
     .    70, 'surface area of each flux surface, 4*pi*pi*R0*hcap*rho*<A
     .BS(grad rho)>', rcode)
          call NCAPTC (niterdb, id_sfareanpsi, 'units', NCCHAR,
     .    8, 'meters^2', rcode)
          call NCENDF (niterdb, rcode)  ! leave define mode
        end if
        call tspline (psival,grho1npsi,npsi,bpar,cs2spline,kpsi,
     .                ier,tension,aspline,bspline,cspline,dspline,
     .                espline,fspline,tmax,psir,work,nj)
        do j=1, nj
          work(j) = 39.478 * rmajor * hcap(j) * r(j) * ABS (work(j))
        end do
        const = 1.0e-4
        call multpl1 (work, nj, const)   ! convert to meter^2
        call NCVPT (niterdb, id_sfareanpsi, c1t, cnt, work, rcode)
      else
        read (niterdb,  7)  starflag
        read (niterdb,  7)  starflag
        read (niterdb, 10) (sfareanpsi(j), j=1,nj)
      end if
c
c --- cross-sectional area
c
      if (irwflag .eq. 0) then
        if (icall .eq. 0) then          ! first time, so define dim, var
          call NCREDF (niterdb, rcode)  ! enter define mode
          id_cxareanpsi = NCVDEF (niterdb, 'cxareanpsi', NCDOUBLE,
     .                            2, dim_njtime, rcode)
          call NCAPTC (niterdb, id_cxareanpsi, 'long_name', NCCHAR,
     .    41, 'cross-sectional area of each flux surface', rcode)
          call NCAPTC (niterdb, id_cxareanpsi, 'units', NCCHAR,
     .    8, 'meters^2', rcode)
          call NCENDF (niterdb, rcode)  ! leave define mode
        end if
        call tspline (psival,cxareanpsi,npsi,bpar,cs2spline,kpsi,
     .                ier,tension,aspline,bspline,cspline,dspline,
     .                espline,fspline,tmax,psir,work,nj)
        const = 1.0e-4
        call multpl1 (work, nj, const)   ! convert to meter^2
        call NCVPT (niterdb, id_cxareanpsi, c1t, cnt, work, rcode)
      else
        read (niterdb,  7)  starflag
        read (niterdb, 10) (cxareanpsi(j), j=1,nj)
      end if
c
c --- flux surface average grad rho
c
      if (irwflag .eq. 0) then
        if (icall .eq. 0) then          ! first time, so define dim, var
          call NCREDF (niterdb, rcode)  ! enter define mode
          id_grho1npsi = NCVDEF (niterdb, 'grho1npsi', NCDOUBLE,
     .                           2, dim_njtime, rcode)
          call NCAPTC (niterdb, id_grho1npsi, 'long_name', NCCHAR,
     .    38, 'flux surface average absolute grad rho', rcode)
          call NCENDF (niterdb, rcode)  ! leave define mode
        end if
        call tspline (psival,grho1npsi,npsi,bpar,cs2spline,kpsi,
     .                ier,tension,aspline,bspline,cspline,dspline,
     .                espline,fspline,tmax,psir,work,nj)
        call NCVPT (niterdb, id_grho1npsi, c1t, cnt, work, rcode)
      else
        read (niterdb,  7)  starflag
        read (niterdb, 10) (grho1npsi(j), j=1,nj)
      end if
c
c --- flux surface average (grad rho)**2
c
      if (irwflag .eq. 0) then
        if (icall .eq. 0) then          ! first time, so define dim, var
          call NCREDF (niterdb, rcode)  ! enter define mode
          id_grho2npsi = NCVDEF (niterdb, 'grho2npsi', NCDOUBLE,
     .                           2, dim_njtime, rcode)
          call NCAPTC (niterdb, id_grho2npsi, 'long_name', NCCHAR,
     .    34, 'flux surface average (grad rho)**2', rcode)
          call NCENDF (niterdb, rcode)  ! leave define mode
        end if
        call tspline (psival,grho2npsi,npsi,bpar,cs2spline,
     .                kpsi,ier,tension,aspline,bspline,cspline,
     .                dspline,espline,fspline,tmax,psir,work,nj)
        call NCVPT (niterdb, id_grho2npsi, c1t, cnt, work, rcode)
      else
        read (niterdb,  7)  starflag
        read (niterdb, 10) (grho2npsi(j), j=1,nj)
      end if
c
c --- plasma boundary
c
      if (irwflag .eq. 0) then
        if (icall .eq. 0) then          ! first time, so define dim, var
          call NCREDF (niterdb, rcode)  ! enter define mode
          id_nplasbdry = NCVDEF (niterdb, 'nplasbdry', NCLONG, 1,
     .                           idim_time, rcode)
          call NCAPTC (niterdb, id_nplasbdry, 'long_name', NCCHAR,
     .    35, 'number of points on plasma boundary', rcode)
          idim_plasbdry = NCDDEF (niterdb, 'dim_plasbdry',
     .    nplasbdry, rcode)             ! define demision
c
          dim_plasbdrytime(1) = idim_plasbdry
          dim_plasbdrytime(2) = idim_time
          id_rplasbdry = NCVDEF (niterdb, 'rplasbdry', NCDOUBLE, 2,
     .                           dim_plasbdrytime, rcode)
          call NCAPTC (niterdb, id_rplasbdry, 'long_name', NCCHAR,
     .    28, 'r points for plasma boundary', rcode)
          call NCAPTC (niterdb, id_rplasbdry, 'units', NCCHAR,
     .    6, 'meters', rcode)
c
          id_zplasbdry = NCVDEF (niterdb, 'zplasbdry', NCDOUBLE, 2,
     .                           dim_plasbdrytime, rcode)
          call NCAPTC (niterdb, id_zplasbdry, 'long_name', NCCHAR,
     .    28, 'z points for plasma boundary', rcode)
          call NCAPTC (niterdb, id_zplasbdry, 'units', NCCHAR,
     .    6, 'meters', rcode)
          call NCENDF (niterdb, rcode)  ! leave define mode
        end if
c
        call NCVPT1 (niterdb, id_nplasbdry, icall+1, nplasbdry, rcode)
        cbt(1) = nplasbdry
        cbt(2) = icall + 1
        call NCVPT (niterdb, id_rplasbdry, c1t, cbt, rplasbdry, rcode)
        call NCVPT (niterdb, id_zplasbdry, c1t, cbt, zplasbdry, rcode)
      else
        read (niterdb,  7)  starflag
        read (niterdb,  9)  nplasbdry
        read (niterdb,  7)  starflag
        read (niterdb, 10) (rplasbdry(j), j=1,nplasbdry)
        read (niterdb,  7)  starflag
        read (niterdb, 10) (zplasbdry(j), j=1,nplasbdry)
      end if
c
c --- close file
c
  100 call NCCLOS (niterdb, rcode)
c
c     subsequent calls will use existing file since ICALL is SAVEd
c
      icall = icall + 1
      return
c
      end
