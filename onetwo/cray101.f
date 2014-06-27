
      subroutine add_spaces (line)
c
      implicit none
c
      integer   INDEX,i,LEN_TRIM
      character *(*) line
      character line_new*134, delimiter*1

      data      delimiter / '='/
       
        i = INDEX(line,delimiter)
        if(i .eq. 1)then !error, line starts with =
                   print *,'error in parsing input,'
                   print *,'offending line is :'
                   print *,line
                   call STOP('add_spaces',1)
         endif
         line_new(1:i-1) =line(1:i-1)
         line_new(i:i+2) = ' = '
         line_new(i+3:) = line(i+1:)
         line(1:LEN_TRIM(line_new)) = 
     .                line_new(1:LEN_TRIM(line_new))

      return
c
      end




      subroutine bc_print(input)
      USE param
      USe soln
      USE numbrs
      USE bd_condtn,only : bc,bctime,ub,fluxb,
     .    ub_save,ub_rho_edge,bctime_zone
      implicit  integer (i-n), real*8 (a-h, o-z)
       character(len =*) :: input
c      include 'bcon.i'
c      include 'numbrs.i'
       nkkk = nk
       nbctimkk = nbctim
       if(nbctim .le. 0)nbctimkk = 1
       if(nk .le. 0)nkkk = kk
       print *,'nbctim,nk,input =',nbctim,nk,input
       do i = 1, nbctimkk
       do j=1,nkkk
          print *,'i,j,bctime(i),bc(i,j) =',i,j,bctime(i),bc(i,j)
       enddo
       enddo
      return
      end

      subroutine bc_zone(profin,knots,rprofin,
     .      bpar,iprofnbr,ksplin, kbctim,bctime,
     .     nbctim, nz,kj,nj,nout,
     .     name_profile,roa,variable_edge)
c---------------------------------------------------------------------------
c     subroutine determines the edge zone values of profiles at the
c     times given in bctime(implicitly through the second index
c     in profin and rprofin).  The zone extends from j =nz to j=nj
c     Note that different profiles can have different values of nz
c     and for a given profile nz may  vary with time.
c     for a given profile nz is in normalized flux (0,1.]
c     or in grid points
c     Spline the results onto a [0,1] r grid
c INPUT

c OUTPUT:
c variable_grid   =0  no zonal boundary condition for tis profile was specified
c                 =1  zonal boundary conditions given in terms of grid points
c                 =-1                                             normalized flux
c---------------------------------------------------------------------------
      implicit none
      integer knots(*),ksplin,kbctim,nj,j,i,l,
     .        iprofnbr,nout,istat,nbctim,kj,i1,i2,
     .        variable_edge,iflux,igrid,index
      real *8
     .      profin(ksplin,kbctim), rprofin(ksplin,kbctim),
     .            bpar(4,kbctim),roa(*),nz(*),fun,
     .            bctime(*),fnj,diff1,diff2
      real *8 ,dimension(:),allocatable :: yt
      character *(*) name_profile
    
      fnj = Dfloat(nj)
      allocate (yt(1:nj),STAT = istat)
      if(istat .ne. 0)
     .          call allocate_error("yt - bc_zone",0,istat) 
      variable_edge = 0
      iflux = 0
      igrid = 0
      do j=1,nbctim ! check if user set any of the nz values 
         if( nz(j) .lt. fnj )then
              variable_edge = 1
              if(nz(j) .le. 2.) iflux = 1
              if(nz(j) .gt. 2.) igrid = 1
         endif
      enddo
      !not both iflux and igrid can be set:
       if(iflux*igrid .gt. 0)then
          write(*,'(" Error in zone edge boundary input ",/,
     .                   " for profile ",a)')name_profile
          write(*,'("cant mix normalized flux and grid point",/,
     .           " boundary zone information for a given profile")')
          write(*,'(" Error in zone edge boundary input ",/,
     .                   " for profile ",a)')name_profile
          call STOP ('subroutine bc_zone:fix_edge inpt prblm',0)
       endif
      if(variable_edge .eq. 0)return  ! no boundary zone for this case

c       at least one nz(j) was set by the user so we have to process the data at all times.
      do j=1,nbctim
         if(j .eq. 1)then  !first value MUST be set,we will use it to determine
                           !if input is in normalized flux or in grid points
            if(2. .le. nz(j).and. nz(j) .lt. fnj)then 
                  variable_edge = 1         !input is in grid points 
            else if(0.0 .lt. nz(j).and. nz(j)  .lt. 2.0) then  !input is in flux
                  variable_edge = -1       !boundary zone in normalized flux
            else 
               write(*,'(" Error in zone edge boundary input ",/,
     .                   " for profile ",a)')name_profile
               call STOP ('subroutine bc_zone:fix_edge inpt prblm',0)
            endif
         elseif( j .gt. 1 .and. nz(j) .gt. kj )then
            !nz(j) contains default value because it was not
            !set in inone. This means use value at time j-1 which could be grid point or flux
            !as determined above by variable_edge
            nz(j) = nz(j-1)
         elseif( j .gt. 1)then  ! input is set by user for j > 1
            if(variable_edge .eq. 1)then
               if(nz(j) .lt. 2.D0 .or. nz(j) .gt. fnj)then
               write(*,'(" Input error for profile ",a)')name_profile 
               write(*,'(" At bctime point # ",i5)')j
               call STOP ('subroutine bc_zone:fix_edge inpt prblm',0)
               endif
            else
               if(nz(j) .le. 0.0D0 .or. nz(j) .gt. 1.D0)then
                write(*,'(" Input error for profile ",a)')name_profile 
                write(*,'(" At bctime point # ",i5)')j
               call STOP ('subroutine bc_zone:fix_edge inpt prblm',0)
               endif
            endif
         endif
         if(variable_edge .eq. 1) then
               nz(j) = MIN(fnj,nz(j))     !gridpoint input
         else  !variable_edge .eq. -1
               nz(j)  = MIN(1.0D0,nz(j))   !normalized flux input
         endif

         fun = nz(j)
         if(variable_edge  .eq. 1)then
            index = INT(fun)
         else                   ! variable_edge  .eq. -1
             !find the index of the r value closest to the normalized flux value fun:
            call find1 (i1, i2, fun, roa, nj)
            if(i1*i2 .eq. 0)then
               write(*,'(" Error in finding ",1pe12.5,/,
     .               " within range of roa")')fun
               call STOP ('subroutine bc_zone',2)
            endif
              !pick closest value:
            diff1 =ABS(fun -roa(i1))
            diff2 =ABS(fun -roa(i2))
            if(diff1 .ge. diff2)then
               index = i2
            else
               index = i1
            endif
         endif
         call intrp (0,-1, rprofin(1,j), profin(1,j),knots(j), 
     .                                               roa, yt, nj) 
         write(nout,2)bctime(j),name_profile
         write(nout,1)(roa(i),yt(i),i=index,nj)
      enddo
 1    format(2(2x,1pe14.4))
 2    format(' rho and boundary zone  values  At Time ',1pe14.6,
     .       '  for profile   ',a)
      deallocate (yt, STAT = istat)
      if(istat .ne. 0)
     .          call deallocate_error("yt,bc_zone",0,istat)

      !print *,'jmp.den: variable_edge = ',variable_edge !jmp.den
      !pause
      return
      end


cjmp.ibm.start
      subroutine cpu_time_12(seconds)
      implicit none
      real*8 seconds
      call cpu_time(seconds)
      end





c
c
      subroutine driver
c
c
      USE param
      USE fusion
      USE io
      USE neut
      USE transp
      USE solcon
      USE soln
      USE ions
      USE mhdpar
      USE nub  
      USE ext_prog_info
      USE rf
      USE extra
      USE yoka
      USE numbrs
      USE mesh
      USE sourc
      USE machin
      USE geom
      USE events
      USE flags
      USE soln2d
      USE bd_condtn
      USE mixcom
      USE mcgo
      USE flxav
      USE iterdbmd
      USE     io_gcnmp,       ONLY : nout_ncd => nout,ncrt_ncd=>ncrt
      !USE platform !jmp.ibm.par
      USE island
      USE pelcom
      USE  pellet_mod
      USE statistics,         ONLY : print_stats,o12_mon_set,   
     .                               o12_index,start_timer,     
     .                               stop_timer,last_mon_index, 
     .                               descrip_mon,collect_stats

      implicit  integer (i-n), real*8 (a-h, o-z)
c
c --- plotcodeid is passed to graphics postprocessor TRPLOT and checked
c --- there to ensure we don't use incompatible versions of ONETWO and
c --- plot codes. plotcodeid does NOT have to match versid, but
c --- plotcodeid should be updated whenever any changes are made in file
c --- trpltfil (which is written by ONETWO).  versid changes whenever
c --- ONETWO is updated. However if the update does not involve the
c --- plotting package, then plotcodeid should not be changed (if it is
c --- then PREPLT and TRPLOT will have to get the corresponding change
c --- and be recompiled.) If the update of ONETWO involves the plotting
c --- package, then both versid and plotcodeid should be changed to the
c --- same new value, and plot programs must be changed correspondingly.
c
c
c ----------------------------------------------------------------------
c --- main driver routine for ONETWO 1-1/2-D plasma transport code
c ----------------------------------------------------------------------
c




      include 'quit.i'

      include 'pckcom.i'

c
      logical      run_preplt
      real*8 timestart,timeend,timer,tot_cpu_time
      character *132 line
      character*64  inname, outname, qikname, yokname, options
      character*10 intfile

      integer len_preplt_to_run
c
      ten_to_the_50th =  1.0d+50
      iterdbfilename  = 'iterdb'  ! name of file to create if iterdb = 1
c
c ----------------------------------------------------------------------
c     assign I/O unit numbers
c     note: all file open statements should use variable unit numbers and
c           should be preceded by a call getioun. getioun may overwrite
c           the value assigned here in order to find an available unit.
c           Note that the unit number passed on calls to myopen(nunit,...)
c           (and other sub.'s which call open) may be changed on return.
c           All file closes should be preceded by a call giveupus(unit).
c     note: unit numbers nout--nitre must remain sequential,
c           because some DO loops over these unit numbers are used!
c           Do not call getioun on these units; reserve them at build time.
c ----------------------------------------------------------------------
c



      !iounit_weiland = 0
      ntemp   =  4
      nin     =  5
c      ncrt    =  6 ! standard output, special case. do not call getioun
      ncrt    = ncrt_ncd
      nout    =  7 ! 7--9 must remain reserved. do not call getioun
      nout_ncd = nout
      nqik    =  8
      nitre   =  9
      nitrex  = 10
      niterdb = 11 ! unit number for ITER database file (iterfilename)
      n_rf    = 14 ! file rf_interfacew,sub rf_mhd_interface
      !iguess  = 19
      ntrplt  = 20
      neqplt  = 21
      ngreen  = 22
      nbplt   = 23
      nscr    = 24
      nfw     = 25
      nsavsol = 26
      nwguess = 27
      nrguess = 27
      neq     = 28
      nbdep   = 29
      nb_strt  = 30    ! used by time dependent beam module as restart file
      nmix    = 31
      ntweak  = 32
      nupel   = 35
      ncorin  = 36
      nunadas = 38
      nunncl  = 39
      !ioc = 41
      !irtrace = 42
      nouthx  =  0 ! if gt 0, call getioun(nouthx,nouthx) in sub. nbsgxn
      nmcgo1  = 50 ! used if MCGO is run
      nmcgo2  = 51 ! used if MCGO is run
      nmcgo3  = 52 ! used if MCGO is run
      ntcc    = 53 ! used if TESTING_NTCC is set
      !unit    = 54
      ! n54    = 54
      !ntemp   = 61 ! used in  get_beam_data (ufiles_12.f90)
      !nlun    = 64 ! used in  get_beam_data (ufiles_12.f90)
      ndset   = 65
      n66     = 66 ! used in sub.'s fw and rfcd
      !iounit   = 66
      !isecur   = 67
      !iozeff   = 69
      !n71      = 71
      !n72      = 72
      n77     = 77 ! getioun in dump_values; used in sub.'s cycsoltn,
                   ! sorpicard, freebdry, fluxav12, and interface_dump_psi
      !miunit = 77
      nmhddat = 88 ! used in fixed bdry calc, file "test_trnspt_mhd.txt"
      lunnbx  = 76 ! nubeam output of lost particle data
      lunres  = 74 ! nubeam required unit no.
      !iodbgprt = 87
      lun_nubeam = 99 ! for nubeam version 201112 and beyond !HSJ

       io_temp  = 149 ! use this for temporary opening and closing of files
                      ! eqw files that do not  need to be attached otuside of the
                      ! local scope where the file  is defined.
       io_toray_hist = 788
c --- local unit numbers:
c --- unit n7 [never opened!] is used in sub. delsol and in cray401.f
c --- unit iguess [=19, with getioun in myopen] is used in subroutine EQDSKRT
c --- unit ioc [getioun(ioc,41)] is used in subroutine printang
c --- unit irtrace [getioun(irtrace,42)] is used in sub. RAYTRACE
c --- unit igenray [getioun(igenray,43)] is used in sub. call_genray
c --- unit n42 [getioun(n42,42)] is used in sub. init
c --- unit unit [getioun(unit,54)] is used in sub. dumper
c --- unit n54 [getioun(n54,54)] is used in sub. ip_chi2
c --- unit iounit [getioun(iounit,66)] is used in sub. WRTEQDSK
c --- unit isecur [getioun(isecur,67)] is used in subroutine INIT
c --- unit iozeff [69, then getioun(iozeff,iozeff)] is used in sub GET_DENSITY
c --- unit n71 [getioun(n71,71)] is used in subroutine nclmetrics
c --- unit n72 [getioun(n72,72)] is used in subroutine nclmetrics
c --- unit iodbgprt [getioun(iodbgprt,87)] is used in subroutine diffus
c --- unit miunit [getioun(miunit,77)] is used in subroutine cntour2
c --- unit iounit_weiland, 0 in weiland_setup, is passed to etawn8_12 [getioun]
c
c ----------------------------------------------------------------------
c --- process (up to four) command-line parameters (in CTSS-like format)
c ----------------------------------------------------------------------
c
      ! start collecting time info here
       IF(.NOT. o12_mon_set)THEN  
         o12_index = last_mon_index+1
         last_mon_index = o12_index
         o12_mon_set = .TRUE.
         descrip_mon(o12_index) = "Onetwo  elapsed time"
c       write(888,FMT='("o12_index c101 start =",i5)')o12_index ! 88888999
      ENDIF
      start_timer(o12_index) =.TRUE.
      stop_timer(o12_index)  =.FALSE.
      CALL collect_stats(o12_index)
      start_timer(o12_index) =.FALSE.


      call LINKETTE ( nin,  inname,  nout, outname,
     .               nqik, qikname,  9999, yokname,  options)
      call TIMELEFT (isecrem0)
      timestart =timer()
c
      if (yokname(1:7) .ne. 'yokfil ') then
        write  (ncrt, 10) TRIM(yokname)
   10   format (/ ' WARNING: You have requested the "yokfil" output',
     .                     ' file be named "', a, '".'                /
     .            '          However, the "yokfil" file is obsolete',
     .                     ' and is no longer being written.'         /
     .            '          Therefore your request will be ignored.' /)
      end if
c
c ----------------------------------------------------------------------
c --- destroy output files (to be created below)
c ----------------------------------------------------------------------
c
      call DESTROY (outname)
      call DESTROY (qikname)
c
c ----------------------------------------------------------------------
c --- open input file, then read in type of run and machine
c ----------------------------------------------------------------------
c
      codeid   = ' '
      machinei = ' '
c
      call getioun(nin,nin)
      open (unit = nin, file = inname, status = 'OLD', iostat = iostat)
c
      if (iostat .ne. 0) then
        write  (ncrt, 7000)  TRIM(inname), nin
 7000   format (/ ' ---- ERROR:  file "', a, '" on logical unit', i3   /
     .                     14x, 'cannot be opened, therefore ',
     .                          'execution of ONETWO cannot continue;' /
     .                     14x, 'check that this file exists ',
     .                          'and that it is readable and writeable')
        call giveupus(nin)
        call STOP ('subroutine DRIVER: bad or missing "inone" file', 1)
      end if
c
c --- copy file inone to temp file INONE, stripping comments from inone
c
      call getioun(ntemp,ntemp)
      open   (unit = ntemp, file = 'INONE', status = 'UNKNOWN')
      call strip_comments (ncrt, nin, ntemp)
      call giveupus(ntemp)
      close  (unit = ntemp)
      close  (unit = nin  )
      open   (unit = nin  , file = 'INONE', status = 'OLD'    )
c
      read   (nin, '(a)')  line
      rewind (unit = nin)
c
c --- build character variable eqgrdsze (equilibrium grid size)
c
      write (intfile, '(2i3)')  nw, nh
      read  (intfile,   '(a)')  eqgrdsze
c
c ----------------------------------------------------------------------
c --- decode the codeid (code id) and machinei (machine id) values
c ----------------------------------------------------------------------
c
      do j=1,125
        if  (codeid .ne. ' ')  go to 7105
c        if ((line(j  ) .eq. 'd') .and.
c     .      (line(j+1) .eq. 'e') .and.
c     .      (line(j+2) .eq. 'e')) codeid = 'dee'
        if(line(j:j+2) == 'dee') codeid ='dee'
c
c        if ((line(j  ) .eq. 'o') .and.
c     .      (line(j+1) .eq. 'n') .and.
c     .      (line(j+2) .eq. 'e') .and.
c     .      (line(j+3) .eq. 'd') .and.
c     .      (line(j+4) .eq. 'e') .and.
c     .      (line(j+5) .eq. 'e')) codeid = 'onedee'
        if(line(j:j+5) == 'onedee') codeid ='onedee'
c
 7105   if (machinei .ne. ' ')  go to 7115
c        if ((line(j  ) .eq. 'd') .and.
c     .      (line(j+1) .eq. 'o') .and.
c     .      (line(j+2) .eq. 'u') .and.
c     .      (line(j+3) .eq. 'b') .and.
c     .      (line(j+4) .eq. '-') .and.
c     .      (line(j+5) .eq. 'i') .and.
c     .      (line(j+6) .eq. 'i') .and.
c     .      (line(j+7) .eq. 'i')) machinei = 'doub-iii'
        if(line(j:j+7) == 'doub-iii') machinei = 'doub-iii'
c
c        if ((line(j  ) .eq. 'd') .and.
c     .      (line(j+1) .eq. 'i') .and.
c     .      (line(j+2) .eq. 'i') .and.
c     .      (line(j+3) .eq. 'i') .and.
c     .      (line(j+4) .eq. '-') .and.
c     .      (line(j+5) .eq. 'd')) machinei = 'diii-d'
        if(line(j:j+5) == 'diii-d') machinei = 'diii-d'
      end do
c
 7115 if (codeid .eq. ' ') then
        codeid = 'onedee'
        write (ncrt, '(a)')
     .       ' WARNING: no codeid found, defaulting to "onedee"'
      end if
c
      if (machinei .eq. ' ') then
        machinei = 'diii-d'
        write (ncrt, '(a)')
     .       ' WARNING: no machine id found, defaulting to "diii-d"'
      end if
    
      runid = 'Onetwo_run'
      
c
c ----------------------------------------------------------------------
c --- create various output files
c ----------------------------------------------------------------------
c
      open (unit = nqik  , file =  qikname  , status = 'UNKNOWN')
      open (unit = nout  , file =  outname  , status = 'UNKNOWN')
      call getioun(ntrplt,ntrplt)
      open (unit = ntrplt, file = 'trpltfil', status = 'UNKNOWN')
      open (unit = nitre , file = 'runlog'  , status = 'UNKNOWN')
      call getioun(nitrex,nitrex)
      open (unit = nitrex, file = 'isllog'  , status = 'UNKNOWN')
      call getioun(nmix,nmix)
      open (unit = nmix  , file = 'mixfil'  , status = 'UNKNOWN')
      call getioun(ntweak,ntweak)
      open (unit = ntweak, file = 'tweaks'  , status = 'UNKNOWN')
      call getioun(nbplt,nbplt)
      open (unit = nbplt , file = 'bpltfil' , status = 'UNKNOWN')
c
c ----------------------------------------------------------------------
c --- create files used only by twodee version
c           'scratch'  is used for FREYA-type plots for beam
c           'scratch1' is used by the 2D RF fastwave module
c ----------------------------------------------------------------------
c
      if (codeid .ne. 'onedee') then
        call getioun (nscr,nscr)
        open  (unit = nscr   , file = 'scratch' , status = 'UNKNOWN')
        close (unit = nscr)
        open  (unit = nscr   , file = 'scratch1', status = 'UNKNOWN')
        call giveupus (nscr)
        close (unit = nscr)
        call getioun (neqplt,neqplt)
        open  (unit = neqplt , file = 'eqpltfil', status = 'UNKNOWN')
        call getioun(nsavsol,nsavsol)
        open  (unit = nsavsol, file = 'savsol'  , status = 'UNKNOWN',
     .         form = 'UNFORMATTED')
      end if
c
c ----------------------------------------------------------------------
c --- collect some system I/O information to be used later for restart
c --- capability and resetting file pointers if the 1-1/2-d version is
c --- being run
c ----------------------------------------------------------------------
c

      call getio

c
c
c ----------------------------------------------------------------------
c --- get the name of the host on which this program is running
c ----------------------------------------------------------------------
c
cibm.jmp.par      host_name_ext = 'undefined' !jmp.ibm  
cibm.jmp.par      nchars = IHOST_NAME (host_name) !jmp.ibm 
cibm.jmp.par      if (nchars .le. 0) then
cibm.jmp.par        host_name = 'indeterminate' 
cibm.jmp.par        nchars    =  13
cibm.jmp.par      end if
cibm.jmp.par      write (ncrt, '(/ a)') ' host name: ' // host_name(1:nchars)
c
c ----------------------------------------------------------------------
c --- write out, to several files, the version of ONETWO being used for
c --- this run
c ----------------------------------------------------------------------
c
      do i=1,nunits
        iunit = ioftn(i)
        write  (iunit, 8010)  versid, codeid, machinei, nw, nh, kj
 8010   format (' ONETWO code version : ', a           /
     .          ' run type: ', a8, '  machine id: ', a8           /
     .          ' equilibrium grid size, nw = ', i3, ' nh = ', i3 /
     .          ' transport   grid size, kj = ', i3)
        write  (iunit, 8011)
 8011   format (//)
      end do
c
      write (ncrt  , 8010)  versid, codeid, machinei, nw, nh, kj
      write (nitre , 8010)  versid, codeid, machinei, nw, nh, kj
      write (nitre , 8011)
      write (nmix  , 8010)  versid, codeid, machinei, nw, nh, kj
      write (nmix  , 8011)
      write (ntweak, 8010)  versid, codeid, machinei, nw, nh, kj
      write (ntweak, 8011)
c
c ----------------------------------------------------------------------
c --- read and edit input file, then initialize variables
c ----------------------------------------------------------------------
c
      write (ncrt, '(a)')  ' initialization routine called'
      call init (tohmwant, run_preplt)
      write (ncrt, '(a)')  ' initialization completed'




c
c ----------------------------------------------------------------------
c --- destroy files if not needed
c ----------------------------------------------------------------------
c
      if (wmix   .eq. 0.0)  call DESTROY ('mixfil')
      if (ttweak .eq. 0.0)  call DESTROY ('tweaks')
c
c ----------------------------------------------------------------------
c --- create output file for TORAY code (NEVER ACTUALLY IMPLEMENTED)
c ----------------------------------------------------------------------
c
****  if (rf_output .eq. 1)
**** .  open (unit = iotoray, file = 'toray_log', status = 'UNKNOWN')
c
c ----------------------------------------------------------------------
c --- create output file for fastwave module, if FASTWAVE is called
c ----------------------------------------------------------------------
c
      do 3025 model=1,krf
        if (irfmodel(model) .ne. 'fastwave') then
          go to 3025
        else
          call getioun(nfw,nfw)
          open (unit = nfw, file = 'fwout', status = 'UNKNOWN')
          go to 3026
        end if
 3025 continue
c
c ----------------------------------------------------------------------
c --- write out case description to trpltfil, eqpltfil, and bpltfil
c ----------------------------------------------------------------------
c
 3026 call outcas (plotcodeid)
      if (mhdonly .eq. 0)  call pelfil (nupel, nin, versid)
c
c ----------------------------------------------------------------------
c --- time is in units of seconds.  timnew is used by subroutine TPORT.
c --- n counts the number of transport steps
c ----------------------------------------------------------------------
c
      time   = time0
      timnew = time0
      n      = 0
c
c ----------------------------------------------------------------------
c this is a 1-d run (i.e., transport only)
c no equilibrium calculations will be done.
c ----------------------------------------------------------------------
c
      if (codeid .eq. 'onedee') then
        ieq = 0
        call tport
      else
        call runtwo             ! this is a 1-1/2-d run (codeid = 'dee')
      end if
c
c ----------------------------------------------------------------------
c the end; terminate plot files with a flag
c ----------------------------------------------------------------------
c
 
        call finquire(ntrplt) ! 8888999


      do i=3,nunits
        iunit = ioftn(i)
        write  (iunit, 8040)
 8040   format (' the end**')
      end do
c
c     write solution monitor info to plot file:
c
      call write_sol (ntrplt)
c
c     destroy beam plotting file if timbplt was never set:
      if (timbplt(1) .eq. 1.0e+06)  call DESTROY ('bpltfil')

      do 3030 k=1,krf
 3030 if (      codeid .eq. 'onedee' .and. wisl .eq. 0.0 .and.
     .     irfmodel(k) .ne. 'ech')  call DESTROY ('runlog')
c
c --- dump current profile data to file outone if requested
c
      if (tohmwant .gt. -1.e30) then
        call sscurdrv (curden, curboot, curdri, curohm,currf,
     .                 etap, hcap, r, tohmwant, nj, nout)
      end if
      if(write_profiles)then
          call Profiles_for_inone
      endif
      call TIMELEFT (isecrem1)
      timeend =timer()
      elapsed_cpu_time = isecrem0 - isecrem1
                            write (  ncrt, 8051) elapsed_cpu_time
                            write (  nout, 8050) elapsed_cpu_time
                            write (  nqik, 8050) elapsed_cpu_time
      if (  wmix .ne. 0.0)  write (  nmix, 8050) elapsed_cpu_time
      if (ttweak .ne. 0.0)  write (ntweak, 8050) elapsed_cpu_time
      tot_cpu_time  = timeend - timestart
                            write (  ncrt, 8051) tot_cpu_time
                            write (  nout, 8050) tot_cpu_time
                            write (  nqik, 8050) tot_cpu_time
      if (  wmix .ne. 0.0)  write (  nmix, 8050) tot_cpu_time
      if (ttweak .ne. 0.0)  write (ntweak, 8050) elapsed_cpu_time
 8050 format (// ' > elapsed CPU time =', f15.1, ' seconds')
 8051 format (/  ' > elapsed CPU time =', f15.1, ' seconds')
      write  (nqik, 8060)  ineucg, ifreya, (irfcnt(k), k=1,krf)
      write  (nout, 8060)  ineucg, ifreya, (irfcnt(k), k=1,krf)
 8060 format (i6, ' NEUCG2 calls', i5, ' FREYA calls', i5, ' RF calls')
      if (iyoka .eq. 2)  call psumry (nqik)
c
      call giveupus (nupel)
      close (unit = nupel )
      close (unit = nqik  )
      close (unit = nout  )
      close (unit = nitre )
      call giveupus (nitrex)
      close (unit = nitrex)
      call giveupus (ntrplt)
      close (unit = ntrplt)
      call giveupus (nunncl)
      close (unit = nunncl)
c
      if (ABS(timpel(1)- HUGE(1.e0)) .lt. 1.0e20)
     .call DESTROY ('peldat'  )
      call DESTROY ('scratch' )
      if(save_scratch1 .eq. 0)call DESTROY ('scratch1')
      call DESTROY ('scratch2')
      call DESTROY ('savsol'  )
c
c     translate from trpltfil to trpltout,
c     for subsequent use by the TRPLOT code
c

c     get right version of preplt to run: preplt_to_run
      if (run_preplt) then
        call get_preplt(kj,ncrt,nout,
     .          preplt_to_run,len_preplt_to_run,preplt_py)
        write (ncrt, '(   / a)')
     .       ' spawning PREPLT to convert from trpltfil to trpltout...'
        irun_preplt = ISHELL (preplt_to_run(1:len_preplt_to_run))
        if(irun_preplt .eq.0 )call DESTROY('trpltfil')
c       let code finish even if preplt fails.
c        if (ISHELL (preplt_to_run(1:len_preplt_to_run)) .ne. 0)
c     .  call STOP ('subroutine DRIVER: failure of spawned PREPLT', 224)
        write (ncrt, '(53x, a)')
     .       '...returned from PREPLT'
c        call DESTROY ('trpltfil')
      end if

      stop_timer(o12_index)  =.TRUE.
      CALL collect_stats(o12_index)
      stop_timer(o12_index)  =.FALSE.

       CALL print_stats



      if(time_step_too_small)
     . write(ncrt,'("quitting because time step is too small")')
      return
c
      end

      subroutine expchk (a, nra, nca, nr, nc, work)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c --- check matrix a for exp'tal data which is 0 (meaning that it is not
c --- to be used).  Because the data is interpolated in time we require
c --- that a value be zero at all times if it is zero at any one time,
c --- so that the interpolation in time won't produce bogus values.
c
c  input
c  a          data matrix
c  nra        exact row    dim of a
c  nca        exact column dim of a
c  nr         actual number of rows    in a
c  nc         actual number of columns in a
c  work       work vector min length nr
c
c  output
c  a(nr,nc)   data matrix modified for zero data elements as
c             explained above
c ------------------------------------------------------------------ HSJ
c
      dimension  a(nra,nca), work(*)
c
      do i=1,nr
        work(i) = 1.0
      end do
c
      do   j=1,nc
        do i=1,nr
          if (a(i,j) .eq. 0.0)  work(i) = 0.0
        end do
      end do
c
      do   j=1,nc
        do i=1,nr
          a(i,j) = a(i,j)*work(i)
        end do
      end do
c
      return
c
      end

      subroutine getio
c
      USE param
      USE geom
      USE  io
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c subroutine GETIO collects some system I/O information.  ioftn holds
c the Fortran unit number.  iocnum holds the ioc number.  iopntr is used
c later to hold disk pointer values.  the order of files is as
c follows:  outone, qikone, trpltfil, dbugfl (if it is used) and,
c if codeid .ne. 'onedee', eqpltfil and convfl files
c ----------------------------------------------------------------------
c
c      include 'param.i'
c      include 'geom.i'
c      include 'io.i'
c
      nunits   = 4
      ioftn(1) = nout
      ioftn(2) = nqik
      ioftn(3) = ntrplt
      ioftn(4) = nbplt
c


      if (codeid .ne. 'onedee') then
        nunits   = 5
        ioftn(5) = neqplt
      end if
      do index=1,nunits
        iopntr(index) = 0
      end do
      return
c
      end

      subroutine getioun (unitv, unitreq)
c  getioun sets io-unit variable unitv to the value unitret.
c  Depending on freeus, unitret either = input variable unitreq or is output.
      integer unitv, unitreq, unitret
      unitret = unitreq
      call freeus(unitret) ! dummy does nothing( unitret not changed) HSJ
      unitv = unitret
      return
      end

      subroutine knotnum (rspln, n, ksplin)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c --- return the number of knots n (.ge.3) if rspln is set up ok
c --- otherwise return n = 0 (meaning rspln is not correctly set up)
c
      dimension rspln(*)
c
      if (rspln(1) .ne. 0.0)  go to 30
      do j=3,ksplin
        n = j
        if ( ABS(rspln(j)-1.D0) .lt. 2*SPACING(1.d0) )  go to 20
      end do
   30 n = 0
   20 return
c
      end

      subroutine myopen (nunit, name, itype, length, iread)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      character*(*) name
c
      if (itype .ne. 2 .and. itype .ne. 3) then
        call STOP ('subroutine MYOPEN: only ITYPE = 2 is allowed', 4)
      end if
c
      call getioun(nunit,nunit)
      open (unit = nunit, file = name, status = 'OLD', iostat = iostat)
      if (iostat .eq. 0) then
        iread = 1
      else
        call giveupus(nunit)
        iread = 0
      end if
      return
c
      end

      subroutine outcas (plotcodeid)
c
      USE param
      USE fusion
      USE io
      USE ions
      USE neut
      USE solcon
      USE limiter
      USE mhdpar
      USE mhdgrid
      USE nub  
      USE nub2
      USE rf
      USE yoka
      USE numbrs
      USE mesh
      USE sourc
      USE machin
      USE tfact
      USE geom
      USE transp,only:beam_data,use_nubeam
      USE flags
      USE tordlrot
      USE ext_prog_info, only : preplt_py,get_preplt,preplt_to_run
      USE mhdcom
      USE ifs
      USE staebler
      USE flxav
      USE oloss_dat
!      USE  io,                        ONLY : ioftn, iopntr, nunits
      USE weiland
      USE island
      USE zeffcom
      USE mhdbcdtn
      USE pelcom
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c --- output case description to trpltfil, eqpltfil, and bpltfil
c ----------------------------------------------------------------------
c
c     Revisions: 08/29/94 G.M. Staebler - added IANGROT, IFSFLAG to trpltfil
c
c      include 'param.i'
c      include 'mhdpar.i'
c      include 'mhdgrid.i'
c      include 'mhdbcdtn.i'
      include 'co2.i'
c      include 'geom.i'
c      include 'io.i'
c      include 'ions.i'
c      include 'island.i'
c      include 'flags.i'
c      include 'flxav.i'
c      include 'fusion.i'
c      include 'limiter.i'
c      include 'machin.i'
c      include 'mesh.i'
c      include 'mhdcom.i'
c      include 'neut.i'
c      include 'nub.i'
c      include 'nub2.i'
c      include 'numbrs.i'
c      include 'pelcom.i'
c      include 'rf.i'
c      include 'sourc.i'
c      include 'solcon.i'
      include 'sxrcom.i'
c      include 'tfact.i'
c      include 'tordlrot.i'
c      include 'yoka.i'
c      include 'zeffcom.i'
      include 'shay.i'
c      include 'staebler.i'
c      include 'weiland.i'
c      include 'ifs.i'
c
      character  irfmod*8, title(10)*8, plotcodeid*(*)
c
      title(1) = ' '
      prado    = 0
      jsep     = 0
c
c --- edit the input file onto the output files ntrplt, nbplt, neqplt
c
      rewind (unit = nin)
c
 2000 read   (unit = nin, fmt = 8000, end = 2015)  (title(i), i=1,9)
      do k=3,nunits
        kunit = ioftn(k)
        write (kunit, '(1x, 9a8)')  (title(i), i=1,9)
      end do
      go to 2000
c
 2015 do k=3,nunits
        kunit = ioftn(k)
        write (kunit, '(''stop'')')
      end do
c
      nbeam = nbeams
      nubeamo =0
      if(use_nubeam)then
         nubeamo= 1
         nbeam = beam_data%nbeam !beams are defined differently
      endif                      !in nubeam
      !print  out both nbeams and nbeam for preplt
      write (ntrplt,'(a)')  plotcodeid
      write (ntrplt, 8000)  codeid


c     get full executeable name preplt_py
        call get_preplt(kj,ncrt,nout,
     .          preplt_to_run,len_preplt_to_run,preplt_py)
      write (ntrplt, '(a)')  preplt_py(1:len_trim(preplt_py))
      write (ntrplt, 8010)  nj,nion,nprim,nimp,npsi,jsep,ifus,iaslow
      write (ntrplt, 8010)  nbeam,nbeams,iborb ,mf,nubeamo
      write (ntrplt, 8010)  nengn,nterow,jsxr,narray,jco2,jzeff,jhirsh
      write (ntrplt, 8011)  wshay, scsmult, maxtimes, iangrot, ifsflag
      write (ntrplt, 8012)  wweiland,include_weiland
      write (ntrplt, 8012)  dorl_kotch, include_ifs
      write (ntrplt, 8020)  theta
      write (ntrplt, 8010)  (itran(i),i=1,kk)
      write (ntrplt, 8010)  (spec_profiles(i),i=1,kj)
      if (jsxr .eq. 0)  go to 2022
      write (ntrplt, 8000)  (namar(k),k=1,narray)
      write (ntrplt, 8010)  (ndiode(k),k=1,narray)
      write (ntrplt, 8010) ((idiode(i,k),i=1,ndiode(k)),k=1,narray)
 2022 if (jco2 .eq. 0)  go to 2100 
      write (ntrplt, 8010)  nco2
      write (ntrplt, 8020) (rtco2(i),i=1,nco2)
 2100 if (jzeff .eq. 0)  go to 2110
      write (ntrplt, 8010) nzeff
      write (ntrplt, 8020) (rtzeff(i),i=1,nzeff)
 2110 irfmod = no_rf
      do k=1,krf
        if (irfmodel(k) .ne. no_rf)  irfmod = 'rf heat'
      end do
      write (ntrplt, '(a8,i2)')  irfmod, irfech
      if (nterow .eq. 0)  go to 2023
      write (ntrplt, '(10i6 )')  (jterow(i), i=1,nterow)
 2023 write (ntrplt, 8000)  (namep(i),i=1,nprim), (namei(i),i=1,nimp)
****  write (*, '("namep and namei should follow")')
****  write (*, 8000)  (namep(i),i=1,nprim), (namei(i),i=1,nimp)
      rfon1 = 1.0e6
      do k=1,krf
        rfon1 = MIN (rfon1, rfon(k))
      end do
      write (ntrplt, 8020)  prado, rfon1, beamon(1), wisl, timmax
      if ( rfon1 .gt. timmax)  go to 2024
      if (irfech .eq. 0     )  go to 2024
      write (ntrplt, 8010) nrfrad
      write (ntrplt, 8020) (rfrad(j),j=1,nrfrad)
      write (ntrplt, 8020) ylaunch
 2024 if (codeid .eq. 'onedee')  go to 2040
      cconst = 0.01
      call multpl1 (rmhdgrid, nw, cconst)
      call multpl1 (zmhdgrid, nh, cconst)
      xdim  = 0.01*xdim
      ydim  = 0.01*ydim
      redge = 0.01*redge
      write (neqplt, 8030) codeid
      write (neqplt, 8010) nj,nion,nprim,nimp,npsi,ishot
      write (neqplt, 8020) xdim,ydim,redge
      write (neqplt, 8020) (rmhdgrid(ix),ix=1,nw)
      write (neqplt, 8020) (zmhdgrid(iy),iy=1,nh)
      write (neqplt, 8000) mhdmode
      write (neqplt, 7998) use_stark
 7998 format (l10)
      write (neqplt, 8010)                  nlimiter
      write (neqplt, 8020) (xlimiter(j),j=1,nlimiter)
      write (neqplt, 8020) (ylimiter(j),j=1,nlimiter)
      write (neqplt, 8010) nfcoil
c
      cconst = 100.0
      call multpl1 (rmhdgrid, nw, cconst)
      call multpl1 (zmhdgrid, nh, cconst)
      xdim  = 100.0 * xdim
      ydim  = 100.0 * ydim
      redge = 100.0 * redge
c
      
 2040 continue
!         write (nbplt, 8000)  codeid
         write (nbplt, 8001)  codeid,versid ! versid added as of v3.94
         write (nbplt, 8010)  nbeams,mfm1,nj,nw,nh,iborb
         write (nbplt, 8010)  npitch
         write (nbplt, 8020)  xdim,ydim,zax
         write (nbplt, 8020) (rmhdgrid(ix),ix=1,nw)
         write (nbplt, 8020) (zmhdgrid(iy),iy=1,nh)
         write (nbplt, 8020)  rmajor,rminor,rout,rin,atw_beam,kappa
c
      call EMPTY (nbplt )
      call EMPTY (neqplt)
      call EMPTY (ntrplt)
      return
c 
 8001 format(a8,2x,'version =',a)
 8000 format (10a8)
 8010 format (6i10)
 8011 format (2x, 1e12.4, 2x, l5, 2x, i5, 2x, i5, 2x, i5)
 8012 format (2x, 1e12.4, 2x, i5)
 8020 format (6(1pe14.6))  !this format must match other formats for ntrplt
 8030 format (a8,'***') 
c
      end

      subroutine prblcin (time, bctime, nbctim, center, edge, alpha,
     .                    gamma, ierr, nout, ncrt)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c --- get parabolic profile specifiers at t = time by interpolation
c --- from the time-dependent input
c
      dimension  bctime(*), center(*), edge(*), alpha(*), gamma(*)
c
      if (nbctim    .eq. 1  )  return
      if (center(2) .eq. 0.0)  return ! => profile not input as parabola
      do i=1,nbctim-1
        dtestl = time-bctime(i)
        dtestu = time-bctime(i+1)
        if (dtestl*dtestu .le. 0.0) then
          dtimei    = 1.0 / (bctime(i+1)-bctime(i))
          slope     = (center(i+1)-center(i))*dtimei
          center(1) = center(i)+slope*dtestl
          slope     = (edge(i+1)-edge(i))*dtimei
          edge(1)   = edge(i)+slope*dtestl
          slope     = (alpha(i+1)-alpha(i))*dtimei
          alpha(1)  = alpha(i)+slope*dtestl
          slope     = (gamma(i+1)-gamma(i))*dtimei
          gamma(1)  = gamma(i)+slope*dtestl
          go to 10
        end if
      end do
      ierr = 1
      write  (nout, 20)  time, bctime(1), bctime(nbctim)
      write  (ncrt, 20)  time, bctime(1), bctime(nbctim)
   20 format ('  ERROR: time,bctime(1),bctime(nbctim) = ',3(2x,f12.5) /
     .        '  can''t get values for parabolic profile input')
   10 return
c
      end

      subroutine Profiles_for_inone
c-----------------------------------------------------------------------
c     create a file of profiles that can be easily copied into inone
c
c------------------------------------------------------------HSJ-------
      USE param 
      USE solcon
      USE numbrs 
      USE soln, only: ene,te,ti,en,curden
      USE ions
      USE mesh
      USE sourc
      USE tordlrot
      USE bd_condtn
      implicit  integer (i-n), real*8 (a-h, o-z)
c      include 'param.i'
c      include 'invers.i'            !njene,etc.,qradr,qradin
c      include 'mesh.i'              !roa(j)
c      include 'numbrs.i'            !nj,nprim,nimp
c      include 'sourc.i'             !qrad
c      include 'tordlrot.i'          !angrot
      ninone = 101


      if(ksplin .lt. nj)
     .           call STOP('Profiles_for_inone,kslpin too small',0)
      

      call getioun(ninone,ninone)
      open (unit = ninone, file = 'inone_newdat', status = 'UNKNOWN', 
     .                                          iostat = iostat)
 1    format(5(2x,1pe14.4))
        do i=1,nion
          write(ninone,'(2x,"njenp =",i5)')nj
          write(ninone,'(2x,"renpin(1,1) =")')
          write(ninone,1)(roa(j),j=1,nj)
          write(ninone,'(2x,"enp(1,1) =")')
          write(ninone,1)(en(j,i),j=1,nj)
          write(ninone,'(//)')
        enddo
        write(ninone,'(2x,"njene =",i5)')nj
        write(ninone,'(2x,"renein(1,1) =")')
        write(ninone,1)(roa(j),j=1,nj)
        write(ninone,'(2x,"enein(1,1) =")')
        write(ninone,1)(ene(j),j=1,nj)


        write(ninone,'(2x,"njte =",i5)')nj
        write(ninone,'(2x,"rtein(1,1) =")')
        write(ninone,1)(roa(j),j=1,nj)
        write(ninone,'(2x,"tein(1,1) =")')
        write(ninone,1)(te(j),j=1,nj)


        write(ninone,'(2x,"njti =",i5)')nj
        write(ninone,'(2x,"rtiin(1,1) =")')
        write(ninone,1)(roa(j),j=1,nj)
        write(ninone,'(2x,"tiin(1,1) =")')
        write(ninone,1)(ti(j),j=1,nj)


        write(ninone,'(2x,"njcur =",i5)')nj
        write(ninone,'(2x,"rcurdein(1,1) =")')
        write(ninone,1)(roa(j),j=1,nj)
        write(ninone,'(2x,"curdenin(1,1) =")')
        write(ninone,1)(curden(j),j=1,nj)


        write(ninone,'(2x,"rangrot(1,1) =")')
        write(ninone,1)(roa(j),j=1,nj)
        write(ninone,'(2x,"angrotin(1,1) =")')
        write(ninone,1)(angrot(j),j=1,nj)


        write(ninone,'(2x,"njzef =",i5)')nj
        write(ninone,'(2x,"rzeffin(1,1) =")')
        write(ninone,1)(roa(j),j=1,nj)
        write(ninone,'(2x,"zeffin(1,1) =")')
        write(ninone,1)(zeff(j),j=1,nj)


        write(ninone,'(2x,"nqrad =",i5)')nj
        write(ninone,'(2x,"qradr  =")')
        write(ninone,1)(roa(j),j=1,nj)
        write(ninone,'(2x,"qradin(1,1) =")')
        write(ninone,1)(qrad(j),j=1,nj)


        close(unit = ninone)

        return
        end

      subroutine prtspln (nbctim, profin, rprofin, bctime, bpar,
     .                    ksplin, kbctim, nout, knots, profile)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c --- prints input profile summary for splninpt = 'new' option
c
      character  profile*8
      dimension  profin(ksplin,kbctim), rprofin(ksplin,kbctim),
     .           bctime(*), knots(*), bpar(4,kbctim)
c
      write  (nout, 10) profile
   10 format (20x, 'SUMMARY FOR INPUT SPLINE PROFILE ', 2x, a)
c
      do m=1,nbctim
        write  (nout, 50)  profile, bctime(m)
   50   format ('  PROFILE   ', a, ' AT TIME =', 1pe10.3)
        write  (nout, 20)  (rprofin(j,m),j=1,knots(m))
   20   format ('  RHO =', (10(2x, 1pe10.3)))
        write  (nout, 30)  profile, (profin(j,m), j=1,knots(m))
   30   format ('  VALUES AT KNOTS FOR   ', a / (10(2x, 1pe10.3)))
        write  (nout, 40)  (bpar(j,m),j=1,4)
   40 format   ('  BPAR =', 4(2x, 1pe10.3))
      end do
c
      write (nout, '(///)')
      return
c
      end





      subroutine readgren (ierr)
c
c
c***********************************************************************
c**                                                                   **
c**     MAIN PROGRAM:  MHD FITTING CODE                               **
c**                                                                   **
c**     SUBPROGRAM DESCRIPTION:                                       **
c**       subroutine READGREN reads a Green's table for use in ONETWO **
c**                                                                   **
c**     CALLING ARGUMENTS:                                            **
c**       ierr      error flag:  if ierr = 0 no error occured         **
c**       if ierr = 1 on return then the Green's table was not read   **
c**                                                                   **
c**     RECORD OF MODIFICATION:                                       **
c**          10/09/89..........first created                          **
c**                                                                   **
c******************************************************************* HSJ
c
      USE param 
      USE io
      USE mhdpar
      USE mhdgrid
      USE ext_prog_info ,only : nchars_12,onetwo_xsct
      USE mhdcom
      implicit  integer (i-n), real*8 (a-h, o-z)
c      include 'param.i'
c      include 'mhdpar.i'
c      include 'io.i'
c      include 'mhdcom.i'
c      include 'mhdgrid.i'
c
      character  eqgrdsze*7, grname*8, intfile*10, greenfilename*64
c
c
c ----------------------------------------------------------------------
c --- generate the file name for the Green's table
c ----------------------------------------------------------------------
c
      if (greentab .eq. 'default') then
         if (nw .lt. 100) then
           if (nh .lt. 100) then
             write  (intfile, 10)  nw, nh
   10        format ('t', i2, i2)
             read   (intfile, 11)  eqgrdsze
   11        format (a5)
             len = 5
           else if (nh .lt. 1000) then
             write  (intfile, 12)  nw, nh
   12        format ('t', i2, i3)
             read   (intfile, 13)  eqgrdsze
   13        format (a6)
             len = 6
           else
             go to 650
           end if
         else if (nw .lt. 1000) then
           if (nh .lt. 100) then
             write (intfile, 14)  nw, nh
   14        format ('t', i3, i2)
             read (intfile, 13) eqgrdsze
             len = 6
           else if (nh .lt. 1000) then
             write (intfile, 15) nw,nh
   15        format ('t', i3, i3)
             read (intfile, 16) eqgrdsze
   16        format ('t', i3, i3)
             len = 7
           else
             go to 650
           end if
         else
           go to 650
         end if
c
         if      (len .eq. 5) then
           grname = eqgrdsze(1:5) // 'd3d'
         else if (len .eq. 6) then
           grname = eqgrdsze(1:6) // 'd3'
         else if (len .eq. 7) then
           grname = eqgrdsze(1:7) // 'd'
         end if
c
         write  (intfile, 17) grname
         read   (intfile, 17) greentab
      else
         write  (intfile, 17) greentab
         read   (intfile, 17) grname
   17    format (a8)
      end if
c
c ----------------------------------------------------------------------
c --- open the file containing the Green's table and read it.
c --- the file's location is found from the environment variable named
c --- ONETWO, which was obtained earlier (in the beginning of
c --- subroutine INIT)
c ----------------------------------------------------------------------
c
      greenfilename = onetwo_xsct(1:nchars_12) // '/data/' // grname
c
      call getioun(ngreen,ngreen)
      open (unit = ngreen, file = greenfilename, status = 'OLD',
     .      form = 'UNFORMATTED', access = 'SEQUENTIAL', err = 100)
c
c --- basic parameters to be verified in ONETWO
c
      ier = 0
      read (ngreen, end=500, err=600) mw, mh, isize, islpfc
      read (ngreen, end=500, err=600) mfcoil, msilop, nagpr2, mrogow,
     .                                mecoil, mesum, mvesel
      if (    mw .ne. nw     )  ier = ier + 1
      if (    mh .ne. nh     )  ier = ier + 1
      if (mfcoil .ne. nfcoil )  ier = ier + 1
      if (msilop .ne. nsilop )  ier = ier + 1
      if (nagpr2 .ne. magpr2 )  ier = ier + 1
      if (mrogow .ne. nrogow )  ier = ier + 1
      if (mecoil .ne. necoil )  ier = ier + 1
      if ( mesum .ne. nesum  )  ier = ier + 1
      if (mvesel .ne. nvessel)  ier = ier + 1
      if (ier .gt. 0)  go to 200
c
c --- grid vectors
c
      read (ngreen, end=500, err=600)  rmhdgrid, zmhdgrid
c
c --- fcoil and plasma current coupling to grid
c
      read (ngreen, end=500, err=600) gridfc
      read (ngreen, end=500, err=600) gridpc
c
c --- fcoil to fcoil and plasma to fcoil coupling
c
      read (ngreen, end=500, err=600) rfcfc
      read (ngreen, end=500, err=600) rfcpc
c
c --- fcoil to psi loop, mag probe, and Rogowski coupling
c
      read (ngreen, end=500, err=600) rsilfc
      read (ngreen, end=500, err=600) rmp2fc
      read (ngreen, end=500, err=600) rgowfc
c
c --- plasma current to psi loop, mag probe, Rogowski coupling
c
      read (ngreen, end=500, err=600) rsilpc
      read (ngreen, end=500, err=600) rmp2pc
      read (ngreen, end=500, err=600) rgowpc
c
c --- e-coil to psi loop,mag probe,fcoil,ecoil and grid coupling
c
      read (ngreen, end=500, err=600) rsilec
      read (ngreen, end=500, err=600) rmp2ec
      read (ngreen, end=500, err=600) rfcec
      read (ngreen, end=500, err=600) recec
      read (ngreen, end=500, err=600) rsisec
      read (ngreen, end=500, err=600) gridec
c
c --- vessel to psi loop,mag probe,fcoil coupling
c
      read (ngreen, end=500, err=600) rsilvs
      read (ngreen, end=500, err=600) rmp2vs
      read (ngreen, end=500, err=600) rfcvs
c
c --- fcoil, ecoil, and vessel to vessel coupling
c
      read (ngreen, end=500, err=600) rvsfc
      read (ngreen, end=500, err=600) rvsec
      read (ngreen, end=500, err=600) rvsvs
      read (ngreen, end=500, err=600) gridvs
c
c --- mhdgrid is used in transport section as well,
c --- therefore we convert it to CGS here
c
      do 50 j=1,nw
   50 rmhdgrid(j) = 100.0 * rmhdgrid(j)
      do 60 j=1,nh
   60 zmhdgrid(j) = 100.0 * zmhdgrid(j)
      hotw        = rmhdgrid(nw) - rmhdgrid(1)
      hoth        = zmhdgrid(nh) - zmhdgrid(1)
      hotrad      =                rmhdgrid(1)
      call giveupus(ngreen)
      close (unit = ngreen)
c
      return
c
  100 write  (ncrt, 110) grname
      write  (nout, 110) grname
  110 format (' FATAL ERROR: subroutine READGREN reports ' /
     .          14x, 'could not open Green''s table ', a8  /
     .          14x, 'ONETWO will be terminated')
      ierr = 1
      return
  200 ierr = 1
      write  (ncrt, 210)  grname, ier
  210 format (' FATAL ERROR: The Green''s table ',a8,' cannot be used'/
     .        ' with this version of ONETWO due to different'         /
     .        ' parameter settings. There are', i5, ' discrepancies:' )
      write (nout, 210)  grname, ier
      write (ncrt, 220)
  220 format ('  ------ ONETWO -----   ------- Green''s table -----')
      write (nout, 220)
      write (ncrt, 230) nw,mw
      write (ncrt, 240) nh,mh
      write (ncrt, 250) nfcoil,mfcoil
      write (ncrt, 260) nsilop,msilop
      write (ncrt, 270) magpr2,nagpr2
      write (ncrt, 280) nrogow,mrogow
      write (ncrt, 290) necoil,mecoil
      write (ncrt, 300) nesum,mesum
      write (ncrt, 310) nvessel,mvesel
      write (nout, 220)
      write (nout, 230) nw,mw
      write (nout, 240) nh,mh
      write (nout, 250) nfcoil,mfcoil
      write (nout, 260) nsilop,msilop
      write (nout, 270) magpr2,nagpr2
      write (nout, 280) nrogow,mrogow
      write (nout, 290) necoil,mecoil
      write (nout, 300) nesum,mesum
      write (nout, 310) nvessel,mvesel
c
  230 format ('  nw =',i5,10x,i5)
  240 format ('  nh =',i5,10x,i5)
  250 format ('  nfcoil =',i5,10x,i5)
  260 format ('  nsilop =',i5,10x,i5)
  270 format ('  magpr2 =',i5,10x,i5)
  280 format ('  nrogow =',i5,10x,i5)
  290 format ('  necoil =',i5,10x,i5)
  300 format ('  nesum =',i5,10x,i5)
  310 format ('  nvessel =',i5,10x,i5)
c
      return
c
c --- error during read (should not happen)
c
  500 ierr = 1
      write (nout, 320) grname
      write (ncrt, 320) grname
  320 format (/' FATAL ERROR: subroutine READGREN reports'         /
     .           14x, 'end of Green''s table ', a8, ' encountered' /
     .           14x, 'ONETWO will be terminated')
      return
  600 ierr = 1
      write (nout, 330) grname
      write (ncrt, 330) grname
  330 format (/ ' ERROR during reading of Green''s table'             /
     .          '       this error should only occur if the file was' /
     .          '       somehow corrupted.  Get a fresh copy of the'  /
     .          '       Green''s table ', a8, ' and try again')
      return
  650 ierr = 1
      write (nout, 400)  nw, nh
      write (ncrt, 400)  nw, nh
  400 format (/ ' FATAL ERROR in READGREN, nw, nh =', 2(2x,i5))
      return
c
      end

      subroutine sscurdrv (curden, curboot, curdrive, curohm,currf,
     .                      eta,hcap, rho, tohmwant, nj, nout)
c
      USE param 
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c SSCURDRV determines the steady state auxiliary driven current profile
c required to match a given current density and bootstrap current
c Normal usage is to have curbeam and curboot given and then determine
c the required RF current profile shape, assuming that the ohmic
c current is in steady state with total ohmic curren given by
c tohmwant (which can be zero, positive or negative as desired).
c
c ------------------------------------------------------------------ HSJ
c
c      include 'param.i'
      include 'storage.i'
c
      real*8 curden(*),curboot(*),eta(*),hcap(*),rho(*),curdrive(*),
     .       curohm(*), currf(*), ssohm(kj), ssdrive(kj), ssdr(kj),
     .       ssor(kj), tohmwant, curblocal(kj), currfwant(kj)
c
      integer      nj, j, nout
      character*32 label
c
      equivalence (  ssohm  (1), xdum(1))
      equivalence (ssdrive  (1), ydum(1))
      equivalence (   ssdr  (1), zdum(1))
      equivalence (   ssor  (1), sdum(1))
      equivalence (curblocal(1), rdum(1))
      equivalence (currfwant(1), wdum(1))
c
c --- get steady state ohmic current under given profile conditions
c --- convert etap to ohm cm using 8.98755179e11
c
      do j=1,nj
        ssohm(j) = 1.0 / (8.98755179e11 * eta(j) * hcap(j))
      end do
c
c --- get corresponding total ohmic current in steady state
c
      call trapv (rho, ssohm, hcap, nj, totssohm)
      ohmfctr = tohmwant / totssohm
      do j=1,nj
        ssohm(j) = ssohm(j) * ohmfctr
      end do
      call trapv (rho, ssohm, hcap, nj, totssohm)
c
      do j=1,nj
        ssdrive  (j) = curden  (j) - ssohm(j) - curboot(j)
        curblocal(j) = curdrive(j) - currf(j)       ! normally currf = 0
        currfwant(j) = ssdrive (j) - curblocal(j)
        if (curdrive(j) .ne. 0.0) then
          ssdr(j) = ssdrive(j) / curdrive(j)
        else
          ssdr(j) = 1000.0
        end if
        if (curohm(j) .ne. 0.0) then
          ssor(j) = ssohm(j) / curohm(j)
        else
          ssor(j) = 1000.0
        end if
      end do
c
      call trapv (rho, currfwant, hcap, nj, totrfwant)
      totrfwant = totrfwant * 6.28
c
      label = 'curden'
      write (nout, 2)   label
      write (nout, 1)  (rho(j), curden(j), j=1,nj)
      label = 'curohm'
      write (nout, 2)   label
      write (nout, 3)   totssohm
      write (nout, 1)  (rho(j), curohm(j), j=1,nj)
      label = 'curdrive'
      write (nout, 2)   label
      write (nout, 1)  (rho(j), curdrive(j), j=1,nj)
      label = 'curboot'
      write (nout, 2)   label
      write (nout, 1)  (rho(j), curboot(j), j=1,nj)
      label = 'ssohm'
      write (nout, 2)   label
      write (nout, 1)  (rho(j), ssohm(j), j=1,nj)
      label = 'ssdrive'
      write (nout, 2)   label
      write (nout, 1)  (rho(j), ssdrive(j), j=1,nj)
      label = 'ssdr01'
      write (nout, 2)   label
      write (nout, 1)  (rho(j), ssdr(j), j=1,nj)
      label = 'ssor01'
      write (nout, 2)   label
      write (nout, 1)  (rho(j), ssor(j), j=1,nj)
      label = 'currfwant'
      write (nout, 2)   label
      write (nout, 1)  (rho(j), currfwant(j), j=1,nj)
c
c     to use in inone need output in rows:
c
      write (nout, 4)  totrfwant
      write (nout, 5) (currfwant(j), j=1,nj)
      write (nout, 5) (      rho(j), j=1,nj)
c
    1 format (2(2x, 1pe14.4))
    2 format (2x, a)
    3 format (' total steady state ohmic current, amps =', 1pe14.4)
    4 format (          ' total RF current wanted, amps ', 1pe12.4)
    5 format (6(2x, 1pe14.4))
c
      return
c
      end



      subroutine strip_comments (ncrt, input, output)
c
      implicit none
c
      integer   INDEX, ncrt, input, output, nchars, last, i,
     .          tab_loc,         number_of_delimiters
      parameter                  (number_of_delimiters = 3)
      character comment_delimiter(number_of_delimiters)*1, line*134
      character  TAB 
      data      comment_delimiter / ';', '!', '#' /
      TAB = char ( 9 )
c
      do while (.true.)
        read (input, '(a)', end = 10) line
        nchars = LEN_TRIM(line)
        line   = 'X' // line(1:nchars)
        nchars = nchars+1
        if (nchars .gt. 133) then
          write  (ncrt, 20)
   20     format (/ ' ---- ERROR: Your "inone" file contains',
     .                   ' one or more lines' /
     .              12x, ' that are longer than 132 characters.' /
     .              12x, ' Please edit the file so no line exceeds',
     .                   ' 132 characters in length.')
          call STOP ('subroutine STRIP_COMMENTS: long input line', 269)
        else if (nchars .gt. 1) then
          last   = nchars
          i      = 0
          do while (i .lt. number_of_delimiters .and. last .gt. 1)
            i    = i + 1
            last = INDEX (line, comment_delimiter(i)) - 1
            if (last .eq. -1)
     .      last = nchars
            line = line(1:last)
          end do
c         replace  tabs with single blanks:
          tab_loc = 1
          do while(tab_loc .gt. 0)
             tab_loc = INDEX(line,TAB)
             if(tab_loc .ne. 0)line(tab_loc:tab_loc) = ' '
          enddo
c          call add_spaces(line)
          write (output, '(a)') line(2:last)
        end if
      end do
   10 return
c
      end









      subroutine time_grad_avg (nbctim, bctime, profin, rprofin, ipoint,
     .                          prof_initial, ksplin, kbctim, nj, nunit,
     .                          knots, time_start, time_end, time0,
     .                          timmax, roa, profout)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c     subroutine determines the time average of the gradient of the
c     input profile at each rho location for the time interval
c     (time_start,time_end)
c
c --- input
c     nbctim     number of time values in bctime vector
c     bctime(j)  j=1,2,..nbctim the time(s) at which the
c                boundary condition or profiles are given
c                (boundary condition for simulation,profile
c                 for analysis mode)
c     profin(i,j) i=1,2..,knots(j),j=1,2..nbctim,is the profile
c     rprofin(i,j) gives the knot location
c     ipoint       is a pointer to the appropriate bpar vector
c                 see variabel iprofnbr in
c                 subroutine intrp for explanation (pretty ugly
c                 but to bothersome to change).
c     prof_initital(i)  i=1,2,..nj is the initial profile,on the
c                  rho grid
c     nj           number of elements in rho grid vector
c     ksplin       parameter sets size of arrays (see above)
c     kbctim       parameter sets size of arrays (see above
c     nunit        fortran unit number for diagnostic i/o
c                  (set to zero for supression of this output)
c     knots(j)     specifies the number of knots used to model
c                  the profile at time bctime(j)c
c
c     time_start
c     time_end     the time interval over which the averaging is done.
c                  Note that this interval is adjusted if necessary
c                  so that it is a subset of the interval (time0,timmax)
c     time0
c     timmax       the code actually runs in simulation and/or
c                  analysis mode over the time interval (time0,timmax)
c     roa(i)       i=1,2,..nj the normalized rho grid vectorc
c
c  output
c
c     profout(i)   i=1,2,..nj the average value of the gradient of
c                  profin over the time interval specified above
c                  The gradient returned is wrt independent variable roa.
c
c ------------------------------------------------------------------ HSJ
c
      include 'storage.i'
c
      integer nbctim, ksplin, kbctim, knots(*), ipoint, nunit,
     .        nj, ierr, i, j, js, je
c
      real*8  profin(ksplin,kbctim), rprofin(ksplin,kbctim),
     .        bctime(*), prof_initial(*),
     .        profs(kstore), profe(kstore), xprof(kstore),
     .        yprof(kstore), grad_pe(kstore), grad_ps(kstore),
     .        profout(*), roa(*), time_start, time_end, time0, timmax,
     .        dt, dtsq, slope, bintcpt
c
c     set up some local names and temporary storage:
c
      equivalence (xdum(1), xprof  (1))
      equivalence (ydum(1), yprof  (1))
      equivalence (zdum(1), profs  (1))
      equivalence (wdum(1), profe  (1))
      equivalence (sdum(1), grad_pe(1))
      equivalence (tdum(1), grad_ps(1))
c
      ierr = 0
c
c     if rgc is input in inone file or we have only a single time point
c
      if (nbctim .eq. 1) then
        call difydx (roa, prof_initial, profout, nj)
        return
      end if
c
c     first some checks to make sure the following calculations
c     can actually be done:
c
      if ((time_start .lt. bctime(1) .or. time_start .ge.
     .     bctime(nbctim)).or.( time_end .lt. time_start
     .                    .or. time_end .gt. bctime(nbctim))) then
        ierr = 1
        if (nunit .gt. 0) then
          write  (nunit, 1) time_start, time_end,
     .                      bctime(1), bctime(nbctim)
    1     format (' subroutine TIME_GRAD_AVG reports an error',
     .            ' in input specification of start, end times' /
     .            ' for the time average gradient calculation'  /
     .            ' time_start     = ', 1pe12.6                 /
     .            ' time_end       = ', 1pe12.6                 /
     .            ' bctime(     1) = ', 1pe12.6                 /
     .            ' bctime(nbctim) = ', 1pe12.6)
        end if
        cconst = 0.0
        call multpl1 (profout, nj, cconst)
        return    ! return with error flag, ierr = 1
      end if
c
c     next search for the indecies js,je
c     js is infimum  such that bctime(js) .le. time_start
c     je is supremum such that bctime(je) .lt. time_end
c
      js = 0
      call tableintrp (bctime, nbctim, time_start, js)
      je = js
      call tableintrp (bctime, nbctim, time_end, je)
c
      dt   = time_end    - time_start
      if (dt .eq. 0.0)  dt = 1.0e30  ! will produce correct result below
      dtsq = time_end**2 - time_start**2
c
      if (je .eq. js) then ! entire average is within 1 bctime interval
c
c       get the profile at time time_start,profs:
c
        call tspldrive (profin,rprofin,profs,knots,ipoint,xprof,
     .                  yprof,time_start)
c
c       get the gradient of the profile:
c
        call difydx(roa,profs,grad_ps,nj)
c
c       get the profile at time time_end,profe:
c
        call tspldrive (profin,rprofin,profe,knots,ipoint,xprof,
     .                  yprof,time_end)
c
c       get the gradient of the profile:
c
        call difydx(roa,profe,grad_pe,nj)
c
c       integrate from time_start to time_end assuming linear
c       variation in time over this interval.
c
        do i=1,nj
          slope      = (grad_pe(i)-grad_ps(i))/dt
          bintcpt    =  grad_ps(i)-slope*time_start
          profout(i) = (bintcpt*dt+0.5*slope*dtsq)/dt
        end do
c
      else ! time interval for average crosses one or more bctime values
c
c       pick up first part,from start time to next bctime value
c
c       get the profile at time time_start,profs:
c
        call tspldrive (profin,rprofin,profs,knots,ipoint,xprof,
     .                  yprof,time_start)
c
c       get the gradient of the profile:
c
        call difydx (roa, profs, grad_ps, nj)
c
c       get the profile at time bctime(je), profe:
c
        call tspldrive (profin,rprofin,profe,knots,ipoint,xprof,
     .                  yprof,bctime(je))
c
c       get the gradient of the profile:
c
        call difydx(roa,profe,grad_pe,nj)
c
c       integrate from time_start to bctime(je) assuming linear
c       variation in time over this interval.
c
        dt   = bctime(js+1)-time_start
        dtsq = bctime(js+1)**2-time_start**2
        do i=1,nj
          slope      = (grad_pe(i)-grad_ps(i))/dt
          bintcpt    =  grad_ps(i)-slope*time_start
          profout(i) = bintcpt*dt+0.5*slope*dtsq
        end do
c
c       now pick up intervals from bctime(js+1) to bctime(je)
c
        j = js + 1
        do while (bctime(j+1) .le. time_end)
c
c         get the gradient at time bctime(j),profs:
c
          call copya(grad_pe,grad_ps,nj)
c
c         get the profile at time bctime(j+1),profe:
c
          call tspldrive (profin,rprofin,profe,knots,ipoint,xprof,
     .                    yprof,bctime(j+1))
c
c         get the gradient of the profile:
c
          call difydx(roa,profe,grad_pe,nj)
c
          dt   = bctime(j+1)-bctime(j)
          dtsq = bctime(j)**2-bctime(j-1)**2
          do i=1,nj
            slope      = (grad_pe(i)-grad_ps(i))/dt
            bintcpt    =  grad_ps(i)-slope*time_start
            profout(i) = bintcpt*dt+0.5*slope*dtsq+profout(i)
          end do
          j = j + 1
        end do
c
c       pick up the last interval from bctime(i) to time_end
c       profe contains the values at the left end of the interval from above
c
        if (je .ne. j)
     .  call STOP ('subroutine TIME_GRAD_AVG: unspecified problem', 63)
        call copya (grad_pe, grad_ps, nj)
c
c       get the profile at time time_end,profe:
c
        call tspldrive (profin, rprofin, profe, knots, ipoint, xprof,
     .                  yprof, time_end)
c
c       get the gradient of the profile:
c
        call difydx(roa,profe,grad_pe,nj)
        dt   = time_end    - bctime(j)
        dtsq = time_end**2 - bctime(j)**2
        do i=1,nj
          if (dt .gt. 0.0) then
            slope      = (grad_pe(i) - grad_ps(i))/dt
            bintcpt    =  grad_ps(i) - slope*time_start
            profout(i) =  bintcpt*dt + 0.5*slope*dtsq + profout(i)
          end if
          profout(i) = profout(i) / (time_end - time_start)
        end do
      end if
      return
c
      end

      real*8 function timer()
cjmp.ibm       real*4 etime, tarray(2)
cjmp.ibm            timer = etime(tarray)
      call cpu_time_12(timer) !jmp.ibm
      end
cjmp.ibm.end

      subroutine tspldrive (profin, rprofin, profout, knots, ipoint,
     .                      xprof, yprof, avgtime)
c
      USE solcon
c
c ----------------------------------------------------------------------
c serves only to get the time value into common block solcon
c (avoids having to include the common blocks in calling routine)
c ------------------------------------------------------------------ HSJ
c
      USE param 
      implicit  integer (i-n), real*8 (a-h, o-z)
c      include 'param.i'
c      include 'solcon.i'
c
      dimension profin(ksplin,kbctim), rprofin(ksplin,kbctim),
     .          profout(*), knots(*), xprof(*), yprof(*)
c
      time_save = time
      time = avgtime    ! because tsplinew uses parameter time in solcon
      call tsplinew (profin,rprofin,profout,knots,ipoint,xprof,yprof)
      time = time_save
      return
c
      end


       SUBROUTINE finquire(iounit)

        integer iounit
        character* 256 filen

         INQUIRE(UNIT = iounit, NAME = filen) 
         print *,'sub finquire ,filen =',iounit,TRIM(filen)

        return
        end


