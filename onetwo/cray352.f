
      subroutine bilinfw (capr,zz,r1,r2,z1,z2,rmaj,z,prfe,prfi,cur,
     .                    kjr,njr,nvba,nvb,nprim,qerz,qirz,qcrz)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      dimension  rmaj(*),z(*),prfe(kjr,*),prfi(kjr,nvba,*),cur(kjr,*)
c
c     bilinear interpolation to find heat and current source at (r,z)
c
      qerz = 0.0
      qirz = 0.0
      qcrz = 0.0
c
c     return if outside RF channel.
c
      if (.not.(zz .lt. z2 .and. zz .gt. z1 .and.
     .        capr .lt. r1 .and. capr .gt. r2))  return
      do 10 ii=2,njr
        if (capr .ge. rmaj(ii))  go to 20
   10 continue
      ii = njr
   20 i1 = ii-1
      i2 = ii
      do 30 jj=2,nvb
      if (z(jj) .ge. zz)  go to 40
   30 continue
      jj = nvb
   40 j1 = jj-1
      j2 = jj
      a1 = (rmaj(i2)-capr)*(z(j2)-zz)
      a2 = (rmaj(i2)-capr)*(zz-z(j1))
      a3 = (capr-rmaj(i1))*(z(j2)-zz)
      a4 = (capr-rmaj(i1))*(zz-z(j1))
      denom = (rmaj(i2)-rmaj(i1))*(z(j2)-z(j1))
      qerz = (a1*prfe(i1,j1)+a2*prfe(i1,j2)
     .      + a3*prfe(i2,j1)+a4*prfe(i2,j2))/denom
      do 50 k=1,nprim
   50 qirz = (a1*prfi(i1,j1,k)+a2*prfi(i1,j2,k)
     . +a3*prfi(i2,j1,k)+a4*prfi(i2,j2,k))/denom  +qirz
      qcrz = (a1*cur(i1,j1)+a2*cur(i1,j2)
     . +a3*cur(i2,j1)+a4*cur(i2,j2))/denom
      return
c
      end

      real*8 function elip (i, x, ier)
c
      implicit none
c
      integer         i, ier, ks
      real*8          elfn, f, s, x, x1, x2, x3
      common /ellipf/ elfn(201)
c
      s  = 200.0 * x
      ks = s + 1
c
      if (ks .ge. 200) then
        call interfw (f, x, elfn(201), elfn(200), elfn(199),
     .                            1.0,     0.995,     0.99)
      else
        x1 =      0.005 * (ks-1.0)
        x2 = x1 + 0.005
        x3 = x2 + 0.005
        call interfw (f, x, elfn(ks), elfn(ks+1), elfn(ks+2),
     .                            x1,         x2,        x3)
      end if
c
      elip = f
      return
c
      end

      subroutine ell
      USE replace_imsl,                    ONLY : my_mmdele
c
      implicit none
c
      integer         i, ier
      real*8          elfn, s
      common /ellipf/ elfn(201)
c
      elfn(201) = 1.0
      do i=1,200
        s       = (i-1.0) * 0.005
        elfn(i) = my_mmdele (1, s, ier)
      end do
      return
c
      end
      subroutine fastcd (codeid, nscr, eqdskin, freq, pin, nk, ntor,
     .                   pant, gamloss, zion, imass, nj, rho, ne, te,
     .                   zeff, rmajor, hcap, qrfe, qrfi, currf, pefcd,
     .                   pifcd, curfcd)
c
      USE param
      USE ename
      USE solcon, only : time
      USE ext_prog_info ,only : fastcd_path
      implicit none
c
      character*(*) codeid, eqdskin
      integer nscr
      real*8  freq,pin
      integer nk
      real*8  ntor(*),pant(*),gamloss
      real*8  zion, imass
      integer nj
      real*8  rho(*),ne(*),te(*),zeff(*),rmajor,hcap(*)
      real*8  qrfe(*),qrfi(*),currf(*)
      real*8  pefcd,pifcd,curfcd
c
****   integer kprim,kimp,kion,kk,kj,kjm1,kbctim,kb,ke,kbe,kcm,kcmp1,
****  .        ksge,kf,kz,krf,kzrf,krt,kevents,nap
****   integer maxp,kjp,ksplin,kar,kdi,ksymbp
****   integer ktab
c
c      include 'param.i'
c
****  character*64  eqdskfilename, eqdskoldname, eqfile, intfl
c
c      include 'ename.i'
c
c ......................................................................
c
c   local variables:
c
c      character*16  infile
      character*18  infile        !02/03/05 HSJ change infcd file name
                                  !to just save last file  
      character*256  command      
      character*80 line
      real*8  pi
      integer n,j
      real*8  rhon(kj),ne20(kj)
      real*8  fploss,dum(kj)
      real*8  psum,csum,pnorm,cnorm
      integer ISHELL
      integer ncall,strleng,strleng2,strleng1,lenge
      logical ex
      data    ncall /0/
      save    ncall,infile
c
      if (codeid .eq. 'onedee') then
        write (6, '(a/a)') ' FASTCD not implemented for "onedee" mode.',
     .                     ' Create a circular EQDSK file.'
        call STOP ('subroutine FASTCD: mismatch of mode and EQDSK', 223)
      end if

      do j=1,nj
        rhon(j) = rho(j) / rho(nj)
        ne20(j) = ne (j) * 1.0e-14
      end do
c

        if(eqdsk_tdem .ne. 'tdem' ) then
             lenge = LEN_TRIM(eqdskfilename)
             eqfile  = eqdskfilename(1:lenge)
        else
           call wrt_tdem_eqdsk(time,eqfile)
           print *,'tdem eqdsk file created in sub fastcd'
        endif
c
      if (eqfile .eq. 'none')  eqfile = eqdskin
c

      ncall = ncall + 1
      IF(ncall .gt.  1)then
c       write  (infile, '(a6, i4.4)')  'infcd_', ncall
        INQUIRE(FILE = infile, EXIST  = ex)
        IF( ex ) CALL DESTROY (infile)       !keep most recent file only
      ENDIF
      write(infile,'(f12.4)')time 
      infile = 'infcd_'//ADJUSTL(infile)
      call getioun(nscr,nscr)
      open   (unit = nscr, file = infile, status = 'UNKNOWN')
      write  (nscr,'(a)')'eqdsk file name'
      write  (nscr,'(a)')eqfile
      write  (nscr,'(a)')'frequency (Hz)'
      write  (nscr,100)freq
      write  (nscr,'(a)')'total input power (W)'
      write  (nscr,100)pin
      write  (nscr,'(a)')'number of toroidal modes in antenna spectrum'
      write  (nscr,400)nk
      write  (nscr,'(a)')'toroidal mode numbers excited by antenna'
      write  (nscr,100)(ntor(n),n=1,nk)
      write  (nscr,'(a)')
     .       'relative power on each toroidal mode number (A.U.)'
      write  (nscr,100)(pant(n),n=1,nk)
      write  (nscr,'(a)')'parasitic damping rate'
      write  (nscr,100)gamloss
      write  (nscr,'(a)')'charge and mass numbers of primary ions'
      write  (nscr,100)zion,imass
      write  (nscr,'(a)')'number of radial grid points'
      write  (nscr,400)nj
      write  (nscr,'(a)')'normalized rho'
      write  (nscr,100)(rhon(j),j=1,nj)
      write  (nscr,'(a)')'electron density (10^20/m^3)'
      write  (nscr,100)(ne20(j),j=1,nj)
      write  (nscr,'(a)')'electron temperature (keV)'
      write  (nscr,100)(te(j),j=1,nj)
      write  (nscr,'(a)')'zeff'
      write  (nscr,100)(zeff(j),j=1,nj)
      write  (nscr,'(a)')'output format'
      write  (nscr,'(a)')'ONETWO'
  100 format (5e16.9)
  400 format (i5)

      close  (unit = nscr)
      call giveupus(nscr)
c
c      command = 'fastcd -f ' // infile
      strleng1 =  len_trim(fastcd_path)
      print *,'fcd path =',fastcd_path
      if(strleng1 .le. 0)
     .    call STOP('fastcd, no path to code found',1)
      strleng = strleng1 + 4                     ! +3 for ' -f '
      strleng = strleng + len(infile)
      strleng2 = len(command)
      if(strleng2 .ge. strleng)then
            command =  fastcd_path(1:strleng1)//' -f ' // infile
            strleng =  len_trim(command)
      else
            call STOP('fastcd, string length inadequate',1)
      endif
      write (6, '(/a/a/3a)')
     .          ' A FASTCD process will now be spawned.',
     .          ' The exact command is shown between the arrows:',
     .          ' --->', command(1:strleng), '<---'
      if (ISHELL ( command(1:strleng))  .ne. 0)
     .  call STOP ('subroutine FASTCD: failure of spawned FASTCD', 221)
c
      call getioun(nscr,nscr)
      open  (unit = nscr, file = infile, status = 'OLD')
c
  200 read  (nscr, '(a)', err = 500)  line
      if (line(1:7) .ne. '$outfcd')  go to 200
      read  (nscr,'(a)')
      read  (nscr,100)pin
      read  (nscr,'(a)')
      read  (nscr,100)fploss
      read  (nscr,'(a)')
      read  (nscr,100)curfcd
      read  (nscr,'(a)')
      read  (nscr,400)nj
      read  (nscr,'(a)')
      read  (nscr,100)(dum(j),j=1,nj)
      read  (nscr,'(a)')
      read  (nscr,100)(qrfe(j),j=1,nj)
      read  (nscr,'(a)')
      read  (nscr,100)(currf(j),j=1,nj)
      call giveupus(nscr)
      close (unit = nscr)
c
      pefcd = (1.0 - fploss) * pin
      pifcd =  0.0
c
c     get right units, and renormalize power and current densities
c
      do j=1,kj
        qrfe (j) = qrfe (j) * 1.0e-6
        qrfi (j) = 0.0
        currf(j) = currf(j) * 1.0e-4
      end do
      call trapv (rho, qrfe , hcap, nj, psum)
      call trapv (rho, currf, hcap, nj, csum)
      pi   = ACOS (-1.0)
      psum = 4.0 * pi**2 * rmajor * psum
      csum = 2.0 * pi             * csum
      if (psum .eq. 0.0)  psum = 1.0e-12
      if (csum .eq. 0.0)  csum = 1.0e-12
      pnorm =  pefcd / psum
      cnorm = curfcd / csum
      write (6, '(a, 2e15.6 /)') ' pnorm, cnorm = ', pnorm, cnorm
      do j=1,nj
        qrfe (j) = qrfe (j) * pnorm
        currf(j) = currf(j) * cnorm
      end do
      return
c
  500 call giveupus(nscr)
      close (unit = nscr)
      call STOP ('subroutine FASTCD: errors in spawned FASTCD', 222)

      return
c
      end

      subroutine fastwave (codeid, nscr,nfw,rmajor,kappa,b0,rbp,r,
     .                     ra,nj,nspfw,te,ti,ene,en,zeff,zii,atw,
     .                     freq,rnp,rfpow,impath,iswch,li,nih,hts,
     .                     zrf,nzrf,pzrf,npsi,x,nx,y,ny,cspln,k2cspln,
     .                     qrfe,qrfi,currf,totfwpe,totfwpi,totfwc)
c
      USE param
      USE contour
        USE mhdpar ,only : kpsi
      USE replace_imsl,                 ONLY : my_dbcevl1
      implicit  integer (i-n), real*8 (a-h, o-z)
c
****  real*8        kappa, codeid
      real*8        kappa                     ! Y. R. Lin-liu
      character*(*) codeid                    ! Y. R. Lin-liu
c
c    Fast wave energy is incident at the outboard side of the tokamak
c    and is restricted to up/down symmetric channels about zrf(i),i = 1,nzrf.
c    Each channel injects fraction pzrf(i),i = 1,nzrf, of RF power rfpow(watts).
c    The channel has full height hts(1) at the outboard edge,
c    hts(2) at r = rmajor, and hts(3) at the inboard edge
c    of the plasma.
c    Heights are the same for each channel.
c    k-parallel varies as 1/R
c
c    Output is:
c      qrfe(j), j = 1,nj,  watts/cm**3 to electrons from all channels.
c      qrfi                             ions
c      currf          ,  amps/cm**2 flux surface averaged RF current.
c                                    from all channels.
c
c     totfwpe(k), k=1,nzrf, power (watts) absorbed by electrons in each channel
c     totfwpi                   ....                  ions      ....
c     totfwc (k), k=1,nzrf, current (amps) driven by each channel
c
c     Use parameters          from ONETWO:  kj, kprim, kf, kzrf
c     kf changed to kpsi 2/10/03 HSJ
c     Use parameter 'nconmax' from ONETWO INCLUDE file 'contour.i'
c
c      include 'param.i'
c      include 'contour.i'
c
      parameter (kjr = 2*(kj-1))
      dimension r(*),ra(*),te(*),ti(*),ene(*),en(kj,kprim),zeff(*),
     .  rbp(*),li(*),nih(*),atw(*),hts(3),qrfe(*),qrfi(*),currf(*)
      dimension rmaj(kjr),ener(kjr),enir(kjr,kprim),ter(kjr),tir(kjr),
     .  zeffr(kjr),zii(kj,kprim),htslb(kjr),hrf(kjr,2),psirr(kj),
     .  ziir(kjr,kprim)
      dimension zrf(nzrf),pzrf(nzrf),freq(nzrf),rnp(nzrf),
     .          totfwpe(nzrf),totfwpi(nzrf),totfwc(nzrf)
      dimension x(nx),y(ny),cspln(k2cspln,nx,*),pds(6)
      dimension bpinv(nconmax),arclen(nconmax),xp(nconmax),yp(nconmax),
     .          work1(nconmax),work2(nconmax),work3(nconmax)
      dimension rho(kpsi),qep(kpsi),qip(kpsi),qcp(kpsi)
c
c     Hard code the number of vertical blocks (use an odd number)
c
      parameter (nvba = 17)
      dimension  z(nvba),cur(kjr,nvba),prfe(kjr,nvba),
     .           prfi(kjr,nvba,kprim)
c
c --- introducing two new arrays, sgnrnp and rnp12, to take care of
c --- the case when rnp is negative
c
      dimension  sgnrnp(kzrf), rnp12(kzrf)
c
      one = 1.0
c
      do nnzrf=1,nzrf
        rnp12 (nnzrf) =            rnp(nnzrf)
        sgnrnp(nnzrf) = SIGN (one, rnp(nnzrf))
        rnp   (nnzrf) =  ABS (     rnp(nnzrf))
      end do
c
      sgnb0 = SIGN (one, b0)
      b0    =  ABS (     b0)
c
c number of major radius bins is twice number of radial bins in
c transport code
c
      njr  = 2*(nj-1)
      njra = kjr
      if (nj .gt. kj)  call STOP ('subroutine FASTWAVE: NJ > KJ', 72)
c
c number of vertical bins
c
      nvb = nvba
      if (nzrf .gt. kzrf)
     .  call STOP ('subroutine FASTWAVE: NZRF > KZRF', 73)
c
      do j=1,nj
        qrfe (j) = 0.0
        qrfi (j) = 0.0
        currf(j) = 0.0
        currf(j) = 0.0
      end do
c
      twopi = 4.0 * ATAN2 (1.0, 0.0)
      b02   = b0**2
c
c     Write heading in "fwout" file
c
      write (nfw, 8131)
c
      if (codeid .eq. 'onedee') then    ! big onedee IF block
        rootkap = SQRT (kappa)
        capg    = 0.5 * (kappa+1.0/kappa)
c
c       Loop over channels
c
        do nnzrf=1,nzrf
c
        wkz   = rnp(nnzrf)*twopi*freq(nnzrf)/3.0e10
        if (zrf(nnzrf) .ge. rootkap*r(nj))
     .  call STOP ('subroutine FASTWAVE: ZRF(NNZRF) outside plasma', 74)
        theta = SIN (zrf(nnzrf)/(r(nj)*rootkap))
        xrf   = r(nj) * COS (thet) / rootkap
        rrf1  = rmajor-xrf
        rrf2  = rmajor+xrf
        dr    = (rrf1-rrf2)/njr
c
        do j=1,nj-1
          rmaj(j) = rrf1-(j-0.5)*dr
          xrf     = rmaj(j)-rmajor
          rr      = SQRT (kappa*xrf**2+zrf(nnzrf)**2/kappa)
          call find (j1, j2, rr, r, nj)
          if (j1 .ne. j2) then
            del = (rr-r(j1))/(r(j2)-r(j1))
          else
            del = 0.0
          end if
          ener(j) = ene(j1)+(ene(j2)-ene(j1))*del
          ter(j) = (te(j1)+(te(j2)-te(j1))*del) * 1000.0
          tir(j) = (ti(j1)+(ti(j2)-ti(j1))*del) * 1000.0
          zeffr(j) = zeff(j1)+(zeff(j2)-zeff(j1))*del
          do k=1,nspfw
            ziir(j,k) = zii(j1,k)+(zii(j2,k)-zii(j1,k))*del
            enir(j,k) = en(j1,k)+(en(j2,k)-en(j1,k))*del
          end do
          rmaj(njr+1-j) = rmaj(j)
          ener(njr+1-j) = ener(j)
          ter(njr+1-j) = ter(j)
          tir(njr+1-j) = tir(j)
          zeffr(njr+1-j) = zeffr(j)
          do k=1,nspfw
            enir(njr+1-j,k) = enir(j,k)
            ziir(njr+1-j,k) = ziir(j,k)
          end do
        end do
c
c       Set up vertical bins for fw
c
        htmax = MAX (hts(1),hts(2),hts(3))
        do i=1,nj-1
          htslb(i) = htmax/nvb
          hrf(i,1) = 0.5*(hts(2)+(rmaj(i)-rmajor)*(hts(1)-hts(2))/
     .             (rmaj(1)-rmajor))
          hrf(i,2) = -hrf(i,1)
        end do
        do i=nj,njr
          htslb(i) = htmax/nvb
          hrf(i,1) = 0.5*(hts(3)+(rmaj(i)-rmaj(njr))*(hts(2)-hts(3))/
     .                            (rmajor-rmaj(njr)))
          hrf(i,2) = -hrf(i,1)
        end do
c
        call fw(njr,njra,rmaj,pzrf(nnzrf)*rfpow,ener,enir,nspfw,ter,tir,
     .   zeffr,iswch,rmajor,b0,freq(nnzrf),wkz,li,nih,ziir,atw,nvb,nvba,
     .   htslb,hrf,prfe,prfi,cur,totfwpe(nnzrf),totfwpi(nnzrf),
     .   totfwc(nnzrf))
c
c       onedee case --- flux surface averages of prfe,prfi,cur
c       set number of points for flux surface averaging
c
        ntheta = 120
        dthet  = twopi/ntheta
        do k=1,nvb
          z(k) = zrf(nnzrf)-0.5*htmax+(k-0.5)*htmax/nvb
        end do
c
        do j=1,nj
          qe = 0.0
          qi = 0.0
          qc = 0.0
          qt = 0.0
          bp02 = (rbp(j)/(r(j)*capg))**2
c
          do i=1,ntheta
            thet = (i-0.5)*dthet
            c = COS (thet)
            s = SIN (thet)
            xx = r(j)*c/rootkap
            zz = r(j)*s*rootkap
            capr = rmajor+xx
            call bilinfw (capr, zz, rrf1, rrf2, z(1), z(nvb), rmaj, z,
     .                    prfe, prfi, cur, kjr, njr, nvba, nvb, nspfw,
     .                    qerz, qirz, qcrz)
            qe = qe + capr*qerz
            qi = qi + capr*qirz
            qc = qc + rmajor
     .              * SQRT (1.0+bp02/b02*(kappa*c*c+s*s/kappa))*qcrz
            qt = qt + capr
          end do
c
          qrfe(j) = qe / qt + qrfe(j)
          qrfi(j) = qi / qt + qrfi(j)
c
c         multiply contribution of the particular ray to currf by sgnrnp
c
          currf(j) = qc/qt*sgnrnp(nnzrf)+currf(j)
        end do
c
c       Write output to file "fwout"
c
        rat0 = pzrf(nnzrf)*rfpow
c
        if (rat0 .ne. 0.0) then
          rat1 = 100.0*totfwpe(nnzrf)/rat0
        else
          rat1 = 0.0
        end if
c
        if (totfwpe(nnzrf) .ne. 0.0) then
           rat2 = totfwc(nnzrf)/totfwpe(nnzrf)
        else
           rat2 = 0.0
        end if
c
        if (rat0 .ne. 0.0) then
          rat3 = totfwc(nnzrf)/rat0
        else
          rat3 = 0.0
        end if
c
****    write (nfw, 8130)  nnzrf, freq(nnzrf), rnp  (nnzrf), rat0,
        write (nfw, 8130)  nnzrf, freq(nnzrf), rnp12(nnzrf), rat0,
     .                     rat1, rat2, rat3
c
        end do
c
c       reset rnp back to the input value
c
        do nnzrf=1,nzrf
          rnp(nnzrf) = rnp12(nnzrf)
        end do
c
        return
      end if    ! end of big onedee IF block
c
 8130 format (i3, 2x, 6(1pe12.2))
 8131 format ('nnzrf','    freq(Hz)  ','  n_para   ','pwr_in(Mw) ',
     .        ' % damping ',' I/P_abs(A/W)',  ' I/P_in(A/W)  ')
c
c     *************  2D CASE  **********************************
c
      do j=1,npsi
        qep(j) = 0.0
        qip(j) = 0.0
        qcp(j) = 0.0
      end do
c
c     Set up FW major radial bins, for 2D case
c     open 'scratch1' file written in FLUXAV and ROWSET
c
      call getioun(nscr,nscr)
      open (unit = nscr, file = 'scratch1', status = 'OLD')
c
c     read through 'scratch1' file to get rho and psir data:
c
      do 106 j=1,npsi
        read (nscr, 8100) mp
        if (mp .eq. 0)  go to 106
        if (mp .gt. nconmax) then
          call giveupus(nscr)
          close (unit = nscr)
          call STOP ('subroutine FASTWAVE: MP > NCONMAX', 75)
        end if
        read (nscr, 8120)  (xp(i),yp(i), i=1,mp)
        read (nscr, 8120)  (bpinv(i)   , i=1,mp)
        read (nscr, 8120)  (arclen(i)  , i=1,mp)
  106 continue
c
      read (nscr, 8100)  nnpsi,nnj          !just space down 
      read (nscr, 8120)   (rho(j),j=1,npsi) !just spaces down
      read (nscr, 8120)     rho(1),rho(2)    !just spaces down
      read (nscr, 8120)  (rho(j),j=1,npsi)
c
c     read unadjusted psir-array, i.e., psi values referenced to eqdsk
c     and evaluated at the r(j),j = 1,nj-values.
c
      read (nscr, 8120)  (psirr(j),j=1,nj)
c
c  Loop over RF channels
c
      do nnzrf=1,nzrf
c
      wkz = rnp(nnzrf)*twopi*freq(nnzrf)/3.0e10
      rewind (unit = nscr)
c
c     Read in xp,yp for outer flux surface
c
      read (nscr, 8100) mp
      if (mp .gt. nconmax) then
        call giveupus(nscr)
        close (unit = nscr)
        call STOP ('subroutine FASTWAVE: MP > NCONMAX', 76)
      end if
      read   (nscr, 8120)  (xp(i),yp(i), i=1,mp)
 8100 format (i6,2x,i6)
 8120 format (6e12.5)
c
c     Find intersection of boundary with z = zrf(nnzrf)
c
      rrf1 = 0.0
      rrf2 = 1.0d+100
c
      do i=1,mp-1
        if ((yp(i+1)-zrf(nnzrf))*(yp(i)-zrf(nnzrf)) .le. 0.0) then
          rrf  = xp(i) + (xp(i+1) - xp(i)) / (yp(i+1) - yp(i))
     .                                  * (zrf(nnzrf) - yp(i))
          rrf1 = MAX (rrf1, rrf)
          rrf2 = MIN (rrf2, rrf)
        end if
      end do
c
c     Divide (rrf1,rrf2) into njr bins, and obtain plasma parameters therein
c
      dr      = (rrf1-rrf2)/njr
      ig      = 0.5*nj
c
      do j=1,njr
        rmaj(j) = rrf1-(j-0.5)*dr
c
c       Find corresponding psi value (dps(1))
c
        call my_dbcevl1 (x,nx,y,ny,cspln,nx,rmaj(j),zrf(nnzrf),
     .                  pds,ier,6)
        if (ier .ne. 0)
     .    call STOP ('subroutine FASTWAVE: error in DBCEVL', 77)
        call find2 (psirr, pds(1), nj, ig, j1, j2, 6, 6)
c
c       Interpolate density, etc., linearly in psi between tabulated points
c
        del      = (pds(1)-psirr(j1))/(psirr(j2)-psirr(j1))
        ener(j)  = ene(j1)+(ene(j2)-ene(j1))*del
        ter(j)   = (te(j1)+(te(j2)-te(j1))*del)*1000.0
        tir(j)   = (ti(j1)+(ti(j2)-ti(j1))*del)*1000.0
        zeffr(j) = zeff(j1)+(zeff(j2)-zeff(j1))*del
        do k=1,nspfw
          ziir(j,k) = zii(j1,k)+(zii(j2,k)-zii(j1,k))*del
          enir(j,k) = en(j1,k)+(en(j2,k)-en(j1,k))*del
        end do
      end do
c
c     Set up vertical heights of RF channel
c
      htmax = MAX (hts(1),hts(2),hts(3))
      do i=1,nj-1
        htslb(i) = htmax/nvb
        hrf(i,1) = 0.5*(hts(2)+(rmaj(i)-rmajor)*(hts(1)-hts(2))/
     .                         (rmaj(1)-rmajor))
        hrf(i,2) = -hrf(i,1)
      end do
c
      do i=nj,njr
        htslb(i) = htmax/nvb
        hrf(i,1) = 0.5*(hts(3)+(rmaj(i)-rmaj(njr))*(hts(2)-hts(3))/
     .                          (rmajor-rmaj(njr)))
        hrf(i,2) = -hrf(i,1)
      end do
c
      do k=1,nvb
        z(k) = zrf(nnzrf) - 0.5 * htmax + (k - 0.5) * htmax / nvb
      end do
c
      call fw (njr,njra,rmaj,pzrf(nnzrf)*rfpow,ener,enir,nspfw,ter,tir,
     .  zeffr,iswch,rmajor,b0,freq(nnzrf),wkz,li,nih,ziir,atw,nvb,nvba,
     .  htslb,hrf,prfe,prfi,cur,totfwpe(nnzrf),totfwpi(nnzrf),
     .  totfwc(nnzrf))
c
c --- account for multiple path absorption
c --- implemented for 2d case only (YRLL)
c
      pwrabs = totfwpe(nnzrf)+totfwpi(nnzrf)
      cfmpow = 1.0
      cfmcur = 1.0
      if ((impath .ne. 1) .and. (pwrabs .ne. 0.0)) then
        rcorr = 1.0 - 1.0/impath
        xcorr = pwrabs/(pzrf(nnzrf)*rfpow)
        cfmpow = 1.0 / xcorr
        cfmcur = 1.0 / ((1.0-rcorr)+rcorr*xcorr)
      end if
c
c     2d flux surface averages
c
c     Read data from 'scratch1' as written by FLUXAV and RHOSET
c
      do 130 j=1,npsi
      if (j .gt. 1) then
        read (nscr, 8100) mp
        if (mp .eq. 0)  go to 130
        if (mp .gt. nconmax) then
          call giveupus(nscr)
          close (unit = nscr)
          call STOP ('subroutine FASTWAVE: MP > NCONMAX', 78)
        end if
        read (nscr, 8120)  (xp(i),yp(i), i=1,mp)
      end if
c
      read (nscr, 8120)  (bpinv (i), i=1,mp)
      read (nscr, 8120)  (arclen(i), i=1,mp)
      call integ (arclen,bpinv,mp,sum1)
      do i=1,mp
        call bilinfw (xp(i), yp(i), rrf1, rrf2, z(1), z(nvb), rmaj, z,
     .                prfe, prfi, cur, kjr, njr, nvba, nvb, nspfw,
     .                work1(i), work2(i), work3(i))
        work1(i) = work1(i)*bpinv(i)
        work2(i) = work2(i)*bpinv(i)
        work3(i) = work3(i)*bpinv(i)
      end do
      call integ (arclen,work1,mp,sum)
      qep(j) = sum/sum1*cfmpow + qep(j)   ! cfmpow: correction factor
      call integ (arclen, work2, mp, sum) ! for multiple path absorption
      qip(j) = sum/sum1*cfmpow+qip(j)
      call integ (arclen,work3,mp,sum)
c
c --- multiply contribution of the particular ray to currf by sgnrnp
c
      qcp(j) = sum/sum1*sgnrnp(nnzrf)*cfmcur+qcp(j)
c
  130 continue
c
c     Write output to file "fwout"
c
      rat0 = pzrf(nnzrf)*rfpow
      if (rat0 .ne. 0.0) then
        rat1 = 100.0 * totfwpe(nnzrf)/rat0
      else
        rat1 = 0.0
      end if
      if (totfwpe(nnzrf) .ne. 0.0) then
        rat2 = totfwc(nnzrf)/totfwpe(nnzrf)
      else
        rat2 = 0.0
      end if
      if (rat0 .ne. 0.0) then
        rat3 = totfwc(nnzrf)/rat0
      else
        rat3 = 0.0
      end if
****  write (nfw, 8130)  nnzrf, freq(nnzrf), rnp  (nnzrf), rat0,
      write (nfw, 8130)  nnzrf, freq(nnzrf), rnp12(nnzrf), rat0,
     .                   rat1, rat2, rat3
      end do
c
      call giveupus(nscr)
      close (unit = nscr)
c
c     linear extrapolation for values at magnetic axis (rho(npsi) = 0.0)
c
      call extrap (rho(npsi-2),rho(npsi-1),rho(npsi),
     .             qep(npsi-2),qep(npsi-1),qep(npsi))
      call extrap (rho(npsi-2),rho(npsi-1),rho(npsi),
     .             qip(npsi-2),qip(npsi-1),qip(npsi))
      call extrap (rho(npsi-2),rho(npsi-1),rho(npsi),
     .             qcp(npsi-2),qcp(npsi-1),qcp(npsi))
c
c     Also use linear interpolation to put qep,... on transport r-mesh.
c     Re-order rho into xp, and qep,qip,qcp into work1,....
c
c     The next statement was commented out by YRLL and HSJ (10/16/92)
c
****  if (kf .gt. 400)  call STOP ('subroutine FASTWAVE: kf > 400', 975)
c
      do j=1,npsi
        xp(j)    = rho(npsi+1-j)
        work1(j) = qep(npsi+1-j)
        work2(j) = qip(npsi+1-j)
        work3(j) = qcp(npsi+1-j)
      end do
c
      qrfe (1) = work1(1)
      qrfi (1) = work2(1)
      currf(1) = work3(1)
      do j=2,nj
        if (r(j) .le. xp(npsi)) then
          call find(j1,j2,r(j),xp,npsi)
          if (j1 .ne. j2) then
            del = (r(j)-xp(j1))/(xp(j2)-xp(j1))
          else
            del = 0.0
          end if
          qrfe (j) = work1(j1)+(work1(j2)-work1(j1))*del
          qrfi (j) = work2(j1)+(work2(j2)-work2(j1))*del
          currf(j) = work3(j1)+(work3(j2)-work3(j1))*del
        else
          call extrap (xp(npsi-1), xp(npsi), r(j), work1(npsi-1),
     .                 work1(npsi), qrfe (j))
          call extrap (xp(npsi-1), xp(npsi), r(j), work2(npsi-1),
     .                 work2(npsi), qrfi (j))
          call extrap (xp(npsi-1), xp(npsi), r(j), work3(npsi-1),
     .                 work3(npsi), currf(j))
        end if
      end do
c
c --- reset rnp back to the input value
c
      do nnzrf=1,nzrf
        rnp(nnzrf) = rnp12(nnzrf)
      end do
      b0 = sgnb0*b0
c
      return
c
      end

      real*8 function funrfcd(ezrin, eziin, omega, j1, nsp, anz2, ier)
      USE param,only : kj
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c      parameter (ni = 3, ndr = 101)
      parameter (ni = 3, ndr = 2*(kj-1))
c     ndr has to be big enough so that 2(nj-1) .ge. ndr   HSJ 
c
      common /blk2/ x1(ndr),xend(ndr),dnit(ni,ndr),tit(ni,ndr),
     .              tet(ndr),dnet(ndr),b0t(ndr),ze0(ndr),fci(ni,ndr),
     .              vi(ni,ndr),rhoi(ni,ndr),fpi(ni,ndr),amu(ni,ndr),
     .              fce(ndr),fpe(ndr),ve(ndr),fci2(ni,ndr),fpi2(ni,ndr),
     .              zef(ndr)
      common /blk4/ ds,dd,p
      data          c/2.9979e10/
c
      ier = 0
      ds  = 1.0
      dd  = 0.0
      om2 = omega**2
      do 100 n=1,nsp
        ds = ds-fpi2(n,j1)/(om2-fci2(n,j1))
        dd = dd+fpi2(n,j1)*omega/((om2-fci2(n,j1))*fci(n,j1))
  100 continue
      ds   = ds+(fpe(j1)/fce(j1))**2
      p    = (fpe(j1)/omega)**2
      anx2 = (dd**2-(anz2-ds)**2)/(anz2-ds+(dd**2/p)*ezrin*anz2
     .      /(anz2-ds))
      if (anx2 .lt. 0.0)  go to 12
      go to 11
   11 continue
      funrfcd = SQRT (anx2*om2/c**2)
      return
   12 ier = 99
      funrfcd = 0.0
      return
c
****   v00 = SQRT (v0)
****   if (x1(j) .ge. 0.0) then
****   sbr = 2.0 * epsi / (1.0+epsi)
****   else
****   sbr = 0.0
****   end if
****   call fwcd(tempe,dens,epsi,sbr,-cer(j),-cei(j),zef(j),v00,eff(j))
c
c      tempe = temperature of electron in keV
c      dens  = density     of electron in /cm^3.
c
****   current(j) = pabe(j)*eff(j)
c
      end

      subroutine fwcd (tempe, dens, epsil, alt2, alphr, alphi, zeff,
     .                 xe, xjop)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c     A subroutine to calculate the local current density
c     driven by RF waves.
c
c --- Uses CD efficiency empirical formaula based on numerical
c     Fokker-Planck bounce-averaged calculations
c --- Ref.: D.A. Ehst and C.F.F. Karney, ANL/FPP/TM-247(1990).
c
c Use (Comraylh)
c
c Use (Comnew)
c
c Use (Comiorlh)
c
c $$$$$$  added by TK Mau on 11/21/91
      alph0 = 1.0 + alphr*alphr+alphi*alphi
c $$$$$$
      akcd  = 11.91/(0.678+zeff)
      c0cd  = 4.13/zeff**0.707
c
      amcd  = 1.38 + 1.10/alph0
      ccd   = 0.389 - 0.2903/alph0
      acd   = 12.3/alph0
c
      alt1  = SQRT (alt2)
      wte   = 1.414*xe
      cprof = 1.0
      if (alt2 .ne. 0.0) then
        ytt = (1.0-alt2)*wte**2/alt2
        arg = (ccd*ytt)**amcd
        cprof = 1.0 - EXP (-arg)
      end if
      rprof   = 1.0 - (epsil**0.77 * SQRT (12.25+wte**2)) /
     .                (3.5*epsil**0.77+wte)
      amprof  = 1.0 + acd*(alt1/wte)**3
c
c ---------------------------------------------------------- SCC 3/12/92
c
      fact = EXP (-xe*xe)
      awt  = ABS (wte)
      aci  = 4.0*awt**2/(5.0 + zeff)
     .     + 6.0 * (6.0 + zeff)/((5.0 + zeff)*(3.0 + zeff))
     .     + 6.0 * (2.0*alphr+2.-alphi**2)/((5.0 + zeff)*alph0)
     .    + 42.0/(awt*(2.0+zeff)*(1.0+zeff))
      aci2 = 3.759942*(1.0-2.0*alphr/alph0)/zeff
     .    +(16.8465/(0.32922+zeff))*(alphr-1.0)/alph0
     .    +(110.85/(0.62879+zeff))/(4.0 * alph0)
     .    +awt**2*((11.5-1.5/alph0)+(7.5+4.5/alph0)*zeff)/
     .     ((2.0+zeff)*(1.0 + zeff))
     .    -awt**4*8.0 * alphr/(4.0 + zeff)/alph0
     .    +awt**2*55.0 * alphi**2/(45.0 + zeff)/alph0
      aci2 = aci2/awt
      eta0 = fact*aci2+(1.0-fact)*aci
c
c $$$$$$$$$$$$$$$$$$$$$$$$ added by TK Mau on 11/21/91
c
      arg1   = 1.0e-3 / tempe * SQRT (dens)
      algame = 24.0 - LOG (arg1)
      eta    = cprof*amprof*eta0*rprof
c
c $$$$$$$$$$$$$$$$$$$$$$$$
c
****  xjop = 3.81e15/algame*tempe/dens*eta
      xjop = 3.84e13/algame*tempe/dens*eta
c
c $$$$$$$$$$$$$$$$$$$$$$$$
c
****  pwce = -SIGN (1.0,pkpar) * ABS (btor/btot)*btsign*pde*xjop
      return
c
      end

      subroutine fw (nj,nja,r,rfpowr,ene,eni,nspi,te,ti,zeff,iswch,
     .               rmajor,b0,freq,wkz,li,nih,ziir,gmii,nvb,nvba,htslb,
     .               hrf,prfe,prfi,cur,totfwpe,totfwpi,totfwc)
c
      USE param,only : kj
      USE io,only : n66,n77
      implicit  integer (i-n), real*8 (a-h, o-z)

c
c ----------------------------------------------------------------------
c     Apr.28,1988.  upgraded June 23,1988.
c ----------------------------------------------------------------------
c  include ion cyclotron damping and electron Landau and transit time damping
c  trapped particle effect included
c  no mode conversion
c
c  input:
c  iswch = 1--ttmp+landau, for ICRF; iswch=2--landau damping, for LHRF
c  b0 = toroidal magnetic field at the plasma axis  (guass)
c  nsp = number of ion species
c  zii(n) = charge of n-th species
c  gmii(n) = atomic mass number of n-th ion species
c  eni(nj,nsp),ene(nj) are densities at center of each slab  (/cm**3)
c  te(nj),ti(nj) are temperature profiles  (eV)
c  freq = frequency in mega-hertz  (Hz)
c  rmajor = major radius at plasma axis=r0  (cms)
c  wkz = parallel wave-number of RF  (/cm)
c  li(ni) = starting harmonic number for damping calculations
c  nih(ni) = number of harmonics to be calculated
c  nvb = number of blocks in vertical direction
c  htslb(j),j = 1,nj: height of each block in jth slab  (cm)
c  hrf(j,2) delimits vertical extent of RF in each slab;  (cm)
c                 (j,1) the upper bound;(j,2) the lower bound;
c  rfpowr = rfpower (watts)
c   r(nj) = major radial coordinates of centers of slabs  (cm)
c   nj = number of slabs on equatorial plane of tokamak
c  zeff(nj) = z-effective at nj
c
c  output:
c  prfi(nsb,nvb,nspi) = rf-power absorbed by each ion species in each block
c  prfe(nsb,nvb) = rf-power absorbed by electrons in each block
c    (watts/cm**3)
c  cur(nsb,nvb) = current density driven by RF in each block
c    (amps/cm**2)
c
c  ndr = number of slabs
c  ni  = number of ion species = nsp
c  il  = number of ion harmonics
c ----------------------------------------------------------------------
c
c      parameter (ndr = 101, ni = 3, ne = 1, il = 35, ndv = 101)
      parameter (ndr = 2*(kj-1), ni = 3, ne = 1, il = 35, ndv = 101)
c     ndr has to be big enough so that 2(nj-1) .ge. ndr   HSJ 
      dimension r(nja),ene(nja),eni(nja,nspi),te(nja),ti(nja),
     .          zeff(nja),li(nspi),nih(nspi),
     .          ziir(nja,nspi),gmii(nspi),
     .          cur(nja,nvba),prfe(nja,nvba),prfi(nja,nvba,nspi),
     .          htslb(nja),hrf(nja,2)
      dimension xv(ndv),vol(ndr),pai(ni)
c
      common /blk1/ zi(ni,ndr), gmi(ni), bt0, wkpl, r0, rfpow, frq,
     .              area(ndr),
     .              nsp, nsb, li0(ni), nli(ni)
c
      common /blk2/x1(ndr),xend(ndr),dnit(ni,ndr),tit(ni,ndr),
     .            tet(ndr),dnet(ndr),b0t(ndr),ze0(ndr),fci(ni,ndr),
     .            vi(ni,ndr),rhoi(ni,ndr),fpi(ni,ndr),amu(ni,ndr),
     .            fce(ndr),fpe(ndr),ve(ndr),fci2(ni,ndr),fpi2(ni,ndr),
     .            zef(ndr)
c
c ----------------------------------------------------------------------
c blk2 --- profiles
c ----------------------------------------------------------------------
c
      real*8                           kx
      common /blk3/ wok(ndr),wok1(ndr),kx(ndr),pabi(ni,ndr),
     .           pabe(ndr),current(ndr),pabsi(ndr),pabse(ndr),totcur_fw,
     .              bessi(il),aiki(ni,ndr),tik(ndr),denom(ndr),aike(ndr)
      common /blk4/ ds,dd,p
      common /blk5/ wk1(il,ni),wk2(il,ni),wk3(il,ni)
      common /blk6/ eff(ndr)
      data          pi/3.14159265/
c
      ierr = 0
      if (ni .lt. nspi .or. ndr .lt. nj .or. ndv .lt. nvb)  ierr = 99
      if(ierr .eq. 99)print *,'ni,nspi,ndr,nj,ndv,nvb =',
     .                         ni,nspi,ndr,nj,ndv,nvb
c
      do n=1,nspi
        if (il .lt. nih(n))  ierr = 99
        print *,'il,nih(n) =',il,nih(n)
      end do
c
      if (ierr .eq. 99) then
      write  (n66, 1111)
 1111 format (' WARNING: you should set ni .ge. nspi, ndr .ge. nj,'
     .        ' il .ge. nih(n) and ndv .ge. nvb' /
     .        ' Subroutines where these parameter statements',
     .        ' occur are FW, RFCD, PARAMFW and FUN.')
      call giveupus(n66)
      call STOP ('subroutine FW: unspecified problem', 80)
      end if
c
      call ell
      nsp = nspi
      frq = 2.0 * pi * freq
      r0 = rmajor
      bt0 = b0
      rfpow = rfpowr
      wkpl = wkz
      nsb = nj
      do 50 n=1,nsp
        li0(n) = li(n)
        nli(n) = nih(n)
   50   gmi(n) = gmii(n)
      do i=1,nj
        x1(i) = r(i)-r0
        do n=1,nsp
          zi(n,i)   = ziir(i,n)
          dnit(n,i) = eni(i,n)
          tit(n,i)  = ti(i)
          if (tit(n,i) .le. 0.0)  tit(n,i) = 1.0
        end do
        dnet(i) = ene(i)
        zef(i) = zeff(i)
        tet(i) = te(i)
        if (tet(i) .le. 0.0) tet(i) = 1.0
      end do
      do i=1,nj
        rmaj = r0+x1(i)
        area(i) = (hrf(i,1)-hrf(i,2))*2.0 * pi * rmaj
      end do
      call rfcd (iswch)
c
      nv = (nvb-1)/2
      nsbt = nsb
      do 300 i=1,nj
      hhtslb = 0.5*htslb(i)
c
      if (i .eq. 1 .or. i .eq. nsbt) then
        if (i .eq. 1) then
          wslb = x1(1)-x1(2)
        else
          wslb = x1(nsbt-1)-x1(nsbt)
        end if
      else
        wslb = 0.5*(x1(i-1)-x1(i+1))
      end if
c
      vol(i) = area(i)*wslb
      area(i) = wslb*(hrf(i,1)-hrf(i,2))
      pe = pabe(i)/vol(i)
      do 310 n=1,nsp
  310 pai(n) = pabi(n,i)/vol(i)
      curr   = current(i)/area(i)
      do 400 j=1,nvb
      xv(j)  = htslb(i)*(j-nv-1)
      if ((xv(j)+hhtslb) .le. hrf(i,2)  .or.
     .    (xv(j)-hhtslb) .ge. hrf(i,1)) then
      prfe(i,j) = 0.0
      do 410 n=1,nsp
        prfi(i,j,n) = 0.0
  410 continue
      cur(i,j) = 0.0
      else
      x3 = xv(j)+hhtslb
      x2 = xv(j)-hhtslb
      if (x3 .ge. hrf(i,1)) x3 = hrf(i,1)
      if (x2 .le. hrf(i,2)) x2 = hrf(i,2)
      fac = (x3-x2)/htslb(i)
      do 420 n=1,nsp
  420 prfi(i,j,n) = pai(n)*fac
      prfe(i,j) = pe*fac
      cur(i,j) = curr*fac
      end if
  400 continue
  300 continue
      totfwpe = pabse(nsbt)
      totfwpi = pabsi(nsbt)
      totfwc = totcur_fw
      return
c
      end

      real*8 function gausin (a, b)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      external gausm
c
      call intgrl (a, b, gausm, sum, 50)
      gausin = sum
      return
c
      end

      real*8 function gausm (x)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      gausm = 0.8862269 * EXP (-x**2)
      return
c
      end

      subroutine interfw (f, z, y1, y2, y3, x1, x2, x3)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c --- interfw interpolates between 3-pts. xi and find f at z
c
      a1 = z-x1
      a2 = z-x2
      a3 = z-x3
      b1 = x1-x2
      b2 = x1-x3
      b3 = x2-x3
      f  = a1*a2*y3/(b2*b3)
      f  = f-a1*a3*y2/(b1*b3)
      f  = f+a2*a3*y1/(b1*b2)
      return
c
      end

      subroutine intgrl (a, b, f, simps, n)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      external F
c
      an     = n
      twoh   = (b-a)/an
      h      = twoh/2.0
      sumend = 0.0
      summid = 0.0
      do k=1,n
        ck1    = k-1.0
        x      = a+ck1*twoh
        sumend = sumend+f(x)
        summid = summid+f(x+h)
      end do
      simps = (2.0 * sumend + 4.0 * summid - f(a) + f(b)) * h / 3.0
      return
c
      end

      subroutine paramfw
      USE param,only : kj
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c  param calculates the plasma parameters necessary for evaluation
c  of the dispersion relation
c
c  dni(i,j)  = ion density of ith species at jth slab
c  dne(j)    = electron density
c  zi(i,j)   = charge state of ith species at the jth slab
c  gmi(i)    = mass number of ith species
c  nsp       = number of species
c  ti(i,j)   = temperature of ith species at jth slab
c  te(j)     = electron temperature
c  frq       = frequency,  bot(j) = toroidal b-field in jth slab
c  fci(i,j)  = ion cyclotron frequency of ith species
c  fce(j)    = electron cyclotron frequency
c  vi(i,j)   = ith species thermal velocity
c  ve(j)     = electron thermal velocity
c  rhoi(i,j) = gyroradius of ith species
c  fpi(i,j)  = plasma frequency of ith ion species
c  fpe(j)    = electron plasma frequency
c  nslb      = number of slabs
c
c      parameter    (ndr = 101, ni = 3, ne = 1, il = 35)
      parameter    (ndr = 2*(kj-1), ni = 3, ne = 1, il = 35)
c     ndr has to be big enough so that 2(nj-1) .ge. ndr   HSJ 
      real*8        kx
c
      common /blk1/ zi(ni,ndr), gmi(ni), bt0, wkpl, r0, rfpow, frq,
     .              area(ndr),
     .              nsp, nsb, li0(ni), nli(ni)
c
c ----------------------------------------------------------------------
c  blk1---input
c ----------------------------------------------------------------------
c
      common /blk2/x1(ndr),xend(ndr),dnit(ni,ndr),tit(ni,ndr),
     .   tet(ndr),dnet(ndr),b0t(ndr),ze0(ndr),fci(ni,ndr),
     .   vi(ni,ndr),rhoi(ni,ndr),fpi(ni,ndr),amu(ni,ndr),
     .   fce(ndr),fpe(ndr),ve(ndr),fci2(ni,ndr),fpi2(ni,ndr),zef(ndr)
c
c ----------------------------------------------------------------------
c  blk2---profiles
c ----------------------------------------------------------------------
c
      common /blk3/wok(ndr),wok1(ndr),kx(ndr),pabi(ni,ndr),
     .          pabe(ndr),current(ndr),pabsi(ndr),pabse(ndr),totcur_fw,
     .          bessi(il),aiki(ni,ndr),tik(ndr),denom(ndr),aike(ndr)
      common /blk4/ds,dd,p
      common /blk5/wk1(il,ni),wk2(il,ni),wk3(il,ni)
      common /blk6/eff(ndr)
c
      data fec/1.76e7/, fic/9.58e3/, fip/1.32e3/, fep/5.64e4/,
     .     fte/5.9255e7/, fti/1.3845e6/
c
      do j=1,nsb
        do i=1,nsp
          fpi (i,j) = fip * zi(i,j) * SQRT (dnit(i,j)/gmi(i))
          fpi2(i,j) = fpi(i,j)**2
          fci (i,j) = fic*zi(i,j)*b0t(j)/gmi(i)
          fci2(i,j) = fci(i,j)**2
          vi  (i,j) = fti * SQRT (tit(i,j)/gmi(i))
          rhoi(i,j) = vi(i,j)/fci(i,j)
        end do
        fce(j) = fec *       b0t (j)
        ve (j) = fte * SQRT (tet (j))
        fpe(j) = fep * SQRT (dnet(j))
      end do
c
      return
c
      end

      subroutine rfcd (iswch)
      USE param,only :kj
      USE io,                     ONLY  : n66,n77
      USE replace_imsl,           ONLY  : my_ mmbsir
c
      implicit  integer (i-n), real*8 (a-h, o-z)

c
c      parameter (ndr = 101, ni = 3, ne = 1, il = 35)
       parameter (ndr = 2*(kj-1), ni = 3, ne = 1, il = 35)
c     ndr has to be big enough so that 2(nj-1) .ge. ndr   HSJ 
      real*8     kx
      complex*16 ai, cz0, cz1
c
      common /blk1/ zi(ni,ndr), gmi(ni), bt0, wkpl, r0, rfpow, frq,
     .              area(ndr),
     .              nsp, nsb, li0(ni), nli(ni)
c
c ----------------------------------------------------------------------
c  blk1---input
c ----------------------------------------------------------------------
c
      common /blk2/x1(ndr),xend(ndr),dnit(ni,ndr),tit(ni,ndr),
     .            tet(ndr),dnet(ndr),b0t(ndr),ze0(ndr),fci(ni,ndr),
     .            vi(ni,ndr),rhoi(ni,ndr),fpi(ni,ndr),amu(ni,ndr),
     .            fce(ndr),fpe(ndr),ve(ndr),fci2(ni,ndr),fpi2(ni,ndr),
     .            zef(ndr)
c
c ----------------------------------------------------------------------
c  blk2---profiles
c ----------------------------------------------------------------------
c
      common /blk3/wok(ndr),wok1(ndr),kx(ndr),pabi(ni,ndr),
     .         pabe(ndr),current(ndr),pabsi(ndr),pabse(ndr),totcur_fw,
     .            bessi(il),aiki(ni,ndr),tik(ndr),denom(ndr),aike(ndr)
      common /blk4/ds,dd,p
      common /blk5/wk1(il,ni),wk2(il,ni),wk3(il,ni)
      common /blk6/eff(ndr)
c
      dimension work(50)
      dimension ierr(ndr)
      dimension wk4(il,ni)
      dimension cer(ndr),cei(ndr)
      real*8    bes_order      ! HSJ
      data      c/2.9979e10/, pi/3.14159265/
c
c   input:
c   nsb = original number of slabs
c   deni(ni,l) = ion density of species ni at l-th slab
c   zi(ni,l) = charge state of ni-th species at the l-th slab
c   gmi(ni) = mass number of ni-th species
c   ti(ni,l) = temperture of ni-th species at l-th slab
c   te(l) = electron temperture
c   frq = frequency*2*pi
c   bt0 = magnetic field at plasma axis
c   r0 = major radius
c   xl0 = coordinate of center of first slab (low field side) measured from
c       the plasma axis
c   nsp = number of ion species
c   wslb = width of original slabs
c        note:slabs of low field and high field parts will also have
c             widths wslb, but slabs in center part will be further
c             subdivided to m*(n2-n1+1) slabs
c   li = number of      ion harmonics
c   le = number of electron harmonics
c   wkpl = parallel wave number
c   nl,nh = slab number (nh .gt. nl) delineating the center part of plasma
c         (will use mode conversion code for this part).
c   pflxin = actual input rfpow flux
c    m = number into which the original slabs are subdivided in the central
c   x(l) = coordinate of center of l-th original slab as measured
c        from the plasma axis
c
c  can output:
c  flxabp = rfpow going into the plasma
c  flxrf is rf-flux and should be normalized to pflxin
c  pabi(ni,j) = integrated rfpow absorbed by ni-th species up to j-th slab
c  pabe(j) = local rfpow absorbed by electron
c  pabsi(j) = integrated rfpow absorbed by ions up to j-th slab
c  pabse(j) = integrated rfpow absorbed by electrons up to j-th slab
c   pabs = total rfpow absorbed per unit area
c   flxabs = total flux incident from the low field side
c     and equals flxin-reflected flux
c ----------------------------------------------------------------------
c
      bes_order = 0.0
      alen = r0 / 100.0
      nsbt = nsb
      a0 = ABS (x1(1))
      do l=1,nsbt
        b0t(l) = bt0*r0/(r0+x1(l))
      end do
      ai = (0.0, 1.0)
      call paramfw
      rtpi = 1.7724539
      do 200 j=1,nsbt
      wkpl1 = wkpl*(r0+a0)/(r0+x1(j))
      vph = frq/wkpl1
      anz = c/vph
      anz2 = anz**2
c
      if (j .eq. 1 .or. j .eq. nsbt) then
        if (j .eq. 1) then
          wslb = x1(1)-x1(2)
        else
          wslb = x1(nsbt-1)-x1(nsbt)
        end if
      else
        wslb = 0.5*(x1(j-1)-x1(j+1))
      end if
c
      ze0(j) = vph/ve(j)
      call zeta(ze0(j),cz0,cz1)
      zabs  =    ABS (cz1)*ze0(j)**2
      ezrin =  -REAL (cz1)*ze0(j)**2/zabs**2
      eziin = -PIMAG (cz1)*ze0(j)**2/zabs**2
      kx(j) = funrfcd(ezrin,eziin,frq,j,nsp,anz2,ierr(j))
      if (ierr(j) .eq. 99)  go to 200
      denom(j) = anz2-ds+(dd**2/p)*ezrin*(anz2/(anz2-ds))
      work(1)  = dd**2*eziin*anz2/p
****  work(1)  = (1.0+dd**2*ezrin/(anz2-ds))*dd**2*eziin*anz2/p
c
c This expression is slightly more accurate, but QL operator has to change also
c
      work(2) = 0.5*rtpi*((fpe(j)/c)*ve(j)/fce(j))**2
     .             *ze0(j)* EXP (-ze0(j)**2)*(anz2-ds)
      aike(j) = kx(j)*(work(1)/denom(j)+work(2))
     .                        /(denom(j)+ds*dd**2*ezrin/(p*(anz2-ds)))
      work(1) = 2.0 * (c/ve(j))**2*(fce(j)/frq)*dd/(ds-anz2)
      cer(j) = work(1)*ezrin/p
      cei(j) = work(1)*eziin/p
      do n=1,nsp
        amu(n,j) = 0.5*(kx(j)*vi(n,j)/fci(n,j))**2
      end do
      wok(j) = 0.0
      tik(j) = 0.0
      do n=1,nsp
      li0(n) = 0
      l1 = li0(n)
      l2 = li0(n)+nli(n)-1
      gl = l1/10.0
      aiki(n,j) = 0.0
c
      if (0.25*amu(n,j)**2 .gt. gl) then
        call my_ mmbsir (amu(n,j), bes_order, l2+1, 2, work, ier)
        noutloc = n66
        if (ier .gt. 0)  write (noutloc, 2030) j, l, ier
        do 222 l=l1+1,l2+1
          l3 = l-l1
  222     bessi(l3) = work(l)
      else
        do l=l1,l2
          l3 = l-l1+1
          sum = 1.0
          do 230 m=1,l
  230       sum = sum*0.5*amu(n,j)/m
          bessi(l3) =
     .      EXP (-amu(n,j)) * sum * (1.0+(0.5*amu(n,j))**2/(l+1.0))
        end do
      end if
c
      do l=1,nli(n)
        wk1(l,n) = (frq-(l1+l-1)*fci(n,j))/(wkpl1*vi(n,j))
        n1 = l1+l-1
        if (n1 .ne. 0) then
          c1 = n1*fci(n,j)*wslb*0.5/(wkpl1*vi(n,j)*(r0+x1(j)))
          c2 = wk1(l,n)+c1
          c3 = wk1(l,n)-c1
          wk2(l,n) =
     .      (0.5 * pi * frq*(r0+x1(j))/(n1*fci(n,j)*wslb))*gausin(c3,c2)
        else
          wk2(l,n) = rtpi*wk1(l,n) * EXP (-wk1(l,n)**2)
        end if
      end do
c
      wok(j) = (fpi(n,j)/frq)**2
      anx2 = (c*kx(j)/frq)**2
      do l=1,nli(n)-1
        wk1(l,n) = (l1+l-1.0)**2*bessi(l)/amu(n,j)
        wk1(l,n) = wk1(l,n)*wk2(l,n)*wok(j)
        if (l .eq. 1) then
          wbsp = 2.0 * bessi(l+1)
        else
          wbsp = bessi(l-1)+bessi(l+1)
        end if
        wk4(l,n) = amu(n,j)*(bessi(l)-0.5*wbsp)
        wk4(l,n) = wk4(l,n)*wk2(l,n)*wok(j)
        wk1(l,n) = wk1(l,n)+wk4(l,n)
        wk3(l,n) = (l1+l-1.0)*(bessi(l)-0.5*wbsp)
     .           * wk2(l,n)*wok(j)
        aiki(n,j) = aiki(n,j)+(((ds-anz2)*wk1(l,n)-dd*wk3(l,n))
     .     *denom(j)/((ds-anz2)**2-dd**2)+0.5*((wk1(l,n)-wk4(l,n))
     .     -(2.0*dd*wk3(l,n)-2.0*ds*wk1(l,n)+(anx2+anz2)*(wk1(l,n)
     .     -wk4(l,n)))*ezrin/p))/(denom(j)+ds*dd**2*ezrin/(p*(anz2-ds)))
      end do
c
      aiki(n,j) = 2.0 * aiki(n,j)*kx(j)
      tik(j) = tik(j)+aiki(n,j)
      end do
      tik(j) = tik(j)+aike(j)
  200 continue
c
 2030 format ('  j = ',i5, '  l=',i5, '  ier=',i5)
c
c ----------------------------------------------------------------------
c     calculate rfpow absorption
c ----------------------------------------------------------------------
c
      totcur_fw = 0.0
      isw = 0
      do j=1,nsbt
        wok (j) = 0.0
        wok1(j) = 0.0
      end do
      do j=1,nsbt
        iflag = 0
        if (j .eq. 1 .or. j .eq. nsbt) then
          if (j .eq. 1) then
            wslb = x1(1)-x1(2)
          else
            wslb = x1(nsbt-1)-x1(nsbt)
          end if
        else
          wslb = 0.5*(x1(j-1)-x1(j+1))
        end if
        if (ierr(j) .eq. 99) tik(j) = 0.0
        wok(j) = tik(j)*wslb
        if (j .eq. 1)  go to 310
        wok(j) = wok(j)+wok(j-1)
  310   wok1(j) = EXP (-wok(j))
        if (j .eq. 1) then
          pabe(j) = 1.0 - wok1(1)
        else
          pabe(j) = wok1(j-1)-wok1(j)
        end if
        pabsi(j) = 0.0
        do 320 n=1,nsp
        if (ierr(j) .eq. 99 .or. tik(j) .le. 1.0e-20) then
          pabi(n,j) = 0.0
          go to 320
        end if
        pabi(n,j) = pabe(j)*aiki(n,j)/tik(j)
        pabi(n,j) = rfpow*pabi(n,j)
  320   pabsi(j) = pabsi(j)+pabi(n,j)
        if (ierr(j) .eq. 99 .or. tik(j) .le. 1.0e-20)  go to 330
        pabe(j) = pabe(j)*rfpow*aike(j)/tik(j)
  330   if (ierr(j) .eq. 99 .or. tik(j) .le. 1.0e-20) then
          pabe(j) = 0.0
          pabsi(j) = 0.0
          isw = 1
        end if
        if (pabe(j)/rfpow .le. 1.0e-8) iflag = 1
        if (j .eq. 1) then
          pabse(1) = pabe (1  )
        else
          pabse(j) = pabse(j-1) + pabe (j  )
          pabsi(j) = pabsi(j  ) + pabsi(j-1)
        end if
        epsi = ABS (x1(j)/r0)
        if (epsi .le. 1.0e-3)  epsi = 1.0e-3
        wkpl1 = wkpl*(r0+a0)/(r0+x1(j))
        v0 = (frq/(wkpl1*ve(j)))**2
        v00 = SQRT (v0)
        if (x1(j) .ge. 0.0) then
          sbr = 2.0 * epsi/(1.0 + epsi)
        else
          sbr = 0.0
        end if
        if (iflag .eq. 0) then
          tempe = tet(j)*1.0e-3    ! electron temperature in keV
          dens  = dnet(j)          ! electron density per cc
          call fwcd (tempe, dens, epsi, sbr, -cer(j), -cei(j),
     .               zef(j), v00, eff(j))
c
c         tempe = electron temperature in keV
c         dens  = electron density per cc.
c
****      alglam = 15.0
****      dnej   = dnet(j)*1.0e-14
****      tej    = tet(j)*1.0e-4
****      eff(j) = eff(j)*1.21*tej/(alglam*dnej*alen)
c
          eff(j) = eff(j)/(6.283185*alen)    ! eff/(2.0 * pi * r0)
          if (epsi .le. 1.0e-3) then
            call interfw (eff(j),x1(j),eff(j-1),eff(j-2),eff(j-3),
     .                    x1(j-1),x1(j-2),x1(j-3))
          end if
        else
          eff(j) = 0.0
        end if
        current(j) = pabe(j)*eff(j)
        totcur_fw     = totcur_fw + current(j)
      end do
      return
c
      end

      subroutine simson (a, b, F, simps, n)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      external F
c
      an     = n
      twoh   = (b-a) / an
      h      =  twoh / 2.0
      sumend = 0.0
      summid = 0.0
      do k=1,n
        ck1    = k - 1.0
        x      = a + ck1 * twoh
        sumend = sumend + f(x)
        summid = summid + f(x+h)
      end do
      simps = (2.0 * sumend + 4.0 * summid - f(a) + f(b)) * h / 3.0
      return
c
      end

      subroutine zeta (x, cz0, cz1)
c
c     This subroutine evaluates the plasma dispersion function
c     and its derivative for real argument x.
c     Originally developed at PPPL.
c
      implicit  integer (i-n), real*8 (a-b, d-h, o-z), complex*16 (c)
c
      data zsqpi /1.7724538509/, error /1.0e-7/
c
      x2 = x * x
      if (ABS (x) .gt. 5.7)  go to 40
c
c     power series
c
      hold  = EXP (-x2)
      sumzr = 1.0e0
      terme = 1.0e0
      fn    = 0.0
   10 fn    = fn + 1.0e0
      terme = terme*x2/fn
      term  = terme/(fn+fn+1.0e0)
      sumzr = sumzr + term
      if (ABS (term/sumzr) .gt. error)  go to 10
      cz0 = hold * CMPLX (-2.0 * x * sumzr, zsqpi)
      cz1 =               -2.0 * (1.0 + x*cz0)
      return
c
c     asymptotic series
c
   40 fn    = 1.0e0
      fact  = 1.0e0
      termz = 1.0e0
      dzr   = 1.0e0
      fmult = 0.5e0/x2
   30 fn    = fn+2.0e0
      fact  = fn*fmult
      if (ABS (fact) .gt. 1.0e0)  go to 20
      termz = termz*fact
      dzr   = dzr + termz
      if (ABS (termz) .gt. error)  go to 30
c
   20 cz1 = dzr / x2
      cz0 = -(1.0 + 0.5*cz1) / x
      return
c
      end
