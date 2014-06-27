
      real*8 function epslon (x)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      real*8  x, a, b, c, eps
c
c     estimate unit roundoff in quantities of size x.
c
c     this program should function properly on all systems
c     satisfying the following two assumptions,
c        1.  the base used in representing floating point
c            numbers is not a power of three.
c        2.  the quantity  a  in statement 10 is represented to
c            the accuracy used in floating point variables
c            that are stored in memory.
c     the statement number 10 and the go to 10 are intended to
c     force optimizing compilers to generate code satisfying
c     assumption 2.
c     under these assumptions, it should be true that,
c            a  is not exactly equal to four-thirds,
c            b  has a zero for its last bit or digit,
c            c  is not exactly equal to one,
c            eps  measures the separation of 1.0 from
c                 the next larger floating point number.
c     the developers of EISPACK would appreciate being informed
c     about any systems where these assumptions do not hold.
c
c     this version dated 4/6/83.
c
      a   = 4.0 / 3.0
   10 b   = a - 1.0
      c   = b + b + b
      eps = ABS (c-1.0)
      if (eps .eq. 0.0)  go to 10
      epslon = eps * ABS (x)
      return
c
      end

      subroutine etawn8_12 (letain_wn6, letain_wn7, cetain_wn30,
     .                      cetain_wn32, lprintin, neq, epsnhin,
     .                      epsnzin, epsnsin, epstein, epsthin,
     .                      epstzin, tauhin, tauzin, fnzin, czin, azin,
     .                      fnsin, ftrapein, ekyrhoin, ekparlin, ndim,
     .                      iounit,weiland_output,omega, gamma, difthi,
     .                      velthi, chieff, nmodes, perform)
c
c ... This is a modified version of the ETAWN8 subroutine as provided
c ... by G. Bateman. Some changes have been made to incorporate the
c ... routine into the ONETWO transport code (thus the "_12" appended
c ... to the name).
c
c ... Changes:  Joop Konings
c
c ... Header: - Arrays letain, cetain changed to the necessary
c ...           parameters letain_wn6, letain_wn7 and
c ...           cetain_wn32, cetain_wn30 respectively.
c ...         - Removed: betahin, betazin
c ... For the eigenvalue solver either IMSL, or NAG or EISPACK can be
c ... used.
c
c ... This version of ETAWN8 is intended for use on NERSC CRAY computers
c ... however, use of the EISPACK eigenvalue solver (subroutine RGG)
c ... makes it a stand-alone subroutine.
c
c ... The following external subprograms need to be provided:
c
c ... GVCRG  IMSL subroutine to compute eigenvalues and eigenvectors
c ... GPIRG  IMSL function   to compute the performance index
c ... F02BJE NAG  subroutine to compute eigenvalues and eigenvectors
c ... F01AAE NAG  subroutine to compute inverse of matrix
c ... F02AKE NAG  subroutine to compute complex eigenvalues of a matrix.
c ... STOP        subroutine to terminate, issue message, set exit value
c
c ... Note that the EISPACK routines RGG etc. are provided with this
c ... source code.
c
c ... JAK All I/O can be excluded for lprintin < 0 (to speed up)
c
      implicit none
c
      logical       inital         , lmatv
      data          inital /.true./, lmatv /.false./
      character*(*) weiland_output
c
      integer    idp
      parameter (idp = 9)
c
      integer    letain(32), letain_wn6, letain_wn7,
     .           lprintin, lprint, neq, idim_l, ndim, nmodes,
     .           ieq, imatrx, j1, j2, j, jd, i, k, iounit
c
      real*8     PIMAG, zero, one,
     .           cetain(32), cetain_wn30, cetain_wn32,
     .           epsnhin, epsnzin, epsnsin,
     .           epstein, epsthin, epstzin, tauhin, tauzin,
     .           fnzin, czin, azin, fnsin, ftrapein,
     .           ekyrhoin, ekparlin, omega(*), gamma(*),
     .           difthi(ndim,*), velthi(*), chieff(*), perform(*)
c
      real*8     betahin      , betazin          ! obsolete?
      data       betahin /0.0/, betazin /0.0/    ! obsolete?
c
c ndim   = first dimension of the 2-D array difthi
c          and the maximum number of unstable modes allowed
c nmodes = number of unstable modes
c ieq    = number of equations
c imatrx = the number of elements computed along each row and
c          column of the transport matrix
c Note: Normally the transport matrix is ieq-1 by ieq-1.
c       However, if there are 6 equations, compute only a 4 by 4 matrix
c       since the impurity temperature equation is not used by most
c       transport codes including BALDUR
c
c ... added for EISPACK routine RG
c
      integer matz, ierr
      real*8  dzam(idp,idp), dzbm(idp,idp), dzalfr(idp), dzalfi(idp),
     .        dzbeta(idp), dzeigenv(idp,idp)
c
      integer iter(idp), ifail
c
      real*8  zam(idp,idp), zbm(idp,idp), zamt(idp,idp), zbmt(idp,idp),
     .        zalfr(idp), zalfi(idp), zbeta(idp), zeigenv(idp,idp),
     .        zomega(idp), zgamma(idp), ztol
c
c zam(i,j) = matrix A
c zbm(i,j) = matrix B
c   Note that the eigenvalues are
c zomega(j) = zalfr(j) / zbeta(j)
c zgamma(j) = zalfi(j) / zbeta(j)
c where beta(j) will be 0.0 in the case of an infinite eigenvalue
c zeigenv(j) = eigenvector
c
****  real*8      gpirg
      complex*16  zalpha(idp), zevec(idp,idp)
c
c zalpha = (zalfr(j), zalfi(j))  from IMSL routine GVCRG
c zevec  = real + imaginary part of eigenvectors
c
c ... NAG; to be sure about dimensions of matrices usews parameter idp
c
      integer    iar      , iai      , ivr      , ivi
      parameter (iar = idp, iai = idp, ivr = idp, ivi = idp)
****  integer    intger(idp)
      real*8     ar(iar,idp), ai(iar,idp), rr(idp), ri(idp),
     .           vr(iar,idp), vi(iar,idp)
      real*8     bori(idp,idp), binv(idp,idp)
cjak  real*8     wkspce(neq)
c
      real*8 zepsmach, zepsqrt,
     .       zetah, zetaz, zetae,
     .       zepsne, zepsnh, zepsnz, zepsns, zepste, zepsth, zepstz,
     .       ztauh, ztauz, zep2nh, zep2nz, zep2ne, zft,
     .       zimp, zfnz, zmass, zfns, zflh, zflz, zgamax
c
      real*8 zgrdth, zgrdte, zgrdtz, zgrdnh, zgrdne, zgrdnz, zgrdns,
     .       zgrdph, zgrdpe, zgrdpz
c
      real*8 zflxph, zflxnh, zflxpe, zflxnz, zflxpz,
     .       ztemp1, ztemp2, zreal, zimag,
     .       ztempa(idp), ztempb(idp)
c
c  These are the effective thermal and particle diffusivities
c  computed in different ways for comparison.
c
      real*8 zflxm(idp,idp), zchim(idp,idp), zgm(idp,idp), zdg
c
c  zflxm(j,jd) are the normalized transport fluxes of:
c  zflxm(1,jd)   n_H T_H
c  zflxm(2,jd)   n_H
c  zflxm(3,jd)   n_e T_e
c  zflxm(4,jd)   n_Z
c  zflxm(5,jd)   n_Z T_Z
c
c  zgm(j,jd)   is a matrix of gradient scale lengths
c  zgm(1,jd) = zgrdth
c  zgm(2,jd) = zgrdnh
c  zgm(3,jd) = zgrdte
c  zgm(4,jd) = zgrdnz
c
c  Here jd=1 uses the input values zgrdph, zgrdnh, zgrdpe, zgrdnz
c  jd=2 uses zgrdph + zdg, zgrdnh, zgrdpe, zgrdnz
c  jd=3 uses zgrdph, zgrdnh + zdg, zgrdpe, zgrdnz
c  jd=4 uses zgrdph, zgrdnh, zgrdpe + zdg, zgrdnz
c  jd=5 uses zgrdph, zgrdnh, zgrdpe, zgrdnz + zdg
c  That is, the gradients are imcemented one at a time.
c
      save idim_l, zepsmach, zepsqrt, inital
c
c ... JAK 4/11/95 fill local arrays letain and cetain
c ... No call to NAG routine eigenvalue solver. Use IMSL
c
      zero = 0.0
      one  = 1.0
c
      do i=1,32
        letain(i) = 0
        cetain(i) = 0.0
      end do
c
      letain( 6) = letain_wn6
      letain( 7) = letain_wn7
      cetain(30) = cetain_wn30
      cetain(32) = cetain_wn32
c
      lprint = lprintin
c
c ... JAK 3/17/95 ... output to file
c
      if (lprint .ge. 0 .and. iounit .gt. 0) then
        call getioun(iounit,iounit)
        open (unit = iounit, file = weiland_output, status = 'UNKNOWN')
      end if
c
c ... initialize variables
c
      ieq = MAX (2, neq)
c
      if (inital) then
c
        idim_l = idp
c
        zepsmach = 0.5
  2     if (0.5 * zepsmach + 1.0 .gt. 1.0) then
          zepsmach = 0.5 * zepsmach
          go to 2
        end if
c
        zepsqrt = SQRT (zepsmach)
c
c ... Print Header
c
        if (lprint .ge. 0) then
c
        if (iounit .gt. 0)  write (iounit, *)
        if (iounit .gt. 0)  write (iounit, *)
     .    'Weiland-Nordman eigenvalue equations, subroutine ETAWN8_12'
        if (iounit .gt. 0)  write (iounit, *) zepsqrt
c
        if      (letain(6) .eq. 2) then
          if (iounit .gt. 0)  write (iounit, *)
     .              ' Eigenvalues computed using IMSL routine gvcrg'
        else if (letain(6) .eq. 1) then
          if (iounit .gt. 0)  write (iounit, *)
     .              ' Eigenvalues computed using NAG14 routine F02BJE'
        else if (letain(6) .eq. 0) then
          if (iounit .gt. 0)  write (iounit, *)
     .              ' Eigenvalues computed using NAG14 routine F02AKE'
        else if (letain(6) .eq. 3) then
          if (iounit .gt. 0)  write (iounit, *)
     .           ' Eigenvalues computed using EISPACK'
        else
          if (iounit .gt. 0)  write (iounit, *)
     .           ' Illegal value letain(6) = ',letain(6)
        end if
c
        if (ieq .eq. 6) then
          if (iounit .gt. 0)  write (iounit, *)
          if (iounit .gt. 0)  write (iounit, *)
     .     'Six eqns for e phi/T_e, T_H, n_H, T_{et}, n_Z, and T_Z'
        else if (ieq .eq. 4) then
          if (iounit .gt. 0)  write (iounit, *)
          if (iounit .gt. 0)  write (iounit, *)
     .      ' Four eqns for e phi/T_e, T_H, n_H, and T_e'
        else if (ieq .eq. 2) then
          if (iounit .gt. 0)  write (iounit, *)
          if (iounit .gt. 0)  write (iounit, *)
     .              ' Two eqns for e phi/T_e and T_H'
        else
          if (iounit .gt. 0)  write (iounit, *)
          if (iounit .gt. 0)  write (iounit, *)
     .               ieq,' = ieq in subroutine ETAWN8_12'
          call giveupus(iounit)
          close (iounit)
          call STOP ('subroutine ETAWN8_12: wrong value #1 of IEQ', 206)
        end if
c
        end if    ! lprint < 0
c
        inital = .false.
c
      end if
c
c ... end of initialization
c
      ieq = MAX (2, neq)
c
      imatrx = MIN (ieq - 1, 4)
c
c ... print header
c
      if (lprint .gt. 0 .and. iounit .gt. 0) then
              write (iounit, *)
              write (iounit, *)
     .    'Weiland-Nordman eigenvalue equations, subroutine ETAWN8_12'
              write (iounit, *)
     .                '(all frequencies normalized by omega_{De})'
              write (iounit, *) '(all diffusivities normalized by '
     .                                      ,'omega_{De} / k_y^2)'
              write (iounit,108) 'epsnhin' , epsnhin
              write (iounit,108) 'epsnzin' , epsnzin
              write (iounit,108) 'epsnsin' , epsnsin
              write (iounit,108) 'epstein' , epstein
              write (iounit,108) 'epsthin' , epsthin
              write (iounit,108) 'epstzin' , epstzin
              write (iounit,108) 'tauhin'  , tauhin
              write (iounit,108) 'tauzin'  , tauzin
              write (iounit,108) 'fnzin'   , fnzin
              write (iounit,108) 'czin'    , czin
              write (iounit,108) 'azin'    , azin
              write (iounit,108) 'fnsin'   , fnsin
              write (iounit,108) 'betahin' , betahin
              write (iounit,108) 'betazin' , betazin
              write (iounit,108) 'ftrapein', ftrapein
              write (iounit,108) 'ekyrhoin', ekyrhoin
              write (iounit,108) 'ekparlin', ekparlin
  108         format (1x, a8, ' = ', 1pe14.6, ',')
      end if
c
      zepsnh = epsnhin
      zepsnz = epsnzin
      zepsns = epsnsin
      zepsth = epsthin
      zepste = epstein
      zepstz = epstzin
c
c ... check validity of input data
c
      if (neq .lt. 2) then
        call giveupus(iounit)
        close (iounit)
        call STOP ('subroutine ETAWN8_12: NEQ < 2', 207)
      end if
c
      if (ndim .gt. idim_l) then
        call giveupus(iounit)
        close (iounit)
        call STOP ('subroutine ETAWN8_12: NDIM > IDIM_L', 208)
      end if
c
      if (ABS (zepste) .lt. zepsqrt)
     .   zepste = SIGN (zepsqrt, zepste)
c
      if (ABS (zepsth) .lt. zepsqrt)
     .   zepsth = SIGN (zepsqrt, zepsth)
c
      if (ABS (zepsnh) .lt. zepsqrt)
     .   zepsnh = SIGN (zepsqrt, zepsnh)
c
      if (ABS (zepstz) .lt. zepsqrt)
     .   zepstz = SIGN (zepsqrt , zepstz)
c
      if (ABS (zepsnz) .lt. zepsqrt)
     .   zepsnz = SIGN (zepsqrt , zepsnz)
c
      if (ABS (zepsns) .lt. zepsqrt)
     .   zepsns = SIGN (zepsqrt , zepsns)
c
c ... initialize arrays
c
      do j1=1,ndim
        omega  (j1) = 0.0
        gamma  (j1) = 0.0
        chieff (j1) = 0.0
        velthi (j1) = 0.0
        perform(j1) = 0.0
        do j2=1,ndim
          difthi(j1,j2) = 0.0
          zchim (j1,j2) = 0.0
          zflxm (j1,j2) = 0.0
        end do
      end do
c
c ... set up initial gradients
c
      zgrdth = 1.0 / zepsth
      zgrdnh = 1.0 / zepsnh
      zgrdte = 1.0 / zepste
      zgrdnz = 1.0 / zepsnz
      zgrdtz = 1.0 / zepstz
      zgrdns = 1.0 / zepsns
c
      zimp   = czin
      zmass  = azin
      zfnz   = fnzin * zimp
      zfns   = fnsin
c
      zgrdne = zgrdnh * ( 1.0 - zfnz - zfns ) + zgrdnz * zfnz
     .  + zgrdns * zfns
      zgrdph = zgrdth + zgrdnh
      zgrdpe = zgrdte + zgrdne
      zgrdpz = zgrdtz + zgrdnz
c
c ... set up gradient matrix zgm(j1,jd)
c
c  zdg = the finite difference used to construct the zgm matrix
c
      zdg = cetain(30)
      if (zdg .lt. zepsqrt ) zdg = 0.01
c
c ... JAK 10/12/95 Error fixed: reset zgm to zero
c
      do k=1,5
      do jd=1,imatrx+1
        zgm(k,jd) = 0.0
      end do
      end do
c
c  input values for zgm matrix
c
      do jd=1,imatrx+1
        zgm(1,jd) = zgrdth
        zgm(2,jd) = zgrdnh
        zgm(3,jd) = zgrdte
        if (ieq .gt. 4) then
          zgm(4,jd) = zgrdnz
          zgm(5,jd) = zgrdtz
        end if
      end do
c
c  incremental values
c  Note:  if zgm(jd,jd) is negative, increment in negative direction
c
      do jd=1,imatrx+1
        zgm(jd,jd+1) = zgm(jd,jd+1) + SIGN (zdg, zgm(jd,jd+1))
      end do
c
      if (lprint .gt. 8) then
c
        if (iounit .gt. 0)  write (iounit, *)
        if (iounit .gt. 0)
     .      write (iounit, *) '  zgm(j1,jd) matrix,  -> jd'
        do j1=1,ieq-1
          if (iounit .gt. 0)  write (iounit, 132) (zgm(j1,jd),jd=1,ieq)
        end do
c
      end if
c
c ... loop over gradient matrix
c
      do jd=1,imatrx+1
c
c ... set variables
c
      do j1=1,ndim
        zomega(j1)  =  0.0
        zgamma(j1)  =  0.0
        zalpha(j1)  = (0.0, 0.0)
        zalfr(j1)   =  0.0
        zalfi(j1)   =  0.0
        zbeta(j1)   =  0.0
        do j2=1,ndim
          zevec  (j1,j2) = (0.0, 0.0)
          zeigenv(j1,j2) =  0.0
        end do
      end do
c
c ... compute gradient scale lengths from the zgm(j1,jd) matrix
c
      zgrdth = zgm(1,jd)
      zgrdnh = zgm(2,jd)
      zgrdte = zgm(3,jd)
      zgrdnz = zgm(4,jd)
      zgrdtz = zgm(5,jd)
c
      zgrdne = zgrdnh * ( 1.0 - zfnz - zfns ) + zgrdnz * zfnz
     .  + zgrdns * zfns
      zgrdph = zgrdph + zgrdnh
      zgrdpe = zgrdpe + zgrdne
      zgrdpz = zgrdpz + zgrdnz
c
      if (ABS (zgrdth) .lt. zepsqrt)
     .   zgrdth = SIGN (zepsqrt , zgrdth)
      if (ABS (zgrdnh) .lt. zepsqrt)
     .   zgrdnh = SIGN (zepsqrt , zgrdnh)
      if (ABS (zgrdte) .lt. zepsqrt)
     .   zgrdte = SIGN (zepsqrt , zgrdte)
      if (ABS (zgrdnz) .lt. zepsqrt)
     .   zgrdnz = SIGN (zepsqrt , zgrdnz)
      if (ABS (zgrdtz) .lt. zepsqrt)
     .   zgrdtz = SIGN (zepsqrt , zgrdtz)
      if (ABS (zgrdne) .lt. zepsqrt)
     .   zgrdne = SIGN (zepsqrt , zgrdne)
c
      zepsth = 1.0 / zgrdth
      zepsnh = 1.0 / zgrdnh
      zepste = 1.0 / zgrdte
      zepsne = 1.0 / zgrdne
      if (ieq .gt. 4) then
        zepsnz = 1.0 / zgrdnz
        zepstz = 1.0 / zgrdtz
      else
        zepsnz = 0.0
        zepstz = 0.0
      end if
c
c ... compute the rest of the dimensionless variables needed
c
      zetah  = zepsnh / zepsth
      zetae  = zepsne / zepste
c
      ztauh  = tauhin
c
      zep2nh = 2.0 * zepsnh
      zep2ne = 2.0 * zepsne
c
      zft    = ftrapein
      zflh   = ekyrhoin**2
c
      if (neq .gt. 4) then
        zetaz  = zepsnz / zepstz
        ztauz  = tauzin / czin
        zep2nz = 2.0 * zepsnz
        zflz   = zmass * zflh / zimp**2
      else
        zetaz  = 0.0
        ztauz  = 1.0
        zep2nz = 0.0
        zflz   = 0.0
      end if
c
c ... Note:
c
c     ztauz = T_Z / ( Z T_e )
c     zfnz  = Z n_Z / n_e
c     zimp  = Z
c     zmass = m_Z / m_H
c     zflh  = k_y^2 \rho_{sH}^2
c     zflz  = k_y^2 \rho_{sZ}^2
c
c ... diagnostic output
c
      if (lprint .gt. 2 .and. iounit .gt. 0) then
        write (iounit, *)
        write (iounit, *) '--------------------------------------'
        write (iounit, *)
        write (iounit, *) ' jd = ',jd
        if (lprint .gt. 2 .or. (lprint .gt. 0 .and. jd .eq. 1)) then
          write (iounit, *)
          write (iounit, *) zgrdph,' = zgrdph'
          write (iounit, *) zgrdpe,' = zgrdpe'
          write (iounit, *) zgrdpz,' = zgrdpz'
          write (iounit, *) zgrdnh,' = zgrdnh'
          write (iounit, *) zgrdne,' = zgrdne'
          write (iounit, *) zgrdnz,' = zgrdnz'
          write (iounit, *) zgrdth,' = zgrdth'
          write (iounit, *) zgrdte,' = zgrdte'
          write (iounit, *) zgrdtz,' = zgrdtz'
          write (iounit, *)
          write (iounit, *) zepsnh,' = zepsnh'
          write (iounit, *) zepsne,' = zepsne'
          write (iounit, *) zepsnz,' = zepsnz'
          write (iounit, *) zepsth,' = zepsth'
          write (iounit, *) zepste,' = zepste'
          write (iounit, *) zepstz,' = zepstz'
          write (iounit, *)
          write (iounit, *) zetah,' = zetah'
          write (iounit, *) zetaz,' = zetaz'
          write (iounit, *) zetae,' = zetae'
          write (iounit, *) ztauh,' = ztauh'
          write (iounit, *) ztauz,' = ztauz'
          write (iounit, *) zep2nh,' = zep2nh'
          write (iounit, *) zep2nz,' = zep2nz'
          write (iounit, *) zep2ne,' = zep2ne'
          write (iounit, *)
          write (iounit, *) zft,' = zft'
          write (iounit, *) zimp,' = zimp'
          write (iounit, *) zmass,' = zmass'
          write (iounit, *) zfnz,' = zfnz'
          write (iounit, *) zflh,' = zflh'
          write (iounit, *) zflz,' = zflz'
          write (iounit, *)
          write (iounit, *) zepsqrt,' = zepsqrt'
          write (iounit, *) zepsmach,' = zepsmach'
        end if
      end if
c
c ... set matrices for eigenvalue equation
c
      do j1=1,idim_l
        zalfr(j1) = 0.0
        zalfi(j1) = 0.0
        zbeta(j1) = 0.0
        do j2=1,idim_l
          zam(j1,j2) = 0.0
          zbm(j1,j2) = 0.0
          zeigenv(j1,j2) = 0.0
        end do
      end do
c
      if (ieq .eq. 6) then
c
c ... Six equations with impurities, trapped electrons, and FLR
c
c  equations for e \phi / T_e, T_H, n_H, T_{et}, n_Z, and T_Z
c
      if (lprint .gt. 1 .and. jd .eq. 1) then
        if (iounit .gt. 0)  write (iounit, *)
        if (iounit .gt. 0)  write (iounit, *)
     .   'Six eqns for e phi/T_e, T_H, n_H, T_{et}, n_Z, and T_Z'
      end if
c
c  hydrogen density
c
        zam(1,1) = - 1.0
     .   + ( 1.0 - zflh * ztauh * ( 1.0 + zetah ) ) / zep2nh
        zam(1,2) = - ztauh
        zam(1,3) = - ztauh
c
        zbm(1,1) = zflh
        zbm(1,3) = 1.0
c
c  hydrogen energy
c
        zam(2,1) = ( zetah - 2.0 / 3.0 ) / zep2nh
        zam(2,2) = - ztauh * 5.0 / 3.0
c
        zbm(2,2) = 1.0
        zbm(2,3) = - 2.0 / 3.0
c
c  trapped electron density
c
        zam(3,1) = - 1.0 + zft / zep2ne
        zam(3,3) = 1.0 - zfnz - zfns
        zam(3,4) = zft
        zam(3,5) = zfnz
c
        zbm(3,1) = zft - 1.0
        zbm(3,3) = 1.0 - zfnz - zfns
        zbm(3,5) = zfnz
c
c  trapped electron energy
c
        zam(4,1) = zft * ( zetae - 2.0 / 3.0 ) / zep2ne
        zam(4,4) = zft * 5.0 / 3.0
c
        zbm(4,1) = ( 1.0 - zft ) * 2.0 / 3.0
        zbm(4,3) = - ( 1.0 - zfnz - zfns ) * 2.0 / 3.0
        zbm(4,4) = zft
        zbm(4,5) = - zfnz * 2.0 / 3.0
c
c  impurity density
c
        zam(5,1) = - 1.0
     .    + ( 1.0 - zflz * ztauz * ( 1.0 + zetaz ) ) / zep2nz
        zam(5,5) = - ztauz
        zam(5,6) = - ztauz
c
        zbm(5,1) = zflz
        zbm(5,5) = 1.0
c
c  impurity energy
c
        zam(6,1) = ( zetaz - 2.0 / 3.0 ) / zep2nz
        zam(6,6) = - ztauz * 5.0 / 3.0
c
        zbm(6,5) = - 2.0 / 3.0
        zbm(6,6) = 1.0
c
c ... 4 equations with trapped electrons and FLR effects
c
      else if (ieq .eq. 4) then
c
c  equations for e phi/T_e, T_H, n_i, and T_e
c
       if (lprint .gt. 1) then
         if (iounit .gt. 0)  write (iounit, *)
         if (iounit .gt. 0)  write (iounit, *)
     .              ' Four eqns for e phi/T_e, T_H, n_H, and T_e'
       end if
c
c  ion continuity
c
        zam(1,1) = 1.0 - zep2nh - zflh * ztauh * (1.0 + zetah)
        zam(1,2) =     - zep2nh        * ztauh
        zam(1,3) =     - zep2nh        * ztauh
c
        zbm(1,1) = zflh * zep2nh
        zbm(1,3) = zep2nh
c
c  ion energy
c
        zam(2,1) = zetah - 2.0 / 3.0
        zam(2,2) = - zep2nh * ztauh * 5.0 / 3.0
c
        zbm(2,2) =   zep2nh
        zbm(2,3) = - zep2nh * 2.0 / 3.0
c
c  trapped electron continuity
c
        zam(3,1) = zft - zep2ne
        zam(3,3) = zep2ne
        zam(3,4) = zft * zep2ne
c
        zbm(3,1) = ( zft - 1.0 ) * zep2ne
        zbm(3,3) = zep2ne
c
c  trapped electron energy
c
        zam(4,1) = zft * ( zetae - 2.0 / 3.0)
        zam(4,4) = zft * zep2ne * 5.0 / 3.0
c
        zbm(4,1) = ( 1.0 - zft ) * zep2ne * 2.0 / 3.0
        zbm(4,3) = - zep2ne * 2.0 / 3.0
        zbm(4,4) = zft * zep2ne
c
c ... two equations when trapped particles and FLR effects omitted
c
      else if (ieq .eq. 2) then
c
c  equations for e phi/T_e and T_H
c
       if (lprint .gt. 1 .and. jd .eq. 1) then
         if (iounit .gt. 0)  write (iounit, *)
         if (iounit .gt. 0)  write (iounit, *)
     .                          ' Two eqns for e phi/T_e and T_H'
       end if
c
        zam(1,1) = ( 1.0 / zep2nh ) - ztauh - 1.0
        zam(1,2) = - ztauh
        zam(2,1) = ( zetah - 2.0 / 3.0 ) / zep2nh
        zam(2,2) = - ztauh * 5.0 / 3.0
c
        zbm(1,1) = 1.0
        zbm(1,2) = 0.0
        zbm(2,1) = - 2.0 / 3.0
        zbm(2,2) = 1.0
c
      else
c
        write (6, *)
        write (6, *)  ieq, ' = ieq in subroutine ETAWN8_12'
        call giveupus(iounit)
        close (iounit)
        call STOP ('subroutine ETAWN8_12: wrong value #2 of IEQ', 210)
c
      end if
c
c ... save copy of matrix which is over-written by subroutine F02BJE
c
      do j2=1,ieq
        do j1=1,ieq
          zamt(j1,j2) = zam(j1,j2)
          zbmt(j1,j2) = zbm(j1,j2)
        end do
      end do
c
c ... diagnostic output
c
      if (lprint .gt. 6) then
        if (iounit .gt. 0)  write (iounit, *)
        if (iounit .gt. 0)  write (iounit, *) ' zam(j1,j2)  j2 ->'
        do j1=1,ieq
          if (iounit .gt. 0)  write (iounit, 192) (zam(j1,j2), j2=1,ieq)
        end do
c
        if (iounit .gt. 0)  write (iounit, *)
        if (iounit .gt. 0)  write (iounit, *) ' zbm(j1,j2)  j2->'
        do j1=1,ieq
          if (iounit .gt. 0)
     .        write  (iounit, 192)  (zbm(j1,j2), j2=1,ieq)
  192         format (1p6e13.5)
        end do
      end if
c
      if (letain(6) .eq. 2) then
c
c ... find the eigenvalues and eigenvectors using IMSL routine gvcrg
c
cjak    call gvcrg (ieq, zam, idim_l, zbm, idim_l, zalpha, zbeta
cjak &    , zevec, idim_l)
c
cjak    perform(jd) = gpirg (ieq, ieq, zam, idim_l, zbm, idim_l
cjak &    , zalpha, zbeta, zevec, idim_l)
c
        if (lprint .gt. 1) then
          if (iounit .gt. 0)  write (iounit, *)
          if (iounit .gt. 0)  write (iounit, *)
     .            ' Eigenvalues computed using IMSL routine GVCRG'
          if (iounit .gt. 0)  write (iounit, *)
          if (iounit .gt. 0)  write (iounit, *)
     .            perform(jd),' = performance index from IMSL GPIRG'
        end if
c
        ifail = 0
c
      else if (letain(6) .eq. 1) then
c
c ... find the eigenvalues and eigenvectors using NAG14 routine F02BJE
c
        ztol  = MAX (zero, cetain(32))
        ifail = -1
        lmatv = .true.
c
cjak    call F02BJE ( ieq, zam, idim_l, zbm, idim_l, ztol
cjak  .  , zalfr, zalfi, zbeta, lmatv, zeigenv, idim_l, iter, ifail )
c
****    if (ifail .gt. 0) call abortb ( 6
****  .  , 'ifail .gt. 0 after call F02BJE in sbrtn ETAWN8_12')
c
        if (lprint .gt. 1) then
          if (iounit .gt. 0)  write (iounit, *)
          if (iounit .gt. 0)  write (iounit, *)
     .               ' Eigenvalues computed using NAG routine F02BJE'
          if (iounit .gt. 0)  write (iounit, *)
        end if
c
      else if (letain(6) .eq. 0) then
c
c ... find the eigenvalues and eigenvectors using NAG14 routine F02AKE
c
      if (ieq.gt.idp) then
        ieq = idp
        write (6, *)  ' WARNING from ETAWN8_12: IEQ set to = ', ieq
      end if
      do j2=1,ieq
        do j1=1,ieq
          zamt(j1,j2) = zam(j1,j2)
          zbmt(j1,j2) = zbm(j1,j2)
          bori(j1,j2) = zbm(j1,j2)
****      zbm (j1,j2) = 0.0
****      if (j1 .eq. j2)  zbm(j1,j2) = j1
          binv(j1,j2) = 0.0
        end do
      end do
c
c ... first invert zbm
c
      ifail = 0
c
****  write (6,*)
****  write (6,*) ' zbm(j1,j2)  j2->'
****  do j1=1,ieq
****    write (6, 192)  (zbm(j1,j2),j2=1,ieq)
****  end do
c
cjak  call F01AAE (bori,idp,ieq,binv,idp,wkspce,ifail)
c
      if (ifail .ne. 0) then
        write (6, *)  'WARNING from ETAWN8_12: F01AAE ifail = ', ifail
        return
      end if
c
c ... then calculate eigenvalues of Binv.A
c
      do j2=1,ieq
        do j1=1,ieq
          ar(j1,j2) = 0.0
          ai(j1,j2) = 0.0
          do k=1,ieq
            ar(j1,j2) = ar(j1,j2) + binv(j1,k)*zam(k,j2)
          end do
        end do
      end do
c
      ifail = 0
c
cjak  call F02AKE ( ar,iar,ai,iai,ieq,rr,ri,vr,ivr,vi,ivi,
cjak .              intger, ifail )
c
      if (ifail.ne.0) then
        write (6, *)  'WARNING from ETAWN8_12: F01AKE ifail = ', ifail
        return
      end if
c
        if (lprint .gt. 1) then
          if (iounit .gt. 0)  write (iounit, *)
          if (iounit .gt. 0)  write (iounit, *)
     .           ' Eigenvalues computed using NAAG routine F02AKE'
          if (iounit .gt. 0)  write (iounit, *)
        end if
c
      else if (letain(6) .eq. 3) then
c
c ... find the eigenvalues and eigenvectors using EISPACK routine RGG
c
        matz = 1
        ierr = 0
c
c ...   convert to double precision ; double precision is not used JAK
c
        do i=1,idp
          do j=1,idp
cjak        dzam(i,j) = DBLE (zam(j,i))
cjak        dzbm(i,j) = DBLE (zbm(j,i))
            dzam(i,j) =       zam(i,j)
            dzbm(i,j) =       zbm(i,j)
          end do
        end do
c
        ierr = 0
        call rgg (idim_l, ieq, dzam, dzbm, dzalfr, dzalfi, dzbeta,
     .            matz, dzeigenv, ierr)
c
c ...   convert to single precision
c
        do i=1,idp
          zalfr(i) = dzalfr(i)
          zalfi(i) = dzalfi(i)
          zbeta(i) = dzbeta(i)
          do j=1,idp
            zeigenv(i,j) = dzeigenv(i,j)
          end do
        end do
c
      end if
c
c ... compute the complex eigenvalues
c
      if (letain(6) .eq. 2) then
c
        do j=1,ieq
          ztemp1 = zbeta(j)
          if (ABS (zbeta(j)) .lt. zepsqrt) ztemp1 = zepsqrt
          zomega(j) =  REAL (zalpha(j)) / ztemp1
          zgamma(j) = PIMAG (zalpha(j)) / ztemp1
          zalfr (j) =  REAL (zalpha(j))
          zalfi (j) = PIMAG (zalpha(j))
        end do
c
      else if (letain(6) .eq. 1) then
c
        do j=1,ieq
          ztemp1 = zbeta(j)
          if (ABS (zbeta(j)) .lt. zepsqrt)  ztemp1 = zepsqrt
          zalpha(j) = CMPLX (zalfr(j), zalfi(j))
          zomega(j) = zalfr(j) / ztemp1
          zgamma(j) = zalfi(j) / ztemp1
        end do
c
        do j2=1,ieq
****      if (ABS (zomega(j2)) .lt. zepsqrt) then
          if (ABS (zgamma(j2)) .lt. zepsqrt) then
            do j1=1,ieq
              zevec(j1,j2) = CMPLX (zeigenv(j2,j1), zero)
            end do
          else if (ABS ( zomega(j2+1) - zomega(j2) ) .lt. zepsqrt) then
            do j1=1,ieq
              zevec(j1,j2) = CMPLX (zeigenv(j2,j1), zeigenv(j2,j1+1))
            end do
          else if (ABS ( zomega(j2-1) - zomega(j2) ) .lt. zepsqrt) then
            do j1=1,ieq
              zevec(j1,j2) = CMPLX (zeigenv(j2,j1-1), -zeigenv(j2,j1))
            end do
          end if
        end do
c
      else if (letain(6) .eq. 0) then
c
        do j=1,ieq
          zomega(j) = rr(j)
          zgamma(j) = ri(j)
          do j1=1,ieq
            zevec(j1,j) = CMPLX (vr(j1,j), vi(j1,j))
          end do
        end do
c
      else if (letain(6) .eq. 3 ) then
c
c ...   JAK 95/11/17 First calculate zomega and zgamma
c
        do j=1,ieq
          ztemp1 = zbeta(j)
          if (ABS (zbeta(j)) .lt. zepsqrt ) ztemp1 = zepsqrt
          zalpha(j) = CMPLX (zalfr(j), zalfi(j))
          zomega(j) = zalfr(j) / ztemp1
          zgamma(j) = zalfi(j) / ztemp1
        end do
        do j=1,ieq
          ztemp1 = zbeta(j)
          if (ABS (zbeta(j)) .lt. zepsqrt ) ztemp1 = zepsqrt
          do i=1,ieq
            if (ABS (zgamma(j)) .lt. zepsqrt ) then
              zevec(i,j) = CMPLX (zeigenv(i,j), zero)
            else if (ABS ( zomega(j+1) - zomega(j) ) .lt. zepsqrt) then
              zevec(i,j) = CMPLX (zeigenv(i,j), zeigenv(i,j+1))
            else if (ABS ( zomega(j-1) - zomega(j) ) .lt. zepsqrt) then
              zevec(i,j) = CMPLX (zeigenv(i,j-1), -zeigenv(i,j))
            end if
c
c ...       JAK 95/11/16 note that this normalization is not essential
c ...       here, as it is taken care of later
c
c           zevec(i,j) = zevec(i,j)/ztemp1
          end do
        end do
c
      end if
c
c ... save growth rates and frequencies during first element of jd loop
c
      if (jd .eq. 1) then
        zgamax = 0.0
        do j=1,ieq
          omega(j) = zomega(j)
          gamma(j) = zgamma(j)
          zgamax = MAX (zgamax, zgamma(j))
        end do
        if (zgamax .lt. zepsqrt)  go to 90
      end if
c
c ... check the eigenfunctions
c
      if (lprint .gt. 12) then
c
        if (iounit .gt. 0)  write (iounit, *)
        if (iounit .gt. 0)  write (iounit, *) ' Checking eigenfunctions'
c
c  Real and imaginary parts
c
        do j=1,ieq-1
c
          do j1=1,ieq
c
            ztempa(j1) = 0.0
c
c ... error fixed: ztmpa should be ztempb JAK 9/28/95
c
            ztempb(j1) = 0.0
c
            do j2=1,ieq
c
              ztempa(j1) = ztempa(j1)
     .          + zamt(j1,j2) *  REAL (zevec(j2,j))
     .          - zbmt(j1,j2) *  REAL (zevec(j2,j)) * zomega(j)
     .          + zbmt(j1,j2) * PIMAG (zevec(j2,j)) * zgamma(j)
              ztempb(j1) = ztempb(j1)
     .          + zamt(j1,j2) * PIMAG (zevec(j2,j))
     .          - zbmt(j1,j2) * PIMAG (zevec(j2,j)) * zomega(j)
     .          - zbmt(j1,j2) *  REAL (zevec(j2,j)) * zgamma(j)
c
             end do
c
           end do
c
           if (iounit .gt. 0)  write (iounit, *)
           if (iounit .gt. 0)  write (iounit, *)
     .                              ' LHS - RHS for j =  ',j
           do j1=1,ieq
             if (iounit .gt. 0)
     .         write (iounit,142)  ztempa(j1), ztempb(j1)
           end do
c
        end do
c
      end if
c
 142  format (1p6e12.4)
c
c ... compute effective diffusivities directly from eigenvalues
c     assume eigenvectors are arranged in the order of
c     e\phi/T_e, T_H, n_H, T_{et}, n_Z, T_Z
c
        nmodes = 0
c
        do j=1,ieq-1
c
          ztemp1 =  REAL (zevec(1,j)) *  REAL (zevec(1,j))
     .           + PIMAG (zevec(1,j)) * PIMAG (zevec(1,j))
c
          if (zgamma(j) .gt. zepsqrt .and. ztemp1 .gt. zepsqrt) then
c
            nmodes = nmodes + 1
c
            zreal  =  REAL (zevec(2,j)) +  REAL (zevec(3,j))
            zimag  = PIMAG (zevec(2,j)) + PIMAG (zevec(3,j))
c
            zflxph = - 2.0 * (ABS (zgamma(j)))**2
     .        * ( zimag *  REAL (zevec(1,j))
     .          - zreal * PIMAG (zevec(1,j)) ) / ztemp1
c
            zflxm(1,jd) = zflxm(1,jd) + zflxph
c
c            if (lprint .gt. 7) then
c              if (iounit .gt. 0)  write (iounit, *)
c              if (iounit .gt. 0)  write (iounit, *)
c    .                             zflxph,' = zflxph for j = ',j
c            end if
c
            zflxnh = - 2.0 * ( abs ( zgamma(j) ) )**2
     .        * ( PIMAG (zevec(3,j)) *  REAL (zevec(1,j))
     .          -  REAL (zevec(3,j)) * PIMAG (zevec(1,j)) ) / ztemp1
c
            zflxm(2,jd) = zflxm(2,jd) + zflxnh
c
c            if (lprint .gt. 7) then
c              if (iounit .gt. 0)  write (iounit, *)
c              if (iounit .gt. 0)  write (iounit, *)
c     .                            zflxnh,' = zflxnh for j = ',j
c            end if
c
            if (ieq .gt. 3) then
c
              zreal =  REAL (zevec(4,j))
     .          + ( 1.0 - zfnz - zfns ) *  REAL (zevec(2,j))
     .          + zfnz * REAL (zevec(5,j))
              zimag = PIMAG (zevec(4,j))
     .          + ( 1.0 - zfnz - zfns ) * PIMAG (zevec(2,j))
     .          + zfnz * PIMAG (zevec(5,j))
c
c  Note, the electron heat flux is reduced by the fraction of
c    trapped electrons
c
              zflxpe =
     .          - 2.0 * zft * (ABS ( zgamma(j) ))**2
     .        * ( zimag *  REAL (zevec(1,j))
     .          - zreal * PIMAG (zevec(1,j)) ) / ztemp1
c
              zflxm(3,jd) = zflxm(3,jd) + zflxpe
c
c              if (lprint .gt. 7) then
c                if (iounit .gt. 0)  write (iounit, *)
c                if (iounit .gt. 0)  write (iounit, *)
c     .                              zflxpe,' = zflxpe for j = ',j
c              end if
            end if
c
            if (ieq .gt. 4) then
              zflxnz = - 2.0 * (ABS ( zgamma(j) ))**2
     .        * ( PIMAG (zevec(5,j)) *  REAL (zevec(1,j))
     .          -  REAL (zevec(5,j)) * PIMAG (zevec(1,j)) ) / ztemp1
c
              zflxm(4,jd) = zflxm(4,jd) + zflxnz
c
c              if (lprint .gt. 7) then
c                if (iounit .gt. 0)  write (iounit, *)
c                if (iounit .gt. 0)  write (iounit, *)
c     .                              zflxnz,' = zflxnz for j = ',j
c              end if
            end if
c
            if (ieq .gt. 5) then
c
              zreal =  REAL (zevec(5,j)) +  REAL (zevec(6,j))
              zimag = PIMAG (zevec(5,j)) + PIMAG (zevec(6,j))
c
              zflxpz = - 2.0 * (ABS ( zgamma(j) ))**2
     .        * ( zimag *  REAL (zevec(1,j))
     .          - zreal * PIMAG (zevec(1,j)) ) / ztemp1
c
              zflxm(5,jd) = zflxm(5,jd) + zflxpz
c
****          if (lprint .gt. 7) then
****            if (iounit .gt. 0)  write (iounit, *)
****            if (iounit .gt. 0)  write (iounit, *)
****  .                             zflxpz,' = zflxpz for j = ', j
****          end if
c
            end if
c
          end if
c
        end do
c
c ... compute effective total diffusivities
c
        zchim(1,jd) = zflxm(1,jd) * zepsth
        zchim(2,jd) = zflxm(2,jd) * zepsnh
        zchim(3,jd) = zflxm(3,jd) * zepste
        zchim(4,jd) = zflxm(4,jd) * zepsnz
        zchim(5,jd) = zflxm(5,jd) * zepstz
c
      if (lprint .gt. 2) then
c
c ... print eigenvalues and eigenfunctions
c
        if (iounit .gt. 0)  write (iounit,121)
        do j=1,ieq
          if (iounit .gt. 0)  write (iounit, 122) zomega(j), zgamma(j),
     .                        zalfr(j), zalfi(j), zbeta(j), iter(j)
        end do
 121    format (/ ' Solution of the eigenvalue equations' /
     .     t4,'zomega',t18,'zgamma',t32,'zalfr',t46,'zalfi',t60,'zbeta',
     .     t74,'iter')
 122    format (1p5e14.5,i5)
c
        if (iounit .gt. 0)  write (iounit, *)
        if (iounit .gt. 0)  write (iounit, *)
     .    ' Effective diffusivities normalized by omega_{De} / k_y^2'
c
        if (iounit .gt. 0)  write (iounit, 130)
        if (iounit .gt. 0)  write (iounit, 132)
     .                            (zchim(j1,jd), j1=1,ieq-1)
c
      end if
c
      if (lprint .gt. 99) then
c
        if (iounit .gt. 0)  write (iounit, *)
        if (iounit .gt. 0)  write (iounit, *)
     .                    ' Eigenvectors zevec(j1,j2) j2->'
        do j1=1,ieq
          if (iounit .gt. 0)  write (iounit,124) (zevec(j1,j2),j2=1,ieq)
  124     format (1p6e11.3)
        end do
c
        if (iounit .gt. 0)  write (iounit, *)
        if (iounit .gt. 0)  write (iounit, *)
     .                           ' Eigenvectors zeigenv(j1,j2) j2->'
        do j1=1,ieq
          if (iounit .gt. 0)  write (iounit,124)
     .                              (zeigenv(j1,j2),j2=1,ieq)
        end do
c
        if (iounit .gt. 0)  write (iounit, *)
        if (iounit .gt. 0)  write (iounit, *)
     .            ztol, ' = tolerance used in subroutine F02BJE'
c
      end if
c
      if (ifail .gt. 0) then
        call giveupus(iounit)
        close (iounit)
        call STOP ('subroutine ETAWN8_12: IFAIL > 0 after F02BJE', 211)
      end if
c
c ... end of loop over diffusivity matrix
c
      end do
c
c ... computation of difthi(j1,j2) from zflxm(j1,jd) and zgm(j1,jd)
c
      do j1=1,imatrx
c
        ztemp2 = 0.0
c
        do jd=1,imatrx
c
          ztemp1 = zgm(jd,jd) - zgm(jd,1)
          if (ABS (ztemp1) .lt. zepsqrt)
     .       ztemp1 = zdg * SIGN (one, zgm(jd,1))
c
          difthi(j1,jd) = ( zflxm(j1,jd+1) - zflxm(j1,1) ) / ztemp1
          ztemp2 = ztemp2 + difthi(j1,jd) * zgm(jd,1)
c
        end do
c
        if (letain(7) .eq. 0) then
c
c ... compute convective velocities
c
          velthi(j1) = zflxm(j1,1) - ztemp2
c
        else
c
c ... alternatively, set the convective velocities to zero
c    and rescale the diffusivity matrix
c
          velthi(j1) = 0.0
          if (ABS (ztemp2) .lt. zepsqrt)
     .      ztemp2 = SIGN (zepsqrt, ztemp2)
          do jd=1,imatrx
            difthi(j1,jd) = difthi(j1,jd) * zflxm(j1,1) / ztemp2
          end do
c
        end if
c
c ... save effective diffusivities
c
        chieff(j1) = zflxm(j1,1) / zgm(j1,1)
c
      end do
c
c ... diagnostic printout
c
      if (lprint .gt. 6) then
c
        if (iounit .gt. 0)  write (iounit, *)
        if (iounit .gt. 0)  write (iounit, *)
     .                    ' matrix zflxm(j1,jd)  -> jd'
        do j1=1,ieq-1
          if (iounit .gt. 0)  write (iounit,132)
     .                       (zflxm(j1,jd),jd=1,ieq-1)
        end do
c
      end if
c
  90  if (lprint .gt. 0) then
c
        if (iounit .gt. 0)  write (iounit, *)
        if (iounit .gt. 0)  write (iounit, *)
     .          '---------------------------------------',
     .          'Final results from subroutine ETAWN8_12'
c
        if (iounit .gt. 0)  write (iounit, *)
        if (iounit .gt. 0)  write (iounit, *)
     .     ' diffusion matrix from eigenvectors'
     .    ,' normalized by omega_{De} / k_y^2'
        if (iounit .gt. 0)  write (iounit, *)
        if (iounit .gt. 0)  write (iounit, *) ' difthi(j1,j2)  -> j2'
        do j1=1,imatrx
          if (iounit .gt. 0)  write (iounit,132)
     .                              (difthi(j1,j2),j2=1,imatrx)
        end do
c
        if (iounit .gt. 0)  write (iounit, *)
        if (iounit .gt. 0)  write (iounit, *) ' convective velocities'
     .    ,' normalized by omega_{De} / k_y'
        do j1=1,imatrx
          if (iounit .gt. 0)  write (iounit,132) velthi(j1)
        end do
c
        if (iounit .gt. 0)  write (iounit, *)
        if (iounit .gt. 0)  write (iounit, *)
     .     ' Effective diffusivities from eigenvectors'
     .    ,' normalized by omega_{De} / k_y^2'
c
        if (iounit .gt. 0)  write (iounit, 130)
        if (iounit .gt. 0)  write (iounit, 132) (chieff(j1),j1=1,imatrx)
c
      end if
c
c ... close output file
c
      if (lprint .ge. 0 .and. iounit .gt. 0) then
        call giveupus(iounit)
        close (iounit)
      end if
c
      return
c
  130 format (/ t3,'chi_i',t16,'dif_h',t29,'chi_e',
     .                     t42,'dif_Z',t55,'chi_Z')
  132 format (1p6e13.5)
c
      end

      subroutine ip_chi2 (itest_chi, RLT, RLN, RLNe, q, kappa, shat,
     .                    zth, nbeam, tau, eps, gnu, g_perp, rmajor,
     .                    rho_i, v_ti, RLTcrit, RLTcritz, chi_0, g,
     .                    counter,
     .                    gamma, chi_i, chi_e)
c
c     modified to use REAL*8  HSJ 2/11/96
c     modified to use tshat everywhere instead of shat 2/20/97 HSJ
c
c     Note: Watch out for apostrophes (single quotes) in comments.
c
c The formulas embodied in this subroutine are documented in the Physics
c of Plasmas article entitled Quantitative Predictions of Tokamak Energy
c Confinement from First-Principles Simulations with Kinetic Effects,
c by M. Kotschenreuther, W. Dorland, M.A. Beer, and G.W. Hammett,
c Vol. 2, p. 2381, (1995). Extensions to non-circular cross-sections are
c described below.
c
c There is a significant typographical error in that paper.  R/Ln* is
c defined to be max(6,R/Ln); it should be min(6,R/Ln) as defined in
c this subroutine.
c
c Also, note that in deriving these formulas, we assumed that the
c density gradient scale lengths for the different species were equal.
c This is an approximation that needs to be relaxed.
c
c As emphasized in the paper, these formulas were derived numerically and
c are therefore not trustworthy outside a particular region of parameter
c space.  For example, we did not parameterize the heat flux in the weak
c magnetic shear limit; thus, one should not use the model in this limit.
c I have attempted to reduce related strange numerical behaviors by
c limiting some inputs to be roughly within their range of validity.
c
c Questions, problems, errors, etc. should be reported to
c bdorland@zonker.ph.utexas.edu or bdorland@pppl.gov.
c
c I will reply as quickly as possible.
c
c Stiffness:
c **********
c For many cases that we have simulated, the transport equations
c tend to be very stiff.  That is, the plasma temperature gradient scale
c length tends to adjust itself to be close to the critical gradient scale
c length over some region of the plasma, because chi becomes very large
c very fast for shorter temperature gradient scale lengths.  Typically,
c we have had to be very careful with the numerical algorithm used in the
c transport equation solver with this experience in mind.  The details
c of our implementation are available to anyone that is interested.
c
c Geometry:
c ********
c
c The nonlinear simulations that were done to obtain these formulas were
c mostly done in a simplified geometry, using a shifted circle, low beta,
c high aspect ratio expansion.  Some modifications due to more sophisticated
c geometrical models have been calculated and have been included here,
c but should be considered preliminary.  There are two important issues
c that must be noted.  First, we derived our formulas using a different
c radial coordinate.  Second, since we are actually calculating the
c transport coefficients in general geometry, we require less assumptions
c for the form of the transport equation to be solved.
c
c Let me describe the <|grad rho|> issue first:
c
c     The database standard modeling assumptions that were agreed upon for
c     this exercise include the assumption that the anomalous fluxes for
c     non-circular cross-sections are simply related to the anomalous
c     fluxes for related circular cross-section plasmas.  That is, in order
c     to get the factor of < |grad rho|**2 > that appears as a coefficient
c     of chi in the energy transport equations, one assumes that
c
c     chi_anom_general = chi_anom_circular * (grad rho).
c
c     One need not make this assumption; one can just calculate the quantity
c     chi_anom_general directly.  One would then have a transport equation
c     of the form
c
c   (3/2) d(n T)/dt = (1/V') d/drho V' <|grad rho|> n chi d/drho(T)] + ...
c
c     in which (grad rho) appears to the first power, rather than the second,
c     and chi is the thermal diffusivity from a general geometry theory.
c
c     This is arguably the better way to proceed, since
c
c         Vprime <|grad rho|> = A
c
c     where A is the surface area.  In this form, the quantity
c
c          -n chi dT/drho
c
c     can be identified as the heat flux per unit area, a natural
c     quantity from a theoretical/simulation point of view.  This chi is
c     the quantity returned by this subroutine.
c
c     If you are solving the transport equations in the form that
c     the ITER Expert Group agreed upon, e.g.,
c
c (3/2) d(n T)/dt = (1/V') d/drho V' <|grad rho|**2> n chi d/drho(T)] + ...
c
c     then you need to multiply the chi_i and chi_e reported by this
c     subroutine by the factor <|grad rho|>/<|grad rho|**2>.  This should
c     result in only small corrections to the predicted profiles.
c
c The choice of radial coordinate is more difficult to resolve:
c
c     We did not use the sqrt(toroidal flux) radial coordinate in our
c     non-circular cross-section simulations.  Instead, we used "rho_d",
c     where rho_d is defined to be the average horizontal minor radius
c     at the elevation of the magnetic axis, normalized to the value of
c     this quantity at the LCFS.
c
c     In other words, denote by "d" the horizontal diameter of a given flux
c     surface measured at the elevation of the magnetic axis.  Denote by
c     "D" the horizontal diameter of the last closed flux surface at the
c     elevation of the magnetic axis.  Then rho_d = d/D.  I believe this
c     is variable number 67 (RMINOR) in the ITER Profile Database Standard
c     List.
c
c     It is not difficult to allow for an arbitrary radial coordinate
c     in a transport code.  One must obtain all of the radial
c     quantities as functions of rho_d rather than rho via interpolation.
c
c     However, I do not expect everyone to go to this length to test our
c     model, since you agreed to use the sqrt(toroidal flux) definition
c     of the radial coordinate.  Thus, I suggest the following alternative:
c     Simply use the rho_d coordinate to define the scale lengths that
c     appear in the formulas below.  For most quantities (such as R/LT),
c     this simply amounts to including an additional factor d rho/d rho_d
c     in the expressions passed to this subroutine.  While not completely
c     correct, this workaround captures the dominant effect, related to
c     evaluating the flux near the critical gradient.
c
c ****** Summary of comments:  ***********
c
c (1) The general geometry extensions to the IFS/PPPL model were derived
c     using rho = d/D = RMINOR as the radial coordinate.  To be
c     most accurate, the transport equation should be solved using d/D
c     as the radial coordinate.
c
c     If you use rho proportional to sqrt(toroidal flux) instead of rho=d/D
c     as your radial coordinate, you should at least carefully define the
c     scale lengths as indicated below (using rho=d/D).
c
c (2) This routine should be used to return the thermal transport
c     coefficients (chi_i, chi_e) for energy transport equations of the form
c
c (3/2) d(n T)/dt = (1/V') d/drho V' <|grad rho|> n chi d/drho(T)] + ...
c
c     Note that <|grad rho|> only appears to the first power according
c     to this definition of chi.  If your code is hardwired to solve an
c     equation of the form
c
c (3/2) d(n T)/dt = (1/V') d/drho V' <|grad rho|**2> n chi d/drho(T)] + ...
c
c     then multiply the chi_i and chi_e obtained from this routine by
c     the factor <|grad rho|>/<|grad rho|**2>.
c
c *****************************************************************
c
c     RLT  R/L_Ti, where R = the major radius and
c                  1/L_Ti = -1/T_i dT_i/drho_d
c     RLN  R/L_ni, where 1/L_ni = -1/n_i dn_i/drho_d
c     RLNe R/L_ne, where 1/L_ne = -1/n_e dn_e/drho_d
c     q    The safety factor.
c     kappa  The elongation, defined here to be the
c          kappa = maximum height/maximum diameter
c     shat == rho_d/q dq/drho_d
c     zth  Thermal Z_eff.  The simulations that were carried out to
c          generate the formulae in this subroutine assumed the plasma
c          was composed of a thermal hydrogenic species, thermal carbon,
c          a hydrogenic beam species, and electrons.  We found that low-Z
c          impurities primarily act to dilute the main ion concentration,
c          and can accounted for to first order by modifying the
c          definition of Z_eff.  Some of the more important effects of
c          the fast ions in the plasma are also partially accounted
c          for by this parameter, which is:
c          zth == (n_i + 36 n_C)/(n_e - n_beam)
c
c     nbeam == local fast ion (beam) density normalized to the electron
c          density.
c     tau  == T_i/T_e.  Note that this is opposite to a widely
c                       used convention.
c     eps  == rho_d/R, the local minor radius normalized to the major radius.
c     gnu   Dimensionless collisionality parameter.
c          gnu == 2.5e-7 * n_e / (T_e**1.5 T_i**0.5) * rmajor
c          where n_e is in units of cm**-3, T_e and T_i are in eV, and
c          rmajor is in units of m.  For an R = 2.4 m, 100 eV, 1e13 plasma,
c          gnu=600.
c     beta  Local total beta; presently ignored.
c     g_perp == velocity shear parameter.  Use g_perp=0 (actual value
c          discussed in the Waltz paper cited below).
c     rmajor == major radius of the plasma
c     rho_i == local thermal gyroradius of thermal hydrogenic species.
c     v_ti == sqrt(T_i/m_i) where T_i and m_i are the local thermal
c          hydrogenic temperature and average thermal hydrogenic mass.
c
c     rho_e,v_te,etg,rlte: dummy parameters at present. Ignore.
c
c     Units: The only dimensional parameters in the inputs are the major
c     radius, rho_i, and v_t.  Their units should be consistent; the chis
c     that are returned will be in units of rho_i**2 v_ti / rmajor.
c
c OUTPUT:
c     RLTcrit: R/L_Tcrit for ITG mode
c     RLTcritz: R/L_Tcrit for carbon branch
c     chi_0: normalized chi (ignore)
c     g: L_Tc/L_T, where L_Tc is the critical temperature gradient
c        scale length for the deuterium branch of the ITG mode.
c     chi_i: Anomalous ion thermal diffusivity from toroidal ITG mode.
c     chi_e: Anomalous electron thermal diffusivity from toroidal ITG mode.
c
c     This parameterization of chi is not complete.  There are significant
c     neglected physical processes that are known to be important in
c     many operational regimes.
c
c     The most significant problems are:
c
c     (1) Trapped ion/long wavelength ITG modes.  These modes are known
c     to be unstable for typical edge tokamak parameters.  However, until
c     we have nonlinear estimates of the associated thermal diffusivity,
c     these modes are ignored, leading to overly optimistic predictions of
c     edge thermal confinement.
c     (2) Trapped electron modes, which can alter the stability boundary
c     significantly for low collisionality.  At high collisionality these
c     modes are generally stable and thus largely irrelevant.  When present,
c     they are associated most strongly with particle transport, although
c     there is also an associated heat transport.
c     (3) Minority ion density gradients, which can strongly change chi
c     and LT_crit.
c     (4) Sheared flows, which are stabilizing.  This includes diamagnetic
c     and ExB shear flows.
c     (5) Finite beta effects, generally stabilizing.
c
c jek 1/21/97 added gamma to argument list and compute rot shear
c             outside of routine (gperp = 0)
c
      implicit none
c
      integer itest_chi
c
c     replaced REAL with REAL*8 in the following HSJ 2/11/96
c
      real*8 RLT, RLN, RLNe, shat, zth, tau, rmajor, chi_i, q,
     .       nbeam, taub, rho_i, v_ti, eps, chi_e, nu, chi_0, g,
     .       g_perp, gnu
      real*8 gamma, rot
      real*8 RLTcrit, RLTcritz, chi0, a_0, b_0
      real*8 f_0, f_z, chiz, g_facz, c_0
      real*8 c1, trln, trlne, tshat, kappa
      real*8 chie1, chie2
      real*8 g_fac1
      real*8 zero, sixth, fourth, half, one, two, six  ! for portability
c
      logical opened
      data    opened  /.false./
      integer counter, n54
****  data    counter /0/
      data    n54 /0/
c
      data a_0/ 0.0/
      data b_0/ 0.0/
      data c_0/ 1.0/
****  data a_0/-1.0/
****  data b_0/-1.0/
c
      save a_0, b_0, c_0 ! not needed since we always compile "-static"
c
      zero   = 0.0
      sixth  = 1.0 / 6.0
      fourth = 0.25
      half   = 0.5
      one    = 1.0
      two    = 2.0
      six    = 6.0
c
      if (c_0 .eq. -1) then
        write (*, *) 'C_0 multiplier'
        read  (*, *)  c_0
      end if
c
      if (b_0 .eq. -1) then
        write (*, *) 'beta=?'
        read  (*, *)  b_0
      end if
      taub  = tau / (1.0-nbeam)
      nu    = gnu * 0.84
c
      tRLN  = MIN (ABS (RLN ), six) * SIGN (one, RLN)
      tRLNe = MIN (ABS (RLNe), six)
c
c Formula is not applicable for shat<0.5; doesn't matter most places.
c
      tshat = MAX (shat, half)
c
c     critical ion temperature gradient:
c
      RLTcrit = 2.46*(1.+2.78/q**2)**0.26*(zth/2.)**0.7*taub**0.52
     .          *( (0.671+0.570*tshat-0.189*tRLN)**2
     .          +0.335*tRLN+0.392-0.779*tshat+0.210*tshat**2)
     .          *( 1.-0.942*(2.95*eps**1.257/nu**0.235-0.2126)
     .          *zth**0.516 / ABS (tshat)**0.671)
c
      RLTcritz = 0.75 * (1.0+taub) * (1.0+tshat)
     .                * MAX (one, 3.0 - 2.0 * tRLNe / 3.0)
     .                * (1.0 + 6.0 * MAX (zero, 2.9-zth))
c
      c1 = 1.0
      if (zth .gt. 3.0)  c1 = (3.0/zth)**1.8
c
      f_0 = 11.8*c1*q**1.13/(1.+tshat**0.84)/taub**1.07
     .      *(1.+6.72*eps/nu**0.26/q**0.96)
     .      /(1.+((kappa-1.)*q/3.6)**2)
c
      f_z = 7.88/(1.+tshat) * MAX (fourth, zth-3) / taub**0.8
     .      /(1.+((kappa-1.)*q/3.6)**2)
c
      chi0 = f_0 * rho_i**2 * v_ti / rmajor
      chiz = f_z * rho_i**2 * v_ti / rmajor
c
      if (RLT-RLTcrit .gt. 0.0) then
        g_fac1 = MIN ((RLT-RLTcrit)**0.5, (RLT-RLTcrit))
      else
        g_fac1 = 0.0
      end if
c
      if (RLT-RLTcritz.gt.0.) then
        g_facz = MIN ((RLT-RLTcritz)**0.5, (RLT-RLTcritz))
      else
        g_facz = 0.0
      end if
c
      chi_i = MAX (chi0*g_fac1, chiz*g_facz)
c
      g     = RLT / RLTcrit
      chi_0 = chi0 * SQRT (ABS (RLTcrit))
c
      chie1 = chi0*g_fac1*1.44*tau**0.4*(q/tshat)**0.3*nu**0.14
     .      * MAX (sixth, eps)
      chie1 = 0.5 * chie1*(1.+tRLNe/3.0)
c
      chie2 = 0.5 * MAX (two, (1.0+RLNe/3.0))
      chie2 = chie2*0.526*tau*nu**0.22
      chie2 = chie2*chiz*g_facz
c
c Correction for n_i/n_e and ratio of heat fluxes rather than chi_s:
c
      chi_e = MAX (chie1, chie2) * (7.0-zth)/6.0
c
c     **Preliminary** model of rotational stabilization.  Based on
c     work described in Waltz, et al., Phys. of Plasmas, Vol. 1,
c     p. 3138 (1992).  Here, the ExB shearing rate is denoted
c     as g_perp and the linear growth rate gamma_l is parameterized as:
c
      gamma = 0.25/(1+0.5*tshat**2)/tau
     .      * (1 + 3.0 * MAX (zero, eps-0.16667))
     .      / (1 + MAX (zero, q-3)/15.0)
     .      * (RLT-RLTcrit)*v_ti / rmajor
      gamma = MAX (zero, gamma) ! GS 12/18/98
      rot   = MIN (one, MAX (zero, (one - ABS (g_perp)/(gamma+1.e-20))))
      chi_i = MAX (zero, c_0*chi_i*rot)
      chi_e = MAX (zero, c_0*chi_e*rot)
c
c test
c
      if (itest_chi .eq. 1) then
        chi_i = 1.0
        chi_e = 1.0
      end if
c
      counter = counter + 1
      if (counter .ge. 1000000000) then    ! formerly 2951 for debugging
        write (6, '(a, i7)') ' IP_CHI2 COUNTER =', counter
        if (.not. opened) then
          call getioun(n54,54)
          open  (unit = n54, file = 'ip_chi2.args', status = 'UNKNOWN')
          opened = .true.
        end if
        write   (unit = n54, fmt = '(i3 / 3(6e13.5/), 4e13.5)')
     .           itest_chi,
     .                 rlt,   rln,  rlne,       q,    kappa,   shat,
     .                 zth, nbeam,   tau,     eps,      gnu, g_perp,
     .              rmajor, rho_i,  v_ti, rltcrit, rltcritz,  chi_0,
     .                   g, gamma, chi_i,   chi_e
****    call STOP ('subroutine IP_CHI2: COUNTER >= 2951', 996)!
      end if
c
      return
c
      end

      subroutine qzhes (nm, n, a, b, matz, z)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      logical matz
      integer i, j, k, l, n, lb, l1, nm, nk1, nm1, nm2
      real*8  a(nm,n), b(nm,n), z(nm,n)
      real*8  r, s, t, u1, u2, v1, v2, rho
c
c     this subroutine is the first step of the qz algorithm
c     for solving generalized matrix eigenvalue problems,
c     siam j. numer. anal. 10, 241-256(1973) by moler and stewart.
c
c     this subroutine accepts a pair of real general matrices and
c     reduces one of them to upper hessenberg form and the other
c     to upper triangular form using orthogonal transformations.
c     it is usually followed by qzit, qzval and, possibly, qzvec.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrices.
c
c        a contains a real general matrix.
c
c        b contains a real general matrix.
c
c        matz should be set to .true. if the right hand transformations
c          are to be accumulated for later use in computing
c          eigenvectors, and to .false. otherwise.
c
c     on output
c
c        a has been reduced to upper hessenberg form.  the elements
c          below the first subdiagonal have been set to zero.
c
c        b has been reduced to upper triangular form.  the elements
c          below the main diagonal have been set to zero.
c
c        z contains the product of the right hand transformations if
c          matz has been set to .true.  otherwise, z is not referenced.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
c     .......... initialize z ..........
c
      if (.not. matz)  go to 10
c
      do 3 j = 1, n
c
         do 2 i = 1, n
            z(i,j) = 0.0
    2    continue
c
         z(j,j) = 1.0
    3 continue
c
c     .......... reduce b to upper triangular form ..........
c
   10 if (n .le. 1)  go to 170
      nm1 = n - 1
c
      do 100 l = 1, nm1
         l1 = l + 1
         s = 0.0
c
         do 20 i = l1, n
            s = s + ABS (b(i,l))
   20    continue
c
         if (s .eq. 0.0)  go to 100
         s = s + ABS (b(l,l))
         r = 0.0
c
         do 25 i = l, n
            b(i,l) = b(i,l) / s
            r = r + b(i,l)**2
   25    continue
c
         r      = SIGN (SQRT (r),b(l,l))
         b(l,l) = b(l,l) + r
         rho    = r * b(l,l)
c
         do 50 j = l1, n
            t = 0.0
c
            do 30 i = l, n
               t = t + b(i,l) * b(i,j)
   30       continue
c
            t = -t / rho
c
            do 40 i = l, n
               b(i,j) = b(i,j) + t * b(i,l)
   40       continue
c
   50    continue
c
         do 80 j = 1, n
            t = 0.0
c
            do 60 i = l, n
               t = t + b(i,l) * a(i,j)
   60       continue
c
            t = -t / rho
c
            do 70 i = l, n
               a(i,j) = a(i,j) + t * b(i,l)
   70       continue
c
   80    continue
c
         b(l,l) = -s * r
c
         do 90 i = l1, n
            b(i,l) = 0.0
   90    continue
c
  100 continue
c
c     .......... reduce a to upper hessenberg form, while
c                keeping b triangular ..........
c
      if (n .eq. 2)  go to 170
      nm2 = n - 2
c
      do 160 k = 1, nm2
         nk1 = nm1 - k
c
c     .......... for l=n-1 step -1 until k+1 do -- ..........
c
         do 150 lb = 1, nk1
            l = n - lb
            l1 = l + 1
c
c     .......... zero a(l+1,k) ..........
c
            s  = ABS (a(l,k)) + ABS (a(l1,k))
            if (s .eq. 0.0)  go to 150
            u1 = a(l,k) / s
            u2 = a(l1,k) / s
            r  = SIGN (SQRT (u1*u1+u2*u2),u1)
            v1 =  -(u1 + r) / r
            v2 = -u2 / r
            u2 = v2 / v1
c
            do 110 j = k, n
               t = a(l,j) + u2 * a(l1,j)
               a(l,j) = a(l,j) + t * v1
               a(l1,j) = a(l1,j) + t * v2
  110       continue
c
            a(l1,k) = 0.0
c
            do 120 j = l, n
               t = b(l,j) + u2 * b(l1,j)
               b(l,j) = b(l,j) + t * v1
               b(l1,j) = b(l1,j) + t * v2
  120       continue
c
c     .......... zero b(l+1,l) ..........
c
            s  = ABS (b(l1,l1)) + ABS (b(l1,l))
            if (s .eq. 0.0)  go to 150
            u1 = b(l1,l1) / s
            u2 = b(l1,l) / s
            r  = SIGN (SQRT (u1*u1+u2*u2),u1)
            v1 =  -(u1 + r) / r
            v2 = -u2 / r
            u2 = v2 / v1
c
            do 130 i = 1, l1
               t = b(i,l1) + u2 * b(i,l)
               b(i,l1) = b(i,l1) + t * v1
               b(i,l) = b(i,l) + t * v2
  130       continue
c
            b(l1,l) = 0.0
c
            do 140 i = 1, n
               t = a(i,l1) + u2 * a(i,l)
               a(i,l1) = a(i,l1) + t * v1
               a(i,l) = a(i,l) + t * v2
  140       continue
c
            if (.not. matz)  go to 150
c
            do 145 i = 1, n
               t = z(i,l1) + u2 * z(i,l)
               z(i,l1) = z(i,l1) + t * v1
               z(i,l) = z(i,l) + t * v2
  145       continue
c
  150    continue
c
  160 continue
c
  170 return
c
      end

      subroutine qzit (nm, n, a, b, eps1, matz, z, ierr)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      logical matz, notlas
      integer i,j,k,l,n,en,k1,k2,ld,ll,l1,na,nm,ish,itn,its,km1,lm1,
     .        enm2,ierr,lor1,enorn
      real*8  a(nm,n),b(nm,n),z(nm,n)
      real*8  r,s,t,a1,a2,a3,ep,sh,u1,u2,u3,v1,v2,v3,ani,a11,
     .        a12,a21,a22,a33,a34,a43,a44,bni,b11,b12,b22,b33,b34,
     .        b44,epsa,epsb,eps1,anorm,bnorm,epslon
c
c     this subroutine is the second step of the qz algorithm
c     for solving generalized matrix eigenvalue problems,
c     siam j. numer. anal. 10, 241-256(1973) by moler and stewart,
c     as modified in technical note nasa tn d-7305(1973) by ward.
c
c     this subroutine accepts a pair of real matrices, one of them
c     in upper hessenberg form and the other in upper triangular form.
c     it reduces the hessenberg matrix to quasi-triangular form using
c     orthogonal transformations while maintaining the triangular form
c     of the other matrix.  it is usually preceded by  qzhes  and
c     followed by  qzval  and, possibly,  qzvec.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrices.
c
c        a contains a real upper hessenberg matrix.
c
c        b contains a real upper triangular matrix.
c
c        eps1 is a tolerance used to determine negligible elements.
c          eps1 = 0.0 (or negative) may be input, in which case an
c          element will be neglected only if it is less than roundoff
c          error times the norm of its matrix.  if the input eps1 is
c          positive, then an element will be considered negligible
c          if it is less than eps1 times the norm of its matrix.  a
c          positive value of eps1 may result in faster execution,
c          but less accurate results.
c
c        matz should be set to .true. if the right hand transformations
c          are to be accumulated for later use in computing
c          eigenvectors, and to .false. otherwise.
c
c        z contains, if matz has been set to .true., the
c          transformation matrix produced in the reduction
c          by  qzhes, if performed, or else the identity matrix.
c          if matz has been set to .false., z is not referenced.
c
c     on output
c
c        a has been reduced to quasi-triangular form.  the elements
c          below the first subdiagonal are still zero and no two
c          consecutive subdiagonal elements are nonzero.
c
c        b is still in upper triangular form, although its elements
c          have been altered.  the location b(n,1) is used to store
c          eps1 times the norm of b for later use by  qzval  and  qzvec.
c
c        z contains the product of the right hand transformations
c          (for both steps) if matz has been set to .true..
c
c        ierr is set to
c          zero       for normal return,
c          j          if the limit of 30*n iterations is exhausted
c                     while the j-th eigenvalue is being sought.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      ierr = 0
c
c     .......... compute epsa,epsb ..........
c
      anorm = 0.0
      bnorm = 0.0
c
      do 30 i = 1, n
         ani = 0.0
         if (i .ne. 1) ani = ABS (a(i,i-1))
         bni = 0.0
c
         do 20 j = i, n
            ani = ani + ABS (a(i,j))
            bni = bni + ABS (b(i,j))
   20    continue
c
         if (ani .gt. anorm) anorm = ani
         if (bni .gt. bnorm) bnorm = bni
   30 continue
c
      if (anorm .eq. 0.0) anorm = 1.0
      if (bnorm .eq. 0.0) bnorm = 1.0
      ep = eps1
      if (ep .gt. 0.0)  go to 50
c
c     .......... use roundoff level if eps1 is zero ..........
c
      ep   = epslon(1.0d0)
   50 epsa = ep * anorm
      epsb = ep * bnorm
c
c     .......... reduce a to quasi-triangular form, while
c                keeping b triangular ..........
c
      lor1 = 1
      enorn = n
      en = n
      itn = 30*n
c
c     .......... begin qz step ..........
c
   60 if (en .le. 2 )  go to 1001
      if (.not. matz)  enorn = en
      its  = 0
      na   = en - 1
      enm2 = na - 1
   70 ish  = 2
c
c     .......... check for convergence or reducibility.
c                for l=en step -1 until 1 do -- ..........
c
      do 80 ll = 1, en
         lm1 = en - ll
         l   = lm1 + 1
         if (l .eq. 1)  go to 95
         if (ABS (a(l,lm1)) .le. epsa)  go to 90
   80 continue
c
   90 a(l,lm1) = 0.0
      if (l .lt. na)  go to 95
c
c     .......... 1-by-1 or 2-by-2 block isolated ..........
c
      en = lm1
      go to 60
c
c     .......... check for small top of b ..........
c
   95 ld     = l
  100 l1     = l + 1
      b11    = b(l,l)
      if (ABS (b11) .gt. epsb)  go to 120
      b(l,l) = 0.0
      s      = ABS (a(l,l)) + ABS (a(l1,l))
      u1     = a(l,l) / s
      u2     = a(l1,l) / s
      r      = SIGN (SQRT (u1*u1+u2*u2),u1)
      v1     = -(u1 + r) / r
      v2     = -u2 / r
      u2     = v2 / v1
c
      do 110 j=l,enorn
         t       = a(l,j) + u2 * a(l1,j)
         a(l,j)  = a(l,j) + t * v1
         a(l1,j) = a(l1,j) + t * v2
         t       = b(l,j) + u2 * b(l1,j)
         b(l,j)  = b(l,j) + t * v1
         b(l1,j) = b(l1,j) + t * v2
  110 continue
c
      if (l .ne. 1)  a(l,lm1) = -a(l,lm1)
      lm1 = l
      l   = l1
      go to 90
  120 a11 = a(l,l) / b11
      a21 = a(l1,l) / b11
      if (ish .eq. 1)  go to 140
c
c     .......... iteration strategy ..........
c
      if (itn .eq. 0)  go to 1000
      if (its .eq. 10)  go to 155
c
c     .......... determine type of shift ..........
c
      b22 = b(l1,l1)
      if (ABS (b22) .lt. epsb) b22 = epsb
      b33 = b(na,na)
      if (ABS (b33) .lt. epsb) b33 = epsb
      b44 = b(en,en)
      if (ABS (b44) .lt. epsb) b44 = epsb
      a33 = a(na,na) / b33
      a34 = a(na,en) / b44
      a43 = a(en,na) / b33
      a44 = a(en,en) / b44
      b34 = b(na,en) / b44
      t = 0.5 * (a43 * b34 - a33 - a44)
      r = t * t + a34 * a43 - a33 * a44
      if (r .lt. 0.0)  go to 150
c
c     .......... determine single shift zeroth column of a ..........
c
      ish = 1
      r = SQRT (r)
      sh = -t + r
      s = -t - r
      if (ABS (s-a44) .lt. ABS (sh-a44)) sh = s
c
c     .......... look for two consecutive small
c                sub-diagonal elements of a.
c                for l=en-2 step -1 until ld do -- ..........
c
      do 130 ll = ld, enm2
         l = enm2 + ld - ll
         if (l .eq. ld)  go to 140
         lm1 = l - 1
         l1 = l + 1
         t = a(l,l)
         if (ABS (b(l,l)) .gt. epsb) t = t - sh * b(l,l)
         if (ABS (a(l,lm1)) .le. ABS (t/a(l1,l)) * epsa)  go to 100
  130 continue
c
  140 a1 = a11 - sh
      a2 = a21
      if (l .ne. ld) a(l,lm1) = -a(l,lm1)
      go to 160
c
c     .......... determine double shift zeroth column of a ..........
c
  150 a12 = a(l,l1) / b22
      a22 = a(l1,l1) / b22
      b12 = b(l,l1) / b22
      a1 = ((a33 - a11) * (a44 - a11) - a34 * a43 + a43 * b34 * a11)
     .     / a21 + a12 - a11 * b12
      a2 = (a22 - a11) - a21 * b12 - (a33 - a11) - (a44 - a11)
     .     + a43 * b34
      a3 = a(l1+1,l1) / b22
      go to 160
c
c     .......... ad hoc shift ..........
c
  155 a1 = 0.0
      a2 = 1.0
      a3 = 1.1605
  160 its = its + 1
      itn = itn - 1
      if (.not. matz) lor1 = ld
c
c     .......... main loop ..........
c
      do 260 k = l, na
         notlas = k .ne. na .and. ish .eq. 2
         k1 = k + 1
         k2 = k + 2
         km1 = max0(k-1,l)
         ll = min0(en,k1+ish)
         if (notlas)  go to 190
c
c     .......... zero a(k+1,k-1) ..........
c
         if (k .eq. l)  go to 170
         a1 = a(k,km1)
         a2 = a(k1,km1)
  170    s = ABS (a1) + ABS (a2)
         if (s .eq. 0.0)  go to 70
         u1 = a1 / s
         u2 = a2 / s
         r = SIGN (SQRT (u1*u1+u2*u2),u1)
         v1 = -(u1 + r) / r
         v2 = -u2 / r
         u2 = v2 / v1
c
         do 180 j = km1, enorn
            t = a(k,j) + u2 * a(k1,j)
            a(k,j) = a(k,j) + t * v1
            a(k1,j) = a(k1,j) + t * v2
            t = b(k,j) + u2 * b(k1,j)
            b(k,j) = b(k,j) + t * v1
            b(k1,j) = b(k1,j) + t * v2
  180    continue
c
         if (k .ne. l) a(k1,km1) = 0.0
         go to 240
c
c     .......... zero a(k+1,k-1) and a(k+2,k-1) ..........
c
  190    if (k .eq. l)  go to 200
         a1 = a(k,km1)
         a2 = a(k1,km1)
         a3 = a(k2,km1)
  200    s = ABS (a1) + ABS (a2) + ABS (a3)
         if (s .eq. 0.0)  go to 260
         u1 = a1 / s
         u2 = a2 / s
         u3 = a3 / s
         r = SIGN (SQRT (u1*u1+u2*u2+u3*u3),u1)
         v1 = -(u1 + r) / r
         v2 = -u2 / r
         v3 = -u3 / r
         u2 = v2 / v1
         u3 = v3 / v1
c
         do 210 j = km1, enorn
            t = a(k,j) + u2 * a(k1,j) + u3 * a(k2,j)
            a(k,j) = a(k,j) + t * v1
            a(k1,j) = a(k1,j) + t * v2
            a(k2,j) = a(k2,j) + t * v3
            t = b(k,j) + u2 * b(k1,j) + u3 * b(k2,j)
            b(k,j) = b(k,j) + t * v1
            b(k1,j) = b(k1,j) + t * v2
            b(k2,j) = b(k2,j) + t * v3
  210    continue
c
         if (k .eq. l)  go to 220
         a(k1,km1) = 0.0
         a(k2,km1) = 0.0
c
c     .......... zero b(k+2,k+1) and b(k+2,k) ..........
c
  220    s = ABS (b(k2,k2)) + ABS (b(k2,k1)) + ABS (b(k2,k))
         if (s .eq. 0.0)  go to 240
         u1 = b(k2,k2) / s
         u2 = b(k2,k1) / s
         u3 = b(k2,k) / s
         r = SIGN (SQRT (u1*u1+u2*u2+u3*u3),u1)
         v1 = -(u1 + r) / r
         v2 = -u2 / r
         v3 = -u3 / r
         u2 = v2 / v1
         u3 = v3 / v1
c
         do 230 i = lor1, ll
            t = a(i,k2) + u2 * a(i,k1) + u3 * a(i,k)
            a(i,k2) = a(i,k2) + t * v1
            a(i,k1) = a(i,k1) + t * v2
            a(i,k) = a(i,k) + t * v3
            t = b(i,k2) + u2 * b(i,k1) + u3 * b(i,k)
            b(i,k2) = b(i,k2) + t * v1
            b(i,k1) = b(i,k1) + t * v2
            b(i,k) = b(i,k) + t * v3
  230    continue
c
         b(k2,k) = 0.0
         b(k2,k1) = 0.0
         if (.not. matz)  go to 240
c
         do 235 i = 1, n
            t = z(i,k2) + u2 * z(i,k1) + u3 * z(i,k)
            z(i,k2) = z(i,k2) + t * v1
            z(i,k1) = z(i,k1) + t * v2
            z(i,k) = z(i,k) + t * v3
  235    continue
c
c     .......... zero b(k+1,k) ..........
c
  240    s = ABS (b(k1,k1)) + ABS (b(k1,k))
         if (s .eq. 0.0)  go to 260
         u1 = b(k1,k1) / s
         u2 = b(k1,k) / s
         r = SIGN (SQRT (u1*u1+u2*u2),u1)
         v1 = -(u1 + r) / r
         v2 = -u2 / r
         u2 = v2 / v1
c
         do 250 i = lor1, ll
            t = a(i,k1) + u2 * a(i,k)
            a(i,k1) = a(i,k1) + t * v1
            a(i,k) = a(i,k) + t * v2
            t = b(i,k1) + u2 * b(i,k)
            b(i,k1) = b(i,k1) + t * v1
            b(i,k) = b(i,k) + t * v2
  250    continue
c
         b(k1,k) = 0.0
         if (.not. matz)  go to 260
c
         do 255 i = 1, n
            t = z(i,k1) + u2 * z(i,k)
            z(i,k1) = z(i,k1) + t * v1
            z(i,k) = z(i,k) + t * v2
  255    continue
c
  260 continue
c
c     .......... end qz step ..........
c
      go to 70
c
c     .......... set error -- all eigenvalues have not
c                converged after 30*n iterations ..........
c
 1000 ierr = en
c
c     .......... save epsb for use by qzval and qzvec ..........
c
 1001 if (n .gt. 1) b(n,1) = epsb
      return
c
      end

      subroutine qzval (nm, n, a, b, alfr, alfi, beta, matz, z)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      logical matz
      integer i, j, n, en, na, nm, nn, isw
      real*8  a(nm,n), b(nm,n), alfr(n), alfi(n), beta(n), z(nm,n)
      real*8  c,d,e,r,s,t,an,a1,a2,bn,cq,cz,di,dr,ei,ti,tr,u1,
     .        u2,v1,v2,a1i,a11,a12,a2i,a21,a22,b11,b12,b22,sqi,sqr,
     .        ssi,ssr,szi,szr,a11i,a11r,a12i,a12r,a22i,a22r,epsb
c
c     this subroutine is the third step of the qz algorithm
c     for solving generalized matrix eigenvalue problems,
c     siam j. numer. anal. 10, 241-256(1973) by moler and stewart.
c
c     this subroutine accepts a pair of real matrices, one of them
c     in quasi-triangular form and the other in upper triangular form.
c     it reduces the quasi-triangular matrix further, so that any
c     remaining 2-by-2 blocks correspond to pairs of complex
c     eigenvalues, and returns quantities whose ratios give the
c     generalized eigenvalues.  it is usually preceded by  qzhes
c     and  qzit  and may be followed by  qzvec.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrices.
c
c        a contains a real upper quasi-triangular matrix.
c
c        b contains a real upper triangular matrix.  in addition,
c          location b(n,1) contains the tolerance quantity (epsb)
c          computed and saved in  qzit.
c
c        matz should be set to .true. if the right hand transformations
c          are to be accumulated for later use in computing
c          eigenvectors, and to .false. otherwise.
c
c        z contains, if matz has been set to .true., the
c          transformation matrix produced in the reductions by qzhes
c          and qzit, if performed, or else the identity matrix.
c          if matz has been set to .false., z is not referenced.
c
c     on output
c
c        a has been reduced further to a quasi-triangular matrix
c          in which all nonzero subdiagonal elements correspond to
c          pairs of complex eigenvalues.
c
c        b is still in upper triangular form, although its elements
c          have been altered.  b(n,1) is unaltered.
c
c        alfr and alfi contain the real and imaginary parts of the
c          diagonal elements of the triangular matrix that would be
c          obtained if a were reduced completely to triangular form
c          by unitary transformations.  non-zero values of alfi occur
c          in pairs, the first member positive and the second negative.
c
c        beta contains the diagonal elements of the corresponding b,
c          normalized to be real and non-negative.  the generalized
c          eigenvalues are then the ratios ((alfr+i*alfi)/beta).
c
c        z contains the product of the right hand transformations
c          (for all three steps) if matz has been set to .true.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      epsb = b(n,1)
      isw  = 1
c
c     .......... find eigenvalues of quasi-triangular matrices.
c                for en=n step -1 until 1 do -- ..........
c
      do 510 nn = 1, n
         en = n + 1 - nn
         na = en - 1
         if (isw .eq. 2)  go to 505
         if (en .eq. 1)  go to 410
         if (a(en,na) .ne. 0.0)  go to 420
c
c     .......... 1-by-1 block, one real root ..........
c
  410    alfr(en) = a(en,en)
         if (b(en,en) .lt. 0.0) alfr(en) = -alfr(en)
         beta(en) = ABS (b(en,en))
         alfi(en) = 0.0
         go to 510
c
c     .......... 2-by-2 block ..........
c
  420    if (ABS (b(na,na)) .le. epsb)  go to 455
         if (ABS (b(en,en)) .gt. epsb)  go to 430
         a1 = a(en,en)
         a2 = a(en,na)
         bn = 0.0
         go to 435
  430    an = ABS (a(na,na)) + ABS (a(na,en)) + ABS (a(en,na))
     .      + ABS (a(en,en))
         bn = ABS (b(na,na)) + ABS (b(na,en)) + ABS (b(en,en))
         a11 = a(na,na) / an
         a12 = a(na,en) / an
         a21 = a(en,na) / an
         a22 = a(en,en) / an
         b11 = b(na,na) / bn
         b12 = b(na,en) / bn
         b22 = b(en,en) / bn
         e = a11 / b11
         ei = a22 / b22
         s = a21 / (b11 * b22)
         t = (a22 - e * b22) / b22
         if (ABS (e) .le. ABS (ei))  go to 431
         e = ei
         t = (a11 - e * b11) / b11
  431    c = 0.5 * (t - s * b12)
         d = c * c + s * (a12 - e * b12)
         if (d .lt. 0.0)  go to 480
c
c     .......... two real roots.
c                zero both a(en,na) and b(en,na) ..........
c
         e = e + (c + SIGN (SQRT (d),c))
         a11 = a11 - e * b11
         a12 = a12 - e * b12
         a22 = a22 - e * b22
         if (ABS (a11) + ABS (a12) .lt.
     .       ABS (a21) + ABS (a22))  go to 432
         a1 = a12
         a2 = a11
         go to 435
  432    a1 = a22
         a2 = a21
c
c     .......... choose and apply real z ..........
c
  435    s = ABS (a1) + ABS (a2)
         u1 = a1 / s
         u2 = a2 / s
         r = SIGN (SQRT (u1*u1+u2*u2),u1)
         v1 = -(u1 + r) / r
         v2 = -u2 / r
         u2 = v2 / v1
c
         do 440 i = 1, en
            t = a(i,en) + u2 * a(i,na)
            a(i,en) = a(i,en) + t * v1
            a(i,na) = a(i,na) + t * v2
            t = b(i,en) + u2 * b(i,na)
            b(i,en) = b(i,en) + t * v1
            b(i,na) = b(i,na) + t * v2
  440    continue
c
         if (.not. matz)  go to 450
c
         do 445 i = 1, n
            t = z(i,en) + u2 * z(i,na)
            z(i,en) = z(i,en) + t * v1
            z(i,na) = z(i,na) + t * v2
  445    continue
c
  450    if (bn .eq. 0.0)  go to 475
         if (an .lt. ABS (e) * bn)  go to 455
         a1 = b(na,na)
         a2 = b(en,na)
         go to 460
  455    a1 = a(na,na)
         a2 = a(en,na)
c
c     .......... choose and apply real q ..........
c
  460    s = ABS (a1) + ABS (a2)
         if (s .eq. 0.0)  go to 475
         u1 = a1 / s
         u2 = a2 / s
         r = SIGN (SQRT (u1*u1+u2*u2),u1)
         v1 = -(u1 + r) / r
         v2 = -u2 / r
         u2 = v2 / v1
c
         do 470 j = na, n
            t = a(na,j) + u2 * a(en,j)
            a(na,j) = a(na,j) + t * v1
            a(en,j) = a(en,j) + t * v2
            t = b(na,j) + u2 * b(en,j)
            b(na,j) = b(na,j) + t * v1
            b(en,j) = b(en,j) + t * v2
  470    continue
c
  475    a(en,na) = 0.0
         b(en,na) = 0.0
         alfr(na) = a(na,na)
         alfr(en) = a(en,en)
         if (b(na,na) .lt. 0.0) alfr(na) = -alfr(na)
         if (b(en,en) .lt. 0.0) alfr(en) = -alfr(en)
         beta(na) = ABS (b(na,na))
         beta(en) = ABS (b(en,en))
         alfi(en) = 0.0
         alfi(na) = 0.0
         go to 505
c
c     .......... two complex roots ..........
c
  480    e = e + c
         ei = SQRT (-d)
         a11r = a11 - e * b11
         a11i = ei * b11
         a12r = a12 - e * b12
         a12i = ei * b12
         a22r = a22 - e * b22
         a22i = ei * b22
         if (ABS (a11r) + ABS (a11i) + ABS (a12r) + ABS (a12i) .lt.
     .       ABS (a21) + ABS (a22r) + ABS (a22i))  go to 482
         a1 = a12r
         a1i = a12i
         a2 = -a11r
         a2i = -a11i
         go to 485
  482    a1 = a22r
         a1i = a22i
         a2 = -a21
         a2i = 0.0
c
c     .......... choose complex z ..........
c
  485    cz = SQRT (a1*a1+a1i*a1i)
         if (cz .eq. 0.0)  go to 487
         szr = (a1 * a2 + a1i * a2i) / cz
         szi = (a1 * a2i - a1i * a2) / cz
         r = SQRT (cz*cz+szr*szr+szi*szi)
         cz = cz / r
         szr = szr / r
         szi = szi / r
         go to 490
  487    szr = 1.0
         szi = 0.0
  490    if (an .lt. (ABS (e) + ei) * bn)  go to 492
         a1 = cz * b11 + szr * b12
         a1i = szi * b12
         a2 = szr * b22
         a2i = szi * b22
         go to 495
  492    a1 = cz * a11 + szr * a12
         a1i = szi * a12
         a2 = cz * a21 + szr * a22
         a2i = szi * a22
c
c     .......... choose complex q ..........
c
  495    cq = SQRT (a1*a1+a1i*a1i)
         if (cq .eq. 0.0)  go to 497
         sqr = (a1 * a2 + a1i * a2i) / cq
         sqi = (a1 * a2i - a1i * a2) / cq
         r = SQRT (cq*cq+sqr*sqr+sqi*sqi)
         cq = cq / r
         sqr = sqr / r
         sqi = sqi / r
         go to 500
  497    sqr = 1.0
         sqi = 0.0
c
c     .......... compute diagonal elements that would result
c                if transformations were applied ..........
c
  500    ssr = sqr * szr + sqi * szi
         ssi = sqr * szi - sqi * szr
         i = 1
         tr = cq * cz * a11 + cq * szr * a12 + sqr * cz * a21
     .      + ssr * a22
         ti = cq * szi * a12 - sqi * cz * a21 + ssi * a22
         dr = cq * cz * b11 + cq * szr * b12 + ssr * b22
         di = cq * szi * b12 + ssi * b22
         go to 503
  502    i = 2
         tr = ssr * a11 - sqr * cz * a12 - cq * szr * a21
     .      + cq * cz * a22
         ti = -ssi * a11 - sqi * cz * a12 + cq * szi * a21
         dr = ssr * b11 - sqr * cz * b12 + cq * cz * b22
         di = -ssi * b11 - sqi * cz * b12
  503    t = ti * dr - tr * di
         j = na
         if (t .lt. 0.0) j = en
         r = SQRT (dr*dr+di*di)
         beta(j) = bn * r
         alfr(j) = an * (tr * dr + ti * di) / r
         alfi(j) = an * t / r
         if (i .eq. 1)  go to 502
  505    isw = 3 - isw
  510 continue
      b(n,1) = epsb
c
      return
c
      end

      subroutine qzvec (nm, n, a, b, alfr, alfi, beta, z)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      integer i, j, k, m, n, en, ii, jj, na, nm, nn, isw, enm2
      real*8  a(nm,n), b(nm,n), alfr(n), alfi(n), beta(n), z(nm,n)
      real*8  d,q,r,s,t,w,x,y,di,dr,ra,rr,sa,ti,tr,t1,t2,w1,x1,
     .        zz,z1,alfm,almi,almr,betm,epsb
c
c     this subroutine is the optional fourth step of the qz algorithm
c     for solving generalized matrix eigenvalue problems,
c     siam j. numer. anal. 10, 241-256(1973) by moler and stewart.
c
c     this subroutine accepts a pair of real matrices, one of them in
c     quasi-triangular form (in which each 2-by-2 block corresponds to
c     a pair of complex eigenvalues) and the other in upper triangular
c     form.  it computes the eigenvectors of the triangular problem and
c     transforms the results back to the original coordinate system.
c     it is usually preceded by  qzhes,  qzit, and  qzval.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrices.
c
c        a contains a real upper quasi-triangular matrix.
c
c        b contains a real upper triangular matrix.  in addition,
c          location b(n,1) contains the tolerance quantity (epsb)
c          computed and saved in  qzit.
c
c        alfr, alfi, and beta  are vectors with components whose
c          ratios ((alfr+i*alfi)/beta) are the generalized
c          eigenvalues.  they are usually obtained from  qzval.
c
c        z contains the transformation matrix produced in the
c          reductions by  qzhes,  qzit, and  qzval, if performed.
c          if the eigenvectors of the triangular problem are
c          desired, z must contain the identity matrix.
c
c     on output
c
c        a is unaltered.  its subdiagonal elements provide information
c           about the storage of the complex eigenvectors.
c
c        b has been destroyed.
c
c        alfr, alfi, and beta are unaltered.
c
c        z contains the real and imaginary parts of the eigenvectors.
c          if alfi(i) .eq. 0.0, the i-th eigenvalue is real and
c            the i-th column of z contains its eigenvector.
c          if alfi(i) .ne. 0.0, the i-th eigenvalue is complex.
c            if alfi(i) .gt. 0.0, the eigenvalue is the first of
c              a complex pair and the i-th and (i+1)-th columns
c              of z contain its eigenvector.
c            if alfi(i) .lt. 0.0, the eigenvalue is the second of
c              a complex pair and the (i-1)-th and i-th columns
c              of z contain the conjugate of its eigenvector.
c          each eigenvector is normalized so that the modulus
c          of its largest component is 1.0 .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      epsb = b(n,1)
      isw = 1
c
c     .......... for en=n step -1 until 1 do -- ..........
c
      do 800 nn = 1, n
         en = n + 1 - nn
         na = en - 1
         if (isw .eq. 2)  go to 795
         if (alfi(en) .ne. 0.0)  go to 710
c
c     .......... real vector ..........
c
         m = en
         b(en,en) = 1.0
         if (na .eq. 0)  go to 800
         alfm = alfr(m)
         betm = beta(m)
c
c     .......... for i=en-1 step -1 until 1 do -- ..........
c
         do 700 ii = 1, na
            i = en - ii
            w = betm * a(i,i) - alfm * b(i,i)
            r = 0.0
c
            do 610 j = m, en
  610       r = r + (betm * a(i,j) - alfm * b(i,j)) * b(j,en)
c
            if (i .eq. 1 .or. isw .eq. 2)  go to 630
            if (betm * a(i,i-1) .eq. 0.0)  go to 630
            zz = w
            s = r
            go to 690
  630       m = i
            if (isw .eq. 2)  go to 640
c
c     .......... real 1-by-1 block ..........
c
            t = w
            if (w .eq. 0.0) t = epsb
            b(i,en) = -r / t
            go to 700
c
c     .......... real 2-by-2 block ..........
c
  640       x = betm * a(i,i+1) - alfm * b(i,i+1)
            y = betm * a(i+1,i)
            q = w * zz - x * y
            t = (x * s - zz * r) / q
            b(i,en) = t
            if (ABS (x) .le. ABS (zz))  go to 650
            b(i+1,en) = (-r - w * t) / x
            go to 690
  650       b(i+1,en) = (-s - y * t) / zz
  690       isw = 3 - isw
  700    continue
c
c     .......... end real vector ..........
c
         go to 800
c
c     .......... complex vector ..........
c
  710    m = na
         almr = alfr(m)
         almi = alfi(m)
         betm = beta(m)
c
c     .......... last vector component chosen imaginary so that
c                eigenvector matrix is triangular ..........
c
         y = betm * a(en,na)
         b(na,na) = -almi * b(en,en) / y
         b(na,en) = (almr * b(en,en) - betm * a(en,en)) / y
         b(en,na) = 0.0
         b(en,en) = 1.0
         enm2 = na - 1
         if (enm2 .eq. 0)  go to 795
c
c     .......... for i=en-2 step -1 until 1 do -- ..........
c
         do 790 ii = 1, enm2
            i = na - ii
            w = betm * a(i,i) - almr * b(i,i)
            w1 = -almi * b(i,i)
            ra = 0.0
            sa = 0.0
c
            do 760 j = m, en
               x = betm * a(i,j) - almr * b(i,j)
               x1 = -almi * b(i,j)
               ra = ra + x * b(j,na) - x1 * b(j,en)
               sa = sa + x * b(j,en) + x1 * b(j,na)
  760       continue
c
            if (i .eq. 1 .or. isw .eq. 2)  go to 770
            if (betm * a(i,i-1) .eq. 0.0)  go to 770
            zz = w
            z1 = w1
            r = ra
            s = sa
            isw = 2
            go to 790
  770       m = i
            if (isw .eq. 2)  go to 780
c
c     .......... complex 1-by-1 block ..........
c
            tr = -ra
            ti = -sa
  773       dr = w
            di = w1
c
c     .......... complex divide (t1,t2) = (tr,ti) / (dr,di) ..........
c
  775       if (ABS (di) .gt. ABS (dr))  go to 777
            rr = di / dr
            d = dr + di * rr
            t1 = (tr + ti * rr) / d
            t2 = (ti - tr * rr) / d
            go to (787,782), isw
  777       rr = dr / di
            d = dr * rr + di
            t1 = (tr * rr + ti) / d
            t2 = (ti * rr - tr) / d
            go to (787,782), isw
c
c     .......... complex 2-by-2 block ..........
c
  780       x = betm * a(i,i+1) - almr * b(i,i+1)
            x1 = -almi * b(i,i+1)
            y = betm * a(i+1,i)
            tr = y * ra - w * r + w1 * s
            ti = y * sa - w * s - w1 * r
            dr = w * zz - w1 * z1 - x * y
            di = w * z1 + w1 * zz - x1 * y
            if (dr .eq. 0.0 .and. di .eq. 0.0) dr = epsb
            go to 775
  782       b(i+1,na) = t1
            b(i+1,en) = t2
            isw = 1
            if (ABS (y) .gt. ABS (w) + ABS (w1))  go to 785
            tr = -ra - x * b(i+1,na) + x1 * b(i+1,en)
            ti = -sa - x * b(i+1,en) - x1 * b(i+1,na)
            go to 773
  785       t1 = (-r - zz * b(i+1,na) + z1 * b(i+1,en)) / y
            t2 = (-s - zz * b(i+1,en) - z1 * b(i+1,na)) / y
  787       b(i,na) = t1
            b(i,en) = t2
  790    continue
c
c     .......... end complex vector ..........
c
  795    isw = 3 - isw
  800 continue
c
c     .......... end back substitution.
c                transform to original coordinate system.
c                for j=n step -1 until 1 do -- ..........
c
      do 880 jj = 1, n
         j = n + 1 - jj
c
         do 880 i = 1, n
            zz = 0.0
c
            do 860 k = 1, j
  860       zz = zz + z(i,k) * b(k,j)
c
            z(i,j) = zz
  880 continue
c
c     .......... normalize so that modulus of largest
c                component of each vector is 1.
c                (isw is 1 initially from before) ..........
c
      do 950 j = 1, n
         d = 0.0
         if (isw .eq. 2)  go to 920
         if (alfi(j) .ne. 0.0)  go to 945
c
         do 890 i = 1, n
            if (ABS (z(i,j)) .gt. d) d = ABS (z(i,j))
  890    continue
c
         do 900 i = 1, n
  900    z(i,j) = z(i,j) / d
c
         go to 950
c
  920    do 930 i = 1, n
            r = ABS (z(i,j-1)) + ABS (z(i,j))
            if (r .ne. 0.0) r = r * SQRT ((z(i,j-1)/r)**2
     .                                     +(z(i,j)/r)**2)
            if (r .gt. d) d = r
  930    continue
c
         do 940 i = 1, n
            z(i,j-1) = z(i,j-1) / d
            z(i,j) = z(i,j) / d
  940    continue
c
  945    isw = 3 - isw
  950 continue
c
      return
c
      end

      subroutine rgg (nm, n, a, b, alfr, alfi, beta, matz, z, ierr)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      logical tf
      integer n, nm, ierr, matz
      real*8  a(nm,n), b(nm,n), alfr(n), alfi(n), beta(n), z(nm,n)
      real*8  zero
c
c     this subroutine calls the recommended sequence of
c     subroutines from the eigensystem subroutine package (EISPACK)
c     to find the eigenvalues and eigenvectors (if desired)
c     for the real general generalized eigenproblem  ax = (lambda)bx.
c
c     on input
c
c        nm  must be set to the row dimension of the two-dimensional
c        array parameters as declared in the calling program
c        dimension statement.
c
c        n  is the order of the matrices  a  and  b.
c
c        a  contains a real general matrix.
c
c        b  contains a real general matrix.
c
c        matz  is an integer variable set equal to zero if
c        only eigenvalues are desired.  otherwise it is set to
c        any non-zero integer for both eigenvalues and eigenvectors.
c
c     on output
c
c        alfr  and  alfi  contain the real and imaginary parts,
c        respectively, of the numerators of the eigenvalues.
c
c        beta  contains the denominators of the eigenvalues,
c        which are thus given by the ratios  (alfr+i*alfi)/beta.
c        complex conjugate pairs of eigenvalues appear consecutively
c        with the eigenvalue having the positive imaginary part first.
c
c        z  contains the real and imaginary parts of the eigenvectors
c        if matz is not zero.  if the j-th eigenvalue is real, the
c        j-th column of  z  contains its eigenvector.  if the j-th
c        eigenvalue is complex with positive imaginary part, the
c        j-th and (j+1)-th columns of  z  contain the real and
c        imaginary parts of its eigenvector.  the conjugate of this
c        vector is the eigenvector for the conjugate eigenvalue.
c
c        ierr  is an integer output variable set equal to an error
c           completion code described in the documentation for qzit.
c           the normal completion code is zero.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      zero = 0.0
      if (n .le. nm)  go to 10
      ierr = 10 * n
      go to 50
c
   10 if (matz .ne. 0)  go to 20
c
c     .......... find eigenvalues only ..........
c
      tf = .false.
      call qzhes (nm, n, a, b, tf, z)
      call  qzit (nm, n, a, b, zero, tf, z, ierr)
      call qzval (nm, n, a, b, alfr, alfi, beta, tf, z)
      go to 50
c
c     .......... find both eigenvalues and eigenvectors ..........
c
   20 tf = .true.
      call qzhes (nm, n, a, b, tf, z)
      call  qzit (nm, n, a, b, zero, tf, z, ierr)
      call qzval (nm, n, a, b, alfr, alfi, beta, tf, z)
      if (ierr .ne. 0)  go to 50
      call qzvec (nm, n, a, b, alfr, alfi, beta, z)
   50 return
c
      end

      subroutine weiland_12 (
     .    arho_exp,        ! minor radius [m]
     .    rmajor_exp,      ! major radius [m]
     .    rho_exp,         ! local r/a
     .    elonga_exp,      ! local plasma elongation
     .    bt_exp,          ! toriodal magnetic field [Gauss]
     .        te_exp,      !     Te [keV]
     .    gradte_exp,      ! gradTe [keV/m]
     .        ti_exp,      !     Ti [keV]
     .    gradti_exp,      ! gradTi [keV/m]
     .        tz_exp,      !     Tz [keV]   main impurity species
     .    gradtz_exp,      ! gradTz [keV/m] main impurity species
     .        ne_exp,      !     ne [*10^19 m-3]
     .    gradne_exp,      ! gradne [*10^19 m-4]
     .        ni_exp,      !     ni [*10^19 m-3]
     .    gradni_exp,      ! gradni [*10^19 m-4]
     .        nz_exp,      !     nz [*10^19 m-3] impurity
     .    gradnz_exp,      ! gradnz [*10^19 m-4] impurity
     .    ftrapein_exp,    ! fraction trapped electrons
     .    azin_exp,        ! mZ/mH impurity mass to H isotope mass
     .    czin_exp,        ! Z impurity charge number
     .    iounit,          ! fortran io unit # for output
     .    weiland_output,  ! i character variable contains file name..
c                          ..for diagnostic printout if iounit > 0
     .    chie_wn,         ! q_e driven diffusion coeff. [m^2/s]
     .    chii_wn,         ! q_i driven
     .    diff_wn,         ! particle flux driven
     .    powem,           ! q_e heat flux [MW/m3]
     .    powim,           ! q_i heat flux [MW/m3]
     .    flowm)           ! particle flux [m-2 s-1]
c
c ... This subroutine calculates transport coefficients according
c ... to the model of Weiland and Nordman.
c
c ...  Input (i): local values of Te, Ti, ne, ni and the gradients at
c ...             a specific radial position.
c ... Output (o): the transport coefficients, including off-diagonal
c ...             terms, as well as the resulting power flows q_e, q_i
c ...             and particle flux.
c
c ...  Input: see annotated argument list above (all but last 6 lines)
c
c ... Output: see annotated argument list above (        last 6 lines)
c
c ... This routine makes basically use of the routine provided by
c ... Bateman (ETAWN8), for which in the present use some changes
c ... in I/O have been made.
c
c ... Three external subprograms need to be provided:
c
c     GVCRG  IMSL subroutine to compute eigenvalues and eigenvectors
c     GPIRG  IMSL function   to compute the performance index
c     STOP        subroutine to terminate, issue message, set exit value
c
c ... References:
c ...    H. Nordman, J. Weiland, and A. Jarmen, "Simulation of
c ...    toroidal drift mode turbulence driven by temperature
c ...    gradients and electron trapping," Nucl. Fusion
c ...    {\bf 30} (1990) 983--996.
c
c ...    J. Weiland and H. Nordman, "Drift wave model for inward
c ...    energy transport in tokamak plasmas," Institute for
c ...    Electromagnetic Field Theory and Plasma Physics,
c ...    Gothenburg, Sweden, (1992) CTH-IEFT/PP-1992-13 ISSN.
c
c ...    J. Weiland, A.B. Jarm\'{e}n, and H. Nordman, "Diffusive
c ...    particle and heat pinch effects in toroidal plasmas,"
c ...    Nucl. Fusion {\bf 29} (1989) 1810--1814.
c
c ...    and many more (see notes of Batemen)
c
c ... This routine has been written by Joop Konings
c                                      619/455-2261
c                                      konings@gav.gat.com
c ... Date created: 4/11/95
c
c ... ===========================================================
c ... Parameter declaration
c ... ===========================================================
c
      implicit none
c
      integer    nmaxwn
      real*8     kevdsecpmw, pi
      character  weiland_output*(*)
c
      parameter (kevdsecpmw = 1.6022e-19*1.0e3*1.0e-6) ! keV/sec per MW
      parameter (pi     = 3.1415926)
      parameter (nmaxwn = 5)
c
      integer    letain_wn6, letain_wn7, lprintin_wn, neq_wn,
     .           ndim_wn, nmodes_wn, iounit
c
      real*8     cetain_wn32, cetain_wn30, sqrtelong_wn,
     .           epsnhin_wn, epsnzin_wn, epsnsin_wn, epstein_wn,
     .           epsthin_wn, epstzin_wn, tauhin_wn, tauzin_wn,
     .           fnzin_wn, czin_wn,   azin_wn, fnsin_wn,
     .           ftrapein_wn, ekyrhoin_wn,
     .           ekparlin_wn, zpmnh_wn, zpmnz_wn, zpmns_wn, zpmte_wn,
     .           zpmne_wn, zpmth_wn, zpmtz_wn,
     .           ky_wn, rhos_wn, omega_de_wn, cs_wn,
     .           gradrhosq_wn, sfactor_wn
c
      real*8     omega_wn(1:nmaxwn), gamma_wn(1:nmaxwn),
     .           difthi_wn(1:nmaxwn,1:nmaxwn), velthi_wn(1:nmaxwn),
     .           chieff_wn(1:nmaxwn), perform_wn(1:nmaxwn)
c
      real*8     arho_exp,
     .           bt_exp, rmajor_exp, rho_exp, elonga_exp,
     .           te_exp, ti_exp, ne_exp, ni_exp, ftrapein_exp,
     .           gradte_exp, gradti_exp, gradne_exp, gradni_exp,
     .           tz_exp, nz_exp, gradtz_exp, gradnz_exp,
     .           azin_exp, czin_exp,
     .           chie_wn, chii_wn, diff_wn,
     .           powem, powim, flowm
c
c     NOTE regarding azin_exp: it should be mZ/mH, impurity mass to
c     hydrogen isotope mass (according to comments in subroutine),
c     but is calculated as mZ/m(species_1). this is being checked.
c
      real*8     zepsmach, zepsqrt
      data       zepsmach /0.0 /
c
      ndim_wn = nmaxwn                ! HSJ 2/26/96
c
c ... ===========================================================
c ... Find smallreal on this machine
c ... ===========================================================
c
      if (zepsmach .eq. 0.0) then
        zepsmach = 0.5
    2   if (0.5 * zepsmach + 1.0 .gt. 1.0) then
          zepsmach = 0.5 * zepsmach
          go to 2
        end if
        zepsqrt = SQRT (zepsmach)
      end if
c
c ... ===========================================================
c ... Prepare the correct quantities for subroutine ETAWN8_12
c ... ===========================================================
c
c ... Quantities like zpmte etc. are not essential, but used to
c ... keep further consistency with Ron Waltz' code.
c ... Assume for the moment that:
c ...      - two ion species with Z=1 and Z=czin_wn respectively
c ...      - a flat Zeff profile
c ...      - Tz = Ti
c ...      - no suprathermal hydrogenic ions
c
      gradrhosq_wn = (1.0 + elonga_exp**2) / 2.0
     .                    / elonga_exp     ! noncircular geometry factor
        sfactor_wn =  2.0 * pi * arho_exp * rho_exp
     .              * 2.0 * pi * rmajor_exp   ! flux surface area [m**2]
      sqrtelong_wn =  SQRT (elonga_exp)
c
c ... check validity of input data
c
      if (ABS (te_exp) .lt. zepsqrt)  te_exp = SIGN (zepsqrt,te_exp)
      if (ABS (ti_exp) .lt. zepsqrt)  ti_exp = SIGN (zepsqrt,ti_exp)
      if (ABS (ne_exp) .lt. zepsqrt)  ne_exp = SIGN (zepsqrt,ne_exp)
      if (ABS (ni_exp) .lt. zepsqrt)  ni_exp = SIGN (zepsqrt,ni_exp)
      if (ABS (tz_exp) .lt. zepsqrt)  tz_exp = SIGN (zepsqrt,tz_exp)
      if (ABS (nz_exp) .lt. zepsqrt)  nz_exp = SIGN (zepsqrt,nz_exp)
c
      zpmte_wn = ABS (gradte_exp / te_exp * arho_exp) ! a/LTe
      zpmth_wn = ABS (gradti_exp / ti_exp * arho_exp) ! a/LTi
      zpmne_wn = ABS (gradne_exp / ne_exp * arho_exp) ! a/Lne
      zpmnh_wn = ABS (gradni_exp / ni_exp * arho_exp) ! a/Lni
      zpmns_wn = 1.0e+6        ! suprathermal hydrogenic ions
      zpmtz_wn = ABS (gradtz_exp / tz_exp * arho_exp) ! a/LTi
      zpmnz_wn = ABS (gradnz_exp / nz_exp * arho_exp) ! a/Lni
c
c ... check validity of input data
c
      if (ABS (zpmte_wn).lt.zepsqrt)  zpmte_wn = SIGN (zepsqrt,zpmte_wn)
      if (ABS (zpmth_wn).lt.zepsqrt)  zpmth_wn = SIGN (zepsqrt,zpmth_wn)
      if (ABS (zpmne_wn).lt.zepsqrt)  zpmne_wn = SIGN (zepsqrt,zpmne_wn)
      if (ABS (zpmnh_wn).lt.zepsqrt)  zpmnh_wn = SIGN (zepsqrt,zpmnh_wn)
      if (ABS (zpmns_wn).lt.zepsqrt)  zpmns_wn = SIGN (zepsqrt,zpmns_wn)
      if (ABS (zpmtz_wn).lt.zepsqrt)  zpmtz_wn = SIGN (zepsqrt,zpmtz_wn)
      if (ABS (zpmnz_wn).lt.zepsqrt)  zpmnz_wn = SIGN (zepsqrt,zpmnz_wn)
c
      epsnhin_wn  = arho_exp/(sqrtelong_wn*rmajor_exp*zpmnh_wn) ! L_nH/R
      epsnzin_wn  = arho_exp/(sqrtelong_wn*rmajor_exp*zpmnz_wn) ! L_nZ/R
      epsnsin_wn  = arho_exp/(sqrtelong_wn*rmajor_exp*zpmns_wn) ! L_ns/R
      epstein_wn  = arho_exp/(sqrtelong_wn*rmajor_exp*zpmte_wn) ! L_Te/R
      epsthin_wn  = arho_exp/(sqrtelong_wn*rmajor_exp*zpmth_wn) ! L_Th/R
      epstzin_wn  = arho_exp/(sqrtelong_wn*rmajor_exp*zpmtz_wn) ! L_Tz/R
c
        tauhin_wn = ti_exp / te_exp      ! Th/Te
        tauzin_wn = tz_exp / te_exp      ! Tz/Te
         fnzin_wn = nz_exp / ne_exp      ! nZ/ne
         fnsin_wn = 0.0 ! ns/ne fraction of superthermal hydrogenic ions
      ftrapein_wn = ftrapein_exp         ! fraction of trapped electrons
          azin_wn = azin_exp
          czin_wn = czin_exp
c
c ... The following parameters control the calculation of
c ... subroutine ETAWN8_12:
c
       letain_wn6  =  3        ! use EISPACK as eigenvalue solver
       letain_wn7  =  0        ! set the convective velocities to zero..
c                              ..and rescale the diffusivity matrix
       cetain_wn30 =  0.01     ! finite difference used to construct..
c                              ..transport matrix
       cetain_wn32 =  0.0001   ! tolerance used in NAG eigenvalue solver
       lprintin_wn = -1        ! no diagnostic output
       neq_wn      =  4        ! number of equations
       ekyrhoin_wn =  0.316    ! normalized poloidal wave number
       ekparlin_wn =  0.1      ! kparr Ln
c
c ... ===========================================================
c ... Call the actual routine of Bateman
c ... ===========================================================
c
c ... JAK 95/11/20
c ... Note that letain(6) = 0, 2 and 3 all give the same solution for
c ... the eigenvalue problem. For some (unknown) reason option 1 gives
c ... the correct eigenvalues but the wrong eigenvectors.
c
      call etawn8_12 (
     .   letain_wn6,  !  letain(6) = 2 use IMSL to compute eigenvalues
c                                  = 1 use NAG (F02BJE)
c                                               to compute eigenvalues
c                                  = 0 use NAG (F02AKE)
c                                               to compute eigenvalues
c                                  = 3 use EISPACK (RGG)
c                                               to compute eigenvalues
     .   letain_wn7,  !  letain(7) = 0 compute convective velocities
c                                      alternatively, set the convective
c                                      velocities to zero and rescale
c                                      the diffusivity matrix
     .   cetain_wn30, !  cetain(30)  finite difference used to construct
c                                    transport matrix
c                                    (see zgm(j1,jd) matrix)
     .   cetain_wn32, !  cetain(32) tolerance used in eigenvalue solver
     .   lprintin_wn, ! i controls printout
c                          higher values => more printout
c                       lprintin_wn = -1 => no output at all
     .   neq_wn,      ! i number of equations
     .   epsnhin_wn,  ! i L_nH/R
     .   epsnzin_wn,  ! i L_nZ/R
     .   epsnsin_wn,  ! i L_ns/R
     .   epstein_wn,  ! i L_Te/R
     .   epsthin_wn,  ! i L_Th/R
     .   epstzin_wn,  ! i L_Tz/R
     .   tauhin_wn,   ! i Th/Te
     .   tauzin_wn,   ! i Tz/Te
     .   fnzin_wn,    ! i nZ/ne
     .   czin_wn,     ! i Z impurity charge number
     .   azin_wn,     ! i mZ/mH impurity mass to hydrogen isotope mass
     .   fnsin_wn,    ! i ns/ne fraction of superthermal hydrogenic ions
     .   ftrapein_wn, ! i fraction of trapped electrons
     .   ekyrhoin_wn, ! i normalized poloidal wave number
     .   ekparlin_wn, ! i kparr Ln
     .   ndim_wn,     ! i first dimension of the 2-D array difthi
c                       and maximum number of unstable modes allowed
     .   iounit,      ! fortran unit for output
     .   weiland_output, ! i character variable contains file name
c                          for diagnostic printout if iounit > 0
     .   omega_wn,    ! o real part of the frequencies
c                         normalized by $ \omega_{De} $
     .   gamma_wn,    ! o growth rates (normalized)
     .   difthi_wn,   !   diffusivity matrix
c                         normalized by $ k_y^2 / \omega_{De} $
     .   velthi_wn,   ! o convective velocities
c                         normalized by $ k_y / \omega_{De} $
     .   chieff_wn,   ! o effective total diffusivities for nhTh,
c                         nh, neTe, nZ, nZTZ
c                         normalized by $ k_y^2 / \omega_{De} $
     .   nmodes_wn,   ! o number of unstable modes
     .   perform_wn)  ! o performance index from IMSL GPIRG
c
c ... ===========================================================
c ... Output
c ... ===========================================================
c
c      Note that here rhos_wn = rhos /Amass**0.5
c                 and   cs_wn = cs_wn*Amass**0.5
c      (so Amass cancels in omega_de_wn.)
c
****  rhos_wn = ((1.02e2*(te_exp*1.0e3)**.5)/bt_exp/1.0e4) * 1.0e-2
      rhos_wn = 32.3 * SQRT (te_exp)/bt_exp  ! = rho_s/sqrt(a_i)
      ky_wn   = ekyrhoin_wn/rhos_wn
      cs_wn   = 9.79e5*(te_exp*1.0e3)**.5 * 1.0e-2
      omega_de_wn = (2.0*ky_wn*rhos_wn*cs_wn)/(rmajor_exp)
c
c ... Note that in the definitions of Bateman:
c
c ...    d( (n~ +T~)v e )/dt = chieff_wn * gradTe
c ...    where n~ is the perturbed n, etc.
c
c ...    But the heat flux is 3/2 d(nT)/dt = 3/2 d( (n~+T~)v~e )/dt so
c
c ...      chie = 3/2 * chieff_wn
c
c ...    and his actual power balance equation reads:
c
c ...      d(nT)/dt = 2/3(nabla . q_e) + 2/3 sourceterms
c
c ... Renormalize the elements in the transport matrix
c ... Note that off-diagonal coupling is implicitly done in the
c ... calculation of the heat fluxes (from which the diffusion
c ... coefficients were obtained).
c
      chie_wn = chieff_wn(3) * omega_de_wn/(ky_wn**2.0) * 3.0/2.0
      chii_wn = chieff_wn(1) * omega_de_wn/(ky_wn**2.0) * 3.0/2.0
      diff_wn = chieff_wn(2) * omega_de_wn/(ky_wn**2.0) * 1.0
c
c ... Calculate the resulting fluxes
c
      powem = kevdsecpmw*te_exp*ne_exp*1.0e19/arho_exp*gradrhosq_wn*
     .        sfactor_wn*(chie_wn*zpmte_wn)
      powim = kevdsecpmw*ti_exp*ni_exp*1.0e19/arho_exp*gradrhosq_wn*
     .        sfactor_wn*(chii_wn*zpmth_wn)
      flowm = kevdsecpmw*1.0*ne_exp*1.0e19/arho_exp*gradrhosq_wn*
     .        sfactor_wn*(diff_wn*zpmne_wn)
      return
c
      end
      subroutine weiland_setup (j, azin_exp, arho_exp, rmajor_exp)
c

c
c ----------------------------------------------------------------------
c    ONETWO interface for Weiland routine                    HSJ 2/16/96
c    map onetwo input variables into the following input
c                               for WEILAND_12
c
c    arho_exp,        ! minor radius [m]
c    rmajor_exp,      ! major radius [m]
c    rho_exp,         ! local r/a
c    elonga_exp,      ! local plasma elongation
c    bt_exp,          ! toriodal magnetic field [Gauss]
c        te_exp,      !     Te [keV]
c    gradte_exp,      ! gradTe [keV/m]
c        ti_exp,      !     Ti [keV]
c    gradti_exp,      ! gradTi [keV/m]
c        tz_exp,      !     Tz [keV]   main impurity species
c    gradtz_exp,      ! gradTz [keV/m] main impurity species
c        ne_exp,      !     ne [*10^19 m-3]
c    gradne_exp,      ! gradne [*10^19 m-4]
c        ni_exp,      !     ni [*10^19 m-3]
c    gradni_exp,      ! gradni [*10^19 m-4]
c        nz_exp,      !     nz [*10^19 m-3] impurity
c    gradnz_exp,      ! gradnz [*10^19 m-4] impurity
c    ftrapein_exp,    ! fraction trapped electrons
c    azin_exp,        ! mZ/mH impurity mass to H isotope mass
c    czin_exp,        ! Z impurity charge number
c
c    map the following weiland_12 output to onetwo variables:
c
c    chie_wn,         ! q_e driven diffusion coeff. [m^2/s]
c    chii_wn,         ! q_i driven
c    diff_wn,         ! particle flux driven
c    powem,           ! q_e heat flux [MW/m3]
c    powim,           ! q_i heat flux [MW/m3]
c    flowm            ! particle flux [m-2 s-1]
c
c ----------------------------------------------------------------------
c
      USE param
      USE ions
      USE soln
      USE numbrs
      USE mesh
      USE machin
      USE geom
      USE flx
      USE neo2d
      USE weiland
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      real*8       ne_exp, ni_exp, nz_exp
      character*32 weiland_output
      data         weiland_output /'weiland_output'/
c
c      include 'param.i'
c      include 'flx.i'     ! flux
c      include 'geom.i'    ! fcap,rcap,etc
c      include 'ions.i'    ! atw,etc
c      include 'numbrs.i'  ! nprim,nj,etc
c      include 'machin.i'  ! btor,rmajor,etc.
c      include 'mesh.i'    ! r,dr
c      include 'neo2d.i'   ! elong_r,ftncl
c      include 'soln.i'    ! te,ti,ene,en
c      include 'weiland.i' ! chie,chii,etc
c
      iounit_weiland = 0  ! set to a Fortran unit not used in ONETWO
c                           to get debug output in "weiland_output" file
c                           note that this file is recreated each time
c                           that weiland_12 is called.
c                           If non-zero, returned by getioun before file open.
c
c     on the half grid (i.e., at (r(j+1)+r(j))/2.0
c
      k = 1               ! use first primary ion ??
c
      rho_exp    = 0.5*(r(j)+r(j+1))/100.
      elonga_exp = 0.5*(elong_r(j)+elong_r(j+1))
      fcapa      = 0.5*(fcap(j)+fcap(j+1))
      rcapa      = 0.5*(rcap(j)+rcap(j+1))
      bt_exp     = ABS (btor*rmajor / (fcapa*rcapa)) ! flux surf avg B_T
      te_exp     = 0.5*(te(j)+te(j+1))
      ti_exp     = 0.5*(ti(j)+ti(j+1))
      tz_exp     = ti_exp
      gradte_exp = 100.0*(te(j+1)-te(j))/dr(j)
      gradti_exp = 100.0*(ti(j+1)-ti(j))/dr(j)
      gradtz_exp = 100.0*(ti(j+1)-ti(j))/dr(j)
      ne_exp     = 0.5e-13*(ene(j)+ene(j+1))
      ni_exp     = 0.5e-13*(en(j,k)+en(j+1,k))
      nz_exp     = 0.5e-13*(en(j,nprim+1)+en(j+1,nprim+1))
      gradne_exp = 1.0e-11*(ene(j+1)-ene(j))/dr(j)
      gradni_exp = 1.0e-11*(en(j+1,1)-en(j,1))/dr(j)
      gradnz_exp = 1.0e-11*(en(j+1,nprim+1)-en(j,nprim+1))/dr(j)
      czin_exp   = 0.5*(z(j,nprim+1)+z(j+1,nprim+1))
      ftrapein_exp =  0.5*(ftncl(j)+ftncl(j+1))
c
c --- the folowing are test values used for debug
c
****  arho_exp    = 0.88
****  rmajor_exp  = 1.7
****  rho_exp     = 0.5
****  elonga_exp  = 1.6
****  bt_exp      = 1.0
****  te_exp      = 0.7
****  gradte_exp  = 0.5
****  ti_exp      = te_exp
****  gradti_exp  = gradte_exp
****  tz_exp      = ti_exp
****  gradtz_exp  = gradti_exp
****  ne_exp      = 2.0
****  gradne_exp  = 6.0
****  ni_exp      = 1.5
****  gradni_exp  = 3.0
****  azin_exp    = 12
****  czin_exp    = 6
****  nz_exp      = (ne_exp-ni_exp)/czin_exp
****  gradnz_exp  = gradni_exp
****  ftrapein_exp = 0.5
c
      call weiland_12 (
     .                arho_exp,        ! minor radius [m]
     .                rmajor_exp,      ! major radius [m]
     .                rho_exp,         ! local r/a
     .                elonga_exp,      ! local plasma elongation
     .                bt_exp,          ! toriodal magnetic field [Gauss]
     .                    te_exp,      !     Te [keV]
     .                gradte_exp,      ! gradTe [keV/m]
     .                    ti_exp,      !     Ti [keV]
     .                gradti_exp,      ! gradTi [keV/m]
     .                    tz_exp, !      Tz [keV]  main impurity species
     .                gradtz_exp, ! gradTz [keV/m] main impurity species
     .                    ne_exp,      !     ne [*10^19 m-3]
     .                gradne_exp,      ! gradne [*10^19 m-4]
     .                    ni_exp,      !     ni [*10^19 m-3]
     .                gradni_exp,      ! gradni [*10^19 m-4]
     .                    nz_exp,      !     nz [*10^19 m-3] impurity
     .                gradnz_exp,      ! gradnz [*10^19 m-4] impurity
     .                ftrapein_exp,    ! fraction trapped electrons
     .                azin_exp,    ! mZ/mH impurity mass to isotope mass
     .                czin_exp,        ! Z impurity charge number
     .          iounit_weiland,        ! fortran unit # for i/o
     .          weiland_output,        ! i character variable contains
     .                                 ! file name for diagnostic
     .                                 ! printout if iounit .gt. 0
     .                chie_wn,     ! q_e driven diffusion coeff. [m^2/s]
     .                chii_wn,         ! q_i driven
     .                diff_wn,         ! particle flux driven
     .                powem,           ! q_e heat flux [MW/m3]
     .                powim,           ! q_i heat flux [MW/m3]
     .                flowm)           ! particle flux [m-2 s-1]
c
      xchie_weiland(j) = 1.0e4*chie_wn
      xchii_weiland(j) = 1.0e4*chii_wn
      d_weiland(j)     = 1.0e4*diff_wn
      qe_weiland(j)    = powem
      qi_weiland(j)    = powim
****  flux(k,j)        = 1.0e-4*flowm               ! [cm-2 s-1 ]
      return
c
      end

      subroutine write_sol (ntrplt)
c
      USE param
      USE monitr
      implicit none
c
c ----------------------------------------------------------------------
c     output info for monitoring the solution
c     into file associated with ntrplt
c     INPUT
c       kk,kj     parameters from param.i
c       all of monitr.i
c     OUTPUT
c       normalized values to file unit ntrplt
c ----------------------------------------------- 2-20-97 ----- HSJ ----
c
      integer ntrplt, j, k
      real*8  totmnn, totl
c
c      include 'param.i'
c      include 'monitr.i'
c
      totmnn = FLOAT (imonn )
      totl   = FLOAT (imontf)
      if (totl .gt. 0.0) then
        do   k=1,kk
          do j=1,kj
            rmonu(j,k) = rmonu(j,k) / totl
          end do
        end do
      end if
      write  (ntrplt, '(6(1pe14.6))') ((rmonu(j,k), j=1,kj), k=1,kk)
      write  (ntrplt, '(6(1pe14.6))')   totl, totmnn  !these formats must be compatible with ntrplt
      return                                          !formats appearing in preplt
c
      end
