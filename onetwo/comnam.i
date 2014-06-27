! -*-f90-*-
!23456789-123456789-123456789-123456789-123456789-123456789-123456789-23
! -*-f90-*-
!
!********************
!       input namelist 'edata' is contained on ecrhdat file
!********************
!
      character*8 raypatt
!
      integer mray(16), inang, gauszone
!
      REAL(p_) thetar(256), phir(256), cr(16), tempie, theta,               &
     &   tangl, bsratio, bhalf, ar, br
!
      common/edata_char/ raypatt
      common/edata_int/ mray, inang, gauszone
      common/edat_real/ br, ar, bhalf, bsratio, tangl, theta, tempie, cr    &
     &   , phir, thetar
      namelist /edata/ thetar, phir, fmu, sgn, cfe0, pfe0, cr, smax,        &
     &   phimin, phimax, thtmin, thtmax, ds, fdout, tout, relerr, abserr    &
     &   , angrid, ezeff, tolcur, tempie, thtc, powinc, bwidth, pdfrac,     &
     &   rea10, rea20, theta, tangl, frctal, bsratio, bhalf, ar, br,        &
     &   alphp, alphc, alphb, rhop, gamalp, tperpp, tperpc, tperpb,         &
     &   gamperp, tparfc, tparfp, tparfb, gamparf, tparbc, tparbp,          &
     &   tparbb, gamparb, cqldat, pwrfmin, scrapef, mray, idamp, igafit,    &
     &   io1, io2, io3, io4, io5, io6, io7, io8, io9, io0, inang, ima10,    &
     &   ima20, iparam, mmzons, mn0, nlmesh, modelc, nbfld, nlcart,         &
     &   nlout, nray, nlpol, ntheta, numphi, netcdfdat, nharm, gauszone,    &
     &   raypatt, dsmin
!
! ----------------------------------------------------------------------
!   YRLL(02/08/00)
!       new variable "scrapef"" is added to this namelist
!       see comments in comeqb.i
!   YRLL(11/10/92)
!       new variable "pwrfmin" is added to this namelist
!       see comments in comend.i
!   YRLL(10/06/93)
!       tprof is a dead variable.  It is removed from namelist.
!   JEK(08/10/01) added raypatt, gauszone
! to namelist for ray pattern generator
! ----------------------------------------------------------------------
!
