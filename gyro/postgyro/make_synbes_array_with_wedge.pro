FUNCTION make_synbes_array_with_wedge, data, wedge, $
  ITMIN = itmin, ITMAX = itmax, OMEGA0=omega0, $
  NR_BESchannels = NR_BESchannels, $
  NZ_BESchannels = NZ_BESchannels, BES_Z00 = BES_Z00, BES_DR=BES_DR, $
  BES_DZ=BES_DZ, SHOW_WINDOWS=show_windows, $
  Aphys = Aphys, A_over_Cs = a_over_cs
;
; C. Holland, UCSD
; v1.0: 3.26.2012: takes postgyro data file plus HDF5 wedge data file
; and generates synthetic BES array using wedge output
;
; use
; IDL> data = GET_GYRO_DATA(simdir, /HDF5)
; IDL> wedge = GET_GYRO_HDF5_WEDGE(data)
; IDL> results = MAKE_SYNBES_ARRAY(data, wedge, ...)
;
; 
; SF is intepolation in theta (set equal to
; data.profile_data.n_theta_mult to use GYRO-caluclated nu_geo)
; ITMIN/ITMAX specify starting/stopping time indexes
; INTERP does time interpolation needed b/c calc done in lab frame
; generally spacing = 0.25 should be ok.
; OMEGA0 is equilibrium ExB rotation frequency at mid-radius norm point)
;  n_tor_Frac = # of equally spaced toroidal locations for 2D
;  calculations
;  BES parameters specifcy # of (R,Z) points, and spacing between them
;

  IF KEYWORD_SET(show_windows) THEN windows_flag=1 ELSE windows_flag=0
  PRINT, 'windows_flag = ', windows_flag

  DEFAULT, omega0, data.w0[data.n_r/2]
  PRINT, 'omega0 (c_s/a): ', omega0
  DEFAULT, n_tor_frac, wedge.n_tor_frac
  PRINT, '# toroidal locations: ', n_tor_frac
  delta_tor_frac = 1./n_tor_frac;
  DEFAULT, itmin, 0
  PRINT, 'itmin = ', itmin
  DEFAULT, itmax, wedge.n_time-1
  PRINT, 'itmax = ', itmax
  DEFAULT, NR_BESchannels, 5
  PRINT, 'NR_BESchannels = ', NR_BESchannels
  DEFAULT, NZ_BESchannels, 6
  PRINT, 'NZ_BESchannels = ', NZ_BESchannels
  DEFAULT, BES_Z00, -4. ;cm
  PRINT, 'BES Z00 (cm): ', BES_Z00
  DEFAULT, BES_DR, 0.9 ;cm
  PRINT, 'BES DR (cm): ', BES_DR
  DEFAULT, BES_DZ, 1.2 ;cm
  PRINT, 'BES DZ (cm): ', BES_DZ

  DEFAULT, Aphys, data.Aphys
  PRINT, 'Aphys (cm) = ', Aphys
  DEFAULT, a_over_cs, 1.e3/data.csda  ;data.csda in kHz
  PRINT, 'a_over_cs (microsec) = ', a_over_cs

  ;get GYRO (R,Z) coords from wedge file
  GYRO_R = wedge.R_GYRO
  GYRO_Z = wedge.Z_GYRO

  ;start with BES
  RESTORE, GETENV('GYRO_DIR') + '/sim/' + data.simdir + '/psf.sav'
  psf.r = (psf.r - psf.stats.r)/Aphys
  psf.p = (psf.p - psf.stats.z)/Aphys

  psf_NR= N_ELEMENTS(psf.r)
  psf_NZ = N_ELEMENTS(psf.p)
  psf_dx = (psf.r[psf_NR-1] - psf.r[0])/(psf_NR-1)
  psf_dz = (psf.p[psf_NZ-1] - psf.p[0])/(psf_NZ-1)

  n_r = wedge.n_r
  n_y = wedge.n_y
  psf_interp = FLTARR(n_r, n_y, NR_BESchannels, NZ_BESchannels)
  psf_interp_norm = FLTARR(NR_BESchannels, NZ_BESchannels)
  BES_xloc = FLTARR(NR_BESchannels, NZ_BESchannels)
  BES_zloc = FLTARR(NR_BESchannels, NZ_BESchannels)
  BES_gyroidx = LONARR(NR_BESchannels, NZ_BESchannels)

  ;R00 = major radius/a of bes channel 9 is at middle of shot
  ir = data.n_r/2
  R00 =  (data.R0[ir] + data.r[ir])
  Z00 = BES_Z00/Aphys

  ;create array BES psfs interpolated onto GYRO grid
  FOR ir = 0, NR_BESchannels-1 DO FOR iz = 0, NZ_BESchannels-1 DO BEGIN
      BES_xloc[ir,iz] = R00 + (ir-NR_BESchannels/2)*BES_DR/Aphys
      BES_zloc[ir,iz] = Z00 - iz*BES_DZ/Aphys

      ;for now, use minidx to find closest GYRO point; maybe interpolate later
      mind = MIN((GYRO_R - BES_xloc[ir,iz])^2 + (GYRO_Z - BES_zloc[ir,iz])^2, $
                 minidx)
      BES_gyroidx[ir,iz] = minidx

      ;regrid PSF onto GYRO coords
      ir_GYRO = (GYRO_R - psf.R[0] - BES_xloc[ir,iz])/psf_dx
      iz_GYRO = (GYRO_Z - psf.P[0] - BES_zloc[ir,iz])/psf_dz
      psf_interp[*,*,ir,iz] = BILINEAR(psf.psf, ir_GYRO, iz_GYRO, MISSING=0)
      psf_interp_norm[ir,iz] = TOTAL(psf_interp[*,*,ir,iz])
  ENDFOR

  it = itmin
  ne_RZ = REFORM(wedge.n_rz[data.n_kinetic-1,0,*,*,0])

  IF (windows_flag) THEN BEGIN ; plot BES locations
      WINDOW, 0
      xr =BES_XLOC[NR_BESchannels/2,0]*Aphys + [-6,6]
      GYRO_RZ_COLOR_CONTOUR, TRANSPOSE(ne_RZ), GYRO_R, GYRO_Z, $
        XRANGE=xr, YRANGE=[-10,2]+BES_Z00, Aphys=Aphys, $
        TITLE = '!4d!Xn!De!N (t = ' + NUMTOSTRING(data.t[it]) + ')'
      
      FOR ir = 0, NR_BESchannels-1 DO FOR iz=0, NZ_BESchannels-1 DO BEGIN
          PLOTS, BES_xloc[ir,iz]*Aphys, BES_zloc[ir,iz]*Aphys, psym=4
          PLOTS, GYRO_R(BES_gyroidx[ir,iz])*Aphys, $
                 GYRO_Z(BES_gyroidx[ir,iz])*Aphys, $
                 psym=2, color=50
      ENDFOR
      CONTOUR, psf.psf, (psf.r+BES_xloc[2,0])*Aphys, $
               (psf.p+BES_zloc[2,0])*Aphys, $
               LEVELS=[0.1,0.5,0.9], $
               C_LINESTYLE=[4,2,0], /OVERPLOT
      CONTOUR, psf_interp[*,*,2,0], GYRO_R*Aphys, GYRO_Z*Aphys, $
               LEVELS=[0.1,0.5,0.9], COLOR=50, $
               C_LINESTYLE=[5,3,1], /OVERPLOT
  ENDIF

  ;set up time axis and result storage arrays
  NT = itmax - itmin + 1
  t = wedge.time_skip*(itmin + FINDGEN(NT))

  GYRO_ne_sig = FLTARR(NR_BESchannels, NZ_BESchannels, N_tor_frac, NT)
  syn_ne_sig = FLTARR(NR_BESchannels, NZ_BESchannels, N_tor_frac, NT)
  C_I = COMPLEX(0,1)

  ;begin iterations
  IY = FLTARR(n_y) + 1
  FOR it = 0,NT-1 DO BEGIN
      PRINT, it,  '/', NT-1

      FOR i_tf = 0, n_tor_frac-1 DO BEGIN
          ne_RZ = REFORM(wedge.n_RZ[data.n_kinetic-1, i_tf,*,*,it])/$
                  (REFORM(data.n_eq[data.n_spec-1,*])#IY)

          FOR ir=0,NR_BESchannels-1 DO FOR iz=0,NZ_BESchannels-1 DO BEGIN
              GYRO_ne_sig[ir,iz,i_tf,it] = ne_RZ[BES_gyroidx[ir,iz]]
              syn_ne_sig[ir,iz,i_tf,it] = TOTAL(ne_RZ*psf_interp[*,*,ir,iz])/$
                psf_interp_norm[ir,iz]
          ENDFOR
      ENDFOR
  ENDFOR

  IF (windows_flag) THEN BEGIN
  	WINDOW, 1
  	PLOT, t, GYRO_ne_sig[0,0,0,*], /XS, title='GYRO ne, syn. BES (red)'
  	OPLOT, t, syn_ne_sig[0,0,0,*], COLOR=100
  ENDIF

  results = {NT:nt, T:t, $
                 ITMIN:itmin, ITMAX:itmax, $
                 Aphys:Aphys, omega0:omega0, $
                 N_tor_frac: n_tor_frac, $
                 a_over_cs: a_over_cs, $
                 NR_BESCHANNELS:NR_BESchannels, $
                 NZ_BESCHANNELS:NZ_BESchannels, $
                 GYRO_ne:GYRO_ne_sig, syn_ne:syn_ne_sig, $
                 syn_BES_XLOC:BES_xloc, syn_BES_ZLOC: BES_zloc, $
                 GYRO_BES_XLOC:GYRO_R[BES_gyroidx], $
                 GYRO_BES_ZLOC:GYRO_Z[BES_gyroidx]}

  RETURN, results
END ;make_synarrays
