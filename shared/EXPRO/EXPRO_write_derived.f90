!--------------------------------------------------------
! EXPRO_write_derived.f90
!
! PURPOSE:
!  Write all derived quantities to input.profiles.extra
!--------------------------------------------------------

subroutine EXPRO_write_derived(io,datafile)

  use EXPRO_globals
  use EXPRO_interface

  integer, intent(in) :: io
  character (len=*), intent(in) :: datafile

  open(unit=io,file=trim(path)//trim(datafile),status='replace')

  write(io,'(a)') '# input.profiles.extra'
  write(io,'(a)') '#'
  write(io,'(a)') '# 01    EXPRO_bunit(:)        B_unit (T)'
  write(io,'(a)') '# 02    EXPRO_s(:)            (r/q)(dq/dr)'
  write(io,'(a)') '# 03    EXPRO_drmaj(:)        dR_0/dr'
  write(io,'(a)') '# 04    EXPRO_dzmag(:)        dZ_0/dr'
  write(io,'(a)') '# 05    EXPRO_sdelta(:)       r d(delta)/dr'
  write(io,'(a)') '# 06    EXPRO_skappa(:)       (r/kappa) d(kappa)/dr'
  write(io,'(a)') '# 07    EXPRO_szeta(:)        r d(zeta)/dr'
  write(io,'(a)') '# 08    EXPRO_dlnnedr(:)      -dln(ne)/dr (1/m)'
  write(io,'(a)') '# 09    EXPRO_dlntedr(:)      -dln(Te)/dr (1/m)'
  write(io,'(a)') '# 10-19 EXPRO_dlnnidr(1:10,:) -dln(ni)/dr (1/m)'
  write(io,'(a)') '# 20-29 EXPRO_dlntidr(1:10,:) -dln(ti)/dr (1/m)'
  write(io,'(a)') '# 30    EXPRO_dlnptotdr(:)    -dln(ptot)/dr (1/m)'
  write(io,'(a)') '# 31    EXPRO_drdrho(:)        dr/d(rho) [DIMENSIONLESS]'
  write(io,'(a)') '# 32    EXPRO_w0p(:)           d(w0)/dr (1/s/m)'
  write(io,'(a)') '# 33    EXPRO_vol(:)           V (m^3)'
  write(io,'(a)') '# 34    EXPRO_volp(:)          dV/dr (m^2)'
  write(io,'(a)') '# 35    EXPRO_cs(:)            cs (m/s)'
  write(io,'(a)') '# 36    EXPRO_rhos(:)          rhos (m)'
  write(io,'(a)') '# 37    EXPRO_ni_new(:)        ni [Corrected for quasin.]'
  write(io,'(a)') '# 38    EXPRO_dlnnidr_new(:)   -dln(ni)/dr (1/m) [Corrected for quasin.]' 
  write(io,'(a)') '# 39    EXPRO_grad_r0(:)       |grad r| at theta=0 '
  write(io,'(a)') '# 40    EXPRO_ave_grad_r(:)    Flux-surface average <|grad r|>' 
  write(io,'(a)') '# 41    EXPRO_bp0(:)           B_pol at theta=0 (T)'
  write(io,'(a)') '# 42    EXPRO_bt0(:)           B_tor at theta=0 (T)'
  write(io,'(a)') '# 43    EXPRO_gamma_e(:)       r/q d(w0)/dr (1/s)'
  write(io,'(a)') '# 44    EXPRO_gamma_p(:)       R_0 d(w0)/dr (1/s)'
  write(io,'(a)') '# 45    EXPRO_mach(:)          R_0 w0/cs'
  write(io,'(a)') '# 46    EXPRO_ip(:)            Plasma current (A) = mu Int[Bp dl]'
  write(io,'(a)') "# 47    EXPRO_adlnnedr(:)      -ne''/ne (1/m^2)"
  write(io,'(a)') "# 48    EXPRO_adlntedr(:)      -Te''/Te (1/m^2)"
  write(io,'(a)') "# 49-58 EXPRO_sdlnnidr(1:10,:) -ni''/ni (1/m^2)"
  write(io,'(a)') "# 59-68 EXPRO_sdlntidr(1:10,:) -Ti''/Ti (1/m^2)"
  write(io,'(a)') '# 69    EXPRO_nuee(:)          nuee (1/s)'
  write(io,'(a)') '#'
  write(io,'(a,i3)') '# Each vector has length ',EXPRO_n_exp
  write(io,'(a)') '#'

  write(io,10) EXPRO_bunit(:)    ! 1
  write(io,10) EXPRO_s(:)        ! 2
  write(io,10) EXPRO_drmaj(:)    ! 3
  write(io,10) EXPRO_dzmag(:)    ! 4 
  write(io,10) EXPRO_sdelta(:)   ! 5  
  write(io,10) EXPRO_skappa(:)   ! 6
  write(io,10) EXPRO_szeta(:)    ! 7
  write(io,10) EXPRO_dlnnedr(:)  ! 8
  write(io,10) EXPRO_dlntedr(:)  ! 9 
  write(io,10) transpose(EXPRO_dlnnidr(:,:)) ! 10-19
  write(io,10) transpose(EXPRO_dlntidr(:,:)) ! 20-29
  write(io,10) EXPRO_dlnptotdr(:)! 30 
  write(io,10) EXPRO_drdrho(:)   ! 31
  write(io,10) EXPRO_w0p(:)      ! 32
  write(io,10) EXPRO_vol(:)      ! 33
  write(io,10) EXPRO_volp(:)     ! 34
  write(io,10) EXPRO_cs(:)       ! 35
  write(io,10) EXPRO_rhos(:)     ! 36
  write(io,10) EXPRO_ni_new(:)   ! 37
  write(io,10) EXPRO_dlnnidr_new(:) ! 38 
  write(io,10) EXPRO_grad_r0(:)     ! 39
  write(io,10) EXPRO_ave_grad_r(:)  ! 40 
  write(io,10) EXPRO_bp0(:)      ! 41
  write(io,10) EXPRO_bt0(:)      ! 42
  write(io,10) EXPRO_gamma_e(:)  ! 43
  write(io,10) EXPRO_gamma_p(:)  ! 44
  write(io,10) EXPRO_mach(:)     ! 45
  write(io,10) EXPRO_ip(:)       ! 46
  write(io,10) EXPRO_sdlnnedr(:) ! 47
  write(io,10) EXPRO_sdlntedr(:) ! 48
  write(io,10) transpose(EXPRO_sdlnnidr(:,:)) ! 49-58
  write(io,10) transpose(EXPRO_sdlntidr(:,:)) ! 59-68
  write(io,10) EXPRO_nuee(:)     ! 69

  close(io)

10 format(t2,1pe14.7)

end subroutine EXPRO_write_derived
