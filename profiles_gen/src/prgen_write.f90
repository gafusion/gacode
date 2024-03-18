!---------------------------------------------------------
! prgen_write.f90
!
! PURPOSE:
!  Write input.gacode 
!----------------------------------------------------------

subroutine prgen_write

  use prgen_globals
  use expro

  implicit none
  integer :: i

  print '(a)','INFO: (prgen_write) Created these species'

  do i=1,expro_n_ion
     print '(t6,i2,2(1x,a))',i,trim(expro_name(i)),trim(expro_type(i))
  enddo

  !-------------------------------------------------------------------------------------
  ! Map (common) data from EFIT analysis [skip if reading input.* without gfile]
  !
  if (format_type < 7) then 
     expro_torfluxa = -btccw*abs(torfluxa)
     expro_rcentr   = rcentr
     expro_bcentr   = -btccw*abs(bcentr)
     expro_current  = -ipccw*abs(current)

     expro_q       = ipccw*btccw*abs(q)
     expro_polflux = -ipccw*abs(dpsi)
     expro_kappa   = kappa
     expro_delta   = delta
     expro_zeta    = zeta
     expro_zmag    = zmag

     expro_shape_sin3 = shape_sin(3,:)
     expro_shape_sin4 = shape_sin(4,:)
     expro_shape_sin5 = shape_sin(5,:)
     expro_shape_sin6 = shape_sin(6,:)
     expro_shape_cos0 = shape_cos(0,:)
     expro_shape_cos1 = shape_cos(1,:)
     expro_shape_cos2 = shape_cos(2,:)
     expro_shape_cos3 = shape_cos(3,:)
     expro_shape_cos4 = shape_cos(4,:)
     expro_shape_cos5 = shape_cos(5,:)
     expro_shape_cos6 = shape_cos(6,:)

     expro_z_eff = zeff

     ! EFIT passthrough functions
     expro_ptot = p_tot
     expro_fpol = fpol
  endif
  !-------------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------------
  ! Now we need to write the profile data
  !
  expro_head_original  = '#  *original : '//trim(date)
  expro_head_statefile = '# *statefile : '//trim(file_state)
  expro_head_gfile     = '#     *gfile : '//trim(file_g)
  expro_head_cerfile   = '#   *cerfile : '//trim(file_cer)

  ! Run check for Jacobian
  call expro_compute_derived

  ! Use unit=2 since expro use unit=1
  open(unit=2,file='success',status='replace')
  if (expro_jerr == 0) then
     call expro_write('input.gacode')
     print '(a)','INFO: (prgen_write) Wrote input.gacode.'
     ! Successful completion
     write(2,*) 1
  else
     write(2,*) 2
  endif
  close(2)

end subroutine prgen_write

