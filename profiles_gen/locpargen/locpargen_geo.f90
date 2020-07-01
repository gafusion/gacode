subroutine locpargen_geo

  use locpargen_globals
  use expro_locsim_interface

  implicit none

  integer :: j1,j2
 
  open(unit=1,file='input.geo',status='replace')
  write(1,'(a)') '# input.geo'
  write(1,'(a)') '#'
  write(1,'(a,f10.6)') '# NOTE: Derived from input.gacode.geo at r/a=',rmin_loc
  write(1,'(a)') '# Lengths normalized to a' 
  write(1,'(a,f10.6)') '# ASPECT_RATIO=',rmaj_loc
  write(1,'(a,f10.6)') '# SAFETY_FACTOR=',q_loc
  write(1,'(a,f10.6)') '# SHEAR=',s_loc
  write(1,'(a,f4.1)') '# BTCCW=',btccw
  write(1,'(a,f4.1)') '# IPCCW=',ipccw
  write(1,'(a)') '#'
  write(1,'(a)') '# File format:'
  write(1,'(a)') '#-------------------'
  write(1,'(a)') '# nfourier'
  write(1,'(a)') '# a[8,0:nfourier]'    
  write(1,'(a)') '#-------------------'
  write(1,'(a)') '#'
  write(1,*) geo_ny_loc

  do j2=0,geo_ny_loc
     do j1=1,8
        write(1,'(1pe20.13)') geo_yin_loc(j1,j2)
     enddo
  enddo

  print '(a)','INFO: (locpargen_geo) Wrote input.geo.'

  close(1)
  
end subroutine locpargen_geo
