subroutine locpargen_geo

  use EXPRO_interface
  use EXPRO_locsim_interface

  implicit none

  integer :: j1,j2
  real, dimension(1) :: x,y
  real, dimension(:,:,:), allocatable :: geo_p
  real, dimension(:), allocatable :: x_vec

  allocate(x_vec(EXPRO_n_exp))
  allocate(geo_p(8,0:EXPRO_nfourier,EXPRO_n_exp))

  geo_p(1:4,:,:) = EXPRO_geo(:,:,:)/EXPRO_rmin(EXPRO_n_exp)
  geo_p(5:8,:,:) = EXPRO_dgeo(:,:,:)

  open(unit=1,file='input.geo',status='replace')
  write(1,'(a)') '# input.geo'
  write(1,'(a)') '#'
  write(1,'(a)') '# See https://fusion.gat.com/theory/input.geo for complete documentation.'
  write(1,'(a)') '#'
  write(1,'(a,f10.6)') '# NOTE: Derived from input.profiles.geo at r/a=',x
  write(1,'(a)') '# Lengths normalized to a' 
  write(1,'(a,f10.6)') '# ASPECT_RATIO=',rmaj_loc
  write(1,'(a,f10.6)') '# SAFETY_FACTOR=',q_loc
  write(1,'(a,f10.6)') '# SHEAR=',s_loc
  write(1,'(a,i3)') '# BTCCW=',-EXPRO_signb
  write(1,'(a,i3)') '# IPCCW=',-EXPRO_signb*EXPRO_signq
  write(1,'(a)') '#'
  write(1,'(a)') '# File format:'
  write(1,'(a)') '#-------------------'
  write(1,'(a)') '# nfourier'
  write(1,'(a)') '# a[8,0:nfourier]'    
  write(1,'(a)') '#-------------------'
  write(1,'(a)') '#'
  write(1,*) EXPRO_nfourier

  do j2=0,EXPRO_nfourier
     do j1=1,8
        call cub_spline(x_vec,geo_p(j1,j2,:),EXPRO_n_exp,x,y,1)
        write(1,'(1pe20.13)') y(1)
     enddo
  enddo

  print *
  print '(a)','INFO: (locpargen) Wrote input.geo.'

  close(1)

  deallocate(x_vec)
  deallocate(geo_p)
  
end subroutine locpargen_geo
