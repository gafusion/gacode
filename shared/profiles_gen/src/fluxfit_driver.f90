! model (1=parameterized, 2=Fourier)
! ns (=number of Fouier harmonics)
! npsi (number of flux gridpoints)
! nd (number of arclength datapoints)
! i_print (0=quiet,1=print)

subroutine fluxfit_driver(model_in,ns_in,npsi_in,nd_in,rd_in,zd_in,i_print)

  use prgen_fluxfit_globals

  implicit none

  ! Input parameters
  integer, intent(in) :: model_in 
  integer, intent(in) :: ns_in 
  integer, intent(in) :: npsi_in 
  integer, intent(in) :: nd_in 
  real, dimension(nd_in,npsi_in) :: rd_in
  real, dimension(nd_in,npsi_in) :: zd_in
  integer, intent(in) :: i_print

  ! Internal variables
  integer :: i,j,p
  integer :: npsi
  real :: t,r,z

  integer :: im,ip,imax
  real :: z0,z2
  real :: r0,r2
  real :: rm,rp
  real :: zm,zp
  real :: dl0,l_tot,dl,rs
  real :: err

  real, dimension(3) :: s

  model = model_in
  ns = ns_in
  npsi = npsi_in
  nd = nd_in

  !------------------------------------------------------
  ! Define fit model (1=parameterized, 2=Fourier) 
  ! and parameter tags
  !
  select case (model)

  case (1)

     nc = 6

     allocate(tag(6))
     tag(1) = 'r'
     tag(2) = 'z_mag'
     tag(3) = 'rmaj'
     tag(4) = 'kappa'
     tag(5) = 'delta'
     tag(6) = 'sqr'

     allocate(c(nc)) 

  case (2)

     allocate(tag(5))
     tag(1) = 'n'
     tag(2) = 'a_R'
     tag(3) = 'b_R'
     tag(4) = 'a_Z'
     tag(5) = 'b_Z'

     allocate(ar(0:ns))
     allocate(br(0:ns))
     allocate(az(0:ns))
     allocate(bz(0:ns))

  case default

     print '(a)','ERROR: model invalid in fluxfit.'
     stop

  end select
  !------------------------------------------------------

  nd   = nd_in
  npsi = npsi_in

  allocate(rd(nd))
  allocate(zd(nd))
  allocate(theta(nd))

  p = 0 

  if (model == 1) then
     open(unit=2,file='fluxfit.profile',status='replace')
     write(2,20) tag(:) 
  else
     open(unit=2,file='fluxfit.geo',status='replace')
     write(2,*) ns
  endif

  open(unit=3,file='fluxfit.rz',status='replace')
  write(3,*) nd
  write(3,*) npsi

  open(unit=4,file='fluxfit.error',status='replace')

  do j=1,npsi

     rd(:) = rd_in(:,j)
     zd(:) = zd_in(:,j)

     ! Find the centroid (z_c,r_c)

     call fluxfit_moments(rd,zd,nd,s)
     r_c = s(2)/s(1)
     z_c = s(3)/s(1)

     ! Find major radii (rp,rm) at the height of 
     ! the centroid.

     call fluxfit_intersect(rd,zd,nd,z_c,rp,rm)

     rmin = 0.5*(rp-rm)
     rmaj = 0.5*(rp+rm)

     if (model == 2) then

        !-------------------
        ! Fourier (fr_model)
        !-------------------

        ! Find rightmost R (r0)

        call fluxfit_minmax(rd,zd,nd,r0,z0,'max')

        imax = maxloc(rd,1) 
        r2 = rd(imax)
        z2 = zd(imax)

        if (z2 > z0) then
           dl0 = sqrt((r2-r0)**2+(z2-z0)**2)
        else
           dl0 = -sqrt((r2-r0)**2+(z2-z0)**2)
        endif

        l_tot = 0.0
        do i=1,nd-1
           dl = sqrt((rd(i+1)-rd(i))**2+(zd(i+1)-zd(i))**2)
           l_tot = l_tot+dl
        enddo

        ! Construct poloidal angle starting at r0 and 
        ! winding counter-clockwise.

        theta(imax) = 2*pi*(dl0)/l_tot
        do i=imax,imax+nd-2
           ip = modulo(i,nd)+1
           im = modulo(i-1,nd)+1
           dl = sqrt((rd(ip)-rd(im))**2+(zd(ip)-zd(im))**2)
           theta(ip) = theta(im) + 2*pi*dl/l_tot
        enddo

        call fluxfit_fourier()

        do i=0,ns
           write(2,10) ar(i),br(i),az(i),bz(i) 
        enddo

        ! Get fit error
        call fluxfit_error(err)
        write(4,40) rmin,err

        if (i_print == 1) then
           print '(t3,a,i4,2x,3(a,1pe12.6,3x))',&
                'Surface',p,'rmin = ',rmin,'error = ',err
           print 20,tag(:)
           do i=0,ns
              print 30,i,ar(i),br(i),az(i),bz(i)
           enddo
           print *
        endif
     else

        !----------------------------
        ! Parameterized (see f_model)
        !----------------------------

        ! Minor radius

        c(1) = rmin

        ! Elevation

        c(2) = z_c

        ! Major radius

        c(3) = rmaj

        ! Elongation

        call fluxfit_minmax(zd,rd,nd,zp,rp,'max')
        call fluxfit_minmax(zd,rd,nd,zm,rm,'min')

        c(4) = (zp-zm)/(2*rmin)

        ! Triangularity (average of upper and lower)

        c(5) = (rmaj-0.5*(rp+rm))/rmin

        ! Outer squareness (average of outer/upper and outer/lower)

        rs = rmaj+rmin*cos(0.25*pi+asin(c(5))/sqrt(2.0))

        call fluxfit_intersect(zd,rd,nd,rs,zp,zm)

        c(6) = -0.25*pi+0.5*(&
             asin((zp-z_c)/(c(4)*rmin))+asin((z_c-zm)/(c(4)*rmin)))

        ! Circle
        !c(4) = 1.0
        !c(5:6) = 0.0

        write(2,10) c(:)

        ! Get fit error
        call fluxfit_error(err)
        write(4,40) rmin,err

        if (i_print == 1) then
           print '(t3,a,i4,2x,3(a,1pe12.6,3x))',&
                'Surface',p,'rmin = ',rmin,'error = ',err
           print 20,tag(:)
           print 10,c(:)
           print *
        endif

     endif

     ! Write (R,Z) contours for plotting

     do i=1,nd
        t = (i-1.0)*2.0*pi/nd
        call fluxfit_f_model(t,r,z)
        write(3,40) r,z,rd(i),zd(i)
     enddo

  enddo

  close(1)
  close(2)
  close(3)
  close(4)


  deallocate(rd)
  deallocate(zd)
  deallocate(theta)
  deallocate(tag)
  if (model == 1) then
     deallocate(c)
  else
     deallocate(ar)
     deallocate(br)
     deallocate(az)
     deallocate(bz)
  endif

10 format(t2,10(1pe15.8,1x))
20 format(t2,10(1x,a,4x))
30 format((t2,i2,11x,4(1pe15.8,1x)))
40 format(4(1pe15.8,1x))

end subroutine fluxfit_driver
