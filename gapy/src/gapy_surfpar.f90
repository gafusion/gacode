module surfpar

  integer :: iselect=16
  integer :: nf=2
  double precision, dimension(:), allocatable :: r,z
  double precision, dimension(:), allocatable :: vr,vz
  double precision, dimension(:), allocatable :: pr,pz
  double precision, dimension(:), allocatable :: x
  double precision, dimension(:), allocatable :: rp,zp

contains

  ! ** Main callable routine **
  
  subroutine surfpar_do(n_psi,n_arc,r_in,z_in)

    !-------------------------------------
    implicit none
    !
    integer, intent(in) :: n_psi,n_arc
    double precision, intent(in), dimension(n_psi,n_arc) :: r_in,z_in


    integer :: i,ix,p
    integer :: n1,m1,n2,m2,s
    double precision, dimension(:), allocatable :: dl,ur,uz,l
    double precision, dimension(:), allocatable :: cr,cz,sr,sz
    double precision :: rmin,rmaj,zmin,zmaj,eps,pi
    !-------------------------------------

    pi = 4.0*atan(1.0)

    allocate(r(n_arc))
    allocate(z(n_arc))
    allocate(dl(n_arc))
    allocate(ur(n_arc))
    allocate(uz(n_arc))
    allocate(vr(n_arc))
    allocate(vz(n_arc))
    allocate(pr(n_arc))
    allocate(pz(n_arc))
    allocate(rp(n_arc))
    allocate(zp(n_arc))
    allocate(l(n_arc))
    allocate(x(n_arc))
    allocate(sr(nf))
    allocate(sz(nf))
    allocate(cr(0:nf))
    allocate(cz(0:nf))

    r = r_in(iselect,:)
    z = z_in(iselect,:)

    ! Pointwise Extrema ------------------------------
    n1 = maxloc(z,1) ; m1 = minloc(z,1)
    zmaj = 0.5*(z(n1)+z(m1))
    zmin = 0.5*(z(n1)-z(m1))

    n2 = maxloc(r,1) ; m2 = minloc(r,1)
    rmaj = 0.5*(r(n2)+r(m2))
    rmin = 0.5*(r(n2)-r(m2))

    s = n2

    ! Shift elements so that first index at max(R).
    r(1:n_arc-1) = cshift(r(1:n_arc-1),s)
    z(1:n_arc-1) = cshift(z(1:n_arc-1),s)
    r(n_arc) = r(1)
    z(n_arc) = z(1)

    if (z(2) < z(1)) then
       ! Reverse order (may be needed)
       r = r(n_arc:1:-1)
       z = z(n_arc:1:-1)
    endif

    ! Compute arc lengths
    do i=1,n_arc-1
       dl(i) = sqrt((r(i+1)-r(i))**2+(z(i+1)-z(i))**2)
    enddo
    dl(n_arc) = dl(1)

    !---------------------------------------------------
    ! Arc length summation 
    l(1) = 0.0
    do i=1,n_arc-1
       l(i+1) = l(i)+dl(i)
    enddo
    x = 2*pi*l/l(n_arc)

    ! Compute generalized angles ---------------------
    eps = 1.0-1e-9
    uz = asin(eps*(z-zmaj)/zmin)
    ur = acos(eps*(r-rmaj)/rmin)

    ! Sort out proper branches
    vr = ur ; vz = uz
    do i=1,n_arc-2
       if (ur(i+1) < ur(i)) vr(i+1) = 2*pi-ur(i+1)
       if (uz(i+1) <= uz(i)) vz(i+1) = pi-uz(i+1)
       if (uz(i+1) >= uz(i) .and. x(i) > pi) then
          vz(i+1) = 2*pi+uz(i+1)
       endif
    enddo

    vr = vr-x
    vz = vz-x
    vr(n_arc) = vr(1)
    vz(n_arc) = vz(1)

    cr = 0.0 ; cz = 0.0
    cr(0) = moment(n_arc,vr,x*0.0+1,dl)
    cz(0) = moment(n_arc,vz,x*0.0+1,dl)

    do p=1,nf
       cr(p) = moment(n_arc,vr,cos(p*x),dl)
       cz(p) = moment(n_arc,vz,cos(p*x),dl)
       sr(p) = moment(n_arc,vr,sin(p*x),dl)
       sz(p) = moment(n_arc,vz,sin(p*x),dl)
    enddo

    ! Regenerage based on shape param
    pr = cr(0) ; pz = cz(0)
    do p=1,nf
       pr = pr+cr(p)*cos(p*x)+sr(p)*sin(p*x)
       pz = pz+cz(p)*cos(p*x)+sz(p)*sin(p*x)
    enddo

    rp = rmaj+rmin*cos(x+pr)
    zp = zmaj+zmin*sin(x+pz)

  end subroutine surfpar_do

  double precision function moment(n,f,w,d)

    implicit none

    integer, intent(in) :: n
    double precision, intent(in) :: f(n),w(n),d(n)
    integer :: i
    double precision s0,s1

    s0 = 0.0
    s1 = 0.0
    do i=1,n-1
      s0 = s0+0.5*(f(i)*w(i)+f(i+1)*w(i+1))*d(i)
      s1 = s1+0.5*(w(i)**2+w(i+1)**2)*d(i)
    enddo

   moment = s0/s1

 end function moment
  
end module surfpar
