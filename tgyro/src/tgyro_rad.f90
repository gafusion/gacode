!---------------------------------------------------------
! Portable subroutines for calculation of radiation:
!
! rad_sync -> synchrotron 
! rad_ion  -> Bremsstrahlung and impurity line 
!---------------------------------------------------------

subroutine rad_sync(aspect_rat,r_min,b_ref,ne,te,qsync,n)

  ! IN:
  !
  ! aspect_rat = R/a [-]
  ! r_min      = a   [cm]
  ! b_ref      = B   [G]
  !
  ! ne [1/cm^3]
  ! te [eV]
  !
  ! OUT:
  !
  ! qsync [erg/cm^3/s]

  implicit none

  integer, intent(in) :: n
  real, intent(in) :: aspect_rat,r_min
  real, intent(in), dimension(n) :: b_ref,ne,te
  real, intent(inout), dimension(n) :: qsync

  integer :: i
  real :: pi,e,k,me,c
  real :: g,phi,wpe,wce

  ! Reflection coefficient (Rosenbluth)
  real, parameter :: r_coeff=0.8

  !------------------------------------------------------
  ! PHYSICAL CONSTANTS
  !
  pi = 4.0*atan(1.0)
  e  = 4.8032e-10 ! statcoul
  k  = 1.6022e-12 ! erg/eV
  me = 9.1094e-28 ! g
  c  = 2.9979e10  ! cm/s
  !------------------------------------------------------

  !-------------------------------------------------------
  ! Synchrotron radiation
  ! - Trubnikov, JETP Lett. 16 (1972) 25.
  do i=1,n   
     wpe = sqrt(4*pi*ne(i)*e**2/me)
     wce = e*abs(b_ref(i))/(me*c)
     g   = k*te(i)/(me*c**2)
     phi = 60*g**1.5*sqrt((1.0-r_coeff)*(1+1/aspect_rat/sqrt(g))/(r_min*wpe**2/c/wce))

     qsync(i) = me/(3*pi*c)*g*(wpe*wce)**2*phi
  enddo
  !-------------------------------------------------------

end subroutine rad_sync

subroutine rad_ion_adas(&
     te,&
     ne,&
     ni,&
     z,&
     name,&
     qbrem,&
     qline,&
     nion,&
     n)

  ! IN:
  !
  ! te    [eV]
  ! ne    [1/cm^3]
  ! ni    [1/cm^3]
  ! z (charge) [-]
  ! name (ion names) [string] 
  ! nion (number of ions)
  ! n    (number of radii)
  !
  ! OUT:
  !
  ! qbrem  [erg/cm^3/s]
  ! qline  [erg/cm^3/s]

  implicit none

  integer, intent(in) :: nion,n
  real, intent(in), dimension(n) :: te,ne
  real, intent(in), dimension(nion,n) :: ni
  real, intent(in), dimension(nion) :: z
  character(len=3), intent(in), dimension(nion) :: name
  real, intent(inout), dimension(n) :: qbrem,qline

  ! Internals

  integer :: i
  real, dimension(n) :: lz
  real, dimension(nion,n) :: qbremi,qcool

  qbrem = 0.0
  qline = 0.0

  do i=1,nion

     ! *Approximate* NRL Bremsstrahlung radiation [erg/cm^3/s]
     ! - 1 W/cm^3 = 1e7 erg/cm^3/s

     qbremi(i,:) = 1e7*1.69e-32*ne*ni(i,:)*z(i)**2*sqrt(te)

     ! Chebyshev interpolation of ADAS data
     !  (see notes in adas21 routine)
     call adas21(te/1e3,name(i),lz,n)

     ! Total ion radiation (continuum plus line, *includes* Bremsstrahlung)
     qcool(i,:) = ne*ni(i,:)*lz

     ! Split into approximate Brem (qbrem) and approximate line (qline)
     ! where qbrem+qline is precise ADAS value
     qbrem = qbrem+qbremi(i,:)
     qline = qline+qcool(i,:)-qbremi(i,:)         

  enddo

end subroutine rad_ion_adas

!-----------------------------------------------------------------------------------
! 12-term Chebyshev polynomial fits to ADAS data
!
!    ln[ Lz(x) ] = sum c_n T_n(x)
!
! where 
!
!    Lz = cooling rate in erg/cm^3/s
!
!    T_n(x) = cos[n*arccos(x)] (Chebyshev polynomials)
!
!    T = electron temperature
! 
!                 ln(T/T_min)
!    x = -1 + 2 --------------- 
!               ln(T_max/T_min)
!
!    c_n = tabulated polynomial coefficients for each ion
!
! Acknowledgements:
! - F. Sciortino for providing access to ADAS data via Aurora 
! - T. Pütterich for up-to-data ADAS data
! - T. Odstrčil for help/checking of 2025 updates
! References:
! - Open ADAS: https://open.adas.ac.uk
! - T. Pütterich et al 2019 Nucl. Fusion 59 056013
! Notes:
! - Lz = Lz_line + Lz_continuum 
! - Aurora follows the radiation nomenclature of ADAS (as described here), separating
!   "line" and "continuum" radiation. Line radiation basically comes from ADF11 PLT
!   files and continuum radiation comes from ADF11 PRB files. Bremsstrahlung is
!   included in the continuum term.
! - For generation of fit coefficients, see tgyro/tools/radiation
!-----------------------------------------------------------------------------------

subroutine adas21(te,name,lz,n)

  implicit none

  integer, intent(in) :: n

  character(len=3), intent(in) :: name
  real, intent(in), dimension(n) :: te
  real, intent(inout), dimension(n) :: lz

  ! Number of terms in Chebyshev series
  integer, parameter :: nc=12
  ! Min and max values of Te 
  real, parameter :: t0=0.05 ! keV
  real, parameter :: t1=50.0 ! keV

  real, dimension(nc) :: c

  real :: s,x
  integer :: i,j

  if (name == 'W') then
     c(:) = (/-4.098179995282e+01,-7.936786403815e-01,-4.619044787359e-01,-1.456232613687e-01,&
              +2.937985817197e-01,+5.063695264214e-02,-1.062685541668e-01,+2.661063713322e-02,&
              +8.700125448336e-02,-6.988848399826e-02,+4.169797142892e-03,+9.826847754665e-03/)
  else if (name == 'B') then
     c(:) = (/-4.824500817201e+01,-9.583639579115e-01,+1.121730883301e+00,-1.860489734455e-01,&
              -1.336561575251e-01,+1.485164070054e-01,-9.072467667403e-02,+4.806968067742e-02,&
              -2.497690268099e-02,+1.137544241417e-02,-2.413605670955e-03,-1.124244591834e-03/)
  else if (name == 'He') then
     c(:) = (/-5.134700413297e+01,+8.218593504938e-01,+4.572258948229e-01,-1.864060260034e-01,&
              +6.186336465363e-02,-1.425022648657e-02,-3.228599626925e-04,+2.723193135411e-03,&
              -1.770156311492e-03,+5.783337251866e-04,+2.308975845735e-05,-4.276172509745e-04/)
  else if (name == 'Be') then
     c(:) = (/-4.901561102918e+01,-5.746662236007e-01,+1.102477316209e+00,-3.421385730883e-01,&
              +3.464880456112e-02,+3.087936176685e-02,-1.824676177229e-02,+2.267469723172e-03,&
              +2.551848339714e-03,-1.268897504553e-03,+2.514022413007e-04,-8.449017995079e-04/)
  else if (name == 'C') then
     c(:) = (/-4.783829635152e+01,-8.537270849071e-01,+7.710123846565e-01,+2.354004692170e-01,&
              -4.440438456161e-01,+2.927614297369e-01,-1.156983192437e-01,+1.733993408112e-02,&
              +2.170627066061e-02,-3.195887243365e-02,+2.481800173124e-02,-8.949718256775e-03/)
  else if (name == 'O') then
     c(:) = (/-4.709052058267e+01,-7.634895015203e-01,+2.965651964690e-01,+4.801687007883e-01,&
              -1.986647439460e-01,-2.151688086625e-01,+3.313680950846e-01,-2.113807929523e-01,&
              +4.535099018906e-02,+5.682622389415e-02,-7.724691172118e-02,+4.671388518804e-02/)
  else if (name == 'N') then
     c(:) = (/-4.748406078178e+01,-6.895926848662e-01,+3.403024330483e-01,+5.832113228037e-01,&
              -5.085976351090e-01,+1.385618342393e-01,+1.094185690421e-01,-1.741954140729e-01,&
              +1.367404249744e-01,-7.076042601889e-02,+1.377707920484e-02,+1.681956225002e-02/)
  else if (name == 'F') then
     c(:) = (/-4.608334126395e+01,-2.032569469519e+00,+1.174720067392e+00,-1.601584693147e-01,&
              +2.522457754350e-01,-4.105666049867e-01,+2.639488083059e-01,-1.240057037307e-04,&
              -1.483369633643e-01,+1.362062267100e-01,-4.689893298094e-02,-2.644300012089e-02/)
  else if (name == 'Ne') then
     c(:) = (/-4.616289844396e+01,-1.476875032140e+00,+8.751481882578e-01,-1.554907576477e-01,&
              +3.453234274569e-01,-4.448210670608e-01,+2.235707626168e-01,+6.309160185099e-02,&
              -1.804901948539e-01,+1.138691735006e-01,+1.238414608039e-02,-7.946295188381e-02/)
  else if (name == 'Al') then
     c(:) = (/-4.486524782614e+01,-2.322475619179e+00,+9.328274092632e-01,+4.744804675756e-02,&
              -6.858095854371e-02,+2.345581654349e-01,-3.716843250349e-01,+1.795673265946e-01,&
              +1.194090655176e-01,-2.087103454071e-01,+8.992619526224e-02,+3.923399787447e-02/)
  else if (name == 'Si') then
     c(:) = (/-4.473392580102e+01,-2.118859872224e+00,+7.535819784959e-01,+1.199722562015e-01,&
              -1.300605991542e-01,+3.223304318098e-01,-3.674984537405e-01,+5.351426696879e-02,&
              +2.108242428513e-01,-1.792043287525e-01,+9.913893641685e-03,+7.490321841247e-02/)
  else if (name == 'Ar') then
     c(:) = (/-4.425505195350e+01,-1.616269823924e+00,+8.194549090984e-02,+4.830273331631e-01,&
              -2.675106132694e-01,+1.174875472786e-01,+1.676830849068e-01,-2.891215850568e-01,&
              +4.842130798317e-02,+1.715110727683e-01,-1.077751857577e-01,-4.574007342660e-02/)
  else if (name == 'Ca') then
     c(:) = (/-4.410322851491e+01,-1.396392282770e+00,+6.228605606294e-02,+2.603343415878e-01,&
              +6.868991709327e-02,-3.101682260197e-01,+5.258665807121e-01,-3.702755461767e-01,&
              -1.771137529772e-02,+1.469773755999e-01,-1.213876050533e-02,-5.491593884657e-02/)
  else if (name == 'Fe') then
     c(:) = (/-4.289975750516e+01,-2.074724550183e+00,+2.592709540821e-01,+2.670554581981e-01,&
              -3.245193180198e-02,-3.447451297355e-02,-1.365795512615e-01,+2.145851896399e-01,&
              +2.140105004146e-02,-1.785057925001e-01,+5.091860663164e-02,+6.171044344444e-02/)
  else if (name == 'Ni') then
     c(:) = (/-4.281682927312e+01,-1.989503354777e+00,+4.039460832605e-01,+2.224196803176e-01,&
              -1.451323786293e-01,+7.644152992866e-02,-1.717626198140e-01,+1.431632992224e-01,&
              +1.513611521059e-01,-2.308718720128e-01,-9.023390843523e-03,+1.781391433306e-01/)
  else if (name == 'Kr') then
     c(:) = (/-4.248982974133e+01,-1.299186536746e+00,-4.512694457002e-01,+6.804566946335e-01,&
              -6.310168823628e-03,-1.033700187711e-01,+2.624047807761e-02,-1.531016488962e-01,&
              +1.229092533413e-01,+5.721604110683e-02,-3.858573093198e-02,-7.656826000830e-02/)
  else if (name == 'Mo') then
     c(:) = (/-4.192185652052e+01,-1.792215448077e+00,+2.732017007262e-03,+1.000290040334e-01,&
              +4.133083599928e-01,-1.495307361845e-01,-1.026532480504e-01,+1.359298888695e-01,&
              -1.749938253055e-01,+7.007588723552e-02,+8.557768114101e-02,-7.082562416821e-02/)
  else if (name == 'Xe') then
     c(:) = (/-4.136087494072e+01,-1.651105892275e+00,-2.928373643947e-01,+2.928170785526e-01,&
              -5.910859658247e-02,+1.267994598292e-01,+2.256208498121e-02,-4.369518838603e-02,&
              +1.738552444143e-02,+3.245918863202e-03,-6.535620998184e-02,+3.866947213111e-03/)
  else if (name == 'Li') then
     c(:) = (/-5.007779179318e+01,+1.379180885582e-01,+7.886186847887e-01,-2.784155208154e-01,&
              +5.393512924579e-02,+8.568540178773e-03,-8.339965761965e-03,-3.235331808761e-03,&
              +7.058675207572e-03,-3.451108596954e-03,-2.663068572132e-03,+6.045007428451e-03/)
  else if (name == 'H' .or. name == 'D' .or. name == 'T') then
     c(:) = (/-5.313111260147e+01,+1.418007193691e+00,+1.261973276523e-01,-4.362701135592e-02,&
              +1.071044131666e-02,-3.046777668538e-03,+7.410284109890e-04,+5.926084113166e-05,&
              +4.909076663983e-05,-5.820357621765e-05,+8.674684934388e-05,-9.027712584843e-05/)
  else 

     ! ion not found
     lz(:) = 0.0
     return

  endif
  
  do j=1,n
     ! Convert Te to Chebyshev grid x = [-1,1]
     x = -1.0+2*log(te(j)/t0)/log(t1/t0)
     ! If outside domain use boundary value
     if (x < -1.0) x = -1.0
     if (x > +1.0) x = +1.0
     ! Sum the Chebyshev series where T_n(x) = cos[n*arccos(x)]
     s = 0.0
     do i=1,nc
        s = s+c(i)*cos((i-1)*acos(x))
     enddo
     ! Lz (cooling rate) in erg/cm^3/s
     lz(j) = exp(s)
  enddo
  
end subroutine adas21
