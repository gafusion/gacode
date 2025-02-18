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
! Chebyshev polynomial fits to ADAS data
!
!    ln[ Lz(x) ] = sum c_n T_n(x)
!
!  where
! 
!                 ln(T/T_min)
!    x = -1 + 2 ---------------
!               ln(T_max/T_min)
!
! Acknowledgements:
!  - F. Sciortino for providing access to ADAS data via Aurora 
!  - T. Pütterich for up-to-data ADAS data
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
!-------------------------------------------------------------------------------------

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

  ! Chebyshev expansion coefficients: c_n 
  
  if (name == 'W') then
     c(:) = (/-4.093426327035e+01,-8.887660631564e-01,-3.780990284830e-01,-1.950023337795e-01,&
              +3.138290691843e-01,+4.782989513315e-02,-9.942946187466e-02,+8.845089763161e-03,&
              +9.069526573697e-02,-5.245048352825e-02,-1.487683353273e-02,+1.917578018825e-02/)
  else if (name == 'Xe') then
     c(:) = (/-4.126366679797e+01,-1.789569183388e+00,-2.380331458294e-01,+2.916911530426e-01,&
              -6.217313390606e-02,+1.177929596352e-01,+3.114580325620e-02,-3.551020007260e-02,&
              -4.850122964780e-03,+1.132323304719e-02,-5.275312157892e-02,-9.051568201374e-03/)
  else if (name == 'Mo') then
     c(:) = (/-4.178151951275e+01,-1.977018529373e+00,+5.339155696054e-02,+1.164267551804e-01,&
              +3.697881990263e-01,-9.594816048640e-02,-1.392054581553e-01,+1.272648056277e-01,&
              -1.336366483240e-01,+3.666060293888e-02,+9.586025795242e-02,-7.210209944439e-02/)
  else if (name == 'Kr') then
     c(:) = (/-4.235332287815e+01,-1.508707679199e+00,-3.300772886398e-01,+6.166385849657e-01,&
              +1.752687990068e-02,-1.004626261246e-01,+5.175682671490e-03,-1.275380183939e-01,&
              +1.087790584052e-01,+6.846942959545e-02,-5.558980841419e-02,-6.669294912560e-02/)
  else if (name == 'Ni') then
     c(:) = (/-4.269403899818e+01,-2.138567547684e+00,+4.165648766103e-01,+2.507972619622e-01,&
              -1.454986877598e-01,+4.044612562765e-02,-1.231313167536e-01,+1.307076922327e-01,&
              +1.176971646853e-01,-1.997449027896e-01,-8.027057678386e-03,+1.583614529900e-01/)
  else if (name == 'Fe') then
     c(:) = (/-4.277490044241e+01,-2.232798257858e+00,+2.871183684045e-01,+2.903760139426e-01,&
              -4.662374777924e-02,-4.436273974526e-02,-1.004882554335e-01,+1.794710746088e-01,&
              +3.168699330882e-02,-1.813266337535e-01,+5.762415716395e-02,+6.379542965373e-02/)
  else if (name == 'Ca') then
     c(:) = (/-4.390083075521e+01,-1.692920511934e+00,+1.896825846094e-01,+2.333977195162e-01,&
              +5.307786998918e-02,-2.559420140904e-01,+4.733492400000e-01,-3.788430571182e-01,&
              +3.375702537147e-02,+1.030183684347e-01,+1.523656115806e-02,-7.482021324342e-02/)
  else if (name == 'Ar') then
     c(:) = (/-4.412345259739e+01,-1.788450950589e+00,+1.322515262175e-01,+4.876947389538e-01,&
              -2.869002749245e-01,+1.699452914498e-01,+9.950501421570e-02,-2.674585184275e-01,&
              +7.451345261250e-02,+1.495713760953e-01,-1.089524173155e-01,-4.191575231760e-02/)
  else if (name == 'Si') then
     c(:) = (/-4.459983387390e+01,-2.279998599897e+00,+7.703525425589e-01,+1.494919348709e-01,&
              -1.136851457700e-01,+2.767894295326e-01,-3.577491771736e-01,+7.013841334798e-02,&
              +2.151919651291e-01,-2.052895326141e-01,+2.210085804088e-02,+9.270982150548e-02/)
  else if (name == 'Al') then
     c(:) = (/-4.475065090279e+01,-2.455868594007e+00,+9.468903008039e-01,+6.944445017599e-02,&
              -4.550919134508e-02,+1.804382546971e-01,-3.573462505157e-01,+2.075274089736e-01,&
              +1.024482383310e-01,-2.254367207993e-01,+1.150695613575e-01,+3.414328980459e-02/)
  else if (name == 'Ne') then
     c(:) = (/-4.599844680574e+01,-1.684860164232e+00,+9.039325377493e-01,-7.544604235334e-02,&
              +2.849631706915e-01,-4.827471944126e-01,+3.138177972060e-01,+2.876874062690e-03,&
              -1.809607030192e-01,+1.510609882754e-01,-2.475867654255e-02,-6.269602018004e-02/)
  else if (name == 'F') then
     c(:) = (/-4.595870691474e+01,-2.176917325041e+00,+1.176783264877e+00,-7.712313240060e-02,&
              +1.847534287214e-01,-4.297192280031e-01,+3.374503944631e-01,-5.862051731844e-02,&
              -1.363051725174e-01,+1.580531615737e-01,-7.677594113938e-02,-5.498186771891e-03/)
  else if (name == 'N') then
     c(:) = (/-4.719917668483e+01,-1.128938430123e+00,+5.686617156868e-01,+5.565647850806e-01,&
              -6.103105546858e-01,+2.559496676285e-01,+3.204394187397e-02,-1.347036917773e-01,&
              +1.166192946931e-01,-6.001774708924e-02,+1.078186024405e-02,+1.336864982060e-02/)
  else if (name == 'O') then
     c(:) = (/-4.688092238361e+01,-1.045540847894e+00,+3.574644442831e-01,+6.007860794100e-01,&
              -3.812470436912e-01,-9.944716626912e-02,+3.141455586422e-01,-2.520592337580e-01,&
              +9.745206757309e-02,+1.606664371633e-02,-5.269687016804e-02,+3.726780755484e-02/)
  else if (name == 'C') then
     c(:) = (/-4.752370087442e+01,-1.370806613078e+00,+1.119762977201e+00,+6.244262441360e-02,&
              -4.172077577493e-01,+3.237504483005e-01,-1.421660253114e-01,+2.526893756273e-02,&
              +2.320010310338e-02,-3.487271688767e-02,+2.758311539699e-02,-1.063180164276e-02/)
  else if (name == 'Be') then
     c(:) = (/-4.883447566291e+01,-8.543314577695e-01,+1.305444973614e+00,-4.830394934711e-01,&
              +1.005512839480e-01,+1.392590190604e-02,-1.980609625444e-02,+5.342857189984e-03,&
              +2.324970825974e-03,-2.466382923947e-03,+1.073116177574e-03,-9.834117466066e-04/)
  else if (name == 'Li') then
     c(:) = (/-4.978562496154e+01,-2.505216545881e-01,+8.650665756334e-01,-1.291711692636e-01,&
              -1.599934526332e-01,+1.322662928235e-01,-5.487930945808e-03,-6.479708903897e-02,&
              +4.332749716030e-02,+1.391112355350e-02,-3.887843175798e-02,+1.377942187934e-02/)
  else if (name == 'He') then
     c(:) = (/-5.128490291648e+01,+7.743125302555e-01,+4.674917416545e-01,-2.087203609904e-01,&
              +7.996303682551e-02,-2.450841492530e-02,+4.177032799848e-03,+1.109529527611e-03,&
              -1.080271138220e-03,+1.914061606095e-04,+2.501544833223e-04,-3.856698155759e-04/)
  else if (name == 'H' .or. name == 'D' .or. name == 'T') then
     ! Hydrogen-like ions (H,D,T)
     c(:) = (/-5.307012989032e+01,+1.382271913121e+00,+1.111772196884e-01,-3.989144654893e-02,&
              +1.043427394534e-02,-3.038480967797e-03,+5.851591993347e-04,+3.472228652286e-04,&
              -8.418918897927e-05,+3.973067124523e-05,-3.853620366361e-05,+2.005063821667e-04/)
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
