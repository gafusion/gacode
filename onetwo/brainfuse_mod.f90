MODULE BRAINFUSE_MOD

  INTEGER include_brainfuse
  
CONTAINS
!========================
!========================
  FUNCTION to_lower(strIn) result(strOut)
    IMPLICIT NONE
    CHARACTER(len=*), intent(in) :: strIn
    CHARACTER(len=len(strIn)) :: strOut
    INTEGER i,j
    do i = 1, len(strIn)
       j = iachar(strIn(i:i))
       if (j>= iachar("A") .and. j<=iachar("Z") ) then
          strOut(i:i) = achar(iachar(strIn(i:i))+32)
       else
          strOut(i:i) = strIn(i:i)
       end if
    end do
  END FUNCTION to_lower
!========================
!========================
  SUBROUTINE GRADIENT(nn, y, yy)
    IMPLICIT NONE
    INTEGER*4 nn, j
    REAL*8 y(nn), yy(nn)
    yy(1) = y(2)-y(1)
    yy(2:nn-1) = (y(3:nn)-y(1:nn-2))/2.
    yy(nn) = y(nn)-y(nn-1)
  END SUBROUTINE GRADIENT
!========================
!========================
  SUBROUTINE brainfuse(nn, bt, ip, r, rmaj, kappa, ne, ni, te, ti, q, vol, wt, press, Qe, Qi)
    IMPLICIT NONE
    INTEGER*4 debug, nn, j, i, dummy
    PARAMETER (debug=1)
    REAL*8  bt, ip, r(nn), rmaj(nn), kappa(nn), delta(nn), ne(nn), ni(nn), te(nn), ti(nn), q(nn), vol(nn), wt(nn), press(nn), Qe(nn), Qi(nn)
    REAL*8, DIMENSION (:), ALLOCATABLE :: cs, bunit, dte, dti, dne, dni, dr, dq, dkappa, ddelta, dvol, dwt, rhos, dpress
    REAL*8 a, mD, eV, pi
    PARAMETER (mD = 3.3475E-27)
    PARAMETER (eV = 1.60217646E-19)
    PARAMETER (pi = 3.14159265359)
!        _mE = 9.10938291E-31
!        _c = 3E8
!        _e0 = 8.8541878176E-12
!        _Kb = 1.3806505E-23
    CHARACTER input_names(20)*255
    CHARACTER output_names(20)*255

    REAL*4, DIMENSION (:,:), ALLOCATABLE :: input, output
    INTEGER num_data, num_input, num_output

    IF (debug) WRITE(*,*)'Running brainfuse!'

!========================
    input_names(1)='rmin_loc'
    input_names(2)='rmaj_loc'
    input_names(3)='kappa'
    input_names(4)='taus'
    input_names(5)='aus'
    input_names(6)='q'

    output_names(1)='Qe'
    output_names(2)='Qi'
!========================

    a=r(nn)

    ALLOCATE( cs(nn) )
    ALLOCATE( rhos(nn) )
    ALLOCATE( bunit(nn) )
    ALLOCATE( dr(nn) )
    ALLOCATE( dte(nn) )
    ALLOCATE( dti(nn) )
    ALLOCATE( dne(nn) )
    ALLOCATE( dni(nn) )
    ALLOCATE( dq(nn) )
    ALLOCATE( dkappa(nn) )
    ALLOCATE( ddelta(nn) )
    ALLOCATE( dvol(nn) )
    ALLOCATE( dwt(nn) )
    ALLOCATE( dpress(nn) )

    cs = SQRT(1E3*eV*te/mD)
    rhos=1.022e-4*SQRT(1E3*eV*ti)/mD/Bt/a
    bunit = cs*mD/(eV*rhos*a)

    call GRADIENT(nn,r,dr)
    call GRADIENT(nn,te,dte)
    dte=dte/dr
    call GRADIENT(nn,ti,dti)
    dti=dti/dr
    call GRADIENT(nn,ne,dne)
    dne=dne/dr
    call GRADIENT(nn,ni,dni)
    dni=dni/dr
    call GRADIENT(nn,q,dq)
    dq=dq/dr
    call GRADIENT(nn,kappa,dkappa)
    dkappa=dkappa/dr
    call GRADIENT(nn,delta,ddelta)
    ddelta=ddelta/dr
    call GRADIENT(nn,vol,dvol)
    dvol=dvol/dr
    call GRADIENT(nn,press,dpress)
    dpress=dpress/dr

!========================

    num_input=6
    num_output=2
    ALLOCATE( input(nn,num_input)   )
    ALLOCATE( output(nn,num_output) )
    DO j=1,num_input
       dummy=1
       SELECT CASE (to_lower(trim(input_names(j))))
       CASE('r')
          input(:,j) = r                       !minor radius
       CASE('rmaj')
          input(:,j) = rmaj                    !major radius
       CASE('te')
          input(:,j) = te                      !electron temperature
       CASE('dte')
          input(:,j) = dte                     !electron temperature gradient
       CASE('ti')
          input(:,j) = ti                      !main ion temperature
       CASE('dti')
          input(:,j) = dti                     !main ion temperature gradient
       CASE('ne')
          input(:,j) = ne                      !electron density
       CASE('dne')
          input(:,j) = dne                      !electron density gradient
       CASE('ni')
          input(:,j) = ni                      !main ion density
       CASE('dni')
          input(:,j) = dni                      !main ion density gradient
       CASE('cs')
          input(:,j) = cs                      !ion sound speed
       CASE('vol')
          input(:,j) = vol                     !volume
       CASE('dvol')
          input(:,j) = dvol                    !volume gradient
       CASE('press')
          input(:,j) = press                   !total pressure
       CASE('dpress')
          input(:,j) = dpress                  !total pressure gradient
       CASE('rmin_loc')
          input(:,j) = r/a                     !normalized minor radius
       CASE('rmaj_loc')
          input(:,j) = rmaj/a                  !normalized major radius
       CASE('kappa')
          input(:,j) = kappa                   !elongation
       CASE('dkappa')
          input(:,j) = dkappa                  !elongation gradient
       CASE('s_kappa_loc')
          input(:,j) = r*dkappa                !normalized elongation shear
       CASE('delta')
          input(:,j) = delta                   !triangularity
       CASE('ddelta')
          input(:,j) = ddelta                  !triangularity gradient
       CASE('s_delta_loc')
          input(:,j) = r*ddelta                !normalized triangularity shear
       CASE('taus')
          input(:,j) = ti/te                   !temperature ratio
       CASE('aus')
          input(:,j) = ni/ne                   !density ratio
       CASE('q')
          input(:,j) = q                       !safety factor
       CASE('dq')
          input(:,j) = dq                      !safety factor shear
       CASE('q_prime_loc')
          input(:,j) = r*dq                    !normalized safety factor shear
       CASE('gamma_p')
          input(:,j) = -rmaj*dwt               !gamma P
       CASE('gamma_e')
          input(:,j) = -dwt*r/q                !gamma E
       CASE('vpar')
          input(:,j) = -ABS(Ip)/Ip*rmaj*wt/cs  !vpar
       CASE('vpar_shear')
          input(:,j) = -ABS(Ip)/Ip*(-rmaj*dwt)*r/cs  !vpar shear
       CASE('vexb_shear')
          input(:,j) = -ABS(Ip)/Ip*(-dwt*r/q)*r/cs   !vper shear
       CASE('betae')
          input(:,j) = 8*pi*ne*te/bunit**2     !electron beta
       CASE('xnue')
          input(:,j) = 1.33E5*(ne/1E20)/te**1.5 * a/cs     !normalized ei collisionality
       CASE('debye')
          input(:,j) = 2.35E-5*SQRT(te/(ne/1E20)) / (rhos*a)     !normalized Debye length
       CASE('qgb')
          input(:,j) = ne*cs*(te*1E3*eV)*(rhos/a)**2     !gyrobohm flux
       CASE('lte','rlts_1')
          input(:,j) = -a*dte/te               !electron temperature scale length
       CASE('lti','rlts_2')
          input(:,j) = -a*dti/ti               !ion temperature scale length
       CASE('lne','rlns_1')
          input(:,j) = -a*dne/ne               !electron density scale length
       CASE('lni','rlns_2')
          input(:,j) = -a*dni/ni               !ion density scale length
       CASE('p_prime_loc','lp','rlps')
          input(:,j) = -a*dpress/press         !total pressure scale length
       CASE DEFAULT
          dummy=0
       END SELECT
       IF (dummy) THEN
          IF (debug) WRITE(*,*)'BRAINFUSE input variable: ',trim(input_names(j))
       ELSE
          WRITE(*,*)'ERROR in BRAINFUSE: input variable`',trim(input_names(j)),'` is not defined!'
          STOP
       ENDIF
    ENDDO
!========================

    CALL run_net_on_data(nn, num_input, num_output, input, output, debug)

!========================

    DO j=1,num_output
       dummy=1
       SELECT CASE (to_lower(trim(output_names(j))))
       CASE('qe')
          Qe=output(:,j)
       CASE('qe_qnorm')
          output(:,j)=output(:,j)*dvol
          Qe=output(:,j)
       CASE('qi')
          Qi=output(:,j)
       CASE('qi_qnorm')
          output(:,j)=output(:,j)*dvol
          Qi=output(:,j)
       CASE DEFAULT
          dummy=0
       END SELECT
       IF (dummy) THEN
          IF (debug) WRITE(*,*)'BRAINFUSE output variable: ',trim(output_names(j))
       ELSE
          WRITE(*,*)'ERROR in BRAINFUSE: output variable`',trim(output_names(j)),'` is not defined!'
          STOP
       ENDIF
    ENDDO

!========================

    IF (debug) THEN
       OPEN(unit=17,file='brainfuse_in.dat',FORM='FORMATTED',STATUS='REPLACE')
       WRITE(17,'(100(2X,A15))')(input_names(j),j=1,num_input)
       DO i=1,nn
          WRITE(17,'(100(1X,e16.9))')(input(i,j),j=1,num_input)
       ENDDO
       CLOSE(17)
       OPEN(unit=17,file='brainfuse_out.dat',FORM='FORMATTED',STATUS='REPLACE')
       WRITE(17,'(100(2X,A15))')(output_names(j),j=1,num_output)
       DO i=1,nn
          WRITE(17,'(100(1X,e16.9))')(output(i,j),j=1,num_output)
       ENDDO
       CLOSE(17)
    ENDIF

!========================

    DEALLOCATE( cs )
    DEALLOCATE( rhos )
    DEALLOCATE( bunit )
    DEALLOCATE( dr )
    DEALLOCATE( dte )
    DEALLOCATE( dti )
    DEALLOCATE( dne )
    DEALLOCATE( dni )
    DEALLOCATE( dq )
    DEALLOCATE( dkappa )
    DEALLOCATE( ddelta )
    DEALLOCATE( dvol )
    DEALLOCATE( dwt )
    DEALLOCATE( dpress )
    DEALLOCATE( input )
    DEALLOCATE( output )
    
  end subroutine brainfuse
  
end module

