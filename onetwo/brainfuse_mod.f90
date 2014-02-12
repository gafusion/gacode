MODULE BRAINFUSE_MOD

  INTEGER include_brainfuse
  CHARACTER*256 brainfuse_path
  
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
  SUBROUTINE brainfuse(nn, bt, ip, r, rmaj, kappa, ne, ni, te, ti, q, vol, wt, press, Qe, Qi, Ge, Gi, Pr)
    IMPLICIT NONE
    INTEGER*4 debug, nn, j, i, dummy, num_extra
    PARAMETER (debug=1)
    REAL*8  bt, ip, r(nn), rmaj(nn), kappa(nn), delta(nn), ne(nn), ni(nn), te(nn), ti(nn), q(nn), vol(nn), wt(nn), press(nn), Qe(nn), Qi(nn), Ge(nn), Gi(nn), Pr(nn)
    REAL*8, DIMENSION (:), ALLOCATABLE :: cs, bunit, dte, dti, dne, dni, dr, dq, dkappa, ddelta, dvol, dwt, rg, dpress
    REAL*8 a, mD, eV, pi
    PARAMETER (mD = 3.3475E-27)
    PARAMETER (eV = 1.60217646E-19)
    PARAMETER (pi = 3.14159265359)

    INTEGER slogi, slogo, logo
    CHARACTER input_names(50)*255
    CHARACTER output_names(50)*255
    CHARACTER dbgtbl_names(100)*255
    REAL*8, DIMENSION (:,:), ALLOCATABLE :: input, output, dbgtbl
    INTEGER num_input, num_output

!========================

    IF (debug) WRITE(*,*)'Running brainfuse!'

    num_input=19
    num_output=5
    input_names(1 )='RMIN_LOC'
    input_names(2 )='RMAJ_LOC'
    input_names(3 )='kappa'
    input_names(4 )='S_KAPPA_LOC'
    input_names(5 )='RLNS_1'
    input_names(6 )='RLNS_2'
    input_names(7 )='RLTS_1'
    input_names(8 )='RLTS_2'
    input_names(9 )='P_PRIME_LOC'
    input_names(10)='XNUE'
    input_names(11)='DEBYE'
    input_names(12)='VPAR'
    input_names(13)='VEXB_SHEAR'
    input_names(14)='VPAR_SHEAR'
    input_names(15)='TAUS'
    input_names(16)='AS'
    input_names(17)='q'
    input_names(18)='Q_PRIME_LOC'
    input_names(19)='RHOS'
    output_names(1)='Qe_norm'
    output_names(2)='Qi_norm'
    output_names(3)='Ge_norm'
    output_names(4)='Gi_norm'
    output_names(5)='Pr_norm'
    slogi = 0
    slogo = 0
    logo  = 0

!========================

    a=r(nn)
    !q=q*SIGN(Ip*0+1,Ip)*SIGN(Bt*0+1,Bt)

    ALLOCATE( cs(nn) )
    cs = SQRT(1E3*eV*te/mD)

    ALLOCATE( rg(nn) )
    rg=mD*SQRT(2.*1E3*eV*ti/mD)/(eV*ABS(Bt))

    ALLOCATE( bunit(nn) )
    bunit = cs*mD/(eV*rg)

    ALLOCATE( dr(nn) )
    call GRADIENT(nn,r,dr)

    ALLOCATE( dte(nn) )
    call GRADIENT(nn,te,dte)
    dte=dte/dr

    ALLOCATE( dti(nn) )
    call GRADIENT(nn,ti,dti)
    dti=dti/dr

    ALLOCATE( dne(nn) )
    call GRADIENT(nn,ne,dne)
    dne=dne/dr

    ALLOCATE( dni(nn) )
    call GRADIENT(nn,ni,dni)
    dni=dni/dr

    ALLOCATE( dq(nn) )
    call GRADIENT(nn,q,dq)
    dq=dq/dr

    ALLOCATE( dkappa(nn) )
    call GRADIENT(nn,kappa,dkappa)
    dkappa=dkappa/dr

    ALLOCATE( ddelta(nn) )
    call GRADIENT(nn,delta,ddelta)
    ddelta=ddelta/dr

    ALLOCATE( dvol(nn) )
    call GRADIENT(nn,vol,dvol)
    dvol=dvol/dr

    ALLOCATE( dwt(nn) )
    call GRADIENT(nn,wt,dwt)
    dwt=dwt/dr

    ALLOCATE( dpress(nn) )
    call GRADIENT(nn,press,dpress)
    dpress=dpress/dr

!========================

    ALLOCATE( input(nn,num_input)   )
    input=0.
    ALLOCATE( output(nn,num_output) )
    output=0.

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
          input(:,j) = dne                     !electron density gradient
       CASE('ni')
          input(:,j) = ni                      !main ion density
       CASE('dni')
          input(:,j) = dni                     !main ion density gradient
       CASE('wt')
          input(:,j) = wt                      !angular velocity
       CASE('dwt')
          input(:,j) = dwt                     !angular velocity gradient
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
       CASE('as')
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
          input(:,j) = -SIGN(Ip*0+1,Ip)*rmaj*wt/cs       !vpar
       CASE('vpar_shear')
          input(:,j) = -SIGN(Ip*0+1,Ip)*(-rmaj*dwt)*r/cs !vpar shear
       CASE('vexb_shear')
          input(:,j) = -SIGN(Ip*0+1,Ip)*(-dwt*r/q)*r/cs  !vper shear
       CASE('betae')
          input(:,j) = 8*pi*ne*te/bunit**2               !electron beta
       CASE('xnue')
          input(:,j) = 1.33E5*(ne/1E20)/te**1.5*a/cs     !normalized ei collisionality
       CASE('debye')
          input(:,j) = 2.35E-5*SQRT(te/(ne/1E20))/rg     !normalized Debye length
       CASE('qgb')
          input(:,j) = ne*cs*(te*1E3*eV)*rg**2           !gyrobohm flux
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
       CASE('rhos')
          input(:,j) = rg/a                    !rhostar
       CASE('rhos_loc')
          input(:,j) = rg/r                    !local rhostar
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
    
    DO j=1,num_input
       IF (slogi.eq.1) THEN
          input(:,j)=SIGN(input(:,j)*0+1,input(:,j))*(LOG10(ABS(input(:,j))+1)-ALOG10(1.))
       ENDIF
    ENDDO

    CALL run_net_on_data(nn, num_input, num_output, input, output, trim(brainfuse_path)//CHAR(0), debug)

    DO j=1,num_input
       IF (slogi.eq.1) THEN
          input(:,j)=SIGN(input(:,j)*0+1,input(:,j))*(10**(ABS(input(:,j))+LOG10(1.))-1)
       ENDIF
    ENDDO
    DO j=1,num_output
       IF (slogo.eq.1) THEN
          output(:,j)=SIGN(output(:,j)*0+1,output(:,j))*(10**(ABS(output(:,j))+LOG10(1.))-1)
       ENDIF
       IF (logo) THEN 
          output(:,j)=10**output(:,j)
       ENDIF
    ENDDO

!========================

    DO j=1,num_output
       dummy=1
       SELECT CASE (to_lower(trim(output_names(j))))

       CASE('qe')
          Qe=output(:,j)
       CASE('qe_norm')
          Qe=output(:,j)*dvol

       CASE('qi')
          Qi=output(:,j)
       CASE('qi_norm')
          Qi=output(:,j)*dvol

       CASE('ge')
          Ge=output(:,j)
       CASE('ge_norm')
          Ge=output(:,j)*dvol

       CASE('gi')
          Gi=output(:,j)
       CASE('gi_norm')
          Gi=output(:,j)*dvol

       CASE('pr')
          Pr=output(:,j)
       CASE('pr_norm')
          Pr=output(:,j)*dvol

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
       
       num_extra=2
       ALLOCATE( dbgtbl(nn,num_input+num_output+num_extra) ) !input+output+dvol+bt

       !input
       DO j=1,num_input
          dbgtbl_names(j)=input_names(j)
          dbgtbl(:,j)=input(:,j)
       ENDDO
       !output
       DO j=1,num_output
          dbgtbl_names(num_input+j)=output_names(j)
          dbgtbl(:,num_input+j)=output(:,j)
       ENDDO
       !extra
       dbgtbl_names(num_input+num_output+1)='dvol'
       dbgtbl(:,num_input+num_output+1)=dvol
       dbgtbl_names(num_input+num_output+2)='Bt'
       dbgtbl(:,num_input+num_output+2)=bt

       OPEN(unit=17,file='brainfuse.dat',FORM='FORMATTED',STATUS='REPLACE')
       WRITE(17,'(100(2X,A15))')(dbgtbl_names(j),j=1,num_input+num_output+num_extra)
       DO i=1,nn
          WRITE(17,'(100(1X,e16.9))')(dbgtbl(i,j),j=1,num_input+num_output+num_extra)
       ENDDO
       CLOSE(17)

       DEALLOCATE( dbgtbl )

    ENDIF

!========================

    DEALLOCATE( cs )
    DEALLOCATE( rg )
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

