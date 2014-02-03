module brainfuse_mod

  integer include_brainfuse
  
contains
  
  SUBROUTINE GRADIENT(nn,y,yy)
    IMPLICIT NONE
    INTEGER*4 nn,j
    REAL*8 y(nn),yy(nn)

    yy(1)=y(2)-y(1)
    yy(2:nn-1)=(y(3:nn)-y(1:nn-2))/2.
    yy(nn)=y(nn)-y(nn-1)

  END SUBROUTINE GRADIENT

  subroutine brainfuse(nn,bt,ip,r,rmaj,kappa,ne,ni,te,ti,q,vol,wt)
    IMPLICIT NONE
    INTEGER*4 debug, nn, j, i, dummy
    PARAMETER (debug=1)
    REAL*8  bt, ip, &
         r(nn), rmaj(nn), kappa(nn), &
         ne(nn), ni(nn), te(nn), ti(nn), &
         q(nn), vol(nn), wt(nn)
    REAL*8, DIMENSION (:), ALLOCATABLE :: cs, dte, dti, dne, dni, dr, dq, dkappa, dvol
    REAL*8 a, mD, eV
    PARAMETER (mD= 3.3475E-27)
    PARAMETER (eV= 1.60217646E-19)
!        _mE = 9.10938291E-31
!        _c = 3E8
!        _e0 = 8.8541878176E-12
!        _Kb = 1.3806505E-23
    CHARACTER input_names(20)*255
    CHARACTER output_names(20)*255

    REAL*4, DIMENSION (:,:), ALLOCATABLE :: input, output
    INTEGER num_data, num_input, num_output

    IF (debug) WRITE(*,*)'Running brainfuse!'

!===============================
    input_names(1)='rmin_loc'
    input_names(2)='rmaj_loc'
    input_names(3)='kappa'
    input_names(4)='taus'
    input_names(5)='aus'
    input_names(6)='q'

    output_names(1)='Qe'
    output_names(2)='Qi'
!===============================

    a=r(nn)

    ALLOCATE( cs(nn) )
    ALLOCATE( dr(nn) )
    ALLOCATE( dte(nn) )
    ALLOCATE( dti(nn) )
    ALLOCATE( dne(nn) )
    ALLOCATE( dni(nn) )
    ALLOCATE( dq(nn) )
    ALLOCATE( dkappa(nn) )
    ALLOCATE( dvol(nn) )

    cs = SQRT(1E3*eV*te/mD)
    call GRADIENT(nn,r,dr)
    call GRADIENT(nn,te,dte)
    call GRADIENT(nn,ti,dti)
    call GRADIENT(nn,ne,dne)
    call GRADIENT(nn,ni,dni)
    call GRADIENT(nn,q,dq)
    call GRADIENT(nn,kappa,dkappa)
    call GRADIENT(nn,vol,dvol)

!===============================

    num_input=6
    num_output=2
    ALLOCATE( input(nn,num_input)   )
    ALLOCATE( output(nn,num_output) )
    DO j=1,num_input
       dummy=1
       SELECT CASE (trim(input_names(j)))
       CASE('r')
          input(:,j) = r                       !minor radius
       CASE('rmaj')
          input(:,j) = rmaj                    !major radius
       CASE('te')
          input(:,j) = te                      !main electron temperature
       CASE('ti')
          input(:,j) = ti                      !main ion temperature
       CASE('ne')
          input(:,j) = ne                      !electron density
       CASE('ni')
          input(:,j) = ni                      !main ion density
       CASE('cs')
          input(:,j) = cs                      !electron sound speed
       CASE('vol')
          input(:,j) = vol                     !volume
       CASE('dvol')
          input(:,j) = dvol                    !differential volume
       CASE('rmin_loc')
          input(:,j) = r/a                     !normalized minor radius
       CASE('rmaj_loc')
          input(:,j) = rmaj/a                  !normalized major radius
       CASE('kappa')
          input(:,j) = kappa                   !elongation
       CASE('dkappa')
          input(:,j) = dkappa                  !elongation shear
       CASE('taus')
          input(:,j) = ti/te                   !temperature ratio
       CASE('aus')
          input(:,j) = ni/ne                   !density ratio
       CASE('q')
          input(:,j) = q                       !safety factor
       CASE('dq')
          input(:,j) = dq                      !safety factor shear
       CASE('vpar')
          input(:,j) = -ABS(Ip)/Ip*rmaj*wt/cs  !vpar
       CASE('lte')
          input(:,j) = -a*dte/te               !electron temperature scale length
       CASE('lti')
          input(:,j) = -a*dti/ti               !ion temperature scale length
       CASE('lne')
          input(:,j) = -a*dne/ne               !electron density scale length
       CASE('lni')
          input(:,j) = -a*dni/ni               !ion density scale length
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
!===============================

    CALL run_net_on_data(nn, num_input, num_output, input, output, debug)

!===============================

    IF (debug) THEN
!       WRITE(6,10) nn, num_input, num_output
!10     FORMAT('nn= ',i4,' num_input= ',i4,' num_output= ',i4)
!       write(*,*) 'r=',r
!       write(*,*) 'Bt=',bt
!       write(*,*) 'Ip=',ip
!       write(*,*) 'R=',rmaj
!       write(*,*) 'kappa=',kappa
!       write(*,*) 'ne=',ne
!       write(*,*) 'ni=',ni
!       write(*,*) 'te=',te
!       write(*,*) 'ti=',ti
!       write(*,*) 'q=',q
!       write(*,*) 'vol=',vol
!       write(*,*) 'wt=',wt
!       WRITE(6,*)'dr=',dr
!       WRITE(6,*)'Qe=',output(:,1)
!       WRITE(6,*)'Qi=',output(:,2)

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
!===============================

    DEALLOCATE( cs )
    DEALLOCATE( dr )
    DEALLOCATE( dte )
    DEALLOCATE( dti )
    DEALLOCATE( dne )
    DEALLOCATE( dni )
    DEALLOCATE( dq )
    DEALLOCATE( dkappa )
    DEALLOCATE( dvol )
    DEALLOCATE( input )
    DEALLOCATE( output )
    
  end subroutine brainfuse
  
end module

