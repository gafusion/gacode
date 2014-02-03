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
    INTEGER*4 debug, nn, j
    PARAMETER (debug=1)
    REAL*8  bt, ip, &
         r(nn), rmaj(nn), kappa(nn), &
         ne(nn), ni(nn), te(nn), ti(nn), &
         q(nn), vol(nn), wt(nn)
    REAL*8, DIMENSION (:), ALLOCATABLE :: cs, dte, dti, dne, dni, dr
    REAL*8 a, mD, eV
    PARAMETER (mD= 3.3475E-27)
    PARAMETER (eV= 1.60217646E-19)
!        _mE = 9.10938291E-31
!        _c = 3E8
!        _e0 = 8.8541878176E-12
!        _Kb = 1.3806505E-23

    REAL*4, DIMENSION (:,:), ALLOCATABLE :: input, output
    INTEGER num_data, num_input, num_output

!===============================

    a=r(nn)

    ALLOCATE( cs(nn) )
    ALLOCATE( dr(nn) )
    ALLOCATE( dte(nn) )
    ALLOCATE( dti(nn) )
    ALLOCATE( dne(nn) )
    ALLOCATE( dni(nn) )

    cs = SQRT(1E3*eV*te/mD)
    call GRADIENT(nn,r,dr)
    call GRADIENT(nn,te,dte)
    call GRADIENT(nn,ti,dti)
    call GRADIENT(nn,ne,dne)
    call GRADIENT(nn,ni,dni)

!===============================

    num_input=6
    num_output=2
    ALLOCATE( input(nn,num_input)   )
    ALLOCATE( output(nn,num_output) )
    input(:,1) = r/a                     !normalized minor radius
    input(:,2) = rmaj/a                  !normalized major radius
    input(:,3) = kappa                   !elongation
    input(:,4) = ti/te                   !temperature ratio
    input(:,5) = ni/ne                   !density ratio
    input(:,6) = q                       !safety factor
    !input(:,7) = -sign(Ip)*rmaj*wt/cs    !vpar
    !input(:,8) = -a*dte/te               !electron temperature scale length
    !input(:,9) = -a*dti/ti               !ion temperature scale length
    !input(:,10) = -a*dne/ne               !electron density scale length
    !input(:,11) = -a*dni/ni               !ion density scale length

!===============================

    CALL run_net_on_data(nn, num_input, num_output, input, output, debug)

!===============================

    IF (debug) THEN
       write(*,*)'Running brainfuse!'
       WRITE(6,10) nn, num_input, num_output
10     FORMAT('nn= ',i4,' num_input= ',i4,' num_output= ',i4)
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
       WRITE(6,*)'dr=',dr
!       WRITE(6,*)'Qe=',output(:,1)
!       WRITE(6,*)'Qi=',output(:,2)
    ENDIF

!===============================

    DEALLOCATE( cs )
    DEALLOCATE( dr )
    DEALLOCATE( dte )
    DEALLOCATE( dti )
    DEALLOCATE( dne )
    DEALLOCATE( dni )
    DEALLOCATE( input )
    DEALLOCATE( output )
    
  end subroutine brainfuse
  
end module

