module brainfuse_mod

  integer include_brainfuse
  
contains
  
  subroutine brainfuse (nn,a,bt,r,rmaj,kappa,ne,ni,te,ti,q,vol)
    implicit none
    integer*4 debug, nn, j
    parameter (debug=1)
    real*8   a, bt, &
         r(nn), rmaj(nn), kappa(nn), &
         ne(nn), ni(nn), te(nn), ti(nn), &
         q(nn), vol(nn)

    INTEGER num_data, num_input, num_output
    REAL*4, DIMENSION (:,:), ALLOCATABLE :: input, output

!===============================
!         densh(*),    densimp(*),  densfe(*), &
!         xzeff(*),   tekev(*),    tikev(*),    q(*),       btor(*), &
!         avezimp(*), amassimp(*), amasshyd(*), aimass(*),  wexbs(*), &
!         grdne(*),   grdni(*),    grdnh(*),    grdnz(*), &
!         grdte(*),   grdti(*),    grdq(*)

    write(*,*)'Running brainfuse!',nn
    IF (debug) THEN
!       write(*,*) 'r/a=',r/a
!       write(*,*) 'Bt=',bt
!       write(*,*) 'R/a=',rmaj/a
!       write(*,*) 'kappa=',kappa
!       write(*,*) 'ne=',ne
!       write(*,*) 'ni=',ni
!       write(*,*) 'te=',te
!       write(*,*) 'ti=',ti
!       write(*,*) 'q=',q
       write(*,*) 'vol=',vol
    ENDIF
!===============================

    num_input=6
    num_output=2
    
    ALLOCATE( input(nn,num_input)   )
    ALLOCATE( output(nn,num_output) )

    WRITE(6,10) nn, num_input, num_output
10  FORMAT('nn= ',i4,' num_input= ',i4,' num_output= ',i4)
    
    input(:,1) = r/a
    input(:,2) = rmaj/a
    input(:,3) = kappa
    input(:,4) = ti/te
    input(:,5) = ni/ne
    input(:,6) = q
    
    CALL run_net_on_data(nn, num_input, num_output, input, output)
    
    WRITE(6,*)'Qe=',output(:,1)
    WRITE(6,*)'Qi=',output(:,2)
    
    DEALLOCATE(input)
    DEALLOCATE(output)
    
  end subroutine brainfuse
  
end module

