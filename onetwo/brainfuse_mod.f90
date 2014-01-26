module brainfuse_mod

  integer include_brainfuse
  
contains
  
  subroutine brainfuse (nn,a,bt,r,rmaj,kappa,ene,te,ti,q)
    implicit none
    integer*4 debug, nn, j
    parameter (debug=1)
    real*8   a, bt, &
         r(nn), rmaj(nn), kappa(nn), &
         ene(nn), te(nn), ti(nn), &
         q(nn)

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
       write(*,*) 'r/a=',r/a
       write(*,*) 'Bt=',bt
       write(*,*) 'R/a=',rmaj/a
       write(*,*) 'kappa=',kappa
       write(*,*) 'ne=',ene
       write(*,*) 'te=',te
       write(*,*) 'ti=',ti
       write(*,*) 'q=',q
    ENDIF
!===============================

    num_data=4
    num_input=2
    num_output=1
    
    ALLOCATE( input(num_data,num_input)   )
    ALLOCATE( output(num_data,num_output) )
    
    WRITE(6,10) num_data, num_input, num_output
10  FORMAT('num_data= ',i2,' num_input= ',i2,' num_output= ',i2)
    
    input(1,1) =-1.
    input(1,2) =-1.
    output(1,1)=-1.
    
    input(2,1) =-1.
    input(2,2) =1.
    output(2,1)=1.
    
    input(3,1) =1.
    input(3,2) =-1.
    output(3,1)=1.
    
    input(4,1) =1.
    input(4,2) =1.
    output(4,1)=-1.
    
    CALL run_net_on_data(num_data, num_input, num_output, input, output)
    
    WRITE(6,*)output
    
    DEALLOCATE(input)
    DEALLOCATE(output)
    
  end subroutine brainfuse
  
end module

