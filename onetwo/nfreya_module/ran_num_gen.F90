    MODULE ran_num_gen
       USE nrtype,                                            ONLY:  SP,DP,I4B

       USE common_constants,                                  ONLY  : izero,zeroc

       CONTAINS 
         SUBROUTINE set_random_number_seed_orig
           !-----------------------------------------------------------------------------------
           ! -- Allow expereimentationwith seed, random number generators
           !-----------------------------------------------------------------------------------

!#if defined (USEMPI)
           USE MPI_data,                                            ONLY : myid
!#endif


           IMPLICIT  NONE
           REAL(DP) rn
           REAL(SP) ranseed
           ranseed = 7**7 + myid            ! different sequence on each cpu
           ranseed =  myid
           rn  = random12 (ranseed)         ! ranseed is saved in random12

         END SUBROUTINE set_random_number_seed_orig


         SUBROUTINE set_random_number_seed
!-------------------------------------------------------------------------
! -- use F90 intrinsics to set seed for each processor
! -- want a different seed for each process:
!--------------------------------------------------------------HSJ--------


  USE nrtype,                                        ONLY : SP,DP,I4B

  USE MPI_data,                                      ONLY : myid

  USE neutral_beams,                                 ONLY : randomize_seed

  USE io_gcnmp,                                      ONLY : ncrt

      IMPLICIT NONE
 
      INTEGER(I4B) i,n,values(8),start_val
      INTEGER(I4B), ALLOCATABLE,DIMENSION(:) :: seed
      LOGICAL randomize


      CALL RANDOM_SEED( SIZE = n )  ! need n integers to define seed
      ALLOCATE(seed(n))
      IF(randomize_seed)THEN ! change it each time we run the code
           CALL DATE_AND_TIME ( VALUES = values )
           start_val = myid*1000*values(7) + values(8) ! number based on msec time
           SEED = (/ ( i * start_val, i = 1, n )/)
           !WRITE(myid+500,'(" seed  init=",10(i5,x))')seed
           !CALL RANDOM_SEED( PUT = (/ ( i * start_val, i = 1, n ) /) )   !
           CALL RANDOM_SEED( PUT = seed)
      ELSE   ! keep the seed the same each time we run the code
            start_val = myid + 100
            SEED = (/ ( i * start_val, i = 1, n )/)
            CALL RANDOM_SEED( PUT = seed )   !
      ENDIF
      CALL  RANDOM_SEED(GET = seed)
      !write to fort.myid+500:
      !WRITE(myid+500,'("process =",i5," n ints req =",i5)')myid,n
      !WRITE(myid+500,'(" seed =",10(i10,x))')seed


      RETURN

      END SUBROUTINE set_random_number_seed


      FUNCTION random12_old (seed)
!-------------------------------------------------------------------------
! -- basic uniform (0,1) pseudo-random number generator (portable) ----
!-------------------------------------------------------------------------
      IMPLICIT NONE
      REAL(SP) saved_seed,seed
      SAVE         saved_seed
      REAL(DP) random12_old,  big, big_minus_1
      DATA   big /2147483648.0_DP/, big_minus_1 /2147483647.0_DP/

      IF (seed .NE. 0.0 )saved_seed = seed ! seed is alwayas zero except on first call
                                           ! in sub set_random_number_seed
                                           ! so saved_seed is alwyas changed by next statement
      saved_seed = MOD (16807.0_DP*saved_seed, big_minus_1)
      random12_old     = saved_seed / big

      RETURN

      END FUNCTION random12_old

      FUNCTION random12 (seed)
!-----------------------------------------------------------------------
! -- This routine  takes place of random12_old and should be used
! -- with set_random_number_seed. The input value of seed is retained
! -- for compatibility but is not used. The seed issue was taken care of in
! -- set_random_number_seed.
!-------------------------------------------------------------------------
      IMPLICIT NONE
      REAL(SP) ran_num,seed
      REAL(DP) random12
      
           CALL RANDOM_NUMBER(ran_num)

           random12 = ran_num ! converts to double

      RETURN

      END FUNCTION random12


      FUNCTION ranorm_old ( )
!--------------------------------------------------------------------------
!     generate one normal (0,1) pseudo-random number using a method
!     based on central limit theorem and power residue method.
!     output becomes truly normal as k (below) goes to infinity.
!     general formula for the deviate is:
!         y  =  (sum of x(i),i=1,k) -k/2.0) / SQRT (k/12.0)
!     where x(i) are are uniformly distributed on 0,1.
!     method is borrowed from IBM; they use k = 12
!--------------------------------------------------------------------------

      IMPLICIT  NONE
      INTEGER(I4B),PARAMETER :: k = 12
      REAL(DP),    PARAMETER :: rtkd12 = 1.0_DP
      REAL(DP) ranorm_old,a,y
      REAL(SP) seed0
      INTEGER(I4B) i


 
      DATA       seed0 /0.0_sP/

      a      = zeroc

! --- note:  this loop is vector hazard

      DO i=1,k
        y = random12 (seed0)
        a = a + y
      END DO

      ranorm_old = (a - 0.5 * k) / rtkd12

      RETURN

      END FUNCTION ranorm_old


      FUNCTION ranorm ( ) RESULT( ranormd )
!--------------------------------------------------------------------------
! -- generate N(0,1) distributed variable
!--------------------------------------------------------------------------

!   Generate a random normal deviate using the polar method.
!   Reference: Marsaglia,G. & Bray,T.A. 'A convenient method for generating
!              normal variables', Siam Rev., vol.6, 260-264, 1964.

        IMPLICIT NONE

        REAL(DP) ranormd



        ! Local variables
        REAL(SP)            :: fn_val
        REAL(SP)            :: u, sum
        REAL(SP), SAVE      :: v, sln
        LOGICAL, SAVE       :: second = .FALSE.
        REAL(SP), PARAMETER :: one = 1.0, vsmall = TINY( one )

        IF (second) THEN
           ! If second, use the second random number generated on last call

           second = .false.
           fn_val = v*sln

        ELSE
           ! First call; generate a pair of random normals

           second = .true.
           DO
              CALL RANDOM_NUMBER( u )
              CALL RANDOM_NUMBER( v )
              u = SCALE( u, 1 ) - one
              v = SCALE( v, 1 ) - one
              sum = u*u + v*v + vsmall ! vsmall added to prevent LOG(zero) / zero
              IF(sum < one) EXIT
           END DO
           sln = SQRT(- SCALE( LOG(sum), 1 ) / sum)
           fn_val = u*sln
        END IF
        ranormd = fn_val ! converst to DP

        RETURN

      END FUNCTION ranorm


    END MODULE ran_num_gen
