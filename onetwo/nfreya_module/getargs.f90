

 SUBROUTINE my_getarg(carg,numarg,unit)

!---------------------------------------------------------------------
!get command line arguments 
!carg(0) = name of program
!carg(1) = first argument
!....
!....
!carg(numargs) = last argument
!----------------------------------------------------------------HSJ---



  USE nrtype,  ONLY : DP,I4B,I2B
  USE io_gcnmp, ONLY : ncrt,nlog
  USE error_handler

  IMPLICIT NONE
  INTEGER(I2B),INTENT(IN),OPTIONAL :: unit
  CHARACTER (len=256) arg
  CHARACTER ( len = *),DIMENSION(0:),INTENT(OUT) ::  carg
  INTEGER(I4B) i, iargc,ierr
  INTEGER(I2B)   iounit
  INTEGER(I2B),INTENT(OUT) :: numarg
  INTEGER(I4B) rl,rlc
  ierr =0
  IF(PRESENT( unit))THEN
     iounit = unit
  ELSE
     iounit =ncrt
  ENDIF
  CALL timestamp ( iounit)

  numarg = iargc ( )

  IF (numarg == 0 )THEN
     lerrno = 200 + iomaxerr
     CALL terminate(lerrno,nlog)
  ENDIF



  numarg = MIN(numarg,3)  ! if this is an mpirun numarg > 3 but
                          ! we only need the i = 0,1,2,3 slots of carg(i)
  DO i = 0, numarg
    CALL getarg ( i, arg )
    carg(i)(:) =' '
    rlc = LEN(carg(i))
    rl  = LEN_TRIM(arg)
    IF(rl .LE. rlc)THEN
         carg(i)(1:rl) = arg(1:rl)
    ELSE
       lerrno = 2
       CALL terminate(lerrno,nlog)
    ENDIF
    WRITE ( *, '(2x,i3,2x,a)' ) i, carg(i)

  END DO

  CALL timestamp (iounit )

  RETURN

END SUBROUTINE my_getarg





SUBROUTINE timestamp (iounit )

!*******************************************************************************
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    May 31 2001   9:45:54.872 AM
!
!  Modified:
!
!    31 May 2001
!
!  Author:
!
!    John Burkardt (original, customized by HSJ)
!
!  Parameters:
!
!    None
!



  USE nrtype,  ONLY : DP,I4B,I2B
  USE nfreya_version, ONLY : p_nfreya_ver

  IMPLICIT NONE

  CHARACTER ( len = 8 ) ampm
  INTEGER(I4B) d
  CHARACTER ( len = 8 ) date
  INTEGER(I4B) h
  INTEGER(I2B) iounit
  INTEGER(I4B) m
  INTEGER(I4B) mm
  CHARACTER ( len = 9 ), PARAMETER, DIMENSION(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  INTEGER(I4B) n
  INTEGER(I4B) s
  CHARACTER ( len = 10 )  time
  INTEGER(I4B) values(8)
  INTEGER(I4B) y
  CHARACTER ( len = 5 ) zone

  CALL DATE_AND_TIME ( date, time, zone, values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  IF ( h < 12 ) THEN
    ampm = 'AM'
  ELSE IF ( h == 12 ) THEN
    IF ( n == 0 .AND. s == 0 ) THEN
      ampm = 'Noon'
    ELSE
      ampm = 'PM'
    END IF
  ELSE
    h = h - 12
    IF ( h < 12 ) THEN
      ampm = 'PM'
    ELSE IF ( h == 12 ) THEN
      IF ( n == 0 .AND. s == 0 ) THEN
        ampm = 'Midnight'
      ELSE
        ampm = 'AM'
      END IF
    END IF
  END IF

  WRITE ( iounit, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a,5x,a)' ) &
    TRIM ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, TRIM ( ampm ),TRIM(P_nfreya_ver)

  RETURN
END SUBROUTINE timestamp
