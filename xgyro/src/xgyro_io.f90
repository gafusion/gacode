module xgyro_io 

contains

  !-----------------------------------------------------------
  ! xgyro_info.f90
  !
  ! PURPOSE:
  !  Routine to write line to run file.
  !-----------------------------------------------------------

  subroutine xgyro_info(message)

    use xgyro_globals, only : xgyro_runfile_info, xgyro_i_proc
    use cgyro_globals, only : silent_flag, io, path

    implicit none

    character (len=*), intent(in) :: message

    if (silent_flag == 0 .and. xgyro_i_proc == 0) then
       open(unit=io,file=trim(path)//xgyro_runfile_info,status='old',position='append')
       write(io,'(a)') 'INFO: (XGYRO) '//message
       close(io)
    endif

  end subroutine xgyro_info

  !-----------------------------------------------------------
  ! xgyro_error.f90
  !
  ! PURPOSE:
  !  Routine to write line to error file.
  !-----------------------------------------------------------

  subroutine xgyro_error(message)

    use xgyro_globals, only : xgyro_runfile_info, xgyro_i_proc
    use cgyro_globals, only : error_status, error_message, &
         silent_flag, io, path

    implicit none

    character (len=*), intent(in) :: message

    error_status  = 1
    error_message = message

    if (silent_flag == 0 .and. xgyro_i_proc == 0) then
       open(unit=io,file=trim(path)//xgyro_runfile_info,status='old',position='append')
       write(io,'(a)') 'ERROR: (XGYRO) '//message
       close(io)
    endif

  end subroutine xgyro_error

end module xgyro_io
