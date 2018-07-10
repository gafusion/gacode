subroutine EXPRO_skip_header(io)

  implicit none

  integer, intent(in) :: io
  character (len=1) :: cdummy

  do
     read(io,'(a)') cdummy     
     if (cdummy /= '#') exit
  enddo
  backspace io 

end subroutine EXPRO_skip_header
