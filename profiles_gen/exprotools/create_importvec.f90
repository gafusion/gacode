subroutine create_importvec(datafile,x,y,n)

  use create_globals

  implicit none

  character (len=*), intent(in) :: datafile
  integer, intent(in) :: n
  real, intent(inout), dimension(n) :: x
  real, intent(inout), dimension(n) :: y

  real :: a
  integer :: ierr
  integer :: n_data
  integer :: eof
  real, dimension(:,:), allocatable :: data


  open(unit=1,file=trim(datafile),status='old',iostat=ierr) 

  if (ierr == 0) then
     n_data = 0
     do 
        read(1,*,iostat=eof) a
        if (eof /= 0) exit
        n_data = n_data+1
     enddo
     rewind(1)
     allocate(data(2,n_data)) 
     read(1,*) data(:,:)
     close(1)
  else
     print '(a)','ERROR: (create) '//datafile//' does not exist'
     stop
  endif
  call cub_spline(data(1,:),data(2,:),n_data,x,y,n)

  print '(a)','INFO: (create) Imported data from '//datafile

end subroutine create_importvec
