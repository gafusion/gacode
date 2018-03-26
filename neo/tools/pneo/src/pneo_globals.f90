module pneo_globals
  
  character(len=80) :: path
  
  integer :: i_proc
  integer :: n_proc
  integer :: i_err

  integer :: ntot
  integer, dimension(:), allocatable :: ic1,ic2,ic3,ic4,ic5,ic6,ic7,ic8,ic9
  real, dimension(:,:), allocatable :: indata_loc
  real, dimension(:,:), allocatable :: indata
  real, dimension(:,:), allocatable :: ingeodata_loc
  real, dimension(:,:), allocatable :: ingeodata
  real, dimension(:,:), allocatable :: outdata_loc
  real, dimension(:,:), allocatable :: outdata
  real, dimension(:,:), allocatable :: outsaudata_loc
  real, dimension(:,:), allocatable :: outsaudata


end module pneo_globals
