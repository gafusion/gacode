module ptglf_globals
  
  character(len=80) :: path
  
  integer :: i_proc
  integer :: n_proc
  integer :: i_err

  integer :: ntot
  integer :: loc_comm
  integer :: color,splitkey
  
  real, dimension(:,:), allocatable :: indata_loc
  real, dimension(:,:), allocatable :: indata
  real, dimension(:,:), allocatable :: outdata_loc
  real, dimension(:,:), allocatable :: outdata
  
end module ptglf_globals
