module parallel_lib

  implicit none

  ! lib

  integer, private :: nproc,iproc
  integer, private :: ni,nj
  integer, private :: ni_loc
  integer, private :: nj_loc
  integer, private :: lib_comm
  integer, private :: nsend

  complex, dimension(:,:,:), allocatable, private :: fsendf
  complex, dimension(:,:,:), allocatable, private :: fsendr
  real, dimension(:,:,:), allocatable, private :: fsendr_real

  ! slib

  integer, private :: nsproc,isproc
  integer, private :: nn
  integer, private :: nexch
  integer, private :: nkeep
  integer, private :: nsplit
  integer, private :: slib_comm

  !integer, private, parameter :: default_size = 1024*1024*32
  integer, private, parameter :: default_size = 1

contains

  !=========================================================

  !  parallel_lib_f -> f(ni_loc,nj) -> g(nj_loc,ni) 
  !  parallel_lib_r -> g(nj_loc,ni) -> f(ni_loc,nj)

  subroutine parallel_lib_init(ni_in,nj_in,ni_loc_out,nj_loc_out,comm)

    use mpi

    implicit none

    integer, intent(in) :: ni_in,nj_in
    integer, intent(in) :: comm
    integer, intent(inout) :: ni_loc_out,nj_loc_out
    integer, external :: parallel_dim
    integer :: ierr

    lib_comm = comm

    call MPI_COMM_RANK(lib_comm,iproc,ierr)
    call MPI_COMM_SIZE(lib_comm,nproc,ierr)

    ni = ni_in
    nj = nj_in

    ni_loc = parallel_dim(ni,nproc)
    nj_loc = parallel_dim(nj,nproc)

    ni_loc_out = ni_loc
    nj_loc_out = nj_loc

    nsend = ni*nj/nproc**2

    allocate(fsendf(nj_loc,ni_loc,nproc))
    allocate(fsendr(ni_loc,nj_loc,nproc))
    allocate(fsendr_real(ni_loc,nj_loc,nproc))

  end subroutine parallel_lib_init

  !=========================================================

  subroutine parallel_lib_f(f,ft)

    use mpi

    implicit none

    complex, intent(in), dimension(ni_loc,nj) :: f
    complex, intent(inout), dimension(nj_loc,ni) :: ft
    integer :: ierr,i_loc,i,j,k,i1,i2

    i1 = 1+iproc*ni_loc
    i2 = (1+iproc)*ni_loc

!$omp parallel do if (size(fsendf) >= default_size) default(none) &
!$omp& shared(nproc,i1,i2,nj_loc) &
!$omp& private(i,i_loc,j) &
!$omp& shared(f,fsendf)
    do k=1,nproc
       do i=i1,i2
          i_loc = i-i1+1 
          do j=1,nj_loc
             fsendf(j,i_loc,k) = f(i_loc,j+(k-1)*nj_loc) 
          enddo
       enddo
    enddo

    call MPI_ALLTOALL(fsendf, &
         nsend, &
         MPI_DOUBLE_COMPLEX, &
         ft, &
         nsend, &
         MPI_DOUBLE_COMPLEX, &
         lib_comm, &
         ierr)

  end subroutine parallel_lib_f

  ! inline version, set one row
  subroutine parallel_lib_f_i_set(i_loc,bv)

    implicit none

    integer, intent(in) :: i_loc
    complex, intent(in), dimension(nj) :: bv
    integer :: j,k

    do k=1,nproc
          do j=1,nj_loc
             fsendf(j,i_loc,k) = bv(j+(k-1)*nj_loc)
          enddo
    enddo

  end subroutine parallel_lib_f_i_set

  ! inline vesion, just execute
  subroutine parallel_lib_f_i_do(ft)

    use mpi

    implicit none

    complex, intent(inout), dimension(nj_loc,ni) :: ft
    integer :: ierr

    call MPI_ALLTOALL(fsendf, &
         nsend, &
         MPI_DOUBLE_COMPLEX, &
         ft, &
         nsend, &
         MPI_DOUBLE_COMPLEX, &
         lib_comm, &
         ierr)

  end subroutine parallel_lib_f_i_do

  !=========================================================

  subroutine parallel_lib_r(ft,f)

    use mpi

    implicit none

    complex, intent(in), dimension(nj_loc,ni) :: ft
    complex, intent(inout), dimension(ni_loc,nj) :: f
    integer :: ierr,j_loc,i,j,k,j1,j2

    j1 = 1+iproc*nj_loc
    j2 = (1+iproc)*nj_loc

!$omp parallel do if (size(fsendr) >= default_size) default(none) &
!$omp& shared(nproc,j1,j2,ni_loc) &
!$omp& private(j,j_loc,i) &
!$omp& shared(ft,fsendr)
    do k=1,nproc
       do j=j1,j2
          j_loc = j-j1+1 
          do i=1,ni_loc
             fsendr(i,j_loc,k) = ft(j_loc,i+(k-1)*ni_loc) 
          enddo
       enddo
    enddo

    call MPI_ALLTOALL(fsendr, &
         nsend, &
         MPI_DOUBLE_COMPLEX, &
         f, &
         nsend, &
         MPI_DOUBLE_COMPLEX, &
         lib_comm, &
         ierr)

  end subroutine parallel_lib_r

  subroutine parallel_lib_rtrans(fin,f)
! -----------------------------------------
! transpose version of parallel_lib_r(fin,f)
! -----------------------------------------
    use mpi

    implicit none

    complex, intent(in), dimension(:,:) :: fin
    complex, intent(inout), dimension(ni_loc,nj) :: f
    integer :: ierr,j_loc,i,j,k,j1,j2

    j1 = 1+iproc*nj_loc
    j2 = (1+iproc)*nj_loc

!$omp parallel do if (size(fsendr) >= default_size) default(none) &
!$omp& shared(nproc,j1,j2,ni_loc) &
!$omp& private(j,j_loc,i) &
!$omp& shared(fin,fsendr)
    do k=1,nproc
       do j=j1,j2
          j_loc = j-j1+1 
          do i=1,ni_loc
             fsendr(i,j_loc,k) = fin(i+(k-1)*ni_loc,j_loc) 
          enddo
       enddo
    enddo

    call MPI_ALLTOALL(fsendr, &
         nsend, &
         MPI_DOUBLE_COMPLEX, &
         f, &
         nsend, &
         MPI_DOUBLE_COMPLEX, &
         lib_comm, &
         ierr)

  end subroutine parallel_lib_rtrans

  !=========================================================

  subroutine parallel_lib_rtrans_real(fin,f)

    use mpi

    implicit none

    real, intent(in), dimension(:,:) :: fin
    real, intent(inout), dimension(ni_loc,nj) :: f
    integer :: ierr,j_loc,i,j,k,j1,j2

    j1 = 1+iproc*nj_loc
    j2 = (1+iproc)*nj_loc

!$omp parallel do if (size(fsendr_real) >= default_size) default(none) &
!$omp& shared(nproc,j1,j2,ni_loc) &
!$omp& private(j,j_loc,i) &
!$omp& shared(fin,fsendr_real)
    do k=1,nproc
       do j=j1,j2
          j_loc = j-j1+1
          do i=1,ni_loc
             fsendr_real(i,j_loc,k) = fin(i+(k-1)*ni_loc,j_loc) 
          enddo
       enddo
    enddo

    call MPI_ALLTOALL(fsendr_real, &
         nsend, &
         MPI_DOUBLE_PRECISION,&
         f, &
         nsend, &
         MPI_DOUBLE_PRECISION, &
         lib_comm, &
         ierr)

  end subroutine parallel_lib_rtrans_real

  !=========================================================

  !  parallel_slib_f -> f(ni_loc,nj) -> g(nj_loc,ni) 
  !  parallel_slib_r -> g(nj_loc,ni) -> f(ni_loc,nj)

  subroutine parallel_slib_init(nn_in,nexch_in,nkeep_in,nsplit_out,comm)

    use mpi

    !-------------------------------------------
    implicit none
    !
    integer, intent(in) :: nn_in
    integer, intent(in) :: nexch_in
    integer, intent(in) :: nkeep_in
    integer, intent(inout) :: nsplit_out
    integer, intent(in) :: comm
    integer :: ierr
    !-------------------------------------------

    slib_comm = comm

    !-------------------------------------------------
    ! Get rank and number of processors from slib_comm 
    ! (which is assumed to already exist).
    !
    call MPI_COMM_RANK(slib_comm,isproc,ierr)
    call MPI_COMM_SIZE(slib_comm,nsproc,ierr)
    !-----------------------------------------------

    nn  = nn_in
    nexch = nexch_in
    nkeep = nkeep_in

    nsplit = 1+(nexch-1)/nn

    nsplit_out = nsplit

  end subroutine parallel_slib_init

!=========================================================

  subroutine parallel_slib_f_nc(x,xt)

    use mpi

    !-------------------------------------------------------
    implicit none
    !
    complex, intent(in), dimension(nkeep,nsplit*nn) :: x
    complex, intent(inout), dimension(nkeep,nsplit,nn) :: xt
    !
    integer :: ierr
    !-------------------------------------------------------

    call MPI_ALLTOALL(x, &
         nkeep*nsplit, &
         MPI_DOUBLE_COMPLEX, &
         xt, &
         nkeep*nsplit, &
         MPI_DOUBLE_COMPLEX, &
         slib_comm, &
         ierr)

  end subroutine parallel_slib_f_nc

  subroutine parallel_slib_f(x_in,xt)

    use mpi

    !-------------------------------------------------------
    implicit none
    !
    complex, intent(in), dimension(nkeep,nexch) :: x_in
    complex, dimension(nkeep,nsplit*nn) :: x
    complex, intent(inout), dimension(nkeep,nsplit,nn) :: xt
    !
    integer :: j
    !-------------------------------------------------------

!$omp parallel do private(j)
    do j=1,nexch
       x(:,j) = x_in(:,j)
    enddo

    do j=nexch+1,nsplit*nn
       x(:,j) = (0.0,0.0)
    enddo

    call parallel_slib_f_nc(x,xt)

  end subroutine parallel_slib_f

 !=========================================================

  subroutine parallel_slib_r_nc (xt,x)

    use mpi

    !-------------------------------------------------------
    implicit none
    !
    complex, intent(in), dimension(nkeep,nsplit,nn) :: xt
    complex, intent(inout), dimension(nkeep,nsplit*nn) :: x
    !
    integer :: ierr
    !-------------------------------------------------------

    call MPI_ALLTOALL(xt, &
         nkeep*nsplit, &
         MPI_DOUBLE_COMPLEX, &
         x, &
         nkeep*nsplit, &
         MPI_DOUBLE_COMPLEX, &
         slib_comm, &
         ierr)

  end subroutine parallel_slib_r_nc

  subroutine parallel_slib_r(xt,x_out)

    use mpi

    !-------------------------------------------------------
    implicit none
    !
    complex, intent(in), dimension(nkeep,nsplit,nn) :: xt
    complex, intent(inout), dimension(nkeep,nexch) :: x_out
    complex, dimension(nkeep,nsplit*nn) :: x
    !
    integer :: j
    !-------------------------------------------------------

    call parallel_slib_r_nc(xt,x)

!$omp parallel do private(j)
    do j=1,nexch
       x_out(:,j) = x(:,j)
    enddo

  end subroutine parallel_slib_r

end module parallel_lib
