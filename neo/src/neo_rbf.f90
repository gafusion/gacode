subroutine neo_rbf(x0,c0)

  use neo_globals

  implicit none

  ! Number of inputs (hypercube dimensions)
  integer, parameter :: n_in = 6
  ! NUmber of outputs (gradient coefficients)
  integer, parameter :: n_out = 6

  real, dimension(n_in), intent(inout) :: x0
  real, dimension(n_out), intent(out) :: c0
  real, dimension(n_in) :: xs0
  real, dimension(n_in) :: xmin, xmax
  real :: cscale_nu

  integer :: p,k
  integer :: stat

  integer :: ntot
  real, dimension(:,:), allocatable :: x
  real, dimension(:,:), allocatable :: b

  real, dimension(2,n_in) :: s

  character(len=255) :: data
  character(len=218) :: root

  character(len=8) :: rbf_type
  real :: rbf_eps

  call get_environment_variable('GACODE_ROOT',root)
  data = trim(root)//'/neo/tools/pneo/data/'//trim(rbf_dir)

  open(unit=1,file=trim(data)//'/out.pneo.scale',status='old',iostat=stat)
  if (stat /= 0) then
     call neo_error('ERROR: (neo_rbf) out.pneo.scale not available')
     return
  endif
  read(1,*) rbf_type
  read(1,*) rbf_eps
  read(1,*) ntot
  do k=1,n_in
     read(1,*) xmin(k), xmax(k)
  enddo
  do k=1,n_in
     read(1,*) s(:,k)  
  enddo
  close(1)

  ! Re-scale x0 if out-of-range of training data
  ! For now, do not re-scale for epsilon and ftrap
  cscale_nu = 1.0
  do k=2,n_in-1
     if(x0(k) > xmax(k))  then
        x0(k) = xmax(k)
        !print *, 'max rescaling ', k
        if(k == 3) then
           cscale_nu = 10.0**xmax(k)/10.0**x0(k)
        endif
     endif
     if(x0(k) < xmin(k))  then
        x0(k) = xmin(k)
        !print *, 'min rescaling ', k
     endif
  enddo

  allocate(x(n_in,ntot))
  allocate(b(ntot,n_out))

  open(unit=1,file=trim(data)//'/out.pneo.rbf',status='old',iostat=stat)
  if (stat /= 0) then
     call neo_error('ERROR: (neo_rbf) out.pneo.rbf not available')
     return
  endif
  read(1,*) x(:,:)
  read(1,*) b(:,:)
  close(1)

  ! Rescale data
  do k=1,n_in
     xs0(k) = s(1,k)*x0(k)+s(2,k)
  enddo

  select case (trim(rbf_type))

  case ('linear') 
     do k=1,n_out
        c0(k) = 0.0
        do p=1,ntot
           c0(k) = c0(k)+b(p,k)*sqrt( sum((xs0(:)-x(:,p))**2) )
        enddo
     enddo

  case ('cubic') 
     do k=1,n_out
        c0(k) = 0.0
        do p=1,ntot
           c0(k) = c0(k)+b(p,k)*sqrt( sum((xs0(:)-x(:,p))**2) )**3
        enddo
     enddo

  case ('gaussian')
     do k=1,n_out
        c0(k) = 0.0
        do p=1,ntot
           c0(k) = c0(k)+b(p,k)*exp(-sum((xs0(:)-x(:,p))**2)/rbf_eps**2)
        enddo
     enddo

  end select

  ! Rescaling
  c0(:) = c0(:) * cscale_nu

end subroutine neo_rbf
