program pneo_rbf

  implicit none

  ! Number of inputs (hypercube dimensions)
  integer, parameter :: n_in = 6
  ! NUmber of outputs (gradient coefficients)
  integer, parameter :: n_out = 6

  integer :: p,q,k
  integer, dimension(9) :: n

  integer :: ntot
  real, dimension(:,:), allocatable :: indata
  real, dimension(:,:), allocatable :: ingeodata
  real, dimension(:,:), allocatable :: outdata
  real, dimension(:,:), allocatable :: x
  real, dimension(16) :: dum
  real, dimension(6) :: xmin,xmax
  integer :: stat

  real, dimension(:,:), allocatable :: a
  real, dimension(:,:), allocatable :: b
  integer, dimension(:), allocatable :: ipiv
  integer :: info

  character(len=8), parameter :: rbf_type='cubic'
  !character(len=8), parameter :: rbf_type='gaussian'
  !character(len=8), parameter :: rbf_type='linear'
  real, parameter :: rbf_eps=1.0

  open(unit=1,file='input.pneo',status='old',iostat=stat)
  if (stat /= 0) then
     print '(a)','ERROR: (pneo_rbf) input.pneo not available'
     stop
  endif
  do k=1,9
     read(1,*) n(k)
     read(1,*) dum(1:n(k))
  enddo
  close(1)

  ntot = n(1)*n(2)*n(3)*n(4)*n(5)*n(6)*n(7)*n(8)*n(9)
  allocate(indata(9,ntot))
  allocate(ingeodata(12,ntot))
  allocate(outdata(6,ntot))
  allocate(x(n_in,ntot))

  open(unit=1,file='out.pneo.indata',status='old',iostat=stat)
  if (stat /= 0) then
     print '(a)','ERROR: (neo_rbf) out.pneo.indata not available'
     stop
  endif
  open(unit=2,file='out.pneo.ingeo',status='old',iostat=stat)
  if (stat /= 0) then
     print '(a)','ERROR: (neo_rbf) out.pneo.ingeo not available'
     stop
  endif
  open(unit=3,file='out.pneo.c',status='old',iostat=stat)
  if (stat /= 0) then
     print '(a)','ERROR: (neo_rbf) out.pneo.c not available'
     stop
  endif

  do p=1,ntot

     do k=1,9
        read(1,*) indata(k,p)
     enddo
     do k=1,6
        read(3,*) outdata(k,p)
     enddo
     do k=1,12
        read(2,*) ingeodata(k,p)
     enddo

  enddo

  close(1)
  close(2)
  close(3)

  x(1:5,:) = indata(1:5,:)
  x(6,:)   = ingeodata(2,:)

  ! Magic geo is eps
  n(6) = n(1)

  do k=1,n_in
     xmin(k) = minval(x(k,:))
     xmax(k) = maxval(x(k,:))
  enddo

  ! Rescale data
  do k=1,n_in
     if(xmax(k)-xmin(k) > 0.0) then
        x(k,:) = (x(k,:)-xmin(k))/(xmax(k)-xmin(k))*n(k)
     endif
  enddo

  allocate(ipiv(ntot))
  allocate(b(ntot,n_out))
  allocate(a(ntot,ntot))

  select case (trim(rbf_type))

  case ('linear') 

     do p=1,ntot
        do q=1,ntot
           a(p,q) = sqrt( sum((x(:,q)-x(:,p))**2) )
        enddo
     enddo

  case ('cubic') 

     do p=1,ntot
        do q=1,ntot
           a(p,q) = sqrt( sum((x(:,q)-x(:,p))**2) )**3
        enddo
     enddo

  case ('gaussian')

     do p=1,ntot
        do q=1,ntot
           a(p,q) = exp(-sum((x(:,q)-x(:,p))**2)/rbf_eps**2 )
        enddo
     enddo

  end select

  b(:,:) = transpose(outdata(:,:))
  call DGESV(ntot,n_out,a,ntot,ipiv,b,ntot,info) 

  open(unit=1,file='out.pneo.rbf',status='replace')
  write(1,10) x(:,:)
  write(1,10) b(:,:)
  close(1)

  open(unit=1,file='out.pneo.scale',status='replace')
  ! Write (a,b) such that xscale = a*x0+b
  write(1,*) rbf_type
  write(1,*) rbf_eps
  write(1,*) ntot
  do k=1,n_in
     write(1,20) xmin(k), xmax(k)
  enddo
  do k=1,n_in
     if(xmax(k)-xmin(k) > 0.0) then
        write(1,20) n(k)/(xmax(k)-xmin(k)),-n(k)*xmin(k)/(xmax(k)-xmin(k))
     else
        write(1,20) 1.0,0.0
     endif
  enddo
  close(1)

10 format(1pe17.10)
20 format(2(1pe17.10,1x))

end program pneo_rbf
