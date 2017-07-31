!------------------------------------------------------------------------
! globalmap.f90
!
! PURPOSE: 
!------------------------------------------------------------------------

program globalmap

  use mpi
  use globalmap_globals
  use cgyro_globals
  use cgyro_io
  use GEO_interface

  implicit none

  integer :: supported
  integer :: ir,it
  real :: rho_star,dr,ri,lr

  !-----------------------------------------------------------------
  ! Initialize MPI_COMM_WORLD communicator
  call MPI_INIT_THREAD(MPI_THREAD_SINGLE,supported,i_err)
  !-----------------------------------------------------------------

  ! Path is cwd:
  path= './'

  !-----------------------------------------------------------------
  ! Set the world MPI communicator
  !
  CGYRO_COMM_WORLD = MPI_COMM_WORLD
  !
  ! Query rank and size
  !
  call MPI_COMM_RANK(CGYRO_COMM_WORLD,i_proc,i_err)
  call MPI_COMM_SIZE(CGYRO_COMM_WORLD,n_proc,i_err)
  !-----------------------------------------------------------------

  call cgyro_read_input
  open(unit=1,file='input.globalmap',status='old')
  read(1,*) rho_star
  close(1)

  test_flag   = 1
  silent_flag = 1
  call cgyro_kernel

  lr = rho_star/rho*length
 
  d_theta = (2*pi/n_theta)
  do it=1,n_theta
     theta(it) = -pi+(it-1)*d_theta
  enddo

  print '(a,1pe12.5)','INFO: (globalmap) L_x   = ',lr
  print '(a,1pe12.5)','INFO: (globalmap) rho_* = ',rho_star
  print *
  print '(a)',' r/a          q            Rmaj'

  open(unit=1,file='out.globalmap',status='replace')
  do ir=1,n_radial

     ri = rmin + (ir-(n_radial+1.0)/2)*lr/(n_radial-1)
     dr = ri-rmin

     GEO_rmin_in  = ri
     GEO_rmaj_in  = rmaj+shift*dr
     GEO_drmaj_in = shift
     GEO_zmag_in  = zmag+dzmag*dr
     GEO_dzmag_in = dzmag
     GEO_q_in     = q*(1+s*dr/rmin)
     GEO_s_in     = s
     GEO_kappa_in = kappa*(1+s_kappa*dr/rmin)
     GEO_s_kappa_in = s_kappa
     GEO_delta_in   = delta*(1+s_delta*dr/rmin)
     GEO_s_delta_in = s_delta
     GEO_zeta_in    = zeta*(1+s_zeta*dr/rmin)
     GEO_s_zeta_in  = s_zeta
     GEO_beta_star_in = beta_star(0)

     print '(3(1pe12.5,1x))',ri,GEO_q_in,GEO_rmaj_in

     call GEO_do()  

     do it=1,n_theta

        call GEO_interp(theta(it))     

        write(1,10) GEO_bigr
        write(1,10) GEO_nu

     enddo
  enddo
  close(1)
  print '(a,1pe12.5)','INFO: (globalmap) Wrote out.globalmap'
  
  call MPI_FINALIZE(i_err)

10 format(1pe12.5)
  
end program globalmap
