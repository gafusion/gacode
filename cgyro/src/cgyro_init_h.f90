subroutine cgyro_init_h

  use mpi
  use cgyro_globals
  use cgyro_io

  implicit none

  integer :: ir,it,is,ie,ix
  real :: arg, ang

  !---------------------------------------------------------------------------
  ! Check to see if we have restart data available
  !
  if (i_proc == 0) then

     open(unit=io,&
          file=trim(path)//runfile_restart_tag,&
          status='old',iostat=i_err)
     close(io)

     if (i_err == 0) then
        restart_flag = 1
     else
        restart_flag = 0
     endif

  endif

  call MPI_BCAST(restart_flag,1,MPI_DOUBLE_PRECISION,0,CGYRO_COMM_WORLD,i_err)
  !---------------------------------------------------------------------------

  if (restart_flag == 1) then

     call cgyro_info('Restart data found.')
     call cgyro_read_restart
     gtime = 0.0

  else

     i_current = 0
     t_current = 0.0
     gtime = 0.0

     !-------------------------------------------------------------------------
     ! Generate analytic initial conditions
     !-------------------------------------------------------------------------

     h_x(:,:) = (0.0,0.0)

     if (zf_test_flag == 1) then

        ! 1. ZONAL-FLOW TEST

        do iv=nv1,nv2

           iv_loc = iv-nv1+1
           is = is_v(iv)
           ix = ix_v(iv)
           ie = ie_v(iv)

           do ic=1,nc

              ir = ir_c(ic) 
              it = it_c(ic)

              if (is == 1 .and. px(ir) /= 0) then
                 arg = k_perp(ic)*rho*vth(is)*mass(is)/(z(is)*bmag(it)) &
                      *sqrt(2.0*energy(ie))*sqrt(1.0-xi(ix)**2)
                 h_x(ic,iv_loc) = 1e-6*bessel_j0(abs(arg))

                 ! J0 here for the ions is equivalent to having
                 ! the electrons deviate in density.
                 ! Alternatively this is the result of instantaneous
                 ! gyroaveraging after the deposition of particles in
                 ! a certain k_radial mode.

              endif
           enddo
        enddo

     else if (n_toroidal == 1 .and. n > 0) then

        ! 2. LINEAR n>0 SIMULATION

        do iv=nv1,nv2

           iv_loc = iv-nv1+1
           is = is_v(iv)

           if (is == 1) then
              do ic=1,nc
                 ir = ir_c(ic) 
                 it = it_c(ic)
                 ang = theta(it)+2*pi*px(ir)
                 h_x(ic,iv_loc) = rho/(1.0+ang**4) 
              enddo
           endif

        enddo

     else

        ! 3. GENERAL CASE

        do ic=1,nc

           ir = ir_c(ic) 
           it = it_c(ic)

           if (n == 0) then

              ! Zonal-flow initial condition

              arg = abs(px(ir))/real(n_radial)
              h_x(ic,:) = amp0*rho*exp(-arg)

           else 

              ! Finite-n initial condition

              if (amp > 0.0) then
                 h_x(ic,:) = amp*rho
              else
                 h_x(ic,:) = amp*rho/n**2
              endif

           endif

        enddo

     endif
  endif

  call cgyro_filter
  call cgyro_field_c

  ! Initialize time-history of fields (-3,-2,-1) to initial field.
  field_old  = field
  field_old2 = field
  field_old3 = field

end subroutine cgyro_init_h
