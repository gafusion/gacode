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

     ! First, assume not a restart
     restart_flag = 0

     ! Check for tag file
     open(unit=io,&
          file=trim(path)//runfile_restart_tag,&
          status='old',iostat=i_err)
     close(io)

     if (i_err == 0) then
        ! Tag file exists: do restart
        restart_flag = 1
     else
        ! Check for restart file with no tagfile
        open(unit=io,&
             file=trim(path)//runfile_restart,&
             status='old',iostat=i_err)
        close(io)
        if (i_err == 0) then
           ! Use restart data as initial h
           restart_flag = 2
        endif
     endif

  endif

  call MPI_BCAST(restart_flag,1,MPI_INTEGER,0,CGYRO_COMM_WORLD,i_err)
  !---------------------------------------------------------------------------

  select case (restart_flag)

  case (1)

     call cgyro_info('Restart data found.')
     call cgyro_read_restart
     if (error_status /=0 ) return
     gtime = 0.0

     ! Rescale to prevent overflow
     if (nonlinear_flag == 0) then
        h_x = h_x/sum(abs(h_x))
     endif
     
  case (2)

     call cgyro_info('Initializing with restart data.')
     call cgyro_read_restart
     if (error_status /=0 ) return
     i_current = 0
     t_current = 0.0
     gtime = 0.0

  case (0)

     i_current = 0
     t_current = 0.0
     gtime = 0.0

     !-------------------------------------------------------------------------
     ! Generate analytic initial conditions
     !-------------------------------------------------------------------------

     h_x(:,:) = (0.0,0.0)

     if (zf_test_mode == 1) then

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

     else if (zf_test_mode >= 2) then

        ! 2. ELECTROMAGNETIC ZONAL-FLOW TEST

        call cgyro_zftest_em

     else if (n_toroidal == 1 .and. my_toroidal > 0) then

        ! 3. LINEAR n>0 SIMULATION

        do iv=nv1,nv2

           iv_loc = iv-nv1+1
           is = is_v(iv)

           if (is == 1) then
              do ic=1,nc
                 ir = ir_c(ic) 
                 it = it_c(ic)
                 ang = theta(it)+2*pi*px(ir)
                 if (amp >  0.0) then
                    h_x(ic,iv_loc) = rho/(1.0+ang**4)
                 else
                    h_x(ic,iv_loc) = rho*ang/(1.0+ang**4)
                 endif
              enddo
           endif

        enddo

     else

        ! 4. GENERAL CASE

        do ic=1,nc

           ir = ir_c(ic) 
           it = it_c(ic)

           if (my_toroidal == 0) then

              ! Zonal-flow initial condition

              arg = abs(px(ir))/real(n_radial)
              h_x(ic,:) = amp0*rho*exp(-arg)
              if (ir == 1 .or. px(ir) == 0) then
                 h_x(ic,:) = 0.0
              endif

           else 

              ! Finite-n initial condition

              if (amp > 0.0) then
                 h_x(ic,:) = amp*rho
              else
                 h_x(ic,:) = amp*rho/my_toroidal**2
              endif

           endif

        enddo

     endif
  end select

  call cgyro_field_c_cpu

  ! Initialize time-history of fields (-3,-2,-1) to initial field.
  field_old  = field
  field_old2 = field
  field_old3 = field

end subroutine cgyro_init_h
