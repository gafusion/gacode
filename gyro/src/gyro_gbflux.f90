!------------------------------------------------------------------
! gyro_gbflux.f90
!
! PURPOSE:
!  GyroBohm normalized forms of fluxes computed in 
!  gyro_nonlinear_flux.  Also compute averaged fluxes.
!
! NOTES:
!  Index 1: Particle flux (gyroBohm units)
!  Index 2: Energy flux (gyroBohm units)
!  Index 3: Total momentum flux (gyroBohm units)
!  Index 4: Turbulent exchange flux/length (gyroBohm units)
!------------------------------------------------------------------

subroutine gyro_gbflux

  use mpi
  use gyro_globals
  use math_constants

  !-------------------------------------------------------------
  implicit none
  !
  real, dimension(n_kinetic,n_field,p_moment,n_x) :: gbflux_i_loc
  real, dimension(n_kinetic,3) :: gbflux_mom_loc
  real, dimension(n_kinetic,2) :: gbflux_exc_loc
  real, dimension(n_kinetic,n_field,p_moment) :: gbflux_loc
  real :: gbflux_norm
  !-------------------------------------------------------------

  !-----------------------------------------------------
  ! Coefficients for spectral sum over non-negative n
  !
  if (n_1(in_1) == 0) then
     j = 1
  else
     j = 2
  endif
  !-----------------------------------------------------

  !----------------------------------------------------------
  ! Compute particle/momentum/energy diffusivities:
  !
  do ix=1,n_field
     do is=1,n_kinetic
        do i=1,n_x

           gbflux_i_loc(is,ix,:,i) = j*nonlinear_flux(i,is,ix,:)/rhos_norm**2

        enddo ! i
     enddo ! is
  enddo ! ix

  gbflux_mom_loc(:,:) = j*nonlinear_flux_momparts(:,:)/rhos_norm**2
  gbflux_exc_loc(:,:) = j*nonlinear_flux_excparts(:,:)/rhos_norm**2
  !----------------------------------------------------------

  !-----------------------------------------------------------------
  ! Define factors useful for quasilinear normalizations
  !
  call gyro_phi_kp_squared
  !-----------------------------------------------------------------

  !-----------------------------------------------------------------
  ! Renomalize linear diffusivities
  ! 
  if (nonlinear_flag == 0 .and. n_1(in_1) /= 0) then

     select case (lindiff_method)

     case (3)

        ! Gamma_{a,n} / (k_theta rho_s |phi_n|^2)
        !     Q_{a,n} / (k_theta rho_s |phi_n|^2)
        !    Pi_{a,n} / (k_theta rho_s |phi_n|^2)
        !     S_{a,n} / (k_theta rho_s |phi_n|^2)

        do i=1,n_x

           gbflux_i_loc(:,:,:,i) = nonlinear_flux(i,:,:,:)/sum(phi_squared(:)/n_x)

        enddo ! i

        gbflux_mom_loc(:,:) = nonlinear_flux_momparts(:,:)/sum(phi_squared(:)/n_x)
        gbflux_exc_loc(:,:) = nonlinear_flux_excparts(:,:)/sum(phi_squared(:)/n_x)

     end select

  endif
  !-----------------------------------------------------------------

  !----------------------------------------------------------
  ! Get (crude) averages over radius.  Do not include the 
  ! region inside the explicit damper.
  !
  gbflux_loc(:,:,:) = 0.0
  !
  do i=1+n_explicit_damp,n_x-n_explicit_damp
     gbflux_loc(:,:,:) = gbflux_loc(:,:,:)+gbflux_i_loc(:,:,:,i)
  enddo ! i
  gbflux_loc(:,:,:) = gbflux_loc(:,:,:)/(n_x-2*n_explicit_damp)
  !----------------------------------------------------------

  !----------------------------------------------------------
  ! Store flux_loc as flux_n for n-dependent diffusivity:
  !
  gbflux_n(:,:,:) = gbflux_loc(:,:,:)
  !----------------------------------------------------------

  !----------------------------------------------------------
  ! At this point, all "loc" quantities are local to 
  ! subgroups.  We need to sum these over n (i.e.,
  ! over NEW_COMM_2):
  !
  call MPI_ALLREDUCE(gbflux_mom_loc, &
       gbflux_mom, &
       size(gbflux_mom), &
       MPI_DOUBLE_PRECISION, &
       MPI_SUM, &
       NEW_COMM_2, &
       i_err)

  call MPI_ALLREDUCE(gbflux_exc_loc, &
       gbflux_exc, &
       size(gbflux_exc), &
       MPI_DOUBLE_PRECISION, &
       MPI_SUM, &
       NEW_COMM_2, &
       i_err)

  call MPI_ALLREDUCE(gbflux_loc, &
       gbflux, &
       size(gbflux), &
       MPI_DOUBLE_PRECISION, &
       MPI_SUM, &
       NEW_COMM_2, &
       i_err)

  call MPI_ALLREDUCE(gbflux_i_loc, &
       gbflux_i, &
       size(gbflux_i), &
       MPI_DOUBLE_PRECISION, &
       MPI_SUM, &
       NEW_COMM_2, &
       i_err)
  !----------------------------------------------------------

  if (nonlinear_flag == 0 .and. n_1(in_1) /= 0) then
     if (lindiff_method == 2) then

        ! [Gamma_a, Q_a, Pi_a, S_a] / Q_i normalization for linear diffusivities.

        gbflux_norm = gbflux(1,1,2)

        gbflux_mom(:,:)           = gbflux_mom(:,:)/gbflux_norm
        gbflux(:,:,:)             = gbflux(:,:,:)/gbflux_norm
        gbflux_i(:,:,:,:)         = gbflux_i(:,:,:,:)/gbflux_norm

     endif
  endif

  ! Compute running time-average (factor 0.9) of radially-dependent fluxes
  if (step > int((1.0-0.9)*nstep)) then
     p_ave = p_ave+1
     gbflux_vec(:,:,:,:) = ((p_ave-1)*gbflux_vec(:,:,:,:)+gbflux_i(:,:,:,:))/p_ave
  else
     p_ave = 0
     gbflux_vec(:,:,:,:) = 0.0
  endif

  if (i_proc == 0) then
     if (debug_flag == 1) print *,'[gyro_gbflux called]'
  endif

end subroutine gyro_gbflux
