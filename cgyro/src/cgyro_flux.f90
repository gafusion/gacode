!-----------------------------------------------------------------
! cgyro_flux.f90
!
! PURPOSE:
!  Compute flux and mean-square fluctuation amplitudes as
!  functions of (kx,ky).
!-----------------------------------------------------------------

subroutine cgyro_flux

  use mpi
  use cgyro_globals

  implicit none

  integer :: ie,ix,is,it,ir,l
  real :: dv,cn
  real :: vpar
  real :: prod
  real, dimension(n_field) :: fprod,fprod2
  real :: dvr
  real :: erot
  real :: flux_norm
  complex :: cprod, cprod2

  !-----------------------------------------------------
  ! 1. Compute kx-ky fluxes (no field breakdown)
  !-----------------------------------------------------

  flux_loc(:,:) = 0.0
  moment_loc(:,:,:,:) = 0.0
  gflux_loc(:,:,:) = 0.0


  iv_loc = 0
  do iv=nv1,nv2

     iv_loc = iv_loc+1

     is = is_v(iv)
     ix = ix_v(iv)
     ie = ie_v(iv)

     ! Integration weight
     dv = w_xi(ix)*w_e(ie)

     ! Parallel velocity
     vpar = vth(is)*sqrt(2.0)*vel(ie)*xi(ix)

     do ic=1,nc

        ir = ir_c(ic)
        it = it_c(ic)

        prod =  aimag(cap_h_c(ic,iv_loc)*conjg(psi(ic,iv_loc)))

        dvr   = w_theta(it)*dens_rot(it,is)*dens(is)*dv
        erot  = (energy(ie)+lambda_rot(it,is))*temp(is)

        ! Energy flux : Q_a
        flux_loc(ir,is) = flux_loc(ir,is)-prod*dvr*erot

        if (itp(it) > 0) then
           cprod = cap_h_c(ic,iv_loc)*dvjvec_c(1,ic,iv_loc)/z(is)
           cn    = dv*z(is)*dens(is)*dens_rot(it,is)/temp(is)

           ! Density moment: (delta n_a)/(n_norm rho_norm)
           moment_loc(ir,itp(it),is,1) = moment_loc(ir,itp(it),is,1)-(cn*field(1,ic)-cprod)

           ! Energy moment : (delta E_a)/(n_norm T_norm rho_norm)
           moment_loc(ir,itp(it),is,2) = moment_loc(ir,itp(it),is,2)-(cn*field(1,ic)-cprod)*erot
        endif

        ! Global fluxes (complex)
        do l=0,n_global

           cprod  = 0.0 
           cprod2 = 0.0

           ! i H * J0 phi^* - i H^* J0 phi

           if (ir-l > 0) then
              cprod  = cprod +i_c*cap_h_c(ic,iv_loc)*conjg(psi(ic_c(ir-l,it),iv_loc))
              cprod2 = cprod2+i_c*cap_h_c(ic,iv_loc)*conjg(i_c*chi(ic_c(ir-l,it),iv_loc))
           endif
           if (ir+l <= n_radial) then
              cprod  = cprod -i_c*conjg(cap_h_c(ic,iv_loc))*psi(ic_c(ir+l,it),iv_loc)
              cprod2 = cprod2-i_c*conjg(cap_h_c(ic,iv_loc))*(i_c*chi(ic_c(ir+l,it),iv_loc))
           endif

           gflux_loc(l,is,1) = gflux_loc(l,is,1)+cprod*dvr
           gflux_loc(l,is,2) = gflux_loc(l,is,2)+cprod*dvr*erot

           cprod = cprod*(mach*bigr(it)/rmaj+btor(it)/bmag(it)*vpar)+cprod2

           gflux_loc(l,is,3) = gflux_loc(l,is,3)+cprod*dvr*bigr(it)*mass(is)

        enddo

     enddo
  enddo

  !-----------------------------------------------------
  ! 2. Compute ky energy fluxes (with field breakdown)
  !~----------------------------------------------------

  fflux_loc(:,:,:) = 0.0

  iv_loc = 0
  do iv=nv1,nv2

     iv_loc = iv_loc+1

     is = is_v(iv)
     ix = ix_v(iv)
     ie = ie_v(iv)

     ! Integration weight
     dv = w_xi(ix)*w_e(ie)

     ! Parallel velocity
     vpar = vth(is)*sqrt(2.0)*vel(ie)*xi(ix)

     do ic=1,nc

        it = it_c(ic)

        ! Im [ H * J0 phi^* ]

        fprod(:)  = aimag(cap_h_c(ic,iv_loc)*conjg(jvec_c(:,ic,iv_loc)*field(:,ic)))
        fprod2(:) = aimag(cap_h_c(ic,iv_loc)*conjg(i_c*jxvec_c(:,ic,iv_loc)*field(:,ic)))

        dvr  = w_theta(it)*dens_rot(it,is)*dens(is)*dv
        erot = (energy(ie)+lambda_rot(it,is))*temp(is)

        ! Density flux: Gamma_a
        fflux_loc(is,1,:) = fflux_loc(is,1,:)-fprod(:)*dvr

        ! Energy flux : Q_a
        fflux_loc(is,2,:) = fflux_loc(is,2,:)-fprod(:)*dvr*erot

        fprod(:) = fprod(:)*(mach*bigr(it)/rmaj+btor(it)/bmag(it)*vpar)+fprod2(:)

        ! Momentum flux: Pi_a
        fflux_loc(is,3,:) = fflux_loc(is,3,:)-fprod(:)*dvr*bigr(it)*mass(is)

     enddo

  enddo

  !-----------------------------------------------------
  ! 3. Renormalize fluxes to GB or quasilinear forms
  !~----------------------------------------------------

  if (nonlinear_flag == 0 .and. n > 0) then

     ! Quasilinear normalization (divide by |phi|^2)
     flux_norm = 0.0
     do ir=1,n_radial
        flux_norm = flux_norm+sum(abs(field(1,ic_c(ir,:)))**2*w_theta(:))
     enddo

     ! Correct for sign of q
     flux_norm = flux_norm*q/abs(q)

     flux_loc  = flux_loc/flux_norm
     gflux_loc = gflux_loc/flux_norm 
     fflux_loc = fflux_loc/flux_norm 

  else

     ! Complete definition of fluxes (NOTE: no factor of 2 in gflux; see above)
     flux_loc  =  flux_loc*(2*k_theta*rho)
     gflux_loc = gflux_loc*(k_theta*rho)
     fflux_loc = fflux_loc*(2*k_theta*rho)

     ! GyroBohm normalizations
     flux_loc   =  flux_loc/rho**2
     gflux_loc  = gflux_loc/rho**2
     fflux_loc  = fflux_loc/rho**2

  endif

  moment_loc = moment_loc/rho

  ! Reduced real flux(kx,ky), below, is still distributed over n 

  call MPI_ALLREDUCE(flux_loc(:,:), &
       flux(:,:), &
       size(flux), &
       MPI_DOUBLE_PRECISION, &
       MPI_SUM, &
       NEW_COMM_1, &
       i_err)

  ! Reduced complex moment(kx,ky), below, is still distributed over n 

  call MPI_ALLREDUCE(moment_loc(:,:,:,:), &
       moment(:,:,:,:), &
       size(moment), &
       MPI_DOUBLE_COMPLEX, &
       MPI_SUM, &
       NEW_COMM_1, &
       i_err)

  ! Global complex gflux(l,ky), below, is still distributed over n 

  call MPI_ALLREDUCE(gflux_loc(:,:,:), &
       gflux(:,:,:), &
       size(gflux), &
       MPI_DOUBLE_COMPLEX, &
       MPI_SUM, &
       NEW_COMM_1, &
       i_err)

  ! Reduced real fflux(ky), below, is still distributed over n 

  call MPI_ALLREDUCE(fflux_loc(:,:,:), &
       fflux(:,:,:), &
       size(fflux), &
       MPI_DOUBLE_PRECISION, &
       MPI_SUM, &
       NEW_COMM_1, &
       i_err)

end subroutine cgyro_flux
