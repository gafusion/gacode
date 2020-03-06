!---------------------------------------------------------------------------
! cgyro_flux.f90
!
! PURPOSE:
!  Compute flux and mean-square fluctuation amplitudes as
!  functions of (kx,ky).
!
! NOTES:
!
!  Flux definitions (these are all broken down by ky) 
!
!   gflux(l,ky) = Radial Fourier expansion of flux (l is the Fourier harmonic)
!   gflux(0,ky) = Domain-averaged flux [this is the standard flux-tube flux]
!     cflux(ky) = Interior-averaged flux [this is the "central" flux]
!
!  Global flux logic:
!
!   Retain radial harmonics l = [0,...,N_GLOBAL], where N_GLOBAL=4 by default
!
!   To trigger output of the global flux (viewable by -plot xflux),
!    set GFLUX_PRINT_FLAG=1
!
! NOTE:
!  cflux is only of interest when profile or ExB shear is active
!  cflux(ky) is the average flux in the "positive-shear region"
!  gflux(ky,0)-cflux is the flux in the "negative-shear region"
!---------------------------------------------------------------------------

subroutine cgyro_flux

  use mpi
  use cgyro_globals
  use cgyro_interface

  implicit none

  integer :: ie,ix,is,it,ir
  integer :: l,icl
  real :: dv,cn
  real :: vpar
  complex, dimension(0:n_global,n_field) :: prod1,prod2
  real :: dvr
  real :: erot
  real :: flux_norm
  complex :: cprod
  real, parameter :: x_fraction=0.2
  real :: u

  !-----------------------------------------------------
  ! 1. Compute kx-ky moments (n,E)
  !-----------------------------------------------------

  moment_loc(:,:,:,:) = 0.0

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

        erot  = (energy(ie)+lambda_rot(it,is))*temp(is)

        if (itp(it) > 0) then
           cprod = cap_h_c(ic,iv_loc)*dvjvec_c(1,ic,iv_loc)/z(is)
           cn    = dv*z(is)*dens(is)*dens_rot(it,is)/temp(is)

           ! Density moment: (delta n_a)/(n_norm rho_norm)
           moment_loc(ir,itp(it),is,1) = moment_loc(ir,itp(it),is,1)-(cn*field(1,ic)-cprod)

           ! Energy moment : (delta E_a)/(n_norm T_norm rho_norm)
           moment_loc(ir,itp(it),is,2) = moment_loc(ir,itp(it),is,2)-(cn*field(1,ic)-cprod)*erot

           ! Velocity moment : (delta v_a)/(n_norm v_norm rho_norm)
           moment_loc(ir,itp(it),is,3) = moment_loc(ir,itp(it),is,3)-(cn*field(1,ic)-cprod)*vpar
        endif

     enddo
  enddo

  !-------------------------------------------------------------
  ! 2. Compute global ky-dependent fluxes (with field breakdown)
  !-------------------------------------------------------------

  gflux_loc(:,:,:,:) = 0.0
  cflux_loc(:,:,:) = 0.0

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

        prod1 = 0.0 
        prod2 = 0.0

        ! Global fluxes (complex)
        do l=0,n_global

           ! i H J0 phi^* - i H^* J0 phi

           if (ir-l > 0) then
              icl = ic_c(ir-l,it)
              prod1(l,:) = prod1(l,:)+i_c*cap_h_c(ic,iv_loc)*&
                   conjg(jvec_c(:,icl,iv_loc)*field(:,icl))
              prod2(l,:) = prod2(l,:)+i_c*cap_h_c(ic,iv_loc)*&
                   conjg(i_c*jxvec_c(:,icl,iv_loc)*field(:,icl))
           endif
           if (ir+l <= n_radial) then
              icl = ic_c(ir+l,it)
              prod1(l,:) = prod1(l,:)-i_c*conjg(cap_h_c(ic,iv_loc))*&
                   jvec_c(:,icl,iv_loc)*field(:,icl)
              prod2(l,:) = prod2(l,:)-i_c*conjg(cap_h_c(ic,iv_loc))*&
                   i_c*jxvec_c(:,icl,iv_loc)*field(:,icl)
           endif

        enddo

        dvr  = w_theta(it)*dens_rot(it,is)*dens(is)*dv
        erot = (energy(ie)+lambda_rot(it,is))*temp(is)

        ! Density flux: Gamma_a
        gflux_loc(:,is,1,:) = gflux_loc(:,is,1,:)+prod1(:,:)*dvr

        ! Energy flux : Q_a
        gflux_loc(:,is,2,:) = gflux_loc(:,is,2,:)+prod1(:,:)*dvr*erot

        prod1(:,:) = prod1(:,:)*(mach*bigr(it)/rmaj+btor(it)/bmag(it)*vpar)+prod2(:,:)

        ! Momentum flux: Pi_a
        gflux_loc(:,is,3,:) = gflux_loc(:,is,3,:)+prod1(:,:)*dvr*bigr(it)*mass(is)

        ! Construct "positive/interior" flux:
        cflux_loc = real(gflux_loc(0,:,:,:))
        do l=1,n_global
           u = 2*pi*l*x_fraction
           cflux_loc = cflux_loc+2*sin(u)*real(gflux_loc(l,:,:,:))/u
        enddo

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
     flux_norm = flux_norm*q/abs(q)*2 ! need 2 for regression compatibility

     gflux_loc = gflux_loc/flux_norm 
     cflux_loc = cflux_loc/flux_norm 

  else

     ! Complete definition of fluxes
     gflux_loc = gflux_loc*(k_theta*rho)
     cflux_loc = cflux_loc*(k_theta*rho)

     ! GyroBohm normalizations
     gflux_loc  = gflux_loc/rho**2
     cflux_loc  = cflux_loc/rho**2

  endif

  moment_loc = moment_loc/rho

  ! Reduced complex moment(kx,ky), below, is still distributed over n 

  call MPI_ALLREDUCE(moment_loc(:,:,:,:), &
       moment(:,:,:,:), &
       size(moment), &
       MPI_DOUBLE_COMPLEX, &
       MPI_SUM, &
       NEW_COMM_1, &
       i_err)

  ! Global complex gflux(l,ky), below, is still distributed over n 

  call MPI_ALLREDUCE(gflux_loc, &
       gflux, &
       size(gflux), &
       MPI_DOUBLE_COMPLEX, &
       MPI_SUM, &
       NEW_COMM_1, &
       i_err)

  ! Reduced real cflux(ky), below, is still distributed over n 

  call MPI_ALLREDUCE(cflux_loc, &
       cflux, &
       size(cflux), &
       MPI_DOUBLE_PRECISION, &
       MPI_SUM, &
       NEW_COMM_1, &
       i_err)

  ! Reorder output for QLGYRO
  do l=1,3
     cgyro_gbflux_out(l, :, :) = real(gflux(0, :, l, :))
  end do


end subroutine cgyro_flux
