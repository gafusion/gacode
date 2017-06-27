module neo_transport

  implicit none

  public :: TRANSP_alloc, TRANSP_do, TRANSP_write

  ! transport coefficients
  real, dimension(:), allocatable :: pflux        ! (ns): gamma/(n0*vt0)
  real, dimension(:), allocatable :: eflux        ! (ns): Q/(n0*vt0*T0)
  real, dimension(:), allocatable :: mflux        ! (ns): Pi/(n0*a*T0)     
  real, dimension(:), allocatable :: uparB        ! (ns): upar*B/(vto*B0)
  real, dimension(:), allocatable :: uparBN       ! (ns): upar*B*n/vto*B0*n0
  real, dimension(:), allocatable :: klittle_upar ! (ns): upar coefficient
  real, dimension(:), allocatable :: kbig_upar    ! (ns): upar coefficient
  real, dimension(:), allocatable :: pvisc        ! (ns): B dot del dot Pi
  real                            :: jpar         ! <sum(Z*upar*B*n)>
  real                            :: jtor         ! <sum(Z*utor*n/R)>/<1/R>

  ! Sugama gyro-viscosity "H" fluxes
  real, dimension(:), allocatable :: pflux_gv      ! (ns): gamma/(n0*vt0)
  real, dimension(:), allocatable :: eflux_gv      ! (ns): Q/(n0*vt0*T0)
  real, dimension(:), allocatable :: mflux_gv      ! (ns): Pi/(n0*a*T0)  

  ! poloidal and toroidal velocities
  integer :: m_theta
  real, dimension(:,:), allocatable :: vpol, vtor, upar ! (ns,nt): vel/vt0
  real, dimension(:,:,:), allocatable :: vpol_fourier  ! (ns,mt,2)
  real, dimension(:,:,:), allocatable :: vtor_fourier  ! (ns,mt,2)
  real, dimension(:,:,:), allocatable :: upar_fourier  ! (ns,mt,2)
  real, dimension(:,:), allocatable :: vtor_0order_fourier  ! (mt,2)
  real, dimension(:,:), allocatable :: upar_0order_fourier  ! (mt,2)
  real, dimension(:),   allocatable :: vpol_th0      ! (ns): vel/vt0 at theta=0
  real, dimension(:),   allocatable :: vtor_th0      ! (ns): vel/vt0 at theta=0
  real :: vtor_0order_th0 ! 0th order toroidal(at theta=0) (species indep)
  real :: uparB_0order    ! 0th order <upar B> (species indep)

  ! potential: delta_phi(theta)
  ! NOTE: defined by sum_s Z_s e int f_s = 0
  ! -- i.e. really only valid for first-order potential
  real, dimension(:), allocatable :: d_phi       ! (ntheta)
  real                            :: d_phi_sqavg ! <d_phi^2>_theta 

  integer, parameter, private :: io=1
  character(len=80),private :: runfile_transp = 'out.neo.transport'
  character(len=80),private :: runfile_phi    = 'out.neo.phi'
  character(len=80),private :: runfile_vel    = 'out.neo.vel'
  character(len=80),private :: runfile_vel_fourier    = 'out.neo.vel_fourier'
  character(len=80),private :: runfile_exp    = 'out.neo.transport_exp'
  character(len=80),private :: runfile_gv     = 'out.neo.transport_gv'
  character(len=80),private :: runfile_flux   = 'out.neo.transport_flux'
  character(len=80),private :: runfile_check  = 'out.neo.prec'
  logical, private :: initialized = .false.
  real, private :: check_sum

contains

  subroutine TRANSP_alloc(flag)
    use neo_globals
    implicit none
    integer, intent (in) :: flag  ! flag=1: allocate; else deallocate

    if(flag == 1) then
       if(initialized) return

       allocate(pflux(n_species))
       allocate(eflux(n_species))
       allocate(mflux(n_species))
       allocate(pflux_gv(n_species))
       allocate(eflux_gv(n_species))
       allocate(mflux_gv(n_species))
       allocate(uparB(n_species))
       allocate(uparBN(n_species))
       allocate(klittle_upar(n_species))
       allocate(kbig_upar(n_species))
       allocate(pvisc(n_species))
       allocate(d_phi(n_theta))

       allocate(vpol(n_species,n_theta))
       allocate(vtor(n_species,n_theta))
       allocate(upar(n_species,n_theta))
       m_theta = (n_theta-1)/2-1
       allocate(vpol_fourier(n_species,0:m_theta,2))
       allocate(vtor_fourier(n_species,0:m_theta,2))
       allocate(upar_fourier(n_species,0:m_theta,2))
       allocate(vtor_0order_fourier(0:m_theta,2))
       allocate(upar_0order_fourier(0:m_theta,2))
       allocate(vpol_th0(n_species))
       allocate(vtor_th0(n_species))

       check_sum=0.0

       if(silent_flag == 0 .and. i_proc == 0) then
          open(unit=io,file=trim(path)//runfile_transp,status='replace')
          close(io)
          open(unit=io,file=trim(path)//runfile_phi,status='replace')
          close(io)
          open(unit=io,file=trim(path)//runfile_vel,status='replace')
          close(io)
          open(unit=io,file=trim(path)//runfile_vel_fourier,status='replace')
          close(io)
          open(unit=io,file=trim(path)//runfile_gv,status='replace')
          close(io)
          open(unit=io,file=trim(path)//runfile_flux,status='replace')
          close(io)
          if(profile_model >= 2) then
             open(unit=io,file=trim(path)//runfile_exp,status='replace')
             close(io)
          endif
       endif

       initialized = .true.

    else
       if(.NOT. initialized) return

       deallocate(pflux)
       deallocate(eflux)
       deallocate(mflux)
       deallocate(pflux_gv)
       deallocate(eflux_gv)
       deallocate(mflux_gv)
       deallocate(uparB)
       deallocate(uparBN)
       deallocate(klittle_upar)
       deallocate(kbig_upar)
       deallocate(pvisc)
       deallocate(d_phi)
       deallocate(vpol)
       deallocate(vtor)
       deallocate(upar)
       deallocate(vpol_fourier)
       deallocate(vtor_fourier)
       deallocate(upar_fourier)
       deallocate(vtor_0order_fourier)
       deallocate(upar_0order_fourier)
       deallocate(vpol_th0)
       deallocate(vtor_th0)

       if(silent_flag == 0 .and. i_proc == 0) then
          open(unit=io,file=trim(path)//runfile_check,status='replace')
          write (io,'(e16.8)',advance='no') check_sum
          close(io)
       endif

       initialized = .false.

    endif

  end subroutine TRANSP_alloc

  subroutine TRANSP_do(ir)
    use neo_globals
    use neo_energy_grid
    use neo_theory
    use neo_equilibrium
    use neo_rotation
    implicit none
    integer :: i, is, ie, ix, it
    real :: fac, fac1, fac2, rfac, poisson_F0fac, &
         B2_div_dens, bigR2_avg
    integer, intent (in) :: ir

    ! Compute the neoclassical transport coefficients

    pflux(:)  = 0.0 
    eflux(:)  = 0.0
    mflux(:)  = 0.0
    upar(:,:) = 0.0
    uparB(:)  = 0.0
    uparBN(:) = 0.0
    pvisc(:)  = 0.0
    d_phi(:)  = 0.0

    do i=1,n_row

       is = is_indx(i)
       ie = ie_indx(i)
       ix = ix_indx(i)
       it = it_indx(i)

       ! This is lambda_a (the effective potential energy)
       rfac = Z(is)/temp(is,ir) * (phi_rot(it) - phi_rot_avg) &
            - (omega_rot(ir) * bigR(it) / vth(is,ir))**2 * 0.5
       if(aniso_model(is) == 2) then
          rfac = rfac + lam_rot_aniso(is,it)
       endif

       if (ix == 0) then  
          pflux(is) = pflux(is) + w_theta(it) &
               * 4.0/sqrt(pi) * dens(is,ir) * dens_fac(is,it) * g(i) &
               * (driftx(is,it) * (4.0/3.0) * evec_e1(ie,ix) &
               + driftxrot1(is,it) * evec_e0(ie,ix))

          eflux(is) = eflux(is) + w_theta(it) * temp(is,ir) &
               * 4.0/sqrt(pi) * dens(is,ir) * dens_fac(is,it) * g(i) &
               * (driftx(is,it) * (4.0/3.0) &
               * (evec_e2(ie,ix) + rfac * evec_e1(ie,ix)) &
               + driftxrot1(is,it) * (evec_e1(ie,ix) + rfac * evec_e0(ie,ix)))

          mflux(is) = mflux(is) + w_theta(it) * temp(is,ir) &
               * dens(is,ir) * dens_fac(is,it) &
               * 4.0/sqrt(pi) * g(i) &
               * bigR(it) / vth(is,ir) &
               * (1.0/3.0 * driftxrot2(is,it) &
               * sqrt(2.0) * Btor(it)/Bmag(it) &
               * evec_e1(ie,ix) &
               + 4.0/3.0 * driftx(is,it) &
               * omega_rot(ir) * bigR(it) / vth(is,ir) &
               * evec_e1(ie,ix) &
               + driftxrot1(is,it) &
               * omega_rot(ir) * bigR(it) / vth(is,ir) &
               * evec_e0(ie,ix))

          d_phi(it) = d_phi(it) + Z(is) *  dens(is,ir) &
               * dens_fac(is,it) * g(i) &
               * 4.0/sqrt(pi) * evec_e0(ie,ix)

       else if (ix == 1) then

          pflux(is) = pflux(is) + w_theta(it) &
               * dens(is,ir) * dens_fac(is,it) &
               * 4.0/sqrt(pi) * g(i) &
               * driftxrot2(is,it) * (1.0/3.0) * evec_e05(ie,ix)

          eflux(is) = eflux(is) + w_theta(it) * temp(is,ir) &
               * dens(is,ir) * dens_fac(is,it) &
               * 4.0/sqrt(pi) * g(i) &
               * driftxrot2(is,it) * (1.0/3.0) &
               * (evec_e105(ie,ix) + rfac * evec_e05(ie,ix))

          mflux(is) = mflux(is) + w_theta(it) &
               * 4.0/sqrt(pi) * dens(is,ir) * dens_fac(is,it) * g(i) &
               * temp(is,ir) / vth(is,ir) * bigR(it) &
               * (8.0/15.0 * driftx(is,it) * Btor(it)/Bmag(it) * sqrt(2.0) &
               * evec_e105(ie,ix) &
               + 1.0/3.0 * driftxrot1(is,it) &
               * sqrt(2.0) * Btor(it)/Bmag(it) * evec_e05(ie,ix) &
               + 1.0/3.0 * driftxrot2(is,it) &
               * omega_rot(ir) * bigR(it) / vth(is,ir) * evec_e05(ie,ix) &
               + 2.0/15.0 * driftxrot3(is,it)* evec_e105(ie,ix))

          ! uparB = < B * 1/n * int vpar * (F0 g)>
          uparB(is) = uparB(is) + w_theta(it) &
               * Bmag(it) * sqrt(2.0) * vth(is,ir) &
               * (1.0/3.0) * g(i) &
               * 4.0/sqrt(pi) * evec_e05(ie,ix)

          uparBN(is) = uparBN(is) + w_theta(it) &
               * Bmag(it) * sqrt(2.0) * vth(is,ir) &
               * (1.0/3.0) * g(i) &
               * 4.0/sqrt(pi) * dens(is,ir) * dens_fac(is,it) * evec_e05(ie,ix)

          upar(is,it) = upar(is,it) &
               + sqrt(2.0) * vth(is,ir) &
               * (1.0/3.0) * g(i) &
               * 4.0/sqrt(pi) * evec_e05(ie,ix)

       else if (ix == 2) then
          pflux(is) = pflux(is) + w_theta(it) &
               * 4.0/sqrt(pi) * dens(is,ir) * dens_fac(is,it) * g(i) &
               * driftx(is,it) * (2.0/15.0) * evec_e1(ie,ix)

          eflux(is) = eflux(is) + w_theta(it) * temp(is,ir) &
               * 4.0/sqrt(pi) * dens(is,ir) * dens_fac(is,it) * g(i) &
               * driftx(is,it) * (2.0/15.0) &
               * (evec_e2(ie,ix) + rfac * evec_e1(ie,ix))

          mflux(is) = mflux(is) + w_theta(it) * temp(is,ir) &
               * dens(is,ir) * dens_fac(is,it) &
               * 4.0/sqrt(pi) * g(i) &
               * bigR(it) / vth(is,ir) &
               * (2.0/15.0 * driftxrot2(is,it) &
               * sqrt(2.0) * Btor(it)/Bmag(it) &
               * evec_e1(ie,ix) &
               + 2.0/15.0 * driftx(is,it) &
               * omega_rot(ir) * bigR(it) / vth(is,ir) &
               * evec_e1(ie,ix))

       else if (ix == 3) then
          mflux(is) = mflux(is) + w_theta(it) &
               * 4.0/sqrt(pi) * dens(is,ir) * dens_fac(is,it) * g(i) &
               * temp(is,ir) / vth(is,ir) * bigR(it) &
               * 2.0/35.0 * evec_e105(ie,ix) &
               * (driftx(is,it) * Btor(it)/Bmag(it) * sqrt(2.0) &
               - driftxrot3(is,it))

       endif
    enddo

    do is=1, n_species
       eflux(is) = eflux(is) + omega_rot(ir) * mflux(is)
    enddo

    ! d_phi: sum_s Z_s e int f_s = 0
    ! (really only valid for first-order, i.e. local solution f_1)
    ! sum Z^2 n_0s / T0s * d_phi = sum Z int F0s*g
    if(rotation_model == 2) then
       ! phi_1 not computed for case with rotation effects included
       d_phi(:) = 0.0
       d_phi_sqavg = 0.0
    else
       ! Compute the poisson int F0 factor: sum Zs^2/T0s*(int of F0s)
       poisson_F0fac = 0.0
       do is=1, n_species
          poisson_F0fac = poisson_F0fac &
               + dens(is,ir) * Z(is) * Z(is) / temp(is,ir)
       enddo
       if(adiabatic_ele_model == 1) then
          poisson_F0fac = poisson_F0fac + ne_ade(ir) / te_ade(ir)
       endif
       d_phi_sqavg = 0.0
       do it=1, n_theta
          d_phi(it) = d_phi(it) / poisson_F0fac
          d_phi_sqavg = d_phi_sqavg + d_phi(it) * d_phi(it) * w_theta(it)
       enddo
    endif

    ! Bootstrap current = sum <Z*n*upar B>
    jpar = 0.0
    do is=1, n_species
       jpar = jpar + Z(is) * uparBN(is) 
    enddo

    ! U_parallel coefficient
    do is=1, n_species
       ! <B^2 / n_0a>
       B2_div_dens = 0.0
       bigR2_avg = 0.0
       do it=1,n_theta
          B2_div_dens = B2_div_dens + w_theta(it) &
               * Bmag(it)**2 / (dens(is,ir) * dens_fac(is,it))
          bigR2_avg = bigR2_avg + w_theta(it) * bigR(it)**2
       enddo
       kbig_upar(is) = (1.0 / B2_div_dens) &
            * ( uparB(is) &
            - (I_div_psip * rho(ir) * temp(is,ir) / (z(is)*1.0)) &
            * (dlnndr(is,ir) &
            - (z(is)*1.0)/temp(is,ir) * dphi0dr(ir) &
            + dlntdr(is,ir) & 
            * (1.0 + z(is) / temp(is,ir) * phi_rot_avg &
            + omega_rot(ir)**2 * 0.5/vth(is,ir)**2 &
            * (bigR_th0**2 - bigR2_avg)) &
            + omega_rot_deriv(ir) &
            * omega_rot(ir)/vth(is,ir)**2 * (bigR_th0**2 - bigR2_avg) &
            + omega_rot(ir)**2 * bigR_th0 / vth(is,ir)**2 * bigR_th0_rderiv))
       if(aniso_model(is) == 2) then
          kbig_upar(is) = kbig_upar(is) &
               - (I_div_psip * rho(ir) * temp(is,ir) / (z(is)*1.0)) &
               * (dlntdr(is,ir) * lam_rot_avg_aniso(is))
       endif
       if(abs(dlntdr(is,ir)) > epsilon(0.)) then
          klittle_upar(is) = -kbig_upar(is) *  B2_div_dens &
               / (dlntdr(is,ir) * &
               I_div_psip * rho(ir) * temp(is,ir) / (z(is)*1.0))
       else
          klittle_upar(is) = 0.0
       end if
    enddo

    ! Poloidal and Toroidal Velocity
    call compute_velocity(ir)

    ! Toroidal component of Bootstrap current = sum <Z*n*utor/R>/<1/R>
    jtor = 0.0
    do is=1, n_species
       do it=1,n_theta
          jtor = jtor + Z(is) * dens(is,ir) * dens_fac(is,it) * vtor(is,it) &
               * w_theta(it) / bigR(it)
       enddo
    enddo
    jtor = jtor/bigRinv_avg
    
    if(silent_flag == 0 .and. i_proc == 0) then
       open(unit=io_neoout,file=trim(path)//runfile_neoout,&
            status='old',position='append')
       fac=0.0
       do is=1, n_species
          fac = fac + Z(is) * pflux(is)
       enddo
       write(io_neoout,*)
       write(io_neoout,'(t2,a,t21,e14.5)') 'r/a = ', r(ir)
       write(io_neoout,'(t2,a,t21,e14.5)') 'jpar = ', jpar
       write(io_neoout,'(t2,a,t21,e14.5)') 'sum Z_s Gamma_s = ', fac
       write(io_neoout,'(t3,a,t11,a,t25,a,t39,a)') &
            'Z', 'pflux', 'eflux', 'mflux'
       do is=1,n_species
          write (io_neoout,'(f7.3,3(e14.5))') &
               Z(is),pflux(is),eflux(is),mflux(is)
       enddo
       write(io_neoout,*) '------------------'
       write(io_neoout,*)
       close(io_neoout)
    endif

    ! Sugama gyro-viscosity "H" fluxes
    pflux_gv(:)  = 0.0
    eflux_gv(:)  = 0.0
    mflux_gv(:)  = 0.0
    do is=1,n_species
       fac1 = 0.0
       fac2 = 0.0
       do it=1,n_theta
          rfac = Z(is)/temp(is,ir) * (phi_rot(it) - phi_rot_avg) &
               - (omega_rot(ir) * bigR(it) / vth(is,ir))**2 * 0.5
          if(aniso_model(is) == 2) then
             rfac = rfac + lam_rot_aniso(is,it)
          endif
          fac1 = fac1 + w_theta(it) * dens(is,ir)  * dens_fac(is,it) &
               / Bmag(it)**3 * (2.0*gradr(it) * k_par(it) * gradr_tderiv(it) &
               - 1.0/Bmag(it) * gradr(it)**2 * gradpar_Bmag(it))
          fac2 = fac2 + w_theta(it) * dens(is,ir)  * dens_fac(is,it) &
               / Bmag(it)**3 * (2.0 * gradr(it) * k_par(it) &
               * gradr_tderiv(it) &
               - 1.0/Bmag(it) * gradr(it)**2 * gradpar_Bmag(it)) &
               * (1.0 + rfac)
       enddo
       pflux_gv(is) = -0.5 * rho(ir)**2 * mass(is) &
            * temp(is,ir) / (Z(is)*1.0)**2 * I_div_psip * r(ir) / q(ir) &
            * fac1 * omega_rot_deriv(ir)
       ! EAB: 02/05/14 fixed bug in mflux_gv -- dlntdr bigR_th0**2 term
       ! had wrong sign
       mflux_gv(is) =  -0.5 * temp(is,ir) * rho(ir)**2 * mass(is) &
            * temp(is,ir) / (Z(is)*1.0)**2 * I_div_psip * r(ir) / q(ir) &
            * (fac2 * dlntdr(is,ir) + fac1 * (dlnndr(is,ir) &
            - z(is)/temp(is,ir) * dphi0dr(ir) &
            + omega_rot(ir) * bigR_th0**2 / vth(is,ir)**2 &
            * omega_rot_deriv(ir) &
            + omega_rot(ir)**2 * bigR_th0 / vth(is,ir)**2 &
            * bigR_th0_rderiv &
            + dlntdr(is,ir) * (1.0 + 0.5 * omega_rot(ir)**2 * bigR_th0**2  &
            / vth(is,ir)**2) &
            + dlntdr(is,ir) * z(is)/temp(is,ir) * phi_rot_avg) )
       eflux_gv(is) = -0.5  * temp(is,ir) * rho(ir)**2 * mass(is) &
            * temp(is,ir) / (Z(is)*1.0)**2 * I_div_psip * r(ir) / q(ir) &
            * fac2 * omega_rot_deriv(ir) &
            + 2.5 * pflux_gv(is) * temp(is,ir) &
            + omega_rot(ir) * mflux_gv(is)
    enddo

    do is=1,n_species
       check_sum = check_sum &
            +  (abs(pflux(is)) + abs(eflux(is)) + abs(mflux(is)))/rho(ir)**2 &
            + abs(uparB(is))/rho(ir)
    enddo

  end subroutine TRANSP_do

  subroutine compute_velocity(ir)
    use neo_globals
    use neo_rotation
    use neo_equilibrium
    implicit none
    integer, intent (in) :: ir
    integer:: is, it, jt
    real :: RBt

    ! 0th-order toroidal flow (at theta=0) and <u_par B>
    vtor_0order_th0  = omega_rot(ir) * bigR_th0 
    RBt = 0.0
    do it=1,n_theta
       RBt = RBt + w_theta(it) * bigR(it) * Btor(it)
    enddo
    uparB_0order = omega_rot(ir) * RBt

    do is=1, n_species
       do it=1, n_theta
          vpol(is,it) = kbig_upar(is) &
               * Bpol(it) / (dens(is,ir) * dens_fac(is,it))
          ! EAB: 02/05/14 -- fixed bug in vtor -- missing bigR_th0_rderiv term
          vtor(is,it) = kbig_upar(is) &
               * Btor(it) / (dens(is,ir) * dens_fac(is,it)) &
               + (I_div_psip * rho(ir) * temp(is,ir) / (z(is)*1.0)) &
               * (1.0 / Btor(it)) &
               * (dlnndr(is,ir) &
               - (1.0*z(is))/temp(is,ir) * dphi0dr(ir) &
               + dlntdr(is,ir) & 
               * (1.0 + z(is) / temp(is,ir) * phi_rot(it) &
               + omega_rot(ir)**2 * 0.5/vth(is,ir)**2 &
               * (bigR_th0**2 - bigR(it)**2)) &
               + omega_rot(ir)**2 * bigR_th0 &
                     / vth(is,ir)**2 * bigR_th0_rderiv &
               + omega_rot_deriv(ir) * omega_rot(ir)/vth(is,ir)**2 &
               * (bigR_th0**2 - bigR(it)**2))
          if(aniso_model(is) == 2) then
             vtor(is,it) = vtor(is,it) &
                  + (I_div_psip * rho(ir) * temp(is,ir) / (z(is)*1.0)) &
                  * (1.0 / Btor(it)) &
                  * (dlntdr(is,ir) * lam_rot_aniso(is,it))
          endif
       end do

       vpol_th0(is) = 0.0
       vtor_th0(is) = 0.0
       do jt=0, m_theta
          vpol_fourier(is,jt,:) = 0.0
          vtor_fourier(is,jt,:) = 0.0
          upar_fourier(is,jt,:) = 0.0
          do it=1, n_theta
             vpol_fourier(is,jt,1) = vpol_fourier(is,jt,1) &
                  + vpol(is,it) * cos(jt * theta(it))
             vtor_fourier(is,jt,1) = vtor_fourier(is,jt,1) &
                  + vtor(is,it) * cos(jt * theta(it))
             upar_fourier(is,jt,1) = upar_fourier(is,jt,1) &
                  + upar(is,it) * cos(jt * theta(it))
             if(jt > 0) then
                vpol_fourier(is,jt,2) = vpol_fourier(is,jt,2) &
                     + vpol(is,it) * sin(jt * theta(it))
                vtor_fourier(is,jt,2) = vtor_fourier(is,jt,2) &
                     + vtor(is,it) * sin(jt * theta(it))
                upar_fourier(is,jt,2) = upar_fourier(is,jt,2) &
                     + upar(is,it) * sin(jt * theta(it))
             endif
          enddo
          if(jt == 0) then
             vpol_fourier(is,jt,1) = vpol_fourier(is,jt,1) / (1.0*n_theta)
             vtor_fourier(is,jt,1) = vtor_fourier(is,jt,1) / (1.0*n_theta)
             upar_fourier(is,jt,1) = upar_fourier(is,jt,1) / (1.0*n_theta)
          else
             vpol_fourier(is,jt,:) = vpol_fourier(is,jt,:) / (0.5*n_theta)
             vtor_fourier(is,jt,:) = vtor_fourier(is,jt,:) / (0.5*n_theta)
             upar_fourier(is,jt,:) = upar_fourier(is,jt,:) / (0.5*n_theta)
          end if
          ! vel(theta=0) = sum vel_fourier cos coefficients
          vpol_th0(is) = vpol_th0(is) + vpol_fourier(is,jt,1)
          vtor_th0(is) = vtor_th0(is) + vtor_fourier(is,jt,1)
       enddo

    enddo

    do jt=0, m_theta
       vtor_0order_fourier(jt,:) = 0.0
       upar_0order_fourier(jt,:) = 0.0
       do it=1, n_theta
          vtor_0order_fourier(jt,1) = vtor_0order_fourier(jt,1) &
               + omega_rot(ir) * bigR(it) &
               * cos(jt * theta(it))
          upar_0order_fourier(jt,1)= upar_0order_fourier(jt,1) &
               + omega_rot(ir) * bigR(it) * Btor(it) / Bmag(it) &
               * cos(jt * theta(it))
          if(jt > 0) then
             vtor_0order_fourier(jt,2) = vtor_0order_fourier(jt,2) &
                  + omega_rot(ir) * bigR(it) &
                  * sin(jt * theta(it))
             upar_0order_fourier(jt,2)= upar_0order_fourier(jt,2) &
                  + omega_rot(ir) * bigR(it) * Btor(it) / Bmag(it) &
                  * sin(jt * theta(it))
          end if
       enddo
       if(jt == 0) then
          vtor_0order_fourier(jt,1) = vtor_0order_fourier(jt,1) / (1.0*n_theta)
          upar_0order_fourier(jt,1) = upar_0order_fourier(jt,1) / (1.0*n_theta)
       else
          vtor_0order_fourier(jt,:) = vtor_0order_fourier(jt,:) / (0.5*n_theta)
          upar_0order_fourier(jt,:) = upar_0order_fourier(jt,:) / (0.5*n_theta)
       end if
    enddo

  end subroutine compute_velocity

  subroutine TRANSP_write(ir)
    use neo_globals
    use neo_rotation
    implicit none
    integer, intent (in) :: ir
    integer :: is, jt, is_ele
    real :: pgb, egb, mgb, dens_ele, temp_ele

    if(silent_flag > 0 .or. i_proc > 0) return

    ! transport coefficients (normalized)
    open(io,file=trim(path)//runfile_transp,status='old',position='append')
    write (io,'(e16.8)',advance='no') r(ir)
    write (io,'(e16.8)',advance='no') d_phi_sqavg
    write (io,'(e16.8)',advance='no') jpar
    write (io,'(e16.8)',advance='no') vtor_0order_th0
    write (io,'(e16.8)',advance='no') uparB_0order
    do is=1, n_species
       write (io,'(e16.8)',advance='no') pflux(is)
       write (io,'(e16.8)',advance='no') eflux(is)
       write (io,'(e16.8)',advance='no') mflux(is)
       write (io,'(e16.8)',advance='no') uparB(is)
       write (io,'(e16.8)',advance='no') klittle_upar(is)
       write (io,'(e16.8)',advance='no') kbig_upar(is)
       write (io,'(e16.8)',advance='no') vpol_th0(is)
       write (io,'(e16.8)',advance='no') vtor_th0(is)
    enddo
    write (io,*)
    close(io)

    ! transport coefficients (units)
    if(profile_model >= 2) then
       open(io,file=trim(path)//runfile_exp,status='old',position='append')
       write (io,'(e16.8)',advance='no') r(ir) * a_meters
       ! m
       write (io,'(e16.8)',advance='no') &
            d_phi_sqavg * (temp_norm(ir)*temp_norm_fac/charge_norm_fac)**2
       ! V^2
       write (io,'(e16.8)',advance='no') &
            jpar * (charge_norm_fac*dens_norm(ir)*vth_norm(ir)*a_meters)
       ! A/m^2 (jpar B/Bunit)
       write (io,'(e16.8)',advance='no') &
            vtor_0order_th0 * (vth_norm(ir)*a_meters)
       ! m/s
       write (io,'(e16.8)',advance='no') &
            uparB_0order * (vth_norm(ir)*a_meters)
       ! m/s (upar B/Bunit)
       do is=1, n_species
          write (io,'(e16.8)',advance='no') &
               pflux(is) * (dens_norm(ir)*vth_norm(ir)*a_meters) 
          ! e19 m-2 s-1
          write (io,'(e16.8)',advance='no') &
               eflux(is) * (dens_norm(ir)*vth_norm(ir)*a_meters &
               *temp_norm(ir)*temp_norm_fac)
          ! W/m^2 
          write (io,'(e16.8)',advance='no') & 
               mflux(is) * (dens_norm(ir)*a_meters*temp_norm(ir)*temp_norm_fac)
          ! N/m
          write (io,'(e16.8)',advance='no') &
               uparB(is) * (vth_norm(ir)*a_meters)
          ! m/s (upar B/Bunit)
          write (io,'(e16.8)',advance='no') klittle_upar(is)  ! dimensionless
          write (io,'(e16.8)',advance='no') &
               kbig_upar(is) * (dens_norm(ir)*vth_norm(ir)*a_meters/b_unit(ir))
          ! e19 1/(m^2 s T)
          write (io,'(e16.8)',advance='no') &
               vpol_th0(is) * (vth_norm(ir)*a_meters)
          ! m/s
          write (io,'(e16.8)',advance='no') &
               vtor_th0(is) * (vth_norm(ir)*a_meters)
          ! m/s
       enddo
       write (io,*)
       close(io)
    end if

    ! delta phi(theta)
    open(io,file=trim(path)//runfile_phi,status='old',position='append')
    write(io,*) d_phi(:)
    close(io)

    ! u_par(theta)
    open(io,file=trim(path)//runfile_vel,status='old',position='append')    
    write(io,*) upar(:,:)
    close(io)

    ! u_par, vpol, vtor theta fourier coefficients
    open(io,file=trim(path)//runfile_vel_fourier,status='old',position='append')
    do is=1,n_species
       do jt=0,m_theta
          write (io,'(e16.8)',advance='no') upar_fourier(is,jt,1)
       enddo
       do jt=0,m_theta
          write (io,'(e16.8)',advance='no') upar_fourier(is,jt,2)
       enddo
       do jt=0,m_theta
          write (io,'(e16.8)',advance='no') vpol_fourier(is,jt,1)
       enddo
       do jt=0,m_theta
          write (io,'(e16.8)',advance='no') vpol_fourier(is,jt,2)
       enddo
       do jt=0,m_theta
          write (io,'(e16.8)',advance='no') vtor_fourier(is,jt,1)
       enddo
       do jt=0,m_theta
          write (io,'(e16.8)',advance='no') vtor_fourier(is,jt,2)
       enddo
    enddo
    write (io,*)
    close(io)

    ! gyroviscosity transport coefficients
    open(io,file=trim(path)//runfile_gv,status='old',position='append')
    write (io,'(e16.8)',advance='no') r(ir)
    do is=1, n_species
       write (io,'(e16.8)',advance='no') pflux_gv(is)
       write (io,'(e16.8)',advance='no') eflux_gv(is)
       write (io,'(e16.8)',advance='no') mflux_gv(is)
    enddo
    write (io,*)
    close(io)

    ! fluxes in GB units
    if(adiabatic_ele_model == 0) then
       do is=1, n_species
          if(Z(is) < 0.0) then
             is_ele = is
             exit
          endif
       enddo
       dens_ele = dens(is_ele,ir)
       temp_ele = temp(is_ele,ir)
    else
       dens_ele = ne_ade(ir)
       temp_ele = te_ade(ir)
    endif
    pgb = dens_ele * rho(ir)**2 * temp_ele**1.5
    egb = dens_ele * rho(ir)**2 * temp_ele**2.5
    mgb = dens_ele * rho(ir)**2 * temp_ele**2
    open(io,file=trim(path)//runfile_flux,status='old',position='append')
    write (io,'(a,e16.8)') '# r/a=', r(ir)
    write(io,'(a,t3,a,t11,a,t25,a,t39,a)') &
         '#', 'Z', 'pflux_dke', 'eflux_dke', 'mflux_dke'
    write(io,'(a,t11,a,t25,a,t39,a)') &
         '#','(GB)', '(GB)', '(GB)'
    do is=1,n_species
       write (io,'(f7.3)',advance='no') Z(is)
       write (io,'(e14.5)',advance='no') pflux(is)/pgb
       write (io,'(e14.5)',advance='no') eflux(is)/egb
       write (io,'(e14.5)',advance='no') mflux(is)/mgb
       write (io,*)
    enddo
    write(io,'(a, t3,a,t11,a,t25,a,t39,a)') &
         '#', 'Z', 'pflux_gv', 'eflux_gv', 'mflux_gv'
     write(io,'(a,t11,a,t25,a,t39,a)') &
         '#', '(GB)', '(GB)', '(GB)'
    do is=1,n_species
       write (io,'(f7.3)',advance='no') Z(is)
       write (io,'(e14.5)',advance='no') pflux_gv(is)/pgb
       write (io,'(e14.5)',advance='no') eflux_gv(is)/egb
       write (io,'(e14.5)',advance='no') mflux_gv(is)/mgb
       write(io,*)
    enddo
    write(io,'(a, t3,a,t11,a,t25,a,t39,a)') &
         '#', 'Z', 'pflux_tgyro', 'eflux_tgyro', 'mflux_tgyro'
     write(io,'(a,t11,a,t25,a,t39,a)') &
         '#', '(GB)', '(GB)', '(GB)'
    do is=1,n_species
       write (io,'(f7.3)',advance='no') Z(is)
       write (io,'(e14.5)',advance='no') (pflux(is)+pflux_gv(is))/pgb
       write (io,'(e14.5)',advance='no') (eflux(is)+eflux_gv(is) &
            -omega_rot(ir)*mflux(is)-omega_rot(ir)*mflux_gv(is))/egb
       write (io,'(e14.5)',advance='no') (mflux(is)+mflux_gv(is))/mgb
       write (io,*)
    enddo
    close(io)
    
  end subroutine TRANSP_write

end module neo_transport
