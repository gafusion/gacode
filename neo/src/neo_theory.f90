module neo_theory

  ! Contains routines to compute analytical neoclassical theory values.
  ! NOTE: The theories assume the diamagnetic ordering.  Thus, rotation
  ! effects are not included here.

  implicit none

  public :: THEORY_alloc, THEORY_do
  real :: pflux_HH, efluxi_HH, efluxe_HH, jpar_HH, kpar_HH, uparB_HH, &
       vpol_ion_HH,  efluxi_CH, efluxi_TG, &
       jpar_S, kpar_S, uparB_S, vpol_ion_S, jtor_S, phi_HR, jpar_K
  real, dimension(:), allocatable :: pflux_multi_HS, eflux_multi_HS
  
  real, private :: eps ! r/rmaj
  real, private :: EparB_avg     ! <E_|| B>
  real, private :: nui_HH, nue_HH, nui_star_HH, nue_star_HH
  integer, private :: is_ele, is_ion
  real, private :: dens_ele, temp_ele, zeff
  integer, private :: ir_global, is_global, ietype
  integer, parameter, private :: io=40
  character(len=80),private :: runfile = 'out.neo.theory'
  logical, private :: initialized = .false.

contains

  subroutine THEORY_alloc(flag)
    use neo_globals
    use neo_nclass_dr
    implicit none
    integer, intent (in) :: flag  ! flag=1: allocate; else deallocate
    integer :: is

    if(flag == 1) then
       if(initialized) return
       if(adiabatic_ele_model == 0) then
          do is=1, n_species
             if(Z(is) < 0.0) then
                is_ele = is
                exit
             endif
          enddo
       endif

       allocate(pflux_multi_HS(n_species))
       allocate(eflux_multi_HS(n_species))

       if(silent_flag == 0 .and. i_proc == 0) then
          open(io,file=trim(path)//runfile,status='replace')
          close(io)
       end if

       if(sim_model == 1 .or. sim_model == 3) then
          call NCLASS_DR_alloc(1)
       endif
       
       initialized = .true.

    else
       if(.NOT. initialized) return
       deallocate(pflux_multi_HS)
       deallocate(eflux_multi_HS)
       if(sim_model == 1 .or. sim_model == 3) then
          call NCLASS_DR_alloc(0)
       endif
       initialized = .false.
    end if

  end subroutine THEORY_alloc

  ! Driver routine for computing theory results
  ! Theories assume no rotation -- so use input dens(psi) as dens
  ! (i.e. not with poloidal part)
  subroutine THEORY_do(ir)
    use neo_globals
    use neo_equilibrium
    use neo_nclass_dr
    implicit none
    integer, intent (in) :: ir
    integer :: is
    real :: d_max

    ! inverse aspect ratio
    eps = r(ir) / rmaj(ir)
    
    ! majority ion species and zeff
           
    if(adiabatic_ele_model == 0) then
       dens_ele = dens(is_ele,ir)
       temp_ele = temp(is_ele,ir)
    else
       dens_ele = ne_ade(ir)
       temp_ele = te_ade(ir)
    endif
    is_ion = -1
    d_max = -1.0
    do is=1, n_species
       if(Z(is) > 0 .and. dens(is,ir) > d_max) then
          is_ion = is
          d_max  = dens(is,ir) 
       endif
    enddo
    if(is_ion == -1) then
       if(silent_flag==0 .and. i_proc==0) then
          open(unit=io_neoout,file=trim(path)//runfile_neoout,&
               status='old',position='append')
          write(io_neoout,*) 'Warning: No ion species for neo_theory'
          close(io_neoout)
       endif
       ! assume primary ion species is is=1
       is_ion = 1
    endif
    
    zeff = 0.0
    do is=1,n_species
       if(is /= is_ele) then
          zeff = zeff + dens(is,ir) * z(is)**2 / dens_ele
       endif
    enddo
    
    ! theory coll freqs
    nui_HH = nu(is_ion,ir) * (4.0/3.0) / sqrt(2.0*pi)
    nui_star_HH = nui_HH * rmaj(ir) * abs(q(ir)) &
         / (sqrt(eps)*sqrt(eps)*sqrt(eps) * vth(is_ion,ir))
    if(adiabatic_ele_model == 1) then
       nue_HH = 0.0
       nue_star_HH = 0.0
    else
       nue_HH = nu(is_ion,ir) * (4.0/3.0) / sqrt(pi) &
            * sqrt(mass(is_ion)/mass(is_ele)) &
            * (temp(is_ion,ir) / temp_ele)**1.5 &
            * 1.0 / (1.0* Z(is_ion))**2
       nue_star_HH = nue_HH * rmaj(ir) * abs(q(ir)) &
            / (sqrt(eps)*sqrt(eps)*sqrt(eps) * vth(is_ele,ir))
    endif

    ! inductive electric field
    ! < E_|| B>
    EparB_avg = epar0(ir)

    ! compute the theory
    call compute_HH(ir,pflux_HH, efluxi_HH, efluxe_HH, &
         jpar_HH, kpar_HH, uparB_HH)
    call compute_vpol_ion(ir,kpar_HH, vpol_ion_HH)
    call compute_CH(ir,efluxi_CH)
    call compute_TG(ir,efluxi_TG)
    call compute_Sauter(ir,jpar_S, kpar_S, uparB_S,jtor_S)
    call compute_vpol_ion(ir,kpar_S, vpol_ion_S)
    call compute_HR(ir,phi_HR)
    call compute_HS(ir,pflux_multi_HS, eflux_multi_HS)
    call compute_Koh(ir,jpar_K)

    if(silent_flag == 0 .and. i_proc == 0) then
       open(io,file=trim(path)//runfile,status='old',position='append')
       write(io,'(e16.8)',advance='no') r(ir)
       write(io,'(e16.8)',advance='no') pflux_HH
       write(io,'(e16.8)',advance='no') efluxi_HH
       write(io,'(e16.8)',advance='no') efluxe_HH
       write(io,'(e16.8)',advance='no') jpar_HH
       write(io,'(e16.8)',advance='no') kpar_HH
       write(io,'(e16.8)',advance='no') uparB_HH
       write(io,'(e16.8)',advance='no') vpol_ion_HH
       write(io,'(e16.8)',advance='no') efluxi_CH
       write(io,'(e16.8)',advance='no') efluxi_TG
       write(io,'(e16.8)',advance='no') jpar_S
       write(io,'(e16.8)',advance='no') kpar_S
       write(io,'(e16.8)',advance='no') uparB_S
       write(io,'(e16.8)',advance='no') vpol_ion_S
       write(io,'(e16.8)',advance='no') phi_HR
       do is=1, n_species
          write(io,'(e16.8)',advance='no') pflux_multi_HS(is)
          write(io,'(e16.8)',advance='no') eflux_multi_HS(is)
       enddo
       write(io,'(e16.8)',advance='no') jpar_K
       write (io,*)
       close(io)
    end if

    if(sim_model == 1 .or. sim_model == 3) then
       call NCLASS_DR_do(ir)
    endif

  end subroutine THEORY_do

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Hinton-Rosenbluth potential
  ! Phys. Fluids, vol. 16, 836 (1973)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine compute_HR(ir,phi)
    use neo_globals
    implicit none
    integer, intent (in) :: ir
    real, intent (out) :: phi
    real :: thavg_fac=0.5        ! theta-avg factor <sin^2> = 0.5
    real ::  nu_crit, phi_1, phi_2, fac

    fac = (q(ir)/eps) * eps * dlntdr(is_ion,ir) * rho(ir) &
         * (sqrt(temp(is_ion,ir) * mass(is_ion)) / (1.0 * Z(is_ion))) &
         / (1.0*Z(is_ion) / temp(is_ion,ir) + 1.0/temp_ele)
    
    ! banana regime -- phi_1 = phi_1 * nui_HH
    phi_1 = 1.81 * nui_star_HH / nui_HH
    ! plateau regime
    phi_2 = sqrt(pi/2.0) 
  
    nu_crit = phi_2 / phi_1

    if(nui_HH > nu_crit) then
       phi = phi_2 * phi_2 * thavg_fac * fac**2
    else
       phi = (phi_1 * nui_HH) * (phi_1 * nui_HH) * thavg_fac * fac**2
    endif
    
  end subroutine compute_HR
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Hinton-Hazeltine fluxes and flows
  ! Rev. Mod. Phys., vol. 48, 239 (1976)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine compute_HH(ir,pflux, efluxi, efluxe, jpar, kpar, uparB)
    use neo_globals
    use neo_equilibrium, only : I_div_psip, Bmag2_avg, Bmag2inv_avg
    implicit none
    integer, intent (in) :: ir
    real, intent (out) :: pflux, efluxi, efluxe, jpar, kpar, uparB
    real :: k0, a0, b0, c0, fac, beta1, A1e, K11, K12, K22, K2, K13, K23, K33,&
         gi, qi, sigma_spitzer

    fac = (1.17-0.35*sqrt(nui_star_HH)) / (1+0.7*sqrt(nui_star_HH))
    beta1 = (fac - 2.1*nui_star_HH*nui_star_HH*eps*eps*eps) &
         / (1.0 + nui_star_HH*nui_star_HH*eps*eps*eps)
    
    kpar = beta1
    
    ! Note: only jpar and upar B use I_div_psip rather than q/eps
    uparB  = I_div_psip * rho(ir) * temp(is_ion,ir) / (z(is_ion)*1.0) &
         * (dlnndr(is_ion,ir) &
         - (z(is_ion)*1.0)/temp(is_ion,ir) * dphi0dr(ir) &
         + (1.0 - kpar) * dlntdr(is_ion,ir))

    k0=0.66; a0=1.03; b0=0.31; c0=0.74
    K2 = k0 * ( 1.0/(1 + a0*sqrt(nui_star_HH) + b0*nui_star_HH) &
         + eps*eps*eps * (c0*c0/b0) * nui_star_HH &
         / (1 + c0 * sqrt(eps*eps*eps) * nui_star_HH) )
    
    ! Note: ion Q get additional part below if not adiabatic ele
    qi = (q(ir)/eps)**2 &
         * dens(is_ion,ir) * temp(is_ion,ir) &
         * (2.0 * mass(is_ion) * temp(is_ion,ir) / (1.0 *Z(is_ion)**2)) &
         * rho(ir)**2 * sqrt(eps) * nui_HH * K2 * dlntdr(is_ion,ir)
    
    if(adiabatic_ele_model == 1) then
       pflux  = 0.0
       efluxe = 0.0
       jpar   = 0.0
       efluxi = qi
       
    else
       
       A1e = -(dlnndr(is_ele,ir) + dlntdr(is_ele,ir)) &
            + 2.5 * dlntdr(is_ele,ir) &
            + temp(is_ion,ir) / (temp(is_ele,ir) * Z(is_ion)) &
              * ( -(dlnndr(is_ion,ir) + dlntdr(is_ion,ir)) &
                  + beta1 / (1.0 + nue_star_HH*nue_star_HH*eps*eps) &
                    * dlntdr(is_ion,ir) )

       k0=1.04; a0=2.01; b0=1.53; c0=0.89
       K11 = k0 * ( 1.0/(1 + a0*sqrt(nue_star_HH) + b0*nue_star_HH) &
            + eps*eps*eps * (c0*c0/b0) * nue_star_HH  &
            / (1 + c0 * sqrt(eps*eps*eps) * nue_star_HH) )
       
       k0=1.20; a0=0.76; b0=0.67; c0=0.56
       K12 =  k0 * ( 1.0/(1 + a0*sqrt(nue_star_HH) + b0*nue_star_HH) &
            + eps*eps*eps * (c0*c0/b0) * nue_star_HH &
            / (1 + c0 * sqrt(eps*eps*eps) * nue_star_HH) )

       k0=2.55; a0=0.45; b0=0.43; c0=0.43
       K22 = k0 * ( 1.0/(1 + a0*sqrt(nue_star_HH) + b0*nue_star_HH) &
            + eps*eps*eps * (c0*c0/b0) * nue_star_HH &
            / (1 + c0 * sqrt(eps*eps*eps) * nue_star_HH) )

       k0=2.30; a0=1.02; b0=1.07; c0=1.07
       K13 = k0 * ( 1.0 / (1 + a0*sqrt(nue_star_HH) + b0* nue_star_HH) ) &
            * ( 1.0 / (1 + c0 * sqrt(eps*eps*eps) * nue_star_HH) )

       k0=4.19; a0=0.57; b0=0.61; c0=0.61
       K23 = k0 * ( 1.0 / (1 + a0*sqrt(nue_star_HH) + b0*nue_star_HH) ) &
            * ( 1.0 / (1 + c0  * sqrt(eps*eps*eps) * nue_star_HH) )
       
       k0=1.83; a0=0.68; b0=0.32; c0=0.66
       K33 = k0 * ( 1.0 / (1 + a0*sqrt(nue_star_HH) + b0*nue_star_HH) ) &
            * ( 1.0 / (1 + c0  * sqrt(eps*eps*eps) * nue_star_HH) )

       gi  = -(q(ir)/eps)**2 * dens(is_ele,ir) &
            * (2.0*mass(is_ele)*temp(is_ele,ir)) *rho(ir)**2 &
            * sqrt(eps) * nue_HH * (K11*A1e - K12*dlntdr(is_ele,ir)) &
            - EparB_avg  * dens(is_ele,ir) * sqrt(eps) * rho(ir) * I_div_psip &
            * K13

       pflux  = gi

       efluxi = qi &
            - beta1 * temp(is_ion,ir) / (1.0 *Z(is_ion)) * gi &
            / (1.0 + nue_star_HH*nue_star_HH*eps*eps*eps) &
            + 2.5 * temp(is_ion,ir) * gi

       efluxe = -(q(ir)/eps)**2 &
            * dens(is_ele,ir) * temp(is_ele,ir) &
            * (2.0*mass(is_ele)*temp(is_ele,ir)) *rho(ir)**2 &
            * sqrt(eps) * nue_HH * (K12*A1e - K22*dlntdr(is_ele,ir)) &
            - EparB_avg  * dens(is_ele,ir) * sqrt(eps) * rho(ir) * I_div_psip &
            * temp(is_ele,ir) * K23
       

       ! Note: only jpar and upar B and inductive E-field 
       ! use I_div_psip rather than q/eps

       sigma_spitzer = dens_ele / ( mass(is_ele) * nue_HH ) &
            / (0.29 + 0.46/(1.08+z(is_ion)))

       jpar   = -(I_div_psip) * dens(is_ele,ir) * temp(is_ele,ir) &
            * sqrt(eps) * rho(ir) * (K13 * A1e - K23*dlntdr(is_ele,ir)) &
            - EparB_avg * sigma_spitzer * sqrt(eps) * K33 &
            + EparB_avg * sigma_spitzer * Bmag2_avg * Bmag2inv_avg

    end if

  end subroutine compute_HH


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Chang-Hinton ion heat flux (q_i) (Q = q_i + 2.5 * particle flux)
  ! Phys. Plasmas, vol. 25, 1493 (1982)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine compute_CH(ir,efluxi)
    use neo_globals
    implicit none
    integer, intent (in) :: ir
    real, intent (out) :: efluxi
    real :: k0, a0, b0, c0, K1, K2, F2, &
         CH_Bmag2inv_avg, CH_Bmag2avg_inv, CH_I_div_psip, alpha, mu_star

    k0=0.66; a0=1.03; b0=0.31; c0=0.74

    alpha  = zeff - 1.0
    mu_star = (1.0 + 1.54*alpha) * nui_star_HH

    CH_Bmag2inv_avg = (1.0 + 1.5*(eps*eps + eps*shift(ir)) &
         + 0.375 * eps*eps*eps * shift(ir)) &
         / (1.0 + 0.5*eps*shift(ir))
    CH_Bmag2avg_inv = (sqrt(1.0-eps*eps) * (1.0 + 0.5*eps*shift(ir))) &
         / (1.0 + (shift(ir)/eps) * (sqrt(1.0 - eps*eps) - 1.0))
    CH_I_div_psip   = q(ir)/eps

    F2 = (0.5/sqrt(eps)) * (CH_Bmag2inv_avg - CH_Bmag2avg_inv)

    K1 = -alpha * (0.83 + 0.42*alpha) / (0.58 + alpha) &
         * (c0 * mu_star * sqrt(eps*eps*eps) * F2) &
         / (1.0 + c0 * mu_star *  sqrt(eps*eps*eps))

    K2 = (k0*(1.0+1.54*alpha) &
         + (1.88*sqrt(eps) - 1.54*eps)*(1.0+3.75*alpha)) * CH_Bmag2inv_avg &
         / (1 + a0*sqrt(mu_star) + b0*mu_star) &
         + k0 * sqrt(eps*eps*eps) * (c0*c0/b0) * mu_star * F2 &
         * (1.0 + 1.33*alpha * (1.0+0.6*alpha)/(1.0+1.79*alpha)) &
         / (1.0 + c0 * sqrt(eps*eps*eps) * mu_star)

    efluxi = CH_I_div_psip**2 &
         * dens(is_ion,ir) * temp(is_ion,ir) &
         * (2.0*mass(is_ion)*temp(is_ion,ir) / (1.0 *Z(is_ion)**2)) &
         * rho(ir)**2 * sqrt(eps) * nui_HH &
         * ( (K2 + K1) * dlntdr(is_ion,ir) + K1 * dlnndr(is_ion,ir))
    
  end subroutine compute_CH

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Taguchi ion heat flux (q_i) (Q = q_i + 2.5 * particle flux)
  ! modified with Chang-Hinton collisional interpolation factor
  ! PPCF, vol. 30, 1897 (1988)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine compute_TG(ir,efluxi)
    use neo_globals
    use neo_equilibrium, only: ftrap
    implicit none
    integer, intent (in) :: ir
    real, intent (out) :: efluxi
    real :: k0, a0, b0, c0, K2, K2_star, F2, &
         CH_Bmag2inv_avg, CH_Bmag2avg_inv, CH_I_div_psip, fc

    ! Geo factors     
    ! original Chang-Hinton
    CH_Bmag2inv_avg = (1.0 + 1.5*(eps*eps + eps*shift(ir)) &
         + 0.375 * eps*eps*eps * shift(ir)) &
         / (1.0 + 0.5*eps*shift(ir))
    CH_Bmag2avg_inv = (sqrt(1.0-eps*eps) * (1.0 + 0.5*eps*shift(ir))) &
         / (1.0 + (shift(ir)/eps) * (sqrt(1.0 - eps*eps) - 1.0))
    CH_I_div_psip   = q(ir)/eps

    ! Taguchi thermal conductivity
    fc = 1.0 - ftrap

    K2_star = (CH_Bmag2inv_avg - CH_Bmag2avg_inv) &
         + CH_Bmag2avg_inv * (1.0-fc)/(1.0+1.17*fc)

    ! Chang-Hinton collisional interpolation factor
    k0=0.66; a0=1.03; b0=0.31; c0=0.74

    F2 = (0.5/sqrt(eps)) * (CH_Bmag2inv_avg - CH_Bmag2avg_inv)
    
    K2 = k0 *( (K2_star/k0) / (1 + a0*sqrt(nui_star_HH) + b0*nui_star_HH) &
         + sqrt(eps*eps*eps) * (c0*c0/b0) * nui_star_HH * F2 &
         / (1.0 + c0 * sqrt(eps*eps*eps) * nui_star_HH) )

    ! sqrt(eps) factor is already in K2_star
    efluxi = CH_I_div_psip**2 &
         * dens(is_ion,ir) * temp(is_ion,ir) &
         *   (2.0*mass(is_ion)*temp(is_ion,ir) / (1.0 *Z(is_ion)**2)) &
         * rho(ir)**2 * nui_HH * K2 * dlntdr(is_ion,ir)
    
  end subroutine compute_TG

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Sauter bootstrap current model
  ! Phys. Plasmas, vol. 6, 2834 (1999)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine compute_Sauter(ir,jpar, kpar, uparB, jtor)
    use neo_globals
    use neo_equilibrium, only: I_div_psip, ftrap, Btor, w_theta
    implicit none
    integer, intent (in) :: ir
    real, intent (out) :: jpar, kpar, uparB, jtor
    real :: X31, L31_S, X32e, F32_ee, X32i, F32_ei, L32_S, alpha_0, alpha_S, &
         X34, L34_S, X33, sigma_S, sigma_spitzer
    real :: nue_S, nue_star_S, nui_star_S
    integer :: is, it
    real :: dens_sum, press_sum, fac

    ! EAB: 07/11/2013 changed def of Sauter nui_star_S for impurities
    ! from (Zi**2 * Zeff * ne) to (Zi**4)*(n_all_ions) 
    !nui_star_S = nui_star_HH / (dens(is_ion,ir) * z(is_ion)**2) &
    !     * dens_ele*zeff
    press_sum = 0.0
    dens_sum = 0.0
    do is=1,n_species
       if(is /= is_ele) then
          dens_sum  = dens_sum + dens(is,ir)
          press_sum = press_sum + dens(is,ir)*temp(is,ir)
       endif
    enddo
    nui_star_S = nui_star_HH / dens(is_ion,ir) * dens_sum
    ! Osborne's interpretation
    !nui_star_S = nui_star_HH / (dens(is_ion,ir) * z(is_ion)**4) &
    !     * dens_sum * zeff**4

    alpha_0 = -1.17*(1.0-ftrap) / (1.0 - 0.22*ftrap - 0.19*ftrap*ftrap)
    
    alpha_S = ( (alpha_0 + 0.25*(1.0-ftrap*ftrap) * sqrt(nui_star_S)) &
         / (1.0 + 0.5 * sqrt(nui_star_S)) &
         + (0.315) * nui_star_S * nui_star_S  &
         * ftrap*ftrap*ftrap*ftrap*ftrap*ftrap ) &
         / (1.0 + 0.15 *nui_star_S*nui_star_S &
         * ftrap*ftrap*ftrap*ftrap*ftrap*ftrap)
    
    kpar = -1.0 * alpha_S

    uparB  = I_div_psip * rho(ir) * temp(is_ion,ir) / (z(is_ion)*1.0) &
         * (dlnndr(is_ion,ir) &
         - (z(is_ion)*1.0)/temp(is_ion,ir) * dphi0dr(ir) &
         + (1.0 - kpar) * dlntdr(is_ion,ir))

    if(adiabatic_ele_model == 1) then
       jpar = 0.0
       jtor = 0.0
    else

       nue_S  = nue_HH * (dens_ele/dens(is_ion,ir)) &
            * (1.0*zeff) / (1.0* Z(is_ion))**2 
       nue_star_S = nue_S * rmaj(ir) * abs(q(ir)) &
            / (sqrt(eps)*sqrt(eps)*sqrt(eps) * vth(is_ele,ir))
       
       X31    = ftrap / (1.0 + (1.0-0.1*ftrap) * sqrt(nue_star_S)  &
            + 0.5*(1.0-ftrap) * nue_star_S / (1.0*zeff))
       
       L31_S  = (1.0 + 1.4/(zeff+1)) * X31 &
            - (1.9/(zeff+1)) * X31*X31 &
            + (0.3/(zeff+1)) * X31*X31*X31 &
            + (0.2/(zeff+1)) * X31*X31*X31*X31
       
       X32e   = ftrap / (1.0 + 0.26*(1-ftrap) * sqrt(nue_star_S) &
            + 0.18*(1.0-0.37*ftrap) * nue_star_S / sqrt(1.0*zeff))
       
       F32_ee = (0.05 + 0.62*zeff) &
            / (zeff*(1+0.44*zeff)) * (X32e - X32e*X32e*X32e*X32e) &
            + 1.0/(1+0.22*zeff) * (X32e*X32e - X32e*X32e*X32e*X32e &
            - 1.2*(X32e*X32e*X32e - X32e*X32e*X32e*X32e)) &
            + 1.2/(1+0.5*zeff) * X32e*X32e*X32e*X32e
       
       X32i   = ftrap / (1.0 + (1+0.6*ftrap) * sqrt(nue_star_S)  &
            + 0.85*(1-0.37*ftrap) * nue_star_S*(1.0+zeff))
       
       F32_ei = -(0.56 + 1.93*zeff) & 
            /(zeff*(1+0.44*zeff)) * (X32i - X32i*X32i*X32i*X32i)  &
            + 4.95/(1+2.48*zeff) * (X32i*X32i - X32i*X32i*X32i*X32i &
            - 0.55*(X32i*X32i*X32i - X32i*X32i*X32i*X32i)) &
            + (-1.2)/(1+0.5*zeff) * X32i*X32i*X32i*X32i;
       
       L32_S  = F32_ee + F32_ei
       
       X34   = ftrap / (1.0 + (1.0-0.1*ftrap) * sqrt(nue_star_S) &
            + 0.5*(1.0-0.5*ftrap) * nue_star_S/(1.0*zeff))
       
       L34_S = (1 + 1.4/(zeff+1.0)) * X34 &
            - (1.9/(zeff+1.0)) * X34*X34 &
            + (0.3/(zeff+1.0)) * X34*X34*X34 &
            + (0.2/(zeff+1.0)) * X34*X34*X34*X34
       
       X33 = ftrap / (1.0 + (0.55-0.1*ftrap) * sqrt(nue_star_S) &
            + 0.45*(1.0-ftrap) * nue_star_S/(1.0*zeff**1.5))

       ! EAB note 05/09/11 
       ! -- Sauter spizer is a little different than HH spitzer
       ! HH: L11 = 1.0  / (0.29 + 0.46/(1.08+z(is_ion)))
       ! S : L11 = 0.58 *32 / (3.0*pi) / (0.58+0.74/(0.76+z(is_ion)))
       
       sigma_spitzer = dens_ele / ( mass(is_ele) * nue_S ) &
            * 0.58 * 32 / (3.0*pi) / (0.58 + 0.74/(0.76+zeff))

       sigma_S = sigma_spitzer &
            * (1.0 - (1.0+0.36/(1.0*zeff)) * X33 &
            + 0.59/(1.0*zeff) * X33 * X33 &
            - 0.23/(1.0*zeff) * X33 * X33 * X33)

       jpar = sigma_S * EparB_avg &
            + L32_S * I_div_psip * rho(ir) * &
            dens(is_ele,ir) * temp(is_ele,ir) * dlntdr(is_ele,ir)

       do is=1,n_species
          jpar = jpar + L31_S * I_div_psip * rho(ir) &
               * dens(is,ir) * temp(is,ir) &
               * (dlntdr(is,ir) + dlnndr(is,ir))
          !if(is /= is_ele) then
          !   jpar = jpar + L34_S * alpha_S * I_div_psip * rho(ir) &
          !        * dlntdr(is,ir) * dens(is,ir) * temp(is,ir)
          !endif
       enddo
       ! EAB: 08/15/13 changed interpretation for ions with unlike temps
       ! from the above loop to this
       jpar = jpar + L34_S * alpha_S * I_div_psip * rho(ir) &
            * dlntdr(is_ion,ir) * press_sum

       ! Approximation for <jtor/R>/<1/R> ~ <jpar B/Bunit>/<BT/Bunit>
       fac = 0.0
       do it=1,n_theta
          fac = fac + w_theta(it) * Btor(it)
       enddo
       jtor = jpar/fac
       
    end if

  end subroutine compute_Sauter

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Hirshman-Sigmar fluxes
  ! Phys. Fluids, vol. 20, 418 (1977)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine compute_HS(ir,pflux_multi, eflux_multi)
    use neo_globals
    use neo_equilibrium, only: I_div_psip, Bmag2_avg, ftrap
    implicit none
    integer, intent (in) :: ir
    real, intent (inout) :: pflux_multi(n_species), eflux_multi(n_species)
    real, dimension(:), allocatable :: nux0, nux2, nux4
    real :: A1, A2, omega_fac, HS_I_div_psip, &
         L11, L12, L21, L22, L_a, L_b, sum_nm, eii_val
    integer :: Nx=100
    integer, parameter :: integ_order = 64
    integer :: js

    omega_fac = 1.0 / Bmag2_avg
    HS_I_div_psip = I_div_psip

    allocate(nux0(n_species))
    allocate(nux2(n_species))
    allocate(nux4(n_species))

    ir_global = ir
    do is_global=1, n_species
       do ietype=1, 3
          call gauss_integ(-1.0,1.0,myHSenefunc,integ_order,Nx,eii_val)
          if(ietype == 1) then
             nux0(is_global) = eii_val * 4.0 / (3.0 * sqrt(pi))
          else if(ietype == 2) then
             nux2(is_global) = eii_val * 4.0 / (3.0 * sqrt(pi))
          else if(ietype == 3) then
             nux4(is_global) = eii_val * 4.0 / (3.0 * sqrt(pi))
          endif
       enddo
    enddo

    sum_nm = 0.0
    do is_global=1, n_species
       sum_nm = sum_nm + mass(is_global)*dens(is_global,ir)*nux0(is_global)
    enddo

    do is_global=1, n_species
             
       A1 = -dlnndr(is_global,ir) + 1.5*dlntdr(is_global,ir)
       A2 = -dlntdr(is_global,ir)
       pflux_multi(is_global) = 0.0
       eflux_multi(is_global) = 0.0

       L_a = nux0(is_global) * omega_fac * HS_I_div_psip**2 &
            * rho(ir)**2 * ftrap * dens(is_global,ir) &
            * temp(is_global,ir) * mass(is_global) &
            / (Z(is_global) * Z(is_global) * 1.0)

       do js=1, n_species

          L_b = nux0(js) * omega_fac * HS_I_div_psip**2 &
               * rho(ir)**2 * ftrap * dens(js,ir) &
               * temp(js,ir) * mass(js) / (Z(js) * Z(js) * 1.0)

          if(is_global == js) then
             L11 = -L_a * (sum_nm &
                  - mass(is_global) * dens(is_global,ir) * nux0(is_global)) &
                  / sum_nm
             L12 = L11 * (nux2(is_global) / nux0(is_global))
             L22 = (nux2(is_global) / nux0(is_global)) &
                  * (L12 &
                  + L_a * (nux2(is_global) / nux0(is_global) &
                  - nux4(is_global) / nux2(is_global)))
             L21 = L12
             
          else       
             L11 = L_a * (Z(is_global) * temp(js,ir)) &
                  / (Z(js) * temp(is_global,ir)) &
                  * mass(js) * dens(js,ir) * nux0(js) &
                  / sum_nm
             L12 = (nux2(js) / nux0(js)) * L11
             L21 = (nux2(is_global) / nux0(is_global)) &
                  * Z(js)/(1.0*Z(is_global)) &
                  * (mass(is_global) * dens(is_global,ir) * nux0(is_global) &
                  / sum_nm) * L_b
             L22 = (nux2(is_global) / nux0(is_global)) * L12
             
          end if


          pflux_multi(is_global) = pflux_multi(is_global) &
               + L11 * A1 + L12 * A2
          eflux_multi(is_global) = eflux_multi(is_global) &
               + L21 * A1 + L22 * A2
          
       enddo
    enddo

    deallocate(nux0)
    deallocate(nux2)
    deallocate(nux4)

  end subroutine compute_HS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Ion Poloidal Flow at theta=0 (from kpar)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine compute_vpol_ion(ir, kpar, vpol_ion)
    use neo_globals
    use neo_equilibrium, only: I_div_psip, Bmag2_avg, Bpol
    implicit none
    real, intent (out) :: vpol_ion
    real, intent (in) :: kpar
    integer, intent (in) :: ir
    integer :: it, jt, m_theta
    real :: bigk, sum
    real, dimension(:), allocatable :: vpol_ion_th

    bigk = -kpar *  (dlntdr(is_ion,ir) * I_div_psip * rho(ir) &
         * temp(is_ion,ir) / (z(is_ion)*1.0)) &
         / (Bmag2_avg / dens(is_ion,ir))
    
    allocate(vpol_ion_th(n_theta))
    do it=1,n_theta
       vpol_ion_th(it) = bigk * Bpol(it) / dens(is_ion,ir)
    end do

    vpol_ion = 0.0
    m_theta = (n_theta-1)/2-1
    do jt=0, m_theta
       ! compute the fourier cosine series coefficient c_n:
       ! vpol = sum_n c_n * cos(n*theta)
       sum = 0.0
       do it=1, n_theta
          sum = sum +  vpol_ion_th(it) * cos(jt * theta(it))
       enddo
       if(jt == 0) then
          sum = sum / (1.0*n_theta)
       else
          sum = sum / (0.5*n_theta)
       end if
       ! vpol(theta=0) = sum_n c_n
       vpol_ion = vpol_ion + sum
    enddo

    deallocate(vpol_ion_th)

  end subroutine compute_vpol_ion

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Function Defininitions
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! uses ir_global, is_global, ietype, eps
  real function myHSenefunc(x)

    use neo_globals
    use neo_equilibrium, only : ftrap
    implicit none
    integer :: js
    real :: x, xa, xb
    real :: ene, de, val
    real,parameter :: emin=0.0, emax=16.0
    real :: ft_fac, ft_star, nu_d_tot
    real :: nu_d

    ! x grid: -1 -> 1 (equally spaced) for int emin..emax
    ! x = 2 (v/vth)/sqrt(2*emax)-1
    xa = 2.0 / (1.0 - sqrt(emin / emax))
    xb = - (1.0 + sqrt(emin / emax)) &
         / (1.0 - sqrt(emin / emax))
    ene = emax * ( (x-xb) / xa )**2
    de = 2.0 * sqrt(emax) / xa
    val = de * exp(-ene) ! weight w/o ene factor

    ! nu_d * ene 
    nu_d_tot = 0.0
    do js=1, n_species
       call get_coll_freqs(0, ir_global,is_global,js,ene,nu_d)
       nu_d_tot = nu_d_tot + nu_d * nu(is_global,ir_global)
    enddo

    ft_star = (3.0*pi/16.0) * eps**2 * vth(is_global,ir_global) &
         * sqrt(2.0) * ene**1.5 &
         / (rmaj(ir_global) * q(ir_global) * nu_d_tot)
    ft_fac  = 1.0 / (1.0 + ftrap / ft_star)

    ! interpolated
    if(ietype == 1) then
       myHSenefunc = val * nu_d_tot * ene * ft_fac
    else if(ietype == 2) then
       myHSenefunc = val * nu_d_tot * ene * ene * ft_fac
    else if(ietype == 3) then
       myHSenefunc = val * nu_d_tot * ene * ene * ene * ft_fac

    endif

  end function myHSenefunc

  ! Compute the energy-dependent coll freqs for a given ir, is, js
  ! returns nu -> (nu * ene)
  subroutine get_coll_freqs(flag, ir_loc,is_loc,js_loc,ene,nu_d)
    use neo_globals
    implicit none
    integer, intent(in)  :: flag
    ! 0->use HS0 model; else->use collision_model
    integer, intent(in) :: ir_loc, is_loc, js_loc
    real, intent(in)    :: ene
    real, intent(inout) :: nu_d
    real :: xa, xb, Hd_coll, Xd_coll, Hs_coll, Xs_coll, fac

    ! (Note: nu(is,ir) and pol part of dens from rotation will be added 
    !  in coll term in kinetic equation)
    fac = (1.0*Z(js_loc))**2 / (1.0*Z(is_loc))**2 &
         * dens(js_loc,ir_loc)/dens(is_loc,ir_loc)
    xa = sqrt(ene)
    xb = xa * (vth(is_loc,ir_loc)/vth(js_loc,ir_loc))

    if (flag .ne. 0 .and. &         
         (collision_model == 1 .or. &
         (spitzer_model==1 .and. js_loc .ne. is_loc)) ) then
       ! Connor model: nu_d = nu_s, form dependent on (is,js)
       if ((is_loc == js_loc) .or. &
            (abs(mass(is_loc) - mass(js_loc)) < epsilon(0.))) then
          ! case 1: like-species/same mass collisions
          if(xb < 1e-4) then
             ! special case limit ene->0 
             Hd_coll = (1.0/sqrt(pi)) * &
                  (4.0/3.0    *  (vth(is_loc,ir_loc)/vth(js_loc,ir_loc)) &
                  - 4.0/15.0  *  (vth(is_loc,ir_loc)/vth(js_loc,ir_loc))**3 &
                  * ene &
                  + 2.0/35.0  *  (vth(is_loc,ir_loc)/vth(js_loc,ir_loc))**5 &
                  * ene**2 &
                  - 2.0/189.0 *  (vth(is_loc,ir_loc)/vth(js_loc,ir_loc))**7 &
                  * ene**3)
             Xd_coll = 1.0
          else
             Hd_coll = exp(-xb*xb) / (xb*sqrt(pi)) &
                  + (1.0-1.0/(2.0*xb*xb)) * erf(xb)
             Xd_coll = 1.0 / (xa)
             endif
       else if(mass(is_loc) < mass(js_loc)) then
          ! case 2: ele-ion and ion-imp(heavy) collisions
          Hd_coll = 1.0
          Xd_coll = 1.0 / (xa)
       else
          ! case 3: ion-ele and imp(heavy)-ion collisions
          Hd_coll = 1.0
          Xd_coll = 1.0 * (xa**2)
          fac = fac * 4.0/(3.0*sqrt(pi)) &
               * sqrt(mass(js_loc) / mass(is_loc)) &
               * (temp(is_loc,ir_loc) / temp(js_loc,ir_loc))**1.5
       endif
       nu_d = fac * Hd_coll * Xd_coll

    else if(xb < 1e-4) then
       ! Hirshman-Sigmar model: nu_d != nu_s
       ! special case limit ene->0 
       ! return nu_d * ene (finite); nu_s * ene (zero)

       nu_d = fac * (1.0/sqrt(pi)) * &
            (4.0/3.0    * (vth(is_loc,ir_loc)/vth(js_loc,ir_loc)) &
            - 4.0/15.0  * (vth(is_loc,ir_loc)/vth(js_loc,ir_loc))**3 * ene & 
            + 2.0/35.0  * (vth(is_loc,ir_loc)/vth(js_loc,ir_loc))**5 * ene**2 &
            - 2.0/189.0 * (vth(is_loc,ir_loc)/vth(js_loc,ir_loc))**7 * ene**3)
    else
       ! Hirshman-Sigmar model: nu_d != nu_s
       Hd_coll = exp(-xb*xb) / (xb*sqrt(pi)) & 
            + (1.0-1.0/(2.0*xb*xb)) * erf(xb)
       Xd_coll = 1.0 / (xa)
       Hs_coll = -exp(-xb*xb) / (xb*sqrt(pi)) & 
            + (1.0/(2.0*xb*xb)) * erf(xb)
       Xs_coll = (xa) * (2.0 * temp(is_loc,ir_loc)/temp(js_loc,ir_loc)) &
            * (1.0 + mass(js_loc) / mass(is_loc))
       nu_d = fac * Hd_coll * Xd_coll
    endif

  end subroutine get_coll_freqs

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Koh modification to
  ! Sauter bootstrap current model
  ! Phys. Plasmas, vol. 19, 072505 (2012)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine compute_Koh(ir,jpar)
    use neo_globals
    use neo_equilibrium, only: I_div_psip, ftrap
    implicit none
    integer, intent (in) :: ir
    real, intent (out) :: jpar
    real :: X31, L31_S, X32e, F32_ee, X32i, F32_ei, L32_S, alpha_0, alpha_S, &
         X34, L34_S, X33, sigma_S, sigma_spitzer
    real :: nue_S, nue_star_S, nui_star_S
    real :: ftrap_new, delta_param, alpha_param, beta_param, h_param, pfac
    integer :: h_flag = 2
    integer :: is
    real :: dens_sum, press_sum

    if(adiabatic_ele_model == 1) then
       jpar = 0.0

    else

       ! alpha unchanged from Sauter
       ! EAB: 07/11/2013 changed def of Sauter nui_star_S for impurities
       ! from (Zi**2 * Zeff * ne) to (Zi**4)*(n_all_ions) 
       !nui_star_S = nui_star_HH / (dens(is_ion,ir) * z(is_ion)**2) &
       !     * dens_ele*zeff
       dens_sum = 0.0
       press_sum = 0.0
       do is=1,n_species
          if(is /= is_ele) then
             dens_sum = dens_sum + dens(is,ir)
             press_sum = press_sum + dens(is,ir)*temp(is,ir)
          endif
       enddo
       nui_star_S = nui_star_HH / dens(is_ion,ir) * dens_sum
       
       alpha_0 = -1.17*(1.0-ftrap) / (1.0 - 0.22*ftrap - 0.19*ftrap*ftrap)
       
       alpha_S = ( (alpha_0 + 0.25*(1.0-ftrap*ftrap) * sqrt(nui_star_S)) &
            / (1.0 + 0.5 * sqrt(nui_star_S)) &
            + (0.315) * nui_star_S * nui_star_S  &
            * ftrap*ftrap*ftrap*ftrap*ftrap*ftrap ) &
            / (1.0 + 0.15 *nui_star_S*nui_star_S &
            * ftrap*ftrap*ftrap*ftrap*ftrap*ftrap)

       nue_S  = nue_HH * (dens_ele/dens(is_ion,ir)) &
            * (1.0*zeff) / (1.0* Z(is_ion))**2 
       nue_star_S = nue_S * rmaj(ir) * abs(q(ir)) &
            / (sqrt(eps)*sqrt(eps)*sqrt(eps) * vth(is_ele,ir))

       if(zeff > 5.0) then
          alpha_param = 0.0
       else
          alpha_param = (-zeff**2 + 5.998*zeff - 4.981) & 
               / (4.294*zeff**2 - 14.07*zeff + 12.61)
       endif

       if(eps < 0.44) then
          beta_param = (abs(eps-0.44))**0.7 * cos(pi*0.7)
       else
          beta_param = (eps - 0.44)**0.7
       endif

       delta_param = 0.55 * zeff**0.2 * (tanh(3.2*beta_param &
            *(eps**1.5 * nue_star_S)**1.4 / zeff**alpha_param) &
            + (1.0 - exp(-nue_star_S/0.1)) &
            * tanh(2.2*beta_param &
            * eps**2.8 * nue_star_S**0.1 / zeff**alpha_param))


        if(profile_model < 2 .or. h_flag == 2) then
           h_param = 1.0
        else
           pfac = 1.0/(rho(ir)*mass(is_ele)*vth(is_ele,ir)*sqrt(2.0)&
                *rmaj(ir)*sqrt(eps)) &
                * psiN_polflux_a * (1.0 - psiN_polflux(ir)) &
                /(a_meters**2 * b_unit(ir))
           if(h_flag == 1) then
              ! single null
              h_param = 1.0 - (0.2/zeff**4) &
                   * exp(-abs(pfac/(2.7*log(eps**1.5 *nue_star_S/3.2 + 3.0))))
           else
              ! double null
              h_param = 1.0 - (0.6/zeff**4) &
                   * exp(-abs(pfac/(3.3*log(eps**1.5 *nue_star_S/3.2 + 2.0))))
           endif
        endif

       ftrap_new = ftrap * h_param


       X31    = ftrap_new * (1.0 + delta_param) &
            / (1.0 + (1.0-0.1*ftrap_new) * sqrt(nue_star_S)  &
            + 0.5*(1.0-ftrap_new) * nue_star_S / (1.0*zeff))
       
       L31_S  = (1.0 + 1.4/(zeff+1)) * X31 &
            - (1.9/(zeff+1)) * X31*X31 &
            + (0.3/(zeff+1)) * X31*X31*X31 &
            + (0.2/(zeff+1)) * X31*X31*X31*X31
       
       X32e   = ftrap_new * (1.0 + delta_param) & 
            / (1.0 + 0.26*(1-ftrap_new) * sqrt(nue_star_S) &
            + 0.18*(1.0-0.37*ftrap_new) * nue_star_S / sqrt(1.0*zeff))
       
       F32_ee = (0.05 + 0.62*zeff) &
            / (zeff*(1+0.44*zeff)) * (X32e - X32e*X32e*X32e*X32e) &
            + 1.0/(1+0.22*zeff) * (X32e*X32e - X32e*X32e*X32e*X32e &
            - 1.2*(X32e*X32e*X32e - X32e*X32e*X32e*X32e)) &
            + 1.2/(1+0.5*zeff) * X32e*X32e*X32e*X32e
       
       X32i   = ftrap_new * (1.0 + delta_param) &
            / (1.0 + (1+0.6*ftrap_new) * sqrt(nue_star_S)  &
            + 0.85*(1-0.37*ftrap_new) * nue_star_S*(1.0+zeff))
       
       F32_ei = -(0.56 + 1.93*zeff) & 
            /(zeff*(1+0.44*zeff)) * (X32i - X32i*X32i*X32i*X32i)  &
            + 4.95/(1+2.48*zeff) * (X32i*X32i - X32i*X32i*X32i*X32i &
            - 0.55*(X32i*X32i*X32i - X32i*X32i*X32i*X32i)) &
            + (-1.2)/(1+0.5*zeff) * X32i*X32i*X32i*X32i;
       
       L32_S  = F32_ee + F32_ei
       
       X34   = ftrap_new * (1.0 + delta_param) &
            / (1.0 + (1.0-0.1*ftrap_new) * sqrt(nue_star_S) &
            + 0.5*(1.0-0.5*ftrap_new) * nue_star_S/(1.0*zeff))
       
       L34_S = (1 + 1.4/(zeff+1.0)) * X34 &
            - (1.9/(zeff+1.0)) * X34*X34 &
            + (0.3/(zeff+1.0)) * X34*X34*X34 &
            + (0.2/(zeff+1.0)) * X34*X34*X34*X34
       
       !!!!!!!!!!!!
       ! Epar components unchanged from Sauter
       ! but EAB modified X33 to match other X_ modifications
       
       X33 = ftrap_new * (1.0 + delta_param) & 
            / (1.0 + (0.55-0.1*ftrap_new) * sqrt(nue_star_S) &
            + 0.45*(1.0-ftrap_new) * nue_star_S/(1.0*zeff**1.5))


       sigma_spitzer = dens_ele / ( mass(is_ele) * nue_S ) &
            * 0.58 * 32 / (3.0*pi) / (0.58 + 0.74/(0.76+zeff))

       sigma_S = sigma_spitzer &
            * (1.0 - (1.0+0.36/(1.0*zeff)) * X33 &
            + 0.59/(1.0*zeff) * X33 * X33 &
            - 0.23/(1.0*zeff) * X33 * X33 * X33)
       !!!!!!!!!!!!

       jpar = sigma_S * EparB_avg &
            + L32_S * I_div_psip * rho(ir) * &
            dens(is_ele,ir) * temp(is_ele,ir) * dlntdr(is_ele,ir)

       do is=1,n_species
          jpar = jpar + L31_S * I_div_psip * rho(ir) &
               * dens(is,ir) * temp(is,ir) &
               * (dlntdr(is,ir) + dlnndr(is,ir))
          !if(is /= is_ele) then
          !   jpar = jpar + L34_S * alpha_S * I_div_psip * rho(ir) &
          !        * dlntdr(is,ir) * dens(is,ir) * temp(is,ir)
          !endif
       enddo

       ! EAB: 08/15/13 changed interpretation for ions with unlike temps
       ! from the above loop to this
       jpar = jpar + L34_S * alpha_S * I_div_psip * rho(ir) &
            * dlntdr(is_ion,ir) * press_sum

    end if
    
  end subroutine compute_Koh

end module neo_theory
