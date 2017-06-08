module neo_nclass_dr
  
  ! Driver for NCLASS calculation
  
  implicit none
  
  public :: NCLASS_DR_alloc, NCLASS_DR_do
  
  integer, parameter, private :: mx_mi = 9
  integer, parameter, private :: mx_ms = 40
  integer, parameter, private :: mx_mz = 100 ! increased this and pamx_mz.inc
					     ! for higher Z
  
  integer, parameter, private :: io_nc = 41
  character(len=80),private :: runfile = 'out.neo.theory_nclass'
  logical, private :: initialized = .false.
  real, dimension(:), allocatable :: pflux_nc, eflux_nc
  real, dimension(:), allocatable :: uparB_nc, vpol_nc, vtor_nc
  real :: jbs_nc
  real :: cc, loglam
  integer :: z_indx ! NOTE: NCLASS Z's are NEO Z's cast to integers
  
contains
  
  subroutine NCLASS_DR_alloc(flag)
    use neo_globals
    implicit none
    integer, intent (in) :: flag  ! flag=1: allocate; else deallocate

    if (n_species < 2) then
       ! NCLASS driver requires at least 2 species
       return
    endif

    if(flag == 1) then
       if(initialized) return
       allocate(pflux_nc(n_species))
       allocate(eflux_nc(n_species))
       allocate(uparB_nc(n_species))
       allocate(vpol_nc(n_species))
       allocate(vtor_nc(n_species))
       if(silent_flag == 0 .and. i_proc == 0) then
          open(io_nc,file=trim(path)//runfile,status='replace')
          close(io_nc)
       end if
       initialized = .true.

    else
       if(.NOT. initialized) return
       deallocate(pflux_nc)
       deallocate(eflux_nc)
       deallocate(uparB_nc)
       deallocate(vpol_nc)
       deallocate(vtor_nc)
       initialized = .false.
    end if


  end subroutine NCLASS_DR_alloc

  subroutine NCLASS_DR_do(ir)
    use neo_globals
    use neo_equilibrium
    implicit none
    integer, intent (in) :: ir
    real, parameter :: potato_flag=0
    real, parameter :: squeeze_flag=0
    real, parameter :: z_protonmass=1.6726
    real, parameter :: z_coulomb=1.6022e-19
    real, parameter :: z_j7kv=1.6022e-16
    integer :: is, it, i, j, k, im, iz, iza
    real    :: ppr, uthai
    real, dimension(8) :: rdum, dum, edum
    ! NCLASS input
    integer :: k_order, k_potato
    real    :: c_potb, c_potl
    real    :: p_grphi, p_gr2phi, p_eb
    real    :: p_eps, p_q, p_b2, p_bm2, p_grbm2, p_ngrth, p_ft, p_fhat
    real, dimension(3) :: p_fm
    real    :: btor0_nc, bpol0_nc
    integer :: m_i, m_z
    real    :: c_den
    real, dimension(3,mx_mi,mx_mz) :: fex_iz
    real, dimension(mx_mi) :: amu_i, temp_i, grt_i
    real, dimension(mx_mi,mx_mz) :: den_iz, grp_iz
    ! conversions
    real    :: pflux_nc_norm, eflux_nc_norm, jbs_nc_norm, v_nc_norm
    real    :: sum1, sum2, sum3, sum4, fac1
    ! NCLASS output
    integer :: iflag, m_s
    real    :: p_bsjb, p_etap, p_exjb
    integer, dimension(mx_ms)    :: jm_s, jz_s
    real, dimension(3,3,mx_mi)   ::  calm_i
    real, dimension(3,3,mx_mi, mx_mi) :: caln_ii,capm_ii, capn_ii
    real, dimension(mx_ms)       :: bsjbp_s, bsjbt_s, dn_s, sqz_s
    real, dimension(5,mx_ms)     :: gfl_s,qfl_s
    real, dimension(3,3,mx_ms)   :: upar_s, utheta_s, ymu_s
    real, dimension(mx_ms)       :: vn_s, veb_s, qeb_s, xi_s
    real, dimension(mx_ms,mx_ms) :: chip_ss, chit_ss, dp_ss, dt_ss
    ! declaration of functions
    real  RARRAY_SUM

    if (n_species < 2) then
       ! NCLASS driver requires at least 2 species
       if(silent_flag==0 .and. i_proc==0) then
          open(unit=io_neoout,file=trim(path)//runfile_neoout,&
               status='old',position='append')
          write(io_neoout,*) 'NCLASS not computed: Requires at least 2 species'
          close(io_neoout)
       endif
       return
    endif

    if (profile_model /= 2) then
       ! For local mode, need to set normalizations
       ! Assume normalizing mass is deuterium
       ! Set T_norm = 1kEV and a_meters=1m
       ! Then determine B_unit from input rho_star
       ! Dens in units and norm will be determined from input collision freqs
       temp_norm_fac   = 1.6022*1000
       charge_norm_fac = 1.6022
       a_meters        = 1.0
       temp_norm(ir)   = 1.0
       vth_norm(ir)    = sqrt(temp_norm(ir) * temp_norm_fac &
            / (mass_deuterium)) &
            * 1.0e4 / a_meters
       b_unit(ir)      =  sqrt(temp_norm(ir) * temp_norm_fac &
            * mass_deuterium) &
            / (charge_norm_fac * rho(ir)) &
            * 1.0e-4 / a_meters
       do is=2,n_species
          if(abs(nu(is,ir)-nu(1,ir)*(1.0*z(is))**4/(1.0*z(1))**4 &
               * dens(is,ir) / dens(1,ir) &
               * sqrt(mass(1)/mass(is)) * (temp(1,ir)/temp(is,ir))**1.5) &
               > 1e-6) then
             if(silent_flag==0 .and. i_proc==0) then
                open(unit=io_neoout,file=trim(path)//runfile_neoout,&
                     status='old',position='append')
                write(io_neoout,*) 'Warning: NCLASS requires self-consistent collision frequencies'
                close(io_neoout)
                exit
             endif
          endif
       enddo
    endif
    
    !  k_order-order of v moments to be solved
    !         = 2 u and p_q
    !         = 3 u, p_q, and u2
    !         = else error
    k_order=3
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Potato-orbit parameters
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  k_potato-option to include potato orbits (-)
    !          =0 off
    !          =else on
    if(potato_flag == 1) then
       k_potato=1
       !  c_potb-kappa(0)*Bt(0)/[2*q(0)**2] (T)
       c_potb = kappa(ir) * Btor_th0 * b_unit(ir) / (2*q(ir)**2)
       !  c_potl-q(0)*R(0) (m)
       c_potl = q(ir) * (rmaj(ir) + r(ir)) * a_meters
    else
       k_potato=0
       c_potb = 0.0
       c_potl = 0.0
    end if
       

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Er0
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  p_grphi- dPhi/drho (V/rho)
    ! This gets added to the pressure gradient, grp_iz, in computing the flow
    p_grphi = dphi0dr(ir) * (1000 * temp_norm(ir)) / a_meters
    !  p_gr2phi-radial electric field gradient Psi'(Phi'/Psi')' (V/rho**2)
    if(squeeze_flag == 1) then
       p_gr2phi = -r(ir) / q(ir) * b_unit(ir) &
            * omega_rot_deriv(ir) * vth_norm(ir)
    else
       p_gr2phi = 0.0
    end if

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Epar
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  p_eb-<E.B> (V*T/m)
    p_eb = epar0(ir) * (1000 * temp_norm(ir)) / a_meters * b_unit(ir)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Geometry
    ! rho -> r
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! inverse aspect ratio
    p_eps = r(ir) / rmaj(ir)
    
    ! safety factor
    p_q = q(ir)
    
    !  p_ft-trapped fraction
    p_ft = ftrap
    
    !  p_fhat-mu_0*F/(dPsi/dr) (rho/m)
    p_fhat = I_div_psip
    
    ! B at outboard midplane (T)
    btor0_nc = Btor_th0 * b_unit(ir)
    bpol0_nc = Bpol_th0 * b_unit(ir)
    
    !  p_b2-<B**2> (T**2)
    p_b2 = Bmag2_avg * b_unit(ir)**2
    
    !  p_bm2-<1/B**2> (/T**2)
    p_bm2 = Bmag2inv_avg / b_unit(ir)**2
    
    !  p_grbm2-<grad(rho)**2/B**2> (rho**2/m**2/T**2)
    p_grbm2 = 0.0
    do it=1,n_theta
       p_grbm2 = p_grbm2 + w_theta(it) * gradr(it)**2 &
            / (Bmag(it) * b_unit(ir))**2
    enddo
    
    !  p_ngrth-<n.grad(Theta)> (1/m)
    !  This refers to the nclass theta
    !  inverse of theta-avg of 1/(bhat dot grad theta)
    p_ngrth = 0.0
    do it=1,n_theta
       p_ngrth = p_ngrth + 1.0/(2.0*pi)  * 1.0 / (k_par(it) / a_meters) &
            * 2.0 * pi / n_theta
    enddo
    p_ngrth = 1.0 / p_ngrth

    !  p_fm(3)-poloidal moments of geometric factor for PS viscosity (-)
    ! first get <B>
    fac1 = 0.0
    do it=1,n_theta
       fac1 = fac1 + w_theta(it) * Bmag(it) * b_unit(ir)
    enddo
    do i=1,3
       sum1 = 0.0
       sum2 = 0.0
       sum3 = 0.0
       sum4 = 0.0
       do it=1,n_theta
          sum1 = sum1 + w_theta(it) * sin(i*theta_nc(it)) &
               * gradpar_Bmag(it) * b_unit(ir) / a_meters
          sum2 = sum2 + w_theta(it) * sin(i*theta_nc(it)) &
               * gradpar_Bmag(it) * b_unit(ir) / a_meters &
               * Bmag(it) * b_unit(ir) * p_ngrth
          sum3 = sum3 + w_theta(it) * cos(i*theta_nc(it)) &
               * gradpar_Bmag(it) * b_unit(ir) / a_meters
          sum4 = sum4 + w_theta(it) * cos(i*theta_nc(it)) &
               * gradpar_Bmag(it) * b_unit(ir) / a_meters &
               * Bmag(it) * b_unit(ir) * p_ngrth
       enddo
       p_fm(i) = 2.0 /  ( p_b2 * p_ngrth * fac1) &
            * (sum1 * sum2 + sum3 * sum4)
     enddo

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Species-dependent paramters
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !  m_i-number of isotopes (1<mi<mx_mi+1)
    m_i = n_species
    
    !  m_z-highest charge state of all species (0<mz<mx_mz+1)
    m_z = 1
    do is=1, n_species
       if(z(is) > m_z) then
          m_z = int(z(is))
       endif
    enddo

    !  amu_i(i)-atomic mass number of i
    !  EAB: Note that this presently assumes mass is relative to deuterium
    amu_i(:) = 0.0
    do is=1,n_species
       amu_i(is) = mass(is) * mass_deuterium / z_protonmass
    enddo
    
    !  temp_i(i)-temperature of i (keV)
    temp_i(:) = 0.0
    do is=1,n_species
       temp_i(is) = temp(is,ir) * temp_norm(ir)
    enddo
    
    !  grt_i(i)-temperature gradient of i (keV/rho)
    grt_i(:) = 0.0 
    do is=1,n_species
       grt_i(is) = -dlntdr(is,ir) * temp_i(is) / a_meters
    enddo

    !  den_iz(i,z)-density of i,z (/m**3) 
    den_iz(:,:) = 0.0
    if(profile_model /=2) then
       ! For local mode, determine dens in units from collision freqs
       cc = sqrt(2.0) * pi * charge_norm_fac**4 &
            * 1.0 / (4.0 * pi * 8.8542)**2 &
            * 1.0 / (sqrt(mass_deuterium) * temp_norm_fac**1.5) &
            * 1e9
       ! use an approx for loglam
       loglam = 24.0 - log(sqrt(1.0*1e13)/(1.0*1000))
       is=1
       z_indx = int(abs(z(is)))
       den_iz(is,z_indx) = (nu(is,ir)*vth_norm(ir)) &
            / (cc * loglam * z(is)**4) &
            * (sqrt(mass(is)) * (temp(is,ir)*temp_norm(ir))**1.5) * 1e19
       dens_norm(ir) = den_iz(is,z_indx) * 1e-19 / dens(is,ir)
       do is=2,n_species
          z_indx = int(abs(z(is)))
          den_iz(is,z_indx) = dens(is,ir)  * dens_norm(ir) * 1e19
       enddo
    else
       do is=1,n_species
          z_indx = int(abs(z(is)))
          den_iz(is,z_indx) = dens(is,ir) * dens_norm(ir) * 1e19
       enddo
    endif

    !  c_den-density cutoff below which species is ignored (/m**3)
    c_den=1.0e10
    do is=1,n_species
       z_indx = int(abs(z(is)))
       if(den_iz(is,z_indx) < c_den) then
          c_den = den_iz(is,z_indx) 
       endif
    enddo
    c_den = c_den - 1.0

    !  grp_iz(i,z)-pressure gradient of i,z (keV/m**3/rho)
    grp_iz(:,:) = 0.0
    do is=1,n_species
       z_indx = int(abs(z(is)))
       grp_iz(is,z_indx) = -den_iz(is,z_indx) *  temp_i(is) &
            * (dlnndr(is,ir) + dlntdr(is,ir)) / a_meters
    enddo

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Sources
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! fex_iz(3,i,z)-moments of external parallel force on i,z (T*j/m**3) 
    !do is=1,n_species
    !   z_indx = int(abs(z(is)))
    !   fex_iz(1,is,z_indx) = source_nclass(is) &
    !        * dens_norm(ir) * 1e19 * b_unit(ir) &
    !        * (vth_norm(ir) * a_meters)**2 * (mass_deuterium*1e-27)
    !enddo
    fex_iz(:,:,:) = 0.0

    ! Normalizations
    pflux_nc_norm = dens_norm(ir) * vth_norm(ir) * a_meters * 1e19
    eflux_nc_norm = dens_norm(ir)*vth_norm(ir)*a_meters &
         * temp_norm(ir)*temp_norm_fac
    jbs_nc_norm   = charge_norm_fac*dens_norm(ir)*vth_norm(ir)*a_meters
    v_nc_norm   = vth_norm(ir) * a_meters

    ! compute friction and viscosity coefficients and flows
    CALL NCLASS(k_order,k_potato,m_i,m_z,c_den,c_potb,c_potl,p_b2, &
         p_bm2,p_eb,p_fhat,p_fm,p_ft,p_grbm2,p_grphi,p_gr2phi, &
         p_ngrth,amu_i,grt_i,temp_i,den_iz,fex_iz,grp_iz,m_s, &
         jm_s,jz_s,p_bsjb,p_etap,p_exjb,calm_i,caln_ii,capm_ii, &
         capn_ii,bsjbp_s,bsjbt_s,dn_s,gfl_s,qfl_s,sqz_s,upar_s, &
         utheta_s,vn_s,veb_s,qeb_s,xi_s,ymu_s,chip_ss,chit_ss, &
         dp_ss,dt_ss,iflag)

    ! check warning flags
    if(silent_flag == 0 .and. i_proc == 0) then
       open(unit=io_neoout,file=trim(path)//runfile_neoout,&
            status='old',position='append')
       if(iflag == -1) then
          ! write(io_neoout,*)  'WARNING: (NEO) NCLASS - no potato orbit viscosity'
       else if(iflag == -2) then
          write(io_neoout,*) 'WARNING: (NEO) NCLASS - no Pfirsch-Schluter viscosity'
       else if(iflag == -3) then
          write(io_neoout,*)  'WARNING: (NEO) NCLASS - no banana viscosity'
       else if(iflag == -4) then
          write(io_neoout,*)  'WARNING: (NEO) NCLASS - no viscosity'
       end if
       close(io_neoout)
    end if

    ! check error flags
    if(iflag > 0) then
       if(silent_flag == 0 .and. i_proc == 0) then
          open(unit=io_neoout,file=trim(path)//runfile_neoout,&
               status='old',position='append')
          if(iflag == 1) then
             write(io_neoout,*) 'WARNING: (NEO) NCLASS - k_order must be 2 or 3'
          elseif(iflag == 2) then
             write(io_neoout,*) 'WARNING: (NEO) NCLASS - require 1<m_i<mx_mi'
          else if(iflag == 3) then
             write(io_neoout,*) 'WARNING: (NEO) NCLASS - require 0<m_z<mx_mz'
          else if(iflag == 4) then
             write(io_neoout,*) 'WARNING: (NEO) NCLASS - require 0<m_s<mx_ms'
          else if(iflag == 5) then
             write(io_neoout,*) 'WARNING: (NEO) NCLASS - inversion of flow matrix failed'
          else if(iflag == 6) then
             write(io_neoout,*) 'WARNING: (NEO) NCLASS - invalid trapped fraction'
          endif
          close(io_neoout)
       endif
       return
    endif

    !  Parallel flows: upar_ij (T*m/s) for each species
    do i=1,m_s
       do j=1,k_order
          do k=1,k_order
             rdum(k)=upar_s(j,k,i)
          enddo
          if(j==1) then
             uparB_nc(i) = (rdum(1)+rdum(2)+rdum(3)) / (v_nc_norm * b_unit(ir))
          endif
       enddo
    enddo
    
    ! Bootstrap current: <J_bs.B>/Bunit (A/m**2)
    rdum(1)=p_bsjb/b_unit(ir)
    if(abs(p_eb) > 0.0) then
       rdum(2) = (p_eb/p_etap)/b_unit(ir)
    else
       rdum(2) = 0.0
    endif
    if(abs(p_exjb) > 0.0) then
       rdum(3) = p_exjb/b_unit(ir)
    else
       rdum(3) = 0.0
    endif
    jbs_nc = (rdum(1)+rdum(2)+rdum(3)) / jbs_nc_norm

    !  Flow velocities on outside midplane (m/s)                      
    do i=1,m_s
       im=jm_s(i)
       iz=jz_s(i)
       iza=IABS(iz)
       ppr=p_fhat*grp_iz(im,iza)*z_j7kv &
            /(z_coulomb*iz*den_iz(im,iza))+p_fhat*p_grphi
       uthai=utheta_s(1,1,i)+utheta_s(1,2,i)+utheta_s(1,3,i)
       ! Toroidal
       rdum(1)=uthai*btor0_nc-ppr/btor0_nc
       vtor_nc(i) = rdum(1) / v_nc_norm
       ! Poloidal
       rdum(2)=uthai*bpol0_nc
       vpol_nc(i) = rdum(2) / v_nc_norm
    enddo

    ! Radial particle fluxes (/m**2/s)
    ! BP, PS, CL, <E.B>, src, total
    ! Load into rdum the five flux components and total
    ! Load into dum the z-weighted + charge (ion) components
    ! Load into edum the z-weighted - charge (electron) components 
    call RARRAY_ZERO(6,dum)
    call RARRAY_ZERO(6,edum)
    do i=1,m_s
       iz=jz_s(i)
       call RARRAY_COPY(5,gfl_s(1,i),1,rdum,1)
       rdum(6)=RARRAY_SUM(5,rdum,1)
       pflux_nc(i) = (rdum(1) + rdum(2) + rdum(4) + rdum(5)) / pflux_nc_norm
       do k=1,5
          if(iz > 0) then
             dum(k)=dum(k)+iz*gfl_s(k,i)
          else
             edum(k)=edum(k)+iz*gfl_s(k,i)
          end if
       enddo
       if(iz > 0) THEN
          dum(6)=dum(6)+iz*rdum(6)
       else
          edum(6)=edum(6)+iz*rdum(6)
       end if
    enddo

    ! Ambipolarity check
    do k=1,6
       rdum(k)=edum(k)+dum(k)
    enddo

    !  Radial conduction fluxes (W/m**2)
    ! BP, PS, CL, <E.B>, src, total
    call RARRAY_ZERO(6,dum)
    call RARRAY_ZERO(6,edum)
    do i=1,m_s
       iz=jz_s(i)
       call RARRAY_COPY(5,qfl_s(1,i),1,rdum,1)
       rdum(6)=RARRAY_SUM(5,qfl_s(1,i),1)         
       do k=1,5
          if(iz > 0) then
             dum(k)=dum(k)+qfl_s(k,i)
          else
             edum(k)=edum(k)+qfl_s(k,i)
          end if
       end do
       if(iz > 0) then
          dum(6)=dum(6)+rdum(6)
       else
          edum(6)=edum(6)+rdum(6)
       end if
    end do
    call RARRAY_ZERO(6,dum)
    do i=1,m_s
       im=jm_s(i)
       iz=jz_s(i)
       do k=1,5
          rdum(k)=qfl_s(k,i)+2.5*gfl_s(k,i)*temp_i(im)*z_j7kv
       end do
       rdum(6)=RARRAY_SUM(5,rdum,1)
       eflux_nc(i) = (rdum(1) + rdum(2) + rdum(4) + rdum(5)) / eflux_nc_norm
    enddo

    if(silent_flag == 0 .and. i_proc == 0) then
       open(io_nc,file=trim(path)//runfile,status='old',position='append')
       write (io_nc,'(e16.8)',advance='no') r(ir)
       write(io_nc,'(e16.8)',advance='no') jbs_nc 
       do is=1,n_species
          write(io_nc,'(e16.8)',advance='no') pflux_nc(is)
          write(io_nc,'(e16.8)',advance='no') eflux_nc(is) 
          write(io_nc,'(e16.8)',advance='no') uparB_nc(is) 
          write(io_nc,'(e16.8)',advance='no') vpol_nc(is)
          write(io_nc,'(e16.8)',advance='no') vtor_nc(is)
       enddo
       write(io_nc,*)
       close(io_nc)
    end if

    ! Store the local nclass transport values at ir=1 in neo_nclass_out
    if(ir == 1) then
       neo_nclass_out(:,:) = 0.0
       neo_nclass_1d_out   = 0.0
       do is=1,n_species
          neo_nclass_out(is,1) = pflux_nc(is)
          neo_nclass_out(is,2) = eflux_nc(is)
          neo_nclass_out(is,3) = vpol_nc(is)
          neo_nclass_out(is,4) = vtor_nc(is)
          neo_nclass_out(is,5) = ymu_s(1,1,is) / (mass_deuterium * 1e-27 &
            * dens_norm(ir) * 1e19 * vth_norm(ir))
       enddo
       neo_nclass_1d_out = jbs_nc 
    endif
    
  end subroutine NCLASS_DR_do

end module neo_nclass_dr
