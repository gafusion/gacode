subroutine cgyro_init_collision

  use timer_lib

  use cgyro_globals

  implicit none

  real, dimension(:,:,:), allocatable :: nu_d, nu_par
  real, dimension(:,:), allocatable :: rs
  real, dimension(:,:,:,:), allocatable :: rsvec, rsvect0, rsvect1
  real, dimension(:,:), allocatable :: klor_fac, kdiff_fac

  real :: arg
  real :: xa, xb, tauinv_ab
  real :: rval
  integer :: jv
  integer :: is,ir,it,ix,ie,js,je,jx,ks
  integer :: dv
  ! parameters for matrix solve
  real, dimension(:,:), allocatable :: amat,cmat_loc
  real, dimension(:,:,:,:,:,:), allocatable :: ctest
  real, dimension(:,:,:,:,:), allocatable :: bessel
  
  if (collision_model == 5) then
     call cgyro_init_collision_simple
     return
  endif

  allocate(nu_d(n_energy,n_species,n_species))
  allocate(nu_par(n_energy,n_species,n_species))
  allocate(klor_fac(n_species,n_species))
  allocate(kdiff_fac(n_species,n_species))
  nu_d(:,:,:) = 0.0
  nu_par(:,:,:) = 0.0
  klor_fac(:,:) = 0.0
  kdiff_fac(:,:) = 0.0
  
  do ie=1,n_energy
     do is=1,n_species
        do js=1,n_species

           xa = vel(ie)
           xb = xa * vth(is) / vth(js)
           tauinv_ab = nu(is) * z(js)**2 / z(is)**2 &
                * dens(js)/dens(is)
           ! re-scale only electron collisions (ee,ei,ie)
           if(is == is_ele .or. js == is_ele) then
              tauinv_ab = tauinv_ab * collision_ele_scale
           endif
              
           select case (collision_model)

           case (1)

              ! Only ee,ei Connor-like Lorentz
              if (is == is_ele) then
                 if (is == js) then
                    ! e-e
                    nu_d(ie,is,js) = tauinv_ab * (1.0/xa**3) &
                         * (exp(-xb*xb)/(xb*sqrt(pi)) &
                         + (1.0-1.0/(2.0*xb*xb)) * erf(xb))
                 else
                    ! e-i
                    if(z_eff_method == 2) then
                       nu_d(ie,is,js) = tauinv_ab * (1.0/xa**3)
                    else
                       nu_d(ie,is,js) = tauinv_ab * (1.0/xa**3) &
                            * z(is)**2 / z(js)**2 &
                            * dens(is)/dens(js) * z_eff/(n_species-1)
                    endif
                 endif
              endif

           case (2)

              ! Connor model
              if (is == js .or. &
                   (abs(mass(is) - mass(js)) < epsilon(0.0))) then
                 ! case 1: like-species/same mass collisions
                 nu_d(ie,is,js) = tauinv_ab * (1.0/xa**3) &
                      * (exp(-xb*xb)/(xb*sqrt(pi)) &
                      + (1.0-1.0/(2.0*xb*xb)) * erf(xb))

              else if (mass(is) < mass(js)) then
                 ! case 2: ele-ion and ion-imp(heavy) collisions
                 nu_d(ie,is,js) = tauinv_ab * (1.0/xa**3)

              else
                 ! case 3: ion-ele and imp(heavy)-ion collisions
                 nu_d(ie,is,js) = tauinv_ab * 4.0/(3.0*sqrt(pi)) &
                      * sqrt(mass(js)/mass(is)) * (temp(is)/temp(js))**1.5
              endif

           case(4)

              ! Ad hoc op
              ! (Fix for underflow)
              nu_d(ie,is,js) = tauinv_ab * (1.0/xa**3) &
                   * (exp(-xb*xb)/(xb*sqrt(pi)) &
                   + (1.0-1.0/(2.0*xb*xb)) * erf(xb))
              ! No i-e Lorentz
              !if (is /= is_ele .and. js == is_ele) then
              !   nu_d(ie,is,js) = 0.0
              !endif
              if(collision_kperp == 1) then
                 klor_fac(is,js) = 1.0
              endif

              ! Diffusion 
              nu_par(ie,is,js) = tauinv_ab * (2.0/xa**3) &
                   * (-exp(-xb*xb)/(xb*sqrt(pi)) &
                   + (1.0/(2.0*xb*xb)) * erf(xb))
              if(collision_kperp == 1) then
                 kdiff_fac(is,js) = 1.0
              endif
              
           end select

        enddo
     enddo
  enddo

  if (collision_ion_model == 1) then
     do is=1,n_species
        if(is /= is_ele) then
           do js=1,n_species
              nu_d(:,is,js) = 0.0
              nu_par(:,is,js) = 0.0
           enddo
        endif
     enddo
  endif
  
  ! Printout used for CGYRO paper.
  !if (i_proc == 0) then
  !   do ie=1,n_energy
  !      print '(i1,3(" &{\tt ",f9.5,"}")," \\")',&
  !           ie,vel(ie),w_e(ie),nu_d(ie,n_species,n_species)
  !   enddo
  !   print *,sum(w_e)
  !endif

  if ( collision_model == 4 .and. collision_kperp == 1 .and. &
       (collision_mom_restore == 1 .or. collision_ene_restore == 1)) then
     allocate(bessel(n_species,n_xi,n_energy,nc_loc,0:1))
!$omp parallel do private(ic_loc,it,ie,ix,is,arg)
     do ic=nc1,nc2
        ic_loc = ic-nc1+1
        it = it_c(ic)
        do ie=1,n_energy
           do ix=1,n_xi
              do is=1,n_species   
                 arg = k_perp(ic)*rho*vth(is)*mass(is)&
                      /(z(is)*bmag(it)) *sqrt(2.0*energy(ie)) &
                      *sqrt(1.0-xi(ix)**2)
                 bessel(is,ix,ie,ic_loc,0) = bessel_j0(arg)
                 bessel(is,ix,ie,ic_loc,1) = bessel_j1(arg)
              enddo
           enddo
        enddo
     enddo
  endif

  allocate(ctest(n_species,n_species,n_xi,n_xi,n_energy,n_energy))
  allocate(rs(n_species,n_species))
  allocate(rsvec(n_species,n_species,n_xi,n_energy))
  allocate(rsvect0(n_species,n_species,n_xi,n_energy))
  allocate(rsvect1(n_species,n_species,n_xi,n_energy))

  allocate(amat(nv,nv))
  allocate(cmat_loc(nv,nv))

  ! Collision test particle component
  ctest = 0.0

  ! Lorentz
  do is=1,n_species
     do ix=1,n_xi
        do ie=1,n_energy
           do js=1,n_species
              do jx=1,n_xi
                 je = ie
                 ctest(is,js,ix,jx,ie,je) &
                      = ctest(is,js,ix,jx,ie,je) &
                      + xi_lor_mat(ix,jx) *0.5*nu_d(ie,is,js)
              enddo
           enddo
        enddo
     enddo
  enddo

  ! Diffusion
  if (collision_model == 4 .and. collision_ene_diffusion == 1) then
!$omp parallel do collapse(5) private(is,ix,ie,js,je,jx)     
     do is=1,n_species 
        do ix=1,n_xi
           do ie=1,n_energy
              do js=1,n_species
                 do je=1,n_energy
                    jx = ix
                    ! From K. Hallatschek
                    ! self-adjoint part of ctest written self-adjointly
                    ctest(is,js,ix,jx,ie,je) = ctest(is,js,ix,jx,ie,je) &
                         -0.5 / w_e(ie) &
                         *sum(w_e(:)*e_deriv1_mat(:,ie)*energy(:) &
                         *nu_par(:,is,js) *e_deriv1_mat(:,je))/(1.0*e_max)
                    ! non-self-adjoint part proportional 1-Ta/Tb written
                    ! in a way that supports inherent particle number 
                    ! conservation for small kperp
                    ctest(is,js,ix,jx,ie,je) = ctest(is,js,ix,jx,ie,je) &
                         + (1-temp(is)/temp(js)) / sqrt(1.0*e_max)/w_e(ie) &
                         * w_e(je)*e_deriv1_mat(je,ie) &
                         * nu_par(je,is,js)*energy(je)**1.5
                 enddo
              enddo
           enddo
        enddo
     enddo
  endif

  ! matrix solve parameters
  allocate(i_piv(nv))

  ! Construct the collision matrix

!$omp  parallel do  default(none) &
!$omp& shared(nc1,nc2,nv,my_toroidal,delta_t,n_species,rho,is_ele,n_field,n_energy,n_xi) &
!$omp& shared(collision_kperp,collision_field_model,explicit_trap_flag) &
!$omp& firstprivate(collision_model,collision_mom_restore,collision_ene_restore) &
!$omp& shared(ae_flag,lambda_debye,dens_ele,temp_ele,dens_rot) &
!$omp& shared(betae_unit,sum_den_h) &
!$omp& shared(it_c,ir_c,px,is_v,ix_v,ie_v,ctest,xi_deriv_mat) &
!$omp& shared(temp,jvec_v,omega_trap,dens,energy,vel) &
!$omp& shared(omega_rot_trap,omega_rot_u,e_deriv1_mat,e_deriv1_rot_mat,e_max,bessel) &
!$omp& shared(xi_lor_mat) &
!$omp& shared(k_perp,vth,mass,z,bmag,nu_d,xi,nu_par,w_e,w_xi) &
!$omp& shared(klor_fac,kdiff_fac) &
!$omp& private(ic,ic_loc,it,ir,info,rval) &
!$omp& private(iv,is,ix,ie,jv,js,jx,je,ks) &
!$omp& private(amat,cmat_loc,i_piv,rs,rsvec,rsvect0,rsvect1) &
!$omp& private(dv) firstprivate(collision_precision_mode, collision_full_stripes) &
!$omp& shared(cmat,cmat_fp32,cmat_stripes)
  do ic=nc1,nc2
   
     ic_loc = ic-nc1+1

     it = it_c(ic)
     ir = ir_c(ic)

     ! Collision field particle component
     amat(:,:)   = 0.0
     cmat_loc(:,:) = 0.0

     select case (collision_model)

     case(2)
     if (collision_mom_restore == 1) then
        do is=1,n_species
           do js=1,n_species
              rs(is,js) = 0.0
              do ie=1,n_energy
                 rs(is,js) = rs(is,js) + w_e(ie)*nu_d(ie,is,js)*energy(ie)
              enddo
           enddo
        enddo

           do iv=1,nv  
              is = is_v(iv)
              ix = ix_v(iv)
              ie = ie_v(iv)

              do jv=1,nv
                 js = is_v(jv)
                 jx = ix_v(jv)
                 je = ie_v(jv)

                 if (abs(rs(is,js)) > epsilon(0.0)) then
                    cmat_loc(iv,jv) = &
                         cmat_loc(iv,jv) &
                         + 3.0 * (mass(js)/mass(is)) &
                         * (dens(js)/dens(is)) * dens_rot(it,js) &
                         * (vth(js)/vth(is)) * nu_d(ie,is,js) &
                         * vel(ie) * xi(ix) &
                         * nu_d(je,js,is) * sqrt(energy(je)) &
                         * xi(jx) * w_e(je) * w_xi(jx) / rs(is,js)
                 endif
              enddo
           enddo

     endif

     case(4)

     ! Momentum Restoring

     if (collision_mom_restore == 1) then

        ! C_test_ab(v_par f0a,f0b) and w_e v_par C_test_ab(H_a)
        rsvec = 0.0
        rsvect0 = 0.0
        rs(:,:) = 0.0
        
        do is=1,n_species
           do js=1,n_species
              do ix=1,n_xi
                 do jx=1,n_xi
                    do ie=1,n_energy
                       do je=1,n_energy
                          rsvec(is,js,ix,ie) = rsvec(is,js,ix,ie) &
                               + ctest(is,js,ix,jx,ie,je) &
                               * sqrt(2.0*energy(je)) * xi(jx) * vth(is)
                          rsvect0(is,js,ix,ie) = rsvect0(is,js,ix,ie) &
                               + ctest(is,js,jx,ix,je,ie) &
                               * sqrt(2.0*energy(je)) * xi(jx) &
                               * w_e(je)*w_xi(jx) * vth(is)
                       enddo
                    enddo
                 enddo
              enddo

              ! int v_par C_test_ab(v_par f0a,f0b) / n_0a
              do ix=1,n_xi
                 do ie=1,n_energy
                    rs(is,js) = rs(is,js) + w_e(ie)*w_xi(ix) * dens(is) &
                         * rsvec(is,js,ix,ie) * sqrt(2.0*energy(ie)) * xi(ix) &
                         * vth(is)
                 enddo
              enddo
           enddo
        enddo

        if (collision_kperp == 0) then
              do iv=1,nv  
                 is = is_v(iv)
                 ix = ix_v(iv)
                 ie = ie_v(iv)

                 do jv=1,nv
                    js = is_v(jv)
                    jx = ix_v(jv)
                    je = ie_v(jv)

                    if (abs(rs(is,js))>epsilon(0.0)) then
                       cmat_loc(iv,jv) &
                            = cmat_loc(iv,jv) &
                            - mass(js)/mass(is) * dens(js) * dens_rot(it,js) &
                            * rsvec(is,js,ix,ie) / rs(is,js) &
                            * rsvect0(js,is,jx,je)
                    endif
                 enddo
              enddo

        else
              rsvect0(:,:,:,:) = 0.0
              rsvect1(:,:,:,:) = 0.0
              do is=1,n_species
                 do js=1, n_species
                    do ix=1,n_xi
                       do jx=1,n_xi
                          do ie=1,n_energy
                             do je=1,n_energy
                                rsvect0(is,js,ix,ie) = rsvect0(is,js,ix,ie) &
                                     + ctest(is,js,jx,ix,je,ie) &
                                     * sqrt(2.0*energy(je)) * xi(jx) &
                                     * w_e(je)*w_xi(jx) * vth(is) &
                                     * bessel(is,ix,ie,ic_loc,0)
                                rsvect1(is,js,ix,ie) = rsvect1(is,js,ix,ie) &
                                     + ctest(is,js,jx,ix,je,ie) &
                                     * sqrt(2.0*energy(je)) * xi(jx) &
                                     * w_e(je)*w_xi(jx) * vth(is) &
                                     * bessel(is,ix,ie,ic_loc,1) &
                                     * sqrt(1.0-xi(ix)**2)/xi(ix)
                             enddo
                          enddo
                       enddo
                    enddo
                 enddo
              enddo
              
              do iv=1,nv  
                 is = is_v(iv)
                 ix = ix_v(iv)
                 ie = ie_v(iv)

                 do jv=1,nv
                    js = is_v(jv)
                    jx = ix_v(jv)
                    je = ie_v(jv)

                    if (abs(rs(is,js)) > epsilon(0.)) then 
                       cmat_loc(iv,jv) &
                            = cmat_loc(iv,jv) &
                            - mass(js)/mass(is) &
                            * dens(js) * dens_rot(it,js) &
                            * rsvec(is,js,ix,ie) &
                            * bessel(is,ix,ie,ic_loc,0) / rs(is,js) &
                            * rsvect0(js,is,jx,je)
                       cmat_loc(iv,jv) &
                            = cmat_loc(iv,jv) &
                            - mass(js)/mass(is) &
                            * dens(js) * dens_rot(it,js) &
                            * rsvec(is,js,ix,ie) / rs(is,js) &
                            * bessel(is,ix,ie,ic_loc,1) &
                            * sqrt(1.0-xi(ix)**2)/xi(ix) &
                            * rsvect1(js,is,jx,je) 
                    endif
                 enddo
              enddo
        endif

     endif

     ! Energy Restoring

     if (collision_ene_restore == 1) then
              
        ! C_test_ab(u_a^2 f0a,f0b) and w_e u_a^2 C_test_ab(H_a)
        rsvec  = 0.0
        rsvect0 = 0.0
        rs(:,:) = 0.0

        do is=1,n_species
           do js=1,n_species
              do ix=1,n_xi
                 do jx=1,n_xi
                    do ie=1,n_energy
                       do je=1,n_energy
                          rsvec(is,js,ix,ie) = rsvec(is,js,ix,ie) &
                               + ctest(is,js,ix,jx,ie,je) * energy(je)
                          rsvect0(is,js,ix,ie) = rsvect0(is,js,ix,ie) &
                               + ctest(is,js,jx,ix,je,ie) * energy(je) &
                               * w_e(je)*w_xi(jx)
                       enddo
                    enddo
                 enddo
              enddo
           
              ! int v^2 C_test_ab(u_a^2 f0a,f0b) 
              
              do ix=1,n_xi
                 do ie=1,n_energy
                    rs(is,js) = rs(is,js) + w_e(ie)*w_xi(ix) &
                         * dens(is) * rsvec(is,js,ix,ie) * energy(ie) 
                 enddo
              enddo
           enddo
        enddo

        if (collision_kperp == 0) then
              do iv=1,nv  
                 is = is_v(iv)
                 ix = ix_v(iv)
                 ie = ie_v(iv)

                 do jv=1,nv
                    js = is_v(jv)
                    jx = ix_v(jv)
                    je = ie_v(jv)

                    if (abs(rs(is,js)) > epsilon(0.0)) then
                       cmat_loc(iv,jv) &
                            = cmat_loc(iv,jv) &
                            - temp(js)/temp(is) * dens(js) * dens_rot(it,js) &
                            * rsvec(is,js,ix,ie) &
                            / rs(is,js) * rsvect0(js,is,jx,je) 
                    endif
                 enddo
              enddo

        else
              rsvect0(:,:,:,:) = 0.0
              do is=1,n_species
                 do js=1,n_species
                    do ix=1,n_xi
                       do jx=1,n_xi
                          do ie=1,n_energy
                             do je=1,n_energy
                                rsvect0(is,js,ix,ie) = rsvect0(is,js,ix,ie) &
                                     + ctest(is,js,jx,ix,je,ie) * energy(je) &
                                     * w_e(je)*w_xi(jx) &
                                     * bessel(is,ix,ie,ic_loc,0) 
                             enddo
                          enddo
                       enddo
                    enddo
                 enddo
              enddo
              
              do iv=1,nv  
                 is = is_v(iv)
                 ix = ix_v(iv)
                 ie = ie_v(iv)
                 
                 do jv=1,nv
                    js = is_v(jv)
                    jx = ix_v(jv)
                    je = ie_v(jv)
                    
                    if (abs(rs(is,js)) > epsilon(0.0)) then
                       cmat_loc(iv,jv) &
                            = cmat_loc(iv,jv) &
                            - temp(js)/temp(is) * dens(js) * dens_rot(it,js) &
                            * rsvec(is,js,ix,ie) &
                            * bessel(is,ix,ie,ic_loc,0) / rs(is,js) &
                            * rsvect0(js,is,jx,je)
                    endif
                 enddo
              enddo
        endif
        
     endif
     
     end select

     ! Avoid singularity of n=0,p=0:
     if (px(ir) == 0 .and. my_toroidal == 0) then

        do iv=1,nv
           cmat_loc(iv,iv) =  1.0
           amat(iv,iv) = 1.0
        enddo

     else

        ! Already has field particle collisions
        do iv=1,nv
           do jv=1,nv
              rval = (0.5*delta_t) * cmat_loc(jv,iv)
              amat(jv,iv)     =  rval
              cmat_loc(jv,iv) = -rval
           enddo
        enddo

        do iv=1,nv

           is = is_v(iv)
           ix = ix_v(iv)
           ie = ie_v(iv)

           do jv=1,nv

              js = is_v(jv)
              jx = ix_v(jv)
              je = ie_v(jv)

              ! Collision component: Test particle
              if (is == js) then
                 do ks=1,n_species
                    rval = (0.5*delta_t) * ctest(is,ks,ix,jx,ie,je) &
                         * dens_rot(it,ks)
                    cmat_loc(iv,jv) = cmat_loc(iv,jv) - rval
                    amat(iv,jv) = amat(iv,jv) + rval
                 enddo
              endif

              ! Trapping 
              ! (not part of collision operator but contains xi-derivative)
              if (explicit_trap_flag == 0 .and. is == js .and. ie == je) then
                 rval = (0.5*delta_t) * (omega_trap(it,is) * vel(ie) &
                      + omega_rot_trap(it,is) / vel(ie)) &
                      * (1.0 - xi(ix)**2) * xi_deriv_mat(ix,jx)
                 cmat_loc(iv,jv) = cmat_loc(iv,jv) + rval
                 amat(iv,jv) = amat(iv,jv) - rval
              endif

              ! Rotation energy derivative
              ! (not part of collision operator but contains e-derivative)
              if (explicit_trap_flag == 0 .and. is == js .and. ix == jx) then
                 rval = (0.5*delta_t) * omega_rot_u(it,is) * xi(ix) &
                      * e_deriv1_rot_mat(ie,je)/sqrt(1.0*e_max)
                 cmat_loc(iv,jv) = cmat_loc(iv,jv) + rval
                 amat(iv,jv) = amat(iv,jv) - rval
              endif

              ! Finite-kperp test particle corrections 
              if (collision_model == 4 .and. collision_kperp == 1) then
                 if (is == js .and. jx == ix .and. je == ie) then
                    do ks=1,n_species
                       rval = (0.5*delta_t) &
                            * (-0.25*(k_perp(ic)*rho*vth(is)*mass(is) &
                            / (z(is)*bmag(it)))**2 * 2.0*energy(ie) &
                            * (klor_fac(is,ks)*nu_d(ie,is,ks) * (1+xi(ix)**2) &
                            + kdiff_fac(is,ks)*nu_par(ie,is,ks)* (1-xi(ix)**2)))
                       cmat_loc(iv,jv) = cmat_loc(iv,jv) - rval
                       amat(iv,jv) = amat(iv,jv) + rval
                    enddo
                 endif
              endif

              if (collision_field_model == 1) then

                 ! Poisson component l
                 if (my_toroidal == 0 .and. ae_flag == 1) then
                    ! Cannot include Poisson in collision matrix
                    ! for n=0 with ade because depends on theta
                    ! i.e. ne0 ~ phi - <phi>
                    !cmat_loc(iv,jv)    = cmat_loc(iv,jv) + 0.0
                    !amat(iv,jv)        = amat(iv,jv) + 0.0
                 else
                    rval =  z(is)/temp(is) * jvec_v(1,ic_loc,iv) &
                         / (k_perp(ic)**2 * lambda_debye**2 &
                         * dens_ele / temp_ele + sum_den_h(it)) &
                         * z(js)*dens(js)*dens_rot(it,js) &
                         * jvec_v(1,ic_loc,jv) * w_e(je) * w_xi(jx) 
                    cmat_loc(iv,jv) = cmat_loc(iv,jv) - rval
                    amat(iv,jv) = amat(iv,jv) - rval
                 endif

                 ! Ampere component
                 if (n_field > 1) then
                    rval =  z(is)/temp(is) * (jvec_v(2,ic_loc,iv) &
                         / (2.0*k_perp(ic)**2 * rho**2 / betae_unit & 
                         * dens_ele * temp_ele)) &
                         * z(js)*dens(js)*dens_rot(it,js) &
                         * jvec_v(2,ic_loc,jv) * w_e(je) * w_xi(jx)  
                    cmat_loc(iv,jv) = cmat_loc(iv,jv) + rval
                    amat(iv,jv) = amat(iv,jv) + rval
                 endif

                 ! Ampere Bpar component
                 if (n_field > 2) then
                    rval = jvec_v(3,ic_loc,iv) &
                         * (-0.5*betae_unit)/(dens_ele*temp_ele) &
                         * w_e(je)*w_xi(jx)*dens(js)*dens_rot(it,js)*temp(js) &
                         * jvec_v(3,ic_loc,jv)/(temp(is)/z(is))/(temp(js)/z(js))
                    cmat_loc(iv,jv) = cmat_loc(iv,jv) - rval
                    amat(iv,jv) = amat(iv,jv) - rval
                 endif

              endif

           enddo
        enddo

        ! constant part
        do iv=1,nv
           cmat_loc(iv,iv) = cmat_loc(iv,iv) + 1.0
           amat(iv,iv) = amat(iv,iv) + 1.0
        enddo

     endif

     ! H_bar = (1 - dt/2 C - Poisson)^(-1) * (1 + dt/2 C + Poisson) H
     ! Lapack factorization and inverse of LHS
     call DGESV(nv,nv,cmat_loc(:,:),size(cmat_loc,1), &
          i_piv,amat,size(amat,1),info)


     ! result in amat, transfer to the right cmat matrix
     if (collision_precision_mode /= 0) then
        do jv=1,nv
           cmat_stripes(:,jv,ic_loc) = 0.0
           do iv=1,nv
              dv = iv-jv
              if (abs(dv) .GT. collision_full_stripes) then
                 ! far from diagonal, keep low precision only
                 cmat_fp32(iv,jv,ic_loc) = amat(iv,jv)
              else
                 ! close to the diagonal, keep full precision
                 cmat_stripes(dv,jv,ic_loc) = amat(iv,jv)
                 ! set main matrix to 0, for ease of compute later
                 cmat_fp32(iv,jv,ic_loc) = 0.0
              endif
           enddo
        enddo
     else
        ! keep all cmat in full precision
        cmat(:,:,ic_loc) = amat(:,:)
     endif

  enddo
  deallocate(cmat_loc)
  deallocate(amat)

  if (collision_model == 4 .and. collision_kperp == 1 .and. &
       (collision_mom_restore == 1 .or. collision_ene_restore == 1)) then
     deallocate(bessel)
  end if


  if (collision_precision_mode /= 0) then
!$acc enter data copyin(cmat_stripes,cmat_fp32) if (gpu_bigmem_flag == 1)
  else
!$acc enter data copyin(cmat) if (gpu_bigmem_flag == 1)
  endif

  deallocate(i_piv)
  deallocate(nu_d)
  deallocate(nu_par)
  deallocate(rs)
  deallocate(rsvec)
  deallocate(rsvect0)
  deallocate(rsvect1)
  deallocate(ctest)
  deallocate(klor_fac)
  deallocate(kdiff_fac)

end subroutine cgyro_init_collision
