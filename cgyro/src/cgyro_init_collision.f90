subroutine cgyro_init_collision

  use timer_lib

  use cgyro_globals

  implicit none

  real, dimension(:,:,:), allocatable :: nu_d, nu_s, nu_par, nu_par_deriv
  real, dimension(:,:), allocatable :: rs
  real, dimension(:,:,:,:), allocatable :: rsvec, rsvect 

  real :: arg
  real :: xa, xb, tauinv_ab
  integer :: jv
  integer :: is,ir,it,ix,ie,js,je,jx,ks
  ! parameters for matrix solve
  real, dimension(:,:), allocatable :: amat
  real, dimension(:,:,:,:,:,:), allocatable :: ctest
  real, dimension(:,:,:,:,:), allocatable :: bessel

  if (collision_model == 5) then
     call cgyro_init_collision_simple
     return
  endif

  allocate(nu_d(n_energy,n_species,n_species))
  allocate(nu_s(n_energy,n_species,n_species))
  allocate(nu_par(n_energy,n_species,n_species))
  allocate(nu_par_deriv(n_energy,n_species,n_species))
  nu_d(:,:,:) = 0.0
  nu_s(:,:,:) = 0.0
  nu_par(:,:,:) = 0.0
  nu_par_deriv(:,:,:) = 0.0

  do ie=1,n_energy
     do is=1,n_species
        do js=1,n_species

           xa = vel(ie)
           xb = xa * vth(is) / vth(js)
           tauinv_ab = nu(is) * (1.0*z(js))**2 / (1.0*z(is))**2 &
                * dens(js)/dens(is)

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
                    nu_d(ie,is,js) = tauinv_ab * (1.0/xa**3)
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
              nu_s(ie,is,js) = nu_d(ie,is,js)

           case (3)

              ! Reduced Hirshman-Sigmar model
              ! (Fix for underflow)
              nu_d(ie,is,js) = tauinv_ab * (1.0/xa**3) &
                   * (exp(-xb*xb)/(xb*sqrt(pi)) &
                   + (1.0-1.0/(2.0*xb*xb)) * erf(xb))
              nu_s(ie,is,js) = tauinv_ab * (1.0/xa) &
                   * (-exp(-xb*xb)/(xb*sqrt(pi)) &
                   + (1.0/(2.0*xb*xb)) * erf(xb)) &
                   * (2.0*temp(is)/temp(js))*(1.0+mass(js)/mass(is))

           case(4)

              ! Ad hoc op
              ! (Fix for underflow)
              nu_d(ie,is,js) = tauinv_ab * (1.0/xa**3) &
                   * (exp(-xb*xb)/(xb*sqrt(pi)) &
                   + (1.0-1.0/(2.0*xb*xb)) * erf(xb))
              ! No i-e Lorentz
              if (is /= is_ele .and. js == is_ele) then
                 nu_d(ie,is,js) = 0.0
              endif

              ! Only ii, ee Diffusion
              if(is == js) then
                 nu_par(ie,is,js) = tauinv_ab * (2.0/xa**3) &
                      * (-exp(-xb*xb)/(xb*sqrt(pi)) &
                      + (1.0/(2.0*xb*xb)) * erf(xb))
                 nu_par_deriv(ie,is,js) = tauinv_ab * (2.0/xa**3) &
                      * (-3/xa * (-exp(-xb*xb)/(xb*sqrt(pi)) &
                      + (1.0/(2.0*xb*xb)) * erf(xb)) + vth(is)/vth(js) &
                      * (2.0*exp(-xb*xb)/(xb**2*sqrt(pi)) &
                      + 2.0*exp(-xb*xb)/sqrt(pi) - erf(xb)/xb**3))
              endif

           case(6)
              ! Const nu(energy)
              nu_d(ie,is,js) = tauinv_ab
              
           end select

        enddo
     enddo
  enddo

  ! Printout used for CGYRO paper.
  !if (i_proc == 0) then
  !   do ie=1,n_energy
  !      print '(i1,3(" &{\tt ",f9.5,"}")," \\")',&
  !           ie,vel(ie),w_e(ie),nu_d(ie,n_species,n_species)
  !   enddo
  !   print *,sum(w_e)
  !endif

  if (collision_model == 4 .and. collision_kperp == 1 .and. &
       (collision_mom_restore == 1 .or. collision_ene_restore == 1)) then
     allocate(bessel(n_species,n_xi,n_energy,nc_loc,0:1))
     ic_loc = 0
     do ic=nc1,nc2
        ic_loc = ic_loc+1
        it = it_c(ic)
        ir = ir_c(ic)
        do is=1,n_species   
           do ix=1,n_xi
              do ie=1,n_energy
                 arg = k_perp(ic)*rho*vth(is)*mass(is)&
                      /(z(is)*bmag(it)) *sqrt(2.0*energy(ie)) &
                      *sqrt(1.0-xi(ix)**2)
                 bessel(is,ix,ie,ic_loc,0) = bessel_j0(abs(arg))
                 bessel(is,ix,ie,ic_loc,1) = bessel_j1(abs(arg))
              enddo
           enddo
        enddo
     enddo
     ! asymptotic limits for electrons
     do is=1,n_species
        if (is == is_ele) then
           bessel(is,:,:,:,0) = 1.0
           bessel(is,:,:,:,1) = 0.0
        endif
     enddo
  endif

  allocate(ctest(n_species,n_species,n_xi,n_xi,n_energy,n_energy))
  allocate(rs(n_species,n_species))
  allocate(rsvec(n_species,n_species,n_xi,n_energy))
  allocate(rsvect(n_species,n_species,n_xi,n_energy))
  allocate(amat(nv,nv))

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

  ! U factor for HS0
  if (collision_model == 3 .and. collision_mom_restore == 1) then
     do is=1,n_species
        do ix=1,n_xi
           do ie=1,n_energy
              do js=1,n_species
                 do jx=1,n_xi
                    je = ie
                    ctest(is,js,ix,jx,ie,je) &
                         = ctest(is,js,ix,jx,ie,je) &
                         + (nu_d(ie,is,js)-nu_s(ie,is,js)) &
                         * 3.0 * xi(ix) * xi(jx) * w_xi(jx)
                 enddo
              enddo
           enddo
        enddo
     enddo
  endif

  ! Diffusion
  if (collision_model == 4 .and. collision_ene_diffusion == 1) then
     do is=1,n_species 
        do ix=1,n_xi
           do ie=1,n_energy
              do js=1,n_species
                 do je=1,n_energy
                    jx = ix                       
                    if (ie == je) then
                       ctest(is,js,ix,jx,ie,je) &
                            = ctest(is,js,ix,jx,ie,je) &
                            + (1.0-temp(is)/temp(js)) &
                            * (-nu_par_deriv(ie,is,js) * energy(ie)**1.5 & 
                            + nu_par(ie,is,js) & 
                            * (2.0*energy(ie)**2 - 5.0*energy(ie)))
                    endif
                    ctest(is,js,ix,jx,ie,je) &
                         = ctest(is,js,ix,jx,ie,je) &
                         + nu_par(ie,is,js) * 0.5 *energy(ie) &
                         * e_deriv2_mat(ie,je) / e_max
                    ctest(is,js,ix,jx,ie,je) &
                         = ctest(is,js,ix,jx,ie,je) &
                         + e_deriv1_mat(ie,je)/sqrt(1.0*e_max) &
                         * (nu_par_deriv(ie,is,js) * 0.5*energy(ie) &
                         + nu_par(ie,is,js) * (2.0*vel(ie) &
                         + (temp(is)/temp(js)-2.0) * energy(ie)**1.5))
                 enddo
              enddo
           enddo
        enddo
     enddo
  endif

  ! Collision field particle component
  amat(:,:)   = 0.0
  cmat(:,:,:) = 0.0

  select case (collision_model)

  case(1,2,3,6)
     if (collision_mom_restore == 1) then
        do is=1,n_species
           do js=1,n_species
              rs(is,js) = 0.0
              do ie=1,n_energy
                 rs(is,js) = rs(is,js) + w_e(ie)*nu_s(ie,is,js)*energy(ie)
              enddo
           enddo
        enddo

        ic_loc = 0
        do ic=nc1,nc2
           ic_loc = ic_loc+1

           do iv=1,nv  
              is = is_v(iv)
              ix = ix_v(iv)
              ie = ie_v(iv)

              do jv=1,nv
                 js = is_v(jv)
                 jx = ix_v(jv)
                 je = ie_v(jv)

                 if (abs(rs(is,js)) > epsilon(0.0)) then
                    cmat(iv,jv,ic_loc) = &
                         cmat(iv,jv,ic_loc) &
                         + 3.0 * (mass(js)/mass(is)) * (dens(js)/dens(is)) &
                         * (vth(js)/vth(is)) * nu_s(ie,is,js) &
                         * vel(ie) * xi(ix) &
                         * nu_s(je,js,is) * sqrt(energy(je)) &
                         * xi(jx) * w_e(je) * w_xi(jx) / rs(is,js)
                 endif
              enddo
           enddo
        enddo

     endif

  case(4)

     ! Momentum Restoring

     if (collision_mom_restore == 1) then

        ! C_test_ab(v_par f0a,f0b) / vth_a
        rsvec = 0.0
        do is=1,n_species
           do js=1,n_species
              do ix=1,n_xi
                 do jx=1,n_xi
                    do ie=1,n_energy
                       do je=1,n_energy
                          rsvec(is,js,ix,ie) = rsvec(is,js,ix,ie) &
                               + ctest(is,js,ix,jx,ie,je) &
                               * sqrt(2.0*energy(je)) * xi(jx) 
                       enddo
                    enddo
                 enddo
              enddo

              ! int v_par C_test_ab(v_par f0a,f0b) / (n_0a vth_a^2)
              rs(is,js) = 0.0
              do ix=1,n_xi
                 do ie=1,n_energy
                    rs(is,js) = rs(is,js) + w_e(ie)*w_xi(ix) &
                         * rsvec(is,js,ix,ie) * sqrt(2.0*energy(ie)) * xi(ix) 
                 enddo
              enddo
           enddo
        enddo

        if (collision_kperp == 0) then
           ic_loc = 0
           do ic=nc1,nc2
              ic_loc = ic_loc+1

              do iv=1,nv  
                 is = is_v(iv)
                 ix = ix_v(iv)
                 ie = ie_v(iv)

                 do jv=1,nv
                    js = is_v(jv)
                    jx = ix_v(jv)
                    je = ie_v(jv)

                    if (abs(rs(is,js))>epsilon(0.0)) then
                       cmat(iv,jv,ic_loc) &
                            = cmat(iv,jv,ic_loc) &
                            - mass(js)/mass(is) * (dens(js)/dens(is)) &
                            * (vth(js)/vth(is)) * rsvec(is,js,ix,ie) &
                            / rs(is,js) * rsvec(js,is,jx,je) &
                            * w_e(je)*w_xi(jx)
                    endif
                 enddo
              enddo
           enddo

        else
           ic_loc = 0
           do ic=nc1,nc2
              ic_loc = ic_loc+1
              it = it_c(ic)
              ir = ir_c(ic)

              do iv=1,nv  
                 is = is_v(iv)
                 ix = ix_v(iv)
                 ie = ie_v(iv)

                 do jv=1,nv
                    js = is_v(jv)
                    jx = ix_v(jv)
                    je = ie_v(jv)

                    ! EAB: 07/08/16 fixed bug, had (vtb/vta)**2
                    if (abs(rs(is,js)) > epsilon(0.)) then 
                       cmat(iv,jv,ic_loc) &
                            = cmat(iv,jv,ic_loc) &
                            - mass(js)/mass(is) * (dens(js)/dens(is)) &
                            * (vth(js)/vth(is)) * rsvec(is,js,ix,ie) &
                            * bessel(is,ix,ie,ic_loc,0) / rs(is,js) &
                            * rsvec(js,is,jx,je) * bessel(js,jx,je,ic_loc,0) &
                            * w_xi(jx)*w_e(je)
                       cmat(iv,jv,ic_loc) &
                            = cmat(iv,jv,ic_loc) &
                            - mass(js)/mass(is) * (dens(js)/dens(is)) &
                            * (vth(js)/vth(is)) * rsvec(is,js,ix,ie) &
                            * bessel(is,ix,ie,ic_loc,1) &
                            * sqrt(1.0-xi(ix)**2)/xi(ix) / rs(is,js) &
                            * rsvec(js,is,jx,je) * bessel(js,jx,je,ic_loc,1) &
                            * sqrt(1.0-xi(jx)**2)/xi(jx) * w_xi(jx)*w_e(je)
                    endif
                 enddo
              enddo
           enddo

        endif

     endif

     ! Energy Restoring

     if (collision_ene_restore == 1) then

        ! Only ii, ee Diffusion
              
        ! C_test_ab(u_a^2 f0a,f0b) and w_e u_a^2 C_test_ab(H_a)
        rsvec  = 0.0
        rsvect = 0.0
        rs(:,:) = 0.0

        do is=1,n_species
           js=is
           do ix=1,n_xi
              do jx=1,n_xi
                 do ie=1,n_energy
                    do je=1,n_energy
                       rsvec(is,js,ix,ie) = rsvec(is,js,ix,ie) &
                            + ctest(is,js,ix,jx,ie,je) * energy(je)
                       rsvect(is,js,ix,ie) = rsvect(is,js,ix,ie) &
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

        if (collision_kperp == 0) then
           ic_loc = 0
           do ic=nc1,nc2
              ic_loc = ic_loc+1

              do iv=1,nv  
                 is = is_v(iv)
                 ix = ix_v(iv)
                 ie = ie_v(iv)

                 do jv=1,nv
                    js = is_v(jv)
                    jx = ix_v(jv)
                    je = ie_v(jv)

                    if (abs(rs(is,js)) > epsilon(0.0)) then
                       !cmat(iv,jv,ic_loc) &
                       !     = cmat(iv,jv,ic_loc) &
                       !     - temp(js)/temp(is) * dens(js) &
                       !     * rsvec(is,js,ix,ie) &
                       !     / rs(is,js) * rsvec(js,is,jx,je) &
                       !     * w_xi(jx)*w_e(je)
                       cmat(iv,jv,ic_loc) &
                            = cmat(iv,jv,ic_loc) &
                            - temp(js)/temp(is) * dens(js) &
                            * rsvec(is,js,ix,ie) &
                            / rs(is,js) * rsvect(js,is,jx,je) 
                    endif
                 enddo
              enddo
           enddo

        else
           ic_loc = 0
           do ic=nc1,nc2
              ic_loc = ic_loc+1
              it = it_c(ic)
              ir = ir_c(ic)

              rsvect(:,:,:,:) = 0.0
              do is=1,n_species
                 js=is
                 do ix=1,n_xi
                    do jx=1,n_xi
                       do ie=1,n_energy
                          do je=1,n_energy
                             rsvect(is,js,ix,ie) = rsvect(is,js,ix,ie) &
                                  + ctest(is,js,jx,ix,je,ie) * energy(je) &
                                  * w_e(je)*w_xi(jx) &
                                  * bessel(is,ix,ie,ic_loc,0) 
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
                       !cmat(iv,jv,ic_loc) &
                       !     = cmat(iv,jv,ic_loc) &
                       !     - temp(js)/temp(is) * dens(js) &
                       !     * rsvec(is,js,ix,ie) &
                       !     * bessel(is,ix,ie,ic_loc,0) / rs(is,js) &
                       !     * rsvec(js,is,jx,je) * bessel(js,jx,je,ic_loc,0) &
                       !     * w_xi(jx)*w_e(je)
                       cmat(iv,jv,ic_loc) &
                            = cmat(iv,jv,ic_loc) &
                            - temp(js)/temp(is) * dens(js) &
                            * rsvec(is,js,ix,ie) &
                            * bessel(is,ix,ie,ic_loc,0) / rs(is,js) &
                            * rsvect(js,is,jx,je)
                    endif
                 enddo
              enddo
           enddo
        endif

     endif

  end select

  if (collision_model == 4 .and. collision_kperp == 1 .and. &
       (collision_mom_restore == 1 .or. collision_ene_restore == 1)) then
     deallocate(bessel)
  end if

  ! matrix solve parameters
  allocate(i_piv(nv))

  ! Construct the collision matrix

!$omp  parallel do  default(none) &
!$omp& shared(nc1,nc2,nv,n,delta_t,n_species,rho,is_ele,n_field) &
!$omp& shared(collision_model,collision_kperp,collision_field_model) &
!$omp& shared(ae_flag,lambda_debye,dens_ele,temp_ele) &
!$omp& shared(betae_unit,sum_den_h) &
!$omp& shared(it_c,ir_c,px,is_v,ix_v,ie_v,ctest,xi_deriv_mat) &
!$omp& shared(temp,jvec_v,omega_trap,dens,energy,vel) &
!$omp& shared(k_perp,vth,mass,z,bmag,nu_d,xi,nu_par,w_e,w_xi) &
!$omp& private(ic,ic_loc,it,ir,info) &
!$omp& private(iv,is,ix,ie,jv,js,jx,je,ks) &
!$omp& private(amat,i_piv) &
!$omp& shared(cmat)
  do ic=nc1,nc2
   
     ic_loc = ic-nc1+1

     it = it_c(ic)
     ir = ir_c(ic)

     ! Initialize work array
     amat(:,:) = 0.0

     ! Avoid singularity of n=0,p=0:
     if (px(ir) == 0 .and. n == 0) then

        do iv=1,nv
           cmat(iv,iv,ic_loc) =  1.0
           amat(iv,iv) = 1.0
        enddo

     else

        ! Already has field particle collisions
        amat(:,:)        = (0.5*delta_t)  * cmat(:,:,ic_loc)
        cmat(:,:,ic_loc) = -(0.5*delta_t) * cmat(:,:,ic_loc)

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
                    cmat(iv,jv,ic_loc) = cmat(iv,jv,ic_loc) &
                         - (0.5*delta_t) * ctest(is,ks,ix,jx,ie,je)
                    amat(iv,jv) = amat(iv,jv) &
                         + (0.5*delta_t) * ctest(is,ks,ix,jx,ie,je)
                 enddo
              endif

              ! Trapping (not part of collision operator but contains xi-derivative)
              if (collision_model /= 6) then
                 if (is == js .and. ie == je) then
                    cmat(iv,jv,ic_loc) = cmat(iv,jv,ic_loc) &
                         + (0.5*delta_t) * omega_trap(it,is) &
                         * vel(ie) * (1.0 - xi(ix)**2) &
                         * xi_deriv_mat(ix,jx) 
                    amat(iv,jv) = amat(iv,jv) &
                         - (0.5*delta_t) * omega_trap(it,is) &
                         * vel(ie) * (1.0 - xi(ix)**2) &
                         * xi_deriv_mat(ix,jx) 
                 endif
              endif

              ! Finite-kperp test particle corrections 
              if (collision_model == 4 .and. collision_kperp == 1) then
                 if (is == js .and. jx == ix .and. je == ie) then
                    do ks=1,n_species
                       cmat(iv,jv,ic_loc) = cmat(iv,jv,ic_loc) &
                            - (0.5*delta_t) &
                            * (-0.25*(k_perp(ic)*rho*vth(is)*mass(is) &
                            / (z(is)*bmag(it)))**2 * 2.0*energy(ie) &
                            * (nu_d(ie,is,ks) * (1+xi(ix)**2) &
                            + nu_par(ie,is,ks) * (1-xi(ix)**2)) )
                       amat(iv,jv) = amat(iv,jv) &
                            + (0.5*delta_t) &
                            * (-0.25*(k_perp(ic)*rho*vth(is)*mass(is) &
                            / (z(is)*bmag(it)))**2 * 2.0*energy(ie) &
                            * (nu_d(ie,is,ks) * (1+xi(ix)**2) &
                            + nu_par(ie,is,ks) * (1-xi(ix)**2)) )
                    enddo
                 endif
              endif

              if (collision_field_model == 1) then

                 ! Poisson component l
                 if (n == 0 .and. ae_flag == 1) then
                    ! Cannot include Poisson in collision matrix
                    ! for n=0 with ade because depends on theta
                    ! i.e. ne0 ~ phi - <phi>
                    cmat(iv,jv,ic_loc) = cmat(iv,jv,ic_loc) + 0.0
                    amat(iv,jv)        = amat(iv,jv) + 0.0
                 else
                    cmat(iv,jv,ic_loc) = cmat(iv,jv,ic_loc) &
                         - z(is)/temp(is) * jvec_v(1,ic_loc,iv) &
                         / (k_perp(ic)**2 * lambda_debye**2 &
                         * dens_ele / temp_ele + sum_den_h) &
                         * z(js)*dens(js) &
                         * jvec_v(1,ic_loc,jv) * w_e(je) * w_xi(jx) 
                    amat(iv,jv) = amat(iv,jv) &
                         - z(is)/temp(is) * jvec_v(1,ic_loc,iv) &
                         / (k_perp(ic)**2 * lambda_debye**2 &
                         * dens_ele / temp_ele + sum_den_h) &
                         * z(js)*dens(js) &
                         * jvec_v(1,ic_loc,jv) * w_e(je) * w_xi(jx) 
                 endif

                 ! Ampere component
                 if (n_field > 1) then
                    cmat(iv,jv,ic_loc) = cmat(iv,jv,ic_loc) &
                         + z(is)/temp(is) * (jvec_v(2,ic_loc,iv) &
                         / (2.0*k_perp(ic)**2 * rho**2 / betae_unit & 
                         * dens_ele * temp_ele)) &
                         * z(js)*dens(js) &
                         * jvec_v(2,ic_loc,jv) * w_e(je) * w_xi(jx)  
                    amat(iv,jv) = amat(iv,jv) &
                         + z(is)/temp(is) * (jvec_v(2,ic_loc,iv) &
                         / (2.0*k_perp(ic)**2 * rho**2 / betae_unit & 
                         * dens_ele * temp_ele)) &
                         * z(js)*dens(js) &
                         * jvec_v(2,ic_loc,jv) * w_e(je) * w_xi(jx) 
                 endif

                 ! Ampere Bpar component
                 if (n_field > 2) then
                    cmat(iv,jv,ic_loc) = cmat(iv,jv,ic_loc) &
                         - jvec_v(3,ic_loc,iv) &
                         * (-0.5*betae_unit)/(dens_ele*temp_ele) &
                         * w_e(je)*w_xi(jx)*dens(js)*temp(js) &
                         * jvec_v(3,ic_loc,jv)/(temp(is)/z(is))/(temp(js)/z(js))
                    amat(iv,jv) = amat(iv,jv) &
                         - jvec_v(3,ic_loc,iv) &
                         * (-0.5*betae_unit)/(dens_ele*temp_ele) &
                         * w_e(je)*w_xi(jx)*dens(js)*temp(js) &
                         * jvec_v(3,ic_loc,jv)/(temp(is)/z(is))/(temp(js)/z(js))
                 endif

              endif

           enddo
        enddo

        ! constant part
        do iv=1,nv
           cmat(iv,iv,ic_loc) = cmat(iv,iv,ic_loc) + 1.0
           amat(iv,iv) = amat(iv,iv) + 1.0
        enddo

     endif

     ! H_bar = (1 - dt/2 C - Poisson)^(-1) * (1 + dt/2 C + Poisson) H
     ! Lapack factorization and inverse of LHS
     call DGESV(nv,nv,cmat(:,:,ic_loc),size(cmat,1), &
          i_piv,amat,size(amat,1),info)
     cmat(:,:,ic_loc) = amat(:,:)

  enddo
!$acc  enter data copyin(cmat)

  deallocate(amat)
  deallocate(i_piv)
  deallocate(nu_d)
  deallocate(nu_s)
  deallocate(nu_par)
  deallocate(nu_par_deriv)
  deallocate(rs)
  deallocate(rsvec)
  deallocate(rsvect)
  deallocate(ctest)

end subroutine cgyro_init_collision
