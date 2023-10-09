subroutine cgyro_init_collision

  use timer_lib

  use cgyro_globals
  use mpi
  use cgyro_io

  implicit none

  real, dimension(:,:,:), allocatable :: nu_d, nu_par
  real, dimension(:,:), allocatable :: rs
  real, dimension(:,:,:,:), allocatable :: rsvec, rsvect0
  real, dimension(:,:), allocatable :: klor_fac, kdiff_fac
  real :: rsvtmp, rsvtmp0
  real :: my_dens2_rot,my_bj0,my_bj1

  character(len=160) :: msg
  real :: arg,xi_s1s,xi_prop
  real :: xa, xb, tauinv_ab
  real :: rval
  integer :: jv
  integer :: is,ir,it,ix,ie,js,je,jx,ks
  integer :: itor
  ! parameters for matrix solve
  real, dimension(:,:), allocatable :: cmat_base1,cmat_base2
  real, dimension(:,:), allocatable :: amat,cmat_loc
  real, dimension(:,:,:,:,:,:), allocatable :: ctest
  real, dimension(:,:,:,:), allocatable :: bessel
  ! diagnostics
  real :: my_cmat_fp32
  real :: amat_sum, cmat_sum, cmat_diff, cmat_rel_diff
  real :: cmat32_sum, cmat32_diff
  real :: cmat_diff_global_loc, cmat32_diff_global_loc
  ! use real as 32-bit int may overflow
  real, dimension(8:19) :: cmap_fp32_error_abs_cnt_loc
  real, dimension(8:18) :: cmap_fp32_error_rel_cnt_loc
  real, dimension(8:19) :: cmap_fp32_error_abs_cnt
  real, dimension(8:18) :: cmap_fp32_error_rel_cnt
  real, dimension(2) :: cmap_fp32_error_sum_loc
  real, dimension(2) :: cmap_fp32_error_sum

  if (collision_model == 5) then
     call cgyro_init_collision_simple
     return
  endif

  cmap_fp32_error_abs_cnt_loc(:) = 0
  cmap_fp32_error_rel_cnt_loc(:) = 0

  allocate(nu_d(n_energy,n_species,n_species))
  allocate(nu_par(n_energy,n_species,n_species))
  allocate(klor_fac(n_species,n_species))
  allocate(kdiff_fac(n_species,n_species))
  nu_d(:,:,:) = 0.0
  nu_par(:,:,:) = 0.0
  klor_fac(:,:) = 0.0
  kdiff_fac(:,:) = 0.0
  
  do js=1,n_species
     do is=1,n_species
        do ie=1,n_energy

           xa = vel(ie)
           xb = xa*vth(is)/vth(js)
           tauinv_ab = nu(is)*z(js)**2/z(is)**2*dens(js)/dens(is)
              
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
     allocate(bessel(0:1,nv,nc_loc,nt1:nt2))
!$omp parallel do collapse(2) private(ic_loc,it,ie,ix,is,iv,arg,xi_s1s,xi_prop)
     do itor=nt1,nt2
      do ic=nc1,nc2
        ic_loc = ic-nc1+1
        it = it_c(ic)
        do iv=1,nv
              is = is_v(iv)
              ix = ix_v(iv)
              ie = ie_v(iv)
              xi_s1s = sqrt(1.0-xi(ix)**2)
              xi_prop = xi_s1s / xi(ix)
                 arg = k_perp(ic,itor)*rho*vth(is)*mass(is)&
                      /(z(is)*bmag(it)) * vel2(ie) * xi_s1s
                 bessel(0,iv,ic_loc,itor) = bessel_j0(arg)
                 ! always used with the correction, so do it here
                 bessel(1,iv,ic_loc,itor) = bessel_j1(arg) * xi_prop
        enddo
      enddo
     enddo
  endif

  allocate(ctest(n_species,n_species,n_xi,n_xi,n_energy,n_energy))

  ! Collision test particle component
  ! ctest

  ! Lorentz
!$omp parallel do collapse(2) private(is,js,ix,jx,ie,je)
  do je=1,n_energy
     do ie=1,n_energy
        if (je == ie ) then
           do jx=1,n_xi
              do ix=1,n_xi
                 do js=1,n_species
                    do is=1,n_species
                       ctest(is,js,ix,jx,ie,je) = &
                             + xi_lor_mat(ix,jx) *0.5*nu_d(ie,is,js)
                    enddo
                 enddo
              enddo
           enddo
        else
           ctest(:,:,:,:,ie,je) = 0
        endif
     enddo
  enddo

  ! Diffusion
  if (collision_model == 4 .and. collision_ene_diffusion == 1) then
!$omp parallel do collapse(2) private(is,js,ix,jx,ie,je)
     do je=1,n_energy
        do ie=1,n_energy
           do ix=1,n_xi
              jx = ix
              do js=1,n_species
                 do is=1,n_species
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

  allocate(cmat_base1(nv,nv))
  allocate(cmat_base2(nv,nv))

  allocate(rs(n_species,n_species))
  allocate(rsvec(n_species,n_species,n_xi,n_energy))
  allocate(rsvect0(n_species,n_species,n_xi,n_energy))

  select case (collision_model)

  case(2)
     if (collision_mom_restore == 1) then
        do js=1,n_species
           do is=1,n_species
              rs(is,js) = 0.0
              do ie=1,n_energy
                 rs(is,js) = rs(is,js) + w_e(ie)*nu_d(ie,is,js)*energy(ie)
              enddo
           enddo
        enddo

!$omp parallel do collapse(2) private(is,js,ix,jx,ie,je)
        do jv=1,nv
           do iv=1,nv
              js = is_v(jv)
              jx = ix_v(jv)
              je = ie_v(jv)

              is = is_v(iv)
              ix = ix_v(iv)
              ie = ie_v(iv)

              if (abs(rs(is,js)) > epsilon(0.0)) then
                 cmat_base1(iv,jv) = &
                         + 3.0 * (mass(js)/mass(is)) &
                         * (vth(js)/vth(is)) * nu_d(ie,is,js) &
                         * vel(ie) * xi(ix) &
                         * nu_d(je,js,is) * vel(je) &
                         * xi(jx) * w_exi(je,jx) / rs(is,js) &
                         / dens(is)
             else
                cmat_base1(iv,jv) = 0
             endif
           enddo
       enddo

     endif

  case(4)

     ! Momentum Restoring

     if (collision_mom_restore == 1) then

        ! C_test_ab(v_par f0a,f0b) and w_e v_par C_test_ab(H_a)
        rs(:,:) = 0.0
        
        do ie=1,n_energy
           do ix=1,n_xi
              do js=1,n_species
                 do is=1,n_species
                    rsvtmp = 0
                    rsvtmp0 = 0
                    do je=1,n_energy
                       do jx=1,n_xi
                          rsvtmp = rsvtmp &
                               + ctest(is,js,ix,jx,ie,je) &
                               * vel2(je) * xi(jx)
                          rsvtmp0 = rsvtmp0 &
                               + ctest(is,js,jx,ix,je,ie) &
                               * vel2(je) * xi(jx) &
                               * w_exi(je,jx)
                       enddo
                    enddo
                    rsvtmp = rsvtmp * vth(is)
                    rsvtmp0 = rsvtmp0 * vth(is)
                    rsvec(is,js,ix,ie) = rsvtmp
                    rsvect0(is,js,ix,ie) = rsvtmp0
                    ! int v_par C_test_ab(v_par f0a,f0b) / n_0a
                    rs(is,js) = rs(is,js) + w_exi(ie,ix) * dens(is) &
                         * rsvtmp * vel2(ie) * xi(ix) &
                         * vth(is)
                 enddo
              enddo
           enddo
        enddo

        if (collision_kperp == 0) then
!$omp parallel do collapse(2) private(is,js,ix,jx,ie,je)
           do jv=1,nv
              do iv=1,nv
                 js = is_v(jv)
                 jx = ix_v(jv)
                 je = ie_v(jv)

                 is = is_v(iv)
                 ix = ix_v(iv)
                 ie = ie_v(iv)

                 if (abs(rs(is,js))>epsilon(0.0)) then
                    cmat_base1(iv,jv) = &
                            - mass(js)/mass(is) &
                            * rsvec(is,js,ix,ie) / rs(is,js) &
                            * rsvect0(js,is,jx,je)
                 else
                    cmat_base1(iv,jv) = 0
                 endif
              enddo
           enddo

        else
!$omp parallel do collapse(2) private(is,js,ix,jx,ie,je)
              do jv=1,nv
                 do iv=1,nv
                    js = is_v(jv)
                    jx = ix_v(jv)
                    je = ie_v(jv)

                    is = is_v(iv)
                    ix = ix_v(iv)
                    ie = ie_v(iv)

                    if (abs(rs(is,js)) > epsilon(0.)) then 
                       cmat_base1(iv,jv) = &
                            - mass(js)/mass(is) &
                            * rsvec(is,js,ix,ie) / rs(is,js) &
                            * rsvect0(js,is,jx,je)
                    else
                       cmat_base1(iv,jv) = 0
                    endif
                 enddo
              enddo
        endif

     endif

     ! Energy Restoring

     if (collision_ene_restore == 1) then
              
        ! C_test_ab(u_a^2 f0a,f0b) and w_e u_a^2 C_test_ab(H_a)
        rs(:,:) = 0.0

        do ie=1,n_energy
           do ix=1,n_xi
              do js=1,n_species
                 do is=1,n_species
                    rsvtmp = 0
                    rsvtmp0 = 0
                    do je=1,n_energy
                       do jx=1,n_xi
                          rsvtmp = rsvtmp &
                               + ctest(is,js,ix,jx,ie,je) * energy(je)
                          rsvtmp0 = rsvtmp0 &
                               + ctest(is,js,jx,ix,je,ie) * energy(je) &
                               * w_exi(je,jx)
                       enddo
                    enddo
                    rsvec(is,js,ix,ie) = rsvtmp
                    rsvect0(is,js,ix,ie) = rsvtmp0
                    ! int v^2 C_test_ab(u_a^2 f0a,f0b) 
                    rs(is,js) = rs(is,js) + w_exi(ie,ix) &
                         * dens(is) * rsvtmp * energy(ie) 
                 enddo
              enddo
           enddo
        enddo

        if (collision_kperp == 0) then
!$omp parallel do collapse(2) private(is,js,ix,jx,ie,je)
              do jv=1,nv
                 do iv=1,nv
                    js = is_v(jv)
                    jx = ix_v(jv)
                    je = ie_v(jv)

                    is = is_v(iv)
                    ix = ix_v(iv)
                    ie = ie_v(iv)

                    if (abs(rs(is,js)) > epsilon(0.0)) then
                       cmat_base2(iv,jv) = &
                            - temp(js)/temp(is) &
                            * rsvec(is,js,ix,ie) / rs(is,js) &
                            * rsvect0(js,is,jx,je)
                    else
                       cmat_base2(iv,jv) = 0
                    endif
                 enddo
              enddo

        else
!$omp parallel do collapse(2) private(is,js,ix,jx,ie,je)
              do jv=1,nv
                 do iv=1,nv
                    js = is_v(jv)
                    jx = ix_v(jv)
                    je = ie_v(jv)

                    is = is_v(iv)
                    ix = ix_v(iv)
                    ie = ie_v(iv)
                 
                    if (abs(rs(is,js)) > epsilon(0.0)) then
                       cmat_base2(iv,jv) = &
                            - temp(js)/temp(is) &
                            * rsvec(is,js,ix,ie) / rs(is,js) &
                            * rsvect0(js,is,jx,je)
                    else
                       cmat_base2(iv,jv) = 0
                    endif
                 enddo
              enddo
        endif
        
     endif
     
  end select
  deallocate(rs)
  deallocate(rsvec)
  deallocate(rsvect0)

  ! matrix solve parameters
  allocate(i_piv(nv))

  allocate(amat(nv,nv))
  allocate(cmat_loc(nv,nv))

  ! Construct the collision matrix
  cmat_diff_global_loc = 0.0
  cmat32_diff_global_loc = 0.0

!$omp  parallel do collapse(2) default(none) &
!$omp& shared(nc1,nc2,nt1,nt2,nv,delta_t,n_species,rho,is_ele,n_field,n_energy,n_xi) &
!$omp& shared(collision_kperp,collision_field_model,explicit_trap_flag) &
!$omp& firstprivate(collision_model,collision_mom_restore,collision_ene_restore) &
!$omp& shared(ae_flag,lambda_debye,dens_ele,temp_ele,dens_rot,dens2_rot) &
!$omp& shared(cmat_base1,cmat_base2,bessel) &
!$omp& shared(betae_unit,sum_den_h) &
!$omp& shared(it_c,ir_c,px,is_v,ix_v,ie_v,ctest,xi_deriv_mat) &
!$omp& shared(temp,jvec_v,omega_trap,dens,energy,vel,vel2) &
!$omp& shared(omega_rot_trap,omega_rot_u,e_deriv1_mat,e_deriv1_rot_mat,e_max) &
!$omp& shared(xi_lor_mat) &
!$omp& shared(k_perp,vth,mass,z,bmag,nu_d,xi,nu_par,w_e,w_exi) &
!$omp& shared(klor_fac,kdiff_fac) &
!$omp& private(ic,ic_loc,it,ir,info,rval,my_dens2_rot,my_bj0,my_bj1) &
!$omp& private(iv,is,ix,ie,jv,js,jx,je,ks) &
!$omp& private(amat_sum,cmat_sum,cmat_diff,cmat_rel_diff,cmat32_sum,cmat32_diff) &
!$omp& private(amat,cmat_loc,my_cmat_fp32,i_piv) &
!$omp& firstprivate(collision_precision_mode,n_low_energy) &
!$omp& shared(cmat,cmat_fp32,cmat_stripes,cmat_e1) &
!$omp& reduction(+:cmat_diff_global_loc,cmat32_diff_global_loc) &
!$omp& reduction(+:cmap_fp32_error_abs_cnt_loc,cmap_fp32_error_rel_cnt_loc)
  do itor=nt1,nt2
   do ic=nc1,nc2
   
     ic_loc = ic-nc1+1

     it = it_c(ic)
     ir = ir_c(ic)

     ! Collision field particle component
     amat(:,:)   = 0.0

     select case (collision_model)

     case(2)
     if (collision_mom_restore == 1) then
        do jv=1,nv
           js = is_v(jv)
           my_dens2_rot = dens2_rot(it,js)

           do iv=1,nv

                 ! we know we are the first one to modify cmat_loc, so no need for +=
                 cmat_loc(iv,jv) = &
                         cmat_base1(iv,jv) &
                         * my_dens2_rot
           enddo
       enddo

     else
       cmat_loc(:,:) = 0.0
     endif

     case(4)

     ! Momentum Restoring

     if (collision_mom_restore == 1) then

        if (collision_kperp == 0) then
           do jv=1,nv
              js = is_v(jv)
              my_dens2_rot = dens2_rot(it,js)

              do iv=1,nv

                    ! we know we are the first one to modify cmat_loc, so no need for +=
                    cmat_loc(iv,jv) = &
                            cmat_base1(iv,jv) &
                            * my_dens2_rot
              enddo
           enddo

        else
              do jv=1,nv
                 js = is_v(jv)
                 my_dens2_rot = dens2_rot(it,js)
                 my_bj0 = bessel(0,jv,ic_loc,itor)
                 my_bj1 = bessel(1,jv,ic_loc,itor)

                 do iv=1,nv

                       ! we know we are the first one to modify cmat_loc, so no need for +=
                       cmat_loc(iv,jv) = &
                            cmat_base1(iv,jv) &
                            * (  bessel(0,iv,ic_loc,itor) * my_bj0 &
                               + bessel(1,iv,ic_loc,itor) * my_bj1 ) &
                            * my_dens2_rot
                 enddo
              enddo
        endif

     else
       cmat_loc(:,:) = 0.0
     endif

     ! Energy Restoring

     if (collision_ene_restore == 1) then

        if (collision_kperp == 0) then
              do jv=1,nv
                 js = is_v(jv)
                 my_dens2_rot = dens2_rot(it,js)

                 do iv=1,nv

                       ! likely not the first, use +=
                       cmat_loc(iv,jv) &
                            = cmat_loc(iv,jv) &
                            + cmat_base2(iv,jv) &
                            * my_dens2_rot
                 enddo
              enddo

        else
              do jv=1,nv
                 js = is_v(jv)
                 my_dens2_rot = dens2_rot(it,js)
                 my_bj0 = bessel(0,jv,ic_loc,itor)
                    
                 do iv=1,nv
                       ! likely not the first, use +=
                       cmat_loc(iv,jv) &
                            = cmat_loc(iv,jv) &
                            + cmat_base2(iv,jv) &
                            * bessel(0,iv,ic_loc,itor) * my_bj0 &
                            * my_dens2_rot
                 enddo
              enddo
        endif
        
     endif

     case default
        cmat_loc(:,:) = 0.0
     
     end select

     ! Avoid singularity of n=0,p=0:
     if (px(ir) == 0 .and. itor == 0) then

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

        do jv=1,nv

           js = is_v(jv)
           jx = ix_v(jv)
           je = ie_v(jv)

           do iv=1,nv

              is = is_v(iv)
              ix = ix_v(iv)
              ie = ie_v(iv)

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
                 rval = (0.5*delta_t) * (omega_trap(it,is,itor) * vel(ie) &
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
                            * (-0.25*(k_perp(ic,itor)*rho*vth(is)*mass(is) &
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
                 if (itor == 0 .and. ae_flag == 1) then
                    ! Cannot include Poisson in collision matrix
                    ! for n=0 with ade because depends on theta
                    ! i.e. ne0 ~ phi - <phi>
                    !cmat_loc(iv,jv)    = cmat_loc(iv,jv) + 0.0
                    !amat(iv,jv)        = amat(iv,jv) + 0.0
                 else
                    rval =  z(is)/temp(is) * jvec_v(1,ic_loc,itor,iv) &
                         / (k_perp(ic,itor)**2 * lambda_debye**2 &
                         * dens_ele / temp_ele + sum_den_h(it)) &
                         * z(js)*dens2_rot(it,js) &
                         * jvec_v(1,ic_loc,itor,jv) * w_exi(je,jx) 
                    cmat_loc(iv,jv) = cmat_loc(iv,jv) - rval
                    amat(iv,jv) = amat(iv,jv) - rval
                 endif

                 ! Ampere component
                 if (n_field > 1) then
                    rval =  z(is)/temp(is) * (jvec_v(2,ic_loc,itor,iv) &
                         / (2.0*k_perp(ic,itor)**2 * rho**2 / betae_unit & 
                         * dens_ele * temp_ele)) &
                         * z(js)*dens2_rot(it,js) &
                         * jvec_v(2,ic_loc,itor,jv) * w_exi(je,jx)  
                    cmat_loc(iv,jv) = cmat_loc(iv,jv) + rval
                    amat(iv,jv) = amat(iv,jv) + rval
                 endif

                 ! Ampere Bpar component
                 if (n_field > 2) then
                    rval = jvec_v(3,ic_loc,itor,iv) &
                         * (-0.5*betae_unit)/(dens_ele*temp_ele) &
                         * w_exi(je,jx)*dens2_rot(it,js)*temp(js) &
                         * jvec_v(3,ic_loc,itor,jv)/(temp(is)/z(is))/(temp(js)/z(js))
                    cmat_loc(iv,jv) = cmat_loc(iv,jv) - rval
                    amat(iv,jv) = amat(iv,jv) - rval
                 endif

              endif

              ! constant part
              if ( iv == jv ) then
                 cmat_loc(iv,iv) = cmat_loc(iv,iv) + 1.0
                 amat(iv,iv) = amat(iv,iv) + 1.0
              endif
           enddo
        enddo


     endif

     ! H_bar = (1 - dt/2 C - Poisson)^(-1) * (1 + dt/2 C + Poisson) H
     ! Lapack factorization and inverse of LHS
     call DGESV(nv,nv,cmat_loc(:,:),size(cmat_loc,1), &
          i_piv,amat,size(amat,1),info)


     ! result in amat, transfer to the right cmat matrix
     if (collision_precision_mode /= 0) then
        ! keep all cmat in fp32 precision
        cmat_fp32(:,:,ic_loc,itor) = amat(:,:)
        ! keep the remaining precision for select elements
        do jv=1,nv
           je = ie_v(jv)
           js = is_v(jv)
           jx = ix_v(jv)
           amat_sum = 0.0
           cmat_sum = 0.0
           cmat32_sum = 0.0
           do iv=1,nv
              ! using abs values, as I am gaugung precision errors, and this avoid symmetry cancellations
              amat_sum = amat_sum + abs(amat(iv,jv))
              cmat32_sum = cmat32_sum + abs(cmat_fp32(iv,jv,ic_loc,itor))
              ie = ie_v(iv)
              is = is_v(iv)
              ix = ix_v(iv)
              my_cmat_fp32 = cmat_fp32(iv,jv,ic_loc,itor)
              if (ie<=n_low_energy) then ! always keep all detail for lowest energy
                 cmat_e1(ix,is,ie,jv,ic_loc,itor) = amat(iv,jv) - my_cmat_fp32
                 ! my_cmat_fp32 not used for the original purpose anymore, reuse to represent the reduced precsion
                 my_cmat_fp32 = my_cmat_fp32 + cmat_e1(ix,is,ie,jv,ic_loc,itor) ! my_cmat_fp32 is actually fp64, so sum OK
              else ! only keep if energy and species the same
                 if ((je == ie) .AND. (js == is)) then
                    cmat_stripes(ix,is,ie,jx,ic_loc,itor) = amat(iv,jv) - my_cmat_fp32
                    ! my_cmat_fp32 not used for the original purpose anymore, reuse to represent the reduced precsion
                    my_cmat_fp32 = my_cmat_fp32 + cmat_stripes(ix,is,ie,jx,ic_loc,itor) ! my_cmat_fp32 is actually fp64, so sum OK
                 endif
              endif
              cmat_sum = cmat_sum + abs(my_cmat_fp32)
           enddo
           cmat_diff = abs(amat_sum-cmat_sum)
           cmat_diff_global_loc = cmat_diff_global_loc + cmat_diff
           cmat32_diff = abs(amat_sum-cmat32_sum)
           cmat32_diff_global_loc = cmat32_diff_global_loc + cmat32_diff
           cmat_rel_diff = cmat_diff/abs(amat_sum)
           ! reuse amt_sum as 10^-(ie), to reduce numbe rof variables in use 
           do ie=8,18
             amat_sum = 10.0 ** (-ie)
             if (cmat_diff>amat_sum) then
               cmap_fp32_error_abs_cnt_loc(ie) = cmap_fp32_error_abs_cnt_loc(ie) + 1
             endif
             if (cmat_rel_diff>amat_sum) then
               cmap_fp32_error_rel_cnt_loc(ie) = cmap_fp32_error_rel_cnt_loc(ie) + 1
             endif
           enddo
           ! use 19 as absolute counter for normalization
           cmap_fp32_error_abs_cnt_loc(19) = cmap_fp32_error_abs_cnt_loc(19) + 1
        enddo
     else
        ! keep all cmat in full precision
        cmat(:,:,ic_loc,itor) = amat(:,:)
     endif

   enddo
  enddo
  deallocate(cmat_loc)
  deallocate(amat)

  if (collision_model == 4 .and. collision_kperp == 1 .and. &
       (collision_mom_restore == 1 .or. collision_ene_restore == 1)) then
     deallocate(bessel)
  end if


  if (collision_precision_mode /= 0) then
#if defined(OMPGPU)
     ! no async for OMPGPU for now
!$omp target enter data map(to:cmat_fp32,cmat_stripes,cmat_e1) if (gpu_bigmem_flag == 1)
#elif defined(_OPENACC)
!$acc enter data copyin(cmat_fp32,cmat_stripes,cmat_e1) async if (gpu_bigmem_flag == 1)
#endif
     call MPI_ALLREDUCE(cmap_fp32_error_abs_cnt_loc,&
          cmap_fp32_error_abs_cnt,&
          12,&
          MPI_DOUBLE_PRECISION,&
          MPI_SUM,&
          CGYRO_COMM_WORLD,&
          i_err)

     call MPI_ALLREDUCE(cmap_fp32_error_rel_cnt_loc,&
          cmap_fp32_error_rel_cnt,&
          11,&
          MPI_DOUBLE_PRECISION,&
          MPI_SUM,&
          CGYRO_COMM_WORLD,&
          i_err)

     cmap_fp32_error_sum_loc(1) = cmat_diff_global_loc
     cmap_fp32_error_sum_loc(2) = cmat32_diff_global_loc
     call MPI_ALLREDUCE(cmap_fp32_error_sum_loc,&
          cmap_fp32_error_sum,&
          2,&
          MPI_DOUBLE_PRECISION,&
          MPI_SUM,&
          CGYRO_COMM_WORLD,&
          i_err)

     if (i_proc==0) then
        ! reuse amat_sum to reduce number of variables
        amat_sum = 0.01 * cmap_fp32_error_abs_cnt(19)
        write (msg, "(A,A,A,A,A,A,A,A,A)") &
                    "                           ", &
                    " >1.e-8", "    e-9", &
                    "   e-10", "   e-11", &
                    "   e-12", "   e-13", &
                    "   e-14", "   e-15"
        call cgyro_info(msg)
        write (msg, "(A,4(F6.2,A),4(F6.1,A))") &
                    "Abs cmat_fp32 error rates: ", cmap_fp32_error_abs_cnt(8)/amat_sum, &
                    "%",  cmap_fp32_error_abs_cnt(9)/amat_sum, &
                    "%", cmap_fp32_error_abs_cnt(10)/amat_sum, &
                    "%", cmap_fp32_error_abs_cnt(11)/amat_sum, &
                    "%", cmap_fp32_error_abs_cnt(12)/amat_sum, &
                    "%", cmap_fp32_error_abs_cnt(13)/amat_sum, &
                    "%", cmap_fp32_error_abs_cnt(14)/amat_sum, &
                    "%", cmap_fp32_error_abs_cnt(15)/amat_sum, "%"
        call cgyro_info(msg)
        write (msg, "(A,4(F6.2,A),4(F6.1,A))") &
                    "Rel cmat_fp32 error rates: ", cmap_fp32_error_rel_cnt(8)/amat_sum, &
                    "%", cmap_fp32_error_rel_cnt(9)/amat_sum, &
                    "%",cmap_fp32_error_rel_cnt(10)/amat_sum, &
                    "%",cmap_fp32_error_rel_cnt(11)/amat_sum, &
                    "%",cmap_fp32_error_rel_cnt(12)/amat_sum, &
                    "%",cmap_fp32_error_rel_cnt(13)/amat_sum, &
                    "%",cmap_fp32_error_rel_cnt(14)/amat_sum, &
                    "%",cmap_fp32_error_rel_cnt(15)/amat_sum, "%"
        call cgyro_info(msg)
        write (msg, "(A,1PE9.2)") &
                    "Abs cmat_fp32 error avg: ", 0.01*cmap_fp32_error_sum(1)/amat_sum
        call cgyro_info(msg)
        ! Printout used for CGYRO paper.
        !write (msg, "(A,1PE9.2)") &
        !            "Abs cmat plain fp32 error avg: ", 0.01*cmap_fp32_error_sum(2)/amat_sum
        !call cgyro_info(msg)
     endif
#if (!defined(OMPGPU)) && defined(_OPENACC)
!$acc wait
#endif
  else
#if defined(OMPGPU)
!$omp target enter data map(to:cmat) if (gpu_bigmem_flag == 1)
#elif defined(_OPENACC)
!$acc enter data copyin(cmat) if (gpu_bigmem_flag == 1)
#endif
  endif

  deallocate(cmat_base2)
  deallocate(cmat_base1)
  deallocate(i_piv)
  deallocate(nu_d)
  deallocate(nu_par)
  deallocate(ctest)
  deallocate(klor_fac)
  deallocate(kdiff_fac)

end subroutine cgyro_init_collision
