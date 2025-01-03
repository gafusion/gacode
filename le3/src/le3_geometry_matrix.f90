subroutine le3_geometry_matrix

  use le3_globals
  implicit none
  integer :: i,j,k,its,ips,kt,kp
  
  real, dimension(:), allocatable :: vec_vdriftx
  real, dimension(:), allocatable :: vec_flux
  real, dimension(:), allocatable :: vec_upar
  real, dimension(:), allocatable :: vec_uparB 
  real, dimension(:), allocatable :: vec_fsa
  real, dimension(:), allocatable :: vec_bmag
  real, dimension(:), allocatable :: vec_thetabar
  real, dimension(:), allocatable :: vec_ntv

  real :: bsq_avg
  real :: J_boozer, I_boozer
  real, dimension(:), allocatable :: vec_boozer_p, vec_boozer_t, vec_boozer
  real, dimension(:,:), allocatable :: mat_boozer, t_boozer, p_boozer
  real, dimension(:,:), allocatable :: xtemp
  real, dimension(:), allocatable :: tmap
  real :: xsum
  integer, dimension(:), allocatable :: i_piv

  ! anorm* bhat dot grad = bdotgrad * (iota d/dt + d/dp)  
  bdotgrad(:,:) = 1.0/(bmag * g)

  ! anorm * (bhat dot grad B)/B
  bdotgradB_overB(:,:) = bdotgrad * (iota * dbdt + dbdp) / bmag

  ! bhat cross grad B dot grad chi / (B^2 * r)
  vdrift_x(:,:) = 1/(rmin*bmag* g**2) &
       * (-dbdt * (gpp + iota * gpt) + dbdp * (gpt + iota*gtt)) / bmag**2

  ! bhat cross grad B dot grad theta / B^2
  vdrift_dt(:,:) = 1/g**2 &
       * ((b1*bmag/chi1 - dtheta*dbdt)*(gpp + iota * gpt) &
       - dbdp*gc) / bmag**3

  ! bhat cross grad B dot grad phi / B^2
  vdrift_dp(:,:) = 1/g**2 &
       * (-(b1*bmag/chi1 - dtheta*dbdt)*(gpt + iota * gtt) &
       + dbdt*gc) / bmag**3
 
  ! (a/cs)*vdrift_gradB = 1/(cs*anorm*Omega_ca_unit)*(vperp^2/2+vpar^2)*[(1) d/d(r/a) + (2) (i k_theta a)]
  ! Remap theta: 0..2pi -> -pi..pi

  allocate(xtemp(nt,np))
  allocate(tmap(nt))
  do i=1,nt/2
     k=nt/2+i
     tmap(i) = t(k)-2*pi
     tmap(k)   = t(i)
  enddo

  vdrift_gk(:,:,1) = vdrift_x(:,:)
  call remap_theta(vdrift_gk(:,:,1),xtemp(:,:),0.0)
  vdrift_gk(:,:,1) = xtemp(:,:)

  vdrift_gk(:,:,2)  = (-iota*vdrift_dp(:,:) + vdrift_dt(:,:))*rmin
  call remap_theta(vdrift_gk(:,:,2),xtemp(:,:),0.0)
  vdrift_gk(:,:,2) = xtemp(:,:)
  do i=1,nt
     do j=1,np
        vdrift_gk(i,j,2) = vdrift_gk(i,j,2) &
             - vdrift_gk(i,j,1) * rmin**2 *(iota_p/iota)*tmap(i)
     enddo
  enddo

  ! (a/cs)*vdrift_gradp =  1/(cs*anorm*Omega_ca_unit)*(vpar^2) [(3) (i k_theta a)]
  vdrift_gk(:,:,3) = (-0.5*beta_star)/bmag(:,:)**2
  call remap_theta(vdrift_gk(:,:,3),xtemp(:,:),0.0)
  vdrift_gk(:,:,3) = xtemp(:,:)
  
  ! grad_perpsq: (1) d^2/d(r/a)^2 + (2) (k_theta a)^2 + (3) (i k_theta a) d/d(r/a)
  grad_cc =  1.0/g**2 * (gtt * gpp - gpt**2)
  grad_tt =  1.0/g**2 * (gpp * gcc - gcp**2)
  grad_pp =  1.0/g**2 * (gtt * gcc - gct**2)
  grad_pt =  1.0/g**2 * (gcp * gct - gpt*gcc)
  grad_ct =  1.0/g**2 * (gpt * gcp - gpp*gct)
  grad_cp =  1.0/g**2 * (gpt * gct - gtt*gcp)
  
  grad_perpsq_gk(:,:,1) = 1.0/rmin**2 * grad_cc(:,:)
  call remap_theta(grad_perpsq_gk(:,:,1),xtemp(:,:),0.0)
  grad_perpsq_gk(:,:,1) = xtemp(:,:)
  
  grad_perpsq_gk(:,:,2) = grad_pp + 1.0/iota**2*grad_tt - 2.0/iota*grad_pt
  call remap_theta(grad_perpsq_gk(:,:,2),xtemp(:,:),0.0)
  grad_perpsq_gk(:,:,2) = xtemp(:,:)
  call remap_theta(grad_cc(:,:),xtemp(:,:),0.0)
  do i=1,nt
     do j=1,np
        grad_perpsq_gk(i,j,2) = grad_perpsq_gk(i,j,2) + (iota_p/iota**2*tmap(i))**2 * xtemp(i,j)
     enddo
  enddo
  call remap_theta(grad_cp(:,:),xtemp(:,:),0.0)
    do i=1,nt
     do j=1,np
        grad_perpsq_gk(i,j,2) = grad_perpsq_gk(i,j,2) + 2.0*iota_p/iota**2*tmap(i)*xtemp(i,j)
     enddo
  enddo
  call remap_theta(grad_ct(:,:),xtemp(:,:),0.0)
    do i=1,nt
     do j=1,np
        grad_perpsq_gk(i,j,2) = grad_perpsq_gk(i,j,2) - 2.0*iota_p/iota**3*tmap(i)*xtemp(i,j)
     enddo
  enddo
  grad_perpsq_gk(:,:,2) = grad_perpsq_gk(:,:,2)*(iota*rmin)**2
  
  grad_perpsq_gk(:,:,3) = grad_cp - 1.0/iota*grad_ct
  call remap_theta(grad_perpsq_gk(:,:,3),xtemp(:,:),0.0)
  grad_perpsq_gk(:,:,3) = xtemp(:,:)
  call remap_theta(grad_cc(:,:),xtemp(:,:),0.0)
  do i=1,nt
     do j=1,np
        grad_perpsq_gk(i,j,3) =  grad_perpsq_gk(i,j,3) + iota_p/iota**2*tmap(i)* xtemp(i,j)
     enddo
  enddo
  grad_perpsq_gk(:,:,3) =  grad_perpsq_gk(:,:,3)*(-2.0*iota)
  
  ! flux-surface d volume / dr
  vprime = 0.0
  do i=1,nt
     do j=1,np
        vprime = vprime + g(i,j)
     enddo
  enddo
  vprime = vprime / (nt*np)

  ! -bhat cross grad f dot grad r / B 
  ! in the "DKES" formulation, replace 1/B with B/<B^2>
  ! < B^2>
  bsq_avg      = 0.0
  do i=1,nt
     do j=1,np
        bsq_avg = bsq_avg + g(i,j) * bmag(i,j)**2
     enddo
  enddo
  bsq_avg      = bsq_avg / (nt*np) / vprime
  vexb_dt(:,:) = -1/(rmin*bmag * g**2) &
       * (-(gpp + iota * gpt)) * bmag/bsq_avg
  vexb_dp(:,:) = -1/(rmin*bmag * g**2) & 
       * (gpt + iota*gtt) * bmag/bsq_avg

  ! construct the geo collocation matices for NEO-3D

  allocate(mat_stream_dt(matsize,matsize))
  allocate(mat_stream_dp(matsize,matsize))
  allocate(mat_trap(matsize,matsize))
  allocate(mat_coll(matsize,matsize))
  allocate(mat_vexb_dt(matsize,matsize))
  allocate(mat_vexb_dp(matsize,matsize))
  allocate(mat_vdrift_dt(matsize,matsize))
  allocate(mat_vdrift_dp(matsize,matsize))
  allocate(vec_vdriftx(matsize))
  allocate(vec_flux(matsize))
  allocate(vec_upar(matsize))
  allocate(vec_uparB(matsize))
  allocate(vec_fsa(matsize))
  allocate(vec_bmag(matsize))
  allocate(vec_thetabar(matsize))
  allocate(vec_ntv(matsize))

  mat_stream_dt(:,:) = 0.0
  mat_stream_dp(:,:) = 0.0
  mat_trap(:,:)      = 0.0
  mat_coll(:,:)      = 0.0
  mat_vexb_dt(:,:)   = 0.0
  mat_vexb_dp(:,:)   = 0.0
  mat_vdrift_dt(:,:) = 0.0
  mat_vdrift_dp(:,:) = 0.0
  vec_vdriftx(:)  = 0.0
  vec_flux(:)     = 0.0
  vec_upar(:)     = 0.0
  vec_uparB(:)    = 0.0
  vec_fsa(:)      = 0.0
  vec_bmag(:)     = 0.0
  vec_thetabar(:) = 0.0
  vec_ntv(:)      = 0.0

  do i=1,matsize
     call le3_basis(itype(i),m_indx(i),n_indx(i),bk,'d0')
     do j=1,matsize
        call le3_basis(itype(j),m_indx(j),n_indx(j),bkp(:,:),'d0')
        call le3_basis(itype(j),m_indx(j),n_indx(j),bkp_t(:,:),'dt')
        call le3_basis(itype(j),m_indx(j),n_indx(j),bkp_p(:,:),'dp')
        do kt=1,nt
           do kp=1,np
              mat_trap(i,j) = mat_trap(i,j) &
                   + bk(kt,kp) * bkp(kt,kp) * bdotgradB_overB(kt,kp) 
              mat_stream_dt(i,j) = mat_stream_dt(i,j) &
                   + bk(kt,kp) * bkp_t(kt,kp) * bdotgrad(kt,kp) * iota
              mat_stream_dp(i,j) = mat_stream_dp(i,j) &
                   + bk(kt,kp) * bkp_p(kt,kp) * bdotgrad(kt,kp) 
              mat_coll(i,j) = mat_coll(i,j) &
                   + bk(kt,kp) * bkp(kt,kp)
              mat_vexb_dt(i,j) = mat_vexb_dt(i,j) &
                   + bk(kt,kp) * bkp_t(kt,kp) * vexb_dt(kt,kp)
              mat_vexb_dp(i,j) = mat_vexb_dp(i,j) &
                   + bk(kt,kp) * bkp_p(kt,kp) * vexb_dp(kt,kp)
              mat_vdrift_dt(i,j) = mat_vdrift_dt(i,j) &
                   + bk(kt,kp) * bkp_t(kt,kp) * vdrift_dt(kt,kp)
              mat_vdrift_dp(i,j) = mat_vdrift_dp(i,j) &
                   + bk(kt,kp) * bkp_p(kt,kp) * vdrift_dp(kt,kp)
           enddo
        enddo
     enddo
     do kt=1,nt
        do kp=1,np
           vec_thetabar(i) = vec_thetabar(i) &
                + bk(kt,kp) * (tb(kt,kp)-t(kt))
           vec_vdriftx(i) = vec_vdriftx(i) &
                + bk(kt,kp) * vdrift_x(kt,kp)
           vec_flux(i)    = vec_flux(i) &
                + bk(kt,kp) * vdrift_x(kt,kp) * g(kt,kp)
           vec_upar(i)    = vec_upar(i) &
                + bk(kt,kp) 
           vec_uparB(i)   = vec_uparB(i) &
                + bk(kt,kp) * bmag(kt,kp) * g(kt,kp)
           vec_fsa(i)   = vec_fsa(i) &
                + bk(kt,kp) * g(kt,kp)
           vec_bmag(i)   = vec_bmag(i) &
                + bk(kt,kp) * bmag(kt,kp) 
           vec_ntv(i)   = vec_ntv(i) &
                + bk(kt,kp) * dgdp(kt,kp)
        enddo
     enddo
  enddo

  mat_stream_dt(:,:) = mat_stream_dt(:,:) / (nt*np)
  mat_stream_dp(:,:) = mat_stream_dp(:,:) / (nt*np)
  mat_trap(:,:)      = mat_trap(:,:)      / (nt*np)
  mat_coll(:,:)      = mat_coll(:,:)      / (nt*np)
  mat_vexb_dt(:,:)   = mat_vexb_dt(:,:)   / (nt*np)
  mat_vexb_dp(:,:)   = mat_vexb_dp(:,:)   / (nt*np)
  mat_vdrift_dt(:,:) = mat_vdrift_dt(:,:) / (nt*np)
  mat_vdrift_dp(:,:) = mat_vdrift_dp(:,:) / (nt*np)

  open(unit=1,file='out.le3.geomatrix',status='replace')
  do i=1,matsize
     do j=1, matsize
        write (1,'(e16.8)',advance='no') mat_trap(i,j)
        write (1,'(e16.8)',advance='no') mat_stream_dt(i,j)
        write (1,'(e16.8)',advance='no') mat_stream_dp(i,j)
        write (1,'(e16.8)',advance='no') mat_coll(i,j)
        write (1,'(e16.8)',advance='no') mat_vexb_dt(i,j)
        write (1,'(e16.8)',advance='no') mat_vexb_dp(i,j)
        write (1,'(e16.8)',advance='no') mat_vdrift_dt(i,j)
        write (1,'(e16.8)',advance='no') mat_vdrift_dp(i,j)
        write (1,*)
     enddo
  enddo
  close(1)

  ! Construct the geo vectors for NEO-3D

  vec_thetabar(:) = vec_thetabar(:) / (nt*np)
  vec_vdriftx(:)  = vec_vdriftx(:) / (nt*np)
  vec_flux(:)     = vec_flux(:)  / (nt*np) / vprime
  vec_uparB(:)    = vec_uparB(:) / (nt*np) / vprime
  vec_upar(:)     = vec_upar(:)  / (nt*np) * (2*pi)**2
  vec_fsa(:)      = vec_fsa (:)  / (nt*np) / vprime
  vec_bmag(:)     = vec_bmag (:)  / (nt*np)
  vec_ntv(:)      = vec_ntv(:) / (nt*np) / vprime

  open(unit=1,file='out.le3.geovector',status='replace')
  do i=1,matsize
     write (1,'(e16.8)',advance='no') vec_thetabar(i)
     write (1,'(e16.8)',advance='no') vec_vdriftx(i)
     write (1,'(e16.8)',advance='no') vec_flux(i)
     write (1,'(e16.8)',advance='no') vec_uparB(i)
     write (1,'(e16.8)',advance='no') vec_upar(i)
     write (1,'(e16.8)',advance='no') vec_fsa(i)
     write (1,'(e16.8)',advance='no') vec_bmag(i)
     write (1,'(e16.8)',advance='no') vec_ntv(i)
     write (1,*)
  enddo
  close(1)

  open(unit=1,file='out.le3.geoscalar',status='replace')
  write (1,'(i3)') nts
  write (1,'(i3)') nps
  write (1,'(i3)') matsize
  write (1,'(i3)') indx_c00
  write (1,'(e16.8)') bsq_avg
  close(1)

  open(unit=1,file='out.le3.geoindx',status='replace')
  do i=1,matsize
     write (1,'(i2,i2,i2)') itype(i),m_indx(i), n_indx(i)
  enddo
  close(1)
  
  ! Map to Boozer coordinates

  ! J: B ~ J(psi) grad phi
  J_boozer = 0.0
  do its=1,nt
     do ips=1,np
        J_boozer = J_boozer + 1/g(its,ips)*(gpp(its,ips)+iota*gpt(its,ips)) 
     enddo
  enddo
  J_boozer = J_boozer/(nt*np)

  ! I: B ~ I(psi) grad theta
  I_boozer = 0.0
  do ips=1,np
     do its=1,nt
        I_boozer = I_boozer + 1/g(its,ips)*(gpt(its,ips)+iota*gtt(its,ips)) 
     enddo
  enddo
  I_boozer = I_boozer/(nt*np)

  ! Compute omega: Solve domega/dphi and domega/dtheta equations for omega
  ! using fourier expansion

  allocate(vec_boozer_p(matsize))
  vec_boozer_p(:) = 0.0
  do i=1,matsize
     call le3_basis(itype(i),m_indx(i),n_indx(i),bk,'dp')
     do kt=1,nt
        do kp=1,np
           if(n_indx(i) /= 0) then
              vec_boozer_p(i) = vec_boozer_p(i) &
                   + bk(kt,kp) * 1.0/(J_boozer + iota*I_boozer) * (J_boozer &
                   - 1/g(kt,kp)*(gpp(kt,kp)+iota*gpt(kt,kp))) &
                   / n_indx(i)**2
           endif
        enddo
     enddo
  enddo
  vec_boozer_p(:) = vec_boozer_p(:)/(nt*np)

  allocate(vec_boozer_t(matsize))
  vec_boozer_t(:) = 0.0
  do i=1,matsize
     call le3_basis(itype(i),m_indx(i),n_indx(i),bk,'dt')
     do kt=1,nt
        do kp=1,np
           if(m_indx(i) /= 0) then
              vec_boozer_t(i) = vec_boozer_t(i) &
                   + bk(kt,kp) * 1.0/(J_boozer + iota*I_boozer) * (I_boozer &
                   - 1/g(kt,kp)*(gpt(kt,kp)+iota*gtt(kt,kp))) &
                   / m_indx(i)**2
           endif
        enddo
     enddo
  enddo
  vec_boozer_t(:) = vec_boozer_t(:)/(nt*np)

  allocate(vec_boozer(matsize))
  do i=1,matsize
     if(itype(i) == 1 .and. n_indx(i) == 0) then
        vec_boozer(i) = vec_boozer_t(i)
     else if(itype(i) == 3 .and. n_indx(i) == 0) then
        vec_boozer(i) = vec_boozer_t(i)
     else if(itype(i) == 3 .and. m_indx(i) == 0) then
        vec_boozer(i) = vec_boozer_p(i)
     else if(itype(i) == 4 .and. m_indx(i) == 0) then
        vec_boozer(i) = vec_boozer_p(i)
     else
        vec_boozer(i) = 0.5*(vec_boozer_p(i)+vec_boozer_t(i))
     endif
  enddo

  ! Form theta_b = theta - iota*omega and phi_b = phi - omega

  allocate(t_boozer(nt,np))
  allocate(p_boozer(nt,np))
  do kt=1,nt
     do kp=1,np
        xsum=0.0
        do i=1, matsize
           call le3_basis(itype(i),m_indx(i),n_indx(i),bk,'d0')
           xsum = xsum + bk(kt,kp)*vec_boozer(i)
        enddo
        t_boozer(kt,kp) = t(kt) - iota*xsum
        p_boozer(kt,kp) = p(kp) - xsum
     enddo
  enddo

  ! Map B(theta,phi) fourier coeffs to B(theta_b,phi_b)

  allocate(mat_boozer(matsize,matsize))

  ! RHS
  mat_boozer(:,:) = 0.0
  do i=1,matsize
     call le3_basis(itype(i),m_indx(i),n_indx(i),bk,'d0')
     do j=1,matsize
        call le3_basis(itype(j),m_indx(j),n_indx(j),bkp,'d0')
        do kt=1,nt
           do kp=1,np
              mat_boozer(i,j) = mat_boozer(i,j) + bk(kt,kp) * bkp(kt,kp)
           enddo
        enddo
     enddo
  enddo
  mat_boozer(:,:) = mat_boozer(:,:)/(nt*np)
  vec_boozer_t(:) = vec_bmag(:)
  call DGEMV ('N', matsize, matsize,1.0,mat_boozer(:,:),matsize, &
       vec_boozer_t,1,0.0,vec_boozer,1)

  mat_boozer(:,:) = 0.0
  do i=1,matsize
     call le3_basis(itype(i),m_indx(i),n_indx(i),bk,'d0')
     do j=1,matsize
        ! EAB not yet implemented
        bkp(:,:) = 0.0
        if(itype(j) == 1) then
           do kt=1,nt
              do kp=1,np
                 bkp(kt,kp) = sin(m_indx(j)*t_boozer(kt,kp)) &
                      * cos(n_indx(j)*p_boozer(kt,kp))
              enddo
           enddo
        else if(itype(j) == 2) then
           do kt=1,nt
              do kp=1,np
                 bkp(kt,kp) = sin(m_indx(j)*t_boozer(kt,kp)) &
                      * sin(n_indx(j)*p_boozer(kt,kp))
              enddo
           enddo
        else if(itype(j) == 3) then
           do kt=1,nt
              do kp=1,np
                 bkp(kt,kp) = cos(m_indx(j)*t_boozer(kt,kp)) &
                      * cos(n_indx(j)*p_boozer(kt,kp))
              enddo
           enddo
        else if(itype(j) == 4) then
           do kt=1,nt
              do kp=1,np
                 bkp(kt,kp) = cos(m_indx(j)*t_boozer(kt,kp)) &
                      * sin(n_indx(j)*p_boozer(kt,kp))
              enddo
           enddo
        endif
        !
        do kt=1,nt
           do kp=1,np
              mat_boozer(i,j) = mat_boozer(i,j) + bk(kt,kp) * bkp(kt,kp)
           enddo
        enddo
        !
     enddo
  enddo
  mat_boozer(:,:) = mat_boozer(:,:)/(nt*np)
  
  ! Solve the matrix problem
  allocate(i_piv(matsize))
  call DGETRF(matsize,matsize,mat_boozer(:,:),matsize,i_piv,info)
  call DGETRS('N',matsize,1,mat_boozer(:,:),matsize,i_piv,vec_boozer(:),&
       matsize,info)
  deallocate(i_piv)

  open(unit=1,file='out.le3.geoboozer',status='replace')
  write (1,'(e16.8)') J_boozer
  write (1,'(e16.8)') I_boozer
  do i=1,matsize
     if(itype(i) == 1) then
        write (1,'(e16.8,e16.8,i2,i2,i2)') vec_boozer(i), vec_bmag(i), &
          itype(i),m_indx(i), n_indx(i)
     endif
  enddo
  do i=1,matsize
     if(itype(i) == 2) then
        write (1,'(e16.8,e16.8,i2,i2,i2)') vec_boozer(i), vec_bmag(i), &
          itype(i),m_indx(i), n_indx(i)
     endif
  enddo
  do i=1,matsize
     if(itype(i) == 3) then
        write (1,'(e16.8,e16.8,i2,i2,i2)') vec_boozer(i), vec_bmag(i), &
          itype(i),m_indx(i), n_indx(i)
     endif
  enddo
  do i=1,matsize
     if(itype(i) == 4) then
        write (1,'(e16.8,e16.8,i2,i2,i2)') vec_boozer(i), vec_bmag(i), &
          itype(i),m_indx(i), n_indx(i)
     endif
  enddo
  close(1)

  deallocate(mat_boozer)
  deallocate(vec_boozer)
  deallocate(vec_boozer_p)
  deallocate(vec_boozer_t)

  call le3_compute_theory

  ! Write out quantitites for GK-3D
  ! remap theta: 0..2pi -> -pi..pi
  
  open(unit=1,file='out.le3.geogk',status='replace')
  write (1,'(i4)') nt
  write (1,'(i4)') np
  write (1,'(e16.8)') tmap
  write (1,'(e16.8)') p

  ! theta_bar(theta,phi) and derivatives
  call remap_theta(tb(:,:),xtemp(:,:),1.0)
  write (1,'(e16.8)') xtemp(:,:)
  call remap_theta(dtbdt(:,:),xtemp(:,:),0.0)
  write (1,'(e16.8)') xtemp(:,:)
  call remap_theta(dtbdp(:,:),xtemp(:,:),0.0)
  write (1,'(e16.8)') xtemp(:,:)
  
  ! B/Bunit
  call remap_theta(bmag(:,:),xtemp(:,:),0.0)
  write (1,'(e16.8)') xtemp(:,:)
  
  ! anorm*(bhat dot grad) = bdotgrad * (iota d/dt + d/dp) = (1) d/dt + (2) d/dp
  call remap_theta(bdotgrad(:,:),xtemp(:,:),0.0)
  write (1,'(e16.8)') xtemp(:,:)*iota
  write (1,'(e16.8)') xtemp(:,:)
  
  ! anorm*(bhat dot grad B)/B
  call remap_theta(bdotgradB_overB(:,:),xtemp(:,:),0.0)
  write (1,'(e16.8)') xtemp(:,:)
  
  ! (a/cs)*vdrift_gradB = 1/(cs*anorm*Omega_ca_unit)*(vperp^2/2+vpar^2)*[(1) d/d(r/a) + (2) (i k_theta a)]
  write (1,'(e16.8)') vdrift_gk(:,:,1)
  write (1,'(e16.8)') vdrift_gk(:,:,2)
  
  ! (a/cs)*vdrift_gradp =  1/(cs*anorm*Omega_ca_unit)*(par^2)*(i k_theta a) [(3) (i ktheta a)]
  write (1,'(e16.8)') vdrift_gk(:,:,3)
  
  ! grad_perpsq: (1) d^2/d(r/a)^2 + (2) (k_theta a)^2 + (3) (i k_theta a) d/d(r/a)
  write (1,'(e16.8)') grad_perpsq_gk(:,:,1)
  write (1,'(e16.8)') grad_perpsq_gk(:,:,2)
  write (1,'(e16.8)') grad_perpsq_gk(:,:,3)

  call remap_theta(b1(:,:),xtemp(:,:),0.0)
  write (1,'(e16.8)') xtemp(:,:)

  call remap_theta(chi1(:,:),xtemp(:,:),0.0)
  write (1,'(e16.8)') xtemp(:,:)

  call remap_theta(dchi(:,:),xtemp(:,:),0.0)
  write (1,'(e16.8)') xtemp(:,:)

  call remap_theta(dtheta(:,:),xtemp(:,:),0.0)
  write (1,'(e16.8)') xtemp(:,:)

  call remap_theta(dtheta_t(:,:),xtemp(:,:),0.0)
  write (1,'(e16.8)') xtemp(:,:)
  
  close(1)
  deallocate(xtemp)
  deallocate(tmap)
  
  deallocate(itype)
  deallocate(m_indx)
  deallocate(n_indx)

  deallocate(gc)
  deallocate(b1)

  deallocate(g)
  deallocate(gpp)
  deallocate(gtt)
  deallocate(gpt)
  deallocate(gct)
  deallocate(gcp)
  deallocate(gcc)
  deallocate(cosu)
  deallocate(btor)
  deallocate(bpol)
  deallocate(bmag)
  deallocate(dbdt)
  deallocate(dbdp)
  deallocate(dgdp)
  deallocate(deriv_t)
  deallocate(deriv_p)

  deallocate(bdotgrad)
  deallocate(bdotgradB_overB)
  deallocate(vdrift_x)
  deallocate(vexb_dt)
  deallocate(vexb_dp)
  deallocate(vdrift_dt)
  deallocate(vdrift_dp)
  deallocate(vdrift_gk)
  deallocate(grad_perpsq_gk)
  deallocate(grad_cc)
  deallocate(grad_tt)
  deallocate(grad_pp)
  deallocate(grad_pt)
  deallocate(grad_ct)
  deallocate(grad_cp)
  deallocate(mat_stream_dt)
  deallocate(mat_stream_dp)
  deallocate(mat_trap)
  deallocate(mat_coll)
  deallocate(mat_vexb_dt)
  deallocate(mat_vexb_dp)
  deallocate(mat_vdrift_dt)
  deallocate(mat_vdrift_dp)
  deallocate(vec_vdriftx)
  deallocate(vec_flux)
  deallocate(vec_upar)
  deallocate(vec_uparB)
  deallocate(vec_fsa)
  deallocate(vec_bmag)
  deallocate(vec_thetabar)
  deallocate(vec_ntv)

end subroutine le3_geometry_matrix

subroutine remap_theta(xold,xnew,fac)
  use le3_globals
  implicit none
  real, intent(inout), dimension(nt,np) :: xold,xnew
  real, intent(in) :: fac
  integer :: i,j,k
  do i=1,nt/2
     k=nt/2+i
     do j=1,np
        xnew(i,j)   = xold(k,j) - fac*2*pi
        xnew(k,j)   = xold(i,j)
     enddo
  enddo
end subroutine remap_theta
