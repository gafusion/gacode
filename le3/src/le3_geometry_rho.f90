subroutine le3_geometry_rho

  use le3_globals
  use geo

  implicit none

  integer :: i,j,k,kp
  integer :: nrow
  real :: c_n0,c_m0
  real :: c_dn,c_dm
  real :: x,y
  real :: ycosuv
  real :: ysinuv
  real :: up,pp,pt
  real :: d0, eta, d1, gpp_1, gtt_1, gpt_1
  real :: sum1, sum2, sum3, sum4
  
  real, dimension(:,:), allocatable :: s1a,s1b
  real, dimension(:,:), allocatable :: s2a,s2b
  real, dimension(:,:), allocatable :: a11a,a11b
  real, dimension(:,:), allocatable :: a12a
  real, dimension(:,:), allocatable :: a21a,a21b

  real, dimension(:,:), allocatable :: sys_m
  real, dimension(:), allocatable :: sys_b
  integer, dimension(:), allocatable :: i_piv

  ! Since this is a 2x2 system, define new matsize
  nrow = 2*matsize

  ! Step 1: Generate required functions on (theta,phi) mesh
  allocate(s1a(nt,np))
  allocate(s1b(nt,np))
  allocate(s2a(nt,np))
  allocate(s2b(nt,np))
  allocate(a11a(nt,np))
  allocate(a11b(nt,np))
  allocate(a12a(nt,np))
  allocate(a21a(nt,np))
  allocate(a21b(nt,np))

  ! Global second order quantities
  allocate(dtheta(nt,np))
  allocate(dchi(nt,np))
  allocate(dtheta_t(nt,np))
  allocate(dtheta_p(nt,np))
  allocate(b1(nt,np))
  allocate(gc(nt,np))
  allocate(gct(nt,np))
  allocate(gcp(nt,np))
  allocate(gcc(nt,np))

  do i=1,nt
     do j=1,np

        !-------------------------------------------------------------
        ! DEFINITION BLOCK (reused below)
        !
        x = sqrt(gtt(i,j))
        y = sqrt(drdp(i,j)**2+dzdp(i,j)**2)

        ycosuv = gpt(i,j)/x
        ysinuv = (dzdt(i,j)*drdp(i,j)-drdt(i,j)*dzdp(i,j))/x
        up = (drdt(i,j)*dzdpt(i,j)-dzdt(i,j)*drdpt(i,j))/gtt(i,j)
        pp =          up*ycosuv-chi1p(i,j)/chi1(i,j)*ysinuv 
        pt = (x/rc(i,j))*ycosuv-chi1t(i,j)/chi1(i,j)*ysinuv

        eta = 1.0/rc(i,j)+cosu(i,j)/r(i,j)
        
        c_m0 = iota*x**2+x*ycosuv
        c_n0 = r(i,j)**2+y**2+iota*x*ycosuv

        c_dm = x*up+pt+iota*x**2*2.0/rc(i,j)+chi1(i,j)*iota_p*x**2-c_m0*eta 
        c_dn = 2*r(i,j)*cosu(i,j)+iota*x*up+2*pp+iota*pt+chi1(i,j)*iota_p*x*ycosuv-c_n0*eta

        d0 = g(i,j)
        !-------------------------------------------------------------

        s1a(i,j) = -0.5*d0*beta_star/rmin+(c_dn+iota*c_dm)/(chi1(i,j)*d0)
        s1b(i,j) = ysinuv/(chi1(i,j)*d0)

        s2a(i,j) = c_dm/(chi1(i,j)*d0)
        s2b(i,j) = c_dn/(chi1(i,j)*d0)

        a11a(i,j) = c_m0/d0
        a11b(i,j) = c_n0/d0

        a12a(i,j) = (c_n0+iota*c_m0)/d0

        a21a(i,j) = (-2*c_m0+iota*x**2)/d0
        a21b(i,j) = x**2/d0

     enddo
  enddo

  allocate(sys_m(nrow,nrow))
  allocate(sys_b(nrow))
  allocate(i_piv(nrow))

  ! Initialize arrays
  sys_m(:,:) = 0.0
  sys_b(:)   = 0.0

  do k=1,matsize

     ! Projection basis
     call le3_basis(itype(k),m_indx(k),n_indx(k),bk(:,:),'d0')
     call le3_basis(itype(k),m_indx(k),n_indx(k),bk_t(:,:),'dt')
     call le3_basis(itype(k),m_indx(k),n_indx(k),bk_p(:,:),'dp')

     ! Construct source
     do i=1,nt
        do j=1,np

           ! (1)
           sys_b(k) = sys_b(k) &
                +bk(i,j)*s1a(i,j) &
                +(bk_p(i,j)+iota*bk_t(i,j))*s1b(i,j)

           ! (2)
           sys_b(k+matsize) = sys_b(k+matsize) &
                -bk_p(i,j)*s2a(i,j) &
                +bk_t(i,j)*s2b(i,j)
        enddo
     enddo

     do kp=1,matsize

        ! Expansion basis
        call le3_basis(itype(kp),m_indx(kp),n_indx(kp),bkp(:,:),'d0')
        call le3_basis(itype(kp),m_indx(kp),n_indx(kp),bkp_t(:,:),'dt')
        call le3_basis(itype(kp),m_indx(kp),n_indx(kp),bkp_p(:,:),'dp')

        ! Construct matrix
        do i=1,nt
           do j=1,np

              ! (1,1)
              sys_m(k,kp) = sys_m(k,kp) &
                   +bk(i,j)*a11a(i,j)*bkp_p(i,j) &
                   -bk(i,j)*a11b(i,j)*bkp_t(i,j)

              ! (1,2)
              sys_m(k,kp+matsize) = sys_m(k,kp+matsize) &
                   -bk(i,j)*a12a(i,j)*bkp(i,j)

              ! (2,1)
              sys_m(k+matsize,kp) = sys_m(k+matsize,kp) &
                   -bk_t(i,j)*a21a(i,j)*(bkp_p(i,j)+iota*bkp_t(i,j)) &
                   -bk_t(i,j)*a12a(i,j)*bkp_t(i,j) &
                   -bk_p(i,j)*a21b(i,j)*(bkp_p(i,j)+iota*bkp_t(i,j)) &
                   +bk_p(i,j)*a11a(i,j)*bkp_t(i,j) &
                   -bk_t(i,j)*a11a(i,j)*bkp_p(i,j)

              ! (2,2)
              sys_m(k+matsize,kp+matsize) = sys_m(k+matsize,kp+matsize) &
                   -bk_t(i,j)*a11b(i,j)*bkp(i,j) &
                   +bk_p(i,j)*a11a(i,j)*bkp(i,j)

           enddo
        enddo

     enddo
  enddo

  ! Replace zero-row (due to zero theta-average) with c00=0 for delta_theta.
  sys_m(matsize+1,:) = 0.0
  sys_m(matsize+1,1) = 1.0
  sys_b(matsize) = 0.0

  call DGESV(nrow,1,sys_m,nrow,i_piv,sys_b,nrow,info)

  if (info /= 0) then
     print '(a)', 'ERROR (le3_geometry_rho) Singular matrix in DGESV.'
     stop
  endif

  open(unit=1,file='out.le3.rho',status='replace')
  do k=1,matsize
     write(1,'(2(1pe13.6,1x))') sys_b(k),sys_b(k+matsize)
  enddo
  close(1)

  dtheta  = 0.0
  dchi    = 0.0
  dtheta_t = 0.0
  dtheta_p = 0.0
  do k=1,matsize
     call le3_basis(itype(k),m_indx(k),n_indx(k),bk(:,:),'d0')
     call le3_basis(itype(k),m_indx(k),n_indx(k),bk_t(:,:),'dt')
     call le3_basis(itype(k),m_indx(k),n_indx(k),bk_p(:,:),'dp')
     do j=1,np
        do i=1,nt
           dchi(i,j)     = dchi(i,j)     + sys_b(k+matsize)*bk(i,j)
           dtheta(i,j)   = dtheta(i,j)   + sys_b(k)*bk(i,j)
           dtheta_t(i,j) = dtheta_t(i,j) + sys_b(k)*bk_t(i,j)
           dtheta_p(i,j) = dtheta_p(i,j) + sys_b(k)*bk_p(i,j)

        enddo
     enddo
  enddo

  ! Compute additional quantities required for GK drifts:
  !
  ! b1   -> B1/Bs
  ! g_c -> g_cp + i g_ct

  do j=1,np
     do i=1,nt

        !-------------------------------------------------------------
        ! DEFINITION BLOCK (reused above)
        !
        x = sqrt(gtt(i,j))
        y = sqrt(drdp(i,j)**2+dzdp(i,j)**2)

        ycosuv = gpt(i,j)/x
        ysinuv = (dzdt(i,j)*drdp(i,j)-drdt(i,j)*dzdp(i,j))/x
        up = (drdt(i,j)*dzdpt(i,j)-dzdt(i,j)*drdpt(i,j))/gtt(i,j)
        pp =          up*ycosuv-chi1p(i,j)/chi1(i,j)*ysinuv 
        pt = (x/rc(i,j))*ycosuv-chi1t(i,j)/chi1(i,j)*ysinuv

        eta = 1.0/rc(i,j)+cosu(i,j)/r(i,j)

        c_m0 = iota*x**2+x*ycosuv
        c_n0 = r(i,j)**2+y**2+iota*x*ycosuv

        c_dm = x*up+pt+iota*x**2*2.0/rc(i,j)+chi1(i,j)*iota_p*x**2-c_m0*eta 
        c_dn = 2*r(i,j)*cosu(i,j)+iota*x*up+2*pp+iota*pt+chi1(i,j)*iota_p*x*ycosuv-c_n0*eta

        ! D0 = sqrt(g_s)
        d0 = g(i,j)

        ! D1/D0
        d1 = eta - chi1(i,j)*(dchi(i,j)+dtheta_t(i,j))

        gpp_1 = 2.0*(r(i,j)*cosu(i,j) - chi1p(i,j)/chi1(i,j)*ysinuv &
             + (up - chi1(i,j)*x*dtheta_p(i,j))*ycosuv)

        gtt_1 = 2.0 * x**2 * (1.0/rc(i,j)- chi1(i,j)*dtheta_t(i,j))

        gpt_1 = x*(up - chi1(i,j)*x*dtheta_p(i,j)) &
             + (x/rc(i,j) - chi1(i,j)*x*dtheta_t(i,j))*ycosuv &
             - chi1t(i,j)/chi1(i,j)*ysinuv
        
        !-------------------------------------------------------------

        ! B1/Bs
        b1(i,j) = -d1 &
             + 0.5/(c_n0 + iota*c_m0) * (gpp_1 + iota**2 * gtt_1 + 2.0*iota*gpt_1 &
             + chi1(i,j) * iota_p * (2*iota*gtt(i,j) + gpt(i,j)) )
        
        ! g_cp + i g_ct
        gc(i,j) = ysinuv/chi1(i,j)-c_m0*dtheta(i,j)

        ! g_ct
        gct(i,j) = -dtheta(i,j) * x**2
        
        ! g_cp
        gcp(i,j) = ysinuv/chi1(i,j) - dtheta(i,j)*x*ycosuv
        
        ! g_cc
        gcc(i,j) = 1.0/chi1(i,j)**2 + dtheta(i,j)**2 * x**2

     enddo
  enddo

  !do i=1,nt
  !   j=1
  !   x = sqrt(gtt(i,j))
     !print *, b1(i,j), -1.0/(r(i,j)**2 + iota**2 * x**2) &
     !     * (r(i,j)*cosu(i,j) + iota**2*x**2/rc(i,j)), &
     !     chi1(i,j)*r(i,j)/sqrt(gtt(i,j)), &
     !     r(i,j)/sqrt(gtt(i,j))&
     !  *(-1/rc(i,j)+cosu(i,j)/r(i,j)+chi1(i,j)*(dchi(i,j)+dtheta_t(i,j)))/iota
  !enddo
  
  !--------------------------------------------------------------------
  ! Check with GEO result
  !
  !if (equilibrium_model == 0) then
  !   GEO_model_in = 0
  !   GEO_rmin_in = rmin
  !   GEO_rmaj_in = rmaj
  !   GEO_drmaj_in = shift
  !   GEO_zmag_in = zmag
  !   GEO_dzmag_in = dzmag
  !   GEO_q_in = q
  !   GEO_s_in = s
  !   GEO_kappa_in = kappa
  !   GEO_s_kappa_in = s_kappa
  !   GEO_delta_in = delta
  !   GEO_s_delta_in = s_delta
  !   GEO_zeta_in = zeta
  !   GEO_s_zeta_in = s_zeta
  !   GEO_beta_star_in = beta_star
  !   call GEO_alloc(1)
  !   call GEO_do()
  !call GEO_interp(nt,tb(:,1),.false.)
  !do i=1,nt
  !   print *, tb(i,1), btor(i,1), GEO_bt(i), bpol(i,1), GEO_bp(i)
  !enddo
  !close(1)
  !call GEO_alloc(0)
  !endif
  !--------------------------------------------------------------------
  
  deallocate(sys_m)
  deallocate(sys_b)
  deallocate(i_piv)

end subroutine le3_geometry_rho
