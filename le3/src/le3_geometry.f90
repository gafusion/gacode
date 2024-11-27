subroutine le3_geometry

  use le3_globals

  implicit none
  integer :: i,j,k,its,ips,ip

  print '(a,1pe12.5)','INFO: (le3) Root accuracy ->',sum(abs(yfunc))/size(yfunc)

  if (nps > 0) then
     print 30,'a_{mn}','b_{mn}','c_{mn}','d_{mn}'
     do ips=0,nps
        do its=0,nts
           print 20,its,ips,as(its,ips),bs(its,ips),cs(its,ips),ds(its,ips)
        enddo
     enddo
  else
     print 30,'a_{m}','c_{m}'
     ips=0
     do its=0,nts
        print 10,its,as(its,ips),cs(its,ips)
     enddo
  endif

  allocate(g(nt,np))
  allocate(gpp(nt,np))
  allocate(gtt(nt,np))
  allocate(gpt(nt,np))
  allocate(cosu(nt,np))
  allocate(chi1(nt,np))
  allocate(chi1p(nt,np))
  allocate(chi1t(nt,np))
  allocate(drdpt(nt,np))
  allocate(dzdpt(nt,np))
  allocate(drdt(nt,np))
  allocate(drdp(nt,np))
  allocate(dzdt(nt,np))
  allocate(dzdp(nt,np))
  allocate(btor(nt,np))
  allocate(bpol(nt,np))
  allocate(bmag(nt,np))
  allocate(dbdt(nt,np))
  allocate(dbdp(nt,np))
  allocate(bdotgrad(nt,np))
  allocate(bdotgradB_overB(nt,np))
  allocate(vdrift_x(nt,np))
  allocate(vdrift_dt(nt,np))
  allocate(vdrift_dp(nt,np))
  allocate(vdrift_gk(nt,np,3))
  allocate(grad_perpsq_gk(nt,np,3))
  allocate(grad_cc(nt,np))
  allocate(grad_tt(nt,np))
  allocate(grad_pp(nt,np))
  allocate(grad_pt(nt,np))
  allocate(grad_ct(nt,np))
  allocate(grad_cp(nt,np))
  allocate(vexb_dt(nt,np))
  allocate(vexb_dp(nt,np))
  allocate(dgdp(nt,np))

  do i=1,nt
     tb(i,:) = t(i)
  enddo
  dtbdt(:,:) = 1.0
  dtbdp(:,:) = 0.0

  do j=1,np
     do i=1,nt

        do ips=0,nps
           do its=0,nts

              tb(i,j) = tb(i,j) + &
                   sinm(i,its)*(bs(its,ips)*sinn(j,ips)+as(its,ips)*cosn(j,ips)) &
                   +cosm(i,its)*(ds(its,ips)*sinn(j,ips)+cs(its,ips)*cosn(j,ips))

              dtbdt(i,j) = dtbdt(i,j) + &
                   its*cosm(i,its)*(bs(its,ips)*sinn(j,ips)+as(its,ips)*cosn(j,ips)) &
                   -its*sinm(i,its)*(ds(its,ips)*sinn(j,ips)+cs(its,ips)*cosn(j,ips))

              dtbdp(i,j) = dtbdp(i,j) + &
                   ips*sinm(i,its)*(bs(its,ips)*cosn(j,ips)-as(its,ips)*sinn(j,ips)) &
                   +ips*cosm(i,its)*(ds(its,ips)*cosn(j,ips)-cs(its,ips)*sinn(j,ips))

           enddo
        enddo

        call le3_rz(tb(i,j),&
             p(j),&
             r(i,j),&
             z(i,j),&
             drdtb(i,j),&
             drdpb(i,j),&
             dzdtb(i,j),&
             dzdpb(i,j),&
             jac(i,j),&
             rc(i,j))

     enddo
  enddo

  ! sqrt(g)
  g(:,:) = dtbdt*jac/rmin

  ! R_t
  drdt(:,:) = drdtb*dtbdt

  ! Z_t
  dzdt(:,:) = dzdtb*dtbdt

  ! R_p
  drdp(:,:) = drdpb+drdtb*dtbdp

  ! Z_p
  dzdp(:,:) = dzdpb+dzdtb*dtbdp

  ! g_pp = R^2 + R_p^2 + Z_p^2
  gpp(:,:) = r**2+drdp**2+dzdp**2

  ! g_tt = R_t^2 + Z_t^2
  gtt(:,:) = drdt**2+dzdt**2

  ! g_pt = R_t R_p + Z_t Z_p
  gpt(:,:) = drdt*drdp+dzdt*dzdp

  ! cos(u) = Z_t/sqrt(g_tt)
  cosu(:,:) = dzdt/sqrt(gtt)

  btor(:,:) = 1.0/(r*g)

  bpol(:,:) = 1.0/g*(gpt/sqrt(gtt)+iota*sqrt(gtt))

  bmag(:,:) = 1.0/g*sqrt(gpp+2.0*iota*gpt+iota**2*gtt)

  ! chi_1 = x R /sqrt(g)
  chi1(:,:) = sqrt(gtt)*r/g

  ! db/dtheta
  allocate(deriv_t(0:nt-1))
  deriv_t(0) = 0.0
  do i=1,nt-1
     deriv_t(i) = -0.5*(-1)**i/tan(0.5*t(i+1))
  enddo
  dbdt(:,:) = 0.0
  chi1t(:,:) = 0.0
  do j=1,np
     do i=1,nt
        do ip=1,nt
           k = ip-i
           if (k < 0) then
              k = k + nt
           endif
           dbdt(i,j)  = dbdt(i,j)+deriv_t(k)*bmag(ip,j)
           chi1t(i,j) = chi1t(i,j)+deriv_t(k)*chi1(ip,j)
        enddo
     enddo
  enddo

  ! db/dphi and dg/dphi
  allocate(deriv_p(0:np-1))
  deriv_p(0) = 0.0
  do i=1,np-1
     deriv_p(i) = -0.5*(-1)**i/tan(0.5*p(i+1))
  enddo
  dbdp(:,:) = 0.0
  dgdp(:,:) = 0.0
  chi1p(:,:) = 0.0
  drdpt(:,:) = 0.0
  dzdpt(:,:) = 0.0
  do j=1,nt
     do i=1,np
        do ip=1,np
           k = ip-i
           if (k < 0) then
              k = k + np
           endif
           dbdp(j,i)  = dbdp(j,i) +deriv_p(k)*bmag(j,ip)
           dgdp(j,i)  = dgdp(j,i) +deriv_p(k)*g(j,ip)
           chi1p(j,i) = chi1p(j,i)+deriv_p(k)*chi1(j,ip)
           drdpt(j,i) = drdpt(j,i)+deriv_p(k)*drdt(j,ip)
           dzdpt(j,i) = dzdpt(j,i)+deriv_p(k)*dzdt(j,ip)
        enddo
     enddo
  enddo

  matsize = 4*nts*nps+2*(nts+nps)+1

  allocate(m_indx(matsize))
  allocate(n_indx(matsize))
  allocate(itype(matsize))
  allocate(bk(nt,np))
  allocate(bkp(nt,np))
  allocate(bk_t(nt,np))
  allocate(bkp_t(nt,np))
  allocate(bk_p(nt,np))
  allocate(bkp_p(nt,np))

  i=0
  do ips=0,nps
     do its=0,nts
        if(its > 0) then 
           ! amn
           i=i+1
           itype(i) = 1
           m_indx(i)     = its
           n_indx(i)     = ips
        endif
        if(ips > 0 .and. its > 0) then 
           !bmn
           i=i+1
           itype(i) = 2
           m_indx(i)     = its
           n_indx(i)     = ips
        endif
        ! cmn
        i=i+1
        itype(i) = 3
        m_indx(i)     = its
        n_indx(i)     = ips
        if(its == 0 .and. ips == 0) then
           indx_c00 = i
        endif
        if(ips  > 0) then
           ! dmn
           i=i+1
           itype(i) = 4
           m_indx(i)     = its
           n_indx(i)     = ips
        endif
     enddo
  enddo

10 format('(',i2,'):',2x,2(1pe14.7,1x))
20 format('(',i2,',',i2,'):',2x,4(1pe14.7,1x))
30 format(t15,4(a,9x))

end subroutine le3_geometry
