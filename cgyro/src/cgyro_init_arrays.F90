subroutine cgyro_init_arrays

  use mpi
  use cgyro_globals
  use parallel_lib

  implicit none

  real :: arg
  real :: efac
  real :: u
  real :: fac
  integer :: ir,it,is,ie,ix
  integer :: itm,itl,itor,mytor,itf
  integer :: it_loc
  integer :: jr,jt,id
  integer :: i_field
  integer :: l,ll
  integer :: iltheta_min
  complex :: thfac,carg
  complex :: jval
  real, dimension(:,:,:), allocatable :: res_loc
  real, dimension(:,:,:), allocatable :: jloc_c
  real, dimension(:,:,:), allocatable :: res_norm
  real, dimension(:,:), allocatable :: res_weight
  real, external :: spectraldiss

  !-------------------------------------------------------------------------
  ! Distributed Bessel-function Gyroaverages

  allocate(jloc_c(2,nc,nt1:nt2))

  do itor=nt1,nt2
   iv_loc = 0
   do iv=nv1,nv2

     iv_loc = iv_loc+1

     is = is_v(iv)
     ix = ix_v(iv)
     ie = ie_v(iv)

     do ic=1,nc

        it = it_c(ic)
        ir = ir_c(ic)

        arg = k_perp(ic,itor)*rho*vth(is)*mass(is)/(z(is)*bmag(it)) &
             *vel2(ie)*sqrt(1.0-xi(ix)**2)

        ! Need this for (Phi, A_parallel) terms in GK and field equations

        jloc_c(1,ic,itor) = bessel_j0(arg)

        ! Needed for B_parallel in GK and field equations

        jloc_c(2,ic,itor) = 0.5*(jloc_c(1,ic,itor) + bessel_jn(2,arg))/bmag(it)
        
     enddo

     ! Psi factors: 

     ! J0 phi
     efac = 1.0
     jvec_c(1,:,iv_loc,itor) = efac*jloc_c(1,:,itor)
     
     if (n_field > 1) then
        ! J0 vpar Apar
        efac = -xi(ix)*vel2(ie)*vth(is)
        jvec_c(2,:,iv_loc,itor) = efac*jloc_c(1,:,itor)
        
        if (n_field > 2) then
           ! J2 bpar
           efac = 2.0*energy(ie)*(1-xi(ix)**2)*temp(is)/z(is)
           jvec_c(3,:,iv_loc,itor) = efac*jloc_c(2,:,itor)
        endif

     endif
     
     ! Chi factors (for momentum flux, not GK equation) 
     do ic=1,nc
        it = it_c(ic)
        fac = rho * temp(is)/(z(is) * bmag(it)) * bpol(it)/bmag(it) &
             * 2.0 * energy(ie)*(1-xi(ix)**2) * k_x(ic,itor)
        
        jxvec_c(1,ic,iv_loc,itor) =  fac * (bmag(it) * jloc_c(2,ic,itor))
        
        if (n_field > 1) then
           efac = -xi(ix)*vel2(ie)*vth(is)
           jxvec_c(2,ic,iv_loc,itor) = efac * fac * (bmag(it) * jloc_c(2,ic,itor))
           
           if (n_field > 2) then
              if(itor == 0) then
                 jxvec_c(3,ic,iv_loc,itor) = 0.0
              else
                 jxvec_c(3,ic,iv_loc,itor) = fac * z(is)*bmag(it)/mass(is) &
                      /(k_perp(ic,itor)*rho)**2 &
                      * (bmag(it) * jloc_c(2,ic,itor) - jloc_c(1,ic,itor))
              endif
           endif
           
        endif
     enddo
   enddo
  enddo
 
  deallocate(jloc_c)
#if defined(OMPGPU)
!$omp target enter data map(to:jvec_c)
#elif defined(_OPENACC)
!$acc enter data copyin(jvec_c)
#endif

  do i_field=1,n_field
     call parallel_lib_rtrans_real(jvec_c(i_field,:,:,:),jvec_v(i_field,:,:,:))
  enddo

  if (nonlinear_flag == 1) then
!
! This is a NL transpose
! Reminder: nc ~= n_radial*n_theta
! First half of the transpose is done locally with sub-sampling
!  from (n_field,n_theta,n_radial,nv_loc,nt_loc) -> (n_field,n_radial,n_jtheta,nv_loc,nt_loc,n_toroidal_procs)
! Then AlltoAll finishes the transpose, keeping
!  (n_field,n_radial,n_jtheta,nv_loc,nt_loc,n_toroidal_procs)
! 

#if defined(OMPGPU)
!$omp target teams distribute parallel do simd collapse(6) &
!$omp&   private(itor,it,iltheta_min,mytor,ir,itf,jval)
#elif defined(_OPENACC)
!$acc parallel loop gang vector independent collapse(6) &
!$acc&         private(itor,it,iltheta_min,mytor,ir,itf,jval) &
!$acc&         present(jvec_c_nl,jvec_c,ic_c) &
!$acc&         present(n_toroidal_procs,nt_loc,nv_loc,n_jtheta,n_radial) &
!$acc&         present(nt1,n_theta,n_field,nsplit) default(none)
#else
!$omp parallel do collapse(4) &
!$omp&         private(it_loc,itor,mytor,it,ir,itf,iltheta_min,jval)
#endif
   do itm=1,n_toroidal_procs
    do itl=1,nt_loc
     do iv_loc=1,nv_loc
      do it_loc=1,n_jtheta
       do ir=1,n_radial
        do itf=1,n_field
          iltheta_min = 1+((itm-1)*nsplit)/nv_loc
          it = it_loc+iltheta_min-1
          itor = itl+(itm-1)*nt_loc
          jval = (0.0,0.0)
          if (it <= n_theta) then
            mytor = nt1+itl-1
            ! ic_c(ir,it) = (ir-1)*n_theta+it
            jval = jvec_c(itf,(ir-1)*n_theta+it,iv_loc,mytor)
          endif
          ! else just padding
          jvec_c_nl(itf,ir,it_loc,iv_loc,itor) = jval
        enddo
       enddo
      enddo
     enddo
    enddo
   enddo
   call parallel_slib_distribute_real(n_field*n_radial*n_jtheta*nv_loc*nt_loc,jvec_c_nl)
  endif

  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  ! Conservative upwind factor
  !
  allocate(res_loc(nc,n_species,nt1:nt2))
  allocate(res_norm(nc,n_species,nt1:nt2))
  allocate(res_weight(n_xi,n_energy))

  if (res_weight_power < 1e-2) then
     res_weight(:,:) = 1.0
  else
     do ix=1,n_xi
        do ie=1,n_energy
           res_weight(ix,ie) = (abs(xi(ix))*vel(ie))**res_weight_power
        enddo
     enddo
  endif

  res_loc(:,:,:) = 0.0

!$omp parallel private(ic,iv_loc,is,ix,ie)
!$omp do collapse(2) reduction(+:res_loc)
  do itor=nt1,nt2
   do iv=nv1,nv2
     iv_loc = iv-nv1+1
     is = is_v(iv)
     ix = ix_v(iv)
     ie = ie_v(iv)
     do ic=1,nc
        res_loc(ic,is,itor) = res_loc(ic,is,itor) + &
                w_exi(ie,ix)*jvec_c(1,ic,iv_loc,itor)**2*res_weight(ix,ie)
     enddo
   enddo
  enddo
!$omp end do
!$omp end parallel

  call MPI_ALLREDUCE(res_loc,&
       res_norm,&
       size(res_norm),&
       MPI_DOUBLE_PRECISION,&
       MPI_SUM,&
       NEW_COMM_1,&
       i_err)

!$omp parallel do collapse(2) private(iv_loc,is,ix,ie,ic)
  do itor=nt1,nt2
   do iv=nv1,nv2
     iv_loc = iv-nv1+1
     is = is_v(iv)
     ix = ix_v(iv)
     ie = ie_v(iv)
     do ic=1,nc
        upfac1(ic,iv_loc,itor) = w_exi(ie,ix)*abs(xi(ix))*vel(ie) * &
                jvec_c(1,ic,iv_loc,itor)
        upfac2(ic,iv_loc,itor) = jvec_c(1,ic,iv_loc,itor)*res_weight(ix,ie)/ &
                res_norm(ic,is,itor)
     enddo
   enddo
  enddo

  deallocate(res_weight)
  deallocate(res_norm)
  deallocate(res_loc)
  
#if defined(OMPGPU)
!$omp target enter data map(to:upfac1,upfac2)
#elif defined(_OPENACC)
!$acc enter data copyin(upfac1,upfac2)
#endif

  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  ! Coefficient setup
  !
  allocate(vfac(nv_loc))
  do iv=nv1,nv2

     iv_loc = iv-nv1+1
     is = is_v(iv)
     ix = ix_v(iv)
     ie = ie_v(iv)

     vfac(iv_loc) = w_exi(ie,ix)*z(is)**2/temp(is)*dens(is)

  enddo

  allocate(sum_den_h(n_theta))
  sum_den_h(:) = 0.0
  do is=1,n_species
     do ie=1,n_energy
        do ix=1,n_xi
           do it=1,n_theta
              sum_den_h(it) = sum_den_h(it) + w_exi(ie,ix) &
                   *z(is)**2/temp(is)*dens2_rot(it,is)
           enddo
        enddo
     enddo
  enddo

  if (ae_flag == 1) then
     sum_den_h(:) = sum_den_h(:) + dens_ele*dens_ele_rot(:)/temp_ele
  endif

  allocate(sum_den_x(nc,nt1:nt2))
  if (n_field > 1) allocate(sum_cur_x(nc,nt1:nt2))

  call cgyro_field_coefficients
  !------------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  ! Zonal flow with adiabatic electrons:
  !
  do itor=nt1,nt2
   if (itor == 0 .and. ae_flag == 1) then
     ! since this applies only to itor == 0, we do not need to extend the matrix
     allocate(hzf(n_radial,n_theta,n_theta))
     hzf(:,:,:) = 0.0      
     do ir=1,n_radial
        do it=1,n_theta
           ! my_toroidal==0
           hzf(ir,it,it) = k_perp(ic_c(ir,it),0)**2 * lambda_debye**2 &
                * dens_ele/temp_ele + sum_den_h(it)
           do jt=1,n_theta
              hzf(ir,it,jt) = hzf(ir,it,jt) &
                   - dens_ele*dens_ele_rot(it)/temp_ele*w_theta(jt)
           enddo
        enddo
     enddo

     allocate(work(n_theta))
     allocate(i_piv(n_theta))
     do ir=1,n_radial
        call DGETRF(n_theta,n_theta,hzf(ir,:,:),n_theta,i_piv,info)
        call DGETRI(n_theta,hzf(ir,:,:),n_theta,i_piv,work,n_theta,info)
     enddo
     deallocate(i_piv)
     deallocate(work)

     allocate(xzf(n_radial,n_theta,n_theta))
     xzf(:,:,:) = 0.0     
     do ir=1,n_radial
        do it=1,n_theta
           ! my_toroidal==0
           xzf(ir,it,it) = k_perp(ic_c(ir,it),0)**2*lambda_debye**2 &
                * dens_ele/temp_ele+sum_den_x(ic_c(ir,it),0)
           do jt=1,n_theta
              xzf(ir,it,jt) = xzf(ir,it,jt) &
                   - dens_ele*dens_ele_rot(it)/temp_ele*w_theta(jt)
           enddo
        enddo
     enddo

     allocate(work(n_theta))
     allocate(i_piv(n_theta))
     do ir=1,n_radial
        call DGETRF(n_theta,n_theta,xzf(ir,:,:),n_theta,i_piv,info)
        call DGETRI(n_theta,xzf(ir,:,:),n_theta,i_piv,work,n_theta,info)
     enddo
     deallocate(i_piv)
     deallocate(work)

   endif
  enddo

  if (n_field > 1) deallocate(sum_cur_x)
  deallocate(sum_den_x)

  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  ! Parallel derivative and dissipation stencils
  !
  allocate(cderiv(-nup_theta:nup_theta))
  allocate(uderiv(-nup_theta:nup_theta))
  call advect_schemes(d_theta,nup_theta,cderiv,uderiv)
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  ! Wavenumber advection stencil (coefficients of triangle wave)
  !
  allocate(c_wave(n_wave))
  do l=1,n_wave
     ll = 2*l-1
     c_wave(l) = 2.0/pi/ll**2*(-1)**(l-1)
  enddo
  source = 0.0
  sa     = 0.0
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  ! Streaming coefficient arrays
  !
  do itor=nt1,nt2
   do ir=1,n_radial
     do it=1,n_theta
        do id=-nup_theta,nup_theta
           jt = modulo(it+id-1,n_theta)+1
           if (it+id < 1) then
              thfac = exp(2*pi*i_c*k_theta_base*itor*rmin)
              jr = modulo(ir-itor*box_size*sign_qs-1,n_radial)+1
           else if (it+id > n_theta) then
              thfac = exp(-2*pi*i_c*k_theta_base*itor*rmin)
              jr = modulo(ir+itor*box_size*sign_qs-1,n_radial)+1
           else
              thfac = (1.0,0.0)
              jr = ir
           endif
           dtheta(id, ic_c(ir,it), itor)    = cderiv(id)*thfac
           dtheta_up(id, ic_c(ir,it), itor) = uderiv(id)*thfac*up_theta
           icd_c(id, ic_c(ir,it), itor)     = ic_c(jr,modulo(it+id-1,n_theta)+1)
        enddo
     enddo
   enddo
  enddo
#if defined(OMPGPU)
!$omp target enter data map(to:dtheta,dtheta_up,icd_c,c_wave)
#elif defined(_OPENACC)
!$acc enter data copyin(dtheta,dtheta_up,icd_c,c_wave)
#endif

  ! Streaming coefficients (for speed optimization)

!$omp parallel do collapse(3) &
!$omp& private(iv,ic,iv_loc,is,ix,ie,ir,it,carg,u)
  do itor=nt1,nt2
   do iv=nv1,nv2
     do ic=1,nc

        iv_loc = iv-nv1+1
        is = is_v(iv)
        ix = ix_v(iv)
        ie = ie_v(iv)
        ir = ir_c(ic) 
        it = it_c(ic)

        u = (pi/n_toroidal)*itor

        ! omega_dalpha
        omega_cap_h(ic,iv_loc,itor) = &
             -omega_adrift(it,is)*energy(ie)*(1.0+xi(ix)**2)*&
             (n_toroidal*q/pi/rmin)*(i_c*u)

        ! omega_dalpha [UPWIND: iu -> spectraldiss]
        omega_h(ic,iv_loc,itor) = &
             -abs(omega_adrift(it,is))*energy(ie)*(1.0+xi(ix)**2)*&
             (n_toroidal*q/pi/rmin)*spectraldiss(u,nup_alpha)*up_alpha

        ! (i ktheta) components from drifts        
        omega_cap_h(ic,iv_loc,itor) = omega_cap_h(ic,iv_loc,itor) &
             - i_c*k_theta_base*itor*(omega_aprdrift(it,is)*energy(ie)*xi(ix)**2 &
             + omega_cdrift(it,is)*vel(ie)*xi(ix) + omega_rot_drift(it,is) &
             + omega_rot_edrift(it))
        
        ! Note that we shift the dissipation with px0 (ballooning angle linear mode)
        u = (2.0*pi/n_radial)*(px(ir)+px0)

        ! (d/dr) components from drifts
        
        omega_cap_h(ic,iv_loc,itor) = omega_cap_h(ic,iv_loc,itor) & 
             - (n_radial/length)*i_c*u &
             * (omega_rdrift(it,is)*energy(ie)*(1.0+xi(ix)**2) &
             + omega_cdrift_r(it,is)*vel(ie)*xi(ix) &
             + omega_rot_drift_r(it,is) &
             + omega_rot_edrift_r(it))
        
        ! (d/dr) upwind components from drifts [UPWIND: iu -> spectraldiss]
        omega_h(ic,iv_loc,itor) = omega_h(ic,iv_loc,itor) &
             - (n_radial/length)*spectraldiss(u,nup_radial)*up_radial &
             * (abs(omega_rdrift(it,is))*energy(ie)*(1.0+xi(ix)**2) &
             + abs(omega_cdrift_r(it,is)*xi(ix))*vel(ie) &
             + abs(omega_rot_drift_r(it,is)) &
             + abs(omega_rot_edrift_r(it)))          
             
        ! omega_star 
        carg = -i_c*k_theta_base*itor*rho*(dlnndr(is)+dlntdr(is)*(energy(ie)-1.5)) &
             -i_c*k_theta_base*itor*rho*(vel2(ie)*xi(ix)/vth(is) &
             *omega_gammap(it)) -i_c*k_theta_base*itor*rho*omega_rot_star(it,is)

        omega_s(:,ic,iv_loc,itor) = carg*jvec_c(:,ic,iv_loc,itor)

        ! Profile curvature via wavenumber advection (ix -> d/dp)
        ! See whiteboard notes.
        ! JC: Re-checked sign and normalization (Oct 2019)
        carg = -k_theta_base*itor*length*(sdlnndr(is)+sdlntdr(is)*(energy(ie)-1.5))/(2*pi)

        omega_ss(:,ic,iv_loc,itor) = carg*jvec_c(:,ic,iv_loc,itor)

     enddo
   enddo
  enddo
#if defined(OMPGPU)
!$omp target enter data map(to:omega_cap_h,omega_h,omega_s,omega_ss)
#elif defined(_OPENACC)
!$acc enter data copyin(omega_cap_h,omega_h,omega_s,omega_ss)
#endif
  !-------------------------------------------------------------------------

end subroutine cgyro_init_arrays

! Spectral dissipation function

real function spectraldiss(u,n)

  implicit none
  real, intent(in) :: u
  integer, intent(in) :: n

  select case(n)

  case (1)

     ! 2nd order spectral dissipation
     spectraldiss = 1-cos(u)

  case (2)

     ! 4th order spectral dissipation
     spectraldiss = (3-4*cos(u)+cos(2*u))/6.0

  case (3)

     ! 6th order spectral dissipation
     spectraldiss = (20-30*cos(u)+12*cos(2*u)-2*cos(3*u))/60.0

  case (4)

     ! 8th order spectral dissipation
     spectraldiss = (70-112*cos(u)+56*cos(2*u)-16*cos(3*u)+2*cos(4*u))/280.0

  case default

     print *,'Order out of range in spectraldiss'
     spectraldiss = 0.0
     stop

  end select

end function spectraldiss

subroutine advect_schemes(dx,n,d,f)

  implicit none
  real, intent(in) :: dx
  integer, intent(in) :: n
  real, intent(inout), dimension(-n:n) :: d,f

  select case(n)

  case (1)

     ! 1st-order UPWIND

     ! 2nd-order centered derivative
     d(-1) = -1.0/(2.0*dx)
     d(0)  =  0.0/(2.0*dx)
     d(1)  =  1.0/(2.0*dx)

     ! 2nd-derivative filter (-1/2 * dx * d^2/dx^2)
     f(-1) = -1.0/(2.0*dx)
     f(0)  =  2.0/(2.0*dx)
     f(1)  = -1.0/(2.0*dx)

  case (2)

     ! 3rd-order UPWIND

     ! 4th-order centered derivative
     d(-2) =  1.0/(12.0*dx)
     d(-1) = -8.0/(12.0*dx)
     d(0)  =  0.0/(12.0*dx)
     d(1)  =  8.0/(12.0*dx)
     d(2)  = -1.0/(12.0*dx)

     ! 4th-derivative filter (1/12 * dx^3 * d^4/dx^4)
     f(-2) =  1.0/(12.0*dx)
     f(-1) = -4.0/(12.0*dx)
     f(0)  =  6.0/(12.0*dx)
     f(1)  = -4.0/(12.0*dx)
     f(2)  =  1.0/(12.0*dx)

  case (3)

     ! 5th-order UPWIND

     ! 6th-order centered derivative
     d(-3) =  -1.0/(60.0*dx)
     d(-2) =   9.0/(60.0*dx)
     d(-1) = -45.0/(60.0*dx)
     d(0)  =   0.0/(60.0*dx)
     d(1)  =  45.0/(60.0*dx)
     d(2)  =  -9.0/(60.0*dx)
     d(3)  =   1.0/(60.0*dx)

     ! 6th-derivative filter (-1/60 * dx^5 * d^6/dx^6)
     f(-3) =  -1.0/(60.0*dx)
     f(-2) =   6.0/(60.0*dx)
     f(-1) = -15.0/(60.0*dx)
     f(0)  =  20.0/(60.0*dx)
     f(1)  = -15.0/(60.0*dx)
     f(2)  =   6.0/(60.0*dx)
     f(3)  =  -1.0/(60.0*dx)

  end select

end subroutine advect_schemes
