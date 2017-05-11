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
  integer :: jr,jt,id,ccw_fac
  integer :: i_field
  integer :: l,ll
  complex :: thfac,carg
  real, dimension(nc,n_species) :: res_loc
  real, dimension(:,:), allocatable :: jloc_c
  real, external :: spectraldiss

  !-------------------------------------------------------------------------
  ! Distributed Bessel-function Gyroaverages

  allocate(jloc_c(2,nc))

  iv_loc = 0
  do iv=nv1,nv2

     iv_loc = iv_loc+1

     is = is_v(iv)
     ix = ix_v(iv)
     ie = ie_v(iv)

     do ic=1,nc

        it = it_c(ic)
        ir = ir_c(ic)

        arg = k_perp(ic)*rho*vth(is)*mass(is)/(z(is)*bmag(it)) &
             *sqrt(2.0*energy(ie))*sqrt(1.0-xi(ix)**2)

        ! Need this for (Phi, A_parallel) terms in GK and field equations

        jloc_c(1,ic) = bessel_j0(abs(arg))

        ! Needed for B_parallel in GK and field equations

        jloc_c(2,ic) = 0.5*(jloc_c(1,ic) + bessel_jn(2,abs(arg)))/bmag(it)

     enddo

     ! psi factors
     efac = 1.0
     jvec_c(1,:,iv_loc) = efac*jloc_c(1,:)
     
     if (n_field > 1) then
        efac = -xi(ix)*sqrt(2.0*energy(ie))*vth(is)
        jvec_c(2,:,iv_loc) = efac*jloc_c(1,:)
        
        if (n_field > 2) then
           efac = 2.0*energy(ie)*(1-xi(ix)**2)*temp(is)/z(is)
           jvec_c(3,:,iv_loc) = efac*jloc_c(2,:)
        endif

     endif
     
     ! chi factors
     do ic=1,nc
        it = it_c(ic)
        fac = rho * temp(is)/(z(is) * bmag(it)) * bpol(it)/bmag(it) &
             * 2.0 * energy(ie)*(1-xi(ix)**2) * k_x(ic)
        
        jxvec_c(1,ic,iv_loc) =  fac * (bmag(it) * jloc_c(2,ic))
        
        if (n_field > 1) then
           efac = -xi(ix)*sqrt(2.0*energy(ie))*vth(is)
           jxvec_c(2,ic,iv_loc) = efac * fac * (bmag(it) * jloc_c(2,ic))
           
           if (n_field > 2) then
              if(n==0) then
                 jxvec_c(3,ic,iv_loc) = 0.0
              else
                 jxvec_c(3,ic,iv_loc) = fac * z(is)*bmag(it)/mass(is) &
                      /(k_perp(ic)*rho)**2 &
                      * (bmag(it) * jloc_c(2,ic) - jloc_c(1,ic))
              endif
           endif
           
        endif
     enddo
  enddo
  
  deallocate(jloc_c)
  do i_field=1,n_field
     call parallel_lib_rtrans_real(jvec_c(i_field,:,:),jvec_v(i_field,:,:))
  enddo
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  ! Conservative upwind factor
  !
  allocate(res_norm(nc,n_species))

  res_loc(:,:) = 0.0

!$omp parallel private(ic,iv_loc,is,ix,ie)
!$omp do reduction(+:res_loc)
  do iv=nv1,nv2
     iv_loc = iv-nv1+1
     is = is_v(iv)
     ix = ix_v(iv)
     ie = ie_v(iv)
     do ic=1,nc
        res_loc(ic,is) = res_loc(ic,is)+w_xi(ix)*w_e(ie)*jvec_c(1,ic,iv_loc)**2 
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

!$omp parallel do private(iv_loc,is,ix,ie,ic)
  do iv=nv1,nv2
     iv_loc = iv-nv1+1
     is = is_v(iv)
     ix = ix_v(iv)
     ie = ie_v(iv)
     do ic=1,nc
        upfac1(ic,iv_loc) = w_e(ie)*w_xi(ix)*abs(xi(ix))*vel(ie)*jvec_c(1,ic,iv_loc)
        upfac2(ic,iv_loc) = jvec_c(1,ic,iv_loc)/res_norm(ic,is)
     enddo
  enddo

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

     vfac(iv_loc) = w_xi(ix)*w_e(ie)*z(is)**2/temp(is)*dens(is)

  enddo

  allocate(sum_den_h(n_theta))
  sum_den_h(:) = 0.0
  do is=1,n_species
     do ie=1,n_energy
        do ix=1,n_xi
           do it=1,n_theta
              sum_den_h(it) = sum_den_h(it) + w_xi(ix)*w_e(ie) &
                   *z(is)**2/temp(is)*dens(is)*dens_rot(it,is)
           enddo
        enddo
     enddo
  enddo

  if (ae_flag == 1) then
     sum_den_h(:) = sum_den_h(:) + dens_ele*dens_ele_rot(:)/temp_ele
  endif

  allocate(sum_den_x(nc))
  if (n_field > 1) allocate(sum_cur_x(nc))

  call cgyro_field_coefficients
  !------------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  ! Zonal flow with adiabatic electrons:
  !
  if (n == 0 .and. ae_flag == 1) then

     allocate(hzf(n_radial,n_theta,n_theta))
     hzf(:,:,:) = 0.0      
     do ir=1,n_radial
        do it=1,n_theta
           hzf(ir,it,it) = k_perp(ic_c(ir,it))**2 * lambda_debye**2 &
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
           xzf(ir,it,it) = k_perp(ic_c(ir,it))**2*lambda_debye**2 &
                * dens_ele/temp_ele+sum_den_x(ic_c(ir,it))
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

     ! Need to allocate these for future use
     allocate(pvec_outr(n_theta))
     allocate(pvec_outi(n_theta))

  endif
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
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  ! Streaming coefficient arrays
  !
  if (ipccw*btccw < 0) then
     ccw_fac = -1
  else
     ccw_fac = 1
  endif

  do ir=1,n_radial
     do it=1,n_theta
        do id=-nup_theta,nup_theta
           jt = modulo(it+id-1,n_theta)+1
           if (it+id < 1) then
              thfac = exp(2*pi*i_c*k_theta*rmin)
              jr = modulo(ir-n*box_size*ccw_fac-1,n_radial)+1
           else if (it+id > n_theta) then
              thfac = exp(-2*pi*i_c*k_theta*rmin)
              jr = modulo(ir+n*box_size*ccw_fac-1,n_radial)+1
           else
              thfac = (1.0,0.0)
              jr = ir
           endif
           dtheta(ic_c(ir,it),id)    = cderiv(id)*thfac
           dtheta_up(ic_c(ir,it),id) = uderiv(id)*thfac*up_theta
           icd_c(ic_c(ir,it),id)     = ic_c(jr,modulo(it+id-1,n_theta)+1)
        enddo
     enddo
  enddo
!$acc enter data copyin(dtheta,dtheta_up,icd_c)

  ! Streaming coefficients (for speed optimization)

!$omp parallel do collapse(2) &
!$omp& private(iv,ic,iv_loc,is,ix,ie,ir,it,carg,u)
  do iv=nv1,nv2
     do ic=1,nc

        iv_loc = iv-nv1+1
        is = is_v(iv)
        ix = ix_v(iv)
        ie = ie_v(iv)
        ir = ir_c(ic) 
        it = it_c(ic)

        u = (pi/n_toroidal)*n

        ! omega_dalpha
        omega_cap_h(ic,iv_loc) = &
             -omega_adrift(it,is)*energy(ie)*(1.0+xi(ix)**2)*&
             (n_toroidal*q/pi/rmin)*(i_c*u)

        ! omega_dalpha [UPWIND: iu -> spectraldiss]
        omega_h(ic,iv_loc) = &
             -abs(omega_adrift(it,is))*energy(ie)*(1.0+xi(ix)**2)*&
             (n_toroidal*q/pi/rmin)*spectraldiss(u,nup_alpha)*up_alpha

        ! (i ktheta) components from drifts
        
        omega_cap_h(ic,iv_loc) = omega_cap_h(ic,iv_loc) &
             - i_c * k_theta * (omega_aprdrift(it,is)*energy(ie)*xi(ix)**2 &
             + omega_cdrift(it,is)*vel(ie)*xi(ix) + omega_rot_drift(it,is) &
             + omega_rot_edrift(it))
        
        u = (2.0*pi/n_radial)*px(ir)

        ! (d/dr) components from drifts
        
        omega_cap_h(ic,iv_loc) = omega_cap_h(ic,iv_loc) & 
             - (n_radial/length)*(i_c*u) &
             * (omega_rdrift(it,is)*energy(ie)*(1.0+xi(ix)**2) &
             + omega_cdrift_r(it,is)*vel(ie)*xi(ix) &
             + omega_rot_drift_r(it,is) &
             + omega_rot_edrift_r(it))
        
        ! (d/dr) upwind components from drifts [UPWIND: iu -> spectraldiss]
        omega_h(ic,iv_loc) = omega_h(ic,iv_loc) &
             - (n_radial/length)*spectraldiss(u,nup_radial)*up_radial &
             * (abs(omega_rdrift(it,is))*energy(ie)*(1.0+xi(ix)**2) &
             + abs(omega_cdrift_r(it,is)*xi(ix))*vel(ie) &
             + abs(omega_rot_drift_r(it,is)) &
             + abs(omega_rot_edrift_r(it)))          
             
        ! omega_star 
        carg = -i_c*k_theta*rho*(dlnndr(is)+dlntdr(is)*(energy(ie)-1.5)) &
             -i_c*k_theta*rho*(sqrt(2.0*energy(ie))*xi(ix)/vth(is) &
             *omega_gammap(it)) -i_c*k_theta*rho*omega_rot_star(it,is)

        omega_s(:,ic,iv_loc) = carg*jvec_c(:,ic,iv_loc)

        ! Profile curvature via wavenumber advection
        carg = k_theta*rho*(sdlnndr(is)+sdlntdr(is)*(energy(ie)-1.5))*length/(2*pi)
        
        omega_ss(:,ic,iv_loc) = carg*jvec_c(:,ic,iv_loc)

     enddo
  enddo
!$acc enter data copyin(omega_cap_h,omega_h,omega_s)
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

     ! 2nd-derivative filter
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

     ! 4th-derivative filter 
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

     ! 6th-derivative filter 
     f(-3) =  -1.0/(60.0*dx)
     f(-2) =   6.0/(60.0*dx)
     f(-1) = -15.0/(60.0*dx)
     f(0)  =  20.0/(60.0*dx)
     f(1)  = -15.0/(60.0*dx)
     f(2)  =   6.0/(60.0*dx)
     f(3)  =  -1.0/(60.0*dx)

  end select

end subroutine advect_schemes
