subroutine cgyro_rhs_comm_async_hx
  use cgyro_nl_comm
  use cgyro_globals

  implicit none

  if (nonlinear_flag == 1) then
       ! prepare and transfer h_x (first half)
       call cgyro_nl_fftw_comm1_async
  endif

end subroutine cgyro_rhs_comm_async_hx

subroutine cgyro_rhs_comm_async_fd
  use cgyro_nl_comm
  use cgyro_globals

  implicit none

  if (nonlinear_flag == 1) then
       ! transfer fields
       call cgyro_nl_fftw_comm2_async
       ! and the other half of h_x
       call cgyro_nl_fftw_comm3_async
  endif

end subroutine cgyro_rhs_comm_async_fd

subroutine cgyro_rhs_r_comm_async(ij)
  use cgyro_nl_comm
  use cgyro_globals

  implicit none

  integer, intent(in) :: ij

  if (nonlinear_flag == 1) then
       call cgyro_nl_fftw_comm1_r(ij)
  endif

end subroutine cgyro_rhs_r_comm_async

! Note: Calling test propagates the async operations in some MPI implementations
subroutine cgyro_rhs_comm_test
  use cgyro_nl_comm
  use cgyro_globals

  implicit none

  if (nonlinear_flag == 1) then
     call cgyro_nl_fftw_comm_test
  endif

end subroutine cgyro_rhs_comm_test

subroutine cgyro_rhs_comp1(ij,update_cap)

  use timer_lib
  use cgyro_globals
  use cgyro_field_mod, only : field

  implicit none

  integer, intent(in) :: ij
  logical, intent(in) :: update_cap
  !--------------------------------
  integer :: is,itor
  complex :: h_el,cap_el,my_psi,rhs_el

  call timer_lib_in('str')

#if (!defined(OMPGPU)) && defined(_OPENACC)
!$acc data present(h_x,rhs,field,cap_h_c)

!$acc data  &
!$acc& present(rhs) &
!$acc& present(cap_h_c,z,temp,jvec_c) &
!$acc& present(is_v,ix_v,ie_v,it_c) &
!$acc& present(omega_cap_h,omega_h,omega_s)

#endif

  ! get the initial rhs initialization
  ! that depends on h_x,cap_h_c and field only
#if defined(OMPGPU)
  ! no async for OMPGPU for now
!$omp target teams distribute parallel do simd collapse(3) &
!$omp&  firstprivate(nv2,nv1,nt2,nt1,nc,update_cap,ij) &
!$omp&  private(iv,ic,iv_loc,rhs_el,h_el,is,cap_el,my_psi)
#elif defined(_OPENACC)
!$acc parallel loop gang vector collapse(3) async(1) &
!$acc&  firstprivate(nv2,nv1,nt2,nt1,nc,update_cap,ij) &
!$acc&  private(iv,ic,iv_loc,rhs_el,h_el,is,cap_el,my_psi)
#else
!$omp parallel do collapse(2) &
!$omp&  firstprivate(nv2,nv1,nt2,nt1,nc,update_cap,ij) &
!$omp&  private(iv,iv_loc,itor,is,ic,rhs_el,h_el,cap_el,my_psi) 
#endif
  do itor=nt1,nt2
   do iv=nv1,nv2
     do ic=1,nc
        iv_loc = iv-nv1+1
        h_el = h_x(ic,iv_loc,itor)
        if (update_cap) then
           is = is_v(iv)
           my_psi = sum( jvec_c(:,ic,iv_loc,itor)*field(:,ic,itor))
           cap_el = h_el+my_psi*z(is)/temp(is)
           cap_h_c(ic,iv_loc,itor) = cap_el
        else
           cap_el = cap_h_c(ic,iv_loc,itor)
        endif
        ! Diagonal terms
        rhs_el = &
             omega_cap_h(ic,iv_loc,itor)*cap_el+&
             omega_h(ic,iv_loc,itor)*h_el

        rhs(ic,iv_loc,itor,ij) = rhs_el + &
             sum(omega_s(:,ic,iv_loc,itor)*field(:,ic,itor))
     enddo
   enddo
  enddo

  call cgyro_rhs_comm_test

#if (!defined(OMPGPU)) && defined(_OPENACC)
  !no async for OMPGPU for now
!$acc wait(1)

!$acc end data    

! h_x
!$acc end data 
#endif

  call timer_lib_out('str')

end subroutine cgyro_rhs_comp1

subroutine cgyro_rhs_comp2(ij)

  use timer_lib
  use cgyro_globals
  use cgyro_field_mod, only : field

  implicit none

  integer, intent(in) :: ij
  !--------------------------------
  integer :: itor,ir,it
  ! ir loop specific
  integer :: itorbox
  !integer :: iv_loc
  integer :: is
  integer :: jr0(0:2)   ! n_theta*(pre-compute jr-1)
  real :: my_vel,my_xi
  ! it loop specific
  !integer :: ic
  integer :: id
  integer :: itd   ! precompute modulo(it+id-1,n_theta)+1, use for iteration
  integer :: itd_class
  integer :: jc
  real :: rval,rval2s
  complex :: thfac
  real :: z_temp
  complex :: g_val
  complex :: rhs_stream
  ! CPU/MPI conservative-upwind path (post-stencil projection)
  complex :: rhs_diss

  call timer_lib_in('str')

#if (!defined(OMPGPU)) && defined(_OPENACC)
!$acc data  &
!$acc& present(rhs) &
!$acc& present(cap_h_c,h_x,field,jvec_c,z,temp) &
!$acc& present(is_v,ix_v,ie_v) &
!$acc& present(omega_stream,xi,vel) &
!$acc& present(thfac_itor,cderiv,uderiv) &
!$acc& present(upwind_flux)

#endif
  ! add stream to rhs

  ! =====================================================================
  ! Conservative-upwind path -- post-stencil projection
  !
  !     D g = -c Q S |v_par| g ,      Q = I - P
  !
  ! P is the ORTHOGONAL projection (in the free-energy metric <a,b>=int f0 a b)
  ! at each theta point onto J0(gamma_a) and, for n_field>1, J0*v_par -- i.e.
  ! it subtracts exactly the J0 / J0*v_par component of the flux, no |v_par|
  ! weight (see upfac_num/upfac_cur).  It is applied to the OUTPUT of the
  ! theta-nonlocal stencil S, using the LOCAL J0.  Because P Q = 0 at each point,
  !
  !     (moment functional) . D g = P Q S g = 0   exactly,
  !
  ! for any g, any grid, any d(J0)/d(theta) -- i.e. number and current are
  ! conserved exactly even where J0(gamma_a) varies rapidly with the
  ! ballooning angle.  (Projecting the INPUT instead, D = S Q, leaves a
  ! residual P S Q ~ [S,P] ~ d(J0)/d(theta): the spurious grid-scale source
  ! behind the slow finite-beta instability.  See Desktop/upwind_projection.)
  !
  ! For n_field>1 Q removes BOTH the number moment (jvec_c(1)=J0, couples to
  ! phi/Poisson) and the parallel-current moment (jvec_c(2)~J0*v_par, couples
  ! to A_par/Ampere).  The two are parity-orthogonal, so they are independent
  ! subtractions.  Restoring current conservation removes the spurious A_par
  ! drive that the 2023 (number-only) scheme had dropped.
  !
  ! g_x logically holds |v_par|*g (unprojected; set in cgyro_upwind).  Here:
  ! Pass A : central streaming (advection) -> rhs, and S[g_x] -> upwind_flux
  ! Reduce : number & current moments of upwind_flux over velocity  (MPI)
  ! Pass B : subtract the local-J0(ic) projection(s); add dissipation to rhs
  ! =====================================================================

  ! ---- Pass A: advection -> rhs; raw dissipation stencil -> upwind_flux
#if defined(OMPGPU)
  ! no async for OMPGPU for now
!$omp target teams distribute parallel do simd collapse(4) &
!$omp&  firstprivate(n_radial,nv2,nv1,nt2,nt1,n_theta) &
!$omp&  firstprivate(sign_qs,nup_theta,ij,box_size,n_field) &
!$omp&  private(itor,iv,ir,it) &
!$omp&  private(itorbox,iv_loc,is,jr0,my_vel,my_xi,g_val,z_temp) &
!$omp&  private(ic,id,itd,itd_class,jc,rval,rval2s,thfac,rhs_stream,rhs_diss)
#elif defined(_OPENACC)
!$acc  parallel loop gang vector collapse(4) async(1) &
!$acc&  firstprivate(n_radial,nv2,nv1,nt2,nt1,n_theta) &
!$acc&  firstprivate(sign_qs,nup_theta,ij,box_size,n_field) &
!$acc&  private(itor,iv,ir,it) &
!$acc&  private(itorbox,iv_loc,is,jr0,my_vel,my_xi,g_val,z_temp) &
!$acc&  private(ic,id,itd,itd_class,jc,rval,rval2s,thfac,rhs_stream,rhs_diss)
#else
!$omp parallel do collapse(3) &
!$omp&  firstprivate(n_radial,nv2,nv1,nt2,nt1,n_theta) &
!$omp&  firstprivate(sign_qs,nup_theta,ij,box_size,n_field) &
!$omp&  private(itor,iv,ir,it) &
!$omp&  private(itorbox,iv_loc,is,jr0,my_vel,my_xi,g_val,z_temp) &
!$omp&  private(ic,id,itd,itd_class,jc,rval,rval2s,thfac,rhs_stream,rhs_diss)
#endif
  do itor=nt1,nt2
   do iv=nv1,nv2
    do ir=1,n_radial
#if defined(OMPGPU) || defined(_OPENACC)
     ! keep loop high for maximal collapse on GPU
     do it=1,n_theta
#endif
        itorbox = itor*box_size*sign_qs
        iv_loc = iv-nv1+1

        is = is_v(iv)
        my_vel = vel(ie_v(iv))
        my_xi = xi(ix_v(iv))
        jr0(0) = n_theta*modulo(ir-itorbox-1,n_radial)
        jr0(1) = n_theta*(ir-1)
        jr0(2) = n_theta*modulo(ir+itorbox-1,n_radial)

        if (n_field > 1) then
          z_temp = z(is)/temp(is)
        else
          z_temp = 0.0
        endif

#if !(defined(OMPGPU) || defined(_OPENACC))
        ! loop as late as possible on CPU, to minimize recompute
        ! compiler will vectorize
        do it=1,n_theta
#endif
          ic = (ir-1)*n_theta + it ! ic_c(ir,it)

          rval2s = omega_stream(it,is,itor)
          rval  = rval2s*my_vel*my_xi

          rhs_stream = 0.0   ! central advection D[cap_h]
          rhs_diss   = 0.0   ! upwind dissipation S[g_x], g_x = |v_par|*(I-P)g

          ! historical reminder
          !icd_c(ic, id, itor)     = ic_c(jr,modulo(it+id-1,n_theta)+1)
          !jc = icd_c(ic, id, itor)
          !dtheta(ic, id, itor)    := cderiv(id)*thfac
          !dtheta_up(ic, id, itor) := uderiv(id)*thfac*up_theta

          itd = n_theta+it-nup_theta
          itd_class = 0
          jc = jr0(itd_class)+itd
          thfac = thfac_itor(itd_class,itor)
          do id=-nup_theta,nup_theta
              if (itd > n_theta) then
                itd = itd - n_theta
                itd_class = itd_class + 1
                jc = jr0(itd_class)+itd
                thfac = thfac_itor(itd_class,itor)
              endif

              ! g_val = g_x(jc,iv_loc,itor)
              g_val = h_x(jc,iv_loc,itor)
              if (n_field > 1) then
                 g_val = g_val + &
                   z_temp*jvec_c(2,jc,iv_loc,itor)*field(2,jc,itor)
              endif
              g_val = abs(my_xi)*my_vel*g_val

              rhs_stream = rhs_stream &
                - thfac*rval*cderiv(id)*cap_h_c(jc,iv_loc,itor)
              rhs_diss = rhs_diss &
                + thfac*uderiv(id)*g_val
              itd = itd + 1
              jc = jc + 1
          enddo

          rhs(ic,iv_loc,itor,ij) = rhs(ic,iv_loc,itor,ij) + rhs_stream
          upwind_flux(ic,iv_loc,itor) = rhs_diss
     enddo
    enddo
   enddo
  enddo

  call cgyro_rhs_comm_test
#if (!defined(OMPGPU)) && defined(_OPENACC)
!$acc wait(1)

!$acc end data    

#endif

  call timer_lib_out('str')

end subroutine cgyro_rhs_comp2

subroutine cgyro_rhs_upwind_comm(ij)

  use timer_lib
  use cgyro_globals
  use parallel_lib
  use cgyro_math

  implicit none

  integer, intent(in) :: ij
  !--------------------------------
  integer :: itor,is
  integer :: ix,ie
  integer :: n_upwind_arr
  integer :: iu
  complex :: res_loc(2) ! keep it fixed for convenience and better optimization


  call timer_lib_in('str')
  n_upwind_arr = 2
  if (n_field<2) n_upwind_arr = 1 ! do not need current if I have a single field

  ! ---- number (J0) and current (J0*v_par) moments of the flux ----
  ! Both moments are reduced over velocity.  The current moment is only
  ! meaningful (and jvec_c(2) only allocated) for n_field>1.
#if defined(OMPGPU)
!$omp target teams distribute parallel do simd collapse(3) &
!$omp&  firstprivate(nt1,nt2,ns1,ns2,nv1,nv2,nc,n_field,n_upwind_arr) &
!$omp&  private(res_loc,iv,iv_loc,ix,ie,iu)
#elif defined(_OPENACC)
!$acc parallel loop gang vector collapse(3) async(1) &
!$acc&  firstprivate(nt1,nt2,ns1,ns2,nv1,nv2,nc,n_field,n_upwind_arr) &
!$acc&  private(res_loc,iv,iv_loc,ix,ie,iu) &
!$acc&  present(upwind_res_loc) &
!$acc&  present(w_exi,jvec_c,upwind_flux,ix_v,ie_v)
#else
!$omp parallel do collapse(3) &
!$omp&  firstprivate(nt1,nt2,ns1,ns2,nv1,nv2,nc,n_field,n_upwind_arr) &
!$omp&  private(res_loc,iv,iv_loc,ix,ie,iu)
#endif
  do itor=nt1,nt2
   do is=ns1,ns2
     do ic=1,nc
       res_loc(:) = (0.0,0.0)
       do iv=nv1,nv2
          iv_loc = iv-nv1+1
          if (is == is_v(iv)) then
             ix = ix_v(iv)
             ie = ie_v(iv)
             do iu=1,n_upwind_arr
               res_loc(iu) = res_loc(iu) &
                 + w_exi(ie,ix)*jvec_c(iu,ic,iv_loc,itor)*upwind_flux(ic,iv_loc,itor)
             enddo
          endif
       enddo
       do iu=1,n_upwind_arr
         upwind_res_loc(ic,iu,is,itor) = res_loc(iu)
       enddo
     enddo
   enddo
  enddo

  call cgyro_rhs_comm_test
#if (!defined(OMPGPU)) && defined(_OPENACC)
!$acc wait(1)
#endif

  if (upwind_single_flag == 0) then
    call timer_lib_out('str')
    call timer_lib_in('str_comm')
    call cgyro_rhs_comm_test
    call parallel_clib_sum_upwind(upwind_res_loc,upwind_res)
    call cgyro_rhs_comm_test
    call timer_lib_out('str_comm')
  else
    ! We will do the upwind AllReduce using fp32 precision
    ! The local buffers are small enough, so computed them using full precision
    ! and we now make a copy to the final comm buffer.
    ! The expensive part is only comm, not local compute, 
    ! so this logic is acceptable.
    call cgyro_cmpl_copy_to_fp32(nt_loc*ns_loc*nc*n_upwind_arr, &
                                 upwind32_res_loc, upwind_res_loc)
    call timer_lib_out('str')
    call timer_lib_in('str_comm')
    call cgyro_rhs_comm_test
    call parallel_clib_sum_upwind32(upwind32_res_loc,upwind32_res)
    call cgyro_rhs_comm_test
    call timer_lib_out('str_comm')
    call timer_lib_in('str')
    ! And we need to copy back, too
    call cgyro_cmpl_copy_from_fp32(nt_loc*ns_loc*nc*n_upwind_arr, &
                                   upwind_res, upwind32_res)
    call timer_lib_out('str')
  endif
end subroutine cgyro_rhs_upwind_comm

subroutine cgyro_rhs_comp3(ij)

  use timer_lib
  use cgyro_globals

  implicit none

  integer, intent(in) :: ij
  !--------------------------------
  integer :: itor,it,is
  integer :: n_upwind_arr
  integer :: iu
  real :: rval2
  complex :: cflux

  call timer_lib_in('str')
  n_upwind_arr = 2
  if (n_field<2) n_upwind_arr = 1 ! do not need current if I have a single field

  ! ---- Pass B: subtract number (and current) projection; add to rhs
  ! cflux = (I - P_number [- P_current]) applied to the dissipation flux.
  ! The current subtraction (n_field>1) restores current conservation,
  ! removing the spurious parallel-current source that drives A_par.
#if defined(OMPGPU)
!$omp target teams distribute parallel do simd collapse(3) &
!$omp&  firstprivate(nv2,nv1,nt2,nt1,ij,up_theta,nc,n_upwind_arr) &
!$omp&  private(iv_loc,is,it,rval2,cflux)
#elif defined(_OPENACC)
!$acc parallel loop gang vector collapse(3) async(1) &
!$acc&  firstprivate(nv2,nv1,nt2,nt1,ij,up_theta,nc,n_upwind_arr) &
!$acc&  private(iv_loc,is,it,rval2,cflux) &
!$acc&  present(omega_stream,rhs) &
!$acc&  present(upwind_flux,upfac_num,upwind_res)
#else
!$omp parallel do collapse(3) &
!$omp&  firstprivate(nv2,nv1,nt2,nt1,ij,up_theta,nc,n_upwind_arr) &
!$omp&  private(iv_loc,is,it,rval2,cflux)
#endif
  do itor=nt1,nt2
   do iv=nv1,nv2
     do ic=1,nc
        iv_loc = iv-nv1+1
        is = is_v(iv)
        it = it_c(ic)
        rval2 = abs(omega_stream(it,is,itor))
        cflux = upwind_flux(ic,iv_loc,itor)
        do iu=1,n_upwind_arr
          cflux = cflux &
                - upfac_num(ic,iu,iv_loc,itor)*upwind_res(ic,iu,is,itor)
        enddo
        rhs(ic,iv_loc,itor,ij) = rhs(ic,iv_loc,itor,ij) &
             - rval2*up_theta*cflux
     enddo
   enddo
  enddo

  call cgyro_rhs_comm_test
#if (!defined(OMPGPU)) && defined(_OPENACC)
!$acc wait(1)
#endif

  call timer_lib_out('str')

end subroutine cgyro_rhs_comp3

#if defined(OMPGPU) || defined(_OPENACC)
! gpu code

subroutine cgyro_rhs_trap(ij)

  use cgyro_globals

  implicit none

  integer, intent(in) :: ij
  !--------------------------------

  if (explicit_trap_flag == 1) then
     ! we should never get in here... should have failed during init
     ! hard abort if we somehow end up here
     call abort
  endif

end subroutine cgyro_rhs_trap

#else
! cpu code

subroutine cgyro_rhs_trap(ij)

  use timer_lib
  use parallel_lib
  use cgyro_globals

  implicit none

  integer, intent(in) :: ij
  !--------------------------------
  integer :: is, ix, ie, js, jx, je, jv, it, j, k, itor
  complex, dimension(:,:), allocatable :: rhs_trap
  complex, dimension(:), allocatable   :: bvec_trap
  integer :: nj_loc

  ! Explicit trapping term
  if (explicit_trap_flag == 1) then
     ! assert(n_sim==1), not supported else
     call timer_lib_in('str')

     allocate(rhs_trap(nc,nv_loc))
     allocate(bvec_trap(nv))
     call parallel_lib_rtrans_pack(cap_h_c)
     call parallel_lib_r_do(cap_h_v)
     call parallel_lib_nj_loc(nj_loc)
     
     do itor=nt1,nt2
      do ic=nc1,nc2
        ic_loc = ic-nc1+1
        it = it_c(ic)
        
        bvec_trap(:) = (0.0,0.0)
        do iv=1,nv
           is = is_v(iv)
           ix = ix_v(iv)
           ie = ie_v(iv)
           do jv=1,nv
              js = is_v(jv)
              jx = ix_v(jv)
              je = ie_v(jv)
              if (is == js) then
                 if (ie == je) then
                    bvec_trap(iv) = bvec_trap(iv) - (omega_trap(it,is,itor) * vel(ie) &
                         + omega_rot_trap(it,is) / vel(ie)) &
                         * (1.0 - xi(ix)**2) * xi_deriv_mat(ix,jx) * cap_h_v(ic_loc,itor,jv,i_sim)
                 endif
                 if (ix == jx) then
                    bvec_trap(iv) = bvec_trap(iv) - omega_rot_u(it,is) * xi(ix) &
                         * e_deriv1_rot_mat(ie,je)/sqrt(1.0*e_max) * cap_h_v(ic_loc,itor,jv,i_sim)
                 endif
              endif
           enddo
        enddo

        do k=1,nproc
           do j=1,nj_loc
              fsendf(j,itor,ic_loc,k) = bvec_trap(j+(k-1)*nj_loc)
           enddo
        enddo
      enddo
     enddo
     
     call parallel_lib_f_i_do(cap_h_ct)
     do itor=nt1,nt2
      do iv=nv1,nv2
        iv_loc = iv-nv1+1
        do ic=1,nc
           rhs_trap(ic,iv_loc) = cap_h_ct(iv_loc,itor,ic)
        enddo
      enddo

      rhs(:,:,itor,ij) = rhs(:,:,itor,ij) +  rhs_trap(:,:)
     enddo
     
     deallocate(rhs_trap)
     deallocate(bvec_trap)
     
     call timer_lib_out('str')
  endif
  

end subroutine cgyro_rhs_trap

#endif

subroutine cgyro_rhs(ij,update_cap)

  use timer_lib
  use cgyro_globals
  use cgyro_nl

  implicit none

  integer, intent(in) :: ij
  logical, intent(in) :: update_cap
  !--------------------------------

  ! fields is ready by now
  call cgyro_rhs_comm_async_fd

  call cgyro_rhs_comp1(ij,update_cap)

  call cgyro_rhs_comm_test

  ! Wavenumber advection (ExB shear and/or global) terms
  if (source_flag == 1) then
     call cgyro_globalshear(ij)
  endif
     
  call cgyro_rhs_comm_test

  ! Nonlinear evaluation [f,g]
  if (nonlinear_flag == 1) then
     ! assumes someone already started the input comm
     ! and will finish the output comm
     call cgyro_nl_fftw()
  endif

  call cgyro_rhs_comm_test

  call cgyro_rhs_comp2(ij)
  call cgyro_rhs_upwind_comm(ij)
  call cgyro_rhs_comp3(ij)

  call cgyro_rhs_trap(ij)

  ! updates rhs
  call cgyro_rhs_r_comm_async(ij)

end subroutine cgyro_rhs
