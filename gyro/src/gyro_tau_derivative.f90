!------------------------------------------------
! gyro_tau_derivative.f90 [caller gyro_rhs_total]
!
! PURPOSE:
!  Compute the tau-derivative terms in the RHS 
!  (i.e., derivatives in the theta directions)
!  using 5-point rules.
!
! NOTES:
!  Upwind dissipation is necessary for Landau
!  damping of zonal flows.  No dissipation will
!  result in "bouncing" instead of damping.
!------------------------------------------------

subroutine gyro_tau_derivative

  use gyro_globals
  use gyro_pointers
  use ompdata

  !---------------------------------------
  implicit none
  !
  integer :: p
  !
  integer :: mff,mf
  integer :: mc
  integer :: mb,mbb
  !
  real, dimension(n_x) :: v0
  real :: upwind
  !
  complex :: pff,pf
  complex :: pc
  complex :: pb,pbb
  !
  complex, dimension(n_x,n_stack) :: fh0
  complex, dimension(n_x,n_stack) :: fh
  !
  complex :: temp1
  complex :: temp2
  !---------------------------------------

  ! We should be calling this for gyrokinetic species only

!$omp parallel private(upwind,p_nek_loc,ie,k,ck,p,mff,mf,mc,mb,mbb,pff,pf,pc,pb,pbb,temp1,temp2) &
!$omp private(is,p_nek,i,m)
  do is=1,n_gk

     ! q has sign ipccw*btccw, but upwind dissipation must have 
     ! fixed sign independent of q:
     upwind = orbit_upwind_vec(is)*abs(q_norm)/q_norm

     p_nek_loc = 0
     do p_nek=1+i_proc_1,n_nek_1,n_proc_1

        p_nek_loc = p_nek_loc+1

        ie = nek_e(p_nek)  
        k  = nek_k(p_nek)   

        ck = class(k)

        do i = ibeg, iend

           fh0(i,:) = h(:,i,p_nek_loc,is) 
           fh(i,:)  = fh0(i,:)+ &
                z(is)*alpha_s(is,i)*gyro_u(:,i,p_nek_loc,is)

           v0(i) = mu(is)*sqrt(tem_s(is,i))*v_theta(i,ie,k,is)/d_tau(ck)

           RHS_dt(:,i,p_nek_loc,is) = 0.0

           do m=1,n_stack

              p = p_phys(ck,m)

              mff = m_cyc(ck,m+2,p)
              mf  = m_cyc(ck,m+1,p)
              mc  = m_cyc(ck,m,p)
              mb  = m_cyc(ck,m-1,p)
              mbb = m_cyc(ck,m-2,p)

              pff = p_cyc(ck,i,m+2,p)
              pf  = p_cyc(ck,i,m+1,p)
              pc  = p_cyc(ck,i,m,p)
              pb  = p_cyc(ck,i,m-1,p)
              pbb = p_cyc(ck,i,m-2,p)

              temp1 = (-1.0/12.0)*pff*fh(i,mff) &
                   +(8.0/12.0)*pf*fh(i,mf) &
                   +(-8.0/12.0)*pb*fh(i,mb) &
                   +(1.0/12.0)*pbb*fh(i,mbb) 

              temp2 = (1.0/12.0)*pff*fh0(i,mff) &
                   +(-4.0/12.0)*pf*fh0(i,mf) &
                   +(6.0/12.0)*pc*fh0(i,mc) &
                   +(-4.0/12.0)*pb*fh0(i,mb) &
                   +(1.0/12.0)*pbb*fh0(i,mbb) 


              !------------------------------------------------
              ! Entropy tracking:
              !
              RHS_dt(m,i,p_nek_loc,is) =  RHS_dt(m,i,p_nek_loc,is) &
                   +real(v0(i)*(-upwind*temp2)*conjg(fh(i,m)))
              !------------------------------------------------

              RHS(m,i,p_nek_loc,is) = RHS(m,i,p_nek_loc,is)+ &
                   v0(i)*(sigma_tau(ck,p)*temp1-upwind*temp2)

           enddo ! m

        enddo ! i

     enddo ! p_nek

  enddo ! is
!$omp end parallel

end subroutine gyro_tau_derivative
