!-----------------------------------------------------------------
! neo_energy_grid.f90
!
! PURPOSE:
!  Energy-basis matrices for NEO, in the basis
!
!    phi_n(x) = x^beta L_n^(alpha)(x^2),
!    alpha = e_lag/2, beta = xi_beta_l.
!
!  Every integral is done by quadrature, with L_n from its
!  three-term recurrence
!
!  - Streaming and source matrices: closed form, or Gauss-Laguerre
!    (exact).
!  - Test-particle operator: composite Gauss-Legendre in x.
!  - Field-particle operator (full FP): composite Gauss-Legendre
!    Rosenbluth inner integrals always carry the positive power
!    and are done in scaled form,
!      int_0^x y^p e^(-c y^2) L_n(c y^2) dy
!        = x^(p+1) int_0^1 t^p e^(-u) L_n(u) dt,  u = c x^2 t^2,
!    by one Gauss-Legendre rule in t whose Laguerre table is
!    shared by every l.
!
!  Supported:  collision_model = 1 (Connor), 4 (full FP),
!                                5 (FP test + ad hoc field)
!              laguerre_method = 1 (Laguerre 1/2+3/2), 2 (Sonine)
!-----------------------------------------------------------------

module neo_energy_grid

  implicit none

  public :: ENERGY_basis_ints_alloc, ENERGY_basis_ints, &
       ENERGY_coll_ints_alloc, ENERGY_coll_ints

  ! energy matrices
  ! (energy, energy, xi, left/right diagonal in xi)
  real, dimension(:,:,:,:), allocatable :: emat_e05, emat_en05, emat_e05de, &
       emat_e0, emat_e1

  ! energy vectors
  real, dimension(:,:), allocatable :: evec_e0, evec_e1, evec_e2, evec_e05, evec_e105

  ! xi vectors -- order of associated Laguerre basis
  integer, dimension(:), allocatable  :: e_lag
  integer, dimension(:), allocatable  :: xi_beta_l
  integer, parameter  :: e_alpha=2

  ! collision matrices
  ! (species,species,energy,energy,xi)
  real, dimension(:,:,:,:,:), allocatable :: emat_coll_test, emat_coll_field

  ! private variables
  logical, private :: initialized_basis = .false.
  logical, private :: initialized_coll  = .false.

  ! monomial basis x^(2n+beta) in place of Laguerre (write_fullcoll_mono only)
  logical, private :: mono = .false.

  ! composite Gauss-Legendre: points per panel, panel width
  integer, parameter, private :: ng = 16
  real, parameter, private    :: hpan = 0.25

  ! outer-node block size for the field-particle inner integrals
  integer, parameter, private :: nblk = 256

contains

  !===============================================================
  ! Allocation
  !===============================================================

  subroutine ENERGY_basis_ints_alloc(flag)

    use neo_globals, only : n_energy, n_xi, laguerre_method
    implicit none
    integer, intent (in) :: flag
    integer :: ix

    if (flag == 1) then
       if(initialized_basis) return

       allocate(emat_e05(0:n_energy,0:n_energy,0:n_xi,2))
       allocate(emat_en05(0:n_energy,0:n_energy,0:n_xi,2))
       allocate(emat_e05de(0:n_energy,0:n_energy,0:n_xi,2))
       allocate(emat_e0(0:n_energy,0:n_energy,0:n_xi,1))
       allocate(emat_e1(0:n_energy,0:n_energy,0:n_xi,3))

       allocate(evec_e0(0:n_energy,0:n_xi))
       allocate(evec_e1(0:n_energy,0:n_xi))
       allocate(evec_e2(0:n_energy,0:n_xi))
       allocate(evec_e05(0:n_energy,0:n_xi))
       allocate(evec_e105(0:n_energy,0:n_xi))

       allocate(e_lag(0:n_xi))
       allocate(xi_beta_l(0:n_xi))

       select case (laguerre_method)

       case (1)

          ! Laguerre 1/2+3/2
          e_lag(:) = 3
          e_lag(0) = 1
          xi_beta_l(:) = 1
          xi_beta_l(0) = 0

       case (2)

          ! Sonine
          do ix=0,n_xi
             e_lag(ix) = 2*ix + 1
             xi_beta_l(ix) = ix
          enddo

       case default

          call neo_error('ERROR: (NEO) laguerre_method must be 1 or 2')
          e_lag(:) = 3
          e_lag(0) = 1
          xi_beta_l(:) = 1
          xi_beta_l(0) = 0

       end select

       initialized_basis = .true.

    else
       if(.NOT. initialized_basis) return

       deallocate(emat_e05)
       deallocate(emat_en05)
       deallocate(emat_e05de)
       deallocate(emat_e0)
       deallocate(emat_e1)
       deallocate(evec_e0)
       deallocate(evec_e1)
       deallocate(evec_e2)
       deallocate(evec_e05)
       deallocate(evec_e105)
       deallocate(e_lag)
       deallocate(xi_beta_l)

       initialized_basis = .false.

    endif

  end subroutine ENERGY_basis_ints_alloc

  subroutine ENERGY_coll_ints_alloc(flag)
    use neo_globals, only : n_species, n_energy, n_xi
    implicit none
    integer, intent (in) :: flag

    if(flag == 1) then
       if(initialized_coll) return

       allocate(emat_coll_test(n_species,n_species,0:n_energy,0:n_energy,0:n_xi))
       allocate(emat_coll_field(n_species,n_species,0:n_energy,0:n_energy,0:n_xi))

       initialized_coll = .true.

    else
       if(.NOT. initialized_coll) return

       deallocate(emat_coll_test)
       deallocate(emat_coll_field)

       initialized_coll = .false.

    endif

  end subroutine ENERGY_coll_ints_alloc

  !===============================================================
  ! Streaming and source integrals
  !
  !   int_0^inf x^q e^{-x^2} phi_i phi_j dx
  !     = (1/2) int_0^inf s^gam e^{-s} L_i L_j ds,
  !   gam = (q+beta_i+beta_j-1)/2.
  !
  ! Vectors and emat_e0 in closed form; the other matrices by
  ! Gauss-Laguerre with n_energy+2 points (exact).
  !===============================================================

  subroutine ENERGY_basis_ints
    use neo_globals, only : n_energy, n_xi
    implicit none
    integer :: ie, ix, jx, kx
    real :: alpha, beta

    ! vectors: int_0^inf x^q e^{-x^2} phi_i dx
    do ix=0,n_xi
       alpha = 0.5*e_lag(ix)
       beta  = 1.0*xi_beta_l(ix)
       do ie=0,n_energy
          evec_e0(ie,ix)   = lag_moment(ie,alpha,0.5*(1+beta))
          evec_e05(ie,ix)  = lag_moment(ie,alpha,0.5*(2+beta))
          evec_e1(ie,ix)   = lag_moment(ie,alpha,0.5*(3+beta))
          evec_e105(ie,ix) = lag_moment(ie,alpha,0.5*(4+beta))
          evec_e2(ie,ix)   = lag_moment(ie,alpha,0.5*(5+beta))
       enddo
    enddo

    emat_e0(:,:,:,:)    = 0.0
    emat_e05(:,:,:,:)   = 0.0
    emat_en05(:,:,:,:)  = 0.0
    emat_e05de(:,:,:,:) = 0.0
    emat_e1(:,:,:,:)    = 0.0

    do ix=0,n_xi

       ! q=2, jx=ix: orthogonality, (1/2) Gamma(n+alpha+1)/n!
       alpha = 0.5*e_lag(ix)
       do ie=0,n_energy
          emat_e0(ie,ie,ix,1) = 0.5*exp(log_gamma(ie+alpha+1.0)-log_gamma(ie+1.0))
       enddo

       ! jx = ix-1, ix+1
       do kx=1,2
          jx = ix+2*kx-3
          if (jx >= 0 .and. jx <= n_xi) then
             call lag_quad(ix,jx,3,.false.,emat_e05(:,:,ix,kx))
             call lag_quad(ix,jx,1,.false.,emat_en05(:,:,ix,kx))
             call lag_quad(ix,jx,1,.true.,emat_e05de(:,:,ix,kx))
          endif
       enddo

       ! jx = ix-2, ix, ix+2
       do kx=1,3
          jx = ix+2*kx-4
          if (jx >= 0 .and. jx <= n_xi) then
             call lag_quad(ix,jx,4,.false.,emat_e1(:,:,ix,kx))
          endif
       enddo

    enddo

  end subroutine ENERGY_basis_ints

  ! m(i,j) = int_0^inf x^q e^{-x^2} phi_i(x) psi_j(x) dx, with psi_j = phi_j,
  ! or psi_j = x d(phi_j)/dx if deriv
  subroutine lag_quad(ix,jx,q,deriv,m)
    use neo_globals, only : n_energy
    implicit none
    integer, intent(in) :: ix, jx, q
    logical, intent(in) :: deriv
    real, dimension(0:n_energy,0:n_energy), intent(out) :: m
    integer :: nq, je
    real :: ai, aj, bj, gam
    real, dimension(:), allocatable :: s, w
    real, dimension(:,:), allocatable :: li, lj, pj

    nq = n_energy+2
    ai = 0.5*e_lag(ix)
    aj = 0.5*e_lag(jx)
    bj = 1.0*xi_beta_l(jx)
    gam = 0.5*(q+xi_beta_l(ix)+xi_beta_l(jx)-1)

    allocate(s(nq),w(nq),li(0:n_energy,nq),lj(0:n_energy,nq),pj(0:n_energy,nq))
    call gauss_laguerre(nq,gam,s,w)
    call lag_eval(n_energy,ai,nq,s,li)
    call lag_eval(n_energy,aj,nq,s,lj)

    if (deriv) then
       ! x d/dx [x^b L_n] = x^b [(b+2n) L_n - 2(n+a) L_{n-1}]
       pj(0,:) = bj*lj(0,:)
       do je=1,n_energy
          pj(je,:) = (bj+2*je)*lj(je,:)-2*(je+aj)*lj(je-1,:)
       enddo
    else
       pj = lj
    endif

    m = 0.5*proj(li,pj,w)

    deallocate(s,w,li,lj,pj)

  end subroutine lag_quad

  !===============================================================
  ! Collision integrals
  !===============================================================

  subroutine ENERGY_coll_ints(ir)
    use neo_globals, only : n_species, n_energy, dens, vth, mass, temp, &
         collision_model, spitzer_model, coll_uncoupledei_model, ae_flag, is_ele, &
         write_cmoments_flag
    implicit none
    integer, intent(in) :: ir
    real :: tauinv_ab, tauinv_ba, lambda
    integer :: is, js, ix, ia, ib
    integer :: model, ccase

    emat_coll_test(:,:,:,:,:)  = 0.0
    emat_coll_field(:,:,:,:,:) = 0.0

    if (collision_model /= 1 .and. collision_model /= 4 .and. &
         collision_model /= 5) then
       call neo_error('ERROR: (NEO) collision_model must be 1, 4 or 5')
       return
    endif

    do is=1, n_species
       do js=1, n_species

          ! (Note: pol part of dens from rotation will be added
          !  in coll term in kinetic equation)
          call pair_rates(ir,is,js,tauinv_ab,tauinv_ba,lambda)

          ccase = 0
          if(collision_model == 1 .or. &
               (spitzer_model==1 .and. js .ne. is)) then
             ! Connor model
             model = 1
             ccase = connor_case(is,js)
          else
             model = collision_model
          endif

          call coll_test_ints(is,js,ir,tauinv_ab,lambda,model,ccase)
          call coll_field_ints(is,js,ir,tauinv_ab,tauinv_ba,lambda,model,ccase)

       enddo
    enddo

    if(collision_model==5) then
       ! Replace actual field-particle op with ad-hoc op
       emat_coll_field(:,:,:,:,:) = 0.0
       do is=1, n_species
          do js=1, n_species

             ix = 1
             ia = (1-xi_beta_l(ix))/e_alpha
             if(abs(emat_coll_test(is,js,ia,ia,ix)) > epsilon(0.)) then
                do ib=0,n_energy
                   emat_coll_field(is,js,:,ib,ix) = &
                        -(mass(js)*dens(js,ir)*vth(js,ir)) &
                        / (mass(is)*dens(is,ir)*vth(is,ir)) &
                        * emat_coll_test(is,js,:,ia,ix) &
                        * emat_coll_test(js,is,ia,ib,ix) &
                        / emat_coll_test(is,js,ia,ia,ix)
                enddo
             endif

             ix = 0
             ia = (2-xi_beta_l(ix))/e_alpha
             if(abs(emat_coll_test(is,js,ia,ia,ix)) > epsilon(0.)) then
                do ib=0,n_energy
                   emat_coll_field(is,js,:,ib,ix) = &
                        -(dens(js,ir)*temp(js,ir)) &
                        / (dens(is,ir)*temp(is,ir)) &
                        * emat_coll_test(is,js,:,ia,ix) &
                        * emat_coll_test(js,is,ia,ib,ix) &
                        / emat_coll_test(is,js,ia,ia,ix)
                enddo
             endif

          enddo
       enddo
    endif

    if(coll_uncoupledei_model == 1 .or. coll_uncoupledei_model == 2) then
       if(ae_flag == 1) then
          call neo_error('ERROR: (NEO) Must have electron species for uncoupled e-i problem')
          return
       endif

       do is=1,n_species
          do js=1,n_species

             call pair_rates(ir,is,js,tauinv_ab,tauinv_ba,lambda)

             if(is .ne. is_ele .and. js == is_ele) then
                emat_coll_test(is,js,:,:,:) = 0.0
                emat_coll_field(is,js,:,:,:) = 0.0
                if(coll_uncoupledei_model == 2) then
                   ! C_ie = Connor case ma>mb
                   call coll_test_ints(is,js,ir,tauinv_ab,lambda,1,3)
                   call coll_field_ints(is,js,ir,tauinv_ab,tauinv_ba,lambda,1,3)
                endif
             endif
             if(is == is_ele .and. js .ne. is) then
                ! C_ei = L(fe) with Connor case ma<mb
                emat_coll_test(is,js,:,:,:) = 0.0
                emat_coll_field(is,js,:,:,:) = 0.0
                call coll_test_ints(is,js,ir,tauinv_ab,lambda,1,2)
                call coll_field_ints(is,js,ir,tauinv_ab,tauinv_ba,lambda,1,2)
             endif
          enddo
       enddo
    endif

    if (write_cmoments_flag == 1) then
       call write_fullcoll_mono(ir)
    endif

  end subroutine ENERGY_coll_ints

  subroutine pair_rates(ir,is,js,tauinv_ab,tauinv_ba,lambda)
    use neo_globals, only : nu, Z, dens, vth
    implicit none
    integer, intent(in) :: ir, is, js
    real, intent(out) :: tauinv_ab, tauinv_ba, lambda

    tauinv_ab = nu(is,ir) * (1.0*Z(js))**2 / (1.0*Z(is))**2 &
         * dens(js,ir)/dens(is,ir)
    tauinv_ba = nu(js,ir) * (1.0*Z(is))**2 / (1.0*Z(js))**2 &
         * dens(is,ir)/dens(js,ir)
    lambda = (vth(is,ir) / vth(js,ir))**2

  end subroutine pair_rates

  ! Connor model: 1 for ma=mb, 2 for ma<mb, 3 for ma>mb
  integer function connor_case(is,js)
    use neo_globals, only : mass
    implicit none
    integer, intent(in) :: is, js

    if(is == js .or. (abs(mass(is)-mass(js)) < epsilon(0.))) then
       connor_case = 1
    else if(mass(is) < mass(js)) then
       connor_case = 2
    else
       connor_case = 3
    endif

  end function connor_case

  ! True if ix starts a new basis (e_lag, xi_beta_l differ from ix-1)
  logical function new_basis(ix)
    implicit none
    integer, intent(in) :: ix
    new_basis = .true.
    if (ix > 0) new_basis = (e_lag(ix) /= e_lag(ix-1) .or. &
         xi_beta_l(ix) /= xi_beta_l(ix-1))
  end function new_basis

  !---------------------------------------------------------------
  ! Test-particle matrices emat_coll_test(is,js,:,:,ix), all ix.
  ! k0, k2 as in kfun; D phi = x dphi/dx.
  !   Connor 1:   -c0 l(l+1) <phi phi (k0-k2)>
  !   Connor 2:   -tau/2 l(l+1) <phi phi /x>
  !   Connor 3:   -c0 (Ta/Tb) (2/3) l(l+1) <phi phi x^2>
  !   FP:         -c0 l(l+1) <phi phi (k0-k2)>
  !               + 4 c0 [(1-Ta/Tb) <D phi phi x^2 k2> - 1/2 <D phi D phi k2>]
  ! with c0 = tau sqrt(lambda/pi) and <..> = int dx e^{-x^2} (..).
  ! The matrices depend on l only through l(l+1) and the basis, so
  ! they are recomputed only when the basis changes.
  !---------------------------------------------------------------
  subroutine coll_test_ints(is,js,ir,tauinv_ab,lambda,model,ccase)
    use neo_globals, only : n_energy, n_xi, pi, temp
    implicit none
    integer, intent(in) :: is, js, ir, model, ccase
    real, intent(in) :: tauinv_ab, lambda
    real, dimension(:), allocatable :: x, w, k0, k2
    real, dimension(:,:), allocatable :: phi, dphi, tpa, ten
    real :: c0, ta
    integer :: nq, ix, q, ne

    ne = n_energy
    call x_grid(lambda,x,w)
    nq = size(x)
    allocate(k0(nq),k2(nq))
    allocate(phi(0:ne,nq),dphi(0:ne,nq),tpa(0:ne,0:ne),ten(0:ne,0:ne))
    do q=1,nq
       k0(q) = kfun(0,lambda*x(q)**2)
       k2(q) = kfun(2,lambda*x(q)**2)
    enddo

    c0 = tauinv_ab*sqrt(lambda/pi)
    ta = temp(is,ir)/temp(js,ir)

    do ix=0,n_xi

       if (new_basis(ix)) then

          call basis_eval(ne,0.5*e_lag(ix),xi_beta_l(ix),nq,x,phi,dphi)
          ten(:,:) = 0.0

          if (model == 1) then
             select case (ccase)
             case (1)
                tpa = c0*proj(phi,phi,w*(k0-k2))
             case (2)
                tpa = 0.5*tauinv_ab*proj(phi,phi,w/x)
             case (3)
                tpa = c0*ta*(2.0/3.0)*proj(phi,phi,w*x**2)
             end select
          else
             tpa = c0*proj(phi,phi,w*(k0-k2))
             ten = 4*c0*((1-ta)*proj(dphi,phi,w*x**2*k2) &
                  -0.5*proj(dphi,dphi,w*k2))
          endif

       endif

       emat_coll_test(is,js,:,:,ix) = -ix*(ix+1)*tpa+ten

    enddo

    deallocate(x,w,k0,k2,phi,dphi,tpa,ten)

  end subroutine coll_test_ints

  !---------------------------------------------------------------
  ! Field-particle matrices emat_coll_field(is,js,:,:,ix), all ix.
  ! x is the test-species (a) speed in units of v_ta; the field
  ! species (b) basis is phi_j(sqrt(lambda) x).
  ! Connor: l=1 only, separable, with kernels K_n(x;lambda) (k0,k2)
  ! and K_n(x;1/lambda) (q0,q2) as in kfun.  Full FP: fp_field.
  !---------------------------------------------------------------
  subroutine coll_field_ints(is,js,ir,tauinv_ab,tauinv_ba,lambda,model,ccase)
    use neo_globals, only : n_energy, n_xi, pi, temp, mass, dens, vth
    implicit none
    integer, intent(in) :: is, js, ir, model, ccase
    real, intent(in) :: tauinv_ab, tauinv_ba, lambda
    real, dimension(:), allocatable :: x, w, wa, k0, k2, q0, q2
    real, dimension(:,:), allocatable :: phi, dphi
    real, dimension(0:n_energy) :: va, vb
    real :: alpha, beta, tb, mm, xmax
    integer :: nq, ie, q, ne

    emat_coll_field(is,js,:,:,:) = 0.0

    if (model == 5) return

    if (model == 4) then
       call fp_field(is,js,tauinv_ab,lambda,mass(is)/mass(js))
       return
    endif

    ! Connor: l=1 only
    if (n_xi < 1) return

    ne    = n_energy
    tb    = temp(js,ir)/temp(is,ir)
    mm    = (mass(js)*dens(js,ir)*vth(js,ir))/(mass(is)*dens(is,ir)*vth(is,ir))
    alpha = 0.5*e_lag(1)
    beta  = 1.0*xi_beta_l(1)

    select case (ccase)

    case (1)

       call field_grid(lambda,x,w)
       nq = size(x)
       xmax = x_extent()
       allocate(wa(nq),k0(nq),k2(nq),q0(nq),q2(nq))
       allocate(phi(0:ne,nq),dphi(0:ne,nq))
       wa = 0.0; k0 = 0.0; k2 = 0.0; q0 = 0.0; q2 = 0.0
       do q=1,nq
          if (x(q) > xmax) cycle
          wa(q) = w(q)*exp(-x(q)**2)
          k0(q) = kfun(0,lambda*x(q)**2)
          k2(q) = kfun(2,lambda*x(q)**2)
          q0(q) = kfun(0,x(q)**2/lambda)
          q2(q) = kfun(2,x(q)**2/lambda)
       enddo
       call basis_eval(ne,alpha,xi_beta_l(1),nq,merge(x,0.0,x <= xmax),phi,dphi)
       va = matmul(phi,wa*x*(k0-k2))
       vb = matmul(phi,wa*x*(q0-q2))
       emat_coll_field(is,js,:,:,1) = mm*2*tauinv_ba/sqrt(lambda*pi) &
            *outer(va,vb)/sum(wa*x**2*(k0-k2))
       deallocate(x,w,wa,k0,k2,q0,q2,phi,dphi)

    case (2)

       do ie=0,ne
          va(ie) = 2*basis_moment(ie,alpha,0.5*(beta-1))
          vb(ie) = 2*basis_moment(ie,alpha,0.5*(beta+2))
       enddo
       emat_coll_field(is,js,:,:,1) = mm*(2.0/3.0)*tauinv_ba*tb/sqrt(lambda*pi) &
            *outer(va,vb)

    case (3)

       do ie=0,ne
          va(ie) = 2*basis_moment(ie,alpha,0.5*(beta+2))
          vb(ie) = 2*basis_moment(ie,alpha,0.5*(beta-1))
       enddo
       emat_coll_field(is,js,:,:,1) = mm*0.5*tauinv_ba*tb/gamma(2.5)*outer(va,vb)

    end select

  end subroutine coll_field_ints

  !---------------------------------------------------------------
  ! Full linearized FP field-particle operator
  !
  !   C = tau_ab 2/sqrt(pi) lambda^{3/2} lambda^{beta/2}
  !       [ m_a/m_b R + w2 D(1,2) + w3 Db(1,2) + w45 (D(3,2) + Db(1,4))
  !                   + w67 (D(1,4) + Db(3,2)) ]
  !   R       = 2 int x^{2 beta+2} e^{-(1+lambda)x^2} L_i(x^2) L_j(lambda x^2) dx
  !   D(A,B)  = int dx x^{2 beta+A+B+1} e^{-x^2} L_i(x^2)
  !               x^{-(p+1)} int_0^x y^p e^{-lambda y^2} L_j(lambda y^2) dy
  !   Db(A,B) = int dy y^{2 beta+A+B+1} e^{-lambda y^2} L_j(lambda y^2)
  !               y^{-(p+1)} int_0^y x^p e^{-x^2} L_i(x^2) dx
  ! with p = l+B+beta and the numbers in D(A,B) giving the power
  ! offsets.  The scaled inner integrals
  !   I_n(x;c,p) = int_0^1 t^p e^{-u} L_n(u) dt,  u = c x^2 t^2,
  ! come from one Gauss-Legendre rule in t: for each outer node the
  ! table e^{-u} L_n(u) is built once and contracted with t^p for
  ! every p the basis group needs.  Beyond u = scut the tail of
  ! int_0^inf is negligible and the closed form is used.
  !---------------------------------------------------------------
  subroutine fp_field(is,js,tauinv_ab,lambda,mab)
    use neo_globals, only : n_energy, n_xi, pi
    implicit none
    integer, intent(in) :: is, js
    real, intent(in) :: tauinv_ab, lambda, mab
    real, dimension(:), allocatable :: x, w, wa, wb, tg, tw
    logical, dimension(:), allocatable :: ma, mb
    real, dimension(:,:), allocatable :: la, lb, tp
    real, dimension(:,:,:), allocatable :: f, ita, itb
    integer, dimension(:), allocatable :: plist
    integer, dimension(:), allocatable :: pidx
    real :: xmax, scut, cf, beta, alpha, el, w2, w3, w45, w67
    integer :: ne, nq, nt, pmax, np, ix, ix0, ix1, b0, b1, nb, k, p, ip1, ip2
    real, dimension(:), allocatable :: va1, va2, vb1, vb2, xb4, xb6

    ne = n_energy

    call field_grid(lambda,x,w)
    nq = size(x)
    xmax = x_extent()

    allocate(ma(nq),mb(nq),wa(nq),wb(nq))
    ma = x <= xmax
    mb = sqrt(lambda)*x <= xmax
    wa = 0.0
    wb = 0.0
    where (ma) wa = w*exp(-x*x)
    where (mb) wb = w*exp(-lambda*x*x)

    ! t-rule: exact for t^p g(t) with p <= pmax whenever g is a
    ! polynomial of degree < 2*(nt - (pmax+1)/2)
    pmax = n_xi+4+maxval(xi_beta_l)
    scut = 10.0*ne+2*pmax+80
    nt   = ne+ceiling(3*sqrt(scut))+12+(pmax+2)/2
    allocate(tg(nt),tw(nt))
    call gauss_legendre01(nt,tg,tw)

    allocate(la(0:ne,nq),lb(0:ne,nq),f(0:ne,0:ne,0:n_xi))
    allocate(pidx(0:pmax),xb4(nq),xb6(nq))
    allocate(va1(nblk),va2(nblk),vb1(nblk),vb2(nblk))
    f(:,:,:) = 0.0

    ix0 = 0
    do while (ix0 <= n_xi)

       ! basis group ix0..ix1 (same alpha, beta)
       ix1 = ix0
       do while (ix1 < n_xi)
          if (new_basis(ix1+1)) exit
          ix1 = ix1+1
       enddo
       alpha = 0.5*e_lag(ix0)
       beta  = 1.0*xi_beta_l(ix0)
       cf    = lambda**(0.5*beta)

       ! test-species basis on x <= xmax, field-species on sqrt(lambda) x <= xmax
       call lag_eval(ne,alpha,nq,merge(x*x,0.0,ma),la)
       call lag_eval(ne,alpha,nq,merge(lambda*x*x,0.0,mb),lb)
       do k=0,ne
          la(k,:) = merge(la(k,:),0.0,ma)
          lb(k,:) = merge(lb(k,:),0.0,mb)
       enddo
       xb4 = x**(4+2*beta)
       xb6 = x**(6+2*beta)

       ! R: local term
       f(:,:,ix0) = mab*2*cf*proj(la,lb, &
            merge(w*exp(-(1+lambda)*x*x),0.0,ma.and.mb)*x**(2*beta+2))
       do ix=ix0+1,ix1
          f(:,:,ix) = f(:,:,ix0)
       enddo

       ! powers p needed by this group: l+2+beta, l+4+beta
       pidx(:) = 0
       do ix=ix0,ix1
          pidx(ix+2+nint(beta)) = 1
          pidx(ix+4+nint(beta)) = 1
       enddo
       np = count(pidx == 1)
       allocate(plist(np),tp(nt,np))
       np = 0
       do p=0,pmax
          if (pidx(p) == 1) then
             np = np+1
             plist(np) = p
             pidx(p) = np
             tp(:,np) = tw*tg**p
          endif
       enddo
       allocate(ita(0:ne,nblk,np),itb(0:ne,nblk,np))

       ! outer nodes in blocks
       do b0=1,nq,nblk
          b1 = min(nq,b0+nblk-1)
          nb = b1-b0+1

          call inner_ints(ne,alpha,lambda,nt,tg,np,plist,tp,scut, &
               nb,x(b0:b1),ma(b0:b1),ita)
          call inner_ints(ne,alpha,1.0,nt,tg,np,plist,tp,scut, &
               nb,x(b0:b1),mb(b0:b1),itb)

          do ix=ix0,ix1
             el  = 1.0*ix
             w2  = -2.0/(el+0.5)*(mab-el*(1-mab))
             w3  = -2.0/(el+0.5)*(1+el*(1-mab))
             w67 = (el+1.0)*(el+2.0)/((el+1.5)*(el+0.5))
             if (ix > 1) then
                w45 = -el*(el-1.0)/(el*el-0.25)
             else
                w45 = 0.0
             endif
             ip1 = pidx(ix+2+nint(beta))
             ip2 = pidx(ix+4+nint(beta))

             va1(1:nb) = wa(b0:b1)*(w2*xb4(b0:b1)+w45*xb6(b0:b1))
             va2(1:nb) = wa(b0:b1)*w67*xb6(b0:b1)
             vb1(1:nb) = wb(b0:b1)*(w3*xb4(b0:b1)+w45*xb6(b0:b1))
             vb2(1:nb) = wb(b0:b1)*w67*xb6(b0:b1)

             f(:,:,ix) = f(:,:,ix) + cf*( &
                  proj(la(:,b0:b1),ita(:,1:nb,ip1),va1(1:nb)) &
                  +proj(la(:,b0:b1),ita(:,1:nb,ip2),va2(1:nb)) &
                  +proj(itb(:,1:nb,ip1),lb(:,b0:b1),vb1(1:nb)) &
                  +proj(itb(:,1:nb,ip2),lb(:,b0:b1),vb2(1:nb)))
          enddo
       enddo

       deallocate(plist,tp,ita,itb)
       ix0 = ix1+1

    enddo

    emat_coll_field(is,js,:,:,:) = tauinv_ab*2/sqrt(pi)*lambda**1.5*f

    deallocate(x,w,wa,wb,ma,mb,tg,tw,la,lb,f,pidx,xb4,xb6,va1,va2,vb1,vb2)

  end subroutine fp_field

  ! it(n,q,k) = int_0^1 t^p e^{-u} L_n^(alpha)(u) dt, u = c x_q^2 t^2,
  ! p = plist(k), for the nodes x_q with mask(q) (zero elsewhere).
  ! tp(:,k) = Gauss-Legendre weights times t^p.  Beyond u = scut:
  !   int_0^inf y^p e^{-c y^2} L_n(c y^2) dy = c^{-(p+1)/2} basis_moment(n,alpha,(p-1)/2)
  subroutine inner_ints(ne,alpha,c,nt,tg,np,plist,tp,scut,nx,x,mask,it)
    implicit none
    integer, intent(in) :: ne, nt, np, nx
    real, intent(in) :: alpha, c, scut
    real, intent(in) :: tg(nt), tp(nt,np), x(nx)
    integer, intent(in) :: plist(np)
    logical, intent(in) :: mask(nx)
    real, intent(inout) :: it(0:ne,nblk,np)
    real :: u(nt), l(0:ne,nt), full(0:ne,np), s
    integer :: q, k, n

    do k=1,np
       do n=0,ne
          full(n,k) = c**(-0.5*(plist(k)+1))*basis_moment(n,alpha,0.5*(plist(k)-1))
       enddo
    enddo

    do q=1,nx
       if (.not. mask(q)) then
          it(:,q,:) = 0.0
          cycle
       endif
       s = c*x(q)**2
       if (s > scut) then
          do k=1,np
             it(:,q,k) = full(:,k)/x(q)**(plist(k)+1)
          enddo
       else
          u = s*tg**2
          call lag_eval(ne,alpha,nt,u,l)
          do n=0,ne
             l(n,:) = l(n,:)*exp(-u)
          enddo
          it(:,q,:) = matmul(l,tp)
       endif
    enddo

  end subroutine inner_ints

  !---------------------------------------------------------------
  ! Diagnostic: full FP operator (test and field parts) in the
  ! monomial basis x^(2n+beta), written to out.neo.diagnostic_coll*
  !---------------------------------------------------------------
  subroutine write_fullcoll_mono(ir)
    use neo_globals, only : n_species, n_energy, n_xi, path
    implicit none
    integer, intent(in) :: ir
    real :: tauinv_ab, tauinv_ba, lambda
    integer :: is, js, ie, je, ix
    real, dimension(:,:,:,:,:), allocatable :: save_test, save_field
    real, dimension(:,:), allocatable :: ctot

    allocate(save_test,source=emat_coll_test)
    allocate(save_field,source=emat_coll_field)
    allocate(ctot(0:n_energy,0:n_energy))

    mono = .true.
    do is=1, n_species
       do js=1, n_species
          call pair_rates(ir,is,js,tauinv_ab,tauinv_ba,lambda)
          call coll_test_ints(is,js,ir,tauinv_ab,lambda,4,0)
          call coll_field_ints(is,js,ir,tauinv_ab,tauinv_ba,lambda,4,0)
       enddo
    enddo
    mono = .false.

    open(unit=1,file=trim(path)//'out.neo.diagnostic_coll',status='replace')
    write(1,'(a)') '# Collision matrix energy-element order [C00,C01,C02,C11,C12,C22,...]'
    write(1,'(a)') '# '
    close(1)
    open(unit=1,file=trim(path)//'out.neo.diagnostic_colltest',status='replace')
    write(1,'(a)') '# Collision matrix energy-element order [C00,C01,C02,C10,C11,C12,C20,C21,C22,...]'
    write(1,'(a)') '# '
    close(1)
    open(unit=1,file=trim(path)//'out.neo.diagnostic_collfield',status='replace')
    write(1,'(a)') '# Collision matrix energy-element order [C00,C01,C02,C10,C11,C12,C20,C21,C22,...]'
    write(1,'(a)') '# '
    close(1)

    do is=1, n_species
       do js=1, n_species

          call pair_rates(ir,is,js,tauinv_ab,tauinv_ba,lambda)

          open(unit=1,file=trim(path)//'out.neo.diagnostic_coll',status='old',position='append')
          do ix=0,n_xi
             ctot(:,:) = (emat_coll_test(is,js,:,:,ix)+emat_coll_field(is,js,:,:,ix))/tauinv_ab
             write(1,'(a,i2,a,i2,a,i2,a)') '# ix = ',ix,'  (is,js)=(',is,',',js,')'
             do ie=0,n_energy
                do je=ie,n_energy
                   write(1,'(1pe16.8,1x)',advance='no') ctot(ie,je)
                enddo
             enddo
             write(1,*)
          enddo
          close(1)
          open(unit=1,file=trim(path)//'out.neo.diagnostic_colltest',status='old',position='append')
          do ix=0,n_xi
             ctot(:,:) = emat_coll_test(is,js,:,:,ix)/tauinv_ab
             write(1,'(a,i2,a,i2,a,i2,a)') '# ix = ',ix,'  (is,js)=(',is,',',js,')'
             do ie=0,n_energy
                do je=0,n_energy
                   write(1,'(1pe16.8,1x)',advance='no') ctot(ie,je)
                enddo
             enddo
             write(1,*)
          enddo
          close(1)
          open(unit=1,file=trim(path)//'out.neo.diagnostic_collfield',status='old',position='append')
          do ix=0,n_xi
             ctot(:,:) = emat_coll_field(is,js,:,:,ix)/tauinv_ab
             write(1,'(a,i2,a,i2,a,i2,a)') '# ix = ',ix,'  (is,js)=(',is,',',js,')'
             do ie=0,n_energy
                do je=0,n_energy
                   write(1,'(1pe16.8,1x)',advance='no') ctot(ie,je)
                enddo
             enddo
             write(1,*)
          enddo
          close(1)

       enddo
    enddo

    emat_coll_test  = save_test
    emat_coll_field = save_field
    deallocate(save_test,save_field,ctot)

  end subroutine write_fullcoll_mono

  !===============================================================
  ! Basis functions
  !===============================================================

  ! L_n^(alpha)(s), n=0..nmax, by the three-term recurrence
  ! (s^n in monomial mode)
  subroutine lag_eval(nmax,alpha,ns,s,l)
    implicit none
    integer, intent(in) :: nmax, ns
    real, intent(in) :: alpha
    real, dimension(ns), intent(in) :: s
    real, dimension(0:nmax,ns), intent(out) :: l
    integer :: n

    l(0,:) = 1.0
    if (mono) then
       do n=1,nmax
          l(n,:) = l(n-1,:)*s
       enddo
       return
    endif
    if (nmax > 0) l(1,:) = 1.0+alpha-s
    do n=1,nmax-1
       l(n+1,:) = ((2*n+1+alpha-s)*l(n,:)-(n+alpha)*l(n-1,:))/(n+1)
    enddo

  end subroutine lag_eval

  ! phi_n = x^beta L_n(x^2) and D phi_n = x dphi_n/dx at the points x
  subroutine basis_eval(ne,alpha,ibeta,nx,x,phi,dphi)
    implicit none
    integer, intent(in) :: ne, ibeta, nx
    real, intent(in) :: alpha
    real, intent(in) :: x(nx)
    real, intent(out) :: phi(0:ne,nx), dphi(0:ne,nx)
    real :: l(0:ne,nx), xb(nx)
    integer :: n

    call lag_eval(ne,alpha,nx,x*x,l)
    xb = x**ibeta
    do n=0,ne
       phi(n,:) = xb*l(n,:)
    enddo
    if (mono) then
       do n=0,ne
          dphi(n,:) = (ibeta+2*n)*phi(n,:)
       enddo
    else
       ! x d/dx [x^b L_n] = x^b [(b+2n) L_n - 2(n+a) L_{n-1}]
       dphi(0,:) = ibeta*phi(0,:)
       do n=1,ne
          dphi(n,:) = xb*((ibeta+2*n)*l(n,:)-2*(n+alpha)*l(n-1,:))
       enddo
    endif

  end subroutine basis_eval

  ! (1/2) int_0^inf s^gam e^{-s} L_n^(alpha)(s) ds
  !   = (1/2) Gamma(gam+1) (alpha-gam)_n / n!
  real function lag_moment(n,alpha,gam)
    implicit none
    integer, intent(in) :: n
    real, intent(in) :: alpha, gam
    integer :: k

    lag_moment = 0.5*gamma(gam+1.0)
    do k=0,n-1
       lag_moment = lag_moment*(alpha-gam+k)/(k+1.0)
    enddo

  end function lag_moment

  ! As lag_moment, for the current basis (Laguerre or monomial)
  real function basis_moment(n,alpha,gam)
    implicit none
    integer, intent(in) :: n
    real, intent(in) :: alpha, gam

    if (mono) then
       basis_moment = 0.5*gamma(gam+n+1.0)
    else
       basis_moment = lag_moment(n,alpha,gam)
    endif

  end function basis_moment

  ! k_n(z) = K_n(x)/x^{n+1}, K_n(x) = int_0^x y^n e^{-lambda y^2} dy,
  ! z = lambda x^2, for n = 0, 2:
  !   k_n = e^{-z}/(n+1) M(1,(n+3)/2,z)    (positive series, z < 2)
  !   k_0 = sqrt(pi) erf(sqrt z)/(2 sqrt z),  k_2 = (k_0-e^{-z})/(2z)
  real function kfun(n,z)
    use neo_globals, only : pi
    implicit none
    integer, intent(in) :: n
    real, intent(in) :: z
    real :: b, term, k0
    integer :: k

    if (z < 2.0) then
       b = 0.5*(n+3)
       term = 1.0
       kfun = 1.0
       k = 0
       do while (term > epsilon(term)*kfun)
          term = term*z/(b+k)
          kfun = kfun+term
          k = k+1
       enddo
       kfun = exp(-z)*kfun/(n+1)
    else
       k0 = sqrt(pi)*erf(sqrt(z))/(2*sqrt(z))
       if (n == 0) then
          kfun = k0
       else
          kfun = (k0-exp(-z))/(2*z)
       endif
    endif

  end function kfun

  !===============================================================
  ! Quadrature
  !===============================================================

  ! Extent in x of e^{-x^2} x^{2 beta} L_n^2 for all basis functions
  real function x_extent()
    use neo_globals, only : n_energy
    implicit none
    x_extent = sqrt(6.0*n_energy+maxval(e_lag)+60.0)
  end function x_extent

  ! Composite Gauss-Legendre grid in x on [0,x_max], with weights times
  ! e^{-x^2}.  For lambda > 1 the panels double from hpan/sqrt(lambda) to
  ! resolve the field-species layer at x ~ lambda^{-1/2}; then width hpan.
  subroutine x_grid(lambda,x,w)
    implicit none
    real, intent(in) :: lambda
    real, dimension(:), allocatable, intent(out) :: x, w
    real, dimension(ng) :: t, wt
    real, dimension(:), allocatable :: b
    real :: xmax, a
    integer :: np, k

    xmax = x_extent()

    allocate(b(0:ceiling(xmax/hpan)+64))
    b(0) = 0.0
    np = 0
    if (lambda > 1.0) then
       a = hpan/sqrt(lambda)
       do while (a < hpan)
          np = np+1
          b(np) = a
          a = 2*a
       enddo
    endif
    do while (b(np) < xmax)
       np = np+1
       b(np) = b(np-1)+hpan
    enddo

    call gauss_legendre01(ng,t,wt)
    allocate(x(np*ng),w(np*ng))
    do k=1,np
       x((k-1)*ng+1:k*ng) = b(k-1)+(b(k)-b(k-1))*t
       w((k-1)*ng+1:k*ng) = (b(k)-b(k-1))*wt
    enddo
    w = w*exp(-x*x)
    deallocate(b)

  end subroutine x_grid

  ! Merged composite Gauss-Legendre grid on [0,inf): panels of width
  ! hpan*s, s = 1, lambda^{-1/2}, lambda^{1/2}, each sequence to x_max*s, so
  ! that functions of x, sqrt(lambda) x and x/sqrt(lambda) are all resolved.
  ! Between the smallest and largest scale the panels also double in width,
  ! for the power-law tails x^{-k} of the Rosenbluth potentials there.
  ! Raw weights (no Maxwellian).
  subroutine field_grid(lambda,x,w)
    implicit none
    real, intent(in) :: lambda
    real, dimension(:), allocatable, intent(out) :: x, w
    real, dimension(ng) :: t, wt
    real, dimension(:), allocatable :: b
    real :: xmax, s(3), v
    integer :: kmax, i, k, np, nb

    xmax = x_extent()
    kmax = ceiling(xmax/hpan)
    s = (/ 1.0, 1.0/sqrt(lambda), sqrt(lambda) /)

    allocate(b(3*(kmax+1)+200))
    nb = 0
    do i=1,3
       do k=0,kmax
          nb = nb+1
          b(nb) = hpan*s(i)*k
       enddo
    enddo
    if (lambda /= 1.0) then
       v = hpan*minval(s)
       do while (v < hpan*maxval(s))
          nb = nb+1
          b(nb) = v
          v = 2*v
       enddo
    endif
    call sort_unique(b,nb)

    call gauss_legendre01(ng,t,wt)
    np = nb-1
    allocate(x(np*ng),w(np*ng))
    do k=1,np
       x((k-1)*ng+1:k*ng) = b(k)+(b(k+1)-b(k))*t
       w((k-1)*ng+1:k*ng) = (b(k+1)-b(k))*wt
    enddo
    deallocate(b)

  contains

    subroutine sort_unique(a,n)
      real, dimension(:), intent(inout) :: a
      integer, intent(inout) :: n
      integer :: i, j, m
      real :: v
      do i=2,n
         v = a(i)
         j = i-1
         do while (j >= 1)
            if (a(j) <= v) exit
            a(j+1) = a(j)
            j = j-1
         enddo
         a(j+1) = v
      enddo
      m = 1
      do i=2,n
         if (a(i)-a(m) > 1e-9*a(i)) then
            m = m+1
            a(m) = a(i)
         endif
      enddo
      n = m
    end subroutine sort_unique

  end subroutine field_grid

  ! Gauss-Laguerre rule for weight s^gam e^{-s} on [0,inf).
  ! Jacobi matrix: diagonal 2k+gam+1, off-diagonal sqrt(k(k+gam)).
  subroutine gauss_laguerre(n,gam,s,w)
    implicit none
    integer, intent(in) :: n
    real, intent(in) :: gam
    real, dimension(n), intent(out) :: s, w
    real, dimension(n) :: e
    integer :: k

    do k=1,n
       s(k) = 2*k-1+gam
       e(k) = 0.0
    enddo
    do k=1,n-1
       e(k) = sqrt(k*(k+gam))
    enddo
    call gauss_jacobi_rule(n,s,e,gamma(gam+1.0),w)

  end subroutine gauss_laguerre

  ! Gauss-Legendre rule on [0,1].
  ! Jacobi matrix on [-1,1]: diagonal 0, off-diagonal k/sqrt(4k^2-1).
  subroutine gauss_legendre01(n,t,w)
    implicit none
    integer, intent(in) :: n
    real, dimension(n), intent(out) :: t, w
    real, dimension(n) :: e
    integer :: k

    t(:) = 0.0
    e(:) = 0.0
    do k=1,n-1
       e(k) = k/sqrt(4.0*k*k-1.0)
    enddo
    call gauss_jacobi_rule(n,t,e,2.0,w)
    t = 0.5*(t+1.0)
    w = 0.5*w

  end subroutine gauss_legendre01

  ! Golub-Welsch: on input s is the diagonal and e the off-diagonal
  ! (e(k) couples k,k+1) of the Jacobi matrix, mu0 the total weight.
  ! On output s holds the nodes (eigenvalues) and w = mu0 v_1^2, with v_1
  ! the first component of each normalized eigenvector.  Implicit QL with
  ! Wilkinson shifts, carrying only the first row of the eigenvector matrix.
  subroutine gauss_jacobi_rule(n,s,e,mu0,w)
    implicit none
    integer, intent(in) :: n
    real, dimension(n), intent(inout) :: s, e
    real, intent(in) :: mu0
    real, dimension(n), intent(out) :: w
    real, dimension(n) :: z
    integer :: i, k, l, m, iter
    real :: b, c, f, g, p, r, sn, dd
    logical :: split

    z(:) = 0.0
    z(1) = 1.0

    do l=1,n
       iter = 0
       do
          m = n
          do k=l,n-1
             dd = abs(s(k))+abs(s(k+1))
             if (abs(e(k)) <= epsilon(dd)*dd) then
                m = k
                exit
             endif
          enddo
          if (m == l) exit
          iter = iter+1
          if (iter > 100) call neo_error('ERROR: (NEO) gauss_jacobi_rule failed')
          g  = (s(l+1)-s(l))/(2*e(l))
          r  = hypot(g,1.0)
          g  = s(m)-s(l)+e(l)/(g+sign(r,g))
          sn = 1.0
          c  = 1.0
          p  = 0.0
          split = .false.
          do i=m-1,l,-1
             f = sn*e(i)
             b = c*e(i)
             r = hypot(f,g)
             e(i+1) = r
             if (r == 0.0) then
                s(i+1) = s(i+1)-p
                e(m) = 0.0
                split = .true.
                exit
             endif
             sn = f/r
             c  = g/r
             g  = s(i+1)-p
             r  = (s(i)-g)*sn+2*c*b
             p  = sn*r
             s(i+1) = g+p
             g  = c*r-b
             f  = z(i+1)
             z(i+1) = sn*z(i)+c*f
             z(i)   = c*z(i)-sn*f
          enddo
          if (split) cycle
          s(l) = s(l)-p
          e(l) = g
          e(m) = 0.0
       enddo
    enddo

    w = mu0*z**2

  end subroutine gauss_jacobi_rule

  !===============================================================
  ! Small linear-algebra helpers
  !===============================================================

  ! P(i,j) = sum_q v_q a(i,q) b(j,q)
  function proj(a,b,v)
    real, dimension(:,:), intent(in) :: a, b
    real, dimension(:), intent(in) :: v
    real, dimension(size(a,1),size(b,1)) :: proj
    real, dimension(:,:), allocatable :: bv
    integer :: q

    allocate(bv(size(b,2),size(b,1)))
    do q=1,size(v)
       bv(q,:) = v(q)*b(:,q)
    enddo
    proj = matmul(a,bv)
    deallocate(bv)

  end function proj

  function outer(a,b)
    real, dimension(:), intent(in) :: a, b
    real, dimension(size(a),size(b)) :: outer
    integer :: i
    do i=1,size(a)
       outer(i,:) = a(i)*b(:)
    enddo
  end function outer

end module neo_energy_grid
