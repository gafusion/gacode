module neo_energy_grid

  implicit none

  public :: ENERGY_basis_ints_alloc, ENERGY_basis_ints, &
       ENERGY_coll_ints_alloc,ENERGY_coll_ints

  ! energy matrices
  ! (energy, xi, left/right diagonal in xi)
  real, dimension(:,:,:,:), allocatable :: emat_e05, emat_en05, emat_e05de, &
       emat_e0, emat_e1

  ! energy vectors
  real, dimension(:,:), allocatable :: evec_e0, evec_e1, evec_e2, evec_e05, evec_e105

  ! xi vectors -- order of associated Laguerre basis
  integer, dimension(:), allocatable  :: e_lag 
  integer, dimension(:), allocatable  :: xi_beta_l
  integer, parameter  :: e_alpha=2  ! should be 1 or 2

  ! collision matrices
  ! (species,species,energy,xi)
  real, dimension(:,:,:,:,:), allocatable :: emat_coll_test, emat_coll_field

  ! private variables
  logical, private :: initialized_basis = .false.
  logical, private :: initialized_coll  = .false.
  real, dimension(:), allocatable :: mygamma2

contains

  subroutine ENERGY_basis_ints_alloc(flag)

    use neo_globals, only : n_energy, n_xi, laguerre_method
    implicit none
    integer, intent (in) :: flag
    integer :: ie,ix
    integer :: xarg

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

       case (3)

          ! Laguerre 1/2
          e_lag(:) = 1
          xi_beta_l(:) = 0

       case (4)

          ! Laguerre 3/2
          e_lag(:) = 3
          xi_beta_l(:) = 1

       end select

       xarg = 4*e_alpha * n_energy + 4*n_xi + 12
       allocate(mygamma2(1:xarg))
       do ie=1,xarg
          mygamma2(ie) = gamma2(ie)
       enddo

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
       deallocate(mygamma2)

       initialized_basis = .false.

    endif

  end subroutine ENERGY_basis_ints_alloc
  
  subroutine ENERGY_basis_ints
    ! Energy integrals of the basis phi_n = x^beta L_n^(alpha)(x^2),
    ! alpha = e_lag/2, beta = xi_beta_l.  With s = x^2,
    !
    !   int_0^inf x^q e^{-x^2} phi_i phi_j dx = (1/2) int_0^inf s^gam e^{-s} L_i L_j ds
    !
    ! with gam = (q+beta_i+beta_j-1)/2.  Vectors and emat_e0 are closed form;
    ! the remaining matrices use Gauss-Laguerre quadrature with L_n evaluated
    ! by recurrence (exact for nq > n_energy).  No power-series expansion of
    ! L_n is used, so there is no cancellation.
    use neo_globals, only : n_energy, n_xi
    implicit none
    integer :: ie, ix, jx, kx
    real :: alpha, beta

    ! vectors: int_0^inf x^q e^{-x^2} phi_i dx
    do ix=0,n_xi
       alpha = 0.5*e_lag(ix)
       beta  = 1.0*xi_beta_l(ix)
       do ie=0,n_energy
          evec_e0(ie,ix)   = lag_moment(ie,alpha,0.5*(2+beta-1))
          evec_e05(ie,ix)  = lag_moment(ie,alpha,0.5*(3+beta-1))
          evec_e1(ie,ix)   = lag_moment(ie,alpha,0.5*(4+beta-1))
          evec_e105(ie,ix) = lag_moment(ie,alpha,0.5*(5+beta-1))
          evec_e2(ie,ix)   = lag_moment(ie,alpha,0.5*(6+beta-1))
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

  ! m(i,j) = int_0^inf x^q e^{-x^2} phi_i(x) psi_j(x) dx, with psi_j = phi_j,
  ! or psi_j = x d(phi_j)/dx if deriv
  subroutine lag_quad(ix,jx,q,deriv,m)
    use neo_globals, only : n_energy
    implicit none
    integer, intent(in) :: ix, jx, q
    logical, intent(in) :: deriv
    real, dimension(0:n_energy,0:n_energy), intent(out) :: m
    integer :: nq, ie, je
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

    do ie=0,n_energy
       do je=0,n_energy
          m(ie,je) = 0.5*sum(w*li(ie,:)*pj(je,:))
       enddo
    enddo

    deallocate(s,w,li,lj,pj)

  end subroutine lag_quad

  ! L_n^(alpha)(s), n=0..nmax, by the three-term recurrence
  subroutine lag_eval(nmax,alpha,ns,s,l)
    implicit none
    integer, intent(in) :: nmax, ns
    real, intent(in) :: alpha
    real, dimension(ns), intent(in) :: s
    real, dimension(0:nmax,ns), intent(out) :: l
    integer :: n

    l(0,:) = 1.0
    if (nmax > 0) l(1,:) = 1.0+alpha-s
    do n=1,nmax-1
       l(n+1,:) = ((2*n+1+alpha-s)*l(n,:)-(n+alpha)*l(n-1,:))/(n+1)
    enddo

  end subroutine lag_eval

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

  ! Composite Gauss-Legendre grid in x on [0,x_max], with weights times
  ! e^{-x^2}.  For lambda > 1 the panels double from 0.25/sqrt(lambda) to
  ! resolve the field-species layer at x ~ lambda^{-1/2}; then width 0.25.
  ! x_max covers e^{-x^2} x^{2 beta} L_n^2 for all basis functions.
  subroutine x_grid(lambda,x,w)
    use neo_globals, only : n_energy
    implicit none
    real, intent(in) :: lambda
    real, dimension(:), allocatable, intent(out) :: x, w
    integer, parameter :: ng=16
    real, parameter :: h=0.25
    real, dimension(ng) :: t, wt
    real, dimension(:), allocatable :: b
    real :: xmax, a
    integer :: np, k

    xmax = sqrt(6.0*n_energy+maxval(e_lag)+60.0)

    ! panel breakpoints
    allocate(b(0:1000))
    b(0) = 0.0
    np = 0
    if (lambda > 1.0) then
       a = h/sqrt(lambda)
       do while (a < h)
          np = np+1
          b(np) = a
          a = 2*a
       enddo
    endif
    do while (b(np) < xmax)
       np = np+1
       b(np) = b(np-1)+h
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

  ! Test-particle matrices emat_coll_test(is,js,:,:,ix), all ix.
  ! With phi_i = x^beta L_i^(alpha)(x^2), D phi_i = x d(phi_i)/dx and
  ! K_n as in kfun, the moment sums of the original formulation map as
  !   sum c c F(p+u,n)        -> int e^{-x^2} x^u phi_i phi_j K_n dx
  !   (2k+beta_i) c ...       -> D phi_i,  (2m+beta_j) c ... -> D phi_j
  !   sum c c Gamma((p+c)/2)  -> 2 int x^{c-1} e^{-x^2} phi_i phi_j dx
  ! model: 1 Connor (ccase 1: ma=mb, 2: ma<mb, 3: ma>mb), 2 HS0, 3 HS, 4/5 FP
  subroutine coll_test_ints(is,js,ir,tauinv_ab,lambda,model,ccase)
    use neo_globals, only : n_energy, n_xi, pi, temp, mass
    implicit none
    integer, intent(in) :: is, js, ir, model, ccase
    real, intent(in) :: tauinv_ab, lambda
    real, dimension(:), allocatable :: x, w, k0, k2, g
    real, dimension(:,:), allocatable :: l, phi, dphi, t
    real, dimension(:), allocatable :: r1
    real :: c0, ta, mr, alpha, beta, r2
    integer :: nq, ix, ie, q

    call x_grid(lambda,x,w)
    nq = size(x)
    allocate(k0(nq),k2(nq),g(nq),r1(0:n_energy))
    allocate(l(0:n_energy,nq),phi(0:n_energy,nq),dphi(0:n_energy,nq))
    allocate(t(0:n_energy,0:n_energy))
    do q=1,nq
       k0(q) = kfun(0,lambda*x(q)**2)
       k2(q) = kfun(2,lambda*x(q)**2)
    enddo

    c0 = tauinv_ab*sqrt(lambda/pi)
    ta = temp(is,ir)/temp(js,ir)
    mr = mass(js)/mass(is)

    do ix=0,n_xi

       alpha = 0.5*e_lag(ix)
       beta  = 1.0*xi_beta_l(ix)
       call lag_eval(n_energy,alpha,nq,x**2,l)
       do ie=0,n_energy
          phi(ie,:) = x**beta*l(ie,:)
          if (ie == 0) then
             dphi(ie,:) = x**beta*beta*l(0,:)
          else
             dphi(ie,:) = x**beta*((beta+2*ie)*l(ie,:)-2*(ie+alpha)*l(ie-1,:))
          endif
       enddo

       select case (model)

       case (1)

          select case (ccase)
          case (1)
             t = -c0*ix*(ix+1)*proj(phi,phi,k0-k2)
          case (2)
             ! 2 int x^{-1} e^{-x^2} phi_i phi_j; for beta=0 the p=0 term
             ! (divergent, omitted in the moment form) is removed
             t = -0.5*tauinv_ab*ix*(ix+1)*proj0(phi,1.0/x)
          case (3)
             t = -c0*ta*(2.0/3.0)*ix*(ix+1)*proj(phi,phi,x**2)
          end select

       case (2,3)

          t = -c0*ix*(ix+1)*proj(phi,phi,k0-k2)

          if (ix == 1) then
             ! u-restoring
             t = t+2*c0*proj(phi,phi,k0-k2-2*ta*(1+mr)*x**2*k2)
             if (model == 3) then
                ! k-heating friction
                g = 12*x**3*k2-8*ta*(1+mr)*x**5*k2+4*(x**3*k0-x**3*k2)
                r2 = sum(w*x**3*g)
                do ie=0,n_energy
                   r1(ie) = sum(w*phi(ie,:)*g)
                enddo
                do ie=0,n_energy
                   t(ie,:) = t(ie,:)+c0*r1(ie)*r1(:)/r2
                enddo
             endif
          endif

          if (model == 3 .and. ix == 2) then
             ! pi-energy restoring
             t = t+c0*proj(phi,phi,-8*ta*(1+mr)*x**2*k2+4*k0)
          endif

          if (model == 3 .and. ix == 0) then
             ! energy diffusion
             t = t-2*c0*proj(dphi,dphi,k2)
          endif

       case default

          ! full linearized FP: pitch-angle + energy scattering
          t = -c0*ix*(ix+1)*proj(phi,phi,k0-k2) &
               +4*c0*((1-ta)*proj(dphi,phi,x**2*k2)-0.5*proj(dphi,dphi,k2))

       end select

       emat_coll_test(is,js,:,:,ix) = t

    enddo

    deallocate(x,w,k0,k2,g,r1,l,phi,dphi,t)

  contains

    ! P(i,j) = sum_q w_q f_q a(i,q) b(j,q)
    function proj(a,b,f)
      real, dimension(0:n_energy,nq), intent(in) :: a, b
      real, dimension(nq), intent(in) :: f
      real, dimension(0:n_energy,0:n_energy) :: proj
      integer :: i, j
      do i=0,n_energy
         do j=0,n_energy
            proj(i,j) = sum(w*f*a(i,:)*b(j,:))
         enddo
      enddo
    end function proj

    ! As proj(a,a,f), but for beta=0 with the constant term a_i(0) a_j(0)
    ! removed (the divergent p=0 moment, omitted in the moment form);
    ! L_n^(alpha)(0) = (alpha+1)_n/n!
    function proj0(a,f)
      real, dimension(0:n_energy,nq), intent(in) :: a
      real, dimension(nq), intent(in) :: f
      real, dimension(0:n_energy,0:n_energy) :: proj0
      real :: a0(0:n_energy)
      integer :: i, j
      a0(:) = 0.0
      if (xi_beta_l(ix) == 0) then
         a0(0) = 1.0
         do i=1,n_energy
            a0(i) = a0(i-1)*(alpha+i)/i
         enddo
      endif
      do i=0,n_energy
         do j=0,n_energy
            proj0(i,j) = sum(w*f*(a(i,:)*a(j,:)-a0(i)*a0(j)))
         enddo
      enddo
    end function proj0

  end subroutine coll_test_ints

  subroutine ENERGY_coll_ints_alloc(flag)
    use neo_globals, only : n_species, n_energy,n_xi
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

  subroutine ENERGY_coll_ints(ir)
    ! Collision matrices in the energy basis, by quadrature
    ! (coll_test_ints, coll_field_ints); no Laguerre power-series sums.
    use neo_globals
    implicit none
    integer, intent(in) :: ir
    real :: tauinv_ab, tauinv_ba, lambda, rs, ru
    integer :: is, js, ie, je, ix
    integer :: model, ccase

    emat_coll_test(:,:,:,:,:)  = 0.0
    emat_coll_field(:,:,:,:,:) = 0.0

    do is=1, n_species
       do js=1, n_species

          ! (Note: pol part of dens from rotation will be added 
          !  in coll term in kinetic equation)
          tauinv_ab = nu(is,ir) * (1.0*Z(js))**2 / (1.0*Z(is))**2 & 
               * dens(js,ir)/dens(is,ir)
          tauinv_ba = nu(js,ir) * (1.0*Z(is))**2 / (1.0*Z(js))**2 & 
               * dens(is,ir)/dens(js,ir)

          lambda = (vth(is,ir) / vth(js,ir))**2

          ccase = 0
          if(collision_model == 1 .or. &
               (spitzer_model==1 .and. js .ne. is)) then
             ! Connor model
             model = 1
             if(is == js .or. &
                  (abs(mass(is)-mass(js)) < epsilon(0.))) then
                ccase = 1
             else if(mass(is) < mass(js)) then
                ccase = 2
             else
                ccase = 3
             endif
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
             do ie=0, n_energy
                do je=0,n_energy

                   ix = 1
                   if(abs(emat_coll_test(is,js,(1-xi_beta_l(ix))/e_alpha,&
                        (1-xi_beta_l(ix))/e_alpha,ix)) &
                        > epsilon(0.)) then
                      rs = -(mass(js)*dens(js,ir)*vth(js,ir)) &
                           / (mass(is)*dens(is,ir)*vth(is,ir)) &
                           * emat_coll_test(is,js,ie,&
                           (1-xi_beta_l(ix))/e_alpha,ix) &
                           * emat_coll_test(js,is,&
                           (1-xi_beta_l(ix))/e_alpha,je,ix) &
                           / emat_coll_test(is,js,(1-xi_beta_l(ix))/e_alpha,&
                           (1-xi_beta_l(ix))/e_alpha,ix)
                   else
                      rs = 0.0
                   endif
                   emat_coll_field(is,js,ie,je,ix) = rs

                   ix = 0
                   if(abs(emat_coll_test(is,js,(2-xi_beta_l(ix))/e_alpha,&
                        (2-xi_beta_l(ix))/e_alpha,ix)) &
                        > epsilon(0.)) then
                      ru = -(dens(js,ir)*temp(js,ir)) &
                           / (dens(is,ir)*temp(is,ir)) &
                           * emat_coll_test(is,js,ie,&
                           (2-xi_beta_l(ix))/e_alpha,ix) &
                           * emat_coll_test(js,is,&
                           (2-xi_beta_l(ix))/e_alpha,je,ix) &
                           / emat_coll_test(is,js,(2-xi_beta_l(ix))/e_alpha,&
                           (2-xi_beta_l(ix))/e_alpha,ix)
                   else
                      ru = 0.0
                   endif
                   emat_coll_field(is,js,ie,je,ix) = ru

                enddo
             enddo
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

             tauinv_ab = nu(is,ir) * (1.0*Z(js))**2 / (1.0*Z(is))**2 & 
                  * dens(js,ir)/dens(is,ir)
             tauinv_ba = nu(js,ir) * (1.0*Z(is))**2 / (1.0*Z(js))**2 & 
                  * dens(is,ir)/dens(js,ir)
             lambda = (vth(is,ir) / vth(js,ir))**2

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

  ! Merged composite Gauss-Legendre grid on [0,inf): panels of width
  ! h*s, s = 1, lambda^{-1/2}, lambda^{1/2}, each sequence to x_max*s, so
  ! that functions of x, sqrt(lambda) x and x/sqrt(lambda) are all resolved.
  ! Between the smallest and largest scale the panels also double in width,
  ! for the power-law tails x^{-k} of the Rosenbluth potentials there.
  ! Raw weights (no Maxwellian).
  subroutine field_grid(lambda,x,w)
    use neo_globals, only : n_energy
    implicit none
    real, intent(in) :: lambda
    real, dimension(:), allocatable, intent(out) :: x, w
    integer, parameter :: ng=16
    real, parameter :: h=0.25
    real, dimension(ng) :: t, wt
    real, dimension(:), allocatable :: b
    real :: xmax, s(3), v
    integer :: kmax, i, k, np, nb

    xmax = sqrt(6.0*n_energy+maxval(e_lag)+60.0)
    kmax = ceiling(xmax/h)
    s = (/ 1.0, 1.0/sqrt(lambda), sqrt(lambda) /)

    allocate(b(3*(kmax+1)+200))
    nb = 0
    do i=1,3
       do k=0,kmax
          nb = nb+1
          b(nb) = h*s(i)*k
       enddo
    enddo
    if (lambda /= 1.0) then
       v = h*minval(s)
       do while (v < h*maxval(s))
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

  ! Gauss-Jacobi rule on [0,1] for weight t^p.
  ! Jacobi weight (1-z)^0 (1+z)^p on [-1,1], mapped by t=(1+z)/2.
  subroutine gauss_jacobi01(n,p,t,w)
    implicit none
    integer, intent(in) :: n, p
    real, dimension(n), intent(out) :: t, w
    real, dimension(n) :: e
    real :: b, c
    integer :: k

    b = 1.0*p
    do k=0,n-1
       c = 2*k+b
       if (k == 0) then
          t(1) = b/(b+2)
       else
          t(k+1) = b*b/(c*(c+2))
       endif
    enddo
    e(:) = 0.0
    do k=1,n-1
       c = 2*k+b
       e(k) = sqrt(4.0*k*k*(k+b)*(k+b)/(c*c*(c+1)*(c-1)))
    enddo
    call gauss_jacobi_rule(n,t,e,2.0**(p+1)/(p+1),w)
    t = 0.5*(t+1.0)
    w = w/2.0**(p+1)

  end subroutine gauss_jacobi01

  ! it(n,q) = int_0^1 t^p e^{-u} L_n^(alpha)(u) dt, u = c x_q^2 t^2,
  ! for x_q with mask(q).  For c x^2 beyond the basis extent the tail of
  ! int_0^inf is negligible and the closed form is used:
  !   int_0^inf y^p e^{-c y^2} L_n(c y^2) dy = c^{-(p+1)/2} lag_moment(n,alpha,(p-1)/2)
  subroutine itilde(nmax,alpha,p,c,nq,x,mask,it)
    implicit none
    integer, intent(in) :: nmax, p, nq
    real, intent(in) :: alpha, c
    real, dimension(nq), intent(in) :: x
    logical, dimension(nq), intent(in) :: mask
    real, dimension(0:nmax,nq), intent(out) :: it
    real, dimension(:), allocatable :: t, w, u
    real, dimension(:,:), allocatable :: l
    real, dimension(0:nmax) :: full
    real :: scut, s
    integer :: nt, n, q

    scut = 10.0*nmax+2*p+80
    nt = nmax+ceiling(3*sqrt(scut))+12
    allocate(t(nt),w(nt),u(nt),l(0:nmax,nt))
    call gauss_jacobi01(nt,p,t,w)
    do n=0,nmax
       full(n) = c**(-0.5*(p+1))*lag_moment(n,alpha,0.5*(p-1))
    enddo

    it(:,:) = 0.0
    do q=1,nq
       if (.not. mask(q)) cycle
       s = c*x(q)**2
       if (s > scut) then
          it(:,q) = full(:)/x(q)**(p+1)
       else
          u = s*t**2
          call lag_eval(nmax,alpha,nt,u,l)
          do n=0,nmax
             it(n,q) = sum(w*exp(-u)*l(n,:))
          enddo
       endif
    enddo

    deallocate(t,w,u,l)

  end subroutine itilde

  ! Field-particle matrices emat_coll_field(is,js,:,:,ix), all ix.
  ! x is the test-species (a) speed, y the field-species (b) speed, both in
  ! units of v_ta; phi_j(sqrt(lambda) y) is the b basis in its own units.
  ! Separable terms (models 1-3) are single integrals with kernels
  ! K_n(x;lambda) (k0,k2) and K_n(x;1/lambda) (q0,q2) as in kfun.
  ! Full FP (model 4):
  !   C = tau_ab 2/sqrt(pi) lambda^{3/2} [ m_a/m_b R + sum_k w_k D_k ]
  !   R      = 2 int x^2 e^{-(1+lambda)x^2} phi_i(x) phi_j(sqrt(lambda) x) dx
  !   D(A,B) = int dx x^A e^{-x^2} phi_i(x) int_0^x dy y^B e^{-lambda y^2} phi_j(sqrt(lambda) y)
  !   Db(A,B)= int dy y^B e^{-lambda y^2} phi_j(sqrt(lambda) y) int_0^y dx x^A e^{-x^2} phi_i(x)
  ! The inner integral always carries the positive power (l+2 or l+4), in
  ! the scaled form itilde; the combined outer power is A+B+2 beta+1 >= 4.
  subroutine coll_field_ints(is,js,ir,tauinv_ab,tauinv_ba,lambda,model,ccase)
    use neo_globals, only : n_energy, n_xi, pi, temp, mass, dens, vth
    implicit none
    integer, intent(in) :: is, js, ir, model, ccase
    real, intent(in) :: tauinv_ab, tauinv_ba, lambda
    real, dimension(:), allocatable :: x, w, wa, wb, k0, k2, q0, q2, g1, g3
    real, dimension(:,:), allocatable :: la, lb, phi, dphi, f
    real, dimension(:,:), allocatable :: it1, it2, ib1, ib2
    real, dimension(0:n_energy) :: va, vb
    logical, dimension(:), allocatable :: ma, mb
    real :: xmax, alpha, beta, ta, tb, mab, mm, r2, r3, cf
    real :: w2, w3, w45, w67
    integer :: nq, ix, ie, q, l

    if (model == 5) return

    call field_grid(lambda,x,w)
    nq = size(x)
    xmax = sqrt(6.0*n_energy+maxval(e_lag)+60.0)

    allocate(ma(nq),mb(nq),wa(nq),wb(nq))
    ma = x <= xmax
    mb = sqrt(lambda)*x <= xmax
    wa = 0.0
    wb = 0.0
    where (ma) wa = w*exp(-x*x)
    where (mb) wb = w*exp(-lambda*x*x)

    allocate(k0(nq),k2(nq),q0(nq),q2(nq),g1(nq),g3(nq))
    k0 = 0.0; k2 = 0.0; q0 = 0.0; q2 = 0.0
    do q=1,nq
       if (.not. ma(q)) cycle
       k0(q) = kfun(0,lambda*x(q)**2)
       k2(q) = kfun(2,lambda*x(q)**2)
       q0(q) = kfun(0,x(q)**2/lambda)
       q2(q) = kfun(2,x(q)**2/lambda)
    enddo

    allocate(la(0:n_energy,nq),lb(0:n_energy,nq))
    allocate(phi(0:n_energy,nq),dphi(0:n_energy,nq),f(0:n_energy,0:n_energy))
    allocate(it1(0:n_energy,nq),it2(0:n_energy,nq),ib1(0:n_energy,nq),ib2(0:n_energy,nq))

    ta  = temp(is,ir)/temp(js,ir)
    tb  = 1.0/ta
    mab = mass(is)/mass(js)
    mm  = (mass(js)*dens(js,ir)*vth(js,ir))/(mass(is)*dens(is,ir)*vth(is,ir))

    do ix=0,n_xi

       l = ix
       alpha = 0.5*e_lag(ix)
       beta  = 1.0*xi_beta_l(ix)

       ! test-species basis on x <= xmax; zero elsewhere
       call lag_eval(n_energy,alpha,nq,merge(x*x,0.0,ma),la)
       do ie=0,n_energy
          phi(ie,:) = merge(x**beta*la(ie,:),0.0,ma)
          if (ie == 0) then
             dphi(ie,:) = beta*phi(0,:)
          else
             dphi(ie,:) = merge(x**beta*((beta+2*ie)*la(ie,:) &
                  -2*(ie+alpha)*la(ie-1,:)),0.0,ma)
          endif
       enddo

       f(:,:) = 0.0

       select case (model)

       case (1)

          if (l == 1) then
             select case (ccase)
             case (1)
                va = vec(phi,x*(k0-k2))
                vb = vec(phi,x*(q0-q2))
                f = mm*2*tauinv_ba/sqrt(lambda*pi)*outer(va,vb)/sum(wa*x**2*(k0-k2))
             case (2)
                do ie=0,n_energy
                   va(ie) = 2*lag_moment(ie,alpha,0.5*(beta-1))
                   vb(ie) = 2*lag_moment(ie,alpha,0.5*(beta+2))
                enddo
                f = mm*(2.0/3.0)*tauinv_ba*tb/sqrt(lambda*pi)*outer(va,vb)
             case (3)
                do ie=0,n_energy
                   va(ie) = 2*lag_moment(ie,alpha,0.5*(beta+2))
                   vb(ie) = 2*lag_moment(ie,alpha,0.5*(beta-1))
                enddo
                f = mm*0.5*tauinv_ba*tb/gamma(2.5)*outer(va,vb)
             end select
          endif

       case (2,3)

          if (l == 1) then
             ! slowing-down
             va = vec(phi,x**3*k2)
             vb = vec(phi,x**3*q2)
             f = mm*4*tauinv_ba*tb*(1+mab)/sqrt(lambda*pi) &
                  *outer(va,vb)/sum(wa*x**4*k2)
             if (model == 3) then
                ! h-heating friction
                g1 = 4*(2*mab-1)*x**3*k2+8*ta*(1+1/mab)*x**5*k2 &
                     -4*(1+mab)*(x**3*k0-x**3*k2)
                g3 = 4*(2/mab-1)*x**3*q2+8*tb*(1+mab)*x**5*q2 &
                     -4*(1+1/mab)*(x**3*q0-x**3*q2)
                r2 = sum(wa*x**3*g1)
                f = f+(mass(js)*dens(js,ir)*vth(js,ir)**3) &
                     /(mass(is)*dens(is,ir)*vth(is,ir)**3) &
                     *1.5*tauinv_ba/sqrt(lambda*pi) &
                     *outer(vec(phi,g1),vec(phi,g3))/r2
             endif
          endif

          ! nup-energy restoring.  For odd beta the moment form used
          ! F(m,n) with m+n even, which neo_compute_fcoll leaves zero, so
          ! the term vanished; that behaviour is kept.
          if (model == 3 .and. l == 2 .and. mod(xi_beta_l(ix),2) == 0) then
             g1 = 8*(1+1.5*mab)*x**2*k2+8*ta*(1+1/mab)*x**4*k2 &
                  -4*(1+1.5*mab)*(x**2*k0-x**2*k2)
             g3 = 8*(1+1.5/mab)*x**2*q2+8*tb*(1+mab)*x**4*q2 &
                  -4*(1+1.5/mab)*(x**2*q0-x**2*q2)
             r2 = sum(wa*x**2*g1)
             f = (temp(js,ir)*dens(js,ir))/(temp(is,ir)*dens(is,ir)) &
                  *tauinv_ba/sqrt(lambda*pi) &
                  *outer(vec(phi,g1),vec(phi,g3))/r2
          endif

          if (model == 3 .and. l == 0) then
             ! energy diffusion (restoring)
             r3 = sum(wa*x**3*(8*ta*(1+1/mab)*x**3*k2-4*x*k0))
             va = vec(dphi,x**2*k2)
             vb = vec(phi,8*tb*(1+mab)*x**4*q2-4*x**2*q0)
             f = 4*sqrt(lambda/pi) &
                  *(temp(js,ir)*dens(js,ir))/(temp(is,ir)*dens(is,ir)) &
                  *tauinv_ba/lambda*outer(va,vb)/r3
          endif

       case (4)

          ! field-species basis on sqrt(lambda) x <= xmax
          call lag_eval(n_energy,alpha,nq,merge(lambda*x*x,0.0,mb),lb)
          do ie=0,n_energy
             lb(ie,:) = merge(lb(ie,:),0.0,mb)
          enddo
          cf = lambda**(0.5*beta)

          ! R: local term
          f = mab*2*cf*proj2(la,lb,merge(w*exp(-(1+lambda)*x*x),0.0,ma.and.mb) &
               *x**(2*beta+2))

          ! inner integrals with positive power p = l+2+beta, l+4+beta
          call itilde(n_energy,alpha,l+2+nint(beta),lambda,nq,x,ma,it1)
          call itilde(n_energy,alpha,l+4+nint(beta),lambda,nq,x,ma,it2)
          call itilde(n_energy,alpha,l+2+nint(beta),1.0,nq,x,mb,ib1)
          call itilde(n_energy,alpha,l+4+nint(beta),1.0,nq,x,mb,ib2)

          w2  = -2.0/(l+0.5)*(mab-l*(1-mab))
          w3  = -2.0/(l+0.5)*(1+l*(1-mab))
          w45 = -l*(l-1.0)/(l*l-0.25)
          w67 = (l+1.0)*(l+2.0)/((l+1.5)*(l+0.5))

          ! D(A,B), outer x:  x^{A+B+2beta+1} e^{-x^2} L_i(x^2) cf it_j
          ! Db(A,B), outer y: y^{A+B+2beta+1} e^{-lambda y^2} cf L_j(lambda y^2) ib_i
          ! A+B = 3 (w2,w3) or 5 (w45,w67)
          f = f+w2*cf*proj2(la,it1,wa*x**(4+2*beta)) &
               +w3*cf*transpose(proj2(lb,ib1,wb*x**(4+2*beta))) &
               +w67*cf*proj2(la,it2,wa*x**(6+2*beta)) &
               +w67*cf*transpose(proj2(lb,ib2,wb*x**(6+2*beta)))
          if (l > 1) then
             f = f+w45*cf*proj2(la,it1,wa*x**(6+2*beta)) &
                  +w45*cf*transpose(proj2(lb,ib1,wb*x**(6+2*beta)))
          endif

          f = tauinv_ab*2/sqrt(pi)*lambda**1.5*f

       end select

       emat_coll_field(is,js,:,:,ix) = f

    enddo

    deallocate(x,w,wa,wb,ma,mb,k0,k2,q0,q2,g1,g3,la,lb,phi,dphi,f,it1,it2,ib1,ib2)

  contains

    ! v(i) = sum_q wa_q g_q a(i,q)
    function vec(a,g)
      real, dimension(0:n_energy,nq), intent(in) :: a
      real, dimension(nq), intent(in) :: g
      real, dimension(0:n_energy) :: vec
      integer :: i
      do i=0,n_energy
         vec(i) = sum(wa*g*a(i,:))
      enddo
    end function vec

    ! P(i,j) = sum_q v_q a(i,q) b(j,q)
    function proj2(a,b,v)
      real, dimension(0:n_energy,nq), intent(in) :: a, b
      real, dimension(nq), intent(in) :: v
      real, dimension(0:n_energy,0:n_energy) :: proj2
      integer :: i, j
      do i=0,n_energy
         do j=0,n_energy
            proj2(i,j) = sum(v*a(i,:)*b(j,:))
         enddo
      enddo
    end function proj2

    function outer(a,b)
      real, dimension(0:n_energy), intent(in) :: a, b
      real, dimension(0:n_energy,0:n_energy) :: outer
      integer :: i
      do i=0,n_energy
         outer(i,:) = a(i)*b(:)
      enddo
    end function outer

  end subroutine coll_field_ints

  subroutine write_fullcoll_mono(ir)
    use neo_globals
    implicit none
    integer, intent(in) :: ir
    real :: tauinv_ab, tauinv_ba, lambda
    integer :: xarg, yarg
    real :: rd, r1, r2, r3, r4, r5, r6, r7
    real, dimension(:,:), allocatable :: fcoll, fcoll_bar, &
         fcollinv, fcollinv_bar
    integer :: fmarg
    integer :: is, js, ie, je, ix, jx
    real, dimension(:,:,:,:,:), allocatable :: test_mono
    real, dimension(:,:,:,:,:), allocatable :: field_mono
    real, dimension(:,:), allocatable :: ctot

    allocate(test_mono(n_species,n_species,0:n_energy,0:n_energy,0:n_xi))
    allocate(field_mono(n_species,n_species,0:n_energy,0:n_energy,0:n_xi))

    test_mono(:,:,:,:,:)  = 0.0
    field_mono(:,:,:,:,:) = 0.0

    fmarg = 2*e_alpha*n_energy + (2*n_xi)+4
    allocate(fcoll(-fmarg:fmarg,-fmarg:fmarg))
    allocate(fcoll_bar(-fmarg:fmarg,-fmarg:fmarg))
    allocate(fcollinv(-fmarg:fmarg,-fmarg:fmarg))
    allocate(fcollinv_bar(-fmarg:fmarg,-fmarg:fmarg))
    allocate(ctot(0:n_energy,0:n_energy))

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

          ! (Note: pol part of dens from rotation will be added 
          !  in coll term in kinetic equation)
          tauinv_ab = nu(is,ir) * (1.0*Z(js))**2 / (1.0*Z(is))**2 & 
               * dens(js,ir)/dens(is,ir)
          tauinv_ba = nu(js,ir) * (1.0*Z(is))**2 / (1.0*Z(js))**2 & 
               * dens(is,ir)/dens(js,ir)

          lambda = (vth(is,ir) / vth(js,ir))**2
          call neo_compute_fcoll(fmarg,lambda,fcoll,fcoll_bar)
          call neo_compute_fcoll(fmarg,1.0/lambda,fcollinv,fcollinv_bar)

          do ie=0, n_energy
             do je=0,n_energy
                do ix=0, n_xi

                   jx = ix

                   ! Full linearized FP op
                   xarg = e_alpha*(ie+je) &
                        + xi_beta_l(ix) + xi_beta_l(jx)
                   test_mono(is,js,ie,je,ix) = &
                        test_mono(is,js,ie,je,ix) &
                        -tauinv_ab &
                        * sqrt(lambda/pi) *ix*(ix+1) &
                        * (fcoll(xarg-1,0) - fcoll(xarg-3,2)) 

                   rd = 4.0*tauinv_ab * sqrt(lambda/pi) &
                        * (e_alpha*ie + xi_beta_l(ix)) &
                        * ( (1.0 - temp(is,ir)/temp(js,ir)) &
                        * fcoll(xarg-1,2) &
                        - 0.5 * (e_alpha*je + xi_beta_l(jx)) &
                        * fcoll(xarg-3,2) )

                   test_mono(is,js,ie,je,ix) = &
                        test_mono(is,js,ie,je,ix) + rd

                   ! Real full field-particle operator

                   r1 = mass(is)/mass(js) &
                        * (1.0 + lambda)**(-0.5*(xarg+3)) &
                        * mygamma2(xarg + 3)

                   xarg = e_alpha*ie + xi_beta_l(ix)
                   yarg = e_alpha*je + xi_beta_l(jx)

                   r2 = -2.0/(ix+0.5) &
                        * (mass(is)/mass(js) &
                        - ix*(1.0-mass(is)/mass(js))) &
                        * fcoll(xarg-ix+1,yarg+jx+2)

                   r3 = -2.0/(ix+0.5) &
                        * (1.0 + ix*(1.0-mass(is)/mass(js))) &
                        * fcoll_bar(xarg+ix+2,yarg-jx+1)

                   r4 = -ix*(ix-1.0)/(ix*ix - 0.25) &
                        * fcoll(xarg-ix+3,yarg+jx+2)

                   r5 = -ix*(ix-1.0)/(ix*ix - 0.25) &
                        * fcoll_bar(xarg+ix+2,yarg-jx+3)

                   r6 = (ix+1.0)*(ix+2.0)/(ix+1.5)/(ix+0.5) &
                        * fcoll(xarg-ix+1,yarg+jx+4)

                   r7 = (ix+1.0)*(ix+2.0)/(ix+1.5)/(ix+0.5) &
                        * fcoll_bar(xarg+ix+4,yarg-jx+1)

                   field_mono(is,js,ie,je,ix) = &
                        field_mono(is,js,ie,je,ix) &
                        + tauinv_ab * 2.0/sqrt(pi) * lambda**1.5 &
                        * lambda**(0.5*(e_alpha*je &
                        + xi_beta_l(ix))) &
                        * (r1 + r2 + r3 + r4 + r5 + r6 + r7)

                enddo
             enddo
          enddo

          open(unit=1,file=trim(path)//'out.neo.diagnostic_coll',status='old',position='append')
          do ix=0,n_xi
             ctot(:,:) = (test_mono(is,js,:,:,ix)+field_mono(is,js,:,:,ix))/tauinv_ab 
             write(1,'(a,i2,a,i2,a,i2,a)') '# ix = ',ix,'  (is,js)=(',is,',',js,')'
             do ie=0,n_energy
                do je=ie,n_energy
                   write(1,'(1pe16.8,1x)',advance='no') ctot(ie,je)
                enddo
             enddo
             write(1,*)
          enddo
          open(unit=1,file=trim(path)//'out.neo.diagnostic_colltest',status='old',position='append')
          do ix=0,n_xi
             ctot(:,:) = (test_mono(is,js,:,:,ix))/tauinv_ab 
             write(1,'(a,i2,a,i2,a,i2,a)') '# ix = ',ix,'  (is,js)=(',is,',',js,')'
             do ie=0,n_energy
                do je=0,n_energy
                   write(1,'(1pe16.8,1x)',advance='no') ctot(ie,je)
                enddo
             enddo
             write(1,*)
          enddo
          open(unit=1,file=trim(path)//'out.neo.diagnostic_collfield',status='old',position='append')
          do ix=0,n_xi
             ctot(:,:) = (field_mono(is,js,:,:,ix))/tauinv_ab
             write(1,'(a,i2,a,i2,a,i2,a)') '# ix = ',ix,'  (is,js)=(',is,',',js,')'
             do ie=0,n_energy
                do je=0,n_energy
                   write(1,'(1pe16.8,1x)',advance='no') ctot(ie,je)
                enddo
             enddo
             write(1,*)
          enddo

       enddo ! js
    enddo ! is

    close(1)

  end subroutine write_fullcoll_mono

  ! returns Gamma(n/2)
  real function gamma2(n)
    use neo_globals, only : pi
    implicit none
    integer, intent (in) :: n
    integer :: i, i1
    
    if(mod(n,2) == 0) then
       ! Gamma of integer
       gamma2 = 1.0
       i1=4
    else 
       ! Gamma of half integer
       gamma2 = sqrt(pi)
       i1=3
    endif

    do i=i1,n,2
       gamma2 = gamma2 * (i/2.0-1.0)
    enddo

  end function gamma2

end module neo_energy_grid
