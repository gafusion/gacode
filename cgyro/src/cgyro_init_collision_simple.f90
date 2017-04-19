subroutine cgyro_init_collision_simple

  use timer_lib

  use cgyro_globals

  implicit none

  real, dimension(:,:,:), allocatable :: nu_d

  real :: xa, xb, tauinv_ab
  integer :: is,it,ix,ie,js,jx
  ! parameters for matrix solve
  real, dimension(:,:), allocatable :: amat
  real, dimension(:,:,:,:), allocatable :: ctest

  allocate(nu_d(n_energy,n_species,n_species))
  nu_d(:,:,:) = 0.0

  do ie=1,n_energy
     do is=1,n_species
        do js=1,n_species

           xa = vel(ie)
           xb = xa * vth(is) / vth(js)
           tauinv_ab = nu(is) * z(js)**2 / z(is)**2 &
                * dens(js)/dens(is)

           ! Only ee,ei Connor-like Lorentz
           if (is == is_ele) then
              if (is == js) then
                 ! e-e
                 nu_d(ie,is,js) = tauinv_ab * (1.0/xa**3) &
                      * (exp(-xb*xb)/(xb*sqrt(pi)) &
                      + (1.0-1.0/(2.0*xb*xb)) * erf(xb))
              else
                 ! e-i
                 nu_d(ie,is,js) = tauinv_ab * (1.0/xa**3)
              endif
           endif

        enddo
     enddo
  enddo

  allocate(ctest(n_species,n_energy,n_xi,n_xi))
  allocate(amat(n_xi,n_xi))

  ! Collision test particle component
  ctest = 0.0

  ! Lorentz
  do is=1,n_species
     do ie=1,n_energy
        do ix=1,n_xi
           do jx=1, n_xi
              do js=1,n_species
                 ctest(is,ie,ix,jx) &
                      = ctest(is,ie,ix,jx) &
                      + xi_lor_mat(ix,jx) *0.5*nu_d(ie,is,js)
              enddo
           enddo
        enddo
     enddo
  enddo

  ! matrix solve parameters
  allocate(i_piv(n_xi))

  ! set-up the collision matrix

  do it=1,n_theta
     do is=1,n_species
        do ie=1,n_energy

           cmat_simple(:,:,ie,is,it) = 0.0
           amat(:,:) = 0.0

                 ! constant part
           do ix=1,n_xi
              cmat_simple(ix,ix,ie,is,it) = 1.0
              amat(ix,ix) = 1.0
           enddo

           do jx=1,n_xi
              do ix=1,n_xi

                 ! Collision component: Test particle
                 cmat_simple(ix,jx,ie,is,it) = cmat_simple(ix,jx,ie,is,it) &
                      - (0.5*delta_t) * ctest(is,ie,ix,jx)
                 amat(ix,jx) = amat(ix,jx) &
                      + (0.5*delta_t) * ctest(is,ie,ix,jx)

                 ! Trapping 
                 cmat_simple(ix,jx,ie,is,it) = cmat_simple(ix,jx,ie,is,it) &
                      + (0.5*delta_t) * (omega_trap(it,is) * vel(ie) &
                      + omega_rot_trap(it,is) / vel(ie)) &
                      * (1.0 - xi(ix)**2) * xi_deriv_mat(ix,jx) 
                       
                 amat(ix,jx) = amat(ix,jx) &
                      - (0.5*delta_t) * (omega_trap(it,is) * vel(ie) &
                      + omega_rot_trap(it,is) / vel(ie)) &
                      * (1.0 - xi(ix)**2) * xi_deriv_mat(ix,jx)
                      

              enddo
           enddo

           ! H_bar = (1 - dt/2 C)^(-1) * (1 + dt/2 C) H
           ! Lapack factorization and inverse of LHS
           call DGESV(n_xi,n_xi,cmat_simple(:,:,ie,is,it),n_xi,&
                i_piv,amat,n_xi,info)
           cmat_simple(:,:,ie,is,it) = amat(:,:)

        enddo
     enddo
  enddo

  deallocate(amat)
  deallocate(i_piv)
  deallocate(nu_d)
  deallocate(ctest)

end subroutine cgyro_init_collision_simple
