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
    use neo_globals, only : n_energy, n_xi
    implicit none
    integer :: xarg
    integer :: ie, je, ke, me, ix, jx, kx
    real :: zarg0, zarg1, zarg2

    ! vectors
    do ie=0, n_energy
       do ix=0, n_xi

          evec_e0(ie,ix)   = 0.0
          evec_e1(ie,ix)   = 0.0
          evec_e2(ie,ix)   = 0.0
          evec_e05(ie,ix)  = 0.0
          evec_e105(ie,ix) = 0.0

          do ke=0, ie

             zarg0 = (-1.0)**ke
             zarg1 = mygamma2(2 + 2*ie + e_lag(ix)) &
                  / mygamma2(2 + 2*(ie-ke)) &
                  / mygamma2(2 + 2*ke + e_lag(ix))
             zarg2 = mygamma2(2 + 2*ke)

             evec_e0(ie,ix) = evec_e0(ie,ix) +  0.5 * zarg0 * zarg1 &
                  * (mygamma2(e_alpha*ke + xi_beta_l(ix) + 3) / zarg2)

             evec_e1(ie,ix) = evec_e1(ie,ix) +  0.5 * zarg0 * zarg1 &
                  * (mygamma2(e_alpha*ke + xi_beta_l(ix) + 5) / zarg2)

             evec_e2(ie,ix) = evec_e2(ie,ix) +  0.5 * zarg0 * zarg1 &
                  * (mygamma2(e_alpha*ke + xi_beta_l(ix) + 7) / zarg2)
             
             evec_e05(ie,ix) = evec_e05(ie,ix) + 0.5 * zarg0 * zarg1 &
                  * (mygamma2(e_alpha*ke + xi_beta_l(ix) + 4) / zarg2)

             evec_e105(ie,ix) = evec_e105(ie,ix) + 0.5 * zarg0 * zarg1 &
                  * (mygamma2(e_alpha*ke + xi_beta_l(ix) + 6) / zarg2)
             
          enddo
       enddo
    enddo

    do ie=0, n_energy
       do je=0,n_energy
          do ix=0, n_xi
 
             emat_e05(ie,je,ix,:)   = 0.0
             emat_en05(ie,je,ix,:)  = 0.0
             emat_e05de(ie,je,ix,:) = 0.0
             emat_e0(ie,je,ix,:)    = 0.0
             emat_e1(ie,je,ix,:)    = 0.0

             do ke=0,ie
                do me=0,je
                   
                   kx=1
                   jx=ix
                   zarg0 = (-1.0)**(ke+me)
                   zarg1 = (mygamma2(2 + 2*ie + e_lag(ix)) &
                        / mygamma2(2 + 2*(ie-ke)) &
                        / mygamma2(2 + 2*ke + e_lag(ix))) &
                        * (mygamma2(2 + 2*je + e_lag(jx)) &
                        / mygamma2(2 + 2*(je-me)) &
                        / mygamma2(2 + 2*me + e_lag(jx)))
                   zarg2 = mygamma2(2 + 2*ke) * mygamma2(2 + 2*me)
                   xarg = e_alpha*(ke+me) + xi_beta_l(ix) &
                        + xi_beta_l(jx) + 3
                   emat_e0(ie,je,ix,kx) = emat_e0(ie,je,ix,kx) &
                        + 0.5*zarg0 * zarg1 * (mygamma2(xarg) / zarg2)

                   do kx=1,2
                      if(kx == 1) then 
                         jx = ix-1
                      else
                         jx = ix+1
                      endif

                      if (jx >= 0 .and. jx <= n_xi) then

                         zarg0 = (-1.0)**(ke+me)
                         zarg1 = (mygamma2(2 + 2*ie + e_lag(ix)) &
                              / mygamma2(2 + 2*(ie-ke)) &
                              / mygamma2(2 + 2*ke + e_lag(ix))) &
                              * (mygamma2(2 + 2*je + e_lag(jx)) &
                              / mygamma2(2 + 2*(je-me)) &
                              / mygamma2(2 + 2*me + e_lag(jx)))
                         zarg2 = mygamma2(2 + 2*ke) * mygamma2(2 + 2*me)

                         xarg = e_alpha*(ke+me) + xi_beta_l(ix) &
                              + xi_beta_l(jx) + 4
                         emat_e05(ie,je,ix,kx) = emat_e05(ie,je,ix,kx) &
                              + 0.5*zarg0 * zarg1 * (mygamma2(xarg) / zarg2)

                         xarg = e_alpha*(ke+me) + xi_beta_l(ix) &
                              + xi_beta_l(jx) + 2
                         emat_en05(ie,je,ix,kx) = emat_en05(ie,je,ix,kx) &
                              + 0.5*zarg0 * zarg1 * (mygamma2(xarg) / zarg2)
                         emat_e05de(ie,je,ix,kx) = emat_e05de(ie,je,ix,kx) &
                              + 0.5*zarg0 * zarg1 * (mygamma2(xarg) / zarg2) &
                              * (e_alpha*me + xi_beta_l(jx))

                      endif
                   enddo

                   do kx=1,3
                      if(kx == 1) then 
                         jx = ix-2
                      else if(kx == 2) then
                         jx = ix
                      else 
                         jx = ix+2
                      endif

                      if (jx >= 0 .and. jx <= n_xi) then

                         zarg0 = (-1.0)**(ke+me)
                         zarg1 = (mygamma2(2 + 2*ie + e_lag(ix)) &
                              / mygamma2(2 + 2*(ie-ke)) &
                              / mygamma2(2 + 2*ke + e_lag(ix))) &
                              * (mygamma2(2 + 2*je + e_lag(jx)) &
                              / mygamma2(2 + 2*(je-me)) &
                              / mygamma2(2 + 2*me + e_lag(jx)))
                         zarg2 = mygamma2(2 + 2*ke) * mygamma2(2 + 2*me)

                         xarg = e_alpha*(ke+me) + xi_beta_l(ix) &
                              + xi_beta_l(jx) + 5
                         emat_e1(ie,je,ix,kx) = emat_e1(ie,je,ix,kx) &
                              + 0.5*zarg0 * zarg1 * (mygamma2(xarg) / zarg2)
                         
                      endif
                   enddo


                enddo
             enddo
          enddo
       enddo
    enddo
    
  end subroutine ENERGY_basis_ints

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
    use neo_globals
    implicit none
    integer, intent(in) :: ir
    real :: tauinv_ab, tauinv_ba, lambda
    integer :: xarg, yarg
    real :: rs, ru, rd, rh, rk, rp, rpi, rv, r1, r2, r3, r4, r5, r6, r7
    real, dimension(:,:), allocatable :: fcoll, fcoll_bar, &
         fcollinv, fcollinv_bar
    integer :: fmarg
    integer :: is, js, ie, je, ix, jx, ke, me
    real :: zarg0, zarg1, zarg2
    integer :: is_ele

    emat_coll_test(:,:,:,:,:)  = 0.0
    emat_coll_field(:,:,:,:,:) = 0.0

    fmarg = 2*e_alpha*n_energy + (2*n_xi)+4
    allocate(fcoll(-fmarg:fmarg,-fmarg:fmarg))
    allocate(fcoll_bar(-fmarg:fmarg,-fmarg:fmarg))
    allocate(fcollinv(-fmarg:fmarg,-fmarg:fmarg))
    allocate(fcollinv_bar(-fmarg:fmarg,-fmarg:fmarg))

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

                   do ke=0,ie
                      do me=0,je

                         jx = ix

                         zarg0 = (-1.0)**(ke+me)
                         zarg1 = (mygamma2(2 + 2*ie + e_lag(ix)) &
                              / mygamma2(2 + 2*(ie-ke)) &
                              / mygamma2(2 + 2*ke + e_lag(ix))) &
                              * (mygamma2(2 + 2*je + e_lag(jx)) &
                              / mygamma2(2 + 2*(je-me)) &
                              / mygamma2(2 + 2*me + e_lag(jx)))
                         zarg2 = mygamma2(2 + 2*ke) * mygamma2(2 + 2*me)

                         if(collision_model == 1 .or. &
                              (spitzer_model==1 .and. js .ne. is)) then
                            ! Connor model

                            if(is == js .or. &
                                 (abs(mass(is)-mass(js)) < epsilon(0.))) then
                               ! case 1: ma = mb
                               xarg = e_alpha*(ke+me) &
                                    + xi_beta_l(ix) + xi_beta_l(jx)
                               emat_coll_test(is,js,ie,je,ix) = &
                                    emat_coll_test(is,js,ie,je,ix) &
                                    -tauinv_ab &
                                    * sqrt(lambda/pi) *ix*(ix+1) &
                                    * zarg0 * zarg1 &
                                    * (fcoll(xarg-1,0) - fcoll(xarg-3,2)) &
                                    / zarg2

                               if(ix == 1) then
                                  xarg = e_alpha*ke + xi_beta_l(ix)
                                  yarg = e_alpha*me + xi_beta_l(jx)
                                  emat_coll_field(is,js,ie,je,ix) = &
                                       emat_coll_field(is,js,ie,je,ix) &
                                       + (mass(js)*dens(js,ir)*vth(js,ir)) &
                                       / (mass(is)*dens(is,ir)*vth(is,ir)) &
                                       * 2.0*tauinv_ba &
                                       / sqrt(lambda*pi) &
                                       * zarg0 * zarg1 &
                                       * (fcoll(xarg,0) - fcoll(xarg-2,2)) &
                                       / zarg2 &
                                       / (fcoll(1,0) - fcoll(-1,2)) &
                                       * (fcollinv(yarg,0) - fcollinv(yarg-2,2))

                               endif

                            else if(mass(is) < mass(js)) then
                               ! case 2: ma < mb
                               ! e-i, i-z
                               xarg = e_alpha*(ke+me) &
                                    + xi_beta_l(ix) + xi_beta_l(jx)
                               if(xarg > 0) then
                                  emat_coll_test(is,js,ie,je,ix) = &
                                       emat_coll_test(is,js,ie,je,ix) &
                                       -tauinv_ab &
                                       * 0.25 *ix*(ix+1) &
                                       * zarg0 * zarg1 &
                                       * mygamma2(xarg) &
                                       / zarg2
                               endif
                               if(ix == 1) then
                                  xarg = e_alpha*ke + xi_beta_l(ix)
                                  yarg = e_alpha*me + xi_beta_l(jx)
                                  emat_coll_field(is,js,ie,je,ix) = &
                                       emat_coll_field(is,js,ie,je,ix) &
                                       + (mass(js)*dens(js,ir)*vth(js,ir)) &
                                       / (mass(is)*dens(is,ir)*vth(is,ir)) &
                                       * (2.0/3.0)*tauinv_ba &
                                       * temp(js,ir) / temp(is,ir) &
                                       / sqrt(lambda*pi) &
                                       * zarg0 * zarg1 &
                                       * mygamma2(xarg+1) / mygamma2(2) &
                                       / zarg2 &
                                       * mygamma2(yarg+4)
                               endif

                            else 
                               ! case 3: ma > mb
                               ! i-e, z-i
                               xarg = e_alpha*(ke+me) &
                                    + xi_beta_l(ix) + xi_beta_l(jx)
                               emat_coll_test(is,js,ie,je,ix) = &
                                    emat_coll_test(is,js,ie,je,ix) &
                                    -tauinv_ab &
                                    * sqrt(lambda/pi) &
                                    * temp(is,ir)/temp(js,ir) &
                                    * (1.0/3.0) *ix*(ix+1) &
                                    * zarg0 * zarg1 &
                                    * mygamma2(xarg+3) &
                                    / zarg2
                               if(ix == 1) then
                                  xarg = e_alpha*ke + xi_beta_l(ix)
                                  yarg = e_alpha*me + xi_beta_l(jx)
                                  emat_coll_field(is,js,ie,je,ix) = &
                                       emat_coll_field(is,js,ie,je,ix) &
                                       + (mass(js)*dens(js,ir)*vth(js,ir)) &
                                       / (mass(is)*dens(is,ir)*vth(is,ir)) &
                                       * 0.5 * tauinv_ba &
                                       * temp(js,ir) / temp(is,ir) &
                                       * zarg0 * zarg1 &
                                       * mygamma2(xarg+4) / mygamma2(5) &
                                       / zarg2 &
                                       * mygamma2(yarg+1)
                               endif

                            endif

                         else if(collision_model == 2) then
                            ! HS0
                            xarg = e_alpha*(ke+me) &
                                 + xi_beta_l(ix) + xi_beta_l(jx)
                            emat_coll_test(is,js,ie,je,ix) = &
                                 emat_coll_test(is,js,ie,je,ix) &
                                 -tauinv_ab &
                                 * sqrt(lambda/pi) *ix*(ix+1) &
                                 * zarg0 * zarg1 &
                                 * (fcoll(xarg-1,0) - fcoll(xarg-3,2)) &
                                 / zarg2

                            if(ix == 1) then
                               xarg = e_alpha*ke + xi_beta_l(ix)
                               yarg = e_alpha*me + xi_beta_l(jx)
                               rs = (mass(js)*dens(js,ir)*vth(js,ir)) &
                                    / (mass(is)*dens(is,ir)*vth(is,ir)) &
                                    * 4.0*tauinv_ba &
                                    * temp(js,ir)/temp(is,ir) &
                                    * (1.0 + mass(is)/mass(js)) &
                                    / sqrt(lambda*pi) &
                                    * zarg0 * zarg1 * fcoll(xarg,2) / zarg2 &
                                    / fcoll(1,2) * fcollinv(yarg,2)
                               xarg = e_alpha*(ke+me) &
                                    + xi_beta_l(ix) + xi_beta_l(jx)
                               ru = 2.0*tauinv_ab * sqrt(lambda/pi) &
                                    * zarg0 * zarg1 &
                                    * (fcoll(xarg-1,0)-fcoll(xarg-3,2) &
                                    - 2.0*temp(is,ir)/temp(js,ir) &
                                    * (1.0 + mass(js)/mass(is)) &
                                    * fcoll(xarg-1,2)) &
                                    / zarg2 
                               emat_coll_field(is,js,ie,je,ix) = &
                                    emat_coll_field(is,js,ie,je,ix) + rs
                               emat_coll_test(is,js,ie,je,ix) = &
                                    emat_coll_test(is,js,ie,je,ix) + ru
                            endif

                         else if(collision_model == 3) then
                            ! Full HS
                            xarg = e_alpha*(ke+me) &
                                 + xi_beta_l(ix) + xi_beta_l(jx)
                            emat_coll_test(is,js,ie,je,ix) = &
                                 emat_coll_test(is,js,ie,je,ix) &
                                 -tauinv_ab &
                                 * sqrt(lambda/pi) *ix*(ix+1) &
                                 * zarg0 * zarg1 &
                                 * (fcoll(xarg-1,0) - fcoll(xarg-3,2)) &
                                 / zarg2

                            if(ix == 1) then
                               ! slowing-down
                               xarg = e_alpha*ke + xi_beta_l(ix)
                               yarg = e_alpha*me + xi_beta_l(jx)
                               rs = (mass(js)*dens(js,ir)*vth(js,ir)) &
                                    / (mass(is)*dens(is,ir)*vth(is,ir)) &
                                    * 4.0*tauinv_ba &
                                    * temp(js,ir)/temp(is,ir) &
                                    * (1.0 + mass(is)/mass(js)) &
                                    / sqrt(lambda*pi) &
                                    * zarg0 * zarg1 * fcoll(xarg,2) / zarg2 &
                                    / fcoll(1,2) * fcollinv(yarg,2)
                               ! u-restoring
                               xarg = e_alpha*(ke+me) &
                                    + xi_beta_l(ix) + xi_beta_l(jx)
                               ru = 2.0*tauinv_ab * sqrt(lambda/pi) &
                                    * zarg0 * zarg1 &
                                    * (fcoll(xarg-1,0)-fcoll(xarg-3,2) &
                                    - 2.0*temp(is,ir)/temp(js,ir) &
                                    * (1.0 + mass(js)/mass(is)) &
                                    * fcoll(xarg-1,2)) &
                                    / zarg2
                               ! h-heating friction
                               xarg = e_alpha*ke + xi_beta_l(ix)
                               r1 = 4.0*(2.0*mass(is)/mass(js)-1.0) &
                                    * fcoll(xarg,2) &
                                    + 8.0*temp(is,ir)/temp(js,ir) &
                                    * (1.0 + mass(js)/mass(is)) &
                                    * fcoll(xarg+2,2) &
                                    - 4.0*(1.0 + mass(is)/mass(js)) &
                                    * (fcoll(xarg+2,0)-fcoll(xarg,2))
                               xarg = 3
                               r2 = 4.0*(2.0*mass(is)/mass(js)-1.0) &
                                    * fcoll(xarg,2) &
                                    + 8.0*temp(is,ir)/temp(js,ir) &
                                    * (1.0 + mass(js)/mass(is)) &
                                    * fcoll(xarg+2,2) &
                                    - 4.0*(1.0 + mass(is)/mass(js)) &
                                    * (fcoll(xarg+2,0)-fcoll(xarg,2))
                               xarg = e_alpha*me + xi_beta_l(jx)
                               r3 = 4.0*(2.0*mass(js)/mass(is)-1.0) &
                                    * fcollinv(xarg,2) &
                                    + 8.0*temp(js,ir)/temp(is,ir) &
                                    * (1.0 + mass(is)/mass(js)) &
                                    * fcollinv(xarg+2,2) &
                                    - 4.0*(1.0 + mass(js)/mass(is)) &
                                    * (fcollinv(xarg+2,0)-fcollinv(xarg,2))
                               rh = (mass(js)*dens(js,ir)*vth(js,ir)**3) &
                                    / (mass(is)*dens(is,ir)*vth(is,ir)**3) &
                                    * 1.5*tauinv_ba / sqrt(lambda*pi) &
                                    * zarg0 * zarg1 * r1 / zarg2 &
                                    / r2 * r3
                               ! k-heating friction
                               xarg = e_alpha*ke + xi_beta_l(ix)
                               r1 = 12.0 * fcoll(xarg,2) &
                                    - 8.0*temp(is,ir)/temp(js,ir) &
                                    * (1.0 + mass(js)/mass(is)) &
                                    * fcoll(xarg+2,2) &
                                    + 4.0*(fcoll(xarg+2,0)-fcoll(xarg,2))
                               xarg = 3
                               r2 = 12.0 * fcoll(xarg,2) &
                                    - 8.0*temp(is,ir)/temp(js,ir) &
                                    * (1.0 + mass(js)/mass(is)) &
                                    * fcoll(xarg+2,2) &
                                    + 4.0*(fcoll(xarg+2,0)-fcoll(xarg,2))
                               xarg = e_alpha*me + xi_beta_l(jx)
                               r3 = 12.0 * fcoll(xarg,2) &
                                    - 8.0*temp(is,ir)/temp(js,ir) &
                                    * (1.0 + mass(js)/mass(is)) &
                                    * fcoll(xarg+2,2) &
                                    + 4.0*(fcoll(xarg+2,0)-fcoll(xarg,2))
                               rk = tauinv_ab * sqrt(lambda/pi) &
                                    * zarg0 * zarg1 * r1 / zarg2 &
                                    / r2 * r3
                               emat_coll_field(is,js,ie,je,ix) = &
                                    emat_coll_field(is,js,ie,je,ix) + rs + rh
                               emat_coll_test(is,js,ie,je,ix) = &
                                    emat_coll_test(is,js,ie,je,ix) + ru + rk
                            endif

                            if(ix == 2) then
                               ! nup-energy restoring
                               xarg = e_alpha*ke + xi_beta_l(ix)
                               r1 = 8.0 * (1.0 + 1.5*mass(is)/mass(js)) &
                                    * fcoll(xarg-1,2) &
                                    + 8.0*temp(is,ir)/temp(js,ir) &
                                    * (1.0 + mass(js)/mass(is)) &
                                    * fcoll(xarg+1,2) &
                                    - 4.0 * (1.0 + 1.5*mass(is)/mass(js)) &
                                    * (fcoll(xarg+1,0)-fcoll(xarg-1,2))
                               xarg = 2
                               r2 = 8.0 * (1.0 + 1.5*mass(is)/mass(js)) &
                                    * fcoll(xarg-1,2) &
                                    + 8.0*temp(is,ir)/temp(js,ir) &
                                    * (1.0 + mass(js)/mass(is)) &
                                    * fcoll(xarg+1,2) &
                                    - 4.0 * (1.0 + 1.5*mass(is)/mass(js)) &
                                    * (fcoll(xarg+1,0)-fcoll(xarg-1,2))
                               xarg = e_alpha*me + xi_beta_l(jx)
                               r3 = 8.0 * (1.0 + 1.5*mass(js)/mass(is)) &
                                    * fcollinv(xarg-1,2) &
                                    + 8.0*temp(js,ir)/temp(is,ir) &
                                    * (1.0 + mass(is)/mass(js)) &
                                    * fcollinv(xarg+1,2) &
                                    - 4.0 * (1.0 + 1.5*mass(js)/mass(is)) &
                                    * (fcollinv(xarg+1,0)-fcollinv(xarg-1,2))
                               rp = (temp(js,ir)*dens(js,ir)) &
                                    / (temp(is,ir)*dens(is,ir)) &
                                    * tauinv_ba / sqrt(lambda*pi) &
                                    * zarg0 * zarg1 * r1 / zarg2 &
                                    / r2 * r3
                               emat_coll_field(is,js,ie,je,ix) = &
                                    emat_coll_field(is,js,ie,je,ix) + rp

                               ! pi-energy restoring
                               xarg = e_alpha*ke + xi_beta_l(ix)
                               yarg = e_alpha*me + xi_beta_l(jx)
                               r1 = -8.0*temp(is,ir)/temp(js,ir) &
                                    * (1.0 + mass(js)/mass(is)) &
                                    * fcoll(xarg+yarg-1,2) &
                                    + 4.0 * fcoll(xarg+yarg-1,0)
                               rpi = tauinv_ab * sqrt(lambda/pi) &
                                    * zarg0 * zarg1 * r1 / zarg2
                               emat_coll_test(is,js,ie,je,ix) = &
                                    emat_coll_test(is,js,ie,je,ix) + rpi
                            endif

                            if(ix == 0) then
                               ! energy diffusion
                               xarg = e_alpha*me + xi_beta_l(jx) 
                               r2 = 8.0*temp(js,ir)/temp(is,ir) &
                                    * (1.0 + mass(is)/mass(js)) &
                                    * fcollinv(xarg+1,2) &
                                    - 4.0 * fcollinv(xarg+1,0)
                               xarg = 2
                               r3 = 8.0*temp(is,ir)/temp(js,ir) &
                                    * (1.0 + mass(js)/mass(is)) &
                                    * fcoll(xarg+1,2) &
                                    - 4.0 * fcoll(xarg+1,0)
                               r1 = (temp(js,ir)*dens(js,ir)) &
                                    / (temp(is,ir)*dens(is,ir)) &
                                    * tauinv_ba / tauinv_ab / lambda * r2 / r3
                               xarg = e_alpha*ke + xi_beta_l(ix) 
                               yarg = e_alpha*me + xi_beta_l(jx)
                               rd =  -4.0*tauinv_ab * sqrt(lambda/pi) &
                                    * xarg * zarg0 * zarg1 &
                                    * 0.5 * yarg &
                                    * fcoll(xarg+yarg-3,2) / zarg2
                               rv =  4.0*tauinv_ab * sqrt(lambda/pi) &
                                    * xarg * zarg0 * zarg1 &
                                    * r1 * fcoll(xarg-1,2) / zarg2
                               emat_coll_test(is,js,ie,je,ix) = &
                                    emat_coll_test(is,js,ie,je,ix) + rd
                               emat_coll_field(is,js,ie,je,ix) = &
                                    emat_coll_field(is,js,ie,je,ix) + rv
                            endif

                         else 
                            ! Full linearized FP op
                            xarg = e_alpha*(ke+me) &
                                 + xi_beta_l(ix) + xi_beta_l(jx)
                            emat_coll_test(is,js,ie,je,ix) = &
                                 emat_coll_test(is,js,ie,je,ix) &
                                 -tauinv_ab &
                                 * sqrt(lambda/pi) *ix*(ix+1) &
                                 * zarg0 * zarg1 &
                                 * (fcoll(xarg-1,0) - fcoll(xarg-3,2)) &
                                 / zarg2

                            rd = 4.0*tauinv_ab * sqrt(lambda/pi) &
                                 * (e_alpha*ke + xi_beta_l(ix)) &
                                 * zarg0 * zarg1 &
                                 * ( (1.0 - temp(is,ir)/temp(js,ir)) &
                                 * fcoll(xarg-1,2) &
                                 - 0.5 * (e_alpha*me + xi_beta_l(jx)) &
                                 * fcoll(xarg-3,2) ) &
                                 / zarg2

                            emat_coll_test(is,js,ie,je,ix) = &
                                 emat_coll_test(is,js,ie,je,ix) + rd

                            ! Real full field-particle operator

                            r1 = mass(is)/mass(js) &
                                 * (1.0 + lambda)**(-0.5*(xarg+3)) &
                                 * mygamma2(xarg + 3)

                            xarg = e_alpha*ke + xi_beta_l(ix)
                            yarg = e_alpha*me + xi_beta_l(jx)

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

                            emat_coll_field(is,js,ie,je,ix) = &
                                 emat_coll_field(is,js,ie,je,ix) &
                                 + tauinv_ab * 2.0/sqrt(pi) * lambda**1.5 &
                                 * lambda**(0.5*(e_alpha*me &
                                 + xi_beta_l(ix))) &
                                 * zarg0 * zarg1 &
                                 * (r1 + r2 + r3 + r4 + r5 + r6 + r7) &
                                 / zarg2
                         endif
                      enddo
                   enddo
                enddo
             enddo
          enddo
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
       is_ele = -1
       do is=1, n_species
          if(Z(is) < 0.0) then
             is_ele = is
             exit
          endif
       enddo
       if(is_ele == -1) then
          call neo_error('ERROR: (NEO) Must have electron species for uncoupled e-i problem')
          return
       endif

       do is=1,n_species
          do js=1,n_species

             if(is .ne. is_ele .and. js == is_ele) then
                if(coll_uncoupledei_model == 1) then
                   ! C_ie = 0
                   emat_coll_test(is,js,:,:,:) = 0.0
                   emat_coll_field(is,js,:,:,:) = 0.0
                else
                   ! C_ie = Connor case ma>mb
                   emat_coll_test(is,js,:,:,:) = 0.0
                   emat_coll_field(is,js,:,:,:) = 0.0
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
                            do ke=0,ie
                               do me=0,je
                                  jx = ix
                                  zarg0 = (-1.0)**(ke+me)
                                  zarg1 = (mygamma2(2 + 2*ie + e_lag(ix)) &
                                       / mygamma2(2 + 2*(ie-ke)) &
                                       / mygamma2(2 + 2*ke + e_lag(ix))) &
                                       * (mygamma2(2 + 2*je + e_lag(jx)) &
                                       / mygamma2(2 + 2*(je-me)) &
                                       / mygamma2(2 + 2*me + e_lag(jx)))
                                  zarg2 = mygamma2(2 + 2*ke) * mygamma2(2 + 2*me)
                                  xarg = e_alpha*(ke+me) &
                                       + xi_beta_l(ix) + xi_beta_l(jx)
                                  emat_coll_test(is,js,ie,je,ix) = &
                                       emat_coll_test(is,js,ie,je,ix) &
                                       -tauinv_ab &
                                       * sqrt(lambda/pi) &
                                       * temp(is,ir)/temp(js,ir) &
                                       * (1.0/3.0) *ix*(ix+1) &
                                       * zarg0 * zarg1 &
                                       * mygamma2(xarg+3) &
                                       / zarg2
                                  if(ix == 1) then
                                     xarg = e_alpha*ke + xi_beta_l(ix)
                                     yarg = e_alpha*me + xi_beta_l(jx)
                                     emat_coll_field(is,js,ie,je,ix) = &
                                          emat_coll_field(is,js,ie,je,ix) &
                                          + (mass(js)*dens(js,ir)*vth(js,ir)) &
                                          / (mass(is)*dens(is,ir)*vth(is,ir)) &
                                          * 0.5 * tauinv_ba &
                                          * temp(js,ir) / temp(is,ir) &
                                          * zarg0 * zarg1 &
                                          * mygamma2(xarg+4) / mygamma2(5) &
                                          / zarg2 &
                                          * mygamma2(yarg+1)
                                  endif
                               enddo
                            enddo
                         enddo
                      enddo
                   enddo
                endif
             endif
             if(is == is_ele .and. js .ne. is) then
                ! C_ei = L(fe) with Connor case ma<mb
                emat_coll_field(is,js,:,:,:) = 0.0
                emat_coll_test(is,js,:,:,:) = 0.0
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
                         do ke=0,ie
                            do me=0,je
                               jx = ix
                               zarg0 = (-1.0)**(ke+me)
                               zarg1 = (mygamma2(2 + 2*ie + e_lag(ix)) &
                                    / mygamma2(2 + 2*(ie-ke)) &
                                    / mygamma2(2 + 2*ke + e_lag(ix))) &
                                    * (mygamma2(2 + 2*je + e_lag(jx)) &
                                    / mygamma2(2 + 2*(je-me)) &
                                    / mygamma2(2 + 2*me + e_lag(jx)))
                               zarg2 = mygamma2(2 + 2*ke) * mygamma2(2 + 2*me)
                               xarg = e_alpha*(ke+me) &
                                    + xi_beta_l(ix) + xi_beta_l(jx)
                               if(xarg > 0) then
                                  emat_coll_test(is,js,ie,je,ix) = &
                                       emat_coll_test(is,js,ie,je,ix) &
                                       -tauinv_ab &
                                       * 0.25 *ix*(ix+1) &
                                       * zarg0 * zarg1 &
                                       * mygamma2(xarg) &
                                       / zarg2
                               endif
                               if(ix == 1) then
                                  xarg = e_alpha*ke + xi_beta_l(ix)
                                  yarg = e_alpha*me + xi_beta_l(jx)
                                  emat_coll_field(is,js,ie,je,ix) = &
                                       emat_coll_field(is,js,ie,je,ix) &
                                       + (mass(js)*dens(js,ir)*vth(js,ir)) &
                                       / (mass(is)*dens(is,ir)*vth(is,ir)) &
                                       * (2.0/3.0)*tauinv_ba &
                                       * temp(js,ir) / temp(is,ir) &
                                       / sqrt(lambda*pi) &
                                       * zarg0 * zarg1 &
                                       * mygamma2(xarg+1) / mygamma2(2) &
                                       / zarg2 &
                                       * mygamma2(yarg+4)
                               endif
                            enddo
                         enddo
                      enddo
                   enddo
                enddo
             endif
          enddo
       enddo
    endif

    deallocate(fcoll)
    deallocate(fcoll_bar)
    deallocate(fcollinv)
    deallocate(fcollinv_bar)

    if (write_cmoments_flag == 1) then
       call write_fullcoll_mono(ir)
    endif

  end subroutine ENERGY_coll_ints

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
