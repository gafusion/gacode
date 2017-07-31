!-----------------------------------------------------------------
! neo_do.f90
!
! PURPOSE:
!  Subroutinized main neo program.  
!
! NOTES:
!  This can be called directly using the driver routine neo 
!  (in which case input data will read from input.dat) or called 
!  as a subroutine using neo_sub.
!-----------------------------------------------------------------

subroutine neo_do

  use neo_globals
  use neo_energy_grid
  use neo_sparse_solve
  use neo_equilibrium
  use neo_rotation
  use neo_transport
  use neo_theory
  use neo_nclass_dr
  use neo_g_velocitygrids
  use neo_allocate_profile
  use neo_3d_driver
  use mpi
  implicit none

  integer :: ir, is, ie, ix, it, js, je, jx, jt, ks
  integer, dimension(:,:,:,:), allocatable :: mindx  ! (ns,ne,nxi+1,nth)
  integer :: i, j, k, id
  integer :: ierr
  integer :: n_elem, asize
  real, dimension(:), allocatable :: a
  integer, dimension(:), allocatable :: a_iindx
  integer, dimension(:), allocatable :: a_jindx
  integer :: ifac, matfac_err, max_ifac=5

  ! kinetic equation terms: 
  real :: stream, trap, rotkin

  integer, parameter :: io_neo=10, io_f=11
  character(len=80)  :: runfile_f = 'out.neo.f'

  ! First, define some grid dimensions for band-diagonal LAPACK
  if (threed_model == 1) then
     nb    = (n_energy+1)*n_species*tpmatsize
     n_row = (n_xi+1)*nb
  else
     nb    = (n_energy+1)*n_species*n_theta
     n_row = (n_xi+1)*nb
  endif

  call neo_make_profiles
  if(error_status > 0) goto 100
  call neo_check
  if(error_status > 0) goto 100

  if(spitzer_model==1) then
     call neo_spitzer
     goto 100
  endif

  if(threed_model==1) then
     call ThreeD_do
     goto 100
  endif

  ! cyclic index (for theta-periodicity)
  allocate(thcyc(1-n_theta:2*n_theta))
  do it=1,n_theta
     thcyc(it-n_theta) = it
     thcyc(it) = it
     thcyc(it+n_theta) = it
  enddo
  ! coefficients for 4th order centered derivative
  cderiv(-2) =  1
  cderiv(-1) = -8
  cderiv(0)  =  0
  cderiv(1)  =  8
  cderiv(2)  = -1

  if(sim_model == 0 .or. sim_model == 3) then
     ! Theory calculation only -- no numerical kinetic calculation
     call EQUIL_alloc(1)
     call THEORY_alloc(1)
     do ir=1, n_radial
        call EQUIL_do(ir)
        call THEORY_do(ir)
        if(error_status > 0) goto 100
        ! Store the local neo transport values at ir=1 in neo_xtheory_out
        ! (n_species_max, transport coeff)
        ! transport coeff: 1-> gamma, 2-> Q, 3->Pi, 4-> upar
        if(ir == 1) then
           neo_dke_out(:,:)  = 0.0
           neo_dke_1d_out(:) = 0.0
           neo_gv_out(:,:)   = 0.0
           neo_th_out(:)     = 0.0
           neo_th_out(1)     = pflux_HH
           neo_th_out(2)     = efluxi_HH
           do is=1, n_species
              if(Z(is) == -1) then
                 neo_th_out(3) = efluxe_HH
              endif
           enddo
           neo_th_out(4)     = efluxi_CH
           neo_th_out(5)     = jpar_S
           neo_th_out(6)     = jpar_K
           neo_th_out(8)     = jtor_S
           neo_thHS_out(:,:) = 0.0
           do is=1, n_species
              neo_thHS_out(is,1) = pflux_multi_HS(is) 
              neo_thHS_out(is,2) = eflux_multi_HS(is) 
           enddo
        end if
     enddo
     call THEORY_alloc(0)
     call EQUIL_alloc(0)
     deallocate(thcyc)
     return
  end if

  ! Set-up the energy grids for basis ints (indep of r)
  call ENERGY_basis_ints_alloc(1)
  call ENERGY_basis_ints
  call ENERGY_coll_ints_alloc(1)

  ! Matrix solve allocations
  n_row = n_species*(n_energy+1)*(n_xi+1)*n_theta
  ! Estimate number of elements
  i = 0
  ! constraint
  if(collision_model == 1 .or. collision_model == 2) then
     i = i + n_species*n_theta * (n_energy+1)
  else
     i = i + n_species*n_theta * 2
  endif
  ! collisions
  i = i + n_species*(n_xi+1)*n_theta*(n_energy+1)**2 * (1 + n_species)
  ! streaming
  i = i + n_species*n_theta*(n_energy+1)**2*((n_xi-1)*4*2 + 2*4)
  ! trapping/rotation
  i = i + n_species*n_theta*(n_energy+1)**2*((n_xi-1)*2 + 2)
  n_max = i
  asize = n_max
  allocate(a(n_max),stat=ierr)
  if(ierr /= 0) then
     call neo_error('ERROR: (NEO) Array allocation failed')
     goto 100
  end if
  allocate(a_iindx(n_max),stat=ierr)
  if(ierr /= 0) then
     call neo_error('ERROR: (NEO) Array allocation failed')
     goto 100
  end if
  allocate(a_jindx(n_max),stat=ierr)
  if(ierr /= 0) then
     call neo_error('ERROR: (NEO) Array allocation failed')
     goto 100
  end if
  allocate(g(n_row))
  allocate(mindx(n_species,0:n_energy,0:n_xi,n_theta))
  allocate(is_indx(n_row))
  allocate(ie_indx(n_row))
  allocate(ix_indx(n_row))
  allocate(it_indx(n_row))

  ! matrix indices
  i = 0
  do is=1,n_species
     do ie=0,n_energy
        do ix=0,n_xi
           do it=1,n_theta
              i = i+1
              mindx(is,ie,ix,it) = i
              is_indx(i) = is
              ie_indx(i) = ie
              ix_indx(i) = ix
              it_indx(i) = it
           enddo
        enddo
     enddo
  enddo

  allocate(driftx(n_species,n_theta))
  allocate(driftxrot1(n_species,n_theta))
  allocate(driftxrot2(n_species,n_theta))
  allocate(driftxrot3(n_species,n_theta))

  call EQUIL_alloc(1)
  call ROT_alloc(1)
  call TRANSP_alloc(1)

  ! Set-up the matrix equation: LHS radially local matrix

  if (silent_flag == 0 .and. i_proc == 0) then
     open(unit=io_neo,file=trim(path)//'out.neo.grid',status='replace')
     write(io_neo,*) n_species
     write(io_neo,*) n_energy
     write(io_neo,*) n_xi
     write(io_neo,*) n_theta
     do it=1, n_theta
        write(io_neo,*) theta(it)
     enddo
     write(io_neo,*) n_radial
     do ir=1, n_radial
        write(io_neo,*) r(ir)
     enddo
     close(io_neo)

     open(unit=io_f,file=trim(path)//runfile_f,status='replace')
     close(io_f)
  end if

  do ir=1, n_radial

     if(collision_model == 3) then
        do is=2, n_species
           if(abs(temp(is,ir)-temp(1,ir)) > 1e-3) then
              if(silent_flag == 0 .and. i_proc == 0) then
                 open(unit=io_neoout,file=trim(path)//runfile_neoout,status='old',position='append')
                 write(io_neoout,*) 'WARNING: (NEO) Full HS collisions with unequal temps'
                 close(io_neoout)
              endif
           endif
        enddo
     endif

     ! Get the equilibrium parameters (th)
     call EQUIL_do(ir)

     ! Get the rotation phi
     call ROT_solve_phi(ir)
     if(error_status > 0) goto 100

     ! Get the energy grids for coll ints
     call ENERGY_coll_ints(ir)

     ! Set the LHS
     a_iindx(:) = 0
     a_jindx(:) = 0
     a(:) = 0.0
     k = 0

     do is=1,n_species
        do ie=0,n_energy
           do ix=0,n_xi
              do it=1,n_theta

                 ! vpar bhat dot grad 
                 ! -> stream * xi * sqrt(ene)
                 stream = sqrt(2.0) * vth(is,ir) * k_par(it) &
                      / (12*d_theta)
                 ! mu bdot grad B d/dvpar 
                 ! -> trap * (1-xi^2)*d/dxi * sqrt(ene)
                 trap   = ( gradpar_Bmag(it) / Bmag(it) ) &
                      * sqrt(0.5) * vth(is,ir)
                 ! vdrift dot grad r 
                 ! -> driftx * (1+xi^2) * d/d(r/a) * ene(ie)
                 driftx(is,it) = v_drift_x(it) &
                      * mass(is)/(1.0*Z(is)) * (vth(is,ir))**2 
                 ! rotation
                 rotkin = 0.5 * sqrt(2.0) * vth(is,ir) *  k_par(it) &
                      * (-Z(is)/temp(is,ir)*phi_rot_deriv(it) &
                      + omega_rot(ir)**2 * bigR(it) / vth(is,ir)**2 &
                      * bigR_tderiv(it))
                 driftxrot1(is,it) = I_div_psip * k_par(it) &
                      * mass(is)/(1.0*Z(is)) * rho(ir) / Bmag(it) &
                      * (vth(is,ir))**2 &
                      * (-Z(is)/temp(is,ir)*phi_rot_deriv(it) &
                      + omega_rot(ir)**2 * bigR(it)/ vth(is,ir)**2 &
                      * bigR_tderiv(it))
                 if(aniso_model(is) == 2) then
                    rotkin = rotkin - 0.5 * sqrt(2.0) * vth(is,ir) &
                         * lam_rot_kpar_aniso(is,it)
                    driftxrot1(is,it) = driftxrot1(is,it) - I_div_psip &
                         * mass(is)/(1.0*Z(is)) * rho(ir) / Bmag(it) &
                         * (vth(is,ir))**2 * lam_rot_kpar_aniso(is,it)
                 endif
                 driftxrot2(is,it) = I_div_psip* k_par(it) / Btor(it) &
                      * mass(is)/(1.0*Z(is)) * rho(ir) &
                      * vth(is,ir) * 2.0 * sqrt(2.0) &
                      * bigR_tderiv(it) * omega_rot(ir)
                 driftxrot3(is,it) = 1.0/sqrt(2.0) &
                      * vth(is,ir)**2 * mass(is)/(1.0*Z(is)) &
                      * rho(ir) / Bmag(it) &
                      * Btor(it)/(Bmag(it)* I_div_psip) &
                      * (2.0*k_par(it) * gradr(it) * gradr_tderiv(it) &
                      - gradpar_Bmag(it) / Bmag(it) * gradr(it)**2)

                 i = mindx(is,ie,ix,it)

                 ! Impose constant for g0
                 ! local  -> at each species and ene 

                 if (ix==0 .and. it == 1 &
                      .and. (collision_model==1 .or. collision_model==2))&
                      then
                    ! <f_ie> = 0
                    js = is; je = ie; jx=0
                    do jt=1, n_theta
                       j = mindx(js,je,jx,jt)
                       k = k+1
                       a(k) = w_theta(jt)
                       a_iindx(k) = i
                       a_jindx(k) = j
                    enddo

                 else if(ix==0 .and. it == 1 &
                      .and. (ie==0 .or. ie==1) &
                      .and. (collision_model==3 .or. &
                      collision_model==4 .or. collision_model==5)) then
                    ! <f_ie> = 0
                    js = is; je = ie; jx=0
                    do jt=1, n_theta
                       j = mindx(js,je,jx,jt)
                       k = k+1
                       a(k) = w_theta(jt)
                       a_iindx(k) = i
                       a_jindx(k) = j
                    enddo

                 else

                    ! Collisions 
                    jt = it; jx = ix

                    ! test particle
                    js = is
                    do je=0, n_energy
                       j = mindx(js,je,jx,jt)
                       k = k+1
                       a(k) = 0.0
                       do ks=1, n_species
                          if (coll_uncoupledaniso_model == 1) then
                             if(is/= ks .and. aniso_model(ks) == 2) then
                                a(k) = a(k) + 0.0
                             else
                                a(k) = a(k) &
                                     - emat_coll_test(is,ks,ie,je,ix) &
                                     * dens_fac(ks,it)
                             endif
                          else
                             a(k) = a(k) &
                               - emat_coll_test(is,ks,ie,je,ix) &
                               * dens_fac(ks,it)
                          endif
                       enddo
                       a_iindx(k) = i
                       a_jindx(k) = j
                    enddo

                    ! field particle                     
                    do je=0, n_energy
                       do js=1,n_species
                          j = mindx(js,je,jx,jt)
                          k = k+1
                          a(k) = -emat_coll_field(is,js,ie,je,ix) &
                               * dens_fac(js,it)
                          if (coll_uncoupledaniso_model == 1) then
                             if (is/= js .and. aniso_model(js) == 2) then
                                a(k) = 0.0
                             endif
                          endif
                          a_iindx(k) = i
                          a_jindx(k) = j
                       enddo
                    enddo


                    ! Streaming
                    js = is
                    do je=0, n_energy
                       do id=-2,2
                          if (id /= 0) then
                             jt = thcyc(it+id)
                             jx = ix-1
                             if (jx >= 0) then
                                j = mindx(js,je,jx,jt)
                                k = k+1
                                a(k) = stream &
                                     * ix/(2*ix-1.0) * cderiv(id) &
                                     * emat_e05(ie,je,ix,1) 
                                a_iindx(k) = i
                                a_jindx(k) = j
                             endif
                             jx = ix+1
                             if (jx <= n_xi) then
                                j = mindx(js,je,jx,jt)
                                k = k+1
                                a(k) = stream &
                                     * (ix+1.0)/(2*ix+3.0) * cderiv(id) &
                                     * emat_e05(ie,je,ix,2) 
                                a_iindx(k) = i
                                a_jindx(k) = j
                             endif
                          endif
                       enddo
                    enddo

                    ! Trapping and Rotation
                    js = is; jt = it
                    do je=0, n_energy
                       jx = ix-1
                       if (jx >= 0) then
                          j = mindx(js,je,jx,jt)
                          k = k+1
                          a(k) = trap &
                               * ix*(ix-1.0)/(2*ix-1.0) &
                               * emat_e05(ie,je,ix,1) &
                               + rotkin &
                               * ix/(2*ix-1.0) * emat_e05de(ie,je,ix,1) &
                               - rotkin &
                               * ix*(ix-1.0)/(2*ix-1.0) &
                               * emat_en05(ie,je,ix,1)
                          a_iindx(k) = i
                          a_jindx(k) = j                      
                       endif
                       jx = ix+1
                       if (jx <= n_xi) then
                          j = mindx(js,je,jx,jt)
                          k = k+1
                          a(k) = -trap &
                               * (ix+1.0)*(ix+2.0)/(2*ix+3.0) &
                               * emat_e05(ie,je,ix,2) &
                               + rotkin &
                               * (ix+1.0)/(2*ix+3.0) &
                               * emat_e05de(ie,je,ix,2) &
                               + rotkin &
                               * (ix+1.0)*(ix+2.0)/(2*ix+3.0) &
                               * emat_en05(ie,je,ix,2)
                          a_iindx(k) = i
                          a_jindx(k) = j
                       endif
                    enddo

                 endif ! bc if/else
              enddo ! it
           enddo ! ix
        enddo ! ie
     enddo ! is

     ! Factor the Matrix -- uses a(:) and a_indx(:)

     n_elem = k
     n_max = n_elem*matsz_scalefac
     matfac_err = 0

     do ifac = 1, max_ifac

        if(silent_flag == 0 .and. i_proc == 0) then
           open(unit=io_neoout,file=trim(path)//runfile_neoout,&
                status='old',position='append')
           write(io_neoout,'(t2,a,e12.5)') 'Estimated memory (GB) = ', 2*8.0*(asize+n_max)/1.0e9
           close(io_neoout)
        endif
        if(allocated(amat))       deallocate(amat)
        allocate(amat(n_max),stat=ierr)
        if(ierr /= 0) then
           call neo_error('ERROR: (NEO) Array allocation failed')
           goto 100
        end if
        if(allocated(amat_indx))  deallocate(amat_indx)
        allocate(amat_indx(2*n_max),stat=ierr)
        if(ierr /= 0) then
           call neo_error('ERROR: (NEO) Array allocation failed')
           goto 100
        end if

        amat(:) = 0.0
        amat_indx(:) = 0
        do k=1,n_elem
           amat(k) = a(k)
           amat_indx(k) = a_iindx(k)
           amat_indx(n_elem+k) = a_jindx(k)
        enddo

        if(silent_flag == 0 .and. i_proc == 0) then
           open(unit=io_neoout,file=trim(path)//runfile_neoout,&
                status='old',position='append')
           if(ifac == 1) then
              write(io_neoout,'(t2,a)') 'Begin matrix factor'
           else
              write(io_neoout,'(t2,a)') 'Re-trying matrix factorization'
           endif
           close(io_neoout)
        endif
        call SOLVE_factor(n_elem)
        if(error_status > 0) then
           error_status = 0
           n_max = n_max * 2
        else
           matfac_err = 1
           exit
        endif
     enddo
     if(matfac_err == 0) then
        call neo_error('ERROR: (NEO) Matrix factorization failed. Try increasing MATSZ_SCALEFAC.')
        goto 100
     endif
     if(silent_flag == 0 .and. i_proc == 0) then
        open(unit=io_neoout,file=trim(path)//runfile_neoout,&
             status='old',position='append')
        write(io_neoout,'(t2,a)') 'Done matrix factor'
        close(io_neoout)
     endif

     ! Set the RHS source 
     call set_RHS_source

     ! Matrix solve -- uses g(:), a(:), and a_indx(:)
     if(silent_flag == 0 .and. i_proc == 0) then
        open(unit=io_neoout,file=trim(path)//runfile_neoout,&
             status='old',position='append')
        write(io_neoout,'(t2,a)') 'Begin matrix solve'
        close(io_neoout)
     endif
     call SOLVE_do
     if(silent_flag == 0 .and. i_proc == 0) then
        open(unit=io_neoout,file=trim(path)//runfile_neoout,&
             status='old',position='append')
        write(io_neoout,'(t2,a)') 'Done matrix solve'
        close(io_neoout)
     endif

     ! Compute the neo transport coefficients
     call TRANSP_do(ir)
     call TRANSP_write(ir)

     ! re-construct the energy dependence
     ! call g_energy(ir)

     ! Write the rotation parameters
     call ROT_write(ir)

     ! Compute the theory transport coefficients
     if(ir == 1) then
        call THEORY_alloc(1)
     end if
     call  THEORY_do(ir)
     if(error_status > 0) goto 100


     ! Store the local neo transport values at ir=1 in neo_x_out
     ! (n_species_max, transport coeff)
     ! transport coeff: 1-> gamma, 2-> Q, 3->Pi, 4-> upar
     if(ir == 1) then
        neo_dke_out(:,:) = 0.0
        neo_dke_1d_out   = 0.0
        neo_gv_out(:,:)  = 0.0
        do is=1, n_species
           neo_dke_out(is,1) = pflux(is)
           neo_dke_out(is,2) = eflux(is)
           neo_dke_out(is,3) = mflux(is)
           neo_dke_out(is,4) = eflux(is) - omega_rot(ir)*mflux(is) 
           neo_dke_out(is,5) = vpol_th0(is)
           neo_dke_out(is,6) = vtor_th0(is) + vtor_0order_th0
           neo_gv_out(is,1)  = pflux_gv(is)
           neo_gv_out(is,2)  = eflux_gv(is)
           neo_gv_out(is,3)  = mflux_gv(is)
           neo_gv_out(is,4)  = eflux_gv(is) - omega_rot(ir)*mflux_gv(is)
        enddo
        neo_dke_1d_out(1)    = jpar
        neo_dke_1d_out(2)    = jtor
        neo_th_out(:) = 0.0
        neo_th_out(1) = pflux_HH
        neo_th_out(2) = efluxi_HH
        do is=1, n_species
           if(Z(is) == -1) then
              neo_th_out(3) = efluxe_HH
           endif
        enddo
        neo_th_out(4) = efluxi_CH
        neo_th_out(5)     = jpar_S
        neo_th_out(6)     = jpar_K
        neo_th_out(8)     = jtor_S
        if(sim_model == 1) then
           neo_th_out(7)     = jbs_nc
        endif
        neo_thHS_out(:,:) = 0.0
        do is=1, n_species
           neo_thHS_out(is,1) = pflux_multi_HS(is) 
           neo_thHS_out(is,2) = eflux_multi_HS(is) 
        enddo
     end if

     if(ir == n_radial) then
        call THEORY_alloc(0)
     end if

     if (silent_flag == 0 .and. i_proc == 0) then
        open(io_f,file=trim(path)//runfile_f,status='old',position='append')
        write(io_f,'(1pe12.5)') g(:)
        close(io_f)
     endif

     call g_xi

  end do ! ir


  ! Clean-up
100 continue
  call ENERGY_basis_ints_alloc(0)
  call ENERGY_coll_ints_alloc(0)
  call EQUIL_alloc(0)
  call ROT_alloc(0)
  call TRANSP_alloc(0)
  call PROFILE_SIM_alloc(0)
  if(allocated(a))          deallocate(a)
  if(allocated(a_iindx))    deallocate(a_iindx)
  if(allocated(a_jindx))    deallocate(a_jindx)
  if(allocated(amat))       deallocate(amat)
  if(allocated(amat_indx))  deallocate(amat_indx)
  if(allocated(g))          deallocate(g)
  if(allocated(mindx))      deallocate(mindx)
  if(allocated(is_indx))    deallocate(is_indx)
  if(allocated(ie_indx))    deallocate(ie_indx)
  if(allocated(ix_indx))    deallocate(ix_indx)
  if(allocated(it_indx))    deallocate(it_indx)
  if(allocated(thcyc))      deallocate(thcyc)
  if(allocated(driftx))     deallocate(driftx)
  if(allocated(driftxrot1)) deallocate(driftxrot1)
  if(allocated(driftxrot2)) deallocate(driftxrot2)
  if(allocated(driftxrot3)) deallocate(driftxrot3)

contains

  subroutine set_RHS_source
    implicit none
    real :: src_F0_Ln, src_F0_Lt, src_P0, src_Rot1, src_Rot2

    ! First-Order Source Term
    g(:) = 0.0
    do is=1,n_species

       ! src = (1/F0) * dF0/dr + Ze/T dPhi0/dr
       src_F0_Ln = -(dlnndr(is,ir) - 1.5*dlntdr(is,ir))  ! ene^0 part
       src_F0_Lt = -dlntdr(is,ir)                        ! ene^1 part
       src_P0    = (1.0*Z(is))/temp(is,ir) * dphi0dr(ir) ! ene^0 Er part

       do ie=0,n_energy
          do ix=0,n_xi
             do it=1,n_theta

                i = mindx(is,ie,ix,it) 

                src_Rot1   = -omega_rot(ir) * bigR_th0**2 &
                     / vth(is,ir)**2 * omega_rot_deriv(ir) &
                     - dlntdr(is,ir) * Z(is) / temp(is,ir) * phi_rot(it) &
                     + dlntdr(is,ir) * (omega_rot(ir)/vth(is,ir))**2 &
                     * 0.5 * (bigR(it)**2 - bigR_th0**2) &
                     - omega_rot(ir)**2 * bigR_th0 &
                     / vth(is,ir)**2 * bigR_th0_rderiv

                if(aniso_model(is) == 2) then
                   src_Rot1 = src_Rot1 - dlntdr(is,ir) * lam_rot_aniso(is,it)
                endif

                src_Rot2 = omega_rot_deriv(ir) * bigR(it) / vth(is,ir)

                ! Impose constant for g0

                if(ix==0 .and. it==1 .and. &
                     (collision_model==1 .or. collision_model==2)) then
                   g(i) = 0.0

                else if(ix==0 .and. it==1 &
                     .and. (ie==0 .or. ie==1) &
                     .and. (collision_model==3 .or. &
                     collision_model==4 .or. collision_model==5)) then
                   g(i) = 0.0

                else 

                   ! Equilibrium source:
                   ! -(vdrift dot grad (F0 + Ze/T Phi0))

                   if (ix == 0) then
                      g(i) = -(4.0/3.0) * driftx(is,it) &
                           * ( (src_F0_Ln + src_P0 + src_Rot1) &
                           * evec_e1(ie,ix) &
                           + src_F0_Lt * evec_e2(ie,ix) ) &
                           - driftxrot1(is,it) &
                           * ( (src_F0_Ln + src_P0 + src_Rot1) &
                           * evec_e0(ie,ix) &
                           + src_F0_Lt * evec_e1(ie,ix) ) &
                           - src_Rot2 * 1.0/3.0* driftxrot2(is,it) &
                           * sqrt(2.0) * Btor(it)/Bmag(it) &
                           * evec_e1(ie,ix) &
                           - src_Rot2 * 4.0/3.0 * driftx(is,it) &
                           * omega_rot(ir) * bigR(it)/vth(is,ir) &
                           * evec_e1(ie,ix) &
                           - src_Rot2 * driftxrot1(is,it) &
                           * omega_rot(ir) * bigR(it)/vth(is,ir) &
                           * evec_e0(ie,ix)

                   else if(ix == 1) then
                      g(i) = sqrt(2.0) * vth(is,ir) * evec_e05(ie,ix) &
                           * (1.0*Z(is))/temp(is,ir) * epar0(ir) * Bmag(it) &
                           / Bmag2_avg &
                           - driftxrot2(is,it) &
                           * ( (src_F0_Ln + src_P0 + src_Rot1) &
                           * evec_e05(ie,ix) &
                           + src_F0_Lt * evec_e105(ie,ix) ) &
                           - src_Rot2 * 8.0/5.0 * driftx(is,it) &
                           * sqrt(2.0) * Btor(it)/Bmag(it) &
                           * evec_e105(ie,ix) &
                           - src_Rot2 * driftxrot1(is,it) &
                           * sqrt(2.0) * Btor(it)/Bmag(it) &
                           * evec_e05(ie,ix) &
                           - src_Rot2 * driftxrot2(is,it) &
                           * omega_rot(ir) * bigR(it)/vth(is,ir) &
                           * evec_e05(ie,ix) &
                           - src_Rot2 * 2.0/5.0* driftxrot3(is,it) &
                           * evec_e105(ie,ix)

                   else if (ix == 2) then
                      g(i) = -(2.0/3.0) * driftx(is,it) &
                           * ( (src_F0_Ln + src_P0 + src_Rot1) &
                           * evec_e1(ie,ix) &
                           + src_F0_Lt * evec_e2(ie,ix) ) &
                           - src_Rot2 * 2.0/3.0* driftxrot2(is,it) &
                           * sqrt(2.0) * Btor(it)/Bmag(it) &
                           * evec_e1(ie,ix) &
                           - src_Rot2 * 2.0/3.0 * driftx(is,it) &
                           * omega_rot(ir) * bigR(it)/vth(is,ir) &
                           * evec_e1(ie,ix)

                   else if(ix == 3) then
                      g(i) = -src_Rot2 * 2.0/5.0 * driftx(is,it) &
                           * sqrt(2.0) * Btor(it)/Bmag(it) &
                           * evec_e105(ie,ix) &
                           + src_Rot2 * 2.0/5.0* driftxrot3(is,it) &
                           * evec_e105(ie,ix)

                   endif


                endif

             enddo
          enddo
       enddo
    enddo

  end subroutine set_RHS_source

end subroutine neo_do
