module neo_3d_driver

  implicit none

  public :: ThreeD_do

  character(len=80),private  :: runfile_f = 'out.neo.f_3d'
  character(len=80),private  :: runfile_transp = 'out.neo.transport_3d'
  character(len=80),private :: runfile_check  = 'out.neo.prec'
  real, private :: check_sum
  integer, dimension(:,:,:,:), allocatable, private :: mindx 
  integer, private :: stat
  real, dimension(:,:), allocatable, private :: tpmat_trap, tpmat_stream_dt, &
       tpmat_stream_dp, tpmat_coll, tpmat_vexb_dt, tpmat_vexb_dp, &
       tpmat_vdrift_dt, tpmat_vdrift_dp
  real, dimension(:), allocatable, private :: tpvec_vdriftx, tpvec_flux, &
       tpvec_uparB, tpvec_upar, tpvec_fsa, tpvec_thetabar, tpvec_bmag, &
       tpvec_ntv
  ! transport coefficients
  real, dimension(:), allocatable, private :: pflux, eflux, uparB, ntv
  real, dimension(:,:), allocatable, private :: upar
  real, private :: jpar
  

contains

  subroutine ThreeD_do
    use neo_globals
    use neo_energy_grid
    use neo_allocate_profile
    implicit none
    integer :: ir, i, j
    integer :: status

    ! Set-up the energy grids for basis ints (indep of r)
    call ENERGY_basis_ints_alloc(1)
    call ENERGY_basis_ints
    call ENERGY_coll_ints_alloc(1)

    allocate(tpmat_trap(tpmatsize,tpmatsize))
    allocate(tpmat_stream_dt(tpmatsize,tpmatsize))
    allocate(tpmat_stream_dp(tpmatsize,tpmatsize))
    allocate(tpmat_coll(tpmatsize,tpmatsize))
    allocate(tpmat_vexb_dt(tpmatsize,tpmatsize))
    allocate(tpmat_vexb_dp(tpmatsize,tpmatsize))
    allocate(tpmat_vdrift_dt(tpmatsize,tpmatsize))
    allocate(tpmat_vdrift_dp(tpmatsize,tpmatsize))
    open(unit=1,file='out.le3.geomatrix',status='old',iostat=stat)
    if(stat .ne. 0) then
       call neo_error('ERROR: (NEO) le3 files not available')
       goto 200
    endif
    do i=1,tpmatsize
       do j=1, tpmatsize
          read(1,*) tpmat_trap(i,j), tpmat_stream_dt(i,j), &
               tpmat_stream_dp(i,j), tpmat_coll(i,j), &
               tpmat_vexb_dt(i,j),   tpmat_vexb_dp(i,j), &
               tpmat_vdrift_dt(i,j), tpmat_vdrift_dp(i,j)
       enddo
    enddo

    allocate(tpvec_thetabar(tpmatsize))
    allocate(tpvec_vdriftx(tpmatsize))
    allocate(tpvec_flux(tpmatsize))
    allocate(tpvec_uparB(tpmatsize))
    allocate(tpvec_upar(tpmatsize))
    allocate(tpvec_fsa(tpmatsize))
    allocate(tpvec_bmag(tpmatsize))
    allocate(tpvec_ntv(tpmatsize))
    open(unit=1,file='out.le3.geovector',status='old',iostat=stat)
    if(stat .ne. 0) then
       call neo_error('ERROR: (NEO) le3 files not available')
       goto 200
    endif
    do i=1,tpmatsize
       read(1,*) tpvec_thetabar(i), tpvec_vdriftx(i), tpvec_flux(i), &
            tpvec_uparB(i), tpvec_upar(i), tpvec_fsa(i), tpvec_bmag(i), &
            tpvec_ntv(i)
    enddo

    if (silent_flag == 0 .and. i_proc == 0) then
       open(unit=1,file=trim(path)//'out.neo.grid_3d',status='replace')
       write(1,*) n_species
       write(1,*) n_energy
       write(1,*) n_xi
       write(1,*) tpmatsize
       close(1)
       
       open(unit=1,file=trim(path)//runfile_f,status='replace')
       close(1)
    end if

    do ir=1,n_radial
       call threed_matrix_solve(ir,status)
       if(status == 1) goto 200
    enddo

200 call ENERGY_basis_ints_alloc(0)
    call ENERGY_coll_ints_alloc(0)
    call PROFILE_SIM_alloc(0)

    if (allocated(tpmat_trap))           deallocate(tpmat_trap)
    if (allocated(tpmat_stream_dt))      deallocate(tpmat_stream_dt)
    if (allocated(tpmat_stream_dp))      deallocate(tpmat_stream_dp)
    if (allocated(tpmat_coll))           deallocate(tpmat_coll)
    if (allocated(tpmat_vexb_dt))        deallocate(tpmat_vexb_dt)
    if (allocated(tpmat_vexb_dp))        deallocate(tpmat_vexb_dp)
    if (allocated(tpmat_vdrift_dt))      deallocate(tpmat_vdrift_dt)
    if (allocated(tpmat_vdrift_dp))      deallocate(tpmat_vdrift_dp)
    if (allocated(tpvec_thetabar))       deallocate(tpvec_thetabar) 
    if (allocated(tpvec_vdriftx))        deallocate(tpvec_vdriftx) 
    if (allocated(tpvec_flux))           deallocate(tpvec_flux)
    if (allocated(tpvec_uparB))          deallocate(tpvec_uparB) 
    if (allocated(tpvec_upar))           deallocate(tpvec_upar)
    if (allocated(tpvec_fsa))            deallocate(tpvec_fsa)
    if (allocated(tpvec_bmag))           deallocate(tpvec_bmag)
    if (allocated(tpvec_ntv))            deallocate(tpvec_ntv)

  end subroutine ThreeD_do

  subroutine threed_matrix_solve(ir,status)

    use neo_globals
    use neo_energy_grid
    use neo_g_velocitygrids

    implicit none
    integer, intent(in)  :: ir
    integer, intent(out) :: status
    integer :: is, ie, ix, it, js, je, jx, jt, ks, iab
    integer :: i,j,j_loc
    integer :: n_piv
    integer :: n_fill
    integer :: ierr
    integer :: info
    real :: fac

    status = 0

    ! Bandwidth = KL = KU = 2*nb-1
    bw  = 2*nb-1

    ! LAPACK Parameters
    
    npb = n_row
    
    n_fill = bw
    
    ! LDAB = 1+2*KL+KU 
    ldab = 1+2*bw+n_fill
    
    n_piv = n_row

    ! Here, npb means 'number of columns in a distributed block'
    allocate(ab(ldab,npb),stat=ierr)
    if(ierr /= 0) then
       call neo_error('ERROR: (NEO) Array allocation failed')
       status=1
       goto 100
    endif

    allocate(ipiv(n_piv),stat=ierr)
    if(ierr /= 0) then
       call neo_error('ERROR: (NEO) Array allocation failed')
       status=1
       goto 100
    endif

    allocate(g_loc(npb))
    allocate(g(n_row))
    allocate(mindx(n_species,0:n_energy,0:n_xi,tpmatsize))
    allocate(is_indx(n_row))
    allocate(ie_indx(n_row))
    allocate(ix_indx(n_row))
    allocate(it_indx(n_row))

    ! matrix indices
    i = 0
    do ix=0,n_xi
       do is=1,n_species
          do ie=0,n_energy
             do it=1,tpmatsize
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

    ! Set-up the matrix equation: LHS radially local matrix

    if(silent_flag == 0 .and. i_proc == 0) then
       open(unit=io_neoout,file=trim(path)//runfile_neoout,status='old',position='append')
       write(io_neoout,'(a,i3)') 'INFO: (NEO) ir =',ir
       write(io_neoout,'(a)') 'INFO: (NEO) Begin matrix setup'
       close(io_neoout)
    endif

    if(collision_model == 3) then
       do is=2, n_species
          if(abs(temp(is,ir)-temp(1,ir)) > 1e-3) then
             if(silent_flag == 0 .and. i_proc == 0) then
                open(unit=io_neoout,file=trim(path)//runfile_neoout,status='old',position='append')
                write(io_neoout,'(a)') 'WARNING: (NEO) Full HS collisions with unequal temps'
                close(io_neoout)
             endif
          endif
       enddo
    endif

    ! Get the energy grids for coll ints
    call ENERGY_coll_ints(ir)

    ! Set the LHS
    ab(:,:) = 0.0

    do i=1,n_row

       is = is_indx(i) 
       ie = ie_indx(i) 
       ix = ix_indx(i) 
       it = it_indx(i) 

       do j_loc=1,npb

          ! Compute global column index j from local index
          j = j_loc+i_proc*npb

          js = is_indx(j) 
          je = ie_indx(j) 
          jx = ix_indx(j) 
          jt = it_indx(j) 

          ! Band-storage index for a(i,j)
          iab = n_fill+bw+1+i-j

          ! Impose constant for g0
          ! local  -> at each species and ene 

          if (ix==0 .and. it == indx_c00 &
               .and. (collision_model==1 .or. collision_model==2))&
               then
             ! <f_ie> = 0
             if(js == is .and. je == ie .and. jx == 0) then
                ab(iab,j_loc) = tpvec_fsa(jt)
             endif

          else if(ix==0 .and. it == indx_c00 &
               .and. (ie==0 .or. ie==1) &
               .and. (collision_model==3 .or. &
               collision_model==4 .or. collision_model==5)) then
             ! <f_ie> = 0
             if(js == is .and. je == ie .and. jx == 0) then
                ab(iab,j_loc) = tpvec_fsa(jt)
             endif

          else

             ! Collisions 
             if(jx == ix) then

                ! test particle
                if(js == is) then
                   if(abs(tpmat_coll(it,jt)) > 1e-12) then
                      ab(iab,j_loc) = 0.0
                      do ks=1,n_species
                         ab(iab,j_loc) = ab(iab,j_loc) &
                              - emat_coll_test(is,ks,ie,je,ix) &
                              * tpmat_coll(it,jt)
                      enddo
                   endif
                endif

                ! field particle 
                if(abs(tpmat_coll(it,jt)) > 1e-12) then
                   ab(iab,j_loc) = ab(iab,j_loc)-emat_coll_field(is,js,ie,je,ix) &
                        * tpmat_coll(it,jt)
                endif
             endif

             ! hypervisc term
             if(is==js .and. ix==jx) then
                if(abs(tpmat_coll(it,jt)) > 1e-12) then
                   ab(iab,j_loc) = ab(iab,j_loc) &
                        + threed_hyperxi*(1.0*ix/(1.0*n_xi))**n_xi &
                        * tpmat_coll(it,jt) * emat_e0(ie,je,ix,1)
                endif
             endif

             ! Streaming -- d/dtheta and d/dphi
             if(js == is) then
                fac = sqrt(2.0) * vth(is,ir)
                if(jx == ix-1 .and. jx >= 0) then
                   ab(iab,j_loc) = ab(iab,j_loc) + fac &
                        * ix/(2*ix-1.0) &
                        * emat_e05(ie,je,ix,1) &
                        * (tpmat_stream_dt(it,jt) &
                        + tpmat_stream_dp(it,jt))  
                endif
                if(jx == ix+1 .and. jx <= n_xi) then
                   ab(iab,j_loc) = ab(iab,j_loc) + fac &
                        * (ix+1.0)/(2*ix+3.0) &
                        * emat_e05(ie,je,ix,2) &
                        * (tpmat_stream_dt(it,jt) &
                        + tpmat_stream_dp(it,jt))
                endif
             endif

             ! Trapping 
             if(js == is) then
                fac  = sqrt(0.5) * vth(is,ir)
                if(jx == ix-1 .and. jx >= 0) then
                   ab(iab,j_loc) = ab(iab,j_loc) + fac &
                        * ix*(ix-1.0)/(2*ix-1.0) &
                        * tpmat_trap(it,jt) &
                        * emat_e05(ie,je,ix,1) 
                endif
                if(jx == ix+1 .and. jx <= n_xi) then
                   ab(iab,j_loc) = ab(iab,j_loc) - fac &
                        * (ix+1.0)*(ix+2.0)/(2*ix+3.0) &
                        * tpmat_trap(it,jt) &
                        * emat_e05(ie,je,ix,2)
                endif
             endif

             ! V_ExB term
             if(threed_exb_model == 1) then
                if(js == is .and. jx == ix) then
                   ab(iab,j_loc) = ab(iab,j_loc) &
                        + rho(ir)*threed_exb_dphi0dr &
                        * emat_e0(ie,je,ix,1) &
                        * (tpmat_vexb_dt(it,jt) &
                        + tpmat_vexb_dp(it,jt))
                endif
             endif

             if(threed_drift_model == 1) then
                if(js == is) then
                   fac = rho(ir)*mass(is)/(1.0*Z(is)) * (vth(is,ir))**2
                   if(jx == ix-2 .and. jx >= 0) then
                      ab(iab,j_loc) = ab(iab,j_loc) + fac &
                           * (ix)*(ix-1.0)/((2*ix-3.0)*(2*ix-1.0)) &
                           * (0.0*tpmat_vdrift_dt(it,jt) + tpmat_vdrift_dp(it,jt))&
                           * emat_e1(ie,je,ix,1)
                   endif
                   if(jx == ix) then
                      ab(iab,j_loc) = ab(iab,j_loc) + fac &
                           * (1.0 + (ix+1.0)**2/((2*ix+1.0)*(2*ix+3.0)) &
                           + ix**2/((2*ix+1.0)*(2*ix-1.0))) &
                           * (0.0*tpmat_vdrift_dt(it,jt) + tpmat_vdrift_dp(it,jt))&
                           * emat_e1(ie,je,ix,2)
                   endif
                   if(jx == ix+2 .and. jx <= n_xi) then
                      ab(iab,j_loc) = ab(iab,j_loc) + fac &
                           * (ix+1.0)*(ix+2.0)/((2*ix+3.0)*(2*ix+5.0)) &
                           * (0.0*tpmat_vdrift_dt(it,jt) + tpmat_vdrift_dp(it,jt))&
                           * emat_e1(ie,je,ix,3)
                   endif
                endif
             endif

          endif ! bc if/else
       enddo
    enddo

    ! Set the RHS source 
    call threed_set_source(ir)

    ! Factor the Matrix and Solve -- uses a(:) and a_indx(:)

    if(silent_flag == 0 .and. i_proc == 0) then
       open(unit=io_neoout,file=trim(path)//runfile_neoout,&
            status='old',position='append')
       write(io_neoout,'(a)') 'INFO: (NEO) Begin matrix factor and solve'
       close(io_neoout)
    endif

    ! Factor and solve with LAPACK band-matrix routine
    g(:) = g_loc(:)
    call DGBSV(n_row,2*nb-1,2*nb-1,1,ab,ldab,ipiv,g,n_row,info)

    if(info /= 0) then
       call neo_error('ERROR: (NEO) Matrix factorization failed')
       status=1
       goto 100
    endif

    if(silent_flag == 0 .and. i_proc == 0) then
       open(unit=io_neoout,file=trim(path)//runfile_neoout,&
            status='old',position='append')
       write(io_neoout,'(a)') 'INFO: (NEO) Done matrix factor and solve'
       close(io_neoout)
    endif

    ! Compute the neo transport coefficients
    call threed_transport(ir)

    if(silent_flag == 0 .and. i_proc == 0) then
       open(1,file=trim(path)//runfile_f,status='old',position='append')
       write(1,'(1pe12.5)') g(:)
       close(1)
    endif

    call g_xi


100 continue
    if(allocated(ab))              deallocate(ab)
    if(allocated(ipiv))            deallocate(ipiv)
    if(allocated(g))               deallocate(g)
    if(allocated(g_loc))           deallocate(g_loc)
    if(allocated(mindx))           deallocate(mindx)
    if(allocated(is_indx))         deallocate(is_indx)
    if(allocated(ie_indx))         deallocate(ie_indx)
    if(allocated(ix_indx))         deallocate(ix_indx)
    if(allocated(it_indx))         deallocate(it_indx)

  end subroutine threed_matrix_solve

  subroutine threed_set_source(ir)

    use neo_globals
    use neo_energy_grid

    implicit none

    integer, intent(in) :: ir
    real :: src_F0_Ln, src_F0_Lt, src_P0
    integer :: i,i_loc
    integer :: is, ie, ix, it
    real :: fac


    ! First-Order Source Term
    g_loc(:) = 0.0

    do i_loc=1,npb

       ! Compute global column index i from local index
       i = i_loc+i_proc*npb

       is = is_indx(i) 
       ie = ie_indx(i) 
       ix = ix_indx(i) 
       it = it_indx(i) 

       fac = rho(ir)*mass(is)/(1.0*Z(is)) * (vth(is,ir))**2

       ! src = (1/F0) * dF0/dr + Ze/T dPhi0/dr
       src_F0_Ln = -(dlnndr(is,ir) - 1.5*dlntdr(is,ir))  ! ene^0 part
       src_F0_Lt = -dlntdr(is,ir)                        ! ene^1 part
       src_P0    = (1.0*Z(is))/temp(is,ir) * dphi0dr(ir) ! ene^0 Er part

       ! Impose constant for g0

       if(ix==0 .and. it== indx_c00.and. &
            (collision_model==1 .or. collision_model==2)) then
          g_loc(i_loc) = 0.0

       else if(ix==0 .and. it== indx_c00 &
            .and. (ie==0 .or. ie==1) &
            .and. (collision_model==3 .or. &
            collision_model==4 .or. collision_model==5)) then
          g_loc(i_loc) = 0.0

       else 

          ! Equilibrium source:
          ! -(vdrift dot grad (F0 + Ze/T Phi0))

          if (ix == 0) then
             g_loc(i_loc) = -(4.0/3.0) * fac * tpvec_vdriftx(it) &
                  * ( (src_F0_Ln + src_P0) &
                  * evec_e1(ie,ix) &
                  + src_F0_Lt * evec_e2(ie,ix) )

          else if (ix == 2) then
             g_loc(i_loc) = -(2.0/3.0) * fac * tpvec_vdriftx(it) &
                  * ( (src_F0_Ln + src_P0) &
                  * evec_e1(ie,ix) &
                  + src_F0_Lt * evec_e2(ie,ix) )

          else if (ix == 1) then
             g_loc(i_loc) = sqrt(2.0) * vth(is,ir) * evec_e05(ie,ix) &
                  * (1.0*Z(is))/temp(is,ir) * epar0(ir) * tpvec_bmag(it) &
                  / threed_bmag2_avg

          endif

       endif

    enddo

  end subroutine threed_set_source

  subroutine threed_transport(ir)
    use neo_globals
    use neo_energy_grid
    implicit none
    integer, intent(in) :: ir
    integer :: is, ie, ix, it, i
    real    :: fac_drift, fac_upar

    allocate(pflux(n_species))
    allocate(eflux(n_species))
    allocate(uparB(n_species))
    allocate(upar(n_species,tpmatsize))
    allocate(ntv(n_species))

    check_sum=0.0

    pflux(:)  = 0.0
    eflux(:)  = 0.0
    uparB(:)  = 0.0
    upar(:,:) = 0.0
    ntv(:)    = 0.0

    ! Compute the transport coefficients here

    do i=1,n_row

       is = is_indx(i)
       ie = ie_indx(i)
       ix = ix_indx(i)
       it = it_indx(i)

       fac_drift = rho(ir)*mass(is)/(1.0*Z(is)) * (vth(is,ir))**2

       fac_upar  = sqrt(2.0) * vth(is,ir)

       if(ix == 0) then

          pflux(is) = pflux(is) + g(i) &
               * 4.0/sqrt(pi) * dens(is,ir) &
               * (4.0/3.0) * evec_e1(ie,ix) &
               * fac_drift * tpvec_flux(it)

          eflux(is) = eflux(is) + g(i) &
               * 4.0/sqrt(pi) * dens(is,ir) &
               * (4.0/3.0) * evec_e2(ie,ix) * temp(is,ir) &
               * fac_drift * tpvec_flux(it)

          ntv(is) = ntv(is) + g(i) &
               * 4.0/sqrt(pi) * dens(is,ir) &
               * (4.0/3.0) * evec_e1(ie,ix) * temp(is,ir) &
               * (-0.5) * tpvec_ntv(it)

       else if(ix == 2) then
          
          pflux(is) = pflux(is) + g(i) &
               * 4.0/sqrt(pi) * dens(is,ir) &
               * (2.0/15.0) * evec_e1(ie,ix) &
               * fac_drift * tpvec_flux(it)

          eflux(is) = eflux(is) + g(i) &
               * 4.0/sqrt(pi) * dens(is,ir) &
               * (2.0/15.0) * evec_e2(ie,ix) * temp(is,ir) &
               * fac_drift * tpvec_flux(it)

          ntv(is) = ntv(is) + g(i) &
               * 4.0/sqrt(pi) * dens(is,ir) &
               * (2.0/15.0) * evec_e1(ie,ix) * temp(is,ir) &
               * (-0.5) * tpvec_ntv(it)

       else if(ix == 1) then

          upar(is,it) = upar(is,it) + g(i) &
               * 4.0/sqrt(pi) * (1.0/3.0) * evec_e05(ie,ix) &
               * fac_upar * tpvec_upar(it)
          
          uparB(is) = uparB(is) + g(i) &
               * 4.0/sqrt(pi) * (1.0/3.0) * evec_e05(ie,ix) &
               * fac_upar * tpvec_uparB(it)

       endif

    enddo

    ! Bootstrap current = sum <Z*n*upar B>
    jpar = 0.0
    do is=1, n_species
       jpar = jpar + Z(is) * uparB(is) * dens(is,ir)
    enddo

    if(silent_flag == 0 .and. i_proc == 0) then
       if(ir==1) then
          open(unit=1,file=trim(path)//runfile_transp,status='replace')
       else
          open(1,file=trim(path)//runfile_transp,status='old',position='append')
       endif
       write (1,'(e16.8)',advance='no') r(ir)
       write (1,'(e16.8)',advance='no') jpar
       do is=1, n_species
          write (1,'(e16.8)',advance='no') pflux(is)
          write (1,'(e16.8)',advance='no') eflux(is)
          write (1,'(e16.8)',advance='no') uparB(is)
          write (1,'(e16.8)',advance='no') ntv(is)
       enddo
       write (1,*)
       close(1)
    endif

    do is=1,n_species
       check_sum = check_sum &
            +  (abs(pflux(is)) + abs(eflux(is)) + abs(ntv(is)))/rho(ir)**2 &
            + abs(uparB(is))/rho(ir)
    enddo
    if(silent_flag == 0 .and. i_proc == 0) then
       open(unit=1,file=trim(path)//runfile_check,status='replace')
       write (1,'(e16.8)',advance='no') check_sum
       close(1)
    endif

    deallocate(pflux)
    deallocate(eflux)
    deallocate(uparB)
    deallocate(upar)
    deallocate(ntv)

  end subroutine threed_transport

end module neo_3d_driver
