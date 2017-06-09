subroutine neo_spitzer

  ! Spitzer problem: 
  ! C_ee + nu_ei * L g_e 
  ! = -v_par*(Z_e e E_par/T_e - grad_par ln p_e 
  !           + (ene-5/2) grad_par ln T_e

  use neo_globals
  use neo_energy_grid
  use neo_sparse_solve
  implicit none
  integer :: io_sp=60
  integer :: is_ele, is_ion
  integer :: ie, is, ir, ix, je, i, j, k, n_elem, ierr
  real :: L0, L11, L12, L21, L22
  real :: nu_ei
  real, dimension(3) :: sp_pflux, sp_eflux, src1, src2

  is_ele = -1
  do is=1, n_species
     if(Z(is) < 0.0) then
        is_ele = is
        exit
     endif
  enddo
  if(is_ele == -1) then
     call neo_error('ERROR: (NEO) Must have electron species for Spitzer problem')
     return
  endif
  is_ion = -1
  do is=1, n_species
     if(Z(is) > 0) then
        is_ion = is
        exit
     endif
  enddo
  if(is_ele == -1) then
     call neo_error('ERROR: (NEO) Must have ion species for Spitzer problem')
     return
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  is = is_ele
  ix = 1
  ir = 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Matrix set-up
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call ENERGY_basis_ints_alloc(1)
  call ENERGY_basis_ints
  call ENERGY_coll_ints_alloc(1)
  call ENERGY_coll_ints(ir)

  n_row = n_energy+1
  n_max = (n_energy+1)**2 * 10
  allocate(amat(n_max),stat=ierr)
  if(ierr /= 0) then
     call neo_error('ERROR: (NEO) Spitzer allocation failed')
     goto 100
  end if
  allocate(amat_indx(2*n_max),stat=ierr)
  if(ierr /= 0) then
     call neo_error('ERROR: (NEO) Spitzer allocation failed')
     goto 100
  end if
  allocate(g(n_row))

  amat_indx(:) = 0
  amat(:) = 0.0
  k = 0
  do ie=0,n_energy
     i = ie+1
     do je=0,n_energy
        j = je+1
        k = k+1
        amat(k) = -(emat_coll_test(is,is,ie,je,ix) &
             + emat_coll_field(is,is,ie,je,ix) &
             + emat_coll_test(is,is_ion,ie,je,ix))
        amat_indx(k) = i
        amat_indx(k+n_max) = j
     enddo
  enddo
  n_elem = k
  do k=1,n_elem
     amat_indx(n_elem+k) = amat_indx(n_max+k)
  enddo

  if(silent_flag == 0 .and. i_proc == 0) then
     open(unit=io_neoout,file=trim(path)//runfile_neoout,&
          status='old',position='append')
     write(io_neoout,*) 'Begin matrix factor'
     close(io_neoout)
  endif
  call SOLVE_factor(n_elem)
  if(error_status > 0) return
  if(silent_flag == 0 .and. i_proc == 0) then
     open(unit=io_neoout,file=trim(path)//runfile_neoout,&
          status='old',position='append')
     write(io_neoout,*) 'Done matrix factor'
     close(io_neoout)
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Matrix solve
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Source term for L11 and L21
  src1(1) =  1.0    
  src2(1) =  0.0
  ! Source term for L21 and L22
  src1(2) = 0.0
  src2(2) = 1.0 
  ! Source term for check
  src1(3) =  (1.0*Z(is))/temp(is,ir) * epar0_spitzer  &
       + dlnndr(is,ir) + dlntdr(is,ir)    
  src2(3) =  dlntdr(is,ir)  

  do j=1,3
     do ie=0,n_energy
        i = ie+1
        g(i) = sqrt(2.0) * vth(is,ir) &
             * ( (src1(j) - 2.5*src2(j)) * evec_e05(ie,ix) &
             + src2(j) * evec_e105(ie,ix) )
     enddo

     if(silent_flag == 0 .and. i_proc == 0) then
        open(unit=io_neoout,file=trim(path)//runfile_neoout,&
             status='old',position='append')
        write(io_neoout,*) 'Begin matrix solve'
        close(io_neoout)
     endif
     call SOLVE_do
     if(silent_flag == 0 .and. i_proc == 0) then
        open(unit=io_neoout,file=trim(path)//runfile_neoout,&
             status='old',position='append')
        write(io_neoout,*) 'Done matrix solve'
        close(io_neoout)
     endif

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Transport coefficients
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     sp_pflux(j) = 0.0
     sp_eflux(j) = 0.0
     do ie=0,n_energy
        i = ie + 1
        sp_pflux(j) = sp_pflux(j) &
             + dens(is,ir) &
             * sqrt(2.0) * vth(is,ir) &
             * (1.0/3.0) *  g(i) &
             * 4.0/sqrt(pi) * evec_e05(ie,ix)

        sp_eflux(j) = sp_eflux(j) &
             +  dens(is,ir) &
             * sqrt(2.0) * vth(is,ir) &
             * (1.0/3.0) * temp(is,ir) *  g(i) &
             * 4.0/sqrt(pi) * (evec_e105(ie,ix) - 2.5 * evec_e05(ie,ix))
     enddo
  enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Resistivity Parameters
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!

  nu_ei = nu(is,ir) * (1.0*Z(is_ion))**2 / (1.0*Z(is))**2 & 
       * dens(is_ion,ir)/dens(is,ir) * (4.0/3.0) / sqrt(pi)
  L0 = dens(is,ir) * temp(is,ir) / mass(is) / nu_ei

  L11 = sp_pflux(1) / src1(1) / L0
  L21 = sp_eflux(1) / src1(1) / L0
  L12 = sp_pflux(2) / src2(2) / L0
  L22 = sp_eflux(2) / src2(2) / L0

  if (silent_flag == 0 .and. i_proc == 0) then
     open(unit=io_sp,file=trim(path)//'out.neo.spitzer',status='replace')
     write (io_sp,'(e16.8)',advance='no') L11
     write (io_sp,'(e16.8)',advance='no') L12
     write (io_sp,'(e16.8)',advance='no') L21
     write (io_sp,'(e16.8)',advance='no') L22
     write (io_sp,'(e16.8)',advance='no') sp_pflux(3)
     write (io_sp,'(e16.8)',advance='no') (L11*src1(3) + L12*src2(3))*L0
     write (io_sp,'(e16.8)',advance='no') sp_eflux(3)
     write (io_sp,'(e16.8)',advance='no') (L21*src1(3) + L22*src2(3))*L0
     close(io_sp)
  endif

  ! Clean-up
100 continue
  call ENERGY_basis_ints_alloc(0)
  call ENERGY_coll_ints_alloc(0)
  if(allocated(amat))      deallocate(amat)
  if(allocated(amat_indx)) deallocate(amat_indx)
  if(allocated(g))      deallocate(g)

end subroutine neo_spitzer

