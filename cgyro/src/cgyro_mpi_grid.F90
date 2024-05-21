!-----------------------------------------------------------------
! cgyro_mpi_grid.f90
!
! PURPOSE:
!  Subroutinized main cgyro program.  
!-----------------------------------------------------------------

subroutine cgyro_mpi_grid

  use timer_lib
  use parallel_lib
  use mpi

  use cgyro_globals
  use cgyro_io

  implicit none

  integer :: ie,ix,is,ir,it,itm
  integer :: iltheta_min,iltheta_max
  integer :: d
  integer :: splitkey
  integer :: nproc_3
  character(len=192) :: msg

  integer, external :: omp_get_max_threads, omp_get_thread_num

  ! Velocity-space (v) and configuration-space (c) dimensions
  nv = n_energy*n_xi*n_species
  nc = n_radial*n_theta

  call gcd(nv,nc,d)

  if (modulo(n_toroidal,nt_loc) /= 0) then
     write (msg, "(A,I3,A,I2,A)") "N_TOROIDAL (",n_toroidal,") not a multiple of N_TOROIDAL_PER_PROCESS (",nt_loc,")"
     call cgyro_error(msg)
     return
  endif

  n_toroidal_procs = n_toroidal/nt_loc

  !-------------------------------------------------------------------------
  ! MPI diagnostics need to come early
  !
  if (silent_flag == 0 .and. i_proc == 0) then
     open(unit=io,file=trim(path)//runfile_mpi,status='replace')
     write(io,*) 'Parallelization and distribution diagnostics'
     write(io,*)
     write(io,'(a,i5)') '         nv: ',nv
     write(io,'(a,i5)') '         nc: ',nc
     write(io,'(a,i5)') ' GCD(nv,nc): ',d
     write(io,'(a,i5)') ' n_toroidal: ',n_toroidal
     write(io,'(a,i5)') '     nt_loc: ',nt_loc
     write(io,*)
     write(io,*) '          [coll]     [str]      [NL]      [NL]      [NL]    [coll]     [str]'
     write(io,*) ' n_MPI    nc_loc    nv_loc   n_split  atoa[MB] atoa proc atoa proc ared proc'
     write(io,*) '------    ------    ------   -------  -------- --------- --------- ---------'
     do it=1,d*n_toroidal_procs
        if (mod(d*n_toroidal_procs,it) == 0 .and. mod(it,n_toroidal_procs) == 0) then
           n_proc_1 = it/n_toroidal_procs
           ! further filter out incompatible multiples for velocity==2
           if ((velocity_order==1) .or. &
               (n_proc_1 == 1) .or. ( modulo(n_proc_1, n_species) == 0 ) ) then
                nc_loc = nc/n_proc_1           
                nv_loc = nv/n_proc_1           
                nsplit = 1+(nv_loc*n_theta-1)/n_toroidal_procs
                nproc_3 = n_proc_1
                if ((n_proc_1 /= 1) .and. (velocity_order==2)) nproc_3 = n_proc_1/3
                write(io,'(t2,4(i6,4x),f6.2,4x,i6,4x,i6,4x,i6)') &
                     it,nc_loc,nv_loc,nsplit,16.0*n_radial*nt_loc*nsplit/1e6,&
                     n_toroidal_procs,n_proc_1,nproc_3
           endif
        endif
     enddo
     close(io)
  endif
  !-------------------------------------------------------------------------

  allocate(ie_v(nv))
  allocate(ix_v(nv))
  allocate(is_v(nv))

  allocate(ir_c(nc))
  allocate(it_c(nc))

  allocate(ic_c(n_radial,n_theta))
  allocate(iv_v(n_energy,n_xi,n_species))

  ! Velocity pointers
  iv = 0
  if (velocity_order==1) then
    call cgyro_info('Velocity order 1')
    !original
    do ie=1,n_energy
      do ix=1,n_xi
        do is=1,n_species
           iv = iv+1
           ie_v(iv) = ie
           ix_v(iv) = ix
           is_v(iv) = is
           iv_v(ie,ix,is) = iv
        enddo
      enddo
    enddo
  else if (velocity_order==2) then
    call cgyro_info('Velocity order 2')
    ! optimized for minimizing species
    do is=1,n_species
      do ie=1,n_energy
        do ix=1,n_xi
           iv = iv+1
           ie_v(iv) = ie
           ix_v(iv) = ix
           is_v(iv) = is
           iv_v(ie,ix,is) = iv
        enddo
      enddo
    enddo
  else
     call cgyro_error('Unknown VELOCITY_ORDER.')
     return
  endif

#if defined(OMPGPU)
!$omp target enter data map(to:ie_v,ix_v,is_v,iv_v)
#elif defined(_OPENACC)
!$acc enter data copyin(ie_v,ix_v,is_v,iv_v)
#endif

  ! Configuration pointers
  ic = 0
  do ir=1,n_radial
     do it=1,n_theta
        ic = ic+1
        ir_c(ic) = ir
        it_c(ic) = it
        ic_c(ir,it) = ic
     enddo
  enddo

#if defined(OMPGPU)
!$omp target enter data map(to:ir_c,it_c,ic_c)
#elif defined(_OPENACC)
!$acc enter data copyin(ir_c,it_c,ic_c)
#endif

  ! Shear pointers
  allocate(ica_c(nc))
  allocate(icb_c(nc))
  ic = 0
  do ir=1,n_radial
     do it=1,n_theta
        ic = ic+1
        if (ir < n_radial) then
           ica_c(ic) = ic_c(ir,it)
           icb_c(ic) = ic_c(ir+1,it)
        else
           ica_c(ic) = ic_c(ir,it)
           icb_c(ic) = ic_c(1,it)
        endif
     enddo
  enddo

  if (test_flag == 1) then
     ! Set dimensions for calculation of memory in test mode
     nv_loc = nv
     nc_loc = nc
     nsplit = nv_loc*n_theta/n_toroidal
     nt_loc = 1
     nt1 = 0
     nt2 = 0
     return
  endif

  if (zf_test_mode > 0) then
     ! Zonal flow (n=0) test
     nt_loc = 1
  endif

  !-------------------------------------------------------------
  ! Check that n_proc is a multiple of n_toroidal_procs
  !
  if (modulo(n_proc,n_toroidal_procs) /= 0) then
     write (msg, "(A,I3,A,I2,A)") "MPI processes (",n_proc,&
                 ") not a multiple of N_TOROIDAL/N_TOROIDAL_PER_PROCESS (", &
                 n_toroidal_procs,")"
     call cgyro_error(msg)
     return
  endif

#if defined(OMPGPU)
!$omp target enter data map(to:nv,nc,n_toroidal,n_toroidal_procs)
!$omp target enter data map(to:n_radial,n_theta,n_field,n_energy,n_xi,n_species)
#elif defined(_OPENACC)
!$acc enter data copyin(nv,nc,n_toroidal,n_toroidal_procs)
!$acc enter data copyin(n_radial,n_theta,n_field,n_energy,n_xi,n_species)
#endif

  ! Assign subgroup dimensions: n_proc = n_proc_1 * n_proc_2

  n_proc_1 = n_proc/n_toroidal_procs
  n_proc_2 = n_toroidal_procs

  ! Check that nv and nc are multiples of toroidal MPI multiplier

  if (modulo(nv,n_proc_1) /= 0) then
     write (msg, "(A,I6,A,I3,A)") "nv (",nv,") not a multiple of coll atoa procs (",n_proc_1,")"
     call cgyro_error(msg)
     return
  endif

  if (modulo(nc,n_proc_1) /= 0) then
     write (msg, "(A,I6,A,I3,A)") "nc (",nc,") not a multiple of coll atoa procs (",n_proc_1,")"
     call cgyro_error(msg)
     return
  endif

  ! Local group indices:

  if (mpi_rank_order == 1) then
     i_group_1 = i_proc/n_proc_1
     i_group_2 = modulo(i_proc,n_proc_1)
     call cgyro_info('MPI rank alignment 1')
  else
     i_group_1 = modulo(i_proc,n_proc_2)
     i_group_2 = i_proc/n_proc_2
     call cgyro_info('MPI rank alignment 2')
  endif
  !------------------------------------------------

  !-----------------------------------------------------------
  ! Split up GYRO_COMM_WORLD into groups and adjoint:
  !
  !             NEW_COMM_1  and  NEW_COMM_2
  !
  splitkey = i_proc
  call MPI_COMM_SPLIT(CGYRO_COMM_WORLD,&
       i_group_1,& 
       splitkey,&
       NEW_COMM_1, &
       i_err)
  if (i_err /= 0) then
     call cgyro_error('NEW_COMM_1 not created')
     return
  endif

  ! Local adjoint Group number

  call MPI_COMM_SPLIT(CGYRO_COMM_WORLD,&
       i_group_2,& 
       splitkey,&
       NEW_COMM_2, &
       i_err)
  if (i_err /= 0) then
     call cgyro_error('NEW_COMM_2 not created')
     return
  endif
  !
  call MPI_COMM_RANK(NEW_COMM_1,i_proc_1,i_err)
  call MPI_COMM_RANK(NEW_COMM_2,i_proc_2,i_err)
  !-----------------------------------------------------------

  ! Define the toroidal range
  ! Note: The same test is in cgyro_make_progiles, too
  !       But we need my_toroidal early
  if (zf_test_mode > 0) then
     ! Zonal flow (n=0) test
     nt1 = 0
     nt2 = 0
  else if (n_toroidal == 1) then
     ! Single linear mode (assume n=1)
     nt1 = 1
     nt2 = 1
  else
     ! Multiple modes (n=0,1,2,...,n_toroidal-1)
     nt1 = i_group_1*nt_loc
     nt2 = nt1 + nt_loc -1
  endif

  ! Linear parallelization dimensions

  ! ni -> nc
  ! nj -> nv  
  call parallel_lib_init(nc,nv,nt1,nt_loc,nc_loc,nv_loc,NEW_COMM_1)

  nv1 = 1+i_proc_1*nv_loc
  nv2 = (1+i_proc_1)*nv_loc

  ns1 = 1
  ns2 = n_species
  i_group_3 = 1
  if (velocity_order==2) then
    ! Paricles are contiguous in this order
    ns1 = is_v(nv1)
    ns2 = is_v(nv2)
    ! We need a clean split, so that all ranks have the same number of species
    ! n_proc_1 == 1 is an exception, since we do not need to split anything
    if ( (n_proc_1 /= 1) .and. ( modulo(n_proc_1, n_species) /= 0 ) ) then
           write (msg, "(A,I3,A,I2,A)") "coll atoa procs (",n_proc_1,") not a multiple of n_species (",n_species,&
                   "), needed for VELOCITY_ORDER=2"
           call cgyro_error(msg)
           return
    endif
    i_group_3 = ns1
  endif

! when exchaning only specific species, we need a dedicated comm
  call MPI_COMM_SPLIT(NEW_COMM_1,&
       i_group_3,&
       splitkey,&
       NEW_COMM_3, &
       i_err)
  if (i_err /= 0) then
     call cgyro_error('NEW_COMM_3 not created')
     return
  endif
  !
  call MPI_COMM_RANK(NEW_COMM_3,i_proc_3,i_err)

  call parallel_clib_init(ns1,ns2,ns_loc,NEW_COMM_3)


  nc1 = 1+i_proc_1*nc_loc
  nc2 = (1+i_proc_1)*nc_loc

  ! Nonlinear parallelization dimensions (returns nsplit)

  call parallel_slib_init(n_toroidal_procs,nv_loc*n_theta,n_radial,nsplit,NEW_COMM_2)
  ! we will move the big fpack in two independent pieces
  nsplitA = (nsplit+1)/2
  nsplitB = nsplit-nsplitA

  n_jtheta = 0
  if (nonlinear_flag == 1) then
     do itm=1,n_toroidal_procs
       ! find max n_jtheta among all processes
       ! since we will need that to have equal number of rows
       ! in all the gpack buffers
       iltheta_min = 1+((itm-1)*nsplit)/nv_loc
       ! use min since we have uneven splitting and compute can go past limit
       iltheta_max = min(n_theta,1+(itm*nsplit-1)/nv_loc)
       n_jtheta = max(n_jtheta,iltheta_max-iltheta_min+1)
     enddo
  endif

#if defined(OMPGPU)
!$omp target enter data map(to:nt1,nt2,nt_loc,nv1,nv2,nv_loc,ns1,ns2,ns_loc)
!$omp target enter data map(to:nc1,nc2,nc_loc,n_jtheta,nsplit,nsplitA,nsplitB)
#elif defined(_OPENACC)
!$acc enter data copyin(nt1,nt2,nt_loc,nv1,nv2,nv_loc,ns1,ns2,ns_loc)
!$acc enter data copyin(nc1,nc2,nc_loc,n_jtheta,nsplit,nsplitA,nsplitB)
#endif

  fA_req_valid = .FALSE.
  fB_req_valid = .FALSE.
  g_req_valid = .FALSE.

  ! OMP code
  n_omp = omp_get_max_threads()

  !----------------------------------------------------------------------------
  ! Restart IO setup

  if (mpiio_stripe_factor > 0) then
     call cgyro_info("MPI IO factor set for restart file")
     write (mpiio_stripe_str,"(I3.3)") mpiio_stripe_factor
  endif

  if (mpiio_small_stripe_factor > 0) then
     call cgyro_info("MPI IO factor set for data files")
     write (mpiio_small_stripe_str,"(I2.2)") mpiio_small_stripe_factor
  endif

  ! save hostname configuration
  call cgyro_write_hosts
  !----------------------------------------------------------------------------

end subroutine cgyro_mpi_grid

subroutine gcd(m,n,d)

  implicit none

  integer, intent(in) :: m,n
  integer, intent(inout) :: d
  integer :: a,b,c

  a = m
  b = n

  if (a < b) then
     c = a
     a = b
     b = c
  endif

  do          
     c = mod(a, b)    
     if (c == 0) exit
     a = b         
     b = c 
  enddo

  d = b

end subroutine gcd
