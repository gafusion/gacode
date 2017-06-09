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

  integer :: ie,ix,is,ir,it
  integer :: d
  integer :: splitkey
  integer, external :: parallel_dim

  integer, external :: omp_get_max_threads, omp_get_thread_num

  ! Velocity-space (v) and configuration-space (c) dimensions
  nv = n_energy*n_xi*n_species
  nc = n_radial*n_theta

  call gcd(nv,nc,d)

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
     write(io,*)
     write(io,*) '          [coll]     [str]      [NL]      [NL]      [NL]'
     write(io,*) ' n_MPI    nc_loc    nv_loc   n_split  atoa[MB] atoa proc'
     write(io,*) '------    ------    ------   -------  -------- ---------'
     do it=1,d*n_toroidal
        if (mod(d*n_toroidal,it) == 0 .and. mod(it,n_toroidal) == 0) then
           n_proc_1 = it/n_toroidal
           nc_loc = nc/n_proc_1           
           nv_loc = nv/n_proc_1           
           nsplit = 1+(nv_loc*n_theta-1)/n_toroidal
           write(io,'(t2,4(i6,4x),f5.2,4x,i6)') &
                it,nc_loc,nv_loc,nsplit,16.0*n_radial*nsplit/1e6,n_toroidal
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
!$acc enter data copyin(ie_v,ix_v,is_v,iv_v)

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
!$acc enter data copyin(ir_c,it_c,ic_c)

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
     return
  endif

  !-------------------------------------------------------------
  ! Check that n_proc is a multiple of n_toroidal
  !
  if (modulo(n_proc,n_toroidal) /= 0) then
     call cgyro_error('Number of processors must be a multiple of N_TOROIDAL.')
     return
  endif

  ! Assign subgroup dimensions: n_proc = n_proc_1 * n_proc_2

  n_proc_1 = n_proc/n_toroidal
  n_proc_2 = n_toroidal

  ! Check that nv and nc are multiples of the local processor count

  if (modulo(nv,n_proc_1) /= 0 .or. modulo(nc,n_proc_1) /= 0) then
     call cgyro_error('nv or nc not a multiple of the local processor count.')
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
     call cgyro_error('ERROR: (CGYRO) NEW_COMM_1 not created')
     return
  endif

  ! Local adjoint Group number

  call MPI_COMM_SPLIT(CGYRO_COMM_WORLD,&
       i_group_2,& 
       splitkey,&
       NEW_COMM_2, &
       i_err)
  if (i_err /= 0) then
     call cgyro_error('ERROR: (CGYRO) NEW_COMM_2 not created')
     return
  endif
  !
  call MPI_COMM_RANK(NEW_COMM_1,i_proc_1,i_err)
  call MPI_COMM_RANK(NEW_COMM_2,i_proc_2,i_err)
  !
  !-----------------------------------------------------------

  ! Linear parallelization dimensions

  ! ni -> nc
  ! nj -> nv  
  call parallel_lib_init(nc,nv,nc_loc,nv_loc,NEW_COMM_1)

  nv1 = 1+i_proc_1*nv_loc
  nv2 = (1+i_proc_1)*nv_loc

  nc1 = 1+i_proc_1*nc_loc
  nc2 = (1+i_proc_1)*nc_loc

  ! Nonlinear parallelization dimensions (returns nsplit)

  if (nonlinear_method == 1) then
     call parallel_slib_init(n_toroidal,nv_loc,nc,nsplit,NEW_COMM_2)
  else
     call parallel_slib_init(n_toroidal,nv_loc*n_theta,n_radial,nsplit,NEW_COMM_2)
  endif

  ! OMP code
  n_omp = omp_get_max_threads()

  ! Restart file chunking logic (max_filesize is 8.0 GB)
  max_filesize = 8.0e9
  n_chunk = 1+int(16.0*n_toroidal*nc*nv/max_filesize)

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
