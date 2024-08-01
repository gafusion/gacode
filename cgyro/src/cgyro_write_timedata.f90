!-----------------------------------------------------------------
! cgyro_write_timedata.f90
!
! PURPOSE:
!  Output of time-dependent data.
!-----------------------------------------------------------------

subroutine cgyro_write_timedata

  use mpi
  use cgyro_globals
  use cgyro_step

  implicit none

  logical :: has_zf, has_balloon
  integer :: i_field,i_moment
  integer :: ir,it
  integer :: p_field
  real :: fvec(2,nt1:nt2)
  real :: vec(4)
  complex :: a_norm
  complex :: ftemp(n_theta,n_radial)
  complex :: field_plot(n_radial,theta_plot,nt1:nt2)

  ! Print this data on print steps only; otherwise exit now
  if (mod(i_time,print_step) /= 0) return

  ! Increment the print counter on actual output steps
  if (io_control == 2) i_current = i_current+1

  !---------------------------------------------------------------------------
  if (n_toroidal == 1 .and. h_print_flag == 1) then
     call write_distribution(trim(path)//binfile_hb)
  endif
  !---------------------------------------------------------------------------

  call cgyro_flux

  ! ky flux for all species with field breakdown
  call cgyro_write_distributed_breal(&
       trim(path)//binfile_ky_flux,&
       size(gflux(0,:,1:nflux,:,:)),&
       real(gflux(0,:,1:nflux,:,:)))

  ! central ky flux for all species with field breakdown
  call cgyro_write_distributed_breal(&
       trim(path)//binfile_ky_cflux,&
       size(cflux(:,1:nflux,:,:)),&
       cflux(:,1:nflux,:,:))

  if (gflux_print_flag == 1) then
     ! Global (n,e,v) fluxes for all species
     do i_moment=1,3
        call cgyro_write_distributed_bcomplex(&
             trim(path)//binfile_lky_flux(i_moment),&
             size(gflux(:,:,i_moment,:,:)),&
             gflux(:,:,i_moment,:,:))
     enddo
  endif

  if (moment_print_flag == 1) then
     ! (n,e,v) moment for all species at selected thetas.
     ! prattq, allow moment printing for linear runs.
     do i_moment=1,3
        call cgyro_write_distributed_bcomplex(&
             trim(path)//binfile_kxky(i_moment),&
             size(moment(:,:,:,:,i_moment)),&
             moment(:,:,:,:,i_moment))
     enddo
  endif

  if (field_print_flag == 1) then
     p_field = n_field
  else
     p_field = 1
  endif

  do i_field=1,p_field
     ! Sort out subset of theta values for plotting
     do ic=1,nc
        ir = ir_c(ic)
        it = it_c(ic)
        if (itp(it) > 0) then
           field_plot(ir,itp(it),nt1:nt2) = field(i_field,ic,nt1:nt2)
        endif
     enddo

     ! Complex potentials at selected thetas
     call cgyro_write_distributed_bcomplex(&
          trim(path)//binfile_kxky_field(i_field),&
          size(field_plot),&
          field_plot)
  enddo

  ! Checksum for regression testing
  ! Note that checksum is a distributed real scalar
  if (zf_test_mode == 0) then
     ! Do not include exchange in precision
     call write_precision(trim(path)//runfile_prec,sum(abs(real(gflux(0,:,1:3,:,:)))))
  else
     call write_precision(trim(path)//runfile_prec,sum(abs(field)))
  endif

  !------------------------------------------------------------------
  ! Ballooning mode (or ZF) output for linear runs with a single mode
  ! (can both be plotted with cgyro_plot -plot ball)
  !
  has_balloon = (n_toroidal == 1) .and. ((nt1 > 0)  .and. (box_size == 1))
  has_zf      = zf_test_mode > 0
  if ( (i_proc==0) .and. (has_zf .or. has_balloon) ) then
     ! NOTE: Only process the first my_toroidal

     do i_field=1,n_field

        do ir=1,n_radial
           do it=1,n_theta
              ftemp(it,ir) = field(i_field,ic_c(ir,it),nt1)
           enddo
        enddo

        if (has_balloon) then
           if (i_field == 1) then
              it = maxloc(abs(ftemp(:,n_radial/2+1)),dim=1)
              a_norm = ftemp(it,n_radial/2+1)
           endif
           if (abs(a_norm) .lt. 1e-20) a_norm = 1.0
           call extended_ang(ftemp)
        else
           a_norm = 1.0
        endif
        if (abs(a_norm) .le. 1e-20) a_norm = 1.0
        call write_binary(trim(path)//binfile_fieldb(i_field),&
             ftemp(:,:)/a_norm,size(ftemp))
     enddo
  endif
  !---------------------------------------------------------------

  ! Linear frequency diagnostics for every value of n
  call cgyro_freq
  fvec(1,:) = real(freq(:)) ; fvec(2,:) = aimag(freq(:))
  if (n_toroidal > 1) then
     ! NOTE: Update counter when having more than one my_toroidal
     call cgyro_write_distributed_breal(trim(path)//binfile_freq,size(fvec),fvec(:,:))
  else 
     call write_ascii(trim(path)//runfile_freq,2,fvec(:,:))
  endif
  
  ! Output to screen
  if (printout) call print_scrdata()

  ! Output to files
  vec(1) = t_current ; vec(2:3) = integration_error(:) ; vec(4) = delta_t_gk
  call write_ascii(trim(path)//runfile_time,4,vec(1:4))

  call MPI_BCAST(signal,1,MPI_INTEGER,0,CGYRO_COMM_WORLD,i_err)

end subroutine cgyro_write_timedata

!===============================================================================
! Individual parallel-safe I/O routines
!===============================================================================

!------------------------------------------------------
! cgyro_write_distributed_bcomplex.f90
!
! PURPOSE:
!  Control merged output of complex distributed array.
!------------------------------------------------------

subroutine cgyro_write_distributed_bcomplex(datafile,n_fn,fn)

  use mpi
  use cgyro_globals
  use cgyro_io, only: cgyro_error

  !------------------------------------------------------
  implicit none
  !
  character (len=*), intent(in) :: datafile
  integer, intent(in) :: n_fn
  complex, intent(in) :: fn(n_fn)
  !
  ! Required for MPI-IO:
  !
  integer :: filemode
  integer :: finfo
  integer :: fh
  integer :: fstatus(MPI_STATUS_SIZE)
  integer(kind=MPI_OFFSET_KIND) :: disp
  integer(kind=MPI_OFFSET_KIND) :: offset1
  character :: cdummy
  !------------------------------------------------------

  if (i_proc_1 /= 0) return

  select case (io_control)

  case (0)

     return

  case (1)

     ! Open

     if (i_proc == 0) then
        open(unit=io,file=datafile,status='replace')
        close(io)
     endif

  case (2)

     ! Append

     ! Write in parallel to the binary datafile
     filemode = MPI_MODE_WRONLY
     disp = i_current-1
     disp = disp*n_proc_2*size(fn)*BYTE*2

     offset1 = i_proc_2*size(fn)

     call MPI_INFO_CREATE(finfo,i_err)

     if (mpiio_small_stripe_factor > 0) then
        ! user asked us to explicitly set the MPI IO striping factor
        call MPI_INFO_SET(finfo,"striping_factor",mpiio_small_stripe_str,i_err)
     endif

     call MPI_FILE_OPEN(NEW_COMM_2,&
          datafile,&
          filemode,&
          finfo,&
          fh,&
          i_err)

     if (BYTE == 4) then

        ! Single (default) 
        call MPI_FILE_SET_VIEW(fh,&
             disp,&
             MPI_COMPLEX8,&
             MPI_COMPLEX8,&
             'native',&
             finfo,&
             i_err)

        call MPI_FILE_WRITE_AT(fh,&
             offset1,&
             cmplx(fn,kind=4),&
             n_fn,&
             MPI_COMPLEX8,&
             fstatus,&
             i_err)
     else

        call MPI_FILE_SET_VIEW(fh,&
             disp,&
             MPI_COMPLEX16,&
             MPI_COMPLEX16,&
             'native',&
             finfo,&
             i_err)

        call MPI_FILE_WRITE_AT(fh,&
             offset1,&
             fn,&
             n_fn,&
             MPI_COMPLEX16,&
             fstatus,&
             i_err)
     endif

     call MPI_FILE_SYNC(fh,i_err)
     call MPI_FILE_CLOSE(fh,i_err)
     call MPI_INFO_FREE(finfo,i_err)

  case (3)

     ! Rewind

     if (i_proc == 0) then

        disp = i_current
        disp = disp*n_proc_2*size(fn)*BYTE*2

        open(unit=io,file=datafile,status='old',access='stream',iostat=i_err)
        if (i_err/=0) then
           call cgyro_error('[REWIND] Failed to open '//datafile)
           return
        endif
        if (disp>0) then
           read(io,pos=disp, iostat=i_err) cdummy
           if (i_err/=0) then
              call cgyro_error('[REWIND] Failed to rewind '//datafile)
              close(io)
              return
           endif
        endif
        endfile(io)
        close(io)

     endif

  end select

end subroutine cgyro_write_distributed_bcomplex

!------------------------------------------------------
! cgyro_write_distributed_breal.f90
!
! PURPOSE:
!  Control merged output of real distributed array.
!------------------------------------------------------

subroutine cgyro_write_distributed_breal(datafile,n_fn,fn)

  use mpi
  use cgyro_globals
  use cgyro_io, only: cgyro_error

  !------------------------------------------------------
  implicit none
  !
  character (len=*), intent(in) :: datafile
  integer, intent(in) :: n_fn
  real, intent(in) :: fn(n_fn)
  !
  ! Required for MPI-IO:
  !
  integer :: filemode
  integer :: finfo
  integer :: fh
  integer :: fstatus(MPI_STATUS_SIZE)
  integer(kind=MPI_OFFSET_KIND) :: disp
  integer(kind=MPI_OFFSET_KIND) :: offset1
  !
  character :: cdummy
  !------------------------------------------------------

  if (i_proc_1 /= 0) return

  select case (io_control)

  case (0)

     return

  case (1)

     ! Open

     if (i_proc == 0) then
        open(unit=io,file=datafile,status='replace')
        close(io)
     endif

  case (2)

     ! Append

     ! Write in parallel to the binary datafile
     filemode = MPI_MODE_WRONLY
     disp = i_current-1
     disp = disp*n_proc_2*size(fn)*BYTE

     offset1 = i_proc_2*size(fn)

     call MPI_INFO_CREATE(finfo,i_err)

     if (mpiio_small_stripe_factor > 0) then
        ! user asked us to explicitly set the MPI IO striping factor
        call MPI_INFO_SET(finfo,"striping_factor",mpiio_small_stripe_str,i_err)
     endif

     call MPI_FILE_OPEN(NEW_COMM_2,&
          datafile,&
          filemode,&
          finfo,&
          fh,&
          i_err)

     if (BYTE == 4) then

        ! Single (default) 
        call MPI_FILE_SET_VIEW(fh,&
             disp,&
             MPI_REAL4,&
             MPI_REAL4,&
             'native',&
             finfo,&
             i_err)

        call MPI_FILE_WRITE_AT(fh,&
             offset1,&
             real(fn,kind=4),&
             n_fn,&
             MPI_REAL4,&
             fstatus,&
             i_err)
     else

        call MPI_FILE_SET_VIEW(fh,&
             disp,&
             MPI_REAL8,&
             MPI_REAL8,&
             'native',&
             finfo,&
             i_err)

        call MPI_FILE_WRITE_AT(fh,&
             offset1,&
             fn,&
             n_fn,&
             MPI_REAL8,&
             fstatus,&
             i_err)
     endif

     call MPI_FILE_SYNC(fh,i_err)
     call MPI_FILE_CLOSE(fh,i_err)
     call MPI_INFO_FREE(finfo,i_err)

  case (3)

     ! Rewind

     if (i_proc == 0) then

        disp = i_current
        disp = disp*n_proc_2*size(fn)*BYTE

        open(unit=io,file=datafile,status='old',access='stream',iostat=i_err)
        if (i_err/=0) then
           call cgyro_error('[REWIND] Failed to open '//datafile)
           return
        endif
        if (disp>0) then
           read(io,pos=disp, iostat=i_err) cdummy
           if (i_err/=0) then
              call cgyro_error('[REWIND] Failed to rewind '//datafile)
              close(io)
              return
           endif
        endif
        endfile(io)
        close(io)

     endif

  end select

end subroutine cgyro_write_distributed_breal

!------------------------------------------------------
! write_precision.f90
!
! PURPOSE:
!  Reduce across n and then write precision scalar.
!------------------------------------------------------

subroutine write_precision(datafile,fn)

  use mpi
  use cgyro_globals

  !------------------------------------------------------
  implicit none
  !
  character (len=*), intent(in) :: datafile
  real, intent(in) :: fn
  real :: fn_sum
  integer :: i
  !------------------------------------------------------

  call MPI_ALLREDUCE(fn, &
       fn_sum, &
       1, &
       MPI_DOUBLE_PRECISION, &
       MPI_SUM, &
       NEW_COMM_2, &
       i_err)

  if (i_proc > 0) return

  select case (io_control)

  case(0)

     return

  case(1)

     ! Open

     open(unit=io,file=datafile,status='replace')
     close(io)

  case(2)

     ! Append

     open(unit=io,file=datafile,status='old',position='append')
     write(io,fmtstr_hi) fn_sum
     close(io)

  case(3)

     ! Rewind

     open(unit=io,file=datafile,status='old')
     do i=1,i_current
        read(io,fmtstr_hi) fn_sum
     enddo
     endfile(io)
     close(io)

  end select

end subroutine write_precision

!----------------------------------------------------------------
! write_ascii.f90
!
! PURPOSE:
!  ASCII output of vector of reals 
!----------------------------------------------------------------

subroutine write_ascii(datafile,n_fn,fn)

  use cgyro_globals

  !------------------------------------------------------
  implicit none
  !
  character (len=*), intent(in) :: datafile
  integer, intent(in) :: n_fn
  real, intent(inout), dimension(n_fn) :: fn
  integer :: i
  !------------------------------------------------------

  if (i_proc > 0) return

  select case (io_control)

  case(0)

     return

  case(1)

     ! Open

     open(unit=io,file=datafile,status='replace')
     close(io)

  case(2)

     ! Append

     open(unit=io,file=datafile,position='append')
     write(io,fmtstrn) fn(:)  
     close(io)

  case(3)

     ! Rewind

     open(unit=io,file=datafile,status='old')
     do i=1,i_current
        read(io,fmtstr) fn(1)
     enddo
     endfile(io)
     close(io)

  end select

end subroutine write_ascii

! NOTE: Can only be called with a single my_toroidal
subroutine write_distribution(datafile)

  use mpi
  use cgyro_globals

  !------------------------------------------------------
  implicit none
  !
  character (len=*), intent(in) :: datafile
  integer :: ir,it
  complex, dimension(:,:), allocatable :: h_x_glob
  complex :: ftemp(n_theta,n_radial)
  !------------------------------------------------------

  select case (io_control)

  case(0)

     return

  case(1)

     ! Open
     if (i_proc == 0) then
        open(unit=io,file=datafile,status='replace')
        close(io)
     endif

  case(2)

     ! Append

     if (i_proc == 0) then
        open(unit=io,file=datafile,status='old',position='append',access='stream')
     endif

     allocate(h_x_glob(nc,nv))
     ! Collect distribution onto process 0
     call MPI_GATHER(cap_h_c(:,:,nt1),&
          size(cap_h_c(:,:,nt1)),&
          MPI_DOUBLE_COMPLEX,&
          h_x_glob(:,:),&
          size(h_x_glob),&
          MPI_DOUBLE_COMPLEX,&
          0,&
          NEW_COMM_1,&
          i_err)

     if (i_proc == 0) then
        do iv=1,nv
           do ir=1,n_radial
              do it=1,n_theta
                 ftemp(it,ir) = h_x_glob(ic_c(ir,it),iv)
              enddo
           enddo
           if (zf_test_mode == 0) call extended_ang(ftemp)
           write(io) cmplx(ftemp,kind=4)
        enddo
        close(io)
     endif

     deallocate(h_x_glob)

  case(3)

     ! Rewind not implemented

  end select

end subroutine write_distribution

!----------------------------------------------------------------
! write_timers.f90
!
! PURPOSE:
!
! Initialization:  stream_init, coll_init
! Runtime: field,stream,nl,nl_comm,coll,coll_field,coll_comm,io
!
!----------------------------------------------------------------

subroutine write_timers(datafile)

  use cgyro_globals
  use timer_lib

  !-----------------------------------------------
  implicit none
  !
  character (len=*), intent(in) :: datafile
  integer :: i_dummy
  real, dimension(11) :: dummy
  character (len=2) :: sdummy
  !-------------------------------------------------

  if (io_control == 1 .or. io_control == 3) then
     ! Timer initialization (starts at timer 6)
     call timer_lib_init('str')
     call timer_lib_init('str_mem')
     call timer_lib_init('str_comm')
     call timer_lib_init('nl')
     call timer_lib_init('nl_mem')
     call timer_lib_init('nl_comm')
     call timer_lib_init('field')
     call timer_lib_init('field_com')
     call timer_lib_init('shear')
     call timer_lib_init('coll')
     call timer_lib_init('coll_mem')
     call timer_lib_init('coll_comm')
     call timer_lib_init('io')
     call timer_lib_init('TOTAL')
  endif

  select case (io_control)

  case(0)

     return

  case (1)

     ! Initial open
     if (i_proc == 0) then
        open(unit=io,file=datafile,status='replace')
        write(io,'(a)') 'Setup time'
        ! align for f10.3, i.e. at least 6 chars, last empty
        write(io,'(1x,9(a10,2x))') &
             '    input ', &
             ' str_init ', &
             '  nl_init ', &
             'coll_init ', &
             '  io_init '
        write(io,'(9(f10.3,2x))') &
             timer_lib_time('input'),&
             timer_lib_time('str_init'),&
             timer_lib_time('nl_init'),&
             timer_lib_time('coll_init'),&
             timer_lib_time('io_init')
        write(io,'(a)') 'Run time'
        ! align for f10.3, i.e. at least 6 chars, last empty
        write(io,'(1x,14(a10,1x))') &
             '    str   ', &
             '  str_mem ', &
             ' str_comm ', &
             '    nl    ', &
             '   nl_mem ', &
             '  nl_comm ', &
             '    field ', &
             'field_com ', &
             '    shear ', &
             '    coll  ', &
             ' coll_mem ', &
             'coll_comm ', &
             '    io    ', &
             '    TOTAL '
        close(io)
     endif

  case (2)

     !---------------------------------------------------------------------------
     ! Print timers

     call MPI_barrier(CGYRO_COMM_WORLD,i_err)

     if (i_proc == 0) then
        open(unit=io,file=datafile,status='old',position='append')
        write(io,'(14(f10.3,1x))') &
             timer_lib_time('str'),& 
             timer_lib_time('str_mem'),&
             timer_lib_time('str_comm'),& 
             timer_lib_time('nl'),& 
             timer_lib_time('nl_mem'),&
             timer_lib_time('nl_comm'),&
             timer_lib_time('field'),&
             timer_lib_time('field_com'),&
             timer_lib_time('shear'),&
             timer_lib_time('coll'),&
             timer_lib_time('coll_mem'),&
             timer_lib_time('coll_comm'),&
             timer_lib_time('io'),& 
             timer_lib_time('TOTAL') 
        close(io)
     endif
     !---------------------------------------------------------------------------

     ! Reset all timers
     timer_cpu = 0.0

  case (3)

     ! Rewind

     if (i_proc == 0) then 
        open(unit=io,file=datafile,status='old')
        read(io,'(a)') sdummy
        read(io,'(a)') sdummy
        read(io,'(a)') sdummy
        read(io,'(a)') sdummy
        read(io,'(a)') sdummy
        do i_dummy=1,i_current
           read(io,*) dummy(1)
        enddo
        endfile(io,iostat=i_dummy)
        close(io)
     endif

     ! Reset all timers
     timer_cpu = 0.0

  end select

end subroutine write_timers


subroutine print_scrdata()

  use cgyro_globals

  !------------------------------------------------------
  implicit none
  !------------------------------------------------------

  ! Lots of immediate return conditions
  if (restart_flag == 1 .and. i_time == 0) return
  if (io_control == 0 .or. i_proc > 0) return

  ! Different output for 1-mode (linear) or multiple-mode runs

  if (n_toroidal > 1) then
     print '(a,1pe9.3,a,1pe9.3,1x,1pe9.3,a)',&
          '[t: ',t_current,&
          '][e: ',integration_error(:),']'
  else
     print '(a,1pe9.3,a,1pe10.3,1x,1pe10.3,a,1pe10.3,a,1pe9.3,1x,1pe9.3,a)',&
          '[t: ',t_current,&
          '][w: ',freq(nt1),&
          '][dw:',abs(freq_err),&
          '][e: ',integration_error(:),']'

  endif

end subroutine print_scrdata

subroutine write_binary(datafile,fn,n_fn)

  use cgyro_globals
  use cgyro_io, only: cgyro_error

  !------------------------------------------------------
  implicit none
  !
  character (len=*), intent(in) :: datafile
  integer, intent(in) :: n_fn
  complex, intent(in) :: fn(n_fn)
  integer :: disp
  character :: cdummy
  !------------------------------------------------------

  if (i_proc > 0) return

  select case (io_control)

  case(0)

     return

  case (1)

     ! Open

     open(unit=io,file=datafile,status='replace')
     close(io)

  case (2)

     ! Append

     open(unit=io,file=datafile,status='old',position='append',access='stream')

     write(io) cmplx(fn,kind=4)
     close(io)

  case (3)

     ! Rewind

     disp = i_current
     disp = disp*size(fn)*BYTE*2

     open(unit=io,file=datafile,status='old',access='stream', iostat=i_err)
     if (i_err/=0) then
        call cgyro_error('[REWIND] Failed to open '//datafile)
        return
     endif
     if (disp>0) then
        read(io,pos=disp, iostat=i_err) cdummy
        if (i_err/=0) then
           call cgyro_error('[REWIND] Failed to rewind '//datafile)
           close(io)
           return
        endif
     endif
     endfile(io)
     close(io)

  end select

end subroutine write_binary

!----------------------------------------------------------------
! extended_ang.f90
!
! PURPOSE:
!
! Map from (r,theta) to extended angle (f2d -> f1d)
!----------------------------------------------------------------
! NOTE: Can only be called with a single my_toroidal

subroutine extended_ang(f2d)

  use cgyro_globals    

  implicit none

  integer :: ir
  complex, intent(inout), dimension(n_theta,n_radial) :: f2d
  complex, dimension(n_theta,n_radial) :: f1d 

  ! Assumption is that box_size=1

  do ir=1,n_radial
     f1d(:,ir) = f2d(:,ir)*exp(-2*pi*i_c*(px(ir)+px0)*k_theta_base*nt1*rmin*sign_qs)
  enddo

  if (sign_qs < 0) then
     ! Reverse output direction
     do ir=1,n_radial
        f2d(:,ir) = f1d(:,n_radial-ir+1)
     enddo
  else
     f2d = f1d
  endif

end subroutine extended_ang
