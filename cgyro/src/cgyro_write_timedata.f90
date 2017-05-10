!-----------------------------------------------------------------
! cgyro_write_timedata.f90
!
! PURPOSE:
!  Output of time-dependent data.
!-----------------------------------------------------------------

subroutine cgyro_write_timedata

  use mpi
  use cgyro_globals

  implicit none

  complex :: a_norm
  integer :: i_field,i_moment
  integer :: ir,it
  logical :: lfe
  real :: vfreq(2)
  complex :: ftemp(n_radial,n_theta)
  complex :: field_plot(n_radial,theta_plot)

  ! Print this data on print steps only; otherwise exit now
  if (mod(i_time,print_step) /= 0) return

  ! Increment the print counter on actual output steps
  if (io_control == 2) i_current = i_current+1

  !---------------------------------------------------------------------------
  if (n_toroidal == 1 .and. h_print_flag == 1) then
     call write_distribution(trim(path)//runfile_hb)
  endif
  !---------------------------------------------------------------------------

  call cgyro_flux

  ! ky flux for all species with field breakdown
   call cgyro_write_distributed_real(&
       trim(path)//runfile_ky_flux,&
       size(fflux(:,:,:)),&
       fflux(:,:,:))
 
  if (nonlinear_flag == 1 .and. kxkyflux_print_flag == 1) then
     ! kxky energy flux for all species
     call cgyro_write_distributed_real(&
          trim(path)//runfile_kxky_flux,&
          size(flux(:,:)),&
          flux(:,:))
  endif

  if (nonlinear_flag == 1 .and. globalflux_print_flag == 1) then
     ! Global (n,e,v) fluxes for all species
     do i_moment=1,3
        call cgyro_write_distributed_complex(&
             trim(path)//runfile_lky_flux(i_moment),&
             size(gflux(:,:,i_moment)),&
             gflux(:,:,i_moment))
     enddo
  endif
  
  if (nonlinear_flag == 1 .and. moment_print_flag == 1) then
     ! (n,e) moment for all species at selected thetas.
     do i_moment=1,2
        call cgyro_write_distributed_complex(&
             trim(path)//runfile_kxky(i_moment),&
             size(moment(:,:,:,i_moment)),&
             moment(:,:,:,i_moment))
     enddo
  endif

  ! Sort out subset of theta values for plotting
  do ic=1,nc
     ir = ir_c(ic)
     it = it_c(ic)
     if (itp(it) > 0) then
        field_plot(ir,itp(it)) = field(1,ic)
     endif
  enddo
        
  ! Complex potential at selected thetas
  call cgyro_write_distributed_complex(&
       trim(path)//runfile_kxky_phi,&
       size(field_plot),&
       field_plot)

  ! Checksum for regression testing
  ! Note that checksum is a distributed real scalar
  call write_precision(trim(path)//runfile_prec,sum(abs(fflux)))

  !------------------------------------------------------------------
  ! Ballooning mode (or ZF) output for linear runs with a single mode
  ! (can both be plotted with cgyro_plot -plot ball)
  if (n_toroidal == 1) then
     do i_field=1,n_field

        do ir=1,n_radial
           do it=1,n_theta
              ftemp(ir,it) = field(i_field,ic_c(ir,it))
           enddo
        enddo

        if (i_field == 1) a_norm = ftemp(n_radial/2+1,n_theta/2+1) 

        if (n == 0) then
           call write_zf(&
                trim(path)//runfile_fieldb(i_field),&
                ftemp(:,:)/a_norm)
        else
           call write_balloon(&
                trim(path)//runfile_fieldb(i_field),&
                ftemp(:,:)/a_norm)
        endif
     enddo
  endif
  !---------------------------------------------------------------

  ! Linear frequency diagnostics for every value of n
  call cgyro_freq
  vfreq(1) = real(freq) 
  vfreq(2) = aimag(freq)
  call cgyro_write_distributed_real(trim(path)//runfile_freq,size(vfreq),vfreq)

  ! Output to screen
  call print_scrdata()

  ! Output to files
  call write_time(trim(path)//runfile_time)
  call write_timers(trim(path)//runfile_timers)

  ! Check for manual halt signal
  if (i_proc == 0) then
     inquire(file=trim(path)//'halt',exist=lfe)
     if (lfe .eqv. .true.) then
        open(unit=io,file='halt',status='old')
        read(io,*) signal
        close(io)
     endif
  endif

  call MPI_BCAST(signal,1,MPI_INTEGER,0,CGYRO_COMM_WORLD,i_err)

end subroutine cgyro_write_timedata

!===============================================================================
! Individual parallel-safe I/O routines
!===============================================================================

!------------------------------------------------------
! cgyro_write_distributed_complex.f90
!
! PURPOSE:
!  Control merged output of complex distributed array.
!------------------------------------------------------

subroutine cgyro_write_distributed_complex(datafile,n_fn,fn)

  use mpi
  use cgyro_globals

  !------------------------------------------------------
  implicit none
  !
  character (len=*), intent(in) :: datafile
  integer, intent(in) :: n_fn
  complex, intent(in) :: fn(n_fn)
  !
  integer :: i_group_send
  integer :: i_send
  integer :: in
  integer :: i_dummy
  !
  complex :: fn_recv(n_fn)
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

     if (i_proc == 0) &
          open(unit=io,file=datafile,status='old',position='append')

     do in=1,n_toroidal

        !-----------------------------------------
        ! Subgroup collector:
        !
        i_group_send = in-1

        if (i_group_send /= 0) then

           i_send = i_group_send*n_proc_1

           if (i_proc == 0) then

              call MPI_RECV(fn_recv,&
                   n_fn,&
                   MPI_DOUBLE_COMPLEX,&
                   i_send,&
                   in,&
                   CGYRO_COMM_WORLD,&
                   recv_status,&
                   i_err)

           else if (i_proc == i_send) then

              call MPI_SEND(fn,&
                   n_fn,&
                   MPI_DOUBLE_COMPLEX,&
                   0,&
                   in,&
                   CGYRO_COMM_WORLD,&
                   i_err)

           endif

        else

           fn_recv(:) = fn(:)

        endif
        !
        !-----------------------------------------

        if (i_proc == 0) then

           write(io,fmtstr) fn_recv(:)

        endif

     enddo ! in

     if (i_proc == 0) close(io)

  case(3)

     ! Rewind

     if (i_proc == 0) then

        open(unit=io,file=datafile,status='old')
        do i_dummy=1,i_current

           do in=1,n_toroidal
              read(io,fmtstr) fn_recv(:)
           enddo

        enddo

        endfile(io)
        close(io)

     endif

  end select

end subroutine cgyro_write_distributed_complex


!------------------------------------------------------
! cgyro_write_distributed_real.f90
!
! PURPOSE:
!  Control merged output of real distributed array.
!------------------------------------------------------

subroutine cgyro_write_distributed_real(datafile,n_fn,fn)

  use mpi
  use cgyro_globals

  !------------------------------------------------------
  implicit none
  !
  character (len=*), intent(in) :: datafile
  integer, intent(in) :: n_fn
  real, intent(in) :: fn(n_fn)
  !
  integer :: i_group_send
  integer :: i_send
  integer :: in
  integer :: i_dummy
  !
  real :: fn_recv(n_fn)
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

     if (i_proc == 0) &
          open(unit=io,file=datafile,status='old',position='append')

     do in=1,n_toroidal

        !-----------------------------------------
        ! Subgroup collector:
        !
        i_group_send = in-1

        if (i_group_send /= 0) then

           i_send = i_group_send*n_proc_1

           if (i_proc == 0) then

              call MPI_RECV(fn_recv,&
                   n_fn,&
                   MPI_DOUBLE_PRECISION,&
                   i_send,&
                   in,&
                   CGYRO_COMM_WORLD,&
                   recv_status,&
                   i_err)

           else if (i_proc == i_send) then

              call MPI_SEND(fn,&
                   n_fn,&
                   MPI_DOUBLE_PRECISION,&
                   0,&
                   in,&
                   CGYRO_COMM_WORLD,&
                   i_err)

           endif

        else

           fn_recv(:) = fn(:)

        endif
        !
        !-----------------------------------------

        if (i_proc == 0) then

           write(io,fmtstr) fn_recv(:)

        endif

     enddo ! in

     if (i_proc == 0) close(io)

  case(3)

     ! Rewind

     if (i_proc == 0) then

        open(unit=io,file=datafile,status='old')
        do i_dummy=1,i_current

           do in=1,n_toroidal
              read(io,fmtstr) fn_recv(:)
           enddo

        enddo

        endfile(io)
        close(io)

     endif

  end select

end subroutine cgyro_write_distributed_real


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
  integer :: i_dummy
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
     do i_dummy=1,i_current
        read(io,fmtstr_hi) fn_sum
     enddo
     endfile(io)
     close(io)

  end select

end subroutine write_precision


!------------------------------------------------------
! write_balloon.f90
!
! PURPOSE:
!  Manage IO of ballooning potentials
!------------------------------------------------------

subroutine write_balloon(datafile,fn)

  use cgyro_globals

  !------------------------------------------------------
  implicit none
  !
  character (len=*), intent(in) :: datafile
  complex, intent(in) :: fn(n_radial,n_theta)
  !
  integer :: i_dummy
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

     open(unit=io,file=datafile,status='old',position='append')

     ! Construct ballooning-space form of field
     call extended_ang(fn,f_balloon)

     write(io,fmtstr) transpose(f_balloon(:,:))
     close(io)

     !-------------------------------------------------------

  case(3)

     ! Rewind

     open(unit=io,file=datafile,status='old')
     do i_dummy=1,i_current
        read(io,fmtstr) f_balloon(:,:)
     enddo
     endfile(io)
     close(io)

  end select

end subroutine write_balloon


!----------------------------------------------------------------
! extended_ang.f90
!
! PURPOSE:
!
! Map from (r,theta) to extended angle (f2d -> f1d)
!----------------------------------------------------------------

subroutine extended_ang(f2d,f1d)

  use cgyro_globals    

  implicit none

  integer :: ir,jr,it,np
  complex, intent(in), dimension(n_radial,n_theta) :: f2d
  complex, intent(inout), dimension(n_radial,n_theta) :: f1d 

  np = n_radial/2/box_size

  do ir=-np,np-1
     do it=1,n_theta
        ! Manage positive/negative q
        if (ipccw*btccw > 0) then
           jr = box_size*ir+n_radial/2+1
        else
           jr = -box_size*ir+n_radial/2
        endif
        f1d(ir+np+1,it) = f2d(jr,it)*exp(-2*pi*i_c*ir*k_theta*rmin)
     enddo
  enddo

  if (ipccw*btccw < 0) then
     f1d = f1d*exp(2*pi*i_c*abs(k_theta)*rmin)
  endif

end subroutine extended_ang


!----------------------------------------------------------------
! write_zf.f90
!
! PURPOSE:
!
!----------------------------------------------------------------

subroutine write_zf(datafile,fn)

  use cgyro_globals

  !------------------------------------------------------
  implicit none
  !
  character (len=*), intent(in) :: datafile
  complex, intent(in) :: fn(n_radial,n_theta)
  complex :: ftmp(n_radial,n_theta)
  integer :: i_dummy
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

     open(unit=io,file=datafile,status='old',position='append')
    
     write(io,fmtstr) fn(:,:)
     close(io)

  case(3)

     ! Rewind

     open(unit=io,file=datafile,status='old')
     do i_dummy=1,i_current
        read(io,fmtstr) ftmp(:,:)
     enddo
     endfile(io)
     close(io)

  end select

end subroutine write_zf

!----------------------------------------------------------------
! write_time.f90
!
! PURPOSE:
!
!----------------------------------------------------------------

subroutine write_time(datafile)

  use cgyro_globals

  !------------------------------------------------------
  implicit none
  !
  character (len=*), intent(in) :: datafile
  integer :: i_dummy
  real :: dummy
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

     open(unit=io,file=datafile,status='old',position='append')
     write(io,fmtstrn) t_current,integration_error(:)
     close(io)

  case(3)

     ! Rewind

     open(unit=io,file=datafile,status='old')
     do i_dummy=1,i_current
        read(io,fmtstr2) dummy,dummy
     enddo
     endfile(io)
     close(io)

  end select

end subroutine write_time


subroutine write_distribution(datafile)

  use mpi
  use cgyro_globals

  !------------------------------------------------------
  implicit none
  !
  character (len=*), intent(in) :: datafile
  integer :: ir,it
  complex, dimension(:,:), allocatable :: h_x_glob
  complex :: ftemp(n_radial,n_theta)
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
        open(unit=io,file=datafile,status='old',position='append')
     endif

     allocate(h_x_glob(nc,nv))

     ! Collect distribution onto process 0
     call MPI_GATHER(cap_h_c(:,:),&
          size(h_x),&
          MPI_DOUBLE_COMPLEX,&
          h_x_glob(:,:),&
          size(h_x),&
          MPI_DOUBLE_COMPLEX,&
          0,&
          NEW_COMM_1,&
          i_err)

     if (i_proc == 0) then
        do iv=1,nv
           do ir=1,n_radial
              do it=1,n_theta
                 ftemp(ir,it) = h_x_glob(ic_c(ir,it),iv)
              enddo
           enddo
           call extended_ang(ftemp,f_balloon)
           write(io,fmtstr) transpose(f_balloon(:,:))
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
     ! Timer initialization (starts at timer 4)
     call timer_lib_init('str')
     call timer_lib_init('str_comm')
     call timer_lib_init('nl')
     call timer_lib_init('nl_comm')
     call timer_lib_init('field_h')
     call timer_lib_init('field_H')
     call timer_lib_init('shear')
     call timer_lib_init('coll')
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
        write(io,'(1x,9(a11,1x))') timer_cpu_tag(1:3)
        write(io,'(9(1pe10.3,2x))') &
             timer_lib_time('str_init'),timer_lib_time('coll_init'),timer_lib_time('io_init')
        write(io,'(a)') 'Run time'
        write(io,'(1x,11(a10,1x))') timer_cpu_tag(4:14)
        close(io)
     endif

  case (2)

     !---------------------------------------------------------------------------
     ! Print timers
     if (i_proc == 0) then
        open(unit=io,file=datafile,status='old',position='append')
        write(io,'(11(1pe10.3,1x))') &
             timer_lib_time('str'),& 
             timer_lib_time('str_comm'),& 
             timer_lib_time('nl'),& 
             timer_lib_time('nl_comm'),&
             timer_lib_time('field_h'),&
             timer_lib_time('field_H'),&
             timer_lib_time('shear'),&
             timer_lib_time('coll'),&
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
          '][w: ',freq,&
          '][dw:',abs(freq_err),&
          '][e: ',integration_error(:),']'

  endif

end subroutine print_scrdata
