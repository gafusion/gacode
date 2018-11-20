!------------------------------------------------------
! gyro_write_timedata.F90
!
! PURPOSE:
!  This is the master file controlling output of
!  data on a per-timestep basis.  This file also 
!  contains the MPI IO routines 
!
!  - write_distributed_real
!  - write_distributed_complex
!  - write_local_real
!-----------------------------------------------------

subroutine gyro_write_timedata

  use gyro_globals
  use mpi

  !---------------------------------------------------
  implicit none
  !
  real, dimension(:,:,:), allocatable :: a3
  !
  complex, dimension(n_theta_plot,n_x,n_kinetic) :: n_plot
  complex, dimension(n_theta_plot,n_x,n_kinetic) :: e_plot
  complex, dimension(n_theta_plot,n_x,n_kinetic) :: v_plot
  !---------------------------------------------------

  !---------------------------------------------------
  ! Timestep data:
  !
  if (i_proc == 0) then
     call gyro_write_step(trim(path)//'out.gyro.t',1)
  endif
  !---------------------------------------------------

  !--------------------------------------------------
  ! Output of field-like quantities:
  !
  if (plot_n_flag+plot_e_flag+plot_v_flag > 0) then
     n_plot(:,:,:) = moments_plot(:,:,:,1)
     e_plot(:,:,:) = moments_plot(:,:,:,2)
     v_plot(:,:,:) = moments_plot(:,:,:,3)
  endif
  !
  if (plot_u_flag == 1) then

     ! POTENTIALS

     call write_distributed_complex(&
          trim(path)//'out.gyro.moment_u',&
          10,&
          size(phi_plot(:,:,1:n_field)),&
          phi_plot(:,:,1:n_field))

  endif !u_flag==1

  if (plot_epar_flag == 1) then

     ! PARALLEL ELECTRIC FIELD

     call write_distributed_complex(&
          trim(path)//'out.gyro.moment_epar',&
          10,&
          size(phi_plot(:,:,n_field+1)),&
          phi_plot(:,:,n_field+1))

  endif !epar_flag==1

  if (plot_n_flag == 1) then

     ! DENSITY

     call write_distributed_complex(&
          trim(path)//'out.gyro.moment_n',&
          10,&
          size(n_plot),&
          n_plot)

  endif !n_flag ==1 

  if (plot_e_flag == 1) then

     ! ENERGY

     call write_distributed_complex(&
          trim(path)//'out.gyro.moment_e',&
          10,&
          size(e_plot),&
          e_plot)

  endif !e_flag==1

  if (plot_v_flag == 1) then

     ! PARALLEL VELOCITY

     call write_distributed_complex(&
          trim(path)//'out.gyro.moment_v',&
          10,&
          size(v_plot),&
          v_plot)

  endif !v_flag==1

  !--------------------------------------------------

  call gyro_kxky_spectrum

  call write_distributed_real(&
       trim(path)//'out.gyro.kxkyspec',&
       10,&
       size(kxkyspec),&
       kxkyspec)

  if (i_proc == 0) then
     call write_local_real(&
          trim(path)//'out.gyro.k_perp_squared',&
          10,&
          size(k_perp_squared),&
          k_perp_squared)
  endif

  call gyro_field_fluxave

  !-------------------------------------------------------------------
  ! Calculation of fundamental nonlinear fluxes
  !
  call gyro_nonlinear_flux
  call gyro_gbflux
  !-------------------------------------------------------------------

  !-------------------------------------------------------------------
  ! Output specific to linear/nonlinear operation:
  !
  if (nonlinear_flag == 0) then

     !=============
     ! BEGIN LINEAR 
     !=============

     call gyro_write_freq(trim(path)//'out.gyro.freq',10)

     if (plot_u_flag == 1) then        

        ! PHI
        call gyro_ballooning_mode(trim(path)//'out.gyro.balloon_phi',10,1,0)

        if (n_field > 1) then
           ! A_PARALLEL 

           call gyro_ballooning_mode(trim(path)//'out.gyro.balloon_a',10,2,0)

        endif

        if (n_field > 2) then
           ! B_PARALLEL 

           call gyro_ballooning_mode(trim(path)//'out.gyro.balloon_aperp',10,3,0)

        endif

        ! E_PARALLEL
        if (eparallel_plot_flag == 1) then

           call gyro_ballooning_mode(trim(path)//'out.gyro.balloon_epar',10,n_field+1,0)

        endif

     endif

     if (plot_n_flag == 1) then

        ! DENSITY
        if (electron_method /= 3) then

           call gyro_ballooning_mode(trim(path)//'out.gyro.balloon_n_ion',10,5,1)

        endif
        if (electron_method > 1) then 

           call gyro_ballooning_mode(trim(path)//'out.gyro.balloon_n_elec',10,5,indx_e)

        endif
     endif

     if (plot_e_flag == 1) then

        ! ENERGY
        if (electron_method /= 3) then

           call gyro_ballooning_mode(trim(path)//'out.gyro.balloon_e_ion',10,6,1)

        endif
        if (electron_method > 1) then 

           call gyro_ballooning_mode(trim(path)//'out.gyro.balloon_e_elec',10,6,indx_e)

        endif
     endif

     if (plot_v_flag == 1) then

        ! ENERGY
        if (electron_method /= 3) then

           call gyro_ballooning_mode(trim(path)//'out.gyro.balloon_v_ion',10,7,1)

        endif
        if (electron_method > 1) then 

           call gyro_ballooning_mode(trim(path)//'out.gyro.balloon_v_elec',10,7,indx_e)

        endif
     endif

     !-----------------------------------------------------------------
     ! Distribution function data:
     !
     if (n_proc == 1 .and. n_n == 1 .and. dist_print == 1) then
        call gyro_write_h(trim(path)//'out.gyro.hp',trim(path)//'out.gyro.ht',10,11)
     endif
     !-----------------------------------------------------------------

     if (i_proc == 0 .and. lindiff_method > 1) then

        call write_local_real( &
             trim(path)//'out.gyro.gbflux',10,size(gbflux),gbflux)
        call write_local_real( &
             trim(path)//'out.gyro.gbflux_mom',10,size(gbflux_mom),gbflux_mom)
        call write_local_real( &
             trim(path)//'out.gyro.gbflux_i',10,size(gbflux_i),gbflux_i)

        if (trapdiff_flag == 1) then
           call write_local_real( &
                trim(path)//'out.gyro.gbflux_trapped',&
                10,size(gbflux_trapped),gbflux_trapped)
           call write_local_real( &
                trim(path)//'out.gyro.gbflux_i_trapped',&
                10,size(gbflux_i_trapped),gbflux_i_trapped)
        endif

     endif !i_proc ==0 and lindiff >1 

     !=============
     ! END LINEAR 
     !=============

  else

     !================
     ! BEGIN NONLINEAR 
     !================

     call write_distributed_real(&
          trim(path)//'out.gyro.gbflux_n',&
          10,&
          size(gbflux_n),&
          gbflux_n)

     if (i_proc == 0) then

        call write_local_real(trim(path)//'out.gyro.field_rms',10,size(ave_phi),ave_phi)

        call write_local_real( &
             trim(path)//'out.gyro.gbflux',10,size(gbflux),gbflux)
        call write_local_real( &
             trim(path)//'out.gyro.gbflux_i',10,size(gbflux_i),gbflux_i)

        call write_local_real( &
             trim(path)//'out.gyro.gbflux_mom',10,size(gbflux_mom),gbflux_mom)
        call write_local_real( &
             trim(path)//'out.gyro.gbflux_exc',10,size(gbflux_exc),gbflux_exc)

        if (trapdiff_flag == 1) then
           call write_local_real( &
                trim(path)//'out.gyro.gbflux_trapped',10,&
                size(gbflux_trapped),gbflux_trapped)
           call write_local_real( &
                trim(path)//'out.gyro.gbflux_i_trapped',10,&
                size(gbflux_i_trapped),gbflux_i_trapped)
        endif !trapdiff_flag == 1

        call write_local_real( &
             trim(path)//'out.gyro.zerobar',10,&
             size(field_fluxave),transpose(field_fluxave))


        allocate(a3(n_kinetic,4,n_x))
        do i=1,n_x
           a3(:,1,i) = h0_n(:,i)
           a3(:,2,i) = h0_e(:,i)
           a3(:,3,i) = source_n(:,i)
           a3(:,4,i) = source_e(:,i)
        enddo

        call write_local_real( &
             trim(path)//'out.gyro.source',10,size(a3),a3)

        call write_local_real( &
             trim(path)//'out.gyro.moment_zero',10,&
             size(moments_zero_plot),moments_zero_plot)

        deallocate(a3)

     endif! i_proc ==0

     !================
     ! END NONLINEAR 
     !================

  endif
  !-------------------------------------------------------------------

  call gyro_write_error(trim(path)//'out.gyro.error',10)

  !------------------------------------------------------------
  ! Entropy diagnostics
  !
  if (entropy_flag == 1) then
     call gyro_entropy 
     if (i_proc == 0) then 
        call write_local_real(&
             trim(path)//'out.gyro.entropy.out',10,size(entropy),entropy)
     endif
  endif
  !------------------------------------------------------------

  !------------------------------------------------------------
  ! Write precision-monitoring data
  !
  call gyro_write_precision(10,sum(abs(gbflux)))
  !------------------------------------------------------------

  !------------------------------------------------------------
  call gyro_write_timers(trim(path)//'out.gyro.timing',10)
  !------------------------------------------------------------

  if (i_proc == 0 .and. debug_flag == 1) print *,'[gyro_write_timedata done]'

end subroutine gyro_write_timedata

!===========================================================================
!------------------------------------------------------
! write_distributed_real.f90
!
! PURPOSE:
!  Control merged output of real distributed array.
!------------------------------------------------------

subroutine write_distributed_real(datafile,io,n_fn,fn)

  use mpi
  use gyro_globals

  !------------------------------------------------------
  implicit none
  !
  character (len=*), intent(in) :: datafile
  integer, intent(in) :: n_fn
  real, intent(in) :: fn(n_fn)
  !
  integer, intent(in) :: io
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
  character(len=fmtstr_len*n_fn) :: fnstr
  character(len=fmtstr_len) :: tmpstr
  character :: c
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

     ! Create human readable string ready to be written
     ! Do it in one shot, to minimize IO
     do in=1,n_fn
        write(tmpstr, fmtstr) fn(in)
        fnstr((in-1)*fmtstr_len+1:in*fmtstr_len-1) = tmpstr(1:11)
        fnstr(in*fmtstr_len:in*fmtstr_len) = NEW_LINE('A')
     enddo

     ! now write in parallel to the common file
     filemode = MPI_MODE_WRONLY
     disp = data_step
     disp = disp*n_proc_2
     disp = disp*fmtstr_len*n_fn

     offset1 = i_proc_2
     offset1 = offset1*fmtstr_len*n_fn

     call MPI_INFO_CREATE(finfo,i_err)

     call MPI_INFO_SET(finfo,"striping_factor","4",i_err)

     call MPI_FILE_OPEN(NEW_COMM_2,&
          datafile,&
          filemode,&
          finfo,&
          fh,&
          i_err)

     call MPI_FILE_SET_VIEW(fh,&
          disp,&
          MPI_CHAR,&
          MPI_CHAR,&
          'native',&
          finfo,&
          i_err)

     call MPI_FILE_WRITE_AT(fh,&
          offset1,&
          fnstr,&
          n_fn*fmtstr_len,&
          MPI_CHAR,&
          fstatus,&
          i_err)

     call MPI_FILE_SYNC(fh,i_err)
     call MPI_FILE_CLOSE(fh,i_err)
     call MPI_INFO_FREE(finfo,i_err)

  case(3)

     ! Rewind

     if (i_proc == 0) then

        disp = data_step+1
        disp = disp*n_proc_2
        disp = disp*fmtstr_len*n_fn

        open(unit=io,file=datafile,status='old',access='STREAM')
        if (disp > 0) then
           read(io,pos=disp) c
        endif
        endfile(io)
        close(io)

     endif

  end select

end subroutine write_distributed_real

!------------------------------------------------------
! write_distributed_complex.f90
!
! PURPOSE:
!  Control merged output of complex distributed array.
!------------------------------------------------------

subroutine write_distributed_complex(datafile,io,n_fn,fn)

  use mpi
  use gyro_globals

  !------------------------------------------------------
  implicit none
  !
  character (len=*), intent(in) :: datafile
  integer, intent(in) :: n_fn
  complex, intent(in) :: fn(n_fn)
  !
  integer, intent(in) :: io
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
  character(len=fmtstr_len*n_fn*2) :: fnstr
  character(len=fmtstr_len) :: tmpstr
  character :: c
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

  case(2)

     ! Append

     ! Create human readable string ready to be written
     ! Do it in one shot, to minimize IO
     do in=1,n_fn
        write(tmpstr, fmtstr) real(fn(in))
        fnstr((in-1)*fmtstr_len*2+1:(in-1)*fmtstr_len*2+fmtstr_len-1) = tmpstr(1:11)
        fnstr((in-1)*fmtstr_len*2+fmtstr_len:(in-1)*fmtstr_len*2+fmtstr_len) = NEW_LINE('A')
        write(tmpstr, fmtstr) aimag(fn(in))
        fnstr((in-1)*fmtstr_len*2+fmtstr_len+1:in*fmtstr_len*2-1) = tmpstr(1:11)
        fnstr(in*fmtstr_len*2:in*fmtstr_len*2) = NEW_LINE('A')
     enddo

     ! now write in parallel to the common file
     filemode = MPI_MODE_WRONLY
     disp = data_step
     disp = disp*n_proc_2
     disp = disp*fmtstr_len*2*n_fn

     offset1 = i_proc_2
     offset1 = offset1*fmtstr_len*2*n_fn

     call MPI_INFO_CREATE(finfo,i_err)

     call MPI_INFO_SET(finfo,"striping_factor","4",i_err)

     call MPI_FILE_OPEN(NEW_COMM_2,&
          datafile,&
          filemode,&
          finfo,&
          fh,&
          i_err)

     call MPI_FILE_SET_VIEW(fh,&
          disp,&
          MPI_CHAR,&
          MPI_CHAR,&
          'native',&
          finfo,&
          i_err)

     call MPI_FILE_WRITE_AT(fh,&
          offset1,&
          fnstr,&
          n_fn*fmtstr_len*2,&
          MPI_CHAR,&
          fstatus,&
          i_err)

     call MPI_FILE_SYNC(fh,i_err)
     call MPI_FILE_CLOSE(fh,i_err)
     call MPI_INFO_FREE(finfo,i_err)

  case(3)

     ! Rewind

     if (i_proc == 0) then

        disp     = data_step+1
        disp     = disp * n_proc_2
        disp = disp * fmtstr_len * 2 * n_fn

        open(unit=io,file=datafile,status='old',access="STREAM")
        if (disp>0) then
           read(io,pos=disp) c
        endif
        endfile(io)
        close(io)

     endif

  end select

end subroutine write_distributed_complex

!------------------------------------------------------
! write_local_real.f90
!
! PURPOSE:
!  This routine write a vector of nondistributed reals.
!------------------------------------------------------

subroutine write_local_real(datafile,io,n_fn,fn)

  use gyro_globals, only : &
       data_step, &
       io_control, &
       fmtstr

  !---------------------------------------------------
  implicit none
  !
  character (len=*), intent(in) :: datafile
  integer, intent(in) :: io
  integer, intent(in) :: n_fn
  real, intent(in) :: fn(n_fn)
  !
  integer :: data_loop
  real :: dummy(n_fn)
  !---------------------------------------------------

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
     write(io,fmtstr)  fn(:)
     close(io)

  case(3)

     ! Rewind

     open(unit=io,file=datafile,status='old')

     do data_loop=0,data_step
        read(io,fmtstr) dummy(:)
     enddo

     endfile(io)
     close(io)

  end select

end subroutine write_local_real

