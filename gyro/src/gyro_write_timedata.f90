!------------------------------------------------------
! gyro_write_timedata.F90
!
! PURPOSE:
!  This is the master file controlling output of
!  data on a per-timestep basis.  This file also 
!  contains the MPI IO routines 
!
!  - write_distributed_real
!  - write_distributed_bcomplex
!  - write_binary
!-----------------------------------------------------

subroutine gyro_write_timedata

  use gyro_globals
  use mpi

  !---------------------------------------------------
  implicit none
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

     call write_distributed_bcomplex(&
          trim(path)//'bin.gyro.moment_u',&
          10,&
          size(phi_plot(:,:,1:n_field)),&
          phi_plot(:,:,1:n_field))

  endif !u_flag==1

  if (plot_epar_flag == 1) then

     ! PARALLEL ELECTRIC FIELD

     call write_distributed_bcomplex(&
          trim(path)//'bin.gyro.moment_epar',&
          10,&
          size(phi_plot(:,:,n_field+1)),&
          phi_plot(:,:,n_field+1))

  endif !epar_flag==1

  if (plot_n_flag == 1) then

     ! DENSITY

     call write_distributed_bcomplex(&
          trim(path)//'bin.gyro.moment_n',&
          10,&
          size(n_plot),&
          n_plot)

  endif !n_flag ==1 

  if (plot_e_flag == 1) then

     ! ENERGY

     call write_distributed_bcomplex(&
          trim(path)//'bin.gyro.moment_e',&
          10,&
          size(e_plot),&
          e_plot)

  endif !e_flag==1

  if (plot_v_flag == 1) then

     ! PARALLEL VELOCITY

     call write_distributed_bcomplex(&
          trim(path)//'bin.gyro.moment_v',&
          10,&
          size(v_plot),&
          v_plot)

  endif !v_flag==1

  !--------------------------------------------------

  call gyro_kxky_spectrum

  call write_distributed_breal(&
       trim(path)//'bin.gyro.kxkyspec',&
       10,&
       size(kxkyspec),&
       kxkyspec)

  if (i_proc == 0 .and. extra_print_flag == 1) then
     call gyro_write_binary(&
          trim(path)//'bin.gyro.k_perp_squared',&
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

        call gyro_write_binary( &
             trim(path)//'bin.gyro.gbflux_i',10,size(gbflux_i),gbflux_i)
        if (extra_print_flag == 1) then
           call gyro_write_binary( &
                trim(path)//'bin.gyro.gbflux_mom',10,size(gbflux_mom),gbflux_mom)
        endif

     endif !i_proc ==0 and lindiff >1 

     !=============
     ! END LINEAR 
     !=============

  else

     !================
     ! BEGIN NONLINEAR 
     !================

     call write_distributed_breal(&
          trim(path)//'bin.gyro.gbflux_n',&
          10,&
          size(gbflux_n),&
          gbflux_n)

     if (i_proc == 0) then

        call gyro_write_binary( &
             trim(path)//'bin.gyro.gbflux_i',10,size(gbflux_i),gbflux_i)

        call gyro_write_binary( &
             trim(path)//'bin.gyro.moment_zero',10,&
             size(moments_zero_plot),moments_zero_plot)

        call gyro_write_binary(&
             trim(path)//'bin.gyro.field_rms',10,size(ave_phi),ave_phi)

        if (extra_print_flag == 1) then
           call gyro_write_binary( &
                trim(path)//'bin.gyro.gbflux_mom',10,size(gbflux_mom),gbflux_mom)
           call gyro_write_binary( &
                trim(path)//'bin.gyro.gbflux_exc',10,size(gbflux_exc),gbflux_exc)
           call gyro_write_binary( &
                trim(path)//'bin.gyro.zerobar',10,&
                size(field_fluxave),transpose(field_fluxave))
        endif

     endif ! i_proc = 0

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
        call gyro_write_binary(&
             trim(path)//'bin.gyro.entropy.out',10,size(entropy),entropy)
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

subroutine gyro_write_binary(datafile,io,n_fn,fn)

  use gyro_globals, only: BYTE,io_control,data_step,i_proc

  !------------------------------------------------------
  implicit none
  !
  character (len=*), intent(in) :: datafile
  integer, intent(in) :: io,n_fn
  real, intent(in) :: fn(n_fn)
  integer :: i_err,disp
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

     write(io) real(fn,kind=4)
     close(io)

  case (3)

     ! Rewind

     disp = data_step+1
     disp = disp*size(fn)*BYTE

     open(unit=io,file=datafile,status='old',access='stream', iostat=i_err)
     if (i_err/=0) then
       call catch_error('ERROR: (CGYRO) [REWIND] Failed to open '//datafile)
       return
     endif
     if (disp>0) then
       read(io,pos=disp, iostat=i_err) cdummy
       if (i_err/=0) then
         call catch_error('ERROR: (CGYRO) [REWIND] Failed to rewind '//datafile)
         close(io)
         return
       endif
     endif
     endfile(io)
     close(io)

  end select

end subroutine gyro_write_binary

subroutine write_distributed_breal(datafile,io,n_fn,fn)

  use mpi
  use gyro_globals

  !------------------------------------------------------
  implicit none
  !
  character (len=*), intent(in) :: datafile
  integer, intent(in) :: io,n_fn
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
     disp = data_step
     disp = disp*n_proc_2*size(fn)*BYTE

     offset1 = i_proc_2*size(fn)

     call MPI_INFO_CREATE(finfo,i_err)

     call MPI_INFO_SET(finfo,"striping_factor",'4',i_err)

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

        disp = data_step+1
        disp = disp*n_proc_2*size(fn)*BYTE

        open(unit=io,file=datafile,status='old',access='stream',iostat=i_err)
        if (i_err/=0) then
          call catch_error('ERROR: (CGYRO) [REWIND] Failed to open '//datafile)
          return
        endif
        if (disp>0) then
          read(io,pos=disp, iostat=i_err) cdummy
          if (i_err/=0) then
            call catch_error('ERROR: (CGYRO) [REWIND] Failed to rewind '//datafile)
            close(io)
            return
          endif
        endif
        endfile(io)
        close(io)

     endif

  end select

end subroutine write_distributed_breal

subroutine write_distributed_bcomplex(datafile,io,n_fn,fn)

  use mpi
  use gyro_globals

  !------------------------------------------------------
  implicit none
  !
  character (len=*), intent(in) :: datafile
  integer, intent(in) :: io,n_fn
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
     disp = data_step
     disp = disp*n_proc_2*size(fn)*BYTE*2

     offset1 = i_proc_2*size(fn)

     call MPI_INFO_CREATE(finfo,i_err)

     call MPI_INFO_SET(finfo,"striping_factor",'4',i_err)

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

        disp = data_step+1
        disp = disp*n_proc_2*size(fn)*BYTE*2

        open(unit=io,file=datafile,status='old',access='stream',iostat=i_err)
        if (i_err/=0) then
          call catch_error('ERROR: (CGYRO) [REWIND] Failed to open '//datafile)
          return
        endif
        if (disp>0) then
          read(io,pos=disp, iostat=i_err) cdummy
          if (i_err/=0) then
            call catch_error('ERROR: (CGYRO) [REWIND] Failed to rewind '//datafile)
            close(io)
            return
          endif
        endif
        endfile(io)
        close(io)

     endif

  end select

end subroutine write_distributed_bcomplex
