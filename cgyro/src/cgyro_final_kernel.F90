!-----------------------------------------------------------------
! cgyro_final_kernel.f90
!  
! PURPOSE:  
!  Finalizations/deallocation for cgyro_kernel
!
! NOTES:
!  This can be called directly using the driver routine cgyro 
!  (in which case input data will read from input.dat) or called 
!  as a subroutine using cgyro_sub.
!-----------------------------------------------------------------

subroutine cgyro_final_kernel
  use timer_lib
  use mpi
  use cgyro_globals
  use cgyro_io

  implicit none

  character(len=30) :: final_msg
  real :: exit_dt
  
  call timer_lib_in('coll_mem')
#if defined(OMPGPU)
!$omp target update from(field,cap_h_c,h_x,source,rhs(:,:,:,1))
#elif defined(_OPENACC)
!$acc update host(field,cap_h_c,h_x,source,rhs(:,:,:,1))
#endif
  call timer_lib_out('coll_mem')

  ! Manage exit message

  if (error_status == 0) then
     if (nonlinear_flag == 1) then
        final_msg = 'Normal'
     else
        if (signal == 1) then
           final_msg = 'Linear converged'
        else
           final_msg = 'Linear terminated at max time'
        endif
     endif
     if (silent_flag == 0 .and. i_proc == 0) then
        open(unit=io,file=trim(path)//runfile_info,status='old',position='append')
        write(io,'(a)') 'EXIT: (CGYRO) '//trim(final_msg)
        close(io)
     endif
  endif
  
  call system_clock(kernel_exit_time,kernel_count_rate,kernel_count_max)
  if (kernel_exit_time .gt. kernel_start_time) then
    exit_dt = (kernel_exit_time-kernel_start_time)/real(kernel_count_rate)
  else
    exit_dt = (kernel_exit_time-kernel_start_time+kernel_count_max)/real(kernel_count_rate)
 endif
 
 call cgyro_cleanup
 
end subroutine cgyro_final_kernel
