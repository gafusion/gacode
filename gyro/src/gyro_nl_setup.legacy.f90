!-----------------------------------------------------
! gyro_nl_setup.f90
!
! PURPOSE:
!  Allocate (and define) selected arrays for use with 
!  direct nonlinear method only.
!-------------------------------------------------------

subroutine gyro_nl_setup

  use gyro_globals
  use gyro_nl_private

  !------------------------------------
  implicit none
  !------------------------------------

  !----------------------------------------------------
  do nn=-n_max,n_max
     i_p(nn) = abs(nn)+1
     if (nn < 0) then 
        n_p(nn) = -n(i_p(nn))
     else
        n_p(nn) = n(i_p(nn))
     endif
  enddo
  !
  ! Nonlinear coupling coefficient:
  !
  c_nl_i(:) = -rhos_norm*q_s(:)/(r_s(:)*b_unit_s(:))
  !
  ! ... and add the nonuniform grid effect:
  !
  c_nl_i(:) = dr_eodr(:)*c_nl_i(:)
  !----------------------------------------------------

end subroutine gyro_nl_setup
