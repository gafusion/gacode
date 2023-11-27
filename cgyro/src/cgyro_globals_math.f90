!-----------------------------------------------------------------
! cgyro_globals_math.f90
!
! PURPOSE:
!  Helper functions for scalable math operations
!  on matrices with known sizes
!-----------------------------------------------------------------

module cgyro_globals_math

    implicit none 

contains

  !=========================================================
  ! Velocity distributed arrays
  !=========================================================

  !-------------------------------------------------------
  ! Copy one or more arrays
  !-------------------------------------------------------

subroutine cgyro_vel_copy(left, r1)
    use cgyro_globals
    use cgyro_math
    !-------------------------------------------------------
    implicit none
    !
    complex, intent(out), dimension(nc,nv_loc,nt_loc) :: left
    complex, intent(in), dimension(nc,nv_loc,nt_loc) :: r1
    !-------------------------------------------------------
    call cgyro_cmpl_copy(nc*nv_loc*nt_loc, left, r1)
end subroutine cgyro_vel_copy

subroutine cgyro_vel_copy2(left1, left2, r1)
    use cgyro_globals
    use cgyro_math
    !-------------------------------------------------------
    implicit none
    !
    complex, intent(out), dimension(nc,nv_loc,nt_loc) :: left1
    complex, intent(out), dimension(nc,nv_loc,nt_loc) :: left2
    complex, intent(in), dimension(nc,nv_loc,nt_loc) :: r1
    !-------------------------------------------------------
    call cgyro_cmpl_copy2(nc*nv_loc*nt_loc, left1, left2, r1)
end subroutine cgyro_vel_copy2

  !-------------------------------------------------------
  ! Multiple-add of array without using the value of left
  !-------------------------------------------------------

subroutine cgyro_vel_fma2(left, r1, c2, r2, abssum)
    use cgyro_globals
    use cgyro_math
    !-------------------------------------------------------
    implicit none
    !
    complex, intent(out), dimension(nc,nv_loc,nt_loc) :: left
    complex, intent(in), dimension(nc,nv_loc,nt_loc) :: r1
    real, intent(in) :: c2
    complex, intent(in), dimension(nc,nv_loc,nt_loc) :: r2
    real, intent(inout), optional :: abssum
    !-------------------------------------------------------
    call cgyro_cmpl_fma2(nc*nv_loc*nt_loc, left,r1,c2,r2,abssum)
end subroutine cgyro_vel_fma2

subroutine cgyro_vel_fma3(left, r1, c2, r2, c3, r3, abssum)
    use cgyro_globals
    use cgyro_math
    !-------------------------------------------------------
    implicit none
    !
    complex, intent(out), dimension(nc,nv_loc,nt_loc) :: left
    complex, intent(in), dimension(nc,nv_loc,nt_loc) :: r1
    real, intent(in) :: c2
    complex, intent(in), dimension(nc,nv_loc,nt_loc) :: r2
    real, intent(in) :: c3
    complex, intent(in), dimension(nc,nv_loc,nt_loc) :: r3
    real, intent(inout), optional :: abssum
    !-------------------------------------------------------
    call cgyro_cmpl_fma3(nc*nv_loc*nt_loc, left,r1,c2,r2,c3,r3,abssum)
end subroutine cgyro_vel_fma3

subroutine cgyro_vel_fma4(left, r1, c2, r2, c3, r3, c4, r4, abssum)
    use cgyro_globals
    use cgyro_math
    !-------------------------------------------------------
    implicit none
    !
    complex, intent(out), dimension(nc,nv_loc,nt_loc) :: left
    complex, intent(in), dimension(nc,nv_loc,nt_loc) :: r1
    real, intent(in) :: c2
    complex, intent(in), dimension(nc,nv_loc,nt_loc) :: r2
    real, intent(in) :: c3
    complex, intent(in), dimension(nc,nv_loc,nt_loc) :: r3
    real, intent(in) :: c4
    complex, intent(in), dimension(nc,nv_loc,nt_loc) :: r4
    real, intent(inout), optional :: abssum
    !-------------------------------------------------------
    call cgyro_cmpl_fma4(nc*nv_loc*nt_loc, left,r1,c2,r2,c3,r3,c4,r4,abssum)
end subroutine cgyro_vel_fma4

subroutine cgyro_vel_fma5(left, r1, c2, r2, c3, r3, c4, r4, c5, r5, abssum)
    use cgyro_globals
    use cgyro_math
    !-------------------------------------------------------
    implicit none
    !
    complex, intent(out), dimension(nc,nv_loc,nt_loc) :: left
    complex, intent(in), dimension(nc,nv_loc,nt_loc) :: r1
    real, intent(in) :: c2
    complex, intent(in), dimension(nc,nv_loc,nt_loc) :: r2
    real, intent(in) :: c3
    complex, intent(in), dimension(nc,nv_loc,nt_loc) :: r3
    real, intent(in) :: c4
    complex, intent(in), dimension(nc,nv_loc,nt_loc) :: r4
    real, intent(in) :: c5
    complex, intent(in), dimension(nc,nv_loc,nt_loc) :: r5
    real, intent(inout), optional :: abssum
    !-------------------------------------------------------
    call cgyro_cmpl_fma5(nc*nv_loc*nt_loc, left,r1,c2,r2,c3,r3,c4,r4,c5,r5,abssum)
end subroutine cgyro_vel_fma5

subroutine cgyro_vel_fma6(left, r1, c2, r2, c3, r3, c4, r4, c5, r5, c6, r6, abssum)
    use cgyro_globals
    use cgyro_math
    !-------------------------------------------------------
    implicit none
    !
    complex, intent(out), dimension(nc,nv_loc,nt_loc) :: left
    complex, intent(in), dimension(nc,nv_loc,nt_loc) :: r1
    real, intent(in) :: c2
    complex, intent(in), dimension(nc,nv_loc,nt_loc) :: r2
    real, intent(in) :: c3
    complex, intent(in), dimension(nc,nv_loc,nt_loc) :: r3
    real, intent(in) :: c4
    complex, intent(in), dimension(nc,nv_loc,nt_loc) :: r4
    real, intent(in) :: c5
    complex, intent(in), dimension(nc,nv_loc,nt_loc) :: r5
    real, intent(in) :: c6
    complex, intent(in), dimension(nc,nv_loc,nt_loc) :: r6
    real, intent(inout), optional :: abssum
    !-------------------------------------------------------
    call cgyro_cmpl_fma6(nc*nv_loc*nt_loc, &
            left,r1,c2,r2,c3,r3,c4,r4,c5,r5,c6,r6,abssum)
end subroutine cgyro_vel_fma6

subroutine cgyro_vel_fmaN(nr, left, r1, cN, rN, abssum)
    use cgyro_globals
    use cgyro_math
    !-------------------------------------------------------
    implicit none
    !
    integer, intent(in) :: nr
    complex, intent(out), dimension(nc,nv_loc,nt_loc) :: left
    complex, intent(in), dimension(nc,nv_loc,nt_loc) :: r1
    real, intent(in), dimension(nr) :: cN
    complex, intent(in), dimension(nc,nv_loc,nt_loc,nr) :: rN
    real, intent(inout), optional :: abssum
    !-------------------------------------------------------
    integer :: real_nr

    ! allow for tail-end zeros
    real_nr = nr
    do while ((real_nr>1) .AND. (cN(real_nr) == 0.d0 ))
      real_nr = real_nr -1
    enddo

    select case (real_nr)
    case(1)
       call cgyro_vel_fma2(left,r1,cN(1),rN(:,:,:,1),abssum)
    case(2)
       call cgyro_vel_fma3(left,r1,cN(1),rN(:,:,:,1),cN(2),rN(:,:,:,2),abssum)
    case(3)
       call cgyro_vel_fma4(left,r1,cN(1),rN(:,:,:,1),cN(2),rN(:,:,:,2), &
                           cN(3),rN(:,:,:,3),abssum)
    case(4)
       call cgyro_vel_fma5(left,r1,cN(1),rN(:,:,:,1),cN(2),rN(:,:,:,2), &
                           cN(3),rN(:,:,:,3),cN(4),rN(:,:,:,4),abssum)
    case(5)
       call cgyro_vel_fma6(left,r1,cN(1),rN(:,:,:,1),cN(2),rN(:,:,:,2), &
                           cN(3),rN(:,:,:,3),cN(4),rN(:,:,:,4), &
                           cN(5),rN(:,:,:,5),abssum)
    case default
       ! use the generic implementation
       call cgyro_cmpl_fmaN(nc*nv_loc*nt_loc, real_nr, left,r1,cN,rN,abssum)
   end select
end subroutine cgyro_vel_fmaN

  !=========================================================
  ! Specialized merge of 2 FMA with abssum used in gk
  !=========================================================

subroutine cgyro_vel_solution_werror(nr, left, r0, c1, &
                m1, cN, rN, ec1, ecN, abssum_left, abssum_m)
    use cgyro_globals
    use cgyro_math
    !-------------------------------------------------------
    implicit none
    !
    integer, intent(in) :: nr
    complex, intent(inout), dimension(nc,nv_loc,nt_loc) :: left
    complex, intent(in), dimension(nc,nv_loc,nt_loc) :: r0
    real, intent(in) :: c1
    complex, intent(inout), dimension(nc,nv_loc,nt_loc) :: m1
    real, intent(in), dimension(nr) :: cN
    complex, intent(in), dimension(nc,nv_loc,nt_loc,nr) :: rN
    real, intent(in) :: ec1
    real, intent(in), dimension(nr) :: ecN
    real, intent(inout), optional :: abssum_left
    real, intent(inout), optional :: abssum_m
    !-------------------------------------------------------
    call cgyro_cmpl_solution_werror(nc*nv_loc*nt_loc ,nr, &
            left,r0,c1,m1,cN,rN,ec1,ecN,abssum_left,abssum_m)

end subroutine cgyro_vel_solution_werror

end module cgyro_globals_math


