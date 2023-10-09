!-----------------------------------------------------------------
! cgyro_math.F90
!
! PURPOSE:
!  Helper functions for scalable math operations.
!-----------------------------------------------------------------

module cgyro_math

    implicit none 

contains

  !=========================================================
  ! Copy one or more arrays
  !=========================================================

subroutine cgyro_cmpl_copy(sz, left, r1)
    !-------------------------------------------------------
    implicit none
    !
    integer, intent(in) :: sz
    complex, intent(out), dimension(*) :: left
    complex, intent(in), dimension(*) :: r1
    !
    integer :: i
    !-------------------------------------------------------
#if defined(OMPGPU)
!$omp target teams distribute parallel do simd &
!$omp&         map(from:left(1:sz)) &
!$omp&         map(to:r1(1:sz))
#elif defined(_OPENACC)
!$acc parallel loop independent gang vector &
!$acc&         present(left,r1)
#else
!$omp parallel do 
#endif
    do i=1,sz
       left(i) = r1(i)
    enddo
end subroutine cgyro_cmpl_copy

subroutine cgyro_cmpl_copy2(sz, left1, left2, r1)
    !-------------------------------------------------------
    implicit none
    !
    integer, intent(in) :: sz
    complex, intent(out), dimension(*) :: left1
    complex, intent(out), dimension(*) :: left2
    complex, intent(in), dimension(*) :: r1
    !
    integer :: i
    complex :: tmp
    !-------------------------------------------------------
#if defined(OMPGPU)
!$omp target teams distribute parallel do simd &
!$omp&         private(tmp) &
!$omp&         map(from:left1(1:sz),left2(1:sz)) &
!$omp&         map(to:r1(1:sz))
#elif defined(_OPENACC)
!$acc parallel loop independent gang vector &
!$acc&         present(left1,left2,r1) private(tmp)
#else
!$omp parallel do private(tmp)
#endif
    do i=1,sz
       tmp = r1(i)
       left1(i) = tmp
       left2(i) = tmp
    enddo
end subroutine cgyro_cmpl_copy2

  !=========================================================
  ! Multiple-add of array without using the value of left
  !=========================================================

subroutine cgyro_cmpl_fma2(sz, left, r1, c2, r2, abssum)
    !-------------------------------------------------------
    implicit none
    !
    integer, intent(in) :: sz
    complex, intent(out), dimension(*) :: left
    complex, intent(in), dimension(*) :: r1
    real, intent(in) :: c2
    complex, intent(in), dimension(*) :: r2
    real, intent(inout), optional :: abssum
    !
    integer :: i
    complex :: tmp
    real :: s
    !-------------------------------------------------------
    if (present(abssum)) then
      s = 0.0
#if defined(OMPGPU)
!$omp target teams distribute parallel do simd &
!$omp&         private(tmp) &
!$omp&         map(from:left(1:sz)) &
!$omp&         map(to:r1(1:sz),r2(1:sz)) &
!$omp&         reduction(+:s)
#elif defined(_OPENACC)
!$acc parallel loop independent gang vector &
!$acc&         present(left,r1,r2) private(tmp) reduction(+:s)
#else
!$omp parallel do private(tmp) reduction(+:s)
#endif
      do i=1,sz
        tmp = r1(i) + c2 * r2(i)
        left(i) = tmp
        s = s + abs(tmp)
      enddo
      abssum = s
    else
#if defined(OMPGPU)
!$omp target teams distribute parallel do simd &
!$omp&         map(from:left(1:sz)) &
!$omp&         map(to:r1(1:sz),r2(1:sz))
#elif defined(_OPENACC)
!$acc parallel loop independent gang vector &
!$acc&         present(left,r1,r2)
#else
!$omp parallel do
#endif
      do i=1,sz
        left(i) = r1(i) + c2 * r2(i)
      enddo
    endif
end subroutine cgyro_cmpl_fma2

subroutine cgyro_cmpl_fma3(sz, left, r1, c2, r2, c3, r3, abssum)
    !-------------------------------------------------------
    implicit none
    !
    integer, intent(in) :: sz
    complex, intent(out), dimension(*) :: left
    complex, intent(in), dimension(*) :: r1
    real, intent(in) :: c2
    complex, intent(in), dimension(*) :: r2
    real, intent(in) :: c3
    complex, intent(in), dimension(*) :: r3
    real, intent(inout), optional :: abssum
    !
    integer :: i
    complex :: tmp
    real :: s
    !-------------------------------------------------------
    if (present(abssum)) then
      s = 0.0
#if defined(OMPGPU)
!$omp target teams distribute parallel do simd &
!$omp&         private(tmp) &
!$omp&         map(from:left(1:sz)) &
!$omp&         map(to:r1(1:sz),r2(1:sz),r3(1:sz)) &
!$omp&         reduction(+:s)
#elif defined(_OPENACC)
!$acc parallel loop independent gang vector &
!$acc&         present(left,r1,r2,r3) private(tmp) reduction(+:s)
#else
!$omp parallel do private(tmp) reduction(+:s)
#endif
      do i=1,sz
        tmp = r1(i) + c2 * r2(i) + c3 * r3(i)
        left(i) = tmp
        s = s + abs(tmp)
      enddo
      abssum = s
    else
#if defined(OMPGPU)
!$omp target teams distribute parallel do simd &
!$omp&         map(from:left(1:sz)) &
!$omp&         map(to:r1(1:sz),r2(1:sz),r3(1:sz)) 
#elif defined(_OPENACC)
!$acc parallel loop independent gang vector &
!$acc&         present(left,r1,r2,r3)
#else
!$omp parallel do
#endif
      do i=1,sz
        left(i) = r1(i) + c2 * r2(i) + c3 * r3(i)
      enddo
    endif
end subroutine cgyro_cmpl_fma3

subroutine cgyro_cmpl_fma4(sz, left, r1, c2, r2, c3, r3, c4, r4, abssum)
    !-------------------------------------------------------
    implicit none
    !
    integer, intent(in) :: sz
    complex, intent(out), dimension(*) :: left
    complex, intent(in), dimension(*) :: r1
    real, intent(in) :: c2
    complex, intent(in), dimension(*) :: r2
    real, intent(in) :: c3
    complex, intent(in), dimension(*) :: r3
    real, intent(in) :: c4
    complex, intent(in), dimension(*) :: r4
    real, intent(inout), optional :: abssum
    !
    integer :: i
    complex :: tmp
    real :: s
    !-------------------------------------------------------
    if (present(abssum)) then
      s = 0.0
#if defined(OMPGPU)
!$omp target teams distribute parallel do simd &
!$omp&         private(tmp) &
!$omp&         map(from:left(1:sz)) &
!$omp&         map(to:r1(1:sz),r2(1:sz),r3(1:sz)) &
!$omp&         map(to:r4(1:sz)) &
!$omp&         reduction(+:s)
#elif defined(_OPENACC)
!$acc parallel loop independent gang vector &
!$acc&         present(left,r1,r2,r3,r4) private(tmp) reduction(+:s)
#else
!$omp parallel do private(tmp) reduction(+:s)
#endif
      do i=1,sz
        tmp = r1(i) + c2 * r2(i) + c3 * r3(i)+ c4 * r4(i)
        left(i) = tmp
        s = s + abs(tmp)
      enddo
      abssum = s
    else
#if defined(OMPGPU)
!$omp target teams distribute parallel do simd &
!$omp&         map(from:left(1:sz)) &
!$omp&         map(to:r1(1:sz),r2(1:sz),r3(1:sz)) &
!$omp&         map(to:r4(1:sz)) 
#elif defined(_OPENACC)
!$acc parallel loop independent gang vector &
!$acc&         present(left,r1,r2,r3,r4)
#else
!$omp parallel do
#endif
      do i=1,sz
        left(i) = r1(i) + c2 * r2(i) + c3 * r3(i)+ c4 * r4(i)
      enddo
    endif
end subroutine cgyro_cmpl_fma4

subroutine cgyro_cmpl_fma5(sz, left, r1, c2, r2, c3, r3, c4, r4, c5, r5, abssum)
    !-------------------------------------------------------
    implicit none
    !
    integer, intent(in) :: sz
    complex, intent(out), dimension(*) :: left
    complex, intent(in), dimension(*) :: r1
    real, intent(in) :: c2
    complex, intent(in), dimension(*) :: r2
    real, intent(in) :: c3
    complex, intent(in), dimension(*) :: r3
    real, intent(in) :: c4
    complex, intent(in), dimension(*) :: r4
    real, intent(in) :: c5
    complex, intent(in), dimension(*) :: r5
    real, intent(inout), optional :: abssum
    !
    integer :: i
    complex :: tmp
    real :: s
    !-------------------------------------------------------
    if (present(abssum)) then
      s = 0.0
#if defined(OMPGPU)
!$omp target teams distribute parallel do simd &
!$omp&         private(tmp) &
!$omp&         map(from:left(1:sz)) &
!$omp&         map(to:r1(1:sz),r2(1:sz),r3(1:sz)) &
!$omp&         map(to:r4(1:sz),r5(1:sz)) &
!$omp&         reduction(+:s)
#elif defined(_OPENACC)
!$acc parallel loop independent gang vector &
!$acc&         present(left,r1,r2,r3,r4,r5) private(tmp) reduction(+:s)
#else
!$omp parallel do private(tmp) reduction(+:s)
#endif
      do i=1,sz
        tmp = r1(i) + c2 * r2(i) + c3 * r3(i)+ c4 * r4(i) + c5 * r5(i)
        left(i) = tmp
        s = s + abs(tmp)
      enddo
      abssum = s
    else
#if defined(OMPGPU)
!$omp target teams distribute parallel do simd &
!$omp&         map(from:left(1:sz)) &
!$omp&         map(to:r1(1:sz),r2(1:sz),r3(1:sz)) &
!$omp&         map(to:r4(1:sz),r5(1:sz))
#elif defined(_OPENACC)
!$acc parallel loop independent gang vector &
!$acc&         present(left,r1,r2,r3,r4,r5)
#else
!$omp parallel do
#endif
      do i=1,sz
        left(i) = r1(i) + c2 * r2(i) + c3 * r3(i)+ c4 * r4(i) + c5 * r5(i)
      enddo
    endif
end subroutine cgyro_cmpl_fma5

subroutine cgyro_cmpl_fma6(sz, left, r1, c2, r2, c3, r3, c4, r4, c5, r5, c6, r6, abssum)
    !-------------------------------------------------------
    implicit none
    !
    integer, intent(in) :: sz
    complex, intent(out), dimension(*) :: left
    complex, intent(in), dimension(*) :: r1
    real, intent(in) :: c2
    complex, intent(in), dimension(*) :: r2
    real, intent(in) :: c3
    complex, intent(in), dimension(*) :: r3
    real, intent(in) :: c4
    complex, intent(in), dimension(*) :: r4
    real, intent(in) :: c5
    complex, intent(in), dimension(*) :: r5
    real, intent(in) :: c6
    complex, intent(in), dimension(*) :: r6
    real, intent(inout), optional :: abssum
    !
    integer :: i
    complex :: tmp
    real :: s
    !-------------------------------------------------------
    if (present(abssum)) then
      s = 0.0
#if defined(OMPGPU)
!$omp target teams distribute parallel do simd &
!$omp&         private(tmp) &
!$omp&         map(from:left(1:sz)) &
!$omp&         map(to:r1(1:sz),r2(1:sz),r3(1:sz)) &
!$omp&         map(to:r4(1:sz),r5(1:sz),r6(1:sz)) &
!$omp&         reduction(+:s)
#elif defined(_OPENACC)
!$acc parallel loop independent gang vector &
!$acc&         present(left,r1,r2,r3,r4,r5,r6) private(tmp) reduction(+:s)
#else
!$omp parallel do private(tmp) reduction(+:s)
#endif
      do i=1,sz
        tmp = r1(i) + c2 * r2(i) + c3 * r3(i)+ c4 * r4(i) + c5 * r5(i) + c6 * r6(i)
        left(i) = tmp
        s = s + abs(tmp)
      enddo
      abssum = s
    else
#if defined(OMPGPU)
!$omp target teams distribute parallel do simd &
!$omp&         map(from:left(1:sz)) &
!$omp&         map(to:r1(1:sz),r2(1:sz),r3(1:sz)) &
!$omp&         map(to:r4(1:sz),r5(1:sz),r6(1:sz))
#elif defined(_OPENACC)
!$acc parallel loop independent gang vector &
!$acc&         present(left,r1,r2,r3,r4,r5,r6)
#else
!$omp parallel do
#endif
      do i=1,sz
        left(i) = r1(i) + c2 * r2(i) + c3 * r3(i)+ c4 * r4(i) + c5 * r5(i) + c6 * r6(i)
      enddo
    endif
end subroutine cgyro_cmpl_fma6

! rN should logically be (sz,nr) in size,
subroutine cgyro_cmpl_fmaN(sz, nr, left, r1, cN, rN, abssum)
    !-------------------------------------------------------
    implicit none
    !
    integer, intent(in) :: sz
    integer, intent(in) :: nr
    complex, intent(out), dimension(*) :: left
    complex, intent(in), dimension(*) :: r1
    real, intent(in), dimension(nr) :: cN
    complex, intent(in), dimension(*) :: rN
    real, intent(inout), optional :: abssum
    !
    integer :: i,j
    complex :: tmp
    real :: s
    !-------------------------------------------------------
    if (present(abssum)) then
      s = 0.0
#if defined(OMPGPU)
!$omp target teams distribute parallel do simd &
!$omp&         private(tmp,j) &
!$omp&         map(from:left(1:sz))  map(to:r1(1:sz)) &
!$omp&         map(to:rN(1:sz*nr),cN(1:nr)) &
!$omp&         reduction(+:s)
#elif defined(_OPENACC)
!$acc parallel loop independent gang vector &
!$acc&         present(left,r1,rN) copyin(cN) private(tmp,j) reduction(+:s)
#else
!$omp parallel do private(tmp,j) reduction(+:s)
#endif
      do i=1,sz
        tmp = r1(i)
#if (!defined(OMPGPU)) && defined(_OPENACC)
!$acc loop seq
#endif
        do j=1,nr
           tmp = tmp + cN(j) * rN((j-1)*sz+i)
        enddo
        left(i) = tmp
        s = s + abs(tmp)
      enddo
      abssum = s
    else
#if defined(OMPGPU)
!$omp target teams distribute parallel do simd &
!$omp&         private(tmp,j) &
!$omp&         map(from:left(1:sz))  map(to:r1(1:sz)) &
!$omp&         map(to:rN(1:sz*nr),cN(1:nr))
#elif defined(_OPENACC)
!$acc parallel loop independent gang vector &
!$acc&         present(left,r1,rN) copyin(cN) private(tmp,j)
#else
!$omp parallel do private(tmp,j)
#endif
      do i=1,sz
        tmp = r1(i)
#if (!defined(OMPGPU)) && defined(_OPENACC)
!$acc loop seq
#endif
        do j=1,nr
           tmp = tmp + cN(j) * rN((j-1)*sz+i)
        enddo
        left(i) = tmp
      enddo
    endif
end subroutine cgyro_cmpl_fmaN

  !=========================================================
  ! Specialized merge of 2 FMA with abssum used in gk
  !=========================================================

! rN should logically be (sz,nr) in size,
subroutine cgyro_cmpl_solution_werror(sz, nr, left, r0, c1, m1, cN, rN, ec1, ecN, abssum_left, abssum_m)
    !-------------------------------------------------------
    implicit none
    !
    integer, intent(in) :: sz
    integer, intent(in) :: nr
    complex, intent(inout), dimension(*) :: left
    complex, intent(in), dimension(*) :: r0
    real, intent(in) :: c1
    complex, intent(inout), dimension(*) :: m1
    real, intent(in), dimension(nr) :: cN
    complex, intent(in), dimension(*) :: rN
    real, intent(in) :: ec1
    real, intent(in), dimension(nr) :: ecN
    real, intent(inout) :: abssum_left
    real, intent(inout) :: abssum_m
    !
    integer :: i,j
    complex :: tmp, tmpl, tmpm
    real :: sl,sm
    !-------------------------------------------------------
    sl = 0.0
    sm = 0.0
#if defined(OMPGPU)
!$omp target teams distribute parallel do simd &
!$omp&         private(tmp,tmpl,tmpm,j) &
!$omp&         map(from:left(1:sz)) map(tofrom:m1(1:sz)) &
!$omp&         map(to:rN(1:sz*nr),r0(1:sz),cN(1:nr),ecN(1:nr)) &
!$omp&         reduction(+:sl,sm)
#elif defined(_OPENACC)
!$acc parallel loop independent gang vector &
!$acc&         present(left,r0,m1,rN) copyin(cN,ecN) private(tmp,tmpl,tmpm,j) &
!$acc&         reduction(+:sl,sm)
#else
!$omp parallel do private(tmp,tmpl,tmpm,j) reduction(+:sl,sm)
#endif
    do i=1,sz
       ! compute solution using FMA of r0,m1 and rN using c -> left
       ! also FMA of m1 and rN using ec -> m1
       tmp = m1(i)
       tmpl = r0(i) + c1 * tmp
       tmpm = ec1*tmp
#if (!defined(OMPGPU)) && defined(_OPENACC)
!$acc loop seq private(tmp)
#endif
       do j=1,nr
          tmp = rN((j-1)*sz+i)
          tmpl = tmpl +  cN(j) * tmp
          tmpm = tmpm + ecN(j) * tmp
       enddo
       left(i) = tmpl
       m1(i) = tmpm
       sl = sl + abs(tmpl)
       sm = sm + abs(tmpm)
    enddo
    abssum_left = sl
    abssum_m = sm
end subroutine cgyro_cmpl_solution_werror

end module cgyro_math

