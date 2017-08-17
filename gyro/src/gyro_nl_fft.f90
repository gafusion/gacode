!----------------------------------------------------------------
! gyro_nl_fft.fftw3.f90
!
! PURPOSE:
!  This routine evaluates the ExB nonlinearity with periodic or
!  nonperiodic boundary conditions using the (F,G)-conservative 
!  difference scheme with FFT in the toroidal direction.
!
! NOTES:
!  FFTW3 version.
!----------------------------------------------------------------

subroutine gyro_nl_fft

  use gyro_globals
  use gyro_pointers
  use gyro_nl_private
  use math_constants
  use ompdata

  !--------------------------------------------
  implicit none
  !
  complex, dimension(-n_max:n_max,i1_buffer:i2_buffer) :: fn
  complex, dimension(-n_max:n_max,i1_buffer:i2_buffer) :: gn
  complex, dimension(0:n_max,n_x) :: nl
  complex, dimension(0:n_max,i1_buffer:i2_buffer) :: fgp
  !
  complex :: fn_p
  complex :: gn_p
  complex :: fn_r
  complex :: gn_r
  !
  complex :: fgr
  complex :: fgr_p
  complex :: fgp_r
  complex :: fg2
  !
  real :: fg_r
  real :: gf_r
  real :: gf_p
  real :: fg_p
  real :: f_pg_r
  real :: f_rg_p
  !
  complex :: fg_r_c
  complex :: gf_r_c
  complex :: f_pg_r_c
  complex :: f_rg_p_c
  !
  !--------------------------------------------

  do is=1,n_kinetic
     do i_split=1,msplit_SSUB

        ! fn, gn and fgp have ip indices; 
        ! must be zeroed.

!$omp parallel private(fn_r,gn_r,fn_p,gn_p,i,nn,i_diff)
!$omp single
        do i=i1_buffer,0
           fn(:,i)  = (0.0,0.0)
           gn(:,i)  = (0.0,0.0)
           fgp(:,i) = (0.0,0.0)
        enddo
        do i=n_x+1,i2_buffer
           fn(:,i)  = (0.0,0.0)
           gn(:,i)  = (0.0,0.0)
           fgp(:,i) = (0.0,0.0)
        enddo
!$omp end single
        do i = ibeg, iend
           do nn=0,n_max
              gn(nn,i) = h_tran(i,i_split,i_p(nn),is)
              fn(nn,i) = gyro_u_tran(i,i_split,i_p(nn),is)
           enddo ! nn
        enddo ! i
!$omp barrier  ! wait for all fn, gn to be available

        !------------------------------------------------
        !
        !
        do i = ibeg, iend
           v_fft3(:,i,:) = 0.0
           do nn=0,n_max
              fn_r = (0.0,0.0)
              gn_r = (0.0,0.0)
              !---------------------------------------------------------------
              ! df/dr, dg/dr
              !
              do i_diff=-m_dx,m_dx-i_dx
                 fn_r = fn_r+w_d1(i_diff)*fn(nn,i_loop(i+i_diff))
                 gn_r = gn_r+w_d1(i_diff)*gn(nn,i_loop(i+i_diff))
              enddo ! i_diff
              !------------------------------------------------
              ! df/dp, dg/dp
              !
              fn_p = -i_c*n_p(nn)*fn(nn,i)
              gn_p = -i_c*n_p(nn)*gn(nn,i)
              !---------------------------------------------------------------
              ! Dealiasing and wrap-around of 
              ! f, g, df/dp, dg/dp, df/dr, dg/dr
              !--------------------------------------------------
              ! Supervector loading I
              !
              v_fft3(nn,i,1) = real(fn(nn,i))
              v_fft3(nn,i,2) = real(gn(nn,i))
              v_fft3(nn,i,3) = real(fn_p)
              v_fft3(nn,i,4) = real(gn_p)
              v_fft3(nn,i,5) = real(fn_r)
              v_fft3(nn,i,6) = real(gn_r)
              if (nn /= 0) then
                 v_fft3(n_fft-nn,i,1) = aimag(fn(nn,i))
                 v_fft3(n_fft-nn,i,2) = aimag(gn(nn,i))
                 v_fft3(n_fft-nn,i,3) = aimag(fn_p)
                 v_fft3(n_fft-nn,i,4) = aimag(gn_p)
                 v_fft3(n_fft-nn,i,5) = aimag(fn_r)
                 v_fft3(n_fft-nn,i,6) = aimag(gn_r)
              endif

           enddo ! nn
        enddo ! i
!$omp end parallel

        !---------------------------------------------------
        ! Backward FFT
        !
        call dfftw_execute(plan_b)

        !---------------------------------------------------------------
        ! Real space multiplications
        !
!$omp parallel private(fg_r,gf_r,gf_p,fg_p,f_pg_r,f_rg_p)
        do i = ibeg, iend
           do nn=0,n_fft-1
              !f dg/dr
              fg_r   =  vt_fft3(nn,i,1)*vt_fft3(nn,i,6)

              !g df/dr
              gf_r   =  vt_fft3(nn,i,2)*vt_fft3(nn,i,5)

              !g df/dp
              gf_p   =  vt_fft3(nn,i,2)*vt_fft3(nn,i,3)

              !f dg/dp
              fg_p   =  vt_fft3(nn,i,1)*vt_fft3(nn,i,4)

              ! df/dp dg/dr
              f_pg_r =  vt_fft3(nn,i,3)*vt_fft3(nn,i,6)

              ! df/dr dg/dp
              f_rg_p =  vt_fft3(nn,i,5)*vt_fft3(nn,i,4)

              !---------------------------------------------------------------
              ! Supervector loading II
              ! 
              vt_fft3(nn,i,1) = fg_r
              vt_fft3(nn,i,2) = gf_r
              vt_fft3(nn,i,3) = gf_p
              vt_fft3(nn,i,4) = fg_p
              vt_fft3(nn,i,5) = f_pg_r
              vt_fft3(nn,i,6) = f_rg_p   
           enddo !nn
        enddo
!$omp end parallel
        !---------------------------------------------------------------

        !---------------------------------------------------------------
        ! Forward FFT
        !
        call dfftw_execute(plan_f)

        !--------------------------------------------------
        ! Re-Construction of complex arrays
        ! 
        !--------------------------------------------------
        !  g df/dp - f dg/dp, 
!$omp parallel private(fgp_r,fg_r_c,gf_r_c,f_pg_r_c,f_rg_p_c,fgr,fg2,fgr_p)
        do i = ibeg, iend
           fgp(0,i) = (cmplx(v_fft3(0,i,3),0.0) - &
                       cmplx(v_fft3(0,i,4),0.0))/n_fft
           do nn=1,n_max
              fgp(nn,i) = (cmplx(v_fft3(nn,i,3),v_fft3(n_fft-nn,i,3))- &
                           cmplx(v_fft3(nn,i,4),v_fft3(n_fft-nn,i,4)))/n_fft
           enddo ! nn 
        enddo ! i
!$omp barrier   ! ensure all fgp values are available
        !--------------------------------------------------

        !---------------------------------------------------------------
        do i = ibeg, iend
           do nn=0,n_max
              !---------------------------------------------------------------
              ! d/dr (g df/dp - f dg/dp)
              !
              fgp_r = (0.0,0.0)
              do i_diff=-m_dx,m_dx-i_dx
                 fgp_r = fgp_r+w_d1(i_diff)*fgp(nn,i_loop(i+i_diff))
              enddo ! i_diff

              if (nn == 0) then
                 fg_r_c   = cmplx(v_fft3(0,i,1),0.0)
                 gf_r_c   = cmplx(v_fft3(0,i,2),0.0)
                 f_pg_r_c = cmplx(v_fft3(0,i,5),0.0)
                 f_rg_p_c = cmplx(v_fft3(0,i,6),0.0)   
              else
                 fg_r_c   = cmplx(v_fft3(nn,i,1),v_fft3(n_fft-nn,i,1))
                 gf_r_c   = cmplx(v_fft3(nn,i,2),v_fft3(n_fft-nn,i,2))
                 f_pg_r_c = cmplx(v_fft3(nn,i,5),v_fft3(n_fft-nn,i,5))
                 f_rg_p_c = cmplx(v_fft3(nn,i,6),v_fft3(n_fft-nn,i,6)) 
              end if

              !----------------------------------------------------------------
              ! f dg/dr - g df/dr, 
              !  g df/dp - f dg/dp,   (above, before this loop)
              !   df/dp dg/dr - df/dr dg/dp
              !

              !f dg/dr - g df/dr
              fgr = (fg_r_c - gf_r_c)/n_fft

              !df/dp dg/dr - df/dr dg/dp
              fg2 = (f_pg_r_c - f_rg_p_c)/n_fft

              ! d/dp (f dg/dr - g df/dr)
              fgr_p = -i_c*n_p(nn)*fgr

              !------------------------------------------------
              ! Arakawa scheme:
              !
              ! d/dp (f dg/dr - g df/dr)
              !   + d/dr (g df/dp - f dg/dp)
              !       + df/dp dg/dr - df/dr dg/dp
              !
              nl(nn,i) = fgr_p + fgp_r + fg2
              !------------------------------------------------
              ! Finally, update global RHS (use h_tran for efficiency):
              h_tran(i,i_split,i_p(nn),is) = (c_nl_i(i)/3.0)*nl(nn,i)

           enddo
        enddo
!$omp end parallel

     enddo ! i_split 
  enddo ! is

end subroutine gyro_nl_fft

