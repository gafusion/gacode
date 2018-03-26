!------------------------------------------------------------------------------
! cgyro_step_collision.F90
!
! PURPOSE:
!  Take an implicit collision step using the pre-computed collision 
!  matrix.  Effectively, we compute the new collisional cap_H: 
!
!                       H = h + ze/T G phi
!------------------------------------------------------------------------------

subroutine cgyro_step_collision

  use parallel_lib
  use timer_lib

  use cgyro_globals

  implicit none

  integer :: is,ivp,it
  complex, dimension(nv) :: bvec,cvec
  real :: cvec_re,cvec_im
  real :: cval

  !----------------------------------------------------------------
  ! Perform data tranpose from _c to _v data layouts:
  call timer_lib_in('coll_comm')
  call parallel_lib_rtrans(cap_h_c,cap_h_v)
  call timer_lib_out('coll_comm')
  !----------------------------------------------------------------

  call timer_lib_in('coll')

!$omp parallel do private(ic,ic_loc, it, iv,ivp,cvec,bvec,cvec_re,cvec_im,cval) firstprivate(collision_model)
  do ic=nc1,nc2

     ic_loc = ic-nc1+1

     it = it_c(ic)

     ! Set-up the RHS: H = f + ze/T G phi

     do iv=1,nv
        cvec(iv) = cap_h_v(ic_loc,iv)
     enddo

     bvec(:) = (0.0,0.0)

     ! This is a key loop for performance
     if (collision_model == 6) then
        do ivp=1,nv
           cvec_re = real(cvec(ivp))
           cvec_im = aimag(cvec(ivp))
           do iv=1,nv
             cval = cmat_base(iv,ivp,it) + cmat_diff(iv,ivp,ic_loc)
             bvec(iv) = bvec(iv)+ cmplx(cval*cvec_re, cval*cvec_im)
           enddo
        enddo
     else
       do ivp=1,nv
           cvec_re = real(cvec(ivp))
           cvec_im = aimag(cvec(ivp))
           do iv=1,nv
             cval = cmat(iv,ivp,ic_loc)
             bvec(iv) = bvec(iv)+ cmplx(cval*cvec_re, cval*cvec_im)
           enddo
        enddo
     endif

     call parallel_lib_f_i_set(ic_loc, bvec)
     if (collision_field_model == 1) then
       ! cap_h_v not re-used else
       do iv=1,nv
          cap_h_v(ic_loc,iv) = bvec(iv)
        enddo
     endif
  enddo

  call timer_lib_out('coll')

  ! Compute the new phi
  if (collision_field_model == 1) then
     call cgyro_field_v
  endif

  call timer_lib_in('coll_comm')
  call parallel_lib_f_i_do(cap_h_ct)
  call timer_lib_out('coll_comm')

  call timer_lib_in('coll')

  ! Compute H given h and [phi(h), apar(h)]

!$omp parallel do private(iv_loc,is,ic,iv)
  do iv=nv1,nv2
     iv_loc = iv-nv1+1
     do ic=1,nc
        psi(ic,iv_loc) = sum(jvec_c(:,ic,iv_loc)*field(:,ic))
        chi(ic,iv_loc) = sum(jxvec_c(:,ic,iv_loc)*field(:,ic))
     enddo
     do ic=1,nc
        cap_h_c(ic,iv_loc) = cap_h_ct(iv_loc,ic)
     enddo
     is = is_v(iv)
     do ic=1,nc
        h_x(ic,iv_loc) = cap_h_c(ic,iv_loc)-psi(ic,iv_loc)*(z(is)/temp(is))
     enddo
  enddo

  call timer_lib_out('coll')

  if (collision_field_model == 0) then
     call cgyro_field_c
  endif

end subroutine cgyro_step_collision
