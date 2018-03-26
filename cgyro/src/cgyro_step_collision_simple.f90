!------------------------------------------------------------------------------
! cgyro_step_collision_simple.F90
!
! PURPOSE:
!  Take an implicit collision step using the pre-computed collision 
!  matrix.  Effectively, we compute the new collisional cap_H: 
!
!                       H = h + ze/T G phi
!------------------------------------------------------------------------------

subroutine cgyro_step_collision_simple

  use parallel_lib
  use timer_lib

  use cgyro_globals

  implicit none

  integer :: is,ie,ix,jx,it,ir
  integer :: ivp
  complex, dimension(:,:,:),allocatable :: bvec,cvec
  complex :: bvec_flat(nv)
  real :: cvec_re,cvec_im

  !----------------------------------------------------------------
  ! Perform data tranpose from _c to _v data layouts:
  call timer_lib_in('coll_comm')
  call parallel_lib_rtrans(cap_h_c,cap_h_v)
  call timer_lib_out('coll_comm')
  !----------------------------------------------------------------

  call timer_lib_in('coll')

  allocate(bvec(n_xi,n_energy,n_species))
  allocate(cvec(n_xi,n_energy,n_species))

!$omp parallel do private(ic_loc,ivp,iv,is,ix,jx,ie,ir,it,cvec_re,cvec_im,bvec,cvec,bvec_flat)
  do ic=nc1,nc2
     ic_loc = ic-nc1+1
     ir = ir_c(ic)
     it = it_c(ic)

     ! Set-up the RHS: H = f + ze/T G phi

     do iv=1,nv
        cvec(ix_v(iv),ie_v(iv),is_v(iv)) = cap_h_v(ic_loc,iv)
     enddo

     ! Avoid singularity of n=0,p=0:
     if (px(ir) == 0 .and. n == 0) then
        bvec = cvec
     else

        bvec = 0.0

        do is=1,n_species
           do ie=1,n_energy              
              do jx=1,n_xi

                 cvec_re = real(cvec(jx,ie,is))
                 cvec_im = aimag(cvec(jx,ie,is))

                 do ix=1,n_xi
                    bvec(ix,ie,is) = bvec(ix,ie,is)+ &
                         cmplx(cmat_simple(ix,jx,ie,is,it)*cvec_re, &
                         cmat_simple(ix,jx,ie,is,it)*cvec_im)
                 enddo
              enddo
           enddo
        enddo
     endif

     do iv=1,nv
        bvec_flat(iv) = bvec(ix_v(iv),ie_v(iv),is_v(iv))
     enddo
    call parallel_lib_f_i_set(ic_loc, bvec_flat) 

  enddo

  deallocate(bvec,cvec)

  call timer_lib_out('coll')

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

  call cgyro_field_c
  
end subroutine cgyro_step_collision_simple
