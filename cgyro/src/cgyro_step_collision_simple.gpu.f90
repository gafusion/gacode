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

  integer :: is,ie,ix,jx,it,ir,j,k
  integer :: ivp,nj_loc
  complex, dimension(:,:,:),allocatable :: bvec,cvec
  complex :: bvec_flat(nv)
  real :: cvec_re,cvec_im

  !----------------------------------------------------------------
  ! Perform data tranpose from _c to _v data layouts:
  call timer_lib_in('coll_mem')
  call parallel_lib_rtrans_pack_gpu(cap_h_c)
  call timer_lib_out('coll_mem')
  call timer_lib_in('coll_comm')
  call parallel_lib_r_do_gpu(cap_h_v)
  call timer_lib_out('coll_comm')
  !----------------------------------------------------------------

  call timer_lib_in('coll')

  allocate(bvec(n_xi,n_energy,n_species))
  allocate(cvec(n_xi,n_energy,n_species))

  call parallel_lib_nj_loc(nj_loc)

!$acc parallel loop gang num_workers(4) vector_length(32) &
!$acc&         present(ix_v,ie_v,is_v,ir_c,it_c,px,cap_h_v,cmat_simple,fsendf) &
!$acc&         private(bvec,cvec,bvec_flat) &
!$acc&         private(ic_loc,ir,it) default(none)
  do ic=nc1,nc2

     ic_loc = ic-nc1+1
     ir = ir_c(ic)
     it = it_c(ic)

     ! Set-up the RHS: H = f + ze/T G phi

!$acc loop worker vector
     do iv=1,nv
        cvec(ix_v(iv),ie_v(iv),is_v(iv)) = cap_h_v(ic_loc,iv)
     enddo

     ! Avoid singularity of n=0,p=0:
     if (px(ir) == 0 .and. n == 0) then

!$acc loop worker collapse(2)
        do is=1,n_species
           do ie=1,n_energy              
!$acc loop vector
             do ix=1,n_xi
               bvec(ix,ie,is) = cvec(ix,ie,is)
             enddo
           enddo
        enddo
     else

!$acc loop worker collapse(2)
        do is=1,n_species
           do ie=1,n_energy              

!$acc loop vector
              do ix=1,n_xi
                 bvec(ix,ie,is) = 0.0
              enddo

!$acc loop seq
              do jx=1,n_xi

                 cvec_re = real(cvec(jx,ie,is))
                 cvec_im = aimag(cvec(jx,ie,is))

!$acc loop vector
                 do ix=1,n_xi
                    bvec(ix,ie,is) = bvec(ix,ie,is)+ &
                         cmplx(cmat_simple(ix,jx,ie,is,it)*cvec_re, &
                         cmat_simple(ix,jx,ie,is,it)*cvec_im)
                 enddo
              enddo
           enddo
        enddo
     endif

!$acc loop worker vector
     do iv=1,nv
        bvec_flat(iv) = bvec(ix_v(iv),ie_v(iv),is_v(iv))
     enddo

!$acc loop worker
     do k=1,nproc
!$acc loop vector
        do j=1,nj_loc
           fsendf(j,ic_loc,k) = bvec_flat(j+(k-1)*nj_loc)
        enddo
     enddo

  enddo

  deallocate(bvec,cvec)

  call timer_lib_out('coll')

  call timer_lib_in('coll_comm')
  call parallel_lib_f_i_do_gpu(cap_h_ct)
  call timer_lib_out('coll_comm')

  call timer_lib_in('coll')

  ! Compute H given h and [phi(h), apar(h)]

!$acc parallel loop collapse(2) gang vector private(iv_loc,is) &
!$acc&         present(is_v,psi,cap_h_c,cap_h_ct,cap_h_c,jvec_c,field,z,temp,h_x) &
!$acc&         default(none)
  do iv=nv1,nv2
     do ic=1,nc
        iv_loc = iv-nv1+1
        is = is_v(iv)
        psi(ic,iv_loc) = sum(jvec_c(:,ic,iv_loc)*field(:,ic))
        cap_h_c(ic,iv_loc) = cap_h_ct(iv_loc,ic)
        h_x(ic,iv_loc) = cap_h_c(ic,iv_loc)-psi(ic,iv_loc)*(z(is)/temp(is))
     enddo
  enddo


  call timer_lib_out('coll')

  call cgyro_field_c_gpu

end subroutine cgyro_step_collision_simple
