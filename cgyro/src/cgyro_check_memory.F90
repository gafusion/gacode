subroutine cgyro_check_memory(datafile)

  use cgyro_globals

  implicit none

  character (len=*), intent(in) :: datafile
  integer :: mult
  
  if (i_proc == 0) then

     open(unit=io,file=datafile,status='replace')

     total_memory = 0
     write(io,*) ' == Non-Distributed =='

     write(io,*)
     write(io,*) 'Pointers'
     write(io,*)    
     call cgyro_alloc_add(io,nv*4.0,'ie_v')
     call cgyro_alloc_add(io,nv*4.0,'ix_v')
     call cgyro_alloc_add(io,nv*4.0,'is_v')
     call cgyro_alloc_add(io,nc*4.0,'ir_c')
     call cgyro_alloc_add(io,nc*4.0,'it_c')
     call cgyro_alloc_add(io,n_radial*n_theta*4.0,'ic_c')
     call cgyro_alloc_add(io,n_energy*n_xi*n_species*4.0,'iv_v')

     write(io,*)
     write(io,*) 'Fields and field solve'
     write(io,*)
     call cgyro_alloc_add_3d(io,n_field,nc,nt_loc,16,'field')
     call cgyro_alloc_add_3d(io,n_field,nc,nt_loc,16,'field_dot')
     call cgyro_alloc_add_3d(io,n_field,nc,nt_loc,16,'field_loc')
     call cgyro_alloc_add_3d(io,n_field,nc,nt_loc,16,'field_old')
     call cgyro_alloc_add_3d(io,n_field,nc,nt_loc,16,'field_old2')
     call cgyro_alloc_add_3d(io,n_field,nc,nt_loc,16,'field_old3')
     call cgyro_alloc_add_3d(io,n_field,nc,nt_loc,8,'fcoef')
     if (n_field < 3) then
        call cgyro_alloc_add_3d(io,n_field,nc,nt_loc,8,'gcoef')
     else
        call cgyro_alloc_add_3d(io,5,nc,nt_loc,8,'gcoef')
     endif
     if ((collision_model /= 5) .AND. (collision_field_model == 1)) then
       call cgyro_alloc_add_3d(io,n_field,nt_loc,nc,16,'field_v')
       call cgyro_alloc_add_3d(io,n_field,nt_loc,nc_loc,16,'field_loc_v')
       call cgyro_alloc_add_4d(io,n_field,nv,nt_loc,nc_loc,8,'dvjvec_v')
     endif
     call cgyro_alloc_add_4d(io,n_field,nc,nv_loc,nt_loc,8,'dvjvec_c')

     if (nonlinear_flag == 1) then
        write(io,*)
        write(io,*) 'Nonlinear'
        write(io,*)
        ! nsplit * n_toroidal = nv_loc * n_theta
        if (nl_single_flag > 1) then
          ! single precision nl
#if !(defined(OMPGPU) || defined(_OPENACC))
        call cgyro_alloc_add_s_3d(io,(ny/2+1),nx,n_omp,8,'fx')
        call cgyro_alloc_add_s_3d(io,(ny/2+1),nx,n_omp,8,'gx')
        call cgyro_alloc_add_s_3d(io,(ny/2+1),nx,n_omp,8,'fy')
        call cgyro_alloc_add_s_3d(io,(ny/2+1),nx,n_omp,8,'gy')
        call cgyro_alloc_add_s_3d(io,ny,nx,nsplitA,4,'ux')
        call cgyro_alloc_add_s_3d(io,ny,nx,nsplit,4,'vx')
        call cgyro_alloc_add_s_3d(io,ny,nx,nsplitA,4,'uy')
        call cgyro_alloc_add_s_3d(io,ny,nx,nsplit,4,'vy')
        call cgyro_alloc_add_s_3d(io,ny,nx,n_omp,4,'uv')
#else
        call cgyro_alloc_add_s_3d(io,(ny/2+1),nx,nsplitA,8,'fx')
        call cgyro_alloc_add_s_3d(io,(ny/2+1),nx,nsplit,8,'gx')
        call cgyro_alloc_add_s_3d(io,(ny/2+1),nx,nsplitA,8,'fy')
        call cgyro_alloc_add_s_3d(io,(ny/2+1),nx,nsplit,8,'gy')
        call cgyro_alloc_add_s_3d(io,ny,nx,nsplitA,4,'ux')
        call cgyro_alloc_add_s_3d(io,ny,nx,nsplit,4,'vx')
        call cgyro_alloc_add_s_3d(io,ny,nx,nsplitA,4,'uy')
        call cgyro_alloc_add_s_3d(io,ny,nx,nsplit,4,'vy')
        call cgyro_alloc_add_s_3d(io,ny,nx,nsplitA,4,'uv')
#endif
        else
#if !(defined(OMPGPU) || defined(_OPENACC))
        call cgyro_alloc_add_3d(io,(ny/2+1),nx,n_omp,16,'fx')
        call cgyro_alloc_add_3d(io,(ny/2+1),nx,n_omp,16,'gx')
        call cgyro_alloc_add_3d(io,(ny/2+1),nx,n_omp,16,'fy')
        call cgyro_alloc_add_3d(io,(ny/2+1),nx,n_omp,16,'gy')
        call cgyro_alloc_add_3d(io,ny,nx,nsplitA,8,'ux')
        call cgyro_alloc_add_3d(io,ny,nx,nsplit,8,'vx')
        call cgyro_alloc_add_3d(io,ny,nx,nsplitA,8,'uy')
        call cgyro_alloc_add_3d(io,ny,nx,nsplit,8,'vy')
        call cgyro_alloc_add_3d(io,ny,nx,n_omp,8,'uv')
#else
        call cgyro_alloc_add_3d(io,(ny/2+1),nx,nsplitA,16,'fx')
        call cgyro_alloc_add_3d(io,(ny/2+1),nx,nsplit,16,'gx')
        call cgyro_alloc_add_3d(io,(ny/2+1),nx,nsplitA,16,'fy')
        call cgyro_alloc_add_3d(io,(ny/2+1),nx,nsplit,16,'gy')
        call cgyro_alloc_add_3d(io,ny,nx,nsplitA,8,'ux')
        call cgyro_alloc_add_3d(io,ny,nx,nsplit,8,'vx')
        call cgyro_alloc_add_3d(io,ny,nx,nsplitA,8,'uy')
        call cgyro_alloc_add_3d(io,ny,nx,nsplit,8,'vy')
        call cgyro_alloc_add_3d(io,ny,nx,nsplitA,8,'uv')
#endif
        endif
     endif

     write(io,*)
     write(io,*) 'TOTAL'
     write(io,10) total_memory/1e6,' MB [persistent memory per MPI process]'
     write(io,*)

     total_memory = 0
     write(io,*) ' == Distributed =='

     write(io,*)
     write(io,*) 'Distribution-like arrays'
     write(io,*)

     select case(delta_t_method)
     case(1)
        call cgyro_alloc_add_3d(io,nc,nv_loc,nt_loc,16,'h0_old')
        call cgyro_alloc_add_4d(io,nc,nv_loc,nt_loc,6,16,'rhs')
     case(2)
        call cgyro_alloc_add_3d(io,nc,nv_loc,nt_loc,16,'h0_old')
        call cgyro_alloc_add_4d(io,nc,nv_loc,nt_loc,7,16,'rhs')
     case(3)
        call cgyro_alloc_add_3d(io,nc,nv_loc,nt_loc,16,'h0_old')
        call cgyro_alloc_add_4d(io,nc,nv_loc,nt_loc,9,16,'rhs')
     case default
        ! Normal timestep
        call cgyro_alloc_add_4d(io,nc,nv_loc,nt_loc,4,16,'rhs')
     end select

     call cgyro_alloc_add_4d(io,n_field,nc,nv_loc,nt_loc,16,'omega_s')
     call cgyro_alloc_add_4d(io,n_field,nc,nv_loc,nt_loc,16,'omega_ss')
     call cgyro_alloc_add_3d(io,nc,nv_loc,nt_loc,16,'omega_cap_h')
     call cgyro_alloc_add_3d(io,nc,nv_loc,nt_loc,16,'omega_h')
     if (triad_print_flag == 1) then
       call cgyro_alloc_add_3d(io,nc,nv_loc,nt_loc,16,'diss_r')
     endif
     call cgyro_alloc_add_3d(io,nc,nv_loc,nt_loc,16,'h_x')
     call cgyro_alloc_add_3d(io,nc,nv_loc,nt_loc,16,'h0_x')
     call cgyro_alloc_add_3d(io,nc,nv_loc,nt_loc,16,'g_x')
     call cgyro_alloc_add_3d(io,nc,nv_loc,nt_loc,16,'cap_h_c')
     call cgyro_alloc_add_3d(io,nc,nv_loc,nt_loc,16,'cap_h_ct')
     call cgyro_alloc_add_3d(io,nc,nv_loc,nt_loc,16,'cap_h_c_dot')
     call cgyro_alloc_add_3d(io,nc,nv_loc,nt_loc,16,'cap_h_c_old')
     call cgyro_alloc_add_3d(io,nc,nv_loc,nt_loc,16,'cap_h_c_old2')
     call cgyro_alloc_add_3d(io,nc_loc,nv,nt_loc,16,'cap_h_v')
     call cgyro_alloc_add_4d(io,n_field,nc,nv_loc,nt_loc,8,'jvec_c')
     if (nonlinear_flag == 1) call cgyro_alloc_add(io,n_field*n_radial*n_jtheta*nv_loc*n_toroidal*8.0,'jvec_c_nl')
     call cgyro_alloc_add_4d(io,n_field,nc_loc,nv,nt_loc,8,'jvec_v')
     call cgyro_alloc_add_4d(io,n_field,nc,nv_loc,nt_loc,8,'jxvec_c')
     call cgyro_alloc_add_3d(io,nc,nv_loc,nt_loc,8,'upfac1')
     call cgyro_alloc_add_3d(io,nc,nv_loc,nt_loc,8,'upfac2')

     if (nonlinear_flag == 1) then
        write(io,*)
        write(io,*) 'Nonlinear bracket'
        write(io,*)
        ! nsplit * n_toroidal = nv_loc * n_theta
        if (nl_single_flag > 1) then
          ! single precision nl
        call cgyro_alloc_add_s_4d(io,n_radial,nt_loc,nsplitA,n_toroidal_procs,8,'fA_nl')
        if (nsplitB > 0) call cgyro_alloc_add_s_4d(io,n_radial,nt_loc,nsplitB,n_toroidal_procs,8,'fB_nl')
        call cgyro_alloc_add_s_4d(io,n_field,n_radial,n_jtheta,n_toroidal,8,'g_nl')
        call cgyro_alloc_add_s_3d(io,n_radial,nt_loc,nsplitA*n_toroidal_procs,8,'fpackA')
        if (nsplitB > 0) call cgyro_alloc_add_s_3d(io,n_radial,nt_loc,nsplitB*n_toroidal_procs,8,'fpackB')
        call cgyro_alloc_add_s_4d(io,n_field,n_radial,n_jtheta,n_toroidal,8,'gpack')
        else
                ! double precision nl
        call cgyro_alloc_add_4d(io,n_radial,nt_loc,nsplitA,n_toroidal_procs,16,'fA_nl')
        if (nl_single_flag > 0) call cgyro_alloc_add_s_4d(io,n_radial,nt_loc,nsplitA,n_toroidal_procs,8,'fA_nl32')
        if (nsplitB > 0) then ! nsplitB can be zero at large MPI
        call cgyro_alloc_add_4d(io,n_radial,nt_loc,nsplitB,n_toroidal_procs,16,'fB_nl')
        if (nl_single_flag > 0) call cgyro_alloc_add_s_4d(io,n_radial,nt_loc,nsplitB,n_toroidal_procs,16,'fB_nl32')
        endif
        call cgyro_alloc_add_4d(io,n_field,n_radial,n_jtheta,n_toroidal,16,'g_nl')
        call cgyro_alloc_add_3d(io,n_radial,nt_loc,nsplitA*n_toroidal_procs,16,'fpackA')
        if (nl_single_flag > 0) call cgyro_alloc_add_s_3d(io,n_radial,nt_loc,nsplitA*n_toroidal_procs,8,'fpackA32')
        if (nsplitB > 0) then ! nsplitB can be zero at large MPI
        call cgyro_alloc_add_3d(io,n_radial,nt_loc,nsplitB*n_toroidal_procs,16,'fpackB')
        if (nl_single_flag > 0) call cgyro_alloc_add_s_3d(io,n_radial,nt_loc,nsplitB*n_toroidal_procs,16,'fpackB32')
        endif
        call cgyro_alloc_add_4d(io,n_field,n_radial,n_jtheta,n_toroidal,16,'gpack')
        endif
        if (triad_print_flag == 1) then ! always double precision
          call cgyro_alloc_add_4d(io,n_radial,nt_loc,nsplitA,n_toroidal_procs,16,'eA_nl')
          call cgyro_alloc_add_4d(io,n_radial,nt_loc,nsplitB,n_toroidal_procs,16,'eB_nl')
          call cgyro_alloc_add_3d(io,n_radial,nt_loc,nsplitA*n_toroidal_procs,16,'epackA')
          call cgyro_alloc_add_3d(io,n_radial,nt_loc,nsplitB*n_toroidal_procs,16,'epackB')
        endif
     endif

     write(io,*)
     write(io,*) 'Collision operator'
     write(io,*)

     if(collision_model == 5) then
        call cgyro_alloc_add(io,(8.0*n_xi)*n_xi*n_species*n_energy*n_theta*nt_loc,'cmat')
     else
        if (collision_precision_mode == 1) then
           call cgyro_alloc_add_s_4d(io,nv,nv,nc_loc,nt_loc,4,'cmat_fp32')
           call cgyro_alloc_add(io,4.0*n_xi*n_species*(n_energy-n_low_energy)*n_xi*nc_loc*nt_loc,'cmat_stripes')
           call cgyro_alloc_add(io,4.0*n_xi*n_species*nv*nc_loc*n_low_energy*nt_loc,'cmat_e1')
        else if (collision_precision_mode == 32) then
           call cgyro_alloc_add_s_4d(io,nv,nv,nc_loc,nt_loc,4,'cmat_fp32')
        else
           call cgyro_alloc_add_4d(io,nv,nv,nc_loc,nt_loc,8,'cmat')
        endif
#if defined(OMPGPU) || defined(_OPENACC)
        if (gpu_bigmem_flag /= 1) then
           write(io,*) 'Note: cmat is not in GPU memory'
        endif
#endif
     endif

     write(io,*)
     write(io,*) 'TOTAL'
     if (test_flag == 1) then
        write(io,10) total_memory/1e9,&
             ' GB [per toroidal mode ; halved with every doubling of MPI processes]'
        write(io,*) ' '
        do mult=2,16,2
           write(io,20) total_memory/1e9/mult,' GB ',mult*n_toroidal,' MPI processes'
        enddo           
     else
        write(io,10) total_memory/1e9,' GB [per MPI process]'
     endif    
     close(io)

  end if

10 format(t2,f8.3,a,3x,a)
20 format(t2,f8.3,a,i5,a)

end subroutine cgyro_check_memory

!------------------------------------------------
! cgyro_alloc_add.f90
!
! PURPOSE:
!  Primitive allocation addition routine.
!------------------------------------------------

subroutine cgyro_alloc_add(my_io,bytes,name)

  use cgyro_globals

  implicit none
  !
  integer, intent(in) :: my_io
  real, intent(in) :: bytes
  character (len=*), intent(in) :: name
 
  if (i_proc == 0) then
     
     total_memory = total_memory+bytes

     if (bytes < 1e8) then 
        write(my_io,10) bytes/1e6,' MB',name
     else
        write(my_io,10) bytes/1e9,' GB',name
     endif
     
  endif

10 format(t2,f8.3,a,3x,a)
  
end subroutine cgyro_alloc_add

subroutine cgyro_alloc_add_3d(my_io,d1,d2,d3,elsize,name)

  use cgyro_globals

  implicit none
  !
  integer, intent(in) :: my_io
  integer, intent(in) :: d1,d2,d3,elsize
  character (len=*), intent(in) :: name
  !
  real :: bytes
  character (15) :: name15
  character (14) :: typestr
  character(len=50) :: fstr
 
  if (i_proc == 0) then

     bytes = elsize
     bytes = ((bytes*d1)*d2)*d3
     
     total_memory = total_memory+bytes

     name15 = name
     if ( (d1<=9999) .and. (d2<=9999) .and.(d3<=9999) ) then
       fstr = '(t2,f8.3,a,3x,a11,3x,a,i4,a,i4,a,i4,a,a)'
     else if ( (d1<=99999) .and. (d2<=99999) .and.(d3<=99999)) then
       fstr = '(t2,f8.3,a,3x,a11,3x,a,i5,a,i5,a,i5,a,a)'
     else
       fstr = '(t2,f8.3,a,3x,a11,3x,a,i7,a,i7,a,i7,a,a)'
     endif

     if (elsize<5) then
        typestr = 'integer'
     else if (elsize<9) then
        typestr = 'real'
     else
        typestr = 'complex'
     endif

     if (bytes < 1e8) then 
        write(my_io,fstr) bytes/1e6,' MB',name15,' (',d1,',',d2,',',d3,') ',typestr
     else
        write(my_io,fstr) bytes/1e9,' GB',name15,' (',d1,',',d2,',',d3,') ',typestr
     endif
     
  endif
  
end subroutine cgyro_alloc_add_3d

subroutine cgyro_alloc_add_4d(my_io,d1,d2,d3,d4,elsize,name)

  use cgyro_globals

  implicit none
  !
  integer, intent(in) :: my_io
  integer, intent(in) :: d1,d2,d3,d4,elsize
  character (len=*), intent(in) :: name
  !
  real :: bytes
  character (15) :: name15
  character (14) :: typestr
  character(len=50) :: fstr
 
  if (i_proc == 0) then

     bytes = elsize
     bytes = (((bytes*d1)*d2)*d3)*d4
     
     total_memory = total_memory+bytes

     name15 = name
     if ( (d1<=9999) .and. (d2<=9999) .and.(d3<=9999) .and.(d4<=9999) ) then
        fstr = '(t2,f8.3,a,3x,a11,3x,a,i4,a,i4,a,i4,a,i4,a,a)'
     else if ( (d1<=99999) .and. (d2<=99999) .and.(d3<=99999) .and.(d4<=99999)) then
        fstr = '(t2,f8.3,a,3x,a11,3x,a,i5,a,i5,a,i5,a,i5,a,a)'
     else
        fstr = '(t2,f8.3,a,3x,a11,3x,a,i7,a,i7,a,i7,a,i7,a,a)'
     endif

     if (elsize<5) then
        typestr = 'integer'
     else if (elsize<9) then
        typestr = 'real'
     else
        typestr = 'complex'
     endif

     if (bytes < 1e8) then 
        write(my_io,fstr) bytes/1e6,' MB',name15,' (',d1,',',d2,',',d3,',',d4,') ',typestr
     else
        write(my_io,fstr) bytes/1e9,' GB',name15,' (',d1,',',d2,',',d3,',',d4,') ',typestr
     endif
     
  endif

end subroutine cgyro_alloc_add_4d

subroutine cgyro_alloc_add_s_3d(my_io,d1,d2,d3,elsize,name)

  use cgyro_globals

  implicit none
  !
  integer, intent(in) :: my_io
  integer, intent(in) :: d1,d2,d3,elsize
  character (len=*), intent(in) :: name
  !
  real :: bytes
  character (15) :: name15
  character (14) :: typestr
  character(len=50) :: fstr
 
  if (i_proc == 0) then

     bytes = elsize
     bytes = ((bytes*d1)*d2)*d3
     
     total_memory = total_memory+bytes

     name15 = name
     if ( (d1<=9999) .and. (d2<=9999) .and.(d3<=9999) ) then
       fstr = '(t2,f8.3,a,3x,a11,3x,a,i4,a,i4,a,i4,a,a)'
     else if ( (d1<=99999) .and. (d2<=99999) .and.(d3<=99999)) then
       fstr = '(t2,f8.3,a,3x,a11,3x,a,i5,a,i5,a,i5,a,a)'
     else
       fstr = '(t2,f8.3,a,3x,a11,3x,a,i7,a,i7,a,i7,a,a)'
     endif

     if (elsize<5) then
        typestr = 'real(fp32)'
     else
        typestr = 'complex(fp32)'
     endif

     if (bytes < 1e8) then 
        write(my_io,fstr) bytes/1e6,' MB',name15,' (',d1,',',d2,',',d3,') ',typestr
     else
        write(my_io,fstr) bytes/1e9,' GB',name15,' (',d1,',',d2,',',d3,') ',typestr
     endif
     
  endif
  
end subroutine cgyro_alloc_add_s_3d

subroutine cgyro_alloc_add_s_4d(my_io,d1,d2,d3,d4,elsize,name)

  use cgyro_globals

  implicit none
  !
  integer, intent(in) :: my_io
  integer, intent(in) :: d1,d2,d3,d4,elsize
  character (len=*), intent(in) :: name
  !
  real :: bytes
  character (15) :: name15
  character (14) :: typestr
  character(len=50) :: fstr
 
  if (i_proc == 0) then

     bytes = elsize
     bytes = (((bytes*d1)*d2)*d3)*d4
     
     total_memory = total_memory+bytes

     name15 = name
     if ( (d1<=9999) .and. (d2<=9999) .and.(d3<=9999) .and.(d4<=9999) ) then
        fstr = '(t2,f8.3,a,3x,a11,3x,a,i4,a,i4,a,i4,a,i4,a,a)'
     else if ( (d1<=99999) .and. (d2<=99999) .and.(d3<=99999) .and.(d4<=99999)) then
        fstr = '(t2,f8.3,a,3x,a11,3x,a,i5,a,i5,a,i5,a,i5,a,a)'
     else
        fstr = '(t2,f8.3,a,3x,a11,3x,a,i7,a,i7,a,i7,a,i7,a,a)'
     endif

     if (elsize<5) then
        typestr = 'real(fp32)'
     else
        typestr = 'complex(fp32)'
     endif

     if (bytes < 1e8) then 
        write(my_io,fstr) bytes/1e6,' MB',name15,' (',d1,',',d2,',',d3,',',d4,') ',typestr
     else
        write(my_io,fstr) bytes/1e9,' GB',name15,' (',d1,',',d2,',',d3,',',d4,') ',typestr
     endif
     
  endif

end subroutine cgyro_alloc_add_s_4d
