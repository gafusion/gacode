! This is called after cgyro_kernel to perform clean-ups from CPU and GPU (deallocate)

subroutine cgyro_cleanup
  use cgyro_globals
  use cgyro_coll_data, only : cgyro_coll_data_cleanup
  use parallel_lib
  use cgyro_field_mod, only : cgyro_field_c_cleanup, cgyro_field_e_cleanup, &
       cgyro_field_v_cleanup
  use cgyro_flux_mod

#if defined(_OPENACC) || defined(OMPGPU)
#define CGYRO_GPU_FFT
#endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!
  ! From cgyro_init_manager
  !!!!!!!!!!!!!!!!!!!!!!!!!!
  
  if(allocated(energy))        deallocate(energy)
  if(allocated(vel))        then
#if defined(OMPGPU)
!$omp target exit data map(release:vel)
#elif defined(_OPENACC)
!$acc exit data delete(vel)
#endif
     deallocate(vel)
  endif
  if(allocated(vel2))        then
#if defined(OMPGPU)
!$omp target exit data map(release:vel2)
#elif defined(_OPENACC)
!$acc exit data delete(vel2)
#endif
     deallocate(vel2)
  endif
  if(allocated(w_e))           deallocate(w_e)
  if(allocated(e_deriv1_mat))  deallocate(e_deriv1_mat)
  if(allocated(e_deriv1_rot_mat))  deallocate(e_deriv1_rot_mat)
  if(allocated(xi))            then
#if defined(OMPGPU)
!$omp target exit data map(release:xi)
#elif defined(_OPENACC)
!$acc exit data delete(xi)
#endif
     deallocate(xi)
  endif
  if(allocated(w_xi))          deallocate(w_xi)
  if(allocated(w_exi))         deallocate(w_exi)
  if(allocated(xi_lor_mat))    deallocate(xi_lor_mat)
  if(allocated(xi_deriv_mat))  deallocate(xi_deriv_mat)
  
  if(allocated(theta))          deallocate(theta)
  if(allocated(thetab))         deallocate(thetab)
  if(allocated(w_theta))        deallocate(w_theta)
  if(allocated(g_theta))        deallocate(g_theta)
  if(allocated(g_theta_geo))    deallocate(g_theta_geo)
  if(allocated(bmag))           deallocate(bmag)
  if(allocated(btor))           deallocate(btor)
  if(allocated(bpol))           deallocate(bpol)
  if(allocated(k_perp))         deallocate(k_perp)
  if(allocated(k_x))            deallocate(k_x)
  if(allocated(bigr))           deallocate(bigr)
  if(allocated(bigr_r))         deallocate(bigr_r)
  if(allocated(captheta))       deallocate(captheta)
  if(allocated(itp))            deallocate(itp)
  if(allocated(omega_stream))   then
#if defined(OMPGPU)
!$omp target exit data map(release:omega_stream)
#elif defined(_OPENACC)
!$acc exit data delete(omega_stream)
#endif
     deallocate(omega_stream)
  endif
  if(allocated(omega_trap))          deallocate(omega_trap)
  if(allocated(omega_rdrift))        deallocate(omega_rdrift)
  if(allocated(omega_adrift))        deallocate(omega_adrift)
  if(allocated(omega_aprdrift))      deallocate(omega_aprdrift)
  if(allocated(omega_cdrift))        deallocate(omega_cdrift)
  if(allocated(omega_cdrift_r))      deallocate(omega_cdrift_r)
  if(allocated(omega_gammap))        deallocate(omega_gammap)
  if(allocated(lambda_rot))          deallocate(lambda_rot)
  if(allocated(dlambda_rot))         deallocate(dlambda_rot)
  if(allocated(dens_rot))            deallocate(dens_rot)
  if(allocated(dens2_rot))           deallocate(dens2_rot)
  if(allocated(dens_ele_rot))        deallocate(dens_ele_rot)
  if(allocated(dens_avg_rot))        deallocate(dens_avg_rot)
  if(allocated(dlnndr_avg_rot))      deallocate(dlnndr_avg_rot)
  if(allocated(omega_rot_trap))      deallocate(omega_rot_trap)
  if(allocated(omega_rot_u))         deallocate(omega_rot_u)
  if(allocated(omega_rot_drift))     deallocate(omega_rot_drift)
  if(allocated(omega_rot_drift_r))   deallocate(omega_rot_drift_r)
  if(allocated(omega_rot_edrift))    deallocate(omega_rot_edrift)
  if(allocated(omega_rot_edrift_r))  deallocate(omega_rot_edrift_r)
  if(allocated(omega_rot_star))      deallocate(omega_rot_star)
  if(allocated(gtime))               deallocate(gtime)
  if(allocated(freq))                deallocate(freq)
  call cgyro_field_e_cleanup
  call cgyro_field_v_cleanup
  call cgyro_field_c_cleanup
  call cgyro_flux_cleanup
  if(allocated(epar))                deallocate(epar)
  if(allocated(recv_status))         deallocate(recv_status)
  if(allocated(source)) then
#if defined(OMPGPU)
!$omp target exit data map(release:source)
#elif defined(_OPENACC)
!$acc exit data delete(source)
#endif
     deallocate(source)
  endif
  if(allocated(thfac_itor)) then
#if defined(OMPGPU)
!$omp target exit data map(release:thfac_itor)
#elif defined(_OPENACC)
!$acc exit data delete(thfac_itor)
#endif
     deallocate(thfac_itor)
  endif
  if(allocated(h0_old)) then
#if defined(OMPGPU)
!$omp target exit data map(release:h0_old)
#elif defined(_OPENACC)
!$acc exit data delete(h0_old)
#endif
     deallocate(h0_old)
  endif
  if(allocated(rhs)) then
#if defined(OMPGPU)
!$omp target exit data map(release:rhs)
#elif defined(_OPENACC)
!$acc exit data delete(rhs)
#endif
     deallocate(rhs)
  endif
  if(allocated(h_x)) then
#if defined(OMPGPU)
!$omp target exit data map(release:h_x)
#elif defined(_OPENACC)
!$acc exit data delete(h_x)
#endif
     deallocate(h_x)
  endif
  if(allocated(g_x)) then
#if defined(OMPGPU)
!$omp target exit data map(release:g_x)
#elif defined(_OPENACC)
!$acc exit data delete(g_x)
#endif
     deallocate(g_x)
  endif
  if(allocated(h0_x)) then
#if defined(OMPGPU)
!$omp target exit data map(release:h0_x)
#elif defined(_OPENACC)
!$acc exit data delete(h0_x)
#endif
     deallocate(h0_x)
  endif
  if(allocated(cap_h_c)) then
#if defined(OMPGPU)
!$omp target exit data map(release:cap_h_c)
#elif defined(_OPENACC)
!$acc exit data delete(cap_h_c)
#endif
     deallocate(cap_h_c)
  endif
  if(allocated(cap_h_ct))  then
#if defined(OMPGPU)
!$omp target exit data map(release:cap_h_ct)
#elif defined(_OPENACC)
!$acc exit data delete(cap_h_ct)
#endif
     deallocate(cap_h_ct)
  endif
  if(allocated(omega_cap_h)) then
#if defined(OMPGPU)
!$omp target exit data map(release:omega_cap_h)
#elif defined(_OPENACC)
!$acc exit data delete(omega_cap_h)
#endif
     deallocate(omega_cap_h)
  endif
  if(allocated(omega_h)) then
#if defined(OMPGPU)
!$omp target exit data map(release:omega_h)
#elif defined(_OPENACC)
!$acc exit data delete(omega_h)
#endif
     deallocate(omega_h)
  endif
  if(allocated(omega_s)) then
#if defined(OMPGPU)
!$omp target exit data map(release:omega_s)
#elif defined(_OPENACC)
!$acc exit data delete(omega_s)
#endif
     deallocate(omega_s)
  endif
  if(allocated(omega_ss)) then
#if defined(OMPGPU)
!$omp target exit data map(release:omega_ss)
#elif defined(_OPENACC)
!$acc exit data delete(omega_ss)
#endif
     deallocate(omega_ss)
  endif
  if(allocated(omega_sbeta)) then
#if defined(OMPGPU)
!$omp target exit data map(release:omega_sbeta)
#elif defined(_OPENACC)
!$acc exit data delete(omega_sbeta)
#endif
     deallocate(omega_sbeta)
  endif
  if(allocated(jvec_c))  then
#if defined(OMPGPU)
!$omp target exit data map(release:jvec_c)
#elif defined(_OPENACC)
!$acc exit data delete(jvec_c)
#endif
     deallocate(jvec_c)
  endif
  if(allocated(jvec_c_nl))  then
#if defined(OMPGPU)
!$omp target exit data map(release:jvec_c_nl)
#elif defined(_OPENACC)
!$acc exit data delete(jvec_c_nl)
#endif
     deallocate(jvec_c_nl)
  endif
  if(allocated(jvec_v))              deallocate(jvec_v)
  if(allocated(jxvec_c))  then
     deallocate(jxvec_c)
  endif
  if(allocated(upfac1))   then
#if defined(OMPGPU)
!$omp target exit data map(release:upfac1)
#elif defined(_OPENACC)
!$acc exit data delete(upfac1)
#endif
     deallocate(upfac1)
  endif
  if(allocated(upfac2))   then
#if defined(OMPGPU)
!$omp target exit data map(release:upfac2)
#elif defined(_OPENACC)
!$acc exit data delete(upfac2)
#endif
     deallocate(upfac2)
  endif
  if(allocated(cap_h_v))  then
#if defined(OMPGPU)
!$omp target exit data map(release:cap_h_v)
#elif defined(_OPENACC)
!$acc exit data delete(cap_h_v)
#endif
     deallocate(cap_h_v)
  endif
  if(allocated(upwind_res_loc))   then
#if defined(OMPGPU)
!$omp target exit data map(release:upwind_res_loc)
#elif defined(_OPENACC)
!$acc exit data delete(upwind_res_loc)
#endif
     deallocate(upwind_res_loc)
  endif
  if(allocated(upwind_res))   then
#if defined(OMPGPU)
!$omp target exit data map(release:upwind_res)
#elif defined(_OPENACC)
!$acc exit data delete(upwind_res)
#endif
     deallocate(upwind_res)
  endif
  if(allocated(upwind32_res_loc))   then
#if defined(OMPGPU)
!$omp target exit data map(release:upwind32_res_loc)
#elif defined(_OPENACC)
!$acc exit data delete(upwind32_res_loc)
#endif
     deallocate(upwind32_res_loc)
  endif
  if(allocated(upwind32_res))   then
#if defined(OMPGPU)
!$omp target exit data map(release:upwind32_res)
#elif defined(_OPENACC)
!$acc exit data delete(upwind32_res)
#endif
     deallocate(upwind32_res)
  endif
  if(allocated(fA_nl))   then
#if defined(OMPGPU)
!$omp target exit data map(release:fA_nl)
#elif defined(_OPENACC)
!$acc exit data delete(fA_nl)
#endif
     deallocate(fA_nl)
  endif
  if(allocated(fB_nl))   then
#if defined(OMPGPU)
!$omp target exit data map(release:fB_nl)
#elif defined(_OPENACC)
!$acc exit data delete(fB_nl)
#endif
     deallocate(fB_nl)
  endif
  if(allocated(fA_nl32))   then
#if defined(OMPGPU)
!$omp target exit data map(release:fA_nl32)
#elif defined(_OPENACC)
!$acc exit data delete(fA_nl32)
#endif
     deallocate(fA_nl32)
  endif
  if(allocated(fB_nl32))   then
#if defined(OMPGPU)
!$omp target exit data map(release:fB_nl32)
#elif defined(_OPENACC)
!$acc exit data delete(fB_nl32)
#endif
     deallocate(fB_nl32)
  endif
  if(allocated(g_nl))   then
#if defined(OMPGPU)
!$omp target exit data map(release:g_nl)
#elif defined(_OPENACC)
!$acc exit data delete(g_nl)
#endif
     deallocate(g_nl)
  endif
  if(allocated(g_nl32))   then
#if defined(OMPGPU)
!$omp target exit data map(release:g_nl32)
#elif defined(_OPENACC)
!$acc exit data delete(g_nl32)
#endif
     deallocate(g_nl32)
  endif
  if(allocated(fpackA))   then
#if defined(OMPGPU)
!$omp target exit data map(release:fpackA)
#elif defined(_OPENACC)
!$acc exit data delete(fpackA)
#endif
     deallocate(fpackA)
  endif
  if(allocated(fpackB))   then
#if defined(OMPGPU)
!$omp target exit data map(release:fpackB)
#elif defined(_OPENACC)
!$acc exit data delete(fpackB)
#endif
     deallocate(fpackB)
  endif
  if(allocated(fpackA32))   then
#if defined(OMPGPU)
!$omp target exit data map(release:fpackA32)
#elif defined(_OPENACC)
!$acc exit data delete(fpackA32)
#endif
     deallocate(fpackA32)
  endif
  if(allocated(fpackB32))   then
#if defined(OMPGPU)
!$omp target exit data map(release:fpackB32)
#elif defined(_OPENACC)
!$acc exit data delete(fpackB32)
#endif
     deallocate(fpackB32)
  endif
  if(allocated(gpack))   then
#if defined(OMPGPU)
!$omp target exit data map(release:gpack)
#elif defined(_OPENACC)
!$acc exit data delete(gpack)
#endif
     deallocate(gpack)
  endif
  if(allocated(gpack32))   then
#if defined(OMPGPU)
!$omp target exit data map(release:gpack32)
#elif defined(_OPENACC)
!$acc exit data delete(gpack32)
#endif
     deallocate(gpack32)
  endif
  call cgyro_coll_data_cleanup(gpu_bigmem_flag)

#ifndef CGYRO_GPU_FFT
  if(allocated(fx32))                deallocate(fx32)
  if(allocated(gx32))                deallocate(gx32)
  if(allocated(fy32))                deallocate(fy32)
  if(allocated(gy32))                deallocate(gy32)
  if(allocated(vxmany32))            deallocate(vxmany32)
  if(allocated(vymany32))            deallocate(vymany32)
  if(allocated(uxmany32))            deallocate(uxmany32)
  if(allocated(uymany32))            deallocate(uymany32)
  if(allocated(uv32))                deallocate(uv32)

  if(allocated(fx))                deallocate(fx)
  if(allocated(gx))                deallocate(gx)
  if(allocated(fy))                deallocate(fy)
  if(allocated(gy))                deallocate(gy)
  if(allocated(vxmany))            deallocate(vxmany)
  if(allocated(vymany))            deallocate(vymany)
  if(allocated(uxmany))            deallocate(uxmany)
  if(allocated(uymany))            deallocate(uymany)
  if(allocated(uv))                deallocate(uv)
#endif  

#ifdef CGYRO_GPU_FFT
  if(allocated(fxmany))    then
#if defined(OMPGPU)
!$omp target exit data map(release:fxmany)
#elif defined(_OPENACC)
!$acc exit data delete(fxmany)
#endif
     deallocate(fxmany)
  endif
  if(allocated(gxmany))    then
#if defined(OMPGPU)
!$omp target exit data map(release:gxmany)
#elif defined(_OPENACC)
!$acc exit data delete(gxmany)
#endif
     deallocate(gxmany)
  endif
  if(allocated(fymany))    then
#if defined(OMPGPU)
!$omp target exit data map(release:fymany)
#elif defined(_OPENACC)
!$acc exit data delete(fymany)
#endif
     deallocate(fymany)
  endif
  if(allocated(gymany))    then
#if defined(OMPGPU)
!$omp target exit data map(release:gymany)
#elif defined(_OPENACC)
!$acc exit data delete(gymany)
#endif
     deallocate(gymany)
  endif
  if(allocated(uxmany))    then
#if defined(OMPGPU)
!$omp target exit data map(release:uxmany)
#elif defined(_OPENACC)
!$acc exit data delete(uxmany)
#endif
     deallocate(uxmany)
  endif
  if(allocated(uymany))    then
#if defined(OMPGPU)
!$omp target exit data map(release:uymany)
#elif defined(_OPENACC)
!$acc exit data delete(uymany)
#endif
     deallocate(uymany)
  endif
  if(allocated(vxmany))     then
#if defined(OMPGPU)
!$omp target exit data map(release:vxmany)
#elif defined(_OPENACC)
!$acc exit data delete(vxmany)
#endif
     deallocate(vxmany)
  endif
  if(allocated(vymany))     then
#if defined(OMPGPU)
!$omp target exit data map(release:vymany)
#elif defined(_OPENACC)
!$acc exit data delete(vymany)
#endif
     deallocate(vymany)
  endif
  if(allocated(uvmany))     then
#if defined(OMPGPU)
!$omp target exit data map(release:uvmany)
#elif defined(_OPENACC)
!$acc exit data delete(uvmany)
#endif
     deallocate(uvmany)
  endif
  if(allocated(fxmany32))    then
#if defined(OMPGPU)
!$omp target exit data map(release:fxmany32)
#elif defined(_OPENACC)
!$acc exit data delete(fxmany32)
#endif
     deallocate(fxmany32)
  endif
  if(allocated(gxmany32))    then
#if defined(OMPGPU)
!$omp target exit data map(release:gxmany32)
#elif defined(_OPENACC)
!$acc exit data delete(gxmany32)
#endif
     deallocate(gxmany32)
  endif
  if(allocated(fymany32))    then
#if defined(OMPGPU)
!$omp target exit data map(release:fymany32)
#elif defined(_OPENACC)
!$acc exit data delete(fymany32)
#endif
     deallocate(fymany32)
  endif
  if(allocated(gymany32))    then
#if defined(OMPGPU)
!$omp target exit data map(release:gymany32)
#elif defined(_OPENACC)
!$acc exit data delete(gymany32)
#endif
     deallocate(gymany32)
  endif
  if(allocated(uxmany32))    then
#if defined(OMPGPU)
!$omp target exit data map(release:uxmany32)
#elif defined(_OPENACC)
!$acc exit data delete(uxmany32)
#endif
     deallocate(uxmany32)
  endif
  if(allocated(uymany32))    then
#if defined(OMPGPU)
!$omp target exit data map(release:uymany32)
#elif defined(_OPENACC)
!$acc exit data delete(uymany32)
#endif
     deallocate(uymany32)
  endif
  if(allocated(vxmany32))     then
#if defined(OMPGPU)
!$omp target exit data map(release:vxmany32)
#elif defined(_OPENACC)
!$acc exit data delete(vxmany32)
#endif
     deallocate(vxmany32)
  endif
  if(allocated(vymany32))     then
#if defined(OMPGPU)
!$omp target exit data map(release:vymany32)
#elif defined(_OPENACC)
!$acc exit data delete(vymany32)
#endif
     deallocate(vymany32)
  endif
  if(allocated(uvmany32))     then
#if defined(OMPGPU)
!$omp target exit data map(release:uvmany32)
#elif defined(_OPENACC)
!$acc exit data delete(uvmany32)
#endif
     deallocate(uvmany32)
  endif
#endif    

  !!!!!!!!!!!!!!!!!!!!!!!!!!
  ! From cgyro_make_profiles
  !!!!!!!!!!!!!!!!!!!!!!!!!!
  
  if(allocated(px))  then
#if defined(OMPGPU)
!$omp target exit data map(release:px)
#elif defined(_OPENACC)
!$acc exit data delete(px)
#endif
     deallocate(px)
  endif
  
#if defined(OMPGPU)
!$omp target exit data map(release:z)
#elif defined(_OPENACC)
!$acc exit data delete(z)
#endif

#if defined(OMPGPU)
!$omp target exit data map(release:temp)
#elif defined(_OPENACC)
!$acc exit data delete(temp)
#endif
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!
  ! From cgyro_init_arrays
  !!!!!!!!!!!!!!!!!!!!!!!!!!

  if(allocated(sum_den_h))        deallocate(sum_den_h)
  if(allocated(cderiv))           deallocate(cderiv)
  if(allocated(uderiv))           deallocate(uderiv)
  if(allocated(c_wave)) then
#if defined(OMPGPU)
!$omp target exit data map(release:c_wave)
#elif defined(_OPENACC)
!$acc exit data delete(c_wave)
#endif
     deallocate(c_wave)
  endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!
  ! From cgyro_mpi_grid
  !!!!!!!!!!!!!!!!!!!!!!!!!!

  if(allocated(ie_v)) then
#if defined(OMPGPU)
!$omp target exit data map(release:ie_v)
#elif defined(_OPENACC)
!$acc exit data delete(ie_v)
#endif
     deallocate(ie_v)
  endif

  if(allocated(ix_v)) then
#if defined(OMPGPU)
!$omp target exit data map(release:ix_v)
#elif defined(_OPENACC)
!$acc exit data delete(ix_v)
#endif
     deallocate(ix_v)
  endif

  if(allocated(is_v)) then
#if defined(OMPGPU)
!$omp target exit data map(release:is_v)
#elif defined(_OPENACC)
!$acc exit data delete(is_v)
#endif
     deallocate(is_v)
  endif

  if(allocated(iv_v)) then
#if defined(OMPGPU)
!$omp target exit data map(release:iv_v)
#elif defined(_OPENACC)
!$acc exit data delete(iv_v)
#endif
     deallocate(iv_v)
  endif

    if(allocated(ir_c)) then
#if defined(OMPGPU)
!$omp target exit data map(release:ir_c)
#elif defined(_OPENACC)
!$acc exit data delete(ir_c)
#endif
     deallocate(ir_c)
  endif

  if(allocated(it_c)) then
#if defined(OMPGPU)
!$omp target exit data map(release:it_c)
#elif defined(_OPENACC)
!$acc exit data delete(it_c)
#endif
     deallocate(it_c)
  endif

  if(allocated(ic_c)) then
#if defined(OMPGPU)
!$omp target exit data map(release:ic_c)
#elif defined(_OPENACC)
!$acc exit data delete(ic_c)
#endif
     deallocate(ic_c)
  endif

  if(allocated(ica_c))            deallocate(ica_c)
  if(allocated(icb_c))            deallocate(icb_c)

  !!!!!!!!!!!!!!!!!!!!!!!!!!
  ! From cgyro_parallel_lib
  !!!!!!!!!!!!!!!!!!!!!!!!!!

  call parallel_lib_clean
  
end subroutine cgyro_cleanup
 
