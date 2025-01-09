! This is called after cgyro_kernel to perform clean-ups from CPU and GPU (deallocate)

subroutine cgyro_cleanup
  use cgyro_globals
  use parallel_lib

#if defined(_OPENACC) || defined(OMPGPU)
#define CGYRO_GPU_FFT
#endif

#if defined(OMPGPU)

#define ccl_del_device(x) \
!$omp target exit data map(release:x)
#define ccl_del_bigdevice(x) \
!$omp target exit data map(release:x) if (gpu_bigmem_flag == 1)

#elif defined(_OPENACC)

#define ccl_del_device(x) \
!$acc exit data delete(x)
#define ccl_del_bigdevice(x) \
!$acc exit data delete(x) if (gpu_bigmem_flag == 1)

#else

  ! nothing to do
#define ccl_del_device(x)
#define ccl_del_bigdevice(x)

#endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!
  ! From cgyro_init_manager
  !!!!!!!!!!!!!!!!!!!!!!!!!!
  
  if(allocated(energy))        deallocate(energy)
  if(allocated(vel))        then
     ccl_del_device(vel)
     deallocate(vel)
  endif
  if(allocated(vel2))        then
     ccl_del_device(vel2)
     deallocate(vel2)
  endif
  if(allocated(w_e))           deallocate(w_e)
  if(allocated(e_deriv1_mat))  deallocate(e_deriv1_mat)
  if(allocated(e_deriv1_rot_mat))  deallocate(e_deriv1_rot_mat)
  if(allocated(xi))            then
     ccl_del_device(xi)
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
     ccl_del_device(omega_stream)
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
  if(allocated(fcoef))  then
     ccl_del_device(fcoef)     
     deallocate(fcoef)
  endif
  if(allocated(gcoef))  then
     ccl_del_device(gcoef)     
     deallocate(gcoef)
  endif
  if(allocated(field_v))  then
     ccl_del_device(field_v)     
     deallocate(field_v)
  endif
  if(allocated(field_loc_v))  then
     ccl_del_device(field_loc_v)     
     deallocate(field_loc_v)
  endif
  if(allocated(field))  then
     ccl_del_device(field)     
     deallocate(field)
  endif
  if(allocated(field_loc))  then
     ccl_del_device(field_loc)     
     deallocate(field_loc)
  endif
  if(allocated(field_dot))           deallocate(field_dot)
  if(allocated(field_old))           deallocate(field_old)
  if(allocated(field_old2))          deallocate(field_old2)
  if(allocated(field_old3))          deallocate(field_old3)
  if(allocated(moment))              deallocate(moment)
  if(allocated(moment_loc))          deallocate(moment_loc)
  if(allocated(cflux))               deallocate(cflux)
  if(allocated(cflux_loc))           deallocate(cflux_loc)
  if(allocated(gflux))               deallocate(gflux)
  if(allocated(gflux_loc))           deallocate(gflux_loc)
  if(allocated(cflux_tave))          deallocate(cflux_tave)
  if(allocated(gflux_tave))          deallocate(gflux_tave)
  if(allocated(recv_status))         deallocate(recv_status)
  if(allocated(source)) then
      ccl_del_device(source)      
     deallocate(source)
  endif
  if(allocated(triad)) then
      ccl_del_device(triad)      
     deallocate(triad)
  endif
  if(allocated(triad_loc)) then
      ccl_del_device(triad_loc)      
     deallocate(triad_loc)
  endif
  if(allocated(triad_loc_old)) then
      ccl_del_device(triad_loc_old)      
     deallocate(triad_loc_old)
  endif
  if(allocated(thfac_itor)) then
      ccl_del_device(thfac_itor)
     deallocate(thfac_itor)
  endif
  if(allocated(h0_old)) then
     ccl_del_device(h0_old)      
     deallocate(h0_old)
  endif
  if(allocated(rhs)) then
     ccl_del_device(rhs)       
     deallocate(rhs)
  endif
  if(allocated(h_x)) then
     ccl_del_device(h_x)       
     deallocate(h_x)
  endif
  if(allocated(g_x)) then
     ccl_del_device(g_x)       
     deallocate(g_x)
  endif
  if(allocated(h0_x)) then
     ccl_del_device(h0_x)        
     deallocate(h0_x)
  endif
  if(allocated(cap_h_c)) then
     ccl_del_device(cap_h_c)       
     deallocate(cap_h_c)
  endif
  if(allocated(cap_h_ct))  then
     ccl_del_device(cap_h_ct)       
     deallocate(cap_h_ct)
  endif
  if(allocated(cap_h_c_dot)) then
     ccl_del_device(cap_h_c_dot)
     deallocate(cap_h_c_dot)
  endif
  if(allocated(cap_h_c_old)) then
     ccl_del_device(cap_h_c_old)
     deallocate(cap_h_c_old)
  endif
  if(allocated(cap_h_c_old2)) then
     ccl_del_device(cap_h_c_old2)
     deallocate(cap_h_c_old2)
  endif
  if(allocated(omega_cap_h)) then
     ccl_del_device(omega_cap_h)        
     deallocate(omega_cap_h)
  endif
  if(allocated(omega_h)) then
     ccl_del_device(omega_h)        
     deallocate(omega_h)
  endif
  if(allocated(omega_s)) then
     ccl_del_device(omega_s)        
     deallocate(omega_s)
  endif
  if(allocated(omega_ss)) then
     ccl_del_device(omega_ss)        
     deallocate(omega_ss)
  endif
  if(allocated(omega_sbeta)) then
     ccl_del_device(omega_sbeta)        
     deallocate(omega_sbeta)
  endif
  if(allocated(diss_r)) then
     ccl_del_device(diss_r)        
     deallocate(diss_r)
  endif
  if(allocated(jvec_c))  then
     ccl_del_device(jvec_c)     
     deallocate(jvec_c)
  endif
  if(allocated(jvec_c_nl))  then
     ccl_del_device(jvec_c_nl)
     deallocate(jvec_c_nl)
  endif
  if(allocated(jvec_v))              deallocate(jvec_v)
  if(allocated(dvjvec_c)) then
     ccl_del_device(dvjvec_c)     
     deallocate(dvjvec_c)
  endif
  if(allocated(dvjvec_v)) then
     ccl_del_device(dvjvec_v)     
     deallocate(dvjvec_v)
  endif
  if(allocated(jxvec_c))  then
     deallocate(jxvec_c)
  endif
  if(allocated(upfac1))   then
     ccl_del_device(upfac1)
     deallocate(upfac1)
  endif
  if(allocated(upfac2))   then
     ccl_del_device(upfac2)
     deallocate(upfac2)
  endif
  if(allocated(cap_h_v))  then
     ccl_del_device(cap_h_v)     
     deallocate(cap_h_v)
  endif
  if(allocated(upwind_res_loc))   then
     ccl_del_device(upwind_res_loc)
     deallocate(upwind_res_loc)
  endif
  if(allocated(upwind_res))   then
     ccl_del_device(upwind_res)
     deallocate(upwind_res)
  endif
  if(allocated(upwind32_res_loc))   then
     ccl_del_device(upwind32_res_loc)
     deallocate(upwind32_res_loc)
  endif
  if(allocated(upwind32_res))   then
     ccl_del_device(upwind32_res)
     deallocate(upwind32_res)
  endif
  if(allocated(fA_nl))   then
     ccl_del_device(fA_nl)
     deallocate(fA_nl)
  endif
  if(allocated(fB_nl))   then
     ccl_del_device(fB_nl)
     deallocate(fB_nl)
  endif
  if(allocated(g_nl))   then
     ccl_del_device(g_nl)
     deallocate(g_nl)
  endif
  if(allocated(eA_nl))   then
     ccl_del_device(eA_nl)
     deallocate(eA_nl)
  endif
  if(allocated(eB_nl))   then
     ccl_del_device(eB_nl)
     deallocate(eB_nl)
  endif
  if(allocated(fpackA))   then
     ccl_del_device(fpackA)
     deallocate(fpackA)
  endif
  if(allocated(fpackB))   then
     ccl_del_device(fpackB)
     deallocate(fpackB)
  endif
  if(allocated(gpack))   then
     ccl_del_device(gpack)
     deallocate(gpack)
  endif
  if(allocated(epackA))   then
     ccl_del_device(epackA)
     deallocate(epackA)
  endif
  if(allocated(epackB))   then
     ccl_del_device(epackB)
     deallocate(epackB)
  endif
  if (allocated(cmat)) then
     ccl_del_bigdevice(cmat)
     deallocate(cmat)
  endif
  if (allocated(cmat_fp32)) then
     ccl_del_bigdevice(cmat_fp32)
     deallocate(cmat_fp32)
  endif
  if (allocated(cmat_stripes)) then
     ccl_del_bigdevice(cmat_stripes)
     deallocate(cmat_stripes)
  endif
    if (allocated(cmat_simple)) then
     ccl_del_device(cmat_simple)     
     deallocate(cmat_simple)
  endif

#ifndef CGYRO_GPU_FFT
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
     ccl_del_device(fxmany)     
     deallocate(fxmany)
  endif
  if(allocated(gxmany))    then
     ccl_del_device(gxmany)     
     deallocate(gxmany)
  endif
  if(allocated(fymany))    then
     ccl_del_device(fymany)     
     deallocate(fymany)
  endif
  if(allocated(gymany))    then
     ccl_del_device(gymany)     
     deallocate(gymany)
  endif
  if(allocated(uxmany))    then
      ccl_del_device(uxmany)    
     deallocate(uxmany)
  endif
  if(allocated(uymany))    then
     ccl_del_device(uymany)     
     deallocate(uymany)
  endif
  if(allocated(vxmany))     then
     ccl_del_device(vxmany)     
     deallocate(vxmany)
  endif
  if(allocated(vymany))     then
     ccl_del_device(vymany)     
     deallocate(vymany)
  endif
  if(allocated(uvmany))     then
     ccl_del_device(uvmany)     
     deallocate(uvmany)
  endif
#endif    

  !!!!!!!!!!!!!!!!!!!!!!!!!!
  ! From cgyro_make_profiles
  !!!!!!!!!!!!!!!!!!!!!!!!!!
  
  if(allocated(px))  then
     ccl_del_device(px)      
     deallocate(px)
  endif

  if(allocated(geo_yin))          deallocate(geo_yin)
  
     ccl_del_device(z)

     ccl_del_device(temp) 
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!
  ! From cgyro_init_arrays
  !!!!!!!!!!!!!!!!!!!!!!!!!!

  if(allocated(vfac))             deallocate(vfac)
  if(allocated(sum_den_h))        deallocate(sum_den_h)
  if(allocated(cderiv))           deallocate(cderiv)
  if(allocated(uderiv))           deallocate(uderiv)
  if(allocated(c_wave)) then
     ccl_del_device(c_wave)     
     deallocate(c_wave)
  endif
  if(allocated(hzf))              deallocate(hzf)
  if(allocated(xzf))              deallocate(xzf)

  !!!!!!!!!!!!!!!!!!!!!!!!!!
  ! From cgyro_mpi_grid
  !!!!!!!!!!!!!!!!!!!!!!!!!!

  if(allocated(ie_v)) then
     ccl_del_device(ie_v)     
     deallocate(ie_v)
  endif

  if(allocated(ix_v)) then
     ccl_del_device(ix_v)     
     deallocate(ix_v)
  endif

  if(allocated(is_v)) then
     ccl_del_device(is_v)     
     deallocate(is_v)
  endif

  if(allocated(iv_v)) then
     ccl_del_device(iv_v)     
     deallocate(iv_v)
  endif

    if(allocated(ir_c)) then
     ccl_del_device(ir_c)     
     deallocate(ir_c)
  endif

  if(allocated(it_c)) then
     ccl_del_device(it_c)     
     deallocate(it_c)
  endif

  if(allocated(ic_c)) then
     ccl_del_device(ic_c)     
     deallocate(ic_c)
  endif

  if(allocated(ica_c))            deallocate(ica_c)
  if(allocated(icb_c))            deallocate(icb_c)

  !!!!!!!!!!!!!!!!!!!!!!!!!!
  ! From cgyro_parallel_lib
  !!!!!!!!!!!!!!!!!!!!!!!!!!

  call parallel_lib_clean
  
end subroutine cgyro_cleanup
 
