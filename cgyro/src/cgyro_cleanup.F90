! This is called after cgyro_kernel to perform clean-ups from CPU and GPU (deallocate)

subroutine cgyro_cleanup
  use cgyro_globals
  use parallel_lib

  !!!!!!!!!!!!!!!!!!!!!!!!!!
  ! From cgyro_init_manager
  !!!!!!!!!!!!!!!!!!!!!!!!!!
  
  if(allocated(energy))        deallocate(energy)
  if(allocated(vel))        then
!$acc exit data delete(vel)
     deallocate(vel)
  endif
  if(allocated(w_e))           deallocate(w_e)
  if(allocated(e_deriv1_mat))  deallocate(e_deriv1_mat)
  if(allocated(e_deriv1_rot_mat))  deallocate(e_deriv1_rot_mat)
  if(allocated(xi))            then
!$acc exit data delete(xi)
     deallocate(xi)
  endif
  if(allocated(w_xi))          deallocate(w_xi)
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
!$acc exit data delete(omega_stream)
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
  if(allocated(freq_err))            deallocate(freq_err)
  if(allocated(fcoef))  then
!$acc exit data delete(fcoef)     
     deallocate(fcoef)
  endif
  if(allocated(gcoef))  then
!$acc exit data delete(gcoef)     
     deallocate(gcoef)
  endif
  if(allocated(field))  then
!$acc exit data delete(field)     
     deallocate(field)
  endif
  if(allocated(field_loc))  then
!$acc exit data delete(field_loc)     
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
  if(allocated(icd_c))  then
!$acc exit data delete(icd_c)     
     deallocate(icd_c)
  endif
  if(allocated(dtheta)) then
!$acc exit data delete(dtheta)     
     deallocate(dtheta)
  endif
  if(allocated(dtheta_up))  then
!$acc exit data delete(dtheta_up)     
     deallocate(dtheta_up)
  endif
  if(allocated(source)) then
 !$acc exit data delete(source)      
     deallocate(source)
  endif
  if(allocated(h0_old)) then
!$acc exit data delete(h0_old)      
     deallocate(h0_old)
  endif
  if(allocated(rhs)) then
!$acc exit data delete(rhs)       
     deallocate(rhs)
  endif
  if(allocated(h_x)) then
!$acc exit data delete(h_x)       
     deallocate(h_x)
  endif
  if(allocated(g_x)) then
!$acc exit data delete(g_x)       
     deallocate(g_x)
  endif
  if(allocated(h0_x)) then
!$acc exit data delete(h0_x)        
     deallocate(h0_x)
  endif
  if(allocated(cap_h_c)) then
!$acc exit data delete(cap_h_c)       
     deallocate(cap_h_c)
  endif
  if(allocated(cap_h_ct))  then
!$acc exit data delete(cap_h_ct)       
     deallocate(cap_h_ct)
  endif
  if(allocated(cap_h_c_dot)) then
!$acc exit data delete(cap_h_c_dot)
     deallocate(cap_h_c_dot)
  endif
  if(allocated(cap_h_c_old)) then
!$acc exit data delete(cap_h_c_old)
     deallocate(cap_h_c_old)
  endif
  if(allocated(cap_h_c_old2)) then
!$acc exit data delete(cap_h_c_old2)
     deallocate(cap_h_c_old2)
  endif
  if(allocated(omega_cap_h)) then
!$acc exit data delete(omega_cap_h)        
     deallocate(omega_cap_h)
  endif
  if(allocated(omega_h)) then
!$acc exit data delete(omega_h)        
     deallocate(omega_h)
  endif
  if(allocated(omega_s)) then
!$acc exit data delete(omega_s)        
     deallocate(omega_s)
  endif
  if(allocated(omega_ss)) then
!$acc exit data delete(omega_ss)        
     deallocate(omega_ss)
  endif
  if(allocated(jvec_c))  then
!$acc exit data delete(jvec_c)     
     deallocate(jvec_c)
  endif
  if(allocated(jvec_c_nl))  then
!$acc exit data delete(jvec_c_nl)
     deallocate(jvec_c_nl)
  endif
  if(allocated(jvec_v))              deallocate(jvec_v)
  if(allocated(dvjvec_c)) then
!$acc exit data delete(dvjvec_c)     
     deallocate(dvjvec_c)
  endif
  if(allocated(dvjvec_v)) then
!$acc exit data delete(dvjvec_v)     
     deallocate(dvjvec_v)
  endif
  if(allocated(jxvec_c))  then
     deallocate(jxvec_c)
  endif
  if(allocated(upfac1))   then
!$acc exit data delete(upfac1)
     deallocate(upfac1)
  endif
  if(allocated(upfac2))   then
!$acc exit data delete(upfac2)
     deallocate(upfac2)
  endif
  if(allocated(cap_h_v))  then
!$acc exit data delete(cap_h_v)     
     deallocate(cap_h_v)
  endif
  if(allocated(upwind_res_loc))   then
!$acc exit data delete(upwind_res_loc)
     deallocate(upwind_res_loc)
  endif
  if(allocated(upwind_res))   then
!$acc exit data delete(upwind_res)
     deallocate(upwind_res)
  endif
  if(allocated(upwind32_res_loc))   then
!$acc exit data delete(upwind32_res_loc)
     deallocate(upwind32_res_loc)
  endif
  if(allocated(upwind32_res))   then
!$acc exit data delete(upwind32_res)
     deallocate(upwind32_res)
  endif
  if(allocated(f_nl))   then
!$acc exit data delete(f_nl)
     deallocate(f_nl)
  endif
  if(allocated(g_nl))   then
!$acc exit data delete(g_nl)
     deallocate(g_nl)
  endif
  if(allocated(fpack))   then
!$acc exit data delete(fpack)
     deallocate(fpack)
  endif
  if(allocated(gpack))   then
!$acc exit data delete(gpack)
     deallocate(gpack)
  endif
  if (allocated(cmat)) then
!$acc exit data delete(cmat) if (gpu_bigmem_flag == 1)
     deallocate(cmat)
  endif
  if (allocated(cmat_fp32)) then
!$acc exit data delete(cmat_fp32) if (gpu_bigmem_flag == 1)
     deallocate(cmat_fp32)
  endif
  if (allocated(cmat_stripes)) then
!$acc exit data delete(cmat_stripes) if (gpu_bigmem_flag == 1)
     deallocate(cmat_stripes)
  endif
    if (allocated(cmat_simple)) then
!$acc exit data delete(cmat_simple)     
     deallocate(cmat_simple)
  endif

#ifndef _OPENACC
  if(allocated(fx))                deallocate(fx)
  if(allocated(gx))                deallocate(gx)
  if(allocated(fy))                deallocate(fy)
  if(allocated(gy))                deallocate(gy)
  if(allocated(uxmany))            deallocate(uxmany)
  if(allocated(uymany))            deallocate(uymany)
  if(allocated(vx))                deallocate(vx)
  if(allocated(vy))                deallocate(vy)
  if(allocated(uv))                deallocate(uv)
#endif  

#ifdef _OPENACC
  if(allocated(fxmany))    then
!$acc exit data delete(fxmany)     
     deallocate(fxmany)
  endif
  if(allocated(gxmany))    then
!$acc exit data delete(gxmany)     
     deallocate(gxmany)
  endif
  if(allocated(fymany))    then
!$acc exit data delete(fymany)     
     deallocate(fymany)
  endif
  if(allocated(gymany))    then
!$acc exit data delete(gymany)     
     deallocate(gymany)
  endif
  if(allocated(uxmany))    then
 !$acc exit data delete(uxmany)    
     deallocate(uxmany)
  endif
  if(allocated(uymany))    then
!$acc exit data delete(uymany)     
     deallocate(uymany)
  endif
  if(allocated(vxmany))     then
!$acc exit data delete(vxmany)     
     deallocate(vxmany)
  endif
  if(allocated(vymany))     then
!$acc exit data delete(vymany)     
     deallocate(vymany)
  endif
  if(allocated(uvmany))     then
!$acc exit data delete(uvmany)     
     deallocate(uvmany)
  endif
#endif    

  !!!!!!!!!!!!!!!!!!!!!!!!!!
  ! From cgyro_make_profiles
  !!!!!!!!!!!!!!!!!!!!!!!!!!
  
  if(allocated(px))  then
!$acc exit data delete(px)      
     deallocate(px)
  endif

  if(allocated(geo_yin))          deallocate(geo_yin)
  
!$acc exit data delete(z)

!$acc exit data delete(temp) 
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!
  ! From cgyro_init_arrays
  !!!!!!!!!!!!!!!!!!!!!!!!!!

  if(allocated(vfac))             deallocate(vfac)
  if(allocated(sum_den_h))        deallocate(sum_den_h)
  if(allocated(cderiv))           deallocate(cderiv)
  if(allocated(uderiv))           deallocate(uderiv)
  if(allocated(c_wave)) then
!$acc exit data delete(c_wave)     
     deallocate(c_wave)
  endif
  if(allocated(hzf))              deallocate(hzf)
  if(allocated(xzf))              deallocate(xzf)

  !!!!!!!!!!!!!!!!!!!!!!!!!!
  ! From cgyro_mpi_grid
  !!!!!!!!!!!!!!!!!!!!!!!!!!

  if(allocated(ie_v)) then
!$acc exit data delete(ie_v)     
     deallocate(ie_v)
  endif

  if(allocated(ix_v)) then
!$acc exit data delete(ix_v)     
     deallocate(ix_v)
  endif

  if(allocated(is_v)) then
!$acc exit data delete(is_v)     
     deallocate(is_v)
  endif

  if(allocated(iv_v)) then
!$acc exit data delete(ie_v)     
     deallocate(iv_v)
  endif

    if(allocated(ir_c)) then
!$acc exit data delete(ir_c)     
     deallocate(ir_c)
  endif

  if(allocated(it_c)) then
!$acc exit data delete(it_c)     
     deallocate(it_c)
  endif

  if(allocated(ic_c)) then
!$acc exit data delete(ic_c)     
     deallocate(ic_c)
  endif

  if(allocated(ica_c))            deallocate(ica_c)
  if(allocated(icb_c))            deallocate(icb_c)

  !!!!!!!!!!!!!!!!!!!!!!!!!!
  ! From cgyro_parallel_lib
  !!!!!!!!!!!!!!!!!!!!!!!!!!

  call parallel_lib_clean
  
end subroutine cgyro_cleanup
 
