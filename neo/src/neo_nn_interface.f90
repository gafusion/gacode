

module neo_nn_interface

  implicit none 

  real :: nn_rmin_in
  real :: nn_q_in
  real :: nn_nuee_in
  real :: nn_ni1_ne_in
  real :: nn_ti1_te_in

  real :: nn_charge_norm_fac
  real :: nn_dens_norm_f    
  real :: nn_vnorm
  real :: nn_anorm
  real :: nn_jbs_norm
  real :: nn_I_over_phi_prime
  real :: nn_rho_star
  real :: nn_1_over_Lte
  real :: nn_1_over_Lne
  real :: nn_ni1_dens
  real :: nn_ni2_dens
  real :: nn_1_over_LtD
  real :: nn_1_over_LnD
  real :: nn_1_over_LtC
  real :: nn_1_over_LnC
  real :: nn_NEO_jbs_in_A_m2
  real :: nn_SAU_jbs_in_A_m2
  real :: nn_rho_in
  real :: nn_delta_in
  real :: nn_kappa_in


  real :: OUT_CNEO_CTD
  real :: OUT_CNEO_CTe
  real :: OUT_CNEO_CTC
  real :: OUT_CNEO_CnD
  real :: OUT_CNEO_Cne
  real :: OUT_CNEO_CnC
  real :: OUT_CSAU_CTD
  real :: OUT_CSAU_CTe
  real :: OUT_CSAU_CTC
  real :: OUT_CSAU_CnD
  real :: OUT_CSAU_Cne
  real :: OUT_CSAU_CnC

 
  
end module neo_nn_interface
