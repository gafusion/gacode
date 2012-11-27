!------------------------------------------------------------
! prgen_map_peqdsk.f90
!
! PURPOSE:
!  Map native peqdsk data onto input.profiles standard.  
!------------------------------------------------------------

subroutine prgen_map_peqdsk

  use prgen_globals

  implicit none


  ! Compute rho, bref and arho:
  call prgen_get_chi(nx,q_gato,kappa,rmin,dpsi,rho,peqdsk_bref,peqdsk_arho)

  !---------------------------------------------------------
  ! Map profile data onto single array:
  !
  allocate(vec(n_indx,peqdsk_nj))
  vec(:,:) = 0.0
  !
  vec(1,:)  = rho(:)
  vec(2,:)  = rmin(:)
  vec(3,:)  = rmaj(:)
  vec(4,:)  = q(:)
  vec(5,:)  = kappa(:)
  vec(6,:)  = delta(:)
  vec(7,:)  = peqdsk_te(:)
  vec(8,:)  = peqdsk_ne(:)*10
  vec(9,:)  = 0.0      ! zeff
  vec(10,:) = -1e3*peqdsk_omegat(:) ! Omega_tor
  vec(11,:) = 0.0      ! flow_mom
  vec(12,:) = 0.0      ! pow_e
  vec(13,:) = 0.0      ! pow_i 
  vec(14,:) = 0.0      ! pow_ei_exp
  vec(15,:) = zeta(:)
  vec(16,:) = 0.0      ! flow_beam
  vec(17,:) = 0.0      ! flow_wall_exp
  vec(18,:) = zmag(:)  
  vec(19,:) = 0.0      
  vec(20,:) = dpsi(:)

  ! ni
  vec(21,:) = peqdsk_ni(:)*10

  ! ti
  vec(26,:) = peqdsk_ti(:)

  ! vphi
  vec(31,:) = 0.0
  vec(32,:) = peqdsk_omegat(:)*1000*(rmaj(:)+rmin(:))
  vec(33,:) = 0.0
  vec(34,:) = 0.0
  vec(35,:) = 0.0

  ! vpol
  vec(36,:) = 0.0
  vec(37,:) = 0.0
  vec(38,:) = 0.0
  vec(39,:) = 0.0
  vec(40,:) = 0.0

  signpsi = abs(dpsi(peqdsk_nj)-dpsi(1))/&
       (dpsi(peqdsk_nj)-dpsi(1))

  ! Construct impurity temperature and density profiles assuming 
  ! Z2 given by pfile_z2:
  if (pfile_z2 > 0.0) then
     vec(27,:) = vec(26,:)
     vec(22,:) = (vec(8,:)-vec(21,:))/pfile_z2
  endif

end subroutine prgen_map_peqdsk
