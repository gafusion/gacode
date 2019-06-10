      SUBROUTINE tglf_startup
!*********************************************
!
!     initalization of tglf 
!
!*********************************************
      USE tglf_dimensions
      USE tglf_global
!
      IMPLICIT NONE
      INTEGER :: i
!
! initialized trace_path that records the flow of tglf internal calls
!
      do i=1,7
        trace_path(i)=0
      enddo
!
! set global constants
!
      pi = atan2( 0.0, -1.0 )
      pi_2 = 2.0*pi
      sqrt_pi = SQRT(pi)
      sqrt_two = SQRT(2.0)
!
!
! set dimensions of allocatable arrays
!
      nx = 2*nxgrid_in-1 
      nbasis = nbasis_max_in
      nbasis_max = nbasis
      ns0 = 1
      if(adiabatic_elec_in) then
        ns0 = 2
      endif
      if(use_default_species)ns_in=2
      ns = ns_in
! SAT_RULE=0 has different fits for nmodes=2 and nmodes=4
      if(sat_rule_in.eq.0)then
        if(nmodes_in.gt.2)nmodes_in=4
      endif
! debug
!      write(*,*)"nx=",nx
!      write(*,*)"nbasis =",nbasis    
!      write(*,*)"ns0 = ",ns0
!      write(*,*)"ns = ",ns
!      write(*,*)"nxgrid_in=",nxgrid_in
!      write(*,*)"nky_in=",nky_in
!      write(*,*)"nmodes_in=",nmodes_in
!      write(*,*)"nbasis_max_in=",nbasis_max_in
!      write(*,*)"nbasis_min_in=",nbasis_min_in
!      write(*,*)"nwidth_in=",nwidth_in
! 
! allocate internal matricies
!
      CALL tglf_allocate
!
! fill species arrays
!
        trace_path(3)=1
      CALL get_species
!      write(*,*)"nxgrid_in=",nxgrid_in
!      write(*,*)"nky_in=",nky_in
!      write(*,*)"nmodes_in=",nmodes_in
!      write(*,*)"nbasis_max_in=",nbasis_max_in
!      write(*,*)"nbasis_min_in=",nbasis_min_in
!      write(*,*)"nwidth_in=",nwidth_in
!
      if(new_start)then
        trace_path(2)=1
        new_width=.TRUE.
      endif
      new_start = .FALSE.
!
!   load the xgrid and gauss-hermite weights
!
      if(gauher_uncalled)CALL gauher
!
! load the hermite basis functions on the xgrid
!  
      if(gauss_hermite_uncalled)CALL gauss_hermite
!
      END SUBROUTINE tglf_startup
!
!
      SUBROUTINE get_species
!*********************************************
!
!*********************************************
      USE tglf_dimensions
      USE tglf_global
      USE tglf_species
!
      IMPLICIT NONE
      INTEGER :: is
      REAL :: xnu
!
! electrons=1, ions =2,...
!
      if(use_default_species)then
       mass_in(1) = 0.0002723
       mass_in(2) = 1.0
       zs_in(1) = -1.0
       zs_in(2) = 1.0
       as_in(1) = 1.0
       as_in(2) = 1.0
       taus_in(1)=1.0
       taus_in(2)=1.0
      endif
!
      nfields_out = 1
      if(use_bper_in)nfields_out = nfields_out + 1
      if(use_bpar_in)nfields_out = nfields_out + 1
!
      pol = 0.0
      U0 = 0.0
      do is=1,ns
        rlns(is) = rlns_in(is)
        rlts(is) = rlts_in(is)
        taus(is) = taus_in(is)
        as(is) = as_in(is)
        vpar_s(is)=0.0
        if(vpar_model_in.eq.0)vpar_s(is)=alpha_mach_in*vpar_in(is)
        vpar_shear_s(is)=alpha_p_in*vpar_shear_in(is)
!        if(nbasis_min_in.eq.1.and.(vpar_shear_s(is).ne.0.0.or.vpar_s(is).ne.0.0))then
!          nbasis_min_in = 2      
!        endif
        zs(is) = zs_in(is)
        mass(is) = mass_in(is)
        vs(is) = SQRT(taus(is)/mass(is))
        pol = pol +  zs(is)*zs(is)*as(is)/taus(is)
        U0 = U0 + as(is)*vpar_s(is)*zs(is)*zs(is)/taus(is)
!        write(*,*)"species",is
!        write(*,*)rlns(is),rlts(is)
!        write(*,*)taus(is),as(is)
!        write(*,*)zs(is),mass(is)
      enddo
!
      vexb_shear_s = vexb_shear_in*sign_It_in
!      if(vpar_shear_model_in.eq.1)then  
!        vexb_shear_s = sign_Bt_in*vexb_shear_in
!      endif
!
!
      xnu = 0.0
! energy exchange
      ei_exch(1,1) = -3.0*xnu*taus(1)*mass(1)/mass(2)
      ei_exch(1,2) =  3.0*xnu*mass(1)/mass(2)
      ei_exch(2,1) = 3.0*xnu*taus(2)*mass(1)/mass(2)
      ei_exch(2,2) = -3.0*xnu*mass(1)/mass(2)
! resistivity
!      xnu = xnue_in*4.0*0.513/(3.0*sqrt_pi)
      resist(1,1) = - xnu
      resist(1,2) = xnu*as(2)*vs(2)/vs(1)
      resist(2,1) = xnu*(vs(1)/vs(2))*(mass(1)/mass(2))*as(1)/as(2) 
      resist(2,2) = - xnu*mass(1)/mass(2)
!
      END SUBROUTINE get_species
!
!
