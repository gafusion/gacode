subroutine EXPRO_alloc_control(i_proc,flag)

  use EXPRO_globals
  use EXPRO_interface

  implicit none

  integer, intent(in) :: i_proc
  integer, intent(in) :: flag

  integer :: ierr
  integer, parameter :: io=1


  if (flag == 1) then

     if (i_proc == 0 .and. EXPRO_n_exp == 0) then

        !--------------------------------------------------------------
        ! Read input.profiles.gen if it exists and EXPRO_n_exp not set
        ! 
        open(unit=io,&
             file=trim(path)//'input.profiles.gen',&
             status='old',&
             iostat=ierr)

        if (ierr == 0) then
           read(io,*) EXPRO_shot
           read(io,*) EXPRO_n_ion
           read(io,*) EXPRO_n_exp
           close(io)
           open(unit=io,file=trim(path)//trim(runfile),status='replace')
           print('(a)'),'INFO: (EXPRO) input.profiles.gen found'
           close(io)
        else
           close(io)
           open(unit=io,file=trim(path)//trim(runfile),status='replace')
           print('(a)'),'ERROR: (EXPRO) input.profiles.gen does not exist'
           close(io)
           stop
        endif
        !--------------------------------------------------------------

     endif

     if (i_proc == 0) then

        !--------------------------------------------------------------
        ! Read geometry Fourier harmonics
        ! 
        open(unit=io,&
             file=trim(path)//'input.profiles.geo',&
             status='old',&
             iostat=ierr)

        if (ierr /= 0) then
           EXPRO_nfourier = -1
           close(io)
        else
           call EXPRO_skip_header(io)
           read(io,*) EXPRO_nfourier
           close(io)
           open(unit=io,file=trim(path)//trim(runfile),status='replace')
           write(io,'(a)') 'INFO: (EXPRO) input.profiles.geo found'
           close(io)
        endif

        !--------------------------------------------------------------

     endif

     allocate(EXPRO_rho(EXPRO_n_exp))
     EXPRO_rho(:)=0.0
     allocate(EXPRO_rmin(EXPRO_n_exp))
     EXPRO_rmin(:)=0.0
     allocate(EXPRO_rmaj(EXPRO_n_exp))
     EXPRO_rmaj(:)=0.0
     allocate(EXPRO_q(EXPRO_n_exp))
     EXPRO_q(:)=0.0
     allocate(EXPRO_kappa(EXPRO_n_exp))
     EXPRO_kappa(:)=0.0

     allocate(EXPRO_delta(EXPRO_n_exp))
     EXPRO_delta(:)=0.0
     allocate(EXPRO_te(EXPRO_n_exp))
     EXPRO_te(:)=0.0
     allocate(EXPRO_ne(EXPRO_n_exp))
     EXPRO_ne(:)=0.0
     allocate(EXPRO_z_eff(EXPRO_n_exp))
     EXPRO_z_eff(:)=0.0
     allocate(EXPRO_w0(EXPRO_n_exp))
     EXPRO_w0(:)=0.0

     allocate(EXPRO_flow_mom(EXPRO_n_exp))
     EXPRO_flow_mom(:)=0.0
     allocate(EXPRO_pow_e(EXPRO_n_exp))
     EXPRO_pow_e(:)=0.0
     allocate(EXPRO_pow_i(EXPRO_n_exp))
     EXPRO_pow_i(:)=0.0
     allocate(EXPRO_pow_ei(EXPRO_n_exp))
     EXPRO_pow_ei(:)=0.0
     allocate(EXPRO_zeta(EXPRO_n_exp))
     EXPRO_zeta(:)=0.0

     allocate(EXPRO_sbeame(EXPRO_n_exp))
     EXPRO_sbeame(:)=0.0
     allocate(EXPRO_sbcx(EXPRO_n_exp))
     EXPRO_sbcx(:)=0.0
     allocate(EXPRO_sscxl(EXPRO_n_exp))
     EXPRO_sscxl(:)=0.0

     allocate(EXPRO_flow_beam(EXPRO_n_exp))
     EXPRO_flow_beam(:)=0.0
     allocate(EXPRO_flow_wall(EXPRO_n_exp))
     EXPRO_flow_wall(:)=0.0
     allocate(EXPRO_zmag(EXPRO_n_exp))
     EXPRO_zmag(:)=0.0
     allocate(EXPRO_ptot(EXPRO_n_exp))
     EXPRO_ptot(:)=0.0
     allocate(EXPRO_dlnptotdr(EXPRO_n_exp))
     EXPRO_dlnptotdr(:)=0.0

     allocate(EXPRO_ni(EXPRO_n_ion_max,EXPRO_n_exp))
     EXPRO_ni(:,:)=0.0
     allocate(EXPRO_ti(EXPRO_n_ion_max,EXPRO_n_exp))
     EXPRO_ti(:,:)=0.0
     allocate(EXPRO_vtor(EXPRO_n_ion_max,EXPRO_n_exp))
     EXPRO_vtor(:,:)=0.0
     allocate(EXPRO_vpol(EXPRO_n_ion_max,EXPRO_n_exp))
     EXPRO_vpol(:,:)=0.0

     allocate(EXPRO_pow_e_fus(EXPRO_n_exp))
     EXPRO_pow_e_fus(:)=0.0
     allocate(EXPRO_pow_i_fus(EXPRO_n_exp))
     EXPRO_pow_i_fus(:)=0.0
     allocate(EXPRO_pow_e_sync(EXPRO_n_exp))
     EXPRO_pow_e_sync(:)=0.0
     allocate(EXPRO_pow_e_brem(EXPRO_n_exp))
     EXPRO_pow_e_brem(:)=0.0
     allocate(EXPRO_pow_e_line(EXPRO_n_exp))
     EXPRO_pow_e_line(:)=0.0

     allocate(EXPRO_pow_e_aux(EXPRO_n_exp))
     EXPRO_pow_e_aux(:)=0.0
     allocate(EXPRO_pow_i_aux(EXPRO_n_exp))
     EXPRO_pow_i_aux(:)=0.0

     allocate(EXPRO_bunit(EXPRO_n_exp))
     EXPRO_bunit(:)=0.0
     allocate(EXPRO_s(EXPRO_n_exp))
     EXPRO_s(:)=0.0
     allocate(EXPRO_drmaj(EXPRO_n_exp))
     EXPRO_drmaj(:)=0.0
     allocate(EXPRO_dzmag(EXPRO_n_exp))
     EXPRO_dzmag(:)=0.0
     allocate(EXPRO_sdelta(EXPRO_n_exp))
     EXPRO_sdelta(:)=0.0
     allocate(EXPRO_skappa(EXPRO_n_exp))
     EXPRO_skappa(:)=0.0
     allocate(EXPRO_szeta(EXPRO_n_exp))
     EXPRO_szeta(:)=0.0
     allocate(EXPRO_dlnnedr(EXPRO_n_exp))
     EXPRO_dlnnedr(:)=0.0
     allocate(EXPRO_dlntedr(EXPRO_n_exp))
     EXPRO_dlntedr(:)=0.0
     allocate(EXPRO_sdlnnedr(EXPRO_n_exp))
     EXPRO_sdlnnedr(:)=0.0
     allocate(EXPRO_sdlntedr(EXPRO_n_exp))
     EXPRO_sdlntedr(:)=0.0

     allocate(EXPRO_dlnnidr(EXPRO_n_ion_max,EXPRO_n_exp))
     EXPRO_dlnnidr(:,:)=0.0
     allocate(EXPRO_dlntidr(EXPRO_n_ion_max,EXPRO_n_exp))
     EXPRO_dlntidr(:,:)=0.0
     allocate(EXPRO_sdlnnidr(EXPRO_n_ion_max,EXPRO_n_exp))
     EXPRO_sdlnnidr(:,:)=0.0
     allocate(EXPRO_sdlntidr(EXPRO_n_ion_max,EXPRO_n_exp))
     EXPRO_sdlntidr(:,:)=0.0

     allocate(EXPRO_vol(EXPRO_n_exp))
     EXPRO_vol(:)=0.0
     allocate(EXPRO_volp(EXPRO_n_exp))
     EXPRO_volp(:)=0.0

     allocate(EXPRO_cs(EXPRO_n_exp))
     EXPRO_cs(:)=0.0
     allocate(EXPRO_rhos(EXPRO_n_exp))
     EXPRO_rhos(:)=0.0

     allocate(EXPRO_nuee(EXPRO_n_exp))
     EXPRO_nuee(:)=0.0
     
     allocate(EXPRO_w0p(EXPRO_n_exp))
     EXPRO_w0p(:)=0.0
     allocate(EXPRO_gamma_e(EXPRO_n_exp))
     EXPRO_gamma_e(:)=0.0
     allocate(EXPRO_gamma_p(EXPRO_n_exp))
     EXPRO_gamma_p(:)=0.0
     allocate(EXPRO_mach(EXPRO_n_exp))
     EXPRO_mach(:)=0.0
     allocate(EXPRO_thetascale(EXPRO_n_exp))
     EXPRO_thetascale(:)=1.0

     if (EXPRO_nfourier > 0) then  
        allocate(EXPRO_geo(4,0:EXPRO_nfourier,EXPRO_n_exp))
        EXPRO_geo(:,:,:)=0.0
        allocate(EXPRO_dgeo(4,0:EXPRO_nfourier,EXPRO_n_exp))
        EXPRO_dgeo(:,:,:)=0.0
     endif

     allocate(EXPRO_ni_new(EXPRO_n_exp))
     EXPRO_ni_new(:)=0.0
     allocate(EXPRO_dlnnidr_new(EXPRO_n_exp))
     EXPRO_dlnnidr_new(:)=0.0
     allocate(EXPRO_sdlnnidr_new(EXPRO_n_exp))
     EXPRO_sdlnnidr_new(:)=0.0
     allocate(EXPRO_grad_r0(EXPRO_n_exp))
     EXPRO_grad_r0(:)=0.0
     allocate(EXPRO_ave_grad_r(EXPRO_n_exp))
     EXPRO_ave_grad_r(:)=0.0
     allocate(EXPRO_drdrho(EXPRO_n_exp))
     EXPRO_drdrho(:)=0.0
     allocate(EXPRO_bp0(EXPRO_n_exp))
     EXPRO_bp0(:)=0.0
     allocate(EXPRO_bt0(EXPRO_n_exp))
     EXPRO_bt0(:)=0.0
     allocate(EXPRO_polflux(EXPRO_n_exp))
     EXPRO_polflux(:)=0.0
     allocate(EXPRO_ip(EXPRO_n_exp))
     EXPRO_ip(:)=0.0

  else

     deallocate(EXPRO_rho)
     deallocate(EXPRO_rmin)
     deallocate(EXPRO_rmaj)
     deallocate(EXPRO_q)
     deallocate(EXPRO_kappa)

     deallocate(EXPRO_delta)
     deallocate(EXPRO_te)
     deallocate(EXPRO_ne)
     deallocate(EXPRO_z_eff)
     deallocate(EXPRO_w0)

     deallocate(EXPRO_flow_mom)
     deallocate(EXPRO_pow_e)
     deallocate(EXPRO_pow_i)
     deallocate(EXPRO_pow_ei)
     deallocate(EXPRO_zeta)

     deallocate(EXPRO_sbeame)
     deallocate(EXPRO_sbcx)
     deallocate(EXPRO_sscxl)

     deallocate(EXPRO_flow_beam)
     deallocate(EXPRO_flow_wall)
     deallocate(EXPRO_zmag)
     deallocate(EXPRO_ptot)
     deallocate(EXPRO_polflux)

     deallocate(EXPRO_ni)
     deallocate(EXPRO_ti)
     deallocate(EXPRO_vtor)
     deallocate(EXPRO_vpol)

     deallocate(EXPRO_pow_e_fus)
     deallocate(EXPRO_pow_i_fus)
     deallocate(EXPRO_pow_e_sync)
     deallocate(EXPRO_pow_e_brem)
     deallocate(EXPRO_pow_e_line)

     deallocate(EXPRO_pow_e_aux)
     deallocate(EXPRO_pow_i_aux)

     deallocate(EXPRO_bunit)
     deallocate(EXPRO_s)
     deallocate(EXPRO_drmaj)
     deallocate(EXPRO_dzmag)
     deallocate(EXPRO_sdelta)
     deallocate(EXPRO_skappa)
     deallocate(EXPRO_szeta)
     deallocate(EXPRO_dlnnedr)
     deallocate(EXPRO_dlntedr)
     deallocate(EXPRO_sdlnnedr)
     deallocate(EXPRO_sdlntedr)

     deallocate(EXPRO_dlnnidr)
     deallocate(EXPRO_dlntidr)
     deallocate(EXPRO_sdlnnidr)
     deallocate(EXPRO_sdlntidr)

     deallocate(EXPRO_vol)
     deallocate(EXPRO_volp)

     deallocate(EXPRO_cs)
     deallocate(EXPRO_rhos)

     deallocate(EXPRO_nuee)

     deallocate(EXPRO_w0p)
     deallocate(EXPRO_gamma_e)
     deallocate(EXPRO_gamma_p)
     deallocate(EXPRO_mach)
     deallocate(EXPRO_thetascale)

     if (allocated(EXPRO_geo)) deallocate(EXPRO_geo)
     if (allocated(EXPRO_dgeo)) deallocate(EXPRO_dgeo)

     deallocate(EXPRO_ni_new)
     deallocate(EXPRO_dlnnidr_new)
     deallocate(EXPRO_sdlnnidr_new)
     deallocate(EXPRO_grad_r0)
     deallocate(EXPRO_ave_grad_r)
     deallocate(EXPRO_drdrho)
     deallocate(EXPRO_bp0)
     deallocate(EXPRO_bt0)
     deallocate(EXPRO_dlnptotdr)
     deallocate(EXPRO_ip)

  endif

end subroutine EXPRO_alloc_control
