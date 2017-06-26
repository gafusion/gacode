!
      SUBROUTINE tglf_max
!
      USE tglf_global
      USE tglf_dimensions
      USE tglf_species
      USE tglf_pkg
!
      IMPLICIT NONE
      LOGICAL :: save_iflux
      INTEGER :: nt,i,is,imax
      INTEGER :: save_nbasis
      INTEGER :: branch
      REAL :: width_max=1.65
      REAL :: width_min=0.15
      REAL :: tp,dt
      REAL :: t1,t2,tm,g1,g2,gm
      REAL :: tmax,tmin,gamma_max,gmax
      REAL :: save_width,dtmin
      REAL :: gamma_n(nt0),freq_n(nt0),width_n(nt0)
      REAL :: save_vexb_shear
      REAL :: save_alpha_kx_p
      REAL :: wgp_max,width_p_max 
      REAL :: kyi,ft2
!
      CALL tglf_setup_geometry
!
      do i=1,nmodes_in
        gamma_reference_kx0(i)=0.0
        freq_reference_kx0(i)=0.0
      enddo
      save_iflux = iflux_in
      save_nbasis = nbasis_max_in
      save_width = width_in
      save_vexb_shear = vexb_shear_in
      if(alpha_quench_in.eq.0.0)vexb_shear_in = 0.0
      save_alpha_kx_p = alpha_kx_p_in
      alpha_kx_p_in=0.0
      width_min = width_min_in
      width_max = ABS(width_in)
!
!      write(*,*)"R_unit=",R_unit,"q_unit=",q_unit
!      write(*,*)ns0,ns,ky_in
      do is=ns0,ns
        kyi = ky_in*SQRT(taus(is)*mass(is))/ABS(zs(is))
        wgp_max = ABS((taus(is)/zs(is))*vpar_shear_in(is)/vs(is))*ky_in/(1+kyi**2)
        width_p_max = 3.6*vs(is)/(sqrt_two*R_unit*q_unit*MAX(wgp_max,0.001))
        width_p_max=MAX(width_p_max,0.01)
         if(width_p_max.lt.width_min_in)then
          width_min = width_p_max
        endif
      enddo
!        kyi = ky_in*SQRT(taus_in(2)*mass_in(2))/ABS(zs_in(2))
!        wgp_max = ABS(vpar_shear_in(2)/vs(2))*kyi/(1+kyi**2)
!        width_p_max = 3.6/(sqrt_two*R_unit*q_unit*MAX(wgp_max,0.001))
!        width_p_max=MAX(width_p_max,0.01)
!         if(width_p_max.lt.width_min_in)then
!          width_min = width_p_max
!        endif
!      write(*,*)ky," width_p_max = ", width_p_max,width_min
!
! for ibranch_in > 0 the most unstable positive frequency mode is stored 
! in gamma_out(1) and the most unstable negative frequency mode 
! (ion diamagnetic drift direction) is stored in gamma_out(2)
! for ibranch_in = -1 the unstable modes in rank order are stored in
! gamma_out(i) i=1,nmodes_in
!
      if(ibranch_in.gt.0)branch = ibranch_in
      if(ibranch_in.eq.-1)branch = 1
      iflux_in=.FALSE.
      if(nbasis_min_in.ne.0)then
        nbasis = nbasis_min_in
      endif
!       write(*,*)"nbasis = ",nbasis
      tmin=LOG10(width_min)
      tmax=LOG10(width_max)
      nt=nwidth_in
      dtmin=(tmax-tmin)/REAL(nt-1)
      if(use_bisection_in)nt=5
!
!      open(2,file = 'width.dat',status='unknown')
      dt = (tmax-tmin)/REAL(nt-1)
      tp = tmin
      do i=1,nt
       gamma_n(i)=0.0
       freq_n(i)=0.0
       width_n(i)=0.0
       tp = tmin + REAL(i-1)*dt
       width_in = 10.0**tp
!       write(*,*)"width_in = ",width_in
       new_width = .TRUE.
       call tglf_LS
       if(ibranch_in.eq.0)then
        branch = 1
        if(gamma_out(2).gt.gamma_out(1))branch = 2
       endif
!       write(*,*)i,width_in,gamma_out(branch),freq_out(branch),ft
       width_n(i)=width_in
       gamma_n(i) = gamma_out(branch)
       freq_n(i) = freq_out(branch)
      enddo
!      close(2)
! find the global maximum
      gamma_max=gamma_n(nt)
      imax=nt
      do i=nt-1,1,-1
        if(gamma_n(i).gt.gamma_max)then
          gamma_max=gamma_n(i)
          imax=i
        endif
!        write(*,*)"debug",i,imax,gamma_max,gamma_n(i)
      enddo
      width_in = width_n(imax)
!
      if(use_bisection_in.and.gamma_max.gt.0.0)then
!
! use bounded bisection search to refine width
!
         if(imax.eq.1)then
! maximum is against bottom width
           g1 = gamma_n(1)
           t1 = tmin
           g2 = gamma_n(2)
           t2 = LOG10(width_n(2))
           tp = (t2+t1)/2.0    
           width_in = 10.0**tp
           new_width = .TRUE.
           call tglf_LS
           if(ibranch_in.eq.0)then
            branch = 1
            if(gamma_out(2).gt.gamma_out(1))branch = 2
           endif
           gm = gamma_out(branch)
           tm = tp
         elseif(imax.eq.nt)then
! maximum is against top width
           g1 = gamma_n(nt-1)
           t1 = LOG10(width_n(nt-1))
           g2 = gamma_n(nt)
           t2 = tmax           
           tp = (t2+t1)/2.0    
           width_in = 10.0**tp
           new_width = .TRUE.
           call tglf_LS
           if(ibranch_in.eq.0)then
             branch = 1
            if(gamma_out(2).gt.gamma_out(1))branch = 2
           endif
           gm = gamma_out(branch)
           tm = tp
         else
! maximum is away from boundaries
           g1 = gamma_n(imax-1)
           t1 = LOG10(width_n(imax-1))
           g2 = gamma_n(imax+1)
           t2 = LOG10(width_n(imax+1))
           gm = gamma_n(imax)
           tm = LOG10(width_n(imax))
         endif
! start bisection search
         dt=(t2-t1)/2.0
!         write(*,*)"dtmin=",dtmin,"tmin=",tmin,"tmax=",tmax
         do while(dt.gt.dtmin)
           dt=dt/2.0
           gmax = MAX(gm,MAX(g1,g2))
!           write(*,*)"dt=",dt,gmax
!           write(*,*)g1,gm,g2
!           write(*,*)t1,tm,t2
!
           if(g1.eq.gmax)then  
             if(t1.gt.tmin)then
! shift past t1 and compute new g1,t1 
               tp = t1-dt
               width_in = 10.0**tp
               new_width = .TRUE.
               call tglf_LS
               if(ibranch_in.eq.0)then
                 branch = 1
                 if(gamma_out(2).gt.gamma_out(1))branch = 2
               endif
               tm = t1
               gm = g1
               g1 = gamma_out(branch)
               t1 = tp
! compute new g2,t2
               tp = tm+dt
               width_in = 10.0**tp
               new_width = .TRUE.
               call tglf_LS
               if(ibranch_in.eq.0)then
                 branch = 1
                 if(gamma_out(2).gt.gamma_out(1))branch = 2
               endif
               g2 = gamma_out(branch)
               t2 = tp
             else   ! t1 at tmin
! shrink towards t1
               tp = t1 + dt
               width_in = 10.0**tp
               new_width = .TRUE.
               call tglf_LS
               if(ibranch_in.eq.0)then
                 branch = 1
                 if(gamma_out(2).gt.gamma_out(1))branch = 2
               endif
               g2 = gm
               t2 = tm
               gm = gamma_out(branch)
               tm = tp               
             endif
           elseif(g2.eq.gmax)then
             if(t2.lt.tmax)then
! shift past t2 and compute new g2,t2
               tp = t2 + dt
               width_in = 10.0**tp
               new_width = .TRUE.
               call tglf_LS
               if(ibranch_in.eq.0)then
                 branch = 1
                 if(gamma_out(2).gt.gamma_out(1))branch = 2
               endif
               gm = g2
               tm = t2
               g2 = gamma_out(branch)
               t2 = tp
! compute new g1,t1
               tp = tm - dt
               width_in = 10.0**tp
               new_width = .TRUE.
               call tglf_LS
               if(ibranch_in.eq.0)then
                 branch = 1
                 if(gamma_out(2).gt.gamma_out(1))branch = 2
               endif
               g1 = gamma_out(branch)
               t1 = tp
             else  ! t2 at tmax 
! shrink towards t2
               tp = t2 - dt
               width_in = 10.0**tp
               new_width = .TRUE.
               call tglf_LS
               if(ibranch_in.eq.0)then
                 branch = 1
                 if(gamma_out(2).gt.gamma_out(1))branch = 2
               endif
               g1 = gm
               t1 = tm
               gm = gamma_out(branch)
               tm = tp               
             endif
           else  ! gm.eq.gmax
! compute new g1,t1 and g2,t2 closer to gm,tm
             tp = tm - dt
             width_in = 10.0**tp
             new_width = .TRUE.
             call tglf_LS
             if(ibranch_in.le.0)then
               branch = 1
               if(gamma_out(2).gt.gamma_out(1))branch = 2
             endif
             g1 = gamma_out(branch)
             t1 = tp
!
             tp = tm + dt
             width_in = 10.0**tp
             new_width = .TRUE.
             call tglf_LS
             if(ibranch_in.eq.0)then
               branch = 1
               if(gamma_out(2).gt.gamma_out(1))branch = 2
             endif
             g2 = gamma_out(branch)
             t2 = tp
           endif           
         enddo  ! end of bisection search main loop
! find final maximum
         gmax=gm
         tp=tm
         if(g1.gt.gmax)then
          gmax=g1
          tp=t1
         endif
         if(g2.gt.gmax)then
          gmax=g2
          tp=t2
         endif
         gamma_max = gmax
         width_in = 10.0**tp        
      endif ! done with bisection search
!      write(*,*)"gamma_max=",gamma_max,"width_in=",width_in
       gamma_nb_min_out = gamma_max
!
       if(gamma_max.ne.0.0)then
! refine eigenvalue with more basis functions
         nbasis = save_nbasis
!         write(*,*)"nbasis=",nbasis
         iflux_in=save_iflux
         new_width=.TRUE.
         call tglf_LS
! check for inward ballooning modes
!         write(*,*)"modB_test = ",modB_test
         if(inboard_detrapped_in.ne.0.and.ft_test.gt.modB_test)then
           ft2 = 1.0-ft_test*(1.0-ft*ft)
           if(ft2.lt.0.0)then
             ft = 0.0
           else
!             ft = SQRT(ft2)
             ft =  0.0
           endif
!           write(*,*)"changed ft",ft
           new_geometry = .FALSE.
           new_width = .FALSE.
           new_matrix = .TRUE.
           call tglf_LS
         endif
         if(alpha_quench_in.eq.0.0)then
           if(save_vexb_shear.ne.0.0.or.wgp_max.ne.0.0)then
             do i=1,nmodes_out
               gamma_reference_kx0(i) = gamma_out(i)
               freq_reference_kx0(i) = freq_out(i)
             enddo
             vexb_shear_in = save_vexb_shear
             alpha_kx_p_in = save_alpha_kx_p
             iflux_in=save_iflux
             new_width=.TRUE.
             call tglf_LS
           endif
         endif
         if(ibranch_in.eq.0)then
            branch = 1
            if(gamma_out(2).gt.gamma_out(1))branch = 2
         endif
         gamma_max = gamma_out(branch)
!
!        write(*,*)"width_p_max = ", width_p_max,width_in,gamma_max
!        write(*,*)ky,width_in,gamma_out(1),freq_out(1)
!        write(*,*)" maximum gamma for nbasis = ",nbasis_max_in
!        write(*,*)"gamma_out(1) = ",gamma_out(1)
!        write(*,*)"freq_out(1) = ",freq_out(1)
!        write(*,*)"gamma_out(2) = ",gamma_out(2)
!        write(*,*)"freq_out(2) = ",freq_out(2)
!        gamma_n(igamma)=gamma_out(branch)
!        freq_n(igamma) = freq_out(branch)
       endif
!
       if(gamma_max.eq.0.0)then
         width_in=save_width
         do i=1,nmodes_in
          gamma_out(i)=0.0
          freq_out(i)=0.0
          gamma_reference_kx0(i)=0.0
          freq_reference_kx0(i)=0.0
         enddo
       endif
!
       nbasis=save_nbasis
       iflux_in=save_iflux
       vexb_shear_in = save_vexb_shear
       alpha_kx_p_in = save_alpha_kx_p
!
      END SUBROUTINE tglf_max   
