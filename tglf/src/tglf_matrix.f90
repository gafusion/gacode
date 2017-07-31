!         
      SUBROUTINE get_matrix
!
      USE tglf_dimensions
      USE tglf_global
      USE tglf_coeff
!
      IMPLICIT NONE
      INTEGER :: i,j,k 
!
      new_matrix=.FALSE.
!      write(*,*)"get_matrix"
!
      call FLR_xgrid
!
      call ave_theta
!
!      write(*,*)"ave_bp(1,1)=",ave_bp(1,1)
      if(vpar_model_in.eq.0.and.use_bper_in)then
        call ave_inv0(ave_bp,ave_bpinv)
        do i=1,nbasis
        do j=1,nbasis
          ave_p0inv(i,j) = 0.0
          ave_b0inv(i,j) = 0.0
          do k=1,nbasis
            ave_p0inv(i,j) = ave_p0inv(i,j)+ave_bpinv(i,k)*ave_b0(k,j)
            ave_b0inv(i,j) = ave_b0inv(i,j)+ave_bpinv(i,k)*ave_p0(k,j)
          enddo
        enddo
        enddo
!          write(*,*)"ave_b0inv(i,j)=",ave_b0inv(1,1)
      else
        call ave_inv0(ave_p0,ave_p0inv)
        if(use_bper_in)call ave_inv0(ave_b0,ave_b0inv)
!          write(*,*)"ave_b0inv(i,j)=",ave_b0inv(1,1)
      endif
!
! debug
!      do i=1,nbasis
!      do j=1,nbasis
!        p0 = 0.0
!        do k=1,nbasis
!          p0 = p0 + ave_b0(i,k)*ave_b0inv(k,j)
!        enddo
!        write(*,*)"check b0/b0",i,j,p0
!      enddo
!      enddo
! debug
!
!    compute ave_modwd
!
      call modwd
!       write(*,*)"ave_modwd(1,1)=",ave_modwd(1,1)
!
!    compute ave_modkpar
!
      call modkpar
!
!   compute the  average bessel functions
!
      call ave_h
! debug
!
!       write(*,*)"check hp1"
!       do is = ns0,ns
!        do i=1,nbasis
!        do j=1,nbasis
!         hp1 = 0.0
!         do k=1,nbasis
!           hp1 = hp1 + ave_ht1(is,i,k)*ave_hn(is,k,j)
!         enddo
!         write(*,*)is,i,j,ave_hp1(is,i,j)
!        enddo
!        enddo
!       enddo
!
       call ave_inv(ave_hn,ave_hninv)
       call ave_inv(ave_hp1,ave_hp1inv)
       call ave_inv(ave_hp3,ave_hp3inv)
!
       call h_ratios
!
       if(Linsker_factor_in.ne.0.0)call grad_ave_h
! debug
!       write(*,*)"gradhp1 = "
!       do is=ns0,ns
!       do i=1,nbasis
!       do j=1,nbasis
!         write(*,*)is,i,j,ave_gradhp1(is,i,j)
!       enddo
!       enddo
!       enddo
!
       call ave_hp0
       if(betae_in.gt.0.0)then
         call ave_hb0
         if(vpar_model_in.eq.0)call ave_hbp
       endif
!
       call wd_h
!
       call kpar_h
!
       if(gradB_factor_in.ne.0.0)call gradB_h
!
       if(ft.gt.ft_min)then
         call ave_g
! debug
!       write(*,*)"check gp1"
!       do is = 1,ns
!        do i=1,nbasis
!        do j=1,nbasis
!         gp1 = 0.0
!         do k=1,nbasis
!           gp1 = gp1 + ave_gt1(is,i,k)*ave_gn(is,k,j)
!         enddo
!         write(*,*)is,i,j,ave_gp1(is,i,j)-gp1
!        enddo
!        enddo
!       enddo
!
         call ave_inv(ave_gn,ave_gninv)
         call ave_inv(ave_gp1,ave_gp1inv) 
         call ave_inv(ave_gp3,ave_gp3inv)
!
         call g_ratios
!
         if(Linsker_factor_in.ne.0.0)call grad_ave_g
!
         call ave_gp0
         if(betae_in.gt.0.0)then
           call ave_gb0
           if(vpar_model_in.eq.0)call ave_gbp
         endif
!
         call wd_g
!
         call kpar_g
!
         if(gradB_factor_in.ne.0.0)call gradB_g
!
      endif
!
      END SUBROUTINE get_matrix
!
      SUBROUTINE ave_h
!***************************************************************
!
!   compute the average h-bessel functions
!
!***************************************************************
      USE tglf_dimensions
      USE tglf_hermite
      USE tglf_xgrid
      USE tglf_coeff
!
      IMPLICIT NONE
      INTEGER :: i,j,k,is


      REAL :: ww, zero_cut
!
      zero_cut = 1.E-12
!
      do is = ns0,ns   
       do i=1,nbasis
       do j=i,nbasis
        ave_hn(is,i,j) = 0.0
        ave_hp1(is,i,j) = 0.0
        ave_hp3(is,i,j) = 0.0
        ave_hr11(is,i,j) = 0.0
        ave_hr13(is,i,j) = 0.0
        ave_hr33(is,i,j) = 0.0
        ave_hw113(is,i,j) = 0.0
        ave_hw133(is,i,j) = 0.0
        ave_hw333(is,i,j) = 0.0
!
        do k=1,nx
         hn = hxn(is,k)
         hp1 = hxp1(is,k)
         hp3 = hxp3(is,k)
         hr11 = hxr11(is,k)
         hr13 = hxr13(is,k)
         hr33 = hxr33(is,k)
         hw113 = hxw113(is,k)
         hw133 = hxw133(is,k)
         hw333 = hxw333(is,k)
!
         ww = wx(k)*h(i,k)*h(j,k)
!        write(*,*)i,j,k,ww,wx(k),h(i,k),h(j,k)
         ave_hn(is,i,j)   = ave_hn(is,i,j)  + hn*ww
         ave_hp1(is,i,j)  = ave_hp1(is,i,j)  + hp1*ww
         ave_hp3(is,i,j)  = ave_hp3(is,i,j)  + hp3*ww
         ave_hr11(is,i,j) = ave_hr11(is,i,j) + hr11*ww
         ave_hr13(is,i,j) = ave_hr13(is,i,j) + hr13*ww
         ave_hr33(is,i,j) = ave_hr33(is,i,j) + hr33*ww
         ave_hw113(is,i,j) = ave_hw113(is,i,j) + hw113*ww
         ave_hw133(is,i,j) = ave_hw133(is,i,j) + hw133*ww
         ave_hw333(is,i,j) = ave_hw333(is,i,j) + hw333*ww
        enddo
        if(ABS(ave_hn(is,i,j)).lt.zero_cut)ave_hn(is,i,j)=0.0
        if(ABS(ave_hp1(is,i,j)).lt.zero_cut)ave_hp1(is,i,j)=0.0
        if(ABS(ave_hp3(is,i,j)).lt.zero_cut)ave_hp3(is,i,j)=0.0
        if(ABS(ave_hr11(is,i,j)).lt.zero_cut)ave_hr11(is,i,j)=0.0
        if(ABS(ave_hr13(is,i,j)).lt.zero_cut)ave_hr13(is,i,j)=0.0
        if(ABS(ave_hr33(is,i,j)).lt.zero_cut)ave_hr33(is,i,j)=0.0
        if(ABS(ave_hw113(is,i,j)).lt.zero_cut)ave_hw113(is,i,j)=0.0
        if(ABS(ave_hw133(is,i,j)).lt.zero_cut)ave_hw133(is,i,j)=0.0
        if(ABS(ave_hw333(is,i,j)).lt.zero_cut)ave_hw333(is,i,j)=0.0
! symmetrize
         ave_hn(is,j,i)   = ave_hn(is,i,j) 
         ave_hp1(is,j,i)  = ave_hp1(is,i,j)
         ave_hp3(is,j,i)  = ave_hp3(is,i,j) 
         ave_hr11(is,j,i) = ave_hr11(is,i,j) 
         ave_hr13(is,j,i) = ave_hr13(is,i,j)
         ave_hr33(is,j,i) = ave_hr33(is,i,j)
         ave_hw113(is,j,i) = ave_hw113(is,i,j) 
         ave_hw133(is,j,i) = ave_hw133(is,i,j)
         ave_hw333(is,j,i) = ave_hw333(is,i,j)
       enddo
       enddo
      enddo
!
      END SUBROUTINE ave_h
!
      SUBROUTINE grad_ave_h
!***************************************************************
!
!   compute the parallel gradient of ave_m matrix
!
!***************************************************************
      USE tglf_dimensions
      USE tglf_coeff
!
      IMPLICIT NONE
      INTEGER :: i,j,k,is
!
      do is = ns0,ns
       do i=1,nbasis
       do j=1,nbasis
        gradhp1 = 0.0
!        gradhr11 = 0.0
!        gradhr13 = 0.0
        do k=1,nbasis
         gradhp1 = gradhp1 + &
         ave_kpar(i,k)*ave_hp1(is,k,j) - ave_hp1(is,i,k)*ave_kpar(k,j) 
!         gradhr11 = gradhr11 + &
!         ave_kpar(i,k)*ave_hr11(is,k,j) - ave_hr11(is,i,k)*ave_kpar(k,j)
!         gradhr13 = gradhr13 + &
!         ave_kpar(i,k)*ave_hr13(is,k,j) - ave_hr13(is,i,k)*ave_kpar(k,j)
        enddo
        ave_gradhp1(is,i,j) = gradhp1    
!        ave_gradhr11(is,i,j) = gradhr11
!        ave_gradhr13(is,i,j) = gradhr13
       enddo
       enddo
       do i=1,nbasis
       do j=1,nbasis
        gradhr11 = 0.0
        gradhr13 = 0.0
        do k=1,nbasis
         gradhr11 = gradhr11 + &
         ave_hu1(is,i,k)*ave_gradhp1(is,k,j) 
         gradhr13 = gradhr13 + &
         ave_hu3(is,i,k)*ave_gradhp1(is,k,j)
        enddo
        ave_gradhr11(is,i,j) = gradhr11    
        ave_gradhr13(is,i,j) = gradhr13 
       enddo
       enddo
       do i=1,nbasis
       do j=1,nbasis
         gradhp1p1 = 0.0
         gradhr11p1 = 0.0
         gradhr13p1 = 0.0
         gradhp1 = 0.0
         gradhr11 = 0.0
         gradhr13 = 0.0
         do k=1,nbasis
            gradhp1p1 = gradhp1p1 &
            + ave_gradhp1(is,i,k)*ave_hp1inv(is,k,j)
            gradhr11p1 = gradhr11p1 &
            + ave_gradhr11(is,i,k)*ave_hp1inv(is,k,j)
            gradhr13p1 = gradhr13p1 &
            + ave_gradhr13(is,i,k)*ave_hp1inv(is,k,j)
            gradhp1 = gradhp1 + ave_gradhp1(is,i,k)*ave_p0inv(k,j)
            gradhr11 = gradhr11 + ave_gradhr11(is,i,k)*ave_p0inv(k,j)
            gradhr13 = gradhr13 + ave_gradhr13(is,i,k)*ave_p0inv(k,j)
         enddo
         ave_gradhp1p1(is,i,j) = gradhp1p1
         ave_gradhr11p1(is,i,j) = gradhr11p1
         ave_gradhr13p1(is,i,j) = gradhr13p1
         ave_gradhp1p0(is,i,j) = gradhp1
         ave_gradhr11p0(is,i,j) = gradhr11
         ave_gradhr13p0(is,i,j) = gradhr13
       enddo
       enddo
      enddo   
!
      END SUBROUTINE grad_ave_h
!
      SUBROUTINE h_ratios
!***************************************************************
!
!   compute the h_ratios
!
!***************************************************************
      USE tglf_dimensions
      USE tglf_coeff
!
      IMPLICIT NONE
      INTEGER :: i,j,k,is
!
! compute matrix ratios
!
      do is=ns0,ns
       do i=1,nbasis
       do j=1,nbasis
          ht1 = 0.0
          ht3 = 0.0
          hu1 = 0.0
          hu3 = 0.0
          hu33 = 0.0
          do k=1,nbasis
            ht1   = ht1  + ave_hp1(is,i,k)*ave_hninv(is,k,j)
            ht3   = ht3  + ave_hp3(is,i,k)*ave_hninv(is,k,j)
            hu1   = hu1  + ave_hr11(is,i,k)*ave_hp1inv(is,k,j)
            hu3   = hu3  + ave_hr13(is,i,k)*ave_hp1inv(is,k,j)
            hu33  = hu33 + ave_hr33(is,i,k)*ave_hp3inv(is,k,j)
          enddo
          ave_ht1(is,i,j)  = ht1
          ave_ht3(is,i,j)  = ht3
          ave_hu1(is,i,j)  = hu1
          ave_hu3(is,i,j)  = hu3
          ave_hu33(is,i,j) = hu33
       enddo
       enddo
      enddo          
      do is=ns0,ns
       do i=1,nbasis
       do j=1,nbasis
          hu3ht1 = 0.0
          hu3ht3 = 0.0
          hu33ht1 = 0.0
          hu33ht3 = 0.0
          do k=1,nbasis
            hu3ht1   = hu3ht1  + ave_hu3(is,i,k)*ave_ht1(is,k,j)
            hu33ht1  = hu33ht1 + ave_hu33(is,i,k)*ave_ht1(is,k,j)
            hu3ht3   = hu3ht3  + ave_hu3(is,i,k)*ave_ht3(is,k,j)
            hu33ht3  = hu33ht3 + ave_hu33(is,i,k)*ave_ht3(is,k,j)
          enddo
          ave_hu3ht1(is,i,j)  = hu3ht1
          ave_hu3ht3(is,i,j)  = hu3ht3
          ave_hu33ht1(is,i,j)  = hu33ht1
          ave_hu33ht3(is,i,j)  = hu33ht3
       enddo
       enddo
      enddo          
!
      END SUBROUTINE h_ratios
!
      SUBROUTINE ave_hp0
!***************************************************************
!
!   compute the products ave_h*ave_p0inv
!
!***************************************************************
      USE tglf_dimensions
      USE tglf_coeff
!
      IMPLICIT NONE
      INTEGER :: i,j,k,is
!
! compute matrix products
!
      do is=ns0,ns
       do i=1,nbasis
       do j=1,nbasis
          hn = 0.0
          hp1 = 0.0
          hp3 = 0.0
          hr11 = 0.0
          hr13 = 0.0
          hr33 = 0.0
          do k=1,nbasis
            hn    = hn    + ave_hn(is,i,k)*ave_p0inv(k,j)
            hp1   = hp1   + ave_hp1(is,i,k)*ave_p0inv(k,j)
            hp3   = hp3   + ave_hp3(is,i,k)*ave_p0inv(k,j)
            hr11  = hr11  + ave_hr11(is,i,k)*ave_p0inv(k,j)
            hr13  = hr13  + ave_hr13(is,i,k)*ave_p0inv(k,j)
            hr33  = hr33  + ave_hr33(is,i,k)*ave_p0inv(k,j)
          enddo
          ave_hnp0(is,i,j)      = hn
          ave_hp1p0(is,i,j)     = hp1
          ave_hp3p0(is,i,j)     = hp3
          ave_hr11p0(is,i,j)    = hr11
          ave_hr13p0(is,i,j)    = hr13
          ave_hr33p0(is,i,j)    = hr33
       enddo
       enddo
      enddo          
!
      do is = ns0,ns
        do i=1,nbasis
        do j=1,nbasis
        hp1 = 0.0
        hr11 = 0.0
        hr13 = 0.0
        do k=1,nbasis
          hp1 = hp1 + ave_c_tor_par(i,k)*ave_hp1p0(is,k,j)
          hr11 = hr11 + ave_c_tor_par(i,k)*ave_hr11p0(is,k,j)
          hr13 = hr13 + ave_c_tor_par(i,k)*ave_hr13p0(is,k,j)
        enddo
        ave_c_tor_par_hp1p0(is,i,j) = hp1
        ave_c_tor_par_hr11p0(is,i,j) = hr11        
        ave_c_tor_par_hr13p0(is,i,j) = hr13        
        enddo
        enddo
      enddo
!
      END SUBROUTINE ave_hp0
!
      SUBROUTINE ave_hb0
!***************************************************************
!
!   compute the products ave_h*ave_b0inv
!
!***************************************************************
      USE tglf_dimensions
      USE tglf_coeff
!
      IMPLICIT NONE
      INTEGER :: i,j,k,is
!
! compute matrix products
!
      do is=ns0,ns
       do i=1,nbasis
       do j=1,nbasis
          hn    = 0.0
          hp1   = 0.0
          hp3   = 0.0
          hr11  = 0.0
          hr13  = 0.0
          hr33  = 0.0
          hw113 = 0.0
          hw133 = 0.0
          hw333 = 0.0
          do k=1,nbasis
            hn     = hn     + ave_hn(is,i,k)*ave_b0inv(k,j)
            hp1    = hp1    + ave_hp1(is,i,k)*ave_b0inv(k,j)
            hp3    = hp3    + ave_hp3(is,i,k)*ave_b0inv(k,j)
            hr11   = hr11   + ave_hr11(is,i,k)*ave_b0inv(k,j)
            hr13   = hr13   + ave_hr13(is,i,k)*ave_b0inv(k,j)
            hr33   = hr33   + ave_hr33(is,i,k)*ave_b0inv(k,j)
            hw113  = hw113  + ave_hw113(is,i,k)*ave_b0inv(k,j)
            hw133  = hw133  + ave_hw133(is,i,k)*ave_b0inv(k,j)
            hw333  = hw333  + ave_hw333(is,i,k)*ave_b0inv(k,j)
          enddo
          ave_hnb0(is,i,j)     = hn
          ave_hp1b0(is,i,j)    = hp1
          ave_hp3b0(is,i,j)    = hp3
          ave_hr11b0(is,i,j)   = hr11
          ave_hr13b0(is,i,j)   = hr13
          ave_hr33b0(is,i,j)   = hr33
          ave_hw113b0(is,i,j)  = hw113
          ave_hw133b0(is,i,j)  = hw133
          ave_hw333b0(is,i,j)  = hw333
       enddo
       enddo
      enddo          
!
      END SUBROUTINE ave_hb0
!
      SUBROUTINE ave_hbp
!***************************************************************
!
!   compute the products ave_h*ave_bpinv
!
!***************************************************************
      USE tglf_dimensions
      USE tglf_coeff
!
      IMPLICIT NONE
      INTEGER :: i,j,k,is
!
! compute matrix products
!
      do is=ns0,ns
       do i=1,nbasis
       do j=1,nbasis
          hn    = 0.0
          hp1   = 0.0
          hp3   = 0.0
          hr11  = 0.0
          hr13  = 0.0
          hr33  = 0.0
          hw113 = 0.0
          hw133 = 0.0
          hw333 = 0.0
          do k=1,nbasis
            hn     = hn     + ave_hn(is,i,k)*ave_bpinv(k,j)
            hp1    = hp1    + ave_hp1(is,i,k)*ave_bpinv(k,j)
            hp3    = hp3    + ave_hp3(is,i,k)*ave_bpinv(k,j)
            hr11   = hr11   + ave_hr11(is,i,k)*ave_bpinv(k,j)
            hr13   = hr13   + ave_hr13(is,i,k)*ave_bpinv(k,j)
            hr33   = hr33   + ave_hr33(is,i,k)*ave_bpinv(k,j)
            hw113  = hw113  + ave_hw113(is,i,k)*ave_bpinv(k,j)
            hw133  = hw133  + ave_hw133(is,i,k)*ave_bpinv(k,j)
            hw333  = hw333  + ave_hw333(is,i,k)*ave_bpinv(k,j)
          enddo
          ave_hnbp(is,i,j)     = hn
          ave_hp1bp(is,i,j)    = hp1
          ave_hp3bp(is,i,j)    = hp3
          ave_hr11bp(is,i,j)   = hr11
          ave_hr13bp(is,i,j)   = hr13
          ave_hr33bp(is,i,j)   = hr33
          ave_hw113bp(is,i,j)  = hw113
          ave_hw133bp(is,i,j)  = hw133
          ave_hw333bp(is,i,j)  = hw333
       enddo
       enddo
      enddo          
!
      END SUBROUTINE ave_hbp
!
      SUBROUTINE gradB_h
!***************************************************************
!
!   compute the products gradB*h
!
!***************************************************************
      USE tglf_dimensions
      USE tglf_coeff
!
      IMPLICIT NONE
      INTEGER :: i,j,k,is
!
! compute matrix ratios
!
      do is=ns0,ns
       do i=1,nbasis
       do j=1,nbasis
          gradBhp1 = 0.0
          gradBhp3 = 0.0
          gradBhr11 = 0.0
          gradBhr13 = 0.0
          gradBhr33 = 0.0
          gradBhu1 = 0.0
          gradBhu3 = 0.0
          gradBhu33 = 0.0
          do k=1,nbasis
            gradBhp1  = gradBhp1 + ave_gradB(i,k)*ave_hp1p0(is,k,j)
            gradBhp3  = gradBhp3 + ave_gradB(i,k)*ave_hp3p0(is,k,j)
            gradBhr11  = gradBhr11 + ave_gradB(i,k)*ave_hr11p0(is,k,j)
            gradBhr13  = gradBhr13 + ave_gradB(i,k)*ave_hr13p0(is,k,j)
            gradBhr33  = gradBhr33 + ave_gradB(i,k)*ave_hr33p0(is,k,j)
            gradBhu1  = gradBhu1 + ave_gradB(i,k)*ave_hu1(is,k,j)
            gradBhu3  = gradBhu3 + ave_gradB(i,k)*ave_hu3(is,k,j)
            gradBhu33  = gradBhu33 + ave_gradB(i,k)*ave_hu33(is,k,j)
          enddo
          ave_gradBhp1(is,i,j) = gradBhp1
          ave_gradBhp3(is,i,j) = gradBhp3
          ave_gradBhr11(is,i,j) = gradBhr11
          ave_gradBhr13(is,i,j) = gradBhr13
          ave_gradBhr33(is,i,j) = gradBhr33
          ave_gradBhu1(is,i,j) = gradBhu1
          ave_gradBhu3(is,i,j) = gradBhu3
          ave_gradBhu33(is,i,j) = gradBhu33
       enddo
       enddo
      enddo          
!
      END SUBROUTINE gradB_h
!
      SUBROUTINE ave_g
!***************************************************************
!
!   compute the average g- bessel functions
!
!***************************************************************
      USE tglf_dimensions
      USE tglf_hermite
      USE tglf_xgrid
      USE tglf_coeff
!
      IMPLICIT NONE
      INTEGER :: i,j,k,is
      REAL :: ww, zero_cut
!
      zero_cut=1.E-12
!
      do is=ns0,ns
       do i=1,nbasis
       do j=i,nbasis
!  initialize the bessel function averages
        ave_gn(is,i,j) = 0.0
        ave_gp1(is,i,j) = 0.0
        ave_gp3(is,i,j) = 0.0
        ave_gr11(is,i,j) = 0.0
        ave_gr13(is,i,j) = 0.0
        ave_gr33(is,i,j) = 0.0
        ave_gw113(is,i,j) = 0.0
        ave_gw133(is,i,j) = 0.0
        ave_gw333(is,i,j) = 0.0
!
        do k=1,nx
         gn = gxn(is,k)
         gp1 = gxp1(is,k)
         gp3 = gxp3(is,k)
         gr11 = gxr11(is,k)
         gr13 = gxr13(is,k)
         gr33 = gxr33(is,k)
         gw113 = gxw113(is,k)
         gw133 = gxw133(is,k)
         gw333 = gxw333(is,k)
!
         ww = wx(k)*h(i,k)*h(j,k)
         ave_gn(is,i,j) = ave_gn(is,i,j) + gn*ww
         ave_gp1(is,i,j) = ave_gp1(is,i,j) + gp1*ww
         ave_gp3(is,i,j) = ave_gp3(is,i,j) + gp3*ww
         ave_gr11(is,i,j) = ave_gr11(is,i,j) + gr11*ww
         ave_gr13(is,i,j) = ave_gr13(is,i,j) + gr13*ww
         ave_gr33(is,i,j) = ave_gr33(is,i,j) + gr33*ww
         ave_gw113(is,i,j) = ave_gw113(is,i,j) + gw113*ww
         ave_gw133(is,i,j) = ave_gw133(is,i,j) + gw133*ww
         ave_gw333(is,i,j) = ave_gw333(is,i,j) + gw333*ww
        enddo
        if(ABS(ave_gn(is,i,j)).lt.zero_cut)ave_gn(is,i,j)=0.0
        if(ABS(ave_gp1(is,i,j)).lt.zero_cut)ave_gp1(is,i,j)=0.0
        if(ABS(ave_gp3(is,i,j)).lt.zero_cut)ave_gp3(is,i,j)=0.0
        if(ABS(ave_gr11(is,i,j)).lt.zero_cut)ave_gr11(is,i,j)=0.0
        if(ABS(ave_gr13(is,i,j)).lt.zero_cut)ave_gr13(is,i,j)=0.0
        if(ABS(ave_gr33(is,i,j)).lt.zero_cut)ave_gr33(is,i,j)=0.0
        if(ABS(ave_gw113(is,i,j)).lt.zero_cut)ave_gw113(is,i,j)=0.0
        if(ABS(ave_gw133(is,i,j)).lt.zero_cut)ave_gw133(is,i,j)=0.0
        if(ABS(ave_gw333(is,i,j)).lt.zero_cut)ave_gw333(is,i,j)=0.0
! symmetrize
         ave_gn(is,j,i) = ave_gn(is,i,j) 
         ave_gp1(is,j,i) = ave_gp1(is,i,j)
         ave_gp3(is,j,i) = ave_gp3(is,i,j) 
         ave_gr11(is,j,i) = ave_gr11(is,i,j) 
         ave_gr13(is,j,i) = ave_gr13(is,i,j)
         ave_gr33(is,j,i) = ave_gr33(is,i,j)
         ave_gw113(is,j,i) = ave_gw113(is,i,j) 
         ave_gw133(is,j,i) = ave_gw133(is,i,j)
         ave_gw333(is,j,i) = ave_gw333(is,i,j)
       enddo
       enddo
      enddo
!
      END SUBROUTINE ave_g
!
      SUBROUTINE grad_ave_g
!***************************************************************
!
!   compute the parallel gradient of ave_g matricies
!
!***************************************************************
      USE tglf_dimensions
      USE tglf_coeff
!
      IMPLICIT NONE
      INTEGER :: i,j,k,is
!
      do is = ns0,ns
       do i=1,nbasis
       do j=1,nbasis
        gradgp1 = 0.0
!        gradgr11 = 0.0
!        gradgr13 = 0.0
        do k=1,nbasis
         gradgp1 = gradgp1 + &
         ave_kpar(i,k)*ave_gp1(is,k,j) - ave_gp1(is,i,k)*ave_kpar(k,j) 
!         gradgr11 = gradgr11 + &
!         ave_kpar(i,k)*ave_gr11(is,k,j) - ave_gr11(is,i,k)*ave_kpar(k,j)
!         gradgr13 = gradgr13 + &
!         ave_kpar(i,k)*ave_gr13(is,k,j) - ave_gr13(is,i,k)*ave_kpar(k,j)
        enddo
        ave_gradgp1(is,i,j) = gradgp1    
        ave_gradgr11(is,i,j) = gradgr11    
        ave_gradgr13(is,i,j) = gradgr13 
       enddo
       enddo
       do i=1,nbasis
       do j=1,nbasis
        gradgr11 = 0.0
        gradgr13 = 0.0
        do k=1,nbasis
         gradgr11 = gradgr11 + &
         ave_gu1(is,i,k)*ave_gradgp1(is,k,j) 
         gradgr13 = gradgr13 + &
         ave_gu3(is,i,k)*ave_gradgp1(is,k,j)
        enddo
        ave_gradgr11(is,i,j) = gradgr11    
        ave_gradgr13(is,i,j) = gradgr13 
       enddo
       enddo
       do i=1,nbasis
       do j=1,nbasis
         gradgp1p1 = 0.0
         gradgr11p1 = 0.0
         gradgr13p1 = 0.0
         gradgp1 = 0.0
         gradgr11 = 0.0
         gradgr13 = 0.0
         do k=1,nbasis
            gradgp1p1 = gradgp1p1 &
            + ave_gradgp1(is,i,k)*ave_gp1inv(is,k,j)
            gradgr11p1 = gradgr11p1 &
            + ave_gradgr11(is,i,k)*ave_gp1inv(is,k,j)
            gradgr13p1 = gradgr13p1 &
            + ave_gradgr13(is,i,k)*ave_gp1inv(is,k,j)
            gradgp1  = gradgp1  + ave_gradgp1(is,i,k)*ave_p0inv(k,j)
            gradgr11  = gradgr11  + ave_gradgr11(is,i,k)*ave_p0inv(k,j)
            gradgr13  = gradgr13  + ave_gradgr13(is,i,k)*ave_p0inv(k,j)
         enddo
         ave_gradgp1p1(is,i,j) = gradgp1p1
         ave_gradgr11p1(is,i,j) = gradgr11p1
         ave_gradgr13p1(is,i,j) = gradgr13p1
         ave_gradgp1p0(is,i,j) = gradgp1
         ave_gradgr11p0(is,i,j) = gradgr11
         ave_gradgr13p0(is,i,j) = gradgr13
       enddo
       enddo
      enddo   
      return
      END SUBROUTINE grad_ave_g
!
      SUBROUTINE g_ratios
!***************************************************************
!
!   compute the g_ratios
!
!***************************************************************
      USE tglf_dimensions
      USE tglf_coeff
!
      IMPLICIT NONE
      INTEGER :: i,j,k,is
!
! compute matrix ratios
!
      do is=ns0,ns
       do i=1,nbasis
       do j=1,nbasis
          gt1 = 0.0
          gt3 = 0.0
          gu1 = 0.0
          gu3 = 0.0
          gu33 = 0.0
          do k=1,nbasis
            gt1   = gt1  + ave_gp1(is,i,k)*ave_gninv(is,k,j)
            gt3   = gt3  + ave_gp3(is,i,k)*ave_gninv(is,k,j)
            gu1   = gu1  + ave_gr11(is,i,k)*ave_gp1inv(is,k,j)
            gu3   = gu3  + ave_gr13(is,i,k)*ave_gp1inv(is,k,j)
            gu33  = gu33 + ave_gr33(is,i,k)*ave_gp3inv(is,k,j)
          enddo
          ave_gt1(is,i,j)  = gt1
          ave_gt3(is,i,j)  = gt3
          ave_gu1(is,i,j)  = gu1
          ave_gu3(is,i,j)  = gu3
          ave_gu33(is,i,j) = gu33
       enddo
       enddo
      enddo  
!
      do is=ns0,ns
       do i=1,nbasis
       do j=1,nbasis
          gu3gt1 = 0.0
          gu3gt3 = 0.0
          gu33gt1 = 0.0
          gu33gt3 = 0.0
          do k=1,nbasis
            gu3gt1   = gu3gt1  + ave_gu3(is,i,k)*ave_gt1(is,k,j)
            gu33gt1  = gu33gt1 + ave_gu33(is,i,k)*ave_gt1(is,k,j)
            gu3gt3   = gu3gt3  + ave_gu3(is,i,k)*ave_gt3(is,k,j)
            gu33gt3  = gu33gt3 + ave_gu33(is,i,k)*ave_gt3(is,k,j)
          enddo
          ave_gu3gt1(is,i,j)  = gu3gt1
          ave_gu3gt3(is,i,j)  = gu3gt3
          ave_gu33gt1(is,i,j)  = gu33gt1
          ave_gu33gt3(is,i,j)  = gu33gt3
       enddo
       enddo
      enddo                  
!
      END SUBROUTINE g_ratios
!
      SUBROUTINE ave_gp0
!***************************************************************
!
!   compute the average g- bessel functions
!
!***************************************************************
      USE tglf_dimensions
      USE tglf_coeff
!
      IMPLICIT NONE
      INTEGER :: i,j,k,is
!
!   compute matrix products
!
      do is=ns0,ns
       do i=1,nbasis
       do j=1,nbasis
          gn = 0.0
          gp1 = 0.0
          gp3 = 0.0
          gr11 = 0.0
          gr13 = 0.0
          gr33 = 0.0
          do k=1,nbasis
            gn    = gn    + ave_gn(is,i,k)*ave_p0inv(k,j)
            gp1   = gp1   + ave_gp1(is,i,k)*ave_p0inv(k,j)
            gp3   = gp3   + ave_gp3(is,i,k)*ave_p0inv(k,j)
            gr11  = gr11  + ave_gr11(is,i,k)*ave_p0inv(k,j)
            gr13  = gr13  + ave_gr13(is,i,k)*ave_p0inv(k,j)
            gr33  = gr33  + ave_gr33(is,i,k)*ave_p0inv(k,j)
          enddo
          ave_gnp0(is,i,j) = gn
          ave_gp1p0(is,i,j) = gp1
          ave_gp3p0(is,i,j) = gp3
          ave_gr11p0(is,i,j) = gr11
          ave_gr13p0(is,i,j) = gr13
          ave_gr33p0(is,i,j) = gr33
       enddo
       enddo
      enddo
!
      do is = ns0,ns
        do i=1,nbasis
        do j=1,nbasis
        gp1 = 0.0
        gr11 = 0.0
        gr13 = 0.0
        do k=1,nbasis
          gp1 = gp1 + ave_c_tor_par(i,k)*ave_gp1p0(is,k,j)
          gr11 = gr11 + ave_c_tor_par(i,k)*ave_gr11p0(is,k,j)
          gr13 = gr13 + ave_c_tor_par(i,k)*ave_gr13p0(is,k,j)
        enddo
        ave_c_tor_par_gp1p0(is,i,j) = gp1
        ave_c_tor_par_gr11p0(is,i,j) = gr11        
        ave_c_tor_par_gr13p0(is,i,j) = gr13        
        enddo
        enddo
      enddo
!
      END SUBROUTINE ave_gp0
!
      SUBROUTINE ave_gb0
!***************************************************************
!
!   compute the products ave_g*ave_b0inv
!
!***************************************************************
      USE tglf_dimensions
      USE tglf_coeff
!
      IMPLICIT NONE
      INTEGER :: i,j,k,is
!
!   compute matrix products
!
      do is=ns0,ns
       do i=1,nbasis
       do j=1,nbasis
          gn = 0.0
          gp1 = 0.0
          gp3 = 0.0
          gr11 = 0.0
          gr13 = 0.0
          gr33 = 0.0
          gw113 = 0.0
          gw133 = 0.0
          gw333 = 0.0
          do k=1,nbasis
            gn    = gn    + ave_gn(is,i,k)*ave_b0inv(k,j)
            gp1   = gp1   + ave_gp1(is,i,k)*ave_b0inv(k,j)
            gp3   = gp3   + ave_gp3(is,i,k)*ave_b0inv(k,j)
            gr11  = gr11  + ave_gr11(is,i,k)*ave_b0inv(k,j)
            gr13  = gr13  + ave_gr13(is,i,k)*ave_b0inv(k,j)
            gr33  = gr33  + ave_gr33(is,i,k)*ave_b0inv(k,j)
            gw113 = gw113 + ave_gw113(is,i,k)*ave_b0inv(k,j)
            gw133 = gw133 + ave_gw133(is,i,k)*ave_b0inv(k,j)
            gw333 = gw333 + ave_gw333(is,i,k)*ave_b0inv(k,j)
          enddo
          ave_gnb0(is,i,j)    = gn
          ave_gp1b0(is,i,j)   = gp1
          ave_gp3b0(is,i,j)   = gp3
          ave_gr11b0(is,i,j)  = gr11
          ave_gr13b0(is,i,j)  = gr13
          ave_gr33b0(is,i,j)  = gr33
          ave_gw113b0(is,i,j) = gw113
          ave_gw133b0(is,i,j) = gw133
          ave_gw333b0(is,i,j) = gw333
       enddo
       enddo
      enddo
!
      END SUBROUTINE ave_gb0
!!
      SUBROUTINE ave_gbp
!***************************************************************
!
!   compute the products ave_g*ave_bpinv
!
!***************************************************************
      USE tglf_dimensions
      USE tglf_coeff
!
      IMPLICIT NONE
      INTEGER :: i,j,k,is
!
!   compute matrix products
!
      do is=ns0,ns
       do i=1,nbasis
       do j=1,nbasis
          gn = 0.0
          gp1 = 0.0
          gp3 = 0.0
          gr11 = 0.0
          gr13 = 0.0
          gr33 = 0.0
          gw113 = 0.0
          gw133 = 0.0
          gw333 = 0.0
          do k=1,nbasis
            gn    = gn    + ave_gn(is,i,k)*ave_bpinv(k,j)
            gp1   = gp1   + ave_gp1(is,i,k)*ave_bpinv(k,j)
            gp3   = gp3   + ave_gp3(is,i,k)*ave_bpinv(k,j)
            gr11  = gr11  + ave_gr11(is,i,k)*ave_bpinv(k,j)
            gr13  = gr13  + ave_gr13(is,i,k)*ave_bpinv(k,j)
            gr33  = gr33  + ave_gr33(is,i,k)*ave_bpinv(k,j)
            gw113 = gw113 + ave_gw113(is,i,k)*ave_bpinv(k,j)
            gw133 = gw133 + ave_gw133(is,i,k)*ave_bpinv(k,j)
            gw333 = gw333 + ave_gw333(is,i,k)*ave_bpinv(k,j)
          enddo
          ave_gnbp(is,i,j)    = gn
          ave_gp1bp(is,i,j)   = gp1
          ave_gp3bp(is,i,j)   = gp3
          ave_gr11bp(is,i,j)  = gr11
          ave_gr13bp(is,i,j)  = gr13
          ave_gr33bp(is,i,j)  = gr33
          ave_gw113bp(is,i,j) = gw113
          ave_gw133bp(is,i,j) = gw133
          ave_gw333bp(is,i,j) = gw333
       enddo
       enddo
      enddo
!
      END SUBROUTINE ave_gbp
!
      SUBROUTINE gradB_g
!***************************************************************
!
!   compute the products gradB*ave_g
!
!***************************************************************
      USE tglf_dimensions
      USE tglf_coeff
!
      IMPLICIT NONE
      INTEGER :: i,j,k,is

!
! compute gradB*ave_g
!
      do is=ns0,ns
       do i=1,nbasis
       do j=1,nbasis
          gradBgp1 = 0.0
          gradBgp3 = 0.0
          gradBgr11 = 0.0
          gradBgr13 = 0.0
          gradBgr33 = 0.0
          gradBgu1 = 0.0
          gradBgu3 = 0.0
          gradBgu33 = 0.0
          do k=1,nbasis
            gradBgp1  = gradBgp1 + ave_gradB(i,k)*ave_gp1p0(is,k,j)
            gradBgp3  = gradBgp3 + ave_gradB(i,k)*ave_gp3p0(is,k,j)
            gradBgr11  = gradBgr11 + ave_gradB(i,k)*ave_gr11p0(is,k,j)
            gradBgr13  = gradBgr13 + ave_gradB(i,k)*ave_gr13p0(is,k,j)
            gradBgr33  = gradBgr33 + ave_gradB(i,k)*ave_gr33p0(is,k,j)
            gradBgu1  = gradBgu1 + ave_gradB(i,k)*ave_gu1(is,k,j)
            gradBgu3  = gradBgu3 + ave_gradB(i,k)*ave_gu3(is,k,j)
            gradBgu33  = gradBgu33 + ave_gradB(i,k)*ave_gu33(is,k,j)
          enddo
          ave_gradBgp1(is,i,j) = gradBgp1
          ave_gradBgp3(is,i,j) = gradBgp3
          ave_gradBgr11(is,i,j) = gradBgr11
          ave_gradBgr13(is,i,j) = gradBgr13
          ave_gradBgr33(is,i,j) = gradBgr33
          ave_gradBgu1(is,i,j) = gradBgu1
          ave_gradBgu3(is,i,j) = gradBgu3
          ave_gradBgu33(is,i,j) = gradBgu33
       enddo
       enddo
      enddo          
!
      END SUBROUTINE gradB_g
!
      SUBROUTINE ave_theta
!***********************************************************
!  compute  k-independent hermite basis averages
!***********************************************************
      USE tglf_dimensions
      USE tglf_global
      USE tglf_species
      USE tglf_hermite
      USE tglf_xgrid
      USE tglf_coeff
!
      IMPLICIT NONE
      INTEGER :: i,j,k,is
      REAL :: ww
      REAL :: zero_cut
      REAL :: debye,betaU
      REAL :: ft2
!
      zero_cut = 1.E-12
      ft2 = ft*ft
!
! fill the Poisson equation phi multiplier px0 x-grid array
! note that pol is set in get_species
!
      do i=1,nx
        debye = debye_factor_in*b0x(i)*(ky*debye_in)**2  ! (k_per*debye length unit)^2
        p0x(i) = debye + pol
      enddo
!debug
!      write(*,*)"check p0",nx
!      do i=1,nx
!        write(*,*)i,p0x(i)
!      enddo
!debug
!
!  compute the guass-hermite intregrals
!
       do i=1,nbasis
       do j=1,nbasis
         ave_kpar(i,j) = 0.0
       enddo
       enddo
       do i=1,nbasis        
        if(i.lt.nbasis)then
         ave_kpar(i,i+1)=SQRT(REAL(i)/2.0)
! kpar is an odd function
         ave_kpar(i+1,i)=-ave_kpar(i,i+1)
        endif
       do j=i,nbasis
!   initialize the averages
        ave_wdh(i,j) = 0.0
        ave_wdg(i,j)=0.0
        ave_b0(i,j) = 0.0
        ave_lnB(i,j) = 0.0
        ave_p0inv(i,j) = 0.0
        ave_p0(i,j) = 0.0
        ave_kx(i,j) = 0.0
        ave_c_tor_par(i,j) = 0.0
        ave_c_tor_per(i,j) = 0.0
        ave_c_par_par(i,j) = 0.0
!
        do k=1,nx
         ww=wx(k)*h(i,k)*h(j,k)     
         ave_wdh(i,j)    = ave_wdh(i,j)  + ww*wdx(k)
         ave_wdg(i,j)   = ave_wdg(i,j) + ww*(wdx(k)+wdpx(k)*(1.0-ft2)/(1.0+ft2))
         ave_b0(i,j)    = ave_b0(i,j)  + ww*b0x(k)
!         ave_lnB(i,j)   = ave_lnB(i,j) + ww*DLOG(Bx(k))
         ave_p0inv(i,j) = ave_p0inv(i,j) + ww/p0x(k)
         ave_p0(i,j) = ave_p0(i,j) + ww*p0x(k)
         ave_kx(i,j) = ave_kx(i,j) + ww*kxx(k)
         ave_c_tor_par(i,j) = ave_c_tor_par(i,j) + ww*cx_tor_par(k)
         ave_c_tor_per(i,j) = ave_c_tor_per(i,j) + ww*cx_tor_per(k)
         ave_c_par_par(i,j) = ave_c_par_par(i,j) + ww*cx_par_par(k)
        enddo
        if(ABS(ave_wdh(i,j)).lt.zero_cut)ave_wdh(i,j) = 0.0
        if(ABS(ave_wdg(i,j)).lt.zero_cut)ave_wdg(i,j) = 0.0
        if(ABS(ave_b0(i,j)).lt.zero_cut)ave_b0(i,j) = 0.0
        if(ABS(ave_lnB(i,j)).lt.zero_cut)ave_lnB(i,j) = 0.0
        if(ABS(ave_p0inv(i,j)).lt.zero_cut)ave_p0inv(i,j) = 0.0
        if(ABS(ave_p0(i,j)).lt.zero_cut)ave_p0(i,j) = 0.0
        if(ABS(ave_kx(i,j)).lt.zero_cut)ave_kx(i,j) = 0.0
        if(ABS(ave_c_tor_par(i,j)).lt.zero_cut)ave_c_tor_par(i,j) = 0.0
        if(ABS(ave_c_tor_per(i,j)).lt.zero_cut)ave_c_tor_per(i,j) = 0.0
        if(ABS(ave_c_par_par(i,j)).lt.zero_cut)ave_c_par_par(i,j) = 0.0
! symmetrize
        ave_wdh(j,i)    = ave_wdh(i,j)
        ave_wdg(j,i)    = ave_wdg(i,j)
        ave_b0(j,i)    = ave_b0(i,j)
        ave_lnB(j,i)   = ave_lnB(i,j)
        ave_p0inv(j,i) = ave_p0inv(i,j)
        ave_p0(j,i) = ave_p0(i,j)
        ave_kx(j,i) = ave_kx(i,j)
        ave_c_tor_par(j,i) = ave_c_tor_par(i,j)
        ave_c_tor_per(j,i) = ave_c_tor_per(i,j)
        ave_c_par_par(j,i) = ave_c_par_par(i,j)
        enddo
       enddo
!       write(*,*)"ave_wd(1,1)=",ave_wd(1,1)
!       write(*,*)"ave_p0inv(1,1) = ",ave_p0inv(1,1)
!       write(*,*)"ave_p0(1,1) = ",ave_p0(1,1)
!       write(*,*)"ave_b0(1,1) = ",ave_b0(1,1)
!       write(*,*)"ave_kx(1,1) = ",ave_kx(1,1)
!       write(*,*)"ave_tor_par(1,1) = ",ave_c_tor_par(1,1)
!       write(*,*)"ave_tor_per(1,1) = ",ave_c_tor_per(1,1)
!       write(*,*)"ave_par_par(1,1) = ",ave_c_par_par(1,1)
       do i=1,nbasis
       do j=1,nbasis
         gradB = 0.0
         do k=1,nbasis
           gradB = gradB + ave_kpar(i,k)*ave_lnB(k,j) &
           - ave_lnB(i,k)*ave_kpar(k,j)
         enddo
         ave_gradB(i,j) = gradB
!       write(*,*)i,j,"gradB= ",gradB
       enddo
       enddo
!
       do is=ns0,ns
         if(nbasis.eq.1)then
           ave_kpar_eff(is,1,1) = -xi/sqrt_two
           if(vpar_model_in.eq.1)then
             ave_kpar_eff(is,1,1) = ave_kpar_eff(is,1,1)   &
             +xi*ABS(alpha_mach_in*vpar_in(is))*ky*R_unit*q_unit*width_in*mass(is)/zs(is)
           endif
         else
           do i=1,nbasis
           do j=1,nbasis
             ave_kpar_eff(is,i,j) = ave_kpar(i,j)     
             if(vpar_model_in.eq.1.and.i.eq.j)then
               ave_kpar_eff(is,i,j) = ave_kpar_eff(is,i,j)  &
               -xi*alpha_mach_in*vpar_in(is)*ky*R_unit*q_unit*width_in*mass(is)/zs(is)
             endif
           enddo
           enddo
         endif
       enddo
!
       if(vpar_model_in.eq.0.and.use_bper_in)then
         betaU = 0.5*betae_in/(ky*ky)
         betaU = betaU*U0*U0
         do i=1,nbasis
         do j=1,nbasis
           ave_bp(i,j) = 0.0
           do k=1,nbasis
             ave_bp(i,j) = ave_bp(i,j) + ave_b0(i,k)*ave_p0(k,j)
           enddo
           if(i.eq.j)ave_bp(i,j) = ave_bp(i,j) + betaU
!           write(*,*)i,j,"ave_bp(i,j)=",ave_bp(i,j),betaU
         enddo
         enddo
       endif
! 
!
!  debug
!       write(*,*)"check ave_p0inv"
!       do i=1,nbasis
!         write(*,*)"ave_p0inv(i,j) = ",(ave_p0inv(i,j),j=1,nbasis)
!       enddo
!       write(*,*)"check ave_kpar "
!       do i=1,nbasis
!       do j=1,nbasis
!        write(*,*)i,j,ave_kpar(i,j)
!       enddo
!       enddo
!       write(*,*)"check ave_kpar_eff "
!       do is=ns0,ns
!       do i=1,nbasis
!       do j=1,nbasis
!        write(*,*)i,j,ave_kpar_eff(is,i,j)
!       enddo
!       enddo
!       enddo
!       write(*,*)"check ave_wdh "
!       do i=1,nbasis
!       do j=1,nbasis
!        write(*,*)i,j,ave_wdh(i,j)
!       enddo
!       enddo
!       write(*,*)"check ave_wdg "
!       do i=1,nbasis
!       do j=1,nbasis
!        write(*,*)i,j,ave_wdg(i,j)
!       enddo
!       enddo
!       write(*,*)"check ave_gradB "
!       do i=1,nbasis
!       do j=1,nbasis
!        write(*,*)i,j,ave_gradB(i,j)
!       enddo
!       enddo
!
      END SUBROUTINE ave_theta
!
      SUBROUTINE wd_h
!***************************************************************
!
!   compute the products ave_wd*ave_h
!
!***************************************************************
      USE tglf_dimensions
      USE tglf_coeff
      USE tglf_global
!
      IMPLICIT NONE
      INTEGER :: is,ib,jb,k
!
      do is=ns0,ns
      do ib=1,nbasis
      do jb=1,nbasis
       wdhp1p0 = 0.0
       wdhr11p0 = 0.0
       wdhr13p0 = 0.0
       wdhp1b0 = 0.0
       wdhr11b0 = 0.0
       wdhr13b0 =0.0
       wdhp1bp = 0.0
       wdhr11bp = 0.0
       wdhr13bp =0.0
       wdht1 = 0.0
       wdht3 = 0.0
       wdhu1 = 0.0
       wdhu3 = 0.0
       wdhu3ht1 = 0.0
       wdhu3ht3 = 0.0
       wdhu33 = 0.0
       wdhu33ht1 = 0.0
       wdhu33ht3 = 0.0
       modwdhu1 = 0.0
       modwdhu3 = 0.0
       modwdhu33 = 0.0
       modwdht1 = 0.0
       modwdht3 = 0.0
       modwdhu3ht1 = 0.0
       modwdhu3ht3 = 0.0
       modwdhu33ht1 = 0.0
       modwdhu33ht3 = 0.0
       do k=1,nbasis
        wdhp1p0   = wdhp1p0   + ave_wdh(ib,k)*ave_hp1p0(is,k,jb)
        wdhr11p0  = wdhr11p0  + ave_wdh(ib,k)*ave_hr11p0(is,k,jb)
        wdhr13p0  = wdhr13p0  + ave_wdh(ib,k)*ave_hr13p0(is,k,jb)
        wdht1  = wdht1  + ave_wdh(ib,k)*ave_ht1(is,k,jb)
        wdht3  = wdht3  + ave_wdh(ib,k)*ave_ht3(is,k,jb)
        wdhu1   = wdhu1   + ave_wdh(ib,k)*ave_hu1(is,k,jb)
        wdhu3   = wdhu3   + ave_wdh(ib,k)*ave_hu3(is,k,jb)
        wdhu3ht1   = wdhu3ht1 &
          + ave_wdh(ib,k)*ave_hu3ht1(is,k,jb)
        wdhu3ht3   = wdhu3ht3 &
          + ave_wdh(ib,k)*ave_hu3ht3(is,k,jb)
        wdhu33  = wdhu33  + ave_wdh(ib,k)*ave_hu33(is,k,jb)
        wdhu33ht1   = wdhu33ht1 &
          + ave_wdh(ib,k)*ave_hu33ht1(is,k,jb)
        wdhu33ht3   = wdhu33ht3 &
          + ave_wdh(ib,k)*ave_hu33ht3(is,k,jb)
        modwdht1  = modwdht1  + ave_modwdh(ib,k)*ave_ht1(is,k,jb)
        modwdht3  = modwdht3  + ave_modwdh(ib,k)*ave_ht3(is,k,jb)
        modwdhu1  = modwdhu1  + ave_modwdh(ib,k)*ave_hu1(is,k,jb)
        modwdhu3  = modwdhu3  + ave_modwdh(ib,k)*ave_hu3(is,k,jb)
        modwdhu3ht1  = modwdhu3ht1 &
          + ave_modwdh(ib,k)*ave_hu3ht1(is,k,jb)
        modwdhu3ht3  = modwdhu3ht3 & 
          + ave_modwdh(ib,k)*ave_hu3ht3(is,k,jb)
        modwdhu33  = modwdhu33 & 
          + ave_modwdh(ib,k)*ave_hu33(is,k,jb)
        modwdhu33ht1  = modwdhu33ht1 & 
          + ave_modwdh(ib,k)*ave_hu33ht1(is,k,jb)
        modwdhu33ht3  = modwdhu33ht3 & 
          + ave_modwdh(ib,k)*ave_hu33ht3(is,k,jb)
       enddo
       if(use_bper_in)then
         do k=1,nbasis
           wdhp1b0   = wdhp1b0   + ave_wdh(ib,k)*ave_hp1b0(is,k,jb)
           wdhr11b0  = wdhr11b0  + ave_wdh(ib,k)*ave_hr11b0(is,k,jb)
           wdhr13b0  = wdhr13b0  + ave_wdh(ib,k)*ave_hr13b0(is,k,jb)
         enddo
         if(vpar_model_in.eq.0)then
           do k=1,nbasis
             wdhp1bp   = wdhp1bp   + ave_wdh(ib,k)*ave_hp1bp(is,k,jb)
             wdhr11bp  = wdhr11bp  + ave_wdh(ib,k)*ave_hr11bp(is,k,jb)
             wdhr13bp  = wdhr13bp  + ave_wdh(ib,k)*ave_hr13bp(is,k,jb)
           enddo
         endif
       endif
       ave_wdhp1p0(is,ib,jb) = wdhp1p0
       ave_wdhr11p0(is,ib,jb) = wdhr11p0
       ave_wdhr13p0(is,ib,jb) = wdhr13p0
       ave_wdhp1b0(is,ib,jb) = wdhp1b0
       ave_wdhr11b0(is,ib,jb) = wdhr11b0
       ave_wdhr13b0(is,ib,jb) = wdhr13bp
       ave_wdhp1bp(is,ib,jb) = wdhp1bp
       ave_wdhr11bp(is,ib,jb) = wdhr11bp
       ave_wdhr13bp(is,ib,jb) = wdhr13bp
       ave_wdht1(is,ib,jb)= wdht1
       ave_wdht3(is,ib,jb)= wdht3
       ave_wdhu1(is,ib,jb) = wdhu1
       ave_wdhu3(is,ib,jb) = wdhu3
       ave_wdhu3ht1(is,ib,jb) = wdhu3ht1
       ave_wdhu3ht3(is,ib,jb) = wdhu3ht3
       ave_wdhu33(is,ib,jb) = wdhu33
       ave_wdhu33ht1(is,ib,jb) = wdhu33ht1
       ave_wdhu33ht3(is,ib,jb) = wdhu33ht3
       ave_modwdht1(is,ib,jb)= modwdht1
       ave_modwdht3(is,ib,jb)= modwdht3
       ave_modwdhu1(is,ib,jb)= modwdhu1
       ave_modwdhu3(is,ib,jb)= modwdhu3
       ave_modwdhu3ht1(is,ib,jb)= modwdhu3ht1
       ave_modwdhu3ht3(is,ib,jb)= modwdhu3ht3
       ave_modwdhu33(is,ib,jb)= modwdhu33
       ave_modwdhu33ht1(is,ib,jb)= modwdhu33ht1
       ave_modwdhu33ht3(is,ib,jb)= modwdhu33ht3
      enddo
      enddo
      enddo
!
      END SUBROUTINE wd_h
!
      SUBROUTINE wd_g
!***************************************************************
!
!   compute the products ave_wd*ave_g
!
!***************************************************************
      USE tglf_dimensions
      USE tglf_coeff
      USE tglf_global
!
      IMPLICIT NONE
      INTEGER :: is,ib,jb,k
!
      do is=ns0,ns
      do ib=1,nbasis
      do jb=1,nbasis
       wdgp1p0 = 0.0
       wdgr11p0 = 0.0
       wdgr13p0 = 0.0
       wdgp1b0 = 0.0
       wdgr11b0 = 0.0
       wdgr13b0 = 0.0
       wdgp1bp = 0.0
       wdgr11bp = 0.0
       wdgr13bp = 0.0
       wdgu1 = 0.0
       wdgu3 = 0.0
       wdgu33 = 0.0
       wdgt1 = 0.0
       wdgt3 = 0.0
       wdgu3gt1 = 0.0
       wdgu3gt3 = 0.0
       wdgu33gt1 = 0.0
       wdgu33gt3 = 0.0
       modwdgu1 = 0.0
       modwdgu3 = 0.0
       modwdgu33 = 0.0
       modwdgt1 = 0.0
       modwdgt3 = 0.0
       modwdgu3gt1 = 0.0
       modwdgu3gt3 = 0.0
       modwdgu33gt1 = 0.0
       modwdgu33gt3 = 0.0
       do k=1,nbasis
         wdgp1p0   = wdgp1p0   + ave_wdg(ib,k)*ave_gp1p0(is,k,jb)
         wdgr11p0  = wdgr11p0  + ave_wdg(ib,k)*ave_gr11p0(is,k,jb)
         wdgr13p0  = wdgr13p0  + ave_wdg(ib,k)*ave_gr13p0(is,k,jb)
         wdgu1   = wdgu1   + ave_wdg(ib,k)*ave_gu1(is,k,jb)
         wdgu3   = wdgu3   + ave_wdg(ib,k)*ave_gu3(is,k,jb)
         wdgu33  = wdgu33  + ave_wdg(ib,k)*ave_gu33(is,k,jb)
         wdgt1  = wdgt1  + ave_wdg(ib,k)*ave_gt1(is,k,jb)
         wdgt3  = wdgt3  + ave_wdg(ib,k)*ave_gt3(is,k,jb)
         wdgu3gt1   = wdgu3gt1 &  
          + ave_wdg(ib,k)*ave_gu3gt1(is,k,jb)
         wdgu33gt1  = wdgu33gt1 &
          + ave_wdg(ib,k)*ave_gu33gt1(is,k,jb)
         wdgu3gt3   = wdgu3gt3 & 
          + ave_wdg(ib,k)*ave_gu3gt3(is,k,jb)
         wdgu33gt3  = wdgu33gt3 &  
          + ave_wdg(ib,k)*ave_gu33gt3(is,k,jb)
         modwdgu1  = modwdgu1 & 
          + ave_modwdg(ib,k)*ave_gu1(is,k,jb)
         modwdgu3  = modwdgu3 &  
          + ave_modwdg(ib,k)*ave_gu3(is,k,jb)
         modwdgu33  = modwdgu33 & 
          + ave_modwdg(ib,k)*ave_gu33(is,k,jb)
         modwdgt1  = modwdgt1 & 
          + ave_modwdg(ib,k)*ave_gt1(is,k,jb)
         modwdgt3  = modwdgt3 & 
          + ave_modwdg(ib,k)*ave_gt3(is,k,jb)
         modwdgu3gt1  = modwdgu3gt1 & 
          + ave_modwdg(ib,k)*ave_gu3gt1(is,k,jb)
         modwdgu33gt1  = modwdgu33gt1 &  
          + ave_modwdg(ib,k)*ave_gu33gt1(is,k,jb)
         modwdgu3gt3  = modwdgu3gt3 &
          + ave_modwdg(ib,k)*ave_gu3gt3(is,k,jb)
         modwdgu33gt3  = modwdgu33gt3 &  
          + ave_modwdg(ib,k)*ave_gu33gt3(is,k,jb)
       enddo
       if(use_bper_in)then
         do k=1,nbasis
           wdgp1b0   = wdgp1b0   + ave_wdg(ib,k)*ave_gp1b0(is,k,jb)
           wdgr11b0  = wdgr11b0  + ave_wdg(ib,k)*ave_gr11b0(is,k,jb)
           wdgr13b0  = wdgr13b0  + ave_wdg(ib,k)*ave_gr13b0(is,k,jb)
         enddo
         if(vpar_model_in.eq.0)then
           do k=1,nbasis
             wdgp1bp   = wdgp1bp   + ave_wdg(ib,k)*ave_gp1bp(is,k,jb)
             wdgr11bp  = wdgr11bp  + ave_wdg(ib,k)*ave_gr11bp(is,k,jb)
             wdgr13bp  = wdgr13bp  + ave_wdg(ib,k)*ave_gr13bp(is,k,jb)
           enddo
         endif
       endif
       ave_wdgp1p0(is,ib,jb) = wdgp1p0
       ave_wdgr11p0(is,ib,jb) = wdgr11p0
       ave_wdgr13p0(is,ib,jb) = wdgr13p0
       ave_wdgp1b0(is,ib,jb) = wdgp1b0
       ave_wdgr11b0(is,ib,jb) = wdgr11b0
       ave_wdgr13b0(is,ib,jb) = wdgr13b0
       ave_wdgp1bp(is,ib,jb) = wdgp1bp
       ave_wdgr11bp(is,ib,jb) = wdgr11bp
       ave_wdgr13bp(is,ib,jb) = wdgr13bp
       ave_wdgu1(is,ib,jb) = wdgu1
       ave_wdgu3(is,ib,jb) = wdgu3
       ave_wdgu33(is,ib,jb) = wdgu33
       ave_wdgt1(is,ib,jb)= wdgt1
       ave_wdgt3(is,ib,jb)= wdgt3
       ave_wdgu3gt1(is,ib,jb)= wdgu3gt1
       ave_wdgu3gt3(is,ib,jb)= wdgu3gt3
       ave_wdgu33gt1(is,ib,jb)= wdgu33gt1
       ave_wdgu33gt3(is,ib,jb)= wdgu33gt3
       ave_modwdgu1(is,ib,jb)= modwdgu1
       ave_modwdgu3(is,ib,jb)= modwdgu3
       ave_modwdgu33(is,ib,jb)= modwdgu33
       ave_modwdgt1(is,ib,jb)= modwdgt1
       ave_modwdgt3(is,ib,jb)= modwdgt3
       ave_modwdgu3gt1(is,ib,jb)= modwdgu3gt1
       ave_modwdgu3gt3(is,ib,jb)= modwdgu3gt3
       ave_modwdgu33gt1(is,ib,jb)= modwdgu33gt1
       ave_modwdgu33gt3(is,ib,jb)= modwdgu33gt3
      enddo
      enddo
      enddo
!
      END SUBROUTINE wd_g
!
      SUBROUTINE kpar_h
!***************************************************************
!
!   compute the products ave_kpar*ave_h and ave_modkpar*ave_h
!
!***************************************************************
      USE tglf_dimensions
      USE tglf_coeff
      USE tglf_global
!
      IMPLICIT NONE
      INTEGER :: is,ib,jb,k
      COMPLEX :: kp_hnp0,kp_hp1p0,kp_hp3p0
      COMPLEX :: kp_hp1b0,kp_hr11b0,kp_hr13b0
      COMPLEX :: kp_hnbp,kp_hp3bp
      COMPLEX :: kp_hp1bp,kp_hr11bp,kp_hr13bp
      COMPLEX :: kp_hu1,kp_hu3
      COMPLEX :: kp_ht1,kp_ht3
!
      do is=ns0,ns
       do ib=1,nbasis
       do jb=1,nbasis
        kp_hnp0 = 0.0
        kp_hp1p0 = 0.0
        kp_hp3p0 = 0.0
        kp_hp1b0 = 0.0
        kp_hr11b0 = 0.0
        kp_hr13b0 = 0.0
        kp_hnbp = 0.0
        kp_hp3bp = 0.0
        kp_hp1bp = 0.0
        kp_hr11bp = 0.0
        kp_hr13bp = 0.0
        kp_hu1 = 0.0
        kp_hu3 = 0.0
        kp_ht1 = 0.0
        kp_ht3 = 0.0
        modkpar_hu1 = 0.0
        modkpar_hu3 = 0.0
        do k=1,nbasis
         kp_hnp0        = kp_hnp0 &
            +ave_hnp0(is,ib,k)*ave_kpar(k,jb)
         kp_hp1p0        = kp_hp1p0 &
            +ave_hp1p0(is,ib,k)*ave_kpar(k,jb)
         kp_hp3p0        = kp_hp3p0 &
            +ave_hp3p0(is,ib,k)*ave_kpar(k,jb)
         kp_hu1     = kp_hu1  &
            +ave_hu1(is,ib,k)*ave_kpar_eff(is,k,jb)
         kp_hu3     = kp_hu3 &
            +ave_hu3(is,ib,k)*ave_kpar_eff(is,k,jb)
         kp_ht1     = kp_ht1 &
            +ave_kpar_eff(is,ib,k)*ave_ht1(is,k,jb)
         kp_ht3     = kp_ht3 &
            +ave_kpar_eff(is,ib,k)*ave_ht3(is,k,jb)
         modkpar_hu1        = modkpar_hu1 &
            +ave_modkpar_eff(is,ib,k)*ave_hu1(is,k,jb)
         modkpar_hu3        = modkpar_hu3 &
            +ave_modkpar_eff(is,ib,k)*ave_hu3(is,k,jb)
        enddo
        if(use_bper_in)then
          do k=1,nbasis
            kp_hp1b0        = kp_hp1b0 &
            +ave_hp1b0(is,ib,k)*ave_kpar(k,jb)
            kp_hr11b0       = kp_hr11b0 &
            +ave_hr11b0(is,ib,k)*ave_kpar(k,jb)
            kp_hr13b0       = kp_hr13b0 &
            +ave_hr13b0(is,ib,k)*ave_kpar(k,jb)
          enddo
          if(vpar_model_in.eq.0)then
            do k=1,nbasis
             kp_hnbp        = kp_hnbp &
             +ave_hnbp(is,ib,k)*ave_kpar(k,jb)
             kp_hp3bp        = kp_hp3bp &
             +ave_hp3bp(is,ib,k)*ave_kpar(k,jb)
             kp_hp1bp        = kp_hp1bp &
             +ave_hp1bp(is,ib,k)*ave_kpar(k,jb)
             kp_hr11bp       = kp_hr11bp &
             +ave_hr11bp(is,ib,k)*ave_kpar(k,jb)
             kp_hr13bp       = kp_hr13bp &
             +ave_hr13bp(is,ib,k)*ave_kpar(k,jb)
            enddo
          endif
        endif
        ave_kparhnp0(is,ib,jb) = kp_hnp0
        ave_kparhp1p0(is,ib,jb) = kp_hp1p0
        ave_kparhp3p0(is,ib,jb) = kp_hp3p0
        ave_kparhp1b0(is,ib,jb) = kp_hp1b0
        ave_kparhr11b0(is,ib,jb) = kp_hr11b0
        ave_kparhr13b0(is,ib,jb) = kp_hr13b0
        ave_kparhnbp(is,ib,jb) = kp_hnbp
        ave_kparhp3bp(is,ib,jb) = kp_hp3bp
        ave_kparhp1bp(is,ib,jb) = kp_hp1bp
        ave_kparhr11bp(is,ib,jb) = kp_hr11bp
        ave_kparhr13bp(is,ib,jb) = kp_hr13b0
        ave_kparhu1(is,ib,jb) = kp_hu1
        ave_kparhu3(is,ib,jb) = kp_hu3
        ave_kparht1(is,ib,jb) = kp_ht1
        ave_kparht3(is,ib,jb) = kp_ht3
        ave_modkparhu1(is,ib,jb) = modkpar_hu1
        ave_modkparhu3(is,ib,jb) = modkpar_hu3
       enddo
       enddo
      enddo
!
      END SUBROUTINE kpar_h
!
      SUBROUTINE kpar_g
!***************************************************************
!
!   compute the products ave_kpar*ave_h and ave_modkpar*ave_h
!
!***************************************************************
      USE tglf_dimensions
      USE tglf_coeff
      USE tglf_global
!
      IMPLICIT NONE
      INTEGER :: is,ib,jb,k
      COMPLEX :: kp_gnp0,kp_gp1p0,kp_gp3p0
      COMPLEX :: kp_gp1b0,kp_gr11b0,kp_gr13b0
      COMPLEX :: kp_gnbp,kp_gp3bp
      COMPLEX :: kp_gp1bp,kp_gr11bp,kp_gr13bp
      COMPLEX :: kp_gu1,kp_gu3
      COMPLEX :: kp_gt1,kp_gt3
!
      do is=ns0,ns
       do ib=1,nbasis
       do jb=1,nbasis
        kp_gnp0 = 0.0
        kp_gp1p0 = 0.0
        kp_gp3p0 = 0.0
        kp_gp1b0 = 0.0
        kp_gr11b0 = 0.0
        kp_gr13b0 = 0.0
        kp_gnbp = 0.0
        kp_gp3bp = 0.0
        kp_gp1bp = 0.0
        kp_gr11bp = 0.0
        kp_gr13bp = 0.0
        kp_gu1 = 0.0
        kp_gu3 = 0.0
        kp_gt1 = 0.0
        kp_gt3 = 0.0
        modkpar_gu1 = 0.0
        modkpar_gu3 = 0.0
        do k=1,nbasis
         kp_gnp0        = kp_gnp0 & 
            +ave_gnp0(is,ib,k)*ave_kpar(k,jb)
         kp_gp1p0        = kp_gp1p0 & 
            +ave_gp1p0(is,ib,k)*ave_kpar(k,jb)
         kp_gp3p0        = kp_gp3p0 & 
            +ave_gp3p0(is,ib,k)*ave_kpar(k,jb)
         kp_gu1        = kp_gu1 & 
            +ave_gu1(is,ib,k)*ave_kpar_eff(is,k,jb)
         kp_gu3        = kp_gu3 & 
            +ave_gu3(is,ib,k)*ave_kpar_eff(is,k,jb)
         kp_gt1        = kp_gt1 &
            +ave_kpar_eff(is,ib,k)*ave_gt1(is,k,jb)
         kp_gt3        = kp_gt3 &
            +ave_kpar_eff(is,ib,k)*ave_gt3(is,k,jb)
         modkpar_gu1  = modkpar_gu1 &
         +ave_modkpar_eff(is,ib,k)*ave_gu1(is,k,jb)
         modkpar_gu3  = modkpar_gu3 &
         +ave_modkpar_eff(is,ib,k)*ave_gu3(is,k,jb)
        enddo
        if(use_bper_in)then
          do k=1,nbasis
            kp_gp1b0        = kp_gp1b0 & 
            +ave_gp1b0(is,ib,k)*ave_kpar(k,jb)
            kp_gr11b0       = kp_gr11b0 &
            +ave_gr11b0(is,ib,k)*ave_kpar(k,jb)
            kp_gr13b0       = kp_gr13b0 &
            +ave_gr13b0(is,ib,k)*ave_kpar(k,jb)
          enddo
          if(vpar_model_in.eq.0)then
           do k=1,nbasis
            kp_gnbp        = kp_gnbp & 
            +ave_gnbp(is,ib,k)*ave_kpar(k,jb)
            kp_gp3bp        = kp_gp3bp & 
            +ave_gp3bp(is,ib,k)*ave_kpar(k,jb)
            kp_gp1bp        = kp_gp1bp & 
            +ave_gp1bp(is,ib,k)*ave_kpar(k,jb)
            kp_gr11bp       = kp_gr11bp &
            +ave_gr11bp(is,ib,k)*ave_kpar(k,jb)
            kp_gr13bp       = kp_gr13bp &
            +ave_gr13bp(is,ib,k)*ave_kpar(k,jb)
           enddo
          endif
        endif
        ave_kpargnp0(is,ib,jb) = kp_gnp0
        ave_kpargp1p0(is,ib,jb) = kp_gp1p0
        ave_kpargp3p0(is,ib,jb) = kp_gp3p0
        ave_kpargp1b0(is,ib,jb) = kp_gp1b0
        ave_kpargr11b0(is,ib,jb) = kp_gr11b0
        ave_kpargr13b0(is,ib,jb) = kp_gr13b0
        ave_kpargnbp(is,ib,jb) = kp_gnbp
        ave_kpargp3bp(is,ib,jb) = kp_gp3bp
        ave_kpargp1bp(is,ib,jb) = kp_gp1bp
        ave_kpargr11bp(is,ib,jb) = kp_gr11bp
        ave_kpargr13bp(is,ib,jb) = kp_gr13bp
        ave_kpargu1(is,ib,jb) = kp_gu1
        ave_kpargu3(is,ib,jb) = kp_gu3
        ave_kpargt1(is,ib,jb) = kp_gt1
        ave_kpargt3(is,ib,jb) = kp_gt3
        ave_modkpargu1(is,ib,jb) = modkpar_gu1
        ave_modkpargu3(is,ib,jb) = modkpar_gu3
       enddo
       enddo
      enddo
!
      END SUBROUTINE kpar_g
!
      SUBROUTINE modwd
!***************************************************************
!
!   compute the matricies modwdh and modwdg
!
!***************************************************************
      USE tglf_dimensions
      USE tglf_global
      USE tglf_coeff
!
      IMPLICIT NONE
      INTEGER :: nm,i,j,k
      INTEGER :: lwork,info
      REAL :: a(nbasis,nbasis)
      REAL :: w(nbasis)
      REAL :: work(34*nbasis)
!
      nm = nbasis
!
      if(nm.eq.1)then
        ave_modwdh(1,1) = ABS(ave_wdh(1,1))
        ave_modwdg(1,1) = ABS(ave_wdg(1,1))
        RETURN
      endif
!
!  find the eigenvalues of ave_wdh
!
       do i=1,nm
       do j=i,nm
         a(i,j) = ave_wdh(i,j)
       enddo
       enddo
       lwork=34*nm
! call LAPACK routine for symmetric real eigenvalue problem DSYEV
       call DSYEV('V','U',nm,a,nm,w,work,lwork,info)
       if(info.ne.0)CALL tglf_error(1,"DSYEV failed in modwd")
! debug
!       write(*,*)"ave_wdh eigenvalues"
!       do i=1,nm
!       write(*,*)i,w(i)
!       enddo
!       write(*,*)"ave_wdh eigenvectors"
!       do i=1,nm
!       write(*,*)"******",i
!       do j=1,nm
!       write(*,*)j,a(j,i)
!       enddo
!       enddo
! 
! regularize the eigenvalues
       do k=1,nm
           if(ABS(w(k)).lt.wd_zero_in)then
             if(w(k).ge.0.0)then
               w(k)=wd_zero_in
             else
               w(k)=-wd_zero_in
             endif
           endif
       enddo
! compute ave_modwd and recompute ave_wd with regularized eigenvalues
! note that the DSYEV normalized eigenvectors are now in a(i,j)
       do i=1,nm
       do j=1,nm
         ave_modwdh(i,j) = 0.0
         ave_wdh(i,j) = 0.0
         do k=1,nm
           ave_modwdh(i,j) = ave_modwdh(i,j) + ABS(w(k))*a(i,k)*a(j,k)
           ave_wdh(i,j) = ave_wdh(i,j) + w(k)*a(i,k)*a(j,k)
         enddo
       enddo
       enddo
!
! debug check modwdh eigenvalues
!       do i=1,nm
!       do j=i,nm
!         a(i,j) = ave_modwdh(i,j)
!       enddo
!       enddo
!       call DSYEV('V','U',nm,a,nb,w,work,lwork,info)
!       write(*,*)"modwdh eigenvalues"
!       do i=1,nm
!       write(*,*)i,w(i)
!       enddo
!
!
!  find the eigenvalues of ave_wdh
!
       do i=1,nm
       do j=i,nm
         a(i,j) = ave_wdg(i,j)
       enddo
       enddo
       lwork=34*nm
! call LAPACK routine for symmetric real eigenvalue problem DSYEV
       call DSYEV('V','U',nm,a,nm,w,work,lwork,info)
       if(info.ne.0)CALL tglf_error(1,"DSYEV failed in modwd")
! debug
!       write(*,*)"ave_wd eigenvalues"
!       do i=1,nm
!       write(*,*)i,w(i)
!       enddo
!       write(*,*)"ave_wd eigenvectors"
!       do i=1,nm
!       write(*,*)"******",i
!       do j=1,nm
!       write(*,*)j,a(j,i)
!       enddo
!       enddo
! 
! regularize the eigenvalues
       do k=1,nm
           if(ABS(w(k)).lt.wd_zero_in)then
             if(w(k).ge.0.0)then
               w(k)=wd_zero_in
             else
               w(k)=-wd_zero_in
             endif
           endif
       enddo
! compute ave_modwd and recompute ave_wd with regularized eigenvalues
! note that the DSYEV normalized eigenvectors are now in a(i,j)
       do i=1,nm
       do j=1,nm
         ave_modwdg(i,j) = 0.0
         ave_wdg(i,j) = 0.0
         do k=1,nm
           ave_modwdg(i,j) = ave_modwdg(i,j) + ABS(w(k))*a(i,k)*a(j,k)
           ave_wdg(i,j) = ave_wdg(i,j) + w(k)*a(i,k)*a(j,k)
         enddo
       enddo
       enddo
!
! debug check modwdg eigenvalues
!       do i=1,nm
!       do j=i,nm
!         a(i,j) = ave_modwdg(i,j)
!       enddo
!       enddo
!       call DSYEV('V','U',nm,a,nb,w,work,lwork,info)
!       write(*,*)"modwdg eigenvalues"
!       do i=1,nm
!       write(*,*)i,w(i)
!       enddo
!
!
      END SUBROUTINE modwd
!
      SUBROUTINE modkpar
!***************************************************************
!
!   compute the matrix mod_kpar
!
!***************************************************************
      USE tglf_dimensions
      USE tglf_global
      USE tglf_coeff
!
      IMPLICIT NONE
      INTEGER :: nm,i,j,k,is
      INTEGER :: lwork,info
      REAL :: w(nbasis)
      REAL :: rwork(3*nbasis-2)
      COMPLEX :: work(34*nbasis)
      COMPLEX :: a(nbasis,nbasis),b(nbasis,nbasis)
!      COMPLEX :: xi
!
      nm = nbasis
      if(nm.eq.1)then
        do is=ns0,ns
         ave_modkpar_eff(is,1,1) = xi*ave_kpar_eff(is,1,1)
        enddo
      else
!
! find the eigenvalues and eigenvectors
!
      if(vpar_model_in.ne.1)then     
        do i=1,nm
        do j=i,nm
          a(i,j) = xi*ave_kpar(i,j)
        enddo
        enddo
        lwork=34*nm
        call ZHEEV('V','U',nm,a,nm,w,work,lwork,rwork,info)
        if(info.ne.0)CALL tglf_error(1,"ZHEEV failed in modkpar")
!       write(*,*)"kpar eigenvalues"
!       do i=1,nm
!       write(*,*)i,w(i)
!       enddo
!       write(*,*)"kpar eigenvectors"
!       do i=1,nm
!       write(*,*)"******",i
!       do j=1,nm
!       write(*,*)j,a(i,j)
!       enddo
!       enddo
! construct mod_kpar
! note modkpar is a real symmetric matrix
!       write(*,*)"modkpar"
         do i=1,nm
         do j=1,nm
           b(i,j) = 0.0
           do k=1,nm
             b(i,j) = b(i,j)+ ABS(w(k))*a(i,k)*CONJG(a(j,k))
           enddo
!         write(*,*)i,j,b(i,j)
           ave_modkpar(i,j)=REAl(b(i,j))
           do is=ns0,ns
             ave_modkpar_eff(is,i,j) = b(i,j)
           enddo
         enddo
         enddo
       else ! vpar_model_in = 1
!
         do is=ns0,ns
           do i=1,nm
           do j=1,nm
             a(i,j) = xi*ave_kpar_eff(is,i,j)
!            write(*,*)"a",is,i,j,a(i,j)
           enddo
           enddo
           lwork=34*nm
           call ZHEEV('V','U',nm,a,nm,w,work,lwork,rwork,info)
           if(info.ne.0)CALL tglf_error(1,"ZHEEV failed in modkpar_eff")
           do i=1,nm
!            write(*,*)i,"w = ",w(i)
           do j=1,nm
             b(i,j) = 0.0
             do k=1,nm
               b(i,j) = b(i,j)+ ABS(w(k))*a(i,k)*CONJG(a(j,k))
             enddo
!         write(*,*)i,j,b(i,j)
             ave_modkpar_eff(is,i,j) = b(i,j)
           enddo
           enddo
         enddo 
       endif
! debug check modk eigenvalues
!      do i=1,nm
!      do j=i,nm
!        a(i,j) = ave_modkpar(i,j)
!      enddo
!      enddo
!      call ZHEEV('V','U',nm,a,nb,w,work,lwork,rwork,info)
!      if(info.ne.0)write(*,*)"modkpar info =",info
!      write(*,*)"kpar eigenvalues"
!      do i=1,nm
!        write(*,*)i,w(i)
!      enddo
!
      endif
!
      END SUBROUTINE modkpar
!
      SUBROUTINE ave_inv(ave_m,ave_minv)
!***************************************************************
!
!   compute the inverse matrix ave_minv
!   of the symmetric real matrix ave_m
!
!***************************************************************
!
      USE tglf_dimensions
! 
      IMPLICIT NONE
      INTEGER :: nm,is,i,j,k
      INTEGER :: lwork,info
      REAL,INTENT(IN),DIMENSION(ns,nbasis_max,nbasis_max) :: ave_m
      REAL,INTENT(OUT),DIMENSION(ns,nbasis_max,nbasis_max) :: ave_minv
      REAL :: detm,zero
!      REAL :: check
      REAL,DIMENSION(nbasis,nbasis) :: a
      REAL,DIMENSION(nbasis) :: w
      REAL,DIMENSION(34*nbasis) :: work

      nm = nbasis
!
      if(nm.eq.1)then
        do is=ns0,ns
          ave_minv(is,1,1) = 1.0/ave_m(is,1,1)
        enddo
        go to 100
      endif
!
      if(nm.eq.2)then
        do is=ns0,ns
         detm = ave_m(is,1,1)*ave_m(is,2,2)-ave_m(is,1,2)*ave_m(is,2,1)
         if(detm.eq.0.0)detm = 1.E-12
         ave_minv(is,1,1) = ave_m(is,2,2)/detm
         ave_minv(is,2,2) = ave_m(is,1,1)/detm
         ave_minv(is,1,2) = -ave_m(is,1,2)/detm
         ave_minv(is,2,1) = -ave_m(is,2,1)/detm
        enddo
        go to 100
      endif
!
!  for nm > 2 need to use eigensytem solver
!
!  find the eigenvalues of ave_m
!
      zero = 1.0E-10
      do is = ns0,ns
       do i=1,nm
       do j=i,nm
         a(i,j) = ave_m(is,i,j)
       enddo
       enddo
       lwork=34*nm
! call LAPACK routine for symmetric real eigenvalue problem DSYEV
       call DSYEV('V','U',nm,a,nm,w,work,lwork,info)
       if(info.ne.0)CALL tglf_error(1,"DSYEV failed in ave_inv")
!
! debug
!       write(*,*)"ave_m eigenvalues"
!       do i=1,nm
!       write(*,*)i,rr(i),ri(i)
!       enddo
!       write(*,*)"ave_m eigenvectors"
!       do i=1,nm
!       write(*,*)"******",i
!       do j=1,nm
!       write(*,*)j,vr(j,i),vi(j,i)
!       enddo
!       enddo
!
! regularize the eigenvalues
       do k=1,nm
           if(ABS(w(k)).lt.zero)then
             if(w(k).ge.0.0)then
               w(k)=zero
             else
               w(k)=-zero
             endif
           endif
       enddo
! compute ave_inv
       do i=1,nm
       do j=1,nm
         ave_minv(is,i,j) = 0.0
         do k=1,nm
           ave_minv(is,i,j) = ave_minv(is,i,j) + a(i,k)*a(j,k)/w(k)
         enddo
       enddo
       enddo
      enddo ! is loop
 100  continue
! debug check minv
!      do is=ns0,ns
!       write(*,*)"check minv",is
!        do i=1,nm
!        do j=1,nm
!         check = 0.0
!         do k=1,nm
!           check = check + ave_m(is,i,k)*ave_minv(is,k,j)
!         enddo
!         write(*,*)i,j,check
!        enddo
!        enddo
!      enddo
!
      END SUBROUTINE ave_inv
!
      SUBROUTINE ave_inv0(ave_m,ave_minv)
!***************************************************************
!
!   compute the inverse matrix ave_minv
!   of the symmetric real matrix ave_m
!
!***************************************************************
      USE tglf_dimensions
!
      IMPLICIT NONE
      INTEGER :: nm,i,j,k
      INTEGER :: lwork,info
      REAL :: detm,zero
!      REAL :: check
      REAL :: a(nbasis,nbasis)
      REAL :: w(nbasis)
      REAL :: work(34*nbasis)
      REAL,INTENT(IN),DIMENSION(nbasis_max,nbasis_max) :: ave_m
      REAL,INTENT(OUT),DIMENSION(nbasis_max,nbasis_max) :: ave_minv
!
      nm = nbasis
      zero = 1.0E-12
!
      if(nm.eq.1)then
        ave_minv(1,1) = 1.0/ave_m(1,1)
        go to 100
      endif
!
      if(nm.eq.2)then
        detm = ave_m(1,1)*ave_m(2,2)-ave_m(1,2)*ave_m(2,1)
        if(detm.eq.0.0)detm = zero
        ave_minv(1,1) = ave_m(2,2)/detm
        ave_minv(2,2) = ave_m(1,1)/detm
        ave_minv(1,2) = -ave_m(1,2)/detm
        ave_minv(2,1) = -ave_m(2,1)/detm
        go to 100
      endif
!
!  for nm > 2 need to use eigensytem solver
!
!
!  find the eigenvalues of ave_m
!
       do i=1,nm
       do j=i,nm
         a(i,j) = ave_m(i,j)
       enddo
       enddo
       lwork=34*nm
! call LAPACK routine for symmetric real eigenvalue problem DSYEV
       call DSYEV('V','U',nm,a,nm,w,work,lwork,info)
       if(info.ne.0)CALL tglf_error(1,"DSYEV failed in ave_inv0")
!
! debug
!       write(*,*)"ave_m eigenvalues"
!       do i=1,nm
!       write(*,*)i,rr(i),ri(i)
!       enddo
!       write(*,*)"ave_m eigenvectors"
!       do i=1,nm
!       write(*,*)"******",i
!       do j=1,nm
!       write(*,*)j,vr(j,i),vi(j,i)
!       enddo
!       enddo
!
! regularize the eigenvalues
       do k=1,nm
           if(ABS(w(k)).lt.zero)then
             if(w(k).ge.0.0)then
               w(k)=zero
             else
               w(k)=-zero
             endif
           endif
       enddo
! compute ave_inv
       do i=1,nm
       do j=1,nm
         ave_minv(i,j) = 0.0
         do k=1,nm
           ave_minv(i,j) = ave_minv(i,j) + a(i,k)*a(j,k)/w(k)
         enddo
       enddo
       enddo
 100   continue
! debug check minv
!       write(*,*)"check minv"
!        do i=1,nm
!        do j=1,nm
!         check = 0.0
!         do k=1,nm
!           check = check + ave_m(i,k)*ave_minv(k,j)
!         enddo
!         write(*,*)i,j,check
!        enddo
!        enddo
!
      END SUBROUTINE ave_inv0
!
      SUBROUTINE FLR_xgrid
!
!*************************************************************
!  begin computation of the FLR integrals
!*************************************************************
!  compute the FLR integrals at the hermite nodes
!
      USE tglf_dimensions
      USE tglf_global
      USE tglf_species
      USE tglf_xgrid
!
      IMPLICIT NONE
      INTEGER :: i,is
!      
      REAL :: b,bb,ftx,fth,ft2
!      REAL :: FLR_Hn,FLR_dHp1,FLR_dHp3
!      REAL :: FLR_dHr11,FLR_dHr13,FLR_dHr33
!
! debug
!        b = ky_in
!        ft = rmin_sa
!        write(*,*)"b = ",b,"ft = ",ft
!        write(*,*)"hn = ",FLR_Hn(ft,b)
!        write(*,*)"dhp1 = ",FLR_dHp1(ft,b)
!        write(*,*)"dhp3 = ",FLR_dHp3(ft,b)
!        write(*,*)"dhr11 = ",FLR_dHr11(ft,b)
!        write(*,*)"dhr13 = ",FLR_dHr13(ft,b)
!        write(*,*)"dhr33 = ",FLR_dHr33(ft,b)
!        write(*,*)"dhw113= ",FLR_dHw113(ft,b)
!        write(*,*)"dhw133= ",FLR_dHw133(ft,b)
!        write(*,*)"dhw333= ",FLR_dHw333(ft,b)
!        write(*,*)"FLR ft=",ft,ft_min

! debug
!
      fth=1.0
      do i=1,nx
      do is =ns0,ns
        bb=taus(is)*mass(is)*(ky/zs(is))**2 
        b=bb*b0x(i)/b2x(i)
        hxn(is,i)    = FLR_Hn(fth,b)
        hxp1(is,i)   = hxn(is,i)
        hxp3(is,i)   = FLR_dHp3(fth,b)+hxn(is,i)
        hxr11(is,i)  = 3.0*hxp1(is,i)
        hxr13(is,i)  = FLR_dHr13(fth,b)+(5.0/3.0)*hxp1(is,i)
        hxr33(is,i)  = FLR_dHr33(fth,b)+(5.0/3.0)*hxp3(is,i)
        hxw113(is,i) = FLR_dHw113(fth,b)+(7.0/3.0)*hxr11(is,i)
        hxw133(is,i) = FLR_dHw133(fth,b)+(7.0/3.0)*hxr13(is,i)
        hxw333(is,i) = FLR_dHw333(fth,b)+(7.0/3.0)*hxr33(is,i)
! debug
!          write(*,*)is,i,"b = ",b
!          write(*,*)"hxn = ",hxn(is,i)
!          write(*,*)"hxp1 = ",hxp1(is,i)
!          write(*,*)"hxp3 = ",hxp3(is,i)
!          write(*,*)"hxr11 = ",hxr11(is,i)
!          write(*,*)"hxr13 = ",hxr13(is,i)
!          write(*,*)"hxr33 = ",hxr33(is,i)
!          write(*,*)"hxw113 = ",hxw113(is,i)
!          write(*,*)"hxw133 = ",hxw133(is,i)
!          write(*,*)"hxw333 = ",hxw333(is,i)
!
          ftx=ft
          if(ftx.ge.ft_min)then
           ft2 = ftx**2
           gxn(is,i)    = FLR_Hn(ftx,b)
           gxp1(is,i)   = FLR_dHp1(ftx,b)  + ft2*gxn(is,i)
           gxp3(is,i)   = FLR_dHp3(ftx,b)  + gxn(is,i)
           gxr11(is,i)  = FLR_dHr11(ftx,b) + 3.0*ft2*gxp1(is,i)
           gxr13(is,i)  = FLR_dHr13(ftx,b) + (5.0/3.0)*gxp1(is,i)
           gxr33(is,i)  = FLR_dHr33(ftx,b) + (5.0/3.0)*gxp3(is,i)
           gxw113(is,i) = FLR_dHw113(ftx,b)+(7.0/3.0)*gxr11(is,i)
           gxw133(is,i) = FLR_dHw133(ftx,b)+(7.0/3.0)*gxr13(is,i)
           gxw333(is,i) = FLR_dHw333(ftx,b)+(7.0/3.0)*gxr33(is,i)
! debug
!          write(*,*)is,i,"b = ",b
!          write(*,*)"gxn = ",gxn(is,i)
!          write(*,*)"gxp1 = ",gxp1(is,i)
!          write(*,*)"gxp3 = ",gxp3(is,i)
!          write(*,*)"gxr11 = ",gxr11(is,i)
!          write(*,*)"gxr13 = ",gxr13(is,i)
!          write(*,*)"gxr33 = ",gxr33(is,i)
!          write(*,*)"gxw113 = ",gxw113(is,i)
!          write(*,*)"gxw133 = ",gxw133(is,i)
!          write(*,*)"gxw333 = ",gxw333(is,i)
          endif
        enddo
      enddo
!
      CONTAINS
!
      REAL FUNCTION FLR_Hn(ft,b)
!******************************************************************
!     Approxmation to the integral of (J0^2)*Fmaxwellian over a
!     wedge of velocity space -ft < (v_par/v) < ft.
!     b = (k_per*sqrt(T/m)/(eB/mc))**2
!     The approximation has the form
!     Hn = ft*(h0+a1*h1+a2*h2+a3*h3+a4*h4+a5*h5+a6*h6
!             +a7*h7+a8*h8+a9*h9+a10*h10+a11*h11++a12*h12+a13*h13)
!     where a1 - a13 are fit coefficients which are
!     functions of ft and
!     xa1 = (0.25)**4
!     xa2 = (0.5)**4
!     xa3 = (0.75)**4
!     xa4 = (1.0)**4
!     xa5 = 0.25*(1.5)**5
!     xa6 = 0.25*(2.0)**5
!     xa7 = 0.25*(2.5)**5
!     xa8 = 0.25*(3.0)**5
!     xa9 = 0.25*(4.0)**5
!     xa10 = 0.25*(6.0)**5
!     xa11 = 0.25*(9.0)**5
!     xa12 = 0.25*(15.0)**5
!     xa13 = 0.25*(24.0)**5
!     h0 = 1/(1+b)
!     h1 = b/(xa1+b**2)
!     h2 = b/(xa2+b**2)
!     h3 = b/(xa3+b**2)
!     h4 = b/(xa4+b**2)
!     h5 = b**2/(xa5+b**2.5)
!     h6 = b**2/(xa6+b**2.5)
!     h7 = b**2/(xa7+b**2.5)
!     h8 = b**2/(xa8+b**2.5)
!     h9 = b**2/(xa9+b**2.5)
!     h10 = b**2/(xa10+b**2.5)
!     h11 = b**2/(xa11+b**2.5)
!     h12 = b**2/(xa12+b**2.5)
!     h13 = b**2/(xa13+b**2.5)
!******************************************************************
      IMPLICIT NONE
      INTEGER :: i,j,k
      INTEGER,PARAMETER :: nf=40,na=13
      REAL :: y(na),h(na)
      REAL :: a(nf,na)
      REAL :: g(nf)
      REAL :: b2,b25
      REAL :: gt,dg,ca
      REAL :: h0,hs
! inputs
      REAL,INTENT(IN) :: ft,b
!
      data y &
       / 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, &
         6.0, 9.0, 15.0, 24.0 /
      data g &
       / 0.0, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, &
         0.225, 0.25, 0.275, 0.3, 0.325, 0.35, 0.375, 0.4, &
         0.425, 0.45, 0.475, 0.5, 0.525, 0.55, 0.575, 0.6, &
         0.625, 0.65, 0.675, 0.7, 0.725, 0.75, 0.775, 0.8, &
         0.825, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975, 0.995 / 
!
      data a(1:nf,1) &
      /0.00006827619289334661,0.00006739908529205798, &
       0.00006642340225559538,0.00006521725089427256, &
       0.0000633200944587625,0.00006094658676147767, &
       0.00005846319823717128,0.00005589429651048847, &
       0.0000530120714860781,0.00004959640561852401, &
       0.00004557299901515817,0.00004100465479354848, &
       0.00003601632697866836,0.00003072437701496192, &
       0.00002519950920812322,0.00001946243688380542, &
       0.00001349894106825044,7.280449497902018E-6, &
       7.81467538640711E-7,-6.009281872077906E-6, &
       -0.00001308574283124015,-0.00002042487253389656, &
       -0.00002799015843716406,-0.00003573558423659366, &
       -0.00004360896680832663,-0.00005155411928868095, &
       -0.00005951179518009602,-0.00006741949596189416, &
       -0.00007521130873048697,-0.0000828162567642308, &
       -0.0000901578406383946,-0.0000971533116929947, &
       -0.0001037134062249893,-0.0001097425459503292, &
       -0.0001151396744077512,-0.0001197986909616581, &
       -0.0001236110865185249,-0.0001264666625018391, &
       -0.0001282566808816975,-0.0001288505286545868 /
!
      data a(1:nf,2) &
       /0.001917837998461118,0.001894831420943146, &
       0.001863229869005367,0.001815272475878094, &
       0.001761957193383365,0.001698266272655728, &
       0.001615575531500948,0.001513410116041394, &
       0.00139731379569662,0.001272663694123142, &
       0.001141382167906945,0.001002149288014151, &
       0.000852197090052292,0.0006890069730533658, &
       0.0005112080261794595,0.0003186989974295953, &
       0.0001123128948300837,-0.0001066484461681042, &
       -0.0003367983984318343,-0.0005769001774191032, &
       -0.000825906772387593,-0.001082918951414185, &
       -0.00134706776423393,-0.001617387963595768, &
       -0.001892709032091519,-0.002171577917563381, &
       -0.002452215603816502,-0.002732507214908166, &
       -0.003009997786640783,-0.003281932260499747, &
       -0.003545276142994971,-0.003796752493922634, &
       -0.004032878793634551,-0.004250004200407228, &
       -0.004444343075483848,-0.004612034155290956, &
       -0.004749155747978272,-0.004851797135527148, &
       -0.004916102097701667,-0.004937428976052244 /
!
      data a(1:nf,3) &
       /0.004416936495960673,0.004198733934450849, &
       0.003971086739384624,0.003708926908667967, &
       0.003246708666177354,0.002662336696735139, &
       0.002090151211288818,0.001539229712702826, &
       0.000925055186225536,0.0001654075826281448, &
       -0.0007678497544560263,-0.001851503118857133, &
       -0.003038775870953285,-0.004285898148816364, &
       -0.005566240742194187,-0.006872174787164498, &
       -0.00820958235439169,-0.00959015165104458, &
       -0.01102487565531641,-0.01251981890794384, &
       -0.01407503433860473,-0.01568474990971705, &
       -0.01733882935018275,-0.01902442479016823, &
       -0.02072740520430564,-0.02243333915151355, &
       -0.0241280023381104,-0.0257974198545245, &
       -0.0274279395300456,-0.02900560713734305, &
       -0.03051604924365631,-0.03194415363020924, &
       -0.03327388215571598,-0.03448819645817093, &
       -0.03556915300848081,-0.03649777861334273, &
       -0.0372546390596217,-0.03781983035637242, &
       -0.03817337756393268,-0.03829054341125325 /
!
      data a(1:nf,4) &
      /-0.06434740889079478,-0.06555002710632141, &
       -0.06716091272446522,-0.06954428615096775, &
       -0.07226395352344956,-0.07549890029449482, &
       -0.079568296276747,-0.0844668744919215, &
       -0.0899450417206399,-0.0957542938024496, &
       -0.1017803397879442,-0.1080365500818661, &
       -0.1145929327871927,-0.1215080461724494, &
       -0.1287927933777406,-0.1364058246606678, &
       -0.14426796169964,-0.1522824579899923, &
       -0.1603522527399797,-0.1683914667030247, &
       -0.1763287626439948,-0.1841074578501872, &
       -0.1916823025860401,-0.1990157211922125, &
       -0.2060745743641509,-0.2128280027713059, &
       -0.2192464162408039,-0.2253015883567525, &
       -0.2309665351919952,-0.2362171252884824, &
       -0.2410321706414794,-0.2453939214965719, &
       -0.2492880428797856,-0.2527031253719028, &
       -0.255629588807376,-0.2580589003271746, &
       -0.2599809348080443,-0.2613824761385884, &
       -0.2622445622423197,-0.2625277808506313 /
!
      data a(1:nf,5) &
      /-0.03624280424486432,-0.03723296874210282, &
       -0.03831724944381009,-0.03964501491717279, &
       -0.04171279574879952,-0.04425427975868798, &
       -0.04681804294494354,-0.04934326946419526, &
       -0.05206380537642876,-0.05520350809595954, &
       -0.05880904776156679,-0.06275731278122044, &
       -0.06684507477683072,-0.07087651702787476, &
       -0.07471155461418671,-0.07827452456930438, &
       -0.0815386491009324,-0.0845027765070332, &
       -0.0871716712709526,-0.0895433805919317, &
       -0.0916069153850389,-0.0933438444668102, &
       -0.0947339689063293,-0.0957613533025944, &
       -0.0964192638900892,-0.0967132188886513, &
       -0.0966620301579191,-0.0962968840936858, &
       -0.0956603433776841,-0.094802444777618, &
       -0.0937787264303954,-0.0926473428405014, &
       -0.0914667262922187,-0.0902937373207098, &
       -0.0891824397036881,-0.0881827421345534, &
       -0.0873405788974251,-0.0866971358537724, &
       -0.0862887933912129,-0.0861525246521266 /
!
      data a(1:nf,6) &
      /0.01128655852507744,0.01020872249706667, &
       0.00851651326216379,0.005772464044507234, &
       0.003347846287051066,0.0006932624041000129, &
       -0.003193044807520767,-0.00826498174976763, &
       -0.01366754423319679,-0.01853967034316094, &
       -0.02247035312789522,-0.02549524395997627, &
       -0.02786890651122131,-0.02983641653922668, &
       -0.03150722151204166,-0.03283650291946563, &
       -0.03367458801809213,-0.03383949884464882, &
       -0.03318022427354093,-0.03161903395956742, &
       -0.02916207463310429,-0.02589509750521636, &
       -0.02196374572987713,-0.01754939504207028, &
       -0.0128456485993702,-0.00803891030763563, &
       -0.003294512355188495,0.001250765601437556, &
       0.005494870164759264,0.00936581810109262, &
       0.01282136297566723,0.01584373349342064, &
       0.0184342894837606,0.02060808299382928, &
       0.02238881061745491,0.02380347751967073, &
       0.02487892836818135,0.02563731365940614, &
       0.02609254062011278,0.02624017922376688/
!
      data a(1:nf,7) &
      /0.07716778195454137,0.07577796814735106, &
       0.07487044993891371,0.07456467557427752, &
       0.07228028431565159,0.06942333403698992, &
       0.06839594867212106,0.06946098783699406, &
       0.07119768466824401,0.07217735457826062, &
       0.07194152183425922,0.07100130369363078, &
       0.07032971762937827,0.07082545646004733, &
       0.07298752917554011,0.07683284670564355, &
       0.0819847429704652,0.0878383235022717, &
       0.0937275500858581,0.0990598580406921, &
       0.1033857852503181,0.1064324713021665, &
       0.1080951392106117,0.1084082680038694, &
       0.1075075725715335,0.1055916433148604, &
       0.1028885640650713,0.0996315236035448, &
       0.0960347086041697,0.0922912915764859, &
       0.0885624966343097,0.0849782014931663, &
       0.081639144247468,0.07862090109280684, &
       0.07597868877063493,0.07375318266839013, &
       0.07197497005936664,0.07067008929205547, &
       0.06986444678430698,0.06959931338237935 /
!
      data a(1:nf,8) &
      /0.03252198432252359,0.03187560768572029, &
       0.03028260743034252,0.02725806944799758, &
       0.02589483612780368,0.02495262373469141, &
       0.02213580845248175,0.01759867251178233, &
       0.01343197340274228,0.01169969887580031, &
       0.01321209741178049,0.0174288000038641, &
       0.02300798838619348,0.0284624280925517, &
       0.03262346720408684,0.03484156301372182, &
       0.03497747005449,0.03327617228209583, &
       0.03020581830019144,0.02630536286853378, &
       0.02208622098388894,0.01796422932222187, &
       0.01423599512727149,0.01108031695812263, &
       0.00857513749841089,0.006721258502168342, &
       0.005466786629886941,0.004727319060838244, &
       0.00440878789680525,0.004411592611411729, &
       0.004644132057690319,0.005026076917536618, &
       0.005490209639766908,0.00598245202641033, &
       0.006460797970312576,0.006893101077564756, &
       0.007255556791773701,0.007530200916588709, &
       0.007703156974403311,0.007760615257692991 /
!
      data a(1:nf,9) &
      /0.0911779056244026,0.0881253590690804, &
       0.0853208653531631,0.0827684736944843, &
       0.0791486707825628,0.07625944801228456, &
       0.07659806636179881,0.0805814745902774, &
       0.0866492189583986,0.0926684663335451, &
       0.0970527221721862,0.0991469655426323, &
       0.0990847660430412,0.097444235439649, &
       0.0949245606030262,0.0921342180016143, &
       0.0894963195444687,0.0872389713566548, &
       0.0854311012813858,0.0840393319756445, &
       0.0829719843149831,0.0821214919909893, &
       0.0813882903552908,0.080693826633367, &
       0.07998482247566951,0.0792319905977627, &
       0.07842584418448374,0.07757204590779291, &
       0.0766841544079564,0.07578258460537013, &
       0.0748892444965231,0.0740257750738367, &
       0.07321221977707108,0.07246636730824236, &
       0.07180360905053374,0.07123724272575758, &
       0.07077883562412619,0.07043888330053032, &
       0.07022740961583604,0.07015754151538306/

!
      data a(1:nf,10) &
      /0.06765638260432305,0.06406703169875727, &
       0.06034672675706307,0.05705375385345814, &
       0.05821703315375604,0.063314399750297, &
       0.06899189691530054,0.07267290696161772, &
       0.07369130484563384,0.07268452864881212, &
       0.07067671321051321,0.06850481214800995, &
       0.06663505277387536,0.0652099899437958, &
       0.0641740883509613,0.06339110949150478, &
       0.06272192669260517,0.06206243660420492, &
       0.0613543088335812,0.06057429717681843, &
       0.05972932298989056,0.05883985172638258, &
       0.05793161680300605,0.05702944948409843, &
       0.05615399182853526,0.05532052754557671, &
       0.05453914486145863,0.05381536630118985, &
       0.05315212089313132,0.05254953724630415, &
       0.05200666932318899,0.05152203345183076, &
       0.0510940808406292,0.05072153894530461, &
       0.05040363930566202,0.05014017958500987, &
       0.0499317088414802,0.04977952912351625, &
       0.04968580520521201,0.04965498975246665/
!
      data a(1:nf,11) &
      /0.05699370071348491,0.05149192664499243, &
       0.05065942319058143,0.05579104611023448, &
       0.06101239212577761,0.06269404285014086, &
       0.06147918738660241,0.05932768192334556, &
       0.05756129058611158,0.05655639877010695, &
       0.05613262067648053,0.05595963312614589, &
       0.05577673199589012,0.0554521660364829, &
       0.05496030814801516,0.05433562003366362, &
       0.0536329244956911,0.05290289972570367, &
       0.05217936683016035,0.0514862114387255, &
       0.05082826446135012,0.05020419824440915, &
       0.04960880211228827,0.04903614759550817, &
       0.04848149731707176,0.04794216839334134, &
       0.04741766503005173,0.04690948638857824, &
       0.04642022875428433,0.04595366075598298, &
       0.04551391014910902,0.04510515782097256, &
       0.044731450473408,0.04439657803398839, &
       0.04410405631717926,0.04385717974542047, &
       0.04365913763521101,0.04351316083460527, &
       0.04342269696911762,0.04339286277019279 /
!
      data a(1:nf,12) &
      /0.03664939498385491,0.03237561518152475, &
       0.03813476878412953,0.04030660953309315, &
       0.03887188802507055,0.03730243844922979, &
       0.03679802474932086,0.03693425942632053, &
       0.03708479461328912,0.03696189241002144, &
       0.03657275414972018,0.03604058147594843, &
       0.03548204718394607,0.03496309985507814, &
       0.03450174939043032,0.03408721802067922, &
       0.03369841700236624,0.03331581242012639, &
       0.0329323792962798,0.03253169017150431, &
       0.03212006873120008,0.03170202110919428, &
       0.03128372693263086,0.03087155033826579, &
       0.0304711538031869,0.03008710153749724, &
       0.02972279803982096,0.02938055073186501, &
       0.02906204419282494,0.02876820705144789, &
       0.0284995935552234,0.02825660499916138, &
       0.02803954705347184,0.02784874667744107, &
       0.0276846335155976,0.02754777567441373, &
       0.02743897150030604,0.02735927996606062, &
       0.02731009442136884,0.02729390561702728/
!
      data a(1:nf,13) &
      /0.02975846135747427,0.03236122711191199, &
       0.03175957255928162,0.03023615839610341, &
       0.03024285964529651,0.03054126825305659, &
       0.03040821483874899,0.02995175957230299, &
       0.02946249782742549,0.02909263806845435, &
       0.02884648465085155,0.02866313229854521, &
       0.02848087627563395,0.02826339910374234, &
       0.02800056111964719,0.02769915739103121, &
       0.02737300394614417,0.02703606615031493, &
       0.02669307514727589,0.02636242071503557, &
       0.02604141834126838,0.02573083707347553, &
       0.02543023233784817,0.02513882924118646, &
       0.02485608114950881,0.02458194111914954, &
       0.02431692212366228,0.024062057734766, &
       0.02381865631293933,0.02358829463268834, &
       0.02337270302971941,0.02317351248508457, &
       0.02299231519769485,0.02283061055901999, &
       0.02268981329467284,0.02257127952887872, &
       0.02247636740274565,0.02240649781968627, &
       0.02236322889550153,0.02234896652916776 /
!
      b2 = b*b
      b25 = b**2.5
      h0 = 1.0/(1.0+b)
      do i=1,4
       h(i) = b/( y(i)**4 + b2)
      enddo
      do i=5,na
       h(i) = b2/(0.25*y(i)**5 + b25)       
      enddo
! transform to gt grid
      gt = SQRT(1-ft)
      if(gt.gt.g(nf))gt = g(nf)
! find grid position
      if(gt.le.g(2))then
        i=1
      else
        i = INT(gt/0.025)
      endif
      j = i+1
!      write(*,*)"FLR_Hn",i,gt,ft
! interpolate the coefficients
      dg = (gt - g(i))/(g(j)-g(i))
      hs = h0
      do k=1,na
        ca = a(i,k) + (a(j,k)-a(i,k))*dg
! sum up the terms
        hs = hs + ca*h(k)
      enddo
! final answer
      FLR_Hn = ft*hs
!
      END FUNCTION FLR_Hn
!
!
      REAL FUNCTION FLR_dHp1(ft,b)
!******************************************************************
!     Approxmation to the integral of (J0^2)*Fmaxwellian over a
!     wedge of velocity space -ft < (v_par/v) < ft.
!     b = (k_per*sqrt(T/m)/(eB/mc))**2
!     The approximation has the form
!     Hn = ft*(h0+a1*h1+a2*h2+a3*h3+a4*h4+a5*h5+a6*h6
!             +a7*h7+a8*h8+a9*h9+a10*h10+a11*h11++a12*h12+a13*h13)
!     where a1 - a13 are fit coefficients which are
!     functions of ft and
!     xa1 = (0.25)**4
!     xa2 = (0.5)**4
!     xa3 = (0.75)**4
!     xa4 = (1.0)**4
!     xa5 = 0.25*(1.5)**5
!     xa6 = 0.25*(2.0)**5
!     xa7 = 0.25*(2.5)**5
!     xa8 = 0.25*(3.0)**5
!     xa9 = 0.25*(4.0)**5
!     xa10 = 0.25*(6.0)**5
!     xa11 = 0.25*(9.0)**5
!     xa12 = 0.25*(15.0)**5
!     xa13 = 0.25*(24.0)**5
!     h0 = 1/(1+b)
!     h1 = b/(xa1+b**2)
!     h2 = b/(xa2+b**2)
!     h3 = b/(xa3+b**2)
!     h4 = b/(xa4+b**2)
!     h5 = b**2/(xa5+b**2.5)
!     h6 = b**2/(xa6+b**2.5)
!     h7 = b**2/(xa7+b**2.5)
!     h8 = b**2/(xa8+b**2.5)
!     h9 = b**2/(xa9+b**2.5)
!     h10 = b**2/(xa10+b**2.5)
!     h11 = b**2/(xa11+b**2.5)
!     h12 = b**2/(xa12+b**2.5)
!     h13 = b**2/(xa13+b**2.5)
!******************************************************************
      IMPLICIT NONE
      INTEGER,PARAMETER :: nf=40,na=13
      INTEGER :: i,j,k
      REAL :: y(na),h(na)
      REAL :: a(nf,na)
      REAL :: g(nf)
      REAL :: b2,b25
      REAL :: gt,dg,ca
      REAL :: h0,hs
! inputs
      REAL,INTENT(IN) :: ft,b
!
      data y &
      / 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, &
        6.0, 9.0, 15.0, 24.0 /
      data g &
      / 0.0, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, &
        0.225, 0.25, 0.275, 0.3, 0.325, 0.35, 0.375, 0.4, &
        0.425, 0.45, 0.475, 0.5, 0.525, 0.55, 0.575, 0.6, &
        0.625, 0.65, 0.675, 0.7, 0.725, 0.75, 0.775, 0.8, &
        0.825, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975, 0.995 /
!
      data a(1:nf,1) &
       /0,-1.930999091095923E-6,-2.960991528348877E-6, &
       -6.154174511454453E-6,-0.00001114642641835047, &
       -0.00001496341695015972,-0.00001777711412125118, &
       -0.00002189234582485312,-0.00002885391045732016, &
       -0.00003844640850068675,-0.00004939871893680908, &
       -0.00006045236037147072,-0.000070990727933951, &
       -0.0000811127772586047,-0.0000913730276388475, &
       -0.000102444966057082,-0.0001148663051755398, &
       -0.0001289201953026221,-0.0001446331835174276, &
       -0.0001618429427948575,-0.0001802880130045098, &
       -0.0001996876099350452,-0.0002197944806753826, &
       -0.0002404183476049577,-0.0002614251248907947, &
       -0.0002827204205832545,-0.0003042255829558891, &
       -0.000325854892310157,-0.0003474920599192017, &
       -0.0003689807865177399,-0.0003901066234334019, &
       -0.0004106169137956872,-0.0004301838415965056, &
       -0.0004484654074302117,-0.0004650668595640494, &
       -0.000479575280909874,-0.0004915659226616756, &
       -0.0005006175919946134,-0.0005063189593678531, &
       -0.0005082157906882695 /
!
      data a(1:nf,2) &
        /0,-0.00004194525666284548, &
        -0.0001273156308883975,-0.0002056141090789176, &
        -0.0002857559910660781,-0.0004389428374474458, &
        -0.0006617588108732073,-0.000900134130104115, &
        -0.001118233049555858,-0.001322327045180353, &
        -0.001544024892708917,-0.001814796975413207, &
        -0.002150943267801318,-0.002551783869108248, &
        -0.003005781915370591,-0.003498463624428219, &
        -0.004018368587377641,-0.004559734264079984, &
        -0.005122388308748703,-0.005709968162943081, &
        -0.006327610674374817,-0.006979882323169755, &
        -0.007669355681763603,-0.00839589743212053, &
        -0.0091565494387851,-0.00994580575745742, &
        -0.01075609567750988,-0.0115782655678188, &
        -0.01240213088038841,-0.01321673786472565, &
        -0.01401074445648378,-0.01477237535347542, &
        -0.01548985492815352,-0.01615100847176631, &
        -0.0167436906236944,-0.01725562758025234, &
        -0.01767455366934878,-0.01798831362495496, &
        -0.01818491482121334,-0.01825013519611904 /
!
      data a(1:nf,3) &
       / 0,-0.0004949458659930586, &
        -0.0006223389545283466,-0.001425278848220574, &
        -0.002777762953477727,-0.003592632754394858, &
        -0.003924601161759841,-0.004621076696713783, &
        -0.006256088990015495,-0.00875017365645226, &
        -0.01162549543224132,-0.0144032952638076, &
        -0.01684202469464682,-0.01896846946059326, &
        -0.02098241046601086,-0.02313062220828303, &
        -0.02561028533866318,-0.02852321735637834, &
        -0.0318743399300858,-0.03559723690450516, &
        -0.03958861135502096,-0.04373956572105996, &
        -0.04795689127872948,-0.05217324721801679, &
        -0.05634803211077842,-0.06046210119442085, &
        -0.06450945249074329,-0.0684891545861474, &
        -0.07239683993958543,-0.07622086652002834, &
        -0.07993823320236293,-0.0835150186803319, &
        -0.0869050660668986,-0.0900543242751263, &
        -0.092900899412087,-0.0953788527732187, &
        -0.097420047849393,-0.0989564183116274, &
        -0.0999230523070371,-0.1002441867837462 /
!
      data a(1:nf,4) &
       / 0,-0.002257777688720037,-0.006305363728668376, &
        -0.01034798017113135,-0.01466436256998603, &
        -0.0219793722487998,-0.03209517213917592, &
        -0.04276949973637403,-0.05244246580218026, &
        -0.06121622935485103,-0.07020559130958571, &
        -0.080519292080036,-0.0926472400237077, &
        -0.1063790403546331,-0.1210487879360195, &
        -0.135860537033621,-0.1501394639530332, &
        -0.1634530299978431,-0.1756183677292047, &
        -0.1866395914432249,-0.1966219016087455, &
        -0.2056934177153338,-0.2139525405581872, &
        -0.221443738531069,-0.2281571285321903, &
        -0.2340436664618902,-0.2390377750439819, &
        -0.2430789850704455,-0.2461338545152663, &
        -0.2482055663196422,-0.2493428478278852, &
        -0.2496376780529274,-0.2492246313952511, &
        -0.2482680973602374,-0.2469565590532087, &
        -0.2454888837408556,-0.2440639664115824, &
        -0.2428687887653431,-0.2420639893766472, &
        -0.2417879747636104 /
!
      data a(1:nf,5) &
       / 0,-0.002144312832523405,-0.003282843034891315, &
        -0.006741411445268095,-0.0120774728563809, &
        -0.01578560102424071,-0.01791575569806934, &
        -0.02103521175655153,-0.02684614094396873, &
        -0.03494085352248321,-0.04359567240931105, &
        -0.05104908266612131,-0.0562885330954339, &
        -0.05917467864730624,-0.06015117398261926, &
        -0.05984564077307062,-0.05875976743154406, &
        -0.05712158212172106,-0.05488256493166144, &
        -0.05180545830350876,-0.04758219094596531, &
        -0.04194204832100806,-0.03472512233475206, &
        -0.02591626376634247,-0.015644683580604, &
        -0.004159620274591441,0.00820683006251159, &
        0.02107876801389288,0.03407582136073201, &
        0.04683903727868152,0.0590533499257193, &
        0.07045668579590236,0.0808458131119755, &
        0.0900697306378407,0.098024779618702, &
        0.1046416985177536,0.1098735098435689, &
        0.1136813486258382,0.1160178052559541, &
        0.1167845253426408 /
!
      data a(1:nf,6) &
       / 0,-0.001733236253032295,-0.007594073959997902, &
        -0.01068792943290476,-0.01173058924414864, &
        -0.01927159543066403,-0.03257681082882676, &
        -0.04413538432869502,-0.04850613282741194, &
        -0.04574402483725617,-0.03944287876935876, &
        -0.03339043151790365,-0.02944655693792947, &
        -0.02717792694860531,-0.02464215853843571, &
        -0.01951717274411607,-0.01002146884643456, &
        0.004602913165110195,0.02402815898572942, &
        0.04706879795934987,0.07203598980710533, &
        0.0971025256697242,0.120594334820713, &
        0.1411804546140524,0.15796222501998, &
        0.1704810270734963,0.178670379140822, &
        0.1827794420597866,0.1832792857413, &
        0.1807843764157766,0.1759749689355863, &
        0.169542902741839,0.1621476910360063, &
        0.1543952652268083,0.1468212677682953, &
        0.1398898636187544,0.1339956510623273, &
        0.1294701463104257,0.1265934481629084, &
        0.1256324484135063 /
!
      data a(1:nf,7) &
       / 0,-0.003527983222841474, &
        -0.00004831085369737476,-0.004756774305300199, &
        -0.01541520139419636,-0.01334101436680315, &
        0.001218078785539171,0.01379655191457437, &
        0.01413077975140365,0.003667368224953139, &
        -0.00837867150763902,-0.01230460416243886, &
        -0.002903730724170366,0.01927092283480719, &
        0.04930697836054414,0.0803705152020502, &
        0.1059511817375929,0.1213152536294332, &
        0.1240790969285497,0.114089159552476, &
        0.0929112004202932,0.06316328529183107, &
        0.02789090820620711,-0.00992522289453935, &
        -0.04769584609522803,-0.0833949235688121, &
        -0.115608675256831,-0.1435026774642829, &
        -0.1667290111458289,-0.1853311035321232, &
        -0.1996278447543471,-0.2101191646227736, &
        -0.2174006427959979,-0.2221039119533069, &
        -0.2248457626736612,-0.2261976054603467, &
        -0.2266622953907129,-0.226658125500399, &
        -0.2265132290226014,-0.226440981014999 /
!
      data a(1:nf,8) &
       / 0,-0.0004069673455518104,-0.00900496080525841, &
        -0.00951171692779561,-0.003369577324746711, &
        -0.01008732325445694,-0.02793147683269029, &
        -0.0388696504512379,-0.02964969295642721, &
        -0.001243576427339135,0.03486111183554177, &
        0.06490112340939137,0.07915217961432494, &
        0.07431478617202438,0.05277793522926171, &
        0.02037900572399001,-0.01599334230597047, &
        -0.05020302004671966,-0.07785553394622869, &
        -0.0965629487693895,-0.1057379603699174, &
        -0.1061301802107192,-0.0993110132341835, &
        -0.087214445237185,-0.07178899786719991, &
        -0.05477141525768695,-0.03756734730057455, &
        -0.02121154797874673,-0.006392097308036287, &
        0.006510923726594563,0.01735631306708534, &
        0.02618086388428061,0.03313849993138422, &
        0.03845460474644287,0.04238423274854166, &
        0.04518438418894524,0.04709145785845059, &
        0.0483040704093707,0.04897677093826909, &
        0.04918462614440588 /
!
      data a(1:nf,9) &
       / 0,-0.006351126139988073,-0.00701003729386106, &
        -0.01252056869487079,-0.01891579902576641, &
        -0.00945435616320531,0.01633297079349021, &
        0.04395566414018533,0.05893419348768764, &
        0.05569695313187913,0.03746920138705813, &
        0.01189411288153407,-0.01331381643501007, &
        -0.03286107285662678,-0.04450288667237289, &
        -0.04850108724749244,-0.04663401809397132, &
        -0.0412444492036532,-0.03455195156344903, &
        -0.02827899485009187,-0.02352099413983975, &
        -0.0207929854789008,-0.02015605257406988, &
        -0.02136898590164682,-0.02402731184503589, &
        -0.02767220122661359,-0.0318652307154578, &
        -0.03623282457453952,-0.04048337750422562, &
        -0.04441339857303793,-0.04789728625573435, &
        -0.05087489555514178,-0.05333584077219736, &
        -0.05530540650204516,-0.05683071729203423, &
        -0.05796989587650346,-0.05878271363041665, &
        -0.05932234624298428,-0.05963264967096353, &
        -0.05973045210817956 /
!
      data a(1:nf,10) &
       / 0,-0.005645165914549386,-0.01303234076222609, &
        -0.005724039695679197,0.01815741807803584, &
        0.03619860407391771,0.03516855255739328, &
        0.01893821783208634,-0.001470537132467038, &
        -0.01729803744204529,-0.02518246546732237, &
        -0.0260398522601473,-0.02271866050273372, &
        -0.01813877954860406,-0.01435383282610346, &
        -0.01235500198075294,-0.01228017954061346, &
        -0.0137528470519161,-0.01619142210692707, &
        -0.01901800996943862,-0.0217794513175185, &
        -0.0241800503637666,-0.02607234911039552, &
        -0.02742380307223777,-0.028277729692998, &
        -0.02871841758195068,-0.02884444597229968, &
        -0.02875074811272316,-0.02851986032101957, &
        -0.02821632429994077,-0.02788757243237246, &
        -0.02756550352801054,-0.02726995936712024, &
        -0.02701136094151338,-0.0267941587215538, &
        -0.02661896270299162,-0.02648456922987775, &
        -0.02638971853872839,-0.02633188086790837, &
        -0.02631323881761514 /
!
      data a(1:nf,11) &
       / 0,-0.007757811267608873,0.00842653731688701, &
        0.02499866208531534,0.01822303046942917, &
        -0.0005013209021285298,-0.01413767476888905, &
        -0.01708647660151259,-0.01280756385345587, &
        -0.006833211684983275,-0.002807705766008619, &
        -0.001830204708260426,-0.003334982601312177, &
        -0.006150412369837154,-0.00920004184066017, &
        -0.01179686839345535,-0.01366255970052648, &
        -0.01481764368265418,-0.01544195932675034, &
        -0.01576981293099308,-0.01600519868900191, &
        -0.01629563779234062,-0.01672304681597176, &
        -0.01731267057872737,-0.01804888424073965, &
        -0.01889172826835179,-0.01979117010639083, &
        -0.02069751177924051,-0.02156758022983927, &
        -0.0223686481992779,-0.02307900609638849, &
        -0.02368733664487918,-0.02419112765781742, &
        -0.02459465187466544,-0.02490682531806957, &
        -0.02513918797020957,-0.02530406930214788, &
        -0.02541277265781793,-0.02547488944478554, &
        -0.02549438289422755 /
!
      data a(1:nf,12) &
       / 0,0.00836432994430984,0.01474526217669311, &
        -0.001158550730979778,-0.00968842798062101, &
        -0.006326200318437962,-0.0004503779639845268, &
        0.001827553302003046,0.0002935478441003649, &
        -0.002809883286850672,-0.005527719248852609, &
        -0.007059073509546865,-0.007499989413092357, &
        -0.007333101221443564,-0.007048157025783568, &
        -0.006967211502022494,-0.007219386104956538, &
        -0.007789137152399589,-0.00858832575606895, &
        -0.00948906363222313,-0.01039271646635876, &
        -0.01122379025109202,-0.01193844115822858, &
        -0.01251980775068381,-0.0129703468783664, &
        -0.0133043085245006,-0.01354120818902937, &
        -0.01370144335737454,-0.01380394794462716, &
        -0.01386445834477761,-0.01389546651374912, &
        -0.01390642452150214,-0.01390437260199559, &
        -0.01389434631951948,-0.01388011056440104, &
        -0.01386448497541121,-0.01384970187289114, &
        -0.01383765078695345,-0.01382951221230162, &
        -0.01382677022888323 /
!
      data a(1:nf,13) &
       / 0,0.005067318298867268,-0.006386882621518983, &
        -0.002545894837045421,0.001862106998487167, &
        0.0006745433991216182,-0.002461288343219566, &
        -0.0041932474756031,-0.004072253513993561, &
        -0.003110581064337257,-0.002331015073283393, &
        -0.002201308116790514,-0.002708040950352473, &
        -0.003604975663801907,-0.004617916654148874, &
        -0.005549883923172879,-0.006305894020107905, &
        -0.006873934570591084,-0.007285627859546998, &
        -0.007609277747774316,-0.007893101760097637, &
        -0.00817754302711777,-0.00848590987813958, &
        -0.00882590095988645,-0.00919361074942149, &
        -0.00957771060129377,-0.00996395267157994, &
        -0.01033798749492546,-0.01068738469464519, &
        -0.01100284553441664,-0.01127859625905964, &
        -0.01151206535500384,-0.01170355950440385, &
        -0.01185557243281633,-0.01197212477309358, &
        -0.01205806679914275,-0.01211845123050876, &
        -0.01215790597504746,-0.01218021096394671, &
        -0.01218718760671441 /
!
      b2 = b*b
      b25 = b**2.5
      h0 = 0.0
      do i=1,4
       h(i) = b/( y(i)**4 + b2)
      enddo
      do i=5,na
       h(i) = b2/(0.25*y(i)**5 + b25)       
      enddo
! transform to gt grid
      gt = SQRT(1-ft)
      if(gt.gt.g(nf))gt = g(nf)
! find grid position
      if(gt.le.g(2))then
        i=1
      else
        i = INT(gt/0.025)
      endif
      j = i+1
! interpolate the coefficients
      dg = (gt - g(i))/(g(j)-g(i))
      hs = h0
      do k=1,na
        ca = a(i,k) + (a(j,k)-a(i,k))*dg
! sum up the terms
        hs = hs + ca*h(k)
      enddo
! final answer
      FLR_dHp1 = (ft**3)*hs
!
      END FUNCTION FLR_dHp1
!
!
      REAL FUNCTION FLR_dHp3(ft,b)
!******************************************************************
!     Approxmation to the integral of (J0^2)*Fmaxwellian over a
!     wedge of velocity space -ft < (v_par/v) < ft.
!     b = (k_per*sqrt(T/m)/(eB/mc))**2
!     The approximation has the form
!     Hn = ft*(h0+a1*h1+a2*h2+a3*h3+a4*h4+a5*h5+a6*h6
!             +a7*h7+a8*h8+a9*h9+a10*h10+a11*h11++a12*h12+a13*h13)
!     where a1 - a13 are fit coefficients which are
!     functions of ft and
!     xa1 = (0.25)**4
!     xa2 = (0.5)**4
!     xa3 = (0.75)**4
!     xa4 = (1.0)**4
!     xa5 = 0.25*(1.5)**5
!     xa6 = 0.25*(2.0)**5
!     xa7 = 0.25*(2.5)**5
!     xa8 = 0.25*(3.0)**5
!     xa9 = 0.25*(4.0)**5
!     xa10 = 0.25*(6.0)**5
!     xa11 = 0.25*(9.0)**5
!     xa12 = 0.25*(15.0)**5
!     xa13 = 0.25*(24.0)**5
!     h0 = 1/(1+b)
!     h1 = b/(xa1+b**2)
!     h2 = b/(xa2+b**2)
!     h3 = b/(xa3+b**2)
!     h4 = b/(xa4+b**2)
!     h5 = b**2/(xa5+b**2.5)
!     h6 = b**2/(xa6+b**2.5)
!     h7 = b**2/(xa7+b**2.5)
!     h8 = b**2/(xa8+b**2.5)
!     h9 = b**2/(xa9+b**2.5)
!     h10 = b**2/(xa10+b**2.5)
!     h11 = b**2/(xa11+b**2.5)
!     h12 = b**2/(xa12+b**2.5)
!     h13 = b**2/(xa13+b**2.5)
!******************************************************************
      IMPLICIT NONE
      INTEGER,PARAMETER :: nf=40,na=13
      INTEGER :: i,j,k
      REAL :: y(na),h(na)
      REAL :: a(nf,na)
      REAL :: g(nf)
      REAL :: b2,b25
      REAL :: gt,dg,ca
      REAL :: h0,hs
! inputs
      REAL,INTENT(IN) :: ft,b
!
      data y &
       / 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, &
         6.0, 9.0, 15.0, 24.0 /
      data g &
       / 0.0, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, &
         0.225, 0.25, 0.275, 0.3, 0.325, 0.35, 0.375, 0.4, &
         0.425, 0.45, 0.475, 0.5, 0.525, 0.55, 0.575, 0.6, &
         0.625, 0.65, 0.675, 0.7, 0.725, 0.75, 0.775, 0.8, &
         0.825, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975, 0.995 /
!
      data a(1:nf,1) &
       /-0.0002868056294331331,-0.0002875811537591722, &
        -0.0002881695050265826,-0.0002896843391151882, &
        -0.0002916904269459931,-0.0002933380065348055, &
        -0.0002949384306718349,-0.0002972833835205653, &
        -0.0003007291595217119,-0.0003050393271926225, &
        -0.0003097220956871584,-0.0003143899995580001, &
        -0.0003189217612157989,-0.0003234357443871498, &
        -0.0003281719785246126,-0.0003333708778457172, &
        -0.0003391951158999734,-0.0003457031689108731, &
        -0.0003528619389890245,-0.0003605795765590993, &
        -0.0003687399655396059,-0.0003772310913434102, &
        -0.0003859611117457418,-0.0003948633104879826, &
        -0.000403892603720068,-0.0004130171354556294, &
        -0.0004222081447110765,-0.0004314305973136029, &
        -0.0004406351655552077,-0.0004497540155887716, &
        -0.0004586980234852225,-0.0004673568091365604, &
        -0.0004756003459592098,-0.0004832819693103312, &
        -0.0004902419334322565,-0.0004963123797404506, &
        -0.0005013221108640335,-0.0005050961588005414, &
        -0.0005074728023367473,-0.0005082629513683144/
!
      data a(1:nf,2) &
       /-0.0101208580437598,-0.01014481030575884, &
        -0.01018383000601481,-0.01022219869439501, &
        -0.01027035924740222,-0.01034870160364673, &
        -0.01045006191146344,-0.0105558696455702, &
        -0.01065795607833994,-0.01076231550021279, &
        -0.01088106625271674,-0.01102386729336503, &
        -0.01119404219185928,-0.01138918771674483, &
        -0.0116039544329392,-0.01183289804024278, &
        -0.01207229568563148,-0.0123207229211148, &
        -0.01257869871355472,-0.01284784422238216, &
        -0.0131300035969187,-0.01342651056979119, &
        -0.01373775106687112,-0.01406299480311878, &
        -0.01440043408780798,-0.01474734746446646, &
        -0.01510031432141053,-0.0154554223024956, &
        -0.01580845626257836,-0.01615500998828404, &
        -0.01649057879407347,-0.0168106010988151, &
        -0.01711047653157854,-0.01738556873158428, &
        -0.01763121028584308,-0.01784268780360466, &
        -0.01801530419113404,-0.01814428117225634, &
        -0.01822501268990321,-0.01825177303289219/
!
      data a(1:nf,3) &
       /-0.06017478726489287,-0.06034456260662489, &
        -0.06042293158534254,-0.06077873337316765, &
        -0.06125702160568059,-0.06154670069267421, &
        -0.06175900055964754,-0.06218494189356277, &
        -0.0629568508855498,-0.06398660265535944, &
        -0.0650890707142704,-0.06611636344459421, &
        -0.06701952257324696,-0.0678397735652097, &
        -0.06866489285578081,-0.06958366964443525, &
        -0.07065618018334854,-0.07190349621608494, &
        -0.07331221618575131,-0.07484696916433946, &
        -0.0764636855730334,-0.07812078629995476, &
        -0.07978574291611052,-0.0814374186105153, &
        -0.0830651466618845,-0.0846658686700284, &
        -0.0862405411037517,-0.0877907779045641, &
        -0.0893158988477271,-0.090811466221972, &
        -0.0922682071078722,-0.0936719612726272, &
        -0.095004157323956,-0.0962426344953322, &
        -0.0973625676046369,-0.0983377276087636, &
        -0.0991408352116349,-0.0997457957429468, &
        -0.1001261102915091,-0.1002524836142409/
!
      data a(1:nf,4) &
       /-0.2309847438764876,-0.2315114977755848, &
        -0.232501425640584,-0.2332495936659627, &
        -0.2341358206878437,-0.2359373047551776, &
        -0.2383502261623568,-0.2406051149760865, &
        -0.2423333394159475,-0.243727566836771, &
        -0.245226734582312,-0.2471717998527675, &
        -0.2496468336944673,-0.2525010906994787, &
        -0.2554612640574402,-0.2582492566332228, &
        -0.2606597109545827,-0.2625876997036404, &
        -0.2640182425274357,-0.2649951368464562, &
        -0.2655877006440267,-0.26586271651991, &
        -0.265868234068929,-0.2656281788507019, &
        -0.265145305623525,-0.2644090647734432, &
        -0.2634052280740651,-0.2621247374706666, &
        -0.2605713285221448,-0.2587650287080486, &
        -0.2567444869450332,-0.2545663803011062, &
        -0.2523032279127653,-0.2500400953147265, &
        -0.2478708008680447,-0.2458931119783777, &
        -0.2442058563725805,-0.2429015990773533, &
        -0.2420675550572284,-0.2417880531051465/
!
      data a(1:nf,5) &
       /-0.01487702637412047,-0.01496647668479415, &
        -0.01466094692896006,-0.014972540790499, &
        -0.01540527517197185,-0.01495444670013606, &
        -0.01393468082868539,-0.01323458980930192, &
        -0.01324003850835957,-0.01362599380725006, &
        -0.01374565628822912,-0.01306407787666781, &
        -0.01136441432478322,-0.00872750659979418, &
        -0.005395978690390215,-0.001628911985786109, &
        0.002394276290798997,0.006607348199373491, &
        0.01104205484460605,0.01578607200154645, &
        0.02093808641853327,0.02656966435475682, &
        0.03270308095431584,0.03930413904016015, &
        0.04628702444861925,0.05352675545436241, &
        0.06087494182059941,0.0681752307233167, &
        0.07527760357099771,0.0820470199997801, &
        0.0883704591079797,0.0941594107094766, &
        0.0993495811874999,0.1038983319380035, &
        0.1077805843062816,0.1109830484338663, &
        0.1134994083284368,0.1153212809340351, &
        0.116436311780189,0.1168015764157358/
!
      data a(1:nf,6) &
       /0.0634113482479667,0.06371375051423698, &
        0.06310152741164138,0.06422251048125092, &
        0.06592539439511852,0.06574870049813625, &
        0.06466463521299261,0.06520053765250568, &
        0.06861111868177238,0.07427716749202572, &
        0.0806926528344103,0.0866075701349628, &
        0.091590067345508,0.0959830177626426, &
        0.1005356409514118,0.1059931313416857, &
        0.1128107607192332,0.1210402588025594, &
        0.1303598640614494,0.1401930190702618, &
        0.1498494143588391,0.1586583157163637, &
        0.1660634156253317,0.1716765258391989, &
        0.1752939207266033,0.176884811640769, &
        0.1765625679518662,0.1745488190985015, &
        0.1711331887341624,0.1666439218206136, &
        0.161417904304662,0.1557812062441614, &
        0.150035904275458,0.1444528491687018, &
        0.1392691145312419,0.1346896977024388, &
        0.1308893346508469,0.12802187106561, &
        0.1262167974565677,0.1256170094524324/
!
      data a(1:nf,7) &
       /-0.0615291252932581,-0.06192923753744736, &
        -0.06034976793978348,-0.06197805854545935, &
        -0.06430831851941375,-0.06204730874632855, &
        -0.0569524829601291,-0.05408729979183575, &
        -0.05590947533676593,-0.06090284005180366, &
        -0.06561660275035353,-0.06716646901420011, &
        -0.06458070396981354,-0.0588395912698829, &
        -0.05213795219643817,-0.04695930854344886, &
        -0.04534270143077448,-0.04848698865156153, &
        -0.05667163378730945,-0.06940353155569939, &
        -0.0856595657877482,-0.1041595311611984, &
        -0.1235941363470313,-0.1427895966031381, &
        -0.1608042500034685,-0.1769669855727497, &
        -0.1908727884749794,-0.202352705776142, &
        -0.21142247436486,-0.2182423626156722, &
        -0.2230651068954094,-0.2261975413760589, &
        -0.2279688905766291,-0.2287069907272046, &
        -0.228721523617018,-0.2282936440506118, &
        -0.2276670828933931,-0.2270510518176625, &
        -0.2266030821011039,-0.2264445072399632/
!
      data a(1:nf,8) &
       /-0.01460015606988496,-0.01433994477025534, &
        -0.01617779535043274,-0.01435257338140886, &
        -0.01143452984619481,-0.01309085907327989, &
        -0.01712944945092585,-0.01758557490912519, &
        -0.01150106705415536,-0.0007494675431253071, &
        0.01010620972120545,0.01672565356329256, &
        0.01678036805511051,0.01030724892916775, &
        -0.000890753596719823,-0.01418466039694179, &
        -0.02696600230061974,-0.03719980974196325, &
        -0.04366190060614404,-0.04592290346814614, &
        -0.04419963834406588,-0.03913494151528692, &
        -0.03158926291283939,-0.02247021380638297, &
        -0.01261381917439763,-0.002715558874454864, &
        0.006698073102144609,0.01527058337689374, &
        0.02279152861393773,0.02917754713056737, &
        0.03443710130904834,0.0386434621974366, &
        0.04191024085660712,0.04437149086995909, &
        0.04616635348017084,0.04742820968557335, &
        0.04827467460065945,0.04880867547921657, &
        0.04909997144919772,0.04918957605413321/
!
      data a(1:nf,9) &
       /-0.03268211293774046,-0.03283837471710354, &
        -0.03126702352807229,-0.03149936613845066, &
        -0.03130396623553432,-0.02631830836425632, &
        -0.01815783454354737,-0.01204979796709452, &
        -0.01181881861940296,-0.01783228333061884, &
        -0.02781926750173951,-0.03865121531702853, &
        -0.04774469553230983,-0.05368487977541645, &
        -0.0562152831547101,-0.05589632544639403, &
        -0.05368796707793209,-0.05060470343196925, &
        -0.04749591764475936,-0.04495237024191775, &
        -0.04329164625700513,-0.0426055967323814, &
        -0.04282414794825351,-0.0437807534884101, &
        -0.04526722161663855,-0.04707367403803413, &
        -0.04901381278900674,-0.05093849982524341, &
        -0.05273754396795742,-0.05434099142534901, &
        -0.05571171232753943,-0.05683910746033794, &
        -0.05773216335541091,-0.05841301124955212, &
        -0.05891138265540217,-0.05926012913657658, &
        -0.05949074770768033,-0.05963405040886491, &
        -0.05971043052267028,-0.05973360020159963/
!
      data a(1:nf,10) &
       /-0.02298570615364132,-0.02252332134664494, &
        -0.02254613881864249,-0.01795650552692745, &
        -0.0108873631669562,-0.00853739412224544, &
        -0.01306363173516015,-0.02136381191849219, &
        -0.02922413698772504,-0.03410392380466562, &
        -0.03557863949657546,-0.03455937833500515, &
        -0.03236027087296684,-0.03008898804014205, &
        -0.02840967266777042,-0.02755737285766853, &
        -0.02746638112662963,-0.0279164571965106, &
        -0.02864979999785831,-0.02943785304900157, &
        -0.03012038492668657,-0.03060548541042286, &
        -0.03085997227785798,-0.03089263810553368, &
        -0.03073742664610525,-0.03043958188286178, &
        -0.03004578791622565,-0.02959789685902671, &
        -0.0291308260692537,-0.0286704809442081, &
        -0.02823513129496223,-0.02783651880128571, &
        -0.02748135150607498,-0.02717280566771195, &
        -0.02691184676159109,-0.02669822256787656, &
        -0.02653180221019968,-0.02641100548279871, &
        -0.02633751850511063,-0.02631346979493314/
!
      data a(1:nf,11) &
       /-0.01828230664388564,-0.01723292093175429, &
        -0.01132504465009096,-0.00934679892988526, &
        -0.01528716112300812,-0.0228464071031182, &
        -0.02672281587162796,-0.02637200856749093, &
        -0.02387783279196622,-0.02137798118732101, &
        -0.01997943083276751,-0.01981634158945404, &
        -0.02050222909606981,-0.02153503452203458, &
        -0.02252292451092742,-0.02325362921230601, &
        -0.02367268399873643,-0.02382723968323626, &
        -0.02380778165088634,-0.02371529536301569, &
        -0.02362549778538901,-0.02358820664512237, &
        -0.02362603565351737,-0.0237400929241458, &
        -0.023917287918038,-0.02413720041210286, &
        -0.02437751903509851,-0.02461787695447153, &
        -0.02484160140465608,-0.02503721677775493, &
        -0.02519811109604487,-0.02532204816719052, &
        -0.02541038927994028,-0.02546716097076922, &
        -0.02549815086004078,-0.02551007360410793, &
        -0.0255096887462481,-0.02550392842589301, &
        -0.02549776666451141,-0.02549531042372639/
!
      data a(1:nf,12) &
       /-0.0117718769470323,-0.006155963076068777, &
        -0.007933657063456124,-0.01476426117029363, &
        -0.01670502607753369,-0.0145883642535117, &
        -0.01237851342740976,-0.01183684476062507, &
        -0.01259977055658306,-0.01370224520534855, &
        -0.014487865077719,-0.01477301854980721, &
        -0.0146756335811828,-0.01440932514231363, &
        -0.01415479486507491,-0.01401312317819495, &
        -0.01401036777614579,-0.01412312543731965, &
        -0.014311464826639,-0.01451648732303002, &
        -0.01470543962730434,-0.01485356278690819, &
        -0.01494913480895602,-0.01499064979056821, &
        -0.01498334658876554,-0.01493605767107214, &
        -0.01485880813248459,-0.01476117351736615, &
        -0.01465168090084777,-0.01453712043036627, &
        -0.01442274110022205,-0.01431251545183906, &
        -0.01420935464297278,-0.0141154250869964, &
        -0.01403241006644194,-0.01396168995244758, &
        -0.01390461205661277,-0.01386213296284805, &
        -0.01383574596415517,-0.01382702374402716/
!
      data a(1:nf,13) &
       /-0.00944680887541073,-0.00952560731960679, &
        -0.01297313737190309,-0.01070729761586587, &
        -0.00929922025626993,-0.00997218167017799, &
        -0.01100870712104808,-0.01135687257110129, &
        -0.0110664087726277,-0.01059014997624408, &
        -0.01027529465371409,-0.01023562861771975, &
        -0.01041987147434059,-0.01071524633071631, &
        -0.01101860759977515,-0.01126623125467858, &
        -0.01143526970585649,-0.01153207451056483, &
        -0.01157195930591625,-0.01159141951029263, &
        -0.01160565275266543,-0.01162876814204084, &
        -0.01166761863620499,-0.01172304459613098, &
        -0.0117917594389806,-0.01186824989721547, &
        -0.01194635277515668,-0.01202041351364881, &
        -0.01208585772844124,-0.01213964646024652, &
        -0.01218033272913241,-0.01220776899214369, &
        -0.01222298976357256,-0.01222789618071573, &
        -0.01222496998758532,-0.01221699239780094, &
        -0.01220677657512377,-0.01219706719147791, &
        -0.01219006069959786,-0.0121875871031617/
!
      b2 = b*b
      b25 = b**2.5
      h0 = 0.0
      do i=1,4
       h(i) = b/( y(i)**4 + b2)
      enddo
      do i=5,na
       h(i) = b2/(0.25*y(i)**5 + b25)       
      enddo
! transform to gt grid
      gt = SQRT(1-ft)
      if(gt.gt.g(nf))gt = g(nf)
! find grid position
      if(gt.le.g(2))then
        i=1
      else
        i = INT(gt/0.025)
      endif
      j = i+1
! interpolate the coefficients
      dg = (gt - g(i))/(g(j)-g(i))
      hs = h0
      do k=1,na
        ca = a(i,k) + (a(j,k)-a(i,k))*dg
! sum up the terms
        hs = hs + ca*h(k)
      enddo
! final answer
      FLR_dHp3 = ft*hs
!
      END  FUNCTION FLR_dHp3
!
!
      REAL FUNCTION FLR_dHr11(ft,b)
!******************************************************************
!     Approxmation to the integral of (J0^2)*Fmaxwellian over a
!     wedge of velocity space -ft < (v_par/v) < ft.
!     b = (k_per*sqrt(T/m)/(eB/mc))**2
!     The approximation has the form
!     Hn = ft*(h0+a1*h1+a2*h2+a3*h3+a4*h4+a5*h5+a6*h6
!             +a7*h7+a8*h8+a9*h9+a10*h10+a11*h11++a12*h12+a13*h13)
!     where a1 - a13 are fit coefficients which are
!     functions of ft and
!     xa1 = (0.25)**4
!     xa2 = (0.5)**4
!     xa3 = (0.75)**4
!     xa4 = (1.0)**4
!     xa5 = 0.25*(1.5)**5
!     xa6 = 0.25*(2.0)**5
!     xa7 = 0.25*(2.5)**5
!     xa8 = 0.25*(3.0)**5
!     xa9 = 0.25*(4.0)**5
!     xa10 = 0.25*(6.0)**5
!     xa11 = 0.25*(9.0)**5
!     xa12 = 0.25*(15.0)**5
!     xa13 = 0.25*(24.0)**5
!     h0 = 1/(1+b)
!     h1 = b/(xa1+b**2)
!     h2 = b/(xa2+b**2)
!     h3 = b/(xa3+b**2)
!     h4 = b/(xa4+b**2)
!     h5 = b**2/(xa5+b**2.5)
!     h6 = b**2/(xa6+b**2.5)
!     h7 = b**2/(xa7+b**2.5)
!     h8 = b**2/(xa8+b**2.5)
!     h9 = b**2/(xa9+b**2.5)
!     h10 = b**2/(xa10+b**2.5)
!     h11 = b**2/(xa11+b**2.5)
!     h12 = b**2/(xa12+b**2.5)
!     h13 = b**2/(xa13+b**2.5)
!******************************************************************
      IMPLICIT NONE
      INTEGER,PARAMETER :: nf=40,na=13
      INTEGER :: i,j,k
      REAL :: y(na),h(na)
      REAL :: a(nf,na)
      REAL :: g(nf)
      REAL :: b2,b25
      REAL :: gt,dg,ca
      REAL :: h0,hs
! inputs
      REAL,INTENT(IN) :: ft,b
!
      data y &
       / 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, &
         6.0, 9.0, 15.0, 24.0 /
      data g &
       / 0.0, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, &
         0.225, 0.25, 0.275, 0.3, 0.325, 0.35, 0.375, 0.4, &
         0.425, 0.45, 0.475, 0.5, 0.525, 0.55, 0.575, 0.6, &
         0.625, 0.65, 0.675, 0.7, 0.725, 0.75, 0.775, 0.8, &
         0.825, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975, 0.995 /
!
      data a(1:nf,1) &
       / 0,-1.61480787693818E-6,-2.537113820017782E-6, &
        -8.31157970830145E-6,-0.0000116132537675013, &
        -0.00001105333686085922,-0.00001320279577416003, &
        -0.00002208552741141314,-0.00003561033662777171, &
        -0.00004903172912028999,-0.00005909800963140778, &
        -0.00006562977470166618,-0.00007076817724970686, &
        -0.00007735214336528201,-0.0000875719640745126, &
        -0.000102354028160256,-0.0001214168229530987, &
        -0.0001437081678255171,-0.000167932818129024, &
        -0.0001929743882021291,-0.0002181292918781592, &
        -0.0002431594974547346,-0.0002682142224565121, &
        -0.0002936842381292472,-0.0003200427305939506, &
        -0.000347708679393561,-0.000376950548726963, &
        -0.0004078318261086756,-0.0004401981653489624, &
        -0.0004736855405431877,-0.0005077599850560801, &
        -0.0005417340151742959,-0.0005748463324263308, &
        -0.0006062444728986654,-0.0006350644175051447, &
        -0.0006604341617535589,-0.0006815001415549538, &
        -0.000697444547864125,-0.0007075092504000766, &
        -0.0007108565792670304/
!
      data a(1:nf,2) &
       / 0,-0.00004954363511903542, &
        -0.0001379389885596075,-0.0001561890383860288, &
        -0.0002792055098167485,-0.0005416620765113432, &
        -0.0007890483882218802,-0.000927759730441187, &
        -0.00101039166686448,-0.00115341406067672, &
        -0.001438296911230777,-0.00187383543032476, &
        -0.002414000265058838,-0.002996507066376157, &
        -0.003574732210492949,-0.004132142598963165, &
        -0.0046806734381118,-0.005249998894107652, &
        -0.005874633519063111,-0.006583546164506688, &
        -0.007394257479688847,-0.00831125627963979, &
        -0.00932755437786761,-0.0104278718762007, &
        -0.01159218188948583,-0.01279877235534703, &
        -0.01402641282688522,-0.01525560663476376, &
        -0.01646891827435403,-0.01765087962457441, &
        -0.01878735743531674,-0.01986511781887736, &
        -0.02087087902178184,-0.0217912885623709, &
        -0.02261224548429878,-0.02331898349002739, &
        -0.02389610662368788,-0.02432778657782713, &
        -0.02459819163920405,-0.02468786623946138/
!
      data a(1:nf,3) &
       / 0,-0.0003827112301546564, &
        -0.0004688908980583522,-0.002199708150226318, &
        -0.002931744846976895,-0.002123908974616216, &
        -0.002181210474289071,-0.004605423333166559, &
        -0.00863061286159717,-0.01247896120027687, &
        -0.01490803911088084,-0.01581755291376305, &
        -0.01598292633724478,-0.0164523414444748, &
        -0.01803788686089472,-0.02107702937146363, &
        -0.02544729448860753,-0.03072806761086058, &
        -0.03640046263277656,-0.04201020772964428, &
        -0.04726080332521798,-0.05203804633666803, &
        -0.05638399506338718,-0.06044420080935465, &
        -0.06440881490663079,-0.06846159555453574, &
        -0.07274396852879252,-0.07733503396189518, &
        -0.082247477842852,-0.0874321563337234, &
        -0.0927912300671821,-0.098191823665854, &
        -0.1034828835253662,-0.1085053953857707, &
        -0.1131043124106145,-0.1171338311267498, &
        -0.1204607288723746,-0.1229651387401242, &
        -0.124538253005698,-0.1250605883422854/
!
      data a(1:nf,4) &
       / 0,-0.002540372889083718,-0.006685919799796655, &
        -0.00834605758419367,-0.0141986539147988, &
        -0.02557734925105977,-0.03625459630117534, &
        -0.04227288758573232,-0.04546522299937311, &
        -0.05025556459114397,-0.05968923457552303, &
        -0.07386429667506597,-0.0906062894421602, &
        -0.1070154817640277,-0.1207803021598348, &
        -0.1308000231528375,-0.1371567694371842, &
        -0.1407069992863893,-0.1425725471992045, &
        -0.143724987229947,-0.1447486148288815, &
        -0.1457802871667912,-0.1465799352529385, &
        -0.1466703773913419,-0.1454929393471232, &
        -0.1425421896473747,-0.1374608373642871, &
        -0.1300919670260973,-0.1204889118613249, &
        -0.1089006835862105,-0.0957345305934422, &
        -0.081514682275815,-0.06683406995292605, &
        -0.0523193297548779,-0.03859424084886148, &
        -0.02625769163187086,-0.01586716608291115, &
        -0.007928273209733072,-0.002890900945523495, &
        -0.001209815958154486/
!
      data a(1:nf,5) &
       / 0,-0.001792688758988115,-0.002769334624159409, &
        -0.00903072715878006,-0.01224035602468057, &
        -0.0104701528677431,-0.01114062916707894, &
        -0.01888232132085385,-0.03104849726511882, &
        -0.04162522715725507,-0.04620798991681785, &
        -0.04403877610602334,-0.03722985591944204, &
        -0.02885743804932722,-0.02131431397769581, &
        -0.01551683154802667,-0.01093526435562682, &
        -0.006119618855481249,0.0006331807302609925, &
        0.0107225461892253,0.0249150821131089, &
        0.0432498972338039,0.06512425323975039, &
        0.0894829208256717,0.1150410994136554, &
        0.140490090635192,0.1646568038324851, &
        0.1866080807744422,0.2056986862087899, &
        0.2215793258033935,0.2341697500122426, &
        0.2436161409108589,0.2502350705545203, &
        0.2544621919634341,0.2567980838069289, &
        0.257764233383958,0.2578637759587643, &
        0.2575473666618887,0.2571849149164662, &
        0.2570371538312529/
!
      data a(1:nf,6) &
       / 0,-0.002579894294451777,-0.00868836512863523, &
        -0.004151714235597528,-0.00966327723673087, &
        -0.02955908048197177,-0.04344166745514717, &
        -0.03760821193173628,-0.01735374864978274, &
        0.002829990164965488,0.01237312778333665, &
        0.01029521003725834,0.003360774042383296, &
        0.000975020729876641,0.01052580679073009, &
        0.03491100358406264,0.07233107728016348, &
        0.1175796072828099,0.1639080533928105, &
        0.2047533285326446,0.2349451944210703, &
        0.2512992231332235,0.2526825518381504, &
        0.2397215392774724,0.2143266619692717, &
        0.1791760436826848,0.1372512712269478, &
        0.0914728018704398,0.04445382504143827, &
        -0.001643225776121706,-0.0451576704959326, &
        -0.0849210613485667,-0.1201992709515226, &
        -0.1506178886800969,-0.1760692800704513, &
        -0.1966235656113909,-0.2124415460954997, &
        -0.223693089822455,-0.2304864048226563, &
        -0.2326955189416534/
!
      data a(1:nf,7) &
       / 0,-0.001720989358719734,0.002655002069413889, &
        -0.01655636097378585,-0.01571428023441065, &
        0.01620579676303291,0.0386404470553015, &
        0.02419253117004278,-0.01437636752769674, &
        -0.043761910891743,-0.03840241904021477, &
        0.005821170954926033,0.07392011139124511, &
        0.1420611940080034,0.1880419549288359, &
        0.1978934809136671,0.1678815766955462, &
        0.1030663656774289,0.01412979214867896, &
        -0.0860649302446496,-0.1852722570397587, &
        -0.2736551747036851,-0.3445846087654706, &
        -0.394668620597397,-0.4232936643631845, &
        -0.4319340552376819,-0.4234203952235916, &
        -0.4012850752494226,-0.3692475655913782, &
        -0.3308469968631464,-0.2892193543136855, &
        -0.2469855198888891,-0.2062315390250497, &
        -0.1685400368686189,-0.1350659365687666, &
        -0.1066256315408952,-0.0837947127652663, &
        -0.06700534499962984,-0.05662879235820103, &
        -0.05321355209684428/
!
      data a(1:nf,8) &
       / 0,-0.002262588440434705,-0.01133103758999956, &
        0.005512409887285708,0.001785778821827228, &
        -0.03304305504582141,-0.051006152957975, &
        -0.01819593868311817,0.05206573654128066, &
        0.118873182275197,0.1458711885892983, &
        0.1192208487612893,0.04805182121679226, &
        -0.04491012880114669,-0.1344114209439556, &
        -0.2006754049283574,-0.2330004029865321, &
        -0.229693338144006,-0.1958454688123706, &
        -0.140379719653087,-0.07334941596504991, &
        -0.003971923526596438,0.06049924349470718, &
        0.1151849033084929,0.1574853654008633, &
        0.1866956326265865,0.2035043052558006, &
        0.2094942509527074,0.2067157898879144, &
        0.1973556530359046,0.1835110804397204, &
        0.1670506126272552,0.1495510737820277, &
        0.1322808818588104,0.1162244558669098, &
        0.1021233145888178,0.0905287643365276, &
        0.0818581882985255,0.07643987726667087, &
        0.07464723169343828/
!
      data a(1:nf,9) &
       / 0,-0.004663514246367259,-0.003302175033624517, &
        -0.01626271728424343,-0.006927029818041818, &
        0.03583464411147573,0.07557557699775736, &
        0.07751853170360747,0.03797820324277154, &
        -0.02209863375705659,-0.07643705021221824, &
        -0.1073759244208064,-0.1099963284678398, &
        -0.0894827668931036,-0.05602515152837479, &
        -0.02022064342625347,0.0097852838227128, &
        0.02937002559026081,0.03722594065626922, &
        0.03448998286381039,0.02370357547977392, &
        0.00790347926706084,-0.01002688547163994, &
        -0.02773222972442279,-0.04353463945576764, &
        -0.05642664024075585,-0.06597396014680557, &
        -0.07217983187557496,-0.07534615491225473, &
        -0.07595069878721227,-0.07455052918585156, &
        -0.07171196749373892,-0.06796611332581359, &
        -0.0637831603992775,-0.05956303811825382, &
        -0.05563557139134839,-0.05226727969739592, &
        -0.04967202235479381,-0.04801758034904477, &
        -0.04746498730350564/
!
      data a(1:nf,10) &
       / 0,-0.00568634139924752,-0.00822582447647544, &
        0.02301045775197675,0.04757759823427576, &
        0.03043585706859193,-0.0126430864202438, &
        -0.04861377172319868,-0.05943371303779723, &
        -0.04704424085954308,-0.02358667786744889, &
        -0.001295517140360498,0.01253448605457893, &
        0.01608750804979808,0.01130859442886369, &
        0.00175411033141797,-0.00904973399308471, &
        -0.01849313713869238,-0.02513445389912654, &
        -0.0285585300760157,-0.0290644761800976, &
        -0.02733605920068332,-0.02417505669254302, &
        -0.02032410105390975,-0.01637324749690578, &
        -0.01273012603019102,-0.00963099875292301, &
        -0.007172985167553786,-0.005352864892608785, &
        -0.004104470621596094,-0.003329035084428863, &
        -0.002918088771334962,-0.002767887709929662, &
        -0.002788448514355451,-0.002906630629599076, &
        -0.003066553790213365,-0.003227724651528923, &
        -0.003361694537510787,-0.003451104539742693, &
        -0.003481310629500706/
!
      data a(1:nf,11) &
       / 0,-0.0007589724709777822,0.02647893420387985, &
        0.01529645073380253,-0.02104322587836903, &
        -0.03698561695551365,-0.02420375518276133, &
        -0.001293962140512122,0.01408050798582458, &
        0.01624830996718007,0.00863429166725772, &
        -0.002453259101887982,-0.01195690918015689, &
        -0.01744373455883946,-0.01875896303187782, &
        -0.01706311886120541,-0.01392168716257888, &
        -0.01070792123884202,-0.00833297734721454, &
        -0.007213129393777075,-0.007367248459544207, &
        -0.00855929815650642,-0.01043371471683491, &
        -0.01261868016223015,-0.01479150167431653, &
        -0.01671078422542534,-0.0182238580834682, &
        -0.01925897286693473,-0.01980955136825273, &
        -0.01991610023117437,-0.01964914343855488, &
        -0.01909495507087581,-0.01834475909104508, &
        -0.01748725007653101,-0.01660401631152724, &
        -0.01576723420477974,-0.01503889457376441, &
        -0.01447119336692177,-0.01410640590075366, &
        -0.01398407501675658/
!
      data a(1:nf,12) &
       / 0,0.01772082440177327,-0.005399105420125728, &
        -0.02102292777675349,-0.005486020255254523, &
        0.01043558450340497,0.010560795864701, &
        0.000808093398786703,-0.00830688615598142, &
        -0.01176145807470193,-0.01011240588867037, &
        -0.00622438462877651,-0.002717351155977317, &
        -0.000967158754682295,-0.001163484210631195, &
        -0.002782541549688099,-0.005065457752769284, &
        -0.007337100363630739,-0.00914952481545349, &
        -0.01029871500447182,-0.01077310386303381, &
        -0.01068092607925349,-0.01018324730250431, &
        -0.00944519498641531,-0.00860717564879396, &
        -0.007772448715207865,-0.007006522714800609, &
        -0.006342695816834707,-0.00579030724532338, &
        -0.005343302736600071,-0.004987363248153474, &
        -0.004705429401740818,-0.004481132584773598, &
        -0.004301034057383054,-0.004155297951008572, &
        -0.004037737279447664,-0.00394535735828707, &
        -0.003877552443268055,-0.003835600902513023, &
        -0.003821704040399243/
!
      data a(1:nf,13) &
       / 0,-0.00863238583154102,-0.004612900286640256, &
        0.00795792730289846,0.002940046272011631, &
        -0.005753532302313752,-0.007291610791487279, &
        -0.003279625353781592,0.00094606330672321, &
        0.002404592426712804,0.001075361710552384, &
        -0.001649120953695526,-0.004344175923882975, &
        -0.006172355292474261,-0.006927748546083111, &
        -0.006823381803957064,-0.006243544356524337, &
        -0.005561882455387374,-0.005046711180469677, &
        -0.004835792671766108,-0.004954117900239117, &
        -0.00534881080096494,-0.005926412470291415, &
        -0.006583157101910044,-0.007225392895838522, &
        -0.007781085325131736,-0.00820318660756901, &
        -0.00846883567761993,-0.00857538846382309, &
        -0.00853544523410952,-0.00837182208519905, &
        -0.00811318942239733,-0.007790571492123561, &
        -0.007434835050358763,-0.007075017623380862, &
        -0.006737402894926134,-0.006445057721582979, &
        -0.006217752770321989,-0.006071927902178915, &
        -0.006023036419830553/
!
      b2 = b*b
      b25 = b**2.5
      h0 = 0.0
      do i=1,4
       h(i) = b/( y(i)**4 + b2)
      enddo
      do i=5,na
       h(i) = b2/(0.25*y(i)**5 + b25)       
      enddo
! transform to gt grid
      gt = SQRT(1-ft)
      if(gt.gt.g(nf))gt = g(nf)
! find grid position
      if(gt.le.g(2))then
        i=1
      else
        i = INT(gt/0.025)
      endif
      j = i+1
! interpolate the coefficients
      dg = (gt - g(i))/(g(j)-g(i))
      hs = h0
      do k=1,na
        ca = a(i,k) + (a(j,k)-a(i,k))*dg
! sum up the terms
        hs = hs + ca*h(k)
      enddo
! final answer
      FLR_dHr11 = (3.0*ft**5)*hs
! 
      END FUNCTION FLR_dHr11
!
!
      REAL FUNCTION FLR_dHr13(ft,b)
!******************************************************************
!     Approxmation to the integral of (J0^2)*Fmaxwellian over a
!     wedge of velocity space -ft < (v_par/v) < ft.
!     b = (k_per*sqrt(T/m)/(eB/mc))**2
!     The approximation has the form
!     Hn = ft*(h0+a1*h1+a2*h2+a3*h3+a4*h4+a5*h5+a6*h6
!             +a7*h7+a8*h8+a9*h9+a10*h10+a11*h11++a12*h12+a13*h13)
!     where a1 - a13 are fit coefficients which are
!     functions of ft and
!     xa1 = (0.25)**4
!     xa2 = (0.5)**4
!     xa3 = (0.75)**4
!     xa4 = (1.0)**4
!     xa5 = 0.25*(1.5)**5
!     xa6 = 0.25*(2.0)**5
!     xa7 = 0.25*(2.5)**5
!     xa8 = 0.25*(3.0)**5
!     xa9 = 0.25*(4.0)**5
!     xa10 = 0.25*(6.0)**5
!     xa11 = 0.25*(9.0)**5
!     xa12 = 0.25*(15.0)**5
!     xa13 = 0.25*(24.0)**5
!     h0 = 1/(1+b)
!     h1 = b/(xa1+b**2)
!     h2 = b/(xa2+b**2)
!     h3 = b/(xa3+b**2)
!     h4 = b/(xa4+b**2)
!     h5 = b**2/(xa5+b**2.5)
!     h6 = b**2/(xa6+b**2.5)
!     h7 = b**2/(xa7+b**2.5)
!     h8 = b**2/(xa8+b**2.5)
!     h9 = b**2/(xa9+b**2.5)
!     h10 = b**2/(xa10+b**2.5)
!     h11 = b**2/(xa11+b**2.5)
!     h12 = b**2/(xa12+b**2.5)
!     h13 = b**2/(xa13+b**2.5)
!******************************************************************
      IMPLICIT NONE
      INTEGER,PARAMETER :: nf=40,na=13
      INTEGER :: i,j,k
      REAL :: y(na),h(na)
      REAL :: a(nf,na)
      REAL :: g(nf)
      REAL :: b2,b25
      REAL :: gt,dg,ca
      REAL :: h0,hs
! inputs
      REAL,INTENT(IN) :: ft,b
!
      data y &
       / 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, &
         6.0, 9.0, 15.0, 24.0 /
      data g &
       / 0.0, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, &
         0.225, 0.25, 0.275, 0.3, 0.325, 0.35, 0.375, 0.4, &
         0.425, 0.45, 0.475, 0.5, 0.525, 0.55, 0.575, 0.6, &
         0.625, 0.65, 0.675, 0.7, 0.725, 0.75, 0.775, 0.8, &
         0.825, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975, 0.995 /
!
      data a(1:nf,1) &
       /-0.0001720852901934805,-0.0001732210791781909, &
        -0.0001745908453220229,-0.0001785611038002381, &
        -0.0001806824667585374,-0.000181436420050485, &
        -0.0001848940191357485,-0.0001925445226894828, &
        -0.0002023545910054129,-0.0002114153882053508, &
        -0.0002182577778816039,-0.0002233709461574561, &
        -0.0002284335384018101,-0.0002352356710411683, &
        -0.0002449286297895609,-0.0002577811318946121, &
        -0.0002733333996055133,-0.0002907387962370997, &
        -0.0003091117281991792,-0.0003277748361715801, &
        -0.0003463714608366705,-0.0003648626200444751, &
        -0.0003834481175709106,-0.0004024505830790349, &
        -0.0004222216122482294,-0.0004430376930850746, &
        -0.000465046025643263,-0.000488241670004778, &
        -0.0005124655902783491,-0.0005374153098267218, &
        -0.000562679442471772,-0.0005877444866355384, &
        -0.0006120626413953809,-0.0006350225504770473, &
        -0.0006560175941168399,-0.0006744394266359633, &
        -0.0006896955507491604,-0.0007012186771003332, &
        -0.0007084829291577841,-0.0007108987177533797/
!
      data a(1:nf,2) &
       /-0.00607248941980851,-0.006121842246581287, &
        -0.006185104721909846,-0.006225664343816172, &
        -0.006349950844885546,-0.006547697707938205, &
        -0.006723329264828326,-0.006842913158970105, &
        -0.006956756613763591,-0.007136384736828071, &
        -0.007419326924405844,-0.007796701511558622, &
        -0.00823148091121215,-0.00868409957162872, &
        -0.00913021255405388,-0.00956634910570389, &
        -0.01000605693840316,-0.01047156763679631, &
        -0.01098530708173541,-0.011563803939349, &
        -0.01221481442801715,-0.01293720010341187, &
        -0.01372262569351523,-0.01455808113066192, &
        -0.0154282683090352,-0.01631755736516912, &
        -0.01721148182993951,-0.01809720875394026, &
        -0.01896360811806563,-0.01980106818821257, &
        -0.02060092226024048,-0.02135514484741672, &
        -0.02205555956851759,-0.02269394240300271, &
        -0.0232614286708232,-0.02374862367189322, &
        -0.02414561013298049,-0.02444206650734732, &
        -0.02462757801053459,-0.02468905824616768/
!
      data a(1:nf,3) &
       /-0.03610514817180309,-0.03631989906846428, &
        -0.03656905631329484,-0.03765525106956956, &
        -0.03796962377259699,-0.03766628355445286, &
        -0.03823686590133487,-0.04024004226278083, &
        -0.0429241388784395,-0.04519574992804575, &
        -0.04648677346161791,-0.04695960329474374, &
        -0.0472276494678239,-0.04795303749841143, &
        -0.04956180542031306,-0.05214754200198028, &
        -0.05552580933831855,-0.0593617336712144, &
        -0.06330244880460282,-0.06707313127839791, &
        -0.07052267517214655,-0.07362538680258994, &
        -0.07645321613722542,-0.07913473749285917, &
        -0.0818142441025618,-0.0846211503089763, &
        -0.0876445767547864,-0.0909261040406512, &
        -0.0944596334911643,-0.098196620610693, &
        -0.1020569655744038,-0.1059388653512883, &
        -0.1097312616454955,-0.1133199231231983, &
        -0.1165959024251741,-0.119458274569542, &
        -0.1218158799104798,-0.1235873255659947, &
        -0.1246984475806911,-0.125067198374595/
!
      data a(1:nf,4) &
       /-0.1385901922759657,-0.139766145728718, &
        -0.1412683267711806,-0.1414634122318736, &
        -0.1444149642405343,-0.1496983163145244, &
        -0.153447924999075,-0.1541622173165233, &
        -0.1536878747736368,-0.15474703211027, &
        -0.1587108907896043,-0.1650658316430402, &
        -0.1721263166779406,-0.1780687405191892, &
        -0.1816679144253255,-0.1825510271287887, &
        -0.1810635556519464,-0.1779453603383395, &
        -0.1739926754415925,-0.1698127593489506, &
        -0.1657077694335004,-0.1616721210528721, &
        -0.1574660730854656,-0.1527235764503118, &
        -0.147060486722944,-0.1401546916002475, &
        -0.1318129952494658,-0.1219900157874474, &
        -0.1107878642103149,-0.0984405883756617, &
        -0.0852838212431927,-0.07172520321363982, &
        -0.05820976016898415,-0.04519825231321475, &
        -0.03314250183571823,-0.02247308784515265, &
        -0.01359035967359817,-0.006859649833472442, &
        -0.002612103937521253,-0.001198252815861577/
!
      data a(1:nf,5) &
       /-0.00892692407679801,-0.00880823634961514, &
        -0.00856663905689591,-0.01045811880961684, &
        -0.00940553237565339,-0.005760940689889478, &
        -0.004143204726761469,-0.006251669434260353, &
        -0.00957242356828436,-0.01042494027457343, &
        -0.006765279524783684,0.001097854366635039, &
        0.01135813722401534,0.02199366305433781, &
        0.03169510529196419,0.04019279321566679, &
        0.04809025118497229,0.05644795662269263, &
        0.06634080916306947,0.07852954393492596, &
        0.0932980601106741,0.1104422785161131, &
        0.1293652407629442,0.1492251584722129, &
        0.1690926185857799,0.1880756938476721, &
        0.2054300831967992,0.2206072041970011, &
        0.2332751821155011,0.2433134875727382, &
        0.2507853884673869,0.2559033594280402, &
        0.2589863153551677,0.260423903058823, &
        0.2606382515182517,0.2600549351187595, &
        0.2590771943601224,0.2580637462762444, &
        0.2573114307293284,0.2570421479913989/
!
      data a(1:nf,6) &
       /0.03804841377048915,0.0379189520509829, &
        0.03768722760561733,0.04336948074526167, &
        0.04228822136888659,0.03565387607441993, &
        0.03602639992988004,0.04861914381882945, &
        0.06755415340373426,0.0836340816873745, &
        0.0917838969490696,0.0931238032369501, &
        0.0928067431954716,0.096556161138451, &
        0.108016102629155,0.1276756722544397, &
        0.1531751322983258,0.1803962931720272, &
        0.2047407962358801,0.222188216900156, &
        0.2299484691503001,0.2266998732040241, &
        0.2125022863290836,0.1885104871776007, &
        0.1566001888027653,0.1190264621464491, &
        0.07808278072603833,0.03589724612113653, &
        -0.005712742060272189,-0.04531148786097683, &
        -0.0818546042849582,-0.1146627601617537, &
        -0.1433650853353979,-0.1678378142981397, &
        -0.1881312172354966,-0.2044034836843851, &
        -0.2168569883688929,-0.2256788844296539, &
        -0.2309910015827104,-0.2327159213488404/
!
      data a(1:nf,7) &
       /-0.03692053533845384,-0.03625229233735828, &
        -0.03500148450632328,-0.0449668368428443, &
        -0.03970636576829314,-0.0208644296404824, &
        -0.01396160148453779,-0.02968318958793911, &
        -0.05525507316342408,-0.07034008553818321, &
        -0.06323719092091001,-0.03626147084551928, &
        -0.001915273790302761,0.02436740197791111, &
        0.0304035935073635,0.01045739886177917, &
        -0.03448596802822405,-0.0983552717574589, &
        -0.1723827939436131,-0.2473715920050457, &
        -0.3153613396086348,-0.3705583416168511, &
        -0.4096139977286206,-0.4314376954806513, &
        -0.4367343928673433,-0.4275028172951602, &
        -0.4064371290767913,-0.3765266002647731, &
        -0.3407143517364498,-0.3016753809206991, &
        -0.2616966249879707,-0.2226297179055868, &
        -0.1859048614913397,-0.1525724134123183, &
        -0.1233732126028972,-0.0988125977913461, &
        -0.07923811007844654,-0.06491563715733523, &
        -0.05609142376348755,-0.05319176145033045/
!
      data a(1:nf,8) &
       /-0.00875677401339384,-0.00975595113007712, &
        -0.01114777054970348,0.0003473822110466612, &
        -0.003918827228190719,-0.02170828237349894, &
        -0.02370281041802966,0.002943045418759251, &
        0.04332487366172104,0.07166450568737245, &
        0.07014411649664129,0.03663121694941117, &
        -0.01794874706898569,-0.07681692396596423, &
        -0.1242438809641733,-0.149790751387856, &
        -0.1496005978417584,-0.1255426810287697, &
        -0.083287210253144,-0.03020422765607875, &
        0.02636588935324343,0.0802503340649529, &
        0.1269961861239374,0.1639629209508007, &
        0.19010523575555,0.2056693392345858, &
        0.2117253462721215,0.2098537244016842, &
        0.2018415356466485,0.1894666691438577, &
        0.1743602249470348,0.1579287115236018, &
        0.1413289337181229,0.1254703989523378, &
        0.1110455234248186,0.098567500950907, &
        0.0884145859781989,0.0808759171705099, &
        0.07618535928370845,0.07463693366785633/
!
      data a(1:nf,9) &
       /-0.0196115031539788,-0.01877247063399729, &
        -0.01670289952690608,-0.02137842176950579, &
        -0.01182932239728334,0.01090883014914241, &
        0.0233201482975775,0.0100466790537178, &
        -0.02431553568473666,-0.06345739842201368, &
        -0.0917097668138074,-0.1010509091687145, &
        -0.0918268433439754,-0.0699189151001772, &
        -0.04313381918904742,-0.01842431095306482, &
        -0.0004365818365540958,0.00878763514917603, &
        0.00937977778068444,0.002904198934665961, &
        -0.00839036927748457,-0.02216507270911849, &
        -0.03637773979221481,-0.04947678036168672, &
        -0.06043978576449562,-0.0687605285008166, &
        -0.07430460643575127,-0.07723102770742746, &
        -0.07788087652575094,-0.07669045883020041, &
        -0.07412620286506364,-0.07063950963392594, &
        -0.06664029499133834,-0.06248303642334688, &
        -0.05846424058079481,-0.05482567237063669, &
        -0.05176173897160026,-0.04942946052287771, &
        -0.04795354480608085,-0.04746238092042564/
!
      data a(1:nf,10) &
       /-0.01378907486603525,-0.01358715656656082, &
        -0.01070614754682926,0.006320194546080384, &
        0.01056455011446667,-0.00982744322023507, &
        -0.03837997616599219,-0.05540388526553331, &
        -0.05413324686595755,-0.03969696910126052, &
        -0.02144950574376002,-0.006947415198097316, &
        0.0002915035702351631,0.0003437633953838147, &
        -0.004631007426789141,-0.01187127413063315, &
        -0.01899785169033723,-0.02445362196688425, &
        -0.02754178965868374,-0.02824477653585585, &
        -0.02697091211501547,-0.02432000434665262, &
        -0.02090933859971328,-0.01726807391569666, &
        -0.01379827545711129,-0.0107345650358524, &
        -0.00821096573152233,-0.006257997140070079, &
        -0.004840730867236192,-0.003886262547089658, &
        -0.003304851118438372,-0.003005353496246069, &
        -0.002904207022886703,-0.002930974366737139, &
        -0.003029442340981522,-0.00315726000469797, &
        -0.003284085406968878,-0.00338854891090074, &
        -0.003458156626376737,-0.003481585360045994/
!
      data a(1:nf,11) &
       /-0.01097253337919979,-0.006172154696782852, &
        0.003898751178140989,-0.01185106226478099, &
        -0.03337873743213152,-0.03622613229751872, &
        -0.02256147890871899,-0.006825237382246531, &
        0.001075879534249901,-0.0001832951384483741, &
        -0.006880473500531914,-0.01450719812840289, &
        -0.02005464699330575,-0.02245793407424188, &
        -0.02208329587871056,-0.01998431044231326, &
        -0.01730356614628235,-0.01492930710085854, &
        -0.01337587035405614,-0.0128090867809197, &
        -0.01314061252573559,-0.01413686867854607, &
        -0.01551246164417715,-0.01699601073067774, &
        -0.018356944893546,-0.0194634249623887, &
        -0.02022017414805451,-0.02059890229310948, &
        -0.02061522749840841,-0.02031480805519981, &
        -0.01976102115710788,-0.01902508853973488, &
        -0.01817895891483604,-0.01729057443005808, &
        -0.01642126509213705,-0.01562466334716363, &
        -0.01494656118747039,-0.0144258498752461, &
        -0.01409427464320224,-0.01398357257390258/
!
      data a(1:nf,12) &
       /-0.007057952931431471,0.00188180644779207, &
        -0.01699593740851121,-0.02093760182610706, &
        -0.007630792102594064,0.001078818368279755, &
        -0.001274006898779019,-0.00840531859105285, &
        -0.01356307780528856,-0.01448852791706003, &
        -0.0123547633976167,-0.00931915048785103, &
        -0.007007905894443318,-0.00608903837754338, &
        -0.00646626249037237,-0.007655087578056301, &
        -0.00910198051678586,-0.01037081775195304, &
        -0.01120867549977393,-0.01153359743147452, &
        -0.01138621009814308,-0.01087488941296659, &
        -0.0101300372644464,-0.00927340368030586, &
        -0.00841686833479723,-0.007599094861185431, &
        -0.006871638216755282,-0.006250208109854893, &
        -0.005734486891982527,-0.005314552987916608, &
        -0.00497589564542442,-0.004703239393670884, &
        -0.004482682955137855,-0.004303169598002921, &
        -0.004156694645857662,-0.004038189605381585, &
        -0.003945308408091114,-0.003877374796793711, &
        -0.003835527538492145,-0.003821694312494478/
!
      data a(1:nf,13) &
       /-0.005670146637748142,-0.01398120796873224, &
        -0.00676642424165872,-0.0002236440852413532, &
        -0.005149345648505287,-0.01013995909406717, &
        -0.00983497897604507,-0.006634065815399415, &
        -0.004088742793615588,-0.003663698889886223, &
        -0.004933611052919992,-0.006811193179665199, &
        -0.00839152949402182,-0.00924434269509694, &
        -0.00935727980956315,-0.00895542343024548, &
        -0.00833078375020208,-0.007732201922198258, &
        -0.007317473576082806,-0.007149838191399584, &
        -0.007218479939464472,-0.007466540448954832, &
        -0.007817497552333055,-0.00819492518987236, &
        -0.00853327020648865,-0.00878948420670655, &
        -0.00893562772551707,-0.00896233044654812, &
        -0.00887396696179394,-0.00868476998137084, &
        -0.00841517959618544,-0.00808873655812051, &
        -0.007729786737521298,-0.007361876006765772, &
        -0.007006734854415985,-0.006683891373024692, &
        -0.006410256291920601,-0.006200609889542852, &
        -0.006067333609944024,-0.00602284666988543/
!
      b2 = b*b
      b25 = b**2.5
      h0 = 0.0
      do i=1,4
       h(i) = b/( y(i)**4 + b2)
      enddo
      do i=5,na
       h(i) = b2/(0.25*y(i)**5 + b25)       
      enddo
! transform to gt grid
      gt = SQRT(1-ft)
      if(gt.gt.g(nf))gt = g(nf)
! find grid position
      if(gt.le.g(2))then
        i=1
      else
        i = INT(gt/0.025)
      endif
      j = i+1
! interpolate the coefficients
      dg = (gt - g(i))/(g(j)-g(i))
      hs = h0
      do k=1,na
        ca = a(i,k) + (a(j,k)-a(i,k))*dg
! sum up the terms
        hs = hs + ca*h(k)
      enddo
! final answer
      FLR_dHr13 = (5.0/3.0)*(ft**3)*hs
! 
      END FUNCTION FLR_dHr13
!
!
      REAL FUNCTION FLR_dHr33(ft,b)
!******************************************************************
!     Approxmation to the integral of (J0^2)*Fmaxwellian over a
!     wedge of velocity space -ft < (v_par/v) < ft.
!     b = (k_per*sqrt(T/m)/(eB/mc))**2
!     The approximation has the form
!     Hn = ft*(h0+a1*h1+a2*h2+a3*h3+a4*h4+a5*h5+a6*h6
!             +a7*h7+a8*h8+a9*h9+a10*h10+a11*h11++a12*h12+a13*h13)
!     where a1 - a13 are fit coefficients which are
!     functions of ft and
!     xa1 = (0.25)**4
!     xa2 = (0.5)**4
!     xa3 = (0.75)**4
!     xa4 = (1.0)**4
!     xa5 = 0.25*(1.5)**5
!     xa6 = 0.25*(2.0)**5
!     xa7 = 0.25*(2.5)**5
!     xa8 = 0.25*(3.0)**5
!     xa9 = 0.25*(4.0)**5
!     xa10 = 0.25*(6.0)**5
!     xa11 = 0.25*(9.0)**5
!     xa12 = 0.25*(15.0)**5
!     xa13 = 0.25*(24.0)**5
!     h0 = 1/(1+b)
!     h1 = b/(xa1+b**2)
!     h2 = b/(xa2+b**2)
!     h3 = b/(xa3+b**2)
!     h4 = b/(xa4+b**2)
!     h5 = b**2/(xa5+b**2.5)
!     h6 = b**2/(xa6+b**2.5)
!     h7 = b**2/(xa7+b**2.5)
!     h8 = b**2/(xa8+b**2.5)
!     h9 = b**2/(xa9+b**2.5)
!     h10 = b**2/(xa10+b**2.5)
!     h11 = b**2/(xa11+b**2.5)
!     h12 = b**2/(xa12+b**2.5)
!     h13 = b**2/(xa13+b**2.5)
!******************************************************************
      IMPLICIT NONE
      INTEGER,PARAMETER :: nf=40,na=13
      INTEGER :: i,j,k
      REAL :: y(na),h(na)
      REAL :: a(nf,na)
      REAL :: g(nf)
      REAL :: b2,b25
      REAL :: gt,dg,ca
      REAL :: h0,hs
! inputs
      REAL,INTENT(IN) :: ft,b
!
      data y &
       / 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, &
         6.0, 9.0, 15.0, 24.0 /
      data g &
       / 0.0, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, &
         0.225, 0.25, 0.275, 0.3, 0.325, 0.35, 0.375, 0.4, &
         0.425, 0.45, 0.475, 0.5, 0.525, 0.55, 0.575, 0.6, &
         0.625, 0.65, 0.675, 0.7, 0.725, 0.75, 0.775, 0.8, &
         0.825, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975, 0.995 /
!
      data a(1:nf,1) &
       / -0.0003775984382690199,-0.0003784909214937589, &
        -0.0003795916973605756,-0.0003818149907637014, &
        -0.0003836784414571823,-0.0003853559443972812, &
        -0.0003882097321739775,-0.0003927219118005234, &
        -0.0003981953250677991,-0.0004036559734204095, &
        -0.0004086316467407606,-0.0004133151223572498, &
        -0.0004182920476007317,-0.0004241690695009534, &
        -0.0004313196185939467,-0.0004398069924043014, &
        -0.0004494443915202362,-0.0004599182822971585, &
        -0.0004709112101420126,-0.000482186907640068, &
        -0.0004936268174467098,-0.000505224943403366, &
        -0.0005170560369479108,-0.0005292328610916521, &
        -0.0005418647927014806,-0.0005550251939125082, &
        -0.0005687305322148416,-0.0005829310005787, &
        -0.0005975105154485622,-0.0006122932698482241, &
        -0.0006270540292924884,-0.0006415298623842052, &
        -0.0006554317658925765,-0.000668454969992573, &
        -0.0006802877506793563,-0.0006906185388511559, &
        -0.0006991404526695681,-0.0007055626243545941, &
        -0.0007096020520976943,-0.0007109444669523535/
!
      data a(1:nf,2) &
       /-0.01348913689051852,-0.01352414027264805, &
        -0.013568412169703,-0.01361450586444158, &
        -0.01369794607297215,-0.01381507377093354, &
        -0.01393400027506674,-0.01404385782540596, &
        -0.01416198192962586,-0.01431241890062295, &
        -0.0145073474359888,-0.01474315937983791, &
        -0.01500691177657436,-0.01528514572443071, &
        -0.01556988621174899,-0.01586040271585612, &
        -0.01616169641878973,-0.01648147318641083, &
        -0.01682712761247745,-0.01720362545751941, &
        -0.01761254486734975,-0.01805211345268071, &
        -0.01851788575349568,-0.01900368923618344, &
        -0.0195025502483146,-0.02000742644496875, &
        -0.02051167770574439,-0.02100928372259071, &
        -0.02149485983617412,-0.02196353867929394, &
        -0.02241078392656314,-0.02283218921383102, &
        -0.02322329965127276,-0.02357947980071302, &
        -0.02389583506594167,-0.02416718912992774, &
        -0.02438808912498339,-0.02455297332675244, &
        -0.02465606301267123,-0.02469022149849027/
!
      data a(1:nf,3) &
       /-0.0714340485045578,-0.07159391471857771, &
        -0.07178774647878488,-0.07230377377456472, &
        -0.07260556088199886,-0.07274883985689242, &
        -0.07323212216066693,-0.07423559683903874, &
        -0.07550062516477435,-0.07666010361627115, &
        -0.07753018946529523,-0.07817517028720911, &
        -0.07880842950827824,-0.07965362730748939, &
        -0.080848191208916,-0.0824128300241397, &
        -0.0842730077940696,-0.0863051801643067, &
        -0.0883837005355458,-0.0904140499421565, &
        -0.09234788526751,-0.094182203711143, &
        -0.0959481408131762,-0.0976953357321773, &
        -0.099476577575556,-0.1013356538753044, &
        -0.1032996266722729,-0.105375499291938, &
        -0.1075504816204462,-0.1097947482041486, &
        -0.1120655578751981,-0.1143117986392525, &
        -0.1164782298361088,-0.1185089901968578, &
        -0.1203501442010024,-0.1219512124567196, &
        -0.1232660285045255,-0.1242514975251616, &
        -0.1248692905671082,-0.1250741770435721/
!
      data a(1:nf,4) &
       /-0.1140330903982963,-0.1143627500803044, &
        -0.1147815955364487,-0.1147294686719825, &
        -0.1155559248341305,-0.1171014492244224, &
        -0.1180659488907664,-0.1179547633598834, &
        -0.1173933622129277,-0.1172839824655971, &
        -0.1180566125492014,-0.1195002340034679, &
        -0.1210157589042285,-0.1219730812309402, &
        -0.1219609734839107,-0.1208679190240416, &
        -0.1188291865204107,-0.1161096622214874, &
        -0.1129843925065605,-0.1096539578833283, &
        -0.1062065331552173,-0.1026209030026535, &
        -0.0987962906320349,-0.0945936696334355, &
        -0.0898762966499645,-0.0845417960637506, &
        -0.07854251306609659,-0.07189414530143539, &
        -0.06467466153932774,-0.05701637111979519, &
        -0.04909411044501502,-0.04111202914525327, &
        -0.03329093625483218,-0.0258573764382468, &
        -0.01903508597409919,-0.01303898541577865, &
        -0.00807090524935842,-0.0043200857632364, &
        -0.001957147187660228,-0.001171474245520463/
!
      data a(1:nf,5) &
       /0.1120065210674691,0.1123493252008408, &
        0.1128098021524352,0.1127187275372962, &
        0.1137727412093214,0.1158452410425273, &
        0.1173973469746181,0.1178868235314183, &
        0.1181822452951641,0.1195227802566951, &
        0.1225739157733163,0.1271991179308292, &
        0.1327629695771262,0.1385733262946243, &
        0.1441954616069611,0.1495555444773096, &
        0.1548737819133407,0.1605131280045995, &
        0.1668223486417832,0.1740223071688268, &
        0.1821527045078547,0.1910737286092492, &
        0.2005055052813683,0.2100858153141898, &
        0.2194297326224316,0.2281802828155573, &
        0.2360447002647248,0.2428151751174638, &
        0.2483757821664405,0.2526986764479607, &
        0.2558330590715641,0.2578900444656525, &
        0.2590260625590885,0.2594265357636098, &
        0.2592909513438653,0.2588198279312485, &
        0.2582030188958099,0.2576121666700541, &
        0.2571870419908934,0.257036957597013/
!
      data a(1:nf,6) &
       /0.006383272857812912,0.006260261985697357, &
        0.006085217459442507,0.007830752487486956, &
        0.007266174862488754,0.004824725841576648, &
        0.004695022601964442,0.00855056958047666, &
        0.01434287666604675,0.01895499536150693, &
        0.02070741707592249,0.02002004494687525, &
        0.01864044962183274,0.01844653897221966, &
        0.02054762362803171,0.02494008876019121, &
        0.03064333701665674,0.03610389710436661, &
        0.03965685109507544,0.03990135329280295, &
        0.03592652316111244,0.02738595153406965, &
        0.01445460912590438,-0.002285772530618868, &
        -0.0219900719476506,-0.04368990430968335, &
        -0.06641103042118135,-0.0892594127506793, &
        -0.1114751112829822,-0.1324578786696552, &
        -0.1517706125179691,-0.1691270285605245, &
        -0.1843695887668296,-0.1974422075073752, &
        -0.2083611532495079,-0.2171863628244425, &
        -0.2239932264935012,-0.2288511365928624, &
        -0.2317908700470113,-0.232748385330183/
!
      data a(1:nf,7) &
       /-0.1399180612585296,-0.1399527770192664, &
        -0.1398659147279096,-0.1436352417354466, &
        -0.142448788096095,-0.1369359121423814, &
        -0.1356379836134975,-0.1420250000237716, &
        -0.1517088542682499,-0.1578805294817769, &
        -0.156771488648278,-0.1493882045783507, &
        -0.1401189896476811,-0.1342026561637165, &
        -0.1356065053131723,-0.1459634450292836, &
        -0.1645327401079661,-0.1888226513567488, &
        -0.2154656370007496,-0.2410367009082698, &
        -0.2626472807163433,-0.2782688427159729, &
        -0.286818802037827,-0.2880773884656281, &
        -0.2825090326426683,-0.2710496541754656, &
        -0.2549025441981624,-0.2353673240486916, &
        -0.2137119702073841,-0.1910881760754131, &
        -0.1684844591953555,-0.1467093821352162, &
        -0.1263962950880066,-0.1080225261065884, &
        -0.0919370251786104,-0.07839188337252965, &
        -0.0675765040034599,-0.05964254399746428, &
        -0.05474747429519171,-0.05313701974272243/
!
      data a(1:nf,8) &
       /0.05494146246597086,0.05476757208751299, &
        0.05451105698340681,0.05861486430138322, &
        0.0575266314524935,0.05207229762679846, &
        0.05204456303275735,0.06159124608888851, &
        0.07555550481366912,0.0852716240791904, &
        0.0849386357110087,0.07415555446470105, &
        0.05692724029902103,0.03913105778741518, &
        0.02612569784087279,0.02130898651549851, &
        0.0257135518833241,0.03835237179467515, &
        0.05693178375234079,0.07861592238552704, &
        0.1006499060714908,0.1207649119687487, &
        0.1373691066368469,0.1495710432258705, &
        0.1570948685498103,0.1601414902576158, &
        0.1592364778949157,0.1550907257673508, &
        0.1484872654016621,0.1401985597832289, &
        0.1309326060229386,0.1213036099385823, &
        0.1118212605551676,0.1028934329337832, &
        0.0948375607200834,0.0878967385117478, &
        0.0822596020206743,0.07807196843530644, &
        0.07546782548441131,0.07460766649563456/
!
      data a(1:nf,9) &
       /-0.0398451488497834,-0.03961591188359483, &
        -0.0389957040923425,-0.04064464403843992, &
        -0.03758406651666704,-0.03023979849759407, &
        -0.02649324119765648,-0.0313745208395301, &
        -0.04317034507169728,-0.05628528536114965, &
        -0.06545198920626315,-0.06808339151985381, &
        -0.06447137303911115,-0.05677541120425746, &
        -0.04775907426661572,-0.03982754336960694, &
        -0.03453430996322867,-0.03249376020584908, &
        -0.03354844409205691,-0.03704223049115427, &
        -0.04209305312288392,-0.04780786665049304, &
        -0.05342130029848845,-0.05836373757996456, &
        -0.06227573084101287,-0.0649880463552563, &
        -0.0664840907258233,-0.06685705714556745, &
        -0.06626960597163369,-0.06492025854182747, &
        -0.06301786976990095,-0.06076417384515747, &
        -0.05834307420144191,-0.05591542286841405, &
        -0.0536178133904425,-0.05156392890448625, &
        -0.04984837822568631,-0.04854663733176499, &
        -0.04772601854970666,-0.04745307775023967/
!
      data a(1:nf,10) &
       /-0.01161133931023969,-0.01153938064611452, &
        -0.01057285820012358,-0.0049223271202109, &
        -0.003588251155853147,-0.01044346905972372, &
        -0.01988273674147974,-0.02531049461919185, &
        -0.02454392106385011,-0.01941308233312685, &
        -0.01313083616566695,-0.00824837513590935, &
        -0.005908378770655955,-0.006017293456034769, &
        -0.007778602790189302,-0.01020957040362717, &
        -0.01248138503932321,-0.01406757306810569, &
        -0.01475167403820347,-0.01455860473271645, &
        -0.01366164787490814,-0.01229724284499061, &
        -0.01070221285236652,-0.00907602816744122, &
        -0.007564294252473362,-0.006257234909331375, &
        -0.005197035277564422,-0.004389137963683196, &
        -0.00381415561483539,-0.003438390482107981, &
        -0.003222182438620295,-0.003125716389222503, &
        -0.003112767740640059,-0.003152621772492611, &
        -0.003220646391017606,-0.003298095798317274, &
        -0.003370788802488007,-0.003430389378159032, &
        -0.00346900263127742,-0.00348203091077785/
!
      data a(1:nf,11) &
       /-0.01475488097845215,-0.0131659579956464, &
        -0.00984530987779518,-0.01513735127789423, &
        -0.02228592493738319,-0.02314160299699469, &
        -0.0185117836593792,-0.01328166069492467, &
        -0.01075703158102055,-0.01132667559915125, &
        -0.01368747563581602,-0.01629567312440655, &
        -0.01814047602813843,-0.01888579560729661, &
        -0.01868553199253987,-0.01792236671361118, &
        -0.01699950437264116,-0.016222625939717, &
        -0.015760631845118,-0.01565717215672738, &
        -0.01586605904287814,-0.01629138577017035, &
        -0.01682158948314361,-0.01735331218790162, &
        -0.01780494562670211,-0.01812183749868778, &
        -0.0182757449959151,-0.01826096079271996, &
        -0.01808897996091699,-0.017783067225668, &
        -0.01737333768796418,-0.01689296241701561, &
        -0.01637525021394273,-0.0158517700075011, &
        -0.01535132537740526,-0.01489943715840474, &
        -0.01451859438722707,-0.01422748225800201, &
        -0.01404300000371334,-0.01398147400601609/
!
      data a(1:nf,12) &
       /-0.006838632233982654,-0.003867366377080605, &
        -0.0101582215475142,-0.01142922646116672, &
        -0.006964029077651357,-0.004081630393476146, &
        -0.004906577329531932,-0.007295796077075815, &
        -0.00897884709760874,-0.00922246511907712, &
        -0.00844514582875874,-0.007387802617258163, &
        -0.006598480752564865,-0.006291829842488817, &
        -0.006420444298881697,-0.00680709611944442, &
        -0.007255794018551719,-0.007615871115776418, &
        -0.007803641218360879,-0.007796807511151531, &
        -0.007616407455299416,-0.007306690784178742, &
        -0.006918592206117059,-0.00649867652769282, &
        -0.006083394283866631,-0.005697384311009696, &
        -0.005354375251937292,-0.0050594072678861, &
        -0.004811483692306592,-0.004605969343503996, &
        -0.004436684642927429,-0.00429716652445709, &
        -0.004181701046353341,-0.00408579512581202, &
        -0.004006299773789828,-0.003941323661963146, &
        -0.003889911981985543,-0.003852401291763883, &
        -0.003829117182666242,-0.003821432865102725/
!
      data a(1:nf,13) &
       /-0.006920509825410892,-0.00969183149312854, &
        -0.007278276009691722,-0.005113151110433332, &
        -0.006775129741594955,-0.00843624714694329, &
        -0.00831891531707461,-0.007247331121387556, &
        -0.006418359905497736,-0.006309510672439072, &
        -0.006765199413237366,-0.007410034624458625, &
        -0.007937363920071598,-0.00820651717389971, &
        -0.00822070603409865,-0.00806318954690526, &
        -0.007837694577168687,-0.007630355200276639, &
        -0.007493620414041794,-0.007445765449390644, &
        -0.007478710682293188,-0.007568519013239796, &
        -0.007684986373830725,-0.007798876915496117, &
        -0.007886351130596064,-0.007930997146500806, &
        -0.007924093198801465,-0.007863766031389614, &
        -0.007753580296661654,-0.007601017859814874, &
        -0.007415949571695458,-0.007209537896678934, &
        -0.006993230917857934,-0.006778096523200723, &
        -0.00657441467855574,-0.006391592662227044, &
        -0.00623800806451757,-0.006120893611210529, &
        -0.006046727902490881,-0.006022003432101048/
!
      b2 = b*b
      b25 = b**2.5
      h0 = 0.0
      do i=1,4
       h(i) = b/( y(i)**4 + b2)
      enddo
      do i=5,na
       h(i) = b2/(0.25*y(i)**5 + b25)       
      enddo
! transform to gt grid
      gt = SQRT(1-ft)
      if(gt.gt.g(nf))gt = g(nf)
! find grid position
      if(gt.le.g(2))then
        i=1
      else
        i = INT(gt/0.025)
      endif
      j = i+1
! interpolate the coefficients
      dg = (gt - g(i))/(g(j)-g(i))
      hs = h0
      do k=1,na
        ca = a(i,k) + (a(j,k)-a(i,k))*dg
! sum up the terms
        hs = hs + ca*h(k)
      enddo
! final answer
      FLR_dHr33 = (5.0/3.0)*ft*hs
! 
      END FUNCTION FLR_dHr33
!
      REAL FUNCTION FLR_dHw113(ft,b)
!******************************************************************
!     Approxmation to the integral of (J0^2)*Fmaxwellian over a
!     wedge of velocity space -ft < (v_par/v) < ft.
!     b = (k_per*sqrt(T/m)/(eB/mc))**2
!     The approximation has the form
!     Hn = ft*(h0+a1*h1+a2*h2+a3*h3+a4*h4+a5*h5+a6*h6
!             +a7*h7+a8*h8+a9*h9+a10*h10+a11*h11++a12*h12+a13*h13)
!     where a1 - a13 are fit coefficients which are
!     functions of ft and
!     xa1 = (0.25)**4
!     xa2 = (0.5)**4
!     xa3 = (0.75)**4
!     xa4 = (1.0)**4
!     xa5 = 0.25*(1.5)**5
!     xa6 = 0.25*(2.0)**5
!     xa7 = 0.25*(2.5)**5
!     xa8 = 0.25*(3.0)**5
!     xa9 = 0.25*(4.0)**5
!     xa10 = 0.25*(6.0)**5
!     xa11 = 0.25*(9.0)**5
!     xa12 = 0.25*(15.0)**5
!     xa13 = 0.25*(24.0)**5
!     h0 = 1/(1+b)
!     h1 = b/(xa1+b**2)
!     h2 = b/(xa2+b**2)
!     h3 = b/(xa3+b**2)
!     h4 = b/(xa4+b**2)
!     h5 = b**2/(xa5+b**2.5)
!     h6 = b**2/(xa6+b**2.5)
!     h7 = b**2/(xa7+b**2.5)
!     h8 = b**2/(xa8+b**2.5)
!     h9 = b**2/(xa9+b**2.5)
!     h10 = b**2/(xa10+b**2.5)
!     h11 = b**2/(xa11+b**2.5)
!     h12 = b**2/(xa12+b**2.5)
!     h13 = b**2/(xa13+b**2.5)
!******************************************************************
      IMPLICIT NONE
      INTEGER,PARAMETER :: nf=40,na=13
      INTEGER :: i,j,k
      REAL :: y(na),h(na)
      REAL :: a(nf,na)
      REAL :: g(nf)
      REAL :: b2,b25
      REAL :: gt,dg,ca
      REAL :: h0,hs
! inputs
      REAL,INTENT(IN) :: ft,b
!
      data y &
       / 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, &
         6.0, 9.0, 15.0, 24.0 /
      data g &
       / 0.0, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, &
         0.225, 0.25, 0.275, 0.3, 0.325, 0.35, 0.375, 0.4, &
         0.425, 0.45, 0.475, 0.5, 0.525, 0.55, 0.575, 0.6, &
         0.625, 0.65, 0.675, 0.7, 0.725, 0.75, 0.775, 0.8, &
         0.825, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975, 0.995 /
!
      data a(1:nf,1) &
      / -0.00012291818154852585,-0.00012375263084951627, &
        -0.00012705548124603006,-0.00013098651928378224, &
        -0.0001300239666075509,-0.00013177749861248422, &
        -0.0001416508163252228,-0.00015550803169839256, &
        -0.00016620392704025377,-0.00017074566009234182, &
        -0.00017159930353227715,-0.00017383782191315866, &
        -0.00018175346650273205,-0.00019701140785954863, &
        -0.00021859900611786998,-0.000243889853534629, &
        -0.0002699815290798657,-0.0002947261685595226, &
        -0.00031722035541981164,-0.000337792202805387, &
        -0.00035766055972796806,-0.0003784658823848819, &
        -0.0004018290813234682,-0.00042902777021169847, &
        -0.0004608180506222581,-0.000497387812614436, &
        -0.000538406255968793,-0.000583129436013019, &
        -0.0006305267041906149,-0.000679402758511749, &
        -0.0007284999586893992,-0.000776551952746507, &
        -0.0008224421940399296,-0.0008650106003864622, &
        -0.0009032790398664226,-0.0009363362607852643, &
        -0.0009633464765768096,-0.0009835318689397354, &
        -0.0009961581541146902,-0.0010003403782987697 /
!
      data a(1:nf,2) &
      / -0.004337490438419467,-0.004408593486916912, &
        -0.004444732893871262,-0.004513084920870769, &
        -0.004746716224472536,-0.004966905189926152, &
        -0.0050474815226666525,-0.005090215286763122, &
        -0.00526942385712359,-0.005661262378387377, &
        -0.006212451332793689,-0.006808196295868374, &
        -0.007352835643260391,-0.007813653125623654, &
        -0.008221799233683269,-0.008646592106579254, &
        -0.009163283140971812,-0.009828142697037467, &
        -0.010666433178700307,-0.011672358593866283, &
        -0.012816820800042539,-0.014058236645157152, &
        -0.015352710681727944,-0.016661457616571118, &
        -0.017954829105449654,-0.019213298822106317, &
        -0.020426258365389804,-0.021589584129613204, &
        -0.02270280514345191,-0.023766462455146353, &
        -0.024780010400852648,-0.025740854056276685, &
        -0.026641335452965764,-0.027473109041528776, &
        -0.028222892422324453,-0.0288753343510586, &
        -0.029413354607232156,-0.029819040736136504, &
        -0.03007457841145844,-0.03015958147317277 /
!
      data a(1:nf,3) &
      / -0.025789416540476395,-0.025863671916492414, &
        -0.026752033343753423,-0.02776285010137891, &
        -0.0268528163535402,-0.026777820587090995, &
        -0.02955528498350568,-0.033673506613657755, &
        -0.03645613813116211,-0.03673869769667215, &
        -0.035395183334922864,-0.03429677800247921, &
        -0.03504395519063408,-0.03825854020503772, &
        -0.04354700092219477,-0.04989002421580002, &
        -0.05614561815707608,-0.061444719056882535, &
        -0.06538654474137906,-0.06804351162676547, &
        -0.0698390675808378,-0.0713731415263339, &
        -0.07325510396397994,-0.07597958209111222, &
        -0.07985731453905554,-0.08499685144457025, &
        -0.0913242428176313,-0.0986254232853867, &
        -0.106597513712263,-0.11489872950221636, &
        -0.12319028091078232,-0.13116420569679899, &
        -0.1385767458004059,-0.1452272176309921, &
        -0.15098577275477298,-0.1557705550296975, &
        -0.15953716233820825,-0.16226276560602115, &
        -0.16392684460284235,-0.16447100550571703 /
!
      data a(1:nf,4) &
      / -0.09899293384649144,-0.10084946294926933, &
        -0.1011302345174534,-0.10205538450119378, &
        -0.1088043627695936,-0.11428737382708065, &
        -0.1132409786928783,-0.10941771762362462, &
        -0.10956627048746492,-0.11656133481170383, &
        -0.12802655812453567,-0.13899099079728106, &
        -0.14514680123880308,-0.14468281373946246, &
        -0.13839656320743265,-0.12870614932231228, &
        -0.11836297686876596,-0.10943536099478379, &
        -0.10280497357206575,-0.09815413251447813, &
        -0.094282732421193,-0.08956277078150854, &
        -0.08237530193949194,-0.07143726470766931, &
        -0.0559852634434137,-0.035825963542865935, &
        -0.011285648139688131,0.016901680837224875, &
        0.047731729899838204,0.08006560789079664, &
        0.11274593826724322,0.14467903891931744, &
        0.1749189408310664,0.20264034778176632, &
        0.22719625442186364,0.2480742405743276, &
        0.2648713768482409,0.2772571559940007, &
        0.28492665268140294,0.28745394505792565 /
!
      data a(1:nf,5) &
      / -0.006376438143248286,-0.005678703446335254, &
        -0.007181414270206465,-0.008507180417327342, &
        -0.0030587441824594386,0.0006837746342349682, &
        -0.003603714846206696,-0.010986064023135178, &
        -0.012610429834335024,-0.004322314309491038, &
        0.01147327536594922,0.029033446667665053, &
        0.04334096052149583,0.05244466490401989, &
        0.05764587021879841,0.062269453796300445, &
        0.07001457297002311,0.08361609614926602, &
        0.10414575356235867,0.1309488498638237, &
        0.16202779131495504,0.19463388527655678, &
        0.22586639273901976,0.25315153437607774, &
        0.2745478854122989,0.2888796411883875, &
        0.2957310136847081,0.29534700453521834, &
        0.28848448632249624,0.2762489169055351, &
        0.2599416114678207,0.24093559859516578, &
        0.22056029982770697,0.200073858767419, &
        0.18058915594392255,0.16307899744588872, &
        0.14837629590025692,0.1371861806229382, &
        0.1301087453944052,0.12775157039472929 /
!
      data a(1:nf,6) &
      / 0.027177570664762007,0.025660708555256706, &
        0.030019069202118165,0.034780580419145335, &
        0.023320137904761618,0.01817172647548837, &
        0.037079002425348784,0.06844564789619578, &
        0.08974929620717964,0.09029991374164348, &
        0.07675836283085558,0.06499550453404435, &
        0.06914834574169637,0.09482386270215615, &
        0.1379794637879581,0.18783617939276331, &
        0.2313128627924559,0.2569871646789491, &
        0.2575565093754758,0.23063525011110664, &
        0.17825163066590288,0.10560141337121365, &
        0.01957505245077118,-0.07257663622752686, &
        -0.16424535314961952,-0.25008372936185164, &
        -0.3262765411701469,-0.3905460154948117, &
        -0.4419829042948522,-0.48078482669577016, &
        -0.5079654459657901,-0.5250811103553854, &
        -0.533968296299798,-0.5365968807075228, &
        -0.5348823617886482,-0.5306091732121783, &
        -0.5253617203540606,-0.5204788174681856, &
        -0.5170289758707942,-0.5158210314604874 /
!
      data a(1:nf,7) &
      / -0.026372063603654183,-0.02279738648994778, &
        -0.030477397857253452,-0.037794455199759325, &
        -0.00939980517958805,0.009403815519082537, &
        -0.018262663726081207,-0.0680292795978279, &
        -0.09089953648681193,-0.06195155056397739, &
        0.005267972306630231,0.07375742911626848, &
        0.10607920078327226,0.08054572994722053, &
        -0.0037416917029386987,-0.13066104130892753, &
        -0.2750429500969378,-0.41081291792853647, &
        -0.5169218876480672,-0.5803207592929356, &
        -0.5963392649702908,-0.5673339915198976, &
        -0.5005196294106693,-0.4056994907480056, &
        -0.293336018214041,-0.1731546737880798, &
        -0.053299059767304646,0.06004487236804623, &
        0.1626780757255154,0.2521944294437535, &
        0.32766885362617326,0.3893090864225345, &
        0.4380909891912348,0.475501904284213, &
        0.5032448280184814,0.5230635316340198, &
        0.5365860646977645,0.545201480352514, &
        0.5499677169652563,0.5514437316786911 /
!
      data a(1:nf,8) &
      / -0.00625452482978095,-0.010416059191566518, &
        -0.001837935506956212,0.0074448059579057535, &
        -0.02114644284628575,-0.03645015952019959, &
        0.004726101489112122,0.07532648626612293, &
        0.11557282717997908,0.08838736835677219, &
        -0.000215820077887674,-0.11361432793391568, &
        -0.20828156159367683,-0.252490968921111, &
        -0.2343238009491586,-0.16036367103460591, &
        -0.04923305210270379,0.07603021407758648, &
        0.19416763599623532,0.28941698163847396, &
        0.3528196974926264,0.38176272264811706, &
        0.37853625228503907,0.3485698691497733, &
        0.29879096227258295,0.23633494367827068, &
        0.16767439318924515,0.09813397367500065, &
        0.031709885855513154,-0.028898797279542254, &
        -0.08213176697345648,-0.12735127914181832, &
        -0.16459048291547096,-0.19436840361777086, &
        -0.21746401051399267,-0.23477675909638607, &
        -0.24719798343637045,-0.2555074745988286, &
        -0.26029535462350917,-0.2618135305648628 /
!
      data a(1:nf,9) &
      / -0.014008491750317154,-0.010858175579017582, &
        -0.014690660926245824,-0.015761620614778815, &
        0.01464447972770333,0.03910597699451668, &
        0.017355992072813797,-0.04339672372302816, &
        -0.10486212034317854,-0.13337007660116607, &
        -0.11833310862022395,-0.07042592050882823, &
        -0.010451770597382337,0.04151937789201676, &
        0.07250837919265096,0.07821919335872429, &
        0.06144968616200275,0.029176814055095868, &
        -0.010307945322355128,-0.049496048796219494, &
        -0.08281126920765991,-0.10695048986674838, &
        -0.12068142913632585,-0.1243622088058931, &
        -0.11938555224841174,-0.10767388546887524, &
        -0.09128581316653617,-0.07214835131302855, &
        -0.051902456974487254,-0.031837591871321624, &
        -0.012888376560811876,0.004330913752435994, &
        0.0194676674054024,0.03237615024360607, &
        0.04304968282263044,0.05156948665906291, &
        0.05805638992308726,0.06262932414316424, &
        0.06537178272882738,0.06626058438614202 /
!
      data a(1:nf,10) &
      / -0.00984897811154295,-0.011213861929878277, &
        0.0012627451468573136,0.01941374484061742, &
        -0.003406601605029813,-0.048003219521157714, &
        -0.07005618536760494,-0.05470962561924697, &
        -0.018155545739173773,0.015651831815350792, &
        0.03188101641918051,0.02836082749762614, &
        0.011364129277122603,-0.010148131909476943, &
        -0.028652566233441146,-0.03985093338495993, &
        -0.04267261453849436,-0.03835268498579425, &
        -0.029289165056639144,-0.018089834102037905, &
        -0.006962454864817502,0.002562340739394875, &
        0.009661844996713187,0.014107576054367299, &
        0.016092146732905155,0.016054256814645118, &
        0.014533169216485287,0.012064290565305558, &
        0.009115043843360704,0.006054004665941082, &
        0.0031439570816537143,0.0005495812643497855, &
        -0.0016428522137874602,-0.0034124813988241076, &
        -0.004778238574853821,-0.005785197916573281, &
        -0.006489628654333224,-0.006946983561648912, &
        -0.007203437489834874,-0.007283427867099679 /
!
      data a(1:nf,11) &
      / -0.007838136118266448,0.0023283823846186325, &
        0.001350735558413696,-0.03198178458521528, &
        -0.040596949464847265,-0.014980284228098828, &
        0.012791580189536633,0.019819783531897617, &
        0.007271275441677716,-0.011886331165359865, &
        -0.0262525339412496,-0.03110547627515836, &
        -0.02746227703920945,-0.019147853285635952, &
        -0.010223201311571284,-0.0035626382723261685, &
        -0.0004672150329625424,-0.0009112845231917799, &
        -0.00403890933760076,-0.008656764137227269, &
        -0.01360137615094259,-0.017952400198493423, &
        -0.021113970947071126,-0.02280533204267643, &
        -0.023001538474963334,-0.021855988348116417, &
        -0.019625406994920214,-0.016608135320441697, &
        -0.013099497826239093,-0.00936373068355234, &
        -0.0056195452451353045,-0.00203576662275351, &
        0.0012658540684111763,0.004204433896718651, &
        0.006731066917082362,0.008820075004633776, &
        0.010459947060637909,0.011645075021106455, &
        0.01236858575520039,0.012605261987329808 /
!
      data a(1:nf,12) &
      / -0.005039917187131149,-0.0029013239412580827, &
        -0.024496509753129114,-0.006786436881478286, &
        0.01106983509522319,0.004712881019976245, &
        -0.011312489664947023,-0.019457127209444147, &
        -0.01638408334685515,-0.008003856477868387, &
        -0.0008415088360402834,0.0019291226677922158, &
        0.0004012449126942208,-0.0036058000357095388, &
        -0.00801660331642573,-0.011332947267893578, &
        -0.012870604210423853,-0.012631556175140224, &
        -0.01103644643052104,-0.0086631699346853, &
        -0.006059339539222552,-0.0036405843282927286, &
        -0.0016594081760095136,-0.00021950624531985596, &
        0.0006876199576112407,0.0011406916704034264, &
        0.0012514544748930545,0.001136719835680533, &
        0.0009004434494527258,0.0006245666075586076, &
        0.00036636584786098325,0.00015932496829951237, &
        0.00001877796219751371,-0.00005610040197390198, &
        -0.00007504145200787349,-0.00005390788669945579, &
        -0.000010858216070053217,0.000036449164767748954, &
        0.00007251477797254324,0.00008554944171834666 /
!
      data a(1:nf,13) &
      / -0.004051448637862887,-0.012710130801553365, & 
        0.004123873000763445,-0.00079524839806977, &
        -0.011795890485173155,-0.010660987764235186, &
        -0.002896230830155222,0.0017138601566251133, &
        0.00044455948714605675,-0.004059851259818226, &
        -0.008315848655032099,-0.010388721340253149, &
        -0.010085141914889073,-0.008238370611403312, &
        -0.005940061549884135,-0.004056071559860097, &
        -0.0030460609094442592,-0.0029882442919335617, &
        -0.003698817078671146,-0.004867744899458737, &
        -0.006168002777779735,-0.007325451991562426, &
        -0.008150961417197289,-0.008544807514787145, &
        -0.008484275269092567,-0.008003669996159202, &
        -0.007173265967389542,-0.00608077451028044, &
        -0.004817263337921096,-0.0034675722223658223, &
        -0.0021048852864626877,-0.0007880084210414307, &
        0.00043828047126626135,0.0015418095470033366, &
        0.002501029414641698,0.0033022891471159443, &
        0.003937044506494347,0.00439923872933079, &
        0.004682915960920964,0.0047759730286914825 /
!
      b2 = b*b
      b25 = b**2.5
      h0 = 0.0
      do i=1,4
       h(i) = b/( y(i)**4 + b2)
      enddo
      do i=5,na
       h(i) = b2/(0.25*y(i)**5 + b25)       
      enddo
! transform to gt grid
      gt = SQRT(1-ft)
      if(gt.gt.g(nf))gt = g(nf)
! find grid position
      if(gt.le.g(2))then
        i=1
      else
        i = INT(gt/0.025)
      endif
      j = i+1
! interpolate the coefficients
      dg = (gt - g(i))/(g(j)-g(i))
      hs = h0
      do k=1,na
        ca = a(i,k) + (a(j,k)-a(i,k))*dg
! sum up the terms
        hs = hs + ca*h(k)
      enddo
! final answer
      FLR_dHw113 = hs*7.0*ft**5
! 
      END FUNCTION FLR_dHw113
!
      REAL FUNCTION FLR_dHw133(ft,b)
!******************************************************************
!     Approxmation to the integral of (J0^2)*Fmaxwellian over a
!     wedge of velocity space -ft < (v_par/v) < ft.
!     b = (k_per*sqrt(T/m)/(eB/mc))**2
!     The approximation has the form
!     Hn = ft*(h0+a1*h1+a2*h2+a3*h3+a4*h4+a5*h5+a6*h6
!             +a7*h7+a8*h8+a9*h9+a10*h10+a11*h11++a12*h12+a13*h13)
!     where a1 - a13 are fit coefficients which are
!     functions of ft and
!     xa1 = (0.25)**4
!     xa2 = (0.5)**4
!     xa3 = (0.75)**4
!     xa4 = (1.0)**4
!     xa5 = 0.25*(1.5)**5
!     xa6 = 0.25*(2.0)**5
!     xa7 = 0.25*(2.5)**5
!     xa8 = 0.25*(3.0)**5
!     xa9 = 0.25*(4.0)**5
!     xa10 = 0.25*(6.0)**5
!     xa11 = 0.25*(9.0)**5
!     xa12 = 0.25*(15.0)**5
!     xa13 = 0.25*(24.0)**5
!     h0 = 1/(1+b)
!     h1 = b/(xa1+b**2)
!     h2 = b/(xa2+b**2)
!     h3 = b/(xa3+b**2)
!     h4 = b/(xa4+b**2)
!     h5 = b**2/(xa5+b**2.5)
!     h6 = b**2/(xa6+b**2.5)
!     h7 = b**2/(xa7+b**2.5)
!     h8 = b**2/(xa8+b**2.5)
!     h9 = b**2/(xa9+b**2.5)
!     h10 = b**2/(xa10+b**2.5)
!     h11 = b**2/(xa11+b**2.5)
!     h12 = b**2/(xa12+b**2.5)
!     h13 = b**2/(xa13+b**2.5)
!******************************************************************
      IMPLICIT NONE
      INTEGER,PARAMETER :: nf=40,na=13
      INTEGER :: i,j,k
      REAL :: y(na),h(na)
      REAL :: a(nf,na)
      REAL :: g(nf)
      REAL :: b2,b25
      REAL :: gt,dg,ca
      REAL :: h0,hs
! inputs
      REAL,INTENT(IN) :: ft,b
!
      data y &
       / 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, &
         6.0, 9.0, 15.0, 24.0 /
      data g &
       / 0.0, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, &
         0.225, 0.25, 0.275, 0.3, 0.325, 0.35, 0.375, 0.4, &
         0.425, 0.45, 0.475, 0.5, 0.525, 0.55, 0.575, 0.6, &
         0.625, 0.65, 0.675, 0.7, 0.725, 0.75, 0.775, 0.8, &
         0.825, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975, 0.995 /
!
      data a(1:nf,1) &
      / -0.00021099192118256838,-0.00021215432152454916, &
        -0.00021496794078679216,-0.0002184834618284004, &
        -0.00021941654558710333,-0.00022236560345456624, &
        -0.00023054234474129132,-0.00024140332656386043, &
        -0.00025065651787525844,-0.0002565888453816001, &
        -0.0002608008638638193,-0.0002664239352132558, &
        -0.00027603751579952096,-0.0002905707088558529, &
        -0.0003092997186611157,-0.00033054345666383433, &
        -0.0003525021614725343,-0.00037390291791483676, &
        -0.0003942847702524954,-0.0004139702307464077, &
        -0.0004338300882209545,-0.0004549712564694275, &
        -0.0004784477632904238,-0.0005050564419257775, &
        -0.0005352091135584214,-0.0005689241641615439, &
        -0.0006058623700062604,-0.0006453965891872748, &
        -0.0006867023311796228,-0.0007288423617886064, &
        -0.0007708375233062004,-0.0008117191175689431, &
        -0.0008505628778718854,-0.0008865065057619415, &
        -0.0009187544439011042,-0.000946573693707986, &
        -0.0009692844353736962,-0.000986248091849884, &
        -0.000996855899839666,-0.0010003701997520764 /
!
      data a(1:nf,2) &
      / -0.007516118124619053,-0.007582653238938919, &
        -0.0076342482696716285,-0.007717619946425058, &
        -0.007912275240496758,-0.008110604773897245, &
        -0.008238136660917217,-0.008357959075613097, &
        -0.00857542931809041,-0.008934964749197194, &
        -0.009402441872275125,-0.009907613537582405, &
        -0.010393740847795874,-0.010843457621503294, &
        -0.011278773768073291,-0.011744144219876684, &
        -0.01228631426532667,-0.012938552618297683, &
        -0.013713667488640668,-0.014604396509227957, &
        -0.015588738798927082,-0.016637126131197943, &
        -0.017719061869290398,-0.01880785442289601, &
        -0.019883216640246504,-0.02093179613428897, &
        -0.02194603978368015,-0.022922659607983453, &
        -0.023860610571888596,-0.02475927288492752, &
        -0.02561701354887691,-0.026430228044694637, &
        -0.027192842147565788,-0.027896229073599343, &
        -0.028529436907107186,-0.02907964016534592, &
        -0.029532728426256316,-0.029873977513902528, &
        -0.030088748064251103,-0.030160149517984625 /
!
      data a(1:nf,3) &
      / -0.04092965887698519,-0.04108791365041209, &
        -0.041763798914349115,-0.0425664958148011, &
        -0.042278363337207914,-0.04256324420863211, &
        -0.04461782800418379,-0.047509554123032816, &
        -0.04962694813025692,-0.0502981180991402, &
        -0.05008861592961211,-0.050149991249548, &
        -0.05143946941445865,-0.05429915348307257, &
        -0.05844359709436264,-0.06321572470635073, &
        -0.06790231607299124,-0.07198280269361823, &
        -0.07524345444452774,-0.07777308816600903, &
        -0.07987848551630261,-0.0819682367379857, &
        -0.08444341412682776,-0.08761775738184019, &
        -0.0916744799274587,-0.09665451715586393, &
        -0.10247622032857273,-0.10896266845117664, &
        -0.11587842040592,-0.12296421003080438, &
        -0.1299660012955417,-0.1366561446309582, &
        -0.1428463246014633,-0.14839273386785035, &
        -0.15319480246229555,-0.15718889908364564, &
        -0.1603384260088896,-0.16262154429987885, &
        -0.16401773974421308,-0.1644748074882163 /
!
      data a(1:nf,4) &
      / -0.08847008213185137,-0.08950266057354117, &
        -0.08956185900412872,-0.08996295761725803, &
        -0.09378179777065743,-0.09672447409018736, &
        -0.09565159344782037,-0.0928608899057608, &
        -0.09239446228466153,-0.0958843333953816, &
        -0.1017727372178554,-0.10697509329509525, &
        -0.10888489003418633,-0.10646645546812472, &
        -0.10029654359164386,-0.09191637778103923, &
        -0.08302373128178075,-0.07483562439194014, &
        -0.0677936479278527,-0.061575909426923925, &
        -0.055317050127185574,-0.047911243168033746, &
        -0.03829861846785132,-0.02567582409925162, &
        -0.009610567767488254,0.009925448413904636, &
        0.03260752800772515,0.057830581929129155, &
        0.08480769084824225,0.1126647160849239, &
        0.14052111495654585,0.16755062498063844, &
        0.19302051955117072,0.21631032366319747, &
        0.23691314734665203,0.25442331729582923, &
        0.26851403894306414,0.27890857729656493, &
        0.2853484322805304,0.28747135871112917 /
!
      data a(1:nf,5) &
      / 0.045453688373416234,0.04626290143985179, &
        0.04585367945437757,0.04577111442596715, &
        0.0499589016652473,0.05330966706851292, &
        0.05208762198549877,0.04935155804204483, &
        0.05044740361272304,0.057794683933097946, &
        0.06982437678311504,0.08299637982910524, &
        0.09429715873010602,0.1026383982728325, &
        0.10892961993416905,0.11527147675679583, &
        0.1239146997798235,0.13642253622273448, &
        0.15325643130012767,0.17376356002468496, &
        0.1964473317261186,0.2193656620521789, &
        0.24052762546858075,0.25820676470338244, &
        0.2711363746739077,0.2786000026369567, &
        0.280407411908587,0.27683130318548554, &
        0.26849973901749846,0.2562818191500492, &
        0.24118083185905492,0.2242449310287703, &
        0.20649930448060227,0.1889004118981693, &
        0.17231010099450683,0.15748585485200262, &
        0.14508309002105757,0.13566507695607855, &
        0.12971578300635223,0.12773533987413643 /
!
      data a(1:nf,6) &
      / 0.013603933077689256,0.012591474509046074, &
        0.015079012475345044,0.017716598924231275, &
        0.010581749781845007,0.007259092335027417, &
        0.018248038775213082,0.03632319309268495, &
        0.047845675319290404,0.04653361035465098, &
        0.036627082512225684,0.02776866231010633, &
        0.02826599415022646,0.041004639312099744, &
        0.06288802887621081,0.08677855628285569, &
        0.10433223567017424,0.10850111132573803, &
        0.09504472730320379,0.06297632126508756, &
        0.014178569506898775,-0.04744856825517485, &
        -0.11697402519166666,-0.1893205368009927, &
        -0.2599203016796996,-0.3251409560240426, &
        -0.38242253093690515,-0.4302941182747253, &
        -0.46824573921545687,-0.49654398632720564, &
        -0.5160279201073619,-0.5279145091664055, &
        -0.5336296737802027,-0.5346724705739401, &
        -0.5325133419956618,-0.528522020749775, &
        -0.5239192120489788,-0.5197440963231221, &
        -0.5168285606668845,-0.5158125215880243 /
!
      data a(1:nf,7) &
      / -0.07050900771183866,-0.0686975577860276, &
        -0.07373385238843966,-0.07866448446474927, &
        -0.06242096129883623,-0.05237518033556077, &
        -0.07043540472773335,-0.10153706526429085, &
        -0.11618472431056714,-0.09994946353896172, &
        -0.06181656281371972,-0.024650018414544395, &
        -0.010941871352922128,-0.032963310527291956, &
        -0.08989064133836466,-0.17041318237790737, &
        -0.25789054820827384,-0.33554691071613396, &
        -0.3901927617015346,-0.41403975247359126, &
        -0.40484819105767644,-0.364971289777271, &
        -0.2998876714125851,-0.21668467065511643, &
        -0.12278708351305312,-0.024986710734075768, &
        0.07104398063872697,0.1610385498868735, &
        0.24213512599639442,0.3127362290975515, &
        0.3722883517760024,0.42103463757413195, &
        0.45977409743540043,0.48964758964309485, &
        0.5119597863176426,0.5280361378086809, &
        0.5391110170919904,0.546238573808289, &
        0.550215585719926,0.5514534563665242 /
!
      data a(1:nf,8) &
      / 0.02103968515659871,0.018749726693712776, &
        0.02416594695083818,0.030039876044564595, &
        0.013363681266622685,0.005014909438857185, &
        0.030604783735602048,0.07330614924915957, &
        0.0970420702127327,0.08013790795081499, &
        0.02732821326222068,-0.03828562118749379, &
        -0.08998446200538379,-0.1089120233478913, &
        -0.08883641275171583,-0.035095649581409205, &
        0.039624910198758356,0.12018528129120487, &
        0.19291409333782994,0.2479270425793969, &
        0.2798934726627267,0.2876724675326394, &
        0.2733245745173756,0.24092719367772997, &
        0.19549218584558548,0.14205819672395614, &
        0.08518864395550296,0.028604635884938667, &
        -0.024917898760627555,-0.07352083974121615, &
        -0.11615020091855266,-0.1524017856018114, &
        -0.1823560916626552,-0.20642259638960136, &
        -0.2252043588140924,-0.23938475913918444, &
        -0.24963608694360828,-0.25654522717870565, &
        -0.26054968111939036,-0.26182364162586635 /
!
      data a(1:nf,9) &
      / -0.02267690620648155,-0.020853142754175824, &
        -0.023246021122883098,-0.02398655869104238, &
        -0.005980840914769692,0.008127357531870771, &
        -0.0056389618983803125,-0.042427349745674636, &
        -0.07882659085443677,-0.09465226101777025, &
        -0.08403049565035081,-0.054063270014842235, &
        -0.017848556590575493,0.012252622790507428, &
        0.02849947340466963,0.028669755413650022, &
        0.014981256784400099,-0.0077980630459040245, &
        -0.03418288668465741,-0.059300688974529114, &
        -0.07960860555352034,-0.09308989363006326, &
        -0.09910239631491424,-0.09805028667273619, &
        -0.09101876607690973,-0.07940238951241763, &
        -0.06469142843808795,-0.04826165784184866, &
        -0.03127959321588003,-0.014659766686568876, &
        0.000935731765266401,0.015072129484134389, &
        0.02750275149427417,0.03812495846413677, &
        0.04693728272866049,0.05399983101759953, &
        0.05939977729285373,0.06322154527160595, &
        0.06552035019571101,0.06626658759955795 /
!
      data a(1:nf,10) &
      / -0.008918951998630131,-0.009731508499073227, &
        -0.002246677334367675,0.008534718066663394, &
        -0.005288410945825994,-0.03188947369196704, &
        -0.04456670555932152,-0.03466355885865724, &
        -0.012293645829144273,0.00793659909365163, &
        0.01717037366725427,0.014356135519150826, &
        0.003577728137044933,-0.009533464680362869, &
        -0.02034367688825789,-0.026295397720626612, &
        -0.02687415641950297,-0.023002209882333036, &
        -0.016304062467434433,-0.008498812608818351, &
        -0.0010154439450344999,0.0051774491554359425, &
        0.009583877958216247,0.012100455274225319, &
        0.01290774458565025,0.012319130757485364, &
        0.01073268694619689,0.008529947576679908, &
        0.006041107103704646,0.003524793265988757, &
        0.0011642528144465425,-0.0009267228000346961, &
        -0.002691871993917294,-0.004118560691564799, &
        -0.0052244307541455814,-0.006045057923909525, &
        -0.00662365351408245,-0.00700271947315767, &
        -0.0072168584828951055,-0.007283940735384631 /
!
      data a(1:nf,11) &
      / -0.009454838622372286,-0.0033747396034542465, &
        -0.004028638604561596,-0.024007904457737483, &
        -0.028984475859608283,-0.013465325799545624, &
        0.0030802974517105675,0.006942564627486547, &
        -0.000952378635537432,-0.012613818172834579, &
        -0.0211289828336112,-0.023731197425967965, &
        -0.021179015590600425,-0.01591431678508548, &
        -0.01046884925933722,-0.006592072704649121, &
        -0.005027106945724391,-0.005680711113942083, &
        -0.007946275624467036,-0.011020610139283349, &
        -0.014138197475893843,-0.016705466141295977, &
        -0.018349914793863054,-0.018910759120774445, &
        -0.01840869879242457,-0.01695058460079446, &
        -0.014740430786835956,-0.011998047575626725, &
        -0.008938606090050039,-0.005753732969049352, &
        -0.002602140005407705,0.0003934534332633577, &
        0.0031441625290369313,0.00559013095463301, &
        0.007694632459363082,0.009437439374961332, &
        0.010808282477788977,0.011801063192635963, &
        0.012408110760064872,0.012606868908186408 /
!
      data a(1:nf,12) &
      / -0.0049530754615332295,-0.0036860483572920266, &
        -0.016596872750482605,-0.005903058059625987, &
        0.004749609789748477,0.000814230403579197, &
        -0.008815553974327539,-0.013565406529131963, &
        -0.011531272818870564,-0.006379686296766662, &
        -0.0020774128507643486,-0.0005091528046776705, &
        -0.0015458529611884764,-0.004019888273323113, &
        -0.006634793247972892,-0.00847585966331299, &
        -0.009148512794365127,-0.008689807128358096, &
        -0.007395042952766051,-0.005650441393761696, &
        -0.0038140398945734066,-0.0021523039005730826, &
        -0.0008219820258978938,0.00011888193146225579, &
        0.0007008990719956643,0.0009570629160070565, &
        0.0009818501387315326,0.0008584488527589862, &
        0.00065985815140579,0.0004429907489466034, &
        0.0002471255867746436,0.00009490272095102892, &
        -4.806822015979151e-6,-0.00005346540342010542, &
        -0.000059632250945185206,-0.000036041703672440306, &
        3.189653592583497e-6,0.00004408404303835045, &
        0.00007467182178100273,0.00008564840951297281 /
!
      data a(1:nf,13) &
      / -0.004584398039556792,-0.009770773553784196, &
        0.0003216386554683037,-0.0026766173509775726, &
        -0.00926520047976287,-0.00852685252302643, &
        -0.0038534968915098222,-0.0011496150762506652, &
        -0.002008686623485856,-0.004772308961799254, &
        -0.007317442559895859,-0.008488859002033777, &
        -0.008204978691192855,-0.007006262348884906, &
        -0.005577623351493033,-0.004450019369705904, &
        -0.0038926995216486127,-0.003933108381547923, &
        -0.004435510899232664,-0.005188535392230742, &
        -0.005975393841914922,-0.006618037713559399, &
        -0.006997080711682813,-0.007053657147542136, &
        -0.00678202537435979,-0.006210314562479224, &
        -0.005394149228184153,-0.004399105587299168, &
        -0.0032921697210512857,-0.0021354956080661402, &
        -0.0009825451665379514,0.00012336222228015448, &
        0.0011491198000975755,0.002070754785264861, &
        0.002871825427154917,0.0035416343652103777, &
        0.004072924551817958,0.004460369506028794, &
        0.0046984524800082195,0.004776607760963558 /
!
      b2 = b*b
      b25 = b**2.5
      h0 = 0.0
      do i=1,4
       h(i) = b/( y(i)**4 + b2)
      enddo
      do i=5,na
       h(i) = b2/(0.25*y(i)**5 + b25)       
      enddo
! transform to gt grid
      gt = SQRT(1-ft)
      if(gt.gt.g(nf))gt = g(nf)
! find grid position
      if(gt.le.g(2))then
        i=1
      else
        i = INT(gt/0.025)
      endif
      j = i+1
! interpolate the coefficients
      dg = (gt - g(i))/(g(j)-g(i))
      hs = h0
      do k=1,na
        ca = a(i,k) + (a(j,k)-a(i,k))*dg
! sum up the terms
        hs = hs + ca*h(k)
      enddo
! final answer
      FLR_dHw133 = hs*(35.0/9.0)*ft**3
! 
      END FUNCTION FLR_dHw133
!
      REAL FUNCTION FLR_dHw333(ft,b)
!******************************************************************
!     Approxmation to the integral of (J0^2)*Fmaxwellian over a
!     wedge of velocity space -ft < (v_par/v) < ft.
!     b = (k_per*sqrt(T/m)/(eB/mc))**2
!     The approximation has the form
!     Hn = ft*(h0+a1*h1+a2*h2+a3*h3+a4*h4+a5*h5+a6*h6
!             +a7*h7+a8*h8+a9*h9+a10*h10+a11*h11++a12*h12+a13*h13)
!     where a1 - a13 are fit coefficients which are
!     functions of ft and
!     xa1 = (0.25)**4
!     xa2 = (0.5)**4
!     xa3 = (0.75)**4
!     xa4 = (1.0)**4
!     xa5 = 0.25*(1.5)**5
!     xa6 = 0.25*(2.0)**5
!     xa7 = 0.25*(2.5)**5
!     xa8 = 0.25*(3.0)**5
!     xa9 = 0.25*(4.0)**5
!     xa10 = 0.25*(6.0)**5
!     xa11 = 0.25*(9.0)**5
!     xa12 = 0.25*(15.0)**5
!     xa13 = 0.25*(24.0)**5
!     h0 = 1/(1+b)
!     h1 = b/(xa1+b**2)
!     h2 = b/(xa2+b**2)
!     h3 = b/(xa3+b**2)
!     h4 = b/(xa4+b**2)
!     h5 = b**2/(xa5+b**2.5)
!     h6 = b**2/(xa6+b**2.5)
!     h7 = b**2/(xa7+b**2.5)
!     h8 = b**2/(xa8+b**2.5)
!     h9 = b**2/(xa9+b**2.5)
!     h10 = b**2/(xa10+b**2.5)
!     h11 = b**2/(xa11+b**2.5)
!     h12 = b**2/(xa12+b**2.5)
!     h13 = b**2/(xa13+b**2.5)
!******************************************************************
      IMPLICIT NONE
      INTEGER,PARAMETER :: nf=40,na=13
      INTEGER :: i,j,k
      REAL :: y(na),h(na)
      REAL :: a(nf,na)
      REAL :: g(nf)
      REAL :: b2,b25
      REAL :: gt,dg,ca
      REAL :: h0,hs
! inputs
      REAL,INTENT(IN) :: ft,b
!
      data y &
       / 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, &
         6.0, 9.0, 15.0, 24.0 /
      data g &
       / 0.0, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, &
         0.225, 0.25, 0.275, 0.3, 0.325, 0.35, 0.375, 0.4, &
         0.425, 0.45, 0.475, 0.5, 0.525, 0.55, 0.575, 0.6, &
         0.625, 0.65, 0.675, 0.7, 0.725, 0.75, 0.775, 0.8, &
         0.825, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975, 0.995 /
!
      data a(1:nf,1) &
      / -0.0005021661447917047,-0.0005032821076609295, &
        -0.0005051324505362942,-0.0005075804931204964, &
        -0.0005095413350433287,-0.0005125595986448239, &
        -0.0005176962985093958,-0.0005240849070395903, &
        -0.0005302957952606491,-0.0005357877379057463, &
        -0.000541135601926257,-0.0005474105442366817, &
        -0.0005554722584792948,-0.0005656012287806877, &
        -0.0005775157319480861,-0.0005906200014285101, &
        -0.0006042999299704643,-0.0006181410123060616, &
        -0.0006320224600919744,-0.000646095073041586, &
        -0.0006606918344762569,-0.0006762102842052808, &
        -0.0006930068206507656,-0.0007113214829623793, &
        -0.0007312392295752979,-0.0007526836441615892, &
        -0.0007754341036469992,-0.0007991563497466991, &
        -0.0008234375311415576,-0.0008478200414279158, &
        -0.0008718292501548319,-0.0008949944200309368, &
        -0.0009168621249662223,-0.0009370032252962179, &
        -0.0009550148502899258,-0.000970518748484972, &
        -0.0009831573778547043,-0.0009925892401509984, &
        -0.0009984841801797817,-0.0010004366837319534 /
!
      data a(1:nf,2) &
      / -0.01676626001942383,-0.016811573104707822, &
        -0.01685770925135932,-0.016926159234465388, &
        -0.017043303927763986,-0.017173132044786987, &
        -0.017291208861396345,-0.017419186780930795, &
        -0.0175923438434058,-0.017824954508363522, &
        -0.018104818261202915,-0.0184079951677234, &
        -0.018715624661998165,-0.0190226145707173, &
        -0.019337136922951137,-0.019674600793532226, &
        -0.02005048702772938,-0.02047503209976509, &
        -0.020950880208749467,-0.021473459807623607, &
        -0.022033023988279776,-0.02261732212763956, &
        -0.023214015544019473,-0.023812380893777796, &
        -0.0244041678075535,-0.02498371260650023, &
        -0.025547526272024168,-0.026093594461652492, &
        -0.026620598975469312,-0.02712720351248743, &
        -0.027611496563078042,-0.028070622971594572, &
        -0.028500602234052597,-0.028896310974700024, &
        -0.029251591318232073,-0.029559451941255643, &
        -0.029812330949715138,-0.030002389499496296, &
        -0.0301218256080793,-0.030161501544106573 /
!
      data a(1:nf,3) &
      / -0.08747743211009329,-0.08764661004136487, &
        -0.08801769477193111,-0.08848815330301929, &
        -0.08865491552792015,-0.08907598149043983, &
        -0.09014556474079405,-0.09154443771031895, &
        -0.09273415125973661,-0.09350006964875601, &
        -0.09404439541772325,-0.09475947266688478, &
        -0.09596266686374944,-0.09775518112593096, &
        -0.10002488187723135,-0.10253773460922555, &
        -0.10504945570672564,-0.10738969053652747, &
        -0.10949981573102652,-0.11142762427791064, &
        -0.11329448908519368,-0.11525197684642308, &
        -0.11744190504322227,-0.11996769897926418, &
        -0.12287961287805627,-0.1261725331748953, &
        -0.12979302792760739,-0.1336518441341279, &
        -0.1376383428805159,-0.14163430723543247, &
        -0.14552541242352235,-0.1492095558210862, &
        -0.152601852241544,-0.15563651705041526, &
        -0.15826611535976598,-0.16045872690880336, &
        -0.162193560555135,-0.16345559350336014, &
        -0.16422961889047993,-0.16448342735248977 /
!
      data a(1:nf,4) &
      / 0.024018069193986302,0.02395640723199588, &
        0.02429376059764643,0.024662974659645287, &
        0.02405194043789649,0.023908268847112613, &
        0.025275639725348853,0.027377258762718015, &
        0.02887425126351484,0.029259459617459127, &
        0.029101124361500252,0.029463803852806447, &
        0.03122645066301133,0.034715933533281396, &
        0.03970996850831854,0.045668834959039395, &
        0.05202061269766811,0.05837693391365162, &
        0.06462980898529103,0.07093737952129009, &
        0.07763810363200407,0.08513733368265408, &
        0.09380242254639803,0.1038870380145529, &
        0.11549165627346702,0.12855720590155784, &
        0.14288338912277132,0.15816188684137222, &
        0.17401527457350074,0.1900348350490485, &
        0.2058126543461093,0.22096571985144753, &
        0.23515133222891238,0.2480742929207138, &
        0.2594870085491582,0.2691839363530637, &
        0.2769917929127569,0.2827571011244836, &
        0.2863322523680466,0.2875114682855431 /
!
      data a(1:nf,5) &
      / 0.13589385571168056,0.13638989924690303, &
        0.13653661253908334,0.13691100437636483, &
        0.13882029822763464,0.14055255582248272, &
        0.14087805687913524,0.14084296975959132, &
        0.14223165945446115,0.14581605232428818, &
        0.1510202902134128,0.156622892674401, &
        0.1616064285677098,0.16562411284499934, &
        0.1690035652087566,0.17245576105764593, &
        0.17670756689219103,0.1822160634178689, &
        0.1890327365646377,0.19681261234798164, &
        0.20492202952002314,0.21258969112755677, &
        0.21905410974066053,0.22367851037594877, &
        0.22602134776420169,0.22586355105936073, &
        0.223201223993312,0.21821502417844452, &
        0.21122744546480265,0.20265685342225304, &
        0.1929746997223829,0.18266959804582616, &
        0.1722199909970672,0.16207555171684207, &
        0.1526465134756384,0.14429951099960814, &
        0.13735834823038406,0.13210774210036103, &
        0.12879836059531724,0.12769782901179627 /
!
      data a(1:nf,6) &
      / -0.1932302990971495,-0.1940860897918999, &
        -0.19390760222574643,-0.19395842028268384, &
        -0.19753391778371254,-0.20009235967157335, &
        -0.19819256192646506,-0.19434142198774573, &
        -0.193133267571731,-0.19662605079484274, &
        -0.20333639774752732,-0.2100275756764174, &
        -0.2140071923843836,-0.21449027433910972, &
        -0.21273295343504506,-0.21131564643126, &
        -0.21313108690247762,-0.22051234928312624, &
        -0.234718934296255,-0.25580807555233265, &
        -0.2828005854034217,-0.3140104595293831, &
        -0.3474171342336454,-0.3809960338531384, &
        -0.41296305396740096,-0.44192158113087343, &
        -0.46692236100594586,-0.48745615354295935, &
        -0.5034025687296395,-0.5149553539339857, &
        -0.5225404112024763,-0.5267372760145097, &
        -0.5282106774677433,-0.5276549384123927, &
        -0.5257515902136038,-0.5231386456373428, &
        -0.5203892408474124,-0.517996111673019, &
        -0.5163587470878386,-0.5157930593711408 /
!
      data a(1:nf,7) &
      / 0.07723733423632273,0.07821149623368762, &
        0.07699682316161716,0.0760320483314939, &
        0.0823025195413386,0.08663542893397325, &
        0.08181992984932984,0.0730443802587265, &
        0.0701889141648272,0.07791286474511405, &
        0.09297510196658776,0.10759909243666443, &
        0.11436667254925925,0.10953867397976347, &
        0.09393720620964752,0.07194243614624884, &
        0.049627752648504386,0.032920568831193764, &
        0.026304003629301853,0.03221106974617194, &
        0.051015934519029305,0.0814174930803877, &
        0.12100051791260047,0.1668077414405127, &
        0.2158201654248445,0.26530101428281305, &
        0.3130009819818025,0.35724446322466585, &
        0.3969282161682238,0.43146275698254155, &
        0.4606833622553761,0.4847499946955147, &
        0.504049830834671,0.5191098040702664, &
        0.5305227239472774,0.538886896234495, &
        0.5447577929083265,0.5486077793393012, &
        0.5507906398776696,0.5514764904469697 /
!
      data a(1:nf,8) &
      / -0.01782350432979496,-0.01868398679465244, &
        -0.01699868194522569,-0.015237060541358094, &
        -0.02102577931987426,-0.02401927391760572, &
        -0.015790832658060694,-0.0021528170777491695, &
        0.004802834465018835,-0.0019577485281614404, &
        -0.020493532815351267,-0.04276890035835046, &
        -0.059750141857255024,-0.06531416686431091, &
        -0.057800863071744146,-0.03955322245359527, &
        -0.015373058881005486,0.00922585599912118, &
        0.02940219430343549,0.041786324037824896, &
        0.04472173881018038,0.03810218551396061, &
        0.022988350458599593,0.0011584081975084537, &
        -0.02530380381856201,-0.05432949193608971, &
        -0.08408011480050326,-0.11306457993849084, &
        -0.1401796739177925,-0.16469620200812618, &
        -0.18621254277102373,-0.20459184851579737, &
        -0.21989528230797162,-0.23231854023758536, &
        -0.24213588665879282,-0.24965269445022464, &
        -0.2551664381093768,-0.2589336008882255, &
        -0.2611409463700376,-0.2618474880844568 /
!
      data a(1:nf,9) &
      / -0.013457684760015187,-0.012826965986272754, &
        -0.01359888427350131,-0.013797626961156162, &
        -0.007769060255302218,-0.003120907604621248, &
        -0.007778891189063497,-0.019946949225184918, &
        -0.03168206896910619,-0.036273607232296445, &
        -0.03193568050689177,-0.021272147750057036, &
        -0.00884197675948828,0.0011578399452906174, &
        0.006205256757224831,0.005726516746732413, &
        0.0006857299378207893,-0.007098877697572004, &
        -0.015609528617340995,-0.023098063523295798, &
        -0.028333306292091187,-0.03066113582686203, &
        -0.029940084627988695,-0.026414251807115097, &
        -0.020571128652042248,-0.013014048418883428, &
        -0.004362926516056298,0.0048126426904693265, &
        0.014031776566069354,0.0229177776749242, &
        0.031197737403721226,0.03869185774024542, &
        0.0452975124498356,0.05097124481284787, &
        0.055711028308481514,0.05953974492002001, &
        0.06249063646204622,0.06459414660362628, &
        0.06586639707884057,0.06628063054827749 /
!
      data a(1:nf,10) &
      / -0.004759988925198799,-0.0050206337941498536, &
        -0.0025159378099177676,0.0010560261054035625, &
        -0.0035752793527406013,-0.01236515760519974, &
        -0.016380345711647593,-0.012826291137524315, &
        -0.005207460053341606,0.0015331363049428948, &
        0.004465622342074038,0.003332543783961084, &
        -0.00039421669445927066,-0.004755568152973844, &
        -0.008176651928457401,-0.009826737712752331, &
        -0.009593057365947288,-0.007858497634363415, &
        -0.005237608740767952,-0.0023609136870439418, &
        0.00025995858369476377,0.002289778396791786, &
        0.0035706510837660854,0.0040900110363771836, &
        0.003935668967742334,0.003251496170100565, &
        0.0022005944700800306,0.00093940751143895, &
        -0.0003982053975384314,-0.0017082505416627924, &
        -0.0029179380261392662,-0.003982914920275693, &
        -0.0048824210710805005,-0.005613621468617325, &
        -0.006186063659147534,-0.0066165569027401006, &
        -0.0069249552654939706,-0.007130329853453743, &
        -0.007248013375757534,-0.007285191067162744 /
!
      data a(1:nf,11) &
      / -0.006015840408904127,-0.003983048921647914, &
        -0.004207298222644873,-0.01083726786203254, &
        -0.012402968413892212,-0.007146109452279115, &
        -0.0016294063688471194,-0.00040868754079315117, &
        -0.0030976083314366853,-0.006961600302001436, &
        -0.009674558142780176,-0.010341835147023026, &
        -0.00926826610611875,-0.007317034623558927, &
        -0.00536030884990546,-0.003983063007437115, &
        -0.0034123644197303626,-0.0035829210875433226, &
        -0.004254736495853839,-0.005126667102812288, &
        -0.005919733403644933,-0.006423945207871151, &
        -0.006514474855318975,-0.006146721044604897, &
        -0.005339922043664291,-0.00415694939167377, &
        -0.0026846795671993107,-0.0010183686608464804, &
        0.0007499110648079288,0.0025383327769405994, &
        0.004279054387655412,0.00591916840945883, &
        0.007419743643912291,0.008753878122665038, &
        0.009904198616689452,0.010860058852415366, &
        0.011614949429821553,0.012163701339359512, &
        0.012500246517798463,0.012610617542319846 /
!
      data a(1:nf,12) &
      / -0.0027357353070964763,-0.002313408560751684, &
        -0.006594875859577964,-0.0029986295579278455, &
        0.0005436406641105007,-0.0007930995954059661, &
        -0.003988445566829313,-0.005504737136141857, &
        -0.0047408227741938536,-0.002959138508456638, &
        -0.0014990108051475826,-0.0009781124880630099, &
        -0.001327810389155415,-0.0021330397485244124, &
        -0.0029458019595418605,-0.0034587473683276304, &
        -0.0035492118838007514,-0.0032454138885380512, &
        -0.0026631647896235844,-0.0019464520548864161, &
        -0.0012252340952456642,-0.0005936527265095748, &
        -0.00010436774921496195,0.0002266690909606961, &
        0.00041045585547472285,0.0004742471963796735, &
        0.00045198639490429127,0.00037749414408505544, &
        0.000279954263553428,0.00018177692407628232, &
        0.0000980335101594676,0.00003687140950203727, &
        8.621660265717423e-7,-0.000011672476503132856, &
        -5.297826586980747e-6,0.000013745594682812445, &
        0.00003853067813097688,0.00006245571738905387, &
        0.00007974151994749601,0.00008586015594624907 /
!
      data a(1:nf,13) &
      / -0.0028216933336548955,-0.004541787322088764, &
        -0.0011731247151034463,-0.002177822771711968, &
        -0.004356562830916788,-0.004074820577288252, &
        -0.0024929298806266598,-0.0015893652578591322, &
        -0.0018793312522045191,-0.0027875366874603458, &
        -0.003594583498623072,-0.003918600579148723, &
        -0.0037445928472816803,-0.0032664219604783273, &
        -0.002721585936949822,-0.0022892325575789574, &
        -0.002054734617617626,-0.0020182579275327406, &
        -0.002124103461157212,-0.0022920061293874028, &
        -0.002442353177621648,-0.002511670876561789, &
        -0.0024594282464167616,-0.002268289479339103, &
        -0.0019404744305679955,-0.00149242540489114, &
        -0.0009492281596781282,-0.0003398764561750589, &
        0.0003064280012662346,0.0009627422416386144, &
        0.0016057721698065075,0.0022165748703497634, &
        0.0027803981878150585,0.003286262622743341, &
        0.003726270102743068,0.004094926736218141, &
        0.004388185990909066,0.0046026164064695285, &
        0.00473467800986993,0.004778082785167448 /
!
      b2 = b*b
      b25 = b**2.5
      h0 = 0.0
      do i=1,4
       h(i) = b/( y(i)**4 + b2)
      enddo
      do i=5,na
       h(i) = b2/(0.25*y(i)**5 + b25)       
      enddo
! transform to gt grid
      gt = SQRT(1-ft)
      if(gt.gt.g(nf))gt = g(nf)
! find grid position
      if(gt.le.g(2))then
        i=1
      else
        i = INT(gt/0.025)
      endif
      j = i+1
! interpolate the coefficients
      dg = (gt - g(i))/(g(j)-g(i))
      hs = h0
      do k=1,na
        ca = a(i,k) + (a(j,k)-a(i,k))*dg
! sum up the terms
        hs = hs + ca*h(k)
      enddo
! final answer
      FLR_dHw333 = (35.0/9.0)*ft*hs
! 
      END FUNCTION FLR_dHw333
!
      END SUBROUTINE FLR_xgrid

