   MODULE zonal_data

      USE nrtype,                        ONLY : Dp,I4B

      IMPLICIT NONE

      REAL(DP),DIMENSION(:),  ALLOCATABLE :: wnoperm,                     &
                                             zone_volume,                 &
                                             zone_area,potsid,rinsid,     &
                                             rotsid,pinsid,b1ins,b2ins,   &
                                             b1ots,b2ots,fpsio,fpsii,zne, &
                                             zte,zti,zangrot,zenbeam,     &
                                             zenbeamold,psif,psivol,      &
                                             rho_zone

      REAL(DP),DIMENSION(:,:),ALLOCATABLE :: zni,zzi,r_zone_cntr,         &
                                             z_zone_cntr

      INTEGER(I4B) nw,nh,nwh,nh2,nwork,nplasbdry,nlimiter,ncontour,       &
                   mf, mfm1,max_ctr_pts
      INTEGER(I4B),DIMENSION(:), ALLOCATABLE :: n_zone_cntr


   END MODULE zonal_data
