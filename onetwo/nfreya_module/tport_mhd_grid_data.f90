   MODULE tport_mhd_grid_data

      USE nrtype,                        ONLY : Dp,I4B

      IMPLICIT NONE

      REAL(DP),DIMENSION(:),  ALLOCATABLE :: psir,  psival,rmhdgrid,       &
                                             zmhdgrid,xlimiter,ylimiter,   &
                                             rplasbdry, zplasbdry, bpmag,  &
                                             rcontour,zcontour,fpsi,ene,te,&
                                             ti,angrot,psivolp      

      !rcontour(ncontour),zcontour(ncontour) hold values for last flux surface traced

      REAL(DP),DIMENSION(:,:),ALLOCATABLE :: psi_g_cm2,en  ! psi in gauss cm**2

      REAL(DP),DIMENSION(:,:,:),ALLOCATABLE :: hibr,hdep,angmpf 

      REAL(DP) xmagn1,ymagn1,rcmin,rcmax,zcmin,                           &
               zcmax,zrcmin,zrcmax,rzcmin,rzcmax,volume,                  &
               xmin_lim,xmax_lim,ymin_lim,ymax_lim,                       &
               psiax,psilim,                                              &
               rplasmin,rplasmax,zplasmin,zplasmax

      INTEGER(I4B) ncontour

   END MODULE tport_mhd_grid_data
