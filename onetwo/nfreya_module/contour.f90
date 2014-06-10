   MODULE CONTOUR
! --- INCLUDE file contour.i
!
!
! --- rplasbdry, zplasbdry holds the given fixed boundary
! --- for the fixed boundary equilibrium calculations
!
      integer ,parameter :: nconmax = 2000, nrplas = 1500 ! see kstore (storage.i)
!
      integer *4 ::        nplasbdry, ncontr, nboundpts
      real *8, dimension(:) ::                                   &
                      rcontr(nconmax), zcontr(nconmax),          &
                      rplasbdry(nrplas), zplasbdry(nrplas),      &
                      boundrcntr(nrplas), boundzcntr(nrplas)
      real *8  rplasmin, rplasmax, zplasmin, zplasmax,rsep, zsep
!
   END MODULE CONTOUR
