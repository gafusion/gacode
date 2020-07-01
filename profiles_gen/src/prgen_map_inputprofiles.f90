!------------------------------------------------------------
! prgen_map_inputprofiles.f90
!
! PURPOSE:
!  Map input profiles plus gfile data.
!------------------------------------------------------------

subroutine prgen_map_inputprofiles

  use prgen_globals
  use expro

  implicit none

  integer :: i
  real :: xnew(nx)
  
  if (efit_method == 1) then

     call prgen_get_chi(nx,q,dpsi,rho,torfluxa)

     ! Align with new rho grid, not old rho grid, expro_rho:
     call cub_spline(expro_rho,expro_te,nx,rho,xnew,nx)
     expro_te = xnew
     call cub_spline(expro_rho,expro_ne,nx,rho,xnew,nx)
     expro_ne = xnew
     call cub_spline(expro_rho,expro_z_eff,nx,rho,xnew,nx)
     expro_z_eff = xnew
     call cub_spline(expro_rho,expro_w0,nx,rho,xnew,nx)
     expro_w0 = xnew
     call cub_spline(expro_rho,expro_ptot,nx,rho,xnew,nx)
     expro_ptot = xnew
     do i=1,10
        call cub_spline(expro_rho,expro_ni(i,:),nx,rho,xnew,nx)
        expro_ni(i,:) = xnew(:)
        call cub_spline(expro_rho,expro_ti(i,:),nx,rho,xnew,nx)
        expro_ti(i,:) = xnew(:)
        call cub_spline(expro_rho,expro_vtor(i,:),nx,rho,xnew,nx)
        expro_vtor(i,:) = xnew(:)
        call cub_spline(expro_rho,expro_vpol(i,:),nx,rho,xnew,nx)
        expro_vpol(i,:) = xnew(:)
     enddo

     call cub_spline(expro_rho,expro_qohme,nx,rho,xnew,nx)  ; expro_qohme = xnew
     call cub_spline(expro_rho,expro_qbeame,nx,rho,xnew,nx) ; expro_qbeame = xnew
     call cub_spline(expro_rho,expro_qbeami,nx,rho,xnew,nx) ; expro_qbeami = xnew
     call cub_spline(expro_rho,expro_qrfe,nx,rho,xnew,nx)   ; expro_qrfe = xnew
     call cub_spline(expro_rho,expro_qrfi,nx,rho,xnew,nx)   ; expro_qrfi = xnew
     call cub_spline(expro_rho,expro_qfuse,nx,rho,xnew,nx) ; expro_qfuse = xnew
     call cub_spline(expro_rho,expro_qfusi,nx,rho,xnew,nx) ; expro_qfusi = xnew
     call cub_spline(expro_rho,expro_qbrem,nx,rho,xnew,nx) ; expro_qbrem = xnew
     call cub_spline(expro_rho,expro_qsync,nx,rho,xnew,nx) ; expro_qsync = xnew
     call cub_spline(expro_rho,expro_qline,nx,rho,xnew,nx) ; expro_qline = xnew
     call cub_spline(expro_rho,expro_qei,nx,rho,xnew,nx)   ; expro_qei = xnew
     call cub_spline(expro_rho,expro_qione,nx,rho,xnew,nx) ; expro_qione = xnew
     call cub_spline(expro_rho,expro_qioni,nx,rho,xnew,nx) ; expro_qioni = xnew
     call cub_spline(expro_rho,expro_qcxi,nx,rho,xnew,nx)  ; expro_qcxi = xnew
     call cub_spline(expro_rho,expro_qmom,nx,rho,xnew,nx)  ; expro_qmom = xnew
     call cub_spline(expro_rho,expro_qpar_wall,nx,rho,xnew,nx)  ; expro_qpar_wall = xnew
     call cub_spline(expro_rho,expro_qpar_beam,nx,rho,xnew,nx)  ; expro_qpar_beam = xnew

     !---------------------------------------------------------
     ! Map profile data into expro interface.
     ! NOTE: expro_alloc already called in prgen_read_inputprofiles
     !
     expro_rho(:)  = rho
     expro_rmin(:) = rmin
     expro_rmaj(:) = rmaj
     !---------------------------------------------------------
 
  endif

  !---------------------------------------------------------
  ! Notification about CER file
  !
  if (file_cer /= "null") then
    print '(a)','INFO: (prgen_map_inputprofiles) IGNORING cer file'
  endif
  !---------------------------------------------------------

end subroutine prgen_map_inputprofiles

