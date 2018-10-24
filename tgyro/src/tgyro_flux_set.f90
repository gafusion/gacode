!-----------------------------------------------------------------------
! tgyro_flux_set.f90
!
! PURPOSE:
!  Manage common operation of mapping fluxes to internal array.
!-----------------------------------------------------------------------

subroutine tgyro_flux_set(f_vec)

  use tgyro_globals

  implicit none

  integer :: i
  integer :: p
  integer :: i_ion

  real, intent(inout), dimension(p_max) :: f_vec

  p = 0
  do i=2,n_r

     if (loc_ti_feedback_flag == 1) then
        p = p+1
        f_vec(p) = eflux_i_tot(i)
     endif

     if (loc_te_feedback_flag == 1) then
        p = p+1
        f_vec(p) = eflux_e_tot(i)
     endif

     if (loc_er_feedback_flag == 1) then
        p = p+1
        f_vec(p) = mflux_tot(i)
     endif

     if (evo_e(0) == 1) then
        p = p+1
        f_vec(p) = pflux_e_tot(i)
     endif

     do i_ion=1,loc_n_ion
        if (evo_e(i_ion) >= 1) then
           p = p+1
           f_vec(p) = pflux_i_tot(i_ion,i)
        endif
     enddo

  enddo

end subroutine tgyro_flux_set
