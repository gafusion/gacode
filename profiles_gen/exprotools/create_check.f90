subroutine create_check

  use create_globals

  implicit none

  integer :: i

  if (nx < 2) then
     print *, 'ERROR: (create) User must set NX > 1.'
     stop
  endif

  select case(set_exm_b_ref)
  case(0)
  case(1)
  case default
     print *, 'ERROR: SET parameters must be 0 or 1'
     stop
  end select

  select case(set_exm_arho)
  case(0)
  case(1)
  case default
     print *, 'ERROR: SET parameters must be 0 or 1'
     stop
  end select

  select case(set_exm_rho)
  case(0)
  case(1)
  case default
     print *, 'ERROR: SET parameters must be 0 or 1'
     stop
  end select

  select case(set_exm_rmin)
  case(0)
  case(1)
     if(exm_a_meters <= 0.0) then
        print *, 'ERROR: A_METERS must be positive'
     endif
  case default
     print *, 'ERROR: SET parameters must be 0 or 1'
     stop
  end select

  select case(set_exm_rmaj)
  case(0)
  case(1)
     if(exm_rmaj < 0.0) then
        print *, 'ERROR: RMAJ must be positive'
        stop
     endif
  case default
     print *, 'ERROR: SET parameters must be 0 or 1'
     stop
  end select

  select case(set_exm_q)
  case(0)
  case(1)
  case default  
     print *, 'ERROR: SET parameters must be 0 or 1'
     stop
  end select

  select case(set_exm_kappa)
  case(0)
  case(1)
     if(exm_kappa < 1.0) then
        print *, 'ERROR: KAPPA must be >= 1'
        stop
     endif
  case default  
     print *, 'ERROR: SET parameters must be 0 or 1'
     stop
  end select

  select case(set_exm_delta)
  case(0)
  case(1)
  case default  
     print *, 'ERROR: SET parameters must be 0 or 1'
     stop
  end select

  select case(set_exm_te)
  case(0)
  case(1)
     if(exm_te_axis < 0.0) then
        print *, 'ERROR: TE_AXIS must be positive'
        stop
     endif
  case default  
     print *, 'ERROR: SET parameters must be 0 or 1'
     stop
  end select

  select case(set_exm_z_eff)
  case(0)
  case(1)
     if(exm_z_eff < 1.0) then
        print *, 'ERROR: Z_EFF must be >= 1.0'
        stop
     endif
  case default  
     print *, 'ERROR: SET parameters must be 0 or 1'
     stop
  end select

  select case(set_exm_w0)
  case(0)
  case(1)
  case default  
     print *, 'ERROR: SET parameters must be 0 or 1'
     stop
  end select

  select case(set_exm_zeta)
  case(0)
  case(1)
  case default  
     print *, 'ERROR: SET parameters must be 0 or 1'
     stop
  end select

  select case(set_exm_zmag)
  case(0)
  case(1)
  case default  
     print *, 'ERROR: SET parameters must be 0 or 1'
     stop
  end select

  do i=1,nions_max
     select case(set_exm_ti(i))
     case(0)
     case(1)
        if(exm_ti_axis(i) < 0.0) then
           print *, 'ERROR: TI_AXIS must be positive'
           stop
        endif
     case default  
        print *, 'ERROR: SET parameters must be 0 or 1'
        stop
     end select
  end do

end subroutine create_check
