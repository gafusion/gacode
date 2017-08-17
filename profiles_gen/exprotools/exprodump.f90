!------------------------------------------------------------
! exprodump.f90
!
! PURPOSE:
!  Dump a given EXPRO interface variable as a function of 
!  rmin and rho.
!------------------------------------------------------------

program exprodump

  use EXPRO_interface

  implicit none

  character(len=20) :: var,x
  integer :: i
  real, dimension(5) :: z
  integer :: density_method


  open(unit=1,file='input.exprodump',status='old')
  read(1,*) var
  read(1,*) x
  read(1,*) z(1)
  read(1,*) z(2)
  read(1,*) z(3)
  read(1,*) z(4)
  read(1,*) z(5)
  read(1,*) density_method
  close(1)

  EXPRO_ctrl_n_ion = 5
  EXPRO_ctrl_quasineutral_flag = density_method-1
  EXPRO_ctrl_z(1:5) = z(1:5)
  EXPRO_ctrl_numeq_flag = 0 

  call EXPRO_alloc('./',1) 
  call EXPRO_read

  if (trim(var) == 'list') then

     print *,'rho'
     print *,'rmin'
     print *,'rmaj'
     print *,'q'
     print *,'kappa'

     print *,'delta'
     print *,'te'
     print *,'ne'
     print *,'z_eff'
     print *,'w0'

     print *,'flow_mom'
     print *,'pow_e'
     print *,'pow_i'
     print *,'pow_ei'
     print *,'zeta'

     print *,'flow_beam'
     print *,'flow_wall'
     print *,'zmag'
     print *,'ptot'
     print *,'polflux'

     print *,'ni_1'
     print *,'ni_2'
     print *,'ni_3'
     print *,'ni_4'
     print *,'ni_5'

     print *,'Ti_1'
     print *,'Ti_2'
     print *,'Ti_3'
     print *,'Ti_4'
     print *,'Ti_5'

     print *,'vtor_1'
     print *,'vtor_2'
     print *,'vtor_3'
     print *,'vtor_4'
     print *,'vtor_5'

     print *,'vpol_1'
     print *,'vpol_2'
     print *,'vpol_3'
     print *,'vpol_4'
     print *,'vpol_5'

     print *,'bunit'
     print *,'s'
     print *,'drmaj'
     print *,'dzmag'
     print *,'sdelta'

     print *,'vol'
     print *,'volp'
     print *,'cs'
     print *,'rhos'
     print *,'drdrho'

     print *,'grad_r0'
     print *,'ave_grad_r'

  else   

     do i=1,EXPRO_n_exp

        ! Column 1
        select case (trim(x)) 
        case ('none')
           write(*,'(a)',advance='no') ' ' 
        case ('r')
           write(*,10,advance='no') EXPRO_rmin(i)/EXPRO_rmin(EXPRO_n_exp)
        case ('rho')
           write(*,10,advance='no') EXPRO_rho(i)
        case ('psi')
           write(*,10,advance='no') EXPRO_polflux(i)/EXPRO_polflux(EXPRO_n_exp)
        end select

        ! Column 2
        select case (trim(var))

        case ('rho')
           print 20,EXPRO_rho(i)
        case ('rmin')
           print 20,EXPRO_rmin(i)
        case ('rmaj')
           print 20,EXPRO_rmaj(i)
        case ('q')
           print 20,EXPRO_q(i)
        case ('kappa')
           print 20,EXPRO_kappa(i)
        case ('delta')
           print 20,EXPRO_delta(i)
        case ('te')
           print 20,EXPRO_te(i)
        case ('ne')
           print 20,EXPRO_ne(i)
        case ('z_eff')
           print 20,EXPRO_z_eff(i)
        case ('w0')
           print 20,EXPRO_w0(i)
        case ('flow_mom')
           print 20,EXPRO_flow_mom(i)
        case ('pow_e')
           print 20,EXPRO_pow_e(i)
        case ('pow_i')
           print 20,EXPRO_pow_i(i)
        case ('pow_ei')
           print 20,EXPRO_pow_ei(i)
        case ('zeta')
           print 20,EXPRO_zeta(i)
        case ('flow_beam')
           print 20,EXPRO_flow_beam(i)
        case ('flow_wall')
           print 20,EXPRO_flow_wall(i)
        case ('zmag')
           print 20,EXPRO_zmag(i)
        case ('ptot')
           print 20,EXPRO_ptot(i)
        case ('polflux')
           print 20,EXPRO_polflux(i)
        case ('ni_1')
           if (EXPRO_ctrl_quasineutral_flag == 1) then
              print 20,EXPRO_ni_new(i)
           else
              print 20,EXPRO_ni(1,i)
           endif
        case ('ni_2')
           print 20,EXPRO_ni(2,i)
        case ('ni_3')
           print 20,EXPRO_ni(3,i)
        case ('ni_4')
           print 20,EXPRO_ni(4,i)
        case ('ni_5')
           print 20,EXPRO_ni(5,i)
        case ('Ti_1')
           print 20,EXPRO_Ti(1,i)
        case ('Ti_2')
           print 20,EXPRO_Ti(2,i)
        case ('Ti_3')
           print 20,EXPRO_Ti(3,i)
        case ('Ti_4')
           print 20,EXPRO_Ti(4,i)
        case ('Ti_5')
           print 20,EXPRO_Ti(5,i)   
        case ('vtor_1')
           print 20,EXPRO_vtor(1,i)
        case ('vtor_2')
           print 20,EXPRO_vtor(2,i)
        case ('vtor_3')
           print 20,EXPRO_vtor(3,i)
        case ('vtor_4')
           print 20,EXPRO_vtor(4,i)
        case ('vtor_5')
           print 20,EXPRO_vtor(5,i)   
        case ('vpol_1')
           print 20,EXPRO_vpol(1,i)
        case ('vpol_2')
           print 20,EXPRO_vpol(2,i)
        case ('vpol_3')
           print 20,EXPRO_vpol(3,i)
        case ('vpol_4')
           print 20,EXPRO_vpol(4,i)
        case ('vpol_5')
           print 20,EXPRO_vpol(5,i)
        case ('bunit')
           print 20,EXPRO_bunit(i)
        case ('s')
           print 20,EXPRO_s(i)
        case ('drmaj')
           print 20,EXPRO_drmaj(i)
        case ('dzmag')
           print 20,EXPRO_dzmag(i)
        case ('sdelta')
           print 20,EXPRO_sdelta(i)
        case ('vol')
           print 20,EXPRO_vol(i)
        case ('volp')
           print 20,EXPRO_volp(i)
        case ('cs')
           print 20,EXPRO_cs(i)
        case ('rhos')
           print 20,EXPRO_rhos(i)
        case ('drdrho')
           print 20,EXPRO_drdrho(i)
        case ('grad_r0')
           print 20,EXPRO_grad_r0(i)
        case ('ave_grad_r')
           print 20,EXPRO_ave_grad_r(i)
        case default
           print 20,0.0
        end select

     enddo

  endif

  call EXPRO_alloc('./',0) 

10 format(t2,1pe12.5)
20 format(1x,1pe12.5)

end program exprodump
