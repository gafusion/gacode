!-----------------------------------------------------------------
! cgyro_write_timedata.f90
!
! PURPOSE:
!  Output of time-independent (initial) data.
!-----------------------------------------------------------------

subroutine cgyro_write_initdata

  use mpi
  use cgyro_globals
  use cgyro_experimental_globals

  implicit none

  integer :: p,in,is
  real :: kymax,z_eff
  real, external ::spectraldiss

  !----------------------------------------------------------------------------
  ! Runfile to give complete summary to user
  ! 
  ! NOTE: Some of this data is reproduced in other data files.
  !
  if (silent_flag == 0 .and. i_proc == 0) then

     open(unit=io,file=trim(path)//runfile_info,status='old',position='append')


     write(io,*)
     write(io,'(a)') ' n_theta | n_species | n_energy | n_xi '
     write(io,'(t4,i3,t16,i1,t26,i2,t36,i2)') n_theta,n_species,n_energy,n_xi
     if (test_flag == 0) then
        write(io,*) 
        write(io,'(a)') ' nc_loc | nv_loc | nsplit | n_MPI | n_OMP'
        write(io,'(t3,i4,t12,i4,t21,i4,t29,i4,t36,i4)') nc_loc,nv_loc,nsplit,n_proc,n_omp
     endif

     if (zf_test_flag == 0) then

        ! Compute kymax
        if (n_toroidal == 1) then
           kymax = q/rmin*rho
        else
           kymax = q/rmin*(n_toroidal-1)*rho
        endif

        write(io,*)
        write(io,*) '          n    Delta      Max     L/rho'
        write(io,'(a,i4,2x,2(f7.3,2x),2x,f6.2)') ' kx*rho:',&
             n_radial,2*pi*rho/length,2*pi*rho*(n_radial/2-1)/length,length/rho
        write(io,'(a,i4,2x,2(f7.3,2x),2x,f6.2)') ' ky*rho:',&
             n_toroidal,q/rmin*rho,kymax,2*pi/ky

     else
        if (n_radial==1) then
           write(io,*) ' kx*rho:',2*pi*rho/length
        else
           write(io,*) '          n          Delta            Max           L/rho'
           write(io,'(a,i4,2x,2(g0.8,2x),2x,g0.8)') ' kx*rho:',&
                n_radial,2*pi*rho/length,2*pi*rho*(n_radial/2-1)/length,length/rho
        endif

     endif

     ! Compute Z_eff (diagnostic only)
     z_eff = 0.0
     do is=1,n_species
        if (z(is) > 0.0) then 
           z_eff = z_eff+dens(is)*z(is)**2/dens_ele
        endif
     enddo

     ! Dissipation information
     write(io,*)
     write(io,'(t2,3(a,1x,i2,3x),t48,a)') &
          ' D-theta:',2*nup_theta, ' D-radial:',2*nup_radial,' D-alpha:',2*nup_alpha,&
          '[dissipation order]'
     write(io,'(t2,3(a,1x,f3.1,2x),t48,a)') &
          'up_theta:',up_theta , 'up_radial:',up_radial ,'up_alpha:',up_alpha,&
          '[dissipation strength]'

     write(io,*) 
     write(io,20) '    r/a:',rmin
     write(io,20) '    R/a:',rmaj, '  shift:',shift,  '  betae:',betae_unit
     write(io,20) '      q:',q,    '      s:',s,      ' beta_*:',beta_star(0)
     write(io,20) '  kappa:',kappa,'s_kappa:',s_kappa,' lamb_*:',lambda_star
     write(io,20) '  delta:',delta,'s_delta:',s_delta,'gamma_e:',gamma_e
     write(io,20) '   zeta:',zeta, ' s_zeta:',s_zeta, 'gamma_p:',gamma_p
     write(io,20) '   zmag:',zmag, '  dzmag:',dzmag,  '   mach:',mach

     write(io,*)
     write(io,20) '[rho/a]:',rho,'[z_eff]:',z_eff,'[w_E*dt]',(n_toroidal-1)*q/rmin*length*gamma_e/(2*pi)*delta_t
     write(io,*)
     write(io,'(a)') &
          ' i  z  n/n_norm   T/T_norm   m/m_norm     a/Ln       a/Lt       nu         s_n        s_t'
     do is=1,n_species
        write(io,'(t1,i2,1x,i2,3(2x,1pe9.3),2(1x,1pe10.3),(2x,1pe9.3),2(1x,1pe10.3))') &
             is,int(z(is)),dens(is),temp(is),mass(is),dlnndr(is),dlntdr(is),nu(is),sdlnndr(is),sdlntdr(is)
     enddo

     if (profile_model == 2) then
        write(io,*)
        write(io,10) '           a[m]:',a_meters, '  b_unit[T]:',b_unit
        write(io,10) 'n_norm[e19/m^3]:',dens_norm,'v_norm[m/s]:',vth_norm,'T_norm[keV]:',temp_norm
     endif
     write(io,*)
    
     close(io)

  endif
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  ! Write the initial equilibrium data
  !
  if (silent_flag == 0 .and. i_proc == 0) then

     open(unit=io,file=trim(path)//'out.cgyro.equilibrium',status='replace')
     write (io,fmtstr) rmin
     write (io,fmtstr) rmaj
     write (io,fmtstr) q
     write (io,fmtstr) s
     write (io,fmtstr) shift
     write (io,fmtstr) kappa
     write (io,fmtstr) s_kappa
     write (io,fmtstr) delta
     write (io,fmtstr) s_delta
     write (io,fmtstr) zeta
     write (io,fmtstr) s_zeta
     write (io,fmtstr) zmag
     write (io,fmtstr) dzmag
     write (io,fmtstr) rho
     write (io,fmtstr) ky
     write (io,fmtstr) betae_unit
     write (io,fmtstr) beta_star(0)
     write (io,fmtstr) lambda_star
     write (io,fmtstr) gamma_e
     write (io,fmtstr) gamma_p
     write (io,fmtstr) mach
     write (io,fmtstr) a_meters
     write (io,fmtstr) b_unit
     write (io,fmtstr) dens_norm
     write (io,fmtstr) temp_norm
     write (io,fmtstr) vth_norm
     do is=1,n_species
        write (io,fmtstr) z(is)
        write (io,fmtstr) mass(is)
        write (io,fmtstr) dens(is)
        write (io,fmtstr) temp(is)
        write (io,fmtstr) dlnndr(is)
        write (io,fmtstr) dlntdr(is)
        write (io,fmtstr) nu(is)
     enddo
     close(io)

  endif

  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  ! Write the initial equilibrium geometry data
  !
  if (silent_flag == 0 .and. i_proc == 0) then

     open(unit=io,file=trim(path)//'out.cgyro.geo',status='replace')
     write(io,fmtstr) theta(:)
     write(io,fmtstr) w_theta(:)
     write(io,fmtstr) bmag(:)
     write(io,fmtstr) omega_stream(:,1)
     write(io,fmtstr) omega_trap(:,1)
     write(io,fmtstr) omega_rdrift(:,1)
     write(io,fmtstr) omega_adrift(:,1)
     write(io,fmtstr) omega_aprdrift(:,1)
     write(io,fmtstr) omega_cdrift(:,1)
     write(io,fmtstr) omega_cdrift_r(:,1)
     write(io,fmtstr) omega_gammap(:)
     write(io,fmtstr) k_perp(ic_c(n_radial/2+1,:))
     close(io)

  endif
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  ! Write the initial grid data 
  !
  if (silent_flag == 0 .and. i_proc == 0) then

     open(unit=io,file=trim(path)//runfile_grids,status='replace')
     write(io,'(i4)') n_toroidal
     write(io,'(i4)') n_species
     write(io,'(i4)') n_field
     write(io,'(i4)') n_radial
     write(io,'(i4)') n_theta
     write(io,'(i4)') n_energy
     write(io,'(i4)') n_xi
     write(io,'(i4)') box_size
     write(io,'(1pe12.5)') length/rho
     write(io,'(i4)') n_global
     write(io,'(i4)') px(:)
     write(io,'(1pe12.5)') theta(:)
     write(io,'(1pe12.5)') energy(:)
     write(io,'(1pe12.5)') xi(:)
     write(io,'(1pe12.5)') transpose(thetab(:,:))
     write(io,'(1pe12.5)') (rho*q/rmin*in,in=0,n_toroidal-1)
     write(io,'(1pe12.5)') (spectraldiss((pi/n_toroidal)*in,nup_alpha),in=0,n_toroidal-1)
     write(io,'(1pe12.5)') (spectraldiss((2*pi/n_radial)*p,nup_radial),p=-n_radial/2,n_radial/2-1)
     close(io)

  endif
  !----------------------------------------------------------------------------

10 format(t2,4(a,1x,1pe9.3,2x))  
20 format(t2,4(a,1x,1pe10.3,2x)) 

end subroutine cgyro_write_initdata
