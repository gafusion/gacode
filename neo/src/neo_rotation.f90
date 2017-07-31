module neo_rotation

  implicit none

  public :: ROT_alloc, ROT_solve_phi, ROT_write
  real, dimension(:), allocatable   :: phi_rot        ! phi(r)-phi(r,theta0)
  real, dimension(:), allocatable   :: phi_rot_deriv  ! theta derivative
  real, dimension(:), allocatable   :: phi_rot_rderiv ! radial derivative
  real                              :: phi_rot_avg    ! flux surface avg
  real, dimension(:,:), allocatable :: dens_fac       ! (ns,nth):
                                                      ! n=n(th0)*dens_fac
  real, dimension(:,:), allocatable :: lam_rot_aniso  ! lambda_aniso
  real, dimension(:,:), allocatable :: lam_rot_kpar_aniso ! bhatdotgrad lam
  real, dimension(:), allocatable   :: lam_rot_avg_aniso ! <lambda_aniso>
  real, dimension(:), allocatable   :: rotavg_e0, rotavg_e1, &
       rotavg_e2, rotavg_e3, rotavg_e4, rotavg_e5a, rotavg_e5b, rotavg_e6

  logical, private :: initialized = .false.
  integer, private :: io_rot=50
  character(len=80),private :: runfile = 'out.neo.rotation'

  contains
    subroutine ROT_alloc(flag)
      use neo_globals
      implicit none
      integer, intent (in) :: flag  ! flag=1: allocate; else deallocate

      if(flag == 1) then
         if(initialized) return
         allocate(phi_rot(n_theta))
         allocate(phi_rot_deriv(n_theta))
         allocate(phi_rot_rderiv(n_theta))
         allocate(dens_fac(n_species,n_theta))
         allocate(lam_rot_aniso(n_species,n_theta))
         allocate(lam_rot_kpar_aniso(n_species,n_theta))
         allocate(lam_rot_avg_aniso(n_species))
         allocate(rotavg_e0(n_species))
         allocate(rotavg_e1(n_species))
         allocate(rotavg_e2(n_species))
         allocate(rotavg_e3(n_species))
         allocate(rotavg_e4(n_species))
         allocate(rotavg_e5a(n_species))
         allocate(rotavg_e5b(n_species))
         allocate(rotavg_e6(n_species))
         if(silent_flag == 0 .and. i_proc == 0 .and. rotation_model == 2) then
            open(unit=io_rot,file=trim(path)//runfile,status='replace')
            close(io_rot)
         end if
         initialized = .true.

      else
         if(.NOT. initialized) return
         deallocate(phi_rot)
         deallocate(phi_rot_deriv)
         deallocate(phi_rot_rderiv)
         deallocate(dens_fac)
         deallocate(lam_rot_aniso)
         deallocate(lam_rot_kpar_aniso)
         deallocate(lam_rot_avg_aniso)
         deallocate(rotavg_e0)
         deallocate(rotavg_e1)
         deallocate(rotavg_e2)
         deallocate(rotavg_e3)
         deallocate(rotavg_e4)
         deallocate(rotavg_e5a)
         deallocate(rotavg_e5b)
         deallocate(rotavg_e6)
         initialized = .false.

      endif

    end subroutine ROT_alloc

    ! solve the quasi-neutrality relation for poloidal part of phi
    subroutine ROT_solve_phi(ir)
      use neo_globals
      use neo_equilibrium
      implicit none
      integer, intent (in) :: ir
      integer, parameter :: nmax = 200
      integer :: it, is, jt, id, n, is_ele
      real :: x, x0, sum_zn, dsum_zn, fac
      real, dimension(:,:), allocatable :: dlnndr_fac, dens_fac_new, &
           lam_rot_rderiv_aniso
      real, dimension(:), allocatable :: dens_avg, dens_avg_new
      real :: dens_new

      if(rotation_model == 1 .or. spitzer_model==1) then
         phi_rot(:) = 0.0
         phi_rot_deriv(:) = 0.0
         phi_rot_rderiv(:) = 0.0
         dens_fac(:,:) = 1.0
         lam_rot_aniso(:,:) = 0.0
         lam_rot_kpar_aniso(:,:) = 0.0
         lam_rot_avg_aniso(:) = 0.0
         phi_rot_avg = 0.0
         omega_rot(ir) = 0.0
         omega_rot_deriv(ir) = 0.0
         rotavg_e0(:)  = 0.0
         rotavg_e1(:)  = 0.0
         rotavg_e2(:)  = 0.0
         rotavg_e3(:)  = 0.0
         rotavg_e4(:)  = 0.0
         rotavg_e5a(:) = 0.0
         rotavg_e5b(:) = 0.0
         rotavg_e6(:)  = 0.0

      else 
         
         ! check for equilibrium-scale QN at theta=0
         sum_zn  = 0.0
         dsum_zn = 0.0
         do is=1, n_species
            sum_zn  = sum_zn  + z(is) * dens(is,ir)
            dsum_zn = dsum_zn + z(is) * dens(is,ir) * dlnndr(is,ir)
         enddo
         if(adiabatic_ele_model == 1) then
            sum_zn  = sum_zn - ne_ade(ir)
            sum_zn  = sum_zn / ne_ade(ir)
            dsum_zn = dsum_zn - ne_ade(ir) * dlnndre_ade(ir)
            dsum_zn = dsum_zn / ne_ade(ir)
         else
            do is=1, n_species
               if(Z(is) < 0.0) then
                  is_ele = is
                  exit
               endif
            enddo
            sum_zn   = sum_zn / dens(is_ele,ir)
            dsum_zn  = dsum_zn / dens(is_ele,ir)
         endif

         if(abs(sum_zn) > 1.0e-3) then
            if(silent_flag == 0 .and. i_proc == 0) then
               open(unit=io_neoout,file=trim(path)//runfile_neoout,&
                    status='old',position='append')
               write(io_neoout,*) &
                    'WARNING: (NEO) ROTATION IS BEING RUN WITHOUT QN DENSITIES'
               close(io_neoout)
            endif
         endif

         allocate(dens_avg(n_species))
         allocate(dlnndr_fac(n_species,n_theta))
         allocate(dens_avg_new(n_species))
         allocate(dens_fac_new(n_species,n_theta))
         allocate(lam_rot_rderiv_aniso(n_species,n_theta))

         lam_rot_aniso(:,:) = 0.0
         lam_rot_kpar_aniso(:,:) = 0.0
         lam_rot_avg_aniso(:) = 0.0
         lam_rot_rderiv_aniso(:,:) = 0.0
         do is=1,n_species
            if(aniso_model(is) == 2) then
               lam_rot_avg_aniso(is) = 0.0
               do it=1,n_theta
                  !!!!! first model
                  !lam_rot_aniso(is,it) = &
                  !     (temp_perp(is,ir)/temp_para(is,ir) - 1.0)**0.75 &
                  !     * log(Bmag(it)/Bmag_th0)
                  !lam_rot_kpar_aniso(is,it) = &
                  !     (temp_perp(is,ir)/temp_para(is,ir) - 1.0)**0.75 &
                  !     * gradpar_Bmag(it) / Bmag(it)
                  !lam_rot_rderiv_aniso(is,it) = log(Bmag(it)/Bmag_th0) &
                  !     * 0.75 &
                  !     / (temp_perp(is,ir)/temp_para(is,ir) - 1.0)**0.25 &
                  !     * temp_perp(is,ir)/temp_para(is,ir) &
                  !     * (-dlntdr_perp(is,ir) + dlntdr_para(is,ir)) &
                  !     + (temp_perp(is,ir)/temp_para(is,ir) - 1.0)**0.75 &
                  !     * (1.0/Bmag(it) * Bmag_rderiv(it) &
                  !     - 1.0/Bmag_th0 * Bmag_th0_rderiv)

                  !!!!! second model
                  lam_rot_aniso(is,it) = &
                       log(temp_perp(is,ir)/temp_para(is,ir) &
                       - Bmag_th0/Bmag(it) &
                       * (temp_perp(is,ir)/temp_para(is,ir) - 1.0))
                  lam_rot_kpar_aniso(is,it) = &
                       1.0/(temp_perp(is,ir)/temp_para(is,ir) &
                       - Bmag_th0/Bmag(it) &
                       * (temp_perp(is,ir)/temp_para(is,ir) - 1.0)) &
                       * (temp_perp(is,ir)/temp_para(is,ir) - 1.0) &
                       *  Bmag_th0/Bmag(it)**2 * gradpar_Bmag(it)
                  lam_rot_rderiv_aniso(is,it) = &
                       1.0/(temp_perp(is,ir)/temp_para(is,ir) &
                       - Bmag_th0/Bmag(it) &
                       * (temp_perp(is,ir)/temp_para(is,ir) - 1.0)) &
                       * ( (1.0 - Bmag_th0/Bmag(it)) &
                       * temp_perp(is,ir)/temp_para(is,ir) &
                       * (-dlntdr_perp(is,ir) + dlntdr_para(is,ir)) &
                       - (temp_perp(is,ir)/temp_para(is,ir) - 1.0) &
                       * (1.0/Bmag_th0 * Bmag_th0_rderiv &
                       - Bmag_th0/Bmag(it)**2 * Bmag_rderiv(it)) )

                  lam_rot_avg_aniso(is) = lam_rot_avg_aniso(is) &
                       + w_theta(it) * lam_rot_aniso(is,it)
               enddo
            endif
         enddo

         ! partial component of n/n(theta0) 
         ! -- will add the phi_rot component after the QN solve
         do is=1, n_species
            do it=1,n_theta
               dens_fac(is,it) = &
                    exp(omega_rot(ir)**2 * 0.5/vth_para(is,ir)**2 &
                    * (bigR(it)**2 - bigR_th0**2) - lam_rot_aniso(is,it)) 
            enddo
         enddo

         phi_rot_avg = 0.0
         x = 0.05          ! initial guess for phi_rot(1)
         do it=1,n_theta
            
            n=1

            do 
               ! use Newton's method to solve the quasi-neutrality relation
               
               sum_zn  = 0.0
               dsum_zn = 0.0
               do is=1, n_species
                  fac = z(is) * dens(is,ir) * dens_fac(is,it) &
                       * exp(-z(is) / temp_para(is,ir) * x)
                  sum_zn  = sum_zn  + fac
                  dsum_zn = dsum_zn - z(is) / temp_para(is,ir) * fac
               enddo
            
               if(adiabatic_ele_model == 1) then
                  fac = -ne_ade(ir) * exp(1.0/te_ade(ir) * x)
                  sum_zn  = sum_zn  + fac
                  dsum_zn = dsum_zn + 1.0 / te_ade(ir) * fac
               endif

               x0 = x
               x  = x0 - sum_zn / dsum_zn
               
               if (abs(x-x0) < 1.0e-12) exit
               
               n = n + 1
               if(n > nmax) exit
               
            enddo
            
            if(n > nmax) then
               call neo_error('ERROR: (NEO) Rotation density computation failed to converge')
               return
            endif
            
            ! phi_rot is total phi(r) - total phi(r,theta0)
            phi_rot(it) = x  
            
            ! flux surface avg
            phi_rot_avg = phi_rot_avg &
                 + w_theta(it) * phi_rot(it)  
            
         enddo

         ! n(theta)/n(0)
         do is=1, n_species
            do it=1,n_theta
               dens_fac(is,it) = dens_fac(is,it) &
                    * exp(-z(is) / temp_para(is,ir) * phi_rot(it))
               
            enddo
         enddo

         ! <n>
         dens_avg(:)   = 0.0
         do is=1, n_species
            do it=1,n_theta
               dens_avg(is) = dens_avg(is) + &
                    w_theta(it) * dens(is,ir) * dens_fac(is,it)
            enddo
         enddo

         ! dln n(theta)/dr
         do is=1, n_species
            do it=1,n_theta
               dlnndr_fac(is,it) = -(-dlnndr(is,ir) + omega_rot_deriv(ir) &
                    * omega_rot(ir)/vth_para(is,ir)**2 &
                    * (bigR(it)**2 - bigR_th0**2) &
                    + omega_rot(ir)**2 / vth_para(is,ir)**2 &
                    * (bigR(it)*bigR_rderiv(it) &
                    - bigR_th0*bigR_th0_rderiv) &
                    - phi_rot(it) * z(is)/temp_para(is,ir) &
                    * dlntdr_para(is,ir) &
                    + 0.5 * omega_rot(ir)**2 / vth_para(is,ir)**2 &
                    * dlntdr_para(is,ir) * (bigR(it)**2 - bigR_th0**2) &
                    - lam_rot_rderiv_aniso(is,it))
            enddo
         enddo

         ! solve the radial derivative of QN for dphi*/dr
         do it=1,n_theta
            phi_rot_rderiv(it) = 0.0
            sum_zn  = 0.0
            do is=1, n_species
               fac = z(is) * dens(is,ir) * dens_fac(is,it)
               sum_zn = sum_zn + z(is) / temp_para(is,ir) * fac
            enddo
            if(adiabatic_ele_model == 1) then
               fac = -ne_ade(ir) * exp(1.0/te_ade(ir) * phi_rot(it))
               sum_zn = sum_zn - 1.0/te_ade(ir) * fac
            endif
            
            do is=1,n_species
               fac = z(is) * dens(is,ir) * dens_fac(is,it)
               phi_rot_rderiv(it) = phi_rot_rderiv(it) &
                    + fac * (-dlnndr_fac(is,it))
            enddo
            if(adiabatic_ele_model == 1) then
               fac = -ne_ade(ir) * exp(1.0/te_ade(ir) * phi_rot(it))
               phi_rot_rderiv(it) = phi_rot_rderiv(it) &
                    + fac * (-dlnndre_ade(ir) &
                    + phi_rot(it)/te_ade(ir) * dlntdre_ade(ir))
            endif
            phi_rot_rderiv(it) = phi_rot_rderiv(it) / sum_zn
         enddo

         ! d phi*/d theta
         do it=1,n_theta
            phi_rot_deriv(it) = 0.0
            do id=-2,2
               if (id /= 0) then
                  jt = thcyc(it+id)
                  phi_rot_deriv(it) = phi_rot_deriv(it) &
                       + phi_rot(jt) * cderiv(id) / (12.0*d_theta)
               endif
            enddo
         enddo

         do is=1,n_species
            rotavg_e0(is)  = 0.0
            rotavg_e1(is)  = 0.0
            rotavg_e2(is)  = 0.0
            rotavg_e3(is)  = 0.0
            rotavg_e4(is)  = 0.0
            rotavg_e5a(is) = 0.0
            rotavg_e5b(is) = 0.0
            rotavg_e6(is)  = 0.0
            do it=1,n_theta
               fac = exp(omega_rot(ir)**2 * 0.5 / vth_para(is,ir)**2 &
                    * (bigR(it)**2 - bigR_th0**2) &
                    - z(is) / temp_para(is,ir) * phi_rot(it) &
                    - lam_rot_aniso(is,it))
               rotavg_e0(is) = rotavg_e0(is) + w_theta(it) * fac
               rotavg_e1(is) = rotavg_e1(is) + w_theta(it) * fac &
                    * z(is)/temp_para(is,ir) * phi_rot(it)
               rotavg_e2(is) = rotavg_e2(is) + w_theta(it) * fac &
                    * z(is)/temp_para(is,ir) * phi_rot_rderiv(it)
               rotavg_e3(is) = rotavg_e3(is) + w_theta(it) * fac &
                    * (bigR(it)**2 - bigR_th0**2)
               rotavg_e4(is) = rotavg_e4(is) + w_theta(it) * fac &
                    * 2.0*(bigR(it)*bigR_rderiv(it) &
                    - bigR_th0*bigR_th0_rderiv)
               rotavg_e5a(is) = rotavg_e5a(is) + w_theta(it) &
                    * jacobln_rderiv(it)
               rotavg_e5b(is) = rotavg_e5b(is) + w_theta(it) * fac &
                    * jacobln_rderiv(it)
               if(aniso_model(is) == 2) then
                  rotavg_e6(is) = rotavg_e6(is) + w_theta(it) * fac &
                       * (-lam_rot_rderiv_aniso(is,it))
               endif
            enddo
            rotavg_e3(is) = rotavg_e3(is) &
                 * (omega_rot(ir)/vth_para(is,ir)**2 &
                 * omega_rot_deriv(ir) - 0.5 * omega_rot(ir)**2 &
                 / vth_para(is,ir)**2 * (-dlntdr_para(is,ir)) )
            rotavg_e4(is) = rotavg_e4(is) * 0.5 * omega_rot(ir)**2 &
                 / vth_para(is,ir)**2
            rotavg_e1(is) = rotavg_e1(is) * (-dlntdr_para(is,ir))
            rotavg_e5a(is) = -rotavg_e5a(is) * rotavg_e0(is)
         enddo
         
         ! redefine nm wrt Teff and keeping <nm> fixed
         do is=1,n_species
            if(aniso_model(is) == 2) then
               dens_new = dens(is,ir)
               dens_fac_new(is,:) = exp(omega_rot(ir)**2 * 0.5/vth(is,ir)**2 &
                    * (bigR(:)**2 - bigR_th0**2) &
                    - z(is) / temp(is,ir) * phi_rot(:) &
                    - lam_rot_aniso(is,:))
               dens_avg_new(is)   = 0.0
               do it=1,n_theta
                  dens_avg_new(is) = dens_avg_new(is) &
                       + w_theta(it) * dens_fac_new(is,it)
               enddo
               dens_new   = dens_avg(is) / dens_avg_new(is)
               nu(is,ir)  = nu(is,ir) / dens(is,ir)
               dens(is,ir) = dens_new
               dens_fac(is,:) = dens_fac_new(is,:)
               nu(is,ir)  = nu(is,ir) * dens(is,ir)
            endif
         enddo
         
         deallocate(dens_avg)
         deallocate(dlnndr_fac)
         deallocate(dens_avg_new)
         deallocate(dens_fac_new)
         deallocate(lam_rot_rderiv_aniso)
         
      endif
      
    end subroutine ROT_solve_phi
    
    subroutine ROT_write(ir)
      use neo_globals
      implicit none
      integer, intent (in) :: ir
      integer :: is, it

      if(silent_flag == 0 .and. i_proc == 0 .and. rotation_model == 2) then
         open(io_rot,file=trim(path)//runfile,status='old',position='append')
         write (io_rot,'(e16.8)',advance='no') r(ir)
         write (io_rot,'(e16.8)',advance='no') phi_rot_avg
         do is=1, n_species
            write (io_rot,'(e16.8)',advance='no') 1.0/rotavg_e0(is)
            write (io_rot,'(e16.8)',advance='no') ( -rotavg_e2(is) &
                 + rotavg_e3(is) + rotavg_e4(is) + rotavg_e1(is) &
                 + rotavg_e5a(is) + rotavg_e5b(is) + rotavg_e6(is)) &
                 / rotavg_e0(is) 
         enddo
         do it=1, n_theta
            write (io_rot,'(e16.8)',advance='no') phi_rot(it)
         enddo
         do is=1, n_species
            do it=1, n_theta
                write (io_rot,'(e16.8)',advance='no') dens_fac(is,it)
             enddo
          enddo
         write(io_rot,*)
         close(io_rot)
      endif

      if(silent_flag == 0 .and. i_proc == 0 .and. rotation_model == 2) then
            open(unit=io_rot,file=trim(path)//'out.neo.diagnostic_rot',&
                 status='replace')
            write(io_rot,'(a,i3)') "# n_theta           = ",n_theta
            write(io_rot,'(a,i3)') "# n_species         = ",n_species
            write(io_rot,'(a)') "# Functions:"
            write(io_rot,'(a)') "#   theta(:)"
            write(io_rot,'(a)') "#   phi_rot (theta)"
            write(io_rot,'(a)') "#   dphi_rot/dr (theta)"
            write(io_rot,'(a)') "#   do is=1,n_species: total e0,e1,e2,e3,e4,e5,eeta terms"
            write(io_rot,'(1pe16.8)') theta(:)
            write(io_rot,'(1pe16.8)') phi_rot(:)
            write(io_rot,'(1pe16.8)') phi_rot_rderiv(:)
            do is=1,n_species
               write(io_rot,'(1pe16.8)') rotavg_e0(is)
               write(io_rot,'(1pe16.8)') rotavg_e1(is)
               write(io_rot,'(1pe16.8)') rotavg_e2(is)
               write(io_rot,'(1pe16.8)') rotavg_e3(is)
               write(io_rot,'(1pe16.8)') rotavg_e4(is)
               write(io_rot,'(1pe16.8)') rotavg_e5a(is)+rotavg_e5b(is)
               write(io_rot,'(1pe16.8)') rotavg_e6(is)
            enddo
            close(io_rot)
         end if


    end subroutine ROT_write

end module neo_rotation
