program pneo

  use mpi
  use pneo_globals
  use neo_interface
  use EXPRO_interface

  implicit none

  integer :: p,k, is
  integer :: i1,i2,i3,i4,i5,i6,i7,i8,i9
  integer :: n1,n2,n3,n4,n5,n6,n7,n8,n9
  real, dimension(:), allocatable :: rmin_over_rmaj
  real, dimension(:), allocatable :: q
  real, dimension(:), allocatable :: nu_ee
  real, dimension(:), allocatable :: ni_over_ne
  real, dimension(:), allocatable :: ti_over_te
  real, dimension(:), allocatable :: delta
  real, dimension(:), allocatable :: s_delta
  real, dimension(:), allocatable :: kappa
  real, dimension(:), allocatable :: s_kappa
  real, dimension(:), allocatable :: jfac
  
  allocate(rmin_over_rmaj(9))
  allocate(q(9))
  allocate(nu_ee(9))
  allocate(ni_over_ne(9))
  allocate(ti_over_te(9))
  allocate(delta(9))
  allocate(s_delta(9))
  allocate(kappa(9))
  allocate(s_kappa(9))

  open(unit=1,file='input.pneo',status='old')
  read(1,*) n1
  read(1,*) rmin_over_rmaj(1:n1)
  read(1,*) n2
  read(1,*) q(1:n2)
  read(1,*) n3
  read(1,*) nu_ee(1:n3)
  read(1,*) n4
  read(1,*) ni_over_ne(1:n4)
  read(1,*) n5
  read(1,*) ti_over_te(1:n5)
  read(1,*) n6
  read(1,*) delta(1:n6)
  read(1,*) n7
  read(1,*) s_delta(1:n7)
  read(1,*) n8
  read(1,*) kappa(1:n8)
  read(1,*) n9
  read(1,*) s_kappa(1:n9)
  close(1)

  !---------------------------------------------------------------------
  ! Initialize MPI_COMM_WORLD communicator.
  !
  call MPI_INIT(i_err)
  call MPI_COMM_RANK(MPI_COMM_WORLD,i_proc,i_err)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,n_proc,i_err)
  !---------------------------------------------------------------------

  ! Path is cwd:
  path= './'

  call neo_init_serial(path)

  ! pointers

  ntot = n1*n2*n3*n4*n5*n6*n7*n8*n9
  allocate(ic1(ntot))
  allocate(ic2(ntot))
  allocate(ic3(ntot))
  allocate(ic4(ntot))
  allocate(ic5(ntot))
  allocate(ic6(ntot))
  allocate(ic7(ntot))
  allocate(ic8(ntot))
  allocate(ic9(ntot))
  
  allocate(indata_loc(7,ntot))
  allocate(indata(7,ntot))

  allocate(outdata_j_loc(6,ntot))
  allocate(outdata_j(6,ntot))
  allocate(outdata_u_loc(18,ntot))
  allocate(outdata_u(18,ntot))
  allocate(outdata_g_loc(18,ntot))
  allocate(outdata_g(18,ntot))
  allocate(outdata_q_loc(18,ntot))
  allocate(outdata_q(18,ntot))
  
  p = 0
  do i1=1,n1
     do i2=1,n2
        do i3=1,n3
           do i4=1,n4
              do i5=1,n5
                 do i6=1,n6
                    do i7=1,n7
                       do i8=1,n8
                          do i9=1,n9

                             p = p+1
                             ic1(p) = i1
                             ic2(p) = i2
                             ic3(p) = i3
                             ic4(p) = i4
                             ic5(p) = i5
                             ic6(p) = i6
                             ic7(p) = i7
                             ic8(p) = i8
                             ic9(p) = i9
                          enddo
                       enddo
                    enddo
                 enddo
              enddo
           enddo
        enddo
     enddo
  enddo

  ! Fixed NEO subroutine inputs
  neo_silent_flag_in = 1
  neo_n_energy_in = 6
  neo_n_xi_in     = 17
  neo_n_theta_in  = 17
  neo_collision_model_in = 4   ! Full FP collisions
  neo_rho_star_in = 0.001      ! arbitrary

  ! geometry
  neo_equilibrium_model_in = 2 ! Miller equilibrium
  neo_rmaj_over_a_in = 1.0     ! anorm = rmaj

  ! for shaping test:
  !neo_n_xi_in     = 27
  !neo_n_theta_in  = 57
  !neo_shift_in = -0.26850
  !neo_zeta_in  = -4.2304e-2
  !neo_s_zeta_in = -2.2832e-1
  !neo_zmag_over_a_in = -0.0293479
  !neo_s_zmag_in = -8.8566e-1

  ! 3 species: ele + D + C
  neo_n_species_in = 3
  ! electrons
  neo_z_in(1)     = -1
  neo_mass_in(1)  = 0.0002724486
  neo_temp_in(1)  = 1.0        ! Tnorm = Te
  neo_dens_in(1)  = 1.0        ! nnorm = ne
  ! D
  neo_z_in(2)     = 1
  neo_mass_in(2)  = 1.0
  ! C
  neo_z_in(3)     = 6
  neo_mass_in(3)  = 6.0

  allocate(jfac(neo_n_species_in))
  
  ! For testing, use THEORY sim_model=0;
  ! For nn, use sim_model=4;
  ! else use NEO sim_model=1
  !neo_sim_model_in = 4
  !!!!!!

  if (i_proc == 0) print '(a,i5)','NTOT = ',ntot

  do p=1+i_proc,ntot,n_proc

     if (i_proc == 0) print '(i5,a,i5)',p,' - ',p+n_proc-1

     i1 = ic1(p) 
     i2 = ic2(p)
     i3 = ic3(p) 
     i4 = ic4(p)
     i5 = ic5(p) 
     i6 = ic6(p)
     i7 = ic7(p) 
     i8 = ic8(p)
     i9 = ic9(p)

     neo_rmin_over_a_in = rmin_over_rmaj(i1)
     neo_q_in           = abs(q(i2))
     neo_nu_1_in        = nu_ee(i3)
     neo_dens_in(2)     = ni_over_ne(i4)
     neo_dens_in(3)     = (1.0-neo_z_in(2)*neo_dens_in(2))/(1.0*neo_z_in(3))
     neo_temp_in(2)     = ti_over_te(i5)
     neo_temp_in(3)     = ti_over_te(i5)
     neo_delta_in       = delta(i6)
     neo_s_delta_in     = s_delta(i7)
     neo_kappa_in       = kappa(i8)  
     neo_s_kappa_in     = s_kappa(i9)

     do is=1,neo_n_species_in
        jfac(is) = neo_dens_in(is) * abs(neo_z_in(is))
     enddo
     
     ! Cne
     neo_dlnndr_in(:) = 0.0; neo_dlntdr_in(:) = 0.0; neo_dlnndr_in(1) = 1.0
     call neo_run()
     outdata_j_loc(1,p) = neo_jpar_dke_out/jfac(1)
     do is=1,neo_n_species_in
        outdata_u_loc(1+6*(is-1),p)  = neo_vpol_dke_out(is)
        outdata_g_loc(1+6*(is-1),p)  = neo_pflux_dke_out(is)
        outdata_q_loc(1+6*(is-1),p)  = neo_efluxtot_dke_out(is)
     enddo
     
     ! CTe
     neo_dlnndr_in(:) = 0.0;  neo_dlntdr_in(:) = 0.0; neo_dlntdr_in(1) = 1.0 
     call neo_run()
     outdata_j_loc(2,p) = neo_jpar_dke_out/jfac(1)
     do is=1,neo_n_species_in
        outdata_u_loc(2+6*(is-1),p)  = neo_vpol_dke_out(is)
        outdata_g_loc(2+6*(is-1),p)  = neo_pflux_dke_out(is)
        outdata_q_loc(2+6*(is-1),p)  = neo_efluxtot_dke_out(is)
     enddo
     
     ! Cni1
     neo_dlnndr_in(:) = 0.0;  neo_dlntdr_in(:) = 0.0; neo_dlnndr_in(2) = 1.0 
     call neo_run()
     outdata_j_loc(3,p) = neo_jpar_dke_out/jfac(2)
     do is=1,neo_n_species_in
        outdata_u_loc(3+6*(is-1),p)  = neo_vpol_dke_out(is)
        outdata_g_loc(3+6*(is-1),p)  = neo_pflux_dke_out(is)
        outdata_q_loc(3+6*(is-1),p)  = neo_efluxtot_dke_out(is)
     enddo
     
     ! CTi1
     neo_dlnndr_in(:) = 0.0;  neo_dlntdr_in(:) = 0.0; neo_dlntdr_in(2) = 1.0  
     call neo_run()
     outdata_j_loc(4,p) = neo_jpar_dke_out/jfac(2)
     do is=1,neo_n_species_in
        outdata_u_loc(4+6*(is-1),p)  = neo_vpol_dke_out(is)
        outdata_g_loc(4+6*(is-1),p)  = neo_pflux_dke_out(is)
        outdata_q_loc(4+6*(is-1),p)  = neo_efluxtot_dke_out(is)
     enddo
     
     ! Cni2
     neo_dlnndr_in(:) = 0.0;  neo_dlntdr_in(:) = 0.0; neo_dlnndr_in(3) = 1.0 
     call neo_run()
     outdata_j_loc(5,p) = neo_jpar_dke_out/jfac(3)
     do is=1,neo_n_species_in
        outdata_u_loc(5+6*(is-1),p)  = neo_vpol_dke_out(is)
        outdata_g_loc(5+6*(is-1),p)  = neo_pflux_dke_out(is)
        outdata_q_loc(5+6*(is-1),p)  = neo_efluxtot_dke_out(is)
     enddo
     
     ! CTi2 
     neo_dlnndr_in(:) = 0.0;  neo_dlntdr_in(:) = 0.0; neo_dlntdr_in(3) = 1.0
     call neo_run()
     outdata_j_loc(6,p) = neo_jpar_dke_out/jfac(3)
     do is=1,neo_n_species_in
        outdata_u_loc(6+6*(is-1),p)  = neo_vpol_dke_out(is)
        outdata_g_loc(6+6*(is-1),p)  = neo_pflux_dke_out(is)
        outdata_q_loc(6+6*(is-1),p)  = neo_efluxtot_dke_out(is)
     enddo

     ! <jpar B>/Bunit j_s ~ rho (I/psip) sum_s |z_a| n_a C_a 1/L_a
     outdata_j_loc(:,p) = outdata_j_loc(:,p) &
           / (neo_rho_star_in * neo_geoparams_out(1))
          
     ! K_a/n_a ~ Vpol / Bpol
     ! <B^2/Bunit^2> K_a Bunit/(n_a c_s) ~ rho (I/psip) C 1/L
     outdata_u_loc(:,p) = outdata_u_loc(:,p) &
          * neo_geoparams_out(3) / neo_geoparams_out(4) &
          / (neo_rho_star_in * neo_geoparams_out(1))
  
     ! Gamma_a/(n_e c_s) ~ nu_ee rho_s^2 (I/psip)^2 / <B^2/Bunit^2> C 1/L
     outdata_g_loc(:,p) = outdata_g_loc(:,p) &
          * neo_geoparams_out(3) &
          / (neo_rho_star_in * neo_geoparams_out(1))**2 / neo_nu_1_in

     ! Q_a/(n_e c_s T_e) ~ nu_ee rho_a^2 (I/psip)^2  <B^2/Bunit^2> C 1/L
     outdata_q_loc(:,p) = outdata_q_loc(:,p) &
          * neo_geoparams_out(3) &
          / (neo_rho_star_in * neo_geoparams_out(1))**2 / neo_nu_1_in

     do is=1,neo_n_species_in
        outdata_u_loc(2+6*(is-1),p) = outdata_u_loc(2+6*(is-1),p) &
             + 1.5*outdata_u_loc(1+6*(is-1),p)
        outdata_u_loc(4+6*(is-1),p) = outdata_u_loc(4+6*(is-1),p) &
             + 1.5*outdata_u_loc(3+6*(is-1),p)
        outdata_u_loc(6+6*(is-1),p) = outdata_u_loc(6+6*(is-1),p) &
             + 1.5*outdata_u_loc(5+6*(is-1),p)
        outdata_g_loc(2+6*(is-1),p) = outdata_g_loc(2+6*(is-1),p) &
             + 1.5*outdata_g_loc(1+6*(is-1),p)
        outdata_g_loc(4+6*(is-1),p) = outdata_g_loc(4+6*(is-1),p) &
             + 1.5*outdata_g_loc(3+6*(is-1),p)
        outdata_g_loc(6+6*(is-1),p) = outdata_g_loc(6+6*(is-1),p) &
             + 1.5*outdata_g_loc(5+6*(is-1),p)
        outdata_q_loc(2+6*(is-1),p) = outdata_q_loc(2+6*(is-1),p) &
             + 1.5*outdata_q_loc(1+6*(is-1),p)
        outdata_q_loc(4+6*(is-1),p) = outdata_q_loc(4+6*(is-1),p) &
             + 1.5*outdata_q_loc(3+6*(is-1),p)
        outdata_q_loc(6+6*(is-1),p) = outdata_q_loc(6+6*(is-1),p) &
             + 1.5*outdata_q_loc(5+6*(is-1),p)
     enddo
     
     ! 7 inputs: eps,ft,q,log10(nuee),ni,Ti, <B^2><1/B^2>-1
     indata_loc(1,p) = neo_rmin_over_a_in
     indata_loc(2,p) = neo_geoparams_out(2)
     indata_loc(3,p) = neo_q_in
     indata_loc(4,p) = log10(neo_nu_1_in)
     indata_loc(5,p) = neo_dens_in(2)
     indata_loc(6,p) = neo_temp_in(2)
     indata_loc(7,p) = neo_geoparams_out(5)
     
  enddo

  ! Collect all data 
  call MPI_ALLREDUCE(indata_loc,indata,size(indata), &
       MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,i_err)
  call MPI_ALLREDUCE(outdata_j_loc,outdata_j,size(outdata_j), &
       MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,i_err)
  call MPI_ALLREDUCE(outdata_u_loc,outdata_u,size(outdata_u), &
       MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,i_err)
  call MPI_ALLREDUCE(outdata_g_loc,outdata_g,size(outdata_g), &
       MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,i_err)
  call MPI_ALLREDUCE(outdata_q_loc,outdata_q,size(outdata_q), &
       MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,i_err)
  
  if (i_proc == 0) then

     open(unit=1,file='out.pneo.indata',status='replace')
     write(1,10) indata(:,:)
     close(1)

     open(unit=1,file='out.pneo.c_j',status='replace')
     write(1,10) outdata_j(:,:)
     close(1)

     open(unit=1,file='out.pneo.c_u',status='replace')
     write(1,10) outdata_u(:,:)
     close(1)

     open(unit=1,file='out.pneo.c_g',status='replace')
     write(1,10) outdata_g(:,:)
     close(1)

     open(unit=1,file='out.pneo.c_q',status='replace')
     write(1,10) outdata_q(:,:)
     close(1)
     
  endif

  call MPI_finalize(i_err)

10 format(1pe17.10)

end program pneo
