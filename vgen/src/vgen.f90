! NOTES
! - negative ion density check (part of sanity checker)


!------------------------------------------------------------------------
! vgen.f90
!
! PURPOSE: 
!  Driver for the vgen (velocity-generation) capability of NEO.  This 
!  will write a new input.gacode with NEO-computed electric field 
!  and/or velocities. 
!------------------------------------------------------------------------

program vgen

  use mpi
  use vgen_globals
  use neo_interface
  use expro

  implicit none

  integer :: i
  integer :: j
  integer :: ix
  integer :: rotation_model 
  real :: grad_p
  real :: ya
  real :: yb
  real :: vtor_diff
  real :: er0
  real :: omega
  real :: omega_deriv
  integer :: simntheta
  integer :: iteration_flag
  real :: cpu_tot_in, cpu_tot_out
  character(len=14) :: er_tag
  character(len=17) :: vel_tag
  character(len=7)  :: ix_tag
  character(len=15) :: j_tag
  
  real, dimension(:), allocatable :: er_exp

  !---------------------------------------------------------------------
  ! Initialize MPI_COMM_WORLD communicator.
  !
  call MPI_INIT(i_err)
  call MPI_COMM_RANK(MPI_COMM_WORLD,i_proc,i_err)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,n_proc,i_err)
  !---------------------------------------------------------------------

  cpu_tot_in = MPI_Wtime()

  ! Path is cwd:
  path= './'

  ! Obscure definition of number tags 
  do ix=1,100
     write (tag(ix),fmt) ix-1
  enddo

  ! Read vgen control parameters

  open(unit=1,file='vgen.dat',status='old')
  read(1,*) er_method
  read(1,*) vel_method
  read(1,*) erspecies_indx
  read(1,*) epar_flag
  read(1,*) nth_min
  read(1,*) nth_max
  read(1,*) nn_flag
  close(1)

  select case(er_method)
  case(1)
     if(i_proc == 0) then
        print '(a)','INFO: (VGEN) Computing omega0 (Er)  from force balance'
        er_tag=' -er 1 (FB)'
        write(ix_tag,'(i0)') erspecies_indx
        ix_tag=' -ix ' // trim(ix_tag)
     endif
  case(2)
     if(i_proc == 0) then
        print '(a)', 'INFO: (VGEN) Computing omega0 (Er) from NEO (weak rotation limit)'
        er_tag=' -er 2 (NEO)'
        write(ix_tag,'(i0)') erspecies_indx
        ix_tag=' -ix ' // trim(ix_tag)
     endif
  case(4)
     if(i_proc == 0) then
        print '(a)','INFO: (VGEN) Returning given omega0 (Er)'
        er_tag=' -er 4 (given)'
        ix_tag=''
     endif
  case default
     if(i_proc == 0) then
        print '(a)','ERROR: Invalid er_method'
     endif
     call MPI_finalize(i_err)
     stop
  end select

  select case(vel_method)
  case(1)
     if (i_proc == 0) then
        print '(a)','INFO: (VGEN) Computing velocities from NEO (weak rotation limit)'
        vel_tag=' -vel 1 (weak)'
     endif
  case(2)
     if (i_proc == 0) then
        print '(a)','INFO: (VGEN) Computing velocities from NEO (strong rotation limit)'
        vel_tag=' -vel 2 (strong)'
     endif
  case default
     if (i_proc == 0) then
        print '(a)','ERROR: Invalid vel_method'
     endif
     call MPI_finalize(i_err)
     stop
  end select

  !---------------------------------------------------------------------
  ! Initialize vgen parameters
  !
  call vgen_init
  allocate(er_exp(EXPRO_n_exp))
  if (er_method /= 4) then
     er_exp(:)    = 0.0
     EXPRO_w0(:)  = 0.0
     EXPRO_w0p(:) = 0.0
  endif
  !---------------------------------------------------------------------
  
  select case(neo_sim_model_in)
  case(0)
     if (i_proc == 0) then
        print '(a)','INFO: (VGEN) Bootstrap current from Sauter'
        j_tag=' -jbs (Sauter)'
     endif
  case(1)
     if (i_proc == 0) then
        print '(a)','INFO: (VGEN) Bootstrap current from NEO'
        j_tag=' -jbs (NEO)'
     endif
  case(2)
     if (i_proc == 0) then
        print '(a)','INFO: (VGEN) Bootstrap current from NEO'
        j_tag=' -jbs (NEO)'
     endif
  case(3)
     if (i_proc == 0) then
        print '(a)','INFO: (VGEN) Bootstrap current from Sauter'
        j_tag=' -jbs (Sauter)'
     endif
  case(4)
     if (i_proc == 0) then
        print '(a)','INFO: (VGEN) Bootstrap current from NEO NN'
        j_tag=' -jbs (NEO NN)'
     endif
  case(5)
     if (i_proc == 0) then
        print '(a)','INFO: (VGEN) Bootstrap current from modified Sauter'
        j_tag=' -jbs (Sauter mod)'
     endif
  case default
     if (i_proc == 0) then
        print '(a)','ERROR: Invalid neo_sim_model'
     endif
     call MPI_finalize(i_err)
     stop
  end select
  
  select case(epar_flag)
  case(0)
     if (i_proc == 0) then
        print '(a)','INFO: (VGEN) Do not include conductivity calculation'
     endif
  case(1)
     if (i_proc == 0) then
        print '(a)','INFO: (VGEN) Include conductivity calculation'
     endif
  case default
     if (i_proc == 0) then
        print '(a)','ERROR: Invalid epar_flag'
     endif
     call MPI_finalize(i_err)
     stop
  end select
  
  if (i_proc == 0) then
     print '(a,i2,a,i2)','INFO: (VGEN) Using NEO Theta Resolution: ',nth_min,'-',nth_max
     print '(a,i2)','INFO: (VGEN) MPI tasks: ',n_proc
  endif

  !---------------------------------------------------------------------
  ! Distribution scheme.
  !
  n_loc = EXPRO_n_exp/n_proc+1
  allocate(i_glob(n_loc))

  i_loc = 0
  do i=2+i_proc,EXPRO_n_exp-1,n_proc
     i_loc = i_loc+1
     i_glob(i_loc) = i
  enddo
  n_loc = i_loc
  !---------------------------------------------------------------------

  !======================================================================
  ! Four alternatives for Er calculation:
  !
  ! 1. Compute Er from force balance using measured vpol and vtor 
  ! 2. Compute Er by matching vtor_measured with vtor_neo at theta=0 assuming
  !    weak rotation
  ! 3. Compute Er by matching vtor_measured with vtor_neo at theta=0 assuming
  !    strong rotation
  ! 4. Return the given Er

  select case (er_method) 

  case (1,4)

     if (er_method == 1) then

        ! Compute Er from force balance using measured vpol and vtor

        ! Er calculation first

        do i_loc=1,n_loc
           i = i_glob(i_loc)
           if (erspecies_indx == 1) then
              grad_p = -(EXPRO_dlnnidr_new(i) + EXPRO_dlntidr(1,i))
           else
              grad_p = -(EXPRO_dlnnidr(erspecies_indx,i) &
                   + EXPRO_dlntidr(erspecies_indx,i))
           endif
           er_exp(i) = (grad_p * EXPRO_grad_r0(i) &
                * EXPRO_ti(erspecies_indx,i)*temp_norm_fac &
                / (EXPRO_z(erspecies_indx) * charge_norm_fac) & 
                + EXPRO_vtor(erspecies_indx,i) * EXPRO_bp0(i) &
                - EXPRO_vpol(erspecies_indx,i) * EXPRO_bt0(i)) &
                / 1000
        enddo

        call vgen_reduce(er_exp(2:EXPRO_n_exp-1),EXPRO_n_exp-2)

        ! Compute omega and omega_deriv from newly-generated Er:

        do i=2,EXPRO_n_exp-1
           EXPRO_w0(i) = 2.9979e10*EXPRO_q(i)*(er_exp(i)/30.0)/ &
                ((1e4*EXPRO_bunit(i))*(1e2*EXPRO_rmin(i))*EXPRO_grad_r0(i))
        enddo

        ! Compute w0p
        call bound_deriv(EXPRO_w0p(2:EXPRO_n_exp-1),EXPRO_w0(2:EXPRO_n_exp-1),&
             EXPRO_rmin,EXPRO_n_exp-2)

     endif

     if (er_method == 4) then

        do i=1,EXPRO_n_exp
           er_exp(i) = 1.0/2.9979e10/EXPRO_q(i)*EXPRO_w0(i)*30.0* &
                ((1e4*EXPRO_bunit(i))*(1e2*EXPRO_rmin(i))*EXPRO_grad_r0(i))
        enddo

     endif

     ! Flow calculation based on existing EXPRO_w0 and EXPRO_w0p

     do i_loc=1,n_loc

        i = i_glob(i_loc)
        if (vel_method == 1) then
           rotation_model = 1   ! weak rotation
        else
           rotation_model = 2   ! strong rotation
        endif
        er0 = er_exp(i)
        omega = EXPRO_w0(i) 
        omega_deriv = EXPRO_w0p(i) 

        iteration_flag = 1
        call vgen_compute_neo(i,vtor_diff,rotation_model,er0,omega,omega_deriv, simntheta,iteration_flag)

        print 10,EXPRO_rho(i),&
             er_exp(i),EXPRO_vpol(1,i)/1e3,simntheta,i_proc

     enddo

     ! Reduce vpol,vtor
     do j=1,n_ions
        call vgen_reduce(EXPRO_vpol(j,2:EXPRO_n_exp-1),EXPRO_n_exp-2)
        call vgen_reduce(EXPRO_vtor(j,2:EXPRO_n_exp-1),EXPRO_n_exp-2)
     enddo

  case (2)

     ! Compute Er using NEO (weak rotation limit) 
     ! by matching vtor_measured with vtor_neo at theta=0 

     do i_loc=1,n_loc

        i = i_glob(i_loc)

        rotation_model = 1
        er0 = 0.0
        omega = 0.0
        omega_deriv = 0.0

        if (vel_method == 2 .and. epar_flag == 1) then
           iteration_flag=2
        else
           iteration_flag=1
        endif
        
        call vgen_compute_neo(i,vtor_diff, rotation_model, er0, omega, &
             omega_deriv, simntheta,iteration_flag)
        iteration_flag=1

        ! omega = (vtor_measured - vtor_neo_ater0) / R

        er_exp(i) = (vtor_diff/(vth_norm * EXPRO_rmin(EXPRO_n_exp)))  &
             / (neo_rmaj_over_a_in + neo_rmin_over_a_in) &
             * neo_rmin_over_a_in / (abs(neo_q_in) * EXPRO_signq) &
             / (abs(neo_rho_star_in) * EXPRO_signb) &
             * EXPRO_grad_r0(i) &
             * (temp_norm*temp_norm_fac / charge_norm_fac &
             / EXPRO_rmin(EXPRO_n_exp)) / 1000

        ! Part of vtor due to Er
        do j=1,n_ions
           EXPRO_vtor(j,i) = EXPRO_vtor(j,i) + vtor_diff
        enddo

        print 10,EXPRO_rho(i),&
             er_exp(i),EXPRO_vpol(1,i)/1e3,simntheta, i_proc

     enddo

     ! Reduce er,vpol,vtor
     call vgen_reduce(er_exp(2:EXPRO_n_exp-1),EXPRO_n_exp-2)

     if (vel_method == 1) then

        do j=1,n_ions
           call vgen_reduce(EXPRO_vpol(j,2:EXPRO_n_exp-1),EXPRO_n_exp-2)
           call vgen_reduce(EXPRO_vtor(j,2:EXPRO_n_exp-1),EXPRO_n_exp-2)
        enddo

     else

        EXPRO_vpol(:,:) = 0.0
        EXPRO_vtor(:,:) = 0.0

        if (i_proc == 0) print '(a)', 'INFO: (VGEN) Recomputing flows using NEO in the strong rotation limit.'

        ! Re-compute the flows using strong rotation
        ! omega and omega_deriv 
        do i=2,EXPRO_n_exp-1
           EXPRO_w0(i) = 2.9979e10*EXPRO_q(i)*(er_exp(i)/30.0)/ &
                ((1e4*EXPRO_bunit(i))*(1e2*EXPRO_rmin(i))*EXPRO_grad_r0(i))
        enddo
        ! Compute w0p
        call bound_deriv(EXPRO_w0p(2:EXPRO_n_exp-1),EXPRO_w0(2:EXPRO_n_exp-1),&
             EXPRO_rmin,EXPRO_n_exp-2)

        do i_loc=1,n_loc

           i = i_glob(i_loc)

           rotation_model = 2  
           er0 = er_exp(i)
           omega = EXPRO_w0(i) 
           omega_deriv = EXPRO_w0p(i) 
           call vgen_compute_neo(i,vtor_diff, rotation_model, er0, omega, &
                omega_deriv, simntheta,iteration_flag)

           print 10,EXPRO_rho(i),&
                er_exp(i),EXPRO_vpol(1,i)/1e3,simntheta,i_proc

        enddo

        ! Reduce vpol,vtor
        do j=1,n_ions
           call vgen_reduce(EXPRO_vpol(j,2:EXPRO_n_exp-1),EXPRO_n_exp-2)
           call vgen_reduce(EXPRO_vtor(j,2:EXPRO_n_exp-1),EXPRO_n_exp-2)
        enddo

     endif

  end select
  !======================================================================

  ! Additional reductions
  call vgen_reduce(pflux_sum(2:EXPRO_n_exp-1),EXPRO_n_exp-2)
  call vgen_reduce(jbs_neo(2:EXPRO_n_exp-1),EXPRO_n_exp-2)
  call vgen_reduce(jsigma_neo(2:EXPRO_n_exp-1),EXPRO_n_exp-2)
  call vgen_reduce(jtor_neo(2:EXPRO_n_exp-1),EXPRO_n_exp-2)
  call vgen_reduce(jbs_sauter(2:EXPRO_n_exp-1),EXPRO_n_exp-2)
  call vgen_reduce(jsigma_sauter(2:EXPRO_n_exp-1),EXPRO_n_exp-2)
  call vgen_reduce(jtor_sauter(2:EXPRO_n_exp-1),EXPRO_n_exp-2)
  call vgen_reduce(jbs_sauter_mod(2:EXPRO_n_exp-1),EXPRO_n_exp-2)
  call vgen_reduce(jsigma_sauter_mod(2:EXPRO_n_exp-1),EXPRO_n_exp-2)
  call vgen_reduce(jtor_sauter_mod(2:EXPRO_n_exp-1),EXPRO_n_exp-2)
  
  !------------------------------------------------------------------------
  ! Extrapolation for r=0 and r=n_exp boundary points
  !
  if (er_method /= 4) then

     call bound_extrap(ya,yb,er_exp,EXPRO_rmin,EXPRO_n_exp)
     er_exp(1) = ya
     er_exp(EXPRO_n_exp) = yb

     do i=2,EXPRO_n_exp-1
        EXPRO_w0(i) = 2.9979e10*EXPRO_q(i)*(er_exp(i)/30.0)/ &
             ((1e4*EXPRO_bunit(i))*(1e2*EXPRO_rmin(i))*EXPRO_grad_r0(i))
     enddo
     call bound_extrap(ya,yb,EXPRO_w0,EXPRO_rmin,EXPRO_n_exp)
     EXPRO_w0(1) = ya
     EXPRO_w0(EXPRO_n_exp) = yb

     call bound_deriv(EXPRO_w0p(2:EXPRO_n_exp-1),EXPRO_w0(2:EXPRO_n_exp-1),&
          EXPRO_rmin,EXPRO_n_exp-2)
     call bound_extrap(ya,yb,EXPRO_w0p,EXPRO_rmin,EXPRO_n_exp)
     EXPRO_w0p(1) = ya
     EXPRO_w0p(EXPRO_n_exp) = yb

  endif
  !
  do j=1,n_ions
     call bound_extrap(ya,yb,EXPRO_vpol(j,:),EXPRO_rmin,EXPRO_n_exp)
     EXPRO_vpol(j,1) = ya
     EXPRO_vpol(j,EXPRO_n_exp) = yb
     call bound_extrap(ya,yb,EXPRO_vtor(j,:),EXPRO_rmin,EXPRO_n_exp)
     EXPRO_vtor(j,1) = ya
     EXPRO_vtor(j,EXPRO_n_exp) = yb
  enddo
  !
  call bound_extrap(ya,yb,jbs_neo,EXPRO_rmin,EXPRO_n_exp)
  jbs_neo(1)           = ya
  jbs_neo(EXPRO_n_exp) = yb
  call bound_extrap(ya,yb,jsigma_neo,EXPRO_rmin,EXPRO_n_exp)
  jsigma_neo(1)           = ya
  jsigma_neo(EXPRO_n_exp) = yb
  call bound_extrap(ya,yb,jtor_neo,EXPRO_rmin,EXPRO_n_exp)
  jtor_neo(1)           = ya
  jtor_neo(EXPRO_n_exp) = yb
  call bound_extrap(ya,yb,jbs_sauter,EXPRO_rmin,EXPRO_n_exp)
  jbs_sauter(1)           = ya
  jbs_sauter(EXPRO_n_exp) = yb
  call bound_extrap(ya,yb,jsigma_sauter,EXPRO_rmin,EXPRO_n_exp)
  jsigma_sauter(1)           = ya
  jsigma_sauter(EXPRO_n_exp) = yb
  call bound_extrap(ya,yb,jtor_sauter,EXPRO_rmin,EXPRO_n_exp)
  jtor_sauter(1)           = ya
  jtor_sauter(EXPRO_n_exp) = yb
  call bound_extrap(ya,yb,jbs_sauter_mod,EXPRO_rmin,EXPRO_n_exp)
  jbs_sauter_mod(1)           = ya
  jbs_sauter_mod(EXPRO_n_exp) = yb
  call bound_extrap(ya,yb,jsigma_sauter_mod,EXPRO_rmin,EXPRO_n_exp)
  jsigma_sauter_mod(1)           = ya
  jsigma_sauter_mod(EXPRO_n_exp) = yb
  call bound_extrap(ya,yb,jtor_sauter_mod,EXPRO_rmin,EXPRO_n_exp)
  jtor_sauter_mod(1)           = ya
  jtor_sauter_mod(EXPRO_n_exp) = yb
  call bound_extrap(ya,yb,pflux_sum,EXPRO_rmin,EXPRO_n_exp)
  pflux_sum(1)           = ya
  pflux_sum(EXPRO_n_exp) = yb
  !------------------------------------------------------------------------
  
  if(neo_sim_model_in == 0 .or. neo_sim_model_in == 3) then
     EXPRO_jbs(:)      = jbs_sauter(:)
     EXPRO_jbstor(:)   = jtor_sauter(:)
     EXPRO_sigmapar(:) = jsigma_sauter(:)
  else if(neo_sim_model_in == 5) then
     EXPRO_jbs(:)      = jbs_sauter_mod(:)
     EXPRO_jbstor(:)   = jtor_sauter_mod(:)
     EXPRO_sigmapar(:) = jsigma_sauter_mod(:)
  else
     EXPRO_jbs(:)      = jbs_neo(:)
     EXPRO_jbstor(:)   = jtor_neo(:)
     EXPRO_sigmapar(:) = jsigma_neo(:)
  endif
  
  !------------------------------------------------------------------------
  
  ! Write output on processor 0

  if (i_proc == 0) then

     open(unit=1,file='out.vgen.vel',status='replace')
     do i=1,EXPRO_n_exp
        write(1,'(e16.8)',advance='no') EXPRO_rho(i)
        write(1,'(e16.8)',advance='no') er_exp(i)
        write(1,'(e16.8)',advance='no') EXPRO_w0(i)
        write(1,'(e16.8)',advance='no') EXPRO_w0p(i)
        do j=1,n_ions
           write(1,'(e16.8)',advance='no') EXPRO_vpol(j,i)
           write(1,'(e16.8)',advance='no') EXPRO_vtor(j,i)
           write(1,'(e16.8)',advance='no') EXPRO_vpol(j,i) / EXPRO_bp0(i)
           write(1,'(e16.8)',advance='no') (EXPRO_vtor(j,i) &
                - EXPRO_vpol(j,i)*EXPRO_bt0(i)/ EXPRO_bp0(i)) &
                / (EXPRO_rmaj(i)+EXPRO_rmin(i))
        enddo
        write(1,*)
     enddo
     close(1)
     
     open(unit=1,file='out.vgen.ercomp',status='replace')
     do i=1,EXPRO_n_exp
        write(1,'(e16.8)',advance='no') EXPRO_rho(i)
        write(1,'(e16.8)',advance='no') EXPRO_w0(i)
        do j=1,n_ions
           if(j == 1) then
              grad_p = -(EXPRO_dlnnidr_new(i) + EXPRO_dlntidr(j,i))
           else
              grad_p = -(EXPRO_dlnnidr(j,i) + EXPRO_dlntidr(j,i))
           endif
           write(1,'(e16.8)',advance='no') grad_p * EXPRO_grad_r0(i) &
                * EXPRO_ti(j,i)*temp_norm_fac &
                / (EXPRO_z(j) * charge_norm_fac) / (EXPRO_rmaj(i) + EXPRO_rmin(i)) &
                / EXPRO_bp0(i)
           write(1,'(e16.8)',advance='no') EXPRO_vtor(j,i) / (EXPRO_rmaj(i) + EXPRO_rmin(i))
           write(1,'(e16.8)',advance='no') -EXPRO_vpol(j,i) * EXPRO_bt0(i)/EXPRO_bp0(i) &
                / (EXPRO_rmaj(i) + EXPRO_rmin(i))
        enddo
        write(1,*)
     enddo
     close(1)
     
     !----------------------------------------------------------------------
     ! Generate new input.gacode.* files

     ! 1. input.gacode
     !call expro_write_original('input.gacode','input.gacode.new',' ')

     ! NEW APPROACH
     expro_head_vgen = '#      vgen : ' //  trim(er_tag) // trim(ix_tag) // trim(vel_tag) // trim(j_tag)
     call expro_write('input.gacode')
     print '(a)', 'INFO: (VGEN) Created vgen/input.gacode [new data]'

     ! 3. input.profiles.jbs
     open(unit=1,file='out.vgen.jbs',status='replace')
     if(neo_sim_model_in == 5) then
        write(1,'(a)') '#'
        write(1,'(a)') '# expro_rho'
        write(1,'(a)') '# sum z*pflux_neo/(c_s n_e) /(rho_s/a_norm)**2'
        write(1,'(a)') '# jbs_neo    (MA/m^2)'
        write(1,'(a)') '# jbs_sauter_mod (MA/m^2)'
        write(1,'(a)') '# jsigma_neo    (MS/m)'
        write(1,'(a)') '# jsigma_sauter_mod (MS/m)'
        write(1,'(a)') '# jtor_neo    (MA/m^2)'
        write(1,'(a)') '# jtor_sauter_mod (MA/m^2)' 
        write(1,'(a)') '# where jbs = < j_parallel B > / B_unit'
        write(1,'(a)') '# where jtor = < j_tor/R > / <1/R>'
        write(1,'(a)') '#'
        do i=1,EXPRO_n_exp
           write(1,'(8(1pe14.7,2x))') EXPRO_rho(i), pflux_sum(i), &
                jbs_neo(i), jbs_sauter_mod(i), jsigma_neo(i), jsigma_sauter_mod(i), &
                jtor_neo(i), jtor_sauter_mod(i)
        enddo
     else
        write(1,'(a)') '#'
        write(1,'(a)') '# expro_rho'
        write(1,'(a)') '# sum z*pflux_neo/(c_s n_e) /(rho_s/a_norm)**2'
        write(1,'(a)') '# jbs_neo    (MA/m^2)'
        write(1,'(a)') '# jbs_sauter (MA/m^2)'
        write(1,'(a)') '# jsigma_neo    (MS/m)'
        write(1,'(a)') '# jsigma_sauter (MS/m)'
        write(1,'(a)') '# jtor_neo    (MA/m^2)'
        write(1,'(a)') '# jtor_sauter (MA/m^2)' 
        write(1,'(a)') '# where jbs = < j_parallel B > / B_unit'
        write(1,'(a)') '# where jtor = < j_tor/R > / <1/R>'
        write(1,'(a)') '#'
        do i=1,EXPRO_n_exp
           write(1,'(8(1pe14.7,2x))') EXPRO_rho(i), pflux_sum(i), &
                jbs_neo(i), jbs_sauter(i), jsigma_neo(i), jsigma_sauter(i), &
                jtor_neo(i), jtor_sauter(i)
        enddo
     endif
     close(1)
     !----------------------------------------------------------------------

     call vgen_getgeo()

  endif

  deallocate(er_exp)
  deallocate(vtor_measured)
  deallocate(jbs_neo)
  deallocate(jbs_sauter)
  deallocate(jsigma_neo)
  deallocate(jsigma_sauter)
  deallocate(jtor_neo)
  deallocate(jtor_sauter)
  deallocate(pflux_sum)
  deallocate(i_glob)

  cpu_tot_out = MPI_Wtime()

  if (timing_flag == 1) then
     print *, 'total cpu time=', cpu_tot_out-cpu_tot_in
  endif

  call MPI_finalize(i_err)

10 format(&
       'rho=',f6.4,2x,&
       'Er_0(kV/m)=',1pe9.2,2x,&
       'vpol_1(km/s)=',1pe9.2,2x,&
       'nth=',i2,2x,'[',i2,']')

end program vgen
