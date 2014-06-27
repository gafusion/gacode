
  SUBROUTINE nfreya_namelist_defaults
! ----------------------------------------------------------------------
! --- set defaults before reading namelist:
! ------------------------------------------------------------HSJ-------
    USE nrtype,                ONLY  : DP,I4B,I2B

    USE common_constants,      ONLY  : izero,zeroc

 
    USE MPI_data,              ONLY  : parallel_model

    USE solcon_gcnmp,          ONLY  : time0,time_max,dt


    USE  io_gcnmp,             ONLY  : rpc,iterdb_12_input_file,               &
                                       switch_iterdb_output,                   & ! this
                                       switch_statefile_output,                & ! and this are equivalent
                                       statefile_output_name

    USE ions_gcnmp,            ONLY :  ni_sc

    USE source_terms_gcnmp,    ONLY  : beam_th_mult, stfuse_mult, sbfuse_mult, &
                                       qrad_mult

    USE grid_class,            ONLY  : nj,nj_out,nj_max,nj_start,allow_regrid, &
                                       use_compact_schemes

    USE neutral_beams,         ONLY  :                                         &
                                       timbplt, beam_on, beam_time, nameb, relnub,  &
                                        anglev, angleh,nashape, aheigh, &
                                       awidth, bcur, bptor, blenp, nbshape,    &
                                       bleni,bheigh, bwidth, bhfoc, bvfoc,     &
                                       bhdiv, bvdiv, ebkev, fbcur, nbeams,     &
                                       naptr, alen, bvofset, bhofset, nsourc,  &
                                       sfrac1, npart, npskip, rpivot,          &
                                       zpivot, randomize_seed, fionx,fd_beam,mstate,  &
                                       ilorent,iexcit,ngl, kdeni, kdenz, ksvi, &
                                       ksvz, ksve, krad, ngh,ncont,kdene,ne_tk,&
                                       fe_tk,ds_tk,izstrp,iz,                  &
                                       iborb,iterate_beam,npart_all_beamlines, &
                                       no_injectors,beam_sim_time_start,       &
                                       beam_sim_time_end,                      &
                                       nfreya_plot_file_name,namelist_nameb,   &
                                       namelist_nbion,split_injectors,         &
                                       write_performance_data,nfreya_vb,       &
                                       calc_variance,use_ufile

    USE grid_class,            ONLY :  fit_grid,set_cap1
                                       
    USE zonal_data,            ONLY :  mf

    USE xsct,                  ONLY :  adas_xsct_path
 
    USE nf_param,              ONLY :  kb,kprim,kimp

    USE P_nfreya_interface,    ONLY : beam_data_namelist,beam_data_ufile,d_fast_ion

    USE nub,                   ONLY : hdepsmth,fidiff_on,bfr_neutrlz

    IMPLICIT NONE
    INTEGER(I4B) i,j
  
    IF(.NOT. ASSOCIATED(namelist_nameb))ALLOCATE(namelist_nameb(1))
    namelist_nameb(1) = 'none'
    namelist_nbion =1 

    time0    = -1._DP  ! a value <0.0 means that tGCNMs given in iterdb file
                       ! will be used for the start time of the analysis
    time_max = -1._DP  ! a value <0.0 means that tGCNMf given in iterdb file
                       ! will be used for the end  time of the analysis


    !---------------------------------------------------------------------------------
    ! we set defaults for fast ion diffusion here instead if in nubeam defaults section
    ! because nubeam namelist may not be used but default settings are required either way
    !----------------------------------------------------------------------------------
     d_fast_ion%adiff_a(:,:) = -1._DP ; d_fast_ion%adiff_0(:,:) = -1._DP 
     d_fast_ion%adiff_xpin(:,:) = HUGE(1._DP) ; d_fast_ion%adiff_xpout(:,:) =  HUGE(1._DP)  
     d_fast_ion%fidif_on = -1  ! will be set to 1 (sub check_d_fast_ion_inpt) if activated
                               ! in nubeam namelist


    split_injectors = .TRUE. !if true nad ncpus > nbeams then create additional
                             !beams with reduced number of particles
                             ! if false create additional beams with same  number of particles
                             ! as beamline that was copied
    beam_data_namelist  = 'nubeam type namelist file not given'
    beam_data_ufile     = 'ufile name not given'
    use_ufile           = .FALSE.
    use_compact_schemes = izero
    write_performance_data = 'NONE'
    switch_iterdb_output = 0
    switch_statefile_output =  .FALSE.
    statefile_output_name   =  'P_nfreya_statefile'
    nj_out              = nj
    nj_start            = nj	
    allow_regrid        = .TRUE.
    nj_max              = 257   ! must be 2^n+1 for integer n
    mf                  = 20    ! number of grid points defining zones
                                ! for Nfreya zonal scheme

     fit_grid             = .FALSE.
     iterate_beam         = .FALSE.
     calc_variance        = .FALSE.
     iterdb_12_input_file =' '
     set_cap1 = izero
     nfreya_vb = .FALSE.       ! turn off copious output to screen
     hdepsmth = -1.            ! turn of smoothing
     fidiff_on = .FALSE.       ! turn of smoothing by fast ion diffusion
     bfr_neutrlz = .TRUE.
                                          
! ----------------------------------------------------------------------        
! --- NEUTRAL BEAM HEATING PARAMETERS                                           
! ----------------------------------------------------------------------        
!     The default is for co-injection. To get ctr injection                     
!     set angleh to a negative value (usually -19.5 degrees)                    
!     and switch the sources (sfrac(1) becomes sfrac(2), etc.)                  
!     Be sure to check the graphi! output from program NUBPLT if you            
!     change the default parameters!                                            
!                                                                               
! iterate_beam  logical, if .true., allow up to imaxbeamconvg iterations       
!               default: iterate_beam  = .false.                               
!               Iterate_beam = false is set internally (and user input value is  
!               ignored) if time_dep_beam = 1 is set or if nubeam is selected.   
! imaxbeamconvg iterations are done to get consistency                         
!                in the thermal and fast ion densities (ONLY).                  
!                default: imaxbeamconvg =  5                                    
! relaxden      relaxation parameter for beam iteration,must be in (0,1]       
!                with 1 meaning no relaxation.                                  
! relaxden_err  relative error for beam convergence (The relative error is     
!                measured in the change in electron density if icenez=0         
!                and as a relative change in the thermal ion species            
!                corresponding to the beam if icenez=1)                         
! timbplt       Times (up to ntimbplt) to produce data for Freya-like plots           
!                of beam deposition. output is processed by the NUBPLT code.    
!                Defaulted to OFF.                                              
!                timbplt(1) .le. time0 .and. beam_on .lt. time0 gives o/p        
!                for initial time. 

!                If  time_dep_beam = 0 only index i=1 will be used in the       
!                following!!
! beam_on(i)     Time (sec) at which beam heating is turned on for all sources
!                of beam line i                                                 
! beam_time(i)      Time interval (sec) during which beam heating is on
!                this is the pulse length for each pulse from all sources       
!                associated with beam line i
!  external_beam_cur = 0,1, use this switch to input an external beam driven
!                current as a function of space and time.
!                default =0 , which  calculates the beam driven current internally
!                if external_beam_cur = 1 then the electron drag,the
!                trapped electron correction to the electron drag and the
!                pure ion current  (curbe,curbet and curbi) are not calculated
!                individually and hence will appear as 0.0 in the
!                output.
!                The time and space dependent input quantities for the transp
!                beam drive current are given by
!                         knotscurb_external(1...kbctim)
!                         bparcurb_external(kbctim)
!                         rcurb_external(ksplin,kbctim)
!                         curbeam_external(ksplin,kbctim) 
!                where the first index,ksplin, ranges over knots and
!                second index, kbctime ranges over times in bctime (this
!                structure is identical to the other profiles inputs such as
!                tein, described above)
!                
!   -----            NEW INPUTS USED ONLY IF time_dep_beam =1 --------
! beamoff(i)     the length of time                                            
!                all the sources of beam i are off before they turn on again    
!                Note that all sources of a given beam line are controlled by   
!                beam_on(i),beamoff(i) ,beam_time(i) ,and beam_end(i)                
!                Obviously single pulse operation is achieved by setting either 
!                beam_time or beamoff to a value larger than timmax                 
!                (HSJ fix sub getpow for pulsed beams )                         ###
!  beam_end(i)   the time beyond which beam i is shutdown (e.g. never comes on again)
!  beam_cycles(i)An alternative to inputting beam_end is to give the number of
!                cycles that the beam goes through. This is not an integer,
!                fractionla parts of a cycle are allowed. A cycle is defined as the
!                sum of the on and off time (eg. beam_time(i)+beamoff(i) )
!                if both beam_end and beam_cycles is set then beam_end is used
! source2_phase(i)  This is the time lead (if source2_phase(i) is input as      
!                a negative number) or time lag ((if source2_phase(i) is input  
!                as a positive number). Of source 2, beamline i, relative to    
!                source 1,beamline i in seconds.                                
!                (if source2_phase(1) = -0.050 for example it means that        
!                 source 2 of beam line 1 is started 50 mse! before source 1    
!                 of this beam line starts. The start time of source 1 is       
!                 given by beam_on(1)                                            
!                 if the input nsour! is set to 1 instead of 2 then this has    
!                 no effect. By making source2_phase(i) a number larger than    
!                 timmax you can effectively eliminate the second source of     
!                 any beam line.
!  beam_thermal_cutoff flag for determination of lower limit on slowing down    
!                 distribution of fast ions, set to -1,0 or 1, 1 is default     
!                 If   beam_thermal_cutoff = 1 then the thermal speed           
!                of the ions is used as the lower limit of integration          
!                (see beam_pulse_control below).                                
!                if beam_thermal_cutoff =0 then the lower limit of integration  
!                over the fast ion     
!                distribution function is taken as  therm_frac* vthermal ion    
!                therm_fra! is DEFAULTED TO  1.0                                
!                Both beam_thermal_cutoff=0 and beam_thermal_cutoff=1 will      
!                result in a lower limit that depends on Ti and hence rho.      
!                NOTE that in either case the Fokker-Planck equation used to    
!                determine the slowing down distribution is not valid    at     
!                these energies and hence this is an approximation
!                To get a cutoff independent of rho set beam_thermal_cutoff =   
!                -1 and input the value of the thermal ion speed to be used     
!                in  beam_thermal_speed ( see below)                            
!                The other factors that go into the integration limit (ne,te)   
!                are averaged over the grid, using a straight average at        
!                present. see tau0_avg subroutine                               
!  beam_mode     character variable set to "test" or "run"                      
!                ( "run" is default) in test mode Freya is not called. this     
!                 mode is intended primarily to build up waveforms of beam      
!                input before doing a transport run for developers only.        
! beam_restart_file     A file, if present in the working directory, 
!                       will be used to set beam initial conditions, 
!                       provided time0 matches the  input values in inone 
!                       and beam_init_restart_file = 0.
!                       ( beam_restart_file can hold up to 256 characters 
!                         - no check or error exit at this time )
!                       Note that this is NOT the restart file for Nubeam !!!!!
! beam_init_restart_file   
!                Set to 1 to create file beam_restart_file in this run.
!                (default 0) if file exists it will be replaced by a new one.


! ngauss         the number of gaussian quadrature points to be used in         
!                carrying out the integrations. This is a accuracy versus       
!                execution speed issue. This value probably does not have to    
!                be changed by the causual user. The maximum value is           
!                ngauss_max (see  gauss_info.i)                                 

! beam_pulse_control  This integer valued switch determines how the upper and   
!                lower limits of integration on the fast ion distribution       
!                function (and all required moments thereof) are handled in     
!                determination of the fundamental time step dt.
!                if beam_thermal_cutoff = -1 then beam_pulse_control is forced  
!                to be 0 ( no matter what the user set it to)

!                Besides the     
!                obvious fact that all pulse on and off times must be           
!                reflected in the time                                          
!                steps taken by the code we also have to consider those times   
!                at which the upper and lower limits of integration reach       
!                the saturated value of the thermal cutoff.
!                But this thermal cutoff normally depends on the local          
!                temperature and hence is different at every grid point.        
!                Thus we see that dt will have to be modified many times to     
!                account for each beam,injection energy, source,grid point and  
!                pulse number. This means that many time steps will be          
!                generated internally. This is what will happen if              
!                beam_pulse_control is left at its default value of 0           
!                and beam_thermal_cutoff = 0 or 1.

!                Set  beam_pulse_control = -1 if you      want to skip the      
!                modification of the time step due to this effect altogether    
!                In this case there will be some error in the calculated        
!                beam integrals (because we wont hit the thermalization time    
!                exactly)The amount of this error is then controlled            
!                by the maximum allowed time step, dtmax. If you select         
!                beam_pulse_control = -1 then you should cut dtmax  down from   
!                its default value appropriately.                               

!                beam_pulse_control = n, where n >0 and n < nj, ( nj            
!                is the number of rho grid points)  means that the plasma is    
!                broken up into n regions, with the beam thermalization         
!                time calulated at the            WARNING not yet do
!                mid point of each region. This means that instead of nj mods   
!                of the time step for each beam,energy,source,pulse, there will 
!                only be n such interruptions. What we have in mind here is that 
!                n=1,2,3,4,5, et! but any n up to nj-1 can be used.             
!                For n=1 for example the slowing down parmeters will be         
!                evaluated at nj/2 and held constant at that value at all       
!                grid points assuming beam_thermal_cutoff =0 or 1                

!                Note that  beam_thermal_cutoff = -1 is the same as setting     
!                beam_pulse-control = 1 but  using beam_thermal_speed for the   
!                single value to be used as the lower limit of integration      
!                at all the rho grid points.

!    
! beam_thermal_speed          in cm/sec, used when beam_pulse_control =0        
!                the default value corresponds to an ion temperature of 50ev  

! -----            END of       NEW INPUTS USED ONLY IF time_dep_beam =1 ----      
!           
!                
!  ibcur         Flag for neutral beam driven current                           
!                  1, include beam driven current                               
!                  0, neglect beam driven current                               
!   use_Callen    integer,has meaning only if mcgo run is done.                 
!                 Intended for checking out some model details.                 
!                 Users should leave it set at the default value of 0.          
!  ibcx          Flag for secondary charge exchange between fast ions           
!                  and thermal neutrals                                         
!                  1, include secondary charge exchange (default)               
!                  0, neglect secondary charge exchange                         
!                  NOTE: thermal ion particle source due to secondary charge    
!                        exchange is not included in the ion density equation,  
!                        independent of the setting of ibcx.                    
!                                                                               
!  nbeamtcx    switch used to determine the form of the beam torque.            
!              if nbeamtcx = 0 (default) then the loss of beam torque due       
!              to secondary charge exchange of fast ions with thermal           
!              neutrals is neglected. if nbeamtcx = 1 the transfer of angular   
!              momentum from fast ions to thermal ions and electrons is         
!              modified to account for the cx losses (done with array ssprcxl)  
!  ibslow      Flag for neutral beam slowing down                               
!                1, include neutral beam slowing down                           
!                0, neglect neutral beam slowing down                           
!                If the new multiple beam pulse model is selected               
!                (time_dep_beam =1) then ibslow =1 is enforced by the code      
!  fionx       Allows testing for non-classical slowing down                    
!                  (see subroutine SLOW1).                                      
!                                                                               
!  fast_ion_target  integer,if set to 1 will let the incoming neutral beam      
!                   see the stored fast ion density                             
!                   no attempt to modify the stopping rates due to the          
!                   fact that the fast ion distribution is not Maxwellian       
!                   has been made!!!!!!!!!!!!!!!!!!!!!!!                        
!                   Note that if the stored fast ion density is large enough    
!                   that it is of concern then it is also large enough so that  
!                   a simple one pass linearization doesnt make sense. So you   
!                   should use this option in conjunction with the iterate      
!                   beam option.                                                
!                   For Mcgo coupled runs the birth points of the fast ions     
!                   are determined by the Freya in Onetwo. Hence this switch    
!                   will affect the Mcgo results as well. Note that Mcgo        
!                   itself does not let the Monte Carlo particles see the       
!                   stored fast ion density. As an approximate remedy for       
!                   this situation you can set fast_ion_target = -1 when        
!                   running a Mcgo coupled case. This will add the stored fast  
!                   ion density determined in Freya to the thermal ion density  
!                   used in Mcgo, so that the Mcgo particles see a correct      
!                   total density.                                              
!                   (If Mcgo is run with fast_ion_target = 1 then               
!                   the birth points feed to Mcgo will account for the stored   
!                   fast ions but subsequent slowing down of the Monte Carlo    
!                   particles  will see only the thermal ions)                  
!  rtstcx       A factor between zero and 1 used to fudge secondary charge      
!               of fast ions. ie the original fast ion charge exchanges         
!               with a thermal neutral. The resulting fast neutral              
!               will however most likely be reionized. To model this            
!               reionization in Onetwo we don't do the secondary charge         
!               exchange calculation, which requires a Monte Carlo approach.    
!               Instead we assume that the probability against charge exchange  
!               of a fast ion, normally taken as                                
!               {(vbeam**3+vcrit**3)/(v**3+vcrit**3)}**A                        
!               where   A=(-taus/(3.*taucx)) is the ratio of slowing down       
!               to charge exchange lifetime is given instead by                 
!               {(vbeam**3+vcrit**3)/(v**3+vcrit**3)}**B                        
!               where B=rtstcx*A. To decrease charge exchange losses            
!               set rtstcx less than 1. Note that this can be interpreted       
!               as assuming that the neutral density, enn, is replaced by       
!               rtstcx*enn for the purposes of charge exchange probability      
!               calculations ONLY. The actual neutral density is not changed!   
!                                                                               
!               Set rtstcx to a negative number to use sigma*v, rather          
!               than <sigma*v> (i.e., Maxwellian average), to determine         
!               the mean lifetime of a fast ion against charge exchange.        
!               The energy and velocity used will be that of the beam           
!               corrected for bulk rotation of ions (note the assumption        
!               that neutrals rotate as ions do).                               
!               If rtstcx is set to a number less than -5. then                 
!               the average fast ion energy (again corrected for rotation)      
!               for each beam slowing down distribution is used.                
!                                                                               
!  nameb        Name of neutral species in beam                                 
!              'h', 'd', 't', 'dt'                                              
!               MUST BE PRIMARY ION SPECIES  #1 or #2.                          
!               if nameb = 'h' then fdbeam defaults to natural isotopi!         
!                  content of d in 'h'.                                         
!               if nameb =  't' fdbeam is explicitly set to 0                   
!               if nameb = 'dt' fdbeam is set to fd, (the fraction              
!               of d in thermal dt mixture)  USER HAS NO CHOICE                 
!                                            UNLESS IFUS = -1 is selected       
!               NOTE that nameb=dt selects a SINGLE fast ion fluid              
!               with a fictitious mass !!!!                                     
!  fdbeam      Fraction of deuterium atoms in neutral beam.                     
!  relnub      Minimum relative change in ion density or electron               
!              temperature between neutral beam calculations                    
!              (0.10 is suggested)                                              
!  tfusbb      'thermal fusion balance for beams', fraction by which            
!              the net energy gain from thermal fusion must exceed              
!              the net energy loss for automati! beam turnoff.                  
!              If tfusbb = 0, automati! beam turnoff is not done.               
!  iddcal     Flag controlling treatment of beam effects on fuscal,             
!             the calculated fusion neutron rate:                               
!             0 = do not include knock-on or beam-d neutrons in fuscal          
!             1 = include only knock-on neutrons in fuscal                      
!             2 = include only beam-d neutrons in fuscal                        
!             3 = include both knock-on and beam-d neutrons in fuscal           
!  randomize_seed  set to true to give a romdom starting seed to each processor
!             used. The seed will be different on each processor.
!             if false the seed will still be different on each processor but
!             it will not be randomized. Instead the same starting seed will
!             will be used each time the code is run.   NOTE: the
!             option randomize_seed = .TRUE. is not available if the old
!             nfreya random number gneraor is used! In that case each
!             processor gets adifferent  seed but no randomization is done.

!  npart      Number of particles followed into plasma (suggest 10000)          
!  npskip     Ratio of number of particles followed into plasma                 
!               to number of source particles (suggest 1)                       
!             npskip =1 is set in Nfreya if this is a parallel run with more    
!             than one processor                                                
!  iborb      Flag for modeling orbit effects on beam-generated fast ions       
             
!             1 => do     model orbit effects (this is the original method)     
!             0 => do not model orbit effects                                   
!  itrapfi    Flag for trapped fast ions.  If iborb = 1 then                    
!             setting itrapfi = 1 will modify the beam driven                   
!             current by the initial trapped ion fraction.                      
!             (Subsequent pitch angle diffusion is NOT taken into account).     
!             itrapfi = 0 neglects this effect.  If iborb = 0,                  
!             itrapfi has no effect.  itrapfi = 0 is default.                   
!  ds_tk     maximum trajectory increment (cm) used in subroutine INJECT to     
!            calculate psi(s), where s is the neutral trajectory                
!            pathlength from the first closed flux surface it                   
!            encounters.  used for non-zero toroidal rotation, where            
!            mean free path as a function of path length is required.           
!  fe_tk     factor (>1) to set upper limit of energy in n*sigma array.         
!            max(ebins) = max(ebkev)*fe_tk.  required for nonzero               
!            rotation cases. fe_tk=1.3 is default                               
!  ne_tk     number of equi-width energy bins used in forming n*sigma           
!            array.  required for nonzero rotation cases.  is internally        
!            reset to zero (used as flag to turn off rotational effects         
!            on neutral stopping) if angular rotation is not present            
!            (iangrot = 0).                                                     
!  iexcit    Selection switch for atomi! cross section data and model:          
!            0,  Use the fundamental atomi! data of Freeman & Jones (1972).     
!                  This option is not advised since this atomi! data is         
!                  considered outdated and excited state effects are            
!                  not considered.                                              
!            1,  Use hexnb package but do not include excitations in            
!                  its calculation of cross sections.                           
!            2,  Use hexnb package, include excitations in calculation          
!                  of cross sections.                                           
!                  NOTE: Options 1 and 2 are not advised because the hexnb      
!                  model has been found to be flawed and is based on            
!                  atomi! data that is considered outdated.                     
!                  (Boley et al, Nuclear Fusion 29, 1984)                       
!            5,  Use the JET-ADAS effective stopping cross sections.            
!                  This model is preferred since it is based on the most        
!                  recent atomi! data available and includes multi-step         
!                  ionization processes due to excited states. This is an       
!                  important effect for considering neutral beam penetration    
!                  into fusion-grade plasmas. For a detailed description        
!                  of the ADAS Atomi! Data and Analysis Structure (ADAS)        
!                  developed by JET, see Finkenthal, 1994.                      
!                  Note: All ions are considered fully stripped.                
!                  (Daniel Finkenthal Ph.D. Thesis, U! Berkeley, 1994)          
!                                                                               
!  neg_ion_source(i)  i=1,...,nbeams                                            
!                     an integer array, set to 1 to indicate                    
!                     that a NEGATIVE ion source is used for                    
!                     neutral beam line i. At this time, the only               
!                     effect of this switch is to set the neutralization        
!                     efficiency, arbitrarily, to 98%, independent of the       
!                     negative ions energy and to eliminate the second and      
!                     third energy components from the beam.                    
!  time_dep_beam      set to 1 to indicate that new time dependent beam input   
!                     will be used. (default = 0).(NOTE: if this is a snapshot  
!                     run then the code will set time_dep_beam =0 even if user  
!                     selected time_dep_beam =1)
!                                                                               
!  izstrp(i) i = 1,2..nimp set to 0 for coronal equilibrium                     
!            values of <z> and <zsq>. set to 1 for fully stripped impurities    
!  mstate    principal quantum number n above which excitations                 
!            count as ionizations.                                              
!  inubpat   two-dimensional beam deposition option                             
!            0,  do not calculate beam deposition in (r,z) coordinates          
!            1,  calculate beam deposition on (r,z) grid (default = eqdsk       
!                grid).  write deposition array and n = 3 excited state         
!                fraction to file 'beamdep' for standalone analysis.            
!  npat      neutral beam deposition (r,z) grid dimensions, used if             
!            inubpat = 1.0                                                      
!            npat(1)  =  number of elements in 'r' (<=2*nw)                     
!            npat(2)  =  number of elements in 'z' (<=2*nh)                     
!            defaults, npat(1) = nw, npat(2) = nh                               
!  mf        Number of flux zones plus 1                             
!                                                                               
!  In the following list the index ib designates the beam injector,             
!  while ie refers to one of the three energy components.                       
!  iap refers to one of the apertures.                                          
!                                                                               
!  nbeams          Number of neutral beam injectors (1 to 2 )               
!  nsour!          Number of sources per beamline.                              
!                    If 1, source is centered on beamline axis.                 
!                    If nsour! = 2, distinguish between the beamline            
!                    axis and the source centerline (optical axis).             
!                    The two sources are assumed to be mirror images            
!                    through the beamline axis.                                 
!                    In either case, the exit grid plane is perpendicular       
!                    to the beamline axis, and contains the source              
!                    exit grid center(s).                                       
!                    If nsour! = 2, the alignment of the sources w.r.t.         
!                    the beamline axis is specified through bhofset,            
!                    bvofset, and bleni (described further below).              
!  naptr           Total number of apertures encountered by a particle          
!                    as is moves from the source into the plasma chamber.       
!                    Maximum is specified by parameter nap ( = 10).             
!                    First set of apertures encountered by the particles        
!                    are assumed centered on the source axis, and subsequent    
!                    apertures are centered on the beamline axis;               
!                    the distinction is made through ashape.                    
!  anglev(ib)      Vertical angle (degrees) between optical axis                
!                    and horizontal plane; a positive value indicates           
!                    particles move upward                                      
!  angleh(ib)      Horizontal angle (degrees) between optical axis and          
!                    vertical plane passing through pivot point and             
!                    toroidal axis; a zero value denotes perpendicular          
!                    injection, while a positive value indicates par-           
!                    ticles move in the co-current direction                    
!  bvofset(ib)     Vertical offset from beamline axis to center                 
!                    of each source (cm; used only for nsour! = 2)              
!  bhofset(ib)     Horizontal offset from beamline axis to center               
!                    of each source (cm; used only for nsour! = 2)              
!  bleni(ib)       Length along source centerline (source optical axis) from    
!                    source to intersection point with the beamline axis.       
!  sfrac1(ib)      Fraction of source current per beamline coming               
!                    from upper source (used only for nsour! = 2)               
!                    (the upper source is the more perpendicular,               
!                      or right source normally)                                
!  bcur(ib)        Total current (a) in ion beam (used only if bptor            
!                    is zero)                                                   
!  bptor(ib)       Total power (w) through aperture into torus; when            
!                    nonzero, bptor takes precedence over bcur                  
!  nbshape(ib)      Beam shape                                                   
!                    'circ':  circular                                          
!                    'rect':  rectangular                                       
!                'rect-lps':  rect. long pulse source (DIII-D only)             
!                             a choice of short or long pulse sources is        
!                             available by injector (not by source).  one       
!                             or both injectors may be long pulse by            
!                             setting one or both to 'rect-lps'                 
!                  CAUTION:  DIII-D sources are defaulted to lps specs.         
!                  It is the user's responsibility to overide these for         
!                   sps configuration(s).                                       
!  bheigh(ib)      Height of source (cm)                                        
!                             Default is bshape(1) = bshape(2) = 'rect-lps'     
!  bwidth(ib)      Width of source (cm); diameter for circular source.          
!  bhfo! (ib)      Horizontal focal length of source (cm)                       
!  bvfo! (ib)      Vertical focal length of source (cm)                         
!  bhdiv (ib)      Horizontal divergence of source (degrees)                    
!  bvdiv (ib)      Vertical divergence of source (degrees)                      
!  ebkev (ib)      Maximum particle energy in source (keV)                      
!  fbcur (ie,ib)   Fraction of current at energy ebkev/ie                       
!                  Note that this is the current fraction at the source,        
!                  before it enters the neutralizer.                            
!  ashape(iap,ib)  Aperture shape.                                              
!                   Prefix 's-' indicates source axis centered.                 
!                   Prefix 'b-' indicates beamline axis centered.               
!                     's-circ'          'b-circ'                                
!                     's-rect'          'b-rect'                                
!                     's-vert'          'b-vert'                                
!                     's-horiz'         'b-horiz'                               
!                                       'b-d3d'                                 
!                    cir! = circular aperture,                                  
!                    rect = rectangular,                                        
!                    vert = limits vertical height of source particles,         
!                   horiz = limits horizontal height of source particles,       
!                     d3d = special DIII-D polygonal aperture                   
!  aheigh(iap,ib)  Height of aperture (cm)                                      
!  awidth(iap,ib)  Width  of aperture (cm); diameter for circular               
!                    aperture                                                   
!  alen(iap,ib)    Length from source to aperture for 's-type' apertures,       
!                    and from exit grid plane along beamline axis for           
!                    'b-type' apertures.                                        
!  blenp(ib)       Distance along beamline axis from source exit                
!                    plane to the fiducial "pivot" point.                       
!  rpivot(ib)      Radial position of pivot (cm)                                
!  zpivot(ib)      Axial position of pivot (cm)                                 
!                                                                               
!  hdepsmth    set this param to a positive value (gt.0.0 and .le.  10) to turn off          
!              the smoothing of hibrz and hdepz in subroutine POSTNUB.          
!              if this option is used then enough zones must be specified       
!              for adequate resolution (zones = number of radial grid points)   
!              and enough injected neutrals must be followed to minimize        
!              the statistical noise enough so that no greatly uneven           
!              profiles result! this option was added because the smoothing     
!              of the profiles by subroutine SMOOTH can lead to unphysical      
!              peaking of the birth and deposition profiles.                    
!              Matching of Monte Carlo results for fast ion distributions,      
!              especially in the presence of mhd activity, indicates that       
!              smoothing with hdepsmth may be required. Hence                   
!              if hdepsmth .gt. 10 then ipass = hdepsmth -10 passes are made    
!              over the data to get a smooth profile. Each pass averages        
!              nhdep = 6 (not adjustable) grid points together. Note in         
!              particular that if the grid is coarse then 6 points cover        
!              a wider range in rho than if the grid were finer. hence varying  
!              the grid size  probably also means adjusting hdepsmth to         
!              maintain a more or less                                          
!              constant deposition profile. Note that the smoothed profiles     
!              are renormalized to the plasma volume.                           
!              It is possible to account for the fixed nhdep by adjusting       
!              hdepsmth. As an example it was observed that with 51 grid        
!              points and hdepsmth =45 a fast ion deposition profile resulted   
!              which  was subsequently matched using 201 grid points by         
!              changing hdepsmth from 45 to 600. (mf = 12 and typical DIII-D    
!              densities for h mode shots)                                      
!                             
! beam_sim_time_start
! beam_sim_time_end  start and end times of beam simulation for this run.
!                                                  
! nfreya_vb     logical  switch used for debug purposes                           
!               if  nfreya_vb is TRUE  then some FREYA-related diagnostic           
!               output will be written to the screen                             
!                                                                               
!     DIII-D beam input.  DEFAULT IS LONG-PULSE-SOURCE SPECIFICATIONS.          
!                                                                               
      naptr = 4                                                                 
!                                                                               
      do i=1,kb                                                                 
        anglev(i)    =   0.0                                                    
        angleh(i)    =  19.5                                                    
        bvofset(i)   =   0.0                                                    
        bhofset(i)   =  42.074                                                  
        bleni(i)     = 556.808                                                  
        bcur(i)      = 110.0   ! dont change,see logic in freya
        bptor(i)     =  0.001                                                     
        nbshape(i)   = 'rect-lps'                                               
        bheigh(i)    =  48.0                                                    
        bwidth(i)    =  12.0                                                    
        bhdiv(i)     =   0.50       !degrees                                    
        bvdiv(i)     =   1.3        !degrees                                    
        fbcur(1,i)   =   0.7                                                    
        fbcur(2,i)   =   0.2                                                    
        fbcur(3,i)   =   0.1                                                    
        bhfoc(i)     =   1.0d100                                                
        bvfoc(i)     =   1.0e3                                                  
        ebkev(i)     =  75.0                                                    
        ebkev(i)     =   0.0001         ! changed HSJ                           
        sfrac1(i)    =   0.5                                                    
        nashape(1,i) = 's-rect'                                                 
        nashape(2,i) = 's-rect'                                                 
        nashape(3,i) = 'b-d3d'                                                   
        nashape(4,i) = 'b-circ'                                                 
        aheigh(1,i)  =  47.8                                                    
        awidth(1,i)  =  13.8                                                    
        alen(1,i)    = 186.1                                                    
        aheigh(2,i)  =  48.0                                                    
        awidth(2,i)  =  17.7                                                    
        alen(2,i)    = 346.0                                                    
        alen(3,i)    = 449.0                                                    
        awidth(4,i)  =  50.9                                                    
        alen(4,i)    = 500.0                                                    
        blenp(i)     = 539.0                                                    
        rpivot(i)    = 286.6                                                    
        zpivot(i)    =   0.0                                                    
      end do                                                                    
!                                                                               
!  parameters for including rotation in neutral stopping                        
!                                                                               
 9030 ds_tk = 5.0                                                               
      fe_tk = 1.4  !passed to nbsgxn where it is used as ebfac                  
      ne_tk = 20   !passed to nbsgxn where it is used as nebin   
                   ! but check ksge (nf_param) for consistency.
                   ! currently not done.               
!                                                                               
! --- following parameters are used in subroutine HEXNB                         
!                   
      ilorent =  0 
      mstate  =  4
      ncont   =  30
      kdene   =  1                                                              
      kdeni   =  1                                                              
      kdenz   =  1                                                              
      ksvi    =  0                                                              
      ksvz    =  0                                                              
      ksve    =  0                                                              
      krad    =  1                                                              
      ngh     = 10                                                              
      ngl     = 10                                                              
      iexcit  =  5 ! note this is 0 in onetwo                                          
      adas_xsct_path = '/usc-data/p2/linux/onetwo/' 
      iborb   = 1
      ilorent =  0                                                              
      mstate  =  4                                                              
      ncont   = 30                                                              
                             
      npart = -1  ;  npart_all_beamlines = npart ! total # pseudo neutrals
      nbeams = -1 ;  no_injectors = nbeams        ! total # injectors(eg beamlines)
      DO  j=1,kimp                                                          
         iz(j)     = izero                                                          
         izstrp(j) = izero                                                             
!                                                                               
!      note izstrp = 0 implies coronal equilibrium for impurity j in hexnb      
!                                                                               

      ENDDO                                                      
      beam_sim_time_start = -HUGE(1._DP) ; beam_sim_time_end = beam_sim_time_start



      nfreya_plot_file_name = 'P_NF_bpltfil'
      IF(.NOT. ASSOCIATED(ni_sc)) ALLOCATE(ni_sc(4))

    RETURN


  END SUBROUTINE nfreya_namelist_defaults
