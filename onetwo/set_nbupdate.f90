subroutine set_nbupdate

  use nbi_types
  USE transp
  USE solcon,only : time,dt
  implicit NONE



  !---------------------
  ! update NBI input items 
  ! 
  ! this covers the minority of items that can be changed from timestep
  ! to timestep:  beam powers, voltages, plasma parameter profiles, etc.
  !
  ! see "nbspec.dat" which defines which input quantities are updatable,
  ! and which are not.
  !
  ! generally, the call "nbi_get_xxx" fetches the current settings
  !            and call "nbi_update_xxx" applies the new values
  !
  ! only "updatable" quantities will be modified.  If some updatable
  ! quantities are *not* to be modified or reset explicitly, then, the
  ! "nbi_get_xxx" call is required, to be sure that the correct values
  ! of such quantities are not lost during the "nbi_update_xxx" call.

  type(nbitype_times)  :: ztimes       ! basic system
  type(nbitype_grid)   :: zgrid        ! computational grids
  type(nbitype_powers)  :: zpowers     ! beam powers, voltages, energy fractns
  type(nbitype_fusion) :: zfusion      ! fusion product species
  type(nbitype_num)    :: znum         ! numerical controls
  type(nbitype_saw)    :: zsaw         ! sawtooth controls
  type(nbitype_adif)   :: zadif        ! anomolous diffusion
  type(nbitype_atomic) :: zatomic      ! atomic physics
  type(nbitype_collid) :: zcollid      ! collision operator
  type(nbitype_minority) :: zminority  ! RF minority specie
  type(nbitype_ripple) :: zripple      ! TF rippple model controls
  type(nbitype_fishbone) :: zfbone     ! "fishbone" model controls
  type(nbitype_misc)   :: zmisc        ! miscellany

  !--------------------------------------------------------------

  integer id_ripple,id_mcgrid,nbeam_tr
  real *8 time_nubeam
  !--------------------------------------------------------------
  !  timestep information
  !  (nbi_get_times call not necessary)
  !  tbm1,tbm2 are set in nubeam_12 to time and time + dt
  ztimes%tbm1 = tbm1  ! start of timestep (time, seconds)
  ztimes%tbm2 = tbm2  ! end of timestep (time, seconds)

  call nbi_update_times(ztimes)
     print *,'sub set_nbupdate,tbm1,tbm2 =',tbm1,tbm2

  !--------------------------------------------------------------
  !  update xplasma id:  Monte Carlo 2d grid
  !
  call xpload_mcgrid_gen      !added 2/13/04 HSJ
  call nbi_get_grid(zgrid)    ! previous settings

  call xpload_mcgrid(id_mcgrid)
  zgrid%id_mcgrid = id_mcgrid ! current xplasma id for 2d MC irregular grid

  call nbi_update_grid(zgrid) ! apply update (id_mcgrid)
     print *,'sub set_nbupdate,done nbi_update_grid'





  !--------------------------------------------------------------
  !  update beam powers, voltages, and energy fractions
  !  (nbi_get_powers not used; all powers, etc. are specified explicitly).

  nbeam_tr = beam_data%nbeam
   time_nubeam = 0.5*(tbm2+tbm1)      !use values at midpoint
   time_nubeam = time + dt            ! = tbm2 to plot correctly need tbm1

   call beam_power_interp(time_nubeam)

  zpowers%einja(1:nbeam_tr) = beam_data%einja(1:nbeam_tr)     ! voltages
  zpowers%pinja(1:nbeam_tr) = beam_data%pinja(1:nbeam_tr)     ! powers
  zpowers%ffulla(1:nbeam_tr) = beam_data%ffulla(1:nbeam_tr)   ! full energy fractions
  zpowers%fhalfa(1:nbeam_tr) = beam_data%fhalfa(1:nbeam_tr)   ! half energy fractions

  !  for details on definitions of ffulla & fhalfa, see nbi_com_mod.f90
  !  comments...

  call nbi_update_powers(zpowers) ! apply update in NUBEAM
     print *,'no beams ',nbeam_tr
     print *,'beam  powers',beam_data%pinja(1:nbeam_tr)
     print *,'sub set_nbupdate,done nbi_update_powers'


  !--------------------------------------------------------------
  !  Fusion product species... He3, He4(alphas), Tritons, Protons
  !
  !  update source rates, where requested.

  call nbi_get_fusion(zfusion) ! values currently in NUBEAM

  ! (because nbi_get_fusion is used, some of the following could be
  ! omitted safely...)

  zfusion%nlusf3 = nlusf3      ! switches to supply birth rate data
  zfusion%nlusfa = nlusfa      ! otherwise use internal fusion rate calculation
  zfusion%nlusft = nlusft
  zfusion%nlusfp = nlusfp

  zfusion%xdatsf3 = xdatsf3    ! the birthrate data given here
  zfusion%xdatsfa = xdatsfa    ! (ignored unless corresponding switch is set)
  zfusion%xdatsft = xdatsft
  zfusion%xdatsfp = xdatsfp

  call nbi_update_fusion(zfusion) ! apply update in NUBEAM
     print *,'sub set_nbupdate,done nbi_update_fusion'

  !--------------------------------------------------------------
  !  "num" block -- numerical controls and adjustments
  !  this block allows TRANSP to modify some numerical controls
  !  time dependently...
  !

  call nbi_get_num(znum)    ! *** fetch previous settings ***
                            ! This call is important, because not all
                            ! updatable values are being explicitly set.
                            ! This call makes sure that values from the
                            ! previous timestep are reapplied for the
                            ! next timestep.

  znum%nptcls = nptcls      ! *** number of MC model ions for BEAM SPECIES ***
  znum%nptclf = nptclf      ! *** number of MC ions for FUSION PRODUCTS ***

  call nbi_update_num(znum) ! apply update...
     print *,'no MC ions for beam species ',nptcls
     print *,'no MC ions for fusion  products ',nptclf
     print *,'sub set_nbupdate,done nbi_update_num'

  !--------------------------------------------------------------
  !  atomic physics operator

  call nbi_get_atomic(zatomic)

  zatomic%dn0out = dn0out   ! "external" neutral density

  call nbi_update_atomic(zatomic)
     print *,'sub set_nbupdate,done nbi_update_atomic'

  !--------------------------------------------------------------
  !  collision operator

  call nbi_get_collid(zcollid)

!  zcollid%nlorbo = nlorbo   ! "orbit only" switch ONLY FOR DEBUGGING...

  call nbi_update_collid(zcollid)

  !--------------------------------------------------------------
  !  RF minority

  call nbi_get_minority(zminority)

  zminority%nthrf = nthrf   ! =1 if RF minority is thermalized.

  call nbi_update_minority(zminority)

  !--------------------------------------------------------------
  !  sawtooth options

  call nbi_get_saw(zsaw)  ! previous settings...

  if ((tsaw-tbm2).le.dtmint) then
     zsaw%nlsawb = nlsawb      ! .TRUE. for sawtooth mixing of beam ions
     zsaw%nlsawf = nlsawf      ! .TRUE. for sawtooth mixing of fusion products
  else
     zsaw%nlsawb = .FALSE.     ! no sawtooth expected on this step...
     zsaw%nlsawf = .FALSE.
  endif

  call nbi_update_saw(zsaw)    ! apply updates...
     print *,'sub set_nbupdate,done nbi_update_saw'

  !--------------------------------------------------------------
  !  anomolous diffusion options (to activate the model, provide an
  !  anomolous diffusivity profile... see set_nbprofiles(...))

  call nbi_get_adif(zadif)   ! previous settings...

  zadif%nkdifb = nkdifb      ! =3: apply to all Monte Carlo fast ions
                             ! =1: beam ions only; =2: fusion products only.

  zadif%ndifbe = ndifbe      ! energy dependent adjustment (see nbspec.dat).
  zadif%fdifbe(1:ndifbe) = fdifbe(1:ndifbe)
  zadif%edifbe(1:ndifbe) = edifbe(1:ndifbe)

  call nbi_update_adif(zadif)  ! apply updates...

     print *,'sub set_nbupdate,done nbi_update_adif'

  !--------------------------------------------------------------
  !  ripple model:
  !   (a) taurip, loss time for simple model, can be updated.
  !   (b) id_ripple, xplasma id for log(TF_RIPPLE_FIELD) = f(R,Z)
  !  see nbspec.dat or nbi_com_mod.f90 for details
  !

  call nbi_get_ripple(zripple)  ! previous settings...

  call eq_gfnum('LOG_TF_RIPPLE',id_ripple)

  zripple%taurip = taurip
  zripple%id_ripple = id_ripple

  call nbi_update_ripple(zripple)  ! apply update...
      print *,'sub set_nbupdate,done nbi_update_ripple'


  !--------------------------------------------------------------
  ! fishbone model:  can be activated/de-activated over time

  call nbi_get_fishbone(zfbone)    ! previous settings

  zfbone%nlfbon = nlfbon

  call nbi_update_fishbone(zfbone) ! apply update
     print *,'sub set_nbupdate,done nbi_update_fishbone'






  !--------------------------------------------------------------
  ! misc block -- update plasma elevation, if this is an updown
  ! symmetric run

  if(nlsym2b) then

     call nbi_get_misc(zmisc)      ! previous settings

     zmisc%plhgt = plhgt

     call nbi_update_misc(zmisc)   ! apply updates

  endif
     print *,'sub set_nbupdate,done nbi_update_misc'
  !---------------------

end subroutine set_nbupdate
