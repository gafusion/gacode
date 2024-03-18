module expro

  ! List of all useful interface objects
  character*12, dimension(120) :: expro_list 

  character(len=2) :: ident='# '
  double precision :: expro_mass_deuterium=3.34358e-24  ! md (g)
  
  !-----------------------------------------------------------
  ! DATA:
  !
  != Header entries

  logical :: hasmpi
  integer :: expro_comm
  integer :: expro_n_exp
  integer :: expro_n_ion
  integer :: expro_shot=0
  integer :: expro_time=0
  integer, parameter :: expro_ion_max=200
  character*10, dimension(expro_ion_max) :: expro_name
  character*10, dimension(expro_ion_max) :: expro_type

  double precision :: expro_masse=5.44887413e-4   ! me/m_H
  double precision, dimension(:), allocatable :: expro_mass
  double precision :: expro_ze=-1.0
  double precision, dimension(:), allocatable :: expro_z
  double precision :: expro_torfluxa=0.0
  double precision :: expro_rcentr=0.0
  double precision :: expro_bcentr=0.0
  double precision :: expro_current=0.0

  != 1D and 2D profile arrays contained in input.gacode

  double precision, dimension(:), allocatable :: expro_rho
  double precision, dimension(:), allocatable :: expro_rmin
  double precision, dimension(:), allocatable :: expro_polflux
  double precision, dimension(:), allocatable :: expro_q
  double precision, dimension(:), allocatable :: expro_w0
  double precision, dimension(:), allocatable :: expro_rmaj
  double precision, dimension(:), allocatable :: expro_zmag
  double precision, dimension(:), allocatable :: expro_kappa
  double precision, dimension(:), allocatable :: expro_delta
  double precision, dimension(:), allocatable :: expro_zeta
  double precision, dimension(:), allocatable :: expro_shape_cos0
  double precision, dimension(:), allocatable :: expro_shape_cos1
  double precision, dimension(:), allocatable :: expro_shape_cos2
  double precision, dimension(:), allocatable :: expro_shape_cos3
  double precision, dimension(:), allocatable :: expro_shape_cos4
  double precision, dimension(:), allocatable :: expro_shape_cos5
  double precision, dimension(:), allocatable :: expro_shape_cos6
  double precision, dimension(:), allocatable :: expro_shape_sin3
  double precision, dimension(:), allocatable :: expro_shape_sin4
  double precision, dimension(:), allocatable :: expro_shape_sin5
  double precision, dimension(:), allocatable :: expro_shape_sin6
  double precision, dimension(:), allocatable :: expro_ne
  double precision, dimension(:,:), allocatable :: expro_ni
  double precision, dimension(:), allocatable :: expro_te
  double precision, dimension(:,:), allocatable :: expro_ti
  double precision, dimension(:), allocatable :: expro_ptot
  double precision, dimension(:), allocatable :: expro_fpol
  double precision, dimension(:), allocatable :: expro_johm
  double precision, dimension(:), allocatable :: expro_jbs
  double precision, dimension(:), allocatable :: expro_jrf
  double precision, dimension(:), allocatable :: expro_jnb
  double precision, dimension(:), allocatable :: expro_jbstor
  double precision, dimension(:), allocatable :: expro_sigmapar
  double precision, dimension(:), allocatable :: expro_z_eff
  double precision, dimension(:,:), allocatable :: expro_vpol
  double precision, dimension(:,:), allocatable :: expro_vtor
  double precision, dimension(:), allocatable :: expro_qohme
  double precision, dimension(:), allocatable :: expro_qbeame
  double precision, dimension(:), allocatable :: expro_qbeami
  double precision, dimension(:), allocatable :: expro_qrfe
  double precision, dimension(:), allocatable :: expro_qrfi
  double precision, dimension(:), allocatable :: expro_qfuse
  double precision, dimension(:), allocatable :: expro_qfusi
  double precision, dimension(:), allocatable :: expro_qbrem
  double precision, dimension(:), allocatable :: expro_qsync
  double precision, dimension(:), allocatable :: expro_qline
  double precision, dimension(:), allocatable :: expro_qei
  double precision, dimension(:), allocatable :: expro_qione
  double precision, dimension(:), allocatable :: expro_qioni
  double precision, dimension(:), allocatable :: expro_qcxi
  double precision, dimension(:), allocatable :: expro_qpar_beam
  double precision, dimension(:), allocatable :: expro_qpar_wall
  double precision, dimension(:), allocatable :: expro_qmom

  != 0D Derived quantities
  double precision :: expro_betap
  double precision :: expro_betat
  double precision :: expro_betan
  double precision :: expro_greenwald
  double precision :: expro_ptransp
  double precision :: expro_tau
  double precision :: expro_tau98y2 

  != 1D Derived quantities
  double precision, dimension(:), allocatable :: expro_bunit
  double precision, dimension(:), allocatable :: expro_gamma_e
  double precision, dimension(:), allocatable :: expro_gamma_p
  double precision, dimension(:), allocatable :: &
       expro_s,&
       expro_drmaj,&
       expro_dzmag,&
       expro_skappa,&
       expro_sdelta,&
       expro_szeta,&
       expro_shape_scos0,&
       expro_shape_scos1,&
       expro_shape_scos2,&
       expro_shape_scos3,&
       expro_shape_scos4,&
       expro_shape_scos5,&
       expro_shape_scos6,&
       expro_shape_ssin3,&
       expro_shape_ssin4,&
       expro_shape_ssin5,&
       expro_shape_ssin6,&
       expro_dlnnedr,&
       expro_dlntedr,&
       expro_sdlnnedr,&
       expro_sdlntedr,&
       expro_dlnptotdr,&
       expro_w0p,&
       expro_surf,&
       expro_vol,&
       expro_volp,&
       expro_cs,&
       expro_rhos,&
       expro_nuee,&
       expro_ni_new,&
       expro_dlnnidr_new,&
       expro_sdlnnidr_new,&
       expro_grad_r0,&
       expro_ave_grad_r,&
       expro_bp0,&
       expro_bt0,&
       expro_bp2,&
       expro_bt2,&
       expro_mach,&
       expro_thetascale,&
       expro_flow_beam,&
       expro_flow_wall,&
       expro_flow_mom,&
       expro_pow_e,&
       expro_pow_i,&
       expro_pow_ei,&
       expro_pow_e_aux,&
       expro_pow_i_aux,&
       expro_pow_e_fus,&
       expro_pow_i_fus,&
       expro_pow_e_sync,&
       expro_pow_e_brem,&
       expro_pow_e_line,&
       expro_pow_e_ohmic

  != 2D Derived quantities
  double precision, dimension(:,:), allocatable :: &
       expro_dlnnidr,&
       expro_dlntidr,&
       expro_sdlnnidr,&
       expro_sdlntidr
  !-----------------------------------------------------------

  ! input.gacode.geo dimension and arrays

  integer :: expro_nfourier
  double precision, dimension(:,:,:), allocatable :: expro_geo
  double precision, dimension(:,:,:), allocatable :: expro_dgeo

  ! Field orientation parameters

  integer :: expro_signb
  integer :: expro_signq

  ! Control parameters

  integer :: expro_ctrl_n_ion = -1
  integer :: expro_ctrl_quasineutral_flag 
  integer :: expro_ctrl_numeq_flag
  integer :: expro_error=0
  integer :: expro_jerr=0

  ! Header information
  character(len=70) :: expro_head_original =  '#  *original : null'
  character(len=70) :: expro_head_statefile = '# *statefile : null'
  character(len=70) :: expro_head_gfile =     '#     *gfile : null'
  character(len=70) :: expro_head_cerfile =   '#   *cerfile : null'
  character(len=70) :: expro_head_vgen =      '#      *vgen : null'
  character(len=70) :: expro_head_tgyro =     '#     *tgyro : null'

contains

  subroutine expro_init(flag)

    implicit none

    integer :: nexp,nion
    integer, intent(in) :: flag
        
    if (flag == 1) then

       nexp = expro_n_exp
       nion = expro_n_ion

       allocate(expro_mass(nion))    ; expro_mass = 1.0
       allocate(expro_z(nion))       ; expro_z = 1.0

       allocate(expro_rho(nexp))     ; expro_rho = 0.0
       allocate(expro_rmin(nexp))    ; expro_rmin = 0.0
       allocate(expro_q(nexp))       ; expro_q = 0.0
       allocate(expro_polflux(nexp)) ; expro_polflux = 0.0
       allocate(expro_w0(nexp))      ; expro_w0 = 0.0
       allocate(expro_rmaj(nexp))    ; expro_rmaj = 0.0

       allocate(expro_zmag(nexp))    ; expro_zmag = 0.0
       allocate(expro_kappa(nexp))   ; expro_kappa = 0.0
       allocate(expro_delta(nexp))   ; expro_delta = 0.0
       allocate(expro_zeta(nexp))    ; expro_zeta = 0.0
       allocate(expro_shape_cos0(nexp))   ; expro_shape_cos0 = 0.0
       allocate(expro_shape_cos1(nexp))   ; expro_shape_cos1 = 0.0
       allocate(expro_shape_cos2(nexp))   ; expro_shape_cos2 = 0.0
       allocate(expro_shape_cos3(nexp))   ; expro_shape_cos3 = 0.0
       allocate(expro_shape_cos4(nexp))   ; expro_shape_cos4 = 0.0
       allocate(expro_shape_cos5(nexp))   ; expro_shape_cos5 = 0.0
       allocate(expro_shape_cos6(nexp))   ; expro_shape_cos6 = 0.0
       allocate(expro_shape_sin3(nexp))   ; expro_shape_sin3 = 0.0
       allocate(expro_shape_sin4(nexp))   ; expro_shape_sin4 = 0.0
       allocate(expro_shape_sin5(nexp))   ; expro_shape_sin5 = 0.0
       allocate(expro_shape_sin6(nexp))   ; expro_shape_sin6 = 0.0
       allocate(expro_ne(nexp))      ; expro_ne = 0.0
       allocate(expro_te(nexp))      ; expro_te = 0.0
       allocate(expro_ptot(nexp))    ; expro_ptot = 0.0
       allocate(expro_johm(nexp))    ; expro_johm = 0.0
       allocate(expro_jbs(nexp))     ; expro_jbs = 0.0
       allocate(expro_jrf(nexp))     ; expro_jrf = 0.0
       allocate(expro_jnb(nexp))     ; expro_jnb = 0.0
       allocate(expro_jbstor(nexp))  ; expro_jbstor = 0.0
       allocate(expro_sigmapar(nexp)); expro_sigmapar = 0.0
       allocate(expro_z_eff(nexp))   ; expro_z_eff = 0.0

       allocate(expro_qohme(nexp))  ; expro_qohme = 0.0
       allocate(expro_qbeame(nexp)) ; expro_qbeame = 0.0
       allocate(expro_qbeami(nexp)) ; expro_qbeami = 0.0
       allocate(expro_qrfe(nexp))   ; expro_qrfe = 0.0
       allocate(expro_qrfi(nexp))   ; expro_qrfi = 0.0
       allocate(expro_qfuse(nexp))  ; expro_qfuse = 0.0
       allocate(expro_qfusi(nexp))  ; expro_qfusi = 0.0
       allocate(expro_qbrem(nexp))  ; expro_qbrem = 0.0
       allocate(expro_qsync(nexp))  ; expro_qsync = 0.0
       allocate(expro_qline(nexp))  ; expro_qline = 0.0
       allocate(expro_qei(nexp))    ; expro_qei = 0.0
       allocate(expro_qione(nexp))  ; expro_qione = 0.0
       allocate(expro_qioni(nexp))  ; expro_qioni = 0.0
       allocate(expro_qcxi(nexp))   ; expro_qcxi = 0.0
       allocate(expro_qpar_beam(nexp))   ; expro_qpar_beam = 0.0
       allocate(expro_qpar_wall(nexp))   ; expro_qpar_wall = 0.0
       allocate(expro_qmom(nexp))   ; expro_qmom = 0.0

       allocate(expro_ni(nion,nexp))   ; expro_ni = 0.0
       allocate(expro_ti(nion,nexp))   ; expro_ti = 0.0
       allocate(expro_vpol(nion,nexp)) ; expro_vpol = 0.0
       allocate(expro_vtor(nion,nexp)) ; expro_vtor = 0.0

       ! Derived quantities

       allocate(expro_bunit(nexp))        ; expro_bunit = 0.0
       allocate(expro_s(nexp))            ; expro_s = 0.0
       allocate(expro_drmaj(nexp))        ; expro_drmaj = 0.0
       allocate(expro_dzmag(nexp))        ; expro_dzmag = 0.0
       allocate(expro_skappa(nexp))       ; expro_skappa = 0.0
       allocate(expro_sdelta(nexp))       ; expro_sdelta = 0.0
       allocate(expro_szeta(nexp))        ; expro_szeta = 0.0
       allocate(expro_shape_scos0(nexp))  ; expro_shape_scos0 = 0.0
       allocate(expro_shape_scos1(nexp))  ; expro_shape_scos1 = 0.0
       allocate(expro_shape_scos2(nexp))  ; expro_shape_scos2 = 0.0
       allocate(expro_shape_scos3(nexp))  ; expro_shape_scos3 = 0.0
       allocate(expro_shape_scos4(nexp))  ; expro_shape_scos4 = 0.0
       allocate(expro_shape_scos5(nexp))  ; expro_shape_scos5 = 0.0
       allocate(expro_shape_scos6(nexp))  ; expro_shape_scos6 = 0.0
       allocate(expro_shape_ssin3(nexp))  ; expro_shape_ssin3 = 0.0
       allocate(expro_shape_ssin4(nexp))  ; expro_shape_ssin4 = 0.0
       allocate(expro_shape_ssin5(nexp))  ; expro_shape_ssin5 = 0.0
       allocate(expro_shape_ssin6(nexp))  ; expro_shape_ssin6 = 0.0
       allocate(expro_dlnnedr(nexp))      ; expro_dlnnedr = 0.0
       allocate(expro_dlntedr(nexp))      ; expro_dlntedr = 0.0
       allocate(expro_sdlnnedr(nexp))     ; expro_sdlnnedr = 0.0
       allocate(expro_sdlntedr(nexp))     ; expro_sdlntedr = 0.0
       allocate(expro_dlnptotdr(nexp))    ; expro_dlnptotdr = 0.0
       allocate(expro_w0p(nexp))          ; expro_w0p = 0.0
       allocate(expro_surf(nexp))         ; expro_surf = 0.0
       allocate(expro_vol(nexp))          ; expro_vol = 0.0
       allocate(expro_volp(nexp))         ; expro_volp = 0.0
       allocate(expro_cs(nexp))           ; expro_cs = 0.0
       allocate(expro_rhos(nexp))         ; expro_rhos = 0.0
       allocate(expro_nuee(nexp))         ; expro_nuee = 0.0
       allocate(expro_ni_new(nexp))       ; expro_ni_new = 0.0
       allocate(expro_dlnnidr_new(nexp))  ; expro_dlnnidr_new = 0.0
       allocate(expro_sdlnnidr_new(nexp)) ; expro_sdlnnidr_new = 0.0
       allocate(expro_grad_r0(nexp))      ; expro_grad_r0 = 0.0
       allocate(expro_ave_grad_r(nexp))   ; expro_ave_grad_r = 0.0
       allocate(expro_bp0(nexp))          ; expro_bp0 = 0.0
       allocate(expro_bt0(nexp))          ; expro_bt0 = 0.0
       allocate(expro_bp2(nexp))          ; expro_bp2 = 0.0
       allocate(expro_bt2(nexp))          ; expro_bt2 = 0.0
       allocate(expro_fpol(nexp))         ; expro_fpol = 0.0
       allocate(expro_gamma_e(nexp))      ; expro_gamma_e = 0.0
       allocate(expro_gamma_p(nexp))      ; expro_gamma_p = 0.0
       allocate(expro_mach(nexp))         ; expro_mach = 0.0
       allocate(expro_thetascale(nexp))   ; expro_thetascale = 0.0
       allocate(expro_flow_beam(nexp))  ; expro_flow_beam = 0.0
       allocate(expro_flow_wall(nexp))  ; expro_flow_wall = 0.0
       allocate(expro_flow_mom(nexp))   ; expro_flow_mom = 0.0
       allocate(expro_pow_e(nexp))      ; expro_pow_e = 0.0
       allocate(expro_pow_i(nexp))      ; expro_pow_i = 0.0
       allocate(expro_pow_ei(nexp))     ; expro_pow_ei = 0.0
       allocate(expro_pow_e_aux(nexp))  ; expro_pow_e_aux = 0.0
       allocate(expro_pow_i_aux(nexp))  ; expro_pow_i_aux = 0.0
       allocate(expro_pow_e_fus(nexp))  ; expro_pow_e_fus = 0.0
       allocate(expro_pow_i_fus(nexp))  ; expro_pow_i_fus = 0.0
       allocate(expro_pow_e_sync(nexp)) ; expro_pow_e_sync = 0.0
       allocate(expro_pow_e_brem(nexp)) ; expro_pow_e_brem = 0.0
       allocate(expro_pow_e_line(nexp)) ; expro_pow_e_line = 0.0
       allocate(expro_pow_e_ohmic(nexp)) ; expro_pow_e_ohmic = 0.0

       allocate(expro_dlnnidr(nion,nexp))  ; expro_dlnnidr = 0.0
       allocate(expro_dlntidr(nion,nexp))  ; expro_dlntidr = 0.0
       allocate(expro_sdlnnidr(nion,nexp)) ; expro_sdlnnidr = 0.0
       allocate(expro_sdlntidr(nion,nexp)) ; expro_sdlntidr = 0.0

    else

       deallocate(expro_mass) 
       deallocate(expro_z) 

       deallocate(expro_rho)
       deallocate(expro_rmin)
       deallocate(expro_q)
       deallocate(expro_polflux)
       deallocate(expro_w0)
       deallocate(expro_rmaj)
       deallocate(expro_zmag)
       deallocate(expro_kappa)
       deallocate(expro_delta)
       deallocate(expro_zeta)
       deallocate(expro_shape_cos0)
       deallocate(expro_shape_cos1)
       deallocate(expro_shape_cos2)
       deallocate(expro_shape_cos3)
       deallocate(expro_shape_cos4)
       deallocate(expro_shape_cos5)
       deallocate(expro_shape_cos6)
       deallocate(expro_shape_sin3)
       deallocate(expro_shape_sin4)
       deallocate(expro_shape_sin5)
       deallocate(expro_shape_sin6)
       deallocate(expro_ne)
       deallocate(expro_te)
       deallocate(expro_ptot)
       deallocate(expro_johm)
       deallocate(expro_jbs)
       deallocate(expro_jrf)
       deallocate(expro_jnb)
       deallocate(expro_jbstor)
       deallocate(expro_sigmapar)
       deallocate(expro_z_eff)

       deallocate(expro_qohme)
       deallocate(expro_qbeame)
       deallocate(expro_qbeami)
       deallocate(expro_qrfe)
       deallocate(expro_qrfi)
       deallocate(expro_qfuse)
       deallocate(expro_qfusi)
       deallocate(expro_qbrem)
       deallocate(expro_qsync)
       deallocate(expro_qline)
       deallocate(expro_qei)
       deallocate(expro_qione)
       deallocate(expro_qioni)
       deallocate(expro_qcxi)
       deallocate(expro_qpar_beam)
       deallocate(expro_qpar_wall)
       deallocate(expro_qmom)

       deallocate(expro_ni,expro_ti,expro_vpol,expro_vtor)

       ! Derived

       deallocate(expro_bunit)
       deallocate(expro_s)
       deallocate(expro_drmaj)
       deallocate(expro_dzmag)       
       deallocate(expro_skappa)
       deallocate(expro_sdelta)
       deallocate(expro_szeta)
       deallocate(expro_shape_scos0)
       deallocate(expro_shape_scos1)
       deallocate(expro_shape_scos2)
       deallocate(expro_shape_scos3)
       deallocate(expro_shape_scos4)
       deallocate(expro_shape_scos5)
       deallocate(expro_shape_scos6)
       deallocate(expro_shape_ssin3)
       deallocate(expro_shape_ssin4)
       deallocate(expro_shape_ssin5)
       deallocate(expro_shape_ssin6)
       deallocate(expro_dlnnedr)      
       deallocate(expro_dlntedr)      
       deallocate(expro_sdlnnedr)      
       deallocate(expro_sdlntedr)      
       deallocate(expro_dlnptotdr)    
       deallocate(expro_w0p)          
       deallocate(expro_surf)          
       deallocate(expro_vol)          
       deallocate(expro_volp)         
       deallocate(expro_cs)           
       deallocate(expro_rhos)         
       deallocate(expro_nuee)
       deallocate(expro_ni_new)       
       deallocate(expro_dlnnidr_new)  
       deallocate(expro_sdlnnidr_new)
       deallocate(expro_grad_r0)      
       deallocate(expro_ave_grad_r)   
       deallocate(expro_bp0)          
       deallocate(expro_bt0)          
       deallocate(expro_bp2)          
       deallocate(expro_bt2)          
       deallocate(expro_fpol)           
       deallocate(expro_gamma_e)      
       deallocate(expro_gamma_p)      
       deallocate(expro_mach)   
       deallocate(expro_thetascale)  
       deallocate(expro_flow_beam)
       deallocate(expro_flow_wall)
       deallocate(expro_flow_mom)
       deallocate(expro_pow_e)
       deallocate(expro_pow_i)
       deallocate(expro_pow_ei)
       deallocate(expro_pow_e_aux)
       deallocate(expro_pow_i_aux)
       deallocate(expro_pow_e_fus)
       deallocate(expro_pow_i_fus)
       deallocate(expro_pow_e_sync)
       deallocate(expro_pow_e_brem)
       deallocate(expro_pow_e_line)
       deallocate(expro_pow_e_ohmic)

       deallocate(expro_dlnnidr)
       deallocate(expro_dlntidr)  
       deallocate(expro_sdlnnidr) 
       deallocate(expro_sdlntidr)  

    endif

  end subroutine expro_init

  subroutine expro_read(thisinfile,comm)

    implicit none

    character(len=*), intent(in) :: thisinfile
    integer, intent(in) :: comm
    integer :: nexp,nion,ierr,i,nd
    character(len=70) :: ytag
    character(len=22) :: c

    if (comm == 0) then
       hasmpi = .false.
    else
       hasmpi = .true.
       expro_comm = comm
    endif

    ! ORDERING NOTE: nexp should appear before any profile arrays
    
    open(unit=1,file=trim(thisinfile),status='old')

    do

       read(1,'(a)',end=99) ytag

       if (index(ytag,'*original') > 0) then
          expro_head_original=ytag ; cycle
       else if (index(ytag,'*statefile') > 0) then
          expro_head_statefile=ytag ; cycle
       else if (index(ytag,'*gfile') > 0) then
          expro_head_gfile=ytag ; cycle
       else if (index(ytag,'*cerfile') > 0) then
          expro_head_cerfile=ytag ; cycle
       else if (index(ytag,'*vgen') > 0) then
          expro_head_vgen=ytag ; cycle
       else if (index(ytag,'*tgyro') > 0) then
          expro_head_tgyro=ytag ; cycle
       endif

       nd = scan(ytag,'|')
       if (nd == 0) then
          ! no units field so trim all whitespace
          c = trim(ytag(3:))
       else
          ! trim units and whitespace
          c = trim(ytag(3:nd-1))
       endif

       select case (c)
       case ('nexp')
          call expro_icomm(expro_n_exp) 
          nexp = expro_n_exp
       case ('nion')
          call expro_icomm(expro_n_ion)
          nion = expro_n_ion
          if (allocated(expro_rho)) call expro_init(0)
          call expro_init(1) 
       case ('shot')
          call expro_icomm(expro_shot) 
       case ('time')
          call expro_icomm(expro_time) 
       case ('name')
          call expro_tcomm(expro_name(1:nion),nion)
       case ('type')
          call expro_tcomm(expro_type(1:nion),nion)
       case ('masse')
          call expro_rcomm(expro_masse) 
       case ('mass')
          call expro_lcomm(expro_mass,nion)
       case ('ze')
          call expro_rcomm(expro_ze) 
       case ('z')
          call expro_lcomm(expro_z,nion)
       case ('torfluxa')
          call expro_rcomm(expro_torfluxa) 
       case ('rcentr')
          call expro_rcomm(expro_rcentr) 
       case ('bcentr')
          call expro_rcomm(expro_bcentr) 
       case ('current')
          call expro_rcomm(expro_current) 
       case ('rho')
          call expro_vcomm(expro_rho,nexp)  
       case ('rmin')
          call expro_vcomm(expro_rmin,nexp)  
       case ('polflux')
          call expro_vcomm(expro_polflux,nexp)  
       case ('q')
          call expro_vcomm(expro_q,nexp)  
       case ('w0')
          call expro_vcomm(expro_w0,nexp)  
       case ('rmaj')
          call expro_vcomm(expro_rmaj,nexp)  
       case ('zmag')
          call expro_vcomm(expro_zmag,nexp)  
       case ('kappa')
          call expro_vcomm(expro_kappa,nexp)  
       case ('delta')
          call expro_vcomm(expro_delta,nexp)  
       case ('zeta')
          call expro_vcomm(expro_zeta,nexp)
       case ('shape_cos0')
          call expro_vcomm(expro_shape_cos0,nexp)
       case ('shape_cos1')
          call expro_vcomm(expro_shape_cos1,nexp)
       case ('shape_cos2')
          call expro_vcomm(expro_shape_cos2,nexp)
       case ('shape_cos3')
          call expro_vcomm(expro_shape_cos3,nexp)
       case ('shape_cos4')
          call expro_vcomm(expro_shape_cos4,nexp)
       case ('shape_cos5')
          call expro_vcomm(expro_shape_cos5,nexp)
       case ('shape_cos6')
          call expro_vcomm(expro_shape_cos6,nexp)
       case ('shape_sin3')
          call expro_vcomm(expro_shape_sin3,nexp)
       case ('shape_sin4')
          call expro_vcomm(expro_shape_sin4,nexp)
       case ('shape_sin5')
          call expro_vcomm(expro_shape_sin5,nexp)
       case ('shape_sin6')
          call expro_vcomm(expro_shape_sin6,nexp) 
       case ('ne')
          call expro_vcomm(expro_ne,nexp) 
       case ('te')
          call expro_vcomm(expro_te,nexp)  
       case ('ptot')
          call expro_vcomm(expro_ptot,nexp)  
       case ('fpol')
          call expro_vcomm(expro_fpol,nexp)  
       case ('johm')
          call expro_vcomm(expro_johm,nexp)  
       case ('jbs')
          call expro_vcomm(expro_jbs,nexp)  
       case ('jrf')
          call expro_vcomm(expro_jrf,nexp)  
       case ('jnb')
          call expro_vcomm(expro_jnb,nexp)  
       case ('jbstor')
          call expro_vcomm(expro_jbstor,nexp)  
       case ('sigmapar')
          call expro_vcomm(expro_sigmapar,nexp)  
       case ('z_eff')
          call expro_vcomm(expro_z_eff,nexp) 
       case ('ni')
          call expro_acomm(expro_ni(:,:),nion,nexp) 
       case ('ti')
          call expro_acomm(expro_ti(:,:),nion,nexp) 
       case ('vpol')
          call expro_acomm(expro_vpol(:,:),nion,nexp) 
       case ('vtor')
          call expro_acomm(expro_vtor(:,:),nion,nexp) 
       case ('qohme')
          call expro_vcomm(expro_qohme,nexp) 
       case ('qbeame')
          call expro_vcomm(expro_qbeame,nexp) 
       case ('qbeami')
          call expro_vcomm(expro_qbeami,nexp) 
       case ('qrfe')
          call expro_vcomm(expro_qrfe,nexp) 
       case ('qrfi')
          call expro_vcomm(expro_qrfi,nexp) 
       case ('qfuse')
          call expro_vcomm(expro_qfuse,nexp) 
       case ('qfusi')
          call expro_vcomm(expro_qfusi,nexp) 
       case ('qbrem')
          call expro_vcomm(expro_qbrem,nexp) 
       case ('qsync')
          call expro_vcomm(expro_qsync,nexp) 
       case ('qline')
          call expro_vcomm(expro_qline,nexp) 
       case ('qei')
          call expro_vcomm(expro_qei,nexp) 
       case ('qione')
          call expro_vcomm(expro_qione,nexp) 
       case ('qioni')
          call expro_vcomm(expro_qioni,nexp) 
       case ('qcxi')
          call expro_vcomm(expro_qcxi,nexp)
       case ('qpar') ! for backward compatibility
          call expro_vcomm(expro_qpar_beam,nexp)
       case ('qpar_beam')
          call expro_vcomm(expro_qpar_beam,nexp)
       case ('qpar_wall')
          call expro_vcomm(expro_qpar_wall,nexp)
       case ('qmom')
          call expro_vcomm(expro_qmom,nexp) 
       end select

    enddo

99  close(1)

    ! ** input.gacode.geo **

    nexp = expro_n_exp
    open(unit=1,file=trim(thisinfile)//'.geo',status='old',iostat=ierr)
    if (ierr == 0) then
       call expro_skip_header(1)
       call expro_icomm(expro_nfourier)
       if (allocated(expro_geo)) deallocate(expro_geo)
       if (allocated(expro_dgeo)) deallocate(expro_dgeo)
       allocate(expro_geo(4,0:expro_nfourier,nexp)) ; expro_geo(:,:,:)=0.0
       allocate(expro_dgeo(4,0:expro_nfourier,nexp)) ; expro_dgeo(:,:,:)=0.0
       do i=1,nexp
          call expro_scomm(expro_geo(:,:,i),4*(expro_nfourier+1))
       enddo
    else
       expro_nfourier = -1
    endif
    close(1)

    if (expro_ctrl_n_ion <= expro_n_ion) then
       call expro_compute_derived
       call expro_list_set
    else
       expro_error = 1
    endif

  end subroutine expro_read

  subroutine expro_write(thisinfile)

    implicit none

    integer :: i,nexp,nion
    character(len=*), intent(in) :: thisinfile 

    nexp = expro_n_exp
    nion = expro_n_ion

    ! Write header
    open(unit=1,file=trim(thisinfile),status='replace')
    write(1,'(a)') expro_head_original
    write(1,'(a)') expro_head_statefile 
    write(1,'(a)') expro_head_gfile
    write(1,'(a)') expro_head_cerfile
    write(1,'(a)') expro_head_vgen
    write(1,'(a)') expro_head_tgyro
    write(1,'(a)') '#'

    ! Write data
    call expro_writei(nexp,'nexp')
    call expro_writei(nion,'nion')
    call expro_writei(expro_shot,'shot')
    call expro_writei(expro_time,'time')
    write(1,'(a)') ident//'name' ; write(1,'(20(a,1x))') (trim(expro_name(i)),i=1,nion)
    write(1,'(a)') ident//'type' ; write(1,'(20(a,1x))') (trim(expro_type(i)),i=1,nion)
    write(1,'(a)') ident//'masse'; write(1,30) expro_masse
    write(1,'(a)') ident//'mass' ; write(1,40) expro_mass
    write(1,'(a)') ident//'ze'   ; write(1,30) expro_ze
    write(1,'(a)') ident//'z'    ; write(1,40) expro_z

    ! Write vector/array data, skipping objects that are 0.0
    call expro_writes(expro_torfluxa,'torfluxa','Wb/radian')
    call expro_writes(expro_rcentr,'rcentr','m')
    call expro_writes(expro_bcentr,'bcentr','T')
    call expro_writes(expro_current,'current','MA')
    call expro_writev(expro_rho,nexp,'rho','-')
    call expro_writev(expro_rmin,nexp,'rmin','m')
    call expro_writev(expro_polflux,nexp,'polflux','Wb/radian')
    call expro_writev(expro_q,nexp,'q','-')
    call expro_writev(expro_w0,nexp,'w0','rad/s')
    call expro_writev(expro_rmaj,nexp,'rmaj','m')
    call expro_writev(expro_zmag,nexp,'zmag','m')
    call expro_writev(expro_kappa,nexp,'kappa','-')
    call expro_writev(expro_delta,nexp,'delta','-')
    call expro_writev(expro_zeta,nexp,'zeta','-')
    call expro_writev(expro_shape_cos0,nexp,'shape_cos0','-')
    call expro_writev(expro_shape_cos1,nexp,'shape_cos1','-')
    call expro_writev(expro_shape_cos2,nexp,'shape_cos2','-')
    call expro_writev(expro_shape_cos3,nexp,'shape_cos3','-')
    call expro_writev(expro_shape_cos4,nexp,'shape_cos4','-')
    call expro_writev(expro_shape_cos5,nexp,'shape_cos5','-')
    call expro_writev(expro_shape_cos6,nexp,'shape_cos6','-')
    call expro_writev(expro_shape_sin3,nexp,'shape_sin3','-')
    call expro_writev(expro_shape_sin4,nexp,'shape_sin4','-')
    call expro_writev(expro_shape_sin5,nexp,'shape_sin5','-')
    call expro_writev(expro_shape_sin6,nexp,'shape_sin6','-')
    call expro_writev(expro_ne,nexp,'ne','10^19/m^3')
    call expro_writea(expro_ni(:,:),nion,nexp,'ni','10^19/m^3')
    call expro_writev(expro_te,nexp,'te','keV')
    call expro_writea(expro_ti(:,:),nion,nexp,'ti','keV')
    call expro_writev(expro_ptot,nexp,'ptot','Pa')
    call expro_writev(expro_fpol,nexp,'fpol','T-m')
    call expro_writev(expro_johm,nexp,'johm','MA/m^2')
    call expro_writev(expro_jbs,nexp,'jbs','MA/m^2')
    call expro_writev(expro_jrf,nexp,'jrf','MA/m^2')
    call expro_writev(expro_jnb,nexp,'jnb','MA/m^2')
    call expro_writev(expro_jbstor,nexp,'jbstor','MA/m^2')
    call expro_writev(expro_sigmapar,nexp,'sigmapar','MSiemens/m')
    call expro_writev(expro_z_eff,nexp,'z_eff','-')
    call expro_writea(expro_vpol(:,:),nion,nexp,'vpol','m/s')
    call expro_writea(expro_vtor(:,:),nion,nexp,'vtor','m/s')
    call expro_writev(expro_qohme,nexp,'qohme','MW/m^3')
    call expro_writev(expro_qbeame,nexp,'qbeame','MW/m^3')
    call expro_writev(expro_qbeami,nexp,'qbeami','MW/m^3')
    call expro_writev(expro_qrfe,nexp,'qrfe','MW/m^3')
    call expro_writev(expro_qrfi,nexp,'qrfi','MW/m^3')
    call expro_writev(expro_qfuse,nexp,'qfuse','MW/m^3')
    call expro_writev(expro_qfusi,nexp,'qfusi','MW/m^3')
    call expro_writev(expro_qbrem,nexp,'qbrem','MW/m^3')
    call expro_writev(expro_qsync,nexp,'qsync','MW/m^3')
    call expro_writev(expro_qline,nexp,'qline','MW/m^3')
    call expro_writev(expro_qei,nexp,'qei','MW/m^3')
    call expro_writev(expro_qione,nexp,'qione','MW/m^3')
    call expro_writev(expro_qioni,nexp,'qioni','MW/m^3')
    call expro_writev(expro_qcxi,nexp,'qcxi','MW/m^3')
    call expro_writev(expro_qpar_beam,nexp,'qpar_beam','1/m^3/s') 
    call expro_writev(expro_qpar_wall,nexp,'qpar_wall','1/m^3/s')
    call expro_writev(expro_qmom,nexp,'qmom','N/m^2')

    close(1)

30  format(1pe14.7)
40  format(10(1pe14.7))

  end subroutine expro_write

! This is the full list of user variable for the expro interface
 
subroutine expro_list_set
  
  expro_list(1) = 'n_exp'
  expro_list(2) = 'n_ion'
  expro_list(3) = 'mass'
  expro_list(4) = 'z'
  expro_list(5) = 'torfluxa'
  expro_list(6) = 'rcentr'
  expro_list(7) = 'bcentr'
  expro_list(8) = 'current'
  expro_list(9) = 'rho'
  expro_list(10) = 'rmin'
  expro_list(11) = 'q'
  expro_list(12) = 'w0'
  expro_list(13) = 'rmaj'
  expro_list(14) = 'zmag'
  expro_list(15) = 'kappa'
  expro_list(16) = 'delta'
  expro_list(17) = 'zeta'
  expro_list(18) = 'shape_cos0'
  expro_list(19) = 'shape_cos1'
  expro_list(20) = 'shape_cos2'
  expro_list(21) = 'shape_cos3'
  expro_list(22) = 'shape_cos4'
  expro_list(23) = 'shape_cos5'
  expro_list(24) = 'shape_cos6'
  expro_list(25) = 'shape_sin3'
  expro_list(26) = 'shape_sin4'
  expro_list(27) = 'shape_sin5'
  expro_list(28) = 'shape_sin6'
  expro_list(29) = 'ne'
  expro_list(30) = 'ni'
  expro_list(31) = 'te'
  expro_list(32) = 'ti'
  expro_list(33) = 'ptot'
  expro_list(34) = 'fpol'
  expro_list(35) = 'johm'
  expro_list(36) = 'jbs'
  expro_list(37) = 'jrf'
  expro_list(38) = 'jnb'
  expro_list(39) = 'jbstor'
  expro_list(40) = 'sigmapar'
  expro_list(41) = 'z_eff'
  expro_list(42) = 'vpol'
  expro_list(43) = 'vtor'
  expro_list(44) = 'qohme'
  expro_list(45) = 'qbeame'
  expro_list(46) = 'qbeami'
  expro_list(47) = 'qrfe'
  expro_list(48) = 'qrfi'
  expro_list(49) = 'qfuse'
  expro_list(50) = 'qfusi'
  expro_list(51) = 'qbrem'
  expro_list(52) = 'qsync'
  expro_list(53) = 'qline'
  expro_list(54) = 'qei'
  expro_list(55) = 'qione'
  expro_list(56) = 'qioni'
  expro_list(57) = 'qcxi'
  expro_list(58) = 'qpar_beam'
  expro_list(59) = 'qpar_wall'
  expro_list(60) = 'qmom'
  expro_list(61) = 'bunit'
  expro_list(62) = 'gamma_e'
  expro_list(63) = 'gamma_p'
  expro_list(64) = 's'
  expro_list(65) = 'drmaj'
  expro_list(66) = 'dzmag'
  expro_list(67) = 'skappa'
  expro_list(68) = 'sdelta'
  expro_list(69) = 'szeta'
  expro_list(70) = 'shape_scos0'
  expro_list(71) = 'shape_scos1'
  expro_list(72) = 'shape_scos2'
  expro_list(73) = 'shape_scos3'
  expro_list(74) = 'shape_scos4'
  expro_list(75) = 'shape_scos5'
  expro_list(76) = 'shape_scos6'
  expro_list(77) = 'shape_ssin3'
  expro_list(78) = 'shape_ssin4'
  expro_list(79) = 'shape_ssin5'
  expro_list(80) = 'shape_ssin6'
  expro_list(81) = 'dlnnedr'
  expro_list(82) = 'dlntedr'
  expro_list(83) = 'dlnnidr'
  expro_list(84) = 'dlntidr'
  expro_list(85) = 'w0p'
  expro_list(86) = 'surf'
  expro_list(87) = 'vol'
  expro_list(88) = 'volp'
  expro_list(89) = 'cs'
  expro_list(90) = 'rhos'
  expro_list(91) = 'nuee'
  expro_list(92) = 'grad_r0'
  expro_list(93) = 'ave_grad_r'
  expro_list(94) = 'bp0'
  expro_list(95) = 'bt0'
  expro_list(96) = 'mach'
  expro_list(97) = 'flow_beam'
  expro_list(98) = 'flow_wall'
  expro_list(99) = 'flow_mom'
  expro_list(100) = 'pow_e'
  expro_list(101) = 'pow_i'
  expro_list(102) = 'pow_ei'
  expro_list(103) = 'pow_e_aux'
  expro_list(104) = 'pow_i_aux'
  expro_list(105) = 'pow_e_fus'
  expro_list(106) = 'pow_i_fus'
  expro_list(107) = 'pow_e_sync'
  expro_list(108) = 'pow_e_brem'
  expro_list(109) = 'pow_e_line'
  expro_list(110) = 'pow_e_ohmic'
  expro_list(111) = 'polflux'
  expro_list(112) = 'shot'
  expro_list(113) = 'time'

  ! new scalars
  
  expro_list(114) = 'betap'
  expro_list(115) = 'betat'
  expro_list(116) = 'betan'
  expro_list(117) = 'greenwald'
  expro_list(118) = 'ptransp'
  expro_list(119) = 'tau'
  expro_list(120) = 'tau98y2'
  
end subroutine expro_list_set

end module expro
