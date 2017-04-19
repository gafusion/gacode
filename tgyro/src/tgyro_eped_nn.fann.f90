 module toq_profiles_mod

   implicit none
   private
   public :: toq_profiles

 contains

   subroutine toq_profiles( &
        psin, npsi, widthp, &
        nped13, tpedEV, &
        ncore13, tcoreEV, &
        nedge13, tedgeEV, &
        nexpin, nexpout, texpin, texpout, &
        p_0, n13)

     implicit none
     real, dimension (:) :: psin, p_0, n13
     integer :: npsi, i
     real :: widthp, nped13, tpedEV, ncore13, tcoreEV, nedge13, tedgeEV
     real :: nexpin, nexpout, texpin, texpout
     real :: xphalf, pconst, a_n, a_t, xped, ncoretanh, tcoretanh
     real :: xpsi, nval, nvalp, tval, tvalp, xtoped

     ! From model 127 of psetup.f in TOQ code

     xphalf=1.-widthp
     pconst=1.-tanh((1.0-xphalf)/widthp)
     a_n=2.*(nped13-nedge13)/(1.+tanh(1.)-pconst)
     a_t=2.*(tpedEV-tedgeEV)/(1.+tanh(1.)-pconst)
     xped=xphalf-widthp
     ncoretanh=0.5*a_n*(1.-tanh(-xphalf/widthp)-pconst)+nedge13
     tcoretanh=0.5*a_t*(1.-tanh(-xphalf/widthp)-pconst)+tedgeEV
     do i=1,npsi
        xpsi=psin(i)
        nval=0.5*a_n*(1.-tanh((xpsi-xphalf)/widthp)-pconst)+nedge13
        nvalp=-0.5*a_n/widthp*(1./cosh((xpsi-xphalf)/widthp)**2)
        tval=0.5*a_t*(1.-tanh((xpsi-xphalf)/widthp)-pconst)+tedgeEV
        tvalp=-0.5*a_t/widthp*(1./cosh((xpsi-xphalf)/widthp)**2)
        if (ncore13.gt.0. .and. xpsi.lt.xped) then
           xtoped=xpsi/xped
           nval=nval+(ncore13-ncoretanh)*(1.-xtoped**nexpin)**nexpout
           nvalp=nvalp-(ncore13-ncoretanh)*nexpin*nexpout*xtoped**(nexpin-1.)*&
                (1.-xtoped**nexpin)**(nexpout-1.)
        endif
        if (tcoreEV.gt.0. .and. xpsi.lt.xped) then
           xtoped=xpsi/xped
           tval=tval+(tcoreEV-tcoretanh)*(1.-xtoped**texpin)**texpout
           tvalp=tvalp-(tcoreEV-tcoretanh)*texpin*texpout*xtoped**(texpin-1.)*&
                (1.-xtoped**texpin)**(texpout-1.)
        endif
        p_0(i)=2.*(nval*1.e13)*(tval*1.6022e-12)
        n13(i)=nval
     enddo
   end subroutine toq_profiles

 end module toq_profiles_mod

 subroutine tgyro_eped_nn

   use toq_profiles_mod
   use tgyro_ped
   use tgyro_globals

   implicit none

   character(len=1000) :: nn_executable
   character(len=1000) :: nn_files
   character(len=1) :: dummy

   character(len=1000) :: epednn_model

   real(4) :: INPUT_PARAMETERS(10)
   real(4) :: OUTPUT_PARAMETERS(5)

   real :: nn_p_ped, nn_t_ped
   real :: nn_t_edg, nn_n_edg
   real :: nn_t_cor, nn_n_cor

   integer :: i

   real :: nexpin=1.1, nexpout=1.1, texpin=1.2, texpout=1.4

   CHARACTER NUL
   PARAMETER(NUL = CHAR(0))

   include 'brainfuse_lib.inc'

   INPUT_PARAMETERS( 1) = a_in
   INPUT_PARAMETERS( 2) = betan_in
   INPUT_PARAMETERS( 3) = bt_in
   INPUT_PARAMETERS( 4) = delta_in
   INPUT_PARAMETERS( 5) = ip_in
   INPUT_PARAMETERS( 6) = kappa_in
   INPUT_PARAMETERS( 7) = m_in
   INPUT_PARAMETERS( 8) = nped_in
   INPUT_PARAMETERS( 9) = r_in
   INPUT_PARAMETERS(10) = zeffped_in

   call get_environment_variable('EPEDNN_MODEL_DIR',epednn_model)

   ierr = load_anns(1, TRIM(epednn_model)//NUL,'brainfuse'//NUL)
   ierr = load_anns_inputs(INPUT_PARAMETERS)
   ierr = run_anns()
   ierr = get_anns_avg_array(OUTPUT_PARAMETERS)

   nn_w_ped = OUTPUT_PARAMETERS(3)

   ! nn_p* -> Pa
   ! nn_t* -> eV
   ! nn_n* -> 10^13/cm^3 

   nn_p_ped = OUTPUT_PARAMETERS(1)*1e6
   nn_t_ped = (10*nn_p_ped)/(2*(1e13*nped_in)*k)

   nn_n_cor = nped_in*1.5
   nn_t_cor = t_axis

   nn_n_edg = nped_in*0.25
   nn_t_edg = 75.0

   nn_vec=0.0
   do i=0,nx_nn-1
      nn_vec(i+1,1) = i/(nx_nn-1.0)
   enddo

   call toq_profiles( &
        nn_vec(:,1), nx_nn, nn_w_ped/2.0, &
        nped_in, nn_t_ped,                &
        nn_n_cor, nn_t_cor,               &
        nn_n_edg, nn_t_edg,               &
        nexpin, nexpout, texpin, texpout, &
        nn_vec(:,3), nn_vec(:,2))

   nn_vec(:,2) = nn_vec(:,2)*1e13
   nn_vec(:,3) = nn_vec(:,3)/10.0

   ! At this point: 
   !   nn_vec(:,1) -> psi_norm
   !   nn_vec(:,2) -> <n> [1/cm^3]
   !   nn_vec(:,3) -> P   [Pa]

 end subroutine tgyro_eped_nn
