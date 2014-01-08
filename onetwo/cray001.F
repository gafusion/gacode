
      block data initialize
c
      USE param
      USE adaptive
      USE constnts
      USE mcgo
      USE staebler
      implicit  integer (i-n), real*8 (a-h, o-z)
c

      include 'quit.i'
c      include 'staebler.i'


      data time_step_too_small /.false./    ! quit.i

      end




#ifdef PARALLEL_12
      subroutine onetwo
#else
      program onetwo
#endif
      USE ext_prog_info
      USE io,   ONLY   : versid
c
      implicit none
c
      integer*4 mpiIerror,myid,ierr,omp_get_thread_num
      character rcs_id*63 
      save      rcs_id
      data      rcs_id /
     ."$Id: cray001.F,v 1.20 2013/12/12 18:37:01 stjohn Exp $"/


      include 'quit.i'
D      include 'mpif.h'
D     call MPI_Init(mpiIerror)
D      call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr) !get processor id
D      print *,'proc. id initalization',myid 
D      call MPI_Barrier( MPI_COMM_WORLD,ierr)


!------------------------------------------------------------
! -- NOTE: VERSID IS DECODED SO DO NOT CHANGEGE THE FORMAT:
! --       'Onetwo VX.XX' (or Vx.x)
! -- wih x.xx some number
!------------------------------------------------------------------
!      versid = 'Onetwo V3.2 '  ! 2/13/04  HSJ
!      versid = 'Onetwo V3.4 '  ! 6/21/04  HSJ(introduces vloop_bc)
!      versid = 'Onetwo V3.5 '  ! 03/04/05 HSJ(paleoclassical and nubeam
!                                               restart refinement)
!      versid = 'Onetwo V3.6 '  ! changed nubeam restart logic
!      versid = 'Onetwo V3.93'  ! 
!      versid = 'Onetwo V3.94'  ! add some nubplt output and new toray interface
!      versid = 'Onetwo V4.10'  ! merge with JMP 
!      versid = 'Onetwo V4.10'  ! merge with JMP 
!      versid = 'Onetwo V4.11'  ! merge with JMP 
!      versid = 'Onetwo V4.20'  ! parallel version 
!      versid = 'Onetwo V4.21'  ! merged with JMP ,include STAR1-STAR6, HSJ 01/08/06 
!      versid = 'Onetwo V4.22'  ! fixed h98y2 and H_Petty to include
                                ! thermal energy only
!      versid = 'Onetwo V4.23'  ! modify paleoclassical transport formulation
                                ! paleocl.f90 to account for new development.2/13/07  
!       versid = 'Onetwo V4.24'  ! include option to use old version 1.29 of iterdb
                                ! file, see cray102.f for info.

!       versid = 'Onetwo V4.25'  ! Correction to paloecl.f90
!       versid = 'Onetwo V4.26'   ! update ext_prog_info.f90 to point to latest preplt for 51 frid point version (201 was ok)
       versid = 'Onetwo V5.0'   ! this version adds GCNMP to Onetwo solver
       versid = 'Onetwo V5.1'   ! this version updates The nubeam interface slightly
       versid = 'Onetwo V5.2'   ! this version adds more glf23 diagnostics 
       versid = 'Onetwo V5.3'   ! added density_ult,te_mult,ti_mult,etc.,new
                                ! statefile output
       versid = 'Onetwo V5.43'  ! started on 9/25/10  continuation of transport 
                                ! statefile development
       versid = 'Onetwo V5.44'  ! started on 9/25/10  continuation of transport 
                                ! statefile development

       versid = 'Onetwo V5.45'  ! update for ech and neg IP in tdem mode
       versid = 'Onetwo V5.46'  ! update for venus,star,P_Nfreya
       versid = 'Onetwo V5.47'  ! add timing routines for perfomance study 5/7/13  
       versid ='Onetwo  V5.48'  ! removed dbg calls
       versid ='Onetwo  V5.49'  ! added he3 fusion calcs 8/29/13 HSJ
       versid ='Onetwo  V5.50'  ! added VIDAR as running host 11/4/13 working yet in this verson 12/12/13   HSJ
                                ! robin bc not 
!$OMP parallel
!$      print *,"OMP  thread number =",omp_get_thread_num()
!$OMP end parallel


      call driver
      call DESTROY ('INONE')
D      print *,'waiting for processes to finalize'
D      call MPI_Barrier( MPI_COMM_WORLD,ierr)
D      print *,'done processes to finalize'
D      call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr) !get processor id
D      print *,'proc. id calling finalize',myid
D     call MPI_Finalize(mpiIerror)
        if(time_step_too_small)then
#ifdef PARALLEL_12
	     print *,'program ONETWO: ABNORMAL  termination' !jmp.ibm.par
#else
         call STOP ('program ONETWO: ABNORMAL  termination', 0)
#endif
      else
#ifdef PARALLEL_12
         print *,'program ONETWO: normal termination'
#else
         call STOP ('program ONETWO: normal termination', 0)
#endif
      endif
c
      end
