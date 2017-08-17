!-----------------------------------------------------
! GKEIGEN_do.F
!
! PURPOSE:
!  Determine the linear time-derivative matrix 
!  time_derivative_matrix and find the gkeigen_n_values
!  eigenvectors and eigenvalues satisfying the condition
!  specified in GKEIGEN_METHOD.
!  Only works in linear mode.
!-----------------------------------------------------

      subroutine GKEIGEN_do

      use gyro_globals
      use gyro_pointers
      use GKEIGEN_globals

      implicit none

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#include "finclude/petsc.h"
!#include "finclude/petscvec.h"
!#include "finclude/petscmat.h"
!#include "finclude/petscis.h"
!#include "finclude/petscsys.h"
!#include "finclude/petscviewer.h"
#include "finclude/slepc.h"
!#include "finclude/slepceps.h"

!--------------------------------------------------
! PETSc, SLEPc, and other local variable declarations. These variables
! are called in PETSc/SLEPc functions so are declared this way to
! conform to the specific PETSc/SLEPc build being used.
!
      Integer        n_converged,                                       &
     &               n_iter,                                            &
     &               irows(h_length_block),                             &
     &               irows_t(h_length_block_t),                         &
     &               irows_vec(0:h_length_block-1),                     &
     &               icols(h_width_loc),                                &
     &               icols_t(h_length_loc),                             &
     &               mrow,                                              &
     &               ncol

      Vec            dummy_vec
!     -------------------------------------------------------
!     declare dummy_vec2 to avoid using PETSC_NULL in Fortran
!     -------------------------------------------------------
      Vec            dummy_vec2
      PetscReal      eigr
      PetscReal      eigi

      PetscScalar    dummy_scalar,dummy_scalar2,                        &
     &               gkeigen_target,                                    &
     &               one_col(h_length_loc),                             &
     &               M_seq(seq_length),                                 &
     &               M_seq_t(seq_length_t),                             &
     &               ev_vec_loc(0:h_length_block-1)

      PetscReal      residual

      Mat                                                               &
     &               time_derivative_matrix

      PetscErrorCode ierr, ierr0

      EPS            epsolve

      ST             trans

      KSP            ksub

      PC             precon

      PetscViewer    vecfile
!
!-------------------------------------------------

      integer :: irestart
      integer :: istrln
      integer :: i_init
      complex, dimension(500,h_length_loc) :: restart_cols
      complex, dimension(:,:), allocatable :: M_loc
      complex, dimension(:), allocatable :: M_pass_1
      complex, dimension(:), allocatable :: M_pass_2
      complex, dimension(0:h_length_loc-1) :: ev_vec
      character (len=80) :: gkeigen_message

10    Format(I2, ':', I2, ':', I2, ':' I3)
20    Format('Mode: ', I3, '    Omega: ', es18.9, '    Gamma: ', es18.9,&
     &       '    Residulal: ', es18.9) 
30    Format(es20.9, es20.9, es20.9) 
40    Format(I2, ',', I2, ' owns rows', I5, '-',I5, ' and columns ', I5,&
     &       '-', I5, ' with h_length_loc = ', I5)
50    Format(I2, ',', I2, ' jstate_start = ', I5, ' Local jsets: ',     &
     &       I5, '-', I5)
60    Format('Subset 1: ', I5, ' Subset 2: ', I5, ' Send from: ', I5,   &
     &    ' Send to: ', I5)
70    Format('Set: ', I2, ',', ' Proc: ', I2, ',',' Loc: ', I4, ',',    &
     &       ' Glob: ',I4,',',' Val: ',g12.5, ' + ', g12.5, ' i')
80    Format(I5, ' eigenvalues found')


      l_print = ( (gkeigen_j_set == 0) .AND. (i_proc == 0) )

      allocate(M_loc(0:h_length_loc-1,0:h_width_loc-1))     
      allocate(M_pass_1(h_length_loc*h_width_loc))

      call GKEIGEN_catch_error

!-------------------------------------------------
! Initialize the SLEPc operating environment,
! placeholder vector dummy_vec, and matrix
! time_derivative_matrix.
!
      call SlepcInitialize(PETSC_NULL_CHARACTER,ierr)
!
      call VecCreate(MPI_COMM_WORLD,dummy_vec,ierr)
      call VecCreate(MPI_COMM_WORLD,dummy_vec2,ierr)
      call VecSetSizes(dummy_vec,PETSC_DECIDE,h_length,ierr)
      call VecSetSizes(dummy_vec2,PETSC_DECIDE,h_length,ierr)
      call VecSetFromOptions(dummy_vec,ierr)
      call VecSetFromOptions(dummy_vec2,ierr)
!
!
      call MatCreate(MPI_COMM_WORLD,time_derivative_matrix,ierr)
      call MatSetSizes(time_derivative_matrix,PETSC_DECIDE,             &
     &                 PETSC_DECIDE,h_length,h_length,ierr)
      call MatSetType(time_derivative_matrix,MATDENSE,ierr)

!-------------------------------------------------


!------------------------------------------------
! Get all elements of time_derivative_matrix column
! by column with get_RHS.  Grid points owned by each
! processor are grouped together in the total state
! vector that defines time_derivative_matrix. 
!

      call MatGetOwnershipRange(time_derivative_matrix,jstate_start,    &
     &                          jstate_end,ierr)
      call MatGetSize(time_derivative_matrix,mrow,ncol,ierr) 

!      Print *, j_proc_tot, 'PETSc owns ', jstate_start, jstate_end
!     &         'size: ', mrow, ' x ', ncol

      istate_start = gkeigen_j_set * h_width_loc
      istate_end   = istate_start + h_width_loc - 1

      jstate_start = 0
      Do i_proc_e = 0, i_proc - 1
         jstate_start = jstate_start + n_nek_1 / n_proc_1
         If ( i_proc_e < Mod(n_nek_1,n_proc_1) )                        &
     &                  jstate_start = jstate_start + 1
      EndDo
      jstate_start = jstate_start * n_kinetic * n_x * n_stack
      jstate_end = jstate_start + h_length_loc - 1

!      Do jelem = 1, h_length_loc
!        irows(jelem) = jstate_start + jelem - 1
!      EndDo

      call date_and_time(values=time_array)
      If (l_print) Then
       Print 10, time_array(5),time_array(6),time_array(7),time_array(8)
       Print *, 'Starting value assignments.'
      EndIf

!      If (eigensolve_restart_flag == 1) Then
!        call GKEIGEN_matrix_read(n_last_col)
!        If (l_print) print *, " Restarting from column ", n_last_col
!      Else
!        n_last_col = 0
!        If ( (l_print) .AND. (gkeigen_mwrite_flag==1) ) Then
!          open(unit=1,file=trim(path)//file_eigen_restart,              &
!     &         status='replace')
!          write (1,*) "!  Matrix elements for eigensolver."
!          close(unit=1)
!        EndIf
!      EndIf

!      do step = 1, nstep
!        call gyro_field_solve_explicit
!        call gyro_rhs_total
!        h = h + dt*RHS
!      enddo
!      jset_start = (i_proc*gkeigen_proc_mult+gkeigen_j_set) *           &
!     &                                             h_length_block
      jset_start = j_proc_tot*h_length_block
      jset_end = jset_start + h_length_block - 1
      do jelem = jset_start, jset_end
        jelem0 = jelem - jset_start
        irows_vec(jelem0) = jelem
      enddo
!      print *, j_proc_tot, i_proc, gkeigen_j_set, jset_start, jset_end, &
!     &         jstate_start, '*****************'
      jstate = jstate_start
      jelem0 = 0
      Do is2 = 1, n_kinetic
       p_nek_loc2 = 0
       Do p_nek2 = 1+i_proc, n_nek_1, n_proc_1
        p_nek_loc2 = p_nek_loc2 + 1
        Do i2 = 1, n_x
         Do m2 = 1, n_stack
           if ((jstate>=jset_start).AND.(jstate<=jset_end)) then
             ev_vec_loc(jelem0) = h(m2,i2,p_nek_loc2,is2)
             jelem0 = jelem0 + 1
           endif
           jstate = jstate + 1
         EndDo ! m2
        EndDo ! i2
       EndDo ! p_nek2
      EndDo ! is2
      call VecSetValues(dummy_vec,h_length_block,irows_vec,             &
     &                  ev_vec_loc,INSERT_VALUES,ierr)
      call VecAssemblyBegin(dummy_vec,ierr)
      call VecAssemblyEnd(dummy_vec,ierr)

!      call VecView(dummy_vec,PETSC_VIEWER_STDOUT_WORLD,ierr)

      restart_cols(:,:) = (0.0,0.0)
      istate = 0

      Do i_proc_e = 0, n_proc_1 - 1

       Do is1 = 1, n_kinetic

        p_nek_loc1 = 0
        Do p_nek1 = 1+i_proc_e, n_nek_1, n_proc_1

         If (i_proc_e == i_proc) p_nek_loc1 = p_nek_loc1 + 1
         Do i1 = 1, n_x

          Do m1 = 1, n_stack

           l_in_i_block = ((istate >= istate_start) .AND.               &
     &                   (istate <= istate_end))
!           If ((eigensolve_restart_flag == 1) .AND.                     &
!     &                 (istate < n_last_col)) Then
!
!                   one_col(1:h_length_loc) =                            &
!     &                       eigensolve_matrix(1:h_length_loc,istate+1)
!
           If (l_in_i_block) Then

             h(:,:,:,:) = (0.0,0.0)  
             If (i_proc_e == i_proc)                                    &
     &                h(m1,i1,p_nek_loc1,is1) = (1.0,0.0)
               call gyro_field_solve_explicit
               call gyro_rhs_total
               ielem = istate - istate_start
               jstate = jstate_start
               Do is2 = 1, n_kinetic
                p_nek_loc2 = 0
                Do p_nek2 = 1+i_proc, n_nek_1, n_proc_1
                 p_nek_loc2 = p_nek_loc2 + 1
                 Do i2 = 1, n_x
                  Do m2 = 1, n_stack
                    jelem = jstate-jstate_start
                    one_col(jelem+1) = RHS(m2,i2,p_nek_loc2,is2)
                    M_loc(jelem,ielem)= RHS(m2,i2,p_nek_loc2,is2)  
!                    M_loc(jelem,ielem)=(1.,0.)*jelem+(0.,1.)*ielem
!                    if (l_print) print *, M_loc(jelem,ielem)
                    jstate = jstate + 1
                  EndDo ! m2
                 EndDo ! i2
                EndDo ! p_nek2
               EndDo ! is2

             EndIf

!             If (gkeigen_matrixonly == 0)                               &
!     &        call MatSetValues(time_derivative_matrix,h_length_loc,    &
!     &                     irows,1,istate,one_col,INSERT_VALUES,ierr)
             istate = istate + 1
!             irestart = istate-n_last_col
!             If (irestart >= 1) restart_cols(irestart,:)=one_col(:)
!             If (irestart == 500) Then
!               If (gkeigen_mwrite_flag==1)                              &
!     &           call GKEIGEN_matrix_write(restart_cols,n_last_col)
!               n_last_col = istate
!               restart_cols(:,:) = (0.0,0.0)
!             EndIf

          EndDo ! m1

         EndDo ! i1

        EndDo ! p_nek1

       EndDo ! is1

      EndDo ! i_proc_e
!
!---------------------------------------------

!      If (gkeigen_mwrite_flag==1)                                       &
!     &       call GKEIGEN_matrix_write(restart_cols,n_last_col)

      If (gkeigen_matrixonly == 1) Then
        call date_and_time(values=time_array)
        If (l_print) Then
          Print 10, time_array(5),time_array(6),time_array(7),          &
     &              time_array(8)
          Print *, 'Finished matrix construction without solving.'
        EndIf
        Goto 1000
      EndIf

      call date_and_time(values=time_array)
      If (l_print) Then
       Print 10, time_array(5),time_array(6),time_array(7),time_array(8)
       Print *, 'Starting message passing phase.'
      EndIf


!--------------------------------------------------------------------
! Use MPI_ALLTOALL to pass matrix elements to the process that owns
! the matrix elements in PETSc. Then assign values in PETSc.
! Memory allocation and deallocation here is necessary to keep memory 
! in bounds for large cases.
!
      Do jelem = 0, h_length_loc-1
        Do ielem = 0, h_width_loc-1
           iseq = jelem*h_width_loc + ielem + 1
           M_pass_1(iseq) = M_loc(jelem,ielem)
        EndDo
      EndDo
!
      deallocate(M_loc)
      allocate(M_pass_2(h_length_loc*h_width_loc))
!
       call MPI_ALLTOALL(M_pass_1,                                      &
     &                   seq_length,                                    &
     &                   MPI_DOUBLE_COMPLEX,                            &
     &                   M_pass_2,                                      &
     &                   seq_length,                                    &
     &                   MPI_DOUBLE_COMPLEX,                            &
     &                   GYRO_COMM_UNIPROC,                             &
     &                   ierr)
!
      deallocate(M_pass_1)
!
! Assign the values in M_pass_2 to the appropriate
! PETSc matrix elements
!
      Do gkeigen_j_set_e0 = 0, gkeigen_proc_mult-1
        iseq_start = gkeigen_j_set_e0*h_width_loc*h_length_block + 1
        iseq_end = iseq_start - 1 + h_width_loc*h_length_block
        Do iseq =  iseq_start, iseq_end
          iseq0 = iseq - iseq_start + 1
          M_seq(iseq0) = M_pass_2(iseq)
        EndDo
        Do ielem = 0, h_width_loc-1
          icols(ielem+1) = gkeigen_j_set_e0*h_width_loc + ielem
        EndDo
        jelem_start = gkeigen_j_set*h_length_block
        jelem_end = jelem_start + h_length_block - 1
        Do jelem = jelem_start, jelem_end
          i_index = jelem - jelem_start + 1
          irows(i_index) = jstate_start + jelem
        EndDo
        call MatSetValues(time_derivative_matrix,h_length_block,irows,  &
     &       h_width_loc,icols,M_seq,INSERT_VALUES,ierr)

      EndDo  ! gkeigen_j_set_e0
!
      deallocate(M_pass_2)
!
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!
!-----------------------------------------------------------------------

      call date_and_time(values=time_array)
      If (l_print) Then
       Print 10, time_array(5),time_array(6),time_array(7),time_array(8)
       Print *, 'Starting matrix assembly.'
      EndIf
      call MatAssemblyBegin(time_derivative_matrix,MAT_FINAL_ASSEMBLY,  &
     &                      ierr)
      call MatAssemblyEnd(time_derivative_matrix,MAT_FINAL_ASSEMBLY,    &
     &                    ierr)

!--------------------------------------------------------------------
! This creates the EigenProbelmSolver (EPS) framework and solves for the 
! gkeigen_n_values eigenvalues and vectors with the largest real part using a 
! Krylov subspace of maximum dimension gkeigen_kspace_dim.

!      gkeigen_kspace_dim = Max ( gkeigen_n_values+1,gkeigen_kspace_dim )
!      gkeigen_kspace_dim = Min ( gkeigen_kspace_dim,Int(h_length)   )

!      call MatView(time_derivative_matrix,PETSC_VIEWER_STDOUT_WORLD,    &
!     &             ierr)

      call date_and_time(values=time_array)
      If (l_print) Then
       Print 10, time_array(5),time_array(6),time_array(7),time_array(8)
       Print *, 'Starting eigensolver.'
       Print *, "Using Krylov subspace of dimension ",gkeigen_kspace_dim
      EndIf

      call EPSCreate(MPI_COMM_WORLD,epsolve,ierr)
!      call EPSSetType(epsolve,EPSARPACK,ierr)
      call EPSSetOperators(epsolve,time_derivative_matrix,              &
     &                     PETSC_NULL_OBJECT,ierr)

!      call EPSSetDimensions(epsolve,gkeigen_n_values,                   &
!     &                      gkeigen_kspace_dim,h_length,ierr)
      call EPSSetDimensions(epsolve,gkeigen_n_values,                   &
     &                      gkeigen_kspace_dim,PETSC_DECIDE,ierr)
!      call EPSSetTrueResidual(epsolve,PETSC_TRUE,ierr)
!      call EPSSetBalance(epsolve,EPS_BALANCE_TWOSIDE,200,PETSC_DECIDE,  &
!     &                   ierr)

!      call EPSGetST(epsolve,trans,ierr)
!      call STSetType(trans,STSHIFT,ierr)

!      call STGetKSP(trans,ksub,ierr)
!      call KSPSetType(ksub,KSPGMRES,ierr)
      
!      call KSPGetPC(ksub,precon,ierr)
!      call PCSetType(precon,PCASM,ierr)

!      call STSetType(trans,STPRECOND,ierr)
!      call KSPSetType(ksub,KSPCG,ierr)
!      call PCSetType(precon,PCJACOBI,ierr)

!      call KSPSetPC(ksub,precon,ierr)
!      call STSetKSP(trans,ksub,ierr)
!      call EPSSetST(epsolve,trans,ierr)

!--------------------------------------------------------------
! Search for eigenvalues in the region of the complex plane
! specified by gkeigen_method. GA frequency sign convention
! is opposite to sign of the imagainary part returned by SLEPc.
! GA convention is used in descriptions.
!
! (1) Largest real part / Largest growth rate
! (2) Largest imaginary part / Largest frequency in ion direction
! (3) Smallest imaginary part / Largest frequency in electron direction
! (4) Smallest real part / Most stable modes
!
!
      gkeigen_target = (1.0,0.0)*gkeigen_gamma_target -                 &
     &                                 (0.0,1.0)*gkeigen_omega_target
      call EPSSetTarget(epsolve,gkeigen_target,ierr)

      select case (gkeigen_method)

      case (1)

        call send_line(' Seeking largest growth rates.')
        If (l_print) print *, ' Seeking largest growth rates.'
        call EPSSetWhichEigenpairs(epsolve,EPS_LARGEST_REAL,ierr)

      case (2)

        call send_line(' Seeking largest negative frequency.')
        If (l_print) print *, ' Seeking largest negative frequency.'
        call EPSSetWhichEigenpairs(epsolve,EPS_LARGEST_IMAGINARY,ierr)

      case (3)

        call send_line(' Seeking largest positive frequency.')
        If (l_print) print *, ' Seeking largest positive frequency.'
        call EPSSetWhichEigenpairs(epsolve,EPS_SMALLEST_IMAGINARY,ierr)

      case (4)

        call send_line(' Seeking smallest growth rates.')
        If (l_print) print *, ' Seeking smallest growth rates.'
        call EPSSetWhichEigenpairs(epsolve,EPS_SMALLEST_REAL,ierr)

      case (5)

        write(gkeigen_message,*)                                        &
     & ' Seeking eigenvalues with frequency near ', gkeigen_omega_target
        call send_line(gkeigen_message)
        If (l_print) print *, gkeigen_message
        call EPSSetWhichEigenpairs(epsolve,EPS_TARGET_IMAGINARY,ierr)

      case (6)

        write(gkeigen_message,*)                                        &
     & ' Seeking eigenvalues with growth rate near ',                   &
     &        gkeigen_gamma_target
        call send_line(gkeigen_message)
        If (l_print) print *, gkeigen_message
        call EPSSetWhichEigenpairs(epsolve,EPS_TARGET_REAL,ierr)

      case default

        call catch_error('INVALID: gkeigen_method')

      end select
!--------------------------------------------------------------------

!      call PetscViewerBinaryOpen(MPI_COMM_WORLD,'vector.out',           &
!     &                           FILE_MODE_READ,vecfile,ierr)
!      call VecLoad(vecfile,PETSC_NULL,dummy_vec,ierr)

      call EPSSetInitialSpace(epsolve,1,dummy_vec,ierr)

      call EPSSetTolerances(epsolve,gkeigen_tol,gkeigen_iter,ierr)
      call EPSSetFromOptions(epsolve,ierr)
      call EPSSolve(epsolve,ierr)

      call date_and_time(values=time_array)
      If (l_print) Then
       Print 10, time_array(5),time_array(6),time_array(7),time_array(8)
       Print *, 'Finished eigensolve.'
      EndIf
!-------------------------------------------------------------------

!----------------------------------------------------------------
! All converged solutions are returned. Eigenvalues are stored
! in eigensolve_lambda, vectors in eigensolve_vector.

      call EPSGetConverged(epsolve,n_converged,ierr)
      If (l_print) Print *, n_converged, ' converged solution(s)'
      step = 0
      Do iev = 0, n_converged-1
        call EPSGetEigenpair(epsolve,iev,eigr,eigi,                     &
     &          dummy_vec,dummy_vec2,ierr)
        ! ev_omega = -aimag(dummy_scalar)
        ! ev_gamma = real(dummy_scalar)
        ev_omega = -eigi
        ev_gamma = eigr


        ! Interface change in SLEPC
        ! call EPSComputeResidualNorm(epsolve,iev,residual,ierr)
        call  EPSComputeError(epsolve,iev,EPS_ERROR_RELATIVE,           &
     &                        residual,ierr)
        If (l_print) Then
          print 20, iev+1, ev_omega, ev_gamma, residual
          if (iev==0) then
           open(unit=1,file=trim(path)//file_eigen_freq,                &
     &          status='replace')
          else
           open(unit=1,file=trim(path)//file_eigen_freq,                &
     &          status='old',position='append')
          endif
          write (1,30) ev_omega, ev_gamma, residual
          close (unit=1)
          call gyro_write_precision(10,abs(ev_omega))
          step = step + 1
          call gyro_write_precision(10,abs(ev_gamma))
        EndIf
        call VecAssemblyBegin(dummy_vec,ierr)
        call VecAssemblyEnd(dummy_vec,ierr)
        jset_start = j_proc_tot*h_length_block
        jset_end = jset_start + h_length_block - 1
        do jelem = jset_start, jset_end
          jelem0 = jelem - jset_start
          irows_vec(jelem0) = jelem
        enddo
        call VecGetValues(dummy_vec,h_length_block,irows_vec,           &
     &                    ev_vec_loc,ierr)

!        call VecView(dummy_vec,PETSC_VIEWER_STDOUT_WORLD,ierr)
!        do j_proc_tot_e1 = 0, n_proc_tot-1
!          if (j_proc_tot==j_proc_tot_e1) print *, ev_vec_loc
!          call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!        enddo

        if (gkeigen_j_set==0) then
          do jelem = 0, h_length_block-1
            ev_vec(jelem) = ev_vec_loc(jelem)
          enddo
        endif

        do gkeigen_j_set_e1 = 1, gkeigen_proc_mult-1
          j_proc_tot_e1 = i_proc*gkeigen_proc_mult + gkeigen_j_set_e1
          j_proc_tot_0  = i_proc*gkeigen_proc_mult
          if (gkeigen_j_set==0) then
            call MPI_RECV(ev_vec_loc,                                   &
     &                    size(ev_vec_loc),                             &
     &                    MPIU_SCALAR,                                  &
     &                    j_proc_tot_e1,                                &
     &                    j_proc_tot_e1,                                &
     &                    MPI_COMM_WORLD,                               &
     &                    recv_status,                                  &
     &                    ierr)
            jelem_start = gkeigen_j_set_e1*h_length_block
            jelem_end = jelem_start + h_length_block - 1
            do jelem = jelem_start, jelem_end
              jelem0 = jelem-jelem_start
              ev_vec(jelem) = ev_vec_loc(jelem0)
            enddo
          endif

          if (gkeigen_j_set==gkeigen_j_set_e1)                          &
     &      call MPI_SEND(ev_vec_loc,                                   &
     &                    size(ev_vec_loc),                             &
     &                    MPIU_SCALAR,                                  &
     &                    j_proc_tot_0,                                 &
     &                    j_proc_tot_e1,                                &
     &                    MPI_COMM_WORLD,                               &
     &                    ierr)
           call MPI_BARRIER(MPI_COMM_WORLD,ierr)
         enddo

!        call PetscViewerBinaryOpen(MPI_COMM_WORLD,'vector.out',         &
!     &                             FILE_MODE_WRITE,vecfile,ierr)
!        call VecView(dummy_vec,vecfile,ierr)

!        call VecView(dummy_vec,PETSC_VIEWER_STDOUT_WORLD,ierr)
!        do i_proc_e1 = 0, n_proc-1
!          if ((i_proc==i_proc_e1).AND.(gkeigen_j_set==0)) Then
!            print *, "Process ", i_proc_e1
!            print *, ev_vec
!          endif
!          call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!        enddo

        ev_vec = (1.*h_length)**0.5 * ev_vec
        alltime_index = 0
        time_skip = 1
        plot_filter = 0.0

        if (gkeigen_j_set==0) then
          t_current = t_current + 1.0
          data_step = data_step + 1
          jelem = 0
          Do is2 = 1, n_kinetic
            p_nek_loc2 = 0
            Do p_nek2 = 1+i_proc, n_nek_1, n_proc_1
              p_nek_loc2 = p_nek_loc2 + 1
              Do i2 = 1, n_x
                Do m2 = 1, n_stack
                  h(m2,i2,p_nek_loc2,is2)=ev_vec(jelem)
                  jelem = jelem + 1
                EndDo ! m2
              EndDo ! i2
            EndDo ! p_nek2
          EndDo ! is2

          call gyro_field_solve_explicit

          call gyro_field_interpolation
!          call gyro_field_fluxave
!          call gyro_field_time_derivative
          call gyro_field_plot
          call gyro_moments_plot
          call gyro_write_timedata

        endif

      EndDo

!---------------------------------------------------------------------

      If (l_print) Then
        call EPSGetIterationNumber(epsolve,n_iter,ierr)
        Print *, n_iter, ' iterations performed'
      EndIf

      write(gkeigen_message,80) n_converged
      call gyro_set_exit_status(gkeigen_message,2)

      call VecDestroy(dummy_vec,ierr)
      call VecDestroy(dummy_vec2,ierr)
      call MatDestroy(time_derivative_matrix,ierr)
      call EPSDestroy(epsolve,ierr)
!      call PetscViewerDestroy(vecfile,ierr)

      call SlepcFinalize(ierr)

1000  end subroutine

