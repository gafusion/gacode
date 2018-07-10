subroutine SOLVE_petsc(m_size, m_iindx, m_jindx, m, m_asize, m_rows, rhs, sol)
    use neo_globals
#include <petsc/finclude/petscksp.h>
    use petscksp
    implicit none

    integer, intent (in) :: m_size
    integer, dimension(*), intent (in) :: m_iindx
    integer, dimension(*), intent (in) :: m_jindx
    real, dimension(*), intent (in) :: m
    integer, intent (in) :: m_asize
    integer, intent (in) :: m_rows
    real, dimension(*), intent (in) :: rhs
    real, dimension(*), intent (inout) :: sol

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!                   Variable declarations
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      PetscInt         n   ! (n x n matrix)
      Mat              A
      Vec              x,b,u
      KSP              ksp
      PC               pc
      PetscReal        norm,tol
      PetscErrorCode   ierr
      PetscInt         row(1), col(1)
      PetscScalar      el(1)
      PetscScalar      neg_one
      PetscInt         one
      PetscInt         i,k,its
      PetscBool        solve

      PetscInt,    ALLOCATABLE :: ix(:), nnz(:)
      PetscScalar, ALLOCATABLE :: s(:)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!                 Beginning of program
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
      if (ierr .ne. 0) then
        print*,'Unable to initialize PETSc'
        stop
      endif
      n       = m_rows
      one     = 1
      neg_one = -1.0
      solve = PETSC_FALSE

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Compute the matrix and right-hand-side vector that define
!      the linear system, Ax = b.
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      allocate(nnz(0:n-1))
      nnz(:) = 0
      do i=1,m_size
        k = m_iindx(i)-1
        nnz(k) = nnz(k) + 1
      enddo
      call MatCreateSeqAIJ(PETSC_COMM_SELF,n,n,PETSC_DEFAULT_INTEGER,nnz,A,ierr)
      call MatSetFromOptions(A,ierr)
      call MatSetUp(A,ierr)

      call MatZeroEntries(A,ierr)
      do i=1,m_size
        row(1) = m_iindx(i)-1
        col(1) = m_jindx(i)-1
        el(1) = m(i)
        call MatSetValues(A,one,row,one,col,el,ADD_VALUES,ierr)
      enddo

!  Assemble matrix, using the 2-step process:
      call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
      call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
      !call MatView(A,PETSC_VIEWER_STDOUT_SELF,ierr)

!  Set right-hand-side vector and create solution vectors space
      call VecCreateSeq(PETSC_COMM_SELF,n,x,ierr)
      call VecSetFromOptions(x,ierr)
      call VecDuplicate(x,b,ierr)
      call VecDuplicate(x,u,ierr)

      allocate(ix(0:n-1))
      do i=0,n-1
          ix(i) = i
      enddo
      call VecSetValues(b,n,ix,rhs,INSERT_VALUES,ierr)

      call VecAssemblyBegin(b, ierr)
      call VecAssemblyEnd(b, ierr)
      !call VecView(b,PETSC_VIEWER_STDOUT_SELF,ierr)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!         Create the linear solver and set various options
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      call KSPCreate(PETSC_COMM_SELF,ksp,ierr)

!  Set operators. Here the matrix that defines the linear system
!  also serves as the preconditioning matrix.
      call KSPSetOperators(ksp,A,A,ierr)


!  Set solver context and factorize
     if(silent_flag == 0 .and. i_proc == 0) then
        open(unit=io_neoout,file=trim(path)//runfile_neoout,&
             status='old',position='append')
        write(io_neoout,'(t2,a)') 'Begin matrix factor -- PETSc'
     endif

    if (use_slu == 1) then
#ifdef PETSC_HAVE_SUPERLU
     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,'(t2,a)') 'PETSc matrix solver -- SuperLU'
     endif
      call KSPSetType(ksp,KSPPREONLY,ierr)
      call KSPGetPC(ksp,pc,ierr)
      call PCSetType(pc,PCLU,ierr)
      call PCFactorSetMatSolverType(pc,MATSOLVERSUPERLU,ierr)
      call PCFactorSetUpMatSolverType(pc,ierr)
      solve = PETSC_TRUE
#endif
    endif

    if (solve .eqv. PETSC_FALSE) then
#if defined (PETSC_HAVE_UMFPACK) || defined (PETSC_HAVE_SUITESPARSE)
     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,'(t2,a)') 'PETSc matrix solver -- UMFPACK'
     endif
      call KSPSetType(ksp,KSPPREONLY,ierr)
      call KSPGetPC(ksp,pc,ierr)
      call PCSetType(pc,PCLU,ierr)
      call PCFactorSetMatSolverType(pc,MATSOLVERUMFPACK,ierr)
      call PCFactorSetUpMatSolverType(pc,ierr)
      solve = PETSC_TRUE
#endif
    endif

    if (solve .eqv. PETSC_FALSE) then
#ifdef PETSC_HAVE_MUMPS
     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,'(t2,a)') 'PETSc matrix solver -- MUMPS'
     endif
      call KSPSetType(ksp,KSPPREONLY,ierr)
      call KSPGetPC(ksp,pc,ierr)
      call PCSetType(pc,PCLU,ierr)
      call PCFactorSetMatSolverType(pc,MATSOLVERMUMPS,ierr)
      call PCFactorSetUpMatSolverType(pc,ierr)
      solve = PETSC_TRUE
#endif
    endif

    if (solve .eqv. PETSC_FALSE) then
      if(silent_flag == 0 .and. i_proc == 0) then
        close(io_neoout)
      endif
      call neo_error('ERROR: (NEO) PETSc not configured to use a sparse solver')
      return
    endif

      call KSPSetFromOptions(ksp,ierr)
      call KSPSetUp(ksp,ierr)

     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,'(t2,a)') 'Done matrix factor -- PETSc'
        close(io_neoout)
     endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!                      Solve the linear system
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     if(silent_flag == 0 .and. i_proc == 0) then
        open(unit=io_neoout,file=trim(path)//runfile_neoout,&
             status='old',position='append')
        write(io_neoout,'(t2,a)') 'Begin matrix solve -- PETSc'
        close(io_neoout)
     endif

      call KSPSolve(ksp,b,x,ierr)

      !View solver info; we could instead use the option -ksp_view
      !call KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD,ierr)
      !call VecView(x,PETSC_VIEWER_STDOUT_SELF,ierr)

     if(silent_flag == 0 .and. i_proc == 0) then
        open(unit=io_neoout,file=trim(path)//runfile_neoout,&
             status='old',position='append')
        write(io_neoout,'(t2,a)') 'Done matrix solve -- PETSc'
        close(io_neoout)
     endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!                     Verify solution
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      call MatMult(A,x,u,ierr)
!      call VecAXPY(u,neg_one,b,ierr)
!      call VecNorm(u,NORM_2,norm,ierr)
!      if (norm .gt. 1.e-10) then
!        call neo_error('ERROR: (NEO) PETSc verification failed')
!        return
!      endif

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!                     Copy solution and clean up
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      allocate(s(0:n-1))
      s(:) = 0.0
      call VecGetValues(x,n,ix,s,ierr)
      sol(1:n) = s(0:n-1)

!  Free work space.  All PETSc objects should be destroyed when they
!  are no longer needed.
      call VecDestroy(x,ierr)
      call VecDestroy(b,ierr)
      call VecDestroy(u,ierr)
      call MatDestroy(A,ierr)
      call KSPDestroy(ksp,ierr)
      if(allocated(s))     deallocate(s)
      if(allocated(ix))    deallocate(ix)
      if(allocated(nnz))   deallocate(nnz)

!  Always call PetscFinalize() before exiting a program.
      call PetscFinalize(ierr)
end subroutine SOLVE_petsc
