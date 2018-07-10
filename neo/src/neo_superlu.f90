subroutine SOLVE_superlu(m_size, m_iindx, m_jindx, m, m_asize, m_rows, rhs, sol)
    use neo_globals
    implicit none
    integer, intent (in) :: m_size
    integer, dimension(*), intent (in) :: m_iindx
    integer, dimension(*), intent (in) :: m_jindx
    real, dimension(*), intent (in) :: m
    integer, intent (in) :: m_asize
    integer, intent (in) :: m_rows
    real, dimension(*), intent (in) :: rhs
    real, dimension(*), intent (inout) :: sol

    real, dimension(:), allocatable :: acc
    real, dimension(:), allocatable :: b
    integer, dimension(:), allocatable :: rowindx
    integer, dimension(:), allocatable :: colindx

    integer :: iopt, n, nnz, nrhs=1
    integer(KIND=8):: factors

     !Get the number of non-zero elements in compressed column format
     call st_to_cc_size ( m_size, m_iindx, m_jindx, nnz )

     !Allocate arrays to hold compressed column format
     n = m_rows
     allocate(acc(1:nnz))
     allocate(rowindx(1:nnz))
     allocate(colindx(1:n+1))
     allocate(b(1:n))
     b(1:n) = rhs(1:n)

     !Convert sparse triplet to compressed column format 
     call st_to_cc_index ( m_size, m_iindx, m_jindx, nnz, n, rowindx, colindx )
     call st_to_cc_values ( m_size, m_iindx, m_jindx, m, nnz, n, rowindx, colindx, acc )

     ! Factor the Matrix -- uses acc(:) and [row/col]indx(:)
     if(silent_flag == 0 .and. i_proc == 0) then
        open(unit=io_neoout,file=trim(path)//runfile_neoout,&
             status='old',position='append')
        write(io_neoout,'(t2,a)') 'Begin matrix factor -- SuperLU'
        close(io_neoout)
     endif

     !First, factorize the matrix. The factors are stored in *factors* handle.
     iopt = 1
     call c_fortran_dgssv( iopt, n, nnz, nrhs, acc, rowindx, colindx, b, n, factors, error_status )
     if (error_status .ne. 0) return

     if(silent_flag == 0 .and. i_proc == 0) then
        open(unit=io_neoout,file=trim(path)//runfile_neoout,&
             status='old',position='append')
        write(io_neoout,'(t2,a)') 'Done matrix factor -- SuperLU'
        close(io_neoout)
     endif

     ! Matrix solve -- uses rhs(:), acc(:), and [row/col]indx(:)
     if(silent_flag == 0 .and. i_proc == 0) then
        open(unit=io_neoout,file=trim(path)//runfile_neoout,&
             status='old',position='append')
        write(io_neoout,'(t2,a)') 'Begin matrix solve -- SuperLU'
        close(io_neoout)
     endif

     !Second, solve the system using the existing factors.
     iopt = 2
     call c_fortran_dgssv( iopt, n, nnz, nrhs, acc, rowindx, colindx, b, n, factors, error_status )
     if (error_status .ne. 0) return

     if(silent_flag == 0 .and. i_proc == 0) then
        open(unit=io_neoout,file=trim(path)//runfile_neoout,&
             status='old',position='append')
        write(io_neoout,'(t2,a)') 'Done matrix solve -- SuperLU'
        close(io_neoout)
     endif

     sol(1:n) = b(1:n)

     !Last, free the storage allocated inside SuperLU
     iopt = 3
     call c_fortran_dgssv( iopt, n, nnz, nrhs, acc, rowindx, colindx, b, n, factors, error_status )
     deallocate(acc)
     deallocate(b)
     deallocate(rowindx)
     deallocate(colindx)
end subroutine SOLVE_superlu
