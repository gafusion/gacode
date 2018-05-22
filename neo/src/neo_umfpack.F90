module neo_umfpack

  implicit none

  public :: SOLVE_umfpack, SOLVE_factor, SOLVE_do

  real,    dimension(10), private ::  cntl
  integer, dimension(20), private ::  icntl
  real,    dimension(20), private ::  rinfo
  integer, dimension(40), private ::  uinfo
  integer, dimension(20), private ::  keep
  character (len=40), private :: mystr
  
contains

  subroutine SOLVE_umfpack(m_size, m_iindx, m_jindx, m, m_asize, m_rows, rhs, sol)
    use neo_globals
    implicit none
    integer, intent (in) :: m_size
    integer, dimension(:), intent (in) :: m_iindx
    integer, dimension(:), intent (in) :: m_jindx
    real, dimension(:), intent (in) :: m
    integer, intent (in) :: m_asize
    integer, intent (in) :: m_rows
    real, dimension(:), intent (in) :: rhs
    real, dimension(:), intent (inout) :: sol
    integer :: k, ierr, ifac, matfac_err, max_ifac=5

     n_max = m_size*matsz_scalefac
     matfac_err = 0

     do ifac = 1, max_ifac

        if(silent_flag == 0 .and. i_proc == 0) then
           open(unit=io_neoout,file=trim(path)//runfile_neoout,&
                status='old',position='append')
           write(io_neoout,'(t2,a,e12.5)') 'Estimated memory (GB) = ', 2*8.0*(m_asize+n_max)/1.0e9
           close(io_neoout)
        endif
        if(allocated(amat))       deallocate(amat)
        allocate(amat(n_max),stat=ierr)
        if(ierr /= 0) then
           call neo_error('ERROR: (NEO) Array allocation failed')
           error_status = 100
           return
        end if
        if(allocated(amat_indx))  deallocate(amat_indx)
        allocate(amat_indx(2*n_max),stat=ierr)
        if(ierr /= 0) then
           call neo_error('ERROR: (NEO) Array allocation failed')
           error_status = 100
           return
        end if

        amat(:) = 0.0
        amat_indx(:) = 0
        do k=1,m_size
           amat(k) = m(k)
           amat_indx(k) = m_iindx(k)
           amat_indx(m_size+k) = m_jindx(k)
        enddo

        if(silent_flag == 0 .and. i_proc == 0) then
           open(unit=io_neoout,file=trim(path)//runfile_neoout,&
                status='old',position='append')
           if(ifac == 1) then
              write(io_neoout,'(t2,a)') 'Begin matrix factor -- UMFPACK'
           else
              write(io_neoout,'(t2,a)') 'Re-trying matrix factorization -- UMFPACK'
           endif
           close(io_neoout)
        endif
        call SOLVE_factor(m_size)
        if(error_status > 0) then
           error_status = 0
           n_max = n_max * 2
        else
           matfac_err = 1
           exit
        endif
     enddo
     if(matfac_err == 0) then
        call neo_error('ERROR: (NEO) Matrix factorization failed. Try increasing MATSZ_SCALEFAC.')
        error_status = 100
        return
     endif
     if(silent_flag == 0 .and. i_proc == 0) then
        open(unit=io_neoout,file=trim(path)//runfile_neoout,&
             status='old',position='append')
        write(io_neoout,'(t2,a)') 'Done matrix factor -- UMFPACK'
        close(io_neoout)
     endif

     ! Matrix solve -- uses g(:), a(:), and a_indx(:)
     if(silent_flag == 0 .and. i_proc == 0) then
        open(unit=io_neoout,file=trim(path)//runfile_neoout,&
             status='old',position='append')
        write(io_neoout,'(t2,a)') 'Begin matrix solve -- UMFPACK'
        close(io_neoout)
     endif
     call SOLVE_do
     if(silent_flag == 0 .and. i_proc == 0) then
        open(unit=io_neoout,file=trim(path)//runfile_neoout,&
             status='old',position='append')
        write(io_neoout,'(t2,a)') 'Done matrix solve -- UMFPACK'
        close(io_neoout)
     endif
  end subroutine SOLVE_umfpack

  subroutine SOLVE_factor(n_elem)
    use neo_globals
    implicit none
    integer, intent (in) :: n_elem
    
    call UMD21I(keep,cntl,icntl)
 
    !--------------------------------------------
    ! Error and diagnostic messages from UMFPACK,
    ! to unit=5.
    !
    icntl(1) = 6
    icntl(2) = 6
    !
    ! Use icntl(3) = 2 for more verbosity.
    icntl(3) = 0
    !
    ! Use icntl(8)=n for n levels of iterative 
    ! refinement (set job=1).
    icntl(8) = 0
    !--------------------------------------------

    ! Factor Maxwell matrices
    call UMD2FA(n_row,&
         n_elem,&
         0,&
         .false.,&
         n_max,&
         2*n_max,&
         amat,&
         amat_indx,&
         keep,&
         cntl,&
         icntl,&
         uinfo,&
         rinfo)

    if(uinfo(1) < 0) then
       call neo_error('ERROR: (NEO) Matrix factorization failed in neo_sparse_solve')
       write(mystr,'(A6,I5)') "uinfo=", uinfo(1)
       call neo_error(mystr)
       return
    endif


  end subroutine SOLVE_factor

  subroutine SOLVE_do
    use neo_globals
    implicit none
    real, dimension(:), allocatable :: w, gs

    allocate(gs(n_row))
    allocate(w(2*n_row))

    call UMD2SO(n_row,&
         0,&
         .false.,&
         n_max,&
         2*n_max,&
         amat,&
         amat_indx,&
         keep,&
         g,&
         gs,&
         w,&
         cntl,&
         icntl,&
         uinfo,&
         rinfo)
    
    g = gs

    deallocate(gs)
    deallocate(w)

  end subroutine SOLVE_do

end module neo_umfpack
