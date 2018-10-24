!-----------------------------------------------------------
! tgyro_restart.f90
!
! PURPOSE:
!  Manage restart.
!----------------------------------------------------------

subroutine tgyro_restart

  use mpi
  use tgyro_globals

  implicit none

  integer :: i,ip
  integer :: j
  integer :: ioerr
  integer :: p
  integer :: is,i_ion

  character(len=1) :: dummy
  character(len=2) :: itag
  real, dimension(6) :: x_read
  real, dimension(2:n_r,loc_n_ion+4) :: res2,relax2


  if (i_proc_global == 0) then

     open(unit=1,&
          file='out.tgyro.control',&
          status='old',&
          iostat=ioerr)

     if (ioerr /= 0) then

        close(1)
        loc_restart_flag = 0
        stop

     else

        read(1,*) n_r
        read(1,*) n_evolve
        read(1,*) i_tran
        close(1)

        ! 1. Gradients

        open(unit=1,file='out.tgyro.profile',status='old')
        do j=0,i_tran
           read(1,'(a)') dummy
           read(1,'(a)') dummy
           do i=1,n_r
              read(1,*) x_read(1:6)
              f_rot(i) = x_read(6)/r_min
           enddo
        enddo
        close(1)

        open(unit=1,file='out.tgyro.profile_e',status='old')
        do j=0,i_tran
           read(1,'(a)') dummy
           read(1,'(a)') dummy
           do i=1,n_r
              read(1,*) x_read(1:6)
              dlnnedr(i) = x_read(3)/r_min
              dlntedr(i) = x_read(5)/r_min
           enddo
        enddo
        close(1)

        do i_ion=1,loc_n_ion
           itag = 'i'//trim(ion_tag(i_ion))     
           open(unit=1,file='out.tgyro.profile_'//itag,status='old')
           do j=0,i_tran
              read(1,'(a)') dummy
              read(1,'(a)') dummy
              do i=1,n_r
                 read(1,*) x_read(1:6)
                 dlnnidr(i_ion,i) = x_read(3)/r_min
                 dlntidr(i_ion,i) = x_read(5)/r_min
              enddo
           enddo
           close(1)
        enddo

        ! 2. Total fluxes

        open(unit=1,file='out.tgyro.evo_er',status='old')
        do j=0,i_tran
           read(1,'(a)') dummy
           read(1,'(a)') dummy
           do i=1,n_r
              read(1,*) x_read(1:3)
              mflux_tot(i) = x_read(2)
           enddo
        enddo
        close(1)
        open(unit=1,file='out.tgyro.evo_ne',status='old')
        do j=0,i_tran
           read(1,'(a)') dummy
           read(1,'(a)') dummy
           do i=1,n_r
              read(1,*) x_read(1:3)
              pflux_e_tot(i) = x_read(2)
           enddo
        enddo
        close(1)
        open(unit=1,file='out.tgyro.evo_te',status='old')
        do j=0,i_tran
           read(1,'(a)') dummy
           read(1,'(a)') dummy
           do i=1,n_r
              read(1,*) x_read(1:3)
              eflux_e_tot(i) = x_read(2)
           enddo
        enddo
        close(1)
        open(unit=1,file='out.tgyro.evo_ti',status='old')
        do j=0,i_tran
           read(1,'(a)') dummy
           read(1,'(a)') dummy
           do i=1,n_r
              read(1,*) x_read(1:3)
              eflux_i_tot(i) = x_read(2)
           enddo
        enddo
        close(1)
        do i_ion=1,loc_n_ion
           open(unit=1,file='out.tgyro.evo_n'//trim(ion_tag(i_ion)),status='old')
           do j=0,i_tran
              read(1,'(a)') dummy
              read(1,'(a)') dummy
              do i=1,n_r
                 read(1,*) x_read(1:3)
                 pflux_i_tot(i_ion,i) = x_read(2)
              enddo
           enddo
           close(1)
        enddo

        ! 3. Residuals

        open(unit=1,file='out.tgyro.residual',status='old')
        do j=0,i_tran 
           read(1,'(a)') dummy
           read(1,'(a)') dummy
           p = 0
           read(1,'(a)') dummy
           do i=2,n_r
              read(1,*) &
                   x_read(1),(res2(i,ip),relax2(i,ip),ip=1,4+loc_n_ion)
              if (loc_ti_feedback_flag == 1) then
                 p = p+1
                 res(p) = res2(i,1) 
                 relax(p) = relax2(i,1) 
              endif
              if (loc_te_feedback_flag == 1) then
                 p  = p+1
                 res(p) = res2(i,2) 
                 relax(p) = relax2(i,2)
              endif
              if (loc_er_feedback_flag == 1) then
                 p  = p+1
                 res(p) = res2(i,3) 
                 relax(p) = relax2(i,3)
              endif
              do is=0,loc_n_ion
                 if (evo_e(is) >= 1) then
                    p  = p+1
                    res(p) = res2(i,4+is)
                    relax(p) = relax2(i,4+is)
                 endif
              enddo
           enddo
        enddo
        close(1)

     endif

  endif

  ! Reset n_evolve in case it was changed from value stored in control.out.

  n_evolve = loc_ti_feedback_flag+loc_te_feedback_flag+loc_er_feedback_flag
  do is=0,loc_n_ion
     if (evo_e(is) >= 1) n_evolve=n_evolve+1
  enddo

  call MPI_BCAST(loc_restart_flag,&
       1,&
       MPI_INTEGER,&
       0,&
       MPI_COMM_WORLD,&
       ierr)

  call MPI_BCAST(i_tran,&
       1,&
       MPI_INTEGER,&
       0,&
       MPI_COMM_WORLD,&
       ierr)

  if (loc_restart_flag == 1) then

     call MPI_BCAST(dlnnidr,&
          size(dlnnidr),&
          MPI_DOUBLE_PRECISION,&
          0,&
          MPI_COMM_WORLD,&
          ierr)

     call MPI_BCAST(dlnnedr,&
          size(dlnnedr),&
          MPI_DOUBLE_PRECISION,&
          0,&
          MPI_COMM_WORLD,&
          ierr)

     call MPI_BCAST(dlntidr,&
          size(dlntidr),&
          MPI_DOUBLE_PRECISION,&
          0,&
          MPI_COMM_WORLD,&
          ierr)

     call MPI_BCAST(dlntedr,&
          size(dlntedr),&
          MPI_DOUBLE_PRECISION,&
          0,&
          MPI_COMM_WORLD,&
          ierr)

     call MPI_BCAST(f_rot,&
          size(f_rot),&
          MPI_DOUBLE_PRECISION,&
          0,&
          MPI_COMM_WORLD,&
          ierr)

     call MPI_BCAST(eflux_i_tot,&
          size(eflux_i_tot),&
          MPI_DOUBLE_PRECISION,&
          0,&
          MPI_COMM_WORLD,&
          ierr)

     call MPI_BCAST(eflux_e_tot,&
          size(eflux_e_tot),&
          MPI_DOUBLE_PRECISION,&
          0,&
          MPI_COMM_WORLD,&
          ierr)

     call MPI_BCAST(pflux_e_tot,&
          size(pflux_e_tot),&
          MPI_DOUBLE_PRECISION,&
          0,&
          MPI_COMM_WORLD,&
          ierr)

     call MPI_BCAST(pflux_i_tot,&
          size(pflux_i_tot),&
          MPI_DOUBLE_PRECISION,&
          0,&
          MPI_COMM_WORLD,&
          ierr)

     call MPI_BCAST(mflux_tot,&
          size(mflux_tot),&
          MPI_DOUBLE_PRECISION,&
          0,&
          MPI_COMM_WORLD,&
          ierr)

     call MPI_BCAST(relax,&
          size(relax),&
          MPI_DOUBLE_PRECISION,&
          0,&
          MPI_COMM_WORLD,&
          ierr)

     call MPI_BCAST(res,&
          size(res),&
          MPI_DOUBLE_PRECISION,&
          0,&
          MPI_COMM_WORLD,&
          ierr)

  endif

end subroutine tgyro_restart
