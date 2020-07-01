program ptglf

  use mpi
  use ptglf_globals
  use tglf_interface

  implicit none

  integer :: p

  !---------------------------------------------------------------------
  ! Initialize MPI_COMM_WORLD communicator.
  !
  call MPI_INIT(i_err)
  call MPI_COMM_RANK(MPI_COMM_WORLD,i_proc,i_err)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,n_proc,i_err)
  !---------------------------------------------------------------------

  !---------------------------------------------------------------------
  ! Split MPI_COMM_WORLD into sub-communicators of size 1: loc_comm
  ! NOTE: this is trivial, and would be unnecessary if calling a
  !       serial version of TGLF
  ntot = n_proc
  splitkey = i_proc
  color = 0
  call MPI_COMM_SPLIT(MPI_COMM_WORLD,&
       color,&
       splitkey,&
       loc_comm,&
       i_err)
  !---------------------------------------------------------------------

  ! Path is cwd:
  path= './'

  !---------------------------------------------------------------------
  ! Allocate input and output data containers
  ! NOTE: will need to be expanded to cover all inputs/outputs
  allocate(indata_loc(2,ntot))
  allocate(outdata_loc(3,ntot))
  allocate(indata(2,ntot))
  allocate(outdata(3,ntot))
  !---------------------------------------------------------------------

  call tglf_init(path,loc_comm)
  call ptglf_init()

  if (i_proc == 0) print '(a,i0)','ntot = ',ntot

  !---------------------------------------------------------------------
  ! Set inputs depending in task rank
  p = 1+i_proc
  tglf_rlts_in(:) = 3.0+real(p)/ntot

  ! Pack inputs into array
  indata_loc(1,p) = tglf_rlts_in(1)
  indata_loc(2,p) = tglf_rlts_in(2)
  !---------------------------------------------------------------------

  !---------------------------------------------------------------------
  ! Run TGLF
  call tglf_run()

  ! Capture outputs
  outdata_loc(1,p) = tglf_elec_pflux_out
  outdata_loc(2,p) = tglf_elec_eflux_out
  outdata_loc(3,p) = tglf_elec_mflux_out
  !---------------------------------------------------------------------

  ! Collect all data 
  call MPI_ALLREDUCE(indata_loc,indata,size(indata), &
       MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,i_err)
  call MPI_ALLREDUCE(outdata_loc,outdata,size(outdata), &
       MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,i_err)

  if (i_proc == 0) then

     open(unit=1,file='out.ptglf.indata',status='replace')
     do p=1,ntot
        write(1,10) indata(:,p)
     enddo
     close(1)

     open(unit=1,file='out.ptglf.outdata',status='replace')
     do p=1,ntot
        write(1,10) outdata(:,p)
     enddo
     close(1)

  endif

  call MPI_finalize(i_err)

10 format(20(1pe17.10,1x))

end program ptglf
