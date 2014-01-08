  MODULE error_handler
  USE nrtype ,                         ONLY : DP,I4B,I2B
#ifdef GCNMP
  USE io_gcnmp,                        ONLY : ncrt 
#elif defined  NFREYA
  USE io_gcnmp,                        ONLY : ncrt 
#else
  USE  io,                             ONLY : ncrt
#endif
  USE MPI_data,                        ONLY : mpiierr,myid,master,initialized_mpi
#if defined (USEMPI)
      USE mpi
#endif
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: terminate, load_errno,iomaxerr,dbg_print
  INTEGER(I4B),PARAMETER   :: iomaxerr = 100
  INTEGER(I4B),PARAMETER   :: istart =1
  INTEGER(I4B),PARAMETER   :: mmaxerr = iomaxerr + 400
  INTEGER(I4B),PARAMETER   :: iend = iomaxerr + mmaxerr
  INTEGER(I4B),PUBLIC      :: lerrno, nfreya_offset,mhd_offset,nfreya_end
  INTEGER(I4B) imhd
  LOGICAL dbg_print


  TYPE errdat
    CHARACTER(len = 132) :: descrip
    CHARACTER(len = 132) :: attrib
    INTEGER(I2B)  errno
  END TYPE errdat
  TYPE (errdat) :: error(istart:iend)

 
  CONTAINS

    SUBROUTINE terminate(errno,nlog)
      INTEGER(I4B), INTENT(IN) ::  errno
      INTEGER(i2B),INTENT(IN)  ::  nlog
      
      IF(myid == master)THEN
        WRITE(ncrt,FMT='(a)')error(errno)%descrip(1:LEN_TRIM(error(errno)%descrip))
        IF(nlog .ne. ncrt) &
        WRITE(nlog,FMT='(a)')error(errno)%descrip(1:LEN_TRIM(error(errno)%descrip))
        IF(LEN_TRIM(error(errno)%attrib ) > 0)   THEN
           WRITE(ncrt,FMT='(a)')error(errno)%attrib(1:LEN_TRIM(error(errno)%attrib))
        IF(nlog .ne. ncrt) &
           WRITE(nlog,FMT='(a)')error(errno)%attrib(1:LEN_TRIM(error(errno)%attrib))
        ENDIF
        WRITE(ncrt,FMT='("error number =",i5)')errno
        IF(nlog .ne. ncrt) &
        WRITE(nlog,FMT='("error number =",i5)')errno
      ENDIF
#if defined (USEMPI)    
      CALL MPI_initialized(initialized_mpi,mpiierr)   ! if mpi is not initialized dont call MPI_ABORT
      IF(initialized_mpi ) &
           CALL MPI_ABORT(MPI_COMM_WORLD,lerrno,mpiierr)
#endif

      STOP
    END SUBROUTINE terminate



 
    SUBROUTINE load_errno
      IMPLICIT NONE
      INTEGER j
      DO j=istart,iend
         error(j)%descrip ='undefined'
         error(j)%attrib =''
         error(j)%errno  = j
      ENDDO





     error(1)%descrip = 'insufficient  command line arguments found'
     error(1)%attrib = 'usage: gcnmp  iterdb_file namelist_file run_log_file'
     error(1)%errno  = 1

     error(2)%descrip  = 'arg,carg, size mismatch in sub my_getarg '
     error(2)%attrib   = 'increase size of carg in program GCNM'
     error(2)%errno     = 2


     error(3)%descrip  = 'could not open iterdb input or output file'
     error(3)%errno     = 3


     error(4)%descrip  = 'could not open namelist  input file'
     error(4)%errno     = 4

     error(5)%descrip  = 'could not open output  file'
     error(5)%errno     = 5

     error(6)%descrip  = 'ERROR: max unit exceeded in get_next_io_unit'
     error(6)%errno     = 6

     error(7)%descrip  = 'ERROR: number of dependent variables does not agree '
     error(7)%errno     = 7

     error(8)%descrip  = 'ERROR: NAMELIST input: use_glf23 and use_glf23_flux not consistent'
     error(8)%errno     = 8

!     error(9)%descrip  = 'ERROR: NAMELIST input: use_glf23_flux and use_hyperbolic_den_eq not consistent'
     error(9)%descrip  = 'ERROR: NAMELIST input: use_mmm_flux and use_mmm not consistent'
     error(9)%errno     = 9


     error(10)%descrip  = 'ERROR: namelist lines are  longer than 255 characters'
     error(10)%errno     = 10

     error(11)%descrip  = 'ERROR: NAMELIST input: use_tglf(4) must be zero (eg Faradays law)'
     error(11)%errno     = 11

     error(12)%descrip  = 'ERROR: NAMELIST input: use_tglf(5) must be zero ( no toridal rotation in tglf at this time)'
     error(12)%errno     = 12

     error(13)%descrip  = 'ERROR: Iterdb file input: psir grid is not monotonic'
     error(13)%errno     = 13

     error(14)%descrip  = 'ERROR: NAMELIST  input: select_solver must be set explicitely for this case'
     error(14)%errno     = 14

     error(15)%descrip  = 'ERROR: NAMELIST  input: use_glf23(1) or use_tglf(1) or use_mmm(1) requies itenpd(1) is set'  
     error(15)%errno     = 15

     error(16)%descrip  = 'ERROR: NAMELIST  input: use_glf23(2) or use_tglf(2) or use_mmm(1) requies itte is set'
     error(16)%errno     = 16

     error(17)%descrip  = 'ERROR: NAMELIST  input: use_glf23(3) or use_tglf(3) or use_mmm(1) requies itti is set'
     error(17)%errno     = 17

     error(18)%descrip  = 'ERROR: NAMELIST  input: use_glf23(5) or use_tglf(5) or use_mmm(1) requies itw is set'
     error(18)%errno     = 18

     error(19)%descrip  = 'ERROR: NAMELIST  input: Maccormack solver requires that use_glf23_flux(1) is set '
     error(19)%errno     = 19

     error(20)%descrip  = 'ERROR: NAMELIST  input: Maccormack solver requires that use_glf23_flux(2) is set '
     error(20)%errno     = 20

     error(21)%descrip  = 'ERROR: NAMELIST  input: Maccormack solver requires that use_glf23_flux(3) is set '
     error(21)%errno     = 21

     error(22)%descrip  = 'ERROR: NAMELIST  input: Maccormack solver requires that use_glf23_flux(5) is set '
     error(22)%errno     = 22

     error(23)%descrip  = 'ERROR: NAMELIST  input: Only one turbulent confinement model may be on at a time'
     error(23)%errno     = 23

     error(24)%descrip  = 'ERROR: NAMELIST  input: Maccormac solver requires steady_state=1'
     error(24)%errno     = 24

     error(25)%descrip   =  'subroutine bc_conditions te error 1'
     error(25)%errno     = 25

     error(26)%descrip   =  'subroutine bc_conditions te error 2'
     error(26)%errno     = 26


     error(27)%descrip   =  'subroutine bc_conditions ti error 1'
     error(27)%errno     = 27


     error(28)%descrip   =  'subroutine bc_conditions ti error 2'
     error(28)%errno     = 28


     error(29)%descrip   =  'subroutine bc_conditions angrot  error 1'
     error(29)%errno     = 29


     error(30)%descrip   =  'subroutine bc_conditions angrot  error 2'
     error(30)%errno     =  30


     error(31)%descrip   =  'subroutine time_interp_ic, t_w out of range of  bctime'
     error(31)%errno     =  31



     error(32)%descrip   =  'subroutine RW_ITERDB_NETCDF, NETCDF ERROR'
     error(32)%errno     =  32

     error(33)%descrip   =  'subroutine RW_ITERDB_NETCDF, NETCDF FILE DIMENSION  ERROR'
     error(33)%errno     =  33

     error(34)%descrip   =  'subroutine RW_ITERDB_NETCDF, NETCDF CHARACTER DIMENSION  ERROR'
     error(34)%errno     =  34

     error(35)%descrip   =  'subroutine  SET_LABEL, ERROR in determination of indicator'
     error(35)%errno     =  35

     error(36)%descrip   =  'subroutine  READ_NAMELIST,time dependent Newton solver only valid for glf23'
     error(36)%errno     =  36

     error(37)%descrip   =  'subroutine  READ_NAMELIST,neoclassical transport,jneo,not set correctly'
     error(37)%errno     =  37

     error(38)%descrip   =  'subroutine  READ_NAMELIST,Maccormack method not useable for this case'
     error(38)%errno     =  38

     error(39)%descrip   =  'subroutine  READ_NAMELIST,use_flow_eq must be explicitely set for this case'
     error(39)%errno     =  39


     error(40)%descrip   =  'subroutine  set_itran,steady state ambiguous equation selection'
     error(40)%errno     =  40
     error(40)%attrib    =  ' only itran = +1 or itran = -1 is allowed in steady state'

     error(41)%descrip   =  'subroutine iter_dbase_txt,(file rw_iterdb.F90) io error'
     error(41)%errno     =  41
     error(41)%attrib    =  ' Error in processing iterdb type text file input'

     error(42)%descrip   =  'subroutine swim_gcnmp_map,inconsistent nprim specification'
     error(42)%errno     =  42

     error(43)%descrip   =  'subroutine swim_gcnmp_map,inconsistent nimp specification'
     error(43)%errno     =  43


     error(44)%errno     =  44
     error(44)%descrip   =  'subroutine swim_gcnmp_map,inconsistent neutral specification'
     
     error(45)%errno     =  45
     error(45)%descrip   =  'subroutine ion_source_rate,number of neutrals no allowed' 

     error(46)%errno     =  46
     error(46)%descrip   =  'subroutine cxrv argument out of range' 

     error(47)%errno     =  47
     error(47)%descrip   =  'subroutine swim_gcnmp_map edge neutral source does not match  any primary ion species' 

     error(48)%errno     =  48
     error(48)%descrip   =  'subroutine curr_convert (in statefile) problem with dimensions' 

     error(49)%errno     =  49
     error(49)%descrip   =  'subroutine readnet_swim_state problem with character arrays' 



     error(50)%errno     =  50
     error(50)%descrip   =  'subroutine read_namelist nion_max too small'
 
     error(51)%errno     =  51
     error(51)%descrip   =  'subroutine set_forcebal_glf_input:  option not implemented' 

     error(52)%errno     =  52
     error(52)%descrip   =  'subroutine set_forcebal_tglf_input:  option not implemented' 





     error(55)%errno     =  55
     error(55)%descrip   =  'subroutine forcebal detected error,see file msg_**' 
!    

     error(56)%errno     =  56
     error(56)%descrip   =  'subroutine forcebal_input_for_gcnmp,nprim error' 


     error(57)%errno     =  57
     error(57)%descrip   =  'subroutine mmm7_1 returned error' 

     error(58)%errno     =  58
     error(58)%descrip   =  'Error:  tglf grid pt paralleization must be on as well as wave no'

     error(59)%errno     =  59
     error(59)%descrip   =  'Error:  sub test_split_comm,communicator split did not check out'

     error(60)%errno     =  60
     error(60)%descrip   =  'Error:  sub gcnmp_init,test_tglf requires # processors GE # grid points' 

     error(61)%errno     =  61
     error(61)%descrip   =  'Error:  sub collect_tglf_outputv193,3 species allowed: electron,primary ion and impurity' 

     error(62)%errno     =  62
     error(62)%descrip   =  'Error:  sub mmm_load_statefile_vectors,dimension bounds problem' 

     error(63)%errno     =  63
     error(63)%descrip   =  'Error:  sub tport,n=nmax or time =time_max,no steps possible' 

     error(64)%errno     =  64
     error(64)%descrip   =  'Error:  sub set_fixed_transpt option currentluy not available: frz_tglf,frz_glf' 

     error(65)%errno     =  65
     error(65)%descrip   =  'Error:  sub set_fixed_transpt options not consistent with simulation flags'

     error(66)%errno     =  66
     error(66)%descrip   =  'Error: nml input,ion_torrot,ion_den,curden : fixed_transpt not implemented'

     error(67)%errno     =  67
     error(67)%descrip   =  'Error: read_efit_eqdsk,forcebal.f90 mhd data not in state file'

     error(68)%errno     =  68
     error(68)%descrip   =  'Error: sub robinbc requires common boundary'

     error(69)%errno     =  69
     error(69)%descrip   =  'Error: sub ss_boundary_flux,equation not known'

     error(iomaxerr)%descrip   =  'Generic error exit'
     error(iomaxerr)%errno     =  iomaxerr
     error(iomaxerr)%attrib    =  'used to exit code at various places for testing'


!    other errors
     error(iomaxerr+1)%descrip = 'time_max .lt. time0'
     error(iomaxerr+1)%errno   = iomaxerr + 1_I4B

     error(iomaxerr+2)%descrip = 'namelist run_data did not encounter end of namelist'
     error(iomaxerr+2)%errno   = iomaxerr + 2_I4B

     error(iomaxerr+3)%descrip = 'namelist run_data failed in read'
     error(iomaxerr+3)%errno   = iomaxerr + 3_I4B

     error(iomaxerr+4)%descrip = 'itenid dimension too small'
     error(iomaxerr+4)%errno   = iomaxerr + 4_I4B

     error(iomaxerr+5)%descrip = 'itenpd dimension too small'
     error(iomaxerr+5)%errno   = iomaxerr + 5_I4B


     error(iomaxerr+6)%descrip = '0.0 .le. fix_edge .le. 1.0 required'
     error(iomaxerr+6)%errno   = iomaxerr + 6_I4B


     error(iomaxerr+7)%descrip = 'sub get_vars: no. equations to solve = 0'
     error(iomaxerr+7)%errno   = iomaxerr + 7_I4B

 
     error(iomaxerr+8)%descrip = 'sub read_namelist: jbootstrap not recognized'
     error(iomaxerr+8)%errno   = iomaxerr + 8_I4B

     error(iomaxerr+9)%descrip = 'sub read_namelist: mul_den,mul_flux dimension  too small'
     error(iomaxerr+9)%errno   = iomaxerr + 9_I4B

     error(iomaxerr+10)%descrip = 'sub gcnm_solve: no. equations to solve = 0'  !not used
     error(iomaxerr+10)%errno   = iomaxerr + 10_I4B


     error(iomaxerr+12)%descrip = 'sub gcnm_solve: no solution method selected'
     error(iomaxerr+12)%errno   = iomaxerr + 12_I4B


     error(iomaxerr+14)%descrip = 'sub gcnm_solve: global set incorrectly '
     error(iomaxerr+14)%errno   = iomaxerr + 14_I4B

     error(iomaxerr+14)%descrip = 'sub gcnm_solve:  error return '
     error(iomaxerr+14)%errno   = iomaxerr + 14_I4B


     error(iomaxerr+15)%descrip = 'sub gcnm_solve:  coding err, termcode =0 '
     error(iomaxerr+15)%errno   = iomaxerr + 15_I4B



     error(iomaxerr+16)%descrip = 'sub gcnm_solve:  case 3 '
     error(iomaxerr+16)%errno   = iomaxerr + 16_I4B


     error(iomaxerr+17)%descrip = 'sub gcnm_solve: no freeze possible 17 '
     error(iomaxerr+17)%errno   = iomaxerr + 17_I4B


     error(iomaxerr+18)%descrip = 'sub read_namelist: nbctim=0 not allowed7 '
     error(iomaxerr+18)%errno   = iomaxerr + 18_I4B


     error(iomaxerr+20)%descrip = 'sub set_itran: neocl_mult(3,3) = 1 is required '
     error(iomaxerr+20)%errno   = iomaxerr + 20_I4B

     error(iomaxerr+21)%descrip = 'sub read_namelist: mul_den,mul_fluxor gamma_bc  error'
     error(iomaxerr+21)%errno   = iomaxerr + 21_I4B

     error(iomaxerr+30)%descrip = 'sub : set_ion_prop: of primary ion  not recognized '
     error(iomaxerr+30)%errno = iomaxerr + 30_I4B




     error(iomaxerr+31)%descrip = 'sub : set_ion_prop: name of  impurity ion  not recognized '
     error(iomaxerr+31)%errno   = iomaxerr + 31_I4B



     error(iomaxerr+32)%descrip = 'sub :     method not implemented '
     error(iomaxerr+32)%errno   = iomaxerr + 32_I4B



     error(iomaxerr+33)%descrip = 'sub :fdjac testing section'
     error(iomaxerr+33)%errno   = iomaxerr + 33_I4B



     error(iomaxerr+34)%descrip = 'sub nemodel :dormqr error'
     error(iomaxerr+34)%errno   = iomaxerr + 34_I4B

     error(iomaxerr+35)%descrip = 'sub nemodel : R is singular  '
     error(iomaxerr+35)%errno   = iomaxerr + 35_I4B

     error(iomaxerr+36)%descrip = 'sub nl_driver : input error  '
     error(iomaxerr+36)%errno   = iomaxerr + 36_I4B

     error(iomaxerr+37)%descrip = 'sub set_ion_prop: fast- thermal  species dont match'
     error(iomaxerr+37)%errno   = iomaxerr + 37_I4B

     error(iomaxerr+40)%descrip = 'subroutine bc_zone:fix_edge inpt prblm'  
     error(iomaxerr+40)%errno   = iomaxerr + 40_I4B

     error(iomaxerr+45)%descrip = 'subroutine glf_driver_pert'  
     error(iomaxerr+45)%errno   = iomaxerr + 45_I4B
     error(iomaxerr+45)%attrib  = 'index jmm > nj-1 should not have occured'

     error(iomaxerr+50)%descrip = 'diffuse,iptrpa = 3 (dv method)  not implemented'
     error(iomaxerr+50)%errno   = iomaxerr + 50_I4B


     error(iomaxerr+51)%descrip = 'diffuse,iptrpa (dv method)  not implemented'
     error(iomaxerr+51)%errno   = iomaxerr + 51_I4B


     error(iomaxerr+52)%descrip = 'dv method  not implemented for rotation'
     error(iomaxerr+52)%errno   = iomaxerr + 52_I4B


     error(iomaxerr+53)%descrip = 'subroutine neoresist: failed resistivity '
     error(iomaxerr+53)%errno   = iomaxerr + 53_I4B

     error(iomaxerr+54)%descrip = 'ERROR in sub simple_advance'
     error(iomaxerr+54)%errno   =  iomaxerr + 54_I4B
     error(iomaxerr+54)%attrib  = 'solution method not specified or not avaiable for this case'

     error(iomaxerr+55)%descrip = 'ERROR in sub advance_set1'
     error(iomaxerr+55)%errno   =  iomaxerr + 55_I4B
     error(iomaxerr+55)%attrib  = 'solution method not specified'


     error(iomaxerr+56)%descrip = 'ERROR in sub advance_set2'
     error(iomaxerr+56)%errno   =  iomaxerr + 56_I4B
     error(iomaxerr+56)%attrib  = 'solution method not specified'


     error(iomaxerr+57)%descrip = 'ERROR in subroutine mmm7_1: nerr .NE. 0 was returned '
     error(iomaxerr+57)%errno   =  iomaxerr + 57_I4B
     error(iomaxerr+57)%attrib  = 'mmm7_1 trapped error'

     error(iomaxerr+75)%descrip = 'subroutine setup_solve_equations,nkt =0  '
     error(iomaxerr+75)%errno   = iomaxerr + 75_I4B

     error(iomaxerr+80)%descrip = 'subroutine check_zero_bc error  '
     error(iomaxerr+80)%errno   = iomaxerr + 80_I4B

     error(iomaxerr+81)%descrip = 'subroutine check_zero_bc error  '
     error(iomaxerr+81)%errno   = iomaxerr + 81_I4B


     error(iomaxerr+85)%descrip = 'Inquiry failure in sub ohacd  '
     error(iomaxerr+85)%errno   = iomaxerr + 85_I4B


     error(iomaxerr+86)%descrip = 'vloop_monitor.txt  open failure  in sub ohacd  '
     error(iomaxerr+86)%errno   = iomaxerr + 86_I4B


     error(iomaxerr+87)%descrip = 'ERROR in sub ohacd writting file '
     error(iomaxerr+87)%errno   = iomaxerr + 87_I4B


     error(iomaxerr+88)%descrip = 'ERROR in sub pturb_vars,isym out of range '
     error(iomaxerr+88)%errno   = iomaxerr + 88_I4B
     error(iomaxerr+88)%attrib  = 'isym currently rangesfrom 1 to itran_max'

     error(iomaxerr+89)%descrip = 'ERROR in sub select_par_model,solver selection '
     error(iomaxerr+89)%errno   = iomaxerr + 89_I4B
     error(iomaxerr+89)%attrib  = 'mmm, ncpu=1 case, solver not defined ?? '


     error(iomaxerr+90)%descrip = 'ERROR in sub fdaysource,jboot value not implemented'
     error(iomaxerr+90)%errno   = iomaxerr + 90_I4B

     error(iomaxerr+91)%descrip = 'ERROR in sub zfit/zsqfit,max temperature range exceeded'
     error(iomaxerr+91)%errno   = iomaxerr + 91_I4B

     error(iomaxerr+92)%descrip = 'ERROR in sub zfit/zsqfit,impurity not recognized'
     error(iomaxerr+92)%errno   = iomaxerr + 92_I4B

     error(iomaxerr+95)%descrip = 'ERROR in sub dgbsv,pred-cor sol, info ne 0'
     error(iomaxerr+95)%errno   = iomaxerr + 95_I4B

     error(iomaxerr+96)%descrip = 'ERROR in sub nwt_pred_cor incorrect input settings'
     error(iomaxerr+96)%errno   = iomaxerr + 96_I4B

     error(iomaxerr+97)%descrip = 'ERROR in sub zsqfit, impurity names not consistent'
     error(iomaxerr+97)%errno   = iomaxerr + 97_I4B

     error(iomaxerr+98)%descrip = 'ERROR in sub zfit, impurity names not consistent'
     error(iomaxerr+98)%errno   = iomaxerr + 98_I4B

     error(iomaxerr+99)%descrip = 'ERROR  sub set_mmm_vars_zct is not available'
     error(iomaxerr+99)%errno   = iomaxerr + 99_I4B

     error(iomaxerr+100)%descrip = 'ERROR non existent error '
     error(iomaxerr+100)%errno   = iomaxerr + 100_I4B

     error(iomaxerr+101)%descrip = 'ERROR in sub fdjacp_sparse,numprocs error'
     error(iomaxerr+101)%errno   = iomaxerr + 101_I4B
     error(iomaxerr+101)%attrib  = 'internal error'

     error(iomaxerr+102)%descrip = 'ERROR in sub fdjacp_sparse,numprocs error'
     error(iomaxerr+102)%errno   =  iomaxerr + 102_I4B
     error(iomaxerr+102)%attrib  = 'must run either 1 process or at least enough processors to cover the Jacobian'


     error(iomaxerr+103)%descrip = 'ERROR in sub read_namelist,parallel_model error'
     error(iomaxerr+103)%errno   =  iomaxerr + 102_I4B
     error(iomaxerr+103)%attrib  = 'jacobian_type = 1 is implemented only with parallel_model = 0'

     error(iomaxerr+104)%descrip = 'ERROR in sub finite_elements , k profile index problem'
     error(iomaxerr+104)%errno   =  iomaxerr + 104_I4B
     error(iomaxerr+104)%attrib  = 'K index is used because boundary of profile is not at rho = a'

     error(iomaxerr+105)%descrip = 'ERROR in sub abcg,selected boundary condition not implemented'
     error(iomaxerr+105)%errno   =  iomaxerr + 105_I4B

     error(iomaxerr+106)%descrip = 'ERROR in sub scale_densities,too many ions densities selected for simulation'
     error(iomaxerr+106)%attrib  = 'only one density can be selected for simulation in this case.'
     error(iomaxerr+106)%errno   =  iomaxerr + 106_I4B

     error(iomaxerr+107)%descrip = 'ERROR in sub eval_diffusive_eq , eq_split option not implemented' 
     error(iomaxerr+107)%attrib  = 'inconsistent input is most likely cause of this error'
     error(iomaxerr+107)%errno   =  iomaxerr + 107_I4B

     error(iomaxerr+108)%descrip = 'ERROR in sub regrid , unexpected size in r' 
     error(iomaxerr+108)%attrib  = 'this is a programming error'
     error(iomaxerr+108)%errno   =  iomaxerr + 108_I4B

     error(iomaxerr+109)%descrip = 'ERROR in sub scale_densities ,array reallocation not allowed' 
     error(iomaxerr+109)%attrib  = 'encountered when dynaic gridding is in effect'
     error(iomaxerr+109)%errno   =  iomaxerr + 109_I4B

     error(iomaxerr+110)%descrip = 'ERROR in sub com_group_init ,' 
     error(iomaxerr+110)%attrib  = 'this is a programming error, should not be encountered'
     error(iomaxerr+110)%errno   =  iomaxerr + 110_I4B

     error(iomaxerr+115)%descrip = 'ERROR in sub set_tglf_vars,' 
     error(iomaxerr+115)%attrib  = 'detected profile input error: all values equal to zero ?'
     error(iomaxerr+115)%errno   =  iomaxerr + 115_I4B

     error(iomaxerr+111)%descrip = 'ERROR in sub com_group_init ,' 
     error(iomaxerr+111)%attrib  = 'Number of processors not commensurate with no of communicator groups'
     error(iomaxerr+111)%errno   =  iomaxerr + 111_I4B

     error(iomaxerr+112)%descrip = 'ERROR in sub distribute_namelist ,' 
     error(iomaxerr+112)%attrib  = 'character arry size problem gcnmp versus forecebal'
     error(iomaxerr+112)%errno   =  iomaxerr + 112_I4B


     error(iomaxerr+120)%descrip = 'ERROR in sub set_tglf_vars,' 
     error(iomaxerr+120)%attrib  = 'internal error,zeff not calculated correctly'
     error(iomaxerr+120)%errno   =  iomaxerr + 120_I4B

     error(iomaxerr+122)%descrip = 'ERROR in sub set_tglf_vars or setup_models,' 
     error(iomaxerr+122)%attrib  = 'igeo_tg geometry switch = 2,3(Fourier or Elite) not implemented'
     error(iomaxerr+122)%errno   =  iomaxerr + 122_I4B

     error(iomaxerr+123)%descrip = 'ERROR in sub set_parallel_model option  branch not defined,' 
     error(iomaxerr+123)%attrib  = 'no turbulent model selected'
     error(iomaxerr+123)%errno   =  iomaxerr + 123_I4B

     error(iomaxerr+124)%descrip = 'ERROR in sub tglf_driver,test_tglf runs only with 1 cpu,' 
     error(iomaxerr+124)%attrib  = 'Multiple cpu option not implemented'
     error(iomaxerr+124)%errno   =  iomaxerr + 124_I4B

     error(iomaxerr+125)%descrip = 'ERROR in sub tglf_driver,namelist tglfin file problem' 
     error(iomaxerr+125)%attrib  = 'namelist file read error,' 
     error(iomaxerr+125)%errno   =  iomaxerr + 125_I4B

     error(iomaxerr+126)%descrip = 'ERROR in sub tglf_driver,namelist tglfin namelist read problem' 
     error(iomaxerr+126)%attrib  = 'namelist file read error,' 
     error(iomaxerr+126)%errno   =  iomaxerr + 126_I4B

     error(iomaxerr+127)%descrip = 'ERROR in sub compact_deriv' 
     error(iomaxerr+127)%attrib  = 'non uniform mesh requires that order = 2 or 4,(received -1)'
     error(iomaxerr+127)%errno   =  iomaxerr + 127_I4B

     error(iomaxerr+128)%descrip = 'ERROR in sub compact_deriv' 
     error(iomaxerr+128)%attrib  = 'non uniform mesh requires that order = 2 or 4,(received -2)' 
     error(iomaxerr+128)%errno   =  iomaxerr + 128_I4B


     error(iomaxerr+129)%descrip = 'ERROR in sub compact_deriv' 
     error(iomaxerr+129)%attrib  = 'Lapack dgtsv returned non zero value for info,order =4 ' 
     error(iomaxerr+129)%errno   =  iomaxerr + 129_I4B   

     error(iomaxerr+130)%descrip = 'ERROR in sub compact_deriv' 
     error(iomaxerr+130)%attrib  = 'unifrom grid selected for 4th order method, (must use nonuniform option)' 
     error(iomaxerr+130)%errno   =  iomaxerr + 130_I4B   

     error(iomaxerr+131)%descrip = 'ERROR in sub compact_deriv' 
     error(iomaxerr+131)%attrib  = 'non-unifrom grid selected for 6th order method, (must use uniform option)' 
     error(iomaxerr+131)%errno   =  iomaxerr + 131_I4B   


     error(iomaxerr+132)%descrip = 'ERROR in sub compact_deriv' 
     error(iomaxerr+132)%attrib  = 'unifrm zctr grid selected use unifrm fctr' 
     error(iomaxerr+132)%errno   =  iomaxerr + 132_I4B   

     error(iomaxerr+133)%descrip = 'ERROR in sub compact_interp' 
     error(iomaxerr+133)%attrib  = 'lapack dgbsv returned non zero error code' 
     error(iomaxerr+133)%errno   =  iomaxerr + 133_I4B   


     error(iomaxerr+134)%descrip = 'ERROR in sub tglf_pack_data' 
     error(iomaxerr+134)%attrib  = 'Amount of data to pack  does  not match parameter tglf_max_pack' 
     error(iomaxerr+134)%errno   =  iomaxerr + 134_I4B   


     error(iomaxerr+135)%descrip = 'ERROR in sub single_comm_group' 
     error(iomaxerr+135)%attrib  = 'Amount of data received  does  not match parameter tglf_max_pack' 
     error(iomaxerr+135)%errno   =  iomaxerr + 135_I4B   


     error(iomaxerr+136)%descrip = 'ERROR in sub tglf_store_packed' 
     error(iomaxerr+136)%attrib  = 'Amount of data stored  does  not match parameter tglf_max_pack' 
     error(iomaxerr+136)%errno   =  iomaxerr + 136_I4B   

     error(iomaxerr+137)%descrip = 'ERROR in sub get_M' 
     error(iomaxerr+137)%attrib  = 'rho grid index j=1 or j =nj not allowed' 
     error(iomaxerr+137)%errno   =  iomaxerr + 137_I4B  
 
     error(iomaxerr+138)%descrip = 'ERROR in sub get_W' 
     error(iomaxerr+138)%attrib  = 'rho grid index j=1 or j =nj not allowed' 
     error(iomaxerr+138)%errno   =  iomaxerr + 138_I4B   

     error(iomaxerr+139)%descrip = 'ERROR in sub set_bdry_index' 
     error(iomaxerr+139)%attrib  = 'all bdry index values must be the same for this case' 
     error(iomaxerr+139)%errno   =  iomaxerr + 139_I4B   

     error(iomaxerr+140)%descrip = 'ERROR in sub set_bdry_index' 
     error(iomaxerr+140)%attrib  = 'all bdry index values must be the same for this case' 
     error(iomaxerr+140)%errno   =  iomaxerr + 140_I4B   

     error(iomaxerr+141)%descrip = 'ERROR in sub refactor_profiles/tglf_flux_processing' 
     error(iomaxerr+141)%attrib  = 'tglf_p,e_flux/glf_p,e,m_flux array size problem' 
     error(iomaxerr+141)%errno   =  iomaxerr + 141_I4B 
  
     error(iomaxerr+142)%descrip = 'ERROR in sub single_comm_group_domain_decomp' 
     error(iomaxerr+142)%attrib  = 'subroutine requires that MPI2_EXT is used in compilation' 
     error(iomaxerr+142)%errno   =  iomaxerr + 142_I4B 

     error(iomaxerr+145)%descrip = 'ERROR in sub compact_deriv' 
     error(iomaxerr+145)%attrib  = 'dgtsv returned non xero info value,order =6 section' 
     error(iomaxerr+145)%errno   =  iomaxerr + 145_I4B 

     error(iomaxerr+146)%descrip = 'ERROR in sub midpt_interp' 
     error(iomaxerr+146)%attrib  = 'dgbsv returned non zero info value' 
     error(iomaxerr+146)%errno   =  iomaxerr + 146_I4B 

     error(iomaxerr+148)%descrip = 'ERROR in sub compact_deriv' 
     error(iomaxerr+148)%attrib  = 'order=4,zctr_deriv = 1 option not implemented yet' 
     error(iomaxerr+148)%errno   =  iomaxerr + 148_I4B 

     error(iomaxerr+149)%descrip = 'ERROR in sub tglf/mmm _single_grid_point' 
     error(iomaxerr+149)%attrib  = 'pertrb  option not implemented' 
     error(iomaxerr+149)%errno   =  iomaxerr + 149_I4B 



     error(iomaxerr+150)%descrip = 'ERROR in sub single_comm_group_ms' 
     error(iomaxerr+150)%attrib  = 'number of processors must be <= rho grid size -1,(eg nj-1) ' 
     error(iomaxerr+150)%errno   =  iomaxerr + 150_I4B 



     error(iomaxerr+151)%descrip = 'ERROR in sub glf_pack_data' 
     error(iomaxerr+151)%attrib  = 'Amount of data to pack  does  not match parameter glf_max_pack' 
     error(iomaxerr+151)%errno   =  iomaxerr + 151_I4B   


     error(iomaxerr+152)%descrip = 'ERROR in sub single_comm_group'
     error(iomaxerr+152)%attrib  = 'Amount of data received  does  not match parameter glf_max_pack' 
     error(iomaxerr+152)%errno   =  iomaxerr + 152_I4B   


     error(iomaxerr+153)%descrip = 'ERROR in sub glf_store_packed' 
     error(iomaxerr+153)%attrib  = 'Amount of data stored  does  not match parameter glf_max_pack' 
     error(iomaxerr+153)%errno   =  iomaxerr + 153_I4B   



     error(iomaxerr+154)%descrip = 'ERROR in sub get_tot_pressure' 
     error(iomaxerr+154)%attrib  = 'profile%ene,profile%te,etc not alloacated' 
     error(iomaxerr+154)%errno   =  iomaxerr + 154_I4B 


     error(iomaxerr+155)%descrip = 'ERROR in sub glf_driver_pert' 
     error(iomaxerr+155)%attrib  = 'itport_pt(4) < 0 or itport_pt(5) .NE. 0 is not implemented' 
     error(iomaxerr+155)%errno   =  iomaxerr + 155_I4B 

     error(iomaxerr+156)%descrip = 'ERROR in sub process_shape' 
     error(iomaxerr+156)%attrib  = 'rmin at zelev1  and rmax at zelev2 not consistent' 
     error(iomaxerr+156)%errno   =  iomaxerr + 156_I4B 

     error(iomaxerr+157)%descrip = 'ERROR in sub scale_betan_ped' 
     error(iomaxerr+157)%attrib  = 'rohn_betan_ped is out of range (0,1)' 
     error(iomaxerr+157)%errno   =  iomaxerr + 157_I4B 

     error(iomaxerr+158)%descrip = 'ERROR in sub scale_betan_ped' 
     error(iomaxerr+158)%attrib  = 'negative density in betan ped scaling' 
     error(iomaxerr+158)%errno   =  iomaxerr + 158_I4B 

     error(iomaxerr+159)%descrip = 'ERROR in sub get_den_ptrb' 
     error(iomaxerr+159)%attrib  = 'info from DGETRF or DGETRS .NE. 0' 
     error(iomaxerr+159)%errno   =  iomaxerr + 159_I4B 

     error(iomaxerr+160)%descrip = 'ERROR in sub eval_flow_eq , eq_split option not implemented' 
     error(iomaxerr+160)%attrib  = 'equation set splitting is not yet ready for use in this context'
     error(iomaxerr+160)%errno   =  iomaxerr + 160_I4B


     error(iomaxerr+161)%descrip = 'ERROR in sub tglf_allocate , itran_max .NE. itprt_max' 
     error(iomaxerr+161)%attrib  = 'dimensional alignment problem - fatal error'
     error(iomaxerr+161)%errno   =  iomaxerr + 161_I4B

     error(iomaxerr+162)%descrip = 'ERROR in sub calc_densities , nspline too small' 
     error(iomaxerr+162)%attrib  = '      nspline  must be .GE.  3'
     error(iomaxerr+162)%errno   =  iomaxerr + 162_I4B


     error(iomaxerr+163)%descrip = 'ERROR in sub single_comm_group' 
     error(iomaxerr+163)%attrib  = 'Amount of data received  does  not match parameter mmm_max_pack' 
     error(iomaxerr+163)%errno   =  iomaxerr + 163_I4B   


     error(iomaxerr+165)%descrip = 'ERROR in sub flow_eq_resid, grid index out of range' 
     error(iomaxerr+165)%attrib  = 'valid grid indecies vary from 1 to nj depending on dependent variable'
     error(iomaxerr+165)%errno   =  iomaxerr + 165_I4B


     error(iomaxerr+166)%descrip = 'ERROR in sub flow_eq_resid, grid index 1 not valid for rbp' 
     error(iomaxerr+166)%attrib  = 'valid grid indecies vary from 2 to nj for rbp'
     error(iomaxerr+166)%errno   =  iomaxerr + 166_I4B


     error(iomaxerr+168)%descrip = 'ERROR in sub Jac_ss_flow, grid index ot of range' 
     error(iomaxerr+168)%errno   =  iomaxerr + 168_I4B

     error(iomaxerr+169)%errno   =  iomaxerr + 169_I4B
     error(iomaxerr+169)%descrip = 'subroutine source: stfuse not set in statefile  but required' 
     error(iomaxerr+169)%attrib  = 'use "internal_thermal_fusion = .TRUE."  for this case'


     error(iomaxerr+170)%descrip = 'ERROR in sub tlf_ssfdjac_par, compact schemes not implemented yet ' 
     error(iomaxerr+170)%errno   =  iomaxerr + 170_I4B


     error(iomaxerr+171)%descrip = 'ERROR in namelist input for tglf species, must be > = 2' 
     error(iomaxerr+171)%errno   =  iomaxerr + 171_I4B



     error(iomaxerr+175)%descrip = 'ERROR in sub check_scale_factrs, zero scaling encountered' 
     error(iomaxerr+175)%errno   =  iomaxerr + 175_I4B


     error(iomaxerr+176)%descrip = 'ERROR in sub terminate_iterations, situaton should not happen' 
     error(iomaxerr+176)%errno   =  iomaxerr + 176_I4B

     error(iomaxerr+177)%errno   =  iomaxerr + 177_I4B
     error(iomaxerr+177)%descrip = 'subroutine source: ion species of type dt' 
     error(iomaxerr+177)%attrib  = 'use "internal_thermal_fusion = .TRUE."  for this case'
     error(iomaxerr+178)%errno   =  iomaxerr + 178_I4B
     error(iomaxerr+178)%descrip = 'subroutine :collect_stats, start and stop problem' 
     error(iomaxerr+178)%attrib  = 'cant call with start and stop set simultaneously'

     error(iomaxerr+179)%errno   =  iomaxerr + 179_I4B
     error(iomaxerr+179)%descrip = 'subroutine :collect_stats, stat_index out of range' 
     error(iomaxerr+179)%attrib  = 'require 1 LE stat_index LE nstats'



     error(iomaxerr+180)%descrip = 'ERROR in sub ss_fdjac_par, master id must be 0' 
     error(iomaxerr+180)%attrib  = 'requires redefining master internally and recompiling'
     error(iomaxerr+180)%errno   =  iomaxerr + 180_I4B

     error(iomaxerr+181)%descrip = 'ERROR in sub simple_advance, time dependent option not available' 
     error(iomaxerr+181)%attrib  = 'Newton method combined with steady_state = 1 not allowed here'
     error(iomaxerr+181)%errno   =  iomaxerr + 181_I4B

     error(iomaxerr+182)%descrip = 'Error in sub get_tglf_flux_pert, must be compiled with mpi'
     error(iomaxerr+182)%attrib  = 'Tglflux perturbations are done only in parallel mode'
     error(iomaxerr+182)%errno   =  iomaxerr + 182_I4B

     error(iomaxerr+183)%descrip = 'Error in sub set_event_error'
     error(iomaxerr+183)%attrib  = '  event in cur_events not recognized'
     error(iomaxerr+183)%errno   =  iomaxerr + 183_I4B


     error(iomaxerr+184)%descrip = 'Error in sub spl2bc'
     error(iomaxerr+184)%attrib  = 'spli2d returned iflag .NE. 0 '
     error(iomaxerr+184)%errno   =  iomaxerr + 184_I4B

     error(iomaxerr+185)%descrip = 'Error in sub mmm_check_input'
     error(iomaxerr+185)%attrib  = 'Weights for Wieland,DRBM,ETG not set correctly in multimode model '
     error(iomaxerr+185)%errno   =  iomaxerr + 185_I4B

! temporary errors:
     error(iomaxerr+186)%descrip = 'ERROR in sub simple_advance' 
     error(iomaxerr+186)%attrib  = 'Full Newton solver selected for mmm  implies ss not yet implemented'
     error(iomaxerr+186)%errno   =  iomaxerr + 186_I4B

     error(iomaxerr+187)%descrip = 'ERROR in sub set_mmm_vars' 
     error(iomaxerr+187)%attrib  = 'use_compact_schemes not yet implemented'
     error(iomaxerr+187)%errno   =  iomaxerr + 187_I4B

     error(iomaxerr+188)%descrip = 'ERROR in sub psi_contours' 
     error(iomaxerr+188)%attrib  = '2d spline failure'
     error(iomaxerr+188)%errno   =  iomaxerr + 188_I4B


     error(iomaxerr+189)%descrip = 'ERROR in sub mmm_single_grid_point' 
     error(iomaxerr+189)%attrib  = 'grid point error: must have 1 .LE. grid point .LE. nj-1'
     error(iomaxerr+189)%errno   =  iomaxerr + 189_I4B


     error(iomaxerr+190)%descrip = 'ERROR in sub parabolic_prof_set' 
     error(iomaxerr+190)%attrib  = 'n_r does not match nj'
     error(iomaxerr+190)%errno   =  iomaxerr + 190_I4B

     error(iomaxerr+192)%descrip = 'ERROR deliberate termination,' 
     error(iomaxerr+192)%attrib  = 'for testing '
     error(iomaxerr+192)%errno   =  iomaxerr + 192_I4B

     error(iomaxerr+193)%descrip = 'ERROR in sub tglf_single_grid_point,' 
     error(iomaxerr+193)%attrib  = 'value of prtrb is out of range '
     error(iomaxerr+193)%errno   =  iomaxerr + 193_I4B

     error(iomaxerr+194)%descrip = 'ERROR in sub jac_grad_u. branch not yet defined,' 
     error(iomaxerr+194)%attrib  = 'temporary error during code build'
     error(iomaxerr+194)%errno   =  iomaxerr + 194_I4B

     error(iomaxerr+195)%descrip = 'ERROR in sub jac_grad_u.include_glf branch not yet defined,' 
     error(iomaxerr+195)%attrib  = 'temporary error during code build'
     error(iomaxerr+195)%errno   =  iomaxerr + 195_I4B

     error(iomaxerr+196)%descrip = 'ERROR in sub jac_u. branch not yet defined,' 
     error(iomaxerr+196)%attrib  = 'temporary error during code build'
     error(iomaxerr+196)%errno   =  iomaxerr + 196_I4B

     error(iomaxerr+197)%descrip = 'ERROR in sub jac_u.include_tglf branch not yet defined,' 
     error(iomaxerr+197)%attrib  = 'temporary error during code build'
     error(iomaxerr+197)%errno   =  iomaxerr + 197_I4B

     error(iomaxerr+198)%descrip = 'ERROR in sub get_fluxes.include_glf branch not yet defined,' 
     error(iomaxerr+198)%attrib  = 'temporary error during code build'
     error(iomaxerr+198)%errno   =  iomaxerr + 198_I4B


     error(iomaxerr+199)%descrip = 'ERROR in sub linearized_solver sub set_glf_vars not yet defined,' 
     error(iomaxerr+199)%attrib  = 'temporary error during code build'
     error(iomaxerr+199)%errno   =  iomaxerr + 199_I4B



!---------------------------------------------------------------------------------
!   stand alone nfreya specific errors:
!---------------------------------------------------------------------------------

    nfreya_offset = 200 + iomaxerr
 
 
    error(nfreya_offset)%descrip  = 'insufficient  command line arguments found'
    error(nfreya_offset)%attrib   = 'usage: P_Nfreya statefile  namelist_file run_log_file'
    error(nfreya_offset)%errno    =  iomaxerr + 200_I4B

    nfreya_offset = 201 + iomaxerr
    error(nfreya_offset)%descrip  = 'NFREYA ERROR: sub nbsgxn: ' 
    error(nfreya_offset)%attrib   = 'Only iexcit =5 is programmed'
    error(nfreya_offset)%errno    =  iomaxerr + 201_I4B

    nfreya_offset = 202 + iomaxerr
    error(nfreya_offset)%descrip  = 'NFREYA ERROR: sub adassgxn : ' 
    error(nfreya_offset)%attrib   = 'local Adas only has  H, He, B, Be, C, O, N,Ne ions'
    error(nfreya_offset)%errno    =  iomaxerr + 202_I4B

    nfreya_offset = 203 + iomaxerr
    error(nfreya_offset)%descrip  = 'NFREYA ERROR: sub tozone: ' 
    error(nfreya_offset)%attrib   = 'zonal grid generation failed'
    error(nfreya_offset)%errno    =  iomaxerr + 203_I4B


    nfreya_offset = 204 + iomaxerr
    error(nfreya_offset)%descrip  = 'NFREYA ERROR: sub prenub : ' 
    error(nfreya_offset)%attrib   = 'nr,nz do not mach values in param module'
    error(nfreya_offset)%errno    =  iomaxerr + 204_I4B

    nfreya_offset = 205 + iomaxerr
    error(nfreya_offset)%descrip  = 'NFREYA ERROR: sub prenub : ' 
    error(nfreya_offset)%attrib   = 'psi axis psi lim values not consistent'
    error(nfreya_offset)%errno    =  iomaxerr + 205_I4B

    nfreya_offset = 206 + iomaxerr
    error(nfreya_offset)%descrip  = 'NFREYA ERROR: sub cntour1 : ' 
    error(nfreya_offset)%attrib   = 'dbcevl1 bicubic spline eval error'
    error(nfreya_offset)%errno    =  iomaxerr + 206_I4B


    nfreya_offset = 207 + iomaxerr
    error(nfreya_offset)%descrip  = 'NFREYA ERROR: sub cntour1 : ' 
    error(nfreya_offset)%attrib   = 'requires max psi on axis'
    error(nfreya_offset)%errno    =  iomaxerr + 207_I4B

    nfreya_offset = 206 + iomaxerr
    error(nfreya_offset)%descrip  = 'NFREYA ERROR: sub cntour2 : ' 
    error(nfreya_offset)%attrib   = 'dbcevl1 bicubic spline eval error'
    error(nfreya_offset)%errno    =  iomaxerr + 208_I4B


    nfreya_offset = 209 + iomaxerr
    error(nfreya_offset)%descrip  = 'NFREYA ERROR: sub cntour2 : ' 
    error(nfreya_offset)%attrib   = 'requires max psi on axis'
    error(nfreya_offset)%errno    =  iomaxerr + 209_I4B


    nfreya_offset = 210 + iomaxerr
    error(nfreya_offset)%descrip  = 'NFREYA ERROR: sub CNTOUR: use_cnt problem: ' 
    error(nfreya_offset)%attrib   = 'one of (use_cnt1,use_cnt2) must be true'
    error(nfreya_offset)%errno    =  iomaxerr + 210_I4B


    nfreya_offset = 211 + iomaxerr
    error(nfreya_offset)%descrip  = 'NFREYA ERROR: sub cntour2 : ' 
    error(nfreya_offset)%attrib   = 'requires max psi on axis'
    error(nfreya_offset)%errno    =  iomaxerr + 211_I4B


    nfreya_offset = 212 + iomaxerr
    error(nfreya_offset)%descrip  = 'NFREYA ERROR: sub cntour : ' 
    error(nfreya_offset)%attrib   = 'cntour1 or 2 returned ierr ne 0'
    error(nfreya_offset)%errno    =  iomaxerr + 212_I4B


    nfreya_offset = 213 + iomaxerr
    error(nfreya_offset)%descrip  = 'NFREYA ERROR: sub prenub : ' 
    error(nfreya_offset)%attrib   = 'cntour  or returned ierr ne 0'
    error(nfreya_offset)%errno    =  iomaxerr + 213_I4B


    nfreya_offset = 214 + iomaxerr
    error(nfreya_offset)%descrip  = 'NFREYA ERROR: sub read_namelist : ' 
    error(nfreya_offset)%attrib   = 'nameb not given,must be one of h,d,dt,t'
    error(nfreya_offset)%errno    =  iomaxerr + 214_I4B


    nfreya_offset = 215 + iomaxerr
    error(nfreya_offset)%descrip  = 'NFREYA ERROR: sub set_ion_prop : ' 
    error(nfreya_offset)%attrib   = 'primary ion name not found in master list'
    error(nfreya_offset)%errno    =  iomaxerr + 215_I4B


    nfreya_offset = 216 + iomaxerr
    error(nfreya_offset)%descrip  = 'NFREYA ERROR: sub cntour1 : ' 
    error(nfreya_offset)%attrib   = 'primary ion name not found in master list'
    error(nfreya_offset)%errno    =  iomaxerr + 216_I4B

    nfreya_offset = 217 + iomaxerr
    error(nfreya_offset)%descrip  = 'NFREYA ERROR: sub cntour1 : ' 
    error(nfreya_offset)%attrib   = 'primary ion name not found in master list'
    error(nfreya_offset)%errno    =  iomaxerr + 217_I4B

    nfreya_offset = 218 + iomaxerr
    error(nfreya_offset)%descrip  = 'NFREYA ERROR: sub cntour2 : newti problem,dxx .LE. dxmin ' 
    error(nfreya_offset)%attrib   = 'primary ion name not found in master list'
    error(nfreya_offset)%errno    =  iomaxerr + 218_I4B

    nfreya_offset = 219 + iomaxerr
    error(nfreya_offset)%descrip  = 'NFREYA ERROR: sub nbsgxn : master impurity list error' 
    error(nfreya_offset)%attrib   = 'impurity ion name not found in master list'
    error(nfreya_offset)%errno    =  iomaxerr + 219_I4B

    nfreya_offset = 220 + iomaxerr
    error(nfreya_offset)%descrip  = 'NFREYA ERROR: sub nbsgxn : hexnb error' 
    error(nfreya_offset)%attrib   = 'file coronb not found'
    error(nfreya_offset)%errno    =  iomaxerr + 220_I4B

    nfreya_offset = 221 + iomaxerr
    error(nfreya_offset)%descrip  = 'NFREYA ERROR: subroutine ADASSGXN' 
    error(nfreya_offset)%attrib   = 'unallowed impurity'
    error(nfreya_offset)%errno    =  iomaxerr + 221_I4B

    nfreya_offset = 222 + iomaxerr
    error(nfreya_offset)%descrip  = 'NFREYA ERROR: subroutine ADASQH6' 
    error(nfreya_offset)%attrib   = 'ADASSGXN called adasqh6 but got  problem #1'
    error(nfreya_offset)%errno    =  iomaxerr + 222_I4B

    nfreya_offset = 223 + iomaxerr
    error(nfreya_offset)%descrip  = 'NFREYA ERROR: subroutine ADASQH6' 
    error(nfreya_offset)%attrib   = 'ADASSGXN called adasqh6 but got  problem #2'
    error(nfreya_offset)%errno    =  iomaxerr + 223_I4B

    nfreya_offset = 224 + iomaxerr
    error(nfreya_offset)%descrip  = 'TSPLINE module error: subroutine seval' 
    error(nfreya_offset)%attrib   = 'of bounds interpolation is occuring'
    error(nfreya_offset)%errno    =  iomaxerr + 224_I4B


    nfreya_offset = 225 + iomaxerr
    error(nfreya_offset)%descrip  = 'NFREYA ERROR: sub set_ion_prop : ' 
    error(nfreya_offset)%attrib   = 'impurity ion name not found in master list'
    error(nfreya_offset)%errno    =  iomaxerr + 225_I4B


    nfreya_offset = 226 + iomaxerr
    error(nfreya_offset)%descrip  = 'NFREYA ERROR: sub rotate : ' 
    error(nfreya_offset)%attrib   = 'geometry specification error'
    error(nfreya_offset)%errno    =  iomaxerr + 226_I4B


    nfreya_offset = 227 + iomaxerr
    error(nfreya_offset)%descrip  = 'NFREYA ERROR: sub nbdep2 : ' 
    error(nfreya_offset)%attrib   = 'ndum >  2*nh not possible'
    error(nfreya_offset)%errno    =  iomaxerr + 227_I4B


    nfreya_offset = 228 + iomaxerr
    error(nfreya_offset)%descrip  = 'NFREYA ERROR: sub getsgxn : ' 
    error(nfreya_offset)%attrib   = 'maximum rotational energy bin exceeded'
    error(nfreya_offset)%errno    =  iomaxerr + 228_I4B


    nfreya_offset = 229 + iomaxerr
    error(nfreya_offset)%descrip  = 'NFREYA ERROR: sub read_namelist : ' 
    error(nfreya_offset)%attrib   = 'total number of pseudo neutrals to follow into plasma not specified '
    error(nfreya_offset)%errno    =  iomaxerr + 229_I4B

    nfreya_offset = 230 + iomaxerr
    error(nfreya_offset)%descrip  = 'NFREYA ERROR: sub read_namelist : ' 
    error(nfreya_offset)%attrib   = 'total number of injectors(eg beamlines)  not specified '
    error(nfreya_offset)%errno    =  iomaxerr + 230_I4B


    nfreya_offset = 231 + iomaxerr
    error(nfreya_offset)%descrip  = 'NFREYA ERROR: sub read_namelist : ' 
    error(nfreya_offset)%attrib   = 'start and or end time of simulation not specified'
    error(nfreya_offset)%errno    =  iomaxerr + 231_I4B


    nfreya_offset = 232 + iomaxerr
    error(nfreya_offset)%descrip  = 'NFREYA ERROR: sub read_namelist : ' 
    error(nfreya_offset)%attrib   = 'Number of beamlines  parameter kb too small for this case'
    error(nfreya_offset)%errno    =  iomaxerr + 232_I4B


    nfreya_offset = 233 + iomaxerr
    error(nfreya_offset)%descrip  = 'NFREYA ERROR: sub launch_nfreya : ' 
    error(nfreya_offset)%attrib   = 'pseudo injector count error'
    error(nfreya_offset)%errno    =  iomaxerr + 233_I4B


    nfreya_offset = 234 + iomaxerr
    error(nfreya_offset)%descrip  = 'NFREYA ERROR: sub nfreya_store_packed : ' 
    error(nfreya_offset)%attrib   = 'error in number of 0di  items to be unpacked'
    error(nfreya_offset)%errno    =  iomaxerr + 234_I4B


    nfreya_offset = 235 + iomaxerr
    error(nfreya_offset)%descrip  = 'NFREYA ERROR: sub nfreya_store_packed : ' 
    error(nfreya_offset)%attrib   = 'error in number of 1di items to be unpacked'
    error(nfreya_offset)%errno    =  iomaxerr + 235_I4B

    nfreya_offset = 236 + iomaxerr
    error(nfreya_offset)%descrip  = 'NFREYA ERROR: sub nfreya_store_packed : ' 
    error(nfreya_offset)%attrib   = 'error in number of 1dr items to be unpacked'
    error(nfreya_offset)%errno    =  iomaxerr + 236_I4B

    nfreya_offset = 237 + iomaxerr
    error(nfreya_offset)%descrip  = 'NFREYA ERROR: sub nfreya_store_packed : ' 
    error(nfreya_offset)%attrib   = 'error in number of 2dr items to be unpacked'
    error(nfreya_offset)%errno    =  iomaxerr + 237_I4B

    nfreya_offset = 238 + iomaxerr
    error(nfreya_offset)%descrip  = 'NFREYA ERROR: sub nfreya_store_packed : ' 
    error(nfreya_offset)%attrib   = 'error in number of 3dr items to be unpacked'
    error(nfreya_offset)%errno    =  iomaxerr + 238_I4B

    nfreya_offset = 239 + iomaxerr
    error(nfreya_offset)%descrip  = 'NFREYA ERROR: sub nfreya_store_packed : ' 
    error(nfreya_offset)%attrib   = 'error in number of 4dr items to be unpacked'
    error(nfreya_offset)%errno    =  iomaxerr + 239_I4B


    nfreya_offset = 240 + iomaxerr
    error(nfreya_offset)%descrip  = 'NFREYA ERROR: sub write_perf_data : ' 
    error(nfreya_offset)%attrib   = 'fatal error appending to write_performance_data file'
    error(nfreya_offset)%errno    =  iomaxerr + 240_I4B

    nfreya_offset = 241 + iomaxerr
    error(nfreya_offset)%descrip  = 'NFREYA ERROR: sub write_perf_data : ' 
    error(nfreya_offset)%attrib   = 'fatal error creating write_performance_data file'
    error(nfreya_offset)%errno    =  iomaxerr + 241_I4B

    nfreya_offset = 242 + iomaxerr
    error(nfreya_offset)%descrip  = 'NFREYA ERROR: sub setup_plot_file ' 
    error(nfreya_offset)%attrib   = 'fatal error ,task 3 call returned error'
    error(nfreya_offset)%errno    =  iomaxerr + 242_I4B

    nfreya_offset = 243 + iomaxerr
    error(nfreya_offset)%descrip  = 'NFREYA ERROR: sub setup_plot_file ' 
    error(nfreya_offset)%attrib   = 'fatal error ,task= -1  call returned error'
    error(nfreya_offset)%errno    =  iomaxerr + 243_I4B

    nfreya_offset = 244 + iomaxerr
    error(nfreya_offset)%descrip  = 'NFREYA ERROR: sub set_profiles ' 
    error(nfreya_offset)%attrib   = 'error in zetaz interpolation'
    error(nfreya_offset)%errno    =  iomaxerr + 244_I4B

    nfreya_offset = 245 + iomaxerr
    error(nfreya_offset)%descrip  = 'NFREYA ERROR: sub read_namelist ' 
    error(nfreya_offset)%attrib   = 'error in statfile_output_name ==> too long'
    error(nfreya_offset)%errno    =  iomaxerr + 245_I4B

    nfreya_offset = 246 + iomaxerr
    error(nfreya_offset)%descrip  = 'NFREYA ERROR: sub strip_path ' 
    error(nfreya_offset)%attrib   = 'no shot number found'
    error(nfreya_offset)%errno    =  iomaxerr + 246_I4B

    nfreya_offset = 247 + iomaxerr
    error(nfreya_offset)%descrip  = 'NFREYA ERROR: sub strip_path ' 
    error(nfreya_offset)%attrib   = 'Error, ufile 64 character limit exceeded'
    error(nfreya_offset)%errno    =  iomaxerr + 247_I4B

    nfreya_offset = 248 + iomaxerr
    error(nfreya_offset)%descrip  = 'NFREYA ERROR: sub get_beam_data ' 
    error(nfreya_offset)%attrib   = 'Error, reading beam_data_namelist'
    error(nfreya_offset)%errno    =  iomaxerr + 248_I4B

    nfreya_offset = 249 + iomaxerr
    error(nfreya_offset)%descrip  = 'NFREYA ERROR: sub get_beam_data ' 
    error(nfreya_offset)%attrib   = 'Error, namelist read failed'
    error(nfreya_offset)%errno    =  iomaxerr + 249_I4B

    nfreya_offset = 250 + iomaxerr
    error(nfreya_offset)%descrip  = 'NFREYA ERROR: sub get_beam_data ' 
    error(nfreya_offset)%attrib   = 'Error, zero beam case not allowed'
    error(nfreya_offset)%errno    =  iomaxerr + 250_I4B

    nfreya_offset = 251 + iomaxerr
    error(nfreya_offset)%descrip  = 'NFREYA ERROR: sub get_beam_data ' 
    error(nfreya_offset)%attrib   = 'Error, beam_data_ufile not found'
    error(nfreya_offset)%errno    =  iomaxerr + 251_I4B

    nfreya_offset = 252 + iomaxerr
    error(nfreya_offset)%descrip  = 'NFREYA ERROR: sub get_beam_data ' 
    error(nfreya_offset)%attrib   = 'Error, get_beam_data returned error'
    error(nfreya_offset)%errno    =  iomaxerr + 252_I4B

    nfreya_offset = 254 + iomaxerr
    error(nfreya_offset)%descrip  = 'NFREYA ERROR: sub get_beam_data ' 
    error(nfreya_offset)%attrib   = 'Error, Expected to read NSC =1 scalar'
    error(nfreya_offset)%errno    =  iomaxerr + 254_I4B

    nfreya_offset = 255 + iomaxerr
    error(nfreya_offset)%descrip  = 'NFREYA ERROR: sub get_beam_data ' 
    error(nfreya_offset)%attrib   = 'Error,namelist read has wrong nbeams'
    error(nfreya_offset)%errno    =  iomaxerr + 255_I4B

    nfreya_offset = 256 + iomaxerr
    error(nfreya_offset)%descrip  = 'NFREYA ERROR: sub set_integer_value' 
    error(nfreya_offset)%attrib   = 'Error,in fetching integer value from string'
    error(nfreya_offset)%errno    =  iomaxerr + 256_I4B

    nfreya_offset = 257 + iomaxerr
    error(nfreya_offset)%descrip  = 'NFREYA ERROR: sub check_d_fast_ion_inpt' 
    error(nfreya_offset)%attrib   = 'Error,in requested array size'
    error(nfreya_offset)%errno    =  iomaxerr + 257_I4B

    nfreya_offset = 258 + iomaxerr
    error(nfreya_offset)%descrip  = 'NFREYA ERROR: sub nblist_read' 
    error(nfreya_offset)%attrib   = 'Error,nbeam >   nbeamx,tbona,tboffa too small '
    error(nfreya_offset)%errno    =  iomaxerr + 258_I4B

    nfreya_offset = 259 + iomaxerr
    error(nfreya_offset)%descrip  = 'NFREYA ERROR: sub nblist_read_parse' 
    error(nfreya_offset)%attrib   = 'Error,nbeam >   nbeamx,tbona,tboffa too small '
    error(nfreya_offset)%errno    =  iomaxerr + 259_I4B

   nfreya_offset = 260 + iomaxerr
    error(nfreya_offset)%descrip  = 'NFREYA ERROR: sub nblist_read_parse' 
    error(nfreya_offset)%attrib   = 'Error,nbeam ndifbe,ndifbep consistency error'
    error(nfreya_offset)%errno    =  iomaxerr + 260_I4B

   nfreya_offset = 265 + iomaxerr
    error(nfreya_offset)%descrip  = 'STRING_UTIL ERROR : sub strip path'
    error(nfreya_offset)%attrib   = 'Error,could not determine shot no.'
    error(nfreya_offset)%errno    =  iomaxerr + 265_I4B


   nfreya_offset = 266 + iomaxerr
    error(nfreya_offset)%descrip  = 'STRING_UTIL ERROR : sub strip path'
    error(nfreya_offset)%attrib   = 'Error ufile limit is 64 characters'
    error(nfreya_offset)%errno    =  iomaxerr + 266_I4B

   nfreya_offset = 267 + iomaxerr
    error(nfreya_offset)%descrip  = 'STRING_UTIL ERROR : sub set_r4_1darray_value'
    error(nfreya_offset)%attrib   = 'Error, see subroutine err 10'
    error(nfreya_offset)%errno    =  iomaxerr + 267_I4B

   nfreya_offset = 268 + iomaxerr
    error(nfreya_offset)%descrip  = 'STRING_UTIL ERROR : sub set_r4_1darray_value'
    error(nfreya_offset)%attrib   = 'Error,see subroutine error 11'
    error(nfreya_offset)%errno    =  iomaxerr + 268_I4B

   nfreya_offset = 269 + iomaxerr
    error(nfreya_offset)%descrip  = 'STRING_UTIL ERROR : sub set_r4_1darray_value'
    error(nfreya_offset)%attrib   = 'Error,see subroutine error 1'
    error(nfreya_offset)%errno    =  iomaxerr + 269_I4B

  nfreya_offset = 270 + iomaxerr
    error(nfreya_offset)%descrip  = 'STRING_UTIL ERROR : sub set_r4_1darray_value'
    error(nfreya_offset)%attrib   = 'Error,see subroutine error 2'
    error(nfreya_offset)%errno    =  iomaxerr + 270_I4B

  nfreya_offset = 271 + iomaxerr
    error(nfreya_offset)%descrip  = 'STRING_UTIL ERROR : sub set_r8_1darray_value a'
    error(nfreya_offset)%attrib   = 'Error,see subroutine error 1'
    error(nfreya_offset)%errno    =  iomaxerr + 271_I4B

  nfreya_offset = 272 + iomaxerr
    error(nfreya_offset)%descrip  = 'STRING_UTIL ERROR : sub set_i4_1darray_value'
    error(nfreya_offset)%attrib   = 'Error,see subroutine error 1'
    error(nfreya_offset)%errno    =  iomaxerr + 272_I4B

  nfreya_offset = 273 + iomaxerr
    error(nfreya_offset)%descrip  = 'STRING_UTIL ERROR : sub set_L_1darray_value'
    error(nfreya_offset)%attrib   = 'Error,see subroutine error 1'
    error(nfreya_offset)%errno    =  iomaxerr + 273_I4B

  nfreya_offset = 274 + iomaxerr
    error(nfreya_offset)%descrip  = 'STRING_UTIL ERROR : sub set_L_1darray_value'
    error(nfreya_offset)%attrib   = 'Error,see subroutine error 2'
    error(nfreya_offset)%errno    =  iomaxerr + 274_I4B

  nfreya_offset = 275 + iomaxerr
    error(nfreya_offset)%descrip  = 'STRING_UTIL ERROR : sub set_L_1darray_value'
    error(nfreya_offset)%attrib   = 'Error,parse error'
    error(nfreya_offset)%errno    =  iomaxerr + 275_I4B

  nfreya_offset = 276 + iomaxerr
    error(nfreya_offset)%descrip  = 'STRING_UTIL ERROR : sub set_r8_1darray_value'
    error(nfreya_offset)%attrib   = 'Error,see subroutine error 3'
    error(nfreya_offset)%errno    =  iomaxerr + 276_I4B

  nfreya_offset = 277 + iomaxerr
    error(nfreya_offset)%descrip  = 'STRING_UTIL ERROR : sub set_r8_1darray_value'
    error(nfreya_offset)%attrib   = 'Error,see subroutine error 4'
    error(nfreya_offset)%errno    =  iomaxerr + 277_I4B


  nfreya_offset = 278 + iomaxerr
    error(nfreya_offset)%descrip  = 'STRING_UTIL ERROR : sub set_r8_1darray_value'
    error(nfreya_offset)%attrib   = 'Error,see subroutine error 1'
    error(nfreya_offset)%errno    =  iomaxerr + 278_I4B

 nfreya_offset = 279 + iomaxerr
    error(nfreya_offset)%descrip  = 'STRING_UTIL ERROR : sub set_r8_1darray_value'
    error(nfreya_offset)%attrib   = 'Error,see subroutine error 2'
    error(nfreya_offset)%errno    =  iomaxerr + 279_I4B

 nfreya_offset = 280 + iomaxerr
    error(nfreya_offset)%descrip  = 'STRING_UTIL ERROR : sub set_r8_value'
    error(nfreya_offset)%attrib   = 'Error,see subroutine error 1'
    error(nfreya_offset)%errno    =  iomaxerr + 280_I4B

 nfreya_offset = 281 + iomaxerr
    error(nfreya_offset)%descrip  = 'STRING_UTIL ERROR : sub set_integer_value'
    error(nfreya_offset)%attrib   = 'Error,see subroutine error 1'
    error(nfreya_offset)%errno    =  iomaxerr + 281_I4B

 nfreya_offset = 282 + iomaxerr
    error(nfreya_offset)%descrip  = 'STRING_UTIL ERROR : sub readblock'
    error(nfreya_offset)%attrib   = 'Error,see subroutine error 1'
    error(nfreya_offset)%errno    =  iomaxerr + 282_I4B

 nfreya_offset = 283 + iomaxerr
    error(nfreya_offset)%descrip  = 'STRING_UTIL ERROR : sub readblock'
    error(nfreya_offset)%attrib   = 'Error,see subroutine error 2'
    error(nfreya_offset)%errno    =  iomaxerr + 283_I4B

 nfreya_offset = 284 + iomaxerr
    error(nfreya_offset)%descrip  = 'STRING_UTIL ERROR : sub readblock'
    error(nfreya_offset)%attrib   = 'Error,see subroutine error 3'
    error(nfreya_offset)%errno    =  iomaxerr + 284_I4B

 nfreya_offset = 285 + iomaxerr
    error(nfreya_offset)%descrip  = 'STRING_UTIL ERROR : sub readblock'
    error(nfreya_offset)%attrib   = 'Error,see subroutine error 4'
    error(nfreya_offset)%errno    =  iomaxerr + 285_I4B

 nfreya_offset = 286 + iomaxerr
    error(nfreya_offset)%descrip  = 'STRING_UTIL ERROR : sub readblock'
    error(nfreya_offset)%attrib   = 'Error,see subroutine error 6'
    error(nfreya_offset)%errno    =  iomaxerr + 286_I4B

 nfreya_offset = 287 + iomaxerr
    error(nfreya_offset)%descrip  = 'STRING_UTIL ERROR : sub readblock'
    error(nfreya_offset)%attrib   = 'Error,see subroutine error 7'
    error(nfreya_offset)%errno    =  iomaxerr + 287_I4B

 nfreya_offset = 288 + iomaxerr
    error(nfreya_offset)%descrip  = 'STRING_UTIL ERROR : sub read_scalar'
    error(nfreya_offset)%attrib   = 'Error,see subroutine error 1'
    error(nfreya_offset)%errno    =  iomaxerr + 288_I4B

 nfreya_offset = 289 + iomaxerr
    error(nfreya_offset)%descrip  = 'STRING_UTIL ERROR : sub read_scalar'
    error(nfreya_offset)%attrib   = 'Error,see subroutine error 2'
    error(nfreya_offset)%errno    =  iomaxerr + 289_I4B

nfreya_offset = 290 + iomaxerr
    error(nfreya_offset)%descrip  = 'STRING_UTIL ERROR : sub read_scalar'
    error(nfreya_offset)%attrib   = 'Error,see subroutine error 3'
    error(nfreya_offset)%errno    =  iomaxerr + 290_I4B

 nfreya_offset = 291 + iomaxerr
    error(nfreya_offset)%descrip  = 'STRING_UTIL ERROR : sub read_scalar1'
    error(nfreya_offset)%attrib   = 'Error,see subroutine error 1'
    error(nfreya_offset)%errno    =  iomaxerr + 291_I4B

 nfreya_offset = 292 + iomaxerr
    error(nfreya_offset)%descrip  = 'STRING_UTIL ERROR : sub read_scalar1'
    error(nfreya_offset)%attrib   = 'Error,see subroutine error 2'
    error(nfreya_offset)%errno    =  iomaxerr + 292_I4B

 nfreya_offset = 293 + iomaxerr
    error(nfreya_offset)%descrip  = 'STRING_UTIL ERROR : sub read_scalar1'
    error(nfreya_offset)%attrib   = 'Error,see subroutine error 3'
    error(nfreya_offset)%errno    =  iomaxerr + 293_I4B

 nfreya_offset = 294 + iomaxerr
    error(nfreya_offset)%descrip  = 'STRING_UTIL ERROR : sub read_scalar2'
    error(nfreya_offset)%attrib   = 'Error,see subroutine error 1'
    error(nfreya_offset)%errno    =  iomaxerr + 294_I4B

 nfreya_offset = 295 + iomaxerr
    error(nfreya_offset)%descrip  = 'STRING_UTIL ERROR : sub read_scalar2'
    error(nfreya_offset)%attrib   = 'Error,see subroutine error 2'
    error(nfreya_offset)%errno    =  iomaxerr + 295_I4B

 nfreya_offset = 296 + iomaxerr
    error(nfreya_offset)%descrip  = 'STRING_UTIL ERROR : sub convrt_12'
    error(nfreya_offset)%attrib   = 'Error,see subroutine error 1'
    error(nfreya_offset)%errno    =  iomaxerr + 296_I4B

 nfreya_offset = 297 + iomaxerr
    error(nfreya_offset)%descrip  = 'ERROR : sub get_nf_beam_data'
    error(nfreya_offset)%attrib   = 'Need ufile info because beam_data_namelist doesnt exist  but use_ufile = .FALSE. was specified'
    error(nfreya_offset)%errno    =  iomaxerr + 297_I4B


 nfreya_offset = 298 + iomaxerr
    error(nfreya_offset)%descrip  = 'ERROR : sub nfreya_load'
    error(nfreya_offset)%attrib   = 'Only one beam species is currently allowed'
    error(nfreya_offset)%errno    =  iomaxerr + 298_I4B

 nfreya_offset = 299 + iomaxerr
    error(nfreya_offset)%descrip  = 'ERROR : sub nfreya_load'
    error(nfreya_offset)%attrib   = 'Invalid neutral beam species'
    error(nfreya_offset)%errno    =  iomaxerr + 299_I4B

 nfreya_offset = 300 + iomaxerr
    error(nfreya_offset)%descrip  = 'ERROR : sub set_beam_data'
    error(nfreya_offset)%attrib   = 'Error, inconistency in dimensions'
    error(nfreya_offset)%errno    =  iomaxerr + 300_I4B

 nfreya_offset = 301 + iomaxerr
    error(nfreya_offset)%descrip  = 'ERROR : sub nfreya_load'
    error(nfreya_offset)%attrib   = 'Not set up to handle more than 2 tangency radii'
    error(nfreya_offset)%errno    =  iomaxerr + 301_I4B

 nfreya_offset = 302 + iomaxerr
    error(nfreya_offset)%descrip  = 'ERROR : sub nfreya_load'
    error(nfreya_offset)%attrib   = 'dt from nubeam for P_Nfreya not implemented'
    error(nfreya_offset)%errno    =  iomaxerr + 302_I4B

 nfreya_offset = 303 + iomaxerr
    error(nfreya_offset)%descrip  = 'ERROR : sub nfreya_load'
    error(nfreya_offset)%attrib   = 'unknown beam species'
    error(nfreya_offset)%errno    =  iomaxerr + 303_I4B


 nfreya_offset = 304 + iomaxerr
    error(nfreya_offset)%descrip  = 'ERROR : sub nfreya_load'
    error(nfreya_offset)%attrib   = 'number beams requested is too large'
    error(nfreya_offset)%errno    =  iomaxerr + 304_I4B


 nfreya_offset = 305 + iomaxerr
    error(nfreya_offset)%descrip  = 'ERROR : set_freya_beam'
    error(nfreya_offset)%attrib   = 'beam tangency radii  problem'
    error(nfreya_offset)%errno    =  iomaxerr + 305_I4B

 nfreya_offset = 306 + iomaxerr
    error(nfreya_offset)%descrip  = 'ERROR : sub nbdrive_nubeam_checknaml'
    error(nfreya_offset)%attrib   = 'duplicate species in thermal species list'
    error(nfreya_offset)%errno    =  iomaxerr + 306_I4B


 nfreya_offset = 307 + iomaxerr
    error(nfreya_offset)%descrip  = 'ERROR : sub nbdrive_nubeam_checknaml'
    error(nfreya_offset)%attrib   = 'duplicate species in impurity species list'
    error(nfreya_offset)%errno    =  iomaxerr + 307_I4B

 nfreya_offset = 308 + iomaxerr
    error(nfreya_offset)%descrip  = 'ERROR : sub nbdrive_nubeam_checknaml'
    error(nfreya_offset)%attrib   = 'nbi_datchkb detected problem with beam species'
    error(nfreya_offset)%errno    =  iomaxerr + 308_I4B

 nfreya_offset = 309 + iomaxerr
    error(nfreya_offset)%descrip  = 'ERROR : sub nbdrive_nubeam_checknaml'
    error(nfreya_offset)%attrib   = 'checknaml species count invalid'
    error(nfreya_offset)%errno    =  iomaxerr + 309_I4B

 nfreya_offset = 310 + iomaxerr
    error(nfreya_offset)%descrip  = 'ERROR : sub nbdrive_nubeam_checknaml'
    error(nfreya_offset)%attrib   = 'checknaml species type has a > 2'
    error(nfreya_offset)%errno    =  iomaxerr + 309_I4B



 nfreya_offset = 310 + iomaxerr
    error(nfreya_offset)%descrip  = 'ERROR : sub set_nubeam_thermal'
    error(nfreya_offset)%attrib   = 'number of thermal species not consistent'
    error(nfreya_offset)%errno    =  iomaxerr + 310_I4B

 nfreya_offset = 311 + iomaxerr
    error(nfreya_offset)%descrip  = 'ERROR : sub set_nubeam_impure'
    error(nfreya_offset)%attrib   = 'number of impurity  species not consistent'
    error(nfreya_offset)%errno    =  iomaxerr + 311_I4B

 nfreya_offset = 312 + iomaxerr
    error(nfreya_offset)%descrip  = 'ERROR : sub set_nubeam_impure'
    error(nfreya_offset)%attrib   = 'impurity species not reckognized'
    error(nfreya_offset)%errno    =  iomaxerr + 312_I4B


 nfreya_offset = 313 + iomaxerr
    error(nfreya_offset)%descrip  = 'ERROR : sub set_freya_beam'
    error(nfreya_offset)%attrib   = 'number beam speceis not set correctly'
    error(nfreya_offset)%errno    =  iomaxerr + 313_I4B

 nfreya_offset = 314 + iomaxerr
    error(nfreya_offset)%descrip  = 'ERROR : sub set_Pnfreya_beam'
    error(nfreya_offset)%attrib   = 'beam simulation times not set correctly in run directives file'
    error(nfreya_offset)%errno    =  iomaxerr + 314_I4B

    nfreya_offset = 315 + iomaxerr
    error(nfreya_offset)%descrip  = 'ERROR ONETWO : sub write_nfreya_run_directives '
    error(nfreya_offset)%attrib   = 'Onetwo reports error in creating P_nfreya namelist file'
    error(nfreya_offset)%errno    =  iomaxerr + 315_I4B

    nfreya_offset = 316 + iomaxerr
    error(nfreya_offset)%descrip  = 'ERROR P_nfreya_interface : get_nf_beam_data '
    error(nfreya_offset)%attrib   = 'nubeam namelist file not found'
    error(nfreya_offset)%errno    =  iomaxerr + 316_I4B

    nfreya_offset = 320 + iomaxerr
    error(nfreya_offset)%descrip  = 'ERROR P_nfreya_interface : sub UFILE_BYPASS '
    error(nfreya_offset)%attrib   = 'beam is not allowed to be on before time0'
    error(nfreya_offset)%errno    =  iomaxerr + 320_I4B

    nfreya_offset = 325 + iomaxerr
    error(nfreya_offset)%descrip  = 'ERROR file nfreya_load : sub set_Pnfreya_beam '
    error(nfreya_offset)%attrib   = 'name_beam%beam_times must be ge 2'
    error(nfreya_offset)%errno    =  iomaxerr + 325_I4B


    nfreya_offset = 350 + iomaxerr
    error(nfreya_offset)%descrip  = 'ERROR file P_Nfreya_rpc_interface.f90 : sub P_Nfreya_rpc_driver'
    error(nfreya_offset)%attrib   = 'P_nfreya did not reutrn expected statefile'
    error(nfreya_offset)%errno    =  iomaxerr + 350_I4B


    nfreya_offset = 351 + iomaxerr
    error(nfreya_offset)%descrip  = 'ERROR file P_Nfreya_rpc_interface.f90 : sub P_Nfreya_rpc_driver'
    error(nfreya_offset)%attrib   = 'unexpected end in "P_Nfreya_rpc_config"  file namelist read'
    error(nfreya_offset)%errno    =  iomaxerr + 351_I4B


    nfreya_offset = 352 + iomaxerr
    error(nfreya_offset)%descrip  = 'ERROR file P_Nfreya_rpc_interface.f90 : sub setup_run_P_Nfreya'
    error(nfreya_offset)%attrib   = 'Error, character string length inconsistency'
    error(nfreya_offset)%errno    =  iomaxerr + 352_I4B


    nfreya_offset = 353 + iomaxerr
    error(nfreya_offset)%descrip  = 'ERROR file P_Nfreya_rpc_interface.f90 : sub P_Nfreya_rpc_driver'
    error(nfreya_offset)%attrib   = 'Error, input configuration file "P_Nfreya_rpc_config" not found'
    error(nfreya_offset)%errno    =  iomaxerr + 353_I4B


    nfreya_offset = 354 + iomaxerr
    error(nfreya_offset)%descrip  = 'ERROR file P_Nfreya_rpc_interface.f90 : sub P_Nfreya_rpc'
    error(nfreya_offset)%attrib   = 'Error,  nubeam_namelist or nubeam_ufile not set correctly'
    error(nfreya_offset)%errno    =  iomaxerr + 354_I4B



    nfreya_offset = 355 + iomaxerr
    error(nfreya_offset)%descrip  = 'ERROR file thermalization.f90 : function cxrv_m'
    error(nfreya_offset)%attrib   = 'Error, argument x is out of range '
    error(nfreya_offset)%errno    =  iomaxerr + 355_I4B



    nfreya_offset = 356 + iomaxerr
    error(nfreya_offset)%descrip  = 'ERROR file nfreya_namelist.F90 : sub add_quotes_to_namelist'
    error(nfreya_offset)%attrib   = 'Error, processing run directives namelist'
    error(nfreya_offset)%errno    =  iomaxerr + 356_I4B

    nfreya_offset = 357 + iomaxerr
    error(nfreya_offset)%descrip  = 'ERROR file nfreya_routines.f90 : sub beam_prop'
    error(nfreya_offset)%attrib   = 'Error, did not get valid beam energies'
    error(nfreya_offset)%errno    =  iomaxerr + 357_I4B


    nfreya_offset = 358 + iomaxerr
    error(nfreya_offset)%descrip  = 'ERROR file rw_iterdb/netcdf.f90'
    error(nfreya_offset)%attrib   = 'Inconsistent dimension for bptor'
    error(nfreya_offset)%errno    =  iomaxerr + 358_I4B


    nfreya_offset = 359 + iomaxerr
    error(nfreya_offset)%descrip  = 'ERROR file nfreya_load.f90'
    error(nfreya_offset)%attrib   = 'Only beamlet numbers 1-14 are currently recognized'
    error(nfreya_offset)%errno    =  iomaxerr + 359_I4B



    nfreya_offset = 375 + iomaxerr
    error(nfreya_offset)%descrip  = 'ERROR file nfreya_load.f90'
    error(nfreya_offset)%attrib   = 'nmbr beamlets not consistent in Check_P_nfreya_beamlets'
    error(nfreya_offset)%errno    =  iomaxerr + 375_I4B

    nfreya_offset = 376 + iomaxerr
    error(nfreya_offset)%descrip  = 'ERROR file P_nfreya_interface.F90'
    error(nfreya_offset)%attrib   = 'namelist input NLBDAT turns off ufile read but namelistt doesnt contain powers,etc.'
    error(nfreya_offset)%errno    =  iomaxerr + 376_I4B

    nfreya_end = 400 + iomaxerr

    mhd_offset = nfreya_end + 1
    imhd = mhd_offset 
    error(imhd)%descrip = 'ERROR in sub gs_solver'
    error(imhd)%attrib  = '  invalid mhd_dat%equil_type input  in namelist'
    error(imhd)%errno   =  imhd

    imhd = mhd_offset+1
    error(imhd)%descrip = 'ERROR in sub gs_solver'
    error(imhd)%attrib  = '  invalid mhd_dat%equil_type input  in namelist'
    error(imhd)%errno   =  imhd

    imhd = mhd_offset+2
    error(imhd)%descrip = 'ERROR in sub cntour1'
    error(imhd)%attrib  = '  more points on plasma were found than dimension of xcontr,ycontr '
    error(imhd)%errno   =  imhd


    imhd = mhd_offset+3
    error(imhd)%descrip = 'ERROR detected in sub psi_contours'
    error(imhd)%attrib  = '  ierr non zero returned by sub stage_cntour_call'
    error(imhd)%errno   =  imhd

    imhd = mhd_offset+4
    error(imhd)%descrip = 'ERROR detected in sub psi_contours'
    error(imhd)%attrib  = ' rmajor inboard/outboard calcualtion failed'
    error(imhd)%errno   =  imhd


    imhd = mhd_offset+5
    error(imhd)%descrip = 'ERROR detected in sub psi_contours'
    error(imhd)%attrib  = 'value of psi traced was modified from input value but modificstion is not allowed '
    error(imhd)%errno   =  imhd


    imhd = mhd_offset+6
    error(imhd)%descrip = 'ERROR detected in sub stage_cntour_call'
    error(imhd)%attrib  = 'must select cntour1 or cntour2 with use_cnt'
    error(imhd)%errno   =  imhd

    imhd = mhd_offset+7
    error(imhd)%descrip = 'ERROR detected in sub stage_cntour_call'
    error(imhd)%attrib  = 'cntour1 or cntour2 returned  an error'
    error(imhd)%errno   =  imhd


    imhd = mhd_offset+8
    error(imhd)%descrip = 'ERROR detected in sub cntour1'
    error(imhd)%attrib  = 'eval2d spline error (not possible currently)'
    error(imhd)%errno   =  imhd


    imhd = mhd_offset+9
    error(imhd)%descrip = 'ERROR detected in sub cntour2'
    error(imhd)%attrib  = 'eval2d spline error (not possible currently)'
    error(imhd)%errno   =  imhd

    imhd = mhd_offset+10
    error(imhd)%descrip = 'ERROR detected in sub cntour1'
    error(imhd)%attrib  = 'must have max psi on axis'
    error(imhd)%errno   =  imhd

    imhd = mhd_offset+11
    error(imhd)%descrip = 'ERROR detected in sub cntour2'
    error(imhd)%attrib  = 'must have max psi on axis'
    error(imhd)%errno   =  imhd


    imhd = mhd_offset+12
    error(imhd)%descrip = 'ERROR detected in sub cntour1'
    error(imhd)%attrib  = 'unspecified erro'
    error(imhd)%errno   =  imhd

    imhd = mhd_offset+13
    error(imhd)%descrip = 'ERROR detected in sub cntour2'
    error(imhd)%attrib  = 'unspecified erro'
    error(imhd)%errno   =  imhd



    imhd = mhd_offset+14
    error(imhd)%descrip = 'ERROR detected in sub cntour2'
    error(imhd)%attrib  = 'newti problem,dxx .LE. dxmin'
    error(imhd)%errno   =  imhd

   END SUBROUTINE load_errno


  END MODULE error_handler
