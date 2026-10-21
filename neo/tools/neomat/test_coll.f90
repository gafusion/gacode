! Print emat_coll_test, emat_coll_field for 3 species (D, e, W), all pairs.
! Args: n_energy n_xi laguerre_method collision_model [coll_uncoupledei_model]
program test_coll
  use neo_globals
  use neo_energy_grid
  implicit none
  integer :: is, js, ie, je, ix
  character(len=64) :: arg
  call get_command_argument(1,arg); read(arg,*) n_energy
  call get_command_argument(2,arg); read(arg,*) n_xi
  call get_command_argument(3,arg); read(arg,*) laguerre_method
  call get_command_argument(4,arg); read(arg,*) collision_model
  n_species = 3
  spitzer_model = 0
  coll_uncoupledei_model = 0
  if (command_argument_count() > 4) then
     call get_command_argument(5,arg); read(arg,*) coll_uncoupledei_model
  endif
  coll_uncoupledaniso_model = 0
  write_cmoments_flag = 0
  ae_flag = 0
  is_ele = 2
  allocate(mass(3), z(3), dens(3,1), temp(3,1), vth(3,1), nu(3,1))
  mass = (/ 1.0, 2.724486e-4, 92.0 /)
  z    = (/ 1.0, -1.0, 30.0 /)
  dens(:,1) = (/ 1.0, 1.03, 1.0e-3 /)
  temp(:,1) = (/ 1.0, 1.2, 1.0 /)
  nu(:,1)   = (/ 1.0e-3, 6.0e-2, 1.0e-3 /)
  vth(:,1)  = sqrt(temp(:,1)/mass)
  call ENERGY_basis_ints_alloc(1)
  call ENERGY_coll_ints_alloc(1)
  call ENERGY_coll_ints(1)
  do is=1,3
     do js=1,3
        do ix=0,n_xi
           do ie=0,n_energy
              do je=0,n_energy
                 write(*,'(5i4,2es45.35)') is, js, ix, ie, je, emat_coll_test(is,js,ie,je,ix), &
                      emat_coll_field(is,js,ie,je,ix)
              enddo
           enddo
        enddo
     enddo
  enddo
end program test_coll
