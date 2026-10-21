! Print the ENERGY_basis_ints arrays.  Args: n_energy n_xi laguerre_method
program test_basis
  use neo_globals, only : n_energy, n_xi, laguerre_method
  use neo_energy_grid
  implicit none
  integer :: ie, je, ix, kx
  character(len=64) :: arg
  call get_command_argument(1,arg); read(arg,*) n_energy
  call get_command_argument(2,arg); read(arg,*) n_xi
  call get_command_argument(3,arg); read(arg,*) laguerre_method
  call ENERGY_basis_ints_alloc(1)
  call ENERGY_basis_ints
  do ix=0,n_xi
     do ie=0,n_energy
        write(*,'(a,3i4,5es45.35)') 'v', ix, ie, 0, evec_e0(ie,ix), evec_e05(ie,ix), &
             evec_e1(ie,ix), evec_e105(ie,ix), evec_e2(ie,ix)
        do je=0,n_energy
           write(*,'(a,3i4,es45.35)') 'e0', ix, ie, je, emat_e0(ie,je,ix,1)
           do kx=1,2
              write(*,'(a,4i4,3es45.35)') 'k2', ix, ie, je, kx, emat_e05(ie,je,ix,kx), &
                   emat_en05(ie,je,ix,kx), emat_e05de(ie,je,ix,kx)
           enddo
           do kx=1,3
              write(*,'(a,4i4,es45.35)') 'e1', ix, ie, je, kx, emat_e1(ie,je,ix,kx)
           enddo
        enddo
     enddo
  enddo
end program test_basis
