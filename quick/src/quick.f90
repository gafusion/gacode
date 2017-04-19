program quick

  use quick_globals
  use EXPRO_interface

  implicit none

  ! Read QUICK input
  open(unit=1,file='input.quick.gen',status='old')
  read(1,*) n_ion
  read(1,*) z(1)
  read(1,*) z(2)
  read(1,*) z(3)
  close(1)
  
  ! Need to tell EXPRO about the number of ions and charges you want (just like in TGYRO).
  EXPRO_ctrl_n_ion = n_ion
  EXPRO_ctrl_quasineutral_flag = 0
  EXPRO_ctrl_z(1:n_ion) = z(1:n_ion)
  EXPRO_ctrl_numeq_flag = 0

  ! Read input.profiles.gen (i.e., call EXPRO)
  call EXPRO_alloc('./',1)
  call EXPRO_read

  ! Change some data
  EXPRO_ne(:)   = 1.5*EXPRO_ne(:)
  EXPRO_rmin(:) = 0.9*EXPRO_rmin(:)
  
  ! Write new input.profiles
  call EXPRO_write_original(1,'input.profiles',2,'input.profiles.new',&
       'Profiles modified by QUICK')

end program quick
