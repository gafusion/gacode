program test

  use expro
  use expro_locsim_interface

  implicit none

  real :: ipccw,btccw,a_meters

  call expro_read('../data/input.gacode')

  print '(2a)'              ,'name1    ',expro_name(1)
  print '(2a)'              ,'name2    ',expro_name(2)
  print '(a,i0)'            ,'nexp     ',expro_n_exp
  print '(a,1pe12.5)'       ,'torfluxa ',expro_torfluxa
  print '(a,60(1pe12.5,1x))','ni       ',expro_ni(1,:)
  print '(a,60(1pe12.5,1x))','bunit    ',expro_bunit

end program test
