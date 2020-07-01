program test

  use geo

  double precision, dimension(2) :: theta

  theta(:) = (/0.1,0.2/)
  
  call geo_interp(size(theta),theta,.true.)
  
  print *,geo_b
  
end program test
