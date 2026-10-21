! Replicates NEO ENERGY_coll_ints, collision_model=4, for one species pair
! (tauinv_ab=1).  Args: n_energy laguerre_method n_xi m_a/m_b T_a/T_b
program drv
  implicit none
  integer, parameter :: e_alpha=2
  integer :: n_energy, n_xi, meth, ix, jx, ie, je, ke, me, xarg, yarg, fmarg, i
  integer, allocatable :: e_lag(:), xb(:)
  real, allocatable :: g2(:), fc(:,:), fcb(:,:), e0(:,:), ct(:,:), cf(:,:)
  real :: pi, lam, mr, tr, z0, z1, z2, rd, r1, r2, r3, r4, r5, r6, r7
  character(len=64) :: arg
  call get_command_argument(1,arg); read(arg,*) n_energy
  call get_command_argument(2,arg); read(arg,*) meth
  call get_command_argument(3,arg); read(arg,*) n_xi
  call get_command_argument(4,arg); read(arg,*) mr
  call get_command_argument(5,arg); read(arg,*) tr
  pi = 4.0*atan(1.0); lam = tr/mr
  allocate(e_lag(0:n_xi), xb(0:n_xi))
  if (meth == 1) then
     e_lag = 3; e_lag(0) = 1; xb = 1; xb(0) = 0
  else
     do i=0,n_xi; e_lag(i) = 2*i+1; xb(i) = i; enddo
  endif
  allocate(g2(1:4*e_alpha*n_energy+4*n_xi+12))
  do i=1,size(g2); g2(i) = gamma2(i); enddo
  fmarg = 2*e_alpha*n_energy + 2*n_xi + 4
  allocate(fc(-fmarg:fmarg,-fmarg:fmarg), fcb(-fmarg:fmarg,-fmarg:fmarg))
  call neo_compute_fcoll(fmarg,lam,fc,fcb)
  allocate(e0(0:n_energy,0:n_energy), ct(0:n_energy,0:n_energy), cf(0:n_energy,0:n_energy))
  do ix=0,n_xi
     jx = ix
     e0 = 0.0; ct = 0.0; cf = 0.0
     do ie=0,n_energy
        do je=0,n_energy
           do ke=0,ie
              do me=0,je
                 z0 = (-1.0)**(ke+me)
                 z1 = (g2(2+2*ie+e_lag(ix))/g2(2+2*(ie-ke))/g2(2+2*ke+e_lag(ix))) &
                      *(g2(2+2*je+e_lag(jx))/g2(2+2*(je-me))/g2(2+2*me+e_lag(jx)))
                 z2 = g2(2+2*ke)*g2(2+2*me)
                 xarg = e_alpha*(ke+me) + xb(ix) + xb(jx) + 3
                 e0(ie,je) = e0(ie,je) + 0.5*z0*z1*(g2(xarg)/z2)
                 xarg = e_alpha*(ke+me) + xb(ix) + xb(jx)
                 ct(ie,je) = ct(ie,je) - sqrt(lam/pi)*ix*(ix+1)*z0*z1 &
                      *(fc(xarg-1,0)-fc(xarg-3,2))/z2
                 rd = 4.0*sqrt(lam/pi)*(e_alpha*ke+xb(ix))*z0*z1 &
                      *((1.0-tr)*fc(xarg-1,2) - 0.5*(e_alpha*me+xb(jx))*fc(xarg-3,2))/z2
                 ct(ie,je) = ct(ie,je) + rd
                 r1 = mr*(1.0+lam)**(-0.5*(xarg+3))*g2(xarg+3)
                 xarg = e_alpha*ke + xb(ix); yarg = e_alpha*me + xb(jx)
                 r2 = -2.0/(ix+0.5)*(mr-ix*(1.0-mr))*fc(xarg-ix+1,yarg+jx+2)
                 r3 = -2.0/(ix+0.5)*(1.0+ix*(1.0-mr))*fcb(xarg+ix+2,yarg-jx+1)
                 r4 = -ix*(ix-1.0)/(ix*ix-0.25)*fc(xarg-ix+3,yarg+jx+2)
                 r5 = -ix*(ix-1.0)/(ix*ix-0.25)*fcb(xarg+ix+2,yarg-jx+3)
                 r6 = (ix+1.0)*(ix+2.0)/(ix+1.5)/(ix+0.5)*fc(xarg-ix+1,yarg+jx+4)
                 r7 = (ix+1.0)*(ix+2.0)/(ix+1.5)/(ix+0.5)*fcb(xarg+ix+4,yarg-jx+1)
                 cf(ie,je) = cf(ie,je) + 2.0/sqrt(pi)*lam**1.5*lam**(0.5*(e_alpha*me+xb(ix))) &
                      *z0*z1*(r1+r2+r3+r4+r5+r6+r7)/z2
              enddo
           enddo
        enddo
     enddo
     do ie=0,n_energy
        do je=0,n_energy
           write(*,'(3i4,3es45.35)') ix, ie, je, e0(ie,je), ct(ie,je), cf(ie,je)
        enddo
     enddo
  enddo
contains
  real function gamma2(n)
    integer, intent(in) :: n
    integer :: i, i1
    if (mod(n,2) == 0) then
       gamma2 = 1.0; i1 = 4
    else
       gamma2 = sqrt(pi); i1 = 3
    endif
    do i=i1,n,2
       gamma2 = gamma2*(i/2.0-1.0)
    enddo
  end function gamma2
end program drv
