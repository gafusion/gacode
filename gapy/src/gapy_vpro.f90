module vpro

  integer, parameter :: ntag = 38
  integer :: nt

  ! Fundamental input
  
  character*10, dimension(ntag) :: tag = (/&
       'nexp      ',& !1
       'nion      ',& !2
       'bt_exp    ',& !3
       'arho_exp  ',& !4
       'rho       ',& !5
       'rmin      ',& !6
       'polflux   ',& !7
       'q         ',& !8
       'omega0    ',& !9
       'rmaj      ',& !10
       'zmag      ',& !11
       'kappa     ',& !12
       'delta     ',& !13
       'zeta      ',& !14
       'ne        ',& !15
       'Te        ',& !16
       'ptot      ',& !17
       'z_eff     ',& !18
       'flow_beam ',& !19
       'flow_wall ',& !20
       'flow_mom  ',& !21
       'pow_e     ',& !22
       'pow_i     ',& !23
       'pow_ei    ',& !24
       'pow_e_aux ',& !25
       'pow_i_aux ',& !26
       'pow_e_fus ',& !27
       'pow_i_fus ',& !28
       'pow_e_sync',& !29
       'pow_e_brem',& !30
       'pow_e_line',& !31
       'sbeame    ',& !32
       'sbcx      ',& !33
       'sscxl     ',& !34
       'ni        ',& !35
       'ti        ',& !36
       'vpol      ',& !37
       'vtor      '&  !38
       /)

  integer :: &
       nexp,&
       nion

  double precision :: &
       bt_exp,&
       arho_exp

  double precision, dimension(:), allocatable :: &
       rho,&
       rmin,&
       polflux,&
       q,&
       omega0,&
       rmaj,&
       zmag,&
       kappa,&
       delta,&
       zeta,&
       ne,&
       te,&
       ptot,&
       z_eff,&
       flow_beam,&
       flow_wall,&
       flow_mom,&
       pow_e,&
       pow_i,&
       pow_ei,&
       pow_e_aux,&
       pow_i_aux,&
       pow_e_fus,&
       pow_i_fus,&
       pow_e_sync,&
       pow_e_brem,&
       pow_e_line,&
       sbeame,&
       sbcx,&
       sscxl

  double precision, dimension(:,:), allocatable :: &
       ni,&
       ti,&
       vpol,&
       vtor

  ! Derived
  
  double precision, dimension(:), allocatable :: &
       bunit,&
       s
  
contains

  subroutine vpro_read

    implicit none

    integer :: i
    character(len=10) :: ytag

    ! NOTE: nexp should appear before any profile arrays

    open(unit=1,file='dat.bin',status='old',access='stream')
    read(1) nt
    
    do i=1,nt

       read(1) ytag

       select case (trim(ytag))
       case('nexp')
          read(1) nexp
          call vpro_init(1)
       case('nion')
          read(1) nion
       case('bt_exp')
          read(1) bt_exp
       case('arho_exp')
          read(1) arho_exp
       case ('rho')
          read(1) rho 
       case ('rmin')
          read(1) rmin 
       case ('polflux')
          read(1) polflux 
       case ('q')
          read(1) q 
       case ('omega0')
          read(1) omega0 
       case ('rmaj')
          read(1) rmaj 
       case ('zmag')
          read(1) zmag 
       case ('kappa')
          read(1) kappa 
       case ('delta')
          read(1) delta 
       case ('zeta')
          read(1) zeta 
       case ('ne')
          read(1) ne
       case ('Te')
          read(1) te 
       case ('ptot')
          read(1) ptot 
       case ('z_eff')
          read(1) z_eff 
       case ('ni')
          read(1) ni(:,1:nion)
       case ('ti')
          read(1) ti(:,1:nion)
       case ('vpol')
          read(1) vpol(:,1:nion)
       case ('vtor')
          read(1) vtor(:,1:nion)
       case ('flow_beam')
          read(1) flow_beam 
       case ('flow_wall')
          read(1) flow_wall 
       case ('flow_mom')
          read(1) flow_mom 
       case ('pow_e')
          read(1) pow_e 
       case ('pow_i')
          read(1) pow_i 
       case ('pow_ei')
          read(1) pow_ei 
       case ('pow_e_aux')
          read(1) pow_e_aux 
       case ('pow_i_aux')
          read(1) pow_i_aux 
       case ('pow_e_fus')
          read(1) pow_e_fus 
       case ('pow_i_fus')
          read(1) pow_i_fus 
       case ('pow_e_sync')
          read(1) pow_e_sync
       case ('pow_e_brem')
          read(1) pow_e_brem 
       case ('pow_e_line')
          read(1) pow_e_line 
       case ('sbeame')
          read(1) sbeame 
       case ('sbcx')
          read(1) sbcx 
       case ('sscxl')
          read(1) sscxl 
       end select
       
    enddo

  end subroutine vpro_read

  subroutine vpro_init(flag)

    implicit none

    integer, intent(in) :: flag

    if (flag == 1) then
       
       allocate(rho(nexp))
       allocate(rmin(nexp))
       allocate(q(nexp))
       allocate(polflux(nexp))
       allocate(omega0(nexp))
       allocate(rmaj(nexp))
       allocate(zmag(nexp))
       allocate(kappa(nexp))
       allocate(delta(nexp))
       allocate(zeta(nexp))
       allocate(ne(nexp))
       allocate(te(nexp))
       allocate(ptot(nexp))
       allocate(z_eff(nexp))

       allocate(flow_beam(nexp))
       allocate(flow_wall(nexp))
       allocate(flow_mom(nexp))
       allocate(pow_e(nexp))
       allocate(pow_i(nexp))
       allocate(pow_ei(nexp))
       allocate(pow_e_aux(nexp))
       allocate(pow_i_aux(nexp))
       allocate(pow_e_fus(nexp))
       allocate(pow_i_fus(nexp))
       allocate(pow_e_sync(nexp))
       allocate(pow_e_brem(nexp))
       allocate(pow_e_line(nexp))
       allocate(sbeame(nexp))
       allocate(sbcx(nexp))
       allocate(sscxl(nexp))

       rho     = 0.0
       rmin    = 0.0
       q       = 0.0
       polflux = 0.0
       omega0  = 0.0
       rmaj    = 0.0
       zmag    = 0.0
       kappa   = 0.0
       delta   = 0.0
       zeta    = 0.0
       ne      = 0.0
       te      = 0.0
       ptot    = 0.0
       z_eff   = 0.0

       flow_beam = 0.0
       flow_wall = 0.0
       flow_mom  = 0.0
       pow_e     = 0.0
       pow_i     = 0.0
       pow_ei    = 0.0
       pow_e_aux = 0.0
       pow_i_aux = 0.0
       pow_e_fus = 0.0
       pow_i_fus = 0.0
       pow_e_sync = 0.0
       pow_e_brem = 0.0
       pow_e_line = 0.0
       sbeame = 0.0
       sbcx   = 0.0
       sscxl  = 0.0
       
       allocate(ni(nexp,10))
       allocate(ti(nexp,10))
       allocate(vpol(nexp,10))
       allocate(vtor(nexp,10))

       ni = 0.0
       ti = 0.0
       vpol = 0.0
       vtor = 0.0

       ! Derived
       
       allocate(bunit(nexp))
       allocate(s(nexp))

       bunit = 0.0
       s      = 0.0
       
    else

       deallocate(rho)
       deallocate(rmin)
       deallocate(q)
       deallocate(polflux)
       deallocate(omega0)
       deallocate(rmaj)
       deallocate(zmag)
       deallocate(kappa)
       deallocate(delta)
       deallocate(zeta)
       deallocate(ne)
       deallocate(te)
       deallocate(ptot)
       deallocate(z_eff)

       deallocate(flow_beam)
       deallocate(flow_wall)
       deallocate(flow_mom)
       deallocate(pow_e)
       deallocate(pow_i)
       deallocate(pow_ei)
       deallocate(pow_e_aux)
       deallocate(pow_i_aux)
       deallocate(pow_e_fus)
       deallocate(pow_i_fus)
       deallocate(pow_e_sync)
       deallocate(pow_e_brem)
       deallocate(pow_e_line)
       deallocate(sbeame)
       deallocate(sbcx)
       deallocate(sscxl)
      
       deallocate(ni,ti,vpol,vtor)

       ! Derived
       
       deallocate(bunit,s)
       
    endif

  end subroutine vpro_init

  subroutine vpro_write

    implicit none

    open(unit=1,file='dat.bin',status='replace',access='stream')
    write(1) ntag

    write(1) tag(1)
    write(1) nexp

    write(1) tag(2)
    write(1) nion

    write(1) tag(3)
    write(1) bt_exp

    write(1) tag(4)
    write(1) arho_exp
    
    write(1) tag(5)
    write(1) rho

    write(1) tag(6)
    write(1) rmin

    write(1) tag(7)
    write(1) polflux

    write(1) tag(8)
    write(1) q

    write(1) tag(9)
    write(1) omega0

    write(1) tag(10)
    write(1) rmaj

    write(1) tag(11)
    write(1) zmag

    write(1) tag(12)
    write(1) kappa

    write(1) tag(13)
    write(1) delta

    write(1) tag(14)
    write(1) zeta

    write(1) tag(15)
    write(1) ne

    write(1) tag(16)
    write(1) te

    write(1) tag(17)
    write(1) ptot

    write(1) tag(18)
    write(1) z_eff

    write(1) tag(19)
    write(1) flow_beam

    write(1) tag(20)
    write(1) flow_wall

    write(1) tag(21)
    write(1) flow_mom

    write(1) tag(22)
    write(1) pow_e

    write(1) tag(23)
    write(1) pow_i

    write(1) tag(24)
    write(1) pow_ei

    write(1) tag(25)
    write(1) pow_e_aux

    write(1) tag(26)
    write(1) pow_i_aux

    write(1) tag(27)
    write(1) pow_e_fus

    write(1) tag(28)
    write(1) pow_i_fus

    write(1) tag(29)
    write(1) pow_e_sync

    write(1) tag(30)
    write(1) pow_e_brem

    write(1) tag(31)
    write(1) pow_e_line

    write(1) tag(32)
    write(1) sbeame

    write(1) tag(33)
    write(1) sbcx

    write(1) tag(34)
    write(1) sscxl

    write(1) tag(35)
    write(1) ni(:,1:nion)

    write(1) tag(36)
    write(1) ti(:,1:nion)

    write(1) tag(37)
    write(1) vpol(:,1:nion)

    write(1) tag(38)
    write(1) vtor(:,1:nion)

    close(1)

  end subroutine vpro_write

  subroutine vpro_read_legacy

    implicit none

    integer :: i
    character(len=99) :: line
    double precision :: x(5)

    open(unit=1,file='input.profiles',status='old')
    do while (line(1:2) /= '#r')
       read(1,'(a)') line
       if (line(1:5) == 'N_EXP') then
          read(line(7:),*) nexp
       endif
       if (line(1:5) == 'N_ION') then
          read(line(7:),*) nion
       endif
       if (line(1:6) == 'BT_EXP') then
          read(line(8:),*) bt_exp
       endif
       if (line(1:8) == 'ARHO_EXP') then
          read(line(10:),*) arho_exp
       endif
    enddo

    call vpro_init(1)

    ! 1
    do i=1,nexp
       read(1,*) x
       rho(i)     = x(1)
       rmin(i)    = x(2)
       polflux(i) = x(3)
       q(i)       = x(4)
       omega0(i)  = x(5)
    enddo

    read(1,'(a)') line
    read(1,'(a)') line

    ! 2
    do i=1,nexp
       read(1,*) x
       rmaj(i)  = x(1)
       zmag(i)  = x(2)
       kappa(i) = x(3)
       delta(i) = x(4)
       zeta(i)  = x(5)
    enddo

    read(1,'(a)') line
    read(1,'(a)') line

    ! 3
    do i=1,nexp
       read(1,*) x
       ne(i)    = x(1)
       te(i)    = x(2)
       ptot(i)  = x(3)
       z_eff(i) = x(4)
    enddo

    read(1,'(a)') line
    read(1,'(a)') line

    ! 4
    do i=1,nexp
       read(1,*) x
       ni(i,1:nion) = x(1:nion)
    enddo

    read(1,'(a)') line
    read(1,'(a)') line

    ! 5
    do i=1,nexp
       read(1,*) x
    enddo

    read(1,'(a)') line
    read(1,'(a)') line

    ! 6
    do i=1,nexp
       read(1,*) x
       ti(i,1:nion) = x(1:nion)
    enddo

    read(1,'(a)') line
    read(1,'(a)') line

    ! 7
    do i=1,nexp
       read(1,*) x
    enddo
    close(1)

  end subroutine vpro_read_legacy

  subroutine vpro_compute_derived

    use util
    
    implicit none

    double precision :: rhod(nexp)
    
    rhod(:) = arho_exp*rho(:)

    ! b_unit
    call util_bound_deriv(bunit,rhod**2,rmin**2,nexp)
    bunit(:) = bt_exp*bunit

    ! s
    call util_bound_deriv(s,q,rmin,nexp)
    s(:) = (rmin(:)/q(:))*s(:)

  end subroutine vpro_compute_derived
  
end module vpro
