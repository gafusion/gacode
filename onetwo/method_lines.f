

      subroutine addem2(no_neq,t,y_ode,ml,mu,pa,nrowpa)
c-----------------------------------------------------------------------------
c     routine called by dlsoid
c     given values of the dependent variable in y
c     calculate matrix em and add it to pa,pb,pc
c------------------------------------------------------------------HSJ

      USE param
      USE aid_newton
      USE m_lines
      USE     ions
      USE solcon
      USE soln
      USE io
      USE extra
      USE numbrs
      USE mesh
      USE verbose
      USE sourc
      USE machin
      USE geom
      USE flags
      USE tordlrot
      USE soln2d
      USE bd_condtn
      USE rhog
      USE nonlin
      USE flx
      USE ode_pkg
c      USE  io,                        ONLY : ioftn, iopntr, nunits

      implicit  integer (i-n), real*8 (a-h, o-z)


      real *8 t,y_ode(*),pa(nrowpa,*)
      real *8 em(kk,kk)


c      matrix em,set in get_em (and in abcg)
c      uses te,ti, rpb,e tc (NOT VECTOR U(K,J))
c      hence in the following we do not update u:
      do j=1,no_neq
         do i = 1,nprim
              if(nameu_gp(j)(1:2) .eq.  nameu(i))then
                 !primary ion densitites
              !en(grid_point(j),?) implement
           CALL STOP('subroutine addem:density not implemented', 6)
              endif
         enddo
         do i =1 ,nimp
            if(nameu_gp(j)(1:2) .eq.  nameu(nprim+i))then
               !impurity densities
               !eni(grid_point(j),?) implement
           CALL STOP('subroutine addem:density not implemented', 6)
            endif
         enddo 
         if(nameu_gp(j)(1:2) .eq. nameu(nprim+nimp+1))then          !te
            te(grid_point(j)) =  (y_ode(j)*typx(j))**2
         else if(nameu_gp(j)(1:2) .eq. nameu(nprim+nimp+2))then     !ti
            ti(grid_point(j)) =  (y_ode(j)*typx(j))**2 
         else if(nameu_gp(j)(1:2) .eq. nameu(nprim+nimp+3))then     !rbp
            rbp(grid_point(j)) =  (y_ode(j)*typx(j))**2
            if(grid_point(j) .eq. 1)rbp(j) =0.0d0
         else if(nameu_gp(j)(1:2) .eq. nameu(nprim+nimp+4))then     !wr
            angrot(grid_point(j)) =   y_ode(j)*typx(j)
         else
            CALL STOP('subroutine addem2:nameu_gp error', 6)
         endif
       enddo
c
c     

      print *,'addem2 ,ml,mu,nrowpa =',ml,mu,nrowpa
      !note that ml=mu = 2*nkt -1 are the upper and lower
      !bandwidth for the jacobian. (the total bandwith is ml+mu+1)
      !but em has lower and upper bandwith  = nkt-1 (total bamdwith of 2*nkt-1 
       
      do k = 1,nb
         call get_em(k,em)
         if (nkt .lt. nk)
     .    call reduce_em (em, itran, nk, kk)
           print *,'done reduce,k =',k
         call STOP('not finished coding',0)
c         do j = 1,mbb 
c            do i = 1,mbb
               if(k .eq. 1)then
                  print *,'em(i,j)',em(i,j)
               endif
                pa(k,k) = pa(k,k) + em(i,j)*2.*SQRT(te(k))*
     .                                    typx(k)/eqtn_scale(k)
c            enddo
c          enddo
      enddo
      j=0
      i=0
      do nb=1,te_index-1        !nb =1 to number of blocks of size nkt need rbp the same here
         j=j+1                           !grid point index
         call get_em(j,em)
         if (nkt .lt. nk)
     .      call reduce_em (em, itran, nk, kk)


         do nr = 1,nkt                  !loop over rows of block of size nkt
            i=i+1                       !counts number of equations
            jl = i-nr +1
            do nc = jl,jl+nkt-1
               j=i-nc+1
            print *,'i,j=',i,j
            enddo
         enddo
         print *,'nb,j=',nb,j
      enddo
      print *,' no_neq,i = ',no_neq,i  
      print *,'leaving addem2'
      call exit(1)
      return
      end


      subroutine addem(no_neq,t,y_ode,mb,nb,pa,pb,pc)
c-----------------------------------------------------------------------------
c     routine called by dlsoibt
c     given values of the dependent variable in y
c     calculate matrix em and add it to pa,pb,pc
c------------------------------------------------------------------HSJ

      USE param
      USE aid_newton
      USE m_lines
      USE     ions
      USE solcon
      USE soln
      USE io
      USE extra
      USE numbrs
      USE mesh
      USE verbose
      USE sourc
      USE machin
      USE geom
      USE flags
      USE tordlrot
      USE soln2d
      USE bd_condtn
      USE rhog
      USE nonlin
      USE flx
      USE ode_pkg
c      USE  io,                        ONLY : ioftn, iopntr, nunits
      implicit  integer (i-n), real*8 (a-h, o-z)


      real *8 t,y_ode(*),pa(mb,mb,nb), pb(mb,mb,nb), pc(mb,mb,nb)
      real *8 em(kk,kk)


c      matrix em,set in get_em (and in abcg)
c      uses te,ti, rpb,e tc (NOT VECTOR U(K,J))
c      hence in the following we do not update u:
      do j=1,no_neq
         do i = 1,nprim
              if(nameu_gp(j)(1:2) .eq.  nameu(i))then
                 !primary ion densitites
                 !en(grid_point(j),?) implement
           CALL STOP('subroutine addem:density not implemented', 6)
              endif
         enddo
         do i =1 ,nimp
            if(nameu_gp(j)(1:2) .eq.  nameu(nprim+i))then
               !impurity densities
               !eni(grid_point(j),?) implement
           CALL STOP('subroutine addem:density not implemented', 6)
            endif
         enddo 
         if(nameu_gp(j)(1:2) .eq. nameu(nprim+nimp+1))then          !te
            te(grid_point(j)) =  (y_ode(j)*typx(j))**2
            u(nprim+nimp+1,grid_point(j))=te(grid_point(j)) 
         else if(nameu_gp(j)(1:2) .eq. nameu(nprim+nimp+2))then     !ti
            ti(grid_point(j)) =  (y_ode(j)*typx(j))**2 
            u(nprim+nimp+2,grid_point(j))=ti(grid_point(j))
         else if(nameu_gp(j)(1:2) .eq. nameu(nprim+nimp+3))then     !rbp
            rbp(grid_point(j)) =  (y_ode(j)*typx(j))**2
            if(grid_point(j) .eq. 1)rbp(grid_point(j)) =0.0d0
            u(nprim+nimp+3,grid_point(j))=rbp(grid_point(j))
         else if(nameu_gp(j)(1:2) .eq. nameu(nprim+nimp+4))then     !wr
            angrot(grid_point(j)) =   y_ode(j)*typx(j)
            u(nprim+nimp+4,grid_point(j))=angrot(grid_point(j))
         else
            CALL STOP('subroutine addem:nameu_gp error', 6)
         endif
       enddo
c
c
      j=1
      print *,'addem ,nb,mb =',nb,mb
      do k = 1,nb
         call get_em(k,em)  !each block k is associated with a grid point
         if (nkt .lt. nk)
     .    call reduce_em (em, itran, nk, kk)

         do j = 1,mb
            do i = 1,mb
                 ieqrd =(k-1)*nkt + i
                 dumy = u(itranrd(i),k)  !block and gridpoint index the same
                if(k.eq.1 .and. itranrd(i) .eq. nprim+nimp+3
     .            .and. itranrd(j) .eq. nprim+nimp+3 )dumy = 1.d0
                if(itranrd(i) .eq. nprim+nimp+4 .and. itranrd(j)
     .             .eq. nprim+nimp+4)then
                   pa(i,j,k) = pa(i,j,k) + em(i,j)*dumy*
     .                         typx(ieqrd)/eqtn_scale(ieqrd)
                else
                   pa(i,j,k) = pa(i,j,k) + em(i,j)*2.*SQRT(dumy)*
     .                         typx(ieqrd)/eqtn_scale(ieqrd)
                endif
c               pb(i,j,k) = pb(i,j,k) +?
c               pc(i,j,k) = pc(i,j,k) +?
            enddo
          enddo
      enddo


      return
      end




      subroutine eval_ode2 (no_eq, ode_time, y_ode, ode_rhs,
     .                      dumy1,idumy1)
c
c -------------------------------------------------------------HSJ-----
c - This subroutine  is called when the (longitudinal) method of lines
c - is used to solve the diffusion equations.
c - The rhs of
c        M* dy/dt = -(1/hr)*(hr(q+c))+s
c   is returned in ode_rhs
c   the structure of y_ode and ode_rhs is determined
c   by subroutine get_vars. 
c   The assumed sequence is always nprim then nimp then te then ti then rbp
c   then rotation. This sequence must be strictly adhered to(variables that
c   are not run in simulation mode in this sequqnce are simply skipped over)
c   For example suppose we have
c   nkt=3 dependent variables,Te,Ti and rbp.
c   Since rbp at grid point  is laways zero it is skipped.
c     y_ode(1) = sqrt(te(1))/typx(1)
c     y_ode(2) = sqrt(ti(1))/typx(2)
c     y_ode(4) = sqrt(te(2))/typx(3)
c     y_ode(5) = sqrt(ti(2))/typx(4)
c     y_ode(6) = sqrt(rbp(2))/typx(5)
c     y_ode(7) = sqrt(te(3))/typx(6)
c     y_ode(8) = sqrt(ti(3))/typx(7)
c     y_ode(9) = sqrt(rbp(3))/typx(8)
c      and so on all the way to the edge grid points
c      Note that the edge grid point for each variable can be different and
c      the structure in y_ode reflects this.
c     note that for rotation we have
c     y_ode = angrot/typx rather than Sqrt(angrot)/typx since rotation
c              can be negative

c
c  INPUT
c  no_eq         # of equations to be solved 
c  ode_time      the time at which ode_rhs is to be evaluated
c  y_ode         the dependent variable vector at time t
c  from common (aid_newton):
c  typx()        scale factors
c  grid_point()  grid point number
c
c  OUTPUT
c  ode_rhs(no_eq)     this is the rhs of  M*dy/dt = rhs 

c -----------------------------------------------------------------HSJ-----
c
      USE param
      USE aid_newton
      USE ions
      USE io
      USE solcon
      USE soln
      USE extra
      USE numbrs
      USE mesh
      USE verbose
      USE sourc
      USE machin
      USE geom
      USE flags
      USE tordlrot
      USE soln2d
      USE bd_condtn
      USE rhog
      USE nonlin
      USE flx
      USE ode_pkg
c      USE  io,                        ONLY : ioftn, iopntr, nunits

      implicit  integer (i-n), real*8 (a-h, o-z)
c
      dimension y_ode(*), ode_rhs(*)

c
      dimension  y_odeb(kk)
c

      time         = ode_time
      write (6, '("time_ode (radau5) =", 1pe16.8)') ode_time
      print *,'eval_ode2,(radau5), te(1),ti(1) =',te(1),ti(1),iodeflg1

c
c     given the dependent variable values that are run in simulation
c     mode in y_ode load the usual vectors used in SOURCE, etc.:
c
      do j=1,no_eq
         do i = 1,nprim
              if(nameu_gp(j)(1:2) .eq.  nameu(i))then
                 !primary ion densitites
              !en(grid_point(j),?) implement
           CALL STOP('subroutine eval_ode2:density not implemented', 6)
              endif
         enddo
         do i =1 ,nimp
            if(nameu_gp(j)(1:2) .eq.  nameu(nprim+i))then
               !impurity densities
               !eni(grid_point(j),?) implement
           CALL STOP('subroutine eval_ode2:density not implemented', 6)
            endif
         enddo 
         if(nameu_gp(j)(1:2) .eq. nameu(nprim+nimp+1))then          !te
            te(grid_point(j)) =  (y_ode(j)*typx(j))**2
c            print *,'te ,j =',te(grid_point(j)),j,grid_point(j)
            u(nprim+nimp+1,grid_point(j))=(y_ode(j)*typx(j))**2
         else if(nameu_gp(j)(1:2) .eq. nameu(nprim+nimp+2))then     !ti
            ti(grid_point(j)) =  (y_ode(j)*typx(j))**2 
            u(nprim+nimp+2,grid_point(j))=(y_ode(j)*typx(j))**2
c            print *,'ti ,j =',ti(grid_point(j)),j,grid_point(j)
         else if(nameu_gp(j)(1:2) .eq. nameu(nprim+nimp+3))then     !rbp
c            if(grid_point(j) .eq. 1)
c     .      print *,'rbp,j',rbp(grid_point(j)),grid_point(j)
            rbp(grid_point(j)) =  (y_ode(j)*typx(j))**2

            if(grid_point(j) .eq. 1)rbp(grid_point(j))=0.0
             u(nprim+nimp+3,grid_point(j))=rbp(grid_point(j))

         else if(nameu_gp(j)(1:2) .eq. nameu(nprim+nimp+4))then     !wr
            angrot(grid_point(j)) =   y_ode(j)*typx(j)
            u(nprim+nimp+4,grid_point(j))= y_ode(j)*typx(j)
c            print *,'angrot ,j =',angrot(grid_point(j)),j,
c     .                                      grid_point(j)
         else
            CALL STOP('subroutine eval_ode2:nameu_gp error', 6)
         endif
       enddo

       call update (u, en, te, ti, rbp, nk, nj, kj, kk, 
     .                                         iangrot, angrot)
c
c
c
c

      call rhomsh (ode_time)

c
c     take care of bc:
c
      call set_boundary_condition (ode_time)

      call set_error_field(timnew)

      call specify  ! set values of variables run in analysis mode

      call zen      ! get ene and zeff

       call diffus (xi_include)
c      print *,'done diffus'
       call source
c      print *,'done source'
       call fluxx
c      print *,'done fluxx'
       call solve                               !uses u(j,k) not te(j),etc.
c       print *,'done solve'
c      nonlin.i now has a common block with the residuals,resid,
c      stored in it.

       do j=1,no_eq
           ode_rhs(j) = resid(j)
c        if(grid_point(j) .eq. 1 .and. nameu_gp(j)(1:2) .eq.
c     .       nameu(nprim+nimp+3))
c     .            print *,'rbp resid at 0:',resid(j)
       enddo



      return
c
      end

      subroutine eval_ode3 (no_eq, ode_time, y_ode, dydt,
     .                      resid_ode,ires)
c
c -----------------------------------------------------------------------
c - This subroutine  is called when the (longitudinal) method of lines
c - is used to solve the diffusion equations using dlsoibt.
c - The residual  of
c        resid_ode = g(t,y) - A*t,y)*dy/dt
c   is returned in resid_ode
c   the structure of y_ode and resid_ode is determined
c   by subroutine get_vars. 
c   The assumed sequence is always nprim then nimp then te then ti then rbp
c   then rotation. This sequence must be strictly adhered to(variables that
c   are not run in simulation mode in this sequqnce are simply skipped over)
c   For example suppose we have
c   nkt=3 dependent variables,Te,Ti and rbp.
c   Since rbp at grid point  is laways zero it is skipped.
c     y_ode(1) = sqrt(te(1))/typx(1)
c     y_ode(2) = sqrt(ti(1))/typx(2)
c     y_ode(3) = sqrt(rbp(1))/typx(3)  NOTE special equation for rbp at rho = 0 !!!!
c     y_ode(4) = sqrt(te(2))/typx(4)
c     y_ode(5) = sqrt(ti(2))/typx(5)
c     y_ode(6) = sqrt(rbp(2))/typx(6)
c     y_ode(7) = sqrt(te(3))/typx(7)
c     y_ode(8) = sqrt(ti(3))/typx(8)
c     y_ode(9) = sqrt(rbp(3))/typx(9)
c      and so on all the way to the edge grid points
c      Note that the edge grid point for each variable can be different and
c      the structure in y_ode reflects this.
c     note that for rotation we have
c     y_ode = angrot/typx rather than Sqrt(angrot)/typx since rotation
c              can be negative

c
c  INPUT
c  no_eq         # of equations to be solved 
c  ode_time      the time at which ode_rhs is to be evaluated
c  y_ode         the dependent variable vector at time t
c  dydt      derivative of y_ode
c
c  input through common:
c  typx()        scale factors
c  grid_point()  grid point number
c
c  OUTPUT
c  resid_ode (no_eq)   residual of each equation

c -----------------------------------------------------------------HSJ-----
c
      USE param
      USE aid_newton
      USE ions
      USE io
      USE solcon
      USE soln
      USE extra
      USE numbrs
      USE mesh
      USE verbose
      USE sourc
      USE machin
      USE geom
      USE flags
      USE tordlrot
      USE soln2d
      USE bd_condtn
      USE rhog
      USE nonlin
      USE flx
      USE ode_pkg
c      USE  io,                        ONLY : ioftn, iopntr, nunits
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      dimension y_ode(*),resid_ode(*),dydt(*)
c


c
      dimension  y_odeb(kk)
c

      time         = ode_time
      write (6, '("time_ode (dlsodi,dlsodibt) =", 1pe16.8)') ode_time
      print *,'eval_ode3 te(1),ti(1) =',te(1),ti(1)

c
c     given the dependent variable values that are run in simulation
c     mode in y_ode load the usual vectors used in SOURCE, etc.:
c
      do j=1,no_eq
         do i = 1,nprim
              if(nameu_gp(j)(1:2) .eq.  nameu(i))then
                 !primary ion densitites
              !en(grid_point(j),?) implement
           CALL STOP('subroutine eval_ode3:density not implemented', 6)
              endif
         enddo
         do i =1 ,nimp
            if(nameu_gp(j)(1:2) .eq.  nameu(nprim+i))then
               !impurity densities
               !eni(grid_point(j),?) implement
           CALL STOP('subroutine eval_ode3:density not implemented', 6)
            endif
         enddo 
         if(nameu_gp(j)(1:2) .eq. nameu(nprim+nimp+1))then          !te
            te(grid_point(j)) =  (y_ode(j)*typx(j))**2
c            print *,'te ,j =',te(grid_point(j)),j,grid_point(j)
            u(nprim+nimp+1,grid_point(j))=(y_ode(j)*typx(j))**2
         else if(nameu_gp(j)(1:2) .eq. nameu(nprim+nimp+2))then     !ti
            ti(grid_point(j)) =  (y_ode(j)*typx(j))**2 
            u(nprim+nimp+2,grid_point(j))=(y_ode(j)*typx(j))**2
c            print *,'ti ,j =',ti(grid_point(j)),j,grid_point(j)
         else if(nameu_gp(j)(1:2) .eq. nameu(nprim+nimp+3))then     !rbp
c            print *,'typxrbp(j) =',typx(j)
            rbp(grid_point(j)) =  (y_ode(j)*typx(j))**2
            if(grid_point(j) .eq.1)rbp(grid_point(j)) =0.0d0
            u(nprim+nimp+3,grid_point(j))=rbp(grid_point(j))
c            print *,'rbp ,j =',rbp(grid_point(j)),j,grid_point(j)
         else if(nameu_gp(j)(1:2) .eq. nameu(nprim+nimp+4))then     !wr
            angrot(grid_point(j)) =   y_ode(j)*typx(j)
            u(nprim+nimp+4,grid_point(j))= y_ode(j)*typx(j)
c            print *,'angrot ,j =',angrot(grid_point(j)),j,
c     .                                      grid_point(j)
         else
            print *,'nameu_gp(j),eqtn no.=',nameu_gp(j),j
            CALL STOP('subroutine eval_ode3:nameu_gp error', 6)
         endif
       enddo
       call update (u, en, te, ti, rbp, nk, nj, kj, kk, 
     .                                         iangrot, angrot)
c
c
c
c

      call rhomsh (ode_time)

c
c     take care of bc:
c
      call set_boundary_condition (ode_time)

      call set_error_field(timnew)

      call specify  ! set values of variables run in analysis mode

      call zen      ! get ene and zeff

       call diffus (xi_include)
c      print *,'done diffus'

       call source
c      print *,'done source'
       call fluxx
c      print *,'done fluxx'

       call resid_ml(dydt)

c      nonlin.i now has a common block with the residuals,resid,
c      stored in it (done in ml_resid2 which is called from resid_ml)

       do j=1,no_eq
            resid_ode(j) = resid(j)
       enddo



      return
c
      end

      subroutine eval_ode (neq1, ode_time, y_odel, y_odeprime)
c
c ----------------------------------------------------------------------
c - Subroutine evaluates the time derivative of the dependent variables.
c - This subroutien is called when the (longitudinal) method of lines
c - is used to solve the diffusion equations using Dgear.
c - The rhs of
c         dy/dt = -(1/hr)*(hr(q+c))+s
c   is returned in y_odeprime
c   the structure of y_odel and y_odeprime is as follows
c   let nsolve_pde be the number of coupled pdes to be solved
c   then we take y_odel(1)...y_odel(nsolve_pde) to be the dependent variable
c   vector at grid point 1,
c   y_odel(nsolve_pde+1),.....y_odel(nsolve_pde+nsolve_pde)
c                             is the dependent variable vector
c   at grid point 2,
c   and in general we have
c   y_odel((j-1)*nsolve_pde+1),.....y_odel(j*nsolve_pde)
c                                   is the dependent variable vector
c   associated with grid point j.
c   There are nk dependent variables. However only nkt of these are
c   run in simulation mode. The ones that aren't must
c   get values at each time and space point by interpolation from the
c   data in the inone file (sub SPECIFY takes care of this).
c
c  INPUT
c  neq1        # of equations to be solved (= nkt)
c  ode_time    the time at which y_odeprime is to be evaluated
c  y_odel      the dependent variable vector at time t
c              note that y_odel in this routine IS NOT !!!!!!!
c              The same  y_ode that is passed to DGEAR
c
c  OUTPUT
c  y_odeprime(neq1)  the time derivative of y
c -----------------------------------------------------------------HSJ-----
c
      USE param
      USE fusion
      USE ions
      USE solcon
      USE soln
      USE neut
      USE mhdpar
      USE mhdgrid
      USE nub
      USE nub2
      USE rf
      USE io   
      USE extra
      USE numbrs
      USE mesh
      USE verbose
      USE sourc
      USE machin
      USE tfact
      USE geom
      USE flags
      USE tordlrot
      USE soln2d
      USE tcoef
      USE bd_condtn
      USE mixcom
      USE rhog
      USE flx
      USE ode_pkg
      USE oloss_dat
c      USE  io,                        ONLY : ioftn, iopntr, nunits

      USE gpsi
      USE island
      USE pelcom
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      dimension y_odel(*), y_odeprime(*)
c

      include 'imsl.i'


c      include 'island.i'


c      include 'pelcom.i'

      include 'sxrcom.i'

      include 'rebut.i'

c      include 'small.i'

c
      dimension ent(kj), y_odeb(kk)
c
      rather_small = 0.0001 ! for portability .. jf .. 06nov97
      nsolve_pde   = nkt
      time         = ode_time
c
      write (6, '("time_ode 1 =", 1pe16.8)') ode_time
      if (nsolve_pde .eq. 0)
     .  call STOP ('subroutine EVAL_ODE: nsolve_pde = 0', 265)
c
      call rhomsh (ode_time)
c
c     take care of bc:
c

      call set_boundary_condition (ode_time)
      call set_error_field(timnew)
      call specify  ! set values of variables run in analysis mode
      call zen      ! get ene and zeff
c
      do j=1,nj
        ent(j) = 0.0
        do i=1,nion
          ent(j) = ent(j)+en(j,i)           ! total ion density
        end do
      end do
c
      if (rho_edge .lt. 1.0) then
           nj_save = nj
****       nj      = nj_rho_edge
           call copya(ub,ub_save,kk)     ! save ub
           call copya(ub_rho_edge,ub,kk) ! set bc at rho_edge using ub
      end if
c
c
c     load boundary values for analysis and simulation variables:
       do i=1,nk
         y_odeb(i) = 0.0
         if (itran(i) .gt. 0) then
           y_odeb(i) = ub(i)
           if (i .eq. nion+1)
     .       y_odeb(i) = y_odeb(i)*g2c_ode*
     .         ene(nj_rho_edge)/ene_scale ! y_odeb(nion+1)=1.5NeTe, etc
           if (i .eq. nion+2)
     .       y_odeb(i) = y_odeb(i)*g2c_ode*
     .         ent(nj_rho_edge)/eni_scale
c           if (i .eq. nion+3)
              !rbp
c           if (i .eq. nion+4)
              !edge rotation
         end if
       end do


c
c     given the dependent variable values that are run in simulation
c     mode in y_odel load the usual vectors used in SOURCE, etc.:
c     (Vectors not loaded here were loaded in subroutine SPECIFY)
c
      ieq_ode = 0
      do j=1,nj_rho_edge
         do i=1,nk
            if (itran(i) .gt. 0) then
                ieq_ode = ieq_ode+1
                if (i .le. nion ) then
c                     not implemented
                else if (i .eq. nion+1) then
                      te(j) = ene_scale*y_odel(ieq_ode)/(g2c_ode*ene(j))
                      te(j) = MAX (te(j), rather_small)
                      if (j .eq. nj_rho_edge)
     .                te(j) = ene_scale*y_odeb(nion+1)/(g2c_ode*ene(j))
                else if (i .eq. nion+2) then
                      ti(j) = eni_scale*y_odel(ieq_ode)/(g2c_ode*ent(j))
                      ti(j) = MAX (ti(j), rather_small)
                      if (j .eq. nj_rho_edge)
     .                ti(j) = eni_scale*y_odeb(nion+2)/(g2c_ode*ent(j))
                else if (i .eq. nion+3) then
c                     not implemented
c                      if (j .eq. nj_rho_edge)
c                     rbp(j) = ??
                else if (i .eq. nion+4) then
c                     not implemented
c                      if (j .eq. nj_rho_edge)
c                     angrot(j) = 
                end if
             end if
          end do      ! i=1,..nk
       end do
c
c
c
c
       call diffus (xi_include)
       call source
       call fluxx



c-------------at mag axis, r(1) section---------------------------------
c
c     at r(1) = 0.0 we must satisfy the bc of zero flux which leads
c     to slightly different set of equations:
c
         ieq_ode = 0
         hr32    = (hcap(1)+hcap(2))*r(2)*0.25         ! r(1) = 0.0
         dra     = dr(1)*0.5
         coef    = -2.0 / (hcap(1)*dra*dra)
         do i=1,nk  !loop over all (not just simulation) variables
            if (itran(i) .gt. 0) then
                ieq_ode = ieq_ode+1
                if (i .le. nion) then
c                  not implemented
                else if (i .eq. nion+1) then
                   qed32 = 0.0
                   do k=1,nion
                      qed32 = qed32-d(i,k,1)*(en(2,k)-en(1,k))/r(2)
                   end do
                   qed32 = qed32-d(i,i,1)*(te(2)-te(1))/r(2)
                   qed32 = qed32-d(i,i+1,1)*(ti(2)-ti(1))/r(2)
                   qec32 = g1e_ode*fluxe(1)*te(1)
                   y_odeprime(ieq_ode) = coef*(hr32*(qed32+qec32))
     .                      - qdimpl(1)*(te(1)-ti(1))+s(nion+1,1)
                   y_odeprime(ieq_ode) = y_odeprime(ieq_ode)/ene_scale
                else if (i .eq. nion+2) then
                   qid32 = 0.0
                   do k=1,nion
                      qid32 = qid32-d(i,k,1)*(en(2,k)-en(1,k))/r(2)
                   end do
                   qid32=qid32-d(i,i,1)*(ti(2)-ti(1))/r(2)
                   qid32=qid32-d(i,i-1,1)*(te(2)-te(1))/r(2)
                   qic32=g1i_ode*fluxi(1)*ti(1)
                   y_odeprime(ieq_ode)= coef*(hr32*(qid32+qic32))
     .                       +qdimpl(1)*(te(1)-ti(1))+s(nion+2,1)
                   y_odeprime(ieq_ode)=y_odeprime(ieq_ode)/eni_scale
                else if (i .eq. nion+3) then
c                  not implemented
                else if (i .eq. nion+4) then
c                  not implemented
                end if
             end if
          end do  ! i=1,,..nk
c --------------end section at r(1) -----------------------------


c-------interior section,r(2) to r(nj_rho_edge-1) ------------------------
c
c     for the remaining (interior) grid points.
c     The equations for grid point nj-1 will involve
c     grid point nj. The appropriate values of y_odel at nj_rho_edge
c     were loaded above.
c
      ieq_ode=nsolve_pde
      do j=2,nj_rho_edge -1
         hrap=(hcap(j+1)*r(j+1)+hcap(j)*r(j))*0.5
         hram=(hcap(j)*r(j)+hcap(j-1)*r(j-1))*0.5
         hrj=-hcap(j)*r(j)*(ra(j)-ra(j-1))      ! r(j+1/2) -r(j-1/2)
         do i=1,nk
           if (itran(i) .gt. 0) then
               ieq_ode=ieq_ode+1
c
               if (i .le. nion) then
c
c              ion densities, not implemented at this time
c
               else if (i .eq. nion+1) then ! Te
c
c                  electron cond + convct heat flux grid point j+/-.5:
c
                   qerp= 0.0
                   qerm=0.0
                   do k=1,nion
                      qerp=qerp-d(i,k,j)*(en(j+1,k)-en(j,k))/dr(j)
                      qerm=qerm-d(i,k,j-1)*(en(j,k)-en(j-1,k))/dr(j-1)
                   end do
                   qerp=qerp-d(i,i,j)*(te(j+1)-te(j))/dr(j)
                   qerp=qerp+0.5*(te(j+1)+te(j))*fluxe(j)
                   qerm=qerm-d(i,i,j-1)*(te(j)-te(j-1))/dr(j-1)
                   qerm=qerm+0.5*(te(j)+te(j-1))*fluxe(j-1)
                   y_odeprime(ieq_ode)=(hrap*qerp-hram*qerm)/hrj
                   y_odeprime(ieq_ode)=y_odeprime(ieq_ode)+s(nion+1,j)
     .                         - qdimpl(j)*(te(j)-ti(j))
                   y_odeprime(ieq_ode)=y_odeprime(ieq_ode)/ene_scale
               else if (i .eq. nion+2) then ! Ti
c
c                  ion cond + convct heat flux grid point j+/-.5:
c
                   qirp=0.0
                   qirm=0.0
                   do k=1,nion
                      qirp=qirp-d(i,k,j)*(en(j+1,k)-en(j,k))/dr(j)
                      qirm=qirm-d(i,k,j-1)*(en(j,k)-en(j-1,k))/dr(j-1)
                   end do
                   qirp=qirp-d(i,i,j)*(ti(j+1)-ti(j))/dr(j)
                   qirm=qirm-d(i,i,j-1)*(ti(j)-ti(j-1))/dr(j-1)
                   qirp=qirp+0.5*(ti(j+1)+ti(j))*fluxi(j)
                   qirm=qirm+0.5*(ti(j)+ti(j-1))*fluxi(j-1)
                   y_odeprime(ieq_ode)=(hrap*qirp-hram*qirm)/hrj
                   y_odeprime(ieq_ode)=y_odeprime(ieq_ode)+s(nion+2,j)
     .                         + qdimpl(j)*(te(j)-ti(j))
                   y_odeprime(ieq_ode)=y_odeprime(ieq_ode)/eni_scale
               else if (i .eq. nion+3) then ! rbp
c                not implemented
               else if (i .eq. nion+4) then !toroidal rotation
c                not implemented
               else
                 call STOP ('subroutine EVAL_ODE: fell thru IF', 266)
               end if
            end if    ! itran .gt. 0
         end do       ! do loop on i=1,nk
       end do         ! do loop on j=2,nj_rho_edge-1
c
****  if (rho_edge .lt. 1.0) then
****    nj=nj_save
****    call copya(ub_save,ub,kk)           ! restore ub
****  end if
c
      return
c
      end

      subroutine eval_ode_jacob (neq, ode_time, y_odel, pd)
c
c ----------------------------------------------------------------------
c
c     evaluates the jacobian of y_odel wrt y
c     This is a dummy routine. The jacobian is determined
c     numerically for all method of lines solvers
c     (Dgear,Radau5,dlsodi,dlsiobt). But thisd routine
c     must still br present.
c

      USE ode_pkg
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c
      return
c
      end




      subroutine get_em (j,em)
c
      USE param
      USE solcon
      USE soln
      USE     ions
      USE numbrs
      USE mesh
      USE adaptive
      USE sourc
      USE geom
      USE flags
      USE tordlrot
      USE tcoef
      USE bd_condtn
      USE flx
      USE tmpcom
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c     this subroutine generates the matrice em  for mesh point j.
c     the matrices and vectors (all of order nk) are calculated only
c     as needed for a particular mesh point in order to save storage.
c     Note that en,te,ti,rbp,angrot (the quantities at the
c     central,n+theta, time point) are used.
c ----------------------------------------------------------------------
c
c      include 'param.i'
c      include 'bcon.i'
c      include 'flags.i'
c      include 'flx.i'
c      include 'geom.i'
c      include 'invers.i'
c      include 'ions.i'
c      include 'mesh.i'
c      include 'numbrs.i'
c      include 'solcon.i'
c      include 'soln.i'
c      include 'sourc.i'
c      include 'tcoef.i'
c      include 'tmpcom.i'
c      include 'tordlrot.i'
c      include 'adaptive.i'
      
c
      dimension em(kk,kk)

c     note  kk = kion + 4 , used to be kion+3 before rotation was added
c

      real*8          kevperg
      data            kevperg /6.242e8/, xmassp /1.673e-24/

c
c     calculate time weighting factorc
c     calculate inverse time step
c
      dtinv = 1.0 / dtt
      if(diffeq_methd .eq. 1)dtinv =0.0d0 !eliminates contribution
                                          !from em to B and g dtinv is local
                                          !to this routine
c
      zero  = 0.0
      theta = 1.
      onemt = (1.0-theta)/theta
      thetas = theta
      if(ABS(steady_state-1.0) .gt. 1.e-5)then !change HSJ
         onemt = 0.0
         theta =1.0
         dtinv = 1.0   !diffeq_methd = 1 not valid for steady state solutions
      endif

c
c     common multipliers for toroidal rotation
c
      fmult   = angrcple*angrot(j)*kevperg*0.5*iangrot
      r2omega = angrot(j)*r2capi(j)
c
c     calculate hcap ratios
c
      if (j .eq. 1 )  go to 4
      hjm1 = 0.5*(hcap(j-1)+hcap(j))/hcap(j)
      if (j .eq. nj)  go to 6
    4 hjp1 = 0.5*(hcap(j+1)+hcap(j))/hcap(j)
c
c     calculate matrix em (i.e., matrix m in GAA report)
c
    6 do 10 k=1,nk
      do 10 l=1,nk
   10 em(k,l) = 0.0
      if(ABS(steady_state -1.0d0) .lt. 1.e-5)then
        sume = ene(j)
        sumi = 0.0
        do 20 k=1,nk-3-iangrot   !loop over primary and impurity densities
        em(k,k) = 1.0
        em(nk-2-iangrot,k) = 1.5*z(j,k)*te(j)
        em(nk-1-iangrot,k) = 1.5*ti(j)+fmult*r2omega*atw(k)*xmassp
        sume = sume+en(j,k)*te(j)*dzdte(j,k)
   20   sumi = sumi + en(j,k)
        em(nk-2-iangrot,nk-2-iangrot) = 1.5*sume
c        print *,'emj1,j =', em(nk-2-iangrot,nk-2-iangrot),j
        em(nk-1-iangrot,nk-1-iangrot) = 1.5*sumi
        if (iangrot .eq. 1) then
          smassden = 0.0
          do 200 k=1,nk-4
          xmassi   = atw(k)*xmassp
          smassden = smassden+xmassi*en(j,k)
  200     em(nk,k) = r2omega*xmassi
          em(nk-2,nk) = smassden*r2omega*kevperg
          em(nk,nk) = smassden*r2capi(j)
        end if
      endif



      if (j .eq. 1)  go to 30
      em(nk-iangrot,nk-iangrot) = 1.0 /
     .                (fcap(j) * gcap(j) * (hcap(j) * r(j))**2 * twkfar)
c       print *,'em here,nk =',em(nk-iangrot,nk-iangrot),nk
        if(ABS(steady_state -1.0d0) .gt. 1.e-5)
     .  em(nk-iangrot,nk-iangrot) =em(nk-iangrot,nk-iangrot)*
     .                                      steady_state
        go to 40
 30      continue
         em(nk-iangrot,nk-iangrot) = 1.0  !for rbp at rho =0.0
 40      theta =thetas
c      print *,'return from get_em'
      return
      end



        subroutine Mat_M(nqr,am,lmas,rtypx,ipar)
c-----------------------------------------------------------
c       subroutine returns the matrix M in 
c          M*dy/dt =F(t,y)
c       Matrix M is just matrix em in subroutine abcg
c       over the entire grid
c       rtypx is not used. instead typx from aid_newton is used
c------------------------------------------------------------
        USE param
        USE aid_newton
        USE soln
        USE numbrs
        USE mesh
      USE flags
      USE tordlrot
      USE nonlin
        implicit  integer (i-n), real*8 (a-h, o-z)
c        include 'param.i'
c        include 'flags.i'
c        include 'mesh.i'
c        include 'nonlin.i'
c        include 'numbrs.i'  !nj
c        include 'soln.i'
c        include 'tordlrot.i' !iangrot
        dimension a(kk,kk), b(kk,kk), c(kk,kk), g(kk)
        dimension em(kk,kk)
        integer *4 ipar(1),nqr,lmas,neqm,nvar
        real*8 am(lmas,nqr),rtypx(*)



c       M is block tridiagonal.The blocks, em, are  obtained from abcg.
c       the block size is (nkt,nkt) after analysis mode variables are
c       eliminated.



        nrbp = nprim+nimp+3 ! local index for rbp. Recall that em is obtained
c                           ! directly from abcg without reduction to nkt form.

        neqm = 0            !total # equation to be solved(= nqr at end).

        do j=1,nj-1
c           get matrix em, note that em is NOT  reduced to(nkt,nkt)  here
c           (it is of size nk,nk)):
            call abcg (j, a, b, c, g, em)
            nblk = -1       !each grid point j produces  a block matrix em
                            !of size nk by nk . nblk is equation  # relative
                            !to  block em starting with 0

            do kc = 1,nk    !down coloums of block matrix em  at grid point j
               if(itran(kc) .eq. 1)then
c                  if(j .eq. 1 .and. kc .eq.  nrbp)go to 10
c                 skip boundary values:
                  if(j .ge. te_index  .and. kc .eq.
     .                                         nk-iangrot-2)go to 10
                  if(j .ge. ti_index .and. kc .eq. 
     .                                         nk-iangrot-1)go to 10
                  if(j .ge. rot_index .and. kc .eq. 
     .                              nk .and. iangrot .eq. 1)go to 10
                     neqm = neqm +1
                     nvar = -1   !counts variables in siumulation mode for
                                 !grid point j, starting at 0 for firstr one.

                     nblk = nblk + 1  ! = 0,1,2..nkt-1 (or nkt-2 for j =1)

                     do kr =1,nk !across rows of block matrix em at gridpoint j

                        if(itran(kr) .eq.1)then
c                           if(j .eq. 1 .and . kr .eq. nrbp)go to 20
c                          skip boundary values:
                            if(j .ge. te_index  .and. kr .eq.
     .                                         nk-iangrot-2)go to 20
                            if(j .ge. ti_index .and. kr .eq. 
     .                                         nk-iangrot-1)go to 20
                            if(j .ge. rot_index .and. kr .eq. 
     .                              nk .and. iangrot .eq. 1)go to 20
                           nvar = nvar +1  !reduced variable no
                            !full matrix form:
c                           am(neqm,neqm-nblk+nvar ) = em(kc,kr) 
c     .                                  *2.*SQRT(u(kr,j))*typx(neqm)
c                          banded form:
                           if(kr .eq. nk .and. iangrot .eq. 1)then
                           am(nblk-nvar+nkt,neqm-nblk+nvar)= 
     .                               em(kc,kr)*u(kr,
     .                                 grid_point(neqm))*typx(neqm)/
     .                                 eqtn_scale(neqm)
                           else
                           am(nblk-nvar+nkt,neqm-nblk+nvar)=
     .                           em(kc,kr)*2.*SQRT(u(kr,
     .                                    grid_point(neqm)))*typx(neqm)/
     .                                          eqtn_scale(neqm)
                           if(j .eq. 1 .and. kr .eq. nk-iangrot .and. kc
     .                          .eq. nk-iangrot)
     .                           am(nblk-nvar+nkt,neqm-nblk+nvar)=1.
                                !at rho=0 for rbp we change the eqauation
                                !to avoid introduction of what would become
                                !a DAE system (which requires that initial 
                                !derivatives be input)
                                !so at rho = 0 we have
                                ! 1*d/dt(rbp) = 0.0
                                !instead of
                                ! 0*d/dt(rbp) = 0.0
                           endif
 20                        continue
                        endif

                     enddo    ! kr



 10                continue
                endif

            enddo       !kc 

        enddo           !over grid points,j

        if(neqm .ne. nqr)then
           print *,'neqm,nqr =',neqm,nqr
          call STOP("Subroutine Mat_M neqm  wrong",0)
        endif

        return
        end

      subroutine ml_resid2(jrgrid,a,b,c,g,em,dydt)
c---------------------------------------------------HSJ---------------
c      This subroutine calculates the residual for the finite differenced
c      form of the transport equations. This subroutine assumes that the
c      matrices a,b,c,g,em have been reduced to nkt by nkt !!!!
c      Each radial grid point, designated by jrgrid,introduces a set of
c      nkt unknowns where nkt is the number of itran(k) values that
c      are non zero (eg the number of dependent variables run in
c      simulation mode).

C      INPUT
C        THE MATRICES A,B,C,G,EM  FROM resid_ml
C        jrgrid      grid point at which the nkt equations
c                    are to be evaluated




c      OUTPUT
c      nequations    current equation number (the total number of equations
c                    and the total number of unknowns must be equal of course).
c      resid         residual of equation number neq,in nonlin.i
c---------------------------------------------------------------------
c
c
      USE param
      USE aid_newton
      USE soln
      USE numbrs
      USE mesh
      USE flags
      USE tordlrot
      USE nonlin
      implicit  integer (i-n), real*8 (a-h, o-z)
      integer l,k1,k2,jrgrid 
c      include 'param.i'    !kk,
c      include 'flags.i'    !itran(kk)
c      include 'numbrs.i'   !nk,nkt,nj
c      include 'mesh.i'     !te,ti,rot_index
c      include 'nonlin.i'   !resid,nequations
c      include 'soln.i'     !u(kk,kj)
c      include 'tordlrot.i' !iangrot
      real *8 sum, a(kk,kk),b(kk,kk),c(kk,kk),g(kk),dydt(*),
     .             em(kk,kk)
    
      if(jrgrid .eq. nj)return
      
      k1=0
      do k = 1,nk
         if( itran(k) .gt. 0)then
            k1 =k1 +1
            scale_eqn = 1.e0
            sum = 0.0d0          ! nkt is the number of simulation variables
c            if(jrgrid .eq. 1 .and. k .eq. nk-iangrot)
c     .                           go to 30
            if(jrgrid .ge. te_index .and. k .eq. 
     .                         nk-2-iangrot)go to 30
            if(jrgrid .ge. ti_index  .and. k .eq. 
     .                         nk-1-iangrot)go to 30
            if(jrgrid .ge. rot_index .and. k .eq. nk .and. 
     .                                   iangrot .eq. 1)go to 30
            nequations = nequations+1 
            k2 =0
            do l =1,nk
               if(itran(l) .gt. 0)then
                  k2 = k2+1
                  if(jrgrid .eq. 1)go to 10
                  sum =sum +a(k1,k2)*u(l,jrgrid-1)
 10               sum =sum + b(k1,k2)*u(l,jrgrid)
c-------------------------------------

                  if(l .eq. nk .and. iangrot .eq. 1)then
                      dumy = em(k1,k2)*u(l,jrgrid)*typx(nequations)
                  else
                      dumy  = em(k1,k2)*2.*SQRT(u(l,jrgrid))*
     .                                              typx(nequations)
                  endif
c------------------------------------
                  sumsave = sum 
                  sum = sum + dumy*dydt(nequations)
                  if(jrgrid .eq. nj) go to 20
                  sum = sum + c(k1,k2)*u(l,jrgrid+1)
                  sumsave =sumsave+ c(k1,k2)*u(l,jrgrid+1)
 20            endif
            enddo
            resid(nequations)  = ( g(k1) - sum)/eqtn_scale(nequations)
c            if(k1 .eq. 1)print *,'g(k1),j=',g(k1),jrgrid
            sumsave = g(k1) - sumsave
c            if(k1 .eq. 1)print *,'sumsave,j=',sumsave


 30      endif
       enddo
c       print *,'jrgrid,nequations',jrgrid,nequations

       return

       end




      subroutine ml_resid(jrgrid,a,b,c,g)
c---------------------------------------------------HSJ---------------
c      This subroutine calculates the residual for the finite differenced
c      form of the transport equations. This residual can be used in
c      various non linear solver schemes.
c      Each radial grid point, designated by jrgrid,introduces a set of
c      nkt unknowns where nkt is the number of itran(k) values that
c      are non zero (eg the number of dependent variables run in
c      simulation mode).

C      INPUT
C        THE MATRICES A,B,C,G FROM SOLVE
C        jrgrid      grid point at which the nkt equations
c                    are to be evaluated




c      OUTPUT
c      nequations    current equation number (the total number of equations
c                    and the total number of unknowns must be equal of course).
c      resid         residual of equation number neq,in nonlin.i
c---------------------------------------------------------------------
c
c
      USE param
      use aid_newton
      USE soln
      USE numbrs
      USE mesh
      USE flags
      USE tordlrot
      USE nonlin
      implicit  integer (i-n), real*8 (a-h, o-z)
      integer l,k1,k2,jrgrid 
c      include 'param.i'    !kk,
c      include 'flags.i'    !itran(kk)
c      include 'numbrs.i'   !nk,nkt,nj
c      include 'mesh.i'     !te,ti,rot_index
c      include 'nonlin.i'   !resid,nequations
c      include 'soln.i'     !u(kk,kj)
c      include 'tordlrot.i' !iangrot
      real *8 sum, a(kk,kk),b(kk,kk),c(kk,kk),g(kk)

      if(jrgrid .eq. nj)return
      
      k1=0
      do k = 1,nk
         if( itran(k) .gt. 0)then
            k1 =k1 +1
            scale_eqn = 1.e0
            sum = 0.0d0          ! nkt is the number of simulation variables
c            if(jrgrid .eq. 1 .and. k .eq. nk-iangrot)
c     .                           go to 30
            if(jrgrid .ge. te_index .and. k .eq. 
     .                         nk-2-iangrot)go to 30
            if(jrgrid .ge. ti_index  .and. k .eq. 
     .                         nk-1-iangrot)go to 30
            if(jrgrid .ge. rot_index .and. k .eq. nk .and. 
     .                                   iangrot .eq. 1)go to 30
            nequations = nequations+1 
            k2 =0
            do l =1,nk
               if(itran(l) .gt. 0)then
                  k2 = k2+1
                  if(jrgrid .eq. 1)go to 10
                  sum =sum +a(k1,k2)*u(l,jrgrid-1)
 10               sum =sum + b(k1,k2)*u(l,jrgrid)
c                  if(l .eq. nk -1)print *,'b,k1,k2=', b(k1,k2),k1,k2
                  if(jrgrid .eq. nj) go to 20
                  sum = sum + c(k1,k2)*u(l,jrgrid+1)
 20            endif
            enddo
            resid(nequations)  =( -sum + g(k1))/eqtn_scale(nequations)
c           note that resid(i) does not have to be small here
c           because it does not include M*dy/dt

 30      endif
       enddo
c       print *,'jrgrid,nequations',jrgrid,nequations
       return

       end






      subroutine output_sub(nr,told,t,y,cont,lrc,n,
     .                      rpar,ipar,irtrn)
c-------------------------------------------------------------------
c     output subroutine called by radua5
c     not currently used
c-------------------------------------------------------------------
      implicit none
      integer*4 nr,lrc,n,irtrn,ipar
      real *8 told,rpar,cont(lrc),t,y(n)

      return 
      end






      subroutine reduce_em (em, itran, nk, kk)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
      integer diffeq_methd
c
c  this subroutine reduces the order of the matrices a, b, and c
c     and the vector g from nk to nkt.
c
      dimension em(kk,*),itran(*)

c  reduce em  by eliminating all rows and columns
c     corresponding to specified variables
c
      kt = 0
      do 50 k=1,nk
      if (itran(k) .le. 0)  go to 50
      kt = kt + 1
      do 40 l=1,nk
 40      em(kt,l)= em(k,l)
      do 45 l=1,nk
 45   em(l,kt)= em(l,k)
   50 continue
      return
c
      end







      subroutine resid_ml(dydt)
c------------------------------------------------------------------
c
c
c
c  this subroutine solves for the vector uu that satisfies the
c  linear matrix equation
c       aa*uu = gg
c  the matrix aa and the vectors uu and gg are of order nj*nkt.
c  moreover, the matrix aa is block tridiagonal with nj blocks
c  of order nkt.  the jth block-matrix equation can be written as
c       a(j)*u(j-1) + b(j)*u(j) + c(j)*u(j+1) = g(j) .
c
c  NOTE:
c  kk  = leading dimension of u
c  nk  = number of dependent variables in the problem (as opposed to
c        the total number of dependent variables the code allows, which
c        is set to kk)
c        The value of nk depends on how many primary ions
c        and impurities there are in the input (counting both simulation and
c        analysis mode variables)
c        nk = #prim ions  +#impurity ions   + 1 for te + 1 for ti
c                           +1 for rbp  + 1 for rotation
c        Inclusion of rotation is optional,equations for te,ti,rbp
c        and at least one primary ion are always present so min value of
c        nk is 4. Max value of nk is kk ( or kk-1 if rotation is excluded).
c
c  nkt = number of dependent variables for which a transport solution
c        is to be found (i.e., the number of variables for which
c        itran(i),i = 1,2..nk, is equal to 1).
c  the matrices a, b, c,and em and the vector g are calculated by the
c        subroutine ABCG and are of order nk.  if nkt < nk (i.e., if
c        the values of some dependent variables are specified), then
c        the arrays are reduced to order nkt before the matrix equation
c        is solved.  subsequently the solution vector u is expanded
c        back to order nk.
c
c  because aa is block tridiagonal, the solution uu can be obtained
c     efficiently by a direct method.  specifically, gaussian
c     elimination (i.e., forward elimination followed by backward
c     substitution) is used.  the kth element of the solution
c     vector for block j is stored in u(k,j).
c---------------------------------------------------------------HSJ--
c 
      USE param
      USE     io
      USE   m_lines
      USE soln
      USE numbrs
      USE mesh
      
      USE flags
      USE tordlrot
      USE bd_condtn,only : bctime,ub,fluxb,
     .    ub_save,ub_rho_edge,bctime_zone
      USE nonlin
      implicit  integer (i-n), real*8 (a-h, o-z)
c      include 'param.i'
c      include 'bcon.i'
      include 'imsl.i'
c      include 'flags.i'
c      include 'mesh.i'  !te,ti,rot_index
c      include 'numbrs.i'
c      include 'soln.i'
c      include 'tordlrot.i'
c      include 'io.i'
c      include 'nonlin.i'

c
      parameter (kk1 = kk+1)
c

      dimension  a(kk,kk),b(kk,kk),c(kk,kk),g(kk),em(kk,kk)
      dimension dydt(*)
      data       iacc/6/


      nequations  = 0             !counts number of equations  

      do j=1,nj
c
c  calculate the matrices a, b, c,and em, and the vector g
c  the matrix em is required explicitely for the method of
c  lines solution. In that case (specified by diffeq_methd =1),
c  em is not included in matrix b and vector g.
      call abcg (j, a, b, c, g, em)



c     redefine a,b,c,g for those profiles which have the boundary set
c     at a value less than nj:
      if(j .ge. te_index  .and. itran(nk-2-iangrot) .eq. 1)then
         kij = nion+1    !index into matrices for te 
         do ijk =1,nk   ! over row kij of a,b,c,g
            a(kij,ijk)= 0.0
            b(kij,ijk)= 0.0
            c(kij,ijk)= 0.0
            em(kij,kij) =0.0
         enddo 
         b(kij,kij)= 1.
         g(kij) = bctime_zone(j,kion+1)

      endif
      if(j .ge. ti_index .and. itran(nk-1-iangrot) .eq. 1)then
         kij = nion+2    !index into matrices for ti 
         do ijk =1,nk   ! over row kij of a,b,c,g
            a(kij,ijk)= 0.0
            b(kij,ijk)= 0.0
            c(kij,ijk)= 0.0
            em(kij,kij) =0.0
         enddo
         b(kij,kij)= 1.
         g(kij) = bctime_zone(j,kion+2)

      endif
      if(j .ge. rot_index  .and.iangrot .eq. 1
     .                .and.             itran(nk) .eq. 1)then
         kij = nion+4    !index into matrices for toroidal rotation
         do ijk =1,nk   ! over row kij of a,b,c,g
            a(kij,ijk)= 0.0
            b(kij,ijk)= 0.0
            c(kij,ijk)= 0.0
         enddo
         b(kij,kij)= 1.
         g(kij) = bctime_zone(j,kion+4)
      endif


c
c  reduce the arrays to order nkt if nkt < nk
c
      if (nkt .lt. nk)
     .call reduce (j, a, b, c, g,em, u, itran, nk, nj, kk,diffeq_methd)

      call ml_resid2(j,a,b,c,g,em,dydt)

 

      enddo


      return
      end
c method_lines.f



      subroutine solve_ml2
c------------------------------------------------------------------------
c
c     use the method of (longitudinal) lines to solve the equations
c     (a completed and updated version of solve_ml)
c
c
c------------------------------------------------------------------HSJ----

      USE param
      USE aid_newton               !freeze_*,jac_known
      USE m_lines
      USE ions
      USE solcon     !theta
      USE soln
      USE io
      USE extra
      USE numbrs
      USE mesh
      USE verbose
      USE sourc
      USE machin
      USE geom
      USE flags
      USE tordlrot
      USE soln2d
      USE bd_condtn
      USE rhog
      USE nonlin
      USE ode_pkg
      USE replace_imsl,                ONLY : my_dgear
c      USE  io,                        ONLY : ioftn, iopntr, nunits

      implicit  integer (i-n), real*8 (a-h, o-z)



c
      parameter (maxpde = kprim+3)
      parameter (kw_ode = 17*maxpde*kjm1+maxpde*kjm1*maxpde*kj)
      real *8  x(kk)
      real *8 ,save,dimension(:),allocatable ::wk_ode
      real *8 ,save, dimension(:),allocatable ::y_ode,y_ode_old
      real *8 ,save, dimension(:),allocatable ::ydot_ode
      integer *4 ,save,dimension(:),allocatable ::itran_save
      integer *4,save, dimension(:), allocatable :: iwk_ode
      logical maxscale,R5_matsing
      data maxscale,R5_matsing/.false. ,.false. /

c
      data index_ode /1/
c     initialized in block data:
c      data g1e_ode, g1i_ode, g2c_ode! flux => 5/2, energy density => 3/2
c     .    /    2.5,     2.5,     1.5/
c
      external eval_ode,  eval_ode_jacob,Mat_M,
     .         output_sub,eval_ode2,eval_ode3,addem,addem2
c

c


c
c --- set initial condition
c



      if( .not. allocated(itran_save))
     .                   allocate (itran_save(nk),STAT = istat)
      init_nl = 1 !make sure problem  is un-frozen start of cycle
      call set_freeze(freeze_type,itte,itti,itangrot,itran,
     .                freeze_xte,freeze_xti,freeze_xwr,
     .                freeze_xrbp,freeze_xni,itran_save,
     .                itti_save,itte_save,itangrot_save,
     .                init_nl,iodeflg1,nk,iangrot,jac_known)

c       we cannot change the boundary location 
c       while the ode solver is active. hence we define:
        if(index_ode .ne. 1 .and.    !index_ode = 1 above in data statement
     .        (te_index .ne. te_index_save .or. ti_index .ne.
     .       ti_index_save .or. rot_index .ne. rot_index_save))then
          !get mapping and scaling of dependent variables,
          !we will solve neq_ode equations:
          call get_vars(maxscale,xvset,typx,
     .                    nameu_gp,grid_point,x,y_ode,y_ode_old,neq_ode)
          te_index_save = te_index
          ti_index_save = ti_index
          rot_index_save = rot_index
          index_ode = 1    !remap the problem, since boundary location cahnged
        endif
      freeze_te_index =.True.
      freeze_ti_index =.True.
      freeze_rot_index=.True.
 1    timestep_ode=dt 
      time_ode = time
      timeend_ode=time_ode+dt
      if (index_ode .eq. 1) then !very first call to solve_ml or boundary
                                     !location changed:
          index_ode = 0
          if( .not. allocated(itranrd))then
              allocate (itranrd(nkt),STAT = istat)
              if(istat .ne. 0)
     .          call allocate_error("itranrd, solve_ml2",0,istat)
              itranrd(:) = 0
          endif
          if( .not. allocated(typx))then
              allocate (typx(nkt*nj),STAT = istat)
              if(istat .ne. 0)
     .          call allocate_error("typx, solve_ml2",0,istat)
              typx(:) = 0.0d0
          endif
          if( .not. allocated(eqtn_scale))then
              allocate (eqtn_scale(nkt*nj),STAT = istat)
              if(istat .ne. 0)
     .         call allocate_error("eqtn_scale, solve_ml2",0,istat)
              eqtn_scale(:) = 0.0d0
          endif
          if( .not. allocated(nameu_gp))then
              allocate (nameu_gp(nkt*nj),STAT = istat)
              if(istat .ne. 0)
     .        call allocate_error("nameu_gp, solve_ml2",0,istat)
              nameu_gp(:) = ' '
          endif
          if( .not. allocated(grid_point))then
              allocate (grid_point(nkt*nj),STAT = istat)
              if(istat .ne. 0)
     .        call allocate_error("grid_point, solve_ml2",0,istat)
              grid_point(:) = 0
          endif
          if( .not. allocated(y_ode))then
              allocate (y_ode(nkt*nj),STAT = istat)
              if(istat .ne. 0)
     .          call allocate_error("y_ode, solve_ml2",0,istat)
          endif
          if( .not. allocated(y_ode_old))then
              allocate (y_ode_old(nkt*nj),STAT = istat)
              if(istat .ne. 0)
     .          call allocate_error("y_ode_old, solve_ml2",0,istat)
          endif


           nrmax = nj-1 ! values on boundary are known, 
                        ! hence do not form part of the equation set
           theta = 1.   !fully implicit solution always
          !get mapping and scaling of dependent variables,
          !we will solve neq_ode equations:
          call get_vars(maxscale,xvset,typx,
     .                    nameu_gp,grid_point,x,y_ode,y_ode_old,neq_ode)
          te_index_save = te_index     !save current boundary location
          ti_index_save = ti_index
          rot_index_save = rot_index



             do j=1,neq_ode
                  y_ode(j) = y_ode(j)/typx(j)
                 if(nameu_gp(j)(1:2)  .eq. nameu(nprim+nimp+1))
     .               eqtn_scale(j) = s(nprim+nimp+1,nj/2)     !kev/cm**3sec,te
                 if(nameu_gp(j)(1:2)  .eq. nameu(nprim+nimp+2))
     .               eqtn_scale(j) = s(nprim+nimp+2,nj/2)    !kev/cm**3sec,ti
                 if(nameu_gp(j)(1:2)  .eq. nameu(nprim+nimp+3))
     .               eqtn_scale(j) = s(nk-iangrot,nj/2)           !gauss cm/sec
                 if(nameu_gp(j)(1:2)  .eq. nameu(nprim+nimp+4)
     .             .and.    iangrot .eq. 1)
     .               eqtn_scale(j) = s(nk,nj/2)        !

             enddo
             l=0
             do j =1,nk
                if(itran(j) .eq. 1)then
                  l=l+1
                  itranrd(l) = j ! for reduced variables set ,use as index
                 endif           !into u; u(itranrd(l),j) is fit variable l
              enddo              !at gridpoint j


          tol_ode = 1.e-4    !cant make this much bigger
          atol = 1.e-3       !cant make this much bigger
          itol = 0              !tol_ode, atol are scalars



          if(non_lin_method .eq. 1)then    !Radau5
            itol = 0              !tol_ode, atol are scalars
            ijac =0               !compute jacobian in radau5
            mljac = 2*nkt-1       !lower jacobian bandwidth
            mujac = mljac         !upper     
            imas = 1              !matrix M is supplied in subroutine Mat_M
            mlmas = nkt-1         !bandwidth of matrix M (M is a banded matrix)
            mumas = mlmas
            iout =0               ! =0 dont call output_sub,= 1 call it
            lwork = mljac+mujac+1+neq_ode+3*(2*mljac+mujac+1)+12
            lwork = neq_ode*lwork + 20  !dim of wk_ode
            liwork = lwork
          else if(non_lin_method .eq. 2)then  !Dlsoibt
            !check if the input is appropriate for this case:
            ierr = 0
            maxpt =MIN(te_index,ti_index,rot_index,nj)
            if (maxpt .ne. nj .and. itran(nprim+nimp+3)
     .                .eq. 1 )then
               print *,' ERROR, dlsoibt requires that'
               print *,' 1)the boundary conditions are given at '
               print *,'   rho =1.0 for all profiles'
               print *,'            OR'
               print *,' 2) The current profile is not run in'
               print *,'    simulation mode. In this case the'
               print *,'    profiles that are run in simulation'
               print *,'    mode must have a a single common boundary'
               print *,'    which may be at rho < 1.0.'
               ierr = 1
             else if (maxpt .ne. nj .and. itran(nprim+nimp+3)
     .                .eq. 0 )then
               !all boundaries must be the same so block structure
               !of dlsoibt is maintained:
               do j=1,nk
                   !add density cases here
                   if(j .eq. nprim+nimp+1 .and. itran(j) .eq. 1
     .             .and. te_index .ne. maxpt)ierr=1
                   if(j .eq. nprim+nimp+2 .and. itran(j) .eq. 1
     .             .and. ti_index .ne. maxpt)ierr=1
                   if(j .eq. nprim+nimp+4 .and. itran(j) .eq. 1
     .             .and. rot_index .ne. maxpt)ierr=1
               enddo
               if(ierr .eq. 1)then
                 print *,' ERROR, dlsoibt requires that'
                 print *,' all profiles run in simulation mode'
                 print *,' must have the same value of rho'
                 print *,' for the boundary'
               endif
             endif
             if(ierr .eq. 1)  call STOP('solve_ml2,input invalid',0)
            itol = 1             !tol_ode, atol are scalars
            lwork = 22+9*neq_ode + 3*nkt*nkt*(nj-1)   !dim of wk_ode
            liwork = 20 + neq_ode
c            print *,'lwork,liwork =',lwork,liwork
            itask =1     !normal computation of output values
            istate = 0   ! 0 means initial dy/dt not known, 1 means it is supplied
                      !use 0 on frst call and 1 on subsequent calls??
            iopt3 = 0     !no optional inputs
            mf3 = 22     !determine jacobian using finite differences.

          else if(non_lin_method .eq. 4)then     !for dlsodi
            mf3 = 25             !determine banded jacobian using finite differences.
            itol = 1             !tol_ode, atol are scalars
            ml =  2*nkt-1        !lower bandwidth of matrices involved (not just matrix em)
            mu = ml              !upper
            lwork = 22 +10*neq_ode + neq_ode*(2*ml+mu)   !dim of wk_ode for banded jacob
            liwork = 20 + neq_ode
c            print *,'lwork,liwork =',lwork,liwork
            itask =1     !normal computation of output values
            istate = 0   ! 0 means initial dy/dt not known, 1 means it is supplied
            iopt3 = 0     !no optional inputs


          else
               call STOP("SUBROUTINE solve-ml2,no method selected",0)
          endif
      end if          ! index_ode .eq. 1




c
c
c
c

c

c


      if(non_lin_method .eq. 1)then  !Radau5
         if( .not. allocated(wk_ode))then
           allocate (wk_ode(lwork),STAT = istat)
           if(istat .ne. 0)
     .        call allocate_error("wk_ode, method_lines2",0,istat)
           wk_ode(:) =0.0d0
           wk_ode(3) = 0.1  !frequency of calling jacobian( .1 less, .001 more)
         endif
         if( .not. allocated(iwk_ode))then
           allocate(iwk_ode(liwork),STAT = istat)
           if(istat .ne. 0)
     .        call allocate_error("iwk_ode, method_lines2",0,istat)
           iwk_ode(:) = 0.0d0
         endif
c         print *,'lwork,liwork =',lwork,liwork
         call Radau5(neq_ode,eval_ode2,time_ode,y_ode,timeend_ode,
     .               timestep_ode,tol_ode,atol,itol,eval_ode_jacob,
     .               ijac,mljac,mujac,Mat_M,imas,mlmas,mumas,
     .               output_sub,iout,wk_ode,lwork,iwk_ode,liwork,
     .               typx,grid_point,idid)
          if(idid .lt. 0)then
             write(ncrt,'("radau5 reports error,idid =",i5)')idid
             SELECT CASe(idid)
                CASE(-1)
                   write(ncrt,'("ERROR in input")')
                CASE(-2)
                   write(ncrt,'("Nmax too small")')             
                CASE(-3)
                   write(ncrt,'("Step size too small")') 
                     init_nl = 0
                     call set_freeze(freeze_type,itte,itti,itangrot,
     .                       itran,freeze_xte,freeze_xti,freeze_xwr,
     .                       freeze_xrbp,freeze_xni,itran_save,
     .                       itti_save,itte_save,itangrot_save,
     .                       init_nl,iodeflg1,nk,iangrot,jac_known)
                     R5_matsing = .true. !only one pass per call to Radua5
                                         !is possible
                     go to 1  
                CASE(-4)
                   write(ncrt,'("Matrix repeatedly singular")')
                   If( .not. R5_matsing )then
                     init_nl = 0
                     call set_freeze(freeze_type,itte,itti,itangrot,
     .                       itran,freeze_xte,freeze_xti,freeze_xwr,
     .                       freeze_xrbp,freeze_xni,itran_save,
     .                       itti_save,itte_save,itangrot_save,
     .                       init_nl,iodeflg1,nk,iangrot,jac_known)
                     R5_matsing = .true. !only one pass per call to Radua5
                                         !is possible
                     go to 1
                   endif
             END SELECT
             call STOP('Subroutine method_lines 2',0)     
           endif
         !check if starting guess for newton needs adjustment (on next pass):
         nstep_r5 = iwk_ode(16)
         naccpt_r5 = iwk_ode(17)
         nrjct_r5 = iwk_ode(18)
         if(nstep_r5 .gt. naccpt_r5 + nrjct_r5)iwk_ode(4) = 1
         init_nl = -1 !restore problem setup
         call set_freeze(freeze_type,itte,itti,itangrot,itran,
     .                freeze_xte,freeze_xti,freeze_xwr,
     .                freeze_xrbp,freeze_xni,itran_save,
     .                 itti_save,itte_save,itangrot_save,
     .                 init_nl,iodeflg1,nk,iangrot,jac_known)
         init_nl = 0

      else if(non_lin_method .eq. 2)then      !dlsoibt,= from odepack
         if( .not. allocated(wk_ode))then
           allocate (wk_ode(lwork),STAT = istat)
           if(istat .ne. 0)
     .        call allocate_error("wk_ode, method_lines2",0,istat)
           wk_ode(:) =0.0d0
           wk_ode(3) = 0.1  !frequency of calling jacobian( .1 less, .001 more)
         endif
         if( .not. allocated(iwk_ode))then
           allocate(iwk_ode(liwork),STAT = istat)
           if(istat .ne. 0)
     .        call allocate_error("iwk_ode, method_lines2",0,istat)
           iwk_ode(:) = 0.0d0
           iwk_ode(1) = nkt        !block size
           iwk_ode(2) = Min(nj-1,te_index-1)       !blocks in each direction
           if(iwk_ode(1)*iwk_ode(2) .ne. neq_ode)then
             print *,'block size,#blocks ',iwk_ode(1),iwk_ode(2)
            call STOP("SUB solve_ml2, block size incorrect",0)
           endif
         endif
         if( .not. allocated(ydot_ode))then
           allocate(ydot_ode(liwork),STAT = istat)
           if(istat .ne. 0)
     .        call allocate_error("ydot_ode, method_lines2",0,istat)
           ydot_ode(:) = 0.0d0
           istate = 0   !compute initial value of ydot_ode
         endif
           print *,'te_index,neq_ode =',te_index,neq_ode
         call dlsoibt(eval_ode3,addem,eval_ode_jacob,neq_ode,
     .                 y_ode,ydot_ode,time_ode,timeend_ode,itol,
     .                 tol_ode,atol,itask,istate,iopt3,wk_ode,lwork,
     .                iwk_ode,liwork,mf3)    
             if(istate  .lt. 0)then
             write(ncrt,'("dlsoibt reports erro,itask =",i5)')itask
             SELECT CASE(istate)
                CASE(-1)
                   write(ncrt,'("check inputs")')
                CASE(-2)
                   write(ncrt,'("tolerances too small")')             
                CASE(-3)
                   write(ncrt,'("illegal input")')   
                CASE(-4)
                   write(ncrt,'("check inputs")')
                CASE(-5)
                   write(ncrt,'("bad jac or tolerances")')
                CASE(-6)
                   write(ncrt,'("atol problem")')             
                CASE(-7)
                   write(ncrt,'("not supposed to happen")')   
                CASE(-8)
                   write(ncrt,'("Initial dydt not found")')
             END SELECT
             call STOP('Subroutine solve_ml2',0)     
           endif
      else if (non_lin_method .eq. 3) then   !Dgear, implemented only for te,ti
c
        miter_ode=3   ! use 0,2,3 with dummy eval_ode_jac in dgear
        nlc_ode=5
        nlu_ode=nlc_ode
        meth_ode=1
       if( .not. allocated(wk_ode))then
           allocate (wk_ode(kw_ode),STAT = istat)
           if(istat .ne. 0)
     .        call allocate_error("wk_ode, method_lines",0,istat)
        endif
        if( .not. allocated(iwk_ode))then
           allocate(iwk_ode(kk*kj),STAT = istat)
           if(istat .ne. 0)
     .        call allocate_error("iwk_ode, method_lines",0,istat)
        endif
        call my_dgear (neq_ode, eval_ode, eval_ode_jacob, time_ode,
     .            timestep_ode, y_ode, timeend_ode, tol_ode,
     .            meth_ode, miter_ode, index_ode, iwk_ode,
     .            wk_ode, ier_ode)
c


      else if(non_lin_method .eq. 4)then      !dlsodi  from odepack
         !method works with single profile but multiple
         !profiles (eq. te,ti,rbp simultaneously) is not yet
         !implemented(need to modify sub addem2,etc.). Use above 
         !methods for now. HSJ   
         if( .not. allocated(wk_ode))then
           allocate (wk_ode(lwork),STAT = istat)
           if(istat .ne. 0)
     .        call allocate_error("wk_ode, method_lines2",0,istat)
           wk_ode(:) =0.0d0
c           wk_ode(3) = 0.1  !frequency of calling jacobian( .1 less, .001 more)
         endif
         if( .not. allocated(iwk_ode))then
           allocate(iwk_ode(liwork),STAT = istat)
           if(istat .ne. 0)
     .        call allocate_error("iwk_ode, method_lines2",0,istat)
           iwk_ode(:) = 0.0d0
           iwk_ode(1) = ml        !block size
           iwk_ode(2) = mu        !blocks in each direction
         endif
         if( .not. allocated(ydot_ode))then
           allocate(ydot_ode(liwork),STAT = istat)
           if(istat .ne. 0)
     .        call allocate_error("ydot_ode, method_lines2",0,istat)
           ydot_ode(:) = 0.0d0
           istate = 0
         endif

         call dlsodi(eval_ode3,addem2,eval_ode_jacob,neq_ode,
     .                 y_ode,ydot_ode,time_ode,timeend_ode,itol,
     .                 tol_ode,atol,itask,istate,iopt3,wk_ode,lwork,
     .                iwk_ode,liwork,mf3)    
             if(istate  .lt. 0)then
             write(ncrt,'("dlsoibt reports erro,itask =",i5)')itask
             SELECT CASE(istate)
                CASE(-1)
                   write(ncrt,'("check inputs")')
                CASE(-2)
                   write(ncrt,'("tolerances too small")')             
                CASE(-3)
                   write(ncrt,'("illegal input")')   
                CASE(-4)
                   write(ncrt,'("check inputs")')
                CASE(-5)
                   write(ncrt,'("bad jac or tolerances")')
                CASE(-6)
                   write(ncrt,'("atol problem")')             
                CASE(-7)
                   write(ncrt,'("not supposed to happen")')   
                CASE(-8)
                   write(ncrt,'("Initial dydt not found")')
             END SELECT
             call STOP('Subroutine solve_ml2',0)     
           endif


      else

             call STOP('Sub solve_ml2 option not implemented',0)     

      end if

      timnew = timeend_ode
      freeze_te_index =.False.      !indicates that boundary location
      freeze_ti_index =.False.      !may now be changed if necessary
      freeze_rot_index=.False.
      R5_matsing = .False.
      return
c
      end


