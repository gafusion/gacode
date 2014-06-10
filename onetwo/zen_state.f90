    
     MODULE zen_state
         USE param,                                      ONLY :  kj,kion
         USE nrtype,                                     ONLY : I4B,DP

         USE numbrs,                                     ONLY :  nprim, nimp, nion,nj 
         USE ions,                                       ONLY :  namep,namei,nameb
         USE common_constants,                           ONLY :  izero,zeroc  
         REAL(DP),SAVE,dimension(:,:),ALLOCATABLE             :: zavold
         REAL(DP),DIMENSION(:,:),ALLOCATABLE                  :: z_l,zsq_l
         LOGICAL set_dzdt


         CONTAINS

      SUBROUTINE get_charge_state(te_l)
!------------------------------------------------------------------------------
! --- This subroutine sets z_l and zsq_l for impurities given the electron
! --- temperature,te_l over the rho grid
!------------------------------------------------------------------------------
         USE io,                                       ONLY : ncrt,nout
         USE numbrs,                                   ONLY : nj,nion
         USE solcon,                                   ONLY : te_range_check,timzen,   &
                                                              time0,timnew,timzen
         USE ions,                                     ONLY : dzdte,dzdtim,z,zsq


         IMPLICIT NONE
         INTEGER(I4B) i,k,j
         REAL(DP)  dtime,dte
         REAL(DP)  te0(nj),zav0(nj)                  ! temporary arrays
         REAL(DP)  te_l(nj)                          ! input array

         IF(.NOT. ALLOCATED(zavold)) ALLOCATE(zavold(nj,nimp))
         IF(.NOT. ALLOCATED(z_l))ALLOCATE(z_l(kj,kion))
         IF(.NOT. ALLOCATED(zsq_l))ALLOCATE(zsq_l(kj,kion))
         z_l(1:nj,1:nprim) = z(1:nj,1:nprim)
         zsq_l(1:nj,1:nprim) = zsq(1:nj,1:nprim)
         dtime = 0.0_DP
         if (timnew .ne. time0) dtime = timnew - timzen  ! dtime is local to zen
         te0(1:nj) = 0.99 * te_l(1:nj)

      do i=1,nimp
        k = nprim+i
        ! set impurity charge:
        call zfit   (zav0    , te0, namei(i), nj, ncrt, nout,te_range_check)
        call zfit   (z_l(1,k), te_l , namei(i), nj, ncrt, nout,te_range_check)
        call zsqfit (zsq_l(1,k), te_l , namei(i), nj, ncrt, nout, te_range_check)
!



     IF(set_dzdt)THEN 
        do j=1,nj
          dte = 0.01*te_l(j)
          z(j,k) = z_l(j,k)
          zsq(j,k) = zsq_l(j,k)
!
! --- dzdte = 0 (set in INIT) for primary ions, which are assumed
! --- to be fully ionized for all te.
! --- similarly for dzdtim, the value for the primary ions is always zero.
! --- this is implicit in the code for dzdtim and explicit for dzdte.
!


             dzdte(j,k) = (z_l(j,k)-zav0(j))/dte
             IF (dtime .NE.  zeroc)THEN 
                dzdtim(j,i) = (z_l(j,k) - zavold(j,i))/dtime
                if(ABS(dzdtim(j,i)) .gt. 1.d50)then
                   print *,'jzen ,i,k =',j,i,k
                   print *,'z(j,k),zavold(j,i) =',z_l(j,k),zavold(j,i)
                   call exit(1)
                endif
             ENDIF
             zavold(j,i) =  z_l(j,k)

        end do 
     ENDIF
      end do
 


      RETURN
      END SUBROUTINE get_charge_state




       SUBROUTINE en_calc(ene_l,zeff_l,en_l,ld_en,izenstop,set_u)
!----------------------------------------------------------------------------
! --- get primary and impurity densities given beam,alpha and electron densities
! --- and charge state information
! ---------------------------------------------------------------------------


!          inenez 1,-1 valid here
!          compute ion densities from ene and zeff
!

       USE nub2,                            ONLY :  enbeam
       USE fusion,                          ONLY :  enalp
       USE io,                              ONLY :  nout
       USE fusion,                          ONLY :  ifus, ihe,idt,id,it
       USE flags,                           ONLY :  itran,inenez,zfrac
       USE soln,                            ONLY :  u
       USE bd_condtn,                       ONLY :  ub,bc,enpb
       USE ions,                            ONLY :  adjzeff
       USE io,                              ONLY :  ncrt,nout
       USE verbose,                         ONLY :  zenvb
       USE solcon,                          ONLY :  time,time0,timzen

       IMPLICIT NONE
       INTEGER(I4B) k,ier,jer,j,i,izenstop,ntry,jl,ld_en
       REAL(DP) ene_l(kj),zeff_l(kj),en_l(ld_en,*)
       REAL(DP) dn,alpha_mult,zefftry,dzeff,entryp,entryi,zfracp
       LOGICAL set_u


      k   = nprim+1
!
      ier = IZERO
      jer = IZERO
      alpha_mult=1.0_DP
!
!     ne = z1*n1+z2*n2+zb*nb+za*na                with z1=zb=1 & za=2
!     zeff*ne = z1sq*n1+z2sq*n2+zbsq*nb+zasq*na   solve for n1 & n2:
!

 
      do 3200 j=1,nj
!
      if (ifus .ge. 0) then
!
          dn = z_l(j,1)*zsq_l(j,k)-z_l(j,k)*zsq_l(j,1)
          en_l(j,1) = (zsq_l(j,k)*(ene_l(j)-enbeam(j)-2.0*alpha_mult*enalp(j)) &
            -z_l(j,k)*(zeff_l(j)*ene_l(j)-enbeam(j)-4.0 *alpha_mult* &
                                                           enalp(j)))/dn

          en_l(j,k) = (ene_l(j)*(zeff_l(j)-1.0)-2.0*alpha_mult*enalp(j))/dn
          if (en_l(j,1) .le. 0. .or. en_l(j,k) .le. 0.0)  ier = 1
!
! --- try to correct for negative densities
! --- i.e., search for smallest zeff that will work
!
         if (ier .eq. 1) then
             izenstop = 1
             jer     = j
             ier     = 0
             ntry    = 10
             zefftry = zeff_l(jer)
             dzeff   = (zefftry-1.0)/(ntry-1)
             if (dzeff .le. 0.0)  go to 3040
            do jl=1,ntry
               zefftry = zefftry - dzeff
               zefftry = MAX (zefftry, 1._DP)
               entryp  = (zsq_l(j,k)*(ene_l(j)-enbeam(j)-2.0*enalp(j)) &
              -z_l(j,k)*(zefftry*ene_l(j)-enbeam(j)-4.0 * enalp(j)))/dn
               entryi  = (ene_l(j)*(zefftry-1.0)-2.0*enalp(j))/dn
               if ((entryp .gt. 0.0) .and. (entryi .gt. 0.0)) &
                                                         go to 3050
            end do
 3040       continue
!            write (nout, 9010) time,jer,ene_l(jer),enbeam(jer),
!     .                                 enalp(jer), zeff_l(jer)
!            write (nqik, 9010) time,jer,ene_l(jer),enbeam(jer),
!     .                                 enalp(jer), zeff_l(jer)
            write (ncrt, 9010) time,jer,ene_l(jer),enbeam(jer), &
                                       enalp(jer), zeff_l(jer)
 9010       format (' *************** ERROR: ZEFF IMPLIES NEGATIVE', &
                 ' ION DENSITY AT TIME ', f10.6,' *****************' / &
                 ' j, ene(j), enbeam(j), enapl(j), zeff(j) ='        / &
                   2x, i5, 4(1pe14.6, 2x))
!            write  (nout, 9005)
!            write  (nqik, 9005)
            write  (ncrt, 9005)
 9005       format (' even zeff = 1 will not work, ene too', &
                 ' small and/or enbeam, enalp too large' / &
                 ' ****** NO ADJUSTMENT TO ZEFF WAS MADE *****')

            write  (ncrt, 9018)  dn, en_l(j,1), en_l(j,k),k
 9018       format ('  dn      = ', 1pe14.4 / &
              '  en_(j,1) = ', 1pe14.4 / &
              '  en(j,k) = ', 1pe14.4 / &
              '  k       = ', i5)
!            write  (nout, 9019)  z_l(j,1),z_l(j,k),zsq_l(j,1),zsq_l(j,k)
!            write  (nqik, 9019)  z_l(j,1),z_l(j,k),zsq_l(j,1),zsq_l(j,k)
            write  (ncrt, 9019)  z_l(j,1),z_l(j,k),zsq_l(j,1),zsq_l(j,k)
 9019       format ('  z  (j,1) = ', 1pe14.4 / &
              '  z  (j,k) = ', 1pe14.4 / &
              '  zsq(j,1) = ', 1pe14.4 / &
              '  zsq(j,k) = ', 1pe14.4 )
            go to 3060
 3050       continue
!            write  (nout, 9015) time
 9015       format (' ******************* ERROR: ZEFF IMPLIES NEGATIVE', &
                 ' ION DENSITY AT TIME ', f10.6, ' *******************')
            if (adjzeff .eq. 1) then
!                write  (nout, 9020) jer,zeff_l(jer),zefftry,en_l(jer,1),
!     .                        entryp,en_l(jer,k),entryi
!                write  (nqik, 9020) jer,zeff_l(jer),zefftry,en_l(jer,1),
!     .                        entryp,en_l(jer,k),entryi
                write  (ncrt, 9020) jer,zeff_l(jer),zefftry,en_l(jer,1), &
                              entryp,en_l(jer,k),entryi
 9020          format ( &
            ' ********** WARNING: AT J = ', i5,' THE FOLLOWING' &
            ' CHANGES WERE MADE *****************************'        / &
            ' ZEFF (OLD) =  ', F10.4, '  ZEFF (NEW) =', f10.4         / &
            ' EN(1,J) (OLD) =', 1pe12.4, '  EN(1,J) (NEW) =', 1pe12.4 / &
            ' EN(1,K) (OLD) =', 1pe12.4, '  EN(1,K) (NEW) =', 1pe12.4)
              zeff_l(jer) = zefftry
              en_l(jer,1) = entryp
              en_l(jer,k) = entryi
              izenstop = 1 !HSJ 6/16/07   STOP  code per request from pr
           else
!               write (nout, 9025)  zefftry, jer, entryp, entryi
!               write (nqik, 9025)  zefftry, jer, entryp, entryi
               write (ncrt, 9025)  zefftry, jer, entryp, entryi
 9025          format (' must have Zeff < ', f10.4, ' at j = ', i5 / &
                       ' to get primary and impurity densitites of ', &
                         2(1pe12.4, 2x))
               izenstop = 1 !if density is negative we have to stop
           end if ! adjzeff branch
         end if   ! ier = 1 branch
 3060    if (nprim .eq. 1)  go to 3100
         en_l(j,2) = (1.0 - zfrac)*en_l(j,1)
         en_l(j,1) = zfrac*en_l(j,1)
!
      else        ! ifus < 0 branch
           zfracp = -1.0_DP
           if (time .eq. time0)  zfracp = zfrac

           call get_density (enbeam,enalp,ene_l,en_l,zsq_l,zeff_l,z_l, &
                             itran,kj,j,nprim,nimp,ihe,id,it, &
                             ncrt,nout,zfracp,time)

           if (time .eq. time0 .and. j .eq. nj) then
              if (0.0 .lt. zfrac .and. zfrac .lt. 1.0) then
                 do i=1,nprim
                     bc(1,i)   =  en_l(j,i)
                     enpb(1,i) =  bc(1,i)
                     ub(i)     =  bc(1,i)
                     u(i,j)    =  en_l(j,i)    ! set boundary condition
                 end do
              end if
           end if
      end if
 
 3100 IF(set_u)THEN
 
        do  i=1,nion
           if      (inenez .eq. -1) then
              u(i,j) = en_l(j,i)  ! density may be in simulation mode
           else if (itran(i) .eq. 0) then
              u(i,j) = en_l(j,i)  ! density in analysis mode
           else if (itran(i) .eq. 1 .and. time .eq. time0) then
              u(i,j) = en_l(j,i)
           else
              en_l(j,i)= u(i,j)   ! density in simulation mode HSJ 2/14/96
           end if
        ENDDO
      ENDIF

3200  continue ! end nested loops over spatial mesh j




       if (zenvb .gt. 0) then
          write (*,'(/ " on exit from ZEN we have nprim = ", i5)') nprim
          write (*,'(  " primary and impurity ion densities", &
                       " determined from ENE and ZEFF")')
          do i=1,nprim
             write (*, '(" i = ", i3, &
                         " en(1,i),en(nj,i) = ",2(1pe14.4))')i,en_l(1,i), &
                                                               en_l(nj,i)
          end do
          write (*, '(" on exit from ZEN we have nimp =", i5)') nimp
          if (nimp .gt. 0) then
            do i=1,nimp
               write (*,'(" i = ",i3, &
                " en(1,nprim+i),en(nj,nprim+i)  = ",2(1pe14.4))')i, &
                              en_l(1,nprim+i),en_l(nj,nprim+i)
            end do
          end if
          write (*,'(" enalp (1:5) =  ", 5(1pe14.4))') (enalp(j),j=1,5)
          write (*,'(" enbeam(1:5) =  ", 5(1pe14.4))') (enbeam(j),j=1,5)
          write (*,'(" ene(1), ene(nj) =  ",2(1pe14.4))') ene_l(1),ene_l(nj)
       end if



       DEALLOCATE(z_l,zsq_l)
       RETURN
       END SUBROUTINE en_calc


       END MODULE zen_state
