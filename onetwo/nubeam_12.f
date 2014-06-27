      subroutine  beam_power_interp(time)
c ----------------------------------------------------------------
c use data from ufiles (or from transp namelist file) to load
c powers, energy and beam fractions
c --------------------------------------------------------HSJ 12/04/03
      USE transp, only : beam_data

      IMPLICIT none
      INTEGER i,j,l,nt,nchan
      REAL *8 time
      REAL *8  dt1,dt2,slope,intc
 
      i=0
      nt = SIZE(beam_data%beam_times)
      do j=1,nt-1
         dt1 =  time - beam_data%beam_times(j)
         dt2 =  time - beam_data%beam_times(j+1)
         if(dt1*dt2 .le. 0.0)then
            i =j
            exit
         endif
      enddo


!       do j=1,SIZE(beam_data%beam_inject,DIM=1)
!          print*, 'beam_times,inject  ='
!          do l=1,beam_data%nbeam
!            print *,beam_data%beam_times(j),beam_data%beam_inject(j,l)
!          enddo
!       enddo


      if(i .eq.0 .and. time .gt. beam_data%beam_times(nt) )then
         beam_data%beam_times(nt) = time !extend due to roundoff
         i = nt-1
      elseif( i .eq. 0)then
         print *,'dt1,dt2 =',dt1,dt2
         print *,'time =',time
         print *,'beam_data%beam_times =',beam_data%beam_times
         call STOP('beam_power_interp',1)
      endif
c     linearly interpolate between (i) and (i+1)  values to get quantities:
       dt1 = beam_data%beam_times(i+1)-beam_data%beam_times(i)


       do j=1,beam_data%nbeam          !loop over beamlines

           slope = (beam_data%beam_inject(i+1,j) -
     &                    beam_data%beam_inject(i,j) )  / dt1
           intc = beam_data%beam_inject(i+1,j) - 
     &                   slope *beam_data%beam_times(i+1)
           beam_data%pinja(j)  =  slope *time +intc
           print *,'in interp, dt1,time,j ',dt1,time,j
           PRINT *,' beam_data%pinja(j)=', beam_data%pinja(j)
           print *,'injct i,i+1 =',beam_data%beam_inject(i,j),
     .            beam_data%beam_inject(i+1,j)

           slope = (beam_data%beam_inject(i+1,beam_data%nbeam+j) -
     &              beam_data%beam_inject(i,beam_data%nbeam+j) )/dt1
           intc = beam_data%beam_inject(i+1,beam_data%nbeam+j) - 
     &                   slope *beam_data%beam_times(i+1)
           beam_data%einja(j)  =  slope *time +intc 


           slope = (beam_data%beam_inject(i+1,2*beam_data%nbeam+j) -
     &               beam_data%beam_inject(i,2*beam_data%nbeam+j) )/dt1
           intc = beam_data%beam_inject(i+1,2*beam_data%nbeam+j) - 
     &                   slope *beam_data%beam_times(i+1)
           beam_data%ffulla(j) = slope *time +intc 


           slope = (beam_data%beam_inject(i+1,3*beam_data%nbeam+j) -
     &               beam_data%beam_inject(i,3*beam_data%nbeam+j) )/dt1
           intc = beam_data%beam_inject(i+1,3*beam_data%nbeam+j) - 
     &                   slope *beam_data%beam_times(i+1)
           beam_data%fhalfa(j) = slope *time +intc

       enddo
       write(173,fmt='("time,beam_data%nbeam",
     & x,1pe12.6,x,i3)')time,beam_data%nbeam
       write(173,fmt='("beam_data%pinja(j) =",5(1pe12.2,x))')
     & beam_data%pinja(1:beam_data%nbeam)  ! 888888999999
       !call stop('nubeam_12.f l 79',222)

      return
      end subroutine  beam_power_interp



