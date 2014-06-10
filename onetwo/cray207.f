

      subroutine areaclc (psisw, psinw, psine, psise, psibdry, psik,
     .                    drmhdgrd, dzmhdgrd, drhalf, dzhalf, drdzhalf,
     .                    iounit, darea, daream, psictr, rctr)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c --- subroutine calculates the area associated with filament currents
c --- near the plasma boundary. for filament (i,j) the coordinate system
c --- we work with here has its origin centered on the point (i,j).
c --- it is assumed that psi values less than psibdry represent points
c --- inside the plasma and psi values greater than psibdry represent
c --- points outside the plasma.
c --- input:
c  psisw         psi value, south west corner
c  psinw                    north west
c  psine                    north east
c  psise                    south east
c --- note that the above psi values are defined on the half grid!
c --- thus each of these psi values is the average of 4 psi values on
c --- the (rmhdgrid(i),zmhdgrid(j)) grid.
c  psik         psi at cell center (at rmhdgrid(i),zmhdgrid(j))
c  psibdry      boundary value of psi
c  drmhdgrd     spacing of mhdgrid in r
c  dzmhdgrd                           z
c  drhalf       drmhdgrd/2
c  dzhalf       dzmhdgrd/2
c  drdzhalf     drhalf*dzhalf
c  iounit       unit number for writting of error message in case
c               routine fails. if iounit = 0 no error message generated
c
c --- output
c  the output consists of 4 quantities: darea,daream,psictr and rctr.
c  if the point at which psi = psik is inside the plasma (i.e., psibdry
c  is greater than psik) then darea is returned as that part of the
c  area of the cell which is inside the plasma,and psictr is returned
c  as the value of psi corresponding to the geometrical center of the
c  region with this area. rctr is the radial distance of psictr from
c  psik. Note that if the point psik is inside the plasma and is
c  sufficiently far away from the boundary that the cell is
c  entirely inside tha plasma then darea is just the area of the cell,
c  psictr is just psik and rctr is returned as zero. (rctr is measured
c  from the pt at which psi = psik).
c  whenever the pt psi = psik is inside the plasma,daream is returned
c  as 0.0.
c
c  if the point at which psi = psik is outside the plasma then darea
c  is returned as 0.0 . In this case daream is returned as the area
c  of that part of the cell which is inside the plasma.
c  again rctr and psictr correspond to values at the geometric center
c  of that part of the cell which is inside the plasma boundary.
c  if both darea and daream are returned as 0.0 then the cell is
c  totally outside of the plasma. In this case psictr and rctr are not
c  used in subsequent calculations and hence are not set on return
c  from this subroutine.
c  NOTE that rctr and zctr (used below) are only rough approximations
c  to the desired geometric center of the area!
c  if any of the four psi values become close to identical then
c  the interpolation becomes ill conditioned. no check for this is made!
c ------------------------------------------------------------------ HSJ
c
      dimension  rpostn(2), zpostn(2)
c
      fdarea (r1,r2,a,b) = 0.5 * a * (r2 * r2 - r1 * r1) + b * (r2 - r1)
c
c     fdarea integrates z = a*r+b  w.r.t. r from r=r1 to r=r2
c
c --- first step: get the intersection of line with the cell
c
      areat  = drmhdgrd*dzmhdgrd
      rctr   = 0.0
      ipostn = 0
      if ((psibdry-psisw)*(psibdry-psise) .gt. 0.0)  go to 100
        dpsi   = psise-psisw
        a      = drmhdgrd/dpsi
        b      = drhalf-a*psise
        ipostn = 1
        rpostn(ipostn) = a*psibdry+b
        zpostn(ipostn) = -dzhalf
  100 if ((psibdry-psinw)*(psibdry-psisw) .gt. 0.0)  go to 200
        dpsi   = psinw-psisw
        a      = dzmhdgrd/dpsi
        b      = -dzhalf-a*psisw
        ipostn = ipostn+1
        rpostn(ipostn) = -drhalf
        zpostn(ipostn) = a*psibdry+b
        if (ipostn .eq. 2)  go to 700
  200 if ((psibdry-psinw)*(psibdry-psine) .gt. 0.0)  go to 300
        dpsi   = psine-psinw
        a      = drmhdgrd/dpsi
        b      = drhalf-a*psine
        ipostn = ipostn+1
        rpostn(ipostn) = a*psibdry+b
        zpostn(ipostn) = dzhalf
        if (ipostn .eq. 2)  go to 700
  300   if ((psibdry-psine)*(psibdry-psise) .gt. 0.0)  go to 600
        dpsi   = psine-psise
        a      = dzmhdgrd/dpsi
        b      = -dzhalf-a*psise
        ipostn = ipostn+1
        rpostn(ipostn) = drhalf
        zpostn(ipostn) = a*psibdry+b
        if (ipostn .eq. 2)  go to 700
c
c --- next get the cell area inside the plasma associated with the point
c --- case a: cell is totally inside or totally outside the plasma
c
  600 if (ipostn .lt. 2) then          ! didn't find pisbdry in cell
        daream = 0.0
        if (psik .le. psibdry) then    ! pt is  inside plasma
          darea  = areat
          psictr = psik
        else                           ! pt is outside plasma
          darea  = 0.0
        end if
        return
      end if
c
c --- case b:cell intersects plasma boundary
c --- area calc for vertical boundary
c
  700   if (ABS (rpostn(1)-rpostn(2)) .lt. 1.0e-06) then ! vertical bdry
          if (psik .lt. psibdry) then                    ! inside plasma
            daream = 0.0
            if (rpostn(1) .ge. 0.0) then
              darea = (rpostn(1)+drhalf)*dzmhdgrd
              rctr  = (-drhalf+rpostn(1)) * 0.5
            else
              darea = (drhalf-rpostn(1))*dzmhdgrd
              rctr  = (rpostn(1)+drhalf) * 0.5
            end if
          else                                          ! outside plasma
            darea = 0.0
            if (rpostn(1) .ge. 0.0) then
              daream = (drhalf-rpostn(1))*dzmhdgrd
              rctr   = (rpostn(1)+drhalf) * 0.5
            else
              daream = (rpostn(1)+drhalf)*dzmhdgrd
              rctr   = (rpostn(1)-drhalf) * 0.5
            end if
          end if
          zctr = 0.0
          go to 800
        end if
c
c --- area calc for horizontal boundary
c
        if (ABS (zpostn(1)-zpostn(2)) .lt. 1.0e-06) then ! horiz bdry
          if (psik .lt. psibdry) then                    ! inside plasma
            daream = 0.0
            if (zpostn(1) .ge. 0.0) then
              darea = (zpostn(1)+dzhalf)*drmhdgrd
              zctr  = (-dzhalf+zpostn(1)) * 0.5
            else
              darea = (dzhalf-zpostn(1))*drmhdgrd
              zctr  = (zpostn(1)+dzhalf) * 0.5
            end if
          else                                          ! outside plasma
            darea = 0.0
            if (zpostn(1) .ge. 0.0) then
              daream = (dzhalf-zpostn(1))*drmhdgrd
              zctr   = (zpostn(1)+dzhalf) * 0.5
            else
              daream = (zpostn(1)+dzhalf)*drmhdgrd
              zctr   = (zpostn(1)-dzhalf) * 0.5
            end if
          end if
          rctr = 0.0
          go to 800
        end if
c
c --- area calculation for inclined boundary
c
        if (rpostn(1) .gt. rpostn(2)) then
        r1 = rpostn(2)
        z1 = zpostn(2)
        r2 = rpostn(1)
        z2 = zpostn(1)
      else
        r1 = rpostn(1)
        z1 = zpostn(1)
        r2 = rpostn(2)
        z2 = zpostn(2)
      end if
      a = (z2-z1)/(r2-r1)
      b = z1-a*r1
c
c ----------------------------------------------------------------------
c --- given (r1,z1),(r2,z2),a and b,we now determine the associated
c --- area. Here (r1,z1) and (r2,z2) are the start and end of the line
c --- segement,a is the slope and b is the intercept of the line segement.
c ----------------------------------------------------------------------
c --- positive slope line
c
      if (a .gt. 0) then
        if (b .gt. 0.0) then                   ! positive intercept
          if (r1 .eq. -drhalf) then
            if (z1 .ge. 0.0) then              ! case a
              darea = fdarea(r1,r2,a,b)
              if (z2 .eq. dzhalf)  darea = darea+(drhalf-r2)*dzhalf
              darea = darea+2.0*drdzhalf
              if (psik .gt. psibdry) then      ! point is outside
                daream = areat-darea
                darea  = 0.0
                rctr   = (r1+r2) * 0.5
                zctr   = (dzhalf+z1) * 0.5
              else                             ! point is inside
                daream = 0.0
                rctr   = 0.0
                zb     = MIN (b,z2)
                zctr   = (zb-dzhalf) * 0.5
              end if
              go to 800
            else                               ! z1 < 0.0, case b
              r1a   = -b/a
              darea = fdarea(r1a,r2,a,b)
              if (z2 .eq. dzhalf)  darea = darea+(drhalf-r2)*dzhalf
              darea = darea+2.0*drdzhalf+fdarea(r1,r1a,a,b)
              if (psik .gt. psibdry) then      ! point is outside
                daream = areat-darea
                darea  = 0.0
                rctr   = (r1+r2) * 0.5
                zctr   = (dzhalf+z1) * 0.5
              else                             ! point is inside
                daream = 0.0
                rctr   = 0.0
                zb     = MIN (b,z2)
                zctr   = (zb-dzhalf) * 0.5
              end if
              go to 800
            end if
          else                ! r1 > -drhalf, z1 must be = -dzhalf,casec
            r1a   = -b/a
            darea = fdarea(r1a,r2,a,b)
            darea = darea+2.0*drdzhalf+(drhalf-r2)*dzhalf
            darea = darea+fdarea(r1,r1a,a,b)-(r1+drhalf)*dzhalf
              if (psik .gt. psibdry) then        ! point is outside
                daream = areat-darea
                darea  = 0.0
                rctr   = (-drhalf+r2) * 0.5
                zctr   = 0.0
              else                               ! point is inside
                daream = 0.0
                rctr   = (r1+drhalf) * 0.5
                zb     = MIN (b,z2)
                zctr   = (zb-dzhalf) * 0.5
              end if
              go to 800
          end if
        else                               ! b < 0
          if (r1 .eq. -drhalf) then        ! implies z2 < dzhalf, case D
            r1a   = -b/a
            darea = drdzhalf-fdarea(r1,r1a,a,b)
            darea = darea+drdzhalf-fdarea(r1a,r2,a,b)
              if (psik .gt. psibdry) then  ! point is outside
                daream = areat-darea
                darea  = 0.0
                rctr   = 0.0
                zctr   = (-dzhalf+z2) * 0.5
              else                         ! point is inside
                daream = 0.0
                rctr   = 0.0
                zctr   = (z1+dzhalf) * 0.5
              end if
              go to 800
          else                             ! r1 > -drhalf
            if (z2 .eq. dzhalf) then       ! case E
              r1a   = -b/a
              darea = drdzhalf+(r1+drhalf)*dzhalf
              darea = darea-fdarea(r1,r1a,a,b)
              darea = darea
     .              + drdzhalf-fdarea(r1a,r2,a,b)-(drhalf-r2)*dzhalf
              if (psik .gt. psibdry) then  ! point is outside
                daream = areat-darea
                darea  = 0.0
                rctr   = (drhalf+r1) * 0.5
                zctr   = 0.0
              else                         ! point is inside
                daream = 0.0
                rctr   = (-drhalf+r2) * 0.5
                zb     = MAX (b, -dzhalf)
                zctr   = (zb+dzhalf) * 0.5
              end if
              go to 800
            else if (z2 .ge. 0.0) then     ! case f
              r1a   = -b/a
              darea = drdzhalf+(r1+drhalf)*dzhalf-fdarea(r1,r1a,a,b)
              darea = darea+drdzhalf-fdarea(r1a,r2,a,b)
              if (psik .gt. psibdry) then  ! point is outside
                daream = areat-darea
                darea  = 0.0
                rctr   = (r1+r2) * 0.5
                zctr   = (-dzhalf+z2) * 0.5
              else                         ! point is inside
                daream = 0.0
                rctr   = 0.0
                zb     = MAX (-dzhalf,b)
                zctr   = (zb+dzhalf) * 0.5
              end if
              go to 800
            else                           ! z2 < 0.0, case g
              darea = 2.0 * drdzhalf - fdarea(r1,r2,a,b)
     .              + (r1 + drhalf) * dzhalf
              if (psik .gt. psibdry) then  ! point is outside
                daream = areat-darea
                darea  = 0.0
                rctr   = (r1+r2) * 0.5
                zctr   = (-dzhalf+z2) * 0.5
              else                         ! point is inside
                daream = 0.0
                rctr   = 0.0
                zb     = MAX (b,-dzhalf)
                zctr   = (zb+dzhalf) * 0.5
              end if
              go to 800
            end if
          end if                           ! r1 branch
        end if                             ! b < 0 branch
c
c --- second half, slope < 0
c
      else
        if (b .gt. 0.0) then
          if (r1 .eq. -drhalf) then
            if (z2 .ge. 0.0) then              ! case al
              darea = 2.0 * drdzhalf+fdarea(r1,r2,a,b)
              if (psik .gt. psibdry) then      ! point is outside
                daream = areat-darea
                darea  = 0.0
                rctr   = 0.0
                zctr   = (dzhalf+z2) * 0.5
              else                             ! point is inside
                daream = 0.0
                rctr   = 0.0
                zctr   = (b-dzhalf) * 0.5
              end if
              go to 800
            else                               ! z2 < 0, case BL
              r1a   = -b/a
              darea = fdarea(r1,r1a,a,b)+drdzhalf+fdarea(r1a,r2,a,b)
     .              + drdzhalf
              if (psik .gt. psibdry) then      ! point is outside
                daream = areat-darea
                darea  = 0.0
                rctr   = 0.0
                zctr   = (dzhalf+z2) * 0.5
              else                             ! point is inside
                daream = 0.0
                rctr   = 0.0
                zctr   = (b-dzhalf) * 0.5
              end if
              go to 800
            end if
          else                                 ! r1 > -drhalf
            if (z2 .ge. 0.0) then              ! case CL
              darea = (r1+drhalf)*dzhalf+fdarea(r1,r2,a,b)+2.0*drdzhalf
              if (psik .gt. psibdry) then      ! point is outside
                daream = areat-darea
                darea  = 0.0
                rctr   = (r2+r1) * 0.5
                zctr   = (z1+z2) * 0.5
              else                             ! point is inside
                daream = 0.0
                rctr   = 0.0
                zb     = MIN (b,dzhalf)
                zctr   = (zb-dzhalf) * 0.5
              end if
              go to 800
            else                               ! z2 < 0.0, case DL
              r1a = -b/a
              darea = drdzhalf+(r1+drhalf)*dzhalf+fdarea(r1,r1a,a,b)
     .              + drdzhalf+fdarea(r1a,r2,a,b)
              if (psik .gt. psibdry) then      ! point is outside
                daream = areat-darea
                darea  = 0.0
                rctr   = (r2+r1) * 0.5
                zctr   = (z1+z2) * 0.5
              else                             ! point is inside
                daream = 0.0
                rctr   = (r1a-drhalf) * 0.5
                zb     = MIN (b,dzhalf)
                zctr   = (zb-dzhalf) * 0.5
              end if
              go to 800
            end if
          end if
        else                                    ! intercept b < 0.0
          if (r1 .eq. -drhalf) then
            if (z1 .gt. 0.0) then               ! case EL
              r1a   = -b/a
              darea = drdzhalf-fdarea(r1a,r2,a,b)+(drhalf-r2)*dzhalf
     .              + drdzhalf-fdarea(r1,r1a,a,b)
              if (psik .gt. psibdry) then       ! point is outside
                daream = areat-darea
                darea  = 0.0
                rctr   = (r2+r1) * 0.5
                zctr   = (z1+z2) * 0.5
              else                              ! point is inside
                daream = 0.0
                rctr   = (r1a+drhalf) * 0.5
                zb     = MAX (b,-dzhalf)
                zctr   = (zb+dzhalf) * 0.5
              end if
              go to 800
              else                              ! z1 < 0.0
                if (r2 .eq. drhalf) then        ! case FL
                  darea = drdzhalf+drdzhalf-fdarea(r1,r2,a,b)
              if (psik .gt. psibdry) then       ! point is outside
                daream  = areat-darea
                darea   = 0.0
                rctr    = 0.0
                zctr    = (z1-dzhalf) * 0.5
              else                              ! point is inside
                daream  = 0.0
                rctr    = 0.0
                zctr    = (b+dzhalf) * 0.5
              end if
              go to 800
                else                            ! r2 < drhalf, case GL
                  darea = -fdarea(r1,r2,a,b)+(drhalf-r2)*dzhalf
     .                    +drdzhalf+drdzhalf
              if (psik .gt. psibdry) then       ! point is outside
                daream = areat-darea
                darea  = 0.0
                rctr   = (r2+r1) * 0.5
                zctr   = (z1+z2) * 0.5
              else                              ! point is inside
                daream = 0.0
                rctr   = 0.0
                zb     = MAX (b,-dzhalf)
                zctr   = (zb+dzhalf) * 0.5
              end if
              go to 800
                end if
              end if                            ! z1 branch
            else                                ! r1 > -drhalf
              if (r2 .lt. 0.0) then             ! case HL
                r1a = -b/a
                darea = 3.0 * drdzhalf - fdarea(r1,r1a,a,b)
     .                - (r1 + drhalf) * dzhalf
     .                - fdarea(r1a,r2,a,b) - r2 * dzhalf
              if (psik .gt. psibdry) then       ! point is outside
                daream = areat-darea
                darea = 0.0
                rctr = (r2-drhalf) * 0.5
                zctr = 0.0
              else                              ! point is inside
                daream = 0.0
                rctr = (r1a+drhalf) * 0.5
                zctr = 0.0
              end if
              go to 800
              else                              ! r2 > 0.0, case IL
                r1a = -b/a
                darea = drdzhalf-fdarea(r1,r1a,a,b)-(r1+drhalf)*dzhalf
     .           -fdarea(r1a,r2,a,b)+(drhalf-r2)*dzhalf+drdzhalf
              if (psik .gt. psibdry) then       ! point is outside
                daream = areat-darea
                darea = 0.0
                rctr = (r2-drhalf) * 0.5
                zctr = 0.0
              else                              ! point is inside
                daream = 0.0
                rctr = (r1a+drhalf) * 0.5
                zctr = (b+dzhalf) * 0.5
              end if
              go to 800
              end if
            end if                              ! r1 branch
        end if                                  ! b < 0 branch
      end if                                    ! a < 0 branch
c
c --- should not reach this point
c
      if (iounit .ne. 0)
     .write (iounit, 1000)  r1,z1,r2,z2,psisw,psinw,psine,psise,a,b
 1000 format ('  ERROR'                            /
     . '  r1,z1,r2,z2 =',4(2x,1pe12.6)             /
     . '  psisw,psinw,psine,psise =',4(2x,1pe12.6) /
     . '  a,b =',2(2x,1pe12.6))
c
  800 if (ABS (rctr) .gt. drhalf  .or.  ABS (zctr) .gt. dzhalf .or.
     .        daream .gt.  areat  .or.       darea .gt. areat)
     .  call STOP ('subroutine AREACLC: unspecified problem', 5)
      call areaintp (psisw,psinw,psine,psise,rctr,zctr,
     .               areat,drhalf,dzhalf,psictr)
      psictr = MIN (psictr, psibdry)
      return
c
      end

      subroutine areaintp (psisw, psinw, psine, psise, rp, zp,
     .                     areat, drhalf, dzhalf, psival)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      area1  = (drhalf-rp) * (zp+dzhalf)
      area2  = (drhalf-rp) * (dzhalf-zp)
      area3  = (rp+drhalf) * (dzhalf-zp)
      area4  = (rp+drhalf) * (zp+dzhalf)
      psival = (psinw*area1 + psisw*area2  +
     .          psise*area3 + psine*area4) / areat
      return
c
      end

      subroutine bound (psi, nw, nh, nwh, zero, x, y, xctr, yctr, rmin,
     .                  ix, limitr, xlim, ylim, nerr, limfag, radold,
     .                  npoint, xcontr, ycontr, ncontr, dpsi,
     .                  xmin, xmax, ymin, ymax, rymin, rymax,
     .                  zxmin, zxmax, psivl)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c***********************************************************************
c**                                                                   **
c**     main program:  MHD fitting code                               **
c**                                                                   **
c**     subprogram description:                                       **
c**          bound finds the outermost contour of a function          **
c**          given in the array psi. The contoured                    **
c**          function psi is assumed to be a decreasing function      **
c**          of distance away from the magnetic axis.                 **
c**              (but see switch nerr below)                          **
c**          bound uses bilinear interpolation rather than the cubic  **
c**          spline representation of psi used elsewhere.             **
c**          if ix<0, then BOUND just traces the field line           **
c**          starting at (xctr,yctr) either to closure or to the      **
c**          wall. ix = -1, trace clockwise. -2, counter clockwise.   **
c**   to find the plasma surface bound searches along                 **
c**   the line defined by (xctr,yctr),(rmin,yctr),                    **
c**   until it finds the largest surface that does not cross          **
c**   a separatrix or pass outside the limiter box defined            **
c**   by array zero. Hence rmin should be set to the minimum          **
c**   R(major) that the plasma boundary can achieve. This is normally **
c**   taken simply as the smallest R coordiante of the limiter minus, **
c**   say, 0.01 meter or one grid width, it doesn't matter.           **
c**   The benefit you gain for having precise estimates (that is      **
c**   xctr-rmin is small)is that the subroutine has to do fewer       **
c**   iterations to converge.                                         **
c**   dpsi is a flag returned by bound used to determine if the       **
c**   plasma was limited or diverted. Normally dpsi is passed         **
c**   on to other routines (eg xpoint) which then finds the           **
c**   location of the xpoint. The way this works is as follows.       **
c**   during the contour tracing done in bound the subroutine         **
c**   keeps track of whether or not the contour it traced             **
c**   passed through an x,y cell that also contained a segement       **
c**   of the limiter. if so then bound calculates the difference      **
c**   in psi between the value of psi on the plasma contour and       **
c**   the largest value of psi on the limiter segement in the         **
c**   cell and sets this equal to dpsi. thus if dpsi is "small"       **
c**   it signals that the plasma boundary was essentially             **
c**   touching (tangent to) the limiter segement. If dpsi is "large"  **
c**   then the plasma did not touch the limiter segement. If the      **
c**   plasma boundary is so far away from the limiter that the        **
c**   plasma contour and limiter segement never exist inside the      **
c**   same x,y cell (which implies a well diverted plasma)            **
c**   then dpsi is never calculated,its returned value will be        **
c**   1.0e10.                                                         **
c**   On output if dpsi .ne. 1.0e10 then the plasma boundary is close **
c**   enough to the limiter that they pass through the same x,y cell. **
c**   This indicates a marginally diverted/limited plasma. We         **
c**   may search for an x point under these circumstances (failure    **
c**   to find an xpoint because the plasma is limited is o.k in       **
c**   subroutine XPOINT). Note that the determination of whether or   **
c**   not the plasma is limited or diverted does not play a role in   **
c**   the calculations as such. The only important thing is that the  **
c**   boundary of the plasma is found. Hence the terms "small" and    **
c**  "large" used above should not be of paramount concern.           **
c**                                                                   **
c**   routine bound uses cellb,chkcrn,extrapbd,fqlin,maxpsi,minmax    **
c**                                                                   **
c**     input:                                                        **
c**       psi.............function to be contoured                    **
c**       nw..............row dimension of psi                        **
c**       nh..............column dimension of psi                     **
c**       nwh.............nw x nh                                     **
c**       zero............weighting array (found by subroutine ZLIM)  **
c**       x...............r grid                                      **
c**       y...............z grid                                      **
c**       xctr............r guess of contour center                   **
c**       yctr............z guess of contour center                   **
c**                           the (xctr,yctr) guess has to be an      **
c**                           interior point of the plasma. It is     **
c**                           usually sufficient to pick xctr,yctr    **
c**                           as the geometric center of the grid,    **
c**                           assuming that the plasma covers the     **
c**                           grid.                                   **
c**       rmin............see above                                   **
c**       ix..............flag explained above                        **
c**       limitr..........number of limiter points                    **
c**       xlim............r coordinates of limiter points             **
c**       ylim............z coordinates of limiter points             **
c**       nerr............flag used to indicate that psi is           **
c**                       negative on the magnetic axis and           **
c**                       increasing away from the magnetic axis      **
c**                       if nerr = -1 on input then                  **
c**                       psi is changed to -psi,                     **
c**                       the contour is generated,and then psi       **
c**                       is switched back to its original sign.      **
c**       limfag..........flag passed to subroutine zlim              **
c**                       limfag = 1 psi shape is convex                **
c**                       limfag = 2 psi shape is general               **
c**       radold..........guess of inside radius of plasma boundary   **
c**       npoint..........maximum number of contour points            **
c**                    which bound may find. that is,if bound         **
c**                    finds more points than can be stored in        **
c**                    xcontr,ycontr,an error exit (nerr = 3)           **
c**                    must be taken. There is no way to reduce       **
c**                    the number of points bound finds. The          **
c**                    number of points found will vary with the      **
c**                    nw,nh of the grid. if nw,nh are made large     **
c**                    the dimension of xcontr,ycontr,will have to    **
c**                    be increased. The number of points bound       **
c**                    finds is roughly proportional to 2(nw+nh).     **
c**                                                                   **
c**   output:                                                         **
c**       xcontr..........output r coordinates of contour             **
c**       ycontr..........output z coordinates of contour             **
c**       ncontr..........number of contour points found              **
c**       dpsi............an output flag described above              **
c**       xmin............minimum r value on contour                  **
c**       xmax............maximum r value                             **
c**       ymin............minimum z value                             **
c**       ymax............maximum z value                             **
c**                       and on limiter.  dpsi is used to            **
c**                       distinguish between a limited and a         **
c**                       diverted plasma                             **
c**       rymin...........r at ymin                                   **
c**       rymax...........r at ymax                                   **
c**       zxmin...........z at xmin                                   **
c**       zxmax...........z at xmax                                   **
c**       psivl...........psi value on the outermost contour          **
c**       dpsi............explained above                             **
c**                       if ix = -2,dpsi = 1.0e-06 is returned         **
c**       nerr............if an error occurred nerr .gt. 0            **
c**                       is returned                                 **
c**                                                                   **
c**     record of modification:                                       **
c**          26/04/83..........first created                          **
c**          01/09/83..........replaced contouring routines           **
c**                            with new dec10 version                 **
c**          12/03/84..........modifications to allow function as a   **
c**                            one field line tracer routine.         **
c**                            R Stambaugh.                           **
c**                                                                   **
c******************************************************************* HSJ
c
      dimension  psi(*), zero(*), x(*), y(*), xcontr(*), ycontr(*)
      dimension  dist(5), xlim(*), ylim(*)
      data       etolc, etol, nloop/1.0e-08, 1.0e-07, 120/
      data       mecopy/0/
c
c ----------------------------------------------------------------------
c nerr = -1, negative plasma current
c ----------------------------------------------------------------------
c
      nosign = 0
      if (nerr .eq. -1) then
        nosign = 1
        do i=1,nwh
          psi(i) = -psi(i)
        end do
      end if
c
      nerr  = 0
      loop  = 0
      psib0 = -1.0e10
      rad   = radold
      rin   = xctr
      rout  = rmin
c
c ----------------------------------------------------------------------
c field line tracing ix < 0
c ----------------------------------------------------------------------
c
      if (ix .lt. 0)  rad = xctr
      if (mecopy .gt. 0)  go to 20
      dx = x(2) - x(1)
      dy = y(2) - y(1)
      area = dx * dy
      rmid = 1.02 * (x(1)+x(nw)) / 2.0
      mecopy = 1
c
c ----------------------------------------------------------------------
c find starting value of psi
c ----------------------------------------------------------------------
c
   20 loop = loop+1
      i = 1+(rad-x(1))/dx
      if (rad-x(i) .lt. 0.0) i = i-1
      j = 1+(yctr-y(1))/(dy-0.000001)
      jjj = j
c
      if ((ix .eq. -2) .or. (ix .lt. -2 .and. rad .gt. rmid)) then
        j   = j - 1
        jjj = j + 1
      end if
c
      if ((ix .gt. 0) .and. (rad .gt. rmid)) then
        j = j-1
        jjj = j+1
      end if
      kstrt = (i-1)*nh+j
      kk = kstrt
      kold = kk
c
      if (ix .lt. 0 .and. ix .ge. -2)  go to 21
c
      dpsi = 1.0e10
      xmin = 1.0e10
      xmax = -1.0e10
      ymin = 1.0e10
      ymax = -1.0e10
   21 xt   = rad
      yt   = y(jjj)
      ncontr = 1
      xcontr(ncontr) = xt
      ycontr(ncontr) = yt
      a3 = (xt-x(i))*dy
      a4 = area-a3
      psivl = (psi(kk+nh+1)*a3+psi(kk+1)*a4)/area
      if (yt .eq. y(j))  psivl = psi(kk)+(psi(kk+nh)-psi(kk))
     .  *(xt-x(i))/dx
c
   30 f1 = psi(kk)
      f2 = psi(kk+nh)
      f3 = psi(kk+nh+1)
      f4 = psi(kk+1)
      x1 = x(i)
      x2 = x(i+1)
      y1 = y(j)
      y2 = y(j+1)
      if (ncontr .eq. 1)  go to 100
c
c ----------------------------------------------------------------------
c check for proximity to corner
c ----------------------------------------------------------------------
c
      dist(1) = (xt-x1)**2+(yt-y1)**2
      dist(2) = (xt-x2)**2+(yt-y1)**2
      dist(3) = (xt-x2)**2+(yt-y2)**2
      dist(4) = (xt-x1)**2+(yt-y2)**2
      dist(5) = MIN (dist(1),dist(2),dist(3),dist(4))
      if (dist(5) .gt. etolc)  go to 100
      do l=1,4
        kj = l
        if (dist(l) .eq. dist(5))  go to 50
      end do
c
c ----------------------------------------------------------------------
c kj points to appropriate corner
c ----------------------------------------------------------------------
c
   50 call chkcrn(psi,nwh,psivl,kold,knew,kj,kk,nh,i,i1)
      kk = knew
      i  = i1
      j  = kk-(i-1)*nh
   90 f1 = psi(kk)
      f2 = psi(kk+nh)
      f3 = psi(kk+nh+1)
      f4 = psi(kk+1)
      x1 = x(i)
      x2 = x(i+1)
      y1 = y(j)
      y2 = y(j+1)
c
c ----------------------------------------------------------------------
c check for limiter in cell
c ----------------------------------------------------------------------
c
  100 zsum = zero(kk)+zero(kk+1)+zero(kk+nh)+zero(kk+nh+1)
      if (zsum .eq. 0.0)  go to 1005
      if (ABS (zsum-4.0) .lt. 1.0e-03)  go to 540
c
c ----------------------------------------------------------------------
c from one to three corners of cell are inside limiter.  get max
c psi on line segment of limiter in cell and compare this max
c with current value of psilim (or psisep)
c note: do loop index assumes point 'limitr+1' is the same as
c point 'limitr'
c ----------------------------------------------------------------------
c
      psilx = -1.0e10
c
      do 520 k=1,limitr-1
        xc1 = xlim(k)
        yc1 = ylim(k)
        xc2 = xlim(k+1)
        yc2 = ylim(k+1)
        ik1 = 1+(xlim(k)-x(1))/dx
        jk1 = 1+(ylim(k)-y(1))/dy
        ik2 = 1+(xlim(k+1)-x(1))/dx
        jk2 = 1+(ylim(k+1)-x(1))/dy
        kij1 = (ik1-1)*nh+jk1
        kij2 = (ik2-1)*nh+jk2
        if ((kij1 .eq. kk) .and. (kij2 .eq. kk))  go to 510
c
c ----------------------------------------------------------------------
c       at least one limiter point is not in cell.  subroutine CELLB
c       returns intersections of cell boundaries and line defined by
c       points k and k+1 or one cell boundary and one interior point.
c       ifail=1 if points k and k+1 do not intersect the current cell
c       of interest.
c ----------------------------------------------------------------------
c
        ifail = 0
        call cellb(xc1,yc1,xc2,yc2,x1,y1,x2,y2,ifail)
        if (ifail .eq. 1)  go to 520
  510   call maxpsi (xc1,yc1,xc2,yc2,x1,y1,x2,y2,f1,f2,f3,f4,
     .               area,psilm,xtry1,ytry1)
        psilx = MAX (psilm, psilx)
        if (psilx .gt. psilm)  go to 520
        xtry  = xtry1
        ytry  = ytry1
  520 continue
c
      if (psilx .eq. -1.0e10)  go to 1090
      dpsi = MIN (dpsi, ABS (psivl-psilx))
      if (psilx-psivl)  540, 540, 530
  530 call zlim (zerol, 1, 1, limitr, xlim, ylim, xt, yt, limfag)
      if (zerol .le. 0.01)  go to 1005
  540 call extrapbd (f1,f2,f3,f4,x1,y1,x2,y2,xt,yt,xt1,yt1,xt2,yt2,
     .               psivl,area,dx,dy)
c
c ----------------------------------------------------------------------
c decide which intersection (xt1,yt1) or (xt2,yt2) is required
c ----------------------------------------------------------------------
c
      dist1 = (yt1-yt)**2+(xt1-xt)**2
      dist2 = (yt2-yt)**2+(xt2-xt)**2
      if (dist1 .lt. dist2)  go to 560
      yt = yt1
      xt = xt1
      go to 570
  560 yt = yt2
      xt = xt2
  570 ncontr = ncontr+1
      if (ncontr .gt. npoint)  go to 1090
      xcontr(ncontr) = xt
      ycontr(ncontr) = yt
c
c ----------------------------------------------------------------------
c find next cell
c ----------------------------------------------------------------------
c
  600 if (xt .eq. x2)  i = i+1
      if (xt .eq. x1)  i = i-1
      if (yt .eq. y2)  j = j+1
      if (yt .eq. y1)  j = j-1
c3/12/84
      if (ix .lt. 0 .and. ix .ge. -2)  go to 601
c3/12/84
      if (yt .lt. ymin)  rymin = xt
      if (yt .gt. ymax)  rymax = xt
      if (xt .lt. xmin)  zxmin = yt
      if (xt .gt. xmax)  zxmax = yt
      xmin = MIN (xmin,xt)
      xmax = MAX (xmax,xt)
      ymin = MIN (yt,ymin)
      ymax = MAX (yt,ymax)
  601 kold = kk
c
c ----------------------------------------------------------------------
c find new cell index
c ----------------------------------------------------------------------
c
      kk = (i-1)*nh+j
      if (kk .eq. kstrt)  go to 1040
      dis2p = SQRT ((xcontr(1)-xt)**2+(ycontr(1)-yt)**2)
      if ((dis2p .lt. 0.1*dx) .and. (ncontr .gt. 5))  go to 1040
      go to 30
c
c ----------------------------------------------------------------------
c psi on boundary smaller than psi on limiter, decrease rad and try again
c ----------------------------------------------------------------------
c
 1005 psib0 = psivl
c
      if (ix .lt. 0)  go to 2000
c
c ----------------------------------------------------------------------
c --- new code added 10/9/92
c --- for roundoff in limiter calculations
c
      if (ncontr .eq. 1) then
        rout = rin + 1.0e-9
        rin  = rin + 0.001
      end if
c
      if (loop .gt. nloop)  go to 1080
      rout = rad
      rad = (rin+rout) * 0.5
      go to 20
c
c ----------------------------------------------------------------------
c check for convergence of boundary
c converged in psi or in radius, whichever comes first
c ----------------------------------------------------------------------
c
 1040 err=ABS ((psivl-psib0)/psivl)
      if (ABS (rin-rout) .lt. 1.0e-6)  err=etol     ! new, 10/9/92 HSJ
c
        if (ix .ge. 0)  go to 1045
        if (ix .lt. -2) then
          dpsi = 1.0e-06
          go to 2000
        end if
        go to 2000
c
 1045 if (err .le. etol)  go to 1080
      if (loop .gt. nloop)  go to 1080
c
c ----------------------------------------------------------------------
c new rad, psi and try again
c ----------------------------------------------------------------------
c
      psib0 = psivl
      call zlim(zerol,1,1,limitr,xlim,ylim,rad,yctr,limfag)
      if (zerol .le. 0.01) then
        rout = rad
      else
        rin = rad
      end if
      rad = (rin+rout) * 0.5
      go to 20
 1080 radold = rad
      psib0 = psivl
      if (ABS (ycontr(1)-ycontr(ncontr)) .gt. 0.5 * dy)  go to 1090
      if (ABS (xcontr(1)-xcontr(ncontr)) .gt. 0.5 * dx)  go to 1090
      go to 2000
 1090 nerr = 3
c
 2000 if (nosign .eq. 1) then
        do i=1,nwh
          psi(i) = -psi(i)
        end do
        psivl = -psivl
      end if
      return
c
      end

      subroutine calc_torflux (torflux, psivolp, fpsi,
     .                         ratave, npsi, rmajor)
c
      USE constnts
      implicit none
c
c ----------------------------------------------------------------------
c --- subroutine calculates the toroidal flux without using dv/dpsi.
c --- instead we integrate dv/dtorflux, see eq. 2.2-4 gaa
c
c --- input
c  NOTE the following vectors are defined on the psival grid.
c  the convention is that psival(1) = plasma edge
c  and psival(npsi) = magnetic axis
c  psivolp(j)   j = 1,2,..npsi  volume
c  fpsi(j)      j = 1,2,..npsi  f(psi)
c  ratave(j)    j = 1,2,..npsi <Rmajor**2/R**2>
c  npsi         number of psi values in psival grid
c  rmajor
c
c --- output
c  torflux(j)   j = 1,2..npsi,toroidal flux on psival grid
c
c ------------------------------------------------------------------ HSJ
c
      integer j,npsi
      real*8  torflux(*),psivolp(*),fpsi(*),ratave(*),fctr,
     .        rmajor,dv,rmsq
c
c      include 'constnts.i'
c
      rmsq          = rmajor*rmajor * pi * 4.0
      torflux(npsi) = 0.0                     ! magnetic axis
      do j=npsi-1,1,-1                        ! from axis to plasma edge
        dv         = psivolp(j)-psivolp(j+1)  ! dv is positive
        fctr       = fpsi(j+1)*ratave(j+1)+fpsi(j)*ratave(j)
        torflux(j) = torflux(j+1)+fctr*dv/rmsq
      end do
      return
c
      end

      subroutine cellb (xc1, yc1, xc2, yc2, x1, y1, x2, y2, ifail)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c***********************************************************************
c**                                                                   **
c**     MAIN PROGRAM:  MHD FITTING CODE                               **
c**                                                                   **
c**     SUBPROGRAM DESCRIPTION:                                       **
c**          cellb redefines (xc1,yc1) and/or (xc2,yc2) so that       **
c**          they are intersections of cell boundaries unless         **
c**          one of the points is an interior point in which          **
c**          case it is not disturbed.cellb is called only by bound.  **
c**                                                                   **
c**     RECORD OF MODIFICATION:                                       **
c**          26/04/83..........first created                          **
c**                                                                   **
c******************************************************************* HSJ
c
      if ((xc1 .lt. x1) .and. (xc2 .lt. x1))  go to 1000
      if ((xc1 .gt. x2) .and. (xc2 .gt. x2))  go to 1000
      if ((yc1 .lt. y1) .and. (yc2 .lt. y1))  go to 1000
      if ((yc1 .gt. y2) .and. (yc2 .gt. y2))  go to 1000
      dx = xc1-xc2
      dy = yc1-yc2
      if ((dx .eq. 0.0) .and. (dy .eq. 0.0))  ifail = 1
      if (ifail .eq. 1  )  return
      if (   dx .eq. 0.0)  go to 200
      if (   dy .eq. 0.0)  go to 500
c
c ----------------------------------------------------------------------
c line is inclined. get equation y = alpha*x+beta.
c ----------------------------------------------------------------------
c
      alpha = dy/dx
      beta = yc1-alpha*xc1
      y1b = alpha*x1+beta
      if (y1b .lt. y1)  go to 20
      if (y1b .gt. y2)  go to 30
      xs1 = x1
      ys1 = y1b
      go to 40
   20 x1b = (y1-beta)/alpha
      if ((x1b .lt. x1) .or. (x1b .gt. x2))  go to 1000
      xs1 = x1b
      ys1 = y1
      go to 40
   30 x1b = (y2-beta)/alpha
      if ((x1b .lt. x1) .or. (x1b .gt. x2))  go to 1000
      xs1 = x1b
      ys1 = y2
   40 y2b = alpha*x2+beta
      if (y2b .lt. y1)  go to 50
      if (y2b .gt. y2)  go to 60
      xs2 = x2
      ys2 = y2b
      go to 70
   50 x2b = (y1-beta)/alpha
      if ((x2b .gt. x2) .or. (x2b .lt. x1))  go to 1000
      xs2 = x2b
      ys2 = y1
      go to 70
   60 x2b = (y2-beta)/alpha
      if ((x2b .gt. x2) .or. (x2b .lt. x1))  go to 1000
      xs2 = x2b
      ys2 = y2
c
c ----------------------------------------------------------------------
c at this point we have intersection (xs1,ys1) and (xs2,ys2)
c check for interior point
c ----------------------------------------------------------------------
c
   70 if ((x1 .lt. xc1) .and. (xc1 .lt. x2))  go to 80
      if ((x1 .lt. xc2) .and. (xc2 .lt. x2))  go to 90
      go to 100
   80 if ((y1 .lt. yc1) .and. (yc1 .lt. y2))  go to 160
      go to 100
   90 if ((y1 .lt. yc2) .and. (yc2 .lt. y2))  go to 170
  100 xc1 = xs1
      yc1 = ys1
      xc2 = xs2
      yc2 = ys2
      return
c
c ----------------------------------------------------------------------
c point (xc1,yc1) is interior point
c ----------------------------------------------------------------------
c
  160 if (yc2 .gt. yc1)  go to 165
      xc2 = xs1
      yc2 = ys1
      if (ys1 .gt. ys2)  yc2 = ys2
      if (ys1 .gt. ys2)  xc2 = xs2
      return
  165 xc2 = xs2
      yc2 = ys2
      if (ys1 .gt. ys2)  yc2 = ys1
      if (ys1 .gt. ys2)  xc2 = xs1
      return
c
c ----------------------------------------------------------------------
c point (xc2,yc2) is interior point
c ----------------------------------------------------------------------
c
  170 if (yc1 .gt. yc2)  go to 190
      xc1 = xs1
      yc1 = ys1
      if (ys1 .gt. ys2)  yc1 = ys2
      if (ys1 .gt. ys2)  xc1 = xs1
      return
  190 xc1 = xs2
      yc1 = ys2
      if (ys1 .gt. ys2)  yc1 = ys1
      if (ys1 .gt. ys2)  xc1 = xs1
      return
c
c ----------------------------------------------------------------------
c check if intersection exists for vertical line
c ----------------------------------------------------------------------
c
  200 if ((xc1 .lt. x1) .or. (xc1 .gt. x2))  go to 1000
c
c ----------------------------------------------------------------------
c is there an interior point?
c ----------------------------------------------------------------------
c
      if ((y1 .le. yc1) .and. (yc1 .le. y2))  go to 300
      if ((y1 .le. yc2) .and. (yc2 .le. y2))  go to 400
c
c ----------------------------------------------------------------------
c no interior points
c ----------------------------------------------------------------------
c
      yc1 = y1
      yc2 = y2
      return
c
c ----------------------------------------------------------------------
c point (xc1,yc1) is interior point
c ----------------------------------------------------------------------
c
  300 if (yc2 .gt. yc1)  yc2 = y2
      if (yc2 .lt. yc1)  yc2 = y1
      return
c
c ----------------------------------------------------------------------
c point  (xc2,yc2) is interior point
c ----------------------------------------------------------------------
c
  400 if (yc2 .gt. yc1)  yc1 = y1
      if (yc2 .lt. yc1)  yc1 = y2
      return
c
c ----------------------------------------------------------------------
c check if intersection exists for horizontal line
c ----------------------------------------------------------------------
c
  500 if ((yc1 .lt. y1) .or. (yc1 .gt. y2))  go to 1000
c
c ----------------------------------------------------------------------
c is there an interior point?
c ----------------------------------------------------------------------
c
      if ((x1 .lt. xc1) .and. (xc1 .le. x2))  go to 550
      if ((x1 .le. xc2) .and. (xc2 .le. x2))  go to 600
c
c ----------------------------------------------------------------------
c no interior points
c ----------------------------------------------------------------------
c
      xc1 = x1
      xc2 = x2
      return
c
c ----------------------------------------------------------------------
c point (xc1,yc1) is interior point
c ----------------------------------------------------------------------
c
  550 if (xc1 .lt. xc2)  xc2 = x2
      if (xc1 .gt. xc2)  xc2 = x1
      return
c
c ----------------------------------------------------------------------
c point (xc2,yc2) is interior point
c ----------------------------------------------------------------------
c
  600 if (xc2 .gt. xc1)  xc1 = x1
      if (xc2 .lt. xc1)  xc1 = x2
      return
 1000 ifail = 1
      return
c
      end

      subroutine chekeqd (p, nyeqd, xdimeqd, ydimeqd, redeqd, iounit,
     .                    rcontr, zcontr, ncontr, beqd, xaxis, yaxis,
     .                    psimag, psilim, toteqd, reqd,nconmax)
c
c
c ----------------------------------------------------------------------
c
c     We must have a reference point to start from. Unfortunately not
c     all eqdsks are internally self-consistent (some of the early
c     diii-d eqdsks have the wrong ordering,and or wrong sign in some
c     of the required quantities). The object of this subroutine is
c     twofold. First we must make sure that the eqdsk is self consistent
c     and second,we want to define the equilibrium relative to a right
c     hand coordinate system,(R,Z,phi),which is used in transport
c     calculations. The transport formalism embedded in ONETWO assumes
c     that Jphi and btor (or j parallel and b) are in the positive phi
c     direction. Furthermore the poloidal flux function psi is assumed
c     to be equal to (minus poloidal flux)/(2.0*pi),with the poloidal flux
c     defined from the symmetry (z) axis out to some flux surface,(not
c     from the magnetic axis out to some flux surface).
c     In order to accomodate these conditions we must define the positive
c     phi direction in the direction of the current flow. Thus the
c     (R,Z,phi) coordinate system used has the z axis pointing downward
c     for current flow in the counter clockwise direction (viewed from
c     above). This condition is evidenced in the eqdsk by the fact that
c     toteqd is greater than zero. Similarly the z axis points upward
c     for current flow in the clockwise direction (viewed from above).
c     This condition is evidenced by the fact that toteqd is less than
c     zero. Because the current and b field are most often oppositely
c     directed (for dIII-d) we must retain the sign of the toroidal
c     field in the calculations. The directional sense of the b field
c     is given by the sign of beqd in the eqdsk file. Below we change
c     these signs as appropriate and set some flags so that the
c     boundary conditions for the MHD equations can be adjusted in sign
c     as appropriate in another subroutine.
c     The Grad-Shafranov equation used in ONETWO (consistent with
c     the definition of psi given above) is
c
c        delstar psi = (4pi/c)*R*Jphi,with Jphi = -c*R*pprim-(c/4pi)*ffprim/R
c        ( delstar psi = u0*R*Jphi,with Jphi = -R*pprim-ffprim/(u0*R))
c
c        In order to check the eqdsk we have to assume
c     that toteqd and beqd are properly referenced to the phi direction
c     of the diii-d coordinate system (which is (R,phi,Z) with phi
c     positive in the counter clockwise direction,viewed from above).
c     (i.e., toteqd, btor > 0.0 (ON INPUT) means the toroidal current and b
c     field point in the counterclockwise direction and toteqd, btor < 0
c     (ON INPUT) means these quantities point in the clockwise direction
c     (both when viewed from above). Generally beqd and toteqd have
c     opposite signs for DIII-D.
c     ON OUTPUT toteqd is positive (i.e., current flows in positive phi
c     direction of transport grid) and beqd is referenced relative to
c     this direction of toteqd. Note that it is only the coordinate
c     system that is redefined. The actual current flow (clockwise or
c     counter clockwise ) is not changed.
c
c ---   input (argument list):
c           p             psi from eqdsk (in gauss-cm**2)
c           nxeqd
c           nyeqd         size of grid on eqdsk
c           xdimeqd
c           ydimeqd
c           redeqd        define the eqdsk MHD grid  (in cm)
c           iounit        for error messages, set to zero to suppress
c           ncontr
c           rcontr
c           zcontr        define the plasma boundary  (in cm)
c           nconmax       dimension of rcontr,zcontr
c           beqd          btor at reqd in eqdsk (in gauss)
c           reqd          defines location of btor value (in cm)
c           xaxis
c           yaxis         define magnetic axis of eqdsk equilibrium
c           psimag
c           psilim        axis and edge psi values
c           toteqd        total toroidal current in eqdsk (in amps)
c --- input from INCLUDE files:
c           pppsi     - pressure derivative array
c                           (in units of gram/(gauss*cm**3*sec**2)
c           ffppsi    - vacuum toroidal field modification array (gauss)
c           xdum,ydum,zdum,wdum,rdum   temporary storage vectors
c           cspln,wnoperm,pds     bicubic spline work area
c           psival
c           presspsi       erg/cm**3 (or dyne/cm**2)
c           qpsi
c
c --- output:
c          p(i,j)
c          beqd
c          toteqd
c          psimag
c          psilim
c          pppsi
c          ffppsi
c          presspsi
c          qpsi
c          mhdmultp =  1 if (R,Z,phi) with Z axis downward
c                   = -1                   Z      upward
c --- the output has the same numerical values as the input. The only
c --- changes are possibly in the signs and ordering of the vectors
c --- (if they were incorrectly ordered on input). The output is
c --- consistent with the transport convention discussed above.
c --- q and beqd can be negative on output!
c
c ------------------------------------------------------------------ HSJ
c

      USE mhdpar        
      USE constnts
      USE psig
      USE bicube
      USE mhdbcdtn
      USE io,                                  ONLY : nitre
      USE replace_imsl,                        ONLY : my_usmnmx,
     .                                                my_ibcccu,
     .                                                my_dbcevl1

      implicit  integer (i-n), real*8 (a-h, o-z)

      include 'imsl.i'
      include 'storage.i'

c
      dimension p(nw,nh),rcontr(nconmax),zcontr(nconmax)
c
c --- set mhdmultp = 1 for clockwise current, = -1 for counter clockwise
c --- current then set sign of beqd relative to toteqd
c
      mhdmultp = 1
      if (toteqd .lt. 0.0) then
        toteqd = -toteqd
        mhdmultp = -1
      end if
c      beqd = beqd*mhdmultp       ! 88889999 disabled HSJ 1/31/13 
c
c --- first get the grid that corresponds to the eqdsk
c
      if (nyeqd .gt. nh)  go to 1000
      if (nxeqd .gt. nw)  go to 1000
      dxneqd = xdimeqd/(nxeqd-1)
      dyneqd = ydimeqd/(nyeqd-1)
      do 10 i=1,nxeqd
   10 xdum(i) = redeqd+(i-1)*dxneqd
      do 20 j=1,nyeqd
   20 ydum(j) = -ydimeqd * 0.5+(j-1)*dyneqd
c
c --- get psi and gradient on magnetic axis
c
      icalc = 6
      imslmd='1331c207'

      call my_ibcccu (p,xdum,nxeqd,ydum,nyeqd,cspln,nw,wnoperm,ier)
      call my_dbcevl1 (xdum,nxeqd,ydum,nyeqd,cspln,nw,xaxis,yaxis,pds,
     .ier,icalc)
      psimag = pds(1)    ! get exact psimag using bicubic representation
      if (pds(5)*pds(6) .lt. 0.0)  go to 3000
      if (pds(5)        .lt. 0.0) then     ! psi is max on axis, reverse
          do 40 j=1,nyeqd
          do 40 i=1,nxeqd
   40       p(i,j) = -p(i,j)
      end if
c
c --- the four lines below were commented out on 18 Jan 95 by YRLL and HSJ
c
****  if (psimag .gt. 0.0) then
****     psimag = -psimag
****     psilim = -psilim
****  end if
c
c --- get the r,z extremes of the plasma boundary contour
c
      if (ncontr .gt. 0) then
        call my_usmnmx (rcontr, ncontr, 1, xplmin, xplmax)
        call my_usmnmx (zcontr, ncontr, 1, yplmin, yplmax)
      else
c
c --- if plasma boundary was not read in (because some old
c --- eqdsks and eqdsks from other codes didn't write the
c --- complete eqdsk) then get the plasma boundary now:
c
c           need p positive n axis for CNTOUR
c
            cconst = -1.0
            call multpl1 (p, nxeqd*nyeqd, cconst)
c
c           fit bicubic spline
c
            imslmd ='1367c207'
            call my_ibcccu  (p,xdum,nxeqd,ydum,nyeqd,cspln,nw,
     .                       wnoperm,ier)
            call multpl1 (p,nxeqd*nyeqd,cconst)    ! restore p
            psivl  = -psilim
            delta_psi=(-psimag+psilim)/nxeqd
            arcl   = 2.00   ! 2.0 cm arc length
            bperr  = 0.05   ! 5% rel. change in bp along contour allowed
            bperr  = 1.e-3  ! 999,
            dx     = xdum(2)-xdum(1)
            dy     = ydum(2)-ydum(1)
            iauto  = 1
            iconvg = 1
            xmin   = xdum(2)    ! define search box for CNTOUR
            xmax   = xdum(nxeqd-1)
            ymax   = ydum(nyeqd-1)
            ymin   = ydum(1)
c
            call cntour (xaxis,yaxis,psivl,xplmin,xplmax,yplmin,
     .                   yplmax,yxmin,yxmax,xymin,xymax,dang,arcl,
     .                   bperr,dx,dy,xmin,xmax,ymin,ymax,iauto,iautoc,
     .                   rcontr,zcontr,ncontr,xdum,nxeqd,ydum,nyeqd,
     .                   cspln,n2cspln,nh2,iounit,nconmax,ierr,
     .                   rdum,iconvg,delta_psi)
c
            if (ncontr .lt. 25)  go to 2000
            if (  ierr .ne.  0)  go to 2000
            if (iautoc .eq.  1)  psilim = -psivl
      end if
c
c --- next construct a temporary psi grid,psival,such that psival(1) =
c --- plasma edge psi,psival(nxeqd) = mag axis psi
c
      iunfrm = 1
      call psiset (nxeqd, psimag, psilim, psival, iunfrm)
c
c --- check the pressure,presspsi. Here we simply assume that kinetic
c --- pressure on axis is always greater than it is at the plasma edge
c
      if (presspsi(nxeqd) .lt. presspsi(1))
     .  call revers(presspsi,zdum,nxeqd)
c
c --- differentiate the pressure
c
      call difydx (psival, presspsi, zdum, nxeqd)
c
c --- now compare zdum with eqdsk pressure profile gradient,pppsi.
c --- this and other tests below are fuzzy. We cannot look for an
c --- exact match. Hence the indirect approach used below is necessary.
c --- The choices are: either zdum(i) agrees with +/- pppsi(i) or
c --- zdum(i) agrees with +/- pppsi(nxeqd-i+1).
c
      sum1 = 0.0
      sum2 = 0.0
      sum3 = 0.0
      sum4 = 0.0
      do i=2,nxeqd-1
        sum1 = sum1+(zdum(i)-pppsi(i))**2
        sum2 = sum2+(zdum(i)-pppsi(nxeqd-i+1))**2
        sum3 = sum3+(zdum(i)+pppsi(i))**2
        sum4 = sum4+(zdum(i)+pppsi(nxeqd-i+1))**2
      end do
      sum = MIN (sum1, sum2, sum3, sum4)
      IF ((sum .eq. sum2) .or. (sum .eq. sum4))THEN
        call revers (pppsi, zdum, nxeqd)
         WRITE(*,506)
         WRITE(nitre,506)
 506     FORMAT( '*** reversed pppsi in chekeqd  ****')
      ENDIF
      cconst = -1.0
      IF ((sum .eq. sum3) .or. (sum .eq. sum4))THEN
        call multpl1 (pppsi, nxeqd, cconst)
        WRITE(*,507)
        WRITE(nitre,507)
 507    FORMAT('*** changed sign of pppsi  in chekeqd ****')
      ENDIF
c
c --- next get the sign of ffppsi by comparing with
c --- delstar psi = (4pi/c)*R*Jphi = -4pi*R**2*pppsi-ffppsi.
c --- For convenience,
c --- we use the bicubic spline rather than the finite difference
c --- method for evaluating del-star.
c --- (the two methods should agree well enough to allow this).
c
      imslmd='1451c207'
      call my_ibcccu (p,xdum,nxeqd,ydum,nyeqd,cspln,nw,wnoperm,ier)
      call copya  (ffppsi,zdum,nxeqd)
      call revers (zdum, wdum, nxeqd)       ! reversed version of ffppsi
      sum1   = 0.0
      sum2   = 0.0
      sum3   = 0.0
      sum4   = 0.0
      imslmd = 'chekeqd'
c
c --- test inside the plasma only
c
      do 70 j=1,nyeqd
        if (ydum(j) .lt. yplmin)  go to 70
        if (ydum(j) .gt. yplmax)  go to 70
        do 71 i=1,nxeqd
          if (xdum(i) .lt. xplmin)  go to 71
          if (xdum(i) .gt. xplmax)  go to 71
          if (p (i,j) .gt. psilim)  go to 71
      imslmd='1469c207'
          call my_dbcevl1(xdum,nxeqd,ydum,nyeqd,cspln,nw,xdum(i),
     .                    ydum(j),pds,ier,icalc)
          delstar = pds(5)-pds(2)/xdum(i)+pds(6)
c
c         get pppsi,ffppsi,and ffppsi reversed on the psival grid
c
          call intrp (1,1,psival,pppsi,nxeqd,pds(1),pprim1,1)
          call intrp (1,1,psival,ffppsi,nxeqd,pds(1),ffprim1,1)
          call intrp (1,1,psival,zdum,nxeqd,pds(1),ffprim2,1)
          dumy = delstar + pprim1*fourpi*xdum(i)**2    ! in gauss
          sum1 = sum1+(dumy+ffprim1)**2
          sum2 = sum2+(dumy-ffprim1)**2
          sum3 = sum3+(dumy+ffprim2)**2
          sum4 = sum4+(dumy-ffprim2)**2
   71   continue
   70 continue
c
      sum = MIN (sum1,sum2,sum3,sum4)
      IF ((sum .eq. sum3) .or. (sum .eq. sum4))THEN
         call revers (ffppsi, zdum, nxeqd)
         WRITE(*,504)
         WRITE(nitre,504)
 504     FORMAT( '*** reversed ffppsi in chekeqd  ****')
      ENDIF
      cconst = -1.0
      IF ((sum .eq. sum2) .or. (sum .eq. sum4))THEN
        call multpl1 (ffppsi, nxeqd, cconst)
        WRITE(*,505)
        WRITE(nitre,505)
 505    FORMAT('*** changed sign of ffppsi  in chekeqd ****')
      ENDIF
c     
c --- now we can construct f(psi) from ffppsi
c
      isgnbtor = 1
      if (beqd .lt. 0.0)  isgnbtor = -1
      flimsq = (reqd*beqd)**2
      call trap2(psival,ffppsi,zdum,nxeqd)
      do j=1,nxeqd
        wdum(j) = isgnbtor * SQRT (flimsq+2.0*zdum(j))
      end do
c
c --- wdum is f(psi) constructed from the (possibly corrected) ffppsi vector
c
      call copya  (fpsi, zdum, nxeqd)
      call revers (zdum, xdum, nxeqd)    ! reversed version of fpsi
      sum1 = 0.0
      sum2 = 0.0
      sum3 = 0.0
      sum4 = 0.0
      do j=1,nxeqd
        sum1 = sum1 + (wdum(j)-fpsi(j))**2
        sum2 = sum2 + (wdum(j)+fpsi(j))**2
        sum3 = sum3 + (wdum(j)-fpsi(nxeqd-j+1))**2
        sum4 = sum4 + (wdum(j)+fpsi(nxeqd-j+1))**2
      end do
      sum = MIN (sum1,sum2,sum3,sum4)
      IF ((sum .eq. sum4) .or. (sum .eq. sum3))THEN
        call revers(fpsi,xdum,nxeqd)
        WRITE(*,501)
        WRITE(nitre,501)
  501   FORMAT( '*** reversed fpsi in chekeqd  ****')
      ENDIF

      cconst = -1.0
      IF ((sum .eq. sum2) .or. (sum .eq. sum4))THEN
        call multpl1 (fpsi,nxeqd,cconst)
        WRITE(*,502)
        WRITE(nitre,502)
 502    FORMAT('*** changed sign of fpsi  in chekeqd ****')
      ENDIF

c
c --- finally change the sign of qpsi to agree with the sign of fpsi
c
      IF (ABS (qpsi(nxeqd)) .gt. ABS (qpsi(1)))THEN
        call revers(qpsi,xdum,nxeqd)
        WRITE(*,503)
        WRITE(nitre,503)
 503    FORMAT( '*** reversed qpsi in chekeqd  ****')
      ENDIF

      do j=1,nxeqd
        qpsi(j) = isgnbtor * ABS (qpsi(j))
      end do
c
c --- the eqdsk data should now conform to the coordinate
c --- system (r,z,phi)  and the grad shafranov equation
c --- delstar psi = +u0*R*Jphi,with Jphi = -R*pppsi-ffppsi/(u0*R).
c
      return
c
c --- fatal error handling
c
 1000 if (iounit .ne. 0)
     .write  (iounit, 1010)  nxeqd, nyeqd, nw, nh
 1010 format (' subroutine CHEKEQD detects a FATAL ERROR'  /
     .        '   nxeqd,nyeqd =',2(2x,i5)                  /
     .        '   nw,nh =',2(2x,i5)                        /
     .        '   must have nxeqd .le. nw,  nyeqd .le. nh' /
     .        '   ONETWO must stop'                        )
      call STOP ('subroutine CHEKEQD: problem #1', 6)
 2000 if (iounit .ne. 0)
     .write  (iounit, 2010)
 2010 format (' subroutine CHEKEQD detects a FATAL ERROR'        /
     .        '   the eqdsk does not contain the (r,z) points'   /
     .        '   required to check the eqdsk. ONETWO must stop.')
      call STOP ('subroutine CHEKEQD: problem #2', 7)
 3000 if (iounit .ne. 0)
     .write  (iounit, 3010)
 3010 format (' subroutine CHEKEQ detects a FATAL ERROR in eqdsk:' /
     .        '   the magnetic axis setting is not correct.'       /
     .        '   ONETWO must stop.'                               )
      call STOP ('subroutine CHEKEQD: problem #3', 8)
c
      end

      subroutine chekgm

c
c ----------------------------------------------------------------------
c   subroutine to print out geometric quantities and check for
c   convergence. Convergence is assumed achieved if:
c    1) the relative change in rhoa is less than delrho
c       and
c    2) the relative change in any of the cap quantities at all grid
c       points (1...nj) is less than delcap.
c   Note that the relative change is the difference between predicted
c   and calculated quantities at the current time point. The predicted
c   value is determined by linear in time extrapolation.
c   if the above two criteria are met,ind = 0 is returned.If one or both
c   of the above criteria are not met then ind = 1 is returned,unless
c   ibypas = 1 in which case ind=0 is returned even though convergence
c   was not achieved.
c --- input
c --- INCLUDE file param:
c  kj                         dimension of vectors indexed by nj below
c
c --- INCLUDE file geom:
c  rhoa0                      rho at plasma boundary
c  rhoa                       rho at plasma boundary t = time
c  drhoadt_geom               deriv of rhoa
c
c  fcap(1..nj)
c  fcap0(1..nj)
c  dfdt(1..nj)                deriv of fcap
c
c  gcap(1..nj)
c  gcap0(1..nj)
c  dgdt(1..nj)                deriv of gcap
c
c  hcap(1..nj)
c  hcap0(1..nj)
c  dhdt(1..nj)                deriv of hcap
c
c  r2cap(1..nj)
c  r2cap0(1..nj)
c  dr2dt(1..nj)              deriv of r2cap
c
c  r2capi(1..nj)
c  r2capi0(1..nj)
c  dr2idt(1..nj)            deriv of r2capi
c
c  rcap(1..nj)
c  rcap0(1..nj)
c  drcapdt(1..nj)           deriv of rcap
c

c  rcapi(1..nj)
c  rcap0i(1..nj)
c  drcapidt(1..nj)           deriv of rcap
c
c --- INCLUDE file neo2d.i
c  eps(1..nj)
c  eps0(1..nj)
c  depsdt(1..nj)            deriv of eps
c
c  xhm2(1..nj)
c  xhm20(1..nj)
c  dxhm2dt(1..nj)           deriv of xhm2
c
c  xi11(1..nj)
c  xi110(1..nj)
c  dxi11dt(1..nj)           deriv of xi11
c
c  xi33(1..nj)
c  xi330(1..nj)
c  dxi33dt(1..nj)           deriv of xi33
c
c  xips(1..nj)
c  xips0(1..nj)
c  dxipsdt(1..nj)           deriv of xips
c
c --- INCLUDE file io
c  nout                        writes to file outone
c  ncrt                        writes to screen
c  nitre                       writes to file runlog
c  jprt                        print increment
c
c --- INCLUDE file numbrs:
c  nj                      size of transport rho grid
c
c --- INCLUDE file solcon:
c  time                     current time
c
c --- INCLUDE file soln2d:
c  delrho                   convergence criteria for rhoa
c  delcap                                            cap quantitites
c  itre                     equilibrium iteration ctr
c                             (for number of iters required
c                              to converge geometric parameters)
c  maxitr                   max iters allowed
c  ibypas
c
c --- output
c --- INCLUDE file soln2d:
c  icon                       = 0,converged,=1 not converged
c                            if ibypass = 1,icon=0 on output
c
c ----------------------------------------------------------------------
c
c
      USE param
      USE io
      USE solcon
      USE numbrs
      USE geom
      USE soln2d
      USE neo2d
      implicit  integer (i-n), real*8 (a-h, o-z)
c      include 'param.i'
c      include 'geom.i'
c      include 'io.i'
c      include 'numbrs.i'
c      include 'solcon.i'
c      include 'soln2d.i'
c      include 'neo2d.i'
c
      character*8 name(12)
      data        name / 'fcap    ', 'gcap    ', 'hcap    ', 'r2cap   ',
     .                   'r2capi  ', 'eps     ', 'xhm2    ', 'xi11    ',
     .                   'xi33    ', 'xips    ', 'RESERVED', 'UNKNOWN '/
      data        puny / 1.0e-7 /,
     .            id   / 12     / ! added 31 Mar 94 for the UNKNOWN case
c
      icon = 0
c
c ----------------------------------------------------------------------
c check rhoa for convergence
c ----------------------------------------------------------------------
c
      rhoprd = rhoa0 + drhoadt_geom * (time - timcap)
      drho   = ABS ((rhoprd-rhoa)/rhoa)
      adrhot = drhoadt_geom / rhoa
      rhoerr = ABS (rhoprd-rhoa) / (ABS (rhoa-rhoa0) + puny)
      if (drho .gt. delrho)  icon = 1
c
c ----------------------------------------------------------------------
c check cap quantities for convergence
c ----------------------------------------------------------------------
c
      dcap0 = 0.0
      dcap  = 0.0
      df    = 0.0
      dg    = 0.0
      dh    = 0.0
      deps  = 0.0
      dxhm2 = 0.0
      dxi11 = 0.0
      dxi33 = 0.0
      dxips = 0.0
c
      do j=1,nj
        gprd = gcap0(j)+dgdt(j)*(time-timcap)
        if (.not. implicit_fh) then
          fprd = fcap0(j)+dfdt(j)*(time-timcap)
          hprd = hcap0(j)+dhdt(j)*(time-timcap)
        end if
        r2prd = r2cap0(j)+dr2dt(j)*(time-timcap)
        r2iprd = r2capi0(j)+dr2idt(j)*(time-timcap)
        rprd = rcap0(j)+drcapdt(j)*(time-timcap)
        rprdi = rcap0i(j)+drcapidt(j)*(time-timcap)
        epsprd = eps0(j)+depsdt(j)*(time-timcap)
        xhm2prd = xhm20(j)+dxhm2dt(j)*(time-timcap)
        xi11prd = xi110(j)+dxi11dt(j)*(time-timcap)
        xi33prd = xi330(j)+dxi33dt(j)*(time-timcap)
        xipsprd = xips0(j)+dxipsdt(j)*(time-timcap)
        dg = ABS ((gprd-gcap(j))/gcap(j))
        if (.not. implicit_fh) then
          df = ABS ((fprd-fcap(j))/fcap(j))
          dh = ABS ((hprd-hcap(j))/hcap(j))
        end if
        dr2 = ABS ((r2prd-r2cap(j))/r2cap(j))
        dr2i = ABS ((r2iprd-r2capi(j))/r2capi(j))
        dr2r = ABS ((rprd-rcap(j))/rcap(j))
        dr2ri = ABS((rprdi-rcapi(j))/rcapi(j))
        if (j .gt. 1) then
          deps = ABS ((epsprd-eps(j))/eps(j))
          dxhm2 = ABS ((xhm2prd-xhm2(j))/xhm2(j))
          dxi11 = ABS ((xi11prd-xi11(j))/xi11(j))
          dxi33 = ABS ((xi33prd-xi33(j))/xi33(j))
          dxips = ABS ((xipsprd-xips(j))/xips(j))
        end if
        dcap = MAX (dcap,dg)
        if (.not. implicit_fh) then
          dcap = MAX (dcap,df)
          dcap = MAX (dcap,dh)
        end if
        dcap = MAX (dcap,dr2)
        dcap = MAX (dcap,dr2i)
        dcap = MAX (dcap,dr2r)
c        dcap = MAX (dcap,dr2ri)  again, rcapi is for kinetic fit only
        if (xi_include) then
          dcap = MAX (dcap,deps)
          dcap = MAX (dcap,dxhm2)
          dcap = MAX (dcap,dxi11)
          dcap = MAX (dcap,dxi33)
          dcap = MAX (dcap,dxips)
        end if
c
        jd = 0
        if (dcap .ne. dcap0) then
          jd    = j
          dcap0 = dcap
          if      (dcap .eq. df   ) then
            id =  1
          else if (dcap .eq. dg   ) then
            id =  2
          else if (dcap .eq. dh   ) then
            id =  3
          else if (dcap .eq. dr2  ) then
            id =  4
          else if (dcap .eq. dr2i ) then
            id =  5
c          else if (dcap .eq. dr2ri ) then   see rcapi comments
c            id =  5
          else if (dcap .eq. deps ) then
            id =  6
          else if (dcap .eq. dxhm2) then
            id =  7
          else if (dcap .eq. dxi11) then
            id =  8
          else if (dcap .eq. dxi33) then
            id =  9
          else if (dcap .eq. dxips) then
            id = 10
          else
            id = 12
          end if
        end if
c
      end do
c
      if (dcap .gt. delcap)  icon = 1
c
c ----------------------------------------------------------------------
c write out results,
c but results for eps,xhm2,xi11,xi33,xips not written out at present.
c ----------------------------------------------------------------------
c
      timms  = time   * 1.0e3
      timoms = timcap * 1.0e3
      call headeq (nout)
      write  (nout , 30)  timms,timoms,rhoa0,rhoprd,rhoa,adrhot,rhoerr
      write  (nitre, 30)  timms,timoms,rhoa0,rhoprd,rhoa,adrhot,rhoerr
   30 format (' geometry check ' // '  time = ', f16.2,
     .        ' msec', '   timold = ', f16.2, ' msec' // 10x,
     .        ' rho max ' / ' ',
     .        7h old   , 7h pred  , 7h new   , 8h logdt  , 8hrel err /
     .        3f7.3, 2e10.2 // 10x, 'fcap', 50x, 'gcap' /
     .        5x, 6h old  , 6h pred , 6h new  , 8h logdt  , 8hrel err )
c
      do j=1,nj,jprt
        fprd  = fcap0(j)+dfdt(j)*(time-timcap)
        gprd  = gcap0(j)+dgdt(j)*(time-timcap)
        adfdt = dfdt(j)/fcap0(j)
        adgdt = dgdt(j)/gcap0(j)
        ferr  = ABS (fprd-fcap(j))/(ABS (fcap(j)-fcap0(j))+puny)
        gerr  = ABS (gprd-gcap(j))/(ABS (gcap(j)-gcap0(j))+puny)
        write (nout , 40) j,fcap0(j),fprd,fcap(j),adfdt,ferr,
     .                    gcap0(j),gprd,gcap(j),adgdt,gerr
        write (nitre, 40) j,fcap0(j),fprd,fcap(j),adfdt,ferr,
     .                    gcap0(j),gprd,gcap(j),adgdt,gerr
   40   format (i4,3f6.3,2e10.2,10x,3f6.3,2e10.2)
      end do
c
      write  (nout , 130)
      write  (nitre, 130)
  130 format (/ 10x, 'hcap', 46x, 'r2cap' / 5x,
     .           6h old  ,6h pred ,6h new  ,8h logdt  ,8hrel err ,10x,
     .           6h old  ,6h pred ,6h new  ,8h logdt  ,8hrel err )
c
      do j=1,nj,jprt
        hprd  = hcap0(j)+dhdt(j)*(time-timcap)
        r2prd = r2cap0(j)+dr2dt(j)*(time-timcap)
        adhdt = dhdt(j)/hcap0(j)
        adrdt = dr2dt(j)/r2cap0(j)
        herr  = ABS (hprd-hcap(j))/(ABS (hcap(j)-hcap0(j))+puny)
        r2err = ABS (r2prd-r2cap(j))/(ABS (r2cap(j)-r2cap0(j))+puny)
        write (nout , 40) j,hcap0(j),hprd,hcap(j),adhdt,herr,
     .                    r2cap0(j),r2prd,r2cap(j),adrdt,r2err
        write (nitre, 40) j,hcap0(j),hprd,hcap(j),adhdt,herr,
     .                    r2cap0(j),r2prd,r2cap(j),adrdt,r2err
      end do
c
      write  (nout , 60)  drho, delrho, name(id), jd, dcap, delcap
      write  (nitre, 60)  drho, delrho, name(id), jd, dcap, delcap
      write  (ncrt , 60)  drho, delrho, name(id), jd, dcap, delcap
   60 format (/ 10x, '% change : rhoa', 11x, e12.4, '  versus', e12.4 /
     .               21x, a10, '(', i3, ')', e12.4, '  versus', e12.4)

c
c ----------------------------------------------------------------------
c determine if an iteration is to be done
c ----------------------------------------------------------------------
c
      if (icon .eq. 1)  go to 70
      write  (nout , 80)
      write  (nitre, 80)
   80 format (' no iteration is required')
      return
c
   70 if (itre .ge. maxitr .and. ibypas .eq. 1)  go to 100
      write  (nout , 90) itre
      write  (nitre, 90) itre
   90 format (' iteration number', i5, '  will be carried out')
      return
c
  100 write  (nout , 110)
      write  (nitre, 110)
  110 format (' iteration is called for, but itre > maxitr',
     .        ' and ibypas = 1, so no iteration will be done')
      icon = 0
      return
c
      end

      subroutine chkcrn (psi, nwh, psivl, kold, knew, icrnr,
     .                   kk, nh, i, i1)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c***********************************************************************
c**                                                                   **
c**     MAIN PROGRAM:  MHD FITTING CODE                               **
c**                                                                   **
c**     SUBPROGRAM DESCRIPTION:                                       **
c**          chkcrn decides which cell is next when an error in       **
c**          the cell step is likely.                                 **
c**                                                                   **
c**     RECORD OF MODIFICATION:                                       **
c**          07/09/83..........first created                          **
c**                                                                   **
c******************************************************************* HSJ
c
      dimension psi(nwh)
c
      knew = 0
      i1   = 0
c
c ----------------------------------------------------------------------
c cell #1
c ----------------------------------------------------------------------
c
      kn = kk
      if (kn .eq. kold)  go to 10
      call minmax (psi,nwh,nh,kn,psivl,iflag,knew)
      if (iflag .eq. 1)  go to 10
      i1 = i
      return
   10 go to (100,200,300,400),icrnr
c
c ----------------------------------------------------------------------
c corner #1
c cell #2
c ----------------------------------------------------------------------
c
  100 kn = kk-nh
      if (kn .eq. kold)  go to 110
      call minmax (psi,nwh,nh,kn,psivl,iflag,knew)
      if (iflag .eq. 1)  go to 110
      i1 = i-1
      return
c
c ----------------------------------------------------------------------
c cell #3
c ----------------------------------------------------------------------
c
  110 kn = kk-nh-1
      if (kold .eq. kn)  go to 120
      call minmax (psi,nwh,nh,kn,psivl,iflag,knew)
      if (iflag .eq. 1)  go to 120
      i1 = i-1
      return
c
c ----------------------------------------------------------------------
c cell #4
c ----------------------------------------------------------------------
c
  120 kn = kk-1
      call minmax (psi,nwh,nh,kn,psivl,iflag,knew)
      i1 = i-1
      return
c
c ----------------------------------------------------------------------
c corner #2
c cell #2
c ----------------------------------------------------------------------
c
  200 kn = kk-1
      if (kn .eq. kold)  go to 210
      call minmax (psi,nwh,nh,kn,psivl,iflag,knew)
      if (iflag .eq. 1)  go to 210
      i1 = i
      return
c
c ----------------------------------------------------------------------
c cell #3
c ----------------------------------------------------------------------
c
  210 kn = kk+nh-1
      if (kn .eq. kold)  go to 220
      call minmax (psi,nwh,nh,kn,psivl,iflag,knew)
      if (iflag .eq. 1)  go to 220
      i1 = i+1
      return
c
c ----------------------------------------------------------------------
c cell #4
c ----------------------------------------------------------------------
c
  220 kn = kk+nh
      call minmax (psi,nwh,nh,kn,psivl,iflag,knew)
      i1 = i+1
      return
c
c ----------------------------------------------------------------------
c corner #3
c cell #2
c ----------------------------------------------------------------------
c
  300 kn = kk+nh
      if (kn .eq. kold)  go to 310
      call minmax (psi,nwh,nh,kn,psivl,iflag,knew)
      if (iflag .eq. 1)  go to 310
      i1 = i+1
      return
c
c ----------------------------------------------------------------------
c cell #3
c ----------------------------------------------------------------------
c
  310 kn = kk+nh+1
      if (kn .eq. kold)  go to 320
      call minmax (psi,nwh,nh,kn,psivl,iflag,knew)
      if (iflag .eq. 1)  go to 320
      i1 = i+1
      return
c
c ----------------------------------------------------------------------
c cell #4
c ----------------------------------------------------------------------
c
  320 kn = kk+1
      call minmax (psi,nwh,nh,kn,psivl,iflag,knew)
      i1 = i
      return
c
c ----------------------------------------------------------------------
c corner #4
c cell #2
c ----------------------------------------------------------------------
c
  400 kn = kk+1
      if (kn .eq. kold)  go to 410
      call minmax (psi,nwh,nh,kn,psivl,iflag,knew)
      if (iflag .eq. 1)  go to 410
      i1 = i
      return
c
c ----------------------------------------------------------------------
c cell #3
c ----------------------------------------------------------------------
c
  410 kn = kk-nh+1
      if (kn .eq. kold)  go to 420
      call minmax (psi,nwh,nh,kn,psivl,iflag,knew)
      if (iflag .eq. 1)  go to 420
      i1 = i-1
      return
c
c ----------------------------------------------------------------------
c cell #4
c ----------------------------------------------------------------------
c
  420 kn = kk-nh
      call minmax (psi,nwh,nh,kn,psivl,iflag,knew)
      i1 = i-1
      return
c
      end

      subroutine chkinout (rpt, zpt, psirzp, ind, isigncur, psiedge)
c

c
c ----------------------------------------------------------------------
c set ind = 1 if (rpt,zpt) is inside plasma ...
c set ind = 0 if (rpt,zpt) is outside plasma ...
c return psirzp,value of psi at (rpt,zpt) if ind = 1
c --- INCLUDE file variables used:
c  INCLUDE file mhdpar:
c  nw,nh
c  INCLUDE file mhdgrid:
c  rmhdgrid(nw)
c  zmhdgrid(nh)
c  INCLUDE file bicube:
c  cspln                  assumed set up on input
c  pds(6)
c  INCLUDE file contour:
c  rplasmin,rplasmax
c  zplasmin,zplasmax
c ------------------------------------------------------------------ HSJ
c

      USE contour
      USE mhdpar 
      USE mhdgrid       
      USE mhdcom
      USE bicube
      USE replace_imsl,                   ONLY : my_dbcevl1

      implicit  integer (i-n), real*8 (a-h, o-z)

      include 'imsl.i'
c
      psichk = psibdry
      if (psiedge .ne. -1.0e30)  psichk = psiedge
      ind    = 1
      imslmd = 'chkinout'
      if (rpt .lt. rplasmin)  ind = 0
      if (rpt .gt. rplasmax)  ind = 0
      if (zpt .lt. zplasmin)  ind = 0
      if (zpt .gt. zplasmax)  ind = 0
      if (ind .eq. 0       )  return
      call my_dbcevl1 (rmhdgrid,nw,zmhdgrid,nh,cspln,nw,rpt,zpt,
     .                 pds,ier,1)
      psirzp = pds(1)
      if (isigncur .gt. 0 .and. psirzp .lt. psichk)  ind = 0
      if (isigncur .lt. 0 .and. psirzp .gt. psichk)  ind = 0
      return
c
      end

      subroutine chord (nlimiter, xlimiter, ylimiter, xlimcen,
     .                  ychord, x1, x2)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c inputs
c ----------------------------------------------------------------------
c nlimiter: number of limiter points
c xlimiter: x coordinate of limiter point (meters)
c ylimiter: y coordinate of limiter point (meters)
c xlimcen : x coordinate of center of the limiter
c ychord  : y location of chord passing through
c           the limiter (parallel to the xaxis) (meters)
c ----------------------------------------------------------------------
c outputs
c ----------------------------------------------------------------------
c x1,x2   : define left and right end points of chord (meters)
c
      dimension  xlimiter(*), ylimiter(*)
c
      do 210 i=1,nlimiter-1
        del1 = ychord-ylimiter(i)
        del2 = ychord-ylimiter(i+1)
        if (    del1*del2 .gt. 0.0        )  go to 210
        if (xlimiter(i+1) .ne. xlimiter(i))  go to 200
        xint = xlimiter(i)
        yint = ychord
        go to 204
  200   s1   = (ylimiter(i+1)-ylimiter(i))/(xlimiter(i+1)-xlimiter(i))
        b1   = ylimiter(i)-s1*xlimiter(i)
        s2   = 0.0
        b2   = ychord
        xint = (b1-b2)/(s2-s1)
        yint = s1*xint+b1
  204   if (xlimiter(i) .lt. xlimcen) x1 = xint
        if (xlimiter(i) .gt. xlimcen) x2 = xint
  210 continue
c
      return
c
      end

      subroutine cmpltcyr (diagl, diag, diagu, f, nw, nw1, nh, nhpwr,
     .                     phi, phi1, diag1, v, wk, wk2)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c --- subroutine solves the set of three point vector equations
c ---         -y(i,j-1)+c(i,j)*y(i,j)-y(i,j+1) = f(i,j)
c --- for each j,from j = 2 to j=nh-1,by complete cyclic reduction.
c --- all intermediate solutions are generated using
c --- a standard tridiagonal solver. thus the upper limit of index i,
c --- nw, is not restricted in any way. however nh-1 must be a power of 2.
c --- here c is an (nw-2) by (nw-2) tridiagonal matrix with
c --- lower diagonal diagl(i),diagonal diag(i) and upper diagonal
c --- diagu(i),i = 1,2,...nw-2.
c --- input:
c
c  diagl(i)        i = 1,..nw-2 is the lower diagonal of c
c  diag(i)         i = 1,..nw-2 is the        diagonal
c  diagu(i)        i = 1,..nw-2 is the upper diagonal of c
c                  diagl,diag,and diagu are not changed by this routine.
c                  diagl(1),and diagu(nw-2) are not used.
c  f(i,j)          i = 2,..nw-1,j = 2..nh-1 is the inhomogenous term
c                  on input f(i,1),f(i,nh),i = 1,2,...nw must contain
c                  the respective boundary conditions
c  nw              is the number of grid points in the i direction
c  nw1             is the leading dimension of the array f as
c                  declared in the calling program. (nw  .le. nw1)
c  nh              is the number of grid points in the j direction
c                  nh-1 must be a power of 2 (see nhpwr).
c  nhpwr            is defined so that 2**nhpwr = nh-1 (defines the j direction)
c  output---
c  f(i,j)          the solution (found for i = 2,...nw-1 and j=2,..nh-1)
c --- storage required:
c  phi(i)            i = 1,2,...nw-2
c  phi1(i)           i = 1,2....nw-2
c  diag1(i)          i = 1,2....nw-2
c  v(i)              i = 1,2....nw-2
c  wk(i)             i = 1,2....nw-2
c  wk2(i)            i = 1,2....nw-2
c
c  notes:
c              this is a general 3 point vector equation solver.
c              to use it you must specify the diagonals
c              of matrix c and the inhomogeneous term f. these
c              vector equations arise naturally due to finite
c              differencing elliptic equations.
c               to solve the intermediate sets of tridiagonal
c               scalar equations we have a number of choices:
c               a)tridiagonal elimination (eg. subroutine TRIDIAG)
c                 this method is used here.
c               b)tridiagonal elimination with pivoting
c                 (subroutine SGTSL could be used but it is slower
c                  and offers no increase in accuracy if matrix c
c                  is diagonally dominant)
c               c)cyclic reduction of 3 point scalar equations
c                  Takes more operations to complete (~11n
c                  versus ~5n). also this method would restrict nw-1
c                  to be a power of 2.
c              the tridiagonal solver used here assumes that
c              the matrix c* is diagonally dominant. if c* is not
c              diagonally dominant the solution can fail or be in
c              error by a large amount. (for such cases the
c              subroutine tridiag should be replaced with a
c              version that incorporates pivoting).
c              here c* is the matrix c with a modified diagonal
c              term. the modification adds  -2.0*COS (some angle)
c              to the diagonal.
c              in terms of matrix c itself it means that
c              (cy,y)should be .ge. 2(y,y) which makes c*
c              positive definite. if this does not hold then
c              the cyclic reduction method is not valid!
c              the reason for this is the method used to reduce powers of
c              matrix c during the reduction process.
c
c              two trigonometric expressions are required
c              in this subroutine: cos ((2*l-1)*pi/2**k )
c              and SIN ( (2*l-1)*pi/2**k).these are calculated
c              recursively using half angle,double angle
c              and sum of two angle
c              trigonometric identities,for k = 1,2,...nhpwr and
c              l = 1,2,....2**(k-1)
c
c           conversion to CRAY requires only that double precision
c           parameters (used for recursive definition of angles)
c           be made single:
c             a) remove real*8 declaration line
c             b) change all occurrences of 1.0e0,2.0e0 to 1.0e0,2.0e0
c             c)                           SQRT          SQRT
c             d) exclamation-delimited comments are ok (CIVIC compiler)
c      (above conversions have been done here, so this is CRAY version)
c
c    This subroutine was created from first priciples by HSJ, 11/15/89
c    It is designed to allow arbitrary grids in one direction and
c    for easy inclusion of the facr(l) method later.
c ------------------------------------------------------------------ HSJ
c
      dimension  f(nw1,nh),diagl(*),diag(*),diagu(*)
      dimension  phi(*),phi1(*),v(*),diag1(*),wk(*),wk2(*)
c
      root2inv = 1.0e0 / SQRT (2.0e0)
c
c --- the forward reduction steps
c
      cosdj = -1.0e0
      cosdt = 1.0e0
      nred  = nhpwr - 1        ! nred steps are required
      jpow  = 1
      do 10 k=1,nred           ! k is the reduction step index
        jpowm1 = jpow          ! 2**(k-1)
        jpow   = jpow*2        ! 2**k
        jstart = jpow + 1      ! 2**k+1
        jend   = nh - jpow     ! nh-2**k
        jstep  = jpow          ! 2**k
c
c --- set up for trigonometric recursions
c
        ddum1 = cosdt
        cosdt = root2inv * SQRT (1.0e0 + ddum1)    ! half-angle formulae
        sindt = root2inv * SQRT (1.0e0 - ddum1)
        if (k .eq. 1)  cosdt = -cosdt
        ddum1 =  root2inv * SQRT (1.0e0+cosdj)
        ddum2 = -root2inv * SQRT (1.0e0-cosdj)
        cosdj = ddum1
        do 20 j=jstart,jend,jstep
          jdm = j-jpowm1
          jdp = j+jpowm1
          do i=2,nw-1       ! define vector phi
            phi(i-1) = f(i,jdm) + f(i,jdp)
          end do
          cosdi = ddum1     ! reset cosdi,sindi each time l loop is done
          sindi = ddum2
          m1    = -1
          do 40 l=1,jpowm1
          m1    = -m1
c
c --- evaluate COS ((2*l-1)*pi/2**k) = cosdi and SIN ((2*l-1)*pi/2**k)
c --- =sindi
c
          cosdii = cosdi*cosdt-sindi*sindt
          sindi  = sindi*cosdt+cosdi*sindt
          cosdi  = cosdii
          cosdii = -2.0e0*cosdi
          alpha  = m1*sindi/jpowm1
            do i=1,nw-2         ! adjust diagonal of c and vector phi
              diag1(i) = diag(i) + cosdii
              phi1 (i) = phi (i) * alpha
            end do
c
c --- solve the tridiagonal system with modified diagonal terms
c
            call tridiag (diagl, diag1, diagu, phi1, v, wk, nw-2)
            do i=2,nw-1
              f(i,j) = f(i,j) + v(i-1)
            end do
   40     continue                    ! end loop on l
          do i=2,nw-1
            f(i,j) = 0.5 * f(i,j)
          end do
   20   continue                      ! end loop on j
   10 continue                        ! end loop on k
c
c --- reduction to a single set of equations is complete, now backsolve
c
c --- set up the initial cosine and sine terms
c
      ddum1 =  cosdt
      cosdt =  root2inv * SQRT (1.0e0+ddum1)
      sindt =  root2inv * SQRT (1.0e0-ddum1)
      ddum1 =  cosdt
      cosdt =  root2inv * SQRT (1.0e0+ddum1)
      sindt =  root2inv * SQRT (1.0e0-ddum1)
      ddum1 =  cosdt
      cosdj =  root2inv * SQRT (1.0e0+ddum1)
      sindj = -root2inv * SQRT (1.0e0-ddum1)
c
c --- cosdt is now COS (2*pi/2**(nhpwr+1))
c --- sindt is now SIN (2*pi/2**(nhpwr+1))
c --- cosdj is now COS (-pi/2**(nhpwr+1))
c --- sindj is now SIN (-pi/2**(nhpwr+1))
c
      nred = nred+1
      jpow = 2**nred          ! 2**nhpwr
      do 100 k=nred,1,-1
        jstep  = jpow
        jpow   = jpow/2       ! 2**k-1
        jstart = jpow+1
        jend   = nh-jpow
c
c --- set up for trigonometric recursions (angle dt is doubled)
c
        sindt = 2.0e0*sindt*cosdt                ! double-angle formulae
        cosdt = 2.0e0*cosdt*cosdt - 1.0e0
        if (k .eq. 1)  cosdt = -1.0e0
        sindj = 2.0e0*sindj*cosdj
        cosdj = 2.0e0*cosdj*cosdj-1.0e0
        do 110 j=jstart,jend,jstep
          jdm = j-jpow
          jdp = j+jpow
          do i=2,nw-1
            phi(i-1) = f(i,jdm)+f(i,jdp)
            wk2(i-1) = f(i,j)
          end do
          cosdi = cosdj   ! reset cosdi, sindi each time l loop is done
          sindi = sindj
          m1    = -1
          do 130 l=1,jpow
            m1 = -m1
c
c --- evaluate COS ((2*l-1)*pi/2**k) = cosdi and SIN ((2*l-1)*pi/2**k)
c --- =sindi
c
          cosdii = cosdi*cosdt-sindi*sindt
          sindi  = sindi*cosdt+cosdi*sindt
          cosdi  = cosdii
          cosdii = -2.0e0*cosdi
          alpha  = m1*sindi/jpow
            do i=1,nw-2
              diag1(i) = diag(i) + cosdii
              phi1 (i) = wk2 (i) + alpha * phi(i)
            end do
c
c --- solve the tridiagonal system with modified diagonal terms
c
            call tridiag (diagl, diag1, diagu, phi1, v, wk, nw-2)
c
c --- progressively build up the required solution vectors
c
            if (l .eq. 1) then
              do i=2,nw-1
                 f(i,j) = v(i-1)
              end do
            else
              do i=2,nw-1
                f(i,j) = f(i,j)+v(i-1)
              end do
            end if
  130     continue    ! l loop
  110   continue      ! j loop
  100 continue        ! k loop
      return
c
      end

      subroutine cntour1 (xaxd,yaxd,psivl,xemin,xemax,yemin,yemax,yxmin,
     .                    yxmax,xymin,xymax,dang,arcl,bperr,dx,dy,xmin,
     .                    xmax,ymin,ymax,iauto,iautoc,xc,yc,ipts,x,nw,y,
     .                    nh,cspln,n2cspln,nh2,itty,iptsm,ierr,
     .                    bpmag,iconvg,delta_psi)
c

c
c ----------------------------------------------------------------------
c
c                     version 5/17/94
c                     version 7/07/94
c             (same as 5/17/94 but uses IMSL spline routines)
c
c          CNTOUR1 AND ITS OLDER VERSION, CNTOUR2, ARE NOT INTENDED TO
c          TO BE GENERAL-PURPOSE CONTOURING ROUTINES. THESE ROUTINES
c          ARE DESIGNED SPECIFICALLY FOR FINDING PLASMA FLUX
c          SURFACES AND ASSUME THAT:
c
c          1) THE PLASMA SURFACE IS CONVEX (BEAN SHAPES MAY BE
c             PROBLEMATIC)
c          2) THE FUNCTION TO BE CONTOURED IS MONOTONICALLY
c             DECREASING IN THE REGION OF INTEREST
c
c     this routine uses the bicubic spline representation of psi to get
c     psi values at an arbitrary point on the MHD grid. the spline
c     coefficients must exist (in cspln) before this subroutine is called.
c     given (xaxd,yaxd) and a psi value, psivl, generate a contour of
c     ordered points,(xc(i),yc(i)), i = 1,ipts, which has (xaxd,yaxd) as
c     as an interior point. the contour must fully encircle (xaxd,yaxd).
c     the search is limited to a rectangle specified by (xmin,xmax,ymin
c     ymax).  dx and dy determine the basic increment in x and y for
c     the coarse grid search.  they should be picked so that over the
c     distance dx or dy, psi is single valued.  a sliding interval
c     search rather than a binary search is used for the coarse grid so
c     that non-monotonic psi can be handled.
c     arcl is the arclength in meters taken from the current point to get
c     to the next point.
c     bperr is relative change in
c     poloidal b field between (xc(i),yc(i)) and (xc(i+1),yc(i+1)).
c     if the number of elements (xc,yc) exceeds the limit,iptsm,
c     bperr is relaxed (up to twice) by a constant increment
c     dbperr before an error exit is taken.
c     bperr = 0.03 is suggested.
c     (xemin,xemax,yemin,
c     yemax) are min and max values of contour.  (yxmin,yxmax), are y
c     values at x = xemin and x = xemax respectively.  (xymin,xymax) are x
c     values at y = yemin and y = yemax respectively.  iauto is used as an
c     "error recovery switch".  that is, if a flux surface is found to
c     pass outside the search box an error exit is taken by way of label
c     1000, if iauto = 0.  if iauto = 1, the input value of psi, psivl, is
c     increased toward the value of psi at (xaxd,yaxd) until a closed
c     surface inside the search box is found.  obviously this option
c     makes sense only for the limiter flux surface.  it could be used
c     to find the plasma boundary if an appopriate search box were
c     defined.  however its primary use here is to slightly modify
c     psilim so that a closed limiter flux surface is found.  the option
c     is necessary due to the fact some codes generate eqdsks
c     using linear interpolation, which differs from the bicubic inter-
c     polation used here.  the user is informed that psivl was changed
c     by returning iautoc = 1.  if no change occurred iautoc = 0.
c     if no error ocurred ierr = 0 is returned. if an error ocurred then
c     ierr = 1 is returned and the results of this calculation are not valid.
c
c --- input (see description above)
c
c  xaxd
c  yaxd
c  psivl   psi value for which contour is desired
c  arcl    arc length of step (in same units of length as xaxd for example)
c  bperr   arcl,bperr control point density on contour returned
c  dx
c  dy      defines step length taken
c          input dx and dy as 0.0 to set step length in this subroutine
c  xmin
c  xmax
c  ymin
c  ymax    define the search box
c             note: the search box could be defined as the entire
c                   MHD grid (i.e., xmin = x(1),xmax = x(nw),ymin = y(1),
c                   ymax = y(nh)) however this will normally not work
c                   satisfactorily because there are local minima and
c                   maxima around the f coils which would confuse the
c                   contour tracing and require many more calculations.
c                   (these additional calculations are not done here!)
c                   the best way to define the search box is to set it
c                   equal to a rectangle,within which the plasma must
c                   reside under all circumstances. the only such
c                   rectangle to suggest itself is the one defined by
c                   the extremes of the limiter. hence it is suggested
c                   that xmin be set to the min x value of the limiter
c                   points vector,xmax the max x value of limiter points
c                   vector and similarly for ymin,ymax.
c  iaouto   = 1      automatic error recovery, = 0 no recovery
c                  (iauto may be set to 0 in this subroutine
c                  if a failure occurs)
c                   There are three ways to run
c                     a)iauto=0    should be ok for all interior plasma
c                                  surfaces. The value of psi to be traced,
c                                  psivl,is not altered. Either the routine
c                                  returns with the desired contour or
c                                  a failure is indicated with ierr
c                                  (iconvg is not used in this mode)
c
c                        the other two modes are intended for use on the
c                        plasma surface:
c                     b)iauto=1,iconvg=0 . May change the value of psivl!!
c                       (hence use call by reference). If the contour
c                       with the value of psivl could not be traced or
c                       passes outside the "box" (see below) the value
c                       of psivl will be REPEATEDLY  brought closer
c                       to the value
c                       on the magnetic axis. Either by some fraction
c                       of delta_psi if it is set or some fraction of
c                       psiaxis-psilim if it is not. The first contour
c                       succesfully traced will be returned.
c                     c)iauto=1,iconvg=1 as b except that a binary search
c                       is done to find  the plasma surface. FOR THIS
c                       PURPOSE THE "PLASMA SURFACE" IS DEFINED AS THE
c                       LARGEST CLOSED FLUX SURFACE THAT CAN BE FOUND
c                       WHICH DOES NOT PASS OUTSIDE THE SEARCH BOX.
c                       The actual search box used is not critical for
c                       diverted plasmas but is very important for
c                       limited plasmas. It is recommended that the
c                       boundary finder "bound" be used in those cases
c                       where an appropriate psivl for the plasma
c                       boundary is not known a priori.
c
c  x(1..nw)
c  y(1..nh)        the MHD grid vectors
c  cspln(nw,nh,4)  bicubic spline coefficient array
c  cspln(n2cspln,nw,nh2)  for IMSL routines
c  itty            iounit for diagnostic messages. set to 0 to suppress.
c  iptsm           max allowed number of points to be returned in
c                  vectors xc,yc
c  iconvg          used only if iauto = 1 is set . if iauto = 1 and iconvg = 1
c                  the contour is converged to the outermost cntour
c                  (meaning the contour with the smallest psi) that never
c                  leaves the search box given by xmin,xmax,ymin,ymax.
c delta_psi        this is the psi grid spacing. Define delta_psi as
c                  psi(penultimate value)-psi(last(i.e., edge) value)
c                  If the value of psi to be traced leads to an error
c                  and iauto=1 then delta_psi is used to control the
c                  amount by which the value of psi is adjusted. Otherwise
c                  delta_psi is not used,so you set its value accordingly.
c
c --- output:
c
c  xemin
c  xemax
c  yemin
c  yemax
c  yxmin
c  yxmax
c  xymax
c  xymin             limits of contour (see above)
c  xc(1,..ipts)
c  yc(1,..ipts)      the ipts contour points found
c  bpmag(1,...,ipts) magnitude of poloidal bfield at xc, yc
c  ipts
c  ierr     error flag, = 0 no error, = 1 error, output not usable.
c           if ierr > 1 then ierr equals the number of points found
c           on the contour up to this point. since this number will
c           be exceeded if the code continues,the subroutine must
c           stop. the error condition on bperr has allready been relaxed
c           before this condition causes the subroutine to terminate.
c           hence in order to continue this subroutine could be called
c           again with a larger value of arcl.
c           (note: in practice we have not found this to be a problem,
c           provided that iptsm is large enough to begin with to give
c           a reasonable representation of the contour).
c
c ------------------------------------------------------------------ HSJ
c
      USE replace_imsl,                   ONLY : my_dbcevl1

      implicit  integer (i-n), real*8 (a-h, o-z)

      dimension  pds(6), xc(*), yc(*), x(nw), y(nh),
     .           cspln(n2cspln,nw,nh2), bpmag(*)
      data       fpiov4     , spiov4     , twopi      , nitr8
     .         / 3.926990818, 5.497787145, 6.283185307, 10 /
      logical    bptol, reverse
c
      nstepc = 2  ! # steps along ray before tangent line method is used
      bpmintol = 0.05 ! if min bp drops below bpmintol and iauto = 1..
c                     ..the contour is searched for an x point
      itry     = 0
      irestart = 0
      psiin    = 0.0
      ierr     = 0
      iautoc   = 0
      bperrsav = bperr
      dbperr   = bperr/2.0
      arclmin  = 0.10 * arcl ! used to set minimum angle increment below
      ilsave   = 1
      jlsave   = 1
c
c --- get psi at (xaxd,yaxd)
c
****  call eval_bicubic_spline (x,nw,y,nh,cspln,nw,nh,xaxd,yaxd,pds,ier,
****.                           1,ilsave,jlsave)
c

      call my_dbcevl1 (x, nw, y, nh, cspln, nw, xaxd, yaxd, pds, ier, 1)
      if (ier .ne. 0) then
        if (itty .ne. 0)  write (itty, 4) ier
    4   format (/ ' subroutine CNTOUR1,cray207, detects an error:'   /
     .            ' EVAL_BICUBIC_SPLINE returned ier =', i5 /
     .            ' while getting psi on axis'               )
        call STOP ('subroutine CNTOUR1: problem #1', 157)
      end if
      psiaxd = pds(1)
c
 
      if (psiaxd .lt. psivl) then
        if (itty .ne. 0)  write (itty, 3)  psiaxd, psivl,xaxd,yaxd
    3   format (/ ' subroutine CNTOUR1,cray207, detects an ERROR:'/
     .            ' psi on magnetic axis    ', 1pe16.8          /
     .            ' psi contour to be found ', 1pe16.8          /
     .            ' xaxd,yaxd ',2(1pe16.8,x)
     .            ' must have max psi on axis, non-recoverable' /
     .            ' execution must stop'                         )
        call STOP ('subroutine CNTOUR1: problem #2', 158)
      end if
c
    5 xemin    =  1.0e10     ! entry point for generating contour..
      xemax    = -1.0e10     ! ..come back to here if a correctable..
      yemin    =  1.0e10     ! ..failure occured.
      yemax    = -1.0e10
      nothing  = -1
      xymin    =  0.0
      xymax    =  0.0
      yxmin    =  0.0
      yxmax    =  0.0
      ipts     =  0
      irestart = irestart + 1
      tstart   = fpiov4
      tend     = tstart + twopi
      thet     = tstart
      dthetmin = 1.0e-5   ! min. increment in theta from one ray to next
c                           note that this is increased below,
c                           after the first point is found
      dthet = 0.0
      dxx   = dx
      dyy   = dy
      if (dxx .eq. 0.0)  dxx = x(3)-x(1)  ! 2.0 times grid spacing is OK
      if (dyy .eq. 0.0)  dyy = y(3)-y(1)
      dxmin = 0.05*(x(2)-x(1)) ! subroutine stops if dxx becomes < dxmin
      dymin = 0.05*(y(2)-y(1)) ! subroutine stops if dyy becomes < dymin
      dx    = dxx
      dy    = dyy
      if (dxmin .ge. dxx)  dxmin = 0.5 * dxx
      serrt = 3.5e-06
      derrt = 0.5e-07
c
c --- serrt is absolute error convergence criteria for Newton's method below
c
c --- loop over theta from tstart to tend ( = tstart+twopi)
c
   10 thet   = thet + dthet
      if (thet .gt. twopi .and. tstart .ne. 0.0) then
        thet = thet - twopi
        tend = tend - twopi
      end if
      thet   = MIN (thet, tend)
      if (thet .eq. tend)  go to 200   ! normal exit from routine
      dxx    = dx
      dyy    = dy
      iqpt   = 0
c
c --- get equation of ray emanating from (xaxd,yaxd)
c
      call findeqlin (xaxd, yaxd, thet, a, bincp, iflg, isgn)
c
c --- now have y = a*x+bincp    (iflg = 0)   or
c ---          x = a*y+bincp    (iflg = 1)
c
c     start search from axis point
c
      x1   = xaxd
      y1   = yaxd
      bp1  = 0.0D0
      psi1 = psiaxd
      cost = COS (thet)
      sint = SIN (thet)
c
c --- sliding interval search. max width of interval ~1.41*(dx or dy)
c
      nsteps = 0
      cmult  = 1.0
      xsave  = -1.0e20
   40 if (iflg .eq. 1)  go to 50
c
c --- search in x
c
      xmult = 1.0
   41 x2    = x1+isgn*dxx*xmult*cmult    ! increment x
      if (x2 .eq. x1)  go to 1000        ! psivl is not within reach
      if (ABS(x2- x1)  .lt. 1.e-12 ) go to 1000
      x2    = MAX (x2,xmin)
      x2    = MIN (x2,xmax)              ! limit the search
      if (x2 .eq. x1)  go to 1000        ! psivl is not within reach
      y2    = a*x2+bincp                 ! get corresponding value of y
      if (y2 .lt. ymin .or. y2 .gt. ymax) then
        xmult = 0.5 * xmult
        go to 41
      end if
      go to 60
c
c --- search in y
c
   50 ymult = 1.0
   51 y2    = y1+isgn*dyy*ymult*cmult    ! increment 
      if (y2 .eq. y1)  go to 1000        ! psivl is not within reach
      if (ABS(y2- y1)  .lt. 1.e-12)  go to 1000  ! psivl is not within reach
      y2    = MAX (y2,ymin)
      y2    = MIN (y2,ymax)              ! limit the search
      if (y2 .eq. y1)  go to 1000        ! psivl is not within reach
      x2    = a*y2+bincp                 ! get corresponding value of x
      if (x2 .lt. xmin .or. x2 .gt. xmax) then
          ymult = ymult * 0.5
          go to 51
      end if
c
** 60 call eval_bicubic_spline (x,nw,y,nh,cspln,nw,nh,x2,y2,pds,ier,1,
**** .                          ilsave,jlsave)
   60 call my_dbcevl1 (x,nw,y,nh,cspln,nw,x2,y2,pds,ier,1)
      if (ier .ne. 0) then
        write (6, *)  'CNTOUR1 reports, call to EVAL_BI.., ier =', ier
        write (6, *)  'arguments    x2, y2 =', x2  , y2
        write (6, *)  'xgrid(1), xgrid(nw) =', x(1), x(nw)
        write (6, *)  'ygrid(1), ygrid(nh) =', y(1), y(nh)
        call STOP ('subroutine CNTOUR1: problem #3', 159)
      end if
      psi2 =  pds(1)
      dpsi = (psivl-psi1)*(psivl-psi2)
      if (dpsi .le. 0.0)  go to 70
c
c ----------------------------------------------------------------------
c
      if (cmult .lt. 0.99) then
        cmult = cmult * 0.5
        go to 40
      end if
c
c ----------------------------------------------------------------------
c --- if psi starts to increase we are in the vicinity of the saddle
c --- (i.e., x point). in this case we want to be sure that we stay
c --- on the contour which envelopes the plasma. the only robust way
c --- to do this is to first search for the minimum in psi along
c --- the ray. onece this minimum is found we guarantee that the
c --- point we find is on the plasma contour by not allowing the ray
c --- to extend past this minimum.
c --- to find the minimum we use Newton's method to search for the
c --- point along the ray where the directional derivative of pis
c --- along the ray is zero.
c ----------------------------------------------------------------------
c
      reverse = .false.
      if (psi2-psi1 .gt. 0.0) then    ! psi is increasing
         reverse  = .true.
         xsave    = x1   ! save current guess for possible later restore
         ysave    = y1
         psisave  = psi1
         isgnsave = isgn
         xn       = x1
         yn       = y1
         costsq   = cost * cost
         sintsq   = sint * sint
         newti    = 0
         step     = serrt + 2
         do while (ABS (step) .gt. serrt .and. newti .le. nitr8)
****         call eval_bicubic_spline (x,nw,y,nh,cspln,nw,nh,xn,yn,pds,
**** .                                 ier,6,ilsave,jlsave)
             call my_dbcevl1 (x,nw,y,nh,cspln,nw,xn,yn,pds,ier,6)
             step  = -(pds(2)*cost+pds(3)*sint)/(costsq*pds(5)
     .                  + 2.0*cost*sint*pds(4) + sintsq*pds(6))
             xn    = xn+step*cost
             yn    = yn+step*sint
             newti = newti+1
         end do
         if (newti .le. nitr8) then          ! Newton's method converged
             isgn = -isgn
             if (iflg .eq. 1)  dyy = 0.5 * dyy
             if (iflg .eq. 0)  dxx = 0.5 * dxx
             x2   = xn
             y2   = yn
             psi2 = pds(1)
c
c            we have found the minimum in psi along the ray
c            this minimum must be less than or equal to the value for
c            which we are searching,psivl. If this is not the case then
c            we can't find psivl on the ray so abandon the search:
c
             if (psi2 .gt. psivl)  go to 1000
         else
           go to 1500            ! Newton's method failed, try to refine
         end if
       end if
       x1     = x2
       y1     = y2
       psi1   = psi2
       nsteps = nsteps+1
       go to 40
c
c ----------------------------------------------------------------------
c
c --- now have psivl between psi1 and psi2, converge using Newton-Raphson
c
   70 newti = 0
      if (iflg .eq. 1)  go to 75
      xn = x1+isgn*dxx * 0.5
      yn = a*xn+bincp
      go to 80
   75 yn = y1+isgn*dyy * 0.5
      xn = a*yn+bincp
** 80 call eval_bicubic_spline (x,nw,y,nh,cspln,nw,nh,xn,yn,pds,ier,3,
**** .                          ilsave,jlsave)
   80 call my_dbcevl1(x,nw,y,nh,cspln,nw,xn,yn,pds,ier,3)
      if (ier .ne. 0) then
        newti = nitr8
        go to 76
      end if
      dpsids =  pds(2)*cost+pds(3)*sint
      dpsi   =  pds(1)-psivl
      serr   = -dpsi/dpsids
      if (ABS (serr) .lt. serrt)  go to 90
      if (psivl .eq. 0.0 .and. ABS (dpsi) .lt. derrt)  go to 90
      if (psivl .ne. 0.0) then
        if (ABS (dpsi/psivl) .lt. derrt)  go to 90
      end if
      delx  = serr*cost
      dely  = serr*sint
      xn    = xn+delx
      yn    = yn+dely
      newti = newti+1
   76 if (newti .ge. nitr8) then
****       if ( dxx .le. dxmin)  stop 'newti'
****       if ( dyy .le. dymin)  stop 'newti'
           if ( dxx .le. dxmin)  go to 90
           if ( dyy .le. dymin)  go to 90
           if (iflg .eq. 0    )  dxx = 0.5 * dxx
           if (iflg .eq. 1    )  dyy = 0.5 * dyy
           if (iqpt .eq. 0    )  go to 40           ! 6/11/93
           thet = thet - dthet
           go to 10
      end if
      go to 80
c
c ----------------------------------------------------------------------
c --- end of Newton iteration
c --- check for sufficient accuracy in point spacing as determined by thet
c --- accuracy test is based on a relative error in poloidal b field of bperr
c ----------------------------------------------------------------------
c
   90 bp2 = SQRT (pds(2)**2 + pds(3)**2) / xn
      if (thet .eq. tstart)  go to 100          ! first point of contour
c
c --- if relative change in bp is less than bperr store the point just found:
c
      if (ABS (bp2-bp1) / MAX (bp2,bp1) .lt. bperr)  go to 100
c
c --- if dthet is currently at its smallest allowed value of dthetmin
c --- and relative change in bp is still greater than bperr continue anyway:
c
      if (dthet .eq. dthetmin)  go to 100
c
c --- spacing too large for grad psi. decrease theta and try again
c
****  thet = thet-dthet             ! reset to previous value
c
      dthet = dthet * 0.5           ! take half the angle increment
      dthet = MAX (dthet,dthetmin)  ! dthet lower limit must be observed
c
c --- we found an acceptable point. collect it and set up for next point
c
  100 ipts = ipts+1
      if (ipts .le. iptsm-1)  go to 150
c
c --- we ran out of storage for the contour points. assuming iptsm is
c --- set to a reasonable value this should not happen unless bperr is
c --- set unrealistically small. we increase bperr by dbperr and try
c --- to generate the contour again with the relaxed error requirement.
c
      if (bperr-bperrsav .ge. 2.0*dbperr)  go to 110
      bperr = bperr+dbperr
      go to 5
  110 if (itty .ne. 0) then
      write  (itty, 1040)  iptsm,dtarcl,bperr,psivl
 1040 format
     . ('  more points on plasma '                                     /
     .  '  boundary were found than dimension of xcontr,ycontr allows' /
     .  '  allowed number = ',i5                                       /
     .  '  adjustment in point spacing will be made',
     .  '  dtarcl,bperr,psivl =',5(2x,1pe12.4))
      end if
      ierr = ipts
      return
  150 xc(ipts)    = xn  ! found  new point; save and set up for next ray
      yc(ipts)    = yn
      bpmag(ipts) = bp2
      xemin       = MIN (xemin,xn)
      xemax       = MAX (xemax,xn)
      if (xemax .eq. xn)  yxmax = yn
      if (xemin .eq. xn)  yxmin = yn
      yemin       = MIN (yemin, yn)
      yemax       = MAX (yemax, yn)
      if (yemax .eq. yn)  xymax = xn
      if (yemin .eq. yn)  xymin = xn
c
      if ( ipts .eq. 1 ) then  ! get a radius estimate so we can set..
c                              ..the theta increment based on arc length
        estmrad  = SQRT ((xn-xaxd)**2+(yn-yaxd)**2)
        dthetmin = arclmin/estmrad
        dtarcl   = arcl/estmrad
      end if
c
      if (ABS (bp2-bp1) / MAX (bp2,bp1) .le. 0.25 * bperr)
     .dthet = 2.0 * dthet
      bp1   = bp2                    ! needed for error test on next ray
      dthet = MIN (dthet, dtarcl  )  ! don't let dthet get too large
      dthet = MAX (dthet, dthetmin)  ! don't let dthet get too small
c
c ----------------------------------------------------------------------
c --- for contours sufficiently far removed from xaxd,yaxd,that the
c --- search along the ray is expensive we use the following alternative.
c --- an approximation to the new point,(xns,yns)is found by moving along
c --- the tangent line, for a distance arcl.  see subroutine NEW_POINT.
c --- if we are near an x point, as signaled by bpmag < bpmintol,
c --- then we switch back to the ray method, with small
c --- dx,dy increments to avoid crossing the x point.
c ----------------------------------------------------------------------
c
      bptol = bpmag(ipts) .gt. bpmintol
      if (nsteps .gt. nstepc .and. bptol) then
c
c          use tangent line method:
c
           call new_point (pds,arcl,xn,yn,thet,twopi,spiov4,xmin,xmax,
     .                     ymin,ymax,yaxd,xaxd,ier,sint,cost,
     .                     xns,yns,thetnew)
           if (ier .gt. 0)  go to 10    ! failed, go back use ray method
           dthet = thetnew -  thet
           thet  = thet    + dthet      ! do not use thetnew here
           if (thet .gt. twopi .and. tstart .ne. 0.0) then
             thet = thet - twopi
             tend = tend - twopi
           end if
           thet   = MIN (thet, tend)
           if (thet .eq. tend) then
             nothing = 1
             go to 200
           end if
           sint  = SIN (thet)
           cost  = COS (thet)
           xn    = xns
           yn    = yns
           newti = 0
           iqpt  = 1
           go to 80                   ! skip directly to Newton's method
      else                            ! use ray method
           go to 10
      end if
c
c --- done with this contour. close contour
c
  200 ipts        = ipts+1
      xc(ipts)    = xc(1)
      yc(ipts)    = yc(1)
      bpmag(ipts) = bpmag(1)
c
      bperr = bperrsav
      if (itry .ne. 0 .and. iconvg .eq. 1) then
        psiin  = psivl
        psivln = 0.5 * (psiout+psiin)
        if (ABS (psivln-psivl) .lt. 1.0e-6) then
          if (itty .ne. 0) then
            write  (itty, 1045) psivln
 1045       format (' final converged value of psi = ', 1pe16.8)
          end if
          return
        else
          psivl = psivln
          go to 5
        end if
      end if
      return
c
c --- errors
c
*1000 call eval_bicubic_spline (x,nw,y,nh,cspln,nw,nh,x2,y2,pds,ier,1,
**** .                          ilsave,jlsave)
c
 1000 call my_dbcevl1 (x, nw, y, nh, cspln, nw, x2, y2, pds, ier, 1)
      if (itty .ne. 0) then
        write  (itty, 1010) xmin, xmax, x2, ymin, ymax, y2,
     .                      psivl, pds(1), irestart
 1010   format (' open contour in CNTOUR1'        /
     .          ' xmin, xmax, x = ', 3(2x, e16.8) /
     .          ' ymin, ymax, y = ', 3(2x, e16.8) /
     .          ' psivl, pds(1) = ', 2(2x, e16.8), '    irestart =', i4)
      end if
      if (iauto .eq. 1)  go to 1030
      ierr = 1
      if (itty .ne. 0) then
          write  (itty, 2001)  psiaxd, psivl, ipts
 2001     format (' subroutine CNTOUR1 has detected an error:' /
     .          '   psiaxd, psivl, ipts =', 2(2x, 1pe14.6), 2x, i5)
          if (itry .gt. 100) then
            write  (itty, 2002) itry
 2002       format (/
     .      ' The given input value of psi on the plasma boundary' /
     .      ' is so inappropriate that after ', i5, ' tries'       /
     .      ' we still could not find a closed flux surface'       /
     .      ' consistent with the mhdgrid and limiter'             /
     .      ' geometry.  WE ARE GIVING UP.  Please fix your value' /
     .      ' of psilim in the input eqdsk file and try again.'    /)
          end if
 2003     format ( ' ERROR: x1, x2, psivl =', 3(2x, 1pe14.6))
 2004     format ( ' ERROR: y1, y2, psivl =', 3(2x, 1pe14.6))
 2005     format ( ' ERROR: psi2, psivl =',   2(2x, 1pe14.6))
 2006     format ( ' ERROR: netwi, nitr8 =',  2(2x, i5))
 2007     format ( ' ERROR: cmult, xmult, ymult, xsave, ysave, psisave',
     .             ' nsteps =' / 6(2x, 1pe16.8), i5)
 2008     format ( ' ERROR: x1, y1, psi1, x2, y2, psi2, isgn = ' /
     .               6(2x,1pe16.8),i5)
          if (   x1 .eq. x2   )  write (itty, 2003) x1,x2,psivl
          if (   y1 .eq. y2   )  write (itty, 2004) y1,y2,psivl
          if ( psi2 .gt. psivl)  write (itty, 2005) psi2 , psivl
          if (newti .gt. nitr8)  write (itty, 2006) newti, nitr8
          if (cmult .lt. 0.99 )  write (itty, 2007) cmult, xmult, ymult,
     .                                  xsave, ysave, psisave, nsteps
          if (cmult .lt. 0.99 )  write (itty, 2008) x1, y1, psi1, x2,
     .                                  y2, psi2, isgn
      end if
      return
c
 1030 psivl0 = psivl
      psiout = psivl
      dapsi  = psiaxd - psivl0
      dmult  = 0.0005
      if (delta_psi .gt. 0.0) then
        maxtry=20
        dapsi=delta_psi
        dmult=1.0D0 / DFLOAT (maxtry)  ! increments psivl because of psivl0
      end if
      psivl  = psivl0 + dapsi * dmult
      if ((itry .gt. 0) .and. (psiin .ne. 0.0) .and. (iconvg .eq. 1))
     .  psivl = 0.5 * (psiin + psiout)
****  if (iconvg .eq. 0)  psivl = pds(1)  ! HSJ 8/28/96
      iautoc = 1
      if (itty .ne. 0) then
        write  (itty, 1020)  psivl0, psivl
 1020   format (' will try to correct improper setting of psilim' /
     .          ' changing psilim from', e16.8,'  to  ',e16.8)
      end if
      itry = itry + 1
      if (itry .gt. 100) then
        iauto = 0           ! allows error exit return to be taken
        go to 1000          ! can't do it, so quit
      end if
      go to 5               ! go back and try again
c
 1500 cmult = 0.5           ! set cmult < 1 so logic above kicks in
      x1    = xsave
      y1    = ysave
      psi1  = psisave
      isgn  = isgnsave      ! restore search direction
      go to 40
c
      end

      subroutine cntour2 (xaxd,yaxd,psivl,xemin,xemax,yemin,yemax,yxmin,
     .                    yxmax,xymin,xymax,dang,arcl,bperr,dx,dy,xmin,
     .                    xmax,ymin,ymax,iauto,iautoc,xc,yc,ipts,x,nw,y,
     .                    nh,cspln,n2cspln,nh2,itty,iptsm,ierr,bpmag,
     .                    iconvg,delta_psi)
c
 
c
c ----------------------------------------------------------------------
c this routine uses the bicubic spline representation of psi to get
c psi values at an arbitrary point on the MHD grid. the spline
c coefficients must exist (in cspln) before this subroutine is called.
c given (xaxd,yaxd) and a psi value, psivl, generate a contour of
c ordered points,(xc(i),yc(i)), i = 1,ipts, which has (xaxd,yaxd) as
c as an interior point. the contour must fully encircle (xaxd,yaxd).
c the search is limited to a rectangle specified by (xmin,xmax,ymin
c ymax).  dx and dy determine the basic increment in x and y for
c the coarse grid search.  they should be picked so that over the
c distance dx or dy, psi is single valued.  a sliding interval
c search rather than a binary search is used for the coarse grid so
c that non-monotonic psi can be handled.  dx and dy are in meters.
c dx = dy = 0.01 is suggested.  dang is angular step for the rays which
c emmanate from (xaxd,yaxd) in degrees.  dang should be set
c according to shape of equilibrium.  if dang is too small, many
c unnecessary points will be generated in (xc,yc).  if dang is too
c large for bperr(see below), the routine wastes computation time
c reducing the magnitude of dang.  dang = 10 deg. near psilim and
c dang = 30 deg. near psimax is suggested.  for highly elongated
c plasmas control with dang is difficult to set for all contours.
c therefore arcl may be used to set dang internally.  arcl is the
c arclength in meters taken from the current point to get to the
c next point.  set arcl to a large number (e.g. 10 meters) to
c overide this option and use dang only.  the angle increment will
c be the minimum of (arcl/rad,dang).  bperr is relative change in
c poloidal b field between (xc(i),yc(i)) and (xc(i+1),yc(i+1)).
c if the number of elements (xc,yc) exceeds the limit
c bperr is relaxed (up to twice) by a constant increment
c dbperr before an error exit is taken.
c bperr = 0.03 is suggested.  note: if dang yields bperr(computed)<
c bperr(input), then dang is used.  otherwise dang is successively
c reduced by 0.5 until this condition is meet.  (xemin,xemax,yemin,
c yemax) are min and max values of contour.  (yxmin,yxmax), are y
c values at x = xemin and x = xemax respectively.  (xymin,xymax) are x
c values at y = yemin and y = yemax respectively.  iauto is used as an
c "error recovery switch".  that is, if a flux surface is found to
c pass outside the search box an error exit is taken by way of label
c 1000, if iauto = 0.  if iauto = 1, the input value of psi, psivl, is
c increased toward the value of psi at (xaxd,yaxd) until a closed
c surface inside the search box is found.  obviously this option
c makes sense only for the limiter flux surface.  it could be used
c to find the plasma boundary if an appopriate search box were
c defined.  however its primary use here is to slightly modify
c psilim so that a closed limiter flux surface is found.  the option
c is necessary due to the fact some codes generate eqdsks
c using linear interpolation, which differs from the bicubic inter-
c polation used here.  the user is informed that psivl was changed
c by returning iautoc = 1.  if no change occured iautoc = 0.
c if no error ocurred ierr = 0 is returned. if an error ocurred then
c ierr = 1 is returned and the results of this calculation are not valid.
c
c --- input (see description above)
c  xaxd
c  yaxd
c  psivl   psi value for which contour is desired
c  dang
c  arcl   !!!!!!!!!!!!!!!!!!!no longer used
c  bperr   dang,arcl,bperr control point density on contour returned
c  dx
c  dy      defines step length taken
c          input dx and dy as 0.0 to set step length in this subroutine
c  xmin
c  xmax
c  ymin
c  ymax    define the search box
c             NOTE: the search box could be defined as the entire
c                   MHD grid (i.e., xmin = x(1),xmax = x(nw),ymin = y(1),
c                   ymax = y(nh)) However this will normally not work
c                   satisfactorily because there are local minima and
c                   maxima around the f coils which would confuse the
c                   contour tracing and require many more calculations.
c                   (these additional calculations are not done here!)
c                   The best way to define the search box is to set it
c                   equal to a rectangle,within which the plasma must
c                   reside under all circumstances. The only such
c                   rectangle to suggest itself is the one defined by
c                   the extremes of the limiter. Hence it is suggested
c                   that xmin be set to the min x value of the limiter
c                   points vector,xmax the max x value of limiter points
c                   vector and similarly for ymin,ymax.
c  iaouto  = 1 automatic error recovery, = 0 no recovery
c             (iauto may be set to 0 in this subroutine
c              if a failure occurs)
c  x(1..nw)
c  y(1..nh)   the MHD grid vectors
c  cspln(n2cspln,nw,nh2)  bicubic spline coefficient array
c  itty      iounit for diagnostic messages. set to 0 to suppress.
c  iptsm     max allowed number of points to be returned in
c            vectors xc,yc
c
c delta_psi        this is the psi grid spacing. Define delta_psi as
c                  psi(penultimate value)-psi(last(i.e., edge) value)
c                  If the value of psi to be traced leads to an error
c                  and iauto=1 then delta_psi is used to control the
c                  amount by which the value of psi is adjusted. Otherwise
c                  delta_psi is not used,so you set its value accordingly.
c
c --- output:
c  xemin
c  xemax
c  yemin
c  yemax
c  yxmin
c  yxmax
c  xymax
c  xymin    limits of contour (see above)
c  xc(1,..ipts)
c  yc(1,..ipts)   the ipts contour points found
c  bpmag(1,..ipts)    magnitude of poloidal bfield at xc,yc
c  ipts
c  ierr   error flag, = 0 no error, = 1 error, output not usable.
c             if ierr .gt. 1 then ierr equals the number of points found
c             on the contour up to this point. Since this number will
c             be exceeded if the code continues,the subroutine must
c             stop. The error condition on bperr has allready been relaxed
c             before this condition causes the subroutine to terminate.
c             Hence in order to continue this subroutine could be called
c             again with a larger value of arcl. (Note: in practice we have
c             not found this to be a problem,provided that iptsm is large
c             enough to begin with toi give an reasonable representation
c             of the contour).
c ------------------------------------------------------------------ HSJ
c
      USE replace_imsl,                   ONLY : my_dbcevl1

      implicit  integer (i-n), real*8 (a-h, o-z)
      include 'imsl.i'
c
      dimension  pds(6), xc(*), yc(*), x(nw), y(nh),
     .           cspln(n2cspln,nw,nh2), bpmag(*)
      data       pi, piov2, piov4, fpiov4, spiov4
     .  /3.141592654,1.570796327,0.7853981634,3.926990818,5.497787145/
      data       twopi,tpiov4, tpiov2, nitr8
     .  /6.283185307,2.356194491,4.712388904,20/
c
      imslmd   = 'cntour2'
      itry     = 0
      irestart = 0
      psiin    = 0.0
      ierr     = 0
      iautoc   = 0
      bperrsav = bperr
      dbperr   = bperr / 2.0
      iprintit = 77
      iprintito = 0
      iptsmon  = 64
      miunit = 0 !avoids dec f90 complaint
      if (iprintit .eq. 1) then
          iprintito = 1
cLLL: had been:
c     open (unit = iprintit, file = 'cntourdebug.dat',
c    .    status = 'UNKNOWN')
cLLL: note write's to miunit below, but miunit is undefined. propose:
          call getioun(miunit,77)
          open (unit = miunit, file = 'cntourdebug.dat',
     .    status = 'UNKNOWN')
      end if
c
c --- get psi at (xaxd,yaxd)
c
      call my_dbcevl1 (x,nw,y,nh,cspln,nw,xaxd,yaxd,pds,ier,1)
      if (ier .ne. 0) then
          if (itty .ne. 0)  write (itty, 4) ier
    4     format (' subroutine CNTOUR2 detects an error:' /
     .            ' DBCEVL1 returned ier =', i5,
     .            ' while getting psi on axis')
          call STOP ('subroutine CNTOUR2: cnt1 problem', 160)
      end if
      psiaxd = pds(1)
      if (psiaxd .lt. psivl) then
          if (itty .ne. 0)  write (itty, 3)  psiaxd, psivl
    3     format (' subroutine CNTOUR2 detects an error :'     /
     .            ' psi on magnetic axis ',1pe12.6             /
     .            ' psi contour to be found ',1pe12.6          /
     .            ' must have max psi on axis,non recoverable' /
     .            ' ERROR, program must stop')
          call STOP ('subroutine CNTOUR2: cnt2 problem', 161)
      end if
    5 xemin =  1.0e10
      xemax = -1.0e10
      yemin =  1.0e10
      yemax = -1.0e10
      xymin =  0.0
      xymax =  0.0
      yxmin =  0.0
      yxmax =  0.0
      ipts  =  0
      irestart = irestart + 1
      tstart   = fpiov4
      tend     = tstart + twopi
      thet     = tstart
      dthet0   = twopi * dang / 360.0
      dthetmin = 1.0e-3    ! min increment in theta from one ray to next
      if (dthetmin .gt. dthet0)  dthetmin = 0.5 * dthet0
      dthet    = 0.0
      dxx      = dx
      dyy      = dy
      if (dxx .eq. 0.0)  dxx = x(3)-x(1)    ! 2 times grid spacing is ok
      if (dyy .eq. 0.0)  dyy = y(3)-y(1)
      dxmin = 0.05*(x(2)-x(1)) ! subroutine stops if dxx becomes < dxmin
      dx    = dxx
      dy    = dyy
      if (dxmin .ge. dxx)  dxmin = 0.5 * dxx
      serrt = 3.5e-06
      derrt = 0.5e-07
c
c --- serrt is absolute error convergence criteria for Newton's method below
c --- loop over theta from tstart to tend( = tstart+twopi)
c
   10 thet   = thet + dthet
      if (thet .gt. twopi .and. tstart .ne. 0.0) then
        thet = thet - twopi
        tend = tend - twopi
      end if
      thet = MIN (thet, tend)
      if (thet .eq. tend)  go to 200
      dxx = dx
      dyy = dy
c
c --- get equation of ray emanating from (xaxd, yaxd)
c
      if (( piov4 .le. thet) .and. (thet .le. tpiov4))  go to 20
      if ((fpiov4 .le. thet) .and. (thet .le. spiov4))  go to 20
c
c --- y as a function of x            y = a*x+bincp
c
      isgn = -1
      if ((thet .lt. piov4) .or. (thet .gt. spiov4))  isgn = 1
      a = TAN (thet)
      iflg = 0
      bincp = yaxd-a*xaxd
      go to 30
c
c --- x as a function of y            x = a*y+bincp
c
   20 isgn = 1
      if (thet .gt. pi)  isgn = -1
      if (isgn .eq. -1)  go to 22
      thet1 = piov2 - thet
      if (thet .gt. piov2)  thet1 = twopi - ABS (thet1)
      go to 25
   22 thet1 = tpiov2 - thet
      if (thet .gt. tpiov2)  thet1 = pi - ABS (thet1)
   25 a     = TAN (thet1)
      iflg  = 1
      bincp = xaxd - a*yaxd
   30 continue
c
c --- now have y = a*x+bincp    (iflg=0)   or
c --- x = a*y+bincp             (iflg=1)
c
      x1   = xaxd
      y1   = yaxd
      cost = COS (thet)
      sint = SIN (thet)
      psi1 = psiaxd
c
c --- sliding interval search. max width of interval ~1.41*(dx or dy)
c
   40 if (iflg .eq. 1)  go to 50
c
c --- search in x
c
      xmult = 1.0
   41 x2 = x1+isgn*dxx*xmult         ! increment x
      if (iprintit .eq. 1 .and. ipts .eq. iptsmon)
     .    write (miunit, 3001) psivl,x1,y1,x2,y2,xmin,xmax,ymin,ymax
 3001 format (2x,"psivlx,x1,y1,x2,y2,xmin,xmax,ymin,ymax  =" /
     .        5(2x,1pe16.8),4(2x,1pe16.8))
      if (x2 .eq. x1)  go to 1000    ! psivl is not within reach
      x2 = MAX (x2,xmin)
      x2 = MIN (x2,xmax)             ! limit the search
      if (x2 .eq. x1)  go to 1000    ! psivl is not within reach
      y2 = a*x2+bincp                ! get corresponding value of y
      if (y2 .lt. ymin .or. y2 .gt. ymax) then
        xmult = 0.5 * xmult
        go to 41
      end if
      go to 60
c
c --- search in y
c
   50 ymult = 1.0
   51 y2 = y1+isgn*dyy*ymult          ! increment y
      if (iprintit .eq. 1 .and. ipts .eq. iptsmon)
     .        write(miunit,3000)psivl,x1,y1,x2,y2,
     .                            xmin,xmax,ymin,ymax
c
 3000         format(2x,"psivly,x1,y1,x2,y2,xmin,xmax,ymin,ymax  =",/,
     .               5(2x,1pe16.8),4(2x,1pe16.8))
      if (y2 .eq. y1)  go to 1000     ! psivl is not within reach
      y2 = MAX (y2,ymin)
      y2 = MIN (y2,ymax)              ! limit the search
      if (y2 .eq. y1)  go to 1000     ! psivl is not within reach
      x2 = a*y2+bincp                 ! get corresponding value of x
      if (x2 .lt. xmin .or. x2 .gt. xmax) then
          ymult = ymult * 0.5
          go to 51
      end if
c
   60 call my_dbcevl1 (x,nw,y,nh,cspln,nw,x2,y2,pds,ier,1)
      if (ier .ne. 0) then
        write (6, *)  'CNTOUR2 reports, call to DBCEVL1 ,ier =', ier
        write (6, *)  'arguments x2,y2 =',x2,y2
        write (6, *)  'xgrid(1),xgrid(nw) =',x(1),x(nw)
        write (6, *)  'ygrid(1),ygrid(nh) =',y(1),y(nh)
        call STOP ('subroutine CNTOUR2: unspecified problem', 162)
      end if
      psi2 = pds(1)
      dpsi = (psivl-psi1)*(psivl-psi2)
      if (iprintit .eq. 1 .and. ipts .eq. iptsmon)
     .        write(miunit,3002)dpsi,psivl,psi1,psi2
 3002 format(2x,"dpsi,psivl,psi1,psi2 =",/,
     . 4(2x,1pe16.8))
      if (dpsi .le. 0.0)  go to 70
      x1   = x2
      y1   = y2
      psi1 = psi2
      go to 40
c
c --- now have psivl between psi1 and psi2,converge using newton-raphson
c
   70 newti = 0
      if (iflg .eq. 1)  go to 75
      xn = x1+isgn*dxx * 0.5
      yn = a*xn+bincp
      go to 80
   75 yn = y1+isgn*dyy * 0.5
      xn = a*yn+bincp
   80 call my_dbcevl1 (x,nw,y,nh,cspln,nw,xn,yn,pds,ier,3)
      if (ier .ne. 0) then
        newti = nitr8
        go to 76
      end if
      dpsids = pds(2)*cost+pds(3)*sint
      dpsi   = pds(1)-psivl
      serr   = -dpsi/dpsids
      if (ABS (serr) .lt. serrt)  go to 90
      if (psivl .eq. 0.0 .and. ABS (dpsi) .lt. derrt)  go to 90
      if (psivl .ne. 0.0) then
          if (ABS (dpsi/psivl) .lt. derrt)  go to 90
      end if
      delx  = serr*cost
      dely  = serr*sint
      xn    = xn+delx
      yn    = yn+dely
      newti = newti+1
   76 if (newti .ge. nitr8) then
        if (dxx .le. dxmin)
     .    call STOP ('subroutine CNTOUR2: newti problem', 163)
        dxx = 0.5 * dxx
        dyy = 0.5 * dyy
      if (iprintit .eq. 1 .and. ipts .eq. iptsmon)
     .        write(miunit,3004)
 3004 format(2x,"returning to label 40 with dxx,dyy = ",/,
     .       2(2x,1pe16.8))
        go to 40        ! 6/11/93
      end if
c
      go to 80
c
c --- end of newton iteration
c --- check for sufficient accuracy in point spacing as determined by thet
c --- accuracy test is based on a relative error in poloidal b field of bperr
c
   90 bp2 = SQRT (pds(2)**2 + pds(3)**2) / xn
      if (thet .eq. tstart)  go to 100    ! first point of contour
c
c --- if relative change in bp is less than bperr store the point just found:
c
      if (iprintit .eq. 1 .and. ipts .eq. iptsmon)
     .        write(miunit,3005)bp2,bp1,bperr,dthetmin,dthet
 3005 format(2x,"bp2,bp1,bperr,dthetmin,dthet =",/,
     .5(2x,1pe16.6))
      if (ABS (bp2-bp1) / MAX (bp2,bp1) .lt. bperr)  go to 100
c
c --- if dthet is currently at its smallest allowed value of dthetmin
c --- and relative change in bp is still greater than bperr continue anyway:
c
      if (dthet .eq. dthetmin)  go to 100
c
c --- spacing too large for grad psi. decrease theta and try again
c
      thet  = thet - dthet          ! reset to previous value
      dthet = dthet * 0.5             ! take half the angle increment
      dthet = MAX (dthet, dthetmin) ! dthet lower limit must be observed
      go to 10        ! go back and try with a new ray closer to old one
  100 bp1   = bp2
      ipts  = ipts+1
      if (ipts .le. iptsm-1)  go to 150
c
c --- we ran out of storage for the contour points. assuming iptsm is
c --- set to a reasonable value this should not happen unless bperr is
c --- set unrealistically small. we increase bperr by dbperr and try
c --- to generate the contour again with the relaxed error requirement.
c
      if (bperr-bperrsav .ge. 2.0*dbperr)  go to 110
      bperr = bperr+dbperr
      go to 5
  110 if (itty .ne. 0) then
        write (itty, 1040)  iptsm, dthet0, bperr, psivl
 1040   format
     . ('  more points on plasma '                                     /
     .  '  boundary were found than dimension of xcontr,ycontr allows' /
     .  '  allowed number = ',i5                                       /
     .  '  adjustment in point spacing will be made',
     .  '  dthet0,bperr,psivl =',5(2x,1pe12.4))
      end if
c
      ierr = ipts
      return
c
  150 xc(ipts) = xn   ! found new point; save it and set up for next ray
      yc(ipts) = yn
      bpmag(ipts) = bp2
      xemin = MIN (xemin, xn)
      xemax = MAX (xemax, xn)
      if (xemax .eq. xn)  yxmax = yn
      if (xemin .eq. xn)  yxmin = yn
      yemin = MIN (yemin, yn)
      yemax = MAX (yemax, yn)
      if (yemax .eq. yn)  xymax = xn
      if (yemin .eq. yn)  xymin = xn
      if (ABS (bp2-bp1) / MAX (bp2,bp1) .le. 0.25 * bperr)
     .dthet = 2.0 * dthet
      dthet = MIN (dthet, dthet0)        ! don't let dthet get too large
      dthet = MAX (dthet, dthetmin)      ! don't let dthet get too small
      go to 10                           ! continue with the next ray
c
c --- done with this contour. close contour
c
  200 ipts        = ipts+1
      xc(ipts)    = xc(1)
      yc(ipts)    = yc(1)
      bpmag(ipts) = bpmag(1)
      bperr       = bperrsav
c
      if (itry .ne. 0 .and. iconvg .eq. 1) then
        psiin  = psivl
        psivln = 0.5 * (psiout+psiin)
        if (ABS (psivln-psivl) .lt. 1.0e-6) then
          if (itty .ne. 0) then
            write  (itty, 1045) psivln
 1045       format (' final converged value of psi = ', 1pe16.8)
          end if
          return
        else
          psivl = psivln
          go to 5
        end if
      end if
      return
c
c --- errors
c
 1000 call my_dbcevl1 (x, nw, y, nh, cspln, nw, x2, y2, pds, ier, 1)
      if (iprintito .eq. 1) then
        iprintit = miunit
        call giveupus(iprintit)
        close (unit = iprintit)
      end if
      if (itty .ne. 0) then
        write  (itty, 1010) xmin, xmax, x2, ymin, ymax, y2,
     .                      psivl, pds(1), irestart
 1010   format (' open contour in CNTOUR2'       /
     .          ' xmin, xmax, x = ', 3(2x, e16.8) /
     .          ' ymin, ymax, y = ', 3(2x, e16.8) /
     .          ' psivl, pds(1) = ', 2(2x, e16.8), '    irestart =', i4)
      end if
      if (iauto .eq. 1)  go to 1030
      ierr   = 1
      return
 1030 psivl0 = psivl
      psiout = psivl
      dapsi  = psiaxd - psivl0
      dmult  = 0.0005
      if (delta_psi .gt. 0.0) then
        maxtry=99
        dapsi=delta_psi
        dmult = 1.0D0 / DFLOAT (maxtry)! increments psivl because of psivl0
      end if
      psivl = psivl0 + dapsi * dmult
      if ((itry .gt. 0) .and. (psiin .ne. 0.0) .and. (iconvg .eq. 1))
     .  psivl = 0.5 * (psiin+psiout)
****  if (iconvg .eq. 0)  psivl = pds(1)  ! HSJ 8/28/96
      iautoc = 1
      if (itty .ne. 0)
     .write  (itty, 1020)  psivl0, psivl
 1020 format (' boundary search, will change psilim from' /
     .          e16.8, '  to  ', e16.8, '  and try again')
      itry = itry + 1
      if (itry .gt. 50) then
        iauto = 0
        go to 1000        ! can't do it, quit
      end if
      go to 5             ! go back and try again
c
      end

      subroutine concek (p, po, pmax, pmin,
     .                   nx, ny, toleq, ind, s, imax, jmax)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c     check convergence of equilibrium solution (called by FREEBDRY)
c ----------------------------------------------------------------------
c
      dimension  p(nx,ny), po(nx,ny)
c
      s   = 0.0
      del = 0.0
c
      do   i=1,nx
        do j=1,ny
          dif = ABS (p(i,j) - po(i,j))
          del = MAX (del, dif)
          if (del .eq. dif) then
            imax = i
            jmax = j
          end if
        end do
      end do
c
      s   = del / ABS (pmax)
      ind = 0
      if (s .lt. toleq)  ind = 1
      return
c
      end

      subroutine concur (ic, tolc)
c
      USE param
      USE io
      USE solcon
      USE soln
      USE mhdpar   
      USE numbrs      
      USE mesh
      USE soln2d
      USE rhog
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c This subroutine compares current profile quantities before and after
c the last equilibrium calculation and calculates the change in pprim.
c NOTE: for irguess = +99 tolc will equal 1.0 on the first pass through
c       this routine because ppold is not set initially.
c ----------------------------------------------------------------------
c
c      include 'param.i'
c      include 'mhdpar.i'
c      include 'io.i'
c      include 'mesh.i'
c      include 'numbrs.i'
c      include 'rhog.i'
c      include 'solcon.i'
c      include 'soln.i'
c      include 'soln2d.i'
c
      tolc  = 0.0
      denom = 0.0
      do j=1,nj
        tolc  = tolc  + (pprim(j)-ppold(j))**2
        denom = denom +  pprim(j)**2
      end do
      tolc  = tolc/denom
      timms = time * 1.0e3
      icc   = ic - 1
      call headeq (6)
      write  (nout , 50)  timms, icc, tolc
      write  (nitre, 50)  timms, icc, tolc
   50 format (/ ' subroutine CONCUR reports quantities changed by new ',
     .           'equilibrium' //
     .          ' time = ', f10.2, ' msec', 5x, 'ic=', i3, 5x, 'tolc=',
     .                      e12.3 //
     .            5x,  8x, 8h   rho   , 8x, 8x, 8h    psi /
     .            5x, 12h   old       , 12h   new      ,
     .                12h   old       , 12h   new      )
      do j=1,nj,jprt
        write  (nitre, 70)  j, rold(j), r(j), psiold(j), psir(j)
        write  (nout , 70)  j, rold(j), r(j), psiold(j), psir(j)
   70   format (1x, i5, 8e12.4)
      end do
c
      write (nout , 80)
      write (nitre, 80)
   80 format (/ 5x,  8x, 8h   cur   , 8x, 8x, 8h  pprim  , 8x,
     .               8x, 8h ffprim /
     .          5x, 12h   old       , 12h   new      ,
     .              12h   old       , 12h   new      ,
     .              12h   old       , 12h   new      )
c
      do j=1,nj,jprt
        write (nitre, 70)  j,curold(j),curden(j),
     .                     ppold(j),pprim(j),ffpold(j),ffprim(j)
        write (nout , 70)  j,curold(j),curden(j),
     .                     ppold(j),pprim(j),ffpold(j),ffprim(j)
      end do
      return
c
      end

      subroutine curcalc (rbp, fcap, hcap, gcap, r, curden,
     .                    curpar_soln, nj)
c
c
c ----------------------------------------------------------------------
c --- subroutine calculates the current density ,<Jphi*R0/R>,using Ampere's
c --- law in the form:
c           <Jphi*R0/R> = (c/4pi)*(1.0/(h*r)*d(g*h*r*bp)/dr
c --- The parallel current density
c     < Jtotal  dot B/Bt0> = (c/4pi)*(1./(F**2H*r)*d(FGHrbp)/dr
c --- it is assumed that d<Jphi*R0/R>/dr = 0 at the magnetic axis.
c --- We enforce this condition below by extrapolating the cubic polynomial
c --- fit of curden near the magnetic axis (see subroutine CUBICEXTRP)
c
c --- input
c argument list:
c  rbp(j)                      j = 1,..nj,rbp(j) = f(j)*g(j)*h(j)*r(j)*bp(j)
c                              in gauss cm
c
c  fcap(j)                     j = 1,2..nj,fcap(j) = f(psilim)/f(psi)
c
c  hcap(j)                     j = 1,2..nj,hcap(j) = fcap(j)/<R0**2/R**2>
c  gcap(j)                     j = 1,2,..nj,gcap(j) = <(grad rho R0/R)**2>
c
c  r(j)                        j = 1,2..nj,r(j) = rho value (in cm)
c
c  nj                          the grid size
c
c  INCLUDE files:
c  INCLUDE file constnts.i:
c  pi
c
c  INCLUDE file storage:
c  xdum(j)
c  ydum(j)
c  zdum(j)
c  wdum(j)                  j = 1,2..nj temporary work storage
c
c --- output
c  argument list:
c  curden(j)                j = 1,2..nj, the current density,amps/cm**2,
c                           on the r(j) grid,with zero gradient at r = 0
c  curpar_soln(j)           parallel current, see below
c
c ------------------------------------------------------------------ HSJ
c
      USE param
      USE constnts
      implicit  integer (i-n), real*8 (a-h, o-z)
c      include 'param.i'
c      include 'constnts.i'
      include 'storage.i'
c
      dimension    rbp(*), fcap(*), hcap(*), r(*), curden(*),
     .             curpar_soln(*), gcap(*)
      dimension    ghrbp (kj), dfcapdr(kj), dghrbpdr(kj)
      equivalence (ghrbp ( 1), xdum(1)),   (dghrbpdr( 1), ydum(1))
      equivalence (dfcapdr(1), zdum(1)) !note rdum is used in curperp

c
c --- get g*h*rho*bp
c
      do j=1,nj
        ghrbp(j) = rbp(j) / fcap(j)
      end do
c
c --- differentiate ghrbp
c
      call difydx(r,ghrbp,dghrbpdr,nj)
      const = 1.0 / (0.4*pi)    ! converts to amps/cm**2
c
c --- get curden in amps/cm**2
c
      do j=2,nj 
        curden(j) = (dghrbpdr(j)/(hcap(j)*r(j)))*const
      end do


c
c --- get current density at magnetic axis,consistent with
c --- d<Jphi*R0/R>/dr  = 0 at r=0
c
      call cubicextrp (curden(2), curden(3), curden(4),
     .                 r(2), r(3), r(4), curden(1), 2)


c --- integrate curden:
      sum = 0.0
      DO j=1,nj-1
         sum = sum +0.5*(curden(j)*hcap(j)*r(j)
     .         +curden(j+1)*hcap(j+1)*r(j+1))*(r(j+1)-r(j))
      ENDDO




c     make sure curden is normalized to the total current
c     (uses rbp(nj):
      call curnorm

c
c     get the parallel current , curpar_soln,given  by
c     curpar_soln =
c     <J*B/Bt0> = < J*R0/R>/fcap +(c/4pi)*(gcap/fcap**2)*d/dr(fcap)
c
c      call difydx (r, fcap, dfcapdr, nj)
c      do j=1,nj
c        curpar_soln(j) = curden(j)/fcap(j) +
c     .                   const*(gcap(j)/(fcap(j)*fcap(j)))*dfcapdr(j)
c      end do
      call difydx (r, rbp, dfcapdr, nj)
      do j=2,nj
        curpar_soln(j) = const/(fcap(j)*fcap(j)*hcap(j)*r(j))
     .                                                 *dfcapdr(j)
      end do
      curpar_soln(1)=curpar_soln(2)
c     get the contribution to the toroidal current due to perpendicular cur.
      call curperp
      return
c
      end

      subroutine curdenmhd
c
c
c ----------------------------------------------------------------------
c --- load pcurrent
c  input
c  INCLUDE file param.i        for parameters used in other INCLUDE files
c  INCLUDE file numbrs.i
c       nj
c  INCLUDE file mhdpar.i
c       nw,nh
c  INCLUDE file mhdgrid.i
c       rmhdgrid
c       zmhdgrid
c  INCLUDE file zerocom.i
c       wzero(nwh)
c  INCLUDE file mhdcom.i
c       psi(i,j)
c       pcurrent(nwh) (output,amps/m**2)
c  INCLUDE file rhog.i
c       psir(nj)
c       pprim(nj)
c       ffprim(nj)
c  INCLUDE file constnts.i
c       u0
c  INCLUDE file bicube.i
c       wnoperm(nwork)          temporary storage
c ----------------------------------------------------------------------
c
      USE param
      USE mhdpar
      USE mhdgrid 
      USE numbrs        
      USE constnts
      USE rhog
      USE mhdcom
      USE bicube
      implicit  integer (i-n), real*8 (a-h, o-z)

c
      dimension                curden2d(nw,nh)
c      equivalence (wnoperm(1), curden2d( 1, 1))
c
      jsymetric = 1
      call curdencalc2 (curden2d,nw,nh,rmhdgrid,wzero,psi,
     .                  psir,nj,u0,totc,pprim,
     .                  ffprim,jsymetric)
      do j=1,nh
         do i=1,nw
            k0 = (i-1)*nh+j
            pcurrent(k0) = curden2d(i,j)
         end do
      end do
c
      darea = (rmhdgrid(2)-rmhdgrid(1))
      darea = darea*(zmhdgrid(2)-zmhdgrid(1))
      totc  = totc*darea
****  write (6, *)  'curdenmhd total current =', totc
      return
c
      end

      subroutine curnorm
c ----------------------------------------------------------------------------
c     normalize curden  to total toroidal current
c
c-----------------------------------------------------------------------------
      USE param
      USE soln
      USE numbrs
      USE mesh
      USE geom
      USE constnts
      implicit  integer (i-n), real*8 (a-h, o-z)
c      include 'param.i' !kj
c      include 'constnts.i' ! twopi
c      include 'geom.i'  !gcap,hcap,
c      include 'mesh.i'  !r(1..kj)
c      include 'numbrs.i'!nj
c      include 'soln.i'  !rbp, curden (amps/cm**2)

      real *8 ,dimension(:),allocatable :: xdum,ydum
      allocate(xdum(kj),ydum(kj))
      xdum = hcap * r * curden
      call trap2(r,xdum,ydum,nj)
      xnorm = 5.*rbp(nj)/(twopi*ydum(nj))

      curden = curden * xnorm 


      deallocate(xdum)
      deallocate(ydum)
      return
      end

      


      subroutine curperp
c
c ----------------------------------------------------------------------
c --- calculate the contribution to toroidal current due to
c --- the perpendicular and Pfirsh Schluter currents.
c ---------------------------------------------------------HSJ-8/28/98--
c
      USE param
      USE soln
      USE numbrs
      USE mesh
      USE machin
      USE geom
      USE constnts
      USE rhog
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c      include 'param.i'
c      include 'constnts.i'  ! pi
c      include 'geom.i'      ! r2cap =<R0**2/R**2>,fcap,etc,
c      include 'numbrs.i'    ! nj
c      include 'machin.i'
c      include 'mesh.i'      ! r(j)
c      include 'rhog.i'      ! press(j)
c      include 'soln.i'      ! cur_tor_ps_soln(j)
      include 'storage.i'
c
      dimension   dpressdr(kj)
      equivalence(dpressdr(1), rdum(1))
c
      ten = 10.0
c
c     get the pressure (returned in rhog.i::press(1-nj) )
c
      call pressr (0, 1)        ! get total press in units of dyne/cm**2
      call difydx (r, press, dpressdr, nj)
      do j=2,nj
        factor_ = (1.0 - (1.0/fcap(j))**2*r2cap(j)/bsq_avg_cap(j))
        bp0_    = rbp(j)/(r(j)*fcap(j)*gcap(j)*hcap(j))  ! gauss
        cur_tor_ps_soln(j) =
     .          -(dpressdr(j)/bp0_)*factor_ * ten ! amps/cm**2
      end do
      cur_tor_ps_soln(1) = cur_tor_ps_soln(2)
      return
c
      end

      subroutine curpro (ii)
c
c
c ----------------------------------------------------------------------
c This subroutine converts pprim and ffprim to the transport r mesh
c and then calculates the (absolute value of) current density.
c NOTE: curden = flux surface average of toroidal current density
c       time rmajor/r in units of amps/m**2
c ----------------------------------------------------------------------
c
      USE param
      USE soln
      USE mhdpar   
      USE numbrs    
      USE machin
      USE geom
      USE constnts
      USE rhog
      implicit  integer (i-n), real*8 (a-h, o-z)
c      include 'param.i'
c      include 'constnts.i'
c      include 'mhdpar.i'
c      include 'geom.i'
      include 'imsl.i'
c      include 'machin.i'
c      include 'numbrs.i'
c      include 'rhog.i'
c      include 'soln.i'     !njcur
c
      if (njcur .ne. 0)  return
c
      imslmd = 'curpro'
      if (ii .eq. 0) then
        call intrp (1, 1, psiold,  ppold, nj, psir,  pprim, nj)
        call intrp (1, 1, psiold, ffpold, nj, psir, ffprim, nj)
      end if
      do j=1,nj
        curden(j) = -rmajor*pprim(j) - ffprim(j)*r2cap(j)/(u0*rmajor)
      end do
      return
c
      end

      subroutine dbcevl1 (x, nx, y, ny, c, ic, xl, yl, pds, ier, icalc)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c --- this is a special version of dbcevl, used by the ONETWO suite
c --- of subroutines. It differs from the IMSL version in the calling
c --- arguments and the calculations below are reordered for efficiency
c --- That is, many times we don't need all of the pds array. This
c --- version returns only those values of pds specifically requested
c --- using the argument icalc. This leads to about at 50% improvement
c --- in execution time. --- HSJ
c
      integer  nx,ny,ic,ier,i,j,k,km1,kp1,kp2,lxpl,lx,ly,l,lxp1
      real*8   x(*),y(*),c(2,ic,*),xl,yl,pds(6),
     .         hx,hy,sux(2),suy(2),su(2),svx(2),sv(2),sxy(2),
     .         u,v,spln0,spln1,spln2,s0,sh,sp0,sph,h,d
c
      spln0(s0,sh,sp0,sph,h,d) = s0+d*(h*sp0+d*(3.0*(sh-s0)-
     . (sph+2.0*sp0)*h+d*(2.0*(s0-sh)+(sph+sp0)*h)))
c
      spln1(s0,sh,sp0,sph,h,d) = sp0+d*(6.0 * (sh-s0)/h-2.0*
     . (sph+2.0*sp0)+3.0*d*(2.0*(s0-sh)/h+(sph+sp0)))
c
      spln2(s0,sh,sp0,sph,h,d) = 6.0 * (sh-s0)/h**2-2.0*
     . (sph+2.0*sp0)/h+d*(2.0*(s0-sh)/h**2+(sph+sp0)/h)*6.
c
      data  lx, ly /0, 0/
c
      ier = 0
c
c --- correlated table search for xl
c
      call tableintrp (x, nx, xl, lx)
      if (lx  .eq. 0 )  ier = 33
      if (lx  .eq. nx)  ier = 35
      if (ier .ne. 0 )  go to 100
c
c --- correlated table search for yl
c
      call tableintrp (y, ny, yl, ly)
      if ( ly .eq. 0 )  ier = 34
      if ( ly .eq. ny)  ier = 36
      if (ier .ne. 0 )  go to 100
c
c --- replaces the following linear search
c
****  if (xl .lt. x(1))  ier = 33
****  do 5 i=2,nx
****     lx = i-1
****     if (xl .le. x(i))  go to 10
****5 continue
****  ier = 35
***10 if (yl .lt. y(1))  ier = 34
****  do 15 j=2,ny
****     ly = j-1
****     if (yl .le. y(j))  go to 20
***15 continue
****  ier = 36
***20 if (ier .ne. 0)  go to 100
c
      lxp1 = lx+1
      hx   = x(lxp1)-x(lx)
      hy   = y(ly+1)-y(ly)
      u    = (xl-x(lx))/hx
      v    = (yl-y(ly))/hy
      k    = 2*ly
      kp1  = k+1
      kp2  = k+2
      km1  = k-1
c
      do 25 l=1,2
         lxpl = lx-1+l
         i    = 2*(ly-1+l)
         j    = i-1
         sv(l) = spln0(c(1,lxpl,km1),c(1,lxpl,kp1),c(1,lxpl,k),
     .   c(1,lxpl,kp2),hy,v)
         svx(l) = spln0(c(2,lxpl,km1),c(2,lxpl,kp1),c(2,lxpl,k),
     .   c(2,lxpl,kp2),hy,v)
      if (icalc .lt. 3)  go to 25                    ! rearranged by HSJ
         su(l) = spln0(c(1,lx,j),c(1,lxp1,j),c(2,lx,j),
     .   c(2,lxp1,j),hx,u)
         suy(l) = spln0(c(1,lx,i),c(1,lxp1,i),c(2,lx,i),
     .   c(2,lxp1,i),hx,u)
      if (icalc .lt. 4)  go to 25
         sux(l) = spln1(c(1,lx,j),c(1,lxp1,j),c(2,lx,j),
     .   c(2,lxp1,j),hx,u)
         sxy(l) = spln1(c(1,lx,i),c(1,lxp1,i),c(2,lx,i),
     .   c(2,lxp1,i),hx,u)
   25 continue
c
      pds(1) = spln0(sv(1),sv(2),svx(1),svx(2),hx,u)
      if (icalc .eq. 1)  return ! special routine only calculates pds(1)
      pds(2) = spln1(sv(1),sv(2),svx(1),svx(2),hx,u)
      if (icalc .eq. 2)  return
      pds(3) = spln1(su(1),su(2),suy(1),suy(2),hy,v)
      if (icalc .eq. 3)  return
      pds(4) = spln1(sux(1),sux(2),sxy(1),sxy(2),hy,v)
      if (icalc .eq. 4)  return
      pds(5) = spln2(sv(1),sv(2),svx(1),svx(2),hx,u)
      if (icalc .eq. 5)  return
      pds(6) = spln2(su(1),su(2),suy(1),suy(2),hy,v)
  100 if (ier .gt. 0) then
****    call uertst1 (ier, 'dbcevl')
      end if
      return
c
      end

      subroutine difydx (x, y, yprime, nj)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c --- subroutine differentiates y w.r.t. x, for nj values of x
c --- and returns result in yprime.
c --- we use the basic quadratic langrangian interpolation formula
c --- for non-uniform mesh spacing in x to model y and then
c --- differentiate this model analytically. This yields central
c --- differences, accurate to order of the mesh spacing squared,
c --- appropriate forward (at j = 1) and backward (at j = nj) formulae
c --- are used.
c --- input
c  x(1..nj)                   independent variable
c  y(1..nj)                   dependent variable
c  nj                         size of x,y
c --- output
c  yprime(1..nj)              dy/dx
c ------------------------------------------------------------------ HSJ
c
      dimension x(nj), y(nj), yprime(nj)
c
      dxm       = x(2) - x(1)
      dxp       = x(3) - x(2)
      yprime(1) = -y(1)*(2.0*dxm+dxp)/(dxm*(dxm+dxp))
     .           + y(2)*(dxm+dxp)/(dxm*dxp)
     .           - y(3)*(dxm/dxp)/(dxp+dxm)
c

      do j=2,nj-1
        dxr = dxp/dxm
        yprime(j) = -y(j-1)*dxr/(dxp+dxm)
     .    +y(j)*(1.0/dxm-1.0/dxp)
     .    +y(j+1)/(dxr*(dxp+dxm))
        if (j .ne. nj-1) then
          dxm = dxp
          dxp = x(j+2)-x(j+1)
        end if
      end do
c
      yprime(nj) = y(nj-2)/(dxm+dxp)
     .           - y(nj-1)*(dxp+dxm)/(dxp*dxm)
     .           + y(nj)*(2.0*dxp+dxm)/(dxp*(dxp+dxm))
c
c     yprime(1 ) = (y(2)-y(1))/(x(2)-x(1))
      yprime(nj) = (y(nj)-y(nj-1))/(x(nj)-x(nj-1))
      return
c
      end

      subroutine drive (ii, ic)
c
      USE param
      USE neut
      USE io
      USE solcon
      USE soln
      USE contour
      USE mhdpar
      USE mhdgrid
      USE nub  
      USE rf    
      USE extra 
      USE numbrs
      USE machin
      USE geom
      USE constnts
      USE soln2d
      USE bd_condtn
      USE psig
      USE rhog
      USE mhdcom
      USE bicube
      USE shapctr
      USE flxav
      USE neo2dp
      USE etc
      USE iterdbmd,                      ONLY : smth_mhd_parms
      USE smthspline,                    ONLY : smooth_values
      USE gpsi
      USE mhdbcdtn
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c driver for D s
c ----------------------------------------------------------------------
c

c      include 'small.i'

c      include 'zerocom.i'
c      include 'mhdbcdtn.i'
      REAL(DP),ALLOCATABLE,DIMENSION(:) :: smthdata,dsmthdatadx,     
     .                                     d2smthdatadx2,weight,psi_temp
c
      vadjmin = -1.0e30
      vadjmax =  1.0e30
      klimold =  0
      indeq   =  0
      if (ii .ne. 0)  write (nitre, 8200)
 8200 format ('  entering plasma volume iteration')
      do it=1,itvol
        psibdry = psir(nj)
        psiaxis = psir(1)
        if (ii .eq. 0)  go to 2080
        iounit = 0                       ! turn off output from FREEBDRY
        if (ieqprt .eq. 1)  iounit = ncrt
        if (mhdmode .eq. 'coils') then
          if (minchisq .eq. 1 .and. ic .lt. 2)
     .    call fittocur (tocur, iounit)
          call freebdry (iteq, iounit, indeq, relerr)
        else
          call fixdbdry (iteq,toleq,iounit,indeq,relerr,
     .                                   tocur,ieq,omeq,ic)
        end if
c
c --- if we came close to the required error don't quit
c
        if (ABS (relerr-toleq) / toleq .lt. 5.0)  indeq = 1
        indeq = 1                         ! assume converged, always
        if (ieqfail .eq. 0)  go to 2060
        write (ncrt , 8100)  ieqfail
        write (nout , 8100)  ieqfail
        write (nitre, 8100)  ieqfail
 8100   format (' no convergence in freebdry - calculation aborted' /
     .          ' ieqfail in subroutine DRIVE =', i5)
        return
c
 2060   if (indeq .eq. 1)  go to 2080
        ieqfail = 1
        write  (ncrt , 8120)
        write  (nout , 8120)
        write  (nitre, 8120)
 8120   format (' freebdry did not converge; time step will be halved')
        return
c
c ----------------------------------------------------------------------
c psiset constructs a psi-grid ,psival,on the MHD (i.e., npsi) grid
c psival(1) = plasma edge,psival(npsi) = mag axis:
c if iunfrm = 0 spaced in psi square
c          ne 0 spaced in psi
c ----------------------------------------------------------------------
c
 2080   iunfrm = 1
        if(mhdmethd.eq.'toq') goto 2081
        call psiset (npsi, psiaxis, psibdry, psival, iunfrm)
        
c
c ----------------------------------------------------------------------
c --- given ffprim on psir grid integrate ffprim to get f(psi)
c --- on psival grid. ffprim is not changed
c --- HOWEVER: A matching ffprim is returned as ffppsival(psi) - D.F.
c ----------------------------------------------------------------------
c
****    if (fpsi(npsi) .eq. 0.0)  call fofpsi (flim)
        call fofpsi (flim)
c
c ----------------------------------------------------------------------
c FLUXAV12 returns a set of quantities defined on the psival grid;
c     flux is in units of volt-sec; dimensions are in m;
c     bax is in tesla; and curaxx is in amp/m**2
c ----------------------------------------------------------------------
c
        xmagn1 = rma
        ymagn1 = zma
        curaxx = -xax(1)*pprim(1)-ffprim(1)/(xax(1)*u0)
****    if (ii .ne. 0)  curaxx = 1.0e7*curden(1)*xax(1)/rmajor
c        bax    = fpsi(1)/xax(1)    ! tesla
        bax = fpsi(npsi)/xax(1)       ! changed 8/24/06 HSJ
        
c
c     set iscr = 1 for FLUXAV12 write to file 'scratch1'
c
        iscr  = 1
        ixcal = 1 
        ispln = 0
!        nscr = 0
!        If(save_scratch1 > 0)call getioun(nscr,nscr) 
        call SECOND (timeflx)
        write (ncrt, '(a,i5)')  ' calling FLUXAV12,npsi =',npsi

        call fluxav12 (psi,rmhdgrid,nw,zmhdgrid,nh,psival,npsi,fpsi,
     .               ffppsival,qpsival,bax,
     .               curaxx,rma,zma,ixcal,btor,rmajor,nscr,iscr, ncrt,
     .               rplasmin,rplasmax,zplasmin,zplasmax,elongx,epsp,
     .               vprime,ratave,psivolp,xhm2p,xi11p,xi33p,
     .               xipsp,ratavei,ravg,ravgi,bsqinvavg,bsq_avg,b_avg,
     .               h_factr,ierr, widep,hitep,kappa,
     .               rminor,rgeom,btgeom,circum,rvloop,zvloop,psivloop,
     .               tocur,ali,mhdmode,mhdmethd,psifctr,cxareanpsi,
     .               grho1npsi,grho2npsi,rmajavnpsi,rminavnpsi,
     .               triangnpsi,triangnpsi_l,pindentnpsi,sfareanpsi,
     .               torfluxnpsi,ifixshap)
        write (ncrt, '(a)')  ' returned from FLUXAV12'
        if(smth_mhd_parms .gt. 0.0)then
c           smooth the triangularity, triangnpsi, and 
c           elongation,elongx. These quantities are interpoalted 
c           onto the rho grid and then output in the iterdb file.
            no_weight = 1
            IF(.not. allocated(smthdata))ALLOCATE(smthdata(npsi),             
     .                         dsmthdatadx(npsi),d2smthdatadx2(npsi),
     .                         weight(npsi),psi_temp(npsi))

            smthdata(1:npsi) = triangnpsi(1:npsi)
            psi_temp(1:npsi) = ABS(psival(1:npsi))
            derivr = (triangnpsi(npsi) -triangnpsi(npsi-1))/
     .               (psi_temp(npsi) -psi_temp(npsi-1)) 
            derivl = 0.0
            derivr = (triangnpsi(2) -triangnpsi(1))/
     .               (psi_temp(2) -psi_temp(1)) 
            derivr = 0.0
            CALL smooth_values(smthdata,psi_temp,npsi,triangnpsi,
     .                       dsmthdatadx,d2smthdatadx2,weight,
     .                       no_weight,smth_mhd_parms,1,derivl,derivr)

            WRITE(*,FMT='("Spline smooth output upper")')
c            do j = 1, npsi
c             WRITE(*,FMT='(i3,4(1pe12.4,2x))')j,psival(j),
c     .             smthdata(j),triangnpsi(j),weight(j)
c            enddo

            smthdata(1:npsi) = triangnpsi_l(1:npsi)
            psi_temp(1:npsi) = ABS(psival(1:npsi))
            derivr = (triangnpsi_l(npsi) -triangnpsi_l(npsi-1))/
     .               (psi_temp(npsi) -psi_temp(npsi-1)) 
            derivl = 0.0
            CALL smooth_values(smthdata,psi_temp,npsi-1,triangnpsi_l,
     .                       dsmthdatadx,d2smthdatadx2,weight,
     .                       no_weight,smth_mhd_parms,2,derivl,derivr)

c            WRITE(*,FMT='("Spline smooth output lower")')
c            do j = 1, npsi
c             WRITE(*,FMT='(i3,4(1pe12.4,2x))')j,psival(j),
c     .             smthdata(j),triangnpsi_l(j),weight(j)
c            enddo

            smthdata(1:npsi) = elongx(1:npsi)
            CALL smooth_values(smthdata,psi_temp,npsi-1,elongx,
     .                       dsmthdatadx,d2smthdatadx2,weight,
     .                       no_weight,smth_mhd_parms,2,derivl,derivr)

c            WRITE(*,FMT='("Spline smooth output elongation")')
c            do j = 1, npsi
c             WRITE(*,FMT='(i3,4(1pe12.4,2x))')j,psival(j),
c     .                            smthdata(j),elongx(j),weight(j)
c            enddo
        endif

c
c ----------------------------------------------------------------------
c define the toroidal box used in subroutine TIMTOR.
c note that rmin is not used there. instead it is rin, the inside
c radius of the vessel that is used. This distinction is important
c for re-entrant neutrals. hence do not change rin to rmin !!!!!!!!
c rmin,rmax,zmin,zmax are in cm and are stored in file machin.i
c ----------------------------------------------------------------------
c
        rmin    = 100.0 * rplasmin
        rmax    = 100.0 * rplasmax
        zmin    = 100.0 * zplasmin
        zmax    = 100.0 * zplasmax
        psibdry = psival(1)
        call SECOND (timeflx1)
        timeflx = timeflx1 - timeflx
        write  (nitre, 2050)  timeflx
 2050   format (' time in FLUXAV12 = ', 1pe12.6, ' seconds')
 2081   continue
c
        if (raneut .eq. 0.0)
     .    raneut = SQRT (kappa) * rminor * 100.0  ! in cm
        do i=1,npsi
          ratave(i) = rmajor**2*ratave(i)    ! <R0**2/R**2>
        end do
        volume = psivolp(1)
        if (    ii .eq. 0    )  go to 5000
        if (tolvol .ge. 100.0)  go to 5000
        volerr = ABS (volume-volwant) / volwant
        del    = 1.0e30
        do k=1,3,2
          delx = ABS (xlim-xlimpos(k))
          if (delx .le. del) then
            del  = delx
            klim = k
          end if
        end do
c
        write  (ncrt , 8000)  volume, volwant, voladj, klim
        write  (nitre, 8000)  volume, volwant, voladj, klim
        write  (nout , 8000)  volume, volwant, voladj, klim
 8000   format (/ ' volume,volwant,voladj,klim ', 3e15.3, 2x, i5)
        if (limpos .eq. 'left'  .and. klim .ne. 1)  go to 3300
        if (limpos .eq. 'right' .and. klim .ne. 3)  go to 3300
        if (volerr .lt. tolvol)  go to 5000
c
 3300   deladj  = ABS (vadjmax-vadjmin)
        if (deladj .gt. dvadjmax)  go to 3360
        volwant = 0.98 * volume
        do i=1,mxtbcmhd
          if (volaray(i) .gt. volwant) volaray(i) = volwant
        end do
c
 3360   call getadjv (voladj,volume,volwant,volnudge,klim,
     .                klimold,limpos,vadjmin,vadjmax)
      end do
c
      ieqfail = 1
      write  (ncrt , 8130)
      write  (nout , 8130)
      write  (nitre, 8130)
 8130 format (' ERROR: could not find proper volume - run aborted')
      return
c
 5000 if (ii .ne. 0)  write (nitre, 8210)
 8210 format (' done with plasma volume iteration')
c
      write  (ncrt , 8150)  widep, hitep, kappa, volume, xax(1), yax(1),
     .                      elongx(npsi)
      write  (nitre, 8150)  widep, hitep, kappa, volume, xax(1), yax(1),
     .                      elongx(npsi)
      write  (nout , 8150)  widep, hitep, kappa, volume, xax(1), yax(1),
     .                      elongx(npsi)
      write  (nqik , 8150)  widep, hitep, kappa, volume, xax(1), yax(1),
     .                      elongx(npsi)
 8150 format (// ' plasma shape parameters - MKS units'               //
     .           ' width = ', f8.1, '  height =', f8.1, '  h/w =', f8.2,
     .        '  volume = ', e12.3, '  magnetic axis ', f8.2, 1x, f8.2 /
     .         ' kappa0 = ', f10.4 //)


      return
c
      end

      subroutine eqdskext (iounit)
c
      USE param
      USE ions
      USE soln
      USE contour
      USE limiter
      USE mhdpar
      USE nub
      USE nub2 
      USE numbrs   
      USE mesh
      USE rhog
      USE mhdcom
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c --- This subroutine extends the eqdsk to include plasma data
c ----------------------------------------------------------------------
c

c
      parameter (nfast = ke * kb + 1)
      dimension  efion(nfast), bmpwr(nfast), azfst(nfast), atmfst(nfast)
      dimension  dum(kj)
c$$$$    TKM   10/99, 04/00, 05/00
      dimension  denbm(kj), tembm(kj)
c$$$$$$$$$$$$$$$$$
c
      ncontr = 0
      write  (iounit, 8020)  ncontr, nlimiter
 8020 format (5i5)
      write  (iounit, 8010) (rcontr  (i),   zcontr(i), i=1,ncontr  )
      write  (iounit, 8010) (xlimiter(i), ylimiter(i), i=1,nlimiter)
 8010 format (5e16.9)
c
c ----------------------------------------------------------------------
c --- the following output is not part of the "official" eqdsk. It is
c --- what distinguishes eqdsk created in ONETWO from those created by EFIT.
c --- we put out enough information so that the eqdsk could be used as
c --- a restart file
c ----------------------------------------------------------------------
c  Latest update:  TKM    05/00
c
      nti = 1
      write (iounit, 8020)  nj, nprim, nimp, nti
      write (iounit, 8010) (atw(k),k=1,nion)
      do 700 k=1,nion
  700 write (iounit, 8010) (z(j,k),j=1,nj)
      write (iounit, 8010) (zeff(j),j=1,nj)
      write (iounit, 8010) (r(j),j=1,nj)
      write (iounit, 8010) (psir(j),j=1,nj)
      write (iounit, 8010) (te(j),j=1,nj)
      do 710 k=1,nti
  710 write (iounit, 8010) (ti(j),j=1,nj)
      do 720 j=1,nj
  720 dum(j) = 1.0e6 * ene(j)
      write (iounit, 8010)  (dum(j),j=1,nj)
      do 730 k=1,nprim
      do 731 j=1,nj
  731 dum(j) = 1.0e6 * en(j,k)
  730 write (iounit, 8010) (dum(j),j=1,nj)
c
      if (nimp .ne. 0) then
        do 740 k=nprim+1,nion
          do 741 j=1,nj
  741       dum(j) = 1.0e6 * en(j,k)
  740     write (iounit, 8010)  (dum(j), j=1,nj)
      end if
c
      if (ibeam.lt.3) then
         nalfa = 0
         write (iounit, '(2i5)') nalfa,nbeams
      else
          nalfa = ke*nbeams
          do ib=1,nbeams
             do ie=1,ke
                j        = (ib-1) * ke + ie
                efion(j) = ebkev(ib) / (2.0**(ie-1))
                bmpwr(j) = bptor(ib) * fbcur(ie,ib)
             end do
          end do
c
          do ib = 1,nbeams
             efion(ib) = ebkev(ib)
             bmpwr(ib) = 0.
             do ie = 1,ke
               bmpwr(ib) = bmpwr(ib)+bptor(ib)*fbcur(ie,ib)
             end do
          enddo
          do nab=1,2
            if (nameb .eq. namep(nab)) then
              do ib=1,nbeams
                do ie=1,ke
                  j         = (ib-1) * ke + ie
                  azfst (j) = atomno(nab)
                  atmfst(j) = atw(nab)
                end do
c               azfst(ib) = atomno(nab)
c               atmfst(ib) = atw(nab)
               end do
            end if
          end do
          write (iounit, '(2i5)') nalfa,nbeams
          do ib=1,nbeams
            do ie=1,ke
              i = (ib-1) * ke + ie
              write (iounit, '(4e16.9)')
     .                     efion(i), bmpwr(i), azfst(i), atmfst(i)
              write (iounit, 8010) (hdep(j,ie,ib),j=1,kj)
              write (iounit, 8010) (zeta(j,ie,ib),j=1,kj)
            end do
          end do
          do j = 1,kj
             denbm(j) = 0.
             tembm(j) = 0.
             do ib = 1,nbeams
                do ie = 1,ke
                  denbm(j) = denbm(j)+enb(j,ie,ib)
                   tembm(j) = tembm(j)+6.24142e15*wb(j,ie,ib)
                enddo
             enddo
             if(denbm(j).gt.0.) then
                tembm(j) = 0.667*tembm(j)/denbm(j)
             else
                tembm(j) = 1.0
             endif
          enddo
          write (iounit, 8010) (denbm(j),j=1,kj)
          write (iounit, 8010) (tembm(j),j=1,kj)
      endif
c $$$$$$$$$$$$$$$$$
      return
c
      end

      subroutine eqplot (ii)

c ----------------------------------------------------------------------
c this subroutine outputs various MHD quantities to the file eqpltfil
c this file is later accessed by the standalone plot program EQPLOT
c ----------------------------------------------------------------------
c
c
      USE param 
      USE io
      USE solcon
      USE soln
      USE contour
      USE mhdpar    
      USE mhdgrid 
      USE extra 
      USE numbrs
      USE mesh
      USE sourc
      USE machin
      USE geom
      USE tordlrot
      USE constnts
      USE soln2d
      USE psig
      USE rhog
      USE mhdcom
      USE flx,              ONLY : flux,fluxe,fluxi,qieneo
      USE bicube
      USE shapctr
      USE flxav
      USE etc
      USE tmpcom
      USE gpsi
      USE mhdbcdtn
      USE replace_imsl,      ONLY : my_ibcccu

      implicit  integer (i-n), real*8 (a-h, o-z)
c

      include 'imsl.i'


      include 'storage.i'

c
      parameter   (nrad = 2 * nw)
      dimension    rad(nrad), terad(nrad), tirad(nrad), enerad(nrad),
     .             safrad(nrad), rhorad(nrad), psirad(nrad),
     .             currad(nrad), current(nw,nh)
      imslmd = 'eqpl00'
      if (ii .ne. 0) then
      do 110 j=1,nh
        do 110 i=1,nw
          k0 = (i-1) * nh + j
  110     current(i,j) = pcurrent(k0)
      else
        call zeroa (current, nwh)
        do 115 j=1,nh
          if (zmhdgrid(j) .lt. zplasmin)  go to 115
          if (zmhdgrid(j) .gt. zplasmax)  go to 115
          do 116 i=1,nw
            if (rmhdgrid(i) .lt. rplasmin)  go to 116
            if (rmhdgrid(i) .gt. rplasmax)  go to 116
            if (   psi(i,j) .gt. psibdry )  go to 116
            call intrp (1,1,psir, pprim,nj,psi(i,j), pprim1,1)
            call intrp (1,1,psir,ffprim,nj,psi(i,j),ffprim1,1)
            current(i,j) = -rmhdgrid(i)*pprim1-ffprim1/(u0*rmhdgrid(i))
  116     continue
  115   continue
      end if
c
      if (ieqfail .gt. 0)  go to 2040
c
c --- set up bicubic spline representation of psi
c
      imslmd='4654c207'
      call my_ibcccu (psi,rmhdgrid,nw,zmhdgrid,nh,cspln,nw,wnoperm,ier)
c
c --- get major radius vector and corresponding rho values
c
      isigncur = +1
      if (psiaxis .lt. psibdry)
     .isigncur = -1
      psiedge  = -1.0e30
      call getrmaj (zma, nrad, ierr, rad, psirad, isigncur, psiedge)
c
c --- rad(1...nrad) now contains major radius (in meters) and psirad(1..nrad)
c --- contains corresponding psi values (in volt-sec/rad)
c
      imslmd = 'eqpl01'
      call intrp (1,1,psival,rho,npsi,psirad,rhorad,nrad)
      imslmd = 'eqpl02'
      call intrp (0,1,r,te,nj,rhorad,terad,nrad)
      imslmd = 'eqpl03'
      call intrp (0,1,r,ti,nj,rhorad,tirad,nrad)
      imslmd = 'eqpl05'
      call intrp (0,1,r,ene,nj,rhorad,enerad,nrad)
      imslmd = 'eqpl05'
      call intrp (0,1,r,q,nj,rhorad,safrad,nrad)
      imslmd = 'eqpl05'
c
 2040 kfar = nk - iangrot
c
      write (neqplt,'(a)') ' **** continue ****'
      write (neqplt, 8020)  ieqfail,fixfcoil,ncontr,nxeqd,irguess,kfar
      write (neqplt, 8010)  time,psiaxis,psibdry,kappa,volume,tocur
      write (neqplt, 8010)  circum,betapmhd,betatmhd,pvbar,btor
      write (neqplt, 8010)  ((    psi(ix,iy), ix=1,nw), iy=1,nh)
      write (neqplt, 8010)  ((current(ix,iy), ix=1,nw), iy=1,nh)
      write (neqplt, 8010)  (psir(j),pprim(j),ffprim(j), j=1,nj)
      write (neqplt, 8010)  (   press(i),i=1,nj)
      write (neqplt, 8010)  (  pressb(i),i=1,nj)
      write (neqplt, 8010)  (prseqdsk(i),i=1,nxeqd)
      write (neqplt, 8010)  (psieqdsk(i),i=1,nxeqd)
      write (neqplt, 8010)  (  r(i),q(i),i=1,nj)
      write (neqplt, 8010)  (  curden(i),i=1,nj)
      write (neqplt, 8010)  (  curohm(i),i=1,nj)
      write (neqplt, 8010)  (  curdri(i),i=1,nj)
      write (neqplt, 8010)  ( curboot(i),i=1,nj)
      write (neqplt, 8010)  (    gcap(i),i=1,nj)
      write (neqplt, 8010)  (    hcap(i),i=1,nj)
      write (neqplt, 8010)  (    fcap(i),i=1,nj)
      write (neqplt, 8010)  (     rbp(i),i=1,nj)
c
c --- 2d source terms in Faraday's Law
c
      write (neqplt, 8010)  (  fday2d1(i), i=1,nj)
      write (neqplt, 8010)  (  fday2d2(i), i=1,nj)
      write (neqplt, 8010)  (  fday2d3(i), i=1,nj)
      write (neqplt, 8010)  (   dscrip(i), i=1,nj)
      write (neqplt, 8010)  (  scurdri(i), i=1,nj)
      write (neqplt, 8010)  ( scurfast(i), i=1,nj)
      write (neqplt, 8010)  (     etap(i), i=1,nj)
      write (neqplt, 8010)  (   s(kfar,i), i=1,nj)
      write (neqplt, 8010)  (flux(kfar,i), i=1,nj-1)
c
c --- plasma boundary at the current time
c
      write (neqplt, 8010)  (rcontr(i),zcontr(i), i=1,ncontr)
c
c --- MSE data if requested
c
      if (use_stark) then
        write (neqplt, 8010)  (  tstarkexptl(i), i=1,nstark)
        write (neqplt, 8010)  (sigstarkexptl(i), i=1,nstark)
        write (neqplt, 8010)  (  rstarkexptl(i), i=1,nstark)
        write (neqplt, 8010)  (  zstarkexptl(i), i=1,nstark)
        write (neqplt, 8010)  (   tstarkcalc(i), i=1,nstark)
        write (neqplt, 8010)  (   chisqstark(i), i=1,nstark)
        write (neqplt, 8010)      chisqstarktot
      end if
      if (mhdmode .eq. 'coils') then
c
c --- exptl psiloop and probe values
c
        write (neqplt, 8010)  ( psiloop(i), i=1,nsilop)
        write (neqplt, 8010)  (probeval(i), i=1,magpr2)
c
c --- calculated psiloop and probe values
c
        write (neqplt, 8010)  (psilopcl(i), i=1,nsilop)
        write (neqplt, 8010)  (probeclc(i), i=1,magpr2)
c
c --- reduced xchisq and number of observations
c
        write  (neqplt, 8015)  xchisq, nobsrvd
 8015   format (e13.6, 2x, i5)
      end if
c
c --- induced parallel electric field
c
        write (neqplt, 8010)  (etor(i), i=1,nj)
c
c --- loop voltage, inductance
c
        write (neqplt, 8010)  vlopmhdt, vloopcl(ieq+1), ali,
     .                        voltag, voltoh, psivloop
        if (fixfcoil .eq. 0 .and. mhdmode .eq. 'coils') then
c
c --- calculated fcoil currents
c
          write (neqplt, 8010)  (curfcoil(i), i=1,nfcoil)
c
c --- experimental fcoil currents
c
          write (neqplt, 8010)  (expfcoil(i), i=1,nfcoil)
        end if
****  end if             ! changed so that  'no coils'  option will work
c
      if (ieqfail .gt. 0)  return
c
      write (neqplt, 8020)                  nrad
      write (neqplt, 8010)  (   rad(j), j=1,nrad)
      write (neqplt, 8010)  (rhorad(j), j=1,nrad)
      write (neqplt, 8010)  ( terad(j), j=1,nrad)
      write (neqplt, 8010)  ( tirad(j), j=1,nrad)
      write (neqplt, 8010)  (enerad(j), j=1,nrad)
      write (neqplt, 8010)  (safrad(j), j=1,nrad)
c
      imslmd = 'eqpl07'
      do j=1,nrad
        call intrp (1,1,psir, pprim,nj,psirad(j), pprim1,1)
        call intrp (1,1,psir,ffprim,nj,psirad(j),ffprim1,1)
        currad(j) = -rad(j)*pprim1-ffprim1/(u0*rad(j))
      end do
c
      write  (neqplt, 8010)  (currad(j), j=1,nrad)
 8010 format (6e13.6)
 8020 format (6(2x, i5))
      imslmd = 'eqpl08'
      return
c
      end

      subroutine exptlclc (kstart, kend)
c
c
c ----------------------------------------------------------------------
c --- subroutine calculates psiloop and magnetic probe values,so that
c --- these can be compared with experimentally measured quantities.
c    kstart              defines support set of pcurrent,
c    kend                to avoid going over whole grid
c
c --- input(through INCLUDE files):
c --- INCLUDE file mhdpar
c   nfcoil           these are the parameters set in mhdpar
c   nsilop           which are used in this subroutine
c   magpr2
c   nwh
c   nesum
c   nvessel
c
c --- common block mhdcom (INCLUDE file mhdpar):
c   pcurrent(nwh)  toroidal plasma filament current (amps)
c   ifitpsi       =1 use psi loop values in fcoil cur determination
c                 =0 do not use psi loop values in fcoil cur determination
c   ifitprob      =1 use mag probe values in fcoil cur determination
c                 =0 do not use mag probe values in fcoil cur determination
c   iecurr        =1 include ecoil current effects
c                 =0 exclude ecoil current effects
c   ivessel       =1 include vessel current effects
c                 =0 exclude vessel current effects
c
c --- from common blocks associated with Green's table (INCLUDE file mhdcom):
c --- these are the various inductive coupling factors used in this subroutine
c   rsilfc(nsilop,nfcoil)         fcoil to psi loop
c   rmp2fc(magpr2,nfcoil)         fcoil to magnetic probe
c   rsilec(nsilop,nesum)          ecoil to psi loop
c   rmp2ec(magpr2,nesum)          ecoil to magnetic probe
c   rsilvs(nsilop,nvessel)        vessel to psi loop
c   rmp2vs(magpr2,nvessel)        vessel to magnetic probe
c   rsilpc(nsilop,nwh)            plasma to psi loop
c   rmp2pc(magpr2,nwh)            plasma to magnetic probe
c
c --- common block mhdbctdn(INCLUDE file mhdbctdn)
c   psiloop(nsilop)        psi loop values at the current time
c                          (not used if ifitpsi = 0)
c   probeval(magpr2)       probe values at the current time
c                          (not used if ifitprob = 0)
c   ecoilcur(nesum)        e coil current at the current time
c                          (not used if iecurr = 0)
c   vescur(nvessel)        vessel current at the current time
c                          (not used if ivessel = 0)
c   curfcoil(nfcoil)       calculated fcoil currents,amps
c  fwtpsilp(nsilop)
c  fwtmprbe(1..magpri)
c  fwtfcoil(1..nfcoil)
c  fwtecoil(1...nesum)
c  fwtvescr(1..nvessel)       the weighting vectors ( = 1/sigma)
c                             fwtecoil and fwtvescr are currently not used
c  isgngren                   sign to use (+/- 1) see getfcur
c
c --- output
c ------------------------------------------------------------------ HSJ
c
      USE mhdpar      
      USE mhdcom
      USE mhdbcdtn
      implicit  integer (i-n), real*8 (a-h, o-z)

c
        do 500 m=1,nsilop
          psilopcl(m) = 0.0
c
c --- get contribution to psi loop m from all sources (other than fcoil)
c --- a)plasma current
c
            do 100 k=kstart,kend
  100       psilopcl(m) = psilopcl(m)+isgngren*rsilpc(m,k)*pcurrent(k)
            if (iecurr .ne. 0.0) then
c
c --- b)from ecoils
c
              do 200 k=1,nesum
  200         psilopcl(m) = psilopcl(m)+isgngren*rsilec(m,k)*ecoilcur(k)
            end if
            if (ivessel .ne. 0) then
c
c --- c)from vessel currents
c
              do 300 k=1,nvessel
  300         psilopcl(m) = psilopcl(m)+isgngren*rsilvs(m,k)*vescur(k)
            end if
c
c --- d)from fcoils
c
          do 400 k=1,nfcoil
  400        psilopcl(m) = psilopcl(m)+isgngren*curfcoil(k)*rsilfc(m,k)
  500   continue
c
c --- next get the probe values
c
        do 1000 m=1,magpr2
          probeclc(m) = 0.0
c
c --- a)plasma current contribution
c
            do 900 k=kstart,kend
  900       probeclc(m) = probeclc(m)+rmp2pc(m,k)*pcurrent(k)
            if (iecurr .ne. 0) then
c
c --- b)ecoil contribution
c
              do 1200 k=1,nesum
 1200         probeclc(m) = probeclc(m)+rmp2ec(m,k)*ecoilcur(k)
            end if
            if (ivessel .ne. 0) then
c
c --- c)vessel contribution
c
              do 1300 k=1,nvessel
 1300         probeclc(m) = probeclc(m)+rmp2vs(m,k)*vescur(k)
            end if
c
c --- d)fcoil current contribution
c
          do 1500 k=1,nfcoil
 1500        probeclc(m) = probeclc(m)+rmp2fc(m,k)*curfcoil(k)
 1000   continue
c
c --- calculate reduced chi-square
c
      xchisq  = 0.0
      nobsrvd = 0
      if (ifitpsi .eq. 1) then
        do j=1,nsilop
          if (fwtpsilp(j) .ne. 0.0) then
            nobsrvd = nobsrvd + 1
            xchisq  = xchisq + (fwtpsilp(j)*(psiloop(j)-psilopcl(j)))**2
          end if
        end do
      end if
      if (ifitprob .eq. 1) then
        do j=1,magpr2
          if (fwtmprbe(j) .ne. 0.0) then
            nobsrvd = nobsrvd + 1
            xchisq  = xchisq+(fwtmprbe(j)*(probeval(j)-probeclc(j)))**2
          end if
        end do
      end if
      if (fixfcoil .eq. 0 .and. fitfcur .eq. 1) then
        do j=1,nfcoil
          if (fwtfcoil(j) .ne. 0.0) then
            nobsrvd = nobsrvd+1
            xchisq  = xchisq+(fwtfcoil(j)*(expfcoil(j)-curfcoil(j)))**2
          end if
        end do
      end if
      if (nobsrvd .ne. 0)  xchisq = xchisq / nobsrvd
      return
c
      end




      subroutine runtwo
c
      USE param
      USE io
      USE fusion
      USE ions
      USE solcon
      USE soln
      USE transp
      USE contour
      USE limiter
      USE mhdpar
      USE mhdgrid 
      USE nub 
      USE numbrs
      USE rf 
      USE extra
      USE mesh
      USE sourc
      USE machin
      USE tfact
      USE geom
      USE tordlrot
      USE constnts
      USE soln2d
      USE bd_condtn
      USE rhog
      USE mhdcom
      USE bicube
      USE shapctr
      USE flxav
      USE iterdbmd                     ! get irwflag for old iterdb from here

      USE iterdbmd_gcnmp,               ONLY : irwflag_gcnmp => irwflag
      USE etc 
c      USE  io,                         ONLY : ioftn, iopntr, nunits

      USE  set_12_gcnmp_vars,           ONLY : set_gcnmp_vars,     
     .                                          set_onetwo_vars

      USE iterdbmd_gcnmp,               ONLY : iterdb_file_name

      USE gcnmp_input,                  ONLY : write_iterdb_txt,
     .                                         switch_iterdb_output

      USE  kinetic_efit,                ONLY : kine_message


      USE mhdbcdtn
      USE pelcom
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      character rcs_id*63
      save      rcs_id
      data      rcs_id /
     ."$Id: cray207.f,v 1.85 2014/04/29 18:15:45 stjohn Exp $"/
c
c ----------------------------------------------------------------------
c  driver for 1-1/2-d runs (codeid .ne.  'onedee')
c
c  THE PROBLEM WITH UNITS IS HANDLED AS FOLLOWS:
c    on input to RUNTWO the units are gaussian, but
c    psi(i,j), p(i,j), psir(j), psival(i), psiaxis, psibdry, and
c    bp(j) carry units of kilogauss; temperatures are in keV
c
c    subroutine MHD converts all necessary quantitites to MKS,
c    does the required equilibrium calculations in MKS,
c    and then converts back to the above gaussian system.
c ----------------------------------------------------------------------
c


      include 'pckcom.i'

      include 'storage.i'

c
      external FLUSH
c     some temporary stuff to test xplasma:

c
c ----------------------------------------------------------------------
c start up MHD calculations. all MHD calcs are done in MKS units.
c obtain guess for equilibrium flux, pprime, and ffprime
c make initial call to equilibrium
c ----------------------------------------------------------------------
c
      ieq    = 0
      itrfix = itrfix0
      itre   = 1



      eq_init: if (irguess .eq. 0) then 
c
c --- convert curden to a/m**2
c
        do 100 j=1,nj
  100     curden(j) = curden(j)*1.0e4
        tocur = pcurmhdt     ! tocur is defined by rbp in subroutine MHD
c
c --- at present irguess = 0 will not work until tocur is defined!!!!!
c
        btor   = btormhdt
        rmajor = rmajor*0.01
        flim   = btor*rmajor
        call initmhd (rmhdgrid, nw, zmhdgrid, nh, zero, rmagax,
     .                zmagax, rscale, zscale, rmajor, u0, tocur, cspln,
     .                n2cspln, nw, nh2, pds, psi1d,
     .                wnoperm, nlimiter, xlimiter, ylimiter, rcontr,
     .                zcontr, nconmax, curden, nj, ncrt, psi, pprim,
     .                ffprim, zdum, psir, press)
        rma     = rmagax
        zma     = zmagax
        psiaxis = psir(1)
        psibdry = psir(nj)
      else if (mhdmethd .ne. 'tdem') then eq_init
        if(.NOT. initialize_from_statefile )call reqdsk
        call dump_values ('after reqdsk')

      else    eq_init                                ! this is a tdem run
c
c       construct all quantities that appear on eqdsk from input file:
c
        call set_cdf_init
        call dump_values ('after set_cdf_init') !if tdemvb = 0 does nothing

      end if eq_init
      if (ieqfail .gt. 0)  return    ! failure in INITMHD routine






      write  (nitre, 8000)
      write  (ncrt , 8000)
 8000 format (' making initial call to equilibrium')
      call SECOND (tmhdstart)
      write(6,FMT = '(" calling mhd with ieq = ",i5)')ieq

      kine_message = "before_init"
      call mhd   ! defines tocur, total current for the MHD calculations

      call SECOND (tmhdend)
      timemhd = tmhdend - tmhdstart
      write  (6, 42)  ' time in MHD-1 module =', timemhd, ' seconds'

   42 format (a, f15.1, a)
      write  (nitre, '(a)')  ' done with initial equilibrium'




c
c ----------------------------------------------------------------------
c --- return if calculations failed or this was a single equilibrium run
c     (mhdonly = 1)
c ----------------------------------------------------------------------
c
      if (ieqfail .gt. 0)  return
      if (mhdonly .eq. 1)  return


      call getioun(nmhddat,nmhddat)
      open (unit = nmhddat, file = 'test_trnspt_mhd.txt',
     .    status = 'UNKNOWN')
c
c ----------------------------------------------------------------------
c --- set up the inititial profiles in vector u.
c --- subroutine REDATE copies en, te, ti, etc into vector u
c --- This is where the statefile results get into vector u 
c ----------------------------------------------------------------------
c

      call redate (u, en, te, ti, rbp, nk, nj, kj, kk, iangrot, angrot)
c


 2000 ieq = ieq + 1                 !loop over equilibrium time steps





c
c ----------------------------------------------------------------------
c --- write transport solution to binary file "savsol", to be read back
c --- if necessary by subroutine READIT, called from subroutine PREPAR
c ----------------------------------------------------------------------
c
      call savit
      if (ieq .gt. keqmax) then
        write  (nout, 1500)  ieq, keqmax
 1500   format ('  ieq =', i5 /
     .          '  the code must stop because a maximum of  ', i5 /
     .          '  equilibrium steps are allowed')
        return
      end if


c
c ----------------------------------------------------------------------
c start a new equilibrium/transport cycle
c if ibeam = -1 beam was previously on but should now be turned off
c             (accounts for residual beam slowing down,etc).
c if ibeam = 0 beam is currently off (and has never been on)
c if ibeam = 1 turn beam on
c if ibeam = 2 turn beam off
c if ibeam = 3 beam is currently on
c if irf = 1 turn RF on
c if irf = 2 turn RF off
c ----------------------------------------------------------------------
c
      itrfix = itrfix0
      itre   = 0
      if ( ibeam .eq. 1)  ibeam  =  3
      if ( ibeam .eq. 2)  ibeam  = -1
      do 2005 k=1,krf
      if (irf(k) .eq. 1)  irf(k) =  3
 2005 if (irf(k) .eq. 2)  irf(k) = -1
c
c ----------------------------------------------------------------------
c set up dteq, the time increment for the equilibrium calculation
c i.e., the next equilibrium calculation will be done
c dteq seconds after the current time. (unless dteq is cut down by the
c code). dteq is picked so that we land on a time given in deqlst
c ----------------------------------------------------------------------
c
      if (ifixshap .eq. 0 .and. mhdmethd .ne. 'tdem') then
          do j=2,ieqmax
            itindex = j
            if (deqlst(itindex) .gt. time)  go to 110
          end do
  110     if (deqlst(itindex)-time .lt. 1.0e-7)  itindex = itindex + 1
          itindex = MIN0 (itindex, ieqmax)
          dteq    = deqlst(itindex)-time
          if ((time+dteq) .gt. timmax)  dteq = timmax - time
      else                                                ! ifixshap = 1
          dteq    = timmax - time
          if (dteq .lt. 0.0)  dteq = 0.0
      end if
      if (dteq .gt. 0.0 .and. mhdmethd .eq. 'tdem')
     .    dteq = timmax - time + 100.0
c
c ----------------------------------------------------------------------
c start transport/equilibrium iterations
c ----------------------------------------------------------------------
c
 2010 itre  = itre + 1
      ptime = time + dteq
      if (time+dteq .gt. timmax)  ptime = timmax
      write  (nitre, 8010)  ieq, itre, time, ptime
      write  (ncrt , 8010)  ieq, itre, time, ptime
 8010 format (' current equilibrium iteration number =', i5    /
     .        ' current   transport iteration number =', i5    /
     .        ' using the current equilibrium, try to advance' /
     .        ' the transport calculations from time =', f14.6 /
     .                                 30x, 'to time =', f14.6 )
c
c --- flush standard output (unit 6 = ncrt)
c    (non-zero iflush implies that an error occurred in subroutine FLUSH)
c
      call FLUSH (ncrt, iflush)
c
      call SECOND (trsptstart)
c

      call tport                ! do some transport

      rbp_save(:) = rbp(:)

c
      call SECOND (trsptend  )
      trsptime = trsptend - trsptstart
      write (6, 42)  ' time in transport module =', trsptime, ' seconds'
c
c     calculate the psi grid corresponding to the transport rho grid:
c
      call psirho(0)
c
c   istop was set to 1 because of one of the following:
c     n > nmax
c     dt < dtmin
c     timmax < time0
c     nmax < 0
c     time > timmax
c   reset istop to 0 for the last case, so that the final equilibrium
c   will be converged if it isn't
c
      if (istop .eq. 1 .and. time .ge. timmax)  istop = 0
c
c ----------------------------------------------------------------------
c call MHD even if ifixshap or istop = 1, to obtain final o/p to file 'eqdsk'
c this requires that ieq is set to zero if dtime = 0 (as will be the case
c if timmax .le. time0) otherwise subroutine PREPAR will give trouble
c ----------------------------------------------------------------------
c
      if (timmax .le. time0)  ieq = 0
      write  (ncrt , 8020)  ieq, itre, time
      write  (nitre, 8020)  ieq, itre, time
 8020 format (' ieq, itre, time ', i4, i6, 2x, e12.6, ' calling MHD')


      call SECOND (tmhdstart)
      ! kine_messagee ='final' set in sub tport
      if (mhdmethd .ne. 'tdem')  call mhd
      call SECOND (tmhdend  )
      timemhd = tmhdend - tmhdstart
      write (6, 42)  ' time in MHD module =', timemhd, ' seconds'




c
c ------------------------------------------------------------------------------
c write text or netcdf version of iterdb (onetwo state)  file:
c ------------------------------------------------------------------------------
c

      CALL Statefile_proc
 
        go to 1509

c      above call takes the place of the following:



      if (ifixshap .eq. 1 .and. iterdb .eq.  1 .and.               
     .                                create_XPTOR_input )THEN
                      irwflag =0  ! write Xptor  text iterdb file
                      call iter_dbase129
      ENDIF
      irwflag_gcnmp = 0   ! write new statefile
      if (ifixshap .eq. 1 .and. iterdb .eq.  2 .and.               
     .              .NOT. create_GCNMP_input )  call iter_dbase !returns withot writing
      IF (ifixshap .eq. 1 .and. iterdb .eq.  1 .and.              
     .              create_GCNMP_input )THEN
          CALL set_gcnmp_vars
          IF(statefile_type == 0 .AND. switch_iterdb_output ==1)THEN
c          IF(write_iterdb_txt .AND. switch_iterdb_output ==1)THEN
           !here the gcnmp input file was in text mode and the gcnmp
           !output file is in netcdf form
              iterdb_file_name = 'iterdb.nc'
              CALL iter_dbase_nc
          ELSEIF(statefile_type == 0 .AND. switch_iterdb_output ==0)THEN
           !here the gcnmp input file was in text mode and the gcnmp
           !output file is in text  form
              iterdb_file_name = 'iterdb.txt'
              CALL iter_dbase_txt
           ELSEIF(statefile_type == 1 
     .                             .AND. switch_iterdb_output ==1)THEN
c          ELSEIF(.NOT. write_iterdb_txt .AND. 
c     .                                  switch_iterdb_output ==1)THEN
           !here the gcnmp input file was in netcdf mode and the gcnmp
           !output file is in text form
              iterdb_file_name = 'iterdb.txt'
              CALL iter_dbase_txt
          ELSEIF(statefile_type == 1 .AND. switch_iterdb_output ==0)THEN
            ! here the input and output files are in netcdf form
             iterdb_file_name = 'iterdb.nc'
              CALL iter_dbase_nc
          
          ELSEIF(statefile_type == -1)THEN
             ! here there was no input statefile. write the output according to
             ! setting of write_iterdb_txt
             IF(write_iterdb_txt )THEN
               iterdb_file_name = 'iterdb.txt'
               CALL iter_dbase_txt
             ELSE
               iterdb_file_name = 'iterdb.nc'
               CALL iter_dbase_nc
             ENDIF
          ENDIF
      ENDIF
c      Highjacking the iterdb = -1 netcdf write to write the onetwo statefile
c      instead of the old netcdf iterdb file:
c      if (ifixshap .eq. 1 .and. iterdb .eq. -1)  call iter_dbase1
c
c      if (mhdmethd .eq. 'tdem' .and. iterdb .eq. -1 
c     .           .AND. .NOT. create_GCNMP_input) call iter_dbase1



      if (mhdmethd .eq. 'tdem' .and. iterdb .eq.  2
     .          .AND. .NOT. create_GCNMP_input) call iter_dbase !returns withot writing



      IF ( (iterdb == -1 
     .          .AND.  mhdmethd == 'tdem')  .OR.  (iterdb == -1 
     .          .AND.  ifixshap == 1))THEN
          write_iterdb_txt = .FALSE.
          iterdb_file_name = 'iterdb.nc'
          CALL set_gcnmp_vars
          CALL iter_dbase_nc
      ENDIF

      IF (mhdmethd  ==  'tdem' .AND.iterdb  ==   1
     .           .AND.  create_GCNMP_input  )THEN
          CALL set_gcnmp_vars
          IF(write_iterdb_txt)THEN
             iterdb_file_name = 'iterdb_tdem.txt'
             CALL iter_dbase_txt
          ELSE
             iterdb_file_name = 'iterdb_tdem.nc'
             CALL iter_dbase_nc
          ENDIF
 
      ENDIF

 1509 Continue
 

****  if (ifixshap .eq. 1 .or.   istop .eq. 1)  return
      if (ifixshap .eq. 1                    )  return
      if (mhdmethd .eq. 'tdem')  return
c
c --- igoitr was set to 1 in subroutine PREPAR if rhoa and/or
c --- cap parameters are not converged.
c --- if converged then igoitr = 0
c --- if the equilibrium failed cut the time step in half immediately
c --- by setting igoitr = 1
c
      if (ieqfail .gt. 0) then
        itre   = maxitr
        igoitr = 1
      end if
      write  (nitre, 8030)  igoitr, ieqfail
 8030 format (' done with mhd, igoitr = ', i5, '  (0 means converged)' /
     .        ' ieqfail in subroutine RUNTWO =', i5)
      if (igoitr .ne. 0)  go to 2015
c
c ----------------------------------------------------------------------
c equilibrium/transport cycle has converged
c ----------------------------------------------------------------------
c
      write  (nout , 8015) time
      write  (ncrt , 8015) time
      write  (nitre, 8015) time
 8015 format (' transport/equilibrium cycle converged, time = ', f14.6)


c--------------------------------------------------------------------
c      here ifixshap ne 1 ( returns above otherwise)

c      xptor type text file:
      if (ifixshap .ne. 1 .and. iterdb .eq.  1 
     .          .AND. .NOT. create_GCNMP_input)  call iter_dbase129

c     iter_dbase(cray204.f) text file modified xptor not used (obsolete)
      if (ifixshap .ne. 1 .and. iterdb .eq.  2 
     .          .AND. .NOT. create_GCNMP_input)  call iter_dbase

c     text version of onetwo/gcnmp statefile
      if (ifixshap .ne. 1 .and. iterdb .eq.  1
     .  .AND.  create_GCNMP_input .AND. write_iterdb_txt) 
     .                                        call iter_dbase_txt

c     netcdf version of onetwo/gcnmp statefile:
      if (ifixshap .ne. 1 .and. iterdb .eq.  1
     . .AND.  create_GCNMP_input .AND. .NOT. write_iterdb_txt) 
     .                                        call iter_dbase_nc
c     new netcdf version with backwards compatability
      if (ifixshap .ne. 1 .and. iterdb .eq.  1
     . .AND.  .NOT. create_GCNMP_input .AND. .not. create_xptor_input 
     . .AND. .NOT. write_iterdb_txt ) call iter_dbase_nc_compat

c--------------------------------------------------------------------




c--------------------------------------------------------------------
c     these are not active  because ifixshap == 1 returns above

c     xptor type netcdf file:
      if (ifixshap .eq. 1 .and. iterdb .eq. -1
     .          .AND. .NOT. create_GCNMP_input)  call iter_dbase1

c     text version of onetwo/gcnmp statefile
      if (ifixshap .eq. 1 .and. iterdb .eq.  1
     .          .AND.  create_GCNMP_input .AND. write_iterdb_txt) 
     .                                        call iter_dbase_txt
c     netcdf version of onetwo/gcnmp statefile:
      if (ifixshap .eq. 1 .and. iterdb .eq.  1
     . .AND.  create_GCNMP_input .AND. .NOT. write_iterdb_txt) 
     .                                        call iter_dbase_nc
c-----------------------------------------------------------------





      if ( time .ge. timmax)  return
      if (istop .eq. 1)  return ! return even if time < timmax because..
c                               ..istop was set to 1 for some other reason
      go to 2000
c
c ----------------------------------------------------------------------
c --- no convergence:
c --- here we enter a type of predictor/corrector loop . The geometric
c --- parameters (rhoa and the cap quantities) calculated using the transport
c --- solution at t = time+dteq are not sufficiently close to those predicted
c --- to exist (in subroutine PREPAR) at this time.
c --- Consequently we redo the transport calculation,
c --- again advancing the transport dependent variables from t = time
c --- to t = time+dteq,but with a new set of time derivatives for the
c --- geometric parameters (see subroutine PREPAR).
c --- The time derivatives of the geometric parameters,
c --- which are assumed constant in the interval dteq, are used in tport
c --- to generate the proper time-dependent parameters during the dteq interval.
c --- (see subroutine RHOMSH for this.
c --- also, the time derivatives enter into the 2d source terms;
c --- see subroutine SOURCE).
c ----------------------------------------------------------------------
c
 2015 write (nitre, 8100)  itre, maxitr, igoitr
 8100 format ('  itre = ',i5, '  maxitr = ',i5, '  igoitr = ',i5)
      if (itre .lt. maxitr)  write (ncrt , 8105)
      if (itre .lt. maxitr)  write (nitre, 8105)
      if (itre .lt. maxitr)  write (nout , 8105)
 8105 format ('  rhoa and/or cap parameters not converged,' /
     .        '  transport/equilibrium cycle will have to be repeated')
      if (itre .lt. maxitr)  go to 2010
c
c ----------------------------------------------------------------------
c we have exceeded the maximum number of iterations
c (i.e. predictor/corrector steps mentioned above), so we
c will cut the time between equilibrium calculations in half.
c this should cause the change in the geometric parameters to
c become small enough for convergence to be achieved.
c ----------------------------------------------------------------------
c
 2020 dteq = dteq * 0.5
****  ieqmax = ieqmax+1
****  if (ieqmax .gt. keqmax)  go to 3000
****  do 2030 i=ieqmax,ieq+1,-1
****  deqlst(i) = deqlst(i-1)
*2030 continue
****  deqlst(ieq) = dteq
c
      write (ncrt,8050) ieq,itre,dteq
 8050 format (2x / 'ieq = ',i6, '  itre=',i6, '  dteq being halved to',
     .    e15.3 / '*********#########*********#########*******######')
      write (nitre,8050) ieq,itre,dteq
      write (nqik,8050) ieq,itre,dteq
c
c ----------------------------------------------------------------------
c reset iteration counter, itre
c ----------------------------------------------------------------------
c
      itre = 0
      if (dteq .ge. dtmine)  go to 2010
c
c ----------------------------------------------------------------------
c time between equilibrium calculations has been reduced below minimum;
c run will be terminated.
c ----------------------------------------------------------------------
c
 3000 write (ncrt, 8060)
      write (nout, 8060)
      write (nqik, 8060)
 8060 format (' exceeded keqmax or dteq reduced below dtmine ',
     .         '**** run will be aborted *****')
      return
c
      end
