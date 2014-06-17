      subroutine convert_to_fixed_grid (zeta, rho, v_dep, v_dep_out,
     .                                  zeta_loc, nj)
c
c
c ----------------------------------------------------------------------
c     subroutine converts v_dep to the fixed (in time) normalized
c     zeta grid
c
c     INPUT:
c     zeta(j)          j=1,2..nj   the fixed zeta grid
c     rho(j)           the rho grid
c     v_dep(j)         the dependent variable
c     zeta_loc(j)      work storage of length nj
c     nj
c
c     OUTPUT:
c       v_dep_out(j)   j=1,2..nj the values of v_dep interpolated onto
c                      the zeta grid.
c
c ----------------------------------------------------------------------
c

      USE param
      USE mhdpar
      USE bicube
      USE replace_imsl,                 ONLY : my_icsccu,my_icsevu
      implicit  integer (i-n), real*8 (a-h, o-z)

      dimension zeta(*),rho(*),v_dep(*),v_dep_out(*),
     .          zeta_loc(*),csp(kj,3)
c      equivalence (csp(1,1),cspln(1,1,1))
c
      rhomax = rho(nj)
      do j=1,nj-1
        zeta_loc(j)=rho(j)/rhomax
      end do
      zeta_loc(1) = 1.0          ! avoids roundoff problem
      call my_icsccu (zeta_loc,v_dep,nj,cspln,kj,ier)
      call my_icsevu (zeta_loc,v_dep,nj,cspln,kj,zeta,v_dep_out,nj,ier)
      return
c
      end

      subroutine dump_values (input_string)
c

      USE param
      USE io
      USE solcon
      USE soln
      USE contour
      USE limiter
      USE mhdpar
      USE mhdgrid
      USE tdem
      USE extra
      USE numbrs
      USE mesh
      USE geom
      USE flags
      USE tordlrot
      USE constnts
      USE soln2d
      USE bd_condtn,only : bctime,ub,fluxb,
     .    ub_save,ub_rho_edge,bctime_zone
      USE psig
      USE rhog
      USE mhdcom
      USE flxav
      USE neo2d
      USE etc
      USE gpsi
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c --- primarily used for debugging
c ------------------------------------------------------------------ HSJ
c
      logical   exists, opened
c
c      include 'param.i'
c      include 'bcon.i'
c      include 'contour.i'         ! get nconmax for tdem.i
c      include 'constnts.i'        ! psikgaus
c      include 'etc.i'             ! bp(i)
c      include 'extra.i'           ! toteqd,eloc,fact,q
c      include 'flags.i'           ! itran
c      include 'geom.i'            ! r2cap,rhoa,dfdt,dgdt,dhdt
c      include 'io.i'              ! eqdskin,ncrt
c      include 'limiter.i'
c      include 'neo2d.i'           ! eps,xips,xi11,xi33,xhm2
c      include 'numbrs.i'          ! nj,nk
c      include 'mesh.i'            ! r(i)
c      include 'mhdpar.i'          ! kpsi
c      include 'mhdgrid.i'
c      include 'mhdcom.i'          ! psi(i,j)
c      include 'psig.i'            ! rho,psival,qpsi
c      include 'rhog.i'            ! psir,press,pprim,ffprim,qpsir
c      include 'small.i'           ! p(i,j),pmin
c      include 'soln.i'            ! rbp
c      include 'soln2d.i'          ! xi_include
c      include 'solcon.i'          ! time
      include 'storage.i'         ! ydum
c      include 'tdem.i'            ! tdemvb
c      include 'tordlrot.i'        ! iangrot
c      include 'flxav.i'           ! npsi,xmagn1,ymagn1
c
      character*(*) input_string
c
      print *,'in dump_values,tdemvb =',tdemvb
      print *,'called with :',input_string
      if (tdemvb .le. 0)  return  ! normally don't dump the data
c

      iodump = n77          ! shared with DUMP_PSI_VALUES routine
c
      inquire (file = 'debug.data', iostat = iostat,
     .              exist = exists, opened = opened)
           if (iostat .ne. 0) then         ! problem with INQUIRE
             write (ncrt, '(/ 2a, i6)')
     .                    ' ERROR: Fatal INQUIRE failure,',
     .                           ' IOSTAT =', iostat
            call STOP ('subroutine DUMP_VALUES: bad INQUIRE', 12)
           end if
   20      if (exists .and. opened) then   ! file is open
             write (iodump,'(" time = ", 1pe14.6, 1x, a /)')
     .              time, input_string
c
           else                           ! file is not open
            if (exists) then              ! file exists; trash it
               call DESTROY ('debug.data')
             end if
             call getioun(n77,iodump)
             iodump = n77
             open (unit = iodump, file = 'debug.data',
     .             status = 'NEW', iostat = iostat)
             if (iostat .ne. 0) then
               write (ncrt, '(/ 2a, i6)')
     .                      ' ERROR: Fatal OPEN failure,',
     .                             ' IOSTAT =', iostat
               call giveupus(iodump)
               call STOP ('subroutine DUMP_DATA: bad OPEN', 165)
             end if
             call GET_DATE_TIME (timestamp)
             write (iodump, '(2a /)')
     .                      'FILE debug.data CREATED  ', timestamp
             exists = .true.
             opened = .true.
             go to 20
           end if
           write (iodump,'(" fpsi = ",1pe14.6," psivolp = ",1pe14.6,
     .                     " psival = ",1pe14.6)')(fpsi(j),
     .                       psivolp(j),psival(j),j=1,npsi)
           write (iodump,'(" fcap = ",1pe14.6," gcap = ",1pe14.6,
     .                     " hcap = ",1pe14.6)')(fcap(j),
     .                       gcap(j),hcap(j),j=1,nj)
          write (iodump,'(" rcap = ",1pe14.6,
     .                    " rcapi = ",1pe14.6,
     .                    " r2cap = ",1pe14.6,
     .                    " r2capi = ",1pe14.6)')(rcap(j),
     .                      rcapi(j),r2cap(j),r2capi(j),j=1,nj)
c          note rcapi added 1/17/2013 HSJ

c
          write (iodump,'(" roa = ",1pe14.6," r = ",1pe14.6,
     .                    " psir = ",1pe14.6 )')(roa(j),
     .                                 r(j),psir(j),j=1,nj)
c
          write (iodump,'(" drr = ",1pe14.6," rrm = ",1pe14.6,
     .                    " rrp = ",1pe14.6 )')(drr(j),
     .                                 rrm(j),rrp(j),j=1,nj-1)
           write (iodump,'(" eloc = ",1pe14.6," fact = ",1pe14.6,
     .                     " psir = ",1pe14.6)')(eloc(j),
     .                               fact(j),psir(j),j=1,nj)
c
           write (iodump,'(" psival = ",1pe14.6," qpsi = ",1pe14.6,
     .                     " psival = ",1pe14.6)')(psival(j),
     .                               qpsi(j),psival(j),j=1,npsi)
c
           if ((input_string .eq. 'after reqdsk') .or.
     .          (input_string .eq. 'after set_cdf_init') .or.
     .          (input_string .eq. 'from RHOSET')) then
               do j=1,nj
                  zdum(j)=0.0
                  wdum(j)=0.0
               end do
           else
              do j=1,nj
                 zdum(j)=rbp(j)
                 wdum(j)=bp(j)
              end do
           end if
c
       print *,'dump_values'
           write (iodump,'(" q = ",1pe14.6," bp = ",1pe14.6,
     .                     " rbp = ",1pe14.6)')(q(j),
     .                               wdum(j),zdum(j),j=1,nj)
       print *,'returning from dump_values'
      return
c
      end

      subroutine ech_netcdf_interface (pdata_t,nx_t,nz_t,rlim_eqd_t,
     .           zlim_eqd_t,nlim_eqd_t,rp_eqd_t,zp_eqd_t,np_eqd_t,
     .           xdim_t,zdim_t,rmajor_t,redgem_t,zmid_t,
     .           xmaxis_t,zmaxis_t,psimagax_t,psilimitr_t,
     .           btorus_t,time_t)
c ----------------------------------------------------------------------
c load values in calling list with info from netcdf file
c    input
c    nx_t,nz_t     size of mhdgrid (to be checked below)
c    time_t the time at which information is to be take from the newtcdf file
c    all other quantities in the namelist are outputs.
c -------------------------------------------------------------- HSJ ---
c
c
      USE param
      USE io
      USE soln
      USE contour             ! get nconmax for tdem.f90
      USE limiter
      USE mhdpar
      USE mhdgrid
      USE tdem
      USE numbrs
      USE extra
      USE mesh
      USE machin
      USE geom
      USE flags
      USE tordlrot
      USE constnts
      USE soln2d
      USE bd_condtn,only : bctime,ub,fluxb,
     .    ub_save,ub_rho_edge,bctime_zone
      USE psig
      USE rhog
      USE mhdcom
      USE bicube
      USE flxav
      USE neo2d
      USE etc
      USE gpsi
      implicit  integer (i-n), real*8 (a-h, o-z)
      dimension pdata_t(nx_t,nz_t),rp_eqd_t(*),zp_eqd_t(*),
     .          rlim_eqd_t(*),zlim_eqd_t(*)

      include 'storage.i'  ! ydum,zdum

c
c     the following call sets the input to the values at current time
c
      igetcur_t=-1         ! only get relevant stuff at time time_t

      print *,'calling update_tdem_data from ech_necdf_interface'
      print *,'to get tdem data at time ',time_t
      call update_tdem_data (time_t,igetcur_t)
      print *,'done update_tdem_data'


      if (nx_t .ne. nw .or. nz_t .ne. nh)
     .call STOP ('subroutine ECH_NETCDF_INTERFACE: nx, nz problem', 262)
      xdim_t   = xdim         ! lengths in meters
      zdim_t   = ydim         ! meters
      rmajor_t = rmajor/100.  ! meters
      redgem_t = redge        ! meters
      zmid_t   = ymid         ! meters
      xmaxis_t = rma_t/100.0
      zmaxis_t = zma_t/100.0
      psimagax_t = psimag_t
      psilimitr_t = psilim_t
      btorus_t = btor/1.e4    ! tesla
      nlim_eqd_t = nlimiter
      call copya(xlimiter,rlim_eqd_t,nlimiter)
      call copya(ylimiter,zlim_eqd_t,nlimiter)
c       print *,'xlimiter(1) ,ylimiter(1) =',xlimiter(1),ylimiter(1)
c
      np_eqd_t = nplasbdry
      call copya (rplasbdry,rp_eqd_t, nplasbdry)
      call copya (zplasbdry,zp_eqd_t, nplasbdry)
c
c       print *,'rpeqd_t(1),zp-eqd_t(1) =',rp_eqd_t(1),zp_eqd_t(1)
      call copya (p, pdata_t, nw*nh)
      call multpl1 (pdata_t, nw*nh, 1.0/psikgaus) ! convert from..
c                                                 ..kgauss to volt-sec
      return
c
      end




      subroutine get_cdf_data
c
c ----------------------------------------------------------------------
c   load limiter points vectors from tdem mode list
c   get info necessary to construc mhd grid
c
c  INPUT (INCLUDE FILES)
c
c        tdem.i
c             xlim_tdem
c             ylim_tdem
c             rmisc_tdem
c
c  OUTPUT (INCLUDE FILES)
c
c        limiter.i
c             xlimiter
c             ylimiter
c             nlimiter
c
c        mhdgrid.i
c              xdim                       ! meters
c              ydim
c              ymid
c              redge
c              rmajor
c
c ------------------------------------------------------------------ HSJ
c
      USE param
      USE contour
      USE limiter
      USE mhdpar
      USE mhdgrid
      USE tdem
      implicit  integer (i-n), real*8 (a-h, o-z)
c      include 'param.i'
c      include 'contour.i'       ! get nconmax for tdem.i
c      include 'limiter.i'
c      include 'mhdpar.i'
c      include 'mhdgrid.i'
c      include 'tdem.i'
c
      if (nlim_tdem .lt. 5) then   ! 5 some arbitrary number
        call STOP ('subroutine GET_LIMITER: points not loaded', 238)
      else
        do j=1,nlim_tdem
          xlimiter(j) = xlim_tdem(j)             ! in meters
          ylimiter(j) = ylim_tdem(j)
        end do
        nlimiter = nlim_tdem
      end if


c
      rmhdgrid(nw) = -1.0e30
      zmhdgrid(nh) = -1.0e30
      xdim         = real_tdem(1)                       ! meters
      ydim         = real_tdem(2)
      ymid         = real_tdem(3)
      redge        = real_tdem(4)
      rmajor       = real_tdem(5)
 
      return
c
      end

      subroutine load_variable (itype,ncid,id_var,
     .                          start,count,input_file_name,var12name,
     .                          byte,aschar,intg2a,intg4a,real4a,real8a)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c     store the results of the netCDF read
c     in the right kind of variable
c     ONLY THE INTEGER*4 and REAL** WORK FOR NOW!
c ------------------------------------------------------------------ HSJ
c
      integer       start(*), count(*),rcode
      character*(*) input_file_name,var12name
      logical  * 1  byte(*)
      character*(*) aschar(*)  ! need to fix up call
      integer  * 2  intg2a(*)
      integer  * 4  intg4a(*)
      real     * 4  real4a(*)
      real     * 8  real8a(*)
c
      if      (itype .eq. 1) then
        call NCVGT (ncid,id_var,start,count,byte,rcode)
      else if (itype .eq. 2) then
        call NCVGT (ncid,id_var,start,count,aschar,rcode)
      else if (itype .eq. 3) then
        call NCVGT (ncid,id_var,start,count,intg2a,rcode)
      else if (itype .eq. 4) then
        call NCVGT (ncid,id_var,start,count,intg4a,rcode)
      else if (itype .eq. 5) then
        call NCVGT (ncid,id_var,start,count,real4a,rcode)
      else if (itype .eq. 6) then
        call NCVGT (ncid,id_var,start,count,real8a,rcode)
      else
        call STOP ('subroutine LOAD_VARIABLE: ITYPE out of range', 236)
      end if
c
      if (rcode .ne. 0) then
        write (*, '(" Subroutine LOAD_VARIABLE reports:" /
     .              " an error occured while trying to " /
     .              " read variable ", a /
     .              " from the netCDF file ", a,/,
     .              " itype,id_var  = ",i5,i5)')
     .                var12name, input_file_name,itype,id_var
        write(*,'("netcdf routine NCVGT returned rcode =",i5)')rcode
     .          
        call STOP ('subroutine LOAD_VARIABLE: NCVGT problem', 237)
      end if
      return
c
      end

      subroutine read_tdem (itdem_flag, file_name, ishot)
c

c
c ----------------------------------------------------------------------
c
c     subroutine READ_TDEM reads selected parts of the netCDF file
c     whose name is stored in eqdskin. This should be the first

c     subroutine that calls RMEPC_12 so that ncdfile_open is set

c     correctly.
c     The actual data read depends on the setting of itdem_flag.
c
c  INPUT
c     itdem_flag = 1 reads
c
c     itdem_flag = 2 reads
c
c     file_name      character variable contains name of netCDF file
c
c     ishot          shot number of data in netCDF file (not significant
c                    at this time)
c
c  INPUT from include files (i.e., common blocks):
c     param.i
c         ??
c
c     contour.i
c         nconmax        sizes rontr_tdem,etc
c
c     io.i:
c        eqdskin         name of netCDF file
c     mhdpar.i
c        mxtbcmhd        max number of eqdsks, i.e., MaX Time Boundary Cond. MHD
c
c     storage.i
c        xdum            temporary work vector
c        ydum            temporary work vector
c
c     limiter.i
c        maxlimpt        used to dimension xlim_tdem,ylim_tdem
c
c  OUTPUT through  include files (i.e., common blocks)
c
c     tdem.i
c         shot_tdem(i) i=1,2...ntime_tdem list of available shot numbers
c         time_tdem(i) i=1,2...ntime_tdem list of available times
c         xlim_tdem(i) i=1,2... nlim_tdem limiter coords
c         ylim_tdem(i) i=1,2... nlim_tdem
c         time_tdem(1:ntime_tdem)
c         rtime_tdem()
c         timeqbcd()
c         real_misc()
c         int_misc()
c         vol_tdem
c         circum_tdem
c         rma_tdem(),zma_tdem()
c         psilim_tdem()
c         psimag_tdem()
c         rsep_tdem(),zsep_tdem(),psisep_tdem()
c         toteqd_tdem() 
c         totcur()
c         beq_tdem()
c         btorax_tdem()
c         rhoa_tdem()
c         csp_rho_fit_tdem()
c
c ---------------------------------------------------------10/22/96- HSJ
c
      USE param
      USE io
      USE solcon
      USE contour
      USE contour
      USE limiter
      USE mhdpar   
      USE tdem
      USE numbrs  
      USE extra
      USE bd_condtn,only : bctime,ub,fluxb,
     .    ub_save,ub_rho_edge,bctime_zone,totcur
      USE mhdbcdtn
      implicit  integer (i-n), real*8 (a-h, o-z)

      include 'storage.i'

c
      character*(*) file_name
      character*8   slice
c
c     the following are placeholders i the  call to rmepc_12
c
      logical*1    byte
      character*16 aschar  ! need to fix up call
      integer*2    intg2a
      integer*4    intg4a,istart,iend
      real*4       real4a
      real*8       real8a,totcurt_0,totcurt_max

c
      slice = 'none'    ! following are 1d, only time slice is available
      if (itdem_flag .eq. 1) then
c
c         read the vector of shot #s,shot_tdem
c
          call rmepc_12 (iostat,file_name,'shot',mxtbcmhd,
     .               byte,aschar,intg2a,shot_tdem,real4a,real8a,
     .               ishot,itime,itype,ntime_tdem,ncdfile_open,
     .               slice)
c
c         read the vector of available times,time_tdem
          call rmepc_12 (iostat,file_name,'time',mxtbcmhd,
     .               byte,aschar,intg2a,time_tdem,real4a,real8a,
     .               ishot,itime,itype,ntime_tdem,ncdfile_open,
     .               slice)
          do j=1,ntime_tdem
             rtime_tdem(j) = time_tdem (j)*0.001 ! convert to real*8,sec
             timeqbcd  (j) = rtime_tdem(j)
          end do
          itbcmhd = ntime_tdem
c
          ierr_tdem = 0
c
c         CHECK AVAILABLE TIMES:
c
          if (ntime_tdem .lt. 3) then  
            ierr_tdem = ierr_tdem + 1
            write (*, '(" need at least three eqdsks in netCDF file")')
          end if
          ! roundoff issue requires this:
          IF(ABS(time0 -rtime_tdem(1))  .LT. 1.e-6)
     .          rtime_tdem(1) = time0 - 1.e-6
          if ((time0 .lt. rtime_tdem(1))
     .               .or. ( time0 .gt. rtime_tdem(ntime_tdem))
     .               .or. ( timmax .lt. rtime_tdem(1))
     .               .or. ( timmax .gt. rtime_tdem(ntime_tdem))) then
             ierr_tdem = ierr_tdem+1
             write (*,'("error in time settings :"         /
     .             "time0 =",1pe12.5,5x,"timmax =",1pe12.5 /
     .             2x,"time span  in netcdf file:"         /
     .             2x,2(2x,1pe12.5))')time0,timmax,rtime_tdem(1),
     .             rtime_tdem(ntime_tdem)
          end if
c

c
          if (ierr_tdem .gt. 0)
     .      call STOP ('subroutine READ_TDEM: general error', 228)
c
c         read the limiter points,xlim_tdem,ylim_tdem
c
          call rmepc_12 (iostat,file_name,'xlimiter',maxlimpt,
     .                   byte,aschar,intg2a,intg4a,real4a,xlim_tdem,
     .                   ishot,itime,itype,nlim_tdem,ncdfile_open,slice)
c
          call rmepc_12 (iostat,file_name,'ylimiter',maxlimpt,
     .                   byte,aschar,intg2a,intg4a,real4a,ylim_tdem,
     .                   ishot,itime,itype,nlim_tdem,ncdfile_open,slice)
c
c         read the vector of miscelaneous stuff:
c
          call rmepc_12 (iostat,file_name,'real_misc',miscelan,
     .                   byte,aschar,intg2a,intg4a,real4a,real_tdem,
     .                  ishot,itime,itype,nmisc_tdem,ncdfile_open,slice)
c
c         read the integer vector of miscelaneous stuff:
c
          call rmepc_12 (iostat,file_name,'intg_misc',miscelan,
     .               byte,aschar,intg2a,intg_tdem,real4a,real8a,
     .               ishot,itime,itype,nmisc_tdem,ncdfile_open,
     .               slice)
c
c         read total plasma volume:
c
          call rmepc_12 (iostat,file_name,'volume',mxtbcmhd,
     .               byte,aschar,intg2a,intg4a,real4a,vol_tdem,
     .               ishot,itime,itype,ntime_tdem,ncdfile_open,
     .               slice)
c
c         read  plasma circumference:
c
          call rmepc_12 (iostat,file_name,'circum',mxtbcmhd,
     .               byte,aschar,intg2a,intg4a,real4a,circum_tdem,
     .               ishot,itime,itype,ntime_tdem,ncdfile_open,
     .               slice)
c
c         read  plasma magnetic axis:
c
          call rmepc_12 (iostat,file_name,'rma',mxtbcmhd,
     .               byte,aschar,intg2a,intg4a,real4a,rma_tdem,
     .               ishot,itime,itype,ntime_tdem,ncdfile_open,
     .               slice)
c
          call rmepc_12 (iostat,file_name,'zma',mxtbcmhd,
     .               byte,aschar,intg2a,intg4a,real4a,zma_tdem,
     .               ishot,itime,itype,ntime_tdem,ncdfile_open,
     .               slice)
c
c         read  psilim and psimag:
c
          call rmepc_12 (iostat,file_name,'psilim',mxtbcmhd,
     .               byte,aschar,intg2a,intg4a,real4a,psilim_tdem,
     .               ishot,itime,itype,ntime_tdem,ncdfile_open,
     .               slice)

c
          call rmepc_12 (iostat,file_name,'psimag',mxtbcmhd,
     .               byte,aschar,intg2a,intg4a,real4a,psimag_tdem,
     .               ishot,itime,itype,ntime_tdem,ncdfile_open,
     .               slice)

c
c         read rsep and zsep:
c
          call rmepc_12 (iostat,file_name,'rsep',mxtbcmhd,
     .               byte,aschar,intg2a,intg4a,real4a,rsep_tdem,
     .               ishot,itime,itype,ntime_tdem,ncdfile_open,
     .               slice)
c
          call rmepc_12 (iostat,file_name,'zsep',mxtbcmhd,
     .               byte,aschar,intg2a,intg4a,real4a,zsep_tdem,
     .               ishot,itime,itype,ntime_tdem,ncdfile_open,
     .               slice)
c
          call rmepc_12 (iostat,file_name,'psisep',mxtbcmhd,
     .               byte,aschar,intg2a,intg4a,real4a,psisep_tdem,
     .               ishot,itime,itype,ntime_tdem,ncdfile_open,
     .               slice)

c
c         read total plasma current:
c
          call rmepc_12 (iostat,file_name,'toteqd',mxtbcmhd,
     .                   byte,aschar,intg2a,intg4a,real4a,toteqd_tdem,
     .                   ishot,itime,itype,ntime_tdem,ncdfile_open,
     .                   slice)
 
c
c         interpolate toteqd_tdem onto the bctime array:
c
          istart = nbctim
          iend =0
c         at least one value of bctime must be in the interval[time0,timmax]
          do i=1,nbctim
             if(bctime(i) .le. time0)i0 =i
             if(bctime(i) .le. timmax)i0m=i
            if (bctime(i) .ge. time0   .and.                
     .          bctime(i) .le. timmax) then  
                istart = MIN(istart,i)
                iend   = MAX(iend,i)
                call tdem_intrp (toteqd_tdem, ntime_tdem, rtime_tdem,
     .                           bctime(i), totcur(i))
            else                              
             totcur(i) = 0.0                  
            end if                             
          end do
          print *,'istart,iend =',istart,iend
c         new 1/2703

 
          if(i0 .eq. i0m) then   !at least one value of bctime 
                                 !must be in the interval[time0,timmax]
           !get current at time time0:
           call tdem_intrp (toteqd_tdem, ntime_tdem, rtime_tdem,
     .                           time0, totcurt_0)
           !get current at time timmax:
           call tdem_intrp (toteqd_tdem, ntime_tdem, rtime_tdem,
     .                           time0, totcurt_max)
           ! set up bctime to reproduce these currents:
           do i = 1,nbctim     !extend totcur at either end arbitrarily
             totcur(i) = (totcurt_0+totcurt_max)*.5
           enddo
           
          else
           do i = 1,nbctim            !extend totcur at either end arbitrarily
             if(i .lt. istart) totcur(i) = totcur(istart)
             if(i .gt. iend) totcur(i) = totcur(iend)
           enddo
          endif


          totcur(:)  = ABS(totcur(:)) !88888899999 changed to take 
                                      ! absolute value  1/31/13 HSJ
                                      



c
c         read beqd, vac field:
c
          call rmepc_12 (iostat,file_name,'beqd',mxtbcmhd,
     .               byte,aschar,intg2a,intg4a,real4a,beqd_tdem,
     .               ishot,itime,itype,ntime_tdem,ncdfile_open,
     .               slice)
c
c       read btorax, mag axis toroidal field:
c
          call rmepc_12 (iostat,file_name,'btorax',mxtbcmhd,
     .               byte,aschar,intg2a,intg4a,real4a,btorax_tdem,
     .               ishot,itime,itype,ntime_tdem,ncdfile_open,
     .               slice)
c
c         read rho at palsma edge,rhoa:
c
          call rmepc_12 (iostat,file_name,'rhoa',mxtbcmhd,
     .               byte,aschar,intg2a,intg4a,real4a,rhoa_tdem,
     .               ishot,itime,itype,ntime_tdem,ncdfile_open,
     .               slice)
c
c         read smoothed rho at palsma edge,rhoa:
c
          call rmepc_12 (iostat,file_name,'rhomax_smooth',mxtbcmhd,
     .               byte,aschar,intg2a,intg4a,real4a,
     .               rhomax_smooth_tdem,
     .               ishot,itime,itype,ntime_tdem,ncdfile_open,
     .               slice)

c
c         read fitting coefficient vector for rhomax,all fitting models
c
          slice = 'fit2d'
          call rmepc_12 (iostat,file_name,'fit_rhomax',mxtbcmhd,
     .               byte,aschar,intg2a,intg4a,real4a,
     .               csp_rho_fit_tdem,
     .               ishot,itime,itype,ncoldim_csp,ncdfile_open,
     .               slice)
c
      else
c
      end if
      return
c
      end

      subroutine rmepc_12 (iostat,input_file_name,var12name,ndim,
     .                     byte,aschar,intg2a,intg4a,real4a,real8a,
     .                     shot,ifixd,itype,n_length,ncdfile_open,
     .                     slice)
c
c
c ----------------------------------------------------------------------
c
c     RMEPC_12 reads input_file written in netCDF format by subroutine
c     WMEPC_12 in program MEPC.
c
c     input
c             iostat = 0 open netCDF file on entry(if closed)
c                        and close netCDF file on exit (if open)
c                    =-1 close file on exit (if open)
c                    = 1 open file on entry (if closed) and leave open
c                        on exit.
c                    =-1 close file on exit (if open)
c
c             input_file_name
c             var12name  variable name  (must match one of names in
c                        netCDF file exactly)
c
c             shot       shot number (currently not used for anything!)
c             ifixd      subscript in time at which the variable
c                        will be returned as a function of space
c                        or subscript in space at which the variable
c                        will be returned as a function of time
c                        (see defn of slice below).
c             ndim       the DECLARED dimension of the output vector.
c                        The output vector is one of the six quantities
c                        byte,....real8a,see below. ndim is the DECLARED
c                        dimension of whichever one of these six vectors
c                        that is returned (depends on itype,see below).
c                        If ndim is less than the length of the vector
c                        that will be returned by this subroutine then
c                        an error exit is taken!!!!
c
c             slice       = 'space' or 'time', return a vector over
c                         the spatial dimension at a given time or
c                         over the time dimension at a given spatial point.
c     output
c                        output depends on type (i.e., byte, real*8,etc) of
c                        requested variable. the output will be in
c                        one of the following vectors,depending on the
c                        type. The six possible types are
c      integer*2     byte(*)                            returns itype=1
c      character*(*) aschar(*)  ! need to fix up call   returns itype=2
c      integer*2 intg2a(*)                              returns itype=3
c      integer*4 intg4a(*)                              returns itype=4
c      real*4    real4a(*)                              returns itype=5
c      real*8    real8a(*)                              returns itype=6
c                       One way to check where the returned results will
c                       be is to check itype in the calling program.
c                       (for example if itype =3 then the results will be in
c                         vector intg2a(*))
c      n_length          the size of the returned vector
c      ncdfile_open input/output  true/false if file is open/closed
c ---- (ANOTHER!! GREAT SUBROUTINE BY) --- HSJ --- cut the crap --------



      USE param
      USE tdem, only: ncid
      USE flags
      implicit  integer (i-n), real*8 (a-h, o-z)
       include 'netcdf.inc'                ! from the netCDF package..
c                                         ..symlinked to local directory
c
c      include 'param.i'  !KK is all that is needed (so flags.i can be included0
c      include 'flags.i'  !interp_prof
      integer       shot, ifixd, ndim
      logical*1     byte(*)
      logical       ncdfile_open
      character*(*) aschar(*)  ! need to fix up call
      integer*2     intg2a(*)
      integer*4     intg4a(*)
      real*4        real4a(*)
      real*8        real8a(*)
      character*(*) input_file_name , var12name     
      character*(*) slice
c
c     local variables
c
      integer   rcode                               ! error code
      integer   start(3),count(3)                   ! local work space
       integer   vdims(maxvdims),recdim_def
      character*128 dname
c
c      print *,'ncdfile_open,ncid =',ncdfile_open,ncid
      if (iostat .ge. 0 .and. .not. ncdfile_open) then
        call NCPOPT (ncverbos) ! allow netCDF error message
c                                but do not terminate on fatal error
c
c       open input_file, stop if fail
c
        ncid = NCOPN (input_file_name, NCNOWRIT, rcode)


        call save_ncid(ncid)

        if (rcode .ne. 0)
     .  call STOP ('subroutine RMEPC_12: netCDF file open problem', 229)
        ncdfile_open = .true.
      end if
c
c     get some information about this file that will be required
c     when the reads are performed below:
c     number of dimensions defined in input file:          ndims_def
c     number of variables  defined in input file:          nvars_def
c     number of global attrib  "   "   "     "  :          ngatts_def
c     dimension id of unlimited length          :          recdim_def
c     recdim_def returns -1 if no such dimension is defined
c

      call NCINQ (ncid,ndims_def,nvars_def,ngatts_def,recdim_def,rcode)
      if (rcode .ne. 0) then
        call STOP ('subroutine RMEPC_12: NCINQ problem', 230)
      end if
c
c      this puts out an unwanted error message so dont use it HSJ
c      id_var = NCVID (ncid, var12name, rcode) 



c
c      This allows me to control the error: HSJ
       rcode  = NF_INQ_VARID (ncid,var12name,id_var)

       IF(rcode .NE. NF_NOERR)THEN
c        write (*, '("tried to read undefined variable ", a)') var12name
        ! add some logic for rcapi(pre jan17,2013 rcapi not in file)
        IF(var12name .NE. 'rcapi')THEN
           call STOP ('subroutine RMEPC_12: no such variable #1', 231)
        ELSE
           real8a(1:ndim) = .5
           RETURN   
        ENDIF  
      ENDIF

c     now get the id_var (after error check above) HSJ
      id_var = NCVID (ncid, var12name, rcode) 
c
c     given the name of a dimension (here dim_kj) NCDID
c     returns its id no (the number of the dimension
c     defined in the netCDF file):
c
      idim_kj =     NCDID (ncid, 'dim_kj'  , rcode)
      idim_eq =     NCDID (ncid, 'dim_time', rcode)
      idim_nw =     NCDID (ncid, 'dim_nw'  , rcode)
      idim_nh =     NCDID (ncid, 'dim_nh'  , rcode)
      idim_3  =     NCDID (ncid, 'dim_spline_coef'   , rcode)
      idim_neqdsk = NCDID(ncid,'dim_eqdsk',rcode)
c
c     NCDINQ RETURNS dname(='dim_kj  ') and nj:
c     call NCDINQ (ncid, idim_kj, dname, nj, rcode)
c
c     given the id_var of the variable find the following:
c     variable name                        dname  (character)
c     variable type                        itype  (integer  )
c                                              1 = NCBYTE    2 = NCCHAR
c                                              3 = NCSHORT   4 = NCLONG
c                                              5 = NCFLOAT   6 = NCDOUBLE
c     number of dimensions:                nvdims (integer, 0 means scalar, etc)
c     number of assigned attributes        nvatts (integer)
c     vdims(1,2,..nvdims)                         (integer)
c     vdims(i) gives dimension id of dimension i
c
      call NCVINQ (ncid,id_var,dname,itype,nvdims,vdims,nvatts,rcode)
      if (rcode.ne.0) then
         write (*, '("NCVINQ problem while reading variable:", a)')
     .                var12name
         call STOP ('subroutine RMEPC_12: no such variable #2', 232)
      end if
      do i=1,nvdims              ! get the length of each dimension
        call NCDINQ (ncid, vdims(i), dname, count(i), rcode)
      end do
c
 
      if (nvdims .eq. 0) then    ! read scalar
c
c     DO SOMETHING HERE, NO SCALARS IN DATABASE AT PRESENT
c

      else if (nvdims .eq. 1) then ! read vector of values
c
c            the (single) dimension id is vdims(nvdims)
c            use this dimension id to get the length of the vector
c            we want to read.
c            NCDINQ RETURNS dname  and n_length for both fixed and
c            unlimited dimensions given the dimension id.
c
             start(1) = 1
             n_length = count(1)

c
      else if (nvdims .eq. 2) then    ! quantity we want to read is 2-D
c
c           we want to return the vector of values over the space
c           dimension  at time point ifixd (slice = 'space') or
c           we want to return the vector of values over the time
c           dimension at a given spatial point (slice = 'time').
c           The time dimension is identified by 'dim_time' (or idim_eq)
c           the space dimension by 'dim_kj' or idim_kj.
c
         if      (slice .eq. 'space') then
           start(1) = 1        ! start reading of spatial dimension here
           start(2) = ifixd    ! start reading of time dimension here
           count(2) = 1        ! # of values to read in time dimension
           n_length = count(1) ! # of values to read in space dimension
 447    else if (slice .eq. 'time') then
           start(1) = ifixd    ! start reading of spatial dimension here
           start(2) = 1        ! start reading of time dimension here
           count(1) = 1
           n_length = count(2) ! # of values to read in space dimension
         else if (slice .eq. 'fit2d') then
           start(1)=1
           start(2)=1
c           count(1),count(2),the lengths of the two dimensions, were
c           determined above in call to NCDINQ
           n_length=count(1)
         else
           call STOP ('subroutine RMEPC_12: wrong slice setting', 233)
         end if 
c
      else if (nvdims .eq. 3) then
         if (slice .eq. 'fit3d' .and. interp_prof .ne. 'fit_psi') then
           start(1)=1
           start(2)=1
           start(3)=ifixd ! index to rho_12
           count(3)=1  ! return spatial value corresponding to rho_12(j)
         else if (slice .eq. 'fit3d' .and. interp_prof .eq. 'fit_psi')
     .                                                            then
           start(1)=1
           start(2)=1
           start(3)=1
           !count(1..3) is set by NCDINQ above 
         else
c
c           we want to return a two dimensional matrix at a particular
c           value of the third (in this case time) variable:
c           The time dimension is identified by 'dim_time' (or idim_eq)
c           the space dimension by 'dim_nw' or idim_nw
c           and dim_nh or idim_nh .
c           The only 3D
c           quantity is psi (two space ,one time dimension). We do not
c           try to read time slices at one of the nwh space points
c           at present!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
           n_length = count(1)
           start(1) = 1
           start(2) = 1
           start(3) = ifixd     ! start reading of time dimension here
           count(3) = 1         ! only want one value here
        end if
c
      else                      ! more than 3d, not implemented
        call STOP ('subroutine RMEPC_12: only 3D is implemented', 234)
      end if                    ! end nvdims branch
c
        if (n_length .gt. ndim) then
            write (*,'(" subroutine RMEPC_12 reports":                /
     .                 " dimension out of range,"                     /
     .                 " declared dimension =",i5                     /
     .                 " actual size of vector to be returned = ", i5 /
     .                 " variable name = ",a)') ndim,n_length,var12name
            call STOP ('subroutine RMEPC_12: dim out of range', 235)
        end if
c
c       load variable according to its type:
c

        call load_variable (itype,ncid,id_var,
     .                      start,count,input_file_name,var12name,
     .                      byte,aschar,intg2a,intg4a,real4a,real8a)


c
c     close file
c
      if ((iostat .le. 0) .and. ncdfile_open) then
        ncdfile_open = .false.
        call NCCLOS (ncid, rcode)
      end if
      return
c
      end


      subroutine save_ncid(ncid1)
c------------------------------------------------------------------------------------------
c       the only purpose of this subroutine is to get ncid stored in the tdem.i
c       common block
c---------------------------------------------------------------------------HSJ---11-13-00
      USE param
      USE tdem, only: ncid
      implicit  integer (i-n), real*8 (a-h, o-z)

c      include 'param.i'  !
c      include 'contour.i' ! nconmax, used in tdem.i
c      include 'limiter.i' !maxlimpt , used in tdem.i
c      include 'mhdpar.i' ! for parameters in tdem.i
c      include  'tdem.i'  !ncid
     
      ncid = ncid1

      return
      end



      subroutine set_cdf_init
c
c
c --------------------------------------------------------------------
c     this subroutine takes the place of reading the initial eqdsk
c     (reqdsk) in the case of running in the time dependent eqdsk mode.
c     all output is in gaussian units
c
c     INPUT
c        INCLUDE FILES
c           solcon.i
c                 time0                  start of analysis time
c
c ------------------------------------------------------ 6/18/96 --- HSJ
c
      USE param
      USE io
      USE solcon
      USE soln
      USE contour
      USE contour
      USE limiter
      USE mhdpar
      USE mhdgrid
      USE tdem
      USE numbrs
      USE extra
      USE mesh
      USE machin
      USE geom
      USE tordlrot
      USE constnts
      USE bd_condtn,only : bc,bctime,ub,fluxb,
     .    ub_save,ub_rho_edge,bctime_zone,totcur
      USE psig
      USE rhog
      USE mhdcom
      USE etc
      USE gpsi
      USE grid_class,     ONLY : reverse
      implicit  integer (i-n), real*8 (a-h, o-z)
      include 'storage.i'

      dimension    torflux(kstore)
      equivalence (torflux(1), xdum(1))
      character*8  slice
c
c     the following are placeholders i the  call to rmepc_12
c
      logical*1    byte
      character*16 aschar  ! need to fix up call
      integer*2    intg2a
      integer*4    intg4a
      real*4       real4a
      real*8       real8a
c
      xdimeqd = real_tdem(1)                       ! meters
      ydimeqd = real_tdem(2)
      ymideqd = real_tdem(3)
      redeqd  = real_tdem(4)
      reqdsk_box_edge = 100.*redeqd  !to set rin in inject correctly
      reqd    = real_tdem(5)
c
      nxeqd   = intg_tdem(1)
      nyeqd   = intg_tdem(2)
      njeqd   = intg_tdem(3)
c

      if (nj .ne. njeqd) then
        write (ncrt,'(" subroutine SET_CDF_INIT detected an error:" /
     .                " nj in file ",a,2x,i3 /
     .                " nj in Onetwo ",i3 /
     .                " The value of nj in inone must be changed to"
     .                " match the value in file ",a)')
     .                  TRIM(eqdskin),njeqd,nj,eqdskin
        call STOP ('subroutine SET_CDF_INIT: unspecified crud', 239)
      end if
      if (nw .ne. nxeqd) then
         write (ncrt,'(" Subroutine SET_CDF_INIT detected an error:" /
     .                 " nw in file ",a,2x,i3 /
     .                 " nw in Onetwo ",i3 /
     .                 " The value of nw in file",a /
     .                 " must be set to the Onetwo value")')
     .                   eqdskin,nxeqd,nw,eqdskin
         call STOP ('subroutine SET_CDF_INIT: problem city', 240)
      end if
c
      if (nh .ne. nyeqd) then
         write (ncrt,'(" Subroutine SET_CDF_INIT detected an error:" /
     .                 " nh in file ",a,2x,i3 /
     .                 " nh in Onetwo ",i3 /
     .                 " The value of nh in file",a /
     .                 " must be set to the Onetwo value")')
     .                   eqdskin,nyeqd,nh,eqdskin
         call STOP ('subroutine SET_CDF_INIT: bad news', 241)
      end if

      if(intg_tdem(7) .eq. 2 ) then
         write (ncrt,'(" Subroutine SET_CDF_INIT detected an error:" /
     .                 " intg_tdem(7)  in file ",a,2x,i3 /
     .                 " The ICSVKU fitting option is not implemented",
     .                 " in the onetwo code")')
         call STOP ('subroutine SET_CDF_INIT: method not available',241)
      endif
      if(intg_tdem(8) .eq. 2 ) then
         write (ncrt,'(" Subroutine SET_CDF_INIT detected an error:" /
     .                 " intg_tdem(8)  in file ",a,2x,i3 /
     .                 " The ICSVKU fitting option is not implemented",
     .                 " in the onetwo code")')
         call STOP ('subroutine SET_CDF_INIT: method not available',241)
      endif
c
c     find time index
c
      call tableintrp(rtime_tdem,ntime_tdem,time0,il_tdem)
      if (il_tdem .lt. 1 .or. il_tdem .gt. ntime_tdem-1) then
        write (ncrt, '(" subroutine SET_CDF_INIT reports:" /
     .                 " Time interpolation is out of bounds" /
     .                 " time0 = ",1pe12.6 /
     .                 " Start,end times in file ",a,2(2x,1pe14.8))')
     .               time0,eqdskin,rtime_tdem(1),rtime_tdem(ntime_tdem)
        call STOP ('subroutine SET_CDF_INIT: time interp problem', 242)
      end if
c
c ----------------------------------------------------------------------
c     read ntime_tdem time values of the profiles at each of the nw psi
c     values and get the profiles at time0.
c ----------------------------------------------------------------------
c
       slice = 'time'
       do j=1,nw
c
            iostat = 1
c
            call rmepc_12 (iostat,eqdskin,'fpsi',mxtbcmhd,
     .               byte,aschar,intg2a,intg4a,real4a,ydum,
     .               ishot,j,itype,ntime_tdem,ncdfile_open,
     .               slice)

c           given ydum(1:ntime_tdem) on the time grid
c           rtime_tdem(1:ntime_tdem) get the value of each profile
c           at time0 for the nw values that are normally on the eqdsks:


            call tdem_intrp(ydum,ntime_tdem,rtime_tdem,time0,
     .                                            fpsi_tdem(j))
c
            call rmepc_12 (iostat,eqdskin,'presspsi',mxtbcmhd,
     .               byte,aschar,intg2a,intg4a,real4a,ydum,
     .               ishot,j,itype,nw_tdem,ncdfile_open,
     .               slice)
            call tdem_intrp(ydum,ntime_tdem,rtime_tdem,time0,
     .                                          presspsi_tdem(j))
c
           call rmepc_12 (iostat,eqdskin,'ffppsi',mxtbcmhd,
     .               byte,aschar,intg2a,intg4a,real4a,ydum,
     .               ishot,j,itype,nw_tdem,ncdfile_open,
     .               slice)
            call tdem_intrp(ydum,ntime_tdem,rtime_tdem,time0,
     .                                          ffprimpsi_tdem(j))
c
           call rmepc_12 (iostat,eqdskin,'pppsi',mxtbcmhd,
     .               byte,aschar,intg2a,intg4a,real4a,ydum,
     .               ishot,j,itype,nw_tdem,ncdfile_open,
     .               slice)
            call tdem_intrp(ydum,ntime_tdem,rtime_tdem,time0,
     .                                          pprimpsi_tdem(j))
c
           call rmepc_12 (iostat,eqdskin,'qpsi',mxtbcmhd,
     .               byte,aschar,intg2a,intg4a,real4a,ydum,
     .               ishot,j,itype,nw_tdem,ncdfile_open,
     .               slice)
          call tdem_intrp(ydum,ntime_tdem,rtime_tdem,time0,
     .                                          qpsi_tdem(j))
c
      end do     ! end loop over nw psir values


c     for output to statefile
      IF(.NOT. ALLOCATED( pppsi_eqdsk))ALLOCATE(pppsi_eqdsk(nw))
      IF(.NOT. ALLOCATED( presspsi_eqdsk))ALLOCATE(presspsi_eqdsk(nw))
      IF(.NOT. ALLOCATED( fpsi_eqdsk))ALLOCATE(fpsi_eqdsk(nw))
      IF(.NOT. ALLOCATED( ffppsi_eqdsk))ALLOCATE(ffppsi_eqdsk(nw))
      IF(.NOT. ALLOCATED( qpsi_eqdsk))ALLOCATE(qpsi_eqdsk(nw))
      pppsi_eqdsk(1:nw)  = pprimpsi_tdem(1:nw)
!      CALL reverse(nw,pppsi_eqdsk)
      presspsi_eqdsk(1:nw)  = presspsi_tdem(1:nw)
!      CALL reverse(nw,presspsi_eqdsk)
      fpsi_eqdsk(1:nw)   = fpsi_tdem(1:nw)
!      CALL reverse(nw,fpsi_eqdsk)
      ffppsi_eqdsk(1:nw) = ffprimpsi_tdem(1:nw)
!     CALL reverse(nw,ffppsi_eqdsk)
      qpsi_eqdsk(1:nw)   = qpsi_tdem(1:nw)
!      CALL reverse(nw,qpsi_eqdsk)
      nxeqd_eqdsk           = nw


      call tdem_intrp(psilim_tdem,ntime_tdem,rtime_tdem,time0,
     .                                          psilim)
      call tdem_intrp(psimag_tdem,ntime_tdem,rtime_tdem,time0,
     .                                          psimag)
c      psilim_eqdsk          = psilim_tdem(?)
c      psimag_eqdsk          = psimag_tdem(?)

c
c     more time interpolation,this time on the nj grid:
      do j=1,nj

          call rmepc_12 (iostat,eqdskin,'curden',mxtbcmhd,
     .                   byte,aschar,intg2a,intg4a,real4a,ydum,
     .                   ishot,j,itype,nw_tdem,ncdfile_open,
     .                   slice) ! amps/m**2
          call tdem_intrp(ydum,ntime_tdem,rtime_tdem,time0,
     .                                          curden_tdem(j))
c
           call rmepc_12 (iostat,eqdskin,'hcap',mxtbcmhd,
     .               byte,aschar,intg2a,intg4a,real4a,ydum,
     .               ishot,j,itype,nw_tdem,ncdfile_open,
     .               slice)
          call tdem_intrp(ydum,ntime_tdem,rtime_tdem,time0,
     .                                          hcap_tdem(j))
c
           call rmepc_12 (iostat,eqdskin,'r2cap',mxtbcmhd,
     .               byte,aschar,intg2a,intg4a,real4a,ydum,
     .               ishot,j,itype,nw_tdem,ncdfile_open,
     .               slice)
          call tdem_intrp(ydum,ntime_tdem,rtime_tdem,time0,
     .                                          r2cap_tdem(j))
      end do   ! end time interpolation on nj grid


      convrt = 1.0e-4
      call multpl1 (curden_tdem, nj, convrt)    ! amps/cm**2
c
      call copya(curden_tdem,curden,nj)         !curden in analysis mode
      call copya(hcap_tdem,hcap,nj)             !is obtained here
      call copya(r2cap_tdem,r2cap,nj)

c     load the onetwo variables:
c
            call copya (fpsi_tdem,fpsi,nw)
            call copya (ffprimpsi_tdem,ffppsi,nw)
            call copya (presspsi_tdem,presspsi,nw)
            call copya (pprimpsi_tdem,pppsi,nw)
            call copya (qpsi_tdem,qpsi,nw)
            call tdem_intrp(rma_tdem,ntime_tdem,rtime_tdem,
     .                                           time0,rma)
            call tdem_intrp(zma_tdem,ntime_tdem,rtime_tdem,
     .                                           time0,zma)
            call tdem_intrp(psimag_tdem,ntime_tdem,rtime_tdem,
     .                                           time0,psimag)
            call tdem_intrp(psilim_tdem,ntime_tdem,rtime_tdem,
     .                                           time0,psilim)

            call tdem_intrp(beqd_tdem,ntime_tdem,rtime_tdem,
     .                                           time0,beqd)
            call tdem_intrp(toteqd_tdem,ntime_tdem,rtime_tdem,
     .                                           time0,toteqd)
            call tdem_intrp(psisep_tdem,ntime_tdem,rtime_tdem,
     .                                           time0,psimx1)
            call tdem_intrp(rsep_tdem,ntime_tdem,rtime_tdem,
     .                                           time0,rsep)
            call tdem_intrp(zsep_tdem,ntime_tdem,rtime_tdem,
     .                                           time0,zsep)
           call tdem_intrp(vol_tdem,ntime_tdem,rtime_tdem,
     .                                           time0,eqdskvol)

           call tdem_intrp(rhoa_tdem,ntime_tdem,rtime_tdem,
     .                                           time0,rhoa0)
c
            psimx2 = 0.0
            xax1   = rma
            xax2   = 0.0
            zax1   = zma
            zax2   = 0.0
            psisep = psimx1
c
            write (ncrt, '(1x, a, f12.6)')
     .      'volume determined from eqdsk boundary values = ', eqdskvol

            slice = 'none'
c
c          do quick linear interp on psi
c
           call rmepc_12 (iostat,eqdskin,'psi',nw,
     .               byte,aschar,intg2a,intg4a,real4a,psi_tdem,
     .               ishot,il_tdem,itype,nw_tdem,ncdfile_open,
     .               slice)    ! at time point il_tdem
          call rmepc_12 (iostat,eqdskin,'psi',nw,
     .               byte,aschar,intg2a,intg4a,real4a,p,
     .               ishot,il_tdem+1,itype,nw_tdem,ncdfile_open,
     .               slice) ! at time point il_tdem+1, use p temporarily
          delt = time0-rtime_tdem(il_tdem)
          delt = delt/(rtime_tdem(il_tdem+1)-rtime_tdem(il_tdem))
          do j=1,nh
             do i=1,nw
                psi(i,j) = psi_tdem(i,j)+(p(i,j)-psi_tdem(i,j))*delt
             end do
          end do
c
c ----------------------------------------------------------------------
c convert values to gauss cm units
c ----------------------------------------------------------------------
c
      xdimeqd = 100.0 * xdimeqd
      ydimeqd = 100.0 * ydimeqd
      reqd    = 100.0 * reqd
      redeqd  = 100.0 * redeqd
      ymideqd = 100.0 * ymideqd
      rma     = 100.0 * rma
      zma     = 100.0 * zma
      xax(1)  = rma
      yax(1)  = zma
      rsep    = 100.0 * rsep
      zsep    = 100.0 * zsep
      psimag  = psikgaus*1.0e3*psimag ! gauss-cm**2, psikgaus*1000=1.0e8
      psilim  = psikgaus*1.0e3*psilim
      call multpl1 (psi, nwh, psikgaus*1.0e3)
      call copya ( psi,p,nwh)
      beqd    = 1.0e4 * beqd                      ! gauss
      cconst  = 1.0e+6
      call multpl1 (fpsi     , nxeqd    , cconst) ! fpsi in gauss-cm
      cconst  = 1.0e+4
      call multpl1 (ffppsi   , nxeqd    , cconst) ! ffppsi in gauss
      cconst  = 10.0
      call multpl1 (presspsi , nxeqd    , cconst) ! presspsi: erg/cm**3
      cconst  = 1.0e-7
      call multpl1 (pppsi    , nxeqd    , cconst) ! pppsi in..
c                                              ..gram/gauss*cm**3*sec**2
c
c      if time0 is one of the times in rtime_tdem then
c      read in the contours and use them . Otherwise
c      use the time interpolated psi to generate a new
c      plasma contour (interpolation in time of the contours
c      requires that we evaluate the velocity field,
c      not worth the effort here)
c
       cconst  = 100.0
       ircontr=0
       ncontr=0  ! forces contour to be found in chekeqd
c
       if (ABS (rtime_tdem(il_tdem)-time0) .lt. 1.0e-6) then
         ircontr=1
         iuse=il_tdem
       else if (ABS (rtime_tdem(il_tdem+1)-time0) .lt. 1.0e-6) then
         ircontr = 1
         iuse    = il_tdem+1
       else
         ircontr=0
       end if
       if (ircontr .eq. 1) then
           slice = 'space'
           call rmepc_12 (iostat,eqdskin,'rcontr',nconmax,
     .          byte,aschar,intg2a,intg4a,real4a,rcontr_tdem,
     .          ishot,iuse,itype,ncontrr,ncdfile_open,
     .          slice)
           call rmepc_12 (iostat,eqdskin,'zcontr',nconmax,
     .          byte,aschar,intg2a,intg4a,real4a,zcontr_tdem,
     .          ishot,iuse,itype,ncontrr,ncdfile_open,
     .          slice)
          call rmepc_12 (iostat,eqdskin,'ncontr',mxtbcmhd,
     .          byte,aschar,intg2a,ncontr_tdem,real4a,real8a,
     .          ishot,iuse,itype,ntime_tdem,ncdfile_open,
     .          slice)
c
           ncontr=ncontr_tdem(iuse)
           call copya (rcontr_tdem,rcontr, ncontr)
           call copya (zcontr_tdem,zcontr, ncontr)
           call multpl1 (rcontr   , ncontr   , cconst) ! cm
           call multpl1 (zcontr   , ncontr   , cconst) ! cm
           nplasbdry=ncontr
           call copya (rcontr,rplasbdry, nplasbdry)
           call copya (zcontr,zplasbdry, nplasbdry)
c
       end if
c
      rmajor    = reqd
      btor      = beqd
      flim      = rmajor*btor
      totcur(1) = toteqd
      rbp(nj)   = 0.2 * ABS (toteqd)
      if (time .eq. bctime(1))  bc(1,nk-iangrot) = 0.2 * ABS (toteqd)
c
c     get contour with time interpolated psi
c     here is where we are not able to use the eqdsk
c     plasma boundary since that information does not
c     exist at the desired time
c

 
      call chekeqd (psi, nyeqd, xdimeqd, ydimeqd, redeqd, ncrt,
     .              rcontr, zcontr, ncontr, btor, rma, zma,
     .              psimag, psilim, totcur(1), rmajor, nconmax)

       psiaxis = psimag
       psibdry = psilim
       nplasbdry = ncontr
       call copya(rcontr,rplasbdry,ncontr)
       call copya(zcontr,zplasbdry,ncontr)
c
c --- set up the psival grid: psival(1) = edge, psival(nw) = axis
c
      iunfrm = 1
      call psiset (nw, psiaxis, psibdry, psival, iunfrm)
c
c --- save eqdsk pressure and psi for plotting (subroutine EQPLOT)
c
      do j=1,nxeqd
        psieqdsk(j) = psival  (nxeqd-j+1) * 1.0e-08    ! volt-sec
        prseqdsk(j) = presspsi(nxeqd-j+1) * 0.1        ! nt/m**2
      end do
c
c ----------------------------------------------------------------------
c --- calculate the toroidal flux
c ----------------------------------------------------------------------
c
      call trap2 (psival, qpsi, torflux, nw)
      do j=1,nw
        torflux(j) = twopi * (torflux(j) - torflux(nw))
      end do
c
c ----------------------------------------------------------------------
c --- define the rho grid
c ----------------------------------------------------------------------
c
      do i=1,nw
        rho(i) = SQRT (torflux(i)/(pi*btor))
      end do
c
c ----------------------------------------------------------------------
c --- define the r grid
c ----------------------------------------------------------------------
c
****  do i=1,nj
****    r(i) = r(i) * rho(1) / r(nj)
****  end do
c
c     mepc code puts out fixed in time zeta grid values:
c
      write (*,'("rhoa from tdem,rho(1) calc =",2(2x,1pe12.2))')
     .            rhoa0,rho(1)
      dtdemr = 1.0 / (nj - 1)
      do i=1,nj
         r(i)= (i-1)*dtdemr*rhoa0
      end do
c
c ----------------------------------------------------------------------
c --- get q on r grid
c ----------------------------------------------------------------------
c
  605 call intrp (0, 1, rho, qpsi, nw, r, q, nj)
c
c --- calculate bp0
c ----------------------------------------------------------------------
c
      do j=1,nj
        bp(j) = r(j)*btor/(rmajor*q(j))
      end do
c
c ----------------------------------------------------------------------
c --- now define the psir grid(this is the psi grid corresponding to
c --- the r grid. it is defined by eq. 2.1-8)
c ----------------------------------------------------------------------
c
      call trap2 (r, bp, psir, nj)
      do  j=1,nj
         psir(j) = psival(nw) + (psival(1) - psival(nw))*psir(j)
     .                                                   /psir(nj)
      end do
c
c ----------------------------------------------------------------------
c --- get pressure on psir grid
c ----------------------------------------------------------------------
c
      call intrp (1, 1, psival, presspsi, nw, psir, press, nj)
c
c ----------------------------------------------------------------------
c --- get pprim on psir grid
c ----------------------------------------------------------------------
c
      call intrp (1, 1, psival, pppsi, nw, psir, pprim, nj)
****  call difydx(psir,press,pprim,nj)
c
c ----------------------------------------------------------------------
c --- get ffprim on psir grid
c ----------------------------------------------------------------------
c
      call intrp (1, 1, psival, ffppsi, nw, psir, ffprim, nj)
c
c ----------------------------------------------------------------------
c --- get qpsir on psir grid
c ----------------------------------------------------------------------
c
      call intrp (1, 1, psival, qpsi, nw, psir, qpsir, nj)
c
c ----------------------------------------------------------------------
c --- convert psi(i,j),p(i,j),psival(i),psir(j) to kgauss cm**2
c --- convert bp(j) to kgauss
c ----------------------------------------------------------------------
c
      cconst  = 1.0e-3
      call multpl1 (bp, nj, cconst)
      psibdry = psilim  * 1.0e-3
      psiaxis = psiaxis * 1.0e-3
      call multpl1 (psival, nw , cconst)
      call multpl1 (psir  , nj , cconst)
      call multpl1 (psi   , nwh, cconst)
      call copya   (psi   , p  , nwh   ) ! p must be set for irguess < 0
c
c --- SAVCUR copies pprim into ppold, ffprim into ffpold,
c --- curden into curold, psir into psiold, and r into rold,
c --- so that subroutine CONCUR can do its thing
c
      call savcur
      return
c
      end

      subroutine tdem_diffwrt_time (v_dep, nval, v_ind, v_want, v_ret,
     .                              num_derivs, vderiv_ret, vdderiv_ret)
c

c
c ----------------------------------------------------------------------
c     get the value of v_dep at v_want and
c     get the derivative of vdep wrt v_ind at v_want
c     this routine was created to allow easy changing from linear to
c     spline interpolation
c INPUT
c     v_dep(i)  i=1,2...nval dependent variable
c     v_ind(i)               independent
c     v_want    value in range of v_ind at which v_dep is desired
c     num_derivs   # of derivatives wanted ,0,1,or 2
c     MAY use vdum from storage.i !!
c OUTPUT
c     v_ret         interpolated value
c     vderiv_ret    derivative
c     vdderiv_ret   second derivative
c
c --------------------------------------------------- 6/23/96-----HSJ
c
      USE param
      USE mhdpar
      USE flags
      USE bicube
      USE replace_imsl,                 ONLY : my_ibcccu,my_dbcevl1,
     .                                         my_icsicu,my_icsevu,
     .                                         my_icsvku
      implicit  integer (i-n), real*8 (a-h, o-z)
      character*16 interp_type
c
c      include 'param.i'
c      include 'mhdpar.i'   ! nw, mxtbcmhd
c      include 'bicube.i'   ! wnoperm
c      include 'flags.i'    ! interp_prof
      include 'storage.i'
c
****  parameter   (nwt = MAX (nw, mxtbcmhd))       ! illegal in CRAY f77
      parameter   (nwt =          mxtbcmhd)
      parameter   (nxkmax = 100)
c
      dimension xk(nxkmax),yk(nxkmax)
      dimension    v_dep(*), v_ind(*), csp(nwt,3),work_icsvku(1)
c      equivalence (csp(1,1), wnoperm(1))
      equivalence (work_icsvku(1),vdum(1))
      dimension    bpar(4) ,vwa(1),vwr(1)
      data         bpar /0.0, 0.0, 0.0, 0.0/
      data         nxkd /6 /
c
      if (MAX (nw, mxtbcmhd) .ne. nwt) then
        call STOP ('subroutine TDEM_DIFFWRT_TIME: bad parameter', 243)
      end if
c
      interp_type = 'nat spline'             ! because bpar = 0
c
c     all kinetic input data was smoothed outside of Onetwo
c     The poloidal B field must be smoothed here however so that
c     its derivative in time can be taken in a meaningful manner:
c
****  if (interp_prof .eq. 'bp0'   )  interp_type = 'spline_fit'
****  if (interp_prof .eq. 'psilim')  interp_type = 'spline_fit'
c
      if (nval .le. 3)  interp_type = 'nat spline'
c
      if (interp_type .eq. 'linear') then    ! try linear interp. later?
        call STOP ('subroutine TDEM_DIFFWRT_TIME: linear is bad', 10)
      else if (interp_type .eq. 'nat spline') then! NATURAL cubic spline
        call my_icsicu (v_ind,v_dep,nval,bpar,csp,nwt,ier)

        ! for ifort:
        vwa(1) = v_want 
        call my_icsevu (v_ind,v_dep,nval,csp,nwt,vwa,vwr,1,ier)
        !call my_icsevu (v_ind,v_dep,nval,csp,nwt,v_want,v_ret,1,ier)
        v_ret = vwr(1)

      else if (interp_type .eq. 'spline_fit') then
c       dont smooth the data value itself:
        call my_icsicu (v_ind,v_dep,nval,bpar,csp,nwt,ier)
        vwa(1) = v_want 
        call my_icsevu (v_ind,v_dep,nval,csp,nwt,vwa,vwr,1,ier)
        v_ret = vwr(1)
        !call my_icsevu (v_ind,v_dep,nval,csp,nwt,v_want,v_ret,1,ier)

        if (num_derivs .eq. 0)  return
         nxk = nxkd
         nxk = MAX (nxk, 3)                   ! number of knots to use
         if (nxk .gt. nxkmax)
     .   call STOP ('subroutine TDEM_DIFFWRT_TIME: nxkmax small', 256)
         if (nval .le. 3)
     .   call STOP ('subroutine TDEM_DIFFWRT_TIME: nval too small', 257)
         xk(1)=v_ind(1)               ! initial guess for knot locations
         xk(nxk)=v_ind(nval)
         dxk=(xk(nxk)-xk(1))/(nxk-1)
         do j=2,nxk-1
           xk(j)=xk(j-1)+dxk
         end do
        work_dim = nval * (nxk + 6)
        if (work_dim .gt. kstore)
     .  call STOP ('subroutine TDEM_DIFFWRT_TIME: insuff. storage', 258)
c
        call my_icsvku (v_ind,v_dep,nval,xk,nxk,yk,csp,nwt,error,
     .               work_icsvku,ier)
        if (ier .gt. 100)
     .  call STOP ('subroutine TDEM_DIFFWRT_TIME: ISCVKU problem', 259)
c
        call tableintrp(xk,nxk,v_want,i1)
        deltav = v_want-xk(i1)
        v_rets  = ((csp(i1,3)*deltav +csp(i1,2))*deltav+
     .              csp(i1,1))*deltav+yk(i1)
        v_ret=v_rets
        if (num_derivs .eq. 0) return
        vderiv_ret = (3.0*csp(i1,3)*deltav
     .              + 2.0*csp(i1,2))*deltav + csp(i1,1)
        if (num_derivs .eq. 1)  return
        vdderiv_ret = 6.0*csp(i1,3)*deltav+2.0*csp(i1,2)
        return
      else
        call STOP ('subroutine TDEM_DIFFWRT_TIME: bad INTERP_TYPE', 11)
      end if
c
      if (num_derivs .eq. 0)  return
c
      call tableintrp(v_ind,nval,v_want,il)
      if (il .le. 0 .or. il .ge. nval) then
        call STOP ('subroutine tdem_diffwrt_time: problem #1', 200)
      end if
      deltav     = v_want-v_ind(il)
      vderiv_ret = (3.0*csp(il,3)*deltav+ 2.0*csp(il,2))*deltav
     .                                                  +csp(il,1)
      if (num_derivs .eq. 1)  return
c
      vdderiv_ret = 6.0*csp(il,3)*deltav+2.0*csp(il,2)
      return
c
      end

      subroutine tdem_intrp (v_dep, nval, v_ind, v_want, v_ret)
c
      implicit none
c
c ----------------------------------------------------------------------
c     interpolate v_dep
c INPUT
c     v_dep(i)  i=1,2...nval dependent variable
c     v_ind(i)               independent
c     v_want    value in range of v_ind at which v_dep is desired
c OUTPUT
c     v_ret     interpolated value
c
c ------------------------------------------------------------------ HSJ
c
      integer nval, num_derivs
      real*8  v_dep(*), v_ind(*), v_want, v_ret
c
      num_derivs = 0
c
c     TDEM_DIFFWRT_TIME is called with 2 fewer arguments:
c
      call tdem_diffwrt_time (v_dep, nval, v_ind, v_want, v_ret,
     .                        num_derivs)
      return
c
      end

      subroutine update_tdem_data (t_interparg, igetcur)
c
c ----------------------------------------------------------------------
c this suboutine provides the parameters given below at the time
c t_interparg by interplation in the netCDF eqdsk input file.
c derivatives (in time) of required parms are also returned
c Note that the units are returned suitable for use in TPORT!!!
c if this is analysis mode for Faraday's law then:
c       If igetcur =1 then only current related quantities are obtained
c       If igetcur =0            both mhd and current
c       If igetcur =-1           mhd
c ------------------------------------------------------ 6/25/96 --- HSJ
c
      USE param
      USE io
      USE soln
      USE contour
      USE limiter
      USE mhdpar
      USE mhdgrid
      USE tdem
      USE numbrs
      USE extra
      USE mesh
      USE machin
      USE geom
      USE flags
      USE tordlrot
      USE constnts
      USE soln2d
      USE bd_condtn,only : bctime,ub,fluxb,
     .    ub_save,ub_rho_edge,bctime_zone
      USE psig
      USE rhog
      USE mhdcom
      USE bicube
      USE flxav
      USE neo2d
      USE etc
      USE gpsi        ! p(nwmhd,nhmhd) ,*_eqdsk
      USE grid_class,     ONLY : reverse
      USE tension_spline,             ONLY : tspline90,t716_TSPSI,             
     .                                       t716_TSVAL1,t716_TSPSS
      USE statistics,                 ONLY : descrip_mon,tdem_mon_set,
     .                                       tdem_index,last_mon_index,
     .                                       start_timer,stop_timer,
     .                                       collect_stats,elapsed_time,
     .                                       o12_index


      implicit  integer (i-n), real*8 (a-h, o-z)

      include 'storage.i'  ! ydum,zdum
c                            (NOTE: VDUM is used in TDEM_DIFFWRT_TIME)

c
      character*8 slice
c
c     the following are placeholders in the call to RMEPC_12
c
      logical  * 1 byte
      character*16 aschar  ! need to fix up call
      character*128 dname
      integer  * 2 intg2a
      integer  * 4 intg4a , vdims(3),rcode, id_var,itype,
     .             nvdims,nvatts, count(3),nj_ncd,iendc,iflag,
     .             ncd
      real     * 4 real4a
      REAL*8   deriv_bc(4),v_interp(4),dv_interp(4)
      real *8 ,DIMENSION(:,:,:),ALLOCATABLE,save :: csp_psi_fit_tdem
****  real     * 8 real8a  ! not used
c      data tol /1.0e-9/    ! poor, make same as other tols in program
      data tol /1.0e-7/     ! HSJ 6/13/13


           if(o12_index .Gt. 1) THEN
             write(888,fmt='("o12 index =",i5)')o12_index
             call stop("o12 index prob c206,l 1825",1) ! 8888889999
          ENDIF
      ! start collecting time info here
      IF(.NOT. tdem_mon_set)THEN  
         tdem_index = last_mon_index+1
         last_mon_index = tdem_index
         tdem_mon_set = .TRUE.
         descrip_mon(tdem_index) = "tdem  elapsed time"
      ENDIF
          if(o12_index .Gt. 1) THEN
             write(888,fmt='("o12 index =",i5)')o12_index
             call stop("o12 index prob c206,l 1836",1) ! 8888889999
          ENDIF
      start_timer(tdem_index) = .TRUE.
      stop_timer(tdem_index)  = .FALSE.
      CALL collect_stats(tdem_index)
          if(o12_index .Gt. 1) THEN
             write(888,fmt='("o12 index =",i5)')o12_index
             call stop("o12 index prob c206,l 1843",1) ! 8888889999
          ENDIF
      start_timer(tdem_index) = .FALSE. 
c      write(888,FMT='("start tdem 1833",1pe14.6)')
c     .      elapsed_time(tdem_index)  ! 88888889999999
c
      iostat = 1


c       figure out what size csp_psi_fit should be
c
      if( .not. allocated(csp_psi_fit_tdem))then
              !replaces  rmepc functionality for 3d cases:
              !get variable ID:
              id_var = NCVID (ncid, 'fit_psi', rcode)
              !get # dimensions,nvdims,dimension id, vdims(1,..nvdims)
              call NCVINQ (ncid,id_var,dname,itype,nvdims,vdims,
     .                                                  nvatts,rcode)
              ! get the length of each dimension
              do i=1,nvdims              
                   call NCDINQ (ncid, vdims(i), dname, count(i), rcode)
              end do
         allocate( csp_psi_fit_tdem(count(1),count(2),count(3)))
      endif
c
c





c     find time index
c
      t_interp=t_interparg
      call tableintrp(rtime_tdem,ntime_tdem,t_interp,il_tdem)
      if (il_tdem .lt. 1 .or. il_tdem .gt. ntime_tdem-1) then
          if (ABS (t_interp-rtime_tdem(ntime_tdem)) .lt. 3.0*tol) then
            il_tdem=ntime_tdem-1
            t_interp=rtime_tdem(ntime_tdem)
          else if (ABS (t_interp-rtime_tdem(1)) .lt. 3.0*tol) then
            il_tdem=1
            t_interp=rtime_tdem(1)
          else
            write (ncrt, '(" subroutine UPDATE_TDEM_DATA reports:" /
     .                     " Time interpolation is out of bounds" /
     .                     " t_interp = ",1pe12.6 /
     .                   " Start, end times in file ",a,2(2x,1pe14.8))')
     .             t_interp,eqdskin,rtime_tdem(1),rtime_tdem(ntime_tdem)
!        call STOP ('subroutine UPDATE_TDEM_DATA: T_INTERP problem', 9)
            difs= t_interp-rtime_tdem(1) 
            dife = t_interp-rtime_tdem(ntime_tdem)
            IF(ABS(difs) .lt. ABS(dife))THEN
               write(ncrt,FMT='("correcting by adjusting ",/,
     .                           "interpolation time from ",/,
     .                           1pe16.8,/,
     .                           "  to:"
     .                           1pe16.8)')t_interp,rtime_tdem(1)
               il_tdem=1
               t_interp=rtime_tdem(1)
            ELSE
               write(ncrt,FMT='("correcting by adjusting",/
     .                          "interpolation time from ",/,
     .                           1pe16.8,/,
     .                           "  to:"
     .                           1pe16.8)')t_interp,
     .                           rtime_tdem(ntime_tdem)
               il_tdem=ntime_tdem-1
               t_interp=rtime_tdem(ntime_tdem)
            ENDIF
          end if
      end if
c
c ----------------------------------------------------------------------
c     read ntime_tdem time values of the profiles at the fixed (in time)
c     normalized zeta grid point rho_12(j). The actual rho grid is
c     given by rho_12(j)*rhoa(at this time)
c     do  time interpolation at each of  the nj spatial grid points
c     to get the actual desired values:
c ----------------------------------------------------------------------
c
      slice = 'time'
c
c          get a value of  rhomax at this time:
c
                  if (intg_tdem(7) .le. 0) then      ! no smoothing
c
c                     get drhoadt_geom by differentiating the spline
c                     rhoa_tdem was loaded in read_tdem for the times
c                     given in rtime_tdem(1,2,..ntime_tdem):
c
                      interp_prof='drhoadt_geom'
                      call tdem_diffwrt_time (rhoa_tdem, ntime_tdem,
     .                                        rtime_tdem, t_interp,
     .                                        rhoa, 1, drhoadt_geom,
     .                                        dumy)
                      drhoadt_geom = drhoadt_geom* 100.0 ! convert to cm
                      drhomaxdt_tdem = drhoadt_geom
                      rhoa   = rhoa  * 100.0
c
c                     fit  coefficients,csp_rho_fit_tdem,  for all models
c                     were read in read_tdem
c
                  else if (intg_tdem(7) .eq. 1) then ! linear fit
                           rhoa=(csp_rho_fit_tdem(1,1)
     .                             +csp_rho_fit_tdem(2,1)*t_interp)*100.
                  else if (intg_tdem(7) .eq. 2) then ! spline fit icsvku
c
c                 this option found to be too unstable, not implemented
                       call STOP('ICSVKU fit of tdem data not available
     . in Onetwo',0)
c
                  else if (intg_tdem(7) .eq. 3) then ! spline fit icsscv
                     delta_time=t_interp-rtime_tdem(il_tdem)
                   rhoa=100.0*(((csp_rho_fit_tdem(il_tdem,3)*delta_time+
     .                       csp_rho_fit_tdem(il_tdem,2))*delta_time+
     .                       csp_rho_fit_tdem(il_tdem,1))*delta_time+
     .                       rhomax_smooth_tdem(il_tdem))
                     drhomaxdt_tdem = 100*
     .                   ((3.0*csp_rho_fit_tdem(il_tdem,3) *delta_time+
     .                     2.0*csp_rho_fit_tdem(il_tdem,2))*delta_time+
     .                         csp_rho_fit_tdem(il_tdem,1))
                  else
c
                  end if
c

      do j=1,nj 
c
c          if (itran(nk-iangrot) .eq. 0
c     .       .and. igetcur .ge. 0 ) then ! analysis mode for rbp
c
              call rmepc_12 (iostat,eqdskin,'curden',mxtbcmhd,
     .                       byte,aschar,intg2a,intg4a,real4a,ydum,
     .                       ishot,j,itype,nw_tdem,ncdfile_open,
     .                       slice)    ! <J_phi*R0/R>,amps/m**2
              interp_prof='cur'

              call tdem_intrp(ydum,ntime_tdem,rtime_tdem,t_interp,
     .                                             curden_tdem(j))
              curden_tdem(j) = 1.e-4*curden_tdem(j)   !convert to amps/cm**2
c

              call rmepc_12 (iostat,eqdskin,'bp0',mxtbcmhd,
     .               byte,aschar,intg2a,intg4a,real4a,ydum,
     .               ishot,j,itype,nw_tdem,ncdfile_open,
     .               slice) ! tesla
              interp_prof='bp0'
              call tdem_diffwrt_time (ydum, ntime_tdem, rtime_tdem,
     .                                t_interp, bp0_tdem(j), 1,
     .                                dbdt_tdem(j), dumy)
c
              interp_prof='q'
              call rmepc_12 (iostat,eqdskin,'q',mxtbcmhd,
     .               byte,aschar,intg2a,intg4a,real4a,ydum,
     .               ishot,j,itype,nw_tdem,ncdfile_open,
     .               slice)
              call tdem_intrp(ydum,ntime_tdem,rtime_tdem,t_interp,
     .                                                  q_tdem(j))

              interp_prof='psiconz'
              call rmepc_12 (iostat,eqdskin,'psi_con_zeta',mxtbcmhd,
     .               byte,aschar,intg2a,intg4a,real4a,ydum,
     .               ishot,j,itype,nw_tdem,ncdfile_open,slice)
              call tdem_diffwrt_time (ydum, ntime_tdem, rtime_tdem,
     .                                t_interp,
     .                                const_zeta_psigrid_smooth_tdem(j),
     .                                1, dpsidt_fin_dif, dumy)
              slice='fit3d'
              interp_prof='fit_psi'

              call rmepc_12 (iostat,eqdskin,'fit_psi',mxtbcmhd,
     .               byte,aschar,intg2a,intg4a,real4a,
     .               csp_psi_fit_tdem,
     .               ishot,j,itype,nw_tdem,ncdfile_open,slice)

              delta_time=t_interp-rtime_tdem(il_tdem)

              if (intg_tdem(8) .eq. 0.0) then ! finite difference form
                 dpsidt_const_zeta_tdem(j)=dpsidt_fin_dif
              else if (intg_tdem(8) .eq. 1) then ! fitted forms
                 dpsidt_const_zeta_tdem(j)= csp_psi_fit_tdem(1,2,j)
              else if (intg_tdem(8) .eq. 2) then
c                 case not implemented
                  call STOP('ICSVKU fit of tdem data not available in 
     .Onetwo',0)
              else if (intg_tdem(8) .eq. 3) then
                 dpsidt_const_zeta_tdem(j) =
     .                 ((3.0*csp_psi_fit_tdem(il_tdem,3,j) *delta_time +
     .                   2.0*csp_psi_fit_tdem(il_tdem,2,j))*delta_time +
     .                   csp_psi_fit_tdem(il_tdem,1,j)) ! volts/radian
              else
                call STOP ('subroutine UPDATE_TDEM_DATA: error', 260)
              end if

              slice = 'time'
c
              dpsidt_const_rho_tdem(j)=dpsidt_const_zeta_tdem(j)-
     .          dpsidrho_const_t_tdem(j)*(r(j)/rhoa)*
     .          drhomaxdt_tdem    ! volts/radian
c
****      end if
c

          if (igetcur .le. 0) then
c
              call rmepc_12 (iostat,eqdskin,'fcap',mxtbcmhd,
     .                       byte,aschar,intg2a,intg4a,real4a,ydum,
     .                       ishot,j,itype,nw_tdem,ncdfile_open,slice)

              interp_prof='fcap'
              call tdem_diffwrt_time (ydum, ntime_tdem, rtime_tdem,
     .                                t_interp, fcap(j), 1, dfdt(j),
     .                                dumy)
c
              call rmepc_12 (iostat,eqdskin,'gcap',mxtbcmhd,
     .                       byte,aschar,intg2a,intg4a,real4a,ydum,
     .                       ishot,j,itype,nw_tdem,ncdfile_open,slice)

              interp_prof='gcap'
              call tdem_diffwrt_time (ydum, ntime_tdem, rtime_tdem,
     .                                t_interp, gcap(j), 1, dgdt(j),
     .                                dumy)

c
              call rmepc_12 (iostat,eqdskin,'hcap',mxtbcmhd,
     .               byte,aschar,intg2a,intg4a,real4a,ydum,
     .               ishot,j,itype,nw_tdem,ncdfile_open,
     .               slice)
              interp_prof='hcap'
              call tdem_diffwrt_time (ydum, ntime_tdem, rtime_tdem,
     .                                t_interp, hcap(j), 1, dhdt(j),
     .                                dumy)
c
              call rmepc_12 (iostat,eqdskin,'r2cap',mxtbcmhd,
     .               byte,aschar,intg2a,intg4a,real4a,ydum,
     .               ishot,j,itype,nw_tdem,ncdfile_open,
     .               slice)   ! <R0**2/R**2>
              interp_prof='r2cap'
              call tdem_intrp(ydum,ntime_tdem,rtime_tdem,t_interp,
     .                                               r2cap(j))
              call tdem_diffwrt_time (ydum, ntime_tdem, rtime_tdem,
     .                                t_interp, r2cap(j), 1, dr2dt(j),
     .                                dumy)
c
              call rmepc_12 (iostat,eqdskin,'r2capi',mxtbcmhd,
     .               byte,aschar,intg2a,intg4a,real4a,ydum,
     .               ishot,j,itype,nw_tdem,ncdfile_open,
     .               slice)   ! <R**2>
              interp_prof='r2capi'
              call tdem_diffwrt_time (ydum, ntime_tdem, rtime_tdem,
     .                                t_interp, r2capi(j), 1, dr2idt(j),
     .                                dumy)
c
              call rmepc_12 (iostat,eqdskin,'rcap',mxtbcmhd,
     .               byte,aschar,intg2a,intg4a,real4a,ydum,
     .               ishot,j,itype,nw_tdem,ncdfile_open,
     .               slice)  ! <R>
              interp_prof='rcap'
              call tdem_diffwrt_time (ydum, ntime_tdem, rtime_tdem,
     .                                t_interp, rcap(j), 1, drcapdt(j),
     .                                dumy)


c             NOTE rcapi = <1/R> is required only by kinetic efit
c             subrotuine wrt_kin_efit_naml. But <1/R> is currently
c             not computed in the mepc code. hence we do not
c             call rmepc with rcapi. This means that  wrt_kin_efit_naml
c             will not function in tdem mode until this is remedied.
c             HSJ   8/27/04
c             UPDATE 1/16/2013 <1/R> is now in the netcdf tdem file
c         
c
              call rmepc_12 (iostat,eqdskin,'rcapi',mxtbcmhd,
     .               byte,aschar,intg2a,intg4a,real4a,ydum,
     .               ishot,j,itype,nw_tdem,ncdfile_open,
     .               slice)  ! <1/R> ?
c                 NOTE: ydum  now contains rcapi AT GRID POINT J as
c                 a funtion of time (not as a function of rho)
                  interp_prof='rcapi'
                  call tdem_diffwrt_time (ydum, ntime_tdem, rtime_tdem,
     .                          t_interp, rcapi(j), 1, drcapidt(j),
     .                          dumy)
c                 rcapi(j) is now value of rcapi at rho(j) at this time

c
              if (xi_include) then
c
c              time derivs not required for these
c              depsdt,...etc used only to predict eps,..etc
c              in mhdmethd .ne. tdem mode
c
                  call rmepc_12 (iostat,eqdskin,'eps',mxtbcmhd,
     .               byte,aschar,intg2a,intg4a,real4a,ydum,
     .               ishot,j,itype,nw_tdem,ncdfile_open,
     .               slice)
                  interp_prof='eps'
                  call tdem_intrp(ydum,ntime_tdem,rtime_tdem,t_interp,
     .                                                      eps(j))
c
                  call rmepc_12 (iostat,eqdskin,'xhm2',mxtbcmhd,
     .               byte,aschar,intg2a,intg4a,real4a,ydum,
     .               ishot,j,itype,nw_tdem,ncdfile_open,
     .               slice)
                  interp_prof='xhm2'
                  call tdem_intrp(ydum,ntime_tdem,rtime_tdem,t_interp,
     .                                                     xhm2(j))
c
                  call rmepc_12 (iostat,eqdskin,'xi11',mxtbcmhd,
     .               byte,aschar,intg2a,intg4a,real4a,ydum,
     .               ishot,j,itype,nw_tdem,ncdfile_open,
     .               slice)
                  interp_prof='xi11'
                  call tdem_intrp(ydum,ntime_tdem,rtime_tdem,t_interp,
     .                                                     xi11(j))
c
                  call rmepc_12 (iostat,eqdskin,'xi33',mxtbcmhd,
     .               byte,aschar,intg2a,intg4a,real4a,ydum,
     .               ishot,j,itype,nw_tdem,ncdfile_open,
     .               slice)
                  interp_prof='xi33'
                  call tdem_intrp(ydum,ntime_tdem,rtime_tdem,t_interp,
     .                                                     xi33(j))
c
                  call rmepc_12 (iostat,eqdskin,'xips',mxtbcmhd,
     .               byte,aschar,intg2a,intg4a,real4a,ydum,
     .               ishot,j,itype,nw_tdem,ncdfile_open,
     .               slice)
                  interp_prof='xips'
                  call tdem_intrp(ydum,ntime_tdem,rtime_tdem,t_interp,
     .                                                     xips(j))
             end if !xi_include
c
        end if    ! end igetcur .le. 0 branch
c

      end do      ! end time interpolation on nj grid

c
      do j=1,nj
        if (j .eq. 1) then
            dpsidrho_const_t_tdem(j)=
     .           (const_zeta_psigrid_smooth_tdem(j+1)-
     .            const_zeta_psigrid_smooth_tdem(j))/
     .                 (rhoa*(r(2)/r(nj)))
        else if ( j .eq. nj) then
           dpsidrho_const_t_tdem(j)=
     .           (const_zeta_psigrid_smooth_tdem(j)-
     .            const_zeta_psigrid_smooth_tdem(j-1))/
     .                 (rhoa*(r(nj)-r(nj-1))/r(nj))
        else
           dpsidrho_const_t_tdem(j)=
     .           (const_zeta_psigrid_smooth_tdem(j+1)-
     .            const_zeta_psigrid_smooth_tdem(j-1))/
     .                 (rhoa*(r(j+1)-r(j-1))/r(nj))
        end if
        dpsidt_const_rho_tdem(j)=dpsidt_const_zeta_tdem(j)-
     .                dpsidrho_const_t_tdem(j)*(r(j)/rhoa)*
     .                               drhomaxdt_tdem    ! volts/radian
c
      end do
c
       stop_timer(tdem_index)  =.TRUE.
       CALL collect_stats(tdem_index)
       stop_timer(tdem_index)  =.FALSE.
c      write(888,FMT='("stop tdem  l 2193",1pe14.6)')
c     .         elapsed_time(tdem_index) ! 88888889999999
      if (igetcur .eq. 1)  return
           start_timer(tdem_index) = .TRUE.
           stop_timer(tdem_index)  = .FALSE.
           CALL collect_stats(tdem_index)
           start_timer(tdem_index) = .FALSE. 
c      write(888,FMT='("start tdem 2199",1pe14.6)') 
c     .            elapsed_time(tdem_index) ! 88888889999999

c
c           some unit conversions:
c
            cconst = 1.0D4
            call multpl1(r2capi,nj,cconst)
            call multpl1(dr2idt,nj,cconst)
c
            cconst = 1.0D2
            call multpl1(rcap,nj,cconst)
            call multpl1(drcapdt,nj,cconst)

c

           cconst = 0.01D0
           call multpl1(rcapi,nj,cconst)   
           call multpl1(drcapidt,nj,cconst)

c
            interp_prof='rma'
            call tdem_intrp(rma_tdem,ntime_tdem,rtime_tdem,t_interp,rma)
            rma = rma * 100.0D0 ! convert to cm
            rma_t = rma
c
            interp_prof='zma'
            call tdem_intrp(zma_tdem,ntime_tdem,rtime_tdem,t_interp,zma)
            zma = zma * 100.0 ! convert to cm
            zma_t = zma
c
            interp_prof='psimag'
            call tdem_intrp(psimag_tdem,ntime_tdem,rtime_tdem,
     .                                           t_interp,psimag)
            psimag_t = psimag
c
****        call tdem_intrp(psilim_tdem,ntime_tdem,rtime_tdem,
**** .                                           t_interp,psilim)
c
c           psilim_tdem was read in once up front. Here we just
c           get the derivative on the (moving) plasma surface
c           as a function of time:
c
            interp_prof='psilim'
c            print *,'before tdem_diffwrt_time,psilim =',psilim
c            call tdem_diffwrt_time (psilim_tdem, ntime_tdem, rtime_tdem,
c     .                              t_interp, psilim, 1, dpsidt_tdem,
c     .                              dumy)
c            print *,'in update_tdem_data,spline psilim =',psilim
c            print *,'in update_tdem_data,spline dpsidt =',dpsidt_tdem
c            replace above spline inerpolation with tesnsion spline: HSJ 9/1/09
             ncd = 1  ! one continuous derivative at knots
             iendc = 0
             deriv_bc(:) =0.0 ; dv_interp(:) = 0.0
             nj_ncd = 1
            CALL t716_TSPSI (ntime_tdem,rtime_tdem,psilim_tdem,
     .                        ncd ,iendc,.FALSE.,.FALSE.,deriv_bc)
            iflag =    0          ! compute value at time t_interp
            v_interp(:) = t_interp
            CALL  t716_TSVAL1 (ntime_tdem,rtime_tdem,psilim_tdem,    
     .                            IFLAG,nj_ncd,v_interp,dv_interp)
            psilim =  dv_interp(1)
            iflag  =   1          ! compute derivative at time t_interp
            CALL  t716_TSVAL1 (ntime_tdem,rtime_tdem,psilim_tdem ,   
     .                            IFLAG,nj_ncd,v_interp,dv_interp)
            dpsidt_tdem = dv_interp(1)
            psilim_t = psilim

c
            interp_prof='beqd'
            call tdem_intrp(beqd_tdem,ntime_tdem,rtime_tdem,
     .                                           t_interp,beqd)
            btor=beqd*1.e4
            interp_prof='toteqd'
            call tdem_intrp(toteqd_tdem,ntime_tdem,rtime_tdem,
     .                                           t_interp,toteqd)
            interp_prof='psimx1'
            call tdem_intrp(psisep_tdem,ntime_tdem,rtime_tdem,
     .                                           t_interp,psimx1)
c
            interp_prof='rsep'
            call tdem_intrp(rsep_tdem,ntime_tdem,rtime_tdem,
     .                                           t_interp,rsep)
            rsep=rsep*100.0
            interp_prof='zsep'
            call tdem_intrp(zsep_tdem,ntime_tdem,rtime_tdem,
     .                                           t_interp,zsep)
            zsep=zsep*100.0
c
            interp_prof='tdem_vol'
            call tdem_intrp(vol_tdem,ntime_tdem,rtime_tdem,
     .                                           t_interp,tdem_vol)

c
c          get the volume of each flux zone on the psigrid at time
c          t_interp by interpolation in time:
c
           slice = 'time'

           do j=1,nw
c
   
              call rmepc_12 (iostat,eqdskin,'psigrid',mxtbcmhd,
     .               byte,aschar,intg2a,intg4a,real4a,ydum,
     .               ishot,j,itype,ntime_tdem,ncdfile_open,
     .               slice)    ! ydum(1,..ntime_tdem) contains
c                                psi value j at
c                                all the available time points
c
              interp_prof='psigrid'
              call tdem_intrp(ydum,ntime_tdem,rtime_tdem,t_interp,
     .                                             psigrid_tdem(j))
c
              call rmepc_12 (iostat,eqdskin,'volpsi',mxtbcmhd,
     .                       byte,aschar,intg2a,intg4a,real4a,ydum,
     .                       ishot,j,itype,ntime_tdem,ncdfile_open,
     .                       slice) ! ydum(1,..ntime_tdem) contains..
c                                   ..volume inside flux zone j at..
c                                   ..all the available time points
              interp_prof = 'volpsi'
              call tdem_intrp(ydum,ntime_tdem,rtime_tdem,t_interp,
     .                                             volpsi_tdem(j))
c             volpsi_tdem(1,..nw) now contains the volumes inside
c                    the flux zones at time t_interp. Note that
c                    volpsi(1) is plasma edge and volpsi(nw)(=0.0)
c                    is magnetic axis. The values of psi associated
c                    with these volumes are in psigrid_tdem
c
c

c-----------------------------------------------------------------------

            call rmepc_12 (iostat,eqdskin,'fpsi',mxtbcmhd,
     .               byte,aschar,intg2a,intg4a,real4a,ydum,
     .               ishot,j,itype,ntime_tdem,ncdfile_open,
     .               slice)

c           given ydum(1:ntime_tdem) on the time grid
c           rtime_tdem(1:ntime_tdem) get the value of each profile
c           at time t_interp  for the nw values that are normally 
c           on the eqdsks:

            interp_prof = 'fpsi'
            call tdem_intrp(ydum,ntime_tdem,rtime_tdem,t_interp,
     .                                            fpsi_tdem(j))
c
            call rmepc_12 (iostat,eqdskin,'presspsi',mxtbcmhd,
     .               byte,aschar,intg2a,intg4a,real4a,ydum,
     .               ishot,j,itype,nw_tdem,ncdfile_open,
     .               slice)
            interp_prof = 'presspsi'
            call tdem_intrp(ydum,ntime_tdem,rtime_tdem,t_interp,
     .                                          presspsi_tdem(j))
c
           call rmepc_12 (iostat,eqdskin,'ffppsi',mxtbcmhd,
     .               byte,aschar,intg2a,intg4a,real4a,ydum,
     .               ishot,j,itype,nw_tdem,ncdfile_open,
     .               slice)
            interp_prof = 'ffppsi'
            call tdem_intrp(ydum,ntime_tdem,rtime_tdem,t_interp,
     .                                          ffprimpsi_tdem(j))
c
           call rmepc_12 (iostat,eqdskin,'pppsi',mxtbcmhd,
     .               byte,aschar,intg2a,intg4a,real4a,ydum,
     .               ishot,j,itype,nw_tdem,ncdfile_open,
     .               slice)
            interp_prof = 'pppsi'
            call tdem_intrp(ydum,ntime_tdem,rtime_tdem,t_interp,
     .                                          pprimpsi_tdem(j))
c
           call rmepc_12 (iostat,eqdskin,'qpsi',mxtbcmhd,
     .               byte,aschar,intg2a,intg4a,real4a,ydum,
     .               ishot,j,itype,nw_tdem,ncdfile_open,
     .               slice)
          interp_prof = 'qpsi'
          call tdem_intrp(ydum,ntime_tdem,rtime_tdem,t_interp,
     .                                          qpsi_tdem(j))
c

c-----------------------------------------------------------------------
           end do 


c     for output to statefile
      IF(.NOT. ALLOCATED( pppsi_eqdsk))ALLOCATE(pppsi_eqdsk(nw))
      IF(.NOT. ALLOCATED( presspsi_eqdsk))ALLOCATE(presspsi_eqdsk(nw))
      IF(.NOT. ALLOCATED( fpsi_eqdsk))ALLOCATE(fpsi_eqdsk(nw))
      IF(.NOT. ALLOCATED( ffppsi_eqdsk))ALLOCATE(ffppsi_eqdsk(nw))
      IF(.NOT. ALLOCATED( qpsi_eqdsk))ALLOCATE(qpsi_eqdsk(nw))
      pppsi_eqdsk(1:nw)  = pprimpsi_tdem(1:nw)
!      CALL reverse(nw,pppsi_eqdsk)
      presspsi_eqdsk(1:nw)  = presspsi_tdem(1:nw)
!      CALL reverse(nw,presspsi_eqdsk)
      fpsi_eqdsk(1:nw)   = fpsi_tdem(1:nw)
!      CALL reverse(nw,fpsi_eqdsk)
      ffppsi_eqdsk(1:nw) = ffprimpsi_tdem(1:nw)
!      CALL reverse(nw,ffppsi_eqdsk)
      qpsi_eqdsk(1:nw)   = qpsi_tdem(1:nw)
!      CALL reverse(nw,qpsi_eqdsk)
      nxeqd_eqdsk           = nw
      psimag_eqdsk          = psimag
      psilim_eqdsk          = psilim



           cconst = 1.0e6
           call multpl1(volpsi_tdem,nw,cconst)     ! cm^3
           call multpl1(psigrid_tdem,nw,psikgaus)  ! in kgauss cm^2
c
           psisep = psimx1
           slice  = 'none'
c
c          do quick linear interp on psi
c
           call rmepc_12 (iostat,eqdskin,'psi',nw,
     .               byte,aschar,intg2a,intg4a,real4a,psi_tdem,
     .               ishot,il_tdem,itype,nw_tdem,ncdfile_open,
     .               slice)    ! at time point il_tdem
           call rmepc_12 (iostat,eqdskin,'psi',nw,
     .               byte,aschar,intg2a,intg4a,real4a,psi,
     .               ishot,il_tdem+1,itype,nw_tdem,ncdfile_open,
     .               slice) ! at time point il_tdem+1, use p temporarily
           delt = t_interp-rtime_tdem(il_tdem)
           delt = delt/(rtime_tdem(il_tdem+1)-rtime_tdem(il_tdem))
           do j=1,nh 
             do i=1,nw
                p(i,j)=psi_tdem(i,j)+(psi(i,j)-psi_tdem(i,j))*delt
                p(i,j)=p(i,j)*psikgaus       ! convert p(i,j) to kguass cm^2
             end do
           end do
c
c         locate magnetic axis of interpolated equilibrium:
c
          isignn  = -1  ! search for minimum
          iknowax =  1  ! good guess for axis exists
          ispln   =  0  ! bicubic spline cspln must be calculated
          ndwk    =  2*nw*nh+2*nh
c
          call magax (p,rmhdgrid,zmhdgrid,nw,nh,isignn,iknowax,
     .                ncrt,ispln,zero,rma,zma,cspln,n2cspln,nh2,
     .                psimag,elax)
          psiaxis = psimag  ! in kgauss cm2 due to p
          psibdry = psilim * psikgaus   ! new value corresponding to p..
c                             ..not calculated, may have to do in future
          xmagn1  = rma
          xax(1)  = rma
          ymagn1  = zma     ! should come up with a couple more names..
          yax(1)  = zma     ! ..here, just to make it really confusing
          pmin    = psiaxis ! store in small.i, use in subroutine PSISET
          pmax    = psibdry
          call psiset (npsi, psiaxis, psibdry, psival, 1)

c     INPUT to intrp : psigrid_tdem, volpsi_tdem, nw

c     OUTPUTof intrp : psival, psivolp, npsi
      call intrp (1, 1, psigrid_tdem, volpsi_tdem, nw, 
     .            psival, psivolp, npsi)


      interp_prof = '????'
       stop_timer(tdem_index)  =.TRUE.
       CALL collect_stats(tdem_index)
       stop_timer(tdem_index)  =.FALSE.

c       write(888,FMT='("stop tdem  2462"1pe14.6)')
c     .         elapsed_time(tdem_index)  ! 88888889999999
 
      return
c
      end
