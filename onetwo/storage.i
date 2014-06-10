c
c --- INCLUDE file storage.i
c
      integer   kstore,kstorest
      real *8 xdum,ydum,zdum,wdum,vdum,sdum,tdum,udum,rdum
c      parameter (kstore   =  2000)                ! must be .ge. nconmax
      parameter (kstore   =  6000)                ! must be .ge. nconmax
      parameter (kstorest = 10000)
c
c --- the following arrays are volatile. that is, they may be used by
c --- any program segment but are not for permanent storage outside that
c --- segment. (see also INCLUDE file bicube.i for additional arrays.)
c
      common /storage/ xdum(kstore  ), ydum(kstore  ), zdum(kstore  ),
     .                 wdum(kstore  ), vdum(kstore  ), sdum(kstorest),
     .                 tdum(kstorest), udum(kstore  ), rdum(kstore  )
c
