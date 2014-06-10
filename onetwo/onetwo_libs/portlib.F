      real *8  function AMAXAF (array, first, last)
c
c --- replacement for this LIBMATH function in UNICOS ------------------
c
      implicit none

      character what_id*45
      save      what_id
      data      what_id/"@(#)portlib.f      17 Nov 98  General Atomics"/
c
      integer   first, last, index
      real*8    array(last), maximum
c
      maximum = array(first)
c
      do index=first,last
        if (array(index) .gt. maximum)  maximum = array(index)
      end do
c
      AMAXAF = maximum
      return
c
      end

      subroutine AMINMX (a, i1, i2, i3, a1, a2)
c
c --- replacement for this LIBMATH subroutine in UNICOS ----------------
c
      implicit none
c
      integer  i1, i2, i3, i, i0
      real*8   a(*), a1, a2
c
      a1 = a(i1)
      a2 = a(i1)
      i0 = i1 + i3
      do i=i0,i2,i3
        a1 = MIN (a(i), a1)
        a2 = MAX (a(i), a2)
      end do
      return
c
      end

      subroutine COPY_TEXT_FILE (oldfile, newfile)
c
c --- copy one 80-column-wide text file to another ---------------------
c
      implicit none
c
      real     * 8    banktime
      logical                   trace
      character*32                     routine
      common /portee/ banktime, trace, routine
c
      character       oldfile*(*), newfile*(*), record*80
      integer         oldunit    , newunit
      parameter      (oldunit = 1, newunit = 2)
c
      routine = 'COPY_TEXT_FILE'
      if (trace)  write (*, '( / '' ---- entering '', a   )')  routine
c
      open    (unit = oldunit, file = oldfile, status = 'OLD'    )
      open    (unit = newunit, file = newfile, status = 'UNKNOWN')
      do while (.true.)
        read  (unit = oldunit, fmt = '(a)', end = 10) record
        write (unit = newunit, fmt = '(a)'          ) record
      end do
   10 close   (unit = oldunit)
      close   (unit = newunit)
c
      if (trace)  write (*, '(   '' ----  exiting '', a / )')  routine
      return
c
      end

      subroutine DESTROY (file)
c
c --- remove the specified file ----------------------------------------
c
      implicit none
c
      real     * 8    banktime
      logical                   trace
      character*32                     routine
      common /portee/ banktime, trace, routine
c
      integer         status, UNLINK
      character*(*)   file
c
      routine = 'DESTROY'
      if (trace)  write (*, '( / '' ---- entering '', a   )')  routine
      status  =  UNLINK (file)
      if (trace)  write (*, '(   '' ----  exiting '', a / )')  routine
      return
c
      end

      real*8 function D1MACH (item)
c
c --- return machine-dependent floating point constants ----------------
c
      implicit none
c
      integer  item
      real*8   rmach(5)
      save     rmach
c
      data     rmach(1) / z'0010000000000000' /,
     .         rmach(2) / z'7fefffffffffffff' /,
     .         rmach(3) / z'3ca0000000000000' /,
     .         rmach(4) / z'3cb0000000000000' /,
     .         rmach(5) / z'3fd34413509f79ff' /
c
      D1MACH = rmach(item)
      return
c
      end

      subroutine EMPTY (unit)
c
c --- flush the I/O buffer of a file (if unit is open) -----------------
c
      implicit none
c
      real     * 8    banktime
      logical                   trace
      character*32                     routine
      common /portee/ banktime, trace, routine
c
      logical         opened
      integer         unit
      character*128 name
c
      routine = 'EMPTY'
      if (trace)  write (*, '( / '' ---- entering '', a   )')  routine
c
c --- check to see if the file is open
c
      inquire (unit = unit, opened = opened, name = name)

c
c --- flush the buffer by closing and then reopening the file at its end
c
      if (opened) then
        close (unit = unit,              status = 'KEEP')
       open  (UNIT = unit, FILE = name, STATUS = 'OLD',
     .         POSITION  = 'APPEND')
      else
        write  (*, 10)  routine, unit
   10   format (/ ' ---- WARNING from routine ', a         /
     .            '      unit number ', i2, ' is not open' /)
      end if
c
      if (trace)  write (*, '(   '' ----  exiting '', a / )')  routine
      return
c
      end

      subroutine EMPTY_AND_RECORD_POINTER (unit, pointer)
c
c --- flush I/O buffer of a file (if unit is open) and record pointer --
c
      implicit none
c
      real     * 8    banktime
      logical                   trace
      character*32                     routine
      common /portee/ banktime, trace, routine
c
      logical         opened
      integer         unit, pointer, iostat, LENGTH
      character       name*128, record*133
      external                               LENGTH
c
      routine = 'EMPTY_AND_RECORD_POINTER'
      if (trace)  write (*, '( / '' ---- entering '', a   )')  routine
c
c --- check to see if the file is open
c
      inquire ( unit =  unit, opened = opened, name = name)
c
c --- flush buffer by closing and reopening file, then going to the end
c
      if (opened) then
        close (unit = unit,              status = 'KEEP')
        open  (unit = unit, file = name, status = 'OLD' )
        pointer = 0
        read (unit, '(a)', iostat = iostat) record
        do while (iostat .eq. 0)
          pointer = pointer + LENGTH (record)
          read (unit, '(a)', iostat = iostat) record
        end do
        backspace unit
      else
        write  (*, 10)  routine, unit
   10   format (/ ' ---- WARNING from routine ', a         /
     .            '      unit number ', i2, ' is not open' /)
      end if
c
      if (trace)  write (*, '(   '' ----  exiting '', a / )')  routine
      return
c
      end

      subroutine FLUSH (unit, status)
c
c --- flush the I/O buffer of a file (if unit is open) -----------------
c
      implicit none
c
      real     * 8    banktime
      logical                   trace
      character*32                     routine
      common /portee/ banktime, trace, routine
c
      integer         unit, status
c
      routine = 'FLUSH'
      if (trace)  write (*, '( / '' ---- entering '', a   )')  routine
c
c --- no functionality yet; return success status
c
      status = 0
c
      if (trace)  write (*, '(   '' ----  exiting '', a / )')  routine
      return
c
      end

      subroutine GET_DATE_TIME (result)
c
c --- get the current date and time as a 20-character ASCII string -----
c
      implicit none
c
      integer    TIME
      character CTIME*24, result*20, buffer*24
c
      buffer =  CTIME (TIME ())
      result =  buffer(5:24)
      return
c
      end

      subroutine GET_PROGRAM_INFO (name, owner,group,size, date, status)
c
c --- get some information about the current program file --------------
c
      implicit none
c
      real     * 8    banktime
      logical                   trace
      character*32                     routine
      common /portee/ banktime, trace, routine
c
      character       name*(*),     CTIME*24, date*24
      external        STAT, GETARG, CTIME
      integer         STAT, stat_block(12), owner, group, size, status
c
      routine = 'GET_PROGRAM_INFO'
      if (trace)  write (*, '( / '' ---- entering '', a   )')  routine
c
      call GETARG (0, name)
      status  = STAT (name, stat_block)
      owner   =             stat_block( 5)
      group   =             stat_block( 6)
      size    =             stat_block( 8)
      date    =      CTIME (stat_block(10))
c
      if (trace)  write (*, '(   '' ----  exiting '', a / )')  routine
      return
c
      end

      integer function ISHELL (command)
c
c --- replacement for this function in UNICOS --------------------------
c
      implicit none
c
      integer       SYSTEM
      character*(*) command
c
      if (SYSTEM (command) .eq. 0) then
        ISHELL =  0
      else
        ISHELL = -1
      end if
      return
c
      end

      integer function LENGTH (string)
c
c --- return trimmed (no trailing spaces or tabs) length of a string ---
c
      implicit none
c
      character  last*1, string*(*)
c
      LENGTH = LEN (string)
      if (LENGTH .ne. 0) then
        last = string(LENGTH:LENGTH)
        do while (last .eq. ' ' .or. last .eq. '	')
          LENGTH = LENGTH - 1
          if (LENGTH .ne. 0) then
            last = string(LENGTH:LENGTH)
          else
            last = '#'
          end if
        end do
      end if
      return
c
      end

      subroutine LINKETTE (unit1, file1, unit2, file2,
     .                     unit3, file3, unit4, file4, options)
c
c --- special initialization routine for use only with ONETWO ----------
c
      implicit none
c
      real     * 8    banktime
      logical                   trace
      character*32                     routine
      common /portee/ banktime, trace, routine
c
      integer         LENGTH,
     .                unit1, unit2, unit3, unit4, iarg, nchars
      real     * 8    bankmins
      character*(*)   file1, file2, file3, file4
      character* 2    flag1, flag2, flag3, flag4, flag5
      character* 6    flag_t1, flag_t2
      character*64    argument, options
      data            flag1/'i='/,flag2/'o='/,flag3/'q='/,flag4/'y='/,
     .                flag5/'t='/,flag_t1/'-trace'/,flag_t2/'/trace'/,
     .                trace/.false./, bankmins/1.0/
      external        LENGTH, GETARG
c
      routine = 'LINKETTE'
      if (trace)  write (*, '( / '' ---- entering '', a   )')  routine
c
c --- load the default filenames and options
c
      file1   = 'inone'
      file2   = 'outone'
      file3   = 'qikone'
      file4   = 'yokfil'
      options = 'reserved for future command-line options'
c
c --- get all the command-line arguments
c
#ifndef PARALLEL_12 
      iarg   = 1
      call GETARG (iarg, argument)
      nchars =   LENGTH (argument)
      do while (nchars .gt. 0)
        if      (INDEX (argument, flag1) .eq. 1) then
          file1 = argument(3:)
        else if (INDEX (argument, flag2) .eq. 1) then
          file2 = argument(3:)
        else if (INDEX (argument, flag3) .eq. 1) then
          file3 = argument(3:)
        else if (INDEX (argument, flag4) .eq. 1) then
          file4 = argument(3:)
        else if (INDEX (argument, flag5) .eq. 1) then
          read  (argument(3:), '(f8.0)') bankmins
        else if (       argument .eq. flag_t1  ) then
          trace = .true.
        else if (       argument .eq. flag_t2  ) then
          trace = .true.
        else
          write  (*, 10)  routine, argument(1:LENGTH (argument))
   10     format (/ ' ---- WARNING from routine ', a /
     .              '      command-line argument "', a,
     .                  '" was not recognized; it will be ignored' /)
        end if
        iarg = iarg + 1
        if (iarg .le. 8) then
          call GETARG (iarg, argument)
          nchars =   LENGTH (argument)
        else
          nchars = -1
        end if
      end do
#endif
c
      banktime = 60.0 * bankmins
c
      if (trace)  write (*, '(   '' ----  exiting '', a / )')  routine
      return
c
      end

      integer function MAXAF (array, first, last)
c
c --- replacement for this LIBMATH function in UNICOS ------------------
c
      implicit none
c
      integer   first, last, index
      integer   array(last), maximum
c
      maximum = array(first)
c
      do index=first,last
        if (array(index) .gt. maximum)  maximum = array(index)
      end do
c
      MAXAF = maximum
      return
c
      end

      subroutine MYFLUSH (unit, filename)
c
c --- flush an I/O buffer ----------------------------------------------
c
      implicit none
c
      integer       unit
      character*(*) filename
c
      close (unit = unit,                  status = 'KEEP')
      open  (unit = unit, file = filename, status = 'OLD',
cjmp.ibm     .       access = 'APPEND')
     .       position = 'APPEND') !jmp.ibm
      return
c
      end

      subroutine OPEN_APPEND (unit, file)
c
c --- open an existing sequential formatted file for append writing ----
c
      implicit none
c
      real     * 8    banktime
      logical                   trace
      character*32                     routine
      common /portee/ banktime, trace, routine
c
      integer         unit
      character*(*)   file
c
      routine = 'OPEN_APPEND'
      if (trace)  write (*, '( / '' ---- entering '', a   )')  routine
c
c --- open the file
c
      open (unit = unit, file = file, status = 'OLD',
cjmp.ibm     .      access = 'APPEND')
     .      position = 'APPEND') !jmp.ibm
c
      if (trace)  write (*, '(   '' ----  exiting '', a / )')  routine
      return
c
      end

      real*8 function PIMAG (complex_number)
c
c --- portable replacement for non-portable Fortran intrinsic IMAG -----
c
      implicit none
c
      complex*16 complex_number
c
      PIMAG = DIMAG (complex_number)
      return
c
      end

      real*8 function RANDOM12 (seed)
c
c --- basic uniform (0,1) pseudo-random number generator (portable) ----
c
      implicit none
c
      save         saved_seed
      real*8 seed, saved_seed,
     .       big,                  big_minus_1
      data   big /2147483648.0d0/, big_minus_1 /2147483647.0d0/
c
      if (seed .ne. 0.0d0)
     .saved_seed = seed
      saved_seed = MOD (16807.0d0*saved_seed, big_minus_1)
      RANDOM12     = saved_seed / big
      return
c
      end

      real*4 function R1MACH (item)
c
c --- return machine-dependent floating point constants ----------------
c
      implicit none
c
      integer  item
      real*4   rmach(5)
      save     rmach
c
      data     rmach(1) / z'00800000' /,
     .         rmach(2) / z'7f7fffff' /,
     .         rmach(3) / z'33800000' /,
     .         rmach(4) / z'34000000' /,
     .         rmach(5) / z'3e9a209b' /
c
      R1MACH = rmach(item)
      return
c
      end

      subroutine SECOND (seconds)
c
c --- replacement for this subroutine in UNICOS ------------------------
c
      implicit none
c
      real*8  seconds
      real*4  ETIME, times(2)
c
c --- return elapsed CPU time
c
cjmp.ibm     seconds = DBLE (ETIME (times))
cjmp.ibm     call cpu_time_12(seconds) !jmp.ibm      
      call cpu_time(seconds)
      return
c
      end

cjmp.ibm_start
#ifdef PARALLEL_12
      subroutine STOP (message, status)
c
c --- print specified message and then exit with the specified status --
c
      use mpi12
      implicit none
c
      integer   status,ierr
      character message*(*)
      INTEGER iflag_mpistop(1) 

      if (status .eq. 0) then
        write (*, '(/ '' STOP in '', a /)')
     .  message
      endif

      iflag_mpistop(1) = 1
      CALL mpi12_send_int(iflag_mpistop,1,99)
      
      CALL MPI_FINALIZE(ierr)
      
      call EXIT (-1)
      
      end
cjmp.ibm_end
#else
      subroutine STOP (message, status)
c
c --- print specified message and then exit with the specified status --
c
      implicit none
c
      integer   ISHELL, status, position
      character message*(*), program_name*128, command*256
      external  GETARG
c
      if (status .eq. 0) then
        write (*, '(/ '' STOP in '', a /)')
     .  message
      else
        call GETARG (0, program_name)
        position = INDEX (program_name, 'onetwo')
        if (position .ne. 0) then
          position = INDEX (program_name(position:), '/')
          print *,'position =',position
          if (position .eq. 0) then
            command = 'rm -f bpltfil eqpltfil ; '              //
     .                'touch outone qikone runlog '            //
     .                      'trpltfil trpltout trpltout.nc ; ' //
     .                'mv outone _outone ; '                   //
     .                'mv qikone _qikone ; '                   //
     .                'mv runlog _runlog ; '                   //
     .                'mv trpltfil _trpltfil ; '               //
     .                'mv trpltout _trpltout ; '               //
     .                'mv trpltout.nc _trpltout.nc ; '         //
     .                'echo "\n OUTPUT IS INVALID"'
            if (ISHELL (command) .ne. 0) then
              stop 'Last-chance abort in subroutine STOP'
            end if
          end if
        end if
        write (*, '(/ '' STOP in '', a / '' EXIT status is'', i4 /)')
     .  message, status
      end if
      call EXIT (status)
c
      end
#endif

      subroutine TIMELEFT (iseconds)
c
c --- get the time remaining -------------------------------------------
c
      implicit none
c
      real     * 8    banktime
      logical                   trace
      character*32                     routine
      common /portee/ banktime, trace, routine
c
      integer         iseconds
      real*8          seconds
c
      routine = 'TIMELEFT'
      if (trace)  write (*, '( / '' ---- entering '', a   )')  routine
c
      call SECOND (seconds)
      iseconds = banktime - seconds
c
      if (trace)  write (*, '(   '' ----  exiting '', a / )')  routine
      return
c
      end

      subroutine TIMEUSED (cpu)
c
c --- get the time used ------------------------------------------------
c
      implicit none
c
      real     * 8    banktime
      logical                   trace
      character*32                     routine
      common /portee/ banktime, trace, routine
c
      integer         cpu
      real*8          seconds
c
      routine = 'TIMEUSED'
      if (trace)  write (*, '( / '' ---- entering '', a   )')  routine
c
      call SECOND (seconds)
      cpu = seconds
c
      if (trace)  write (*, '(   '' ----  exiting '', a / )')  routine
      return
c
      end
      subroutine freeus(unitno)
c
c     $Id: portlib.F,v 1.2 2007/04/09 20:57:02 parkjm Exp $
c
c dummy replacement for Basis's freeus, which returns a free unit-
c number in its argument.
      return
      end

      subroutine giveupus(unitno)
c dummy replacement for Basis's giveupus, which gives up io unit-no. unitno.
      return
      end

      subroutine kaboom(arg)
c dummy replacement for Basis's kaboom, which returns to the controller
c and clears the stack.
      return
      end

      subroutine ruthere
c dummy replacement for Basis's ruthere, which checks for and handles ctrl-c.
      return
      end

      subroutine execuser(f)
c dummy replacement for Basis's execuser, which executes script-function f.
      character*(*) f
      return
      end

      subroutine parsestr(str)
c dummy replacement for Basis's parsestr, which (recursively) calls the
c parser with command-string str.
      character*(*) str
      return
      end
