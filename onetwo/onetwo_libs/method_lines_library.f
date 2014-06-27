c----------------------start of dlsoibt related routines ----------------------

*DECK DLSOIBT
      SUBROUTINE DLSOIBT (RES, ADDA, JAC, NEQ, Y, YDOTI, T, TOUT, ITOL,
     1  RTOL, ATOL, ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, LIW, MF )
      EXTERNAL RES, ADDA, JAC
      INTEGER NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, IWORK, LIW, MF
      DOUBLE PRECISION Y, YDOTI, T, TOUT, RTOL, ATOL, RWORK
      DIMENSION NEQ(*), Y(*), YDOTI(*), RTOL(*), ATOL(*), RWORK(LRW),
     1          IWORK(LIW)
C-----------------------------------------------------------------------
C This is the 7 May 2001 version of
C DLSOIBT: Livermore Solver for Ordinary differential equations given
C          in Implicit form, with Block-Tridiagonal Jacobian treatment.
C
C This version is in double precision.
C
C DLSOIBT solves the initial value problem for linearly implicit
C systems of first order ODEs,
C     A(t,y) * dy/dt = g(t,y) ,  where A(t,y) is a square matrix,
C or, in component form,
C     ( a   * ( dy / dt ))  + ... +  ( a     * ( dy   / dt ))  =
C        i,1      1                     i,NEQ      NEQ
C
C      =   g ( t, y , y ,..., y    )   ( i = 1,...,NEQ )
C           i      1   2       NEQ
C
C If A is singular, this is a differential-algebraic system.
C
C DLSOIBT is a variant version of the DLSODI package, for the case where
C the matrices A, dg/dy, and d(A*s)/dy are all block-tridiagonal.
C-----------------------------------------------------------------------
C Reference:
C     Alan C. Hindmarsh,  ODEPACK, A Systematized Collection of ODE
C     Solvers, in Scientific Computing,  R. S. Stepleman et al. (Eds.),
C     North-Holland, Amsterdam, 1983, pp. 55-64.
C-----------------------------------------------------------------------
C Authors:       Alan C. Hindmarsh and Jeffrey F. Painter
C                Center for Applied Scientific Computing, L-561
C                Lawrence Livermore National Laboratory
C                Livermore, CA 94551
C and
C                Charles S. Kenney
C formerly at:   Naval Weapons Center
C                China Lake, CA 93555
C-----------------------------------------------------------------------
C Summary of Usage.
C
C Communication between the user and the DLSOIBT package, for normal
C situations, is summarized here.  This summary describes only a subset
C of the full set of options available.  See the full description for
C details, including optional communication, nonstandard options,
C and instructions for special situations.  See also the example
C problem (with program and output) following this summary.
C
C A. First, provide a subroutine of the form:
C               SUBROUTINE RES (NEQ, T, Y, S, R, IRES)
C               DOUBLE PRECISION T, Y(*), S(*), R(*)
C which computes the residual function
C     r = g(t,y)  -  A(t,y) * s ,
C as a function of t and the vectors y and s.  (s is an internally
C generated approximation to dy/dt.)  The arrays Y and S are inputs
C to the RES routine and should not be altered.  The residual
C vector is to be stored in the array R.  The argument IRES should be
C ignored for casual use of DLSOIBT.  (For uses of IRES, see the
C paragraph on RES in the full description below.)
C
C B. Next, identify the block structure of the matrices A = A(t,y) and
C dr/dy.  DLSOIBT must deal internally with a linear combination, P, of
C these two matrices.  The matrix P (hence both A and dr/dy) must have
C a block-tridiagonal form with fixed structure parameters
C     MB = block size, MB .ge. 1, and
C     NB = number of blocks in each direction, NB .ge. 4,
C with MB*NB = NEQ.  In each of the NB block-rows of the matrix P
C (each consisting of MB consecutive rows), the nonzero elements are
C to lie in three consecutive MB by MB blocks.  In block-rows
C 2 through NB - 1, these are centered about the main diagonal.
C in block-rows 1 and NB, they are the diagonal blocks and the two
C blocks adjacent to the diagonal block.  (Thus block positions (1,3)
C and (NB,NB-2) can be nonzero.)
C Alternatively, P (hence A and dr/dy) may be only approximately
C equal to matrices with this form, and DLSOIBT should still succeed.
C The block-tridiagonal matrix P is described by three arrays,
C each of size MB by MB by NB:
C     PA = array of diagonal blocks,
C     PB = array of superdiagonal (and one subdiagonal) blocks, and
C     PC = array of subdiagonal (and one superdiagonal) blocks.
C Specifically, the three MB by MB blocks in the k-th block-row of P
C are stored in (reading across):
C     PC(*,*,k) = block to the left of the diagonal block,
C     PA(*,*,k) = diagonal block, and
C     PB(*,*,k) = block to the right of the diagonal block,
C except for k = 1, where the three blocks (reading across) are
C     PA(*,*,1) (= diagonal block), PB(*,*,1), and PC(*,*,1),
C and k = NB, where they are
C     PB(*,*,NB), PC(*,*,NB), and PA(*,*,NB) (= diagonal block).
C (Each asterisk * stands for an index that ranges from 1 to MB.)
C
C C. You must also provide a subroutine of the form:
C     SUBROUTINE ADDA (NEQ, T, Y, MB, NB, PA, PB, PC)
C     DOUBLE PRECISION T, Y(*), PA(MB,MB,NB), PB(MB,MB,NB), PC(MB,MB,NB)
C which adds the nonzero blocks of the matrix A = A(t,y) to the
C contents of the arrays PA, PB, and PC, following the structure
C description in Paragraph B above.
C T and the Y array are input and should not be altered.
C Thus the affect of ADDA should be the following:
C     DO 30 K = 1,NB
C       DO 20 J = 1,MB
C         DO 10 I = 1,MB
C           PA(I,J,K) = PA(I,J,K) +
C             ( (I,J) element of K-th diagonal block of A)
C           PB(I,J,K) = PB(I,J,K) +
C             ( (I,J) element of block in block position (K,K+1) of A,
C             or in block position (NB,NB-2) if K = NB)
C           PC(I,J,K) = PC(I,J,K) +
C             ( (I,J) element of block in block position (K,K-1) of A,
C             or in block position (1,3) if K = 1)
C 10        CONTINUE
C 20      CONTINUE
C 30    CONTINUE
C
C D. For the sake of efficiency, you are encouraged to supply the
C Jacobian matrix dr/dy in closed form, where r = g(t,y) - A(t,y)*s
C (s = a fixed vector) as above.  If dr/dy is being supplied,
C use MF = 21, and provide a subroutine of the form:
C     SUBROUTINE JAC (NEQ, T, Y, S, MB, NB, PA, PB, PC)
C     DOUBLE PRECISION T, Y(*), S(*), PA(MB,MB,NB), PB(MB,MB,NB),
C    1                 PC(MB,MB,NB)
C which computes dr/dy as a function of t, y, and s.  Here T, Y, and
C S are inputs, and the routine is to load dr/dy into PA, PB, PC,
C according to the structure description in Paragraph B above.
C That is, load the diagonal blocks into PA, the superdiagonal blocks
C (and block (NB,NB-2) ) into PB, and the subdiagonal blocks (and
C block (1,3) ) into PC.  The blocks in block-row k of dr/dy are to
C be loaded into PA(*,*,k), PB(*,*,k), and PC(*,*,k).
C     Only nonzero elements need be loaded, and the indexing
C of PA, PB, and PC is the same as in the ADDA routine.
C     Note that if A is independent of Y (or this dependence
C is weak enough to be ignored) then JAC is to compute dg/dy.
C     If it is not feasible to provide a JAC routine, use
C MF = 22, and DLSOIBT will compute an approximate Jacobian
C internally by difference quotients.
C
C E. Next decide whether or not to provide the initial value of the
C derivative vector dy/dt.  If the initial value of A(t,y) is
C nonsingular (and not too ill-conditioned), you may let DLSOIBT compute
C this vector (ISTATE = 0).  (DLSOIBT will solve the system A*s = g for
C s, with initial values of A and g.)  If A(t,y) is initially
C singular, then the system is a differential-algebraic system, and
C you must make use of the particular form of the system to compute the
C initial values of y and dy/dt.  In that case, use ISTATE = 1 and
C load the initial value of dy/dt into the array YDOTI.
C The input array YDOTI and the initial Y array must be consistent with
C the equations A*dy/dt = g.  This implies that the initial residual
C r = g(t,y) - A(t,y)*YDOTI  must be approximately zero.
C
C F. Write a main program which calls Subroutine DLSOIBT once for
C each point at which answers are desired.  This should also provide
C for possible use of logical unit 6 for output of error messages by
C DLSOIBT.  on the first call to DLSOIBT, supply arguments as follows:
C RES    = name of user subroutine for residual function r.
C ADDA   = name of user subroutine for computing and adding A(t,y).
C JAC    = name of user subroutine for Jacobian matrix dr/dy
C          (MF = 21).  If not used, pass a dummy name.
C Note: the names for the RES and ADDA routines and (if used) the
C        JAC routine must be declared External in the calling program.
C NEQ    = number of scalar equations in the system.
C Y      = array of initial values, of length NEQ.
C YDOTI  = array of length NEQ (containing initial dy/dt if ISTATE = 1).
C T      = the initial value of the independent variable.
C TOUT   = first point where output is desired (.ne. T).
C ITOL   = 1 or 2 according as ATOL (below) is a scalar or array.
C RTOL   = relative tolerance parameter (scalar).
C ATOL   = absolute tolerance parameter (scalar or array).
C          the estimated local error in y(i) will be controlled so as
C          to be roughly less (in magnitude) than
C             EWT(i) = RTOL*ABS(Y(i)) + ATOL     if ITOL = 1, or
C             EWT(i) = RTOL*ABS(Y(i)) + ATOL(i)  if ITOL = 2.
C          Thus the local error test passes if, in each component,
C          either the absolute error is less than ATOL (or ATOL(i)),
C          or the relative error is less than RTOL.
C          Use RTOL = 0.0 for pure absolute error control, and
C          use ATOL = 0.0 (or ATOL(i) = 0.0) for pure relative error
C          control.  Caution: Actual (global) errors may exceed these
C          local tolerances, so choose them conservatively.
C ITASK  = 1 for normal computation of output values of y at t = TOUT.
C ISTATE = integer flag (input and output).  Set ISTATE = 1 if the
C          initial dy/dt is supplied, and 0 otherwise.
C IOPT   = 0 to indicate no optional inputs used.
C RWORK  = real work array of length at least:
C             22 + 9*NEQ + 3*MB*MB*NB        for MF = 21 or 22.
C LRW    = declared length of RWORK (in user's dimension).
C IWORK  = integer work array of length at least 20 + NEQ.
C          Input in IWORK(1) the block size MB and in IWORK(2) the
C          number NB of blocks in each direction along the matrix A.
C          These must satisfy  MB .ge. 1, NB .ge. 4, and MB*NB = NEQ.
C LIW    = declared length of IWORK (in user's dimension).
C MF     = method flag.  Standard values are:
C          21 for a user-supplied Jacobian.
C          22 for an internally generated Jacobian.
C          For other choices of MF, see the paragraph on MF in
C          the full description below.
C Note that the main program must declare arrays Y, YDOTI, RWORK, IWORK,
C and possibly ATOL.
C
C G. The output from the first call (or any call) is:
C      Y = array of computed values of y(t) vector.
C      T = corresponding value of independent variable (normally TOUT).
C ISTATE = 2  if DLSOIBT was successful, negative otherwise.
C          -1 means excess work done on this call (check all inputs).
C          -2 means excess accuracy requested (tolerances too small).
C          -3 means illegal input detected (see printed message).
C          -4 means repeated error test failures (check all inputs).
C          -5 means repeated convergence failures (perhaps bad Jacobian
C             supplied or wrong choice of tolerances).
C          -6 means error weight became zero during problem. (Solution
C             component i vanished, and ATOL or ATOL(i) = 0.)
C          -7 cannot occur in casual use.
C          -8 means DLSOIBT was unable to compute the initial dy/dt.
C             In casual use, this means A(t,y) is initially singular.
C             Supply YDOTI and use ISTATE = 1 on the first call.
C
C  If DLSOIBT returns ISTATE = -1, -4, or -5, then the output of
C  DLSOIBT also includes YDOTI = array containing residual vector
C  r = g - A * dy/dt  evaluated at the current t, y, and dy/dt.
C
C H. To continue the integration after a successful return, simply
C reset TOUT and call DLSOIBT again.  No other parameters need be reset.
C
C-----------------------------------------------------------------------
C Example Problem.
C
C The following is an example problem, with the coding needed
C for its solution by DLSOIBT.  The problem comes from the partial
C differential equation (the Burgers equation)
C   du/dt  =  - u * du/dx  +  eta * d**2 u/dx**2,   eta = .05,
C on -1 .le. x .le. 1.  The boundary conditions are
C   du/dx = 0  at x = -1 and at x = 1.
C The initial profile is a square wave,
C   u = 1 in ABS(x) .lt. .5,  u = .5 at ABS(x) = .5,  u = 0 elsewhere.
C The PDE is discretized in x by a simplified Galerkin method,
C using piecewise linear basis functions, on a grid of 40 intervals.
C The equations at x = -1 and 1 use a 3-point difference approximation
C for the right-hand side.  The result is a system A * dy/dt = g(y),
C of size NEQ = 41, where y(i) is the approximation to u at x = x(i),
C with x(i) = -1 + (i-1)*delx, delx = 2/(NEQ-1) = .05.  The individual
C equations in the system are
C   dy(1)/dt = ( y(3) - 2*y(2) + y(1) ) * eta / delx**2,
C   dy(NEQ)/dt = ( y(NEQ-2) - 2*y(NEQ-1) + y(NEQ) ) * eta / delx**2,
C and for i = 2, 3, ..., NEQ-1,
C   (1/6) dy(i-1)/dt + (4/6) dy(i)/dt + (1/6) dy(i+1)/dt
C       = ( y(i-1)**2 - y(i+1)**2 ) / (4*delx)
C         + ( y(i+1) - 2*y(i) + y(i-1) ) * eta / delx**2.
C The following coding solves the problem with MF = 21, with output
C of solution statistics at t = .1, .2, .3, and .4, and of the
C solution vector at t = .4.  Here the block size is just MB = 1.
C
C     EXTERNAL RESID, ADDABT, JACBT
C     DOUBLE PRECISION ATOL, RTOL, RWORK, T, TOUT, Y, YDOTI
C     DIMENSION Y(41), YDOTI(41), RWORK(514), IWORK(61)
C     NEQ = 41
C     DO 10 I = 1,NEQ
C  10   Y(I) = 0.0
C     Y(11) = 0.5
C     DO 20 I = 12,30
C  20   Y(I) = 1.0
C     Y(31) = 0.5
C     T = 0.0
C     TOUT = 0.1
C     ITOL = 1
C     RTOL = 1.0D-4
C     ATOL = 1.0D-5
C     ITASK = 1
C     ISTATE = 0
C     IOPT = 0
C     LRW = 514
C     LIW = 61
C     IWORK(1) = 1
C     IWORK(2) = NEQ
C     MF = 21
C     DO 40 IO = 1,4
C       CALL DLSOIBT (RESID, ADDABT, JACBT, NEQ, Y, YDOTI, T, TOUT,
C    1     ITOL,RTOL,ATOL, ITASK, ISTATE, IOPT, RWORK,LRW,IWORK,LIW, MF)
C       WRITE (6,30) T, IWORK(11), IWORK(12), IWORK(13)
C  30   FORMAT(' At t =',F5.2,'   No. steps =',I4,'  No. r-s =',I4,
C    1         '  No. J-s =',I3)
C       IF (ISTATE .NE. 2)  GO TO 90
C       TOUT = TOUT + 0.1
C  40   CONTINUE
C     WRITE(6,50) (Y(I),I=1,NEQ)
C  50 FORMAT(/' Final solution values..'/9(5D12.4/))
C     STOP
C  90 WRITE(6,95) ISTATE
C  95 FORMAT(///' Error halt.. ISTATE =',I3)
C     STOP
C     END
C
C     SUBROUTINE RESID (N, T, Y, S, R, IRES)
C     DOUBLE PRECISION T, Y, S, R, ETA, DELX, EODSQ
C     DIMENSION Y(N), S(N), R(N)
C     DATA ETA/0.05/, DELX/0.05/
C     EODSQ = ETA/DELX**2
C     R(1) = EODSQ*(Y(3) - 2.0*Y(2) + Y(1)) - S(1)
C     NM1 = N - 1
C     DO 10 I = 2,NM1
C       R(I) = (Y(I-1)**2 - Y(I+1)**2)/(4.0*DELX)
C    1        + EODSQ*(Y(I+1) - 2.0*Y(I) + Y(I-1))
C    2        - (S(I-1) + 4.0*S(I) + S(I+1))/6.0
C  10   CONTINUE
C     R(N) = EODSQ*(Y(N-2) - 2.0*Y(NM1) + Y(N)) - S(N)
C     RETURN
C     END
C
C     SUBROUTINE ADDABT (N, T, Y, MB, NB, PA, PB, PC)
C     DOUBLE PRECISION T, Y, PA, PB, PC
C     DIMENSION Y(N), PA(MB,MB,NB), PB(MB,MB,NB), PC(MB,MB,NB)
C     PA(1,1,1) = PA(1,1,1) + 1.0
C     NM1 = N - 1
C     DO 10 K = 2,NM1
C       PA(1,1,K) = PA(1,1,K) + (4.0/6.0)
C       PB(1,1,K) = PB(1,1,K) + (1.0/6.0)
C       PC(1,1,K) = PC(1,1,K) + (1.0/6.0)
C  10   CONTINUE
C     PA(1,1,N) = PA(1,1,N) + 1.0
C     RETURN
C     END
C
C     SUBROUTINE JACBT (N, T, Y, S, MB, NB, PA, PB, PC)
C     DOUBLE PRECISION T, Y, S, PA, PB, PC, ETA, DELX, EODSQ
C     DIMENSION Y(N), S(N), PA(MB,MB,NB),PB(MB,MB,NB),PC(MB,MB,NB)
C     DATA ETA/0.05/, DELX/0.05/
C     EODSQ = ETA/DELX**2
C     PA(1,1,1) = EODSQ
C     PB(1,1,1) = -2.0*EODSQ
C     PC(1,1,1) = EODSQ
C     DO 10 K = 2,N
C       PA(1,1,K) = -2.0*EODSQ
C       PB(1,1,K) = -Y(K+1)*(0.5/DELX) + EODSQ
C       PC(1,1,K) = Y(K-1)*(0.5/DELX) + EODSQ
C  10   CONTINUE
C     PB(1,1,N) = EODSQ
C     PC(1,1,N) = -2.0*EODSQ
C     PA(1,1,N) = EODSQ
C     RETURN
C     END
C
C The output of this program (on a CDC-7600 in single precision)
C is as follows:
C
C At t = 0.10   No. steps =  35  No. r-s =  45  No. J-s =  9
C At t = 0.20   No. steps =  43  No. r-s =  54  No. J-s = 10
C At t = 0.30   No. steps =  48  No. r-s =  60  No. J-s = 11
C At t = 0.40   No. steps =  51  No. r-s =  64  No. J-s = 12
C
C Final solution values..
C  1.2747e-02  1.1997e-02  1.5560e-02  2.3767e-02  3.7224e-02
C  5.6646e-02  8.2645e-02  1.1557e-01  1.5541e-01  2.0177e-01
C  2.5397e-01  3.1104e-01  3.7189e-01  4.3530e-01  5.0000e-01
C  5.6472e-01  6.2816e-01  6.8903e-01  7.4612e-01  7.9829e-01
C  8.4460e-01  8.8438e-01  9.1727e-01  9.4330e-01  9.6281e-01
C  9.7632e-01  9.8426e-01  9.8648e-01  9.8162e-01  9.6617e-01
C  9.3374e-01  8.7535e-01  7.8236e-01  6.5321e-01  5.0003e-01
C  3.4709e-01  2.1876e-01  1.2771e-01  7.3671e-02  5.0642e-02
C  5.4496e-02
C
C-----------------------------------------------------------------------
C Full Description of User Interface to DLSOIBT.
C
C The user interface to DLSOIBT consists of the following parts.
C
C 1.   The call sequence to Subroutine DLSOIBT, which is a driver
C      routine for the solver.  This includes descriptions of both
C      the call sequence arguments and of user-supplied routines.
C      Following these descriptions is a description of
C      optional inputs available through the call sequence, and then
C      a description of optional outputs (in the work arrays).
C
C 2.   Descriptions of other routines in the DLSOIBT package that may be
C      (optionally) called by the user.  These provide the ability to
C      alter error message handling, save and restore the internal
C      Common, and obtain specified derivatives of the solution y(t).
C
C 3.   Descriptions of Common blocks to be declared in overlay
C      or similar environments, or to be saved when doing an interrupt
C      of the problem and continued solution later.
C
C 4.   Description of two routines in the DLSOIBT package, either of
C      which the user may replace with his/her own version, if desired.
C      These relate to the measurement of errors.
C
C-----------------------------------------------------------------------
C Part 1.  Call Sequence.
C
C The call sequence parameters used for input only are
C     RES, ADDA, JAC, NEQ, TOUT, ITOL, RTOL, ATOL, ITASK,
C     IOPT, LRW, LIW, MF,
C and those used for both input and output are
C     Y, T, ISTATE, YDOTI.
C The work arrays RWORK and IWORK are also used for additional and
C optional inputs and optional outputs.  (The term output here refers
C to the return from Subroutine DLSOIBT to the user's calling program.)
C
C The legality of input parameters will be thoroughly checked on the
C initial call for the problem, but not checked thereafter unless a
C change in input parameters is flagged by ISTATE = 3 on input.
C
C The descriptions of the call arguments are as follows.
C
C RES    = the name of the user-supplied subroutine which supplies
C          the residual vector for the ODE system, defined by
C            r = g(t,y) - A(t,y) * s
C          as a function of the scalar t and the vectors
C          s and y (s approximates dy/dt).  This subroutine
C          is to have the form
C              SUBROUTINE RES (NEQ, T, Y, S, R, IRES)
C              DOUBLE PRECISION T, Y(*), S(*), R(*)
C          where NEQ, T, Y, S, and IRES are input, and R and
C          IRES are output. Y, S, and R are arrays of length NEQ.
C             On input, IRES indicates how DLSOIBT will use the
C          returned array R, as follows:
C             IRES = 1  means that DLSOIBT needs the full residual,
C                       r = g - A*s, exactly.
C             IRES = -1 means that DLSOIBT is using R only to compute
C                       the Jacobian dr/dy by difference quotients.
C          The RES routine can ignore IRES, or it can omit some terms
C          if IRES = -1.  If A does not depend on y, then RES can
C          just return R = g when IRES = -1.  If g - A*s contains other
C          additive terms that are independent of y, these can also be
C          dropped, if done consistently, when IRES = -1.
C             The subroutine should set the flag IRES if it
C          encounters a halt condition or illegal input.
C          Otherwise, it should not reset IRES.  On output,
C             IRES = 1 or -1 represents a normal return, and
C          DLSOIBT continues integrating the ODE.  Leave IRES
C          unchanged from its input value.
C             IRES = 2 tells DLSOIBT to immediately return control
C          to the calling program, with ISTATE = 3.  This lets
C          the calling program change parameters of the problem
C          if necessary.
C             IRES = 3 represents an error condition (for example, an
C          illegal value of y).  DLSOIBT tries to integrate the system
C          without getting IRES = 3 from RES.  If it cannot, DLSOIBT
C          returns with ISTATE = -7 or -1.
C             On an DLSOIBT return with ISTATE = 3, -1, or -7, the
C          values of T and Y returned correspond to the last point
C          reached successfully without getting the flag IRES = 2 or 3.
C             The flag values IRES = 2 and 3 should not be used to
C          handle switches or root-stop conditions.  This is better
C          done by calling DLSOIBT in a one-step mode and checking the
C          stopping function for a sign change at each step.
C             If quantities computed in the RES routine are needed
C          externally to DLSOIBT, an extra call to RES should be made
C          for this purpose, for consistent and accurate results.
C          To get the current dy/dt for the S argument, use DINTDY.
C             RES must be declared External in the calling
C          program. See note below for more about RES.
C
C ADDA   = the name of the user-supplied subroutine which adds the
C          matrix A = A(t,y) to another matrix, P, stored in
C          block-tridiagonal form.  This routine is to have the form
C               SUBROUTINE ADDA (NEQ, T, Y, MB, NB, PA, PB, PC)
C               DOUBLE PRECISION T, Y(*), PA(MB,MB,NB), PB(MB,MB,NB),
C              1                 PC(MB,MB,NB)
C          where NEQ, T, Y, MB, NB, and the arrays PA, PB, and PC
C          are input, and the arrays PA, PB, and PC are output.
C          Y is an array of length NEQ, and the arrays PA, PB, PC
C          are all MB by MB by NB.
C             Here a block-tridiagonal structure is assumed for A(t,y),
C          and also for the matrix P to which A is added here,
C          as described in Paragraph B of the Summary of Usage above.
C          Thus the affect of ADDA should be the following:
C               DO 30 K = 1,NB
C                 DO 20 J = 1,MB
C                   DO 10 I = 1,MB
C                     PA(I,J,K) = PA(I,J,K) +
C                       ( (I,J) element of K-th diagonal block of A)
C                     PB(I,J,K) = PB(I,J,K) +
C                       ( (I,J) element of block (K,K+1) of A,
C                       or block (NB,NB-2) if K = NB)
C                     PC(I,J,K) = PC(I,J,K) +
C                       ( (I,J) element of block (K,K-1) of A,
C                       or block (1,3) if K = 1)
C           10        CONTINUE
C           20      CONTINUE
C           30    CONTINUE
C             ADDA must be declared External in the calling program.
C          See note below for more information about ADDA.
C
C JAC    = the name of the user-supplied subroutine which supplies
C          the Jacobian matrix, dr/dy, where r = g - A*s.  JAC is
C          required if MITER = 1.  Otherwise a dummy name can be
C          passed.  This subroutine is to have the form
C               SUBROUTINE JAC (NEQ, T, Y, S, MB, NB, PA, PB, PC)
C               DOUBLE PRECISION T, Y(*), S(*), PA(MB,MB,NB),
C              1                 PB(MB,MB,NB), PC(MB,MB,NB)
C          where NEQ, T, Y, S, MB, NB, and the arrays PA, PB, and PC
C          are input, and the arrays PA, PB, and PC are output.
C          Y and S are arrays of length NEQ, and the arrays PA, PB, PC
C          are all MB by MB by NB.
C          PA, PB, and PC are to be loaded with partial derivatives
C          (elements of the Jacobian matrix) on output, in terms of the
C          block-tridiagonal structure assumed, as described
C          in Paragraph B of the Summary of Usage above.
C          That is, load the diagonal blocks into PA, the
C          superdiagonal blocks (and block (NB,NB-2) ) into PB, and
C          the subdiagonal blocks (and block (1,3) ) into PC.
C          The blocks in block-row k of dr/dy are to be loaded into
C          PA(*,*,k), PB(*,*,k), and PC(*,*,k).
C          Thus the affect of JAC should be the following:
C               DO 30 K = 1,NB
C                 DO 20 J = 1,MB
C                   DO 10 I = 1,MB
C                     PA(I,J,K) = ( (I,J) element of
C                       K-th diagonal block of dr/dy)
C                     PB(I,J,K) = ( (I,J) element of block (K,K+1)
C                       of dr/dy, or block (NB,NB-2) if K = NB)
C                     PC(I,J,K) = ( (I,J) element of block (K,K-1)
C                       of dr/dy, or block (1,3) if K = 1)
C           10        CONTINUE
C           20      CONTINUE
C           30    CONTINUE
C               PA, PB, and PC are preset to zero by the solver,
C          so that only the nonzero elements need be loaded by JAC.
C          Each call to JAC is preceded by a call to RES with the same
C          arguments NEQ, T, Y, and S.  Thus to gain some efficiency,
C          intermediate quantities shared by both calculations may be
C          saved in a user Common block by RES and not recomputed by JAC
C          if desired.  Also, JAC may alter the Y array, if desired.
C               JAC need not provide dr/dy exactly.  A crude
C          approximation will do, so that DLSOIBT may be used when
C          A and dr/dy are not really block-tridiagonal, but are close
C          to matrices that are.
C               JAC must be declared External in the calling program.
C               See note below for more about JAC.
C
C    Note on RES, ADDA, and JAC:
C          These subroutines may access user-defined quantities in
C          NEQ(2),... and/or in Y(NEQ(1)+1),... if NEQ is an array
C          (dimensioned in the subroutines) and/or Y has length
C          exceeding NEQ(1).  However, these routines should not alter
C          NEQ(1), Y(1),...,Y(NEQ) or any other input variables.
C          See the descriptions of NEQ and Y below.
C
C NEQ    = the size of the system (number of first order ordinary
C          differential equations or scalar algebraic equations).
C          Used only for input.
C          NEQ may be decreased, but not increased, during the problem.
C          If NEQ is decreased (with ISTATE = 3 on input), the
C          remaining components of Y should be left undisturbed, if
C          these are to be accessed in RES, ADDA, or JAC.
C
C          Normally, NEQ is a scalar, and it is generally referred to
C          as a scalar in this user interface description.  However,
C          NEQ may be an array, with NEQ(1) set to the system size.
C          (The DLSOIBT package accesses only NEQ(1).)  In either case,
C          this parameter is passed as the NEQ argument in all calls
C          to RES, ADDA, and JAC.  Hence, if it is an array,
C          locations NEQ(2),... may be used to store other integer data
C          and pass it to RES, ADDA, or JAC.  Each such subroutine
C          must include NEQ in a Dimension statement in that case.
C
C Y      = a real array for the vector of dependent variables, of
C          length NEQ or more.  Used for both input and output on the
C          first call (ISTATE = 0 or 1), and only for output on other
C          calls.  On the first call, Y must contain the vector of
C          initial values.  On output, Y contains the computed solution
C          vector, evaluated at t.  If desired, the Y array may be used
C          for other purposes between calls to the solver.
C
C          This array is passed as the Y argument in all calls to RES,
C          ADDA, and JAC.  Hence its length may exceed NEQ,
C          and locations Y(NEQ+1),... may be used to store other real
C          data and pass it to RES, ADDA, or JAC.  (The DLSOIBT
C          package accesses only Y(1),...,Y(NEQ). )
C
C YDOTI  = a real array for the initial value of the vector
C          dy/dt and for work space, of dimension at least NEQ.
C
C          On input:
C            If ISTATE = 0 then DLSOIBT will compute the initial value
C          of dy/dt, if A is nonsingular.  Thus YDOTI will
C          serve only as work space and may have any value.
C            If ISTATE = 1 then YDOTI must contain the initial value
C          of dy/dt.
C            If ISTATE = 2 or 3 (continuation calls) then YDOTI
C          may have any value.
C            Note: If the initial value of A is singular, then
C          DLSOIBT cannot compute the initial value of dy/dt, so
C          it must be provided in YDOTI, with ISTATE = 1.
C
C          On output, when DLSOIBT terminates abnormally with ISTATE =
C          -1, -4, or -5, YDOTI will contain the residual
C          r = g(t,y) - A(t,y)*(dy/dt).  If r is large, t is near
C          its initial value, and YDOTI is supplied with ISTATE = 1,
C          there may have been an incorrect input value of
C          YDOTI = dy/dt, or the problem (as given to DLSOIBT)
C          may not have a solution.
C
C          If desired, the YDOTI array may be used for other
C          purposes between calls to the solver.
C
C T      = the independent variable.  On input, T is used only on the
C          first call, as the initial point of the integration.
C          On output, after each call, T is the value at which a
C          computed solution y is evaluated (usually the same as TOUT).
C          On an error return, T is the farthest point reached.
C
C TOUT   = the next value of t at which a computed solution is desired.
C          Used only for input.
C
C          When starting the problem (ISTATE = 0 or 1), TOUT may be
C          equal to T for one call, then should .ne. T for the next
C          call.  For the initial T, an input value of TOUT .ne. T is
C          used in order to determine the direction of the integration
C          (i.e. the algebraic sign of the step sizes) and the rough
C          scale of the problem.  Integration in either direction
C          (forward or backward in t) is permitted.
C
C          If ITASK = 2 or 5 (one-step modes), TOUT is ignored after
C          the first call (i.e. the first call with TOUT .ne. T).
C          Otherwise, TOUT is required on every call.
C
C          If ITASK = 1, 3, or 4, the values of TOUT need not be
C          monotone, but a value of TOUT which backs up is limited
C          to the current internal T interval, whose endpoints are
C          TCUR - HU and TCUR (see optional outputs, below, for
C          TCUR and HU).
C
C ITOL   = an indicator for the type of error control.  See
C          description below under ATOL.  Used only for input.
C
C RTOL   = a relative error tolerance parameter, either a scalar or
C          an array of length NEQ.  See description below under ATOL.
C          Input only.
C
C ATOL   = an absolute error tolerance parameter, either a scalar or
C          an array of length NEQ.  Input only.
C
C             The input parameters ITOL, RTOL, and ATOL determine
C          the error control performed by the solver.  The solver will
C          control the vector E = (E(i)) of estimated local errors
C          in y, according to an inequality of the form
C                      RMS-norm of ( E(i)/EWT(i) )   .le.   1,
C          where       EWT(i) = RTOL(i)*ABS(Y(i)) + ATOL(i),
C          and the RMS-norm (root-mean-square norm) here is
C          RMS-norm(v) = SQRT(sum v(i)**2 / NEQ).  Here EWT = (EWT(i))
C          is a vector of weights which must always be positive, and
C          the values of RTOL and ATOL should all be non-negative.
C          The following table gives the types (scalar/array) of
C          RTOL and ATOL, and the corresponding form of EWT(i).
C
C             ITOL    RTOL       ATOL          EWT(i)
C              1     scalar     scalar     RTOL*ABS(Y(i)) + ATOL
C              2     scalar     array      RTOL*ABS(Y(i)) + ATOL(i)
C              3     array      scalar     RTOL(i)*ABS(Y(i)) + ATOL
C              4     array      scalar     RTOL(i)*ABS(Y(i)) + ATOL(i)
C
C          When either of these parameters is a scalar, it need not
C          be dimensioned in the user's calling program.
C
C          If none of the above choices (with ITOL, RTOL, and ATOL
C          fixed throughout the problem) is suitable, more general
C          error controls can be obtained by substituting
C          user-supplied routines for the setting of EWT and/or for
C          the norm calculation.  See Part 4 below.
C
C          If global errors are to be estimated by making a repeated
C          run on the same problem with smaller tolerances, then all
C          components of RTOL and ATOL (i.e. of EWT) should be scaled
C          down uniformly.
C
C ITASK  = an index specifying the task to be performed.
C          Input only.  ITASK has the following values and meanings.
C          1  means normal computation of output values of y(t) at
C             t = TOUT (by overshooting and interpolating).
C          2  means take one step only and return.
C          3  means stop at the first internal mesh point at or
C             beyond t = TOUT and return.
C          4  means normal computation of output values of y(t) at
C             t = TOUT but without overshooting t = TCRIT.
C             TCRIT must be input as RWORK(1).  TCRIT may be equal to
C             or beyond TOUT, but not behind it in the direction of
C             integration.  This option is useful if the problem
C             has a singularity at or beyond t = TCRIT.
C          5  means take one step, without passing TCRIT, and return.
C             TCRIT must be input as RWORK(1).
C
C          Note:  If ITASK = 4 or 5 and the solver reaches TCRIT
C          (within roundoff), it will return T = TCRIT (exactly) to
C          indicate this (unless ITASK = 4 and TOUT comes before TCRIT,
C          in which case answers at t = TOUT are returned first).
C
C ISTATE = an index used for input and output to specify the
C          state of the calculation.
C
C          On input, the values of ISTATE are as follows.
C          0  means this is the first call for the problem, and
C             DLSOIBT is to compute the initial value of dy/dt
C             (while doing other initializations).  See note below.
C          1  means this is the first call for the problem, and
C             the initial value of dy/dt has been supplied in
C             YDOTI (DLSOIBT will do other initializations).
C             See note below.
C          2  means this is not the first call, and the calculation
C             is to continue normally, with no change in any input
C             parameters except possibly TOUT and ITASK.
C             (If ITOL, RTOL, and/or ATOL are changed between calls
C             with ISTATE = 2, the new values will be used but not
C             tested for legality.)
C          3  means this is not the first call, and the
C             calculation is to continue normally, but with
C             a change in input parameters other than
C             TOUT and ITASK.  Changes are allowed in
C             NEQ, ITOL, RTOL, ATOL, IOPT, LRW, LIW, MF, MB, NB,
C             and any of the optional inputs except H0.
C             (See IWORK description for MB and NB.)
C          Note:  A preliminary call with TOUT = T is not counted
C          as a first call here, as no initialization or checking of
C          input is done.  (Such a call is sometimes useful for the
C          purpose of outputting the initial conditions.)
C          Thus the first call for which TOUT .ne. T requires
C          ISTATE = 0 or 1 on input.
C
C          On output, ISTATE has the following values and meanings.
C           0 or 1  means nothing was done; TOUT = t and
C              ISTATE = 0 or 1 on input.
C           2  means that the integration was performed successfully.
C           3  means that the user-supplied Subroutine RES signalled
C              DLSOIBT to halt the integration and return (IRES = 2).
C              Integration as far as T was achieved with no occurrence
C              of IRES = 2, but this flag was set on attempting the
C              next step.
C          -1  means an excessive amount of work (more than MXSTEP
C              steps) was done on this call, before completing the
C              requested task, but the integration was otherwise
C              successful as far as T.  (MXSTEP is an optional input
C              and is normally 500.)  To continue, the user may
C              simply reset ISTATE to a value .gt. 1 and call again
C              (the excess work step counter will be reset to 0).
C              In addition, the user may increase MXSTEP to avoid
C              this error return (see below on optional inputs).
C          -2  means too much accuracy was requested for the precision
C              of the machine being used.  This was detected before
C              completing the requested task, but the integration
C              was successful as far as T.  To continue, the tolerance
C              parameters must be reset, and ISTATE must be set
C              to 3.  The optional output TOLSF may be used for this
C              purpose.  (Note: If this condition is detected before
C              taking any steps, then an illegal input return
C              (ISTATE = -3) occurs instead.)
C          -3  means illegal input was detected, before taking any
C              integration steps.  See written message for details.
C              Note:  If the solver detects an infinite loop of calls
C              to the solver with illegal input, it will cause
C              the run to stop.
C          -4  means there were repeated error test failures on
C              one attempted step, before completing the requested
C              task, but the integration was successful as far as T.
C              The problem may have a singularity, or the input
C              may be inappropriate.
C          -5  means there were repeated convergence test failures on
C              one attempted step, before completing the requested
C              task, but the integration was successful as far as T.
C              This may be caused by an inaccurate Jacobian matrix.
C          -6  means EWT(i) became zero for some i during the
C              integration.  Pure relative error control (ATOL(i) = 0.0)
C              was requested on a variable which has now vanished.
C              The integration was successful as far as T.
C          -7  means that the user-supplied Subroutine RES set
C              its error flag (IRES = 3) despite repeated tries by
C              DLSOIBT to avoid that condition.
C          -8  means that ISTATE was 0 on input but DLSOIBT was unable
C              to compute the initial value of dy/dt.  See the
C              printed message for details.
C
C          Note:  Since the normal output value of ISTATE is 2,
C          it does not need to be reset for normal continuation.
C          Similarly, ISTATE (= 3) need not be reset if RES told
C          DLSOIBT to return because the calling program must change
C          the parameters of the problem.
C          Also, since a negative input value of ISTATE will be
C          regarded as illegal, a negative output value requires the
C          user to change it, and possibly other inputs, before
C          calling the solver again.
C
C IOPT   = an integer flag to specify whether or not any optional
C          inputs are being used on this call.  Input only.
C          The optional inputs are listed separately below.
C          IOPT = 0 means no optional inputs are being used.
C                   Default values will be used in all cases.
C          IOPT = 1 means one or more optional inputs are being used.
C
C RWORK  = a real working array (double precision).
C          The length of RWORK must be at least
C             20 + NYH*(MAXORD + 1) + 3*NEQ + LENWM    where
C          NYH    = the initial value of NEQ,
C          MAXORD = 12 (if METH = 1) or 5 (if METH = 2) (unless a
C                   smaller value is given as an optional input),
C          LENWM  = 3*MB*MB*NB + 2.
C          (See MF description for the definition of METH.)
C          Thus if MAXORD has its default value and NEQ is constant,
C          this length is
C             22 + 16*NEQ + 3*MB*MB*NB     for MF = 11 or 12,
C             22 + 9*NEQ + 3*MB*MB*NB      for MF = 21 or 22.
C          The first 20 words of RWORK are reserved for conditional
C          and optional inputs and optional outputs.
C
C          The following word in RWORK is a conditional input:
C            RWORK(1) = TCRIT = critical value of t which the solver
C                       is not to overshoot.  Required if ITASK is
C                       4 or 5, and ignored otherwise.  (See ITASK.)
C
C LRW    = the length of the array RWORK, as declared by the user.
C          (This will be checked by the solver.)
C
C IWORK  = an integer work array.  The length of IWORK must be at least
C          20 + NEQ .  The first few words of IWORK are used for
C          additional and optional inputs and optional outputs.
C
C          The following 2 words in IWORK are additional required
C          inputs to DLSOIBT:
C            IWORK(1) = MB = block size
C            IWORK(2) = NB = number of blocks in the main diagonal
C          These must satisfy  MB .ge. 1, NB .ge. 4, and MB*NB = NEQ.
C
C LIW    = the length of the array IWORK, as declared by the user.
C          (This will be checked by the solver.)
C
C Note:  The work arrays must not be altered between calls to DLSOIBT
C for the same problem, except possibly for the additional and
C optional inputs, and except for the last 3*NEQ words of RWORK.
C The latter space is used for internal scratch space, and so is
C available for use by the user outside DLSOIBT between calls, if
C desired (but not for use by RES, ADDA, or JAC).
C
C MF     = the method flag.  used only for input.  The legal values of
C          MF are 11, 12, 21, and 22.
C          MF has decimal digits METH and MITER: MF = 10*METH + MITER.
C            METH indicates the basic linear multistep method:
C              METH = 1 means the implicit Adams method.
C              METH = 2 means the method based on Backward
C                       Differentiation Formulas (BDFS).
C                The BDF method is strongly preferred for stiff
C              problems, while the Adams method is preferred when the
C              problem is not stiff.  If the matrix A(t,y) is
C              nonsingular, stiffness here can be taken to mean that of
C              the explicit ODE system dy/dt = A-inverse * g.  If A is
C              singular, the concept of stiffness is not well defined.
C                If you do not know whether the problem is stiff, we
C              recommend using METH = 2.  If it is stiff, the advantage
C              of METH = 2 over METH = 1 will be great, while if it is
C              not stiff, the advantage of METH = 1 will be slight.
C              If maximum efficiency is important, some experimentation
C              with METH may be necessary.
C            MITER indicates the corrector iteration method:
C              MITER = 1 means chord iteration with a user-supplied
C                        block-tridiagonal Jacobian.
C              MITER = 2 means chord iteration with an internally
C                        generated (difference quotient) block-
C                        tridiagonal Jacobian approximation, using
C                        3*MB+1 extra calls to RES per dr/dy evaluation.
C              If MITER = 1, the user must supply a Subroutine JAC
C              (the name is arbitrary) as described above under JAC.
C              For MITER = 2, a dummy argument can be used.
C-----------------------------------------------------------------------
C Optional Inputs.
C
C The following is a list of the optional inputs provided for in the
C call sequence.  (See also Part 2.)  For each such input variable,
C this table lists its name as used in this documentation, its
C location in the call sequence, its meaning, and the default value.
C The use of any of these inputs requires IOPT = 1, and in that
C case all of these inputs are examined.  A value of zero for any
C of these optional inputs will cause the default value to be used.
C Thus to use a subset of the optional inputs, simply preload
C locations 5 to 10 in RWORK and IWORK to 0.0 and 0 respectively, and
C then set those of interest to nonzero values.
C
C Name    Location      Meaning and Default Value
C
C H0      RWORK(5)  the step size to be attempted on the first step.
C                   The default value is determined by the solver.
C
C HMAX    RWORK(6)  the maximum absolute step size allowed.
C                   The default value is infinite.
C
C HMIN    RWORK(7)  the minimum absolute step size allowed.
C                   The default value is 0.  (This lower bound is not
C                   enforced on the final step before reaching TCRIT
C                   when ITASK = 4 or 5.)
C
C MAXORD  IWORK(5)  the maximum order to be allowed.  The default
C                   value is 12 if METH = 1, and 5 if METH = 2.
C                   If MAXORD exceeds the default value, it will
C                   be reduced to the default value.
C                   If MAXORD is changed during the problem, it may
C                   cause the current order to be reduced.
C
C MXSTEP  IWORK(6)  maximum number of (internally defined) steps
C                   allowed during one call to the solver.
C                   The default value is 500.
C
C MXHNIL  IWORK(7)  maximum number of messages printed (per problem)
C                   warning that T + H = T on a step (H = step size).
C                   This must be positive to result in a non-default
C                   value.  The default value is 10.
C-----------------------------------------------------------------------
C Optional Outputs.
C
C As optional additional output from DLSOIBT, the variables listed
C below are quantities related to the performance of DLSOIBT
C which are available to the user.  These are communicated by way of
C the work arrays, but also have internal mnemonic names as shown.
C Except where stated otherwise, all of these outputs are defined
C on any successful return from DLSOIBT, and on any return with
C ISTATE = -1, -2, -4, -5, -6, or -7.  On a return with -3 (illegal
C input) or -8, they will be unchanged from their existing values
C (if any), except possibly for TOLSF, LENRW, and LENIW.
C On any error return, outputs relevant to the error will be defined,
C as noted below.
C
C Name    Location      Meaning
C
C HU      RWORK(11) the step size in t last used (successfully).
C
C HCUR    RWORK(12) the step size to be attempted on the next step.
C
C TCUR    RWORK(13) the current value of the independent variable
C                   which the solver has actually reached, i.e. the
C                   current internal mesh point in t.  On output, TCUR
C                   will always be at least as far as the argument
C                   T, but may be farther (if interpolation was done).
C
C TOLSF   RWORK(14) a tolerance scale factor, greater than 1.0,
C                   computed when a request for too much accuracy was
C                   detected (ISTATE = -3 if detected at the start of
C                   the problem, ISTATE = -2 otherwise).  If ITOL is
C                   left unaltered but RTOL and ATOL are uniformly
C                   scaled up by a factor of TOLSF for the next call,
C                   then the solver is deemed likely to succeed.
C                   (The user may also ignore TOLSF and alter the
C                   tolerance parameters in any other way appropriate.)
C
C NST     IWORK(11) the number of steps taken for the problem so far.
C
C NRE     IWORK(12) the number of residual evaluations (res calls)
C                   for the problem so far.
C
C NJE     IWORK(13) the number of Jacobian evaluations (each involving
C                   an evaluation of a and dr/dy) for the problem so
C                   far.  This equals the number of calls to ADDA and
C                   (if MITER = 1) to JAC, and the number of matrix
C                   LU decompositions.
C
C NQU     IWORK(14) the method order last used (successfully).
C
C NQCUR   IWORK(15) the order to be attempted on the next step.
C
C IMXER   IWORK(16) the index of the component of largest magnitude in
C                   the weighted local error vector ( E(i)/EWT(i) ),
C                   on an error return with ISTATE = -4 or -5.
C
C LENRW   IWORK(17) the length of RWORK actually required.
C                   This is defined on normal returns and on an illegal
C                   input return for insufficient storage.
C
C LENIW   IWORK(18) the length of IWORK actually required.
C                   This is defined on normal returns and on an illegal
C                   input return for insufficient storage.
C
C
C The following two arrays are segments of the RWORK array which
C may also be of interest to the user as optional outputs.
C For each array, the table below gives its internal name,
C its base address in RWORK, and its description.
C
C Name    Base Address      Description
C
C YH      21             the Nordsieck history array, of size NYH by
C                        (NQCUR + 1), where NYH is the initial value
C                        of NEQ.  For j = 0,1,...,NQCUR, column j+1
C                        of YH contains HCUR**j/factorial(j) times
C                        the j-th derivative of the interpolating
C                        polynomial currently representing the solution,
C                        evaluated at t = TCUR.
C
C ACOR     LENRW-NEQ+1   array of size NEQ used for the accumulated
C                        corrections on each step, scaled on output to
C                        represent the estimated local error in y on
C                        the last step.  This is the vector E in the
C                        description of the error control.  It is
C                        defined only on a return from DLSOIBT with
C                        ISTATE = 2.
C
C-----------------------------------------------------------------------
C Part 2.  Other Routines Callable.
C
C The following are optional calls which the user may make to
C gain additional capabilities in conjunction with DLSOIBT.
C (The routines XSETUN and XSETF are designed to conform to the
C SLATEC error handling package.)
C
C     Form of Call                  Function
C   CALL XSETUN(LUN)          Set the logical unit number, LUN, for
C                             output of messages from DLSOIBT, if
C                             the default is not desired.
C                             The default value of LUN is 6.
C
C   CALL XSETF(MFLAG)         Set a flag to control the printing of
C                             messages by DLSOIBT.
C                             MFLAG = 0 means do not print. (Danger:
C                             This risks losing valuable information.)
C                             MFLAG = 1 means print (the default).
C
C                             Either of the above calls may be made at
C                             any time and will take effect immediately.
C
C   CALL DSRCOM(RSAV,ISAV,JOB) saves and restores the contents of
C                             the internal Common blocks used by
C                             DLSOIBT (see Part 3 below).
C                             RSAV must be a real array of length 9
C                             or more, and ISAV must be an integer
C                             array of length 25 or more.
C                             JOB=1 means save Common into RSAV/ISAV.
C                             JOB=2 means restore Common from RSAV/ISAV.
C                                DSRCOM is useful if one is
C                             interrupting a run and restarting
C                             later, or alternating between two or
C                             more problems solved with DLSOIBT.
C
C   CALL DINTDY(,,,,,)        Provide derivatives of y, of various
C        (see below)          orders, at a specified point t, if
C                             desired.  It may be called only after
C                             a successful return from DLSOIBT.
C
C The detailed instructions for using DINTDY are as follows.
C The form of the call is:
C
C   CALL DINTDY (T, K, RWORK(21), NYH, DKY, IFLAG)
C
C The input parameters are:
C
C T         = value of independent variable where answers are desired
C             (normally the same as the t last returned by DLSOIBT).
C             For valid results, T must lie between TCUR - HU and TCUR.
C             (See optional outputs for TCUR and HU.)
C K         = integer order of the derivative desired.  K must satisfy
C             0 .le. K .le. NQCUR, where NQCUR is the current order
C             (see optional outputs).  The capability corresponding
C             to K = 0, i.e. computing y(t), is already provided
C             by DLSOIBT directly.  Since NQCUR .ge. 1, the first
C             derivative dy/dt is always available with DINTDY.
C RWORK(21) = the base address of the history array YH.
C NYH       = column length of YH, equal to the initial value of NEQ.
C
C The output parameters are:
C
C DKY       = a real array of length NEQ containing the computed value
C             of the K-th derivative of y(t).
C IFLAG     = integer flag, returned as 0 if K and T were legal,
C             -1 if K was illegal, and -2 if T was illegal.
C             On an error return, a message is also written.
C-----------------------------------------------------------------------
C Part 3.  Common Blocks.
C
C If DLSOIBT is to be used in an overlay situation, the user
C must declare, in the primary overlay, the variables in:
C   (1) the call sequence to DLSOIBT, and
C   (2) the internal Common block
C         /DLS001/  of length  34  (9 double precision words
C                      followed by 25 integer words).
C
C If DLSOIBT is used on a system in which the contents of internal
C Common blocks are not preserved between calls, the user should
C declare the above Common block in the calling program to insure
C that their contents are preserved.
C
C If the solution of a given problem by DLSOIBT is to be interrupted
C and then later continued, such as when restarting an interrupted run
C or alternating between two or more problems, the user should save,
C following the return from the last DLSOIBT call prior to the
C interruption, the contents of the call sequence variables and the
C internal Common blocks, and later restore these values before the
C next DLSOIBT call for that problem.  To save and restore the Common
C blocks, use Subroutine DSRCOM (see Part 2 above).
C
C-----------------------------------------------------------------------
C Part 4.  Optionally Replaceable Solver Routines.
C
C Below are descriptions of two routines in the DLSOIBT package which
C relate to the measurement of errors.  Either routine can be
C replaced by a user-supplied version, if desired.  However, since such
C a replacement may have a major impact on performance, it should be
C done only when absolutely necessary, and only with great caution.
C (Note: The means by which the package version of a routine is
C superseded by the user's version may be system-dependent.)
C
C (a) DEWSET.
C The following subroutine is called just before each internal
C integration step, and sets the array of error weights, EWT, as
C described under ITOL/RTOL/ATOL above:
C     SUBROUTINE DEWSET (NEQ, ITOL, RTOL, ATOL, YCUR, EWT)
C where NEQ, ITOL, RTOL, and ATOL are as in the DLSOIBT call sequence,
C YCUR contains the current dependent variable vector, and
C EWT is the array of weights set by DEWSET.
C
C If the user supplies this subroutine, it must return in EWT(i)
C (i = 1,...,NEQ) a positive quantity suitable for comparing errors
C in y(i) to.  The EWT array returned by DEWSET is passed to the DVNORM
C routine (see below), and also used by DLSOIBT in the computation
C of the optional output IMXER, the diagonal Jacobian approximation,
C and the increments for difference quotient Jacobians.
C
C In the user-supplied version of DEWSET, it may be desirable to use
C the current values of derivatives of y.  Derivatives up to order NQ
C are available from the history array YH, described above under
C optional outputs.  In DEWSET, YH is identical to the YCUR array,
C extended to NQ + 1 columns with a column length of NYH and scale
C factors of H**j/factorial(j).  On the first call for the problem,
C given by NST = 0, NQ is 1 and H is temporarily set to 1.0.
C NYH is the initial value of NEQ.  The quantities NQ, H, and NST
C can be obtained by including in DEWSET the statements:
C     DOUBLE PRECISION RLS
C     COMMON /DLS001/ RLS(9),ILS(25)
C     NQ = ILS(21)
C     NST = ILS(22)
C     H = RLS(3)
C Thus, for example, the current value of dy/dt can be obtained as
C YCUR(NYH+i)/H  (i=1,...,NEQ)  (and the division by H is
C unnecessary when NST = 0).
C
C (b) DVNORM.
C The following is a real function routine which computes the weighted
C root-mean-square norm of a vector v:
C     D = DVNORM (N, V, W)
C where:
C   N = the length of the vector,
C   V = real array of length N containing the vector,
C   W = real array of length N containing weights,
C   D = SQRT( (1/N) * sum(V(i)*W(i))**2 ).
C DVNORM is called with N = NEQ and with W(i) = 1.0/EWT(i), where
C EWT is as set by Subroutine DEWSET.
C
C If the user supplies this function, it should return a non-negative
C value of DVNORM suitable for use in the error control in DLSOIBT.
C None of the arguments should be altered by DVNORM.
C For example, a user-supplied DVNORM routine might:
C   -substitute a max-norm of (V(i)*W(i)) for the RMS-norm, or
C   -ignore some components of V in the norm, with the effect of
C    suppressing the error control on those components of y.
C-----------------------------------------------------------------------
C
C***REVISION HISTORY  (YYYYMMDD)
C 19840625  DATE WRITTEN
C 19870330  Major update: corrected comments throughout;
C           removed TRET from Common; rewrote EWSET with 4 loops;
C           fixed t test in INTDY; added Cray directives in STODI;
C           in STODI, fixed DELP init. and logic around PJAC call;
C           combined routines to save/restore Common;
C           passed LEVEL = 0 in error message calls (except run abort).
C 20010425  Major update: convert source lines to upper case;
C           added *DECK lines; changed from 1 to * in dummy dimensions;
C           changed names R1MACH/D1MACH to RUMACH/DUMACH;
C           renamed routines for uniqueness across single/double prec.;
C           converted intrinsic names to generic form;
C           removed ILLIN and NTREP (data loaded) from Common;
C           removed all 'own' variables from Common;
C           changed error messages to quoted strings;
C           replaced XERRWV/XERRWD with 1993 revised version;
C           converted prologues, comments, error messages to mixed case;
C           converted arithmetic IF statements to logical IF statements;
C           numerous corrections to prologues and internal comments.
C 20010507  Converted single precision source to double precision.
C 20020502  Corrected declarations in descriptions of user routines.
C
C-----------------------------------------------------------------------
C Other routines in the DLSOIBT package.
C
C In addition to Subroutine DLSOIBT, the DLSOIBT package includes the
C following subroutines and function routines:
C  DAIGBT   computes the initial value of the vector
C             dy/dt = A-inverse * g
C  DINTDY   computes an interpolated value of the y vector at t = TOUT.
C  DSTODI   is the core integrator, which does one step of the
C           integration and the associated error control.
C  DCFODE   sets all method coefficients and test constants.
C  DEWSET   sets the error weight vector EWT before each step.
C  DVNORM   computes the weighted RMS-norm of a vector.
C  DSRCOM   is a user-callable routine to save and restore
C           the contents of the internal Common blocks.
C  DPJIBT   computes and preprocesses the Jacobian matrix
C           and the Newton iteration matrix P.
C  DSLSBT   manages solution of linear system in chord iteration.
C  DDECBT and DSOLBT   are routines for solving block-tridiagonal
C           systems of linear algebraic equations.
C  DGEFA and DGESL   are routines from LINPACK for solving full
C           systems of linear algebraic equations.
C  DDOT     is one of the basic linear algebra modules (BLAS).
C  DUMACH   computes the unit roundoff in a machine-independent manner.
C  XERRWD, XSETUN, XSETF, IXSAV, and IUMACH  handle the printing of all
C           error messages and warnings.  XERRWD is machine-dependent.
C Note:  DVNORM, DDOT, DUMACH, IXSAV, and IUMACH are function routines.
C All the others are subroutines.
C
C-----------------------------------------------------------------------
      EXTERNAL DPJIBT, DSLSBT
      DOUBLE PRECISION DUMACH, DVNORM
      INTEGER ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, 
     1   LYH, LEWT, LACOR, LSAVR, LWM, LIWM, METH, MITER,
     2   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NRE, NJE, NQU
      INTEGER I, I1, I2, IER, IFLAG, IMXER, IRES, KGO,
     1   LENIW, LENRW, LENWM, LP, LYD0, MB, MORD, MXHNL0, MXSTP0, NB
      INTEGER INIT, MXSTEP, MXHNIL, NHNIL, NSLAST, NYH
      DOUBLE PRECISION CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
      DOUBLE PRECISION ATOLI, AYI, BIG, EWTI, H0, HMAX, HMX, RH, RTOLI,
     1   TCRIT, TDIST, TNEXT, TOL, TOLSF, TP, SIZE, SUM, W0
      DIMENSION MORD(2)
      LOGICAL IHIT
      CHARACTER*60 MSG
      SAVE INIT, MXSTEP, MXHNIL, NHNIL, NSLAST, NYH
C-----------------------------------------------------------------------
C The following internal Common block contains variables which are
C communicated between subroutines.  All real variables are listed
C first, followed by all integers.  The block is declared in
C Subroutines DLSOIBT, DINTDY, DSTODI, DPJIBT, DSLSBT
C-----------------------------------------------------------------------
      COMMON /DLS001/ CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,
     1   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, 
     2   LYH, LEWT, LACOR, LSAVR, LWM, LIWM, METH, MITER,
     3   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NRE, NJE, NQU
C
      DATA MORD(1),MORD(2)/12,5/, MXSTP0/500/, MXHNL0/10/
C-----------------------------------------------------------------------
C Block A.
C This code block is executed on every call.
C It tests ISTATE and ITASK for legality and branches appropriately.
C If ISTATE .gt. 1 but the flag INIT shows that initialization has
C not yet been done, an error return occurs.
C If ISTATE = 0 or 1 and TOUT = T, return immediately.
C-----------------------------------------------------------------------
      IF (ISTATE .LT. 0 .OR. ISTATE .GT. 3) GO TO 601
      IF (ITASK .LT. 1 .OR. ITASK .GT. 5) GO TO 602
      IF (ISTATE .LE. 1) GO TO 10
      IF (INIT .EQ. 0) GO TO 603
      IF (ISTATE .EQ. 2) GO TO 200
      GO TO 20
 10   INIT = 0
      IF (TOUT .EQ. T) RETURN
C-----------------------------------------------------------------------
C Block B.
C The next code block is executed for the initial call (ISTATE = 0 or 1)
C or for a continuation call with parameter changes (ISTATE = 3).
C It contains checking of all inputs and various initializations.
C
C First check legality of the non-optional inputs NEQ, ITOL, IOPT,
C MF, MB, and NB.
C-----------------------------------------------------------------------
 20   IF (NEQ(1) .LE. 0) GO TO 604
      IF (ISTATE .LE. 1) GO TO 25
      IF (NEQ(1) .GT. N) GO TO 605
 25   N = NEQ(1)
      IF (ITOL .LT. 1 .OR. ITOL .GT. 4) GO TO 606
      IF (IOPT .LT. 0 .OR. IOPT .GT. 1) GO TO 607
      METH = MF/10
      MITER = MF - 10*METH
      IF (METH .LT. 1 .OR. METH .GT. 2) GO TO 608
      IF (MITER .LT. 1 .OR. MITER .GT. 2) GO TO 608
      MB = IWORK(1)
      NB = IWORK(2)
      IF (MB .LT. 1 .OR. MB .GT. N) GO TO 609
      IF (NB .LT. 4) GO TO 610
      IF (MB*NB .NE. N) GO TO 609
C Next process and check the optional inputs. --------------------------
      IF (IOPT .EQ. 1) GO TO 40
      MAXORD = MORD(METH)
      MXSTEP = MXSTP0
      MXHNIL = MXHNL0
      IF (ISTATE .LE. 1) H0 = 0.0D0
      HMXI = 0.0D0
      HMIN = 0.0D0
      GO TO 60
 40   MAXORD = IWORK(5)
      IF (MAXORD .LT. 0) GO TO 611
      IF (MAXORD .EQ. 0) MAXORD = 100
      MAXORD = MIN(MAXORD,MORD(METH))
      MXSTEP = IWORK(6)
      IF (MXSTEP .LT. 0) GO TO 612
      IF (MXSTEP .EQ. 0) MXSTEP = MXSTP0
      MXHNIL = IWORK(7)
      IF (MXHNIL .LT. 0) GO TO 613
      IF (MXHNIL .EQ. 0) MXHNIL = MXHNL0
      IF (ISTATE .GT. 1) GO TO 50
      H0 = RWORK(5)
      IF ((TOUT - T)*H0 .LT. 0.0D0) GO TO 614
 50   HMAX = RWORK(6)
      IF (HMAX .LT. 0.0D0) GO TO 615
      HMXI = 0.0D0
      IF (HMAX .GT. 0.0D0) HMXI = 1.0D0/HMAX
      HMIN = RWORK(7)
      IF (HMIN .LT. 0.0D0) GO TO 616
C-----------------------------------------------------------------------
C Set work array pointers and check lengths LRW and LIW.
C Pointers to segments of RWORK and IWORK are named by prefixing L to
C the name of the segment.  E.g., the segment YH starts at RWORK(LYH).
C Segments of RWORK (in order) are denoted YH, WM, EWT, SAVR, ACOR.
C-----------------------------------------------------------------------
 60   LYH = 21
      IF (ISTATE .LE. 1) NYH = N
      LWM = LYH + (MAXORD + 1)*NYH
      LENWM = 3*MB*MB*NB + 2
      LEWT = LWM + LENWM
      LSAVR = LEWT + N
      LACOR = LSAVR + N
      LENRW = LACOR + N - 1
      IWORK(17) = LENRW
      LIWM = 1
      LENIW = 20 + N
      IWORK(18) = LENIW
      IF (LENRW .GT. LRW) GO TO 617
      IF (LENIW .GT. LIW) GO TO 618
C Check RTOL and ATOL for legality. ------------------------------------
      RTOLI = RTOL(1)
      ATOLI = ATOL(1)
      DO 70 I = 1,N
        IF (ITOL .GE. 3) RTOLI = RTOL(I)
        IF (ITOL .EQ. 2 .OR. ITOL .EQ. 4) ATOLI = ATOL(I)
        IF (RTOLI .LT. 0.0D0) GO TO 619
        IF (ATOLI .LT. 0.0D0) GO TO 620
 70     CONTINUE
      IF (ISTATE .LE. 1) GO TO 100
C If ISTATE = 3, set flag to signal parameter changes to DSTODI. -------
      JSTART = -1
      IF (NQ .LE. MAXORD) GO TO 90
C MAXORD was reduced below NQ.  Copy YH(*,MAXORD+2) into YDOTI.---------
      DO 80 I = 1,N
 80     YDOTI(I) = RWORK(I+LWM-1)
C Reload WM(1) = RWORK(lWM), since lWM may have changed. ---------------
 90   RWORK(LWM) = SQRT(UROUND)
      IF (N .EQ. NYH) GO TO 200
C NEQ was reduced.  Zero part of YH to avoid undefined references. -----
      I1 = LYH + L*NYH
      I2 = LYH + (MAXORD + 1)*NYH - 1
      IF (I1 .GT. I2) GO TO 200
      DO 95 I = I1,I2
 95     RWORK(I) = 0.0D0
      GO TO 200
C-----------------------------------------------------------------------
C Block C.
C The next block is for the initial call only (ISTATE = 0 or 1).
C It contains all remaining initializations, the call to DAIGBT
C (if ISTATE = 1), and the calculation of the initial step size.
C The error weights in EWT are inverted after being loaded.
C-----------------------------------------------------------------------
 100  UROUND = DUMACH()
      TN = T
      IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GO TO 105
      TCRIT = RWORK(1)
      IF ((TCRIT - TOUT)*(TOUT - T) .LT. 0.0D0) GO TO 625
      IF (H0 .NE. 0.0D0 .AND. (T + H0 - TCRIT)*H0 .GT. 0.0D0)
     1   H0 = TCRIT - T
 105  JSTART = 0
      RWORK(LWM) = SQRT(UROUND)
      NHNIL = 0
      NST = 0
      NRE = 0
      NJE = 0
      NSLAST = 0
      HU = 0.0D0
      NQU = 0
      CCMAX = 0.3D0
      MAXCOR = 3
      MSBP = 20
      MXNCF = 10
C Compute initial dy/dt, if necessary, and load it and initial Y into YH
      LYD0 = LYH + NYH
      LP = LWM + 1
      IF ( ISTATE .EQ. 1 )  GO TO 120
C DLSOIBT must compute initial dy/dt (LYD0 points to YH(*,2)). ---------
         CALL DAIGBT( RES, ADDA, NEQ, T, Y, RWORK(LYD0),
     1               MB, NB, RWORK(LP), IWORK(21), IER )
         NRE = NRE + 1
         IF (IER .LT. 0) GO TO 560
         IF (IER .GT. 0) GO TO 565
         DO 115  I = 1,N
  115       RWORK(I+LYH-1) = Y(I)
         GO TO 130
C Initial dy/dt was supplied.  Load into YH (LYD0 points to YH(*,2).). -
  120    DO 125  I = 1,N
            RWORK(I+LYH-1) = Y(I)
  125       RWORK(I+LYD0-1) = YDOTI(I)
C Load and invert the EWT array.  (H is temporarily set to 1.0.) -------
  130 CONTINUE
      NQ = 1
      H = 1.0D0
      CALL DEWSET (N, ITOL, RTOL, ATOL, RWORK(LYH), RWORK(LEWT))
      DO 135 I = 1,N
        IF (RWORK(I+LEWT-1) .LE. 0.0D0) GO TO 621
 135    RWORK(I+LEWT-1) = 1.0D0/RWORK(I+LEWT-1)
C-----------------------------------------------------------------------
C The coding below computes the step size, H0, to be attempted on the
C first step, unless the user has supplied a value for this.
C First check that TOUT - T differs significantly from zero.
C A scalar tolerance quantity TOL is computed, as MAX(RTOL(i))
C if this is positive, or MAX(ATOL(i)/ABS(Y(i))) otherwise, adjusted
C so as to be between 100*UROUND and 1.0E-3.
C Then the computed value H0 is given by..
C                                      NEQ
C   H0**2 = TOL / ( w0**-2 + (1/NEQ) * Sum ( YDOT(i)/ywt(i) )**2  )
C                                       1
C where   w0      = MAX ( ABS(T), ABS(TOUT) ),
C         YDOT(i) = i-th component of initial value of dy/dt,
C         ywt(i)  = EWT(i)/TOL  (a weight for y(i)).
C The sign of H0 is inferred from the initial values of TOUT and T.
C-----------------------------------------------------------------------
      IF (H0 .NE. 0.0D0) GO TO 180
      TDIST = ABS(TOUT - T)
      W0 = MAX(ABS(T),ABS(TOUT))
      IF (TDIST .LT. 2.0D0*UROUND*W0) GO TO 622
      TOL = RTOL(1)
      IF (ITOL .LE. 2) GO TO 145
      DO 140 I = 1,N
 140    TOL = MAX(TOL,RTOL(I))
 145  IF (TOL .GT. 0.0D0) GO TO 160
      ATOLI = ATOL(1)
      DO 150 I = 1,N
        IF (ITOL .EQ. 2 .OR. ITOL .EQ. 4) ATOLI = ATOL(I)
        AYI = ABS(Y(I))
        IF (AYI .NE. 0.0D0) TOL = MAX(TOL,ATOLI/AYI)
 150    CONTINUE
 160  TOL = MAX(TOL,100.0D0*UROUND)
      TOL = MIN(TOL,0.001D0)
      SUM = DVNORM (N, RWORK(LYD0), RWORK(LEWT))
      SUM = 1.0D0/(TOL*W0*W0) + TOL*SUM**2
      H0 = 1.0D0/SQRT(SUM)
      H0 = MIN(H0,TDIST)
      H0 = SIGN(H0,TOUT-T)
C Adjust H0 if necessary to meet HMAX bound. ---------------------------
 180  RH = ABS(H0)*HMXI
      IF (RH .GT. 1.0D0) H0 = H0/RH
C Load H with H0 and scale YH(*,2) by H0. ------------------------------
      H = H0
      DO 190 I = 1,N
 190    RWORK(I+LYD0-1) = H0*RWORK(I+LYD0-1)
      GO TO 270
C-----------------------------------------------------------------------
C Block D.
C The next code block is for continuation calls only (ISTATE = 2 or 3)
C and is to check stop conditions before taking a step.
C-----------------------------------------------------------------------
 200  NSLAST = NST
      GO TO (210, 250, 220, 230, 240), ITASK
 210  IF ((TN - TOUT)*H .LT. 0.0D0) GO TO 250
      CALL DINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      IF (IFLAG .NE. 0) GO TO 627
      T = TOUT
      GO TO 420
 220  TP = TN - HU*(1.0D0 + 100.0D0*UROUND)
      IF ((TP - TOUT)*H .GT. 0.0D0) GO TO 623
      IF ((TN - TOUT)*H .LT. 0.0D0) GO TO 250
      GO TO 400
 230  TCRIT = RWORK(1)
      IF ((TN - TCRIT)*H .GT. 0.0D0) GO TO 624
      IF ((TCRIT - TOUT)*H .LT. 0.0D0) GO TO 625
      IF ((TN - TOUT)*H .LT. 0.0D0) GO TO 245
      CALL DINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      IF (IFLAG .NE. 0) GO TO 627
      T = TOUT
      GO TO 420
 240  TCRIT = RWORK(1)
      IF ((TN - TCRIT)*H .GT. 0.0D0) GO TO 624
 245  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. 100.0D0*UROUND*HMX
      IF (IHIT) GO TO 400
      TNEXT = TN + H*(1.0D0 + 4.0D0*UROUND)
      IF ((TNEXT - TCRIT)*H .LE. 0.0D0) GO TO 250
      H = (TCRIT - TN)*(1.0D0 - 4.0D0*UROUND)
      IF (ISTATE .EQ. 2) JSTART = -2
C-----------------------------------------------------------------------
C Block E.
C The next block is normally executed for all calls and contains
C the call to the one-step core integrator DSTODI.
C
C This is a looping point for the integration steps.
C
C First check for too many steps being taken, update EWT (if not at
C start of problem), check for too much accuracy being requested, and
C check for H below the roundoff level in T.
C-----------------------------------------------------------------------
 250  CONTINUE
      IF ((NST-NSLAST) .GE. MXSTEP) GO TO 500
      CALL DEWSET (N, ITOL, RTOL, ATOL, RWORK(LYH), RWORK(LEWT))
      DO 260 I = 1,N
        IF (RWORK(I+LEWT-1) .LE. 0.0D0) GO TO 510
 260    RWORK(I+LEWT-1) = 1.0D0/RWORK(I+LEWT-1)
 270  TOLSF = UROUND*DVNORM (N, RWORK(LYH), RWORK(LEWT))
      IF (TOLSF .LE. 1.0D0) GO TO 280
      TOLSF = TOLSF*2.0D0
      IF (NST .EQ. 0) GO TO 626
      GO TO 520
 280  IF ((TN + H) .NE. TN) GO TO 290
      NHNIL = NHNIL + 1
      IF (NHNIL .GT. MXHNIL) GO TO 290
      MSG = 'DLSOIBT- Warning..Internal T (=R1) and H (=R2) are'
      CALL XERRWD (MSG, 50, 101, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG='      such that in the machine, T + H = T on the next step  '
      CALL XERRWD (MSG, 60, 101, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '     (H = step size). Solver will continue anyway.'
      CALL XERRWD (MSG, 50, 101, 0, 0, 0, 0, 2, TN, H)
      IF (NHNIL .LT. MXHNIL) GO TO 290
      MSG = 'DLSOIBT- Above warning has been issued I1 times.  '
      CALL XERRWD (MSG, 50, 102, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '     It will not be issued again for this problem.'
      CALL XERRWD (MSG, 50, 102, 0, 1, MXHNIL, 0, 0, 0.0D0, 0.0D0)
 290  CONTINUE
C-----------------------------------------------------------------------
C     CALL DSTODI(NEQ,Y,YH,NYH,YH1,EWT,SAVF,SAVR,ACOR,WM,IWM,RES,
C                 ADDA,JAC,DPJIBT,DSLSBT)
C Note: SAVF in DSTODI occupies the same space as YDOTI in DLSOIBT.
C-----------------------------------------------------------------------
      CALL DSTODI (NEQ, Y, RWORK(LYH), NYH, RWORK(LYH), RWORK(LEWT),
     1   YDOTI, RWORK(LSAVR), RWORK(LACOR), RWORK(LWM),
     2   IWORK(LIWM), RES, ADDA, JAC, DPJIBT, DSLSBT )
      KGO = 1 - KFLAG
      GO TO (300, 530, 540, 400, 550), KGO
C
C KGO = 1:success; 2:error test failure; 3:convergence failure;
C       4:RES ordered return; 5:RES returned error.
C-----------------------------------------------------------------------
C Block F.
C The following block handles the case of a successful return from the
C core integrator (KFLAG = 0).  Test for stop conditions.
C-----------------------------------------------------------------------
 300  INIT = 1
      GO TO (310, 400, 330, 340, 350), ITASK
C ITASK = 1.  If TOUT has been reached, interpolate. -------------------
 310  IF ((TN - TOUT)*H .LT. 0.0D0) GO TO 250
      CALL DINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      T = TOUT
      GO TO 420
C ITASK = 3.  Jump to exit if TOUT was reached. ------------------------
 330  IF ((TN - TOUT)*H .GE. 0.0D0) GO TO 400
      GO TO 250
C ITASK = 4.  See if TOUT or TCRIT was reached.  Adjust H if necessary.
 340  IF ((TN - TOUT)*H .LT. 0.0D0) GO TO 345
      CALL DINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      T = TOUT
      GO TO 420
 345  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. 100.0D0*UROUND*HMX
      IF (IHIT) GO TO 400
      TNEXT = TN + H*(1.0D0 + 4.0D0*UROUND)
      IF ((TNEXT - TCRIT)*H .LE. 0.0D0) GO TO 250
      H = (TCRIT - TN)*(1.0D0 - 4.0D0*UROUND)
      JSTART = -2
      GO TO 250
C ITASK = 5.  see if TCRIT was reached and jump to exit. ---------------
 350  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. 100.0D0*UROUND*HMX
C-----------------------------------------------------------------------
C Block G.
C The following block handles all successful returns from DLSOIBT.
C If ITASK .ne. 1, Y is loaded from YH and T is set accordingly.
C ISTATE is set to 2, and the optional outputs are loaded into the
C work arrays before returning.
C-----------------------------------------------------------------------
 400  DO 410 I = 1,N
 410    Y(I) = RWORK(I+LYH-1)
      T = TN
      IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GO TO 420
      IF (IHIT) T = TCRIT
  420 ISTATE = 2
      IF ( KFLAG .EQ. -3 )  ISTATE = 3
      RWORK(11) = HU
      RWORK(12) = H
      RWORK(13) = TN
      IWORK(11) = NST
      IWORK(12) = NRE
      IWORK(13) = NJE
      IWORK(14) = NQU
      IWORK(15) = NQ
      RETURN
C-----------------------------------------------------------------------
C Block H.
C The following block handles all unsuccessful returns other than
C those for illegal input.  First the error message routine is called.
C If there was an error test or convergence test failure, IMXER is set.
C Then Y is loaded from YH and T is set to TN.
C The optional outputs are loaded into the work arrays before returning.
C-----------------------------------------------------------------------
C The maximum number of steps was taken before reaching TOUT. ----------
 500  MSG = 'DLSOIBT- At current T (=R1), MXSTEP (=I1) steps   '
      CALL XERRWD (MSG, 50, 201, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '      taken on this call before reaching TOUT     '
      CALL XERRWD (MSG, 50, 201, 0, 1, MXSTEP, 0, 1, TN, 0.0D0)
      ISTATE = -1
      GO TO 580
C EWT(i) .le. 0.0 for some i (not at start of problem). ----------------
 510  EWTI = RWORK(LEWT+I-1)
      MSG = 'DLSOIBT- At T (=R1), EWT(I1) has become R2 .le. 0.'
      CALL XERRWD (MSG, 50, 202, 0, 1, I, 0, 2, TN, EWTI)
      ISTATE = -6
      GO TO 590
C Too much accuracy requested for machine precision. -------------------
 520  MSG = 'DLSOIBT- At T (=R1), too much accuracy requested  '
      CALL XERRWD (MSG, 50, 203, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '      for precision of machine..  See TOLSF (=R2) '
      CALL XERRWD (MSG, 50, 203, 0, 0, 0, 0, 2, TN, TOLSF)
      RWORK(14) = TOLSF
      ISTATE = -2
      GO TO 590
C KFLAG = -1.  Error test failed repeatedly or with ABS(H) = HMIN. -----
 530  MSG = 'DLSOIBT- At T (=R1) and step size H (=R2), the    '
      CALL XERRWD (MSG, 50, 204, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = 'error test failed repeatedly or with ABS(H) = HMIN'
      CALL XERRWD (MSG, 50, 204, 0, 0, 0, 0, 2, TN, H)
      ISTATE = -4
      GO TO 570
C KFLAG = -2.  Convergence failed repeatedly or with ABS(H) = HMIN. ----
 540  MSG = 'DLSOIBT- At T (=R1) and step size H (=R2), the    '
      CALL XERRWD (MSG, 50, 205, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '      corrector convergence failed repeatedly     '
      CALL XERRWD (MSG, 50, 205, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '      or with ABS(H) = HMIN   '
      CALL XERRWD (MSG, 30, 205, 0, 0, 0, 0, 2, TN, H)
      ISTATE = -5
      GO TO 570
C IRES = 3 returned by RES, despite retries by DSTODI.------------------
 550  MSG = 'DLSOIBT- At T (=R1) residual routine returned     '
      CALL XERRWD (MSG, 50, 206, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '      error IRES = 3 repeatedly.        '
      CALL XERRWD (MSG, 40, 206, 0, 0, 0, 0, 1, TN, 0.0D0)
      ISTATE = -7
      GO TO 590
C DAIGBT failed because a diagonal block of A matrix was singular. -----
 560  IER = -IER
      MSG='DLSOIBT- Attempt to initialize dy/dt failed:  Matrix A has a'
      CALL XERRWD (MSG, 60, 207, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '      singular diagonal block, block no. = (I1)   '
      CALL XERRWD (MSG, 50, 207, 0, 1, IER, 0, 0, 0.0D0, 0.0D0)
      ISTATE = -8
      RETURN
C DAIGBT failed because RES set IRES to 2 or 3. ------------------------
 565  MSG = 'DLSOIBT- Attempt to initialize dy/dt failed       '
      CALL XERRWD (MSG, 50, 208, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '      because residual routine set its error flag '
      CALL XERRWD (MSG, 50, 208, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '      to IRES = (I1)'
      CALL XERRWD (MSG, 20, 208, 0, 1, IER, 0, 0, 0.0D0, 0.0D0)
      ISTATE = -8
      RETURN
C Compute IMXER if relevant. -------------------------------------------
 570  BIG = 0.0D0
      IMXER = 1
      DO 575 I = 1,N
        SIZE = ABS(RWORK(I+LACOR-1)*RWORK(I+LEWT-1))
        IF (BIG .GE. SIZE) GO TO 575
        BIG = SIZE
        IMXER = I
 575    CONTINUE
      IWORK(16) = IMXER
C Compute residual if relevant. ----------------------------------------
 580  LYD0 = LYH + NYH
      DO 585 I = 1,N
         RWORK(I+LSAVR-1) = RWORK(I+LYD0-1)/H
 585     Y(I) = RWORK(I+LYH-1)
      IRES = 1
      CALL RES (NEQ, TN, Y, RWORK(LSAVR), YDOTI, IRES)
      NRE = NRE + 1
      IF (IRES .LE. 1)  GO TO 595
      MSG = 'DLSOIBT- Residual routine set its flag IRES       '
      CALL XERRWD (MSG, 50, 210, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '      to (I1) when called for final output.       '
      CALL XERRWD (MSG, 50, 210, 0, 1, IRES, 0, 0, 0.0D0, 0.0D0)
      GO TO 595
C Set Y vector, T, and optional outputs. -------------------------------
 590  DO 592 I = 1,N
 592    Y(I) = RWORK(I+LYH-1)
 595  T = TN
      RWORK(11) = HU
      RWORK(12) = H
      RWORK(13) = TN
      IWORK(11) = NST
      IWORK(12) = NRE
      IWORK(13) = NJE
      IWORK(14) = NQU
      IWORK(15) = NQ
      RETURN
C-----------------------------------------------------------------------
C Block I.
C The following block handles all error returns due to illegal input
C (ISTATE = -3), as detected before calling the core integrator.
C First the error message routine is called.  If the illegal input
C is a negative ISTATE, the run is aborted (apparent infinite loop).
C-----------------------------------------------------------------------
 601  MSG = 'DLSOIBT- ISTATE (=I1) illegal.'
      CALL XERRWD (MSG, 30, 1, 0, 1, ISTATE, 0, 0, 0.0D0, 0.0D0)
      IF (ISTATE .LT. 0) GO TO 800
      GO TO 700
 602  MSG = 'DLSOIBT- ITASK (=I1) illegal. '
      CALL XERRWD (MSG, 30, 2, 0, 1, ITASK, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 603  MSG = 'DLSOIBT- ISTATE.gt.1 but DLSOIBT not initialized. '
      CALL XERRWD (MSG, 50, 3, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 604  MSG = 'DLSOIBT- NEQ (=I1) .lt. 1     '
      CALL XERRWD (MSG, 30, 4, 0, 1, NEQ(1), 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 605  MSG = 'DLSOIBT- ISTATE = 3 and NEQ increased (I1 to I2). '
      CALL XERRWD (MSG, 50, 5, 0, 2, N, NEQ(1), 0, 0.0D0, 0.0D0)
      GO TO 700
 606  MSG = 'DLSOIBT- ITOL (=I1) illegal.  '
      CALL XERRWD (MSG, 30, 6, 0, 1, ITOL, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 607  MSG = 'DLSOIBT- IOPT (=I1) illegal.  '
      CALL XERRWD (MSG, 30, 7, 0, 1, IOPT, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 608  MSG = 'DLSOIBT- MF (=I1) illegal.    '
      CALL XERRWD (MSG, 30, 8, 0, 1, MF, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 609  MSG = 'DLSOIBT- MB (=I1) or NB (=I2) illegal.  '
      CALL XERRWD (MSG, 40, 9, 0, 2, MB, NB, 0, 0.0D0, 0.0D0)
      GO TO 700
 610  MSG = 'DLSOIBT- NB (=I1) .lt. 4 illegal.       '
      CALL XERRWD (MSG, 40, 10, 0, 1, NB, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 611  MSG = 'DLSOIBT- MAXORD (=I1) .lt. 0  '
      CALL XERRWD (MSG, 30, 11, 0, 1, MAXORD, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 612  MSG = 'DLSOIBT- MXSTEP (=I1) .lt. 0  '
      CALL XERRWD (MSG, 30, 12, 0, 1, MXSTEP, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 613  MSG = 'DLSOIBT- MXHNIL (=I1) .lt. 0  '
      CALL XERRWD (MSG, 30, 13, 0, 1, MXHNIL, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 614  MSG = 'DLSOIBT- TOUT (=R1) behind T (=R2)      '
      CALL XERRWD (MSG, 40, 14, 0, 0, 0, 0, 2, TOUT, T)
      MSG = '      Integration direction is given by H0 (=R1)  '
      CALL XERRWD (MSG, 50, 14, 0, 0, 0, 0, 1, H0, 0.0D0)
      GO TO 700
 615  MSG = 'DLSOIBT- HMAX (=R1) .lt. 0.0  '
      CALL XERRWD (MSG, 30, 15, 0, 0, 0, 0, 1, HMAX, 0.0D0)
      GO TO 700
 616  MSG = 'DLSOIBT- HMIN (=R1) .lt. 0.0  '
      CALL XERRWD (MSG, 30, 16, 0, 0, 0, 0, 1, HMIN, 0.0D0)
      GO TO 700
 617  MSG='DLSOIBT- RWORK length needed, LENRW (=I1), exceeds LRW (=I2)'
      CALL XERRWD (MSG, 60, 17, 0, 2, LENRW, LRW, 0, 0.0D0, 0.0D0)
      GO TO 700
 618  MSG='DLSOIBT- IWORK length needed, LENIW (=I1), exceeds LIW (=I2)'
      CALL XERRWD (MSG, 60, 18, 0, 2, LENIW, LIW, 0, 0.0D0, 0.0D0)
      GO TO 700
 619  MSG = 'DLSOIBT- RTOL(=I1) is R1 .lt. 0.0       '
      CALL XERRWD (MSG, 40, 19, 0, 1, I, 0, 1, RTOLI, 0.0D0)
      GO TO 700
 620  MSG = 'DLSOIBT- ATOL(=I1) is R1 .lt. 0.0       '
      CALL XERRWD (MSG, 40, 20, 0, 1, I, 0, 1, ATOLI, 0.0D0)
      GO TO 700
 621  EWTI = RWORK(LEWT+I-1)
      MSG = 'DLSOIBT- EWT(I1) is R1 .le. 0.0         '
      CALL XERRWD (MSG, 40, 21, 0, 1, I, 0, 1, EWTI, 0.0D0)
      GO TO 700
 622  MSG='DLSOIBT- TOUT(=R1) too close to T(=R2) to start integration.'
      CALL XERRWD (MSG, 60, 22, 0, 0, 0, 0, 2, TOUT, T)
      GO TO 700
 623  MSG='DLSOIBT- ITASK = I1 and TOUT (=R1) behind TCUR - HU (= R2)  '
      CALL XERRWD (MSG, 60, 23, 0, 1, ITASK, 0, 2, TOUT, TP)
      GO TO 700
 624  MSG='DLSOIBT- ITASK = 4 or 5 and TCRIT (=R1) behind TCUR (=R2)   '
      CALL XERRWD (MSG, 60, 24, 0, 0, 0, 0, 2, TCRIT, TN)
      GO TO 700
 625  MSG='DLSOIBT- ITASK = 4 or 5 and TCRIT (=R1) behind TOUT (=R2)   '
      CALL XERRWD (MSG, 60, 25, 0, 0, 0, 0, 2, TCRIT, TOUT)
      GO TO 700
 626  MSG = 'DLSOIBT- At start of problem, too much accuracy   '
      CALL XERRWD (MSG, 50, 26, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG='      requested for precision of machine..  See TOLSF (=R1) '
      CALL XERRWD (MSG, 60, 26, 0, 0, 0, 0, 1, TOLSF, 0.0D0)
      RWORK(14) = TOLSF
      GO TO 700
 627  MSG = 'DLSOIBT- Trouble in DINTDY.  ITASK = I1, TOUT = R1'
      CALL XERRWD (MSG, 50, 27, 0, 1, ITASK, 0, 1, TOUT, 0.0D0)
C
 700  ISTATE = -3
      RETURN
C
 800  MSG = 'DLSOIBT- Run aborted.. apparent infinite loop.    '
      CALL XERRWD (MSG, 50, 303, 2, 0, 0, 0, 0, 0.0D0, 0.0D0)
      RETURN
C----------------------- End of Subroutine DLSOIBT ---------------------
      END



*DECK DUMACH
      DOUBLE PRECISION FUNCTION DUMACH ()
C***BEGIN PROLOGUE  DUMACH
C***PURPOSE  Compute the unit roundoff of the machine.
C***CATEGORY  R1
C***TYPE      DOUBLE PRECISION (RUMACH-S, DUMACH-D)
C***KEYWORDS  MACHINE CONSTANTS
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C***DESCRIPTION
C *Usage:
C        DOUBLE PRECISION  A, DUMACH
C        A = DUMACH()
C
C *Function Return Values:
C     A : the unit roundoff of the machine.
C
C *Description:
C     The unit roundoff is defined as the smallest positive machine
C     number u such that  1.0 + u .ne. 1.0.  This is computed by DUMACH
C     in a machine-independent manner.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   930216  DATE WRITTEN
C   930818  Added SLATEC-format prologue.  (FNF)
C***END PROLOGUE  DUMACH
C
C*Internal Notes:
C-----------------------------------------------------------------------
C Subroutines/functions called by DUMACH.. None
C-----------------------------------------------------------------------
C**End
C
      DOUBLE PRECISION U, COMP
C***FIRST EXECUTABLE STATEMENT  DUMACH
      U = 1.0D0
 10   U = U*0.5D0
      COMP = 1.0D0 + U
      IF (COMP .NE. 1.0D0) GO TO 10
      DUMACH = U*2.0D0
      RETURN
C----------------------- End of Function DUMACH ------------------------
      END
*DECK DCFODE
      SUBROUTINE DCFODE (METH, ELCO, TESCO)
C***BEGIN PROLOGUE  DCFODE
C***SUBSIDIARY
C***PURPOSE  Set ODE integrator coefficients.
C***TYPE      DOUBLE PRECISION (SCFODE-S, DCFODE-D)
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C***DESCRIPTION
C
C  DCFODE is called by the integrator routine to set coefficients
C  needed there.  The coefficients for the current method, as
C  given by the value of METH, are set for all orders and saved.
C  The maximum order assumed here is 12 if METH = 1 and 5 if METH = 2.
C  (A smaller value of the maximum order is also allowed.)
C  DCFODE is called once at the beginning of the problem,
C  and is not called again unless and until METH is changed.
C
C  The ELCO array contains the basic method coefficients.
C  The coefficients el(i), 1 .le. i .le. nq+1, for the method of
C  order nq are stored in ELCO(i,nq).  They are given by a genetrating
C  polynomial, i.e.,
C      l(x) = el(1) + el(2)*x + ... + el(nq+1)*x**nq.
C  For the implicit Adams methods, l(x) is given by
C      dl/dx = (x+1)*(x+2)*...*(x+nq-1)/factorial(nq-1),    l(-1) = 0.
C  For the BDF methods, l(x) is given by
C      l(x) = (x+1)*(x+2)* ... *(x+nq)/K,
C  where         K = factorial(nq)*(1 + 1/2 + ... + 1/nq).
C
C  The TESCO array contains test constants used for the
C  local error test and the selection of step size and/or order.
C  At order nq, TESCO(k,nq) is used for the selection of step
C  size at order nq - 1 if k = 1, at order nq if k = 2, and at order
C  nq + 1 if k = 3.
C
C***SEE ALSO  DLSODE
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791129  DATE WRITTEN
C   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
C   890503  Minor cosmetic changes.  (FNF)
C   930809  Renamed to allow single/double precision versions. (ACH)
C***END PROLOGUE  DCFODE
C**End
      INTEGER METH
      INTEGER I, IB, NQ, NQM1, NQP1
      DOUBLE PRECISION ELCO, TESCO
      DOUBLE PRECISION AGAMQ, FNQ, FNQM1, PC, PINT, RAGQ,
     1   RQFAC, RQ1FAC, TSIGN, XPIN
      DIMENSION ELCO(13,12), TESCO(3,12)
      DIMENSION PC(12)
C
C***FIRST EXECUTABLE STATEMENT  DCFODE
      GO TO (100, 200), METH
C
 100  ELCO(1,1) = 1.0D0
      ELCO(2,1) = 1.0D0
      TESCO(1,1) = 0.0D0
      TESCO(2,1) = 2.0D0
      TESCO(1,2) = 1.0D0
      TESCO(3,12) = 0.0D0
      PC(1) = 1.0D0
      RQFAC = 1.0D0
      DO 140 NQ = 2,12
C-----------------------------------------------------------------------
C The PC array will contain the coefficients of the polynomial
C     p(x) = (x+1)*(x+2)*...*(x+nq-1).
C Initially, p(x) = 1.
C-----------------------------------------------------------------------
        RQ1FAC = RQFAC
        RQFAC = RQFAC/NQ
        NQM1 = NQ - 1
        FNQM1 = NQM1
        NQP1 = NQ + 1
C Form coefficients of p(x)*(x+nq-1). ----------------------------------
        PC(NQ) = 0.0D0
        DO 110 IB = 1,NQM1
          I = NQP1 - IB
 110      PC(I) = PC(I-1) + FNQM1*PC(I)
        PC(1) = FNQM1*PC(1)
C Compute integral, -1 to 0, of p(x) and x*p(x). -----------------------
        PINT = PC(1)
        XPIN = PC(1)/2.0D0
        TSIGN = 1.0D0
        DO 120 I = 2,NQ
          TSIGN = -TSIGN
          PINT = PINT + TSIGN*PC(I)/I
 120      XPIN = XPIN + TSIGN*PC(I)/(I+1)
C Store coefficients in ELCO and TESCO. --------------------------------
        ELCO(1,NQ) = PINT*RQ1FAC
        ELCO(2,NQ) = 1.0D0
        DO 130 I = 2,NQ
 130      ELCO(I+1,NQ) = RQ1FAC*PC(I)/I
        AGAMQ = RQFAC*XPIN
        RAGQ = 1.0D0/AGAMQ
        TESCO(2,NQ) = RAGQ
        IF (NQ .LT. 12) TESCO(1,NQP1) = RAGQ*RQFAC/NQP1
        TESCO(3,NQM1) = RAGQ
 140    CONTINUE
      RETURN
C
 200  PC(1) = 1.0D0
      RQ1FAC = 1.0D0
      DO 230 NQ = 1,5
C-----------------------------------------------------------------------
C The PC array will contain the coefficients of the polynomial
C     p(x) = (x+1)*(x+2)*...*(x+nq).
C Initially, p(x) = 1.
C-----------------------------------------------------------------------
        FNQ = NQ
        NQP1 = NQ + 1
C Form coefficients of p(x)*(x+nq). ------------------------------------
        PC(NQP1) = 0.0D0
        DO 210 IB = 1,NQ
          I = NQ + 2 - IB
 210      PC(I) = PC(I-1) + FNQ*PC(I)
        PC(1) = FNQ*PC(1)
C Store coefficients in ELCO and TESCO. --------------------------------
        DO 220 I = 1,NQP1
 220      ELCO(I,NQ) = PC(I)/PC(2)
        ELCO(2,NQ) = 1.0D0
        TESCO(1,NQ) = RQ1FAC
        TESCO(2,NQ) = NQP1/ELCO(1,NQ)
        TESCO(3,NQ) = (NQ+2)/ELCO(1,NQ)
        RQ1FAC = RQ1FAC/FNQ
 230    CONTINUE
      RETURN
C----------------------- END OF SUBROUTINE DCFODE ----------------------
      END
*DECK DINTDY
      SUBROUTINE DINTDY (T, K, YH, NYH, DKY, IFLAG)
C***BEGIN PROLOGUE  DINTDY
C***SUBSIDIARY
C***PURPOSE  Interpolate solution derivatives.
C***TYPE      DOUBLE PRECISION (SINTDY-S, DINTDY-D)
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C***DESCRIPTION
C
C  DINTDY computes interpolated values of the K-th derivative of the
C  dependent variable vector y, and stores it in DKY.  This routine
C  is called within the package with K = 0 and T = TOUT, but may
C  also be called by the user for any K up to the current order.
C  (See detailed instructions in the usage documentation.)
C
C  The computed values in DKY are gotten by interpolation using the
C  Nordsieck history array YH.  This array corresponds uniquely to a
C  vector-valued polynomial of degree NQCUR or less, and DKY is set
C  to the K-th derivative of this polynomial at T.
C  The formula for DKY is:
C               q
C   DKY(i)  =  sum  c(j,K) * (T - tn)**(j-K) * h**(-j) * YH(i,j+1)
C              j=K
C  where  c(j,K) = j*(j-1)*...*(j-K+1), q = NQCUR, tn = TCUR, h = HCUR.
C  The quantities  nq = NQCUR, l = nq+1, N = NEQ, tn, and h are
C  communicated by COMMON.  The above sum is done in reverse order.
C  IFLAG is returned negative if either K or T is out of bounds.
C
C***SEE ALSO  DLSODE
C***ROUTINES CALLED  XERRWD
C***COMMON BLOCKS    DLS001
C***REVISION HISTORY  (YYMMDD)
C   791129  DATE WRITTEN
C   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
C   890503  Minor cosmetic changes.  (FNF)
C   930809  Renamed to allow single/double precision versions. (ACH)
C   010418  Reduced size of Common block /DLS001/. (ACH)
C***END PROLOGUE  DINTDY
C**End
      INTEGER K, NYH, IFLAG
      INTEGER ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, 
     1   LYH, LEWT, LACOR, LSAVR, LWM, LIWM, METH, MITER,
     2   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      INTEGER I, IC, J, JB, JB2, JJ, JJ1, JP1
      DOUBLE PRECISION T, YH, DKY
      DOUBLE PRECISION CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
      DOUBLE PRECISION C, R, S, TP
      CHARACTER*80 MSG
      DIMENSION YH(NYH,*), DKY(*)
      COMMON /DLS001/ CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,
     1   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, 
     2   LYH, LEWT, LACOR, LSAVR, LWM, LIWM, METH, MITER,
     3   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
C
C***FIRST EXECUTABLE STATEMENT  DINTDY
      IFLAG = 0
      IF (K .LT. 0 .OR. K .GT. NQ) GO TO 80
      TP = TN - HU -  100.0D0*UROUND*(TN + HU)
      IF ((T-TP)*(T-TN) .GT. 0.0D0) GO TO 90
C
      S = (T - TN)/H
      IC = 1
      IF (K .EQ. 0) GO TO 15
      JJ1 = L - K
      DO 10 JJ = JJ1,NQ
 10     IC = IC*JJ
 15   C = IC
      DO 20 I = 1,N
 20     DKY(I) = C*YH(I,L)
      IF (K .EQ. NQ) GO TO 55
      JB2 = NQ - K
      DO 50 JB = 1,JB2
        J = NQ - JB
        JP1 = J + 1
        IC = 1
        IF (K .EQ. 0) GO TO 35
        JJ1 = JP1 - K
        DO 30 JJ = JJ1,J
 30       IC = IC*JJ
 35     C = IC
        DO 40 I = 1,N
 40       DKY(I) = C*YH(I,JP1) + S*DKY(I)
 50     CONTINUE
      IF (K .EQ. 0) RETURN
 55   R = H**(-K)
      DO 60 I = 1,N
 60     DKY(I) = R*DKY(I)
      RETURN
C
 80   MSG = 'DINTDY-  K (=I1) illegal      '
      CALL XERRWD (MSG, 30, 51, 0, 1, K, 0, 0, 0.0D0, 0.0D0)
      IFLAG = -1
      RETURN
 90   MSG = 'DINTDY-  T (=R1) illegal      '
      CALL XERRWD (MSG, 30, 52, 0, 0, 0, 0, 1, T, 0.0D0)
      MSG='      T not in interval TCUR - HU (= R1) to TCUR (=R2)      '
      CALL XERRWD (MSG, 60, 52, 0, 0, 0, 0, 2, TP, TN)
      IFLAG = -2
      RETURN
C----------------------- END OF SUBROUTINE DINTDY ----------------------
      END
*DECK DPREPJ
      SUBROUTINE DPREPJ (NEQ, Y, YH, NYH, EWT, FTEM, SAVF, WM, IWM,
     1   F, JAC)
C***BEGIN PROLOGUE  DPREPJ
C***SUBSIDIARY
C***PURPOSE  Compute and process Newton iteration matrix.
C***TYPE      DOUBLE PRECISION (SPREPJ-S, DPREPJ-D)
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C***DESCRIPTION
C
C  DPREPJ is called by DSTODE to compute and process the matrix
C  P = I - h*el(1)*J , where J is an approximation to the Jacobian.
C  Here J is computed by the user-supplied routine JAC if
C  MITER = 1 or 4, or by finite differencing if MITER = 2, 3, or 5.
C  If MITER = 3, a diagonal approximation to J is used.
C  J is stored in WM and replaced by P.  If MITER .ne. 3, P is then
C  subjected to LU decomposition in preparation for later solution
C  of linear systems with P as coefficient matrix.  This is done
C  by DGEFA if MITER = 1 or 2, and by DGBFA if MITER = 4 or 5.
C
C  In addition to variables described in DSTODE and DLSODE prologues,
C  communication with DPREPJ uses the following:
C  Y     = array containing predicted values on entry.
C  FTEM  = work array of length N (ACOR in DSTODE).
C  SAVF  = array containing f evaluated at predicted y.
C  WM    = real work space for matrices.  On output it contains the
C          inverse diagonal matrix if MITER = 3 and the LU decomposition
C          of P if MITER is 1, 2 , 4, or 5.
C          Storage of matrix elements starts at WM(3).
C          WM also contains the following matrix-related data:
C          WM(1) = SQRT(UROUND), used in numerical Jacobian increments.
C          WM(2) = H*EL0, saved for later use if MITER = 3.
C  IWM   = integer work space containing pivot information, starting at
C          IWM(21), if MITER is 1, 2, 4, or 5.  IWM also contains band
C          parameters ML = IWM(1) and MU = IWM(2) if MITER is 4 or 5.
C  EL0   = EL(1) (input).
C  IERPJ = output error flag,  = 0 if no trouble, .gt. 0 if
C          P matrix found to be singular.
C  JCUR  = output flag = 1 to indicate that the Jacobian matrix
C          (or approximation) is now current.
C  This routine also uses the COMMON variables EL0, H, TN, UROUND,
C  MITER, N, NFE, and NJE.
C
C***SEE ALSO  DLSODE
C***ROUTINES CALLED  DGBFA, DGEFA, DVNORM
C***COMMON BLOCKS    DLS001
C***REVISION HISTORY  (YYMMDD)
C   791129  DATE WRITTEN
C   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
C   890504  Minor cosmetic changes.  (FNF)
C   930809  Renamed to allow single/double precision versions. (ACH)
C   010418  Reduced size of Common block /DLS001/. (ACH)
C***END PROLOGUE  DPREPJ
C**End
      EXTERNAL F, JAC
      INTEGER NEQ, NYH, IWM
      INTEGER ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, 
     1   LYH, LEWT, LACOR, LSAVR, LWM, LIWM, METH, MITER,
     2   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      INTEGER I, I1, I2, IER, II, J, J1, JJ, LENP,
     1   MBA, MBAND, MEB1, MEBAND, ML, ML3, MU, NP1
      DOUBLE PRECISION Y, YH, EWT, FTEM, SAVF, WM
      DOUBLE PRECISION CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
      DOUBLE PRECISION CON, DI, FAC, HL0, R, R0, SRUR, YI, YJ, YJJ,
     1   DVNORM
      DIMENSION NEQ(*), Y(*), YH(NYH,*), EWT(*), FTEM(*), SAVF(*),
     1   WM(*), IWM(*)
      COMMON /DLS001/ CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,
     1   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, 
     2   LYH, LEWT, LACOR, LSAVR, LWM, LIWM, METH, MITER,
     3   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
C
C***FIRST EXECUTABLE STATEMENT  DPREPJ
      NJE = NJE + 1
      IERPJ = 0
      JCUR = 1
      HL0 = H*EL0
      GO TO (100, 200, 300, 400, 500), MITER
C If MITER = 1, call JAC and multiply by scalar. -----------------------
 100  LENP = N*N
      DO 110 I = 1,LENP
 110    WM(I+2) = 0.0D0
      CALL JAC (NEQ, TN, Y, 0, 0, WM(3), N)
      CON = -HL0
      DO 120 I = 1,LENP
 120    WM(I+2) = WM(I+2)*CON
      GO TO 240
C If MITER = 2, make N calls to F to approximate J. --------------------
 200  FAC = DVNORM (N, SAVF, EWT)
      R0 = 1000.0D0*ABS(H)*UROUND*N*FAC
      IF (R0 .EQ. 0.0D0) R0 = 1.0D0
      SRUR = WM(1)
      J1 = 2
      DO 230 J = 1,N
        YJ = Y(J)
        R = MAX(SRUR*ABS(YJ),R0/EWT(J))
        Y(J) = Y(J) + R
        FAC = -HL0/R
        CALL F (NEQ, TN, Y, FTEM)
        DO 220 I = 1,N
 220      WM(I+J1) = (FTEM(I) - SAVF(I))*FAC
        Y(J) = YJ
        J1 = J1 + N
 230    CONTINUE
      NFE = NFE + N
C Add identity matrix. -------------------------------------------------
 240  J = 3
      NP1 = N + 1
      DO 250 I = 1,N
        WM(J) = WM(J) + 1.0D0
 250    J = J + NP1
C Do LU decomposition on P. --------------------------------------------
      CALL DGEFA (WM(3), N, N, IWM(21), IER)
      IF (IER .NE. 0) IERPJ = 1
      RETURN
C If MITER = 3, construct a diagonal approximation to J and P. ---------
 300  WM(2) = HL0
      R = EL0*0.1D0
      DO 310 I = 1,N
 310    Y(I) = Y(I) + R*(H*SAVF(I) - YH(I,2))
      CALL F (NEQ, TN, Y, WM(3))
      NFE = NFE + 1
      DO 320 I = 1,N
        R0 = H*SAVF(I) - YH(I,2)
        DI = 0.1D0*R0 - H*(WM(I+2) - SAVF(I))
        WM(I+2) = 1.0D0
        IF (ABS(R0) .LT. UROUND/EWT(I)) GO TO 320
        IF (ABS(DI) .EQ. 0.0D0) GO TO 330
        WM(I+2) = 0.1D0*R0/DI
 320    CONTINUE
      RETURN
 330  IERPJ = 1
      RETURN
C If MITER = 4, call JAC and multiply by scalar. -----------------------
 400  ML = IWM(1)
      MU = IWM(2)
      ML3 = ML + 3
      MBAND = ML + MU + 1
      MEBAND = MBAND + ML
      LENP = MEBAND*N
      DO 410 I = 1,LENP
 410    WM(I+2) = 0.0D0
      CALL JAC (NEQ, TN, Y, ML, MU, WM(ML3), MEBAND)
      CON = -HL0
      DO 420 I = 1,LENP
 420    WM(I+2) = WM(I+2)*CON
      GO TO 570
C If MITER = 5, make MBAND calls to F to approximate J. ----------------
 500  ML = IWM(1)
      MU = IWM(2)
      MBAND = ML + MU + 1
      MBA = MIN(MBAND,N)
      MEBAND = MBAND + ML
      MEB1 = MEBAND - 1
      SRUR = WM(1)
      FAC = DVNORM (N, SAVF, EWT)
      R0 = 1000.0D0*ABS(H)*UROUND*N*FAC
      IF (R0 .EQ. 0.0D0) R0 = 1.0D0
      DO 560 J = 1,MBA
        DO 530 I = J,N,MBAND
          YI = Y(I)
          R = MAX(SRUR*ABS(YI),R0/EWT(I))
 530      Y(I) = Y(I) + R
        CALL F (NEQ, TN, Y, FTEM)
        DO 550 JJ = J,N,MBAND
          Y(JJ) = YH(JJ,1)
          YJJ = Y(JJ)
          R = MAX(SRUR*ABS(YJJ),R0/EWT(JJ))
          FAC = -HL0/R
          I1 = MAX(JJ-MU,1)
          I2 = MIN(JJ+ML,N)
          II = JJ*MEB1 - ML + 2
          DO 540 I = I1,I2
 540        WM(II+I) = (FTEM(I) - SAVF(I))*FAC
 550      CONTINUE
 560    CONTINUE
      NFE = NFE + MBA
C Add identity matrix. -------------------------------------------------
 570  II = MBAND + 2
      DO 580 I = 1,N
        WM(II) = WM(II) + 1.0D0
 580    II = II + MEBAND
C Do LU decomposition of P. --------------------------------------------
      CALL DGBFA (WM(3), MEBAND, N, ML, MU, IWM(21), IER)
      IF (IER .NE. 0) IERPJ = 1
      RETURN
C----------------------- END OF SUBROUTINE DPREPJ ----------------------
      END
*DECK DSOLSY
      SUBROUTINE DSOLSY (WM, IWM, X, TEM)
C***BEGIN PROLOGUE  DSOLSY
C***SUBSIDIARY
C***PURPOSE  ODEPACK linear system solver.
C***TYPE      DOUBLE PRECISION (SSOLSY-S, DSOLSY-D)
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C***DESCRIPTION
C
C  This routine manages the solution of the linear system arising from
C  a chord iteration.  It is called if MITER .ne. 0.
C  If MITER is 1 or 2, it calls DGESL to accomplish this.
C  If MITER = 3 it updates the coefficient h*EL0 in the diagonal
C  matrix, and then computes the solution.
C  If MITER is 4 or 5, it calls DGBSL.
C  Communication with DSOLSY uses the following variables:
C  WM    = real work space containing the inverse diagonal matrix if
C          MITER = 3 and the LU decomposition of the matrix otherwise.
C          Storage of matrix elements starts at WM(3).
C          WM also contains the following matrix-related data:
C          WM(1) = SQRT(UROUND) (not used here),
C          WM(2) = HL0, the previous value of h*EL0, used if MITER = 3.
C  IWM   = integer work space containing pivot information, starting at
C          IWM(21), if MITER is 1, 2, 4, or 5.  IWM also contains band
C          parameters ML = IWM(1) and MU = IWM(2) if MITER is 4 or 5.
C  X     = the right-hand side vector on input, and the solution vector
C          on output, of length N.
C  TEM   = vector of work space of length N, not used in this version.
C  IERSL = output flag (in COMMON).  IERSL = 0 if no trouble occurred.
C          IERSL = 1 if a singular matrix arose with MITER = 3.
C  This routine also uses the COMMON variables EL0, H, MITER, and N.
C
C***SEE ALSO  DLSODE
C***ROUTINES CALLED  DGBSL, DGESL
C***COMMON BLOCKS    DLS001
C***REVISION HISTORY  (YYMMDD)
C   791129  DATE WRITTEN
C   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
C   890503  Minor cosmetic changes.  (FNF)
C   930809  Renamed to allow single/double precision versions. (ACH)
C   010418  Reduced size of Common block /DLS001/. (ACH)
C***END PROLOGUE  DSOLSY
C**End
      INTEGER IWM
      INTEGER ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, 
     1   LYH, LEWT, LACOR, LSAVR, LWM, LIWM, METH, MITER,
     2   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      INTEGER I, MEBAND, ML, MU
      DOUBLE PRECISION WM, X, TEM
      DOUBLE PRECISION CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
      DOUBLE PRECISION DI, HL0, PHL0, R
      DIMENSION WM(*), IWM(*), X(*), TEM(*)
      COMMON /DLS001/ CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,
     1   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, 
     2   LYH, LEWT, LACOR, LSAVR, LWM, LIWM, METH, MITER,
     3   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
C
C***FIRST EXECUTABLE STATEMENT  DSOLSY
      IERSL = 0
      GO TO (100, 100, 300, 400, 400), MITER
 100  CALL DGESL (WM(3), N, N, IWM(21), X, 0)
      RETURN
C
 300  PHL0 = WM(2)
      HL0 = H*EL0
      WM(2) = HL0
      IF (HL0 .EQ. PHL0) GO TO 330
      R = HL0/PHL0
      DO 320 I = 1,N
        DI = 1.0D0 - R*(1.0D0 - 1.0D0/WM(I+2))
        IF (ABS(DI) .EQ. 0.0D0) GO TO 390
 320    WM(I+2) = 1.0D0/DI
 330  DO 340 I = 1,N
 340    X(I) = WM(I+2)*X(I)
      RETURN
 390  IERSL = 1
      RETURN
C
 400  ML = IWM(1)
      MU = IWM(2)
      MEBAND = 2*ML + MU + 1
      CALL DGBSL (WM(3), MEBAND, N, ML, MU, IWM(21), X, 0)
      RETURN
C----------------------- END OF SUBROUTINE DSOLSY ----------------------
      END
*DECK DSRCOM
      SUBROUTINE DSRCOM (RSAV, ISAV, JOB)
C***BEGIN PROLOGUE  DSRCOM
C***SUBSIDIARY
C***PURPOSE  Save/restore ODEPACK COMMON blocks.
C***TYPE      DOUBLE PRECISION (SSRCOM-S, DSRCOM-D)
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C***DESCRIPTION
C
C  This routine saves or restores (depending on JOB) the contents of
C  the COMMON block DLS001, which is used internally
C  by one or more ODEPACK solvers.
C
C  RSAV = real array of length 9 or more.
C  ISAV = integer array of length 25 or more.
C  JOB  = flag indicating to save or restore the COMMON blocks:
C         JOB  = 1 if COMMON is to be saved (written to RSAV/ISAV)
C         JOB  = 2 if COMMON is to be restored (read from RSAV/ISAV)
C         A call with JOB = 2 presumes a prior call with JOB = 1.
C
C***SEE ALSO  DLSODE
C***ROUTINES CALLED  (NONE)
C***COMMON BLOCKS    DLS001
C***REVISION HISTORY  (YYMMDD)
C   791129  DATE WRITTEN
C   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
C   890503  Minor cosmetic changes.  (FNF)
C   921116  Deleted treatment of block /EH0001/.  (ACH)
C   930801  Reduced Common block length by 2.  (ACH)
C   930809  Renamed to allow single/double precision versions. (ACH)
C   010418  Reduced Common block length by 209+12. (ACH)
C***END PROLOGUE  DSRCOM
C**End
      INTEGER ISAV, JOB
      INTEGER ILS
      INTEGER I, LENILS, LENRLS
      DOUBLE PRECISION RSAV,   RLS
      DIMENSION RSAV(*), ISAV(*)
      COMMON /DLS001/ RLS(9), ILS(25)
      DATA LENRLS/9/, LENILS/25/
C
C***FIRST EXECUTABLE STATEMENT  DSRCOM
      IF (JOB .EQ. 2) GO TO 100
C
      DO 10 I = 1,LENRLS
 10     RSAV(I) = RLS(I)
      DO 20 I = 1,LENILS
 20     ISAV(I) = ILS(I)
      RETURN
C
 100  CONTINUE
      DO 110 I = 1,LENRLS
 110     RLS(I) = RSAV(I)
      DO 120 I = 1,LENILS
 120     ILS(I) = ISAV(I)
      RETURN
C----------------------- END OF SUBROUTINE DSRCOM ----------------------
      END
*DECK DSTODE
      SUBROUTINE DSTODE (NEQ, Y, YH, NYH, YH1, EWT, SAVF, ACOR,
     1   WM, IWM, F, JAC, PJAC, SLVS)
C***BEGIN PROLOGUE  DSTODE
C***SUBSIDIARY
C***PURPOSE  Performs one step of an ODEPACK integration.
C***TYPE      DOUBLE PRECISION (SSTODE-S, DSTODE-D)
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C***DESCRIPTION
C
C  DSTODE performs one step of the integration of an initial value
C  problem for a system of ordinary differential equations.
C  Note:  DSTODE is independent of the value of the iteration method
C  indicator MITER, when this is .ne. 0, and hence is independent
C  of the type of chord method used, or the Jacobian structure.
C  Communication with DSTODE is done with the following variables:
C
C  NEQ    = integer array containing problem size in NEQ(1), and
C           passed as the NEQ argument in all calls to F and JAC.
C  Y      = an array of length .ge. N used as the Y argument in
C           all calls to F and JAC.
C  YH     = an NYH by LMAX array containing the dependent variables
C           and their approximate scaled derivatives, where
C           LMAX = MAXORD + 1.  YH(i,j+1) contains the approximate
C           j-th derivative of y(i), scaled by h**j/factorial(j)
C           (j = 0,1,...,NQ).  on entry for the first step, the first
C           two columns of YH must be set from the initial values.
C  NYH    = a constant integer .ge. N, the first dimension of YH.
C  YH1    = a one-dimensional array occupying the same space as YH.
C  EWT    = an array of length N containing multiplicative weights
C           for local error measurements.  Local errors in Y(i) are
C           compared to 1.0/EWT(i) in various error tests.
C  SAVF   = an array of working storage, of length N.
C           Also used for input of YH(*,MAXORD+2) when JSTART = -1
C           and MAXORD .lt. the current order NQ.
C  ACOR   = a work array of length N, used for the accumulated
C           corrections.  On a successful return, ACOR(i) contains
C           the estimated one-step local error in Y(i).
C  WM,IWM = real and integer work arrays associated with matrix
C           operations in chord iteration (MITER .ne. 0).
C  PJAC   = name of routine to evaluate and preprocess Jacobian matrix
C           and P = I - h*el0*JAC, if a chord method is being used.
C  SLVS   = name of routine to solve linear system in chord iteration.
C  CCMAX  = maximum relative change in h*el0 before PJAC is called.
C  H      = the step size to be attempted on the next step.
C           H is altered by the error control algorithm during the
C           problem.  H can be either positive or negative, but its
C           sign must remain constant throughout the problem.
C  HMIN   = the minimum absolute value of the step size h to be used.
C  HMXI   = inverse of the maximum absolute value of h to be used.
C           HMXI = 0.0 is allowed and corresponds to an infinite hmax.
C           HMIN and HMXI may be changed at any time, but will not
C           take effect until the next change of h is considered.
C  TN     = the independent variable. TN is updated on each step taken.
C  JSTART = an integer used for input only, with the following
C           values and meanings:
C                0  perform the first step.
C            .gt.0  take a new step continuing from the last.
C               -1  take the next step with a new value of H, MAXORD,
C                     N, METH, MITER, and/or matrix parameters.
C               -2  take the next step with a new value of H,
C                     but with other inputs unchanged.
C           On return, JSTART is set to 1 to facilitate continuation.
C  KFLAG  = a completion code with the following meanings:
C                0  the step was succesful.
C               -1  the requested error could not be achieved.
C               -2  corrector convergence could not be achieved.
C               -3  fatal error in PJAC or SLVS.
C           A return with KFLAG = -1 or -2 means either
C           abs(H) = HMIN or 10 consecutive failures occurred.
C           On a return with KFLAG negative, the values of TN and
C           the YH array are as of the beginning of the last
C           step, and H is the last step size attempted.
C  MAXORD = the maximum order of integration method to be allowed.
C  MAXCOR = the maximum number of corrector iterations allowed.
C  MSBP   = maximum number of steps between PJAC calls (MITER .gt. 0).
C  MXNCF  = maximum number of convergence failures allowed.
C  METH/MITER = the method flags.  See description in driver.
C  N      = the number of first-order differential equations.
C  The values of CCMAX, H, HMIN, HMXI, TN, JSTART, KFLAG, MAXORD,
C  MAXCOR, MSBP, MXNCF, METH, MITER, and N are communicated via COMMON.
C
C***SEE ALSO  DLSODE
C***ROUTINES CALLED  DCFODE, DVNORM
C***COMMON BLOCKS    DLS001
C***REVISION HISTORY  (YYMMDD)
C   791129  DATE WRITTEN
C   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
C   890503  Minor cosmetic changes.  (FNF)
C   930809  Renamed to allow single/double precision versions. (ACH)
C   010418  Reduced size of Common block /DLS001/. (ACH)
C***END PROLOGUE  DSTODE
C**End
      EXTERNAL F, JAC, PJAC, SLVS
      INTEGER NEQ, NYH, IWM
      INTEGER ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, 
     1   LYH, LEWT, LACOR, LSAVR, LWM, LIWM, METH, MITER,
     2   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      INTEGER I, I1, IREDO, IRET, J, JB, M, NCF, NEWQ
      INTEGER IALTH, IPUP, LMAX, MEO, NQNYH, NSLP
      DOUBLE PRECISION Y, YH, YH1, EWT, SAVF, ACOR, WM
      DOUBLE PRECISION CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
      DOUBLE PRECISION DCON, DDN, DEL, DELP, DSM, DUP, EXDN, EXSM, EXUP,
     1   R, RH, RHDN, RHSM, RHUP, TOLD, DVNORM,
     2   CONIT, CRATE, EL(13), ELCO(13,12), HOLD, RMAX, TESCO(3,12)
      DIMENSION NEQ(*), Y(*), YH(NYH,*), YH1(*), EWT(*), SAVF(*),
     1   ACOR(*), WM(*), IWM(*)
      SAVE CONIT, CRATE, EL, ELCO, HOLD, RMAX, TESCO,
     1   IALTH, IPUP, LMAX, MEO, NQNYH, NSLP
      COMMON /DLS001/ CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,
     1   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, 
     2   LYH, LEWT, LACOR, LSAVR, LWM, LIWM, METH, MITER,
     3   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
C
C***FIRST EXECUTABLE STATEMENT  DSTODE
      KFLAG = 0
      TOLD = TN
      NCF = 0
      IERPJ = 0
      IERSL = 0
      JCUR = 0
      ICF = 0
      DELP = 0.0D0
      IF (JSTART .GT. 0) GO TO 200
      IF (JSTART .EQ. -1) GO TO 100
      IF (JSTART .EQ. -2) GO TO 160
C-----------------------------------------------------------------------
C On the first call, the order is set to 1, and other variables are
C initialized.  RMAX is the maximum ratio by which H can be increased
C in a single step.  It is initially 1.E4 to compensate for the small
C initial H, but then is normally equal to 10.  If a failure
C occurs (in corrector convergence or error test), RMAX is set to 2
C for the next increase.
C-----------------------------------------------------------------------
      LMAX = MAXORD + 1
      NQ = 1
      L = 2
      IALTH = 2
      RMAX = 10000.0D0
      RC = 0.0D0
      EL0 = 1.0D0
      CRATE = 0.7D0
      HOLD = H
      MEO = METH
      NSLP = 0
      IPUP = MITER
      IRET = 3
      GO TO 140
C-----------------------------------------------------------------------
C The following block handles preliminaries needed when JSTART = -1.
C IPUP is set to MITER to force a matrix update.
C If an order increase is about to be considered (IALTH = 1),
C IALTH is reset to 2 to postpone consideration one more step.
C If the caller has changed METH, DCFODE is called to reset
C the coefficients of the method.
C If the caller has changed MAXORD to a value less than the current
C order NQ, NQ is reduced to MAXORD, and a new H chosen accordingly.
C If H is to be changed, YH must be rescaled.
C If H or METH is being changed, IALTH is reset to L = NQ + 1
C to prevent further changes in H for that many steps.
C-----------------------------------------------------------------------
 100  IPUP = MITER
      LMAX = MAXORD + 1
      IF (IALTH .EQ. 1) IALTH = 2
      IF (METH .EQ. MEO) GO TO 110
      CALL DCFODE (METH, ELCO, TESCO)
      MEO = METH
      IF (NQ .GT. MAXORD) GO TO 120
      IALTH = L
      IRET = 1
      GO TO 150
 110  IF (NQ .LE. MAXORD) GO TO 160
 120  NQ = MAXORD
      L = LMAX
      DO 125 I = 1,L
 125    EL(I) = ELCO(I,NQ)
      NQNYH = NQ*NYH
      RC = RC*EL(1)/EL0
      EL0 = EL(1)
      CONIT = 0.5D0/(NQ+2)
      DDN = DVNORM (N, SAVF, EWT)/TESCO(1,L)
      EXDN = 1.0D0/L
      RHDN = 1.0D0/(1.3D0*DDN**EXDN + 0.0000013D0)
      RH = MIN(RHDN,1.0D0)
      IREDO = 3
      IF (H .EQ. HOLD) GO TO 170
      RH = MIN(RH,ABS(H/HOLD))
      H = HOLD
      GO TO 175
C-----------------------------------------------------------------------
C DCFODE is called to get all the integration coefficients for the
C current METH.  Then the EL vector and related constants are reset
C whenever the order NQ is changed, or at the start of the problem.
C-----------------------------------------------------------------------
 140  CALL DCFODE (METH, ELCO, TESCO)
 150  DO 155 I = 1,L
 155    EL(I) = ELCO(I,NQ)
      NQNYH = NQ*NYH
      RC = RC*EL(1)/EL0
      EL0 = EL(1)
      CONIT = 0.5D0/(NQ+2)
      GO TO (160, 170, 200), IRET
C-----------------------------------------------------------------------
C If H is being changed, the H ratio RH is checked against
C RMAX, HMIN, and HMXI, and the YH array rescaled.  IALTH is set to
C L = NQ + 1 to prevent a change of H for that many steps, unless
C forced by a convergence or error test failure.
C-----------------------------------------------------------------------
 160  IF (H .EQ. HOLD) GO TO 200
      RH = H/HOLD
      H = HOLD
      IREDO = 3
      GO TO 175
 170  RH = MAX(RH,HMIN/ABS(H))
 175  RH = MIN(RH,RMAX)
      RH = RH/MAX(1.0D0,ABS(H)*HMXI*RH)
      R = 1.0D0
      DO 180 J = 2,L
        R = R*RH
        DO 180 I = 1,N
 180      YH(I,J) = YH(I,J)*R
      H = H*RH
      RC = RC*RH
      IALTH = L
      IF (IREDO .EQ. 0) GO TO 690
C-----------------------------------------------------------------------
C This section computes the predicted values by effectively
C multiplying the YH array by the Pascal Triangle matrix.
C RC is the ratio of new to old values of the coefficient  H*EL(1).
C When RC differs from 1 by more than CCMAX, IPUP is set to MITER
C to force PJAC to be called, if a Jacobian is involved.
C In any case, PJAC is called at least every MSBP steps.
C-----------------------------------------------------------------------
 200  IF (ABS(RC-1.0D0) .GT. CCMAX) IPUP = MITER
      IF (NST .GE. NSLP+MSBP) IPUP = MITER
      TN = TN + H
      I1 = NQNYH + 1
      DO 215 JB = 1,NQ
        I1 = I1 - NYH
Cdir$ ivdep
        DO 210 I = I1,NQNYH
 210      YH1(I) = YH1(I) + YH1(I+NYH)
 215    CONTINUE
C-----------------------------------------------------------------------
C Up to MAXCOR corrector iterations are taken.  A convergence test is
C made on the R.M.S. norm of each correction, weighted by the error
C weight vector EWT.  The sum of the corrections is accumulated in the
C vector ACOR(i).  The YH array is not altered in the corrector loop.
C-----------------------------------------------------------------------
 220  M = 0
      DO 230 I = 1,N
 230    Y(I) = YH(I,1)
      CALL F (NEQ, TN, Y, SAVF)
      NFE = NFE + 1
      IF (IPUP .LE. 0) GO TO 250
C-----------------------------------------------------------------------
C If indicated, the matrix P = I - h*el(1)*J is reevaluated and
C preprocessed before starting the corrector iteration.  IPUP is set
C to 0 as an indicator that this has been done.
C-----------------------------------------------------------------------
      CALL PJAC (NEQ, Y, YH, NYH, EWT, ACOR, SAVF, WM, IWM, F, JAC)
      IPUP = 0
      RC = 1.0D0
      NSLP = NST
      CRATE = 0.7D0
      IF (IERPJ .NE. 0) GO TO 430
 250  DO 260 I = 1,N
 260    ACOR(I) = 0.0D0
 270  IF (MITER .NE. 0) GO TO 350
C-----------------------------------------------------------------------
C In the case of functional iteration, update Y directly from
C the result of the last function evaluation.
C-----------------------------------------------------------------------
      DO 290 I = 1,N
        SAVF(I) = H*SAVF(I) - YH(I,2)
 290    Y(I) = SAVF(I) - ACOR(I)
      DEL = DVNORM (N, Y, EWT)
      DO 300 I = 1,N
        Y(I) = YH(I,1) + EL(1)*SAVF(I)
 300    ACOR(I) = SAVF(I)
      GO TO 400
C-----------------------------------------------------------------------
C In the case of the chord method, compute the corrector error,
C and solve the linear system with that as right-hand side and
C P as coefficient matrix.
C-----------------------------------------------------------------------
 350  DO 360 I = 1,N
 360    Y(I) = H*SAVF(I) - (YH(I,2) + ACOR(I))
      CALL SLVS (WM, IWM, Y, SAVF)
      IF (IERSL .LT. 0) GO TO 430
      IF (IERSL .GT. 0) GO TO 410
      DEL = DVNORM (N, Y, EWT)
      DO 380 I = 1,N
        ACOR(I) = ACOR(I) + Y(I)
 380    Y(I) = YH(I,1) + EL(1)*ACOR(I)
C-----------------------------------------------------------------------
C Test for convergence.  If M.gt.0, an estimate of the convergence
C rate constant is stored in CRATE, and this is used in the test.
C-----------------------------------------------------------------------
 400  IF (M .NE. 0) CRATE = MAX(0.2D0*CRATE,DEL/DELP)
      DCON = DEL*MIN(1.0D0,1.5D0*CRATE)/(TESCO(2,NQ)*CONIT)
      IF (DCON .LE. 1.0D0) GO TO 450
      M = M + 1
      IF (M .EQ. MAXCOR) GO TO 410
      IF (M .GE. 2 .AND. DEL .GT. 2.0D0*DELP) GO TO 410
      DELP = DEL
      CALL F (NEQ, TN, Y, SAVF)
      NFE = NFE + 1
      GO TO 270
C-----------------------------------------------------------------------
C The corrector iteration failed to converge.
C If MITER .ne. 0 and the Jacobian is out of date, PJAC is called for
C the next try.  Otherwise the YH array is retracted to its values
C before prediction, and H is reduced, if possible.  If H cannot be
C reduced or MXNCF failures have occurred, exit with KFLAG = -2.
C-----------------------------------------------------------------------
 410  IF (MITER .EQ. 0 .OR. JCUR .EQ. 1) GO TO 430
      ICF = 1
      IPUP = MITER
      GO TO 220
 430  ICF = 2
      NCF = NCF + 1
      RMAX = 2.0D0
      TN = TOLD
      I1 = NQNYH + 1
      DO 445 JB = 1,NQ
        I1 = I1 - NYH
Cdir$ ivdep
        DO 440 I = I1,NQNYH
 440      YH1(I) = YH1(I) - YH1(I+NYH)
 445    CONTINUE
      IF (IERPJ .LT. 0 .OR. IERSL .LT. 0) GO TO 680
      IF (ABS(H) .LE. HMIN*1.00001D0) GO TO 670
      IF (NCF .EQ. MXNCF) GO TO 670
      RH = 0.25D0
      IPUP = MITER
      IREDO = 1
      GO TO 170
C-----------------------------------------------------------------------
C The corrector has converged.  JCUR is set to 0
C to signal that the Jacobian involved may need updating later.
C The local error test is made and control passes to statement 500
C if it fails.
C-----------------------------------------------------------------------
 450  JCUR = 0
      IF (M .EQ. 0) DSM = DEL/TESCO(2,NQ)
      IF (M .GT. 0) DSM = DVNORM (N, ACOR, EWT)/TESCO(2,NQ)
      IF (DSM .GT. 1.0D0) GO TO 500
C-----------------------------------------------------------------------
C After a successful step, update the YH array.
C Consider changing H if IALTH = 1.  Otherwise decrease IALTH by 1.
C If IALTH is then 1 and NQ .lt. MAXORD, then ACOR is saved for
C use in a possible order increase on the next step.
C If a change in H is considered, an increase or decrease in order
C by one is considered also.  A change in H is made only if it is by a
C factor of at least 1.1.  If not, IALTH is set to 3 to prevent
C testing for that many steps.
C-----------------------------------------------------------------------
      KFLAG = 0
      IREDO = 0
      NST = NST + 1
      HU = H
      NQU = NQ
      DO 470 J = 1,L
        DO 470 I = 1,N
 470      YH(I,J) = YH(I,J) + EL(J)*ACOR(I)
      IALTH = IALTH - 1
      IF (IALTH .EQ. 0) GO TO 520
      IF (IALTH .GT. 1) GO TO 700
      IF (L .EQ. LMAX) GO TO 700
      DO 490 I = 1,N
 490    YH(I,LMAX) = ACOR(I)
      GO TO 700
C-----------------------------------------------------------------------
C The error test failed.  KFLAG keeps track of multiple failures.
C Restore TN and the YH array to their previous values, and prepare
C to try the step again.  Compute the optimum step size for this or
C one lower order.  After 2 or more failures, H is forced to decrease
C by a factor of 0.2 or less.
C-----------------------------------------------------------------------
 500  KFLAG = KFLAG - 1
      TN = TOLD
      I1 = NQNYH + 1
      DO 515 JB = 1,NQ
        I1 = I1 - NYH
Cdir$ ivdep
        DO 510 I = I1,NQNYH
 510      YH1(I) = YH1(I) - YH1(I+NYH)
 515    CONTINUE
      RMAX = 2.0D0
      IF (ABS(H) .LE. HMIN*1.00001D0) GO TO 660
      IF (KFLAG .LE. -3) GO TO 640
      IREDO = 2
      RHUP = 0.0D0
      GO TO 540
C-----------------------------------------------------------------------
C Regardless of the success or failure of the step, factors
C RHDN, RHSM, and RHUP are computed, by which H could be multiplied
C at order NQ - 1, order NQ, or order NQ + 1, respectively.
C In the case of failure, RHUP = 0.0 to avoid an order increase.
C The largest of these is determined and the new order chosen
C accordingly.  If the order is to be increased, we compute one
C additional scaled derivative.
C-----------------------------------------------------------------------
 520  RHUP = 0.0D0
      IF (L .EQ. LMAX) GO TO 540
      DO 530 I = 1,N
 530    SAVF(I) = ACOR(I) - YH(I,LMAX)
      DUP = DVNORM (N, SAVF, EWT)/TESCO(3,NQ)
      EXUP = 1.0D0/(L+1)
      RHUP = 1.0D0/(1.4D0*DUP**EXUP + 0.0000014D0)
 540  EXSM = 1.0D0/L
      RHSM = 1.0D0/(1.2D0*DSM**EXSM + 0.0000012D0)
      RHDN = 0.0D0
      IF (NQ .EQ. 1) GO TO 560
      DDN = DVNORM (N, YH(1,L), EWT)/TESCO(1,NQ)
      EXDN = 1.0D0/NQ
      RHDN = 1.0D0/(1.3D0*DDN**EXDN + 0.0000013D0)
 560  IF (RHSM .GE. RHUP) GO TO 570
      IF (RHUP .GT. RHDN) GO TO 590
      GO TO 580
 570  IF (RHSM .LT. RHDN) GO TO 580
      NEWQ = NQ
      RH = RHSM
      GO TO 620
 580  NEWQ = NQ - 1
      RH = RHDN
      IF (KFLAG .LT. 0 .AND. RH .GT. 1.0D0) RH = 1.0D0
      GO TO 620
 590  NEWQ = L
      RH = RHUP
      IF (RH .LT. 1.1D0) GO TO 610
      R = EL(L)/L
      DO 600 I = 1,N
 600    YH(I,NEWQ+1) = ACOR(I)*R
      GO TO 630
 610  IALTH = 3
      GO TO 700
 620  IF ((KFLAG .EQ. 0) .AND. (RH .LT. 1.1D0)) GO TO 610
      IF (KFLAG .LE. -2) RH = MIN(RH,0.2D0)
C-----------------------------------------------------------------------
C If there is a change of order, reset NQ, l, and the coefficients.
C In any case H is reset according to RH and the YH array is rescaled.
C Then exit from 690 if the step was OK, or redo the step otherwise.
C-----------------------------------------------------------------------
      IF (NEWQ .EQ. NQ) GO TO 170
 630  NQ = NEWQ
      L = NQ + 1
      IRET = 2
      GO TO 150
C-----------------------------------------------------------------------
C Control reaches this section if 3 or more failures have occured.
C If 10 failures have occurred, exit with KFLAG = -1.
C It is assumed that the derivatives that have accumulated in the
C YH array have errors of the wrong order.  Hence the first
C derivative is recomputed, and the order is set to 1.  Then
C H is reduced by a factor of 10, and the step is retried,
C until it succeeds or H reaches HMIN.
C-----------------------------------------------------------------------
 640  IF (KFLAG .EQ. -10) GO TO 660
      RH = 0.1D0
      RH = MAX(HMIN/ABS(H),RH)
      H = H*RH
      DO 645 I = 1,N
 645    Y(I) = YH(I,1)
      CALL F (NEQ, TN, Y, SAVF)
      NFE = NFE + 1
      DO 650 I = 1,N
 650    YH(I,2) = H*SAVF(I)
      IPUP = MITER
      IALTH = 5
      IF (NQ .EQ. 1) GO TO 200
      NQ = 1
      L = 2
      IRET = 3
      GO TO 150
C-----------------------------------------------------------------------
C All returns are made through this section.  H is saved in HOLD
C to allow the caller to change H on the next step.
C-----------------------------------------------------------------------
 660  KFLAG = -1
      GO TO 720
 670  KFLAG = -2
      GO TO 720
 680  KFLAG = -3
      GO TO 720
 690  RMAX = 10.0D0
 700  R = 1.0D0/TESCO(2,NQU)
      DO 710 I = 1,N
 710    ACOR(I) = ACOR(I)*R
 720  HOLD = H
      JSTART = 1
      RETURN
C----------------------- END OF SUBROUTINE DSTODE ----------------------
      END
*DECK DEWSET
      SUBROUTINE DEWSET (N, ITOL, RTOL, ATOL, YCUR, EWT)
C***BEGIN PROLOGUE  DEWSET
C***SUBSIDIARY
C***PURPOSE  Set error weight vector.
C***TYPE      DOUBLE PRECISION (SEWSET-S, DEWSET-D)
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C***DESCRIPTION
C
C  This subroutine sets the error weight vector EWT according to
C      EWT(i) = RTOL(i)*ABS(YCUR(i)) + ATOL(i),  i = 1,...,N,
C  with the subscript on RTOL and/or ATOL possibly replaced by 1 above,
C  depending on the value of ITOL.
C
C***SEE ALSO  DLSODE
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791129  DATE WRITTEN
C   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
C   890503  Minor cosmetic changes.  (FNF)
C   930809  Renamed to allow single/double precision versions. (ACH)
C***END PROLOGUE  DEWSET
C**End
      INTEGER N, ITOL
      INTEGER I
      DOUBLE PRECISION RTOL, ATOL, YCUR, EWT
      DIMENSION RTOL(*), ATOL(*), YCUR(N), EWT(N)
C
C***FIRST EXECUTABLE STATEMENT  DEWSET
      GO TO (10, 20, 30, 40), ITOL
 10   CONTINUE
      DO 15 I = 1,N
 15     EWT(I) = RTOL(1)*ABS(YCUR(I)) + ATOL(1)
      RETURN
 20   CONTINUE
      DO 25 I = 1,N
 25     EWT(I) = RTOL(1)*ABS(YCUR(I)) + ATOL(I)
      RETURN
 30   CONTINUE
      DO 35 I = 1,N
 35     EWT(I) = RTOL(I)*ABS(YCUR(I)) + ATOL(1)
      RETURN
 40   CONTINUE
      DO 45 I = 1,N
 45     EWT(I) = RTOL(I)*ABS(YCUR(I)) + ATOL(I)
      RETURN
C----------------------- END OF SUBROUTINE DEWSET ----------------------
      END
*DECK DVNORM
      DOUBLE PRECISION FUNCTION DVNORM (N, V, W)
C***BEGIN PROLOGUE  DVNORM
C***SUBSIDIARY
C***PURPOSE  Weighted root-mean-square vector norm.
C***TYPE      DOUBLE PRECISION (SVNORM-S, DVNORM-D)
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C***DESCRIPTION
C
C  This function routine computes the weighted root-mean-square norm
C  of the vector of length N contained in the array V, with weights
C  contained in the array W of length N:
C    DVNORM = SQRT( (1/N) * SUM( V(i)*W(i) )**2 )
C
C***SEE ALSO  DLSODE
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791129  DATE WRITTEN
C   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
C   890503  Minor cosmetic changes.  (FNF)
C   930809  Renamed to allow single/double precision versions. (ACH)
C***END PROLOGUE  DVNORM
C**End
      INTEGER N,   I
      DOUBLE PRECISION V, W,   SUM
      DIMENSION V(N), W(N)
C
C***FIRST EXECUTABLE STATEMENT  DVNORM
      SUM = 0.0D0
      DO 10 I = 1,N
 10     SUM = SUM + (V(I)*W(I))**2
      DVNORM = SQRT(SUM/N)
      RETURN
C----------------------- END OF FUNCTION DVNORM ------------------------
      END
*DECK DIPREP
      SUBROUTINE DIPREP (NEQ, Y, RWORK, IA, JA, IPFLAG, F, JAC)
      EXTERNAL F, JAC
      INTEGER NEQ, IA, JA, IPFLAG
      INTEGER ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, 
     1   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     2   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      INTEGER IPLOST, IESP, ISTATC, IYS, IBA, IBIAN, IBJAN, IBJGP,
     1   IPIAN, IPJAN, IPJGP, IPIGP, IPR, IPC, IPIC, IPISP, IPRSP, IPA,
     2   LENYH, LENYHM, LENWK, LREQ, LRAT, LREST, LWMIN, MOSS, MSBJ,
     3   NSLJ, NGP, NLU, NNZ, NSP, NZL, NZU
      INTEGER I, IMAX, LEWTN, LYHD, LYHN
      DOUBLE PRECISION Y, RWORK
      DOUBLE PRECISION CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
      DOUBLE PRECISION RLSS
      DIMENSION NEQ(*), Y(*), RWORK(*), IA(*), JA(*)
      COMMON /DLS001/ CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,
     1   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, 
     2   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     3   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      COMMON /DLSS01/ RLSS(6),
     1   IPLOST, IESP, ISTATC, IYS, IBA, IBIAN, IBJAN, IBJGP,
     2   IPIAN, IPJAN, IPJGP, IPIGP, IPR, IPC, IPIC, IPISP, IPRSP, IPA,
     3   LENYH, LENYHM, LENWK, LREQ, LRAT, LREST, LWMIN, MOSS, MSBJ,
     4   NSLJ, NGP, NLU, NNZ, NSP, NZL, NZU
C-----------------------------------------------------------------------
C This routine serves as an interface between the driver and
C Subroutine DPREP.  It is called only if MITER is 1 or 2.
C Tasks performed here are:
C  * call DPREP,
C  * reset the required WM segment length LENWK,
C  * move YH back to its final location (following WM in RWORK),
C  * reset pointers for YH, SAVF, EWT, and ACOR, and
C  * move EWT to its new position if ISTATE = 1.
C IPFLAG is an output error indication flag.  IPFLAG = 0 if there was
C no trouble, and IPFLAG is the value of the DPREP error flag IPPER
C if there was trouble in Subroutine DPREP.
C-----------------------------------------------------------------------
      IPFLAG = 0
C Call DPREP to do matrix preprocessing operations. --------------------
      CALL DPREP (NEQ, Y, RWORK(LYH), RWORK(LSAVF), RWORK(LEWT),
     1   RWORK(LACOR), IA, JA, RWORK(LWM), RWORK(LWM), IPFLAG, F, JAC)
      LENWK = MAX(LREQ,LWMIN)
      IF (IPFLAG .LT. 0) RETURN
C If DPREP was successful, move YH to end of required space for WM. ----
      LYHN = LWM + LENWK
      IF (LYHN .GT. LYH) RETURN
      LYHD = LYH - LYHN
      IF (LYHD .EQ. 0) GO TO 20
      IMAX = LYHN - 1 + LENYHM
      DO 10 I = LYHN,IMAX
 10     RWORK(I) = RWORK(I+LYHD)
      LYH = LYHN
C Reset pointers for SAVF, EWT, and ACOR. ------------------------------
 20   LSAVF = LYH + LENYH
      LEWTN = LSAVF + N
      LACOR = LEWTN + N
      IF (ISTATC .EQ. 3) GO TO 40
C If ISTATE = 1, move EWT (left) to its new position. ------------------
      IF (LEWTN .GT. LEWT) RETURN
      DO 30 I = 1,N
 30     RWORK(I+LEWTN-1) = RWORK(I+LEWT-1)
 40   LEWT = LEWTN
      RETURN
C----------------------- End of Subroutine DIPREP ----------------------
      END
*DECK DPREP
      SUBROUTINE DPREP (NEQ, Y, YH, SAVF, EWT, FTEM, IA, JA,
     1                     WK, IWK, IPPER, F, JAC)
      EXTERNAL F,JAC
      INTEGER NEQ, IA, JA, IWK, IPPER
      INTEGER ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, 
     1   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     2   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      INTEGER IPLOST, IESP, ISTATC, IYS, IBA, IBIAN, IBJAN, IBJGP,
     1   IPIAN, IPJAN, IPJGP, IPIGP, IPR, IPC, IPIC, IPISP, IPRSP, IPA,
     2   LENYH, LENYHM, LENWK, LREQ, LRAT, LREST, LWMIN, MOSS, MSBJ,
     3   NSLJ, NGP, NLU, NNZ, NSP, NZL, NZU
      INTEGER I, IBR, IER, IPIL, IPIU, IPTT1, IPTT2, J, JFOUND, K,
     1   KNEW, KMAX, KMIN, LDIF, LENIGP, LIWK, MAXG, NP1, NZSUT
      DOUBLE PRECISION Y, YH, SAVF, EWT, FTEM, WK
      DOUBLE PRECISION CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
      DOUBLE PRECISION CON0, CONMIN, CCMXJ, PSMALL, RBIG, SETH
      DOUBLE PRECISION DQ, DYJ, ERWT, FAC, YJ
      DIMENSION NEQ(*), Y(*), YH(*), SAVF(*), EWT(*), FTEM(*),
     1   IA(*), JA(*), WK(*), IWK(*)
      COMMON /DLS001/ CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,
     1   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, 
     2   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     3   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      COMMON /DLSS01/ CON0, CONMIN, CCMXJ, PSMALL, RBIG, SETH,
     1   IPLOST, IESP, ISTATC, IYS, IBA, IBIAN, IBJAN, IBJGP,
     2   IPIAN, IPJAN, IPJGP, IPIGP, IPR, IPC, IPIC, IPISP, IPRSP, IPA,
     3   LENYH, LENYHM, LENWK, LREQ, LRAT, LREST, LWMIN, MOSS, MSBJ,
     4   NSLJ, NGP, NLU, NNZ, NSP, NZL, NZU
C-----------------------------------------------------------------------
C This routine performs preprocessing related to the sparse linear
C systems that must be solved if MITER = 1 or 2.
C The operations that are performed here are:
C  * compute sparseness structure of Jacobian according to MOSS,
C  * compute grouping of column indices (MITER = 2),
C  * compute a new ordering of rows and columns of the matrix,
C  * reorder JA corresponding to the new ordering,
C  * perform a symbolic LU factorization of the matrix, and
C  * set pointers for segments of the IWK/WK array.
C In addition to variables described previously, DPREP uses the
C following for communication:
C YH     = the history array.  Only the first column, containing the
C          current Y vector, is used.  Used only if MOSS .ne. 0.
C SAVF   = a work array of length NEQ, used only if MOSS .ne. 0.
C EWT    = array of length NEQ containing (inverted) error weights.
C          Used only if MOSS = 2 or if ISTATE = MOSS = 1.
C FTEM   = a work array of length NEQ, identical to ACOR in the driver,
C          used only if MOSS = 2.
C WK     = a real work array of length LENWK, identical to WM in
C          the driver.
C IWK    = integer work array, assumed to occupy the same space as WK.
C LENWK  = the length of the work arrays WK and IWK.
C ISTATC = a copy of the driver input argument ISTATE (= 1 on the
C          first call, = 3 on a continuation call).
C IYS    = flag value from ODRV or CDRV.
C IPPER  = output error flag with the following values and meanings:
C          0  no error.
C         -1  insufficient storage for internal structure pointers.
C         -2  insufficient storage for JGROUP.
C         -3  insufficient storage for ODRV.
C         -4  other error flag from ODRV (should never occur).
C         -5  insufficient storage for CDRV.
C         -6  other error flag from CDRV.
C-----------------------------------------------------------------------
      IBIAN = LRAT*2
      IPIAN = IBIAN + 1
      NP1 = N + 1
      IPJAN = IPIAN + NP1
      IBJAN = IPJAN - 1
      LIWK = LENWK*LRAT
      IF (IPJAN+N-1 .GT. LIWK) GO TO 210
      IF (MOSS .EQ. 0) GO TO 30
C
      IF (ISTATC .EQ. 3) GO TO 20
C ISTATE = 1 and MOSS .ne. 0.  Perturb Y for structure determination. --
      DO 10 I = 1,N
        ERWT = 1.0D0/EWT(I)
        FAC = 1.0D0 + 1.0D0/(I + 1.0D0)
        Y(I) = Y(I) + FAC*SIGN(ERWT,Y(I))
 10     CONTINUE
      GO TO (70, 100), MOSS
C
 20   CONTINUE
C ISTATE = 3 and MOSS .ne. 0.  Load Y from YH(*,1). --------------------
      DO 25 I = 1,N
 25     Y(I) = YH(I)
      GO TO (70, 100), MOSS
C
C MOSS = 0.  Process user's IA,JA.  Add diagonal entries if necessary. -
 30   KNEW = IPJAN
      KMIN = IA(1)
      IWK(IPIAN) = 1
      DO 60 J = 1,N
        JFOUND = 0
        KMAX = IA(J+1) - 1
        IF (KMIN .GT. KMAX) GO TO 45
        DO 40 K = KMIN,KMAX
          I = JA(K)
          IF (I .EQ. J) JFOUND = 1
          IF (KNEW .GT. LIWK) GO TO 210
          IWK(KNEW) = I
          KNEW = KNEW + 1
 40       CONTINUE
        IF (JFOUND .EQ. 1) GO TO 50
 45     IF (KNEW .GT. LIWK) GO TO 210
        IWK(KNEW) = J
        KNEW = KNEW + 1
 50     IWK(IPIAN+J) = KNEW + 1 - IPJAN
        KMIN = KMAX + 1
 60     CONTINUE
      GO TO 140
C
C MOSS = 1.  Compute structure from user-supplied Jacobian routine JAC.
 70   CONTINUE
C A dummy call to F allows user to create temporaries for use in JAC. --
      CALL F (NEQ, TN, Y, SAVF)
      K = IPJAN
      IWK(IPIAN) = 1
      DO 90 J = 1,N
        IF (K .GT. LIWK) GO TO 210
        IWK(K) = J
        K = K + 1
        DO 75 I = 1,N
 75       SAVF(I) = 0.0D0
        CALL JAC (NEQ, TN, Y, J, IWK(IPIAN), IWK(IPJAN), SAVF)
        DO 80 I = 1,N
          IF (ABS(SAVF(I)) .LE. SETH) GO TO 80
          IF (I .EQ. J) GO TO 80
          IF (K .GT. LIWK) GO TO 210
          IWK(K) = I
          K = K + 1
 80       CONTINUE
        IWK(IPIAN+J) = K + 1 - IPJAN
 90     CONTINUE
      GO TO 140
C
C MOSS = 2.  Compute structure from results of N + 1 calls to F. -------
 100  K = IPJAN
      IWK(IPIAN) = 1
      CALL F (NEQ, TN, Y, SAVF)
      DO 120 J = 1,N
        IF (K .GT. LIWK) GO TO 210
        IWK(K) = J
        K = K + 1
        YJ = Y(J)
        ERWT = 1.0D0/EWT(J)
        DYJ = SIGN(ERWT,YJ)
        Y(J) = YJ + DYJ
        CALL F (NEQ, TN, Y, FTEM)
        Y(J) = YJ
        DO 110 I = 1,N
          DQ = (FTEM(I) - SAVF(I))/DYJ
          IF (ABS(DQ) .LE. SETH) GO TO 110
          IF (I .EQ. J) GO TO 110
          IF (K .GT. LIWK) GO TO 210
          IWK(K) = I
          K = K + 1
 110      CONTINUE
        IWK(IPIAN+J) = K + 1 - IPJAN
 120    CONTINUE
C
 140  CONTINUE
      IF (MOSS .EQ. 0 .OR. ISTATC .NE. 1) GO TO 150
C If ISTATE = 1 and MOSS .ne. 0, restore Y from YH. --------------------
      DO 145 I = 1,N
 145    Y(I) = YH(I)
 150  NNZ = IWK(IPIAN+N) - 1
      LENIGP = 0
      IPIGP = IPJAN + NNZ
      IF (MITER .NE. 2) GO TO 160
C
C Compute grouping of column indices (MITER = 2). ----------------------
      MAXG = NP1
      IPJGP = IPJAN + NNZ
      IBJGP = IPJGP - 1
      IPIGP = IPJGP + N
      IPTT1 = IPIGP + NP1
      IPTT2 = IPTT1 + N
      LREQ = IPTT2 + N - 1
      IF (LREQ .GT. LIWK) GO TO 220
      CALL JGROUP (N, IWK(IPIAN), IWK(IPJAN), MAXG, NGP, IWK(IPIGP),
     1   IWK(IPJGP), IWK(IPTT1), IWK(IPTT2), IER)
      IF (IER .NE. 0) GO TO 220
      LENIGP = NGP + 1
C
C Compute new ordering of rows/columns of Jacobian. --------------------
 160  IPR = IPIGP + LENIGP
      IPC = IPR
      IPIC = IPC + N
      IPISP = IPIC + N
      IPRSP = (IPISP - 2)/LRAT + 2
      IESP = LENWK + 1 - IPRSP
      IF (IESP .LT. 0) GO TO 230
      IBR = IPR - 1
      DO 170 I = 1,N
 170    IWK(IBR+I) = I
      NSP = LIWK + 1 - IPISP
      CALL ODRV (N, IWK(IPIAN), IWK(IPJAN), WK, IWK(IPR), IWK(IPIC),
     1   NSP, IWK(IPISP), 1, IYS)
      IF (IYS .EQ. 11*N+1) GO TO 240
      IF (IYS .NE. 0) GO TO 230
C
C Reorder JAN and do symbolic LU factorization of matrix. --------------
      IPA = LENWK + 1 - NNZ
      NSP = IPA - IPRSP
      LREQ = MAX(12*N/LRAT, 6*N/LRAT+2*N+NNZ) + 3
      LREQ = LREQ + IPRSP - 1 + NNZ
      IF (LREQ .GT. LENWK) GO TO 250
      IBA = IPA - 1
      DO 180 I = 1,NNZ
 180    WK(IBA+I) = 0.0D0
      IPISP = LRAT*(IPRSP - 1) + 1
      CALL CDRV (N,IWK(IPR),IWK(IPC),IWK(IPIC),IWK(IPIAN),IWK(IPJAN),
     1   WK(IPA),WK(IPA),WK(IPA),NSP,IWK(IPISP),WK(IPRSP),IESP,5,IYS)
      LREQ = LENWK - IESP
      IF (IYS .EQ. 10*N+1) GO TO 250
      IF (IYS .NE. 0) GO TO 260
      IPIL = IPISP
      IPIU = IPIL + 2*N + 1
      NZU = IWK(IPIL+N) - IWK(IPIL)
      NZL = IWK(IPIU+N) - IWK(IPIU)
      IF (LRAT .GT. 1) GO TO 190
      CALL ADJLR (N, IWK(IPISP), LDIF)
      LREQ = LREQ + LDIF
 190  CONTINUE
      IF (LRAT .EQ. 2 .AND. NNZ .EQ. N) LREQ = LREQ + 1
      NSP = NSP + LREQ - LENWK
      IPA = LREQ + 1 - NNZ
      IBA = IPA - 1
      IPPER = 0
      RETURN
C
 210  IPPER = -1
      LREQ = 2 + (2*N + 1)/LRAT
      LREQ = MAX(LENWK+1,LREQ)
      RETURN
C
 220  IPPER = -2
      LREQ = (LREQ - 1)/LRAT + 1
      RETURN
C
 230  IPPER = -3
      CALL CNTNZU (N, IWK(IPIAN), IWK(IPJAN), NZSUT)
      LREQ = LENWK - IESP + (3*N + 4*NZSUT - 1)/LRAT + 1
      RETURN
C
 240  IPPER = -4
      RETURN
C
 250  IPPER = -5
      RETURN
C
 260  IPPER = -6
      LREQ = LENWK
      RETURN
C----------------------- End of Subroutine DPREP -----------------------
      END
*DECK JGROUP
      SUBROUTINE JGROUP (N,IA,JA,MAXG,NGRP,IGP,JGP,INCL,JDONE,IER)
      INTEGER N, IA, JA, MAXG, NGRP, IGP, JGP, INCL, JDONE, IER
      DIMENSION IA(*), JA(*), IGP(*), JGP(*), INCL(*), JDONE(*)
C-----------------------------------------------------------------------
C This subroutine constructs groupings of the column indices of
C the Jacobian matrix, used in the numerical evaluation of the
C Jacobian by finite differences.
C
C Input:
C N      = the order of the matrix.
C IA,JA  = sparse structure descriptors of the matrix by rows.
C MAXG   = length of available storage in the IGP array.
C
C Output:
C NGRP   = number of groups.
C JGP    = array of length N containing the column indices by groups.
C IGP    = pointer array of length NGRP + 1 to the locations in JGP
C          of the beginning of each group.
C IER    = error indicator.  IER = 0 if no error occurred, or 1 if
C          MAXG was insufficient.
C
C INCL and JDONE are working arrays of length N.
C-----------------------------------------------------------------------
      INTEGER I, J, K, KMIN, KMAX, NCOL, NG
C
      IER = 0
      DO 10 J = 1,N
 10     JDONE(J) = 0
      NCOL = 1
      DO 60 NG = 1,MAXG
        IGP(NG) = NCOL
        DO 20 I = 1,N
 20       INCL(I) = 0
        DO 50 J = 1,N
C Reject column J if it is already in a group.--------------------------
          IF (JDONE(J) .EQ. 1) GO TO 50
          KMIN = IA(J)
          KMAX = IA(J+1) - 1
          DO 30 K = KMIN,KMAX
C Reject column J if it overlaps any column already in this group.------
            I = JA(K)
            IF (INCL(I) .EQ. 1) GO TO 50
 30         CONTINUE
C Accept column J into group NG.----------------------------------------
          JGP(NCOL) = J
          NCOL = NCOL + 1
          JDONE(J) = 1
          DO 40 K = KMIN,KMAX
            I = JA(K)
 40         INCL(I) = 1
 50       CONTINUE
C Stop if this group is empty (grouping is complete).-------------------
        IF (NCOL .EQ. IGP(NG)) GO TO 70
 60     CONTINUE
C Error return if not all columns were chosen (MAXG too small).---------
      IF (NCOL .LE. N) GO TO 80
      NG = MAXG
 70   NGRP = NG - 1
      RETURN
 80   IER = 1
      RETURN
C----------------------- End of Subroutine JGROUP ----------------------
      END
*DECK ADJLR
      SUBROUTINE ADJLR (N, ISP, LDIF)
      INTEGER N, ISP, LDIF
      DIMENSION ISP(*)
C-----------------------------------------------------------------------
C This routine computes an adjustment, LDIF, to the required
C integer storage space in IWK (sparse matrix work space).
C It is called only if the word length ratio is LRAT = 1.
C This is to account for the possibility that the symbolic LU phase
C may require more storage than the numerical LU and solution phases.
C-----------------------------------------------------------------------
      INTEGER IP, JLMAX, JUMAX, LNFC, LSFC, NZLU
C
      IP = 2*N + 1
C Get JLMAX = IJL(N) and JUMAX = IJU(N) (sizes of JL and JU). ----------
      JLMAX = ISP(IP)
      JUMAX = ISP(IP+IP)
C NZLU = (size of L) + (size of U) = (IL(N+1)-IL(1)) + (IU(N+1)-IU(1)).
      NZLU = ISP(N+1) - ISP(1) + ISP(IP+N+1) - ISP(IP+1)
      LSFC = 12*N + 3 + 2*MAX(JLMAX,JUMAX)
      LNFC = 9*N + 2 + JLMAX + JUMAX + NZLU
      LDIF = MAX(0, LSFC - LNFC)
      RETURN
C----------------------- End of Subroutine ADJLR -----------------------
      END
*DECK CNTNZU
      SUBROUTINE CNTNZU (N, IA, JA, NZSUT)
      INTEGER N, IA, JA, NZSUT
      DIMENSION IA(*), JA(*)
C-----------------------------------------------------------------------
C This routine counts the number of nonzero elements in the strict
C upper triangle of the matrix M + M(transpose), where the sparsity
C structure of M is given by pointer arrays IA and JA.
C This is needed to compute the storage requirements for the
C sparse matrix reordering operation in ODRV.
C-----------------------------------------------------------------------
      INTEGER II, JJ, J, JMIN, JMAX, K, KMIN, KMAX, NUM
C
      NUM = 0
      DO 50 II = 1,N
        JMIN = IA(II)
        JMAX = IA(II+1) - 1
        IF (JMIN .GT. JMAX) GO TO 50
        DO 40 J = JMIN,JMAX
          IF (JA(J) - II) 10, 40, 30
 10       JJ =JA(J)
          KMIN = IA(JJ)
          KMAX = IA(JJ+1) - 1
          IF (KMIN .GT. KMAX) GO TO 30
          DO 20 K = KMIN,KMAX
            IF (JA(K) .EQ. II) GO TO 40
 20         CONTINUE
 30       NUM = NUM + 1
 40       CONTINUE
 50     CONTINUE
      NZSUT = NUM
      RETURN
C----------------------- End of Subroutine CNTNZU ----------------------
      END
*DECK DPRJS
      SUBROUTINE DPRJS (NEQ,Y,YH,NYH,EWT,FTEM,SAVF,WK,IWK,F,JAC)
      EXTERNAL F,JAC
      INTEGER NEQ, NYH, IWK
      INTEGER ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, 
     1   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     2   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      INTEGER IPLOST, IESP, ISTATC, IYS, IBA, IBIAN, IBJAN, IBJGP,
     1   IPIAN, IPJAN, IPJGP, IPIGP, IPR, IPC, IPIC, IPISP, IPRSP, IPA,
     2   LENYH, LENYHM, LENWK, LREQ, LRAT, LREST, LWMIN, MOSS, MSBJ,
     3   NSLJ, NGP, NLU, NNZ, NSP, NZL, NZU
      INTEGER I, IMUL, J, JJ, JOK, JMAX, JMIN, K, KMAX, KMIN, NG
      DOUBLE PRECISION Y, YH, EWT, FTEM, SAVF, WK
      DOUBLE PRECISION CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
      DOUBLE PRECISION CON0, CONMIN, CCMXJ, PSMALL, RBIG, SETH
      DOUBLE PRECISION CON, DI, FAC, HL0, PIJ, R, R0, RCON, RCONT,
     1   SRUR, DVNORM
      DIMENSION NEQ(*), Y(*), YH(NYH,*), EWT(*), FTEM(*), SAVF(*),
     1   WK(*), IWK(*)
      COMMON /DLS001/ CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,
     1   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, 
     2   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     3   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      COMMON /DLSS01/ CON0, CONMIN, CCMXJ, PSMALL, RBIG, SETH,
     1   IPLOST, IESP, ISTATC, IYS, IBA, IBIAN, IBJAN, IBJGP,
     2   IPIAN, IPJAN, IPJGP, IPIGP, IPR, IPC, IPIC, IPISP, IPRSP, IPA,
     3   LENYH, LENYHM, LENWK, LREQ, LRAT, LREST, LWMIN, MOSS, MSBJ,
     4   NSLJ, NGP, NLU, NNZ, NSP, NZL, NZU
C-----------------------------------------------------------------------
C DPRJS is called to compute and process the matrix
C P = I - H*EL(1)*J , where J is an approximation to the Jacobian.
C J is computed by columns, either by the user-supplied routine JAC
C if MITER = 1, or by finite differencing if MITER = 2.
C if MITER = 3, a diagonal approximation to J is used.
C if MITER = 1 or 2, and if the existing value of the Jacobian
C (as contained in P) is considered acceptable, then a new value of
C P is reconstructed from the old value.  In any case, when MITER
C is 1 or 2, the P matrix is subjected to LU decomposition in CDRV.
C P and its LU decomposition are stored (separately) in WK.
C
C In addition to variables described previously, communication
C with DPRJS uses the following:
C Y     = array containing predicted values on entry.
C FTEM  = work array of length N (ACOR in DSTODE).
C SAVF  = array containing f evaluated at predicted y.
C WK    = real work space for matrices.  On output it contains the
C         inverse diagonal matrix if MITER = 3, and P and its sparse
C         LU decomposition if MITER is 1 or 2.
C         Storage of matrix elements starts at WK(3).
C         WK also contains the following matrix-related data:
C         WK(1) = SQRT(UROUND), used in numerical Jacobian increments.
C         WK(2) = H*EL0, saved for later use if MITER = 3.
C IWK   = integer work space for matrix-related data, assumed to
C         be equivalenced to WK.  In addition, WK(IPRSP) and IWK(IPISP)
C         are assumed to have identical locations.
C EL0   = EL(1) (input).
C IERPJ = output error flag (in Common).
C       = 0 if no error.
C       = 1  if zero pivot found in CDRV.
C       = 2  if a singular matrix arose with MITER = 3.
C       = -1 if insufficient storage for CDRV (should not occur here).
C       = -2 if other error found in CDRV (should not occur here).
C JCUR  = output flag showing status of (approximate) Jacobian matrix:
C          = 1 to indicate that the Jacobian is now current, or
C          = 0 to indicate that a saved value was used.
C This routine also uses other variables in Common.
C-----------------------------------------------------------------------
      HL0 = H*EL0
      CON = -HL0
      IF (MITER .EQ. 3) GO TO 300
C See whether J should be reevaluated (JOK = 0) or not (JOK = 1). ------
      JOK = 1
      IF (NST .EQ. 0 .OR. NST .GE. NSLJ+MSBJ) JOK = 0
      IF (ICF .EQ. 1 .AND. ABS(RC - 1.0D0) .LT. CCMXJ) JOK = 0
      IF (ICF .EQ. 2) JOK = 0
      IF (JOK .EQ. 1) GO TO 250
C
C MITER = 1 or 2, and the Jacobian is to be reevaluated. ---------------
 20   JCUR = 1
      NJE = NJE + 1
      NSLJ = NST
      IPLOST = 0
      CONMIN = ABS(CON)
      GO TO (100, 200), MITER
C
C If MITER = 1, call JAC, multiply by scalar, and add identity. --------
 100  CONTINUE
      KMIN = IWK(IPIAN)
      DO 130 J = 1, N
        KMAX = IWK(IPIAN+J) - 1
        DO 110 I = 1,N
 110      FTEM(I) = 0.0D0
        CALL JAC (NEQ, TN, Y, J, IWK(IPIAN), IWK(IPJAN), FTEM)
        DO 120 K = KMIN, KMAX
          I = IWK(IBJAN+K)
          WK(IBA+K) = FTEM(I)*CON
          IF (I .EQ. J) WK(IBA+K) = WK(IBA+K) + 1.0D0
 120      CONTINUE
        KMIN = KMAX + 1
 130    CONTINUE
      GO TO 290
C
C If MITER = 2, make NGP calls to F to approximate J and P. ------------
 200  CONTINUE
      FAC = DVNORM(N, SAVF, EWT)
      R0 = 1000.0D0 * ABS(H) * UROUND * N * FAC
      IF (R0 .EQ. 0.0D0) R0 = 1.0D0
      SRUR = WK(1)
      JMIN = IWK(IPIGP)
      DO 240 NG = 1,NGP
        JMAX = IWK(IPIGP+NG) - 1
        DO 210 J = JMIN,JMAX
          JJ = IWK(IBJGP+J)
          R = MAX(SRUR*ABS(Y(JJ)),R0/EWT(JJ))
 210      Y(JJ) = Y(JJ) + R
        CALL F (NEQ, TN, Y, FTEM)
        DO 230 J = JMIN,JMAX
          JJ = IWK(IBJGP+J)
          Y(JJ) = YH(JJ,1)
          R = MAX(SRUR*ABS(Y(JJ)),R0/EWT(JJ))
          FAC = -HL0/R
          KMIN =IWK(IBIAN+JJ)
          KMAX =IWK(IBIAN+JJ+1) - 1
          DO 220 K = KMIN,KMAX
            I = IWK(IBJAN+K)
            WK(IBA+K) = (FTEM(I) - SAVF(I))*FAC
            IF (I .EQ. JJ) WK(IBA+K) = WK(IBA+K) + 1.0D0
 220        CONTINUE
 230      CONTINUE
        JMIN = JMAX + 1
 240    CONTINUE
      NFE = NFE + NGP
      GO TO 290
C
C If JOK = 1, reconstruct new P from old P. ----------------------------
 250  JCUR = 0
      RCON = CON/CON0
      RCONT = ABS(CON)/CONMIN
      IF (RCONT .GT. RBIG .AND. IPLOST .EQ. 1) GO TO 20
      KMIN = IWK(IPIAN)
      DO 275 J = 1,N
        KMAX = IWK(IPIAN+J) - 1
        DO 270 K = KMIN,KMAX
          I = IWK(IBJAN+K)
          PIJ = WK(IBA+K)
          IF (I .NE. J) GO TO 260
          PIJ = PIJ - 1.0D0
          IF (ABS(PIJ) .GE. PSMALL) GO TO 260
            IPLOST = 1
            CONMIN = MIN(ABS(CON0),CONMIN)
 260      PIJ = PIJ*RCON
          IF (I .EQ. J) PIJ = PIJ + 1.0D0
          WK(IBA+K) = PIJ
 270      CONTINUE
        KMIN = KMAX + 1
 275    CONTINUE
C
C Do numerical factorization of P matrix. ------------------------------
 290  NLU = NLU + 1
      CON0 = CON
      IERPJ = 0
      DO 295 I = 1,N
 295    FTEM(I) = 0.0D0
      CALL CDRV (N,IWK(IPR),IWK(IPC),IWK(IPIC),IWK(IPIAN),IWK(IPJAN),
     1   WK(IPA),FTEM,FTEM,NSP,IWK(IPISP),WK(IPRSP),IESP,2,IYS)
      IF (IYS .EQ. 0) RETURN
      IMUL = (IYS - 1)/N
      IERPJ = -2
      IF (IMUL .EQ. 8) IERPJ = 1
      IF (IMUL .EQ. 10) IERPJ = -1
      RETURN
C
C If MITER = 3, construct a diagonal approximation to J and P. ---------
 300  CONTINUE
      JCUR = 1
      NJE = NJE + 1
      WK(2) = HL0
      IERPJ = 0
      R = EL0*0.1D0
      DO 310 I = 1,N
 310    Y(I) = Y(I) + R*(H*SAVF(I) - YH(I,2))
      CALL F (NEQ, TN, Y, WK(3))
      NFE = NFE + 1
      DO 320 I = 1,N
        R0 = H*SAVF(I) - YH(I,2)
        DI = 0.1D0*R0 - H*(WK(I+2) - SAVF(I))
        WK(I+2) = 1.0D0
        IF (ABS(R0) .LT. UROUND/EWT(I)) GO TO 320
        IF (ABS(DI) .EQ. 0.0D0) GO TO 330
        WK(I+2) = 0.1D0*R0/DI
 320    CONTINUE
      RETURN
 330  IERPJ = 2
      RETURN
C----------------------- End of Subroutine DPRJS -----------------------
      END
*DECK DSOLSS
      SUBROUTINE DSOLSS (WK, IWK, X, TEM)
      INTEGER IWK
      INTEGER ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, 
     1   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     2   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      INTEGER IPLOST, IESP, ISTATC, IYS, IBA, IBIAN, IBJAN, IBJGP,
     1   IPIAN, IPJAN, IPJGP, IPIGP, IPR, IPC, IPIC, IPISP, IPRSP, IPA,
     2   LENYH, LENYHM, LENWK, LREQ, LRAT, LREST, LWMIN, MOSS, MSBJ,
     3   NSLJ, NGP, NLU, NNZ, NSP, NZL, NZU
      INTEGER I
      DOUBLE PRECISION WK, X, TEM
      DOUBLE PRECISION CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
      DOUBLE PRECISION RLSS
      DOUBLE PRECISION DI, HL0, PHL0, R
      DIMENSION WK(*), IWK(*), X(*), TEM(*)
      COMMON /DLS001/ CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,
     1   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, 
     2   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     3   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      COMMON /DLSS01/ RLSS(6),
     1   IPLOST, IESP, ISTATC, IYS, IBA, IBIAN, IBJAN, IBJGP,
     2   IPIAN, IPJAN, IPJGP, IPIGP, IPR, IPC, IPIC, IPISP, IPRSP, IPA,
     3   LENYH, LENYHM, LENWK, LREQ, LRAT, LREST, LWMIN, MOSS, MSBJ,
     4   NSLJ, NGP, NLU, NNZ, NSP, NZL, NZU
C-----------------------------------------------------------------------
C This routine manages the solution of the linear system arising from
C a chord iteration.  It is called if MITER .ne. 0.
C If MITER is 1 or 2, it calls CDRV to accomplish this.
C If MITER = 3 it updates the coefficient H*EL0 in the diagonal
C matrix, and then computes the solution.
C communication with DSOLSS uses the following variables:
C WK    = real work space containing the inverse diagonal matrix if
C         MITER = 3 and the LU decomposition of the matrix otherwise.
C         Storage of matrix elements starts at WK(3).
C         WK also contains the following matrix-related data:
C         WK(1) = SQRT(UROUND) (not used here),
C         WK(2) = HL0, the previous value of H*EL0, used if MITER = 3.
C IWK   = integer work space for matrix-related data, assumed to
C         be equivalenced to WK.  In addition, WK(IPRSP) and IWK(IPISP)
C         are assumed to have identical locations.
C X     = the right-hand side vector on input, and the solution vector
C         on output, of length N.
C TEM   = vector of work space of length N, not used in this version.
C IERSL = output flag (in Common).
C         IERSL = 0  if no trouble occurred.
C         IERSL = -1 if CDRV returned an error flag (MITER = 1 or 2).
C                    This should never occur and is considered fatal.
C         IERSL = 1  if a singular matrix arose with MITER = 3.
C This routine also uses other variables in Common.
C-----------------------------------------------------------------------
      IERSL = 0
      GO TO (100, 100, 300), MITER
 100  CALL CDRV (N,IWK(IPR),IWK(IPC),IWK(IPIC),IWK(IPIAN),IWK(IPJAN),
     1   WK(IPA),X,X,NSP,IWK(IPISP),WK(IPRSP),IESP,4,IERSL)
      IF (IERSL .NE. 0) IERSL = -1
      RETURN
C
 300  PHL0 = WK(2)
      HL0 = H*EL0
      WK(2) = HL0
      IF (HL0 .EQ. PHL0) GO TO 330
      R = HL0/PHL0
      DO 320 I = 1,N
        DI = 1.0D0 - R*(1.0D0 - 1.0D0/WK(I+2))
        IF (ABS(DI) .EQ. 0.0D0) GO TO 390
 320    WK(I+2) = 1.0D0/DI
 330  DO 340 I = 1,N
 340    X(I) = WK(I+2)*X(I)
      RETURN
 390  IERSL = 1
      RETURN
C
C----------------------- End of Subroutine DSOLSS ----------------------
      END
*DECK DSRCMS
      SUBROUTINE DSRCMS (RSAV, ISAV, JOB)
C-----------------------------------------------------------------------
C This routine saves or restores (depending on JOB) the contents of
C the Common blocks DLS001, DLSS01, which are used
C internally by one or more ODEPACK solvers.
C
C RSAV = real array of length 15 or more.
C ISAV = integer array of length 59 or more.
C JOB  = flag indicating to save or restore the Common blocks:
C        JOB  = 1 if Common is to be saved (written to RSAV/ISAV)
C        JOB  = 2 if Common is to be restored (read from RSAV/ISAV)
C        A call with JOB = 2 presumes a prior call with JOB = 1.
C-----------------------------------------------------------------------
      INTEGER ISAV, JOB
      INTEGER ILS, ILSS
      INTEGER I, LENIL, LENISS, LENRL, LENRSS
      DOUBLE PRECISION RSAV,   RLS, RLSS
      DIMENSION RSAV(*), ISAV(*)
      COMMON /DLS001/ RLS(9), ILS(25)
      COMMON /DLSS01/ RLSS(6), ILSS(34)
      DATA LENRL/9/, LENIL/25/, LENRSS/6/, LENISS/34/
C
      IF (JOB .EQ. 2) GO TO 100
      DO 10 I = 1,LENRL
 10     RSAV(I) = RLS(I)
      DO 15 I = 1,LENRSS
 15     RSAV(LENRL+I) = RLSS(I)
C
      DO 20 I = 1,LENIL
 20     ISAV(I) = ILS(I)
      DO 25 I = 1,LENISS
 25     ISAV(LENIL+I) = ILSS(I)
C
      RETURN
C
 100  CONTINUE
      DO 110 I = 1,LENRL
 110     RLS(I) = RSAV(I)
      DO 115 I = 1,LENRSS
 115     RLSS(I) = RSAV(LENRL+I)
C
      DO 120 I = 1,LENIL
 120     ILS(I) = ISAV(I)
      DO 125 I = 1,LENISS
 125     ILSS(I) = ISAV(LENIL+I)
C
      RETURN
C----------------------- End of Subroutine DSRCMS ----------------------
      END
*DECK ODRV
      subroutine odrv
     *     (n, ia,ja,a, p,ip, nsp,isp, path, flag)
c                                                                 5/2/83
c***********************************************************************
c  odrv -- driver for sparse matrix reordering routines
c***********************************************************************
c
c  description
c
c    odrv finds a minimum degree ordering of the rows and columns
c    of a matrix m stored in (ia,ja,a) format (see below).  for the
c    reordered matrix, the work and storage required to perform
c    gaussian elimination is (usually) significantly less.
c
c    note.. odrv and its subordinate routines have been modified to
c    compute orderings for general matrices, not necessarily having any
c    symmetry.  the miminum degree ordering is computed for the
c    structure of the symmetric matrix  m + m-transpose.
c    modifications to the original odrv module have been made in
c    the coding in subroutine mdi, and in the initial comments in
c    subroutines odrv and md.
c
c    if only the nonzero entries in the upper triangle of m are being
c    stored, then odrv symmetrically reorders (ia,ja,a), (optionally)
c    with the diagonal entries placed first in each row.  this is to
c    ensure that if m(i,j) will be in the upper triangle of m with
c    respect to the new ordering, then m(i,j) is stored in row i (and
c    thus m(j,i) is not stored),  whereas if m(i,j) will be in the
c    strict lower triangle of m, then m(j,i) is stored in row j (and
c    thus m(i,j) is not stored).
c
c
c  storage of sparse matrices
c
c    the nonzero entries of the matrix m are stored row-by-row in the
c    array a.  to identify the individual nonzero entries in each row,
c    we need to know in which column each entry lies.  these column
c    indices are stored in the array ja.  i.e., if  a(k) = m(i,j),  then
c    ja(k) = j.  to identify the individual rows, we need to know where
c    each row starts.  these row pointers are stored in the array ia.
c    i.e., if m(i,j) is the first nonzero entry (stored) in the i-th row
c    and  a(k) = m(i,j),  then  ia(i) = k.  moreover, ia(n+1) points to
c    the first location following the last element in the last row.
c    thus, the number of entries in the i-th row is  ia(i+1) - ia(i),
c    the nonzero entries in the i-th row are stored consecutively in
c
c            a(ia(i)),  a(ia(i)+1),  ..., a(ia(i+1)-1),
c
c    and the corresponding column indices are stored consecutively in
c
c            ja(ia(i)), ja(ia(i)+1), ..., ja(ia(i+1)-1).
c
c    when the coefficient matrix is symmetric, only the nonzero entries
c    in the upper triangle need be stored.  for example, the matrix
c
c             ( 1  0  2  3  0 )
c             ( 0  4  0  0  0 )
c         m = ( 2  0  5  6  0 )
c             ( 3  0  6  7  8 )
c             ( 0  0  0  8  9 )
c
c    could be stored as
c
c            - 1  2  3  4  5  6  7  8  9 10 11 12 13
c         ---+--------------------------------------
c         ia - 1  4  5  8 12 14
c         ja - 1  3  4  2  1  3  4  1  3  4  5  4  5
c          a - 1  2  3  4  2  5  6  3  6  7  8  8  9
c
c    or (symmetrically) as
c
c            - 1  2  3  4  5  6  7  8  9
c         ---+--------------------------
c         ia - 1  4  5  7  9 10
c         ja - 1  3  4  2  3  4  4  5  5
c          a - 1  2  3  4  5  6  7  8  9          .
c
c
c  parameters
c
c    n    - order of the matrix
c
c    ia   - integer one-dimensional array containing pointers to delimit
c           rows in ja and a.  dimension = n+1
c
c    ja   - integer one-dimensional array containing the column indices
c           corresponding to the elements of a.  dimension = number of
c           nonzero entries in (the upper triangle of) m
c
c    a    - real one-dimensional array containing the nonzero entries in
c           (the upper triangle of) m, stored by rows.  dimension =
c           number of nonzero entries in (the upper triangle of) m
c
c    p    - integer one-dimensional array used to return the permutation
c           of the rows and columns of m corresponding to the minimum
c           degree ordering.  dimension = n
c
c    ip   - integer one-dimensional array used to return the inverse of
c           the permutation returned in p.  dimension = n
c
c    nsp  - declared dimension of the one-dimensional array isp.  nsp
c           must be at least  3n+4k,  where k is the number of nonzeroes
c           in the strict upper triangle of m
c
c    isp  - integer one-dimensional array used for working storage.
c           dimension = nsp
c
c    path - integer path specification.  values and their meanings are -
c             1  find minimum degree ordering only
c             2  find minimum degree ordering and reorder symmetrically
c                  stored matrix (used when only the nonzero entries in
c                  the upper triangle of m are being stored)
c             3  reorder symmetrically stored matrix as specified by
c                  input permutation (used when an ordering has already
c                  been determined and only the nonzero entries in the
c                  upper triangle of m are being stored)
c             4  same as 2 but put diagonal entries at start of each row
c             5  same as 3 but put diagonal entries at start of each row
c
c    flag - integer error flag.  values and their meanings are -
c               0    no errors detected
c              9n+k  insufficient storage in md
c             10n+1  insufficient storage in odrv
c             11n+1  illegal path specification
c
c
c  conversion from real to double precision
c
c    change the real declarations in odrv and sro to double precision
c    declarations.
c
c-----------------------------------------------------------------------
c
      integer  ia(*), ja(*),  p(*), ip(*),  isp(*),  path,  flag,
     *   v, l, head,  tmp, q
c...  real  a(*)
      double precision  a(*)
      logical  dflag
c
c----initialize error flag and validate path specification
      flag = 0
      if (path.lt.1 .or. 5.lt.path)  go to 111
c
c----allocate storage and find minimum degree ordering
      if ((path-1) * (path-2) * (path-4) .ne. 0)  go to 1
        max = (nsp-n)/2
        v    = 1
        l    = v     +  max
        head = l     +  max
        next = head  +  n
        if (max.lt.n)  go to 110
c
        call  md
     *     (n, ia,ja, max,isp(v),isp(l), isp(head),p,ip, isp(v), flag)
        if (flag.ne.0)  go to 100
c
c----allocate storage and symmetrically reorder matrix
   1  if ((path-2) * (path-3) * (path-4) * (path-5) .ne. 0)  go to 2
        tmp = (nsp+1) -      n
        q   = tmp     - (ia(n+1)-1)
        if (q.lt.1)  go to 110
c
        dflag = path.eq.4 .or. path.eq.5
        call sro
     *     (n,  ip,  ia, ja, a,  isp(tmp),  isp(q),  dflag)
c
   2  return
c
c ** error -- error detected in md
 100  return
c ** error -- insufficient storage
 110  flag = 10*n + 1
      return
c ** error -- illegal path specified
 111  flag = 11*n + 1
      return
      end
      subroutine md
     *     (n, ia,ja, max, v,l, head,last,next, mark, flag)
c***********************************************************************
c  md -- minimum degree algorithm (based on element model)
c***********************************************************************
c
c  description
c
c    md finds a minimum degree ordering of the rows and columns of a
c    general sparse matrix m stored in (ia,ja,a) format.
c    when the structure of m is nonsymmetric, the ordering is that
c    obtained for the symmetric matrix  m + m-transpose.
c
c
c  additional parameters
c
c    max  - declared dimension of the one-dimensional arrays v and l.
c           max must be at least  n+2k,  where k is the number of
c           nonzeroes in the strict upper triangle of m + m-transpose
c
c    v    - integer one-dimensional work array.  dimension = max
c
c    l    - integer one-dimensional work array.  dimension = max
c
c    head - integer one-dimensional work array.  dimension = n
c
c    last - integer one-dimensional array used to return the permutation
c           of the rows and columns of m corresponding to the minimum
c           degree ordering.  dimension = n
c
c    next - integer one-dimensional array used to return the inverse of
c           the permutation returned in last.  dimension = n
c
c    mark - integer one-dimensional work array (may be the same as v).
c           dimension = n
c
c    flag - integer error flag.  values and their meanings are -
c             0     no errors detected
c             9n+k  insufficient storage in md
c
c
c  definitions of internal parameters
c
c    ---------+---------------------------------------------------------
c    v(s)     - value field of list entry
c    ---------+---------------------------------------------------------
c    l(s)     - link field of list entry  (0 =) end of list)
c    ---------+---------------------------------------------------------
c    l(vi)    - pointer to element list of uneliminated vertex vi
c    ---------+---------------------------------------------------------
c    l(ej)    - pointer to boundary list of active element ej
c    ---------+---------------------------------------------------------
c    head(d)  - vj =) vj head of d-list d
c             -  0 =) no vertex in d-list d
c
c
c             -                  vi uneliminated vertex
c             -          vi in ek           -       vi not in ek
c    ---------+-----------------------------+---------------------------
c    next(vi) - undefined but nonnegative   - vj =) vj next in d-list
c             -                             -  0 =) vi tail of d-list
c    ---------+-----------------------------+---------------------------
c    last(vi) - (not set until mdp)         - -d =) vi head of d-list d
c             --vk =) compute degree        - vj =) vj last in d-list
c             - ej =) vi prototype of ej    -  0 =) vi not in any d-list
c             -  0 =) do not compute degree -
c    ---------+-----------------------------+---------------------------
c    mark(vi) - mark(vk)                    - nonneg. tag .lt. mark(vk)
c
c
c             -                   vi eliminated vertex
c             -      ei active element      -           otherwise
c    ---------+-----------------------------+---------------------------
c    next(vi) - -j =) vi was j-th vertex    - -j =) vi was j-th vertex
c             -       to be eliminated      -       to be eliminated
c    ---------+-----------------------------+---------------------------
c    last(vi) -  m =) size of ei = m        - undefined
c    ---------+-----------------------------+---------------------------
c    mark(vi) - -m =) overlap count of ei   - undefined
c             -       with ek = m           -
c             - otherwise nonnegative tag   -
c             -       .lt. mark(vk)         -
c
c-----------------------------------------------------------------------
c
      integer  ia(*), ja(*),  v(*), l(*),  head(*), last(*), next(*),
     *   mark(*),  flag,  tag, dmin, vk,ek, tail
      equivalence  (vk,ek)
c
c----initialization
      tag = 0
      call  mdi
     *   (n, ia,ja, max,v,l, head,last,next, mark,tag, flag)
      if (flag.ne.0)  return
c
      k = 0
      dmin = 1
c
c----while  k .lt. n  do
   1  if (k.ge.n)  go to 4
c
c------search for vertex of minimum degree
   2    if (head(dmin).gt.0)  go to 3
          dmin = dmin + 1
          go to 2
c
c------remove vertex vk of minimum degree from degree list
   3    vk = head(dmin)
        head(dmin) = next(vk)
        if (head(dmin).gt.0)  last(head(dmin)) = -dmin
c
c------number vertex vk, adjust tag, and tag vk
        k = k+1
        next(vk) = -k
        last(ek) = dmin - 1
        tag = tag + last(ek)
        mark(vk) = tag
c
c------form element ek from uneliminated neighbors of vk
        call  mdm
     *     (vk,tail, v,l, last,next, mark)
c
c------purge inactive elements and do mass elimination
        call  mdp
     *     (k,ek,tail, v,l, head,last,next, mark)
c
c------update degrees of uneliminated vertices in ek
        call  mdu
     *     (ek,dmin, v,l, head,last,next, mark)
c
        go to 1
c
c----generate inverse permutation from permutation
   4  do 5 k=1,n
        next(k) = -next(k)
   5    last(next(k)) = k
c
      return
      end
      subroutine mdi
     *     (n, ia,ja, max,v,l, head,last,next, mark,tag, flag)
c***********************************************************************
c  mdi -- initialization
c***********************************************************************
      integer  ia(*), ja(*),  v(*), l(*),  head(*), last(*), next(*),
     *   mark(*), tag,  flag,  sfs, vi,dvi, vj
c
c----initialize degrees, element lists, and degree lists
      do 1 vi=1,n
        mark(vi) = 1
        l(vi) = 0
   1    head(vi) = 0
      sfs = n+1
c
c----create nonzero structure
c----for each nonzero entry a(vi,vj)
      do 6 vi=1,n
        jmin = ia(vi)
        jmax = ia(vi+1) - 1
        if (jmin.gt.jmax)  go to 6
        do 5 j=jmin,jmax
          vj = ja(j)
          if (vj-vi) 2, 5, 4
c
c------if a(vi,vj) is in strict lower triangle
c------check for previous occurrence of a(vj,vi)
   2      lvk = vi
          kmax = mark(vi) - 1
          if (kmax .eq. 0) go to 4
          do 3 k=1,kmax
            lvk = l(lvk)
            if (v(lvk).eq.vj) go to 5
   3        continue
c----for unentered entries a(vi,vj)
   4        if (sfs.ge.max)  go to 101
c
c------enter vj in element list for vi
            mark(vi) = mark(vi) + 1
            v(sfs) = vj
            l(sfs) = l(vi)
            l(vi) = sfs
            sfs = sfs+1
c
c------enter vi in element list for vj
            mark(vj) = mark(vj) + 1
            v(sfs) = vi
            l(sfs) = l(vj)
            l(vj) = sfs
            sfs = sfs+1
   5      continue
   6    continue
c
c----create degree lists and initialize mark vector
      do 7 vi=1,n
        dvi = mark(vi)
        next(vi) = head(dvi)
        head(dvi) = vi
        last(vi) = -dvi
        nextvi = next(vi)
        if (nextvi.gt.0)  last(nextvi) = vi
   7    mark(vi) = tag
c
      return
c
c ** error-  insufficient storage
 101  flag = 9*n + vi
      return
      end
      subroutine mdm
     *     (vk,tail, v,l, last,next, mark)
c***********************************************************************
c  mdm -- form element from uneliminated neighbors of vk
c***********************************************************************
      integer  vk, tail,  v(*), l(*),   last(*), next(*),   mark(*),
     *   tag, s,ls,vs,es, b,lb,vb, blp,blpmax
      equivalence  (vs, es)
c
c----initialize tag and list of uneliminated neighbors
      tag = mark(vk)
      tail = vk
c
c----for each vertex/element vs/es in element list of vk
      ls = l(vk)
   1  s = ls
      if (s.eq.0)  go to 5
        ls = l(s)
        vs = v(s)
        if (next(vs).lt.0)  go to 2
c
c------if vs is uneliminated vertex, then tag and append to list of
c------uneliminated neighbors
          mark(vs) = tag
          l(tail) = s
          tail = s
          go to 4
c
c------if es is active element, then ...
c--------for each vertex vb in boundary list of element es
   2      lb = l(es)
          blpmax = last(es)
          do 3 blp=1,blpmax
            b = lb
            lb = l(b)
            vb = v(b)
c
c----------if vb is untagged vertex, then tag and append to list of
c----------uneliminated neighbors
            if (mark(vb).ge.tag)  go to 3
              mark(vb) = tag
              l(tail) = b
              tail = b
   3        continue
c
c--------mark es inactive
          mark(es) = tag
c
   4    go to 1
c
c----terminate list of uneliminated neighbors
   5  l(tail) = 0
c
      return
      end
      subroutine mdp
     *     (k,ek,tail, v,l, head,last,next, mark)
c***********************************************************************
c  mdp -- purge inactive elements and do mass elimination
c***********************************************************************
      integer  ek, tail,  v(*), l(*),  head(*), last(*), next(*),
     *   mark(*),  tag, free, li,vi,lvi,evi, s,ls,es, ilp,ilpmax
c
c----initialize tag
      tag = mark(ek)
c
c----for each vertex vi in ek
      li = ek
      ilpmax = last(ek)
      if (ilpmax.le.0)  go to 12
      do 11 ilp=1,ilpmax
        i = li
        li = l(i)
        vi = v(li)
c
c------remove vi from degree list
        if (last(vi).eq.0)  go to 3
          if (last(vi).gt.0)  go to 1
            head(-last(vi)) = next(vi)
            go to 2
   1        next(last(vi)) = next(vi)
   2      if (next(vi).gt.0)  last(next(vi)) = last(vi)
c
c------remove inactive items from element list of vi
   3    ls = vi
   4    s = ls
        ls = l(s)
        if (ls.eq.0)  go to 6
          es = v(ls)
          if (mark(es).lt.tag)  go to 5
            free = ls
            l(s) = l(ls)
            ls = s
   5      go to 4
c
c------if vi is interior vertex, then remove from list and eliminate
   6    lvi = l(vi)
        if (lvi.ne.0)  go to 7
          l(i) = l(li)
          li = i
c
          k = k+1
          next(vi) = -k
          last(ek) = last(ek) - 1
          go to 11
c
c------else ...
c--------classify vertex vi
   7      if (l(lvi).ne.0)  go to 9
            evi = v(lvi)
            if (next(evi).ge.0)  go to 9
              if (mark(evi).lt.0)  go to 8
c
c----------if vi is prototype vertex, then mark as such, initialize
c----------overlap count for corresponding element, and move vi to end
c----------of boundary list
                last(vi) = evi
                mark(evi) = -1
                l(tail) = li
                tail = li
                l(i) = l(li)
                li = i
                go to 10
c
c----------else if vi is duplicate vertex, then mark as such and adjust
c----------overlap count for corresponding element
   8            last(vi) = 0
                mark(evi) = mark(evi) - 1
                go to 10
c
c----------else mark vi to compute degree
   9            last(vi) = -ek
c
c--------insert ek in element list of vi
  10      v(free) = ek
          l(free) = l(vi)
          l(vi) = free
  11    continue
c
c----terminate boundary list
  12  l(tail) = 0
c
      return
      end
      subroutine mdu
     *     (ek,dmin, v,l, head,last,next, mark)
c***********************************************************************
c  mdu -- update degrees of uneliminated vertices in ek
c***********************************************************************
      integer  ek, dmin,  v(*), l(*),  head(*), last(*), next(*),
     *   mark(*),  tag, vi,evi,dvi, s,vs,es, b,vb, ilp,ilpmax,
     *   blp,blpmax
      equivalence  (vs, es)
c
c----initialize tag
      tag = mark(ek) - last(ek)
c
c----for each vertex vi in ek
      i = ek
      ilpmax = last(ek)
      if (ilpmax.le.0)  go to 11
      do 10 ilp=1,ilpmax
        i = l(i)
        vi = v(i)
        if (last(vi))  1, 10, 8
c
c------if vi neither prototype nor duplicate vertex, then merge elements
c------to compute degree
   1      tag = tag + 1
          dvi = last(ek)
c
c--------for each vertex/element vs/es in element list of vi
          s = l(vi)
   2      s = l(s)
          if (s.eq.0)  go to 9
            vs = v(s)
            if (next(vs).lt.0)  go to 3
c
c----------if vs is uneliminated vertex, then tag and adjust degree
              mark(vs) = tag
              dvi = dvi + 1
              go to 5
c
c----------if es is active element, then expand
c------------check for outmatched vertex
   3          if (mark(es).lt.0)  go to 6
c
c------------for each vertex vb in es
              b = es
              blpmax = last(es)
              do 4 blp=1,blpmax
                b = l(b)
                vb = v(b)
c
c--------------if vb is untagged, then tag and adjust degree
                if (mark(vb).ge.tag)  go to 4
                  mark(vb) = tag
                  dvi = dvi + 1
   4            continue
c
   5        go to 2
c
c------else if vi is outmatched vertex, then adjust overlaps but do not
c------compute degree
   6      last(vi) = 0
          mark(es) = mark(es) - 1
   7      s = l(s)
          if (s.eq.0)  go to 10
            es = v(s)
            if (mark(es).lt.0)  mark(es) = mark(es) - 1
            go to 7
c
c------else if vi is prototype vertex, then calculate degree by
c------inclusion/exclusion and reset overlap count
   8      evi = last(vi)
          dvi = last(ek) + last(evi) + mark(evi)
          mark(evi) = 0
c
c------insert vi in appropriate degree list
   9    next(vi) = head(dvi)
        head(dvi) = vi
        last(vi) = -dvi
        if (next(vi).gt.0)  last(next(vi)) = vi
        if (dvi.lt.dmin)  dmin = dvi
c
  10    continue
c
  11  return
      end
      subroutine sro
     *     (n, ip, ia,ja,a, q, r, dflag)
c***********************************************************************
c  sro -- symmetric reordering of sparse symmetric matrix
c***********************************************************************
c
c  description
c
c    the nonzero entries of the matrix m are assumed to be stored
c    symmetrically in (ia,ja,a) format (i.e., not both m(i,j) and m(j,i)
c    are stored if i ne j).
c
c    sro does not rearrange the order of the rows, but does move
c    nonzeroes from one row to another to ensure that if m(i,j) will be
c    in the upper triangle of m with respect to the new ordering, then
c    m(i,j) is stored in row i (and thus m(j,i) is not stored),  whereas
c    if m(i,j) will be in the strict lower triangle of m, then m(j,i) is
c    stored in row j (and thus m(i,j) is not stored).
c
c
c  additional parameters
c
c    q     - integer one-dimensional work array.  dimension = n
c
c    r     - integer one-dimensional work array.  dimension = number of
c            nonzero entries in the upper triangle of m
c
c    dflag - logical variable.  if dflag = .true., then store nonzero
c            diagonal elements at the beginning of the row
c
c-----------------------------------------------------------------------
c
      integer  ip(*),  ia(*), ja(*),  q(*), r(*)
c...  real  a(*),  ak
      double precision  a(*),  ak
      logical  dflag
c
c
c--phase 1 -- find row in which to store each nonzero
c----initialize count of nonzeroes to be stored in each row
      do 1 i=1,n
  1     q(i) = 0
c
c----for each nonzero element a(j)
      do 3 i=1,n
        jmin = ia(i)
        jmax = ia(i+1) - 1
        if (jmin.gt.jmax)  go to 3
        do 2 j=jmin,jmax
c
c--------find row (=r(j)) and column (=ja(j)) in which to store a(j) ...
          k = ja(j)
          if (ip(k).lt.ip(i))  ja(j) = i
          if (ip(k).ge.ip(i))  k = i
          r(j) = k
c
c--------... and increment count of nonzeroes (=q(r(j)) in that row
  2       q(k) = q(k) + 1
  3     continue
c
c
c--phase 2 -- find new ia and permutation to apply to (ja,a)
c----determine pointers to delimit rows in permuted (ja,a)
      do 4 i=1,n
        ia(i+1) = ia(i) + q(i)
  4     q(i) = ia(i+1)
c
c----determine where each (ja(j),a(j)) is stored in permuted (ja,a)
c----for each nonzero element (in reverse order)
      ilast = 0
      jmin = ia(1)
      jmax = ia(n+1) - 1
      j = jmax
      do 6 jdummy=jmin,jmax
        i = r(j)
        if (.not.dflag .or. ja(j).ne.i .or. i.eq.ilast)  go to 5
c
c------if dflag, then put diagonal nonzero at beginning of row
          r(j) = ia(i)
          ilast = i
          go to 6
c
c------put (off-diagonal) nonzero in last unused location in row
  5       q(i) = q(i) - 1
          r(j) = q(i)
c
  6     j = j-1
c
c
c--phase 3 -- permute (ja,a) to upper triangular form (wrt new ordering)
      do 8 j=jmin,jmax
  7     if (r(j).eq.j)  go to 8
          k = r(j)
          r(j) = r(k)
          r(k) = k
          jak = ja(k)
          ja(k) = ja(j)
          ja(j) = jak
          ak = a(k)
          a(k) = a(j)
          a(j) = ak
          go to 7
  8     continue
c
      return
      end
*DECK CDRV
      subroutine cdrv
     *     (n, r,c,ic, ia,ja,a, b, z, nsp,isp,rsp,esp, path, flag)
c*** subroutine cdrv
c*** driver for subroutines for solving sparse nonsymmetric systems of
c       linear equations (compressed pointer storage)
c
c
c    parameters
c    class abbreviations are--
c       n - integer variable
c       f - real variable
c       v - supplies a value to the driver
c       r - returns a result from the driver
c       i - used internally by the driver
c       a - array
c
c class - parameter
c ------+----------
c       -
c         the nonzero entries of the coefficient matrix m are stored
c    row-by-row in the array a.  to identify the individual nonzero
c    entries in each row, we need to know in which column each entry
c    lies.  the column indices which correspond to the nonzero entries
c    of m are stored in the array ja.  i.e., if  a(k) = m(i,j),  then
c    ja(k) = j.  in addition, we need to know where each row starts and
c    how long it is.  the index positions in ja and a where the rows of
c    m begin are stored in the array ia.  i.e., if m(i,j) is the first
c    nonzero entry (stored) in the i-th row and a(k) = m(i,j),  then
c    ia(i) = k.  moreover, the index in ja and a of the first location
c    following the last element in the last row is stored in ia(n+1).
c    thus, the number of entries in the i-th row is given by
c    ia(i+1) - ia(i),  the nonzero entries of the i-th row are stored
c    consecutively in
c            a(ia(i)),  a(ia(i)+1),  ..., a(ia(i+1)-1),
c    and the corresponding column indices are stored consecutively in
c            ja(ia(i)), ja(ia(i)+1), ..., ja(ia(i+1)-1).
c    for example, the 5 by 5 matrix
c                ( 1. 0. 2. 0. 0.)
c                ( 0. 3. 0. 0. 0.)
c            m = ( 0. 4. 5. 6. 0.)
c                ( 0. 0. 0. 7. 0.)
c                ( 0. 0. 0. 8. 9.)
c    would be stored as
c               - 1  2  3  4  5  6  7  8  9
c            ---+--------------------------
c            ia - 1  3  4  7  8 10
c            ja - 1  3  2  2  3  4  4  4  5
c             a - 1. 2. 3. 4. 5. 6. 7. 8. 9.         .
c
c nv    - n     - number of variables/equations.
c fva   - a     - nonzero entries of the coefficient matrix m, stored
c       -           by rows.
c       -           size = number of nonzero entries in m.
c nva   - ia    - pointers to delimit the rows in a.
c       -           size = n+1.
c nva   - ja    - column numbers corresponding to the elements of a.
c       -           size = size of a.
c fva   - b     - right-hand side b.  b and z can the same array.
c       -           size = n.
c fra   - z     - solution x.  b and z can be the same array.
c       -           size = n.
c
c         the rows and columns of the original matrix m can be
c    reordered (e.g., to reduce fillin or ensure numerical stability)
c    before calling the driver.  if no reordering is done, then set
c    r(i) = c(i) = ic(i) = i  for i=1,...,n.  the solution z is returned
c    in the original order.
c         if the columns have been reordered (i.e.,  c(i).ne.i  for some
c    i), then the driver will call a subroutine (nroc) which rearranges
c    each row of ja and a, leaving the rows in the original order, but
c    placing the elements of each row in increasing order with respect
c    to the new ordering.  if  path.ne.1,  then nroc is assumed to have
c    been called already.
c
c nva   - r     - ordering of the rows of m.
c       -           size = n.
c nva   - c     - ordering of the columns of m.
c       -           size = n.
c nva   - ic    - inverse of the ordering of the columns of m.  i.e.,
c       -           ic(c(i)) = i  for i=1,...,n.
c       -           size = n.
c
c         the solution of the system of linear equations is divided into
c    three stages --
c      nsfc -- the matrix m is processed symbolically to determine where
c               fillin will occur during the numeric factorization.
c      nnfc -- the matrix m is factored numerically into the product ldu
c               of a unit lower triangular matrix l, a diagonal matrix
c               d, and a unit upper triangular matrix u, and the system
c               mx = b  is solved.
c      nnsc -- the linear system  mx = b  is solved using the ldu
c  or           factorization from nnfc.
c      nntc -- the transposed linear system  mt x = b  is solved using
c               the ldu factorization from nnf.
c    for several systems whose coefficient matrices have the same
c    nonzero structure, nsfc need be done only once (for the first
c    system).  then nnfc is done once for each additional system.  for
c    several systems with the same coefficient matrix, nsfc and nnfc
c    need be done only once (for the first system).  then nnsc or nntc
c    is done once for each additional right-hand side.
c
c nv    - path  - path specification.  values and their meanings are --
c       -           1  perform nroc, nsfc, and nnfc.
c       -           2  perform nnfc only  (nsfc is assumed to have been
c       -               done in a manner compatible with the storage
c       -               allocation used in the driver).
c       -           3  perform nnsc only  (nsfc and nnfc are assumed to
c       -               have been done in a manner compatible with the
c       -               storage allocation used in the driver).
c       -           4  perform nntc only  (nsfc and nnfc are assumed to
c       -               have been done in a manner compatible with the
c       -               storage allocation used in the driver).
c       -           5  perform nroc and nsfc.
c
c         various errors are detected by the driver and the individual
c    subroutines.
c
c nr    - flag  - error flag.  values and their meanings are --
c       -             0     no errors detected
c       -             n+k   null row in a  --  row = k
c       -            2n+k   duplicate entry in a  --  row = k
c       -            3n+k   insufficient storage in nsfc  --  row = k
c       -            4n+1   insufficient storage in nnfc
c       -            5n+k   null pivot  --  row = k
c       -            6n+k   insufficient storage in nsfc  --  row = k
c       -            7n+1   insufficient storage in nnfc
c       -            8n+k   zero pivot  --  row = k
c       -           10n+1   insufficient storage in cdrv
c       -           11n+1   illegal path specification
c
c         working storage is needed for the factored form of the matrix
c    m plus various temporary vectors.  the arrays isp and rsp should be
c    equivalenced.  integer storage is allocated from the beginning of
c    isp and real storage from the end of rsp.
c
c nv    - nsp   - declared dimension of rsp.  nsp generally must
c       -           be larger than  8n+2 + 2k  (where  k = (number of
c       -           nonzero entries in m)).
c nvira - isp   - integer working storage divided up into various arrays
c       -           needed by the subroutines.  isp and rsp should be
c       -           equivalenced.
c       -           size = lratio*nsp.
c fvira - rsp   - real working storage divided up into various arrays
c       -           needed by the subroutines.  isp and rsp should be
c       -           equivalenced.
c       -           size = nsp.
c nr    - esp   - if sufficient storage was available to perform the
c       -           symbolic factorization (nsfc), then esp is set to
c       -           the amount of excess storage provided (negative if
c       -           insufficient storage was available to perform the
c       -           numeric factorization (nnfc)).
c
c
c  conversion to double precision
c
c    to convert these routines for double precision arrays..
c    (1) use the double precision declarations in place of the real
c    declarations in each subprogram, as given in comment cards.
c    (2) change the data-loaded value of the integer  lratio
c    in subroutine cdrv, as indicated below.
c    (3) change e0 to d0 in the constants in statement number 10
c    in subroutine nnfc and the line following that.
c
      integer  r(*), c(*), ic(*),  ia(*), ja(*),  isp(*), esp,  path,
     *   flag,  d, u, q, row, tmp, ar,  umax
c     real  a(*), b(*), z(*), rsp(*)
      double precision  a(*), b(*), z(*), rsp(*)
c
c  set lratio equal to the ratio between the length of floating point
c  and integer array data.  e. g., lratio = 1 for (real, integer),
c  lratio = 2 for (double precision, integer)
c
      data lratio/2/
c
      if (path.lt.1 .or. 5.lt.path)  go to 111
c******initialize and divide up temporary storage  *******************
      il   = 1
      ijl  = il  + (n+1)
      iu   = ijl +   n
      iju  = iu  + (n+1)
      irl  = iju +   n
      jrl  = irl +   n
      jl   = jrl +   n
c
c  ******  reorder a if necessary, call nsfc if flag is set  ***********
      if ((path-1) * (path-5) .ne. 0)  go to 5
        max = (lratio*nsp + 1 - jl) - (n+1) - 5*n
        jlmax = max/2
        q     = jl   + jlmax
        ira   = q    + (n+1)
        jra   = ira  +   n
        irac  = jra  +   n
        iru   = irac +   n
        jru   = iru  +   n
        jutmp = jru  +   n
        jumax = lratio*nsp  + 1 - jutmp
        esp = max/lratio
        if (jlmax.le.0 .or. jumax.le.0)  go to 110
c
        do 1 i=1,n
          if (c(i).ne.i)  go to 2
   1      continue
        go to 3
   2    ar = nsp + 1 - n
        call  nroc
     *     (n, ic, ia,ja,a, isp(il), rsp(ar), isp(iu), flag)
        if (flag.ne.0)  go to 100
c
   3    call  nsfc
     *     (n, r, ic, ia,ja,
     *      jlmax, isp(il), isp(jl), isp(ijl),
     *      jumax, isp(iu), isp(jutmp), isp(iju),
     *      isp(q), isp(ira), isp(jra), isp(irac),
     *      isp(irl), isp(jrl), isp(iru), isp(jru),  flag)
        if(flag .ne. 0)  go to 100
c  ******  move ju next to jl  *****************************************
        jlmax = isp(ijl+n-1)
        ju    = jl + jlmax
        jumax = isp(iju+n-1)
        if (jumax.le.0)  go to 5
        do 4 j=1,jumax
   4      isp(ju+j-1) = isp(jutmp+j-1)
c
c  ******  call remaining subroutines  *********************************
   5  jlmax = isp(ijl+n-1)
      ju    = jl  + jlmax
      jumax = isp(iju+n-1)
      l     = (ju + jumax - 2 + lratio)  /  lratio    +    1
      lmax  = isp(il+n) - 1
      d     = l   + lmax
      u     = d   + n
      row   = nsp + 1 - n
      tmp   = row - n
      umax  = tmp - u
      esp   = umax - (isp(iu+n) - 1)
c
      if ((path-1) * (path-2) .ne. 0)  go to 6
        if (umax.lt.0)  go to 110
        call nnfc
     *     (n,  r, c, ic,  ia, ja, a, z, b,
     *      lmax, isp(il), isp(jl), isp(ijl), rsp(l),  rsp(d),
     *      umax, isp(iu), isp(ju), isp(iju), rsp(u),
     *      rsp(row), rsp(tmp),  isp(irl), isp(jrl),  flag)
        if(flag .ne. 0)  go to 100
c
   6  if ((path-3) .ne. 0)  go to 7
        call nnsc
     *     (n,  r, c,  isp(il), isp(jl), isp(ijl), rsp(l),
     *      rsp(d),    isp(iu), isp(ju), isp(iju), rsp(u),
     *      z, b,  rsp(tmp))
c
   7  if ((path-4) .ne. 0)  go to 8
        call nntc
     *     (n,  r, c,  isp(il), isp(jl), isp(ijl), rsp(l),
     *      rsp(d),    isp(iu), isp(ju), isp(iju), rsp(u),
     *      z, b,  rsp(tmp))
   8  return
c
c ** error.. error detected in nroc, nsfc, nnfc, or nnsc
 100  return
c ** error.. insufficient storage
 110  flag = 10*n + 1
      return
c ** error.. illegal path specification
 111  flag = 11*n + 1
      return
      end
      subroutine nroc (n, ic, ia, ja, a, jar, ar, p, flag)
c
c       ----------------------------------------------------------------
c
c               yale sparse matrix package - nonsymmetric codes
c                    solving the system of equations mx = b
c
c    i.   calling sequences
c         the coefficient matrix can be processed by an ordering routine
c    (e.g., to reduce fillin or ensure numerical stability) before using
c    the remaining subroutines.  if no reordering is done, then set
c    r(i) = c(i) = ic(i) = i  for i=1,...,n.  if an ordering subroutine
c    is used, then nroc should be used to reorder the coefficient matrix
c    the calling sequence is --
c        (       (matrix ordering))
c        (nroc   (matrix reordering))
c         nsfc   (symbolic factorization to determine where fillin will
c                  occur during numeric factorization)
c         nnfc   (numeric factorization into product ldu of unit lower
c                  triangular matrix l, diagonal matrix d, and unit
c                  upper triangular matrix u, and solution of linear
c                  system)
c         nnsc   (solution of linear system for additional right-hand
c                  side using ldu factorization from nnfc)
c    (if only one system of equations is to be solved, then the
c    subroutine trk should be used.)
c
c    ii.  storage of sparse matrices
c         the nonzero entries of the coefficient matrix m are stored
c    row-by-row in the array a.  to identify the individual nonzero
c    entries in each row, we need to know in which column each entry
c    lies.  the column indices which correspond to the nonzero entries
c    of m are stored in the array ja.  i.e., if  a(k) = m(i,j),  then
c    ja(k) = j.  in addition, we need to know where each row starts and
c    how long it is.  the index positions in ja and a where the rows of
c    m begin are stored in the array ia.  i.e., if m(i,j) is the first
c    (leftmost) entry in the i-th row and  a(k) = m(i,j),  then
c    ia(i) = k.  moreover, the index in ja and a of the first location
c    following the last element in the last row is stored in ia(n+1).
c    thus, the number of entries in the i-th row is given by
c    ia(i+1) - ia(i),  the nonzero entries of the i-th row are stored
c    consecutively in
c            a(ia(i)),  a(ia(i)+1),  ..., a(ia(i+1)-1),
c    and the corresponding column indices are stored consecutively in
c            ja(ia(i)), ja(ia(i)+1), ..., ja(ia(i+1)-1).
c    for example, the 5 by 5 matrix
c                ( 1. 0. 2. 0. 0.)
c                ( 0. 3. 0. 0. 0.)
c            m = ( 0. 4. 5. 6. 0.)
c                ( 0. 0. 0. 7. 0.)
c                ( 0. 0. 0. 8. 9.)
c    would be stored as
c               - 1  2  3  4  5  6  7  8  9
c            ---+--------------------------
c            ia - 1  3  4  7  8 10
c            ja - 1  3  2  2  3  4  4  4  5
c             a - 1. 2. 3. 4. 5. 6. 7. 8. 9.         .
c
c         the strict upper (lower) triangular portion of the matrix
c    u (l) is stored in a similar fashion using the arrays  iu, ju, u
c    (il, jl, l)  except that an additional array iju (ijl) is used to
c    compress storage of ju (jl) by allowing some sequences of column
c    (row) indices to used for more than one row (column)  (n.b., l is
c    stored by columns).  iju(k) (ijl(k)) points to the starting
c    location in ju (jl) of entries for the kth row (column).
c    compression in ju (jl) occurs in two ways.  first, if a row
c    (column) i was merged into the current row (column) k, and the
c    number of elements merged in from (the tail portion of) row
c    (column) i is the same as the final length of row (column) k, then
c    the kth row (column) and the tail of row (column) i are identical
c    and iju(k) (ijl(k)) points to the start of the tail.  second, if
c    some tail portion of the (k-1)st row (column) is identical to the
c    head of the kth row (column), then iju(k) (ijl(k)) points to the
c    start of that tail portion.  for example, the nonzero structure of
c    the strict upper triangular part of the matrix
c            d 0 x x x
c            0 d 0 x x
c            0 0 d x 0
c            0 0 0 d x
c            0 0 0 0 d
c    would be represented as
c                - 1 2 3 4 5 6
c            ----+------------
c             iu - 1 4 6 7 8 8
c             ju - 3 4 5 4
c            iju - 1 2 4 3           .
c    the diagonal entries of l and u are assumed to be equal to one and
c    are not stored.  the array d contains the reciprocals of the
c    diagonal entries of the matrix d.
c
c    iii. additional storage savings
c         in nsfc, r and ic can be the same array in the calling
c    sequence if no reordering of the coefficient matrix has been done.
c         in nnfc, r, c, and ic can all be the same array if no
c    reordering has been done.  if only the rows have been reordered,
c    then c and ic can be the same array.  if the row and column
c    orderings are the same, then r and c can be the same array.  z and
c    row can be the same array.
c         in nnsc or nntc, r and c can be the same array if no
c    reordering has been done or if the row and column orderings are the
c    same.  z and b can be the same array.  however, then b will be
c    destroyed.
c
c    iv.  parameters
c         following is a list of parameters to the programs.  names are
c    uniform among the various subroutines.  class abbreviations are --
c       n - integer variable
c       f - real variable
c       v - supplies a value to a subroutine
c       r - returns a result from a subroutine
c       i - used internally by a subroutine
c       a - array
c
c class - parameter
c ------+----------
c fva   - a     - nonzero entries of the coefficient matrix m, stored
c       -           by rows.
c       -           size = number of nonzero entries in m.
c fva   - b     - right-hand side b.
c       -           size = n.
c nva   - c     - ordering of the columns of m.
c       -           size = n.
c fvra  - d     - reciprocals of the diagonal entries of the matrix d.
c       -           size = n.
c nr    - flag  - error flag.  values and their meanings are --
c       -            0     no errors detected
c       -            n+k   null row in a  --  row = k
c       -           2n+k   duplicate entry in a  --  row = k
c       -           3n+k   insufficient storage for jl  --  row = k
c       -           4n+1   insufficient storage for l
c       -           5n+k   null pivot  --  row = k
c       -           6n+k   insufficient storage for ju  --  row = k
c       -           7n+1   insufficient storage for u
c       -           8n+k   zero pivot  --  row = k
c nva   - ia    - pointers to delimit the rows of a.
c       -           size = n+1.
c nvra  - ijl   - pointers to the first element in each column in jl,
c       -           used to compress storage in jl.
c       -           size = n.
c nvra  - iju   - pointers to the first element in each row in ju, used
c       -           to compress storage in ju.
c       -           size = n.
c nvra  - il    - pointers to delimit the columns of l.
c       -           size = n+1.
c nvra  - iu    - pointers to delimit the rows of u.
c       -           size = n+1.
c nva   - ja    - column numbers corresponding to the elements of a.
c       -           size = size of a.
c nvra  - jl    - row numbers corresponding to the elements of l.
c       -           size = jlmax.
c nv    - jlmax - declared dimension of jl.  jlmax must be larger than
c       -           the number of nonzeros in the strict lower triangle
c       -           of m plus fillin minus compression.
c nvra  - ju    - column numbers corresponding to the elements of u.
c       -           size = jumax.
c nv    - jumax - declared dimension of ju.  jumax must be larger than
c       -           the number of nonzeros in the strict upper triangle
c       -           of m plus fillin minus compression.
c fvra  - l     - nonzero entries in the strict lower triangular portion
c       -           of the matrix l, stored by columns.
c       -           size = lmax.
c nv    - lmax  - declared dimension of l.  lmax must be larger than
c       -           the number of nonzeros in the strict lower triangle
c       -           of m plus fillin  (il(n+1)-1 after nsfc).
c nv    - n     - number of variables/equations.
c nva   - r     - ordering of the rows of m.
c       -           size = n.
c fvra  - u     - nonzero entries in the strict upper triangular portion
c       -           of the matrix u, stored by rows.
c       -           size = umax.
c nv    - umax  - declared dimension of u.  umax must be larger than
c       -           the number of nonzeros in the strict upper triangle
c       -           of m plus fillin  (iu(n+1)-1 after nsfc).
c fra   - z     - solution x.
c       -           size = n.
c
c       ----------------------------------------------------------------
c
c*** subroutine nroc
c*** reorders rows of a, leaving row order unchanged
c
c
c       input parameters.. n, ic, ia, ja, a
c       output parameters.. ja, a, flag
c
c       parameters used internally..
c nia   - p     - at the kth step, p is a linked list of the reordered
c       -           column indices of the kth row of a.  p(n+1) points
c       -           to the first entry in the list.
c       -           size = n+1.
c nia   - jar   - at the kth step,jar contains the elements of the
c       -           reordered column indices of a.
c       -           size = n.
c fia   - ar    - at the kth step, ar contains the elements of the
c       -           reordered row of a.
c       -           size = n.
c
      integer  ic(*), ia(*), ja(*), jar(*), p(*), flag
c     real  a(*), ar(*)
      double precision  a(*), ar(*)
c
c  ******  for each nonempty row  *******************************
      do 5 k=1,n
        jmin = ia(k)
        jmax = ia(k+1) - 1
        if(jmin .gt. jmax) go to 5
        p(n+1) = n + 1
c  ******  insert each element in the list  *********************
        do 3 j=jmin,jmax
          newj = ic(ja(j))
          i = n + 1
   1      if(p(i) .ge. newj) go to 2
            i = p(i)
            go to 1
   2      if(p(i) .eq. newj) go to 102
          p(newj) = p(i)
          p(i) = newj
          jar(newj) = ja(j)
          ar(newj) = a(j)
   3      continue
c  ******  replace old row in ja and a  *************************
        i = n + 1
        do 4 j=jmin,jmax
          i = p(i)
          ja(j) = jar(i)
   4      a(j) = ar(i)
   5    continue
      flag = 0
      return
c
c ** error.. duplicate entry in a
 102  flag = n + k
      return
      end
      subroutine nsfc
     *      (n, r, ic, ia,ja, jlmax,il,jl,ijl, jumax,iu,ju,iju,
     *       q, ira,jra, irac, irl,jrl, iru,jru, flag)
c*** subroutine nsfc
c*** symbolic ldu-factorization of nonsymmetric sparse matrix
c      (compressed pointer storage)
c
c
c       input variables.. n, r, ic, ia, ja, jlmax, jumax.
c       output variables.. il, jl, ijl, iu, ju, iju, flag.
c
c       parameters used internally..
c nia   - q     - suppose  m*  is the result of reordering  m.  if
c       -           processing of the ith row of  m*  (hence the ith
c       -           row of  u) is being done,  q(j)  is initially
c       -           nonzero if  m*(i,j) is nonzero (j.ge.i).  since
c       -           values need not be stored, each entry points to the
c       -           next nonzero and  q(n+1)  points to the first.  n+1
c       -           indicates the end of the list.  for example, if n=9
c       -           and the 5th row of  m*  is
c       -              0 x x 0 x 0 0 x 0
c       -           then  q  will initially be
c       -              a a a a 8 a a 10 5           (a - arbitrary).
c       -           as the algorithm proceeds, other elements of  q
c       -           are inserted in the list because of fillin.
c       -           q  is used in an analogous manner to compute the
c       -           ith column of  l.
c       -           size = n+1.
c nia   - ira,  - vectors used to find the columns of  m.  at the kth
c nia   - jra,      step of the factorization,  irac(k)  points to the
c nia   - irac      head of a linked list in  jra  of row indices i
c       -           such that i .ge. k and  m(i,k)  is nonzero.  zero
c       -           indicates the end of the list.  ira(i)  (i.ge.k)
c       -           points to the smallest j such that j .ge. k and
c       -           m(i,j)  is nonzero.
c       -           size of each = n.
c nia   - irl,  - vectors used to find the rows of  l.  at the kth step
c nia   - jrl       of the factorization,  jrl(k)  points to the head
c       -           of a linked list in  jrl  of column indices j
c       -           such j .lt. k and  l(k,j)  is nonzero.  zero
c       -           indicates the end of the list.  irl(j)  (j.lt.k)
c       -           points to the smallest i such that i .ge. k and
c       -           l(i,j)  is nonzero.
c       -           size of each = n.
c nia   - iru,  - vectors used in a manner analogous to  irl and jrl
c nia   - jru       to find the columns of  u.
c       -           size of each = n.
c
c  internal variables..
c    jlptr - points to the last position used in  jl.
c    juptr - points to the last position used in  ju.
c    jmin,jmax - are the indices in  a or u  of the first and last
c                elements to be examined in a given row.
c                for example,  jmin=ia(k), jmax=ia(k+1)-1.
c
      integer cend, qm, rend, rk, vj
      integer ia(*), ja(*), ira(*), jra(*), il(*), jl(*), ijl(*)
      integer iu(*), ju(*), iju(*), irl(*), jrl(*), iru(*), jru(*)
      integer r(*), ic(*), q(*), irac(*), flag
c
c  ******  initialize pointers  ****************************************
      np1 = n + 1
      jlmin = 1
      jlptr = 0
      il(1) = 1
      jumin = 1
      juptr = 0
      iu(1) = 1
      do 1 k=1,n
        irac(k) = 0
        jra(k) = 0
        jrl(k) = 0
   1    jru(k) = 0
c  ******  initialize column pointers for a  ***************************
      do 2 k=1,n
        rk = r(k)
        iak = ia(rk)
        if (iak .ge. ia(rk+1))  go to 101
        jaiak = ic(ja(iak))
        if (jaiak .gt. k)  go to 105
        jra(k) = irac(jaiak)
        irac(jaiak) = k
   2    ira(k) = iak
c
c  ******  for each column of l and row of u  **************************
      do 41 k=1,n
c
c  ******  initialize q for computing kth column of l  *****************
        q(np1) = np1
        luk = -1
c  ******  by filling in kth column of a  ******************************
        vj = irac(k)
        if (vj .eq. 0)  go to 5
   3      qm = np1
   4      m = qm
          qm =  q(m)
          if (qm .lt. vj)  go to 4
          if (qm .eq. vj)  go to 102
            luk = luk + 1
            q(m) = vj
            q(vj) = qm
            vj = jra(vj)
            if (vj .ne. 0)  go to 3
c  ******  link through jru  *******************************************
   5    lastid = 0
        lasti = 0
        ijl(k) = jlptr
        i = k
   6      i = jru(i)
          if (i .eq. 0)  go to 10
          qm = np1
          jmin = irl(i)
          jmax = ijl(i) + il(i+1) - il(i) - 1
          long = jmax - jmin
          if (long .lt. 0)  go to 6
          jtmp = jl(jmin)
          if (jtmp .ne. k)  long = long + 1
          if (jtmp .eq. k)  r(i) = -r(i)
          if (lastid .ge. long)  go to 7
            lasti = i
            lastid = long
c  ******  and merge the corresponding columns into the kth column  ****
   7      do 9 j=jmin,jmax
            vj = jl(j)
   8        m = qm
            qm = q(m)
            if (qm .lt. vj)  go to 8
            if (qm .eq. vj)  go to 9
              luk = luk + 1
              q(m) = vj
              q(vj) = qm
              qm = vj
   9        continue
            go to 6
c  ******  lasti is the longest column merged into the kth  ************
c  ******  see if it equals the entire kth column  *********************
  10    qm = q(np1)
        if (qm .ne. k)  go to 105
        if (luk .eq. 0)  go to 17
        if (lastid .ne. luk)  go to 11
c  ******  if so, jl can be compressed  ********************************
        irll = irl(lasti)
        ijl(k) = irll + 1
        if (jl(irll) .ne. k)  ijl(k) = ijl(k) - 1
        go to 17
c  ******  if not, see if kth column can overlap the previous one  *****
  11    if (jlmin .gt. jlptr)  go to 15
        qm = q(qm)
        do 12 j=jlmin,jlptr
          if (jl(j) - qm)  12, 13, 15
  12      continue
        go to 15
  13    ijl(k) = j
        do 14 i=j,jlptr
          if (jl(i) .ne. qm)  go to 15
          qm = q(qm)
          if (qm .gt. n)  go to 17
  14      continue
        jlptr = j - 1
c  ******  move column indices from q to jl, update vectors  ***********
  15    jlmin = jlptr + 1
        ijl(k) = jlmin
        if (luk .eq. 0)  go to 17
        jlptr = jlptr + luk
        if (jlptr .gt. jlmax)  go to 103
          qm = q(np1)
          do 16 j=jlmin,jlptr
            qm = q(qm)
  16        jl(j) = qm
  17    irl(k) = ijl(k)
        il(k+1) = il(k) + luk
c
c  ******  initialize q for computing kth row of u  ********************
        q(np1) = np1
        luk = -1
c  ******  by filling in kth row of reordered a  ***********************
        rk = r(k)
        jmin = ira(k)
        jmax = ia(rk+1) - 1
        if (jmin .gt. jmax)  go to 20
        do 19 j=jmin,jmax
          vj = ic(ja(j))
          qm = np1
  18      m = qm
          qm = q(m)
          if (qm .lt. vj)  go to 18
          if (qm .eq. vj)  go to 102
            luk = luk + 1
            q(m) = vj
            q(vj) = qm
  19      continue
c  ******  link through jrl,  ******************************************
  20    lastid = 0
        lasti = 0
        iju(k) = juptr
        i = k
        i1 = jrl(k)
  21      i = i1
          if (i .eq. 0)  go to 26
          i1 = jrl(i)
          qm = np1
          jmin = iru(i)
          jmax = iju(i) + iu(i+1) - iu(i) - 1
          long = jmax - jmin
          if (long .lt. 0)  go to 21
          jtmp = ju(jmin)
          if (jtmp .eq. k)  go to 22
c  ******  update irl and jrl, *****************************************
            long = long + 1
            cend = ijl(i) + il(i+1) - il(i)
            irl(i) = irl(i) + 1
            if (irl(i) .ge. cend)  go to 22
              j = jl(irl(i))
              jrl(i) = jrl(j)
              jrl(j) = i
  22      if (lastid .ge. long)  go to 23
            lasti = i
            lastid = long
c  ******  and merge the corresponding rows into the kth row  **********
  23      do 25 j=jmin,jmax
            vj = ju(j)
  24        m = qm
            qm = q(m)
            if (qm .lt. vj)  go to 24
            if (qm .eq. vj)  go to 25
              luk = luk + 1
              q(m) = vj
              q(vj) = qm
              qm = vj
  25        continue
          go to 21
c  ******  update jrl(k) and irl(k)  ***********************************
  26    if (il(k+1) .le. il(k))  go to 27
          j = jl(irl(k))
          jrl(k) = jrl(j)
          jrl(j) = k
c  ******  lasti is the longest row merged into the kth  ***************
c  ******  see if it equals the entire kth row  ************************
  27    qm = q(np1)
        if (qm .ne. k)  go to 105
        if (luk .eq. 0)  go to 34
        if (lastid .ne. luk)  go to 28
c  ******  if so, ju can be compressed  ********************************
        irul = iru(lasti)
        iju(k) = irul + 1
        if (ju(irul) .ne. k)  iju(k) = iju(k) - 1
        go to 34
c  ******  if not, see if kth row can overlap the previous one  ********
  28    if (jumin .gt. juptr)  go to 32
        qm = q(qm)
        do 29 j=jumin,juptr
          if (ju(j) - qm)  29, 30, 32
  29      continue
        go to 32
  30    iju(k) = j
        do 31 i=j,juptr
          if (ju(i) .ne. qm)  go to 32
          qm = q(qm)
          if (qm .gt. n)  go to 34
  31      continue
        juptr = j - 1
c  ******  move row indices from q to ju, update vectors  **************
  32    jumin = juptr + 1
        iju(k) = jumin
        if (luk .eq. 0)  go to 34
        juptr = juptr + luk
        if (juptr .gt. jumax)  go to 106
          qm = q(np1)
          do 33 j=jumin,juptr
            qm = q(qm)
  33        ju(j) = qm
  34    iru(k) = iju(k)
        iu(k+1) = iu(k) + luk
c
c  ******  update iru, jru  ********************************************
        i = k
  35      i1 = jru(i)
          if (r(i) .lt. 0)  go to 36
          rend = iju(i) + iu(i+1) - iu(i)
          if (iru(i) .ge. rend)  go to 37
            j = ju(iru(i))
            jru(i) = jru(j)
            jru(j) = i
            go to 37
  36      r(i) = -r(i)
  37      i = i1
          if (i .eq. 0)  go to 38
          iru(i) = iru(i) + 1
          go to 35
c
c  ******  update ira, jra, irac  **************************************
  38    i = irac(k)
        if (i .eq. 0)  go to 41
  39      i1 = jra(i)
          ira(i) = ira(i) + 1
          if (ira(i) .ge. ia(r(i)+1))  go to 40
          irai = ira(i)
          jairai = ic(ja(irai))
          if (jairai .gt. i)  go to 40
          jra(i) = irac(jairai)
          irac(jairai) = i
  40      i = i1
          if (i .ne. 0)  go to 39
  41    continue
c
      ijl(n) = jlptr
      iju(n) = juptr
      flag = 0
      return
c
c ** error.. null row in a
 101  flag = n + rk
      return
c ** error.. duplicate entry in a
 102  flag = 2*n + rk
      return
c ** error.. insufficient storage for jl
 103  flag = 3*n + k
      return
c ** error.. null pivot
 105  flag = 5*n + k
      return
c ** error.. insufficient storage for ju
 106  flag = 6*n + k
      return
      end
      subroutine nnfc
     *     (n, r,c,ic, ia,ja,a, z, b,
     *      lmax,il,jl,ijl,l, d, umax,iu,ju,iju,u,
     *      row, tmp, irl,jrl, flag)
c*** subroutine nnfc
c*** numerical ldu-factorization of sparse nonsymmetric matrix and
c      solution of system of linear equations (compressed pointer
c      storage)
c
c
c       input variables..  n, r, c, ic, ia, ja, a, b,
c                          il, jl, ijl, lmax, iu, ju, iju, umax
c       output variables.. z, l, d, u, flag
c
c       parameters used internally..
c nia   - irl,  - vectors used to find the rows of  l.  at the kth step
c nia   - jrl       of the factorization,  jrl(k)  points to the head
c       -           of a linked list in  jrl  of column indices j
c       -           such j .lt. k and  l(k,j)  is nonzero.  zero
c       -           indicates the end of the list.  irl(j)  (j.lt.k)
c       -           points to the smallest i such that i .ge. k and
c       -           l(i,j)  is nonzero.
c       -           size of each = n.
c fia   - row   - holds intermediate values in calculation of  u and l.
c       -           size = n.
c fia   - tmp   - holds new right-hand side  b*  for solution of the
c       -           equation ux = b*.
c       -           size = n.
c
c  internal variables..
c    jmin, jmax - indices of the first and last positions in a row to
c      be examined.
c    sum - used in calculating  tmp.
c
      integer rk,umax
      integer  r(*), c(*), ic(*), ia(*), ja(*), il(*), jl(*), ijl(*)
      integer  iu(*), ju(*), iju(*), irl(*), jrl(*), flag
c     real  a(*), l(*), d(*), u(*), z(*), b(*), row(*)
c     real tmp(*), lki, sum, dk
      double precision  a(*), l(*), d(*), u(*), z(*), b(*), row(*)
      double precision  tmp(*), lki, sum, dk
c
c  ******  initialize pointers and test storage  ***********************
      if(il(n+1)-1 .gt. lmax) go to 104
      if(iu(n+1)-1 .gt. umax) go to 107
      do 1 k=1,n
        irl(k) = il(k)
        jrl(k) = 0
   1    continue
c
c  ******  for each row  ***********************************************
      do 19 k=1,n
c  ******  reverse jrl and zero row where kth row of l will fill in  ***
        row(k) = 0
        i1 = 0
        if (jrl(k) .eq. 0) go to 3
        i = jrl(k)
   2    i2 = jrl(i)
        jrl(i) = i1
        i1 = i
        row(i) = 0
        i = i2
        if (i .ne. 0) go to 2
c  ******  set row to zero where u will fill in  ***********************
   3    jmin = iju(k)
        jmax = jmin + iu(k+1) - iu(k) - 1
        if (jmin .gt. jmax) go to 5
        do 4 j=jmin,jmax
   4      row(ju(j)) = 0
c  ******  place kth row of a in row  **********************************
   5    rk = r(k)
        jmin = ia(rk)
        jmax = ia(rk+1) - 1
        do 6 j=jmin,jmax
          row(ic(ja(j))) = a(j)
   6      continue
c  ******  initialize sum, and link through jrl  ***********************
        sum = b(rk)
        i = i1
        if (i .eq. 0) go to 10
c  ******  assign the kth row of l and adjust row, sum  ****************
   7      lki = -row(i)
c  ******  if l is not required, then comment out the following line  **
          l(irl(i)) = -lki
          sum = sum + lki * tmp(i)
          jmin = iu(i)
          jmax = iu(i+1) - 1
          if (jmin .gt. jmax) go to 9
          mu = iju(i) - jmin
          do 8 j=jmin,jmax
   8        row(ju(mu+j)) = row(ju(mu+j)) + lki * u(j)
   9      i = jrl(i)
          if (i .ne. 0) go to 7
c
c  ******  assign kth row of u and diagonal d, set tmp(k)  *************
  10    if (row(k) .eq. 0.0d0) go to 108
        dk = 1.0d0 / row(k)
        d(k) = dk
        tmp(k) = sum * dk
        if (k .eq. n) go to 19
        jmin = iu(k)
        jmax = iu(k+1) - 1
        if (jmin .gt. jmax)  go to 12
        mu = iju(k) - jmin
        do 11 j=jmin,jmax
  11      u(j) = row(ju(mu+j)) * dk
  12    continue
c
c  ******  update irl and jrl, keeping jrl in decreasing order  ********
        i = i1
        if (i .eq. 0) go to 18
  14    irl(i) = irl(i) + 1
        i1 = jrl(i)
        if (irl(i) .ge. il(i+1)) go to 17
        ijlb = irl(i) - il(i) + ijl(i)
        j = jl(ijlb)
  15    if (i .gt. jrl(j)) go to 16
          j = jrl(j)
          go to 15
  16    jrl(i) = jrl(j)
        jrl(j) = i
  17    i = i1
        if (i .ne. 0) go to 14
  18    if (irl(k) .ge. il(k+1)) go to 19
        j = jl(ijl(k))
        jrl(k) = jrl(j)
        jrl(j) = k
  19    continue
c
c  ******  solve  ux = tmp  by back substitution  **********************
      k = n
      do 22 i=1,n
        sum =  tmp(k)
        jmin = iu(k)
        jmax = iu(k+1) - 1
        if (jmin .gt. jmax)  go to 21
        mu = iju(k) - jmin
        do 20 j=jmin,jmax
  20      sum = sum - u(j) * tmp(ju(mu+j))
  21    tmp(k) =  sum
        z(c(k)) =  sum
  22    k = k-1
      flag = 0
      return
c
c ** error.. insufficient storage for l
 104  flag = 4*n + 1
      return
c ** error.. insufficient storage for u
 107  flag = 7*n + 1
      return
c ** error.. zero pivot
 108  flag = 8*n + k
      return
      end
      subroutine nnsc
     *     (n, r, c, il, jl, ijl, l, d, iu, ju, iju, u, z, b, tmp)
c*** subroutine nnsc
c*** numerical solution of sparse nonsymmetric system of linear
c      equations given ldu-factorization (compressed pointer storage)
c
c
c       input variables..  n, r, c, il, jl, ijl, l, d, iu, ju, iju, u, b
c       output variables.. z
c
c       parameters used internally..
c fia   - tmp   - temporary vector which gets result of solving  ly = b.
c       -           size = n.
c
c  internal variables..
c    jmin, jmax - indices of the first and last positions in a row of
c      u or l  to be used.
c
      integer r(*), c(*), il(*), jl(*), ijl(*), iu(*), ju(*), iju(*)
c     real l(*), d(*), u(*), b(*), z(*), tmp(*), tmpk, sum
      double precision  l(*), d(*), u(*), b(*), z(*), tmp(*), tmpk,sum
c
c  ******  set tmp to reordered b  *************************************
      do 1 k=1,n
   1    tmp(k) = b(r(k))
c  ******  solve  ly = b  by forward substitution  *********************
      do 3 k=1,n
        jmin = il(k)
        jmax = il(k+1) - 1
        tmpk = -d(k) * tmp(k)
        tmp(k) = -tmpk
        if (jmin .gt. jmax) go to 3
        ml = ijl(k) - jmin
        do 2 j=jmin,jmax
   2      tmp(jl(ml+j)) = tmp(jl(ml+j)) + tmpk * l(j)
   3    continue
c  ******  solve  ux = y  by back substitution  ************************
      k = n
      do 6 i=1,n
        sum = -tmp(k)
        jmin = iu(k)
        jmax = iu(k+1) - 1
        if (jmin .gt. jmax) go to 5
        mu = iju(k) - jmin
        do 4 j=jmin,jmax
   4      sum = sum + u(j) * tmp(ju(mu+j))
   5    tmp(k) = -sum
        z(c(k)) = -sum
        k = k - 1
   6    continue
      return
      end
      subroutine nntc
     *     (n, r, c, il, jl, ijl, l, d, iu, ju, iju, u, z, b, tmp)
c*** subroutine nntc
c*** numeric solution of the transpose of a sparse nonsymmetric system
c      of linear equations given lu-factorization (compressed pointer
c      storage)
c
c
c       input variables..  n, r, c, il, jl, ijl, l, d, iu, ju, iju, u, b
c       output variables.. z
c
c       parameters used internally..
c fia   - tmp   - temporary vector which gets result of solving ut y = b
c       -           size = n.
c
c  internal variables..
c    jmin, jmax - indices of the first and last positions in a row of
c      u or l  to be used.
c
      integer r(*), c(*), il(*), jl(*), ijl(*), iu(*), ju(*), iju(*)
c     real l(*), d(*), u(*), b(*), z(*), tmp(*), tmpk,sum
      double precision l(*), d(*), u(*), b(*), z(*), tmp(*), tmpk,sum
c
c  ******  set tmp to reordered b  *************************************
      do 1 k=1,n
   1    tmp(k) = b(c(k))
c  ******  solve  ut y = b  by forward substitution  *******************
      do 3 k=1,n
        jmin = iu(k)
        jmax = iu(k+1) - 1
        tmpk = -tmp(k)
        if (jmin .gt. jmax) go to 3
        mu = iju(k) - jmin
        do 2 j=jmin,jmax
   2      tmp(ju(mu+j)) = tmp(ju(mu+j)) + tmpk * u(j)
   3    continue
c  ******  solve  lt x = y  by back substitution  **********************
      k = n
      do 6 i=1,n
        sum = -tmp(k)
        jmin = il(k)
        jmax = il(k+1) - 1
        if (jmin .gt. jmax) go to 5
        ml = ijl(k) - jmin
        do 4 j=jmin,jmax
   4      sum = sum + l(j) * tmp(jl(ml+j))
   5    tmp(k) = -sum * d(k)
        z(r(k)) = tmp(k)
        k = k - 1
   6    continue
      return
      end
*DECK DSTODA
      SUBROUTINE DSTODA (NEQ, Y, YH, NYH, YH1, EWT, SAVF, ACOR,
     1   WM, IWM, F, JAC, PJAC, SLVS)
      EXTERNAL F, JAC, PJAC, SLVS
      INTEGER NEQ, NYH, IWM
      INTEGER ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, 
     1   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     2   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      INTEGER JTYP, MUSED, MXORDN, MXORDS
      INTEGER I, I1, ICOUNT, IRFLAG, IREDO, IRET, J, JB, M, NCF, NEWQ
      INTEGER LM1, LM1P1, LM2, LM2P1, NQM1, NQM2
      INTEGER IALTH, IPUP, LMAX, NQNYH, NSLP
      DOUBLE PRECISION Y, YH, YH1, EWT, SAVF, ACOR, WM
      DOUBLE PRECISION CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
      DOUBLE PRECISION PDNORM
      DOUBLE PRECISION CONIT, CRATE, EL(13), ELCO(13,12), HOLD, RMAX,
     1   TESCO(3,12)
      DOUBLE PRECISION DCON, DDN, DEL, DELP, DSM, DUP, EXDN, EXSM, EXUP,
     1   R, RH, RHDN, RHSM, RHUP, TOLD, DMNORM
      DOUBLE PRECISION ALPHA, CM1(12),CM2(5), DM1,DM2, EXM1,EXM2, PDEST,
     1   PDLAST, PDH, PNORM, RATE, RATIO, RH1, RH1IT, RH2, RM, SM1(12)
      DIMENSION NEQ(*), Y(*), YH(NYH,*), YH1(*), EWT(*), SAVF(*),
     1   ACOR(*), WM(*), IWM(*)
      SAVE CONIT, CRATE, EL, ELCO, HOLD, RMAX, TESCO,
     1   IALTH, IPUP, LMAX, NQNYH, NSLP
      SAVE CM1, CM2, PDEST, PDLAST, RATIO, SM1, ICOUNT, IRFLAG
      COMMON /DLS001/ CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,
     1   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, 
     2   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     3   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      COMMON /DLSA01/ PDNORM, JTYP, MUSED, MXORDN, MXORDS
      DATA SM1/0.5D0, 0.575D0, 0.55D0, 0.45D0, 0.35D0, 0.25D0,
     1   0.20D0, 0.15D0, 0.10D0, 0.075D0, 0.050D0, 0.025D0/
C-----------------------------------------------------------------------
C DSTODA performs one step of the integration of an initial value
C problem for a system of ordinary differential equations.
C Note: DSTODA is independent of the value of the iteration method
C indicator MITER, when this is .ne. 0, and hence is independent
C of the type of chord method used, or the Jacobian structure.
C Communication with DSTODA is done with the following variables:
C
C Y      = an array of length .ge. N used as the Y argument in
C          all calls to F and JAC.
C NEQ    = integer array containing problem size in NEQ(1), and
C          passed as the NEQ argument in all calls to F and JAC.
C YH     = an NYH by LMAX array containing the dependent variables
C          and their approximate scaled derivatives, where
C          LMAX = MAXORD + 1.  YH(i,j+1) contains the approximate
C          j-th derivative of y(i), scaled by H**j/factorial(j)
C          (j = 0,1,...,NQ).  On entry for the first step, the first
C          two columns of YH must be set from the initial values.
C NYH    = a constant integer .ge. N, the first dimension of YH.
C YH1    = a one-dimensional array occupying the same space as YH.
C EWT    = an array of length N containing multiplicative weights
C          for local error measurements.  Local errors in y(i) are
C          compared to 1.0/EWT(i) in various error tests.
C SAVF   = an array of working storage, of length N.
C ACOR   = a work array of length N, used for the accumulated
C          corrections.  On a successful return, ACOR(i) contains
C          the estimated one-step local error in y(i).
C WM,IWM = real and integer work arrays associated with matrix
C          operations in chord iteration (MITER .ne. 0).
C PJAC   = name of routine to evaluate and preprocess Jacobian matrix
C          and P = I - H*EL0*Jac, if a chord method is being used.
C          It also returns an estimate of norm(Jac) in PDNORM.
C SLVS   = name of routine to solve linear system in chord iteration.
C CCMAX  = maximum relative change in H*EL0 before PJAC is called.
C H      = the step size to be attempted on the next step.
C          H is altered by the error control algorithm during the
C          problem.  H can be either positive or negative, but its
C          sign must remain constant throughout the problem.
C HMIN   = the minimum absolute value of the step size H to be used.
C HMXI   = inverse of the maximum absolute value of H to be used.
C          HMXI = 0.0 is allowed and corresponds to an infinite HMAX.
C          HMIN and HMXI may be changed at any time, but will not
C          take effect until the next change of H is considered.
C TN     = the independent variable. TN is updated on each step taken.
C JSTART = an integer used for input only, with the following
C          values and meanings:
C               0  perform the first step.
C           .gt.0  take a new step continuing from the last.
C              -1  take the next step with a new value of H,
C                    N, METH, MITER, and/or matrix parameters.
C              -2  take the next step with a new value of H,
C                    but with other inputs unchanged.
C          On return, JSTART is set to 1 to facilitate continuation.
C KFLAG  = a completion code with the following meanings:
C               0  the step was succesful.
C              -1  the requested error could not be achieved.
C              -2  corrector convergence could not be achieved.
C              -3  fatal error in PJAC or SLVS.
C          A return with KFLAG = -1 or -2 means either
C          ABS(H) = HMIN or 10 consecutive failures occurred.
C          On a return with KFLAG negative, the values of TN and
C          the YH array are as of the beginning of the last
C          step, and H is the last step size attempted.
C MAXORD = the maximum order of integration method to be allowed.
C MAXCOR = the maximum number of corrector iterations allowed.
C MSBP   = maximum number of steps between PJAC calls (MITER .gt. 0).
C MXNCF  = maximum number of convergence failures allowed.
C METH   = current method.
C          METH = 1 means Adams method (nonstiff)
C          METH = 2 means BDF method (stiff)
C          METH may be reset by DSTODA.
C MITER  = corrector iteration method.
C          MITER = 0 means functional iteration.
C          MITER = JT .gt. 0 means a chord iteration corresponding
C          to Jacobian type JT.  (The DLSODA/DLSODAR argument JT is
C          communicated here as JTYP, but is not used in DSTODA
C          except to load MITER following a method switch.)
C          MITER may be reset by DSTODA.
C N      = the number of first-order differential equations.
C-----------------------------------------------------------------------
      KFLAG = 0
      TOLD = TN
      NCF = 0
      IERPJ = 0
      IERSL = 0
      JCUR = 0
      ICF = 0
      DELP = 0.0D0
      IF (JSTART .GT. 0) GO TO 200
      IF (JSTART .EQ. -1) GO TO 100
      IF (JSTART .EQ. -2) GO TO 160
C-----------------------------------------------------------------------
C On the first call, the order is set to 1, and other variables are
C initialized.  RMAX is the maximum ratio by which H can be increased
C in a single step.  It is initially 1.E4 to compensate for the small
C initial H, but then is normally equal to 10.  If a failure
C occurs (in corrector convergence or error test), RMAX is set at 2
C for the next increase.
C DCFODE is called to get the needed coefficients for both methods.
C-----------------------------------------------------------------------
      LMAX = MAXORD + 1
      NQ = 1
      L = 2
      IALTH = 2
      RMAX = 10000.0D0
      RC = 0.0D0
      EL0 = 1.0D0
      CRATE = 0.7D0
      HOLD = H
      NSLP = 0
      IPUP = MITER
      IRET = 3
C Initialize switching parameters.  METH = 1 is assumed initially. -----
      ICOUNT = 20
      IRFLAG = 0
      PDEST = 0.0D0
      PDLAST = 0.0D0
      RATIO = 5.0D0
      CALL DCFODE (2, ELCO, TESCO)
      DO 10 I = 1,5
 10     CM2(I) = TESCO(2,I)*ELCO(I+1,I)
      CALL DCFODE (1, ELCO, TESCO)
      DO 20 I = 1,12
 20     CM1(I) = TESCO(2,I)*ELCO(I+1,I)
      GO TO 150
C-----------------------------------------------------------------------
C The following block handles preliminaries needed when JSTART = -1.
C IPUP is set to MITER to force a matrix update.
C If an order increase is about to be considered (IALTH = 1),
C IALTH is reset to 2 to postpone consideration one more step.
C If the caller has changed METH, DCFODE is called to reset
C the coefficients of the method.
C If H is to be changed, YH must be rescaled.
C If H or METH is being changed, IALTH is reset to L = NQ + 1
C to prevent further changes in H for that many steps.
C-----------------------------------------------------------------------
 100  IPUP = MITER
      LMAX = MAXORD + 1
      IF (IALTH .EQ. 1) IALTH = 2
      IF (METH .EQ. MUSED) GO TO 160
      CALL DCFODE (METH, ELCO, TESCO)
      IALTH = L
      IRET = 1
C-----------------------------------------------------------------------
C The el vector and related constants are reset
C whenever the order NQ is changed, or at the start of the problem.
C-----------------------------------------------------------------------
 150  DO 155 I = 1,L
 155    EL(I) = ELCO(I,NQ)
      NQNYH = NQ*NYH
      RC = RC*EL(1)/EL0
      EL0 = EL(1)
      CONIT = 0.5D0/(NQ+2)
      GO TO (160, 170, 200), IRET
C-----------------------------------------------------------------------
C If H is being changed, the H ratio RH is checked against
C RMAX, HMIN, and HMXI, and the YH array rescaled.  IALTH is set to
C L = NQ + 1 to prevent a change of H for that many steps, unless
C forced by a convergence or error test failure.
C-----------------------------------------------------------------------
 160  IF (H .EQ. HOLD) GO TO 200
      RH = H/HOLD
      H = HOLD
      IREDO = 3
      GO TO 175
 170  RH = MAX(RH,HMIN/ABS(H))
 175  RH = MIN(RH,RMAX)
      RH = RH/MAX(1.0D0,ABS(H)*HMXI*RH)
C-----------------------------------------------------------------------
C If METH = 1, also restrict the new step size by the stability region.
C If this reduces H, set IRFLAG to 1 so that if there are roundoff
C problems later, we can assume that is the cause of the trouble.
C-----------------------------------------------------------------------
      IF (METH .EQ. 2) GO TO 178
      IRFLAG = 0
      PDH = MAX(ABS(H)*PDLAST,0.000001D0)
      IF (RH*PDH*1.00001D0 .LT. SM1(NQ)) GO TO 178
      RH = SM1(NQ)/PDH
      IRFLAG = 1
 178  CONTINUE
      R = 1.0D0
      DO 180 J = 2,L
        R = R*RH
        DO 180 I = 1,N
 180      YH(I,J) = YH(I,J)*R
      H = H*RH
      RC = RC*RH
      IALTH = L
      IF (IREDO .EQ. 0) GO TO 690
C-----------------------------------------------------------------------
C This section computes the predicted values by effectively
C multiplying the YH array by the Pascal triangle matrix.
C RC is the ratio of new to old values of the coefficient  H*EL(1).
C When RC differs from 1 by more than CCMAX, IPUP is set to MITER
C to force PJAC to be called, if a Jacobian is involved.
C In any case, PJAC is called at least every MSBP steps.
C-----------------------------------------------------------------------
 200  IF (ABS(RC-1.0D0) .GT. CCMAX) IPUP = MITER
      IF (NST .GE. NSLP+MSBP) IPUP = MITER
      TN = TN + H
      I1 = NQNYH + 1
      DO 215 JB = 1,NQ
        I1 = I1 - NYH
CDIR$ IVDEP
        DO 210 I = I1,NQNYH
 210      YH1(I) = YH1(I) + YH1(I+NYH)
 215    CONTINUE
      PNORM = DMNORM (N, YH1, EWT)
C-----------------------------------------------------------------------
C Up to MAXCOR corrector iterations are taken.  A convergence test is
C made on the RMS-norm of each correction, weighted by the error
C weight vector EWT.  The sum of the corrections is accumulated in the
C vector ACOR(i).  The YH array is not altered in the corrector loop.
C-----------------------------------------------------------------------
 220  M = 0
      RATE = 0.0D0
      DEL = 0.0D0
      DO 230 I = 1,N
 230    Y(I) = YH(I,1)
      CALL F (NEQ, TN, Y, SAVF)
      NFE = NFE + 1
      IF (IPUP .LE. 0) GO TO 250
C-----------------------------------------------------------------------
C If indicated, the matrix P = I - H*EL(1)*J is reevaluated and
C preprocessed before starting the corrector iteration.  IPUP is set
C to 0 as an indicator that this has been done.
C-----------------------------------------------------------------------
      CALL PJAC (NEQ, Y, YH, NYH, EWT, ACOR, SAVF, WM, IWM, F, JAC)
      IPUP = 0
      RC = 1.0D0
      NSLP = NST
      CRATE = 0.7D0
      IF (IERPJ .NE. 0) GO TO 430
 250  DO 260 I = 1,N
 260    ACOR(I) = 0.0D0
 270  IF (MITER .NE. 0) GO TO 350
C-----------------------------------------------------------------------
C In the case of functional iteration, update Y directly from
C the result of the last function evaluation.
C-----------------------------------------------------------------------
      DO 290 I = 1,N
        SAVF(I) = H*SAVF(I) - YH(I,2)
 290    Y(I) = SAVF(I) - ACOR(I)
      DEL = DMNORM (N, Y, EWT)
      DO 300 I = 1,N
        Y(I) = YH(I,1) + EL(1)*SAVF(I)
 300    ACOR(I) = SAVF(I)
      GO TO 400
C-----------------------------------------------------------------------
C In the case of the chord method, compute the corrector error,
C and solve the linear system with that as right-hand side and
C P as coefficient matrix.
C-----------------------------------------------------------------------
 350  DO 360 I = 1,N
 360    Y(I) = H*SAVF(I) - (YH(I,2) + ACOR(I))
      CALL SLVS (WM, IWM, Y, SAVF)
      IF (IERSL .LT. 0) GO TO 430
      IF (IERSL .GT. 0) GO TO 410
      DEL = DMNORM (N, Y, EWT)
      DO 380 I = 1,N
        ACOR(I) = ACOR(I) + Y(I)
 380    Y(I) = YH(I,1) + EL(1)*ACOR(I)
C-----------------------------------------------------------------------
C Test for convergence.  If M .gt. 0, an estimate of the convergence
C rate constant is stored in CRATE, and this is used in the test.
C
C We first check for a change of iterates that is the size of
C roundoff error.  If this occurs, the iteration has converged, and a
C new rate estimate is not formed.
C In all other cases, force at least two iterations to estimate a
C local Lipschitz constant estimate for Adams methods.
C On convergence, form PDEST = local maximum Lipschitz constant
C estimate.  PDLAST is the most recent nonzero estimate.
C-----------------------------------------------------------------------
 400  CONTINUE
      IF (DEL .LE. 100.0D0*PNORM*UROUND) GO TO 450
      IF (M .EQ. 0 .AND. METH .EQ. 1) GO TO 405
      IF (M .EQ. 0) GO TO 402
      RM = 1024.0D0
      IF (DEL .LE. 1024.0D0*DELP) RM = DEL/DELP
      RATE = MAX(RATE,RM)
      CRATE = MAX(0.2D0*CRATE,RM)
 402  DCON = DEL*MIN(1.0D0,1.5D0*CRATE)/(TESCO(2,NQ)*CONIT)
      IF (DCON .GT. 1.0D0) GO TO 405
      PDEST = MAX(PDEST,RATE/ABS(H*EL(1)))
      IF (PDEST .NE. 0.0D0) PDLAST = PDEST
      GO TO 450
 405  CONTINUE
      M = M + 1
      IF (M .EQ. MAXCOR) GO TO 410
      IF (M .GE. 2 .AND. DEL .GT. 2.0D0*DELP) GO TO 410
      DELP = DEL
      CALL F (NEQ, TN, Y, SAVF)
      NFE = NFE + 1
      GO TO 270
C-----------------------------------------------------------------------
C The corrector iteration failed to converge.
C If MITER .ne. 0 and the Jacobian is out of date, PJAC is called for
C the next try.  Otherwise the YH array is retracted to its values
C before prediction, and H is reduced, if possible.  If H cannot be
C reduced or MXNCF failures have occurred, exit with KFLAG = -2.
C-----------------------------------------------------------------------
 410  IF (MITER .EQ. 0 .OR. JCUR .EQ. 1) GO TO 430
      ICF = 1
      IPUP = MITER
      GO TO 220
 430  ICF = 2
      NCF = NCF + 1
      RMAX = 2.0D0
      TN = TOLD
      I1 = NQNYH + 1
      DO 445 JB = 1,NQ
        I1 = I1 - NYH
CDIR$ IVDEP
        DO 440 I = I1,NQNYH
 440      YH1(I) = YH1(I) - YH1(I+NYH)
 445    CONTINUE
      IF (IERPJ .LT. 0 .OR. IERSL .LT. 0) GO TO 680
      IF (ABS(H) .LE. HMIN*1.00001D0) GO TO 670
      IF (NCF .EQ. MXNCF) GO TO 670
      RH = 0.25D0
      IPUP = MITER
      IREDO = 1
      GO TO 170
C-----------------------------------------------------------------------
C The corrector has converged.  JCUR is set to 0
C to signal that the Jacobian involved may need updating later.
C The local error test is made and control passes to statement 500
C if it fails.
C-----------------------------------------------------------------------
 450  JCUR = 0
      IF (M .EQ. 0) DSM = DEL/TESCO(2,NQ)
      IF (M .GT. 0) DSM = DMNORM (N, ACOR, EWT)/TESCO(2,NQ)
      IF (DSM .GT. 1.0D0) GO TO 500
C-----------------------------------------------------------------------
C After a successful step, update the YH array.
C Decrease ICOUNT by 1, and if it is -1, consider switching methods.
C If a method switch is made, reset various parameters,
C rescale the YH array, and exit.  If there is no switch,
C consider changing H if IALTH = 1.  Otherwise decrease IALTH by 1.
C If IALTH is then 1 and NQ .lt. MAXORD, then ACOR is saved for
C use in a possible order increase on the next step.
C If a change in H is considered, an increase or decrease in order
C by one is considered also.  A change in H is made only if it is by a
C factor of at least 1.1.  If not, IALTH is set to 3 to prevent
C testing for that many steps.
C-----------------------------------------------------------------------
      KFLAG = 0
      IREDO = 0
      NST = NST + 1
      HU = H
      NQU = NQ
      MUSED = METH
      DO 460 J = 1,L
        DO 460 I = 1,N
 460      YH(I,J) = YH(I,J) + EL(J)*ACOR(I)
      ICOUNT = ICOUNT - 1
      IF (ICOUNT .GE. 0) GO TO 488
      IF (METH .EQ. 2) GO TO 480
C-----------------------------------------------------------------------
C We are currently using an Adams method.  Consider switching to BDF.
C If the current order is greater than 5, assume the problem is
C not stiff, and skip this section.
C If the Lipschitz constant and error estimate are not polluted
C by roundoff, go to 470 and perform the usual test.
C Otherwise, switch to the BDF methods if the last step was
C restricted to insure stability (irflag = 1), and stay with Adams
C method if not.  When switching to BDF with polluted error estimates,
C in the absence of other information, double the step size.
C
C When the estimates are OK, we make the usual test by computing
C the step size we could have (ideally) used on this step,
C with the current (Adams) method, and also that for the BDF.
C If NQ .gt. MXORDS, we consider changing to order MXORDS on switching.
C Compare the two step sizes to decide whether to switch.
C The step size advantage must be at least RATIO = 5 to switch.
C-----------------------------------------------------------------------
      IF (NQ .GT. 5) GO TO 488
      IF (DSM .GT. 100.0D0*PNORM*UROUND .AND. PDEST .NE. 0.0D0)
     1   GO TO 470
      IF (IRFLAG .EQ. 0) GO TO 488
      RH2 = 2.0D0
      NQM2 = MIN(NQ,MXORDS)
      GO TO 478
 470  CONTINUE
      EXSM = 1.0D0/L
      RH1 = 1.0D0/(1.2D0*DSM**EXSM + 0.0000012D0)
      RH1IT = 2.0D0*RH1
      PDH = PDLAST*ABS(H)
      IF (PDH*RH1 .GT. 0.00001D0) RH1IT = SM1(NQ)/PDH
      RH1 = MIN(RH1,RH1IT)
      IF (NQ .LE. MXORDS) GO TO 474
         NQM2 = MXORDS
         LM2 = MXORDS + 1
         EXM2 = 1.0D0/LM2
         LM2P1 = LM2 + 1
         DM2 = DMNORM (N, YH(1,LM2P1), EWT)/CM2(MXORDS)
         RH2 = 1.0D0/(1.2D0*DM2**EXM2 + 0.0000012D0)
         GO TO 476
 474  DM2 = DSM*(CM1(NQ)/CM2(NQ))
      RH2 = 1.0D0/(1.2D0*DM2**EXSM + 0.0000012D0)
      NQM2 = NQ
 476  CONTINUE
      IF (RH2 .LT. RATIO*RH1) GO TO 488
C THE SWITCH TEST PASSED.  RESET RELEVANT QUANTITIES FOR BDF. ----------
 478  RH = RH2
      ICOUNT = 20
      METH = 2
      MITER = JTYP
      PDLAST = 0.0D0
      NQ = NQM2
      L = NQ + 1
      GO TO 170
C-----------------------------------------------------------------------
C We are currently using a BDF method.  Consider switching to Adams.
C Compute the step size we could have (ideally) used on this step,
C with the current (BDF) method, and also that for the Adams.
C If NQ .gt. MXORDN, we consider changing to order MXORDN on switching.
C Compare the two step sizes to decide whether to switch.
C The step size advantage must be at least 5/RATIO = 1 to switch.
C If the step size for Adams would be so small as to cause
C roundoff pollution, we stay with BDF.
C-----------------------------------------------------------------------
 480  CONTINUE
      EXSM = 1.0D0/L
      IF (MXORDN .GE. NQ) GO TO 484
         NQM1 = MXORDN
         LM1 = MXORDN + 1
         EXM1 = 1.0D0/LM1
         LM1P1 = LM1 + 1
         DM1 = DMNORM (N, YH(1,LM1P1), EWT)/CM1(MXORDN)
         RH1 = 1.0D0/(1.2D0*DM1**EXM1 + 0.0000012D0)
         GO TO 486
 484  DM1 = DSM*(CM2(NQ)/CM1(NQ))
      RH1 = 1.0D0/(1.2D0*DM1**EXSM + 0.0000012D0)
      NQM1 = NQ
      EXM1 = EXSM
 486  RH1IT = 2.0D0*RH1
      PDH = PDNORM*ABS(H)
      IF (PDH*RH1 .GT. 0.00001D0) RH1IT = SM1(NQM1)/PDH
      RH1 = MIN(RH1,RH1IT)
      RH2 = 1.0D0/(1.2D0*DSM**EXSM + 0.0000012D0)
      IF (RH1*RATIO .LT. 5.0D0*RH2) GO TO 488
      ALPHA = MAX(0.001D0,RH1)
      DM1 = (ALPHA**EXM1)*DM1
      IF (DM1 .LE. 1000.0D0*UROUND*PNORM) GO TO 488
C The switch test passed.  Reset relevant quantities for Adams. --------
      RH = RH1
      ICOUNT = 20
      METH = 1
      MITER = 0
      PDLAST = 0.0D0
      NQ = NQM1
      L = NQ + 1
      GO TO 170
C
C No method switch is being made.  Do the usual step/order selection. --
 488  CONTINUE
      IALTH = IALTH - 1
      IF (IALTH .EQ. 0) GO TO 520
      IF (IALTH .GT. 1) GO TO 700
      IF (L .EQ. LMAX) GO TO 700
      DO 490 I = 1,N
 490    YH(I,LMAX) = ACOR(I)
      GO TO 700
C-----------------------------------------------------------------------
C The error test failed.  KFLAG keeps track of multiple failures.
C Restore TN and the YH array to their previous values, and prepare
C to try the step again.  Compute the optimum step size for this or
C one lower order.  After 2 or more failures, H is forced to decrease
C by a factor of 0.2 or less.
C-----------------------------------------------------------------------
 500  KFLAG = KFLAG - 1
      TN = TOLD
      I1 = NQNYH + 1
      DO 515 JB = 1,NQ
        I1 = I1 - NYH
CDIR$ IVDEP
        DO 510 I = I1,NQNYH
 510      YH1(I) = YH1(I) - YH1(I+NYH)
 515    CONTINUE
      RMAX = 2.0D0
      IF (ABS(H) .LE. HMIN*1.00001D0) GO TO 660
      IF (KFLAG .LE. -3) GO TO 640
      IREDO = 2
      RHUP = 0.0D0
      GO TO 540
C-----------------------------------------------------------------------
C Regardless of the success or failure of the step, factors
C RHDN, RHSM, and RHUP are computed, by which H could be multiplied
C at order NQ - 1, order NQ, or order NQ + 1, respectively.
C In the case of failure, RHUP = 0.0 to avoid an order increase.
C The largest of these is determined and the new order chosen
C accordingly.  If the order is to be increased, we compute one
C additional scaled derivative.
C-----------------------------------------------------------------------
 520  RHUP = 0.0D0
      IF (L .EQ. LMAX) GO TO 540
      DO 530 I = 1,N
 530    SAVF(I) = ACOR(I) - YH(I,LMAX)
      DUP = DMNORM (N, SAVF, EWT)/TESCO(3,NQ)
      EXUP = 1.0D0/(L+1)
      RHUP = 1.0D0/(1.4D0*DUP**EXUP + 0.0000014D0)
 540  EXSM = 1.0D0/L
      RHSM = 1.0D0/(1.2D0*DSM**EXSM + 0.0000012D0)
      RHDN = 0.0D0
      IF (NQ .EQ. 1) GO TO 550
      DDN = DMNORM (N, YH(1,L), EWT)/TESCO(1,NQ)
      EXDN = 1.0D0/NQ
      RHDN = 1.0D0/(1.3D0*DDN**EXDN + 0.0000013D0)
C If METH = 1, limit RH according to the stability region also. --------
 550  IF (METH .EQ. 2) GO TO 560
      PDH = MAX(ABS(H)*PDLAST,0.000001D0)
      IF (L .LT. LMAX) RHUP = MIN(RHUP,SM1(L)/PDH)
      RHSM = MIN(RHSM,SM1(NQ)/PDH)
      IF (NQ .GT. 1) RHDN = MIN(RHDN,SM1(NQ-1)/PDH)
      PDEST = 0.0D0
 560  IF (RHSM .GE. RHUP) GO TO 570
      IF (RHUP .GT. RHDN) GO TO 590
      GO TO 580
 570  IF (RHSM .LT. RHDN) GO TO 580
      NEWQ = NQ
      RH = RHSM
      GO TO 620
 580  NEWQ = NQ - 1
      RH = RHDN
      IF (KFLAG .LT. 0 .AND. RH .GT. 1.0D0) RH = 1.0D0
      GO TO 620
 590  NEWQ = L
      RH = RHUP
      IF (RH .LT. 1.1D0) GO TO 610
      R = EL(L)/L
      DO 600 I = 1,N
 600    YH(I,NEWQ+1) = ACOR(I)*R
      GO TO 630
 610  IALTH = 3
      GO TO 700
C If METH = 1 and H is restricted by stability, bypass 10 percent test.
 620  IF (METH .EQ. 2) GO TO 622
      IF (RH*PDH*1.00001D0 .GE. SM1(NEWQ)) GO TO 625
 622  IF (KFLAG .EQ. 0 .AND. RH .LT. 1.1D0) GO TO 610
 625  IF (KFLAG .LE. -2) RH = MIN(RH,0.2D0)
C-----------------------------------------------------------------------
C If there is a change of order, reset NQ, L, and the coefficients.
C In any case H is reset according to RH and the YH array is rescaled.
C Then exit from 690 if the step was OK, or redo the step otherwise.
C-----------------------------------------------------------------------
      IF (NEWQ .EQ. NQ) GO TO 170
 630  NQ = NEWQ
      L = NQ + 1
      IRET = 2
      GO TO 150
C-----------------------------------------------------------------------
C Control reaches this section if 3 or more failures have occured.
C If 10 failures have occurred, exit with KFLAG = -1.
C It is assumed that the derivatives that have accumulated in the
C YH array have errors of the wrong order.  Hence the first
C derivative is recomputed, and the order is set to 1.  Then
C H is reduced by a factor of 10, and the step is retried,
C until it succeeds or H reaches HMIN.
C-----------------------------------------------------------------------
 640  IF (KFLAG .EQ. -10) GO TO 660
      RH = 0.1D0
      RH = MAX(HMIN/ABS(H),RH)
      H = H*RH
      DO 645 I = 1,N
 645    Y(I) = YH(I,1)
      CALL F (NEQ, TN, Y, SAVF)
      NFE = NFE + 1
      DO 650 I = 1,N
 650    YH(I,2) = H*SAVF(I)
      IPUP = MITER
      IALTH = 5
      IF (NQ .EQ. 1) GO TO 200
      NQ = 1
      L = 2
      IRET = 3
      GO TO 150
C-----------------------------------------------------------------------
C All returns are made through this section.  H is saved in HOLD
C to allow the caller to change H on the next step.
C-----------------------------------------------------------------------
 660  KFLAG = -1
      GO TO 720
 670  KFLAG = -2
      GO TO 720
 680  KFLAG = -3
      GO TO 720
 690  RMAX = 10.0D0
 700  R = 1.0D0/TESCO(2,NQU)
      DO 710 I = 1,N
 710    ACOR(I) = ACOR(I)*R
 720  HOLD = H
      JSTART = 1
      RETURN
C----------------------- End of Subroutine DSTODA ----------------------
      END
*DECK DPRJA
      SUBROUTINE DPRJA (NEQ, Y, YH, NYH, EWT, FTEM, SAVF, WM, IWM,
     1   F, JAC)
      EXTERNAL F, JAC
      INTEGER NEQ, NYH, IWM
      INTEGER ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, 
     1   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     2   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      INTEGER JTYP, MUSED, MXORDN, MXORDS
      INTEGER I, I1, I2, IER, II, J, J1, JJ, LENP,
     1   MBA, MBAND, MEB1, MEBAND, ML, ML3, MU, NP1
      DOUBLE PRECISION Y, YH, EWT, FTEM, SAVF, WM
      DOUBLE PRECISION CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
      DOUBLE PRECISION PDNORM
      DOUBLE PRECISION CON, FAC, HL0, R, R0, SRUR, YI, YJ, YJJ,
     1   DMNORM, DFNORM, DBNORM
      DIMENSION NEQ(*), Y(*), YH(NYH,*), EWT(*), FTEM(*), SAVF(*),
     1   WM(*), IWM(*)
      COMMON /DLS001/ CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,
     1   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, 
     2   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     3   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      COMMON /DLSA01/ PDNORM, JTYP, MUSED, MXORDN, MXORDS
C-----------------------------------------------------------------------
C DPRJA is called by DSTODA to compute and process the matrix
C P = I - H*EL(1)*J , where J is an approximation to the Jacobian.
C Here J is computed by the user-supplied routine JAC if
C MITER = 1 or 4 or by finite differencing if MITER = 2 or 5.
C J, scaled by -H*EL(1), is stored in WM.  Then the norm of J (the
C matrix norm consistent with the weighted max-norm on vectors given
C by DMNORM) is computed, and J is overwritten by P.  P is then
C subjected to LU decomposition in preparation for later solution
C of linear systems with P as coefficient matrix.  This is done
C by DGEFA if MITER = 1 or 2, and by DGBFA if MITER = 4 or 5.
C
C In addition to variables described previously, communication
C with DPRJA uses the following:
C Y     = array containing predicted values on entry.
C FTEM  = work array of length N (ACOR in DSTODA).
C SAVF  = array containing f evaluated at predicted y.
C WM    = real work space for matrices.  On output it contains the
C         LU decomposition of P.
C         Storage of matrix elements starts at WM(3).
C         WM also contains the following matrix-related data:
C         WM(1) = SQRT(UROUND), used in numerical Jacobian increments.
C IWM   = integer work space containing pivot information, starting at
C         IWM(21).   IWM also contains the band parameters
C         ML = IWM(1) and MU = IWM(2) if MITER is 4 or 5.
C EL0   = EL(1) (input).
C PDNORM= norm of Jacobian matrix. (Output).
C IERPJ = output error flag,  = 0 if no trouble, .gt. 0 if
C         P matrix found to be singular.
C JCUR  = output flag = 1 to indicate that the Jacobian matrix
C         (or approximation) is now current.
C This routine also uses the Common variables EL0, H, TN, UROUND,
C MITER, N, NFE, and NJE.
C-----------------------------------------------------------------------
      NJE = NJE + 1
      IERPJ = 0
      JCUR = 1
      HL0 = H*EL0
      GO TO (100, 200, 300, 400, 500), MITER
C If MITER = 1, call JAC and multiply by scalar. -----------------------
 100  LENP = N*N
      DO 110 I = 1,LENP
 110    WM(I+2) = 0.0D0
      CALL JAC (NEQ, TN, Y, 0, 0, WM(3), N)
      CON = -HL0
      DO 120 I = 1,LENP
 120    WM(I+2) = WM(I+2)*CON
      GO TO 240
C If MITER = 2, make N calls to F to approximate J. --------------------
 200  FAC = DMNORM (N, SAVF, EWT)
      R0 = 1000.0D0*ABS(H)*UROUND*N*FAC
      IF (R0 .EQ. 0.0D0) R0 = 1.0D0
      SRUR = WM(1)
      J1 = 2
      DO 230 J = 1,N
        YJ = Y(J)
        R = MAX(SRUR*ABS(YJ),R0/EWT(J))
        Y(J) = Y(J) + R
        FAC = -HL0/R
        CALL F (NEQ, TN, Y, FTEM)
        DO 220 I = 1,N
 220      WM(I+J1) = (FTEM(I) - SAVF(I))*FAC
        Y(J) = YJ
        J1 = J1 + N
 230    CONTINUE
      NFE = NFE + N
 240  CONTINUE
C Compute norm of Jacobian. --------------------------------------------
      PDNORM = DFNORM (N, WM(3), EWT)/ABS(HL0)
C Add identity matrix. -------------------------------------------------
      J = 3
      NP1 = N + 1
      DO 250 I = 1,N
        WM(J) = WM(J) + 1.0D0
 250    J = J + NP1
C Do LU decomposition on P. --------------------------------------------
      CALL DGEFA (WM(3), N, N, IWM(21), IER)
      IF (IER .NE. 0) IERPJ = 1
      RETURN
C Dummy block only, since MITER is never 3 in this routine. ------------
 300  RETURN
C If MITER = 4, call JAC and multiply by scalar. -----------------------
 400  ML = IWM(1)
      MU = IWM(2)
      ML3 = ML + 3
      MBAND = ML + MU + 1
      MEBAND = MBAND + ML
      LENP = MEBAND*N
      DO 410 I = 1,LENP
 410    WM(I+2) = 0.0D0
      CALL JAC (NEQ, TN, Y, ML, MU, WM(ML3), MEBAND)
      CON = -HL0
      DO 420 I = 1,LENP
 420    WM(I+2) = WM(I+2)*CON
      GO TO 570
C If MITER = 5, make MBAND calls to F to approximate J. ----------------
 500  ML = IWM(1)
      MU = IWM(2)
      MBAND = ML + MU + 1
      MBA = MIN(MBAND,N)
      MEBAND = MBAND + ML
      MEB1 = MEBAND - 1
      SRUR = WM(1)
      FAC = DMNORM (N, SAVF, EWT)
      R0 = 1000.0D0*ABS(H)*UROUND*N*FAC
      IF (R0 .EQ. 0.0D0) R0 = 1.0D0
      DO 560 J = 1,MBA
        DO 530 I = J,N,MBAND
          YI = Y(I)
          R = MAX(SRUR*ABS(YI),R0/EWT(I))
 530      Y(I) = Y(I) + R
        CALL F (NEQ, TN, Y, FTEM)
        DO 550 JJ = J,N,MBAND
          Y(JJ) = YH(JJ,1)
          YJJ = Y(JJ)
          R = MAX(SRUR*ABS(YJJ),R0/EWT(JJ))
          FAC = -HL0/R
          I1 = MAX(JJ-MU,1)
          I2 = MIN(JJ+ML,N)
          II = JJ*MEB1 - ML + 2
          DO 540 I = I1,I2
 540        WM(II+I) = (FTEM(I) - SAVF(I))*FAC
 550      CONTINUE
 560    CONTINUE
      NFE = NFE + MBA
 570  CONTINUE
C Compute norm of Jacobian. --------------------------------------------
      PDNORM = DBNORM (N, WM(ML+3), MEBAND, ML, MU, EWT)/ABS(HL0)
C Add identity matrix. -------------------------------------------------
      II = MBAND + 2
      DO 580 I = 1,N
        WM(II) = WM(II) + 1.0D0
 580    II = II + MEBAND
C Do LU decomposition of P. --------------------------------------------
      CALL DGBFA (WM(3), MEBAND, N, ML, MU, IWM(21), IER)
      IF (IER .NE. 0) IERPJ = 1
      RETURN
C----------------------- End of Subroutine DPRJA -----------------------
      END
*DECK DMNORM
      DOUBLE PRECISION FUNCTION DMNORM (N, V, W)
C-----------------------------------------------------------------------
C This function routine computes the weighted max-norm
C of the vector of length N contained in the array V, with weights
C contained in the array w of length N:
C   DMNORM = MAX(i=1,...,N) ABS(V(i))*W(i)
C-----------------------------------------------------------------------
      INTEGER N,   I
      DOUBLE PRECISION V, W,   VM
      DIMENSION V(N), W(N)
      VM = 0.0D0
      DO 10 I = 1,N
 10     VM = MAX(VM,ABS(V(I))*W(I))
      DMNORM = VM
      RETURN
C----------------------- End of Function DMNORM ------------------------
      END
*DECK DFNORM
      DOUBLE PRECISION FUNCTION DFNORM (N, A, W)
C-----------------------------------------------------------------------
C This function computes the norm of a full N by N matrix,
C stored in the array A, that is consistent with the weighted max-norm
C on vectors, with weights stored in the array W:
C   DFNORM = MAX(i=1,...,N) ( W(i) * Sum(j=1,...,N) ABS(a(i,j))/W(j) )
C-----------------------------------------------------------------------
      INTEGER N,   I, J
      DOUBLE PRECISION A,   W, AN, SUM
      DIMENSION A(N,N), W(N)
      AN = 0.0D0
      DO 20 I = 1,N
        SUM = 0.0D0
        DO 10 J = 1,N
 10       SUM = SUM + ABS(A(I,J))/W(J)
        AN = MAX(AN,SUM*W(I))
 20     CONTINUE
      DFNORM = AN
      RETURN
C----------------------- End of Function DFNORM ------------------------
      END
*DECK DBNORM
      DOUBLE PRECISION FUNCTION DBNORM (N, A, NRA, ML, MU, W)
C-----------------------------------------------------------------------
C This function computes the norm of a banded N by N matrix,
C stored in the array A, that is consistent with the weighted max-norm
C on vectors, with weights stored in the array W.
C ML and MU are the lower and upper half-bandwidths of the matrix.
C NRA is the first dimension of the A array, NRA .ge. ML+MU+1.
C In terms of the matrix elements a(i,j), the norm is given by:
C   DBNORM = MAX(i=1,...,N) ( W(i) * Sum(j=1,...,N) ABS(a(i,j))/W(j) )
C-----------------------------------------------------------------------
      INTEGER N, NRA, ML, MU
      INTEGER I, I1, JLO, JHI, J
      DOUBLE PRECISION A, W
      DOUBLE PRECISION AN, SUM
      DIMENSION A(NRA,N), W(N)
      AN = 0.0D0
      DO 20 I = 1,N
        SUM = 0.0D0
        I1 = I + MU + 1
        JLO = MAX(I-ML,1)
        JHI = MIN(I+MU,N)
        DO 10 J = JLO,JHI
 10       SUM = SUM + ABS(A(I1-J,J))/W(J)
        AN = MAX(AN,SUM*W(I))
 20     CONTINUE
      DBNORM = AN
      RETURN
C----------------------- End of Function DBNORM ------------------------
      END
*DECK DSRCMA
      SUBROUTINE DSRCMA (RSAV, ISAV, JOB)
C-----------------------------------------------------------------------
C This routine saves or restores (depending on JOB) the contents of
C the Common blocks DLS001, DLSA01, which are used
C internally by one or more ODEPACK solvers.
C
C RSAV = real array of length 10 or more.
C ISAV = integer array of length 29 or more.
C JOB  = flag indicating to save or restore the Common blocks:
C        JOB  = 1 if Common is to be saved (written to RSAV/ISAV)
C        JOB  = 2 if Common is to be restored (read from RSAV/ISAV)
C        A call with JOB = 2 presumes a prior call with JOB = 1.
C-----------------------------------------------------------------------
      INTEGER ISAV, JOB
      INTEGER ILS, ILSA
      INTEGER I, LENRLS, LENILS, LENILA
      DOUBLE PRECISION RSAV
      DOUBLE PRECISION RLS, RLSA
      DIMENSION RSAV(*), ISAV(*)
      COMMON /DLS001/ RLS(9), ILS(25)
      COMMON /DLSA01/ RLSA, ILSA(4)
      DATA LENRLS/9/, LENILS/25/, LENILA/4/
C
      IF (JOB .EQ. 2) GO TO 100
      DO 10 I = 1,LENRLS
 10     RSAV(I) = RLS(I)
      RSAV(LENRLS+1) = RLSA
C
      DO 20 I = 1,LENILS
 20     ISAV(I) = ILS(I)
      DO 25 I = 1,LENILA
 25     ISAV(LENILS+I) = ILSA(I)
C
      RETURN
C
 100  CONTINUE
      DO 110 I = 1,LENRLS
 110     RLS(I) = RSAV(I)
      RLSA = RSAV(LENRLS+1)
C
      DO 120 I = 1,LENILS
 120     ILS(I) = ISAV(I)
      DO 125 I = 1,LENILA
 125     ILSA(I) = ISAV(LENILS+I)
C
      RETURN
C----------------------- End of Subroutine DSRCMA ----------------------
      END
*DECK DRCHEK
      SUBROUTINE DRCHEK (JOB, G, NEQ, Y, YH,NYH, G0, G1, GX, JROOT, IRT)
      EXTERNAL G
      INTEGER JOB, NEQ, NYH, JROOT, IRT
      DOUBLE PRECISION Y, YH, G0, G1, GX
      DIMENSION NEQ(*), Y(*), YH(NYH,*), G0(*), G1(*), GX(*), JROOT(*)
      INTEGER ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, 
     1   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     2   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      INTEGER IRFND, ITASKC, NGC, NGE
      INTEGER I, IFLAG, JFLAG
      DOUBLE PRECISION CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
      DOUBLE PRECISION T0, TLAST, TOUTC
      DOUBLE PRECISION HMING, T1, TEMP1, TEMP2, X
      LOGICAL ZROOT
      COMMON /DLS001/ CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,
     1   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, 
     2   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     3   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      COMMON /DLSR01/ T0, TLAST, TOUTC, IRFND, ITASKC, NGC, NGE
C-----------------------------------------------------------------------
C This routine checks for the presence of a root in the
C vicinity of the current T, in a manner depending on the
C input flag JOB.  It calls Subroutine DROOTS to locate the root
C as precisely as possible.
C
C In addition to variables described previously, DRCHEK
C uses the following for communication:
C JOB    = integer flag indicating type of call:
C          JOB = 1 means the problem is being initialized, and DRCHEK
C                  is to look for a root at or very near the initial T.
C          JOB = 2 means a continuation call to the solver was just
C                  made, and DRCHEK is to check for a root in the
C                  relevant part of the step last taken.
C          JOB = 3 means a successful step was just taken, and DRCHEK
C                  is to look for a root in the interval of the step.
C G0     = array of length NG, containing the value of g at T = T0.
C          G0 is input for JOB .ge. 2, and output in all cases.
C G1,GX  = arrays of length NG for work space.
C IRT    = completion flag:
C          IRT = 0  means no root was found.
C          IRT = -1 means JOB = 1 and a root was found too near to T.
C          IRT = 1  means a legitimate root was found (JOB = 2 or 3).
C                   On return, T0 is the root location, and Y is the
C                   corresponding solution vector.
C T0     = value of T at one endpoint of interval of interest.  Only
C          roots beyond T0 in the direction of integration are sought.
C          T0 is input if JOB .ge. 2, and output in all cases.
C          T0 is updated by DRCHEK, whether a root is found or not.
C TLAST  = last value of T returned by the solver (input only).
C TOUTC  = copy of TOUT (input only).
C IRFND  = input flag showing whether the last step taken had a root.
C          IRFND = 1 if it did, = 0 if not.
C ITASKC = copy of ITASK (input only).
C NGC    = copy of NG (input only).
C-----------------------------------------------------------------------
C
      IRT = 0
      DO 10 I = 1,NGC
 10     JROOT(I) = 0
      HMING = (ABS(TN) + ABS(H))*UROUND*100.0D0
C
      GO TO (100, 200, 300), JOB
C
C Evaluate g at initial T, and check for zero values. ------------------
 100  CONTINUE
      T0 = TN
      CALL G (NEQ, T0, Y, NGC, G0)
      NGE = 1
      ZROOT = .FALSE.
      DO 110 I = 1,NGC
 110    IF (ABS(G0(I)) .LE. 0.0D0) ZROOT = .TRUE.
      IF (.NOT. ZROOT) GO TO 190
C g has a zero at T.  Look at g at T + (small increment). --------------
      TEMP1 = SIGN(HMING,H)
      T0 = T0 + TEMP1
      TEMP2 = TEMP1/H
      DO 120 I = 1,N
 120    Y(I) = Y(I) + TEMP2*YH(I,2)
      CALL G (NEQ, T0, Y, NGC, G0)
      NGE = NGE + 1
      ZROOT = .FALSE.
      DO 130 I = 1,NGC
 130    IF (ABS(G0(I)) .LE. 0.0D0) ZROOT = .TRUE.
      IF (.NOT. ZROOT) GO TO 190
C g has a zero at T and also close to T.  Take error return. -----------
      IRT = -1
      RETURN
C
 190  CONTINUE
      RETURN
C
C
 200  CONTINUE
      IF (IRFND .EQ. 0) GO TO 260
C If a root was found on the previous step, evaluate G0 = g(T0). -------
      CALL DINTDY (T0, 0, YH, NYH, Y, IFLAG)
      CALL G (NEQ, T0, Y, NGC, G0)
      NGE = NGE + 1
      ZROOT = .FALSE.
      DO 210 I = 1,NGC
 210    IF (ABS(G0(I)) .LE. 0.0D0) ZROOT = .TRUE.
      IF (.NOT. ZROOT) GO TO 260
C g has a zero at T0.  Look at g at T + (small increment). -------------
      TEMP1 = SIGN(HMING,H)
      T0 = T0 + TEMP1
      IF ((T0 - TN)*H .LT. 0.0D0) GO TO 230
      TEMP2 = TEMP1/H
      DO 220 I = 1,N
 220    Y(I) = Y(I) + TEMP2*YH(I,2)
      GO TO 240
 230  CALL DINTDY (T0, 0, YH, NYH, Y, IFLAG)
 240  CALL G (NEQ, T0, Y, NGC, G0)
      NGE = NGE + 1
      ZROOT = .FALSE.
      DO 250 I = 1,NGC
        IF (ABS(G0(I)) .GT. 0.0D0) GO TO 250
        JROOT(I) = 1
        ZROOT = .TRUE.
 250    CONTINUE
      IF (.NOT. ZROOT) GO TO 260
C g has a zero at T0 and also close to T0.  Return root. ---------------
      IRT = 1
      RETURN
C G0 has no zero components.  Proceed to check relevant interval. ------
 260  IF (TN .EQ. TLAST) GO TO 390
C
 300  CONTINUE
C Set T1 to TN or TOUTC, whichever comes first, and get g at T1. -------
      IF (ITASKC.EQ.2 .OR. ITASKC.EQ.3 .OR. ITASKC.EQ.5) GO TO 310
      IF ((TOUTC - TN)*H .GE. 0.0D0) GO TO 310
      T1 = TOUTC
      IF ((T1 - T0)*H .LE. 0.0D0) GO TO 390
      CALL DINTDY (T1, 0, YH, NYH, Y, IFLAG)
      GO TO 330
 310  T1 = TN
      DO 320 I = 1,N
 320    Y(I) = YH(I,1)
 330  CALL G (NEQ, T1, Y, NGC, G1)
      NGE = NGE + 1
C Call DROOTS to search for root in interval from T0 to T1. ------------
      JFLAG = 0
 350  CONTINUE
      CALL DROOTS (NGC, HMING, JFLAG, T0, T1, G0, G1, GX, X, JROOT)
      IF (JFLAG .GT. 1) GO TO 360
      CALL DINTDY (X, 0, YH, NYH, Y, IFLAG)
      CALL G (NEQ, X, Y, NGC, GX)
      NGE = NGE + 1
      GO TO 350
 360  T0 = X
      CALL DCOPY (NGC, GX, 1, G0, 1)
      IF (JFLAG .EQ. 4) GO TO 390
C Found a root.  Interpolate to X and return. --------------------------
      CALL DINTDY (X, 0, YH, NYH, Y, IFLAG)
      IRT = 1
      RETURN
C
 390  CONTINUE
      RETURN
C----------------------- End of Subroutine DRCHEK ----------------------
      END
*DECK DROOTS
      SUBROUTINE DROOTS (NG, HMIN, JFLAG, X0, X1, G0, G1, GX, X, JROOT)
      INTEGER NG, JFLAG, JROOT
      DOUBLE PRECISION HMIN, X0, X1, G0, G1, GX, X
      DIMENSION G0(NG), G1(NG), GX(NG), JROOT(NG)
C-----------------------------------------------------------------------
C This subroutine finds the leftmost root of a set of arbitrary
C functions gi(x) (i = 1,...,NG) in an interval (X0,X1).  Only roots
C of odd multiplicity (i.e. changes of sign of the gi) are found.
C Here the sign of X1 - X0 is arbitrary, but is constant for a given
C problem, and -leftmost- means nearest to X0.
C The values of the vector-valued function g(x) = (gi, i=1...NG)
C are communicated through the call sequence of DROOTS.
C The method used is the Illinois algorithm.
C
C Reference:
C Kathie L. Hiebert and Lawrence F. Shampine, Implicitly Defined
C Output Points for Solutions of ODEs, Sandia Report SAND80-0180,
C February 1980.
C
C Description of parameters.
C
C NG     = number of functions gi, or the number of components of
C          the vector valued function g(x).  Input only.
C
C HMIN   = resolution parameter in X.  Input only.  When a root is
C          found, it is located only to within an error of HMIN in X.
C          Typically, HMIN should be set to something on the order of
C               100 * UROUND * MAX(ABS(X0),ABS(X1)),
C          where UROUND is the unit roundoff of the machine.
C
C JFLAG  = integer flag for input and output communication.
C
C          On input, set JFLAG = 0 on the first call for the problem,
C          and leave it unchanged until the problem is completed.
C          (The problem is completed when JFLAG .ge. 2 on return.)
C
C          On output, JFLAG has the following values and meanings:
C          JFLAG = 1 means DROOTS needs a value of g(x).  Set GX = g(X)
C                    and call DROOTS again.
C          JFLAG = 2 means a root has been found.  The root is
C                    at X, and GX contains g(X).  (Actually, X is the
C                    rightmost approximation to the root on an interval
C                    (X0,X1) of size HMIN or less.)
C          JFLAG = 3 means X = X1 is a root, with one or more of the gi
C                    being zero at X1 and no sign changes in (X0,X1).
C                    GX contains g(X) on output.
C          JFLAG = 4 means no roots (of odd multiplicity) were
C                    found in (X0,X1) (no sign changes).
C
C X0,X1  = endpoints of the interval where roots are sought.
C          X1 and X0 are input when JFLAG = 0 (first call), and
C          must be left unchanged between calls until the problem is
C          completed.  X0 and X1 must be distinct, but X1 - X0 may be
C          of either sign.  However, the notion of -left- and -right-
C          will be used to mean nearer to X0 or X1, respectively.
C          When JFLAG .ge. 2 on return, X0 and X1 are output, and
C          are the endpoints of the relevant interval.
C
C G0,G1  = arrays of length NG containing the vectors g(X0) and g(X1),
C          respectively.  When JFLAG = 0, G0 and G1 are input and
C          none of the G0(i) should be zero.
C          When JFLAG .ge. 2 on return, G0 and G1 are output.
C
C GX     = array of length NG containing g(X).  GX is input
C          when JFLAG = 1, and output when JFLAG .ge. 2.
C
C X      = independent variable value.  Output only.
C          When JFLAG = 1 on output, X is the point at which g(x)
C          is to be evaluated and loaded into GX.
C          When JFLAG = 2 or 3, X is the root.
C          When JFLAG = 4, X is the right endpoint of the interval, X1.
C
C JROOT  = integer array of length NG.  Output only.
C          When JFLAG = 2 or 3, JROOT indicates which components
C          of g(x) have a root at X.  JROOT(i) is 1 if the i-th
C          component has a root, and JROOT(i) = 0 otherwise.
C-----------------------------------------------------------------------
      INTEGER I, IMAX, IMXOLD, LAST, NXLAST
      DOUBLE PRECISION ALPHA, T2, TMAX, X2, ZERO
      LOGICAL ZROOT, SGNCHG, XROOT
      SAVE ALPHA, X2, IMAX, LAST
      DATA ZERO/0.0D0/
C
      IF (JFLAG .EQ. 1) GO TO 200
C JFLAG .ne. 1.  Check for change in sign of g or zero at X1. ----------
      IMAX = 0
      TMAX = ZERO
      ZROOT = .FALSE.
      DO 120 I = 1,NG
        IF (ABS(G1(I)) .GT. ZERO) GO TO 110
        ZROOT = .TRUE.
        GO TO 120
C At this point, G0(i) has been checked and cannot be zero. ------------
 110    IF (SIGN(1.0D0,G0(I)) .EQ. SIGN(1.0D0,G1(I))) GO TO 120
          T2 = ABS(G1(I)/(G1(I)-G0(I)))
          IF (T2 .LE. TMAX) GO TO 120
            TMAX = T2
            IMAX = I
 120    CONTINUE
      IF (IMAX .GT. 0) GO TO 130
      SGNCHG = .FALSE.
      GO TO 140
 130  SGNCHG = .TRUE.
 140  IF (.NOT. SGNCHG) GO TO 400
C There is a sign change.  Find the first root in the interval. --------
      XROOT = .FALSE.
      NXLAST = 0
      LAST = 1
C
C Repeat until the first root in the interval is found.  Loop point. ---
 150  CONTINUE
      IF (XROOT) GO TO 300
      IF (NXLAST .EQ. LAST) GO TO 160
      ALPHA = 1.0D0
      GO TO 180
 160  IF (LAST .EQ. 0) GO TO 170
      ALPHA = 0.5D0*ALPHA
      GO TO 180
 170  ALPHA = 2.0D0*ALPHA
 180  X2 = X1 - (X1-X0)*G1(IMAX)/(G1(IMAX) - ALPHA*G0(IMAX))
      IF ((ABS(X2-X0) .LT. HMIN) .AND.
     1   (ABS(X1-X0) .GT. 10.0D0*HMIN)) X2 = X0 + 0.1D0*(X1-X0)
      JFLAG = 1
      X = X2
C Return to the calling routine to get a value of GX = g(X). -----------
      RETURN
C Check to see in which interval g changes sign. -----------------------
 200  IMXOLD = IMAX
      IMAX = 0
      TMAX = ZERO
      ZROOT = .FALSE.
      DO 220 I = 1,NG
        IF (ABS(GX(I)) .GT. ZERO) GO TO 210
        ZROOT = .TRUE.
        GO TO 220
C Neither G0(i) nor GX(i) can be zero at this point. -------------------
 210    IF (SIGN(1.0D0,G0(I)) .EQ. SIGN(1.0D0,GX(I))) GO TO 220
          T2 = ABS(GX(I)/(GX(I) - G0(I)))
          IF (T2 .LE. TMAX) GO TO 220
            TMAX = T2
            IMAX = I
 220    CONTINUE
      IF (IMAX .GT. 0) GO TO 230
      SGNCHG = .FALSE.
      IMAX = IMXOLD
      GO TO 240
 230  SGNCHG = .TRUE.
 240  NXLAST = LAST
      IF (.NOT. SGNCHG) GO TO 250
C Sign change between X0 and X2, so replace X1 with X2. ----------------
      X1 = X2
      CALL DCOPY (NG, GX, 1, G1, 1)
      LAST = 1
      XROOT = .FALSE.
      GO TO 270
 250  IF (.NOT. ZROOT) GO TO 260
C Zero value at X2 and no sign change in (X0,X2), so X2 is a root. -----
      X1 = X2
      CALL DCOPY (NG, GX, 1, G1, 1)
      XROOT = .TRUE.
      GO TO 270
C No sign change between X0 and X2.  Replace X0 with X2. ---------------
 260  CONTINUE
      CALL DCOPY (NG, GX, 1, G0, 1)
      X0 = X2
      LAST = 0
      XROOT = .FALSE.
 270  IF (ABS(X1-X0) .LE. HMIN) XROOT = .TRUE.
      GO TO 150
C
C Return with X1 as the root.  Set JROOT.  Set X = X1 and GX = G1. -----
 300  JFLAG = 2
      X = X1
      CALL DCOPY (NG, G1, 1, GX, 1)
      DO 320 I = 1,NG
        JROOT(I) = 0
        IF (ABS(G1(I)) .GT. ZERO) GO TO 310
          JROOT(I) = 1
          GO TO 320
 310    IF (SIGN(1.0D0,G0(I)) .NE. SIGN(1.0D0,G1(I))) JROOT(I) = 1
 320    CONTINUE
      RETURN
C
C No sign change in the interval.  Check for zero at right endpoint. ---
 400  IF (.NOT. ZROOT) GO TO 420
C
C Zero value at X1 and no sign change in (X0,X1).  Return JFLAG = 3. ---
      X = X1
      CALL DCOPY (NG, G1, 1, GX, 1)
      DO 410 I = 1,NG
        JROOT(I) = 0
        IF (ABS(G1(I)) .LE. ZERO) JROOT (I) = 1
 410  CONTINUE
      JFLAG = 3
      RETURN
C
C No sign changes in this interval.  Set X = X1, return JFLAG = 4. -----
 420  CALL DCOPY (NG, G1, 1, GX, 1)
      X = X1
      JFLAG = 4
      RETURN
C----------------------- End of Subroutine DROOTS ----------------------
      END
*DECK DSRCAR
      SUBROUTINE DSRCAR (RSAV, ISAV, JOB)
C-----------------------------------------------------------------------
C This routine saves or restores (depending on JOB) the contents of
C the Common blocks DLS001, DLSA01, DLSR01, which are used
C internally by one or more ODEPACK solvers.
C
C RSAV = real array of length 13 or more.
C ISAV = integer array of length 33 or more.
C JOB  = flag indicating to save or restore the Common blocks:
C        JOB  = 1 if Common is to be saved (written to RSAV/ISAV)
C        JOB  = 2 if Common is to be restored (read from RSAV/ISAV)
C        A call with JOB = 2 presumes a prior call with JOB = 1.
C-----------------------------------------------------------------------
      INTEGER ISAV, JOB
      INTEGER ILS, ILSA, ILSR
      INTEGER I, IOFF, LENRLS, LENILS, LENILA, LENRLR, LENILR
      DOUBLE PRECISION RSAV
      DOUBLE PRECISION RLS, RLSA, RLSR
      DIMENSION RSAV(*), ISAV(*)
      COMMON /DLS001/ RLS(9), ILS(25)
      COMMON /DLSA01/ RLSA, ILSA(4)
      COMMON /DLSR01/ RLSR(3), ILSR(4)
      DATA LENRLS/9/, LENILS/25/, LENILA/4/
      DATA LENRLR/3/, LENILR/4/
C
      IF (JOB .EQ. 2) GO TO 100
      DO 10 I = 1,LENRLS
 10     RSAV(I) = RLS(I)
      RSAV(LENRLS+1) = RLSA
      IOFF = LENRLS + 1
      DO 20 I = 1,LENRLR
 20     RSAV(IOFF+I) = RLSR(I)
C
      DO 30 I = 1,LENILS
 30     ISAV(I) = ILS(I)
      DO 35 I = 1,LENILA
 35     ISAV(LENILS+I) = ILSA(I)
      IOFF = LENILS + LENILA
      DO 40 I = 1,LENILR
 40     ISAV(IOFF+I) = ILSR(I)
C
      RETURN
C
 100  CONTINUE
      DO 110 I = 1,LENRLS
 110     RLS(I) = RSAV(I)
      RLSA = RSAV(LENRLS+1)
      IOFF = LENRLS + 1
      DO 120 I = 1,LENRLR
 120     RLSR(I) = RSAV(IOFF+I)
C
      DO 130 I = 1,LENILS
 130     ILS(I) = ISAV(I)
      DO 135 I = 1,LENILA
 135     ILSA(I) = ISAV(LENILS+I)
      IOFF = LENILS + LENILA
      DO 140 I = 1,LENILR
 140     ILSR(I) = ISAV(IOFF+I)
C
      RETURN
C----------------------- End of Subroutine DSRCAR ----------------------
      END
*DECK DSTODPK
      SUBROUTINE DSTODPK (NEQ, Y, YH, NYH, YH1, EWT, SAVF, SAVX, ACOR,
     1   WM, IWM, F, JAC, PSOL)
      EXTERNAL F, JAC, PSOL
      INTEGER NEQ, NYH, IWM
      DOUBLE PRECISION Y, YH, YH1, EWT, SAVF, SAVX, ACOR, WM
      DIMENSION NEQ(*), Y(*), YH(NYH,*), YH1(*), EWT(*), SAVF(*),
     1   SAVX(*), ACOR(*), WM(*), IWM(*)
      INTEGER ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, 
     1   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     2   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NPE, NQU
      INTEGER JPRE, JACFLG, LOCWP, LOCIWP, LSAVX, KMP, MAXL, MNEWT,
     1   NNI, NLI, NPS, NCFN, NCFL
      DOUBLE PRECISION CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
      DOUBLE PRECISION DELT, EPCON, SQRTN, RSQRTN
      COMMON /DLS001/ CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,
     1   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, 
     2   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     3   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NPE, NQU
      COMMON /DLPK01/ DELT, EPCON, SQRTN, RSQRTN,
     1   JPRE, JACFLG, LOCWP, LOCIWP, LSAVX, KMP, MAXL, MNEWT,
     2   NNI, NLI, NPS, NCFN, NCFL
C-----------------------------------------------------------------------
C DSTODPK performs one step of the integration of an initial value
C problem for a system of Ordinary Differential Equations.
C-----------------------------------------------------------------------
C The following changes were made to generate Subroutine DSTODPK
C from Subroutine DSTODE:
C 1. The array SAVX was added to the call sequence.
C 2. PJAC and SLVS were replaced by PSOL in the call sequence.
C 3. The Common block /DLPK01/ was added for communication.
C 4. The test constant EPCON is loaded into Common below statement
C    numbers 125 and 155, and used below statement 400.
C 5. The Newton iteration counter MNEWT is set below 220 and 400.
C 6. The call to PJAC was replaced with a call to DPKSET (fixed name),
C    with a longer call sequence, called depending on JACFLG.
C 7. The corrector residual is stored in SAVX (not Y) at 360,
C    and the solution vector is in SAVX in the 380 loop.
C 8. SLVS was renamed DSOLPK and includes NEQ, SAVX, EWT, F, and JAC.
C    SAVX was added because DSOLPK now needs Y and SAVF undisturbed.
C 9. The nonlinear convergence failure count NCFN is set at 430.
C-----------------------------------------------------------------------
C Note: DSTODPK is independent of the value of the iteration method
C indicator MITER, when this is .ne. 0, and hence is independent
C of the type of chord method used, or the Jacobian structure.
C Communication with DSTODPK is done with the following variables:
C
C NEQ    = integer array containing problem size in NEQ(1), and
C          passed as the NEQ argument in all calls to F and JAC.
C Y      = an array of length .ge. N used as the Y argument in
C          all calls to F and JAC.
C YH     = an NYH by LMAX array containing the dependent variables
C          and their approximate scaled derivatives, where
C          LMAX = MAXORD + 1.  YH(i,j+1) contains the approximate
C          j-th derivative of y(i), scaled by H**j/factorial(j)
C          (j = 0,1,...,NQ).  On entry for the first step, the first
C          two columns of YH must be set from the initial values.
C NYH    = a constant integer .ge. N, the first dimension of YH.
C YH1    = a one-dimensional array occupying the same space as YH.
C EWT    = an array of length N containing multiplicative weights
C          for local error measurements.  Local errors in y(i) are
C          compared to 1.0/EWT(i) in various error tests.
C SAVF   = an array of working storage, of length N.
C          Also used for input of YH(*,MAXORD+2) when JSTART = -1
C          and MAXORD .lt. the current order NQ.
C SAVX   = an array of working storage, of length N.
C ACOR   = a work array of length N, used for the accumulated
C          corrections.  On a successful return, ACOR(i) contains
C          the estimated one-step local error in y(i).
C WM,IWM = real and integer work arrays associated with matrix
C          operations in chord iteration (MITER .ne. 0).
C CCMAX  = maximum relative change in H*EL0 before DPKSET is called.
C H      = the step size to be attempted on the next step.
C          H is altered by the error control algorithm during the
C          problem.  H can be either positive or negative, but its
C          sign must remain constant throughout the problem.
C HMIN   = the minimum absolute value of the step size H to be used.
C HMXI   = inverse of the maximum absolute value of H to be used.
C          HMXI = 0.0 is allowed and corresponds to an infinite HMAX.
C          HMIN and HMXI may be changed at any time, but will not
C          take effect until the next change of H is considered.
C TN     = the independent variable. TN is updated on each step taken.
C JSTART = an integer used for input only, with the following
C          values and meanings:
C               0  perform the first step.
C           .gt.0  take a new step continuing from the last.
C              -1  take the next step with a new value of H, MAXORD,
C                    N, METH, MITER, and/or matrix parameters.
C              -2  take the next step with a new value of H,
C                    but with other inputs unchanged.
C          On return, JSTART is set to 1 to facilitate continuation.
C KFLAG  = a completion code with the following meanings:
C               0  the step was succesful.
C              -1  the requested error could not be achieved.
C              -2  corrector convergence could not be achieved.
C              -3  fatal error in DPKSET or DSOLPK.
C          A return with KFLAG = -1 or -2 means either
C          ABS(H) = HMIN or 10 consecutive failures occurred.
C          On a return with KFLAG negative, the values of TN and
C          the YH array are as of the beginning of the last
C          step, and H is the last step size attempted.
C MAXORD = the maximum order of integration method to be allowed.
C MAXCOR = the maximum number of corrector iterations allowed.
C MSBP   = maximum number of steps between DPKSET calls (MITER .gt. 0).
C MXNCF  = maximum number of convergence failures allowed.
C METH/MITER = the method flags.  See description in driver.
C N      = the number of first-order differential equations.
C-----------------------------------------------------------------------
      INTEGER I, I1, IREDO, IRET, J, JB, M, NCF, NEWQ
      INTEGER IALTH, IPUP, LMAX, MEO, NQNYH, NSLP
      DOUBLE PRECISION DCON, DDN, DEL, DELP, DSM, DUP, EXDN, EXSM, EXUP,
     1   R, RH, RHDN, RHSM, RHUP, TOLD, DVNORM
      DOUBLE PRECISION CONIT, CRATE, EL(13), ELCO(13,12), HOLD, RMAX,    
     1   TESCO(3,12)
      SAVE CONIT, CRATE, EL, ELCO, HOLD, RMAX, TESCO,
     1   IALTH, IPUP, LMAX, MEO, NQNYH, NSLP
C
      KFLAG = 0
      TOLD = TN
      NCF = 0
      IERPJ = 0
      IERSL = 0
      JCUR = 0
      ICF = 0
      DELP = 0.0D0
      IF (JSTART .GT. 0) GO TO 200
      IF (JSTART .EQ. -1) GO TO 100
      IF (JSTART .EQ. -2) GO TO 160
C-----------------------------------------------------------------------
C On the first call, the order is set to 1, and other variables are
C initialized.  RMAX is the maximum ratio by which H can be increased
C in a single step.  It is initially 1.E4 to compensate for the small
C initial H, but then is normally equal to 10.  If a failure
C occurs (in corrector convergence or error test), RMAX is set at 2
C for the next increase.
C-----------------------------------------------------------------------
      LMAX = MAXORD + 1
      NQ = 1
      L = 2
      IALTH = 2
      RMAX = 10000.0D0
      RC = 0.0D0
      EL0 = 1.0D0
      CRATE = 0.7D0
      HOLD = H
      MEO = METH
      NSLP = 0
      IPUP = MITER
      IRET = 3
      GO TO 140
C-----------------------------------------------------------------------
C The following block handles preliminaries needed when JSTART = -1.
C IPUP is set to MITER to force a matrix update.
C If an order increase is about to be considered (IALTH = 1),
C IALTH is reset to 2 to postpone consideration one more step.
C If the caller has changed METH, DCFODE is called to reset
C the coefficients of the method.
C If the caller has changed MAXORD to a value less than the current
C order NQ, NQ is reduced to MAXORD, and a new H chosen accordingly.
C If H is to be changed, YH must be rescaled.
C If H or METH is being changed, IALTH is reset to L = NQ + 1
C to prevent further changes in H for that many steps.
C-----------------------------------------------------------------------
 100  IPUP = MITER
      LMAX = MAXORD + 1
      IF (IALTH .EQ. 1) IALTH = 2
      IF (METH .EQ. MEO) GO TO 110
      CALL DCFODE (METH, ELCO, TESCO)
      MEO = METH
      IF (NQ .GT. MAXORD) GO TO 120
      IALTH = L
      IRET = 1
      GO TO 150
 110  IF (NQ .LE. MAXORD) GO TO 160
 120  NQ = MAXORD
      L = LMAX
      DO 125 I = 1,L
 125    EL(I) = ELCO(I,NQ)
      NQNYH = NQ*NYH
      RC = RC*EL(1)/EL0
      EL0 = EL(1)
      CONIT = 0.5D0/(NQ+2)
      EPCON = CONIT*TESCO(2,NQ)
      DDN = DVNORM (N, SAVF, EWT)/TESCO(1,L)
      EXDN = 1.0D0/L
      RHDN = 1.0D0/(1.3D0*DDN**EXDN + 0.0000013D0)
      RH = MIN(RHDN,1.0D0)
      IREDO = 3
      IF (H .EQ. HOLD) GO TO 170
      RH = MIN(RH,ABS(H/HOLD))
      H = HOLD
      GO TO 175
C-----------------------------------------------------------------------
C DCFODE is called to get all the integration coefficients for the
C current METH.  Then the EL vector and related constants are reset
C whenever the order NQ is changed, or at the start of the problem.
C-----------------------------------------------------------------------
 140  CALL DCFODE (METH, ELCO, TESCO)
 150  DO 155 I = 1,L
 155    EL(I) = ELCO(I,NQ)
      NQNYH = NQ*NYH
      RC = RC*EL(1)/EL0
      EL0 = EL(1)
      CONIT = 0.5D0/(NQ+2)
      EPCON = CONIT*TESCO(2,NQ)
      GO TO (160, 170, 200), IRET
C-----------------------------------------------------------------------
C If H is being changed, the H ratio RH is checked against
C RMAX, HMIN, and HMXI, and the YH array rescaled.  IALTH is set to
C L = NQ + 1 to prevent a change of H for that many steps, unless
C forced by a convergence or error test failure.
C-----------------------------------------------------------------------
 160  IF (H .EQ. HOLD) GO TO 200
      RH = H/HOLD
      H = HOLD
      IREDO = 3
      GO TO 175
 170  RH = MAX(RH,HMIN/ABS(H))
 175  RH = MIN(RH,RMAX)
      RH = RH/MAX(1.0D0,ABS(H)*HMXI*RH)
      R = 1.0D0
      DO 180 J = 2,L
        R = R*RH
        DO 180 I = 1,N
 180      YH(I,J) = YH(I,J)*R
      H = H*RH
      RC = RC*RH
      IALTH = L
      IF (IREDO .EQ. 0) GO TO 690
C-----------------------------------------------------------------------
C This section computes the predicted values by effectively
C multiplying the YH array by the Pascal triangle matrix.
C The flag IPUP is set according to whether matrix data is involved
C (JACFLG .ne. 0) or not (JACFLG = 0), to trigger a call to DPKSET.
C IPUP is set to MITER when RC differs from 1 by more than CCMAX,
C and at least every MSBP steps, when JACFLG = 1.
C RC is the ratio of new to old values of the coefficient  H*EL(1).
C-----------------------------------------------------------------------
 200  IF (JACFLG .NE. 0) GO TO 202
      IPUP = 0
      CRATE = 0.7D0
      GO TO 205
 202  IF (ABS(RC-1.0D0) .GT. CCMAX) IPUP = MITER
      IF (NST .GE. NSLP+MSBP) IPUP = MITER
 205  TN = TN + H
      I1 = NQNYH + 1
      DO 215 JB = 1,NQ
        I1 = I1 - NYH
CDIR$ IVDEP
        DO 210 I = I1,NQNYH
 210      YH1(I) = YH1(I) + YH1(I+NYH)
 215    CONTINUE
C-----------------------------------------------------------------------
C Up to MAXCOR corrector iterations are taken.  A convergence test is
C made on the RMS-norm of each correction, weighted by the error
C weight vector EWT.  The sum of the corrections is accumulated in the
C vector ACOR(i).  The YH array is not altered in the corrector loop.
C-----------------------------------------------------------------------
 220  M = 0
      MNEWT = 0
      DO 230 I = 1,N
 230    Y(I) = YH(I,1)
      CALL F (NEQ, TN, Y, SAVF)
      NFE = NFE + 1
      IF (IPUP .LE. 0) GO TO 250
C-----------------------------------------------------------------------
C If indicated, DPKSET is called to update any matrix data needed,
C before starting the corrector iteration.
C IPUP is set to 0 as an indicator that this has been done.
C-----------------------------------------------------------------------
      CALL DPKSET (NEQ, Y, YH1, EWT, ACOR, SAVF, WM, IWM, F, JAC)
      IPUP = 0
      RC = 1.0D0
      NSLP = NST
      CRATE = 0.7D0
      IF (IERPJ .NE. 0) GO TO 430
 250  DO 260 I = 1,N
 260    ACOR(I) = 0.0D0
 270  IF (MITER .NE. 0) GO TO 350
C-----------------------------------------------------------------------
C In the case of functional iteration, update Y directly from
C the result of the last function evaluation.
C-----------------------------------------------------------------------
      DO 290 I = 1,N
        SAVF(I) = H*SAVF(I) - YH(I,2)
 290    Y(I) = SAVF(I) - ACOR(I)
      DEL = DVNORM (N, Y, EWT)
      DO 300 I = 1,N
        Y(I) = YH(I,1) + EL(1)*SAVF(I)
 300    ACOR(I) = SAVF(I)
      GO TO 400
C-----------------------------------------------------------------------
C In the case of the chord method, compute the corrector error,
C and solve the linear system with that as right-hand side and
C P as coefficient matrix.
C-----------------------------------------------------------------------
 350  DO 360 I = 1,N
 360    SAVX(I) = H*SAVF(I) - (YH(I,2) + ACOR(I))
      CALL DSOLPK (NEQ, Y, SAVF, SAVX, EWT, WM, IWM, F, PSOL)
      IF (IERSL .LT. 0) GO TO 430
      IF (IERSL .GT. 0) GO TO 410
      DEL = DVNORM (N, SAVX, EWT)
      DO 380 I = 1,N
        ACOR(I) = ACOR(I) + SAVX(I)
 380    Y(I) = YH(I,1) + EL(1)*ACOR(I)
C-----------------------------------------------------------------------
C Test for convergence.  If M .gt. 0, an estimate of the convergence
C rate constant is stored in CRATE, and this is used in the test.
C-----------------------------------------------------------------------
 400  IF (M .NE. 0) CRATE = MAX(0.2D0*CRATE,DEL/DELP)
      DCON = DEL*MIN(1.0D0,1.5D0*CRATE)/EPCON
      IF (DCON .LE. 1.0D0) GO TO 450
      M = M + 1
      IF (M .EQ. MAXCOR) GO TO 410
      IF (M .GE. 2 .AND. DEL .GT. 2.0D0*DELP) GO TO 410
      MNEWT = M
      DELP = DEL
      CALL F (NEQ, TN, Y, SAVF)
      NFE = NFE + 1
      GO TO 270
C-----------------------------------------------------------------------
C The corrector iteration failed to converge.
C If MITER .ne. 0 and the Jacobian is out of date, DPKSET is called for
C the next try.  Otherwise the YH array is retracted to its values
C before prediction, and H is reduced, if possible.  If H cannot be
C reduced or MXNCF failures have occurred, exit with KFLAG = -2.
C-----------------------------------------------------------------------
 410  IF (MITER.EQ.0 .OR. JCUR.EQ.1 .OR. JACFLG.EQ.0) GO TO 430
      ICF = 1
      IPUP = MITER
      GO TO 220
 430  ICF = 2
      NCF = NCF + 1
      NCFN = NCFN + 1
      RMAX = 2.0D0
      TN = TOLD
      I1 = NQNYH + 1
      DO 445 JB = 1,NQ
        I1 = I1 - NYH
CDIR$ IVDEP
        DO 440 I = I1,NQNYH
 440      YH1(I) = YH1(I) - YH1(I+NYH)
 445    CONTINUE
      IF (IERPJ .LT. 0 .OR. IERSL .LT. 0) GO TO 680
      IF (ABS(H) .LE. HMIN*1.00001D0) GO TO 670
      IF (NCF .EQ. MXNCF) GO TO 670
      RH = 0.5D0
      IPUP = MITER
      IREDO = 1
      GO TO 170
C-----------------------------------------------------------------------
C The corrector has converged.  JCUR is set to 0
C to signal that the Jacobian involved may need updating later.
C The local error test is made and control passes to statement 500
C if it fails.
C-----------------------------------------------------------------------
 450  JCUR = 0
      IF (M .EQ. 0) DSM = DEL/TESCO(2,NQ)
      IF (M .GT. 0) DSM = DVNORM (N, ACOR, EWT)/TESCO(2,NQ)
      IF (DSM .GT. 1.0D0) GO TO 500
C-----------------------------------------------------------------------
C After a successful step, update the YH array.
C Consider changing H if IALTH = 1.  Otherwise decrease IALTH by 1.
C If IALTH is then 1 and NQ .lt. MAXORD, then ACOR is saved for
C use in a possible order increase on the next step.
C If a change in H is considered, an increase or decrease in order
C by one is considered also.  A change in H is made only if it is by a
C factor of at least 1.1.  If not, IALTH is set to 3 to prevent
C testing for that many steps.
C-----------------------------------------------------------------------
      KFLAG = 0
      IREDO = 0
      NST = NST + 1
      HU = H
      NQU = NQ
      DO 470 J = 1,L
        DO 470 I = 1,N
 470      YH(I,J) = YH(I,J) + EL(J)*ACOR(I)
      IALTH = IALTH - 1
      IF (IALTH .EQ. 0) GO TO 520
      IF (IALTH .GT. 1) GO TO 700
      IF (L .EQ. LMAX) GO TO 700
      DO 490 I = 1,N
 490    YH(I,LMAX) = ACOR(I)
      GO TO 700
C-----------------------------------------------------------------------
C The error test failed.  KFLAG keeps track of multiple failures.
C Restore TN and the YH array to their previous values, and prepare
C to try the step again.  Compute the optimum step size for this or
C one lower order.  After 2 or more failures, H is forced to decrease
C by a factor of 0.2 or less.
C-----------------------------------------------------------------------
 500  KFLAG = KFLAG - 1
      TN = TOLD
      I1 = NQNYH + 1
      DO 515 JB = 1,NQ
        I1 = I1 - NYH
CDIR$ IVDEP
        DO 510 I = I1,NQNYH
 510      YH1(I) = YH1(I) - YH1(I+NYH)
 515    CONTINUE
      RMAX = 2.0D0
      IF (ABS(H) .LE. HMIN*1.00001D0) GO TO 660
      IF (KFLAG .LE. -3) GO TO 640
      IREDO = 2
      RHUP = 0.0D0
      GO TO 540
C-----------------------------------------------------------------------
C Regardless of the success or failure of the step, factors
C RHDN, RHSM, and RHUP are computed, by which H could be multiplied
C at order NQ - 1, order NQ, or order NQ + 1, respectively.
C In the case of failure, RHUP = 0.0 to avoid an order increase.
C the largest of these is determined and the new order chosen
C accordingly.  If the order is to be increased, we compute one
C additional scaled derivative.
C-----------------------------------------------------------------------
 520  RHUP = 0.0D0
      IF (L .EQ. LMAX) GO TO 540
      DO 530 I = 1,N
 530    SAVF(I) = ACOR(I) - YH(I,LMAX)
      DUP = DVNORM (N, SAVF, EWT)/TESCO(3,NQ)
      EXUP = 1.0D0/(L+1)
      RHUP = 1.0D0/(1.4D0*DUP**EXUP + 0.0000014D0)
 540  EXSM = 1.0D0/L
      RHSM = 1.0D0/(1.2D0*DSM**EXSM + 0.0000012D0)
      RHDN = 0.0D0
      IF (NQ .EQ. 1) GO TO 560
      DDN = DVNORM (N, YH(1,L), EWT)/TESCO(1,NQ)
      EXDN = 1.0D0/NQ
      RHDN = 1.0D0/(1.3D0*DDN**EXDN + 0.0000013D0)
 560  IF (RHSM .GE. RHUP) GO TO 570
      IF (RHUP .GT. RHDN) GO TO 590
      GO TO 580
 570  IF (RHSM .LT. RHDN) GO TO 580
      NEWQ = NQ
      RH = RHSM
      GO TO 620
 580  NEWQ = NQ - 1
      RH = RHDN
      IF (KFLAG .LT. 0 .AND. RH .GT. 1.0D0) RH = 1.0D0
      GO TO 620
 590  NEWQ = L
      RH = RHUP
      IF (RH .LT. 1.1D0) GO TO 610
      R = EL(L)/L
      DO 600 I = 1,N
 600    YH(I,NEWQ+1) = ACOR(I)*R
      GO TO 630
 610  IALTH = 3
      GO TO 700
 620  IF ((KFLAG .EQ. 0) .AND. (RH .LT. 1.1D0)) GO TO 610
      IF (KFLAG .LE. -2) RH = MIN(RH,0.2D0)
C-----------------------------------------------------------------------
C If there is a change of order, reset NQ, L, and the coefficients.
C In any case H is reset according to RH and the YH array is rescaled.
C Then exit from 690 if the step was OK, or redo the step otherwise.
C-----------------------------------------------------------------------
      IF (NEWQ .EQ. NQ) GO TO 170
 630  NQ = NEWQ
      L = NQ + 1
      IRET = 2
      GO TO 150
C-----------------------------------------------------------------------
C Control reaches this section if 3 or more failures have occured.
C If 10 failures have occurred, exit with KFLAG = -1.
C It is assumed that the derivatives that have accumulated in the
C YH array have errors of the wrong order.  Hence the first
C derivative is recomputed, and the order is set to 1.  Then
C H is reduced by a factor of 10, and the step is retried,
C until it succeeds or H reaches HMIN.
C-----------------------------------------------------------------------
 640  IF (KFLAG .EQ. -10) GO TO 660
      RH = 0.1D0
      RH = MAX(HMIN/ABS(H),RH)
      H = H*RH
      DO 645 I = 1,N
 645    Y(I) = YH(I,1)
      CALL F (NEQ, TN, Y, SAVF)
      NFE = NFE + 1
      DO 650 I = 1,N
 650    YH(I,2) = H*SAVF(I)
      IPUP = MITER
      IALTH = 5
      IF (NQ .EQ. 1) GO TO 200
      NQ = 1
      L = 2
      IRET = 3
      GO TO 150
C-----------------------------------------------------------------------
C All returns are made through this section.  H is saved in HOLD
C to allow the caller to change H on the next step.
C-----------------------------------------------------------------------
 660  KFLAG = -1
      GO TO 720
 670  KFLAG = -2
      GO TO 720
 680  KFLAG = -3
      GO TO 720
 690  RMAX = 10.0D0
 700  R = 1.0D0/TESCO(2,NQU)
      DO 710 I = 1,N
 710    ACOR(I) = ACOR(I)*R
 720  HOLD = H
      JSTART = 1
      RETURN
C----------------------- End of Subroutine DSTODPK ---------------------
      END
*DECK DPKSET
      SUBROUTINE DPKSET (NEQ, Y, YSV, EWT, FTEM, SAVF, WM, IWM, F, JAC)
      EXTERNAL F, JAC
      INTEGER NEQ, IWM
      DOUBLE PRECISION Y, YSV, EWT, FTEM, SAVF, WM
      DIMENSION NEQ(*), Y(*), YSV(*), EWT(*), FTEM(*), SAVF(*),
     1   WM(*), IWM(*)
      INTEGER ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, 
     1   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     2   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NPE, NQU
      INTEGER JPRE, JACFLG, LOCWP, LOCIWP, LSAVX, KMP, MAXL, MNEWT,
     1   NNI, NLI, NPS, NCFN, NCFL
      DOUBLE PRECISION CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
      DOUBLE PRECISION DELT, EPCON, SQRTN, RSQRTN
      COMMON /DLS001/ CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,
     1   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, 
     2   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     3   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NPE, NQU
      COMMON /DLPK01/ DELT, EPCON, SQRTN, RSQRTN,
     1   JPRE, JACFLG, LOCWP, LOCIWP, LSAVX, KMP, MAXL, MNEWT,
     2   NNI, NLI, NPS, NCFN, NCFL
C-----------------------------------------------------------------------
C DPKSET is called by DSTODPK to interface with the user-supplied
C routine JAC, to compute and process relevant parts of
C the matrix P = I - H*EL(1)*J , where J is the Jacobian df/dy,
C as need for preconditioning matrix operations later.
C
C In addition to variables described previously, communication
C with DPKSET uses the following:
C Y     = array containing predicted values on entry.
C YSV   = array containing predicted y, to be saved (YH1 in DSTODPK).
C FTEM  = work array of length N (ACOR in DSTODPK).
C SAVF  = array containing f evaluated at predicted y.
C WM    = real work space for matrices.
C         Space for preconditioning data starts at WM(LOCWP).
C IWM   = integer work space.
C         Space for preconditioning data starts at IWM(LOCIWP).
C IERPJ = output error flag,  = 0 if no trouble, .gt. 0 if
C         JAC returned an error flag.
C JCUR  = output flag = 1 to indicate that the Jacobian matrix
C         (or approximation) is now current.
C This routine also uses Common variables EL0, H, TN, IERPJ, JCUR, NPE.
C-----------------------------------------------------------------------
      INTEGER IER
      DOUBLE PRECISION HL0
C
      IERPJ = 0
      JCUR = 1
      HL0 = EL0*H
      CALL JAC (F, NEQ, TN, Y, YSV, EWT, SAVF, FTEM, HL0,
     1   WM(LOCWP), IWM(LOCIWP), IER)
      NPE = NPE + 1
      IF (IER .EQ. 0) RETURN
      IERPJ = 1
      RETURN
C----------------------- End of Subroutine DPKSET ----------------------
      END
*DECK DSOLPK
      SUBROUTINE DSOLPK (NEQ, Y, SAVF, X, EWT, WM, IWM, F, PSOL)
      EXTERNAL F, PSOL
      INTEGER NEQ, IWM
      DOUBLE PRECISION Y, SAVF, X, EWT, WM
      DIMENSION NEQ(*), Y(*), SAVF(*), X(*), EWT(*), WM(*), IWM(*)
      INTEGER ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, 
     1   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     2   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NPE, NQU
      INTEGER JPRE, JACFLG, LOCWP, LOCIWP, LSAVX, KMP, MAXL, MNEWT,
     1   NNI, NLI, NPS, NCFN, NCFL
      DOUBLE PRECISION CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
      DOUBLE PRECISION DELT, EPCON, SQRTN, RSQRTN
      COMMON /DLS001/ CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,
     1   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, 
     2   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     3   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NPE, NQU
      COMMON /DLPK01/ DELT, EPCON, SQRTN, RSQRTN,
     1   JPRE, JACFLG, LOCWP, LOCIWP, LSAVX, KMP, MAXL, MNEWT,
     2   NNI, NLI, NPS, NCFN, NCFL
C-----------------------------------------------------------------------
C This routine interfaces to one of DSPIOM, DSPIGMR, DPCG, DPCGS, or
C DUSOL, for the solution of the linear system arising from a Newton
C iteration.  It is called if MITER .ne. 0.
C In addition to variables described elsewhere,
C communication with DSOLPK uses the following variables:
C WM    = real work space containing data for the algorithm
C         (Krylov basis vectors, Hessenberg matrix, etc.)
C IWM   = integer work space containing data for the algorithm
C X     = the right-hand side vector on input, and the solution vector
C         on output, of length N.
C IERSL = output flag (in Common):
C         IERSL =  0 means no trouble occurred.
C         IERSL =  1 means the iterative method failed to converge.
C                    If the preconditioner is out of date, the step
C                    is repeated with a new preconditioner.
C                    Otherwise, the stepsize is reduced (forcing a
C                    new evaluation of the preconditioner) and the
C                    step is repeated.
C         IERSL = -1 means there was a nonrecoverable error in the
C                    iterative solver, and an error exit occurs.
C This routine also uses the Common variables TN, EL0, H, N, MITER,
C DELT, EPCON, SQRTN, RSQRTN, MAXL, KMP, MNEWT, NNI, NLI, NPS, NCFL,
C LOCWP, LOCIWP.
C-----------------------------------------------------------------------
      INTEGER IFLAG, LB, LDL, LHES, LIOM, LGMR, LPCG, LP, LQ, LR,
     1   LV, LW, LWK, LZ, MAXLP1, NPSL
      DOUBLE PRECISION DELTA, HL0
C
      IERSL = 0
      HL0 = H*EL0
      DELTA = DELT*EPCON
      GO TO (100, 200, 300, 400, 900, 900, 900, 900, 900), MITER
C-----------------------------------------------------------------------
C Use the SPIOM algorithm to solve the linear system P*x = -f.
C-----------------------------------------------------------------------
 100  CONTINUE
      LV = 1
      LB = LV + N*MAXL
      LHES = LB + N
      LWK = LHES + MAXL*MAXL
      CALL DCOPY (N, X, 1, WM(LB), 1)
      CALL DSCAL (N, RSQRTN, EWT, 1)
      CALL DSPIOM (NEQ, TN, Y, SAVF, WM(LB), EWT, N, MAXL, KMP, DELTA,
     1   HL0, JPRE, MNEWT, F, PSOL, NPSL, X, WM(LV), WM(LHES), IWM,
     2   LIOM, WM(LOCWP), IWM(LOCIWP), WM(LWK), IFLAG)
      NNI = NNI + 1
      NLI = NLI + LIOM
      NPS = NPS + NPSL
      CALL DSCAL (N, SQRTN, EWT, 1)
      IF (IFLAG .NE. 0) NCFL = NCFL + 1
      IF (IFLAG .GE. 2) IERSL = 1
      IF (IFLAG .LT. 0) IERSL = -1
      RETURN
C-----------------------------------------------------------------------
C Use the SPIGMR algorithm to solve the linear system P*x = -f.
C-----------------------------------------------------------------------
 200  CONTINUE
      MAXLP1 = MAXL + 1
      LV = 1
      LB = LV + N*MAXL
      LHES = LB + N + 1
      LQ = LHES + MAXL*MAXLP1
      LWK = LQ + 2*MAXL
      LDL = LWK + MIN(1,MAXL-KMP)*N
      CALL DCOPY (N, X, 1, WM(LB), 1)
      CALL DSCAL (N, RSQRTN, EWT, 1)
      CALL DSPIGMR (NEQ, TN, Y, SAVF, WM(LB), EWT, N, MAXL, MAXLP1, KMP,
     1   DELTA, HL0, JPRE, MNEWT, F, PSOL, NPSL, X, WM(LV), WM(LHES),
     2   WM(LQ), LGMR, WM(LOCWP), IWM(LOCIWP), WM(LWK), WM(LDL), IFLAG)
      NNI = NNI + 1
      NLI = NLI + LGMR
      NPS = NPS + NPSL
      CALL DSCAL (N, SQRTN, EWT, 1)
      IF (IFLAG .NE. 0) NCFL = NCFL + 1
      IF (IFLAG .GE. 2) IERSL = 1
      IF (IFLAG .LT. 0) IERSL = -1
      RETURN
C-----------------------------------------------------------------------
C Use DPCG to solve the linear system P*x = -f
C-----------------------------------------------------------------------
 300  CONTINUE
      LR = 1
      LP = LR + N
      LW = LP + N
      LZ = LW + N
      LWK = LZ + N
      CALL DCOPY (N, X, 1, WM(LR), 1)
      CALL DPCG (NEQ, TN, Y, SAVF, WM(LR), EWT, N, MAXL, DELTA, HL0,
     1          JPRE, MNEWT, F, PSOL, NPSL, X, WM(LP), WM(LW), WM(LZ),
     2          LPCG, WM(LOCWP), IWM(LOCIWP), WM(LWK), IFLAG)
      NNI = NNI + 1
      NLI = NLI + LPCG
      NPS = NPS + NPSL
      IF (IFLAG .NE. 0) NCFL = NCFL + 1
      IF (IFLAG .GE. 2) IERSL = 1
      IF (IFLAG .LT. 0) IERSL = -1
      RETURN
C-----------------------------------------------------------------------
C Use DPCGS to solve the linear system P*x = -f
C-----------------------------------------------------------------------
 400  CONTINUE
      LR = 1
      LP = LR + N
      LW = LP + N
      LZ = LW + N
      LWK = LZ + N
      CALL DCOPY (N, X, 1, WM(LR), 1)
      CALL DPCGS (NEQ, TN, Y, SAVF, WM(LR), EWT, N, MAXL, DELTA, HL0,
     1           JPRE, MNEWT, F, PSOL, NPSL, X, WM(LP), WM(LW), WM(LZ),
     2           LPCG, WM(LOCWP), IWM(LOCIWP), WM(LWK), IFLAG)
      NNI = NNI + 1
      NLI = NLI + LPCG
      NPS = NPS + NPSL
      IF (IFLAG .NE. 0) NCFL = NCFL + 1
      IF (IFLAG .GE. 2) IERSL = 1
      IF (IFLAG .LT. 0) IERSL = -1
      RETURN
C-----------------------------------------------------------------------
C Use DUSOL, which interfaces to PSOL, to solve the linear system
C (no Krylov iteration).
C-----------------------------------------------------------------------
 900  CONTINUE
      LB = 1
      LWK = LB + N
      CALL DCOPY (N, X, 1, WM(LB), 1)
      CALL DUSOL (NEQ, TN, Y, SAVF, WM(LB), EWT, N, DELTA, HL0, MNEWT,
     1   PSOL, NPSL, X, WM(LOCWP), IWM(LOCIWP), WM(LWK), IFLAG)
      NNI = NNI + 1
      NPS = NPS + NPSL
      IF (IFLAG .NE. 0) NCFL = NCFL + 1
      IF (IFLAG .EQ. 3) IERSL = 1
      IF (IFLAG .LT. 0) IERSL = -1
      RETURN
C----------------------- End of Subroutine DSOLPK ----------------------
      END
*DECK DSPIOM
      SUBROUTINE DSPIOM (NEQ, TN, Y, SAVF, B, WGHT, N, MAXL, KMP, DELTA,
     1            HL0, JPRE, MNEWT, F, PSOL, NPSL, X, V, HES, IPVT,
     2            LIOM, WP, IWP, WK, IFLAG)
      EXTERNAL F, PSOL
      INTEGER NEQ,N,MAXL,KMP,JPRE,MNEWT,NPSL,IPVT,LIOM,IWP,IFLAG
      DOUBLE PRECISION TN,Y,SAVF,B,WGHT,DELTA,HL0,X,V,HES,WP,WK
      DIMENSION NEQ(*), Y(*), SAVF(*), B(*), WGHT(*), X(*), V(N,*),
     1   HES(MAXL,MAXL), IPVT(*), WP(*), IWP(*), WK(*)
C-----------------------------------------------------------------------
C This routine solves the linear system A * x = b using a scaled
C preconditioned version of the Incomplete Orthogonalization Method.
C An initial guess of x = 0 is assumed.
C-----------------------------------------------------------------------
C
C      On entry
C
C          NEQ = problem size, passed to F and PSOL (NEQ(1) = N).
C
C           TN = current value of t.
C
C            Y = array containing current dependent variable vector.
C
C         SAVF = array containing current value of f(t,y).
C
C         B    = the right hand side of the system A*x = b.
C                B is also used as work space when computing the
C                final approximation.
C                (B is the same as V(*,MAXL+1) in the call to DSPIOM.)
C
C         WGHT = array of length N containing scale factors.
C                1/WGHT(i) are the diagonal elements of the diagonal
C                scaling matrix D.
C
C         N    = the order of the matrix A, and the lengths
C                of the vectors Y, SAVF, B, WGHT, and X.
C
C         MAXL = the maximum allowable order of the matrix HES.
C
C          KMP = the number of previous vectors the new vector VNEW
C                must be made orthogonal to.  KMP .le. MAXL.
C
C        DELTA = tolerance on residuals b - A*x in weighted RMS-norm.
C
C          HL0 = current value of (step size h) * (coefficient l0).
C
C         JPRE = preconditioner type flag.
C
C        MNEWT = Newton iteration counter (.ge. 0).
C
C           WK = real work array of length N used by DATV and PSOL.
C
C           WP = real work array used by preconditioner PSOL.
C
C          IWP = integer work array used by preconditioner PSOL.
C
C      On return
C
C         X    = the final computed approximation to the solution
C                of the system A*x = b.
C
C         V    = the N by (LIOM+1) array containing the LIOM
C                orthogonal vectors V(*,1) to V(*,LIOM).
C
C         HES  = the LU factorization of the LIOM by LIOM upper
C                Hessenberg matrix whose entries are the
C                scaled inner products of A*V(*,k) and V(*,i).
C
C         IPVT = an integer array containg pivoting information.
C                It is loaded in DHEFA and used in DHESL.
C
C         LIOM = the number of iterations performed, and current
C                order of the upper Hessenberg matrix HES.
C
C         NPSL = the number of calls to PSOL.
C
C        IFLAG = integer error flag:
C                0 means convergence in LIOM iterations, LIOM.le.MAXL.
C                1 means the convergence test did not pass in MAXL
C                  iterations, but the residual norm is .lt. 1,
C                  or .lt. norm(b) if MNEWT = 0, and so X is computed.
C                2 means the convergence test did not pass in MAXL
C                  iterations, residual .gt. 1, and X is undefined.
C                3 means there was a recoverable error in PSOL
C                  caused by the preconditioner being out of date.
C               -1 means there was a nonrecoverable error in PSOL.
C
C-----------------------------------------------------------------------
      INTEGER I, IER, INFO, J, K, LL, LM1
      DOUBLE PRECISION BNRM, BNRM0, PROD, RHO, SNORMW, DNRM2, TEM
C
      IFLAG = 0
      LIOM = 0
      NPSL = 0
C-----------------------------------------------------------------------
C The initial residual is the vector b.  Apply scaling to b, and test
C for an immediate return with X = 0 or X = b.
C-----------------------------------------------------------------------
      DO 10 I = 1,N
 10     V(I,1) = B(I)*WGHT(I)
      BNRM0 = DNRM2 (N, V, 1)
      BNRM = BNRM0
      IF (BNRM0 .GT. DELTA) GO TO 30
      IF (MNEWT .GT. 0) GO TO 20
      CALL DCOPY (N, B, 1, X, 1)
      RETURN
 20   DO 25 I = 1,N
 25     X(I) = 0.0D0
      RETURN
 30   CONTINUE
C Apply inverse of left preconditioner to vector b. --------------------
      IER = 0
      IF (JPRE .EQ. 0 .OR. JPRE .EQ. 2) GO TO 55
      CALL PSOL (NEQ, TN, Y, SAVF, WK, HL0, WP, IWP, B, 1, IER)
      NPSL = 1
      IF (IER .NE. 0) GO TO 300
C Calculate norm of scaled vector V(*,1) and normalize it. -------------
      DO 50 I = 1,N
 50     V(I,1) = B(I)*WGHT(I)
      BNRM = DNRM2(N, V, 1)
      DELTA = DELTA*(BNRM/BNRM0)
 55   TEM = 1.0D0/BNRM
      CALL DSCAL (N, TEM, V(1,1), 1)
C Zero out the HES array. ----------------------------------------------
      DO 65 J = 1,MAXL
        DO 60 I = 1,MAXL
 60       HES(I,J) = 0.0D0
 65     CONTINUE
C-----------------------------------------------------------------------
C Main loop on LL = l to compute the vectors V(*,2) to V(*,MAXL).
C The running product PROD is needed for the convergence test.
C-----------------------------------------------------------------------
      PROD = 1.0D0
      DO 90 LL = 1,MAXL
        LIOM = LL
C-----------------------------------------------------------------------
C Call routine DATV to compute VNEW = Abar*v(l), where Abar is
C the matrix A with scaling and inverse preconditioner factors applied.
C Call routine DORTHOG to orthogonalize the new vector vnew = V(*,l+1).
C Call routine DHEFA to update the factors of HES.
C-----------------------------------------------------------------------
        CALL DATV (NEQ, Y, SAVF, V(1,LL), WGHT, X, F, PSOL, V(1,LL+1),
     1        WK, WP, IWP, HL0, JPRE, IER, NPSL)
        IF (IER .NE. 0) GO TO 300
        CALL DORTHOG (V(1,LL+1), V, HES, N, LL, MAXL, KMP, SNORMW)
        CALL DHEFA (HES, MAXL, LL, IPVT, INFO, LL)
        LM1 = LL - 1
        IF (LL .GT. 1 .AND. IPVT(LM1) .EQ. LM1) PROD = PROD*HES(LL,LM1)
        IF (INFO .NE. LL) GO TO 70
C-----------------------------------------------------------------------
C The last pivot in HES was found to be zero.
C If vnew = 0 or l = MAXL, take an error return with IFLAG = 2.
C otherwise, continue the iteration without a convergence test.
C-----------------------------------------------------------------------
        IF (SNORMW .EQ. 0.0D0) GO TO 120
        IF (LL .EQ. MAXL) GO TO 120
        GO TO 80
C-----------------------------------------------------------------------
C Update RHO, the estimate of the norm of the residual b - A*x(l).
C test for convergence.  If passed, compute approximation x(l).
C If failed and l .lt. MAXL, then continue iterating.
C-----------------------------------------------------------------------
 70     CONTINUE
        RHO = BNRM*SNORMW*ABS(PROD/HES(LL,LL))
        IF (RHO .LE. DELTA) GO TO 200
        IF (LL .EQ. MAXL) GO TO 100
C If l .lt. MAXL, store HES(l+1,l) and normalize the vector v(*,l+1).
 80     CONTINUE
        HES(LL+1,LL) = SNORMW
        TEM = 1.0D0/SNORMW
        CALL DSCAL (N, TEM, V(1,LL+1), 1)
 90     CONTINUE
C-----------------------------------------------------------------------
C l has reached MAXL without passing the convergence test:
C If RHO is not too large, compute a solution anyway and return with
C IFLAG = 1.  Otherwise return with IFLAG = 2.
C-----------------------------------------------------------------------
 100  CONTINUE
      IF (RHO .LE. 1.0D0) GO TO 150
      IF (RHO .LE. BNRM .AND. MNEWT .EQ. 0) GO TO 150
 120  CONTINUE
      IFLAG = 2
      RETURN
 150  IFLAG = 1
C-----------------------------------------------------------------------
C Compute the approximation x(l) to the solution.
C Since the vector X was used as work space, and the initial guess
C of the Newton correction is zero, X must be reset to zero.
C-----------------------------------------------------------------------
 200  CONTINUE
      LL = LIOM
      DO 210 K = 1,LL
 210    B(K) = 0.0D0
      B(1) = BNRM
      CALL DHESL (HES, MAXL, LL, IPVT, B)
      DO 220 K = 1,N
 220    X(K) = 0.0D0
      DO 230 I = 1,LL
        CALL DAXPY (N, B(I), V(1,I), 1, X, 1)
 230    CONTINUE
      DO 240 I = 1,N
 240    X(I) = X(I)/WGHT(I)
      IF (JPRE .LE. 1) RETURN
      CALL PSOL (NEQ, TN, Y, SAVF, WK, HL0, WP, IWP, X, 2, IER)
      NPSL = NPSL + 1
      IF (IER .NE. 0) GO TO 300
      RETURN
C-----------------------------------------------------------------------
C This block handles error returns forced by routine PSOL.
C-----------------------------------------------------------------------
 300  CONTINUE
      IF (IER .LT. 0) IFLAG = -1
      IF (IER .GT. 0) IFLAG = 3
      RETURN
C----------------------- End of Subroutine DSPIOM ----------------------
      END
*DECK DATV
      SUBROUTINE DATV (NEQ, Y, SAVF, V, WGHT, FTEM, F, PSOL, Z, VTEM,
     1                WP, IWP, HL0, JPRE, IER, NPSL)
      EXTERNAL F, PSOL
      INTEGER NEQ, IWP, JPRE, IER, NPSL
      DOUBLE PRECISION Y, SAVF, V, WGHT, FTEM, Z, VTEM, WP, HL0
      DIMENSION NEQ(*), Y(*), SAVF(*), V(*), WGHT(*), FTEM(*), Z(*),
     1   VTEM(*), WP(*), IWP(*)
      INTEGER ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, 
     1   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     2   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NPE, NQU
      DOUBLE PRECISION CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
      COMMON /DLS001/ CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,
     1   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, 
     2   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     3   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NPE, NQU
C-----------------------------------------------------------------------
C This routine computes the product
C
C   (D-inverse)*(P1-inverse)*(I - hl0*df/dy)*(P2-inverse)*(D*v),
C
C where D is a diagonal scaling matrix, and P1 and P2 are the
C left and right preconditioning matrices, respectively.
C v is assumed to have WRMS norm equal to 1.
C The product is stored in z.  This is computed by a
C difference quotient, a call to F, and two calls to PSOL.
C-----------------------------------------------------------------------
C
C      On entry
C
C          NEQ = problem size, passed to F and PSOL (NEQ(1) = N).
C
C            Y = array containing current dependent variable vector.
C
C         SAVF = array containing current value of f(t,y).
C
C            V = real array of length N (can be the same array as Z).
C
C         WGHT = array of length N containing scale factors.
C                1/WGHT(i) are the diagonal elements of the matrix D.
C
C         FTEM = work array of length N.
C
C         VTEM = work array of length N used to store the
C                unscaled version of V.
C
C           WP = real work array used by preconditioner PSOL.
C
C          IWP = integer work array used by preconditioner PSOL.
C
C          HL0 = current value of (step size h) * (coefficient l0).
C
C         JPRE = preconditioner type flag.
C
C
C      On return
C
C            Z = array of length N containing desired scaled
C                matrix-vector product.
C
C          IER = error flag from PSOL.
C
C         NPSL = the number of calls to PSOL.
C
C In addition, this routine uses the Common variables TN, N, NFE.
C-----------------------------------------------------------------------
      INTEGER I
      DOUBLE PRECISION FAC, RNORM, DNRM2, TEMPN
C
C Set VTEM = D * V.
      DO 10 I = 1,N
 10     VTEM(I) = V(I)/WGHT(I)
      IER = 0
      IF (JPRE .GE. 2) GO TO 30
C
C JPRE = 0 or 1.  Save Y in Z and increment Y by VTEM.
      CALL DCOPY (N, Y, 1, Z, 1)
      DO 20 I = 1,N
 20     Y(I) = Z(I) + VTEM(I)
      FAC = HL0
      GO TO 60
C
C JPRE = 2 or 3.  Apply inverse of right preconditioner to VTEM.
 30   CONTINUE
      CALL PSOL (NEQ, TN, Y, SAVF, FTEM, HL0, WP, IWP, VTEM, 2, IER)
      NPSL = NPSL + 1
      IF (IER .NE. 0) RETURN
C Calculate L-2 norm of (D-inverse) * VTEM.
      DO 40 I = 1,N
 40     Z(I) = VTEM(I)*WGHT(I)
      TEMPN = DNRM2 (N, Z, 1)
      RNORM = 1.0D0/TEMPN
C Save Y in Z and increment Y by VTEM/norm.
      CALL DCOPY (N, Y, 1, Z, 1)
      DO 50 I = 1,N
 50     Y(I) = Z(I) + VTEM(I)*RNORM
      FAC = HL0*TEMPN
C
C For all JPRE, call F with incremented Y argument, and restore Y.
 60   CONTINUE
      CALL F (NEQ, TN, Y, FTEM)
      NFE = NFE + 1
      CALL DCOPY (N, Z, 1, Y, 1)
C Set Z = (identity - hl0*Jacobian) * VTEM, using difference quotient.
      DO 70 I = 1,N
 70     Z(I) = FTEM(I) - SAVF(I)
      DO 80 I = 1,N
 80     Z(I) = VTEM(I) - FAC*Z(I)
C Apply inverse of left preconditioner to Z, if nontrivial.
      IF (JPRE .EQ. 0 .OR. JPRE .EQ. 2) GO TO 85
      CALL PSOL (NEQ, TN, Y, SAVF, FTEM, HL0, WP, IWP, Z, 1, IER)
      NPSL = NPSL + 1
      IF (IER .NE. 0) RETURN
 85   CONTINUE
C Apply D-inverse to Z and return.
      DO 90 I = 1,N
 90     Z(I) = Z(I)*WGHT(I)
      RETURN
C----------------------- End of Subroutine DATV ------------------------
      END
*DECK DORTHOG
      SUBROUTINE DORTHOG (VNEW, V, HES, N, LL, LDHES, KMP, SNORMW)
      INTEGER N, LL, LDHES, KMP
      DOUBLE PRECISION VNEW, V, HES, SNORMW
      DIMENSION VNEW(*), V(N,*), HES(LDHES,*)
C-----------------------------------------------------------------------
C This routine orthogonalizes the vector VNEW against the previous
C KMP vectors in the V array.  It uses a modified Gram-Schmidt
C orthogonalization procedure with conditional reorthogonalization.
C This is the version of 28 may 1986.
C-----------------------------------------------------------------------
C
C      On entry
C
C         VNEW = the vector of length N containing a scaled product
C                of the Jacobian and the vector V(*,LL).
C
C         V    = the N x l array containing the previous LL
C                orthogonal vectors v(*,1) to v(*,LL).
C
C         HES  = an LL x LL upper Hessenberg matrix containing,
C                in HES(i,k), k.lt.LL, scaled inner products of
C                A*V(*,k) and V(*,i).
C
C        LDHES = the leading dimension of the HES array.
C
C         N    = the order of the matrix A, and the length of VNEW.
C
C         LL   = the current order of the matrix HES.
C
C          KMP = the number of previous vectors the new vector VNEW
C                must be made orthogonal to (KMP .le. MAXL).
C
C
C      On return
C
C         VNEW = the new vector orthogonal to V(*,i0) to V(*,LL),
C                where i0 = MAX(1, LL-KMP+1).
C
C         HES  = upper Hessenberg matrix with column LL filled in with
C                scaled inner products of A*V(*,LL) and V(*,i).
C
C       SNORMW = L-2 norm of VNEW.
C
C-----------------------------------------------------------------------
      INTEGER I, I0
      DOUBLE PRECISION ARG, DDOT, DNRM2, SUMDSQ, TEM, VNRM
C
C Get norm of unaltered VNEW for later use. ----------------------------
      VNRM = DNRM2 (N, VNEW, 1)
C-----------------------------------------------------------------------
C Do modified Gram-Schmidt on VNEW = A*v(LL).
C Scaled inner products give new column of HES.
C Projections of earlier vectors are subtracted from VNEW.
C-----------------------------------------------------------------------
      I0 = MAX(1,LL-KMP+1)
      DO 10 I = I0,LL
        HES(I,LL) = DDOT (N, V(1,I), 1, VNEW, 1)
        TEM = -HES(I,LL)
        CALL DAXPY (N, TEM, V(1,I), 1, VNEW, 1)
 10     CONTINUE
C-----------------------------------------------------------------------
C Compute SNORMW = norm of VNEW.
C If VNEW is small compared to its input value (in norm), then
C reorthogonalize VNEW to V(*,1) through V(*,LL).
C Correct if relative correction exceeds 1000*(unit roundoff).
C finally, correct SNORMW using the dot products involved.
C-----------------------------------------------------------------------
      SNORMW = DNRM2 (N, VNEW, 1)
      IF (VNRM + 0.001D0*SNORMW .NE. VNRM) RETURN
      SUMDSQ = 0.0D0
      DO 30 I = I0,LL
        TEM = -DDOT (N, V(1,I), 1, VNEW, 1)
        IF (HES(I,LL) + 0.001D0*TEM .EQ. HES(I,LL)) GO TO 30
        HES(I,LL) = HES(I,LL) - TEM
        CALL DAXPY (N, TEM, V(1,I), 1, VNEW, 1)
        SUMDSQ = SUMDSQ + TEM**2
 30     CONTINUE
      IF (SUMDSQ .EQ. 0.0D0) RETURN
      ARG = MAX(0.0D0,SNORMW**2 - SUMDSQ)
      SNORMW = SQRT(ARG)
C
      RETURN
C----------------------- End of Subroutine DORTHOG ---------------------
      END
*DECK DSPIGMR
      SUBROUTINE DSPIGMR (NEQ, TN, Y, SAVF, B, WGHT, N, MAXL, MAXLP1,
     1  KMP, DELTA, HL0, JPRE, MNEWT, F, PSOL, NPSL, X, V, HES, Q,
     2  LGMR, WP, IWP, WK, DL, IFLAG)
      EXTERNAL F, PSOL
      INTEGER NEQ,N,MAXL,MAXLP1,KMP,JPRE,MNEWT,NPSL,LGMR,IWP,IFLAG
      DOUBLE PRECISION TN,Y,SAVF,B,WGHT,DELTA,HL0,X,V,HES,Q,WP,WK,DL
      DIMENSION NEQ(*), Y(*), SAVF(*), B(*), WGHT(*), X(*), V(N,*),
     1    HES(MAXLP1,*), Q(*), WP(*), IWP(*), WK(*), DL(*)
C-----------------------------------------------------------------------
C This routine solves the linear system A * x = b using a scaled
C preconditioned version of the Generalized Minimal Residual method.
C An initial guess of x = 0 is assumed.
C-----------------------------------------------------------------------
C
C      On entry
C
C          NEQ = problem size, passed to F and PSOL (NEQ(1) = N).
C
C           TN = current value of t.
C
C            Y = array containing current dependent variable vector.
C
C         SAVF = array containing current value of f(t,y).
C
C            B = the right hand side of the system A*x = b.
C                B is also used as work space when computing
C                the final approximation.
C                (B is the same as V(*,MAXL+1) in the call to DSPIGMR.)
C
C         WGHT = the vector of length N containing the nonzero
C                elements of the diagonal scaling matrix.
C
C            N = the order of the matrix A, and the lengths
C                of the vectors WGHT, B and X.
C
C         MAXL = the maximum allowable order of the matrix HES.
C
C       MAXLP1 = MAXL + 1, used for dynamic dimensioning of HES.
C
C          KMP = the number of previous vectors the new vector VNEW
C                must be made orthogonal to.  KMP .le. MAXL.
C
C        DELTA = tolerance on residuals b - A*x in weighted RMS-norm.
C
C          HL0 = current value of (step size h) * (coefficient l0).
C
C         JPRE = preconditioner type flag.
C
C        MNEWT = Newton iteration counter (.ge. 0).
C
C           WK = real work array used by routine DATV and PSOL.
C
C           DL = real work array used for calculation of the residual
C                norm RHO when the method is incomplete (KMP .lt. MAXL).
C                Not needed or referenced in complete case (KMP = MAXL).
C
C           WP = real work array used by preconditioner PSOL.
C
C          IWP = integer work array used by preconditioner PSOL.
C
C      On return
C
C         X    = the final computed approximation to the solution
C                of the system A*x = b.
C
C         LGMR = the number of iterations performed and
C                the current order of the upper Hessenberg
C                matrix HES.
C
C         NPSL = the number of calls to PSOL.
C
C         V    = the N by (LGMR+1) array containing the LGMR
C                orthogonal vectors V(*,1) to V(*,LGMR).
C
C         HES  = the upper triangular factor of the QR decomposition
C                of the (LGMR+1) by lgmr upper Hessenberg matrix whose
C                entries are the scaled inner-products of A*V(*,i)
C                and V(*,k).
C
C         Q    = real array of length 2*MAXL containing the components
C                of the Givens rotations used in the QR decomposition
C                of HES.  It is loaded in DHEQR and used in DHELS.
C
C        IFLAG = integer error flag:
C                0 means convergence in LGMR iterations, LGMR .le. MAXL.
C                1 means the convergence test did not pass in MAXL
C                  iterations, but the residual norm is .lt. 1,
C                  or .lt. norm(b) if MNEWT = 0, and so x is computed.
C                2 means the convergence test did not pass in MAXL
C                  iterations, residual .gt. 1, and X is undefined.
C                3 means there was a recoverable error in PSOL
C                  caused by the preconditioner being out of date.
C               -1 means there was a nonrecoverable error in PSOL.
C
C-----------------------------------------------------------------------
      INTEGER I, IER, INFO, IP1, I2, J, K, LL, LLP1
      DOUBLE PRECISION BNRM,BNRM0,C,DLNRM,PROD,RHO,S,SNORMW,DNRM2,TEM
C
      IFLAG = 0
      LGMR = 0
      NPSL = 0
C-----------------------------------------------------------------------
C The initial residual is the vector b.  Apply scaling to b, and test
C for an immediate return with X = 0 or X = b.
C-----------------------------------------------------------------------
      DO 10 I = 1,N
 10     V(I,1) = B(I)*WGHT(I)
      BNRM0 = DNRM2 (N, V, 1)
      BNRM = BNRM0
      IF (BNRM0 .GT. DELTA) GO TO 30
      IF (MNEWT .GT. 0) GO TO 20
      CALL DCOPY (N, B, 1, X, 1)
      RETURN
 20   DO 25 I = 1,N
 25     X(I) = 0.0D0
      RETURN
 30   CONTINUE
C Apply inverse of left preconditioner to vector b. --------------------
      IER = 0
      IF (JPRE .EQ. 0 .OR. JPRE .EQ. 2) GO TO 55
      CALL PSOL (NEQ, TN, Y, SAVF, WK, HL0, WP, IWP, B, 1, IER)
      NPSL = 1
      IF (IER .NE. 0) GO TO 300
C Calculate norm of scaled vector V(*,1) and normalize it. -------------
      DO 50 I = 1,N
 50     V(I,1) = B(I)*WGHT(I)
      BNRM = DNRM2 (N, V, 1)
      DELTA = DELTA*(BNRM/BNRM0)
 55   TEM = 1.0D0/BNRM
      CALL DSCAL (N, TEM, V(1,1), 1)
C Zero out the HES array. ----------------------------------------------
      DO 65 J = 1,MAXL
        DO 60 I = 1,MAXLP1
 60       HES(I,J) = 0.0D0
 65     CONTINUE
C-----------------------------------------------------------------------
C Main loop to compute the vectors V(*,2) to V(*,MAXL).
C The running product PROD is needed for the convergence test.
C-----------------------------------------------------------------------
      PROD = 1.0D0
      DO 90 LL = 1,MAXL
        LGMR = LL
C-----------------------------------------------------------------------
C Call routine DATV to compute VNEW = Abar*v(ll), where Abar is
C the matrix A with scaling and inverse preconditioner factors applied.
C Call routine DORTHOG to orthogonalize the new vector VNEW = V(*,LL+1).
C Call routine DHEQR to update the factors of HES.
C-----------------------------------------------------------------------
        CALL DATV (NEQ, Y, SAVF, V(1,LL), WGHT, X, F, PSOL, V(1,LL+1),
     1        WK, WP, IWP, HL0, JPRE, IER, NPSL)
        IF (IER .NE. 0) GO TO 300
        CALL DORTHOG (V(1,LL+1), V, HES, N, LL, MAXLP1, KMP, SNORMW)
        HES(LL+1,LL) = SNORMW
        CALL DHEQR (HES, MAXLP1, LL, Q, INFO, LL)
        IF (INFO .EQ. LL) GO TO 120
C-----------------------------------------------------------------------
C Update RHO, the estimate of the norm of the residual b - A*xl.
C If KMP .lt. MAXL, then the vectors V(*,1),...,V(*,LL+1) are not
C necessarily orthogonal for LL .gt. KMP.  The vector DL must then
C be computed, and its norm used in the calculation of RHO.
C-----------------------------------------------------------------------
        PROD = PROD*Q(2*LL)
        RHO = ABS(PROD*BNRM)
        IF ((LL.GT.KMP) .AND. (KMP.LT.MAXL)) THEN
          IF (LL .EQ. KMP+1) THEN
            CALL DCOPY (N, V(1,1), 1, DL, 1)
            DO 75 I = 1,KMP
              IP1 = I + 1
              I2 = I*2
              S = Q(I2)
              C = Q(I2-1)
              DO 70 K = 1,N
 70             DL(K) = S*DL(K) + C*V(K,IP1)
 75           CONTINUE
            ENDIF
          S = Q(2*LL)
          C = Q(2*LL-1)/SNORMW
          LLP1 = LL + 1
          DO 80 K = 1,N
 80         DL(K) = S*DL(K) + C*V(K,LLP1)
          DLNRM = DNRM2 (N, DL, 1)
          RHO = RHO*DLNRM
          ENDIF
C-----------------------------------------------------------------------
C Test for convergence.  If passed, compute approximation xl.
C if failed and LL .lt. MAXL, then continue iterating.
C-----------------------------------------------------------------------
        IF (RHO .LE. DELTA) GO TO 200
        IF (LL .EQ. MAXL) GO TO 100
C-----------------------------------------------------------------------
C Rescale so that the norm of V(1,LL+1) is one.
C-----------------------------------------------------------------------
        TEM = 1.0D0/SNORMW
        CALL DSCAL (N, TEM, V(1,LL+1), 1)
 90     CONTINUE
 100  CONTINUE
      IF (RHO .LE. 1.0D0) GO TO 150
      IF (RHO .LE. BNRM .AND. MNEWT .EQ. 0) GO TO 150
 120  CONTINUE
      IFLAG = 2
      RETURN
 150  IFLAG = 1
C-----------------------------------------------------------------------
C Compute the approximation xl to the solution.
C Since the vector X was used as work space, and the initial guess
C of the Newton correction is zero, X must be reset to zero.
C-----------------------------------------------------------------------
 200  CONTINUE
      LL = LGMR
      LLP1 = LL + 1
      DO 210 K = 1,LLP1
 210    B(K) = 0.0D0
      B(1) = BNRM
      CALL DHELS (HES, MAXLP1, LL, Q, B)
      DO 220 K = 1,N
 220    X(K) = 0.0D0
      DO 230 I = 1,LL
        CALL DAXPY (N, B(I), V(1,I), 1, X, 1)
 230    CONTINUE
      DO 240 I = 1,N
 240    X(I) = X(I)/WGHT(I)
      IF (JPRE .LE. 1) RETURN
      CALL PSOL (NEQ, TN, Y, SAVF, WK, HL0, WP, IWP, X, 2, IER)
      NPSL = NPSL + 1
      IF (IER .NE. 0) GO TO 300
      RETURN
C-----------------------------------------------------------------------
C This block handles error returns forced by routine PSOL.
C-----------------------------------------------------------------------
 300  CONTINUE
      IF (IER .LT. 0) IFLAG = -1
      IF (IER .GT. 0) IFLAG = 3
C
      RETURN
C----------------------- End of Subroutine DSPIGMR ---------------------
      END
*DECK DPCG
      SUBROUTINE DPCG (NEQ, TN, Y, SAVF, R, WGHT, N, MAXL, DELTA, HL0,
     1 JPRE, MNEWT, F, PSOL, NPSL, X, P, W, Z, LPCG, WP, IWP, WK, IFLAG)
      EXTERNAL F, PSOL
      INTEGER NEQ, N, MAXL, JPRE, MNEWT, NPSL, LPCG, IWP, IFLAG
      DOUBLE PRECISION TN,Y,SAVF,R,WGHT,DELTA,HL0,X,P,W,Z,WP,WK
      DIMENSION NEQ(*), Y(*), SAVF(*), R(*), WGHT(*), X(*), P(*), W(*),
     1   Z(*), WP(*), IWP(*), WK(*)
C-----------------------------------------------------------------------
C This routine computes the solution to the system A*x = b using a
C preconditioned version of the Conjugate Gradient algorithm.
C It is assumed here that the matrix A and the preconditioner
C matrix M are symmetric positive definite or nearly so.
C-----------------------------------------------------------------------
C
C      On entry
C
C          NEQ = problem size, passed to F and PSOL (NEQ(1) = N).
C
C           TN = current value of t.
C
C            Y = array containing current dependent variable vector.
C
C         SAVF = array containing current value of f(t,y).
C
C            R = the right hand side of the system A*x = b.
C
C         WGHT = array of length N containing scale factors.
C                1/WGHT(i) are the diagonal elements of the diagonal
C                scaling matrix D.
C
C            N = the order of the matrix A, and the lengths
C                of the vectors Y, SAVF, R, WGHT, P, W, Z, WK, and X.
C
C         MAXL = the maximum allowable number of iterates.
C
C        DELTA = tolerance on residuals b - A*x in weighted RMS-norm.
C
C          HL0 = current value of (step size h) * (coefficient l0).
C
C         JPRE = preconditioner type flag.
C
C        MNEWT = Newton iteration counter (.ge. 0).
C
C           WK = real work array used by routine DATP.
C
C           WP = real work array used by preconditioner PSOL.
C
C          IWP = integer work array used by preconditioner PSOL.
C
C      On return
C
C         X    = the final computed approximation to the solution
C                of the system A*x = b.
C
C         LPCG = the number of iterations performed, and current
C                order of the upper Hessenberg matrix HES.
C
C         NPSL = the number of calls to PSOL.
C
C        IFLAG = integer error flag:
C                0 means convergence in LPCG iterations, LPCG .le. MAXL.
C                1 means the convergence test did not pass in MAXL
C                  iterations, but the residual norm is .lt. 1,
C                  or .lt. norm(b) if MNEWT = 0, and so X is computed.
C                2 means the convergence test did not pass in MAXL
C                  iterations, residual .gt. 1, and X is undefined.
C                3 means there was a recoverable error in PSOL
C                  caused by the preconditioner being out of date.
C                4 means there was a zero denominator in the algorithm.
C                  The system matrix or preconditioner matrix is not
C                  sufficiently close to being symmetric pos. definite.
C               -1 means there was a nonrecoverable error in PSOL.
C
C-----------------------------------------------------------------------
      INTEGER I, IER
      DOUBLE PRECISION ALPHA,BETA,BNRM,PTW,RNRM,DDOT,DVNORM,ZTR,ZTR0
C
      IFLAG = 0
      NPSL = 0
      LPCG = 0
      DO 10 I = 1,N
 10     X(I) = 0.0D0
      BNRM = DVNORM (N, R, WGHT)
C Test for immediate return with X = 0 or X = b. -----------------------
      IF (BNRM .GT. DELTA) GO TO 20
      IF (MNEWT .GT. 0) RETURN
      CALL DCOPY (N, R, 1, X, 1)
      RETURN
C
 20   ZTR = 0.0D0
C Loop point for PCG iterations. ---------------------------------------
 30   CONTINUE
      LPCG = LPCG + 1
      CALL DCOPY (N, R, 1, Z, 1)
      IER = 0
      IF (JPRE .EQ. 0) GO TO 40
      CALL PSOL (NEQ, TN, Y, SAVF, WK, HL0, WP, IWP, Z, 3, IER)
      NPSL = NPSL + 1
      IF (IER .NE. 0) GO TO 100
 40   CONTINUE
      ZTR0 = ZTR
      ZTR = DDOT (N, Z, 1, R, 1)
      IF (LPCG .NE. 1) GO TO 50
      CALL DCOPY (N, Z, 1, P, 1)
      GO TO 70
 50   CONTINUE
      IF (ZTR0 .EQ. 0.0D0) GO TO 200
      BETA = ZTR/ZTR0
      DO 60 I = 1,N
 60     P(I) = Z(I) + BETA*P(I)
 70   CONTINUE
C-----------------------------------------------------------------------
C  Call DATP to compute A*p and return the answer in W.
C-----------------------------------------------------------------------
      CALL DATP (NEQ, Y, SAVF, P, WGHT, HL0, WK, F, W)
C
      PTW = DDOT (N, P, 1, W, 1)
      IF (PTW .EQ. 0.0D0) GO TO 200
      ALPHA = ZTR/PTW
      CALL DAXPY (N, ALPHA, P, 1, X, 1)
      ALPHA = -ALPHA
      CALL DAXPY (N, ALPHA, W, 1, R, 1)
      RNRM = DVNORM (N, R, WGHT)
      IF (RNRM .LE. DELTA) RETURN
      IF (LPCG .LT. MAXL) GO TO 30
      IFLAG = 2
      IF (RNRM .LE. 1.0D0) IFLAG = 1
      IF (RNRM .LE. BNRM .AND. MNEWT .EQ. 0) IFLAG = 1
      RETURN
C-----------------------------------------------------------------------
C This block handles error returns from PSOL.
C-----------------------------------------------------------------------
 100  CONTINUE
      IF (IER .LT. 0) IFLAG = -1
      IF (IER .GT. 0) IFLAG = 3
      RETURN
C-----------------------------------------------------------------------
C This block handles division by zero errors.
C-----------------------------------------------------------------------
 200  CONTINUE
      IFLAG = 4
      RETURN
C----------------------- End of Subroutine DPCG ------------------------
      END
*DECK DPCGS
      SUBROUTINE DPCGS (NEQ, TN, Y, SAVF, R, WGHT, N, MAXL, DELTA, HL0,
     1 JPRE, MNEWT, F, PSOL, NPSL, X, P, W, Z, LPCG, WP, IWP, WK, IFLAG)
      EXTERNAL F, PSOL
      INTEGER NEQ, N, MAXL, JPRE, MNEWT, NPSL, LPCG, IWP, IFLAG
      DOUBLE PRECISION TN,Y,SAVF,R,WGHT,DELTA,HL0,X,P,W,Z,WP,WK
      DIMENSION NEQ(*), Y(*), SAVF(*), R(*), WGHT(*), X(*), P(*), W(*),
     1   Z(*), WP(*), IWP(*), WK(*)
C-----------------------------------------------------------------------
C This routine computes the solution to the system A*x = b using a
C scaled preconditioned version of the Conjugate Gradient algorithm.
C It is assumed here that the scaled matrix D**-1 * A * D and the
C scaled preconditioner D**-1 * M * D are close to being
C symmetric positive definite.
C-----------------------------------------------------------------------
C
C      On entry
C
C          NEQ = problem size, passed to F and PSOL (NEQ(1) = N).
C
C           TN = current value of t.
C
C            Y = array containing current dependent variable vector.
C
C         SAVF = array containing current value of f(t,y).
C
C            R = the right hand side of the system A*x = b.
C
C         WGHT = array of length N containing scale factors.
C                1/WGHT(i) are the diagonal elements of the diagonal
C                scaling matrix D.
C
C            N = the order of the matrix A, and the lengths
C                of the vectors Y, SAVF, R, WGHT, P, W, Z, WK, and X.
C
C         MAXL = the maximum allowable number of iterates.
C
C        DELTA = tolerance on residuals b - A*x in weighted RMS-norm.
C
C          HL0 = current value of (step size h) * (coefficient l0).
C
C         JPRE = preconditioner type flag.
C
C        MNEWT = Newton iteration counter (.ge. 0).
C
C           WK = real work array used by routine DATP.
C
C           WP = real work array used by preconditioner PSOL.
C
C          IWP = integer work array used by preconditioner PSOL.
C
C      On return
C
C         X    = the final computed approximation to the solution
C                of the system A*x = b.
C
C         LPCG = the number of iterations performed, and current
C                order of the upper Hessenberg matrix HES.
C
C         NPSL = the number of calls to PSOL.
C
C        IFLAG = integer error flag:
C                0 means convergence in LPCG iterations, LPCG .le. MAXL.
C                1 means the convergence test did not pass in MAXL
C                  iterations, but the residual norm is .lt. 1,
C                  or .lt. norm(b) if MNEWT = 0, and so X is computed.
C                2 means the convergence test did not pass in MAXL
C                  iterations, residual .gt. 1, and X is undefined.
C                3 means there was a recoverable error in PSOL
C                  caused by the preconditioner being out of date.
C                4 means there was a zero denominator in the algorithm.
C                  the scaled matrix or scaled preconditioner is not
C                  sufficiently close to being symmetric pos. definite.
C               -1 means there was a nonrecoverable error in PSOL.
C
C-----------------------------------------------------------------------
      INTEGER I, IER
      DOUBLE PRECISION ALPHA, BETA, BNRM, PTW, RNRM, DVNORM, ZTR, ZTR0
C
      IFLAG = 0
      NPSL = 0
      LPCG = 0
      DO 10 I = 1,N
 10     X(I) = 0.0D0
      BNRM = DVNORM (N, R, WGHT)
C Test for immediate return with X = 0 or X = b. -----------------------
      IF (BNRM .GT. DELTA) GO TO 20
      IF (MNEWT .GT. 0) RETURN
      CALL DCOPY (N, R, 1, X, 1)
      RETURN
C
 20   ZTR = 0.0D0
C Loop point for PCG iterations. ---------------------------------------
 30   CONTINUE
      LPCG = LPCG + 1
      CALL DCOPY (N, R, 1, Z, 1)
      IER = 0
      IF (JPRE .EQ. 0) GO TO 40
      CALL PSOL (NEQ, TN, Y, SAVF, WK, HL0, WP, IWP, Z, 3, IER)
      NPSL = NPSL + 1
      IF (IER .NE. 0) GO TO 100
 40   CONTINUE
      ZTR0 = ZTR
      ZTR = 0.0D0
      DO 45 I = 1,N
 45     ZTR = ZTR + Z(I)*R(I)*WGHT(I)**2
      IF (LPCG .NE. 1) GO TO 50
      CALL DCOPY (N, Z, 1, P, 1)
      GO TO 70
 50   CONTINUE
      IF (ZTR0 .EQ. 0.0D0) GO TO 200
      BETA = ZTR/ZTR0
      DO 60 I = 1,N
 60     P(I) = Z(I) + BETA*P(I)
 70   CONTINUE
C-----------------------------------------------------------------------
C  Call DATP to compute A*p and return the answer in W.
C-----------------------------------------------------------------------
      CALL DATP (NEQ, Y, SAVF, P, WGHT, HL0, WK, F, W)
C
      PTW = 0.0D0
      DO 80 I = 1,N
 80     PTW = PTW + P(I)*W(I)*WGHT(I)**2
      IF (PTW .EQ. 0.0D0) GO TO 200
      ALPHA = ZTR/PTW
      CALL DAXPY (N, ALPHA, P, 1, X, 1)
      ALPHA = -ALPHA
      CALL DAXPY (N, ALPHA, W, 1, R, 1)
      RNRM = DVNORM (N, R, WGHT)
      IF (RNRM .LE. DELTA) RETURN
      IF (LPCG .LT. MAXL) GO TO 30
      IFLAG = 2
      IF (RNRM .LE. 1.0D0) IFLAG = 1
      IF (RNRM .LE. BNRM .AND. MNEWT .EQ. 0) IFLAG = 1
      RETURN
C-----------------------------------------------------------------------
C This block handles error returns from PSOL.
C-----------------------------------------------------------------------
 100  CONTINUE
      IF (IER .LT. 0) IFLAG = -1
      IF (IER .GT. 0) IFLAG = 3
      RETURN
C-----------------------------------------------------------------------
C This block handles division by zero errors.
C-----------------------------------------------------------------------
 200  CONTINUE
      IFLAG = 4
      RETURN
C----------------------- End of Subroutine DPCGS -----------------------
      END
*DECK DATP
      SUBROUTINE DATP (NEQ, Y, SAVF, P, WGHT, HL0, WK, F, W)
      EXTERNAL F
      INTEGER NEQ
      DOUBLE PRECISION Y, SAVF, P, WGHT, HL0, WK, W
      DIMENSION NEQ(*), Y(*), SAVF(*), P(*), WGHT(*), WK(*), W(*)
      INTEGER ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, 
     1   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     2   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NPE, NQU
      DOUBLE PRECISION CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
      COMMON /DLS001/ CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,
     1   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, 
     2   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     3   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NPE, NQU
C-----------------------------------------------------------------------
C This routine computes the product
C
C              w = (I - hl0*df/dy)*p
C
C This is computed by a call to F and a difference quotient.
C-----------------------------------------------------------------------
C
C      On entry
C
C          NEQ = problem size, passed to F and PSOL (NEQ(1) = N).
C
C            Y = array containing current dependent variable vector.
C
C         SAVF = array containing current value of f(t,y).
C
C            P = real array of length N.
C
C         WGHT = array of length N containing scale factors.
C                1/WGHT(i) are the diagonal elements of the matrix D.
C
C           WK = work array of length N.
C
C      On return
C
C
C            W = array of length N containing desired
C                matrix-vector product.
C
C In addition, this routine uses the Common variables TN, N, NFE.
C-----------------------------------------------------------------------
      INTEGER I
      DOUBLE PRECISION FAC, PNRM, RPNRM, DVNORM
C
      PNRM = DVNORM (N, P, WGHT)
      RPNRM = 1.0D0/PNRM
      CALL DCOPY (N, Y, 1, W, 1)
      DO 20 I = 1,N
 20     Y(I) = W(I) + P(I)*RPNRM
      CALL F (NEQ, TN, Y, WK)
      NFE = NFE + 1
      CALL DCOPY (N, W, 1, Y, 1)
      FAC = HL0*PNRM
      DO 40 I = 1,N
 40     W(I) = P(I) - FAC*(WK(I) - SAVF(I))
      RETURN
C----------------------- End of Subroutine DATP ------------------------
      END
*DECK DUSOL
      SUBROUTINE DUSOL (NEQ, TN, Y, SAVF, B, WGHT, N, DELTA, HL0, MNEWT,
     1   PSOL, NPSL, X, WP, IWP, WK, IFLAG)
      EXTERNAL PSOL
      INTEGER NEQ, N, MNEWT, NPSL, IWP, IFLAG
      DOUBLE PRECISION TN, Y, SAVF, B, WGHT, DELTA, HL0, X, WP, WK
      DIMENSION NEQ(*), Y(*), SAVF(*), B(*), WGHT(*), X(*),
     1   WP(*), IWP(*), WK(*)
C-----------------------------------------------------------------------
C This routine solves the linear system A * x = b using only a call
C to the user-supplied routine PSOL (no Krylov iteration).
C If the norm of the right-hand side vector b is smaller than DELTA,
C the vector X returned is X = b (if MNEWT = 0) or X = 0 otherwise.
C PSOL is called with an LR argument of 0.
C-----------------------------------------------------------------------
C
C      On entry
C
C          NEQ = problem size, passed to F and PSOL (NEQ(1) = N).
C
C           TN = current value of t.
C
C            Y = array containing current dependent variable vector.
C
C         SAVF = array containing current value of f(t,y).
C
C            B = the right hand side of the system A*x = b.
C
C         WGHT = the vector of length N containing the nonzero
C                elements of the diagonal scaling matrix.
C
C            N = the order of the matrix A, and the lengths
C                of the vectors WGHT, B and X.
C
C        DELTA = tolerance on residuals b - A*x in weighted RMS-norm.
C
C          HL0 = current value of (step size h) * (coefficient l0).
C
C        MNEWT = Newton iteration counter (.ge. 0).
C
C           WK = real work array used by PSOL.
C
C           WP = real work array used by preconditioner PSOL.
C
C          IWP = integer work array used by preconditioner PSOL.
C
C      On return
C
C         X    = the final computed approximation to the solution
C                of the system A*x = b.
C
C         NPSL = the number of calls to PSOL.
C
C        IFLAG = integer error flag:
C                0 means no trouble occurred.
C                3 means there was a recoverable error in PSOL
C                  caused by the preconditioner being out of date.
C               -1 means there was a nonrecoverable error in PSOL.
C
C-----------------------------------------------------------------------
      INTEGER I, IER
      DOUBLE PRECISION BNRM, DVNORM
C
      IFLAG = 0
      NPSL = 0
C-----------------------------------------------------------------------
C Test for an immediate return with X = 0 or X = b.
C-----------------------------------------------------------------------
      BNRM = DVNORM (N, B, WGHT)
      IF (BNRM .GT. DELTA) GO TO 30
      IF (MNEWT .GT. 0) GO TO 10
      CALL DCOPY (N, B, 1, X, 1)
      RETURN
 10   DO 20 I = 1,N
 20     X(I) = 0.0D0
      RETURN
C Make call to PSOL and copy result from B to X. -----------------------
 30   IER = 0
      CALL PSOL (NEQ, TN, Y, SAVF, WK, HL0, WP, IWP, B, 0, IER)
      NPSL = 1
      IF (IER .NE. 0) GO TO 100
      CALL DCOPY (N, B, 1, X, 1)
      RETURN
C-----------------------------------------------------------------------
C This block handles error returns forced by routine PSOL.
C-----------------------------------------------------------------------
 100  CONTINUE
      IF (IER .LT. 0) IFLAG = -1
      IF (IER .GT. 0) IFLAG = 3
      RETURN
C----------------------- End of Subroutine DUSOL -----------------------
      END
*DECK DSRCPK
      SUBROUTINE DSRCPK (RSAV, ISAV, JOB)
C-----------------------------------------------------------------------
C This routine saves or restores (depending on JOB) the contents of
C the Common blocks DLS001, DLPK01, which are used
C internally by the DLSODPK solver.
C
C RSAV = real array of length 13 or more.
C ISAV = integer array of length 38 or more.
C JOB  = flag indicating to save or restore the Common blocks:
C        JOB  = 1 if Common is to be saved (written to RSAV/ISAV)
C        JOB  = 2 if Common is to be restored (read from RSAV/ISAV)
C        A call with JOB = 2 presumes a prior call with JOB = 1.
C-----------------------------------------------------------------------
      INTEGER ISAV, JOB
      INTEGER ILS, ILSP
      INTEGER I, LENILP, LENRLP, LENILS, LENRLS
      DOUBLE PRECISION RSAV,   RLS, RLSP
      DIMENSION RSAV(*), ISAV(*)
      COMMON /DLS001/ RLS(9), ILS(25)
      COMMON /DLPK01/ RLSP(4), ILSP(13)
      DATA LENRLS/9/, LENILS/25/, LENRLP/4/, LENILP/13/
C
      IF (JOB .EQ. 2) GO TO 100
      CALL DCOPY (LENRLS, RLS, 1, RSAV, 1)
      CALL DCOPY (LENRLP, RLSP, 1, RSAV(LENRLS+1), 1)
      DO 20 I = 1,LENILS
 20     ISAV(I) = ILS(I)
      DO 40 I = 1,LENILP
 40     ISAV(LENILS+I) = ILSP(I)
      RETURN
C
 100  CONTINUE
      CALL DCOPY (LENRLS, RSAV, 1, RLS, 1)
      CALL DCOPY (LENRLP, RSAV(LENRLS+1), 1, RLSP, 1)
      DO 120 I = 1,LENILS
 120     ILS(I) = ISAV(I)
      DO 140 I = 1,LENILP
 140     ILSP(I) = ISAV(LENILS+I)
      RETURN
C----------------------- End of Subroutine DSRCPK ----------------------
      END
*DECK DHEFA
      SUBROUTINE DHEFA (A, LDA, N, IPVT, INFO, JOB)
      INTEGER LDA, N, IPVT(*), INFO, JOB
      DOUBLE PRECISION A(LDA,*)
C-----------------------------------------------------------------------
C     This routine is a modification of the LINPACK routine DGEFA and
C     performs an LU decomposition of an upper Hessenberg matrix A.
C     There are two options available:
C
C          (1)  performing a fresh factorization
C          (2)  updating the LU factors by adding a row and a
C               column to the matrix A.
C-----------------------------------------------------------------------
C     DHEFA factors an upper Hessenberg matrix by elimination.
C
C     On entry
C
C        A       DOUBLE PRECISION(LDA, N)
C                the matrix to be factored.
C
C        LDA     INTEGER
C                the leading dimension of the array  A .
C
C        N       INTEGER
C                the order of the matrix  A .
C
C        JOB     INTEGER
C                JOB = 1    means that a fresh factorization of the
C                           matrix A is desired.
C                JOB .ge. 2 means that the current factorization of A
C                           will be updated by the addition of a row
C                           and a column.
C
C     On return
C
C        A       an upper triangular matrix and the multipliers
C                which were used to obtain it.
C                The factorization can be written  A = L*U  where
C                L  is a product of permutation and unit lower
C                triangular matrices and  U  is upper triangular.
C
C        IPVT    INTEGER(N)
C                an integer vector of pivot indices.
C
C        INFO    INTEGER
C                = 0  normal value.
C                = k  if  U(k,k) .eq. 0.0 .  This is not an error
C                     condition for this subroutine, but it does
C                     indicate that DHESL will divide by zero if called.
C
C     Modification of LINPACK, by Peter Brown, LLNL.
C     Written 7/20/83.  This version dated 6/20/01.
C    
C     BLAS called: DAXPY, IDAMAX
C-----------------------------------------------------------------------
      INTEGER IDAMAX, J, K, KM1, KP1, L, NM1
      DOUBLE PRECISION T
C
      IF (JOB .GT. 1) GO TO 80
C
C A new facorization is desired.  This is essentially the LINPACK
C code with the exception that we know there is only one nonzero
C element below the main diagonal.
C
C     Gaussian elimination with partial pivoting
C
      INFO = 0
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 70
      DO 60 K = 1, NM1
         KP1 = K + 1
C
C        Find L = pivot index
C
         L = IDAMAX (2, A(K,K), 1) + K - 1
         IPVT(K) = L
C
C        Zero pivot implies this column already triangularized
C
         IF (A(L,K) .EQ. 0.0D0) GO TO 40
C
C           Interchange if necessary
C
            IF (L .EQ. K) GO TO 10
               T = A(L,K)
               A(L,K) = A(K,K)
               A(K,K) = T
   10       CONTINUE
C
C           Compute multipliers
C
            T = -1.0D0/A(K,K)
            A(K+1,K) = A(K+1,K)*T
C
C           Row elimination with column indexing
C
            DO 30 J = KP1, N
               T = A(L,J)
               IF (L .EQ. K) GO TO 20
                  A(L,J) = A(K,J)
                  A(K,J) = T
   20          CONTINUE
               CALL DAXPY (N-K, T, A(K+1,K), 1, A(K+1,J), 1)
   30       CONTINUE
         GO TO 50
   40    CONTINUE
            INFO = K
   50    CONTINUE
   60 CONTINUE
   70 CONTINUE
      IPVT(N) = N
      IF (A(N,N) .EQ. 0.0D0) INFO = N
      RETURN
C
C The old factorization of A will be updated.  A row and a column
C has been added to the matrix A.
C N-1 is now the old order of the matrix.
C
  80  CONTINUE
      NM1 = N - 1
C
C Perform row interchanges on the elements of the new column, and
C perform elimination operations on the elements using the multipliers.
C
      IF (NM1 .LE. 1) GO TO 105
      DO 100 K = 2,NM1
        KM1 = K - 1
        L = IPVT(KM1)
        T = A(L,N)
        IF (L .EQ. KM1) GO TO 90
          A(L,N) = A(KM1,N)
          A(KM1,N) = T
  90    CONTINUE
        A(K,N) = A(K,N) + A(K,KM1)*T
 100    CONTINUE
 105  CONTINUE
C
C Complete update of factorization by decomposing last 2x2 block.
C
      INFO = 0
C
C        Find L = pivot index
C
         L = IDAMAX (2, A(NM1,NM1), 1) + NM1 - 1
         IPVT(NM1) = L
C
C        Zero pivot implies this column already triangularized
C
         IF (A(L,NM1) .EQ. 0.0D0) GO TO 140
C
C           Interchange if necessary
C
            IF (L .EQ. NM1) GO TO 110
               T = A(L,NM1)
               A(L,NM1) = A(NM1,NM1)
               A(NM1,NM1) = T
  110       CONTINUE
C
C           Compute multipliers
C
            T = -1.0D0/A(NM1,NM1)
            A(N,NM1) = A(N,NM1)*T
C
C           Row elimination with column indexing
C
               T = A(L,N)
               IF (L .EQ. NM1) GO TO 120
                  A(L,N) = A(NM1,N)
                  A(NM1,N) = T
  120          CONTINUE
               A(N,N) = A(N,N) + T*A(N,NM1)
         GO TO 150
  140    CONTINUE
            INFO = NM1
  150    CONTINUE
      IPVT(N) = N
      IF (A(N,N) .EQ. 0.0D0) INFO = N
      RETURN
C----------------------- End of Subroutine DHEFA -----------------------
      END
*DECK DHESL
      SUBROUTINE DHESL (A, LDA, N, IPVT, B)
      INTEGER LDA, N, IPVT(*)
      DOUBLE PRECISION A(LDA,*), B(*)
C-----------------------------------------------------------------------
C This is essentially the LINPACK routine DGESL except for changes
C due to the fact that A is an upper Hessenberg matrix.
C-----------------------------------------------------------------------
C     DHESL solves the real system A * x = b
C     using the factors computed by DHEFA.
C
C     On entry
C
C        A       DOUBLE PRECISION(LDA, N)
C                the output from DHEFA.
C
C        LDA     INTEGER
C                the leading dimension of the array  A .
C
C        N       INTEGER
C                the order of the matrix  A .
C
C        IPVT    INTEGER(N)
C                the pivot vector from DHEFA.
C
C        B       DOUBLE PRECISION(N)
C                the right hand side vector.
C
C     On return
C
C        B       the solution vector  x .
C
C     Modification of LINPACK, by Peter Brown, LLNL.
C     Written 7/20/83.  This version dated 6/20/01.
C
C     BLAS called: DAXPY
C-----------------------------------------------------------------------
      INTEGER K, KB, L, NM1
      DOUBLE PRECISION T
C
      NM1 = N - 1
C
C        Solve  A * x = b
C        First solve  L*y = b
C
         IF (NM1 .LT. 1) GO TO 30
         DO 20 K = 1, NM1
            L = IPVT(K)
            T = B(L)
            IF (L .EQ. K) GO TO 10
               B(L) = B(K)
               B(K) = T
   10       CONTINUE
            B(K+1) = B(K+1) + T*A(K+1,K)
   20    CONTINUE
   30    CONTINUE
C
C        Now solve  U*x = y
C
         DO 40 KB = 1, N
            K = N + 1 - KB
            B(K) = B(K)/A(K,K)
            T = -B(K)
            CALL DAXPY (K-1, T, A(1,K), 1, B(1), 1)
   40    CONTINUE
      RETURN
C----------------------- End of Subroutine DHESL -----------------------
      END
*DECK DHEQR
      SUBROUTINE DHEQR (A, LDA, N, Q, INFO, IJOB)
      INTEGER LDA, N, INFO, IJOB
      DOUBLE PRECISION A(LDA,*), Q(*)
C-----------------------------------------------------------------------
C     This routine performs a QR decomposition of an upper
C     Hessenberg matrix A.  There are two options available:
C
C          (1)  performing a fresh decomposition
C          (2)  updating the QR factors by adding a row and a
C               column to the matrix A.
C-----------------------------------------------------------------------
C     DHEQR decomposes an upper Hessenberg matrix by using Givens
C     rotations.
C
C     On entry
C
C        A       DOUBLE PRECISION(LDA, N)
C                the matrix to be decomposed.
C
C        LDA     INTEGER
C                the leading dimension of the array  A .
C
C        N       INTEGER
C                A is an (N+1) by N Hessenberg matrix.
C
C        IJOB    INTEGER
C                = 1     means that a fresh decomposition of the
C                        matrix A is desired.
C                .ge. 2  means that the current decomposition of A
C                        will be updated by the addition of a row
C                        and a column.
C     On return
C
C        A       the upper triangular matrix R.
C                The factorization can be written Q*A = R, where
C                Q is a product of Givens rotations and R is upper
C                triangular.
C
C        Q       DOUBLE PRECISION(2*N)
C                the factors c and s of each Givens rotation used
C                in decomposing A.
C
C        INFO    INTEGER
C                = 0  normal value.
C                = k  if  A(k,k) .eq. 0.0 .  This is not an error
C                     condition for this subroutine, but it does
C                     indicate that DHELS will divide by zero
C                     if called.
C
C     Modification of LINPACK, by Peter Brown, LLNL.
C     Written 1/13/86.  This version dated 6/20/01.
C-----------------------------------------------------------------------
      INTEGER I, IQ, J, K, KM1, KP1, NM1
      DOUBLE PRECISION C, S, T, T1, T2
C
      IF (IJOB .GT. 1) GO TO 70
C
C A new facorization is desired.
C
C     QR decomposition without pivoting
C
      INFO = 0
      DO 60 K = 1, N
         KM1 = K - 1
         KP1 = K + 1
C
C           Compute kth column of R.
C           First, multiply the kth column of A by the previous
C           k-1 Givens rotations.
C
            IF (KM1 .LT. 1) GO TO 20
            DO 10 J = 1, KM1
              I = 2*(J-1) + 1
              T1 = A(J,K)
              T2 = A(J+1,K)
              C = Q(I)
              S = Q(I+1)
              A(J,K) = C*T1 - S*T2
              A(J+1,K) = S*T1 + C*T2
   10         CONTINUE
C
C           Compute Givens components c and s
C
   20       CONTINUE
            IQ = 2*KM1 + 1
            T1 = A(K,K)
            T2 = A(KP1,K)
            IF (T2 .NE. 0.0D0) GO TO 30
              C = 1.0D0
              S = 0.0D0
              GO TO 50
   30       CONTINUE
            IF (ABS(T2) .LT. ABS(T1)) GO TO 40
              T = T1/T2
              S = -1.0D0/SQRT(1.0D0+T*T)
              C = -S*T
              GO TO 50
   40       CONTINUE
              T = T2/T1
              C = 1.0D0/SQRT(1.0D0+T*T)
              S = -C*T
   50       CONTINUE
            Q(IQ) = C
            Q(IQ+1) = S
            A(K,K) = C*T1 - S*T2
            IF (A(K,K) .EQ. 0.0D0) INFO = K
   60 CONTINUE
      RETURN
C
C The old factorization of A will be updated.  A row and a column
C has been added to the matrix A.
C N by N-1 is now the old size of the matrix.
C
  70  CONTINUE
      NM1 = N - 1
C
C Multiply the new column by the N previous Givens rotations.
C
      DO 100 K = 1,NM1
        I = 2*(K-1) + 1
        T1 = A(K,N)
        T2 = A(K+1,N)
        C = Q(I)
        S = Q(I+1)
        A(K,N) = C*T1 - S*T2
        A(K+1,N) = S*T1 + C*T2
 100    CONTINUE
C
C Complete update of decomposition by forming last Givens rotation,
C and multiplying it times the column vector (A(N,N), A(N+1,N)).
C
      INFO = 0
      T1 = A(N,N)
      T2 = A(N+1,N)
      IF (T2 .NE. 0.0D0) GO TO 110
        C = 1.0D0
        S = 0.0D0
        GO TO 130
 110  CONTINUE
      IF (ABS(T2) .LT. ABS(T1)) GO TO 120
        T = T1/T2
        S = -1.0D0/SQRT(1.0D0+T*T)
        C = -S*T
        GO TO 130
 120  CONTINUE
        T = T2/T1
        C = 1.0D0/SQRT(1.0D0+T*T)
        S = -C*T
 130  CONTINUE
      IQ = 2*N - 1
      Q(IQ) = C
      Q(IQ+1) = S
      A(N,N) = C*T1 - S*T2
      IF (A(N,N) .EQ. 0.0D0) INFO = N
      RETURN
C----------------------- End of Subroutine DHEQR -----------------------
      END
*DECK DHELS
      SUBROUTINE DHELS (A, LDA, N, Q, B)
      INTEGER LDA, N
      DOUBLE PRECISION A(LDA,*), B(*), Q(*)
C-----------------------------------------------------------------------
C This is part of the LINPACK routine DGESL with changes
C due to the fact that A is an upper Hessenberg matrix.
C-----------------------------------------------------------------------
C     DHELS solves the least squares problem
C
C           min (b-A*x, b-A*x)
C
C     using the factors computed by DHEQR.
C
C     On entry
C
C        A       DOUBLE PRECISION(LDA, N)
C                the output from DHEQR which contains the upper
C                triangular factor R in the QR decomposition of A.
C
C        LDA     INTEGER
C                the leading dimension of the array  A .
C
C        N       INTEGER
C                A is originally an (N+1) by N matrix.
C
C        Q       DOUBLE PRECISION(2*N)
C                The coefficients of the N givens rotations
C                used in the QR factorization of A.
C
C        B       DOUBLE PRECISION(N+1)
C                the right hand side vector.
C
C     On return
C
C        B       the solution vector  x .
C
C     Modification of LINPACK, by Peter Brown, LLNL.
C     Written 1/13/86.  This version dated 6/20/01.
C
C     BLAS called: DAXPY
C-----------------------------------------------------------------------
      INTEGER IQ, K, KB, KP1
      DOUBLE PRECISION C, S, T, T1, T2
C
C        Minimize (b-A*x, b-A*x)
C        First form Q*b.
C
         DO 20 K = 1, N
            KP1 = K + 1
            IQ = 2*(K-1) + 1
            C = Q(IQ)
            S = Q(IQ+1)
            T1 = B(K)
            T2 = B(KP1)
            B(K) = C*T1 - S*T2
            B(KP1) = S*T1 + C*T2
   20    CONTINUE
C
C        Now solve  R*x = Q*b.
C
         DO 40 KB = 1, N
            K = N + 1 - KB
            B(K) = B(K)/A(K,K)
            T = -B(K)
            CALL DAXPY (K-1, T, A(1,K), 1, B(1), 1)
   40    CONTINUE
      RETURN
C----------------------- End of Subroutine DHELS -----------------------
      END
*DECK DLHIN
      SUBROUTINE DLHIN (NEQ, N, T0, Y0, YDOT, F, TOUT, UROUND,
     1   EWT, ITOL, ATOL, Y, TEMP, H0, NITER, IER)
      EXTERNAL F
      DOUBLE PRECISION T0, Y0, YDOT, TOUT, UROUND, EWT, ATOL, Y,
     1   TEMP, H0
      INTEGER NEQ, N, ITOL, NITER, IER
      DIMENSION NEQ(*), Y0(*), YDOT(*), EWT(*), ATOL(*), Y(*), TEMP(*)
C-----------------------------------------------------------------------
C Call sequence input -- NEQ, N, T0, Y0, YDOT, F, TOUT, UROUND,
C                        EWT, ITOL, ATOL, Y, TEMP
C Call sequence output -- H0, NITER, IER
C Common block variables accessed -- None
C
C Subroutines called by DLHIN: F, DCOPY
C Function routines called by DLHIN: DVNORM
C-----------------------------------------------------------------------
C This routine computes the step size, H0, to be attempted on the
C first step, when the user has not supplied a value for this.
C
C First we check that TOUT - T0 differs significantly from zero.  Then
C an iteration is done to approximate the initial second derivative
C and this is used to define H from WRMS-norm(H**2 * yddot / 2) = 1.
C A bias factor of 1/2 is applied to the resulting h.
C The sign of H0 is inferred from the initial values of TOUT and T0.
C
C Communication with DLHIN is done with the following variables:
C
C NEQ    = NEQ array of solver, passed to F.
C N      = size of ODE system, input.
C T0     = initial value of independent variable, input.
C Y0     = vector of initial conditions, input.
C YDOT   = vector of initial first derivatives, input.
C F      = name of subroutine for right-hand side f(t,y), input.
C TOUT   = first output value of independent variable
C UROUND = machine unit roundoff
C EWT, ITOL, ATOL = error weights and tolerance parameters
C                   as described in the driver routine, input.
C Y, TEMP = work arrays of length N.
C H0     = step size to be attempted, output.
C NITER  = number of iterations (and of f evaluations) to compute H0,
C          output.
C IER    = the error flag, returned with the value
C          IER = 0  if no trouble occurred, or
C          IER = -1 if TOUT and t0 are considered too close to proceed.
C-----------------------------------------------------------------------
C
C Type declarations for local variables --------------------------------
C
      DOUBLE PRECISION AFI, ATOLI, DELYI, HALF, HG, HLB, HNEW, HRAT,
     1     HUB, HUN, PT1, T1, TDIST, TROUND, TWO, DVNORM, YDDNRM
      INTEGER I, ITER
C-----------------------------------------------------------------------
C The following Fortran-77 declaration is to cause the values of the
C listed (local) variables to be saved between calls to this integrator.
C-----------------------------------------------------------------------
      SAVE HALF, HUN, PT1, TWO
      DATA HALF /0.5D0/, HUN /100.0D0/, PT1 /0.1D0/, TWO /2.0D0/
C
      NITER = 0
      TDIST = ABS(TOUT - T0)
      TROUND = UROUND*MAX(ABS(T0),ABS(TOUT))
      IF (TDIST .LT. TWO*TROUND) GO TO 100
C
C Set a lower bound on H based on the roundoff level in T0 and TOUT. ---
      HLB = HUN*TROUND
C Set an upper bound on H based on TOUT-T0 and the initial Y and YDOT. -
      HUB = PT1*TDIST
      ATOLI = ATOL(1)
      DO 10 I = 1,N
        IF (ITOL .EQ. 2 .OR. ITOL .EQ. 4) ATOLI = ATOL(I)
        DELYI = PT1*ABS(Y0(I)) + ATOLI
        AFI = ABS(YDOT(I))
        IF (AFI*HUB .GT. DELYI) HUB = DELYI/AFI
 10     CONTINUE
C
C Set initial guess for H as geometric mean of upper and lower bounds. -
      ITER = 0
      HG = SQRT(HLB*HUB)
C If the bounds have crossed, exit with the mean value. ----------------
      IF (HUB .LT. HLB) THEN
        H0 = HG
        GO TO 90
      ENDIF
C
C Looping point for iteration. -----------------------------------------
 50   CONTINUE
C Estimate the second derivative as a difference quotient in f. --------
      T1 = T0 + HG
      DO 60 I = 1,N
 60     Y(I) = Y0(I) + HG*YDOT(I)
      CALL F (NEQ, T1, Y, TEMP)
      DO 70 I = 1,N
 70     TEMP(I) = (TEMP(I) - YDOT(I))/HG
      YDDNRM = DVNORM (N, TEMP, EWT)
C Get the corresponding new value of H. --------------------------------
      IF (YDDNRM*HUB*HUB .GT. TWO) THEN
        HNEW = SQRT(TWO/YDDNRM)
      ELSE
        HNEW = SQRT(HG*HUB)
      ENDIF
      ITER = ITER + 1
C-----------------------------------------------------------------------
C Test the stopping conditions.
C Stop if the new and previous H values differ by a factor of .lt. 2.
C Stop if four iterations have been done.  Also, stop with previous H
C if hnew/hg .gt. 2 after first iteration, as this probably means that
C the second derivative value is bad because of cancellation error.
C-----------------------------------------------------------------------
      IF (ITER .GE. 4) GO TO 80
      HRAT = HNEW/HG
      IF ( (HRAT .GT. HALF) .AND. (HRAT .LT. TWO) ) GO TO 80
      IF ( (ITER .GE. 2) .AND. (HNEW .GT. TWO*HG) ) THEN
        HNEW = HG
        GO TO 80
      ENDIF
      HG = HNEW
      GO TO 50
C
C Iteration done.  Apply bounds, bias factor, and sign. ----------------
 80   H0 = HNEW*HALF
      IF (H0 .LT. HLB) H0 = HLB
      IF (H0 .GT. HUB) H0 = HUB
 90   H0 = SIGN(H0, TOUT - T0)
C Restore Y array from Y0, then exit. ----------------------------------
      CALL DCOPY (N, Y0, 1, Y, 1)
      NITER = ITER
      IER = 0
      RETURN
C Error return for TOUT - T0 too small. --------------------------------
 100  IER = -1
      RETURN
C----------------------- End of Subroutine DLHIN -----------------------
      END
*DECK DSTOKA
      SUBROUTINE DSTOKA (NEQ, Y, YH, NYH, YH1, EWT, SAVF, SAVX, ACOR,
     1   WM, IWM, F, JAC, PSOL)
      EXTERNAL F, JAC, PSOL
      INTEGER NEQ, NYH, IWM
      DOUBLE PRECISION Y, YH, YH1, EWT, SAVF, SAVX, ACOR, WM
      DIMENSION NEQ(*), Y(*), YH(NYH,*), YH1(*), EWT(*), SAVF(*),
     1   SAVX(*), ACOR(*), WM(*), IWM(*)
      INTEGER ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, 
     1   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     2   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NPE, NQU
      INTEGER NEWT, NSFI, NSLJ, NJEV
      INTEGER JPRE, JACFLG, LOCWP, LOCIWP, LSAVX, KMP, MAXL, MNEWT,
     1   NNI, NLI, NPS, NCFN, NCFL
      DOUBLE PRECISION CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
      DOUBLE PRECISION STIFR
      DOUBLE PRECISION DELT, EPCON, SQRTN, RSQRTN
      COMMON /DLS001/ CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,
     1   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, 
     2   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     3   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NPE, NQU
      COMMON /DLS002/ STIFR, NEWT, NSFI, NSLJ, NJEV
      COMMON /DLPK01/ DELT, EPCON, SQRTN, RSQRTN,
     1   JPRE, JACFLG, LOCWP, LOCIWP, LSAVX, KMP, MAXL, MNEWT,
     2   NNI, NLI, NPS, NCFN, NCFL
C-----------------------------------------------------------------------
C DSTOKA performs one step of the integration of an initial value
C problem for a system of Ordinary Differential Equations.
C
C This routine was derived from Subroutine DSTODPK in the DLSODPK
C package by the addition of automatic functional/Newton iteration
C switching and logic for re-use of Jacobian data.
C-----------------------------------------------------------------------
C Note: DSTOKA is independent of the value of the iteration method
C indicator MITER, when this is .ne. 0, and hence is independent
C of the type of chord method used, or the Jacobian structure.
C Communication with DSTOKA is done with the following variables:
C
C NEQ    = integer array containing problem size in NEQ(1), and
C          passed as the NEQ argument in all calls to F and JAC.
C Y      = an array of length .ge. N used as the Y argument in
C          all calls to F and JAC.
C YH     = an NYH by LMAX array containing the dependent variables
C          and their approximate scaled derivatives, where
C          LMAX = MAXORD + 1.  YH(i,j+1) contains the approximate
C          j-th derivative of y(i), scaled by H**j/factorial(j)
C          (j = 0,1,...,NQ).  On entry for the first step, the first
C          two columns of YH must be set from the initial values.
C NYH    = a constant integer .ge. N, the first dimension of YH.
C YH1    = a one-dimensional array occupying the same space as YH.
C EWT    = an array of length N containing multiplicative weights
C          for local error measurements.  Local errors in y(i) are
C          compared to 1.0/EWT(i) in various error tests.
C SAVF   = an array of working storage, of length N.
C          Also used for input of YH(*,MAXORD+2) when JSTART = -1
C          and MAXORD .lt. the current order NQ.
C SAVX   = an array of working storage, of length N.
C ACOR   = a work array of length N, used for the accumulated
C          corrections.  On a successful return, ACOR(i) contains
C          the estimated one-step local error in y(i).
C WM,IWM = real and integer work arrays associated with matrix
C          operations in chord iteration (MITER .ne. 0).
C CCMAX  = maximum relative change in H*EL0 before DSETPK is called.
C H      = the step size to be attempted on the next step.
C          H is altered by the error control algorithm during the
C          problem.  H can be either positive or negative, but its
C          sign must remain constant throughout the problem.
C HMIN   = the minimum absolute value of the step size H to be used.
C HMXI   = inverse of the maximum absolute value of H to be used.
C          HMXI = 0.0 is allowed and corresponds to an infinite HMAX.
C          HMIN and HMXI may be changed at any time, but will not
C          take effect until the next change of H is considered.
C TN     = the independent variable. TN is updated on each step taken.
C JSTART = an integer used for input only, with the following
C          values and meanings:
C               0  perform the first step.
C           .gt.0  take a new step continuing from the last.
C              -1  take the next step with a new value of H, MAXORD,
C                    N, METH, MITER, and/or matrix parameters.
C              -2  take the next step with a new value of H,
C                    but with other inputs unchanged.
C          On return, JSTART is set to 1 to facilitate continuation.
C KFLAG  = a completion code with the following meanings:
C               0  the step was succesful.
C              -1  the requested error could not be achieved.
C              -2  corrector convergence could not be achieved.
C              -3  fatal error in DSETPK or DSOLPK.
C          A return with KFLAG = -1 or -2 means either
C          ABS(H) = HMIN or 10 consecutive failures occurred.
C          On a return with KFLAG negative, the values of TN and
C          the YH array are as of the beginning of the last
C          step, and H is the last step size attempted.
C MAXORD = the maximum order of integration method to be allowed.
C MAXCOR = the maximum number of corrector iterations allowed.
C MSBP   = maximum number of steps between DSETPK calls (MITER .gt. 0).
C MXNCF  = maximum number of convergence failures allowed.
C METH/MITER = the method flags.  See description in driver.
C N      = the number of first-order differential equations.
C-----------------------------------------------------------------------
      INTEGER I, I1, IREDO, IRET, J, JB, JOK, M, NCF, NEWQ, NSLOW
      INTEGER IALTH, IPUP, LMAX, MEO, NQNYH, NSLP
      DOUBLE PRECISION DCON, DDN, DEL, DELP, DRC, DSM, DUP, EXDN, EXSM,
     1   EXUP, DFNORM, R, RH, RHDN, RHSM, RHUP, ROC, STIFF, TOLD, DVNORM
      DOUBLE PRECISION CONIT, CRATE, EL(13), ELCO(13,12), HOLD, RMAX,
     1   TESCO(3,12)
      SAVE CONIT, CRATE, EL, ELCO, HOLD, RMAX, TESCO,
     1   IALTH, IPUP, LMAX, MEO, NQNYH, NSLP
C
      KFLAG = 0
      TOLD = TN
      NCF = 0
      IERPJ = 0
      IERSL = 0
      JCUR = 0
      ICF = 0
      DELP = 0.0D0
      IF (JSTART .GT. 0) GO TO 200
      IF (JSTART .EQ. -1) GO TO 100
      IF (JSTART .EQ. -2) GO TO 160
C-----------------------------------------------------------------------
C On the first call, the order is set to 1, and other variables are
C initialized.  RMAX is the maximum ratio by which H can be increased
C in a single step.  It is initially 1.E4 to compensate for the small
C initial H, but then is normally equal to 10.  If a failure
C occurs (in corrector convergence or error test), RMAX is set at 2
C for the next increase.
C-----------------------------------------------------------------------
      LMAX = MAXORD + 1
      NQ = 1
      L = 2
      IALTH = 2
      RMAX = 10000.0D0
      RC = 0.0D0
      EL0 = 1.0D0
      CRATE = 0.7D0
      HOLD = H
      MEO = METH
      NSLP = 0
      NSLJ = 0
      IPUP = 0
      IRET = 3
      NEWT = 0
      STIFR = 0.0D0
      GO TO 140
C-----------------------------------------------------------------------
C The following block handles preliminaries needed when JSTART = -1.
C IPUP is set to MITER to force a matrix update.
C If an order increase is about to be considered (IALTH = 1),
C IALTH is reset to 2 to postpone consideration one more step.
C If the caller has changed METH, DCFODE is called to reset
C the coefficients of the method.
C If the caller has changed MAXORD to a value less than the current
C order NQ, NQ is reduced to MAXORD, and a new H chosen accordingly.
C If H is to be changed, YH must be rescaled.
C If H or METH is being changed, IALTH is reset to L = NQ + 1
C to prevent further changes in H for that many steps.
C-----------------------------------------------------------------------
 100  IPUP = MITER
      LMAX = MAXORD + 1
      IF (IALTH .EQ. 1) IALTH = 2
      IF (METH .EQ. MEO) GO TO 110
      CALL DCFODE (METH, ELCO, TESCO)
      MEO = METH
      IF (NQ .GT. MAXORD) GO TO 120
      IALTH = L
      IRET = 1
      GO TO 150
 110  IF (NQ .LE. MAXORD) GO TO 160
 120  NQ = MAXORD
      L = LMAX
      DO 125 I = 1,L
 125    EL(I) = ELCO(I,NQ)
      NQNYH = NQ*NYH
      RC = RC*EL(1)/EL0
      EL0 = EL(1)
      CONIT = 0.5D0/(NQ+2)
      EPCON = CONIT*TESCO(2,NQ)
      DDN = DVNORM (N, SAVF, EWT)/TESCO(1,L)
      EXDN = 1.0D0/L
      RHDN = 1.0D0/(1.3D0*DDN**EXDN + 0.0000013D0)
      RH = MIN(RHDN,1.0D0)
      IREDO = 3
      IF (H .EQ. HOLD) GO TO 170
      RH = MIN(RH,ABS(H/HOLD))
      H = HOLD
      GO TO 175
C-----------------------------------------------------------------------
C DCFODE is called to get all the integration coefficients for the
C current METH.  Then the EL vector and related constants are reset
C whenever the order NQ is changed, or at the start of the problem.
C-----------------------------------------------------------------------
 140  CALL DCFODE (METH, ELCO, TESCO)
 150  DO 155 I = 1,L
 155    EL(I) = ELCO(I,NQ)
      NQNYH = NQ*NYH
      RC = RC*EL(1)/EL0
      EL0 = EL(1)
      CONIT = 0.5D0/(NQ+2)
      EPCON = CONIT*TESCO(2,NQ)
      GO TO (160, 170, 200), IRET
C-----------------------------------------------------------------------
C If H is being changed, the H ratio RH is checked against
C RMAX, HMIN, and HMXI, and the YH array rescaled.  IALTH is set to
C L = NQ + 1 to prevent a change of H for that many steps, unless
C forced by a convergence or error test failure.
C-----------------------------------------------------------------------
 160  IF (H .EQ. HOLD) GO TO 200
      RH = H/HOLD
      H = HOLD
      IREDO = 3
      GO TO 175
 170  RH = MAX(RH,HMIN/ABS(H))
 175  RH = MIN(RH,RMAX)
      RH = RH/MAX(1.0D0,ABS(H)*HMXI*RH)
      R = 1.0D0
      DO 180 J = 2,L
        R = R*RH
        DO 180 I = 1,N
 180      YH(I,J) = YH(I,J)*R
      H = H*RH
      RC = RC*RH
      IALTH = L
      IF (IREDO .EQ. 0) GO TO 690
C-----------------------------------------------------------------------
C This section computes the predicted values by effectively
C multiplying the YH array by the Pascal triangle matrix.
C The flag IPUP is set according to whether matrix data is involved
C (NEWT .gt. 0 .and. JACFLG .ne. 0) or not, to trigger a call to DSETPK.
C IPUP is set to MITER when RC differs from 1 by more than CCMAX,
C and at least every MSBP steps, when JACFLG = 1.
C RC is the ratio of new to old values of the coefficient  H*EL(1).
C-----------------------------------------------------------------------
 200  IF (NEWT .EQ. 0 .OR. JACFLG .EQ. 0) THEN
        DRC = 0.0D0
        IPUP = 0
        CRATE = 0.7D0
      ELSE
        DRC = ABS(RC - 1.0D0)
        IF (DRC .GT. CCMAX) IPUP = MITER
        IF (NST .GE. NSLP+MSBP) IPUP = MITER
        ENDIF
      TN = TN + H
      I1 = NQNYH + 1
      DO 215 JB = 1,NQ
        I1 = I1 - NYH
CDIR$ IVDEP
        DO 210 I = I1,NQNYH
 210      YH1(I) = YH1(I) + YH1(I+NYH)
 215    CONTINUE
C-----------------------------------------------------------------------
C Up to MAXCOR corrector iterations are taken.  A convergence test is
C made on the RMS-norm of each correction, weighted by the error
C weight vector EWT.  The sum of the corrections is accumulated in the
C vector ACOR(i).  The YH array is not altered in the corrector loop.
C Within the corrector loop, an estimated rate of convergence (ROC)
C and a stiffness ratio estimate (STIFF) are kept.  Corresponding
C global estimates are kept as CRATE and stifr.
C-----------------------------------------------------------------------
 220  M = 0
      MNEWT = 0
      STIFF = 0.0D0
      ROC = 0.05D0
      NSLOW = 0
      DO 230 I = 1,N
 230    Y(I) = YH(I,1)
      CALL F (NEQ, TN, Y, SAVF)
      NFE = NFE + 1
      IF (NEWT .EQ. 0 .OR. IPUP .LE. 0) GO TO 250
C-----------------------------------------------------------------------
C If indicated, DSETPK is called to update any matrix data needed,
C before starting the corrector iteration.
C JOK is set to indicate if the matrix data need not be recomputed.
C IPUP is set to 0 as an indicator that the matrix data is up to date.
C-----------------------------------------------------------------------
      JOK = 1
      IF (NST .EQ. 0 .OR. NST .GT. NSLJ+50) JOK = -1
      IF (ICF .EQ. 1 .AND. DRC .LT. 0.2D0) JOK = -1
      IF (ICF .EQ. 2) JOK = -1
      IF (JOK .EQ. -1) THEN
        NSLJ = NST
        NJEV = NJEV + 1
        ENDIF
      CALL DSETPK (NEQ, Y, YH1, EWT, ACOR, SAVF, JOK, WM, IWM, F, JAC)
      IPUP = 0
      RC = 1.0D0
      DRC = 0.0D0
      NSLP = NST
      CRATE = 0.7D0
      IF (IERPJ .NE. 0) GO TO 430
 250  DO 260 I = 1,N
 260    ACOR(I) = 0.0D0
 270  IF (NEWT .NE. 0) GO TO 350
C-----------------------------------------------------------------------
C In the case of functional iteration, update Y directly from
C the result of the last function evaluation, and STIFF is set to 1.0.
C-----------------------------------------------------------------------
      DO 290 I = 1,N
        SAVF(I) = H*SAVF(I) - YH(I,2)
 290    Y(I) = SAVF(I) - ACOR(I)
      DEL = DVNORM (N, Y, EWT)
      DO 300 I = 1,N
        Y(I) = YH(I,1) + EL(1)*SAVF(I)
 300    ACOR(I) = SAVF(I)
      STIFF = 1.0D0
      GO TO 400
C-----------------------------------------------------------------------
C In the case of the chord method, compute the corrector error,
C and solve the linear system with that as right-hand side and
C P as coefficient matrix.  STIFF is set to the ratio of the norms
C of the residual and the correction vector.
C-----------------------------------------------------------------------
 350  DO 360 I = 1,N
 360    SAVX(I) = H*SAVF(I) - (YH(I,2) + ACOR(I))
      DFNORM = DVNORM (N, SAVX, EWT)
      CALL DSOLPK (NEQ, Y, SAVF, SAVX, EWT, WM, IWM, F, PSOL)
      IF (IERSL .LT. 0) GO TO 430
      IF (IERSL .GT. 0) GO TO 410
      DEL = DVNORM (N, SAVX, EWT)
      IF (DEL .GT. 1.0D-8) STIFF = MAX(STIFF, DFNORM/DEL)
      DO 380 I = 1,N
        ACOR(I) = ACOR(I) + SAVX(I)
 380    Y(I) = YH(I,1) + EL(1)*ACOR(I)
C-----------------------------------------------------------------------
C Test for convergence.  If M .gt. 0, an estimate of the convergence
C rate constant is made for the iteration switch, and is also used
C in the convergence test.   If the iteration seems to be diverging or
C converging at a slow rate (.gt. 0.8 more than once), it is stopped.
C-----------------------------------------------------------------------
 400  IF (M .NE. 0) THEN
        ROC = MAX(0.05D0, DEL/DELP)
        CRATE = MAX(0.2D0*CRATE,ROC)
        ENDIF
      DCON = DEL*MIN(1.0D0,1.5D0*CRATE)/EPCON
      IF (DCON .LE. 1.0D0) GO TO 450
      M = M + 1
      IF (M .EQ. MAXCOR) GO TO 410
      IF (M .GE. 2 .AND. DEL .GT. 2.0D0*DELP) GO TO 410
      IF (ROC .GT. 10.0D0) GO TO 410
      IF (ROC .GT. 0.8D0) NSLOW = NSLOW + 1
      IF (NSLOW .GE. 2) GO TO 410
      MNEWT = M
      DELP = DEL
      CALL F (NEQ, TN, Y, SAVF)
      NFE = NFE + 1
      GO TO 270
C-----------------------------------------------------------------------
C The corrector iteration failed to converge.
C If functional iteration is being done (NEWT = 0) and MITER .gt. 0
C (and this is not the first step), then switch to Newton
C (NEWT = MITER), and retry the step.  (Setting STIFR = 1023 insures
C that a switch back will not occur for 10 step attempts.)
C If Newton iteration is being done, but using a preconditioner that
C is out of date (JACFLG .ne. 0 .and. JCUR = 0), then signal for a
C re-evalutation of the preconditioner, and retry the step.
C In all other cases, the YH array is retracted to its values
C before prediction, and H is reduced, if possible.  If H cannot be
C reduced or MXNCF failures have occurred, exit with KFLAG = -2.
C-----------------------------------------------------------------------
 410  ICF = 1
      IF (NEWT .EQ. 0) THEN
        IF (NST .EQ. 0) GO TO 430
        IF (MITER .EQ. 0) GO TO 430
        NEWT = MITER
        STIFR = 1023.0D0
        IPUP = MITER
        GO TO 220
        ENDIF
      IF (JCUR.EQ.1 .OR. JACFLG.EQ.0) GO TO 430
      IPUP = MITER
      GO TO 220
 430  ICF = 2
      NCF = NCF + 1
      NCFN = NCFN + 1
      RMAX = 2.0D0
      TN = TOLD
      I1 = NQNYH + 1
      DO 445 JB = 1,NQ
        I1 = I1 - NYH
CDIR$ IVDEP
        DO 440 I = I1,NQNYH
 440      YH1(I) = YH1(I) - YH1(I+NYH)
 445    CONTINUE
      IF (IERPJ .LT. 0 .OR. IERSL .LT. 0) GO TO 680
      IF (ABS(H) .LE. HMIN*1.00001D0) GO TO 670
      IF (NCF .EQ. MXNCF) GO TO 670
      RH = 0.5D0
      IPUP = MITER
      IREDO = 1
      GO TO 170
C-----------------------------------------------------------------------
C The corrector has converged.  JCUR is set to 0 to signal that the
C preconditioner involved may need updating later.
C The stiffness ratio STIFR is updated using the latest STIFF value.
C The local error test is made and control passes to statement 500
C if it fails.
C-----------------------------------------------------------------------
 450  JCUR = 0
      IF (NEWT .GT. 0) STIFR = 0.5D0*(STIFR + STIFF)
      IF (M .EQ. 0) DSM = DEL/TESCO(2,NQ)
      IF (M .GT. 0) DSM = DVNORM (N, ACOR, EWT)/TESCO(2,NQ)
      IF (DSM .GT. 1.0D0) GO TO 500
C-----------------------------------------------------------------------
C After a successful step, update the YH array.
C If Newton iteration is being done and STIFR is less than 1.5,
C then switch to functional iteration.
C Consider changing H if IALTH = 1.  Otherwise decrease IALTH by 1.
C If IALTH is then 1 and NQ .lt. MAXORD, then ACOR is saved for
C use in a possible order increase on the next step.
C If a change in H is considered, an increase or decrease in order
C by one is considered also.  A change in H is made only if it is by a
C factor of at least 1.1.  If not, IALTH is set to 3 to prevent
C testing for that many steps.
C-----------------------------------------------------------------------
      KFLAG = 0
      IREDO = 0
      NST = NST + 1
      IF (NEWT .EQ. 0) NSFI = NSFI + 1
      IF (NEWT .GT. 0 .AND. STIFR .LT. 1.5D0) NEWT = 0
      HU = H
      NQU = NQ
      DO 470 J = 1,L
        DO 470 I = 1,N
 470      YH(I,J) = YH(I,J) + EL(J)*ACOR(I)
      IALTH = IALTH - 1
      IF (IALTH .EQ. 0) GO TO 520
      IF (IALTH .GT. 1) GO TO 700
      IF (L .EQ. LMAX) GO TO 700
      DO 490 I = 1,N
 490    YH(I,LMAX) = ACOR(I)
      GO TO 700
C-----------------------------------------------------------------------
C The error test failed.  KFLAG keeps track of multiple failures.
C Restore TN and the YH array to their previous values, and prepare
C to try the step again.  Compute the optimum step size for this or
C one lower order.  After 2 or more failures, H is forced to decrease
C by a factor of 0.2 or less.
C-----------------------------------------------------------------------
 500  KFLAG = KFLAG - 1
      TN = TOLD
      I1 = NQNYH + 1
      DO 515 JB = 1,NQ
        I1 = I1 - NYH
CDIR$ IVDEP
        DO 510 I = I1,NQNYH
 510      YH1(I) = YH1(I) - YH1(I+NYH)
 515    CONTINUE
      RMAX = 2.0D0
      IF (ABS(H) .LE. HMIN*1.00001D0) GO TO 660
      IF (KFLAG .LE. -3) GO TO 640
      IREDO = 2
      RHUP = 0.0D0
      GO TO 540
C-----------------------------------------------------------------------
C Regardless of the success or failure of the step, factors
C RHDN, RHSM, and RHUP are computed, by which H could be multiplied
C at order NQ - 1, order NQ, or order NQ + 1, respectively.
C in the case of failure, RHUP = 0.0 to avoid an order increase.
C the largest of these is determined and the new order chosen
C accordingly.  If the order is to be increased, we compute one
C additional scaled derivative.
C-----------------------------------------------------------------------
 520  RHUP = 0.0D0
      IF (L .EQ. LMAX) GO TO 540
      DO 530 I = 1,N
 530    SAVF(I) = ACOR(I) - YH(I,LMAX)
      DUP = DVNORM (N, SAVF, EWT)/TESCO(3,NQ)
      EXUP = 1.0D0/(L+1)
      RHUP = 1.0D0/(1.4D0*DUP**EXUP + 0.0000014D0)
 540  EXSM = 1.0D0/L
      RHSM = 1.0D0/(1.2D0*DSM**EXSM + 0.0000012D0)
      RHDN = 0.0D0
      IF (NQ .EQ. 1) GO TO 560
      DDN = DVNORM (N, YH(1,L), EWT)/TESCO(1,NQ)
      EXDN = 1.0D0/NQ
      RHDN = 1.0D0/(1.3D0*DDN**EXDN + 0.0000013D0)
 560  IF (RHSM .GE. RHUP) GO TO 570
      IF (RHUP .GT. RHDN) GO TO 590
      GO TO 580
 570  IF (RHSM .LT. RHDN) GO TO 580
      NEWQ = NQ
      RH = RHSM
      GO TO 620
 580  NEWQ = NQ - 1
      RH = RHDN
      IF (KFLAG .LT. 0 .AND. RH .GT. 1.0D0) RH = 1.0D0
      GO TO 620
 590  NEWQ = L
      RH = RHUP
      IF (RH .LT. 1.1D0) GO TO 610
      R = EL(L)/L
      DO 600 I = 1,N
 600    YH(I,NEWQ+1) = ACOR(I)*R
      GO TO 630
 610  IALTH = 3
      GO TO 700
 620  IF ((KFLAG .EQ. 0) .AND. (RH .LT. 1.1D0)) GO TO 610
      IF (KFLAG .LE. -2) RH = MIN(RH,0.2D0)
C-----------------------------------------------------------------------
C If there is a change of order, reset NQ, L, and the coefficients.
C In any case H is reset according to RH and the YH array is rescaled.
C Then exit from 690 if the step was OK, or redo the step otherwise.
C-----------------------------------------------------------------------
      IF (NEWQ .EQ. NQ) GO TO 170
 630  NQ = NEWQ
      L = NQ + 1
      IRET = 2
      GO TO 150
C-----------------------------------------------------------------------
C Control reaches this section if 3 or more failures have occured.
C If 10 failures have occurred, exit with KFLAG = -1.
C It is assumed that the derivatives that have accumulated in the
C YH array have errors of the wrong order.  Hence the first
C derivative is recomputed, and the order is set to 1.  Then
C H is reduced by a factor of 10, and the step is retried,
C until it succeeds or H reaches HMIN.
C-----------------------------------------------------------------------
 640  IF (KFLAG .EQ. -10) GO TO 660
      RH = 0.1D0
      RH = MAX(HMIN/ABS(H),RH)
      H = H*RH
      DO 645 I = 1,N
 645    Y(I) = YH(I,1)
      CALL F (NEQ, TN, Y, SAVF)
      NFE = NFE + 1
      DO 650 I = 1,N
 650    YH(I,2) = H*SAVF(I)
      IPUP = MITER
      IALTH = 5
      IF (NQ .EQ. 1) GO TO 200
      NQ = 1
      L = 2
      IRET = 3
      GO TO 150
C-----------------------------------------------------------------------
C All returns are made through this section.  H is saved in HOLD
C to allow the caller to change H on the next step.
C-----------------------------------------------------------------------
 660  KFLAG = -1
      GO TO 720
 670  KFLAG = -2
      GO TO 720
 680  KFLAG = -3
      GO TO 720
 690  RMAX = 10.0D0
 700  R = 1.0D0/TESCO(2,NQU)
      DO 710 I = 1,N
 710    ACOR(I) = ACOR(I)*R
 720  HOLD = H
      JSTART = 1
      RETURN
C----------------------- End of Subroutine DSTOKA ----------------------
      END
*DECK DSETPK
      SUBROUTINE DSETPK (NEQ, Y, YSV, EWT, FTEM, SAVF, JOK, WM, IWM,
     1                  F, JAC)
      EXTERNAL F, JAC
      INTEGER NEQ, JOK, IWM
      DOUBLE PRECISION Y, YSV, EWT, FTEM, SAVF, WM
      DIMENSION NEQ(*), Y(*), YSV(*), EWT(*), FTEM(*), SAVF(*),
     1   WM(*), IWM(*)
      INTEGER ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, 
     1   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     2   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NPE, NQU
      INTEGER JPRE, JACFLG, LOCWP, LOCIWP, LSAVX, KMP, MAXL, MNEWT,
     1   NNI, NLI, NPS, NCFN, NCFL
      DOUBLE PRECISION CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
      DOUBLE PRECISION DELT, EPCON, SQRTN, RSQRTN
      COMMON /DLS001/ CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,
     1   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, 
     2   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     3   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NPE, NQU
      COMMON /DLPK01/ DELT, EPCON, SQRTN, RSQRTN,
     1   JPRE, JACFLG, LOCWP, LOCIWP, LSAVX, KMP, MAXL, MNEWT,
     2   NNI, NLI, NPS, NCFN, NCFL
C-----------------------------------------------------------------------
C DSETPK is called by DSTOKA to interface with the user-supplied
C routine JAC, to compute and process relevant parts of
C the matrix P = I - H*EL(1)*J , where J is the Jacobian df/dy,
C as need for preconditioning matrix operations later.
C
C In addition to variables described previously, communication
C with DSETPK uses the following:
C Y     = array containing predicted values on entry.
C YSV   = array containing predicted y, to be saved (YH1 in DSTOKA).
C FTEM  = work array of length N (ACOR in DSTOKA).
C SAVF  = array containing f evaluated at predicted y.
C JOK   = input flag showing whether it was judged that Jacobian matrix
C         data need not be recomputed (JOK = 1) or needs to be
C         (JOK = -1).
C WM    = real work space for matrices.
C         Space for preconditioning data starts at WM(LOCWP).
C IWM   = integer work space.
C         Space for preconditioning data starts at IWM(LOCIWP).
C IERPJ = output error flag,  = 0 if no trouble, .gt. 0 if
C         JAC returned an error flag.
C JCUR  = output flag to indicate whether the matrix data involved
C         is now current (JCUR = 1) or not (JCUR = 0).
C This routine also uses Common variables EL0, H, TN, IERPJ, JCUR, NPE.
C-----------------------------------------------------------------------
      INTEGER IER
      DOUBLE PRECISION HL0
C
      IERPJ = 0
      JCUR = 0
      IF (JOK .EQ. -1) JCUR = 1
      HL0 = EL0*H
      CALL JAC (F, NEQ, TN, Y, YSV, EWT, SAVF, FTEM, HL0, JOK,
     1   WM(LOCWP), IWM(LOCIWP), IER)
      NPE = NPE + 1
      IF (IER .EQ. 0) RETURN
      IERPJ = 1
      RETURN
C----------------------- End of Subroutine DSETPK ----------------------
      END
*DECK DSRCKR
      SUBROUTINE DSRCKR (RSAV, ISAV, JOB)
C-----------------------------------------------------------------------
C This routine saves or restores (depending on JOB) the contents of
C the Common blocks DLS001, DLS002, DLSR01, DLPK01, which
C are used internally by the DLSODKR solver.
C
C RSAV = real array of length 17 or more.
C ISAV = integer array of length 46 or more.
C JOB  = flag indicating to save or restore the Common blocks:
C        JOB  = 1 if Common is to be saved (written to RSAV/ISAV)
C        JOB  = 2 if Common is to be restored (read from RSAV/ISAV)
C        A call with JOB = 2 presumes a prior call with JOB = 1.
C-----------------------------------------------------------------------
      INTEGER ISAV, JOB
      INTEGER ILS, ILS2, ILSR, ILSP
      INTEGER I, IOFF, LENILP, LENRLP, LENILS, LENRLS, LENILR, LENRLR
      DOUBLE PRECISION RSAV,   RLS, RLS2, RLSR, RLSP
      DIMENSION RSAV(*), ISAV(*)
      COMMON /DLS001/ RLS(9), ILS(25)
      COMMON /DLS002/ RLS2, ILS2(4)
      COMMON /DLSR01/ RLSR(3), ILSR(4)
      COMMON /DLPK01/ RLSP(4), ILSP(13)
      DATA LENRLS/9/, LENILS/25/, LENRLP/4/, LENILP/13/
      DATA LENRLR/3/, LENILR/4/
C
      IF (JOB .EQ. 2) GO TO 100
      CALL DCOPY (LENRLS, RLS, 1, RSAV, 1)
      RSAV(LENRLS+1) = RLS2
      CALL DCOPY (LENRLR, RLSR, 1, RSAV(LENRLS+2), 1)
      CALL DCOPY (LENRLP, RLSP, 1, RSAV(LENRLS+LENRLR+2), 1)
      DO 20 I = 1,LENILS
 20     ISAV(I) = ILS(I)
      ISAV(LENILS+1) = ILS2(1)
      ISAV(LENILS+2) = ILS2(2)
      ISAV(LENILS+3) = ILS2(3)
      ISAV(LENILS+4) = ILS2(4)
      IOFF = LENILS + 2
      DO 30 I = 1,LENILR
 30     ISAV(IOFF+I) = ILSR(I)
      IOFF = IOFF + LENILR
      DO 40 I = 1,LENILP
 40     ISAV(IOFF+I) = ILSP(I)
      RETURN
C
 100  CONTINUE
      CALL DCOPY (LENRLS, RSAV, 1, RLS, 1)
      RLS2 = RSAV(LENRLS+1)
      CALL DCOPY (LENRLR, RSAV(LENRLS+2), 1, RLSR, 1)
      CALL DCOPY (LENRLP, RSAV(LENRLS+LENRLR+2), 1, RLSP, 1)
      DO 120 I = 1,LENILS
 120    ILS(I) = ISAV(I)
      ILS2(1) = ISAV(LENILS+1)
      ILS2(2) = ISAV(LENILS+2)
      ILS2(3) = ISAV(LENILS+3)
      ILS2(4) = ISAV(LENILS+4)
      IOFF = LENILS + 2
      DO 130 I = 1,LENILR
 130    ILSR(I) = ISAV(IOFF+I)
      IOFF = IOFF + LENILR
      DO 140 I = 1,LENILP
 140    ILSP(I) = ISAV(IOFF+I)
      RETURN
C----------------------- End of Subroutine DSRCKR ----------------------
      END
*DECK DAINVG
      SUBROUTINE DAINVG (RES, ADDA, NEQ, T, Y, YDOT, MITER,
     1                   ML, MU, PW, IPVT, IER )
      EXTERNAL RES, ADDA
      INTEGER NEQ, MITER, ML, MU, IPVT, IER
      INTEGER I, LENPW, MLP1, NROWPW
      DOUBLE PRECISION T, Y, YDOT, PW
      DIMENSION Y(*), YDOT(*), PW(*), IPVT(*)
C-----------------------------------------------------------------------
C This subroutine computes the initial value
C of the vector YDOT satisfying
C     A * YDOT = g(t,y)
C when A is nonsingular.  It is called by DLSODI for
C initialization only, when ISTATE = 0 .
C DAINVG returns an error flag IER:
C   IER  =  0  means DAINVG was successful.
C   IER .ge. 2 means RES returned an error flag IRES = IER.
C   IER .lt. 0 means the a-matrix was found to be singular.
C-----------------------------------------------------------------------
C
      IF (MITER .GE. 4)  GO TO 100
C
C Full matrix case -----------------------------------------------------
C
      LENPW = NEQ*NEQ
      DO 10  I = 1, LENPW
   10    PW(I) = 0.0D0
C
      IER = 1
      CALL RES ( NEQ, T, Y, PW, YDOT, IER )
      IF (IER .GT. 1) RETURN
C
      CALL ADDA ( NEQ, T, Y, 0, 0, PW, NEQ )
      CALL DGEFA ( PW, NEQ, NEQ, IPVT, IER )
      IF (IER .EQ. 0) GO TO 20
         IER = -IER
         RETURN
   20 CALL DGESL ( PW, NEQ, NEQ, IPVT, YDOT, 0 )
      RETURN
C
C Band matrix case -----------------------------------------------------
C
  100 CONTINUE
      NROWPW = 2*ML + MU + 1
      LENPW = NEQ * NROWPW
      DO 110  I = 1, LENPW
  110    PW(I) = 0.0D0
C
      IER = 1
      CALL RES ( NEQ, T, Y, PW, YDOT, IER )
      IF (IER .GT. 1) RETURN
C
      MLP1 = ML + 1
      CALL ADDA ( NEQ, T, Y, ML, MU, PW(MLP1), NROWPW )
      CALL DGBFA ( PW, NROWPW, NEQ, ML, MU, IPVT, IER )
      IF (IER .EQ. 0) GO TO 120
         IER = -IER
         RETURN
  120 CALL DGBSL ( PW, NROWPW, NEQ, ML, MU, IPVT, YDOT, 0 )
      RETURN
C----------------------- End of Subroutine DAINVG ----------------------
      END
*DECK DSTODI
      SUBROUTINE DSTODI (NEQ, Y, YH, NYH, YH1, EWT, SAVF, SAVR,
     1   ACOR, WM, IWM, RES, ADDA, JAC, PJAC, SLVS )
      EXTERNAL RES, ADDA, JAC, PJAC, SLVS
      INTEGER NEQ, NYH, IWM
      INTEGER ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, 
     1   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     2   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NRE, NJE, NQU
      INTEGER I, I1, IREDO, IRES, IRET, J, JB, KGO, M, NCF, NEWQ
      INTEGER IALTH, IPUP, LMAX, MEO, NQNYH, NSLP
      DOUBLE PRECISION Y, YH, YH1, EWT, SAVF, SAVR, ACOR, WM
      DOUBLE PRECISION CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
      DOUBLE PRECISION DCON, DDN, DEL, DELP, DSM, DUP,
     1   ELJH, EL1H, EXDN, EXSM, EXUP,
     2   R, RH, RHDN, RHSM, RHUP, TOLD, DVNORM
      DOUBLE PRECISION CONIT, CRATE, EL(13), ELCO(13,12), HOLD, RMAX,
     1   TESCO(3,12)
      DIMENSION NEQ(*), Y(*), YH(NYH,*), YH1(*), EWT(*), SAVF(*),
     1   SAVR(*), ACOR(*), WM(*), IWM(*)
      COMMON /DLS001/ CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,
     1   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, 
     2   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     3   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NRE, NJE, NQU
      SAVE CONIT, CRATE, EL, ELCO, HOLD, RMAX, TESCO,
     1   IALTH, IPUP, LMAX, MEO, NQNYH, NSLP
C-----------------------------------------------------------------------
C DSTODI performs one step of the integration of an initial value
C problem for a system of Ordinary Differential Equations.
C Note: DSTODI is independent of the value of the iteration method
C indicator MITER, and hence is independent
C of the type of chord method used, or the Jacobian structure.
C Communication with DSTODI is done with the following variables:
C
C NEQ    = integer array containing problem size in NEQ(1), and
C          passed as the NEQ argument in all calls to RES, ADDA,
C          and JAC.
C Y      = an array of length .ge. N used as the Y argument in
C          all calls to RES, JAC, and ADDA.
C NEQ    = integer array containing problem size in NEQ(1), and
C          passed as the NEQ argument in all calls tO RES, G, ADDA,
C          and JAC.
C YH     = an NYH by LMAX array containing the dependent variables
C          and their approximate scaled derivatives, where
C          LMAX = MAXORD + 1.  YH(i,j+1) contains the approximate
C          j-th derivative of y(i), scaled by H**j/factorial(j)
C          (j = 0,1,...,NQ).  On entry for the first step, the first
C          two columns of YH must be set from the initial values.
C NYH    = a constant integer .ge. N, the first dimension of YH.
C YH1    = a one-dimensional array occupying the same space as YH.
C EWT    = an array of length N containing multiplicative weights
C          for local error measurements.  Local errors in y(i) are
C          compared to 1.0/EWT(i) in various error tests.
C SAVF   = an array of working storage, of length N. also used for
C          input of YH(*,MAXORD+2) when JSTART = -1 and MAXORD is less
C          than the current order NQ.
C          Same as YDOTI in the driver.
C SAVR   = an array of working storage, of length N.
C ACOR   = a work array of length N used for the accumulated
C          corrections. On a succesful return, ACOR(i) contains
C          the estimated one-step local error in y(i).
C WM,IWM = real and integer work arrays associated with matrix
C          operations in chord iteration.
C PJAC   = name of routine to evaluate and preprocess Jacobian matrix.
C SLVS   = name of routine to solve linear system in chord iteration.
C CCMAX  = maximum relative change in H*EL0 before PJAC is called.
C H      = the step size to be attempted on the next step.
C          H is altered by the error control algorithm during the
C          problem.  H can be either positive or negative, but its
C          sign must remain constant throughout the problem.
C HMIN   = the minimum absolute value of the step size H to be used.
C HMXI   = inverse of the maximum absolute value of H to be used.
C          HMXI = 0.0 is allowed and corresponds to an infinite HMAX.
C          HMIN and HMXI may be changed at any time, but will not
C          take effect until the next change of H is considered.
C TN     = the independent variable. TN is updated on each step taken.
C JSTART = an integer used for input only, with the following
C          values and meanings:
C               0  perform the first step.
C           .gt.0  take a new step continuing from the last.
C              -1  take the next step with a new value of H, MAXORD,
C                    N, METH, MITER, and/or matrix parameters.
C              -2  take the next step with a new value of H,
C                    but with other inputs unchanged.
C          On return, JSTART is set to 1 to facilitate continuation.
C KFLAG  = a completion code with the following meanings:
C               0  the step was succesful.
C              -1  the requested error could not be achieved.
C              -2  corrector convergence could not be achieved.
C              -3  RES ordered immediate return.
C              -4  error condition from RES could not be avoided.
C              -5  fatal error in PJAC or SLVS.
C          A return with KFLAG = -1, -2, or -4 means either
C          ABS(H) = HMIN or 10 consecutive failures occurred.
C          On a return with KFLAG negative, the values of TN and
C          the YH array are as of the beginning of the last
C          step, and H is the last step size attempted.
C MAXORD = the maximum order of integration method to be allowed.
C MAXCOR = the maximum number of corrector iterations allowed.
C MSBP   = maximum number of steps between PJAC calls.
C MXNCF  = maximum number of convergence failures allowed.
C METH/MITER = the method flags.  See description in driver.
C N      = the number of first-order differential equations.
C-----------------------------------------------------------------------
      KFLAG = 0
      TOLD = TN
      NCF = 0
      IERPJ = 0
      IERSL = 0
      JCUR = 0
      ICF = 0
      DELP = 0.0D0
      IF (JSTART .GT. 0) GO TO 200
      IF (JSTART .EQ. -1) GO TO 100
      IF (JSTART .EQ. -2) GO TO 160
C-----------------------------------------------------------------------
C On the first call, the order is set to 1, and other variables are
C initialized.  RMAX is the maximum ratio by which H can be increased
C in a single step.  It is initially 1.E4 to compensate for the small
C initial H, but then is normally equal to 10.  If a failure
C occurs (in corrector convergence or error test), RMAX is set at 2
C for the next increase.
C-----------------------------------------------------------------------
      LMAX = MAXORD + 1
      NQ = 1
      L = 2
      IALTH = 2
      RMAX = 10000.0D0
      RC = 0.0D0
      EL0 = 1.0D0
      CRATE = 0.7D0
      HOLD = H
      MEO = METH
      NSLP = 0
      IPUP = MITER
      IRET = 3
      GO TO 140
C-----------------------------------------------------------------------
C The following block handles preliminaries needed when JSTART = -1.
C IPUP is set to MITER to force a matrix update.
C If an order increase is about to be considered (IALTH = 1),
C IALTH is reset to 2 to postpone consideration one more step.
C If the caller has changed METH, DCFODE is called to reset
C the coefficients of the method.
C If the caller has changed MAXORD to a value less than the current
C order NQ, NQ is reduced to MAXORD, and a new H chosen accordingly.
C If H is to be changed, YH must be rescaled.
C If H or METH is being changed, IALTH is reset to L = NQ + 1
C to prevent further changes in H for that many steps.
C-----------------------------------------------------------------------
 100  IPUP = MITER
      LMAX = MAXORD + 1
      IF (IALTH .EQ. 1) IALTH = 2
      IF (METH .EQ. MEO) GO TO 110
      CALL DCFODE (METH, ELCO, TESCO)
      MEO = METH
      IF (NQ .GT. MAXORD) GO TO 120
      IALTH = L
      IRET = 1
      GO TO 150
 110  IF (NQ .LE. MAXORD) GO TO 160
 120  NQ = MAXORD
      L = LMAX
      DO 125 I = 1,L
 125    EL(I) = ELCO(I,NQ)
      NQNYH = NQ*NYH
      RC = RC*EL(1)/EL0
      EL0 = EL(1)
      CONIT = 0.5D0/(NQ+2)
      DDN = DVNORM (N, SAVF, EWT)/TESCO(1,L)
      EXDN = 1.0D0/L
      RHDN = 1.0D0/(1.3D0*DDN**EXDN + 0.0000013D0)
      RH = MIN(RHDN,1.0D0)
      IREDO = 3
      IF (H .EQ. HOLD) GO TO 170
      RH = MIN(RH,ABS(H/HOLD))
      H = HOLD
      GO TO 175
C-----------------------------------------------------------------------
C DCFODE is called to get all the integration coefficients for the
C current METH.  Then the EL vector and related constants are reset
C whenever the order NQ is changed, or at the start of the problem.
C-----------------------------------------------------------------------
 140  CALL DCFODE (METH, ELCO, TESCO)
 150  DO 155 I = 1,L
 155    EL(I) = ELCO(I,NQ)
      NQNYH = NQ*NYH
      RC = RC*EL(1)/EL0
      EL0 = EL(1)
      CONIT = 0.5D0/(NQ+2)
      GO TO (160, 170, 200), IRET
C-----------------------------------------------------------------------
C If H is being changed, the H ratio RH is checked against
C RMAX, HMIN, and HMXI, and the YH array rescaled.  IALTH is set to
C L = NQ + 1 to prevent a change of H for that many steps, unless
C forced by a convergence or error test failure.
C-----------------------------------------------------------------------
 160  IF (H .EQ. HOLD) GO TO 200
      RH = H/HOLD
      H = HOLD
      IREDO = 3
      GO TO 175
 170  RH = MAX(RH,HMIN/ABS(H))
 175  RH = MIN(RH,RMAX)
      RH = RH/MAX(1.0D0,ABS(H)*HMXI*RH)
      R = 1.0D0
      DO 180 J = 2,L
        R = R*RH
        DO 180 I = 1,N
 180      YH(I,J) = YH(I,J)*R
      H = H*RH
      RC = RC*RH
      IALTH = L
      IF (IREDO .EQ. 0) GO TO 690
C-----------------------------------------------------------------------
C This section computes the predicted values by effectively
C multiplying the YH array by the Pascal triangle matrix.
C RC is the ratio of new to old values of the coefficient  H*EL(1).
C When RC differs from 1 by more than CCMAX, IPUP is set to MITER
C to force PJAC to be called.
C In any case, PJAC is called at least every MSBP steps.
C-----------------------------------------------------------------------
 200  IF (ABS(RC-1.0D0) .GT. CCMAX) IPUP = MITER
      IF (NST .GE. NSLP+MSBP) IPUP = MITER
      TN = TN + H
      I1 = NQNYH + 1
      DO 215 JB = 1,NQ
        I1 = I1 - NYH
CDIR$ IVDEP
        DO 210 I = I1,NQNYH
 210      YH1(I) = YH1(I) + YH1(I+NYH)
 215    CONTINUE
C-----------------------------------------------------------------------
C Up to MAXCOR corrector iterations are taken.  A convergence test is
C made on the RMS-norm of each correction, weighted by H and the
C error weight vector EWT.  The sum of the corrections is accumulated
C in ACOR(i).  The YH array is not altered in the corrector loop.
C-----------------------------------------------------------------------
 220  M = 0
      DO 230 I = 1,N
        SAVF(I) = YH(I,2) / H
 230    Y(I) = YH(I,1)
      IF (IPUP .LE. 0) GO TO 240
C-----------------------------------------------------------------------
C If indicated, the matrix P = A - H*EL(1)*dr/dy is reevaluated and
C preprocessed before starting the corrector iteration.  IPUP is set
C to 0 as an indicator that this has been done.
C-----------------------------------------------------------------------
      CALL PJAC (NEQ, Y, YH, NYH, EWT, ACOR, SAVR, SAVF, WM, IWM,
     1   RES, JAC, ADDA )
      IPUP = 0
      RC = 1.0D0
      NSLP = NST
      CRATE = 0.7D0
      IF (IERPJ .EQ. 0) GO TO 250
      IF (IERPJ .LT. 0) GO TO 435
      IRES = IERPJ
      GO TO (430, 435, 430), IRES
C Get residual at predicted values, if not already done in PJAC. -------
 240  IRES = 1
      CALL RES ( NEQ, TN, Y, SAVF, SAVR, IRES )
      NRE = NRE + 1
      KGO = ABS(IRES)
      GO TO ( 250, 435, 430 ) , KGO
 250  DO 260 I = 1,N
 260    ACOR(I) = 0.0D0
C-----------------------------------------------------------------------
C Solve the linear system with the current residual as
C right-hand side and P as coefficient matrix.
C-----------------------------------------------------------------------
 270  CONTINUE
      CALL SLVS (WM, IWM, SAVR, SAVF)
      IF (IERSL .LT. 0) GO TO 430
      IF (IERSL .GT. 0) GO TO 410
      EL1H = EL(1) * H
      DEL = DVNORM (N, SAVR, EWT) * ABS(H)
      DO 380 I = 1,N
        ACOR(I) = ACOR(I) + SAVR(I)
        SAVF(I) = ACOR(I) + YH(I,2)/H
 380    Y(I) = YH(I,1) + EL1H*ACOR(I)
C-----------------------------------------------------------------------
C Test for convergence.  If M .gt. 0, an estimate of the convergence
C rate constant is stored in CRATE, and this is used in the test.
C-----------------------------------------------------------------------
      IF (M .NE. 0) CRATE = MAX(0.2D0*CRATE,DEL/DELP)
      DCON = DEL*MIN(1.0D0,1.5D0*CRATE)/(TESCO(2,NQ)*CONIT)
      IF (DCON .LE. 1.0D0) GO TO 460
      M = M + 1
      IF (M .EQ. MAXCOR) GO TO 410
      IF (M .GE. 2 .AND. DEL .GT. 2.0D0*DELP) GO TO 410
      DELP = DEL
      IRES = 1
      CALL RES ( NEQ, TN, Y, SAVF, SAVR, IRES )
      NRE = NRE + 1
      KGO = ABS(IRES)
      GO TO ( 270, 435, 410 ) , KGO
C-----------------------------------------------------------------------
C The correctors failed to converge, or RES has returned abnormally.
C on a convergence failure, if the Jacobian is out of date, PJAC is
C called for the next try.  Otherwise the YH array is retracted to its
C values before prediction, and H is reduced, if possible.
C take an error exit if IRES = 2, or H cannot be reduced, or MXNCF
C failures have occurred, or a fatal error occurred in PJAC or SLVS.
C-----------------------------------------------------------------------
 410  ICF = 1
      IF (JCUR .EQ. 1) GO TO 430
      IPUP = MITER
      GO TO 220
 430  ICF = 2
      NCF = NCF + 1
      RMAX = 2.0D0
 435  TN = TOLD
      I1 = NQNYH + 1
      DO 445 JB = 1,NQ
        I1 = I1 - NYH
CDIR$ IVDEP
        DO 440 I = I1,NQNYH
 440      YH1(I) = YH1(I) - YH1(I+NYH)
 445    CONTINUE
      IF (IRES .EQ. 2) GO TO 680
      IF (IERPJ .LT. 0 .OR. IERSL .LT. 0) GO TO 685
      IF (ABS(H) .LE. HMIN*1.00001D0) GO TO 450
      IF (NCF .EQ. MXNCF) GO TO 450
      RH = 0.25D0
      IPUP = MITER
      IREDO = 1
      GO TO 170
 450  IF (IRES .EQ. 3) GO TO 680
      GO TO 670
C-----------------------------------------------------------------------
C The corrector has converged.  JCUR is set to 0
C to signal that the Jacobian involved may need updating later.
C The local error test is made and control passes to statement 500
C if it fails.
C-----------------------------------------------------------------------
 460  JCUR = 0
      IF (M .EQ. 0) DSM = DEL/TESCO(2,NQ)
      IF (M .GT. 0) DSM = ABS(H) * DVNORM (N, ACOR, EWT)/TESCO(2,NQ)
      IF (DSM .GT. 1.0D0) GO TO 500
C-----------------------------------------------------------------------
C After a successful step, update the YH array.
C Consider changing H if IALTH = 1.  Otherwise decrease IALTH by 1.
C If IALTH is then 1 and NQ .lt. MAXORD, then ACOR is saved for
C use in a possible order increase on the next step.
C If a change in H is considered, an increase or decrease in order
C by one is considered also.  A change in H is made only if it is by a
C factor of at least 1.1.  If not, IALTH is set to 3 to prevent
C testing for that many steps.
C-----------------------------------------------------------------------
      KFLAG = 0
      IREDO = 0
      NST = NST + 1
      HU = H
      NQU = NQ
      DO 470 J = 1,L
        ELJH = EL(J)*H
        DO 470 I = 1,N
 470      YH(I,J) = YH(I,J) + ELJH*ACOR(I)
      IALTH = IALTH - 1
      IF (IALTH .EQ. 0) GO TO 520
      IF (IALTH .GT. 1) GO TO 700
      IF (L .EQ. LMAX) GO TO 700
      DO 490 I = 1,N
 490    YH(I,LMAX) = ACOR(I)
      GO TO 700
C-----------------------------------------------------------------------
C The error test failed.  KFLAG keeps track of multiple failures.
C restore TN and the YH array to their previous values, and prepare
C to try the step again.  Compute the optimum step size for this or
C one lower order.  After 2 or more failures, H is forced to decrease
C by a factor of 0.1 or less.
C-----------------------------------------------------------------------
 500  KFLAG = KFLAG - 1
      TN = TOLD
      I1 = NQNYH + 1
      DO 515 JB = 1,NQ
        I1 = I1 - NYH
CDIR$ IVDEP
        DO 510 I = I1,NQNYH
 510      YH1(I) = YH1(I) - YH1(I+NYH)
 515    CONTINUE
      RMAX = 2.0D0
      IF (ABS(H) .LE. HMIN*1.00001D0) GO TO 660
      IF (KFLAG .LE. -7) GO TO 660
      IREDO = 2
      RHUP = 0.0D0
      GO TO 540
C-----------------------------------------------------------------------
C Regardless of the success or failure of the step, factors
C RHDN, RHSM, and RHUP are computed, by which H could be multiplied
C at order NQ - 1, order NQ, or order NQ + 1, respectively.
C In the case of failure, RHUP = 0.0 to avoid an order increase.
C The largest of these is determined and the new order chosen
C accordingly.  If the order is to be increased, we compute one
C additional scaled derivative.
C-----------------------------------------------------------------------
 520  RHUP = 0.0D0
      IF (L .EQ. LMAX) GO TO 540
      DO 530 I = 1,N
 530    SAVF(I) = ACOR(I) - YH(I,LMAX)
      DUP = ABS(H) * DVNORM (N, SAVF, EWT)/TESCO(3,NQ)
      EXUP = 1.0D0/(L+1)
      RHUP = 1.0D0/(1.4D0*DUP**EXUP + 0.0000014D0)
 540  EXSM = 1.0D0/L
      RHSM = 1.0D0/(1.2D0*DSM**EXSM + 0.0000012D0)
      RHDN = 0.0D0
      IF (NQ .EQ. 1) GO TO 560
      DDN = DVNORM (N, YH(1,L), EWT)/TESCO(1,NQ)
      EXDN = 1.0D0/NQ
      RHDN = 1.0D0/(1.3D0*DDN**EXDN + 0.0000013D0)
 560  IF (RHSM .GE. RHUP) GO TO 570
      IF (RHUP .GT. RHDN) GO TO 590
      GO TO 580
 570  IF (RHSM .LT. RHDN) GO TO 580
      NEWQ = NQ
      RH = RHSM
      GO TO 620
 580  NEWQ = NQ - 1
      RH = RHDN
      IF (KFLAG .LT. 0 .AND. RH .GT. 1.0D0) RH = 1.0D0
      GO TO 620
 590  NEWQ = L
      RH = RHUP
      IF (RH .LT. 1.1D0) GO TO 610
      R = H*EL(L)/L
      DO 600 I = 1,N
 600    YH(I,NEWQ+1) = ACOR(I)*R
      GO TO 630
 610  IALTH = 3
      GO TO 700
 620  IF ((KFLAG .EQ. 0) .AND. (RH .LT. 1.1D0)) GO TO 610
      IF (KFLAG .LE. -2) RH = MIN(RH,0.1D0)
C-----------------------------------------------------------------------
C If there is a change of order, reset NQ, L, and the coefficients.
C In any case H is reset according to RH and the YH array is rescaled.
C Then exit from 690 if the step was OK, or redo the step otherwise.
C-----------------------------------------------------------------------
      IF (NEWQ .EQ. NQ) GO TO 170
 630  NQ = NEWQ
      L = NQ + 1
      IRET = 2
      GO TO 150
C-----------------------------------------------------------------------
C All returns are made through this section.  H is saved in HOLD
C to allow the caller to change H on the next step.
C-----------------------------------------------------------------------
 660  KFLAG = -1
      GO TO 720
 670  KFLAG = -2
      GO TO 720
 680  KFLAG = -1 - IRES
      GO TO 720
 685  KFLAG = -5
      GO TO 720
 690  RMAX = 10.0D0
 700  R = H/TESCO(2,NQU)
      DO 710 I = 1,N
 710    ACOR(I) = ACOR(I)*R
 720  HOLD = H
      JSTART = 1
      RETURN
C----------------------- End of Subroutine DSTODI ----------------------
      END
*DECK DPREPJI
      SUBROUTINE DPREPJI (NEQ, Y, YH, NYH, EWT, RTEM, SAVR, S, WM, IWM,
     1   RES, JAC, ADDA)
      EXTERNAL RES, JAC, ADDA
      INTEGER NEQ, NYH, IWM
      INTEGER ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, 
     1   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     2   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NRE, NJE, NQU
      INTEGER I, I1, I2, IER, II, IRES, J, J1, JJ, LENP,
     1   MBA, MBAND, MEB1, MEBAND, ML, ML3, MU
      DOUBLE PRECISION Y, YH, EWT, RTEM, SAVR, S, WM
      DOUBLE PRECISION CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
      DOUBLE PRECISION CON, FAC, HL0, R, SRUR, YI, YJ, YJJ
      DIMENSION NEQ(*), Y(*), YH(NYH,*), EWT(*), RTEM(*),
     1   S(*), SAVR(*), WM(*), IWM(*)
      COMMON /DLS001/ CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,
     1   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, 
     2   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     3   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NRE, NJE, NQU
C-----------------------------------------------------------------------
C DPREPJI is called by DSTODI to compute and process the matrix
C P = A - H*EL(1)*J , where J is an approximation to the Jacobian dr/dy,
C where r = g(t,y) - A(t,y)*s.  Here J is computed by the user-supplied
C routine JAC if MITER = 1 or 4, or by finite differencing if MITER =
C 2 or 5.  J is stored in WM, rescaled, and ADDA is called to generate
C P. P is then subjected to LU decomposition in preparation
C for later solution of linear systems with P as coefficient
C matrix.  This is done by DGEFA if MITER = 1 or 2, and by
C DGBFA if MITER = 4 or 5.
C
C In addition to variables described previously, communication
C with DPREPJI uses the following:
C Y     = array containing predicted values on entry.
C RTEM  = work array of length N (ACOR in DSTODI).
C SAVR  = array used for output only.  On output it contains the
C         residual evaluated at current values of t and y.
C S     = array containing predicted values of dy/dt (SAVF in DSTODI).
C WM    = real work space for matrices.  On output it contains the
C         LU decomposition of P.
C         Storage of matrix elements starts at WM(3).
C         WM also contains the following matrix-related data:
C         WM(1) = SQRT(UROUND), used in numerical Jacobian increments.
C IWM   = integer work space containing pivot information, starting at
C         IWM(21).  IWM also contains the band parameters
C         ML = IWM(1) and MU = IWM(2) if MITER is 4 or 5.
C EL0   = el(1) (input).
C IERPJ = output error flag.
C         = 0 if no trouble occurred,
C         = 1 if the P matrix was found to be singular,
C         = IRES (= 2 or 3) if RES returned IRES = 2 or 3.
C JCUR  = output flag = 1 to indicate that the Jacobian matrix
C         (or approximation) is now current.
C This routine also uses the Common variables EL0, H, TN, UROUND,
C MITER, N, NRE, and NJE.
C-----------------------------------------------------------------------
      NJE = NJE + 1
      HL0 = H*EL0
      IERPJ = 0
      JCUR = 1
      GO TO (100, 200, 300, 400, 500), MITER
C If MITER = 1, call RES, then JAC, and multiply by scalar. ------------
 100  IRES = 1
      CALL RES (NEQ, TN, Y, S, SAVR, IRES)
      NRE = NRE + 1
      IF (IRES .GT. 1) GO TO 600
      LENP = N*N
      DO 110 I = 1,LENP
 110    WM(I+2) = 0.0D0
      CALL JAC ( NEQ, TN, Y, S, 0, 0, WM(3), N )
      CON = -HL0
      DO 120 I = 1,LENP
 120    WM(I+2) = WM(I+2)*CON
      GO TO 240
C If MITER = 2, make N + 1 calls to RES to approximate J. --------------
 200  CONTINUE
      IRES = -1
      CALL RES (NEQ, TN, Y, S, SAVR, IRES)
      NRE = NRE + 1
      IF (IRES .GT. 1) GO TO 600
      SRUR = WM(1)
      J1 = 2
      DO 230 J = 1,N
        YJ = Y(J)
        R = MAX(SRUR*ABS(YJ),0.01D0/EWT(J))
        Y(J) = Y(J) + R
        FAC = -HL0/R
        CALL RES ( NEQ, TN, Y, S, RTEM, IRES )
        NRE = NRE + 1
        IF (IRES .GT. 1) GO TO 600
        DO 220 I = 1,N
 220      WM(I+J1) = (RTEM(I) - SAVR(I))*FAC
        Y(J) = YJ
        J1 = J1 + N
 230    CONTINUE
      IRES = 1
      CALL RES (NEQ, TN, Y, S, SAVR, IRES)
      NRE = NRE + 1
      IF (IRES .GT. 1) GO TO 600
C Add matrix A. --------------------------------------------------------
 240  CONTINUE
      CALL ADDA(NEQ, TN, Y, 0, 0, WM(3), N)
C Do LU decomposition on P. --------------------------------------------
      CALL DGEFA (WM(3), N, N, IWM(21), IER)
      IF (IER .NE. 0) IERPJ = 1
      RETURN
C Dummy section for MITER = 3
 300  RETURN
C If MITER = 4, call RES, then JAC, and multiply by scalar. ------------
 400  IRES = 1
      CALL RES (NEQ, TN, Y, S, SAVR, IRES)
      NRE = NRE + 1
      IF (IRES .GT. 1) GO TO 600
      ML = IWM(1)
      MU = IWM(2)
      ML3 = ML + 3
      MBAND = ML + MU + 1
      MEBAND = MBAND + ML
      LENP = MEBAND*N
      DO 410 I = 1,LENP
 410    WM(I+2) = 0.0D0
      CALL JAC ( NEQ, TN, Y, S, ML, MU, WM(ML3), MEBAND)
      CON = -HL0
      DO 420 I = 1,LENP
 420    WM(I+2) = WM(I+2)*CON
      GO TO 570
C If MITER = 5, make ML + MU + 2 calls to RES to approximate J. --------
 500  CONTINUE
      IRES = -1
      CALL RES (NEQ, TN, Y, S, SAVR, IRES)
      NRE = NRE + 1
      IF (IRES .GT. 1) GO TO 600
      ML = IWM(1)
      MU = IWM(2)
      ML3 = ML + 3
      MBAND = ML + MU + 1
      MBA = MIN(MBAND,N)
      MEBAND = MBAND + ML
      MEB1 = MEBAND - 1
      SRUR = WM(1)
      DO 560 J = 1,MBA
        DO 530 I = J,N,MBAND
          YI = Y(I)
          R = MAX(SRUR*ABS(YI),0.01D0/EWT(I))
 530      Y(I) = Y(I) + R
        CALL RES ( NEQ, TN, Y, S, RTEM, IRES)
        NRE = NRE + 1
        IF (IRES .GT. 1) GO TO 600
        DO 550 JJ = J,N,MBAND
          Y(JJ) = YH(JJ,1)
          YJJ = Y(JJ)
          R = MAX(SRUR*ABS(YJJ),0.01D0/EWT(JJ))
          FAC = -HL0/R
          I1 = MAX(JJ-MU,1)
          I2 = MIN(JJ+ML,N)
          II = JJ*MEB1 - ML + 2
          DO 540 I = I1,I2
 540        WM(II+I) = (RTEM(I) - SAVR(I))*FAC
 550      CONTINUE
 560    CONTINUE
      IRES = 1
      CALL RES (NEQ, TN, Y, S, SAVR, IRES)
      NRE = NRE + 1
      IF (IRES .GT. 1) GO TO 600
C Add matrix A. --------------------------------------------------------
  570 CONTINUE
      CALL ADDA(NEQ, TN, Y, ML, MU, WM(ML3), MEBAND)
C Do LU decomposition of P. --------------------------------------------
      CALL DGBFA (WM(3), MEBAND, N, ML, MU, IWM(21), IER)
      IF (IER .NE. 0) IERPJ = 1
      RETURN
C Error return for IRES = 2 or IRES = 3 return from RES. ---------------
 600  IERPJ = IRES
      RETURN
C----------------------- End of Subroutine DPREPJI ---------------------
      END
*DECK DAIGBT
      SUBROUTINE DAIGBT (RES, ADDA, NEQ, T, Y, YDOT,
     1                   MB, NB, PW, IPVT, IER )
      EXTERNAL RES, ADDA
      INTEGER NEQ, MB, NB, IPVT, IER
      INTEGER I, LENPW, LBLOX, LPB, LPC
      DOUBLE PRECISION T, Y, YDOT, PW
      DIMENSION Y(*), YDOT(*), PW(*), IPVT(*), NEQ(*)
C-----------------------------------------------------------------------
C This subroutine computes the initial value
C of the vector YDOT satisfying
C     A * YDOT = g(t,y)
C when A is nonsingular.  It is called by DLSOIBT for
C initialization only, when ISTATE = 0 .
C DAIGBT returns an error flag IER:
C   IER  =  0  means DAIGBT was successful.
C   IER .ge. 2 means RES returned an error flag IRES = IER.
C   IER .lt. 0 means the A matrix was found to have a singular
C              diagonal block (hence YDOT could not be solved for).
C-----------------------------------------------------------------------
      LBLOX = MB*MB*NB
      LPB = 1 + LBLOX
      LPC = LPB + LBLOX
      LENPW = 3*LBLOX
      DO 10 I = 1,LENPW
 10     PW(I) = 0.0D0
      IER = 1
      CALL RES (NEQ, T, Y, PW, YDOT, IER)
      IF (IER .GT. 1) RETURN
      CALL ADDA (NEQ, T, Y, MB, NB, PW(1), PW(LPB), PW(LPC) )
      CALL DDECBT (MB, NB, PW, PW(LPB), PW(LPC), IPVT, IER)
      IF (IER .EQ. 0) GO TO 20
      IER = -IER
      RETURN
 20   CALL DSOLBT (MB, NB, PW, PW(LPB), PW(LPC), YDOT, IPVT)
      RETURN
C----------------------- End of Subroutine DAIGBT ----------------------
      END
*DECK DPJIBT
      SUBROUTINE DPJIBT (NEQ, Y, YH, NYH, EWT, RTEM, SAVR, S, WM, IWM,
     1   RES, JAC, ADDA)
      EXTERNAL RES, JAC, ADDA
      INTEGER NEQ, NYH, IWM
      INTEGER ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, 
     1   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     2   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NRE, NJE, NQU
      INTEGER I, IER, IIA, IIB, IIC, IPA, IPB, IPC, IRES, J, J1, J2,
     1   K, K1, LENP, LBLOX, LPB, LPC, MB, MBSQ, MWID, NB
      DOUBLE PRECISION Y, YH, EWT, RTEM, SAVR, S, WM
      DOUBLE PRECISION CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
      DOUBLE PRECISION CON, FAC, HL0, R, SRUR
      DIMENSION NEQ(*), Y(*), YH(NYH,*), EWT(*), RTEM(*),
     1   S(*), SAVR(*), WM(*), IWM(*)
      COMMON /DLS001/ CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,
     1   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, 
     2   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     3   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NRE, NJE, NQU
C-----------------------------------------------------------------------
C DPJIBT is called by DSTODI to compute and process the matrix
C P = A - H*EL(1)*J , where J is an approximation to the Jacobian dr/dy,
C and r = g(t,y) - A(t,y)*s.  Here J is computed by the user-supplied
C routine JAC if MITER = 1, or by finite differencing if MITER = 2.
C J is stored in WM, rescaled, and ADDA is called to generate P.
C P is then subjected to LU decomposition by DDECBT in preparation
C for later solution of linear systems with P as coefficient matrix.
C
C In addition to variables described previously, communication
C with DPJIBT uses the following:
C Y     = array containing predicted values on entry.
C RTEM  = work array of length N (ACOR in DSTODI).
C SAVR  = array used for output only.  On output it contains the
C         residual evaluated at current values of t and y.
C S     = array containing predicted values of dy/dt (SAVF in DSTODI).
C WM    = real work space for matrices.  On output it contains the
C         LU decomposition of P.
C         Storage of matrix elements starts at WM(3).
C         WM also contains the following matrix-related data:
C         WM(1) = SQRT(UROUND), used in numerical Jacobian increments.
C IWM   = integer work space containing pivot information, starting at
C         IWM(21).  IWM also contains block structure parameters
C         MB = IWM(1) and NB = IWM(2).
C EL0   = EL(1) (input).
C IERPJ = output error flag.
C         = 0 if no trouble occurred,
C         = 1 if the P matrix was found to be unfactorable,
C         = IRES (= 2 or 3) if RES returned IRES = 2 or 3.
C JCUR  = output flag = 1 to indicate that the Jacobian matrix
C         (or approximation) is now current.
C This routine also uses the Common variables EL0, H, TN, UROUND,
C MITER, N, NRE, and NJE.
C-----------------------------------------------------------------------
      NJE = NJE + 1
      HL0 = H*EL0
      IERPJ = 0
      JCUR = 1
      MB = IWM(1)
      NB = IWM(2)
      MBSQ = MB*MB
      LBLOX = MBSQ*NB
      LPB = 3 + LBLOX
      LPC = LPB + LBLOX
      LENP = 3*LBLOX
      GO TO (100, 200), MITER
C If MITER = 1, call RES, then JAC, and multiply by scalar. ------------
 100  IRES = 1
      CALL RES (NEQ, TN, Y, S, SAVR, IRES)
      NRE = NRE + 1
      IF (IRES .GT. 1) GO TO 600
      DO 110 I = 1,LENP
 110    WM(I+2) = 0.0D0
      CALL JAC (NEQ, TN, Y, S, MB, NB, WM(3), WM(LPB), WM(LPC))
      CON = -HL0
      DO 120 I = 1,LENP
 120    WM(I+2) = WM(I+2)*CON
      GO TO 260
C
C If MITER = 2, make 3*MB + 1 calls to RES to approximate J. -----------
 200  CONTINUE
      IRES = -1
      CALL RES (NEQ, TN, Y, S, SAVR, IRES)
      NRE = NRE + 1
      IF (IRES .GT. 1) GO TO 600
      MWID = 3*MB
      SRUR = WM(1)
      DO 205 I = 1,LENP
 205    WM(2+I) = 0.0D0
      DO 250 K = 1,3
        DO 240 J = 1,MB
C         Increment Y(I) for group of column indices, and call RES. ----
          J1 = J+(K-1)*MB
          DO 210 I = J1,N,MWID
            R = MAX(SRUR*ABS(Y(I)),0.01D0/EWT(I))
            Y(I) = Y(I) + R
 210      CONTINUE
          CALL RES (NEQ, TN, Y, S, RTEM, IRES)
          NRE = NRE + 1
          IF (IRES .GT. 1) GO TO 600
          DO 215 I = 1,N
 215        RTEM(I) = RTEM(I) - SAVR(I)
          K1 = K
          DO 230 I = J1,N,MWID
C           Get Jacobian elements in column I (block-column K1). -------
            Y(I) = YH(I,1)
            R = MAX(SRUR*ABS(Y(I)),0.01D0/EWT(I))
            FAC = -HL0/R
C           Compute and load elements PA(*,J,K1). ----------------------
            IIA = I - J
            IPA = 2 + (J-1)*MB + (K1-1)*MBSQ
            DO 221 J2 = 1,MB
 221          WM(IPA+J2) = RTEM(IIA+J2)*FAC
            IF (K1 .LE. 1) GO TO 223
C           Compute and load elements PB(*,J,K1-1). --------------------
            IIB = IIA - MB
            IPB = IPA + LBLOX - MBSQ
            DO 222 J2 = 1,MB
 222          WM(IPB+J2) = RTEM(IIB+J2)*FAC
 223        CONTINUE
            IF (K1 .GE. NB) GO TO 225
C           Compute and load elements PC(*,J,K1+1). --------------------
            IIC = IIA + MB
            IPC = IPA + 2*LBLOX + MBSQ
            DO 224 J2 = 1,MB
 224          WM(IPC+J2) = RTEM(IIC+J2)*FAC
 225        CONTINUE
            IF (K1 .NE. 3) GO TO 227
C           Compute and load elements PC(*,J,1). -----------------------
            IPC = IPA - 2*MBSQ + 2*LBLOX
            DO 226 J2 = 1,MB
 226          WM(IPC+J2) = RTEM(J2)*FAC
 227        CONTINUE
            IF (K1 .NE. NB-2) GO TO 229
C           Compute and load elements PB(*,J,NB). ----------------------
            IIB = N - MB
            IPB = IPA + 2*MBSQ + LBLOX
            DO 228 J2 = 1,MB
 228          WM(IPB+J2) = RTEM(IIB+J2)*FAC
 229      K1 = K1 + 3
 230      CONTINUE
 240    CONTINUE
 250  CONTINUE
C RES call for first corrector iteration. ------------------------------
      IRES = 1
      CALL RES (NEQ, TN, Y, S, SAVR, IRES)
      NRE = NRE + 1
      IF (IRES .GT. 1) GO TO 600
C Add matrix A. --------------------------------------------------------
 260  CONTINUE
      CALL ADDA (NEQ, TN, Y, MB, NB, WM(3), WM(LPB), WM(LPC))
C Do LU decomposition on P. --------------------------------------------
      CALL DDECBT (MB, NB, WM(3), WM(LPB), WM(LPC), IWM(21), IER)
      IF (IER .NE. 0) IERPJ = 1
      RETURN
C Error return for IRES = 2 or IRES = 3 return from RES. ---------------
 600  IERPJ = IRES
      RETURN
C----------------------- End of Subroutine DPJIBT ----------------------
      END
*DECK DSLSBT
      SUBROUTINE DSLSBT (WM, IWM, X, TEM)
      INTEGER IWM
      INTEGER LBLOX, LPB, LPC, MB, NB
      DOUBLE PRECISION WM, X, TEM
      DIMENSION WM(*), IWM(*), X(*), TEM(*)
C-----------------------------------------------------------------------
C This routine acts as an interface between the core integrator
C routine and the DSOLBT routine for the solution of the linear system
C arising from chord iteration.
C Communication with DSLSBT uses the following variables:
C WM    = real work space containing the LU decomposition,
C         starting at WM(3).
C IWM   = integer work space containing pivot information, starting at
C         IWM(21).  IWM also contains block structure parameters
C         MB = IWM(1) and NB = IWM(2).
C X     = the right-hand side vector on input, and the solution vector
C         on output, of length N.
C TEM   = vector of work space of length N, not used in this version.
C-----------------------------------------------------------------------
      MB = IWM(1)
      NB = IWM(2)
      LBLOX = MB*MB*NB
      LPB = 3 + LBLOX
      LPC = LPB + LBLOX
      CALL DSOLBT (MB, NB, WM(3), WM(LPB), WM(LPC), X, IWM(21))
      RETURN
C----------------------- End of Subroutine DSLSBT ----------------------
      END
*DECK DDECBT
      SUBROUTINE DDECBT (M, N, A, B, C, IP, IER)
      INTEGER M, N, IP(M,N), IER
      DOUBLE PRECISION A(M,M,N), B(M,M,N), C(M,M,N)
C-----------------------------------------------------------------------
C Block-tridiagonal matrix decomposition routine.
C Written by A. C. Hindmarsh.
C Latest revision:  November 10, 1983 (ACH)
C Reference:  UCID-30150
C             Solution of Block-Tridiagonal Systems of Linear
C             Algebraic Equations
C             A.C. Hindmarsh
C             February 1977
C The input matrix contains three blocks of elements in each block-row,
C including blocks in the (1,3) and (N,N-2) block positions.
C DDECBT uses block Gauss elimination and Subroutines DGEFA and DGESL
C for solution of blocks.  Partial pivoting is done within
C block-rows only.
C
C Note: this version uses LINPACK routines DGEFA/DGESL instead of
C of dec/sol for solution of blocks, and it uses the BLAS routine DDOT
C for dot product calculations.
C
C Input:
C     M = order of each block.
C     N = number of blocks in each direction of the matrix.
C         N must be 4 or more.  The complete matrix has order M*N.
C     A = M by M by N array containing diagonal blocks.
C         A(i,j,k) contains the (i,j) element of the k-th block.
C     B = M by M by N array containing the super-diagonal blocks
C         (in B(*,*,k) for k = 1,...,N-1) and the block in the (N,N-2)
C         block position (in B(*,*,N)).
C     C = M by M by N array containing the subdiagonal blocks
C         (in C(*,*,k) for k = 2,3,...,N) and the block in the
C         (1,3) block position (in C(*,*,1)).
C    IP = integer array of length M*N for working storage.
C Output:
C A,B,C = M by M by N arrays containing the block-LU decomposition
C         of the input matrix.
C    IP = M by N array of pivot information.  IP(*,k) contains
C         information for the k-th digonal block.
C   IER = 0  if no trouble occurred, or
C       = -1 if the input value of M or N was illegal, or
C       = k  if a singular matrix was found in the k-th diagonal block.
C Use DSOLBT to solve the associated linear system.
C
C External routines required: DGEFA and DGESL (from LINPACK) and
C DDOT (from the BLAS, or Basic Linear Algebra package).
C-----------------------------------------------------------------------
      INTEGER NM1, NM2, KM1, I, J, K
      DOUBLE PRECISION DP, DDOT
      IF (M .LT. 1 .OR. N .LT. 4) GO TO 210
      NM1 = N - 1
      NM2 = N - 2
C Process the first block-row. -----------------------------------------
      CALL DGEFA (A, M, M, IP, IER)
      K = 1
      IF (IER .NE. 0) GO TO 200
      DO 10 J = 1,M
        CALL DGESL (A, M, M, IP, B(1,J,1), 0)
        CALL DGESL (A, M, M, IP, C(1,J,1), 0)
 10     CONTINUE
C Adjust B(*,*,2). -----------------------------------------------------
      DO 40 J = 1,M
        DO 30 I = 1,M
          DP = DDOT (M, C(I,1,2), M, C(1,J,1), 1)
          B(I,J,2) = B(I,J,2) - DP
 30       CONTINUE
 40     CONTINUE
C Main loop.  Process block-rows 2 to N-1. -----------------------------
      DO 100 K = 2,NM1
        KM1 = K - 1
        DO 70 J = 1,M
          DO 60 I = 1,M
            DP = DDOT (M, C(I,1,K), M, B(1,J,KM1), 1)
            A(I,J,K) = A(I,J,K) - DP
 60         CONTINUE
 70       CONTINUE
        CALL DGEFA (A(1,1,K), M, M, IP(1,K), IER)
        IF (IER .NE. 0) GO TO 200
        DO 80 J = 1,M
 80       CALL DGESL (A(1,1,K), M, M, IP(1,K), B(1,J,K), 0)
 100    CONTINUE
C Process last block-row and return. -----------------------------------
      DO 130 J = 1,M
        DO 120 I = 1,M
          DP = DDOT (M, B(I,1,N), M, B(1,J,NM2), 1)
          C(I,J,N) = C(I,J,N) - DP
 120      CONTINUE
 130    CONTINUE
      DO 160 J = 1,M
        DO 150 I = 1,M
          DP = DDOT (M, C(I,1,N), M, B(1,J,NM1), 1)
          A(I,J,N) = A(I,J,N) - DP
 150      CONTINUE
 160    CONTINUE
      CALL DGEFA (A(1,1,N), M, M, IP(1,N), IER)
      K = N
      IF (IER .NE. 0) GO TO 200
      RETURN
C Error returns. -------------------------------------------------------
 200  IER = K
      RETURN
 210  IER = -1
      RETURN
C----------------------- End of Subroutine DDECBT ----------------------
      END
*DECK DSOLBT
      SUBROUTINE DSOLBT (M, N, A, B, C, Y, IP)
      INTEGER M, N, IP(M,N)
      DOUBLE PRECISION A(M,M,N), B(M,M,N), C(M,M,N), Y(M,N)
C-----------------------------------------------------------------------
C Solution of block-tridiagonal linear system.
C Coefficient matrix must have been previously processed by DDECBT.
C M, N, A,B,C, and IP  must not have been changed since call to DDECBT.
C Written by A. C. Hindmarsh.
C Input:
C     M = order of each block.
C     N = number of blocks in each direction of matrix.
C A,B,C = M by M by N arrays containing block LU decomposition
C         of coefficient matrix from DDECBT.
C    IP = M by N integer array of pivot information from DDECBT.
C     Y = array of length M*N containg the right-hand side vector
C         (treated as an M by N array here).
C Output:
C     Y = solution vector, of length M*N.
C
C External routines required: DGESL (LINPACK) and DDOT (BLAS).
C-----------------------------------------------------------------------
C
      INTEGER NM1, NM2, I, K, KB, KM1, KP1
      DOUBLE PRECISION DP, DDOT
      NM1 = N - 1
      NM2 = N - 2
C Forward solution sweep. ----------------------------------------------
      CALL DGESL (A, M, M, IP, Y, 0)
      DO 30 K = 2,NM1
        KM1 = K - 1
        DO 20 I = 1,M
          DP = DDOT (M, C(I,1,K), M, Y(1,KM1), 1)
          Y(I,K) = Y(I,K) - DP
 20       CONTINUE
        CALL DGESL (A(1,1,K), M, M, IP(1,K), Y(1,K), 0)
 30     CONTINUE
      DO 50 I = 1,M
        DP = DDOT (M, C(I,1,N), M, Y(1,NM1), 1)
     1     + DDOT (M, B(I,1,N), M, Y(1,NM2), 1)
        Y(I,N) = Y(I,N) - DP
 50     CONTINUE
      CALL DGESL (A(1,1,N), M, M, IP(1,N), Y(1,N), 0)
C Backward solution sweep. ---------------------------------------------
      DO 80 KB = 1,NM1
        K = N - KB
        KP1 = K + 1
        DO 70 I = 1,M
          DP = DDOT (M, B(I,1,K), M, Y(1,KP1), 1)
          Y(I,K) = Y(I,K) - DP
 70       CONTINUE
 80     CONTINUE
      DO 100 I = 1,M
        DP = DDOT (M, C(I,1,1), M, Y(1,3), 1)
        Y(I,1) = Y(I,1) - DP
 100    CONTINUE
      RETURN
C----------------------- End of Subroutine DSOLBT ----------------------
      END
*DECK DIPREPI
      SUBROUTINE DIPREPI (NEQ, Y, S, RWORK, IA, JA, IC, JC, IPFLAG,
     1   RES, JAC, ADDA)
      EXTERNAL RES, JAC, ADDA
      INTEGER NEQ, IA, JA, IC, JC, IPFLAG
      INTEGER ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, 
     1   LYH, LEWT, LACOR, LSAVR, LWM, LIWM, METH, MITER,
     2   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NRE, NJE, NQU
      INTEGER IPLOST, IESP, ISTATC, IYS, IBA, IBIAN, IBJAN, IBJGP,
     1   IPIAN, IPJAN, IPJGP, IPIGP, IPR, IPC, IPIC, IPISP, IPRSP, IPA,
     2   LENYH, LENYHM, LENWK, LREQ, LRAT, LREST, LWMIN, MOSS, MSBJ,
     3   NSLJ, NGP, NLU, NNZ, NSP, NZL, NZU
      INTEGER I, IMAX, LEWTN, LYHD, LYHN
      DOUBLE PRECISION Y, S, RWORK
      DOUBLE PRECISION CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
      DOUBLE PRECISION RLSS
      DIMENSION NEQ(*), Y(*), S(*), RWORK(*), IA(*), JA(*), IC(*), JC(*)
      COMMON /DLS001/ CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,
     1   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, 
     2   LYH, LEWT, LACOR, LSAVR, LWM, LIWM, METH, MITER,
     3   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NRE, NJE, NQU
      COMMON /DLSS01/ RLSS(6),
     1   IPLOST, IESP, ISTATC, IYS, IBA, IBIAN, IBJAN, IBJGP,
     2   IPIAN, IPJAN, IPJGP, IPIGP, IPR, IPC, IPIC, IPISP, IPRSP, IPA,
     3   LENYH, LENYHM, LENWK, LREQ, LRAT, LREST, LWMIN, MOSS, MSBJ,
     4   NSLJ, NGP, NLU, NNZ, NSP, NZL, NZU
C-----------------------------------------------------------------------
C This routine serves as an interface between the driver and
C Subroutine DPREPI.  Tasks performed here are:
C  * call DPREPI,
C  * reset the required WM segment length LENWK,
C  * move YH back to its final location (following WM in RWORK),
C  * reset pointers for YH, SAVR, EWT, and ACOR, and
C  * move EWT to its new position if ISTATE = 0 or 1.
C IPFLAG is an output error indication flag.  IPFLAG = 0 if there was
C no trouble, and IPFLAG is the value of the DPREPI error flag IPPER
C if there was trouble in Subroutine DPREPI.
C-----------------------------------------------------------------------
      IPFLAG = 0
C Call DPREPI to do matrix preprocessing operations. -------------------
      CALL DPREPI (NEQ, Y, S, RWORK(LYH), RWORK(LSAVR), RWORK(LEWT),
     1   RWORK(LACOR), IA, JA, IC, JC, RWORK(LWM), RWORK(LWM), IPFLAG,
     2   RES, JAC, ADDA)
      LENWK = MAX(LREQ,LWMIN)
      IF (IPFLAG .LT. 0) RETURN
C If DPREPI was successful, move YH to end of required space for WM. ---
      LYHN = LWM + LENWK
      IF (LYHN .GT. LYH) RETURN
      LYHD = LYH - LYHN
      IF (LYHD .EQ. 0) GO TO 20
      IMAX = LYHN - 1 + LENYHM
      DO 10 I=LYHN,IMAX
 10     RWORK(I) = RWORK(I+LYHD)
      LYH = LYHN
C Reset pointers for SAVR, EWT, and ACOR. ------------------------------
 20   LSAVR = LYH + LENYH
      LEWTN = LSAVR + N
      LACOR = LEWTN + N
      IF (ISTATC .EQ. 3) GO TO 40
C If ISTATE = 1, move EWT (left) to its new position. ------------------
      IF (LEWTN .GT. LEWT) RETURN
      DO 30 I=1,N
 30     RWORK(I+LEWTN-1) = RWORK(I+LEWT-1)
 40   LEWT = LEWTN
      RETURN
C----------------------- End of Subroutine DIPREPI ---------------------
      END
*DECK DPREPI
      SUBROUTINE DPREPI (NEQ, Y, S, YH, SAVR, EWT, RTEM, IA, JA, IC, JC,
     1                   WK, IWK, IPPER, RES, JAC, ADDA)
      EXTERNAL RES, JAC, ADDA
      INTEGER NEQ, IA, JA, IC, JC, IWK, IPPER
      INTEGER ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, 
     1   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     2   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NRE, NJE, NQU
      INTEGER IPLOST, IESP, ISTATC, IYS, IBA, IBIAN, IBJAN, IBJGP,
     1   IPIAN, IPJAN, IPJGP, IPIGP, IPR, IPC, IPIC, IPISP, IPRSP, IPA,
     2   LENYH, LENYHM, LENWK, LREQ, LRAT, LREST, LWMIN, MOSS, MSBJ,
     3   NSLJ, NGP, NLU, NNZ, NSP, NZL, NZU
      INTEGER I, IBR, IER, IPIL, IPIU, IPTT1, IPTT2, J, K, KNEW, KAMAX,
     1   KAMIN, KCMAX, KCMIN, LDIF, LENIGP, LIWK, LJFO, MAXG, NP1, NZSUT
      DOUBLE PRECISION Y, S, YH, SAVR, EWT, RTEM, WK
      DOUBLE PRECISION CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
      DOUBLE PRECISION RLSS
      DOUBLE PRECISION ERWT, FAC, YJ
      DIMENSION NEQ(*), Y(*), S(*), YH(*), SAVR(*), EWT(*), RTEM(*),
     1   IA(*), JA(*), IC(*), JC(*), WK(*), IWK(*)
      COMMON /DLS001/ CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,
     1   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, 
     2   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     3   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NRE, NJE, NQU
      COMMON /DLSS01/ RLSS(6),
     1   IPLOST, IESP, ISTATC, IYS, IBA, IBIAN, IBJAN, IBJGP,
     2   IPIAN, IPJAN, IPJGP, IPIGP, IPR, IPC, IPIC, IPISP, IPRSP, IPA,
     3   LENYH, LENYHM, LENWK, LREQ, LRAT, LREST, LWMIN, MOSS, MSBJ,
     4   NSLJ, NGP, NLU, NNZ, NSP, NZL, NZU
C-----------------------------------------------------------------------
C This routine performs preprocessing related to the sparse linear
C systems that must be solved.
C The operations that are performed here are:
C  * compute sparseness structure of the iteration matrix
C      P = A - con*J  according to MOSS,
C  * compute grouping of column indices (MITER = 2),
C  * compute a new ordering of rows and columns of the matrix,
C  * reorder JA corresponding to the new ordering,
C  * perform a symbolic LU factorization of the matrix, and
C  * set pointers for segments of the IWK/WK array.
C In addition to variables described previously, DPREPI uses the
C following for communication:
C YH     = the history array.  Only the first column, containing the
C          current Y vector, is used.  Used only if MOSS .ne. 0.
C S      = array of length NEQ, identical to YDOTI in the driver, used
C          only if MOSS .ne. 0.
C SAVR   = a work array of length NEQ, used only if MOSS .ne. 0.
C EWT    = array of length NEQ containing (inverted) error weights.
C          Used only if MOSS = 2 or 4 or if ISTATE = MOSS = 1.
C RTEM   = a work array of length NEQ, identical to ACOR in the driver,
C          used only if MOSS = 2 or 4.
C WK     = a real work array of length LENWK, identical to WM in
C          the driver.
C IWK    = integer work array, assumed to occupy the same space as WK.
C LENWK  = the length of the work arrays WK and IWK.
C ISTATC = a copy of the driver input argument ISTATE (= 1 on the
C          first call, = 3 on a continuation call).
C IYS    = flag value from ODRV or CDRV.
C IPPER  = output error flag , with the following values and meanings:
C        =   0  no error.
C        =  -1  insufficient storage for internal structure pointers.
C        =  -2  insufficient storage for JGROUP.
C        =  -3  insufficient storage for ODRV.
C        =  -4  other error flag from ODRV (should never occur).
C        =  -5  insufficient storage for CDRV.
C        =  -6  other error flag from CDRV.
C        =  -7  if the RES routine returned error flag IRES = IER = 2.
C        =  -8  if the RES routine returned error flag IRES = IER = 3.
C-----------------------------------------------------------------------
      IBIAN = LRAT*2
      IPIAN = IBIAN + 1
      NP1 = N + 1
      IPJAN = IPIAN + NP1
      IBJAN = IPJAN - 1
      LIWK = LENWK*LRAT
      IF (MOSS .NE. 3 .AND. MOSS .NE. 4) LIWK = LIWK - N
      IF (IPJAN+N-1 .GT. LIWK) GO TO 310
      IF (MOSS .EQ. 0) GO TO 30
C
      IF (ISTATC .EQ. 3) GO TO 20
C ISTATE = 1 and MOSS .ne. 0.  Perturb Y for structure determination.
C Initialize S with random nonzero elements for structure determination.
      DO 10 I=1,N
        ERWT = 1.0D0/EWT(I)
        FAC = 1.0D0 + 1.0D0/(I + 1.0D0)
        Y(I) = Y(I) + FAC*SIGN(ERWT,Y(I))
        S(I) = 1.0D0 + FAC*ERWT
 10     CONTINUE
      GO TO (70,100,150,200), MOSS
C
 20   CONTINUE
C ISTATE = 3 and MOSS .ne. 0. Load Y from YH(*,1) and S from YH(*,2). --
      DO 25 I = 1,N
         Y(I) = YH(I)
 25      S(I) = YH(N+I)
      GO TO (70, 100, 150, 200),  MOSS
C
C MOSS = 0. Process user's IA,JA and IC,JC. ----------------------------
 30   KNEW = IPJAN
      KAMIN = IA(1)
      KCMIN = IC(1)
      IWK(IPIAN) = 1
      DO 60 J = 1,N
        DO 35 I = 1,N
 35       IWK(LIWK+I) = 0
        KAMAX = IA(J+1) - 1
        IF (KAMIN .GT. KAMAX) GO TO 45
        DO 40 K = KAMIN,KAMAX
          I = JA(K)
          IWK(LIWK+I) = 1
          IF (KNEW .GT. LIWK) GO TO 310
          IWK(KNEW) = I
          KNEW = KNEW + 1
 40       CONTINUE
 45     KAMIN = KAMAX + 1
        KCMAX = IC(J+1) - 1
        IF (KCMIN .GT. KCMAX) GO TO 55
        DO 50 K = KCMIN,KCMAX
          I = JC(K)
          IF (IWK(LIWK+I) .NE. 0) GO TO 50
          IF (KNEW .GT. LIWK) GO TO 310
          IWK(KNEW) = I
          KNEW = KNEW + 1
 50   CONTINUE
 55     IWK(IPIAN+J) = KNEW + 1 - IPJAN
        KCMIN = KCMAX + 1
 60     CONTINUE
      GO TO 240
C
C MOSS = 1. Compute structure from user-supplied Jacobian routine JAC. -
 70   CONTINUE
C A dummy call to RES allows user to create temporaries for use in JAC.
      IER = 1
      CALL RES (NEQ, TN, Y, S, SAVR, IER)
      IF (IER .GT. 1) GO TO 370
      DO 75 I = 1,N
        SAVR(I) = 0.
 75     WK(LIWK+I) = 0.
      K = IPJAN
      IWK(IPIAN) = 1
      DO 95 J = 1,N
        CALL ADDA (NEQ, TN, Y, J, IWK(IPIAN), IWK(IPJAN), WK(LIWK+1))
        CALL JAC (NEQ, TN, Y, S, J, IWK(IPIAN), IWK(IPJAN), SAVR)
        DO 90 I = 1,N
          LJFO = LIWK + I
          IF (WK(LJFO) .EQ. 0.) GO TO 80
          WK(LJFO) = 0.
          SAVR(I) = 0.
          GO TO 85
 80       IF (SAVR(I) .EQ. 0.) GO TO 90
          SAVR(I) = 0.
 85       IF (K .GT. LIWK) GO TO 310
          IWK(K) = I
          K = K+1
 90   CONTINUE
        IWK(IPIAN+J) = K + 1 - IPJAN
 95   CONTINUE
      GO TO 240
C
C MOSS = 2. Compute structure from results of N + 1 calls to RES. ------
 100  DO 105 I = 1,N
 105    WK(LIWK+I) = 0.
      K = IPJAN
      IWK(IPIAN) = 1
      IER = -1
      IF (MITER .EQ. 1) IER = 1
      CALL RES (NEQ, TN, Y, S, SAVR, IER)
      IF (IER .GT. 1) GO TO 370
      DO 130 J = 1,N
        CALL ADDA (NEQ, TN, Y, J, IWK(IPIAN), IWK(IPJAN), WK(LIWK+1))
        YJ = Y(J)
        ERWT = 1.0D0/EWT(J)
        Y(J) = YJ + SIGN(ERWT,YJ)
        CALL RES (NEQ, TN, Y, S, RTEM, IER)
        IF (IER .GT. 1) RETURN
        Y(J) = YJ
        DO 120 I = 1,N
          LJFO = LIWK + I
          IF (WK(LJFO) .EQ. 0.) GO TO 110
          WK(LJFO) = 0.
          GO TO 115
 110      IF (RTEM(I) .EQ. SAVR(I)) GO TO 120
 115      IF (K .GT. LIWK) GO TO 310
          IWK(K) = I
          K = K + 1
 120  CONTINUE
      IWK(IPIAN+J) = K + 1 - IPJAN
 130  CONTINUE
      GO TO 240
C
C MOSS = 3. Compute structure from the user's IA/JA and JAC routine. ---
 150  CONTINUE
C A dummy call to RES allows user to create temporaries for use in JAC.
      IER = 1
      CALL RES (NEQ, TN, Y, S, SAVR, IER)
      IF (IER .GT. 1) GO TO 370
      DO 155 I = 1,N
 155    SAVR(I) = 0.
      KNEW = IPJAN
      KAMIN = IA(1)
      IWK(IPIAN) = 1
      DO 190 J = 1,N
        CALL JAC (NEQ, TN, Y, S, J, IWK(IPIAN), IWK(IPJAN), SAVR)
        KAMAX = IA(J+1) - 1
        IF (KAMIN .GT. KAMAX) GO TO 170
        DO 160 K = KAMIN,KAMAX
          I = JA(K)
          SAVR(I) = 0.
          IF (KNEW .GT. LIWK) GO TO 310
          IWK(KNEW) = I
          KNEW = KNEW + 1
 160  CONTINUE
 170  KAMIN = KAMAX + 1
      DO 180 I = 1,N
        IF (SAVR(I) .EQ. 0.) GO TO 180
        SAVR(I) = 0.
        IF (KNEW .GT. LIWK) GO TO 310
        IWK(KNEW) = I
        KNEW = KNEW + 1
 180  CONTINUE
      IWK(IPIAN+J) = KNEW + 1 - IPJAN
 190  CONTINUE
      GO TO 240
C
C MOSS = 4. Compute structure from user's IA/JA and N + 1 RES calls. ---
 200  KNEW = IPJAN
      KAMIN = IA(1)
      IWK(IPIAN) = 1
      IER = -1
      IF (MITER .EQ. 1) IER = 1
      CALL RES (NEQ, TN, Y, S, SAVR, IER)
      IF (IER .GT. 1) GO TO 370
      DO 235 J = 1,N
        YJ = Y(J)
        ERWT = 1.0D0/EWT(J)
        Y(J) = YJ + SIGN(ERWT,YJ)
        CALL RES (NEQ, TN, Y, S, RTEM, IER)
        IF (IER .GT. 1) RETURN
        Y(J) = YJ
        KAMAX = IA(J+1) - 1
        IF (KAMIN .GT. KAMAX) GO TO 225
        DO 220 K = KAMIN,KAMAX
          I = JA(K)
          RTEM(I) = SAVR(I)
          IF (KNEW .GT. LIWK) GO TO 310
          IWK(KNEW) = I
          KNEW = KNEW + 1
 220  CONTINUE
 225  KAMIN = KAMAX + 1
        DO 230 I = 1,N
          IF (RTEM(I) .EQ. SAVR(I)) GO TO 230
          IF (KNEW .GT. LIWK) GO TO 310
          IWK(KNEW) = I
          KNEW = KNEW + 1
 230  CONTINUE
        IWK(IPIAN+J) = KNEW + 1 - IPJAN
 235  CONTINUE
C
 240  CONTINUE
      IF (MOSS .EQ. 0 .OR. ISTATC .EQ. 3) GO TO 250
C If ISTATE = 0 or 1 and MOSS .ne. 0, restore Y from YH. ---------------
      DO 245 I = 1,N
 245    Y(I) = YH(I)
 250  NNZ = IWK(IPIAN+N) - 1
       IPPER = 0
      NGP = 0
      LENIGP = 0
      IPIGP = IPJAN + NNZ
      IF (MITER .NE. 2) GO TO 260
C
C Compute grouping of column indices (MITER = 2). ----------------------
C
      MAXG = NP1
      IPJGP = IPJAN + NNZ
      IBJGP = IPJGP - 1
      IPIGP = IPJGP + N
      IPTT1 = IPIGP + NP1
      IPTT2 = IPTT1 + N
      LREQ = IPTT2 + N - 1
      IF (LREQ .GT. LIWK) GO TO 320
      CALL JGROUP (N, IWK(IPIAN), IWK(IPJAN), MAXG, NGP, IWK(IPIGP),
     1   IWK(IPJGP), IWK(IPTT1), IWK(IPTT2), IER)
      IF (IER .NE. 0) GO TO 320
      LENIGP = NGP + 1
C
C Compute new ordering of rows/columns of Jacobian. --------------------
 260  IPR = IPIGP + LENIGP
      IPC = IPR
      IPIC = IPC + N
      IPISP = IPIC + N
      IPRSP = (IPISP-2)/LRAT + 2
      IESP = LENWK + 1 - IPRSP
      IF (IESP .LT. 0) GO TO 330
      IBR = IPR - 1
      DO 270 I = 1,N
 270    IWK(IBR+I) = I
      NSP = LIWK + 1 - IPISP
      CALL ODRV(N, IWK(IPIAN), IWK(IPJAN), WK, IWK(IPR), IWK(IPIC), NSP,
     1   IWK(IPISP), 1, IYS)
      IF (IYS .EQ. 11*N+1) GO TO 340
      IF (IYS .NE. 0) GO TO 330
C
C Reorder JAN and do symbolic LU factorization of matrix. --------------
      IPA = LENWK + 1 - NNZ
      NSP = IPA - IPRSP
      LREQ = MAX(12*N/LRAT, 6*N/LRAT+2*N+NNZ) + 3
      LREQ = LREQ + IPRSP - 1 + NNZ
      IF (LREQ .GT. LENWK) GO TO 350
      IBA = IPA - 1
      DO 280 I = 1,NNZ
 280    WK(IBA+I) = 0.0D0
      IPISP = LRAT*(IPRSP - 1) + 1
      CALL CDRV(N,IWK(IPR),IWK(IPC),IWK(IPIC),IWK(IPIAN),IWK(IPJAN),
     1   WK(IPA),WK(IPA),WK(IPA),NSP,IWK(IPISP),WK(IPRSP),IESP,5,IYS)
      LREQ = LENWK - IESP
      IF (IYS .EQ. 10*N+1) GO TO 350
      IF (IYS .NE. 0) GO TO 360
      IPIL = IPISP
      IPIU = IPIL + 2*N + 1
      NZU = IWK(IPIL+N) - IWK(IPIL)
      NZL = IWK(IPIU+N) - IWK(IPIU)
      IF (LRAT .GT. 1) GO TO 290
      CALL ADJLR (N, IWK(IPISP), LDIF)
      LREQ = LREQ + LDIF
 290  CONTINUE
      IF (LRAT .EQ. 2 .AND. NNZ .EQ. N) LREQ = LREQ + 1
      NSP = NSP + LREQ - LENWK
      IPA = LREQ + 1 - NNZ
      IBA = IPA - 1
      IPPER = 0
      RETURN
C
 310  IPPER = -1
      LREQ = 2 + (2*N + 1)/LRAT
      LREQ = MAX(LENWK+1,LREQ)
      RETURN
C
 320  IPPER = -2
      LREQ = (LREQ - 1)/LRAT + 1
      RETURN
C
 330  IPPER = -3
      CALL CNTNZU (N, IWK(IPIAN), IWK(IPJAN), NZSUT)
      LREQ = LENWK - IESP + (3*N + 4*NZSUT - 1)/LRAT + 1
      RETURN
C
 340  IPPER = -4
      RETURN
C
 350  IPPER =  -5
      RETURN
C
 360  IPPER = -6
      LREQ = LENWK
      RETURN
C
 370  IPPER = -IER - 5
      LREQ = 2 + (2*N + 1)/LRAT
      RETURN
C----------------------- End of Subroutine DPREPI ----------------------
      END
*DECK DAINVGS
      SUBROUTINE DAINVGS (NEQ, T, Y, WK, IWK, TEM, YDOT, IER, RES, ADDA)
      EXTERNAL RES, ADDA
      INTEGER NEQ, IWK, IER
      INTEGER IPLOST, IESP, ISTATC, IYS, IBA, IBIAN, IBJAN, IBJGP,
     1   IPIAN, IPJAN, IPJGP, IPIGP, IPR, IPC, IPIC, IPISP, IPRSP, IPA,
     2   LENYH, LENYHM, LENWK, LREQ, LRAT, LREST, LWMIN, MOSS, MSBJ,
     3   NSLJ, NGP, NLU, NNZ, NSP, NZL, NZU
      INTEGER I, IMUL, J, K, KMIN, KMAX
      DOUBLE PRECISION T, Y, WK, TEM, YDOT
      DOUBLE PRECISION RLSS
      DIMENSION Y(*), WK(*), IWK(*), TEM(*), YDOT(*)
      COMMON /DLSS01/ RLSS(6),
     1   IPLOST, IESP, ISTATC, IYS, IBA, IBIAN, IBJAN, IBJGP,
     2   IPIAN, IPJAN, IPJGP, IPIGP, IPR, IPC, IPIC, IPISP, IPRSP, IPA,
     3   LENYH, LENYHM, LENWK, LREQ, LRAT, LREST, LWMIN, MOSS, MSBJ,
     4   NSLJ, NGP, NLU, NNZ, NSP, NZL, NZU
C-----------------------------------------------------------------------
C This subroutine computes the initial value of the vector YDOT
C satisfying
C     A * YDOT = g(t,y)
C when A is nonsingular.  It is called by DLSODIS for initialization
C only, when ISTATE = 0.  The matrix A is subjected to LU
C decomposition in CDRV.  Then the system A*YDOT = g(t,y) is solved
C in CDRV.
C In addition to variables described previously, communication
C with DAINVGS uses the following:
C Y     = array of initial values.
C WK    = real work space for matrices.  On output it contains A and
C         its LU decomposition.  The LU decomposition is not entirely
C         sparse unless the structure of the matrix A is identical to
C         the structure of the Jacobian matrix dr/dy.
C         Storage of matrix elements starts at WK(3).
C         WK(1) = SQRT(UROUND), not used here.
C IWK   = integer work space for matrix-related data, assumed to
C         be equivalenced to WK.  In addition, WK(IPRSP) and WK(IPISP)
C         are assumed to have identical locations.
C TEM   = vector of work space of length N (ACOR in DSTODI).
C YDOT  = output vector containing the initial dy/dt. YDOT(i) contains
C         dy(i)/dt when the matrix A is non-singular.
C IER   = output error flag with the following values and meanings:
C       = 0  if DAINVGS was successful.
C       = 1  if the A-matrix was found to be singular.
C       = 2  if RES returned an error flag IRES = IER = 2.
C       = 3  if RES returned an error flag IRES = IER = 3.
C       = 4  if insufficient storage for CDRV (should not occur here).
C       = 5  if other error found in CDRV (should not occur here).
C-----------------------------------------------------------------------
C
      DO 10 I = 1,NNZ
 10     WK(IBA+I) = 0.0D0
C
      IER = 1
      CALL RES (NEQ, T, Y, WK(IPA), YDOT, IER)
      IF (IER .GT. 1) RETURN
C
      KMIN = IWK(IPIAN)
      DO 30 J = 1,NEQ
        KMAX = IWK(IPIAN+J) - 1
        DO 15 K = KMIN,KMAX
          I = IWK(IBJAN+K)
 15       TEM(I) = 0.0D0
        CALL ADDA (NEQ, T, Y, J, IWK(IPIAN), IWK(IPJAN), TEM)
        DO 20 K = KMIN,KMAX
          I = IWK(IBJAN+K)
 20       WK(IBA+K) = TEM(I)
        KMIN = KMAX + 1
 30   CONTINUE
      NLU = NLU + 1
      IER = 0
      DO 40 I = 1,NEQ
 40     TEM(I) = 0.0D0
C
C Numerical factorization of matrix A. ---------------------------------
      CALL CDRV (NEQ,IWK(IPR),IWK(IPC),IWK(IPIC),IWK(IPIAN),IWK(IPJAN),
     1  WK(IPA),TEM,TEM,NSP,IWK(IPISP),WK(IPRSP),IESP,2,IYS)
      IF (IYS .EQ. 0) GO TO 50
      IMUL = (IYS - 1)/NEQ
      IER = 5
      IF (IMUL .EQ. 8) IER = 1
      IF (IMUL .EQ. 10) IER = 4
      RETURN
C
C Solution of the linear system. ---------------------------------------
 50   CALL CDRV (NEQ,IWK(IPR),IWK(IPC),IWK(IPIC),IWK(IPIAN),IWK(IPJAN),
     1  WK(IPA),YDOT,YDOT,NSP,IWK(IPISP),WK(IPRSP),IESP,4,IYS)
      IF (IYS .NE. 0) IER = 5
      RETURN
C----------------------- End of Subroutine DAINVGS ---------------------
      END
*DECK DPRJIS
      SUBROUTINE DPRJIS (NEQ, Y, YH, NYH, EWT, RTEM, SAVR, S, WK, IWK,
     1   RES, JAC, ADDA)
      EXTERNAL RES, JAC, ADDA
      INTEGER NEQ, NYH, IWK
      INTEGER ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, 
     1   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     2   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NRE, NJE, NQU
      INTEGER IPLOST, IESP, ISTATC, IYS, IBA, IBIAN, IBJAN, IBJGP,
     1   IPIAN, IPJAN, IPJGP, IPIGP, IPR, IPC, IPIC, IPISP, IPRSP, IPA,
     2   LENYH, LENYHM, LENWK, LREQ, LRAT, LREST, LWMIN, MOSS, MSBJ,
     3   NSLJ, NGP, NLU, NNZ, NSP, NZL, NZU
      INTEGER I, IMUL, IRES, J, JJ, JMAX, JMIN, K, KMAX, KMIN, NG
      DOUBLE PRECISION Y, YH, EWT, RTEM, SAVR, S, WK
      DOUBLE PRECISION CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
      DOUBLE PRECISION RLSS
      DOUBLE PRECISION CON, FAC, HL0, R, SRUR
      DIMENSION NEQ(*), Y(*), YH(NYH,*), EWT(*), RTEM(*),
     1   S(*), SAVR(*), WK(*), IWK(*)
      COMMON /DLS001/ CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,
     1   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, 
     2   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     3   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NRE, NJE, NQU
      COMMON /DLSS01/ RLSS(6),
     1   IPLOST, IESP, ISTATC, IYS, IBA, IBIAN, IBJAN, IBJGP,
     2   IPIAN, IPJAN, IPJGP, IPIGP, IPR, IPC, IPIC, IPISP, IPRSP, IPA,
     3   LENYH, LENYHM, LENWK, LREQ, LRAT, LREST, LWMIN, MOSS, MSBJ,
     4   NSLJ, NGP, NLU, NNZ, NSP, NZL, NZU
C-----------------------------------------------------------------------
C DPRJIS is called to compute and process the matrix
C P = A - H*EL(1)*J, where J is an approximation to the Jacobian dr/dy,
C where r = g(t,y) - A(t,y)*s.  J is computed by columns, either by
C the user-supplied routine JAC if MITER = 1, or by finite differencing
C if MITER = 2.  J is stored in WK, rescaled, and ADDA is called to
C generate P.  The matrix P is subjected to LU decomposition in CDRV.
C P and its LU decomposition are stored separately in WK.
C
C In addition to variables described previously, communication
C with DPRJIS uses the following:
C Y     = array containing predicted values on entry.
C RTEM  = work array of length N (ACOR in DSTODI).
C SAVR  = array containing r evaluated at predicted y. On output it
C         contains the residual evaluated at current values of t and y.
C S     = array containing predicted values of dy/dt (SAVF in DSTODI).
C WK    = real work space for matrices.  On output it contains P and
C         its sparse LU decomposition.  Storage of matrix elements
C         starts at WK(3).
C         WK also contains the following matrix-related data.
C         WK(1) = SQRT(UROUND), used in numerical Jacobian increments.
C IWK   = integer work space for matrix-related data, assumed to be
C         equivalenced to WK.  In addition,  WK(IPRSP) and IWK(IPISP)
C         are assumed to have identical locations.
C EL0   = EL(1) (input).
C IERPJ = output error flag (in COMMON).
C         =  0 if no error.
C         =  1 if zero pivot found in CDRV.
C         = IRES (= 2 or 3) if RES returned IRES = 2 or 3.
C         = -1 if insufficient storage for CDRV (should not occur).
C         = -2 if other error found in CDRV (should not occur here).
C JCUR  = output flag = 1 to indicate that the Jacobian matrix
C         (or approximation) is now current.
C This routine also uses other variables in Common.
C-----------------------------------------------------------------------
      HL0 = H*EL0
      CON = -HL0
      JCUR = 1
      NJE = NJE + 1
      GO TO (100, 200), MITER
C
C If MITER = 1, call RES, then call JAC and ADDA for each column. ------
 100  IRES = 1
      CALL RES (NEQ, TN, Y, S, SAVR, IRES)
      NRE = NRE + 1
      IF (IRES .GT. 1) GO TO 600
      KMIN = IWK(IPIAN)
      DO 130 J = 1,N
        KMAX = IWK(IPIAN+J)-1
        DO 110 I = 1,N
 110      RTEM(I) = 0.0D0
        CALL JAC (NEQ, TN, Y, S, J, IWK(IPIAN), IWK(IPJAN), RTEM)
        DO 120 I = 1,N
 120      RTEM(I) = RTEM(I)*CON
        CALL ADDA (NEQ, TN, Y, J, IWK(IPIAN), IWK(IPJAN), RTEM)
        DO 125 K = KMIN,KMAX
          I = IWK(IBJAN+K)
          WK(IBA+K) = RTEM(I)
 125      CONTINUE
        KMIN = KMAX + 1
 130    CONTINUE
      GO TO 290
C
C If MITER = 2, make NGP + 1 calls to RES to approximate J and P. ------
 200  CONTINUE
      IRES = -1
      CALL RES (NEQ, TN, Y, S, SAVR, IRES)
      NRE = NRE + 1
      IF (IRES .GT. 1) GO TO 600
      SRUR = WK(1)
      JMIN = IWK(IPIGP)
      DO 240 NG = 1,NGP
        JMAX = IWK(IPIGP+NG) - 1
        DO 210 J = JMIN,JMAX
          JJ = IWK(IBJGP+J)
          R = MAX(SRUR*ABS(Y(JJ)),0.01D0/EWT(JJ))
 210      Y(JJ) = Y(JJ) + R
        CALL RES (NEQ,TN,Y,S,RTEM,IRES)
        NRE = NRE + 1
        IF (IRES .GT. 1) GO TO 600
        DO 230 J = JMIN,JMAX
          JJ = IWK(IBJGP+J)
          Y(JJ) = YH(JJ,1)
          R = MAX(SRUR*ABS(Y(JJ)),0.01D0/EWT(JJ))
          FAC = -HL0/R
          KMIN = IWK(IBIAN+JJ)
          KMAX = IWK(IBIAN+JJ+1) - 1
          DO 220 K = KMIN,KMAX
            I = IWK(IBJAN+K)
            RTEM(I) = (RTEM(I) - SAVR(I))*FAC
 220        CONTINUE
        CALL ADDA (NEQ, TN, Y, JJ, IWK(IPIAN), IWK(IPJAN), RTEM)
        DO 225 K = KMIN,KMAX
          I = IWK(IBJAN+K)
          WK(IBA+K) = RTEM(I)
 225      CONTINUE
 230      CONTINUE
        JMIN = JMAX + 1
 240    CONTINUE
      IRES = 1
      CALL RES (NEQ, TN, Y, S, SAVR, IRES)
      NRE = NRE + 1
      IF (IRES .GT. 1) GO TO 600
C
C Do numerical factorization of P matrix. ------------------------------
 290  NLU = NLU + 1
      IERPJ = 0
      DO 295 I = 1,N
 295    RTEM(I) = 0.0D0
      CALL CDRV (N,IWK(IPR),IWK(IPC),IWK(IPIC),IWK(IPIAN),IWK(IPJAN),
     1  WK(IPA),RTEM,RTEM,NSP,IWK(IPISP),WK(IPRSP),IESP,2,IYS)
      IF (IYS .EQ. 0) RETURN
      IMUL = (IYS - 1)/N
      IERPJ = -2
      IF (IMUL .EQ. 8) IERPJ = 1
      IF (IMUL .EQ. 10) IERPJ = -1
      RETURN
C Error return for IRES = 2 or IRES = 3 return from RES. ---------------
 600  IERPJ = IRES
      RETURN
C----------------------- End of Subroutine DPRJIS ----------------------
      END

*DECK DGEFA
      SUBROUTINE DGEFA (A, LDA, N, IPVT, INFO)
C***BEGIN PROLOGUE  DGEFA
C***PURPOSE  Factor a matrix using Gaussian elimination.
C***CATEGORY  D2A1
C***TYPE      DOUBLE PRECISION (SGEFA-S, DGEFA-D, CGEFA-C)
C***KEYWORDS  GENERAL MATRIX, LINEAR ALGEBRA, LINPACK,
C             MATRIX FACTORIZATION
C***AUTHOR  Moler, C. B., (U. of New Mexico)
C***DESCRIPTION
C
C     DGEFA factors a double precision matrix by Gaussian elimination.
C
C     DGEFA is usually called by DGECO, but it can be called
C     directly with a saving in time if  RCOND  is not needed.
C     (Time for DGECO) = (1 + 9/N)*(Time for DGEFA) .
C
C     On Entry
C
C        A       DOUBLE PRECISION(LDA, N)
C                the matrix to be factored.
C
C        LDA     INTEGER
C                the leading dimension of the array  A .
C
C        N       INTEGER
C                the order of the matrix  A .
C
C     On Return
C
C        A       an upper triangular matrix and the multipliers
C                which were used to obtain it.
C                The factorization can be written  A = L*U  where
C                L  is a product of permutation and unit lower
C                triangular matrices and  U  is upper triangular.
C
C        IPVT    INTEGER(N)
C                an integer vector of pivot indices.
C
C        INFO    INTEGER
C                = 0  normal value.
C                = K  if  U(K,K) .EQ. 0.0 .  This is not an error
C                     condition for this subroutine, but it does
C                     indicate that DGESL or DGEDI will divide by zero
C                     if called.  Use  RCOND  in DGECO for a reliable
C                     indication of singularity.
C
C***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
C                 Stewart, LINPACK Users' Guide, SIAM, 1979.
C***ROUTINES CALLED  DAXPY, DSCAL, IDAMAX
C***REVISION HISTORY  (YYMMDD)
C   780814  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DGEFA
      INTEGER LDA,N,IPVT(*),INFO
      DOUBLE PRECISION A(LDA,*)
C
      DOUBLE PRECISION T
      INTEGER IDAMAX,J,K,KP1,L,NM1
C
C     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
C
C***FIRST EXECUTABLE STATEMENT  DGEFA
      INFO = 0
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 70
      DO 60 K = 1, NM1
         KP1 = K + 1
C
C        FIND L = PIVOT INDEX
C
         L = IDAMAX(N-K+1,A(K,K),1) + K - 1
         IPVT(K) = L
C
C        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
C
         IF (A(L,K) .EQ. 0.0D0) GO TO 40
C
C           INTERCHANGE IF NECESSARY
C
            IF (L .EQ. K) GO TO 10
               T = A(L,K)
               A(L,K) = A(K,K)
               A(K,K) = T
   10       CONTINUE
C
C           COMPUTE MULTIPLIERS
C
            T = -1.0D0/A(K,K)
            CALL DSCAL(N-K,T,A(K+1,K),1)
C
C           ROW ELIMINATION WITH COLUMN INDEXING
C
            DO 30 J = KP1, N
               T = A(L,J)
               IF (L .EQ. K) GO TO 20
                  A(L,J) = A(K,J)
                  A(K,J) = T
   20          CONTINUE
               CALL DAXPY(N-K,T,A(K+1,K),1,A(K+1,J),1)
   30       CONTINUE
         GO TO 50
   40    CONTINUE
            INFO = K
   50    CONTINUE
   60 CONTINUE
   70 CONTINUE
      IPVT(N) = N
      IF (A(N,N) .EQ. 0.0D0) INFO = N
      RETURN
      END
*DECK DGESL
      SUBROUTINE DGESL (A, LDA, N, IPVT, B, JOB)
C***BEGIN PROLOGUE  DGESL
C***PURPOSE  Solve the real system A*X=B or TRANS(A)*X=B using the
C            factors computed by DGECO or DGEFA.
C***CATEGORY  D2A1
C***TYPE      DOUBLE PRECISION (SGESL-S, DGESL-D, CGESL-C)
C***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX, SOLVE
C***AUTHOR  Moler, C. B., (U. of New Mexico)
C***DESCRIPTION
C
C     DGESL solves the double precision system
C     A * X = B  or  TRANS(A) * X = B
C     using the factors computed by DGECO or DGEFA.
C
C     On Entry
C
C        A       DOUBLE PRECISION(LDA, N)
C                the output from DGECO or DGEFA.
C
C        LDA     INTEGER
C                the leading dimension of the array  A .
C
C        N       INTEGER
C                the order of the matrix  A .
C
C        IPVT    INTEGER(N)
C                the pivot vector from DGECO or DGEFA.
C
C        B       DOUBLE PRECISION(N)
C                the right hand side vector.
C
C        JOB     INTEGER
C                = 0         to solve  A*X = B ,
C                = nonzero   to solve  TRANS(A)*X = B  where
C                            TRANS(A)  is the transpose.
C
C     On Return
C
C        B       the solution vector  X .
C
C     Error Condition
C
C        A division by zero will occur if the input factor contains a
C        zero on the diagonal.  Technically this indicates singularity
C        but it is often caused by improper arguments or improper
C        setting of LDA .  It will not occur if the subroutines are
C        called correctly and if DGECO has set RCOND .GT. 0.0
C        or DGEFA has set INFO .EQ. 0 .
C
C     To compute  INVERSE(A) * C  where  C  is a matrix
C     with  P  columns
C           CALL DGECO(A,LDA,N,IPVT,RCOND,Z)
C           IF (RCOND is too small) GO TO ...
C           DO 10 J = 1, P
C              CALL DGESL(A,LDA,N,IPVT,C(1,J),0)
C        10 CONTINUE
C
C***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
C                 Stewart, LINPACK Users' Guide, SIAM, 1979.
C***ROUTINES CALLED  DAXPY, DDOT
C***REVISION HISTORY  (YYMMDD)
C   780814  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DGESL
      INTEGER LDA,N,IPVT(*),JOB
      DOUBLE PRECISION A(LDA,*),B(*)
C
      DOUBLE PRECISION DDOT,T
      INTEGER K,KB,L,NM1
C***FIRST EXECUTABLE STATEMENT  DGESL
      NM1 = N - 1
      IF (JOB .NE. 0) GO TO 50
C
C        JOB = 0 , SOLVE  A * X = B
C        FIRST SOLVE  L*Y = B
C
         IF (NM1 .LT. 1) GO TO 30
         DO 20 K = 1, NM1
            L = IPVT(K)
            T = B(L)
            IF (L .EQ. K) GO TO 10
               B(L) = B(K)
               B(K) = T
   10       CONTINUE
            CALL DAXPY(N-K,T,A(K+1,K),1,B(K+1),1)
   20    CONTINUE
   30    CONTINUE
C
C        NOW SOLVE  U*X = Y
C
         DO 40 KB = 1, N
            K = N + 1 - KB
            B(K) = B(K)/A(K,K)
            T = -B(K)
            CALL DAXPY(K-1,T,A(1,K),1,B(1),1)
   40    CONTINUE
      GO TO 100
   50 CONTINUE
C
C        JOB = NONZERO, SOLVE  TRANS(A) * X = B
C        FIRST SOLVE  TRANS(U)*Y = B
C
         DO 60 K = 1, N
            T = DDOT(K-1,A(1,K),1,B(1),1)
            B(K) = (B(K) - T)/A(K,K)
   60    CONTINUE
C
C        NOW SOLVE TRANS(L)*X = Y
C
         IF (NM1 .LT. 1) GO TO 90
         DO 80 KB = 1, NM1
            K = N - KB
            B(K) = B(K) + DDOT(N-K,A(K+1,K),1,B(K+1),1)
            L = IPVT(K)
            IF (L .EQ. K) GO TO 70
               T = B(L)
               B(L) = B(K)
               B(K) = T
   70       CONTINUE
   80    CONTINUE
   90    CONTINUE
  100 CONTINUE
      RETURN
      END
*DECK DGBFA
      SUBROUTINE DGBFA (ABD, LDA, N, ML, MU, IPVT, INFO)
C***BEGIN PROLOGUE  DGBFA
C***PURPOSE  Factor a band matrix using Gaussian elimination.
C***CATEGORY  D2A2
C***TYPE      DOUBLE PRECISION (SGBFA-S, DGBFA-D, CGBFA-C)
C***KEYWORDS  BANDED, LINEAR ALGEBRA, LINPACK, MATRIX FACTORIZATION
C***AUTHOR  Moler, C. B., (U. of New Mexico)
C***DESCRIPTION
C
C     DGBFA factors a double precision band matrix by elimination.
C
C     DGBFA is usually called by DGBCO, but it can be called
C     directly with a saving in time if  RCOND  is not needed.
C
C     On Entry
C
C        ABD     DOUBLE PRECISION(LDA, N)
C                contains the matrix in band storage.  The columns
C                of the matrix are stored in the columns of  ABD  and
C                the diagonals of the matrix are stored in rows
C                ML+1 through 2*ML+MU+1 of  ABD .
C                See the comments below for details.
C
C        LDA     INTEGER
C                the leading dimension of the array  ABD .
C                LDA must be .GE. 2*ML + MU + 1 .
C
C        N       INTEGER
C                the order of the original matrix.
C
C        ML      INTEGER
C                number of diagonals below the main diagonal.
C                0 .LE. ML .LT.  N .
C
C        MU      INTEGER
C                number of diagonals above the main diagonal.
C                0 .LE. MU .LT.  N .
C                More efficient if  ML .LE. MU .
C     On Return
C
C        ABD     an upper triangular matrix in band storage and
C                the multipliers which were used to obtain it.
C                The factorization can be written  A = L*U  where
C                L  is a product of permutation and unit lower
C                triangular matrices and  U  is upper triangular.
C
C        IPVT    INTEGER(N)
C                an integer vector of pivot indices.
C
C        INFO    INTEGER
C                = 0  normal value.
C                = K  if  U(K,K) .EQ. 0.0 .  This is not an error
C                     condition for this subroutine, but it does
C                     indicate that DGBSL will divide by zero if
C                     called.  Use  RCOND  in DGBCO for a reliable
C                     indication of singularity.
C
C     Band Storage
C
C           If  A  is a band matrix, the following program segment
C           will set up the input.
C
C                   ML = (band width below the diagonal)
C                   MU = (band width above the diagonal)
C                   M = ML + MU + 1
C                   DO 20 J = 1, N
C                      I1 = MAX(1, J-MU)
C                      I2 = MIN(N, J+ML)
C                      DO 10 I = I1, I2
C                         K = I - J + M
C                         ABD(K,J) = A(I,J)
C                10    CONTINUE
C                20 CONTINUE
C
C           This uses rows  ML+1  through  2*ML+MU+1  of  ABD .
C           In addition, the first  ML  rows in  ABD  are used for
C           elements generated during the triangularization.
C           The total number of rows needed in  ABD  is  2*ML+MU+1 .
C           The  ML+MU by ML+MU  upper left triangle and the
C           ML by ML  lower right triangle are not referenced.
C
C***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
C                 Stewart, LINPACK Users' Guide, SIAM, 1979.
C***ROUTINES CALLED  DAXPY, DSCAL, IDAMAX
C***REVISION HISTORY  (YYMMDD)
C   780814  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DGBFA
      INTEGER LDA,N,ML,MU,IPVT(*),INFO
      DOUBLE PRECISION ABD(LDA,*)
C
      DOUBLE PRECISION T
      INTEGER I,IDAMAX,I0,J,JU,JZ,J0,J1,K,KP1,L,LM,M,MM,NM1
C
C***FIRST EXECUTABLE STATEMENT  DGBFA
      M = ML + MU + 1
      INFO = 0
C
C     ZERO INITIAL FILL-IN COLUMNS
C
      J0 = MU + 2
      J1 = MIN(N,M) - 1
      IF (J1 .LT. J0) GO TO 30
      DO 20 JZ = J0, J1
         I0 = M + 1 - JZ
         DO 10 I = I0, ML
            ABD(I,JZ) = 0.0D0
   10    CONTINUE
   20 CONTINUE
   30 CONTINUE
      JZ = J1
      JU = 0
C
C     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
C
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 130
      DO 120 K = 1, NM1
         KP1 = K + 1
C
C        ZERO NEXT FILL-IN COLUMN
C
         JZ = JZ + 1
         IF (JZ .GT. N) GO TO 50
         IF (ML .LT. 1) GO TO 50
            DO 40 I = 1, ML
               ABD(I,JZ) = 0.0D0
   40       CONTINUE
   50    CONTINUE
C
C        FIND L = PIVOT INDEX
C
         LM = MIN(ML,N-K)
         L = IDAMAX(LM+1,ABD(M,K),1) + M - 1
         IPVT(K) = L + K - M
C
C        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
C
         IF (ABD(L,K) .EQ. 0.0D0) GO TO 100
C
C           INTERCHANGE IF NECESSARY
C
            IF (L .EQ. M) GO TO 60
               T = ABD(L,K)
               ABD(L,K) = ABD(M,K)
               ABD(M,K) = T
   60       CONTINUE
C
C           COMPUTE MULTIPLIERS
C
            T = -1.0D0/ABD(M,K)
            CALL DSCAL(LM,T,ABD(M+1,K),1)
C
C           ROW ELIMINATION WITH COLUMN INDEXING
C
            JU = MIN(MAX(JU,MU+IPVT(K)),N)
            MM = M
            IF (JU .LT. KP1) GO TO 90
            DO 80 J = KP1, JU
               L = L - 1
               MM = MM - 1
               T = ABD(L,J)
               IF (L .EQ. MM) GO TO 70
                  ABD(L,J) = ABD(MM,J)
                  ABD(MM,J) = T
   70          CONTINUE
               CALL DAXPY(LM,T,ABD(M+1,K),1,ABD(MM+1,J),1)
   80       CONTINUE
   90       CONTINUE
         GO TO 110
  100    CONTINUE
            INFO = K
  110    CONTINUE
  120 CONTINUE
  130 CONTINUE
      IPVT(N) = N
      IF (ABD(M,N) .EQ. 0.0D0) INFO = N
      RETURN
      END
*DECK DGBSL
      SUBROUTINE DGBSL (ABD, LDA, N, ML, MU, IPVT, B, JOB)
C***BEGIN PROLOGUE  DGBSL
C***PURPOSE  Solve the real band system A*X=B or TRANS(A)*X=B using
C            the factors computed by DGBCO or DGBFA.
C***CATEGORY  D2A2
C***TYPE      DOUBLE PRECISION (SGBSL-S, DGBSL-D, CGBSL-C)
C***KEYWORDS  BANDED, LINEAR ALGEBRA, LINPACK, MATRIX, SOLVE
C***AUTHOR  Moler, C. B., (U. of New Mexico)
C***DESCRIPTION
C
C     DGBSL solves the double precision band system
C     A * X = B  or  TRANS(A) * X = B
C     using the factors computed by DGBCO or DGBFA.
C
C     On Entry
C
C        ABD     DOUBLE PRECISION(LDA, N)
C                the output from DGBCO or DGBFA.
C
C        LDA     INTEGER
C                the leading dimension of the array  ABD .
C
C        N       INTEGER
C                the order of the original matrix.
C
C        ML      INTEGER
C                number of diagonals below the main diagonal.
C
C        MU      INTEGER
C                number of diagonals above the main diagonal.
C
C        IPVT    INTEGER(N)
C                the pivot vector from DGBCO or DGBFA.
C
C        B       DOUBLE PRECISION(N)
C                the right hand side vector.
C
C        JOB     INTEGER
C                = 0         to solve  A*X = B ,
C                = nonzero   to solve  TRANS(A)*X = B , where
C                            TRANS(A)  is the transpose.
C
C     On Return
C
C        B       the solution vector  X .
C
C     Error Condition
C
C        A division by zero will occur if the input factor contains a
C        zero on the diagonal.  Technically this indicates singularity
C        but it is often caused by improper arguments or improper
C        setting of LDA .  It will not occur if the subroutines are
C        called correctly and if DGBCO has set RCOND .GT. 0.0
C        or DGBFA has set INFO .EQ. 0 .
C
C     To compute  INVERSE(A) * C  where  C  is a matrix
C     with  P  columns
C           CALL DGBCO(ABD,LDA,N,ML,MU,IPVT,RCOND,Z)
C           IF (RCOND is too small) GO TO ...
C           DO 10 J = 1, P
C              CALL DGBSL(ABD,LDA,N,ML,MU,IPVT,C(1,J),0)
C        10 CONTINUE
C
C***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
C                 Stewart, LINPACK Users' Guide, SIAM, 1979.
C***ROUTINES CALLED  DAXPY, DDOT
C***REVISION HISTORY  (YYMMDD)
C   780814  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DGBSL
      INTEGER LDA,N,ML,MU,IPVT(*),JOB
      DOUBLE PRECISION ABD(LDA,*),B(*)
C
      DOUBLE PRECISION DDOT,T
      INTEGER K,KB,L,LA,LB,LM,M,NM1
C***FIRST EXECUTABLE STATEMENT  DGBSL
      M = MU + ML + 1
      NM1 = N - 1
      IF (JOB .NE. 0) GO TO 50
C
C        JOB = 0 , SOLVE  A * X = B
C        FIRST SOLVE L*Y = B
C
         IF (ML .EQ. 0) GO TO 30
         IF (NM1 .LT. 1) GO TO 30
            DO 20 K = 1, NM1
               LM = MIN(ML,N-K)
               L = IPVT(K)
               T = B(L)
               IF (L .EQ. K) GO TO 10
                  B(L) = B(K)
                  B(K) = T
   10          CONTINUE
               CALL DAXPY(LM,T,ABD(M+1,K),1,B(K+1),1)
   20       CONTINUE
   30    CONTINUE
C
C        NOW SOLVE  U*X = Y
C
         DO 40 KB = 1, N
            K = N + 1 - KB
            B(K) = B(K)/ABD(M,K)
            LM = MIN(K,M) - 1
            LA = M - LM
            LB = K - LM
            T = -B(K)
            CALL DAXPY(LM,T,ABD(LA,K),1,B(LB),1)
   40    CONTINUE
      GO TO 100
   50 CONTINUE
C
C        JOB = NONZERO, SOLVE  TRANS(A) * X = B
C        FIRST SOLVE  TRANS(U)*Y = B
C
         DO 60 K = 1, N
            LM = MIN(K,M) - 1
            LA = M - LM
            LB = K - LM
            T = DDOT(LM,ABD(LA,K),1,B(LB),1)
            B(K) = (B(K) - T)/ABD(M,K)
   60    CONTINUE
C
C        NOW SOLVE TRANS(L)*X = Y
C
         IF (ML .EQ. 0) GO TO 90
         IF (NM1 .LT. 1) GO TO 90
            DO 80 KB = 1, NM1
               K = N - KB
               LM = MIN(ML,N-K)
               B(K) = B(K) + DDOT(LM,ABD(M+1,K),1,B(K+1),1)
               L = IPVT(K)
               IF (L .EQ. K) GO TO 70
                  T = B(L)
                  B(L) = B(K)
                  B(K) = T
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
  100 CONTINUE
      RETURN
      END
*DECK DAXPY
      SUBROUTINE DAXPY (N, DA, DX, INCX, DY, INCY)
C***BEGIN PROLOGUE  DAXPY
C***PURPOSE  Compute a constant times a vector plus a vector.
C***CATEGORY  D1A7
C***TYPE      DOUBLE PRECISION (SAXPY-S, DAXPY-D, CAXPY-C)
C***KEYWORDS  BLAS, LINEAR ALGEBRA, TRIAD, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DA  double precision scalar multiplier
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C       DY  double precision vector with N elements
C     INCY  storage spacing between elements of DY
C
C     --Output--
C       DY  double precision result (unchanged if N .LE. 0)
C
C     Overwrite double precision DY with double precision DA*DX + DY.
C     For I = 0 to N-1, replace  DY(LY+I*INCY) with DA*DX(LX+I*INCX) +
C       DY(LY+I*INCY),
C     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
C     defined in a similar way using INCY.
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DAXPY
      DOUBLE PRECISION DX(*), DY(*), DA
C***FIRST EXECUTABLE STATEMENT  DAXPY
      IF (N.LE.0 .OR. DA.EQ.0.0D0) RETURN
      IF (INCX .EQ. INCY) IF (INCX-1) 5,20,60
C
C     Code for unequal or nonpositive increments.
C
    5 IX = 1
      IY = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DY(IY) + DA*DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C     Code for both increments equal to 1.
C
C     Clean-up loop so remaining vector length is a multiple of 4.
C
   20 M = MOD(N,4)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1,M
        DY(I) = DY(I) + DA*DX(I)
   30 CONTINUE
      IF (N .LT. 4) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
        DY(I) = DY(I) + DA*DX(I)
        DY(I+1) = DY(I+1) + DA*DX(I+1)
        DY(I+2) = DY(I+2) + DA*DX(I+2)
        DY(I+3) = DY(I+3) + DA*DX(I+3)
   50 CONTINUE
      RETURN
C
C     Code for equal, positive, non-unit increments.
C
   60 NS = N*INCX
      DO 70 I = 1,NS,INCX
        DY(I) = DA*DX(I) + DY(I)
   70 CONTINUE
      RETURN
      END












*DECK DDOT
      DOUBLE PRECISION FUNCTION DDOT (N, DX, INCX, DY, INCY)
C***BEGIN PROLOGUE  DDOT
C***PURPOSE  Compute the inner product of two vectors.
C***CATEGORY  D1A4
C***TYPE      DOUBLE PRECISION (SDOT-S, DDOT-D, CDOTU-C)
C***KEYWORDS  BLAS, INNER PRODUCT, LINEAR ALGEBRA, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C       DY  double precision vector with N elements
C     INCY  storage spacing between elements of DY
C
C     --Output--
C     DDOT  double precision dot product (zero if N .LE. 0)
C
C     Returns the dot product of double precision DX and DY.
C     DDOT = sum for I = 0 to N-1 of  DX(LX+I*INCX) * DY(LY+I*INCY),
C     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
C     defined in a similar way using INCY.
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DDOT
      DOUBLE PRECISION DX(*), DY(*)
C***FIRST EXECUTABLE STATEMENT  DDOT
      DDOT = 0.0D0
      IF (N .LE. 0) RETURN
      IF (INCX .EQ. INCY) IF (INCX-1) 5,20,60
C
C     Code for unequal or nonpositive increments.
C
    5 IX = 1
      IY = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DDOT = DDOT + DX(IX)*DY(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C     Code for both increments equal to 1.
C
C     Clean-up loop so remaining vector length is a multiple of 5.
C
   20 M = MOD(N,5)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1,M
         DDOT = DDOT + DX(I)*DY(I)
   30 CONTINUE
      IF (N .LT. 5) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
      DDOT = DDOT + DX(I)*DY(I) + DX(I+1)*DY(I+1) + DX(I+2)*DY(I+2) +
     1              DX(I+3)*DY(I+3) + DX(I+4)*DY(I+4)
   50 CONTINUE
      RETURN
C
C     Code for equal, positive, non-unit increments.
C
   60 NS = N*INCX
      DO 70 I = 1,NS,INCX
        DDOT = DDOT + DX(I)*DY(I)
   70 CONTINUE
      RETURN
      END












*DECK XERRWD
      SUBROUTINE XERRWD (MSG, NMES, NERR, LEVEL, NI, I1, I2, NR, R1, R2)
C***BEGIN PROLOGUE  XERRWD
C***SUBSIDIARY
C***PURPOSE  Write error message with values.
C***CATEGORY  R3C
C***TYPE      DOUBLE PRECISION (XERRWV-S, XERRWD-D)
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C***DESCRIPTION
C
C  Subroutines XERRWD, XSETF, XSETUN, and the function routine IXSAV,
C  as given here, constitute a simplified version of the SLATEC error
C  handling package.
C
C  All arguments are input arguments.
C
C  MSG    = The message (character array).
C  NMES   = The length of MSG (number of characters).
C  NERR   = The error number (not used).
C  LEVEL  = The error level..
C           0 or 1 means recoverable (control returns to caller).
C           2 means fatal (run is aborted--see note below).
C  NI     = Number of integers (0, 1, or 2) to be printed with message.
C  I1,I2  = Integers to be printed, depending on NI.
C  NR     = Number of reals (0, 1, or 2) to be printed with message.
C  R1,R2  = Reals to be printed, depending on NR.
C
C  Note..  this routine is machine-dependent and specialized for use
C  in limited context, in the following ways..
C  1. The argument MSG is assumed to be of type CHARACTER, and
C     the message is printed with a format of (1X,A).
C  2. The message is assumed to take only one line.
C     Multi-line messages are generated by repeated calls.
C  3. If LEVEL = 2, control passes to the statement   STOP
C     to abort the run.  This statement may be machine-dependent.
C  4. R1 and R2 are assumed to be in double precision and are printed
C     in D21.13 format.
C
C***ROUTINES CALLED  IXSAV
C***REVISION HISTORY  (YYMMDD)
C   920831  DATE WRITTEN
C   921118  Replaced MFLGSV/LUNSAV by IXSAV. (ACH)
C   930329  Modified prologue to SLATEC format. (FNF)
C   930407  Changed MSG from CHARACTER*1 array to variable. (FNF)
C   930922  Minor cosmetic change. (FNF)
C***END PROLOGUE  XERRWD
C
C*Internal Notes:
C
C For a different default logical unit number, IXSAV (or a subsidiary
C routine that it calls) will need to be modified.
C For a different run-abort command, change the statement following
C statement 100 at the end.
C-----------------------------------------------------------------------
C Subroutines called by XERRWD.. None
C Function routine called by XERRWD.. IXSAV
C-----------------------------------------------------------------------
C**End
C
C  Declare arguments.
C
      DOUBLE PRECISION R1, R2
      INTEGER NMES, NERR, LEVEL, NI, I1, I2, NR
      CHARACTER*(*) MSG
C
C  Declare local variables.
C
      INTEGER LUNIT, IXSAV, MESFLG
C
C  Get logical unit number and message print flag.
C
C***FIRST EXECUTABLE STATEMENT  XERRWD
      LUNIT = IXSAV (1, 0, .FALSE.)
      MESFLG = IXSAV (2, 0, .FALSE.)
      IF (MESFLG .EQ. 0) GO TO 100
C
C  Write the message.
C
      WRITE (LUNIT,10)  MSG
 10   FORMAT(1X,A)
      IF (NI .EQ. 1) WRITE (LUNIT, 20) I1
 20   FORMAT(6X,'In above message,  I1 =',I10)
      IF (NI .EQ. 2) WRITE (LUNIT, 30) I1,I2
 30   FORMAT(6X,'In above message,  I1 =',I10,3X,'I2 =',I10)
      IF (NR .EQ. 1) WRITE (LUNIT, 40) R1
 40   FORMAT(6X,'In above message,  R1 =',D21.13)
      IF (NR .EQ. 2) WRITE (LUNIT, 50) R1,R2
 50   FORMAT(6X,'In above,  R1 =',D21.13,3X,'R2 =',D21.13)
C
C  Abort the run if LEVEL = 2.
C
 100  IF (LEVEL .NE. 2) RETURN
      STOP
C----------------------- End of Subroutine XERRWD ----------------------
      END
*DECK XSETF
      SUBROUTINE XSETF (MFLAG)
C***BEGIN PROLOGUE  XSETF
C***PURPOSE  Reset the error print control flag.
C***CATEGORY  R3A
C***TYPE      ALL (XSETF-A)
C***KEYWORDS  ERROR CONTROL
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C***DESCRIPTION
C
C   XSETF sets the error print control flag to MFLAG:
C      MFLAG=1 means print all messages (the default).
C      MFLAG=0 means no printing.
C
C***SEE ALSO  XERRWD, XERRWV
C***REFERENCES  (NONE)
C***ROUTINES CALLED  IXSAV
C***REVISION HISTORY  (YYMMDD)
C   921118  DATE WRITTEN
C   930329  Added SLATEC format prologue. (FNF)
C   930407  Corrected SEE ALSO section. (FNF)
C   930922  Made user-callable, and other cosmetic changes. (FNF)
C***END PROLOGUE  XSETF
C
C Subroutines called by XSETF.. None
C Function routine called by XSETF.. IXSAV
C-----------------------------------------------------------------------
C**End
      INTEGER MFLAG, JUNK, IXSAV
C
C***FIRST EXECUTABLE STATEMENT  XSETF
      IF (MFLAG .EQ. 0 .OR. MFLAG .EQ. 1) JUNK = IXSAV (2,MFLAG,.TRUE.)
      RETURN
C----------------------- End of Subroutine XSETF -----------------------
      END
*DECK XSETUN
      SUBROUTINE XSETUN (LUN)
C***BEGIN PROLOGUE  XSETUN
C***PURPOSE  Reset the logical unit number for error messages.
C***CATEGORY  R3B
C***TYPE      ALL (XSETUN-A)
C***KEYWORDS  ERROR CONTROL
C***DESCRIPTION
C
C   XSETUN sets the logical unit number for error messages to LUN.
C
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C***SEE ALSO  XERRWD, XERRWV
C***REFERENCES  (NONE)
C***ROUTINES CALLED  IXSAV
C***REVISION HISTORY  (YYMMDD)
C   921118  DATE WRITTEN
C   930329  Added SLATEC format prologue. (FNF)
C   930407  Corrected SEE ALSO section. (FNF)
C   930922  Made user-callable, and other cosmetic changes. (FNF)
C***END PROLOGUE  XSETUN
C
C Subroutines called by XSETUN.. None
C Function routine called by XSETUN.. IXSAV
C-----------------------------------------------------------------------
C**End
      INTEGER LUN, JUNK, IXSAV
C
C***FIRST EXECUTABLE STATEMENT  XSETUN
      IF (LUN .GT. 0) JUNK = IXSAV (1,LUN,.TRUE.)
      RETURN
C----------------------- End of Subroutine XSETUN ----------------------
      END
*DECK IXSAV
      INTEGER FUNCTION IXSAV (IPAR, IVALUE, ISET)
C***BEGIN PROLOGUE  IXSAV
C***SUBSIDIARY
C***PURPOSE  Save and recall error message control parameters.
C***CATEGORY  R3C
C***TYPE      ALL (IXSAV-A)
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C***DESCRIPTION
C
C  IXSAV saves and recalls one of two error message parameters:
C    LUNIT, the logical unit number to which messages are printed, and
C    MESFLG, the message print flag.
C  This is a modification of the SLATEC library routine J4SAVE.
C
C  Saved local variables..
C   LUNIT  = Logical unit number for messages.  The default is obtained
C            by a call to IUMACH (may be machine-dependent).
C   MESFLG = Print control flag..
C            1 means print all messages (the default).
C            0 means no printing.
C
C  On input..
C    IPAR   = Parameter indicator (1 for LUNIT, 2 for MESFLG).
C    IVALUE = The value to be set for the parameter, if ISET = .TRUE.
C    ISET   = Logical flag to indicate whether to read or write.
C             If ISET = .TRUE., the parameter will be given
C             the value IVALUE.  If ISET = .FALSE., the parameter
C             will be unchanged, and IVALUE is a dummy argument.
C
C  On return..
C    IXSAV = The (old) value of the parameter.
C
C***SEE ALSO  XERRWD, XERRWV
C***ROUTINES CALLED  IUMACH
C***REVISION HISTORY  (YYMMDD)
C   921118  DATE WRITTEN
C   930329  Modified prologue to SLATEC format. (FNF)
C   930915  Added IUMACH call to get default output unit.  (ACH)
C   930922  Minor cosmetic changes. (FNF)
C   010425  Type declaration for IUMACH added. (ACH)
C***END PROLOGUE  IXSAV
C
C Subroutines called by IXSAV.. None
C Function routine called by IXSAV.. IUMACH
C-----------------------------------------------------------------------
C**End
      LOGICAL ISET
      INTEGER IPAR, IVALUE
C-----------------------------------------------------------------------
      INTEGER IUMACH, LUNIT, MESFLG
C-----------------------------------------------------------------------
C The following Fortran-77 declaration is to cause the values of the
C listed (local) variables to be saved between calls to this routine.
C-----------------------------------------------------------------------
      SAVE LUNIT, MESFLG
      DATA LUNIT/-1/, MESFLG/1/
C
C***FIRST EXECUTABLE STATEMENT  IXSAV
      IF (IPAR .EQ. 1) THEN
        IF (LUNIT .EQ. -1) LUNIT = IUMACH()
        IXSAV = LUNIT
        IF (ISET) LUNIT = IVALUE
        ENDIF
C
      IF (IPAR .EQ. 2) THEN
        IXSAV = MESFLG
        IF (ISET) MESFLG = IVALUE
        ENDIF
C
      RETURN
C----------------------- End of Function IXSAV -------------------------
      END
*DECK IUMACH
      INTEGER FUNCTION IUMACH()
C***BEGIN PROLOGUE  IUMACH
C***PURPOSE  Provide standard output unit number.
C***CATEGORY  R1
C***TYPE      INTEGER (IUMACH-I)
C***KEYWORDS  MACHINE CONSTANTS
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C***DESCRIPTION
C *Usage:
C        INTEGER  LOUT, IUMACH
C        LOUT = IUMACH()
C
C *Function Return Values:
C     LOUT : the standard logical unit for Fortran output.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   930915  DATE WRITTEN
C   930922  Made user-callable, and other cosmetic changes. (FNF)
C***END PROLOGUE  IUMACH
C
C*Internal Notes:
C  The built-in value of 6 is standard on a wide range of Fortran
C  systems.  This may be machine-dependent.
C**End
C***FIRST EXECUTABLE STATEMENT  IUMACH
      IUMACH = 6
C
      RETURN
C----------------------- End of Function IUMACH ------------------------
      END

c----------------------end of dlsoibt related routines ----------------------

c
c
c
c

C ******************************************
C     VERSION OF SEPTEMBER 18, 1995      
C ******************************************
C
      SUBROUTINE DECOMR(N,FJAC,LDJAC,FMAS,LDMAS,MLMAS,MUMAS,
     &            M1,M2,NM1,FAC1,E1,LDE1,IP1,IER,IJOB,CALHES,IPHES)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION FJAC(LDJAC,N),FMAS(LDMAS,NM1),E1(LDE1,NM1),
     &          IP1(NM1),IPHES(N)
      LOGICAL CALHES
      COMMON/LINAL/MLE,MUE,MBJAC,MBB,MDIAG,MDIFF,MBDIAG
C
      GOTO (1,2,3,4,5,6,7,55,55,55,11,12,13,14,15), IJOB
C
C -----------------------------------------------------------
C
   1  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX
      DO J=1,N
         DO  I=1,N
            E1(I,J)=-FJAC(I,J)
         END DO
         E1(J,J)=E1(J,J)+FAC1
      END DO
      CALL DEC (N,LDE1,E1,IP1,IER)
      RETURN
C
C -----------------------------------------------------------
C
  11  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO J=1,NM1
         JM1=J+M1
         DO I=1,NM1
            E1(I,J)=-FJAC(I,JM1)
         END DO
         E1(J,J)=E1(J,J)+FAC1
      END DO
 45   MM=M1/M2
      DO J=1,M2
         DO I=1,NM1
            SUM=0.D0
            DO K=0,MM-1
               SUM=(SUM+FJAC(I,J+K*M2))/FAC1
            END DO
            E1(I,J)=E1(I,J)-SUM
         END DO
      END DO
      CALL DEC (NM1,LDE1,E1,IP1,IER)
      RETURN
C
C -----------------------------------------------------------
C
   2  CONTINUE
C ---  B=IDENTITY, JACOBIAN A BANDED MATRIX
      DO J=1,N
         DO I=1,MBJAC
            E1(I+MLE,J)=-FJAC(I,J)
         END DO
         E1(MDIAG,J)=E1(MDIAG,J)+FAC1
      END DO
      CALL DECB (N,LDE1,E1,MLE,MUE,IP1,IER)
      RETURN
C
C -----------------------------------------------------------
C
  12  CONTINUE
C ---  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER
      DO J=1,NM1
         JM1=J+M1
         DO I=1,MBJAC
            E1(I+MLE,J)=-FJAC(I,JM1)
         END DO
         E1(MDIAG,J)=E1(MDIAG,J)+FAC1
      END DO
  46  MM=M1/M2
      DO J=1,M2
         DO I=1,MBJAC
            SUM=0.D0
            DO K=0,MM-1
               SUM=(SUM+FJAC(I,J+K*M2))/FAC1
            END DO
            E1(I+MLE,J)=E1(I+MLE,J)-SUM
         END DO
      END DO
      CALL DECB (NM1,LDE1,E1,MLE,MUE,IP1,IER)
      RETURN
C
C -----------------------------------------------------------
C
   3  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX
      DO J=1,N
         DO I=1,N
            E1(I,J)=-FJAC(I,J)
         END DO
         DO I=MAX(1,J-MUMAS),MIN(N,J+MLMAS)
            E1(I,J)=E1(I,J)+FAC1*FMAS(I-J+MBDIAG,J)
         END DO
      END DO
      CALL DEC (N,LDE1,E1,IP1,IER)
      RETURN
C
C -----------------------------------------------------------
C
  13  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO J=1,NM1
         JM1=J+M1
         DO I=1,NM1
            E1(I,J)=-FJAC(I,JM1)
         END DO
         DO I=MAX(1,J-MUMAS),MIN(NM1,J+MLMAS)
            E1(I,J)=E1(I,J)+FAC1*FMAS(I-J+MBDIAG,J)
         END DO
      END DO
      GOTO 45
C
C -----------------------------------------------------------
C
   4  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX
      DO J=1,N
         DO I=1,MBJAC
            E1(I+MLE,J)=-FJAC(I,J)
         END DO
         DO I=1,MBB
            IB=I+MDIFF
            E1(IB,J)=E1(IB,J)+FAC1*FMAS(I,J)
         END DO
      END DO
      CALL DECB (N,LDE1,E1,MLE,MUE,IP1,IER)
      RETURN
C
C -----------------------------------------------------------
C
  14  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX, SECOND ORDER
      DO J=1,NM1
         JM1=J+M1
         DO I=1,MBJAC
            E1(I+MLE,J)=-FJAC(I,JM1)
         END DO
         DO I=1,MBB
            IB=I+MDIFF
            E1(IB,J)=E1(IB,J)+FAC1*FMAS(I,J)
         END DO
      END DO
      GOTO 46
C
C -----------------------------------------------------------
C
   5  CONTINUE
C ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX
      DO J=1,N
         DO I=1,N
            E1(I,J)=FMAS(I,J)*FAC1-FJAC(I,J)
         END DO
      END DO
      CALL DEC (N,LDE1,E1,IP1,IER)
      RETURN
C
C -----------------------------------------------------------
C
  15  CONTINUE
C ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO J=1,NM1
         JM1=J+M1
         DO I=1,NM1
            E1(I,J)=FMAS(I,J)*FAC1-FJAC(I,JM1)
         END DO
      END DO
      GOTO 45
C
C -----------------------------------------------------------
C
   6  CONTINUE
C ---  B IS A FULL MATRIX, JACOBIAN A BANDED MATRIX
C ---  THIS OPTION IS NOT PROVIDED
      RETURN
C
C -----------------------------------------------------------
C
   7  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX, HESSENBERG-OPTION
      IF (CALHES) CALL ELMHES (LDJAC,N,1,N,FJAC,IPHES) 
      CALHES=.FALSE.
      DO J=1,N-1
         J1=J+1
         E1(J1,J)=-FJAC(J1,J)
      END DO
      DO J=1,N
         DO I=1,J
            E1(I,J)=-FJAC(I,J)
         END DO
         E1(J,J)=E1(J,J)+FAC1
      END DO
      CALL DECH(N,LDE1,E1,1,IP1,IER)
      RETURN
C
C -----------------------------------------------------------
C
  55  CONTINUE
      RETURN
      END
C
C     END OF SUBROUTINE DECOMR
C
C ***********************************************************
C
      SUBROUTINE DECOMC(N,FJAC,LDJAC,FMAS,LDMAS,MLMAS,MUMAS,
     &            M1,M2,NM1,ALPHN,BETAN,E2R,E2I,LDE1,IP2,IER,IJOB)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION FJAC(LDJAC,N),FMAS(LDMAS,NM1),
     &          E2R(LDE1,NM1),E2I(LDE1,NM1),IP2(NM1)
      COMMON/LINAL/MLE,MUE,MBJAC,MBB,MDIAG,MDIFF,MBDIAG
C
      GOTO (1,2,3,4,5,6,7,55,55,55,11,12,13,14,15), IJOB
C
C -----------------------------------------------------------
C
   1  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX
      DO J=1,N
         DO I=1,N
            E2R(I,J)=-FJAC(I,J)
            E2I(I,J)=0.D0
         END DO
         E2R(J,J)=E2R(J,J)+ALPHN
         E2I(J,J)=BETAN
      END DO
      CALL DECC (N,LDE1,E2R,E2I,IP2,IER)
      RETURN
C
C -----------------------------------------------------------
C
  11  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO J=1,NM1
         JM1=J+M1
         DO I=1,NM1
            E2R(I,J)=-FJAC(I,JM1)
            E2I(I,J)=0.D0
         END DO
         E2R(J,J)=E2R(J,J)+ALPHN
         E2I(J,J)=BETAN
      END DO
  45  MM=M1/M2
      ABNO=ALPHN**2+BETAN**2
      ALP=ALPHN/ABNO
      BET=BETAN/ABNO
      DO J=1,M2
         DO I=1,NM1
            SUMR=0.D0
            SUMI=0.D0
            DO K=0,MM-1
               SUMS=SUMR+FJAC(I,J+K*M2)
               SUMR=SUMS*ALP+SUMI*BET
               SUMI=SUMI*ALP-SUMS*BET
            END DO
            E2R(I,J)=E2R(I,J)-SUMR
            E2I(I,J)=E2I(I,J)-SUMI
         END DO
      END DO
      CALL DECC (NM1,LDE1,E2R,E2I,IP2,IER)
      RETURN
C
C -----------------------------------------------------------
C
   2  CONTINUE
C ---  B=IDENTITY, JACOBIAN A BANDED MATRIX
      DO J=1,N
         DO I=1,MBJAC
            IMLE=I+MLE
            E2R(IMLE,J)=-FJAC(I,J)
            E2I(IMLE,J)=0.D0
         END DO
         E2R(MDIAG,J)=E2R(MDIAG,J)+ALPHN
         E2I(MDIAG,J)=BETAN
      END DO
      CALL DECBC (N,LDE1,E2R,E2I,MLE,MUE,IP2,IER)
      RETURN
C
C -----------------------------------------------------------
C
  12  CONTINUE
C ---  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER
      DO J=1,NM1
         JM1=J+M1
         DO I=1,MBJAC
            E2R(I+MLE,J)=-FJAC(I,JM1)
            E2I(I+MLE,J)=0.D0
         END DO
         E2R(MDIAG,J)=E2R(MDIAG,J)+ALPHN
         E2I(MDIAG,J)=E2I(MDIAG,J)+BETAN
      END DO
  46  MM=M1/M2
      ABNO=ALPHN**2+BETAN**2
      ALP=ALPHN/ABNO
      BET=BETAN/ABNO
      DO J=1,M2
         DO I=1,MBJAC
            SUMR=0.D0
            SUMI=0.D0
            DO K=0,MM-1
               SUMS=SUMR+FJAC(I,J+K*M2)
               SUMR=SUMS*ALP+SUMI*BET
               SUMI=SUMI*ALP-SUMS*BET
            END DO
            IMLE=I+MLE
            E2R(IMLE,J)=E2R(IMLE,J)-SUMR
            E2I(IMLE,J)=E2I(IMLE,J)-SUMI
         END DO
      END DO
      CALL DECBC (NM1,LDE1,E2R,E2I,MLE,MUE,IP2,IER)
      RETURN
C
C -----------------------------------------------------------
C
   3  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX
      DO  J=1,N
         DO  I=1,N
            E2R(I,J)=-FJAC(I,J)
            E2I(I,J)=0.D0
         END DO
      END DO
      DO J=1,N
         DO I=MAX(1,J-MUMAS),MIN(N,J+MLMAS)
            BB=FMAS(I-J+MBDIAG,J)
            E2R(I,J)=E2R(I,J)+ALPHN*BB
            E2I(I,J)=BETAN*BB
         END DO
      END DO
      CALL DECC(N,LDE1,E2R,E2I,IP2,IER)
      RETURN
C
C -----------------------------------------------------------
C
  13  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO J=1,NM1
         JM1=J+M1
         DO I=1,NM1
            E2R(I,J)=-FJAC(I,JM1)
            E2I(I,J)=0.D0
         END DO
         DO I=MAX(1,J-MUMAS),MIN(NM1,J+MLMAS)
            FFMA=FMAS(I-J+MBDIAG,J)
            E2R(I,J)=E2R(I,J)+ALPHN*FFMA
            E2I(I,J)=E2I(I,J)+BETAN*FFMA
         END DO
      END DO
      GOTO 45
C
C -----------------------------------------------------------
C
   4  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX
      DO J=1,N
         DO I=1,MBJAC
            IMLE=I+MLE
            E2R(IMLE,J)=-FJAC(I,J)
            E2I(IMLE,J)=0.D0
         END DO
         DO I=MAX(1,MUMAS+2-J),MIN(MBB,MUMAS+1-J+N)
            IB=I+MDIFF
            BB=FMAS(I,J)
            E2R(IB,J)=E2R(IB,J)+ALPHN*BB
            E2I(IB,J)=BETAN*BB
         END DO
      END DO
      CALL DECBC (N,LDE1,E2R,E2I,MLE,MUE,IP2,IER)
      RETURN
C
C -----------------------------------------------------------
C
  14  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX, SECOND ORDER
      DO J=1,NM1
         JM1=J+M1
         DO I=1,MBJAC
            E2R(I+MLE,J)=-FJAC(I,JM1)
            E2I(I+MLE,J)=0.D0
         END DO
         DO I=1,MBB
            IB=I+MDIFF
            FFMA=FMAS(I,J)
            E2R(IB,J)=E2R(IB,J)+ALPHN*FFMA
            E2I(IB,J)=E2I(IB,J)+BETAN*FFMA
         END DO
      END DO
      GOTO 46
C
C -----------------------------------------------------------
C
   5  CONTINUE
C ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX
      DO J=1,N
         DO I=1,N
            BB=FMAS(I,J)
            E2R(I,J)=BB*ALPHN-FJAC(I,J)
            E2I(I,J)=BB*BETAN
         END DO
      END DO
      CALL DECC(N,LDE1,E2R,E2I,IP2,IER)
      RETURN
C
C -----------------------------------------------------------
C
  15  CONTINUE
C ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO J=1,NM1
         JM1=J+M1
         DO I=1,NM1
            E2R(I,J)=ALPHN*FMAS(I,J)-FJAC(I,JM1)
            E2I(I,J)=BETAN*FMAS(I,J)
         END DO
      END DO
      GOTO 45
C
C -----------------------------------------------------------
C
   6  CONTINUE
C ---  B IS A FULL MATRIX, JACOBIAN A BANDED MATRIX
C ---  THIS OPTION IS NOT PROVIDED
      RETURN
C
C -----------------------------------------------------------
C
   7  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX, HESSENBERG-OPTION
      DO J=1,N-1
         J1=J+1
         E2R(J1,J)=-FJAC(J1,J)
         E2I(J1,J)=0.D0
      END DO
      DO J=1,N
         DO I=1,J
            E2I(I,J)=0.D0
            E2R(I,J)=-FJAC(I,J)
         END DO
         E2R(J,J)=E2R(J,J)+ALPHN
         E2I(J,J)=BETAN
      END DO
      CALL DECHC(N,LDE1,E2R,E2I,1,IP2,IER)
      RETURN
C
C -----------------------------------------------------------
C
  55  CONTINUE
      RETURN
      END
C
C     END OF SUBROUTINE DECOMC
C
C ***********************************************************
C
      SUBROUTINE SLVRAR(N,FJAC,LDJAC,MLJAC,MUJAC,FMAS,LDMAS,MLMAS,MUMAS,
     &          M1,M2,NM1,FAC1,E1,LDE1,Z1,F1,IP1,IPHES,IER,IJOB)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION FJAC(LDJAC,N),FMAS(LDMAS,NM1),E1(LDE1,NM1),
     &          IP1(NM1),IPHES(N),Z1(N),F1(N)
      COMMON/LINAL/MLE,MUE,MBJAC,MBB,MDIAG,MDIFF,MBDIAG
C
      GOTO (1,2,3,4,5,6,7,55,55,55,11,12,13,13,15), IJOB
C
C -----------------------------------------------------------
C
   1  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX
      DO I=1,N
         Z1(I)=Z1(I)-F1(I)*FAC1
      END DO
      CALL SOL (N,LDE1,E1,Z1,IP1)
      RETURN
C
C -----------------------------------------------------------
C
  11  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO I=1,N
         Z1(I)=Z1(I)-F1(I)*FAC1
      END DO
 48   CONTINUE
      MM=M1/M2
      DO J=1,M2
         SUM1=0.D0
         DO K=MM-1,0,-1
            JKM=J+K*M2
            SUM1=(Z1(JKM)+SUM1)/FAC1
            DO I=1,NM1
               IM1=I+M1
               Z1(IM1)=Z1(IM1)+FJAC(I,JKM)*SUM1
            END DO
         END DO
      END DO
      CALL SOL (NM1,LDE1,E1,Z1(M1+1),IP1)
 49   CONTINUE
      DO I=M1,1,-1
         Z1(I)=(Z1(I)+Z1(M2+I))/FAC1
      END DO
      RETURN
C
C -----------------------------------------------------------
C
   2  CONTINUE
C ---  B=IDENTITY, JACOBIAN A BANDED MATRIX
      DO I=1,N
         Z1(I)=Z1(I)-F1(I)*FAC1
      END DO
      CALL SOLB (N,LDE1,E1,MLE,MUE,Z1,IP1)
      RETURN
C
C -----------------------------------------------------------
C
  12  CONTINUE
C ---  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER
      DO I=1,N
         Z1(I)=Z1(I)-F1(I)*FAC1
      END DO
  45  CONTINUE
      MM=M1/M2
      DO J=1,M2
         SUM1=0.D0
         DO K=MM-1,0,-1
            JKM=J+K*M2
            SUM1=(Z1(JKM)+SUM1)/FAC1
            DO I=MAX(1,J-MUJAC),MIN(NM1,J+MLJAC)
               IM1=I+M1
               Z1(IM1)=Z1(IM1)+FJAC(I+MUJAC+1-J,JKM)*SUM1
            END DO
         END DO
      END DO
      CALL SOLB (NM1,LDE1,E1,MLE,MUE,Z1(M1+1),IP1)
      GOTO 49
C
C -----------------------------------------------------------
C
   3  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX
      DO I=1,N
         S1=0.0D0
         DO J=MAX(1,I-MLMAS),MIN(N,I+MUMAS)
            S1=S1-FMAS(I-J+MBDIAG,J)*F1(J)
         END DO
         Z1(I)=Z1(I)+S1*FAC1
      END DO
      CALL SOL (N,LDE1,E1,Z1,IP1)
      RETURN
C
C -----------------------------------------------------------
C
  13  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO I=1,M1
         Z1(I)=Z1(I)-F1(I)*FAC1
      END DO
      DO I=1,NM1
         IM1=I+M1
         S1=0.0D0
         DO J=MAX(1,I-MLMAS),MIN(NM1,I+MUMAS)
            S1=S1-FMAS(I-J+MBDIAG,J)*F1(J+M1)
         END DO
         Z1(IM1)=Z1(IM1)+S1*FAC1
      END DO
      IF (IJOB.EQ.14) GOTO 45
      GOTO 48
C
C -----------------------------------------------------------
C
   4  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX
      DO I=1,N
         S1=0.0D0
         DO J=MAX(1,I-MLMAS),MIN(N,I+MUMAS)
            S1=S1-FMAS(I-J+MBDIAG,J)*F1(J)
         END DO
         Z1(I)=Z1(I)+S1*FAC1
      END DO
      CALL SOLB (N,LDE1,E1,MLE,MUE,Z1,IP1)
      RETURN
C
C -----------------------------------------------------------
C
   5  CONTINUE
C ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX
      DO I=1,N
         S1=0.0D0
         DO J=1,N
            S1=S1-FMAS(I,J)*F1(J)
         END DO
         Z1(I)=Z1(I)+S1*FAC1
      END DO
      CALL SOL (N,LDE1,E1,Z1,IP1)
      RETURN
C
C -----------------------------------------------------------
C
  15  CONTINUE
C ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO I=1,M1
         Z1(I)=Z1(I)-F1(I)*FAC1
      END DO
      DO I=1,NM1
         IM1=I+M1
         S1=0.0D0
         DO J=1,NM1
            S1=S1-FMAS(I,J)*F1(J+M1)
         END DO
         Z1(IM1)=Z1(IM1)+S1*FAC1
      END DO
      GOTO 48
C
C -----------------------------------------------------------
C
   6  CONTINUE
C ---  B IS A FULL MATRIX, JACOBIAN A BANDED MATRIX
C ---  THIS OPTION IS NOT PROVIDED
      RETURN
C
C -----------------------------------------------------------
C
   7  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX, HESSENBERG-OPTION
      DO I=1,N
         Z1(I)=Z1(I)-F1(I)*FAC1
      END DO
      DO MM=N-2,1,-1
          MP=N-MM
          MP1=MP-1
          I=IPHES(MP)
          IF (I.EQ.MP) GOTO 746
          ZSAFE=Z1(MP)
          Z1(MP)=Z1(I)
          Z1(I)=ZSAFE
 746      CONTINUE
          DO I=MP+1,N 
             Z1(I)=Z1(I)-FJAC(I,MP1)*Z1(MP)
          END DO
       END DO
       CALL SOLH(N,LDE1,E1,1,Z1,IP1)
       DO MM=1,N-2
          MP=N-MM
          MP1=MP-1
          DO I=MP+1,N 
             Z1(I)=Z1(I)+FJAC(I,MP1)*Z1(MP)
          END DO
          I=IPHES(MP)
          IF (I.EQ.MP) GOTO 750
          ZSAFE=Z1(MP)
          Z1(MP)=Z1(I)
          Z1(I)=ZSAFE
 750      CONTINUE
      END DO
      RETURN
C
C -----------------------------------------------------------
C
  55  CONTINUE
      RETURN
      END
C
C     END OF SUBROUTINE SLVRAR
C
C ***********************************************************
C
      SUBROUTINE SLVRAI(N,FJAC,LDJAC,MLJAC,MUJAC,FMAS,LDMAS,MLMAS,MUMAS,
     &          M1,M2,NM1,ALPHN,BETAN,E2R,E2I,LDE1,Z2,Z3,
     &          F2,F3,CONT,IP2,IPHES,IER,IJOB)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION FJAC(LDJAC,N),FMAS(LDMAS,NM1),
     &          IP2(NM1),IPHES(N),Z2(N),Z3(N),F2(N),F3(N)
      DIMENSION E2R(LDE1,NM1),E2I(LDE1,NM1)
      COMMON/LINAL/MLE,MUE,MBJAC,MBB,MDIAG,MDIFF,MBDIAG
C
      GOTO (1,2,3,4,5,6,7,55,55,55,11,12,13,13,15), IJOB
C
C -----------------------------------------------------------
C
   1  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX
      DO I=1,N
         S2=-F2(I)
         S3=-F3(I)
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
      CALL SOLC (N,LDE1,E2R,E2I,Z2,Z3,IP2)
      RETURN
C
C -----------------------------------------------------------
C
  11  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO I=1,N
         S2=-F2(I)
         S3=-F3(I)
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
 48   ABNO=ALPHN**2+BETAN**2
      MM=M1/M2
      DO J=1,M2
         SUM2=0.D0
         SUM3=0.D0
         DO K=MM-1,0,-1
            JKM=J+K*M2
            SUMH=(Z2(JKM)+SUM2)/ABNO
            SUM3=(Z3(JKM)+SUM3)/ABNO
            SUM2=SUMH*ALPHN+SUM3*BETAN
            SUM3=SUM3*ALPHN-SUMH*BETAN
            DO I=1,NM1
               IM1=I+M1
               Z2(IM1)=Z2(IM1)+FJAC(I,JKM)*SUM2
               Z3(IM1)=Z3(IM1)+FJAC(I,JKM)*SUM3
            END DO
         END DO
      END DO
      CALL SOLC (NM1,LDE1,E2R,E2I,Z2(M1+1),Z3(M1+1),IP2)
 49   CONTINUE
      DO I=M1,1,-1
         MPI=M2+I
         Z2I=Z2(I)+Z2(MPI)
         Z3I=Z3(I)+Z3(MPI)
         Z3(I)=(Z3I*ALPHN-Z2I*BETAN)/ABNO
         Z2(I)=(Z2I*ALPHN+Z3I*BETAN)/ABNO
      END DO
      RETURN
C
C -----------------------------------------------------------
C
   2  CONTINUE
C ---  B=IDENTITY, JACOBIAN A BANDED MATRIX
      DO I=1,N
         S2=-F2(I)
         S3=-F3(I)
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
      CALL SOLBC (N,LDE1,E2R,E2I,MLE,MUE,Z2,Z3,IP2)
      RETURN
C
C -----------------------------------------------------------
C
  12  CONTINUE
C ---  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER
      DO I=1,N
         S2=-F2(I)
         S3=-F3(I)
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
  45  ABNO=ALPHN**2+BETAN**2
      MM=M1/M2
      DO J=1,M2
         SUM2=0.D0
         SUM3=0.D0
         DO K=MM-1,0,-1
            JKM=J+K*M2
            SUMH=(Z2(JKM)+SUM2)/ABNO
            SUM3=(Z3(JKM)+SUM3)/ABNO
            SUM2=SUMH*ALPHN+SUM3*BETAN
            SUM3=SUM3*ALPHN-SUMH*BETAN
            DO I=MAX(1,J-MUJAC),MIN(NM1,J+MLJAC)
               IM1=I+M1
               IIMU=I+MUJAC+1-J
               Z2(IM1)=Z2(IM1)+FJAC(IIMU,JKM)*SUM2
               Z3(IM1)=Z3(IM1)+FJAC(IIMU,JKM)*SUM3
            END DO
         END DO
      END DO
      CALL SOLBC (NM1,LDE1,E2R,E2I,MLE,MUE,Z2(M1+1),Z3(M1+1),IP2)
      GOTO 49
C
C -----------------------------------------------------------
C
   3  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX
      DO I=1,N
         S2=0.0D0
         S3=0.0D0
         DO J=MAX(1,I-MLMAS),MIN(N,I+MUMAS)
            BB=FMAS(I-J+MBDIAG,J)
            S2=S2-BB*F2(J)
            S3=S3-BB*F3(J)
         END DO
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
      CALL SOLC(N,LDE1,E2R,E2I,Z2,Z3,IP2)
      RETURN
C
C -----------------------------------------------------------
C
  13  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO I=1,M1
         S2=-F2(I)
         S3=-F3(I)
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
      DO I=1,NM1
         IM1=I+M1
         S2=0.0D0
         S3=0.0D0
         DO J=MAX(1,I-MLMAS),MIN(NM1,I+MUMAS)
            JM1=J+M1
            BB=FMAS(I-J+MBDIAG,J)
            S2=S2-BB*F2(JM1)
            S3=S3-BB*F3(JM1)
         END DO
         Z2(IM1)=Z2(IM1)+S2*ALPHN-S3*BETAN
         Z3(IM1)=Z3(IM1)+S3*ALPHN+S2*BETAN
      END DO
      IF (IJOB.EQ.14) GOTO 45
      GOTO 48
C
C -----------------------------------------------------------
C
   4  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX
      DO I=1,N
         S2=0.0D0
         S3=0.0D0
         DO J=MAX(1,I-MLMAS),MIN(N,I+MUMAS)
            BB=FMAS(I-J+MBDIAG,J)
            S2=S2-BB*F2(J)
            S3=S3-BB*F3(J)
         END DO
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
      CALL SOLBC(N,LDE1,E2R,E2I,MLE,MUE,Z2,Z3,IP2)
      RETURN
C
C -----------------------------------------------------------
C
   5  CONTINUE
C ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX
      DO I=1,N
         S2=0.0D0
         S3=0.0D0
         DO J=1,N
            BB=FMAS(I,J)
            S2=S2-BB*F2(J)
            S3=S3-BB*F3(J)
         END DO
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
      CALL SOLC(N,LDE1,E2R,E2I,Z2,Z3,IP2)
      RETURN
C
C -----------------------------------------------------------
C
  15  CONTINUE
C ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO I=1,M1
         S2=-F2(I)
         S3=-F3(I)
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
      DO I=1,NM1
         IM1=I+M1
         S2=0.0D0
         S3=0.0D0
         DO J=1,NM1
            JM1=J+M1
            BB=FMAS(I,J)
            S2=S2-BB*F2(JM1)
            S3=S3-BB*F3(JM1)
         END DO
         Z2(IM1)=Z2(IM1)+S2*ALPHN-S3*BETAN
         Z3(IM1)=Z3(IM1)+S3*ALPHN+S2*BETAN
      END DO
      GOTO 48
C
C -----------------------------------------------------------
C
   6  CONTINUE
C ---  B IS A FULL MATRIX, JACOBIAN A BANDED MATRIX
C ---  THIS OPTION IS NOT PROVIDED
      RETURN
C
C -----------------------------------------------------------
C
   7  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX, HESSENBERG-OPTION
      DO I=1,N
         S2=-F2(I)
         S3=-F3(I)
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
      DO MM=N-2,1,-1
          MP=N-MM
          MP1=MP-1
          I=IPHES(MP)
          IF (I.EQ.MP) GOTO 746
          ZSAFE=Z2(MP)
          Z2(MP)=Z2(I)
          Z2(I)=ZSAFE 
          ZSAFE=Z3(MP)
          Z3(MP)=Z3(I)
          Z3(I)=ZSAFE
 746      CONTINUE
          DO I=MP+1,N 
             E1IMP=FJAC(I,MP1)
             Z2(I)=Z2(I)-E1IMP*Z2(MP)
             Z3(I)=Z3(I)-E1IMP*Z3(MP)
          END DO
       END DO
       CALL SOLHC(N,LDE1,E2R,E2I,1,Z2,Z3,IP2)
       DO MM=1,N-2
          MP=N-MM
          MP1=MP-1
          DO I=MP+1,N 
             E1IMP=FJAC(I,MP1)
             Z2(I)=Z2(I)+E1IMP*Z2(MP)
             Z3(I)=Z3(I)+E1IMP*Z3(MP)
          END DO
          I=IPHES(MP)
          IF (I.EQ.MP) GOTO 750
          ZSAFE=Z2(MP)
          Z2(MP)=Z2(I)
          Z2(I)=ZSAFE 
          ZSAFE=Z3(MP)
          Z3(MP)=Z3(I)
          Z3(I)=ZSAFE
 750      CONTINUE
      END DO
      RETURN
C
C -----------------------------------------------------------
C
  55  CONTINUE
      RETURN
      END
C
C     END OF SUBROUTINE SLVRAI
C
C ***********************************************************
C
      SUBROUTINE SLVRAD(N,FJAC,LDJAC,MLJAC,MUJAC,FMAS,LDMAS,MLMAS,MUMAS,
     &          M1,M2,NM1,FAC1,ALPHN,BETAN,E1,E2R,E2I,LDE1,Z1,Z2,Z3,
     &          F1,F2,F3,CONT,IP1,IP2,IPHES,IER,IJOB)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION FJAC(LDJAC,N),FMAS(LDMAS,NM1),E1(LDE1,NM1),
     &          E2R(LDE1,NM1),E2I(LDE1,NM1),IP1(NM1),IP2(NM1),
     &          IPHES(N),Z1(N),Z2(N),Z3(N),F1(N),F2(N),F3(N)
      COMMON/LINAL/MLE,MUE,MBJAC,MBB,MDIAG,MDIFF,MBDIAG
C
      GOTO (1,2,3,4,5,6,7,55,55,55,11,12,13,13,15), IJOB
C
C -----------------------------------------------------------
C
   1  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX
      DO I=1,N
         S2=-F2(I)
         S3=-F3(I)
         Z1(I)=Z1(I)-F1(I)*FAC1
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
      CALL SOL (N,LDE1,E1,Z1,IP1)
      CALL SOLC (N,LDE1,E2R,E2I,Z2,Z3,IP2)
      RETURN
C
C -----------------------------------------------------------
C
  11  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO I=1,N
         S2=-F2(I)
         S3=-F3(I)
         Z1(I)=Z1(I)-F1(I)*FAC1
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
 48   ABNO=ALPHN**2+BETAN**2
      MM=M1/M2
      DO J=1,M2
         SUM1=0.D0
         SUM2=0.D0
         SUM3=0.D0
         DO K=MM-1,0,-1
            JKM=J+K*M2
            SUM1=(Z1(JKM)+SUM1)/FAC1
            SUMH=(Z2(JKM)+SUM2)/ABNO
            SUM3=(Z3(JKM)+SUM3)/ABNO
            SUM2=SUMH*ALPHN+SUM3*BETAN
            SUM3=SUM3*ALPHN-SUMH*BETAN
            DO I=1,NM1
               IM1=I+M1
               Z1(IM1)=Z1(IM1)+FJAC(I,JKM)*SUM1
               Z2(IM1)=Z2(IM1)+FJAC(I,JKM)*SUM2
               Z3(IM1)=Z3(IM1)+FJAC(I,JKM)*SUM3
            END DO
         END DO
      END DO
      CALL SOL (NM1,LDE1,E1,Z1(M1+1),IP1)
      CALL SOLC (NM1,LDE1,E2R,E2I,Z2(M1+1),Z3(M1+1),IP2)
 49   CONTINUE
      DO I=M1,1,-1
         MPI=M2+I
         Z1(I)=(Z1(I)+Z1(MPI))/FAC1
         Z2I=Z2(I)+Z2(MPI)
         Z3I=Z3(I)+Z3(MPI)
         Z3(I)=(Z3I*ALPHN-Z2I*BETAN)/ABNO
         Z2(I)=(Z2I*ALPHN+Z3I*BETAN)/ABNO
      END DO
      RETURN
C
C -----------------------------------------------------------
C
   2  CONTINUE
C ---  B=IDENTITY, JACOBIAN A BANDED MATRIX
      DO I=1,N
         S2=-F2(I)
         S3=-F3(I)
         Z1(I)=Z1(I)-F1(I)*FAC1
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
      CALL SOLB (N,LDE1,E1,MLE,MUE,Z1,IP1)
      CALL SOLBC (N,LDE1,E2R,E2I,MLE,MUE,Z2,Z3,IP2)
      RETURN
C
C -----------------------------------------------------------
C
  12  CONTINUE
C ---  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER
      DO I=1,N
         S2=-F2(I)
         S3=-F3(I)
         Z1(I)=Z1(I)-F1(I)*FAC1
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
  45  ABNO=ALPHN**2+BETAN**2
      MM=M1/M2
      DO J=1,M2
         SUM1=0.D0
         SUM2=0.D0
         SUM3=0.D0
         DO K=MM-1,0,-1
            JKM=J+K*M2
            SUM1=(Z1(JKM)+SUM1)/FAC1
            SUMH=(Z2(JKM)+SUM2)/ABNO
            SUM3=(Z3(JKM)+SUM3)/ABNO
            SUM2=SUMH*ALPHN+SUM3*BETAN
            SUM3=SUM3*ALPHN-SUMH*BETAN
            DO I=MAX(1,J-MUJAC),MIN(NM1,J+MLJAC)
               IM1=I+M1
               FFJA=FJAC(I+MUJAC+1-J,JKM)
               Z1(IM1)=Z1(IM1)+FFJA*SUM1
               Z2(IM1)=Z2(IM1)+FFJA*SUM2
               Z3(IM1)=Z3(IM1)+FFJA*SUM3
            END DO
         END DO
      END DO
      CALL SOLB (NM1,LDE1,E1,MLE,MUE,Z1(M1+1),IP1)
      CALL SOLBC (NM1,LDE1,E2R,E2I,MLE,MUE,Z2(M1+1),Z3(M1+1),IP2)
      GOTO 49
C
C -----------------------------------------------------------
C
   3  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX
      DO I=1,N
         S1=0.0D0
         S2=0.0D0
         S3=0.0D0
         DO J=MAX(1,I-MLMAS),MIN(N,I+MUMAS)
            BB=FMAS(I-J+MBDIAG,J)
            S1=S1-BB*F1(J)
            S2=S2-BB*F2(J)
            S3=S3-BB*F3(J)
         END DO
         Z1(I)=Z1(I)+S1*FAC1
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
      CALL SOL (N,LDE1,E1,Z1,IP1)
      CALL SOLC(N,LDE1,E2R,E2I,Z2,Z3,IP2)
      RETURN
C
C -----------------------------------------------------------
C
  13  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO I=1,M1
         S2=-F2(I)
         S3=-F3(I)
         Z1(I)=Z1(I)-F1(I)*FAC1
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
      DO I=1,NM1
         IM1=I+M1
         S1=0.0D0
         S2=0.0D0
         S3=0.0D0
         J1B=MAX(1,I-MLMAS)
         J2B=MIN(NM1,I+MUMAS)
         DO J=J1B,J2B
            JM1=J+M1
            BB=FMAS(I-J+MBDIAG,J)
            S1=S1-BB*F1(JM1)
            S2=S2-BB*F2(JM1)
            S3=S3-BB*F3(JM1)
         END DO
         Z1(IM1)=Z1(IM1)+S1*FAC1
         Z2(IM1)=Z2(IM1)+S2*ALPHN-S3*BETAN
         Z3(IM1)=Z3(IM1)+S3*ALPHN+S2*BETAN
      END DO
      IF (IJOB.EQ.14) GOTO 45
      GOTO 48
C
C -----------------------------------------------------------
C
   4  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX
      DO I=1,N
         S1=0.0D0
         S2=0.0D0
         S3=0.0D0
         DO J=MAX(1,I-MLMAS),MIN(N,I+MUMAS)
            BB=FMAS(I-J+MBDIAG,J)
            S1=S1-BB*F1(J)
            S2=S2-BB*F2(J)
            S3=S3-BB*F3(J)
         END DO
         Z1(I)=Z1(I)+S1*FAC1
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
      CALL SOLB (N,LDE1,E1,MLE,MUE,Z1,IP1)
      CALL SOLBC(N,LDE1,E2R,E2I,MLE,MUE,Z2,Z3,IP2)
      RETURN
C
C -----------------------------------------------------------
C
   5  CONTINUE
C ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX
      DO I=1,N
         S1=0.0D0
         S2=0.0D0
         S3=0.0D0
         DO J=1,N
            BB=FMAS(I,J)
            S1=S1-BB*F1(J)
            S2=S2-BB*F2(J)
            S3=S3-BB*F3(J)
         END DO
         Z1(I)=Z1(I)+S1*FAC1
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
      CALL SOL (N,LDE1,E1,Z1,IP1)
      CALL SOLC(N,LDE1,E2R,E2I,Z2,Z3,IP2)
      RETURN
C
C -----------------------------------------------------------
C
  15  CONTINUE
C ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO I=1,M1
         S2=-F2(I)
         S3=-F3(I)
         Z1(I)=Z1(I)-F1(I)*FAC1
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
      DO I=1,NM1
         IM1=I+M1
         S1=0.0D0
         S2=0.0D0
         S3=0.0D0
         DO J=1,NM1
            JM1=J+M1
            BB=FMAS(I,J)
            S1=S1-BB*F1(JM1)
            S2=S2-BB*F2(JM1)
            S3=S3-BB*F3(JM1)
         END DO
         Z1(IM1)=Z1(IM1)+S1*FAC1
         Z2(IM1)=Z2(IM1)+S2*ALPHN-S3*BETAN
         Z3(IM1)=Z3(IM1)+S3*ALPHN+S2*BETAN
      END DO
      GOTO 48
C
C -----------------------------------------------------------
C
   6  CONTINUE
C ---  B IS A FULL MATRIX, JACOBIAN A BANDED MATRIX
C ---  THIS OPTION IS NOT PROVIDED
      RETURN
C
C -----------------------------------------------------------
C
   7  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX, HESSENBERG-OPTION
      DO I=1,N
         S2=-F2(I)
         S3=-F3(I)
         Z1(I)=Z1(I)-F1(I)*FAC1
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
      DO MM=N-2,1,-1
          MP=N-MM
          MP1=MP-1
          I=IPHES(MP)
          IF (I.EQ.MP) GOTO 746
          ZSAFE=Z1(MP)
          Z1(MP)=Z1(I)
          Z1(I)=ZSAFE
          ZSAFE=Z2(MP)
          Z2(MP)=Z2(I)
          Z2(I)=ZSAFE 
          ZSAFE=Z3(MP)
          Z3(MP)=Z3(I)
          Z3(I)=ZSAFE
 746      CONTINUE
          DO I=MP+1,N 
             E1IMP=FJAC(I,MP1)
             Z1(I)=Z1(I)-E1IMP*Z1(MP)
             Z2(I)=Z2(I)-E1IMP*Z2(MP)
             Z3(I)=Z3(I)-E1IMP*Z3(MP)
          END DO
       END DO
       CALL SOLH(N,LDE1,E1,1,Z1,IP1)
       CALL SOLHC(N,LDE1,E2R,E2I,1,Z2,Z3,IP2)
       DO MM=1,N-2
          MP=N-MM
          MP1=MP-1
          DO I=MP+1,N 
             E1IMP=FJAC(I,MP1)
             Z1(I)=Z1(I)+E1IMP*Z1(MP)
             Z2(I)=Z2(I)+E1IMP*Z2(MP)
             Z3(I)=Z3(I)+E1IMP*Z3(MP)
          END DO
          I=IPHES(MP)
          IF (I.EQ.MP) GOTO 750
          ZSAFE=Z1(MP)
          Z1(MP)=Z1(I)
          Z1(I)=ZSAFE
          ZSAFE=Z2(MP)
          Z2(MP)=Z2(I)
          Z2(I)=ZSAFE 
          ZSAFE=Z3(MP)
          Z3(MP)=Z3(I)
          Z3(I)=ZSAFE
 750      CONTINUE
      END DO
      RETURN
C
C -----------------------------------------------------------
C
  55  CONTINUE
      RETURN
      END
C
C     END OF SUBROUTINE SLVRAD
C
C ***********************************************************
C
      SUBROUTINE ESTRAD(N,FJAC,LDJAC,MLJAC,MUJAC,FMAS,LDMAS,MLMAS,MUMAS,
     &          H,DD1,DD2,DD3,FCN,NFCN,Y0,Y,IJOB,X,M1,M2,NM1,
     &          E1,LDE1,Z1,Z2,Z3,CONT,F1,F2,IP1,IPHES,SCAL,ERR,
     &          FIRST,REJECT,FAC1,RPAR,IPAR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION FJAC(LDJAC,N),FMAS(LDMAS,NM1),E1(LDE1,NM1),IP1(NM1),
     &     SCAL(N),IPHES(N),Z1(N),Z2(N),Z3(N),F1(N),F2(N),Y0(N),Y(N)
      DIMENSION CONT(N),RPAR(1),IPAR(1)
      LOGICAL FIRST,REJECT
      COMMON/LINAL/MLE,MUE,MBJAC,MBB,MDIAG,MDIFF,MBDIAG
      HEE1=DD1/H
      HEE2=DD2/H
      HEE3=DD3/H
      GOTO (1,2,3,4,5,6,7,55,55,55,11,12,13,14,15), IJOB
C
   1  CONTINUE
C ------  B=IDENTITY, JACOBIAN A FULL MATRIX
      DO  I=1,N 
         F2(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I)
         CONT(I)=F2(I)+Y0(I)
      END DO
      CALL SOL (N,LDE1,E1,CONT,IP1) 
      GOTO 77
C
  11  CONTINUE
C ------  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO I=1,N 
         F2(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I)
         CONT(I)=F2(I)+Y0(I)
      END DO
  48  MM=M1/M2
      DO J=1,M2
         SUM1=0.D0
         DO K=MM-1,0,-1
            SUM1=(CONT(J+K*M2)+SUM1)/FAC1
            DO I=1,NM1
               IM1=I+M1
               CONT(IM1)=CONT(IM1)+FJAC(I,J+K*M2)*SUM1
            END DO
         END DO
      END DO
      CALL SOL (NM1,LDE1,E1,CONT(M1+1),IP1) 
      DO I=M1,1,-1
         CONT(I)=(CONT(I)+CONT(M2+I))/FAC1
      END DO
      GOTO 77
C
   2  CONTINUE
C ------  B=IDENTITY, JACOBIAN A BANDED MATRIX
      DO I=1,N 
         F2(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I)
         CONT(I)=F2(I)+Y0(I)
      END DO
      CALL SOLB (N,LDE1,E1,MLE,MUE,CONT,IP1)
      GOTO 77
C
  12  CONTINUE
C ------  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER
      DO I=1,N 
         F2(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I)
         CONT(I)=F2(I)+Y0(I)
      END DO
  45  MM=M1/M2
      DO J=1,M2
         SUM1=0.D0
         DO K=MM-1,0,-1
            SUM1=(CONT(J+K*M2)+SUM1)/FAC1
            DO I=MAX(1,J-MUJAC),MIN(NM1,J+MLJAC)
               IM1=I+M1
               CONT(IM1)=CONT(IM1)+FJAC(I+MUJAC+1-J,J+K*M2)*SUM1
            END DO
         END DO
      END DO
      CALL SOLB (NM1,LDE1,E1,MLE,MUE,CONT(M1+1),IP1)
      DO I=M1,1,-1
         CONT(I)=(CONT(I)+CONT(M2+I))/FAC1
      END DO
      GOTO 77
C
   3  CONTINUE
C ------  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX
      DO I=1,N
         F1(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I)
      END DO
      DO I=1,N
         SUM=0.D0
         DO J=MAX(1,I-MLMAS),MIN(N,I+MUMAS)
            SUM=SUM+FMAS(I-J+MBDIAG,J)*F1(J)
         END DO
         F2(I)=SUM
         CONT(I)=SUM+Y0(I)
      END DO
      CALL SOL (N,LDE1,E1,CONT,IP1) 
      GOTO 77
C
  13  CONTINUE
C ------  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO I=1,M1
         F2(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I)
         CONT(I)=F2(I)+Y0(I)
      END DO
      DO I=M1+1,N
         F1(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I)
      END DO
      DO I=1,NM1
         SUM=0.D0
         DO J=MAX(1,I-MLMAS),MIN(NM1,I+MUMAS)
            SUM=SUM+FMAS(I-J+MBDIAG,J)*F1(J+M1)
         END DO
         IM1=I+M1
         F2(IM1)=SUM
         CONT(IM1)=SUM+Y0(IM1)
      END DO
      GOTO 48
C
   4  CONTINUE
C ------  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX
      DO I=1,N
         F1(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I)
      END DO
      DO I=1,N
         SUM=0.D0
         DO J=MAX(1,I-MLMAS),MIN(N,I+MUMAS)
            SUM=SUM+FMAS(I-J+MBDIAG,J)*F1(J)
         END DO
         F2(I)=SUM
         CONT(I)=SUM+Y0(I)
      END DO
      CALL SOLB (N,LDE1,E1,MLE,MUE,CONT,IP1)
      GOTO 77
C
  14  CONTINUE
C ------  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX, SECOND ORDER
      DO I=1,M1
         F2(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I)
         CONT(I)=F2(I)+Y0(I)
      END DO
      DO I=M1+1,N
         F1(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I)
      END DO
      DO I=1,NM1
         SUM=0.D0
         DO J=MAX(1,I-MLMAS),MIN(NM1,I+MUMAS)
            SUM=SUM+FMAS(I-J+MBDIAG,J)*F1(J+M1)
         END DO
         IM1=I+M1
         F2(IM1)=SUM
         CONT(IM1)=SUM+Y0(IM1)
      END DO
      GOTO 45
C
   5  CONTINUE
C ------  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX
      DO I=1,N
         F1(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I)
      END DO
      DO I=1,N
         SUM=0.D0
         DO J=1,N
            SUM=SUM+FMAS(I,J)*F1(J)
         END DO
         F2(I)=SUM
         CONT(I)=SUM+Y0(I)
      END DO
      CALL SOL (N,LDE1,E1,CONT,IP1) 
      GOTO 77
C
  15  CONTINUE
C ------  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO I=1,M1
         F2(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I)
         CONT(I)=F2(I)+Y0(I)
      END DO
      DO I=M1+1,N
         F1(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I)
      END DO
      DO I=1,NM1
         SUM=0.D0
         DO J=1,NM1
            SUM=SUM+FMAS(I,J)*F1(J+M1)
         END DO
         IM1=I+M1
         F2(IM1)=SUM
         CONT(IM1)=SUM+Y0(IM1)
      END DO
      GOTO 48
C
   6  CONTINUE
C ------  B IS A FULL MATRIX, JACOBIAN A BANDED MATRIX
C ------  THIS OPTION IS NOT PROVIDED
      RETURN
C
   7  CONTINUE
C ------  B=IDENTITY, JACOBIAN A FULL MATRIX, HESSENBERG-OPTION
      DO I=1,N 
         F2(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I)
         CONT(I)=F2(I)+Y0(I)
      END DO
      DO MM=N-2,1,-1
         MP=N-MM
         I=IPHES(MP)
         IF (I.EQ.MP) GOTO 310
         ZSAFE=CONT(MP)
         CONT(MP)=CONT(I)
         CONT(I)=ZSAFE
 310     CONTINUE
         DO I=MP+1,N 
            CONT(I)=CONT(I)-FJAC(I,MP-1)*CONT(MP)
         END DO
      END DO
      CALL SOLH(N,LDE1,E1,1,CONT,IP1)
      DO MM=1,N-2
         MP=N-MM
         DO I=MP+1,N 
            CONT(I)=CONT(I)+FJAC(I,MP-1)*CONT(MP)
         END DO
         I=IPHES(MP)
         IF (I.EQ.MP) GOTO 440
         ZSAFE=CONT(MP)
         CONT(MP)=CONT(I)
         CONT(I)=ZSAFE
 440     CONTINUE
      END DO
C
C --------------------------------------
C
  77  CONTINUE
      ERR=0.D0
      DO  I=1,N
         ERR=ERR+(CONT(I)/SCAL(I))**2
      END DO
      ERR=MAX(SQRT(ERR/N),1.D-10)
C
      IF (ERR.LT.1.D0) RETURN
      IF (FIRST.OR.REJECT) THEN
          DO I=1,N
             CONT(I)=Y(I)+CONT(I)
          END DO
          CALL FCN(N,X,CONT,F1,RPAR,IPAR)
          NFCN=NFCN+1
          DO I=1,N
             CONT(I)=F1(I)+F2(I)
          END DO
          GOTO (31,32,31,32,31,32,33,55,55,55,41,42,41,42,41), IJOB
C ------ FULL MATRIX OPTION
  31      CONTINUE
          CALL SOL(N,LDE1,E1,CONT,IP1) 
          GOTO 88
C ------ FULL MATRIX OPTION, SECOND ORDER
 41      CONTINUE
         DO J=1,M2
            SUM1=0.D0
            DO K=MM-1,0,-1
               SUM1=(CONT(J+K*M2)+SUM1)/FAC1
               DO I=1,NM1
                  IM1=I+M1
                  CONT(IM1)=CONT(IM1)+FJAC(I,J+K*M2)*SUM1
               END DO
            END DO
         END DO
         CALL SOL(NM1,LDE1,E1,CONT(M1+1),IP1) 
         DO I=M1,1,-1
            CONT(I)=(CONT(I)+CONT(M2+I))/FAC1
         END DO
         GOTO 88
C ------ BANDED MATRIX OPTION
 32      CONTINUE
         CALL SOLB (N,LDE1,E1,MLE,MUE,CONT,IP1)
         GOTO 88
C ------ BANDED MATRIX OPTION, SECOND ORDER
 42      CONTINUE
         DO J=1,M2
            SUM1=0.D0
            DO K=MM-1,0,-1
               SUM1=(CONT(J+K*M2)+SUM1)/FAC1
               DO I=MAX(1,J-MUJAC),MIN(NM1,J+MLJAC)
                  IM1=I+M1
                  CONT(IM1)=CONT(IM1)+FJAC(I+MUJAC+1-J,J+K*M2)*SUM1
               END DO
            END DO
         END DO
         CALL SOLB (NM1,LDE1,E1,MLE,MUE,CONT(M1+1),IP1)
         DO I=M1,1,-1
            CONT(I)=(CONT(I)+CONT(M2+I))/FAC1
         END DO
          GOTO 88
C ------ HESSENBERG MATRIX OPTION
  33      CONTINUE
          DO MM=N-2,1,-1
             MP=N-MM
             I=IPHES(MP)
             IF (I.EQ.MP) GOTO 510
             ZSAFE=CONT(MP)
             CONT(MP)=CONT(I)
             CONT(I)=ZSAFE
 510         CONTINUE
             DO I=MP+1,N 
                CONT(I)=CONT(I)-FJAC(I,MP-1)*CONT(MP)
             END DO
          END DO
          CALL SOLH(N,LDE1,E1,1,CONT,IP1)
          DO MM=1,N-2
             MP=N-MM
             DO I=MP+1,N 
                CONT(I)=CONT(I)+FJAC(I,MP-1)*CONT(MP)
             END DO
             I=IPHES(MP)
             IF (I.EQ.MP) GOTO 640
             ZSAFE=CONT(MP)
             CONT(MP)=CONT(I)
             CONT(I)=ZSAFE
 640         CONTINUE
          END DO
C -----------------------------------
   88     CONTINUE
          ERR=0.D0 
          DO I=1,N
             ERR=ERR+(CONT(I)/SCAL(I))**2
          END DO
          ERR=MAX(SQRT(ERR/N),1.D-10)
       END IF
       RETURN
C -----------------------------------------------------------
  55   CONTINUE
       RETURN
       END
C
C     END OF SUBROUTINE ESTRAD
C
C ***********************************************************
C
      SUBROUTINE ESTRAV(N,FJAC,LDJAC,MLJAC,MUJAC,FMAS,LDMAS,MLMAS,MUMAS,
     &          H,DD,FCN,NFCN,Y0,Y,IJOB,X,M1,M2,NM1,NS,NNS,
     &          E1,LDE1,ZZ,CONT,FF,IP1,IPHES,SCAL,ERR,
     &          FIRST,REJECT,FAC1,RPAR,IPAR)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION FJAC(LDJAC,N),FMAS(LDMAS,NM1),E1(LDE1,NM1),IP1(NM1),
     &     SCAL(N),IPHES(N),ZZ(NNS),FF(NNS),Y0(N),Y(N)
      DIMENSION DD(NS),CONT(N),RPAR(1),IPAR(1)
      LOGICAL FIRST,REJECT
      COMMON/LINAL/MLE,MUE,MBJAC,MBB,MDIAG,MDIFF,MBDIAG
      GOTO (1,2,3,4,5,6,7,55,55,55,11,12,13,14,15), IJOB
C
   1  CONTINUE
C ------  B=IDENTITY, JACOBIAN A FULL MATRIX
      DO  I=1,N
         SUM=0.D0
         DO K=1,NS
            SUM=SUM+DD(K)*ZZ(I+(K-1)*N)
         END DO
         FF(I+N)=SUM/H
         CONT(I)=FF(I+N)+Y0(I)
      END DO
      CALL SOL (N,LDE1,E1,CONT,IP1) 
      GOTO 77
C
  11  CONTINUE
C ------  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO  I=1,N
         SUM=0.D0
         DO K=1,NS
            SUM=SUM+DD(K)*ZZ(I+(K-1)*N)
         END DO
         FF(I+N)=SUM/H
         CONT(I)=FF(I+N)+Y0(I)
      END DO
  48  MM=M1/M2
      DO J=1,M2
         SUM1=0.D0
         DO K=MM-1,0,-1
            SUM1=(CONT(J+K*M2)+SUM1)/FAC1
            DO I=1,NM1
               IM1=I+M1
               CONT(IM1)=CONT(IM1)+FJAC(I,J+K*M2)*SUM1
            END DO
         END DO
      END DO
      CALL SOL (NM1,LDE1,E1,CONT(M1+1),IP1) 
      DO I=M1,1,-1
         CONT(I)=(CONT(I)+CONT(M2+I))/FAC1
      END DO
      GOTO 77
C
   2  CONTINUE
C ------  B=IDENTITY, JACOBIAN A BANDED MATRIX
      DO  I=1,N
         SUM=0.D0
         DO K=1,NS
            SUM=SUM+DD(K)*ZZ(I+(K-1)*N)
         END DO
         FF(I+N)=SUM/H
         CONT(I)=FF(I+N)+Y0(I)
      END DO
      CALL SOLB (N,LDE1,E1,MLE,MUE,CONT,IP1)
      GOTO 77
C
  12  CONTINUE
C ------  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER
      DO  I=1,N
         SUM=0.D0
         DO K=1,NS
            SUM=SUM+DD(K)*ZZ(I+(K-1)*N)
         END DO
         FF(I+N)=SUM/H
         CONT(I)=FF(I+N)+Y0(I)
      END DO
  45  MM=M1/M2
      DO J=1,M2
         SUM1=0.D0
         DO K=MM-1,0,-1
            SUM1=(CONT(J+K*M2)+SUM1)/FAC1
            DO I=MAX(1,J-MUJAC),MIN(NM1,J+MLJAC)
               IM1=I+M1
               CONT(IM1)=CONT(IM1)+FJAC(I+MUJAC+1-J,J+K*M2)*SUM1
            END DO
         END DO
      END DO
      CALL SOLB (NM1,LDE1,E1,MLE,MUE,CONT(M1+1),IP1)
      DO I=M1,1,-1
         CONT(I)=(CONT(I)+CONT(M2+I))/FAC1
      END DO
      GOTO 77
C
   3  CONTINUE
C ------  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX
      DO  I=1,N
         SUM=0.D0
         DO K=1,NS
            SUM=SUM+DD(K)*ZZ(I+(K-1)*N)
         END DO
         FF(I)=SUM/H
      END DO
      DO I=1,N
         SUM=0.D0
         DO J=MAX(1,I-MLMAS),MIN(N,I+MUMAS)
            SUM=SUM+FMAS(I-J+MBDIAG,J)*FF(J)
         END DO
         FF(I+N)=SUM
         CONT(I)=SUM+Y0(I)
      END DO
      CALL SOL (N,LDE1,E1,CONT,IP1) 
      GOTO 77
C
  13  CONTINUE
C ------  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO  I=1,M1
         SUM=0.D0
         DO K=1,NS
            SUM=SUM+DD(K)*ZZ(I+(K-1)*N)
         END DO
         FF(I+N)=SUM/H
         CONT(I)=FF(I+N)+Y0(I)
      END DO
      DO I=M1+1,N
         SUM=0.D0
         DO K=1,NS
            SUM=SUM+DD(K)*ZZ(I+(K-1)*N)
         END DO
         FF(I)=SUM/H
      END DO
      DO I=1,NM1
         SUM=0.D0
         DO J=MAX(1,I-MLMAS),MIN(NM1,I+MUMAS)
            SUM=SUM+FMAS(I-J+MBDIAG,J)*FF(J+M1)
         END DO
         IM1=I+M1
         FF(IM1+N)=SUM
         CONT(IM1)=SUM+Y0(IM1)
      END DO
      GOTO 48
C
   4  CONTINUE
C ------  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX
      DO  I=1,N
         SUM=0.D0
         DO K=1,NS
            SUM=SUM+DD(K)*ZZ(I+(K-1)*N)
         END DO
         FF(I)=SUM/H
      END DO
      DO I=1,N
         SUM=0.D0
         DO J=MAX(1,I-MLMAS),MIN(N,I+MUMAS)
            SUM=SUM+FMAS(I-J+MBDIAG,J)*FF(J)
         END DO
         FF(I+N)=SUM
         CONT(I)=SUM+Y0(I)
      END DO
      CALL SOLB (N,LDE1,E1,MLE,MUE,CONT,IP1)
      GOTO 77
C
  14  CONTINUE
C ------  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX, SECOND ORDER
      DO  I=1,M1
         SUM=0.D0
         DO K=1,NS
            SUM=SUM+DD(K)*ZZ(I+(K-1)*N)
         END DO
         FF(I+N)=SUM/H
         CONT(I)=FF(I+N)+Y0(I)
      END DO
      DO I=M1+1,N
         SUM=0.D0
         DO K=1,NS
            SUM=SUM+DD(K)*ZZ(I+(K-1)*N)
         END DO
         FF(I)=SUM/H
      END DO
      DO I=1,NM1
         SUM=0.D0
         DO J=MAX(1,I-MLMAS),MIN(NM1,I+MUMAS)
            SUM=SUM+FMAS(I-J+MBDIAG,J)*FF(J+M1)
         END DO
         IM1=I+M1
         FF(IM1+N)=SUM
         CONT(IM1)=SUM+Y0(IM1)
      END DO
      GOTO 45
C
   5  CONTINUE
C ------  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX
      DO I=1,N
         SUM=0.D0
         DO K=1,NS
            SUM=SUM+DD(K)*ZZ(I+(K-1)*N)
         END DO
         FF(I)=SUM/H
      END DO
      DO I=1,N
         SUM=0.D0
         DO J=1,N
            SUM=SUM+FMAS(I,J)*FF(J)
         END DO
         FF(I+N)=SUM
         CONT(I)=SUM+Y0(I)
      END DO
      CALL SOL (N,LDE1,E1,CONT,IP1) 
      GOTO 77
C
  15  CONTINUE
C ------  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO  I=1,M1
         SUM=0.D0
         DO K=1,NS
            SUM=SUM+DD(K)*ZZ(I+(K-1)*N)
         END DO
         FF(I+N)=SUM/H
         CONT(I)=FF(I+N)+Y0(I)
      END DO
      DO I=M1+1,N
         SUM=0.D0
         DO K=1,NS
            SUM=SUM+DD(K)*ZZ(I+(K-1)*N)
         END DO
         FF(I)=SUM/H
      END DO
      DO I=1,NM1
         SUM=0.D0
         DO J=1,NM1
            SUM=SUM+FMAS(I,J)*FF(J+M1)
         END DO
         IM1=I+M1
         FF(IM1+N)=SUM
         CONT(IM1)=SUM+Y0(IM1)
      END DO
      GOTO 48
C
   6  CONTINUE
C ------  B IS A FULL MATRIX, JACOBIAN A BANDED MATRIX
C ------  THIS OPTION IS NOT PROVIDED
      RETURN
C
   7  CONTINUE
C ------  B=IDENTITY, JACOBIAN A FULL MATRIX, HESSENBERG-OPTION
      DO  I=1,N
         SUM=0.D0
         DO K=1,NS
            SUM=SUM+DD(K)*ZZ(I+(K-1)*N)
         END DO
         FF(I+N)=SUM/H
         CONT(I)=FF(I+N)+Y0(I)
      END DO
      DO MM=N-2,1,-1
         MP=N-MM
         I=IPHES(MP)
         IF (I.EQ.MP) GOTO 310
         ZSAFE=CONT(MP)
         CONT(MP)=CONT(I)
         CONT(I)=ZSAFE
 310     CONTINUE
         DO I=MP+1,N 
            CONT(I)=CONT(I)-FJAC(I,MP-1)*CONT(MP)
         END DO
      END DO
      CALL SOLH(N,LDE1,E1,1,CONT,IP1)
      DO MM=1,N-2
         MP=N-MM
         DO I=MP+1,N 
            CONT(I)=CONT(I)+FJAC(I,MP-1)*CONT(MP)
         END DO
         I=IPHES(MP)
         IF (I.EQ.MP) GOTO 440
         ZSAFE=CONT(MP)
         CONT(MP)=CONT(I)
         CONT(I)=ZSAFE
 440     CONTINUE
      END DO
C
C --------------------------------------
C
  77  CONTINUE
      ERR=0.D0
      DO  I=1,N
         ERR=ERR+(CONT(I)/SCAL(I))**2
      END DO
      ERR=MAX(SQRT(ERR/N),1.D-10)
C
      IF (ERR.LT.1.D0) RETURN
      IF (FIRST.OR.REJECT) THEN
          DO I=1,N
             CONT(I)=Y(I)+CONT(I)
          END DO
          CALL FCN(N,X,CONT,FF,RPAR,IPAR)
          NFCN=NFCN+1
          DO I=1,N
             CONT(I)=FF(I)+FF(I+N)
          END DO
          GOTO (31,32,31,32,31,32,33,55,55,55,41,42,41,42,41), IJOB
C ------ FULL MATRIX OPTION
 31      CONTINUE
         CALL SOL (N,LDE1,E1,CONT,IP1) 
          GOTO 88
C ------ FULL MATRIX OPTION, SECOND ORDER
 41      CONTINUE
         DO J=1,M2
            SUM1=0.D0
            DO K=MM-1,0,-1
               SUM1=(CONT(J+K*M2)+SUM1)/FAC1
               DO I=1,NM1
                  IM1=I+M1
                  CONT(IM1)=CONT(IM1)+FJAC(I,J+K*M2)*SUM1
               END DO
            END DO
         END DO
         CALL SOL (NM1,LDE1,E1,CONT(M1+1),IP1) 
         DO I=M1,1,-1
            CONT(I)=(CONT(I)+CONT(M2+I))/FAC1
         END DO
          GOTO 88
C ------ BANDED MATRIX OPTION
 32      CONTINUE
         CALL SOLB (N,LDE1,E1,MLE,MUE,CONT,IP1)
          GOTO 88
C ------ BANDED MATRIX OPTION, SECOND ORDER
 42      CONTINUE
         DO J=1,M2
            SUM1=0.D0
            DO K=MM-1,0,-1
               SUM1=(CONT(J+K*M2)+SUM1)/FAC1
               DO I=MAX(1,J-MUJAC),MIN(NM1,J+MLJAC)
                  IM1=I+M1
                  CONT(IM1)=CONT(IM1)+FJAC(I+MUJAC+1-J,J+K*M2)*SUM1
               END DO
            END DO
         END DO
         CALL SOLB (NM1,LDE1,E1,MLE,MUE,CONT(M1+1),IP1)
         DO I=M1,1,-1
            CONT(I)=(CONT(I)+CONT(M2+I))/FAC1
         END DO
          GOTO 88
C ------ HESSENBERG MATRIX OPTION
  33      CONTINUE
          DO MM=N-2,1,-1
             MP=N-MM
             I=IPHES(MP)
             IF (I.EQ.MP) GOTO 510
             ZSAFE=CONT(MP)
             CONT(MP)=CONT(I)
             CONT(I)=ZSAFE
 510         CONTINUE
             DO I=MP+1,N 
                CONT(I)=CONT(I)-FJAC(I,MP-1)*CONT(MP)
             END DO
          END DO
          CALL SOLH(N,LDE1,E1,1,CONT,IP1)
          DO MM=1,N-2
             MP=N-MM
             DO I=MP+1,N 
                CONT(I)=CONT(I)+FJAC(I,MP-1)*CONT(MP)
             END DO
             I=IPHES(MP)
             IF (I.EQ.MP) GOTO 640
             ZSAFE=CONT(MP)
             CONT(MP)=CONT(I)
             CONT(I)=ZSAFE
 640         CONTINUE
          END DO
C -----------------------------------
  88      CONTINUE
          ERR=0.D0 
          DO I=1,N
             ERR=ERR+(CONT(I)/SCAL(I))**2
          END DO
          ERR=MAX(SQRT(ERR/N),1.D-10)
       END IF
       RETURN
C
C -----------------------------------------------------------
C
  55  CONTINUE
      RETURN
       END
C
C     END OF SUBROUTINE ESTRAV
C
C ***********************************************************
C
      SUBROUTINE SLVROD(N,FJAC,LDJAC,MLJAC,MUJAC,FMAS,LDMAS,MLMAS,MUMAS,
     &          M1,M2,NM1,FAC1,E,LDE,IP,DY,AK,FX,YNEW,HD,IJOB,STAGE1)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION FJAC(LDJAC,N),FMAS(LDMAS,NM1),E(LDE,NM1),
     &          IP(NM1),DY(N),AK(N),FX(N),YNEW(N)
      LOGICAL STAGE1
      COMMON/LINAL/MLE,MUE,MBJAC,MBB,MDIAG,MDIFF,MBDIAG
C
      IF (HD.EQ.0.D0) THEN
         DO  I=1,N
           AK(I)=DY(I)
         END DO
      ELSE
         DO I=1,N
            AK(I)=DY(I)+HD*FX(I)
         END DO
      END IF
C
      GOTO (1,2,3,4,5,6,55,55,55,55,11,12,13,13,15), IJOB
C
C -----------------------------------------------------------
C
   1  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX
      IF (STAGE1) THEN
         DO I=1,N
            AK(I)=AK(I)+YNEW(I)
         END DO
      END IF
      CALL SOL (N,LDE,E,AK,IP)
      RETURN
C
C -----------------------------------------------------------
C
  11  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER
      IF (STAGE1) THEN
         DO I=1,N
            AK(I)=AK(I)+YNEW(I)
         END DO
      END IF
 48   MM=M1/M2
      DO J=1,M2
         SUM=0.D0
         DO K=MM-1,0,-1
            JKM=J+K*M2
            SUM=(AK(JKM)+SUM)/FAC1
            DO I=1,NM1
               IM1=I+M1
               AK(IM1)=AK(IM1)+FJAC(I,JKM)*SUM
            END DO
         END DO
      END DO
      CALL SOL (NM1,LDE,E,AK(M1+1),IP)
      DO I=M1,1,-1
         AK(I)=(AK(I)+AK(M2+I))/FAC1
      END DO
      RETURN
C
C -----------------------------------------------------------
C
   2  CONTINUE
C ---  B=IDENTITY, JACOBIAN A BANDED MATRIX
      IF (STAGE1) THEN
         DO I=1,N
            AK(I)=AK(I)+YNEW(I)
         END DO
      END IF
      CALL SOLB (N,LDE,E,MLE,MUE,AK,IP)
      RETURN
C
C -----------------------------------------------------------
C
  12  CONTINUE
C ---  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER
      IF (STAGE1) THEN
         DO I=1,N
            AK(I)=AK(I)+YNEW(I)
         END DO
      END IF
  45  MM=M1/M2
      DO J=1,M2
         SUM=0.D0
         DO K=MM-1,0,-1
            JKM=J+K*M2
            SUM=(AK(JKM)+SUM)/FAC1
            DO I=MAX(1,J-MUJAC),MIN(NM1,J+MLJAC)
               IM1=I+M1
               AK(IM1)=AK(IM1)+FJAC(I+MUJAC+1-J,JKM)*SUM
            END DO
         END DO
      END DO
      CALL SOLB (NM1,LDE,E,MLE,MUE,AK(M1+1),IP)
      DO I=M1,1,-1
         AK(I)=(AK(I)+AK(M2+I))/FAC1
      END DO
      RETURN
C
C -----------------------------------------------------------
C
   3  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX
      IF (STAGE1) THEN
      DO  I=1,N
         SUM=0.D0
         DO  J=MAX(1,I-MLMAS),MIN(N,I+MUMAS)
            SUM=SUM+FMAS(I-J+MBDIAG,J)*YNEW(J)
         END DO
         AK(I)=AK(I)+SUM
      END DO
      END IF
      CALL SOL (N,LDE,E,AK,IP)
      RETURN
C
C -----------------------------------------------------------
C
  13  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER
      IF (STAGE1) THEN
         DO I=1,M1
            AK(I)=AK(I)+YNEW(I)
         END DO
         DO I=1,NM1
            SUM=0.D0
            DO J=MAX(1,I-MLMAS),MIN(NM1,I+MUMAS)
                SUM=SUM+FMAS(I-J+MBDIAG,J)*YNEW(J+M1)
            END DO
            IM1=I+M1
            AK(IM1)=AK(IM1)+SUM
         END DO
      END IF
      IF (IJOB.EQ.14) GOTO 45
      GOTO 48
C
C -----------------------------------------------------------
C
   4  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX
      IF (STAGE1) THEN
      DO I=1,N
         SUM=0.D0
         DO J=MAX(1,I-MLMAS),MIN(N,I+MUMAS)
            SUM=SUM+FMAS(I-J+MBDIAG,J)*YNEW(J)
         END DO
         AK(I)=AK(I)+SUM
      END DO
      END IF
      CALL SOLB (N,LDE,E,MLE,MUE,AK,IP)
      RETURN
C
C -----------------------------------------------------------
C
   5  CONTINUE
C ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX
      IF (STAGE1) THEN
      DO I=1,N
         SUM=0.D0
         DO J=1,N
            SUM=SUM+FMAS(I,J)*YNEW(J)
         END DO
         AK(I)=AK(I)+SUM
      END DO
      END IF
      CALL SOL (N,LDE,E,AK,IP)
      RETURN
C
C -----------------------------------------------------------
C
  15  CONTINUE
C ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER
      IF (STAGE1) THEN
         DO I=1,M1
            AK(I)=AK(I)+YNEW(I)
         END DO
         DO I=1,NM1
            SUM=0.D0
            DO J=1,NM1
               SUM=SUM+FMAS(I,J)*YNEW(J+M1)
            END DO
            IM1=I+M1
            AK(IM1)=AK(IM1)+SUM
         END DO
      END IF
      GOTO 48
C
C -----------------------------------------------------------
C
   6  CONTINUE
C ---  B IS A FULL MATRIX, JACOBIAN A BANDED MATRIX
C ---  THIS OPTION IS NOT PROVIDED
      IF (STAGE1) THEN
      DO 624 I=1,N
         SUM=0.D0
         DO 623 J=1,N
  623       SUM=SUM+FMAS(I,J)*YNEW(J)
  624    AK(I)=AK(I)+SUM
      CALL SOLB (N,LDE,E,MLE,MUE,AK,IP)
      END IF
      RETURN
C
C -----------------------------------------------------------
C
  55  CONTINUE
      RETURN
      END
C
C     END OF SUBROUTINE SLVROD
C
C
C ***********************************************************
C
      SUBROUTINE SLVSEU(N,FJAC,LDJAC,MLJAC,MUJAC,FMAS,LDMAS,MLMAS,MUMAS,
     &          M1,M2,NM1,FAC1,E,LDE,IP,IPHES,DEL,IJOB)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION FJAC(LDJAC,N),FMAS(LDMAS,NM1),E(LDE,NM1),DEL(N)
      DIMENSION IP(NM1),IPHES(N)
      COMMON/LINAL/MLE,MUE,MBJAC,MBB,MDIAG,MDIFF,MBDIAG
C
      GOTO (1,2,1,2,1,55,7,55,55,55,11,12,11,12,11), IJOB
C
C -----------------------------------------------------------
C
   1  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX
      CALL SOL (N,LDE,E,DEL,IP)
      RETURN
C
C -----------------------------------------------------------
C
  11  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER
      MM=M1/M2
      DO J=1,M2
         SUM=0.D0
         DO K=MM-1,0,-1
            JKM=J+K*M2
            SUM=(DEL(JKM)+SUM)/FAC1
            DO I=1,NM1
               IM1=I+M1
               DEL(IM1)=DEL(IM1)+FJAC(I,JKM)*SUM
            END DO
         END DO
      END DO
      CALL SOL (NM1,LDE,E,DEL(M1+1),IP)
      DO I=M1,1,-1
         DEL(I)=(DEL(I)+DEL(M2+I))/FAC1
      END DO
      RETURN
C
C -----------------------------------------------------------
C
   2  CONTINUE
C ---  B=IDENTITY, JACOBIAN A BANDED MATRIX
      CALL SOLB (N,LDE,E,MLE,MUE,DEL,IP)
      RETURN
C
C -----------------------------------------------------------
C
  12  CONTINUE
C ---  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER
      MM=M1/M2
      DO J=1,M2
         SUM=0.D0
         DO K=MM-1,0,-1
            JKM=J+K*M2
            SUM=(DEL(JKM)+SUM)/FAC1
            DO I=MAX(1,J-MUJAC),MIN(NM1,J+MLJAC)
               IM1=I+M1
               DEL(IM1)=DEL(IM1)+FJAC(I+MUJAC+1-J,JKM)*SUM
            END DO
         END DO
      END DO
      CALL SOLB (NM1,LDE,E,MLE,MUE,DEL(M1+1),IP)
      DO I=M1,1,-1
         DEL(I)=(DEL(I)+DEL(M2+I))/FAC1
      END DO
      RETURN
C
C -----------------------------------------------------------
C
   7  CONTINUE
C ---  HESSENBERG OPTION
      DO MMM=N-2,1,-1
         MP=N-MMM
         MP1=MP-1
         I=IPHES(MP)
         IF (I.EQ.MP) GOTO 110
         ZSAFE=DEL(MP)
         DEL(MP)=DEL(I)
         DEL(I)=ZSAFE
 110     CONTINUE
         DO I=MP+1,N 
            DEL(I)=DEL(I)-FJAC(I,MP1)*DEL(MP)
         END DO
      END DO
      CALL SOLH(N,LDE,E,1,DEL,IP)
      DO MMM=1,N-2
         MP=N-MMM
         MP1=MP-1
         DO I=MP+1,N 
            DEL(I)=DEL(I)+FJAC(I,MP1)*DEL(MP)
         END DO
         I=IPHES(MP)
         IF (I.EQ.MP) GOTO 240
         ZSAFE=DEL(MP)
         DEL(MP)=DEL(I)
         DEL(I)=ZSAFE
 240     CONTINUE
      END DO
      RETURN
C
C -----------------------------------------------------------
C
  55  CONTINUE
      RETURN
      END
C
C     END OF SUBROUTINE SLVSEU
C

      SUBROUTINE DEC (N, NDIM, A, IP, IER)
C VERSION REAL DOUBLE PRECISION
      INTEGER N,NDIM,IP,IER,NM1,K,KP1,M,I,J
      DOUBLE PRECISION A,T
      DIMENSION A(NDIM,N), IP(N)
C-----------------------------------------------------------------------
C  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION.
C  INPUT..
C     N = ORDER OF MATRIX.
C     NDIM = DECLARED DIMENSION OF ARRAY  A .
C     A = MATRIX TO BE TRIANGULARIZED.
C  OUTPUT..
C     A(I,J), I.LE.J = UPPER TRIANGULAR FACTOR, U .
C     A(I,J), I.GT.J = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I - L.
C     IP(K), K.LT.N = INDEX OF K-TH PIVOT ROW.
C     IP(N) = (-1)**(NUMBER OF INTERCHANGES) OR O .
C     IER = 0 IF MATRIX A IS NONSINGULAR, OR K IF FOUND TO BE
C           SINGULAR AT STAGE K.
C  USE  SOL  TO OBTAIN SOLUTION OF LINEAR SYSTEM.
C  DETERM(A) = IP(N)*A(1,1)*A(2,2)*...*A(N,N).
C  IF IP(N)=O, A IS SINGULAR, SOL WILL DIVIDE BY ZERO.
C
C  REFERENCE..
C     C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER,
C     C.A.C.M. 15 (1972), P. 274.
C-----------------------------------------------------------------------
      IER = 0
      IP(N) = 1
      IF (N .EQ. 1) GO TO 70
      NM1 = N - 1
      DO 60 K = 1,NM1
        KP1 = K + 1
        M = K
        DO 10 I = KP1,N
          IF (DABS(A(I,K)) .GT. DABS(A(M,K))) M = I  
 10     CONTINUE
        IP(K) = M
        T = A(M,K)
        IF (M .EQ. K) GO TO 20
        IP(N) = -IP(N)
        A(M,K) = A(K,K)
        A(K,K) = T
 20     CONTINUE
        IF (T .EQ. 0.D0) GO TO 80
        T = 1.D0/T
        DO 30 I = KP1,N
 30       A(I,K) = -A(I,K)*T
        DO 50 J = KP1,N
          T = A(M,J)
          A(M,J) = A(K,J)
          A(K,J) = T
          IF (T .EQ. 0.D0) GO TO 45
          DO 40 I = KP1,N
 40         A(I,J) = A(I,J) + A(I,K)*T
 45       CONTINUE
 50       CONTINUE
 60     CONTINUE
 70   K = N
      IF (A(N,N) .EQ. 0.D0) GO TO 80
      RETURN
 80   IER = K
      IP(N) = 0
      RETURN
C----------------------- END OF SUBROUTINE DEC -------------------------
      END
C
C
      SUBROUTINE SOL (N, NDIM, A, B, IP)
C VERSION REAL DOUBLE PRECISION
      INTEGER N,NDIM,IP,NM1,K,KP1,M,I,KB,KM1
      DOUBLE PRECISION A,B,T
      DIMENSION A(NDIM,N), B(N), IP(N)
C-----------------------------------------------------------------------
C  SOLUTION OF LINEAR SYSTEM, A*X = B .
C  INPUT..
C    N = ORDER OF MATRIX.
C    NDIM = DECLARED DIMENSION OF ARRAY  A .
C    A = TRIANGULARIZED MATRIX OBTAINED FROM DEC.
C    B = RIGHT HAND SIDE VECTOR.
C    IP = PIVOT VECTOR OBTAINED FROM DEC.
C  DO NOT USE IF DEC HAS SET IER .NE. 0.
C  OUTPUT..
C    B = SOLUTION VECTOR, X .
C-----------------------------------------------------------------------
      IF (N .EQ. 1) GO TO 50
      NM1 = N - 1
      DO 20 K = 1,NM1
        KP1 = K + 1
        M = IP(K)
        T = B(M)
        B(M) = B(K)
        B(K) = T
        DO 10 I = KP1,N
 10       B(I) = B(I) + A(I,K)*T
 20     CONTINUE
      DO 40 KB = 1,NM1
        KM1 = N - KB
        K = KM1 + 1
        B(K) = B(K)/A(K,K)
        T = -B(K)
        DO 30 I = 1,KM1
 30       B(I) = B(I) + A(I,K)*T
 40     CONTINUE
 50   B(1) = B(1)/A(1,1)
      RETURN
C----------------------- END OF SUBROUTINE SOL -------------------------
      END
c
c
      SUBROUTINE DECH (N, NDIM, A, LB, IP, IER)
C VERSION REAL DOUBLE PRECISION
      INTEGER N,NDIM,IP,IER,NM1,K,KP1,M,I,J,LB,NA
      DOUBLE PRECISION A,T
      DIMENSION A(NDIM,N), IP(N)
C-----------------------------------------------------------------------
C  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION OF A HESSENBERG
C  MATRIX WITH LOWER BANDWIDTH LB
C  INPUT..
C     N = ORDER OF MATRIX A.
C     NDIM = DECLARED DIMENSION OF ARRAY  A .
C     A = MATRIX TO BE TRIANGULARIZED.
C     LB = LOWER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED, LB.GE.1).
C  OUTPUT..
C     A(I,J), I.LE.J = UPPER TRIANGULAR FACTOR, U .
C     A(I,J), I.GT.J = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I - L.
C     IP(K), K.LT.N = INDEX OF K-TH PIVOT ROW.
C     IP(N) = (-1)**(NUMBER OF INTERCHANGES) OR O .
C     IER = 0 IF MATRIX A IS NONSINGULAR, OR K IF FOUND TO BE
C           SINGULAR AT STAGE K.
C  USE  SOLH  TO OBTAIN SOLUTION OF LINEAR SYSTEM.
C  DETERM(A) = IP(N)*A(1,1)*A(2,2)*...*A(N,N).
C  IF IP(N)=O, A IS SINGULAR, SOL WILL DIVIDE BY ZERO.
C
C  REFERENCE..
C     THIS IS A SLIGHT MODIFICATION OF
C     C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER,
C     C.A.C.M. 15 (1972), P. 274.
C-----------------------------------------------------------------------
      IER = 0
      IP(N) = 1
      IF (N .EQ. 1) GO TO 70
      NM1 = N - 1
      DO 60 K = 1,NM1
        KP1 = K + 1
        M = K
        NA = MIN0(N,LB+K)
        DO 10 I = KP1,NA
          IF (DABS(A(I,K)) .GT. DABS(A(M,K))) M = I
 10     CONTINUE
        IP(K) = M
        T = A(M,K)
        IF (M .EQ. K) GO TO 20
        IP(N) = -IP(N)
        A(M,K) = A(K,K)
        A(K,K) = T
 20     CONTINUE
        IF (T .EQ. 0.D0) GO TO 80
        T = 1.D0/T
        DO 30 I = KP1,NA
 30       A(I,K) = -A(I,K)*T
        DO 50 J = KP1,N
          T = A(M,J)
          A(M,J) = A(K,J)
          A(K,J) = T
          IF (T .EQ. 0.D0) GO TO 45
          DO 40 I = KP1,NA
 40         A(I,J) = A(I,J) + A(I,K)*T
 45       CONTINUE
 50       CONTINUE
 60     CONTINUE
 70   K = N
      IF (A(N,N) .EQ. 0.D0) GO TO 80
      RETURN
 80   IER = K
      IP(N) = 0
      RETURN
C----------------------- END OF SUBROUTINE DECH ------------------------
      END
C
C
      SUBROUTINE SOLH (N, NDIM, A, LB, B, IP)
C VERSION REAL DOUBLE PRECISION
      INTEGER N,NDIM,IP,NM1,K,KP1,M,I,KB,KM1,LB,NA
      DOUBLE PRECISION A,B,T
      DIMENSION A(NDIM,N), B(N), IP(N)
C-----------------------------------------------------------------------
C  SOLUTION OF LINEAR SYSTEM, A*X = B .
C  INPUT..
C    N = ORDER OF MATRIX A.
C    NDIM = DECLARED DIMENSION OF ARRAY  A .
C    A = TRIANGULARIZED MATRIX OBTAINED FROM DECH.
C    LB = LOWER BANDWIDTH OF A.
C    B = RIGHT HAND SIDE VECTOR.
C    IP = PIVOT VECTOR OBTAINED FROM DEC.
C  DO NOT USE IF DECH HAS SET IER .NE. 0.
C  OUTPUT..
C    B = SOLUTION VECTOR, X .
C-----------------------------------------------------------------------
      IF (N .EQ. 1) GO TO 50
      NM1 = N - 1
      DO 20 K = 1,NM1
        KP1 = K + 1
        M = IP(K)
        T = B(M)
        B(M) = B(K)
        B(K) = T
        NA = MIN0(N,LB+K)
        DO 10 I = KP1,NA
 10       B(I) = B(I) + A(I,K)*T
 20     CONTINUE
      DO 40 KB = 1,NM1
        KM1 = N - KB
        K = KM1 + 1
        B(K) = B(K)/A(K,K)
        T = -B(K)
        DO 30 I = 1,KM1
 30       B(I) = B(I) + A(I,K)*T
 40     CONTINUE
 50   B(1) = B(1)/A(1,1)
      RETURN
C----------------------- END OF SUBROUTINE SOLH ------------------------
      END
C
      SUBROUTINE DECC (N, NDIM, AR, AI, IP, IER)
C VERSION COMPLEX DOUBLE PRECISION
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER N,NDIM,IP,IER,NM1,K,KP1,M,I,J
      DIMENSION AR(NDIM,N), AI(NDIM,N), IP(N)
C-----------------------------------------------------------------------
C  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION
C  ------ MODIFICATION FOR COMPLEX MATRICES --------
C  INPUT..
C     N = ORDER OF MATRIX.
C     NDIM = DECLARED DIMENSION OF ARRAYS  AR AND AI .
C     (AR, AI) = MATRIX TO BE TRIANGULARIZED.
C  OUTPUT..
C     AR(I,J), I.LE.J = UPPER TRIANGULAR FACTOR, U ; REAL PART.
C     AI(I,J), I.LE.J = UPPER TRIANGULAR FACTOR, U ; IMAGINARY PART.
C     AR(I,J), I.GT.J = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I - L.
C                                                    REAL PART.
C     AI(I,J), I.GT.J = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I - L.
C                                                    IMAGINARY PART.
C     IP(K), K.LT.N = INDEX OF K-TH PIVOT ROW.
C     IP(N) = (-1)**(NUMBER OF INTERCHANGES) OR O .
C     IER = 0 IF MATRIX A IS NONSINGULAR, OR K IF FOUND TO BE
C           SINGULAR AT STAGE K.
C  USE  SOL  TO OBTAIN SOLUTION OF LINEAR SYSTEM.
C  IF IP(N)=O, A IS SINGULAR, SOL WILL DIVIDE BY ZERO.
C
C  REFERENCE..
C     C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER,
C     C.A.C.M. 15 (1972), P. 274.
C-----------------------------------------------------------------------
      IER = 0
      IP(N) = 1
      IF (N .EQ. 1) GO TO 70
      NM1 = N - 1
      DO 60 K = 1,NM1
        KP1 = K + 1
        M = K
        DO 10 I = KP1,N
          IF (DABS(AR(I,K))+DABS(AI(I,K)) .GT.
     &          DABS(AR(M,K))+DABS(AI(M,K))) M = I
 10     CONTINUE
        IP(K) = M
        TR = AR(M,K)
        TI = AI(M,K)
        IF (M .EQ. K) GO TO 20
        IP(N) = -IP(N)
        AR(M,K) = AR(K,K)
        AI(M,K) = AI(K,K)
        AR(K,K) = TR
        AI(K,K) = TI
 20     CONTINUE
        IF (DABS(TR)+DABS(TI) .EQ. 0.D0) GO TO 80
        DEN=TR*TR+TI*TI
        TR=TR/DEN
        TI=-TI/DEN
        DO 30 I = KP1,N
          PRODR=AR(I,K)*TR-AI(I,K)*TI
          PRODI=AI(I,K)*TR+AR(I,K)*TI
          AR(I,K)=-PRODR
          AI(I,K)=-PRODI
 30       CONTINUE
        DO 50 J = KP1,N
          TR = AR(M,J)
          TI = AI(M,J)
          AR(M,J) = AR(K,J)
          AI(M,J) = AI(K,J)
          AR(K,J) = TR
          AI(K,J) = TI
          IF (DABS(TR)+DABS(TI) .EQ. 0.D0) GO TO 48
          IF (TI .EQ. 0.D0) THEN
            DO 40 I = KP1,N
            PRODR=AR(I,K)*TR
            PRODI=AI(I,K)*TR
            AR(I,J) = AR(I,J) + PRODR
            AI(I,J) = AI(I,J) + PRODI
 40         CONTINUE
            GO TO 48
          END IF
          IF (TR .EQ. 0.D0) THEN
            DO 45 I = KP1,N
            PRODR=-AI(I,K)*TI
            PRODI=AR(I,K)*TI
            AR(I,J) = AR(I,J) + PRODR
            AI(I,J) = AI(I,J) + PRODI
 45         CONTINUE
            GO TO 48
          END IF
          DO 47 I = KP1,N
            PRODR=AR(I,K)*TR-AI(I,K)*TI
            PRODI=AI(I,K)*TR+AR(I,K)*TI
            AR(I,J) = AR(I,J) + PRODR
            AI(I,J) = AI(I,J) + PRODI
 47         CONTINUE
 48       CONTINUE
 50       CONTINUE
 60     CONTINUE
 70   K = N
      IF (DABS(AR(N,N))+DABS(AI(N,N)) .EQ. 0.D0) GO TO 80
      RETURN
 80   IER = K
      IP(N) = 0
      RETURN
C----------------------- END OF SUBROUTINE DECC ------------------------
      END
C
C
      SUBROUTINE SOLC (N, NDIM, AR, AI, BR, BI, IP)
C VERSION COMPLEX DOUBLE PRECISION
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER N,NDIM,IP,NM1,K,KP1,M,I,KB,KM1
      DIMENSION AR(NDIM,N), AI(NDIM,N), BR(N), BI(N), IP(N)
C-----------------------------------------------------------------------
C  SOLUTION OF LINEAR SYSTEM, A*X = B .
C  INPUT..
C    N = ORDER OF MATRIX.
C    NDIM = DECLARED DIMENSION OF ARRAYS  AR AND AI.
C    (AR,AI) = TRIANGULARIZED MATRIX OBTAINED FROM DEC.
C    (BR,BI) = RIGHT HAND SIDE VECTOR.
C    IP = PIVOT VECTOR OBTAINED FROM DEC.
C  DO NOT USE IF DEC HAS SET IER .NE. 0.
C  OUTPUT..
C    (BR,BI) = SOLUTION VECTOR, X .
C-----------------------------------------------------------------------
      IF (N .EQ. 1) GO TO 50
      NM1 = N - 1
      DO 20 K = 1,NM1
        KP1 = K + 1
        M = IP(K)
        TR = BR(M)
        TI = BI(M)
        BR(M) = BR(K)
        BI(M) = BI(K)
        BR(K) = TR
        BI(K) = TI
        DO 10 I = KP1,N
          PRODR=AR(I,K)*TR-AI(I,K)*TI
          PRODI=AI(I,K)*TR+AR(I,K)*TI
          BR(I) = BR(I) + PRODR
          BI(I) = BI(I) + PRODI
 10       CONTINUE
 20     CONTINUE
      DO 40 KB = 1,NM1
        KM1 = N - KB
        K = KM1 + 1
        DEN=AR(K,K)*AR(K,K)+AI(K,K)*AI(K,K)
        PRODR=BR(K)*AR(K,K)+BI(K)*AI(K,K)
        PRODI=BI(K)*AR(K,K)-BR(K)*AI(K,K)
        BR(K)=PRODR/DEN
        BI(K)=PRODI/DEN
        TR = -BR(K)
        TI = -BI(K)
        DO 30 I = 1,KM1
          PRODR=AR(I,K)*TR-AI(I,K)*TI
          PRODI=AI(I,K)*TR+AR(I,K)*TI
          BR(I) = BR(I) + PRODR
          BI(I) = BI(I) + PRODI
 30       CONTINUE
 40     CONTINUE
 50     CONTINUE
        DEN=AR(1,1)*AR(1,1)+AI(1,1)*AI(1,1)
        PRODR=BR(1)*AR(1,1)+BI(1)*AI(1,1)
        PRODI=BI(1)*AR(1,1)-BR(1)*AI(1,1)
        BR(1)=PRODR/DEN
        BI(1)=PRODI/DEN
      RETURN
C----------------------- END OF SUBROUTINE SOLC ------------------------
      END  
C
C
      SUBROUTINE DECHC (N, NDIM, AR, AI, LB, IP, IER)
C VERSION COMPLEX DOUBLE PRECISION
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER N,NDIM,IP,IER,NM1,K,KP1,M,I,J
      DIMENSION AR(NDIM,N), AI(NDIM,N), IP(N)
C-----------------------------------------------------------------------
C  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION
C  ------ MODIFICATION FOR COMPLEX MATRICES --------
C  INPUT..
C     N = ORDER OF MATRIX.
C     NDIM = DECLARED DIMENSION OF ARRAYS  AR AND AI .
C     (AR, AI) = MATRIX TO BE TRIANGULARIZED.
C  OUTPUT..
C     AR(I,J), I.LE.J = UPPER TRIANGULAR FACTOR, U ; REAL PART.
C     AI(I,J), I.LE.J = UPPER TRIANGULAR FACTOR, U ; IMAGINARY PART.
C     AR(I,J), I.GT.J = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I - L.
C                                                    REAL PART.
C     AI(I,J), I.GT.J = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I - L.
C                                                    IMAGINARY PART.
C     LB = LOWER BANDWIDTH OF A (DIAGONAL NOT COUNTED), LB.GE.1.
C     IP(K), K.LT.N = INDEX OF K-TH PIVOT ROW.
C     IP(N) = (-1)**(NUMBER OF INTERCHANGES) OR O .
C     IER = 0 IF MATRIX A IS NONSINGULAR, OR K IF FOUND TO BE
C           SINGULAR AT STAGE K.
C  USE  SOL  TO OBTAIN SOLUTION OF LINEAR SYSTEM.
C  IF IP(N)=O, A IS SINGULAR, SOL WILL DIVIDE BY ZERO.
C
C  REFERENCE..
C     C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER,
C     C.A.C.M. 15 (1972), P. 274.
C-----------------------------------------------------------------------
      IER = 0
      IP(N) = 1
      IF (LB .EQ. 0) GO TO 70
      IF (N .EQ. 1) GO TO 70
      NM1 = N - 1
      DO 60 K = 1,NM1
        KP1 = K + 1
        M = K 
        NA = MIN0(N,LB+K)
        DO 10 I = KP1,NA
          IF (DABS(AR(I,K))+DABS(AI(I,K)) .GT.
     &          DABS(AR(M,K))+DABS(AI(M,K))) M = I
 10     CONTINUE
        IP(K) = M
        TR = AR(M,K)
        TI = AI(M,K)
        IF (M .EQ. K) GO TO 20
        IP(N) = -IP(N)
        AR(M,K) = AR(K,K)
        AI(M,K) = AI(K,K)
        AR(K,K) = TR
        AI(K,K) = TI
 20     CONTINUE
        IF (DABS(TR)+DABS(TI) .EQ. 0.D0) GO TO 80
        DEN=TR*TR+TI*TI
        TR=TR/DEN
        TI=-TI/DEN
        DO 30 I = KP1,NA
          PRODR=AR(I,K)*TR-AI(I,K)*TI
          PRODI=AI(I,K)*TR+AR(I,K)*TI
          AR(I,K)=-PRODR
          AI(I,K)=-PRODI
 30       CONTINUE
        DO 50 J = KP1,N
          TR = AR(M,J)
          TI = AI(M,J)
          AR(M,J) = AR(K,J)
          AI(M,J) = AI(K,J)
          AR(K,J) = TR
          AI(K,J) = TI
          IF (DABS(TR)+DABS(TI) .EQ. 0.D0) GO TO 48
          IF (TI .EQ. 0.D0) THEN
            DO 40 I = KP1,NA
            PRODR=AR(I,K)*TR
            PRODI=AI(I,K)*TR
            AR(I,J) = AR(I,J) + PRODR
            AI(I,J) = AI(I,J) + PRODI
 40         CONTINUE
            GO TO 48
          END IF
          IF (TR .EQ. 0.D0) THEN
            DO 45 I = KP1,NA
            PRODR=-AI(I,K)*TI
            PRODI=AR(I,K)*TI
            AR(I,J) = AR(I,J) + PRODR
            AI(I,J) = AI(I,J) + PRODI
 45         CONTINUE
            GO TO 48
          END IF
          DO 47 I = KP1,NA
            PRODR=AR(I,K)*TR-AI(I,K)*TI
            PRODI=AI(I,K)*TR+AR(I,K)*TI
            AR(I,J) = AR(I,J) + PRODR
            AI(I,J) = AI(I,J) + PRODI
 47         CONTINUE
 48       CONTINUE
 50       CONTINUE
 60     CONTINUE
 70   K = N
      IF (DABS(AR(N,N))+DABS(AI(N,N)) .EQ. 0.D0) GO TO 80
      RETURN
 80   IER = K
      IP(N) = 0
      RETURN
C----------------------- END OF SUBROUTINE DECHC -----------------------
      END
C
C
      SUBROUTINE SOLHC (N, NDIM, AR, AI, LB, BR, BI, IP)
C VERSION COMPLEX DOUBLE PRECISION
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER N,NDIM,IP,NM1,K,KP1,M,I,KB,KM1
      DIMENSION AR(NDIM,N), AI(NDIM,N), BR(N), BI(N), IP(N)
C-----------------------------------------------------------------------
C  SOLUTION OF LINEAR SYSTEM, A*X = B .
C  INPUT..
C    N = ORDER OF MATRIX.
C    NDIM = DECLARED DIMENSION OF ARRAYS  AR AND AI.
C    (AR,AI) = TRIANGULARIZED MATRIX OBTAINED FROM DEC.
C    (BR,BI) = RIGHT HAND SIDE VECTOR.
C    LB = LOWER BANDWIDTH OF A.
C    IP = PIVOT VECTOR OBTAINED FROM DEC.
C  DO NOT USE IF DEC HAS SET IER .NE. 0.
C  OUTPUT..
C    (BR,BI) = SOLUTION VECTOR, X .
C-----------------------------------------------------------------------
      IF (N .EQ. 1) GO TO 50
      NM1 = N - 1
      IF (LB .EQ. 0) GO TO 25
      DO 20 K = 1,NM1
        KP1 = K + 1
        M = IP(K)
        TR = BR(M)
        TI = BI(M)
        BR(M) = BR(K)
        BI(M) = BI(K)
        BR(K) = TR
        BI(K) = TI
        DO 10 I = KP1,MIN0(N,LB+K)
          PRODR=AR(I,K)*TR-AI(I,K)*TI
          PRODI=AI(I,K)*TR+AR(I,K)*TI
          BR(I) = BR(I) + PRODR
          BI(I) = BI(I) + PRODI
 10       CONTINUE
 20     CONTINUE
 25     CONTINUE
      DO 40 KB = 1,NM1
        KM1 = N - KB
        K = KM1 + 1
        DEN=AR(K,K)*AR(K,K)+AI(K,K)*AI(K,K)
        PRODR=BR(K)*AR(K,K)+BI(K)*AI(K,K)
        PRODI=BI(K)*AR(K,K)-BR(K)*AI(K,K)
        BR(K)=PRODR/DEN
        BI(K)=PRODI/DEN
        TR = -BR(K)
        TI = -BI(K)
        DO 30 I = 1,KM1
          PRODR=AR(I,K)*TR-AI(I,K)*TI
          PRODI=AI(I,K)*TR+AR(I,K)*TI
          BR(I) = BR(I) + PRODR
          BI(I) = BI(I) + PRODI
 30       CONTINUE
 40     CONTINUE
 50     CONTINUE
        DEN=AR(1,1)*AR(1,1)+AI(1,1)*AI(1,1)
        PRODR=BR(1)*AR(1,1)+BI(1)*AI(1,1)
        PRODI=BI(1)*AR(1,1)-BR(1)*AI(1,1)
        BR(1)=PRODR/DEN
        BI(1)=PRODI/DEN
      RETURN
C----------------------- END OF SUBROUTINE SOLHC -----------------------
      END  
C
      SUBROUTINE DECB (N, NDIM, A, ML, MU, IP, IER)
      REAL*8 A,T
      DIMENSION A(NDIM,N), IP(N)
C-----------------------------------------------------------------------
C  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION OF A BANDED
C  MATRIX WITH LOWER BANDWIDTH ML AND UPPER BANDWIDTH MU
C  INPUT..
C     N       ORDER OF THE ORIGINAL MATRIX A.
C     NDIM    DECLARED DIMENSION OF ARRAY  A.
C     A       CONTAINS THE MATRIX IN BAND STORAGE.   THE COLUMNS  
C                OF THE MATRIX ARE STORED IN THE COLUMNS OF  A  AND
C                THE DIAGONALS OF THE MATRIX ARE STORED IN ROWS 
C                ML+1 THROUGH 2*ML+MU+1 OF  A.
C     ML      LOWER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED).
C     MU      UPPER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED).
C  OUTPUT..
C     A       AN UPPER TRIANGULAR MATRIX IN BAND STORAGE AND 
C                THE MULTIPLIERS WHICH WERE USED TO OBTAIN IT.  
C     IP      INDEX VECTOR OF PIVOT INDICES.
C     IP(N)   (-1)**(NUMBER OF INTERCHANGES) OR O .
C     IER     = 0 IF MATRIX A IS NONSINGULAR, OR  = K IF FOUND TO BE
C                SINGULAR AT STAGE K.
C  USE  SOLB  TO OBTAIN SOLUTION OF LINEAR SYSTEM.
C  DETERM(A) = IP(N)*A(MD,1)*A(MD,2)*...*A(MD,N)  WITH MD=ML+MU+1.
C  IF IP(N)=O, A IS SINGULAR, SOLB WILL DIVIDE BY ZERO.
C
C  REFERENCE..
C     THIS IS A MODIFICATION OF
C     C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER,
C     C.A.C.M. 15 (1972), P. 274.
C-----------------------------------------------------------------------
      IER = 0
      IP(N) = 1 
      MD = ML + MU + 1
      MD1 = MD + 1
      JU = 0
      IF (ML .EQ. 0) GO TO 70
      IF (N .EQ. 1) GO TO 70
      IF (N .LT. MU+2) GO TO 7
      DO 5 J = MU+2,N
      DO 5 I = 1,ML
  5   A(I,J) = 0.D0
  7   NM1 = N - 1
      DO 60 K = 1,NM1
        KP1 = K + 1
        M = MD
        MDL = MIN(ML,N-K) + MD
        DO 10 I = MD1,MDL
          IF (DABS(A(I,K)) .GT. DABS(A(M,K))) M = I
 10     CONTINUE
        IP(K) = M + K - MD
        T = A(M,K)
        IF (M .EQ. MD) GO TO 20
        IP(N) = -IP(N)
        A(M,K) = A(MD,K)
        A(MD,K) = T
 20     CONTINUE
        IF (T .EQ. 0.D0) GO TO 80
        T = 1.D0/T
        DO 30 I = MD1,MDL
 30       A(I,K) = -A(I,K)*T 
        JU = MIN0(MAX0(JU,MU+IP(K)),N)
        MM = MD
        IF (JU .LT. KP1) GO TO 55
        DO 50 J = KP1,JU
          M = M - 1
          MM = MM - 1
          T = A(M,J) 
          IF (M .EQ. MM) GO TO 35
          A(M,J) = A(MM,J)
          A(MM,J) = T
 35       CONTINUE
          IF (T .EQ. 0.D0) GO TO 45
          JK = J - K
          DO 40 I = MD1,MDL
            IJK = I - JK
 40         A(IJK,J) = A(IJK,J) + A(I,K)*T
 45       CONTINUE
 50       CONTINUE
 55     CONTINUE
 60     CONTINUE
 70   K = N
      IF (A(MD,N) .EQ. 0.D0) GO TO 80
      RETURN
 80   IER = K
      IP(N) = 0
      RETURN
C----------------------- END OF SUBROUTINE DECB ------------------------
      END
C
C
      SUBROUTINE SOLB (N, NDIM, A, ML, MU, B, IP)
      REAL*8 A,B,T
      DIMENSION A(NDIM,N), B(N), IP(N)
C-----------------------------------------------------------------------
C  SOLUTION OF LINEAR SYSTEM, A*X = B .
C  INPUT..
C    N      ORDER OF MATRIX A.
C    NDIM   DECLARED DIMENSION OF ARRAY  A .
C    A      TRIANGULARIZED MATRIX OBTAINED FROM DECB.
C    ML     LOWER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED).
C    MU     UPPER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED).
C    B      RIGHT HAND SIDE VECTOR.
C    IP     PIVOT VECTOR OBTAINED FROM DECB.
C  DO NOT USE IF DECB HAS SET IER .NE. 0.
C  OUTPUT..
C    B      SOLUTION VECTOR, X .
C-----------------------------------------------------------------------
      MD = ML + MU + 1
      MD1 = MD + 1
      MDM = MD - 1
      NM1 = N - 1
      IF (ML .EQ. 0) GO TO 25
      IF (N .EQ. 1) GO TO 50
      DO 20 K = 1,NM1
        M = IP(K)
        T = B(M)
        B(M) = B(K)
        B(K) = T
        MDL = MIN(ML,N-K) + MD
        DO 10 I = MD1,MDL
          IMD = I + K - MD
 10       B(IMD) = B(IMD) + A(I,K)*T
 20     CONTINUE
 25   CONTINUE
      DO 40 KB = 1,NM1
        K = N + 1 - KB
        B(K) = B(K)/A(MD,K)
        T = -B(K) 
        KMD = MD - K
        LM = MAX0(1,KMD+1)
        DO 30 I = LM,MDM
          IMD = I - KMD
 30       B(IMD) = B(IMD) + A(I,K)*T
 40     CONTINUE
 50   B(1) = B(1)/A(MD,1)
      RETURN
C----------------------- END OF SUBROUTINE SOLB ------------------------
      END
C
      SUBROUTINE DECBC (N, NDIM, AR, AI, ML, MU, IP, IER)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION AR(NDIM,N), AI(NDIM,N), IP(N)
C-----------------------------------------------------------------------
C  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION OF A BANDED COMPLEX
C  MATRIX WITH LOWER BANDWIDTH ML AND UPPER BANDWIDTH MU
C  INPUT..
C     N       ORDER OF THE ORIGINAL MATRIX A.
C     NDIM    DECLARED DIMENSION OF ARRAY  A.
C     AR, AI     CONTAINS THE MATRIX IN BAND STORAGE.   THE COLUMNS  
C                OF THE MATRIX ARE STORED IN THE COLUMNS OF  AR (REAL
C                PART) AND AI (IMAGINARY PART)  AND
C                THE DIAGONALS OF THE MATRIX ARE STORED IN ROWS 
C                ML+1 THROUGH 2*ML+MU+1 OF  AR AND AI.
C     ML      LOWER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED).
C     MU      UPPER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED).
C  OUTPUT..
C     AR, AI  AN UPPER TRIANGULAR MATRIX IN BAND STORAGE AND 
C                THE MULTIPLIERS WHICH WERE USED TO OBTAIN IT.  
C     IP      INDEX VECTOR OF PIVOT INDICES.
C     IP(N)   (-1)**(NUMBER OF INTERCHANGES) OR O .
C     IER     = 0 IF MATRIX A IS NONSINGULAR, OR  = K IF FOUND TO BE
C                SINGULAR AT STAGE K.
C  USE  SOLBC  TO OBTAIN SOLUTION OF LINEAR SYSTEM.
C  DETERM(A) = IP(N)*A(MD,1)*A(MD,2)*...*A(MD,N)  WITH MD=ML+MU+1.
C  IF IP(N)=O, A IS SINGULAR, SOLBC WILL DIVIDE BY ZERO.
C
C  REFERENCE..
C     THIS IS A MODIFICATION OF
C     C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER,
C     C.A.C.M. 15 (1972), P. 274.
C-----------------------------------------------------------------------
      IER = 0
      IP(N) = 1 
      MD = ML + MU + 1
      MD1 = MD + 1
      JU = 0
      IF (ML .EQ. 0) GO TO 70
      IF (N .EQ. 1) GO TO 70
      IF (N .LT. MU+2) GO TO 7
      DO 5 J = MU+2,N
      DO 5 I = 1,ML
      AR(I,J) = 0.D0
      AI(I,J) = 0.D0
  5   CONTINUE
  7   NM1 = N - 1
      DO 60 K = 1,NM1
        KP1 = K + 1
        M = MD
        MDL = MIN(ML,N-K) + MD
        DO 10 I = MD1,MDL
          IF (DABS(AR(I,K))+DABS(AI(I,K)) .GT.
     &          DABS(AR(M,K))+DABS(AI(M,K))) M = I
 10     CONTINUE
        IP(K) = M + K - MD
        TR = AR(M,K)
        TI = AI(M,K)
        IF (M .EQ. MD) GO TO 20
        IP(N) = -IP(N)
        AR(M,K) = AR(MD,K)
        AI(M,K) = AI(MD,K)
        AR(MD,K) = TR
        AI(MD,K) = TI
 20     IF (DABS(TR)+DABS(TI) .EQ. 0.D0) GO TO 80
        DEN=TR*TR+TI*TI
        TR=TR/DEN
        TI=-TI/DEN
        DO 30 I = MD1,MDL
          PRODR=AR(I,K)*TR-AI(I,K)*TI
          PRODI=AI(I,K)*TR+AR(I,K)*TI
          AR(I,K)=-PRODR
          AI(I,K)=-PRODI
 30       CONTINUE
        JU = MIN0(MAX0(JU,MU+IP(K)),N)
        MM = MD
        IF (JU .LT. KP1) GO TO 55
        DO 50 J = KP1,JU
          M = M - 1
          MM = MM - 1
          TR = AR(M,J)
          TI = AI(M,J)
          IF (M .EQ. MM) GO TO 35
          AR(M,J) = AR(MM,J)
          AI(M,J) = AI(MM,J)
          AR(MM,J) = TR
          AI(MM,J) = TI
 35       CONTINUE
          IF (DABS(TR)+DABS(TI) .EQ. 0.D0) GO TO 48
          JK = J - K
          IF (TI .EQ. 0.D0) THEN
            DO 40 I = MD1,MDL
            IJK = I - JK
            PRODR=AR(I,K)*TR
            PRODI=AI(I,K)*TR
            AR(IJK,J) = AR(IJK,J) + PRODR
            AI(IJK,J) = AI(IJK,J) + PRODI
 40         CONTINUE
            GO TO 48
          END IF
          IF (TR .EQ. 0.D0) THEN
            DO 45 I = MD1,MDL
            IJK = I - JK
            PRODR=-AI(I,K)*TI
            PRODI=AR(I,K)*TI
            AR(IJK,J) = AR(IJK,J) + PRODR
            AI(IJK,J) = AI(IJK,J) + PRODI
 45         CONTINUE
            GO TO 48
          END IF
          DO 47 I = MD1,MDL
            IJK = I - JK
            PRODR=AR(I,K)*TR-AI(I,K)*TI
            PRODI=AI(I,K)*TR+AR(I,K)*TI
            AR(IJK,J) = AR(IJK,J) + PRODR
            AI(IJK,J) = AI(IJK,J) + PRODI
 47         CONTINUE
 48       CONTINUE
 50       CONTINUE
 55     CONTINUE
 60     CONTINUE
 70   K = N
      IF (DABS(AR(MD,N))+DABS(AI(MD,N)) .EQ. 0.D0) GO TO 80
      RETURN
 80   IER = K
      IP(N) = 0
      RETURN
C----------------------- END OF SUBROUTINE DECBC ------------------------
      END
C
C
      SUBROUTINE SOLBC (N, NDIM, AR, AI, ML, MU, BR, BI, IP)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION AR(NDIM,N), AI(NDIM,N), BR(N), BI(N), IP(N)
C-----------------------------------------------------------------------
C  SOLUTION OF LINEAR SYSTEM, A*X = B ,
C                  VERSION BANDED AND COMPLEX-DOUBLE PRECISION.
C  INPUT..
C    N      ORDER OF MATRIX A.
C    NDIM   DECLARED DIMENSION OF ARRAY  A .
C    AR, AI TRIANGULARIZED MATRIX OBTAINED FROM DECB (REAL AND IMAG. PART).
C    ML     LOWER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED).
C    MU     UPPER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED).
C    BR, BI RIGHT HAND SIDE VECTOR (REAL AND IMAG. PART).
C    IP     PIVOT VECTOR OBTAINED FROM DECBC.
C  DO NOT USE IF DECB HAS SET IER .NE. 0.
C  OUTPUT..
C    BR, BI SOLUTION VECTOR, X (REAL AND IMAG. PART).
C-----------------------------------------------------------------------
      MD = ML + MU + 1
      MD1 = MD + 1
      MDM = MD - 1
      NM1 = N - 1
      IF (ML .EQ. 0) GO TO 25
      IF (N .EQ. 1) GO TO 50
      DO 20 K = 1,NM1
        M = IP(K)
        TR = BR(M)
        TI = BI(M)
        BR(M) = BR(K)
        BI(M) = BI(K)
        BR(K) = TR
        BI(K) = TI
        MDL = MIN(ML,N-K) + MD
        DO 10 I = MD1,MDL
          IMD = I + K - MD
          PRODR=AR(I,K)*TR-AI(I,K)*TI
          PRODI=AI(I,K)*TR+AR(I,K)*TI
          BR(IMD) = BR(IMD) + PRODR
          BI(IMD) = BI(IMD) + PRODI
 10     CONTINUE
 20     CONTINUE
 25     CONTINUE
      DO 40 KB = 1,NM1
        K = N + 1 - KB
        DEN=AR(MD,K)*AR(MD,K)+AI(MD,K)*AI(MD,K)
        PRODR=BR(K)*AR(MD,K)+BI(K)*AI(MD,K)
        PRODI=BI(K)*AR(MD,K)-BR(K)*AI(MD,K)
        BR(K)=PRODR/DEN
        BI(K)=PRODI/DEN
        TR = -BR(K)
        TI = -BI(K)
        KMD = MD - K
        LM = MAX0(1,KMD+1)
        DO 30 I = LM,MDM
          IMD = I - KMD
          PRODR=AR(I,K)*TR-AI(I,K)*TI
          PRODI=AI(I,K)*TR+AR(I,K)*TI
          BR(IMD) = BR(IMD) + PRODR
          BI(IMD) = BI(IMD) + PRODI
 30       CONTINUE
 40     CONTINUE
        DEN=AR(MD,1)*AR(MD,1)+AI(MD,1)*AI(MD,1)
        PRODR=BR(1)*AR(MD,1)+BI(1)*AI(MD,1)
        PRODI=BI(1)*AR(MD,1)-BR(1)*AI(MD,1)
        BR(1)=PRODR/DEN
        BI(1)=PRODI/DEN
 50   CONTINUE
      RETURN
C----------------------- END OF SUBROUTINE SOLBC ------------------------
      END
c
C
      subroutine elmhes(nm,n,low,igh,a,int)
C
      integer i,j,m,n,la,nm,igh,kp1,low,mm1,mp1
      real*8 a(nm,n)
      real*8 x,y
      real*8 dabs
      integer int(igh)
C
C     this subroutine is a translation of the algol procedure elmhes,
C     num. math. 12, 349-368(1968) by martin and wilkinson.
C     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).
C
C     given a real general matrix, this subroutine
C     reduces a submatrix situated in rows and columns
C     low through igh to upper hessenberg form by
C     stabilized elementary similarity transformations.
C
C     on input:
C
C      nm must be set to the row dimension of two-dimensional
C        array parameters as declared in the calling program
C        dimension statement;
C
C      n is the order of the matrix;
C
C      low and igh are integers determined by the balancing
C        subroutine  balanc.      if  balanc  has not been used,
C        set low=1, igh=n;
C
C      a contains the input matrix.
C
C     on output:
C
C      a contains the hessenberg matrix.  the multipliers
C        which were used in the reduction are stored in the
C        remaining triangle under the hessenberg matrix;
C
C      int contains information on the rows and columns
C        interchanged in the reduction.
C        only elements low through igh are used.
C
C     questions and comments should be directed to b. s. garbow,
C     applied mathematics division, argonne national laboratory
C
C     ------------------------------------------------------------------
C
      la = igh - 1
      kp1 = low + 1
      if (la .lt. kp1) go to 200
C
      do 180 m = kp1, la
       mm1 = m - 1
       x = 0.0d0
       i = m
C
       do 100 j = m, igh
          if (dabs(a(j,mm1)) .le. dabs(x)) go to 100
          x = a(j,mm1)
          i = j
  100   continue
C
       int(m) = i
       if (i .eq. m) go to 130
C    :::::::::: interchange rows and columns of a ::::::::::
       do 110 j = mm1, n
          y = a(i,j)
          a(i,j) = a(m,j)
          a(m,j) = y
  110   continue
C
       do 120 j = 1, igh
          y = a(j,i)
          a(j,i) = a(j,m)
          a(j,m) = y
  120   continue
C    :::::::::: end interchange ::::::::::
  130   if (x .eq. 0.0d0) go to 180
       mp1 = m + 1
C
       do 160 i = mp1, igh
          y = a(i,mm1)
          if (y .eq. 0.0d0) go to 160
          y = y / x
          a(i,mm1) = y
C
          do 140 j = m, n
  140      a(i,j) = a(i,j) - y * a(m,j)
C
          do 150 j = 1, igh
  150      a(j,m) = a(j,m) + y * a(j,i)
C
  160   continue
C
  180 continue
C
  200 return
C    :::::::::: last card of elmhes ::::::::::
      end



c
c
c

C----------------------------start of radua related routines -----------------
      SUBROUTINE RADAU5(N,FCN,X,Y,XEND,H,
     &                  RTOL,ATOL,ITOL,
     &                  JAC ,IJAC,MLJAC,MUJAC,
     &                  MAS ,IMAS,MLMAS,MUMAS,
     &                  SOLOUT,IOUT,
     &                  WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)
C ----------------------------------------------------------
C     NUMERICAL SOLUTION OF A STIFF (OR DIFFERENTIAL ALGEBRAIC)
C     SYSTEM OF FIRST 0RDER ORDINARY DIFFERENTIAL EQUATIONS
C                     M*Y'=F(X,Y).
C     THE SYSTEM CAN BE (LINEARLY) IMPLICIT (MASS-MATRIX M .NE. I)
C     OR EXPLICIT (M=I).
C     THE METHOD USED IS AN IMPLICIT RUNGE-KUTTA METHOD (RADAU IIA)
C     OF ORDER 5 WITH STEP SIZE CONTROL AND CONTINUOUS OUTPUT.
C     CF. SECTION IV.8
C
C     AUTHORS: E. HAIRER AND G. WANNER
C              UNIVERSITE DE GENEVE, DEPT. DE MATHEMATIQUES
C              CH-1211 GENEVE 24, SWITZERLAND 
C              E-MAIL:  Ernst.Hairer@math.unige.ch
C                       Gerhard.Wanner@math.unige.ch
C     
C     THIS CODE IS PART OF THE BOOK:
C         E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
C         EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
C         SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS 14,
C         SPRINGER-VERLAG 1991, SECOND EDITION 1996.
C      
C     VERSION OF JULY 9, 1996
C     (latest small correction: January 18, 2002)
C
C     INPUT PARAMETERS  
C     ----------------  
C     N           DIMENSION OF THE SYSTEM 
C
C     FCN         NAME (EXTERNAL) OF SUBROUTINE COMPUTING THE
C                 VALUE OF F(X,Y):
C                    SUBROUTINE FCN(N,X,Y,F,RPAR,IPAR)
C                    DOUBLE PRECISION X,Y(N),F(N)
C                    F(1)=...   ETC.
C                 RPAR, IPAR (SEE BELOW)
C
C     X           INITIAL X-VALUE
C
C     Y(N)        INITIAL VALUES FOR Y
C
C     XEND        FINAL X-VALUE (XEND-X MAY BE POSITIVE OR NEGATIVE)
C
C     H           INITIAL STEP SIZE GUESS;
C                 FOR STIFF EQUATIONS WITH INITIAL TRANSIENT, 
C                 H=1.D0/(NORM OF F'), USUALLY 1.D-3 OR 1.D-5, IS GOOD.
C                 THIS CHOICE IS NOT VERY IMPORTANT, THE STEP SIZE IS
C                 QUICKLY ADAPTED. (IF H=0.D0, THE CODE PUTS H=1.D-6).
C
C     RTOL,ATOL   RELATIVE AND ABSOLUTE ERROR TOLERANCES. THEY
C                 CAN BE BOTH SCALARS OR ELSE BOTH VECTORS OF LENGTH N.
C
C     ITOL        SWITCH FOR RTOL AND ATOL:
C                   ITOL=0: BOTH RTOL AND ATOL ARE SCALARS.
C                     THE CODE KEEPS, ROUGHLY, THE LOCAL ERROR OF
C                     Y(I) BELOW RTOL*ABS(Y(I))+ATOL
C                   ITOL=1: BOTH RTOL AND ATOL ARE VECTORS.
C                     THE CODE KEEPS THE LOCAL ERROR OF Y(I) BELOW
C                     RTOL(I)*ABS(Y(I))+ATOL(I).
C
C     JAC         NAME (EXTERNAL) OF THE SUBROUTINE WHICH COMPUTES
C                 THE PARTIAL DERIVATIVES OF F(X,Y) WITH RESPECT TO Y
C                 (THIS ROUTINE IS ONLY CALLED IF IJAC=1; SUPPLY
C                 A DUMMY SUBROUTINE IN THE CASE IJAC=0).
C                 FOR IJAC=1, THIS SUBROUTINE MUST HAVE THE FORM
C                    SUBROUTINE JAC(N,X,Y,DFY,LDFY,RPAR,IPAR)
C                    DOUBLE PRECISION X,Y(N),DFY(LDFY,N)
C                    DFY(1,1)= ...
C                 LDFY, THE COLUMN-LENGTH OF THE ARRAY, IS
C                 FURNISHED BY THE CALLING PROGRAM.
C                 IF (MLJAC.EQ.N) THE JACOBIAN IS SUPPOSED TO
C                    BE FULL AND THE PARTIAL DERIVATIVES ARE
C                    STORED IN DFY AS
C                       DFY(I,J) = PARTIAL F(I) / PARTIAL Y(J)
C                 ELSE, THE JACOBIAN IS TAKEN AS BANDED AND
C                    THE PARTIAL DERIVATIVES ARE STORED
C                    DIAGONAL-WISE AS
C                       DFY(I-J+MUJAC+1,J) = PARTIAL F(I) / PARTIAL Y(J).
C
C     IJAC        SWITCH FOR THE COMPUTATION OF THE JACOBIAN:
C                    IJAC=0: JACOBIAN IS COMPUTED INTERNALLY BY FINITE
C                       DIFFERENCES, SUBROUTINE "JAC" IS NEVER CALLED.
C                    IJAC=1: JACOBIAN IS SUPPLIED BY SUBROUTINE JAC.
C
C     MLJAC       SWITCH FOR THE BANDED STRUCTURE OF THE JACOBIAN:
C                    MLJAC=N: JACOBIAN IS A FULL MATRIX. THE LINEAR
C                       ALGEBRA IS DONE BY FULL-MATRIX GAUSS-ELIMINATION.
C                    0<=MLJAC<N: MLJAC IS THE LOWER BANDWITH OF JACOBIAN 
C                       MATRIX (>= NUMBER OF NON-ZERO DIAGONALS BELOW
C                       THE MAIN DIAGONAL).
C
C     MUJAC       UPPER BANDWITH OF JACOBIAN  MATRIX (>= NUMBER OF NON-
C                 ZERO DIAGONALS ABOVE THE MAIN DIAGONAL).
C                 NEED NOT BE DEFINED IF MLJAC=N.
C
C     ----   MAS,IMAS,MLMAS, AND MUMAS HAVE ANALOG MEANINGS      -----
C     ----   FOR THE "MASS MATRIX" (THE MATRIX "M" OF SECTION IV.8): -
C
C     MAS         NAME (EXTERNAL) OF SUBROUTINE COMPUTING THE MASS-
C                 MATRIX M.
C                 IF IMAS=0, THIS MATRIX IS ASSUMED TO BE THE IDENTITY
C                 MATRIX AND NEEDS NOT TO BE DEFINED;
C                 SUPPLY A DUMMY SUBROUTINE IN THIS CASE.
C                 IF IMAS=1, THE SUBROUTINE MAS IS OF THE FORM
C                    SUBROUTINE MAS(N,AM,LMAS,RPAR,IPAR)
C                    DOUBLE PRECISION AM(LMAS,N)
C                    AM(1,1)= ....
C                    IF (MLMAS.EQ.N) THE MASS-MATRIX IS STORED
C                    AS FULL MATRIX LIKE
C                         AM(I,J) = M(I,J)
C                    ELSE, THE MATRIX IS TAKEN AS BANDED AND STORED
C                    DIAGONAL-WISE AS
C                         AM(I-J+MUMAS+1,J) = M(I,J).
C
C     IMAS       GIVES INFORMATION ON THE MASS-MATRIX:
C                    IMAS=0: M IS SUPPOSED TO BE THE IDENTITY
C                       MATRIX, MAS IS NEVER CALLED.
C                    IMAS=1: MASS-MATRIX  IS SUPPLIED.
C
C     MLMAS       SWITCH FOR THE BANDED STRUCTURE OF THE MASS-MATRIX:
C                    MLMAS=N: THE FULL MATRIX CASE. THE LINEAR
C                       ALGEBRA IS DONE BY FULL-MATRIX GAUSS-ELIMINATION.
C                    0<=MLMAS<N: MLMAS IS THE LOWER BANDWITH OF THE
C                       MATRIX (>= NUMBER OF NON-ZERO DIAGONALS BELOW
C                       THE MAIN DIAGONAL).
C                 MLMAS IS SUPPOSED TO BE .LE. MLJAC.
C
C     MUMAS       UPPER BANDWITH OF MASS-MATRIX (>= NUMBER OF NON-
C                 ZERO DIAGONALS ABOVE THE MAIN DIAGONAL).
C                 NEED NOT BE DEFINED IF MLMAS=N.
C                 MUMAS IS SUPPOSED TO BE .LE. MUJAC.
C
C     SOLOUT      NAME (EXTERNAL) OF SUBROUTINE PROVIDING THE
C                 NUMERICAL SOLUTION DURING INTEGRATION. 
C                 IF IOUT=1, IT IS CALLED AFTER EVERY SUCCESSFUL STEP.
C                 SUPPLY A DUMMY SUBROUTINE IF IOUT=0. 
C                 IT MUST HAVE THE FORM
C                    SUBROUTINE SOLOUT (NR,XOLD,X,Y,CONT,LRC,N,
C                                       RPAR,IPAR,IRTRN)
C                    DOUBLE PRECISION X,Y(N),CONT(LRC)
C                    ....  
C                 SOLOUT FURNISHES THE SOLUTION "Y" AT THE NR-TH
C                    GRID-POINT "X" (THEREBY THE INITIAL VALUE IS
C                    THE FIRST GRID-POINT).
C                 "XOLD" IS THE PRECEEDING GRID-POINT.
C                 "IRTRN" SERVES TO INTERRUPT THE INTEGRATION. IF IRTRN
C                    IS SET <0, RADAU5 RETURNS TO THE CALLING PROGRAM.
C           
C          -----  CONTINUOUS OUTPUT: -----
C                 DURING CALLS TO "SOLOUT", A CONTINUOUS SOLUTION
C                 FOR THE INTERVAL [XOLD,X] IS AVAILABLE THROUGH
C                 THE FUNCTION
C                        >>>   CONTR5(I,S,CONT,LRC)   <<<
C                 WHICH PROVIDES AN APPROXIMATION TO THE I-TH
C                 COMPONENT OF THE SOLUTION AT THE POINT S. THE VALUE
C                 S SHOULD LIE IN THE INTERVAL [XOLD,X].
C                 DO NOT CHANGE THE ENTRIES OF CONT(LRC), IF THE
C                 DENSE OUTPUT FUNCTION IS USED.
C
C     IOUT        SWITCH FOR CALLING THE SUBROUTINE SOLOUT:
C                    IOUT=0: SUBROUTINE IS NEVER CALLED
C                    IOUT=1: SUBROUTINE IS AVAILABLE FOR OUTPUT.
C
C     WORK        ARRAY OF WORKING SPACE OF LENGTH "LWORK".
C                 WORK(1), WORK(2),.., WORK(20) SERVE AS PARAMETERS
C                 FOR THE CODE. FOR STANDARD USE OF THE CODE
C                 WORK(1),..,WORK(20) MUST BE SET TO ZERO BEFORE
C                 CALLING. SEE BELOW FOR A MORE SOPHISTICATED USE.
C                 WORK(21),..,WORK(LWORK) SERVE AS WORKING SPACE
C                 FOR ALL VECTORS AND MATRICES.
C                 "LWORK" MUST BE AT LEAST
C                             N*(LJAC+LMAS+3*LE+12)+20
C                 WHERE
C                    LJAC=N              IF MLJAC=N (FULL JACOBIAN)
C                    LJAC=MLJAC+MUJAC+1  IF MLJAC<N (BANDED JAC.)
C                 AND                  
C                    LMAS=0              IF IMAS=0
C                    LMAS=N              IF IMAS=1 AND MLMAS=N (FULL)
C                    LMAS=MLMAS+MUMAS+1  IF MLMAS<N (BANDED MASS-M.)
C                 AND
C                    LE=N               IF MLJAC=N (FULL JACOBIAN)
C                    LE=2*MLJAC+MUJAC+1 IF MLJAC<N (BANDED JAC.)
C
C                 IN THE USUAL CASE WHERE THE JACOBIAN IS FULL AND THE
C                 MASS-MATRIX IS THE INDENTITY (IMAS=0), THE MINIMUM
C                 STORAGE REQUIREMENT IS 
C                             LWORK = 4*N*N+12*N+20.
C                 IF IWORK(9)=M1>0 THEN "LWORK" MUST BE AT LEAST
C                          N*(LJAC+12)+(N-M1)*(LMAS+3*LE)+20
C                 WHERE IN THE DEFINITIONS OF LJAC, LMAS AND LE THE
C                 NUMBER N CAN BE REPLACED BY N-M1.
C
C     LWORK       DECLARED LENGTH OF ARRAY "WORK".
C
C     IWORK       INTEGER WORKING SPACE OF LENGTH "LIWORK".
C                 IWORK(1),IWORK(2),...,IWORK(20) SERVE AS PARAMETERS
C                 FOR THE CODE. FOR STANDARD USE, SET IWORK(1),..,
C                 IWORK(20) TO ZERO BEFORE CALLING.
C                 IWORK(21),...,IWORK(LIWORK) SERVE AS WORKING AREA.
C                 "LIWORK" MUST BE AT LEAST 3*N+20.
C
C     LIWORK      DECLARED LENGTH OF ARRAY "IWORK".
C
C     RPAR, IPAR  REAL AND INTEGER PARAMETERS (OR PARAMETER ARRAYS) WHICH  
C                 CAN BE USED FOR COMMUNICATION BETWEEN YOUR CALLING
C                 PROGRAM AND THE FCN, JAC, MAS, SOLOUT SUBROUTINES. 
C
C ----------------------------------------------------------------------
C 
C     SOPHISTICATED SETTING OF PARAMETERS
C     -----------------------------------
C              SEVERAL PARAMETERS OF THE CODE ARE TUNED TO MAKE IT WORK 
C              WELL. THEY MAY BE DEFINED BY SETTING WORK(1),...
C              AS WELL AS IWORK(1),... DIFFERENT FROM ZERO.
C              FOR ZERO INPUT, THE CODE CHOOSES DEFAULT VALUES:
C
C    IWORK(1)  IF IWORK(1).NE.0, THE CODE TRANSFORMS THE JACOBIAN
C              MATRIX TO HESSENBERG FORM. THIS IS PARTICULARLY
C              ADVANTAGEOUS FOR LARGE SYSTEMS WITH FULL JACOBIAN.
C              IT DOES NOT WORK FOR BANDED JACOBIAN (MLJAC<N)
C              AND NOT FOR IMPLICIT SYSTEMS (IMAS=1).
C
C    IWORK(2)  THIS IS THE MAXIMAL NUMBER OF ALLOWED STEPS.
C              THE DEFAULT VALUE (FOR IWORK(2)=0) IS 100000.
C
C    IWORK(3)  THE MAXIMUM NUMBER OF NEWTON ITERATIONS FOR THE
C              SOLUTION OF THE IMPLICIT SYSTEM IN EACH STEP.
C              THE DEFAULT VALUE (FOR IWORK(3)=0) IS 7.
C
C    IWORK(4)  IF IWORK(4).EQ.0 THE EXTRAPOLATED COLLOCATION SOLUTION
C              IS TAKEN AS STARTING VALUE FOR NEWTON'S METHOD.
C              IF IWORK(4).NE.0 ZERO STARTING VALUES ARE USED.
C              THE LATTER IS RECOMMENDED IF NEWTON'S METHOD HAS
C              DIFFICULTIES WITH CONVERGENCE (THIS IS THE CASE WHEN
C              NSTEP IS LARGER THAN NACCPT + NREJCT; SEE OUTPUT PARAM.).
C              DEFAULT IS IWORK(4)=0.
C
C       THE FOLLOWING 3 PARAMETERS ARE IMPORTANT FOR
C       DIFFERENTIAL-ALGEBRAIC SYSTEMS OF INDEX > 1.
C       THE FUNCTION-SUBROUTINE SHOULD BE WRITTEN SUCH THAT
C       THE INDEX 1,2,3 VARIABLES APPEAR IN THIS ORDER. 
C       IN ESTIMATING THE ERROR THE INDEX 2 VARIABLES ARE
C       MULTIPLIED BY H, THE INDEX 3 VARIABLES BY H**2.
C
C    IWORK(5)  DIMENSION OF THE INDEX 1 VARIABLES (MUST BE > 0). FOR 
C              ODE'S THIS EQUALS THE DIMENSION OF THE SYSTEM.
C              DEFAULT IWORK(5)=N.
C
C    IWORK(6)  DIMENSION OF THE INDEX 2 VARIABLES. DEFAULT IWORK(6)=0.
C
C    IWORK(7)  DIMENSION OF THE INDEX 3 VARIABLES. DEFAULT IWORK(7)=0.
C
C    IWORK(8)  SWITCH FOR STEP SIZE STRATEGY
C              IF IWORK(8).EQ.1  MOD. PREDICTIVE CONTROLLER (GUSTAFSSON)
C              IF IWORK(8).EQ.2  CLASSICAL STEP SIZE CONTROL
C              THE DEFAULT VALUE (FOR IWORK(8)=0) IS IWORK(8)=1.
C              THE CHOICE IWORK(8).EQ.1 SEEMS TO PRODUCE SAFER RESULTS;
C              FOR SIMPLE PROBLEMS, THE CHOICE IWORK(8).EQ.2 PRODUCES
C              OFTEN SLIGHTLY FASTER RUNS
C
C       IF THE DIFFERENTIAL SYSTEM HAS THE SPECIAL STRUCTURE THAT
C            Y(I)' = Y(I+M2)   FOR  I=1,...,M1,
C       WITH M1 A MULTIPLE OF M2, A SUBSTANTIAL GAIN IN COMPUTERTIME
C       CAN BE ACHIEVED BY SETTING THE PARAMETERS IWORK(9) AND IWORK(10).
C       E.G., FOR SECOND ORDER SYSTEMS P'=V, V'=G(P,V), WHERE P AND V ARE 
C       VECTORS OF DIMENSION N/2, ONE HAS TO PUT M1=M2=N/2.
C       FOR M1>0 SOME OF THE INPUT PARAMETERS HAVE DIFFERENT MEANINGS:
C       - JAC: ONLY THE ELEMENTS OF THE NON-TRIVIAL PART OF THE
C              JACOBIAN HAVE TO BE STORED
C              IF (MLJAC.EQ.N-M1) THE JACOBIAN IS SUPPOSED TO BE FULL
C                 DFY(I,J) = PARTIAL F(I+M1) / PARTIAL Y(J)
C                FOR I=1,N-M1 AND J=1,N.
C              ELSE, THE JACOBIAN IS BANDED ( M1 = M2 * MM )
C                 DFY(I-J+MUJAC+1,J+K*M2) = PARTIAL F(I+M1) / PARTIAL Y(J+K*M2)
C                FOR I=1,MLJAC+MUJAC+1 AND J=1,M2 AND K=0,MM.
C       - MLJAC: MLJAC=N-M1: IF THE NON-TRIVIAL PART OF THE JACOBIAN IS FULL
C                0<=MLJAC<N-M1: IF THE (MM+1) SUBMATRICES (FOR K=0,MM)
C                     PARTIAL F(I+M1) / PARTIAL Y(J+K*M2),  I,J=1,M2
C                    ARE BANDED, MLJAC IS THE MAXIMAL LOWER BANDWIDTH
C                    OF THESE MM+1 SUBMATRICES
C       - MUJAC: MAXIMAL UPPER BANDWIDTH OF THESE MM+1 SUBMATRICES
C                NEED NOT BE DEFINED IF MLJAC=N-M1
C       - MAS: IF IMAS=0 THIS MATRIX IS ASSUMED TO BE THE IDENTITY AND
C              NEED NOT BE DEFINED. SUPPLY A DUMMY SUBROUTINE IN THIS CASE.
C              IT IS ASSUMED THAT ONLY THE ELEMENTS OF RIGHT LOWER BLOCK OF
C              DIMENSION N-M1 DIFFER FROM THAT OF THE IDENTITY MATRIX.
C              IF (MLMAS.EQ.N-M1) THIS SUBMATRIX IS SUPPOSED TO BE FULL
C                 AM(I,J) = M(I+M1,J+M1)     FOR I=1,N-M1 AND J=1,N-M1.
C              ELSE, THE MASS MATRIX IS BANDED
C                 AM(I-J+MUMAS+1,J) = M(I+M1,J+M1)
C       - MLMAS: MLMAS=N-M1: IF THE NON-TRIVIAL PART OF M IS FULL
C                0<=MLMAS<N-M1: LOWER BANDWIDTH OF THE MASS MATRIX
C       - MUMAS: UPPER BANDWIDTH OF THE MASS MATRIX
C                NEED NOT BE DEFINED IF MLMAS=N-M1
C
C    IWORK(9)  THE VALUE OF M1.  DEFAULT M1=0.
C
C    IWORK(10) THE VALUE OF M2.  DEFAULT M2=M1.
C
C ----------
C
C    WORK(1)   UROUND, THE ROUNDING UNIT, DEFAULT 1.D-16.
C
C    WORK(2)   THE SAFETY FACTOR IN STEP SIZE PREDICTION,
C              DEFAULT 0.9D0.
C
C    WORK(3)   DECIDES WHETHER THE JACOBIAN SHOULD BE RECOMPUTED;
C              INCREASE WORK(3), TO 0.1 SAY, WHEN JACOBIAN EVALUATIONS
C              ARE COSTLY. FOR SMALL SYSTEMS WORK(3) SHOULD BE SMALLER 
C              (0.001D0, SAY). NEGATIV WORK(3) FORCES THE CODE TO
C              COMPUTE THE JACOBIAN AFTER EVERY ACCEPTED STEP.     
C              DEFAULT 0.001D0.
C
C    WORK(4)   STOPPING CRITERION FOR NEWTON'S METHOD, USUALLY CHOSEN <1.
C              SMALLER VALUES OF WORK(4) MAKE THE CODE SLOWER, BUT SAFER.
C              DEFAULT MIN(0.03D0,RTOL(1)**0.5D0)
C
C    WORK(5) AND WORK(6) : IF WORK(5) < HNEW/HOLD < WORK(6), THEN THE
C              STEP SIZE IS NOT CHANGED. THIS SAVES, TOGETHER WITH A
C              LARGE WORK(3), LU-DECOMPOSITIONS AND COMPUTING TIME FOR
C              LARGE SYSTEMS. FOR SMALL SYSTEMS ONE MAY HAVE
C              WORK(5)=1.D0, WORK(6)=1.2D0, FOR LARGE FULL SYSTEMS
C              WORK(5)=0.99D0, WORK(6)=2.D0 MIGHT BE GOOD.
C              DEFAULTS WORK(5)=1.D0, WORK(6)=1.2D0 .
C
C    WORK(7)   MAXIMAL STEP SIZE, DEFAULT XEND-X.
C
C    WORK(8), WORK(9)   PARAMETERS FOR STEP SIZE SELECTION
C              THE NEW STEP SIZE IS CHOSEN SUBJECT TO THE RESTRICTION
C                 WORK(8) <= HNEW/HOLD <= WORK(9)
C              DEFAULT VALUES: WORK(8)=0.2D0, WORK(9)=8.D0
C
C-----------------------------------------------------------------------
C
C     OUTPUT PARAMETERS 
C     ----------------- 
C     X           X-VALUE FOR WHICH THE SOLUTION HAS BEEN COMPUTED
C                 (AFTER SUCCESSFUL RETURN X=XEND).
C
C     Y(N)        NUMERICAL SOLUTION AT X
C 
C     H           PREDICTED STEP SIZE OF THE LAST ACCEPTED STEP
C
C     IDID        REPORTS ON SUCCESSFULNESS UPON RETURN:
C                   IDID= 1  COMPUTATION SUCCESSFUL,
C                   IDID= 2  COMPUT. SUCCESSFUL (INTERRUPTED BY SOLOUT)
C                   IDID=-1  INPUT IS NOT CONSISTENT,
C                   IDID=-2  LARGER NMAX IS NEEDED,
C                   IDID=-3  STEP SIZE BECOMES TOO SMALL,
C                   IDID=-4  MATRIX IS REPEATEDLY SINGULAR.
C
C   IWORK(14)  NFCN    NUMBER OF FUNCTION EVALUATIONS (THOSE FOR NUMERICAL
C                      EVALUATION OF THE JACOBIAN ARE NOT COUNTED)  
C   IWORK(15)  NJAC    NUMBER OF JACOBIAN EVALUATIONS (EITHER ANALYTICALLY
C                      OR NUMERICALLY)
C   IWORK(16)  NSTEP   NUMBER OF COMPUTED STEPS
C   IWORK(17)  NACCPT  NUMBER OF ACCEPTED STEPS
C   IWORK(18)  NREJCT  NUMBER OF REJECTED STEPS (DUE TO ERROR TEST),
C                      (STEP REJECTIONS IN THE FIRST STEP ARE NOT COUNTED)
C   IWORK(19)  NDEC    NUMBER OF LU-DECOMPOSITIONS OF BOTH MATRICES
C   IWORK(20)  NSOL    NUMBER OF FORWARD-BACKWARD SUBSTITUTIONS, OF BOTH
C                      SYSTEMS; THE NSTEP FORWARD-BACKWARD SUBSTITUTIONS,
C                      NEEDED FOR STEP SIZE SELECTION, ARE NOT COUNTED
C-----------------------------------------------------------------------
C *** *** *** *** *** *** *** *** *** *** *** *** ***
C          DECLARATIONS 
C *** *** *** *** *** *** *** *** *** *** *** *** ***
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Y(N),ATOL(*),RTOL(*),WORK(LWORK),IWORK(LIWORK)
      DIMENSION RPAR(*),IPAR(*)
      LOGICAL IMPLCT,JBAND,ARRET,STARTN,PRED
      EXTERNAL FCN,JAC,MAS,SOLOUT
C *** *** *** *** *** *** ***
C        SETTING THE PARAMETERS 
C *** *** *** *** *** *** ***
       NFCN=0
       NJAC=0
       NSTEP=0
       NACCPT=0
       NREJCT=0
       NDEC=0
       NSOL=0
       ARRET=.FALSE.
C -------- UROUND   SMALLEST NUMBER SATISFYING 1.0D0+UROUND>1.0D0  
      IF (WORK(1).EQ.0.0D0) THEN
         UROUND=1.0D-16
      ELSE
         UROUND=WORK(1)
         IF (UROUND.LE.1.0D-19.OR.UROUND.GE.1.0D0) THEN
            WRITE(6,*)' COEFFICIENTS HAVE 20 DIGITS, UROUND=',WORK(1)
            ARRET=.TRUE.
         END IF
      END IF
C -------- CHECK AND CHANGE THE TOLERANCES
      EXPM=2.0D0/3.0D0
      IF (ITOL.EQ.0) THEN
          IF (ATOL(1).LE.0.D0.OR.RTOL(1).LE.10.D0*UROUND) THEN
              WRITE (6,*) ' TOLERANCES ARE TOO SMALL'
              ARRET=.TRUE.
          ELSE
              QUOT=ATOL(1)/RTOL(1)
              RTOL(1)=0.1D0*RTOL(1)**EXPM
              ATOL(1)=RTOL(1)*QUOT
          END IF
      ELSE
          DO I=1,N
          IF (ATOL(I).LE.0.D0.OR.RTOL(I).LE.10.D0*UROUND) THEN
              WRITE (6,*) ' TOLERANCES(',I,') ARE TOO SMALL'
              ARRET=.TRUE.
          ELSE
              QUOT=ATOL(I)/RTOL(I)
              RTOL(I)=0.1D0*RTOL(I)**EXPM
              ATOL(I)=RTOL(I)*QUOT
          END IF
          END DO
      END IF
C -------- NMAX , THE MAXIMAL NUMBER OF STEPS -----
      IF (IWORK(2).EQ.0) THEN
         NMAX=100000
      ELSE
         NMAX=IWORK(2)
         IF (NMAX.LE.0) THEN
            WRITE(6,*)' WRONG INPUT IWORK(2)=',IWORK(2)
            ARRET=.TRUE.
         END IF
      END IF
C -------- NIT    MAXIMAL NUMBER OF NEWTON ITERATIONS
      IF (IWORK(3).EQ.0) THEN
         NIT=7
      ELSE
         NIT=IWORK(3)
         IF (NIT.LE.0) THEN
            WRITE(6,*)' CURIOUS INPUT IWORK(3)=',IWORK(3)
            ARRET=.TRUE.
         END IF
      END IF
C -------- STARTN  SWITCH FOR STARTING VALUES OF NEWTON ITERATIONS
      IF(IWORK(4).EQ.0)THEN
         STARTN=.FALSE.
      ELSE
         STARTN=.TRUE.
      END IF
C -------- PARAMETER FOR DIFFERENTIAL-ALGEBRAIC COMPONENTS
      NIND1=IWORK(5)
      NIND2=IWORK(6)
      NIND3=IWORK(7)
      IF (NIND1.EQ.0) NIND1=N
      IF (NIND1+NIND2+NIND3.NE.N) THEN
       WRITE(6,*)' CURIOUS INPUT FOR IWORK(5,6,7)=',NIND1,NIND2,NIND3
       ARRET=.TRUE.
      END IF
C -------- PRED   STEP SIZE CONTROL
      IF(IWORK(8).LE.1)THEN
         PRED=.TRUE.
      ELSE
         PRED=.FALSE.
      END IF
C -------- PARAMETER FOR SECOND ORDER EQUATIONS
      M1=IWORK(9)
      M2=IWORK(10)
      NM1=N-M1
      IF (M1.EQ.0) M2=N
      IF (M2.EQ.0) M2=M1
      IF (M1.LT.0.OR.M2.LT.0.OR.M1+M2.GT.N) THEN
       WRITE(6,*)' CURIOUS INPUT FOR IWORK(9,10)=',M1,M2
       ARRET=.TRUE.
      END IF
C --------- SAFE     SAFETY FACTOR IN STEP SIZE PREDICTION
      IF (WORK(2).EQ.0.0D0) THEN
         SAFE=0.9D0
      ELSE
         SAFE=WORK(2)
         IF (SAFE.LE.0.001D0.OR.SAFE.GE.1.0D0) THEN
            WRITE(6,*)' CURIOUS INPUT FOR WORK(2)=',WORK(2)
            ARRET=.TRUE.
         END IF
      END IF
C ------ THET     DECIDES WHETHER THE JACOBIAN SHOULD BE RECOMPUTED;
      IF (WORK(3).EQ.0.D0) THEN
         THET=0.001D0
      ELSE
         THET=WORK(3)
         IF (THET.GE.1.0D0) THEN
            WRITE(6,*)' CURIOUS INPUT FOR WORK(3)=',WORK(3)
            ARRET=.TRUE.
         END IF
      END IF
C --- FNEWT   STOPPING CRITERION FOR NEWTON'S METHOD, USUALLY CHOSEN <1.
      TOLST=RTOL(1)
      IF (WORK(4).EQ.0.D0) THEN
         FNEWT=MAX(10*UROUND/TOLST,MIN(0.03D0,TOLST**0.5D0))
      ELSE
         FNEWT=WORK(4)
         IF (FNEWT.LE.UROUND/TOLST) THEN
            WRITE(6,*)' CURIOUS INPUT FOR WORK(4)=',WORK(4)
            ARRET=.TRUE.
         END IF
      END IF
C --- QUOT1 AND QUOT2: IF QUOT1 < HNEW/HOLD < QUOT2, STEP SIZE = CONST.
      IF (WORK(5).EQ.0.D0) THEN
         QUOT1=1.D0
      ELSE
         QUOT1=WORK(5)
      END IF
      IF (WORK(6).EQ.0.D0) THEN
         QUOT2=1.2D0
      ELSE
         QUOT2=WORK(6)
      END IF
      IF (QUOT1.GT.1.0D0.OR.QUOT2.LT.1.0D0) THEN
         WRITE(6,*)' CURIOUS INPUT FOR WORK(5,6)=',QUOT1,QUOT2
         ARRET=.TRUE.
      END IF
C -------- MAXIMAL STEP SIZE
      IF (WORK(7).EQ.0.D0) THEN
         HMAX=XEND-X
      ELSE
         HMAX=WORK(7)
      END IF 
C -------  FACL,FACR     PARAMETERS FOR STEP SIZE SELECTION
      IF(WORK(8).EQ.0.D0)THEN
         FACL=5.D0
      ELSE
         FACL=1.D0/WORK(8)
      END IF
      IF(WORK(9).EQ.0.D0)THEN
         FACR=1.D0/8.0D0
      ELSE
         FACR=1.D0/WORK(9)
      END IF
      IF (FACL.LT.1.0D0.OR.FACR.GT.1.0D0) THEN
            WRITE(6,*)' CURIOUS INPUT WORK(8,9)=',WORK(8),WORK(9)
            ARRET=.TRUE.
         END IF
C *** *** *** *** *** *** *** *** *** *** *** *** ***
C         COMPUTATION OF ARRAY ENTRIES
C *** *** *** *** *** *** *** *** *** *** *** *** ***
C ---- IMPLICIT, BANDED OR NOT ?
      IMPLCT=IMAS.NE.0
      JBAND=MLJAC.LT.NM1
C -------- COMPUTATION OF THE ROW-DIMENSIONS OF THE 2-ARRAYS ---
C -- JACOBIAN  AND  MATRICES E1, E2
      IF (JBAND) THEN
         LDJAC=MLJAC+MUJAC+1
         LDE1=MLJAC+LDJAC
      ELSE
         MLJAC=NM1
         MUJAC=NM1
         LDJAC=NM1
         LDE1=NM1
      END IF
C -- MASS MATRIX
      IF (IMPLCT) THEN
          IF (MLMAS.NE.NM1) THEN
              LDMAS=MLMAS+MUMAS+1
              IF (JBAND) THEN
                 IJOB=4
              ELSE
                 IJOB=3
              END IF
          ELSE
              MUMAS=NM1
              LDMAS=NM1
              IJOB=5
          END IF
C ------ BANDWITH OF "MAS" NOT SMALLER THAN BANDWITH OF "JAC"
          IF (MLMAS.GT.MLJAC.OR.MUMAS.GT.MUJAC) THEN
             WRITE (6,*) 'BANDWITH OF "MAS" NOT SMALLER THAN BANDWITH OF
     & "JAC"'
            ARRET=.TRUE.
          END IF
      ELSE
          LDMAS=0
          IF (JBAND) THEN
             IJOB=2
          ELSE
             IJOB=1
             IF (N.GT.2.AND.IWORK(1).NE.0) IJOB=7
          END IF
      END IF
      LDMAS2=MAX(1,LDMAS)
C ------ HESSENBERG OPTION ONLY FOR EXPLICIT EQU. WITH FULL JACOBIAN
      IF ((IMPLCT.OR.JBAND).AND.IJOB.EQ.7) THEN
         WRITE(6,*)' HESSENBERG OPTION ONLY FOR EXPLICIT EQUATIONS WITH 
     &FULL JACOBIAN'
         ARRET=.TRUE.
      END IF
C ------- PREPARE THE ENTRY-POINTS FOR THE ARRAYS IN WORK -----
      IEZ1=21
      IEZ2=IEZ1+N
      IEZ3=IEZ2+N
      IEY0=IEZ3+N
      IESCAL=IEY0+N
      IEF1=IESCAL+N
      IEF2=IEF1+N
      IEF3=IEF2+N
      IECON=IEF3+N
      IEJAC=IECON+4*N
      IEMAS=IEJAC+N*LDJAC
      IEE1=IEMAS+NM1*LDMAS
      IEE2R=IEE1+NM1*LDE1
      IEE2I=IEE2R+NM1*LDE1
C ------ TOTAL STORAGE REQUIREMENT -----------
      ISTORE=IEE2I+NM1*LDE1-1
      IF(ISTORE.GT.LWORK)THEN
         WRITE(6,*)' INSUFFICIENT STORAGE FOR WORK, MIN. LWORK=',ISTORE
         ARRET=.TRUE.
      END IF
C ------- ENTRY POINTS FOR INTEGER WORKSPACE -----
      IEIP1=21
      IEIP2=IEIP1+NM1
      IEIPH=IEIP2+NM1
C --------- TOTAL REQUIREMENT ---------------
      ISTORE=IEIPH+NM1-1
      IF (ISTORE.GT.LIWORK) THEN
         WRITE(6,*)' INSUFF. STORAGE FOR IWORK, MIN. LIWORK=',ISTORE
         ARRET=.TRUE.
      END IF
C ------ WHEN A FAIL HAS OCCURED, WE RETURN WITH IDID=-1
      IF (ARRET) THEN
         IDID=-1
         RETURN
      END IF
C -------- CALL TO CORE INTEGRATOR ------------
      CALL RADCOR(N,FCN,X,Y,XEND,HMAX,H,RTOL,ATOL,ITOL,
     &   JAC,IJAC,MLJAC,MUJAC,MAS,MLMAS,MUMAS,SOLOUT,IOUT,IDID,
     &   NMAX,UROUND,SAFE,THET,FNEWT,QUOT1,QUOT2,NIT,IJOB,STARTN,
     &   NIND1,NIND2,NIND3,PRED,FACL,FACR,M1,M2,NM1,
     &   IMPLCT,JBAND,LDJAC,LDE1,LDMAS2,WORK(IEZ1),WORK(IEZ2),
     &   WORK(IEZ3),WORK(IEY0),WORK(IESCAL),WORK(IEF1),WORK(IEF2),
     &   WORK(IEF3),WORK(IEJAC),WORK(IEE1),WORK(IEE2R),WORK(IEE2I),
     &   WORK(IEMAS),IWORK(IEIP1),IWORK(IEIP2),IWORK(IEIPH),
     &   WORK(IECON),NFCN,NJAC,NSTEP,NACCPT,NREJCT,NDEC,NSOL,RPAR,IPAR)
      IWORK(14)=NFCN
      IWORK(15)=NJAC
      IWORK(16)=NSTEP
      IWORK(17)=NACCPT
      IWORK(18)=NREJCT
      IWORK(19)=NDEC
      IWORK(20)=NSOL
C -------- RESTORE TOLERANCES
      EXPM=1.0D0/EXPM
      IF (ITOL.EQ.0) THEN
              QUOT=ATOL(1)/RTOL(1)
              RTOL(1)=(10.0D0*RTOL(1))**EXPM
              ATOL(1)=RTOL(1)*QUOT
      ELSE
          DO I=1,N
              QUOT=ATOL(I)/RTOL(I)
              RTOL(I)=(10.0D0*RTOL(I))**EXPM
              ATOL(I)=RTOL(I)*QUOT
          END DO
      END IF
C ----------- RETURN -----------
      RETURN
      END
C
C     END OF SUBROUTINE RADAU5
C
C ***********************************************************
C
      SUBROUTINE RADCOR(N,FCN,X,Y,XEND,HMAX,H,RTOL,ATOL,ITOL,
     &   JAC,IJAC,MLJAC,MUJAC,MAS,MLMAS,MUMAS,SOLOUT,IOUT,IDID,
     &   NMAX,UROUND,SAFE,THET,FNEWT,QUOT1,QUOT2,NIT,IJOB,STARTN,
     &   NIND1,NIND2,NIND3,PRED,FACL,FACR,M1,M2,NM1,
     &   IMPLCT,BANDED,LDJAC,LDE1,LDMAS,Z1,Z2,Z3,
     &   Y0,SCAL,F1,F2,F3,FJAC,E1,E2R,E2I,FMAS,IP1,IP2,IPHES,
     &   CONT,NFCN,NJAC,NSTEP,NACCPT,NREJCT,NDEC,NSOL,RPAR,IPAR)
C ----------------------------------------------------------
C     CORE INTEGRATOR FOR RADAU5
C     PARAMETERS SAME AS IN RADAU5 WITH WORKSPACE ADDED 
C ---------------------------------------------------------- 
C         DECLARATIONS 
C ---------------------------------------------------------- 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Y(N),Z1(N),Z2(N),Z3(N),Y0(N),SCAL(N),F1(N),F2(N),F3(N)
      DIMENSION FJAC(LDJAC,N),FMAS(LDMAS,NM1),CONT(4*N)
      DIMENSION E1(LDE1,NM1),E2R(LDE1,NM1),E2I(LDE1,NM1)
      DIMENSION ATOL(*),RTOL(*),RPAR(*),IPAR(*)
      INTEGER IP1(NM1),IP2(NM1),IPHES(NM1)
      COMMON /CONRA5/NN,NN2,NN3,NN4,XSOL,HSOL,C2M1,C1M1
      COMMON/LINAL/MLE,MUE,MBJAC,MBB,MDIAG,MDIFF,MBDIAG
      LOGICAL REJECT,FIRST,IMPLCT,BANDED,CALJAC,STARTN,CALHES
      LOGICAL INDEX1,INDEX2,INDEX3,LAST,PRED
      EXTERNAL FCN
C *** *** *** *** *** *** ***
C  INITIALISATIONS
C *** *** *** *** *** *** ***
C --------- DUPLIFY N FOR COMMON BLOCK CONT -----
      NN=N
      NN2=2*N
      NN3=3*N 
      LRC=4*N
C -------- CHECK THE INDEX OF THE PROBLEM ----- 
      INDEX1=NIND1.NE.0
      INDEX2=NIND2.NE.0
      INDEX3=NIND3.NE.0
C ------- COMPUTE MASS MATRIX FOR IMPLICIT CASE ----------
      IF (IMPLCT) CALL MAS(NM1,FMAS,LDMAS,RPAR,IPAR)
C ---------- CONSTANTS ---------
      SQ6=DSQRT(6.D0)
      C1=(4.D0-SQ6)/10.D0
      C2=(4.D0+SQ6)/10.D0
      C1M1=C1-1.D0
      C2M1=C2-1.D0
      C1MC2=C1-C2
      DD1=-(13.D0+7.D0*SQ6)/3.D0
      DD2=(-13.D0+7.D0*SQ6)/3.D0
      DD3=-1.D0/3.D0
      U1=(6.D0+81.D0**(1.D0/3.D0)-9.D0**(1.D0/3.D0))/30.D0
      ALPH=(12.D0-81.D0**(1.D0/3.D0)+9.D0**(1.D0/3.D0))/60.D0
      BETA=(81.D0**(1.D0/3.D0)+9.D0**(1.D0/3.D0))*DSQRT(3.D0)/60.D0
      CNO=ALPH**2+BETA**2
      U1=1.0D0/U1
      ALPH=ALPH/CNO
      BETA=BETA/CNO
      T11=9.1232394870892942792D-02
      T12=-0.14125529502095420843D0
      T13=-3.0029194105147424492D-02
      T21=0.24171793270710701896D0
      T22=0.20412935229379993199D0
      T23=0.38294211275726193779D0
      T31=0.96604818261509293619D0
      TI11=4.3255798900631553510D0
      TI12=0.33919925181580986954D0
      TI13=0.54177053993587487119D0
      TI21=-4.1787185915519047273D0
      TI22=-0.32768282076106238708D0
      TI23=0.47662355450055045196D0
      TI31=-0.50287263494578687595D0
      TI32=2.5719269498556054292D0
      TI33=-0.59603920482822492497D0
      IF (M1.GT.0) IJOB=IJOB+10
      POSNEG=SIGN(1.D0,XEND-X)
      HMAXN=MIN(ABS(HMAX),ABS(XEND-X)) 
      IF (ABS(H).LE.10.D0*UROUND) H=1.0D-6
      H=MIN(ABS(H),HMAXN)
      H=SIGN(H,POSNEG)
      HOLD=H
      REJECT=.FALSE.
      FIRST=.TRUE.
      LAST=.FALSE.
      IF ((X+H*1.0001D0-XEND)*POSNEG.GE.0.D0) THEN
         H=XEND-X
         LAST=.TRUE.
      END IF
      HOPT=H
      FACCON=1.D0
      CFAC=SAFE*(1+2*NIT)
      NSING=0
      XOLD=X
      IF (IOUT.NE.0) THEN
          IRTRN=1
          NRSOL=1
          XOSOL=XOLD
          XSOL=X
          DO I=1,N
             CONT(I)=Y(I)
          END DO
          NSOLU=N
          HSOL=HOLD
          CALL SOLOUT(NRSOL,XOSOL,XSOL,Y,CONT,LRC,NSOLU,
     &                RPAR,IPAR,IRTRN)
          IF (IRTRN.LT.0) GOTO 179
      END IF
      MLE=MLJAC
      MUE=MUJAC
      MBJAC=MLJAC+MUJAC+1
      MBB=MLMAS+MUMAS+1
      MDIAG=MLE+MUE+1
      MDIFF=MLE+MUE-MUMAS
      MBDIAG=MUMAS+1
      N2=2*N
      N3=3*N
      IF (ITOL.EQ.0) THEN
          DO I=1,N
             SCAL(I)=ATOL(1)+RTOL(1)*ABS(Y(I))
          END DO
      ELSE
          DO I=1,N
             SCAL(I)=ATOL(I)+RTOL(I)*ABS(Y(I))
          END DO
      END IF
      HHFAC=H
      CALL FCN(N,X,Y,Y0,RPAR,IPAR)
      NFCN=NFCN+1
C --- BASIC INTEGRATION STEP  
  10  CONTINUE
C *** *** *** *** *** *** ***
C  COMPUTATION OF THE JACOBIAN
C *** *** *** *** *** *** ***
      NJAC=NJAC+1
      IF (IJAC.EQ.0) THEN
C --- COMPUTE JACOBIAN MATRIX NUMERICALLY
         IF (BANDED) THEN
C --- JACOBIAN IS BANDED
            MUJACP=MUJAC+1
            MD=MIN(MBJAC,M2)
            DO MM=1,M1/M2+1
               DO K=1,MD
                  J=K+(MM-1)*M2
 12               F1(J)=Y(J)
                  F2(J)=DSQRT(UROUND*MAX(1.D-5,ABS(Y(J))))
                  Y(J)=Y(J)+F2(J)
                  J=J+MD
                  IF (J.LE.MM*M2) GOTO 12 
                  CALL FCN(N,X,Y,CONT,RPAR,IPAR)
                  J=K+(MM-1)*M2
                  J1=K
                  LBEG=MAX(1,J1-MUJAC)+M1
 14               LEND=MIN(M2,J1+MLJAC)+M1
                  Y(J)=F1(J)
                  MUJACJ=MUJACP-J1-M1
                  DO L=LBEG,LEND
                     FJAC(L+MUJACJ,J)=(CONT(L)-Y0(L))/F2(J) 
                  END DO
                  J=J+MD
                  J1=J1+MD
                  LBEG=LEND+1
                  IF (J.LE.MM*M2) GOTO 14
               END DO
            END DO
         ELSE
C --- JACOBIAN IS FULL
            DO I=1,N
               YSAFE=Y(I)
               DELT=DSQRT(UROUND*MAX(1.D-5,ABS(YSAFE)))
               Y(I)=YSAFE+DELT
               CALL FCN(N,X,Y,CONT,RPAR,IPAR)
               DO J=M1+1,N
                 FJAC(J-M1,I)=(CONT(J)-Y0(J))/DELT
               END DO
               Y(I)=YSAFE
            END DO
         END IF
      ELSE
C --- COMPUTE JACOBIAN MATRIX ANALYTICALLY
         CALL JAC(N,X,Y,FJAC,LDJAC,RPAR,IPAR)
      END IF
      CALJAC=.TRUE.
      CALHES=.TRUE.
  20  CONTINUE
C --- COMPUTE THE MATRICES E1 AND E2 AND THEIR DECOMPOSITIONS
      FAC1=U1/H
      ALPHN=ALPH/H
      BETAN=BETA/H
      CALL DECOMR(N,FJAC,LDJAC,FMAS,LDMAS,MLMAS,MUMAS,
     &            M1,M2,NM1,FAC1,E1,LDE1,IP1,IER,IJOB,CALHES,IPHES)
      IF (IER.NE.0) GOTO 78
      CALL DECOMC(N,FJAC,LDJAC,FMAS,LDMAS,MLMAS,MUMAS,
     &            M1,M2,NM1,ALPHN,BETAN,E2R,E2I,LDE1,IP2,IER,IJOB)
      IF (IER.NE.0) GOTO 78
      NDEC=NDEC+1
  30  CONTINUE
      NSTEP=NSTEP+1
      IF (NSTEP.GT.NMAX) GOTO 178
      IF (0.1D0*ABS(H).LE.ABS(X)*UROUND) GOTO 177
          IF (INDEX2) THEN
             DO I=NIND1+1,NIND1+NIND2
                SCAL(I)=SCAL(I)/HHFAC
             END DO
          END IF
          IF (INDEX3) THEN
             DO I=NIND1+NIND2+1,NIND1+NIND2+NIND3
                SCAL(I)=SCAL(I)/(HHFAC*HHFAC)
             END DO
          END IF
      XPH=X+H
C *** *** *** *** *** *** ***
C  STARTING VALUES FOR NEWTON ITERATION
C *** *** *** *** *** *** ***
      IF (FIRST.OR.STARTN) THEN
         DO I=1,N
            Z1(I)=0.D0
            Z2(I)=0.D0
            Z3(I)=0.D0
            F1(I)=0.D0
            F2(I)=0.D0
            F3(I)=0.D0
         END DO
      ELSE
         C3Q=H/HOLD
         C1Q=C1*C3Q
         C2Q=C2*C3Q
         DO I=1,N
            AK1=CONT(I+N)
            AK2=CONT(I+N2)
            AK3=CONT(I+N3)
            Z1I=C1Q*(AK1+(C1Q-C2M1)*(AK2+(C1Q-C1M1)*AK3))
            Z2I=C2Q*(AK1+(C2Q-C2M1)*(AK2+(C2Q-C1M1)*AK3))
            Z3I=C3Q*(AK1+(C3Q-C2M1)*(AK2+(C3Q-C1M1)*AK3))
            Z1(I)=Z1I
            Z2(I)=Z2I
            Z3(I)=Z3I
            F1(I)=TI11*Z1I+TI12*Z2I+TI13*Z3I
            F2(I)=TI21*Z1I+TI22*Z2I+TI23*Z3I
            F3(I)=TI31*Z1I+TI32*Z2I+TI33*Z3I
         END DO
      END IF
C *** *** *** *** *** *** ***
C  LOOP FOR THE SIMPLIFIED NEWTON ITERATION
C *** *** *** *** *** *** ***
            NEWT=0
            FACCON=MAX(FACCON,UROUND)**0.8D0
            THETA=ABS(THET)
  40        CONTINUE
            IF (NEWT.GE.NIT) GOTO 78
C ---     COMPUTE THE RIGHT-HAND SIDE
            DO I=1,N
               CONT(I)=Y(I)+Z1(I)
            END DO
            CALL FCN(N,X+C1*H,CONT,Z1,RPAR,IPAR)
            DO I=1,N
               CONT(I)=Y(I)+Z2(I)
            END DO
            CALL FCN(N,X+C2*H,CONT,Z2,RPAR,IPAR)
            DO I=1,N
               CONT(I)=Y(I)+Z3(I)
            END DO
            CALL FCN(N,XPH,CONT,Z3,RPAR,IPAR)
            NFCN=NFCN+3
C ---     SOLVE THE LINEAR SYSTEMS
           DO I=1,N
              A1=Z1(I)
              A2=Z2(I)
              A3=Z3(I)
              Z1(I)=TI11*A1+TI12*A2+TI13*A3
              Z2(I)=TI21*A1+TI22*A2+TI23*A3
              Z3(I)=TI31*A1+TI32*A2+TI33*A3
           END DO
        CALL SLVRAD(N,FJAC,LDJAC,MLJAC,MUJAC,FMAS,LDMAS,MLMAS,MUMAS,
     &          M1,M2,NM1,FAC1,ALPHN,BETAN,E1,E2R,E2I,LDE1,Z1,Z2,Z3,
     &          F1,F2,F3,CONT,IP1,IP2,IPHES,IER,IJOB)
            NSOL=NSOL+1
            NEWT=NEWT+1
            DYNO=0.D0
            DO I=1,N
               DENOM=SCAL(I)
               DYNO=DYNO+(Z1(I)/DENOM)**2+(Z2(I)/DENOM)**2
     &          +(Z3(I)/DENOM)**2
            END DO
            DYNO=DSQRT(DYNO/N3)
C ---     BAD CONVERGENCE OR NUMBER OF ITERATIONS TO LARGE
            IF (NEWT.GT.1.AND.NEWT.LT.NIT) THEN
                THQ=DYNO/DYNOLD
                IF (NEWT.EQ.2) THEN
                   THETA=THQ
                ELSE
                   THETA=SQRT(THQ*THQOLD)
                END IF
                THQOLD=THQ
                IF (THETA.LT.0.99D0) THEN
                    FACCON=THETA/(1.0D0-THETA)
                    DYTH=FACCON*DYNO*THETA**(NIT-1-NEWT)/FNEWT
                    IF (DYTH.GE.1.0D0) THEN
                         QNEWT=DMAX1(1.0D-4,DMIN1(20.0D0,DYTH))
                         HHFAC=.8D0*QNEWT**(-1.0D0/(4.0D0+NIT-1-NEWT))
                         H=HHFAC*H
                         REJECT=.TRUE.
                         LAST=.FALSE.
                         IF (CALJAC) GOTO 20
                         GOTO 10
                    END IF
                ELSE
                    GOTO 78
                END IF
            END IF
            DYNOLD=MAX(DYNO,UROUND)
            DO I=1,N
               F1I=F1(I)+Z1(I)
               F2I=F2(I)+Z2(I)
               F3I=F3(I)+Z3(I)
               F1(I)=F1I
               F2(I)=F2I
               F3(I)=F3I
               Z1(I)=T11*F1I+T12*F2I+T13*F3I
               Z2(I)=T21*F1I+T22*F2I+T23*F3I
               Z3(I)=T31*F1I+    F2I
            END DO
            IF (FACCON*DYNO.GT.FNEWT) GOTO 40
C --- ERROR ESTIMATION  
      CALL ESTRAD (N,FJAC,LDJAC,MLJAC,MUJAC,FMAS,LDMAS,MLMAS,MUMAS,
     &          H,DD1,DD2,DD3,FCN,NFCN,Y0,Y,IJOB,X,M1,M2,NM1,
     &          E1,LDE1,Z1,Z2,Z3,CONT,F1,F2,IP1,IPHES,SCAL,ERR,
     &          FIRST,REJECT,FAC1,RPAR,IPAR)
C --- COMPUTATION OF HNEW
C --- WE REQUIRE .2<=HNEW/H<=8.
      FAC=MIN(SAFE,CFAC/(NEWT+2*NIT))
      QUOT=MAX(FACR,MIN(FACL,ERR**.25D0/FAC))
      HNEW=H/QUOT
C *** *** *** *** *** *** ***
C  IS THE ERROR SMALL ENOUGH ?
C *** *** *** *** *** *** ***
      IF (ERR.LT.1.D0) THEN
C --- STEP IS ACCEPTED  
         FIRST=.FALSE.
         NACCPT=NACCPT+1
         IF (PRED) THEN
C       --- PREDICTIVE CONTROLLER OF GUSTAFSSON
            IF (NACCPT.GT.1) THEN
               FACGUS=(HACC/H)*(ERR**2/ERRACC)**0.25D0/SAFE
               FACGUS=MAX(FACR,MIN(FACL,FACGUS))
               QUOT=MAX(QUOT,FACGUS)
               HNEW=H/QUOT
            END IF
            HACC=H
            ERRACC=MAX(1.0D-2,ERR)
         END IF
         XOLD=X
         HOLD=H
         X=XPH 
         DO I=1,N
            Y(I)=Y(I)+Z3(I)  
            Z2I=Z2(I)
            Z1I=Z1(I)
            CONT(I+N)=(Z2I-Z3(I))/C2M1
            AK=(Z1I-Z2I)/C1MC2
            ACONT3=Z1I/C1
            ACONT3=(AK-ACONT3)/C2
            CONT(I+N2)=(AK-CONT(I+N))/C1M1
            CONT(I+N3)=CONT(I+N2)-ACONT3
         END DO
         IF (ITOL.EQ.0) THEN
             DO I=1,N
                SCAL(I)=ATOL(1)+RTOL(1)*ABS(Y(I))
             END DO
         ELSE
             DO I=1,N
                SCAL(I)=ATOL(I)+RTOL(I)*ABS(Y(I))
             END DO
         END IF
         IF (IOUT.NE.0) THEN
             NRSOL=NACCPT+1
             XSOL=X
             XOSOL=XOLD
             DO I=1,N
                CONT(I)=Y(I)
             END DO
             NSOLU=N
             HSOL=HOLD
             CALL SOLOUT(NRSOL,XOSOL,XSOL,Y,CONT,LRC,NSOLU,
     &                   RPAR,IPAR,IRTRN)
             IF (IRTRN.LT.0) GOTO 179
         END IF
         CALJAC=.FALSE.
         IF (LAST) THEN
            H=HOPT
            IDID=1
            RETURN
         END IF
         CALL FCN(N,X,Y,Y0,RPAR,IPAR)
         NFCN=NFCN+1
         HNEW=POSNEG*MIN(ABS(HNEW),HMAXN)
         HOPT=HNEW
         HOPT=MIN(H,HNEW)
         IF (REJECT) HNEW=POSNEG*MIN(ABS(HNEW),ABS(H)) 
         REJECT=.FALSE.
         IF ((X+HNEW/QUOT1-XEND)*POSNEG.GE.0.D0) THEN
            H=XEND-X
            LAST=.TRUE.
         ELSE
            QT=HNEW/H 
            HHFAC=H
            IF (THETA.LE.THET.AND.QT.GE.QUOT1.AND.QT.LE.QUOT2) GOTO 30
            H=HNEW 
         END IF
         HHFAC=H
         IF (THETA.LE.THET) GOTO 20
         GOTO 10
      ELSE
C --- STEP IS REJECTED  
         REJECT=.TRUE.
         LAST=.FALSE.
         IF (FIRST) THEN
             H=H*0.1D0
             HHFAC=0.1D0
         ELSE 
             HHFAC=HNEW/H
             H=HNEW
         END IF
         IF (NACCPT.GE.1) NREJCT=NREJCT+1
         IF (CALJAC) GOTO 20
         GOTO 10
      END IF
C --- UNEXPECTED STEP-REJECTION
  78  CONTINUE
      IF (IER.NE.0) THEN
          NSING=NSING+1
          IF (NSING.GE.5) GOTO 176
      END IF
      H=H*0.5D0 
      HHFAC=0.5D0
      REJECT=.TRUE.
      LAST=.FALSE.
      IF (CALJAC) GOTO 20
      GOTO 10
C --- FAIL EXIT
 176  CONTINUE
      WRITE(6,979)X   
      WRITE(6,*) ' MATRIX IS REPEATEDLY SINGULAR, IER=',IER
      IDID=-4
      RETURN
 177  CONTINUE
      WRITE(6,979)X   
      WRITE(6,*) ' STEP SIZE T0O SMALL, H=',H
      IDID=-3
      RETURN
 178  CONTINUE
      WRITE(6,979)X   
      WRITE(6,*) ' MORE THAN NMAX =',NMAX,'STEPS ARE NEEDED' 
      IDID=-2
      RETURN
C --- EXIT CAUSED BY SOLOUT
 179  CONTINUE
      WRITE(6,979)X
 979  FORMAT(' EXIT OF RADAU5 AT X=',E18.4) 
      IDID=2
      RETURN
      END
C
C     END OF SUBROUTINE RADCOR
C
C ***********************************************************
C
      DOUBLE PRECISION FUNCTION CONTR5(I,X,CONT,LRC) 
C ----------------------------------------------------------
C     THIS FUNCTION CAN BE USED FOR CONINUOUS OUTPUT. IT PROVIDES AN
C     APPROXIMATION TO THE I-TH COMPONENT OF THE SOLUTION AT X.
C     IT GIVES THE VALUE OF THE COLLOCATION POLYNOMIAL, DEFINED FOR
C     THE LAST SUCCESSFULLY COMPUTED STEP (BY RADAU5).
C ----------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION CONT(LRC)
      COMMON /CONRA5/NN,NN2,NN3,NN4,XSOL,HSOL,C2M1,C1M1
      S=(X-XSOL)/HSOL
      CONTR5=CONT(I)+S*(CONT(I+NN)+(S-C2M1)*(CONT(I+NN2)
     &     +(S-C1M1)*CONT(I+NN3)))
      RETURN
      END
C
C     END OF FUNCTION CONTR5
C
C ***********************************************************




c-------------------end radau5 related routines -----------------------------






*DECK DLSODI
      SUBROUTINE DLSODI (RES, ADDA, JAC, NEQ, Y, YDOTI, T, TOUT, ITOL,
     1  RTOL, ATOL, ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, LIW, MF )
      EXTERNAL RES, ADDA, JAC
      INTEGER NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, IWORK, LIW, MF
      DOUBLE PRECISION Y, YDOTI, T, TOUT, RTOL, ATOL, RWORK
      DIMENSION NEQ(*), Y(*), YDOTI(*), RTOL(*), ATOL(*), RWORK(LRW),
     1          IWORK(LIW)
C-----------------------------------------------------------------------
C This is the 7 May 2001 version of
C DLSODI: Livermore Solver for Ordinary Differential Equations
C         (Implicit form).
C
C This version is in double precision.
C
C DLSODI solves the initial value problem for linearly implicit
C systems of first order ODEs,
C     A(t,y) * dy/dt = g(t,y) ,  where A(t,y) is a square matrix,
C or, in component form,
C     ( a   * ( dy / dt ))  + ... +  ( a     * ( dy   / dt ))  =
C        i,1      1                     i,NEQ      NEQ
C
C      =   g ( t, y , y ,..., y    )   ( i = 1,...,NEQ )
C           i      1   2       NEQ
C
C If A is singular, this is a differential-algebraic system.
C
C DLSODI is a variant version of the DLSODE package.
C-----------------------------------------------------------------------
C Reference:
C     Alan C. Hindmarsh,  ODEPACK, A Systematized Collection of ODE
C     Solvers, in Scientific Computing, R. S. Stepleman et al. (Eds.),
C     North-Holland, Amsterdam, 1983, pp. 55-64.
C-----------------------------------------------------------------------
C Authors:       Alan C. Hindmarsh and Jeffrey F. Painter
C                Center for Applied Scientific Computing, L-561
C                Lawrence Livermore National Laboratory
C                Livermore, CA 94551
C-----------------------------------------------------------------------
C Summary of Usage.
C
C Communication between the user and the DLSODI package, for normal
C situations, is summarized here.  This summary describes only a subset
C of the full set of options available.  See the full description for
C details, including optional communication, nonstandard options,
C and instructions for special situations.  See also the example
C problem (with program and output) following this summary.
C
C A. First, provide a subroutine of the form:
C               SUBROUTINE RES (NEQ, T, Y, S, R, IRES)
C               DOUBLE PRECISION T, Y(*), S(*), R(*)
C which computes the residual function
C     r = g(t,y)  -  A(t,y) * s ,
C as a function of t and the vectors y and s.  (s is an internally
C generated approximation to dy/dt.)  The arrays Y and S are inputs
C to the RES routine and should not be altered.  The residual
C vector is to be stored in the array R.  The argument IRES should be
C ignored for casual use of DLSODI.  (For uses of IRES, see the
C paragraph on RES in the full description below.)
C
C B. Next, decide whether full or banded form is more economical
C for the storage of matrices.  DLSODI must deal internally with the
C matrices A and dr/dy, where r is the residual function defined above.
C DLSODI generates a linear combination of these two matrices, and
C this is treated in either full or banded form.
C     The matrix structure is communicated by a method flag MF,
C which is 21 or 22 for the full case, and 24 or 25 in the band case.
C     In the banded case, DLSODI requires two half-bandwidth
C parameters ML and MU.  These are, respectively, the widths of the
C lower and upper parts of the band, excluding the main diagonal.
C Thus the band consists of the locations (i,j) with
C i-ML .le. j .le. i+MU, and the full bandwidth is ML+MU+1.
C Note that the band must accommodate the nonzero elements of
C A(t,y), dg/dy, and d(A*s)/dy (s fixed).  Alternatively, one
C can define a band that encloses only the elements that are relatively
C large in magnitude, and gain some economy in storage and possibly
C also efficiency, although the appropriate threshhold for
C retaining matrix elements is highly problem-dependent.
C
C C. You must also provide a subroutine of the form:
C               SUBROUTINE ADDA (NEQ, T, Y, ML, MU, P, NROWP)
C               DOUBLE PRECISION T, Y(*), P(NROWP,*)
C which adds the matrix A = A(t,y) to the contents of the array P.
C T and the Y array are input and should not be altered.
C     In the full matrix case, this routine should add elements of
C to P in the usual order.  I.e., add A(i,j) to P(i,j).  (Ignore the
C ML and MU arguments in this case.)
C     In the band matrix case, this routine should add element A(i,j)
C to P(i-j+MU+1,j).  I.e., add the diagonal lines of A to the rows of
C P from the top down (the top line of A added to the first row of P).
C
C D. For the sake of efficiency, you are encouraged to supply the
C Jacobian matrix dr/dy in closed form, where r = g(t,y) - A(t,y)*s
C (s = a fixed vector) as above.  If dr/dy is being supplied,
C use MF = 21 or 24, and provide a subroutine of the form:
C               SUBROUTINE JAC (NEQ, T, Y, S, ML, MU, P, NROWP)
C               DOUBLE PRECISION T, Y(*), S(*), P(NROWP,*)
C which computes dr/dy as a function of t, y, and s.  Here T, Y, and
C S are inputs, and the routine is to load dr/dy into P as follows:
C     In the full matrix case (MF = 21), load P(i,j) with dr(i)/dy(j),
C the partial derivative of r(i) with respect to y(j).  (Ignore the
C ML and MU arguments in this case.)
C     In the band matrix case (MF = 24), load P(i-j+mu+1,j) with
C dr(i)/dy(j), i.e. load the diagonal lines of dr/dy into the rows of
C P from the top down.
C     In either case, only nonzero elements need be loaded, and the
C indexing of P is the same as in the ADDA routine.
C     Note that if A is independent of y (or this dependence
C is weak enough to be ignored) then JAC is to compute dg/dy.
C     If it is not feasible to provide a JAC routine, use
C MF = 22 or 25, and DLSODI will compute an approximate Jacobian
C internally by difference quotients.
C
C E. Next decide whether or not to provide the initial value of the
C derivative vector dy/dt.  If the initial value of A(t,y) is
C nonsingular (and not too ill-conditioned), you may let DLSODI compute
C this vector (ISTATE = 0).  (DLSODI will solve the system A*s = g for
C s, with initial values of A and g.)  If A(t,y) is initially
C singular, then the system is a differential-algebraic system, and
C you must make use of the particular form of the system to compute the
C initial values of y and dy/dt.  In that case, use ISTATE = 1 and
C load the initial value of dy/dt into the array YDOTI.
C The input array YDOTI and the initial Y array must be consistent with
C the equations A*dy/dt = g.  This implies that the initial residual
C r = g(t,y) - A(t,y)*YDOTI  must be approximately zero.
C
C F. Write a main program which calls Subroutine DLSODI once for
C each point at which answers are desired.  This should also provide
C for possible use of logical unit 6 for output of error messages
C by DLSODI.  On the first call to DLSODI, supply arguments as follows:
C RES    = name of user subroutine for residual function r.
C ADDA   = name of user subroutine for computing and adding A(t,y).
C JAC    = name of user subroutine for Jacobian matrix dr/dy
C          (MF = 21 or 24).  If not used, pass a dummy name.
C Note: the names for the RES and ADDA routines and (if used) the
C        JAC routine must be declared External in the calling program.
C NEQ    = number of scalar equations in the system.
C Y      = array of initial values, of length NEQ.
C YDOTI  = array of length NEQ (containing initial dy/dt if ISTATE = 1).
C T      = the initial value of the independent variable.
C TOUT   = first point where output is desired (.ne. T).
C ITOL   = 1 or 2 according as ATOL (below) is a scalar or array.
C RTOL   = relative tolerance parameter (scalar).
C ATOL   = absolute tolerance parameter (scalar or array).
C          the estimated local error in y(i) will be controlled so as
C          to be roughly less (in magnitude) than
C             EWT(i) = RTOL*ABS(Y(i)) + ATOL     if ITOL = 1, or
C             EWT(i) = RTOL*ABS(Y(i)) + ATOL(i)  if ITOL = 2.
C          Thus the local error test passes if, in each component,
C          either the absolute error is less than ATOL (or ATOL(i)),
C          or the relative error is less than RTOL.
C          Use RTOL = 0.0 for pure absolute error control, and
C          use ATOL = 0.0 (or ATOL(i) = 0.0) for pure relative error
C          control.  Caution: Actual (global) errors may exceed these
C          local tolerances, so choose them conservatively.
C ITASK  = 1 for normal computation of output values of y at t = TOUT.
C ISTATE = integer flag (input and output).  Set ISTATE = 1 if the
C          initial dy/dt is supplied, and 0 otherwise.
C IOPT   = 0 to indicate no optional inputs used.
C RWORK  = real work array of length at least:
C             22 +  9*NEQ + NEQ**2           for MF = 21 or 22,
C             22 + 10*NEQ + (2*ML + MU)*NEQ  for MF = 24 or 25.
C LRW    = declared length of RWORK (in user's dimension).
C IWORK  = integer work array of length at least 20 + NEQ.
C          If MF = 24 or 25, input in IWORK(1),IWORK(2) the lower
C          and upper half-bandwidths ML,MU.
C LIW    = declared length of IWORK (in user's dimension).
C MF     = method flag.  Standard values are:
C          21 for a user-supplied full Jacobian.
C          22 for an internally generated full Jacobian.
C          24 for a user-supplied banded Jacobian.
C          25 for an internally generated banded Jacobian.
C          for other choices of MF, see the paragraph on MF in
C          the full description below.
C Note that the main program must declare arrays Y, YDOTI, RWORK, IWORK,
C and possibly ATOL.
C
C G. The output from the first call (or any call) is:
C      Y = array of computed values of y(t) vector.
C      T = corresponding value of independent variable (normally TOUT).
C ISTATE = 2  if DLSODI was successful, negative otherwise.
C          -1 means excess work done on this call (check all inputs).
C          -2 means excess accuracy requested (tolerances too small).
C          -3 means illegal input detected (see printed message).
C          -4 means repeated error test failures (check all inputs).
C          -5 means repeated convergence failures (perhaps bad Jacobian
C             supplied or wrong choice of tolerances).
C          -6 means error weight became zero during problem. (Solution
C             component i vanished, and ATOL or ATOL(i) = 0.)
C          -7 cannot occur in casual use.
C          -8 means DLSODI was unable to compute the initial dy/dt.
C             In casual use, this means A(t,y) is initially singular.
C             Supply YDOTI and use ISTATE = 1 on the first call.
C
C  If DLSODI returns ISTATE = -1, -4, or -5, then the output of
C  DLSODI also includes YDOTI = array containing residual vector
C  r = g - A * dy/dt  evaluated at the current t, y, and dy/dt.
C
C H. To continue the integration after a successful return, simply
C reset TOUT and call DLSODI again.  No other parameters need be reset.
C
C-----------------------------------------------------------------------
C Example Problem.
C
C The following is a simple example problem, with the coding
C needed for its solution by DLSODI.  The problem is from chemical
C kinetics, and consists of the following three equations:
C     dy1/dt = -.04*y1 + 1.e4*y2*y3
C     dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e7*y2**2
C       0.   = y1 + y2 + y3 - 1.
C on the interval from t = 0.0 to t = 4.e10, with initial conditions
C y1 = 1.0, y2 = y3 = 0.
C
C The following coding solves this problem with DLSODI, using MF = 21
C and printing results at t = .4, 4., ..., 4.e10.  It uses
C ITOL = 2 and ATOL much smaller for y2 than y1 or y3 because
C y2 has much smaller values.  dy/dt is supplied in YDOTI. We had
C obtained the initial value of dy3/dt by differentiating the
C third equation and evaluating the first two at t = 0.
C At the end of the run, statistical quantities of interest are
C printed (see optional outputs in the full description below).
C
C     EXTERNAL RESID, APLUSP, DGBYDY
C     DOUBLE PRECISION ATOL, RTOL, RWORK, T, TOUT, Y, YDOTI
C     DIMENSION Y(3), YDOTI(3), ATOL(3), RWORK(58), IWORK(23)
C     NEQ = 3
C     Y(1) = 1.
C     Y(2) = 0.
C     Y(3) = 0.
C     YDOTI(1) = -.04
C     YDOTI(2) =  .04
C     YDOTI(3) =  0.
C     T = 0.
C     TOUT = .4
C     ITOL = 2
C     RTOL = 1.D-4
C     ATOL(1) = 1.D-6
C     ATOL(2) = 1.D-10
C     ATOL(3) = 1.D-6
C     ITASK = 1
C     ISTATE = 1
C     IOPT = 0
C     LRW = 58
C     LIW = 23
C     MF = 21
C     DO 40  IOUT = 1,12
C       CALL DLSODI(RESID, APLUSP, DGBYDY, NEQ, Y, YDOTI, T, TOUT, ITOL,
C    1     RTOL, ATOL, ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, LIW, MF)
C       WRITE (6,20)  T, Y(1), Y(2), Y(3)
C  20   FORMAT(' At t =',D12.4,'   Y =',3D14.6)
C       IF (ISTATE .LT. 0 )  GO TO 80
C  40   TOUT = TOUT*10.
C     WRITE (6,60)  IWORK(11), IWORK(12), IWORK(13)
C  60 FORMAT(/' No. steps =',I4,'  No. r-s =',I4,'  No. J-s =',I4)
C     STOP
C  80 WRITE (6,90)  ISTATE
C  90 FORMAT(///' Error halt.. ISTATE =',I3)
C     STOP
C     END
C
C     SUBROUTINE RESID(NEQ, T, Y, S, R, IRES)
C     DOUBLE PRECISION T, Y, S, R
C     DIMENSION Y(3), S(3), R(3)
C     R(1) = -.04*Y(1) + 1.D4*Y(2)*Y(3) - S(1)
C     R(2) = .04*Y(1) - 1.D4*Y(2)*Y(3) - 3.D7*Y(2)*Y(2) - S(2)
C     R(3) = Y(1) + Y(2) + Y(3) - 1.
C     RETURN
C     END
C
C     SUBROUTINE APLUSP(NEQ, T, Y, ML, MU, P, NROWP)
C     DOUBLE PRECISION T, Y, P
C     DIMENSION Y(3), P(NROWP,3)
C     P(1,1) = P(1,1) + 1.
C     P(2,2) = P(2,2) + 1.
C     RETURN
C     END
C
C     SUBROUTINE DGBYDY(NEQ, T, Y, S, ML, MU, P, NROWP)
C     DOUBLE PRECISION T, Y, S, P
C     DIMENSION Y(3), S(3), P(NROWP,3)
C     P(1,1) = -.04
C     P(1,2) = 1.D4*Y(3)
C     P(1,3) = 1.D4*Y(2)
C     P(2,1) = .04
C     P(2,2) = -1.D4*Y(3) - 6.D7*Y(2)
C     P(2,3) = -1.D4*Y(2)
C     P(3,1) = 1.
C     P(3,2) = 1.
C     P(3,3) = 1.
C     RETURN
C     END
C
C The output of this program (on a CDC-7600 in single precision)
C is as follows:
C
C   At t =  4.0000e-01   Y =  9.851726e-01  3.386406e-05  1.479357e-02
C   At t =  4.0000e+00   Y =  9.055142e-01  2.240418e-05  9.446344e-02
C   At t =  4.0000e+01   Y =  7.158050e-01  9.184616e-06  2.841858e-01
C   At t =  4.0000e+02   Y =  4.504846e-01  3.222434e-06  5.495122e-01
C   At t =  4.0000e+03   Y =  1.831701e-01  8.940379e-07  8.168290e-01
C   At t =  4.0000e+04   Y =  3.897016e-02  1.621193e-07  9.610297e-01
C   At t =  4.0000e+05   Y =  4.935213e-03  1.983756e-08  9.950648e-01
C   At t =  4.0000e+06   Y =  5.159269e-04  2.064759e-09  9.994841e-01
C   At t =  4.0000e+07   Y =  5.306413e-05  2.122677e-10  9.999469e-01
C   At t =  4.0000e+08   Y =  5.494532e-06  2.197826e-11  9.999945e-01
C   At t =  4.0000e+09   Y =  5.129457e-07  2.051784e-12  9.999995e-01
C   At t =  4.0000e+10   Y = -7.170472e-08 -2.868188e-13  1.000000e+00
C
C   No. steps = 330  No. r-s = 404  No. J-s =  69
C
C-----------------------------------------------------------------------
C Full Description of User Interface to DLSODI.
C
C The user interface to DLSODI consists of the following parts.
C
C 1.   The call sequence to Subroutine DLSODI, which is a driver
C      routine for the solver.  This includes descriptions of both
C      the call sequence arguments and of user-supplied routines.
C      Following these descriptions is a description of
C      optional inputs available through the call sequence, and then
C      a description of optional outputs (in the work arrays).
C
C 2.   Descriptions of other routines in the DLSODI package that may be
C      (optionally) called by the user.  These provide the ability to
C      alter error message handling, save and restore the internal
C      Common, and obtain specified derivatives of the solution y(t).
C
C 3.   Descriptions of Common blocks to be declared in overlay
C      or similar environments, or to be saved when doing an interrupt
C      of the problem and continued solution later.
C
C 4.   Description of two routines in the DLSODI package, either of
C      which the user may replace with his/her own version, if desired.
C      These relate to the measurement of errors.
C
C-----------------------------------------------------------------------
C Part 1.  Call Sequence.
C
C The call sequence parameters used for input only are
C     RES, ADDA, JAC, NEQ, TOUT, ITOL, RTOL, ATOL, ITASK,
C     IOPT, LRW, LIW, MF,
C and those used for both input and output are
C     Y, T, ISTATE, YDOTI.
C The work arrays RWORK and IWORK are also used for conditional and
C optional inputs and optional outputs.  (The term output here refers
C to the return from Subroutine DLSODI to the user's calling program.)
C
C The legality of input parameters will be thoroughly checked on the
C initial call for the problem, but not checked thereafter unless a
C change in input parameters is flagged by ISTATE = 3 on input.
C
C The descriptions of the call arguments are as follows.
C
C RES    = the name of the user-supplied subroutine which supplies
C          the residual vector for the ODE system, defined by
C            r = g(t,y) - A(t,y) * s
C          as a function of the scalar t and the vectors
C          s and y (s approximates dy/dt).  This subroutine
C          is to have the form
C               SUBROUTINE RES (NEQ, T, Y, S, R, IRES)
C               DOUBLE PRECISION T, Y(*), S(*), R(*)
C          where NEQ, T, Y, S, and IRES are input, and R and
C          IRES are output.  Y, S, and R are arrays of length NEQ.
C             On input, IRES indicates how DLSODI will use the
C          returned array R, as follows:
C             IRES = 1  means that DLSODI needs the full residual,
C                       r = g - A*s, exactly.
C             IRES = -1 means that DLSODI is using R only to compute
C                       the Jacobian dr/dy by difference quotients.
C          The RES routine can ignore IRES, or it can omit some terms
C          if IRES = -1.  If A does not depend on y, then RES can
C          just return R = g when IRES = -1.  If g - A*s contains other
C          additive terms that are independent of y, these can also be
C          dropped, if done consistently, when IRES = -1.
C             The subroutine should set the flag IRES if it
C          encounters a halt condition or illegal input.
C          Otherwise, it should not reset IRES.  On output,
C             IRES = 1 or -1 represents a normal return, and
C          DLSODI continues integrating the ODE.  Leave IRES
C          unchanged from its input value.
C             IRES = 2 tells DLSODI to immediately return control
C          to the calling program, with ISTATE = 3.  This lets
C          the calling program change parameters of the problem,
C          if necessary.
C             IRES = 3 represents an error condition (for example, an
C          illegal value of y).  DLSODI tries to integrate the system
C          without getting IRES = 3 from RES.  If it cannot, DLSODI
C          returns with ISTATE = -7 or -1.
C             On an DLSODI return with ISTATE = 3, -1, or -7, the values
C          of T and Y returned correspond to the last point reached
C          successfully without getting the flag IRES = 2 or 3.
C             The flag values IRES = 2 and 3 should not be used to
C          handle switches or root-stop conditions.  This is better
C          done by calling DLSODI in a one-step mode and checking the
C          stopping function for a sign change at each step.
C             If quantities computed in the RES routine are needed
C          externally to DLSODI, an extra call to RES should be made
C          for this purpose, for consistent and accurate results.
C          To get the current dy/dt for the S argument, use DINTDY.
C             RES must be declared External in the calling
C          program.  See note below for more about RES.
C
C ADDA   = the name of the user-supplied subroutine which adds the
C          matrix A = A(t,y) to another matrix stored in the same form
C          as A.  The storage form is determined by MITER (see MF).
C          This subroutine is to have the form
C               SUBROUTINE ADDA (NEQ, T, Y, ML, MU, P, NROWP)
C               DOUBLE PRECISION T, Y(*), P(NROWP,*)
C          where NEQ, T, Y, ML, MU, and NROWP are input and P is
C          output.  Y is an array of length NEQ, and the matrix P is
C          stored in an NROWP by NEQ array.
C             In the full matrix case ( MITER = 1 or 2) ADDA should
C          add  A    to P(i,j).  ML and MU are ignored.
C                i,j
C             In the band matrix case ( MITER = 4 or 5) ADDA should
C          add  A    to  P(i-j+MU+1,j).
C                i,j
C          See JAC for details on this band storage form.
C             ADDA must be declared External in the calling program.
C          See note below for more information about ADDA.
C
C JAC    = the name of the user-supplied subroutine which supplies the
C          Jacobian matrix, dr/dy, where r = g - A*s.  The form of the
C          Jacobian matrix is determined by MITER.  JAC is required
C          if MITER = 1 or 4 -- otherwise a dummy name can be
C          passed.  This subroutine is to have the form
C               SUBROUTINE JAC ( NEQ, T, Y, S, ML, MU, P, NROWP )
C               DOUBLE PRECISION T, Y(*), S(*), P(NROWP,*)
C          where NEQ, T, Y, S, ML, MU, and NROWP are input and P
C          is output.  Y and S are arrays of length NEQ, and the
C          matrix P is stored in an NROWP by NEQ array.
C          P is to be loaded with partial derivatives (elements
C          of the Jacobian matrix) on output.
C             In the full matrix case (MITER = 1), ML and MU
C          are ignored and the Jacobian is to be loaded into P
C          by columns-- i.e., dr(i)/dy(j) is loaded into P(i,j).
C             In the band matrix case (MITER = 4), the elements
C          within the band are to be loaded into P by columns,
C          with diagonal lines of dr/dy loaded into the
C          rows of P.  Thus dr(i)/dy(j) is to be loaded
C          into P(i-j+MU+1,j).  The locations in P in the two
C          triangular areas which correspond to nonexistent matrix
C          elements can be ignored or loaded arbitrarily, as they
C          they are overwritten by DLSODI.  ML and MU are the
C          half-bandwidth parameters (see IWORK).
C               In either case, P is preset to zero by the solver,
C          so that only the nonzero elements need be loaded by JAC.
C          Each call to JAC is preceded by a call to RES with the same
C          arguments NEQ, T, Y, and S.  Thus to gain some efficiency,
C          intermediate quantities shared by both calculations may be
C          saved in a user Common block by RES and not recomputed by JAC
C          if desired.  Also, JAC may alter the Y array, if desired.
C               JAC need not provide dr/dy exactly.  A crude
C          approximation (possibly with a smaller bandwidth) will do.
C               JAC must be declared External in the calling program.
C               See note below for more about JAC.
C
C    Note on RES, ADDA, and JAC:
C          These subroutines may access user-defined quantities in
C          NEQ(2),... and/or in Y(NEQ(1)+1),... if NEQ is an array
C          (dimensioned in the subroutines) and/or Y has length
C          exceeding NEQ(1).  However, these routines should not alter
C          NEQ(1), Y(1),...,Y(NEQ) or any other input variables.
C          See the descriptions of NEQ and Y below.
C
C NEQ    = the size of the system (number of first order ordinary
C          differential equations or scalar algebraic equations).
C          Used only for input.
C          NEQ may be decreased, but not increased, during the problem.
C          If NEQ is decreased (with ISTATE = 3 on input), the
C          remaining components of Y should be left undisturbed, if
C          these are to be accessed in RES, ADDA, or JAC.
C
C          Normally, NEQ is a scalar, and it is generally referred to
C          as a scalar in this user interface description.  However,
C          NEQ may be an array, with NEQ(1) set to the system size.
C          (The DLSODI package accesses only NEQ(1).)  In either case,
C          this parameter is passed as the NEQ argument in all calls
C          to RES, ADDA, and JAC.  Hence, if it is an array,
C          locations NEQ(2),... may be used to store other integer data
C          and pass it to RES, ADDA, or JAC.  Each such subroutine
C          must include NEQ in a Dimension statement in that case.
C
C Y      = a real array for the vector of dependent variables, of
C          length NEQ or more.  Used for both input and output on the
C          first call (ISTATE = 0 or 1), and only for output on other
C          calls.  On the first call, Y must contain the vector of
C          initial values.  On output, Y contains the computed solution
C          vector, evaluated at T.  If desired, the Y array may be used
C          for other purposes between calls to the solver.
C
C          This array is passed as the Y argument in all calls to RES,
C          ADDA, and JAC.  Hence its length may exceed NEQ,
C          and locations Y(NEQ+1),... may be used to store other real
C          data and pass it to RES, ADDA, or JAC.  (The DLSODI
C          package accesses only Y(1),...,Y(NEQ). )
C
C YDOTI  = a real array for the initial value of the vector
C          dy/dt and for work space, of dimension at least NEQ.
C
C          On input:
C            If ISTATE = 0, then DLSODI will compute the initial value
C          of dy/dt, if A is nonsingular.  Thus YDOTI will
C          serve only as work space and may have any value.
C            If ISTATE = 1, then YDOTI must contain the initial value
C          of dy/dt.
C            If ISTATE = 2 or 3 (continuation calls), then YDOTI
C          may have any value.
C            Note: If the initial value of A is singular, then
C          DLSODI cannot compute the initial value of dy/dt, so
C          it must be provided in YDOTI, with ISTATE = 1.
C
C          On output, when DLSODI terminates abnormally with ISTATE =
C          -1, -4, or -5, YDOTI will contain the residual
C          r = g(t,y) - A(t,y)*(dy/dt).  If r is large, t is near
C          its initial value, and YDOTI is supplied with ISTATE = 1,
C          then there may have been an incorrect input value of
C          YDOTI = dy/dt, or the problem (as given to DLSODI)
C          may not have a solution.
C
C          If desired, the YDOTI array may be used for other
C          purposes between calls to the solver.
C
C T      = the independent variable.  On input, T is used only on the
C          first call, as the initial point of the integration.
C          On output, after each call, T is the value at which a
C          computed solution Y is evaluated (usually the same as TOUT).
C          on an error return, T is the farthest point reached.
C
C TOUT   = the next value of t at which a computed solution is desired.
C          Used only for input.
C
C          When starting the problem (ISTATE = 0 or 1), TOUT may be
C          equal to T for one call, then should .ne. T for the next
C          call.  For the initial T, an input value of TOUT .ne. T is
C          used in order to determine the direction of the integration
C          (i.e. the algebraic sign of the step sizes) and the rough
C          scale of the problem.  Integration in either direction
C          (forward or backward in t) is permitted.
C
C          If ITASK = 2 or 5 (one-step modes), TOUT is ignored after
C          the first call (i.e. the first call with TOUT .ne. T).
C          Otherwise, TOUT is required on every call.
C
C          If ITASK = 1, 3, or 4, the values of TOUT need not be
C          monotone, but a value of TOUT which backs up is limited
C          to the current internal T interval, whose endpoints are
C          TCUR - HU and TCUR (see optional outputs, below, for
C          TCUR and HU).
C
C ITOL   = an indicator for the type of error control.  See
C          description below under ATOL.  Used only for input.
C
C RTOL   = a relative error tolerance parameter, either a scalar or
C          an array of length NEQ.  See description below under ATOL.
C          Input only.
C
C ATOL   = an absolute error tolerance parameter, either a scalar or
C          an array of length NEQ.  Input only.
C
C             The input parameters ITOL, RTOL, and ATOL determine
C          the error control performed by the solver.  The solver will
C          control the vector E = (E(i)) of estimated local errors
C          in y, according to an inequality of the form
C                      RMS-norm of ( E(i)/EWT(i) )   .le.   1,
C          where       EWT(i) = RTOL(i)*ABS(Y(i)) + ATOL(i),
C          and the RMS-norm (root-mean-square norm) here is
C          RMS-norm(v) = SQRT(sum v(i)**2 / NEQ).  Here EWT = (EWT(i))
C          is a vector of weights which must always be positive, and
C          the values of RTOL and ATOL should all be non-negative.
C          The following table gives the types (scalar/array) of
C          RTOL and ATOL, and the corresponding form of EWT(i).
C
C             ITOL    RTOL       ATOL          EWT(i)
C              1     scalar     scalar     RTOL*ABS(Y(i)) + ATOL
C              2     scalar     array      RTOL*ABS(Y(i)) + ATOL(i)
C              3     array      scalar     RTOL(i)*ABS(Y(i)) + ATOL
C              4     array      scalar     RTOL(i)*ABS(Y(i)) + ATOL(i)
C
C          When either of these parameters is a scalar, it need not
C          be dimensioned in the user's calling program.
C
C          If none of the above choices (with ITOL, RTOL, and ATOL
C          fixed throughout the problem) is suitable, more general
C          error controls can be obtained by substituting
C          user-supplied routines for the setting of EWT and/or for
C          the norm calculation.  See Part 4 below.
C
C          If global errors are to be estimated by making a repeated
C          run on the same problem with smaller tolerances, then all
C          components of RTOL and ATOL (i.e. of EWT) should be scaled
C          down uniformly.
C
C ITASK  = an index specifying the task to be performed.
C          Input only.  ITASK has the following values and meanings.
C          1  means normal computation of output values of y(t) at
C             t = TOUT (by overshooting and interpolating).
C          2  means take one step only and return.
C          3  means stop at the first internal mesh point at or
C             beyond t = TOUT and return.
C          4  means normal computation of output values of y(t) at
C             t = TOUT but without overshooting t = TCRIT.
C             TCRIT must be input as RWORK(1).  TCRIT may be equal to
C             or beyond TOUT, but not behind it in the direction of
C             integration.  This option is useful if the problem
C             has a singularity at or beyond t = TCRIT.
C          5  means take one step, without passing TCRIT, and return.
C             TCRIT must be input as RWORK(1).
C
C          Note:  If ITASK = 4 or 5 and the solver reaches TCRIT
C          (within roundoff), it will return T = TCRIT (exactly) to
C          indicate this (unless ITASK = 4 and TOUT comes before TCRIT,
C          in which case answers at t = TOUT are returned first).
C
C ISTATE = an index used for input and output to specify the
C          state of the calculation.
C
C          On input, the values of ISTATE are as follows.
C          0  means this is the first call for the problem, and
C             DLSODI is to compute the initial value of dy/dt
C             (while doing other initializations).  See note below.
C          1  means this is the first call for the problem, and
C             the initial value of dy/dt has been supplied in
C             YDOTI (DLSODI will do other initializations).  See note
C             below.
C          2  means this is not the first call, and the calculation
C             is to continue normally, with no change in any input
C             parameters except possibly TOUT and ITASK.
C             (If ITOL, RTOL, and/or ATOL are changed between calls
C             with ISTATE = 2, the new values will be used but not
C             tested for legality.)
C          3  means this is not the first call, and the
C             calculation is to continue normally, but with
C             a change in input parameters other than
C             TOUT and ITASK.  Changes are allowed in
C             NEQ, ITOL, RTOL, ATOL, IOPT, LRW, LIW, MF, ML, MU,
C             and any of the optional inputs except H0.
C             (See IWORK description for ML and MU.)
C          Note:  A preliminary call with TOUT = T is not counted
C          as a first call here, as no initialization or checking of
C          input is done.  (Such a call is sometimes useful for the
C          purpose of outputting the initial conditions.)
C          Thus the first call for which TOUT .ne. T requires
C          ISTATE = 0 or 1 on input.
C
C          On output, ISTATE has the following values and meanings.
C           0 or 1  means nothing was done; TOUT = t and
C              ISTATE = 0 or 1 on input.
C           2  means that the integration was performed successfully.
C           3  means that the user-supplied Subroutine RES signalled
C              DLSODI to halt the integration and return (IRES = 2).
C              Integration as far as T was achieved with no occurrence
C              of IRES = 2, but this flag was set on attempting the
C              next step.
C          -1  means an excessive amount of work (more than MXSTEP
C              steps) was done on this call, before completing the
C              requested task, but the integration was otherwise
C              successful as far as T.  (MXSTEP is an optional input
C              and is normally 500.)  To continue, the user may
C              simply reset ISTATE to a value .gt. 1 and call again
C              (the excess work step counter will be reset to 0).
C              In addition, the user may increase MXSTEP to avoid
C              this error return (see below on optional inputs).
C          -2  means too much accuracy was requested for the precision
C              of the machine being used.  This was detected before
C              completing the requested task, but the integration
C              was successful as far as T.  To continue, the tolerance
C              parameters must be reset, and ISTATE must be set
C              to 3.  The optional output TOLSF may be used for this
C              purpose.  (Note: If this condition is detected before
C              taking any steps, then an illegal input return
C              (ISTATE = -3) occurs instead.)
C          -3  means illegal input was detected, before taking any
C              integration steps.  See written message for details.
C              Note:  If the solver detects an infinite loop of calls
C              to the solver with illegal input, it will cause
C              the run to stop.
C          -4  means there were repeated error test failures on
C              one attempted step, before completing the requested
C              task, but the integration was successful as far as T.
C              The problem may have a singularity, or the input
C              may be inappropriate.
C          -5  means there were repeated convergence test failures on
C              one attempted step, before completing the requested
C              task, but the integration was successful as far as T.
C              This may be caused by an inaccurate Jacobian matrix.
C          -6  means EWT(i) became zero for some i during the
C              integration.  pure relative error control (ATOL(i)=0.0)
C              was requested on a variable which has now vanished.
C              the integration was successful as far as T.
C          -7  means that the user-supplied Subroutine RES set
C              its error flag (IRES = 3) despite repeated tries by
C              DLSODI to avoid that condition.
C          -8  means that ISTATE was 0 on input but DLSODI was unable
C              to compute the initial value of dy/dt.  See the
C              printed message for details.
C
C          Note:  Since the normal output value of ISTATE is 2,
C          it does not need to be reset for normal continuation.
C          Similarly, ISTATE (= 3) need not be reset if RES told
C          DLSODI to return because the calling program must change
C          the parameters of the problem.
C          Also, since a negative input value of ISTATE will be
C          regarded as illegal, a negative output value requires the
C          user to change it, and possibly other inputs, before
C          calling the solver again.
C
C IOPT   = an integer flag to specify whether or not any optional
C          inputs are being used on this call.  Input only.
C          The optional inputs are listed separately below.
C          IOPT = 0 means no optional inputs are being used.
C                   Default values will be used in all cases.
C          IOPT = 1 means one or more optional inputs are being used.
C
C RWORK  = a real working array (double precision).
C          The length of RWORK must be at least
C             20 + NYH*(MAXORD + 1) + 3*NEQ + LENWM    where
C          NYH    = the initial value of NEQ,
C          MAXORD = 12 (if METH = 1) or 5 (if METH = 2) (unless a
C                   smaller value is given as an optional input),
C          LENWM   = NEQ**2 + 2    if MITER is 1 or 2, and
C          LENWM   = (2*ML+MU+1)*NEQ + 2 if MITER is 4 or 5.
C          (See MF description for the definition of METH and MITER.)
C          Thus if MAXORD has its default value and NEQ is constant,
C          this length is
C             22 + 16*NEQ + NEQ**2         for MF = 11 or 12,
C             22 + 17*NEQ + (2*ML+MU)*NEQ  for MF = 14 or 15,
C             22 +  9*NEQ + NEQ**2         for MF = 21 or 22,
C             22 + 10*NEQ + (2*ML+MU)*NEQ  for MF = 24 or 25.
C          The first 20 words of RWORK are reserved for conditional
C          and optional inputs and optional outputs.
C
C          The following word in RWORK is a conditional input:
C            RWORK(1) = TCRIT = critical value of t which the solver
C                       is not to overshoot.  Required if ITASK is
C                       4 or 5, and ignored otherwise.  (See ITASK.)
C
C LRW    = the length of the array RWORK, as declared by the user.
C          (This will be checked by the solver.)
C
C IWORK  = an integer work array.  The length of IWORK must be at least
C          20 + NEQ .  The first few words of IWORK are used for
C          conditional and optional inputs and optional outputs.
C
C          The following 2 words in IWORK are conditional inputs:
C            IWORK(1) = ML     These are the lower and upper
C            IWORK(2) = MU     half-bandwidths, respectively, of the
C                       matrices in the problem-- the Jacobian dr/dy
C                       and the left-hand side matrix A. These
C                       half-bandwidths exclude the main diagonal,
C                       so the total bandwidth is ML + MU + 1 .
C                       The band is defined by the matrix locations
C                       (i,j) with i-ML .le. j .le. i+MU.  ML and MU
C                       must satisfy  0 .le.  ML,MU  .le. NEQ-1.
C                       These are required if MITER is 4 or 5, and
C                       ignored otherwise.
C                       ML and MU may in fact be the band parameters
C                       for matrices to which dr/dy and A are only
C                       approximately equal.
C
C LIW    = the length of the array IWORK, as declared by the user.
C          (This will be checked by the solver.)
C
C Note:  The work arrays must not be altered between calls to DLSODI
C for the same problem, except possibly for the conditional and
C optional inputs, and except for the last 3*NEQ words of RWORK.
C The latter space is used for internal scratch space, and so is
C available for use by the user outside DLSODI between calls, if
C desired (but not for use by RES, ADDA, or JAC).
C
C MF     = the method flag.  Used only for input.  The legal values of
C          MF are 11, 12, 14, 15, 21, 22, 24, and 25.
C          MF has decimal digits METH and MITER: MF = 10*METH + MITER.
C            METH indicates the basic linear multistep method:
C              METH = 1 means the implicit Adams method.
C              METH = 2 means the method based on Backward
C                       Differentiation Formulas (BDFs).
C                The BDF method is strongly preferred for stiff
C              problems, while the Adams method is preferred when
C              the problem is not stiff.  If the matrix A(t,y) is
C              nonsingular, stiffness here can be taken to mean that of
C              the explicit ODE system dy/dt = A-inverse * g.  If A is
C              singular, the concept of stiffness is not well defined.
C                If you do not know whether the problem is stiff, we
C              recommend using METH = 2.  If it is stiff, the advantage
C              of METH = 2 over METH = 1 will be great, while if it is
C              not stiff, the advantage of METH = 1 will be slight.
C              If maximum efficiency is important, some experimentation
C              with METH may be necessary.
C            MITER indicates the corrector iteration method:
C              MITER = 1 means chord iteration with a user-supplied
C                        full (NEQ by NEQ) Jacobian.
C              MITER = 2 means chord iteration with an internally
C                        generated (difference quotient) full Jacobian.
C                        This uses NEQ+1 extra calls to RES per dr/dy
C                        evaluation.
C              MITER = 4 means chord iteration with a user-supplied
C                        banded Jacobian.
C              MITER = 5 means chord iteration with an internally
C                        generated banded Jacobian (using ML+MU+2
C                        extra calls to RES per dr/dy evaluation).
C              If MITER = 1 or 4, the user must supply a Subroutine JAC
C              (the name is arbitrary) as described above under JAC.
C              For other values of MITER, a dummy argument can be used.
C-----------------------------------------------------------------------
C Optional Inputs.
C
C The following is a list of the optional inputs provided for in the
C call sequence.  (See also Part 2.)  For each such input variable,
C this table lists its name as used in this documentation, its
C location in the call sequence, its meaning, and the default value.
C the use of any of these inputs requires IOPT = 1, and in that
C case all of these inputs are examined.  A value of zero for any
C of these optional inputs will cause the default value to be used.
C Thus to use a subset of the optional inputs, simply preload
C locations 5 to 10 in RWORK and IWORK to 0.0 and 0 respectively, and
C then set those of interest to nonzero values.
C
C Name    Location      Meaning and Default Value
C
C H0      RWORK(5)  the step size to be attempted on the first step.
C                   The default value is determined by the solver.
C
C HMAX    RWORK(6)  the maximum absolute step size allowed.
C                   The default value is infinite.
C
C HMIN    RWORK(7)  the minimum absolute step size allowed.
C                   The default value is 0.  (This lower bound is not
C                   enforced on the final step before reaching TCRIT
C                   when ITASK = 4 or 5.)
C
C MAXORD  IWORK(5)  the maximum order to be allowed.  The default
C                   value is 12 if METH = 1, and 5 if METH = 2.
C                   If MAXORD exceeds the default value, it will
C                   be reduced to the default value.
C                   If MAXORD is changed during the problem, it may
C                   cause the current order to be reduced.
C
C MXSTEP  IWORK(6)  maximum number of (internally defined) steps
C                   allowed during one call to the solver.
C                   The default value is 500.
C
C MXHNIL  IWORK(7)  maximum number of messages printed (per problem)
C                   warning that T + H = T on a step (H = step size).
C                   This must be positive to result in a non-default
C                   value.  The default value is 10.
C-----------------------------------------------------------------------
C Optional Outputs.
C
C As optional additional output from DLSODI, the variables listed
C below are quantities related to the performance of DLSODI
C which are available to the user.  These are communicated by way of
C the work arrays, but also have internal mnemonic names as shown.
C Except where stated otherwise, all of these outputs are defined
C on any successful return from DLSODI, and on any return with
C ISTATE = -1, -2, -4, -5, -6, or -7.  On a return with -3 (illegal
C input) or -8, they will be unchanged from their existing values
C (if any), except possibly for TOLSF, LENRW, and LENIW.
C On any error return, outputs relevant to the error will be defined,
C as noted below.
C
C Name    Location      Meaning
C
C HU      RWORK(11) the step size in t last used (successfully).
C
C HCUR    RWORK(12) the step size to be attempted on the next step.
C
C TCUR    RWORK(13) the current value of the independent variable
C                   which the solver has actually reached, i.e. the
C                   current internal mesh point in t.  On output, TCUR
C                   will always be at least as far as the argument
C                   T, but may be farther (if interpolation was done).
C
C TOLSF   RWORK(14) a tolerance scale factor, greater than 1.0,
C                   computed when a request for too much accuracy was
C                   detected (ISTATE = -3 if detected at the start of
C                   the problem, ISTATE = -2 otherwise).  If ITOL is
C                   left unaltered but RTOL and ATOL are uniformly
C                   scaled up by a factor of TOLSF for the next call,
C                   then the solver is deemed likely to succeed.
C                   (The user may also ignore TOLSF and alter the
C                   tolerance parameters in any other way appropriate.)
C
C NST     IWORK(11) the number of steps taken for the problem so far.
C
C NRE     IWORK(12) the number of residual evaluations (RES calls)
C                   for the problem so far.
C
C NJE     IWORK(13) the number of Jacobian evaluations (each involving
C                   an evaluation of A and dr/dy) for the problem so
C                   far.  This equals the number of calls to ADDA and
C                   (if MITER = 1 or 4) JAC, and the number of matrix
C                   LU decompositions.
C
C NQU     IWORK(14) the method order last used (successfully).
C
C NQCUR   IWORK(15) the order to be attempted on the next step.
C
C IMXER   IWORK(16) the index of the component of largest magnitude in
C                   the weighted local error vector ( E(i)/EWT(i) ),
C                   on an error return with ISTATE = -4 or -5.
C
C LENRW   IWORK(17) the length of RWORK actually required.
C                   This is defined on normal returns and on an illegal
C                   input return for insufficient storage.
C
C LENIW   IWORK(18) the length of IWORK actually required.
C                   This is defined on normal returns and on an illegal
C                   input return for insufficient storage.
C
C
C The following two arrays are segments of the RWORK array which
C may also be of interest to the user as optional outputs.
C For each array, the table below gives its internal name,
C its base address in RWORK, and its description.
C
C Name    Base Address      Description
C
C YH      21             the Nordsieck history array, of size NYH by
C                        (NQCUR + 1), where NYH is the initial value
C                        of NEQ.  For j = 0,1,...,NQCUR, column j+1
C                        of YH contains HCUR**j/factorial(j) times
C                        the j-th derivative of the interpolating
C                        polynomial currently representing the solution,
C                        evaluated at t = TCUR.
C
C ACOR     LENRW-NEQ+1   array of size NEQ used for the accumulated
C                        corrections on each step, scaled on output to
C                        represent the estimated local error in y on the
C                        last step. This is the vector E in the descrip-
C                        tion of the error control.  It is defined only
C                        on a return from DLSODI with ISTATE = 2.
C
C-----------------------------------------------------------------------
C Part 2.  Other Routines Callable.
C
C The following are optional calls which the user may make to
C gain additional capabilities in conjunction with DLSODI.
C (The routines XSETUN and XSETF are designed to conform to the
C SLATEC error handling package.)
C
C     Form of Call                  Function
C   CALL XSETUN(LUN)          Set the logical unit number, LUN, for
C                             output of messages from DLSODI, if
C                             the default is not desired.
C                             The default value of LUN is 6.
C
C   CALL XSETF(MFLAG)         Set a flag to control the printing of
C                             messages by DLSODI.
C                             MFLAG = 0 means do not print. (Danger:
C                             This risks losing valuable information.)
C                             MFLAG = 1 means print (the default).
C
C                             Either of the above calls may be made at
C                             any time and will take effect immediately.
C
C   CALL DSRCOM(RSAV,ISAV,JOB) saves and restores the contents of
C                             the internal Common blocks used by
C                             DLSODI (see Part 3 below).
C                             RSAV must be a real array of length 9
C                             or more, and ISAV must be an integer
C                             array of length 25 or more.
C                             JOB=1 means save Common into RSAV/ISAV.
C                             JOB=2 means restore Common from RSAV/ISAV.
C                                DSRCOM is useful if one is
C                             interrupting a run and restarting
C                             later, or alternating between two or
C                             more problems solved with DLSODI.
C
C   CALL DINTDY(,,,,,)        Provide derivatives of y, of various
C        (see below)          orders, at a specified point t, if
C                             desired.  It may be called only after
C                             a successful return from DLSODI.
C
C The detailed instructions for using DINTDY are as follows.
C The form of the call is:
C
C   CALL DINTDY (T, K, RWORK(21), NYH, DKY, IFLAG)
C
C The input parameters are:
C
C T         = value of independent variable where answers are desired
C             (normally the same as the T last returned by DLSODI).
C             For valid results, T must lie between TCUR - HU and TCUR.
C             (See optional outputs for TCUR and HU.)
C K         = integer order of the derivative desired.  K must satisfy
C             0 .le. K .le. NQCUR, where NQCUR is the current order
C             (see optional outputs).  The capability corresponding
C             to K = 0, i.e. computing y(T), is already provided
C             by DLSODI directly.  Since NQCUR .ge. 1, the first
C             derivative dy/dt is always available with DINTDY.
C RWORK(21) = the base address of the history array YH.
C NYH       = column length of YH, equal to the initial value of NEQ.
C
C The output parameters are:
C
C DKY       = a real array of length NEQ containing the computed value
C             of the K-th derivative of y(t).
C IFLAG     = integer flag, returned as 0 if K and T were legal,
C             -1 if K was illegal, and -2 if T was illegal.
C             On an error return, a message is also written.
C-----------------------------------------------------------------------
C Part 3.  Common Blocks.
C
C If DLSODI is to be used in an overlay situation, the user
C must declare, in the primary overlay, the variables in:
C   (1) the call sequence to DLSODI, and
C   (2) the internal Common block
C         /DLS001/  of length  34  (9 double precision words
C                      followed by 25 integer words).
C
C If DLSODI is used on a system in which the contents of internal
C Common blocks are not preserved between calls, the user should
C declare the above Common block in the calling program to insure
C that their contents are preserved.
C
C If the solution of a given problem by DLSODI is to be interrupted
C and then later continued, such as when restarting an interrupted run
C or alternating between two or more problems, the user should save,
C following the return from the last DLSODI call prior to the
C interruption, the contents of the call sequence variables and the
C internal Common blocks, and later restore these values before the
C next DLSODI call for that problem.  To save and restore the Common
C blocks, use Subroutine DSRCOM (see Part 2 above).
C
C-----------------------------------------------------------------------
C Part 4.  Optionally Replaceable Solver Routines.
C
C Below are descriptions of two routines in the DLSODI package which
C relate to the measurement of errors.  Either routine can be
C replaced by a user-supplied version, if desired.  However, since such
C a replacement may have a major impact on performance, it should be
C done only when absolutely necessary, and only with great caution.
C (Note: The means by which the package version of a routine is
C superseded by the user's version may be system-dependent.)
C
C (a) DEWSET.
C The following subroutine is called just before each internal
C integration step, and sets the array of error weights, EWT, as
C described under ITOL/RTOL/ATOL above:
C     SUBROUTINE DEWSET (NEQ, ITOL, RTOL, ATOL, YCUR, EWT)
C where NEQ, ITOL, RTOL, and ATOL are as in the DLSODI call sequence,
C YCUR contains the current dependent variable vector, and
C EWT is the array of weights set by DEWSET.
C
C If the user supplies this subroutine, it must return in EWT(i)
C (i = 1,...,NEQ) a positive quantity suitable for comparing errors
C in y(i) to.  The EWT array returned by DEWSET is passed to the DVNORM
C routine (see below), and also used by DLSODI in the computation
C of the optional output IMXER, the diagonal Jacobian approximation,
C and the increments for difference quotient Jacobians.
C
C In the user-supplied version of DEWSET, it may be desirable to use
C the current values of derivatives of y.  Derivatives up to order NQ
C are available from the history array YH, described above under
C optional outputs.  In DEWSET, YH is identical to the YCUR array,
C extended to NQ + 1 columns with a column length of NYH and scale
C factors of H**j/factorial(j).  On the first call for the problem,
C given by NST = 0, NQ is 1 and H is temporarily set to 1.0.
C NYH is the initial value of NEQ.  The quantities NQ, H, and NST
C can be obtained by including in DEWSET the statements:
C     DOUBLE PRECISION RLS
C     COMMON /DLS001/ RLS(9),ILS(25)
C     NQ = ILS(21)
C     NST = ILS(22)
C     H = RLS(3)
C Thus, for example, the current value of dy/dt can be obtained as
C YCUR(NYH+i)/H  (i=1,...,NEQ)  (and the division by H is
C unnecessary when NST = 0).
C
C (b) DVNORM.
C The following is a real function routine which computes the weighted
C root-mean-square norm of a vector v:
C     D = DVNORM (N, V, W)
C where:
C   N = the length of the vector,
C   V = real array of length N containing the vector,
C   W = real array of length N containing weights,
C   D = SQRT( (1/N) * sum(V(i)*W(i))**2 ).
C DVNORM is called with N = NEQ and with W(i) = 1.0/EWT(i), where
C EWT is as set by Subroutine DEWSET.
C
C If the user supplies this function, it should return a non-negative
C value of DVNORM suitable for use in the error control in DLSODI.
C None of the arguments should be altered by DVNORM.
C For example, a user-supplied DVNORM routine might:
C   -substitute a max-norm of (V(i)*W(i)) for the RMS-norm, or
C   -ignore some components of V in the norm, with the effect of
C    suppressing the error control on those components of y.
C-----------------------------------------------------------------------
C
C***REVISION HISTORY  (YYYYMMDD)
C 19800424  DATE WRITTEN
C 19800519  Corrected access of YH on forced order reduction;
C           numerous corrections to prologues and other comments.
C 19800617  In main driver, added loading of SQRT(UROUND) in RWORK;
C           minor corrections to main prologue.
C 19800903  Corrected ISTATE logic; minor changes in prologue.
C 19800923  Added zero initialization of HU and NQU.
C 19801028  Reorganized RES calls in AINVG, STODI, and PREPJI;
C           in LSODI, corrected NRE increment and reset LDY0 at 580;
C           numerous corrections to main prologue.
C 19801218  Revised XERRWD routine; minor corrections to main prologue.
C 19810330  Added Common block /LSI001/; use LSODE's INTDY and SOLSY;
C           minor corrections to XERRWD and error message at 604;
C           minor corrections to declarations; corrections to prologues.
C 19810818  Numerous revisions: replaced EWT by 1/EWT; used flags
C           JCUR, ICF, IERPJ, IERSL between STODI and subordinates;
C           added tuning parameters CCMAX, MAXCOR, MSBP, MXNCF;
C           reorganized returns from STODI; reorganized type decls.;
C           fixed message length in XERRWD; changed default LUNIT to 6;
C           changed Common lengths; changed comments throughout.
C 19820906  Corrected use of ABS(H) in STODI; minor comment fixes.
C 19830510  Numerous revisions: revised diff. quotient increment;
C           eliminated block /LSI001/, using IERPJ flag;
C           revised STODI logic after PJAC return;
C           revised tuning of H change and step attempts in STODI;
C           corrections to main prologue and internal comments.
C 19870330  Major update: corrected comments throughout;
C           removed TRET from Common; rewrote EWSET with 4 loops;
C           fixed t test in INTDY; added Cray directives in STODI;
C           in STODI, fixed DELP init. and logic around PJAC call;
C           combined routines to save/restore Common;
C           passed LEVEL = 0 in error message calls (except run abort).
C 20010425  Major update: convert source lines to upper case;
C           added *DECK lines; changed from 1 to * in dummy dimensions;
C           changed names R1MACH/D1MACH to RUMACH/DUMACH;
C           renamed routines for uniqueness across single/double prec.;
C           converted intrinsic names to generic form;
C           removed ILLIN and NTREP (data loaded) from Common;
C           removed all 'own' variables from Common;
C           changed error messages to quoted strings;
C           replaced XERRWV/XERRWD with 1993 revised version;
C           converted prologues, comments, error messages to mixed case;
C           converted arithmetic IF statements to logical IF statements;
C           numerous corrections to prologues and internal comments.
C 20010507  Converted single precision source to double precision.
C 20020502  Corrected declarations in descriptions of user routines.
C
C-----------------------------------------------------------------------
C Other routines in the DLSODI package.
C
C In addition to Subroutine DLSODI, the DLSODI package includes the
C following subroutines and function routines:
C  DAINVG   computes the initial value of the vector
C             dy/dt = A-inverse * g
C  DINTDY   computes an interpolated value of the y vector at t = TOUT.
C  DSTODI   is the core integrator, which does one step of the
C           integration and the associated error control.
C  DCFODE   sets all method coefficients and test constants.
C  DPREPJI  computes and preprocesses the Jacobian matrix
C           and the Newton iteration matrix P.
C  DSOLSY   manages solution of linear system in chord iteration.
C  DEWSET   sets the error weight vector EWT before each step.
C  DVNORM   computes the weighted RMS-norm of a vector.
C  DSRCOM   is a user-callable routine to save and restore
C           the contents of the internal Common blocks.
C  DGEFA and DGESL   are routines from LINPACK for solving full
C           systems of linear algebraic equations.
C  DGBFA and DGBSL   are routines from LINPACK for solving banded
C           linear systems.
C  DUMACH   computes the unit roundoff in a machine-independent manner.
C  XERRWD, XSETUN, XSETF, IXSAV, and IUMACH  handle the printing of all
C           error messages and warnings.  XERRWD is machine-dependent.
C Note:  DVNORM, DUMACH, IXSAV, and IUMACH are function routines.
C All the others are subroutines.
C
C-----------------------------------------------------------------------
      EXTERNAL DPREPJI, DSOLSY
      DOUBLE PRECISION DUMACH, DVNORM
      INTEGER ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, 
     1   LYH, LEWT, LACOR, LSAVR, LWM, LIWM, METH, MITER,
     2   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NRE, NJE, NQU
      INTEGER I, I1, I2, IER, IFLAG, IMXER, IRES, KGO,
     1   LENIW, LENRW, LENWM, LP, LYD0, ML, MORD, MU, MXHNL0, MXSTP0
      INTEGER INIT, MXSTEP, MXHNIL, NHNIL, NSLAST, NYH
      DOUBLE PRECISION CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
      DOUBLE PRECISION ATOLI, AYI, BIG, EWTI, H0, HMAX, HMX, RH, RTOLI,
     1   TCRIT, TDIST, TNEXT, TOL, TOLSF, TP, SIZE, SUM, W0
      DIMENSION MORD(2)
      LOGICAL IHIT
      CHARACTER*60 MSG
      SAVE INIT, MXSTEP, MXHNIL, NHNIL, NSLAST, NYH
C-----------------------------------------------------------------------
C The following internal Common block contains variables which are
C communicated between subroutines.  All real variables are listed
C first, followed by all integers.  The block is declared in
C Subroutines DLSODI, DINTDY, DSTODI, DPREPJI, DSOLSY.
C-----------------------------------------------------------------------
      COMMON /DLS001/ CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,
     1   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L, 
     2   LYH, LEWT, LACOR, LSAVR, LWM, LIWM, METH, MITER,
     3   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NRE, NJE, NQU
C
      DATA MORD(1),MORD(2)/12,5/, MXSTP0/500/, MXHNL0/10/
C-----------------------------------------------------------------------
C Block A.
C This code block is executed on every call.
C It tests ISTATE and ITASK for legality and branches appropriately.
C If ISTATE .gt. 1 but the flag INIT shows that initialization has
C not yet been done, an error return occurs.
C If ISTATE = 0 or 1 and TOUT = T, return immediately.
C-----------------------------------------------------------------------
      IF (ISTATE .LT. 0 .OR. ISTATE .GT. 3) GO TO 601
      IF (ITASK .LT. 1 .OR. ITASK .GT. 5) GO TO 602
      IF (ISTATE .LE. 1) GO TO 10
      IF (INIT .EQ. 0) GO TO 603
      IF (ISTATE .EQ. 2) GO TO 200
      GO TO 20
 10   INIT = 0
      IF (TOUT .EQ. T) RETURN
C-----------------------------------------------------------------------
C Block B.
C The next code block is executed for the initial call (ISTATE = 0 or 1)
C or for a continuation call with parameter changes (ISTATE = 3).
C It contains checking of all inputs and various initializations.
C
C First check legality of the non-optional inputs NEQ, ITOL, IOPT,
C MF, ML, and MU.
C-----------------------------------------------------------------------
 20   IF (NEQ(1) .LE. 0) GO TO 604
      IF (ISTATE .LE. 1) GO TO 25
      IF (NEQ(1) .GT. N) GO TO 605
 25   N = NEQ(1)
      IF (ITOL .LT. 1 .OR. ITOL .GT. 4) GO TO 606
      IF (IOPT .LT. 0 .OR. IOPT .GT. 1) GO TO 607
      METH = MF/10
      MITER = MF - 10*METH
      IF (METH .LT. 1 .OR. METH .GT. 2) GO TO 608
      IF (MITER .LE. 0 .OR. MITER .GT. 5) GO TO 608
      IF (MITER .EQ. 3)  GO TO 608
      IF (MITER .LT. 3) GO TO 30
      ML = IWORK(1)
      MU = IWORK(2)
      IF (ML .LT. 0 .OR. ML .GE. N) GO TO 609
      IF (MU .LT. 0 .OR. MU .GE. N) GO TO 610
 30   CONTINUE
C Next process and check the optional inputs. --------------------------
      IF (IOPT .EQ. 1) GO TO 40
      MAXORD = MORD(METH)
      MXSTEP = MXSTP0
      MXHNIL = MXHNL0
      IF (ISTATE .LE. 1) H0 = 0.0D0
      HMXI = 0.0D0
      HMIN = 0.0D0
      GO TO 60
 40   MAXORD = IWORK(5)
      IF (MAXORD .LT. 0) GO TO 611
      IF (MAXORD .EQ. 0) MAXORD = 100
      MAXORD = MIN(MAXORD,MORD(METH))
      MXSTEP = IWORK(6)
      IF (MXSTEP .LT. 0) GO TO 612
      IF (MXSTEP .EQ. 0) MXSTEP = MXSTP0
      MXHNIL = IWORK(7)
      IF (MXHNIL .LT. 0) GO TO 613
      IF (MXHNIL .EQ. 0) MXHNIL = MXHNL0
      IF (ISTATE .GT. 1) GO TO 50
      H0 = RWORK(5)
      IF ((TOUT - T)*H0 .LT. 0.0D0) GO TO 614
 50   HMAX = RWORK(6)
      IF (HMAX .LT. 0.0D0) GO TO 615
      HMXI = 0.0D0
      IF (HMAX .GT. 0.0D0) HMXI = 1.0D0/HMAX
      HMIN = RWORK(7)
      IF (HMIN .LT. 0.0D0) GO TO 616
C-----------------------------------------------------------------------
C Set work array pointers and check lengths LRW and LIW.
C Pointers to segments of RWORK and IWORK are named by prefixing L to
C the name of the segment.  E.g., the segment YH starts at RWORK(LYH).
C Segments of RWORK (in order) are denoted YH, WM, EWT, SAVR, ACOR.
C-----------------------------------------------------------------------
 60   LYH = 21
      IF (ISTATE .LE. 1) NYH = N
      LWM = LYH + (MAXORD + 1)*NYH
      IF (MITER .LE. 2) LENWM = N*N + 2
      IF (MITER .GE. 4) LENWM = (2*ML + MU + 1)*N + 2
      LEWT = LWM + LENWM
      LSAVR = LEWT + N
      LACOR = LSAVR + N
      LENRW = LACOR + N - 1
      IWORK(17) = LENRW
      LIWM = 1
      LENIW = 20 + N
      IWORK(18) = LENIW
      IF (LENRW .GT. LRW) GO TO 617
      IF (LENIW .GT. LIW) GO TO 618
C Check RTOL and ATOL for legality. ------------------------------------
      RTOLI = RTOL(1)
      ATOLI = ATOL(1)
      DO 70 I = 1,N
        IF (ITOL .GE. 3) RTOLI = RTOL(I)
        IF (ITOL .EQ. 2 .OR. ITOL .EQ. 4) ATOLI = ATOL(I)
        IF (RTOLI .LT. 0.0D0) GO TO 619
        IF (ATOLI .LT. 0.0D0) GO TO 620
 70     CONTINUE
      IF (ISTATE .LE. 1) GO TO 100
C If ISTATE = 3, set flag to signal parameter changes to DSTODI. -------
      JSTART = -1
      IF (NQ .LE. MAXORD) GO TO 90
C MAXORD was reduced below NQ.  Copy YH(*,MAXORD+2) into YDOTI.---------
      DO 80 I = 1,N
 80     YDOTI(I) = RWORK(I+LWM-1)
C Reload WM(1) = RWORK(lWM), since lWM may have changed. ---------------
 90   RWORK(LWM) = SQRT(UROUND)
      IF (N .EQ. NYH) GO TO 200
C NEQ was reduced.  Zero part of YH to avoid undefined references. -----
      I1 = LYH + L*NYH
      I2 = LYH + (MAXORD + 1)*NYH - 1
      IF (I1 .GT. I2) GO TO 200
      DO 95 I = I1,I2
 95     RWORK(I) = 0.0D0
      GO TO 200
C-----------------------------------------------------------------------
C Block C.
C The next block is for the initial call only (ISTATE = 0 or 1).
C It contains all remaining initializations, the call to DAINVG
C (if ISTATE = 1), and the calculation of the initial step size.
C The error weights in EWT are inverted after being loaded.
C-----------------------------------------------------------------------
 100  UROUND = DUMACH()
      TN = T
      IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GO TO 105
      TCRIT = RWORK(1)
      IF ((TCRIT - TOUT)*(TOUT - T) .LT. 0.0D0) GO TO 625
      IF (H0 .NE. 0.0D0 .AND. (T + H0 - TCRIT)*H0 .GT. 0.0D0)
     1   H0 = TCRIT - T
 105  JSTART = 0
      RWORK(LWM) = SQRT(UROUND)
      NHNIL = 0
      NST = 0
      NRE = 0
      NJE = 0
      NSLAST = 0
      HU = 0.0D0
      NQU = 0
      CCMAX = 0.3D0
      MAXCOR = 3
      MSBP = 20
      MXNCF = 10
C Compute initial dy/dt, if necessary, and load it and initial Y into YH
      LYD0 = LYH + NYH
      LP = LWM + 1
      IF (ISTATE .EQ. 1) GO TO 120
C DLSODI must compute initial dy/dt (LYD0 points to YH(*,2)). ----------
         CALL DAINVG( RES, ADDA, NEQ, T, Y, RWORK(LYD0), MITER,
     1                ML, MU, RWORK(LP), IWORK(21), IER )
         NRE = NRE + 1
         IF (IER .LT. 0) GO TO 560
         IF (IER .GT. 0) GO TO 565
         DO 115 I = 1,N
 115        RWORK(I+LYH-1) = Y(I)
         GO TO 130
C Initial dy/dt was supplied.  Load into YH (LYD0 points to YH(*,2).). -
 120     DO 125 I = 1,N
            RWORK(I+LYH-1) = Y(I)
 125        RWORK(I+LYD0-1) = YDOTI(I)
C Load and invert the EWT array.  (H is temporarily set to 1.0.) -------
 130  CONTINUE
      NQ = 1
      H = 1.0D0
      CALL DEWSET (N, ITOL, RTOL, ATOL, RWORK(LYH), RWORK(LEWT))
      DO 135 I = 1,N
        IF (RWORK(I+LEWT-1) .LE. 0.0D0) GO TO 621
 135    RWORK(I+LEWT-1) = 1.0D0/RWORK(I+LEWT-1)
C-----------------------------------------------------------------------
C The coding below computes the step size, H0, to be attempted on the
C first step, unless the user has supplied a value for this.
C First check that TOUT - T differs significantly from zero.
C A scalar tolerance quantity TOL is computed, as MAX(RTOL(i))
C if this is positive, or MAX(ATOL(i)/ABS(Y(i))) otherwise, adjusted
C so as to be between 100*UROUND and 1.0E-3.
C Then the computed value H0 is given by..
C                                      NEQ
C   H0**2 = TOL / ( w0**-2 + (1/NEQ) * Sum ( YDOT(i)/ywt(i) )**2  )
C                                       1
C where   w0      = MAX ( ABS(T), ABS(TOUT) ),
C         YDOT(i) = i-th component of initial value of dy/dt,
C         ywt(i)  = EWT(i)/TOL  (a weight for y(i)).
C The sign of H0 is inferred from the initial values of TOUT and T.
C-----------------------------------------------------------------------
      IF (H0 .NE. 0.0D0) GO TO 180
      TDIST = ABS(TOUT - T)
      W0 = MAX(ABS(T),ABS(TOUT))
      IF (TDIST .LT. 2.0D0*UROUND*W0) GO TO 622
      TOL = RTOL(1)
      IF (ITOL .LE. 2) GO TO 145
      DO 140 I = 1,N
 140    TOL = MAX(TOL,RTOL(I))
 145  IF (TOL .GT. 0.0D0) GO TO 160
      ATOLI = ATOL(1)
      DO 150 I = 1,N
        IF (ITOL .EQ. 2 .OR. ITOL .EQ. 4) ATOLI = ATOL(I)
        AYI = ABS(Y(I))
        IF (AYI .NE. 0.0D0) TOL = MAX(TOL,ATOLI/AYI)
 150    CONTINUE
 160  TOL = MAX(TOL,100.0D0*UROUND)
      TOL = MIN(TOL,0.001D0)
      SUM = DVNORM (N, RWORK(LYD0), RWORK(LEWT))
      SUM = 1.0D0/(TOL*W0*W0) + TOL*SUM**2
      H0 = 1.0D0/SQRT(SUM)
      H0 = MIN(H0,TDIST)
      H0 = SIGN(H0,TOUT-T)
C Adjust H0 if necessary to meet HMAX bound. ---------------------------
 180  RH = ABS(H0)*HMXI
      IF (RH .GT. 1.0D0) H0 = H0/RH
C Load H with H0 and scale YH(*,2) by H0. ------------------------------
      H = H0
      DO 190 I = 1,N
 190    RWORK(I+LYD0-1) = H0*RWORK(I+LYD0-1)
      GO TO 270
C-----------------------------------------------------------------------
C Block D.
C The next code block is for continuation calls only (ISTATE = 2 or 3)
C and is to check stop conditions before taking a step.
C-----------------------------------------------------------------------
 200  NSLAST = NST
      GO TO (210, 250, 220, 230, 240), ITASK
 210  IF ((TN - TOUT)*H .LT. 0.0D0) GO TO 250
      CALL DINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      IF (IFLAG .NE. 0) GO TO 627
      T = TOUT
      GO TO 420
 220  TP = TN - HU*(1.0D0 + 100.0D0*UROUND)
      IF ((TP - TOUT)*H .GT. 0.0D0) GO TO 623
      IF ((TN - TOUT)*H .LT. 0.0D0) GO TO 250
      GO TO 400
 230  TCRIT = RWORK(1)
      IF ((TN - TCRIT)*H .GT. 0.0D0) GO TO 624
      IF ((TCRIT - TOUT)*H .LT. 0.0D0) GO TO 625
      IF ((TN - TOUT)*H .LT. 0.0D0) GO TO 245
      CALL DINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      IF (IFLAG .NE. 0) GO TO 627
      T = TOUT
      GO TO 420
 240  TCRIT = RWORK(1)
      IF ((TN - TCRIT)*H .GT. 0.0D0) GO TO 624
 245  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. 100.0D0*UROUND*HMX
      IF (IHIT) GO TO 400
      TNEXT = TN + H*(1.0D0 + 4.0D0*UROUND)
      IF ((TNEXT - TCRIT)*H .LE. 0.0D0) GO TO 250
      H = (TCRIT - TN)*(1.0D0 - 4.0D0*UROUND)
      IF (ISTATE .EQ. 2) JSTART = -2
C-----------------------------------------------------------------------
C Block E.
C The next block is normally executed for all calls and contains
C the call to the one-step core integrator DSTODI.
C
C This is a looping point for the integration steps.
C
C First check for too many steps being taken, update EWT (if not at
C start of problem), check for too much accuracy being requested, and
C check for H below the roundoff level in T.
C-----------------------------------------------------------------------
 250  CONTINUE
      IF ((NST-NSLAST) .GE. MXSTEP) GO TO 500
      CALL DEWSET (N, ITOL, RTOL, ATOL, RWORK(LYH), RWORK(LEWT))
      DO 260 I = 1,N
        IF (RWORK(I+LEWT-1) .LE. 0.0D0) GO TO 510
 260    RWORK(I+LEWT-1) = 1.0D0/RWORK(I+LEWT-1)
 270  TOLSF = UROUND*DVNORM (N, RWORK(LYH), RWORK(LEWT))
      IF (TOLSF .LE. 1.0D0) GO TO 280
      TOLSF = TOLSF*2.0D0
      IF (NST .EQ. 0) GO TO 626
      GO TO 520
 280  IF ((TN + H) .NE. TN) GO TO 290
      NHNIL = NHNIL + 1
      IF (NHNIL .GT. MXHNIL) GO TO 290
      MSG = 'DLSODI-  Warning..Internal T (=R1) and H (=R2) are'
      CALL XERRWD (MSG, 50, 101, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG='      such that in the machine, T + H = T on the next step  '
      CALL XERRWD (MSG, 60, 101, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '     (H = step size). Solver will continue anyway.'
      CALL XERRWD (MSG, 50, 101, 0, 0, 0, 0, 2, TN, H)
      IF (NHNIL .LT. MXHNIL) GO TO 290
      MSG = 'DLSODI-  Above warning has been issued I1 times.  '
      CALL XERRWD (MSG, 50, 102, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '     It will not be issued again for this problem.'
      CALL XERRWD (MSG, 50, 102, 0, 1, MXHNIL, 0, 0, 0.0D0, 0.0D0)
 290  CONTINUE
C-----------------------------------------------------------------------
C     CALL DSTODI(NEQ,Y,YH,NYH,YH1,EWT,SAVF,SAVR,ACOR,WM,IWM,RES,
C                 ADDA,JAC,DPREPJI,DSOLSY)
C Note: SAVF in DSTODI occupies the same space as YDOTI in DLSODI.
C-----------------------------------------------------------------------
      CALL DSTODI (NEQ, Y, RWORK(LYH), NYH, RWORK(LYH), RWORK(LEWT),
     1   YDOTI, RWORK(LSAVR), RWORK(LACOR), RWORK(LWM),
     2   IWORK(LIWM), RES, ADDA, JAC, DPREPJI, DSOLSY )
      KGO = 1 - KFLAG
      GO TO (300, 530, 540, 400, 550), KGO
C
C KGO = 1:success; 2:error test failure; 3:convergence failure;
C       4:RES ordered return. 5:RES returned error.
C-----------------------------------------------------------------------
C Block F.
C The following block handles the case of a successful return from the
C core integrator (KFLAG = 0).  Test for stop conditions.
C-----------------------------------------------------------------------
 300  INIT = 1
      GO TO (310, 400, 330, 340, 350), ITASK
C ITASK = 1.  If TOUT has been reached, interpolate. -------------------
 310  IF ((TN - TOUT)*H .LT. 0.0D0) GO TO 250
      CALL DINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      T = TOUT
      GO TO 420
C ITASK = 3.  Jump to exit if TOUT was reached. ------------------------
 330  IF ((TN - TOUT)*H .GE. 0.0D0) GO TO 400
      GO TO 250
C ITASK = 4.  see if TOUT or TCRIT was reached.  adjust h if necessary.
 340  IF ((TN - TOUT)*H .LT. 0.0D0) GO TO 345
      CALL DINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      T = TOUT
      GO TO 420
 345  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. 100.0D0*UROUND*HMX
      IF (IHIT) GO TO 400
      TNEXT = TN + H*(1.0D0 + 4.0D0*UROUND)
      IF ((TNEXT - TCRIT)*H .LE. 0.0D0) GO TO 250
      H = (TCRIT - TN)*(1.0D0 - 4.0D0*UROUND)
      JSTART = -2
      GO TO 250
C ITASK = 5.  See if TCRIT was reached and jump to exit. ---------------
 350  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. 100.0D0*UROUND*HMX
C-----------------------------------------------------------------------
C Block G.
C The following block handles all successful returns from DLSODI.
C if ITASK .ne. 1, Y is loaded from YH and T is set accordingly.
C ISTATE is set to 2, and the optional outputs are loaded into the
C work arrays before returning.
C-----------------------------------------------------------------------
 400  DO 410 I = 1,N
 410    Y(I) = RWORK(I+LYH-1)
      T = TN
      IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GO TO 420
      IF (IHIT) T = TCRIT
 420  ISTATE = 2
      IF (KFLAG .EQ. -3) ISTATE = 3
      RWORK(11) = HU
      RWORK(12) = H
      RWORK(13) = TN
      IWORK(11) = NST
      IWORK(12) = NRE
      IWORK(13) = NJE
      IWORK(14) = NQU
      IWORK(15) = NQ
      RETURN
C-----------------------------------------------------------------------
C Block H.
C The following block handles all unsuccessful returns other than
C those for illegal input.  First the error message routine is called.
C If there was an error test or convergence test failure, IMXER is set.
C Then Y is loaded from YH and T is set to TN.
C The optional outputs are loaded into the work arrays before returning.
C-----------------------------------------------------------------------
C The maximum number of steps was taken before reaching TOUT. ----------
 500  MSG = 'DLSODI-  At current T (=R1), MXSTEP (=I1) steps   '
      CALL XERRWD (MSG, 50, 201, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '      taken on this call before reaching TOUT     '
      CALL XERRWD (MSG, 50, 201, 0, 1, MXSTEP, 0, 1, TN, 0.0D0)
      ISTATE = -1
      GO TO 580
C EWT(i) .le. 0.0 for some i (not at start of problem). ----------------
 510  EWTI = RWORK(LEWT+I-1)
      MSG = 'DLSODI-  At T (=R1), EWT(I1) has become R2 .le. 0.'
      CALL XERRWD (MSG, 50, 202, 0, 1, I, 0, 2, TN, EWTI)
      ISTATE = -6
      GO TO 590
C Too much accuracy requested for machine precision. -------------------
 520  MSG = 'DLSODI-  At T (=R1), too much accuracy requested  '
      CALL XERRWD (MSG, 50, 203, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '      for precision of machine..  See TOLSF (=R2) '
      CALL XERRWD (MSG, 50, 203, 0, 0, 0, 0, 2, TN, TOLSF)
      RWORK(14) = TOLSF
      ISTATE = -2
      GO TO 590
C KFLAG = -1.  Error test failed repeatedly or with ABS(H) = HMIN. -----
 530  MSG = 'DLSODI-  At T(=R1) and step size H(=R2), the error'
      CALL XERRWD (MSG, 50, 204, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '      test failed repeatedly or with ABS(H) = HMIN'
      CALL XERRWD (MSG, 50, 204, 0, 0, 0, 0, 2, TN, H)
      ISTATE = -4
      GO TO 570
C KFLAG = -2.  Convergence failed repeatedly or with ABS(H) = HMIN. ----
 540  MSG = 'DLSODI-  At T (=R1) and step size H (=R2), the    '
      CALL XERRWD (MSG, 50, 205, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '      corrector convergence failed repeatedly     '
      CALL XERRWD (MSG, 50, 205, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '      or with ABS(H) = HMIN   '
      CALL XERRWD (MSG, 30, 205, 0, 0, 0, 0, 2, TN, H)
      ISTATE = -5
      GO TO 570
C IRES = 3 returned by RES, despite retries by DSTODI. -----------------
 550  MSG = 'DLSODI-  At T (=R1) residual routine returned     '
      CALL XERRWD (MSG, 50, 206, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '      error IRES = 3 repeatedly.        '
      CALL XERRWD (MSG, 40, 206, 0, 0, 0, 0, 1, TN, 0.0D0)
      ISTATE = -7
      GO TO 590
C DAINVG failed because matrix A was singular. -------------------------
 560  IER = -IER
      MSG='DLSODI- Attempt to initialize dy/dt failed:  Matrix A is    '
      CALL XERRWD (MSG, 60, 207, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '      singular.  DGEFA or DGBFA returned INFO = I1'
      CALL XERRWD (MSG, 50, 207, 0, 1, IER, 0, 0, 0.0D0, 0.0D0)
      ISTATE = -8
      RETURN
C DAINVG failed because RES set IRES to 2 or 3. ------------------------
 565  MSG = 'DLSODI-  Attempt to initialize dy/dt failed       '
      CALL XERRWD (MSG, 50, 208, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '      because residual routine set its error flag '
      CALL XERRWD (MSG, 50, 208, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '      to IRES = (I1)'
      CALL XERRWD (MSG, 20, 208, 0, 1, IER, 0, 0, 0.0D0, 0.0D0)
      ISTATE = -8
      RETURN
C Compute IMXER if relevant. -------------------------------------------
 570  BIG = 0.0D0
      IMXER = 1
      DO 575 I = 1,N
        SIZE = ABS(RWORK(I+LACOR-1)*RWORK(I+LEWT-1))
        IF (BIG .GE. SIZE) GO TO 575
        BIG = SIZE
        IMXER = I
 575    CONTINUE
      IWORK(16) = IMXER
C Compute residual if relevant. ----------------------------------------
 580  LYD0 = LYH + NYH
      DO 585  I = 1,N
         RWORK(I+LSAVR-1) = RWORK(I+LYD0-1)/H
 585     Y(I) = RWORK(I+LYH-1)
      IRES = 1
      CALL RES (NEQ, TN, Y, RWORK(LSAVR), YDOTI, IRES )
      NRE = NRE + 1
      IF (IRES .LE. 1) GO TO 595
      MSG = 'DLSODI-  Residual routine set its flag IRES       '
      CALL XERRWD (MSG, 50, 210, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '      to (I1) when called for final output.       '
      CALL XERRWD (MSG, 50, 210, 0, 1, IRES, 0, 0, 0.0D0, 0.0D0)
      GO TO 595
C Set Y vector, T, and optional outputs. -------------------------------
 590  DO 592 I = 1,N
 592    Y(I) = RWORK(I+LYH-1)
 595  T = TN
      RWORK(11) = HU
      RWORK(12) = H
      RWORK(13) = TN
      IWORK(11) = NST
      IWORK(12) = NRE
      IWORK(13) = NJE
      IWORK(14) = NQU
      IWORK(15) = NQ
      RETURN
C-----------------------------------------------------------------------
C Block I.
C The following block handles all error returns due to illegal input
C (ISTATE = -3), as detected before calling the core integrator.
C First the error message routine is called.  If the illegal input
C is a negative ISTATE, the run is aborted (apparent infinite loop).
C-----------------------------------------------------------------------
 601  MSG = 'DLSODI-  ISTATE (=I1) illegal.'
      CALL XERRWD (MSG, 30, 1, 0, 1, ISTATE, 0, 0, 0.0D0, 0.0D0)
      IF (ISTATE .LT. 0) GO TO 800
      GO TO 700
 602  MSG = 'DLSODI-  ITASK (=I1) illegal. '
      CALL XERRWD (MSG, 30, 2, 0, 1, ITASK, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 603  MSG = 'DLSODI-  ISTATE .gt. 1 but DLSODI not initialized.'
      CALL XERRWD (MSG, 50, 3, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 604  MSG = 'DLSODI-  NEQ (=I1) .lt. 1     '
      CALL XERRWD (MSG, 30, 4, 0, 1, NEQ(1), 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 605  MSG = 'DLSODI-  ISTATE = 3 and NEQ increased (I1 to I2). '
      CALL XERRWD (MSG, 50, 5, 0, 2, N, NEQ(1), 0, 0.0D0, 0.0D0)
      GO TO 700
 606  MSG = 'DLSODI-  ITOL (=I1) illegal.  '
      CALL XERRWD (MSG, 30, 6, 0, 1, ITOL, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 607  MSG = 'DLSODI-  IOPT (=I1) illegal.  '
      CALL XERRWD (MSG, 30, 7, 0, 1, IOPT, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 608  MSG = 'DLSODI-  MF (=I1) illegal.    '
      CALL XERRWD (MSG, 30, 8, 0, 1, MF, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 609  MSG = 'DLSODI-  ML(=I1) illegal: .lt. 0 or .ge. NEQ(=I2) '
      CALL XERRWD (MSG, 50, 9, 0, 2, ML, NEQ(1), 0, 0.0D0, 0.0D0)
      GO TO 700
 610  MSG = 'DLSODI-  MU(=I1) illegal: .lt. 0 or .ge. NEQ(=I2) '
      CALL XERRWD (MSG, 50, 10, 0, 2, MU, NEQ(1), 0, 0.0D0, 0.0D0)
      GO TO 700
 611  MSG = 'DLSODI-  MAXORD (=I1) .lt. 0  '
      CALL XERRWD (MSG, 30, 11, 0, 1, MAXORD, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 612  MSG = 'DLSODI-  MXSTEP (=I1) .lt. 0  '
      CALL XERRWD (MSG, 30, 12, 0, 1, MXSTEP, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 613  MSG = 'DLSODI-  MXHNIL (=I1) .lt. 0  '
      CALL XERRWD (MSG, 30, 13, 0, 1, MXHNIL, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 614  MSG = 'DLSODI-  TOUT (=R1) behind T (=R2)      '
      CALL XERRWD (MSG, 40, 14, 0, 0, 0, 0, 2, TOUT, T)
      MSG = '      Integration direction is given by H0 (=R1)  '
      CALL XERRWD (MSG, 50, 14, 0, 0, 0, 0, 1, H0, 0.0D0)
      GO TO 700
 615  MSG = 'DLSODI-  HMAX (=R1) .lt. 0.0  '
      CALL XERRWD (MSG, 30, 15, 0, 0, 0, 0, 1, HMAX, 0.0D0)
      GO TO 700
 616  MSG = 'DLSODI-  HMIN (=R1) .lt. 0.0  '
      CALL XERRWD (MSG, 30, 16, 0, 0, 0, 0, 1, HMIN, 0.0D0)
      GO TO 700
 617  MSG='DLSODI-  RWORK length needed, LENRW (=I1), exceeds LRW (=I2)'
      CALL XERRWD (MSG, 60, 17, 0, 2, LENRW, LRW, 0, 0.0D0, 0.0D0)
      GO TO 700
 618  MSG='DLSODI-  IWORK length needed, LENIW (=I1), exceeds LIW (=I2)'
      CALL XERRWD (MSG, 60, 18, 0, 2, LENIW, LIW, 0, 0.0D0, 0.0D0)
      GO TO 700
 619  MSG = 'DLSODI-  RTOL(=I1) is R1 .lt. 0.0       '
      CALL XERRWD (MSG, 40, 19, 0, 1, I, 0, 1, RTOLI, 0.0D0)
      GO TO 700
 620  MSG = 'DLSODI-  ATOL(=I1) is R1 .lt. 0.0       '
      CALL XERRWD (MSG, 40, 20, 0, 1, I, 0, 1, ATOLI, 0.0D0)
      GO TO 700
 621  EWTI = RWORK(LEWT+I-1)
      MSG = 'DLSODI-  EWT(I1) is R1 .le. 0.0         '
      CALL XERRWD (MSG, 40, 21, 0, 1, I, 0, 1, EWTI, 0.0D0)
      GO TO 700
 622  MSG='DLSODI-  TOUT(=R1) too close to T(=R2) to start integration.'
      CALL XERRWD (MSG, 60, 22, 0, 0, 0, 0, 2, TOUT, T)
      GO TO 700
 623  MSG='DLSODI-  ITASK = I1 and TOUT (=R1) behind TCUR - HU (= R2)  '
      CALL XERRWD (MSG, 60, 23, 0, 1, ITASK, 0, 2, TOUT, TP)
      GO TO 700
 624  MSG='DLSODI-  ITASK = 4 or 5 and TCRIT (=R1) behind TCUR (=R2)   '
      CALL XERRWD (MSG, 60, 24, 0, 0, 0, 0, 2, TCRIT, TN)
      GO TO 700
 625  MSG='DLSODI-  ITASK = 4 or 5 and TCRIT (=R1) behind TOUT (=R2)   '
      CALL XERRWD (MSG, 60, 25, 0, 0, 0, 0, 2, TCRIT, TOUT)
      GO TO 700
 626  MSG = 'DLSODI-  At start of problem, too much accuracy   '
      CALL XERRWD (MSG, 50, 26, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG='      requested for precision of machine..  See TOLSF (=R1) '
      CALL XERRWD (MSG, 60, 26, 0, 0, 0, 0, 1, TOLSF, 0.0D0)
      RWORK(14) = TOLSF
      GO TO 700
 627  MSG = 'DLSODI-  Trouble in DINTDY.  ITASK = I1, TOUT = R1'
      CALL XERRWD (MSG, 50, 27, 0, 1, ITASK, 0, 1, TOUT, 0.0D0)
C
 700  ISTATE = -3
      RETURN
C
 800  MSG = 'DLSODI-  Run aborted.. apparent infinite loop.    '
      CALL XERRWD (MSG, 50, 303, 2, 0, 0, 0, 0, 0.0D0, 0.0D0)
      RETURN
C----------------------- End of Subroutine DLSODI ----------------------
      END
