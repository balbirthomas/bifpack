CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC
CCC      This is    packON  , version   12 Aug 1997
CCC
CCC      for interactive work, using
CCC      the PACKD and/or the PACKA files.
CCC          -----            -----
CCC
CCC      A SPECIAL VERSION   I N T E R A   FOR CONTINUATION OF
CCC      ALGEBRAIC EQUATIONS IS CURRENTLY ATTACHED TO  P A C K A 1 .
CCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                   
                      

c    BIFPACK -- a package for Bifurcation, Continuation, and Stability Analysis.
c    Copyright (C) 1983, 1999   R. Seydel

c    This program is free software; you can redistribute it and/or modify
c    it under the terms of the GNU General Public License as published by
c    the Free Software Foundation; either version 2 of the License, or
c    (at your option) any later version.

c    This program is distributed in the hope that it will be useful,
c    but WITHOUT ANY WARRANTY; without even the implied warranty of
c    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c    GNU General Public License for more details.

c    You should have received a copy of the GNU General Public License
c    along with this program; if not, write to the Free Software
c    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

c========================================================================                      
                      
      SUBROUTINE BIFLAB (NOB,NOS,N,IFLAGSOL,X,EPS,INTE,JACOBI
     1   ,WINDOW,UNIT2,UNIT3,UNIT31,UNIT4,UNITLST,UNITOUT)
      DOUBLE PRECISION EPS,X,XA,XHOPF,Y,T,S,T1,T2,ETACO,ZETACO,PARCO,
     1   CC,DELPAR,DIST,HINIT,HMAX,PAR,PERIOD,EVERE,EVEIM,
     2   PI2,STEP,STEPMA,STEPMI,STEPMP,TEND,WEIGHT,EVARE,EVAIM,
     3   WINDOW,YFINAL,EPSONE,BOUND2,Y2,Y0BUFF,DECR,EVAIMGUESS,
     4   SECPAR1,SECPARINC,PARSTOR,PERSTOR,STEPBUFF
      CHARACTER*9 Y0TEXT
      CHARACTER*2 JUDGEBUFF
      INTEGER HFLAG,HEAD1BUFF,COUNTRUNS,COUNTRUNMAX
     1        ,UNIT2,UNIT3,UNIT31,UNIT4,UNITLST,UNITOUT
      PARAMETER (NDIML=14,NDIMBUFF=6)
      DIMENSION X(NDIML),XA(NDIML),Y(NDIML,7),T(7),WINDOW(2)
     1      ,Y2(NDIML,7),Y0BUFF(NDIML,NDIMBUFF),Y0TEXT(NDIMBUFF)
     2      ,EVERE(NDIML),EVEIM(NDIML),XHOPF(NDIML)
     3      ,HEAD1BUFF(3,NDIMBUFF),JUDGEBUFF(NDIMBUFF)
     4      ,STEPBUFF(NDIMBUFF)
C   DIMENSIONS AT LEAST  X(N+2),Y(N+2,M),T(M)
      EXTERNAL RABZ,FABZ,INTE,FCN,FCNLIN
      COMMON /CONTRO/ ETACO,ZETACO,PARCO, NCO,N1CO,N2CO,
     1                KCO,IND2CO,NSECFLAG,NFLAG
C
C  THIS ROUTINE MONITORS THE USE OF VARIOUS OPTIONS/TOOLS 
C  FOR CALCULATING PERIODIC SOLUTIONS OR INVESTIGATING STATIONARY
C  SOLUTIONS OF AUTONOMOUS SYSTEMS.
C
C  INPUT AND OUTPUT ARE KEPT ON A MODERATE LEVEL, SOME CONSTANTS ARE
C  SET TO AVERAGE VALUES.
C
C     FEEL FREE TO CHANGE THIS ROUTINE ACCORDING TO YOUR SPECIFIC
C     EXAMPLE / HARDWARE / SOFTWARE / NEEDS
C
C  INPUT DATA ARE READ FROM FILE 'START' (UNIT 3), RESULTS ARE WRITTEN
C  TO THIS FILE (VECTOR OF INITIAL VALUES, PARAMETER, PERIOD).

C     INPUT FROM CALLING SEQUENCE:
C     -------------------------
C      N     : NUMBER OF EQUATIONS (=NCO)
C      MODUS : FIXES THE DESIRED OPTION, SEE BELOW.
C      EPS   : REQUIRED ACCURACY OF SOLUTIONS
C      INTE  : NAME OF INTEGRATOR TO BE USED
C      JACOBI: =1 IF THE PARTIAL DERIVATIVES OF FCN ARE PROVIDED
C              IN ROUTINE FCNLIN.   ELSE =0  .
C
C     INPUT FROM LOGICAL UNIT 3:
C     --------------------------
C      N+2  OR  N+1  REAL NUMBERS (WILL BE UPDATED):
C      X(I),I=1,..,N : INITIAL VALUES OR STEADY-STATE SOLUTION
C      X(N+1)        : PARAMETER VALUE
C      X(N+2)        : PERIOD (NOT NEEDED FOR MODUS=7)
C
C     INTERNAL INPUT-OUTPUT FILE  FOR "INDEX" (=KCO)
C     --------------------------
C      THE INDEX OF THE NORMALIZATION IS STORED ON UNIT 2 (FORMAT I2).
C      THIS FILE MIGHT BE MANIPULATED, IF NECESSARY. IT MAY ALSO BE 
C      EMPTY.
C
C
C     MODUS  CONTROLS THE SELECTED OPTION:
C     ----------------------------------
C     =0   VARIOUS INITIAL PREPARATIONS; SUCH AS:
C          choosing the extent of interaction,  or
C          check the jacobian of subroutine FCNLIN
C     =1   CONTINUATION ALONG PERIODIC BRANCH WITH RESPECT TO PARAMETER
C     =2   CALL OF ROUTINE  HOPF  (TRANSITION FROM STATIONARY
C                                  TO PERIODIC SOLUTIONS)
C     =3   CALCULATION OF DATA FOR A PHASE DIAGRAM  OR  SIMULATION.
C          THIS OPTION IS NOT CONFINED TO PERIODIC SOLUTIONS.
C          HERE THE TERM "PERIOD" REFERS TO THE VALUE IN X(N+2).
C     =4   CHECK/RECONFIRM/RECALCULATE THE PRESENT PERIODIC SOLUTION
C                 IN DETAIL:
C            PARAMETRIZATION WITH RESPECT TO INDEX-TH INITIAL VALUE.
C            INCLUDES STABILITY ANALYSIS AND OUTPUT OF DATA FOR
C            A PHASE DIAGRAM. MAY BE HELPFUL ALSO FOR CALCULATING
C            A SOLUTION AFTER HAVING MANIPULATED THE INITIAL VALUES
C            IN FILE "START" (LOGICAL UNIT 3).
C     =5   simplified jumping on emanating branch,
C           starts from bifurcation point.
C           linearization not needed. WORKS WITH SIMPLE PROBLEMS.
C     =7   CONTINUATION AND STABILITY ANALYSIS OF A STEADY-STATE 
C           BRANCH. A MORE COMFORTABLE SPECIAL VERSION FOR THIS CASE IS
C           I N T E R A  , CURRENTLY ATTACHED TO  P A C K A . 
C           IN CASE YOU WANT TO WORK WITH  I N T E R A  RATHER THAN WITH
C           B I F L A B , YOU MAY LIKE TO PLACE A 'C' INTO THE FIRST 
C           COLUMN OF THE STATEMENT 'CALL CONTA ... ' , AND YOU CAN 
C           DISPENSE HERE WITH  PACKA1 AND NEWTON'S METHOD.
C     NEEDED SOFTWARE:
C     ---------------
C     O.D.E. INTEGRATOR FOR MODUS=3 
C     O.D.E. BVP SOLVER FOR MODUS=1,2,4,5 
C     PACKD  FOR MODUS=1,2,4,5
C     PACKA1  FOR MODUS=7
C     NEWTON-LIKE METHOD FOR MODUS=7
C     LINEARIZATION OF FCN3 (IN SUBROUTINE FCNLIN) FOR
C                        MODUS=2  (FOR JBIF.GE.2)
C*************************************************************************

      WRITE (*,*) ' BIFPACK: entering routine BIFLAB,'
      WRITE (UNITLST,*) ' BIFPACK: entering routine BIFLAB,'
      INTERACT=1
      WRITE (*,*) ' with the default medium level of interactions.'
c  some typical values chosen for the case of INTERACT=0:
      COUNTRUNMAX=5
      MAXCONAIM=10

      EPSONE=1.d-3
      NCO=N
      N1=N+1
      N1CO=N1
      ETACO=0.d0
      ZETACO=0.d0
      N2=N+2
      N2CO=N2
      COUNTRUNS=0
      IF (NSECFLAG.EQ.1) WRITE (*,*) 'second parameter is',PARCO
      IF (NSECFLAG.EQ.1) WRITE (UNITLST,*) 'second parameter is',PARCO
      LAST=0
      ISTART=1
      CALL CHECKDATA (N,X,IFLAGSOL,IFLAGPER,ISTART,INTERACT,
     1               NOB,NOS,NEXTMODUS,STEP,Y0BUFF,NUMBUFF,
     2               HEAD1BUFF,Y0TEXT,STEPBUFF,JUDGEBUFF,JBUFFER
     3               ,UNIT3,UNIT31,UNIT4)
      ISTART=0

c--------------------------------------------------------------------
101   CONTINUE
      IFLAGPER=2
      MODUSST=MODUS
      IF (INTERACT.EQ.0) THEN
         COUNTRUNS=COUNTRUNS+1
         IF (COUNTRUNS.GT.COUNTRUNMAX) THEN
            write (*,*) 'BIFPACK stops, because more than ',
     1        'COUNTRUNMAX= ',COUNTRUNMAX,' runs performed.'
            write (UNITLST,*) 'BIFPACK stops, because more than ',
     1        'COUNTRUNMAX= ',COUNTRUNMAX,' runs performed.'
            go to 999
            endif
         IF (MODUSST.EQ.2) THEN
            MODUS=1
            GO TO 16
            ENDIF
         CALL CHECKDATA (N,X,IFLAGSOL,IFLAGPER,ISTART,INTERACT,
     1               NOB,NOS,NEXTMODUS,STEP,Y0BUFF,NUMBUFF,
     2               HEAD1BUFF,Y0TEXT,STEPBUFF,JUDGEBUFF,JBUFFER
     3               ,UNIT3,UNIT31,UNIT4)
         IF (NEXTMODUS.GT.0) THEN
            MODUS=NEXTMODUS
            GO TO 16
            ENDIF
         ENDIF
      WRITE (*,610)
      WRITE (*,*) ' ENTER CHOICE:         '
      READ (*,*) MODUS
      IF (MODUS.LT.0) GO TO 999
      IF (MODUS.EQ.6) GO TO 6
      IF (MODUS.EQ.8) GO TO 8
      IF (MODUS.EQ.0) THEN
        WRITE (*,*) ' Initial Preparations:'
        WRITE (*,*) ' enter -1  to exit, '
        WRITE (*,*) ' enter 0  for a lower level of interaction, '
        WRITE (*,*) ' enter 1  for a medium level of interaction, '
        WRITE (*,*) ' enter 2  for a high level of interaction, '
        WRITE (*,*) ' enter 9  to check the Jacobian. '
        WRITE (*,*) 'ENTER now:'
        READ (*,*) IPREP
        IF (IPREP.LT.0.OR.IPREP.GT.9) GO TO 101
        IF (IPREP.GT.2.AND.IPREP.LT.9) go to 101
        IF (IPREP.EQ.0) THEN
           INTERACT=0
           WRITE (*,*) 'You have chosen a low level of interactions.'
           GO TO 101 
           ENDIF
        IF (IPREP.EQ.1) THEN
           INTERACT=1
           WRITE (*,*) 'You have chosen a medium level of interactions.'
           GO TO 101 
           ENDIF
        IF (IPREP.EQ.2) THEN
           INTERACT=2
           WRITE (*,*) 'You have chosen a high level of interactions.'
           GO TO 101
           ENDIF
        IF (JACOBI.EQ.0) THEN
           WRITE (*,*) ' ... not available without linearization.'
           GO TO 101
           ENDIF
        ENDIF
      WRITE (UNITLST,602) MODUS
      WRITE (UNITLST,*) 'Chosen error tolerance EPS is',EPS
C     EACH RUN CREATES A NEW PRINT-FILE (YOU MAY LIKE TO CHANGE IT):
C      REWIND UNITLST
      REWIND UNIT3
      REWIND UNIT2
      READ (UNIT2,604, END=12) INDEX
      GO TO 13
12    INDEX=0
13    REWIND UNIT2

16    CONTINUE 
      KCO=INDEX
      IF (KCO.EQ.0) KCO=1
      IND2CO=1
c      CALL CHECKDATA (N,X,IFLAGSOL,IFLAGPER,ISTART,INTERACT,
c     1               NOB,NOS,NEXTMODUS,STEP,Y0BUFF,NUMBUFF,
c     2               HEAD1BUFF,Y0TEXT,STEPBUFF,JUDGEBUFF,JBUFFER
c     3               ,UNIT3,UNIT31,UNIT4)
      IF (INTERACT.EQ.0) WRITE (*,*) 'Next chosen MODUS:',MODUS 
      IF (MODUS.EQ.2) GO TO 2
      IF (MODUS.EQ.3) GO TO 3
      IF (MODUS.EQ.7) GO TO 7
      DO 10 I=1,N2
10    Y(I,1)=X(I)
C--  IMPLEMENTING EQUIDISTANT NODES FOR MULTIPLE SHOOTING:
      IF (INTERACT.LE.1) THEN
          M=7
        ELSE
          WRITE (*,*) 'Enter number of shooting nodes (.le.7):'
          READ (*,*) M
        ENDIF
      T(1)=0.d0
      T1=1.d0/DBLE(M-1)
      DO 11 J=2,M
11    T(J)=T(J-1)+T1


C--  CALCULATING OF Y ARRAY:
      IF (MODUS.EQ.5) THEN
         WRITE (*,*) ' Enter 2 for PERIOD DOUBLING, else enter 1'
         READ (*,*) IPER
         IF (IPER.EQ.2) Y(N2,1)=Y(N2,1)*2.d0
         END IF
      NFLAG=1 
      CALL DATEN (FABZ,INTE,N2,M,T,Y,EPS)
      IF (MODUS.EQ.0) GO TO 19

C   DETERMINATION OF INDEX:
      IF (INDEX.EQ.0) THEN
          CC=0.d0
          DO 14 I=1,N
          IF (ABS(Y(I,1)).GT.CC) THEN
              CC=ABS(Y(I,1))
              INDEX=I
              END IF
14        CONTINUE
          WRITE (UNIT2,604) INDEX
          END IF

      IF (MODUS.EQ.1) GO TO 1
      IF (MODUS.EQ.4) GO TO 4
      RETURN

CHECKING THE JACOBIAN 
19    WRITE (*,*) ' BIFPACK is now checking the Jacobian.'
      WRITE (UNITLST,*) ' BIFPACK is now checking the Jacobian.'
      CALL CHECKD (T1,X,N,N,FCN,FCNLIN,ICALL,UNITLST)
      GO TO 101


C************************************************************
C*
C*    A CONTINUATION WITH RESPECT TO PARAMETER WILL BE CARRIED OUT:
C*    (NO CONTROL OF STEPSIZE)
C
1     CONTINUE
      EPSONE=1.d-3
      IF (INTERACT.EQ.0) THEN
         MAXCON=MAXCONAIM
        ELSE
         WRITE (*,*) ' Enter number of contin. steps:    '
         READ (*,*) MAXCON
        ENDIF 
c      IFLAGPER=2

      IF ((MODUSST.EQ.2.AND.MAX.GT.0).OR.
     1    MODUSST.EQ.1) THEN
          DELPAR=Y(N+1,1)-Y2(n+1,1)
          WRITE (*,*) 'Last increment in lambda was:',DELPAR
          IF (INTERACT.EQ.0) GO TO 15
          ENDIF
      IF (INTERACT.GT.0) THEN
          WRITE (*,*) ' Enter stepsize (real number):     '
          READ (*,*) DELPAR
          ENDIF

15    ITAIM=6
C The following may be activated, if the amount of work for
calculating the monodromy matrix is to be reduced:
Choose a large value for EPSONE.
C      WRITE (*,*) ' enter EPSONE :         '
C      READ (*,*) EPSONE
C      WRITE (*,*) ' enter ITAIM:'
C      READ (*,*) ITAIM
      IF (INTERACT.LE.1) THEN
          LEVEL=2
       ELSE
        WRITE (*,*) ' ENTER 2 for variable stepsize, 0 for fixed steps'
        READ (*,*) LEVEL
       ENDIF
      IF (INTERACT.GT.0) THEN       
        write (*,*) ' last No. branch /solution is:',nob,nos
        write (*,*) 'Enter 0 for continuation of numbering,'
        write (*,*) 'enter 1 for new numbers:'
        read (*,*) nosnob
        if (nosnob.ne.0) then
          write (*,*) ' Enter number of branch:'
          read (*,*) nob
          nos=0
          endif
          ENDIF
      IOLD=0
      IF (MODUS.EQ.MODUSST) IOLD=1
      CALL CONTD (NOB,NOS,N,M,T,Y,Y2,IOLD,DELPAR,
     1     MAXCON,EPS,LEVEL,ITAIM,INTERACT,JACOBI,INTE,INDEX,EPSONE,
     2     Y0BUFF,NUMBUFF,Y0TEXT,HEAD1BUFF,STEPBUFF,JUDGEBUFF,
     3     UNITLST,UNITOUT,UNIT3,UNIT31)
      DO 18 I=1,N+2
18    X(I)=Y(I,1)

      MODUSST=MODUS
      GO TO 101



C***********************************************************
C*
C*    CALCULATION OF HOPF BIFURCATION:
C
2     CONTINUE
      IF (INTERACT.GT.0) THEN
        write (*,*) ' last No. branch /solution is:',nob,nos
        write (*,*) 'Enter 0 for continuation of numbering,'
        write (*,*) 'enter 1 for new numbers:'
        read (*,*) nosnob
        if (nosnob.ne.0) then
          write (*,*) ' Enter number of branch:'
          read (*,*) nob
          nos=0
          endif
c        write (*,*) 'p=',x(n+1),' . T=',x(n+2)
        ELSE
          NOB=NOB+1
          NOS=0
        ENDIF
      IF (ABS(X(N+2)).LT.EPS) THEN
         WRITE (*,*) ' your initial guess has zero period.'
         GO TO 101
         ENDIF
      EVAIMGUESS=6.28318530718d0/X(N+2)
      IFLAGVEC=0
calculating a guess for the eigenvector
      DO 24 I=1,N
      EVERE(I)=1.d0
24    EVEIM(I)=0.2d0
c  better starting data are obtained when inverse iteration is applied.

C  Then the linearization is needed in subroutine  FCNLIN .
      IF (JACOBI.EQ.1) THEN
         DO 23 I=1,N
23       XHOPF(I)=X(I)
         NSECFLAG1=0
         NSECFLAG2=0
         SECPAR1=PARCO
231      HFLAG=1
         PAR=X(N+1)
         NSECFLAG1=NSECFLAG1+1
         IF (NSECFLAG1.EQ.1) WRITE (*,*)
     1                 'BIFPACK is entering routine HOPFA.'
         WRITE (UNITLST,*) ' '
         CALL HOPFA (N,XHOPF,PAR,EVAIMGUESS,EVERE,EVEIM,PERIOD,HFLAG,
     1               EPS,UNITLST)
         IF (HFLAG.EQ.1) THEN
            JBIF=1
            IF (NSECFLAG1.EQ.1) THEN
               DO 232 I=1,N
232            X(I)=XHOPF(I)
               PERIOD=6.28318530718d0/EVAIMGUESS
               PARSTOR=PAR
               PERSTOR=PERIOD
               X(N+1)=PAR
               X(N+2)=PERIOD
               ENDIF
            IF (NSECFLAG.EQ.1.AND.NSECFLAG1.EQ.1) THEN
               WRITE (*,*) 'Control of second parameter.'
               WRITE (*,*) '  Vary second parameter?', 
     1               '  (1=yes, 0=no )  ENTER:'
               READ (*,*) NSECFLAG2
               IF (NSECFLAG2.EQ.1) THEN
                 WRITE (*,*) 'Enter increment in second parameter:'
                 READ (*,*) SECPARINC
                 WRITE (*,*) 'Enter the number of solutions that are to'
                 WRITE (*,*) ' be calculated along Hopf curve:'
                 READ (*,*) NSECFLAG3
                 IF (NSECFLAG3.LT.0) NSECFLAG3=5
                 ENDIF  
               ENDIF
            NOS=1
            IF (NSECFLAG.NE.1) WRITE (UNITOUT,9300)
     1       'R',NOB,NOS,N,1,'d','HB',(X(I),I=1,N),PAR,PERIOD 
            IF (NSECFLAG.EQ.1.AND.nsecflag1.eq.1) 
     1          WRITE (UNITOUT,9300)
     2       'R',NOB,NOS,N,2,'d','HB',(X(I),I=1,N),PAR,PARCO,PERIOD
            IF (NSECFLAG2.EQ.1) WRITE (UNITOUT,9300)
     1       'R',0,0,N,2,'d','HB',(X(I),I=1,N),PAR,PARCO,PERIOD
            IF (NSECFLAG2.EQ.1) WRITE (UNITLST,*)
     1          '      second parameter is',PARCO  
            if (NSECFLAG2.EQ.1.AND.NSECFLAG1.lt.NSECFLAG3) then
                parco=parco+SECPARINC
                go to 231
                endif
            IFLAGVEC=1
            ENDIF 
            IF (NSECFLAG2.EQ.1.AND.NSECFLAG1.LT.NSECFLAG3)
     1          write (*,*) ' Not all solutions calculated.' 
            IF (NSECFLAG2.EQ.1) THEN
              WRITE (*,*) 'Re-load the previous value of sec. par.?'
              WRITE (*,*) ' Enter 1 for yes, 0 for no:'
              READ (*,*) NSECFLAG2
              IF (NSECFLAG2.EQ.0) THEN 
                DO 25 I=1,N
25              X(I)=XHOPF(I)
                PERIOD=6.28318530718d0/EVAIMGUESS
                X(N+1)=PAR
                X(N+2)=PERIOD
                ENDIF
              IF (NSECFLAG2.EQ.1) THEN
                PAR=PARSTOR
                PERIOD=PERSTOR  
                PARCO=SECPAR1 
                ENDIF
              ENDIF
         ENDIF
      IF (NSECFLAG.EQ.1) WRITE (*,*) 'second parameter is',PARCO
      WRITE (UNITLST,*) ' '
      WRITE (UNITLST,*) ' For the following calculations, the '
      IF (NSECFLAG.EQ.1) WRITE (UNITLST,*) ' second parameter is',PARCO
C  JBIF=2  IF BIFURCATION POINT IS TO BE CALCULATED.
C  THEN THE LINEARIZATION IS NEEDED IN SUBROUTINE  FCNLIN .
      IF (INTERACT.LE.1) THEN
         JHOPF=0
        ELSE
      WRITE (*,*) ' STANDARD TREATMENT means CALCULATING FIVE OSCILL.'
      WRITE (*,*) ' For STANDARD enter 0, for SPECIFIC WISHES enter 1'
      READ (*,*) JHOPF
      ENDIF
      IF (JHOPF.NE.1) THEN
         JBIF=0
         IF (N.EQ.2.AND.HFLAG.LT.1) JBIF=2
         IF (JACOBI.EQ.0) JBIF=-1
         IF (IFLAGVEC.EQ.1) JBIF=1
         DECR=0.01d0
         MAX=5
         GO TO 21
         END IF
      WRITE (*,*) ' how many SMALL oscillations to calculate? ENTER:'
      READ (*,*) MAX
      IF (HFLAG.LT.1) THEN
        WRITE (*,*) ' enter modus JBIF for  H O P F :          '
        READ (*,*) JBIF
        IF (JBIF.GT.2.OR.JBIF.LT.-1) JBIF=2
        IF (JACOBI.EQ.0) JBIF=-1
        ENDIF
      WRITE (*,*) 'Enter the decrement in amplitude, a small real'
      WRITE (*,*) '  positive number such as 0.02 .    ENTER:'
      READ (*,*) DECR
C    MAX  small-amplitude oscillations are to be calculated.
                      
21    PAR=X(N+1)
      PERIOD=X(N+2)
      MAXSTO=MAX
      CALL HOPF (NOB,NOS,N,X,PAR,PERIOD,EPS,INTE,T,M,Y,IFLAGVEC,
     1    EVERE,EVEIM,MAX,DECR,Y2,JBIF,INDEX,UNITLST,UNITOUT,UNIT31)
      IF (MAX.EQ.-1) GO TO 991
      IF (MAX.GT.0) THEN
         REWIND UNIT3
         WRITE (UNIT3,606) (Y(I,1),I=1,N2)
         DO 22 i=1,N+2
22       X(I)=Y(I,1)
         WRITE (UNIT2,604) INDEX
         END IF
      GO TO 101


C**************************************************************
C*
C*    DATA FOR PHASE DIAGRAM
C*
3      CONTINUE
C     CHOOSE NUMBER OF POINTS PER PERIOD:
      NUMBER=50
      IF (MODUS.NE.4) THEN
	  WRITE (*,*) ' Current final time / period is',X(N+2)
	  WRITE (*,*) ' Enter  1  for integr. of ONE period,      '
	  WRITE (*,*) ' enter  N  for integr. of N periods:       '
          READ (*,*) NFAC
          ELSE
             NFAC=1
          END IF
      TEND=FLOAT(NFAC)
      HINIT=0.01d0*X(N2)
      TDEL=TEND/DBLE(NUMBER)
      T1=0.d0
      HMAX=TDEL
      WRITE (UNITLST,607) X(N+1),X(N2)
      WRITE (UNITLST,609)
      WRITE (UNITLST,605) T1,(X(I),I=1,N)
      NFLAG=1

      DO 35 I=1,NUMBER
      T2=T1+TDEL
      CALL INTE (N2,FABZ,T1,X,T2,EPS,HMAX,HINIT)
      IF (HINIT.eq.0.d0) THEN
         WRITE (*,*) 'Integration not successful'
         GO TO 101
         ENDIF
35    WRITE (UNITLST,605) T1,(X(K),K=1,N)

      WRITE (*,*) 'Simulation is done.'
      WRITE (*,*) 'Enter 1 if final value is to be saved in START,'
      WRITE (*,*) '   else enter 0:'
      READ (*,*) NFAC
      IF (NFAC.GT.0) THEN
                 REWIND UNIT3
                 WRITE (UNIT3,606) (X(I),I=1,N2)
                 END IF
      GO TO 101

C************************************************************
4     CONTINUE
      DIST=0.d0
      CALL CHECKPHASE (Y,INDEXPH)
      IF (INDEXPH.GT.0 .AND. INDEXPH.NE.INDEX) THEN
         INDEX=INDEXPH
         WRITE (*,*) 'Phase condition adapted to INDEX=',INDEXPH
         ENDIF
      WRITE (*,*) 'Phase condition:  DY(',INDEX,')=0.'
      GO TO 51

C************************************************************
currently not installed:
5     DIST=0.01d0
C     DIST=P*0.01d0   FOR A  P-PERCENT CHANGE IN THE  KCO-TH
C     INITIAL CONDITION
51    NFLAG=2
      PI2=6.2831853071796d0
      KCO=INDEX
      DIST=DIST*ABS(Y(KCO,1))
      DO 52 J=1,M
      Y(KCO,J)=Y(KCO,J)+DIST*COS(PI2*T(J))
      DO 52 I=1,N
52    Y2(I,J)=Y(I,J)
      ETACO=Y(KCO,1)
      IT=20
      CALL SOLVER (FABZ,RABZ,INTE,N2,M,T,Y,EPS,IT)
      WRITE (*,*) 'it=',it
      IF (IT.LT.0) GO TO 991
      WRITE (UNITLST,*) 'it=',it
      WRITE (UNITLST,607) Y(N+1,1),Y(N2,1)
      WRITE (*,607) Y(N+1,1),Y(N2,1)
      WRITE (UNITLST,605) (Y(I,1),I=1,N)
      IF (NSECFLAG.EQ.0) WRITE (UNITOUT,9300)
     1     'R',-1,0,n,1,' ','P ',(Y(I,1),I=1,N+1),Y(N2,1)
      IF (NSECFLAG.EQ.1) WRITE (UNITOUT,9300)
     1     'R',-1,0,n,2,' ','P ',(Y(I,1),I=1,N+1),PARCO,Y(N2,1)
      REWIND UNIT3
      WRITE (UNIT3,606) (Y(I,1),I=1,N2)
      IFLAG=0
      CALL STAPER (NOB,N,M,Y,EPS,Y2,EVARE,EVAIM,IFLAG,IND,
     1             JACOBI,INTE,EPSONE,Y0BUFF,NUMBUFF,Y0TEXT
     2            ,HEAD1BUFF,STEPBUFF,JUDGEBUFF,UNITLST,UNITOUT)
      IF (MODUS.EQ.4) THEN
          DO 53 I=1,N2
53        X(I)=Y(I,1)
          GO TO 3
          END IF
      GO TO 101

C*********************************************************************
C*
6     CALL CHECKDATA (N,X,IFLAGSOL,IFLAGPER,ISTART,INTERACT,
     1               NOB,NOS,NEXTMODUS,STEP,Y0BUFF,NUMBUFF,
     2               HEAD1BUFF,Y0TEXT,STEPBUFF,JUDGEBUFF,JBUFFER
     3               ,UNIT3,UNIT31,UNIT4)
      GO TO 101

C*********************************************************************
C*
C*    STABILITY ANALYSIS OF STEADY-STATE BRANCH
C
7     CONTINUE
      STEPMI=EPS*10.d0
      ITAIM=3
      ITMAX=2*ITAIM
      WEIGHT=10.d0
      IEXTRA=1
      IRUN=0
      ISTAB=1
      IPRINT=0
      BOUND2=0.1d0
      IF (INTERACT.EQ.0) THEN
          MAXCON=MAXCONAIM
        ELSE
          WRITE (*,*) ' ENTER number of contin. steps:       '
          READ (*,*) MAXCON
        ENDIF
      IF (INTERACT.LE.1) THEN
          INDEX=0
        ELSE 
      WRITE (*,*) ' ENTER 0 for contin. with respect to parameter,  '
      WRITE (*,*) ' ENTER K for fixing K-th component:   '
      READ (*,*) INDEX
        ENDIF
      IF (INDEX.GT.N .OR. INDEX.LT.1) INDEX=N+1
      IF (INTERACT.EQ.0.AND.ABS(STEP).GT.EPS) THEN
          S=STEP
        ELSE
          WRITE (*,*) ' ENTER STEPSIZE (real number):'
          READ (*,*) S
        ENDIF 
      IF (MAXCON.EQ.1) LEVEL=0
      IF (MAXCON.GT.1) THEN
         IF (INTERACT.LE.1) THEN
          LEVEL=2
          ELSE
         WRITE (*,*) ' ENTER 2 for variable stepsize, 0 for fixed steps'
          READ (*,*) LEVEL
          ENDIF
         END IF
      IF (LEVEL.GT.0) THEN
         IF (INTERACT.EQ.0) THEN
            STEPMA=ABS(WINDOW(2)-WINDOW(1))*0.1d0
           ELSE 
            WRITE (*,*) ' ENTER MAXIMUM STEPSIZE (REAL,  >0 )'
            READ (*,*) STEPMA
           ENDIF
         STEPMP=STEPMA
         END IF
C  TARGET FACILITIES ARE NOT EXPLOITED HERE:
      IFINAL=0
      IF (INTERACT.GT.0) THEN
        write (*,*) ' last No. branch /solution is:',nob,nos
        write (*,*) 'Enter 0 for continuation of numbering,'
        write (*,*) 'enter 1 for new numbers:'
        read (*,*) nosnob
        if (nosnob.ne.0) then
          write (*,*) ' Enter number of branch:'
          read (*,*) nob
          nos=0
          endif
        ENDIF
      IFLAGPER=1
      CALL CONTA (NOB,NOS,N,FCN,X,XA,LAST,EPS,LEVEL,INDEX,S,
     1     STEPMA,STEPMI,STEPMP,MAXCON,
     2     ITAIM,ITMAX,WEIGHT,IEXTRA,IFINAL,YFINAL,ISTAB,IRUN,
     3     IPRINT,WINDOW,BOUND2,Y0BUFF,NUMBUFF,Y0TEXT,
     4     HEAD1BUFF,STEPBUFF,JUDGEBUFF,UNITLST,UNITOUT,UNIT3,UNIT31)
      WRITE (*,601) INDEX,S
      GO TO 101

C************************************************************
8     continue
      CALL doubling (NOB,NOS,N,PAR,PERIOD,EPS,INTE,T,M,Y
     1   ,INDEX,UNITLST,UNITOUT,UNIT31)
      GO TO 101

991   WRITE (*,*) ' NO SUCCESS'
999   RETURN


601   FORMAT (' LAST STEPSIZE FOR K=',I2,' WAS',E10.3,10X)
602   FORMAT (//' ************  next run of BIFLAB  ( MODUS',
     1        I3,' )  ******************'/)
604   FORMAT (I2)
605   FORMAT ('  ',6E15.7)
606   FORMAT (E20.10)
607   FORMAT (' PARAMETER:',E19.11,' , PERIOD: ',E19.11)
608   FORMAT (2E20.10)
609   FORMAT (7X,'T (norm.)',8X,'Y1',13X,'Y2',13X,'Y3             ')
501   FORMAT (I1)
610   FORMAT (' THIS IS THE MENU YOU CAN CHOOSE FROM:              ',
     1       /' -------------------------------------              ',
     2       /'    0  :  PREPARATIONS (Jacobian check, interactions)',
     3       /'    1  :  CONTINUATION ALONG PERIODIC BRANCH        ',
     4       /'    2  :  HANDLING OF HOPF BIFURCATION              ',
     5       /'    3  :  DATA FOR A PHASE DIAGRAM / SIMULATION     ',
     6       /'    4  :  JUST (RE)CALCULATE ONE SOLUTION, CHECK STAB.',
     7       /'    6  :  check solutions in buffer      ',
     8       /'    7  :  STABILITY ANALYSIS OF STEADY-STATE BRANCH ',
     8       /'    8  :  calculate period doubling point (carry out',
     8         ' 4 first!)',
     9       /'   -1  :  TERMINATING THIS SESSION                  ')
9300  format (a1,i2,i4,2i2,a1,a2,13e20.13)
      END
