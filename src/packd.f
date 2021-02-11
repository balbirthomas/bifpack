C******************************************
C**                                      **
C**  this is  PACKD   version 12 August 1997 
C**  ---------------                     **
C**    INCLUDING ROUTINES FOR HANDLING   **
C**    BIFURCATIONS OF PERIODIC          **
C**    OSCILLATIONS OF O.D.E.            **
C**                                      **
C**    PACKd  is used most conveniently  **
C**    in connection with the driving    **
C**    procedure  B I F L A B   (packON) **
C**                                      **
C**      LATEST UPDATE:       **
C**                                      **
C******************************************


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

c=======================================================================


      SUBROUTINE HOPF (NOB,NOS,N,X,PAR,PERIOD,EPS,INTE,
     1  T,M,Y,IFLAGVEC,EVERE,EVEIM,MAXSOL,DECR,YOLD,JBIF,INDEX
     2  ,UNITLST,UNITOUT,UNIT31)

C    ALGORITHMS FOR CALCULATING BIFURCATIONS OF PERIODIC SOLUTIONS
C    OF O.D.E.  DY=F(Y,PAR) ,
C    CONSISTING OF THE SUBROUTINES:  HOPF, FABZ, RABZ, DOUBLING .
C
C    REFERENCES: R.SEYDEL, APPL.MATH.COMPUT.9 (1981),PP.257-271
C                (SEE ALSO REFERENCE IN SUBR. MONOD)

C  FIRST VERSION: DEZ 1981
C-------------------------------------------------------------
C
C    NEEDED SUBROUTINES:
C    ------------------
C         FCN(T,Y,F) : RIGHT-HAND SIDE OF DY=F(Y,PAR),
C                      T IS DUMMY PARAMETER.
C                      THE PARAMETER IS GIVEN BY Y(N+1)
C         A NUMERICAL PROCEDURE FOR SOLVING GENERAL NONLINEAR TWO-POINT
C                      BOUNDARY-VALUE PROBLEMS. SUCH A METHOD IS CALLED
C                      BY THE INTERFACE ROUTINE  S O L V E R .
C         INTE         an integrator

C--  SUBROUTINE, NEEDED IN CASE  JBIF.GE.0:
C         FCNLIN(T,Y,A) :  CALCULATES PARTIAL DERIV. OF RIGHT-HAND SIDE,
C                          T IS DUMMY PARAMETER,
C                      A(I,J)=DF(I)/DY(J)   1.LE. I,J .LE.N
C                      THE PARAMETER IS GIVEN BY Y(N+1)
C                      (BEFORE FIRST USE, CHECK FCNLIN BY ROUTINE CHECKD)

C*****************************************************************
C
C      IF THE USER ONLY WANTS TO ACCURATELY CALCULATE A HOPF POINT
C      WITHOUT BEING INTERESTED IN OSCILLATIONS HE MAY USE
C      A STEADY-STATE ROUTINE AS HOPFA (INCLUDED IN PACKA3)
C
C      THE PRESENT VERSION OF HOPF  IS NOT PARICULARLY EFFICIENT, 
C      MAINLY SINCE EIGENVECTORS ARE NOT EXPLOITED.
C
C*****************************************************************


C--  INPUT VARIABLES  FOR JBIF.LE.2  (HOPF BIFURCATION)
C    -------------------------------------------------
C         N      : NUMBER OF COMPONENTS OF  F(Y,PAR)
C         X      : VECTOR CONTAINING A STATIONARY SOLUTION NEAR HOPF
C                  BIFURCATION POINT
C         PAR    : VALUE OF CORRESPONDING PARAMETER
C         PERIOD : ROUGH APPROXIMATION OF INITIAL PERIOD
C
C         IFLAGVEC: =1 if eigenvector EVERE/EVEIM is provided, else =0
C         EVERE  : real part of eigenvector
C         EVEIM  : imaginary part of eigenvector
C
C         MAXSOL    : CONTROLS CALCULATION OF PERIODIC ORBITS.
C                     =0  :  NO PERIODIC ORBITS ARE CALCULATED.
C                     >0  :  MAXSOL  PERIODIC SOL. ARE TO BE CALCULATED.
C                            (SMALL-AMPLITUDE PERIODIC OSCILLATIONS
C                            NEAR THE HOPF BIFURCATION)
C         EPS= RELATIVE ERROR TOLERANCE
C
C         DECR   : real number controlling the "smallness" in the
C                  amplitudes and their differences.
C                  in the range, say,   0.01 to 0.05
C
C         INDEX  : INDEX OF COMPONENT WHOSE INITIAL VALUE AND SLOPE
C                  ARE PRESCRIBED FOR NORMALIZATION PURPOSES (SAME AS INDK).
C                  =0  : INDEX IS DETERMINED BY THIS ROUTINE.
C
C         JBIF   : OPTION OF CALCULATING THE HOPF BIFURCATION POINT.
C                  =2  : THE BRANCH POINT AND THE LINEARIZATION (H IN REF.)
C                        ARE APPROXIMATED BY SOLVING THE BRANCHING SYSTEM.
C                        SMALL-AMPLITUDE OSCILLATIONS CAN BE CALCULATED.
C                        SUBROUTINE  F C N L I N  IS NEEDED.
C                        THE VALUES OF X, PAR, PERIOD  ARE UPDATED BY
C                        THE HOPF BIFURCATION POINT.
C                        DIMENSION  Y0(2*N+2,M)
C
C                  =1  : THE BRANCH POINT IS NOT APPROXIMATED, A GUESS OF THE
C                        SOLUTION  H  TO THE LINEARIZATION IS CALCULATED.
C                        BASED ON THIS GUESS, SMALL-AMPLITUDE PERIODIC
C                        OSCILLATIONS ARE APPROXIMATED, PARAMETRIZED BY
C                        THE INITIAL VALUE OF THE K-TH COMPONENT.
C                        SUBROUTINE  F C N L I N  IS NEEDED.
C                        THE COMPONENT H(INDK)=Y(N+2+INDK,.) (USUALLY INDK=1)
C                        IS CORRECTLY APPROXIMATED (THE SOLUTION IS JUST THE
C                        COSINE-FUNCTION). IN CASE N=2, ALSO THE OTHER
C                        COMPONENT IS CALCULATED.
C                        DIMENSION Y0(2*N+2,M) , HOWEVER SOLVER NEEDS ONLY
C                                                N+2  VARIABLES
C
C                  =0  : SAME AS JBIF=1 , EXCEPT THAT NOW THE LOCAL
C                        OSCILLATIONS ARE PARAMETRIZED BY AMPLITUDE.
C                        SEE ADDITIONAL REMARKS ON IMPLEMENTATION BELOW.

C                 =-1  : NO APPROXIMATION OF HOPF BIFURCATION POINT,
C                        NO CALCULATION OF  LINEARIZATION.
C                        SUBROUTINE  FCNLIN  IS NOT NEEDED.
C                        IN THIS CASE, SUBROUTINE HOPF PRESCRIBES A NONZERO
C                        AMPLITUDE AND TRIES TO JUMP ONTO THE PERIODIC
C                        BRANCH WITHOUT EXPLOITING ANY INFORMATION.
C                        (ONE SHOULD NOT BE TOO OPTIMISTIC WITH JBIF=-1)
C
C         =-1 AND =0, REMARK ON THE IMPLEMENTATION:
C                        THESE SPECIAL CASES ARE TAILORED TO BE USED ALONG
C                        WITH MULTIPLE SHOOTING CODE ODEBC1, which needs
C                        THE STATEMENT:     COMMON /INTVAL/ J
C                        DIMENSION NEEDED FOR ODEBC1:  Y0(N+3,M)
C
C
C--  OUTPUT VARIABLES (FOR JBIF.LE.2) 
C    ----------------
C         X, PAR, PERIOD  :  SEE INPUT, JBIF=2 .
C         Y      : CONTAINS FINALLY CALCULATED PERIODIC SOLUTION.
C                  Y(.,J) CONTAINS VECTOR AT J=1,..,M EQUIDISTANT NODES.
C                  Y(N+1,.) CONTAINS CORRESPONDING PARAMETER
C                  Y(N+2,.) CONTAINS CORRESPONDING PERIOD
C                  T        VECTOR OF NODES AT WHICH APPROX.Y IS CALCULATED
C                  M        NUMBER OF NODES IN T
C         YOLD:  a previous solution
C         MAXSOL=-1  :  NOT all PERIODIC SOLUTIONS were FOUND.
C               < input value: not all solutions were calculated

      DOUBLE PRECISION T,X,Y,Y0,A,YA,CC,DEKR,DELTA,DELTAR,ETACO,ZETACO,
     1   PARCO,EPS,EPB,H2,PAR,PERIOD,PI2,SS,T1,YOLD,DECR
     2   ,EVERE,EVEIM,BETA,DELTA2 
      INTEGER UNITLST,UNITOUT,UNIT31
      CHARACTER*1 MODEL
      PARAMETER (NDIM=12,NDIML=14)
      DIMENSION T(7),X(NDIML),Y(NDIML,7),Y0(NDIML,7), A(NDIML,NDIML),
     1  YA(NDIML),YOLD(NDIML,7),EVERE(NDIML),EVEIM(NDIML),YN(NDIML,7)
C  dimension at least : T(M),X(N),Y(N+3,M),Y0(N+3,M) OR Y0(2*N+2,M)
C                       A(N,N),YA(N+1)

      EXTERNAL FABZ,RABZ,INTE
      COMMON /CONTRO/ ETACO,ZETACO,PARCO, NCO,N1CO,N2CO,
     1                KCO,IND2CO,NSECFLAG,NFLAG

      WRITE (*,*) 'BIFPACK: entering routine HOPF.'
      WRITE (UNITLST,*) 'BIFPACK: entering routine HOPF.'
C--  TO BE CHANGED IF NECESSARY:
      IF (INDEX.GT.N) INDEX=0
      IF (JBIF.LE.2) M=7
      MODEL='d'
C     M=NUMBER OF EQUIDISTANT NODES FOR MULTIPLE SHOOTING.
C     EQUIDISTANT NODES ARE SUFFICIENT IN THIS CONTEXT, AS
C     SMALL-AMPLITUDE SOLUTIONS ARE CALCULATED.
      PI2=6.2831853071796d0
      IF (ABS(PERIOD).GT.EPS) THEN
	 BETA=PI2/PERIOD
	else
	 write (*,*) 'your guess for period is: ',period
	 write (*,*) 'Enter 1 for continuing with 2*pi ,'
	 write (*,*) 'enter 0 for returning:'
	 read (*,*) ind2
	 if (ind2.le.0) then
	     maxsol=-1
	     return
	     endif
	 beta=1.d0
	endif
      ETACO=DECR
      ITMAX=20
      ITAIM=7

C*************************************************
C--  INITIAL PREPARATIONS

      NCO=N
      N2=N+2
      N2CO=N2
      KCO=INDEX
      IF (ABS(X(1)).GT.1.d0) ETACO=ETACO*ABS(X(1))
      DEKR=ETACO
      DELTA=ETACO
      N3=N+3
      T(1)=0.d0
      T1=1./DBLE(M-1)
      DO 19 J=3,M
19    T(J-1)=T(J-2)+T1
      T(M)=1.d0

      IF (KCO.EQ.0) THEN
          KCO=1
c          CC=0.d0
c          DO 12 I=1,N
c          IF (ABS(X(I)).GT.CC) THEN
c             CC=ABS(X(I))
c             KCO=I
c             END IF
c12        CONTINUE
          INDEX=KCO
          END IF

      DO 1 J=1,M
      DO 11 I=1,N
      Y0(I,J)=X(I)
11    Y(I,J)=X(I)
      Y(N+1,J)=PAR
      Y0(N+1,J)=PAR
      Y(N2 ,J)=PERIOD
1     Y0(N2,J)=PERIOD
      IF (JBIF.LT.0) GO TO 4


C**************************************
C--  APPROXIMATION OF BIFURCATION POINT AND LINEARIZATION

      NL=N+N+2
      DO 3 I=1,N
      YA(I)=Y0(I,1)
      DO 3 J=1,M
3     Y0(N2+I,J)=0.d0
      YA(N+1)=PAR
      IF (IFLAGVEC.EQ.1) THEN
         DO 33 J=1,M
         SS=SIN(PI2*T(J))
         CC=COS(PI2*T(J))
         DO 33 I=1,N
33       Y0(N2+I,J)=CC*EVERE(I)-SS*EVEIM(I)
         GO TO 302
         ENDIF

      IFLAG=0
      GO TO 303
301   IF (IFLAG.EQ.0) THEN
          IFLAG=1
          KCO=0
          END IF
      KCO=KCO+1

C--  INITIAL GUESS FOR BRANCHING SYSTEM
303   CALL FCNLIN (T1,YA,A)
C     T1 is dummy here
      IND2=0
305   IND2=IND2+1
      IF (IND2.EQ.KCO) GO TO 305
      IF (ABS(A(KCO,IND2)).LT.EPS) GO TO 305
      IF (IND2.GT.N) GO TO 301
      DO 31 J=1,M
      SS=SIN(PI2*T(J))
      CC=COS(PI2*T(J))
      H2=(-BETA*SS-A(KCO,KCO)*CC)/A(KCO,IND2)
      Y0(N2+KCO,J)=CC
31    Y0(N2+IND2,J)=H2
302   IF (JBIF.LT.2) GO TO 34

C--  SOLUTION OF THE BRANCHING SYSTEM, HOPF CASE (DIMENSION 2*N+2)
      NFLAG=-2
      J=-1
      EPB=EPS
      IF (EPS.LT.0.005d0) EPB=0.005d0
C     CONVERGENCE CAN BE EXPECTED TO BE FAST FOR THE PARAMETER AND THE 
C           PERIOD (SUPERCONVERGENCE), BUT SLOWER FOR THE COMPONENTS OF THE 
C           LINEARIZATION. ACCURACIES WILL BEHAVE ACCORDINGLY. HENCE THE 
C           RELATIVELY LARGE VALUE OF ERROR TOLERANCE  EPB .
      IT =10
      CALL SOLVER (FABZ,RABZ,INTE,NL,M,T,Y0,EPB,IT)
      IF (IT.EQ.-1) GO TO 34
      IF (IT.LT.0) GO TO 41
      PAR=Y0(N+1,1)
      PERIOD=Y0(N2,1)
      DO 32 I=1,N
32    X(I)=Y0(I,1)
      NOS=1
      IF (NSECFLAG.EQ.0) WRITE (UNITOUT,9300) 
     1    'R',nob,nos,n,1,MODEL,'HB',(x(i),i=1,n),par,period
      IF (NSECFLAG.EQ.1) WRITE (UNITOUT,9300) 
     1    'R',nob,nos,n,2,MODEL,'HB',(x(i),i=1,n),par,PARCO,period
      WRITE (UNITLST,602)
      WRITE (UNITLST,601) IT,Y0(N+1,1),Y0(N2,1)
      WRITE (*,602)
      WRITE (*,601) IT,Y0(N+1,1),Y0(N2,1)
      WRITE (UNITLST,6051)
      WRITE (UNITLST,605) (Y0(I,1),I=1,N)
      REWIND UNIT31
      WRITE (UNIT31,*) NOB
      WRITE (UNIT31,*) NOS
      WRITE (UNIT31,*) 'HB'
      WRITE (UNIT31,*) 2
      WRITE (UNIT31,*) 0.


C--  INITIAL GUESS FOR small-amplitude PERIODIC SOLUTION
34    DELTAR=Y0(N2+KCO,1)-Y0(N2+KCO,2)
      IF (ABS(DELTAR).LT.1.d-7) GO TO 4
      DELTA=ETACO/DELTAR
      DELTA2=DELTA
      DO 36 J=1,M
      Y(N+1,J)=Y0(N+1,J)
      Y(N2,J)=Y0(N2,J)
      DO 36 I=1,N
36    Y(I,J)=Y0(I,J)+DELTA*Y0(I+N2,J)
      GO TO 41

C**************************************
C     COMPUTATION OF PERIODIC SOLUTIONS (LOCALLY)

4     Y(KCO,1)=Y(KCO,1)+DEKR*0.5d0
      Y(KCO,2)=Y(KCO,2)-DEKR*0.5d0
41    IF (MAXSOL.LE.0) RETURN
      MIST=0
      WRITE (UNITLST,604) KCO
      IF (JBIF.LE.0) GO TO 43

C     PARAMETRIZATION BY INITIAL VALUE:
      NFLAG=2
      N3=N2
      ETACO=Y(KCO,1)
      GO TO 5

43    Y(N3,1)=Y(KCO,1)
C     PARAMETRIZATION BY AMPLITUDE
      Y(N3,2)=Y(KCO,2)
      DO 44 J=3,M
44    Y(N3,J)=Y(N3,J-1)
      NFLAG=0

5     IFLAG=-1
      DO 59 ILOOP=1,MAXSOL
      IT=ITMAX
      IF (ILOOP.GT.2) ALPHA=DEKR/(YN(KCO,1)-YOLD(KCO,1))
      DO 55 J=1,M
      DO 55 I=1,N3
      IF (ILOOP.EQ.2.AND.I.LE.N) Y(I,J)=Y0(I,J)+DELTA*2.*Y0(I+N2,J)
      IF (ILOOP.GT.2) Y(I,J)=(1.+ALPHA)*YN(I,J)-ALPHA*YOLD(I,J)
55    CONTINUE
      ETACO=Y(KCO,1)
      CALL SOLVER (FABZ,RABZ,INTE,N3,M,T,Y,EPS,IT)
      IF (IT.GT.0) IFLAG=1
      IF (IT.GT.0) MIST=-1
      IF (IT.LT.0) WRITE (UNITLST,603)
      WRITE (UNITLST,601) IT,Y(N+1,1),Y(N2,1)
      IF (IT.LT.0) WRITE (*,603)
      WRITE (*,601)  IT,Y(N+1,1),Y(N2,1)
      WRITE (UNITLST,6051)
      WRITE (UNITLST,605) (Y(I,1),I=1,N)
      IF (IT.LT.0.AND.MIST.LT.0) GO TO 591
      IF (IT.LT.0) MIST=MIST+1
      IF (MIST.GT.2) GO TO 591

C     ADJUSTMENT OF STEPSIZE AFTER THE second STEP:
      IF (ILOOP.GT.2 .AND. NFLAG.NE.2)
     1                DEKR=DEKR*DBLE(ITAIM)/DBLE(IT)
      IF (ILOOP.GT.2 .AND. NFLAG.EQ.2)
     1                DELTA=DELTA*DBLE(ITAIM)/DBLE(IT)
      ETACO=ETACO+DEKR
      IF (NFLAG.EQ.2) ETACO=ETACO+DELTA
      IF (IT.LT.0) GO TO 59
      NOS=NOS+1
      IF (NSECFLAG.EQ.0) WRITE (UNITOUT,9300) 
     1  'R',nob,nos,n,1,MODEL,'P ',(y(i,1),i=1,n+1),y(n2,1)
      IF (NSECFLAG.EQ.1) WRITE (UNITOUT,9300) 
     1  'R',nob,nos,n,2,MODEL,'P ',(y(i,1),i=1,n+1),PARCO,y(n2,1)
      REWIND UNIT31
      WRITE (UNIT31,*) NOB
      WRITE (UNIT31,*) NOS
      WRITE (UNIT31,*) 'P '
      WRITE (UNIT31,*) 1
      WRITE (UNIT31,*) Y(N+1,1)-YN(N+1,1)
      DO 58 I=1,N3
      DO 58 J=1,M
      YOLD(i,j)=YN(I,J)
      YN(i,j)=Y(I,J)
58    CONTINUE
59    CONTINUE
      MAXSOL=ILOOP
      RETURN
591   IF (IFLAG.LT.0) MAXSOL=-1
      MAXSOL=ILOOP
      RETURN



601   FORMAT (' After ',I2,' iter. PARAMETER:',E19.11,
     1        '   PERIOD:',E19.11)
602   FORMAT (' **** BIFURCATION OBTAINED: ****')
603   FORMAT (' NO CONVERGENCE, TERMINATES WITH:')
604   FORMAT (' SUBROUTINE HOPF ATTEMPTS TO CALCULATE PERIODIC SOL.S'/,
     1            ' FIXED/NORMALIZED COMPONENT:',I3)
605   FORMAT ('   ',4E19.11)
6051  FORMAT (' initial vector:')
9300  format (a1,i2,i4,2i2,a1,a2,18e20.13)
      END

      SUBROUTINE FABZ (T,Y,DY)       
      DOUBLE PRECISION Y,DY,T,ETACO,ZETACO,PARCO,  A,TIME,SP
      PARAMETER (NDIM=12,NDIML=14)
      DIMENSION Y(NDIML),DY(NDIML),A(NDIML,NDIML)
C   DIMENSIONS AT LEAST:  A(N,N), Y(.) AND DY(.) DEPEND ON OPTION:
C                                 N+2  OR  N+3  OR  2*N+2
      COMMON /INTVAL/ J
      COMMON /CONTRO/ ETACO,ZETACO,PARCO, NCO,N1CO,N2CO,
     1                KCO,IND2CO,NSECFLAG,NFLAG

      IF (NFLAG.EQ.-1) Y(N2CO)=ETACO
      CALL FCN (T,Y,DY)
      IF (NFLAG.EQ.-1) RETURN
      TIME=Y(N2CO)

      DO 1 I=1,NCO
1     DY(I)=DY(I)*TIME
      DY(NCO+1)=0.d0
      DY(N2CO)=0.d0
      IF (NFLAG.EQ.-2 .OR. NFLAG.GE.3) GO TO 3
      IF (NFLAG.NE.0) RETURN
      DY(NCO+3)=DY(KCO)
      IF (J.GE.2) DY(NCO+3)=0.d0
      RETURN

3     CALL FCNLIN (T,Y,A)
      DO 5 I=1,NCO
      I1=N2CO+I
      SP=0.d0
      DO 4 K=1,NCO
4     SP=SP+A(I,K)*Y(N2CO+K)
5     DY(I1)=TIME*SP
      IF (NFLAG.EQ.-2) RETURN
      DO 6 I=1,NCO
6     DY(N2CO+I)=DY(N2CO+I)+Y(N2CO+NCO+1)*DY(I)/TIME
      DY(N2CO+NCO+1)=0.d0
      RETURN
      END

      SUBROUTINE RABZ (YA,YB,W)
      DOUBLE PRECISION YA,YB,W,A,SP,T,ETACO,ZETACO,PARCO
      PARAMETER (NDIM=12,NDIML=14)
      DIMENSION YA(NDIML),YB(NDIML),W(NDIML),A(NDIML,NDIML)
C    DIMENSIONS AT LEAST:  A(N,N), OTHER VARIABLES DEPEND ON OPTION
      COMMON /CONTRO/ ETACO,ZETACO,PARCO, NCO,N1CO,N2CO,
     1                KCO,IND2CO,NSECFLAG,NFLAG
      IF (NFLAG.EQ.-2) GO TO 3
      CALL FCN (T,YA,W)
      IF (NFLAG.GE.3) GO TO 2
      W(NCO+1)=W(KCO)
      DO 1 I=1,NCO
1     W(I)=YA(I)-YB(I)

      IF (NFLAG.EQ.1) W(N2CO)=YA(NCO+1)-ETACO
      IF (NFLAG.EQ.2) W(N2CO)=YA(KCO)-ETACO
      IF (NFLAG.NE.0) RETURN

      W(N2CO)=YA(KCO)-YB(NCO+3)-ETACO
      W(NCO+3)=YA(KCO)-YA(NCO+3)
      RETURN

2     W(N2CO+NCO+1)=W(KCO)
3     CALL FCNLIN (T,YA,A)
      SP=0.d0
      DO 4 I=1,NCO
      W(I)=YA(I)-YB(I)
      I1=N2CO+I
      W(I1)=YA(I1)-YB(I1)
4     SP=SP+A(KCO,I)*YA(I1)
      W(NCO+1)=SP
      W(N2CO)=YA(N2CO+KCO)-1.d0
      RETURN
      END


      SUBROUTINE DOUBLING (NOB,NOS,N,PAR,PERIOD,EPS,INTE,
     1  T,M,Y,INDEX,UNITLST,UNITOUT,UNIT31)

C    ALGORITHMS FOR CALCULATING BIFURCATIONS OF PERIODIC SOLUTIONS
C    OF O.D.E.  DY=F(Y,PAR) ,
C    preliminary version for period doubling
C    REFERENCES: R.SEYDEL, APPL.MATH.COMPUT.9 (1981),PP.257-271,
C       with symmetry-reduction to simple interval length.

C-------------------------------------------------------------
C
C    NEEDED SUBROUTINES:
C    ------------------
C         FCN(T,Y,F) : RIGHT-HAND SIDE OF DY=F(Y,PAR),
C                      T IS DUMMY PARAMETER.
C                      THE PARAMETER IS GIVEN BY Y(N+1)
C         A NUMERICAL PROCEDURE FOR SOLVING GENERAL NONLINEAR TWO-POINT
C                      BOUNDARY-VALUE PROBLEMS. SUCH A METHOD IS CALLED
C                      BY THE INTERFACE ROUTINE  S O L V E R .
C         INTE         an integrator
C
C         FCNLIN(T,Y,A) :  CALCULATES PARTIAL DERIV. OF RIGHT-HAND SIDE,
C                          T IS DUMMY PARAMETER,
C                      A(I,J)=DF(I)/DY(J)   1.LE. I,J .LE.N
C                      THE PARAMETER IS GIVEN BY Y(N+1)
C                      (BEFORE FIRST USE, CHECK FCNLIN BY ROUTINE CHECKD)


      DOUBLE PRECISION T,Y,Y0,ETACO,ZETACO,
     1   PARCO,EPS,EPB,PAR,PERIOD 
      INTEGER UNITLST,UNITOUT,UNIT31
      CHARACTER*1 MODEL
      PARAMETER (NDIM=12,NDIML=14)
      DIMENSION T(7),Y(NDIML,7),Y0(NDIML,7)
C  dimensions:  T(M),Y(N+3,M), Y0(2*N+2,M)

      EXTERNAL FABZPD,RABZPD,INTE
      COMMON /CONTRO/ ETACO,ZETACO,PARCO, NCO,N1CO,N2CO,
     1                KCO,IND2CO,NSECFLAG,NFLAG
      WRITE (*,*) 'BIFPACK: entering routine DOUBLING.'
      WRITE (UNITLST,*) 'BIFPACK: entering routine DOUBLING.'
C--  TO BE CHANGED IF NECESSARY:
      IF (INDEX.GT.N) INDEX=0
      MODEL='d'
      ITMAX=20
      ITAIM=7

C*************************************************
C--  INITIAL PREPARATIONS

      NCO=N
c      write (*,*) 'initial prep...'
       N2=N+2
      N2CO=N2
      KCO=INDEX

      IF (KCO.EQ.0) THEN
          KCO=1
          INDEX=KCO
          END IF

      DO 11 J=1,M
      DO 11 I=1,N+2
      Y0(I,J)=Y(I,J)
11    continue


C**************************************
C--  APPROXIMATION OF BIFURCATION POINT AND LINEARIZATION

c      write (*,*) 'approxim...'
      NL=N+N+2
      DO 3 I=1,N
      DO 3 J=1,M
3     Y0(N2+I,J)=0.1d0

      IFLAG=0
      GO TO 303
301   IF (IFLAG.EQ.0) THEN
          IFLAG=1
          KCO=0
          END IF
      KCO=KCO+1

C--  INITIAL GUESS FOR BRANCHING SYSTEM
303   continue
      write (*,*) 'kco=',kco
c      write (*,*) 'calculating initial H:'
      NFLAG=-2
      etaco=y0(n+1,1)
      y0(n+2+kco,1)=1.d0
      IT =10
      EPB=1.d-3
c      write (*,*) 'initial data: '
c      write (*,*) (y0(i,1),i=1,nl)
c      write (*,*) (y0(i,m),i=1,nl)
      CALL SOLVER (FABZPD,RABZPD,INTE,NL,M,T,Y0,EPB,IT)
      if (it.gt.0) write (*,*) 'successful calc. of h, it=',it
      if (it.lt.0) write (*,*) 'NO successful calc. of h'
c      write (*,*) 'final data: '
c      write (*,*) (y0(i,1),i=1,nl)
c      write (*,*) (y0(i,m),i=1,nl)
      EPB=EPS
      y0(n+2+kco,1)=1.d0
c      IF (EPS.LT.0.005d0) EPB=0.005d0
      write (*,*) 'solving branching system:'
C     CONVERGENCE CAN BE EXPECTED TO BE FAST FOR THE PARAMETER AND THE 
C           PERIOD (SUPERCONVERGENCE), BUT SLOWER FOR THE COMPONENTS OF THE 
C           LINEARIZATION. ACCURACIES WILL BEHAVE ACCORDINGLY. HENCE THE 
C           RELATIVELY LARGE VALUE OF ERROR TOLERANCE  EPB .
      IT =10
      NFLAG=2
      CALL SOLVER (FABZPD,RABZPD,INTE,NL,M,T,Y0,EPB,IT)
      if (it.gt.0) write (*,*) 'sol. of branching system, it=',it
      if (it.lt.0) write (*,*) 'NO success'
c      write (*,*) 'final data: '
c      write (*,*) (y0(i,1),i=1,nl)
c      write (*,*) (y0(i,m),i=1,nl)
      IF (IT.LE.0) RETURN
      PAR=Y0(N+1,1)
      PERIOD=Y0(N2,1)
      NOS=1
      IF (NSECFLAG.EQ.0) WRITE (UNITOUT,9300) 
     1    'R',nob,nos,n,1,MODEL,'PD',(y0(i,1),i=1,n),par,period
      IF (NSECFLAG.EQ.1) WRITE (UNITOUT,9300) 
     1    'R',nob,nos,n,2,MODEL,'PD',(y0(i,1),i=1,n),par,PARCO,period
      WRITE (UNITLST,602) epb
      WRITE (UNITLST,601) IT,Y0(N+1,1),Y0(N2,1)
      WRITE (*,602) epb
      WRITE (*,601) IT,Y0(N+1,1),Y0(N2,1)
      WRITE (UNITLST,6051)
      WRITE (UNITLST,605) (Y0(I,1),I=1,N)
      REWIND UNIT31
      WRITE (UNIT31,*) NOB
      WRITE (UNIT31,*) NOS
      WRITE (UNIT31,*) 'PD'
      WRITE (UNIT31,*) 2
      WRITE (UNIT31,*) 0.

      RETURN

601   FORMAT (' After ',I2,' iter. PARAMETER:',E19.11,
     1        '   PERIOD:',E19.11)
602   FORMAT (' **** PERIOD DOUBLING bifurcation obtained ****',
     1         '  eps=',e10.3)
603   FORMAT (' NO CONVERGENCE, TERMINATES WITH:')
604   FORMAT (' SUBROUTINE HOPF ATTEMPTS TO CALCULATE PERIODIC SOL.S'/,
     1            ' FIXED/NORMALIZED COMPONENT:',I3)
605   FORMAT ('   ',4E19.11)
6051  FORMAT (' initial vector:')
9300  format (a1,i2,i4,2i2,a1,a2,18e20.13)
      END


      SUBROUTINE FABZPD (T,Y,DY)       
      DOUBLE PRECISION Y,DY,T,ETACO,ZETACO,PARCO,  A,TIME,SP
      PARAMETER (NDIM=12,NDIML=14)
      DIMENSION Y(NDIML),DY(NDIML),A(NDIML,NDIML)
C   DIMENSIONS AT LEAST:  A(N,N), Y(.) AND DY(.) DEPEND ON OPTION:
C                                 N+2   OR  2*N+2
      COMMON /INTVAL/ J
      COMMON /CONTRO/ ETACO,ZETACO,PARCO, NCO,N1CO,N2CO,
     1                KCO,IND2CO,NSECFLAG,NFLAG

      CALL FCN (T,Y,DY)
      TIME=Y(N2CO)

      DO 1 I=1,NCO
1     DY(I)=DY(I)*TIME
      DY(NCO+1)=0.d0
      DY(N2CO)=0.d0

      CALL FCNLIN (T,Y,A)
      DO 5 I=1,NCO
      I1=N2CO+I
      SP=0.d0
      DO 4 K=1,NCO
4     SP=SP+A(I,K)*Y(N2CO+K)
5     DY(I1)=TIME*SP
      RETURN
      END


      SUBROUTINE RABZPD (YA,YB,W)
      DOUBLE PRECISION YA,YB,W,T,ETACO,ZETACO,PARCO
      PARAMETER (NDIM=12,NDIML=14)
      DIMENSION YA(NDIML),YB(NDIML),W(NDIML)
      COMMON /CONTRO/ ETACO,ZETACO,PARCO, NCO,N1CO,N2CO,
     1                KCO,IND2CO,NSECFLAG,NFLAG
      CALL FCN (T,YA,W)
      W(NCO+1)=W(KCO)
      DO 4 I=1,NCO
      W(I)=YA(I)-YB(I)
      W(N2CO+I)=YA(N2CO+I)+YB(N2CO+I)
4     CONTINUE
      IF (NFLAG.EQ.-2) W(N2CO+1)=YA(NCO+1)-ETACO
      W(N2CO)=YA(N2CO+KCO)-1.d0
      RETURN
      END


      SUBROUTINE MONOD  (N0,Y,EPS,A,INTE,UNITLST)
      DOUBLE PRECISION EPS,Y,Y1,A,A1,HMAX,T0,T1,ETACO,ZETACO,PARCO,HI
      INTEGER UNITLST
      PARAMETER (NDIM=12,NDIML=14)
      DIMENSION Y(NDIML,7),Y1(NDIML),A(NDIML,NDIML)
C    DIMENSIONS AT LEAST:  A(N,N),Y(N+N+2,M),Y1(N+N+2)
      COMMON /CONTRO/ ETACO,ZETACO,PARCO, NCO,N1CO,N2CO,
     1                KCO,IND2CO,NSECFALG,NFLAG
      EXTERNAL FABZ,INTE

C*************************************************************
C  THIS ROUTINE CALCULATES THE MONODROMY MATRIX A
C  OF A PERIODIC SOLUTION PROVIDED BY  Y
C*************************************************************
C  METHOD: Matrix A is calculated columnwise. Each of these N=N0
C  ------  vectors is obtained by integrating an initial-value
C          problem of dimension 2N+2.
C          The initial values of the first n components are
C          those of the periodic solution, the next two are the
C          parameter and period, the final N initial values are
c          the K-th unit vector (in the K-th integration,
C          K=1,2,..,N).
C          The final values of the latter components are the
c          desired column vectors of the monodromy matrix.
C
C  to be used along with subroutine FABZ of  H O P F.
C  subroutines supplied by the user:  F C N   and   F C N L I N
C**************************************************************

      IPRINT=0
C     WITH  IPRINT=1  THE MONODROMY MATRIX WILL BE PRINTED OUT.

      NCO=N0
      N2=NCO+2
      N2CO=N2
      NL=N2+NCO
      HI=0.01d0
      NFLAG=-2
      IF (IPRINT.GT.0) WRITE (UNITLST,604) Y(NCO+1,1),Y(N2,1)

      DO 8 K=1,NCO
      DO 2 I=1,NCO
      Y1(I)=Y(I,1)
2     Y1(I+N2)=0.d0
      Y1(NCO+1)=Y(NCO+1,1)
      Y1(N2)=Y(N2,1)
      Y1(N2+K)=1.d0
      T0=0.d0
      T1=1.d0
      HMAX=0.1d0

      CALL INTE (NL,FABZ,T0,Y1,T1,EPS,HMAX,HI)
      A1=ABS(T0-T1)
      DO 3 I=1,NCO
3     A1=A1+ABS(Y1(I)-Y(I,1))
      IF (A1.GT.10.*EPS*NCO) WRITE (UNITLST,601)
      DO 4 I=1,NCO
4     A(I,K)=Y1(I+N2)
8     CONTINUE

      DO 81 I=1,NCO
81    IF (IPRINT.GT.0) WRITE (UNITLST,602) (A(I,J),J=1,NCO)
      RETURN

601   FORMAT (' PERIODIC SOLUTION NOT CONFIRMED')
602   FORMAT (' MONOD.M. ',4E15.8)
604   FORMAT (' MONODROMY MATRIX AT PAR.:',E15.7,' , PERIOD:',E15.7)
      END


      SUBROUTINE MONOD2 (N,M,E)
      DOUBLE PRECISION E,HA,G,A
      PARAMETER (NDIM=12,NDIML=14)
      DIMENSION E(NDIML,NDIML),HA(NDIML),G(NDIML,NDIML,6),
     1          A(NDIML,NDIML)
C    DIMENSIONS AS IN  BRODE  OR  ODEBC1

      COMMON /BIF/ G,A

C            M       : NUMBER OF NODES OF SHOOTING
C  THIS ROUTINE CALCULATES THE MONODROMY MATRIX VIA THE ITERATION
C  MATRICES OF MULTIPLE SHOOTING CODE ODEBC1.
C     REF.: SAME AS IN MONOD.

C    REMARKS:
C    -------  THIS VERSION  USES THE DATA OBTAINED BY ODEBC1
C             (ARRAYS G, A). THESE DATA ARE USUALLY UPDATED BY
C             RANK-1 APPROXIMATIONS LEADING TO INACCURATE VALUES
C             OF THE MONODROMY MATRIX.
C             FOR MORE ACCURATE VALUES, PERFORM A SECOND (DOUBLE)
C             CALL OF ODEBC1.
C----------------------------------------------------------------------------

      M1=M-1

C--  CALCULATION OF MATRIX  E
      DO 1 I=1,N
      DO 1 J=1,N
1     E(I,J)=G(I,J,1)
      IF (M.EQ.2) GO TO 21
      DO 2 K=2,M1
      DO 2 II=1,N
      DO 15 I=1,N
      HA(I)=0.d0
      DO 15 J=1,N
15    HA(I)=HA(I)+G(I,J,K)*E(J,II)
      DO 2 I=1,N
2     E(I,II)=HA(I)
21    DO 22 I=1,N
      DO 22 J=1,N
22    E(I,J)=A(I,J)+E(I,J)

      DO 3 I=1,N
3     E(I,I)=E(I,I)-1.
      DO 31 I=1,N
      DO 31 J=1,N
31    E(I,J)=-E(I,J)
      RETURN

      END


      SUBROUTINE STAPER (NOB,N,M,Y,EPS,YOLD,EVARE,EVAIM,IFLAG,
     1           IND,JACOBI,INTE,EPSONE,Y0BUFF,NUMBUFF,Y0TEXT,
     2           HEAD1BUFF,STEPBUFF,JUDGEBUFF,UNITLST,UNITOUT)

C     THIS ROUTINE INVESTIGATES THE STABILITY OF BRANCHES OF PERIODIC
C     SOLUTIONS BY CHECKING THE EIGENVALUES OF THE MONODROMY MATRIX.
C     THE CURRENT EIGENVALUES ARE COMPARED TO THE EIGENVALUES OF
C     THE PRECEDING SOLUTION / MONODROMY MATRIX.
C
C
C     INPUT VALUES:
C     ------------    N, M, Y, EPS     as usual,
C                            IFLAG  only on first call of STAPER,
C                            YOLD  is the previous solution.
C                     JACOBI  =1 if jacobian is implemented in
C                             subroutine FCNLIN,  else  =0  .
      DOUBLE PRECISION Y,A,EPS,EVARE,EVAIM,PARAM,B,EIGRE,EIGIM,
     1   ADIST,AIM,AIMA,AMAX,APER,ARE,AMAXOL,EPSONE,EPSTOL,PAR0,FAK
     2   ,YOLD,YINIT0,Y0BUFF,ETACO,ZETACO,PARCO,STEPBUFF
      CHARACTER*9 Y0TEXT,TEXT
      CHARACTER*2 TEXT2,JUDGEBUFF
      CHARACTER*1 MODEL
      INTEGER IDUMMY,HEAD1BUFF,UNITLST,UNITOUT
      PARAMETER (NDIM=12,NDIML=14,NDIMBUFF=6)
      DIMENSION Y(NDIML,7),YOLD(NDIML,7),A(NDIML,NDIML),IDUMMY(NDIML),
     1          YINIT0(NDIML),Y0BUFF(NDIML,NDIMBUFF),Y0TEXT(NDIMBUFF),
     2          HEAD1BUFF(3,NDIMBUFF),STEPBUFF(NDIMBUFF),
     3          JUDGEBUFF(NDIMBUFF)
C    DIMENSIONS AT LEAST  A(N,N),Y(N+2,M)

C  FOR EIGENVALUE CALCULATION THE DIMENSIONS ARE:
C     B(NROW,N),EIGRE(N),EIGIM(N)
      DIMENSION B(NDIML,NDIM),EIGRE(NDIM),EIGIM(NDIM)
C   THE row-DIMENSION OF B HAS TO BE SPECIFIED BELOW BY  NROW=...

      COMMON /CONTRO/ ETACO,ZETACO,PARCO, NCO,N1CO,N2CO,
     1                KCO,IND2CO,NSECFLAG,NFLAG
      EXTERNAL INTE

C IFLAG IS NEEDED TO COMPARE TWO CONSECUTIVE RUNS OF STAPER.
C -----
C   INPUT: 
C     =0     ONLY ON FIRST CALL OF STAPER !
C   OUTPUT: (should not be changed; used as input on next call)
C     .LT.0  IF PERIODIC SOLUTION IS STABLE
C     .GT.0  IF PERIODIC SOLUTION IS UNSTABLE
C     THE ABSOLUTE VALUE OF IFLAG MAY GIVE DETAILS:
C            =1 : EIGENVALUE PASSES UNIT CIRCLE AT +1
C            =2 : EIGENVALUE PASSES UNIT CIRCLE AT -1
C            =3 : A PAIR PASSES WITH NONZERO IMAGINARY PART
C            =4 : UNDECIDED
C
C
C     DUMMY VARIABLES:
C     ---------------
C             The variables EVARE (REalpart of EigenVAlue),
C             EVAIM (IMaginary part) monitor the current eigenvalue
C             of the monodromy matrix that seems to cross the unit
C             circle of complex plane.  EVARE, EVAIM, IND, 
C             and IFLAG must NOT be changed by the calling routine.
C
C     NEEDED ROUTINES:
C     ---------------
C        1)   A ROUTINE FOR CALCULATING THE EIGENVALUES OF A MATRIX
C             (HERE EISPACK ROUTINES ARE USED)
C        2)   A ROUTINE FOR PROVIDING THE MONODROMY MATRIX (HERE
C             BOTH  MONOD  AND  MONOD2  ARE CALLED FOR COMPARISON
C             PURPOSES;  MONOD2 IS EXTREMELY CHEAP BUT OFTEN NOT
C             ACCURATE ENOUGH.)
C
C     TOLERANCES:
C     -----------
C         In this routine there are some checks needed whether an
C         eigenvalue is =1  (EPSONE)  or whether an eigenvalue
C         has a zero imaginary part  (EPSTOL).  the two
C         tolerances  EPSONE  and  EPSTOL  may be modified.
C         A CORRECT NUMERICAL INTERPRETATION OF THE BIFURCATION
C         BEHAVIOR DEPENDS CRITICALLY ON THESE TOLERANCES.
C==========================================================

      NROW=NDIML
      EPSTOL=0.1d0*EPSONE
      BOUNDM=1000.
      IF (NDIML.LT.N) THEN
         WRITE (UNITLST,*) 'WARNING: in STAPER dimension too small!'
         WRITE (*,*) 'WARNING: in STAPER dimension too small!'
         ENDIF
C  BOUNDM is a BOUND on the largest Modulus such that distinguishes 
C  between weakly unstable and strongly unstable orbits.
C  Only in the former case the more accurate and more expensive 
C  calculation of the stability is performed.
      MODEL='d'
      PARAM=YOLD(N+1,1)

CALCULATING MONODROMY MATRIX AND EIGENVALUES:
      IMON=2
1     IF (IMON.EQ.1) CALL MONOD (N,Y,EPS,A,INTE,UNITLST)
      IF (IMON.EQ.1) WRITE (UNITLST,612)
      IF (IMON.EQ.2) CALL MONOD2 (N,M,A)
      IF (IMON.EQ.2) WRITE (UNITLST,611)
      DO 2 I=1,N
      DO 2 J=1,N
2     B(I,J)=A(I,J)
      CALL EIGVAL (N,NROW,B,EIGRE,EIGIM,IERR,IDUMMY)
      DO 28 I=1,N
28    WRITE (UNITLST,602) EIGRE(I),EIGIM(I)

C=========================================================
C  INVESTIGATION OF STABILITY AND PERIODICITY:

      AMAX=0.d0
      APER=10.d0
      DO 3 I=1,N
      ARE=EIGRE(I)-1.d0
      AIM=EIGIM(I)
      ADIST=SQRT(ARE*ARE+AIM*AIM)
      IF (ADIST.LT.APER) THEN
                APER=ADIST
                IND1=I
                END IF
3     CONTINUE
      DO 31 I=1,N
      IF (I.EQ.IND1) GO TO 31
      ARE=EIGRE(I)
      AIM=EIGIM(I)
      ADIST=SQRT(ARE*ARE+AIM*AIM)
      IF (ADIST.GT.AMAX) THEN
                 AMAX=ADIST
                 INDA=I
                 END IF
31    CONTINUE
C  INDA is the index of the eigenvalue with the largest
C  modulus (=AMAX), IND1 is the index of the eigenvalue =1
C  APER is the defect of the candidate for an 1-eigenvalue.
      IF (APER.GT.EPSONE) IND1=0
      WRITE (UNITLST,604) APER
      IF (IND1.EQ.0) WRITE (UNITLST,*) ' PERIODICITY NOT CONFIRMED.  '
      IF (IND1.GT.0) WRITE (UNITLST,*) ' PERIODICITY CONFIRMED.      '
      IF (AMAX.GT.1.+EPSONE) WRITE (UNITLST,*) ' UNSTABLE SOLUTION   '
      IF (AMAX.GT.1.+EPSONE) WRITE (*,613) 'unst. ',Y(N+1,1),AMAX,APER 
      IF (AMAX.LT.1.-EPSONE) WRITE (UNITLST,*) ' STABLE SOLUTION     '
      IF (AMAX.LT.1.-EPSONE) WRITE (*,613) 'stable ',Y(N+1,1),AMAX,APER 
      IF (AMAX.GE.1.-EPSONE.AND.AMAX.LE.1.+EPSONE)
     1   WRITE (*,613) ' ',Y(N+1,1),AMAX,APER 
      WRITE (UNITLST,606) AMAX

      IF (IND1.GT.0) GO TO 4
      IF (IMON.EQ.1) GO TO 4
      IF (JACOBI.EQ.0) GO TO 4
      IF (NDIML.LT.N+N+2) THEN
         WRITE (UNITLST,*) 'NDIML too small for extensive stablity'
     1                     ,' analysis.'
         GO TO 4
         ENDIF
      IF (AMAX.GT.BOUNDM) GO TO 4
      IMON=1
      GO TO 1

C============================================================
C  COMPARISON WITH PREVIOUS PERIODIC SOLUTION:

4     IF (IFLAG.EQ.0) GO TO 7
      IF (IFLAG.GT.0 .AND. AMAX.GT.1.d0) GO TO 8
      IF (IFLAG.LT.0 .AND. AMAX.LT.1.d0) GO TO 8

C--  CHANGE IN STABILITY:

      AIMA=ABS(EIGIM(INDA))
      IF (AIMA.LT.EPSTOL.AND.ABS(EVAIM).GT.EPSTOL) GO TO 7
      IF (AIMA.GT.EPSTOL.AND.ABS(EVAIM).LT.EPSTOL) GO TO 7
      IF (ABS(EVAIM).LT.EPSTOL) GO TO 41

C     CASE BIFURCATION FROM/TO TORUS:
      IFLAG=3
      GO TO 6

41    ARE=EIGRE(INDA)
      IF (ARE*EVARE.LE.0.d0) GO TO 7

C     NOW BOTH ARE CONSIDERED TO BE REAL WITH SAME SIGN
      IF (ARE.GT.0.d0) IFLAG=1
      IF (ARE.LT.0.d0) IFLAG=2

6     WRITE (UNITLST,607) IFLAG
      WRITE (*,607) IFLAG
      TEXT='      '
      AMAXOL=SQRT(EVARE*EVARE+EVAIM*EVAIM)
      FAK=(1.d0-AMAXOL)/(AMAX-AMAXOL)
      PAR0=PARAM+FAK*(Y(N+1,1)-PARAM)
      PER0=YOLD(N+2,1) + FAK* (Y(N+2,1)-YOLD(N+2,1))
      DO 61 I=1,N
61    YINIT0(I)= YOLD(I,1) + FAK* (Y(I,1)-YOLD(I,1))
      IF (IFLAG.EQ.3) THEN
         WRITE (UNITLST,610) PAR0
         WRITE (*,610)  PAR0
         TEXT='Torus bif'
         TEXT2='PT'
         END IF
      IF (IFLAG.EQ.1) THEN
         WRITE (UNITLST,608) PAR0
         WRITE (*,608)  PAR0
         TEXT='branching'
         TEXT2='PB'
         END IF
      IF (IFLAG.EQ.2) THEN
         WRITE (UNITLST,609) PAR0
         WRITE (*,609)  PAR0
         TEXT='doubling'
         TEXT2='PD'
         END IF
      IF (NSECFLAG.EQ.0) WRITE (UNITOUT,9300) 'r',NOB,0,N,1,
     1      MODEL,TEXT2,(YINIT0(I),I=1,N),PAR0,PER0
      IF (NSECFLAG.EQ.1) WRITE (UNITOUT,9300) 'r',NOB,0,N,2,
     1      MODEL,TEXT2,(YINIT0(I),I=1,N),PAR0,PARCO,PER0
      WRITE (UNITLST,*) ' Guess of critical period: ', PER0
      WRITE (UNITLST,*) ' corresponding initial vector by linear'
     1                  ,' interpolation:'
      WRITE (UNITLST,601) (YINIT0(I),I=1,N)
      IF (NUMBUFF.LT.NDIMBUFF) THEN
         NUMBUFF=NUMBUFF+1
         Y0TEXT(NUMBUFF)=TEXT
         DO 62 I=1,N
62       Y0BUFF(I,NUMBUFF)=YINIT0(I)
         Y0BUFF(N+1,NUMBUFF)=PAR0
         Y0BUFF(N+2,NUMBUFF)=PER0
         HEAD1BUFF(1,NUMBUFF)=NOB
         HEAD1BUFF(2,NUMBUFF)=0
         JUDGEBUFF(NUMBUFF)=TEXT2
         HEAD1BUFF(3,NUMBUFF)=0
         STEPBUFF(NUMBUFF)=0.
         ENDIF
      GO TO 8

7     IFLAG=4
C     (undecidable case)

8     IND=INDA
      EVARE=EIGRE(IND)
      EVAIM=EIGIM(IND)
      IF (AMAX.LT.1.d0) IFLAG=-ABS(IFLAG)
      RETURN

601   FORMAT ('   ',4E19.11)
602   FORMAT ('  EIGENVAL. (FLOQUET MULTIPLIERS):',2E15.7)
604   FORMAT (' DEFECT OF EIGENVALUE 1 IS ',E13.5,' , ')
606   FORMAT ('        LARGEST MODULUS OF MULTIPLIERS: ',E15.7)
607   FORMAT (' *********************************',
     1       /' ***    CHANGE IN STABILITY    ***',I5,
     2       /' *********************************')
608   FORMAT (' BIFURCATION NEAR PARAMETER VALUE',E15.8)
609   FORMAT (' PERIOD DOUBLING NEAR PARAMETER VALUE',E15.8)
610   FORMAT (' TORUS BIFURCATION NEAR PARAMETER VALUE',E15.8)
611   FORMAT (' MONODROMY MATRIX TAKEN FROM SHOOTING ALGOR.')
612   FORMAT (' MONODROMY MATRIX RECALCULATED')
613   FORMAT (' ',a7,'at PAR.',E13.6,' largest modulus',E13.6,
     1   '. (defect:',E9.2,')')
9300  FORMAT (A1,I2,I4,2I2,A1,A2,18E20.13)
      END


      SUBROUTINE CONTD (NOB,NOS,N,M,T,YN,YOLD,IOLD,
     1      PARSTE,MAXCON,EPS,LEVEL,ITAIM,INTERACT,
     2      JACOBI,INTE,INDEX,EPSONE,Y0BUFF,NUMBUFF,Y0TEXT,
     3      HEAD1BUFF,STEPBUFF,JUDGEBUFF,
     4      UNITLST,UNITOUT,UNIT3,UNIT31)
C  THIS IS A SIMPLE CONTINUATION WITH RESPECT TO THE PARAMETER
C  (COMPONENT N+1) TO DEMONSTRATE THE USE OF ROUTINE STAPER.
C
C    INPUT:
C    -----    N,M,T,YN=Y,EPS       AS USUAL,
C                  (YN IS THE SOLUTION ARRAY WHERE THE CONTINUATION
C                  STARTS)
C             LEVEL   :  >0 for variable steps, =0 for fixed steps
C             PARSTE  :  (initial) step of parameter
C             MAXCON  :  NUMBER OF CONTINUATION STEPS TO REACH PAREND
C             JACOBI     =1 WHEN LINEARIZATION IS PROVIDED, ELSE =0 .
C             INDEX   :  INDEX OF NORMALIZATION (=KCO)
C             EPSONE  :  tolerance for deciding that an eigenvalue is =1,
C                        used for STAPER only.
C
      DOUBLE PRECISION T,Y,YN,YOLD,EPS,ALPHA,PARSTE,ETACO,
     1      ZETACO,PARCO,EVARE,EVAIM,EPSONE,Y0BUFF,STEPBUFF
      CHARACTER*9 Y0TEXT
      CHARACTER*2 JUDGEBUFF
      INTEGER HEAD1BUFF,UNITLST,UNITOUT,UNIT3,UNIT31
      PARAMETER (NDIM=12,NDIML=14,NDIMBUFF=6)
      DIMENSION T(7),Y(NDIML,7),YN(NDIML,7),YOLD(NDIML,7)
     1          ,Y0BUFF(NDIML,NDIMBUFF),Y0TEXT(NDIMBUFF),
     2  HEAD1BUFF(3,NDIMBUFF),STEPBUFF(NDIMBUFF),JUDGEBUFF(NDIMBUFF)
   
C    DIMENSIONS AT LEAST  T(M),Y(N+2,M)
      EXTERNAL FABZ,RABZ,INTE
      CHARACTER *2 TEXT
      COMMON /CONTRO/ ETACO,ZETACO,PARCO, NCO,N1CO,N2CO,
     1                KCO,IND2CO,NSECFLAG,NFLAG
c--------------
c  to be changed if necessary:
      maxfail=10
      IF (INTERACT.GT.0) THEN
       write (*,*) 'Check of stability:'
       write (*,*) 'Enter  1  if stability is to be fully checked,'
       write (*,*) 'enter  0  for a fast low-accuracy check,'
       write (*,*) 'enter -1  for NO CHECK.     ENTER:'
       read (*,*) ISTAB
       if (ISTAB.LT.-1.OR.ISTAB.GT.1) ISTAB=0
      ELSE
         ISTAB=1
      ENDIF
      if (ISTAB.EQ.0) EPSONE=0.5
c--------------
      NCO=N
      N2=N+2
      N2CO=N2
      KCO=INDEX
      ETACO=YN(N+1,1)
      IMIST=0
      IFLAG=0

      ETACO1=ETACO
      DO 8 ICON=1,MAXCON
1     ETACO=ETACO1+PARSTE
      DO 2 I=1,N2
      DO 2 J=1,M
      IF (ICON.GT.1.OR.IOLD.GT.0) THEN
	 ALPHA=PARSTE/(YN(N+1,1)-YOLD(N+1,1))
	 Y(I,J)=(1.+ALPHA)*YN(I,J) - ALPHA*YOLD(I,J)
	 ENDIF
      IF (ICON.EQ.1.AND.IOLD.EQ.0) THEN
	 Y(I,J)=YN(I,J)
	 ENDIF
2     CONTINUE
      IT=20
      WRITE (UNITLST,600)
      NFLAG=1
      CALL SOLVER (FABZ,RABZ,INTE,N2,M,T,Y,EPS,IT)
      IF (IT.LT.0) WRITE (UNITLST,601)
      WRITE (UNITLST,602) ICON,IT,Y(N+1,1),Y(N2,1)
      WRITE (*,602) ICON,IT,Y(N+1,1),Y(N2,1)
      IF (IT.GT.0) THEN      	
        WRITE (UNITLST,6031)
        WRITE (UNITLST,603) (Y(I,1),I=1,N)
        WRITE (*,6031)
        WRITE (*,603) (Y(I,1),I=1,N)
        ENDIF
      IF (IT.LT.0) THEN
	 IF (IMIST.GT.MAXFAIL) THEN
	      write (*,*) ' Termination: too many failures.'
	      write (UNITLST,*) ' Termination: too many failures.'
	      RETURN
	      ENDIF
	 IMIST=IMIST+1
	 PARSTE=PARSTE*0.3d0
c        if (PARSTE.lt. ...) return
	 GO TO 1
	 ENDIF
      REWIND UNIT3
      ETACO1=ETACO
      WRITE (UNIT3,605) (Y(I,1),I=1,N2)
      IF (ISTAB.GE.0)
     1 CALL STAPER (NOB,N,M,Y,EPS,YN,EVARE,EVAIM,IFLAG,IND,JACOBI,
     2             INTE,EPSONE,Y0BUFF,NUMBUFF,Y0TEXT,
     3             HEAD1BUFF,STEPBUFF,JUDGEBUFF,UNITLST,UNITOUT)
C     Here, EVARA,EVAIM,IND are used like dummy parameters.
      IF (IFLAG.EQ.0) TEXT='P '
      IF (IFLAG.GT.0) TEXT='PU'
      IF (IFLAG.LT.0) TEXT='PS'
      NOS=NOS+1
      IF (NSECFLAG.EQ.0) WRITE (UNITOUT,9300) 'R',nob,nos,n,1,
     1     'd',TEXT,(y(i,1),i=1,n+1),y(n2,1)
      IF (NSECFLAG.EQ.1) WRITE (UNITOUT,9300) 'R',nob,nos,n,2,
     1     'd',TEXT,(y(i,1),i=1,n+1),PARCO,y(n2,1)
      REWIND UNIT31
      WRITE (UNIT31,*) NOB
      WRITE (UNIT31,*) NOS
      WRITE (UNIT31,*) TEXT
      WRITE (UNIT31,*) 1
      WRITE (UNIT31,*) Y(N+1,1)-YN(N+1,1)
      DO 4 I=1,N2
      DO 4 J=1,M
      YOLD(I,J)=YN(I,J)
4     YN(I,J)=Y(I,J)
      IF (LEVEL.GT.0) PARSTE=PARSTE*DBLE(ITAIM)/DBLE(IT)
8     CONTINUE
      RETURN
600   FORMAT (' ')
601   FORMAT (' NO SUCCESS')
602   FORMAT (/' step ',I3,' . iter.: ',I3,' .  PAR.:',E19.11,
     1  ' . PERIOD:',E19.11)
603   FORMAT ('  ',4E19.11)
6031  FORMAT (' initial vector:')
605   FORMAT (E20.10)
9300  format (a1,i2,i4,2i2,a1,a2,18e20.13)
      END


      SUBROUTINE CHECKPHASE (Y,INDEXPH)
      DOUBLE PRECISION T1,Y,Y1,DY,ETACO,ZETACO,PARCO,S
      PARAMETER (NDIM=12,NDIML=14)
      DIMENSION Y(NDIML,7),Y1(NDIML),DY(NDIML)
      EXTERNAL FABZ
      COMMON /CONTRO/ ETACO,ZETACO,PARCO, NCO,N1CO,N2CO,
     1                KCO,IND2CO,NSECFLAG,NFLAG
c-------------- This routine checks which of the phase conditions
c               dy(i)=f(i)=0 was adopted on the solution Y.
c               Output is the index INDEXPH=i .
      DO 1 I=1,NCO+2
      Y1(i)=Y(i,1)
1     CONTINUE
      CALL FABZ (T1,Y1,DY)
      S=999.d0
      I1=0
      DO 2 I=1,NCO
      IF (ABS(DY(I)).LT.S) THEN
          S=ABS(DY(I))
          I1=I
          ENDIF
2     CONTINUE
      IF (I1.EQ.0 .OR. S.GT.1.D-3) THEN
          INDEXPH=0
        ELSE
          INDEXPH=I1
        ENDIF
      RETURN
      END
