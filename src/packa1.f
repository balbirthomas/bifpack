C*************************************
C**  This is  packA1 of BIFPACK   latest update:  2. Jun 1999  
C**  ---------------                **
C**    including routines for a     **
C**    stability and bifurcation analysis
C**    of systems of 'algebraic' equations
C
C      In addition, this file includes  I N T E R A
C      for convenient interactive path following.
C
C*************************************



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

c============================================================================
      SUBROUTINE CONTA (NOB,NOS,N,FCN,YN,YA,LAST,EPS,LEVEL,INDEX,
     1   STEP,STEPMA,STEPMI,STEPMP,MAXCON,ITAIM,ITMAX,WEIGHT,IEXTRA,
     2   IFINAL,YFINAL,ISTAB,IRUN,IPRINT,WINDOW,BOUND2,Y0BUFF,
     3   NUMBUFF,Y0TEXT,HEAD1BUFF,STEPBUFF,JUDGEBUFF
     4   ,UNITLST,UNITOUT,UNIT3,UNIT31)

C  CONTINUATION ALGORITHM,
C  VERSION "A" for systems of nonlinear equations

C    INPUT VARIABLES:
C    ----------------
C    N     : number of original equations (without auxiliary eq.)
C    FCN   : name of routine where equations are implemented (see below)
C    YN    : solution from which the continuation is to start.
C    YA    : a previous solution, may be empty.
C    LAST  : =1  YA carries another solution previously calculated 
C                on the same branch.
C            =0  YA is empty or void.
C    EPS   : accuracy requirement for solver
C    LEVEL : =2 standard case of variable stepsize: both
C               INDEX and steplength  STEP  may be varied.
C            =1 only  STEP  may vary, INDEX is kept fixed.
C            =0 fixed steps: both  STEP  and INDEX are kept fixed.
C    INDEX : index of component to be fixed in first step.
C    STEP  : initial stepsize  (INDEX and  STEP  may vary.)
C    STEPMA: maximum stepsize (see remarks below).
C            STEPMA is measured along the branch (predictor).
C    STEPMI: minimum stepsize, measured along the branch (predictor).
C    STEPMP: maximum stepsize in the parameter                       
C    MAXCON: maximum number of contin. steps
C    ITAIM : preferred number of iterations for solver.
C             (see table below)
C    ITMAX : maximum number of iterations allowed for the
C            nonlinear equation solver.
C        The following table gives hints on how to choose ITAIM, ITMAX:
C            NONL.SOLVER                     EPS     ITAIM     ITMAX
C            --------------------------------------------------------
C            NEWTON (QUADRATIC CONVERGENCE)  1.E-5     3        5-6
C            QUASI-N. (SUPERLINEAR CONV.)    1.E-5     6        10
C            CHORD METHOD (LINEAR CONV.)     1.E-5    10        15
C        (or similar; allow more iterations for more stringent accuracies.)
C
C    WEIGHT: weight for preferred use of the physical parameter.
C    IEXTRA: =1 for secant extrapolation, =0: no extrapolation.
C    IFINAL, YFINAL: required for the option of target points, see below.
C    ISTAB : =1   if stability analysis is required (see below)
C    IRUN  : number of continuation run; only used for output.
C    IPRINT: =1 : data of step control are printed. (=0: omitted)
C    WINDOW: The two values in WINDOW define the parameter interval of
C            interest. The continuation stops if a solution is
C            calculated with a parameter value ouside the window.
C            WINDOW(1) : left bound,
C            WINDOW(2) : right bound.

C     HOW TO REACH TARGET POINTS  (IFINAL=0 FOR NO TARGET POINT):
C     ---------------------------
C           No target point required: set  IFINAL=0 , YFINAL is dummy then.
C        IFINAL:   index of Y-component of target point.
C        YFINAL :  target value for Y(IFINAL)
C            There is no target point check carried out in the
C            first continuation step.

C     OUTPUT:
C     ------   YN, YA, INDEX, STEP     are updated.
C
C     OUTPUT for graphics etc. on file 'DATOUT' (UNIT 8):
C     ---------------------------------------------------
C        a standard record of each solution is written via FORMAT 9300

C   ONLY IN CASE   ISTAB.GT.0  :
C   ---------------------------
C      in case the function defined in FCN is the right-hand side
C      of a system of ordinary differential equations, a stability
C      analysis may be of interest.
C      this routine calls the routine   S T A S T A   for carrying out a
C      STAbility analysis of STAtionary solutions of ordinary differential
C      equations. (omitted in case  ISTAB.LE.0 )
C
C    NEEDED SUBROUTINE FCN (T,Y,F)
C    -------------------------  defines the function  F(Y,PARAM)
C            FCN requires the following three statements:
C                   COMMON /CONTRO/ ETACO,ZETACO,PARCO, NCO,N1CO,N2CO,
C                  1                KCO,IND2CO,NSECFLAG,NFLAG
C                   PARAM=Y(N1CO)
C                   F(N1CO)=Y(KCO)-ZETACO*Y(IND2CO)-ETACO
C             (standard cont. takes ZETA=0 and IND2CO arbitrary.)
C
C      REMARKS ON STEPMI/MA/MP:
C      -----------------------
C        (all to be chosen nonnegativ)
C           if no limitaion is required, take STEPMI=0 and
C           STEPMA = large number.
C           if constant step is required, use  LEVEL=0  .
C
C     STEP CONTROL: 
C           The step control takes curvature into account. Two internal 
C           constants BOUND1 and BOUND2 are responsible: The step 
C           control aims at:
C                      curvature * stepsize  .le. BOUND1   ,
C                      relative change in the solution  .le. BOUND2 .
C           We recommend BOUND1=0.1, or 0.15 , and BOUND2=0.1 .
C           The latter bound means a careful/pessimistic control
C           where the steps are limited such that the solution 
C           (including parameter) does not change more than 10 percent.
C           It only applies as long as norm of the solution is larger
C           than BOUND3 (=1.). BOUND2=0.1 creates much information
C           along the branch. For some examples and applications such
C           as homotopy larger values of BOUND2 may be advisable.
C           The initial step length is not adjusted.
C
C     TURNING POINTS:
C           An approximation (based on quadratic interpolation) of
C           turning points is calculated and written to the output.
C           The accuracy is not worse than that of a graphically determined
C           approximation that is based on the continuation output.
C**********************************************************************

      EXTERNAL FCN,FCNLIN
      COMMON /CONTRO/ ETACO,ZETACO,PARCO, NCO,N1CO,N2CO,
     1                KCO,IND2CO,NSECFLAG,NFLAG
      DOUBLE PRECISION Y,YA,YN,EPS,STEP,STEPMA,STEPMI,WEIGHT,ETACO,
     1  ZETACO,PARCO,D,D1,D2,YNORM,ALPHA,DUMMY1,DUMMY2,PREDNO,SECNO,
     2  SECNO1,YFINAL,WINDOW,CURV,CURV1,CURVNE,SCALPR,SUM2,C1,C2,CY,
     3  CPAR,BETA,BOUND1,BOUND2,BOUND3,STEPMP,FACT,Y0BUFF,YTEST
     4  ,STEPAR,STEPBUFF
      INTEGER NUMBUFF,HEAD1BUFF,UNITLST,UNITOUT,UNIT3,UNIT31
      CHARACTER*2 TEXT,JUDGEBUFF
      CHARACTER*9 Y0TEXT
      PARAMETER (NDIM=12,NDIML=14,NDIMBUFF=6)
      DIMENSION Y(NDIML),YA(NDIML),YN(NDIML),WINDOW(2),YTEST(NDIML),
     1          DUMMY1(2),DUMMY2(2),IDUMM3(2),Y0BUFF(NDIML,NDIMBUFF)
     2          ,Y0TEXT(NDIMBUFF),HEAD1BUFF(3,NDIMBUFF),
     3          STEPBUFF(NDIMBUFF),JUDGEBUFF(NDIMBUFF)
C     DIMENSIONS ALL  N+1
      WRITE (*,*) ' BIFPACK: entering routine CONTA.'
      WRITE (UNITLST,*) ' BIFPACK: entering routine CONTA.'
      KCO=INDEX
      IF (WEIGHT.LE.0.) WEIGHT=1.
      IF (ITAIM.LE.0) ITAIM=3
      IFLAG=0
      ICON=1
      ZETACO=0.d0
      IND2CO=1
      N1CO=N+1
      BOUND1=0.15d0
C      WRITE (*,*) 'ENTER BOUND FOR ANGLE :  '
C      READ (*,*) BOUND1
      BOUND3=1.d0
      IF (IPRINT.GT.0 .AND.LEVEL.GT.0) THEN
         WRITE (UNITLST,*) 'CHOSEN LIMITATIONS FOR STEP LENGTHS:'
         WRITE (UNITLST,*) '   ALONG THE BRANCH (PREDITOR):', STEPMA
         WRITE (UNITLST,*) '   IN THE PARAMETER:           ', STEPMP
         WRITE (UNITLST,*) '   RELATIVE CHANGE IN SOLUTION:', BOUND2               
         WRITE (UNITLST,*) '                 ' 
         END IF
      WRITE (UNITLST,603) 1
      IF (LAST.EQ.1) THEN
           D=YN(KCO)-YA(KCO)
           IF (0..LE.D .AND. D.LT.EPS) ALPHA=1.d0
           IF (-EPS.LT.D .AND. D.LE.0.) ALPHA=-1.d0
           IF (ABS(D).GE.EPS) ALPHA=STEP/D
           END IF
      N1=N1CO
      K2=0                                                              *
      GO TO 3                                                           *
C************************************************************************

1     ICON=ICON+1
      WRITE (UNITLST,603) ICON
      IF (LEVEL.EQ.0) THEN
           ALPHA=1.d0
           GO TO 3
           END IF
      IF (LEVEL.EQ.1) GO TO 2
C..  DETERMINATION OF INDEX KCO=K:
      KCO=0
      D1=0.d0
      K2=0
      D2=0.d0
      DO 15 I=1,N1CO
      YNORM=ABS(YN(I))
      IF (YNORM.LT.1.) YNORM=1.d0
      D=ABS(YN(I)-YA(I)) /YNORM
      IF (I.EQ.N+1 .AND. D.GE.0.01) D=D*WEIGHT
      IF (D.GT.D1) THEN
         D2=D1
         K2=KCO
         D1=D
         KCO=I
      ELSE
         IF (D.GE.D2) THEN
         D2=D
         K2=I
         END IF
      END IF
15    CONTINUE
      IF (D1.LT.EPS) then
          write (UNITLST,*) 'CONTA: step control fails. D1 too small.'
          GO TO  9
          endif              
C..  ADJUSTING STEPSIZE:
2     STEP=YN(KCO)-YA(KCO)
      ALPHA=FLOAT(ITAIM)/FLOAT(IT)
      if (alpha.gt.2.) alpha=2.d0
      if (alpha.lt.0.5) alpha=0.5d0 
      IF (ICON.LT.4) GO TO 25
      CURVNE = CURV + (SECNO/SECNO1) * (CURV-CURV1)
      IF (CURVNE.LT.EPS) CURVNE=EPS
      IF (CURVNE.LT.CURV) CURVNE=0.25d0*(3.d0*CURVNE+CURV)
      BETA=BOUND1/(CURVNE*ALPHA*SECNO)
      IF (BETA.LT.1.) ALPHA=ALPHA*BETA
25    PREDNO=ALPHA*SECNO
      IF (SUM2.LE.BOUND3) SUM2=-1.d0
      IF (SUM2.GT.BOUND3) SUM2=PREDNO/SUM2
      IF (IPRINT.GT.0) THEN
              WRITE (UNITLST,*) ' PROPOSED STEP LENGTH FACTOR:',ALPHA
              WRITE (UNITLST,*) ' LENGTH OF FREE PREDICTOR:  ',PREDNO
              IF (SUM2.GT.0.) 
     1        WRITE (UNITLST,*) ' RELATIVE CHANGE IN SOLUTION:',SUM2
              END IF
      STEPAR=ALPHA*ABS((YN(N1)-YA(N1)))
      IF (STEPMI.LE.PREDNO.AND.PREDNO.LE.STEPMA.AND.
     1    STEPAR.LE.STEPMP.AND.SUM2.LE.BOUND2.AND.IFINAL.EQ.0) GO TO 295
      FACT=1.
C  limitations of step length:   
      IF (PREDNO.GT.STEPMA) THEN
            FACT=STEPMA/PREDNO
            ALPHA=ALPHA*FACT
            STEPAR=STEPAR*FACT
            SUM2=SUM2*FACT
            PREDNO=STEPMA
            END IF
      IF (PREDNO.LT.STEPMI) THEN
            FACT=STEPMI/PREDNO
            ALPHA=ALPHA*FACT
            STEPAR=STEPAR*FACT
            SUM2=SUM2*FACT
            PREDNO=STEPMI
            END IF
      IF (STEPAR.GT.STEPMP) THEN
            FACT=STEPMP/STEPAR
            ALPHA=ALPHA*FACT
            PREDNO=PREDNO*FACT
            SUM2=SUM2*FACT
            END IF
      IF (SUM2.GT.BOUND2) THEN
            FACT=BOUND2/SUM2
            ALPHA=ALPHA*FACT
            PREDNO=PREDNO*FACT
            END IF
      IF (IFINAL.GT.0) THEN
            Y(IFINAL)=YN(IFINAL)+ALPHA*(YN(IFINAL)-YA(IFINAL))
C..         CHECK FOR FINAL STEP TO REACH TARGET POINT
            IF ( YN(IFINAL).LT.YFINAL.AND.YFINAL.LT.Y(IFINAL) .OR.
     1         Y(IFINAL).LT.YFINAL.AND.YFINAL.LT.YN(IFINAL) ) THEN
                KCO=IFINAL
                STEP=YN(IFINAL)-YA(IFINAL)
                FACT=ABS((YFINAL-YN(IFINAL))/(Y(IFINAL)-YN(IFINAL)))
                ALPHA=ALPHA*FACT
                PREDNO=PREDNO*FACT
                END IF
            END IF
      IF (IPRINT.GT.0 .AND. ABS(FACT-1.).GT.EPS) THEN
            WRITE (UNITLST,*) ' step length limited to:    ',PREDNO
            WRITE (UNITLST,*) ' final step length factor:',ALPHA
            END IF
295   STEP=STEP*ALPHA

3     IFAIL=0
C..  INITIAL GUESS:
      if (iprint.gt.0) write (UNITLST,*) ' final alpha=',alpha 
31    DO 35 I=1,N1CO
      IF (LAST.EQ.1 .AND. IEXTRA.EQ.1)
     1     Y(I)=YN(I)+ALPHA*(YN(I)-YA(I))
      IF (LAST.EQ.0 .OR. IEXTRA.EQ.0) Y(I)=YN(I)
35    CONTINUE
      ETACO=YN(KCO)+STEP
      IT=-1 
      CALL NOLE (N1,Y,FCN,EPS,ITMAX,IT)

      IF (IT.GT.0) THEN
         REWIND UNIT3
         WRITE (UNIT3,680) (Y(I),I=1,N1CO)
         WRITE (UNITLST,600) IT
         WRITE (UNITLST,6001) KCO,STEP
	 NOS=IRUN*100+ICON
         IF (LAST.EQ.0) GO TO 45
         IF ( (Y(N1)-YN(N1))*(YN(N1)-YA(N1)) .LT.0. ) THEN
C             approximation of turning point (parameter and y(IC) ) :
	      CPAR=0.
	      DO 43 IC=1,N
	      IF (ABS(YN(IC)-Y(IC)).LT.EPS.or.
     1            ABS(YA(IC)-Y(IC)).LT.EPS.or.
     2            ABS(YA(IC)-YN(IC)).LT.EPS) THEN
		  YTEST(IC)=0.d0
		  GO TO 43
		  ENDIF
              C1=(YN(N1)-Y(N1))/(YN(IC)-Y(IC))
              C2=(YA(N1)-Y(N1))/(YA(IC)-Y(IC))
              C2=(C2-C1)/(YA(IC)-YN(IC))
              CY=0.5*(Y(IC)+YN(IC)-C1/C2)
	      YTEST(IC)=CY
	      IF (CPAR.LT.EPS) CPAR=Y(N1)+(CY-Y(IC))*(C1+C2*(CY-YN(IC)))
43            CONTINUE
C  if required enter here iterative refinement
              IF (ABS(CPAR).GT.EPS) THEN
	        IF (NSECFLAG.EQ.0) WRITE (UNITOUT,9300) 'r',NOB,0,N,1,
     1             'a','T ',(YTEST(IC),IC=1,N),CPAR
	        IF (NSECFLAG.EQ.1) WRITE (UNITOUT,9300) 'r',NOB,0,N,2,
     1             'a','T ',(YTEST(IC),IC=1,N),CPAR,PARCO
	        IC=1
                WRITE (UNITLST,605)
	        WRITE (UNITLST,604) CPAR,IC,YTEST(IC)
                WRITE (UNITLST,605)
                WRITE (*,605)
	        WRITE (*,604) CPAR,IC,YTEST(IC)
                WRITE (*,605)
                ENDIF 
              END IF 
45       SCALPR=0.d0
         SECNO=0.d0
         SUM2=0.d0
         SECNO1=0.d0         
         STEPAR=Y(n+1)-YN(N+1)
         DO 5 I=1,N1
         SUM2=SUM2+Y(I)*Y(I)
         SECNO=SECNO+(Y(I)-YN(I))*(Y(I)-YN(I))
         IF (LAST.EQ.1) THEN
                SCALPR=SCALPR+(Y(I)-YN(I))*(YN(I)-YA(I))
                SECNO1=SECNO1+(YN(I)-YA(I))*(YN(I)-YA(I))
                END IF
         YA(I)=YN(I)
5        YN(I)=Y(I)
         SECNO=SQRT(SECNO)
         SUM2=SQRT(SUM2)
         IF (IPRINT.GT.0) THEN
               WRITE (UNITLST,*) ' LENGTH OF SECANT:          ', SECNO
               WRITE (UNITLST,*) ' L-2 NORM OF SOLUTION: ',SUM2
               END IF
         IF (LAST.EQ.1) THEN
               SECNO1=SQRT(SECNO1)
               SCALPR=ABS(SCALPR)/(SECNO*SECNO1)
               IF (1.-SCALPR .GT.10.*EPS) THEN
                      SCALPR=ACOS(SCALPR)
                      ELSE
                      SCALPR=0.
                      END IF
               IF (ICON.GT.2) CURV1=CURV
               CURV=SCALPR/SECNO
               IF (IPRINT.GT.0) WRITE (UNITLST,*) ' CURVATURE: ', CURV
               END IF
         LAST=1
	 CALL PRINT (Y,NOS,UNITLST)
C        THE FOLLOWING STATEMENT CALLS A STABILITY ANALYSIS:
	 IF (ISTAB.GT.0) CALL STASTA (NOB,N,Y,YA,
     1        DUMMY1,DUMMY2,IFLAG,IDUMM3,Y0BUFF,NUMBUFF,Y0TEXT,
     2        HEAD1BUFF,STEPBUFF,JUDGEBUFF,UNITLST,UNITOUT)
	 TEXT='  '
         IF (ISTAB.GT.0 .AND. IFLAG.EQ.1)  TEXT='U '
         IF (ISTAB.GT.0 .AND. IFLAG.EQ.-1) TEXT='S '
	 IF (NSECFLAG.EQ.0) WRITE (UNITOUT,9300) 'R',NOB,NOS,N,1,'a',
     1                  TEXT,(Y(I),I=1,N),Y(N1CO)
	 IF (NSECFLAG.EQ.1) WRITE (UNITOUT,9300) 'R',NOB,NOS,N,2,'a',
     1                  TEXT,(Y(I),I=1,N),Y(N1CO),PARCO
      REWIND UNIT31
      WRITE (UNIT31,*) NOB
      WRITE (UNIT31,*) NOS
      WRITE (UNIT31,*) TEXT
      WRITE (UNIT31,*) 7
      WRITE (UNIT31,*) STEPAR
	 IF (IFINAL.EQ.0) GO TO 51
         IF (ABS(YFINAL-YN(IFINAL)).LT.EPS) GO TO 99
51       IF (YN(N1).LT.WINDOW(1) .OR. YN(N1).GT.WINDOW(2)) THEN
               WRITE (*,*) ' Continuation has left the window.'
               WRITE (UNITLST,*) ' Continuation has left the window.'
               GO TO 99
               END IF
       IF (ICON.LT.MAXCON) GO TO 1

C************************************************************************
C..  FAIL- AND END-CONDITIONS:                                          *
         GO TO 99                                                       *
      ELSE
         WRITE (UNITLST,6002) IT
         WRITE (UNITLST,6001) KCO,STEP
         IF (LEVEL.EQ.0) THEN
             write (UNITLST,*) 'CONTA: failure of step control.'
             GO TO 9
             endif
         IF (IFAIL.LT.2) THEN
C..          REDUCTION AFTER FAILURE:
              IFAIL=IFAIL+1
              STEP=STEP*0.2d0
              IF (LAST.EQ.1) ALPHA=ALPHA*0.2d0
              GO TO 31
         ELSE
              IF (KCO.EQ.K2.OR.K2.EQ.0) THEN
                 WRITE (UNITLST,*) 'NO CONVERGENCE'
                 GO TO 9
                 END IF
C..           TRY NEXT INDEX:
              IF (LEVEL.EQ.1) GO TO 9
              KCO=K2
              IT=ITAIM
              GO TO 2
         END IF
      END IF

9     INDEX=-1
      RETURN
99    INDEX=KCO
      RETURN
600   FORMAT ('  Solution obtained after ',I3,' iterations.')
6001  FORMAT (' STEP CONTROL: index of fixed component was ',I3,
     2        ' ,   stepsize: ',E10.4)
6002  FORMAT ('  FAILURE ',I3) 
601   FORMAT (' ****  parameter value: ',E17.9,'   ****',
     1       /                   27X,14('-'))
680   FORMAT (E20.10)
9300  FORMAT (A1,I2,I4,2I2,A1,A2,18E20.13)
602   FORMAT (6X,4E17.9)
603   FORMAT (/' continuation step NO.:',I4)
604   FORMAT (' ***  TURNING OF BRANCH AT PARAMETER: ',E13.6,
     1        ' , Y(',I3,')= ',E13.6)
605   FORMAT ('*************************************************')
      END


      SUBROUTINE PRINT (Y,NOS,UNITLST)
      PARAMETER (NDIM=12,NDIML=14)
      DOUBLE PRECISION Y(NDIML),ETACO,ZETACO,PARCO,CORR,RESID,EPMACH
      INTEGER UNITLST
      COMMON /CONTRO/ ETACO,ZETACO,PARCO, NCO,N1CO,N2CO,
     1                KCO,IND2CO,NSECFLAG,NFLAG
      COMMON /INFONO/ IWEAK,NZAEHL,RESID,CORR
C   This is a printing routine that is specifically designed to be used
c   along with NOLE, since via the COMMON statement /INFONO/ 
C   information on the accuracy of the solution is provided.
      EPMACH =1.d-12
      IF (CORR.LT.EPMACH) CORR=EPMACH
      NDIGIT=-IDINT(LOG10(CORR))
      WRITE (UNITLST,600) NOS,Y(N1CO)
      WRITE (*,600) NOS,Y(N1CO)
      WRITE (UNITLST,601) (Y(I),I=1,NCO)
      WRITE (*,601) (Y(I),I=1,NCO)
      WRITE (UNITLST,602) NDIGIT
      RETURN
600   FORMAT ('(step',i4,
     1    ')   parameter: ',E19.11,' ,     sol. vector:')
601   FORMAT (3X,4E19.11)
602   FORMAT ('   Accuracy of Newton-solution about ',I2,
     1        ' significant digits.')
      END


      SUBROUTINE STASTA (NOB,N,Y,YOLD,EVARE,EVAIM,IFLAG,IND,
     1  Y0BUFF,NUMBUFF,Y0TEXT,HEAD1BUFF,STEPBUFF,JUDGEBUFF
     2  ,UNITLST,UNITOUT)
      DOUBLE PRECISION Y,Z,B,EIGRE,EIGIM,EVARE,EVAIM,REALMA
     1  ,REALMI,EIGIMI,FAK,GUESSI,ZAHL,PERIOD,TWOPI,YOLD,Y0
     2  ,Y0BUFF,ETACO,ZETACO,PARCO,STEPBUFF
      INTEGER IDUMMY,NUMBUFF,HEAD1BUFF,UNITLST,UNITOUT
      CHARACTER*1 LABEL
      CHARACTER*2 TEXT2,JUDGEBUFF
      CHARACTER*9 Y0TEXT
      PARAMETER (NDIM=12,NDIML=14,NDIMBUFF=6)
      DIMENSION Y(NDIML),Z(NDIML,NDIML),B(NDIML,NDIML),EIGRE(NDIML),
     1          EIGIM(NDIML),IDUMMY(NDIML),YOLD(NDIML),Y0(NDIML),
     2          EVARE(2),EVAIM(2),IND(2),HEAD1BUFF(3,NDIMBUFF),
     3          Y0BUFF(NDIML,NDIMBUFF),Y0TEXT(NDIMBUFF),
     4          JUDGEBUFF(NDIMBUFF),STEPBUFF(NDIMBUFF)
C     DIMENSION AT LEAST: Y(N+1),Z(N,N),B(NROW,N),EIGRE(N),EIGIM(N)
C     NROW TO BE SPECIFIED BELOW  (NUMBER OF ROWS IN THE DIMENSIONS)

      COMMON /CONTRO/ ETACO,ZETACO,PARCO, NCO,N1CO,N2CO,
     1                KCO,IND2CO,NSECFLAG,NFLAG 
      COMMON /MATRIX/ Z
CAN BE EXCHANGED BY (TDUMMY to be defines as DOUBLE PRECISION):
C     CALL FCNLIN (TDUMMY,Y,Z)

C   INPUT: The JACOBIAN matrix is needed. Here the matrix from the
C   -----  solution procedure  N O L E  is provided via the above
C          COMMON statement.
C          (possible alternatives: recalculate the JACOBIAN  Z  by means
C          of the program  CHECKD ,  or  call FCNLIN  if available.)
C          Y is the current solution, YOLD is the previous solution.
C
C   IFLAG  is needed to compare two consecutive runs of  stasta.
C   -----  =0 :  on first call of stasta
C          =1 :  previous solution was unstable
C          =-1 : previous solution was stable
C
C   DUMMY VARIABLES: The variables  EVARE  (realpart of eigenvalues)
C   -----  and  EVAIM  (imaginary part) monitor eigenvalues which
C          might cross the imaginary axis. Their indices are stored in
C          IND.    EVARE, EVAIM, IFLAG, IND    must not be changed
C          by the calling routine.
C
C   The internal variable Y0 carries linear approximations of branch
C           points that are detected. This approximation is written to
C           file 'diagram' , with -1 or -2 as the first entry of the record,
C           and to the BUFFER.
C    WARNING:  NOT EVERY OUTPUT INDICATING A BIFURCATION IS TRUSTWORTHY
C    -------   (DUE TO TOO LARGE STEPLENGTHS  AND/OR  PERMUTATIONS OF
C              EIGENVALUES).
C************************************************************

      TWOPI=6.283185307179586d0
      NROW=NDIML
      NEIG=N
C     For printing eigenvalues, set NEIG=N, otherwise NEIG=0 .

      DO 2 I=1,N
      DO 2 J=1,N
2     B(I,J)=Z(I,J)
      CALL EIGVAL (N,NROW,B,EIGRE,EIGIM,IERR,IDUMMY)

Calculation of the maximum realpart (=REALMA)
C---------------          and of the realpart closest to zero:

      REALMA=-1.E10
      REALMI=1.E10
      IF (IERR.GT.0) WRITE (UNITLST,604)
      IF (IERR.GT.0 .AND. NEIG.GT.0) NEIG=IERR-1
      DO 3 I=1,NEIG
      IF (NEIG.LT.9 .OR. ABS(EIGRE(I)).LT.10.d0)
     1       WRITE (UNITLST,613) EIGRE(I),EIGIM(I)
C     (In large systems not all eigenvalues are printed.)
      IF (EIGRE(I).GT.REALMA) THEN
                REALMA=EIGRE(I)
                EIGIMA=EIGIM(I)
                INDMA=I
                END IF
      IF (ABS(EIGRE(I)).LT.ABS(REALMI)) THEN
                REALMI=EIGRE(I)
                EIGIMI=EIGIM(I)
                INDMI=I
                END IF
3     CONTINUE

      IF (REALMA.GT.0.) WRITE (UNITLST,601) REALMA
      IF (REALMA.LT.0.) WRITE (UNITLST,602) REALMA
      WRITE (UNITLST,603) REALMI,EIGIMI

Comparing the stability with that of the previous solution:
C----------------------------------------------------------

      IF (IFLAG.EQ.0) GO TO 6
      ICASE=0
      INDICA=-1
      IF (REALMA*EVARE(1) .LT. 0.d0) THEN
              ICASE=1
              FAK=EVARE(1) / (EVARE(1)-REALMA)
              IF (INDMA.EQ.IND(1)) GO TO 4
              END IF
      IF (REALMI*EVARE(2) .LT. 0.d0) THEN
              FAK=EVARE(2) / (EVARE(2)-REALMI)
              IF (INDMI.NE.IND(2)) INDICA=-2
              ICASE=2
              GO TO 4
              END IF
      IF (ICASE.EQ.0) GO TO 6
      INDICA=-2

4     IF (INDICA.EQ.-2) WRITE (UNITLST,608)
C    Output will be emphasized when stability changes (ICASE=1).
      IF (ICASE.EQ.1) WRITE (UNITLST,605)
      IF (ICASE.EQ.1) WRITE (*,605)
      IF (ICASE.EQ.1) WRITE (UNITLST,605)
      WRITE (UNITLST,606)
      IF (ICASE.EQ.1) WRITE (UNITLST,607)

C    linear interpolation to the zero of the realpart of eigenvalue:
      DO 5 I=1,N+1
5     Y0(I)=YOLD(I) + (Y(I)-YOLD(I)) * FAK
      WRITE (UNITLST,609) Y0(N+1)
      IF (ICASE.EQ.1) ZAHL=EIGIMA
      IF (ICASE.EQ.2) ZAHL=EIGIMI
      GUESSI=EVAIM(ICASE)+(ZAHL-EVAIM(ICASE))*FAK
      IF (ICASE.EQ.1.AND.NUMBUFF.LT.NDIMBUFF) THEN
         NUMBUFF=NUMBUFF+1
         DO 51 I=1,N+1
51       Y0BUFF(I,NUMBUFF)=Y0(I)
         Y0TEXT(NUMBUFF)='        '
         Y0BUFF(N+2,NUMBUFF)=0.d0
         HEAD1BUFF(1,NUMBUFF)=NOB
         HEAD1BUFF(2,NUMBUFF)=0
         STEPBUFF(NUMBUFF)=0.
         WRITE (*,*) ' Critical solution written to BUFFER',NUMBUFF
         ENDIF
      IF (ABS(GUESSI).GT.0.0001d0 .AND. ABS(EVAIM(ICASE)).GT.1.d-6
     1    .AND. ABS(ZAHL).GT.1.E-6) THEN
        LABEL='-'
        TEXT2='  '
        PERIOD=TWOPI/GUESSI
        WRITE (UNITLST,610) PERIOD
        If (ICASE.EQ.1.AND.NUMBUFF.LT.NDIMBUFF)
     1              Y0BUFF(N+2,NUMBUFF)=PERIOD
        IF (ICASE.EQ.1) THEN
          WRITE (*,605)
          WRITE (*,614) Y0(N+1),PERIOD
          TEXT2='HB'
          Y0TEXT(NUMBUFF)='Hopf Bif.'
          LABEL='r'
          JUDGEBUFF(NUMBUFF)=TEXT2
          HEAD1BUFF(3,NUMBUFF)=2
          ENDIF
        IF (NSECFLAG.EQ.0) WRITE (UNITOUT,9300) LABEL,nob,0,n,1,
     1          'a',TEXT2,(Y0(I),I=1,N+1),PERIOD
        IF (NSECFLAG.EQ.1) WRITE (UNITOUT,9300) LABEL,nob,0,n,2,
     1          'a',TEXT2,(Y0(I),I=1,N+1),PARCO,PERIOD
	 ELSE
        IF (ICASE.EQ.1) THEN
          Y0TEXT(NUMBUFF)='real B.P.'
          TEXT2='BT'
          LABEL='r'
          JUDGEBUFF(NUMBUFF)=TEXT2
          HEAD1BUFF(3,NUMBUFF)=0
          ENDIF
        IF (NSECFLAG.EQ.0) WRITE (UNITOUT,9300) LABEL,nob,0,n,1,
     1          'a',TEXT2,(Y0(I),I=1,N+1)
        IF (NSECFLAG.EQ.1) WRITE (UNITOUT,9300) LABEL,nob,0,n,2,
     1          'a',TEXT2,(Y0(I),I=1,N+1),PARCO
	 END IF
      WRITE (UNITLST,*) 
     1    '**   CRITICAL SOLUTION BY LINEAR INTERPOLATION:   '
      WRITE (UNITLST,615) (Y0(I),I=1,N)
      IF (ICASE.EQ.1) WRITE (UNITLST,605)
      IF (ICASE.EQ.1) WRITE (UNITLST,605)
      IF (ICASE.EQ.1) WRITE (*,605)




6     IF (REALMA.GT.0.d0) IFLAG=1
      IF (REALMA.LT.0.d0) IFLAG=-1
      IND(1)=INDMA
      EVARE(1)=REALMA
      EVAIM(1)=EIGIMA
      IND(2)=INDMI
      EVARE(2)=REALMI
      EVAIM(2)=EIGIMI

C--  THE FOLLOWING OUTPUT GIVES A SUMMARY ON THE SCREEN:
      IF (IFLAG.EQ.1) WRITE (*,612) Y(N+1),REALMA,REALMI
      IF (IFLAG.NE.1) WRITE (*,611) Y(N+1),REALMA
      RETURN

601   FORMAT (' UNSTABLE ***** , max.realpart of eigenvalue:',E15.7)
602   FORMAT (' STABLE , max. realpart of eigenvalue:',E15.7)
603   FORMAT (' Eigenvalue with realpart closest to zero:', 2E15.7)
604   FORMAT (' not all eigenvalues calculated. Analysis incomplete.')
605   FORMAT (' ****************************************************')
606   FORMAT (' **   Eigenvalue crosses imaginary axis,           **')
607   FORMAT (' **',48X,'**',
     1       /' **   CHANGE IN STABILITY,                         **',
     2       /' **   ',19('-'),26X,'**')
608   FORMAT (' Eigenvalues may have switched,',
     1        ' the following result is unsafe:',
     2       /' (You may try smaller steps  and/or  ',
     3        'check eigenvalues by inspection.)')
609   FORMAT (' **   close to parameter',E17.9,10X,'**')
610   FORMAT (' **   Periodic branch may emanate with initial     **',
     1       /' **            period:',E17.9,12X,'**')
611   FORMAT ('   PAR:',E13.6,' STABLE, maximum realpart is',E13.6)
612   FORMAT ('   PAR:',E13.6,' UNSTABLE, realpart max/min:',2E13.6)
613   FORMAT (' eigenvalue: ',2E19.11)
614   FORMAT (' ***  HOPF bifurcation near PAR.',E17.9,'  PERIOD',E17.9)
615   FORMAT ('   ',4E19.11)
9300  FORMAT (a1,i2,i4,2i2,a1,a2,18e20.13)
      END


      SUBROUTINE INTERA (N,Y,IFLAGSOL,FCN,EPS,JACOBI,WINDOW
     1           ,UNITLST,UNITOUT,UNIT1,UNIT3,UNIT31,UNIT4)
      DOUBLE PRECISION EPS,Y,YOLD,T1,S,ETACO,ZETACO,PARCO,
     1   CC,YFINAL,STEPMA,STEPMI,STEPMP,WEIGHT,WINDOW
     3   ,BOUND2,Y0BUFF,STEPBUFF,STEP,DELTA
      CHARACTER*9 Y0TEXT
      CHARACTER*2 JUDGEBUFF
      INTEGER HEAD1BUFF,UNITLST,UNITOUT,UNIT1,UNIT3,UNIT31,UNIT4
      PARAMETER (NDIM=12,NDIML=14,NDIMBUFF=6)
      DIMENSION Y(NDIML),YOLD(NDIML),WINDOW(2),Y0BUFF(NDIML,NDIMBUFF)
     1  ,Y0TEXT(NDIMBUFF),HEAD1BUFF(3,NDIMBUFF),STEPBUFF(NDIMBUFF)
     2  ,JUDGEBUFF(NDIMBUFF)
C   DIMENSION AT LEAST  Y(N+1),YOLD(N+1)
      EXTERNAL FCN
      COMMON /CONTRO/ ETACO,ZETACO,PARCO, NCO,N1CO,N2CO,
     1                KCO,IND2CO,NSECFLAG,NFLAG
C
C  This routine is a specific version of  B I F L A B  for  PACK-A.
C  It monitors in an INTERactive way the use of various options/tools
C  for investigating parameter-dependence of solutions of "algebraic"
C  equations. Routines of  PACK-A  are called.
C
C  Input and output are kept on a moderate level, some constants are
C  set to average values.
C
C     Feel free to change this routine according to your specific
C     example / hardware / software / needs.
C
C  Input data are read from file with (UNIT 3), results are written
C  to this file (solution-vector, parameter).

C     INPUT FROM CALLING SEQUENCE:
C     ---------------------------
C      N     : number of equations (=NCO)
C      FCN   : name of subroutine that defines the right-hand side.
C      EPS   : Required accuracy of solutions
C    ( JACOBI: =1 if the partial derivatives of FCN are provided
C              in routine FCNLIN.   ELSE =0  . not required in this version.)
C
C     INPUT FROM LOGICAL UNIT 3  (may be empty, see MODUS=0):
C     --------------------------
C      N+1  real numbers (will be updated):
C      Y(I),I=1,..,N : solution vector
C      Y(N+1)        : parameter value
C
C     INPUT FROM LOGICAL UNIT 1  (may be empty)
C     -------------------------
C      Options that can be used for continuation are stored in a file.
C
C     MODUS  CONTROLS THE SELECTED OPTION:
C     -----------------------------------
C     =0   entering data  into file 'START'
C
C     =1   continuation:  no restrictions, no testfunctions
C     =2   continuation:  setting of options via menu
C     =3   continuation following the current options
C
C  NOT INCLUDED IN THIS TEST VERSION:
C     =4   CHECK/RECONFIRM/RECALCULATE THE PRESENT SOLUTION
C                 IN DETAIL:
C            PARAMETRIZATION WITH RESPECT TO INDEX-TH INITIAL VALUE.
C            INCLUDES STABILITY ANALYSIS.
C            MAY BE HELPFUL ALSO FOR CALCULATING
C            A SOLUTION AFTER HAVING MANIPULATED THE INITIAL VALUES
C            IN FILE "START" (LOGICAL UNIT 3).
C     =6   Check the BUFFER
C
C     NEEDED SOFTWARE:  PACKA1,  LINEAR, NOLE, EISPACK
C     ---------------
C     Linearization of FCN (in SUBROUTINE FCNLIN) may be used for
C           MODUS=6
C
C     DEFAULT VALUES  are set to the standard values :
C     --------------  ITAIM=3, ITMAX=2*ITAIM, WEIGHT=1., IEXTRA=1 .
C             Changing these values enables to tune or to manipulate
C             the continuation.
C             ITAIM should match the prescribed accuracy, and the
C                   chosen "nonlinear solver".
C***********************************************************************

C  If necessary, the following default values for the continuation procedure
c  may be changed in this program:
      ITAIM=3
      ITMAX=2*ITAIM
      INTERACT=1
      WRITE (*,*) 'BIFPACK: entering routine INTERA.'
      WRITE (UNITLST,*) 'BIFPACK: entering routine INTERA.'
      WRITE (UNITLST,*) 'Requested accuracy is ',EPS
      ISTART=1
      CALL CHECKDATA (N,Y,IFLAGSOL,IFLAGPER,ISTART,INTERACT,
     1      NOB,NOS,NEXTMODUS,STEP,Y0BUFF,NUMBUFF,HEAD1BUFF,Y0TEXT,
     2      STEPBUFF,JUDGEBUFF,JBUFFER,UNIT3,UNIT31,UNIT4)
      ISTART=0
      STEPMI=EPS*10.d0
      WEIGHT=1.d0
      IEXTRA=1
      BOUND2=0.1d0

C***********************************************************************
      NCO=N
      N1=N+1
      N1CO=N1
      IF (N1.GT.NDIML) THEN
          WRITE (*,*) ' dimension NDIML should be ',N1
          WRITE (UNITLST,*) ' dimension NDIML should be ',N1
          STOP
          END IF
      LAST=0
      ETACO=0.d0
      ZETACO=0.d0
      IND2CO=1
      IRUN=0
      IFINAL=0
      MODUS=-1
      IFLAGPER=1
      ISTART=0

10    CONTINUE
      CALL CHECKDATA (N,Y,IFLAGSOL,IFLAGPER,ISTART,INTERACT,
     1      NOB,NOS,NEXTMODUS,STEP,Y0BUFF,NUMBUFF,HEAD1BUFF,Y0TEXT,
     2      STEPBUFF,JUDGEBUFF,JBUFFER,UNIT3,UNIT31,UNIT4)
      ISTART=0
      write (*,*) 'Warning: In case your START is no solution, better'
      write (*,*) ' start with ONE step to create a starting solution!'
CALL SETSTART (N,Y,IFLAGSOL,IFLAGPER,INTERACT)
      IF (MODUS.GE.0) GO TO 11


C*    READING / SETTING OF OPTIONS FOR CONTINUATION
C
200   REWIND UNIT1
      READ (UNIT1,*,END=201) S
      READ (UNIT1,*,END=202) INDEX
      READ (UNIT1,*,END=203) MAXCON
      READ (UNIT1,*,END=204) STEPMA
      READ (UNIT1,*,END=204) STEPMP
      READ (UNIT1,*,END=205) LEVEL
      READ (UNIT1,*,END=206) ISTAB
c      READ (UNIT1,*,END=207) WINDOW(1)
c      READ (UNIT1,*) WINDOW(2)
      REWIND UNIT1
      GO TO 11

201   WRITE (*,*) ' ENTER INITIAL STEP LENGTH:'
      READ (*,*) S

C   DETERMINATION OF INDEX:
202   CC=0.
          DO 2021 I=1,N
          IF (ABS(Y(I)).GT.CC) THEN
              CC=ABS(Y(I))
              INDEX=I
              END IF
2021      CONTINUE

C     SETTING SOME DEFAULT VALUES IN CASE THE OPTIONS-FILE IS NOT DEFINED.
203   MAXCON=20
204   STEPMA=1.d10
      STEPMP=STEPMA
205   LEVEL=2
206   ISTAB=0
c207   WINDOW(1)=-987654321.
c      WINDOW(2)=-WINDOW(1)
208   IFINAL=0

C***********************************************************************

11    WRITE (*,610)
      WRITE (*,*) ' Enter choice of MODUS:       '
      READ (*,*) MODUS
      IF (MODUS.LT.0) GO TO 999
      IF (MODUS.EQ.0) THEN
         ISTART=-1
         GO TO 10
         ENDIF
      IF (MODUS.EQ.2) GO TO 2
      IF (MODUS.EQ.6) go to 6
      IF (INDEX.LT.0) THEN
         WRITE (*,*) 'Index is negative indicating previous failure;'
         WRITE (*,*) 'BIFPACK changes MODUS to 2 ',
     1                'requesting check of control data.'
         MODUS=2
         GO TO 2
         ENDIF
      IF (MODUS.EQ.1) GO TO 1
      IF (MODUS.EQ.3) GO TO 3
      IF (MODUS.EQ.8) GO TO 8
      WRITE (*,*) ' Wrong modus.    '
      RETURN

CHECKING THE JACOBIAN :
C     IF (JACOBI.EQ.0) GO TO 99
C     CALL CHECKD (T1,Y,N,N,FCN,FCNLIN,ICALL,UNITLST)
C     GO TO 999


C************************************************************
C*
C*    A CONTINUATION WILL BE CARRIED OUT:
C*    (STEPSIZE NOT RESTRICTED)
C
1     CONTINUE
      STEPMA=999999.*Y(N1)
      STEPMP=STEPMA
      LEVEL=2
      WRITE (*,*) ' ENTER NUMBER OF CONTINUATION STEPS:'
      READ (*,*) MAXCON
      ISTAB=0
      BOUND2=9999.
      GO TO 3

C*********************************************************************
C*
C     EDITING OPTIONS:
C
2     WRITE (*,*) ' CURRENT STARTING DATA AT PARAMETER:',Y(N1)
      WRITE (*,613) S, INDEX, LEVEL, MAXCON, STEPMA, STEPMP, BOUND2,
     1              ISTAB,WINDOW(1), WINDOW(2)
      IF (IFINAL.EQ.0) THEN
                   WRITE (*,614)
              ELSE
                   WRITE (*,615) IFINAL, YFINAL
              END IF
      WRITE (*,616) NOB
299   WRITE (*,*) ' Enter  NUMBER of option you want to change,  or'
      WRITE (*,*) ' enter  0      to run the continuation,  or'
      WRITE (*,*) ' enter  -1     to return to menu:'
      READ (*,*) ICHOIC
      IF (ICHOIC.LT.0) GO TO 11
      IF (ICHOIC.EQ.0 .OR. ICHOIC.GT.9) GO TO 3
      IF (ICHOIC.EQ.1) GO TO 21
      IF (ICHOIC.EQ.2) GO TO 22
      IF (ICHOIC.EQ.3) GO TO 23
      IF (ICHOIC.EQ.4) GO TO 24
      IF (ICHOIC.EQ.5) GO TO 25
      IF (ICHOIC.EQ.6) GO TO 26
      IF (ICHOIC.EQ.7) GO TO 27
      IF (ICHOIC.EQ.8) GO TO 28
      IF (ICHOIC.EQ.9) GO TO 29

24    WRITE (*,*) ' Enter maximum NUMBER OF CONTINUATION STEPS:'
      READ (*,*) MAXCON
      GO TO 299
22    WRITE (*,*) ' Choose index for initial step:'
      WRITE (*,*) '      Enter 0 for contin. with respect to parameter,'
      WRITE (*,*) '      Enter K for fixing K-th component:'
      READ (*,*) INDEX
      GO TO 299
21    WRITE (*,*) ' Enter initial stepsize (in real):'
      READ (*,606) S
      GO TO 299
23    WRITE (*,*) ' Choose among three LEVELs of step control:'
      WRITE (*,*) ' Enter  2  for freely variable steps,'
      WRITE (*,*) '        1  for variable step but fixed index,'
      WRITE (*,*) '        0  for both steplength and index fixed:'
      READ (*,*) LEVEL
      GO TO 299
25    WRITE (*,*) ' Enter maximum step size for predictor (real,  >0 ):'
      READ (*,606) STEPMA
          IF (STEPMA.LE.0.) THEN
                      WRITE (*,*) ' must be positive        '
                      GO TO 25
                      END IF
      WRITE (*,*) ' Enter maximum step size for parameter (real, >0):'
      READ (*,*) STEPMP
      WRITE (*,*) ' Solution changes are limited to  10  percent.'  
      WRITE (*,*) ' Enter  10  if O.K., or  20  for 20%, or similar:'
      READ (*,*) IBOUND
      IF (IBOUND.GT.1) BOUND2=0.01*REAL(IBOUND)
      GO TO 299
26    WRITE (*,*) ' STABILITY ANALYSIS:  1  for YES,  0  for NO:'
      READ (*,*) ISTAB
      GO TO 299
27    WRITE (*,*) ' Resetting window for parameter:'

         WRITE (*,*) ' enter LEFT limit for parameter:'
         READ (*,*) WINDOW(1)
         WRITE (*,*) ' enter RIGHT limit for parameter:'
         READ (*,*) WINDOW(2)
      GO TO 299

28    WRITE (*,*) ' Which component has to reach target value?'
      WRITE (*,*) '   Enter index (or 0 if no target):        '
      READ (*,*) IFINAL
      IF (IFINAL.GT.0) THEN
            WRITE (*,611) IFINAL
            READ (*,*) YFINAL
            END IF
      GO TO 299

29    WRITE (*,*) ' Enter number of branch:'
      READ (*,*) NOB
      GO TO 299

C**************************************************************************

3     CONTINUE

      IF (INDEX.GT.N .OR. INDEX.LT.1) INDEX=N1
      IRUN=IRUN+1
      IPRINT=0
c                Write (*,*) ' detailed output: enter 1 , else 0:'
c                read (*,*) iprint
      write (*,*) 'nob ',nob
      WRITE (UNITLST,602) IRUN
      CALL CONTA (NOB,NOS,N,FCN,Y,YOLD,LAST,EPS,LEVEL,INDEX,S,
     1       STEPMA,STEPMI,STEPMP,MAXCON,ITAIM,ITMAX,WEIGHT,IEXTRA,
     2       IFINAL,YFINAL,ISTAB,IRUN,IPRINT,WINDOW,BOUND2,Y0BUFF,
     3       NUMBUFF,Y0TEXT,HEAD1BUFF,STEPBUFF,JUDGEBUFF
     4       ,UNITLST,UNITOUT,UNIT3,UNIT31)
      IF (INDEX.LT.0) GO TO 992
      WRITE (*,*) ' END OF CONTINUATION. exit CONTA'
      WRITE (UNITLST,*) ' END OF CONTINUATION. exit CONTA'
      Y(N+2)=0.d0
      GO TO 11

6     CONTINUE
      CALL CHECKDATA (N,Y,IFLAGSOL,IFLAGPER,ISTART,INTERACT,
     1      NOB,NOS,NEXTMODUS,STEP,Y0BUFF,NUMBUFF,HEAD1BUFF,Y0TEXT,
     2      STEPBUFF,JUDGEBUFF,JBUFFER,UNIT3,UNIT31,UNIT4)
      GO TO 11

8     CONTINUE
      write (*,*) 'Enter index for fixed continuation:'
      read (*,*) INDEX
      if (index.le.0) index=n+1
      write (*,*) 'Enter step length for fixed continuation:'
      read (*,*) DELTA
      write (*,*) 'Enter number of steps:'
      read (*,*) IHMAX
      write (*,*) 'Enter INDL for choice of test function:'
      read (*,*) INDL
      write (*,*) 'Enter INDK for choice of test function:'
      read (*,*) INDK
      call LOCATE (INDL,INDK,INDEX,DELTA,IHMAX,
     1                   N,UNITLST,UNITOUT,UNIT3,JACOBI)
      GO TO 11

c991   WRITE (*,*) ' not available without linearization'
c      RETURN
992   WRITE (*,*) ' NO SUCCESS'
      GO TO 11
999   CONTINUE
      REWIND UNIT1
      WRITE (UNIT1,*) S
      WRITE (UNIT1,*) INDEX
      WRITE (UNIT1,*) MAXCON
      WRITE (UNIT1,*) STEPMA
      WRITE (UNIT1,*) STEPMP
      WRITE (UNIT1,*) LEVEL
      WRITE (UNIT1,*) ISTAB
c      WRITE (UNIT1,*) WINDOW(1)
c      WRITE (UNIT1,*) WINDOW(2)
      RETURN


601   FORMAT (' Last STEPSIZE for K=',I2,' was',E10.3,10X)
602   FORMAT (/' **********  NEXT CONTINUATION RUN ( run',I3,
     1                 ' )  ************')
603   FORMAT (I1)
604   FORMAT (I2)
605   FORMAT ('   ',4E19.11)
606   FORMAT (E20.10)
607   FORMAT (' Enter component Y(',I2,'):                 ')
609   FORMAT (9X,'T',14X,'Y1',13X,'Y2',13X,'Y3             ')
610   FORMAT (' THIS IS THE MENU YOU CAN CHOOSE FROM:              ',
     1       /' -------------------------------------              ',
     2  /'    0  :  Entering STARTing data                         ',
     3  /'    1  :  CONTINUATION: no restrictions                  ',
     4  /'    2  :  CONTINUATION: setting of options, THEN RUN     ',
     5  /'    3  :  CONTINUATION: following previous options       ',
C     6  /'    4  :  just (re)calculate one solution                ',
     7  /'    6  :  check Buffer of special solutions           ',
     8  /'    8  :  branch switching with menu control             ',
     9  /'   -1  :  ending this session')
611   FORMAT (' Enter target for  ',I3,'-th component:')
612   FORMAT (' currently ',E14.4,' .LE. PAR .LE. ',E14.4)
613   FORMAT (' THESE ARE THE CURRENT OPTIONS:                      ',
     1  /' No: 1 ,  step length :  ',E14.5,
     2  /' No: 2 ,  index of fixed component :  ',I3,
     3  /' No: 3 ,  level of continuation :  ',I2,
     4  /' No: 4 ,  maximum number of continuation steps : ',I4,
     5  /' No: 5 ,  step bounds: Predictor:',E9.3,
     5                ' , Param.:',E9.3,' , rel.change:',F4.2,
     6  /' No: 6 ,  stability analysis (yes=1, no=0) :  ',I2,
     7  /' No: 7 ,  window :',E14.5,' .lt. parameter .lt. ',E14.5)
614   FORMAT (' No: 8 ,  target point:  no final value is set.')
615   FORMAT (' No: 8 ,  target point:  Y(',I3,') = ',E14.5)
616   FORMAT (' No: 9 ,  number of branch:',I4)
      END



      SUBROUTINE CHECKDATA (N,Y,IFLAGSOL,IFLAGPER,ISTART,INTERACT,
     1      NOB,NOS,NEXTMODUS,STEP,Y0BUFF,NUMBUFF,HEAD1BUFF,Y0TEXT,
     2      STEPBUFF,JUDGEBUFF,JBUFFER,UNIT3,UNIT31,UNIT4)
      DOUBLE PRECISION Y0BUFF,Y,SUMPERIOD,STEP,STEPBUFF,EXCHANGE
      INTEGER iflagsol,iflagper,istart,interact,numbuff,iempty
     1        ,nout,numbuff1,nextmodus,head1buff,exint,jbuffer
     2        ,UNIT3,UNIT31,UNIT4
      CHARACTER*9 Y0TEXT
      CHARACTER*2 JUDGED,JUDGEBUFF,EXCHAR
      PARAMETER (NDIM=12,NDIML=14,NDIMBUFF=6)
      DIMENSION Y0BUFF(NDIML,NDIMBUFF),Y0TEXT(NDIMBUFF),Y(NDIML)
     1  ,HEAD1BUFF(3,NDIMBUFF),STEPBUFF(NDIMBUFF),JUDGEBUFF(NDIMBUFF)
C***********************************************************************
C IFLAGPER: indicates periodic solution
C    = 0 : not set on input.
C    = 1 : require N+1 reals (N for vector Y, with value of
C                                parameter in position N+1)
C    = 2 : requires in addition on position N+2 a value
C                       for the period / final time
C  IFLAGPER1: related status of data in file START.
C  ISTART=1    upon starting BIFPACK.
C  ISTART=-1   to enter new data.
C***********************************************************************

      N1=N+1
      IF (ISTART.EQ.-1)  go to 3
Check data file START:
      IEMPTY=0
      REWIND UNIT3
      READ (UNIT3,606,END=10) (Y(K),K=1,N1)
      IFLAGPER1=1
      READ (UNIT3,606,END=11) Y(N+2)
      IFLAGPER1=2
c      IF (ABS(Y(N+2)).LT.0.000001d0) IFLAGPER1=1
      GO TO 11
10    IEMPTY=1
11    CONTINUE
      IF (IFLAGPER1.GT.IFLAGPER) IFLAGPER=IFLAGPER1

      REWIND UNIT31
      READ (UNIT31,*,END=12) NOB1
      IF (INTERACT.GT.0) NOB=NOB1
      GO TO 121
12    IF (INTERACT.EQ.0) THEN
          NOB=0
        ELSE
          WRITE (*,*) 'Enter number of branch:'
          READ (*,*) NOB
        ENDIF
      NOS=0
      JUDGED='  '
      NEXTMODUS=0
      STEP=0.
      REWIND UNIT31
      WRITE (UNIT31,*) NOB
      WRITE (UNIT31,*) NOS
      WRITE (UNIT31,*) JUDGED
      WRITE (UNIT31,*) NEXTMODUS
      WRITE (UNIT31,*) STEP
      GO TO 13
121   READ (UNIT31,*) NOS         
      READ (UNIT31,608) JUDGED         
      READ (UNIT31,*) NEXTMODUS
      READ (UNIT31,*) STEP         
13    CONTINUE

Check data file BUFFER in case ISTART=1:
      IF (ISTART.EQ.1) THEN
      REWIND UNIT4
      READ (UNIT4,*,END=27) NUMBUFF
      IF (NUMBUFF.LE.0) GO TO 28
      DO 24 J=1,NUMBUFF
      DO 22 I=1,N+2
22    READ (UNIT4,*,END=26) Y0BUFF(I,J)
      READ (UNIT4,605,END=26) Y0TEXT(J)
      READ (UNIT4,*) HEAD1BUFF(1,J)
      READ (UNIT4,*) HEAD1BUFF(2,J)
      READ (UNIT4,608) JUDGEBUFF(J)
      READ (UNIT4,*) HEAD1BUFF(3,J)
      READ (UNIT4,*) STEPBUFF(J)
24    CONTINUE
      GO TO 28
26    WRITE (*,*) 'Buffer is not complete.'
27    NUMBUFF=0 
28    CONTINUE
      ENDIF



      IF (IEMPTY.NE.1 .OR. NUMBUFF.NE.0) GO TO 4
         WRITE (*,*) 'Both files START and BUFFER are empty.'
         WRITE (*,*) 'Enter  0  for entering new data,'
         IF (IFLAGSOL.GT.0) THEN
           WRITE (*,*) 'enter  1  for taking data of routine SETPROB,'
           ENDIF
         WRITE (*,*) 'enter  -1 for leaving this routine.'
         READ (*,*) IFLAGSOL1
         IF (IFLAGSOL1.LT.0) RETURN
         IF (IFLAGSOL1.GT.0) GO TO 35

c   enter new data:
3        DO 303 I=1,N
         WRITE (*,607) I
303      READ (*,606) Y(I)
         WRITE (*,*) ' enter PARAMETER:'
         READ (*,606) Y(N1)
         IF (IFLAGPER.EQ.2) THEN
            WRITE (*,*) ' enter period or final time:'
            READ (*,606) Y(N+2)
            ENDIF

35       REWIND UNIT3
         WRITE (UNIT3,606) (Y(K),K=1,N+2)
         WRITE (*,*) 'data written to file START.'
         RETURN

C    show data:      

4     WRITE (*,*) 'Overview of stored data:'
      IF (NUMBUFF.GT.4) WRITE (*,*) NUMBUFF,' buffers are in use!'
c      IF (IFLAGPER.EQ.2) THEN
            NLAST=N+2
c         ELSE
c            NLAST=N+1
c         ENDIF
      WRITE (*,611)
      NUMBUFF1=NUMBUFF
      IF (NUMBUFF1.GT.4) NUMBUFF1=4
      WRITE (*,612) (k,k=1,4)
      WRITE (*,613) 
      NOUT=N
      IF (NOUT.GT.9) NOUT=9
      DO 41 I=1,NOUT
      IF (IEMPTY.EQ.0) WRITE (*,609) 
     1   'Y(',I,'):',Y(I),(Y0BUFF(I,K),K=1,NUMBUFF1)
      IF (IEMPTY.EQ.1) WRITE (*,6091) 
     1   'Y(',I,'):',(Y0BUFF(I,K),K=1,NUMBUFF1)
41    CONTINUE
      IF (N.GT.NOUT) WRITE (*,*) ' ( plus further ',N-NOUT,
     1                           ' components)'
      IF (IEMPTY.EQ.0) WRITE (*,610) 
     1   'PAR.:',Y(N+1),(Y0BUFF(N+1,K),K=1,NUMBUFF1)
      IF (IEMPTY.EQ.1) WRITE (*,6101) 
     1   'PAR.:',(Y0BUFF(N+1,K),K=1,NUMBUFF1)
      IF (IEMPTY.EQ.1) SUMPERIOD=0.d0
      IF (IEMPTY.EQ.0) SUMPERIOD=Y(N+2)*Y(N+2)
      DO 42 K=1,NUMBUFF
42    SUMPERIOD=SUMPERIOD+Y0BUFF(N+2,k)*Y0BUFF(N+2,k)
      IF (SUMPERIOD.GT.0.0001d0.and.IEMPTY.EQ.0)
     1  WRITE(*,610) 'PER.:',Y(N+2),(Y0BUFF(N+2,K),K=1,NUMBUFF1)
      IF (SUMPERIOD.GT.0.0001d0.and.IEMPTY.EQ.1)
     1  WRITE(*,6101) 'PER.:',(Y0BUFF(N+2,K),K=1,NUMBUFF1)
      WRITE (*,613)
      WRITE (*,614) NOB,NOS,
     1              (HEAD1BUFF(1,K),HEAD1BUFF(2,K),K=1,NUMBUFF1)
      WRITE (*,615) JUDGED,' ',
     1              (JUDGEBUFF(K),Y0TEXT(K),K=1,NUMBUFF1)
      WRITE (*,616) NEXTMODUS,(HEAD1BUFF(3,K),K=1,NUMBUFF1)
      WRITE (*,617) STEP,(STEPBUFF(K),K=1,NUMBUFF1)
      WRITE (*,613)

c-------------------------------------------------------------------
c  exchange of data between files START/HEAD and BUFFER:
51    IF (INTERACT.GT.0) THEN
      WRITE (*,*) 'STORE           =======>>   :  Writing the ',
     1           'current SOLUTION INTO BUFFER ?'
      WRITE (*,*) '     Enter  0  for no storing, -1 for exit,'
     1 ,'  else enter number of the buffer:'
      READ (*,*) JBUFFER
      IF (JBUFFER.GT.0) THEN
          if (jbuffer.gt.ndimbuff.or.jbuffer.gt.numbuff+1) then
             write (*,*) ' too large.'
             go to 51
             endif
          if (jbuffer.eq.numbuff+1) numbuff=jbuffer
          DO 52 I=1,N+1
52        Y0BUFF(I,jbuffer)=Y(I)
          IF (IFLAGPER.EQ.2.OR.NEXTMODUS.EQ.1) THEN
                Y0BUFF(N+2,jbuffer)=Y(N+2)
              ELSE
                Y0BUFF(N+2,jbuffer)=0.d0
              ENDIF
          Y0TEXT(jbuffer)=' '
          HEAD1BUFF(1,JBUFFER)=NOB
          HEAD1BUFF(2,JBUFFER)=NOS
          HEAD1BUFF(3,JBUFFER)=NEXTMODUS
          JUDGEBUFF(JBUFFER)=JUDGED
          STEPBUFF(JBUFFER)=STEP
          ENDIF
      ENDIF


      IF (INTERACT.EQ.0) JBUFFER=0
      IF (NUMBUFF.LE.0.OR.JBUFFER.LT.0) GO TO 99
      IF (INTERACT.GT.0) THEN
      WRITE (*,*) 'LOAD          <<==========    :  Writing a '
     1    ,'BUFFER INTO file START ?'
          WRITE (*,*) '     Enter  0  for no loading,'
     1    ,'  else enter number of buffer:'
          READ (*,*) JBUFFER
        ELSE
          JBUFFER=0          
          K=0
54        K=K+1
          IF (K.LE.NUMBUFF) THEN
checking for Hopf points:
            IF (HEAD1BUFF(3,K).EQ.2) JBUFFER=K 
            ENDIF
          IF (K.LT.NUMBUFF.AND.JBUFFER.EQ.0) GO TO 54
        ENDIF 
      IF (JBUFFER.GT.0) THEN
        IF (JBUFFER.GT.NUMBUFF) THEN
           write (*,*) 'too large, not available.'
           go to 9
           endif
        REWIND UNIT3
        DO 53 I=1,NLAST
        Y(I)=Y0BUFF(I,JBUFFER)
53      WRITE (UNIT3,606) Y0BUFF(I,JBUFFER)
        NOB1=HEAD1BUFF(1,JBUFFER)
        IF (INTERACT.GT.0) NOB=NOB1
        NOS=HEAD1BUFF(2,JBUFFER)
        JUDGED=JUDGEBUFF(JBUFFER)
        NEXTMODUS=HEAD1BUFF(3,JBUFFER)
        IF (INTERACT.EQ.0) HEAD1BUFF(3,JBUFFER)=0
        STEP=STEPBUFF(JBUFFER)
        REWIND UNIT31
        WRITE (UNIT31,*) NOB1
        WRITE (UNIT31,*) NOS
        WRITE (UNIT31,608) JUDGED
        WRITE (UNIT31,*) NEXTMODUS
        WRITE (UNIT31,*) STEP
        WRITE (*,*) 'Buffer no. ',jbuffer,
     1            ' written to files START and HEAD.'
        GO TO 9
        ENDIF

      IF (INTERACT.GT.0) THEN
        WRITE (*,*) 'EXCHANGE       <<========>>    :  Exchange '
     1      ,'a BUFFER with file START ?'
        WRITE (*,*) '     Enter  0  for no exchange'
     1 ,',  else enter number of buffer:'
        READ (*,*) JBUFFER
      IF (JBUFFER.GT.0) THEN
        IF (JBUFFER.GT.NUMBUFF) THEN
           write (*,*) 'too large, not available.'
           go to 9
           endif
        REWIND UNIT3
        DO 55 I=1,NLAST
        EXCHANGE=Y(I)
        Y(I)=Y0BUFF(I,JBUFFER)
        Y0BUFF(I,JBUFFER)=EXCHANGE
55      WRITE (UNIT3,606) Y(I)
        Y0TEXT(JBUFFER)=' '
        EXINT=NOB1
        NOB1=HEAD1BUFF(1,JBUFFER)
        IF (INTERACT.GT.0) NOB=NOB1 
        HEAD1BUFF(1,JBUFFER)=EXINT
        EXINT=NOS
        NOS=HEAD1BUFF(2,JBUFFER)
        HEAD1BUFF(2,JBUFFER)=EXINT
        EXCHAR=JUDGED
        JUDGED=JUDGEBUFF(JBUFFER)
        JUDGEBUFF(JBUFFER)=EXCHAR
        EXINT=NEXTMODUS
        NEXTMODUS=HEAD1BUFF(3,JBUFFER)
        HEAD1BUFF(3,JBUFFER)=EXINT
        EXCHANGE=STEP
        STEP=STEPBUFF(JBUFFER)
        STEPBUFF(JBUFFER)=EXCHANGE
        REWIND UNIT31
        WRITE (UNIT31,*) NOB1
        WRITE (UNIT31,*) NOS
        WRITE (UNIT31,608) JUDGED
        WRITE (UNIT31,*) NEXTMODUS
        WRITE (UNIT31,*) STEP
        WRITE (*,*) 'Buffer no. ',jbuffer,
     1            ' exchanged with files START / HEAD.'
        ENDIF
      ENDIF



c  write buffer-data to file BUFFER:

9     continue
      rewind UNIT4
      write (UNIT4,*) numbuff
      do 94 j=1,numbuff
      do 92 i=1,n+1
92    write (UNIT4,*) y0buff(i,j)
      write (UNIT4,*) y0buff(n+2,j)
      write (UNIT4,605) y0text(j)
      write (UNIT4,*) head1buff(1,j)
      write (UNIT4,*) head1buff(2,j)
      write (UNIT4,608) judgebuff(j)
      write (UNIT4,*) head1buff(3,j)
      write (UNIT4,*) stepbuff(j)
94    CONTINUE
      write (*,*) 'BUFFER stored to file.'
99    RETURN

605   FORMAT (A9)
606   FORMAT (E20.10)
607   FORMAT (' ENTER COMPONENT Y(',I2,'):                 ')
608   FORMAT (A2)
609   FORMAT (' ',A2,I1,A2,5E14.7)
6091  FORMAT (' ',A2,I1,A2,'              ',4E14.7)
610   FORMAT (' ',A5,5E14.7)
6101  FORMAT (' ',A5,'              ',4E14.7)
611   FORMAT ('        current     ',4('  buffer      '))
612   FORMAT ('        "START"     ',4('  no. ',i1,'       '))
613   FORMAT ('      ',5('  ------------'))
614   FORMAT (' NOB/S:',5(' ',I2,' / ',I4,'    '))
615   FORMAT (' type: ',5(' ',A2,', ',A9))
616   FORMAT (' nextMo:',5(' ',I1,'            '))
617   FORMAT (' step:',5E14.7)
      END
