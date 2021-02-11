
C*************************************
C**  THIS IS  packB1 of BIFPACK     **
C**  ---------------                **
C**    INCLUDING ROUTINES FOR A     **
C**    STABILITY AND BIFURCATION ANALYSIS
C**    OF SYSTEMS OF O.D.E. BOUNDARY-VALUE PROBLEMS
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

c====================================================================
C     packb1, LATEST UPDATE:  5 March 1996


      SUBROUTINE CONTB (NOB,NOS,N,M,T,FCN,BC,INTE,YN,YA,LAST,EPS,
     1   LEVEL,INDEX,STEP,STEPMA,STEPMI,STEPMP,MAXCON,ITAIM,ITMAX,
     2   WEIGHT,IEXTRA,IFINAL,YFINAL,IRUN,IPRINT,WINDOW,BOUND2,
     3   UNITLST,UNITOUT,UNIT3)

C  CONTINUATION ALGRITHM,
C  VERSION "B"

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
C    IRUN  : number of continuation run; only used for output.
C    IPRINT: =1 : data of step control are printed. (=0: omitted)
C    WINDOW: The two values in WINDOW define the parameter interval of
C            interest. The continuation stops if a solution is
C            calculated with a parameter value ouside the window.
C            WINDOW(1) : left bound,
C            WINDOW(2) : right bound.

C     HOW TO REACH TARGET POINTS  (IFINAL=0 FOR NO TARGET POINT):
C     ---------------------------
C           NO TARGET POINT REQUIRED: SET  IFINAL=0 , YFINAL IS DUMMY THEN.
C        IFINAL:   INDEX OF Y-COMPONENT OF TARGET POINT.
C        YFINAL :  TARGET VALUE FOR Y(IFINAL)
C            THERE IS NO TARGET POINT CHECK CARRIED OUT IN THE
C            FIRST CONTINUATION STEP.

C     OUTPUT:
C     ------   YN, YA, INDEX, STEP     are updated.
C
C     OUTPUT FOR GRAPHICS ETC. ON FILE 'DATOUT' (UNIT 8):
C     ---------------------------------------------------
C        A STANDARD RECORD OF EACH SOLUTION IS WRITTEN VIA FORMAT 9300

C    NEEDED SUBROUTINE FCN (T,Y,F)
C    -------------------------  defines the function  F(Y,PARAM)
C            FCN requires the following three statements:
C                   COMMON /CONTRO/ ETACO,ZETACO,PARCO, NCO,N1CO,N2CO,
C                  1                KCO,IND2CO,NSECFLAG,NFLAG
C                   PARAM=Y(N1CO)
C                   F(N1CO)=Y(KCO)-ZETACO*Y(IND2CO)-ETACO
C             (standard cont. takes ZETA=0 and IND2CO arbitrary.)
C
C    NEEDED SUBROUTINE BC(YA,YB,R)
C    -------------------------  defines boundary conditions
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

      EXTERNAL FCN,BC,INTE
      COMMON /CONTRO/ ETACO,ZETACO,PARCO, NCO,N1CO,N2CO,
     1                KCO,IND2CO,NSECFLAG,NFLAG
      DOUBLE PRECISION Y,YA,YN,EPS,STEP,STEPMA,STEPMI,WEIGHT,ETACO,
     1  ZETACO,PARCO,D,D1,D2,YNORM,ALPHA,PREDNO,SECNO,
     2  SECNO1,YFINAL,WINDOW,CURV,CURV1,CURVNE,SCALPR,SUM2,C1,C2,CY,
     3  CPAR,BETA,BOUND1,BOUND2,BOUND3,STEPMP,FACT,YTEST,T
      CHARACTER *2 TEXT
      INTEGER UNITLST,UNITOUT,UNIT3
      PARAMETER (NDIM=12,NDIML=14)
      DIMENSION Y(NDIML,7),YA(NDIML,7),YN(NDIML,7),T(7),
     1          WINDOW(2),YTEST(NDIML,7)
C     DIMENSIONS ALL  N+1
      WRITE (*,*) ' BIFPACK: entering routine CONTB.'
      WRITE (UNITLST,*) ' BIFPACK: entering routine CONTB.'
      KCO=INDEX
      IF (WEIGHT.LE.0.d0) WEIGHT=1.d0
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
	   D=YN(KCO,1)-YA(KCO,1)
           IF (0.d0.LE.D .AND. D.LT.EPS) ALPHA=1.d0
           IF (-EPS.LT.D .AND. D.LE.0.d0) ALPHA=-1.d0
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
      YNORM=ABS(YN(I,1))
      IF (YNORM.LT.0.0001) YNORM=1.d0
      D=ABS(YN(I,1)-YA(I,1)) /YNORM
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
      IF (D1.LT.EPS*0.1d0) THEN
         write (UNITLST,*) 'CONTB: step control failure. '
     1                   ,'D1 too small. stop'
         GO TO  9
         endif
C..  ADJUSTING STEPSIZE:
2     STEP=YN(KCO,1)-YA(KCO,1)
      ALPHA=DBLE(ITAIM)/DBLE(IT)
      if (alpha.gt.2.d0) alpha=2.d0
      if (alpha.lt.0.5d0) alpha=0.5d0
      IF (ICON.LT.4) GO TO 25
      CURVNE = CURV + (SECNO/SECNO1) * (CURV-CURV1)
      IF (CURVNE.LT.EPS) CURVNE=EPS
      IF (CURVNE.LT.CURV) CURVNE=0.25d0*(3.d0*CURVNE+CURV)
      BETA=BOUND1/(CURVNE*ALPHA*SECNO)
      IF (BETA.LT.1.d0) ALPHA=ALPHA*BETA
25    PREDNO=ALPHA*SECNO
      IF (SUM2.LE.BOUND3) SUM2=-1.d0
      IF (SUM2.GT.BOUND3) SUM2=PREDNO/SUM2
      IF (IPRINT.GT.0) THEN
              WRITE (UNITLST,*) ' PROPOSED STEP LENGTH FACTOR:',ALPHA
              WRITE (UNITLST,*) ' LENGTH OF FREE PREDICTOR:  ',PREDNO
              IF (SUM2.GT.0.) 
     1        WRITE (UNITLST,*) ' RELATIVE CHANGE IN SOLUTION:',SUM2
              END IF
      STEPAR=ALPHA*ABS((YN(N1,1)-YA(N1,1)))
      IF (STEPMI.LE.PREDNO.AND.PREDNO.LE.STEPMA.AND.
     1    STEPAR.LE.STEPMP.AND.SUM2.LE.BOUND2.AND.IFINAL.EQ.0) GO TO 295
      FACT=1.d0
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
	    Y(IFINAL,1)=YN(IFINAL,1)+ALPHA*(YN(IFINAL,1)-YA(IFINAL,1))
C..         CHECK FOR FINAL STEP TO REACH TARGET POINT
	    IF ( YN(IFINAL,1).LT.YFINAL.AND.YFINAL.LT.Y(IFINAL,1) .OR.
     1         Y(IFINAL,1).LT.YFINAL.AND.YFINAL.LT.YN(IFINAL,1) ) THEN
                KCO=IFINAL
		STEP=YN(IFINAL,1)-YA(IFINAL,1)
		FACT=ABS((YFINAL-YN(IFINAL,1))/
     1                (Y(IFINAL,1)-YN(IFINAL,1)))
                ALPHA=ALPHA*FACT
                PREDNO=PREDNO*FACT
                END IF                 
            END IF
      IF (IPRINT.GT.0 .AND. ABS(FACT-1.d0).GT.EPS) THEN
            WRITE (UNITLST,*) ' STEP LENGTH LIMITED TO:    ',PREDNO
            WRITE (UNITLST,*) ' FINAL STEP LENGTH FACTOR:',ALPHA
            END IF
295   STEP=STEP*ALPHA

3     IFAIL=0
C..  INITIAL GUESS:
      if (iprint.gt.0) write (UNITLST,*) ' final alpha=',alpha
31    DO 35 I=1,N1CO
      DO 35 J=1,M
      IF (LAST.EQ.1 .AND. IEXTRA.EQ.1)
     1     Y(I,j)=YN(I,j)+ALPHA*(YN(I,j)-YA(I,j))
      IF (LAST.EQ.0 .OR. IEXTRA.EQ.0) Y(I,j)=YN(I,j)
35    CONTINUE
      ETACO=YN(KCO,1)+STEP
      IT=ITMAX
      if (iprint.gt.0) write (UNITLST,6021) (y(i,1),i=1,n1)
      CALL SOLVER (FCN,BC,INTE,N1,M,T,Y,EPS,IT)

      IF (IT.GT.0) THEN
         REWIND UNIT3
	 WRITE (UNIT3,680) (Y(I,1),I=1,N1CO)
         WRITE (UNITLST,600) IT
         WRITE (UNITLST,6001) KCO,STEP
	 NOS=IRUN*100+ICON
	 NOS=ICON
         IF (LAST.EQ.0) GO TO 45
	 IF ( (Y(N1,1)-YN(N1,1))*(YN(N1,1)-YA(N1,1)) .LT.0.d0) THEN
C             approximation of turning point (parameter and y(IC) ) :
	      CPAR=0.d0
	      DO 43 IC=1,N
	      IF (ABS(YN(IC,1)-Y(IC,1)).LT.EPS.or.
     1            ABS(YA(IC,1)-Y(IC,1)).LT.EPS.or.
     2            ABS(YA(IC,1)-YN(IC,1)).LT.EPS) THEN
		  YTEST(IC,1)=0.d0
		  GO TO 43
		  ENDIF
	      C1=(YN(N1,1)-Y(N1,1))/(YN(IC,1)-Y(IC,1))
	      C2=(YA(N1,1)-Y(N1,1))/(YA(IC,1)-Y(IC,1))
	      C2=(C2-C1)/(YA(IC,1)-YN(IC,1))
	      CY=0.5*(Y(IC,1)+YN(IC,1)-C1/C2)
	      YTEST(IC,1)=CY
	      IF (CPAR.LT.EPS) CPAR=Y(N1,1)+
     1           (CY-Y(IC,1))*(C1+C2*(CY-YN(IC,1)))
43            CONTINUE
C  if required enter here iterative refinement
	      IF (ABS(CPAR).GE.EPS) THEN
		 IF (NSECFLAG.EQ.0) WRITE (UNITOUT,9300) 'r',NOB,0,N,
     1               1,'b','T ',(YTEST(IC,1),IC=1,N),CPAR
		 IF (NSECFLAG.EQ.1) WRITE (UNITOUT,9300) 'r',NOB,0,N,
     1               2,'b','T ',(YTEST(IC,1),IC=1,N),CPAR,PARCO
		 IC=1
                 WRITE (UNITLST,606) 
		 WRITE (UNITLST,604) CPAR,IC,YTEST(IC,1)
                 WRITE (UNITLST,606)
                 WRITE (*,606) 
		 WRITE (*,604) CPAR,IC,YTEST(IC,1)
                 WRITE (*,606) 
		 ENDIF
              END IF                                
45       SCALPR=0.d0
         SECNO=0.d0
         SUM2=0.d0
         SECNO1=0.d0         
         DO 5 I=1,N1
	 SUM2=SUM2+Y(I,1)*Y(I,1)
	 SECNO=SECNO+(Y(I,1)-YN(I,1))*(Y(I,1)-YN(I,1))
         IF (LAST.EQ.1) THEN
		SCALPR=SCALPR+(Y(I,1)-YN(I,1))*(YN(I,1)-YA(I,1))
		SECNO1=SECNO1+(YN(I,1)-YA(I,1))*(YN(I,1)-YA(I,1))
                END IF
	 DO 5 J=1,M
	 YA(I,j)=YN(I,j)
5        YN(I,j)=Y(I,j)
         SECNO=SQRT(SECNO)
         SUM2=SQRT(SUM2)
         IF (IPRINT.GT.0) THEN
	       WRITE (UNITLST,*) ' LENGTH OF SECANT:     ', SECNO
               WRITE (UNITLST,*) ' L-2 NORM OF SOLUTION: ',SUM2
               END IF
         IF (LAST.EQ.1) THEN
               SECNO1=SQRT(SECNO1)
               SCALPR=ABS(SCALPR)/(SECNO*SECNO1)
               IF (1.-SCALPR .GT.10.*EPS) THEN
                      SCALPR=ACOS(SCALPR)
                      ELSE
                      SCALPR=0.d0
                      END IF
               IF (ICON.GT.2) CURV1=CURV
               CURV=SCALPR/SECNO
               IF (IPRINT.GT.0) WRITE (11,*) ' CURVATURE: ', CURV
               END IF
         LAST=1

	 write (UNITLST,601) y(n1,1)
	 write (*,605) nos,kco,step,it,y(n1,1)
	 write (*,602) (y(i,1),i=1,n)
	 write (UNITLST,602) (y(i,1),i=1,n)
	 TEXT='  '
	 IF (NSECFLAG.EQ.0) WRITE (UNITOUT,9300) 'R',NOB,NOS,N,1,'b',
     1                  TEXT,(Y(I,1),I=1,N),Y(N1CO,1)
	 IF (NSECFLAG.EQ.1) WRITE (UNITOUT,9300) 'R',NOB,NOS,N,2,'b',
     1                  TEXT,(Y(I,1),I=1,N),Y(N1CO,1),PARCO
	 IF (IFINAL.EQ.0) GO TO 51
	 IF (ABS(YFINAL-YN(IFINAL,1)).LT.EPS) GO TO 99
51       IF (YN(N1,1).LT.WINDOW(1) .OR. YN(N1,1).GT.WINDOW(2)) THEN
               WRITE (*,*) ' Continuation has left the window.'
               WRITE (UNITLST,*) ' Continuation has left the window.'
               GO TO 99
               END IF
       IF (ICON.LT.MAXCON) GO TO 1

C************************************************************************
C..  FAIL- AND END-CONDITIONS:                                          *
         GO TO 99                                                       *
      ELSE
         WRITE (11,6002) IT
         WRITE (UNITLST,6001) KCO,STEP
	 IF (LEVEL.EQ.0) then
	     write (11,*) 'CONTB: failure of step control'
	     GO TO 9
	     ENDIF
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
600   FORMAT (' SOLUTION OBTAINED AFTER ',I3,' ITERATIONS.')
6001  FORMAT (' STEP CONTROL: INDEX OF FIXED COMPONENT WAS ',I3,
     2        ' ,   STEPSIZE: ',E12.6)
6002  FORMAT ('  FAILURE ',I3) 
601   FORMAT (' ****  PARAMETER VALUE: ',E19.11,'   ****',
     1       /                   27X,14('-'))
680   FORMAT (E20.10)
9300  FORMAT (A1,I2,I4,2I2,A1,A2,18E20.13)
602   FORMAT ('   sol. :',4E17.9)
6021  FORMAT (' guess:',4E17.9)
603   FORMAT (/' CONTINUATION STEP NO.:',I4)
604   FORMAT (' ***  TURNING OF BRANCH AT PARAMETER: ',E13.6,
     1        ' , Y(',I3,')= ',E13.6)
605   FORMAT (' step:',i4,'  size:',i2,'/',e10.3,' Iter:',i2,
     1       '         parameter:',e17.9)
606   FORMAT ('**************************************************')
      END
