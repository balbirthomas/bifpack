C*************************************
C**  THIS IS  packA2 of BIFPACK     **
C**  ---------------                **
C**    INCLUDING ROUTINES FOR A     **
C**    BIFURCATION ANALYSIS         **
C**    OF SYSTEMS OF 'ALGEBRAIC'    **
C**    EQUATIONS                    **
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

c=========================================================================

C     packa2, latest update: 8. June 1999



C   REMARKS CONCERNING STORAGE / REDUCING BIFPACK
C   =============================================
C      IN CASE ONE DISPENSES WITH SOME OPTIONS, THE STORAGE
C      REQUIREMENTS / PACKAGE LENGTH CAN BE SIGNIFICANTLY REDUCED.
C      THIS MEANS IN PARTICULAR:
C
C      TESTA : TWO (N*N) MATRICES ARE REDUNDANT IF THE DERIVATIVE OF
C              THE TESTFUNCTION IS NOT CALCULATED.
C
C      SWITA : THE ATTACHED ROUTINE  FCNLAR  CAN BE DROPPED IF ONE
C              DISPENSES WITH THE EXPENSIVE OPTIONS WHICH NEED TO
C              DOUBLE THE NUMBER OF EQUATIONS. THIS SAVES STORAGE
C              WORTH MORE THAN ONE (N*N) MATRIX.

      SUBROUTINE TESTA (N,Y,J,INDL,INDK,TAU,N1,NDER,TAUDER,UNITLST)
      DOUBLE PRECISION Z,ZLK,Y,RS,H,D,WORK,TDUMMY
     1  , Y2,TAU,TAUDER,TDZ,DELTA,ETA,APPR,EPMACH
      PARAMETER (NDIM=12,NDIML=14)
      DIMENSION Z(NDIML,NDIML),ZLK(NDIML,NDIML),Y(NDIML),RS(NDIML)
     1        ,H(NDIML),D(NDIML),WORK(NDIML)
C    (DIMENSIONS SEE BELOW)
      INTEGER PIVOT(NDIML),UNITLST

C   THE FOLLOWING VARIABLES ARE NEEDED FOR THE OPTION  NDER>0 .
C   (THESE ARRAYS MAY BE PURGED AS WELL AS THE SECOND PART OF TESTA.)
      DOUBLE PRECISION Z2(NDIML,NDIML),ZH(NDIML,NDIML),H2(NDIML),
     1                 YH(NDIML)
      INTEGER PIV2(NDIML)

      COMMON /MATRIX/ Z
      NROW=NDIML
C
C--  ALGORITHM FOR CALCULATING A BRANCHING-TESTFUNCTION,
C    ITS DERIVATIVE, A GUESS OF THE NEXT BIFURCATION POINT,
C    AND AN APPROXIMATION OF THE SOLUTION OF THE LINEARIZATION
C    VERSION OF FEBR.1980 (ISNM 54)
C-----------------------------------------------------------------------------
C  INPUT:
C  -----
C     N   NUMBER OF NONLINEAR EQUATIONS
C     Y   VECTOR CONTAINING A SOLUTION, Y(N+1)=PARAMETER.
C         IN CASE N1=N+1 NEED DIMENSION Y(2*N+1)
C     J   .GT.0:  SOME OUTPUT IS PRINTED
C     INDL, INDK :  INDICES BETWEEN 1 AND N.
C         THE CHOICE (RATHER ARBITRARY) DETERMINES WHICH TESTFUNCTION
C         IS TO BE CALCULATED.
C     N1  CONTROL PARAMETER, CONTROLS THE APPROXIMATION OF A SOLUTION
C         TO THE LINEARIZATION.
C         N1=N+1 :  APPROXIMATION IS CALCULATED. THESE VALUES WILL
C                   BE STORED IN Y(K), K=N+2,..,2*N+1 .
C         N1.LT.N : NO APPROXIMATION
C    NDER IS CONTROL PARAMETER FOR THE CALCULATION OF THE DERIVATIVE
C         OF THE TESTFUNCTION.
C         .LE.0 : NO DERIVATIVE
C         .GT.0 : DERIVATIVE IS CALCULATED. SUBROUTINE FCNLIN IS NEEDED
C                 FOR PROVIDING THE  N*(N+1) PARTIAL DERIVATIVES OF THE
C                 SYSTEM WITH RESPECT TO Y AND THE PARAMETER.
C
C  OUTPUT
C  ------
C     TAU     VALUE OF THE TESTFUNCTION
C     TAUDER  DERIVATIVE OF THE TESTFUNCTION (IF NDER.GT.0)
C     Y(N+2,...)  SEE ABOVE
C     J       =-1  IF MATRIX IS SINGULAR
C
C  FURTHER REMARKS
C  ---------------
C     IN THIS VERSION THE FUNCTIONAL MATRIX  Z  EVALUATED AT THE
C     SOLUTION IS PROVIDED BY THE STATEMENT
C                  COMMON /MATRIX/ Z
C     THIS MAKES THE ITERATION MATRIX OF NEWTON'S METHOD AVAILABLE.
C     THE COMMON STATEMENT CAN BE REPLACED BY
C                  CALL FCNLIN (T,Y,Z)
C     LINEAR EQUATIONS IN TESTA ARE SOLVED BY THE DECOMPOSITION ROUTINES
C                  BGDECO  AND BGSOLV
C     DUE TO BUSINGER AND GOLUB.
C
C  DIMENSIONS
C  ----------
C     FOR DIMENSION OF  Y  SEE ABOVE.
C     ZLK(N,N),Z(N,N+1),Z2(N,N),ZH(N,N+1),YN(N+1) ;
C     ALL THE OTHER VECTORS NEED DIMENSION  N .
C
C===================================================================

      IPR=J
      EPMACH=1.d-12
      TAU=0.d0
      DO 11 I=1,N
      DO 11 J=1,N
11    ZLK(I,J)=Z(I,J)
      IF (INDK.LE.0.OR.INDL.LE.0) GO TO 2

      DO 1 J=1,N
      ZLK(INDL,J)=0.d0
      RS(J)=0.d0
1     RS(INDL)=1.d0
      ZLK(INDL,INDK)=1.d0

2     CALL BGDECO (N,N,ZLK,NROW,D,PIVOT,EPMACH,ISING)
      IF (ISING.LT.0) GO TO 99
      IF (INDK.LE.0.OR.INDL.LE.0) GO TO 22
      CALL BGSOLV (N,N,ZLK,NROW,D,PIVOT,RS,H,WORK)

22    TDZ=D(1)
      DO 27 J=2,N
27    TDZ=TDZ*D(J)
      IF (INDK.GT.0.AND.INDL.GT.0) GO TO 45
      IF (IPR.GT.0) WRITE (UNITLST,602) INDL,INDK,TAU,TDZ
      GO TO 98

45    DO 4 J=1,N
4     TAU=TAU+Z(INDL,J)*H(J)
      IF(IPR.GT.0 ) WRITE(UNITLST,602) INDL,INDK,TAU,TDZ

      IF (N1.LT.N) GO TO 49
C     INITIAL GUESS FOR LINEARIZATION:
      DO 48 I=1,N
      I1=N1+I
48    Y(I1)=H(I)
49    IF (NDER.LE.0) GO TO 98

C============================================================
C    The following part calculates the derivative of the testfunction if
C    NDER.GT.0
C    This part, together with the dimension statement for Z2, ZH, ...
C    may be purged if the space- and time consuming derivative is not required.

C--  NUMERICAL DIFFERENTIATION:
      NP1=N+1
      ETA=1.E-4
      DO 51 I=1,N
      RS(I)=-Z(I,NP1)
      DO 51 J=1,N
      Z2(I,J)=0.d0
51    ZH(I,J)=Z(I,J)
      CALL BGDECO (N,N,ZH,NROW,H2,PIV2,EPMACH,ISING)
      IF (ISING.LT.0) GO TO 99
      CALL BGSOLV (N,N,ZH,NROW,H2,PIV2,RS,YH,WORK)
      YH(NP1)=1.d0

      DO 58 K=1,NP1
      Y2=Y(K)
      DELTA=ETA*Y2
      IF (ABS(DELTA).LT.1.d-15) DELTA=ETA
      Y(K)=Y2+DELTA
      CALL FCNLIN (TDUMMY,Y,ZH)
      DO 54 I=1,N
      DO 54 J=1,N
54    Z2(I,J)=Z2(I,J)+(ZH(I,J)-Z(I,J))*YH(K)/DELTA
58    Y(K)=Y2

C--  DERIVATIVE OF TESTFUNCTION:
      DO 72 I=1,N
      RS(I)=0.d0
      IF (I.EQ.INDL) GO TO 72
      DO 71 J=1,N
71    RS(I)=RS(I)-Z2(I,J)*H(J)
72    CONTINUE
      CALL BGSOLV (N,N,ZLK,NROW,D,PIVOT,RS,H2,WORK)
      TAUDER=0.d0
      DO 75 J=1,N
75    TAUDER=TAUDER+Z2(INDL,J)*H(J)+Z(INDL,J)*H2(J)
      APPR=Y(NP1)-TAU/TAUDER
      IF (IPR.GT.0) WRITE (UNITLST,604) TAUDER,APPR

C============================================================

98    J=1
      RETURN
99    J=-1
      WRITE (UNITLST,601)
      RETURN
601   FORMAT (' MATRIX SINGULAR')
602   FORMAT (' BRANCHING-FUNCTION  (',I2,'/',I2,')',E13.5,
     1 '   DET.',E13.5)
604   FORMAT (' DERIVATIVE OF BR. FUNCTION:',E14.7,
     1    '  GUESS OF BRANCH PT.:',E14.7)
      END


      SUBROUTINE SWITA (N,Y1,TAU1,Y2,TAU2,MODEB,METHOD,
     1                  UNITLST,UNITOUT,UNIT3)
      DOUBLE PRECISION Y1,Y2,Y0,YEM,TAU1,TAU2,FAK1,FAK2,EPMACH,PARCO
     1,COMP,COMP1,DELTA,DELTA1,DELTA2,ETACO,ZETACO,EPS,EPSB,EPSB1
     2 ,SECPARINC,SECPAR1 
      INTEGER UNITLST,UNITOUT,UNIT3
      PARAMETER (NDIM=12,NDIML=14)
      DIMENSION Y1(NDIML),Y2(NDIML),Y0(NDIML),YEM(NDIML)
C  DIMENSIONS:  Y1(N+N+1),Y2(N+N+1),Y0(N+N+1),YEM(N+1)

C----------------------------------------------------------------------------
C    ALGORITHM FOR BRANCH SWITCHING
C    AND FOR CALCULATION OF TURNING POINTS AND BIFURCATION POINTS
C    IN SYSTEMS OF NONLINEAR EQUATIONS.
C    REFERENCE: R.SEYDEL: ISNM  69 (1983/84) P.181-191
C-----------------------------------------------------------------------------

C--  NEEDED USER-SUPPLIED SUBROUTINES:
C    --------------------------------
C         FCN(T,Y,F)   EVALUATES FUNCTION VECTOR F(Y,PAR).
C                      PARAMETER IS CONTAINED IN  Y(N+1), T IS DUMMY.
C         NOLE(NDIML,Y,FCN,EPS,ITMAX,IT)   :
C                   SOLVER OF NONLINEAR EQUATIONS FCN=F(Y).
C                   USE, FOR EXAMPLE, the NEWTON METHOD.
C                   SIZE:    NDIML=N+1   OR   NDIML=N+N+1 ,  SEE BELOW.
C                   Y:       INITIAL GUESS AND SOLUTION VECTOR
C                   EPS:     RELATIVE ACCURACY.
C                   IT.GT.0: PRINTING DESIRED (INPUT),
C                            NUMBER OF ITERATIONS (OUTPUT)
C                   IT.LT.0: NO PRINTING (INPUT), FAILURE (OUTPUT).
C    IN CASE A BRANCH POINT IS TO BE CALCULATED (MODEB.GT.0):
C         FCNLIN(T,Y,A)  CALCULATES THE PARTIAL DERIVATIVES OF FUNCTION FCN:
C                   A(I,J)=DF(I)/DY(J)  .
C                   PARAMETER IS CONTAINED IN Y(N+1), T IS DUMMY.
C-----------------------------------------------------------------------------

C--  SUBROUTINES USED ALONG WITH SWITA:
C    ------------------------------------
C         TESTA      FOR EVALUATING BRANCHING-TESTFUNCTION
C                    AND PROVIDING INITIAL GUESS FOR LINEARIZATION.
C         FCNLAR     INTERNAL ROUTINE, USED BY SWITA FOR THE
C                    OPTIONS:  METHOD=1 , MODEB>0
C-----------------------------------------------------------------------------

C--  INPUT VARIABLES:
C    ---------------
C         N         NUMBER OF VARIABLES.
C         Y1        ONE SOLUTION ON BRANCH, OUTPUT OF TESTA  (WITH N1=N+1).
C         TAU1      CORRESPONDING VALUE OF TESTFUNCTION, OUTPUT OF TESTA.
C         Y2, TAU2  ANOTHER SOLUTION ON THE BRANCH WITH TESTFUNCTION.
C
C                   BASIC ASSUMPTION:
C                   Y1 AND Y2 ARE SITUATED ON A BRANCH THAT CAN BE
C                   PARAMETRIZED BY THE PARAMETER.
C
C         METHOD    CONTROLS THE METHOD FOR BRANCH SWITCHING:
C                   =1  METHOD OF PARALLEL COMPUTATION,
C                       NEEDS THE DOUBLE DIMENSION:  NDIML=N+N+1
C                   =2  'SIMPLE' METHOD
C                   =3  EXPLOITS SYMMETRY-BREAKING.
C                       2 AND 3 NEED SIMPLE DIMENSION OF NDIML=N+1.
C                       IF A COMBINATION SYMMETRY/ASYMMETRY CANNOT BE
C                       FOUND,  METHOD 2  IS PERFORMED.
C
C         MODEB     .GT.0  BRANCHING SYSTEM IS SOLVED TO CALCULATE A
C                       BIFURCATION OR TURNING POINT.
C                       VALUE REQUESTED FOR MODEB:
C                       MODEB = K ,  K IS THE INDEX OF THE TESTA-INPUT.
C                       SUBROUTINE FCNLIN IS NEEDED.
C                       THE DOUBLE DIMENSION IS NEEDED: NDIML=N+N+1 .
C                   =0  BIFURCATION POINT IS OBTAINED BY INTERPOLATION.
C                       FCNLIN IS NOT NEEDED.
C     RECOMMENDED STANDARD USE:
C     ------------------------
C              METHOD=3 , MODEB=0  for branch switching
C              (NEEDS NEWTON'S METHOD FOR ONLY THE SIMPLE
C               DIMENSION  NDIML=N+1 ),
C              METHOD=0 , MODEB=K  for calculating a turning point
C                (output on file DATOUT) 
C
C-----------------------------------------------------------------------------


      EXTERNAL FCNLAR,FCN,FCNLIN
      COMMON /CONTRO/ ETACO,ZETACO,PARCO, NCO,N1CO,N2CO,
     1                KCO,IND2CO,NSECFLAG,NFLAG
      N1=N+1
      EPMACH=1.d-12
      NCO=N
      N1CO=N1
      NL=N1+N
      NOB=0
      WRITE (UNITLST,*) 'BIFPACK: entering SWITA.'
      WRITE (*,*) 'BIFPACK: entering SWITA.'
C    The following values can be changed if necessary:
C    EPS  is error tolerance for NOLE,
C    MAXSOL  is number of emanating solutions to be calculated,
C    DELTA  is distance to emanating solution.
      EPS=1.d-4
      DELTA=0.02d0
      MAXSOL =1


C--  BRANCH POINT Y0 BY INTERPOLATION:
      IF (TAU1*TAU2.LT.0.d0) GO TO 12
      FAK1=ABS(TAU1/TAU2)
      IF (FAK1.GT.1.2d0 .OR. FAK1.LT.0.8d0) GO TO 12
      FAK1=1.d0
      FAK2=0.d0
      GO TO 15
12    FAK2=1./(TAU1-TAU2)
      FAK1=TAU1*FAK2
      FAK2=TAU2*FAK2
15    DO 18 I=1,NL
18    Y0(I)=FAK1*Y2(I)-FAK2*Y1(I)
      WRITE (UNITLST,601)
      WRITE (*,601)
      WRITE (UNITLST,606) Y0(N1)
      WRITE (*,606) Y0(N1)
      IF (NSECFLAG.EQ.0) WRITE (UNITOUT,9300) 
     1     'r',-1,0,n,1,' ','B ',(Y0(I),I=1,N1)
      IF (NSECFLAG.EQ.1) WRITE (UNITOUT,9300) 
     1     'r',-1,0,n,2,' ','B ',(Y0(I),I=1,N1),PARCO
      WRITE (UNITLST,600) (Y0(I),I=1,N)
      WRITE (*,600) (Y0(I),I=1,N)
      WRITE (UNITLST,608)
      WRITE (UNITLST,600) (Y0(N1+I),I=1,N)
      IF (MODEB.EQ.0) GO TO 3


C--  SOLVING THE BRANCHING SYSTEM:                                   
      EPSB=MAX(EPS,EPMACH**0.3)
      ITMAX=10
      NFLAG=-1
      KCO=MODEB
      NSECFLAG1=0
      NSECFLAG2=0
      SECPAR1=PARCO
2     J=-1
      NSECFLAG1=NSECFLAG1+1
      CALL NOLE (NL,Y0,FCNLAR,EPSB,ITMAX,J)
      IF (J.GT.0) THEN
         WRITE (UNITLST,*) ' '
         WRITE (UNITLST,602) J
         WRITE (UNITLST,606) Y0(n+1)
         WRITE (UNITLST,600) (Y0(I),I=1,N)
c        WRITE (UNITLST,608)
c        WRITE (UNITLST,600) (Y0(N1+I),I=1,N)
         if (nsecflag1.eq.1) then
           WRITE (*,602) J
           WRITE (*,606) Y0(N+1)
           WRITE (*,600) (Y0(I),I=1,N)
           if (nsecflag.eq.1.and.nsecflag1.eq.1) then
              write (*,*) 'Control of second parameter.'
              write (*,*) 'current value: ',PARCO 
              write (*,*) 'Vary second parameter? ',
     1                   '  (1=yes, 0=no ) ENTER:'
              read (*,*) nsecflag2
              do 292 i=1,NL
292           yem(i)=y0(i)
              if (nsecflag2.eq.1) then 
                write (*,*) 'Enter the increment in second parameter:'
                read (*,*) secparinc
                write (*,*) 'enter the number of solutions that are to' 
                write (*,*) ' be calculated along bifurcation curve:'
                read (*,*) nsecflag3
                if (nsecflag3.lt.0) nsecflag3=5 
                endif 
            endif
         endif 
         IF (NSECFLAG2.EQ.1) WRITE (*,613) Y0(N+1),PARCO
         IF (NSECFLAG2.EQ.1) WRITE (UNITLST,613) Y0(N+1),PARCO
         if (nsecflag.ne.1) WRITE (UNITOUT,9300)
     1         'R',NOB,0,N,1,' ','B ',(Y0(I),I=1,N+1)
         if (nsecflag.eq.1.AND.nsecflag1.eq.1) WRITE (UNITOUT,9300)
     1         'R',NOB,0,N,2,' ','B ',(Y0(I),I=1,N+1),PARCO
         if (nsecflag2.eq.1) WRITE (UNITOUT,9300)
     1         'R',0,0,N,2,' ','B ',(Y0(I),I=1,N+1),PARCO
         if (nsecflag1.lt.nsecflag3
     1       .and.nsecflag2.eq.1)  then
             PARCO=PARCO+secparinc
             go to 2
             endif
         END IF
      if (nsecflag2.eq.1) then
           if (nsecflag1.lt.nsecflag3) then
              write (*,*) 'The bifurcation curve could not '
     1            ,'be traced farther.'
              write (UNITLST,*) 'The bifurcation curve could not '
     1            ,'be traced farther.'
           endif
         PARCO=secpar1
         do 293 i=1,NL 
293      y0(i)=yem(i)
         endif
      if (nsecflag1.gt.1) write (*,*) 'second parameter'
     1   ,' retained to the old value ',PARCO
      if (nsecflag1.gt.1) write (UNITLST,*) 'second parameter'
     1   ,' retained to the old value ',PARCO

C=================================================================

C--  INITIAL GUESS OF EMANATING SOLUTION YEM
3     IF (METHOD.LE.0) GO TO 8
      ITMAX=15
      NFLAG=METHOD

C--  DETERMINATION OF INDEX  K=KCO :
      COMP=0.d0
      DO 305 I=1,N
      COMP1=ABS(Y0(N1+I))
      IF (ABS(Y0(I)).GT.EPS) COMP1=COMP1/ABS(Y0(I))
      IF (COMP1.GT.COMP) KCO=I
305   IF (COMP1.GT.COMP) COMP=COMP1
      DELTA1=ABS(Y0(KCO))*DELTA
      IF (DELTA1.LT.DELTA) DELTA1=DELTA

31    KNUM=0
      DELTA2=0.d0
32    KNUM=KNUM+1
      DELTA2=DELTA2+DELTA1

      DO 33 I=1,N
33    YEM(I)=Y0(I)+DELTA2*Y0(N1+I)
      YEM(N1)=Y0(N1)
      IF (NFLAG.EQ.0 .OR. NFLAG.GT.3) GO TO 8
      IF (NFLAG.EQ.2) GO TO 5
      IF (NFLAG.EQ.3) GO TO 6


C--  PARALLEL COMPUTATION:
      DO 4 I=1,N
4     Y0(N1+I)=YEM(I)
      ETACO=DELTA2
      J=-1
      CALL NOLE (NL,Y0,FCNLAR,EPS,ITMAX,J)
      GO TO 69


C--  SIMPLE COMPUTATION
5     NFLAG=2
      ZETACO=0.d0
      IND2CO=1
      ETACO=YEM(KCO)
      J=-1
      CALL NOLE (N1,YEM,FCN,EPS,ITMAX,J)
      GO TO 69


C--  EXPLOITING SYMMETRY BREAKING:
C     CHECKING OF SYMMETRY:
6     I1=0
      EPSB=EPS*5.d0
61    I1=I1+1
      IF (I1.GT.N) GO TO 5
      EPSB1=EPSB*ABS(Y0(I1))
      IF (EPSB1.LT.EPS) EPSB1=EPSB

      I2=0
615   I2=I2+1
      IF (I2.GT.N) GO TO 61
      IF (ABS(Y0(I1)-Y0(I2)).GT.EPSB1) GO TO 615
C     SYMMETRY TYPE  Y(I1)=Y(I2)  IS FOUND

C     NOW CHECKING OF ASYMMETRY:
      ETACO=YEM(I1)-YEM(I2)
      IF (ABS(ETACO).LE.10.*EPSB) GO TO 615

C     COMBINATION SYMMETRY/ASYMMETRY IS FOUND
      J=-1
      ZETACO=1.d0
      KCO=I1
      IND2CO=I2
      CALL NOLE (N1,YEM,FCN,EPS,ITMAX,J)



C--  OUTPUT OF RESULTS:
69    IF (J.GT.0) THEN
            WRITE (UNITLST,603) J
            WRITE (*,603) J
            END IF
      IF (J.LT.0) THEN
            WRITE (UNITLST,*) 'No convergence, with flag ',J
            WRITE (*,*) 'No convergence, with flag ',J
            END IF
      IF (NFLAG.LE.2) WRITE (UNITLST,604) NFLAG,KCO
      IF (NFLAG.EQ.3) WRITE (UNITLST,605) NFLAG,KCO,IND2CO
      IF (NFLAG.EQ.1) THEN
           WRITE (UNITLST,606) Y0(N1)
           WRITE (UNITLST,600) (Y0(I),I=N+2,NL)
           WRITE (UNITLST,600) (Y0(I),I=1,N)
           GO TO 7
           END IF
      WRITE (UNITLST,606) YEM(N1)
      WRITE (*,606) YEM(N1)
      WRITE (UNIT3,680) (YEM(K),K=1,N1)
      IF (J.GT.0.AND.NSECFLAG.EQ.0) WRITE (UNITOUT,9300) 
     1      'R',-1,KNUM,n,1,' ','E ',(YEM(I),I=1,N+1)
      IF (J.GT.0.AND.NSECFLAG.EQ.1) WRITE (UNITOUT,9300) 
     1      'R',-1,KNUM,n,2,' ','E ',(YEM(I),I=1,N+1),PARCO
      WRITE (UNITLST,600) (YEM(I),I=1,N)
      WRITE (*,600) (YEM(I),I=1,N)



7     IF (DELTA1.LT.0. .AND. KNUM.GE.MAXSOL) GO TO 8
      IF (KNUM.LT.MAXSOL) GO TO 32
      IF (NFLAG.EQ.1) GO TO 8
      DELTA1=-DELTA1
      GO TO 31


8     WRITE (UNITLST,*) '  BIFPACK:  leaving routine  SWITA. '
      RETURN
600   FORMAT (' ',4e19.11)
601   FORMAT (/' Branch point by interpolation:')
602   FORMAT (' Solution of branching system after',I4,'  iterations:')
603   FORMAT (/' emanating solution after',I3,'  iterations:')
604   FORMAT (' METHOD = ',I3,' , fixed component:',I2)
605   FORMAT (' METHOD = ',I3,' , fixing asymmetry of components ',2I3)
606   FORMAT ('    PARAMETER:',E14.7,
     1       /'    obtained vector:')
608   FORMAT ('    vector of approximation to linearization:')
613   FORMAT (' solution at branching parameter:',e19.11,
     1        ' ,  second par.:',e19.11)
680   FORMAT (E20.10)
9300  FORMAT (a1,i2,i4,2i2,a1,a2,18e20.13)
      END



      SUBROUTINE FCNLAR (TDUMMY,Y,F)        
C--  TO BE USED IN CONNECTION WITH  S W I T A
      DOUBLE PRECISION ETACO,ZETACO,PARCO,Y,F,YZ,FZ,Z,TDUMMY 
      PARAMETER (NDIM=12,NDIML=14)
      DIMENSION Y(NDIML),F(NDIML),YZ(NDIML),FZ(NDIML),Z(NDIML,NDIML)
C  DIMENSIONS AT LEAST:  Y(N+N+1),F(N+N+1),Z(N,N),YZ(N+1),FZ(N+1)
      COMMON /CONTRO/ ETACO,ZETACO,PARCO, NCO,N1CO,N2CO,
     1                KCO,IND2CO,NSECFLAG,NFLAG
      DO 1 I=1,N1CO
1     YZ(I)=Y(I)
      CALL FCN (TDUMMY,YZ,FZ)
      DO 11 I=1,N1CO
11    F(I)=FZ(I)
      IF (NFLAG.EQ.1) GO TO 4

C--  BRANCHING SYSTEM
      CALL FCNLIN (TDUMMY,YZ,Z)
      DO 2 I=1,NCO
      I1=NCO+I
      F(I1)=0.d0
      DO 2 J=1,NCO
2     F(I1)=F(I1)+Z(I,J)*Y(N1CO+J)
      F(2*NCO+1)=Y(N1CO+KCO)-1.d0
      RETURN

C-- PARALLEL COMPUTATION
4     DO 41 I=1,NCO
41    YZ(I)=Y(N1CO+I)
      YZ(N1CO)=Y(N1CO)
C     YZ CONTAINS Z-SOLUTION
      CALL FCN (TDUMMY,YZ,FZ)
      DO 48 I=1,NCO
      I1=NCO+I
48    F(I1)=FZ(I)
      F(2*NCO+1)=Y(N1CO+KCO)-Y(KCO)-ETACO
      RETURN
      END

      subroutine LOCATE (INDL,INDK,INDEX,DELTA,IHMAX,
     1                   N,UNITLST,UNITOUT,UNIT3,JACOBI)
c subroutine set up in June 1999
      PARAMETER (NDIM=12,NDIML=14)
      DOUBLE PRECISION Y(NDIML),YSP(NDIML)
     1   ,VFU,VFUD,VFUSP,VFUALT,EPS,ETACO,ZETACO
     2   ,PARCO,DELTA
      EXTERNAL FCN,FCNLIN
      CHARACTER*35 NAME1
      character*1 model
      INTEGER UNITLST,UNITOUT,UNIT3
      COMMON /CONTRO/ ETACO,ZETACO,PARCO, NCO,N1CO,N2CO,
     1                KCO,IND2CO,NSECFLAG,NFLAG

      write(UNITLST,*)'BIFPACK.a1- branch switching: Enter LOCATE.'
      write (UNITLST,*) 'INDL, INDK, k, steplength, number of steps:'
      write (UNITLST,*) INDL,INDK,INDEX,DELTA,IHMAX
      write (UNITLST,*) ' '
      NCO=N
      N1=N+1
      rewind UNIT3
      do 1 i=1,n1
      read (UNIT3,*) y(i)
1     continue
      NOS=0
      ITMAX=20
C    THE NEXT VARIABLES ARE FOR STANDARD CONTIN.:
      N1CO=N1
      ZETACO=0.d0
      IND2CO=1

C  A SIMPLE CONTINUATION:
      ETA=Y(INDEX)
      DO 7 ICONT=1,IHMAX
      KCO=INDEX
      ETA=ETA+DELTA
      ETACO=ETA
      ZETACO=0.d0
      IND2CO=1
      IF (ICONT.EQ.1) GO TO 31
      DO 3 I=1,N1+N
3     YSP(I)=Y(I)
31    Y(KCO)=ETACO
      J=-1
      CALL NOLE (N1,Y,FCN,EPS,ITMAX,J)
      IF (J.LT.0) STOP
      CALL PRINT (Y,ICONT,UNITLST)
      J=1
      CALL TESTA (N,Y,J,INDL,INDK,VFU,N1,1,VFUD,UNITLST)
      IF (ICONT.LT.3) GO TO 69
      VORZ=VFUSP*VFU
      IF (VORZ.GT.0.) GO TO 69
Check for poles:
      IF (VFUALT.GT.ABS(VFUSP)) GO TO 5
      WRITE (UNITLST,602)
      GO TO 69

5     CONTINUE
C    Now, branch switching is to be carried out.
C    The efficient choice is  MODEB=0, METHOD=3
C    S W I T A   writes data of emanating solutions to file 'START' (UNIT 3)
C    for a comfortable start to trace the new branch. 
      MODEB=0
      METHOD=3
      WRITE (*,*) ' Zero of test function found.'
      if (jacobi.eq.1) then 
C        CALL CHECKD (T1,Y,N,N,FCN,FCNLIN,IFLAG,UNITLST)
        WRITE (*,*) ' Enter 1 for solving branching system, else 0'
        read (*,*) MODEB 
        if (MODEB.gt.0.and.NDIML.lt.N+N+2) then 
             MODEB=0
             write (UNITLST,*) 'Dimension too small.'
             write (*,*) 'Dimension too small.'
             endif
        if (MODEB.gt.0) MODEB=INDk
        endif
C      WRITE (*,*) ' enter METHOD :       ' 
C      READ (*,*) METHOD
C    (For a closer calculation of the bifurcation point choose
C    MODEB=INDK ;  needs double dimension.)
      CALL SWITA (N,Y,VFU,YSP,VFUSP,MODEB,METHOD,UNITLST,UNITOUT
     1            ,UNIT3)
69    WRITE (UNITOUT,9300) 'R',-1,ICONT,n,1,MODEL,'  ',(Y(I),I=1,N+1)
      IF (ICONT.GT.1) VFUALT=ABS(VFUSP)
      VFUSP=VFU
7     CONTINUE
      return
601   FORMAT (' ',6E13.6)
602   FORMAT (/' POLE')
603   FORMAT ('0SOLUTION AFTER',I3,'  ITER.,  PARAMETER:',E14.7)
9300  FORMAT (a1,i2,i4,2i2,a1,a2,13e20.13)
      end
