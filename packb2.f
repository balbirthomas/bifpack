C********************************************** 
C**  THIS IS  packB2  version  23 Jan 1996   **
C**  ---------------                         **
C**    INCLUDING ROUTINES FOR HANDLING THE   **
C**    BRANCHING OF BOUNDARY-VALUE PROBLEMS  **
C**    OF SYSTEMS OF O.D.E.s                 **
C**********************************************


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
C
C     REMARK CONCERNING STORAGE 
C     ==========================
C        IN CASE ONE DISPENSES WITH THE OPTIONS THAT DOUBLE THE NUMBER
C        OF EQUATIONS TO BE SOLVED,
C        THE LENGTH AND STORAGE REQUIREMENTS CAN BE
C        REDUCED SIGNIFICANTLY. 


      SUBROUTINE TESTB  (N,M,Y,J,INDL,INDK,BFU,N1,UNITLST)
      DOUBLE PRECISION Y,E,HA,EM,RS,G,A,D,BH, S,EA1,BFU,BFUD
     1    ,EPMACH
      INTEGER UNITLST
      PARAMETER (NDIM=12,NDIML=14)
      DIMENSION  E(NDIML,NDIML),HA(NDIML),EM(NDIML),RS(NDIML),
     1           G(NDIML,NDIML,6),A(NDIML,NDIML), D(NDIML),
     2           NPIV(NDIML),BH(NDIML),Y(NDIML,7)
C   DIMENSIONS AT LEAST:  E(N,N),G(N,N,M-1),A(N,N), Y(.,M) SEE BELOW ,
C                         H(N),EM(N),RS(N),D(N),NPIV(N),BH(N)

      COMMON /BIF/ G,A

C--  SUBROUTINE  T E S T B  (FORMERLY  B R O D E)
C    FOR CALCULATION OF A BRANCHING TESTFUNCTION
C    AND AN INITIAL GUESS FOR THE BRANCHING SYSTEM
C    NEEDED SOFTWARE:  ODEBC1 , WITH ABOVE COMMON-STATEMENT
C---------------------------------------------------------------
C    INPUT:  N       : NUMBER OF COMPONENTS OF RIGHT-HAND SIDE OF THE
C                      BOUNDARY-VALUE PROBLEM (WITHOUT AUXILIARY
C                      VARIABLES)
C            M       : NUMBER OF NODES OF MULTIPLE SHOOTING
C            Y       : ARRAY CONTAINING AN APPROXIMATION TO A SOLUTION.
C                      Y(I,J) IS VALUE OF I-TH COMPONENT  AT NODE T=T(J).
C            N1      : .GT.N : NUMBER OF EQUATIONS INCLUDING AUXILIARY
C                              VARIABLES
C                      = 0   : NO CALCULATION OF INITIAL GUESS
C            INDL, INDK : INDICES L AND K , SEE NUMER.MATH. 32 (1979), 51-
C                      THE K-TH COMPONENT OF THE INITIAL VALUE OF  Y
C                      MUST NOT BE FIXED BY THE BOUNDARY CONDITIONS
C            J         .GT.0 : VALUE OF THE TESTFUNCTION IS PRINTED
C                      = 0   : NO PRINTING
C    OUTPUT: BFU     : VALUE OF BRANCHING TESTFUNCTION AT Y
C            Y       : CONTAINS INITIAL GUESS FOR BRANCHING SYSTEM IN
C                      CASE N1.GT.N ; DIMENSION AT LEAST  Y(N1+N,M)
C                       THE FIRST N1 COMPONENTS REMAIN UNCHANGED.
C            J         .LT.0 : MATRIX SINGULAR, NO TESTFUNCTION

C    REMARKS:
C    -------  THIS VERSION OF TESTB USES THE DATA OBTAINED BY ODEBC1
C             (ARRAYS G, A). THESE DATA ARE USUALLY UPDATED BY RANK-1
C             APPROXIMATIONS. IN SEVERAL APPLICATIONS (EXTRAPOLATION TO
C             A ZERO OF TESTFUNCTION) MORE ACCURATE VALUES ARE NEEDED.
C             IN SUCH A CASE, RE-CALL ODEBC1 FOR
C             ONE NUMERICAL DIFFERENTIATION IN THE SOLUTION POINT.
C    FIRST VERSION 1980 (PUBLISHED IN ISNM 54 (1980) P.163-175)
C----------------------------------------------------------------------------

      EPMACH=1.d-13
      BFU=0.d0
      IPR=J
      M1=M-1
      IF (INDK.GT.0.OR.INDL.LE.0) GO TO 13
      WRITE (11,*) ' INDK is not determined. TERMINATION in TESTB.'
      RETURN

C--  CALCULATION OF MATRIX  E
13    DO 1 I=1,N
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
      IF (INDK.LE.0.OR.INDL.LE.0) GO TO 28

C--  PREPARING MATRIX  E(L,K)  AND RIGHT-HAND SIDE
      DO 25 J=1,N
      EA1=E(INDL,J)
      E(INDL,J)=0.d0
      EM(J)=EA1
25    RS(J)=0.d0
      RS(INDL)=1.d0
      E(INDL,INDK)=1.d0
28    IRANK=N

      CALL BGDECO (N,N,E,NDIML,D,NPIV,EPMACH,IRANK)
      IF (INDK.GT.0.AND.INDL.GT.0)
     1CALL BGSOLV (N,N,E,NDIML,D,NPIV,RS,HA,BH)
      IF (IRANK.EQ.N) GO TO 3
      J=-1
      IF (IPR.GT.0) WRITE (UNITLST,*) ' Matrix singular (in TESTB).'
      RETURN

C--  CALCULATION OF DETERMINANT
3     BFUD=D(1)
      DO 33 J=2,N
33    BFUD=BFUD*D(J)
      IF (INDK.GT.0.AND.INDL.GT.0)  GO TO 42
      BFU=BFUD
      IF (IPR.GT.0) WRITE (UNITLST,601) INDL,INDK,BFU,BFUD
      GO TO 9


C--  CALCULATION OF TESTFUNCTION
42    DO 4 I=1,N
4     BFU=BFU+EM(I)*HA(I)
      IF (IPR.GT.0) WRITE (UNITLST,601) INDL,INDK,BFU,BFUD
      IF (N1.LT.N) GO TO 9


C--  CALCULATION OF INITIAL GUESS
      DO 52 I=1,N
      I1=N1+I
      Y(I1,M)=0.d0
52    Y(I1,1)=HA(I)
      IF (M.LE.2) GO TO 9
      M2=M-2
      DO 56 J=1,M2
      J1=J+1
      DO 56 I=1,N
      I1=N1+I
      S=0.
      DO 54 K=1,N
      K1=N1+K
54    S=S+Y(K1,J)*G(I,K,J)
      Y(I1,J1)=S
56    CONTINUE
9     J=1
      RETURN

601   FORMAT ('  test function (',I2,'/',I2,') =', E13.5,
     1                      '  DET =',E10.3)
      END


      SUBROUTINE SWITB (NOB,LEVEL,N,Y1,TAU1,Y2,TAU2,INDK,N1,M,T,
     1       EPS,DIFME,MODEBI,METHOD,MAXSOL,DELTA,IPRINT
     2       ,UNITLST,UNITOUT)

C---  ALGORITHM FOR BRANCH SWITCHING
C     IN BOUNDARY-VALUE PROBLEMS OF O.D.E.,
C     WITH THE OPTION OF CALCULATING TURNING AND BIFURCATION POINTS
C
C     CONSISTING OF THE THREE SUBROUTINES   SWITB, FCNLAR, BCLAR
C     THE ROUTINES FCNLAR AND BCLAR ARE INTERNAL PROCEDURES,
C     THE USER ONLY CALLS SWITB.
C
C     REFERENCE: R.SEYDEL, NUMER.MATH. 41 (1983), PP.93-116 ,
C                SIMPLIFIED VERSION USING INITIAL GUESS EQ.(10A)
C-------------------------------------------------------------------

C                 THIS VERSION IS TO BE USED IN CONNECTION WITH A
C                MULTIPLE SHOOTING CODE SUCH AS  ODEBC1  WITH
C                APPROPRIATE INTEGRATOR  D I F M E
C                AND THE BRANCHING-SOFTWARE SUBROUTINE   T E S T B .

C--  LEVEL = 1  :  STANDARD VERSION  FOR BRANCH SWITCHING:
C    ------------------------------
C               - THE BRANCH THAT INCLUDES Y1 AND Y2 CAN BE
C                 PARAMETRIZED BY THE PARAMETER.
C               - THE BIFURCATION POINT IS NOT PRECISELY APPROXIMATED.
C               - ONE EMANATING SOLUTION IS TO BE CALCULATED ON EITHER
C                 SIDE OF THE BIFURCATION POINT.
C               - SWITB  CHOOSES AN APPROPRIATE METHOD.
C               - NO FURTHER SUBROUTINES ARE USED EXCEPT THE
C                 PROBLEM-DEFINING SUBROUTINES  FCN  (FOR RIGHT-HAND SIDE)
C                 AND  BC  (FOR BOUNDARY CONDITIONS).
C                 FCN(T,Y,F)  FOR F=F(T,Y,PAR) ,
C                 THE PARAMETER IS GIVEN BY Y(N+1).
C                 BC(YA,YB,R)  FOR R=R(YA,YB) .
C               - THE CALLING-SEQUENCE PARAMETERS MODEBI, METHOD, DELTA, MAXSOL,C                     IPRINT      NEED NOT BE DEFINED.

C--  ESSENTIAL PARAMETERS OF CALLING SEQUENCE:
C    ----------------------------------------
C        N    : NUMBER OF DIFFERENTIAL EQS WITHOUT AUXILIARY VARIABLES
C        N1   : NUMBER OF DIFFERENTIAL EQS WITH AUXILIARY VARIABLES.
C               TWO CASES ARE ALLOWED:
C               N1=N+1 IN CASE OF FIXED BOUNDARY CONDITIONS.
C               N1=N+2 IN THE SPECIAL CASE OF PERIODIC BOUNDARY
C                      CONDITIONS IF METHOD=6  MIGHT BE USED.
C               IN SUBROUTINE  FCN  THE PARAMETER IS ALWAYS GIVEN BY Y(N+1)
C        Y1   : ARRAY OF ONE SOLUTION, CONTAINS OUTPUT OF  TESTB
C               Y1(N+1,.) CONTAINS PARAMETER
C               Y1(N1+I,.) CONTAINS H-VECTORFUNCTION (I=1,..,N)
C        TAU1 : TESTFUNCTION CORRESPONDING TO Y1 (OUTPUT OF TESTB)
C        Y2, TAU2  : SOLUTION AND TESTFUNCTION AT ANOTHER VALUE OF PARAMETER
C               (Y1 AND Y2 ARE TWO CONSECUTIVE SOLUTIONS OF A
C               CONTINUATION PROCEDURE NOT FAR FROM THE BIFURCATION )
C        INDK : INDEX  K (=KCO) USED FOR TESTB
C        M    : NUMBER OF NODES  AND
C        T    : VECTOR OF NODES FOR MULTIPLE SHOOTING
C        EPS  : RELATIVE ERROR TOLERANCE FOR ODEBC1
C        NOB  : number of branch (delivered by the user)
C-----------------------------------------------------------------------------

C--  LEVEL = 2  :  SPECIAL OPTIONS  (CALLING SEQUENCE HAS TO BE FULLY
C    ----------------------------   DEFINED.  Y1 AND Y2 NEED NOT BE SITUATED
C                                   ON A BRANCH THAT CAN BE PARAMETRIZED
C                                   BY THE PARAMETER.)
C        MODEBI = 1 : BIFURCATION OR TURNING POINT IS CALCULATED BY
C                      SOLVING THE BRANCHING SYSTEM.
C                      TWO ADDITIONAL SUBROUTINES ARE USED:
C                      FCNLIN     FOR PARTIAL DERIVATIVES OF RIGHT-HAND SIDE
C                               FCNLIN(T,Y,A)  ,
C                                  A(I,J)=DF(I)/DY(J) ,  1.LE. I,J .LE.N  ,
C                                  Y(N+1)=PAR.
C                      BCLIN(YA,YB,R)      FOR EVALUATING THE LINEARIZED
C                                          BOUNDARY CONDITIONS (IN R).
C                = 0 : APPROXIMATION OF BRANCH POINT BY INTERPOLATION OR
C                      EXTRAPOLATION.  THE BRANCHING SYSTEM IS NOT SOLVED.
C                      FCNLIN  AND  BCLIN  ARE NOT NEEDED.
C                      THREE DIFFERENT APPROXIMATIONS DEPENDING ON METHOD.
C
C        METHOD  : FOR CHOOSING THE METHOD FOR BRANCH SWITCHING
C                .LE.0 : THE BRANCH CAN NOT BE PARAMETRIZED BY THE PARAMETER:
C                        THE Y-SOLUTIONS ARE CONTAINED IN A BRANCH
C                        PERPENDICULAR TO THE PARAMETER-AXIS.
C                = -1 : FOR CALCULATION OF A TURNING POINT.
C                       NO BRANCH SWITCHING, ONLY  MAXSOL=0  MAKES SENSE.
C                =  0 : CASE OF PITCHFORK BIFURCATION.
C                       IN ORDER TO OBTAIN A UNIQUE APPROXIMATION OF THE
C                       BIFURCATION  POINT, Y1 AND Y2 SHOULD BE SITUATED
C                       ON THE SAME HALF-BRANCH.
C                       IF  MAXSOL.GT.0  THE ALGORITHM TRIES TO CALCULATE TWO
C                       SOLUTIONS ON THE BRANCH THAT CAN BE PARAMETRIZED BY
C                       THE PARAMETER, ONE ON EITHER SIDE OF THE BIF.
C                .GT.0 : THE Y-BRANCH CAN BE PARAMETRIZED BY THE PARAMETER
C                =  1 : PARALLEL COMPUTATION  (SLOW, NOT RECOMMENDED)
C                =  2 : SIMPLE METHOD
C                =  3 : EXPLOITING SYMMETRY BREAKING,
C                       TYPE OF SYMMETRY IS CHECKED BY SWITB.
C                .GT.3 : TYPE OF SYMMETRY IS FIXED BY THE USER.
C                =  4 : SYMMETRY TYPE (6A) OF REF.
C                =  5 : SYMMETRY TYPE (6B) OF REF.
C                =  6 : SYMMETRY TYPE (6C) OF REF., FOR PERIODIC OSCILL.,
C                       THE DERIVATIVE OF THE (N1=N+2)-TH VARIABLE IN
C                       SUBROUTINE  FCN  HAS TO BE SUITABLY DEFINED:
C                       ONE NODE (NUMBER JM) HAS TO BE THE MIDDLE OF THE
C                       INTEGRATION INTERVAL,  T(JM)=(A+B)/2  .   THEN:
C                            DY(N+2)=DY(INDK)  FOR SUBINTERVAL J.LT.JM
C                            DY(N+2)=0         FOR SUBINTERVAL J.GE.JM  .
C        MAXSOL  : NUMBER OF EMANATING SOLUTIONS TO BE CALCULATED
C        DELTA   : DISTANCE DELTA0 IN EQ.(4) OF REF., E.G.  =0.02
C        IPRINT  : =1 FOR PRINTOUT, =0 FOR NO PRINTING
C
C    NOTE: LEVEL=1 USES THE PARAMETERS:
C    ----             MODEBI=0 , METHOD=3, DELTA=0.02, MAXSOL=1, IPRINT=1
C-------------------------------------------------------------------------------


      DOUBLE PRECISION T,Y,Y1,Y2,Y0, TAU1,TAU2,FAK1,FAK2,TJINT,TJINT2
     1  ,EPS,EPMACH,EPSB,ETACO,PARCO,ZETACO,DELTA,DELL,DEL1,
     2  SECPAR1,SECPARINC 
      INTEGER UNITLST,UNITOUT
      PARAMETER (NDIM=12,NDIML=14)
      DIMENSION T(7),Y1(NDIML,7),Y2(NDIML,7),Y0(NDIML,7),Y(NDIML,7)
C   DIMENSIONS AT LEAST:  T(M),Y1(.,M),Y2(.,M),Y0(.,M),Y(.,M),
C                   WHERE  .  STANDS FOR N1+N
      COMMON /CONTRO/ ETACO,ZETACO,PARCO, NCO,N1CO,N2CO,
     1                KCO,IND2CO,NSECFLAG,NFLAG
      EXTERNAL FCNLAR,BCLAR,DIFME

      KCO=INDK
C     TO BE CHANGED IF NECESSARY:
      EPMACH=1.d-12

      NCO=N
      N1CO=N1
C     NO PRINTING BY SOLVER:
      IPRISO=-1
      IF (LEVEL.EQ.1) THEN
          MODEBI=0
          METHOD=3
          DELTA=0.02d0
          MAXSOL=1
          IPRINT=1
          END IF
      IF (IPRINT.GT.0) THEN
            WRITE (UNITLST,603)
            WRITE (UNITLST,*) 'BIFPACK: entering routine SWITB.'
            WRITE (*,*) 'BIFPACK: entering routine SWITB.'
            WRITE (UNITLST,603)
            END IF
      NL=N1+N
      NP1=N+1

C     EQ. (4) OF REF.:
      DELL=ABS(Y1(KCO,1))*DELTA
      IF (DELL.LT.DELTA) DELL=DELTA

C===============================================================
C--  APPROXIMATION OF BRANCH POINT

C--  BRANCH POINT  Y0  BY INTERPOLATION OR EXTRAPOLATION:
      IF (TAU1*TAU2.LE.0.d0.AND. METHOD.GT.0) GO TO 2
      FAK1=ABS(TAU1/TAU2)
      IF (FAK1.GT.1.2 .OR. FAK1.LT.0.8) GO TO 2
C     NO EXTRAPOLATION, USE Y0=Y2
14    IF (IPRINT.GT.0 .AND. METHOD.GT.0) WRITE (UNITLST,606)
      IF (MODEBI.EQ.0) MODEBI=-1
      DO 15 I=1,NL
      DO 15 J=1,M
15    Y0(I,J)=Y2(I,J)
      GO TO 28

2     IF (IPRINT.GT.0) THEN
            WRITE (UNITLST,604)
            WRITE (*,604)
            END IF
      FAK2=1.d0/(TAU1-TAU2)
      FAK1=TAU1*FAK2
      FAK2=TAU2*FAK2
      IF (METHOD.NE.0) GO TO 24
C     APPROXIMATION OF PITCHFORK BIF. ASSUMING PARABOLA
      DO 21 J=1,M
21    Y0(NP1,J)=FAK1*Y2(NP1,J)-FAK2*Y1(NP1,J)
      IF (TAU1*TAU2.LT.0.d0) GO TO 14
      IF (ABS(TAU2).LT.ABS(TAU1)) GO TO 215
      FAK1=SQRT(TAU1/TAU2)
      FAK2=1.d0/(FAK1-1.d0)
      FAK1=FAK1*FAK2
      GO TO 218
215   FAK2=SQRT(TAU2/TAU1)
      FAK1=1.d0/(1.d0-FAK2)
      FAK2=FAK1*FAK2
218   DO 22 J=1,M
      DO 22 I=1,NL
22    IF (I.NE.NP1) Y0(I,J)=FAK1*Y2(I,J)-FAK2*Y1(I,J)
      GO TO 28

24    DO 25 J=1,M
      DO 25 I=1,NL
25    Y0(I,J)=FAK1*Y2(I,J)-FAK2*Y1(I,J)
      IF (METHOD.GT.0) GO TO 28
C     APPROXIMATION OF TURNING POINT ASSUMING PARABOLA
      FAK1= Y2(NP1,1)*TAU1*TAU1
      FAK2=Y1(NP1,1)*TAU2*TAU2
      FAK2=(FAK1-FAK2)/((TAU1+TAU2)*(TAU1-TAU2))
      DO 26 J=1,M
26    Y0(NP1,J)=FAK2
28    IF (IPRINT.GT.0) THEN
           WRITE (UNITLST,601) Y0(NP1,1)
           WRITE (*,601) Y0(NP1,1)
           WRITE (UNITLST,609)
           WRITE (UNITLST,602) (Y0(I,1),I=1,N)
           END IF
      IF (NSECFLAG.EQ.0) WRITE (UNITOUT,9300) 
     1    'r',NOB,0,N,1,'b','B ',(Y0(I,1),I=1,N+1)
      IF (NSECFLAG.EQ.1) WRITE (UNITOUT,9300) 
     1    'r',NOB,0,N,2,'b','B ',(Y0(I,1),I=1,N+1),PARCO
      IF (MODEBI.LE.0) GO TO 29


C--   BIF.POINT BY SOLVING THE BRANCHING SYSTEM (EQ.(9) OF REF.):
      NFLAG=-1
      EPSB=EPS
C     IN THE CASE OF A BIFURCATION POINT, THE REQUIRED ACCURACY SHOULD NOT
C           BE STRINGENT, SINCE THE CONVERGENCE OF A GENERAL COMPONENT MAY
C           BE POORE. DUE TO SUPERCONVERGENCE, THE PARAMETER CAN BE EXPECTED
C           TO BE MORE ACCURATE THAN  EPSB  WILL REQUIRE.
      IF (METHOD.GE.0 .AND. EPS.LT.EPMACH**0.3) EPSB=EPMACH**0.3
      nsecflag1=0
      nsecflag2=0
      secpar1=PARCO
291   J=IPRISO
      nsecflag1=nsecflag1+1
      CALL ODEBC1 (FCNLAR,BCLAR,DIFME,NL,M,T,Y0,EPSB,10,J)
      IF (IPRINT.GT.0) THEN
          write (UNITLST,*) ' '
          IF (J.GT.0) WRITE (UNITLST,605)
          WRITE (UNITLST,601) Y0(NP1,1)
          WRITE (UNITLST,609)
          WRITE (UNITLST,602) (Y0(I,1),I=1,N)
          END IF
      IF (J.GT.0) THEN
         if (nsecflag1.eq.1) then
            WRITE (*,605)
            WRITE (*,601) Y0(NP1,1)
            WRITE (*,609)
            WRITE (*,602) (Y0(I,1),I=1,N)
            if (nsecflag.eq.1.and.nsecflag1.eq.1) then
              write (*,*) 'Control of second parameter.'
              write (*,*) 'current value:', PARCO
              write (*,*) 'Vary second parameter? ',
     1                   '  (1=yes, 0=no ) ENTER:'
              read (*,*) nsecflag2
              do 292 j=1,M
              do 292 i=1,NL
292           y(i,j)=y0(i,j)
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
         IF (NSECFLAG2.EQ.1) WRITE (*,613) Y0(NP1,1),PARCO
         IF (NSECFLAG2.EQ.1) WRITE (UNITLST,613) Y0(NP1,1),PARCO
         if (nsecflag.ne.1) WRITE (UNITOUT,9300)
     1         'R',NOB,0,N,1,'b','B ',(Y0(I,1),I=1,N+1)
         if (nsecflag.eq.1.AND.nsecflag1.eq.1) WRITE (UNITOUT,9300)
     1         'R',NOB,0,N,2,'b','B ',(Y0(I,1),I=1,N+1),PARCO
         if (nsecflag2.eq.1) WRITE (UNITOUT,9300)
     1         'R',0,0,N,2,'b','B ',(Y0(I,1),I=1,N+1),PARCO
         if (nsecflag1.lt.nsecflag3
     1       .and.nsecflag2.eq.1)  then
             PARCO=PARCO+secparinc
             go to 291
             endif
         ENDIF
29    IF (IPRINT.GT.0) WRITE (UNITLST,603)
      IF (NSECFLAG2.EQ.1) then
          if (nsecflag1.lt.nsecflag3) then
            write (*,*) 'The bifurcation curve could not '
     1                  ,'be traced farther.'  
            write (UNITLST,*) 'The bifurcation curve could not'
     1                  ,' be traced farther.'  
            endif 
          PARCO=secpar1
          do 293 j=1,M
          do 293 i=1,NL 
293       y0(i,j)=y(i,j)
          endif
      if (nsecflag1.gt.1) write (*,*) 'second parameter'    
     1    ,' retained to the old value ',PARCO    
      if (nsecflag1.gt.1) write (11,*) 'second parameter'
     1    ,' retained to the old value ',PARCO    
      IF (MAXSOL.LE.0) GO TO 99
      IF (METHOD.LE.0) GO TO 5
      IF (METHOD.EQ.1) GO TO 7

C=================================================
C--   MAIN PART  (METHODS THAT SOLVE THE EQ.S (5) OR (7) OF REF.)

401   DEL1=0.d0
      DO 45 ICON=1,MAXSOL
      DEL1=DEL1+DELL
      DO 4 J=1,M
      DO 41 I=NP1,N1
41    Y(I,J)=Y0(I,J)
      DO 4 I=1,N
      I1=I+N1
4     Y(I,J)=Y0(I,J)+DEL1*Y0(I1,J)
      NFLAG=2
      IF (METHOD.EQ.6 .AND. DELL.LT.0.) GO TO 433
      IF (METHOD.EQ.6 .AND. ICON.GT.1) GO TO 433
      IF (METHOD.GT.3) GO TO 423
      IF (METHOD.EQ.2) GO TO 422

C--  CHECKING SYMMETRIES

      I=0
42    I=I+1
      IF (I.LE.N) GO TO 43
      IF (N1-N.LT.2) GO TO 422
      METHOD=6
      JINT=0
431   JINT=JINT+1
      IF (JINT.EQ.M) GO TO 422
      TJINT=ABS(Y(N+2,JINT+1)-Y(N+2,JINT))
      TJINT2=MAX(1.d0,DBLE(ABS(Y(N+2,JINT))))*EPS*10.d0
      IF (TJINT.GT.TJINT2) GO TO 431
433   DO 435 J=1,M
      Y(N+2,J)=Y(KCO,J)
435   IF (J.GT.JINT) Y(N+2,J)=Y(N+2,JINT)
      ETACO=Y(N+2,1)+Y(N+2,M)
      GO TO 423

43    EPSA=10.d0*EPS*MAX(1.d0,DBLE(ABS(Y0(I,1))))
      IF (ABS(Y0(I,1)-Y0(I,M)).LT.EPSA) GO TO 42
      IF (ABS(Y0(KCO,1)).LT.EPS*10.d0) GO TO 422
      EPSA=EPS*10.d0*MAX(1.d0,DBLE(ABS(Y0(KCO,1))))
      IF (ABS(Y0(KCO,1)-Y0(KCO,M)).LT.EPSA) METHOD=4
      IF (ABS(Y0(KCO,1)+Y0(KCO,M)).LT.EPSA) METHOD=5
      IF (METHOD.GT.3) GO TO 423
422   ETACO=Y(KCO,1)
      GO TO 429
423   IF (METHOD.EQ.4) ETACO=Y(KCO,1)-Y(KCO,M)
      IF (METHOD.EQ.5) ETACO=Y(KCO,1)+Y(KCO,M)
      IF (ABS(ETACO).LT.SQRT(EPS)) GO TO 422
      NFLAG=METHOD
429   IF (METHOD.EQ.NFLAG .AND. IPRINT.GT.0) THEN
            WRITE (UNITLST,608) KCO,NFLAG
            WRITE (*,608) KCO,NFLAG
            END IF
      IF (METHOD.NE.NFLAG .AND. IPRINT.GT.0) THEN
            WRITE (UNITLST,608) KCO,NFLAG,METHOD
            WRITE (*,608) KCO,NFLAG,METHOD
            END IF


      J=IPRISO
      CALL ODEBC1 (FCNLAR,BCLAR,DIFME ,N1,M,T,Y,EPS,20,J)
      IF (J.GT.0.AND.NSECFLAG.EQ.0) WRITE (UNITOUT,9300)
     1            'R',0,1,N,1,'b','E ',(Y(I,1),I=1,N+1)
      IF (J.GT.0.AND.NSECFLAG.EQ.1) WRITE (UNITOUT,9300)
     1     'R',0,1,N,2,'b','E ',(Y(I,1),I=1,N+1),PARCO
      IF (IPRINT.GT.0.AND.J.GT.0) THEN
              WRITE (UNITLST,601) Y(NP1,1)
              WRITE (*,601) Y(NP1,1)
              WRITE (UNITLST,607) J
              WRITE (UNITLST,610)
              DO 48 J=1,M
   48         WRITE (UNITLST,602) (Y(I,J),I=1,N)
              END IF
45    CONTINUE
      IF (IPRINT.GT.0) WRITE (UNITLST,603)
      IF (DELL.LT.0.d0) GO TO 49
      DELL=-DELL
      GO TO 401
49    CONTINUE
      GO TO 99

C==================================================
C--  CALCULATE TWO SOLUTIONS ON THE BRANCH THAT CAN BE
C    PARAMETRIZED BY THE PARAMETER.

5     NFLAG=0
C     THE FOLLOWING DISTANCE  DELL  MIGHT BE TOO LARGE IF A BRANCH CAN BE
C     PARAMETRIZED ONLY IN A VERY SMALL NEIGHBORHOOD:
      DELL=ABS(Y0(NP1,1))*DELTA
      IF (DELL.LT.10.d0*EPS) DELL=10.d0*EPS
      IND=1
      GO TO 531
53    IND=2
      DELL=-DELL
531   DO 54 J=1,M
      DO 541 I=1,N1
541   Y(I,J)=Y0(I,J)
54    Y(NP1,J)=Y0(NP1,J)+DELL
      ETACO=Y(NP1,1)
      J=IPRISO
      CALL ODEBC1 (FCNLAR,BCLAR,DIFME,N1,M,T,Y,EPS,20,J)
      IF (J.LT.0) GO TO 591
      IF (IPRINT.GT.0) WRITE (UNITLST,601) Y(NP1,1)
      IF (IPRINT.GT.0) WRITE (UNITLST,607) J
      DO 59 J=1,M
59    IF (IPRINT.GT.0) WRITE (UNITLST,602) (Y(I,J),I=1,N)
591   CONTINUE
      IF (IPRINT.GT.0) WRITE (UNITLST,603)
      IF (IND.LT.2) GO TO 53
      GO TO 99

C==================================================
C--   PARALLEL COMPUTATION (SOLVING EQ.(3) OF REF.):

7     NFLAG=-2
      IF (IPRINT.GT.0) WRITE (UNITLST,608) KCO,METHOD
71    DEL1=0.d0
      DO 75 ICON=1,MAXSOL
      DEL1=DEL1+DELL
      DO 72 J=1,M
      DO 721 I=NP1,N1
721   Y(I,J)=Y0(I,J)
      DO 72 I=1,N
      Y(I,J)=Y0(I,J)
      I1=I+N1
72    Y(I1,J)=Y0(I1,J)*DEL1
      ETACO=Y(KCO+N1,1)
      J=IPRISO
      CALL ODEBC1 (FCNLAR,BCLAR,DIFME ,NL,M,T,Y,EPS,20,J)
      IF (J.LT.0) GO TO 75
      IF(IPRINT.GT.0) WRITE (UNITLST,601) Y(NP1,1)
      IF (IPRINT.GT.0) WRITE (UNITLST,607) J
      DO 78 J=1,M
      DO 77 I=1,N
77    Y(I,J)=Y(I,J)+Y(N1+I,J)
78    IF (IPRINT.GT.0) WRITE (UNITLST,602) (Y(I,J),I=1,N)
75    CONTINUE
      IF (IPRINT.GT.0) WRITE (UNITLST,603)
      IF (DELL.LT.0.d0) GO TO 99
      DELL=-DELL
      GO TO 71

99    IF (IPRINT.GT.0) THEN
           WRITE (UNITLST,*) 'BIFPACK: leaving routine SWITB.'
           WRITE (*,*) 'BIFPACK: leaving routine SWITB.'
           WRITE (UNITLST,603)
           END IF
      RETURN
601   FORMAT ('       PARAMETER:',e20.11,'      ')
602   FORMAT ('   ',4E19.11)
603   FORMAT (' ************************************')
604   FORMAT (' Branch point by interpolation-extrapolation:         ')
605   FORMAT (' Branch point obtained by solving the branching system')
606   FORMAT (' bad values of testfunction, Y2-solution used:')
607   FORMAT (' ',I3,' iterations    ')
608   FORMAT (' emanating solution, fixing component ',I2,
     1                  ' , METHOD =',2I2,'    ')
609   FORMAT (' vector of initial values:               ')
610   FORMAT (' values of the components at the nodes:  ')
613   FORMAT (' solution at branching PARAMETER:',e15.9,
     1        ' , second par.:',e15.9)  
9300  FORMAT (A1,I2,I4,2I2,A1,A2,18E20.13)
      END

      SUBROUTINE FCNLAR (T,YL,DYL)   
      DOUBLE PRECISION T,YL,DYL,Y,DY,A,ETACO,PARCO,ZETACO
      PARAMETER (NDIM=12,NDIML=14)
      DIMENSION Y(NDIML),DY(NDIML),YL(NDIML),DYL(NDIML),A(NDIML,NDIML)
C    DIMENSIONS AT LEAST:  Y(N1+N),DY(N1+N),YL(N1+N),DYL(N1+N),A(N,N)
      COMMON /CONTRO/ ETACO,ZETACO,PARCO, NCO,N1CO,N2CO,
     1                KCO,IND2CO,NSECFALG,NFLAG
      CALL FCN (T,YL,DYL)
      DYL(NCO+1)=0.d0
      IF (NFLAG.GE.0) RETURN
      IF (NFLAG.EQ.-2) GO TO 4

      CALL FCNLIN (T,YL,A)
      DO 2 I=1,NCO
      I1=N1CO+I
      DYL(I1)=0.d0
      DO 2 J=1,NCO
      J1=N1CO+J
2     DYL(I1)=DYL(I1)+A(I,J)*YL(J1)
      RETURN

4     DO 41 I=1,N1CO
41    Y(I)=YL(I)
      DO 43 I= 1,NCO
      I1=N1CO+I
43    Y(I)=Y(I)+YL(I1)
      CALL FCN (T,Y,DY)
      DO 6 I=1,NCO
      I1=N1CO+I
6     DYL(I1)=DY(I)-DYL(I)
      RETURN
      END


      SUBROUTINE  BCLAR (YAL,YBL,WL)
      DOUBLE PRECISION YAL,YBL,WL,YA,YB,W,ETACO,PARCO,ZETACO
      PARAMETER (NDIM=12,NDIML=14)
      DIMENSION YAL(NDIML),YBL(NDIML),WL(NDIML),
     1    YA(NDIML),YB(NDIML),W(NDIML)
      COMMON /CONTRO/ ETACO,ZETACO,PARCO, NCO,N1CO,N2CO,
     1                KCO,IND2CO,NSECFLAG,NFLAG
      CALL BC (YAL,YBL,WL)
      IF (NFLAG.GE.0) GO TO 4

      DO 2 I=1,NCO
      I1=I+N1CO
      YA(I)=YAL(I1)
      YB(I)=YBL(I1)
      IF (NFLAG.EQ.-2) YA(I)=YA(I)+YAL(I)
2     IF (NFLAG.EQ.-2) YB(I)=YB(I)+YBL(I)
c     IF (NFLAG.GE.0) GO TO 4
      IF (NFLAG.EQ.-1) CALL BCLIN (YA,YB,W)
      IF (NFLAG.NE.-2) GO TO 26
      NP1=NCO+1
      DO 25 I=NP1,N1CO
      YA(I)=0.d0
25    YB(I)=0.d0
      CALL BC (YA,YB,W)
26    DO 3 I=1,NCO
      I1=I+N1CO
3     WL(I1)=W(I)

4     IF (NFLAG.EQ. 0) WL(N1CO)=YAL(NCO+1)-ETACO
      IF (NFLAG.EQ. 2) WL(N1CO)=YAL(KCO)-ETACO
      IF (NFLAG.EQ. 4) WL(N1CO)=YAL(KCO)-YBL(KCO)-ETACO
      IF (NFLAG.EQ. 5) WL(N1CO)=YAL(KCO)+YBL(KCO)-ETACO
      IF (NFLAG.EQ. 6) WL(N1CO)=YAL(NCO+2)+YBL(NCO+2)-ETACO
      IF (NFLAG.EQ.-1) WL(N1CO)=YAL(N1CO+KCO)-1.d0
      IF (NFLAG.EQ.-2) WL(N1CO)=YAL(N1CO+KCO)-ETACO
      IF (N1CO-NCO.EQ.2) WL(NCO+1)=YAL(NCO+1)-YAL(KCO)
      RETURN
      END
