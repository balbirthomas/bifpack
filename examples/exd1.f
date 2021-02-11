      SUBROUTINE SETUP (MODEL,N,M,T,NAME1,JACOBI,EPS,
     1                  IFLAGSOL,WINDOW,YSTART)
      double precision eps,ystart,ETACO,ZETACO,PARCO,T,WINDOW
      parameter (NDIML=14)
      dimension ystart(ndiml),T(7),WINDOW(2)
      CHARACTER*35 name1
      character*1 model 
      COMMON /CONTRO/ ETACO,ZETACO,PARCO, NCO,N1CO,N2CO,
     1                KCO,IND2CO,NSECFLAG,NFLAG

C   Enter a short name of your problem:
      name1='exd1 Brusselator, 2 cells'
      model='d'

C   If a particular solution is known, enter the solution into the
C   the vector YSTART, and set IFLAGSOL=1.
      iflagsol=1
      YSTART(1)=2.972
      YSTART(2)=2.349
      YSTART(3)=1.028
      YSTART(4)=3.078
      YSTART(5)=6.
      NSECFLAG=0
      window(1)=-99999.
      window(2)=99999.
C    FIRST, SOME TECHNICAL DETAILS OF THE PARTICULAR
C    PROBLEM ARE INDICATED:

C  N IS THE NUMBER OF DIFFERENTIAL EQUATIONS:
      N=4
C   JACOBI=1  IF THE LINEARIZATION IS PROVIDED IN SUBROUTINE
C   FCNLIN,   ELSE TAKE  JACOBI=0
      JACOBI=1
C  THE RELATIVE ERROR TOLERANCE FOR ALL CALCULATIONS:
      EPS=1.d-4
C**************************************************************
C
C    THE MAINPROGRAM TOGETHER WITH THE SUBROUTINE  B I F L A B
C    COMBINES THE USE OF VARIOUS BIFPACK-TOOLS IN A CONVENIENT,
C    SEMI-AUTOMATED FASHION. THE PROGRAM IS DESIGNED TO BE USED
C    AT A TERMINAL IN INTERACTIVE CALCULATIONS. THE FINAL RESULTS
C    OF EACH RUN ARE STORED IN LOCAL FILES AND USED FOR ANY RESTART.
C    ESPECIALLY, THE FILE 'START' PROVIDES THE  N+2 OR N+1  STARTING
C    VALUES. ON FIRST RUN, THE USER will be prompted for these
C    NUMBERS:
C      1 THROUGH N   : INITIAL VALUES (OR STEADY-STATE VECTOR)
C          N+1       : PARAMETER
C          N+2       : PERIOD    (NOT REQUIRED FOR MODUS =7)
C
C**************************************************************
C**************************************************************

C*
C*     IN THE FOLLOWING, SOME NUMBERS ARE GIVEN FOR TESTING THE
C*     BRUSSELATOR MODEL, TWO CELLS (FOUR EQUATIONS: N=4).
C*
C*     FIRST EXAMPLE: (CHOOSE MODUS= 1) CONTINUATION ALONG PERIODIC
C*                    BRANCH WITH STABILITY ANALYSIS
C*
C*     SECOND EXAMPLE : (CHOOSE  2) CASE OF HOPF BIFURCATION:
C*                    CALCULATE SMALL-AMPLITUDE PERIODIC OSCILLATIONS,
C*                    STARTING FROM A STATIONARY SOLUTION.
C*
C*     THIRD EXAMPLE: (CHOOSE  3) CALCULATE DATA FOR A PHASE
C*                    DIAGRAM
C*
C*     FOURTH EXAMPLE: (CHOOSE  5) SIMPLIFIED BRANCH SWITCHING
C*
C*     FIFTH EXAMPLE: (CHOOSE  7) STEADY-STATE STABILITY ANALYSIS.
C*                    THE PROGRAM WILL PROMPT YOU FOR SOME DATA.
C*                    ENTER, FOR INSTANCE
C*                     10 (CONTIN.STEPS), 0 (FOR FIXING THE PARAMETER),
C*                     0.5 (INITIAL STEP), 2 (variable steps), 2. (MAX.STEP)
C*
C*
C*    These data are to be entered in file START
C*    for generating some test numbers.
C******************************************************************
C*
C*     THIS IS THE EXAMPLE OF A CONTINUATION:

C     FIST THE INITIAL VALUES OF THE SOLUTION WHERE THE
C     CONTINUATION STARTS:
c      X(1)=4.6369
c      X(2)=1.31932
c      X(3)=4.6369
c      X(4)=1.31932
c      PAR=5.548861
c      PERIOD=4.136796

C*************************************************************
C*
C*    THIS IS THE EXAMPLE FOR HANDLING A HOPF POINT:
C*

C  THE FOLLOWING ROUGH INITIAL GUESS OF THE HOPF BIFURCATION
C  IS OBTAINED IN A VERY SIMPLE WAY:
C  (TAKEN FOR INSTANCE FROM THE OUTPUT OF THE ROUTINES
C   C O N T A   AND   S T A S T A , SEE  P A C K A 1 )
C     1) CALCULATE A CONTINUATION ALONG THE STATIONARY BRANCH
C     2) CALCULATE THE EIGENVALUES OF THE LINEARIZATION AT THE
C        SOLUTION POINTS
C     3) CHECK FOR A CHANGE OF SIGN OF A REAL PART,
C        THIS GIVES YOU NEAR-BY STATIONARY SOLUTION  X(1),..,X(N)
C        AND   PAR   AND THE VALUE OF THE CORRESPONDING
C        IMAGINARY PART.
C     4) THE APPROXIMATE PERIOD IS:
C        PERIOD = 6.28 / IMAGINARYPART
c      X(1)=2.
c      X(2)=2.51
c      X(3)=2.
c      X(4)=2.51
c      PAR=4.99
c      PERIOD=3.14
C (EXACT VALUES ARE: PAR=5 , X(2)=X(4)=2.5 , PERIOD=PI=3.14159..)

C********************************************************
C*
C*    THIS IS THE EXAMPLE FOR CALCULATING A GENERAL BIFURCATION:
C*
C     ENTER VECTOR OF INITIAL VALUES, PARAMETER AND PERIOD TO

c      X(1)=0.9984815
c      X(2)=5.644571
c      X(3)=0.9984815
c      X(4)=5.644571
c      PAR=5.636
c      PERIOD=4.42324

C************************************************************
C*
C*    STEADY-STATE DATA:
c      X(1)=1.028
c      X(2)=3.078
c      X(3)=2.972
c      X(4)=2.349
c      PAR=6.
c      PERIOD=0.
C*******************************************************************

      return
      end

      SUBROUTINE FCN (T,Y,F)
      DOUBLE PRECISION T,Y,F,PAR,ETACO,ZETACO,PARCO,
     1    A,D,B,D2,B1,Y12,Y22
      COMMON /CONTRO/ ETACO,ZETACO,PARCO, NCO,N1CO,N2CO,
     1                KCO,IND2CO,NSECFLAG,NFLAG
      PARAMETER (NDIM=12,NDIML=14)
      DIMENSION Y(NDIML),F(NDIML)
      F(N1CO)=Y(KCO)-ZETACO*Y(IND2CO)-ETACO
C     THE PARAMETER IS GIVEN BY  Y(N+1)=Y(N1CO)
      PAR=Y(N1CO)
C     THIS SUBROUTINE DEFINES THE COMPONENTS OF THE RIGHT-HAND SIDE
C     F(1),...,F(n)   depending on the input   Y(1),...,Y(n),PAR  .
C
C     example: BRUSSELATOR, 2 cells

      A=2.d0
      D=1.d0
      B=PAR
C     (THAT IS TO SAY, 'B' IS THE PARAMETER IN THIS EXAMPLE)
      D2=4.d0

      B1=B+1.d0
      Y12=Y(1)*Y(1)
      Y22=Y(3)*Y(3)
      F(1)=A-B1*Y(1)+Y12*Y(2)+D*(Y(3)-Y(1))
      F(2)=B*Y(1)-Y12*Y(2)+D2*(Y(4)-Y(2))
      F(3)=A-B1*Y(3)+Y22*Y(4)-D*(Y(3)-Y(1))
      F(4)=B*Y(3)-Y22*Y(4)-D2*(Y(4)-Y(2))
      RETURN
      END


      SUBROUTINE FCNLIN (T,Y,A) 
      DOUBLE PRECISION T,Y,A,PAR,ETACO,ZETACO,PARCO, D,B,D2,B1
      PARAMETER (NDIM=12,NDIML=14)
      DIMENSION Y(NDIML),A(NDIML,NDIML)
      COMMON /CONTRO/ ETACO,ZETACO,PARCO, NCO,N1CO,N2CO,
     1                KCO,IND2CO,NSECFLAG,NFLAG
C  THIS ROUTINE ESTABLISHES PARTIAL DERIVATIVES OF RIGHT-HAND SIDE.
C     THE PARAMETER IS GIVEN BY  PAR=Y(N+1) (=Y(N1CO)):
      PAR=Y(N1CO)
C  OUTPUT IS  a(i,j)= partial derivative of f(i) with respect to y(j),
C                           i,j=1,...,n  .                       

      D=1.d0
      B=Y(5)
      D2=4.d0

      B1=1.d0+B

      A(1,1)=-B1+2.d0*Y(1)*Y(2)-D
      A(1,2)=Y(1)*Y(1)
      A(1,3)=D
      A(1,4)=0.d0

      A(2,1)=B-2.d0*Y(1)*Y(2)
      A(2,2)=-Y(1)*Y(1)-D2
      A(2,3)=0.d0
      A(2,4)=D2

      A(3,1)=D
      A(3,2)=0.d0
      A(3,3)=-B1+2.d0*Y(3)*Y(4)-D
      A(3,4)=Y(3)*Y(3)

      A(4,1)=0.d0
      A(4,2)=D2
      A(4,3)=B-2.d0*Y(3)*Y(4)
      A(4,4)=-Y(3)*Y(3)-D2
      RETURN
      END

