      SUBROUTINE SETUP (MODEL,N,M,T,NAME1,JACOBI,EPS,
     1                  IFLAGSOL,WINDOW,YSTART)
      double precision T,eps,ystart,ETACO,ZETACO,PARCO,WINDOW
      parameter (NDIML=14)
      dimension ystart(ndiml),WINDOW(2),T(7)
      CHARACTER*35 name1
      character*1 model
      COMMON /CONTRO/ ETACO,ZETACO,PARCO, NCO,N1CO,N2CO,
     1                KCO,IND2CO,NSECFLAG,NFLAG

C   Enter a short name of your problem:
      NAME1='exa1 Brusselator, 3 cells'
      MODEL='a'

C   If a particular solution is known, enter the solution into the
C   the vector YSTART, and set IFLAGSOL=1.
      IFLAGSOL=1
      NSECFLAG=0
C    First, some technical details of the particular
C    problem are indicated:

C  N is the number of state variables (differential equations):
      N=6
C   JACOBI=1  if the linearization is provided in subroutine
C   FCNLIN,   else take  JACOBI=0
      JACOBI=1
C  The relative error tolerance for all calculations chosen to be:
      EPS=1.d-5

C   Initial solution, for instance (may be entered with MODUS=0):
      YSTART(1)=1.446
      YSTART(2)=2.744
      YSTART(3)=3.11
      YSTART(4)=2.54
      YSTART(5)=1.446
      YSTART(6)=2.744
C      PARAMETER=1.435
      YSTART(7)=1.435
      WINDOW(1)=0.d0
      Window(2)=2.d0
C
C  (BIFPACK may prompt you to confirm the data: enter 0)
C
C     THEN WITH MODUS=2, YOU MAY ENTER CONTINUATION OPTIONS.
C     FOR INSTANCE, CHANGE THE DEFAULT VALUES VIA DIALOG AS FOLLOWS:
C        enter option No.:  1  ,  enter  -0.1   (initial step -0.1)
C                           2            0      (parameterized by PAR)
C                           4            80     (this closes the loop)
C     After entering these options, enter  0  to run the continuation.
C
C     THE CONTINUATION IN THIS EXAMPLE IS SENSITIVE TO UNDESIRED BRANCH
C           SWITCHING (THEREFORE WE SUGGEST THE SMALL MAXIMUM STEP LENGTH)

      return
      end


      SUBROUTINE FCN (TDUMMY,Y,F)
      DOUBLE PRECISION Y,F,PAR,A,ETACO,ZETACO,PARCO,TDUMMY
      PARAMETER (NDIM=12,NDIML=14)
      DIMENSION Y(NDIML),F(NDIML)
C   THE FOLLOWING THREE STATEMENTS ARE FOR CONTINUATION AND
C   BRANCH SWITCHING PURPOSES:
      COMMON /CONTRO/ ETACO,ZETACO,PARCO, NCO,N1CO,N2CO,
     1                KCO,IND2CO,NSECFLAG,NFLAG
      PAR=Y(N1CO)
      F(N1CO)=Y(KCO)-ZETACO*Y(IND2CO)-ETACO
C   TDUMMY is dummy.
C   THIS ROUTINE EVALUATES THE PROBLEM-DEFINING FUNCTION.
C   input is Y(1),...,Y(n), and PAR=parameter.
C   NOW, F(1),...,F(n) of THE FUNCTION OF THE RIGHT-HAND SIDE  F(Y,PAR)
C   are entered:  (BRUSSELATOR, 3 CELLS)
      A=Y(1)*Y(1)*Y(2)
      F(1)=2.-7.*Y(1)+A+PAR*(Y(3)-Y(1))
      F(2)=6.*Y(1)-A+10.*PAR*(Y(4)-Y(2))
      A=Y(3)*Y(3)*Y(4)
      F(3)=2.-7.*Y(3)+A+PAR*(Y(1)+Y(5)-2.*Y(3))
      F(4)=6.*Y(3)-A+10.*PAR*(Y(2)+Y(6)-2.*Y(4))
      A=Y(5)*Y(5)*Y(6)
      F(5)=2.-7.*Y(5)+A+PAR*(Y(3)-Y(5))
      F(6)=6.*Y(5)-A+10.*PAR*(Y(4)-Y(6))
      RETURN
      END


      SUBROUTINE FCNLIN (TDUMMY,Y,A)
      DOUBLE PRECISION TDUMMY,Y,A,PAR,ETACO,ZETACO,PARCO,P10,Y12,Y11
      PARAMETER (NDIM=12,NDIML=14)
      DIMENSION Y(NDIML),A(NDIML,NDIML)
      COMMON /CONTRO/ ETACO,ZETACO,PARCO, NCO,N1CO,N2CO,
     1                KCO,IND2CO,NSECFLAG,NFLAG
C  For input the same remarks apply as in SUBROUTINE FCN.
C  Output are the n*n values of the Jacobian,
C   A(i,j)= partial derivative of F(i) with respect to Y(j):
      DO 1 I=1,6
      DO 1 J=1,7
1     A(I,J)=0.
      PAR=Y(7)
      P10=PAR*10.
      Y12=2.*Y(1)*Y(2)
      Y11=Y(1)*Y(1)
      A(1,1)=-7.+Y12-PAR
      A(1,2)=Y11
      A(1,3)=PAR
      A(2,1)=6.-Y12
      A(2,2)=-Y11-P10
      A(2,4)=P10
      Y12=2.*Y(3)*Y(4)
      Y11=Y(3)*Y(3)
      A(3,1)=PAR
      A(3,3)=-7.+Y12-2.*PAR
      A(3,4)=Y11
      A(3,5)=PAR
      A(4,2)=P10
      A(4,3)=6.-Y12
      A(4,4)=-Y11-2.*P10
      A(4,6)=P10
      Y12=2.*Y(5)*Y(6)
      Y11=Y(5)*Y(5)
      A(5,3)=PAR
      A(5,5)=-7.+Y12-PAR
      A(5,6)=Y11
      A(6,4)=P10
      A(6,5)=6.-Y12
      A(6,6)=-Y11-P10

      A(1,7)=Y(3)-Y(1)
      A(2,7)=10.*(Y(4)-Y(2))
      A(3,7)=Y(1)+Y(5)-2.*Y(3)
      A(4,7)=10.*(Y(2)+Y(6)-2.*Y(4))
      A(5,7)=Y(3)-Y(5)
      A(6,7)=10.*(Y(4)-Y(6))
      RETURN
      END
