      SUBROUTINE SETUP (MODEL,N,M,T,NAME1,JACOBI,EPS,IFLAGSOL,
     1                  WINDOW,YSTART)
      double precision eps,T,ystart,WINDOW,ETACO,ZETACO,PARCO
      parameter (NDIML=14)
      dimension T(7),ystart(ndiml),Window(2)
      CHARACTER*35 name1
      character*1 model 
      COMMON /CONTRO/ ETACO,ZETACO,PARCO, NCO,N1CO,N2CO,
     1                KCO,IND2CO,NSECFLAG,NFLAG

C   Enter a short name of your problem:
      name1='exb2 Superconductivity '
      model='b'
C   If a particular solution is known, enter the solution into the
C   the vector YSTART, and set IFLAGSOL=1.
      iflagsol=0
C    FIRST, SOME TECHNICAL DETAILS OF THE PARTICULAR
C    PROBLEM ARE INDICATED:
      nsecflag=1
      PARCO=1.
      M=7
      T(1)=0.
      T(2)=1.
      T(3)=2
      T(4)=2.5
      T(5)=3.
      t(6)=4. 
      T(M)=5. 
      window(1)=0.
      window(2)=1.2
c     example of a solution to start from:d
      ystart(1)=0.669836
      ystart(2)=0.
      ystart(3)=-1.3
      ystart(4)=1.
      ystart(5)=0.86597
      iflagsol=1 
C  N IS THE NUMBER OF DIFFERENTIAL EQUATIONS:
      N=4
C   JACOBI=1  IF THE LINEARIZATION IS PROVIDED IN SUBROUTINE
C   FCNLIN,   ELSE TAKE  JACOBI=0
      JACOBI=1
C  THE RELATIVE ERROR TOLERANCE FOR ALL CALCULATIONS:
      EPS=1.E-5
      return
      end

      SUBROUTINE FCN (T,Y,F)
      DOUBLE PRECISION T,Y,F,ETACO,ZETACO,PARCO,PAR
      COMMON /CONTRO/ ETACO,ZETACO,PARCO, NCO,N1CO,N2CO,
     1                KCO,IND2CO,NSECFLAG,NFLAG
      PARAMETER (NDIM=12,NDIML=14)
      DIMENSION Y(NDIML),F(NDIML)
      F(N1CO)=0.
C     THE BRANCHING PARAMETER IS GIVEN BY PAR=Y(N+1)=Y(N1CO)
      PAR=Y(5)
C     THIS SUBROUTINE DEFINES THE COMPONENTS OF THE RIGHT-HAND SIDE
C     F(1),...,F(n)   depending on the input   Y(1),...,Y(n),PAR  .
C
C     example: superconductivity
      f(1)=y(2)
      f(2)=y(1)*(y(1)*y(1)-PARCO+par*y(3)*y(3))
      f(3)=y(4)
      f(4)=y(1)*y(1)*y(3)   
      f(5)=0.
      RETURN
      END


      SUBROUTINE BC (YA,YB,R)
C   SUBROUTINE DEFINING BOUNDARY CONDITIONS
      DOUBLE PRECISION YA,YB,R,ETACO,ZETACO,PARCO,PAR
      PARAMETER (NDIM=12,NDIML=14)
      DIMENSION YA(NDIML), YB(NDIML),R(NDIML)
      COMMON /CONTRO/ ETACO,ZETACO,PARCO, NCO,N1CO,N2CO,
     1                KCO,IND2CO,NSECFLAG,NFLAG
C      PAR=YA(N1CO)
C  YA=y(a), YB=y(b), and PAR=parameter are input,
C  R(1),...,R(n) have to define n scalar boundary conditions:

      R(1)=ya(2)
      R(2)=yb(2)
      r(3)=ya(4)-1.
      r(4)=yb(4)-1. 
C  R(N+1)=YA(KCO)-ETA  , THE NEEDED VALUES ARE SUPPLIED VIA COMMON.
      R(5)=YA(KCO)-ETACO
      RETURN
      END


      SUBROUTINE FCNLIN (T,Y,A)     
C  SUBROUTINE CALCULATES PARTIAL DERIVATIVES
C  OF THE RIGHT-HAND SIDE  (N,N)-MATRIX
      DOUBLE PRECISION ETACO,ZETACO,PARCO,T,Y,A,PAR,Y12
      COMMON /CONTRO/ ETACO,ZETACO,PARCO, NCO,N1CO,N2CO,
     1                KCO,IND2CO,MSECFLAG,NFLAG
      PARAMETER (NDIM=12,NDIML=14)
      DIMENSION Y(NDIML),A(NDIML,NDIML)
C     THE PARAMETER IS CONTAINED IN Y(N+1)
      Y12=Y(1)*Y(1)
      PAR=Y(5)
      DO 1 I=1,4
      DO 1 J=1,4
1     A(I,J)=0.
      A(1,2)=1.
      A(2,1)=3.*Y12-PARCO+PAR*Y(3)*Y(3)
      A(2,3)=2.*PAR*Y(1)*Y(3)
      A(3,4)=1.
      A(4,1)=2.*Y(1)*Y(3)
      A(4,3)=Y12
      RETURN
      END

      SUBROUTINE BCLIN (YA,YB,R)  
      COMMON /CONTRO/ ETACO,ZETACO,PARCO, NCO,N1CO,N2CO,
     1                KCO,IND2CO,NSECFLAG,NFLAG
C  CONTAINS LINEARIZATION OF BOUNDARY CONDITIONS
      DOUBLE PRECISION YA,YB,R,ETACO,ZETACO,PARCO
      PARAMETER (NDIM=12,NDIML=14)
      DIMENSION YA(NDIML),YB(NDIML),R(NDIML)
      R(1)=YA(2)
      R(2)=YB(2)
      R(3)=YA(4)
      R(4)=YB(4)
      RETURN
      END
