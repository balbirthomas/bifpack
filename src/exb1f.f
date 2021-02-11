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
      name1='exb1 Catalyst '
      model='b'

C   If a particular solution is known, enter the solution into the
C   the vector YSTART, and set IFLAGSOL=1.
      iflagsol=1
      NSECFLAG=1
      PARCO=20. 
C    FIRST, SOME TECHNICAL DETAILS OF THE PARTICULAR
C    PROBLEM ARE INDICATED:
      M=5
      T(1)=0.
      T(2)=0.25
      T(3)=0.5
      T(4)=0.75
      T(M)=1. 
      window(1)=0.
      window(2)=0.2
      ystart(1)=0.9948478
      ystart(2)=0.
      ystart(3)=0.01
C  N IS THE NUMBER OF DIFFERENTIAL EQUATIONS:
      N=2
C   JACOBI=1  IF THE LINEARIZATION IS PROVIDED IN SUBROUTINE
C   FCNLIN,   ELSE TAKE  JACOBI=0
      JACOBI=1
C  THE RELATIVE ERROR TOLERANCE FOR ALL CALCULATIONS:
      EPS=1.E-4
      return
      end

      SUBROUTINE FCN (T,Y,F)
      DOUBLE PRECISION T,Y,F,ETACO,ZETACO,PARCO,PAR
      COMMON /CONTRO/ ETACO,ZETACO,PARCO, NCO,N1CO,N2CO,
     1                KCO,IND2CO,NSECFLAG,NFLAG
      PARAMETER (NDIM=12,NDIML=14)
      DIMENSION Y(NDIML),F(NDIML)
      F(N1CO)=0.
C     THE BIFURCATION PARAMETER IS GIVEN BY  PAR=Y(N+1)=Y(N1CO)
      PAR=Y(N1CO)
C     THIS SUBROUTINE DEFINES THE COMPONENTS OF THE RIGHT-HAND SIDE
C     F(1),...,F(n)   depending on the input   Y(1),...,Y(n),PAR  .
C
C     example: catalyst
      beta=0.4
c   parco=      gamma=20.
      F(1)=Y(2)
      RS=1.-y(1) 
      RS=beta*PARCO*RS/(1.+beta*RS)
      E=EXP(RS)
      F(2)=PAR*Y(1)*E
      RETURN
      END

      SUBROUTINE FCNLIN (T,Y,A) 
      DOUBLE PRECISION T,Y,A,ETACO,ZETACO,PARCO,PAR
      PARAMETER (NDIM=12,NDIML=14)
      DIMENSION Y(NDIML),A(NDIML,NDIML)
      COMMON /CONTRO/ ETACO,ZETACO,PARCO, NCO,N1CO,N2CO,
     1                KCO,IND2CO,NSECFLAG,NFLAG
C  THIS ROUTINE ESTABLISHES PARTIAL DERIVATIVES OF RIGHT-HAND SIDE.
C     THE PARAMETER IS GIVEN BY  PAR=Y(N+1) (=Y(N1CO)):
      PAR=Y(N1CO)
C  OUTPUT IS  a(i,j)= partial derivative of f(i) with respect to y(j),
C                           i,j=1,...,n  .                       

      beta=0.4
c  parco=     gamma=20.
      RS=1.-y(1) 
      E=EXP(beta*parco*RS/(1.+beta*RS))
      A(1,1)=0.
      A(1,2)=1.
      A(2,1)=PAR*E*(1.-y(1)*parco*BETA/((1.+BETA*RS)**2))
      A(2,2)=0.
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

      R(1)=yb(1)-1.
      R(2)=ya(2)
C  R(N+1)=YA(KCO)-ETA  , THE NEEDED VALUES ARE SUPPLIED VIA COMMON.
      R(N1CO)=YA(KCO)-ETACO
      RETURN
      END


      SUBROUTINE BCLIN (YA,YB,R)
      DOUBLE PRECISION YA,YB,R,ETACO,ZETACO,PARCO
      PARAMETER (NDIM=12,NDIML=14)
      DIMENSION YA(NDIML),YB(NDIML),R(NDIML)
      COMMON /CONTRO/ ETACO,ZETACO,PARCO, NCO,N1CO,N2CO,
     1                KCO,IND2CO,NSECFLAG,NFLAG
      R(1)=YB(1)
      R(2)=YA(2)
      RETURN
      END
