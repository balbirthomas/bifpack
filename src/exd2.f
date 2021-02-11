      SUBROUTINE SETUP (MODEL,N,M,T,NAME1,JACOBI,EPS,
     1                  IFLAGSOL,WINDOW,YSTART)
      double precision eps,ystart,ETACO,ZETACO,PARCO,T,WINDOW
      parameter (NDIML=14)
      COMMON /CONTRO/ ETACO,ZETACO,PARCO, NCO,N1CO,N2CO,
     1                KCO,IND2CO,NSECFLAG,NFLAG
      dimension ystart(ndiml),t(7),window(2)
      CHARACTER*35 name1
      character*1 model 

C   Enter a short name of your problem:
      name1='exd2 CSTR:  c. stirred tank reactor'
      model='d'

C   If a particular solution is known, enter the solution into the
C   the vector YSTART, and set IFLAGSOL=1.
      iflagsol=1
      YSTART(1)=0.05747 
      YSTART(2)=0.23276 
      YSTART(3)=0.0483138 
      parco=16.2d0
      nsecflag=1
      window(1)=0.
      window(2)=100.
C    FIRST, SOME TECHNICAL DETAILS OF THE PARTICULAR
C    PROBLEM ARE INDICATED:

C  N IS THE NUMBER OF DIFFERENTIAL EQUATIONS:
      N=2
C   JACOBI=1  IF THE LINEARIZATION IS PROVIDED IN SUBROUTINE
C   FCNLIN,   ELSE TAKE  JACOBI=0
      JACOBI=1
C  THE RELATIVE ERROR TOLERANCE FOR ALL CALCULATIONS:
      EPS=1.d-4
      return
      end

      SUBROUTINE FCN (T,Y,F)
      DOUBLE PRECISION T,Y,F,PAR,ETACO,ZETACO,PARCO,
     1    B,beta,hilf
      COMMON /CONTRO/ ETACO,ZETACO,PARCO, NCO,N1CO,N2CO,
     1                KCO,IND2CO,NSECFLAG,NFLAG
      PARAMETER (NDIM=12,NDIML=14)
      DIMENSION Y(NDIML),F(NDIML)
      F(N1CO)=Y(KCO)-ZETACO*Y(IND2CO)-ETACO
C     THE PARAMETER IS GIVEN BY  Y(N+1)=Y(N1CO)
      PAR=Y(N1CO)
C     THIS SUBROUTINE DEFINES THE COMPONENTS OF THE RIGHT-HAND SIDE
C     F(1),...,F(n)   depending on the input   Y(1),...,Y(n),PAR  .
C     A second parameter may be defined in PARCO. 
C
C     example: CSTR

C parco=     B=16.2
      beta=3.d0
      hilf=par*(1.-y(1))*exp(y(2))
      F(1)=-y(1)+hilf
      F(2)=-y(2)+parco*hilf-beta*y(2)
      RETURN
      END


      SUBROUTINE FCNLIN (T,Y,A) 
      DOUBLE PRECISION T,Y,A,PAR,ETACO,ZETACO,PARCO,B,BETA,hilf
      PARAMETER (NDIM=12,NDIML=14)
      DIMENSION Y(NDIML),A(NDIML,NDIML)
      COMMON /CONTRO/ ETACO,ZETACO,PARCO, NCO,N1CO,N2CO,
     1                KCO,IND2CO,NSECFLAG,NFLAG
C  THIS ROUTINE ESTABLISHES PARTIAL DERIVATIVES OF RIGHT-HAND SIDE.
C     THE PARAMETER IS GIVEN BY  PAR=Y(N+1) (=Y(N1CO)):
      PAR=Y(N1CO)
C  OUTPUT IS  a(i,j)= partial derivative of f(i) with respect to y(j),
C                           i,j=1,...,n  .                       

      B=16.2d0
      beta=3.d0
      hilf=par*(1.-y(1))*exp(y(2))
      A(1,1)=-1.-par*exp(y(2))
      A(1,2)=hilf

      A(2,1)=-parco*par*exp(y(2))
      A(2,2)=-1.+parco*hilf-beta

      RETURN
      END

