      SUBROUTINE SETUP (MODEL,N,M,T,NAME1,JACOBI,EPS,IFLAGSOL
     1                   ,WINDOW,YSTART) 
      double precision eps,T,ystart,window,ETACO,ZETACO,PARCO
      parameter (NDIML=14)
      dimension T(7),ystart(ndiml),window(2)
      CHARACTER*35 name1
      character*1 model
      COMMON /CONTRO/ ETACO,ZETACO,PARCO, NCO,N1CO,N2CO,
     1                KCO,IND2CO,NSECFLAG,NFLAG

C   Enter a short name of your problem:
      name1='exb4 Brusselator'
      model='b'
C   If a particular solution is known, enter the solution into the
C   the vector YSTART, and set IFLAGSOL=1.
      iflagsol=1
      NSECFLAG=0
C    FIRST, SOME TECHNICAL DETAILS OF THE PARTICULAR
C    PROBLEM ARE INDICATED:

C  Nodes (for shooting)
      M=5
      T(1)= 0.
      T(2)=0.25
      T(3)=0.5  
      T(4)=0.75
      T(5)= 1.
      window(1)=0.
      window(2)=100.
      ystart(1)=2.
      YSTART(2)=15.678
      ystart(3)=2.3
      ystart(4)=-4.8565
      ystart(5)=0.099886 
C   ANOTHER RUN ILLUSTRATING AN "OPEN" CONTINUATION MAY BE TRIED,
C   USING FOR INSTANCE
C          IFINAL=0
C          S=0.01
C          MAXCON=60
c      IF (INTER.EQ.0) GO TO 2C  N IS THE NUMBER OF DIFFERENTIAL EQUATIONS:
      N=4
C   JACOBI=1  IF THE LINEARIZATION IS PROVIDED IN SUBROUTINE
C   FCNLIN,   ELSE TAKE  JACOBI=0
      JACOBI=0
C  THE RELATIVE ERROR TOLERANCE FOR ALL CALCULATIONS:
      EPS=1.E-4
      return
      end



      SUBROUTINE FCN (T,Y,F)
      DOUBLE PRECISION T,Y,F,PAR,ETACO,ZETACO,PARCO
     1                 ,Y11,Y12,Y13,PAR2
      PARAMETER (NDIM=12,NDIML=14)
      COMMON /CONTRO/ ETACO,ZETACO,PARCO, NCO,N1CO,N2CO,
     1                KCO,IND2CO,NSECFLAG,NFLAG
      DIMENSION Y(NDIML),F(NDIML)
C  SUBROUTINE DEFINING THE RIGHT-HAND SIDE
C      THE PARAMETER IS GIVEN BY Y(N+1):
      PAR=Y(N1CO)
C  THIS IS THE TRIVIAL DIFFERENTIAL EQUATION:
      F(N1CO)=0.
C  NOW ENTER THE  N  VALUES OF THE RIGHT-HAND SIDE
C  F(1),...,F(n), depending on Y(1),...,Y(n), and PAR:

      Y11=Y(1)*Y(1)
      Y12=Y(1)*Y(3)
      Y13=Y12*Y(1)
      PAR2=PAR*PAR
      F(1)=Y(2)
      F(2)=-PAR2*(2.+Y13-5.6*Y(1))/0.0016
      F(3)=Y(4)
      F(4)=-PAR2*(4.6*Y(1)-Y13)/0.008

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

      R(1)=YA(1)-2.
      R(2)=YB(1)-2.
      R(3)=YA(3)-2.3
      R(4)=YB(3)-2.3

C  R(N+1)=YA(KCO)-ETA  , THE NEEDED VALUES ARE SUPPLIED VIA COMMON.
      R(N1CO)=YA(KCO)-ETACO
      RETURN
      END
