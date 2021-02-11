      SUBROUTINE ODEBC1 (FCN,BC,METHOD,N,M,X,Y,EPS,ITMAX,J)
C
C   Multiple shooting code for ODE two-point Boundary value problems.
C   References for method, and convergence strategy:
C     R.Bulirsch, Report of the Carl-Cranz-Gesellschaft (1971).
C     P.Deuflhard: A Relaxation Strategy for the Modified Newton Method.
C          in: Lecture Notes of Math. Vol. 477, Springer (1975).
C     J.Stoer, R.Bulirsch: Introd. to Numerical Analysis. Springer (1980).
C     Convergence Strategy: P.Deuflhard (Numer.Math. 22 (1974) 289-315)

      DOUBLE PRECISION X,Y,D,U,UH,V,W,WH,Y1,A,E,DY,HH,HHA,YA,
     1  YU,YW,G,DY1,DY1A,EPDIFS,EPMACH,EPS,ETA,EPDMIN
     2  ,FC,FCMIN,FCMIN2,FCA,H,HX,HR,HMAX,HFIN,HR1
     3  ,S,T,ST,X1,X2
      PARAMETER (NDIML=14,MDIM=7, MDIM1=MDIM-1)
      DIMENSION X(MDIM),Y(NDIML,MDIM)
     1, D(NDIML),U(NDIML),UH(NDIML),V(NDIML),W(NDIML),WH(NDIML),
     2  Y1(NDIML), A(NDIML,NDIML),E(NDIML,NDIML)
     3, DY(NDIML,MDIM1),HH(NDIML,MDIM1),HHA(NDIML,MDIM1),YA(NDIML,MDIM)
     4, YU(NDIML,MDIM),YW(NDIML,MDIM),G(NDIML,NDIML,MDIM1)
     5, DY1(NDIML,MDIM1),DY1A(NDIML,MDIM1)
      INTEGER PIVOT(NDIML),UNITLST
C  Dimensions at least:   MDIM.ge.M ,  NDIML.ge.N
      COMMON /INTVAL/ JINT
      COMMON /BIF/ G,A
      EXTERNAL FCN
C
C  Unit number for printout:  UNITLST
      UNITLST=11
C
C  relative machine precision (to be adapted)
      EPMACH=1.d-13
C  minimum relative precision of integrator (to be adapted)
      EPDMIN=EPMACH*10.
C  prescribed relative precision in ode-integrator
      EPMACH=EPMACH*1.d3
      IF(EPS.LT.EPMACH) EPS=EPMACH
      EPDIFS=EPS*1.d-2
      IF(EPDIFS.LT.EPDMIN) EPDIFS=EPDMIN
C
C  prescribed relative deviation for numerical differentiation
      ETA=1.d-5
C
C  minimum permitted value of relaxation factor
      FCMIN=0.01
C  starting value of relaxation factor
      FC=FCMIN
C
C  decision parameter for rank-1 approximations (SIGMA.GT.1.)
C  no rank-1 approximation, if SIGMA.GT.1./FCMIN
      SIGMA=3.
C
C   initial preparations:
      FCMIN2=FCMIN*FCMIN
      ITER=0
      IFIN=0
      KT=0
      M1=M-1
      M2=M-2
      FCA=FC
      HX=(X(2)-X(1))*0.01
	 DO 101 I=1,N
	 T=0.
	    DO 1011 J=1,M
1011        T=T+ABS(Y(I,J))
	 T=T/FLOAT(M)
	 IF(T.LT.1.d0) T=1.
	 DO 101 J=1,M1
101      YW(I,J)=T
C
C   computation of the trajectories:
2     J=1
      KONV=1
      SABS=0.
      SUMJ=0.
      HR=HX
20    J1=J+1
      X1=X(J)
      X2=X(J1)
      H=HR
      HMAX=ABS(X2-X1)
      DO 201 K=1,N
201   U(K)=Y(K,J)
      JINT=J
      CALL METHOD(N,FCN,X1,U,X2,EPDIFS,HMAX,H)
      HFIN=HR
      HR=H
      IF(H.NE.0) GO TO 211
      IF(ITER.EQ.0) GO TO 93
      FC=FC*0.5
      GO TO 305
C
211   IF(IFIN.EQ.1) GO TO 900
      IF(J.EQ.M1) GO TO 23
C
C   CONTINUITY CONDITIONS:
      DO 22 K=1,N
      T=U(K)
      YU(K,J)=T
      T=T-Y(K,J1)
      HH(K,J)=T
      SABS=SABS+T*T
      T=T/YW(K,J1)
22    SUMJ=SUMJ+T*T
      J=J1
      GO TO 20
C
C   TWO-POINT BOUNDARY CONDITIONS:
23    DO 230 I=1,N
230   Y1(I)=Y(I,1)
      CALL BC(Y1,U,W)
      SUMR=0.
      DO 231 K=1,N
      T=W(K)
      HH(K,M1)=T
      SUMR=SUMR+T*T
231   Y(K,M)=U(K)
C
C  FIRST MONOTONICITY TEST
      SUM1=SUMJ+SUMR
      SABS=SABS+SUMR
      IF(ITER.EQ.0) THEN
		    SUMJA=SUMJ
		    GO TO 5
		    ENDIF
      KT=1
      IF(SUM1.LE.SUM1A) LEVEL=1
c      write (*,*) 'first test:',sum1,sabs
      GO TO 70
C
C   MODIFICATION OF NEWTON METHOD BY UNDERRELAXATION:
30    T=SQRT(SUM3/SUM3A)
      T=(T+FC-1.)*2./FC
      T=SQRT(4.*T+1.)-1.
      FC=FC/T
305   IF(FC.LT.FCMIN) GO TO 92
      DO 307 J=1,M1
      DO 307 I=1,N
307   Y(I,J)=YA(I,J)+FC*DY(I,J)
      GO TO 2
C
C   PREPARATIONS TO START THE FOLLOWING ITERATION STEP:
4     KT=0
      IF(FC.LT.FCA.AND.NEW.GT.0) LEVEL=0
      FCA=FC
      DO 44 I=1,N
      U(I)=Y(I,M)
      DO 44 J=1,M1
      T=(ABS(Y(I,J))+ABS(YA(I,J)))*0.5
      IF(T.LT.ABS(DY(I,J))) GO TO 44
      IF(T.LT.1.d0) T=1.d0
      YW(I,J)=T
44    DY1A(I,J)=DY1(I,J)
      IF(LEVEL.EQ.3) GO TO 6
C
C   DIFFERENCE APPROXIMATION OF THE JACOBIAN MATRIX:
5     NEW=0
      J=1
      HR=HX
      IF(KT.EQ.0) GO TO 50
	 DO 5001 I=1,N
	 W(I)=HH(I,M1)
	 U(I)=Y(I,M)
5001     Y1(I)=Y(I,1)
50    J1=J+1
      X11=X(J)
      X2=X(J1)
      HMAX=ABS(X2-X11)
      DO 500 I=1,N
      IF(J.NE.1.OR.I.LE.N) GO TO 502
      DO 501 K=1,N
      G(K,I,1)=0.
501   A(K,I)=0.
      GO TO 500
502   X1=X11
      H=HR
      DO 503 K=1,N
503   V(K)=Y(K,J)
      T=V(I)
      S=YW(I,J)*ETA
      IF(T.LT.0.) S=-S
      V(I)=T+S
      S=1./S
      IF(J.NE.1) GO TO 51
      CALL BC(V,U,WH)
      DO 504 K=1,N
504   A(K,I)=S*(WH(K)-W(K))
51    JINT=J
      CALL METHOD(N,FCN,X1,V,X2,EPDIFS,HMAX,H)
      HR1=H
      IF(H.EQ.0.) GO TO 931
      IF(J.EQ.M1) GO TO 53
      DO 52 K=1,N
52    G(K,I,J)=S*(V(K)-YU(K,J))
      GO TO 500
53    CALL BC(Y1,V,WH)
      DO 531 K=1,N
531   G(K,I,M1)=S*(WH(K)-W(K))
500   CONTINUE
      HR=HR1
      J=J1
      IF(J.LT.M) GO TO 50
      GO TO 7
C
C   APPROXIMATION OF THE JACOBIAN MATRIX BY RANK-1 CORRECTIONS:
6     NEW=NEW+1
      H=FCA-1.
      KBP=0
      IF(M.EQ.2) GO TO 62
      DO 61 I=1,N
      S=0.
      DO 611 K=1,N
611   S=S+A(I,K)*DY(K,1)
      IF(S.NE.0.) KBP=1
61    V(I)=S
      IF(KBP.EQ.0) GO TO 62
      DO 621 I=1,N
621   Y1(I)=YA(I,1)
      CALL BC(Y1,U,WH)
      DO 622 I=1,N
      T=WH(I)+FCA*V(I)
      WH(I)=T+H*HHA(I,M1)
622   V(I)=W(I)-T
62    DO 60 J=1,M1
      ST=0.
      DO 601 I=1,N
      T=DY(I,J)/YW(I,J)
      UH(I)=T
601   ST=ST+T*T
      ST=ST*FCA
      DO 60 K=1,N
      T=UH(K)/(ST*YW(K,J))
      IF(KBP.EQ.0) GO TO 600
      IF(J.NE.1) GO TO 6001
      DO 6002 I=1,N
      S=A(I,K)
      IF(S.NE.0.) A(I,K)=S+T*V(I)
6002  CONTINUE
6001  IF(J.NE.M1) GO TO 600
      DO 6003 I=1,N
      S=G(I,K,M1)
      IF(S.NE.0.) G(I,K,M1)=S+T*WH(I)
6003  CONTINUE
      GO TO 60
600   DO 6004 I=1,N
      S=G(I,K,J)
      IF(S.NE.0.) G(I,K,J)=S+T*(HH(I,J)+H*HHA(I,J))
6004  CONTINUE
60    CONTINUE
C
C   COMPUTATION OF THE FUNCTIONAL MATRIX E=-(A+BG(M1)*...*G(1)):
7     SABSA=SABS
      SUMRA=SUMR
      LEVEL=0
      KRED=0
	 DO 71 I=1,N
	 HHA(I,M1)=HH(I,M1)
	 DO 71 K=1,N
71       E(I,K)=G(I,K,M1)
      IF(M.EQ.2) GO TO 73
      DO 72 JJ=1,M2
      J=M1-JJ
      DO 72 I=1,N
	 DO 721 K=1,N
	 S=0.
	 DO 722 L=1,N
722      S=S+E(I,L)*G(L,K,J)
721      U(K)=S
      DO 723 K=1,N
723   E(I,K)=U(K)
72    HHA(I,J)=HH(I,J)
73    DO 74 L=1,N
      S=YW(L,1)
      DO 74 I=1,N
74    E(I,L)=-(A(I,L)+E(I,L))*S
C
C   COMPUTATION OF THE RIGHT-HAND SIDE VECTOR U(N):
70    DO 701 I=1,N
701   V(I)=HH(I,1)
      IF(M.EQ.2) GO TO 702
      DO 703 J=2,M1
	   DO 7031 I=1,N
	   S=HH(I,J)
	   DO 7032 L=1,N
7032       S=S+V(L)*G(I,L,J)
7031       U(I)=S
      DO 703 I=1,N
703   V(I)=U(I)
702   DO 704 I=1,N
704   U(I)=V(I)
      IF(KT.GT.0) GO TO 81
C
C   LSQ-SOLUTION OF THE LINEAR (N,N)-SYSTEM:
      CALL BGDECO (N,N,E,NDIML,D,PIVOT,EPMACH,IRANK)
C
81    IF(IRANK.LT.0) THEN
        WRITE (UNITLST,*) ' Iteration matrix is singular.'
        GO TO 92
        ENDIF 
C
      CALL BGSOLV (N,N,E,NDIML,D,PIVOT,U,V,UH)
C
810   SUM3=0.
      IF(KT.GT.0.OR.ITER.EQ.0) GO TO 811
      IF(KT.EQ.0) SFC1=0.
      SFC2=0.
811   DO 819 L=1,N
      T=V(L)
      SUM3=SUM3+T*T
      S=T*YW(L,1)
      IF(KT.LE.0) GO TO 813
      DY1(L,1)=S
      GO TO 819
813   IF(ABS(T).GT.EPS) KONV=0
      IF(ITER.EQ.0) GO TO 815
      T=T-DY1A(L,1)/YW(L,1)
      SFC2=SFC2+T*T
      IF(KT.LT.0) GO TO 815
      T=DY(L,1)/YW(L,1)
      SFC1=SFC1+T*T
815   DY(L,1)=S
      YA(L,1)=Y(L,1)
      YA(L,M)=Y(L,M)
819   U(L)=S
C
C   SECOND MONOTONICITY TEST:
      IF(KT.LE.0) GO TO 830
      SUM2=SUMJ+SUM3
c      write (*,*) 'second test:',sum2
      IF(SUM2.LE.SUM2A) LEVEL=LEVEL+1
      GO TO 83
830   SUM2A=SUM3
      IF(ITER.GT.0.AND.KT.EQ.0) SUMJA=0.
C
C   RECURSIVE SOLUTION OF M2 LINEAR (N,N)-SYSTEMS:
83    IF(M.EQ.2) GO TO 84
      DO 831 J=2,M1
      J1=J-1
      DO 832 I=1,N
      ST=HH(I,J1)
      S=ST
      DO 833 K=1,N
833   S=S+U(K)*G(I,K,J1)
      H=YW(I,J)
      T=S/H
      SUM3=SUM3+T*T
      IF(KT.LE.0) GO TO 8320
      DY1(I,J)=S
      GO TO 832
8320  IF(ABS(T).GT.EPS) KONV=0
      IF(ITER.EQ.0) GO TO 8321
      T=T-DY1A(I,J)/H
      SFC2=SFC2+T*T
      IF(KT.LT.0) GO TO 8321
      T=ST/H
      SUMJA=SUMJA+T*T
      T=DY(I,J)/H
      SFC1=SFC1+T*T
8321  DY(I,J)=S
      YA(I,J)=Y(I,J)
832   V(I)=S
      DO 831 I=1,N
831   U(I)=V(I)
C
C   THIRD MONOTONICITY TEST:
84    IF(KT.LE.0) GO TO 840
      IF(SUM3.LE.SUM3A) LEVEL=LEVEL+1
      IF(LEVEL) 30,30,4
840   SUM1A=SUMR+SUMJA
      SUM2A=SUM2A+SUMJA
      SUM3A=SUM3
c      write (*,*) 'third test:',sum1a,sum2a,sum3a
      SK=1.
      TD=0.
      TD=ABS(D(1))
      SK=TD/ABS(D(N))
C
C   ESTIMATION OF RELAXATION FACTOR:
      IF(ITER.GT.0) GO TO 85
      IF(KT.LT.0) FC=FCMIN
      GO TO 852
C
85    IF(SFC2.LT.SFC1*FCMIN2) SFC2=SFC1*FCMIN2
      FCS=SQRT(SFC1/SFC2)*FCA
C  MONITORING THE USE OF QUASI-NEWTON CORRECTIONS
      IF(NEW.EQ.0.OR.FCS.GE.FCA*SIGMA) GO TO 851
      KT=-1
      GO TO 5
851   IF(FCS.LT.FCMIN) GO TO 92
      FC=1.
      IF(FCS.LT.0.7) FC=FCS
852   ITER=ITER+1
      IF(EPDIFS*SK.GT.1.) EPDIFS=1./SK
      IF(EPDIFS.LT.EPDMIN) EPDIFS=1./TD
      IF(EPDIFS.LT.EPDMIN) EPDIFS=EPDMIN
c      write (*,*) 'fc/iter:',fc,iter
      DO 855 J=1,M1
      DO 855 I=1,N
855   Y(I,J)=YA(I,J)+FC*DY(I,J)
C
      IF(KONV.EQ.1) GO TO 90
      IF(ITER.GT.ITMAX) GO TO 91
      GO TO 2
C
C   FINAL COMPUTATION OF THE TRAJECTORY IN SUBINTERVAL M1:
90    IFIN=1
      J=M1
      HR=HFIN
      GO TO 20

C   SOLUTION EXIT:
900   DO 901 I=1,N
901   Y(I,M)=U(I)
      J=ITER
      if (level.lt.1.and.sum2.gt.eps.and.sum3.gt.eps) go to 94
      RETURN

C   FAIL EXITS:
91    J=-1
      WRITE (UNITLST,*) ' Iteration terminates after ITMAX=',ITMAX,
     , ' ITERATIONS.'
      GO TO 999
C
92    J=-2
      WRITE (UNITLST,*) ' Newton method fails to converge.'
      GO TO 999
C
931   WRITE (UNITLST,*) ' Singular trajectory by difference '
     1                  ,'approximation;'
93    J=-3
      J1=J1-1
      WRITE (UNITLST,*) ' integration terminates in subinterval ',J1,
     , ' AT T= ',X1
      GO TO 999
C
94    write (*,*) 'ODE1: iter,LEVEL,sum2,sum3=',j,level,sum2,sum3
      j=-4
999   RETURN
      END



      SUBROUTINE RKF7 (N,FCN,X,Y,XEND,TOL,HMAX,HI)
      DOUBLE PRECISION X,XEND,HI,HMAX,Y,W,HMIN,HMAG,H,XTRIAL
     1  ,TOL,EST,EST1
      PARAMETER (NDIML=14)
      DIMENSION Y(NDIML),W(NDIML,12)
C**** DIMENSION Y(N), W(N,12)  ***************
      INTEGER N,IND,I
      DATA UROUND/1.E-14/
C  UNIT ROUNDOFF IN REAL*4 FOR CYBER 175
      IF (HMAX.LT.TOL) HMAX=ABS(XEND-X)
      HMIN=UROUND/TOL
      HMIN=TOL*0.01*ABS(HMAX)
      IND=1
C
C        prepare to make calculations
C
100   GOTO (110,120,130), IND
110   CALL FCN(X,Y,W(1,1))
      HMAG=ABS(HI)
      GOTO 131
C     STEPSIZE ACCEPTED
120   CALL FCN(X,Y,W(1,1))
      HMAG=1.5E0*ABS(H)
      IF(EST.GT.(0.0279936E0*TOL)) HMAG=0.9E0*(TOL/EST)**(1./7.)*ABS(H)
      HMAG=MIN(HMAG,HMAX)
      GOTO 131
130   HMAG=0.9*(TOL/EST)**(1./7.)
      IF (HMAG.LT.0.5) HMAG=0.5
      HMAG=HMAG*ABS(H)
131   CONTINUE
      IF(HMAG.LT.ABS(XEND-X)) GOTO 140
      GOTO 141
140   H=SIGN(MIN(HMAG,0.5E0*ABS(XEND-X)),XEND-X)
      XTRIAL=X+H
      GOTO 142
141   H=XEND-X
      XTRIAL=XEND
142   CONTINUE
C
C        CALCULATIONS
C
      DO 200 I=1,N
200   W(I,12)=Y(I)+H*W(I,1)*(2.E0/27.E0)
      CALL FCN(X+H*(2.E0/27.E0),W(1,12),W(1,2))
      DO 201 I=1,N
201   W(I,12)=Y(I)+H*(W(I,1)*(1.E0/36.E0) + W(I,2)*(1.E0/12.E0) )
      CALL FCN(X+H*(1.E0/9.E0),W(1,12),W(1,3))
      DO 202 I=1,N
202   W(I,12)=Y(I)+H*(W(I,1)*(1.E0/24.E0) + W(I,3)*(1.E0/8.E0) )
      CALL FCN(X+H*(1.E0/6.E0),W(1,12),W(1,4))
      DO 203 I=1,N
203   W(I,12)=Y(I)+H*(W(I,1)*(5.E0/12.E0) - W(I,3)*(25.E0/16.E0)
     *                                    + W(I,4)*(25.E0/16.E0) )
      CALL FCN(X+H*(5.E0/12.E0),W(1,12),W(1,5))
      DO 204 I=1,N
204   W(I,12)=Y(I)+H*(W(I,1)*(1.E0/20.E0) + W(I,4)*(1.E0/4.E0)
     *                                    + W(I,5)*(1.E0/5.E0) )
      CALL FCN(X+H*0.5E0,W(1,12),W(1,6))
      DO 205 I=1,N
205   W(I,12)=Y(I)+H*(-W(I,1)*(25.E0/108.E0) + W(I,4)*(125.E0/108.E0)
     *               - W(I,5)*(65.E0/27.E0) + W(I,6)*(125.E0/54.E0) )
      CALL FCN(X+H*(5.E0/6.E0),W(1,12),W(1,7))
      DO 206 I=1,N
206   W(I,12)=Y(I)+H*(W(I,1)*(31.E0/300.E0) + W(I,5)*(61.E0/225.E0)
     *                 - W(I,6)*(2.E0/9.E0) + W(I,7)*(13.E0/900.E0) )
      CALL FCN(X+H*(1.E0/6.E0),W(1,12),W(1,8))
      DO 207 I=1,N
207   W(I,12)=Y(I)+H*(W(I,1)*2.E0 - W(I,4)*(53.E0/6.E0)
     *    + W(I,5)*(704.E0/45.E0) - W(I,6)*(107.E0/9.E0)
     *     + W(I,7)*(67.E0/90.E0) + W(I,8)*3.E0 )
      CALL FCN(X+H*(2.E0/3.E0),W(1,12),W(1,9))
      DO 208 I=1,N
208   W(I,12)=Y(I)+H*(-W(I,1)*(91.E0/108.E0) + W(I,4)*(23.E0/108.E0)
     *              - W(I,5)*(976.E0/135.E0) + W(I,6)*(311.E0/54.E0)
     *                - W(I,7)*(19.E0/60.E0) + W(I,8)*(17.E0/6.E0)
     *                                       - W(I,9)*(1.E0/12.E0) )
      CALL FCN(X+H*(1.E0/3.E0),W(1,12),W(1,10))
      DO 209 I=1,N
209   W(I,12)=Y(I)+H*(W(I,1)*(2383.E0/4100.E0) - W(I,4)*(341.E0/164.E0)
     *              + W(I,5)*(4496.E0/1025.E0) - W(I,6)*(301.E0/82.E0)
     *              + W(I,7)*(2133.E0/4100.E0) + W(I,8)*(45.E0/82.E0)
     *                 + W(I,9)*(45.E0/164.E0) + W(I,10)*(18.E0/41.E0) )
      CALL FCN(XTRIAL,W(1,12),W(1,11))
      DO 210 I=1,N
210   W(I,12)=Y(I)+H*(W(I,1)*(3.E0/205.E0) -W(I,6)*(6.E0/41.E0)
     *             - W(I,7)*(3.E0/205.E0) - W(I,8)*(3.E0/41.E0)
     *              + W(I,9)*(3.E0/41.E0) + W(I,10)*(6.E0/41.E0) )
      CALL FCN(X,W(1,12),W(1,2))
      DO 211 I=1,N
211   W(I,12)=Y(I)+H*(-W(I,1)*(1777.E0/4100.E0) - W(I,4)*(341.E0/164.E0)
     *               + W(I,5)*(4496.E0/1025.E0) - W(I,6)*(289.E0/82.E0)
     *               + W(I,7)*(2193.E0/4100.E0) + W(I,8)*(51.E0/82.E0)
     *                  + W(I,9)*(33.E0/164.E0) + W(I,10)*(12.E0/41.E0)
     *                                          + W(I,2) )
      CALL FCN(XTRIAL,W(1,12),W(1,3))
C
C        ESTIMATE ERROR
C
      EST=0.E0
      DO 300 I=1,N
      EST1= ABS((41.E0/840.E0)*(W(I,1)+W(I,11)-W(I,2)-W(I,3)))
300   EST=MAX(EST,EST1)
      IF (EST.LT.1.E-20) EST =1.E-20
C
C        DECISIONS
C
      IF(EST.GE.TOL) GOTO 402
      IND=2
      X=XTRIAL
      DO 401 I=1,N
401   Y(I)=Y(I)+H*(W(I,1)*(41.E0/840.E0) + W(I,6)*(34.E0/105.E0)
     *             + W(I,7)*(9.E0/35.E0) + W(I,8)*(9.E0/35.E0)
     *            + W(I,9)*(9.E0/280.E0) + W(I,10)*(9.E0/280.E0)
     *                                   + W(I,11)*(41.E0/840.E0) )
      IF(X.EQ.XEND) GOTO 406
      GOTO 405
C     IND=3  IF STEP TOO LARGE
402   IND=3
      IF( ABS(H).LE.HMIN) GOTO 403
      GOTO 404
403   IND=-2
      HI=0.
      RETURN
404   CONTINUE
405   CONTINUE
      GOTO 100
406   HI=SIGN(HMAG,H)
      RETURN
      END


      SUBROUTINE SOLVER (RHS,BC,INTE,N,M,T,Y,EPS,IT)
      EXTERNAL RHS,BC,INTE
      DOUBLE PRECISION T,Y,EPS
      INTEGER UNITLST
      PARAMETER (NDIML=14,MDIM=7)
      DIMENSION T(MDIM),Y(NDIML,MDIM)

C--  THIS SUBROUTINE HAS TO PROVIDE A SOLUTION TO A BOUNDARY-VALUE PROBLEM
C    of ordinary differential equations.
C    THE SUBROUTINE SERVES AS INTERFACE TO USE MULTIPLE SHOOTING.
C--
C    RHS    IS THE NAME OF THE SUBROUTINE DEFINING THE RIGHT-HAND SIDE
C    BC     IS THE NAME OF THE SUBR. DEFINING THE BOUNDARY CONDITIONS
C    N      NUMBER OF DIFFERENTIAL EQUATIONS
C    M      NUMBER OF SHOOTING NODES
C    T      VECTOR CONTAINING THE M NODES
C    Y      ARRAY CONTAINING APPROXIMATIONS TO THE SOLUTION
C           Y(.,J) VECTOR OF APPROXIMATION AT NODE T(J)
C           Y(I,.) VALUES OF I-TH COMPONENT AT THE NODES
C    EPS    (RELATIVE) PRECISION
C    INTE   INTEGRATOR OF INITIAL-VALUE PROBLEMS
C    IT      UPON INPUT:  MAXIMUM NUMBER OF ITERATIONS
C            UPON OUTPUT: IT.GT.0 IS NUMBER OF ITERATIONS PERFORMED
C                           .LT.0 POSSIBLY FAILURE
C                           =-1 : ITMAX ITERATIONS WERE NOT ENOUGH
C--
      ITMAX=IT
      UNITLST=11
      JPRINT=0
      IF (JPRINT.GT.0) THEN 
       write (UNITLST,*)
       write (UNITLST,*) ' initial data: T, Y(1), Y(2),...' 
       do 1 j=1,m
1      write (UNITLST,600) t(j),(y(i,j),i=1,n)
       endif
      J=0
      CALL ODEBC1 (RHS,BC,INTE,N,M,T,Y,EPS,ITMAX,J)
      IT=J
      IF (JPRINT.GT.0) THEN
       write (UNITLST,*) 'it=',it
       write (UNITLST,*) ' final data: T, Y(1), Y(2),...' 
       do 2 j=1,m
2      write (UNITLST,600) t(j),(y(i,j),i=1,n)
       endif
      RETURN
600   FORMAT (' ',5E15.8)
      END


      SUBROUTINE DATEN (F,METHOD,N,M,T,Y,EPS)                                  S
C  THIS ROUTINE CALCULATES THE VALUES OF THE Y-ARRAY,
c  Y(I,2), Y(I,3),...,Y(I,M)  FOR I=1,..,N,
C  IF THE INITIAL VALUES Y(I,1) ARE KNOWN.
C  GOOD IN CONNECTION WITH MULTIPLE SHOOTING.
C  Y(I,J) CONTAINS THE VALUE OF THE I-TH COMPONENT OF A SOLUTION
C  OF AN O.D.E. AT THE NODE T=T(J).

C  INPUT:
C  -----
C    F       NAME OF THE SUBROUTINE DEFINING RIGHT-HAND SIDE OF O.D.E.
C    METHOD  NAME OF INTEGRATION ROUTINE FOR INITIAL-VALUE PROBLEMS
C    N,M,T,Y (CLEAR FROMABOVE EXPLANATIONS)
C    EPS     RELATIVE PRECISION OF INITIAL VALUES Y(.,1)
C--------------------------------------------------------------------

      DOUBLE PRECISION T,Y,Y1,TA,TB,HMAX,HI,EPS,EPS2
      PARAMETER (NDIML=14)
      DIMENSION T(7),Y(NDIML,7),Y1(NDIML)
      EXTERNAL F
      EPS2=EPS*0.01
C    THE FOLLOWING INITIAL STEPSIZE MAY BE TOO LARGE FOR SOME
C    VERY SENSITIVE PROBLEMS:
      HI=(T(2)-T(1))*0.01
      J=0

1     J=J+1
      J1=J+1
      TA=T(J)
      TB=T(J1)
      HMAX=TB-TA
      DO 2 I=1,N
2     Y1(I)=Y(I,J)
      IF (HI.GT.HMAX) HI=HMAX
      
      CALL METHOD (N,F,TA,Y1,TB,EPS2,HMAX,HI)
      DO 5 I=1,N
5     Y(I,J1)=Y1(I)
      IF (J1.LT.M) GO TO 1
      RETURN
      END


      SUBROUTINE CHECKD (T,Y,NF,NY,FCN,FCNLIN,IFLAG,UNITLST)

C   THIS ROUTINE CHECKS (VIA NUMERICAL DIFFERENTIATION) WHETHER AN
C   ANALYTICAL JACOBIAN MATRIX IS CORRECTLY PROGRAMMED.

C   Y IS A NY-DIMENSIONAL INPUT VECTOR, MATRICES ARE EVALUATED AT Y.
C   T MAY BE DUMMY PARAMETER.

C   THE NF-DIM. VECTOR FUNCTION F(T,Y)  IS DEFINED IN
C   SUBROUTINE  FCN(T,Y,F)  
         
C   THE NF*NY PARTIAL DERIVATIVES
C                 A(I,J)=DF(I)/DY(J)
C   FOR  1.LE.I.LE.NF  ,  1.LE.J.LE.NY
C   ARE DEFINED IN SUBROUTINE
C             FCNLIN(T,Y,A) 

           
C   IFLAG (OUTPUT):
C   -------------   0 :  SUBROUTINE  FCNLIN  IS VERY LIKELY TO BE CORRECT.
C                .LT.0 :  SUBROUTINE  FCNLIN  MAY BE INCORRECT;
C                         THIS HOLDS FOR ABS(IFLAG) MANY ELEMENTS.
C                        CHECK THE PRINTED (CHANNEL 61) ELEMENT(S).
C=======================================================================


      DOUBLE PRECISION T,Y,F,A,B,D,D1,YST,DELTA,EPS,EPS1,DELTAI,
     1                 EPMACH 
      INTEGER UNITLST
      PARAMETER (NDIML=14)
      DIMENSION Y(NDIML),F(NDIML),A(NDIML,NDIML),B(NDIML)
      EXTERNAL FCN,FCNLIN

      EPMACH=1.E-12
C     printing of matrices, if IPRINT = 1 (no printing for IPRINT=0)
      IPRINT=1

      IFLAG=0
      EPS=EPMACH**0.3
      EPS1=EPMACH**0.25

      CALL FCN (T,Y,F)
      CALL FCNLIN (T,Y,A)

      IF (IPRINT.GT.0) WRITE (UNITLST,6666)
      DO 5 J=1,NY
      YST=Y(J)
      DELTA=EPS*MAX(1.,REAL(ABS(Y(J))))
      IF (Y(J).LT.0.) DELTA=-DELTA
      Y(J)=Y(J)+DELTA
      CALL FCN (T,Y,B)
      DELTAI=1./DELTA

      DO 2 I=1,NF
      B(I)=DELTAI*(B(I)-F(I))
      D=ABS(A(I,J)-B(I))
      D1=MAX(1.,REAL(ABS(B(J))))
      IF (D.LT.EPS1*D1) GO TO 2
      IFLAG =IFLAG-1
      IF (IFLAG.GE.-50) WRITE (UNITLST,666) I,J,B(I),A(I,J)
2     CONTINUE
      IF (IPRINT.GT.0) WRITE (UNITLST,6664) (B(I),I=1,NF)
5     Y(J)=YST

      IF (IFLAG.EQ.0) WRITE (UNITLST,6661)
      IF (IFLAG.EQ.0) WRITE (*,6661)
      IF (IFLAG.LT.0) THEN
             ICALL=-IFLAG
             WRITE (UNITLST,6662) ICALL
             WRITE (*,6662) ICALL
             END IF
      IF (IFLAG.LT.-50) WRITE (UNITLST,6663)
      IF (IPRINT.LE.0) RETURN
      WRITE (UNITLST,6665)
      DO 6 J=1,NY
6     WRITE (UNITLST,6664) (A(I,J),I=1,NF)
      RETURN

666   FORMAT (' CHECK ELEMENT (row=',I2,'/col=',I2,')  NUMER.:',E10.4,
     1              '  ANAL.:',E10.4)
6661  FORMAT (' ANALYTIC JACOBIAN SEEMS TO BE O.K.')
6662  FORMAT (' ',I5,' ELEMENTS OF ANAL. JACOBIAN MAY BE WRONG.')
6663  FORMAT ('        ONLY THE FIRST 50 ARE PRINTED.')
6664  FORMAT (' ',8E15.8)
6665  FORMAT (/' ANALYTIC JACOBIAN (transposed):')
6666  FORMAT (/' JACOBIAN (transposed) BY NUMERICAL DIFFER.:')
      END
