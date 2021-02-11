
c    BIFPACK -- a package for Bifurcation, Continuation  and Stability Analysis.
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

c========================================================================

c This elementary shooting routine is attached to BIFPACK for first tests.



      SUBROUTINE SHOOT (TDUMMY,Y,F)
c
c    latest update:  6 Jan 1996
C
C     This is a MULTIPLE SHOOTING device
C     for solving O.D.E. boundary-value problems.
C     This subroutine is designed to be called by NEWTON'S METHOD
C     or by a related nonlinear solver.

C     The right-hand side of the differential equation has to be
C     defined in SUBROUTINE  FCN ,  the boundary conditions in
C     SUBROUTINE  BC  (see the conventions in  B I F P A C K ).
C
C     INPUT VIA COMMON STATEMENTS:
C     ---------------------------
C           Nodes:  T(1),...,T(MNODES).
C           Total dimension for NEWTON's method:  N1CO,
C                 NCO1=N*(MNODES-1)+1 ,
C                 N : number of differential equations without parameter.
C           EPS : accuracy for integrator (smaller than the tolerance
C                                          of NEWTON's method)
C
C     STORAGE:
C     -------
C           Y(1),...,Y(N)  :  initial vector of O.D.E. sol. at T(1),
C           Y(N+1),...,Y(2N)           :  O.D.E.sol. at T(2),
C           Y(..),...,Y(N*(MNODES-1))  :  O.D.E. sol. at T(MNODES-1),
C           Y(N1CO)                    :  parameter.
C                      
C     FAIL EXIT (OUTPUT):
C     ------------------
C           IFLAG=1  if function in SHOOT cannot be evaluated /
C                       integration does not complete.
C                      
C      REQUIRED SUBROUTINES:
C      --------------------
C           RKF7   :  any O.D.E. integrator
C           FCN    :  right-hand side of O.D.E.
C           BC     :  boundary conditions
C                                               
C      COMMON- variables that are not mentioned here are dummy.    

      DOUBLE PRECISION Y,F,ETACO,ZETACO,PARCO,TDUMMY
     1  ,TA,TB,EPS,YB,T,HMAX,HI,F1
      PARAMETER (NDIML=14)
      DIMENSION Y(NDIML),F(NDIML),YB(NDIML),T(7),F1(NDIML)
      COMMON /CONTRO/ ETACO,ZETACO,PARCO, NCO,N1CO,N2CO,
     1                KCO,IND2CO,NSECFLAG,NFLAG
      COMMON /COSHOO/ T,EPS,MNODES
      COMMON /FFAIL/ IFLAG
      EXTERNAL FCN
      N=(N1CO-1)/(MNODES-1)
      N1=N+1
      HI=(T(2)-T(1))*0.01
      J=0

1     J=J+1
      TA=T(J)
      TB=T(J+1)
      HMAX=ABS(TB-TA)
      DO 2 I=1,N
2     YB(I)=Y( I+(J-1)*N )
      YB(N1)=Y(N1CO)
      CALL RKF7 (N1,FCN,TA,YB,TB,EPS,HMAX,HI)
      IF (TA.GE.TB) GO TO 3
      WRITE (*,*) ' SHOOTING FAILS. TERMINATED AT X=',TA
      WRITE (*,*) ' CONTINUATION STEP SIZE IS REDUCED.'
      IFLAG=1
      RETURN
3     IF (J.EQ.MNODES-1) GO TO 5
      DO 4 I=1,N
4     F(I+(J-1)*N)=Y(I+J*N)-YB(I)
      GO TO 1

5     CALL BC (Y,YB,F1)
      DO 6 I=1,N1
6     F(I+(J-1)*N)=F1(I)
      RETURN
      END


      SUBROUTINE RKF7 (N,FCN,X,Y,XEND,TOL,HMAX,HI)
      DOUBLE PRECISION X,XEND,HI,HMAX,Y,W,HMIN,HMAG,H,XTRIAL
     1  ,TOL,EST,EST1
      PARAMETER (NDIML=14)
      DIMENSION Y(NDIML),W(NDIML,12)
C**** DIMENSION Y(N), W(N,12)  ***************
      INTEGER N,IND,I
      DATA UROUND/1.E-14/
C  UNIT ROUNDOFF
      IF (HMAX.LT.TOL) HMAX=ABS(XEND-X)
      HMIN=UROUND/TOL
      HMIN=TOL*0.01*ABS(HMAX)
      if (HMIN.LT.0.0001) HMIN=0.0001
      IND=1
C
C        PREPARE TO MAKE CALCULATIONS
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
      if (abs(h).le.hmin) go to 403
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
