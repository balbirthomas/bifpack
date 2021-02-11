
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

c=====================================================================
c here: packB3


      SUBROUTINE PBP (N,M,T,PAR1,PAR2,INDL,INDK,IPR,BRANCH,IHMAX,
     1   FCN,BC,INTE,UNITLST)
      DOUBLE PRECISION T,Y,BRANCH,E,HA,EM,RS,G,A,D,BH,WORK,PAR1,PAR2
     1  , S,T1,TAU,H,T11,T2,HMAX,EA1,VFU1,TESTA,TESTB,TESTA2
     2  , ETA,EPS,ETA1,ETAI,EPMACH
      INTEGER UNITLST 
      PARAMETER (NDIM=12,NDIML=14)
      DIMENSION T(7),Y(NDIML,7),BRANCH(10,NDIML,7)
      DIMENSION E(NDIML,NDIML),HA(NDIML),EM(NDIML),RS(NDIML),
     1          G(NDIML,NDIML,7),A(NDIML,NDIML), D(NDIML),
     2          NPIV(NDIML),BH(NDIML),WORK(NDIML)
      EXTERNAL FCN,BC,INTE
      NROW=NDIML

C---------------------------------------------------------------------------
C--  ALGORITHM FOR CALCULATION OF PRIMARY BIFURCATION POINTS AND OF
C    SOLUTIONS THAT EMANATE FROM THE ZERO SOLUTION
C    IN BOUNDARY-VALUE PROBLEMS OF O.D.E.
C    AND FOR LINEAR EIGENVALUE PROBLEMS.
C
C    METHOD: THE VALUES OF A SPECIAL BIFURCATION TESTFUNCTION ARE CALCULATED
C    ------  FOR PARAMETERS BETWEEN  PAR1 (START) AND PAR2.
C            BASED ON INTERPOLATION OF THE TESTFUNCTION, AN APPROX. OF
C            THE PRIMARY BIFURCATION POINTS (EIGENVALUES) AND OF AN
C            EMANATING SOLUTION (EIGENFUNCTION) ARE EVALUATED.
C            THE ACCURACY DEPENDS ON THE AMPLITUDE OF THE TESTFUNCTION.
C            REFERENCE: R. SEYDEL, HABILITATION, TU MUENCHEN (81).
C
C            THE OUTPUT VALUES OF THIS SUBROUTINE ARE TO BE USED AS INITIAL
C            GUESS FOR SPECIAL BOUNDARY-VALUE PROBLEMS AS USED IN
C            CONTINUATION, OR (IF NEEDED) FOR A VERY PRECISE CALCULATION OF
C            THE BIFURCATION POINTS. (SEE R.SEYDEL: ISNM 48 (79), PP. 161-169)
C-----------------------------------------------------------------------------

C--  USER SUPPLIED SUBROUTINES:
C    -------------------------
C         FCN(T,Y,F)     THE RIGHT-HAND SIDE OF  DY=F(T,Y,PAR) .
C                        THE PARAMETER IS CONTAINED IN Y(N+1) ,
C                        F(N+1)=0  HAS TO BE DEFINED.
C         BC(YA,YB,R)    DEFINES THE BOUNDARY CONDITIONS. VECTOR  R,
C                        R=R(Y(A),Y(B))
C    STANDARD SOFTWARE:
C         INTE(N,FCN,TINIT,YINIT,TEND,TOLER,STEPMAX,STEPINI)
C                        INTEGRATOR FOR INITIAL-VALUE PROBLEMS
C                        DY=F(T,Y) ,  Y(TINIT)=YINIT  .
C         BGDECO, BGSOLV  FOR THE SOLUTION OF MATRIX-EQUATION.
C                        USE, FOR EXAMPLE, THE ALGORITHMS  DECOMP  AND
C                        SOLVE  DUE TO BUSINGER AND GOLUB, SEE
C                        WILKINSON, REINSCH: HANDBOOK

C--  INPUT VARIABLES:
C    ---------------
C         N              NUMBER OF DIFFERENTIAL EQS WITHOUT AUXILIARY VAR.
C         M              NUMBER OF NODES FOR SUBDIVISION OF INTEGRATION
C                        INTERVAL.
C         T              VECTOR OF NODES
C         PAR1, PAR2     DEFINE THE RANGE OF PARAMETER VALUES.  IN THIS
C                         PARAMETER INTERVAL BIFURCATION POINTS ARE SEARCHED.
C         INDL, INDK     TWO INDICES THAT FIX THE TESTFUNCTION.
C                        THE INDK- INITIAL VALUE OF Y(A) MUST NOT BE FIXED
C                        IN  BC  BY AN INITIAL CONDITION.
C         IPR .GT.0      RESULTS PRINTED
C         IHMAX          NUMBER OF PARAMETER VALUES FOR CHECKING THE
C                        BIFURCATION TESTFUNCTION

C--  OUTPUT VARIABLE:
C    ---------------
C         BRANCH    THIS ARRAY CONTAINS THE VALUES AT THE NODES  T(J) OF
C                   UP TO 10 EMANATING SOLUTIONS IN BRANCH(.,I,J), I=1,N .
C                   THESE VALUES APPROXIMATE THE SOLUTION OF THE LINEARIZAT
C                   THEY ARE NORMED TO   BRANCH(.,INDK,1)=0.01  .
C                   THE BIFURCATION VALUE OF THE PARAMETER IS CONTAINED
C                   IN  BRANCH(.,N+1) .
C         IPR=-1, IF PBP FAILS
C=============================================================================


      DDD=0.01d0
      PAR=PAR1
      PARDEL=(PAR2-PAR1)/DBLE(IHMAX)
      NP1=N+1
      IBIF=0
      DO 101 I=1,10
      DO 101 J=1,13
      DO 101 K=1,M
101   BRANCH(I,J,K)=0.d0
      IHOM2=0
      EPMACH=1.d-12
      HKZERO=1.d0
      EPS=1.d-6
      M1=M-1
      ETA=1.d-5
      ETAI=1.d0/ETA
      IHM=IHMAX*2



      DO 8 IHOM1=1,IHM
      IF (IPR.GT.0) WRITE (UNITLST,601) PAR
C--  CALCULATION OF A TESTFUNCTION AT THE BASIC SOLUTION
      IVERS=0
10    CONTINUE

      TAU=0.d0
      DO 11 I=1,NP1
11    HA(I)=0.d0
      D(NP1)=PAR
      HA(NP1)=PAR
      H=(T(M)-T(1))*0.05d0
      DO 1 J=1,M1
      T11=T(J)
      T2=T(J+1)
      HMAX=T2-T11
      DO 1 K=1,N
         DO 18 I=1,N
18       D(I)=0.d0
         D(K)=ETA
         IF (J.NE.1) GO TO 13
	     CALL BC (D,HA,RS)
             DO 12 I=1,N
12           A(I,K)=ETAI*RS(I)
13       T1=T11
         CALL INTE (NP1,FCN,T1,D,T2,EPS,HMAX,H)
         IF (H.NE.0.) GO TO 135
              WRITE (UNITLST,609) H,T1
              IPR=-1
              RETURN
135       IF (J.LT.M1) GO TO 15
	     CALL BC (HA,D,RS)
             DO 14 I=1,N
14           G(I,K,M1)=ETAI*RS(I)
             GO TO 1
15       DO 16 I=1,N
16       G(I,K,J)=ETAI*D(I)
1     CONTINUE



      DO 221 I=1,N
      DO 221 J=1,N
221   E(I,J)=G(I,J,1)
      IF (M.EQ.2) GO TO 21
      DO 2 K=2,M1
      DO 2 II=1,N
      DO 215 I=1,N
      HA(I)=0.d0
      DO 215 J=1,N
215   HA(I)=HA(I)+G(I,J,K)*E(J,II)
      DO 2 I=1,N
2     E(I,II)=HA(I)
21    DO 22 I=1,N
      DO 22 J=1,N
22    E(I,J)=A(I,J)+E(I,J)
      IF (INDK.LE.0.OR.INDL.LE.0) GO TO 28

      DO 25 J=1,N
      EA1=E(INDL,J)
      E(INDL,J)=0.d0
      EM(J)=EA1
25    RS(J)=0.d0
      RS(INDL)=HKZERO
      E(INDL,INDK)=1.d0
28    IRANK=N

      CALL BGDECO (N,N,E,NROW,D,NPIV,EPMACH,IRANK)
      IF (INDK.GT.0.AND.INDL.GT.0)
     1 CALL BGSOLV (N,N,E,NROW,D,NPIV,RS,HA,WORK)
      IF (IRANK.EQ.N) GO TO 3
      IF (IPR.GT.0) WRITE (UNITLST,604)
      IPR=-1
      RETURN

C     CALCULATION OF DETERMINANT
3     VFUD=D(1)
      DO 33 J=2,N
33    VFUD=VFUD*D(J)
      IF (INDK.GT.0.AND.INDL.GT.0)  GO TO 32
      TAU=VFUD
      IF (IPR.GT.0) WRITE (UNITLST,605) INDL,INDK,TAU,VFUD
      GO TO 39


C--  CALCULATION OF THE TESTFUNCTION:
32    DO 31 I=1,N
31    TAU=TAU+EM(I)*HA(I)
      IF (IPR.GT.0) WRITE (UNITLST,605) INDL,INDK,TAU,VFUD
      IF (IVERS.EQ.0) GO TO 39


C--  CALCULATION OF EMANATING SOLUTION
      DO 34 I=1,N
      Y(I,M)=0.d0
34    Y(I ,1)=HA(I)
      IF (M.LE.2) GO TO 39
      M2=M-2
      DO 36 J=1,M2
      J1=J+1
      DO 36 I=1,N
      S=0.d0
      DO 35 K=1,N
35    S=S+Y(K ,J)*G(I,K,J)
      Y(I,J1)=S
36    CONTINUE
39    CONTINUE
      IF (IVERS.EQ.1) GO TO 43
      VFU=TAU


      IHOM2=IHOM2+1
      IF (IHOM2.LE.1) TESTA =VFU
      IF (IHOM2.LE.1) PARALT= PAR
      IF (IHOM2.LE.1) TESTA2=ABS(TESTA)
C      CHECK FOR POLES:
      IF (ABS(TESTA).GT.TESTA2 .OR. IHOM2.LE.2) GO TO 7
      TESTB =TESTA *VFU
      IF  (TESTB.GT.0.d0) GO TO 7


C--  CALCULATION OF BIFURCATION POINT
      IBIF=IBIF+1
      VERZW=PARALT+( PAR-PARALT)*TESTA /(TESTA -VFU )
      WRITE (UNITLST,607) VERZW
      PARSP= PAR
      PAR=VERZW
      IVERS=1
      GO TO 10

43    PAR=PARSP
      VC=1.d-5*MAX(1.,REAL(VERZW))
      VD=VERZW-PAR
      IF (ABS(VD).LT.VC) GO TO 48
      VB=(TAU-VFU)/VD
      IF (ABS(PARALT-VERZW).LT.VC) GO TO 48
      VC=(TESTA-TAU)/(PARALT-VERZW)
      CONDB=ABS(VB+VC*(VERZW-PAR))
      CONDC=ABS(VC-VB)
      VC=(VC-VB)/(PARALT- PAR)
      VR1=-VB-VC*VD
      VR=VR1*VR1-4.d0*VC*TAU
      IF(VR.GT.0.d0) GO TO 44
         IF (IPR.GT.0.d0) WRITE (UNITLST,602)
         GO TO 48
44    VR=SQRT(VR)
      CONDA=VR
      IF (VR1.GT.0.) VR3=(VR1+VR)*0.5d0/VC
      IF (VR1.LE.0.) VR3=(VR1-VR)*0.5d0/VC
      VR1=VERZW+VR3
      IF (PARALT.LE.VR1.AND.VR1.LE.PAR.AND.PARALT.LT.PAR) GO TO 47
      IF (PAR.LE.VR1.AND.VR1.LE.PARALT.AND.PARALT.GT.PAR) GO TO 47
      VR3=TAU/(VC*VR3)
47    VERZW=VERZW+VR3
      WRITE (UNITLST,603) VERZW
      WRITE (*,603) VERZW
      VD=1.d-2
      IF (CONDA.GT.VD.AND.CONDB.GT.VD.AND.CONDC.GT.VD) GO TO 48
      IF (IPR.GT.0) WRITE (UNITLST,606)
48    CONTINUE


      IF (IBIF.GT.10) GO TO 69
      DO 62 J=1,M
      BRANCH(IBIF,NP1,J)=VERZW
      DO 62 I=1,N
62    BRANCH(IBIF,I,J)=Y(I,J)*DDD
69    IHOM2=0

7     TESTA2=ABS(TESTA)

      IF (IHOM2.GE.2) DERIV=(VFU-TESTA)/(PAR-PARALT)
      IF (IHOM2.LE.2) GO TO 71
      TESTNE=TESTNE-VFU
      TESTRE=ABS(TESTNE/VFU)
      IF (TESTRE.LT.0.05d0.AND.ABS(VFU).GT.TESTA2) PARDEL=PARDEL*1.2d0
      IF (TESTRE.GT.0.15d0.AND.ABS(VFU).LT.TESTA2) PARDEL=PARDEL*0.8d0
71    CONTINUE
      TESTA =VFU
      PARALT=PAR
      IF (IHOM2.LE.2) PAR=PAR+0.5d0*PARDEL
      IF (IHOM2.GT.2) PAR=PAR+PARDEL
      IF (IHOM2.GE.2) TESTNE=VFU+DERIV*(PAR-PARALT)
      IF (PAR2.GT.PAR1 .AND. PAR.GT.PAR2) RETURN
      IF (PAR2.LT.PAR1 .AND. PAR.LT.PAR2) RETURN


8      CONTINUE
      RETURN


601   FORMAT (' PARAMETER',E15.7)
602   FORMAT (' NO SIMPLE REAL ZERO, THIS SHOULD NOT HAPPEN')
603   FORMAT (' Primary bifurcation at parameter:',E15.7,10X)
604   FORMAT (' MATRIX SINGULAR')
605   FORMAT (' testfunction (',I2,'/',I2,') =', E13.5,
     1                      '  DET =',E10.3)
606   FORMAT (' BAD CONDITIONED, PROBABLY NO IMPROVEMENT OF GUESS')
607   FORMAT (/20H  bifurcation guess:,E15.7)
609   FORMAT (' H =',E15.5,'  AT  T =',E15.7)
      END
