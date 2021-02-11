
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

c This nonlinear-equation solver is attached to BIFPACK for first tests.

      SUBROUTINE NOLE (N,X,FCN,EPS,ITMAX,J)
      PARAMETER (NDIML=14)
      INTEGER PIVOT(NDIML),UNITLST
      DOUBLE PRECISION X,Y,F,XY,A,D,Z,H,XNEU,FH,FH1,WORK
     1 , S,TEST,TESTA,X2,RESID,CORR,ALAM,EPS,EPMACH,ETA,CONST,
     2   DELTA,ACO,TDUMMY,LAMDA,LAMDAA,LAMDAM
      DIMENSION X(NDIML),Y(NDIML),F(NDIML),XY(NDIML),A(NDIML,NDIML),
     1          D(NDIML),Z(NDIML,NDIML),H(NDIML),
     2          XNEU(NDIML),FH(NDIML),FH1(NDIML),WORK(NDIML)
      COMMON /MATRIX/ Z
      COMMON /FFAIL/ IFLAG
      COMMON /INFONO/ IWEAK,NZAEHL,RESID,CORR
      NROW=NDIML
      IFLAG=0
C
C   J: FAIL EXIT
C   ------------
C        ON INPUT: J=-1 for no printing 
C        ON OUTPUT:
C           J=-1 : terminated after ITMAX iterations
C           J=-2 : no convergence: damping coefficient too small
C           J=-3 : function can not be evaluated. This flag is set if
C                  in the 'right-hand side' FCN the variable IFLAG
C                  (see COMMON /FFAIL/ ) is set =1.
C           J=-4 : matrix singular
C           J>0  : number of iterations
C=====================================================================
C    IT IS NOT WORTH LOOKING INTO THIS ROUTINE|
C    (THIS IS A PRELIMINARY, UNFINISHED EQUATION SOLVER)
C    YOU CERTAINLY PREFER YOUR OWN NEWTON SOLVER.
C=====================================================================

C      KT=1
      UNITLST=11
      NN=0
      EPMACH =1.d-13
      EPMACH=EPMACH*10.d0
      IF (EPS.LT.EPMACH) EPS=EPMACH
      ETA=1.d-6
      CONST=1.d-4
      KPRINT=J
      ITER=0
      LAMDA=1.d0
      TESTA=0.d0
      LAMDAM=0.02d0
      ALAM=0.5d0
      IF (KPRINT.GE.0) WRITE (UNITLST,604) (X(I),I=1,N)
      IF (KPRINT.GE.0) WRITE (UNITLST,601)
      NZAEHL=1
      CALL FCN (TDUMMY,X,F)
      IF (IFLAG.EQ.1) GO TO 86
      DO 25 I=1,N
      TESTA=TESTA+F(I)*F(I)
25    XY(I)=X(I) 
      IF (KPRINT.GE.0) WRITE (UNITLST,603) ITER,TESTA

C--------------------

3     KONV=0
      ITER=ITER+1
      IF (N.EQ.NN) CALL FCNLIN (TDUMMY,XY,A)
      IF (N.EQ.NN) GO TO 313

      NZAEHL=NZAEHL+1+N
      CALL  FCN (TDUMMY,XY,FH)
      IF (IFLAG.EQ.1) GO TO 86
      DO 31 J=1,N
      X2=XY(J)
      DELTA=ETA*X2
      IF (ABS(DELTA).LT.EPMACH) DELTA=ETA
      IF (X2.LT.0.) DELTA=-DELTA
      XY(J)=X2+DELTA
      CALL FCN (TDUMMY,XY,FH1)
      IF (IFLAG.EQ.1) GO TO 86
      DO 312 I=1,N
312   A(I,J)=(FH1(I)-FH(I))/DELTA
31    XY(J)=X2

313   CONTINUE
      DO 35 I=1,N
      DO 35 J=1,N
35    Z(I,J)=A(I,J)
      CALL BGDECO (N,N,A,NROW,D,PIVOT,EPMACH,ISING)
      IF (ISING.LT.0) THEN
         J=-4
         IF (KPRINT.GE.0) WRITE (UNITLST,*) ' NOLE: matrix is singular.'
         GO TO 9
         ENDIF
      ACO=ABS(D(1)/D(N))
      CALL BGSOLV (N,N,A,NROW,D,PIVOT,F,H,WORK)
      S=0.d0
      CORR=0.d0
      IWEAK=0
      DO 38 I=1,N
      S=S+XY(I)*XY(I)
38    CORR=CORR+H(I)*H(I)
      S=SQRT(S)
      CORR=SQRT(CORR)
      TOL=EPS*(1.d0+S)
      RESID=SQRT(TESTA)
      IF (RESID.LT.EPMACH .AND. CORR.LT.EPMACH) GO TO 801
      IF (RESID.LT.0.5*EPS .AND. CORR.LT.10.*TOL) IWEAK=1
      IF (RESID.LT.10.*EPS .AND. CORR.LT.0.5*TOL) IWEAK=1
      IF (RESID.LT.SQRT(EPS*EPMACH)) IWEAK=1
      
C     LAMDA-STRATEGIE

4     DO 41 I=1,N
41    XNEU(I)=X(I)-LAMDA*H(I)
      NZAEHL=NZAEHL+1
      CALL FCN (TDUMMY,XNEU,F)
      IF (IFLAG.EQ.1) GO TO 86
      TEST=0.d0
      DO 42 I=1,N
42    TEST=TEST+F(I)*F(I)
      IF (TEST.LT.TESTA) GO TO 43
      IF (KONV.EQ.1) GO TO 47
      GO TO 45
43    IF (KONV.LT.0.OR.LAMDA.GT.0.99) GO TO 48

      KONV=1
      TESTA =TEST
      LAMDAA=LAMDA
      LAMDA=LAMDA*2.
      DO 441 I=1,N
      XY(I)=XNEU(I)
441   FH1(I)=F(I)
      GO TO 4

45    KONV=-1
      IF (KPRINT.GE.0) WRITE (UNITLST,603) ITER,TEST ,LAMDA,CORR
      IF (LAMDA.LT.LAMDAM) GO TO 84
      LAMDA=LAMDA*0.5d0
      GO TO 4

47    DO 471 I=1,N
      XNEU(I)=XY(I)
471   F(I)=FH1(I)
      TEST=TESTA
      LAMDA=LAMDAA
48    TESTA=TEST
      KONV=0

C----------------

      DO 53 I=1,N
53    FH(I)=F(I)
      CALL BGSOLV (N,N,A,NROW,D,PIVOT,FH,H,WORK)
      
      CORR=0.d0
      DO 54 I=1,N
54    CORR=CORR+H(I)*H(I)
      RESID=SQRT(TEST)
      CORR=SQRT(CORR)
      IF (CORR.LT.TOL) KONV=KONV+1
      ALAM=0.5d0+RESID*CONST
      IF (ALAM.GT.1.d0) ALAM=1.d0


C     VORBER. DES NAECHSTEN ITER.-SCHRITTES

      DO 58 I=1,N
      X(I)=XNEU(I)
      Y(I)=XNEU(I)-H(I)*LAMDA
58    XY(I)=ALAM*X(I)+(1.d0-ALAM)*Y(I)
      IF (KPRINT.GE.0) WRITE(UNITLST,603) ITER,TEST,LAMDA,CORR,ALAM,ACO
      IF  (KPRINT.GT.0) WRITE (UNITLST,607) (X(I),Y(I),I=1,N)
      IF (KONV.GT.0 .AND. RESID.LT.EPS) GO TO 8
      IF (ITER.LE.ITMAX) GO TO 3
C-------------------

      J=-1
      IF (KPRINT.GE.0) WRITE (UNITLST,*) ' TOO MANY ITERATIONS'
      GO TO 9
84    IF (IWEAK.EQ.1) THEN
           IWEAK=-1
           IF (KPRINT.GE.0) WRITE (UNITLST,*) ' WEAK SOLUTION'
           IF (ITER.GT.1) GO TO 8
           GO TO 801
           END IF
      J=-2
      IF (KPRINT.GE.0) WRITE (UNITLST,*) 
     1                ' damping coefficient too small' 
      GO TO 9
86    J=-3
      IF (KPRINT.GE.0) WRITE (UNITLST,*) 
     1                ' function can not be evaluated.'
      GO TO 9

C     LOESUNG
8     DO 87 I=1,N
87    X(I)=Y(I)
801   IF (KPRINT.GE.0) WRITE (UNITLST,605)
      NZAEHL=NZAEHL+1
      CALL FCN (TDUMMY,X,F)
      RESID=0.d0
      S=0.d0
      DO 88 I=1,N
      S=S+X(I)*X(I)
88    RESID=RESID+F(I)*F(I)
      RESID=SQRT(RESID)
      S=SQRT(S)
      IF (S.GT.EPMACH) CORR=CORR/S
      IF (KPRINT.GE.0) WRITE(UNITLST,603) ITER,RESID,LAMDA,CORR,ALAM,ACO
      J=ITER
      IF (IWEAK.GE.0) IWEAK=0
      IF (IWEAK.LT.0) IWEAK=1


9     IF (KPRINT.GE.0) WRITE (UNITLST,604) (X(I),I=1,N)
      IF (KPRINT.GE.0) WRITE (UNITLST,608) NZAEHL
      RETURN
601   FORMAT (' ',4X,'ITER.', 8X,'TEST', 5X,'FACTOR', 8X,
     1    'KORR.',7X ,'MODIF.', 6X,'KOND.')
603   FORMAT (' ',I6,E18.8,F8.3, E15.5,F11.5,E12.3)
604   FORMAT  (' ',E25.13)
605   FORMAT (' SOLUTION', 7X,'RESIDUUM')
607   FORMAT (' ',2E25.13)
608   FORMAT (' ',I5,'  function calls')
      END
