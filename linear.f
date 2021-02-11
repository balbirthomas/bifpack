      SUBROUTINE BGDECO (M,N,A,NROW,D,PIVOT,ETA,ISING)
      INTEGER M,N,PIVOT(N),JBAR,I,J,K,K1
      DOUBLE PRECISION A,D,S,T,ETA,H
      DIMENSION A(NROW,N),D(N)
C
C--  BGDECO : Householder transformation of A(M,N), M.GE.N.
C             follows DECOMP by BUSINGER and GOLUB 
C                     see: WILKINSON, REINSCH: HANDBOOK ... p.113-114.
C                     Matrix A is destroyed.
C                     D: diagonal elements of reduced matrix.
C
      ISING=N
      DO     1  J=1,N
             PIVOT(J) = J
    1        CONTINUE
      DO     2  K=1,N
             IF  (K.EQ.N)  GO TO 22
             K1 = K+1

C     COLUMN SUMS:
             DO   210  J=K,N
                  S = 0.d0
                  DO  2101  I=K,M
                      S = S+A(I,J)*A(I,J)
 2101                 CONTINUE
                  D(J) = S
  210             CONTINUE
C
C     COLUMN PIVOTING:
C
             H = D(K)
             JBAR = K
             DO   212  J=K1,N
                  S = D(J)
                  IF  (S.LE.H)  GO TO 212
                  H = S
                  JBAR = J
  212             CONTINUE
C
C     COLUMN INTERCHANGE:
C
             IF (JBAR.NE.K) THEN
                I = PIVOT(K)
                PIVOT(K) = PIVOT(JBAR)
                PIVOT(JBAR) = I
                D(JBAR) = D(K)
                D(K)=H
                DO   215  I=1,M
                   T = A(I,K)
                   A(I,K) = A(I,JBAR)
                   A(I,JBAR) = T
  215              CONTINUE
                END IF

   22        H = 0.d0
             DO   221  I=K,M
                  H = H+A(I,K)*A(I,K)
  221             CONTINUE
             T = SQRT(H)
             IF (T.LE.ETA) THEN
                ISING=-1
                RETURN
                END IF
             S = A(K,K)
             IF  (S.GT.0.d0)  T = -T
             D(K) = T
             A(K,K) = S-T
             IF  (K.EQ.N)  RETURN
             T = 1.d0/(H-S*T)
             DO    24  J=K1,N
                   S = 0.d0
                   DO   241  I=K,M
                        S = S+A(I,K)*A(I,J)
  241                   CONTINUE
                   S = S*T
                   DO   242  I=K,M
                        A(I,J) = A(I,J)-A(I,K)*S
  242                   CONTINUE
                   D(J) = D(J)-A(K,J)*A(K,J)
   24              CONTINUE
    2        CONTINUE
      RETURN
      END


      SUBROUTINE BGSOLV (M,N,A,NROW,D,PIVOT,B,X,V)
C
C--  BGSOLV :  Subroutine SOLVE by BUSINGER and GOLUB .
C              Solution is stored in X[1:n] in correct order.
C
      INTEGER M,N,PIVOT(N)
      DOUBLE PRECISION A,B,D,X,V,S
      DIMENSION A(NROW,N),B(N),D(N),X(N),V(N)
C
C     HOUSEHOLDER TRANSFORMATIONS OF THE RIGHT-HAND SIDE:
C
      DO    31  J=1,N
            S = 0.d0
            DO   311  I=J,M
                 S = S+A(I,J)*B(I)
  311            CONTINUE
            S = S/(D(J)*A(J,J))
            DO   312  I=J,M
                 B(I) = B(I)+A(I,J)*S
  312            CONTINUE
   31       CONTINUE
C
C     SOLUTION OF THE UPPER TRIANGULAR SYSTEM:
C
      DO    41  II=1,N
            I = N+1-II
            S = B(I)
            IF  (II.EQ.1)  GO TO 410
            DO  4111  JJ=I1,N
                S = S-A(I,JJ)*V(JJ)
 4111           CONTINUE
  410       I1 = I
            V(I) = S/D(I)
   41       CONTINUE
C
C     BACK-PERMUTATION OF SOLUTION COMPONENTS:
C
      DO    50  J=1,N
            X(PIVOT(J)) = V(J)
   50       CONTINUE
      RETURN
      END
