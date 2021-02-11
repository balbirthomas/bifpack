      SUBROUTINE EIGVAL (N,NM,B,WR,WI,IERR,INT) 
C   THIS SUBROUTINE CALLS THREE EISPACK-ROUTINES FOR
C   CALCULATING THE EIGENVALUES OF THE MATRIX  B .
C   Reference: B.T.Smith et al.: Matrix Eigensystem ... Springer 1974.
C              Lecture Notes in Computer Science 6,
C              in particular pages 150, 220, 240.
C   FOR COMMENTS SEE THE EISPACK-ROUTINES, ESPECIALLY HQR
C
C  INPUT:
C  -----     B  : MATRIX (WILL BE DESTROYED)
C            N  : DIMENSION
C            NM : NUMBER OF ROWS IN THE DIMENSION OF B
C
C  OUTPUT:
C  ------    WR : VECTOR CONTAINING REAL PARTS OF EIGENVALUES
C            WI : VECTOR CONTAINING IMAG.PARTS OF EIGENVALUES,
C                     BOTH DIMENSION  N
C************************************************************

      DOUBLE PRECISION B,WR,WI
      INTEGER INT,UNITLST
      DIMENSION B(NM,N),WR(N),WI(N),INT(NM)
      UNITLST=11
      CALL BALANC (NM,N,B,LOW,IGH,WI)
      CALL ELMHES (NM,N,LOW,IGH,B,INT)
      CALL HQR (NM,N,LOW,IGH,B,WR,WI,IERR)
      IF (IERR.GT.0) WRITE (UNITLST,601) IERR
      RETURN
601   FORMAT ('0EISPACK: The J-th eigenvalue is not determined, J='
     1       ,I3)
      END
      SUBROUTINE BALANC(NM,N,A,LOW,IGH,SCALE)                            EISPACK
C                                                                        EISPACK
      INTEGER I,J,K,L,M,N,JJ,NM,IGH,LOW,IEXC                             EISPACK
      DOUBLE PRECISION A,SCALE,C,F,G,R,S,B2,RADIX                        EISPACK
      DIMENSION A(NM,N),SCALE(N)                                         EISPACK
C     REAL ABS                                                           EISPACK
      LOGICAL NOCONV                                                     EISPACK
C                                                                        EISPACK
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE BALANCE,   EISPACK
C     NUM. MATH. 13, 293-304(1969) BY PARLETT AND REINSCH.               EISPACK
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 315-326(1971).    EISPACK
C                                                                        EISPACK
C     THIS SUBROUTINE BALANCES A REAL MATRIX AND ISOLATES                EISPACK
C     EIGENVALUES WHENEVER POSSIBLE.                                     EISPACK
C                                                                        EISPACK
C     ON INPUT-                                                          EISPACK
C                                                                        EISPACK
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL          EISPACK
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM           EISPACK
C          DIMENSION STATEMENT,                                          EISPACK
C                                                                        EISPACK
C        N IS THE ORDER OF THE MATRIX,                                   EISPACK
C                                                                        EISPACK
C        A CONTAINS THE INPUT MATRIX TO BE BALANCED.                     EISPACK
C                                                                        EISPACK
C     ON OUTPUT-                                                         EISPACK
C                                                                        EISPACK
C        A CONTAINS THE BALANCED MATRIX,                                 EISPACK
C                                                                        EISPACK
C        LOW AND IGH ARE TWO INTEGERS SUCH THAT A(I,J)                   EISPACK
C          IS EQUAL TO ZERO IF                                           EISPACK
C           (1) I IS GREATER THAN J AND                                  EISPACK
C           (2) J=1,...,LOW-1 OR I=IGH+1,...,N,                          EISPACK
C                                                                        EISPACK
C        SCALE CONTAINS INFORMATION DETERMINING THE                      EISPACK
C           PERMUTATIONS AND SCALING FACTORS USED.                       EISPACK
C                                                                        EISPACK
C     SUPPOSE THAT THE PRINCIPAL SUBMATRIX IN ROWS LOW THROUGH IGH       EISPACK
C     HAS BEEN BALANCED, THAT P(J) DENOTES THE INDEX INTERCHANGED        EISPACK
C     WITH J DURING THE PERMUTATION STEP, AND THAT THE ELEMENTS          EISPACK
C     OF THE DIAGONAL MATRIX USED ARE DENOTED BY D(I,J).  THEN           EISPACK
C        SCALE(J) = P(J),    FOR J = 1,...,LOW-1                         EISPACK
C                 = D(J,J),      J = LOW,...,IGH                         EISPACK
C                 = P(J)         J = IGH+1,...,N.                        EISPACK
C     THE ORDER IN WHICH THE INTERCHANGES ARE MADE IS N TO IGH+1,        EISPACK
C     THEN 1 TO LOW-1.                                                   EISPACK
C                                                                        EISPACK
C     NOTE THAT 1 IS RETURNED FOR IGH IF IGH IS ZERO FORMALLY.           EISPACK
C                                                                        EISPACK
C     THE ALGOL PROCEDURE EXC CONTAINED IN BALANCE APPEARS IN            EISPACK
C     BALANC  IN LINE.  (NOTE THAT THE ALGOL ROLES OF IDENTIFIERS        EISPACK
C     K,L HAVE BEEN REVERSED.)                                           EISPACK
C                                                                        EISPACK
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,         EISPACK
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY          EISPACK
C                                                                        EISPACK
C     ------------------------------------------------------------------ EISPACK
C                                                                        EISPACK
C     ********** RADIX IS A MACHINE DEPENDENT PARAMETER SPECIFYING       EISPACK
C                THE BASE OF THE MACHINE FLOATING POINT REPRESENTATION.  EISPACK
C                                                                        EISPACK
C                **********                                              EISPACK
      RADIX = 2.                                                         EISPACK
C                                                                        EISPACK
      B2 = RADIX * RADIX                                                 EISPACK
      K = 1                                                              EISPACK
      L = N                                                              EISPACK
      GO TO 100                                                          EISPACK
C     ********** IN-LINE PROCEDURE FOR ROW AND                           EISPACK
C                COLUMN EXCHANGE **********                              EISPACK
   20 SCALE(M) = J                                                       EISPACK
      IF (J .EQ. M) GO TO 50                                             EISPACK
C                                                                        EISPACK
      DO 30 I = 1, L                                                     EISPACK
         F = A(I,J)                                                      EISPACK
         A(I,J) = A(I,M)                                                 EISPACK
         A(I,M) = F                                                      EISPACK
   30 CONTINUE                                                           EISPACK
C                                                                        EISPACK
      DO 40 I = K, N                                                     EISPACK
         F = A(J,I)                                                      EISPACK
         A(J,I) = A(M,I)                                                 EISPACK
         A(M,I) = F                                                      EISPACK
   40 CONTINUE                                                           EISPACK
C                                                                        EISPACK
   50 GO TO (80,130), IEXC                                               EISPACK
C     ********** SEARCH FOR ROWS ISOLATING AN EIGENVALUE                 EISPACK
C                AND PUSH THEM DOWN **********                           EISPACK
   80 IF (L .EQ. 1) GO TO 280                                            EISPACK
      L = L - 1                                                          EISPACK
C     ********** FOR J=L STEP -1 UNTIL 1 DO -- **********                EISPACK
  100 DO 120 JJ = 1, L                                                   EISPACK
         J = L + 1 - JJ                                                  EISPACK
C                                                                        EISPACK
         DO 110 I = 1, L                                                 EISPACK
            IF (I .EQ. J) GO TO 110                                      EISPACK
            IF (A(J,I) .NE. 0.0) GO TO 120                               EISPACK
  110    CONTINUE                                                        EISPACK
C                                                                        EISPACK
         M = L                                                           EISPACK
         IEXC = 1                                                        EISPACK
         GO TO 20                                                        EISPACK
  120 CONTINUE                                                           EISPACK
C                                                                        EISPACK
      GO TO 140                                                          EISPACK
C     ********** SEARCH FOR COLUMNS ISOLATING AN EIGENVALUE              EISPACK
C                AND PUSH THEM LEFT **********                           EISPACK
  130 K = K + 1                                                          EISPACK
C                                                                        EISPACK
  140 DO 170 J = K, L                                                    EISPACK
C                                                                        EISPACK
         DO 150 I = K, L                                                 EISPACK
            IF (I .EQ. J) GO TO 150                                      EISPACK
            IF (A(I,J) .NE. 0.0) GO TO 170                               EISPACK
  150    CONTINUE                                                        EISPACK
C                                                                        EISPACK
         M = K                                                           EISPACK
         IEXC = 2                                                        EISPACK
         GO TO 20                                                        EISPACK
  170 CONTINUE                                                           EISPACK
C     ********** NOW BALANCE THE SUBMATRIX IN ROWS K TO L **********     EISPACK
      DO 180 I = K, L                                                    EISPACK
  180 SCALE(I) = 1.0                                                     EISPACK
C     ********** ITERATIVE LOOP FOR NORM REDUCTION **********            EISPACK
  190 NOCONV = .FALSE.                                                   EISPACK
C                                                                        EISPACK
      DO 270 I = K, L                                                    EISPACK
         C = 0.0                                                         EISPACK
         R = 0.0                                                         EISPACK
C                                                                        EISPACK
         DO 200 J = K, L                                                 EISPACK
            IF (J .EQ. I) GO TO 200                                      EISPACK
            C = C + ABS(A(J,I))                                          EISPACK
            R = R + ABS(A(I,J))                                          EISPACK
  200    CONTINUE                                                        EISPACK
C     ********** GUARD AGAINST ZERO C OR R DUE TO UNDERFLOW **********   EISPACK
         IF (C .EQ. 0.0 .OR. R .EQ. 0.0) GO TO 270                       EISPACK
         G = R / RADIX                                                   EISPACK
         F = 1.0                                                         EISPACK
         S = C + R                                                       EISPACK
  210    IF (C .GE. G) GO TO 220                                         EISPACK
         F = F * RADIX                                                   EISPACK
         C = C * B2                                                      EISPACK
         GO TO 210                                                       EISPACK
  220    G = R * RADIX                                                   EISPACK
  230    IF (C .LT. G) GO TO 240                                         EISPACK
         F = F / RADIX                                                   EISPACK
         C = C / B2                                                      EISPACK
         GO TO 230                                                       EISPACK
C     ********** NOW BALANCE **********                                  EISPACK
  240    IF ((C + R) / F .GE. 0.95 * S) GO TO 270                        EISPACK
         G = 1.0 / F                                                     EISPACK
         SCALE(I) = SCALE(I) * F                                         EISPACK
         NOCONV = .TRUE.                                                 EISPACK
C                                                                        EISPACK
         DO 250 J = K, N                                                 EISPACK
  250    A(I,J) = A(I,J) * G                                             EISPACK
C                                                                        EISPACK
         DO 260 J = 1, L                                                 EISPACK
  260    A(J,I) = A(J,I) * F                                             EISPACK
C                                                                        EISPACK
  270 CONTINUE                                                           EISPACK
C                                                                        EISPACK
      IF (NOCONV) GO TO 190                                              EISPACK
C                                                                        EISPACK
  280 LOW = K                                                            EISPACK
      IGH = L                                                            EISPACK
      RETURN                                                             EISPACK
C     ********** LAST CARD OF BALANC **********                          EISPACK
      END                                                                EISPACK


      SUBROUTINE ELMHES(NM,N,LOW,IGH,A,INT)                              EISPACK
C                                                                        EISPACK
      INTEGER I,J,M,N,LA,NM,IGH,KP1,LOW,MM1,MP1                          EISPACK
      DOUBLE PRECISION X,Y,A(NM,N)                                       EISPACK
C     REAL ABS                                                           EISPACK
      INTEGER INT(IGH)                                                   EISPACK
C                                                                        EISPACK
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE ELMHES,    EISPACK
C     NUM. MATH. 12, 349-368(1968) BY MARTIN AND WILKINSON.              EISPACK
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 339-358(1971).    EISPACK
C                                                                        EISPACK
C     GIVEN A REAL GENERAL MATRIX, THIS SUBROUTINE                       EISPACK
C     REDUCES A SUBMATRIX SITUATED IN ROWS AND COLUMNS                   EISPACK
C     LOW THROUGH IGH TO UPPER HESSENBERG FORM BY                        EISPACK
C     STABILIZED ELEMENTARY SIMILARITY TRANSFORMATIONS.                  EISPACK
C                                                                        EISPACK
C     ON INPUT-                                                          EISPACK
C                                                                        EISPACK
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL          EISPACK
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM           EISPACK
C          DIMENSION STATEMENT,                                          EISPACK
C                                                                        EISPACK
C        N IS THE ORDER OF THE MATRIX,                                   EISPACK
C                                                                        EISPACK
C        LOW AND IGH ARE INTEGERS DETERMINED BY THE BALANCING            EISPACK
C          SUBROUTINE  BALANC.  IF  BALANC  HAS NOT BEEN USED,           EISPACK
C          SET LOW=1, IGH=N,                                             EISPACK
C                                                                        EISPACK
C        A CONTAINS THE INPUT MATRIX.                                    EISPACK
C                                                                        EISPACK
C     ON OUTPUT-                                                         EISPACK
C                                                                        EISPACK
C        A CONTAINS THE HESSENBERG MATRIX.  THE MULTIPLIERS              EISPACK
C          WHICH WERE USED IN THE REDUCTION ARE STORED IN THE            EISPACK
C          REMAINING TRIANGLE UNDER THE HESSENBERG MATRIX,               EISPACK
C                                                                        EISPACK
C        INT CONTAINS INFORMATION ON THE ROWS AND COLUMNS                EISPACK
C          INTERCHANGED IN THE REDUCTION.                                EISPACK
C          ONLY ELEMENTS LOW THROUGH IGH ARE USED.                       EISPACK
C                                                                        EISPACK
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,         EISPACK
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY          EISPACK
C                                                                        EISPACK
C     ------------------------------------------------------------------ EISPACK
C                                                                        EISPACK
      LA = IGH - 1                                                       EISPACK
      KP1 = LOW + 1                                                      EISPACK
      IF (LA .LT. KP1) GO TO 200                                         EISPACK
C                                                                        EISPACK
      DO 180 M = KP1, LA                                                 EISPACK
         MM1 = M - 1                                                     EISPACK
         X = 0.0                                                         EISPACK
         I = M                                                           EISPACK
C                                                                        EISPACK
         DO 100 J = M, IGH                                               EISPACK
            IF (ABS(A(J,MM1)) .LE. ABS(X)) GO TO 100                     EISPACK
            X = A(J,MM1)                                                 EISPACK
            I = J                                                        EISPACK
  100    CONTINUE                                                        EISPACK
C                                                                        EISPACK
         INT(M) = I                                                      EISPACK
         IF (I .EQ. M) GO TO 130                                         EISPACK
C    ********** INTERCHANGE ROWS AND COLUMNS OF A **********             EISPACK
         DO 110 J = MM1, N                                               EISPACK
            Y = A(I,J)                                                   EISPACK
            A(I,J) = A(M,J)                                              EISPACK
            A(M,J) = Y                                                   EISPACK
  110    CONTINUE                                                        EISPACK
C                                                                        EISPACK
         DO 120 J = 1, IGH                                               EISPACK
            Y = A(J,I)                                                   EISPACK
            A(J,I) = A(J,M)                                              EISPACK
            A(J,M) = Y                                                   EISPACK
  120    CONTINUE                                                        EISPACK
C    ********** END INTERCHANGE **********                               EISPACK
  130    IF (X .EQ. 0.0) GO TO 180                                       EISPACK
         MP1 = M + 1                                                     EISPACK
C                                                                        EISPACK
         DO 160 I = MP1, IGH                                             EISPACK
            Y = A(I,MM1)                                                 EISPACK
            IF (Y .EQ. 0.0) GO TO 160                                    EISPACK
            Y = Y / X                                                    EISPACK
            A(I,MM1) = Y                                                 EISPACK
C                                                                        EISPACK
            DO 140 J = M, N                                              EISPACK
  140       A(I,J) = A(I,J) - Y * A(M,J)                                 EISPACK
C                                                                        EISPACK
            DO 150 J = 1, IGH                                            EISPACK
  150       A(J,M) = A(J,M) + Y * A(J,I)                                 EISPACK
C                                                                        EISPACK
  160    CONTINUE                                                        EISPACK
C                                                                        EISPACK
  180 CONTINUE                                                           EISPACK
C                                                                        EISPACK
  200 RETURN                                                             EISPACK
C    ********** LAST CARD OF ELMHES **********                           EISPACK
      END                                                                EISPACK
      SUBROUTINE HQR(NM,N,LOW,IGH,H,WR,WI,IERR)                          EISPACK
C                                                                        EISPACK
      INTEGER I,J,K,L,M,N,EN,LL,MM,NA,NM,IGH,ITS,LOW,MP2,ENM2,IERR       EISPACK
      DOUBLE PRECISION H,WR,WI,P,Q,R,S,T,W,X,Y,ZZ,NORM,MACHEP            EISPACK
      DIMENSION H(NM,N),WR(N),WI(N)                                      EISPACK
C     REAL SQRT,ABS,SIGN                                                 EISPACK
C     INTEGER MIN0                                                       EISPACK
      LOGICAL NOTLAS                                                     EISPACK
C                                                                        EISPACK
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE HQR,       EISPACK
C     NUM. MATH. 14, 219-231(1970) BY MARTIN, PETERS, AND WILKINSON.     EISPACK
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 359-371(1971).    EISPACK
C                                                                        EISPACK
C     THIS SUBROUTINE FINDS THE EIGENVALUES OF A REAL                    EISPACK
C     UPPER HESSENBERG MATRIX BY THE QR METHOD.                          EISPACK
C                                                                        EISPACK
C     ON INPUT-                                                          EISPACK
C                                                                        EISPACK
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL          EISPACK
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM           EISPACK
C          DIMENSION STATEMENT,                                          EISPACK
C                                                                        EISPACK
C        N IS THE ORDER OF THE MATRIX,                                   EISPACK
C                                                                        EISPACK
C        LOW AND IGH ARE INTEGERS DETERMINED BY THE BALANCING            EISPACK
C          SUBROUTINE  BALANC.  IF  BALANC  HAS NOT BEEN USED,           EISPACK
C          SET LOW=1, IGH=N,                                             EISPACK
C                                                                        EISPACK
C        H CONTAINS THE UPPER HESSENBERG MATRIX.  INFORMATION ABOUT      EISPACK
C          THE TRANSFORMATIONS USED IN THE REDUCTION TO HESSENBERG       EISPACK
C          FORM BY  ELMHES  OR  ORTHES, IF PERFORMED, IS STORED          EISPACK
C          IN THE REMAINING TRIANGLE UNDER THE HESSENBERG MATRIX.        EISPACK
C                                                                        EISPACK
C     ON OUTPUT-                                                         EISPACK
C                                                                        EISPACK
C        H HAS BEEN DESTROYED.  THEREFORE, IT MUST BE SAVED              EISPACK
C          BEFORE CALLING  HQR  IF SUBSEQUENT CALCULATION AND            EISPACK
C          BACK TRANSFORMATION OF EIGENVECTORS IS TO BE PERFORMED,       EISPACK
C                                                                        EISPACK
C        WR AND WI CONTAIN THE REAL AND IMAGINARY PARTS,                 EISPACK
C          RESPECTIVELY, OF THE EIGENVALUES.  THE EIGENVALUES            EISPACK
C          ARE UNORDERED EXCEPT THAT COMPLEX CONJUGATE PAIRS             EISPACK
C          OF VALUES APPEAR CONSECUTIVELY WITH THE EIGENVALUE            EISPACK
C          HAVING THE POSITIVE IMAGINARY PART FIRST.  IF AN              EISPACK
C          ERROR EXIT IS MADE, THE EIGENVALUES SHOULD BE CORRECT         EISPACK
C          FOR INDICES IERR+1,...,N,                                     EISPACK
C                                                                        EISPACK
C        IERR IS SET TO                                                  EISPACK
C          ZERO       FOR NORMAL RETURN,                                 EISPACK
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN                EISPACK
C                     DETERMINED AFTER 30 ITERATIONS.                    EISPACK
C                                                                        EISPACK
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,         EISPACK
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY          EISPACK
C                                                                        EISPACK
C     ------------------------------------------------------------------ EISPACK
C                                                                        EISPACK
C     ********** MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING      EISPACK
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.    EISPACK
C                                                                        EISPACK
C                **********                                              EISPACK
      MACHEP = 2.**(-53)                                                 EISPACK
C                                                                        EISPACK
      IERR = 0                                                           EISPACK
      NORM = 0.0                                                         EISPACK
      K = 1                                                              EISPACK
C     ********** STORE ROOTS ISOLATED BY BALANC                          EISPACK
C                AND COMPUTE MATRIX NORM **********                      EISPACK
      DO 50 I = 1, N                                                     EISPACK
C                                                                        EISPACK
         DO 40 J = K, N                                                  EISPACK
   40    NORM = NORM + ABS(H(I,J))                                       EISPACK
C                                                                        EISPACK
         K = I                                                           EISPACK
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 50                       EISPACK
         WR(I) = H(I,I)                                                  EISPACK
         WI(I) = 0.0                                                     EISPACK
   50 CONTINUE                                                           EISPACK
C                                                                        EISPACK
      EN = IGH                                                           EISPACK
      T = 0.0                                                            EISPACK
C     ********** SEARCH FOR NEXT EIGENVALUES **********                  EISPACK
   60 IF (EN .LT. LOW) GO TO 1001                                        EISPACK
      ITS = 0                                                            EISPACK
      NA = EN - 1                                                        EISPACK
      ENM2 = NA - 1                                                      EISPACK
C     ********** LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT              EISPACK
C                FOR L=EN STEP -1 UNTIL LOW DO -- **********             EISPACK
   70 DO 80 LL = LOW, EN                                                 EISPACK
         L = EN + LOW - LL                                               EISPACK
         IF (L .EQ. LOW) GO TO 100                                       EISPACK
         S = ABS(H(L-1,L-1)) + ABS(H(L,L))                               EISPACK
         IF (S .EQ. 0.0) S = NORM                                        EISPACK
         IF (ABS(H(L,L-1)) .LE. MACHEP * S) GO TO 100                    EISPACK
   80 CONTINUE                                                           EISPACK
C     ********** FORM SHIFT **********                                   EISPACK
  100 X = H(EN,EN)                                                       EISPACK
      IF (L .EQ. EN) GO TO 270                                           EISPACK
      Y = H(NA,NA)                                                       EISPACK
      W = H(EN,NA) * H(NA,EN)                                            EISPACK
      IF (L .EQ. NA) GO TO 280                                           EISPACK
      IF (ITS .EQ. 30) GO TO 1000                                        EISPACK
      IF (ITS .NE. 10 .AND. ITS .NE. 20) GO TO 130                       EISPACK
C     ********** FORM EXCEPTIONAL SHIFT **********                       EISPACK
      T = T + X                                                          EISPACK
C                                                                        EISPACK
      DO 120 I = LOW, EN                                                 EISPACK
  120 H(I,I) = H(I,I) - X                                                EISPACK
C                                                                        EISPACK
      S = ABS(H(EN,NA)) + ABS(H(NA,ENM2))                                EISPACK
      X = 0.75 * S                                                       EISPACK
      Y = X                                                              EISPACK
      W = -0.4375 * S * S                                                EISPACK
  130 ITS = ITS + 1                                                      EISPACK
C     ********** LOOK FOR TWO CONSECUTIVE SMALL                          EISPACK
C                SUB-DIAGONAL ELEMENTS.                                  EISPACK
C                FOR M=EN-2 STEP -1 UNTIL L DO -- **********             EISPACK
      DO 140 MM = L, ENM2                                                EISPACK
         M = ENM2 + L - MM                                               EISPACK
         ZZ = H(M,M)                                                     EISPACK
         R = X - ZZ                                                      EISPACK
         S = Y - ZZ                                                      EISPACK
         P = (R * S - W) / H(M+1,M) + H(M,M+1)                           EISPACK
         Q = H(M+1,M+1) - ZZ - R - S                                     EISPACK
         R = H(M+2,M+1)                                                  EISPACK
         S = ABS(P) + ABS(Q) + ABS(R)                                    EISPACK
         P = P / S                                                       EISPACK
         Q = Q / S                                                       EISPACK
         R = R / S                                                       EISPACK
         IF (M .EQ. L) GO TO 150                                         EISPACK
         IF (ABS(H(M,M-1)) * (ABS(Q) + ABS(R)) .LE. MACHEP * ABS(P)      EISPACK
     X    * (ABS(H(M-1,M-1)) + ABS(ZZ) + ABS(H(M+1,M+1)))) GO TO 150     EISPACK
  140 CONTINUE                                                           EISPACK
C                                                                        EISPACK
  150 MP2 = M + 2                                                        EISPACK
C                                                                        EISPACK
      DO 160 I = MP2, EN                                                 EISPACK
         H(I,I-2) = 0.0                                                  EISPACK
         IF (I .EQ. MP2) GO TO 160                                       EISPACK
         H(I,I-3) = 0.0                                                  EISPACK
  160 CONTINUE                                                           EISPACK
C     ********** DOUBLE QR STEP INVOLVING ROWS L TO EN AND               EISPACK
C                COLUMNS M TO EN **********                              EISPACK
      DO 260 K = M, NA                                                   EISPACK
         NOTLAS = K .NE. NA                                              EISPACK
         IF (K .EQ. M) GO TO 170                                         EISPACK
         P = H(K,K-1)                                                    EISPACK
         Q = H(K+1,K-1)                                                  EISPACK
         R = 0.0                                                         EISPACK
         IF (NOTLAS) R = H(K+2,K-1)                                      EISPACK
         X = ABS(P) + ABS(Q) + ABS(R)                                    EISPACK
         IF (X .EQ. 0.0) GO TO 260                                       EISPACK
         P = P / X                                                       EISPACK
         Q = Q / X                                                       EISPACK
         R = R / X                                                       EISPACK
  170    S = SIGN(SQRT(P*P+Q*Q+R*R),P)                                   EISPACK
         IF (K .EQ. M) GO TO 180                                         EISPACK
         H(K,K-1) = -S * X                                               EISPACK
         GO TO 190                                                       EISPACK
  180    IF (L .NE. M) H(K,K-1) = -H(K,K-1)                              EISPACK
  190    P = P + S                                                       EISPACK
         X = P / S                                                       EISPACK
         Y = Q / S                                                       EISPACK
         ZZ = R / S                                                      EISPACK
         Q = Q / P                                                       EISPACK
         R = R / P                                                       EISPACK
C     ********** ROW MODIFICATION **********                             EISPACK
         DO 210 J = K, EN                                                EISPACK
            P = H(K,J) + Q * H(K+1,J)                                    EISPACK
            IF (.NOT. NOTLAS) GO TO 200                                  EISPACK
            P = P + R * H(K+2,J)                                         EISPACK
            H(K+2,J) = H(K+2,J) - P * ZZ                                 EISPACK
  200       H(K+1,J) = H(K+1,J) - P * Y                                  EISPACK
            H(K,J) = H(K,J) - P * X                                      EISPACK
  210    CONTINUE                                                        EISPACK
C                                                                        EISPACK
         J = MIN0(EN,K+3)                                                EISPACK
C     ********** COLUMN MODIFICATION **********                          EISPACK
         DO 230 I = L, J                                                 EISPACK
            P = X * H(I,K) + Y * H(I,K+1)                                EISPACK
            IF (.NOT. NOTLAS) GO TO 220                                  EISPACK
            P = P + ZZ * H(I,K+2)                                        EISPACK
            H(I,K+2) = H(I,K+2) - P * R                                  EISPACK
  220       H(I,K+1) = H(I,K+1) - P * Q                                  EISPACK
            H(I,K) = H(I,K) - P                                          EISPACK
  230    CONTINUE                                                        EISPACK
C                                                                        EISPACK
  260 CONTINUE                                                           EISPACK
C                                                                        EISPACK
      GO TO 70                                                           EISPACK
C     ********** ONE ROOT FOUND **********                               EISPACK
  270 WR(EN) = X + T                                                     EISPACK
      WI(EN) = 0.0                                                       EISPACK
      EN = NA                                                            EISPACK
      GO TO 60                                                           EISPACK
C     ********** TWO ROOTS FOUND **********                              EISPACK
  280 P = (Y - X) / 2.0                                                  EISPACK
      Q = P * P + W                                                      EISPACK
      ZZ = SQRT(ABS(Q))                                                  EISPACK
      X = X + T                                                          EISPACK
      IF (Q .LT. 0.0) GO TO 320                                          EISPACK
C     ********** REAL PAIR **********                                    EISPACK
      ZZ = P + SIGN(ZZ,P)                                                EISPACK
      WR(NA) = X + ZZ                                                    EISPACK
      WR(EN) = WR(NA)                                                    EISPACK
      IF (ZZ .NE. 0.0) WR(EN) = X - W / ZZ                               EISPACK
      WI(NA) = 0.0                                                       EISPACK
      WI(EN) = 0.0                                                       EISPACK
      GO TO 330                                                          EISPACK
C     ********** COMPLEX PAIR **********                                 EISPACK
  320 WR(NA) = X + P                                                     EISPACK
      WR(EN) = X + P                                                     EISPACK
      WI(NA) = ZZ                                                        EISPACK
      WI(EN) = -ZZ                                                       EISPACK
  330 EN = ENM2                                                          EISPACK
      GO TO 60                                                           EISPACK
C     ********** SET ERROR -- NO CONVERGENCE TO AN                       EISPACK
C                EIGENVALUE AFTER 30 ITERATIONS **********               EISPACK
 1000 IERR = EN                                                          EISPACK
 1001 RETURN                                                             EISPACK
C     ********** LAST CARD OF HQR **********                             EISPACK
      END                                                                EISPACK
