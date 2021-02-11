      PROGRAM MAIN
      DOUBLE PRECISION T,Y,Y1,TAU,TAU1,EPS,ETACO,ZETACO,PARCO,ETA,DE,
     1                 WINDOW,YSTART
      INTEGER UNITLST,UNITOUT
      PARAMETER (NDIM=12,NDIML=14)
      DIMENSION T(7),Y(NDIML,7),Y1(NDIML,7),WINDOW(2),YSTART(NDIML)
      CHARACTER*35 NAME1
      character*1 model
      EXTERNAL FCN,BC,RKF7
      COMMON /CONTRO/ ETACO,ZETACO,PARCO, NCO,N1CO,N2CO,
     1                KCO,IND2CO,NSECFLAG,NFLAG
      UNITLST=11
      UNITOUT=8
      OPEN (UNITLST, FILE='DATLIST')
      OPEN (UNITOUT, FILE='DATOUT')
      CALL SETUP (MODEL,N,M,T,NAME1,JACOBI,EPS,IFLAGSOL,WINDOW,YSTART)
      WRITE (11,*) name1
      WRITE (11,*) ' -------------------------'
      WRITE (*,*) name1
      WRITE (*,*) ' -------------------------'


C  TWO SPECIAL INDICES DEFINING THE TESTFUNCTION (FOR  T E S T B ):
      INDK=1
      INDL=1

C  DATA OF STARTING SOLUTION:
c  (example superconduction)   
      Y(2,1)=0.
      Y(4,1)=1.
      Y(1,1)=0.669836
      Y(3,1)=-1.3
      Y(5,1)=0.86596
      N1=N+1
      N1CO=n1
C      SUBROUTINE DATEN STARTS FROM THE INITIAL VALUES AND CALCULATES
C      INITIAL GUESS FOR  Y(I,J) , J.GE.2
      CALL DATEN (FCN,RKF7,N1,M,T,Y,EPS)

C     NUMBER of branch NOB, for instance:
      NOB=1

C--  THIS IS A VERY SIMPLE CONTINUATION, CHECKING THE SIGN OF THE
C    TESTFUNCTION THAT IS CALCULATED BY   T E S T B .
C    INDEX PICKS THE COMPONENT TO BE FIXED DURING CONTINUATION.
      DELTA=-0.02
      INDEX=1

      ETA=Y(INDEX,1)
      TAU=0.

c   the following five parameters need not be defined in case LEVEL=1:
      MO=1
      IF (JACOBI.EQ.0) MO=0
      ME=0
      MA=0
      DE=0.02
      IPRINT=1



      DO 6 ICON=1,6
      ETACO=ETA
      KCO=INDEX
      IF (ICON.EQ.1) GO TO 31 
      DO 3 J=1,M
      DO 3 I=1,N1+N
3     Y1(I,J)=Y(I,J)
31    TAU1=TAU
      J=0
      CALL ODEBC1 (FCN,BC,RKF7,N1,M,T,Y,EPS,20,J)
      IT=J
      IF (IT.LT.0) THEN
          WRITE (*,*) 'no convergence'
          stop
          endif
      WRITE (11,601) IT,Y(N1,1)
      WRITE (*,601) IT,Y(N1,1)
      WRITE (11,602)
      WRITE (11,603) (Y(I,1),I=1,N)

      J=1
      CALL TESTB (N,M,Y,J,INDL,INDK,TAU,N1,UNITLST)
      IF (ICON.EQ.1) GO TO 61
      IF (TAU*TAU1.GT.0.) GO TO 61

      LEVEL=1
      CALL SWITB (NOB,LEVEL,N,Y1,TAU1,Y,TAU,INDK,N1,M,T,EPS,RKF7,
     *            MO,ME,MA,DE,IPRINT,UNITLST,UNITOUT)

61    IF (IT.GT.0.and.NSECFLAG.eq.0) WRITE (8,9300)
     1    'R',NOB,ICON,N,1,'b','  ',(Y(I,1),I=1,N+1)
      IF (IT.GT.0.and.NSECFLAG.eq.1) WRITE (8,9300)
     1    'R',NOB,ICON,N,2,'b','  ',(Y(I,1),I=1,N+1),PARCO
6     ETA=ETA+DELTA
      WRITE (*,*) ' Print-out in file DATLIST, summary in DATOUT.   '
      STOP
601   FORMAT (' solution within',I3,' iterations.    PARAMETER:',
     1         E15.7,10X)
602   FORMAT (' vector of initial values:')
603   FORMAT (' ',5E15.7)
9300  FORMAT (A1,I2,I4,2I2,A1,A2,13E20.13)
      END
