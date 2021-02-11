      PROGRAM MAIN  
      PARAMETER (NDIM=12,NDIML=14)
      DIMENSION T(7),Y(NDIML,7),YA(NDIML,7),WINDOW(2),YSTART(NDIML)
      DOUBLE PRECISION T,Y,YA,EPS,S,STEPMA,STEPMI,YFINAL,WEIGHT
     1 ,WINDOW,STEPMP,ETACO,ZETACO,PARCO,BOUND2,YSTART
      CHARACTER*35 NAME1
      character*1 model
      COMMON /CONTRO/ ETACO,ZETACO,PARCO, NCO,N1CO,N2CO,
     1                KCO,IND2CO,NSECFLAG,NFLAG
      EXTERNAL FCN,BC,RKF7
      character*10 fileout,filelst
      integer UNIT1,UNIT2,UNIT3,UNIT31,UNIT4,UNITOUT,UNITLST
      CALL SETUP (model,n,m,t,name1,jacobi,eps,iflagsol,
     1            window,YSTART)
      CALL IOCONTROL (NAME1,UNIT1,UNIT2,UNIT3,UNIT31,UNIT4,
     *           UNITLST,UNITOUT,FILEOUT,FILELST)

C   INITIAL VECTOR must be in file START, together with the parameter
C   N+1 reals. 
C     TEST OF THE CONTINUATION ROUTINE  C O N T B   IN   P A C K B 1 .
      IF (name1.gt.' ') then
         WRITE (UNITOUT,9100) 'P',N,1,model,name1
	WRITE (UNITLST,*) '---------------------------'
	WRITE (UNITLST,*) name1
	WRITE (UNITLST,*) '---------------------------'
	WRITE (*,*) '---------------------------'
	WRITE (*,*) name1
	WRITE (*,*) '---------------------------'
	endif
C NOB, number of branch 
      WRITE (*,*) 'Enter number of branch:'
      read (*,*) NOB


C  DATA OF AN INITIAL SOLUTION TO START FROM: N+1 VALUES OF INITIAL VECTOR.  
C  IN GENERAL READ FROM FILE 'START' (UNIT 3).
C  (HERE DETERMINED BY INTEGRATING THE INITIAL-VALUE PROBLEM IN DATEN):
      N1=N+1
      N1CO=N1

      read (UNIT3,*,end=111)  (y(i,1),I=1,N1)
      GO TO 11
111   write (*,*) 'you need in file START ',n1,' reals.'
      if (iflagsol.eq.0) stop
      do 112 i=1,n1
112   y(i,1)=ystart(i)
      write (*,*) ' Data from YSTART in SETUPB are taken.'
11    CALL DATEN (FCN,RKF7,N1,M,T,Y,EPS)
      write (*,*) ' data array is set up.'
C  INITIAL STEP:  Y (A) := Y (A) + S
C                  K        K

C  VARIABLE STEPSIZES REQUIRED: LEVEL=2  (FIXED STEPSIZE: LEVEL=0)
      LEVEL=2

C  OPTIONS CONTROLLING THE CONTINUATION:
      MAXCON=100
      WEIGHT=1.
      ITAIM=5
      IEXTRA=1
      STEPMA=99999.
      STEPMI=EPS
      STEPMP=10000.
      LAST=0
      IPRINT=0
      ITMAX=10
      BOUND2=0.1
C  IN SEVERAL EXAMPLES YOU MAY WISH TO CHANGE THE MAXIMUM STEP LENGTH
C     STEPMA , OR YOU MAY LIKE TO SET A WINDOW FOR THE PARAMETER.
C     THAT IS, THE ABOVE OPTIONS ARE MORE OR LESS PROBLEM DEPENDENT.

C  EXAMPLE OF A FINAL CONDITION:  Y(2,1)=Y(IFINAL,1)=5.
c      IFINAL=2
c      YFINAL=5.

      IRUN=0

C***********************************************************************

C     THIS IS A SIMPLE DIALOG FOR SETTING BASIC OPTIONS
C     FOR INTERACTIVE WORK:

1     write (*,*) 'Enter  1  for a changed new continuation'
      WRITE (*,*) 'Enter  0  for continuation with previous step, OR '
      WRITE (*,*) 'Enter -1  for TERMINATING.    ENTER:           '
      READ (*,*) KCHOIC
      IF (KCHOIC.LT.0) GO TO 9
      IF (ABS(S).LT.EPS) KCHOIC=1 
      IF (KCHOIC.GT.0) THEN
           WRITE (*,*) ' Current parameter is ', Y(N1,1)
           WRITE (*,*) ' Parameter is stored in component No ', N1
           WRITE (*,*) ' Enter INDEX for next continuation step:'
           READ (*,*) K
           IF (K.LE.0) K=N1                                    
           WRITE (*,*) ' Enter STEP LENGTH:    '
           READ (*,*) S
           END IF
      WRITE (*,*) ' ENTER INDEX FOR TARGET POINT (0 FOR NO TARGET):  '
      READ (*,*) IFINAL
      IF (IFINAL.GT.0) THEN
           WRITE (*,*) ' ENTER TARGET VALUE    '
           READ (*,*) YFINAL
           END IF
      IRUN=IRUN+1

C***********************************************************************


2     CALL CONTB (NOB,NOS,N,M,T,FCN,BC,RKF7,Y,YA,LAST,EPS,
     1           LEVEL,K,S,STEPMA,STEPMI,STEPMP,MAXCON,
     2           ITAIM,ITMAX,WEIGHT,IEXTRA,
     3           IFINAL,YFINAL,IRUN,IPRINT,WINDOW,BOUND2,
     4           UNITLST,UNITOUT,UNIT3)
      IF (K.LT.0) WRITE (*,*) 'no success.'
      GO TO 1

9     WRITE (*,*) ' Last trajectory at parameter ', Y(N1,1)
      WRITE (*,601) 
601   FORMAT ('    t             y(1)               y(2)    ...  ')
      DO 91 J=1,M
91    WRITE (*,602) T(J),(Y(I,J),I=1,N)
602   FORMAT (' ',5E15.7)
      WRITE (*,*) ' Print data are in file ',filelst,' , '
      WRITE (*,*) ' summary of data in file ',fileout,' .  '
      STOP
9100  FORMAT (A1,2I2,A1,A35,2E20.13)
      END


      SUBROUTINE IOCONTROL (NAME1,UNIT1,UNIT2,UNIT3,UNIT31,UNIT4,
     *           UNITLST,UNITOUT,FILEOUT,FILELST)
C   This subroutine controls for input/output  the filenames,
C   and the unitnumbers.
C   This routine checks the leading characters of the string NAME1
C   that is provided by subroutine SETUP. In case a set of characters is
C   found with length at most 6, this substring is used as name for
C   output files, together with extension   .LST  (for a listing) ,
C   and extension   .OUT  (for the standard record of each solution).
C   If no such substring is found, the name will be "RESULT.LST", or
C   "RESULT.OUT", respectively.
C   For example, suppose your example is the Lorenz equation, and you
C   want to store all related files under the name  ex7  . Then begin
C   the string "name1" in SETUP with    ex7 Lorenz .....

      CHARACTER*35 name1
      character*10 fileout,filelst
      character*6 filename
      character*1 file1(10)
      integer UNIT1,UNIT2,UNIT3,UNIT31,UNIT4,UNITOUT,UNITLST
      read (name1,'(6a1)') (file1(i),i=1,6)
      if (file1(1).eq.' ') then
         filename='RESULT'
         read (filename,'(6a1)') (file1(i),i=1,6)
         endif
      logo=1
123   if (logo.eq.6 .or. file1(logo+1).eq.' ') go to 1234
      logo=logo+1
      go to 123
1234  write (filelst,'(10a1)') (file1(i),i=1,logo),'.','L','S','T'
      write (fileout,'(10a1)') (file1(i),i=1,logo),'.','O','U','T'
      UNITLST=11
      UNITOUT=8
      UNIT1=1
      UNIT2=2
      UNIT3=3
      UNIT31=31
      UNIT4=4
       OPEN (UNITLST, FILE = filelst)
       OPEN (UNITOUT, FILE = fileout)
       OPEN (UNIT1 ,  FILE = 'OPTIONS')
       OPEN (UNIT2,   FILE = 'INDFILE')
       OPEN (UNIT3 ,  FILE = 'START')
       OPEN (UNIT31,  FILE = 'HEAD')
       OPEN (UNIT4,   FILE = 'BUFFER')
      RETURN
      END
