      PROGRAM MAIN
      DOUBLE PRECISION EPS,EPSCO,T,Y,YSTART,WINDOW
      PARAMETER (NDIML=14)
      DIMENSION YSTART(NDIML),T(7),Y(NDIML),WINDOW(2)
      character*35 name1
      character*1 model
      INTEGER UNITLST,UNITOUT,UNIT1,UNIT3,UNIT31,UNIT4
      COMMON /COSHOO/ T,EPSCO,MNODES
      EXTERNAL SHOOT
      EXTERNAL FCN
      UNITLST=11
      UNITOUT=8
      UNIT1=1
      UNIT3=3
      UNIT31=31
      UNIT4=4
      OPEN (3 , FILE = 'START')
      OPEN (31, FILE = 'HEAD')
      OPEN (4,  FILE = 'BUFFER')
      OPEN (11, FILE = 'DATLIST')
      OPEN (8 , FILE = 'DATOUT')
      OPEN (1 , FILE = 'OPTIONS')

c     For convenience, we call the data provided by SETUP.
c     But some of the data will be overwritten :

      CALL SETUP (MODEL,N,M,T,NAME1,JACOBI,EPS,
     1           IFLAGSOL,WINDOW,YSTART)
      IF (NAME1.GT.' ') THEN
	WRITE (11,*) ' ============================'
	write (11,*) name1
	WRITE (11,*) ' ============================'
	WRITE (*,*) ' ============================'
	write (*,*) name1
	WRITE (*,*) ' ============================'
	endif

C     INTERACTIVE WORK WITH THE  PACK-A  SUBPACKAGE OF BIFPACK,
C     HERE VERSION THAT CALLS MULTIPLE SHOOTING.
C
C     EXAMPLE SUPERCONDUCTIVITY
C

C  PROBLEM DEPENDENT VARIABLES:
      EPSCO=EPS*0.01d0
C     ND= number of equations of the differential equation,
C     MDODES= number of shooting nodes, including the ends of interval,
C     NS= number of algebraic equations.
      ND=4
      MNODES=4
      NS=(MNODES-1)*ND
      T(1)=0.d0
      T(2)=2.5d0
      T(3)=4.d0
      T(MNODES)=5.d0

C     Initial solution, for instance (may be entered with MODUS=0):

C             Y(1),...,Y(4)    :  components of o.d.e. at T(1)=0.
      Y(1)=0.98237
      Y(2)=0.
      Y(3)=-1.
      Y(4)=1.

c             Y(5),...,Y(8)    :  components at T(2)=2.5
      Y(5)=0.997536
      Y(6)=0.
      Y(7)=0.
      Y(8)=0.16726

c             Y(9),...,Y(12)   :  components at T(3)=4.
      Y(9)=0.991186
      Y(10)=-0.0091779
      Y(11)=0.35515
      Y(12)=0.3906

c             y(13) = y(ns)    :  parameter
      Y(13)=0.0858193
      IFLAGSOL=1      

C     THEN WITH MODUS=2, YOU MAY ENTER / CHANGE CONTINUATION OPTIONS.
C     FOR INSTANCE, ENTER :
C       enter option No.:  4  ,  enter   30     (number of cont.steps)
C                          2             0      (parametrize by parameter)
C                          1             0.1    (initial step)
C                          3             2      (to vary steps freely) 
C                          6             0      (no stability analysis)  
C                          8  ,  enter   3  enter -2.4  (target value)
C      THEN run the continuation by entering  0 .

c**********************************************************
c*
c*    careful:  if you take the file PROBLEM.FOR to get FCN and BC,
c*              make sure to enter the    digit  5  for n+1 , and not N1CO,
c*              in  SUBROUTINE FCN :  PAR=Y(5)  ,
c*              in  SUBROUTINE BC  :  R(5)=YA(KCO)-ETACO  .
C*    This is important because n+1 and N1CO differ in subroutine SHOOT.
c*
c*********************************************************** 


      WRITE (*,*) ' Previous branching data to be purged?            '
      WRITE (*,*) '    Enter  0  for purging, enter  1  for appending: '
      READ (*,*) NPURGE
      IF (NPURGE.EQ.0) REWIND 8
      IF (NPURGE.EQ.1) THEN
                 DO 23 I=1,9999
23               READ (8,'( )',END=29)
29               CONTINUE
                 END IF

      CALL INTERA (NS,Y,IFLAGSOL,SHOOT,EPS,JACOBI,WINDOW,
     1              UNITLST,UNITOUT,UNIT1,UNIT3,UNIT31,UNIT4)
      write (*,*) 'Output on files DATOUT and DATLIST'
      STOP

      END
