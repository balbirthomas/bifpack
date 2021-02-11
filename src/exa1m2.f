      PROGRAM MAIN
      PARAMETER (NDIM=12,NDIML=14)
      DOUBLE PRECISION T(7),WINDOW(2),T1,Y(NDIML),F(NDIML),YSP(NDIML)
     1   ,YSTART(NDIML),VFU,VFUD,VFUSP,VFUALT,EPS,ETACO,ZETACO
     2   ,PARCO,DELTA
      EXTERNAL FCN,FCNLIN
      CHARACTER*35 NAME1
      character*1 model
      INTEGER UNITLST,UNITOUT,UNIT3
      COMMON /CONTRO/ ETACO,ZETACO,PARCO, NCO,N1CO,N2CO,
     1                KCO,IND2CO,NSECFLAG,NFLAG

C     SWITCHING BRANCHES IN NONLINEAR EQUATIONS
C       THIS IS TO TEST BRANCH SWITCHING ABILITIES OF  P A C K A 2 .
C       
      CALL SETUP (MODEL,N,M,T,NAME1,JACOBI,EPS,
     1           IFLAGSOL,WINDOW,YSTART)
      CALL IOCONTROL (NAME1,UNIT1,UNIT2,UNIT3,UNIT31,UNIT4,
     *           UNITLST,UNITOUT,FILEOUT,FILELST)
      OPEN (UNITLST , FILE = 'exa1m2.LST')
      OPEN (UNIT3  , FILE = 'START')
      OPEN (UNITOUT  , FILE = 'exa1m2.OUT')
      WRITE (UNITLST,*) ' Testing Branch switching, EXA1 '
      WRITE (UNITLST,*) ' -------------------------------'



C----  Define a SIMPLE CONTINUATION:
C    (INDEX DETERMINES THE COMPONENT TO CONTROL THE CONTINUATION)
      DELTA=0.03
      IHMAX=5
      INDEX=7
      INDL=5
      INDK=1

C     INITIAL SOLUTION:
      Y(1)=1.4458
      Y(2)=2.7436
      Y(3)=3.1084
      Y(4)=2.539
      Y(5)=1.446
      Y(6)=2.744
      Y(7)=1.435
      rewind UNIT3
      do 1 i=1,n+1
      write (UNIT3,*) y(i)
1     continue
      rewind UNIT3
c      ETA=Y(INDEX)
c----------------------------------------------------
      call locate (INDL,INDK,INDEX,DELTA,IHMAX,
     1             N,UNITLST,UNITOUT,UNIT3,JACOBI)

      write (*,*) 'Printout on file  exa1m2.LST, and ...OUT'
      STOP
601   FORMAT (' ',6E13.6)
602   FORMAT (/' POLE')
603   FORMAT ('0SOLUTION AFTER',I3,'  ITER.,  PARAMETER:',E14.7)
9300  FORMAT (a1,i2,i4,2i2,a1,a2,13e20.13)
      end


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
