      PROGRAM MAIN
      EXTERNAL RKF7
      DOUBLE PRECISION EPS,YSTART,T,WINDOW
      PARAMETER (NDIML=14)
      DIMENSION YSTART(NDIML),T(7),WINDOW(2)
      CHARACTER*35 name1
      character*10 fileout,filelst
      character*1 model
      integer UNIT1,UNIT2,UNIT3,UNIT31,UNIT4,UNITOUT,UNITLST
      CALL SETUP (model,n,m,t,name1,jacobi,eps,iflagsol,
     1            window,YSTART)
      CALL IOCONTROL (NAME1,UNIT1,UNIT2,UNIT3,UNIT31,UNIT4,
     *           UNITLST,UNITOUT,FILEOUT,FILELST)
      IF (name1.gt.' ') then
	WRITE (UNITLST,*) '---------------------------'
	WRITE (UNITLST,*) name1
	WRITE (UNITLST,*) '---------------------------'
	WRITE (*,*) '---------------------------'
	WRITE (*,*) name1
	WRITE (*,*) '---------------------------'
	endif

C   The following positions the file  ...OUT to its end in case
C   such a file is found.
      DO 21 I=1,20000
21    READ (UNITOUT,'( )',END=22) 
22    IF (I.GT.1) THEN
         BACKSPACE UNITOUT
         WRITE (*,*) 'Data will be appended to a '
     1              ,'previous file with name ',fileout,' .'
         ENDIF
      IF (I.LE.1) THEN
         REWIND UNITOUT         
         WRITE (UNITOUT,9100) 'P',n,1,model,name1
         ENDIF
      CALL BIFLAB (nob,nos,N,IFLAGSOL,YSTART,EPS,RKF7,JACOBI,WINDOW,
     1     UNIT2,UNIT3,UNIT31,UNIT4,UNITLST,UNITOUT )
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
