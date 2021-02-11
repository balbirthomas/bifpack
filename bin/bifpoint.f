      PROGRAM MAIN
c     Program   bifpoint.f   , version 20.1.1996
c     to output all bifurcation points in file DATBASE

       OPEN (3,FILE ='DATBASE')
       OPEN (9,FILE= 'BIFPOINTS')
      call extract
      stop
      end


      SUBROUTINE EXTRACT
      DOUBLE PRECISION Y1,P,PERIOD
      character*300 buffer
c      character*35 name1,text1
      character*1 kezi1,label
      CHARACTER*2 KEZI2
C      integer lastb,ncountin,ncountout
      parameter (NDIML=14)
      dimension y1(ndiml),p(ndiml)

      do 2 loop=1,10000
      read (3,9000,end=9) buffer

      if (buffer(1:1).eq.'P') go to 2

      if (buffer(1:1).eq.'B') go to 2

      if (buffer(1:1).eq.'R'.or.buffer(1:1).eq.'r') then
         iperiod=0
         if (buffer(13:13).eq.'P'.or.
     1       buffer(13:13).eq.'H') iperiod=1
         if (iperiod.eq.0) read (buffer,9300) label,
     1       nob,nos,n,m,kezi1,kezi2,(y1(i),i=1,n),(p(i),i=1,m)
         if (iperiod.eq.1)
     1       read (buffer,9300) label,nob,nos,n,m,kezi1,kezi2,
     2       (y1(i),i=1,n),(p(i),i=1,m),period
  
      if (KEZI2.EQ.'HB' .OR. KEZI2.EQ.'PT'
     1    .OR. KEZI2.EQ.'T ' .OR. KEZI2.EQ.'PT'
     2    .OR. KEZI2.EQ.'PD' .OR. KEZI2.EQ.'PB'
     3    .OR. KEZI2.EQ.'BT') THEN
         IF (IPERIOD.EQ.0) WRITE (9,9300) LABEL,
     1      NOB,NOS,N,M,KEZI1,KEZI2,(Y1(I),I=1,N),(P(I),I=1,M)
         if (iperiod.eq.1)
     1      WRITE (9,9300) LABEL,NOB,NOS,N,M,KEZI1,KEZI2,
     2           (Y1(I),I=1,N),(P(I),I=1,M),PERIOD
         ENDIF
      ENDIF
2     continue
9     write (*,*) 'Bifpoints are in file BIFPOINTS. '
      return
9000  format (a300)
9300  format (a1,i2,i4,2i2,a1,a2,18e20.13)
      end

