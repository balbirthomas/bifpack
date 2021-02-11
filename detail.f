      PROGRAM MAIN
c     Program   detail.f   , version 20.1.1996.
c     "data manager" for calculating a detail:
c     Only those data are written to DATBASENEW that match the interval
c     PARLEFT <= PAR  <= PARRIGHT

       OPEN (3,FILE ='DATBASE')
       OPEN (14,FILE='DATBASEN')
      write (*,*) 'Enter left end of PAR-interval:'
      read (*,*) parleft
      write (*,*) 'enter right end of PAR-interval:'
      read (*,*) parright
      call detail (parleft,parright)
      stop
      end


      SUBROUTINE DETAIL (PARLEFT,PARRIGHT)
      DOUBLE PRECISION Y1,P,PERIOD
      character*300 buffer
      character*1 kezi1,label
      character*2 kezi2
      parameter (NDIML=14)
      dimension y1(ndiml),p(ndiml)
      Iraus=0
      irein=0
      do 2 loop=1,10000
      read (3,9000,end=9) buffer

      if (buffer(1:1).eq.'P') go to 2
      if (buffer(1:1).eq.'B') go to 2

      if (buffer(1:1).eq.'R'.or.buffer(1:1).eq.'r') then
          irein=irein+1
	 iperiod=0
	 if (buffer(13:13).eq.'P'.or.
     1       buffer(13:13).eq.'H') iperiod=1
	 if (iperiod.eq.0) read (buffer,9300) label,
     1      nob,nos,n,m,kezi1,kezi2,(y1(i),i=1,n),(p(i),i=1,m)
	 if (iperiod.eq.1)
     1      read (buffer,9300) label,nob,nos,n,m,kezi1,kezi2,
     2           (y1(i),i=1,n),(p(i),i=1,m),period
         if (p(1).gt.parright) go to 2
         if (p(1).lt.parleft) go to 2
         if (iperiod.eq.0) write (14,9300) label,
     1      nob,nos,n,m,kezi1,kezi2,(y1(i),i=1,n),(p(i),i=1,m)
	 if (iperiod.eq.1)
     1      write (14,9300) label,nob,nos,n,m,kezi1,kezi2,
     2           (y1(i),i=1,n),(p(i),i=1,m),period
            iraus=iraus+1
          endif
2     continue
9     continue
      write (*,*) 'The required DETAIL is written to DATBASEN ;'
      write (*,*) iraus,' records selected out of ',irein
      return
9000  format (a300)
9001  format (60('='))
9002  format ('at parameter ',e11.5,a20)
9003  format ('at parameter ',e11.5,a20,'. Period:',e11.5)
9100  format (a1,2i2,a1,a35,2e20.13)
9101  format ('Problem : ',a35)
9102  format (' Dimensions: ',i2,' eqs. and ',i2,' parameters')
9103  format ('   ',a35)
9104  format (' problem type: ',a35)
9200  format (a1,3i2,a35)
9201  format ('Branch ',i2,'.   sol.type: ',a20)
9300  format (a1,i2,i4,2i2,a1,a2,18e20.13)
9301  format ('Solution No.',i2,'.',i4)
9302  format (' sol.vector:',5e13.5)
9303  format (' parameters:',5e13.5)
9304  format (' period:    ',5e13.5)
9606  format (e20.10)
      end

