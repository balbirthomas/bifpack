      PROGRAM MAIN
c     Program   summary.f   , version 23.2.1995
c     "data manager"  to summarize bifurcation results

       OPEN (3,FILE ='DATBASE')
       OPEN (9,FILE= 'PRINTLIST')
      call overview
      stop
      end


      SUBROUTINE OVERVIEW
      DOUBLE PRECISION Y1,P,PERIOD
      character*300 buffer
      character*35 name1,text1
      character*1 kezi1,label
      character*2 kezi2
c      integer lastb,ncountin,ncountout
      parameter (NDIML=14)
      dimension y1(ndiml),p(ndiml)
      last=0
      last1=0
      last2=0
      do 2 loop=1,10000
      read (3,9000,end=9) buffer

      if (buffer(1:1).eq.'P') then
	 read (buffer,9100,end=2) label,n,m,kezi1,name1
	     write (*,9001)
	     write (*,9101) name1
	     write (*,9102) n,m
	     if (kezi1.eq.'a') write (*,9104) 'algebraic eq.'
	     if (kezi1.eq.'d') write (*,9104) 'ODE, autonomous'
	     if (kezi1.eq.'b') then
		 write (*,9104) 'ODE, boundary value problem'
		 endif
	     write (*,9001)
	 write (9,9001)
	 write (9,9101) name1
	 write (9,9102) n,m
	 if (kezi1.eq.'a') write (9,9104) 'algebraic eq.'
	 if (kezi1.eq.'d') write (9,9104) 'ODE, autonomous'
	 if (kezi1.eq.'b') then
		 write (9,9104) 'ODE, boundary value problem'
		 endif
	 write (9,9001)
	 endif

      if (buffer(1:1).eq.'B') then
	 read (buffer,9200,end=2) label,nob,n,m,text1
	 write (9,*) ' '
	 write (*,*) ' '
	 write (*,9201) nob,text1
	 write (9,9201) nob,text1
	 endif

      if (buffer(1:1).eq.'R'.or.buffer(1:1).eq.'r') then
	 iperiod=0
	 if (buffer(13:13).eq.'P'.or.
     1       buffer(13:13).eq.'H') iperiod=1
	 if (iperiod.eq.0) read (buffer,9300) label,
     1      nob,nos,n,m,kezi1,kezi2,(y1(i),i=1,n),(p(i),i=1,m)
	 if (iperiod.eq.1)
     1      read (buffer,9300) label,nob,nos,n,m,kezi1,kezi2,
     2           (y1(i),i=1,n),(p(i),i=1,m),period
	 last2=last1
	 last1=last
	 if (kezi2.eq.'S ') last=1
	 if (kezi2.eq.'U ') last=-1
	 if (kezi2.ne.'S '.and.kezi2.ne.'U ') last=0
	 if (last.eq.1.and.last2.eq.-1) write (*,*)
     1 'Stability is gained (relative to the numbering of solutions).'
	 if (last.eq.-1.and.last2.eq.1) write (*,*)
     1 'Stability is lost (relative to the numbering of solutions).'
	 if (kezi2.eq.'T ') write (*,9002) p(1),'turning point'
	 if (kezi2.eq.'HB'.and.iperiod.eq.1)
     1       write (*,9003) p(1),'Hopf bifurcation',period
	 if (kezi2.eq.'HB'.and.iperiod.eq.0)
     1       write (*,9002) p(1),'Hopf bifurcation'
	 if (kezi2.eq.'PT'.and.iperiod.eq.1)
     1       write (*,9003) p(1),'torus bifurcation',period
	 if (kezi2.eq.'PD'.and.iperiod.eq.1)
     1       write (*,9003) p(1),'period doubling',period

	 if (last.eq.1.and.last2.eq.-1) write (9,*)
     1       'Stability is gained.'
	 if (last.eq.-1.and.last2.eq.1) write (9,*)
     1       'Stability is lost.'
	 if (kezi2.eq.'T ') write (9,9002) p(1),'turning point'
	 if (kezi2.eq.'HB'.and.iperiod.eq.1)
     1       write (9,9003) p(1),'Hopf bifurcation',period
	 if (kezi2.eq.'HB'.and.iperiod.eq.0)
     1       write (9,9002) p(1),'Hopf bifurcation'
	 if (kezi2.eq.'PT'.and.iperiod.eq.1)
     1       write (9,9003) p(1),'torus bifurcation',period
	 if (kezi2.eq.'PD'.and.iperiod.eq.1)
     1       write (9,9003) p(1),'period doubling',period
	 endif

2     continue
9     continue
      write (*,*) ' '
      write (*,*) '(This is also written to PRINTLIST.)'
      return
9000  format (a300)
9001  format (60('='))
9002  format (' at parameter ',e11.5,a20)
9003  format (' at parameter ',e11.5,a20,'. Period:',e11.5)
9100  format (a1,2i2,a1,a35,2e20.13)
9101  format (' Problem : ',a35)
9102  format (' Dimensions: ',i2,' eqs. and ',i2,' parameters')
9103  format ('   ',a35)
9104  format (' problem type: ',a35)
9200  format (a1,3i2,a35)
9201  format (' Branch ',i2,'.   sol.type: ',a20)
9300  format (a1,i2,i4,2i2,a1,a2,13e20.13)
9301  format (' Solution No.',i2,'.',i4)
9302  format (' sol.vector:',5e13.5)
9303  format (' parameters:',5e13.5)
9304  format (' period:    ',5e13.5)
9606  format (e20.10)
      end

