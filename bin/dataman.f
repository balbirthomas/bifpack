
c    BIFPACK -- a package for Bifurcation, Continuation  and Stability Analysis.
c    Copyright (C) 1983, 1999   R. Seydel

c    This program is free software; you can redistribute it and/or modify
c    it under the terms of the GNU General Public License as published by
c    the Free Software Foundation; either version 2 of the License, or
c    (at your option) any later version.

c    This program is distributed in the hope that it will be useful,
c    but WITHOUT ANY WARRANTY; without even the implied warranty of
c    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c    GNU General Public License for more details.

c    You should have received a copy of the GNU General Public License
c    along with this program; if not, write to the Free Software
c    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

c=========================================================================

      PROGRAM MAIN
c     Program   dataman.f   , version 15. Aug 1998
c     "data manager" for conversion and distribution of various data types,
c     and for the print output of solution files.

c  options:  see list of write -statement below

c  FORMAT of 9300 as well as dimensions NDIML :
c  to be adapted if necessary!

c  required software:
c     ODE-integrator (here RKF.; it may be necessary to use stiff integrator)

      parameter (NDIML=14)
      DOUBLE PRECISION Y1,P,PERIOD,ymin,ymax,zetaco,etaco,parco,
     1 t0,t1,t2,tdel,tend,hinit,hmax,eps,y2,y3,y,t,window,ystart
      character*300 buffer
      character*35 name1,text1
      character*1 MODEL,label
      character*2 judged
      Dimension y1(ndiml),p(ndiml),ymin(ndiml),ymax(ndiml)
     1   ,y(ndiml),y2(ndiml),y3(ndiml),t(7),window(2),ystart(ndiml)
      COMMON /CONTRO/ ETACO,ZETACO,PARCO, NCO,N1CO,N2CO,
     1                KCO,IND2CO,NSECFLAG,NFLAG
      external fcnper
       OPEN (3,FILE ='DATBASE')
       OPEN (9,FILE= 'PRINTLIST')
       OPEN (11,FILE='DATBRANCH')
       OPEN (12,FILE='DATONE')
       OPEN (13,FILE='DATIME')
       OPEN (14,FILE='DATBASENEW')
       OPEN (15,FILE='DATFIELD')
       OPEN (16,FILE='START')
      istrict=1
      call setup (model,n,m,t,name1,jacobi,eps,iflgsol,window,ystart)
      nscreen=0
      ymin(1)=99999.
      ymax(1)=-99999.

      write (*,*) 'This is the DATA MANAGER, working on DATBASE.'
11    write (*,*) 'SEARCHING  data:  enter'
      write (*,*) '1  for searching a single SOLUTION,'
      write (*,*) '2  for searching an entire BRANCH,'
      write (*,*) 'EDITING the file:  enter'
      write (*,*) '3  for BROWSING the entire data base,'
      write (*,*) '4  for RENUMBERING within each branch,'
      write (*,*) '5  for entering text records,'
      write (*,*) '6  for renumbering branches,'
      write (*,*) 'working on a single solution:   enter'
      write (*,*) '7  for calculating TIME DEPENDENCE,'
      write (*,*) '8  for setting up a TANGENT FIELD,'
      write (*,*) '9  for SIMULATION,'
      write (*,*) '-1 for terminating.'
      write (*,*) 'ENTER:'
      read (*,*) nsearch
      if (nsearch.eq.-1) stop
      if (nsearch.lt.1.or.nsearch.gt.9) go to 11
      if (nsearch.eq.5) then
	 call texts
	 go to 9
	 endif
      if (nsearch.eq.6) then
	 call branchnu
	 go to 9
	 endif
      if (nsearch.eq.7 .or. nsearch.eq.8
     1     .or.nsearch.eq.9) go to 2
      if (nsearch.eq.1) then
	   istrict=0
	 else
	  write (*,*) 'output of ONLY SOLUTION (label  R ) records:'
     1    ,'  enter  1'
	  write (*,*) 'output of "ALL" DATA (labels R and r) enter'
     1     , '  0           ENTER:'
	  read (*,*) istrict
	 endif
      if (nsearch.eq.4) go to 1
      if (nsearch.eq.3) go to 15
      if (nsearch.le.2) then
	 write (*,*) 'enter no. of branch:'
	 read (*,*) nbranch
	 endif
      if (nsearch.eq.1) then
	 write (*,*) 'enter No. of solution on that branch:'
	 read (*,*) nsolu
	 nscreen=1
	 go to 1
	 endif
15    write (*,*) '1 for output also on screen, 0 for no screen:'
      read (*,*) nscreen
1     call browse (nsearch,nscreen,nbranch,nsolu,nfound,istrict)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (nfound.ne.1.or.nsearch.ne.1) go to 9
	 rewind 12
2        read (12,9000,end=25) buffer
	 write (*,*) 'Data from file DATONE are read.'
         write (*,*) 'Hope this file is up to date.' 
	 if (nsearch.eq.8) go to 53
c         if (buffer(12:12).eq.'0') go to 9
	 if (nsearch.ne.7.and.nsearch.ne.9) then
	   write (*,*) ' '
	   write (*,*) 'Enter  1  if you want to create a file of'
	   write (*,*) 'the TIME DEPENDENCE; else enter  0 :'
	   read (*,*) nwish2
	   if (nwish2.eq.0) go to 5
	   endif
	 if (buffer(13:13).ne.'P'.and.buffer(13:13).ne.'H') 
     1      read (buffer,9300) label,
     2      nob,nos,n,npar,MODEL,judged,(y1(i),i=1,n),(p(i),i=1,npar)
	 if (buffer(13:13).eq.'P'.or.buffer(13:13).eq.'H') then
	    read (buffer,9300) label,nob,nos,n,npar,MODEL,judged,
     1           (y1(i),i=1,n),(p(i),i=1,npar),period
	    tend=period
	    endif
	 if (nsearch.eq.9) go to 26
	 if (buffer(13:13).ne.'P') then
	     if (MODEL.eq.'d') then
		 tinit=0.
		 write (*,*) 'Enter FINAL TIME:'
		 read (*,*) tend
		 endif
	     if (MODEL.eq.'b') then
		 write (*,*) 'Enter the initial  t :'
		 read (*,*) tinit
		 write (*,*) 'Enter end of interval:'
		 read (*,*) tend
		 endif
	     endif
      go to 27

25    write (*,*) 'You have no data in DATONE.'
c      write (*,*) 'Enter N:'
c      read (*,*) n
26    write (*,*) 'Starting vector:'
      write (*,*) '  Enter  1  if you want to continue previous'
      write (*,*) '  simulation, for new data enter  0 :'
      read (*,*) nsimul
      if (nsimul.eq.1) then
	    rewind 16
	    do 262 i=1,n+1
262         read (16,*,end=263) y1(i)
	    go to 264
	    endif
263   do 261 i=1,n
      write (*,*) 'Enter y component #',i
261   read (*,*) y1(i)
      write (*,*) 'Enter parameter:'
      read (*,*) y1(n+1)
264   p(1)=y1(n+1)
      write (*,*) 'Enter final time for simulation:'
      read (*,*) tend



27    t0=tinit
      number=1001
      write (*,*) 'Calculating the TIME-DEPENDENCE:'
      write (*,*) '"Standard" means integration of ONE period,'
      write (*,*) 'calculating  1001  points for output.'
      write (*,*) 'If you accept standard, enter 1 , else enter 0:'
      read (*,*) nwisht
      if (nwisht.eq.1) go to 3
      write (*,*) 'current final time is:',tend
      write (*,*) 'Enter  1  for accepting this,'
      write (*,*) 'enter  2  for doubling,'
      write (*,*) 'enter  3  for tripling and so on... ,'
      write (*,*) 'or enter  0  for entering a new final time:'
      read (*,*) nwisht
      if (nwisht.gt.1) tend=tend*nwisht
      if (nwisht.le.0.or.tend.le.t0) then
	  write (*,*) 'now enter new final time:'
	  read (*,*) tend
	  endif
      write (*,*) 'Enter number of points:'
      read (*,*) number
3     eps=1.e-5
      n1=n+1
      y1(n1)=p(1)
      do 31 j=1,n1
      ymin(j)=9999999.
      ymax(j)=-9999999.
      y2(j)=y1(j)
      y3(j)=y1(j)
31    y(j)=y1(j)
      if (TEND.le.T0) then
          write (*,*) ' final time is not set.'
          go to 27
          endif
      TDEL=(TEND-T0)/FLOAT(NUMBER-1)
      HINIT=0.01*TDEL
      T1=T0
      nper=2
      perguess=0.
      HMAX=TDEL
      zetaco=0.
      ind2co=0
      n1co=n1
      y1(n1co)=p(1)
      kco=1
      etaco=y1(kco)
      WRITE (13,9001) y1(n+1),name1
      WRITE (13,9400) nob,nos,n,number,1,t1,(Y1(j),j=1,n)

      DO 4 I=1,NUMBER-1
      T2=T1+TDEL
      CALL RKF7 (N1,FCNPER,T1,Y1,T2,EPS,HMAX,HINIT)
      if (Hinit.lt.(tend-t0)*1.d-10) then
          write (*,*) 'stepsize too small. h=',hinit
          write (*,*) 'final value of  t  is: ',t1
          write (*,*) 'You may require a stiff integrator.'
          stop
          endif
      if (i.gt.1 .and. nsearch.eq.9 .and. perguess.lt.eps) then
c         rough guess of possible period, based on NPER components:
c         (y2 is initial vector, y3 the previous solution point)
	  iflag =1
	  do 42 j=1,nper
	  if (y1(j).le.y3(j) .and. y2(j).lt.y1(j)) iflag=0
	  if (y1(j).le.y3(j) .and. y2(j).gt.y3(j)) iflag=0
	  if (y1(j).ge.y3(j) .and. y2(j).lt.y3(j)) iflag=0
	  if (y1(j).ge.y3(j) .and. y2(j).gt.y1(j)) iflag=0
42        continue
	  if (iflag.eq.1) perguess=t1
	  endif
      do 41 j=1,n
      y3(j)=y1(j)
      if (y1(j).gt.ymax(j)) ymax(j)=y1(j)
41    if (y1(j).lt.ymin(j)) ymin(j)=y1(j)
      WRITE (13,9400) nob,nos,n,number,i+1,t1,(Y1(j),j=1,n)
4     CONTINUE
      write (*,*) 'Minima and maxima (in low precision); par:',y1(n+1)
      write (*,*) '              y1          y2         ...'
      write (*,9306) 'minimum: ',(ymin(j),j=1,n)
      write (*,9306) 'maximum: ',(ymax(j),j=1,n)
      write (*,9306) 'amplitude',(ymax(j)-ymin(j),j=1,n)
9306  format (a10,5e12.5)
      if (nsearch.eq.9) then
	   rewind 16
	   do 43 j=1,n+1
43         write (16,*) y1(j)
	   endif
      if (perguess.gt.eps) write (*,*)
     1    'Possible period may be close to ',perguess
      write (*,*) 'Time dependence written to file DATIME.'
      if (MODEL.eq.'2') go to 9
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
5     write (*,*) ' '
      write (*,*) 'Next Option: TANGENT FIELD:'
      write (*,*) 'Enter 1 to set up a tangent field, else enter 0:'
      read (*,*) nwisht
      if (nwisht.ne.1) go to 9
      go to 55
53    read (buffer,9300) label,
     1      nob,nos,n,npar,MODEL,judged,(y1(i),i=1,n),(p(i),i=1,npar)
      if (MODEL.eq.'b') then
	write (*,*) '..  is not available for B.V.problems.'
	go to 9
	endif
55    continue
      if (n.gt.2) write (*,*) 'WARNING: tangent field may not make',
     1            ' sense for dimensions > 2  ************'
      WRITE (15,9001) y1(n+1),name1
      call field (y,nob,nos,ymin,ymax)
      write (*,*) 'Tangent field written to file DATFIELD.'
9     stop
9000  format (a300)
9001  format ('par=',e13.7,' ',a35)
9300  format (a1,i2,i4,2i2,a1,a2,18e20.13)
9400  format ('T',i2,i4,i2,2i4,14e20.13)
      end



      SUBROUTINE FIELD (Y,NOB,NOS,YMIN,YMAX)
      parameter (NDIML=14)
      double precision t0,
     2   xleft,xright,yleft,yright,y(ndiml),f(ndiml),param(ndiml)
     3   ,xdelta,ydelta,facmax,window(4)
     4   ,ymin(ndiml),ymax(ndiml)
c    The following may be changed:
c    number of points NPOINTS per axis:
      write (*,*) '"Standard" means to choose indices 1 and 2, ...'
      write (*,*) 'Enter  1  for STANDARD treatment, else  0 :'
      read (*,*) nwish
      if (nwish.eq.1) then
	 ind1=1
	 ind2=2
	 npoints=21
	 endif
      if (nwish.eq.0) then
	 write (*,*) 'Now enter number of points per direction :'
	 read (*,*) npoints
	 write (*,*) 'Enter index for horizontal axis:'
	 read (*,*) ind1
	 write (*,*) 'Enter index for vertical axis:'
	 read (*,*) ind2
	 endif
      if (ymin(1).le.ymax(1)) then
	 window(1)=ymin(ind1)-0.1*dmax1(1.0d0,abs(ymin(ind1)))
	 window(2)=ymax(ind1)+0.1*dmax1(1.0d0,abs(ymax(ind1)))
	 window(3)=ymin(ind2)-0.1*dmax1(1.0d0,abs(ymin(ind2)))
	 window(4)=ymax(ind2)+0.1*dmax1(1.0d0,abs(ymax(ind2)))
	    write (*,*) 'suggested window is:'
	    write (*,*) 'horizontal:',(window(j),j=1,2)
	    write (*,*) 'vertical:',(window(j),j=3,4)
	    xleft=window(1)
	    xright=window(2)
	    yleft=window(3)
	    yright=window(4)
	    write (*,*) 'Window O.K.?  Enter 1 for accepting,'
	    write (*,*) '              enter 0 for new window:'
	    read (*,*) nwish
	    if (nwish.eq.1) go to 3
	 endif

2     write (*,*) 'Enter WINDOW: horizontal axis: enter left end:'
      read (*,*) xleft
      write (*,*) 'horizontal axis: enter right end:'
      read (*,*) xright
      write (*,*) 'vertical axis: enter lower end:'
      read (*,*) yleft
      write (*,*) 'vertical axis: enter upper end:'
      read (*,*) yright
      if (xleft.ge.xright .or. yleft.ge.yright) then
	 write (*,*) 'not monotonous!  re-enter!'
	 go to 2
	 endif
3     xdelta=(xright-xleft)/(npoints-1)
      ydelta=(yright-yleft)/(npoints-1)
      write (15,9500) 'W',nob,nos,xleft,xright,yleft,yright

      facmax=99999999999.
      y(IND2)=yleft-ydelta
      do 5 j=1,npoints
      y(IND2)=y(IND2)+ydelta
      y(IND1)=xleft-xdelta
      do 5 i=1,npoints
      y(IND1)=y(IND1)+xdelta
      call fcn (t0,y,f)
      if (facmax.gt.abs(xdelta/f(ind1)))
     1    facmax=abs(xdelta/f(ind1))
      if (facmax.gt.abs(ydelta/f(ind2)))
     1    facmax=abs(ydelta/f(ind2))
5     continue

      y(IND2)=yleft-ydelta
      do 55 j=1,npoints
      y(IND2)=y(IND2)+ydelta
      y(IND1)=xleft-xdelta
      do 55 i=1,npoints
      y(IND1)=y(IND1)+xdelta
      call fcn (t0,y,f)
      write (15,9500) 'F',nob,nos,
     1       y(IND1),y(IND2),y(IND1)+f(IND1)*facmax,
     2       y(IND2)+facmax*f(IND2)
55    continue
      return
9500  format (a1,i2,i4,4e20.13)
      end


      SUBROUTINE BROWSE (nsearch,nscreen,nbranch,nsolu,nfound
     1                    ,istrict)
      DOUBLE PRECISION Y1,P,PERIOD
      character*300 buffer
      character*35 name1,text1
      character*1 MODEL,label
      character*2 judged
      integer lastb,ncountin,ncountout
      parameter (NDIML=14)
      dimension y1(ndiml),p(ndiml)

c  istrict=1 :  only solution records with label  R  are written to output.
c  istrict=0 :  also data with label  r  are written to output.
      nfound=0
      lastb=-1
      ncountin=0
      ncountout=0
      ncountsec=0
      do 2 loop=1,10000
      read (3,9000,end=9) buffer

      if (buffer(1:1).eq.'P') then
	 read (buffer,9100,end=2) label,n,npar,MODEL,name1
	 if (nsearch.eq.4) then
	     write (14,9100) label,n,npar,MODEL,name1 
	     go to 2
	     endif
	 if (nscreen.eq.1) then
	     write (*,9001)
	     write (*,9101) name1
	     write (*,9102) n,npar
	     if (MODEL.eq.'a') write (*,9104) 'algebraic eq.'
	     if (MODEL.eq.'d') write (*,9104) 'ODE, autonomous'
	     if (MODEL.eq.'b') then
		 write (*,9104) 'ODE, boundary value problem'
		 endif
	     write (*,*) ' '
	     endif
	 write (9,9001)
	 write (9,9101) name1
	 write (9,9102) n,npar
	 if (MODEL.eq.'a') write (9,9104) 'algebraic eq.'
	 if (MODEL.eq.'d') write (9,9104) 'ODE, autonomous'
	 if (MODEL.eq.'b') then
		 write (9,9104) 'ODE, boundary value problem'
		 endif
	 write (9,*) ' '
	 endif

      if (buffer(1:1).eq.'B') then
	 read (buffer,9200,end=2) label,nob,n,npar,text1
	 if (nsearch.eq.4) then
	     write (14,9200) label,nob,n,npar,text1
	     go to 2
	     endif
	 if (nsearch.eq.1) go to 2
	 if (nsearch.eq.2.and.nob.ne.nbranch) go to 2
	 if (nscreen.eq.1) then
	     write (*,9001)
	     write (*,9201) nob
	     write (*,9103) text1
	     write (*,*) ' '
	     endif
	 write (9,9001)
	 write (9,9201) nob
	 write (9,9103) text1
	 write (9,*) ' '
	 endif

      if (buffer(1:1).eq.'R'.or.buffer(1:1).eq.'r') then
	 ncountin=ncountin+1
	 iprint=0
	 if (buffer(1:1).eq.'R'.or. istrict.eq.0) iprint=1
	 iperiod=0
	 if (buffer(13:13).eq.'P'.or.
     1   buffer(13:13).eq.'H') iperiod=1
	 if (iperiod.eq.0) read (buffer,9300) label,
     1      nob,nos,n,npar,MODEL,judged,(y1(i),i=1,n),(p(i),i=1,npar)
	 if (iperiod.eq.1)
     1      read (buffer,9300) label,nob,nos,n,npar,MODEL,judged,
     2           (y1(i),i=1,n),(p(i),i=1,npar),period
	 if (ncountin.eq.1) then
	    if ((n+npar+1)*20+14.gt.300) then
	       write (*,*) 'Format for buffer too short. STOP'
	       stop
	       endif
	    endif
	 if (nsearch.eq.4) then
	   if (nob.ne.lastb) then
	       lastb=nob
	       nosnew=0
	       endif
	   nosnew=nosnew+1
	   if (iprint.eq.0) ncountsec=ncountsec+1
	   if (iprint.eq.1) ncountout=ncountout+1
	   if (iperiod.eq.1.and.iprint.eq.1)
     1           write (14,9300) label,nob,nosnew,n,npar,
     2           MODEL,judged,(y1(i),i=1,n),(p(i),i=1,npar),period
	   if (iperiod.ne.1.and.iprint.eq.1)
     1           write (14,9300) label,nob,nosnew,n,npar,
     2           MODEL,judged,(y1(i),i=1,n),(p(i),i=1,npar)
	   go to 2
	   endif
	 if (nsearch.eq.2.and.nob.ne.nbranch) go to 2
	 if (nsearch.eq.1) then
	    if (nob.ne.nbranch) go to 2
	    if (nos.ne.nsolu) go to 2
	    endif
	 if (nsearch.le.2) nfound=1
	 if (iprint.eq.1) ncountout=ncountout+1
	 if (iprint.eq.0) ncountsec=ncountsec+1
	 if (nsearch.eq.1) then
	    if (iperiod.eq.0.and.iprint.eq.1) then
	       write (12,9300) label,nob,nos,n,npar,
     2         MODEL,judged,(y1(i),i=1,n),(p(i),i=1,npar)
	       write (16,9606) (y1(i),i=1,n),p(1)
	       endif
	    if (iperiod.eq.1.and.iprint.eq.1) then
	       write (12,9300) label,nob,nos,n,npar,
     2         MODEL,judged,(y1(i),i=1,n),(p(i),i=1,npar),period
	       write (16,9606) (y1(i),i=1,n),p(1),period
	       endif
	    endif
	 if (nsearch.eq.2) then
	    if (iperiod.eq.0.and.iprint.eq.1)
     1         write (11,9300) label,nob,nos,n,npar,
     2         MODEL,judged,(y1(i),i=1,n),(p(i),i=1,npar)
	    if (iperiod.eq.1.and.iprint.eq.1)
     1         write (11,9300) label,nob,nos,n,npar,
     2         MODEL,judged,(y1(i),i=1,n),(p(i),i=1,npar),period
	    endif
	 if (nscreen.eq.1.or.nsearch.eq.1) then
	     write (*,9301) nob,nos
	     if (npar.gt.0) write (*,9303) (p(i),i=1,npar)
	     write (*,9302) (y1(i),i=1,n)
	     if (MODEL.eq.'1') write (*,9304) period
	     endif
	 write (9,9301) nob,nos
	 if (npar.gt.0) write (9,9303) (p(i),i=1,npar)
	 if (buffer(1:1).eq.'R') write (9,9302) (y1(i),i=1,n)
	 if (buffer(1:1).eq.'r') write (9,9305) (y1(i),i=1,n)
	 if (MODEL.eq.'1') write (9,9304) period
	 endif

2     continue
9     continue
      write (*,9001)
      write (*,*) ncountin,' solutions in input,'
      write (*,*) ncountout,' solutions in output.'
      if (ncountsec.gt.0) write (*,*) ncountsec, ' interpolated',
     1   ' approximations with label  r  are not in output.'
      write (9,9001)
      write (9,*) 'The file DATBASE includes in total ',ncountin,
     1            ' solutions.'
      if (nsearch.eq.4) then
	 write (*,*) 'renumbered data on file DATBASENEW'
	 return
	 endif
      if (nfound.eq.0.and.nsearch.le.2)
     1   write (*,*) 'not found.'
      if (nfound.eq.1.and.nsearch.eq.2)
     1   write (*,*) 'Branch No.',nbranch,' written to DATBRANCH'
      if (nfound.eq.1.and.nsearch.eq.1)
     1   write (*,*) 'Solution written to DATONE and START'
      write (*,*) 'Print-output written to PRINTLIST'
      return
9000  format (a300)
9001  format (60('='))
9100  format (a1,2i2,a1,a35,2e20.13)
9101  format ('Problem : ',a35)
9102  format (' Dimensions: ',i2,' eqs. and ',i2,' parameters')
9103  format ('   ',a35)
9104  format (' problem type: ',a35)
9200  format (a1,3i2,a35)
9201  format ('Branch ',i2)
9300  format (a1,i2,i4,2i2,a1,a2,18e20.13)
9301  format ('Solution No.   ',i2,'.',i4)
9302  format (' sol.vector:   ',4e15.7)
9303  format (' parameters:   ',4e15.7)
9304  format (' period:       ',4e15.7)
9305  format (' approximation:',4e15.7)
9606  format (e20.10)
      end

      SUBROUTINE TEXTS
      DOUBLE PRECISION Y1,P,PERIOD,T,WINDOW,YSTART,eps
      character*300 buffer
      character*35 name1,text1
      character*1 MODEL,label,last
      character*2 judged
      integer lastb,ncountin,ncountout
      parameter (NDIML=14)
      dimension y1(ndiml),p(ndiml),t(7),window(2),ystart(ndiml)
      last=' '
      lastb=-1
      ncountin=0
      ncountout=0

      do 2 loop=1,10000
      miss=0
      read (3,9000,end=9) buffer
      if (loop.eq.1 .and. buffer(1:1).ne.'P') miss=1
      if (miss.eq.1) then
         call setup (model,n,m,t,name1,jacobi,eps,
     1               iflagsol,window,ystart)
c	 write (*,*) 'Enter dimension of state vector:'
c	 read (*,*) n
	 write (*,*) 'Enter dimension of parameter vector (mostly 1):'
	 read (*,*) npar
c	 write (*,*) 'Enter type of equation:'
c	 write (*,*) '  a for algebraic eq.,'
c	 write (*,*) '  d for autonous ODE,'
c	 write (*,*) '  b for ODE boundary val.problem. ENTER:'
c	 read (*,9100) MODEL
c	 write (*,*) 'Enter a brief text to characterize problem:'
c	 read (*,9105) name1
	 if (MODEL.ne.'b') write (14,9100) 'P',n,npar,MODEL,name1
	 if (MODEL.eq.'b') then
c	    write (*,*) 'Enter left end of interval:'
c	    read (*,*) tinit
c	    write (*,*) 'Enter right end of interval:'
c	    read (*,*) tfinal
            tinit=t(1)
            tfinal=t(m)
	    write (14,9100) 'P',n,npar,MODEL,name1,tinit,tfinal
	    endif
	 endif

      if (buffer(1:1).eq.'P') then
	 read (buffer,9100,end=2) label,n,npar,MODEL,name1
	 write (14,9100) label,n,npar,MODEL,name1
	 go to 2
	 endif

      if (buffer(1:1).eq.'B') then
	 read (buffer,9200,end=2) label,nob,n,npar,text1
	 write (14,9200) label,nob,n,npar,text1
	 endif

      if (buffer(1:1).eq.'-') write (14,9000) buffer

      if (buffer(1:1).eq.'R'.or.buffer(1:1).eq.'r') then
	 ncountin=ncountin+1
	 iperiod=0
	 if (buffer(13:13).eq.'P'.or.buffer(13:13).eq.'H') iperiod=1
	 if (iperiod.eq.0) read (buffer,9300) label,
     1     nob,nos,n,npar,MODEL,judged,(y1(i),i=1,n),(p(i),i=1,npar)
	 if (iperiod.eq.1)
     1      read (buffer,9300) label,nob,nos,n,npar,MODEL,judged,
     2           (y1(i),i=1,n),(p(i),i=1,npar),period
	 if (ncountin.eq.1) then
	    if ((n+npar+1)*20+14.gt.300) then
	       write (*,*) 'Format for buffer too short. STOP'
	       stop
	       endif
	    endif
	 if (nob.ne.lastb) then
	      if (last.ne.'B') then
		 write (*,*) 'Enter text to characterize branch',nob
		 read (*,9105) text1
		 write (14,9200) 'B',nob,n,npar,text1
		 endif
	      lastb=nob
	      endif
	 ncountout=ncountout+1
	 if (iperiod.eq.1) write (14,9300) label,nob,nos,n,npar,
     1        MODEL,judged,(y1(i),i=1,n),(p(i),i=1,npar),period
	 if (iperiod.ne.1) write (14,9300) label,nob,nos,n,npar,
     1        MODEL,judged,(y1(i),i=1,n),(p(i),i=1,npar)
	 go to 2
	 endif

2     last=label
9     continue
      write (*,9001)
      write (*,*) ncountin,' solutions in input,'
      write (*,*) ncountout,' solutions in output, plus text records.'
      write (*,*) 'Output written to  DATBASENEW .'
      return
9000  format (a300)
9001  format (60('='))
9100  format (a1,2i2,a1,a35,2e20.13)
9105  format (a35)
9200  format (a1,3i2,a35)
9300  format (a1,i2,i4,2i2,a1,a2,18e20.13)
      end

      SUBROUTINE BRANCHNU
      DOUBLE PRECISION Y1,P,PERIOD
      character*300 buffer
      character*35 name1,text1
      character*1 MODEL,label,last
      character*2 judged
      integer lastb,ncountin,ncountout
      parameter (NDIML=14)
      dimension y1(ndiml),p(ndiml)
      last=' '
      lastb=-1
      ncountin=0
      ncountout=0

      do 2 loop=1,10000
      read (3,9000,end=9) buffer

      if (buffer(1:1).eq.'P') then
	 read (buffer,9100,end=2) label,n,npar,MODEL,name1
	 write (14,9100) label,n,npar,MODEL,name1
	 go to 2
	 endif

      if (buffer(1:1).eq.'B') then
	 read (buffer,9200,end=2) label,nob,n,npar,text1
	 write (14,9200) label,nob,n,npar,text1
	 endif

      if (buffer(1:1).eq.'-') write (14,9000) buffer

      if (buffer(1:1).eq.'R'.or.buffer(1:1).eq.'r') then
	 ncountin=ncountin+1
	 iperiod=0
	 if (buffer(13:13).eq.'P'.or.buffer(13:13).eq.'H') iperiod=1
	 if (iperiod.eq.0) read (buffer,9300) label,
     1      nob,nos,n,npar,MODEL,judged,(y1(i),i=1,n),(p(i),i=1,npar)
	 if (iperiod.eq.1)
     1      read (buffer,9300) label,nob,nos,n,npar,MODEL,judged,
     2           (y1(i),i=1,n),(p(i),i=1,npar),period
	 if (ncountin.eq.1) then
	    if ((n+npar+1)*20+14.gt.300) then
	       write (*,*) 'Format for buffer too short. STOP'
	       stop
	       endif
	    endif
	 if (nob.ne.lastb) then
	      write (*,*) 'next branch has number:',nob
	      write (*,*) 'Enter new number for that branch:'
	      read (*,*) nobnew
	      lastb=nob
	      endif
	 ncountout=ncountout+1
	 if (iperiod.eq.1) write (14,9300) label,nobnew,nos,n,npar,
     1           MODEL,judged,(y1(i),i=1,n),(p(i),i=1,npar),period
	 if (iperiod.ne.1) write (14,9300) label,nobnew,nos,n,npar,
     1           MODEL,judged,(y1(i),i=1,n),(p(i),i=1,npar)
	 go to 2
	 endif

2     last=label
9     continue
      write (*,9001)
      write (*,*) ncountin,' solutions in input,'
      write (*,*) ncountout,' solutions in output, plus text records.'
      write (*,*) 'Output written to  DATBASENEW .'
      return
9000  format (a300)
9001  format (60('='))
9100  format (a1,2i2,a1,a35,2e20.13)
9105  format (a35)
9200  format (a1,3i2,a35)
9300  format (a1,i2,i4,2i2,a1,a2,18e20.13)
      end

      SUBROUTINE FCNPER (T,Y,F)
      DOUBLE PRECISION T,Y,F,ETACO,ZETACO,PARCO
      COMMON /CONTRO/ ETACO,ZETACO,PARCO, NCO,N1CO,N2CO,
     1                KCO,IND2CO,NSECFLAG,NFLAG
      PARAMETER (NDIM=12,NDIML=14)
      DIMENSION Y(NDIML),F(NDIML)
C  This is more or less a dummy routine with the only purpose to
C  take care of  F(N+1)=0
      CALL FCN (T,Y,F)
      F(N1CO)=0.
      RETURN
      END

