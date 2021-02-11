
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

c========================================================================

      PROGRAM MAIN
c     Program   datextrB.f   , version 26.Feb.99
c  input:  DATBASENEW,   output: DATARRAY
C                               (may be used as input for GNUplot)

      parameter (NDIML=14)
      DOUBLE PRECISION Y1,P,PERIOD,zetaco,etaco,parco,
     1 t0,t1,t2,t3,tdel,tend,hinit,hmax,eps,t,window,ystart
      character*300 buffer
      character*35 name1
      character*1 MODEL,label
      character*2 judged
      Dimension y1(ndiml),p(ndiml),t(7),window(2),ystart(ndiml)
      COMMON /CONTRO/ ETACO,ZETACO,PARCO, NCO,N1CO,N2CO,
     1                KCO,IND2CO,NSECFLAG,NFLAG
       OPEN (11,FILE='DATBASENEW')
       OPEN (16,FILE='DATARRAY')
      index=0

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
1     read (11,9000,end=8) buffer
      if (buffer(1:1).eq.'P'.or.buffer(1:1).eq.'B') go to 1
      if (buffer(13:13).ne.'P'.and.buffer(13:13).ne.'H') 
     1      read (buffer,9300) label,
     2      nob,nos,n,npar,MODEL,judged,(y1(i),i=1,n),(p(i),i=1,npar)
      if (buffer(13:13).eq.'P'.or.buffer(13:13).eq.'H') 
     1    read (buffer,9300) label,nob,nos,n,npar,MODEL,judged,
     1           (y1(i),i=1,n),(p(i),i=1,npar),period
c      if (label.ne.'R') go to 1
      if (index.eq.0) then
         write (*,*) 'input: data of file DATBASENEW'
         write (*,*) 'You may like to copy your DATBASE into DATBASENEW'
         write (*,*)'Dimension ist n=',n
         write (*,*)'Output will be pairs:  parameter and y(k)'
         write (*,*)'Enter index k of component y(k):'
         read (*,*) index
         if (index.lt.1.or. index.gt.n) then
            write (*,*) 'We choose k=1'
            index=1
            endif
         endif
      write (16,*) p(1),y1(index)
      go to 1

8     write (*,*) 'output written to file DATARRAY'
      write (*,*) 'In case of several branches,'
      write (*,*) '   you may insert empty lines.'
      stop
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
9000  format (a300)
9300  format (a1,i2,i4,2i2,a1,a2,18e20.13)
9601  format ()
9606  format (3e20.10)
      end
