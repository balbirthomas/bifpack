
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
c     Program   datextrT.f   , version 18 Jun 1999
c     "data extracting" to set up simple files for graphics.
c     Version for phase portraits / time-dependence.
c     output file is named GINPUT, input file is DATIME  .
c     FORMAT 9500 for output file GINPUT.

c  FORMAT of 9300 as well as dimensions NDIML currently dimension  14
c  to be adapted if necessary!

      parameter (NDIML=14)
      DOUBLE PRECISION Y1(NDIML),T
c     1        ,P(NDIML),PERIOD
      character*300 buffer
c      character*35 name1,text1
      character*53 text2
      character*1 label
c     1           ,model
c      character*2 judged
       OPEN (3, FILE='DATIME')
       OPEN (13,FILE='GINPUT')
      write (*,*) 'This is routine  datextrT  working on file DATIME.'
      index=2
      text2=' '
      istart=0      
      jold=-1
      write (*,*) 'Enter  1  for PHASE portrait,'
      write (*,*) 'enter  0  for y(t) -dependence:'
      read (*,*) iwish

      do 2 loop=1,10000
      read (3,9000,end=9) buffer

      if (buffer(1:1).eq.'p') text2=buffer

      if (buffer(1:1).eq.'T') then
         read (buffer,9400) label,nob,nos,n,m,j,t,(y1(i),i=1,n)
         if (istart.eq.0) then
            istart=1

            if (iwish.eq.1) then
             if (n.eq.2) then
               index1=1
               index2=2
               endif
             if (n.gt.2) then
               write (*,*) 'Phase space is ',n,'-dimensional.'
c               write (*,*) ' y(1) will be projected on horizontal axis.'
               write (*,*) 'Choose index for component of  y  to be',
     1                    ' projected to HORIZONTAL AXIS:'
               read (*,*) index1
               write (*,*) 'Choose index for component of  y  to be',
     1                    ' projected to VERTICAL AXIS:'
               read (*,*) index2
               if (index1.lt.1) index1=1
               if (index2.gt.n) index1=n
               endif
             write (13,9502) text2,index2,index1
             endif 

            if (iwish.ne.1) then
             write (*,*) 'Which component of vector  Y  is to be shown?'
             write (*,*) 'Enter index (.le.',n,') :'
             read (*,*) index
             if (index.lt.1.or.index.gt.n) index=1
             write (13,9501) text2,index
            endif                                    
           endif 

         if (j.lt.jold) then                         
            write (13,9500) 'B',0.,0.
            endif
         jold=j
         If (iwish.eq.1) write (13,9500) y1(index1),y1(index2)
         if (iwish.ne.1) write (13,9500) t,y1(index)
         endif
2     continue      
9     write (*,*) 'File for plot is GINPUT; call gnuplot'
      write (*,*) 'Remove first line first!'
      stop
9000  format (a300)
9001  format (a76)
9100  format (a1,2i2,a1,a35,2e20.13)
9101  format ('Problem : ',a35)
9102  format (' Dimensions: ',i2,' eqs. and ',i2,' parameters')
9103  format ('   ',a35)
9104  format (' problem type: ',a35)
9200  format (a1,3i2,a35)
9201  format ('Branch ',i2)
9300  format (a1,i2,i4,2i2,a1,a2,18e20.13)
9301  format ('Solution No.',i2,'.',i4)
9302  format (' sol.vector:',5e13.5)
9303  format (' parameters:',5e13.5)
9304  format (' period:    ',5e13.5)
9400  format (a1,i2,i4,i2,2i4,14e20.13) 
9500  format (2e15.7)
9501  format (a53,' time dependence of y',i2)
9502  format (a50,' phase diagram: y',i2,' vs. y',i2)
9606  format (e20.10)
      end

