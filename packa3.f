ccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccc
cccc
cccc     this is packA3  of BIFPACK
cccc
cccc     STEADY-STATE METHOD FOR CALCUL. HOPF POINTS
ccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccc


c    BIFPACK -- a package for Bifurcation, Continuation, and Stability Analysis.
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

c=======================================================================

c   packa3,  version 17. Jan 1996


      SUBROUTINE HOPFA (N,Y,PAR,EVAIM,EVERE,EVEIM,PERIOD,HFLAG,EPS
     1                  ,UNITLST)
      DOUBLE PRECISION Y,EVERE,EVEIM,Z,PAR,PERIOD,EPS,PI2,
     1     EVAIM,ETACO,ZETACO,PARCO
      PARAMETER (NDIM=12,NDIML=14)
      DIMENSION Y(NDIML),EVERE(NDIML),EVEIM(NDIML),Z(NDIML)
      INTEGER HFLAG,UNITLST

C  THIS SUBROUTINE CALCULATES HOPF BIFURCATION POINTS BY A STEADY-STATE METHOD
C  EXPLOITING THE EIGENVALUE PROPERTIES IN THE HOPF POINT.
C  REFERENCES: A.Jepson: Phd thesis 1981, D.Roose, several papers
C              A.Griewank, G.Reddien: IMA J.NUMER.ANAL. 3 (1983) 295
C              M.Holodniok, M.Kubicek: APPL.MATH.COMPUT. 15 (1984) 261-274
C--------------------------------------------------------------------------
C  THIS ROUTINE DOES NOT CALCULATE SMALL AMPLITUDE OSCILLATIONS.
C  FOR THE CALCULATION OF THE DYNAMICS USE FOR INSTANCE THE ODE-ROUTINE
C  H O P F  RATHER THAN THE STEADY-STATE ROUTINE  H O P F A .
C  This routine is not optimized for storage.
C--------------------------------------------------------------------------

C  NEEDED ROUTINES (FYHOPF IS INCLUDED):
C  ---------------
C       NOLE   A NONLINEAR EQUATION SOLVER (NEWTON'S METHOD)
C              CALLING SEQUENCE SEE IN ROUTINE  S W I T A
C       FCN    SUBROUTINE DEFINING THE FUNCTION  F(Y,PARAMETER)
C              CALLING SEQUENCE:  FCN(DUMMY,Y,F)
C              THE PARAMETER IS GIVEN BY Y(N+1).
C       FCNLIN   A SUBROUTINE CALCULATING THE N*N PARTIAL DERIVATIVES
C              A(I,J)=DF(I)/DY(J)
C              CALLING SEQUENCE:  FCNLIN(DUMMY,Y,A)
C              THE PARAMETER IS GIVEN BY Y(N+1).

C  INPUT DATA:
C  ----------
C       N      DIMENSION= NUMBER OF COMPONENTS OF FUNCTION F
C       Y      VECTOR OF STEADY-STATE CLOSE TO HOPF POINT
C       PAR    CORRESPONDING PARAMETER VALUE
C       EVAIM  GUESS OF IMAGINARY PART OF EIGENVALUE
C       EVERE  GUESS OF REAL PART OF EIGENVECTOR
C       EVEIM  GUESS OF IMAGINARY PART OF EIGENVECTOR
C       HFLAG  =1: OUTPUT IS PRINTED , =0: NO PRINTING

C  OUTPUT DATA:
C  -----------
C       HFLAG  =1 in case of success,  =0 in case of no success,
C              =-1 in case of a zero eigenvalue (turning point)
C       PERIOD : INITIAL PERIOD OF EMANATING PERIODIC SOL. BRANCH
C       Y, EVAIM, EVERE, EVEIM  ARE UPDATED IN CASE OF CONVERGENCE

C  HOW TO USE OUTPUT DATA:
C  ----------------------
C       FOR SMALL VALUES OF ABS(DELTA) AN APPROXIMATION TO A SMALL-AMPLITUDE
C       OSCILLATION  Y(T)  FOR  0.LE.T.LE.1.  IS GIVEN BY
C       Y(T) = Y + (EVERE * COS(2*PI*T) - EVEIM * SIN(2*PI*T))*DELTA
C-------------------------------------------------------------------------

      EXTERNAL FYHOPF
      COMMON /CONTRO/ ETACO,ZETACO,PARCO, NCO,N1CO,N2CO,
     1                KCO,IND2CO,NSECFLAG,NFLAG
      WRITE (*,*) 'BIFPACK: now entering HOPFA:'
      WRITE (UNITLST,*) 'BIFPACK: now entering HOPFA:'
      NCO=N
      PI2=6.28318530717959d0
      N1=N+2
      DO 2 I=1,N
      Z(I)=Y(I)
      Z(I+N1)=EVERE(I)
2     Z(I+N1+N)=EVEIM(I)
      Z(N+1)=PAR
      Z(N1)=EVAIM
      N3=3*N+2
      IF (N3.GT.NDIML) THEN
          WRITE (*,*) 'Dimensions in HOPFA to small.'
          RETURN
          ENDIF 
      ITMAX=20
      J=-1
      CALL NOLE (N3,Z,FYHOPF,EPS,ITMAX,J)
      IF (J.LT.0) THEN
             WRITE (UNITLST,*) ' NO SUCCESS OF  H O P F A'
             HFLAG=0 
             RETURN
             ENDIF
      IF (ABS(Z(N+2)).LT.EPS) THEN
          WRITE (UNITLST,*) ' Convergence to a zero eigenvalue:'
          WRITE (UNITLST,*) ' may be a turning point, '
     1                      ,'but no Hopf point.'
          HFLAG=-1
          RETURN
          ENDIF  
      HFLAG=1
      PERIOD=PI2/Z(N+2)
      DO 4 I=1,N
      Y(I)=Z(I)
      EVERE(I)=Z(I+N1)
4     EVEIM(I)=Z(I+N+N1)
      PAR=Z(N+1)
      EVAIM=Z(N1)
      WRITE (UNITLST,601) J,PAR,PERIOD
      WRITE (*,601) J,PAR,PERIOD
      WRITE (UNITLST,602)
      DO 5 I=1,N
5     WRITE (UNITLST,604) Y(I),EVERE(I),EVEIM(I)
      WRITE (*,604) (Y(I),I=1,N)
      RETURN

601   FORMAT (' HOPF point calculated (within',I3,' iterations)',
     1       /'    PARAMETER:',E20.11,' ,    PERIOD:',E20.11,10X) 
602   FORMAT (/8X,'solution          eigenvector       eigenvector',
     1        /8X,'components        real part         imaginary p.')
604   FORMAT ('   ',4E19.11)
      END


      SUBROUTINE FYHOPF (DUMMY,Y,F)  
      DOUBLE PRECISION Y,A,F,S1,S2,DUMMY,ETACO,ZETACO,PARCO
      PARAMETER (NDIM=12,NDIML=14)
      DIMENSION Y(NDIML),A(NDIML,NDIML),F(NDIML)
C    DIMENSIONS AT LEAST:  Y(3*N+2),F(3*N+2),A(N,N)
      COMMON /CONTRO/ ETACO,ZETACO,PARCO, N,N1CO,N2CO,
     1                KCO,IND2CO,NSECFLAG,NFLAG
      CALL FCN (DUMMY,Y,F)
      CALL FCNLIN (DUMMY,Y,A)
C--  STORAGE:
C    POSITION 1 THROUGH N         CONTAINS THE PRESENT STEADY STATE Y
C    POSITION N+1                 CONTAINS THE PARAMETER
C    POSITION N+2                 CONTAINS THE COMPLEX PART OF EIGENVALUE
C    POSITION N+3 THROUGH 2N+2    CONTAINS THE REALPART OF EIGENVECTOR
C    POSITION 2N+3 THROUGH 3N+2   CONTAINS THE IMAGINARY PART

      DO 5 I=1,N
      S1=0.d0
      S2=0.d0
      DO 4 J=1,N
      S1=S1+A(I,J)*Y(N+2+J)
4     S2=S2+A(I,J)*Y(N+N+2+J)
      F(N+I)=S1+Y(N+2)*Y(N+N+2+I)
5     F(N+N+I)=S2-Y(N+2)*Y(N+2+I)
      F(3*N+1)=Y(N+3)-1.
      F(3*N+2)=Y(N+N+3)
      RETURN
      END
