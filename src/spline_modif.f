C-*-fortran-*-
c
c
c   VERSION 2.2
c
c   This library contains routines for B-spline interpolation in
c   one, two, and three dimensions. Part of the routines are based
c   on the book by Carl de Boor: A practical guide to Splines (Springer,
c   New-York 1978) and have the same calling sequence and names as
c   the corresponding routines from the IMSL library. For documen-
c   tation see the additional files. NOTE: The results in the demo
c   routines may vary slightly on different architectures.
c
c   by W. Schadow 07/19/98
c   last changed by W. Schadow 07/28/2000
c
c
c   Wolfgang Schadow
c   TRIUMF
c   4004 Wesbrook Mall
c   Vancouver, B.C. V6T 2A3
c   Canada
c
c   email: schadow@triumf.ca  or  schadow@physik.uni-bonn.de
c
c   www  : http://www.triumf.ca/people/schadow
c
c
c  ------------------------------------------------------------------
c
c
c   Copyright (C) 2000 Wolfgang Schadow
c
c   This library is free software; you can redistribute it and/or
c   modify it under the terms of the GNU Library General Public
c   License as published by the Free Software Foundation; either
c   version 2 of the License, or (at your option) any later version.
c
c   This library is distributed in the hope that it will be useful,
c   but WITHOUT ANY WARRANTY; without even the implied warranty of
c   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
c   Library General Public License for more details.
c
c   You should have received a copy of the GNU Library General Public
c   License along with this library; if not, write to the
c   Free Software Foundation, Inc., 59 Temple Place - Suite 330,
c   Boston, MA  02111-1307, USA.
c
c
c  ------------------------------------------------------------------
c
c
c   The following routines are included:
c
c            dbsnak
c
c            dbsint
c            dbsval
c            dbsder
c            dbs1gd
c
c            dbs2in
c            dbs2dr
c            dbs2vl
c            dbs2gd
c
c            dbs3in
c            dbs3vl
c            dbs3dr
c            dbs3gd
c
c  ------------------------------------------------------------------
c
c  NEW: corrected some error messages
c       some changes in the checks of dbs3dg to find a possible
c       error earlier.
c
c  ------------------------------------------------------------------
c
c  NEW: documentation included, changed some comments
c
c  ------------------------------------------------------------------
c

c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      SUBROUTINE DBSNAK(NX,XVEC,KXORD,XKNOT)

c
c  Compute the `not-a-knot' spline knot sequence.
c  (see de Boor p. 167)
c
c   nx     - number of data points.  (input)
c   xvec   - array of length ndata containing the location of the
c            data points.  (input)
c   kxord  - order of the spline.  (input)
c   xknot  - array of length ndata+korder containing the knot
c            sequence.  (output)
c

      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)

      DIMENSION XVEC(NX),XKNOT(NX+KXORD)

      LOGICAL FIRST

      SAVE FIRST,EPS

      DATA FIRST/.TRUE./

      IF (FIRST) THEN
         FIRST=.FALSE.
         EPS = DLAMCH('PRECISION')
C         WRITE(6,*) 'SUBROUTINE DBSNAK: '
C         WRITE(6,*) 'EPS = ',EPS
      ENDIF

      IF((KXORD .LT. 0) .OR. (KXORD .GT. NX)) THEN
         WRITE(6,*) 'SUBROUTINE DBSNAK: ERROR'
         WRITE(6,*) '0 <= KXORD <= NX IS REQUIRED.'
         WRITE(6,*) 'KXORD = ', KXORD, ' AND NX = ', NX,  ' IS GIVEN.'
         STOP
      ENDIF

      DO 30 I = 1, KXORD
         XKNOT(I) = XVEC(1)
 30   CONTINUE

      IF(MOD(KXORD,2) .EQ. 0) THEN
         DO 40 I = KXORD+1,NX
            XKNOT(I) = XVEC(I-KXORD/2)
 40      CONTINUE
      ELSE
         DO 50 I = KXORD+1,NX
            XKNOT(I) = 0.5D0 * (XVEC(I-KXORD/2) + XVEC(I-KXORD/2-1))
 50      CONTINUE
      ENDIF

      DO 60 I = NX+1,NX+KXORD
         XKNOT(I) = XVEC(NX) * (1.0D0 + EPS)
 60   CONTINUE

      RETURN

      END


c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      SUBROUTINE DBSINT(NX, XVEC, XDATA, KX, XKNOT, BCOEF)

c
c  Computes the spline interpolant, returning the B-spline coefficients.
c  (see de Boor p. 204)
c
c   nx     - number of data points.  (input)
c   xvec   - array of length nx containing the data point
c            abscissas.  (input)
c   xdata  - array of length ndata containing the data point
c            ordinates.  (input)
c   kx     - order of the spline.  (input)
c            korder must be less than or equal to ndata.
c   xknot  - array of length nx+kx containing the knot
c            sequence.  (input)
c            xknot must be nondecreasing.
c   bscoef - array of length ndata containing the B-spline
c            coefficients.  (output)
c

      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)

      PARAMETER(KXMAX=8, NXMAX=4096)

      DIMENSION XDATA(NX),XVEC(NX),XKNOT(NX+KX),BCOEF(NX)
      DIMENSION WORK((2*KXMAX-1)*NXMAX)

      IF ((KX .GT. KXMAX) .OR. (NX .GT. NXMAX)) THEN
         WRITE(6,*) 'SUBROUTINE DBSINT: ERROR'
         WRITE(6,*) 'KX > KXMAX OR NX > NXMAX'
         WRITE(6,*) 'KX = ',KX,'  KXMAX = ',KXMAX
         WRITE(6,*) 'NX = ',NX,'  NXMAX = ',NXMAX
         STOP
      ENDIF

      NXP1  = NX + 1
      KXM1  = KX - 1
      KPKM2 = 2 * KXM1
      LEFTX = KX
      LENQ  = NX * (KX + KXM1)

      DO 10 I = 1, LENQ
         WORK(I) = 0.D0
 10   CONTINUE

      DO 20 IX = 1,NX
         XVECI  = XVEC(IX)
         ILP1MX = MIN0(IX+KX,NXP1)
         LEFTX   = MAX0(LEFTX,IX)
         IF (XVECI .LT. XKNOT(LEFTX)) GOTO 998
 30      IF (XVECI .LT. XKNOT(LEFTX+1)) GO TO 40
         LEFTX = LEFTX + 1
         IF (LEFTX .LT. ILP1MX) GO TO 30
         LEFTX = LEFTX - 1
         IF (XVECI .GT. XKNOT(LEFTX+1)) GOTO 998
 40      CALL BSPLVB (XKNOT,NX+KX,KX,1,XVECI,LEFTX,BCOEF)
         JJ = IX - LEFTX + 1 + (LEFTX - KX) * (KX + KXM1)
         DO 50 IK = 1,KX
            JJ       = JJ + KPKM2
            WORK(JJ) = BCOEF(IK)
 50      CONTINUE
 20   CONTINUE

      CALL BANFAC(WORK,KX+KXM1,NX,KXM1,KXM1,IFLAG)
      GO TO (60,999), IFLAG

 60   DO 70 IX = 1,NX
         BCOEF(IX) = XDATA(IX)
 70   CONTINUE

      CALL BANSLV(WORK,KX+KXM1,NX,KXM1,KXM1,BCOEF)

      RETURN

 998  WRITE(6,*) 'SUBROUTINE DBSINT:'
      WRITE(6,*) 'XKNOT(IX) <= XKNOT(IX+1) REQUIRED.'
      WRITE(6,*) IX,XKNOT(IX),XKNOT(IX+1)

      STOP

 999  WRITE(6,*) 'SUBROUTINE DBSINT: ERROR'
      WRITE(6,*) 'NO SOLUTION OF LINEAR EQUATION SYSTEM !!!'

      STOP

      END


c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      FUNCTION DBSVAL(X,KX,XKNOT,NX,BCOEF)

c
c  Evaluates a spline, given its B-spline representation.
c
c   x      - point at which the spline is to be evaluated.  (input)
c   kx     - order of the spline.  (input)
c   xknot  - array of length nx+kx containing the knot
c            sequence.  (input)
c            xknot must be nondecreasing.
c   nx     - number of B-spline coefficients.  (input)
c   bcoef  - array of length nx containing the B-spline
c            coefficients.  (input)
c   dbsval - value of the spline at x.  (output)
c

      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)

      PARAMETER(KXMAX=8)

      DIMENSION XKNOT(NX+KX),BCOEF(NX)
      DIMENSION WORK(KXMAX),DL(KXMAX),DR(KXMAX)


c
c     check if kx <= kxmax
c

      IF(KX .GT. KXMAX) THEN
         WRITE(6,*) 'SUBROUTINE DBSVAL:'
         WRITE(6,*) 'KX <= KXMAX REQUIRED.'
         STOP
      ENDIF

c
c     check if xknot(i) <= xknot(i+1) and calculation of i so that
c     xknot(i) <= x < xknot(i+1)
c

      LEFTX = 0

      DO 10 IX = 1,NX+KX-1
         IF (XKNOT(IX) .GT. XKNOT(IX+1)) THEN
             WRITE(6,*) 'SUBROUTINE DBSVAL:'
             WRITE(6,*) 'XKNOT(IX) <= XKNOT(IX+1) REQUIRED.'
             WRITE(6,*) IX,XKNOT(IX),XKNOT(IX+1)
          STOP
          ENDIF
         IF((XKNOT(IX) .LE. X) .AND. (X .LT. XKNOT(IX+1))) LEFTX = IX
10    CONTINUE

      IF(LEFTX .EQ. 0) THEN
         WRITE(6,*) 'SUBROUTINE DBSVAL:'
         WRITE(6,*) 'IX WITH XKNOT(IX) <= X < XKNOT(IX+1) REQUIRED.'
         WRITE(6,*) 'X = ', X
         STOP
      ENDIF

      DO 20 IK = 1,KX-1
         WORK(IK) = BCOEF(LEFTX+IK-KX)
         DL(IK)   = X - XKNOT(LEFTX+IK-KX)
         DR(IK)   = XKNOT(LEFTX+IK) - X
 20   CONTINUE

      WORK(KX)  = BCOEF(LEFTX)
      DL(KX)    = X - XKNOT(LEFTX)

      DO 30 IK = 1,KX-1
         SAVE2 = WORK(IK)
         DO 40  IL = IK+1,KX
            SAVE1 = WORK(IL)
            WORK(IL) = (DL(IL) * WORK(IL) + DR(IL-IK) * SAVE2)
     .           / (DL(IL) + DR(IL - IK))
            SAVE2 = SAVE1
 40      CONTINUE
 30   CONTINUE

      DBSVAL = WORK(KX)

      RETURN

      END


c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      FUNCTION DBSDER(IDERX,X,KX,XKNOT,NX,BCOEF)

c
c  Evaluates the derivative of a spline, given its B-spline representation.
c
c
c   iderx  - order of the derivative to be evaluated.  (input)
c            in particular, iderx = 0 returns the value of the
c            spline.
c   x      - point at which the spline is to be evaluated.  (input)
c   kx     - order of the spline.  (input)
c   xknot  - array of length nx+kx containing the knot
c            sequence.  (input)
c            xknot must be nondecreasing.
c   nx     - number of B-spline coefficients.  (input)
c   bcoef  - array of length nx containing the B-spline
c            coefficients.  (input)
c   dbsder - value of the iderx-th derivative of the spline at x.
c            (output)
c

      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)

      PARAMETER(KXMAX=8)

      DIMENSION XKNOT(NX+KX),BCOEF(NX)
      DIMENSION WORK(KXMAX),DL(KXMAX),DR(KXMAX),BSP(KXMAX)

c
c     check if <= kxmax
c

      IF(KX .GT. KXMAX) THEN
         WRITE(6,*) 'SUBROUTINE DBSDER:'
         WRITE(6,*) 'KX <= KXMAX REQUIRED.'
         STOP
      ENDIF

c
c     check if xknot(i) <= xknot(i+1) and calculation of i so that
c     xknot(i) <= x < xknot(i+1)
c

      LEFTX = 0
      DO 10 IX = 1,NX+KX-1
         IF (XKNOT(IX) .GT. XKNOT(IX+1)) THEN
            WRITE(6,*) 'SUBROUTINE DBSDER:'
            WRITE(6,*) 'XKNOT(IX) <= XKNOT(IX+1) REQUIRED.'
            STOP
         ENDIF
         IF ((XKNOT(IX) .LE. X) .AND. (X .LT. XKNOT(IX+1))) LEFTX = IX
 10   CONTINUE

      IF (LEFTX .EQ. 0) THEN
         WRITE(6,*) 'SUBROUTINE DBSDER:'
         WRITE(6,*) 'IX WITH XKNOT(IX) <= X < XKNOT(IX+1) REQUIRED.'
         WRITE(6,*) 'XKNOT(1)     = ', XKNOT(1)
         WRITE(6,*) 'XKNOT(NX+KX) = ', XKNOT(NX+KX)
         WRITE(6,*) '         X   = ', X
         STOP
      ENDIF

      IF (IDERX .EQ. 0) THEN

         DO 20 IK = 1,KX-1
            WORK(IK) = BCOEF(LEFTX+IK-KX)
            DL(IK)   = X - XKNOT(LEFTX+IK-KX)
            DR(IK)   = XKNOT(LEFTX+IK) - X
 20      CONTINUE

         WORK(KX)  = BCOEF(LEFTX)
         DL(KX)    = X - XKNOT(LEFTX)

         DO 30 IK = 1,KX-1
            SAVE2 = WORK(IK)
            DO 40  IL = IK+1,KX
               SAVE1 = WORK(IL)
               WORK(IL) = (DL(IL) * WORK(IL) + DR(IL-IK) * SAVE2)
     .              / (DL(IL) + DR(IL - IK))
               SAVE2 = SAVE1
 40         CONTINUE
 30      CONTINUE

         DBSDER = WORK(KX)

      ELSEIF ((IDERX .GE. 1) .AND. (IDERX .LT. KX)) THEN
         BSP(1) = 1.0D0
         DO 50 IK = 1,KX-IDERX-1
            DR(IK) = XKNOT(LEFTX+IK) - X
            DL(IK) = X - XKNOT(LEFTX+1-IK)
            SAVE   = BSP(1)
            BSP(1) = 0.0D0
            DO 60 IL = 1,IK
               Y         = SAVE / (DR(IL) + DL(IK+1-IL))
               BSP(IL)   = BSP(IL) + DR(IL) * Y
               SAVE      = BSP(IL+1)
               BSP(IL+1) = DL(IK+1-IL) * Y
 60         CONTINUE
 50      CONTINUE

         DO 70 IK = 1,KX
            WORK(IK) = BCOEF(LEFTX+IK-KX)
            DR(IK)   = XKNOT(LEFTX+IK) - X
            DL(IK)   = X - XKNOT(LEFTX+IK-KX)
 70      CONTINUE

         DO 80 IK = 1,IDERX
            DIK   = DBLE(KX - IK)
            SAVE2 = WORK(IK)
            DO 90  IL = IK+1,KX
               SAVE1    = WORK(IL)
               WORK(IL) = DIK * (WORK(IL) - SAVE2) /(DL(IL) + DR(IL-IK))
               SAVE2    = SAVE1
 90         CONTINUE
 80      CONTINUE

         SUM = 0.0D0

         DO 100 I = 1,KX-IDERX
            SUM = SUM + BSP(I) * WORK(IDERX+I)
 100     CONTINUE

         DBSDER = SUM

      ELSE
         DBSDER = 0.0D0
      ENDIF

      RETURN

      END


c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      SUBROUTINE DBS1GD(IDERX,NXVEC,XVEC,KX,XKNOT,NX,BCOEF,VAL)

c
c  Evaluates the derivative of a spline on a grid, given its B-spline
c  representation.
c
c   iderx  - order of the derivative to be evaluated.  (input)
c            in particular, iderx = 0 returns the value of the
c            spline.
c   nxvec  - length of vector xvec.  (input)
c   xvec   - array of length nxvec containing the points at which the
c            spline is to be evaluated.  (input)
c            xvec should be strictly increasing.
c   kx     - order of the spline.  (input)
c   xknot  - array of length nx+kx containing the knot
c            sequence.  (input)
c            xknot must be nondecreasing.
c   nx     - number of B-spline coefficients.  (input)
c   bcoef  - array of length nx containing the B-spline
c            coefficients.  (input)
c   val    - array of length nxvec containing the values of the
c            iderx-th derivative of the spline at the points in
c            xvec.  (output)
c

      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)

      PARAMETER(KXMAX=8,NXMAX=350)

      DIMENSION XVEC(NXVEC),XKNOT(NX+KX),BCOEF(NX),VAL(NXVEC)
      DIMENSION DL(NXMAX,KXMAX),DR(NXMAX,KXMAX),SAVE1(NXMAX)
      DIMENSION BIATX(NXMAX,KXMAX),TERM(NXMAX),SAVE2(NXMAX)
      DIMENSION LEFTX(NXMAX),WORK(NXMAX,KXMAX)

      LOGICAL SAME,NEXT

c
c     check if kx <= kxmax
c

      IF(KX .GT. KXMAX) THEN
         WRITE(6,*) 'SUBROUTINE DBS1GD:'
         WRITE(6,*) 'KX <= KXMAX REQUIRED.'
         STOP
      ENDIF

      LEFTX(1) = 0

      CALL HUNTN(XKNOT,NX+KX,KX,XVEC(1),LEFTX(1))

      DO 10 IX = 2,NXVEC
         LEFTX(IX) = LEFTX(IX-1)
         SAME = (XKNOT(LEFTX(IX)) .LE. XVEC(IX))
     .        .AND. (XVEC(IX) .LE. XKNOT(LEFTX(IX)+1))
         IF(.NOT. SAME ) THEN
            LEFTX(IX) = LEFTX(IX) + 1
            NEXT      = (XKNOT(LEFTX(IX)) .LE. XVEC(IX))
     .           .AND. (XVEC(IX) .LE. XKNOT(LEFTX(IX)+1))
            IF (.NOT. NEXT)
     .           CALL HUNTN(XKNOT,NX+KX,KX,XVEC(IX),LEFTX(IX))
         ENDIF
 10   CONTINUE

      DO 20 I = 1,NX+KX-1
         IF (XKNOT(I) .GT. XKNOT(I+1)) THEN
            WRITE(6,*) 'SUBROUTINE DBS1GD:'
            WRITE(6,*) 'XKNOT(I) <= XKNOT(I+1) REQUIRED.'
            WRITE(6,*) I, XKNOT(I), XKNOT(I+1)
            WRITE(6,*)
            WRITE(6,*) XKNOT
            STOP
         ENDIF
 20   CONTINUE

      DO 30 I = 1,NXVEC
         IF ((XVEC(I).LT.XKNOT(1)).OR.(XVEC(I).GT.XKNOT(NX+KX))) THEN
            WRITE(6,*) 'SUBROUTINE DBS1GD:'
            WRITE(6,*) 'IX WITH XKNOT(IX) <= X < XKNOT(IX+1) REQUIRED.'
            WRITE(6,*) 'X = ', XVEC(I)
            STOP
         ENDIF
 30   CONTINUE

      IF (IDERX .EQ. 0) THEN

         DO 40 IX = 1,NXVEC
            BIATX(IX,1) = 1.D0
            VAL(IX)     = 0.D0
 40      CONTINUE

         DO 50 IK = 1,KX-1
            DO 60 IX = 1,NXVEC
               DR(IX,IK) = XKNOT(LEFTX(IX)+IK) - XVEC(IX)
               DL(IX,IK) = XVEC(IX) - XKNOT(LEFTX(IX)+1-IK)
               SAVE1(IX) = 0.D0
 60         CONTINUE

            DO 70 IL = 1,IK
               DO 80 IX = 1,NXVEC
                  TERM(IX)     = BIATX(IX,IL)
     .                 / (DR(IX,IL) + DL(IX,IK+1-IL))
                  BIATX(IX,IL) = SAVE1(IX) + DR(IX,IL) * TERM(IX)
                  SAVE1(IX)    = DL(IX,IK+1-IL) * TERM(IX)
 80            CONTINUE
 70         CONTINUE

            DO 90 IX = 1,NXVEC
               BIATX(IX,IK+1) = SAVE1(IX)
 90         CONTINUE
 50      CONTINUE

         DO 100 IK = 1,KX
            DO 110 IX = 1,NXVEC
               VAL(IX) = VAL(IX) + BIATX(IX,IK) * BCOEF(LEFTX(IX)-KX+IK)
 110        CONTINUE
 100     CONTINUE

      ELSEIF ((IDERX .GE. 1) .AND. (IDERX .LT. KX)) THEN

         DO 120 IX = 1,NXVEC
            BIATX(IX,1) = 1.D0
            VAL(IX)     = 0.D0
 120     CONTINUE

         DO 130 IK = 1,KX-IDERX-1
            DO 140 IX = 1,NXVEC
               DR(IX,IK)   = XKNOT(LEFTX(IX)+IK) - XVEC(IX)
               DL(IX,IK)   = XVEC(IX) - XKNOT(LEFTX(IX)+1-IK)
               SAVE1(IX)    = BIATX(IX,1)
               BIATX(IX,1) = 0.0D0
               DO 150 IL = 1,IK
                  TERM(IX)       = SAVE1(IX)
     .                 / (DR(IX,IL) + DL(IX,IK+1-IL))
                  BIATX(IX,IL)   = BIATX(IX,IL) + DR(IX,IL) * TERM(IX)
                  SAVE1(IX)      = BIATX(IX,IL+1)
                  BIATX(IX,IL+1) = DL(IX,IK+1-IL) * TERM(IX)
 150           CONTINUE
 140        CONTINUE
 130     CONTINUE

         DO 160 IK = 1,KX
            DO 170 IX = 1,NXVEC
               WORK(IX,IK) = BCOEF(LEFTX(IX)+IK-KX)
               DR(IX,IK)   = XKNOT(LEFTX(IX)+IK) - XVEC(IX)
               DL(IX,IK)   = XVEC(IX) - XKNOT(LEFTX(IX)+IK-KX)
 170        CONTINUE
 160     CONTINUE

         DO 180 IK = 1,IDERX
            DIK   = DBLE(KX - IK)
            DO 190 IX = 1,NXVEC
               SAVE2(IX) = WORK(IX,IK)
               DO 200  IL = IK+1,KX
                  SAVE1(IX)   = WORK(IX,IL)
                  WORK(IX,IL) = DIK * (WORK(IX,IL) - SAVE2(IX))
     .                 /(DL(IX,IL) + DR(IX,IL-IK))
                  SAVE2(IX)   = SAVE1(IX)
 200           CONTINUE
 190        CONTINUE
 180     CONTINUE

         DO 210 I = 1,KX-IDERX
            DO 220 IX = 1,NXVEC
               VAL(IX) = VAL(IX) + BIATX(IX,I) * WORK(IX,IDERX+I)
 220        CONTINUE
 210     CONTINUE

      ELSE

         DO 230 IX = 1,NXVEC
            VAL(IX) = 0.0D0
 230     CONTINUE

      ENDIF

      RETURN

      END


c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      FUNCTION DBSDCA(IDERX,X,KX,XKNOT,NX,BCOEF,LEFTX)

c
c This routine is equivalent to the routine dbsder, but it does not
c check the parameters!!c
c
c Evaluates the derivative of a spline, given its B-spline representation.
c
c
c   iderx  - order of the derivative to be evaluated.  (input)
c            in particular, iderx = 0 returns the value of the
c            spline.
c   x      - point at which the spline is to be evaluated.  (input)
c   kx     - order of the spline.  (input)
c   xknot  - array of length nx+kx containing the knot
c            sequence.  (input)
c            xknot must be nondecreasing.
c   nx     - number of B-spline coefficients.  (input)
c   bcoef  - array of length nx containing the B-spline
c            coefficients.  (input)
c   leftx  - number of the intervall of xknot that includes x
c   dbsdca - value of the ideriv-th derivative of the spline at x.
c            (output)
c

      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)

      PARAMETER(KXMAX=8)

      DIMENSION XKNOT(NX+KX),BCOEF(NX)
      DIMENSION WORK(KXMAX),DL(KXMAX),DR(KXMAX),BSP(KXMAX)


      IF (IDERX .EQ. 0) THEN

         DO 20 IK = 1,KX-1
            WORK(IK) = BCOEF(LEFTX+IK-KX)
            DL(IK)   = X - XKNOT(LEFTX+IK-KX)
            DR(IK)   = XKNOT(LEFTX+IK) - X
 20      CONTINUE

         WORK(KX)  = BCOEF(LEFTX)
         DL(KX)    = X - XKNOT(LEFTX)

         DO 30 IK = 1,KX-1
            SAVE2 = WORK(IK)
            DO 40  IL = IK+1,KX
               SAVE1 = WORK(IL)
               WORK(IL) = (DL(IL) * WORK(IL) + DR(IL-IK) * SAVE2)
     .              / (DL(IL) + DR(IL - IK))
               SAVE2 = SAVE1
 40         CONTINUE
 30      CONTINUE

         DBSDCA = WORK(KX)

      ELSEIF ((IDERX .GE. 1) .AND. (IDERX .LT. KX)) THEN
         BSP(1) = 1.0D0
         DO 50 IK = 1,KX-IDERX-1
            DR(IK) = XKNOT(LEFTX+IK) - X
            DL(IK) = X - XKNOT(LEFTX+1-IK)
            SAVE   = BSP(1)
            BSP(1) = 0.0D0
            DO 60 IL = 1,IK
               Y         = SAVE / (DR(IL) + DL(IK+1-IL))
               BSP(IL)   = BSP(IL) + DR(IL) * Y
               SAVE      = BSP(IL+1)
               BSP(IL+1) = DL(IK+1-IL) * Y
 60         CONTINUE
 50      CONTINUE

         DO 70 IK = 1,KX
            WORK(IK) = BCOEF(LEFTX+IK-KX)
            DR(IK)   = XKNOT(LEFTX+IK) - X
            DL(IK)   = X - XKNOT(LEFTX+IK-KX)
 70      CONTINUE

         DO 80 IK = 1,IDERX
            DIK   = DBLE(KX - IK)
            SAVE2 = WORK(IK)
            DO 90  IL = IK+1,KX
               SAVE1    = WORK(IL)
               WORK(IL) = DIK * (WORK(IL) - SAVE2) /(DL(IL) + DR(IL-IK))
               SAVE2    = SAVE1
 90         CONTINUE
 80      CONTINUE

         SUM = 0.0D0

         DO 100 I = 1,KX-IDERX
            SUM = SUM + BSP(I) * WORK(IDERX+I)
 100     CONTINUE

         DBSDCA = SUM

      ELSE
         DBSDCA = 0.0D0
      ENDIF

      RETURN

      END


c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      SUBROUTINE DBS2IN(NX,XVEC,NY,YVEC,XYDATA,LDF,KX,KY,XKNOT,
     .     YKNOT,BCOEF)

c
c  Computes a two-dimensional tensor-product spline interpolant,
c  returning the tensor-product B-spline coefficients.
c
c    nx     - number of data points in the x-direction.  (input)
c    xvec   - array of length nx containing the data points in
c             the x-direction.  (input)
c             xdata must be strictly increasing.
c    ny     - number of data points in the y-direction.  (input)
c    yvec   - array of length ny containing the data points in
c             the y-direction.  (input)
c             ydata must be strictly increasing.
c    xydata - array of size nx by ny containing the values to
c             be interpolated.  (input)
c             fdata(i,j) is the value at (xdata(i),ydata(j)).
c    ldf    - the leading dimension of fdata exactly as specified in
c             the dimension statement of the calling program.
c             (input)
c    kx     - order of the spline in the x-direction.  (input)
c             kxord must be less than or equal to nxdata.
c    ky     - order of the spline in the y-direction.  (input)
c             kyord must be less than or equal to nydata.
c    xknot  - array of length nx+kx containing the knot
c             sequence in the x-direction.  (input)
c             xknot must be nondecreasing.
c    yknot  - array of length ny+ky containing the knot
c             sequence in the y-direction.  (input)
c             yknot must be nondecreasing.
c    bcoef  - array of length nx*ny containing the
c             tensor-product B-spline coefficients.  (output)
c             bscoef is treated internally as a matrix of size nxdata
c             by nydata.
c

      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)

      PARAMETER(NXMAX=128,KXMAX=8,NYMAX=128,KYMAX=8)

c
c     dimensions should be
c                  work1(max(nx,ny),max(nx,ny))
c                  work2(max(nx,ny))
c                  work3(max((2*kx-1)*nx,(2*ky-1)*ny))
c

      DIMENSION WORK1(NXMAX,NYMAX),WORK2(NXMAX)
      DIMENSION WORK3((2*KXMAX-1)*NXMAX)

      DIMENSION XVEC(NX),XKNOT(NX+KX),YVEC(NY),YKNOT(NY+KY)
      DIMENSION XYDATA(LDF,*),BCOEF(NX,NY)


      IF ((KX .GT. KXMAX) .OR. (NX .GT. NXMAX)) THEN
         WRITE(6,*) 'SUBROUTINE DBS2IN: ERROR'
         WRITE(6,*) 'KX > KXMAX OR NX > NXMAX'
         WRITE(6,*) 'KX = ',KX,'  KXMAX = ',KXMAX
         WRITE(6,*) 'NX = ',NX,'  NXMAX = ',NXMAX
         STOP
      ENDIF

      IF ((KY .GT. KYMAX) .OR. (NY .GT. NYMAX)) THEN
         WRITE(6,*) 'SUBROUTINE DBS2IN: ERROR'
         WRITE(6,*) 'KY > KYMAX OR NY > NYMAX'
         WRITE(6,*) 'KY = ',KY,'  KYMAX = ',KYMAX
         WRITE(6,*) 'NY = ',NY,'  NYMAX = ',NYMAX
         STOP
      ENDIF

      CALL SPLI2D(XVEC,LDF,XYDATA,XKNOT,NX,KX,NY,WORK2,WORK3,WORK1)
      CALL SPLI2D(YVEC,NY, WORK1, YKNOT,NY,KY,NX,WORK2,WORK3,BCOEF)

      RETURN
      END


C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      SUBROUTINE SPLI2D(XYVEC,LD,XYDATA,XYKNOT,N,K,M,WORK2,WORK3,BCOEF)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)

      DIMENSION XYVEC(N),XYKNOT(N+K),XYDATA(LD,M),BCOEF(M,N)
      DIMENSION WORK2(N),WORK3((2*K-1)*N)

      NP1   = N + 1
      KM1   = K - 1
      KPKM2 = 2 * KM1
      LEFT  = K
      LENQ  = N * (K + KM1)

      DO 10 I = 1,LENQ
         WORK3(I) = 0.0D0
 10   CONTINUE

      DO 20 I = 1,N
         XYVECI  = XYVEC(I)
         ILP1MX = MIN0(I+K,NP1)
         LEFT   = MAX0(LEFT,I)
         IF (XYVECI .LT. XYKNOT(LEFT)) GO TO 998
 30      IF (XYVECI .LT. XYKNOT(LEFT+1)) GO TO 40
         LEFT = LEFT + 1
         IF (LEFT .LT. ILP1MX) GO TO 30
         LEFT = LEFT - 1
         IF (XYVECI .GT. XYKNOT(LEFT+1)) GO TO 998
 40      CALL BSPLVB(XYKNOT,N+K,K,1,XYVECI,LEFT,WORK2)
         JJ = I - LEFT + 1 + (LEFT - K) * (K + KM1)
         DO 50 J = 1,K
            JJ        = JJ + KPKM2
            WORK3(JJ) = WORK2(J)
 50      CONTINUE
 20   CONTINUE

      CALL BANFAC(WORK3,K+KM1,N,KM1,KM1,IFLAG )

      GO TO (60,999), IFLAG

 60   DO 70 J = 1,M
         DO 80 I = 1,N
            WORK2(I) = XYDATA(I,J)
 80      CONTINUE

         CALL BANSLV(WORK3,K+KM1,N,KM1,KM1,WORK2)

         DO 90 I = 1,N
            BCOEF(J,I) = WORK2(I)
 90      CONTINUE
 70   CONTINUE

      RETURN

 998  WRITE(6,*) 'SUBROUTINE DB2IN:'
      WRITE(6,*) 'I WITH KNOT(I) <= X/Y < KNOT(I+1) REQUIRED.'
      WRITE(6,*) 'KNOT(1)   = ', XYKNOT(1)
      WRITE(6,*) 'KNOT(N+K) = ', XYKNOT(N+K)
      WRITE(6,*) '      X/Y = ', XYVECI

      STOP

 999  WRITE(6,*) 'SUBROUTINE DBS2IN: ERROR'
      WRITE(6,*) 'NO SOLUTION OF LINEAR EQUATION SYSTEM !!!'

      STOP

      END


c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      function dbs2vl(x,y,kx,ky,xknot,yknot,nx,ny,bcoef)

c
c  evaluates a two-dimensional tensor-product spline, given its
c  tensor-product B-spline representation.    use numeric
c
c   x      - x-coordinate of the point at which the spline is to be
c            evaluated.  (input)
c   y      - y-coordinate of the point at which the spline is to be
c            evaluated.  (input)
c   kx     - order of the spline in the x-direction.  (input)
c   ky     - order of the spline in the y-direction.  (input)
c   xknot  - array of length nx+kx containing the knot
c            sequence in the x-direction.  (input)
c            xknot must be nondecreasing.
c   yknot  - array of length ny+ky containing the knot
c            sequence in the y-direction.  (input)
c            yknot must be nondecreasing.
c   nx     - number of B-spline coefficients in the x-direction.
c            (input)
c   ny     - number of B-spline coefficients in the y-direction.
c            (input)
c   bcoef  - array of length nx*ny containing the
c            tensor-product B-spline coefficients.  (input)
c            bscoef is treated internally as a matrix of size nx
c            by ny.
c   dbs2vl - value of the spline at (x,y).  (output)
c

      implicit double precision (a-h,o-z), integer (i-n)

      parameter(kxmax=8,kymax=8)

      dimension xknot(nx+kx),yknot(ny+ky),bcoef(nx,ny)
      dimension work(kymax)

c
c     check if k <= kmax
c

      if(kx .gt. kxmax) then
         write(6,*) 'subroutine dbs2vl:'
         write(6,*) 'kx <= kxmax required.'
         stop
      endif

      if(ky .gt. kymax) then
         write(6,*) 'subroutine dbs2vl:'
         write(6,*) 'ky <= kymax required.'
         stop
      endif

c
c     check if knot(i) <= knot(i+1) and calculation of i so that
c     knot(i) <= x < knot(i+1)
c

      leftx = 0

      do 10 i = 1,nx+kx-1
         if (xknot(i) .gt. xknot(i+1)) then
             write(6,*) 'subroutine dbs2vl:'
             write(6,*) 'xknot(i) <= xknot(i+1) required.'
             write(6,*) i, xknot(i), xknot(i+1)
             write(6,*)
             write(6,*) xknot
          stop
          endif
         if((xknot(i) .le. x) .and. (x .lt. xknot(i+1))) leftx = i
10    continue

      if(leftx .eq. 0) then
         write(6,*) 'subroutine dbs2vl:'
         write(6,*) 'i with xknot(i) <= x < xknot(i+1) required.'
         write(6,*) 'x = ', x
         write(6,*)
         write(6,*) xknot
         stop
      endif

      lefty = 0

      do 20 i = 1,ny+ky-1
         if (yknot(i) .gt. yknot(i+1)) then
             write(6,*) 'subroutine dbs2vl:'
             write(6,*) 'yknot(i) <= yknot(i+1) required.'
             write(6,*) i, yknot(i), yknot(i+1)
          stop
          endif
         if((yknot(i) .le. y) .and. (y .lt. yknot(i+1))) lefty = i
 20   continue

      if(lefty .eq. 0) then
         write(6,*) 'subroutine dbs2vl:'
         write(6,*) 'i with yknot(i) <= y < yknot(i+1) required.'
         write(6,*) 'yknot(i)   = ', yknot(i)
         write(6,*) '  y        = ', y
         write(6,*) 'yknot(i+1) = ', yknot(i+1)
         stop
      endif

      do 30 iky = 1,ky
         work(iky) = dbsdca(0,x,kx,xknot,nx,bcoef(1,lefty-ky+iky),leftx)
 30   continue

      dbs2vl = dbsval(y,ky,yknot(lefty-ky+1),ky,work)

      return

      end


c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      function dbs2dr(iderx,idery,x,y,kx,ky,xknot,yknot,nx,ny,bcoef)

c
c  Evaluates the derivative of a two-dimensional tensor-product spline,
c  given its tensor-product B-spline representation.
c
c   iderx  - order of the derivative in the x-direction.  (input)
c   idery  - order of the derivative in the y-direction.  (input)
c   x      - x-coordinate of the point at which the spline is to be
c            evaluated.  (input)
c   y      - y-coordinate of the point at which the spline is to be
c            evaluated.  (input)
c   kx     - order of the spline in the x-direction.  (input)
c   ky     - order of the spline in the y-direction.  (input)
c   xknot  - array of length nx+kx containing the knot
c            sequence in the x-direction.  (input)
c            xknot must be nondecreasing.
c   yknot  - array of length ny+ky containing the knot
c            sequence in the y-direction.  (input)
c            yknot must be nondecreasing.
c   nx     - number of B-spline coefficients in the x-direction.
c            (input)
c   ny     - number of B-spline coefficients in the y-direction.
c            (input)
c   bcoef  - array of length nx*ny containing the
c            tensor-product B-spline coefficients.  (input)
c            bscoef is treated internally as a matrix of size nx
c            by ny.
c   dbs2dr  - value of the (iderx,idery) derivative of the spline at
c            (x,y).  (output)
c

      implicit double precision (a-h,o-z), integer (i-n)

      parameter(kxmax=8,kymax=8)

      dimension xknot(nx+kx),yknot(ny+ky),bcoef(nx,ny)
      dimension work(kymax)

c
c     check if k <= kmax
c
      if(kx .gt. kxmax) then
         write(6,*) 'subroutine dbs2dr:'
         write(6,*) 'kx <= kxmax required.'
         stop
      endif

      if(ky .gt. kymax) then
         write(6,*) 'subroutine dbs2dr:'
         write(6,*) 'ky <= kymax required.'
         stop
      endif

c
c     check if knot(i) <= knot(i+1) and calculation of i so that
c     knot(i) <= x < knot(i+1)
c

      nintx = 0

      do 10 i = 1,nx+kx-1
         if (xknot(i) .gt. xknot(i+1)) then
            write(6,*) 'subroutine dbs2dr:'
            write(6,*) 'xknot(i) <= xknot(i+1) required.'
            write(6,*) i, xknot(i), xknot(i+1)
            stop
         endif
         if((xknot(i) .le. x) .and. (x .lt. xknot(i+1))) nintx = i
 10   continue

      if(nintx .eq. 0) then
         write(6,*) 'subroutine dbs2dr:'
         write(6,*) 'i with xknot(i) <= x < xknot(i+1) required.'
         write(6,*) 'x = ', x
         stop
      endif

      ninty = 0

      do 20 i = 1,ny+ky-1
         if (yknot(i) .gt. yknot(i+1)) then
             write(6,*) 'subroutine dbs2dr:'
             write(6,*) 'yknot(i) <= yknot(i+1) required.'
             write(6,*) i, yknot(i), yknot(i+1)
          stop
          endif
         if((yknot(i) .le. y) .and. (y .lt. yknot(i+1))) ninty = i
 20   continue

      if(ninty .eq. 0) then
         write(6,*) 'subroutine dbs2dr:'
         write(6,*) 'i with yknot(i) <= y < yknot(i+1) required.'
         write(6,*) 'y = ', y
         stop
      endif

      do 30 iy = 1, ky
         work(iy) =
     .        dbsdca(iderx,x,kx,xknot,nx,bcoef(1,ninty-ky+iy),nintx)
 30   continue

      dbs2dr = dbsder(idery,y,ky,yknot(ninty-ky+1),ky,work)

      return

      end


c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dbs2gd

      subroutine dbs2gd(iderx,idery,nxvec,xvec,nyvec,yvec,kx,
     .     ky,xknot,yknot,nx,ny,bcoef,val,ldvalue)

c
c  Evaluates the derivative of a two-dimensional tensor-product spline,
c  given its tensor-product B-spline representation on a grid.
c
c   iderx   - order of the derivative in the x-direction.  (input)
c   idery   - order of the derivative in the y-direction.  (input)
c   nxvec   - number of grid points in the x-direction.  (input)
c   xvec    - array of length nx containing the x-coordinates at
c             which the spline is to be evaluated.  (input)
c             the points in xvec should be strictly increasing.
c   nyvec   - number of grid points in the y-direction.  (input)
c   yvec    - array of length ny containing the y-coordinates at
c             which the spline is to be evaluated.  (input)
c             the points in yvec should be strictly increasing.
c   kx      - order of the spline in the x-direction.  (input)
c   ky      - order of the spline in the y-direction.  (input)
c   xknot   - array of length nx+kx containing the knot
c             sequence in the x-direction.  (input)
c             xknot must be nondecreasing.
c   yknot   - array of length ny+ky containing the knot
c             sequence in the y-direction.  (input)
c             yknot must be nondecreasing.
c   nx      - number of B-spline coefficients in the x-direction.
c             (input)
c   ny      - number of B-spline coefficients in the y-direction.
c             (input)
c   bcoef   - array of length nx*ny containing the
c             tensor-product B-spline coefficients.  (input)
c             bscoef is treated internally as a matrix of size nx
c             by ny.
c   val     - value of the (iderx,idery) derivative of the spline on
c             the nx by ny grid.  (output)
c             value(i,j) contains the derivative of the spline at the
c             point (xvec(i),yvec(j)).
c   ldf     - leading dimension of value exactly as specified in the
c             dimension statement of the calling program.  (input)
c

      implicit double precision (a-h,o-z), integer (i-n)

      parameter(kxmax=8,kymax=8)
      parameter(nxmax=500,nymax=500)

      dimension xvec(nxvec),xknot(nx+kx)
      dimension yvec(nyvec),yknot(ny+ky)
      dimension bcoef(nx,ny),val(ldvalue,*)

      dimension dl(nxmax,kxmax),dr(nxmax,kxmax),save1(nxmax)
      dimension biatx(nxmax,kxmax),biaty(nymax,kymax)
      dimension leftx(nxmax),lefty(nymax),term(nxmax)
      dimension work(kymax)

      logical same,next

c
c     check if kx <= kxmax
c

      if(kx .gt. kxmax) then
         write(6,*) 'subroutine dbs2gd:'
         write(6,*) 'kx <= kxmax required.'
         stop
      endif

      leftx(1) = 0

      call huntn(xknot,nx+kx,kx,xvec(1),leftx(1))

      do 10 ix = 2,nxvec
         leftx(ix) = leftx(ix-1)
         same = (xknot(leftx(ix)) .le. xvec(ix))
     .        .and. (xvec(ix) .le. xknot(leftx(ix)+1))
         if(.not. same ) then
            leftx(ix) = leftx(ix) + 1
            next      = (xknot(leftx(ix)) .le. xvec(ix))
     .           .and. (xvec(ix) .le. xknot(leftx(ix)+1))
            if (.not. next)
     .           call huntn(xknot,nx+kx,kx,xvec(ix),leftx(ix))
         endif
 10   continue

      do 20 i = 1,nx+kx-1
         if (xknot(i) .gt. xknot(i+1)) then
            write(6,*) 'subroutine dbs2gd:'
            write(6,*) 'xknot(i) <= xknot(i+1) required.'
            write(6,*) i, xknot(i), xknot(i+1)
            write(6,*)
            write(6,*) xknot
            stop
         endif
 20   continue

      do 30 i = 1,nxvec
         if ((xvec(i).lt.xknot(1)).or.(xvec(i).gt.xknot(nx+kx))) then
            write(6,*) 'subroutine dbs2gd:'
            write(6,*) 'ix with xknot(ix) <= x < xknot(ix+1) required.'
            write(6,*) 'x = ', xvec(i)
            stop
         endif
 30   continue

c
c     check if ky <= kymax
c

      if(ky .gt. kymax) then
         write(6,*) 'subroutine dbs2gd:'
         write(6,*) 'ky <= kymax required.'
         stop
      endif

      lefty(1) = 0

      call huntn(yknot,ny+ky,ky,yvec(1),lefty(1))

      do 40 iy = 2,nyvec
         lefty(iy) = lefty(iy-1)
         same = (yknot(lefty(iy)) .le. yvec(iy))
     .        .and. (yvec(iy) .le. yknot(lefty(iy)+1))
         if(.not. same ) then
            lefty(iy) = lefty(iy) + 1
            next      = (yknot(lefty(iy)) .le. yvec(iy))
     .           .and. (yvec(iy) .le. yknot(lefty(iy)+1))
            if (.not. next)
     .           call huntn(yknot,ny+ky,ky,yvec(iy),lefty(iy))
         endif
 40   continue

      do 50 i = 1,ny+ky-1
         if (yknot(i) .gt. yknot(i+1)) then
            write(6,*) 'subroutine dbs2gd:'
            write(6,*) 'yknot(i) <= yknot(i+1) required.'
            write(6,*) i, yknot(i), yknot(i+1)
            write(6,*)
            write(6,*) yknot
            stop
         endif
 50   continue

      do 60 i = 1,nyvec
         if ((yvec(i).lt.yknot(1)).or.(yvec(i).gt.yknot(ny+ky))) then
            write(6,*) 'subroutine dbs2gd:'
            write(6,*) 'iy with yknot(iy) <= y < yknot(iy+1) required.'
            write(6,*) 'y = ', yvec(i)
            stop
         endif
 60   continue

      if ((iderx .eq. 0) .and. (idery .eq. 0)) then

         do 70 ix = 1,nxvec
            biatx(ix,1) = 1.d0
 70      continue

         do 80 ik = 1,kx-1
            do 90 ix = 1,nxvec
               dr(ix,ik) = xknot(leftx(ix)+ik) - xvec(ix)
               dl(ix,ik) = xvec(ix) - xknot(leftx(ix)+1-ik)
               save1(ix) = 0.d0
 90         continue

            do 100 il = 1,ik
               do 110 ix = 1,nxvec
                  term(ix)     = biatx(ix,il)
     .                 / (dr(ix,il) + dl(ix,ik+1-il))
                  biatx(ix,il) = save1(ix) + dr(ix,il) * term(ix)
                  save1(ix)    = dl(ix,ik+1-il) * term(ix)
 110           continue
 100        continue

            do 120 ix = 1,nxvec
               biatx(ix,ik+1) = save1(ix)
 120        continue
 80      continue

         do 130 iy = 1,nyvec
            biaty(iy,1) = 1.d0
 130      continue

         do 140 ik = 1,ky-1
            do 150 iy = 1,nyvec
               dr(iy,ik) = yknot(lefty(iy)+ik) - yvec(iy)
               dl(iy,ik) = yvec(iy) - yknot(lefty(iy)+1-ik)
               save1(iy) = 0.d0
 150         continue

            do 160 il = 1,ik
               do 170 iy = 1,nyvec
                  term(iy)     = biaty(iy,il)
     .                 / (dr(iy,il) + dl(iy,ik+1-il))
                  biaty(iy,il) = save1(iy) + dr(iy,il) * term(iy)
                  save1(iy)    = dl(iy,ik+1-il) * term(iy)
 170           continue
 160        continue

            do 180 iy = 1,nyvec
               biaty(iy,ik+1) = save1(iy)
 180        continue
 140     continue

         do 190 iy = 1,nyvec
            do 200 ix = 1,nxvec
               val(ix,iy) = 0.0d0
 200        continue
 190     continue

         do 210 iky = 1,ky
            do 220 ikx = 1,kx
               do 230 iy = 1,nyvec
                  do 240 ix = 1,nxvec
                     val(ix,iy) = val(ix,iy)
     .                    + biatx(ix,ikx) * biaty(iy,iky)
     .                    * bcoef(leftx(ix)-kx+ikx,lefty(iy)-ky+iky)
 240              continue
 230           continue
 220        continue
 210     continue


      elseif (((iderx .ge. 1) .or. (idery .ge. 1))
     .     .and. ( (iderx .lt. kx) .and. (idery .lt. ky))) then

         do 250 iy = 1,nyvec
            do 260 ix = 1,nxvec
               do 270 iky = 1, ky
                  work(iky) = dbsdca(iderx,xvec(ix),kx,xknot,nx,
     .                 bcoef(1,lefty(iy)-ky+iky),leftx(ix))
 270           continue
               val(ix,iy) = dbsder(idery,yvec(iy),ky,
     .              yknot(lefty(iy)-ky+1),ky,work)
 260        continue
 250     continue

c               val(ix,iy) = dbs2dr(iderx,idery,xvec(ix),yvec(iy),
c     .              kx,ky,xknot,yknot,nx,ny,bcoef)

      else

         do 280 iy = 1,nyvec
            do 290 ix = 1,nxvec
               val(ix,iy) = 0.0d0
 290        continue
 280     continue

      endif

      return

      end


c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine dbs3in(nx,xvec,ny,yvec,nz,zvec,xyzdata,ldf,mdf,
     .     kx,ky,kz,xknot,yknot,zknot,bcoef)

c
c  Computes a three-dimensional tensor-product spline interpolant,
c  returning the tensor-product B-spline coefficients.
c
c   nx      - number of data points in the x-direction.  (input)
c   xvec    - array of length nxdata containing the data points in
c             the x-direction.  (input)
c             xdata must be increasing.
c   ny      - number of data points in the y-direction.  (input)
c   yvec    - array of length nydata containing the data points in
c             the y-direction.  (input)
c             ydata must be increasing.
c   nz      - number of data points in the z-direction.  (input)
c   zvec    - array of length nzdata containing the data points in
c             the z-direction.  (input)
c             zdata must be increasing.
c   xyzdata - array of size nx by ny by nz containing the
c             values to be interpolated.  (input)
c             xyzdata(i,j,k) contains the value at
c             (xvec(i),yvec(j),zvec(k)).
c   ldf     - leading dimension of fdata exactly as specified in the
c             dimension statement of the calling program.  (input)
c   mdf     - middle dimension of fdata exactly as specified in the
c             dimension statement of the calling program.  (input)
c   kx      - order of the spline in the x-direction.  (input)
c             kxord must be less than or equal to nxdata.
c   ky      - order of the spline in the y-direction.  (input)
c             kyord must be less than or equal to nydata.
c   kz      - order of the spline in the z-direction.  (input)
c             kzord must be less than or equal to nzdata.
c   xknot   - array of length nx+kx containing the knot
c             sequence in the x-direction.  (input)
c             xknot must be nondecreasing.
c   yknot   - array of length ny+ky containing the knot
c             sequence in the y-direction.  (input)
c             yknot must be nondecreasing.
c   zknot   - array of length nz+kz containing the knot
c             sequence in the z-direction.  (input)
c             zknot must be nondecreasing.
c   bcoef   - array of length nx*ny*nz containing the
c             tensor-product B-spline coefficients.  (output)
c             bscoef is treated internally as a matrix of size nx
c             by ny by nz.
c

      implicit double precision (a-h,o-z), integer (i-n)

      parameter(nxmax=56,kxmax=8,nymax=56,kymax=8,nzmax=56,kzmax=8)

c
c     dimensions should be
c              work1(nx,ny,nz)
c              work2(nz)
c              work3((2*kz-1)*nz)
c

      dimension work2(nzmax),work3((2*kzmax-1)*nzmax)
      dimension work1(nxmax,nymax,nzmax)

      dimension bcoef(nx,ny,nz),xvec(nx),yvec(ny),zvec(nz)
      dimension xknot(nx+kx),yknot(ny+ky),zknot(nz+kz)
      dimension xyzdata(ldf,mdf,nz)

      if ((kx .gt. kxmax) .or. (nx .gt. nxmax)) then
         write(6,*) 'subroutine dbs3in: error'
         write(6,*) 'kx > kxmax or nx > nxmax'
         write(6,*) 'kx = ',kx,'  kxmax = ',kxmax
         write(6,*) 'nx = ',nx,'  nxmax = ',nxmax
         stop
      endif

      if ((ky .gt. kymax) .or. (ny .gt. nymax)) then
         write(6,*) 'subroutine dbs3in: error'
         write(6,*) 'ky > kymax or ny > nymax'
         write(6,*) 'ky = ',ky,'  kymax = ',kymax
         write(6,*) 'ny = ',ny,'  nymax = ',nymax
         stop
      endif

      if ((kz .gt. kzmax) .or. (nz .gt. nzmax)) then
         write(6,*) 'subroutine dbs3in: error'
         write(6,*) 'kz > kzmax or nz > nzmax'
         write(6,*) 'kz = ',kz,'  kymax = ',kzmax
         write(6,*) 'nz = ',nz,'  nymax = ',nzmax
         stop
      endif

      call spli3d(zvec,ldf,mdf,xyzdata,zknot,nz,kz,nx,ny,work2,work3,
     .   work1,nxmax,nymax,nzmax)

      do 10 iz = 1,nz
         call dbs2in(nx,xvec,ny,yvec,work1(1,1,iz),nxmax,kx,ky,xknot,
     .        yknot,bcoef(1,1,iz))
 10   continue

      return

      end


c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      subroutine spli3d(xyzvec,ldf,mdf,xyzdata,xyzknot,n,k,m,l,work2,
     .		work3,bcoef,nxmax,nymax,nzmax)

      implicit double precision (a-h,o-z), integer (i-n)

      dimension xyzvec(n),xyzknot(n+k),xyzdata(ldf,mdf,*)
      dimension bcoef(nxmax,nymax,nzmax), work2(n),work3((2*k-1)*n)

      np1   = n + 1
      km1   = k - 1
      kpkm2 = 2 * km1
      left  = k
      lenq  = n * (k + km1)

      do 10 i = 1,lenq
         work3(i) = 0.d0
 10   continue

      do 20 i = 1,n
         xyzveci = xyzvec(i)
         ilp1mx  = min0(i+k,np1)
         left    = max0(left,i)
         if (xyzveci .lt. xyzknot(left)) go to 998
 30      if (xyzveci .lt. xyzknot(left+1)) go to 40
         left = left + 1
         if (left .lt. ilp1mx) go to 30
         left = left - 1
         if (xyzveci .gt. xyzknot(left+1)) go to 998
 40      call bsplvb(xyzknot,n+k,k,1,xyzveci,left,work2)
         jj = i - left + 1 + (left - k) * (k + km1)
         do 50 j = 1,k
            jj    = jj + kpkm2
            work3(jj) = work2(j)
 50      continue
 20   continue

      call banfac(work3,k+km1,n,km1,km1,iflag)

      go to (60,999), iflag

 60   do 70 j = 1,l
         do 80 i = 1,m
            do 90 in = 1,n
               work2(in) = xyzdata(i,j,in)
 90         continue

            call banslv(work3,k+km1,n,km1,km1,work2)

            do 100 in = 1,n
               bcoef(i,j,in) = work2(in)
 100        continue
 80      continue
 70   continue

      return

 998  write(6,*) 'subroutine db3in:'
      write(6,*) 'i with knot(i) <= x/y/z < knot(i+1) required.'
      write(6,*) 'knot(1)   = ', xyzknot(1)
      write(6,*) 'knot(n+k) = ', xyzknot(n+k)
      write(6,*) '    x/y/z = ', xyzveci

      stop

 999  write(6,*) 'subroutine dbs3in: error'
      write(6,*) 'no solution of linear equation system !!!'

      stop

      end


c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      function dbs3vl(x,y,z,kx,ky,kz,xknot,yknot,zknot,nx,ny,nz,bcoef)

c
c  Evaluates a three-dimensional tensor-product spline, given its
c  tensor-product B-spline representation.
c
c   x      - x-coordinate of the point at which the spline is to be
c            evaluated.  (input)
c   y      - y-coordinate of the point at which the spline is to be
c            evaluated.  (input)
c   z      - z-coordinate of the point at which the spline is to be
c            evaluated.  (input)
c   kx     - order of the spline in the x-direction.  (input)
c   ky     - order of the spline in the y-direction.  (input)
c   kz     - order of the spline in the z-direction.  (input)
c   xknot  - array of length nx+kx containing the knot
c            sequence in the x-direction.  (input)
c            xknot must be nondecreasing.
c   yknot  - array of length ny+ky containing the knot
c            sequence in the y-direction.  (input)
c            yknot must be nondecreasing.
c   zknot  - array of length nz+kz containing the knot
c            sequence in the z-direction.  (input)
c            zknot must be nondecreasing.
c   nx     - number of B-spline coefficients in the x-direction.
c            (input)
c   ny     - number of B-spline coefficients in the y-direction.
c            (input)
c   nz     - number of B-spline coefficients in the z-direction.
c            (input)
c   bcoef  - array of length nx*ny*nz containing the
c            tensor-product B-spline coefficients.  (input)
c            bscoef is treated internally as a matrix of size nx
c            by ny by nz.
c   dbs3vl - value of the spline at (x,y,z).  (output)
c

      implicit double precision (a-h,o-z), integer (i-n)

      parameter(kxmax=8,kymax=8,kzmax=8)

c
c     dimension should be
c            dimension work(kz)
c

      dimension work(kzmax)

      dimension xknot(nx+kx),yknot(ny+ky),zknot(nz+kz)
      dimension bcoef(nx,ny,nz)

c     check if k <= kmax

      if(kx .gt. kxmax) then
         write(6,*) 'subroutine dbs3vl:'
         write(6,*) 'kx <= kxmax required.'
         stop
      endif

      if(ky .gt. kymax) then
         write(6,*) 'subroutine dbs3vl:'
         write(6,*) 'ky <= kymax required.'
         stop
      endif

      if(kz .gt. kzmax) then
         write(6,*) 'subroutine dbs3vl:'
         write(6,*) 'kz <= kzmax required.'
         stop
      endif

c
c     check if knot(i) <= knot(i+1) and calculation of i so that
c     knot(i) <= x < knot(i+1)
c


      nintz = 0

      do 10 i = 1,nz+kz-1
         if (zknot(i) .gt. zknot(i + 1)) then
             write(6,*) 'subroutine dbs3vl:'
             write(6,*) 'zknot(i) <= zknot(i+1) required.'
             write(6,*) i, zknot(i), zknot(i+1)
          stop
          endif
         if((zknot(i) .le. z) .and. (z .lt. zknot(i + 1))) nintz = i
 10   continue

      if(nintz .eq. 0) then
         write(6,*) 'subroutine dbs3vl:'
         write(6,*) 'i with zknot(i) <= z < zknot(i+1) required.'
         write(6,*) 'zknot(i)   = ', zknot(i)
         write(6,*) '  z        = ', z
         write(6,*) 'zknot(i+1) = ', zknot(i+1)
         stop
      endif

      do 40 jz = 1,kz
         work(jz) = dbs2vl(x,y,kx,ky,xknot,yknot,nx,ny,
     .        bcoef(1,1,nintz-kz+jz))
 40   continue

      dbs3vl = dbsval(z,kz,zknot(nintz-kz+1),kz,work)

      return

      end


c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      function dbs3dr(iderx,idery,iderz,x,y,z,kx,ky,kz,
     .     xknot,yknot,zknot,nx,ny,nz,bcoef)

c
c  Evaluates the derivative of a three-dimensional tensor-product spline,
c  given its tensor-product B-spline representation.
c
c   iderx  - order of the x-derivative.  (input)
c   idery  - order of the y-derivative.  (input)
c   iderz  - order of the z-derivative.  (input)
c   x      - x-coordinate of the point at which the spline is to be
c            evaluated.  (input)
c   y      - y-coordinate of the point at which the spline is to be
c            evaluated.  (input)
c   z      - z-coordinate of the point at which the spline is to be
c            evaluated.  (input)
c   kx     - order of the spline in the x-direction.  (input)
c   ky     - order of the spline in the y-direction.  (input)
c   kz     - order of the spline in the z-direction.  (input)
c   xknot  - array of length nx+kx containing the knot
c            sequence in the x-direction.  (input)
c            xknot must be nondecreasing.
c   yknot  - array of length ny+ky containing the knot
c            sequence in the y-direction.  (input)
c            yknot must be nondecreasing.
c   zknot  - array of length nz+kz containing the knot
c            sequence in the z-direction.  (input)
c            zknot must be nondecreasing.
c   nx     - number of B-spline coefficients in the x-direction.
c            (input)
c   ny     - number of B-spline coefficients in the y-direction.
c            (input)
c   nz     - number of B-spline coefficients in the z-direction.
c            (input)
c   bcoef  - array of length nx*ny*nz containing the
c            tensor-product B-spline coefficients.  (input)
c            bscoef is treated internally as a matrix of size nx
c            by ny by nz.
c   dbs3dr - value of the (iderx,idery,iderz) derivative of the
c            spline at (x,y,z).  (output)
c

      implicit double precision (a-h,o-z), integer (i-n)

      parameter(kxmax=8,kymax=8,kzmax=8)

      dimension xknot(nx+kx)
      dimension yknot(ny+ky)
      dimension zknot(nz+kz)
      dimension bcoef(nx,ny,nz),work(kzmax)

c
c     check if k <= kmax
c
      if(kx .gt. kxmax) then
         write(6,*) 'subroutine dbs3dr:'
         write(6,*) 'kx <= kxmax required.'
         stop
      endif

      if(ky .gt. kymax) then
         write(6,*) 'subroutine dbs3dr:'
         write(6,*) 'ky <= kymax required.'
         stop
      endif

      if(kz .gt. kzmax) then
         write(6,*) 'subroutine dbs3dr:'
         write(6,*) 'kz <= kzmax required.'
         stop
      endif

c
c     check if knot(i) <= knot(i+1) and calculation of i so that
c     knot(i) <= x < knot(i+1)
c

      nintz = 0

      do 10 i = 1,nz+kz-1
         if (zknot(i) .gt. zknot(i + 1)) then
             write(6,*) 'subroutine dbs3vl:'
             write(6,*) 'zknot(i) <= zknot(i+1) required.'
             write(6,*) i, zknot(i), zknot(i+1)
          stop
          endif
         if((zknot(i) .le. z) .and. (z .lt. zknot(i + 1))) nintz = i
 10   continue

      if(nintz .eq. 0) then
         write(6,*) 'subroutine dbs3dr:'
         write(6,*) 'i with zknot(i) <= z < zknot(i+1) required.'
         write(6,*) 'zknot(i)   = ', zknot(i)
         write(6,*) '  z        = ', z
         write(6,*) 'zknot(i+1) = ', zknot(i+1)
         stop
      endif

      do 20 jz = 1,kz
         work(jz) = dbs2dr(iderx,idery,x,y,kx,ky,xknot,yknot,nx,ny,
     .        bcoef(1,1,nintz-kz+jz))
 20   continue

      dbs3dr = dbsder(iderz,z,kz,zknot(nintz-kz+1),kz,work)

      return

      end


c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      subroutine dbs3gd(iderx,idery,iderz,nxvec,xvec,nyvec,yvec,nzvec,
     .     zvec,kx,ky,kz,xknot,yknot,zknot,nx,ny,
     .     nz,bcoef,value,ldvalue,mdvalue)

c
c  Evaluates the derivative of a three-dimensional tensor-product spline,
c  given its tensor-product B-spline representation on a grid.
c
c   iderx  - order of the x-derivative.  (input)
c   idery  - order of the y-derivative.  (input)
c   iderz  - order of the z-derivative.  (input)
c   nx     - number of grid points in the x-direction.  (input)
c   xvec   - array of length nx containing the x-coordinates at
c            which the spline is to be evaluated.  (input)
c            the points in xvec should be strictly increasing.
c   ny     - number of grid points in the y-direction.  (input)
c   yvec   - array of length ny containing the y-coordinates at
c            which the spline is to be evaluated.  (input)
c            the points in yvec should be strictly increasing.
c   nz     - number of grid points in the z-direction.  (input)
c   zvec   - array of length nz containing the z-coordinates at
c            which the spline is to be evaluated.  (input)
c            the points in yvec should be strictly increasing.
c   kx     - order of the spline in the x-direction.  (input)
c   ky     - order of the spline in the y-direction.  (input)
c   kz     - order of the spline in the z-direction.  (input)
c   xknot  - array of length nx+kx containing the knot
c            sequence in the x-direction.  (input)
c            xknot must be nondecreasing.
c   yknot  - array of length ny+ky containing the knot
c            sequence in the y-direction.  (input)
c            yknot must be nondecreasing.
c   zknot  - array of length nz+kz containing the knot
c            sequence in the z-direction.  (input)
c            zknot must be nondecreasing.
c   nx     - number of B-spline coefficients in the x-direction.
c            (input)
c   ny     - number of B-spline coefficients in the y-direction.
c            (input)
c   nz     - number of B-spline coefficients in the z-direction.
c            (input)
c   bcoef  - array of length nx*ny*nz containing the
c            tensor-product B-spline coefficients.  (input)
c            bscoef is treated internally as a matrix of size nx
c            by ny by nz.
c   val    - array of size nx by ny by nz containing the values of
c            the (iderx,idery,iderz) derivative of the spline on the
c            nx by ny by nz grid.  (output)
c            value(i,j,k) contains the derivative of the spline at
c            the point (xvec(i), yvec(j), zvec(k)).
c   ldf    - leading dimension of value exactly as specified in the
c            dimension statement of the calling program.  (input)
c   mdf    - middle dimension of value exactly as specified in the
c            dimension statement of the calling program.  (input)
c

      implicit double precision (a-h,o-z), integer (i-n)

      parameter(kxmax=8,kymax=8,kzmax=8)
      parameter(nxmax=56,nymax=56,nzmax=56)

      dimension xvec(nxvec),xknot(nx+kx)
      dimension yvec(nyvec),yknot(ny+ky)
      dimension zvec(nzvec),zknot(nz+kz)
      dimension bcoef(nx,ny,nz)
      dimension value(ldvalue,mdvalue,*)

      dimension dl(nxmax,kxmax),dr(nxmax,kxmax),save1(nxmax)
      dimension biatx(nxmax,kxmax),biaty(nymax,kymax)
      dimension biatz(nzmax,kzmax)
      dimension leftx(nxmax),lefty(nymax),leftz(nzmax)
      dimension term(nxmax)

      logical same,next

c
c     check if kx <= kxmax
c

      if(kx .gt. kxmax) then
         write(6,*) 'subroutine dbs3gd:'
         write(6,*) 'kx <= kxmax required.'
         stop
      endif

      do 10 i = 1,nx+kx-1
         if (xknot(i) .gt. xknot(i+1)) then
            write(6,*) 'subroutine dbs3gd:'
            write(6,*) 'xknot(i) <= xknot(i+1) required.'
            write(6,*) i, xknot(i), xknot(i+1)
            write(6,*)
            write(6,*) xknot
            stop
         endif
 10   continue

      do 20 i = 1,nxvec
         if ((xvec(i).lt.xknot(1)).or.(xvec(i).gt.xknot(nx+kx))) then
            write(6,*) 'subroutine dbs3gd:'
            write(6,*) 'ix with xknot(ix) <= x < xknot(ix+1) required.'
            write(6,*) 'x = ', xvec(i)
            stop
         endif
 20   continue

      leftx(1) = 0

      call huntn(xknot,nx+kx,kx,xvec(1),leftx(1))

      do 30 ix = 2,nxvec
         leftx(ix) = leftx(ix-1)
         same = (xknot(leftx(ix)) .le. xvec(ix))
     .        .and. (xvec(ix) .le. xknot(leftx(ix)+1))
         if(.not. same ) then
            leftx(ix) = leftx(ix) + 1
            next      = (xknot(leftx(ix)) .le. xvec(ix))
     .           .and. (xvec(ix) .le. xknot(leftx(ix)+1))
            if (.not. next)
     .           call huntn(xknot,nx+kx,kx,xvec(ix),leftx(ix))
         endif
 30   continue

c
c     check if ky <= kymax
c

      if(ky .gt. kymax) then
         write(6,*) 'subroutine dbs3gd:'
         write(6,*) 'ky <= kymax required.'
         stop
      endif

      do 40 i = 1,ny+ky-1
         if (yknot(i) .gt. yknot(i+1)) then
            write(6,*) 'subroutine dbs3gd:'
            write(6,*) 'yknot(i) <= yknot(i+1) required.'
            write(6,*) i, yknot(i), yknot(i+1)
            write(6,*)
            write(6,*) yknot
            stop
         endif
 40   continue

      do 50 i = 1,nyvec
         if ((yvec(i).lt.yknot(1)).or.(yvec(i).gt.yknot(ny+ky))) then
            write(6,*) 'subroutine dbs3gd:'
            write(6,*) 'iy with yknot(iy) <= y < yknot(iy+1) required.'
            write(6,*) 'y = ', yvec(i)
            stop
         endif
 50   continue

      lefty(1) = 0

      call huntn(yknot,ny+ky,ky,yvec(1),lefty(1))

      do 60 iy = 2,nyvec
         lefty(iy) = lefty(iy-1)
         same = (yknot(lefty(iy)) .le. yvec(iy))
     .        .and. (yvec(iy) .le. yknot(lefty(iy)+1))
         if(.not. same ) then
            lefty(iy) = lefty(iy) + 1
            next      = (yknot(lefty(iy)) .le. yvec(iy))
     .           .and. (yvec(iy) .le. yknot(lefty(iy)+1))
            if (.not. next)
     .           call huntn(yknot,ny+ky,ky,yvec(iy),lefty(iy))
         endif
 60   continue

c
c     check if kz <= kzmax
c

      if(kz .gt. kzmax) then
         write(6,*) 'subroutine dbs3gd:'
         write(6,*) 'kz <= kzmax required.'
         stop
      endif

      do 70 i = 1,nz+kz-1
         if (zknot(i) .gt. zknot(i+1)) then
            write(6,*) 'subroutine dbs3gd:'
            write(6,*) 'zknot(i) <= zknot(i+1) required.'
            write(6,*) i, zknot(i), zknot(i+1)
            write(6,*)
            write(6,*) zknot
            stop
         endif
 70   continue

      do 80 i = 1,nzvec
         if ((zvec(i).lt.zknot(1)).or.(zvec(i).gt.zknot(nz+kz))) then
            write(6,*) 'subroutine dbs3gd:'
            write(6,*) 'iz with zknot(iz) <= z < zknot(iz+1) required.'
            write(6,*) 'z = ', zvec(i)
            stop
         endif
 80   continue

      leftz(1) = 0

      call huntn(zknot,nz+kz,kz,zvec(1),leftz(1))

      do 90 iz = 2,nzvec
         leftz(iz) = leftz(iz-1)
         same = (zknot(leftz(iz)) .le. zvec(iz))
     .        .and. (zvec(iz) .le. zknot(leftz(iz)+1))
         if(.not. same ) then
            leftz(iz) = leftz(iz) + 1
            next      = (zknot(leftz(iz)) .le. zvec(iz))
     .           .and. (zvec(iz) .le. zknot(leftz(iz)+1))
            if (.not. next)
     .           call huntn(zknot,nz+kz,kz,zvec(iz),leftz(iz))
         endif
 90   continue

      if ((iderx .eq. 0) .and. (idery .eq. 0) .and. (iderz .eq.0)) then

         do 100 ix = 1,nxvec
            biatx(ix,1) = 1.d0
 100      continue

         do 110 ik = 1,kx-1
            do 120 ix = 1,nxvec
               dr(ix,ik) = xknot(leftx(ix)+ik) - xvec(ix)
               dl(ix,ik) = xvec(ix) - xknot(leftx(ix)+1-ik)
               save1(ix) = 0.d0
 120         continue

            do 130 il = 1,ik
               do 140 ix = 1,nxvec
                  term(ix)     = biatx(ix,il)
     .                 / (dr(ix,il) + dl(ix,ik+1-il))
                  biatx(ix,il) = save1(ix) + dr(ix,il) * term(ix)
                  save1(ix)    = dl(ix,ik+1-il) * term(ix)
 140           continue
 130        continue

            do 150 ix = 1,nxvec
               biatx(ix,ik+1) = save1(ix)
 150        continue
 110     continue

         do 160 iy = 1,nyvec
            biaty(iy,1) = 1.d0
 160     continue

         do 170 ik = 1,ky-1
            do 180 iy = 1,nyvec
               dr(iy,ik) = yknot(lefty(iy)+ik) - yvec(iy)
               dl(iy,ik) = yvec(iy) - yknot(lefty(iy)+1-ik)
               save1(iy) = 0.d0
 180        continue

            do 190 il = 1,ik
               do 200 iy = 1,nyvec
                  term(iy)     = biaty(iy,il)
     .                 / (dr(iy,il) + dl(iy,ik+1-il))
                  biaty(iy,il) = save1(iy) + dr(iy,il) * term(iy)
                  save1(iy)    = dl(iy,ik+1-il) * term(iy)
 200           continue
 190        continue

            do 210 iy = 1,nyvec
               biaty(iy,ik+1) = save1(iy)
 210        continue
 170     continue

         do 220 iz = 1,nzvec
            biatz(iz,1) = 1.d0
 220     continue

         do 230 ik = 1,kz-1
            do 240 iz = 1,nzvec
               dr(iz,ik) = zknot(leftz(iz)+ik) - zvec(iz)
               dl(iz,ik) = zvec(iz) - zknot(leftz(iz)+1-ik)
               save1(iz) = 0.d0
 240        continue

            do 250 il = 1,ik
               do 260 iz = 1,nzvec
                  term(iz)     = biatz(iz,il)
     .                 / (dr(iz,il) + dl(iz,ik+1-il))
                  biatz(iz,il) = save1(iz) + dr(iz,il) * term(iz)
                  save1(iz)    = dl(iz,ik+1-il) * term(iz)
 260           continue
 250        continue

            do 270 iz = 1,nzvec
               biatz(iz,ik+1) = save1(iz)
 270        continue
 230     continue


         do 280 iz = 1,nzvec
            do 290 iy = 1,nyvec
               do 300 ix = 1,nxvec
                  value(ix,iy,iz) = 0.0d0
 300           continue
 290        continue
 280     continue

         do 310 ikz = 1,kz
            do 320 iky = 1,ky
               do 330 ikx = 1,kx
                  do 340 iz = 1,nzvec
                     do 350 iy = 1,nyvec
                        do 360 ix = 1,nxvec
                           value(ix,iy,iz) = value(ix,iy,iz)
     .                          + biatx(ix,ikx) * biaty(iy,iky)
     .                          * biatz(iz,ikz)
     .                          * bcoef(leftx(ix)-kx+ikx,
     .                          lefty(iy)-ky+iky,leftz(iz)-kz+ikz)
 360                    continue
 350                 continue
 340              continue
 330           continue
 320        continue
 310     continue

      else

         do 370 iz = 1,nzvec
            do 380 iy = 1,nyvec
               do 390 ix = 1,nxvec
                  value(ix,iy,iz) = dbs3dr(iderx,idery,iderz,xvec(ix),
     .                 yvec(iy),zvec(iz),kx,ky,kz,xknot,yknot,
     .                 zknot,nx,ny,nz,bcoef)
 390           continue
 380         continue
 370      continue

      endif

      return

      end


c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      subroutine bsplvb(t,n,jhigh,index,x,left,biatx)

      implicit double precision (a-h,o-z), integer (i-n)

      parameter(kmax=8)

      dimension biatx(jhigh),t(n),dl(kmax),dr(kmax)

      data j/1/

      go to (10,20), index

 10   j = 1

      biatx(1) = 1.d0

      if (j .ge. jhigh) go to 99

 20   jp1 = j + 1

      dr(j) = t(left+j) - x
      dl(j) = x - t(left+1-j)
      saved = 0.d0

      do 30 i = 1, j
         term     = biatx(i) / (dr(i) + dl(jp1-i))
         biatx(i) = saved + dr(i) * term
         saved    = dl(jp1-i) * term
 30   continue

      biatx(jp1) = saved
      j          = jp1

      if (j .lt. jhigh) go to 20

 99   return

      end


c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      subroutine banfac(w,nroww,nrow,nbandl,nbandu,iflag)

      implicit double precision (a-h,o-z), integer (i-n)

      dimension w(nroww,nrow)

      iflag  = 1
      middle = nbandu + 1
      nrowm1 = nrow - 1

      if (nrowm1) 999,900,10

 10   if (nbandl .gt. 0) go to 30

      do 20 i = 1,nrowm1
         if (w(middle,i) .eq. 0.d0) go to 999
 20   continue

      go to 900

 30   if (nbandu .gt. 0) go to 60
      do 40 i = 1,nrowm1
         pivot = w(middle,i)
         if(pivot .eq. 0.d0) go to 999
         jmax = min0(nbandl, nrow - i)
         do 50 j = 1,jmax
            w(middle+j,i) = w(middle+j,i) / pivot
 50      continue
 40   continue

      return

 60   do 70 i = 1,nrowm1
         pivot = w(middle,i)
         if (pivot .eq. 0.d0) go to 999
         jmax = min0(nbandl,nrow - i)
         do 80 j = 1,jmax
            w(middle+j,i) = w(middle+j,i) / pivot
 80      continue

         kmax = min0(nbandu,nrow - i)

         do 90 k = 1,kmax
            ipk    = i + k
            midmk  = middle - k
            factor = w(midmk,ipk)
            do 100 j = 1,jmax
               w(midmk+j,ipk) = w(midmk+j,ipk) - w(middle+j,i)
     .              * factor
 100         continue
 90       continue
 70   continue

 900  if (w(middle,nrow) .ne. 0.d0) return
 999  iflag = 2

      return

      end


c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      subroutine banslv(w,nroww,nrow,nbandl,nbandu,b)

      implicit double precision (a-h,o-z), integer (i-n)

      dimension w(nroww,nrow),b(nrow)

      middle = nbandu + 1
      if (nrow .eq. 1) goto 99
      nrowm1 = nrow - 1
      if (nbandl .eq. 0) goto 30

      do 10 i = 1, nrowm1
         jmax = min0(nbandl, nrow - i)
         do 20 j = 1, jmax
            b(i+j) = b(i+j) - b(i) * w(middle+j,i)
 20      continue
 10   continue

 30   if (nbandu .gt. 0)  goto 50
      do 40 i = 1, nrow
         b(i) = b(i) / w(1,i)
 40   continue

      return

 50   do 60 i = nrow, 2, -1
         b(i) = b(i)/w(middle,i)
         jmax = min0(nbandu,i-1)
         do 70 j = 1,jmax
            b(i-j) = b(i-j) - b(i) * w(middle-j,i)
 70      continue
 60   continue

 99   b(1) = b(1) / w(middle,1)

      return

      end


c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      subroutine huntn(xx,n,kord,x,jlo)

      implicit double precision (a-h,o-z), integer (i-n)

      dimension xx(n)

c
c     works only for B-Splines (order n)
c

      max  = n - kord
      null = kord

      if (jlo.le.null.or.jlo.gt.max) then
         jlo = null
         jhi = max+1
         goto 30
      endif

      inc = 1

      if (x .ge. xx(jlo)) then
 10      jhi = jlo + inc
         if (jhi .gt. max) then
            jhi = max + 1
         else if (x .ge. xx(jhi)) then
            jlo = jhi
            inc = inc + inc
            goto 10
         endif
      else
         jhi = jlo
 20      jlo = jhi - inc
         if (jlo .le. null) then
            jlo = null
         else if (x .lt. xx(jlo)) then
            jhi = jlo
            inc = inc + inc
            goto 20
         endif
      endif

 30   if (jhi-jlo.eq.1) return

      jm = (jhi + jlo) / 2
      if (x .gt. xx(jm)) then
         jlo = jm
      else
         jhi = jm
      endif

      goto 30

      end


c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      FUNCTION ICHLENGTH(NAME, MSP)
      IMPLICIT NONE
c     **
c     ** Scalar arguments 
      CHARACTER*(*) NAME     
      INTEGER       MSP      
c     **
*     ..
*     Purpose
*     =======
*
*     ..
*     Arguments
*     =========
*
*     ..
*     Authors
*     =======
*     Freddy Fernandes Guimaraes
*
*     ..
*     Historic
*     ========
*     (10/03/2003) First version ICHLENGTH written by Freddy 
**
*
c     **
c     ** Parameters 
      REAL*8        ZERO, ONE, FOUR
      PARAMETER     (ZERO = +0.0D+0, ONE = +1.0D+0, FOUR =  +4.0D+0)
c     **
c     ** Local scalars 
      CHARACTER*1   LGTH
      INTEGER       I, J, ICHLENGTH
c     .. Start program
c     ..
      I = ZERO
 10   I = I + ONE
      LGTH = NAME(I:I)
      IF(LGTH.NE.' ')GOTO 10
      DO J=1,MSP,1
         I = I + ONE
         LGTH = NAME(I:I)
         IF(LGTH.NE.' ')GOTO 10
      ENDDO
c
      ICHLENGTH = I - (MSP + ONE)
c     .. End function
      RETURN
      END

c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      double precision function dlamch( cmach )
*
*  -- lapack auxiliary routine (version 2.0) --
*     univ. of tennessee, univ. of california berkeley, nag ltd.,
*     courant institute, argonne national lab, and rice university
*     october 31, 1992
*
*     .. scalar arguments ..
      character          cmach
*     ..
*
*  purpose
*  =======
*
*  dlamch determines double precision machine parameters.
*
*  arguments
*  =========
*
*  cmach   (input) character*1
*          specifies the value to be returned by dlamch:
*          = 'e' or 'e',   dlamch := eps
*          = 's' or 's ,   dlamch := sfmin
*          = 'b' or 'b',   dlamch := base
*          = 'p' or 'p',   dlamch := eps*base
*          = 'n' or 'n',   dlamch := t
*          = 'r' or 'r',   dlamch := rnd
*          = 'm' or 'm',   dlamch := emin
*          = 'u' or 'u',   dlamch := rmin
*          = 'l' or 'l',   dlamch := emax
*          = 'o' or 'o',   dlamch := rmax
*
*          where
*
*          eps   = relative machine precision
*          sfmin = safe minimum, such that 1/sfmin does not overflow
*          base  = base of the machine
*          prec  = eps*base
*          t     = number of (base) digits in the mantissa
*          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
*          emin  = minimum exponent before (gradual) underflow
*          rmin  = underflow threshold - base**(emin-1)
*          emax  = largest exponent before overflow
*          rmax  = overflow threshold  - (base**emax)*(1-eps)
*
* =====================================================================
*
*     .. parameters ..
      double precision   one, zero
      parameter          ( one = 1.0d+0, zero = 0.0d+0 )
*     ..
*     .. local scalars ..
      logical            first, lrnd
      integer            beta, imax, imin, it
      double precision   base, emax, emin, eps, prec, rmach, rmax, rmin,
     $                   rnd, sfmin, small, t
*     ..
*     .. external functions ..
      logical            lsame
      external           lsame
*     ..
*     .. external subroutines ..
      external           dlamc2
*     ..
*     .. save statement ..
      save               first, eps, sfmin, base, t, rnd, emin, rmin,
     $                   emax, rmax, prec
*     ..
*     .. data statements ..
      data               first / .true. /
*     ..
*     .. executable statements ..
*
      if( first ) then
         first = .false.
         call dlamc2( beta, it, lrnd, eps, imin, rmin, imax, rmax )
         base = beta
         t = it
         if( lrnd ) then
            rnd = one
            eps = ( base**( 1-it ) ) / 2
         else
            rnd = zero
            eps = base**( 1-it )
         end if
         prec = eps*base
         emin = imin
         emax = imax
         sfmin = rmin
         small = one / rmax
         if( small.ge.sfmin ) then
*
*           use small plus a bit, to avoid the possibility of rounding
*           causing overflow when computing  1/sfmin.
*
            sfmin = small*( one+eps )
         end if
      end if
*
      if( lsame( cmach, 'e' ) ) then
         rmach = eps
      else if( lsame( cmach, 's' ) ) then
         rmach = sfmin
      else if( lsame( cmach, 'b' ) ) then
         rmach = base
      else if( lsame( cmach, 'p' ) ) then
         rmach = prec
      else if( lsame( cmach, 'n' ) ) then
         rmach = t
      else if( lsame( cmach, 'r' ) ) then
         rmach = rnd
      else if( lsame( cmach, 'm' ) ) then
         rmach = emin
      else if( lsame( cmach, 'u' ) ) then
         rmach = rmin
      else if( lsame( cmach, 'l' ) ) then
         rmach = emax
      else if( lsame( cmach, 'o' ) ) then
         rmach = rmax
      end if
*
      dlamch = rmach
      return
*
*     end of dlamch
*
      end
*
************************************************************************
*
      subroutine dlamc1( beta, t, rnd, ieee1 )
*
*  -- lapack auxiliary routine (version 2.0) --
*     univ. of tennessee, univ. of california berkeley, nag ltd.,
*     courant institute, argonne national lab, and rice university
*     october 31, 1992
*
*     .. scalar arguments ..
      logical            ieee1, rnd
      integer            beta, t
*     ..
*
*  purpose
*  =======
*
*  dlamc1 determines the machine parameters given by beta, t, rnd, and
*  ieee1.
*
*  arguments
*  =========
*
*  beta    (output) integer
*          the base of the machine.
*
*  t       (output) integer
*          the number of ( beta ) digits in the mantissa.
*
*  rnd     (output) logical
*          specifies whether proper rounding  ( rnd = .true. )  or
*          chopping  ( rnd = .false. )  occurs in addition. this may not
*          be a reliable guide to the way in which the machine performs
*          its arithmetic.
*
*  ieee1   (output) logical
*          specifies whether rounding appears to be done in the ieee
*          'round to nearest' style.
*
*  further details
*  ===============
*
*  the routine is based on the routine  envron  by malcolm and
*  incorporates suggestions by gentleman and marovich. see
*
*     malcolm m. a. (1972) algorithms to reveal properties of
*        floating-point arithmetic. comms. of the acm, 15, 949-951.
*
*     gentleman w. m. and marovich s. b. (1974) more on algorithms
*        that reveal properties of floating point arithmetic units.
*        comms. of the acm, 17, 276-277.
*
* =====================================================================
*
*     .. local scalars ..
      logical            first, lieee1, lrnd
      integer            lbeta, lt
      double precision   a, b, c, f, one, qtr, savec, t1, t2
*     ..
*     .. external functions ..
      double precision   dlamc3
      external           dlamc3
*     ..
*     .. save statement ..
      save               first, lieee1, lbeta, lrnd, lt
*     ..
*     .. data statements ..
      data               first / .true. /
*     ..
*     .. executable statements ..
*
      if( first ) then
         first = .false.
         one = 1
*
*        lbeta,  lieee1,  lt and  lrnd  are the  local values  of  beta,
*        ieee1, t and rnd.
*
*        throughout this routine  we use the function  dlamc3  to ensure
*        that relevant values are  stored and not held in registers,  or
*        are not affected by optimizers.
*
*        compute  a = 2.0**m  with the  smallest positive integer m such
*        that
*
*           fl( a + 1.0 ) = a.
*
         a = 1
         c = 1
*
*+       while( c.eq.one )loop
   10    continue
         if( c.eq.one ) then
            a = 2*a
            c = dlamc3( a, one )
            c = dlamc3( c, -a )
            go to 10
         end if
*+       end while
*
*        now compute  b = 2.0**m  with the smallest positive integer m
*        such that
*
*           fl( a + b ) .gt. a.
*
         b = 1
         c = dlamc3( a, b )
*
*+       while( c.eq.a )loop
   20    continue
         if( c.eq.a ) then
            b = 2*b
            c = dlamc3( a, b )
            go to 20
         end if
*+       end while
*
*        now compute the base.  a and c  are neighbouring floating point
*        numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and so
*        their difference is beta. adding 0.25 to c is to ensure that it
*        is truncated to beta and not ( beta - 1 ).
*
         qtr = one / 4
         savec = c
         c = dlamc3( c, -a )
         lbeta = c + qtr
*
*        now determine whether rounding or chopping occurs,  by adding a
*        bit  less  than  beta/2  and a  bit  more  than  beta/2  to  a.
*
         b = lbeta
         f = dlamc3( b / 2, -b / 100 )
         c = dlamc3( f, a )
         if( c.eq.a ) then
            lrnd = .true.
         else
            lrnd = .false.
         end if
         f = dlamc3( b / 2, b / 100 )
         c = dlamc3( f, a )
         if( ( lrnd ) .and. ( c.eq.a ) )
     $      lrnd = .false.
*
*        try and decide whether rounding is done in the  ieee  'round to
*        nearest' style. b/2 is half a unit in the last place of the two
*        numbers a and savec. furthermore, a is even, i.e. has last  bit
*        zero, and savec is odd. thus adding b/2 to a should not  change
*        a, but adding b/2 to savec should change savec.
*
         t1 = dlamc3( b / 2, a )
         t2 = dlamc3( b / 2, savec )
         lieee1 = ( t1.eq.a ) .and. ( t2.gt.savec ) .and. lrnd
*
*        now find  the  mantissa, t.  it should  be the  integer part of
*        log to the base beta of a,  however it is safer to determine  t
*        by powering.  so we find t as the smallest positive integer for
*        which
*
*           fl( beta**t + 1.0 ) = 1.0.
*
         lt = 0
         a = 1
         c = 1
*
*+       while( c.eq.one )loop
   30    continue
         if( c.eq.one ) then
            lt = lt + 1
            a = a*lbeta
            c = dlamc3( a, one )
            c = dlamc3( c, -a )
            go to 30
         end if
*+       end while
*
      end if
*
      beta = lbeta
      t = lt
      rnd = lrnd
      ieee1 = lieee1
      return
*
*     end of dlamc1
*
      end
*
************************************************************************
*
      subroutine dlamc2( beta, t, rnd, eps, emin, rmin, emax, rmax )
*
*  -- lapack auxiliary routine (version 2.0) --
*     univ. of tennessee, univ. of california berkeley, nag ltd.,
*     courant institute, argonne national lab, and rice university
*     october 31, 1992
*
*     .. scalar arguments ..
      logical            rnd
      integer            beta, emax, emin, t
      double precision   eps, rmax, rmin
*     ..
*
*  purpose
*  =======
*
*  dlamc2 determines the machine parameters specified in its argument
*  list.
*
*  arguments
*  =========
*
*  beta    (output) integer
*          the base of the machine.
*
*  t       (output) integer
*          the number of ( beta ) digits in the mantissa.
*
*  rnd     (output) logical
*          specifies whether proper rounding  ( rnd = .true. )  or
*          chopping  ( rnd = .false. )  occurs in addition. this may not
*          be a reliable guide to the way in which the machine performs
*          its arithmetic.
*
*  eps     (output) double precision
*          the smallest positive number such that
*
*             fl( 1.0 - eps ) .lt. 1.0,
*
*          where fl denotes the computed value.
*
*  emin    (output) integer
*          the minimum exponent before (gradual) underflow occurs.
*
*  rmin    (output) double precision
*          the smallest normalized number for the machine, given by
*          base**( emin - 1 ), where  base  is the floating point value
*          of beta.
*
*  emax    (output) integer
*          the maximum exponent before overflow occurs.
*
*  rmax    (output) double precision
*          the largest positive number for the machine, given by
*          base**emax * ( 1 - eps ), where  base  is the floating point
*          value of beta.
*
*  further details
*  ===============
*
*  the computation of  eps  is based on a routine paranoia by
*  w. kahan of the university of california at berkeley.
*
* =====================================================================
*
*     .. local scalars ..
      logical            first, ieee, iwarn, lieee1, lrnd
      integer            gnmin, gpmin, i, lbeta, lemax, lemin, lt,
     $                   ngnmin, ngpmin
      double precision   a, b, c, half, leps, lrmax, lrmin, one, rbase,
     $                   sixth, small, third, two, zero
*     ..
*     .. external functions ..
      double precision   dlamc3
      external           dlamc3
*     ..
*     .. external subroutines ..
      external           dlamc1, dlamc4, dlamc5
*     ..
*     .. intrinsic functions ..
      intrinsic          abs, max, min
*     ..
*     .. save statement ..
      save               first, iwarn, lbeta, lemax, lemin, leps, lrmax,
     $                   lrmin, lt
*     ..
*     .. data statements ..
      data               first / .true. / , iwarn / .false. /
*     ..
*     .. executable statements ..
*
      if( first ) then
         first = .false.
         zero = 0
         one = 1
         two = 2
*
*        lbeta, lt, lrnd, leps, lemin and lrmin  are the local values of
*        beta, t, rnd, eps, emin and rmin.
*
*        throughout this routine  we use the function  dlamc3  to ensure
*        that relevant values are stored  and not held in registers,  or
*        are not affected by optimizers.
*
*        dlamc1 returns the parameters  lbeta, lt, lrnd and lieee1.
*
         call dlamc1( lbeta, lt, lrnd, lieee1 )
*
*        start to find eps.
*
         b = lbeta
         a = b**( -lt )
         leps = a
*
*        try some tricks to see whether or not this is the correct  eps.
*
         b = two / 3
         half = one / 2
         sixth = dlamc3( b, -half )
         third = dlamc3( sixth, sixth )
         b = dlamc3( third, -half )
         b = dlamc3( b, sixth )
         b = abs( b )
         if( b.lt.leps )
     $      b = leps
*
         leps = 1
*
*+       while( ( leps.gt.b ).and.( b.gt.zero ) )loop
   10    continue
         if( ( leps.gt.b ) .and. ( b.gt.zero ) ) then
            leps = b
            c = dlamc3( half*leps, ( two**5 )*( leps**2 ) )
            c = dlamc3( half, -c )
            b = dlamc3( half, c )
            c = dlamc3( half, -b )
            b = dlamc3( half, c )
            go to 10
         end if
*+       end while
*
         if( a.lt.leps )
     $      leps = a
*
*        computation of eps complete.
*
*        now find  emin.  let a = + or - 1, and + or - (1 + base**(-3)).
*        keep dividing  a by beta until (gradual) underflow occurs. this
*        is detected when we cannot recover the previous a.
*
         rbase = one / lbeta
         small = one
         do 20 i = 1, 3
            small = dlamc3( small*rbase, zero )
   20    continue
         a = dlamc3( one, small )
         call dlamc4( ngpmin, one, lbeta )
         call dlamc4( ngnmin, -one, lbeta )
         call dlamc4( gpmin, a, lbeta )
         call dlamc4( gnmin, -a, lbeta )
         ieee = .false.
*
         if( ( ngpmin.eq.ngnmin ) .and. ( gpmin.eq.gnmin ) ) then
            if( ngpmin.eq.gpmin ) then
               lemin = ngpmin
*            ( non twos-complement machines, no gradual underflow;
*              e.g.,  vax )
            else if( ( gpmin-ngpmin ).eq.3 ) then
               lemin = ngpmin - 1 + lt
               ieee = .true.
*            ( non twos-complement machines, with gradual underflow;
*              e.g., ieee standard followers )
            else
               lemin = min( ngpmin, gpmin )
*            ( a guess; no known machine )
               iwarn = .true.
            end if
*
         else if( ( ngpmin.eq.gpmin ) .and. ( ngnmin.eq.gnmin ) ) then
            if( abs( ngpmin-ngnmin ).eq.1 ) then
               lemin = max( ngpmin, ngnmin )
*            ( twos-complement machines, no gradual underflow;
*              e.g., cyber 205 )
            else
               lemin = min( ngpmin, ngnmin )
*            ( a guess; no known machine )
               iwarn = .true.
            end if
*
         else if( ( abs( ngpmin-ngnmin ).eq.1 ) .and.
     $            ( gpmin.eq.gnmin ) ) then
            if( ( gpmin-min( ngpmin, ngnmin ) ).eq.3 ) then
               lemin = max( ngpmin, ngnmin ) - 1 + lt
*            ( twos-complement machines with gradual underflow;
*              no known machine )
            else
               lemin = min( ngpmin, ngnmin )
*            ( a guess; no known machine )
               iwarn = .true.
            end if
*
         else
            lemin = min( ngpmin, ngnmin, gpmin, gnmin )
*         ( a guess; no known machine )
            iwarn = .true.
         end if
***
* comment out this if block if emin is ok
         if( iwarn ) then
            first = .true.
            write( 6, fmt = 9999 )lemin
         end if
***
*
*        assume ieee arithmetic if we found denormalised  numbers above,
*        or if arithmetic seems to round in the  ieee style,  determined
*        in routine dlamc1. a true ieee machine should have both  things
*        true; however, faulty machines may have one or the other.
*
         ieee = ieee .or. lieee1
*
*        compute  rmin by successive division by  beta. we could compute
*        rmin as base**( emin - 1 ),  but some machines underflow during
*        this computation.
*
         lrmin = 1
         do 30 i = 1, 1 - lemin
            lrmin = dlamc3( lrmin*rbase, zero )
   30    continue
*
*        finally, call dlamc5 to compute emax and rmax.
*
         call dlamc5( lbeta, lt, lemin, ieee, lemax, lrmax )
      end if
*
      beta = lbeta
      t = lt
      rnd = lrnd
      eps = leps
      emin = lemin
      rmin = lrmin
      emax = lemax
      rmax = lrmax
*
      return
*
 9999 format( / / ' warning. the value emin may be incorrect:-',
     $      '  emin = ', i8, /
     $      ' if, after inspection, the value emin looks',
     $      ' acceptable please comment out ',
     $      / ' the if block as marked within the code of routine',
     $      ' dlamc2,', / ' otherwise supply emin explicitly.', / )
*
*     end of dlamc2
*
      end
*
************************************************************************
*
      double precision function dlamc3( a, b )
*
*  -- lapack auxiliary routine (version 2.0) --
*     univ. of tennessee, univ. of california berkeley, nag ltd.,
*     courant institute, argonne national lab, and rice university
*     october 31, 1992
*
*     .. scalar arguments ..
      double precision   a, b
*     ..
*
*  purpose
*  =======
*
*  dlamc3  is intended to force  a  and  b  to be stored prior to doing
*  the addition of  a  and  b ,  for use in situations where optimizers
*  might hold one of these in a register.
*
*  arguments
*  =========
*
*  a, b    (input) double precision
*          the values a and b.
*
* =====================================================================
*
*     .. executable statements ..
*
      dlamc3 = a + b
*
      return
*
*     end of dlamc3
*
      end
*
************************************************************************
*
      subroutine dlamc4( emin, start, base )
*
*  -- lapack auxiliary routine (version 2.0) --
*     univ. of tennessee, univ. of california berkeley, nag ltd.,
*     courant institute, argonne national lab, and rice university
*     october 31, 1992
*
*     .. scalar arguments ..
      integer            base, emin
      double precision   start
*     ..
*
*  purpose
*  =======
*
*  dlamc4 is a service routine for dlamc2.
*
*  arguments
*  =========
*
*  emin    (output) emin
*          the minimum exponent before (gradual) underflow, computed by
*          setting a = start and dividing by base until the previous a
*          can not be recovered.
*
*  start   (input) double precision
*          the starting point for determining emin.
*
*  base    (input) integer
*          the base of the machine.
*
* =====================================================================
*
*     .. local scalars ..
      integer            i
      double precision   a, b1, b2, c1, c2, d1, d2, one, rbase, zero
*     ..
*     .. external functions ..
      double precision   dlamc3
      external           dlamc3
*     ..
*     .. executable statements ..
*
      a = start
      one = 1
      rbase = one / base
      zero = 0
      emin = 1
      b1 = dlamc3( a*rbase, zero )
      c1 = a
      c2 = a
      d1 = a
      d2 = a
*+    while( ( c1.eq.a ).and.( c2.eq.a ).and.
*    $       ( d1.eq.a ).and.( d2.eq.a )      )loop
   10 continue
      if( ( c1.eq.a ) .and. ( c2.eq.a ) .and. ( d1.eq.a ) .and.
     $    ( d2.eq.a ) ) then
         emin = emin - 1
         a = b1
         b1 = dlamc3( a / base, zero )
         c1 = dlamc3( b1*base, zero )
         d1 = zero
         do 20 i = 1, base
            d1 = d1 + b1
   20    continue
         b2 = dlamc3( a*rbase, zero )
         c2 = dlamc3( b2 / rbase, zero )
         d2 = zero
         do 30 i = 1, base
            d2 = d2 + b2
   30    continue
         go to 10
      end if
*+    end while
*
      return
*
*     end of dlamc4
*
      end
*
************************************************************************
*
      subroutine dlamc5( beta, p, emin, ieee, emax, rmax )
*
*  -- lapack auxiliary routine (version 2.0) --
*     univ. of tennessee, univ. of california berkeley, nag ltd.,
*     courant institute, argonne national lab, and rice university
*     october 31, 1992
*
*     .. scalar arguments ..
      logical            ieee
      integer            beta, emax, emin, p
      double precision   rmax
*     ..
*
*  purpose
*  =======
*
*  dlamc5 attempts to compute rmax, the largest machine floating-point
*  number, without overflow.  it assumes that emax + abs(emin) sum
*  approximately to a power of 2.  it will fail on machines where this
*  assumption does not hold, for example, the cyber 205 (emin = -28625,
*  emax = 28718).  it will also fail if the value supplied for emin is
*  too large (i.e. too close to zero), probably with overflow.
*
*  arguments
*  =========
*
*  beta    (input) integer
*          the base of floating-point arithmetic.
*
*  p       (input) integer
*          the number of base beta digits in the mantissa of a
*          floating-point value.
*
*  emin    (input) integer
*          the minimum exponent before (gradual) underflow.
*
*  ieee    (input) logical
*          a logical flag specifying whether or not the arithmetic
*          system is thought to comply with the ieee standard.
*
*  emax    (output) integer
*          the largest exponent before overflow
*
*  rmax    (output) double precision
*          the largest machine floating-point number.
*
* =====================================================================
*
*     .. parameters ..
      double precision   zero, one
      parameter          ( zero = 0.0d0, one = 1.0d0 )
*     ..
*     .. local scalars ..
      integer            exbits, expsum, i, lexp, nbits, try, uexp
      double precision   oldy, recbas, y, z
*     ..
*     .. external functions ..
      double precision   dlamc3
      external           dlamc3
*     ..
*     .. intrinsic functions ..
      intrinsic          mod
*     ..
*     .. executable statements ..
*
*     first compute lexp and uexp, two powers of 2 that bound
*     abs(emin). we then assume that emax + abs(emin) will sum
*     approximately to the bound that is closest to abs(emin).
*     (emax is the exponent of the required number rmax).
*
      lexp = 1
      exbits = 1
   10 continue
      try = lexp*2
      if( try.le.( -emin ) ) then
         lexp = try
         exbits = exbits + 1
         go to 10
      end if
      if( lexp.eq.-emin ) then
         uexp = lexp
      else
         uexp = try
         exbits = exbits + 1
      end if
*
*     now -lexp is less than or equal to emin, and -uexp is greater
*     than or equal to emin. exbits is the number of bits needed to
*     store the exponent.
*
      if( ( uexp+emin ).gt.( -lexp-emin ) ) then
         expsum = 2*lexp
      else
         expsum = 2*uexp
      end if
*
*     expsum is the exponent range, approximately equal to
*     emax - emin + 1 .
*
      emax = expsum + emin - 1
      nbits = 1 + exbits + p
*
*     nbits is the total number of bits needed to store a
*     floating-point number.
*
      if( ( mod( nbits, 2 ).eq.1 ) .and. ( beta.eq.2 ) ) then
*
*        either there are an odd number of bits used to store a
*        floating-point number, which is unlikely, or some bits are
*        not used in the representation of numbers, which is possible,
*        (e.g. cray machines) or the mantissa has an implicit bit,
*        (e.g. ieee machines, dec vax machines), which is perhaps the
*        most likely. we have to assume the last alternative.
*        if this is true, then we need to reduce emax by one because
*        there must be some way of representing zero in an implicit-bit
*        system. on machines like cray, we are reducing emax by one
*        unnecessarily.
*
         emax = emax - 1
      end if
*
      if( ieee ) then
*
*        assume we are on an ieee machine which reserves one exponent
*        for infinity and nan.
*
         emax = emax - 1
      end if
*
*     now create rmax, the largest machine number, which should
*     be equal to (1.0 - beta**(-p)) * beta**emax .
*
*     first compute 1.0 - beta**(-p), being careful that the
*     result is less than 1.0 .
*
      recbas = one / beta
      z = beta - one
      y = zero
      do 20 i = 1, p
         z = z*recbas
         if( y.lt.one )
     $      oldy = y
         y = dlamc3( y, z )
   20 continue
      if( y.ge.one )
     $   y = oldy
*
*     now multiply by beta**emax to get rmax.
*
      do 30 i = 1, emax
         y = dlamc3( y*beta, zero )
   30 continue
*
      rmax = y
      return
*
*     end of dlamc5
*
      end


      logical          function lsame( ca, cb )
*
*  -- lapack auxiliary routine (version 3.0) --
*     univ. of tennessee, univ. of california berkeley, nag ltd.,
*     courant institute, argonne national lab, and rice university
*     september 30, 1994
*
*     .. scalar arguments ..
      character          ca, cb
*     ..
*
*  purpose
*  =======
*
*  lsame returns .true. if ca is the same letter as cb regardless of
*  case.
*
*  arguments
*  =========
*
*  ca      (input) character*1
*  cb      (input) character*1
*          ca and cb specify the single characters to be compared.
*
* =====================================================================
*
*     .. intrinsic functions ..
      intrinsic          ichar
*     ..
*     .. local scalars ..
      integer            inta, intb, zcode
*     ..
*     .. executable statements ..
*
*     test if the characters are equal
*
      lsame = ca.eq.cb
      if( lsame )
     $   return
*
*     now test for equivalence if both characters are alphabetic.
*
      zcode = ichar( 'z' )
*
*     use 'z' rather than 'a' so that ascii can be detected on prime
*     machines, on which ichar returns a value with bit 8 set.
*     ichar('a') on prime machines returns 193 which is the same as
*     ichar('a') on an ebcdic machine.
*
      inta = ichar( ca )
      intb = ichar( cb )
*
      if( zcode.eq.90 .or. zcode.eq.122 ) then
*
*        ascii is assumed - zcode is the ascii code of either lower or
*        upper case 'z'.
*
         if( inta.ge.97 .and. inta.le.122 ) inta = inta - 32
         if( intb.ge.97 .and. intb.le.122 ) intb = intb - 32
*
      else if( zcode.eq.233 .or. zcode.eq.169 ) then
*
*        ebcdic is assumed - zcode is the ebcdic code of either lower or
*        upper case 'z'.
*
         if( inta.ge.129 .and. inta.le.137 .or.
     $       inta.ge.145 .and. inta.le.153 .or.
     $       inta.ge.162 .and. inta.le.169 ) inta = inta + 64
         if( intb.ge.129 .and. intb.le.137 .or.
     $       intb.ge.145 .and. intb.le.153 .or.
     $       intb.ge.162 .and. intb.le.169 ) intb = intb + 64
*
      else if( zcode.eq.218 .or. zcode.eq.250 ) then
*
*        ascii is assumed, on prime machines - zcode is the ascii code
*        plus 128 of either lower or upper case 'z'.
*
         if( inta.ge.225 .and. inta.le.250 ) inta = inta - 32
         if( intb.ge.225 .and. intb.le.250 ) intb = intb - 32
      end if
      lsame = inta.eq.intb
*
*     return
*
*     end of lsame
*
      end
