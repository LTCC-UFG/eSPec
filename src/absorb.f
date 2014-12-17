      SUBROUTINE ABSORBINGBC(DIM, TPABSOR, N, NP, AL, AR, SH, XI, XF, 
     &     VOI, VABC, ALPHA)
      IMPLICIT NONE
c     **
c     ** Scalar arguments 
cdel      LOGICAL       
      CHARACTER*(*) DIM, TPABSOR
      INTEGER       N
      REAL*8        ALPHA
c     **
c     ** Array arguments
cdel      LOGICAL       
cdel      CHARACTER*(*) 
      INTEGER       NP(*)
      REAL*8        XI(*), XF(*), SH(*), AL(*), AR(*),  VOI(*), VABC(*)
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
*     () First version written by 
*
c     **
c     ** Parameters 
      REAL*8        ZERO, ONE, TWO
      PARAMETER     (ZERO = +0.0D+0, ONE = +1.0D+0, TWO = +2.0D+0)
c     **
c     ** Local scalars 
cdel      LOGICAL       
cdel      CHARACTER*1   
cdel      INTEGER       
cdel      REAL*8        
c     **
c     ** External functions 
cdel      LOGICAL       
cdel      CHARACTER*(*) 
      INTEGER       I, J, L, NPAUX, K
      REAL*8        X1, X2!, X3
c     **
c     ** External subroutines 
cdel      EXTERNAL      
c     **
c     ** Intrinsic functions 
      INTRINSIC     DEXP, COSH 
c     .. Start program       
d      write(*,*)'al(1),ar(1),voi(1)',al(1),ar(1),voi(1)
d      write(*,*)'al(2),ar(2),voi(2)',al(2),ar(2),voi(2)
d      write(*,*)'xf(1),xf(2)',xf(1),xf(2)
d      read(*,*)
      IF(DIM(1:3).EQ.'.1D')THEN
         IF(TPABSOR(1:8).EQ.'.SMOOTHW')THEN
            ALPHA = 1 - 0.1*(XF(1) - AR(1))
            DO I=1,N,1
               X1 = XI(1) + (I - ONE)*SH(1) 
               IF(X1.LT.AL(1))THEN
                  VABC(I) = DEXP(VOI(1)*(X1 - AL(1))**VOI(2)
     &                 /(X1 + SH(1) - XI(1)))
               ELSEIF(X1.GT.AR(1))THEN
                  VABC(I) = DEXP(-VOI(1)*(X1 - AR(1))
     &                 /(XF(1) + SH(1) - X1)**VOI(2))
               ELSE
                  VABC(I) = ONE
               ENDIF
c               write(1,*)X1, VABC(I), AL(1), AR(1)
c               read(*,*)
            ENDDO
d            write(*,*)'VABC(1)',VABC(1),VABC(2)
d            read(*,*)
         ELSEIF(TPABSOR(1:5).EQ.'.SMRS')THEN
            K = 0
            DO I=1,N,1
               X1 = XI(1) + (I - ONE)*SH(1) 
               IF(X1.LT.AL(1))THEN
d                  VABC(I) = DEXP(VOI(1)*(X1 - AL(1))
d     &                 /(X1 + SH(1) - XI(1)))
                  VABC(I) = ONE - VOI(1)*(X1 - AL(1))**2
d               ELSEIF(X1.GT.AR(1) .AND. VABC(I-1).LE.0.99999)THEN
d                  VABC(I) = 0.99999
d               ELSEIF(X1.GT.AR(1) .AND. I.GE.N-9)THEN
d                  VABC(I) = DEXP(-VOI(1)*(X1 - AR(1))
d     &                 /(XF(1) + K**2*SH(1) - X1))
d                  K = K + 1
d               ELSEIF(X1.GT.AR(1)+(XF(1)-AR(1))/2)THEN
d                  print *, 'oi'
d                  VABC(I) = DEXP(-VOI(1)*(X1 - AR(1)-(XF(1)-AR(1))/2)
d     &                 /(XF(1) + SH(1) - X1))
d                  print *, 'oi',x1,VABC(I),(X1 - AR(1)-(XF(1)-AR(1))/2)
d               ELSEIF(I.GT.N-3)THEN
d                  VABC(I) = ONE
               ELSEIF(X1.GT.AR(1))THEN
                  VABC(I) = ONE - VOI(1)*(X1 - AR(1))**2
d                  VABC(I) = DEXP(-VOI(1)*(X1 - AR(1))
d     &                 /(XF(1) + SH(1) - X1))
                  print *, 'oi2',x1,VABC(I),(X1 - AR(1))
               ELSE
                  VABC(I) = ONE
               ENDIF
c               write(1,*)X1, VABC(I), AL(1), AR(1)
c               read(*,*)
            ENDDO
d            write(*,*)'VABC(1)',VABC(1),VABC(2)
d            read(*,*)
         ELSEIF(TPABSOR(1:14).EQ.'.VOPTIC_LINEAR')THEN
c            write(*,*)'oi543',AL(1),AR(1)
            DO I=1,N,1
               X1 = XI(1) + (I - ONE)*SH(1) 
               IF(X1.LT.AL(1))THEN
c                  write(*,*)'AL(1)'
                  VABC(I) = -VOI(1)*(X1 - AL(1))/(X1 + SH(1) - XI(1))
               ELSEIF(X1.GT.AR(1))THEN
c                  write(*,*)'AR(1)'
                  VABC(I) = VOI(1)*(X1 - AR(1))/(XF(1) + SH(1) - X1)
               ELSE
                  
                  VABC(I) = ZERO
               ENDIF
c            write(*,*)'oiI',I,X1,VABC(I)
            ENDDO
         ELSEIF(TPABSOR(1:15).EQ.'.VOPTIC_KOSLOFF')THEN
            DO I=1,N,1
               X1 = XI(1) + (I - ONE)*SH(1) 
               IF(X1.LT.AL(1))THEN
                  VABC(I) = VOI(1)/COSH(XI(1) - AL(1)/VOI(2))
               ELSEIF(X1.GT.AR(1))THEN
                  VABC(I) = VOI(1)/COSH(XF(1) - AL(1)/VOI(2))
               ELSE
                  VABC(I) = ZERO
               ENDIF
            ENDDO
         ELSEIF(TPABSOR(1:13).EQ.'.VOPTIC_POLIN')THEN
            DO I=1,N,1
               X1 = XI(1) + (I - ONE)*SH(1) 
               IF(X1.LT.AL(1))THEN
                  VABC(I) = VOI(1)*(X1 - AL(1))**VOI(2)
               ELSEIF(X1.GT.AR(1))THEN
                  VABC(I) = VOI(1)*(X1 - AR(1))**VOI(2)
               ELSE
                  VABC(I) = ZERO
               ENDIF
            ENDDO  
         ELSE
            WRITE(*,*)
            STOP
         ENDIF
      ELSEIF(DIM(1:3).EQ.'.2D')THEN
         IF(TPABSOR(1:8).EQ.'.SMOOTHW')THEN
d            write(*,*)'oi'
d            read(*,*)
            DO J=1,NP(2),1
               X2 = XI(2) + (J - ONE)*SH(2)
               NPAUX = (J - ONE)*NP(1)
               DO I=1,NP(1),1
                  X1 = XI(1) + (I - ONE)*SH(1)
                  L = I + NPAUX
c                  IF(J.GE.(NP(2) - 5) .OR. I.GE.(NP(1) - 5)
c     &                 )THEN
c                     VABC(L) = ZERO
                  IF(X1.LT.AL(1))THEN
                     VABC(L) = DEXP(VOI(1)*(X1 - AL(1))
     &                    /(X1 + SH(1) - XI(1)))
                  ELSEIF(X2.LT.AL(2))THEN
                     VABC(L) = DEXP(VOI(2)*(X2 - AL(2))
     &                    /(X2 + SH(2) - XI(2)))
                  ELSEIF(X1.GT.AR(1))THEN
                     VABC(L) = DEXP(-VOI(1)*(X1 - AR(1))
     &                    /(XF(1) + SH(1) - X1))
d                  ELSEIF(X2.GT.AR(2) .AND. X1.LT.4.5)THEN
                  ELSEIF(X2.GT.AR(2))THEN
                     VABC(L) = DEXP(-VOI(2)*(X2 - AR(2))
     &                    /(XF(2) + SH(2) - X2))
                  ELSE
                     VABC(L) = ONE
                  ENDIF
D                  if(vabc(l).gt.one)then
D                     write(*,*)'x1,x2',x1,x2,ONE
D                     write(*,*)'voi(1),voi(2)',voi(1),voi(2)
D                     write(*,*)'al(1),ar(1)',al(1),ar(1)
D                     write(*,*)'al(2),ar(2)',al(2),ar(2)
D                     write(*,*)'vabc(l)',vabc(l)
D                     stop
D                  endif
               ENDDO
            ENDDO
         ELSEIF(TPABSOR(1:7).EQ.'.VOPTIC')THEN
            WRITE(*,*)'<<<>>> Absorbing bondary conditions .VOPTIC ' 
            WRITE(*,*)'          is not implemente in the code <<<>>>' 
         STOP
            STOP
         ENDIF
      ELSEIF(DIM(1:3).EQ.'.3D')THEN
         WRITE(*,*)'<<<>>> Absorbing bondary conditions is not ' 
         WRITE(*,*)'      implemente for the 3-dimensional code <<<>>>' 
         STOP
      ELSE
         WRITE(*,*)'<<<>>> Dimension error <<<>>>' 
         STOP
      ENDIF
c     ..
      RETURN
      END


d
d                  IF(X1.LT.AL(1))THEN
d                     VABC(L) = DEXP(VOI(1)*(X1 - AL(1))
d     &                    /(X1 - XI(1)))
d                  ELSEIF(X1.GT.AR(1) .AND. X2.GT.AR(2))THEN
d                     VABC(L) = DEXP(-VOI(1)*(X1 - AR(1))
d     &                    /(XF(1) + SH(1) - X1))
d     &                    * DEXP(-VOI(2)*(X2 - AR(2))
d     &                    /(XF(2) + SH(2) - X2))
d                  ELSEIF(X1.LT.AL(1) .AND. X2.GT.AR(2))THEN
d                     VABC(L) = DEXP(VOI(1)*(X1 - AL(1))
d     &                    /(X1 - XI(1)))
d     &                    * DEXP(-VOI(2)*(X2 - AR(2))
d     &                    /(XF(2) + SH(2) - X2))
d                  ELSEIF(X1.GT.AR(1) .AND. X2.LT.AL(2))THEN
d                     VABC(L) = DEXP(-VOI(1)*(X1 - AR(1))
d     &                    /(XF(1) + SH(1) - X1))
d     &                    * DEXP(VOI(2)*(X2 - AL(2))
d     &                    /(X2 - XI(2)))
d                  ELSEIF(X1.LT.AL(1))THEN
d                     VABC(L) = DEXP(VOI(1)*(X1 - AL(1))
d     &                    /(X1 - XI(1)))
d                  ELSEIF(X1.GT.AR(1))THEN
d                     VABC(L) = DEXP(-VOI(1)*(X1 - AR(1))
d     &                    /(XF(1) + SH(1) - X1))
d                  ELSEIF(X2.LT.AL(2))THEN
d                     VABC(L) = DEXP(VOI(2)*(X2 - AL(2))
d     &                    /(X2 - XI(2)))
d                  ELSEIF(X2.GT.AR(2))THEN
d                     VABC(L) = DEXP(-VOI(2)*(X2 - AR(2))
d     &                    /(XF(2) + SH(2) - X2))
d                  ELSE
d                     VABC(L) = ONE
d                  ENDIF
