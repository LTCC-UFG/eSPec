      SUBROUTINE LANCZSG(CHANGE, DIM, IL, IU, INFO, NREORT, LMTREORT, N, 
     &     NSEED, MXDCT, ABSTOL, IWK, MAXINT, NP, AL, BT, EIGVL, SHM, 
     &     V0, VI, VPOT, WK1, WK2, VAR, LNZVC, EIGVC)
      IMPLICIT NONE
c     **
c     ** Scalar arguments 
      LOGICAL       CHANGE
      CHARACTER*(*) DIM     
      INTEGER       IL, IU, INFO, NREORT, LMTREORT, N, NSEED, MXDCT 
      INTEGER       MAXINT
      REAL*8        ABSTOL
c     **
c     ** Array arguments
cdel      LOGICAL
      INTEGER       IWK(*), NP(*)
      REAL*8        AL(*), BT(*),  EIGVL(*), SHM(*), V0(*), VI(*)
      REAL*8        VPOT(*), WK1(*), WK2(*), VAR(*)
      REAL*8        LNZVC(MXDCT,*), EIGVC(MXDCT,*)
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
*     (18/03/2003) First version LANCZS written by Freddy
*
c     **
c     ** Parameters 
      REAL*8        ZERO, ONE, TWO
      PARAMETER     (ZERO = +0.0D+0, ONE = +1.0D+0, TWO = +2.0D+0)
c     **
c     ** Local scalars 
cdel      LOGICAL
cdel      CHARACTER*1
      INTEGER       I, J, K, M, KK
      REAL*8        VL, VU, SUM, FATN
c     **
c     ** Local arrays 
cdel      LOGICAL
cdel      CHARACTER*1
c      INTEGER       IFAIL(IU-IL+1)
       REAL*8        VIAUX(MXDCT)
D      REAL*8        EIGVC(MXDCT,IU-IL+1), VIAUX(MXDCT)
c     **
c     ** External functions 
cdel      LOGICAL
cdel      CHARACTER*2   CHAUX
cdel      INTEGER       
      REAL*8        ECNORM, RAN
c     **
c     ** External subroutines 
      EXTERNAL      LNZ, REORT, AU, DSTEVX
c     **
c     ** Intrinsic functions 
cdel      INTRINSIC     
c     .. Start program 
c     ..
c     .. init and test parameters 
      IF(NREORT.GT.LMTREORT)THEN
         WRITE(*,*)'Error NREORT > LMTREORT'
         WRITE(*,*)'NREORT=',NREORT,'  LMTREORT=',LMTREORT 
         NREORT = LMTREORT
         WRITE(*,*)'Changing NREORT to correct value!'
         WRITE(*,*)'NREORT=',NREORT,'  LMTREORT=',LMTREORT
      ENDIF
c     
      NSEED = 1
      NREORT = 40
      DO I=1,N,1
         VI(I) = RAN(NSEED) 
      ENDDO
c      
      DO K=1,IU-IL+1,1
         DO I=1,N,1
            VI(I) = RAN(NSEED) 
         ENDDO
c
         KK = ZERO
 100     CONTINUE
         KK = KK + ONE
            DO I=1,N,1
               VIAUX(I) = VI(I)
            ENDDO
c     
            CALL LNZ(CHANGE, DIM, N, NREORT, NREORT, LMTREORT, MXDCT, 
     &           NP, SHM, VPOT, V0, VI, WK1, AL, BT, VAR, LNZVC)
c     
            CALL DSTEVX('V', 'I', NREORT, AL, BT, VL, VU, IL, IU, 
     &           ABSTOL, M, EIGVL, EIGVC, MXDCT, WK2, IWK, INFO)
c
            DO I=1,N,1
               VI(I) = ZERO
            ENDDO
c
            SUM = ZERO
            DO J=1,N,1
               DO I=1,NREORT,1
                  VI(J) = VI(J) + EIGVC(I,K)*LNZVC(J,I)
               ENDDO
               SUM = SUM + ABS(VI(I) - VIAUX(I))
            ENDDO
c
            IF(KK.LT.10)THEN
               WRITE(*,1031)KK,EIGVL(K),SUM
            ELSEIF(KK.LT.100)THEN
               WRITE(*,1032)KK,EIGVL(K),SUM
            ELSEIF(KK.LT.1000)THEN
               WRITE(*,1033)KK,EIGVL(K),SUM
            ELSE
               WRITE(*,*)'<<<>>> Too many interactions <<<>>>'
               WRITE(*,*)'Stopping program!'
               STOP
            ENDIF
c
         IF(SUM.GT.ABSTOL .AND. KK.LE.MAXINT)GOTO 100
         IF(KK.GE.MAXINT)THEN
            WRITE(*,*)'<<<>>> Eigenvector didn´t achieve the ',
     &           'convergence <<<>>>'
            WRITE(*,*)'Stopping eSPec!'
            STOP
         ENDIF
c     
         FATN = ECNORM(N, VI)
         DO J=1,N,1
            EIGVC(J,K) = VI(J)/FATN
         ENDDO
      ENDDO
c
      DO J=1,IU-IL+1,1
         DO I=1,N,1
            LNZVC(I,J) = EIGVC(I,J)
         ENDDO
      ENDDO
c     ..
 1031 FORMAT(5X,'Interaction',1X,I1,';',1X,'Energy =',E16.8,'/a.u.; ',
     &     'WF converg. =',E16.8,';') 
 1032 FORMAT(5X,'Interaction',1X,I2,';',1X,'Energy =',E16.8,'/a.u.; ',
     &     'WF converg. =',E16.8,';') 
 1033 FORMAT(5X,'Interaction',1X,I3,';',1X,'Energy =',E16.8,'/a.u.; ',
     &     'WF converg. =',E16.8,';')
      RETURN
      END
