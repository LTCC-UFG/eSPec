      SUBROUTINE LANCZS(CHANGE, DIM, IL, IU, INFO, NREORT, LMTREORT, N, 
     &     NSEED, MXDCT, ABSTOL, IWK, NP, AL, BT, EIGVL, SHM, V0, VI,
     &     VPOT, WK1, WK2, VAR, LNZVC, EIGVC)
      IMPLICIT NONE
c     **
c     ** Scalar arguments 
      LOGICAL       CHANGE
      CHARACTER*(*) DIM     
      INTEGER       IL, IU, INFO, NREORT, LMTREORT, N, NSEED, MXDCT 
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
      CHARACTER*2   CHNC
      INTEGER       I, IT, J, K, M, NC, NSEEDAUX
      REAL*8        ALPH, BT0, SUM1, VL, VU
c     **
c     ** Local arrays 
cdel      LOGICAL
cdel      CHARACTER*1
c      INTEGER       IFAIL(IU-IL+1)
D      REAL*8        EIGVC(MXDCT,IU-IL+1)
c     **
c     ** External functions 
cdel      LOGICAL
      CHARACTER*2   CHAUX
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
      IF(NREORT.GT.LMTREORT .OR. NREORT.GT.N)THEN
         WRITE(*,*)'Error NREORT > LMTREORT or NREORT > N'
         WRITE(*,*)'NREORT=',NREORT,'  LMTREORT=',LMTREORT,'  N=',N  
         IF(CHANGE)THEN
            IF(NREORT.GT.LMTREORT)THEN
               NREORT = LMTREORT
            ELSEIF(NREORT.GT.N)THEN
               NREORT = N
            ENDIF
            WRITE(*,*)'Modifing NREORT to correct value!'
            WRITE(*,*)'NREORT=',NREORT,'  LMTREORT=',LMTREORT,'  N=',N
         ENDIF       
      ENDIF
           
D      write(*,*)'*****>',n,nreort
c 
      NSEED = 1
      NSEEDAUX = NSEED
      DO I=1,N,1
         VI(I) = RAN(NSEEDAUX) 
      ENDDO
c     
      CALL LNZ(CHANGE, DIM, N, N, NREORT, LMTREORT, MXDCT, NP, SHM, 
     &     VPOT, V0, VI, WK1, AL, BT, VAR, LNZVC)
c
      CALL DSTEVX('V', 'I', N, AL, BT, VL, VU, IL, IU, ABSTOL,
     &        M, EIGVL, EIGVC, MXDCT, WK2, IWK, INFO)
c     
      NC = IL - TWO
c
      DO K=1,IU-IL+1,1
         NC = NC + ONE
         WRITE(CHNC, '(I2)')NC
         IF(CHNC(2:2).EQ.'1')THEN
            CHAUX = 'st' 
         ELSEIF(CHNC(2:2).EQ.'2')THEN
            CHAUX = 'nd'
         ELSEIF(CHNC(2:2).EQ.'3')THEN
            CHAUX = 'rd'
         ELSE
            CHAUX = 'th'
         ENDIF
c
c         IF(IFAIL(K).EQ.ZERO)THEN 
          IF(NC.LT.10 .AND. INFO.EQ.0)THEN
             WRITE(*,1031)NC,CHAUX,EIGVL(K) !,INFO
          ELSEIF(INFO.EQ.0)THEN
             WRITE(*,1032)NC,CHAUX,EIGVL(K) !,INFO
          ELSE
             WRITE(*,1033)NC,CHAUX
          ENDIF
          
c         ENDIF
c     
         DO J=1,N,1
            VI(J) = ZERO
            WK2(J) = ZERO
            DO I=1,NREORT,1
               LNZVC(J,I) = ZERO
            ENDDO
         ENDDO 
c     
         NSEEDAUX = NSEED
c         NSEEDAUX = 1d0
         DO I=1,N,1
            WK1(I) = RAN(NSEEDAUX)
         ENDDO 
c     
         BT0 = ECNORM(N, WK1)
c     
         DO J=1,NREORT,1
c           write(*,*)'j',j,n, it, nreort, mxdct
            DO I=1,N,1
               V0(I) = WK2(I)
               WK2(I) = WK1(I)/BT0
               IF(J.LE.NREORT)THEN
                  LNZVC(I,J) = WK2(I)
c*m1              LNZVC(I,1) = WK2(I)
                  IT = J
               ENDIF
               VI(I) = VI(I) + EIGVC(J,K)*WK2(I)
            ENDDO
c
            CALL AU(DIM, N, NP, SHM, VPOT, WK2, WK1, VAR)
c
            SUM1 = ZERO
            DO I=1,N,1
               WK1(I) = WK1(I) - BT0*V0(I)
               SUM1 = SUM1 + WK2(I)*WK1(I)
            ENDDO 
            ALPH = SUM1
c     
D            write(*,*)'# [DEBUG: LANCZS] i, j, WK1(I)',
            DO I=1,N,1
               WK1(I) = WK1(I) - ALPH*WK2(I)
D               write(*,*)i,j,WK1(I)
            ENDDO  
c
            CALL REORT(N, IT, NREORT, MXDCT, WK1, LNZVC)
c     
            BT0 = ECNORM(N, WK1)
c     
         ENDDO
c     
         DO J=1,N,1
            EIGVC(J,K) = VI(J)
         ENDDO
      ENDDO
c     
      DO J=1,IU-IL+1,1
         DO I=1,N,1
            LNZVC(I,J) = EIGVC(I,J)
         ENDDO
      ENDDO
c      
 1031 FORMAT(4X,I1,A2,1X,'state has been calculating =',E16.8,'/a.u..') 
 1032 FORMAT(3X,I2,A2,1X,'state has been calculating =',E16.8,'/a.u..')
 1033 FORMAT(7X,I2,A2,1X,'state hasn`t been calculating correctly!')  
c     &     ' INFO =',I1)
      RETURN
      END
