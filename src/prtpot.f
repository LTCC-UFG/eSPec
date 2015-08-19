      SUBROUTINE PRTPT(NAME, MC, ND, T, NP, XP, XI, SH, V)
      IMPLICIT NONE
c     **
c     ** Scalar arguments 
cdel      LOGICAL       
      CHARACTER*(*) NAME
      INTEGER       MC, ND
      REAL*8        T
c     **
c     ** Array arguments
cdel      LOGICAL       
cdel      CHARACTER*(*) 
      INTEGER       NP(*)
      REAL*8        XP(*) , XI(*) , SH(*), V(*) 
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
      REAL*8        ZERO, ONE
      PARAMETER     (ZERO = +0.0D+0, ONE = +1.0D+0)
c     **
c     ** Local scalars 
cdel      LOGICAL       
cdel      CHARACTER*1   
      INTEGER       I, J, K, L, NC
      REAL*8        WRD

c     **
c     ** External functions 
cdel      LOGICAL       
cdel      CHARACTER*(*) 
      INTEGER       ICHLENGTH
cdel      REAL*8        
c     **
c     ** External subroutines 
cdel      EXTERNAL      
c     **
c     ** Intrinsic functions 
cdel      INTRINSIC     
c     .. Start program 
D      write(*,*)sh(1),sh(2),np(1),np(2)
D      read(*,*)
      IF(MC.EQ.9)THEN
         NC = ICHLENGTH(NAME, 0)
         OPEN(UNIT=9,STATUS='UNKNOWN',FILE=NAME(1:NC))
         WRITE(MC,1011)('#R',I,'#', I=ND,1,-1),'#V',0,'#'
      ELSE
         WRITE(MC,1011)('*R',I,'*', I=ND,1,-1),'*V',0,'*'
      ENDIF
            
      DO L=1,NP(3),1
         XP(3) = XI(3) + (L - 1)*SH(3)
         DO K=1,NP(2),1
            XP(2) = XI(2) + (K - 1)*SH(2)
            DO J=1,NP(1),1
               XP(1) = XI(1) + (J - 1)*SH(1)
               NC = J + (K - ONE)*NP(1) + (L - ONE)*NP(1)*NP(2)
               IF(ABS(V(NC)).LT.1D-99)THEN
                  WRD = ZERO
               ELSE
                  WRD = V(NC)
               ENDIF
               WRITE(MC,1012)(XP(I), I=ND,1,-1),WRD
            ENDDO
            IF(MC.EQ.9)WRITE(MC,*)
         ENDDO
         IF(MC.EQ.9)WRITE(MC,*)
      ENDDO
c
      IF(MC.EQ.9)THEN
         WRITE(MC,*)'# Time/fs = ',T
         CLOSE(9)
      ENDIF
c
 1011 FORMAT(5X,100(A2,I1,A1,11X)) 
 1012 FORMAT(1X,100(E12.6,3X))
c
      RETURN
      END



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC prints only a few regions


      SUBROUTINE PRTPTREG(NAME, MC, ND, T, NP, XP, XI, SH, V, NREG,
     &     RANGE)
      IMPLICIT NONE
c     **
c     ** Scalar arguments 
cdel      LOGICAL       
      CHARACTER*(*) NAME
      INTEGER       II,MC, ND,CK(3)
      REAL*8        T
c     **
c     ** Array arguments
cdel      LOGICAL       
cdel      CHARACTER*(*) 
      INTEGER       NP(*)
      REAL*8        XP(*) , XI(*) , SH(*), V(*) 
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
      REAL*8        ZERO, ONE
      PARAMETER     (ZERO = +0.0D+0, ONE = +1.0D+0)
c     **
c     ** Local scalars 
cdel      LOGICAL       
cdel      CHARACTER*1   
      INTEGER       I, J, K, L, NC,NREG
      REAL*8        WRD,RANGE(5,7)

c     **
c     ** External functions 
cdel      LOGICAL       
cdel      CHARACTER*(*) 
      INTEGER       ICHLENGTH
cdel      REAL*8        
c     **
c     ** External subroutines 
cdel      EXTERNAL      
c     **
c     ** Intrinsic functions 
cdel      INTRINSIC     
c     .. Start program 
D      write(*,*)sh(1),sh(2),np(1),np(2)
D      read(*,*)
      IF(MC.EQ.9)THEN
         NC = ICHLENGTH(NAME, 0)
         OPEN(UNIT=9,STATUS='UNKNOWN',FILE=NAME(1:NC))
         WRITE(MC,1011)('#R',I,'#', I=ND,1,-1),'#V',0,'#'
      ELSE
         WRITE(MC,1011)('*R',I,'*', I=ND,1,-1),'*V',0,'*'
      ENDIF
      DO II=1,NREG,1
         WRITE(MC,*)'# region ',II
c
         DO L=1,NP(3),1
            XP(3) = XI(3) + (L - 1)*SH(3)
            CK(2) = 0.0D+0;
            IF(XP(2).GE.RANGE(II,3).AND.XP(2).LE.RANGE(II,4).AND.ND.EQ.3
     &           )THEN
               CK(2) = 1.0D+0
            ENDIF
            DO K=1,NP(2),1
               XP(2) = XI(2) + (K - 1)*SH(2)
               CK(1) = 0.0e+0;
               IF(XP(2).GE.RANGE(II,3).AND.XP(2).LE.RANGE(II,4))THEN
                  CK(1) = 1.0D+0
               ENDIF
               DO J=1,NP(1),1
                  XP(1) = XI(1) + (J - 1)*SH(1)
                  NC = J + (K - ONE)*NP(1) + (L - ONE)*NP(1)*NP(2)
                  IF(ABS(V(NC)).LT.1D-99)THEN
                     WRD = ZERO
                  ELSE
                     WRD = V(NC)
                  ENDIF
                  IF(XP(1).GE.RANGE(II,1).AND.XP(1).LE.RANGE(II,2) 
     &                 .AND. XP(2).GE.RANGE(II,3)
     &                 .AND. XP(2).LE.RANGE(II,4)
     &                 .AND. XP(3).GE.RANGE(II,5)
     &                 .AND.XP(3).LE.RANGE(II,6))THEN
                     WRITE(MC,1012)(XP(I), I=ND,1,-1),WRD
                  ENDIF
               ENDDO
               IF(MC.EQ.9.AND.CK(1).EQ.1)WRITE(MC,*)
            ENDDO
            IF(MC.EQ.9.AND.CK(2).EQ.1)WRITE(MC,*)
         ENDDO
c
         WRITE(MC,*)
      ENDDO
c
      IF(MC.EQ.9)THEN
         WRITE(MC,*)'# Time/fs = ',T
         CLOSE(9)
      ENDIF
c
 1011 FORMAT(5X,100(A2,I1,A1,11X)) 
 1012 FORMAT(1X,100(E12.6,3X))
c
      RETURN
      END
