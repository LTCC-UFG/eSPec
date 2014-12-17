      SUBROUTINE DIPLMMT(TPDIPL, DMFILE, N, ND, CSTI, CSTF, CSTDM, XI, 
     &     XF, SH, DM)
      IMPLICIT NONE
c     **
c     ** Scalar arguments 
cdel      LOGICAL       
      CHARACTER*(*) TPDIPL, DMFILE
      INTEGER       N, ND    
cdel      REAL*8      
c     **
c     ** Array arguments
cdel      LOGICAL       
cdel      CHARACTER*(*) 
cdel      INTEGER       
      REAL*8        CSTI(*), CSTF(*), CSTDM(*), XI(*), XF(*), SH(*)
      REAL*8        DM(*)
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
      INTEGER       I
      REAL*8        RAB, FATDM, XMIN

c     **
c     ** External functions 
cdel      LOGICAL       
cdel      CHARACTER*(*) 
cdel      INTEGER       
cdel      REAL*8        
c     **
c     ** External subroutines 
cdel      EXTERNAL      
c     **
c     ** Intrinsic functions 
cdel      INTRINSIC     
c     .. Start program    
c      FATDM = +1.0D-18*2.99792458D+9*1.602176462D-19 
      FATDM = 3.335664D-28
      IF(TPDIPL(1:6).EQ.'.CONST')THEN
         DO I=1,N,1
            RAB = XI(1) + (I - 1)*SH(1)
            DM(I) = CSTDM(1)/RAB
         ENDDO
      ELSEIF(TPDIPL(1:5).EQ.'.EXP+')THEN
         DO I=1,N,1
            RAB = XI(1) + (I - 1)*SH(1)
            DM(I) = CSTDM(1)*(EXP(-CSTDM(2)*RAB) - EXP(-3*CSTDM(2)*RAB))
     &           /(EXP(-CSTDM(2)*CSTF(1)) - EXP(-3*CSTDM(2)*CSTF(1)))
            write(21,*)RAB, DM(I)
         ENDDO
      ELSEIF(TPDIPL(1:4).EQ.'.EXP')THEN
         DO I=1,N,1
            RAB = XI(1) + (I - 1)*SH(1)
            DM(I) = CSTDM(1)*RAB*EXP(-CSTDM(2)*RAB)
         ENDDO
      ELSEIF(TPDIPL(1:5).EQ.'.READ')THEN
         CALL RDPTE(DMFILE, N, ND, XI, XF, XMIN, DM)
      ELSEIF(TPDIPL(1:5).EQ.'.NONE' .OR. TPDIPL(1:5).EQ.'.NULL')THEN 
         CONTINUE
      ELSE
         WRITE(*,1001)
         STOP
      ENDIF
c
c      DO I=1,N,1
c         DM(I) = DM(I)          ! debye
c*FATDM ! debye -> C cm 
c      ENDDO
c
 1001 FORMAT('<<<>>> Desired diplole momentum wasn`t found <<<>>>')
c
      RETURN
      END
