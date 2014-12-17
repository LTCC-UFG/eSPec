      SUBROUTINE GETCORR(NAME, N, DTF, TF, U, V)
      IMPLICIT NONE
c     **
c     ** Scalar arguments 
cdel      LOGICAL       
      CHARACTER*(*) NAME
      INTEGER       N
      REAL*8        DTF, TF
c     **
c     ** Array arguments
cdel      LOGICAL       
cdel      CHARACTER*(*) 
cdel      INTEGER       
      REAL*8        U(*), V(*)
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
cdel      REAL*8        ZERO, ONE
cdel      PARAMETER     (ZERO = +0.0D+0, ONE = +1.0D+0)
c     **
c     ** Local scalars 
cdel      LOGICAL       
      CHARACTER*72  READAUX   
      INTEGER       I, NC
cdel      REAL*8        

c     **
c     ** External functions 
cdel      LOGICAL       
cdel      CHARACTER*(*) 
      INTEGER       ICHLENGTH
      REAL*8        T, T0
c     **
c     ** External subroutines 
cdel      EXTERNAL      
c     **
c     ** Intrinsic functions 
cdel      INTRINSIC     
c     .. Start program      
      NC = ICHLENGTH(NAME, 0)
 2    FORMAT(A72)
      OPEN(1, STATUS='OLD', FILE=NAME(1:NC), ERR=996)
c
 10   CONTINUE
      READ(1,2, ERR=998, END= 997)READAUX
      IF(READAUX(1:22).EQ.' Correlation Function:')THEN
         READ(1, 2, ERR=998, END=999)READAUX
         READ(1, 2, ERR=998, END=999)READAUX
c         write(*,*)2**N
         READ(1,*, ERR=998, END=999)T0, U(1), V(1)
         READ(1,*, ERR=998, END=999)T, U(2), V(2)
c         U(1) = U(1)*1.99964686236410649879
c         U(2) = U(2)*1.99964686236410649879
c         V(1) = -V(1)*1.99964686236410649879
c         V(2) = -V(2)*1.99964686236410649879
         DTF =  T - T0
         DO I=3,2**N/2,1
            READ(1,*, ERR=998, END=999)T, U(I), V(I)
c            U(I) = U(I)*1.99964686236410649879
c            V(I) = -V(I)*1.99964686236410649879
c            write(*,*)I, T, U(I), V(I)
         ENDDO
         DO I=2**N/2+1,2**N,1
            READ(1,*, ERR=998, END=999)T, U(I), V(I)
c            U(I) = U(I)*1.99964686236410649879
c            V(I) = V(I)*1.99964686236410649879
c            write(*,*)I, T, U(I), V(I)
         ENDDO
         TF = T + DTF
         RETURN
      ENDIF
      GOTO 10
c
 996  WRITE(*,*)'<<<>>> File ',NAME(1:NC),' cannot be open <<<>>>'
      WRITE(*,*)'Please check and resubimit '
      STOP
 997  WRITE(*,*)'<<<>>> Correlation function wasn`t found <<<>>>'
      STOP
 998  WRITE(*,*)'<<<>>> Error correlation function file <<<>>>'
      STOP
 999  WRITE(*,*)'<<<>>> Error reading correlation function <<<>>>'
      STOP
c
      RETURN
      END
