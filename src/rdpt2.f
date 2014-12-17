c     .. Read potential from file called INFILE
c     =========================================
      SUBROUTINE RDPT2(INFILE, NT, ND, XI, XF, U, V)
      IMPLICIT NONE
c     **
c     ** Scalar arguments 
cdel      LOGICAL
      CHARACTER*(*) INFILE
      INTEGER       INFLENGTH, NT, ND
cdel       REAL*8      
c     **
c     ** Array arguments
cdel       LOGICAL
cdel       CHARACTER      
cdel       INTEGER  
      REAL*8        XI(*), XF(*), U(*), V(*) 
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
*     Historic
*     ========
*     (14/03/2003) First version RDPT written by Freddy. 
*
c     **
c     ** Parameters 
cdel      REAL*8        ZERO, ONE
cdel      PARAMETER     (ZERO = +0.0D+0, ONE = +1.0D+0)
c     **
c     ** Local scalars 
cdel       LOGICAL
c      CHARACTER*1   
      CHARACTER*72  TITLE
      INTEGER       I, J, IOTEST
c      REAL*8        VAUX

c     **
c     ** External functions 
cdel       LOGICAL
cdel       CHARACTER
      INTEGER       ICHLENGTH    
cdel       REAL*8
c     **
c     ** External subroutines 
cdel       EXTERNAL      
c     **
c     ** Intrinsic functions 
cdel       INTRINSIC     
c     .. Start program
c     ..
      INFLENGTH = ICHLENGTH(INFILE, 0)
c
      OPEN(69, STATUS='OLD', FILE=INFILE(1:INFLENGTH), IOSTAT=IOTEST)
      IF(IOTEST.NE.0)THEN
         WRITE(*,*)'<<<>>> Input file "', INFILE(1:INFLENGTH), 
     &        '" cannot be opened <<<>>>'
         STOP
      ENDIF
      READ(69,'(A)',END=999,ERR=999)TITLE
c
      READ(69,*,END=999,ERR=999)(XI(I), I=ND,1,-1), U(1), V(1)
      DO J=2,NT-1,1
         READ(69,*,END=999,ERR=999)(XF(I), I=ND,1,-1), U(J),V(J)
      ENDDO
      READ(69,*,END=999,ERR=999)(XF(I), I=ND,1,-1), U(NT), V(NT)
c
c 992  FORMAT(1X,'End of file, reading finish!')
c      WRITE(*,992)
c     .. 
      RETURN
 999  WRITE(*,*)'<<<>>> Reading file error, check and resubmit <<<>>>'
      WRITE(*,*)INFILE
      STOP
      END


	
	
