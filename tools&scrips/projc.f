      PROGRAM FCfactors
      IMPLICIT NONE
c     **
*     ..
*     Purpose
*     =======
*     Computation of the Frack-Condon factors between two 
*     different electronic states (Y and W)
*
*     FC_ij = <y_i|w_j>; i=0,1,... and j=0,1,...
*
*     ..
*     Description
*     ===========
*     This program computes the Frack-Condon factors between two 
*     different electronic states (Y and W). It is necessary to give 
*     the r dependece (r(1) ... r(N)) of the Mth eigenvectors
*     of each electronic state (Y_M and W_M), which can be calculated with 
*     eSPec. For this, it is necessary to create a input file
*     called FC.dat which has a specific format: 
*
*     =================================================
*     #Comentary line
*     r(1)  y_0(1)  y_1(1)  y_2(1)  .   .   .  y_M(1)  *END*
*     r(2)  y_0(2)  y_1(2)  y_2(2)  .   .   .  y_M(2) 
*     .     .      .      .
*     .     .      .      .
*     r(N)  y_0(N)  y_1(N)  y_2(N)  .   .   .  y_M(N)
*     #Comentary line
*     r(1)  W_0(1)  W_1(1)  W_2(1)  .   .   .  W_M(1)  *END*
*     r(2)  W_0(2)  W_1(2)  W_2(2)  .   .   .  W_M(2) 
*     .     .      .      .
*     .     .      .      .
*     r(N)  W_0(N)  W_1(N)  W_2(N)  .   .   .  W_M(N)
*     
*     ================================================= 
*     OBS: It is necessary to write '*END*' in the end of the
*     first line where is each eigen-vectors are given, as showed
*     above. The last line of the FC.dat file must to be a white line.
*
*     1) compile it giving the comand: f77 -o FC.x FC.f
*     2) run it giving the comand: ./FC.x > FC.out
*
*     ..
*     Authors
*     =======
*     Freddy Fernandes Guimaraes and Viviane Costa Felissicimo.
*
*     ..
*     Historic
*     ========
*     (05/2004) First version written by Freddy 
*     (10/12/2004) Final version written by Viviane
*
c     **
c     ** Parameters 
      INTEGER      ZERO, ONE, MAXR, MAXC
      PARAMETER    (
     &     ZERO = 0.0D+0, 
     &     ONE = 1.0D+0, 
     &     MAXR = 1.0D+6, 
     &     MAXC = 4.0D+1
     &     )
c     **
c     ** Scalars 
      CHARACTER*4  CHNUM
      CHARACTER*30 FILEIN
      CHARACTER*40 INFILE
      CHARACTER*72 LIXO
      INTEGER      MC, NC1, NC2, NL, I, J, K, IOTEST
      REAL*8       SK, CR, CI
c     **
c     ** Arrays
      REAL*8       US(MAXR,MAXC), U1(MAXR), V1(MAXR), PROJ(MAXC)
c     **
c     ** External functions 
cdel      LOGICAL
cdel      CHARACTER*
cdel      INTEGER       
cdel      REAL*8
c     **
c     ** External subroutines 
cdel      EXTERNAL      
c     **
c     ** Intrinsic functions 
cdel      INTRINSIC     
c     .. Start program 
      FILEIN = 'ReImB'    
      MC = 5
      OPEN(1, STATUS='OLD', FILE='vibsts.dat', IOSTAT=IOTEST)
      IF(IOTEST.NE.ZERO)THEN
         WRITE(*,*)'<<<>>> Input file "vibsts.dat" ',
     &        'cannot be opened <<<>>>'
         STOP
      ENDIF
c
      READ(1,102,ERR=999,END=999)LIXO
      READ(1,*,ERR=10,END=999)SK, (US(I,J), J=1,MAXC,1)
 10   CONTINUE
      NC1 = J - ONE
      WRITE(*,*) '# Number of vibrational states:',NC1
      DO I=2,MAXR,1
         READ(1,*,ERR=999,END=20)SK, (US(I,J), J=1,NC1,1)
      ENDDO
 20   CONTINUE
      NL = I - ONE
      WRITE(*,*) '# Number of points of one vibrational state:',NL
c
      DO I=1,10000,1
         WRITE(CHNUM,'(I4)')I
         IF(I.LT.10)THEN
            INFILE = FILEIN(1:MC)//'_000'//CHNUM(4:4)//'.dat'
         ELSEIF(I.LT.100)THEN
            INFILE = FILEIN(1:MC)//'_00'//CHNUM(3:4)//'.dat'
         ELSEIF(I.LT.1000)THEN
            INFILE = FILEIN(1:MC)//'_0'//CHNUM(2:4)//'.dat'
         ELSEIF(I.LT.10000)THEN
            INFILE = FILEIN(1:MC)//'_'//CHNUM(1:4)//'.dat'
         ELSE
            WRITE(*,*)'Too many files to be read. Stopping!'
            STOP
         ENDIF 
         OPEN(2, STATUS='OLD', FILE=INFILE(1:MC+9), IOSTAT=IOTEST)
         IF(IOTEST.NE.ZERO)THEN
            WRITE(*,*)'<<<>>> Input file "',INFILE(1:MC+9),'" ',
     &           'cannot be opened <<<>>>'
            STOP
         ENDIF
         READ(2,102,ERR=998,END=998)LIXO
         DO J=1,NL,1         
            READ(2,*,ERR=998,END=998)SK, U1(J), V1(J) 
         ENDDO
         CLOSE(2)
c
         DO K=1,NC1,1
            CR = ZERO
            CI = ZERO
            DO J=1,NL,1
c               write(*,*)US(J,K),U1(J),V1(J)
               CR = CR + US(J,K)*U1(J)
               CI = CI + US(J,K)*V1(J)
            ENDDO
            PROJ(K) = CR**2 + CI**2
         ENDDO
         WRITE(*,106)LIXO(8:15),(PROJ(K), K=1,NC1,1)
      ENDDO
c
      STOP
 998  WRITE(*,*)'Reading error in file "',INFILE(1:MC+9),'"!'
      STOP
 999  WRITE(*,*)'Reading error in file "FC.dat"!'
      STOP
 102  FORMAT(A72)
 104  FORMAT(1X,'y_i',1X,'->',1X,'w_j',9X, 'FCF')
 105  FORMAT(1X,I3,1X,A2,1X,I3,3X,E15.8)
 106  FORMAT(A8,40(3X,E15.8))
c     
      END
