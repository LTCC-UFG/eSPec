      FUNCTION EF(EFC, E0, T0, TD, TP, OMG, SNI, KL, T)
      IMPLICIT NONE
c     **
c     ** Scalar arguments 
cdel      LOGICAL       
      CHARACTER*(*) EFC
      INTEGER       KL  
      REAL*8        E0, T0, TD, TP, OMG, SNI, T    
c     **
c     ** Array arguments
cdel      LOGICAL       
cdel      CHARACTER*(*) 
cdel      INTEGER       
cdel      REAL*8        
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
*     (24/06/2003) First version EF written by Freddy.
*
c     **
c     ** Parameters 
      REAL*8        ZERO, ONE, TWO, PI, TWOPI
      PARAMETER     (ZERO = +0.0D+0, ONE = +1.0D+0, TWO = +2.0D+0,
     &     PI = +3.14159265358979323846D+0,
     &     TWOPI = +6.2831853071795864769252867663D+0)
c     **
c     ** Local scalars 
cdel      LOGICAL       
cdel      CHARACTER*1   
cdel      INTEGER       
      REAL*8        EF, FAT, FATE, FATT, A1, TPA
c     **
c     ** External functions 
cdel      LOGICAL       
cdel      CHARACTER*(*) 
cdel      INTEGER       
c      REAL*8        GAMMA
c     **
c     ** External subroutines 
cdel      EXTERNAL      
c     **
c     ** Intrinsic functions 
cdel      INTRINSIC     
c     .. Start program
      FAT = +2.1132D-9
      FATE = +2.29371276D+17*3.335664D-28 
      FATT = +1.5192676D+0 
      IF(EFC(1:5).EQ.'.GAUS')THEN
c         IF(T.GT.848.595)E0=0.0D0
c [1] T. Joseph and J. Mans, Mol. Phys., 1986, Vol. 58, No. 6, 1149-1169
         EF = E0*EXP(-4.0D+0*0.693147181D+0*(T - T0)
     &        *(T - T0)/(TP*TP))*COS(OMG*T*FATT + SNI) ![1] 
         EF = EF*FATE ! V/cm = J/(C*cm) -> a.u./D
      ELSEIF(EFC(1:6).EQ.'.GGAUS')THEN
c         write(*,*)'T,E0,OMG,T0,TP,KL',T,E0,OMG,T0,TP,KL
c         read(*,*)
         A1 = ONE/(TWO*KL)
         TPA = TP/(0.693147181**A1)
         EF = FAT*SQRT(E0*EXP(-((T - T0)/TPA)**(TWO*KL)))
     &        *COS(OMG*T*FATT + SNI)
      ELSEIF(EFC(1:5).EQ.'.SIN2')THEN
         EF = E0*COS(PI*(T - TD)/TP)*COS(PI*(T - TD)/TP)
     &        *COS(OMG*(T - TD)*FATT + SNI)
         EF = EF*FATE ! V/cm = J/(C*cm) -> a.u./(D)
      ELSEIF(EFC(1:5).EQ.'.NONE' .OR. EFC(1:5).EQ.'.NULL')THEN
         EF = ZERO
      ELSE
         WRITE(*,1001)
         STOP
      ENDIF
c
 1001 FORMAT('<<<>>> Desired laser pulse wasn´t found. <<<>>>')
c
      RETURN
      END


      FUNCTION GAMMA(X)
      IMPLICIT NONE
c     **
c     ** Scalar arguments 
      REAL*8        X
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
*     (18/03/2004) First version GAMMA written by Freddy.
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
      REAL*8        GAMMA, XPN
c     **
c     ** Local arrays 
      REAL*8        B(8)
c     .. Starting real arrays values
      DATA B / -0.577191652D+0, 0.988205891D+0, -0.897056937D+0, 
     &     0.918206857D+0, -0.756704078D+0, 0.482199394D+0,
     &     -0.193527818D+0, 0.035868343D+0 /
c      
      GAMMA = ONE 
      XPN = ONE
      IF(X.EQ.ZERO)THEN
         RETURN
      ELSE
         DO I=1,8,1
            XPN = XPN*X
            GAMMA = GAMMA + B(I)*XPN
         ENDDO
         GAMMA = GAMMA/X
      ENDIF
c
      RETURN
      END
