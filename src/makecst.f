      SUBROUTINE MAKECST(DIM, POTCH, SM, CST) 
      IMPLICIT NONE
c     **
c     ** Scalar arguments 
cdel      LOGICAL       
      CHARACTER*(*) DIM, POTCH
cdel      INTEGER       
cdel      REAL*8        
c     **
c     ** Array arguments
cdel      LOGICAL       
cdel      CHARACTER*(*) 
cdel      INTEGER       
      REAL*8        SM(*), CST(*)
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
*     (01/04/2003) First version MAKECST written by Freddy 
*
c     **
c     ** Parameters 
      REAL*8        ONE, TWO, FOUR, A0A, A0CGS, CCGS, EEV, ECM, TAU
      REAL*8        TWOPI, FFREQ, PI
      PARAMETER     (
     &     ONE = +1.0D+0, 
     &     TWO = +2.0D+0,
     &     FOUR = +4.0D+0,  
     &     FFREQ = +4.556335252750D-6,
     &     A0A = +5.291772083D-1, 
     &     A0CGS = +5.291772083D-9, 
     &     CCGS = +2.99792458D+10,
     &     EEV = +27.211383411D+0,
     &     ECM = +2.19474631371017D+5,
     &     TAU = +2.4188843243D-17,
     &     PI = +3.14159265358979323846D+0,
     &     TWOPI = +6.2831853071795864769252867663D+0
     &     )
c     **
c     ** Local scalars 
cdel      LOGICAL       
cdel      CHARACTER*1   
cdel      INTEGER       
      REAL*8        Al1, AL1A, AL1B, AL2, D1, D1A, D1B, D2, FK1, FK2
      REAL*8        FK3, Re1, Re2, Re3, AL3, D3
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
 1020 FORMAT(1X,'These parameters were already inputed in a.u.')
 1031 FORMAT(1X,'Re_',I1,'/a.u. =',1X,E12.6,';',1X,'K_',I1,'/a.u. =',1X,
     &     E12.6,';')
 1041 FORMAT(1X,'Re_',I1,'/a.u. =',1X,E12.6,';',1X,'alpha_',I1,
     &     '/a.u. =',1X,E12.6,';',1X,'D_',I1,'/a.u. =',1X,E12.6)
 1051 FORMAT(1X,'R(min)_',I1,'/a.u. =',1X,E12.6,';',1X,'R(max)_',I1,
     &     '/a.u. =',1X,E12.6,';',1X,'Const_',I1,'/a.u. =',1X,E12.6)
 1061 FORMAT(1X,'Const_',I1,'/a.u. =',1X,E12.6)    
 1071 FORMAT(1X,'a',1X,'=',1X,E10.4,';',1X,'b',1X,'=',1X,
     &     E10.4,';',1X,'c',1X,'=',1X,E10.4)
      IF(DIM(1:3).EQ.'.1D')THEN 
         IF(POTCH(1:4).EQ.'.OHS')THEN
c     .. ({cm^-1}*cm/s*s/au(t))**2*au(m) = au(m)/au(t)^2 
c     .. The twopi factor cames from the convertion of the units.
            FK1 = (CST(2)*CCGS*TAU*TWOPI)**2*SM(1)
c            write(*,*)'w =',SQRT(FK1/SM)
c     .. {A}/(A/au(L)) = au(L) 
            Re1 = CST(1)/A0A
c     .. rewrite constats in a.u. units
            CST(1) = Re1
            CST(2) = FK1
            WRITE(*,1031)1, CST(1), 1, CST(2)
         ELSEIF(POTCH(1:5).EQ.'.DOHS')THEN
c     .. ({cm^-1}*cm/s*s/au(t))**2*au(m) = au(m)/au(t)^2 
            FK1 = (CST(2)*FFREQ)**2*SM(1)
            FK2 = (CST(5)*CCGS*TAU*TWOPI)**2*SM(1)/(A0A*1.0D-10)**2
c     .. {A}/(A/au(L)) = au(L) 
            Re1 = CST(1)/A0A 
            Re2 = CST(3)/A0A
            Re3 = CST(4)/A0A
c     .. rewrite constats in a.u. units
            CST(1) = Re1
            CST(2) = FK1
            CST(3) = Re2
            CST(4) = Re3
            CST(5) = FK2
            WRITE(*,1031)1, CST(1), 1, CST(2)
            WRITE(*,1031)2, CST(3), 2, CST(4)
         ELSEIF(POTCH(1:16).EQ.'.MORSE_Dalpha_au' .OR.
     &           POTCH(1:21).EQ.'.ANTI-MORSE_Dalpha_au' .OR.
     &           POTCH(1:17).EQ.'.MORSEM_Dalpha_au' .OR. 
     &           POTCH(1:22).EQ.'.ANTI-MORSEM_Dalpha_au' .OR.
     &           POTCH(1:20).EQ.'.MOD-MORSE_Dalpha_au')THEN
            WRITE(*,1020)
            WRITE(*,1041)1, CST(1), 1, CST(2), 1, CST(3)
         ELSEIF(POTCH(1:13).EQ.'.MORSE_Dalpha' .OR.
     &           POTCH(1:18).EQ.'.ANTI-MORSE_Dalpha' .OR.
     &           POTCH(1:14).EQ.'.MORSEM_Dalpha' .OR. 
     &           POTCH(1:19).EQ.'.ANTI-MORSEM_Dalpha' .OR.
     &           POTCH(1:17).EQ.'.MOD-MORSE_Dalpha')THEN
            Re1 = CST(1)/A0A
            Al1 = CST(2)*A0A
            D1 = CST(3)/EEV
c     .. rewrite constats in a.u. units
            CST(1) = Re1
            CST(2) = Al1
            CST(3) = D1
            WRITE(*,1041)1, CST(1), 1, CST(2), 1, CST(3)
         ELSEIF(POTCH(1:6).EQ.'.MORSE' .OR.
     &           POTCH(1:11).EQ.'.ANTI-MORSE' .OR.
     &           POTCH(1:7).EQ.'.MORSEM' .OR. 
     &           POTCH(1:12).EQ.'.ANTI-MORSEM' .OR.
     &           POTCH(1:10).EQ.'.MOD-MORSE' .OR.
     &           POTCH(1:11).EQ.'.MORSE_CM^1' .OR.
     &           POTCH(1:16).EQ.'.ANTI-MORSE_CM^1' .OR.
     &           POTCH(1:12).EQ.'.MORSEM_CM^1' .OR. 
     &           POTCH(1:17).EQ.'.ANTI-MORSEM_CM^1' .OR.
     &           POTCH(1:15).EQ.'.MOD-MORSE_CM^1')THEN
c     .. {A}/(A/au(L)) = au(L) 
            Re1 = CST(1)/A0A
c     .. ({cm^-1}*cm/s*s/au(t))**2*au(m) = au(m)/au(t)^2 
            FK1 = (CST(2)*CCGS*TAU*TWOPI)**2*SM(1)
c     .. 
c            D1 = CST(3)/EEV
            D1 = TWOPI*CCGS*TAU*CST(2)**2/(FOUR*CST(3))
c     .. 
            Al1 = SQRT(FK1/(TWO*D1))
c            write(*,*)al1,d1,TWOPI*CCGS*TAU*(1580.19)**2/(4*11.98)
c            read(*,*)
c     .. rewrite constats in a.u. units
            CST(1) = Re1
            CST(2) = Al1
            CST(3) = D1
            WRITE(*,1041)1, CST(1), 1, CST(2), 1, CST(3)
         ELSEIF(POTCH(1:11).EQ.'.BARRIER_au' .OR. 
     &           POTCH(1:8).EQ.'.WELL_au')THEN
            WRITE(*,1020)
            WRITE(*,1051)1, CST(1), 1, CST(2), 1, CST(3)
         ELSEIF(POTCH(1:9).EQ.'.CONST_au')THEN
            WRITE(*,1020)
            WRITE(*,1061)1, CST(1)
         ELSEIF(POTCH(1:8).EQ.'.BARRIER' .OR. 
     &           POTCH(1:5).EQ.'.WELL')THEN
            CST(1) = CST(1)/A0A
            CST(2) = CST(2)/A0A
            CST(3) = CST(3)/EEV
            WRITE(*,1051)1, CST(1), 1, CST(2), 1, CST(3)
         ELSEIF(POTCH(1:6).EQ.'.CONST' .OR. POTCH(1:5).EQ.'.FREE')THEN
            CST(1) = CST(1)/EEV
            WRITE(*,1061)1, CST(1)
         ELSE
            WRITE(*,1011)
            STOP
         ENDIF  
      ELSEIF(DIM(1:3).EQ.'.2D')THEN
         IF(POTCH(1:4).EQ.'.OHS')THEN
c     .. ({cm^-1}*cm/s*s/au(t))**2*au(m) = au(m)/au(t)^2 
            FK1 = (CST(2)*CCGS*TAU*TWOPI)**2*SM(1)
            FK2 = (CST(4)*CCGS*TAU*TWOPI)**2*SM(2)
c     .. {A}/(A/au(L)) = au(L) 
            Re1 = CST(1)/A0A
            Re2 = CST(3)/A0A
c     .. rewrite constats in a.u. units
            CST(1) = Re1
            CST(3) = Re2
            CST(2) = FK1
            CST(4) = FK2
            WRITE(*,1031)1, CST(1), 1, CST(2)
            WRITE(*,1031)2, CST(3), 2, CST(4)
         ELSEIF(POTCH(1:7).EQ.'.LSM_au')THEN
            WRITE(*,1020)
            WRITE(*,1041)1, CST(1), 1, CST(2), 1, CST(3)
            WRITE(*,1041)2, CST(1), 2, CST(4), 2, CST(5)
            WRITE(*,1041)3, CST(6), 3, CST(7), 3, CST(8)
         ELSEIF(POTCH(1:4).EQ.'.LSM')THEN
c     .. {A}/(A/au(L)) = au(L) 
            Re1 = CST(1)/A0A
            Re2 = CST(6)/A0A
c     .. {A^-1}*(A/au(L)) = au(L^{-1})     
            AL1A = CST(2)*A0A 
            AL1B = CST(4)*A0A 
            AL2 = CST(7)*A0A
c     .. {cm^{-1}}*(au(E)*cm) = au(E)              
            D1A = CST(3)/ECM
            D1B = CST(5)/ECM
            D2 = CST(8)/ECM
c     .. rewrite constats in a.u. units
            CST(1) = Re1
            CST(2) = AL1A
            CST(3) = D1A
            CST(4) = AL1B
            CST(5) = D1B
            CST(6) = Re2 
            CST(7) = AL2
            CST(8) = D2 
            WRITE(*,1041)1, CST(1), 1, CST(2), 1, CST(3)
            WRITE(*,1041)2, CST(1), 2, CST(4), 2, CST(5)
            WRITE(*,1041)3, CST(6), 3, CST(7), 3, CST(8)
         ELSEIF(POTCH(1:6).EQ.'.LS_au')THEN
            WRITE(*,1020)
            WRITE(*,1041)1, CST(1), 1, CST(2), 1, CST(3)
            WRITE(*,1041)2, CST(1), 2, CST(4), 2, CST(5)
         ELSEIF(POTCH(1:3).EQ.'.LS')THEN
c     .. {A}/(A/au(L)) = au(L) 
            Re1 = CST(1)/A0A
c     .. {A^-1}*(A/au(L)) = au(L^{-1})     
            AL1A = CST(2)*A0A 
            AL1B = CST(4)*A0A 
c     .. {cm^{-1}}*(au(E)*cm) = au(E)              
            D1A = CST(3)/ECM
            D1B = CST(5)/ECM
c     .. rewrite constats in a.u. units
            CST(1) = Re1
            CST(2) = AL1A
            CST(3) = D1A
            CST(4) = AL1B
            CST(5) = D1B
            WRITE(*,1041)1, CST(1), 1, CST(2), 1, CST(3)
            WRITE(*,1041)2, CST(1), 2, CST(4), 2, CST(5)
         ELSEIF(POTCH(1:7).EQ.'.LBH_CC')THEN
            WRITE(*,*) 'Not implemented!'
            STOP
         ELSEIF(POTCH(1:4).EQ.'.LBH')THEN
            WRITE(*,*) 'Not implemented!'
            STOP
         ELSEIF(POTCH(1:18).EQ.'.LEPS_CC_Dalpha_au' .OR.
     &           POTCH(1:15).EQ.'.LEPS_Dalpha_au')THEN
            WRITE(*,1020)
            WRITE(*,1041)1, CST(1), 1, CST(2), 1, CST(3)
            WRITE(*,1041)2, CST(4), 2, CST(5), 2, CST(6)
            WRITE(*,1041)3, CST(7), 3, CST(8), 3, CST(9)
            WRITE(*,1071)CST(10), CST(11), CST(12)
            CST(13) = PI*CST(13)/1.8D+2
            WRITE(*,*)'Angle (rad) =',CST(13)
         ELSEIF(POTCH(1:15).EQ.'.LEPS_CC_Dalpha' .OR.
     &           POTCH(1:12).EQ.'.LEPS_Dalpha')THEN
            Re1 = CST(1)/A0A
            Al1 = CST(2)*A0A
            D1 = CST(3)/EEV
            Re2 = CST(4)/A0A
            Al2 = CST(5)*A0A
            D2 = CST(6)/EEV
            Re3 = CST(7)/A0A
            Al3 = CST(8)*A0A
            D3 = CST(9)/EEV
            CST(1) = Re1
            CST(2) = Al1
            CST(3) = D1
            CST(4) = Re2
            CST(5) = Al2
            CST(6) = D2
            CST(7) = Re3
            CST(8) = Al3
            CST(9) = D3
            WRITE(*,1041)1, CST(1), 1, CST(2), 1, CST(3)
            WRITE(*,1041)2, CST(4), 2, CST(5), 2, CST(6)
            WRITE(*,1041)3, CST(7), 3, CST(8), 3, CST(9)
            WRITE(*,1071)CST(10), CST(11), CST(12)
         ELSEIF(POTCH(1:13).EQ.'.BARRIERXY_au' .OR. 
     &           POTCH(1:10).EQ.'.WELLXY_au')THEN
            WRITE(*,1020)
            WRITE(*,1051)1, CST(1), 1, CST(2), 1, CST(3)
            WRITE(*,1051)2, CST(4), 2, CST(5), 2, CST(6)
         ELSEIF(POTCH(1:12).EQ.'.BARRIERX_au' .OR. 
     &           POTCH(1:12).EQ.'.BARRIERY_au' .OR.
     &           POTCH(1:9).EQ.'.WELLX_au' .OR.
     &           POTCH(1:9).EQ.'.WELLY_au')THEN
            WRITE(*,1020)
            IF(POTCH(1:9).EQ.'.BARRIERX' .OR. 
     &           POTCH(1:6).EQ.'.WELLX')THEN
               WRITE(*,1051)1, CST(1), 1, CST(2), 1, CST(3)
            ELSE
               WRITE(*,1051)2, CST(1), 2, CST(2), 2, CST(3)
            ENDIF
         ELSEIF(POTCH(1:11).EQ.'.BARRIER_au' .OR. 
     &           POTCH(1:8).EQ.'.WELL_au')THEN
            WRITE(*,1020)
            WRITE(*,1051)1, CST(1), 1, CST(2), 1, CST(5)
            WRITE(*,1051)2, CST(3), 2, CST(4), 2, CST(5)
         ELSEIF(POTCH(1:9).EQ.'.CONST_au')THEN
            WRITE(*,1020)
            WRITE(*,1061)1, CST(1)
         ELSEIF(POTCH(1:10).EQ.'.BARRIERXY' .OR. 
     &           POTCH(1:7).EQ.'.WELLXY')THEN
            CST(1) = CST(1)/A0A
            CST(2) = CST(2)/A0A
            CST(3) = CST(3)/EEV 
            CST(4) = CST(4)/A0A
            CST(5) = CST(5)/A0A
            CST(6) = CST(6)/EEV 
            WRITE(*,1051)1, CST(1), 1, CST(2), 1, CST(3)
            WRITE(*,1051)2, CST(4), 2, CST(5), 2, CST(6)
         ELSEIF(POTCH(1:9).EQ.'.BARRIERX' .OR. 
     &           POTCH(1:9).EQ.'.BARRIERY' .OR.
     &           POTCH(1:6).EQ.'.WELLX' .OR.
     &           POTCH(1:6).EQ.'.WELLY')THEN
            CST(1) = CST(1)/A0A
            CST(2) = CST(2)/A0A
            CST(3) = CST(3)/EEV
            IF(POTCH(1:9).EQ.'.BARRIERX' .OR. 
     &           POTCH(1:6).EQ.'.WELLX')THEN
               WRITE(*,1051)1, CST(1), 1, CST(2), 1, CST(3)
            ELSE
               WRITE(*,1051)2, CST(1), 2, CST(2), 2, CST(3)
            ENDIF
         ELSEIF(POTCH(1:8).EQ.'.BARRIER' .OR. 
     &           POTCH(1:5).EQ.'.WELL')THEN
            CST(1) = CST(1)/A0A
            CST(2) = CST(2)/A0A
            CST(3) = CST(3)/A0A
            CST(4) = CST(4)/A0A
            CST(5) = CST(5)/EEV
            WRITE(*,1051)1, CST(1), 1, CST(2), 1, CST(5)
            WRITE(*,1051)2, CST(3), 2, CST(4), 2, CST(5)
         ELSEIF(POTCH(1:6).EQ.'.CONST' .OR. POTCH(1:6).EQ.'.FREE')THEN
            CST(1) = CST(1)/EEV
            WRITE(*,1061)1, CST(1)
         ELSE
            WRITE(*,1011)
         ENDIF  
      ELSEIF(DIM(1:3).EQ.'.3D')THEN 
         IF(POTCH(1:4).EQ.'.OHS')THEN
c     .. ({cm^-1}*cm/s*s/au(t))**2*au(m) = au(m)/au(t)^2 
         FK1 = (CST(2)*CCGS*TAU*TWOPI)**2*SM(1)
         FK2 = (CST(4)*CCGS*TAU*TWOPI)**2*SM(2)
         FK3 = (CST(6)*CCGS*TAU*TWOPI)**2*SM(3)
c     .. {A}/(A/au(L)) = au(L) 
         Re1 = CST(1)/A0A
         Re2 = CST(3)/A0A
         Re3 = CST(5)/A0A
c     .. rewrite constats in a.u. units
         CST(1) = Re1
         CST(3) = Re2
         CST(5) = Re3
         CST(2) = FK1
         CST(4) = FK2
         CST(6) = FK3  
         WRITE(*,1031)1, CST(1), 1, CST(2)
         WRITE(*,1031)2, CST(3), 2, CST(4)
         WRITE(*,1031)3, CST(5), 3, CST(6)
         ELSEIF(POTCH(1:11).EQ.'.BARRIER_au' .OR. 
     &           POTCH(1:8).EQ.'.WELL_au')THEN
            WRITE(*,1020)
            WRITE(*,1051)1, CST(1), 1, CST(2), 1, CST(7)
            WRITE(*,1051)2, CST(3), 2, CST(4), 2, CST(7)
            WRITE(*,1051)2, CST(5), 2, CST(6), 2, CST(7)
         ELSEIF(POTCH(1:9).EQ.'.CONST_au')THEN
            WRITE(*,1020)
            WRITE(*,1061)1, CST(1)
         ELSEIF(POTCH(1:8).EQ.'.BARRIER' .OR. 
     &           POTCH(1:5).EQ.'.WELL')THEN
            CST(1) = CST(1)/A0A
            CST(2) = CST(2)/A0A
            CST(3) = CST(3)/A0A
            CST(4) = CST(4)/A0A
            CST(5) = CST(5)/A0A
            CST(6) = CST(6)/A0A
            CST(7) = CST(7)/EEV
            WRITE(*,1051)1, CST(1), 1, CST(2), 1, CST(7)
            WRITE(*,1051)2, CST(3), 2, CST(4), 2, CST(7)
            WRITE(*,1051)3, CST(5), 3, CST(6), 3, CST(7)
         ELSEIF(POTCH(1:18).EQ.'.LEPS_CC_Dalpha_au')THEN
            WRITE(*,1020)
            WRITE(*,1041)1, CST(1), 1, CST(2), 1, CST(3)
            WRITE(*,1041)2, CST(4), 2, CST(5), 2, CST(6)
            WRITE(*,1041)3, CST(7), 3, CST(8), 3, CST(9)
            WRITE(*,*)'a', CST(10), 'b', CST(11), 'c', CST(12)
         ELSEIF(POTCH(1:15).EQ.'.LEPS_CC_Dalpha')THEN              
            Re1 = CST(1)/A0A
            Al1 = CST(2)*A0A
            D1 = CST(3)/EEV
            Re2 = CST(4)/A0A
            Al2 = CST(5)*A0A
            D2 = CST(6)/EEV
            Re3 = CST(7)/A0A
            Al3 = CST(8)*A0A
            D3 = CST(9)/EEV
            CST(1) = Re1
            CST(2) = Al1
            CST(3) = D1
            CST(4) = Re2
            CST(5) = Al2
            CST(6) = D2
            CST(7) = Re3
            CST(8) = Al3
            CST(9) = D3
            WRITE(*,1041)1, CST(1), 1, CST(2), 1, CST(3)
            WRITE(*,1041)2, CST(4), 2, CST(5), 2, CST(6)
            WRITE(*,1041)3, CST(7), 3, CST(8), 3, CST(9)
            WRITE(*,*)'a', CST(10), 'b', CST(11), 'c', CST(12)
         ELSEIF(POTCH(1:6).EQ.'.CONST' .OR. POTCH(1:5).EQ.'.FREE')THEN
            CST(1) = CST(1)/EEV
            WRITE(*,1061)1, CST(1)
         ELSE    
            WRITE(*,1011)
         ENDIF   
      ELSE
         WRITE(*,1011)
         WRITE(*,*)'<<<>>> Dimesion error!  <<<>>>'
         STOP
      ENDIF
 1011 FORMAT('<<<>>> Potential parameter converting error <<<>>>') 
c     ..
      RETURN
      END
