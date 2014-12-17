c     subroutine PSODAUX
c     ==================
      SUBROUTINE PABMAUX(ABSORB, TPABSOR, DIM, N, CST, THREECST, 
     &     CSGAMMA1, CSGAMMA2, NP, SHM, U1, U2, V1, V2, HU1, HV1,
     &     VPOT, VABC, WORK, VAR)
      IMPLICIT NONE
c     **
c     ** Scalar arguments 
      LOGICAL       ABSORB
      CHARACTER*(*) TPABSOR, DIM
      INTEGER       N
      REAL*8        CST, THREECST, CSGAMMA1, CSGAMMA2
c     **
c     ** Array arguments
cdel      LOGICAL
cdel      CHARACTER*(*)      
      INTEGER       NP(*)
      REAL*8        SHM(*), U1(*), U2(*), V1(*), V2(*), VPOT(*), VABC(*)
      REAL*8        HU1(*), HV1(*), WORK(*), VAR(*)
c     **
*     ..
*     Purpose
*     =======
*     Solve:
*     |p^(2) = |p^(0) + 2*i*dt*|H* |p^(1)/hbar.
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
*     (22/03/2005) First version PABMAUX written by Freddy based on sod.f.
*
c     **
c     ** Parameters 
      REAL*8        ZERO
      PARAMETER     (ZERO = +0.0D+0)
c     **
c     ** Local scalars 
cdel      LOGICAL
cdel      CHARACTER*1
      INTEGER       I
      REAL*8        UAUX, VAUX 
c     **
c     ** External functions 
cdel      LOGICAL
cdel      CHARACTER*(*)
cdel      INTEGER      
cdel      REAL*8
c     **
c     ** External subroutines 
      EXTERNAL      AU
c     **
c     ** Intrinsic functions 
cdel      INTRINSIC     
c     .. Start program  
c      write(*,*)u1(1), u2(1), v1(1), v2(1), vpot(1), vabc(1)
      CALL AU(DIM, N, NP, SHM, VPOT, V2, WORK, VAR)    
c     
      IF(ABSORB)THEN
         IF(TPABSOR(1:7).EQ.'.VOPTIC')THEN
            DO I=1,N,1  
               UAUX = U2(I)
               U2(I) = U1(I) + CST*WORK(I) + CST*VABC(I)*UAUX ! WORK = H*U
               U1(I) = UAUX
            ENDDO
         ELSEIF(TPABSOR(1:8).EQ.'.SMOOTHW')THEN
            DO I=1,N,1  
               UAUX = U2(I)
               U2(I) = VABC(I)*(U1(I) + CST*WORK(I)) ! WORK = H*U
               U1(I) = UAUX
            ENDDO
         ELSE
            WRITE(*,*)'<<<>>> ABM propagation error <<<>>>'
            STOP
         ENDIF
      ELSE
         IF(CSGAMMA1.EQ.ZERO .AND. CSGAMMA2.EQ.ZERO)THEN
            DO I=1,N,1  
               UAUX  = U2(I)
               U2(I) = UAUX + THREECST*WORK(I) - CST*HU1(I) 
               HU1(I) = WORK(I)
               U1(I) = UAUX
            ENDDO
         ELSE
            DO I=1,N,1  
               UAUX  = U2(I)
               U2(I) = UAUX + THREECST*WORK(I) - CSGAMMA1*UAUX
     &              - CST*HU1(I) + CSGAMMA2*U1(I)
               HU1(I) = WORK(I)
               U1(I) = UAUX
            ENDDO
         ENDIF
      ENDIF
c     
      CALL AU(DIM, N, NP, SHM, VPOT, U1, WORK, VAR)
c   
      IF(ABSORB)THEN 
         IF(TPABSOR(1:7).EQ.'.VOPTIC')THEN
            DO I=1,N,1  
               VAUX = V2(I)
               V2(I) = V1(I) - CST*WORK(I) + CST*VABC(I)*VAUX ! WORK = H*U
               V1(I) = VAUX
            ENDDO   
         ELSEIF(TPABSOR(1:8).EQ.'.SMOOTHW')THEN
            DO I=1,N,1  
               VAUX = V2(I)
               V2(I) = VABC(I)*(V1(I) - CST*WORK(I))
               V1(I) = VAUX   
            ENDDO
         ELSE
            WRITE(*,*)'<<<>>> ABM propagation error <<<>>>'
            STOP
         ENDIF
      ELSE
         IF(CSGAMMA1.EQ.ZERO .AND. CSGAMMA2.EQ.ZERO)THEN
            DO I=1,N,1
               VAUX = V2(I)
               V2(I) = VAUX - THREECST*WORK(I) + CST*HV1(I) 
               HV1(I) = WORK(I)
               V1(I) = VAUX
            ENDDO   
         ELSE
            DO I=1,N,1
               VAUX = V2(I)
               V2(I) = VAUX - THREECST*WORK(I) - CSGAMMA1*VAUX 
     &              + CST*HV1(I) + CSGAMMA2*V1(I)
               HV1(I) = WORK(I)
               V1(I) = VAUX
            ENDDO   
         ENDIF
      ENDIF
c
      RETURN
      END
