C      PROGRAM MAIN
C      IMPLICIT NONE
C      DOUBLE PRECISION, DIMENSION(3) :: R
C      DOUBLE PRECISION :: V,ANGS2AU
C      
C      ANGS2AU=0.529177249D+00
C      
C      R(1)=1.270543D+00/ANGS2AU 
C      R(2)=1.843761D+00/ANGS2AU
C      R(3)=1.843761D+00/ANGS2AU
C      
C      CALL QFF(R,V) 
C      
C      !PRINT*, V
C      
C      END PROGRAM MAIN
C
C##############################################################################
C QFF MODEL FUNCTION FOR AB2 MOLECULES
C##############################################################################
C R1=B-B IN ANGS 
C R2=A-B IN ANGS
C R3=A-B IN ANGS
C V0 IS IN Eh
C RREF IS IN ANGS
C THETAREF IS IN DEGS
C POTENTIAL IN Eh
C##############################################################################

      SUBROUTINE QFF(RAU,POTEH)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (N=3)
      DIMENSION RAU(N),R(N),S(N)
      DIMENSION FIJ(N,N)
      DIMENSION FIJK(N,N,N)
      DIMENSION FIJKL(N,N,N,N)
      DIMENSION C(22)
      
       C(  1)=  0.860183731159355958D+01
       C(  2)=  0.370148748695581409D+01
       C(  3)=  0.609465226419907768D+01
       C(  4)=  0.541584890384679341D+01
       C(  5)= -0.372384504948115662D+02
       C(  6)= -0.116436894944383234D+02
       C(  7)= -0.233909471039851589D+02
       C(  8)= -0.267039947528282404D+02
       C(  9)= -0.465351430906939498D+02
       C( 10)= -0.202609953719755698D+01
       C( 11)=  0.133994860135151299D+03
       C( 12)=  0.292795021731294689D+02
       C( 13)=  0.544000520240434113D+02
       C( 14)=  0.103141780873540341D+03
       C( 15)=  0.119504077440915211D+03
       C( 16)= -0.189583683303783062D+02
       C( 17)=  0.316235674797988679D+03
       C( 18)= -0.415885347649561581D+02
       C( 19)=  0.740690985067943757D+02
       C( 20)=  0.135666380120678376D+01
       C( 21)=  0.647456862504041055D+02
       C( 22)= -0.130889856789944417D+03
     
      RREF=C(20)
      THETAREF=C(21)
      V0=C(22)
       
      POTEH=0.000000D+00
      PI=ACOS(-1.0D+00)
      ANGS2AU=0.529177249D+00
      FAC=0.22937104486906D+00 ! CONVERSION FROM aJ TO Eh
      
      R(1)=RAU(1)!*ANGS2AU
      R(2)=RAU(2)!*ANGS2AU
      R(3)=RAU(3)!*ANGS2AU
      
      THETAREF=THETAREF*PI/180.0D+00
      
      COSTHETA=(R(1)**2-R(2)**2-R(3)**2)/(-2.0D+00*R(2)*R(3))
      IF (COSTHETA>1.00D+00) THEN 
        COSTHETA=1.00D+00
      ELSE IF (COSTHETA<-1.00D+00) THEN 
        COSTHETA=-1.00D+00
      END IF 

      THETA=ACOS(COSTHETA)

      S(1)=(1.0D+00/SQRT(2.0D+00))*((R(2)-RREF)+(R(3)-RREF))
      S(2)=THETA-THETAREF
      S(3)=(1.0D+00/SQRT(2.0D+00))*((R(2)-RREF)-(R(3)-RREF)) 
C
C QUADRATIC PART
C      
      DO I=1,N     
        DO J=1,N
          FIJ(I,J)=0.0D+00
          !WRITE(*,FMT='(A4,I1,A,I1,A2)') "FIJ(",I,",",J,")="
        END DO
      END DO
      
      FIJ(1,1)=C(1)
      FIJ(1,2)=C(2)
      FIJ(1,3)=0.0D+00
      FIJ(2,1)=FIJ(1,2)
      FIJ(2,2)=C(3)
      FIJ(2,3)=0.0D+00
      FIJ(3,1)=FIJ(1,3)
      FIJ(3,2)=FIJ(2,3)
      FIJ(3,3)=C(4)
      
      V2=0.000000D+00
      DO I=1,N
        DO J=1,N
          V2=V2+(1.0D+00/2.0D+00)*FIJ(I,J)*S(I)*S(J)
          !PRINT*, FIJ(I,J),S(I),S(J),I,J,V2
        END DO    
      END DO
C
C CUBIC PART
C
      DO I=1,N      
        DO J=1,N
          DO K=1,N
          FIJK(I,J,K)=0.0D+00
          !WRITE(*,FMT='(A5,I1,A,I1,A,I1,A2)') "FIJK(",I,",",J,",",K,")="
          END DO
        END DO
      END DO
      
      FIJK(1,1,1)=C(5)
      FIJK(1,1,2)=C(6)
      FIJK(1,1,3)=0.0D+00
      FIJK(1,2,1)=FIJK(1,1,2)
      FIJK(1,2,2)=C(7)
      FIJK(1,2,3)=0.0D+00
      FIJK(1,3,1)=FIJK(1,1,3)
      FIJK(1,3,2)=FIJK(1,2,3)
      FIJK(1,3,3)=C(8)
      FIJK(2,1,1)=FIJK(1,1,2)
      FIJK(2,1,2)=FIJK(1,2,2)
      FIJK(2,1,3)=FIJK(1,2,3)
      FIJK(2,2,1)=FIJK(1,2,2)
      FIJK(2,2,2)=C(9)
      FIJK(2,2,3)=0.0D+00
      FIJK(2,3,1)=FIJK(1,2,3)
      FIJK(2,3,2)=FIJK(2,2,3)
      FIJK(2,3,3)=C(10)
      FIJK(3,1,1)=FIJK(1,1,3)
      FIJK(3,1,2)=FIJK(1,2,3)
      FIJK(3,1,3)=FIJK(1,3,3)
      FIJK(3,2,1)=FIJK(1,2,3)
      FIJK(3,2,2)=FIJK(2,2,3)
      FIJK(3,2,3)=FIJK(2,3,3)
      FIJK(3,3,1)=FIJK(1,3,3)
      FIJK(3,3,2)=FIJK(2,3,3)
      FIJK(3,3,3)=0.0D+00
      
      V3=0.000000D+00
      DO I=1,N
        DO J=1,N
          DO K=1,N
          V3=V3+(1.0D+00/6.0D+00)*FIJK(I,J,K)*S(I)*S(J)*S(K)
          !PRINT*, FIJK(I,J,K)*S(I)*S(J)*S(K),I,J,K,V3
          END DO
        END DO    
      END DO
C
C QUARTIC PART
C
      DO I=1,N      
        DO J=1,N
          DO K=1,N
            DO L=1,N
              FIJKL(I,J,K,L)=0.0D+00
              !WRITE(*,FMT='(A6,I1,A,I1,A,I1,A,I1,A2)') 
!     *        "FIJKL(",I,",",J,",",K,",",L,")="
            END DO
          END DO
        END DO
      END DO
      
      FIJKL(1,1,1,1)=C(11)
      FIJKL(1,1,1,2)=C(12)
      FIJKL(1,1,1,3)=0.0D+00
      FIJKL(1,1,2,1)=FIJKL(1,1,1,2)
      FIJKL(1,1,2,2)=C(13)
      FIJKL(1,1,2,3)=0.0D+00
      FIJKL(1,1,3,1)=FIJKL(1,1,1,3)
      FIJKL(1,1,3,2)=FIJKL(1,1,2,3)
      FIJKL(1,1,3,3)=C(14)
      FIJKL(1,2,1,1)=FIJKL(1,1,1,2)
      FIJKL(1,2,1,2)=FIJKL(1,1,2,2)
      FIJKL(1,2,1,3)=FIJKL(1,1,2,3)
      FIJKL(1,2,2,1)=FIJKL(1,1,2,2)
      FIJKL(1,2,2,2)=C(15)
      FIJKL(1,2,2,3)=0.0D+00
      FIJKL(1,2,3,1)=FIJKL(1,1,2,3)
      FIJKL(1,2,3,2)=FIJKL(1,2,2,3)
      FIJKL(1,2,3,3)=C(16)
      FIJKL(1,3,1,1)=FIJKL(1,1,1,3)
      FIJKL(1,3,1,2)=FIJKL(1,1,2,3)
      FIJKL(1,3,1,3)=FIJKL(1,1,3,3)
      FIJKL(1,3,2,1)=FIJKL(1,1,2,3)
      FIJKL(1,3,2,2)=FIJKL(1,2,2,3)
      FIJKL(1,3,2,3)=FIJKL(1,2,3,3)
      FIJKL(1,3,3,1)=FIJKL(1,1,3,3)
      FIJKL(1,3,3,2)=FIJKL(1,2,3,3)
      FIJKL(1,3,3,3)=0.0D+00
      FIJKL(2,1,1,1)=FIJKL(1,1,1,2)
      FIJKL(2,1,1,2)=FIJKL(1,1,2,2)
      FIJKL(2,1,1,3)=FIJKL(1,1,2,3)
      FIJKL(2,1,2,1)=FIJKL(1,1,2,2)
      FIJKL(2,1,2,2)=FIJKL(1,2,2,2)
      FIJKL(2,1,2,3)=FIJKL(1,2,2,3)
      FIJKL(2,1,3,1)=FIJKL(1,1,2,3)
      FIJKL(2,1,3,2)=FIJKL(1,2,2,3)
      FIJKL(2,1,3,3)=FIJKL(1,2,3,3)
      FIJKL(2,2,1,1)=FIJKL(1,1,2,2)
      FIJKL(2,2,1,2)=FIJKL(1,2,2,2)
      FIJKL(2,2,1,3)=FIJKL(1,2,2,3)
      FIJKL(2,2,2,1)=FIJKL(1,2,2,2)
      FIJKL(2,2,2,2)=C(17)
      FIJKL(2,2,2,3)=0.0D+00
      FIJKL(2,2,3,1)=FIJKL(1,2,2,3)
      FIJKL(2,2,3,2)=FIJKL(2,2,2,3)
      FIJKL(2,2,3,3)=C(18)
      FIJKL(2,3,1,1)=FIJKL(1,1,2,3)
      FIJKL(2,3,1,2)=FIJKL(1,2,2,3)
      FIJKL(2,3,1,3)=FIJKL(1,2,3,3)
      FIJKL(2,3,2,1)=FIJKL(1,2,2,3)
      FIJKL(2,3,2,2)=FIJKL(2,2,2,3)
      FIJKL(2,3,2,3)=FIJKL(2,2,3,3)
      FIJKL(2,3,3,1)=FIJKL(1,2,3,3)
      FIJKL(2,3,3,2)=FIJKL(2,2,3,3)
      FIJKL(2,3,3,3)=0.0D+00
      FIJKL(3,1,1,1)=FIJKL(1,1,1,3)
      FIJKL(3,1,1,2)=FIJKL(1,1,2,3)
      FIJKL(3,1,1,3)=FIJKL(1,1,3,3)
      FIJKL(3,1,2,1)=FIJKL(1,1,2,3)
      FIJKL(3,1,2,2)=FIJKL(1,2,2,3)
      FIJKL(3,1,2,3)=FIJKL(1,2,3,3)
      FIJKL(3,1,3,1)=FIJKL(1,1,3,3)
      FIJKL(3,1,3,2)=FIJKL(1,2,3,3)
      FIJKL(3,1,3,3)=FIJKL(1,3,3,3)
      FIJKL(3,2,1,1)=FIJKL(1,1,2,3)
      FIJKL(3,2,1,2)=FIJKL(1,2,2,3)
      FIJKL(3,2,1,3)=FIJKL(1,2,3,3)
      FIJKL(3,2,2,1)=FIJKL(1,2,2,3)
      FIJKL(3,2,2,2)=FIJKL(2,2,2,3)
      FIJKL(3,2,2,3)=FIJKL(2,2,3,3)
      FIJKL(3,2,3,1)=FIJKL(1,2,3,3)
      FIJKL(3,2,3,2)=FIJKL(2,2,3,3)
      FIJKL(3,2,3,3)=FIJKL(2,3,3,3)
      FIJKL(3,3,1,1)=FIJKL(1,1,3,3)
      FIJKL(3,3,1,2)=FIJKL(1,2,3,3)
      FIJKL(3,3,1,3)=FIJKL(1,3,3,3)
      FIJKL(3,3,2,1)=FIJKL(1,2,3,3)
      FIJKL(3,3,2,2)=FIJKL(2,2,3,3)
      FIJKL(3,3,2,3)=FIJKL(2,3,3,3)
      FIJKL(3,3,3,1)=FIJKL(1,3,3,3)
      FIJKL(3,3,3,2)=FIJKL(2,3,3,3)
      FIJKL(3,3,3,3)=C(19)      

      V4=0.000000D+00
      DO I=1,N
        DO J=1,N
          DO K=1,N
            DO L=1,N
              V4=V4+(1.0D+00/24.0D+00)*
     *              FIJKL(I,J,K,L)*S(I)*S(J)*S(K)*S(L)
            END DO
          END DO
        END DO    
      END DO
 
      !POTEH=V0+(V2+V3+V4)*FAC
      
      POTEH=(V2+V3+V4)*FAC   
      
!      DO I=1,N     
!        DO J=I,N
!          WRITE(*,FMT='(2I5,10X,F20.10)') J,I,FIJ(J,I)
!        END DO
!      END DO
!      WRITE(*,*) "   0"
!      DO I=1,N      
!        DO J=I,N
!          DO K=J,N
!            WRITE(*,FMT='(3I5,5X,F20.10)') K,J,I,FIJK(K,J,I)
!          END DO
!        END DO
!      END DO      
!      WRITE(*,*) "   0"
!      DO I=1,N      
!        DO J=I,N
!          DO K=J,N
!            DO L=K,N
!              WRITE(*,FMT='(4I5,F20.10)') L,K,J,I,FIJKL(L,K,J,I)
!            END DO
!          END DO
!        END DO
!      END DO
!      WRITE(*,*) "   0"       
           
      RETURN
      END         
