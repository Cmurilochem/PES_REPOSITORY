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
      PARAMETER (N=3,NN=3)
      DIMENSION RAU(N),R(N),S(N)
      DIMENSION FIJ(N,N)
      DIMENSION FIJK(N,N,N)
      DIMENSION FIJKL(N,N,N,N)
      DIMENSION C(22)
      DIMENSION BLBAFIJ(N,N),BLBAFIJK(N,N,N),BLBAFIJKL(N,N,N,N)
      DIMENSION AMFIJ(N,N),AMFIJK(N,N,N),AMFIJKL(N,N,N,N)
      DIMENSION DSDs(NN),DDSDDs(NN,NN),DDDSDDDs(NN,NN,NN)
      DIMENSION D(NN,NN),DD(NN,NN,NN)
C
C FORCE CONSTANTS IN SYMMETRY COORDINATES
C
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
C
C DERIVATIVES OF THE MORSE-SINE VS INTERNAL VS SYMMETRY
C
      DO I=1,NN  
        DSDs(I)=0.0D+00
        DO J=1,NN
          DDSDDs(I,J)=0.0D+00
          D(I,J)=0.0D+00
          DO K=1,NN
            DDDSDDDs(I,J,K)=0.0D+00 
            DD(I,J,K)=0.0D+00 
          END DO
        END DO
      END DO
C
C JACOBIANS 
C FIRST BLOCK - SYMMETRY TO BLBA
C
      D(1,1)=1.0D+00/SQRT(2.0D+00)
      D(1,3)=1.0D+00/SQRT(2.0D+00)
      D(2,2)=1.0D+00
      D(3,1)=1.0D+00/SQRT(2.0D+00)
      D(3,3)=-1.0D+00/SQRT(2.0D+00)      
C
C QUADRATIC PART
C      
      DO I=1,N     
        DO J=1,N
          FIJ(I,J)=0.0D+00
          BLBAFIJ(I,J)=0.0D+00
          AMFIJ(I,J)=0.0D+00
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
C
C CUBIC PART
C
      DO I=1,N      
        DO J=1,N
          DO K=1,N
          FIJK(I,J,K)=0.0D+00
          BLBAFIJK(I,J,K)=0.0D+00
          AMFIJK(I,J,K)=0.0D+00
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
C
C CONVERTING QUADRATIC AND CUBIC FORCE CONSTANTS FROM
C SYMMETRY COORDINATES TO BLBA 
C
      CALL CONV_SYM2BLBA_QUADRATIC(D,FIJ,BLBAFIJ)      
      CALL CONV_SYM2BLBA_CUBIC(D,FIJK,BLBAFIJK)
      
      SCL=2.05D+00
      
      ALPHA1=-SCL*(BLBAFIJK(1,1,1))/(3.0D+00*BLBAFIJ(1,1))  
      ALPHA2=-SCL*(BLBAFIJK(3,3,3))/(3.0D+00*BLBAFIJ(3,3)) 
C      
C JACOBIANS
C SECOND BLOCK - BLBA TO MORSE
C
      DSDs(1)=ALPHA1 
      DSDs(2)=-SIN(THETAREF)
      DSDs(3)=ALPHA2
      DDSDDs(1,1)=-(ALPHA1**2)
      DDSDDs(2,2)=-COS(THETAREF)
      DDSDDs(3,3)=-(ALPHA2**2)
      DDDSDDDs(1,1,1)=(ALPHA1**3)
      DDDSDDDs(2,2,2)=SIN(THETAREF)
      DDDSDDDs(3,3,3)=(ALPHA2**3)     
      
      S(1)=1.0D+00-EXP(-ALPHA1*(R(2)-RREF)) !(R(2)-RREF)!(1.0D+00/SQRT(2.0D+00))*((R(2)-RREF)+(R(3)-RREF))
      S(2)=COS(THETA)-COS(THETAREF) !THETA-THETAREF
      S(3)=1.0D+00-EXP(-ALPHA2*(R(3)-RREF)) !(R(3)-RREF)!(1.0D+00/SQRT(2.0D+00))*((R(2)-RREF)-(R(3)-RREF)) 
C
C CONVERTING QUADRATIC AND CUBIC FORCE CONSTANTS FROM
C BLBA COORDINATES TO MORSE-COSINE 
C      
      CALL CONV_BLBA2MORSE_QUADRATIC(ALPHA1,ALPHA2,THETAREF,
     *                    DSDs,BLBAFIJ,AMFIJ)
     
      CALL CONV_BLBA2MORSE_CUBIC(ALPHA1,ALPHA2,THETAREF,
     *                DSDs,DDSDDs,AMFIJ,BLBAFIJK,AMFIJK)     
      
      V2=0.000000D+00
      DO I=1,N
        DO J=1,N
          V2=V2+(1.0D+00/2.0D+00)*AMFIJ(I,J)*S(I)*S(J)
          !PRINT*, FIJ(I,J),S(I),S(J),I,J,V2
        END DO    
      END DO      
      
      V3=0.000000D+00
      DO I=1,N
        DO J=1,N
          DO K=1,N
          V3=V3+(1.0D+00/6.0D+00)*AMFIJK(I,J,K)*S(I)*S(J)*S(K)
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
              BLBAFIJKL(I,J,K,L)=0.0D+00
              AMFIJKL(I,J,K,L)=0.0D+00
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
C
C CONVERTING QUARTIC FORCE CONSTANTS FROM
C SYMMETRY COORDINATES TO BLBA 
C      
      CALL CONV_SYM2BLBA_QUARTIC(D,FIJKL,BLBAFIJKL)
C
C CONVERTING QUARTIC FORCE CONSTANTS FROM
C BLBA COORDINATES TO MORSE-COSINE
C      
      CALL CONV_BLBA2MORSE_QUARTIC(ALPHA1,ALPHA2,THETAREF,
     *     DSDs,DDSDDs,DDDSDDDs,AMFIJ,AMFIJK,BLBAFIJKL,AMFIJKL)
     
      V4=0.000000D+00
      DO I=1,N
        DO J=1,N
          DO K=1,N
            DO L=1,N
              V4=V4+(1.0D+00/24.0D+00)*
     *        AMFIJKL(I,J,K,L)*S(I)*S(J)*S(K)*S(L)
            END DO
          END DO
        END DO    
      END DO
 
      !POTEH=V0+(V2+V3+V4)*FAC
      
      POTEH=(V2+V3+V4)*FAC  
      
!      DO I=1,N     
!        DO J=I,N
!          WRITE(*,FMT='(2I5,10X,F20.10)') J,I,AMFIJ(J,I)
!        END DO
!      END DO
!      WRITE(*,*) "   0"
!      DO I=1,N      
!        DO J=I,N
!          DO K=J,N
!            WRITE(*,FMT='(3I5,5X,F20.10)') K,J,I,AMFIJK(K,J,I)
!          END DO
!        END DO
!      END DO      
!      WRITE(*,*) "   0"
!      DO I=1,N      
!        DO J=I,N
!          DO K=J,N
!            DO L=K,N
!              WRITE(*,FMT='(4I5,F20.10)') L,K,J,I,AMFIJKL(L,K,J,I)
!            END DO
!          END DO
!        END DO
!      END DO
!      WRITE(*,*) "   0"       
           
      RETURN
      END

C##############################################################################      
C CONVERSION OF THE FORCE CONSTANT MATRIX FROM SYMMETRY TO BLBA
C##############################################################################
      SUBROUTINE CONV_SYM2BLBA_QUADRATIC(D,FIJ,FFIJ)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (N=3,NN=3)
      DIMENSION FIJ(NN,NN),FFIJ(NN,NN)
      DIMENSION D(NN,NN)
      DO I=1,NN     
        DO J=1,NN
          CALL MULT_QUADRATIC(I,J,D,FIJ,RES)
          FFIJ(I,J)=RES
          !WRITE(*,FMT='(2I5,10X,F20.10)') I,J,FFIJ(I,J)
        END DO
      END DO     
      RETURN
      END
C##############################################################################
C##############################################################################
      SUBROUTINE CONV_SYM2BLBA_CUBIC(D,FIJK,FFIJK)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (N=3,NN=3)
      DIMENSION FIJK(NN,NN,NN),FFIJK(NN,NN,NN)
      DIMENSION D(NN,NN)
      DO I=1,NN      
        DO J=1,NN
          DO K=1,NN
            CALL MULT_CUBIC(I,J,K,D,FIJK,RES)
            FFIJK(I,J,K)=RES
            !WRITE(*,FMT='(3I5,5X,F20.10)') I,J,K,FFIJK(I,J,K)
          END DO
        END DO
      END DO     
      RETURN
      END
C##############################################################################
C##############################################################################
      SUBROUTINE CONV_SYM2BLBA_QUARTIC(D,FIJKL,FFIJKL)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (N=3,NN=3)
      DIMENSION FIJKL(NN,NN,NN,NN),FFIJKL(NN,NN,NN,NN)
      DIMENSION D(NN,NN)
      DO I=1,NN      
        DO J=1,NN
          DO K=1,NN
            DO L=1,NN
              CALL MULT_QUARTIC(I,J,K,L,D,FIJKL,RES)
              FFIJKL(I,J,K,L)=RES
              !WRITE(*,FMT='(4I5,5X,F20.10)') I,J,K,L,FFIJKL(I,J,K,L)
            END DO
          END DO
        END DO
      END DO    
      RETURN
      END      
C##############################################################################
C##############################################################################
      SUBROUTINE MULT_QUADRATIC(IFIX,JFIX,D,FIJ,RES)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (N=3,NN=3) 
      DIMENSION FIJ(NN,NN),D(NN,NN)
      SUM=0.0D+00
      DO I=1,NN     
        DO J=1,NN
          SUM=SUM+D(I,IFIX)*FIJ(I,J)*D(J,JFIX)
        END DO
      END DO
      RES=SUM
      RETURN
      END
      
      SUBROUTINE MULT_CUBIC(IFIX,JFIX,KFIX,D,FIJK,RES)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (N=3,NN=3) 
      DIMENSION FIJK(NN,NN,NN),D(NN,NN)
      SUM=0.0D+00
      DO I=1,NN     
        DO J=1,NN
          DO K=1,NN
            SUM=SUM+D(I,IFIX)*D(J,JFIX)*D(K,KFIX)*FIJK(I,J,K)
          END DO
        END DO
      END DO
      RES=SUM
      RETURN
      END
      
      SUBROUTINE MULT_QUARTIC(IFIX,JFIX,KFIX,LFIX,D,FIJKL,RES)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (N=3,NN=3) 
      DIMENSION FIJKL(NN,NN,NN,NN),D(NN,NN)
      SUM=0.0D+00
      DO I=1,NN      
        DO J=1,NN
          DO K=1,NN
            DO L=1,NN
              SUM=SUM+D(I,IFIX)*D(J,JFIX)*
     *        D(K,KFIX)*D(L,LFIX)*FIJKL(I,J,K,L)
            END DO
          END DO
        END DO
      END DO
      RES=SUM
      RETURN
      END      
C##############################################################################
C CONVERSION OF THE FORCE CONSTANT MATRIX FROM INTERNAL TO MORSE-COSSINE
C##############################################################################      
      SUBROUTINE CONV_BLBA2MORSE_QUADRATIC(A1,A2,ANG0,DSDs,FIJ,FFIJ)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (N=3,NN=3)
      DIMENSION FIJ(NN,NN),FFIJ(NN,NN)
      DIMENSION DSDs(NN)           
      DO I=1,NN     
        DO J=1,NN
          FFIJ(I,J)=FIJ(I,J)/(DSDs(I)*DSDs(J))
          !WRITE(*,FMT='(2I5,10X,F20.10)') I,J,FFIJ(I,J)
        END DO
      END DO
C CHECK      
C      DO I=1,NN     
C        DO J=1,NN
C          WRITE(*,FMT='(2I5,10X,F20.10)') J,I,FFIJ(J,I)
C        END DO
C      END DO     
      RETURN 
      END
C##############################################################################
C##############################################################################      
      SUBROUTINE CONV_BLBA2MORSE_CUBIC(A1,A2,ANG0,
     * DSDs,DDSDDs,FFIJ,FIJK,FFIJK)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (N=3,NN=3) 
      DIMENSION FFIJ(NN,NN)
      DIMENSION FIJK(NN,NN,NN),FFIJK(NN,NN,NN)
      DIMENSION DSDs(NN),DDSDDs(NN,NN)
      DO I=1,NN      
        DO J=1,NN
          DO K=1,NN
            P1=-FFIJ(I,J)*DSDs(J)*DDSDDs(K,I)
            P2=-FFIJ(I,J)*DSDs(I)*DDSDDs(K,J)
            P3=-FFIJ(I,K)*DSDs(K)*DDSDDs(J,I)
            P4=FIJK(I,J,K)+P1+P2+P3
            P5=P4/(DSDs(I)*DSDs(J)*DSDs(K))
            FFIJK(I,J,K)=P5
            !WRITE(*,FMT='(3I5,5X,3F20.10)') I,J,K,P1,P2,FFIJK(I,J,K)
          END DO
        END DO
      END DO
C      DO I=1,NN      
C        DO J=1,NN
C          DO K=1,NN
C            WRITE(*,FMT='(3I5,5X,F20.10)') I,J,K,FFIJK(I,J,K)
C          END DO
C        END DO
C      END DO 
      RETURN 
      END
C##############################################################################
C##############################################################################      
      SUBROUTINE CONV_BLBA2MORSE_QUARTIC(A1,A2,ANG0,
     * DSDs,DDSDDs,DDDSDDDs,FFIJ,FFIJK,FIJKL,FFIJKL)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (N=3,NN=3)
      DIMENSION FFIJ(NN,NN)
      DIMENSION FFIJK(NN,NN,NN)
      DIMENSION FIJKL(NN,NN,NN,NN),FFIJKL(NN,NN,NN,NN)
      DIMENSION DSDs(NN),DDSDDs(NN,NN),DDDSDDDs(NN,NN,NN)      
      DO I=1,NN      
        DO J=1,NN
          DO K=1,NN
            DO L=1,NN
              P1=-FFIJK(I,J,K)*DSDs(J)*DSDs(K)*DDSDDs(L,I)
              P2=-FFIJK(I,J,K)*DSDs(I)*DSDs(K)*DDSDDs(L,J)
              P3=-FFIJK(I,J,K)*DSDs(I)*DSDs(J)*DDSDDs(L,K) 
              P4=-FFIJK(I,J,L)*DSDs(J)*DSDs(L)*DDSDDs(K,I)
              P5=-FFIJ(I,J)*DDSDDs(K,I)*DDSDDs(L,J)
              P6=-FFIJ(I,J)*DSDs(J)*DDDSDDDs(L,K,I)
              P7=-FFIJK(I,K,L)*DSDs(K)*DSDs(L)*DDSDDs(J,I)
              P8=-FFIJ(I,K)*DDSDDs(J,I)*DDSDDs(L,K)
              P9=-FFIJ(I,K)*DSDs(K)*DDDSDDDs(L,J,I)
              P10=-FFIJK(I,J,L)*DSDs(I)*DSDs(L)*DDSDDs(K,J)
              P11=-FFIJ(I,J)*DDSDDs(K,J)*DDSDDs(L,I)
              P12=-FFIJ(I,J)*DSDs(I)*DDDSDDDs(L,K,J)
              P13=-FFIJ(I,L)*DSDs(L)*DDDSDDDs(K,J,I)
              P14=FIJKL(I,J,K,L)+P1+P2+P3+P4+P5+P6+P7+P8+P9+P10+
     *                           P11+P12+P13
              P15=P14/(DSDs(I)*DSDs(J)*DSDs(K)*DSDs(L))
              FFIJKL(I,J,K,L)=P15                       
              !WRITE(*,FMT='(4I5,5X,3F20.10)') I,J,K,L,FFIJKL(I,J,K,L)
            END DO
          END DO
        END DO
      END DO
C      DO I=1,NN      
C        DO J=I,NN
C          DO K=J,NN
C            DO L=K,NN
C              WRITE(*,FMT='(4I5,F20.10)') L,K,J,I,FFIJKL(L,K,J,I)
C            END DO
C          END DO
C        END DO
C      END DO      
      RETURN 
      END      
