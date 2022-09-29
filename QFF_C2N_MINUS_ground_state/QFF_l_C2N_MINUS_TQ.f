C      PROGRAM MAIN
C      IMPLICIT NONE
C      DOUBLE PRECISION, DIMENSION(3) :: R
C      DOUBLE PRECISION :: V,ANGS2AU
C      
C      ANGS2AU=0.529177249D+00
C      
C      R(1)=1.270543D+00 
C      R(2)=1.843761D+00
C      R(3)=1.843761D+00
C      
C      CALL QFF(R,V) 
C      
C      !PRINT*, V
C      
C      END PROGRAM MAIN
C
C##############################################################################
C QFF MODEL FUNCTION FOR LINEAR AB2 MOLECULES
C##############################################################################
C R1=A-B IN ANGS 
C R2=A-B IN ANGS
C R3=B-B IN ANGS
C V0 IS IN Eh
C RREF IS IN ANGS
C THETAREF IS IN DEGS
C POTENTIAL IN Eh
C##############################################################################

      SUBROUTINE QFF(RAU,POTEH)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (N=3,NN=4)
      DIMENSION RAU(N),R(N),S(NN)
      DIMENSION FIJ(NN,NN)
      DIMENSION FIJK(NN,NN,NN)
      DIMENSION FIJKL(NN,NN,NN,NN)
      DIMENSION C(23)

       C(  1)=  0.123643019994067238D+02
       C(  2)=  0.208392850543503405D+01
       C(  3)=  0.549136992698635673D+01
       C(  4)=  0.414145228140827604D+00
       C(  5)= -0.943113195116454222D+02
       C(  6)=  0.247607570843894287D+01
       C(  7)= -0.690869559052359516D+01
       C(  8)= -0.874918702192963971D+00
       C(  9)= -0.386440217450307699D+02
       C( 10)= -0.484421565385030561D+00
       C( 11)=  0.497668308883988800D+03
       C( 12)=  0.211112698245120782D+02
       C( 13)= -0.177900240512611951D+02
       C( 14)= -0.302205981438632920D+00
       C( 15)=  0.163592106504568200D+02
       C( 16)=  0.242813903680661181D+01
       C( 17)=  0.207871439751362857D+03
       C( 18)= -0.746732717348490227D+00
       C( 19)=  0.216428977765001473D+01
       C( 20)=  0.120370612208376171D+01
       C( 21)=  0.135938671502126041D+01
       C( 22)=  0.180000000000000000D+03
       C( 23)= -0.130914535769027225D+03
                 
      RREF1=C(20)
      RREF2=C(21)
      THETAREF=C(22)
      V0=C(23)
       
      POTEH=0.000000D+00     
      PI=ACOS(-1.0D+00)
      ANGS2AU=0.529177249D+00
      FAC=0.22937104486906D+00 ! CONVERSION FROM aJ TO Eh
      
      R(1)=RAU(1)
      R(2)=RAU(2)
      R(3)=RAU(3)
      
      THETAREF=THETAREF*PI/180.0D+00
      
      COSTHETA=(R(2)**2-R(1)**2-R(3)**2)/(-2.0D+00*R(1)*R(3))
      IF (COSTHETA>1.00D+00) THEN 
        COSTHETA=1.00D+00
      ELSE IF (COSTHETA<-1.00D+00) THEN 
        COSTHETA=-1.00D+00
      END IF 

      THETA=ACOS(COSTHETA)

      S(1)=(R(1)-RREF1)
      S(2)=(R(3)-RREF2)
      S(3)=SIN(THETA-THETAREF) 
      S(4)=SIN(THETA-THETAREF)      
C
C QUADRATIC PART
C      
      DO I=1,NN     
        DO J=1,NN
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
C
C ADD LINEAR CONSTRAINTS
C
      FIJ(4,4)=FIJ(3,3)      
      
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
      DO I=1,NN      
        DO J=1,NN
          DO K=1,NN
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
C
C ADD LINEAR CONSTRAINTS
C      
      FIJK(1,4,4)=FIJK(1,3,3)
      FIJK(4,1,4)=FIJK(1,3,3)
      FIJK(4,4,1)=FIJK(1,3,3)
C      
      FIJK(2,4,4)=FIJK(2,3,3)
      FIJK(4,2,4)=FIJK(2,3,3)
      FIJK(4,4,2)=FIJK(2,3,3)      
      
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
      DO I=1,NN      
        DO J=1,NN
          DO K=1,NN
            DO L=1,NN
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
C
C ADD LINEAR CONSTRAINTS
C
      FIJKL(1,1,4,4)=FIJKL(1,1,3,3)
      FIJKL(1,4,1,4)=FIJKL(1,1,3,3)
      FIJKL(1,4,4,1)=FIJKL(1,1,3,3)
      FIJKL(4,1,1,4)=FIJKL(1,1,3,3)
      FIJKL(4,1,4,1)=FIJKL(1,1,3,3)
      FIJKL(4,4,1,1)=FIJKL(1,1,3,3)
C      
      FIJKL(1,2,4,4)=FIJKL(1,2,3,3)
      FIJKL(1,4,2,4)=FIJKL(1,2,3,3)
      FIJKL(1,4,4,2)=FIJKL(1,2,3,3)
      FIJKL(2,1,4,4)=FIJKL(1,2,3,3)
      FIJKL(2,4,1,4)=FIJKL(1,2,3,3)
      FIJKL(2,4,4,1)=FIJKL(1,2,3,3)
      FIJKL(4,1,2,4)=FIJKL(1,2,3,3)
      FIJKL(4,1,4,2)=FIJKL(1,2,3,3)
      FIJKL(4,2,1,4)=FIJKL(1,2,3,3)
      FIJKL(4,2,4,1)=FIJKL(1,2,3,3)
      FIJKL(4,4,1,2)=FIJKL(1,2,3,3)
      FIJKL(4,4,2,1)=FIJKL(1,2,3,3)      
C
      FIJKL(2,2,4,4)=FIJKL(2,2,3,3)
      FIJKL(2,4,2,4)=FIJKL(2,2,3,3)
      FIJKL(2,4,4,2)=FIJKL(2,2,3,3)
      FIJKL(4,2,2,4)=FIJKL(2,2,3,3)
      FIJKL(4,2,4,2)=FIJKL(2,2,3,3)
      FIJKL(4,4,2,2)=FIJKL(2,2,3,3)
C      
      FIJKL(4,4,4,4)=FIJKL(3,3,3,3)
C     
      FIJKL(3,3,4,4)=(FIJKL(3,3,3,3)+
     *                4.0D+00*FIJ(3,3))/3.0D+00
      FIJKL(3,4,3,4)=FIJKL(3,3,4,4)
      FIJKL(3,4,4,3)=FIJKL(3,3,4,4)
      FIJKL(4,3,3,4)=FIJKL(3,3,4,4)
      FIJKL(4,3,4,3)=FIJKL(3,3,4,4)
      FIJKL(4,4,3,3)=FIJKL(3,3,4,4)
      
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
      
C      DO I=1,NN     
C        DO J=I,NN
C          WRITE(*,FMT='(2I5,10X,F20.10)') J,I,FIJ(J,I)
C        END DO
C      END DO
C      WRITE(*,*) "   0"
C      DO I=1,NN      
C        DO J=I,NN
C          DO K=J,NN
C            WRITE(*,FMT='(3I5,5X,F20.10)') K,J,I,FIJK(K,J,I)
C          END DO
C        END DO
C      END DO      
C      WRITE(*,*) "   0"
C      DO I=1,NN      
C        DO J=I,NN
C          DO K=J,NN
C            DO L=K,NN
C              WRITE(*,FMT='(4I5,F20.10)') L,K,J,I,FIJKL(L,K,J,I)
C            END DO
C          END DO
C        END DO
C      END DO
C      WRITE(*,*) "   0"       
           
      RETURN
      END         
