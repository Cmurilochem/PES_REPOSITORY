      PROGRAM TEST
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(3) :: R,DVDR 
      DOUBLE PRECISION :: V
      LOGICAL :: GRAD
      DOUBLE PRECISION, DIMENSION(9) :: X,DVDX
C
C     SET GRAD=.TRUE. FOR ANALYTIC GRADIENTS
C     WITH RESPECT TO INTERNAL (R_1, R_2, R_3) COORDINATES
C
      GRAD=.TRUE.
C
C     C2V CYCLIC GLOBAL MINIMUM 
C 
      R(1)=2.40133452D+00
      R(2)=3.46798418D+00
      R(3)=3.46798418D+00

      CALL SiC2_PES(R,V,GRAD,DVDR)

      WRITE(*,*)"MIN  C2V ENERGY (IN HARTREE):", V,DVDR
C
C     CINFV TS
C 
      R(1)=2.43009593D+00
      R(2)=5.63243414D+00
      R(3)=3.20233821D+00

      CALL SiC2_PES(R,V,GRAD,DVDR)

      WRITE(*,*)"TS CINFV ENERGY (IN HARTREE):", V
C
C     DINFH TS
C 
      R(1)=6.81873273D+00
      R(2)=3.40936637D+00
      R(3)=3.40936637D+00

      CALL SiC2_PES(R,V,GRAD,DVDR)

      WRITE(*,*)"TS DINFH ENERGY (IN HARTREE):", V
C
C     TO CALCULATE ANALYTIC GRADIENTS
C     IN CARTESIAN COORDINATES USE:
C
      X(1)=0.00000000D+00 !C1-X
      X(2)=0.00000000D+00 !C1-Y
      X(3)=0.00000000D+00 !C1-Z 
      X(4)=0.00000000D+00 !C2-X
      X(5)=0.00000000D+00 !C2-Y
      X(6)=6.81873273D+00 !C2-Z
      X(7)=0.00000000D+00 !Si-X
      X(8)=0.00000000D+00 !Si-Y
      X(9)=3.40936637D+00 !Si-Z
            
      CALL CARTDERPES(X,V,DVDX)
      
      WRITE(*,*)"TS DINFH ENERGY (IN HARTREE):", V,DVDX

      END PROGRAM TEST
