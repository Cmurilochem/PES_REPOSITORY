      PROGRAM C3
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(3) :: R 
      DOUBLE PRECISION :: VC3
  
      R(1)=2.445D+00
      R(2)=4.890D+00
      R(3)=2.445D+00

      CALL POTC3SPEC(R,VC3)
 
C      WRITE(*,*)"COORDINATES (IN BOHR):", R

      WRITE(*,*)"ENERGY (IN HARTRE):", VC3

      END PROGRAM C3
