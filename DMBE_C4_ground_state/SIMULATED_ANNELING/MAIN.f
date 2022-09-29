      PROGRAM SA
      IMPLICIT NONE
      INTEGER :: I, J, ERR, NATOM, INTDIST
      INTEGER, PARAMETER :: NATOMMAX=100
      DOUBLE PRECISION, DIMENSION(3*NATOMMAX) :: X, R
      DOUBLE PRECISION, DIMENSION(3*NATOMMAX) :: XOPT, ROPT
      INTEGER :: NCARTCOORD, I1, I2, I3
      CHARACTER(1), DIMENSION(NATOMMAX) :: ATOM
      CHARACTER(5) :: POINTG
      DOUBLE PRECISION :: POT0, FOPT, BOHR2ANGS
      INTEGER :: IER

C INICIALIZING VECTORS

      DO I=1,3*NATOMMAX
       X(I)=0.0D+00
       R(I)=0.0D+00
      END DO

C###############################################################################

      OPEN(UNIT=1,FILE='mol_data.txt',IOSTAT=ERR,
     1 FORM='FORMATTED',STATUS='OLD',ACTION='READ')
        IF (ERR /= 0) THEN
          PRINT*, "PROBLEMS IN OPENING THE FILE MOL_DATA.TXT"
        END IF

      OPEN(UNIT=2,FILE='SOURCE/symm.dat',IOSTAT=ERR,
     1 FORM='FORMATTED',STATUS='REPLACE',ACTION='WRITE')
        IF (ERR /= 0) THEN
          PRINT*, "PROBLEMS IN OPENING THE FILE SYMM.DAT"
        END IF

      WRITE(*,50)
      WRITE(*,75)
      WRITE(*,100)

C###############################################################################

C READ # OF ATOMS

      READ(1,*) NATOM

C CALCULATING # OF INTERPARTICLE DISTANCES

      INTDIST=(NATOM*(NATOM-1))/2 

C # OF CARTESIAN COORDINATES

      NCARTCOORD=3*NATOM 

      WRITE(*,101) NATOM
      WRITE(*,100)
      WRITE(*,102) INTDIST
      WRITE(*,100)

C READING THE ATOMIC SYMBOLS

        DO I=1,NATOM

          READ(1,*) ATOM(I) 

        END DO

C READING THE CARTESIAN COORDINATES

      READ(1,*) (X(I),I=1,NCARTCOORD)

      CLOSE(UNIT=1)

C DETERMINING THE MOLECULAR POINT GROUP      

      WRITE(2,*) NATOM 

      DO I=1,NATOM

        I1=3*(I-1)+1

        I2=3*(I-1)+2

        I3=3*(I-1)+3

        WRITE(2,*) ICHAR(ATOM(I)),(X(3*(I-1)+J),J=1,3)
   
      END DO
     
      CLOSE(UNIT=2)

C EXECUTING THE PROGRAM FOR SYMMETRY DETERMINATION

      CALL SYSTEM("gcc SOURCE/symm.c -o SOURCE/symm.x -lm")

      CALL SYSTEM("./SOURCE/symm.x < SOURCE/symm.dat > SOURCE/symm.out")

      CALL SYSTEM("mv symm.res SOURCE/.")

C###############################################################################

      OPEN(UNIT=3,FILE='SOURCE/symm.res',IOSTAT=ERR,
     1 FORM='FORMATTED',STATUS='OLD',ACTION='READ')
        IF (ERR /= 0) THEN
          PRINT*, "PROBLEMS IN OPENING THE FILE SYMM.RES"
        END IF

C###############################################################################

C DETERMINING THE MOLECULAR POINT GROUP OF THE INITIAL CONFIGURATION

      READ(3,*,IOSTAT=ERR) POINTG

        IF ( ERR /= 0 ) THEN

        POINTG="NAN"

        END IF

      CLOSE(UNIT=3)

      WRITE(*,103) POINTG
      WRITE(*,100)

      WRITE(*,104)
      WRITE(*,100)

      CALL FCN(NCARTCOORD,X,POT0)

      WRITE(*,105) POT0
      WRITE(*,100)

C############################################################################### STARTING OPTIMIZATION

      CALL SIMANN(NCARTCOORD,X,XOPT,FOPT,IER)

C###############################################################################

      BOHR2ANGS=0.529177208D+00

      IF (IER .EQ. 0) THEN 

      WRITE(*,50)
      WRITE(*,106)
      WRITE(*,50)
      WRITE(*,107)
      WRITE(*,100)

C############################################################################### WRITTING THE FINAL STRUCTURE TO A FILE
     
      OPEN(UNIT=4,FILE='min.molden',IOSTAT=ERR,
     1 FORM='FORMATTED',STATUS='REPLACE',ACTION='WRITE')
        IF (ERR /= 0) THEN
          PRINT*, "PROBLEMS IN OPENING THE FILE MIN.MOLDEN"
        END IF

C############################################################################### WRITTING THE FINAL STRUCTURE TO A FILE

      WRITE(4,FMT="(I3)") NATOM
      WRITE(4,FMT="(A)")

      DO I=1,NATOM

        WRITE(4,*) ATOM(I),(XOPT(3*(I-1)+J)*BOHR2ANGS,J=1,3)
   
      END DO

      CLOSE(UNIT=4)

      END IF 

C############################################################################### FORMATS

   50 FORMAT ("****************************************************")
   75 FORMAT (6X,"SIMULATED ANNEALING OPTIMIZATION PROGRAM")
  100 FORMAT ("****************************************************",/)
  101 FORMAT ("1). NUMBER OF ATOMS IN THE MOLECULE",2X,I3,/)
  102 FORMAT ("2). INTERNAL DEGREES OF FREEDOM",2X,I3,/)
  103 FORMAT ("3). SYMMETRY OF THE MOLECULAR GEOMETRY",2X,A5,/)
  104 FORMAT (X,"THE CARTEZIAN/INTERNAL COORDINATES ARE IN BOHRS",/)
  105 FORMAT (X,"THE ENERGY IS (IN HARTREES):",F15.10/)
  106 FORMAT (12X,"OPTIMIZATION COMPLETED")
  107 FORMAT (X,"WRITING THE OPT GEOMETRY IN THE FILE MIN.MOLDEN")
  108 FORMAT ("5). SYMMETRY OF THE OPT MOLECULAR GEOMETRY",2X,A5,/)

C###############################################################################

      END PROGRAM

C###############################################################################

      SUBROUTINE CARTINT(X,CINT,NAT)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NA=100,NT=3*NA,DEGREE=1.0)
      DIMENSION GEO(3,NA),COORD(3,NA)
      DIMENSION X(NT),CINT(NT)
      COMMON/CONECT/NNA(NA),NNB(NA),NNC(NA)
      DO I=1,NAT
        COORD(1,I)=X(3*(I-1)+1)
        COORD(2,I)=X(3*(I-1)+2)
        COORD(3,I)=X(3*(I-1)+3)
      ENDDO
      CALL XYZINT(COORD,NNA,NNB,NNC,DEGREE,GEO,NAT)
      CINT(1)=GEO(1,2)
      CINT(2)=GEO(1,3)
      CINT(3)=GEO(2,3)
      DO I=4,NAT
        CINT(3*(I-3)+1)=GEO(1,I)
        CINT(3*(I-3)+2)=GEO(2,I)
        CINT(3*(I-3)+3)=GEO(3,I)
      ENDDO
      END

C###############################################################################

      SUBROUTINE XYZINT(XYZ,NNA,NNB,NNC,DEGREE,GEO,NUMAT)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NA=100)
      DIMENSION XYZ(3,NA),NNA(NA),NNB(NA),NNC(NA),GEO(3,NA)
***********************************************************************
*
* XYZINT WORKS OUT THE INTERNAL COORDINATES OF A MOLECULE.
*        THE "RULES" FOR THE CONNECTIVITY ARE AS FOLLOWS:
*        ATOM I IS DEFINED AS BEING AT A DISTANCE FROM THE NEAREST
*        ATOM J, ATOM J ALREADY HAVING BEEN DEFINED.
*        ATOM I MAKES AN ANGLE WITH ATOM J AND THE ATOM K, WHICH HAS
*        ALREADY BEEN DEFINED, AND IS THE NEAREST ATOM TO J
*        ATOM I MAKES A DIHEDRAL ANGLE WITH ATOMS J, K, AND L. L HAVING
*        BEEN DEFINED AND IS THE NEAREST ATOM TO K, AND J, K AND L
*        HAVE A CONTAINED ANGLE IN THE RANGE 15 TO 165 DEGREES,
*        IF POSSIBLE.
*
*        IF(NA(2).EQ.1 THEN THE ORIGINAL CONNECTIVITY IS USED.
*
*        NOTE THAT GEO AND XYZ MUST NOT BE THE SAME IN THE CALL.
*
*   ON INPUT XYZ    = CARTESIAN ARRAY OF NUMAT ATOMS
*            DEGREE = 1 IF ANGLES ARE TO BE IN RADIANS
*            DEGREE = 57.29578 IF ANGLES ARE TO BE IN DEGREES
*
***********************************************************************
         DO 30 I=1,NUMAT
            NNA(I)=2
            NNB(I)=3
            NNC(I)=4
            IM1=I-1
            IF(IM1.EQ.0)GOTO 30
            SUM=1.D30
            DO 20 J=1,IM1
               R=(XYZ(1,I)-XYZ(1,J))**2+
     1          (XYZ(2,I)-XYZ(2,J))**2+
     2          (XYZ(3,I)-XYZ(3,J))**2
               IF(R.LT.SUM.AND.NNA(J).NE.J.AND.NNB(J).NE.J) THEN
                  SUM=R
                  K=J
               ENDIF
   20       CONTINUE
C
C   ATOM I IS NEAREST TO ATOM K
C
            NNA(I)=K
            IF(I.GT.2)NNB(I)=NNA(K)
            IF(I.GT.3)NNC(I)=NNB(K)
C
C   FIND ANY ATOM TO RELATE TO NA(I)
C
30    CONTINUE
      NNA(1)=0
      NNB(1)=0
      NNC(1)=0
      NNB(2)=0
      NNC(2)=0
      NNC(3)=0
      CALL XYZGEO(XYZ,NUMAT,NNA,NNB,NNC,DEGREE,GEO)
      RETURN
      END

C###############################################################################
      
      SUBROUTINE XYZGEO(XYZ,NUMAT,NNA,NNB,NNC,DEGREE,GEO)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NA=100)
      DIMENSION XYZ(3,NA), NNA(NA), NNB(NA), NNC(NA), GEO(3,NA)
***********************************************************************
*
*   XYZGEO CONVERTS COORDINATES FROM CARTESIAN TO INTERNAL.
*
*     ON INPUT XYZ  = ARRAY OF CARTESIAN COORDINATES
*              NUMAT= NUMBER OF ATOMS
*              NNA   = NUMBERS OF ATOM TO WHICH ATOMS ARE RELATED
*                     BY DISTANCE
*              NNB   = NUMBERS OF ATOM TO WHICH ATOMS ARE RELATED
*                     BY ANGLE
*              NNC   = NUMBERS OF ATOM TO WHICH ATOMS ARE RELATED
*                     BY DIHEDRAL
*
*    ON OUTPUT GEO  = INTERNAL COORDINATES IN ANGSTROMS, RADIANS,
*                     AND RADIANS
*
***********************************************************************
      DO 30 I=2,NUMAT
         J=NNA(I)
         K=NNB(I)
         L=NNC(I)
         IF(I.LT.3) GOTO 30
         II=I
         CALL BANGLE(XYZ,II,J,K,GEO(2,I))
         GEO(2,I)=GEO(2,I)*DEGREE
         IF(I.LT.4) GOTO 30
C
C   MAKE SURE DIHEDRAL IS MEANINGLFUL
C
         CALL BANGLE(XYZ,J,K,L,ANGL)
         TOL=0.2617994D0
         IF(ANGL.GT.3.1415926D0-TOL.OR.ANGL.LT.TOL)THEN
C
C  ANGLE IS UNSATISFACTORY, LET'S SEARCH FOR ANOTHER ATOM FOR
C  DEFINING THE DIHEDRAL.
   10       SUM=100.D0
            DO 20 I1=1,II-1
               R=(XYZ(1,I1)-XYZ(1,K))**2+
     1          (XYZ(2,I1)-XYZ(2,K))**2+
     2          (XYZ(3,I1)-XYZ(3,K))**2
               IF(R.LT.SUM.AND.I1.NE.J.AND.I1.NE.K) THEN
                  CALL BANGLE(XYZ,J,K,I1,ANGL)
                  IF(ANGL.LT.3.1415926D0-TOL.AND.ANGL.GT.TOL)THEN
                     SUM=R
                     L=I1
                     NNC(II)=L
                  ENDIF
               ENDIF
   20       CONTINUE
            IF(SUM.GT.99.D0.AND.TOL.GT.0.1D0)THEN
C
C ANYTHING WITHIN 5 DEGREES?
C
               TOL=0.087266D0
               GOTO 10
            ENDIF
         ENDIF
         CALL DIHED(XYZ,II,J,K,L,GEO(3,I))
         GEO(3,I)=GEO(3,I)*DEGREE
   30 GEO(1,I)= SQRT((XYZ(1,I)-XYZ(1,J))**2+
     1                   (XYZ(2,I)-XYZ(2,J))**2+
     2                   (XYZ(3,I)-XYZ(3,J))**2)
      GEO(1,1)=0.D0
      GEO(2,1)=0.D0
      GEO(3,1)=0.D0
      GEO(2,2)=0.D0
      GEO(3,2)=0.D0
      GEO(3,3)=0.D0
      RETURN
      END

C###############################################################################

      SUBROUTINE BANGLE(XYZ,I,J,K,ANGLE)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NA=100)
      DIMENSION XYZ(3,NA)
*********************************************************************
*
* BANGLE CALCULATES THE ANGLE BETWEEN ATOMS I,J, AND K. THE
*        CARTESIAN COORDINATES ARE IN XYZ.
*
*********************************************************************
      D2IJ = (XYZ(1,I)-XYZ(1,J))**2+
     1       (XYZ(2,I)-XYZ(2,J))**2+
     2       (XYZ(3,I)-XYZ(3,J))**2
      D2JK = (XYZ(1,J)-XYZ(1,K))**2+
     1       (XYZ(2,J)-XYZ(2,K))**2+
     2       (XYZ(3,J)-XYZ(3,K))**2
      D2IK = (XYZ(1,I)-XYZ(1,K))**2+
     1       (XYZ(2,I)-XYZ(2,K))**2+
     2       (XYZ(3,I)-XYZ(3,K))**2
      XY = SQRT(D2IJ*D2JK)
      TEMP = 0.5D0 * (D2IJ+D2JK-D2IK) / XY
      IF (TEMP .GT. 1.0D0) TEMP=1.0D0
      IF (TEMP .LT. -1.0D0) TEMP=-1.0D0
      ANGLE = ACOS( TEMP )
      RETURN
      END

C###############################################################################

      SUBROUTINE DIHED(XYZ,I,J,K,L,ANGLE)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NA=100)
      DIMENSION XYZ(3,NA)
*********************************************************************
*
*      DIHED CALCULATES THE DIHEDRAL ANGLE BETWEEN ATOMS I, J, K,
*            AND L.  THE CARTESIAN COORDINATES OF THESE ATOMS
*            ARE IN ARRAY XYZ.
*
*     DIHED IS A MODIFIED VERSION OF A SUBROUTINE OF THE SAME NAME
*           WHICH WAS WRITTEN BY DR. W. THEIL IN 1973.
*
*********************************************************************
      XI1=XYZ(1,I)-XYZ(1,K)
      XJ1=XYZ(1,J)-XYZ(1,K)
      XL1=XYZ(1,L)-XYZ(1,K)
      YI1=XYZ(2,I)-XYZ(2,K)
      YJ1=XYZ(2,J)-XYZ(2,K)
      YL1=XYZ(2,L)-XYZ(2,K)
      ZI1=XYZ(3,I)-XYZ(3,K)
      ZJ1=XYZ(3,J)-XYZ(3,K)
      ZL1=XYZ(3,L)-XYZ(3,K)
C      ROTATE AROUND Z AXIS TO PUT KJ ALONG Y AXIS
      DIST= SQRT(XJ1**2+YJ1**2+ZJ1**2)
      COSA=ZJ1/DIST
      IF(COSA.GT.1.0D0) COSA=1.0D0
      IF(COSA.LT.-1.0D0) COSA=-1.0D0
      DDD=1.0D0-COSA**2
      IF(DDD.LE.0.0) GO TO 10
      YXDIST=DIST* SQRT(DDD)
      IF(YXDIST.GT.1.0D-6) GO TO 20
   10 CONTINUE
      XI2=XI1
      XL2=XL1
      YI2=YI1
      YL2=YL1
      COSTH=COSA
      SINTH=0.D0
      GO TO 30
   20 COSPH=YJ1/YXDIST
      SINPH=XJ1/YXDIST
      XI2=XI1*COSPH-YI1*SINPH
      XL2=XL1*COSPH-YL1*SINPH
      YI2=XI1*SINPH+YI1*COSPH
      YJ2=XJ1*SINPH+YJ1*COSPH
      YL2=XL1*SINPH+YL1*COSPH
C      ROTATE KJ AROUND THE X AXIS SO KJ LIES ALONG THE Z AXIS
      COSTH=COSA
      SINTH=YJ2/DIST
   30 CONTINUE
      YI3=YI2*COSTH-ZI1*SINTH
      YL3=YL2*COSTH-ZL1*SINTH
      CALL DANG(XL2,YL3,XI2,YI3,ANGLE)
      IF (ANGLE .LT. 0.) ANGLE=4.0D0* ASIN(1.0D00)+ANGLE
      IF (ANGLE .GE. 6.2831853D0 ) ANGLE=0.D0
      RETURN
      END

C###############################################################################

      SUBROUTINE DANG(A1,A2,B1,B2,RCOS)
      IMPLICIT REAL*8 (A-H,O-Z)
**********************************************************************
*
*    DANG  DETERMINES THE ANGLE BETWEEN THE POINTS (A1,A2), (0,0),
*          AND (B1,B2).  THE RESULT IS PUT IN RCOS.
*
**********************************************************************
      ZERO=1.0D-6
      IF( ABS(A1).LT.ZERO.AND. ABS(A2).LT.ZERO) GO TO 10
      IF( ABS(B1).LT.ZERO.AND. ABS(B2).LT.ZERO) GO TO 10
      ANORM=1.0D0/ SQRT(A1**2+A2**2)
      BNORM=1.0D0/ SQRT(B1**2+B2**2)
      A1=A1*ANORM
      A2=A2*ANORM
      B1=B1*BNORM
      B2=B2*BNORM
      SINTH=(A1*B2)-(A2*B1)
      COSTH=A1*B1+A2*B2
      IF(COSTH.GT.1.0D0) COSTH=1.0D0
      IF(COSTH.LT.-1.0D0) COSTH=-1.0D0
      RCOS= ACOS(COSTH)
      IF( ABS(RCOS).LT.4.0D-4) GO TO 10
      IF(SINTH.GT.0.D0) RCOS=4.0D0* ASIN(1.0D00)-RCOS
      RCOS=-RCOS
      RETURN
   10 RCOS=0.0D0
      RETURN
      END

C###############################################################################
