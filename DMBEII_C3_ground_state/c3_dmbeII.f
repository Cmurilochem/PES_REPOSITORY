C################################################################################
C##############  GROUND-STATE SINGLET GLOBAL PES OF C3 - DMBEII #################
C################################################################################
C####### C.M.R. Rocha and A.J.C. Varandas, Phys. Chem. Chem. Phys (2018) ######## 
C################################################################################
C                                                                              ##
C USE CALL POTC3(R,POT)                                                        ##
C                                                                              ##
C INPUT COORDINATES IN BOHR                                                    ##
C                                                                              ## 
C OUTPUT ENERGIES IN HARTREE                                                   ##
C AND GIVEN WITH RESPECT TO INFINITELY SEPARATED C(3P)+C2(3PIu) FRAGMENTS      ##
C                                                                              ## 
C COORDINATES FOR THE DINFH GLOBAL MINIMUM (IN BOHR)                           ##
C                                                                              ##
C R(1)=2.444D+00                                                               ##
C                                                                              ##  
C R(2)=4.889D+00                                                               ##
C                                                                              ##
C R(3)=2.444D+00                                                               ##
C                                                                              ##
C ENERGY (IN HARTRE)=-0.290427D+00                                             ##
C                                                                              ##
C HARMONIC VIBRATIONAL FREQUENCIES (IN CM-1)                                   ##
C                                                                              ## 
C W1=1203.9                                                                    ##
C                                                                              ## 
C W2=61.0                                                                      ##
C                                                                              ##
C W3=2125.5                                                                    ##
C                                                                              ##
C COORDINATES FOR THE C2V ISOMERIZATION TRANSITION STATE (IN BOHR)             ##                                                               
C                                                                              ##
C R(1)=2.771D+00                                                               ##
C                                                                              ##  
C R(2)=2.401D+00                                                               ##
C                                                                              ##
C R(3)=2.771D+00                                                               ##
C                                                                              ##
C ENERGY (IN HARTRE)=-0.256353D+00                                             ##
C                                                                              ##
C HARMONIC VIBRATIONAL FREQUENCIES (IN CM-1)                                   ##
C                                                                              ## 
C W1=1295.2                                                                    ##
C                                                                              ## 
C W2=1840.5                                                                    ##
C                                                                              ##
C W3=1047.3i                                                                   ##
C                                                                              ## 
C COORDINATES FOR THE C2V VdW (LONG-RANGE) TRANSITION STATE (IN BOHR)          ##                                                               
C                                                                              ##
C R(1)=7.249D+00                                                               ##
C                                                                              ##  
C R(2)=7.249D+00                                                               ##
C                                                                              ##
C R(3)=2.470D+00                                                               ##
C                                                                              ##
C ENERGY (IN HARTRE)=-0.002558D+00                                             ##
C                                                                              ##
C HARMONIC VIBRATIONAL FREQUENCIES (IN CM-1)                                   ##
C                                                                              ## 
C W1=1618.1                                                                    ##
C                                                                              ## 
C W2=160.1                                                                     ##
C                                                                              ##
C W3=129.7i                                                                    ##
C                                                                              ##
C COORDINATES FOR THE C2V (LONG-RANGE) 2ND ORDER SADDLE POINT (IN BOHR)        ##                                                               
C                                                                              ##
C R(1)=5.511D+00                                                               ##
C                                                                              ##  
C R(2)=5.511D+00                                                               ##
C                                                                              ##
C R(3)=2.478D+00                                                               ##
C                                                                              ##
C ENERGY (IN HARTRE)=0.012504D+00                                              ##
C                                                                              ##
C HARMONIC VIBRATIONAL FREQUENCIES (IN CM-1)                                   ##
C                                                                              ## 
C W1=1345.3                                                                    ##
C                                                                              ## 
C W2=459.5i                                                                    ##
C                                                                              ##
C W3=546.5i                                                                    ##
C                                                                              ##
C INTEL AND GNU FORTRAN CONPILERS TESTED ON X86-64 FEDORA PLATFORM             ##
C                                                                              ##
C################################################################################

      SUBROUTINE POTC3(R,POT)
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(3) :: R
      DOUBLE PRECISION :: EHFD3H, EHFC2V, EHFTOT 
      DOUBLE PRECISION :: SCALEDPOT23, POT
      DOUBLE PRECISION, PARAMETER :: DD=0.2323414039409390507984199D+00
      DOUBLE PRECISION :: ONEBODY

      CALL SCALEDDYNPOT23(R(1),R(2),R(3),SCALEDPOT23)

      CALL THRBODYTOT(R(1),R(2),R(3),EHFD3H,EHFC2V,EHFTOT)

      ONEBODY=DD

      POT=SCALEDPOT23+EHFTOT+ONEBODY

      END SUBROUTINE POTC3

c  !!!*********************************************************************************!!!

      SUBROUTINE THRBODYTOT(D1,D2,D3,EHFD3H,EHFC2V,EHFTOT)
      IMPLICIT NONE
      DOUBLE PRECISION :: D1, D2, D3
      DOUBLE PRECISION, DIMENSION(144) :: C 
      DOUBLE PRECISION :: EHFC2V1, EHFC2V2
      DOUBLE PRECISION :: EHFD3H, EHFC2V, EHFTOT

       C(  1)=  0.939419642867667104D-01
       C(  2)=  0.814094356965121627D-02
       C(  3)= -0.137991765583622972D+00
       C(  4)= -0.120752320819536626D-01
       C(  5)= -0.526762995195253561D-01
       C(  6)=  0.344140515403429492D+00
       C(  7)= -0.144177778289186331D-01
       C(  8)= -0.216556741282060983D+00
       C(  9)= -0.118516399062917804D+00
       C( 10)=  0.194334439638676221D+00
       C( 11)= -0.672855804010018006D-02
       C( 12)=  0.445640961467239749D+00
       C( 13)=  0.210270114680647513D+01
       C( 14)= -0.473015437845383593D-01
       C( 15)=  0.424675612328713759D-01
       C( 16)= -0.142313896757015378D-02
       C( 17)=  0.922426072159836935D-02
       C( 18)=  0.779262046319446761D-01
       C( 19)= -0.430756618506507479D-03
       C( 20)=  0.101106392051225646D+01
       C( 21)=  0.168776997019681810D-01
       C( 22)=  0.310076807579411050D-01
       C( 23)= -0.555669706556694617D-03
       C( 24)=  0.189073348055313062D+00
       C( 25)=  0.531839785809717306D-02
       C( 26)=  0.659753281861308993D+00
       C( 27)= -0.405365137471814449D-01
       C( 28)=  0.187354019661505999D+00
       C( 29)=  0.813305566580927765D-02
       C( 30)=  0.782028657391569890D-02
       C( 31)= -0.102320710237646184D-03
       C( 32)=  0.139389147527451004D-02
       C( 33)=  0.551339890164715286D-02
       C( 34)=  0.871095286063372325D-02
       C( 35)=  0.142558452764864302D-02
       C( 36)=  0.982099724137223407D-01
       C( 37)=  0.766319032385726507D-02
       C( 38)=  0.253278266783274787D-01
       C( 39)=  0.815312358197797056D-03
       C( 40)=  0.920746721066789567D-03
       C( 41)=  0.437428583519401120D-04
       C( 42)= -0.381045127804672845D-04
       C( 43)=  0.175798202807236079D-02
       C( 44)= -0.466517441974516019D-03
       C( 45)=  0.351449350549256755D-02
       C( 46)= -0.235648639269382502D-02
       C( 47)=  0.672025193130964838D-04
       C( 48)= -0.212710137098015605D-02
       C( 49)=  0.215021725683087460D-03
       C( 50)=  0.134146357895935923D-03
       C( 51)=  0.641835269147592122D-04
       C( 52)= -0.113967961247292232D-04
       C( 53)=  0.102605859070874463D-04
       C( 54)=  0.552705799921816565D-01
       C( 55)=  0.867056281944996593D-01
       C( 56)= -0.292061509407157527D+00
       C( 57)=  0.340031286778664996D-01
       C( 58)=  0.127723520264307916D+00
       C( 59)=  0.126702812674413323D+01
       C( 60)= -0.153752772318713638D-01
       C( 61)= -0.253399962637563180D-02
       C( 62)= -0.121847965474670805D+00
       C( 63)=  0.762417017494457006D+00
       C( 64)= -0.102626123630331217D-01
       C( 65)=  0.439870788773157273D+00
       C( 66)=  0.168735339849976040D+01
       C( 67)= -0.593511612131693914D-01
       C( 68)=  0.207750399009005826D+00
       C( 69)= -0.884966701418230949D-03
       C( 70)=  0.722406868036933741D-02
       C( 71)=  0.380312893499337584D-01
       C( 72)=  0.340641951855910896D-01
       C( 73)=  0.535047684436933180D+00
       C( 74)=  0.262522419438343813D-01
       C( 75)=  0.643720125996914139D-01
       C( 76)=  0.101519265922605704D-02
       C( 77)=  0.341658750722979229D-01
       C( 78)=  0.149953120015734655D-02
       C( 79)=  0.108240178322483449D+00
       C( 80)= -0.157367215202812889D-01
       C( 81)=  0.438820536029984748D-01
       C( 82)=  0.583518008253752166D-02
       C( 83)=  0.679892780345060874D-02
       C( 84)=  0.138470992549592043D-03
       C( 85)=  0.326796638933332995D+01
       C( 86)=  0.609999999999999987D+00
       C( 87)= -0.818307840867775037D-03
       C( 88)= -0.763963412105609719D-03
       C( 89)= -0.585126332145601002D-03
       C( 90)= -0.254858222399380373D-02
       C( 91)= -0.507353537340084482D-01
       C( 92)=  0.246576863041835619D-01
       C( 93)= -0.891309408429813476D-03
       C( 94)=  0.850094949972518971D-01
       C( 95)=  0.400507858455017435D-01
       C( 96)=  0.272964768595807573D-01
       C( 97)= -0.901669218460573541D-03
       C( 98)= -0.117434181088866413D+00
       C( 99)= -0.157822962081838825D+00
       C(100)= -0.144728632687294745D+01
       C(101)= -0.354159905365007022D+00
       C(102)= -0.278247474532728756D+01
       C(103)= -0.221459200604135731D+01
       C(104)= -0.107676338335859406D+00
       C(105)= -0.136264888251007368D-01
       C(106)= -0.160230674075789725D-01
       C(107)=  0.789401227026172853D-01
       C(108)= -0.469425271984457021D-01
       C(109)=  0.892272218467323429D-02
       C(110)=  0.354883893697421487D-01
       C(111)= -0.139408714066350751D-01
       C(112)=  0.260960999999999999D+01
       C(113)=  0.261999637999999990D+01
       C(114)=  0.260440127999999982D+01
       C(115)=  0.329751996671108705D+01
       C(116)=  0.323232707276617127D-05
       C(117)=  0.486143040388596260D-04
       C(118)= -0.212510693796359815D-02
       C(119)=  0.213367046321904380D-03
       C(120)= -0.128324510130825532D-02
       C(121)= -0.573522106229466290D-02
       C(122)=  0.357905162584196482D-03
       C(123)=  0.197085359974491720D-02
       C(124)= -0.178597719776261250D-02
       C(125)=  0.215027561383702992D-03
       C(126)=  0.219562598118572434D-03
       C(127)= -0.135147884296237010D-01
       C(128)= -0.541527773071893911D-01
       C(129)=  0.197689938834732071D+00
       C(130)= -0.851429167298263928D-01
       C(131)=  0.323164050737253100D-01
       C(132)=  0.288432405228098654D+00
       C(133)= -0.478502565254514192D-01
       C(134)=  0.512262982168152151D-02
       C(135)=  0.117985605912901449D-01
       C(136)= -0.244890509588525018D-02
       C(137)=  0.745602598096316554D-02
       C(138)=  0.641217471119934646D-04
       C(139)= -0.782019372756585958D-02
       C(140)= -0.534177225704332440D-03
       C(141)=  0.331976404999999986D+01
       C(142)=  0.331120610000000015D+01
       C(143)=  0.332404301999999996D+01
       C(144)=  0.430430404981779002D+01

      CALL THRBODYD3H(D1,D2,D3,C,EHFD3H)

      CALL THRBODYC2V1(D1,D2,D3,C,EHFC2V1)

      CALL THRBODYC2V2(D1,D2,D3,C,EHFC2V2)

      EHFC2V=EHFC2V1+EHFC2V2
 
      EHFTOT=EHFD3H+EHFC2V

      END SUBROUTINE THRBODYTOT

c  !!!*********************************************************************************!!!

      SUBROUTINE THRBODYD3H(D1,D2,D3,C,EHFD3H)
      IMPLICIT NONE
      DOUBLE PRECISION :: D1, D2, D3
      INTEGER :: I, J, K, L, NUM
      DOUBLE PRECISION, DIMENSION(144) :: C
      DOUBLE PRECISION, DIMENSION(3) :: R  
      DOUBLE PRECISION, DIMENSION(3) :: Q 
      DOUBLE PRECISION, DIMENSION(3,3) :: A
      DOUBLE PRECISION :: TAU1, TAU2, TAU3
      DOUBLE PRECISION :: G1P1, G2P1, G3P1
      DOUBLE PRECISION :: G1P2, G2P2, G3P2
      INTEGER, PARAMETER :: ORDER1=9 
      INTEGER, PARAMETER :: ORDER2=7 
      DOUBLE PRECISION :: RDISP 
      DOUBLE PRECISION :: GAM 
      DOUBLE PRECISION :: POL1, POL2, EHFD3H 
      DOUBLE PRECISION :: DISP, RANG
      DOUBLE PRECISION :: LIN, RANGJT

      RDISP=C(85)
      GAM=C(86)
 
      R(1)=D1
      R(2)=D2
      R(3)=D3
      
      A(1,1)=SQRT(1.000D+00/3.00D+00)
      A(1,2)=A(1,1)
      A(1,3)=A(1,1)
      A(2,1)=0.000D+00
      A(2,2)=SQRT(1.000D+00/2.000D+00)
      A(2,3)=-A(2,2)
      A(3,1)=SQRT(2.000D+00/3.000D+00)
      A(3,2)=-SQRT(1.000D+00/6.000D+00)
      A(3,3)=A(3,2)
      
      NUM=0  
      POL1=0.000D+00
      POL2=0.000D+00
      EHFD3H=0.000D+00
        
      DO I=1,3
        Q(I)=0.00D+00
          DO J=1,3
            Q(I)=Q(I)+A(I,J)*DISP(R(J),RDISP)
          END DO
      END DO     
          
      TAU1=Q(1)
      TAU2=Q(2)**2+Q(3)**2
      TAU3=Q(3)**3-3.00D+00*Q(3)*Q(2)**2 
               
      DO L=0,ORDER1
        DO I=0,L
         IF (I==0) THEN
           G1P1=1.0000D+00
         ELSE 
           G1P1=TAU1**I
         END IF
         DO J=0,(L-I),2
           IF (J==0) THEN
             G2P1=1.0000D+00
           ELSE
             G2P1=TAU2**(J/2)
           END IF
           K=L-I-J
           IF (MOD(K,3)==0) THEN
             IF (K==0) THEN
               G3P1=1.0000D+00
             ELSE
               G3P1=TAU3**(K/3)
             END IF
             NUM=NUM+1
             POL1=POL1+C(NUM)*G1P1*G2P1*G3P1 
           END IF
         END DO
        END DO
      END DO

      NUM=53
  
      DO L=0,ORDER2
        DO I=0,L
          IF (I==0) THEN
            G1P2=1.0000D+00
          ELSE 
            G1P2=TAU1**I
          END IF
          DO J=0,(L-I),2
            IF (J==0) THEN
              G2P2=1.0000D+00
            ELSE
              G2P2=TAU2**(J/2)
            END IF
            K=L-I-J
            IF (MOD(K,3)==0) THEN
              IF (K==0) THEN
                G3P2=1.0000D+00
              ELSE
                G3P2=TAU3**(K/3)
              END IF
              NUM=NUM+1
              POL2=POL2+C(NUM)*G1P2*G2P2*G3P2 
            END IF
          END DO
        END DO
      END DO

      LIN=SQRT(TAU2)*RANGJT(D1,D2,D3)

      EHFD3H=(POL1-LIN*POL2)*RANG(D1,D2,D3,RDISP,GAM)
             
      END SUBROUTINE THRBODYD3H

c  !!!*********************************************************************************!!!

      SUBROUTINE THRBODYC2V1(D1,D2,D3,C,EHFC2V1)
      IMPLICIT NONE
      DOUBLE PRECISION :: D1, D2, D3
      INTEGER :: I, J, K, L, NUMC2V
      DOUBLE PRECISION, DIMENSION(144) :: C
      DOUBLE PRECISION, DIMENSION(3) :: RC2V  
      DOUBLE PRECISION, DIMENSION(3) :: QC2V 
      DOUBLE PRECISION, DIMENSION(3,3) :: A
      DOUBLE PRECISION :: TAU1C2V, TAU2C2V, TAU3C2V
      DOUBLE PRECISION :: G1P1C2V, G2P1C2V, G3P1C2V
      DOUBLE PRECISION :: G1P2C2V, G2P2C2V, G3P2C2V
      DOUBLE PRECISION :: G1P3C2V, G2P3C2V, G3P3C2V
      INTEGER, PARAMETER :: ORDER1C2V=4 
      INTEGER, PARAMETER :: ORDER2C2V=3
      INTEGER, PARAMETER :: ORDER3C2V=3
      DOUBLE PRECISION :: POL1C2V, POL2C2V, POL3C2V, EHFC2V1 
      DOUBLE PRECISION :: DISP, RDISPC2V, X1, X2, RANGC2V
      DOUBLE PRECISION :: GAM1
      DOUBLE PRECISION :: S1, S2, S3
      DOUBLE PRECISION :: DELTAS,  Q31, RHO0, Q1ABS, Q2ABS, Q3ABS
      DOUBLE PRECISION, DIMENSION(3) :: DELTA, THETA, NORM
      DOUBLE PRECISION, DIMENSION(3) :: Q3SEAM, Q2SEAM
      DOUBLE PRECISION :: PI, ZERO
      INTEGER :: PHASE
      DOUBLE PRECISION :: LIN, RANGJT

      RDISPC2V=C(112)
      X1=C(113)
      X2=C(114)
      GAM1=C(115)

      RC2V(1)=D1
      RC2V(2)=D2
      RC2V(3)=D3
      
      A(1,1)=SQRT(1.000D+00/3.00D+00)
      A(1,2)=A(1,1)
      A(1,3)=A(1,1)
      A(2,1)=0.000D+00
      A(2,2)=SQRT(1.000D+00/2.000D+00)
      A(2,3)=-A(2,2)
      A(3,1)=SQRT(2.000D+00/3.000D+00)
      A(3,2)=-SQRT(1.000D+00/6.000D+00)
      A(3,3)=A(3,2)
      
      NUMC2V=86
      POL1C2V=0.000D+00
      POL2C2V=0.000D+00
      POL3C2V=0.000D+00
      EHFC2V1=0.000D+00
        
      DO I=1,3
        QC2V(I)=0.00D+00
          DO J=1,3
            QC2V(I)=QC2V(I)+A(I,J)*DISP(RC2V(J),RDISPC2V)
          END DO
      END DO   
          
      TAU1C2V=QC2V(1)
      TAU2C2V=QC2V(2)**2+QC2V(3)**2
      TAU3C2V=QC2V(3)**3-3.00D+00*QC2V(3)*QC2V(2)**2 
               
      DO L=0,ORDER1C2V
        DO I=0,L
         IF (I==0) THEN
           G1P1C2V=1.0000D+00
         ELSE 
           G1P1C2V=TAU1C2V**I
         END IF
         DO J=0,(L-I),2
           IF (J==0) THEN
             G2P1C2V=1.0000D+00
           ELSE
             G2P1C2V=TAU2C2V**(J/2)
           END IF
           K=L-I-J
           IF (MOD(K,3)==0) THEN
             IF (K==0) THEN
               G3P1C2V=1.0000D+00
             ELSE
               G3P1C2V=TAU3C2V**(K/3)
             END IF
             NUMC2V=NUMC2V+1
             POL1C2V=POL1C2V+C(NUMC2V)*G1P1C2V*G2P1C2V*G3P1C2V 
           END IF
         END DO
        END DO
      END DO 

      NUMC2V=97
  
      DO L=0,ORDER2C2V
        DO I=0,L
          IF (I==0) THEN
            G1P2C2V=1.0000D+00
          ELSE 
            G1P2C2V=TAU1C2V**I
          END IF
          DO J=0,(L-I),2
            IF (J==0) THEN
              G2P2C2V=1.0000D+00
            ELSE
              G2P2C2V=TAU2C2V**(J/2)
            END IF
            K=L-I-J
            IF (MOD(K,3)==0) THEN
              IF (K==0) THEN
                G3P2C2V=1.0000D+00
              ELSE
                G3P2C2V=TAU3C2V**(K/3)
              END IF
              NUMC2V=NUMC2V+1
              POL2C2V=POL2C2V+C(NUMC2V)*G1P2C2V*G2P2C2V*G3P2C2V 
            END IF
          END DO
        END DO
      END DO

      NUMC2V=104

      DO L=0,ORDER3C2V
        DO I=0,L
          IF (I==0) THEN
            G1P3C2V=1.0000D+00
          ELSE 
            G1P3C2V=TAU1C2V**I
          END IF
          DO J=0,(L-I),2
            IF (J==0) THEN
              G2P3C2V=1.0000D+00
            ELSE
              G2P3C2V=TAU2C2V**(J/2)
            END IF
            K=L-I-J
            IF (MOD(K,3)==0) THEN
              IF (K==0) THEN
                G3P3C2V=1.0000D+00
              ELSE
                G3P3C2V=TAU3C2V**(K/3)
              END IF
              NUMC2V=NUMC2V+1
              POL3C2V=POL3C2V+C(NUMC2V)*G1P3C2V*G2P3C2V*G3P3C2V 
            END IF
          END DO
        END DO
      END DO

      Q1ABS=A(1,1)*(D1+D2+D3)
      Q2ABS=A(2,2)*D2+A(2,3)*D3
      Q3ABS=A(3,1)*D1+A(3,2)*D2+A(3,3)*D3

      PI=3.1415926535897932D+00
 
      ZERO=0.00000000000000D+00

      THETA=(/(PI/2.00D+00),(7.00D+00*PI/6.00D+00),
     &        (11.00D+00*PI/6.00D+00)/)

      RHO0=ABS(Q31(Q1ABS))

      IF (Q31(Q1ABS)>=ZERO) THEN
        PHASE=0
      ELSE IF (Q31(Q1ABS)<ZERO) THEN
        PHASE=1
      ENDIF

      Q3SEAM(1)=RHO0*SIN(THETA(1)+DBLE(PHASE)*PI)
      Q2SEAM(1)=RHO0*COS(THETA(1)+DBLE(PHASE)*PI)
      Q3SEAM(2)=RHO0*SIN(THETA(2)+DBLE(PHASE)*PI)
      Q2SEAM(2)=RHO0*COS(THETA(2)+DBLE(PHASE)*PI)
      Q3SEAM(3)=RHO0*SIN(THETA(3)+DBLE(PHASE)*PI)
      Q2SEAM(3)=RHO0*COS(THETA(3)+DBLE(PHASE)*PI)

      DELTA(1)=SQRT((Q2ABS-Q2SEAM(1))**2+
     &              (Q3ABS-Q3SEAM(1))**2)
      DELTA(2)=SQRT((Q2ABS-Q2SEAM(2))**2+
     &              (Q3ABS-Q3SEAM(2))**2)
      DELTA(3)=SQRT((Q2ABS-Q2SEAM(3))**2+
     &              (Q3ABS-Q3SEAM(3))**2)

      S1=A(1,1)*(DELTA(1)+DELTA(2)+DELTA(3))
      S2=A(2,2)*DELTA(2)+A(2,3)*DELTA(3)
      S3=A(3,1)*DELTA(1)+A(3,2)*DELTA(2)+A(3,3)*DELTA(3)

      DELTAS=SQRT(S2**2+S3**2)

      LIN=SQRT(TAU2C2V)*RANGJT(D1,D2,D3)

      EHFC2V1=(POL1C2V-DELTAS*POL2C2V-
     &         LIN*POL3C2V)*
     &         RANGC2V(D1,D2,D3,X1,X2,GAM1)

      END SUBROUTINE THRBODYC2V1

c  !!!*********************************************************************************!!!

      SUBROUTINE THRBODYC2V2(D1,D2,D3,C,EHFC2V2)
      IMPLICIT NONE
      DOUBLE PRECISION :: D1, D2, D3
      INTEGER :: I, J, K, L, NUMC2V
      DOUBLE PRECISION, DIMENSION(144) :: C
      DOUBLE PRECISION, DIMENSION(3) :: RC2V  
      DOUBLE PRECISION, DIMENSION(3) :: QC2V 
      DOUBLE PRECISION, DIMENSION(3,3) :: A
      DOUBLE PRECISION :: TAU1C2V, TAU2C2V, TAU3C2V
      DOUBLE PRECISION :: G1P1C2V, G2P1C2V, G3P1C2V
      DOUBLE PRECISION :: G1P2C2V, G2P2C2V, G3P2C2V
      DOUBLE PRECISION :: G1P3C2V, G2P3C2V, G3P3C2V
      INTEGER, PARAMETER :: ORDER1C2V=4 
      INTEGER, PARAMETER :: ORDER2C2V=3
      INTEGER, PARAMETER :: ORDER3C2V=3
      DOUBLE PRECISION :: POL1C2V, POL2C2V, POL3C2V, EHFC2V2
      DOUBLE PRECISION :: DISP, RDISPC2V, X1, X2, RANGC2V
      DOUBLE PRECISION :: GAM1
      DOUBLE PRECISION :: S1, S2, S3
      DOUBLE PRECISION :: DELTAS, Q31, RHO0, Q1ABS, Q2ABS, Q3ABS
      DOUBLE PRECISION, DIMENSION(3) :: DELTA, THETA, NORM
      DOUBLE PRECISION, DIMENSION(3) :: Q3SEAM, Q2SEAM
      DOUBLE PRECISION :: PI, ZERO
      INTEGER :: PHASE
      DOUBLE PRECISION :: LIN, RANGJT

      RDISPC2V=C(141)
      X1=C(142)
      X2=C(143)
      GAM1=C(144)

      RC2V(1)=D1
      RC2V(2)=D2
      RC2V(3)=D3
      
      A(1,1)=SQRT(1.000D+00/3.00D+00)
      A(1,2)=A(1,1)
      A(1,3)=A(1,1)
      A(2,1)=0.000D+00
      A(2,2)=SQRT(1.000D+00/2.000D+00)
      A(2,3)=-A(2,2)
      A(3,1)=SQRT(2.000D+00/3.000D+00)
      A(3,2)=-SQRT(1.000D+00/6.000D+00)
      A(3,3)=A(3,2)
      
      NUMC2V=115
      POL1C2V=0.000D+00
      POL2C2V=0.000D+00
      POL3C2V=0.000D+00
      EHFC2V2=0.000D+00
        
      DO I=1,3
        QC2V(I)=0.00D+00
          DO J=1,3
            QC2V(I)=QC2V(I)+A(I,J)*DISP(RC2V(J),RDISPC2V)
          END DO
      END DO   
          
      TAU1C2V=QC2V(1)
      TAU2C2V=QC2V(2)**2+QC2V(3)**2
      TAU3C2V=QC2V(3)**3-3.00D+00*QC2V(3)*QC2V(2)**2 
               
      DO L=0,ORDER1C2V
        DO I=0,L
         IF (I==0) THEN
           G1P1C2V=1.0000D+00
         ELSE 
           G1P1C2V=TAU1C2V**I
         END IF
         DO J=0,(L-I),2
           IF (J==0) THEN
             G2P1C2V=1.0000D+00
           ELSE
             G2P1C2V=TAU2C2V**(J/2)
           END IF
           K=L-I-J
           IF (MOD(K,3)==0) THEN
             IF (K==0) THEN
               G3P1C2V=1.0000D+00
             ELSE
               G3P1C2V=TAU3C2V**(K/3)
             END IF
             NUMC2V=NUMC2V+1
             POL1C2V=POL1C2V+C(NUMC2V)*G1P1C2V*G2P1C2V*G3P1C2V 
           END IF
         END DO
        END DO
      END DO 

      NUMC2V=126
  
      DO L=0,ORDER2C2V
        DO I=0,L
          IF (I==0) THEN
            G1P2C2V=1.0000D+00
          ELSE 
            G1P2C2V=TAU1C2V**I
          END IF
          DO J=0,(L-I),2
            IF (J==0) THEN
              G2P2C2V=1.0000D+00
            ELSE
              G2P2C2V=TAU2C2V**(J/2)
            END IF
            K=L-I-J
            IF (MOD(K,3)==0) THEN
              IF (K==0) THEN
                G3P2C2V=1.0000D+00
              ELSE
                G3P2C2V=TAU3C2V**(K/3)
              END IF
              NUMC2V=NUMC2V+1
              POL2C2V=POL2C2V+C(NUMC2V)*G1P2C2V*G2P2C2V*G3P2C2V 
            END IF
          END DO
        END DO
      END DO

      NUMC2V=133

      DO L=0,ORDER3C2V
        DO I=0,L
          IF (I==0) THEN
            G1P3C2V=1.0000D+00
          ELSE 
            G1P3C2V=TAU1C2V**I
          END IF
          DO J=0,(L-I),2
            IF (J==0) THEN
              G2P3C2V=1.0000D+00
            ELSE
              G2P3C2V=TAU2C2V**(J/2)
            END IF
            K=L-I-J
            IF (MOD(K,3)==0) THEN
              IF (K==0) THEN
                G3P3C2V=1.0000D+00
              ELSE
                G3P3C2V=TAU3C2V**(K/3)
              END IF
              NUMC2V=NUMC2V+1
              POL3C2V=POL3C2V+C(NUMC2V)*G1P3C2V*G2P3C2V*G3P3C2V 
            END IF
          END DO
        END DO
      END DO

      Q1ABS=A(1,1)*(D1+D2+D3)
      Q2ABS=A(2,2)*D2+A(2,3)*D3
      Q3ABS=A(3,1)*D1+A(3,2)*D2+A(3,3)*D3

      PI=3.1415926535897932D+00

      THETA=(/(PI/2.00D+00),(7.00D+00*PI/6.00D+00),
     &        (11.00D+00*PI/6.00D+00)/)

      RHO0=ABS(Q31(Q1ABS))

      IF (Q31(Q1ABS)>=ZERO) THEN
        PHASE=0
      ELSE IF (Q31(Q1ABS)<ZERO) THEN
        PHASE=1
      ENDIF

      Q3SEAM(1)=RHO0*SIN(THETA(1)+DBLE(PHASE)*PI)
      Q2SEAM(1)=RHO0*COS(THETA(1)+DBLE(PHASE)*PI)
      Q3SEAM(2)=RHO0*SIN(THETA(2)+DBLE(PHASE)*PI)
      Q2SEAM(2)=RHO0*COS(THETA(2)+DBLE(PHASE)*PI)
      Q3SEAM(3)=RHO0*SIN(THETA(3)+DBLE(PHASE)*PI)
      Q2SEAM(3)=RHO0*COS(THETA(3)+DBLE(PHASE)*PI)

      DELTA(1)=SQRT((Q2ABS-Q2SEAM(1))**2+
     &              (Q3ABS-Q3SEAM(1))**2)
      DELTA(2)=SQRT((Q2ABS-Q2SEAM(2))**2+
     &              (Q3ABS-Q3SEAM(2))**2)
      DELTA(3)=SQRT((Q2ABS-Q2SEAM(3))**2+
     &              (Q3ABS-Q3SEAM(3))**2)

      S1=A(1,1)*(DELTA(1)+DELTA(2)+DELTA(3))
      S2=A(2,2)*DELTA(2)+A(2,3)*DELTA(3)
      S3=A(3,1)*DELTA(1)+A(3,2)*DELTA(2)+A(3,3)*DELTA(3)

      DELTAS=SQRT(S2**2+S3**2)

      LIN=SQRT(TAU2C2V)*RANGJT(D1,D2,D3)

      EHFC2V2=(POL1C2V-DELTAS*POL2C2V-
     &         LIN*POL3C2V)*
     &         RANGC2V(D1,D2,D3,X1,X2,GAM1)

      END SUBROUTINE THRBODYC2V2

c  !!!*********************************************************************************!!!

      SUBROUTINE SCALEDDYNPOT23(R1,R2,R3,SCALEDPOT23)
      IMPLICIT NONE
      DOUBLE PRECISION :: R1, R2, R3, SCALEDPOT23 
      DOUBLE PRECISION :: TWOBODYDYN, TWOBODYHF, VDC3
      DOUBLE PRECISION :: DEE_DIATOMIC, SWIT
      DOUBLE PRECISION, DIMENSION(3) :: EHF, DYN, V2, TWOBODY 

      DEE_DIATOMIC=0.000000000D+00

      CALL POT2BODY(R1,EHF(1),DYN(1),V2(1))
      CALL POT2BODY(R2,EHF(2),DYN(2),V2(2))
      CALL POT2BODY(R3,EHF(3),DYN(3),V2(3))
 
      TWOBODY(1)=DYN(1)*(1-SWIT(R2,R3,R1))*(1-SWIT(R3,R1,R2))

      TWOBODY(2)=DYN(2)*(1-SWIT(R1,R2,R3))*(1-SWIT(R3,R1,R2))

      TWOBODY(3)=DYN(3)*(1-SWIT(R1,R2,R3))*(1-SWIT(R2,R3,R1))

      TWOBODYDYN=SUM(TWOBODY)

      TWOBODYHF=SUM(EHF)

      CALL POT3DC(R1,R2,R3,VDC3)
  
      SCALEDPOT23=TWOBODYHF+TWOBODYDYN+3*DEE_DIATOMIC+VDC3
  
      END SUBROUTINE SCALEDDYNPOT23

c  !!!*********************************************************************************!!!

      SUBROUTINE POT3DC(R1,R2,R3,VDC3)
      IMPLICIT NONE
      DOUBLE PRECISION :: R1, R2, R3, DC3, VDC3
      DOUBLE PRECISION, DIMENSION(3) :: RE, rjac, cosgamma, RB, RC
      DOUBLE PRECISION, DIMENSION(10,3) :: DM, RM
      DOUBLE PRECISION, DIMENSION(3,10,3) :: a, b
      DOUBLE PRECISION, DIMENSION(3,3) :: PLegend
      DOUBLE PRECISION, DIMENSION(10) :: C_DISP_DIATOM
      DOUBLE PRECISION, DIMENSION(10,3,3) :: C_RE
      DOUBLE PRECISION, DIMENSION(10,3) :: C_RE_GAMMA
      INTEGER :: lmax
      INTEGER :: i, n, l, channel, poli
      DOUBLE PRECISION :: r
      DOUBLE PRECISION, DIMENSION(10,3) :: POLI_A, POLI_B
      DOUBLE PRECISION, PARAMETER :: R0_ATOM_DIATOM = 10.983323302739D+00   !! LEE ROY RADIUS (R0) OF C-Mg
      DOUBLE PRECISION, PARAMETER :: R0 = 7.89107343775231D+00 
      DOUBLE PRECISION :: JABOBI_RADIUS_EXP, DAMP_PART, SWIT_PART 
      DOUBLE PRECISION :: DAMP, SWIT

      CALL TRIANG2JACOBI(R1,R2,R3,RE,rjac,cosgamma,RB,RC)

      RM(6,1) = 4.5000D+00;       DM(6,1) = 27.5279D+00;
      a(1,6,1) = 0.99499999D+00;  a(2,6,1) = 0.25778783D+00;
      a(3,6,1) = 0.00139584D+00;  b(1,6,1) = a(1,6,1);
      b(2,6,1) = 0.27435311D+00;  b(3,6,1) = 0.02546246D+00;

      RM(6,2) = 4.5000D+00;       DM(6,2) = 28.5233D+00;
      a(1,6,2) = 0.76362610D+00;  a(2,6,2) = 0.21435356D+00;
      a(3,6,2) = -0.00008571D+00; b(1,6,2) = a(1,6,2);
      b(2,6,2) = 0.24930236D+00;  b(3,6,2) = 0.00676850D+00;
 
      RM(8,1) = 4.4711D+00;       DM(8,1) = 1037.6370D+00
      a(1,8,1) = 0.93135973D+00;  a(2,8,1) = 0.23731725D+00;
      a(3,8,1) = -0.00013049D+00; b(1,8,1) = a(1,8,1);
      b(2,8,1) = 0.24399061D+00;  b(3,8,1)= 0.02486858;
 
      RM(8,2) = 4.4823D+00;       DM(8,2) = 2941.9456D+00;
      a(1,8,2) = 0.79819039D+00;  a(2,8,2) = 0.22589941D+00;
      a(3,8,2) = -0.00082011D+00; b(1,8,2) = a(1,8,2);
      b(2,8,2) = 0.29644898D+00;  b(3,8,2) = 0.02109793D+00;
 
      RM(8,3) = 4.4873D+00;       DM(8,3) = 452.6876D+00;
      a(1,8,3) = 1.23539671D+00;  a(2,8,3) = 0.57644149D+00;
      a(3,8,3) = 0.07460865D+00;  b(1,8,3) = a(1,8,3);
      b(2,8,3) = 0.58644002D+00;  b(3,8,3) = 0.07326621D+00;
  
      RM(10,1) = 4.4525D+00;      DM(10,1) = 46263.4796D+00;
      a(1,10,1) = 0.89618614D+00; a(2,10,1) = 0.22325495D+00;
      a(3,10,1)= -0.00182731D+00; b(1,10,1) = a(1,10,1);
      b(2,10,1)= 0.23140031D+00;  b(3,10,1) = 0.02681374D+00;

      DO channel=1,3  !! LOOP FOR EACH CHANNEL

      PLegend(1,channel) = 1.0D+00      

      PLegend(2,channel) = 0.5D+00*(3.0D+00*cosgamma(channel)**2
c23456
     + - 1.0D+00)                                     

      PLegend(3,channel) = (35.0D+00*cosgamma(channel)**4 
c23456
     + - 30.0D+00*cosgamma(channel)**2 + 3.0D+00)/8.0D+00   

      END DO

      C_DISP_DIATOM(6) = 40.9D+00

      C_DISP_DIATOM(8) = C_DISP_DIATOM(6)*R0**(1.54)  

      C_DISP_DIATOM(10) = C_DISP_DIATOM(6)*1.31D+00*R0**(2.0*1.54)

      DO channel=1,3         
      DO n=6,10,2            

      IF (n == 6) THEN  

      lmax = 2           
 
      ELSE IF (n == 8) THEN

      lmax = 3             

      ELSE IF (n == 10) THEN 

      lmax = 1             

      END IF 

      DO l=1,lmax

      r = RE(channel) - RM(n,l)

      POLI_A(n,l) = 0.0D+00

      POLI_B(n,l) = 0.0D+00

      DO poli = 1,3
     
      POLI_A(n,l) = POLI_A(n,l) + a(poli,n,l)*r**(poli)

      POLI_B(n,l) = POLI_B(n,l) - b(poli,n,l)*r**(poli)
           
      END DO

      IF (l == 1) THEN

      C_RE(n,l,channel) = C_DISP_DIATOM(n) + C_DISP_DIATOM(n)
c23456
     + + DM(n,l)*(1.0D+00 + POLI_A(n,l))*EXP(POLI_B(n,l))		 

      ELSE

      C_RE(n,l,channel) =
c23456
     +  DM(n,l)*(1.0D+00 + POLI_A(n,l))*EXP(POLI_B(n,l)) 
                     
      END IF

      END DO
      END DO                  
      END DO                  

      DO channel=1,3         
      DO n=6,10,2            

      IF (n == 6) THEN  

      lmax = 2            
 
      ELSE IF (n == 8) THEN

      lmax = 3             

      ELSE IF (n == 10) THEN 

      lmax = 1             

      END IF 

      C_RE_GAMMA(n,channel) = 0.0D+00

      DO l=1,lmax

      C_RE_GAMMA(n,channel) = C_RE_GAMMA(n,channel) 
c23456
     + + C_RE(n,l,channel)*PLegend(l,channel)   
     
      END DO

      END DO                
      END DO                 

      DC3 = 0.0D+00

      JABOBI_RADIUS_EXP = 0.0D+00    

      DAMP_PART = 0.0D+00            

      SWIT_PART = 0.0D+00            

      DO channel = 1,3

      DO n=6,10,2
  
      JABOBI_RADIUS_EXP = rjac(channel)**(-n)

      DAMP_PART = DAMP(rjac(channel),R0_ATOM_DIATOM,n)

      SWIT_PART = SWIT(RE(channel),RB(channel),RC(channel))

      DC3 = DC3 + 
c23456
     + SWIT_PART*C_RE_GAMMA(n,channel)*DAMP_PART*JABOBI_RADIUS_EXP   

      END DO

      END DO

      VDC3 = - DC3

      END SUBROUTINE POT3DC


  !!!*********************************************************************************!!!

      SUBROUTINE TRIANG2JACOBI(R1,R2,R3,RE,rjac,cosgamma,RB,RC)
      IMPLICIT NONE
      DOUBLE PRECISION :: R1, R2, R3
      DOUBLE PRECISION, PARAMETER  :: ZERO = 1.0D-12
      DOUBLE PRECISION, DIMENSION(3) :: RE
      DOUBLE PRECISION, DIMENSION(3) :: rjac
      DOUBLE PRECISION, DIMENSION(3) :: RB, RC
      DOUBLE PRECISION, DIMENSION(3) :: squaredrjac
      DOUBLE PRECISION, DIMENSION(3) :: cosgamma
      DOUBLE PRECISION, PARAMETER  :: pi = 3.1415926535897932384626433832795D+00
      INTEGER :: i

      DO i = 1,3      !!! CYCLIC PERMUTATIONS OF 1,2,3

      IF (i == 1) THEN        !! CHANNEL 1 
 
      RE(i) = R1
      RB(i) = R2
      RC(i) = R3

      ELSE IF (i == 2) THEN    !! CHANNEL 2

      RE(i) = R2
      RB(i) = R3
      RC(i) = R1

      ELSE IF (i == 3) THEN    !! CHANNEL 3

      RE(i) = R3
      RB(i) = R1
      RC(i) = R2

      END IF

      rjac(i) = ZERO

      squaredrjac(i) = 0.5D+00*(RB(i)**2 + RC(i)**2 - 0.5D+00*RE(i)**2)

      IF (squaredrjac(i) > 0) THEN 

      rjac(i) = sqrt(squaredrjac(i)) 

      END IF
    
      cosgamma(i)=(0.5D+00*(RC(i)**2-RB(i)**2))/(rjac(i)*RE(i))

      IF (cosgamma(i) > 1.0D+00) THEN 

      cosgamma(i) = 1.0D+00 

      ELSE IF (cosgamma(i) < -1.0D+00) THEN 
   
      cosgamma(i) = -1.0D+00

      END IF 

      END DO

      END SUBROUTINE TRIANG2JACOBI

c  !!!*********************************************************************************!!!

      SUBROUTINE POT2BODY(R,VEHF,VDC,VCC)
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(4) :: A, B, QUI, CN, DC
      INTEGER, DIMENSION(4) :: ID = (/5,6,8,10/)
      DOUBLE PRECISION, DIMENSION(10) :: COEF, POLI 
      INTEGER :: I
      DOUBLE PRECISION :: R, RE, X, R0, RHO, RHO_LINHA  
      DOUBLE PRECISION :: VDC, D, GAMMA, VEHF, G0, G1, G2, VCC, DEE
      DOUBLE PRECISION :: AN, BN
  
      RE = 2.47932D+00

      X = R - RE

      D = 0.467489D+00
 
      G0 = 0.286063D+00

      G1 = 15.5532D+00

      G2 = 0.0736335D+00	

      DEE = 0.000000000D+00!(0.2323414039409390507984199D+00)

      GAMMA =  G0*(1.0D+00 + G1*TANH(G2*X))

      COEF(1) = 0.850571D+00
      COEF(2) = -1.07115D+00
      COEF(3) = 0.969122D+00
      COEF(4) = -0.612554D+00
      COEF(5) = 0.310833D+00
      COEF(6) = -0.0944244D+00
      COEF(7) = -0.00736976D+00
      COEF(8) = 0.0137852D+00
      COEF(9) = -0.00333289D+00
      COEF(10) = 0.000257898D+00

      DO I=1,10
      POLI(I)=COEF(I)*X**I	
      END DO

      VEHF = -(D/R)*(1.0D+00 + SUM(POLI))*EXP(-GAMMA*X) 

      R0 = 7.89107343775231D+00   

      RHO = 5.5D+00 + 1.25D+00*R0 

      RHO_LINHA = R/RHO

      CN(1) = 14.5404D+00
      CN(2) = 40.9D+00
      CN(3) = CN(2)*R0**(1.54)
      CN(4) = CN(2)*1.31D+00*R0**(2.0*1.54)

      DO I=1,4
      A(I) = AN(ID(I))
      END DO

      DO I=1,4
      B(I) = BN(ID(I))
      END DO

      DO I=1,4
      QUI(I)=(1.0D+00-EXP(-A(I)*RHO_LINHA-B(I)*RHO_LINHA**2))**(ID(I))
      END DO 

      DO I=1,4
      DC(I)= - CN(I)*QUI(I)/R**(ID(I)) 
      END DO 
                            
      VDC = SUM(DC)                           
                                          
      VCC = VEHF + VDC + DEE

      END SUBROUTINE POT2BODY

c  !!!*********************************************************************************!!!

      DOUBLE PRECISION FUNCTION Q31(Q1)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN):: Q1
      DOUBLE PRECISION :: A, B, C, D, E, F, G, H 
      A=0.00454872897741226D+00
      B=0.02509103210095290D+00
      C=1.05938652067574000D+00
      D=4.78291330706927000D+00
      E=-0.9261734635796440D+00
      F=-0.2164776500990330D+00
      G=0.17221936619583700D+00
      H=0.68281911615943100D+00
      Q31=A-B*TANH(C*(Q1-D)+E*(Q1-D)**2+
     &     F*(Q1-D)**3+G*(Q1-D)**4+H*(Q1-D)**5)
      END FUNCTION Q31

c  !!!*********************************************************************************!!!

      DOUBLE PRECISION FUNCTION RANGJT(R1,R2,R3)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN):: R1,R2,R3
      DOUBLE PRECISION, PARAMETER :: X0=2.885169410D+00
      DOUBLE PRECISION, PARAMETER :: GAM1=10000.00D+00
      RANGJT=(1.00D+00-EXP(-GAM1*((R1-X0)**2+(R2-X0)**2+(R3-X0)**2)))
      END FUNCTION RANGJT

c  !!!*********************************************************************************!!!

      DOUBLE PRECISION FUNCTION DISP(R,RDISP)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN):: R,RDISP
      DISP=(R-RDISP)
      END FUNCTION DISP 

c  !!!*********************************************************************************!!!

      DOUBLE PRECISION FUNCTION RANG(R1,R2,R3,X0,GAM1)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN):: R1,R2,R3,X0,GAM1
      RANG=(1.00D+00-TANH(GAM1*(R1-X0)))*
     &     (1.00D+00-TANH(GAM1*(R2-X0)))*
     &     (1.00D+00-TANH(GAM1*(R3-X0)))
      END FUNCTION RANG

c  !!!*********************************************************************************!!!

      DOUBLE PRECISION FUNCTION RANGC2V(R1,R2,R3,X1,X2,GAM1)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN):: R1,R2,R3,X1,X2,GAM1
      DOUBLE PRECISION, DIMENSION(3) :: RANG
      RANG(1)=(1.00D+00-TANH(GAM1*(R1-X1)))*
     &        (1.00D+00-TANH(GAM1*(R2-X2)))*
     &        (1.00D+00-TANH(GAM1*(R3-X2)))
      RANG(2)=(1.00D+00-TANH(GAM1*(R3-X1)))*
     &        (1.00D+00-TANH(GAM1*(R1-X2)))*
     &        (1.00D+00-TANH(GAM1*(R2-X2)))
      RANG(3)=(1.00D+00-TANH(GAM1*(R2-X1)))*
     &        (1.00D+00-TANH(GAM1*(R3-X2)))*
     &        (1.00D+00-TANH(GAM1*(R1-X2)))
      RANGC2V=RANG(1)+RANG(2)+RANG(3)
      END FUNCTION RANGC2V

c  !!!*********************************************************************************!!!

      DOUBLE PRECISION FUNCTION DAMP(rjac,R0,N)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN):: rjac, R0
      INTEGER, INTENT(IN):: N
      DOUBLE PRECISION :: RHO, RHO_LINHA, AN, BN
      DOUBLE PRECISION, DIMENSION(10) :: A, B
      RHO = 5.5D+00 + 1.25D+00*R0 
      RHO_LINHA = rjac/RHO
      A(N) = AN(N)
      B(N) = BN(N)
      DAMP = (1.0D+00-EXP(-A(N)*RHO_LINHA-B(N)*RHO_LINHA**2))**N
      END FUNCTION DAMP

c  !!!*********************************************************************************!!!

      DOUBLE PRECISION FUNCTION AN(N)
      IMPLICIT NONE
      DOUBLE PRECISION, PARAMETER :: ALFA0 = 16.36606000D+00
      DOUBLE PRECISION, PARAMETER :: ALFA1 = 0.70172000D+00
      INTEGER, INTENT(IN):: N
      AN=ALFA0/(DBLE(N))**ALFA1
      END FUNCTION AN 

c  !!!*********************************************************************************!!!

      DOUBLE PRECISION FUNCTION BN(N)
      IMPLICIT NONE
      DOUBLE PRECISION, PARAMETER :: BETA0 = 17.19338D+00
      DOUBLE PRECISION, PARAMETER :: BETA1 = 0.09574D+00
      INTEGER, INTENT(IN):: N
      BN=BETA0*DEXP(-BETA1*DBLE(N))
      END FUNCTION BN

c  !!!*********************************************************************************!!!

      DOUBLE PRECISION FUNCTION SWIT(RA,RB,RC)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: RA, RB, RC
      DOUBLE PRECISION, PARAMETER :: eta = 6.0D+00, csi = 1.0D+00
      SWIT = 0.5D+00*(1.0D+00-TANH(csi*(eta*RA-RB-RC)))
      END FUNCTION SWIT

c  !!!*********************************************************************************!!!
