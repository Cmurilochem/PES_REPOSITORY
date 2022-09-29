      PROGRAM TEST
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(6) :: R
      DOUBLE PRECISION, DIMENSION(3) :: W
      DOUBLE PRECISION :: VCHC3,VCHC3REL
      DOUBLE PRECISION :: VLHC3,VLHC3REL
      DOUBLE PRECISION :: VLLHC3,VLLHC3REL
      DOUBLE PRECISION :: VCCHC3,VCCHC3REL
      DOUBLE PRECISION :: VLC1HC3,VLC1HC3REL
      DOUBLE PRECISION :: VLC2HC3,VLC2HC3REL
      DOUBLE PRECISION :: VLC3H,VLC3HREL
      DOUBLE PRECISION :: VCC3H,VCC3HREL
      DOUBLE PRECISION :: VLC2HC,VLC2HCREL
      DOUBLE PRECISION :: VC2CH,VC2CHREL
      DOUBLE PRECISION :: VC2HCMIN,VC2HCMINREL
      DOUBLE PRECISION :: VC2HCTS,VC2HCTSREL
      DOUBLE PRECISION :: VC2CPLUSH,VC2CPLUSHREL
      DOUBLE PRECISION :: VCPLUSCCH,VCPLUSCCHREL
      DOUBLE PRECISION :: VHC2CMIN,VHC2CMINREL 
      DOUBLE PRECISION :: VHC2CTS,VHC2CTSREL 
      DOUBLE PRECISION :: VDALL,VDALLREL
      DOUBLE PRECISION :: Eh2J
      LOGICAL :: DER
      DOUBLE PRECISION, DIMENSION(6) :: DVDR
      
      DER=.TRUE.
C
C #################################################
C
C     A.U. TO KJ MOL-1 CONVERSION FACTOR
C
      Eh2J=627.5095D+00*4.184D+00
C
C #################################################
C
C     EVALUATING STATIONARY POINTS
C
C          COORDINATES IN BOHR
C          ENERGIES IN HARTREES
C
C #################################################
C
      WRITE(*,*)
      WRITE(*,*) "#################################################"
      WRITE(*,*) "                 VALENCE REGION"
      WRITE(*,*) "#################################################"
      WRITE(*,*)
C
C     C-HC3 GLOBAL MINIMUM - C2V - STRUCTURE #1
C 
      R(1)=4.47004D+00
  
      R(2)=2.03390D+00

      R(3)=4.47004D+00

      R(4)=2.58733D+00

      R(5)=2.55689D+00

      R(6)=2.58733D+00

      CALL POTHC3(R,VCHC3,DER,DVDR)

      VCHC3REL=(VCHC3-VCHC3)*Eh2J

      WRITE(*,FMT='(A,F17.12,4X,A,F17.12)')"V(C-HC3/Eh):", VCHC3,
     *     "RELATIVE ENERGY (kJ mol-1):",VCHC3REL
C
C     L-HC3 LOCAL MINIMUM - CINFV - STRUCTURE #51
C
      R(1)=4.34015D+00
  
      R(2)=6.85893D+00

      R(3)=1.99428D+00

      R(4)=2.51878D+00

      R(5)=2.34587D+00

      R(6)=4.86465D+00

      CALL POTHC3(R,VLHC3,DER,DVDR)

      VLHC3REL=(VLHC3-VCHC3)*Eh2J

      WRITE(*,FMT='(A,F17.12,4X,A,F17.12)')"V(L-HC3/Eh):", VLHC3,
     *     "RELATIVE ENERGY (kJ mol-1):",VLHC3REL
C
C     CC-HC3 ISOMERIZATION TS - C2V - STRUCTURE #39
C
      R(1)=3.85202D+00
  
      R(2)=2.36153D+00

      R(3)=2.36153D+00

      R(4)=2.54009D+00

      R(5)=2.54009D+00

      R(6)=3.02789D+00

      CALL POTHC3(R,VCCHC3,DER,DVDR)

      VCCHC3REL=(VCCHC3-VCHC3)*Eh2J

      WRITE(*,FMT='(A,F17.12,4X,A,F17.12)')"V(CC-HC3/Eh):", VCCHC3,
     *     "RELATIVE ENERGY (kJ mol-1):",VCCHC3REL
C
C     LL-HC3 ISOMERIZATION TS - C2V - STRUCTURE #35
C
      R(1)=2.22565D+00
  
      R(2)=3.34579D+00

      R(3)=3.34579D+00

      R(4)=2.53297D+00

      R(5)=2.53297D+00

      R(6)=5.06532D+00

      CALL POTHC3(R,VLLHC3,DER,DVDR)

      VLLHC3REL=(VLLHC3-VCHC3)*Eh2J

      WRITE(*,FMT='(A,F17.12,4X,A,F17.12)')"V(LL-HC3/Eh):", VLLHC3,
     *     "RELATIVE ENERGY (kJ mol-1):",VLLHC3REL      
C
C     LC1-HC3 ISOMERIZATION TS - C1 - STRUCTURE #43
C
      R(1)=4.47526D+00
  
      R(2)=2.03643D+00

      R(3)=5.21169D+00

      R(4)=2.44583D+00

      R(5)=2.58646D+00

      R(6)=3.68223D+00

      CALL POTHC3(R,VLC1HC3,DER,DVDR)

      VLC1HC3REL=(VLC1HC3-VCHC3)*Eh2J

      WRITE(*,FMT='(A,F17.12,4X,A,F17.12)')"V(LC-HC3/Eh):", VLC1HC3,
     *     "RELATIVE ENERGY (kJ mol-1):",VLC1HC3REL 

      WRITE(*,*)
      WRITE(*,*) "#################################################"
      WRITE(*,*) "              LONG RANGE REGION"
      WRITE(*,*) "#################################################"
      WRITE(*,*)
C
C     HC2...C LR MIN - C1 - STRUCTURE #6
C
      R(1)=6.63161D+00
  
      R(2)=2.00832D+00

      R(3)=4.28940D+00

      R(4)=6.45948D+00

      R(5)=6.96017D+00

      R(6)=2.28124D+00

      CALL POTHC3(R,VHC2CMIN,DER,DVDR)

      VHC2CMINREL=(VHC2CMIN-VCHC3)*Eh2J

      WRITE(*,FMT='(A,F17.12,4X,A,F17.12)') 
     *      "V(HC2...C[MIN]/Eh):", VHC2CMIN,
     *      "RELATIVE ENERGY (kJ mol-1):",VHC2CMINREL
C
C     HC2...C LR TS - CS - STRUCTURE #24
C
      R(1)=2.00420D+00
  
      R(2)=5.98403D+00

      R(3)=4.28529D+00

      R(4)=5.77283D+00

      R(5)=2.31696D+00

      R(6)=6.87261D+00

      CALL POTHC3(R,VHC2CTS,DER,DVDR)

      VHC2CTSREL=(VHC2CTS-VCHC3)*Eh2J

      WRITE(*,FMT='(A,F17.12,4X,A,F17.12)') 
     *      "V(HC2...C[TS]/Eh):", VHC2CTS,
     *      "RELATIVE ENERGY (kJ mol-1):",VHC2CTSREL
C
C     C2...HC LR MIN - C2V - STRUCTURE #64
C
      R(1)=5.46681D+00
  
      R(2)=5.46681D+00

      R(3)=2.10413D+00

      R(4)=2.43883D+00

      R(5)=7.53257D+00

      R(6)=7.53257D+00

      CALL POTHC3(R,VC2HCMIN,DER,DVDR)

      VC2HCMINREL=(VC2HCMIN-VCHC3)*Eh2J

      WRITE(*,FMT='(A,F17.12,4X,A,F17.12)') 
     *      "V(C2...HC[MIN]/Eh):", VC2HCMIN,
     *      "RELATIVE ENERGY (kJ mol-1):",VC2HCMINREL
C
C     C2...HC LR TS - C2V - STRUCTURE #59
C
      R(1)=6.35446D+00
  
      R(2)=6.35446D+00

      R(3)=2.11772D+00

      R(4)=2.44293D+00

      R(5)=7.10808D+00

      R(6)=7.10808D+00

      CALL POTHC3(R,VC2HCTS,DER,DVDR)

      VC2HCTSREL=(VC2HCTS-VCHC3)*Eh2J

      WRITE(*,FMT='(A,F17.12,4X,A,F17.12)') 
     *      "V(C2...HC[TS]/Eh):", VC2HCTS,
     *      "RELATIVE ENERGY (kJ mol-1):",VC2HCTSREL

      WRITE(*,*)
      WRITE(*,*) "#################################################"
      WRITE(*,*) "              ASYMPTOTIC CHANNELS"
      WRITE(*,*) "#################################################"
      WRITE(*,*)
C
C     L-C3 + H ASYMPTOTIC CHANNEL
C
      R(1)=1000.0000D+00
  
      R(2)=1002.445D+00

      R(3)=1004.890D+00

      R(4)=2.445D+00

      R(5)=4.890D+00

      R(6)=2.445D+00

      CALL POTHC3(R,VLC3H,DER,DVDR)

      VLC3HREL=(VLC3H-VCHC3)*Eh2J

      WRITE(*,FMT='(A,F17.12,4X,A,F17.12)')"V(L-C3+H/Eh):", VLC3H,
     *     "RELATIVE ENERGY (kJ mol-1):",VLC3HREL
C
C     C-C3 + H ASYMPTOTIC CHANNEL
C
      R(1)=1000.0000D+00
  
      R(2)=1000.0000D+00

      R(3)=1000.0000D+00

      R(4)=2.580D+00

      R(5)=2.580D+00

      R(6)=2.580D+00

      CALL POTHC3(R,VCC3H,DER,DVDR)

      VCC3HREL=(VCC3H-VCHC3)*Eh2J

      WRITE(*,FMT='(A,F17.12,4X,A,F17.12)')"V(C-C3+H/Eh):", VCC3H,
     *     "RELATIVE ENERGY (kJ mol-1):",VCC3HREL
C
C     L-C2H + C ASYMPTOTIC CHANNEL
C
      R(1)=2.009D+00
  
      R(2)=4.294D+00

      R(3)=1004.294D+00

      R(4)=2.285D+00

      R(5)=1002.285D+00

      R(6)=1000.00000D+00

      CALL POTHC3(R,VLC2HC,DER,DVDR)

      VLC2HCREL=(VLC2HC-VCHC3)*Eh2J

      WRITE(*,FMT='(A,F17.12,4X,A,F17.12)')"V(L-C2H+C/Eh):", VLC2HC,
     *     "RELATIVE ENERGY (kJ mol-1):",VLC2HCREL
C
C     C2 + CH ASYMPTOTIC CHANNEL
C
      R(1)=2.12021D+00
  
      R(2)=1002.12021D+00

      R(3)=1004.47679D+00

      R(4)=1000.0D+00

      R(5)=1002.35658D+00

      R(6)=2.35658D+00

      CALL POTHC3(R,VC2CH,DER,DVDR)

      VC2CHREL=(VC2CH-VCHC3)*Eh2J

      WRITE(*,FMT='(A,F17.12,4X,A,F17.12)')"V(C2+CH/Eh):", VC2CH,
     *     "RELATIVE ENERGY (kJ mol-1):",VC2CHREL
C
C     C2 + C + H ASYMPTOTIC CHANNEL
C
      R(1)=1000.0D+00
  
      R(2)=2000.0D+00

      R(3)=2002.35658D+00

      R(4)=1000.0D+00

      R(5)=1002.35658D+00

      R(6)=2.35658D+00

      CALL POTHC3(R,VC2CPLUSH,DER,DVDR)

      VC2CPLUSHREL=(VC2CPLUSH-VCHC3)*Eh2J

      WRITE(*,FMT='(A,F17.12,4X,A,F17.12)')"V(C2+C+H/Eh):", VC2CPLUSH,
     *     "RELATIVE ENERGY (kJ mol-1):",VC2CPLUSHREL
C
C     C + C + CH ASYMPTOTIC CHANNEL
C
      R(1)=2.12021D+00
  
      R(2)=1002.12021D+00

      R(3)=2002.12021D+00

      R(4)=1000.0D+00

      R(5)=2000.0D+00

      R(6)=1000.0D+00

      CALL POTHC3(R,VCPLUSCCH,DER,DVDR)

      VCPLUSCCHREL=(VCPLUSCCH-VCHC3)*Eh2J

      WRITE(*,FMT='(A,F17.12,4X,A,F17.12)')"V(CH+C+C/Eh):", VCPLUSCCH,
     *     "RELATIVE ENERGY (kJ mol-1):",VCPLUSCCHREL

      WRITE(*,*)
      WRITE(*,*) "#################################################"

      END PROGRAM TEST
