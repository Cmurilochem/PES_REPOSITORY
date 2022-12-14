C ** General program to find stationary points using the Newton-Raphson method.
C ** Single calculation or random search in projected cartesian coordinates.
C ** Internal (Z-matrix) coordinates are obtained as well as the harmonic 
C ** frequencies. LaTex and Molden files are created with the relevant data.
C **  S?rgio Rodrigues 25/5/2002
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(NA=20,NT=3*NA,NFM=1000000,NTES=40000)
      PARAMETER(PREC=1.0D-4,PRECE=1.0D-6,PRECA=1.0D-3)
      DIMENSION X(NT),R(NT),RT(NT),SOL(NT,NTES),VALP(NT)
      DIMENSION VECP(NT,NT),RMM(NA),XU(NT),XT(NT),VDER(NT),RO(NT)
      CHARACTER*2 ATOM(NA)
      DIMENSION NB(0:10,10),LI(0:20),NSUM(0:20),N11(0:NT),N12(0:NT)
      DIMENSION N13(0:NT),N14(0:NT),N15(0:NT)
      DIMENSION NFACT(NT),NPER(NT,NFM),N9(0:NT),N10(0:NT),XO(NT)
      LOGICAL SOBREP,LOG1(NT),LOG2(NT),LOG1P(NT),LOG2P(NT)
      LOGICAL SOMEA,TLINE,TPLAN,CONTIN,MONKEY
C ********************************************************************
C ** The point group is obtained with a system calling to the program 
C **  of Serguei Patchkovskii (http://www.cobalt.chem.ucalgary.ca/ps/)
C ********************************************************************
      character*5 pointg,GS(58)
      INTEGER H(58)
      DATA H/1,2,2,2,3,4,5,6,7,8,4,
     &       6,8,10,12,14,16,4,6,8,10, ! verificar D7 e D8
     &       12,14,16,4,6,8,10,12,14,  ! verificar C7v, C8v e C7h
     &       16,8,12,16,20,24,30,34,8, ! verificar D7h e D8h 
     &       12,16,20,24,28,32,4,6,8,     ! verificar D6d e D7d
     &       12,24,24,24,48,60,120,0,0,0/
      DATA GS/"C1","Cs","Ci","C2","C3","C4","C5","C6","C7","C8","D2",
     &        "D3","D4","D5","D6","D7","D8","C2v","C3v","C4v","C5v",
     &        "C6v","C7v","C8v","C2h","C3h","C4h","C5h","C6h","C7h",
     &        "C8h","D2h","D3h","D4h","D5h","D6h","D7h","D8h","D2d",
     &        "D3d","D4d","D5d","D6d","D7d","D8d","S4","S6","S8", 
     &        "T","Th","Td","O","Oh","I","Ih","Cinfv","Dinfh","Kh"/
C ********************************************************************
      DATA N1,N2,N3,N4,N5,N6,N7,N8/8*0/
      DATA N9,N10,N11,N12/244*0/
      DATA N13,N14,N15/183*0/
      MONKEY=.FALSE.
      OPEN(5,FILE='min.dat',STATUS='unknown')
C ** reading the number of atoms
      READ(5,*)NATOM
C ** number of internuclear distances
      NRD=NATOM*(NATOM-1)/2
      NTA=3*NATOM
      NINT=NTA-6
      RNA=REAL(NATOM)
      DO I=1,NATOM
C ** read atomic symbol (character*2) and atomic mass
      print*,I
      READ(5,*)ATOM(I),RMM(I)
      ENDDO
C ** read DRMAX (maximum Delta X)
C ** NTEST (number of initial tentatives)
C ** IWRIT>0 writes a lot ...
      READ(5,*)NTEST,NFAILM,IWRIT,IENER,NCICLE
      READ(5,*)DRMAX,RSMIN,RSMAX,TEMP,DRMIN,RMAX,RMIN
      READ(5,*)(X(I),I=1,NTA)
      CONTIN=.FALSE.    
C ** channel with previous data if exists!
      OPEN(77,ERR=1451,STATUS='OLD')
      CONTIN=.TRUE.
 1451 IF(CONTIN)THEN
      OPEN(33,STATUS='old',ACCESS='APPEND')
      OPEN(78,STATUS='old',ACCESS='APPEND')
      OPEN(43,FILE='min.res',STATUS='old',ACCESS='APPEND')
      OPEN(44,FILE='min.tex',STATUS='old',ACCESS='APPEND')
      OPEN(45,FILE='min.tab',STATUS='old',ACCESS='APPEND')
      OPEN(98,FILE='min.molden',STATUS='old',ACCESS='APPEND')
      DO I=0,NTA-6
      OPEN(50+I,STATUS='old',ERR=1452,ACCESS='APPEND')
      GOTO 1453
 1452 OPEN(50+I,STATUS='new')
 1453 ENDDO
      READ(77,*)NSTEPS,KTES,N1,N2,N3,N4,N5,N6,N7,N8
      READ(77,*)N9,N10,N11,N12,N13,N14,N15
      READ(77,*)((SOL(J,K),J=1,NRD),K=1,KTES)
      ELSE
      OPEN(43,FILE='min.res',STATUS='unknown')
      OPEN(44,FILE='min.tex',STATUS='unknown')
      OPEN(45,FILE='min.tab',STATUS='unknown')
      OPEN(98,FILE='min.molden',STATUS='unknown')
      KTES=0
      ENDIF
      WRITE(43,*)' ***************************************************'
      WRITE(43,*)'  NUMBER OF ATOMS=',NATOM
      WRITE(43,*)'  NUMBER OF MONTE CARLO CYCLES=',NTEST  
      WRITE(43,*)'  NUMBER OF ACCEPTED FAILURES=',NFAILM  
      WRITE(43,*)'  NUMBER OF STARTING COORDINATES=',NCICLE
      WRITE(43,*)'  IWRIT=',IWRIT  
      WRITE(43,*)'  IENER=',IENER  
      WRITE(43,'(A,F10.4)')'  TEMPERATURE=',TEMP  
      WRITE(43,'(A,F10.4)')'  DRMAX=',DRMAX  
      WRITE(43,'(A,F10.4)')'  RSMIN=',RSMIN  
      WRITE(43,'(A,F10.4)')'  RSMAX=',RSMAX  
      WRITE(43,'(A,F10.4)')'  RMIN =',RMIN  
      WRITE(43,'(A,F10.4)')'  RMAX =',RMAX  
      WRITE(43,'(A,F10.4)')' INITIAL CARTESIAN COORDINATES'  
      WRITE(43,'(6F12.6)')(X(I),I=1,NTA)
      CALL TRANS(X,R,NTA)          
      E=POTEN(R,NRD)
      print*,(R(I),I=1,3),E
      WRITE(43,*)(R(I),I=1,3)
      WRITE(43,'(A,F20.4)')'  INITIAL ENERGY =',E  
      WRITE(43,*)'                                                    '
      WRITE(43,*)' ***************************************************'
      IF(IWRIT.GT.0)THEN
      WRITE(99,'(A,12F10.4)')'initial R:',(R(I),I=1,NRD)
      WRITE(99,'(A,F20.4)')'initial E:',E
      ENDIF
      DO I=1,NTA
      XU(I)=X(I)
      ENDDO
      EANT=E
C ** continue a previous calculation ?
      IF(.NOT.CONTIN)NSTEPS=0
c ** calculate factorials and permutations
      CALL FACTOR(NATOM,NFACT)
      CALL PERMUT(NATOM,NPER,NFM)
c ** calculate valid permutations for the problem
c ** (can be improved if included in the permutations calculation)
      NPMAX=0
      DO J=1,NFACT(NATOM)
      SOMEA=.TRUE.
      DO I=1,NATOM
      IF(ATOM(I).NE.ATOM(NPER(I,J)))THEN
       SOMEA=.FALSE.
       GOTO 106
      ENDIF
      ENDDO
 106  IF(SOMEA)THEN
      NPMAX=NPMAX+1
      DO I=1,NATOM
      NPER(I,NPMAX)=NPER(I,J)
      ENDDO
      ENDIF
      ENDDO
      WRITE(43,*)' NUMBER OF VALID PERMUTATIONS',NPMAX
C **
C ** use some different starting positions on the system
C **
      DO KIC=1,NCICLE
C **
      NACEPT=0
      NFAIL=0
      IF(KIC.EQ.1)GOTO 168
      WRITE(43,*)' CICLE', KIC,' INITIAL CARTESIAN COORDINATES'  
      READ(5,*)(XU(IKJL),IKJL=1,NTA)
      WRITE(43,'(6F12.6)')(XU(IKLJ),IKLJ=1,NTA)
      CALL TRANS(XU,R,NTA)          
      EANT=POTEN(R,NRD)
      WRITE(43,'(A,F20.4)')'  INITIAL ENERGY =',EANT  
      WRITE(43,*)'                                                    '
      WRITE(43,*)' ***************************************************'
      IF(IWRIT.GT.0)THEN
      WRITE(99,'(A,12F10.4)')'initial R:',(R(IKLJ),IKLJ=1,NRD)
      WRITE(99,'(A,F20.4)')'initial E:',E
      ENDIF
 168  CONTINUE
C **
C ** start random walk steps
C **
      DO K=1,NTEST
C ** random walk on the cartesian space
C **
      DO KR=1,NTA
C **
       XUANT=XU(KR)
C ** generate random movement of a coordinate
       DELTAX=DRMAX*(1.0D0-2.0D0*RANFF(DUM))
       XU(KR)=XU(KR)+DELTAX
       CALL TRANS(XU,R,NTA)
C ** try to confine the particles to a bag using an elastic wall...
       DO IKR=1,NRD
       IF(R(IKR).GT.RMAX)THEN
       XU(KR)=XUANT-DELTAX
       GOTO 193
       ENDIF
       ENDDO
 193   E=POTEN(R,NRD)
       DELTAE=(E-EANT)/TEMP
       IF(E.LT.EANT)THEN
       NACEPT=NACEPT+1
       NSTEPS=NSTEPS+1
       EANT=E
       ELSE IF(EXP(-DELTAE).GT.RANFF(DUM))THEN    
       NACEPT=NACEPT+1
       NSTEPS=NSTEPS+1
       EANT=E
       ELSE
       NSTEPS=NSTEPS+1
       XU(KR)=XUANT
       IF(IWRIT.GT.0)THEN
       WRITE(99,'(2I5,A,F20.5)')K,KR,' GEOMETRY NOT ACCEPTED',E
       ENDIF
       GOTO 10
       ENDIF
       IF(IWRIT.GT.0)THEN
       WRITE(99,'(2I5,12F8.4)')K,KR,(R(KL),KL=1,NRD),E
      ENDIF
C ** verify sobreposition with previous stationary points
      DO KT=1,KTES
      SOBREP=.TRUE. 
C ** test for nearby configuration
      DO KTR=1,NRD
C ** This is beeing improved in order to become semething as
C ** a taboo search... must be taken out the zone!!!
      IF(ABS(R(KTR)-SOL(KTR,KT)).GT.DRMIN)THEN
      SOBREP=.FALSE.
      GOTO 102
      ENDIF
      ENDDO
      KTSOBRE=KT
      GOTO 103
 102  ENDDO
 103  IF(SOBREP)THEN
      IF(IWRIT.GT.0)THEN
      WRITE(99,*)' SOBREPOSITION WITH ',KTSOBRE
      ENDIF
      NFAIL=NFAIL+1
      GOTO 10
      ENDIF
      DO I=1,NTA
      X(I)=XU(I)
      ENDDO
C ** keeps starting configuration to actualize the restrictions
      DO I=1,NRD
      RT(I)=R(I)
      ENDDO
C ** start geometry optimization
      N1=N1+1
      CALL NR(X,IND,VALP,VECP,RMM,NATOM,IWRIT,N2,N4,N5,N6,RMIN,RMAX)
C ** if fails jumps to new a new cicle
      IF(IND.EQ.100) THEN
      NFAIL=NFAIL+1
      GOTO 10
      ENDIF
      N7=N7+1
      CALL TRANS(X,R,NTA)
      E=POTEN(R,NRD)
      IF(IWRIT.GT.0)THEN
       WRITE(99,'(A,12F10.4)')'found R:',(R(I),I=1,NRD)
       WRITE(99,'(A,F20.4)')'  with E:',E
      ENDIF
C ** test for repetition of configuration
      DO KT=1,KTES
      DO I=1,NRD
      LOG1(I)=.TRUE.
      LOG2(I)=.TRUE.
      LOG1P(I)=.TRUE.
      LOG2P(I)=.TRUE.
      ENDDO
      NDIST=0
      NDISP=0
      DO I=1,NRD
      DO J=1,NRD
      ADELR=ABS(R(I)-SOL(J,KT))
      IF(ADELR.LT.RSMIN.AND.ADELR.GT.PREC.AND.LOG1P(I).AND.LOG2P(J))THEN 
      LOG1P(I)=.FALSE.
      LOG2P(J)=.FALSE.
      NDISTP=NDISTP+1
      ENDIF
      IF(ADELR.LT.PREC.AND.LOG1(I).AND.LOG2(J))THEN 
      LOG1(I)=.FALSE.
      LOG2(J)=.FALSE.
      NDIST=NDIST+1
      GOTO 105
      ENDIF
      ENDDO
 105  ENDDO
C **
      IF(NDIST.EQ.NRD)THEN
      IF(IWRIT.GT.0)THEN
      WRITE(99,*)' REPETED STATIONARY ',KT
      ENDIF
      NFAIL=NFAIL+1
      GOTO 10
      ENDIF
      ENDDO
C **
c      DO I=1,NTA
c      XU(I)=X(I)
c      ENDDO

      IF(NDISTP.EQ.NRD)THEN
      WRITE(43,*)'                                                  '
      WRITE(43,*)' *************************************************'
      WRITE(43,'(A,6F10.4)')' WARNING:',(R(IKR),IKR=1,NRD)
      WRITE(43,*)' IS LIKELY TO BE A MONKEY SADDLE POINT'
      MONKEY=.TRUE.
      GOTO 10
      ENDIF
C ** if no repetition keeps the stationary
      IF(NDIST.LT.NRD)THEN
      KTES=KTES+1
      IF(IWRIT.GT.0)THEN
      WRITE(99,'(2I5)')NFAIL,KTES
      ENDIF
      NFAIL=0
      NN=0
      NP=0
      DO KTR=1,NTA
      IF(ABS(VALP(KTR)).GT.PRECE)THEN
       IF(VALP(KTR).LT.0.0D0)THEN
        NN=NN+1
       ELSE
        NP=NP+1
       ENDIF
      ENDIF
      ENDDO
      DO KTR=1,NRD
       SOL(KTR,KTES)=R(KTR)
      ENDDO
C **
C ** test for linear configuration
C **
      TLINE=.TRUE.
      DO I=1,NATOM
      DO J=1,NATOM
      DO L=1,NATOM
      IF(I.NE.J.AND.I.NE.L.AND.J.NE.L)THEN
      I1=3*(I-1)+1
      I2=3*(I-1)+2
      I3=3*(I-1)+3
      J1=3*(J-1)+1
      J2=3*(J-1)+2
      J3=3*(J-1)+3
      K1=3*(L-1)+1
      K2=3*(L-1)+2
      K3=3*(L-1)+3
      RAZ1=(X(J1)-X(I1))
      RAD1=(X(K1)-X(I1))
      RAZ2=(X(J2)-X(I2))
      RAD2=(X(K2)-X(I2))
      RAZ3=(X(J3)-X(I3))
      RAD3=(X(K3)-X(I3))
      IF  (ABS(RAZ1*RAD2-RAZ2*RAD1).GT.2*PREC
     & .OR.ABS(RAZ2*RAD3-RAZ3*RAD2).GT.2*PREC
     & .OR.ABS(RAZ1*RAD3-RAZ3*RAD1).GT.2*PREC)THEN
      TLINE=.FALSE.
C **
C ** keeps 3 points not in a line for in-plane search      
C **
      IPL=I
      JPL=J
      KPL=L
      GOTO 109 
      ENDIF
      ENDIF
      ENDDO
      ENDDO
      ENDDO
C **
C ** if not linear test for in-plane configuration
C **
 109  TPLAN=.FALSE.
      IF(.NOT.TLINE.AND.NATOM.GT.3)THEN
      X1=X(3*(IPL-1)+1)
      Y1=X(3*(IPL-1)+2)
      Z1=X(3*(IPL-1)+3)
      X2=X(3*(JPL-1)+1)
      Y2=X(3*(JPL-1)+2)
      Z2=X(3*(JPL-1)+3)
      X3=X(3*(KPL-1)+1)
      Y3=X(3*(KPL-1)+2)
      Z3=X(3*(KPL-1)+3)
C ** parameters of the plane      
      APL=(Y2-Y1)*(Z3-Z1)-(Y3-Y1)*(Z2-Z1)
      BPL=(Z2-Z1)*(X3-X1)-(Z3-Z1)*(X2-X1)
      CPL=(X2-X1)*(Y3-Y1)-(X3-X1)*(Y2-Y1)
      DPL=-APL*X1-BPL*Y1-CPL*Z1
      RNORM=SQRT(APL*APL+BPL*BPL+CPL*CPL)
C ** calculate the distances from the plane
      TPLAN=.TRUE.
      DO I=1,NATOM
      IF(I.NE.IPL.AND.I.NE.JPL.AND.I.NE.KPL)THEN
      I1=3*(I-1)+1
      I2=3*(I-1)+2
      I3=3*(I-1)+3
      DIST=(APL*X(I1)+BPL*X(I2)+CPL*X(I3)+DPL)/RNORM
      IF(ABS(DIST).GT.PREC)THEN
      TPLAN=.FALSE.
      GOTO 110
      ENDIF
      ENDIF
      ENDDO
 110  CONTINUE
      ENDIF
C **
      N8=N8+1
      N9(NN)=N9(NN)+1
      NOPTI=1
C **
      IF(TLINE)THEN
      WRITE(43,*)' ***************************************************'
      WRITE(43,*)' LINEAR CONFIGURATION'
      ELSE IF(TPLAN.AND.NATOM.GT.3) THEN
      WRITE(43,*)' **************************************************'
      WRITE(43,*)' PLANAR CONFIGURATION'
      ELSE IF(NATOM.GT.3)THEN
      WRITE(43,*)' **************************************************'
      WRITE(43,*)' SPATIAL CONFIGURATION'
      ELSE
c ** non-linear 3-atom case
      WRITE(43,*)' *************************************************'
      WRITE(43,*)' NON-LINEAR TRIATOM'
      ENDIF
      WRITE(43,*)' STRUCTURE NUMBER=',KTES,' DRMAX=',DRMAX
      WRITE(43,*)' *************************************************'
C ** for obtaing the point group uses the C program of Serguei Patchkovskii
C ** slightly arranjed in order to use files
      open(34,file='SOURCE/symm.dat',status='unknown')
      WRITE(34,*)NATOM
      DO I=1,NATOM
      I1=3*(I-1)+1
      I2=3*(I-1)+2
      I3=3*(I-1)+3
      WRITE(34,*)ICHAR(ATOM(I)),(X(3*(I-1)+J),J=1,3)
      ENDDO
      close(34)
C ** calling of Serguei Patchkovskii symmetry program
      VPOINTG=0
      CALL SYSTEM("./SOURCE/symm.x < SOURCE/symm.dat")
      open(17,file='symm.res',status='unknown')
      read(17,*,END=2014)pointg
      VPOINTG=1
      close(17)
      DO I=1,58
      IF(GS(I).EQ.POINTG)THEN
      IV=I
      GOTO 134
      ENDIF
      ENDDO
 134  IF(IV.LT.56)THEN
      WRITE(43,*)pointg,'symmetry (finite group)    h=',H(IV)
      write(43,*)' Symmetry number (point-group based)=',2*NPMAX/H(IV) 
      ELSE
      WRITE(43,*) pointg,' symmetry  (continuous group!)'
      ENDIF
 2014 IF(VPOINTG.LT.1)THEN
      WRITE(43,*)'failed to find the point-group'
      close(17)
      ENDIF
C **
      IF(IWRIT.GT.0)THEN
       WRITE(99,*)'SYMMETRY   ',pointg
      ENDIF
C **
      WRITE(43,*)' **************************************************'
      NND=NN+NP
      IF((TLINE.AND.NND.LT.NTA-5).OR.(.NOT.TLINE.AND.NND.LT.NTA-6))THEN
      WRITE(43,*)'                                                  '
      WRITE(43,*)' *************************************************'
      WRITE(43,*)' WARNING: THIS IS LIKELY A DEGENERATE SADDLE POINT'
      MONKEY=.TRUE.
      ENDIF
      CALL WRIT(VALP,ATOM,X,VECP,E,NN,NP,NATOM,KTES,IENER)
      ENDIF
 10   IF(MOD(NSTEPS,100).EQ.0)THEN
      ACEPT=DBLE(NACEPT)/DBLE(NSTEPS)
      IF(ACEPT.GT.0.5D0)THEN
      DRMAX=MIN(1.2D0*DRMAX,RSMAX)
      ELSE
      DRMAX=MAX(0.8D0*DRMAX,RSMIN)
      ENDIF
      ENDIF
      IF((MONKEY.AND.KTES.GT.5000).OR.KTES.GT.40000)GOTO 104
      IF(NFAIL.GT.NFAILM)THEN
      WRITE(43,'(A,F10.4)')'  FINAL OF CYCLE CARTESIAN COORDINATES'  
      WRITE(43,'(6F12.6)')(XU(I),I=1,NTA)
      CALL TRANS(XU,R,NTA)          
      E=POTEN(R,NRD)
      WRITE(43,'(A,F20.4)')'  FINAL OF CYCLE ENERGY =',E  
      GOTO 174
      ENDIF
      ENDDO
      ENDDO
 174  ENDDO
 104  RNT=DBLE(NSTEPS)
      RNR=DBLE(N1)
      WRITE(43,'(A,F10.4)')'  FINAL CARTESIAN COORDINATES'  
      WRITE(43,'(6F12.6)')(XU(I),I=1,NTA)
      CALL TRANS(XU,R,NTA)          
      E=POTEN(R,NRD)
      WRITE(43,'(A,F20.4)')'  FINAL ENERGY =',E  
       WRITE(43,*)'***************************************************'
       WRITE(43,*)'** Statistical analysis **'
       WRITE(43,'(A,I8,F10.4,A)')'Monte Carlo steps',NSTEPS
       WRITE(43,'(A,I8,F10.4,A)')'NR tries ',N1,N1/RNT*100.0,'%'
       WRITE(43,'(A,I8,F10.4,A)')'NR fail ] R [ ',N2,N2/RNR*100.0,'%'
       WRITE(43,'(A,I8,F10.4,A)')'NR fail Numeric',N4,N4/RNR*100.0,'%'
       WRITE(43,'(A,I8,F10.4,A)')'NR fail converg',N5,N5/RNR*100.0,'%'
       WRITE(43,'(A,I8,F10.4,A)')'Stationary found',N7,N7/RNR*100.0,'%'  
      
      CLOSE(43)
      CLOSE(44)
      CLOSE(45)
      CLOSE(78)
      CLOSE(98)

      IF(CONTIN)REWIND 77
      WRITE(77,*)NSTEPS,KTES,N1,N2,N3,N4,N5,N6,N7,N8
      WRITE(77,*)N9,N10,N11,N12,N13,N14,N15
      WRITE(77,*)((SOL(J,K),J=1,NRD),K=1,KTES)
     
      CLOSE(77)

      END

      SUBROUTINE FACTOR(N,NF)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(NA=20,NT=3*NA)
      INTEGER NF(NT)
      NF(1)=1
      DO I=2,N
      NF(I)=NF(I-1)*I
      ENDDO
      RETURN
      END

      SUBROUTINE PERMUT(N,NPER,NFM)
      PARAMETER(NA=20,NT=3*NA)
      DIMENSION NPER(NT,NFM)
      integer a(N)
      logical mtc,even
      mtc=.false.
      J=1
  10  call nexper(n,a,mtc,even)
      DO I=1,N
      NPER(I,J)=a(i)
      ENDDO
      J=J+1
      if(mtc)goto 10
      RETURN
      end

      subroutine nexper(n,a,mtc,even)
c next permutation of {1,...,n}. Ref. NW p. 59
C ** Nijenhuis and Wilf Combinatorial Algorithms, Academic Press, 1978
C ** http://www.cs.sunysb.edu/~algorith/implement/wilf/implement.shtml
C **
      integer a(n),s,d
      logical mtc,even
      if(mtc)goto 10
      nm3=n-3
      do 1 i=1,n
    1 a(i)=i
      mtc=.true.
    5 even=.true.
      if(n.eq.1)goto 8
    6 if(a(n).ne.1.or.a(1).ne.2+mod(n,2))return
      if(n.le.3)goto 8
      do 7 i=1,nm3
      if(a(i+1).ne.a(i)+1)return
    7 continue
    8 mtc=.false.
      return
   10 if(n.eq.1)goto 27
      if(.not.even)goto 20
      ia=a(1)
      a(1)=a(2)
      a(2)=ia
      even=.false.
      goto 6
   20 s=0
      do 26 i1=2,n
   25 ia=a(i1)
      i=i1-1
      d=0
      do 30 j=1,i
   30 if(a(j).gt.ia) d=d+1
      s=d+s
      if(d.ne.i*mod(s,2)) goto 35
   26 continue
   27 a(1)=0
      goto 8
   35 m=mod(s+1,2)*(n+1)
      do 40 j=1,i
      if(isign(1,a(j)-ia).eq.isign(1,a(j)-m))goto 40
      m=a(j)
      l=j
   40 continue
      a(l)=ia
      a(i1)=m
      even=.true.
      return
      end


      SUBROUTINE WRIT(VALP,ATOM,X,VECP,E,NN,NP,NATOM,KTES,IENER)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(NA=20,NT=3*NA,NTES=40000)
      PARAMETER(PREC=1.0D-4,PRECE=1.0D-6,PRECM=0.1)
      PARAMETER(ECM1=219474.62492D0,ATMAS=1.82288853D3)
      DIMENSION X(NT),R(NT),SOL(NT,NTES),VAL(NT),VALP(NT)
      DIMENSION CINT(NT),VECP(NT,NT),RMM(NA),FREQ(NT)
      CHARACTER*1 IFR(NT)
      CHARACTER*2 ATOM(NA)
      NRD=NATOM*(NATOM-1)/2
      NTA=3*NATOM
      IF(IENER.EQ.0)THEN
C ** frequencies in cm-1
      CONEN=ECM1/SQRT(ATMAS)
      ELSE
C ** frequencies in reduced units 
      CONEN=1.0
      ENDIF
C ** writing in MOLDEN format for visualization of normal modes
      WRITE(98,*)'[Molden Format]                 ',KTES
      WRITE(98,*)'               '
      WRITE(98,*)'[FREQ]'
      DO K=1,NTA
      WRITE(98,*)SQRT(ABS(VALP(K)))*CONEN*SIGN(1.0,VALP(K))
      ENDDO
      WRITE(98,*)'               '
      WRITE(98,*)'[FR-COORD]'
      DO IN=1,NATOM
      NK=3*(IN-1)
      IF(IENER.EQ.0)THEN
      WRITE(98,*)ATOM(IN),(X(K),K=NK+1,NK+3)
      ELSE
C ** the factor 2.5 is for cheating MOLDEN!...
      WRITE(98,*)ATOM(IN),(X(K)*2.5,K=NK+1,NK+3)
      ENDIF
      ENDDO
      WRITE(98,*)'               '
      WRITE(98,*)'[FR-NORM-COORD]'
      DO J=1,NTA
      WRITE(98,*)'vibration',J
      WRITE(98,'(3F20.7)')(VECP(K,J),K=1,NTA)
      ENDDO
      WRITE(98,*)'               '
      WRITE(98,*)'               '
      WRITE(98,*)'               '
      WRITE(98,*)'               '
c ** writing the long file of results
        WRITE(43,1005)E,NN,NP
        WRITE(43,*)' INTERNUCLEAR DISTANCES:'
        CALL TRANS(X,R,NTA)
        WRITE(43,1000)(R(K),K=1,NRD)
        RMEDI=0.0
        RMSQR=0.0
        DO IK=1,NRD
        RMEDI=RMEDI+R(IK)/DBLE(NRD)
        RMSQR=RMSQR+R(IK)*R(IK)/DBLE(NRD)
        ENDDO
        WRITE(50+NN,*)KTES,RMEDI,RMSQR-RMEDI**2
        WRITE(33,*)KTES,NN,(R(K),K=1,NRD)
        WRITE(43,*)' HESSIAN EIGENVALUES:'        
        WRITE(43,1000)(VALP(K),K=1,NTA)
        WRITE(43,*)' FREQUENCIES '
        NVIB=0
        DO K=1,NTA
        FTEST=SQRT(ABS(VALP(K)))*CONEN
        IF(FTEST.GT.PRECM.AND.VALP(K).GT.0.0)THEN
        NVIB=NVIB+1
        FREQ(NVIB)=FTEST
        IFR(NVIB)=' '
        ELSE IF(FTEST.GT.PRECM.AND.VALP(K).LT.0.0)THEN
        NVIB=NVIB+1
        FREQ(NVIB)=FTEST
        IFR(NVIB)='i'
        ENDIF
       ENDDO
        WRITE(43,1011)(SQRT(ABS(VALP(K)))*CONEN,K=1,NTA)
        WRITE(43,*)' HESSIAN EIGENVECTORS:'
        DO K=1,NTA        
          WRITE(43,1000)(VECP(K,J),J=1,NTA)
        ENDDO
        WRITE(43,*)' CARTESIAN COORDINATES:'        
        WRITE(43,1000)(X(K),K=1,NTA)
        CALL DERR(X,VAL,CINT,NATOM)
        WRITE(43,*)' INTERNAL COORDINATES:'        
        WRITE(43,1000)(CINT(K),K=1,NTA-6)
        WRITE(43,*)' FIRST DERIVATIVE INTERNAL COORDINATES:'        
        WRITE(43,1000)(VAL(K),K=1,NRD)
C ** file for minimum energy path following
        IF(NN.EQ.1)THEN
          WRITE(78,1000)(X(K),K=1,NTA)
        NNEG=0
        DO K=1,NTA
         NNEG=NNEG+1
         IF(VALP(K).LT.-PRECE)GOTO 124
        ENDDO
 124    WRITE(78,*)NNEG
        DO K=1,NTA
          WRITE(78,1000)(VECP(K,J),J=1,NTA)
        ENDDO
      ENDIF
C ** writing the tex and table files
      IF(NRD.EQ.6)THEN
      WRITE(44,2006)(R(K),K=1,NRD),E,(FREQ(K),IFR(K),K=1,NVIB)
      WRITE(45,3006)(R(K),K=1,NRD),E,(FREQ(K),IFR(K),K=1,NVIB)
      ELSE IF(NRD.EQ.3)THEN
      WRITE(44,2003)(R(K),K=1,NRD),E,(FREQ(K),IFR(K),K=1,NVIB)
      WRITE(45,3003)(R(K),K=1,NRD),E,(FREQ(K),IFR(K),K=1,NVIB)
      ELSE IF(NRD.EQ.10)THEN
      WRITE(44,2010)(R(K),K=1,NRD),E,(FREQ(K),IFR(K),K=1,NVIB)
      WRITE(45,3010)(R(K),K=1,NRD),E,(FREQ(K),IFR(K),K=1,NVIB)
      ENDIF
 1000 FORMAT(12F12.5)
 1011 FORMAT(12F12.2)
 1005 FORMAT(2X,'ENER=',F14.6,3X,'NEG=',I5,3X,'POS=',I5)
 2003 FORMAT(2X,3(F12.3,'&'),F12.4,'&',4(F10.2,A1,'&'))
 3003 FORMAT(2X,3F12.3,F14.6,4(F10.2,A1))
 2006 FORMAT(2X,6(F12.3,'&'),F12.4,'&',7(F10.2,A1,'&'))
 3006 FORMAT(2X,6F12.3,F12.4,7(F10.2,A1))
 2010 FORMAT(2X,10(F12.3,'&'),F12.4,'&',10(F10.2,A1,'&'))
 3010 FORMAT(2X,10F12.3,F12.4,10(F10.2,A1))
      RETURN
      END

      SUBROUTINE NR(X,IND,VALP,VV,RMM,NATOM,IWR,N2,N4,N5,N6,RMIN,RMAX)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(NA=20,NT=3*NA,PREC=1.0D-4,NITER=100)
      DIMENSION X(NT),R(NT),P(NT,NT),VAL(NT),VV(NT,NT),RMAT(NT,NT)
      DIMENSION VALP(NT),V(NT,NT),CINT(NT),AUX1(NT,NT),VDER(NT)
      DIMENSION RAZ(NT),RM(NT,NT),RMM(NA)
      IND=0
      NTA=3*NATOM
      NRD=NATOM*(NATOM-1)/2
C ** start iterations
      DO KK=1,NITER
C ** calculating the projector of the coordinates
      CALL PROJECT(X,P,IND,NTA)
      IF(IND.EQ.100)THEN
      IF(IWR.GT.1)THEN
      WRITE(99,*)' NR FAIL PROJECTION'
      ENDIF
      N4=N4+1
      RETURN
      ENDIF      
      CALL DER(X,VDER,NTA)
      CALL DDER(X,V,NTA)
      CALL MATP(P,V,AUX1,NTA,NTA,NTA,NT)
      CALL MATP(AUX1,P,V,NTA,NTA,NTA,NT)
      CALL TRANS(X,R,NTA)
      DO I=1,NRD
      IF(R(I).GT.RMAX.OR.R(I).LT.RMIN)THEN
      IF(IWR.GT.1)THEN
      WRITE(99,*)' NR FAIL RMAX OR RMIN ATTAINED'
      WRITE(99,'(I5,12F8.4)')KK,(R(KL),KL=1,NRD)
      ENDIF
      N2=N2+1
      IND=100
      RETURN
      ENDIF
      ENDDO
      DO I=1,NTA
       RMAT(I,NTA+1)=-VDER(I)
       DO J=1,NTA
        RMAT(I,J)=V(I,J)
       ENDDO
      ENDDO
      CALL SVD(RMAT,NTA,RAZ,IND)
      IF(IND.EQ.100)THEN
      IF(IWR.GT.1)THEN
      WRITE(99,*)' NR FAIL SVD FAILED'
      ENDIF
      N4=N4+1
      RETURN
      ENDIF
      NZ1=0
      NZ2=0
      DO I=1,NTA
      IF(ABS(RAZ(I)).LT.PREC)NZ1=NZ1+1
      IF(ABS(VDER(I)).LT.PREC)NZ2=NZ2+1
      X(I)=X(I)+RAZ(I)
      ENDDO
      IF(NZ1.EQ.NTA.AND.NZ2.EQ.NTA)GOTO 20
      ENDDO
      N5=N5+1
      IF(IWR.GT.1)THEN
      WRITE(99,*)' NR FAILED MAXIT ATTAINED'
      ENDIF
      IND=100
      RETURN
 20   IND=0
      N6=N6+KK
      CALL PROJECT(X,P,IND,NTA)
      IF(IND.EQ.100)THEN
      IF(IWR.GT.1)THEN
      WRITE(99,*)' NR FAIL PROJECTION FAILED'
      ENDIF
      N4=N4+1
      RETURN
      ENDIF      
      CALL DDER(X,V,NTA)
      CALL DER(X,VDER,NTA)
      CALL MATP(P,V,AUX1,NTA,NTA,NTA,NT)
      CALL MATP(AUX1,P,V,NTA,NTA,NTA,NT)
      DO I=1,NTA
      DO J=1,NTA
       RM(I,J)=0.0
      ENDDO
      ENDDO      
      DO I=1,NATOM
      DO K=1,3
      J=3*(I-1)+K
      RM(J,J)=1.0/SQRT(RMM(I))
      ENDDO
      ENDDO
      CALL MATP(RM,V,AUX1,NTA,NTA,NTA,NT)
      CALL MATP(AUX1,RM,V,NTA,NTA,NTA,NT)
      CALL JACOBIX(V,NTA,NT,VALP,VV,NROT)
      RETURN
      END

      SUBROUTINE DERR(X,VAL,CINT,NATOM)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NA=20,NT=3*NA,PREC=1.0D-4)
      DIMENSION X(NT),VAL(NT),CINT(NT)
      NRD=NATOM*(NATOM-1)/2
      CALL CARTINT(X,CINT,NATOM)
      DO I=1,NRD
      VAL(I)=DVLR(CINT,PREC,I,NATOM)
      ENDDO
      RETURN
      END

      FUNCTION DVLR(CINT,STEP,I,NATOM)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NA=20,NT=3*NA )
      DIMENSION XA(NT),CINT(NT),XD(NT),R(NT)
      NR=NATOM*(NATOM-1)/2
      NTA=3*NATOM
      DO IA=1,NR
      XD(IA)=CINT(IA)
      ENDDO
      XD(I)=CINT(I)-STEP
      CALL INTCART(XD,XA,NATOM)
      CALL TRANS(XA,R,NTA)  
      FMIN1=POTEN(R,NR)
      XD(I)=CINT(I)-2*STEP
      CALL INTCART(XD,XA,NATOM)
      CALL TRANS(XA,R,NTA)  
      FMIN2=POTEN(R,NR)
      XD(I)=CINT(I)+STEP
      CALL INTCART(XD,XA,NATOM)
      CALL TRANS(XA,R,NTA)
      FMAX1=POTEN(R,NR)
      XD(I)=CINT(I)+2*STEP
      CALL INTCART(XD,XA,NATOM)
      CALL TRANS(XA,R,NTA)
      FMAX2=POTEN(R,NR)
      DVLR=(FMIN2-8.0*FMIN1+8.0*FMAX1-FMAX2)/(12.D0*STEP)
      RETURN
      END

      SUBROUTINE DER(X,VAL,NTA)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NA=20, NT=3*NA,PREC=1.0D-4)
      DIMENSION X(NT),VAL(NT)
      DO I=1,NTA
      VAL(I)=DVL(X,PREC,I,NTA)
      ENDDO
      RETURN
      END

      SUBROUTINE DDER(X,VAL,NTA)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NA=20, NT=3*NA,PREC=1.0D-4)
      DIMENSION X(NT),VAL(NT,NT)
      DO I=1,NTA
      DO J=I,NTA
      VAL(I,J)=D2VL(X,PREC,I,J,NTA)
      ENDDO
      ENDDO
      DO I=1,NTA
      DO J=I,NTA
      VAL(J,I)=VAL(I,J)
      ENDDO
      ENDDO
      RETURN
      END

      FUNCTION DVL(X,STEP,I,NTA)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NA=20, NT=3*NA)
      DIMENSION X(NT),XD(NT)
      DO IN=1,NTA
      XD(IN)=X(IN)
      ENDDO
      XD(I)=X(I)-STEP
      FMIN1=VL(XD,NTA)
      XD(I)=X(I)+STEP
      FMAX1=VL(XD,NTA)
      XD(I)=X(I)-2*STEP
      FMIN2=VL(XD,NTA)
      XD(I)=X(I)+2*STEP
      FMAX2=VL(XD,NTA)
      DVL=(FMIN2-8.0*FMIN1+8.0*FMAX1-FMAX2)/(12.D0*STEP)
      RETURN
      END

      FUNCTION D2VL(X,STEP,I,J,NTA)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NA=20, NT=3*NA)
      DIMENSION X(NT),XD(NT)
      DO IN=1,NTA
      XD(IN)=X(IN)
      ENDDO
      XD(J)=X(J)-STEP
      FMIN1=DVL(XD,STEP,I,NTA)
      XD(J)=X(J)+STEP
      FMAX1=DVL(XD,STEP,I,NTA)
      XD(J)=X(J)-2*STEP
      FMIN2=DVL(XD,STEP,I,NTA)
      XD(J)=X(J)+2*STEP
      FMAX2=DVL(XD,STEP,I,NTA)
      D2VL=(FMIN2-8.0*FMIN1+8.0*FMAX1-FMAX2)/(12.D0*STEP)
      RETURN
      END

      FUNCTION VL(X,NTA)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NA=20,NT=3*NA)
      DIMENSION X(NT),R(NT)
      NR=NTA*(NTA/3-1)/6
      CALL TRANS(X,R,NTA)
      VL=POTEN(R,NR)
      RETURN
      END

      SUBROUTINE TRANS(X,R,NTA)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(NA=20,NT=3*NA)
      DIMENSION X(NT),R(NT)
      K=0
      DO I=1,NTA,3
       DO J=I,NTA,3
        IF(I.NE.J)THEN
         K=K+1
         R(K)=SQRT((X(I)-X(J))**2+(X(I+1)-X(J+1))**2+
     &   (X(I+2)-X(J+2))**2)
        ENDIF
       ENDDO
      ENDDO
      RETURN
      END

      SUBROUTINE CARTINT(X,CINT,NAT)
      IMPLICIT REAL*8 (A-H,O-Z)
C      PARAMETER (NA=20,NT=3*NA,DEGREE=57.29578D0)
      PARAMETER (NA=20,NT=3*NA,DEGREE=1.0)
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
      
      SUBROUTINE INTCART(CINT,X,NAT)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NA=20,NT=3*NA)
      DIMENSION GEO(3,NA),COORD(3,NA)
      DIMENSION X(NT),CINT(NT)
      COMMON /CONECT/NNA(NA),NNB(NA),NNC(NA)
      GEO(1,2)=CINT(1)
      GEO(1,3)=CINT(2)
      GEO(2,3)=CINT(3)
      DO I=4,NAT
        GEO(1,I)=CINT(3*(I-3)+1)
        GEO(2,I)=CINT(3*(I-3)+2)
        GEO(3,I)=CINT(3*(I-3)+3)
      ENDDO
      CALL GMETRY(GEO,COORD,NNA,NNB,NNC,NAT)
      DO I=1,NAT
        X(3*(I-1)+1)=COORD(1,I)
        X(3*(I-1)+2)=COORD(2,I)
        X(3*(I-1)+3)=COORD(3,I)
      ENDDO
      END
      
      SUBROUTINE GMETRY(GEO,COORD,NNA,NNB,NNC,NATOMS)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NA=20)
      DIMENSION GEO(3,NA),COORD(3,NA),NNA(NA),NNB(NA),NNC(NA)
C***********************************************************************
C
C    GMETRY  COMPUTES COORDINATES FROM BOND-ANGLES AND LENGTHS.
C *** IT IS ADAPTED FROM THE PROGRAM WRITTEN BY M.J.S. DEWAR.
C    (C) NORMAL CONVERSION FROM INTERNAL TO CARTESIAN COORDINATESISDONE.
C
C  ON INPUT:
C         GEO    = ARRAY OF INTERNAL COORDINATES.
C         NATOMS = NUMBER OF ATOMS, INCLUDING DUMMIES.
C         NA     = ARRAY OF ATOM LABELS FOR BOND LENGTHS.
C
C  ON OUTPUT:
C         COORD  = ARRAY OF CARTESIAN COORDINATES
C
C***********************************************************************
      COORD(1,1)=0.0D00
      COORD(2,1)=0.0D00
      COORD(3,1)=0.0D00
      COORD(1,2)=GEO(1,2)
      COORD(2,2)=0.0D00
      COORD(3,2)=0.0D00
      IF(NATOMS.EQ.2) RETURN
      CCOS=COS(GEO(2,3))
      IF(NNA(3).EQ.1)THEN
         COORD(1,3)=COORD(1,1)+GEO(1,3)*CCOS
      ELSE
         COORD(1,3)=COORD(1,2)-GEO(1,3)*CCOS
      ENDIF
      COORD(2,3)=GEO(1,3)*SIN(GEO(2,3))
      COORD(3,3)=0.0D00
      DO 90 I=4,NATOMS
         COSA=COS(GEO(2,I))
         MB=NNB(I)
         MC=NNA(I)
         XB=COORD(1,MB)-COORD(1,MC)
         YB=COORD(2,MB)-COORD(2,MC)
         ZB=COORD(3,MB)-COORD(3,MC)
         RBC=XB*XB+YB*YB+ZB*ZB
         IF(RBC.LT.1.D-16)THEN
C
C     TWO ATOMS ARE COINCIDENT.  A FATAL ERROR.
C
            WRITE(6,'(A,I4,A,I4,A)')' ATOMS',MB,' AND',MC,' ARE COINCIDE
     1NT'
            WRITE(6,'(A)')' THIS IS A FATAL ERROR, RUN STOPPED IN GMETRY
     1'
            STOP
         ELSE
            RBC=1.0D00/SQRT(RBC)
         ENDIF
         MA=NNC(I)
         XA=COORD(1,MA)-COORD(1,MC)
         YA=COORD(2,MA)-COORD(2,MC)
         ZA=COORD(3,MA)-COORD(3,MC)
C
C     ROTATE ABOUT THE Z-AXIS TO MAKE YB=0, AND XB POSITIVE.  IF XYB IS
C     TOO SMALL, FIRST ROTATE THE Y-AXIS BY 90 DEGREES.
C
         XYB=SQRT(XB*XB+YB*YB)
         K=-1
         IF (XYB.GT.0.1D00) GO TO 40
         XPA=ZA
         ZA=-XA
         XA=XPA
         XPB=ZB
         ZB=-XB
         XB=XPB
         XYB=SQRT(XB*XB+YB*YB)
         K=+1
C
C     ROTATE ABOUT THE Y-AXIS TO MAKE ZB VANISH
C
   40    COSTH=XB/XYB
         SINTH=YB/XYB
         XPA=XA*COSTH+YA*SINTH
         YPA=YA*COSTH-XA*SINTH
         SINPH=ZB*RBC
         COSPH=SQRT(ABS(1.D00-SINPH*SINPH))
         ZQA=ZA*COSPH-XPA*SINPH
C
C     ROTATE ABOUT THE X-AXIS TO MAKE ZA=0, AND YA POSITIVE.
C
         YZA=SQRT(YPA**2+ZQA**2)
         IF(YZA.LT.1.D-4)GOTO 60
c         IF(YZA.LT.2.D-2)THEN
c            WRITE(6,'(//20X,'' CALCULATION ABANDONED AT THIS POINT'')')
c            WRITE(6,'(//10X,'' THREE ATOMS BEING USED TO DEFINE THE'',/
c     110X,'' COORDINATES OF A FOURTH ATOM, WHOSE BOND-ANGLE IS'')')
c            WRITE(6,'(10X,'' NOT ZERO OR 180 DEGREEES, ARE '',
c     1''IN AN ALMOST STRAIGHT'')')
c            WRITE(6,'(10X,'' LINE.  THERE IS A HIGH PROBABILITY THAT THE
c     1'',/10X,'' COORDINATES OF THE ATOM WILL BE INCORRECT.'')')
c            WRITE(6,'(//20X,''THE FAULTY ATOM IS ATOM NUMBER'',I4)')I
c            STOP
c            WRITE(6,'(//20X,''CARTESIAN COORDINATES UP TO FAULTY ATOM'')
c     1')
c            WRITE(6,'(//5X,''I'',12X,''X'',12X,''Y'',12X,''Z'')')
c            DO 50 J=1,I
c   50       WRITE(6,'(I6,F16.5,2F13.5)')J,(COORD(K,J),K=1,3)
c            WRITE(6,'(//6X,'' ATOMS'',I3,'','',I3,'', AND'',I3,
c     1'' ARE WITHIN'',F7.4,'' ANGSTROMS OF A STRAIGHT LINE'')')
c     2MC,MB,MA,YZA
c            STOP
c         ENDIF
         COSKH=YPA/YZA
         SINKH=ZQA/YZA
         GOTO 70
   60    CONTINUE
C
C   ANGLE TOO SMALL TO BE IMPORTANT
C
         COSKH=1.D0
         SINKH=0.D0
   70    CONTINUE
C
C     COORDINATES :-   A=(???,YZA,0),   B=(RBC,0,0),  C=(0,0,0)
C     NONE ARE NEGATIVE.
C     THE COORDINATES OF I ARE EVALUATED IN THE NEW FRAME.
C
         SINA=SIN(GEO(2,I))
         SIND=-SIN(GEO(3,I))
         COSD=COS(GEO(3,I))
         XD=GEO(1,I)*COSA
         YD=GEO(1,I)*SINA*COSD
         ZD=GEO(1,I)*SINA*SIND
C
C     TRANSFORM THE COORDINATES BACK TO THE ORIGINAL SYSTEM.
C
         YPD=YD*COSKH-ZD*SINKH
         ZPD=ZD*COSKH+YD*SINKH
         XPD=XD*COSPH-ZPD*SINPH
         ZQD=ZPD*COSPH+XD*SINPH
         XQD=XPD*COSTH-YPD*SINTH
         YQD=YPD*COSTH+XPD*SINTH
         IF (K.LT.1) GO TO 80
         XRD=-ZQD
         ZQD=XQD
         XQD=XRD
   80    COORD(1,I)=XQD+COORD(1,MC)
         COORD(2,I)=YQD+COORD(2,MC)
         COORD(3,I)=ZQD+COORD(3,MC)
   90 CONTINUE
      RETURN
      END

      SUBROUTINE XYZINT(XYZ,NNA,NNB,NNC,DEGREE,GEO,NUMAT)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NA=20)
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
      
      SUBROUTINE XYZGEO(XYZ,NUMAT,NNA,NNB,NNC,DEGREE,GEO)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NA=20)
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
      
      SUBROUTINE BANGLE(XYZ,I,J,K,ANGLE)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NA=20)
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

      SUBROUTINE DIHED(XYZ,I,J,K,L,ANGLE)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NA=20)
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
      
      SUBROUTINE MATP(A,B,C,N,L,M,ND)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(ND,L),B(ND,M),C(ND,M)
      DO I=1,N
         DO J=1,M
            C(I,J)=0.0
            DO K=1,L
               C(I,J)=C(I,J)+A(I,K)*B(K,J)
            ENDDO
         ENDDO
      ENDDO
      RETURN
      END


      SUBROUTINE PROJECT(X,P,IND,NTA)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NA=20,NT=3*NA)
      DIMENSION X(NT),B(NT,NT),S(NT,NT),SINV(NT,NT),
     &  AUX1(NT,NT),BTR(NT,NT),P(NT,NT),U(NT,NT)
      DO I=1,NTA
      DO J=1,6
      B(I,J)=0.0
      ENDDO
      ENDDO
      DO I=1,NTA
      DO J=1,NTA
      U(I,J)=0.0
      ENDDO
      U(I,I)=1.0
      ENDDO
      DO I=1,NTA,3
      B(I,  1)=1.0
      B(I+1,2)=1.0
      B(I+2,3)=1.0
      B(I+1,4)= X(I+2)
      B(I+2,4)=-X(I+1)
      B(I,  5)=-X(I+2)
      B(I+2,5)= X(I)      
      B(I,  6)= X(I+1)
      B(I+1,6)=-X(I)
      ENDDO
      RNT=REAL(NTA/3)
      DO I=1,6
      DO J=1,6
      S(I,J)=0.0
      ENDDO
      ENDDO
      S(1,1)=RNT
      S(2,2)=RNT
      S(3,3)=RNT
      DO I=1,NTA,3
      S(1,5)=S(1,5)-X(I+2)
      S(1,6)=S(1,6)+X(I+1)
      S(2,4)=S(2,4)+X(I+2)
      S(2,6)=S(2,6)-X(I)
      S(3,4)=S(3,4)-X(I+1)
      S(3,5)=S(3,5)+X(I)
      S(4,5)=S(4,5)-X(I)*X(I+1)
      S(4,6)=S(4,6)-X(I)*X(I+2)
      S(5,6)=S(5,6)-X(I+1)*X(I+2)
      S(4,4)=S(4,4)+X(I+1)**2+X(I+2)**2
      S(5,5)=S(5,5)+X(I)**2  +X(I+2)**2
      S(6,6)=S(6,6)+X(I)**2  +X(I+1)**2
      ENDDO
      S(6,1)=S(1,6)
      S(6,2)=S(2,6)
      S(6,3)=S(3,6)
      S(6,4)=S(4,6)
      S(6,5)=S(5,6)
      S(5,1)=S(1,5)
      S(5,2)=S(2,5)
      S(5,3)=S(3,5)
      S(5,4)=S(4,5)
      S(4,1)=S(1,4)
      S(4,2)=S(2,4)
      S(4,3)=S(3,4)
      S(3,1)=S(1,3)
      S(3,2)=S(2,3)
      S(2,1)=S(1,2)
      CALL MINVER(S,SINV,IND)
      IF(IND.EQ.100)RETURN
      CALL MATP(B,SINV,AUX1,NTA,6,6,NT)
      DO I=1,6
        DO J=1,NTA
          BTR(I,J)=B(J,I)
        ENDDO
      ENDDO
      CALL MATP(AUX1,BTR,P,NTA,6,NTA,NT)
      DO I=1,NTA
      DO J=1,NTA
      P(I,J)=U(I,J)-P(I,J)
      ENDDO
      ENDDO
      RETURN
      END

      SUBROUTINE JACOBIX(AA,N,NP,D,V,NROT)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NMAX=200)
      DIMENSION AA(NP,NP),A(NP,NP),D(NP),V(NP,NP),
     &  B(NMAX),Z(NMAX)
      DO 12 IP=1,N
        DO 11 IQ=1,N
          V(IP,IQ)=0.0
          A(IP,IQ)=AA(IP,IQ)
11      CONTINUE
        V(IP,IP)=1.0
12    CONTINUE
      DO 13 IP=1,N
        B(IP)=A(IP,IP)
        D(IP)=B(IP)
        Z(IP)=0.
13    CONTINUE
      NROT=0
      DO 24 I=1,50
        SM=0.
        DO 15 IP=1,N-1
          DO 14 IQ=IP+1,N
            SM=SM+ABS(A(IP,IQ))
14        CONTINUE
15      CONTINUE
        IF(SM.EQ.0.)RETURN
        IF(I.LT.4)THEN
          TRESH=0.2*SM/N**2
        ELSE
          TRESH=0.
        ENDIF
        DO 22 IP=1,N-1
          DO 21 IQ=IP+1,N
            G=100.*ABS(A(IP,IQ))
            IF((I.GT.4).AND.(ABS(D(IP))+G.EQ.ABS(D(IP)))
     *         .AND.(ABS(D(IQ))+G.EQ.ABS(D(IQ))))THEN
              A(IP,IQ)=0.
            ELSE IF(ABS(A(IP,IQ)).GT.TRESH)THEN
              H=D(IQ)-D(IP)
              IF(ABS(H)+G.EQ.ABS(H))THEN
                T=A(IP,IQ)/H
              ELSE
                THETA=0.5*H/A(IP,IQ)
                T=1./(ABS(THETA)+SQRT(1.+THETA**2))
                IF(THETA.LT.0.)T=-T
              ENDIF
              C=1./SQRT(1+T**2)
              S=T*C
              TAU=S/(1.+C)
              H=T*A(IP,IQ)
              Z(IP)=Z(IP)-H
              Z(IQ)=Z(IQ)+H
              D(IP)=D(IP)-H
              D(IQ)=D(IQ)+H
              A(IP,IQ)=0.
              DO 16 J=1,IP-1
                G=A(J,IP)
                H=A(J,IQ)
                A(J,IP)=G-S*(H+G*TAU)
                A(J,IQ)=H+S*(G-H*TAU)
16            CONTINUE
              DO 17 J=IP+1,IQ-1
                G=A(IP,J)
                H=A(J,IQ)
                A(IP,J)=G-S*(H+G*TAU)
                A(J,IQ)=H+S*(G-H*TAU)
17            CONTINUE
              DO 18 J=IQ+1,N
                G=A(IP,J)
                H=A(IQ,J)
                A(IP,J)=G-S*(H+G*TAU)
                A(IQ,J)=H+S*(G-H*TAU)
18            CONTINUE
              DO 19 J=1,N
                G=V(J,IP)
                H=V(J,IQ)
                V(J,IP)=G-S*(H+G*TAU)
                V(J,IQ)=H+S*(G-H*TAU)
19            CONTINUE
              NROT=NROT+1
            ENDIF
21        CONTINUE
22      CONTINUE
        DO 23 IP=1,N
          B(IP)=B(IP)+Z(IP)
          D(IP)=B(IP)
          Z(IP)=0.
23      CONTINUE
24    CONTINUE
      PAUSE '50 iterations should never happen'
      RETURN
      END


      SUBROUTINE MINVER(A,AINV,IND)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NA=20, NT=3*NA,ND=6)
      DIMENSION RMAT(NT,NT+1),A(NT,NT),AINV(NT,NT),X(NT)
      DO I=1,ND
        DO J=1,ND
          RMAT(I,J)=A(I,J)
        ENDDO
      ENDDO
      DO I=1,ND
      DO J=1,ND
        RMAT(J,ND+1)=0.0
        ENDDO
        RMAT(I,ND+1)=1.0
        CALL SVD(RMAT,ND,X,IND)
        IF(IND.EQ.100)RETURN
        DO J=1,ND
        AINV(J,I)=X(J)
        ENDDO
      ENDDO
      RETURN
      END

      SUBROUTINE SVD(RMAT,N,A,IND)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NA=20, NT=3*NA)
      PARAMETER(NMAX=1000,MMAX=100,TOL=1.0D-9)
      DIMENSION RMAT(NT,N+1),A(NT),V(NT,NT),
     *    U(NT,NT),W(NT),B(NMAX)
      DO 12 I=1,N
        DO 11 J=1,N
          U(I,J)=RMAT(I,J)
11      CONTINUE
        B(I)=RMAT(I,N+1)
12    CONTINUE
      CALL SVDCMP(U,N,N,W,V,IND)
      IF(IND.EQ.100)RETURN
      WMAX=0.0
      DO 13 J=1,N
        IF(W(J).GT.WMAX)WMAX=W(J)
13    CONTINUE
      THRESH=TOL*WMAX
      DO 14 J=1,N
        IF(W(J).LT.THRESH)W(J)=0.0
14    CONTINUE
      CALL SVBKSB(U,W,V,N,N,B,A)
      RETURN
      END

      SUBROUTINE SVDCMP(A,M,N,W,V,IND)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NA=20, NT=3*NA)
      PARAMETER (NMAX=1000)
      DIMENSION A(NT,NT),W(NT),V(NT,NT),RV1(NMAX)
      G=0.0
      SCALE=0.0
      ANORM=0.0
      DO 25 I=1,N
        L=I+1
        RV1(I)=SCALE*G
        G=0.0
        S=0.0
        SCALE=0.0
        IF (I.LE.M) THEN
          DO 11 K=I,M
            SCALE=SCALE+ABS(A(K,I))
11        CONTINUE
          IF (SCALE.NE.0.0) THEN
            DO 12 K=I,M
              A(K,I)=A(K,I)/SCALE
              S=S+A(K,I)*A(K,I)
12          CONTINUE
            F=A(I,I)
            G=-SIGN(SQRT(S),F)
            H=F*G-S
            A(I,I)=F-G
            IF (I.NE.N) THEN
              DO 15 J=L,N
                S=0.0
                DO 13 K=I,M
                  S=S+A(K,I)*A(K,J)
13              CONTINUE
                F=S/H
                DO 14 K=I,M
                  A(K,J)=A(K,J)+F*A(K,I)
14              CONTINUE
15            CONTINUE
            ENDIF
            DO 16 K= I,M
              A(K,I)=SCALE*A(K,I)
16          CONTINUE
          ENDIF
        ENDIF
        W(I)=SCALE *G
        G=0.0
        S=0.0
        SCALE=0.0
        IF ((I.LE.M).AND.(I.NE.N)) THEN
          DO 17 K=L,N
            SCALE=SCALE+ABS(A(I,K))
17        CONTINUE
          IF (SCALE.NE.0.0) THEN
            DO 18 K=L,N
              A(I,K)=A(I,K)/SCALE
              S=S+A(I,K)*A(I,K)
18          CONTINUE
            F=A(I,L)
            G=-SIGN(SQRT(S),F)
            H=F*G-S
            A(I,L)=F-G
            DO 19 K=L,N
              RV1(K)=A(I,K)/H
19          CONTINUE
            IF (I.NE.M) THEN
              DO 23 J=L,M
                S=0.0
                DO 21 K=L,N
                  S=S+A(J,K)*A(I,K)
21              CONTINUE
                DO 22 K=L,N
                  A(J,K)=A(J,K)+S*RV1(K)
22              CONTINUE
23            CONTINUE
            ENDIF
            DO 24 K=L,N
              A(I,K)=SCALE*A(I,K)
24          CONTINUE
          ENDIF
        ENDIF
        ANORM=MAX(ANORM,(ABS(W(I))+ABS(RV1(I))))
25    CONTINUE
      DO 32 I=N,1,-1
        IF (I.LT.N) THEN
          IF (G.NE.0.0) THEN
            DO 26 J=L,N
              V(J,I)=(A(I,J)/A(I,L))/G
26          CONTINUE
            DO 29 J=L,N
              S=0.0
              DO 27 K=L,N
                S=S+A(I,K)*V(K,J)
27            CONTINUE
              DO 28 K=L,N
                V(K,J)=V(K,J)+S*V(K,I)
28            CONTINUE
29          CONTINUE
          ENDIF
          DO 31 J=L,N
            V(I,J)=0.0
            V(J,I)=0.0
31        CONTINUE
        ENDIF
        V(I,I)=1.0
        G=RV1(I)
        L=I
32    CONTINUE
      DO 39 I=N,1,-1
        L=I+1
        G=W(I)
        IF (I.LT.N) THEN
          DO 33 J=L,N
            A(I,J)=0.0
33        CONTINUE
        ENDIF
        IF (G.NE.0.0) THEN
          G=1.0/G
          IF (I.NE.N) THEN
            DO 36 J=L,N
              S=0.0
              DO 34 K=L,M
                S=S+A(K,I)*A(K,J)
34            CONTINUE
              F=(S/A(I,I))*G
              DO 35 K=I,M
                A(K,J)=A(K,J)+F*A(K,I)
35            CONTINUE
36          CONTINUE
          ENDIF
          DO 37 J=I,M
            A(J,I)=A(J,I)*G
37        CONTINUE
        ELSE
          DO 38 J= I,M
            A(J,I)=0.0
38        CONTINUE
        ENDIF
        A(I,I)=A(I,I)+1.0
39    CONTINUE
      DO 49 K=N,1,-1
        DO 48 ITS=1,30
          DO 41 L=K,1,-1
            NM=L-1
            IF ((ABS(RV1(L))+ANORM).EQ.ANORM)  GO TO 2
            IF ((ABS(W(NM))+ANORM).EQ.ANORM)  GO TO 1
41        CONTINUE
1         C=0.0
          S=1.0
          DO 43 I=L,K
            F=S*RV1(I)
            IF ((ABS(F)+ANORM).NE.ANORM) THEN
              G=W(I)
              H=SQRT(F*F+G*G)
              W(I)=H
              H=1.0/H
              C= (G*H)
              S=-(F*H)
              DO 42 J=1,M
                Y=A(J,NM)
                Z=A(J,I)
                A(J,NM)=(Y*C)+(Z*S)
                A(J,I)=-(Y*S)+(Z*C)
42            CONTINUE
            ENDIF
43        CONTINUE
2         Z=W(K)
          IF (L.EQ.K) THEN
            IF (Z.LT.0.0) THEN
              W(K)=-Z
              DO 44 J=1,N
                V(J,K)=-V(J,K)
44            CONTINUE
            ENDIF
            GO TO 3
          ENDIF
          IF (ITS.EQ.30)THEN
          IND=100
          RETURN
          ENDIF
          X=W(L)
          NM=K-1
          Y=W(NM)
          G=RV1(NM)
          H=RV1(K)
          F=((Y-Z)*(Y+Z)+(G-H)*(G+H))/(2.0*H*Y)
          G=SQRT(F*F+1.0)
          F=((X-Z)*(X+Z)+H*((Y/(F+SIGN(G,F)))-H))/X
          C=1.0
          S=1.0
          DO 47 J=L,NM
            I=J+1
            G=RV1(I)
            Y=W(I)
            H=S*G
            G=C*G
            Z=SQRT(F*F+H*H)
            RV1(J)=Z
            C=F/Z
            S=H/Z
            F= (X*C)+(G*S)
            G=-(X*S)+(G*C)
            H=Y*S
            Y=Y*C
            DO 45 NM=1,N
              X=V(NM,J)
              Z=V(NM,I)
              V(NM,J)= (X*C)+(Z*S)
              V(NM,I)=-(X*S)+(Z*C)
45          CONTINUE
            Z=SQRT(F*F+H*H)
            W(J)=Z
            IF (Z.NE.0.0) THEN
              Z=1.0/Z
              C=F*Z
              S=H*Z
            ENDIF
            F= (C*G)+(S*Y)
            X=-(S*G)+(C*Y)
            DO 46 NM=1,M
              Y=A(NM,J)
              Z=A(NM,I)
              A(NM,J)= (Y*C)+(Z*S)
              A(NM,I)=-(Y*S)+(Z*C)
46          CONTINUE
47        CONTINUE
          RV1(L)=0.0
          RV1(K)=F
          W(K)=X
48      CONTINUE
3       CONTINUE
49    CONTINUE
      RETURN
      END


      SUBROUTINE SVBKSB(U,W,V,M,N,B,X)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NMAX=1000)
      PARAMETER (NA=20,NT=3*NA)
      DIMENSION U(NT,NT),W(NT),V(NT,NT),B(NMAX),X(NT),TMP(NMAX)
      DO 12 J=1,N
        S=0.
        IF(W(J).NE.0.)THEN
          DO 11 I=1,M
            S=S+U(I,J)*B(I)
11        CONTINUE
          S=S/W(J)
        ENDIF
        TMP(J)=S
12    CONTINUE
      DO 14 J=1,N
        S=0.
        DO 13 JJ=1,N
          S=S+V(J,JJ)*TMP(JJ)
13      CONTINUE
        X(J)=S
14    CONTINUE
      RETURN
      END

        FUNCTION RANFF ( DUMMY )
        IMPLICIT REAL*8 (A-H,O-Z)
C    *******************************************************************
C    ** RETURNS A UNIFORM RANDOM VARIATE IN THE RANGE 0 TO 1.         **
C    *******************************************************************
        INTEGER     L, C, M
        PARAMETER ( L = 1029, C = 221591, M = 1048576 )
        INTEGER     SEED
        SAVE        SEED
        DATA        SEED / 0 /
C    *******************************************************************
        SEED = MOD ( SEED * L + C, M )
        RANFF = REAL ( SEED ) / M
        RETURN
        END
