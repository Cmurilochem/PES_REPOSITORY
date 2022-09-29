C234567
      function POTEN(R,NRD)
      implicit real*8(A-H,O-Z)
      dimension R(6)
      dimension W(6)
      W(1)=R(1)
      W(2)=R(2)
      W(3)=R(3)
      W(4)=R(6)
      W(5)=R(5)
      W(6)=R(4)
      CALL POTC4(W,POT)
      POTEN=POT
      RETURN
      END

c  !!!****************************************************************!!!
