C--------------------------------------   Fourier Transform
      SUBROUTINE SETF67(IBC1,IQ1,ibc,ip,isl,l1,n2,n3,n4,n7,
     &                  si,indx,Zf,Yf)
c     -----------------------------------------------------
      include 'a_tree.h'
      dimension Zf(narr) , Yf(narr)
      dimension si(nf67) , indx(nf67)

      integer  ibc , ip, isl , l1 , n2 , n3 , n4 , n7
c      INCLUDE 'Control.param'
c      PI=DATAN(1.D0)*4.
       IBC=IBC1
       IQ=IQ1
       IF(IBC.LT.3) GO TO 101
       IQ=IQ-1
  101      N3=2**IQ
       N7=N3/2
       N5=N3/4
       I=1
       INDX(I)=N5
       SI(I)=0.5*SQRT(2.)
       K=I
       I=I+1
  102     IL=I
       IF(I.EQ.N7) GO TO 104
  103     K1=INDX(K)/2
       INDX(I)=K1
       SI(I)=SIN(PI*FLOAT(K1)/FLOAT(N3))
       K1=N7-K1
       I=I+1
       INDX(I)=K1
       SI(I)=SIN(PI*FLOAT(K1)/FLOAT(N3))
       K=K+1
       I=I+1
       IF(K.EQ.IL) GO TO 102
       GO TO 103
  104     RETURN
       END
C-----------------------------------------
       SUBROUTINE TFOLD(IS,L,ZZZ,n2)
       include 'a_tree.h'

       integer  n2

       DIMENSION ZZZ(NARR)
       IH2=N2/2-1
       DO 100 I=IS,IH2
       I1=I+L
       I2=N2-I+L
       A=ZZZ(I1)
       ZZZ(I1)=A-ZZZ(I2)
       ZZZ(I2)=A+ZZZ(I2)
  100     CONTINUE
       RETURN
       END
C------------------------------------
       SUBROUTINE NEG(I1,I3,I2,ibc,ip,isl,l1,n2,n3,n4,n7,si,indx,Zf,Yf)
       include 'a_tree.h'
       dimension Zf(narr) , Yf(narr)
       dimension si(nf67) , indx(nf67)
 
       integer  ibc , ip, isl , l1 , n2 , n3 , n4 , n7

       DO 100 K=I1,I3,I2
       Yf(K)=-Yf(K)
  100     CONTINUE
       RETURN
       END
C-------------------------------------
       SUBROUTINE REVNEG(ibc,ip,isl,l1,n2,n3,n4,n7,si,indx,Zf,Yf)
       include 'a_tree.h'
       dimension Zf(narr) , Yf(narr)
       dimension si(nf67) , indx(nf67)

       integer  ibc , ip, isl , l1 , n2 , n3 , n4 , n7

       DO 100 I=1,N7
       J=N3+1+I
       K=N4+1-I
       A=Yf(J)
       Yf(J)=-Yf(K)
       Yf(K)=-A
  100     CONTINUE
       RETURN
       END
C-------------------------------------
       SUBROUTINE ZEERO(L,ibc,ip,isl,l1,n2,n3,n4,n7,si,indx,Zf,Yf)
       include 'a_tree.h'
       dimension Zf(narr) , Yf(narr)
       dimension si(nf67) , indx(nf67)

       integer  ibc , ip, isl , l1 , n2 , n3 , n4 , n7

       DO 100 I=1,N2
       LI=L+I
       Zf(LI-1)=0.0
  100     CONTINUE
       RETURN
       END
C------------------------------------------
       SUBROUTINE TFOLD1(ibc,ip,isl,l1,n2,n3,n4,n7,si,indx,Zf,Yf)
       include 'a_tree.h'
       dimension Zf(narr) , Yf(narr)
       dimension si(nf67) , indx(nf67)

       integer  ibc , ip, isl , l1 , n2 , n3 , n4 , n7

       II=N2-1
       DO 100 I=1,II
       I1=ISL+I
       I2=L1-I
       A=Zf(I1)
       Zf(I1)=A+Zf(I2)
       Zf(I2)=A-Zf(I2)
  100     CONTINUE
       RETURN
       END
C-------------------------------------
       SUBROUTINE KFOLD(ibc,ip,isl,l1,n2,n3,n4,n7,si,indx,Zf,Yf)
       include 'a_tree.h'
       dimension Zf(narr) , Yf(narr)
       dimension si(nf67) , indx(nf67)

       integer  ibc , ip, isl , l1 , n2 , n3 , n4 , n7

       JS1=N2
       I=1
       J5=ISL+N2
       IS1=ISL
       IC1=L1
       JS1=JS1/2
       IF(JS1.NE.1) GO TO 200
       K1=INDX(I)
       SN=SI(I)
       IS0=IS1
       IS1=IS1+JS1
       IC0=IC1
       IC1=IC1+JS1
       ODD1=SN*(Zf(IC1)-Zf(IS1))
       Yf(K1+1)=Zf(IC0)+ODD1
       N3MK1=N3-K1
       Yf(N3MK1+1)=Zf(IC0)-ODD1
       IF(IBC.LT.3) GO TO 110
       ODD2=SN*(Zf(IC1)+Zf(IS1))
       N3PK1=N3+K1
       Yf(N3PK1+1)=Zf(IS0)+ODD2
       N4MK1=N4-K1
       Yf(N4MK1+1)=-Zf(IS0)+ODD2
  110     RETURN
  200     SN=SI(I)
       IS1=IS1+JS1
       IC1=IC1+JS1
       J3=IS1+JS1
  210     IS0=IS1-JS1
       IC0=IC1-JS1
       ODD1=SN*(Zf(IC1)-Zf(IS1))
       ODD2=SN*(Zf(IC1)+Zf(IS1))
       Zf(IC1)=Zf(IC0)-ODD1
       Zf(IS1)=-Zf(IS0)+ODD2
       Zf(IC0)=Zf(IC0)+ODD1
       Zf(IS0)=Zf(IS0)+ODD2
       IS1=IS1+1
       IC1=IC1+1
       IF(IS1.NE.J3) GO TO 210
       I=I+1
  300     IS1=ISL
       IC1=L1
       JS1=JS1/2
       IF(JS1.EQ.1) GO TO 400
  310     SN=SI(I)
       I=I+1
       CS=SI(I)
       IS1=IS1+JS1
       IC1=IC1+JS1
       J3=IS1+JS1
  320     IS0=IS1-JS1
       IC0=IC1-JS1
       ODD1=CS*Zf(IC1)-SN*Zf(IS1)
       ODD2=SN*Zf(IC1)+CS*Zf(IS1)
       Zf(IC1)=Zf(IC0)-ODD1
       Zf(IC0)=Zf(IC0)+ODD1
       Zf(IS1)=-Zf(IS0)+ODD2
       Zf(IS0)=Zf(IS0)+ODD2
       IS1=IS1+1
       IC1=IC1+1
       IF(IS1.NE.J3) GO TO 320
       IS1=IS1+JS1
       IC1=IC1+JS1
       J3=IS1+JS1
  330     IS0=IS1-JS1
       IC0=IC1-JS1
       ODD1=SN*Zf(IC1)-CS*Zf(IS1)
       ODD2=CS*Zf(IC1)+SN*Zf(IS1)
       Zf(IC1)=Zf(IC0)-ODD1
       Zf(IC0)=Zf(IC0)+ODD1
       Zf(IS1)=-Zf(IS0)+ODD2
       Zf(IS0)=Zf(IS0)+ODD2
       IS1=IS1+1
       IC1=IC1+1
       IF(IS1.NE.J3) GO TO 330
       I=I+1
       IF(IS1.EQ.J5) GO TO 300
       GO TO 310
  400     K1=INDX(I)
       SN=SI(I)
       I=I+1
       CS=SI(I)
       IS0=IS1
       IS1=IS1+JS1
       IC0=IC1
       IC1=IC1+JS1
       ODD1=CS*Zf(IC1)-SN*Zf(IS1)
       Yf(K1+1)=Zf(IC0)+ODD1
       N3MK1=N3-K1
       Yf(N3MK1+1)=Zf(IC0)-ODD1
       IF(IBC.LT.3) GO TO 410
       ODD2=SN*Zf(IC1)+CS*Zf(IS1)
       N3PK1=N3+K1
       Yf(N3PK1+1)=Zf(IS0)+ODD2
       N4MK1=N4-K1
       Yf(N4MK1+1)=-Zf(IS0)+ODD2
  410     IS1=IS1+1
       IC1=IC1+1
       K1=INDX(I)
       IS0=IS1
       IS1=IS1+JS1
       IC0=IC1
       IC1=IC1+JS1
       ODD1=SN*Zf(IC1)-CS*Zf(IS1)
       Yf(K1+1)=Zf(IC0)+ODD1
       N3MK1=N3-K1
       Yf(N3MK1+1)=Zf(IC0)-ODD1
       IF(IBC.LT.3) GO TO 420
       ODD2=CS*Zf(IC1)+SN*Zf(IS1)
       N3PK1=N3+K1
       Yf(N3PK1+1)=Zf(IS0)+ODD2
       N4MK1=N4-K1
       Yf(N4MK1+1)=-Zf(IS0)+ODD2
  420     IS1=IS1+1
       IC1=IC1+1
       I=I+1
       IF(IS1.NE.J5) GO TO 400
       RETURN
       END

c      ------------------------------------------------------
       SUBROUTINE FOUR67(IBC1,IQ1,ip1,isl1,l11,n21,n71,
     &                   si,indx,Zf,Yf)
c      ------------------------------------------------------
       include 'a_tree.h'
       dimension Zf(narr) , Yf(narr)
       dimension si(nf67) , indx(nf67)

       integer  ibc , ip, isl , l1 , n2 , n3 , n4 , n7

c       INCLUDE 'Control.param'
       IBC=IBC1
       IQ=IQ1
       ip = ip1
       isl = isl1
       l1 = l11
       n2 = n21
       n7 = n71
       A5=0.5*SQRT(2.0)
       N4=2**IQ
       N3=N4
       GO TO (103,103,101,102),IBC
  101     Zf(1)=Zf(1)/2.0
       Zf(N3+1)=Zf(1)
       N2=N3

       CALL TFOLD(0,1,Zf,n2)

       N3=N3/2
       IQ=IQ-1
       GO TO 103
  102     N3=N3/2
       IQ=IQ-1
  103     N5=N3/4
       N7=N3/2
       N11=3*N7
       N31=N3+1
       GO TO(300,400,500,600),IBC
  300     Zf(N31)=0.0
       Zf(1)=0.0
       N2=N3
       DO 301 I=2,IQ

       CALL TFOLD(1,1,Zf,n2)

  301     N2=N2/2
       Yf(N7+1)=Zf(2)
       JF=N5
       ISL=1
       DO 302 IP=2,IQ
       L1=N2+1

       CALL ZEERO(1,ibc,ip,isl,l1,n2,n3,n4,n7,si,indx,Zf,Yf)
       CALL KFOLD(ibc,ip,isl,l1,n2,n3,n4,n7,si,indx,Zf,Yf)

       I1=3*JF+1
       I2=4*JF
       I3=I1+(N2/2-1)*I2

       CALL NEG(I1,I3,I2,ibc,ip,isl,l1,n2,n3,n4,n7,si,indx,Zf,Yf)

       N2=N2+N2
  302     JF=JF/2
       RETURN
  400     Zf(1)=Zf(1)/2.0
       Zf(N31)=Zf(N31)/2.0
       N2=N3
       DO 401 I=2,IQ

       CALL TFOLD(0,N31-N2,Zf,n2)

  401     N2=N2/2
       L1=N31-N2
       A=Zf(L1)+Zf(L1+2)
       Yf(1)=-A-Zf(L1+1)
       Yf(N31)=-A+Zf(L1+1)
       Yf(N7+1)=Zf(L1)-Zf(L1+2)
       DO 402 IP=2,IQ
       ISL=N31-N2
       L1=ISL-N2
       idum = isl

       CALL ZEERO(idum,ibc,ip,isl,l1,n2,n3,n4,n7,si,indx,Zf,Yf)
       CALL KFOLD(ibc,ip,isl,l1,n2,n3,n4,n7,si,indx,Zf,Yf)

  402     N2=N2+N2

       CALL NEG(1,N31,2,ibc,ip,isl,l1,n2,n3,n4,n7,si,indx,Zf,Yf)

       RETURN
  500     N2=N3
       L2=N4
       DO 501 IP=2,IQ

       CALL TFOLD(1,1,Zf,n2)
       CALL TFOLD(0,L2-N2+1,Zf,n2)

  501     N2=N2/2
       L1=L2-N2+1
       A=Zf(L1)+Zf(L1+2)
       Yf(N7+1)=2.0*(-Zf(L1)+Zf(L1+2))
       Yf(1)=2.0*(A+Zf(L1+1))
       Yf(N31)=2.0*(A-Zf(L1+1))
       Yf(N11+1)=2.0*Zf(2)
       DO 502 IP=2,IQ
       Zf(N2+1)=2.0*Zf(N2+1)
       ISL=N2+1

       CALL TFOLD1(ibc,ip,isl,l1,n2,n3,n4,n7,si,indx,Zf,Yf)

       L1=L1-N2
       Zf(L1)=-2.0*Zf(L1)

       CALL KFOLD(ibc,ip,isl,l1,n2,n3,n4,n7,si,indx,Zf,Yf)

  502     N2=N2+N2
       Yf(1)=Yf(1)*A5
       Yf(N31)=Yf(N31)*A5
       RETURN
  600     Zf(1)=Zf(1)*A5
       Zf(N31)=Zf(N31)*A5
       N2=N3
       DO 601 IP=2,IQ

       CALL TFOLD(0,N31-N2,Zf,n2)
       CALL TFOLD(1,N31,Zf,n2)

  601     N2=N2/2
       L1=N31-N2
       A=Zf(L1)+Zf(L1+2)
       Yf(1)=2.0*(A+Zf(L1+1))
       Yf(N31)=2.0*(A-Zf(L1+1))
       Yf(N7+1)=2.0*(-Zf(L1)+Zf(L1+2))
       Yf(N11+1)=2.0*Zf(N3+2)
       DO 602 IP=2,IQ
       ISL=N31+N2
       Zf(ISL)=2.0*Zf(ISL)

       CALL TFOLD1(ibc,ip,isl,l1,n2,n3,n4,n7,si,indx,Zf,Yf)

       L1=L1-N2
       Zf(L1)=-2.0*Zf(L1)

       CALL KFOLD(ibc,ip,isl,l1,n2,n3,n4,n7,si,indx,Zf,Yf)

  602     N2=N2+N2

       CALL REVNEG(ibc,ip,isl,l1,n2,n3,n4,n7,si,indx,Zf,Yf)

       idum = n3

       CALL NEG(2,idum,2,ibc,ip,isl,l1,n2,n3,n4,n7,si,indx,Zf,Yf)

       N2=N4

       CALL TFOLD(1,1,Yf,n2)

       Yf(N4+1)=0.0

       RETURN
       END

