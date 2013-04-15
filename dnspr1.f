CA********************************************************************
C                SUBROUTINE DERMAT                   DATE : 20.11.96 *
C--------------------------------------------------------------------*
C                                                                    *
C     BUT : CALCUL DES POINTS DE LA COORDONNEE PHYSIQUE Y            *
C     =====        TABLEAUX :  Y(NY)                                 *
C                                                                    *
C           CALCUL DES MATRICES DERIVEES PAR DIFFERENCES FINIES      *
C            D ORDRE QUATRE MINIMUM SUIVANT X                        *
C                  PREMIERE DX(i,j), SECONDE DDX(i,j)                *
C                  QUATRIEME DDDDX(i,j)                              *
C                                                                    *
C           CALCUL DES MATRICES DERIVEES DE CHEBYSCHEV SUIVANT Y     *
C                  PREMIERE DY(i,j), SECONDE DDY(i,j)                *
C                                                                    *
C     PARAMETRES :                                                   *
C     ============                                                   *
C     Y(i)      (SORTIE) : COORDONNEE PHYSIQUE Y                     *
C     D..(i,j)  (SORTIE) : DERIVEE DE LA MATRICE DANS (X,Y)          *
C     NX,NY     (ENTREE) : NOMBRE DE POINTS DE DISCRETISATION        *
C                          SUIVANT X ET Y                            *
C     YSL       (ENTREE) : POSITIONNEMENT DES POINTS DE COLLOCATION  *
C     PASX      (ENTREE) : PAS D ESPACE SUIVANT LA DIRECTION X       *
CE********************************************************************
C
      SUBROUTINE DERMAT (YSL,PASX,Y,DX,DDX,DDDDX,DY,DDY,NX,NY,YMAX,tr)
C     =============================================================
C
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
C     TRANSFORMATIONS DE VARIABLES ALGEBRIQUE :
C     Y=(YMAX*YSL*(1+KSI))/(2*YSL+YMAX*(1-KSI))
      DIMENSION
     &     PCOLY(NY),Y(NY),DAUX(NY,NY),
     &     DX(7,8),DDX(0:4,5),DDDDX(10,10),DY(NY,NY),DDY(NY,NY)
      dimension tr(ny,2)
C
C     ***************************************************
C     1.) CALCUL DES POINTS DE COLLOCATION DANS (-1,+1) *
C     ***************************************************
C
      CALL CHEBCOL (PCOLY,NY)
C
C     CALCUL DE LA COORDONNEE NORMALE
C
      DO 20 i=1,NY
      Y(NY-i+1)=(YMAX*YSL*(1.D0+PCOLY(i)))/
     &(2.D0*YSL+YMAX*(1.D0-PCOLY(i)))
c      WRITE(6,*) 'Coordonnee normale : ',NY-i+1,Y(NY-i+1)
      write(56,*)i,pcoly(i),y(ny-i+1)
 20   CONTINUE
      close(56)
C
C     *************************************************
C     2.) CALCUL DE LA MATRICE DERIVEE D DANS (-1,+1) *
C     *************************************************
C
      CALL CHEBMAT (DY(1,1),PCOLY(1),NY)
C
C     ****************************************************************
C     3.) CALCUL DES MATRICES DERIVEES DX, DDX, DDDDX, DY, DDY       *
C         DANS (-1,+1) AVEC LES POINTS LIMITES A Y = + 1 AND Y = - 1 *
C     ****************************************************************
C
C CALCUL DES MATRICES DERIVEES PHYSIQUES SUIVANT LA DIRECTION Y
C
      DO 29 i=1,NY
         DO 29 j=1,NY
            DDY(i,j)=0.d0
            DO 29 k=1,NY
               DDY(i,j)=DDY(i,j)+DY(i,k)*DY(k,j)
 29            continue
C
      DO 30 i=1,NY
        FACD=(2.D0*YSL+YMAX*(1.D0-PCOLY(i)))**2/
     &(2.D0*YSL*YMAX*(YSL+YMAX))
        tr(ny+1-i,1)=facd
        FACDD=2.D0*(2.D0*YSL+YMAX*(1.D0-PCOLY(i)))*YMAX/
     &(2.D0*YSL*YMAX*(YSL+YMAX)) 
        tr(ny+1-i,2)=facdd
        DO 30 j=1,NY
          DY(i,j)=DY(i,j)*FACD
          DDY(i,j)=DDY(i,j)*(FACD**2)-FACDD*DY(i,j)
 30   CONTINUE
C
C CHANGEMENT DE LA NOTATION : i = 1 ,PAROI ; i = NY ,INFINI.
C
      DO 40 i=1,NY
      ii=NY-i+1
        DO 40 j=1,NY
        jj=NY-j+1
         DAUX(ii,jj)=DY(i,j)
 40   CONTINUE
      DO 41 i=1,NY
        DO 41 j=1,NY
         DY(i,j)=DAUX(i,j)
 41   CONTINUE
      DO 42 i=1,NY
      ii=NY-i+1
        DO 42 j=1,NY
        jj=NY-j+1
         DAUX(ii,jj)=DDY(i,j)
 42   CONTINUE
      DO 43 i=1,NY
        DO 43 j=1,NY
         DDY(i,j)=DAUX(i,j)
 43   CONTINUE
C
C CALCUL DES MATRICES DERIVEES PHYSIQUES SUIVANT LA DIRECTION X
C
      DO 45 i=1,7
         DO 45 j=1,8
            DX(i,j)=0.0D0
 45   CONTINUE
      DO 46 i=0,4
         DO 46 j=1,5
            DDX(i,j)=0.0D0
 46   CONTINUE
C
         AUX=1.0D0/(840.0D0*PASX)
         DX(5,1)=3.0D0*AUX
         DX(5,2)=-32.0D0*AUX
         DX(5,3)=168.0D0*AUX
         DX(5,4)=-672.0D0*AUX
         DX(5,5)=672.0D0*AUX
         DX(5,6)=-168.0D0*AUX
         DX(5,7)=32.0D0*AUX
         DX(5,8)=-3.0D0*AUX
C
 51      AUX=1.0D0/(60.0D0*PASX)
         DX(4,1)=-1.0D0*AUX
         DX(4,2)=9.0D0*AUX
         DX(4,3)=-45.0D0*AUX
         DX(4,4)=45.0D0*AUX
         DX(4,5)=-9.0D0*AUX
         DX(4,6)=1.0D0*AUX
C
 52      AUX=1.0D0/(12.0D0*PASX)
         DX(3,1)=1.0D0*AUX
         DX(3,2)=-8.0D0*AUX
         DX(3,3)=8.0D0*AUX
         DX(3,4)=-1.0D0*AUX
C
 53      AUX=1.0D0/(12.0D0*PASX)
         DX(2,1)=-3.0D0*AUX
         DX(2,2)=-10.0D0*AUX
         DX(2,3)=18.0D0*AUX
         DX(2,4)=-6.0D0*AUX
         DX(2,5)=1.0D0*AUX
C
 54      AUX=1.0D0/(12.0D0*PASX)
         DX(1,1)=-25.0D0*AUX
         DX(1,2)=48.0D0*AUX
         DX(1,3)=-36.0D0*AUX
         DX(1,4)=16.0D0*AUX
         DX(1,5)=-3.0D0*AUX
C
 55      AUX=1.0D0/(12.0D0*PASX)
         DX(6,5)=3.0D0*AUX
         DX(6,4)=10.0D0*AUX
         DX(6,3)=-18.0D0*AUX
         DX(6,2)=6.0D0*AUX
         DX(6,1)=-1.0D0*AUX
C
 56      AUX=1.0D0/(12.0D0*PASX)
         DX(7,5)=25.0D0*AUX
         DX(7,4)=-48.0D0*AUX
         DX(7,3)=36.0D0*AUX
         DX(7,2)=-16.0D0*AUX
         DX(7,1)=3.0D0*AUX
C
C
         AUX=1.0D0/(12.0D0*(PASX**2))
C
         DDX(0,1)=35.0D0*AUX
         DDX(0,2)=-104.0D0*AUX
         DDX(0,3)=114.0D0*AUX
         DDX(0,4)=-56.0D0*AUX
         DDX(0,5)=11.0D0*AUX
C
         DDX(1,1)=11.0D0*AUX
         DDX(1,2)=-20.0D0*AUX
         DDX(1,3)=6.0D0*AUX
         DDX(1,4)=4.0D0*AUX
         DDX(1,5)=-1.0D0*AUX
C
         DDX(2,1)=-1.0D0*AUX
         DDX(2,2)=16.0D0*AUX
         DDX(2,3)=-30.0D0*AUX
         DDX(2,4)=16.0D0*AUX
         DDX(2,5)=-1.0D0*AUX
C
         DDX(3,5)=11.0D0*AUX
         DDX(3,4)=-20.0D0*AUX
         DDX(3,3)=6.0D0*AUX
         DDX(3,2)=4.0D0*AUX
         DDX(3,1)=-1.0D0*AUX
C
         DDX(4,5)=35.0D0*AUX
         DDX(4,4)=-104.0D0*AUX
         DDX(4,3)=114.0D0*AUX
         DDX(4,2)=-56.0D0*AUX
         DDX(4,1)=11.0D0*AUX
C
C
C
      RETURN
      END
C
CA********************************************************************
C                SUBROUTINE CHEBCOL                   DATE : 20.11.96 
C---------------------------------------------------------------------
C
C     BUT : CALCUL DES POINTS DE COLLOCATION DE CHEBYSHEV DANS (-1,+1)
C     =====        P(i) = COS(PI * i /(IP - 1)),  i = 0,...,IP-1
C
C     PARAMETRES :
C     ============
C
C     P(i)   (SORTIE) : POINTS DE COLLOCATION
C     IP     (ENTREE) : TAILLE DE P ---> = NOMBRE DE POINTS DE
C                                          COLLOCATION IP = NX OU NY
CE********************************************************************
C
      SUBROUTINE CHEBCOL (P, IP)
C     ==========================
C
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DIMENSION P(IP)
C
      NEND  = IP
      N1    = IP-1
      PI    = 2.0D0*DACOS(0.0D0)
C
      DO 100 i = 1,NEND
      P(i) = DCOS(PI*(i-1)/N1)
 100  CONTINUE
C
      RETURN
      END
C
CA********************************************************************
C                SUBROUTINE CHEBMAT                   DATE : 20.11.96
C-------------------------------------------------------------------- 
C
C     BUT : CALCUL DE LA MATRICE DERIVEE PREMIERE CHEBYSHEV D(i,j)
C     =====   DANS (-1,+1) : F(P(i))  = F(i)
C                            F'(P(i)) = F'(i) = D(i,j) * F(j)
C                                     i = 0,...,IP-1;  j = 0,...,IP-1
C
C     PARAMETRES :
C     ============
C
C     D(i,j) (SORTIE) : DERIVEE PREMIERE DE LA MATRICE DANS (-1,+1)
C     P(I)   (ENTREE) : POINTS DE COLLOCATION DE CHEBYSHEV DAN (-1,+1)
C     IP     (ENTREE) : TAILLE DE  P = NOMBRE DE PTS COLL. IP=NX OU NY
CE********************************************************************
C
      SUBROUTINE CHEBMAT (D,  P, IP)
C     ==============================
C
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DIMENSION D(IP,IP),P(IP)
C
      NEND     = IP
      N1       = IP-1
      PI       = 4.0D0 * DATAN(1.0D0)
C
C     CALCUL DE LA MATRICE DERIVEE D(i,j) DANS (-1,+1) :
C     --------------------------------------------------
      DO 500 i = 1,NEND
      DO 500 j = 1,NEND
      I1       = i - 1
      J1       = j - 1
      IF (i.EQ.j) GOTO 100
      CI       = 1.0D0
      IF ((i.EQ.1). OR .(i.EQ.NEND)) CI = 2.0D0
      CJ       = 1.0D0
      IF ((j.EQ.1). OR .(j.EQ.NEND)) CJ = 2.0D0
      D(i,j)   = -(DFLOAT(MOD(I1+J1,2))*2.0D0-1.0D0) * CI/
     & (CJ * (P(i) - P(j)))
      GOTO 500
C     TERMES DIAGONAUX
  100 IF((i.EQ.1). OR .(i.EQ.NEND)) GOTO 200
      PD       = P(i)
      D(i,j)   = -PD/(2.0D0*(1.0D0-PD*PD))
      GOTO 500
  200 IF(i.EQ.1) D(1,1) = (2.0D0*N1*N1+1.D0)/6.0D0
      IF(i.EQ.NEND) D(NEND,NEND) = -(2.0D0*N1*N1+1.0D0)/6.0D0
  500 CONTINUE
C
      RETURN
      END
C
CA********************************************************************
C                 SUBROUTINE BLAS1                   DATE : 20.11.96 *
C--------------------------------------------------------------------*
C                                                                    *
C     BUT : CALCUL DU PROFIL DES VITESSES ADIMENSIONNEES DE BLASIUS  *
C     =====        TABLEAUX :  U(N), V(N) ou N=NCOLY                 *
C                                                                    *
C     PARAMETRES :                                                   *
C     ============                                                   *
C     U(i),V(i) (SORTIE) : VITESSES ADIMENSIONNEES SUIVANT X ET Y    *
C     XA        (ENTREE) : CALCUL DU PROFIL A LA STATION XA          *
C     YCL(I)    (ENTREE) : DETERMINATION DES VITESSES SUIVANT Y      *
C     N         (ENTREE) : NOMBRE DE POINTS DE COLLOCATION SUIVANT Y *
CE********************************************************************
C 
c      SUBROUTINE BLAS1(YCL,um1,U,DU,DDU,DDDU,DDDDU,N)
c
      SUBROUTINE BLAS1(Y,U,V,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(DELT1=1.D-11)
      PARAMETER(DELT2=1.D-11)
      PARAMETER(YMAX=30.0D0)
C
      DIMENSION          U(N),DU(N),DDU(N),DDDU(N),DDDDU(N)
      DIMENSION          Y(N),V(N),um1(N)
      DIMENSION          YCL(N)
      DIMENSION          F(3),DF(3)
      DIMENSION          PARAM(100)
      dimension          iparam(100)
C
      EXTERNAL           BLAS
C
c      DO 10,I  = 1,50
c          param(i)=0.
c   10 CONTINUE
      lrw=100
      liw=100
c      WRITE(6,*) ' CALCUL DU PROFIL DE BLASIUS : BLAS1'
      DO 20 i=1,N
         YCL(N-i+1)=Y(i)*1.7207678D0
c         WRITE(6,*) i,YCL(N-i+1)
 20   CONTINUE
c      PAUSE
C
C     ***********************************************************
C     1.) INITIALISATION DU CALCUL DE F(i)=um1(i) ET DF(i)=U(i) *
C     ***********************************************************
C
      um1  (N) =  0.d0
      U    (N) =  0.D0
      DU   (N) =  0.332057339D0
      DDU  (N) =  0.D0
      DDDU (N) =  0.D0
      DDDDU(N) =  -DU(N)*DU(N)/2.D0
C
      F    (1) =  0.D0
      F    (2) =  U  (N)
      F    (3) =  DU (N)
C
      XANF = YCL(N)
c      WRITE(6,*) ' XANF ',XANF
      IER1 = 1
C
C     *************************************************
C     2.) INTEGRATION AVEC RUNGE-KUTTA (IMSL-ROUTINE) *
C     *************************************************
C
      DO 100 i=1,N-1
      YEND=YCL(N-i)
c      WRITE(6,*) ' YEND = ',YEND
      IF (YEND.LT.YMAX) THEN
        IF(i.EQ.(N-1)) THEN
            IER1=3
        ELSE
            IF(YCL(N-i).GE.YMAX) IER1=3
        END IF
C
c       CALL DIVPRK (IER1, 3, BLAS, XANF, YEND, DELT1, PARAM, F)
        istate=1
        call lsode (blas, 3, f, xanf, yend, 1,delt1, delt2, 1,
     &            istate, 0, param, lrw, iparam, liw, jac, 10)
C
c        write(6,*) ' istate ',istate
        IF (IER1.GT.3) THEN
          write(6,*) ' ---> IER1 > 3'
          STOP
          IER1 = 0
        END IF
C
        um1  (N-i) = F(1)
        VSET     = F(1)
        U    (N-i) = F(2)
        USET     = F(2)
        DU   (N-i) = F(3)
        DDU  (N-i) = -F(1)*F(3)/2.D0
        DDDU (N-i) = (-F(2)*F(3)-F(1)*DDU(N-i))/2.D0
        DDDDU(N-i) = (-F(3)*F(3)-2.D0*F(2)*DDU(N-i)-
     &                F(1)*DDDU(N-i))/2.D0
      ELSE
        U    (N-i) = USET
        DU   (N-i) = 0.D0
        DDU  (N-i) = 0.D0
        DDDU (N-i) = 0.D0
        DDDDU(N-i) = 0.D0
      END IF
 100  CONTINUE
      U(1)=1.D0
      DU(1)=0.D0
      DDU(1)=0.D0
      DDDU(1)=0.D0
      DDDDU(1)=0.D0
C
      DO 200 i=N,1,-1
      IF (YCL(i).gt.YMAX) GOTO 202
      V(i)=YCL(i)*U(i)-um1(i)
      GOTO 200
 202  V(i)=V(i+1)
 200  CONTINUE
C
      DO 300 i=1,N/2
         aux=U(N-i+1)
         U(N-i+1)=U(i)
         U(i)=aux
         aux=DU(N-i+1)
         DU(N-i+1)=DU(i)
         DU(i)=aux
         aux=DDU(N-i+1)
         DDU(N-i+1)=DDU(i)
         DDU(i)=aux
         aux=V(N-i+1)
         V(N-i+1)=V(i)
         V(i)=aux
 300     CONTINUE
      RETURN
      END
C
C
      SUBROUTINE BLAS( M , X , Y , DY )
C
      IMPLICIT REAL * 8 (A-H,O-Z)
      DIMENSION         Y(3), DY(3)
C
C     ********************************************
C     BERECHNUNG DY GEMAESS BLASIUS-GLEICHUNG
C     ********************************************
C
C     DIMENSIONSLOS, WIE IM BLASIUS-FALL:
C
      DY (1) = Y(2)
      DY (2) = Y(3)
      DY (3) = -Y(1)*Y(3)/2.D0
C
      RETURN
      END
C
C *** RESOLUTION DE [A]*[B]=[C] ***
C
      SUBROUTINE RESOL(NI,NJ,A,B,C,IR,IC)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(NI,NI),B(NI,NJ),C(NI,NJ)
      DIMENSION W1(NI),W2(NI)
      DIMENSION IR(NI),IC(NI)
c
      CALL LUPIV(A,NI,NI,IR,IC,IPER,W2)
c
      DO 1 J=1,NJ
         DO 2 I=1,NI
            W2(I)=C(I,J)
   2     CONTINUE
         CALL SYSTEM1(NI,NI,A,W1,W2,IR,IC)
         DO 3 I=1,NI
C        WRITE(6,*) ' W1',W1(I)
            B(I,J)=W1(I)
   3     CONTINUE
   1  CONTINUE
c
      RETURN
      END
C
C *** RESOLUTION DE [A]*[B]=[C] ***
C
      SUBROUTINE RESOLQ(NI,NJ,A,B,C,IR,IC)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(NI,NI),B(NI,NJ),C(NI,NJ)
      DIMENSION W1(NI),W2(NI)
      DIMENSION IR(NI),IC(NI)
C
      DO 1 J=1,NJ
         DO 2 I=1,NI
            W2(I)=C(I,J)
   2     CONTINUE
         CALL SYSTEM1(NI,NI,A,W1,W2,IR,IC)
         DO 3 I=1,NI
C        WRITE(6,*) ' W1',W1(I)
            B(I,J)=W1(I)
   3     CONTINUE
   1  CONTINUE
c
      RETURN
      END
C                                                                       00252700
CA********************************************************************* 00243300
C                SUBROUTINE LUPIV3                     DATE : 18.12.96  00243400
C---------------------------------------------------------------------- 00243500
C                                                                       00243600
C     PURPOSE : LU DECOMPOSITION OF A GENERAL REAL MATRIX A OF ORDER N  00243700
C     =======   USING NOT ROW AND COLUMN PIVOTING
C               A MATRICE DE TYPE BANDE
C                                                                       00243900
C     COMMENTS : PROGRAMMED BY U. EHRENSTEIN                            00244000
C     ========== NEEDS FUNCTION IAMAX AND SUBROUTINE RSCALU             00244100
C                                                                       00244200
C     CHANGES :                                                         00244300
C     =========                                                         00244400
C                                                                       00244500
C     PARAMETERS :                                                      00244600
C     ============                                                      00244700
C                                                                       00244800
C     A      : INPUT : N BY N REAL MATRIX                               00244900
C              OUTPUT: LU DECOMPOSITION OF A                            00245000
C     N      : NUMBER OF ROWS OF MATRIX A                               00245100
C     IR     : INTEGER VECTOR CONTAINING THE ROW PERMUTATIONS (N-DIM)   00245200
C     IC     : INTEGER VECTOR CONTAINING THE COLUMN PERMUTATIONS (N-DIM)00245300
C     IPER   : TOTAL NUMBER OF ALL PERMUTATIONS                         00245400
C     B      : WORKING VECTOR OF LENGTH N                               00245500
CE********************************************************************* 00245600
C                                                                       00245700
      SUBROUTINE LUPIV3 (A, N, NDIM, IR, IC, IPER, B)
C     ===============================================
      IMPLICIT REAL*8(A-H,O-Y), COMPLEX*16(Z)                           00246000
C                                                                       00246100
      DIMENSION A(*),B(*)                                               00246200
      DIMENSION IR(*),IC(*)
C
      IPER=0                                                            00246400
C                                                                       00246500
      DO 1 I=1,N                                                        00246600
      IR(I)=I                                                           00246700
      IC(I)=I                                                           00246800
    1 CONTINUE                                                          00246900
C     ********************************************                      00247000
      DO 10 K=1,N-1                                                     00247200
      C=1./A((K-1)*NDIM+K)                                              00251300
      K1=K+1                                                            00251400
C     RSCALU                                                            00251500
C     PERFORMS Y=A*Y WHERE A REAL NUMBER AND Y REAL VECTOR              00251600
      CALL RSCALU(N-K,C,A((K-1)*NDIM+K1),1)
      DO 5 J=K1,N
         IF ((J.ge.K1).and.((J-K1).le.3)) THEN
            KI=K1
            NI=J+2
            ENDIF
         IF (((J-K1).gt.3).and.(J.lt.(N-3))) THEN
            KI=K1+J-5
            NI=J+2
            ENDIF
         IF (J.ge.(N-3)) THEN
            KI=K1
            NI=N
            ENDIF
      DO 6 I=KI,NI
      A((J-1)*NDIM+I)=A((J-1)*NDIM+I)-A((K-1)*NDIM+I)*A((J-1)*NDIM+K)
    6 CONTINUE                                                          00252100
    5 CONTINUE                                                          00252200
   10 CONTINUE                                                          00252400
      RETURN                                                            00252500
      END                                                               00252600
C                                                                       00252700
CA********************************************************************* 00243300
C                SUBROUTINE LUPIV2                     DATE : 18.12.96  00243400
C---------------------------------------------------------------------- 00243500
C                                                                       00243600
C     PURPOSE : LU DECOMPOSITION OF A GENERAL REAL MATRIX A OF ORDER N  00243700
C     =======   USING NOT ROW AND COLUMN PIVOTING
C               A MATRICE DE TYPE BANDE
C                                                                       00243900
C     COMMENTS : PROGRAMMED BY U. EHRENSTEIN                            00244000
C     ========== NEEDS FUNCTION IAMAX AND SUBROUTINE RSCALU             00244100
C                                                                       00244200
C     CHANGES :                                                         00244300
C     =========                                                         00244400
C                                                                       00244500
C     PARAMETERS :                                                      00244600
C     ============                                                      00244700
C                                                                       00244800
C     A      : INPUT : N BY N REAL MATRIX                               00244900
C              OUTPUT: LU DECOMPOSITION OF A                            00245000
C     N      : NUMBER OF ROWS OF MATRIX A                               00245100
C     IR     : INTEGER VECTOR CONTAINING THE ROW PERMUTATIONS (N-DIM)   00245200
C     IC     : INTEGER VECTOR CONTAINING THE COLUMN PERMUTATIONS (N-DIM)00245300
C     IPER   : TOTAL NUMBER OF ALL PERMUTATIONS                         00245400
C     B      : WORKING VECTOR OF LENGTH N                               00245500
CE********************************************************************* 00245600
C                                                                       00245700
      SUBROUTINE LUPIV2 (A, N, NDIM, IR, IC, IPER, B)
C     ===============================================
      IMPLICIT REAL*8(A-H,O-Y), COMPLEX*16(Z)                           00246000
C                                                                       00246100
      DIMENSION A(*),B(*)                                               00246200
      DIMENSION IR(*),IC(*)
C
      IPER=0                                                            00246400
C                                                                       00246500
      DO 1 I=1,N                                                        00246600
      IR(I)=I                                                           00246700
      IC(I)=I                                                           00246800
    1 CONTINUE                                                          00246900
C     ********************************************                      00247000
      DO 10 K=1,N-1                                                     00247200
      C=1./A((K-1)*NDIM+K)                                              00251300
      K1=K+1                                                            00251400
C     RSCALU                                                            00251500
C     PERFORMS Y=A*Y WHERE A REAL NUMBER AND Y REAL VECTOR              00251600
      CALL RSCALU(N-K,C,A((K-1)*NDIM+K1),1)
      DO 5 J=K1,N
         IF ((J.ge.K1).and.((J-K1).le.4)) THEN
            KI=K1
            NI=J+3
            ENDIF
         IF (((J-K1).gt.4).and.(J.lt.(N-3))) THEN
            KI=K1+J-6
            NI=J+3
            ENDIF
         IF (J.ge.(N-3)) THEN
            KI=K1
            NI=N
            ENDIF
      DO 6 I=KI,NI
      A((J-1)*NDIM+I)=A((J-1)*NDIM+I)-A((K-1)*NDIM+I)*A((J-1)*NDIM+K)
    6 CONTINUE                                                          00252100
    5 CONTINUE                                                          00252200
   10 CONTINUE                                                          00252400
      RETURN                                                            00252500
      END                                                               00252600
C                                                                       00252700
CA********************************************************************* 00243300
C                SUBROUTINE LUPIV                      DATE : 23.11.87  00243400
C---------------------------------------------------------------------- 00243500
C                                                                       00243600
C     PURPOSE : LU DECOMPOSITION OF A GENERAL REAL MATRIX A OF ORDER N  00243700
C     =======   USING ROW AND COLUMN PIVOTING                           00243800
C                                                                       00243900
C     COMMENTS : PROGRAMMED BY U. EHRENSTEIN                            00244000
C     ========== NEEDS FUNCTION IAMAX AND SUBROUTINE RSCALU             00244100
C                                                                       00244200
C     CHANGES :                                                         00244300
C     =========                                                         00244400
C                                                                       00244500
C     PARAMETERS :                                                      00244600
C     ============                                                      00244700
C                                                                       00244800
C     A      : INPUT : N BY N REAL MATRIX                               00244900
C              OUTPUT: LU DECOMPOSITION OF A                            00245000
C     N      : NUMBER OF ROWS OF MATRIX A                               00245100
C     IR     : INTEGER VECTOR CONTAINING THE ROW PERMUTATIONS (N-DIM)   00245200
C     IC     : INTEGER VECTOR CONTAINING THE COLUMN PERMUTATIONS (N-DIM)00245300
C     IPER   : TOTAL NUMBER OF ALL PERMUTATIONS                         00245400
C     B      : WORKING VECTOR OF LENGTH N                               00245500
CE********************************************************************* 00245600
C                                                                       00245700
      SUBROUTINE LUPIV (A, N, NDIM, IR, IC, IPER, B)                    00245800
C     ==============================================                    00245900
      IMPLICIT REAL*8(A-H,O-Y), COMPLEX*16(Z)                           00246000
C                                                                       00246100
      DIMENSION A(*),B(*)                                               00246200
      DIMENSION IR(*),IC(*)                                             00246300
      IPER=0                                                            00246400
C                                                                       00246500
      DO 1 I=1,N                                                        00246600
      IR(I)=I                                                           00246700
      IC(I)=I                                                           00246800
    1 CONTINUE                                                          00246900
C     ********************************************                      00247000
      XMAX=DABS(A(1))                                                   00247100
C
      DO 10 K=1,N-1                                                     00247200
      L=K                                                               00247300
      M=K                                                               00247400
          DO 2 J=K,N                                                    00247500
C    IAMAX                                                              00247600
C    DETERMINES FIRST INDEX I SUCH THAT                                 00247700
C    DABS(XI)=MAX DABS(XJ)  J=1,N                                       00247800
          IMAX=IAMAX(N-K+1,A((J-1)*NDIM+K),1)                           00247900
          YMAX=DABS(A((J-1)*NDIM+K+IMAX-1))                             00248000
          IF (YMAX.GT.XMAX) THEN                                        00248100
          XMAX=YMAX                                                     00248200
          L=IMAX+K-1                                                    00248300
          M=J                                                           00248400
          END IF                                                        00248500
    2     CONTINUE                                                      00248600
C
C
C     ********************************************                      00248700
      IRL=IR(L)                                                         00248800
      IR(L)=IR(K)                                                       00248900
      IR(K)=IRL                                                         00249000
      ICM=IC(M)                                                         00249100
      IC(M)=IC(K)                                                       00249200
      IC(K)=ICM                                                         00249300
C     *********************************************                     00249400
      IF (L.NE.K) THEN                                                  00249500
      IPER=IPER+1                                                       00249600
      DO 3 J=1,N                                                        00249700
      B(J)=A((J-1)*NDIM+K)                                              00249800
      A((J-1)*NDIM+K)=A((J-1)*NDIM+L)                                   00249900
      A((J-1)*NDIM+L)=B(J)                                              00250000
    3 CONTINUE                                                          00250100
      END IF                                                            00250200
      IF (M.NE.K) THEN                                                  00250300
      IPER=IPER+1                                                       00250400
      DO 4 I=1,N                                                        00250500
      B(I)=A((K-1)*NDIM+I)                                              00250600
      A((K-1)*NDIM+I)=A((M-1)*NDIM+I)                                   00250700
      A((M-1)*NDIM+I)=B(I)                                              00250800
    4 CONTINUE                                                          00250900
      END IF                                                            00251000
C     ************************************                              00251100
C                                                                       00251200
      C=1./A((K-1)*NDIM+K)                                              00251300
      K1=K+1                                                            00251400
C     RSCALU                                                            00251500
C     PERFORMS Y=A*Y WHERE A REAL NUMBER AND Y REAL VECTOR              00251600
      CALL RSCALU(N-K,C,A((K-1)*NDIM+K1),1)                             00251700
      DO 5 J=K1,N                                                       00251800
      DO 6 I=K1,N                                                       00251900
      A((J-1)*NDIM+I)=A((J-1)*NDIM+I)-A((K-1)*NDIM+I)*A((J-1)*NDIM+K)   00252000
    6 CONTINUE                                                          00252100
    5 CONTINUE                                                          00252200
      XMAX=DABS(A(K*NDIM+K1))                                           00252300
   10 CONTINUE                                                          00252400
C
      RETURN                                                            00252500
      END                                                               00252600
C                                                                       00252700
CA********************************************************************* 00252800
C                FUNCTION IAMAX                        DATE : 23.11.87  00252900
C---------------------------------------------------------------------- 00253000
C                                                                       00253100
C     PURPOSE : THIS FUNCTION DETERMINES THE FIRST INDEX I SUCH THAT    00253200
C     =======   DABS(CX(I)) IS A MAXIMUM;                               00253300
C               CX IS A REAL ARRAY; SKIP ELEMENT =INCX                  00253400
C                                                                       00253500
C     COMMENTS : PROGRAMMED BY U. EHRENSTEIN                            00253600
C     ==========                                                        00253700
C                                                                       00253800
CE********************************************************************* 00253900
C                                                                       00254000
      FUNCTION IAMAX(N, CX, INCX)                                       00254100
C     ===========================                                       00254200
      IMPLICIT REAL*8(A-H,O-Y), COMPLEX*16(Z)                           00254300
C                                                                       00254400
      DIMENSION CX(*)                                                   00254500
C                                                                       00254600
      EM=DABS(CX(1))                                                    00254700
      DO 1 I=1,N                                                        00254800
      EM1=DABS(CX(1+(I-1)*INCX))                                        00254900
      IF(EM1.GE.EM) THEN                                                00255000
      EM=EM1                                                            00255100
      IM=I                                                              00255200
      END IF                                                            00255300
    1 CONTINUE                                                          00255400
      IAMAX=IM                                                          00255500
      RETURN                                                            00255600
      END                                                               00255700
C                                                                       00255800
CA********************************************************************* 00255900
C                SUBROUTINE RSCALU                     DATE : 23.11.87  00256000
C---------------------------------------------------------------------- 00256100
C                                                                       00256200
C     PURPOSE : PERFORMS THE PRODUCT OF A REAL VECTOR C TIMES A REAL    00256300
C     =======   NUMBER C ---> C * A                                     00256400
C                                                                       00256500
C     COMMENTS : PROGRAMMED BY U. EHRENSTEIN                            00256600
C     ==========                                                        00256700
C                                                                       00256800
CE********************************************************************* 00256900
C                                                                       00257000
      SUBROUTINE RSCALU (N, C, A, INCX)                                 00257100
C     =================================                                 00257200
      IMPLICIT REAL*8(A-H,O-Y), COMPLEX*16(Z)                           00257300
C                                                                       00257400
      DIMENSION A(*)                                                    00257500
C                                                                       00257600
      N1=N*INCX                                                         00257700
C                                                                       00257800
      DO 1 I=1,N1,INCX                                                  00257900
      A(I)=A(I)*C                                                       00258000
    1 CONTINUE                                                          00258100
      RETURN                                                            00258200
      END                                                               00258300
C                                                                       00258400
C                                                                       00252700
CA********************************************************************* 00273100
C                FUNCTION SDOTU                        DATE : 23.11.87  00273200
C---------------------------------------------------------------------- 00273300
C                                                                       00273400
C     PURPOSE : THIS FUNCTION COMPUTES THE SCALAR PRODUCT OF TWO        00273500
C     =======   REAL ARRAYS SX AND SY ; SKIP ELEMENT FOR SX =INCX       00273600
C               SKIP ELEMENT FOR SY = INCY                              00273700
C                                                                       00273800
C     COMMENTS : PROGRAMMED BY U. EHRENSTEIN                            00273900
C     ==========                                                        00274000
C                                                                       00274100
CE********************************************************************* 00274200
C                                                                       00274300
      FUNCTION SDOTU(N, SX, INCX, SY, INCY)                             00274400
C     =====================================                             00274500
      IMPLICIT REAL*8(A-H,O-Y), COMPLEX*16(Z)                           00274600
C                                                                       00274700
      DIMENSION SX(*),SY(*)                                             00274800
C                                                                       00274900
      DOT=0.D0                                                          00275000
      DO 1 I=1,N                                                        00275100
      DOT=DOT+SX(1+I*INCX-INCX)*SY(1+I*INCY-INCY)                       00275200
    1 CONTINUE                                                          00275300
      SDOTU=DOT                                                         00275400
      RETURN                                                            00275500
      END                                                               00275600
C                                                                       00275700
CA********************************************************************* 00270600
C                SUBROUTINE PRVEVE                     DATE : 23.11.87  00270700
C---------------------------------------------------------------------- 00270800
C                                                                       00270900
C     PURPOSE : THIS ROUTINE PERFORMS THE SUM CX(I)*CY(I) ,I=1,N        00271000
C     =======   WHERE CX AND CY ARE REAL VECTORS OF LENGTH N            00271100
C               IN THE MEAN PROGRAM ; SKIP ELEMENT FOR CX =IX           00271200
C               SKIP ELEMENT FOR CY = IY ;                              00271300
C                                                                       00271400
C     COMMENTS : PROGRAMMED BY U. EHRENSTEIN                            00271500
C     ==========                                                        00271600
C                                                                       00271700
CE********************************************************************* 00271800
C                                                                       00271900
      SUBROUTINE PRVEVE(N, CX, IX, CY, IY, CDOT)                        00272000
C     ==========================================                        00272100
      IMPLICIT REAL*8(A-H,O-Y), COMPLEX*16(Z)                           00272200
C                                                                       00272300
      DIMENSION CX(*),CY(*)                                             00272400
C                                                                       00272500
      DOT1=SDOTU(N,CX(1),IX,CY(1),IY)                                   00272600
      CDOT=DOT1                                                         00272700
      RETURN                                                            00272800
      END                                                               00272900
C                                                                       00273000
CA********************************************************************* 00258500
C                SUBROUTINE SYSTEM                     DATE : 23.11.87  00258600
C---------------------------------------------------------------------- 00258700
C                                                                       00258800
C     PURPOSE : SOLUTION OF A GENERAL LINEAR REAL SYSTEM OF EQUATIONS   00258900
C     =======                    A * X = B                              00259000
C               PREVIOUS CALL OF PROGRAM LUPIV NECESSARY                00259100
C               ROW AND COLUMN PERMUTATIONS OF LU DECOMPOSITION GIVEN BY00259200
C               THE VECTORS IR AND IC                                   00259300
C               A CONTAINS ALREADY THE LU DECOMPOSITION                 00259400
C                                                                       00259500
C     COMMENTS : PROGRAMMED BY U. EHRENSTEIN                            00259600
C     ========== NEEDS SUBROUTINE PRVEVE AND FUNCTION SDOTU             00259700
C                                                                       00259800
C     CHANGES :                                                         00259900
C     =========                                                         00260000
C                                                                       00260100
C     PARAMETERS :                                                      00260200
C     ============                                                      00260300
C                                                                       00260400
C     A      : LU DECOMPOSITION OF N BY N MATRIX A                      00260500
C     N      : NUMBER OF ROWS OF MATRIX A                               00260600
C     IR     : INTEGER VECTOR CONTAINING THE ROW PERMUTATIONS (N-DIM)   00260700
C     IC     : INTEGER VECTOR CONTAINING THE COLUMN PERMUTATIONS (N-DIM)00260800
C     B      : RHS INPUT VECTOR (N - DIM)                               00260900
CE********************************************************************* 00261000
C                                                                       00261100
      SUBROUTINE SYSTEM1 (N, NDIM, A, X, B, IR, IC)                      00261200
C     ============================================                      00261300
      IMPLICIT REAL*8(A-H,O-Y), COMPLEX*16(Z)                           00261400
C                                                                       00261500
      DIMENSION A(*),X(*),B(*)                                          00261600
      INTEGER IR(*),IC(*)                                               00261700
C                                                                       00261800
      DO 1 I=1,N                                                        00261900
      X(I)=B(IR(I))                                                     00262000
    1 CONTINUE                                                          00262100
      B(1)=X(1)                                                         00262200
      DO 2 I=2,N                                                        00262300
c
      cdot=SDOTU(i-1,a(i),ndim,b(1),1)
c
c      CALL PRVEVE(I-1,A(I),NDIM,B(1),1,CDOT)                            00262400
      B(I)=X(I)-CDOT                                                    00262500
    2 CONTINUE                                                          00262600
      X(N)=1./A((N-1)*NDIM+N)*B(N)                                      00262700
      DO 3 I=1,N-1                                                      00262800
c
      cdot=SDOTU(i,a((n-i)*ndim+n-i),ndim,x(n-i+1),1)
c
c      CALL PRVEVE(I,A((N-I)*NDIM+N-I),NDIM,X(N-I+1),1,CDOT)             00262900
      X(N-I)=1./A((N-I-1)*NDIM+N-I)*(B(N-I)-CDOT)                       00263000
    3 CONTINUE                                                          00263100
      DO 4 I=1,N                                                        00263200
      B(IC(I))=X(I)                                                     00263300
    4 CONTINUE                                                          00263400
      DO 5 I=1,N                                                        00263500
      X(I)=B(I)                                                         00263600
    5 CONTINUE                                                          00263700
      RETURN                                                            00263800
      END                                                               00263900
C                                                                       00264000
C     
      subroutine ene(nmode,letha,lops,lvec,k,dy,bv,rx,zmode,
     &               zdmode,alpha,rsum)
      implicit real * 8 (a-h,o-z)
      complex * 16 zmode,zdmode
      dimension dy(k,k),bv(k),rx(lvec)
      dimension zmode(-nmode:nmode),zdmode(-nmode:nmode)
      rsum=0.d0
      do 1 l=2,k
         zmode(0)=dcmplx(rx(letha+lops+l),0.d0)
         dpsr=0.d0
         do 2 j=1,k
            dpsr=dpsr+dy(l,j)*rx(letha+lops+j)
    2    continue
         zdmode(0)=dcmplx(dpsr,0.d0)
         do 3 n=1,nmode
            zmode(n)=dcmplx(rx(letha+lops+k+2*(n-1)*k+2*l-1),
     &                      rx(letha+lops+k+2*(n-1)*k+2*l))
            dpsr=0.d0
            dpsi=0.d0
            do 4 j=1,k
               dpsr=dpsr+dy(l,j)*rx(letha+lops+k+2*(n-1)*k+2*j-1)
               dpsi=dpsi+dy(l,j)*rx(letha+lops+k+2*(n-1)*k+2*j)
    4       continue
            zdmode(n)=dcmplx(dpsr,dpsi)
            rsum=rsum+dfloat(n)**2*alpha**2*bv(l)*
     &           (dreal(zmode(n))**2+dimag(zmode(n))**2)
     &           +bv(l)*(dreal(zdmode(n))**2+dimag(zdmode(n))**2)
    3    continue     
    1 continue
      rsum=rsum*2.d0*rx(lvec-1)**2
      return
      end 
C
      subroutine lsode (f, neq, y, t, tout, itol, rtol, atol, itask,            
     1            istate, iopt, rwork, lrw, iwork, liw, jac, mf)                
      external f, jac                                                           
      integer neq, itol, itask, istate, iopt, lrw, iwork, liw, mf               
      double precision y, t, tout, rtol, atol, rwork                            
      dimension neq(1), y(1), rtol(1), atol(1), rwork(lrw), iwork(liw)          
c-----------------------------------------------------------------------        
c this is the march 30, 1987 version of                                         
c lsode.. livermore solver for ordinary differential equations.                 
c this version is in double precision.                                          
c                                                                               
c lsode solves the initial value problem for stiff or nonstiff                  
c systems of first order ode-s,                                                 
c     dy/dt = f(t,y) ,  or, in component form,                                  
c     dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(neq)) (i = 1,...,neq).            
c lsode is a package based on the gear and gearb packages, and on the           
c october 23, 1978 version of the tentative odepack user interface              
c standard, with minor modifications.                                           
c-----------------------------------------------------------------------        
c reference..                                                                   
c     alan c. hindmarsh,  odepack, a systematized collection of ode             
c     solvers, in scientific computing, r. s. stepleman et al. (eds.),          
c     north-holland, amsterdam, 1983, pp. 55-64.                                
c-----------------------------------------------------------------------        
c author and contact.. alan c. hindmarsh,                                       
c                      computing and mathematics research div., l-316           
c                      lawrence livermore national laboratory                   
c                      livermore, ca 94550.                                     
c-----------------------------------------------------------------------        
c summary of usage.                                                             
c                                                                               
c communication between the user and the lsode package, for normal              
c situations, is summarized here.  this summary describes only a subset         
c of the full set of options available.  see the full description for           
c details, including optional communication, nonstandard options,               
c and instructions for special situations.  see also the example                
c problem (with program and output) following this summary.                     
c                                                                               
c a. first provide a subroutine of the form..                                   
c               subroutine f (neq, t, y, ydot)                                  
c               dimension y(neq), ydot(neq)                                     
c which supplies the vector function f by loading ydot(i) with f(i).            
c                                                                               
c b. next determine (or guess) whether or not the problem is stiff.             
c stiffness occurs when the jacobian matrix df/dy has an eigenvalue             
c whose real part is negative and large in magnitude, compared to the           
c reciprocal of the t span of interest.  if the problem is nonstiff,            
c use a method flag mf = 10.  if it is stiff, there are four standard           
c choices for mf, and lsode requires the jacobian matrix in some form.          
c this matrix is regarded either as full (mf = 21 or 22),                       
c or banded (mf = 24 or 25).  in the banded case, lsode requires two            
c half-bandwidth parameters ml and mu.  these are, respectively, the            
c widths of the lower and upper parts of the band, excluding the main           
c diagonal.  thus the band consists of the locations (i,j) with                 
c i-ml .le. j .le. i+mu, and the full bandwidth is ml+mu+1.                     
c                                                                               
c c. if the problem is stiff, you are encouraged to supply the jacobian         
c directly (mf = 21 or 24), but if this is not feasible, lsode will             
c compute it internally by difference quotients (mf = 22 or 25).                
c if you are supplying the jacobian, provide a subroutine of the form..         
c               subroutine jac (neq, t, y, ml, mu, pd, nrowpd)                  
c               dimension y(neq), pd(nrowpd,neq)                                
c which supplies df/dy by loading pd as follows..                               
c     for a full jacobian (mf = 21), load pd(i,j) with df(i)/dy(j),             
c the partial derivative of f(i) with respect to y(j).  (ignore the             
c ml and mu arguments in this case.)                                            
c     for a banded jacobian (mf = 24), load pd(i-j+mu+1,j) with                 
c df(i)/dy(j), i.e. load the diagonal lines of df/dy into the rows of           
c pd from the top down.                                                         
c     in either case, only nonzero elements need be loaded.                     
c                                                                               
c d. write a main program which calls subroutine lsode once for                 
c each point at which answers are desired.  this should also provide            
c for possible use of logical unit 6 for output of error messages               
c by lsode.  on the first call to lsode, supply arguments as follows..          
c f      = name of subroutine for right-hand side vector f.                     
c          this name must be declared external in calling program.              
c neq    = number of first order ode-s.                                         
c y      = array of initial values, of length neq.                              
c t      = the initial value of the independent variable.                       
c tout   = first point where output is desired (.ne. t).                        
c itol   = 1 or 2 according as atol (below) is a scalar or array.               
c rtol   = relative tolerance parameter (scalar).                               
c atol   = absolute tolerance parameter (scalar or array).                      
c          the estimated local error in y(i) will be controlled so as           
c          to be roughly less (in magnitude) than                               
c             ewt(i) = rtol*abs(y(i)) + atol     if itol = 1, or                
c             ewt(i) = rtol*abs(y(i)) + atol(i)  if itol = 2.                   
c          thus the local error test passes if, in each component,              
c          either the absolute error is less than atol (or atol(i)),            
c          or the relative error is less than rtol.                             
c          use rtol = 0.0 for pure absolute error control, and                  
c          use atol = 0.0 (or atol(i) = 0.0) for pure relative error            
c          control.  caution.. actual (global) errors may exceed these          
c          local tolerances, so choose them conservatively.                     
c itask  = 1 for normal computation of output values of y at t = tout.          
c istate = integer flag (input and output).  set istate = 1.                    
c iopt   = 0 to indicate no optional inputs used.                               
c rwork  = real work array of length at least..                                 
c             20 + 16*neq                    for mf = 10,                       
c             22 +  9*neq + neq**2           for mf = 21 or 22,                 
c             22 + 10*neq + (2*ml + mu)*neq  for mf = 24 or 25.                 
c lrw    = declared length of rwork (in user-s dimension).                      
c iwork  = integer work array of length at least..                              
c             20        for mf = 10,                                            
c             20 + neq  for mf = 21, 22, 24, or 25.                             
c          if mf = 24 or 25, input in iwork(1),iwork(2) the lower               
c          and upper half-bandwidths ml,mu.                                     
c liw    = declared length of iwork (in user-s dimension).                      
c jac    = name of subroutine for jacobian matrix (mf = 21 or 24).              
c          if used, this name must be declared external in calling              
c          program.  if not used, pass a dummy name.                            
c mf     = method flag.  standard values are..                                  
c          10 for nonstiff (adams) method, no jacobian used.                    
c          21 for stiff (bdf) method, user-supplied full jacobian.              
c          22 for stiff method, internally generated full jacobian.             
c          24 for stiff method, user-supplied banded jacobian.                  
c          25 for stiff method, internally generated banded jacobian.           
c note that the main program must declare arrays y, rwork, iwork,               
c and possibly atol.                                                            
c                                                                               
c e. the output from the first call (or any call) is..                          
c      y = array of computed values of y(t) vector.                             
c      t = corresponding value of independent variable (normally tout).         
c istate = 2  if lsode was successful, negative otherwise.                      
c          -1 means excess work done on this call (perhaps wrong mf).           
c          -2 means excess accuracy requested (tolerances too small).           
c          -3 means illegal input detected (see printed message).               
c          -4 means repeated error test failures (check all inputs).            
c          -5 means repeated convergence failures (perhaps bad jacobian         
c             supplied or wrong choice of mf or tolerances).                    
c          -6 means error weight became zero during problem. (solution          
c             component i vanished, and atol or atol(i) = 0.)                   
c                                                                               
c f. to continue the integration after a successful return, simply              
c reset tout and call lsode again.  no other parameters need be reset.          
c                                                                               
c-----------------------------------------------------------------------        
c example problem.                                                              
c                                                                               
c the following is a simple example problem, with the coding                    
c needed for its solution by lsode.  the problem is from chemical               
c kinetics, and consists of the following three rate equations..                
c     dy1/dt = -.04*y1 + 1.e4*y2*y3                                             
c     dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e7*y2**2                                 
c     dy3/dt = 3.e7*y2**2                                                       
c on the interval from t = 0.0 to t = 4.e10, with initial conditions            
c y1 = 1.0, y2 = y3 = 0.  the problem is stiff.                                 
c                                                                               
c the following coding solves this problem with lsode, using mf = 21            
c and printing results at t = .4, 4., ..., 4.e10.  it uses                      
c itol = 2 and atol much smaller for y2 than y1 or y3 because                   
c y2 has much smaller values.                                                   
c at the end of the run, statistical quantities of interest are                 
c printed (see optional outputs in the full description below).                 
c                                                                               
c     external fex, jex                                                         
c     double precision atol, rtol, rwork, t, tout, y                            
c     dimension y(3), atol(3), rwork(58), iwork(23)                             
c     neq = 3                                                                   
c     y(1) = 1.d0                                                               
c     y(2) = 0.d0                                                               
c     y(3) = 0.d0                                                               
c     t = 0.d0                                                                  
c     tout = .4d0                                                               
c     itol = 2                                                                  
c     rtol = 1.d-4                                                              
c     atol(1) = 1.d-6                                                           
c     atol(2) = 1.d-10                                                          
c     atol(3) = 1.d-6                                                           
c     itask = 1                                                                 
c     istate = 1                                                                
c     iopt = 0                                                                  
c     lrw = 58                                                                  
c     liw = 23                                                                  
c     mf = 21                                                                   
c     do 40 iout = 1,12                                                         
c       call lsode(fex,neq,y,t,tout,itol,rtol,atol,itask,istate,                
c    1     iopt,rwork,lrw,iwork,liw,jex,mf)                                     
c       write(6,20)t,y(1),y(2),y(3)                                             
c 20    format(7h at t =,e12.4,6h   y =,3e14.6)                                 
c       if (istate .lt. 0) go to 80                                             
c 40    tout = tout*10.d0                                                       
c     write(6,60)iwork(11),iwork(12),iwork(13)                                  
c 60  format(/12h no. steps =,i4,11h  no. f-s =,i4,11h  no. j-s =,i4)           
c     stop                                                                      
c 80  write(6,90)istate                                                         
c 90  format(///22h error halt.. istate =,i3)                                   
c     stop                                                                      
c     end                                                                       
c                                                                               
c     subroutine fex (neq, t, y, ydot)                                          
c     double precision t, y, ydot                                               
c     dimension y(3), ydot(3)                                                   
c     ydot(1) = -.04d0*y(1) + 1.d4*y(2)*y(3)                                    
c     ydot(3) = 3.d7*y(2)*y(2)                                                  
c     ydot(2) = -ydot(1) - ydot(3)                                              
c     return                                                                    
c     end                                                                       
c                                                                               
c     subroutine jex (neq, t, y, ml, mu, pd, nrpd)                              
c     double precision pd, t, y                                                 
c     dimension y(3), pd(nrpd,3)                                                
c     pd(1,1) = -.04d0                                                          
c     pd(1,2) = 1.d4*y(3)                                                       
c     pd(1,3) = 1.d4*y(2)                                                       
c     pd(2,1) = .04d0                                                           
c     pd(2,3) = -pd(1,3)                                                        
c     pd(3,2) = 6.d7*y(2)                                                       
c     pd(2,2) = -pd(1,2) - pd(3,2)                                              
c     return                                                                    
c     end                                                                       
c                                                                               
c the output of this program (on a cdc-7600 in single precision)                
c is as follows..                                                               
c                                                                               
c   at t =  4.0000e-01   y =  9.851726e-01  3.386406e-05  1.479357e-02          
c   at t =  4.0000e+00   y =  9.055142e-01  2.240418e-05  9.446344e-02          
c   at t =  4.0000e+01   y =  7.158050e-01  9.184616e-06  2.841858e-01          
c   at t =  4.0000e+02   y =  4.504846e-01  3.222434e-06  5.495122e-01          
c   at t =  4.0000e+03   y =  1.831701e-01  8.940379e-07  8.168290e-01          
c   at t =  4.0000e+04   y =  3.897016e-02  1.621193e-07  9.610297e-01          
c   at t =  4.0000e+05   y =  4.935213e-03  1.983756e-08  9.950648e-01          
c   at t =  4.0000e+06   y =  5.159269e-04  2.064759e-09  9.994841e-01          
c   at t =  4.0000e+07   y =  5.306413e-05  2.122677e-10  9.999469e-01          
c   at t =  4.0000e+08   y =  5.494529e-06  2.197824e-11  9.999945e-01          
c   at t =  4.0000e+09   y =  5.129458e-07  2.051784e-12  9.999995e-01          
c   at t =  4.0000e+10   y = -7.170586e-08 -2.868234e-13  1.000000e+00          
c                                                                               
c   no. steps = 330  no. f-s = 405  no. j-s =  69                               
c-----------------------------------------------------------------------        
c full description of user interface to lsode.                                  
c                                                                               
c the user interface to lsode consists of the following parts.                  
c                                                                               
c i.   the call sequence to subroutine lsode, which is a driver                 
c      routine for the solver.  this includes descriptions of both              
c      the call sequence arguments and of user-supplied routines.               
c      following these descriptions is a description of                         
c      optional inputs available through the call sequence, and then            
c      a description of optional outputs (in the work arrays).                  
c                                                                               
c ii.  descriptions of other routines in the lsode package that may be          
c      (optionally) called by the user.  these provide the ability to           
c      alter error message handling, save and restore the internal              
c      common, and obtain specified derivatives of the solution y(t).           
c                                                                               
c iii. descriptions of common blocks to be declared in overlay                  
c      or similar environments, or to be saved when doing an interrupt          
c      of the problem and continued solution later.                             
c                                                                               
c iv.  description of two routines in the lsode package, either of              
c      which the user may replace with his own version, if desired.             
c      these relate to the measurement of errors.                               
c                                                                               
c-----------------------------------------------------------------------        
c part i.  call sequence.                                                       
c                                                                               
c the call sequence parameters used for input only are                          
c     f, neq, tout, itol, rtol, atol, itask, iopt, lrw, liw, jac, mf,           
c and those used for both input and output are                                  
c     y, t, istate.                                                             
c the work arrays rwork and iwork are also used for conditional and             
c optional inputs and optional outputs.  (the term output here refers           
c to the return from subroutine lsode to the user-s calling program.)           
c                                                                               
c the legality of input parameters will be thoroughly checked on the            
c initial call for the problem, but not checked thereafter unless a             
c change in input parameters is flagged by istate = 3 on input.                 
c                                                                               
c the descriptions of the call arguments are as follows.                        
c                                                                               
c f      = the name of the user-supplied subroutine defining the                
c          ode system.  the system must be put in the first-order               
c          form dy/dt = f(t,y), where f is a vector-valued function             
c          of the scalar t and the vector y.  subroutine f is to                
c          compute the function f.  it is to have the form                      
c               subroutine f (neq, t, y, ydot)                                  
c               dimension y(1), ydot(1)                                         
c          where neq, t, and y are input, and the array ydot = f(t,y)           
c          is output.  y and ydot are arrays of length neq.                     
c          (in the dimension statement above, 1 is a dummy                      
c          dimension.. it can be replaced by any value.)                        
c          subroutine f should not alter y(1),...,y(neq).                       
c          f must be declared external in the calling program.                  
c                                                                               
c          subroutine f may access user-defined quantities in                   
c          neq(2),... and/or in y(neq(1)+1),... if neq is an array              
c          (dimensioned in f) and/or y has length exceeding neq(1).             
c          see the descriptions of neq and y below.                             
c                                                                               
c          if quantities computed in the f routine are needed                   
c          externally to lsode, an extra call to f should be made               
c          for this purpose, for consistent and accurate results.               
c          if only the derivative dy/dt is needed, use intdy instead.           
c                                                                               
c neq    = the size of the ode system (number of first order                    
c          ordinary differential equations).  used only for input.              
c          neq may be decreased, but not increased, during the problem.         
c          if neq is decreased (with istate = 3 on input), the                  
c          remaining components of y should be left undisturbed, if             
c          these are to be accessed in f and/or jac.                            
c                                                                               
c          normally, neq is a scalar, and it is generally referred to           
c          as a scalar in this user interface description.  however,            
c          neq may be an array, with neq(1) set to the system size.             
c          (the lsode package accesses only neq(1).)  in either case,           
c          this parameter is passed as the neq argument in all calls            
c          to f and jac.  hence, if it is an array, locations                   
c          neq(2),... may be used to store other integer data and pass          
c          it to f and/or jac.  subroutines f and/or jac must include           
c          neq in a dimension statement in that case.                           
c                                                                               
c y      = a real array for the vector of dependent variables, of               
c          length neq or more.  used for both input and output on the           
c          first call (istate = 1), and only for output on other calls.         
c          on the first call, y must contain the vector of initial              
c          values.  on output, y contains the computed solution vector,         
c          evaluated at t.  if desired, the y array may be used                 
c          for other purposes between calls to the solver.                      
c                                                                               
c          this array is passed as the y argument in all calls to               
c          f and jac.  hence its length may exceed neq, and locations           
c          y(neq+1),... may be used to store other real data and                
c          pass it to f and/or jac.  (the lsode package accesses only           
c          y(1),...,y(neq).)                                                    
c                                                                               
c t      = the independent variable.  on input, t is used only on the           
c          first call, as the initial point of the integration.                 
c          on output, after each call, t is the value at which a                
c          computed solution y is evaluated (usually the same as tout).         
c          on an error return, t is the farthest point reached.                 
c                                                                               
c tout   = the next value of t at which a computed solution is desired.         
c          used only for input.                                                 
c                                                                               
c          when starting the problem (istate = 1), tout may be equal            
c          to t for one call, then should .ne. t for the next call.             
c          for the initial t, an input value of tout .ne. t is used             
c          in order to determine the direction of the integration               
c          (i.e. the algebraic sign of the step sizes) and the rough            
c          scale of the problem.  integration in either direction               
c          (forward or backward in t) is permitted.                             
c                                                                               
c          if itask = 2 or 5 (one-step modes), tout is ignored after            
c          the first call (i.e. the first call with tout .ne. t).               
c          otherwise, tout is required on every call.                           
c                                                                               
c          if itask = 1, 3, or 4, the values of tout need not be                
c          monotone, but a value of tout which backs up is limited              
c          to the current internal t interval, whose endpoints are              
c          tcur - hu and tcur (see optional outputs, below, for                 
c          tcur and hu).                                                        
c                                                                               
c itol   = an indicator for the type of error control.  see                     
c          description below under atol.  used only for input.                  
c                                                                               
c rtol   = a relative error tolerance parameter, either a scalar or             
c          an array of length neq.  see description below under atol.           
c          input only.                                                          
c                                                                               
c atol   = an absolute error tolerance parameter, either a scalar or            
c          an array of length neq.  input only.                                 
c                                                                               
c             the input parameters itol, rtol, and atol determine               
c          the error control performed by the solver.  the solver will          
c          control the vector e = (e(i)) of estimated local errors              
c          in y, according to an inequality of the form                         
c                      rms-norm of ( e(i)/ewt(i) )   .le.   1,                  
c          where       ewt(i) = rtol(i)*abs(y(i)) + atol(i),                    
c          and the rms-norm (root-mean-square norm) here is                     
c          rms-norm(v) = sqrt(sum v(i)**2 / neq).  here ewt = (ewt(i))          
c          is a vector of weights which must always be positive, and            
c          the values of rtol and atol should all be non-negative.              
c          the following table gives the types (scalar/array) of                
c          rtol and atol, and the corresponding form of ewt(i).                 
c                                                                               
c             itol    rtol       atol          ewt(i)                           
c              1     scalar     scalar     rtol*abs(y(i)) + atol                
c              2     scalar     array      rtol*abs(y(i)) + atol(i)             
c              3     array      scalar     rtol(i)*abs(y(i)) + atol             
c              4     array      array      rtol(i)*abs(y(i)) + atol(i)          
c                                                                               
c          when either of these parameters is a scalar, it need not             
c          be dimensioned in the user-s calling program.                        
c                                                                               
c          if none of the above choices (with itol, rtol, and atol              
c          fixed throughout the problem) is suitable, more general              
c          error controls can be obtained by substituting                       
c          user-supplied routines for the setting of ewt and/or for             
c          the norm calculation.  see part iv below.                            
c                                                                               
c          if global errors are to be estimated by making a repeated            
c          run on the same problem with smaller tolerances, then all            
c          components of rtol and atol (i.e. of ewt) should be scaled           
c          down uniformly.                                                      
c                                                                               
c itask  = an index specifying the task to be performed.                        
c          input only.  itask has the following values and meanings.            
c          1  means normal computation of output values of y(t) at              
c             t = tout (by overshooting and interpolating).                     
c          2  means take one step only and return.                              
c          3  means stop at the first internal mesh point at or                 
c             beyond t = tout and return.                                       
c          4  means normal computation of output values of y(t) at              
c             t = tout but without overshooting t = tcrit.                      
c             tcrit must be input as rwork(1).  tcrit may be equal to           
c             or beyond tout, but not behind it in the direction of             
c             integration.  this option is useful if the problem                
c             has a singularity at or beyond t = tcrit.                         
c          5  means take one step, without passing tcrit, and return.           
c             tcrit must be input as rwork(1).                                  
c                                                                               
c          note..  if itask = 4 or 5 and the solver reaches tcrit               
c          (within roundoff), it will return t = tcrit (exactly) to             
c          indicate this (unless itask = 4 and tout comes before tcrit,         
c          in which case answers at t = tout are returned first).               
c                                                                               
c istate = an index used for input and output to specify the                    
c          the state of the calculation.                                        
c                                                                               
c          on input, the values of istate are as follows.                       
c          1  means this is the first call for the problem                      
c             (initializations will be done).  see note below.                  
c          2  means this is not the first call, and the calculation             
c             is to continue normally, with no change in any input              
c             parameters except possibly tout and itask.                        
c             (if itol, rtol, and/or atol are changed between calls             
c             with istate = 2, the new values will be used but not              
c             tested for legality.)                                             
c          3  means this is not the first call, and the                         
c             calculation is to continue normally, but with                     
c             a change in input parameters other than                           
c             tout and itask.  changes are allowed in                           
c             neq, itol, rtol, atol, iopt, lrw, liw, mf, ml, mu,                
c             and any of the optional inputs except h0.                         
c             (see iwork description for ml and mu.)                            
c          note..  a preliminary call with tout = t is not counted              
c          as a first call here, as no initialization or checking of            
c          input is done.  (such a call is sometimes useful for the             
c          purpose of outputting the initial conditions.)                       
c          thus the first call for which tout .ne. t requires                   
c          istate = 1 on input.                                                 
c                                                                               
c          on output, istate has the following values and meanings.             
c           1  means nothing was done, as tout was equal to t with              
c              istate = 1 on input.  (however, an internal counter was          
c              set to detect and prevent repeated calls of this type.)          
c           2  means the integration was performed successfully.                
c          -1  means an excessive amount of work (more than mxstep              
c              steps) was done on this call, before completing the              
c              requested task, but the integration was otherwise                
c              successful as far as t.  (mxstep is an optional input            
c              and is normally 500.)  to continue, the user may                 
c              simply reset istate to a value .gt. 1 and call again             
c              (the excess work step counter will be reset to 0).               
c              in addition, the user may increase mxstep to avoid               
c              this error return (see below on optional inputs).                
c          -2  means too much accuracy was requested for the precision          
c              of the machine being used.  this was detected before             
c              completing the requested task, but the integration               
c              was successful as far as t.  to continue, the tolerance          
c              parameters must be reset, and istate must be set                 
c              to 3.  the optional output tolsf may be used for this            
c              purpose.  (note.. if this condition is detected before           
c              taking any steps, then an illegal input return                   
c              (istate = -3) occurs instead.)                                   
c          -3  means illegal input was detected, before taking any              
c              integration steps.  see written message for details.             
c              note..  if the solver detects an infinite loop of calls          
c              to the solver with illegal input, it will cause                  
c              the run to stop.                                                 
c          -4  means there were repeated error test failures on                 
c              one attempted step, before completing the requested              
c              task, but the integration was successful as far as t.            
c              the problem may have a singularity, or the input                 
c              may be inappropriate.                                            
c          -5  means there were repeated convergence test failures on           
c              one attempted step, before completing the requested              
c              task, but the integration was successful as far as t.            
c              this may be caused by an inaccurate jacobian matrix,             
c              if one is being used.                                            
c          -6  means ewt(i) became zero for some i during the                   
c              integration.  pure relative error control (atol(i)=0.0)          
c              was requested on a variable which has now vanished.              
c              the integration was successful as far as t.                      
c                                                                               
c          note..  since the normal output value of istate is 2,                
c          it does not need to be reset for normal continuation.                
c          also, since a negative input value of istate will be                 
c          regarded as illegal, a negative output value requires the            
c          user to change it, and possibly other inputs, before                 
c          calling the solver again.                                            
c                                                                               
c iopt   = an integer flag to specify whether or not any optional               
c          inputs are being used on this call.  input only.                     
c          the optional inputs are listed separately below.                     
c          iopt = 0 means no optional inputs are being used.                    
c                   default values will be used in all cases.                   
c          iopt = 1 means one or more optional inputs are being used.           
c                                                                               
c rwork  = a real working array (double precision).                             
c          the length of rwork must be at least                                 
c             20 + nyh*(maxord + 1) + 3*neq + lwm    where                      
c          nyh    = the initial value of neq,                                   
c          maxord = 12 (if meth = 1) or 5 (if meth = 2) (unless a               
c                   smaller value is given as an optional input),               
c          lwm   = 0             if miter = 0,                                  
c          lwm   = neq**2 + 2    if miter is 1 or 2,                            
c          lwm   = neq + 2       if miter = 3, and                              
c          lwm   = (2*ml+mu+1)*neq + 2 if miter is 4 or 5.                      
c          (see the mf description for meth and miter.)                         
c          thus if maxord has its default value and neq is constant,            
c          this length is..                                                     
c             20 + 16*neq                  for mf = 10,                         
c             22 + 16*neq + neq**2         for mf = 11 or 12,                   
c             22 + 17*neq                  for mf = 13,                         
c             22 + 17*neq + (2*ml+mu)*neq  for mf = 14 or 15,                   
c             20 +  9*neq                  for mf = 20,                         
c             22 +  9*neq + neq**2         for mf = 21 or 22,                   
c             22 + 10*neq                  for mf = 23,                         
c             22 + 10*neq + (2*ml+mu)*neq  for mf = 24 or 25.                   
c          the first 20 words of rwork are reserved for conditional             
c          and optional inputs and optional outputs.                            
c                                                                               
c          the following word in rwork is a conditional input..                 
c            rwork(1) = tcrit = critical value of t which the solver            
c                       is not to overshoot.  required if itask is              
c                       4 or 5, and ignored otherwise.  (see itask.)            
c                                                                               
c lrw    = the length of the array rwork, as declared by the user.              
c          (this will be checked by the solver.)                                
c                                                                               
c iwork  = an integer work array.  the length of iwork must be at least         
c             20        if miter = 0 or 3 (mf = 10, 13, 20, 23), or             
c             20 + neq  otherwise (mf = 11, 12, 14, 15, 21, 22, 24, 25).        
c          the first few words of iwork are used for conditional and            
c          optional inputs and optional outputs.                                
c                                                                               
c          the following 2 words in iwork are conditional inputs..              
c            iwork(1) = ml     these are the lower and upper                    
c            iwork(2) = mu     half-bandwidths, respectively, of the            
c                       banded jacobian, excluding the main diagonal.           
c                       the band is defined by the matrix locations             
c                       (i,j) with i-ml .le. j .le. i+mu.  ml and mu            
c                       must satisfy  0 .le.  ml,mu  .le. neq-1.                
c                       these are required if miter is 4 or 5, and              
c                       ignored otherwise.  ml and mu may in fact be            
c                       the band parameters for a matrix to which               
c                       df/dy is only approximately equal.                      
c                                                                               
c liw    = the length of the array iwork, as declared by the user.              
c          (this will be checked by the solver.)                                
c                                                                               
c note..  the work arrays must not be altered between calls to lsode            
c for the same problem, except possibly for the conditional and                 
c optional inputs, and except for the last 3*neq words of rwork.                
c the latter space is used for internal scratch space, and so is                
c available for use by the user outside lsode between calls, if                 
c desired (but not for use by f or jac).                                        
c                                                                               
c jac    = the name of the user-supplied routine (miter = 1 or 4) to            
c          compute the jacobian matrix, df/dy, as a function of                 
c          the scalar t and the vector y.  it is to have the form               
c               subroutine jac (neq, t, y, ml, mu, pd, nrowpd)                  
c               dimension y(1), pd(nrowpd,1)                                    
c          where neq, t, y, ml, mu, and nrowpd are input and the array          
c          pd is to be loaded with partial derivatives (elements of             
c          the jacobian matrix) on output.  pd must be given a first            
c          dimension of nrowpd.  t and y have the same meaning as in            
c          subroutine f.  (in the dimension statement above, 1 is a             
c          dummy dimension.. it can be replaced by any value.)                  
c               in the full matrix case (miter = 1), ml and mu are              
c          ignored, and the jacobian is to be loaded into pd in                 
c          columnwise manner, with df(i)/dy(j) loaded into pd(i,j).             
c               in the band matrix case (miter = 4), the elements               
c          within the band are to be loaded into pd in columnwise               
c          manner, with diagonal lines of df/dy loaded into the rows            
c          of pd.  thus df(i)/dy(j) is to be loaded into pd(i-j+mu+1,j).        
c          ml and mu are the half-bandwidth parameters (see iwork).             
c          the locations in pd in the two triangular areas which                
c          correspond to nonexistent matrix elements can be ignored             
c          or loaded arbitrarily, as they are overwritten by lsode.             
c               jac need not provide df/dy exactly.  a crude                    
c          approximation (possibly with a smaller bandwidth) will do.           
c               in either case, pd is preset to zero by the solver,             
c          so that only the nonzero elements need be loaded by jac.             
c          each call to jac is preceded by a call to f with the same            
c          arguments neq, t, and y.  thus to gain some efficiency,              
c          intermediate quantities shared by both calculations may be           
c          saved in a user common block by f and not recomputed by jac,         
c          if desired.  also, jac may alter the y array, if desired.            
c          jac must be declared external in the calling program.                
c               subroutine jac may access user-defined quantities in            
c          neq(2),... and/or in y(neq(1)+1),... if neq is an array              
c          (dimensioned in jac) and/or y has length exceeding neq(1).           
c          see the descriptions of neq and y above.                             
c                                                                               
c mf     = the method flag.  used only for input.  the legal values of          
c          mf are 10, 11, 12, 13, 14, 15, 20, 21, 22, 23, 24, and 25.           
c          mf has decimal digits meth and miter.. mf = 10*meth + miter.         
c          meth indicates the basic linear multistep method..                   
c            meth = 1 means the implicit adams method.                          
c            meth = 2 means the method based on backward                        
c                     differentiation formulas (bdf-s).                         
c          miter indicates the corrector iteration method..                     
c            miter = 0 means functional iteration (no jacobian matrix           
c                      is involved).                                            
c            miter = 1 means chord iteration with a user-supplied               
c                      full (neq by neq) jacobian.                              
c            miter = 2 means chord iteration with an internally                 
c                      generated (difference quotient) full jacobian            
c                      (using neq extra calls to f per df/dy value).            
c            miter = 3 means chord iteration with an internally                 
c                      generated diagonal jacobian approximation.               
c                      (using 1 extra call to f per df/dy evaluation).          
c            miter = 4 means chord iteration with a user-supplied               
c                      banded jacobian.                                         
c            miter = 5 means chord iteration with an internally                 
c                      generated banded jacobian (using ml+mu+1 extra           
c                      calls to f per df/dy evaluation).                        
c          if miter = 1 or 4, the user must supply a subroutine jac             
c          (the name is arbitrary) as described above under jac.                
c          for other values of miter, a dummy argument can be used.             
c-----------------------------------------------------------------------        
c optional inputs.                                                              
c                                                                               
c the following is a list of the optional inputs provided for in the            
c call sequence.  (see also part ii.)  for each such input variable,            
c this table lists its name as used in this documentation, its                  
c location in the call sequence, its meaning, and the default value.            
c the use of any of these inputs requires iopt = 1, and in that                 
c case all of these inputs are examined.  a value of zero for any               
c of these optional inputs will cause the default value to be used.             
c thus to use a subset of the optional inputs, simply preload                   
c locations 5 to 10 in rwork and iwork to 0.0 and 0 respectively, and           
c then set those of interest to nonzero values.                                 
c                                                                               
c name    location      meaning and default value                               
c                                                                               
c h0      rwork(5)  the step size to be attempted on the first step.            
c                   the default value is determined by the solver.              
c                                                                               
c hmax    rwork(6)  the maximum absolute step size allowed.                     
c                   the default value is infinite.                              
c                                                                               
c hmin    rwork(7)  the minimum absolute step size allowed.                     
c                   the default value is 0.  (this lower bound is not           
c                   enforced on the final step before reaching tcrit            
c                   when itask = 4 or 5.)                                       
c                                                                               
c maxord  iwork(5)  the maximum order to be allowed.  the default               
c                   value is 12 if meth = 1, and 5 if meth = 2.                 
c                   if maxord exceeds the default value, it will                
c                   be reduced to the default value.                            
c                   if maxord is changed during the problem, it may             
c                   cause the current order to be reduced.                      
c                                                                               
c mxstep  iwork(6)  maximum number of (internally defined) steps                
c                   allowed during one call to the solver.                      
c                   the default value is 500.                                   
c                                                                               
c mxhnil  iwork(7)  maximum number of messages printed (per problem)            
c                   warning that t + h = t on a step (h = step size).           
c                   this must be positive to result in a non-default            
c                   value.  the default value is 10.                            
c-----------------------------------------------------------------------        
c optional outputs.                                                             
c                                                                               
c as optional additional output from lsode, the variables listed                
c below are quantities related to the performance of lsode                      
c which are available to the user.  these are communicated by way of            
c the work arrays, but also have internal mnemonic names as shown.              
c except where stated otherwise, all of these outputs are defined               
c on any successful return from lsode, and on any return with                   
c istate = -1, -2, -4, -5, or -6.  on an illegal input return                   
c (istate = -3), they will be unchanged from their existing values              
c (if any), except possibly for tolsf, lenrw, and leniw.                        
c on any error return, outputs relevant to the error will be defined,           
c as noted below.                                                               
c                                                                               
c name    location      meaning                                                 
c                                                                               
c hu      rwork(11) the step size in t last used (successfully).                
c                                                                               
c hcur    rwork(12) the step size to be attempted on the next step.             
c                                                                               
c tcur    rwork(13) the current value of the independent variable               
c                   which the solver has actually reached, i.e. the             
c                   current internal mesh point in t.  on output, tcur          
c                   will always be at least as far as the argument              
c                   t, but may be farther (if interpolation was done).          
c                                                                               
c tolsf   rwork(14) a tolerance scale factor, greater than 1.0,                 
c                   computed when a request for too much accuracy was           
c                   detected (istate = -3 if detected at the start of           
c                   the problem, istate = -2 otherwise).  if itol is            
c                   left unaltered but rtol and atol are uniformly              
c                   scaled up by a factor of tolsf for the next call,           
c                   then the solver is deemed likely to succeed.                
c                   (the user may also ignore tolsf and alter the               
c                   tolerance parameters in any other way appropriate.)         
c                                                                               
c nst     iwork(11) the number of steps taken for the problem so far.           
c                                                                               
c nfe     iwork(12) the number of f evaluations for the problem so far.         
c                                                                               
c nje     iwork(13) the number of jacobian evaluations (and of matrix           
c                   lu decompositions) for the problem so far.                  
c                                                                               
c nqu     iwork(14) the method order last used (successfully).                  
c                                                                               
c nqcur   iwork(15) the order to be attempted on the next step.                 
c                                                                               
c imxer   iwork(16) the index of the component of largest magnitude in          
c                   the weighted local error vector ( e(i)/ewt(i) ),            
c                   on an error return with istate = -4 or -5.                  
c                                                                               
c lenrw   iwork(17) the length of rwork actually required.                      
c                   this is defined on normal returns and on an illegal         
c                   input return for insufficient storage.                      
c                                                                               
c leniw   iwork(18) the length of iwork actually required.                      
c                   this is defined on normal returns and on an illegal         
c                   input return for insufficient storage.                      
c                                                                               
c the following two arrays are segments of the rwork array which                
c may also be of interest to the user as optional outputs.                      
c for each array, the table below gives its internal name,                      
c its base address in rwork, and its description.                               
c                                                                               
c name    base address      description                                         
c                                                                               
c yh      21             the nordsieck history array, of size nyh by            
c                        (nqcur + 1), where nyh is the initial value            
c                        of neq.  for j = 0,1,...,nqcur, column j+1             
c                        of yh contains hcur**j/factorial(j) times              
c                        the j-th derivative of the interpolating               
c                        polynomial currently representing the solution,        
c                        evaluated at t = tcur.                                 
c                                                                               
c acor     lenrw-neq+1   array of size neq used for the accumulated             
c                        corrections on each step, scaled on output             
c                        to represent the estimated local error in y            
c                        on the last step.  this is the vector e in             
c                        the description of the error control.  it is           
c                        defined only on a successful return from lsode.        
c                                                                               
c-----------------------------------------------------------------------        
c part ii.  other routines callable.                                            
c                                                                               
c the following are optional calls which the user may make to                   
c gain additional capabilities in conjunction with lsode.                       
c (the routines xsetun and xsetf are designed to conform to the                 
c slatec error handling package.)                                               
c                                                                               
c     form of call                  function                                    
c   call xsetun(lun)          set the logical unit number, lun, for             
c                             output of messages from lsode, if                 
c                             the default is not desired.                       
c                             the default value of lun is 6.                    
c                                                                               
c   call xsetf(mflag)         set a flag to control the printing of             
c                             messages by lsode.                                
c                             mflag = 0 means do not print. (danger..           
c                             this risks losing valuable information.)          
c                             mflag = 1 means print (the default).              
c                                                                               
c                             either of the above calls may be made at          
c                             any time and will take effect immediately.        
c                                                                               
c   call srcom(rsav,isav,job) saves and restores the contents of                
c                             the internal common blocks used by                
c                             lsode (see part iii below).                       
c                             rsav must be a real array of length 218           
c                             or more, and isav must be an integer              
c                             array of length 41 or more.                       
c                             job=1 means save common into rsav/isav.           
c                             job=2 means restore common from rsav/isav.        
c                                srcom is useful if one is                      
c                             interrupting a run and restarting                 
c                             later, or alternating between two or              
c                             more problems solved with lsode.                  
c                                                                               
c   call intdy(,,,,,)         provide derivatives of y, of various              
c        (see below)          orders, at a specified point t, if                
c                             desired.  it may be called only after             
c                             a successful return from lsode.                   
c                                                                               
c the detailed instructions for using intdy are as follows.                     
c the form of the call is..                                                     
c                                                                               
c   call intdy (t, k, rwork(21), nyh, dky, iflag)                               
c                                                                               
c the input parameters are..                                                    
c                                                                               
c t         = value of independent variable where answers are desired           
c             (normally the same as the t last returned by lsode).              
c             for valid results, t must lie between tcur - hu and tcur.         
c             (see optional outputs for tcur and hu.)                           
c k         = integer order of the derivative desired.  k must satisfy          
c             0 .le. k .le. nqcur, where nqcur is the current order             
c             (see optional outputs).  the capability corresponding             
c             to k = 0, i.e. computing y(t), is already provided                
c             by lsode directly.  since nqcur .ge. 1, the first                 
c             derivative dy/dt is always available with intdy.                  
c rwork(21) = the base address of the history array yh.                         
c nyh       = column length of yh, equal to the initial value of neq.           
c                                                                               
c the output parameters are..                                                   
c                                                                               
c dky       = a real array of length neq containing the computed value          
c             of the k-th derivative of y(t).                                   
c iflag     = integer flag, returned as 0 if k and t were legal,                
c             -1 if k was illegal, and -2 if t was illegal.                     
c             on an error return, a message is also written.                    
c-----------------------------------------------------------------------        
c part iii.  common blocks.                                                     
c                                                                               
c if lsode is to be used in an overlay situation, the user                      
c must declare, in the primary overlay, the variables in..                      
c   (1) the call sequence to lsode,                                             
c   (2) the two internal common blocks                                          
c         /ls0001/  of length  257  (218 double precision words                 
c                         followed by 39 integer words),                        
c         /eh0001/  of length  2 (integer words).                               
c                                                                               
c if lsode is used on a system in which the contents of internal                
c common blocks are not preserved between calls, the user should                
c declare the above two common blocks in his main program to insure             
c that their contents are preserved.                                            
c                                                                               
c if the solution of a given problem by lsode is to be interrupted              
c and then later continued, such as when restarting an interrupted run          
c or alternating between two or more problems, the user should save,            
c following the return from the last lsode call prior to the                    
c interruption, the contents of the call sequence variables and the             
c internal common blocks, and later restore these values before the             
c next lsode call for that problem.  to save and restore the common             
c blocks, use subroutine srcom (see part ii above).                             
c                                                                               
c-----------------------------------------------------------------------        
c part iv.  optionally replaceable solver routines.                             
c                                                                               
c below are descriptions of two routines in the lsode package which             
c relate to the measurement of errors.  either routine can be                   
c replaced by a user-supplied version, if desired.  however, since such         
c a replacement may have a major impact on performance, it should be            
c done only when absolutely necessary, and only with great caution.             
c (note.. the means by which the package version of a routine is                
c superseded by the user-s version may be system-dependent.)                    
c                                                                               
c (a) ewset.                                                                    
c the following subroutine is called just before each internal                  
c integration step, and sets the array of error weights, ewt, as                
c described under itol/rtol/atol above..                                        
c     subroutine ewset (neq, itol, rtol, atol, ycur, ewt)                       
c where neq, itol, rtol, and atol are as in the lsode call sequence,            
c ycur contains the current dependent variable vector, and                      
c ewt is the array of weights set by ewset.                                     
c                                                                               
c if the user supplies this subroutine, it must return in ewt(i)                
c (i = 1,...,neq) a positive quantity suitable for comparing errors             
c in y(i) to.  the ewt array returned by ewset is passed to the                 
c vnorm routine (see below), and also used by lsode in the computation          
c of the optional output imxer, the diagonal jacobian approximation,            
c and the increments for difference quotient jacobians.                         
c                                                                               
c in the user-supplied version of ewset, it may be desirable to use             
c the current values of derivatives of y.  derivatives up to order nq           
c are available from the history array yh, described above under                
c optional outputs.  in ewset, yh is identical to the ycur array,               
c extended to nq + 1 columns with a column length of nyh and scale              
c factors of h**j/factorial(j).  on the first call for the problem,             
c given by nst = 0, nq is 1 and h is temporarily set to 1.0.                    
c the quantities nq, nyh, h, and nst can be obtained by including               
c in ewset the statements..                                                     
c     double precision h, rls                                                   
c     common /ls0001/ rls(218),ils(39)                                          
c     nq = ils(35)                                                              
c     nyh = ils(14)                                                             
c     nst = ils(36)                                                             
c     h = rls(212)                                                              
c thus, for example, the current value of dy/dt can be obtained as              
c ycur(nyh+i)/h  (i=1,...,neq)  (and the division by h is                       
c unnecessary when nst = 0).                                                    
c                                                                               
c (b) vnorm.                                                                    
c the following is a real function routine which computes the weighted          
c root-mean-square norm of a vector v..                                         
c     d = vnorm (n, v, w)                                                       
c where..                                                                       
c   n = the length of the vector,                                               
c   v = real array of length n containing the vector,                           
c   w = real array of length n containing weights,                              
c   d = sqrt( (1/n) * sum(v(i)*w(i))**2 ).                                      
c vnorm is called with n = neq and with w(i) = 1.0/ewt(i), where                
c ewt is as set by subroutine ewset.                                            
c                                                                               
c if the user supplies this function, it should return a non-negative           
c value of vnorm suitable for use in the error control in lsode.                
c none of the arguments should be altered by vnorm.                             
c for example, a user-supplied vnorm routine might..                            
c   -substitute a max-norm of (v(i)*w(i)) for the rms-norm, or                  
c   -ignore some components of v in the norm, with the effect of                
c    suppressing the error control on those components of y.                    
c-----------------------------------------------------------------------        
c-----------------------------------------------------------------------        
c other routines in the lsode package.                                          
c                                                                               
c in addition to subroutine lsode, the lsode package includes the               
c following subroutines and function routines..                                 
c  intdy    computes an interpolated value of the y vector at t = tout.         
c  stode    is the core integrator, which does one step of the                  
c           integration and the associated error control.                       
c  cfode    sets all method coefficients and test constants.                    
c  prepj    computes and preprocesses the jacobian matrix j = df/dy             
c           and the newton iteration matrix p = i - h*l0*j.                     
c  solsy    manages solution of linear system in chord iteration.               
c  ewset    sets the error weight vector ewt before each step.                  
c  vnorm    computes the weighted r.m.s. norm of a vector.                      
c  srcom    is a user-callable routine to save and restore                      
c           the contents of the internal common blocks.                         
c  dgefa and dgesl   are routines from linpack for solving full                 
c           systems of linear algebraic equations.                              
c  dgbfa and dgbsl   are routines from linpack for solving banded               
c           linear systems.                                                     
c  daxpy, dscal, idamax, and ddot   are basic linear algebra modules            
c           (blas) used by the above linpack routines.                          
c  d1mach   computes the unit roundoff in a machine-independent manner.         
c  xerrwv, xsetun, and xsetf   handle the printing of all error                 
c           messages and warnings.  xerrwv is machine-dependent.                
c note..  vnorm, idamax, ddot, and d1mach are function routines.                
c all the others are subroutines.                                               
c                                                                               
c the intrinsic and external routines used by lsode are..                       
c dabs, dmax1, dmin1, dfloat, max0, min0, mod, dsign, dsqrt, and write.         
c                                                                               
c a block data subprogram is also included with the package,                    
c for loading some of the variables in internal common.                         
c                                                                               
c-----------------------------------------------------------------------        
c the following card is for optimized compilation on llnl compilers.            
clll. optimize                                                                  
c-----------------------------------------------------------------------        
      external prepj, solsy                                                     
      integer illin, init, lyh, lewt, lacor, lsavf, lwm, liwm,                  
     1   mxstep, mxhnil, nhnil, ntrep, nslast, nyh, iowns                       
      integer icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,           
     1   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu                 
      integer i, i1, i2, iflag, imxer, kgo, lf0,                                
     1   leniw, lenrw, lenwm, ml, mord, mu, mxhnl0, mxstp0                      
      double precision rowns,                                                   
     1   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround                          
      double precision atoli, ayi, big, ewti, h0, hmax, hmx, rh, rtoli,         
     1   tcrit, tdist, tnext, tol, tolsf, tp, size, sum, w0,                    
     2   d1mach, vnorm                                                          
      dimension mord(2)                                                         
      logical ihit                                                              
c-----------------------------------------------------------------------        
c the following internal common block contains                                  
c (a) variables which are local to any subroutine but whose values must         
c     be preserved between calls to the routine (own variables), and            
c (b) variables which are communicated between subroutines.                     
c the structure of the block is as follows..  all real variables are            
c listed first, followed by all integers.  within each type, the                
c variables are grouped with those local to subroutine lsode first,             
c then those local to subroutine stode, and finally those used                  
c for communication.  the block is declared in subroutines                      
c lsode, intdy, stode, prepj, and solsy.  groups of variables are               
c replaced by dummy arrays in the common declarations in routines               
c where those variables are not used.                                           
c-----------------------------------------------------------------------        
      common /ls0001/ rowns(209),                                               
     1   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround,                         
     2   illin, init, lyh, lewt, lacor, lsavf, lwm, liwm,                       
     3   mxstep, mxhnil, nhnil, ntrep, nslast, nyh, iowns(6),                   
     4   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,                
     5   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu                 
c                                                                               
      data  mord(1),mord(2)/12,5/, mxstp0/500/, mxhnl0/10/                      
c-----------------------------------------------------------------------        
c block a.                                                                      
c this code block is executed on every call.                                    
c it tests istate and itask for legality and branches appropriately.            
c if istate .gt. 1 but the flag init shows that initialization has              
c not yet been done, an error return occurs.                                    
c if istate = 1 and tout = t, jump to block g and return immediately.           
c-----------------------------------------------------------------------        
      if (istate .lt. 1 .or. istate .gt. 3) go to 601                           
      if (itask .lt. 1 .or. itask .gt. 5) go to 602                             
      if (istate .eq. 1) go to 10                                               
      if (init .eq. 0) go to 603                                                
      if (istate .eq. 2) go to 200                                              
      go to 20                                                                  
 10   init = 0                                                                  
      if (tout .eq. t) go to 430                                                
 20   ntrep = 0                                                                 
c-----------------------------------------------------------------------        
c block b.                                                                      
c the next code block is executed for the initial call (istate = 1),            
c or for a continuation call with parameter changes (istate = 3).               
c it contains checking of all inputs and various initializations.               
c                                                                               
c first check legality of the non-optional inputs neq, itol, iopt,              
c mf, ml, and mu.                                                               
c-----------------------------------------------------------------------        
      if (neq(1) .le. 0) go to 604                                              
      if (istate .eq. 1) go to 25                                               
      if (neq(1) .gt. n) go to 605                                              
 25   n = neq(1)                                                                
      if (itol .lt. 1 .or. itol .gt. 4) go to 606                               
      if (iopt .lt. 0 .or. iopt .gt. 1) go to 607                               
      meth = mf/10                                                              
      miter = mf - 10*meth                                                      
      if (meth .lt. 1 .or. meth .gt. 2) go to 608                               
      if (miter .lt. 0 .or. miter .gt. 5) go to 608                             
      if (miter .le. 3) go to 30                                                
      ml = iwork(1)                                                             
      mu = iwork(2)                                                             
      if (ml .lt. 0 .or. ml .ge. n) go to 609                                   
      if (mu .lt. 0 .or. mu .ge. n) go to 610                                   
 30   continue                                                                  
c next process and check the optional inputs. --------------------------        
      if (iopt .eq. 1) go to 40                                                 
      maxord = mord(meth)                                                       
      mxstep = mxstp0                                                           
      mxhnil = mxhnl0                                                           
      if (istate .eq. 1) h0 = 0.0d0                                             
      hmxi = 0.0d0                                                              
      hmin = 0.0d0                                                              
      go to 60                                                                  
 40   maxord = iwork(5)                                                         
      if (maxord .lt. 0) go to 611                                              
      if (maxord .eq. 0) maxord = 100                                           
      maxord = min0(maxord,mord(meth))                                          
      mxstep = iwork(6)                                                         
      if (mxstep .lt. 0) go to 612                                              
      if (mxstep .eq. 0) mxstep = mxstp0                                        
      mxhnil = iwork(7)                                                         
      if (mxhnil .lt. 0) go to 613                                              
      if (mxhnil .eq. 0) mxhnil = mxhnl0                                        
      if (istate .ne. 1) go to 50                                               
      h0 = rwork(5)                                                             
      if ((tout - t)*h0 .lt. 0.0d0) go to 614                                   
 50   hmax = rwork(6)                                                           
      if (hmax .lt. 0.0d0) go to 615                                            
      hmxi = 0.0d0                                                              
      if (hmax .gt. 0.0d0) hmxi = 1.0d0/hmax                                    
      hmin = rwork(7)                                                           
      if (hmin .lt. 0.0d0) go to 616                                            
c-----------------------------------------------------------------------        
c set work array pointers and check lengths lrw and liw.                        
c pointers to segments of rwork and iwork are named by prefixing l to           
c the name of the segment.  e.g., the segment yh starts at rwork(lyh).          
c segments of rwork (in order) are denoted  yh, wm, ewt, savf, acor.            
c-----------------------------------------------------------------------        
 60   lyh = 21                                                                  
      if (istate .eq. 1) nyh = n                                                
      lwm = lyh + (maxord + 1)*nyh                                              
      if (miter .eq. 0) lenwm = 0                                               
      if (miter .eq. 1 .or. miter .eq. 2) lenwm = n*n + 2                       
      if (miter .eq. 3) lenwm = n + 2                                           
      if (miter .ge. 4) lenwm = (2*ml + mu + 1)*n + 2                           
      lewt = lwm + lenwm                                                        
      lsavf = lewt + n                                                          
      lacor = lsavf + n                                                         
      lenrw = lacor + n - 1                                                     
      iwork(17) = lenrw                                                         
      liwm = 1                                                                  
      leniw = 20 + n                                                            
      if (miter .eq. 0 .or. miter .eq. 3) leniw = 20                            
      iwork(18) = leniw                                                         
      if (lenrw .gt. lrw) go to 617                                             
      if (leniw .gt. liw) go to 618                                             
c check rtol and atol for legality. ------------------------------------        
      rtoli = rtol(1)                                                           
      atoli = atol(1)                                                           
      do 70 i = 1,n                                                             
        if (itol .ge. 3) rtoli = rtol(i)                                        
        if (itol .eq. 2 .or. itol .eq. 4) atoli = atol(i)                       
        if (rtoli .lt. 0.0d0) go to 619                                         
        if (atoli .lt. 0.0d0) go to 620                                         
 70     continue                                                                
      if (istate .eq. 1) go to 100                                              
c if istate = 3, set flag to signal parameter changes to stode. --------        
      jstart = -1                                                               
      if (nq .le. maxord) go to 90                                              
c maxord was reduced below nq.  copy yh(*,maxord+2) into savf. ---------        
      do 80 i = 1,n                                                             
 80     rwork(i+lsavf-1) = rwork(i+lwm-1)                                       
c reload wm(1) = rwork(lwm), since lwm may have changed. ---------------        
 90   if (miter .gt. 0) rwork(lwm) = dsqrt(uround)                              
      if (n .eq. nyh) go to 200                                                 
c neq was reduced.  zero part of yh to avoid undefined references. -----        
      i1 = lyh + l*nyh                                                          
      i2 = lyh + (maxord + 1)*nyh - 1                                           
      if (i1 .gt. i2) go to 200                                                 
      do 95 i = i1,i2                                                           
 95     rwork(i) = 0.0d0                                                        
      go to 200                                                                 
c-----------------------------------------------------------------------        
c block c.                                                                      
c the next block is for the initial call only (istate = 1).                     
c it contains all remaining initializations, the initial call to f,             
c and the calculation of the initial step size.                                 
c the error weights in ewt are inverted after being loaded.                     
c-----------------------------------------------------------------------        
c100  uround = d1mach(4)                                                        
 100  uround = 1.d-14                                                           
      tn = t                                                                    
      if (itask .ne. 4 .and. itask .ne. 5) go to 110                            
      tcrit = rwork(1)                                                          
      if ((tcrit - tout)*(tout - t) .lt. 0.0d0) go to 625                       
      if (h0 .ne. 0.0d0 .and. (t + h0 - tcrit)*h0 .gt. 0.0d0)                   
     1   h0 = tcrit - t                                                         
 110  jstart = 0                                                                
      if (miter .gt. 0) rwork(lwm) = dsqrt(uround)                              
      nhnil = 0                                                                 
      nst = 0                                                                   
      nje = 0                                                                   
      nslast = 0                                                                
      hu = 0.0d0                                                                
      nqu = 0                                                                   
      ccmax = 0.3d0                                                             
      maxcor = 3                                                                
      msbp = 20                                                                 
      mxncf = 10                                                                
c initial call to f.  (lf0 points to yh(*,2).) -------------------------        
      lf0 = lyh + nyh                                                           
      call f (neq, t, y, rwork(lf0))                                            
      nfe = 1                                                                   
c load the initial value vector in yh. ---------------------------------        
      do 115 i = 1,n                                                            
 115    rwork(i+lyh-1) = y(i)                                                   
c load and invert the ewt array.  (h is temporarily set to 1.0.) -------        
      nq = 1                                                                    
      h = 1.0d0                                                                 
      call ewset (n, itol, rtol, atol, rwork(lyh), rwork(lewt))                 
      do 120 i = 1,n                                                            
        if (rwork(i+lewt-1) .le. 0.0d0) go to 621                               
 120    rwork(i+lewt-1) = 1.0d0/rwork(i+lewt-1)                                 
c-----------------------------------------------------------------------        
c the coding below computes the step size, h0, to be attempted on the           
c first step, unless the user has supplied a value for this.                    
c first check that tout - t differs significantly from zero.                    
c a scalar tolerance quantity tol is computed, as max(rtol(i))                  
c if this is positive, or max(atol(i)/abs(y(i))) otherwise, adjusted            
c so as to be between 100*uround and 1.0e-3.                                    
c then the computed value h0 is given by..                                      
c                                      neq                                      
c   h0**2 = tol / ( w0**-2 + (1/neq) * sum ( f(i)/ywt(i) )**2  )                
c                                       1                                       
c where   w0     = max ( abs(t), abs(tout) ),                                   
c         f(i)   = i-th component of initial value of f,                        
c         ywt(i) = ewt(i)/tol  (a weight for y(i)).                             
c the sign of h0 is inferred from the initial values of tout and t.             
c-----------------------------------------------------------------------        
      if (h0 .ne. 0.0d0) go to 180                                              
      tdist = dabs(tout - t)                                                    
      w0 = dmax1(dabs(t),dabs(tout))                                            
      if (tdist .lt. 2.0d0*uround*w0) go to 622                                 
      tol = rtol(1)                                                             
      if (itol .le. 2) go to 140                                                
      do 130 i = 1,n                                                            
 130    tol = dmax1(tol,rtol(i))                                                
 140  if (tol .gt. 0.0d0) go to 160                                             
      atoli = atol(1)                                                           
      do 150 i = 1,n                                                            
        if (itol .eq. 2 .or. itol .eq. 4) atoli = atol(i)                       
        ayi = dabs(y(i))                                                        
        if (ayi .ne. 0.0d0) tol = dmax1(tol,atoli/ayi)                          
 150    continue                                                                
 160  tol = dmax1(tol,100.0d0*uround)                                           
      tol = dmin1(tol,0.001d0)                                                  
      sum = vnorm (n, rwork(lf0), rwork(lewt))                                  
      sum = 1.0d0/(tol*w0*w0) + tol*sum**2                                      
      h0 = 1.0d0/dsqrt(sum)                                                     
      h0 = dmin1(h0,tdist)                                                      
      h0 = dsign(h0,tout-t)                                                     
c adjust h0 if necessary to meet hmax bound. ---------------------------        
 180  rh = dabs(h0)*hmxi                                                        
      if (rh .gt. 1.0d0) h0 = h0/rh                                             
c load h with h0 and scale yh(*,2) by h0. ------------------------------        
      h = h0                                                                    
      do 190 i = 1,n                                                            
 190    rwork(i+lf0-1) = h0*rwork(i+lf0-1)                                      
      go to 270                                                                 
c-----------------------------------------------------------------------        
c block d.                                                                      
c the next code block is for continuation calls only (istate = 2 or 3)          
c and is to check stop conditions before taking a step.                         
c-----------------------------------------------------------------------        
 200  nslast = nst                                                              
      go to (210, 250, 220, 230, 240), itask                                    
 210  if ((tn - tout)*h .lt. 0.0d0) go to 250                                   
      call intdy (tout, 0, rwork(lyh), nyh, y, iflag)                           
      if (iflag .ne. 0) go to 627                                               
      t = tout                                                                  
      go to 420                                                                 
 220  tp = tn - hu*(1.0d0 + 100.0d0*uround)                                     
      if ((tp - tout)*h .gt. 0.0d0) go to 623                                   
      if ((tn - tout)*h .lt. 0.0d0) go to 250                                   
      go to 400                                                                 
 230  tcrit = rwork(1)                                                          
      if ((tn - tcrit)*h .gt. 0.0d0) go to 624                                  
      if ((tcrit - tout)*h .lt. 0.0d0) go to 625                                
      if ((tn - tout)*h .lt. 0.0d0) go to 245                                   
      call intdy (tout, 0, rwork(lyh), nyh, y, iflag)                           
      if (iflag .ne. 0) go to 627                                               
      t = tout                                                                  
      go to 420                                                                 
 240  tcrit = rwork(1)                                                          
      if ((tn - tcrit)*h .gt. 0.0d0) go to 624                                  
 245  hmx = dabs(tn) + dabs(h)                                                  
      ihit = dabs(tn - tcrit) .le. 100.0d0*uround*hmx                           
      if (ihit) go to 400                                                       
      tnext = tn + h*(1.0d0 + 4.0d0*uround)                                     
      if ((tnext - tcrit)*h .le. 0.0d0) go to 250                               
      h = (tcrit - tn)*(1.0d0 - 4.0d0*uround)                                   
      if (istate .eq. 2) jstart = -2                                            
c-----------------------------------------------------------------------        
c block e.                                                                      
c the next block is normally executed for all calls and contains                
c the call to the one-step core integrator stode.                               
c                                                                               
c this is a looping point for the integration steps.                            
c                                                                               
c first check for too many steps being taken, update ewt (if not at             
c start of problem), check for too much accuracy being requested, and           
c check for h below the roundoff level in t.                                    
c-----------------------------------------------------------------------        
 250  continue                                                                  
      if ((nst-nslast) .ge. mxstep) go to 500                                   
      call ewset (n, itol, rtol, atol, rwork(lyh), rwork(lewt))                 
      do 260 i = 1,n                                                            
        if (rwork(i+lewt-1) .le. 0.0d0) go to 510                               
 260    rwork(i+lewt-1) = 1.0d0/rwork(i+lewt-1)                                 
 270  tolsf = uround*vnorm (n, rwork(lyh), rwork(lewt))                         
      if (tolsf .le. 1.0d0) go to 280                                           
      tolsf = tolsf*2.0d0                                                       
      if (nst .eq. 0) go to 626                                                 
      go to 520                                                                 
 280  if ((tn + h) .ne. tn) go to 290                                           
      nhnil = nhnil + 1                                                         
      if (nhnil .gt. mxhnil) go to 290                                          
      call xerrwv(50hlsode--  warning..internal t (=r1) and h (=r2) are,        
     1   50, 101, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)                                  
      call xerrwv(                                                              
     1  60h      such that in the machine, t + h = t on the next step  ,        
     1   60, 101, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)                                  
      call xerrwv(50h      (h = step size). solver will continue anyway,        
     1   50, 101, 0, 0, 0, 0, 2, tn, h)                                         
      if (nhnil .lt. mxhnil) go to 290                                          
      call xerrwv(50hlsode--  above warning has been issued i1 times.  ,        
     1   50, 102, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)                                  
      call xerrwv(50h      it will not be issued again for this problem,        
     1   50, 102, 0, 1, mxhnil, 0, 0, 0.0d0, 0.0d0)                             
 290  continue                                                                  
c-----------------------------------------------------------------------        
c     call stode(neq,y,yh,nyh,yh,ewt,savf,acor,wm,iwm,f,jac,prepj,solsy)        
c-----------------------------------------------------------------------        
      call stode (neq, y, rwork(lyh), nyh, rwork(lyh), rwork(lewt),             
     1   rwork(lsavf), rwork(lacor), rwork(lwm), iwork(liwm),                   
     2   f, jac, prepj, solsy)                                                  
      kgo = 1 - kflag                                                           
      go to (300, 530, 540), kgo                                                
c-----------------------------------------------------------------------        
c block f.                                                                      
c the following block handles the case of a successful return from the          
c core integrator (kflag = 0).  test for stop conditions.                       
c-----------------------------------------------------------------------        
 300  init = 1                                                                  
      go to (310, 400, 330, 340, 350), itask                                    
c itask = 1.  if tout has been reached, interpolate. -------------------        
 310  if ((tn - tout)*h .lt. 0.0d0) go to 250                                   
      call intdy (tout, 0, rwork(lyh), nyh, y, iflag)                           
      t = tout                                                                  
      go to 420                                                                 
c itask = 3.  jump to exit if tout was reached. ------------------------        
 330  if ((tn - tout)*h .ge. 0.0d0) go to 400                                   
      go to 250                                                                 
c itask = 4.  see if tout or tcrit was reached.  adjust h if necessary.         
 340  if ((tn - tout)*h .lt. 0.0d0) go to 345                                   
      call intdy (tout, 0, rwork(lyh), nyh, y, iflag)                           
      t = tout                                                                  
      go to 420                                                                 
 345  hmx = dabs(tn) + dabs(h)                                                  
      ihit = dabs(tn - tcrit) .le. 100.0d0*uround*hmx                           
      if (ihit) go to 400                                                       
      tnext = tn + h*(1.0d0 + 4.0d0*uround)                                     
      if ((tnext - tcrit)*h .le. 0.0d0) go to 250                               
      h = (tcrit - tn)*(1.0d0 - 4.0d0*uround)                                   
      jstart = -2                                                               
      go to 250                                                                 
c itask = 5.  see if tcrit was reached and jump to exit. ---------------        
 350  hmx = dabs(tn) + dabs(h)                                                  
      ihit = dabs(tn - tcrit) .le. 100.0d0*uround*hmx                           
c-----------------------------------------------------------------------        
c block g.                                                                      
c the following block handles all successful returns from lsode.                
c if itask .ne. 1, y is loaded from yh and t is set accordingly.                
c istate is set to 2, the illegal input counter is zeroed, and the              
c optional outputs are loaded into the work arrays before returning.            
c if istate = 1 and tout = t, there is a return with no action taken,           
c except that if this has happened repeatedly, the run is terminated.           
c-----------------------------------------------------------------------        
 400  do 410 i = 1,n                                                            
 410    y(i) = rwork(i+lyh-1)                                                   
      t = tn                                                                    
      if (itask .ne. 4 .and. itask .ne. 5) go to 420                            
      if (ihit) t = tcrit                                                       
 420  istate = 2                                                                
      illin = 0                                                                 
      rwork(11) = hu                                                            
      rwork(12) = h                                                             
      rwork(13) = tn                                                            
      iwork(11) = nst                                                           
      iwork(12) = nfe                                                           
      iwork(13) = nje                                                           
      iwork(14) = nqu                                                           
      iwork(15) = nq                                                            
      return                                                                    
c                                                                               
 430  ntrep = ntrep + 1                                                         
      if (ntrep .lt. 5) return                                                  
      call xerrwv(                                                              
     1  60hlsode--  repeated calls with istate = 1 and tout = t (=r1)  ,        
     1   60, 301, 0, 0, 0, 0, 1, t, 0.0d0)                                      
      go to 800                                                                 
c-----------------------------------------------------------------------        
c block h.                                                                      
c the following block handles all unsuccessful returns other than               
c those for illegal input.  first the error message routine is called.          
c if there was an error test or convergence test failure, imxer is set.         
c then y is loaded from yh, t is set to tn, and the illegal input               
c counter illin is set to 0.  the optional outputs are loaded into              
c the work arrays before returning.                                             
c-----------------------------------------------------------------------        
c the maximum number of steps was taken before reaching tout. ----------        
 500  call xerrwv(50hlsode--  at current t (=r1), mxstep (=i1) steps   ,        
     1   50, 201, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)                                  
      call xerrwv(50h      taken on this call before reaching tout     ,        
     1   50, 201, 0, 1, mxstep, 0, 1, tn, 0.0d0)                                
      istate = -1                                                               
      go to 580                                                                 
c ewt(i) .le. 0.0 for some i (not at start of problem). ----------------        
 510  ewti = rwork(lewt+i-1)                                                    
      call xerrwv(50hlsode--  at t (=r1), ewt(i1) has become r2 .le. 0.,        
     1   50, 202, 0, 1, i, 0, 2, tn, ewti)                                      
      istate = -6                                                               
      go to 580                                                                 
c too much accuracy requested for machine precision. -------------------        
 520  call xerrwv(50hlsode--  at t (=r1), too much accuracy requested  ,        
     1   50, 203, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)                                  
      call xerrwv(50h      for precision of machine..  see tolsf (=r2) ,        
     1   50, 203, 0, 0, 0, 0, 2, tn, tolsf)                                     
      rwork(14) = tolsf                                                         
      istate = -2                                                               
      go to 580                                                                 
c kflag = -1.  error test failed repeatedly or with abs(h) = hmin. -----        
 530  call xerrwv(50hlsode--  at t(=r1) and step size h(=r2), the error,        
     1   50, 204, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)                                  
      call xerrwv(50h      test failed repeatedly or with abs(h) = hmin,        
     1   50, 204, 0, 0, 0, 0, 2, tn, h)                                         
      istate = -4                                                               
      go to 560                                                                 
c kflag = -2.  convergence failed repeatedly or with abs(h) = hmin. ----        
 540  call xerrwv(50hlsode--  at t (=r1) and step size h (=r2), the    ,        
     1   50, 205, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)                                  
      call xerrwv(50h      corrector convergence failed repeatedly     ,        
     1   50, 205, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)                                  
      call xerrwv(30h      or with abs(h) = hmin   ,                            
     1   30, 205, 0, 0, 0, 0, 2, tn, h)                                         
      istate = -5                                                               
c compute imxer if relevant. -------------------------------------------        
 560  big = 0.0d0                                                               
      imxer = 1                                                                 
      do 570 i = 1,n                                                            
        size = dabs(rwork(i+lacor-1)*rwork(i+lewt-1))                           
        if (big .ge. size) go to 570                                            
        big = size                                                              
        imxer = i                                                               
 570    continue                                                                
      iwork(16) = imxer                                                         
c set y vector, t, illin, and optional outputs. ------------------------        
 580  do 590 i = 1,n                                                            
 590    y(i) = rwork(i+lyh-1)                                                   
      t = tn                                                                    
      illin = 0                                                                 
      rwork(11) = hu                                                            
      rwork(12) = h                                                             
      rwork(13) = tn                                                            
      iwork(11) = nst                                                           
      iwork(12) = nfe                                                           
      iwork(13) = nje                                                           
      iwork(14) = nqu                                                           
      iwork(15) = nq                                                            
      return                                                                    
c-----------------------------------------------------------------------        
c block i.                                                                      
c the following block handles all error returns due to illegal input            
c (istate = -3), as detected before calling the core integrator.                
c first the error message routine is called.  then if there have been           
c 5 consecutive such returns just before this call to the solver,               
c the run is halted.                                                            
c-----------------------------------------------------------------------        
 601  call xerrwv(30hlsode--  istate (=i1) illegal ,                            
     1   30, 1, 0, 1, istate, 0, 0, 0.0d0, 0.0d0)                               
      go to 700                                                                 
 602  call xerrwv(30hlsode--  itask (=i1) illegal  ,                            
     1   30, 2, 0, 1, itask, 0, 0, 0.0d0, 0.0d0)                                
      go to 700                                                                 
 603  call xerrwv(50hlsode--  istate .gt. 1 but lsode not initialized  ,        
     1   50, 3, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)                                    
      go to 700                                                                 
 604  call xerrwv(30hlsode--  neq (=i1) .lt. 1     ,                            
     1   30, 4, 0, 1, neq(1), 0, 0, 0.0d0, 0.0d0)                               
      go to 700                                                                 
 605  call xerrwv(50hlsode--  istate = 3 and neq increased (i1 to i2)  ,        
     1   50, 5, 0, 2, n, neq(1), 0, 0.0d0, 0.0d0)                               
      go to 700                                                                 
 606  call xerrwv(30hlsode--  itol (=i1) illegal   ,                            
     1   30, 6, 0, 1, itol, 0, 0, 0.0d0, 0.0d0)                                 
      go to 700                                                                 
 607  call xerrwv(30hlsode--  iopt (=i1) illegal   ,                            
     1   30, 7, 0, 1, iopt, 0, 0, 0.0d0, 0.0d0)                                 
      go to 700                                                                 
 608  call xerrwv(30hlsode--  mf (=i1) illegal     ,                            
     1   30, 8, 0, 1, mf, 0, 0, 0.0d0, 0.0d0)                                   
      go to 700                                                                 
 609  call xerrwv(50hlsode--  ml (=i1) illegal.. .lt.0 or .ge.neq (=i2),        
     1   50, 9, 0, 2, ml, neq(1), 0, 0.0d0, 0.0d0)                              
      go to 700                                                                 
 610  call xerrwv(50hlsode--  mu (=i1) illegal.. .lt.0 or .ge.neq (=i2),        
     1   50, 10, 0, 2, mu, neq(1), 0, 0.0d0, 0.0d0)                             
      go to 700                                                                 
 611  call xerrwv(30hlsode--  maxord (=i1) .lt. 0  ,                            
     1   30, 11, 0, 1, maxord, 0, 0, 0.0d0, 0.0d0)                              
      go to 700                                                                 
 612  call xerrwv(30hlsode--  mxstep (=i1) .lt. 0  ,                            
     1   30, 12, 0, 1, mxstep, 0, 0, 0.0d0, 0.0d0)                              
      go to 700                                                                 
 613  call xerrwv(30hlsode--  mxhnil (=i1) .lt. 0  ,                            
     1   30, 13, 0, 1, mxhnil, 0, 0, 0.0d0, 0.0d0)                              
      go to 700                                                                 
 614  call xerrwv(40hlsode--  tout (=r1) behind t (=r2)      ,                  
     1   40, 14, 0, 0, 0, 0, 2, tout, t)                                        
      call xerrwv(50h      integration direction is given by h0 (=r1)  ,        
     1   50, 14, 0, 0, 0, 0, 1, h0, 0.0d0)                                      
      go to 700                                                                 
 615  call xerrwv(30hlsode--  hmax (=r1) .lt. 0.0  ,                            
     1   30, 15, 0, 0, 0, 0, 1, hmax, 0.0d0)                                    
      go to 700                                                                 
 616  call xerrwv(30hlsode--  hmin (=r1) .lt. 0.0  ,                            
     1   30, 16, 0, 0, 0, 0, 1, hmin, 0.0d0)                                    
      go to 700                                                                 
 617  call xerrwv(                                                              
     1  60hlsode--  rwork length needed, lenrw (=i1), exceeds lrw (=i2),        
     1   60, 17, 0, 2, lenrw, lrw, 0, 0.0d0, 0.0d0)                             
      go to 700                                                                 
 618  call xerrwv(                                                              
     1  60hlsode--  iwork length needed, leniw (=i1), exceeds liw (=i2),        
     1   60, 18, 0, 2, leniw, liw, 0, 0.0d0, 0.0d0)                             
      go to 700                                                                 
 619  call xerrwv(40hlsode--  rtol(i1) is r1 .lt. 0.0        ,                  
     1   40, 19, 0, 1, i, 0, 1, rtoli, 0.0d0)                                   
      go to 700                                                                 
 620  call xerrwv(40hlsode--  atol(i1) is r1 .lt. 0.0        ,                  
     1   40, 20, 0, 1, i, 0, 1, atoli, 0.0d0)                                   
      go to 700                                                                 
 621  ewti = rwork(lewt+i-1)                                                    
      call xerrwv(40hlsode--  ewt(i1) is r1 .le. 0.0         ,                  
     1   40, 21, 0, 1, i, 0, 1, ewti, 0.0d0)                                    
      go to 700                                                                 
 622  call xerrwv(                                                              
     1  60hlsode--  tout (=r1) too close to t(=r2) to start integration,        
     1   60, 22, 0, 0, 0, 0, 2, tout, t)                                        
      go to 700                                                                 
 623  call xerrwv(                                                              
     1  60hlsode--  itask = i1 and tout (=r1) behind tcur - hu (= r2)  ,        
     1   60, 23, 0, 1, itask, 0, 2, tout, tp)                                   
      go to 700                                                                 
 624  call xerrwv(                                                              
     1  60hlsode--  itask = 4 or 5 and tcrit (=r1) behind tcur (=r2)   ,        
     1   60, 24, 0, 0, 0, 0, 2, tcrit, tn)                                      
      go to 700                                                                 
 625  call xerrwv(                                                              
     1  60hlsode--  itask = 4 or 5 and tcrit (=r1) behind tout (=r2)   ,        
     1   60, 25, 0, 0, 0, 0, 2, tcrit, tout)                                    
      go to 700                                                                 
 626  call xerrwv(50hlsode--  at start of problem, too much accuracy   ,        
     1   50, 26, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)                                   
      call xerrwv(                                                              
     1  60h      requested for precision of machine..  see tolsf (=r1) ,        
     1   60, 26, 0, 0, 0, 0, 1, tolsf, 0.0d0)                                   
      rwork(14) = tolsf                                                         
      go to 700                                                                 
 627  call xerrwv(50hlsode--  trouble from intdy. itask = i1, tout = r1,        
     1   50, 27, 0, 1, itask, 0, 1, tout, 0.0d0)                                
c                                                                               
 700  if (illin .eq. 5) go to 710                                               
      illin = illin + 1                                                         
      istate = -3                                                               
      return                                                                    
 710  call xerrwv(50hlsode--  repeated occurrences of illegal input    ,        
     1   50, 302, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)                                  
c                                                                               
 800  call xerrwv(50hlsode--  run aborted.. apparent infinite loop     ,        
     1   50, 303, 2, 0, 0, 0, 0, 0.0d0, 0.0d0)                                  
      return                                                                    
c----------------------- end of subroutine lsode -----------------------        
      end                                                                       
      subroutine solsy (wm, iwm, x, tem)                                        
clll. optimize                                                                  
      integer iwm                                                               
      integer iownd, iowns,                                                     
     1   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,                
     2   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu                 
      integer i, meband, ml, mu                                                 
      double precision wm, x, tem                                               
      double precision rowns,                                                   
     1   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround                          
      double precision di, hl0, phl0, r                                         
      dimension wm(100), iwm(100), x(1), tem(1)                                     
      common /ls0001/ rowns(209),                                               
     2   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround,                         
     3   iownd(14), iowns(6),                                                   
     4   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,                
     5   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu                 
c-----------------------------------------------------------------------        
c this routine manages the solution of the linear system arising from           
c a chord iteration.  it is called if miter .ne. 0.                             
c if miter is 1 or 2, it calls dgesl to accomplish this.                        
c if miter = 3 it updates the coefficient h*el0 in the diagonal                 
c matrix, and then computes the solution.                                       
c if miter is 4 or 5, it calls dgbsl.                                           
c communication with solsy uses the following variables..                       
c wm    = real work space containing the inverse diagonal matrix if             
c         miter = 3 and the lu decomposition of the matrix otherwise.           
c         storage of matrix elements starts at wm(3).                           
c         wm also contains the following matrix-related data..                  
c         wm(1) = sqrt(uround) (not used here),                                 
c         wm(2) = hl0, the previous value of h*el0, used if miter = 3.          
c iwm   = integer work space containing pivot information, starting at          
c         iwm(21), if miter is 1, 2, 4, or 5.  iwm also contains band           
c         parameters ml = iwm(1) and mu = iwm(2) if miter is 4 or 5.            
c x     = the right-hand side vector on input, and the solution vector          
c         on output, of length n.                                               
c tem   = vector of work space of length n, not used in this version.           
c iersl = output flag (in common).  iersl = 0 if no trouble occurred.           
c         iersl = 1 if a singular matrix arose with miter = 3.                  
c this routine also uses the common variables el0, h, miter, and n.             
c-----------------------------------------------------------------------        
      iersl = 0                                                                 
      go to (100, 100, 300, 400, 400), miter                                    
 100  call dgesl (wm(3), n, n, iwm(21), x, 0)                                   
      return                                                                    
c                                                                               
 300  phl0 = wm(2)                                                              
      hl0 = h*el0                                                               
      wm(2) = hl0                                                               
      if (hl0 .eq. phl0) go to 330                                              
      r = hl0/phl0                                                              
      do 320 i = 1,n                                                            
        di = 1.0d0 - r*(1.0d0 - 1.0d0/wm(i+2))                                  
        if (dabs(di) .eq. 0.0d0) go to 390                                      
 320    wm(i+2) = 1.0d0/di                                                      
 330  do 340 i = 1,n                                                            
 340    x(i) = wm(i+2)*x(i)                                                     
      return                                                                    
 390  iersl = 1                                                                 
      return                                                                    
c                                                                               
 400  ml = iwm(1)                                                               
      mu = iwm(2)                                                               
      meband = 2*ml + mu + 1                                                    
C     call dgbsl (wm(3), meband, n, ml, mu, iwm(21), x, 0)                      
      return                                                                    
c----------------------- end of subroutine solsy -----------------------        
      end                                                                       
      subroutine prepj (neq, y, yh, nyh, ewt, ftem, savf, wm, iwm,              
     1   f, jac)                                                                
clll. optimize                                                                  
      external f, jac                                                           
      integer neq, nyh, iwm                                                     
      integer iownd, iowns,                                                     
     1   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,                
     2   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu                 
      integer i, i1, i2, ier, ii, j, j1, jj, lenp,                              
     1   mba, mband, meb1, meband, ml, ml3, mu, np1                             
      double precision y, yh, ewt, ftem, savf, wm                               
      double precision rowns,                                                   
     1   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround                          
      double precision con, di, fac, hl0, r, r0, srur, yi, yj, yjj,             
     1   vnorm                                                                  
      dimension neq(1), y(1), yh(nyh,100), ewt(1), ftem(1), savf(100),              
     1   wm(100), iwm(100)                                                          
      common /ls0001/ rowns(209),                                               
     2   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround,                         
     3   iownd(14), iowns(6),                                                   
     4   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,                
     5   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu                 
c-----------------------------------------------------------------------        
c prepj is called by stode to compute and process the matrix                    
c p = i - h*el(1)*j , where j is an approximation to the jacobian.              
c here j is computed by the user-supplied routine jac if                        
c miter = 1 or 4, or by finite differencing if miter = 2, 3, or 5.              
c if miter = 3, a diagonal approximation to j is used.                          
c j is stored in wm and replaced by p.  if miter .ne. 3, p is then              
c subjected to lu decomposition in preparation for later solution               
c of linear systems with p as coefficient matrix. this is done                  
c by dgefa if miter = 1 or 2, and by dgbfa if miter = 4 or 5.                   
c                                                                               
c in addition to variables described previously, communication                  
c with prepj uses the following..                                               
c y     = array containing predicted values on entry.                           
c ftem  = work array of length n (acor in stode).                               
c savf  = array containing f evaluated at predicted y.                          
c wm    = real work space for matrices.  on output it contains the              
c         inverse diagonal matrix if miter = 3 and the lu decomposition         
c         of p if miter is 1, 2 , 4, or 5.                                      
c         storage of matrix elements starts at wm(3).                           
c         wm also contains the following matrix-related data..                  
c         wm(1) = sqrt(uround), used in numerical jacobian increments.          
c         wm(2) = h*el0, saved for later use if miter = 3.                      
c iwm   = integer work space containing pivot information, starting at          
c         iwm(21), if miter is 1, 2, 4, or 5.  iwm also contains band           
c         parameters ml = iwm(1) and mu = iwm(2) if miter is 4 or 5.            
c el0   = el(1) (input).                                                        
c ierpj = output error flag,  = 0 if no trouble, .gt. 0 if                      
c         p matrix found to be singular.                                        
c jcur  = output flag = 1 to indicate that the jacobian matrix                  
c         (or approximation) is now current.                                    
c this routine also uses the common variables el0, h, tn, uround,               
c miter, n, nfe, and nje.                                                       
c-----------------------------------------------------------------------        
      nje = nje + 1                                                             
      ierpj = 0                                                                 
      jcur = 1                                                                  
      hl0 = h*el0                                                               
      go to (100, 200, 300, 400, 500), miter                                    
c if miter = 1, call jac and multiply by scalar. -----------------------        
 100  lenp = n*n                                                                
      do 110 i = 1,lenp                                                         
 110    wm(i+2) = 0.0d0                                                         
      call jac (neq, tn, y, 0, 0, wm(3), n)                                     
      con = -hl0                                                                
      do 120 i = 1,lenp                                                         
 120    wm(i+2) = wm(i+2)*con                                                   
      go to 240                                                                 
c if miter = 2, make n calls to f to approximate j. --------------------        
 200  fac = vnorm (n, savf, ewt)                                                
      r0 = 1000.0d0*dabs(h)*uround*dfloat(n)*fac                                
      if (r0 .eq. 0.0d0) r0 = 1.0d0                                             
      srur = wm(1)                                                              
      j1 = 2                                                                    
      do 230 j = 1,n                                                            
        yj = y(j)                                                               
        r = dmax1(srur*dabs(yj),r0/ewt(j))                                      
        y(j) = y(j) + r                                                         
        fac = -hl0/r                                                            
        call f (neq, tn, y, ftem)                                               
        do 220 i = 1,n                                                          
 220      wm(i+j1) = (ftem(i) - savf(i))*fac                                    
        y(j) = yj                                                               
        j1 = j1 + n                                                             
 230    continue                                                                
      nfe = nfe + n                                                             
c add identity matrix. -------------------------------------------------        
 240  j = 3                                                                     
      np1 = n + 1                                                               
      do 250 i = 1,n                                                            
        wm(j) = wm(j) + 1.0d0                                                   
 250    j = j + np1                                                             
c do lu decomposition on p. --------------------------------------------        
      call dgefa (wm(3), n, n, iwm(21), ier)                                    
      if (ier .ne. 0) ierpj = 1                                                 
      return                                                                    
c if miter = 3, construct a diagonal approximation to j and p. ---------        
 300  wm(2) = hl0                                                               
      r = el0*0.1d0                                                             
      do 310 i = 1,n                                                            
 310    y(i) = y(i) + r*(h*savf(i) - yh(i,2))                                   
      call f (neq, tn, y, wm(3))                                                
      nfe = nfe + 1                                                             
      do 320 i = 1,n                                                            
        r0 = h*savf(i) - yh(i,2)                                                
        di = 0.1d0*r0 - h*(wm(i+2) - savf(i))                                   
        wm(i+2) = 1.0d0                                                         
        if (dabs(r0) .lt. uround/ewt(i)) go to 320                              
        if (dabs(di) .eq. 0.0d0) go to 330                                      
        wm(i+2) = 0.1d0*r0/di                                                   
 320    continue                                                                
      return                                                                    
 330  ierpj = 1                                                                 
      return                                                                    
c if miter = 4, call jac and multiply by scalar. -----------------------        
 400  ml = iwm(1)                                                               
      mu = iwm(2)                                                               
      ml3 = ml + 3                                                              
      mband = ml + mu + 1                                                       
      meband = mband + ml                                                       
      lenp = meband*n                                                           
      do 410 i = 1,lenp                                                         
 410    wm(i+2) = 0.0d0                                                         
      call jac (neq, tn, y, ml, mu, wm(ml3), meband)                            
      con = -hl0                                                                
      do 420 i = 1,lenp                                                         
 420    wm(i+2) = wm(i+2)*con                                                   
      go to 570                                                                 
c if miter = 5, make mband calls to f to approximate j. ----------------        
 500  ml = iwm(1)                                                               
      mu = iwm(2)                                                               
      mband = ml + mu + 1                                                       
      mba = min0(mband,n)                                                       
      meband = mband + ml                                                       
      meb1 = meband - 1                                                         
      srur = wm(1)                                                              
      fac = vnorm (n, savf, ewt)                                                
      r0 = 1000.0d0*dabs(h)*uround*dfloat(n)*fac                                
      if (r0 .eq. 0.0d0) r0 = 1.0d0                                             
      do 560 j = 1,mba                                                          
        do 530 i = j,n,mband                                                    
          yi = y(i)                                                             
          r = dmax1(srur*dabs(yi),r0/ewt(i))                                    
 530      y(i) = y(i) + r                                                       
        call f (neq, tn, y, ftem)                                               
        do 550 jj = j,n,mband                                                   
          y(jj) = yh(jj,1)                                                      
          yjj = y(jj)                                                           
          r = dmax1(srur*dabs(yjj),r0/ewt(jj))                                  
          fac = -hl0/r                                                          
          i1 = max0(jj-mu,1)                                                    
          i2 = min0(jj+ml,n)                                                    
          ii = jj*meb1 - ml + 2                                                 
          do 540 i = i1,i2                                                      
 540        wm(ii+i) = (ftem(i) - savf(i))*fac                                  
 550      continue                                                              
 560    continue                                                                
      nfe = nfe + mba                                                           
c add identity matrix. -------------------------------------------------        
 570  ii = mband + 2                                                            
      do 580 i = 1,n                                                            
        wm(ii) = wm(ii) + 1.0d0                                                 
 580    ii = ii + meband                                                        
c do lu decomposition of p. --------------------------------------------        
C     call dgbfa (wm(3), meband, n, ml, mu, iwm(21), ier)                       
      if (ier .ne. 0) ierpj = 1                                                 
      return                                                                    
c----------------------- end of subroutine prepj -----------------------        
      end                                                                       
      subroutine stode (neq, y, yh, nyh, yh1, ewt, savf, acor,                  
     1   wm, iwm, f, jac, pjac, slvs)                                           
clll. optimize                                                                  
      external f, jac, pjac, slvs                                               
      integer neq, nyh, iwm                                                     
      integer iownd, ialth, ipup, lmax, meo, nqnyh, nslp,                       
     1   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,                
     2   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu                 
      integer i, i1, iredo, iret, j, jb, m, ncf, newq                           
      double precision y, yh, yh1, ewt, savf, acor, wm                          
      double precision conit, crate, el, elco, hold, rmax, tesco,               
     2   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround                          
      double precision dcon, ddn, del, delp, dsm, dup, exdn, exsm, exup,        
     1   r, rh, rhdn, rhsm, rhup, told, vnorm                                   
      dimension neq(1), y(1), yh(nyh,100), yh1(1), ewt(1), savf(100),               
     1   acor(100), wm(100), iwm(100)                                                 
      common /ls0001/ conit, crate, el(13), elco(13,12),                        
     1   hold, rmax, tesco(3,12),                                               
     2   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround, iownd(14),              
     3   ialth, ipup, lmax, meo, nqnyh, nslp,                                   
     4   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,                
     5   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu                 
c-----------------------------------------------------------------------        
c stode performs one step of the integration of an initial value                
c problem for a system of ordinary differential equations.                      
c note.. stode is independent of the value of the iteration method              
c indicator miter, when this is .ne. 0, and hence is independent                
c of the type of chord method used, or the jacobian structure.                  
c communication with stode is done with the following variables..               
c                                                                               
c neq    = integer array containing problem size in neq(1), and                 
c          passed as the neq argument in all calls to f and jac.                
c y      = an array of length .ge. n used as the y argument in                  
c          all calls to f and jac.                                              
c yh     = an nyh by lmax array containing the dependent variables              
c          and their approximate scaled derivatives, where                      
c          lmax = maxord + 1.  yh(i,j+1) contains the approximate               
c          j-th derivative of y(i), scaled by h**j/factorial(j)                 
c          (j = 0,1,...,nq).  on entry for the first step, the first            
c          two columns of yh must be set from the initial values.               
c nyh    = a constant integer .ge. n, the first dimension of yh.                
c yh1    = a one-dimensional array occupying the same space as yh.              
c ewt    = an array of length n containing multiplicative weights               
c          for local error measurements.  local errors in y(i) are              
c          compared to 1.0/ewt(i) in various error tests.                       
c savf   = an array of working storage, of length n.                            
c          also used for input of yh(*,maxord+2) when jstart = -1               
c          and maxord .lt. the current order nq.                                
c acor   = a work array of length n, used for the accumulated                   
c          corrections.  on a successful return, acor(i) contains               
c          the estimated one-step local error in y(i).                          
c wm,iwm = real and integer work arrays associated with matrix                  
c          operations in chord iteration (miter .ne. 0).                        
c pjac   = name of routine to evaluate and preprocess jacobian matrix           
c          and p = i - h*el0*jac, if a chord method is being used.              
c slvs   = name of routine to solve linear system in chord iteration.           
c ccmax  = maximum relative change in h*el0 before pjac is called.              
c h      = the step size to be attempted on the next step.                      
c          h is altered by the error control algorithm during the               
c          problem.  h can be either positive or negative, but its              
c          sign must remain constant throughout the problem.                    
c hmin   = the minimum absolute value of the step size h to be used.            
c hmxi   = inverse of the maximum absolute value of h to be used.               
c          hmxi = 0.0 is allowed and corresponds to an infinite hmax.           
c          hmin and hmxi may be changed at any time, but will not               
c          take effect until the next change of h is considered.                
c tn     = the independent variable. tn is updated on each step taken.          
c jstart = an integer used for input only, with the following                   
c          values and meanings..                                                
c               0  perform the first step.                                      
c           .gt.0  take a new step continuing from the last.                    
c              -1  take the next step with a new value of h, maxord,            
c                    n, meth, miter, and/or matrix parameters.                  
c              -2  take the next step with a new value of h,                    
c                    but with other inputs unchanged.                           
c          on return, jstart is set to 1 to facilitate continuation.            
c kflag  = a completion code with the following meanings..                      
c               0  the step was succesful.                                      
c              -1  the requested error could not be achieved.                   
c              -2  corrector convergence could not be achieved.                 
c              -3  fatal error in pjac or slvs.                                 
c          a return with kflag = -1 or -2 means either                          
c          abs(h) = hmin or 10 consecutive failures occurred.                   
c          on a return with kflag negative, the values of tn and                
c          the yh array are as of the beginning of the last                     
c          step, and h is the last step size attempted.                         
c maxord = the maximum order of integration method to be allowed.               
c maxcor = the maximum number of corrector iterations allowed.                  
c msbp   = maximum number of steps between pjac calls (miter .gt. 0).           
c mxncf  = maximum number of convergence failures allowed.                      
c meth/miter = the method flags.  see description in driver.                    
c n      = the number of first-order differential equations.                    
c-----------------------------------------------------------------------        
      kflag = 0                                                                 
      told = tn                                                                 
      ncf = 0                                                                   
      ierpj = 0                                                                 
      iersl = 0                                                                 
      jcur = 0                                                                  
      icf = 0                                                                   
      delp = 0.0d0                                                              
      if (jstart .gt. 0) go to 200                                              
      if (jstart .eq. -1) go to 100                                             
      if (jstart .eq. -2) go to 160                                             
c-----------------------------------------------------------------------        
c on the first call, the order is set to 1, and other variables are             
c initialized.  rmax is the maximum ratio by which h can be increased           
c in a single step.  it is initially 1.e4 to compensate for the small           
c initial h, but then is normally equal to 10.  if a failure                    
c occurs (in corrector convergence or error test), rmax is set at 2             
c for the next increase.                                                        
c-----------------------------------------------------------------------        
      lmax = maxord + 1                                                         
      nq = 1                                                                    
      l = 2                                                                     
      ialth = 2                                                                 
      rmax = 10000.0d0                                                          
      rc = 0.0d0                                                                
      el0 = 1.0d0                                                               
      crate = 0.7d0                                                             
      hold = h                                                                  
      meo = meth                                                                
      nslp = 0                                                                  
      ipup = miter                                                              
      iret = 3                                                                  
      go to 140                                                                 
c-----------------------------------------------------------------------        
c the following block handles preliminaries needed when jstart = -1.            
c ipup is set to miter to force a matrix update.                                
c if an order increase is about to be considered (ialth = 1),                   
c ialth is reset to 2 to postpone consideration one more step.                  
c if the caller has changed meth, cfode is called to reset                      
c the coefficients of the method.                                               
c if the caller has changed maxord to a value less than the current             
c order nq, nq is reduced to maxord, and a new h chosen accordingly.            
c if h is to be changed, yh must be rescaled.                                   
c if h or meth is being changed, ialth is reset to l = nq + 1                   
c to prevent further changes in h for that many steps.                          
c-----------------------------------------------------------------------        
 100  ipup = miter                                                              
      lmax = maxord + 1                                                         
      if (ialth .eq. 1) ialth = 2                                               
      if (meth .eq. meo) go to 110                                              
      call cfode (meth, elco, tesco)                                            
      meo = meth                                                                
      if (nq .gt. maxord) go to 120                                             
      ialth = l                                                                 
      iret = 1                                                                  
      go to 150                                                                 
 110  if (nq .le. maxord) go to 160                                             
 120  nq = maxord                                                               
      l = lmax                                                                  
      do 125 i = 1,l                                                            
 125    el(i) = elco(i,nq)                                                      
      nqnyh = nq*nyh                                                            
      rc = rc*el(1)/el0                                                         
      el0 = el(1)                                                               
      conit = 0.5d0/dfloat(nq+2)                                                
      ddn = vnorm (n, savf, ewt)/tesco(1,l)                                     
      exdn = 1.0d0/dfloat(l)                                                    
      rhdn = 1.0d0/(1.3d0*ddn**exdn + 0.0000013d0)                              
      rh = dmin1(rhdn,1.0d0)                                                    
      iredo = 3                                                                 
      if (h .eq. hold) go to 170                                                
      rh = dmin1(rh,dabs(h/hold))                                               
      h = hold                                                                  
      go to 175                                                                 
c-----------------------------------------------------------------------        
c cfode is called to get all the integration coefficients for the               
c current meth.  then the el vector and related constants are reset             
c whenever the order nq is changed, or at the start of the problem.             
c-----------------------------------------------------------------------        
 140  call cfode (meth, elco, tesco)                                            
 150  do 155 i = 1,l                                                            
 155    el(i) = elco(i,nq)                                                      
      nqnyh = nq*nyh                                                            
      rc = rc*el(1)/el0                                                         
      el0 = el(1)                                                               
      conit = 0.5d0/dfloat(nq+2)                                                
      go to (160, 170, 200), iret                                               
c-----------------------------------------------------------------------        
c if h is being changed, the h ratio rh is checked against                      
c rmax, hmin, and hmxi, and the yh array rescaled.  ialth is set to             
c l = nq + 1 to prevent a change of h for that many steps, unless               
c forced by a convergence or error test failure.                                
c-----------------------------------------------------------------------        
 160  if (h .eq. hold) go to 200                                                
      rh = h/hold                                                               
      h = hold                                                                  
      iredo = 3                                                                 
      go to 175                                                                 
 170  rh = dmax1(rh,hmin/dabs(h))                                               
 175  rh = dmin1(rh,rmax)                                                       
      rh = rh/dmax1(1.0d0,dabs(h)*hmxi*rh)                                      
      r = 1.0d0                                                                 
      do 180 j = 2,l                                                            
        r = r*rh                                                                
        do 180 i = 1,n                                                          
 180      yh(i,j) = yh(i,j)*r                                                   
      h = h*rh                                                                  
      rc = rc*rh                                                                
      ialth = l                                                                 
      if (iredo .eq. 0) go to 690                                               
c-----------------------------------------------------------------------        
c this section computes the predicted values by effectively                     
c multiplying the yh array by the pascal triangle matrix.                       
c rc is the ratio of new to old values of the coefficient  h*el(1).             
c when rc differs from 1 by more than ccmax, ipup is set to miter               
c to force pjac to be called, if a jacobian is involved.                        
c in any case, pjac is called at least every msbp steps.                        
c-----------------------------------------------------------------------        
 200  if (dabs(rc-1.0d0) .gt. ccmax) ipup = miter                               
      if (nst .ge. nslp+msbp) ipup = miter                                      
      tn = tn + h                                                               
      i1 = nqnyh + 1                                                            
      do 215 jb = 1,nq                                                          
        i1 = i1 - nyh                                                           
cdir$ ivdep                                                                     
        do 210 i = i1,nqnyh                                                     
 210      yh1(i) = yh1(i) + yh1(i+nyh)                                          
 215    continue                                                                
c-----------------------------------------------------------------------        
c up to maxcor corrector iterations are taken.  a convergence test is           
c made on the r.m.s. norm of each correction, weighted by the error             
c weight vector ewt.  the sum of the corrections is accumulated in the          
c vector acor(i).  the yh array is not altered in the corrector loop.           
c-----------------------------------------------------------------------        
 220  m = 0                                                                     
      do 230 i = 1,n                                                            
 230    y(i) = yh(i,1)                                                          
      call f (neq, tn, y, savf)                                                 
      nfe = nfe + 1                                                             
      if (ipup .le. 0) go to 250                                                
c-----------------------------------------------------------------------        
c if indicated, the matrix p = i - h*el(1)*j is reevaluated and                 
c preprocessed before starting the corrector iteration.  ipup is set            
c to 0 as an indicator that this has been done.                                 
c-----------------------------------------------------------------------        
      call pjac (neq, y, yh, nyh, ewt, acor, savf, wm, iwm, f, jac)             
      ipup = 0                                                                  
      rc = 1.0d0                                                                
      nslp = nst                                                                
      crate = 0.7d0                                                             
      if (ierpj .ne. 0) go to 430                                               
 250  do 260 i = 1,n                                                            
 260    acor(i) = 0.0d0                                                         
 270  if (miter .ne. 0) go to 350                                               
c-----------------------------------------------------------------------        
c in the case of functional iteration, update y directly from                   
c the result of the last function evaluation.                                   
c-----------------------------------------------------------------------        
      do 290 i = 1,n                                                            
        savf(i) = h*savf(i) - yh(i,2)                                           
 290    y(i) = savf(i) - acor(i)                                                
      del = vnorm (n, y, ewt)                                                   
      do 300 i = 1,n                                                            
        y(i) = yh(i,1) + el(1)*savf(i)                                          
 300    acor(i) = savf(i)                                                       
      go to 400                                                                 
c-----------------------------------------------------------------------        
c in the case of the chord method, compute the corrector error,                 
c and solve the linear system with that as right-hand side and                  
c p as coefficient matrix.                                                      
c-----------------------------------------------------------------------        
 350  do 360 i = 1,n                                                            
 360    y(i) = h*savf(i) - (yh(i,2) + acor(i))                                  
      call slvs (wm, iwm, y, savf)                                              
      if (iersl .lt. 0) go to 430                                               
      if (iersl .gt. 0) go to 410                                               
      del = vnorm (n, y, ewt)                                                   
      do 380 i = 1,n                                                            
        acor(i) = acor(i) + y(i)                                                
 380    y(i) = yh(i,1) + el(1)*acor(i)                                          
c-----------------------------------------------------------------------        
c test for convergence.  if m.gt.0, an estimate of the convergence              
c rate constant is stored in crate, and this is used in the test.               
c-----------------------------------------------------------------------        
 400  if (m .ne. 0) crate = dmax1(0.2d0*crate,del/delp)                         
      dcon = del*dmin1(1.0d0,1.5d0*crate)/(tesco(2,nq)*conit)                   
      if (dcon .le. 1.0d0) go to 450                                            
      m = m + 1                                                                 
      if (m .eq. maxcor) go to 410                                              
      if (m .ge. 2 .and. del .gt. 2.0d0*delp) go to 410                         
      delp = del                                                                
      call f (neq, tn, y, savf)                                                 
      nfe = nfe + 1                                                             
      go to 270                                                                 
c-----------------------------------------------------------------------        
c the corrector iteration failed to converge.                                   
c if miter .ne. 0 and the jacobian is out of date, pjac is called for           
c the next try.  otherwise the yh array is retracted to its values              
c before prediction, and h is reduced, if possible.  if h cannot be             
c reduced or mxncf failures have occurred, exit with kflag = -2.                
c-----------------------------------------------------------------------        
 410  if (miter .eq. 0 .or. jcur .eq. 1) go to 430                              
      icf = 1                                                                   
      ipup = miter                                                              
      go to 220                                                                 
 430  icf = 2                                                                   
      ncf = ncf + 1                                                             
      rmax = 2.0d0                                                              
      tn = told                                                                 
      i1 = nqnyh + 1                                                            
      do 445 jb = 1,nq                                                          
        i1 = i1 - nyh                                                           
cdir$ ivdep                                                                     
        do 440 i = i1,nqnyh                                                     
 440      yh1(i) = yh1(i) - yh1(i+nyh)                                          
 445    continue                                                                
      if (ierpj .lt. 0 .or. iersl .lt. 0) go to 680                             
      if (dabs(h) .le. hmin*1.00001d0) go to 670                                
      if (ncf .eq. mxncf) go to 670                                             
      rh = 0.25d0                                                               
      ipup = miter                                                              
      iredo = 1                                                                 
      go to 170                                                                 
c-----------------------------------------------------------------------        
c the corrector has converged.  jcur is set to 0                                
c to signal that the jacobian involved may need updating later.                 
c the local error test is made and control passes to statement 500              
c if it fails.                                                                  
c-----------------------------------------------------------------------        
 450  jcur = 0                                                                  
      if (m .eq. 0) dsm = del/tesco(2,nq)                                       
      if (m .gt. 0) dsm = vnorm (n, acor, ewt)/tesco(2,nq)                      
      if (dsm .gt. 1.0d0) go to 500                                             
c-----------------------------------------------------------------------        
c after a successful step, update the yh array.                                 
c consider changing h if ialth = 1.  otherwise decrease ialth by 1.             
c if ialth is then 1 and nq .lt. maxord, then acor is saved for                 
c use in a possible order increase on the next step.                            
c if a change in h is considered, an increase or decrease in order              
c by one is considered also.  a change in h is made only if it is by a          
c factor of at least 1.1.  if not, ialth is set to 3 to prevent                 
c testing for that many steps.                                                  
c-----------------------------------------------------------------------        
      kflag = 0                                                                 
      iredo = 0                                                                 
      nst = nst + 1                                                             
      hu = h                                                                    
      nqu = nq                                                                  
      do 470 j = 1,l                                                            
        do 470 i = 1,n                                                          
 470      yh(i,j) = yh(i,j) + el(j)*acor(i)                                     
      ialth = ialth - 1                                                         
      if (ialth .eq. 0) go to 520                                               
      if (ialth .gt. 1) go to 700                                               
      if (l .eq. lmax) go to 700                                                
      do 490 i = 1,n                                                            
 490    yh(i,lmax) = acor(i)                                                    
      go to 700                                                                 
c-----------------------------------------------------------------------        
c the error test failed.  kflag keeps track of multiple failures.               
c restore tn and the yh array to their previous values, and prepare             
c to try the step again.  compute the optimum step size for this or             
c one lower order.  after 2 or more failures, h is forced to decrease           
c by a factor of 0.2 or less.                                                   
c-----------------------------------------------------------------------        
 500  kflag = kflag - 1                                                         
      tn = told                                                                 
      i1 = nqnyh + 1                                                            
      do 515 jb = 1,nq                                                          
        i1 = i1 - nyh                                                           
cdir$ ivdep                                                                     
        do 510 i = i1,nqnyh                                                     
 510      yh1(i) = yh1(i) - yh1(i+nyh)                                          
 515    continue                                                                
      rmax = 2.0d0                                                              
      if (dabs(h) .le. hmin*1.00001d0) go to 660                                
      if (kflag .le. -3) go to 640                                              
      iredo = 2                                                                 
      rhup = 0.0d0                                                              
      go to 540                                                                 
c-----------------------------------------------------------------------        
c regardless of the success or failure of the step, factors                     
c rhdn, rhsm, and rhup are computed, by which h could be multiplied             
c at order nq - 1, order nq, or order nq + 1, respectively.                     
c in the case of failure, rhup = 0.0 to avoid an order increase.                
c the largest of these is determined and the new order chosen                   
c accordingly.  if the order is to be increased, we compute one                 
c additional scaled derivative.                                                 
c-----------------------------------------------------------------------        
 520  rhup = 0.0d0                                                              
      if (l .eq. lmax) go to 540                                                
      do 530 i = 1,n                                                            
 530    savf(i) = acor(i) - yh(i,lmax)                                          
      dup = vnorm (n, savf, ewt)/tesco(3,nq)                                    
      exup = 1.0d0/dfloat(l+1)                                                  
      rhup = 1.0d0/(1.4d0*dup**exup + 0.0000014d0)                              
 540  exsm = 1.0d0/dfloat(l)                                                    
      rhsm = 1.0d0/(1.2d0*dsm**exsm + 0.0000012d0)                              
      rhdn = 0.0d0                                                              
      if (nq .eq. 1) go to 560                                                  
      ddn = vnorm (n, yh(1,l), ewt)/tesco(1,nq)                                 
      exdn = 1.0d0/dfloat(nq)                                                   
      rhdn = 1.0d0/(1.3d0*ddn**exdn + 0.0000013d0)                              
 560  if (rhsm .ge. rhup) go to 570                                             
      if (rhup .gt. rhdn) go to 590                                             
      go to 580                                                                 
 570  if (rhsm .lt. rhdn) go to 580                                             
      newq = nq                                                                 
      rh = rhsm                                                                 
      go to 620                                                                 
 580  newq = nq - 1                                                             
      rh = rhdn                                                                 
      if (kflag .lt. 0 .and. rh .gt. 1.0d0) rh = 1.0d0                          
      go to 620                                                                 
 590  newq = l                                                                  
      rh = rhup                                                                 
      if (rh .lt. 1.1d0) go to 610                                              
      r = el(l)/dfloat(l)                                                       
      do 600 i = 1,n                                                            
 600    yh(i,newq+1) = acor(i)*r                                                
      go to 630                                                                 
 610  ialth = 3                                                                 
      go to 700                                                                 
 620  if ((kflag .eq. 0) .and. (rh .lt. 1.1d0)) go to 610                       
      if (kflag .le. -2) rh = dmin1(rh,0.2d0)                                   
c-----------------------------------------------------------------------        
c if there is a change of order, reset nq, l, and the coefficients.             
c in any case h is reset according to rh and the yh array is rescaled.          
c then exit from 690 if the step was ok, or redo the step otherwise.            
c-----------------------------------------------------------------------        
      if (newq .eq. nq) go to 170                                               
 630  nq = newq                                                                 
      l = nq + 1                                                                
      iret = 2                                                                  
      go to 150                                                                 
c-----------------------------------------------------------------------        
c control reaches this section if 3 or more failures have occured.              
c if 10 failures have occurred, exit with kflag = -1.                           
c it is assumed that the derivatives that have accumulated in the               
c yh array have errors of the wrong order.  hence the first                     
c derivative is recomputed, and the order is set to 1.  then                    
c h is reduced by a factor of 10, and the step is retried,                      
c until it succeeds or h reaches hmin.                                          
c-----------------------------------------------------------------------        
 640  if (kflag .eq. -10) go to 660                                             
      rh = 0.1d0                                                                
      rh = dmax1(hmin/dabs(h),rh)                                               
      h = h*rh                                                                  
      do 645 i = 1,n                                                            
 645    y(i) = yh(i,1)                                                          
      call f (neq, tn, y, savf)                                                 
      nfe = nfe + 1                                                             
      do 650 i = 1,n                                                            
 650    yh(i,2) = h*savf(i)                                                     
      ipup = miter                                                              
      ialth = 5                                                                 
      if (nq .eq. 1) go to 200                                                  
      nq = 1                                                                    
      l = 2                                                                     
      iret = 3                                                                  
      go to 150                                                                 
c-----------------------------------------------------------------------        
c all returns are made through this section.  h is saved in hold                
c to allow the caller to change h on the next step.                             
c-----------------------------------------------------------------------        
 660  kflag = -1                                                                
      go to 720                                                                 
 670  kflag = -2                                                                
      go to 720                                                                 
 680  kflag = -3                                                                
      go to 720                                                                 
 690  rmax = 10.0d0                                                             
 700  r = 1.0d0/tesco(2,nqu)                                                    
      do 710 i = 1,n                                                            
 710    acor(i) = acor(i)*r                                                     
 720  hold = h                                                                  
      jstart = 1                                                                
      return                                                                    
c----------------------- end of subroutine stode -----------------------        
      end                                                                       
      double precision function vnorm (n, v, w)                                 
clll. optimize                                                                  
c-----------------------------------------------------------------------        
c this function routine computes the weighted root-mean-square norm             
c of the vector of length n contained in the array v, with weights              
c contained in the array w of length n..                                        
c   vnorm = sqrt( (1/n) * sum( v(i)*w(i) )**2 )                                 
c-----------------------------------------------------------------------        
      integer n,   i                                                            
      double precision v, w,   sum                                              
      dimension v(n), w(n)                                                      
      sum = 0.0d0                                                               
      do 10 i = 1,n                                                             
 10     sum = sum + (v(i)*w(i))**2                                              
      vnorm = dsqrt(sum/dfloat(n))                                              
      return                                                                    
c----------------------- end of function vnorm -------------------------        
      end                                                                       
      subroutine cfode (meth, elco, tesco)                                      
clll. optimize                                                                  
      integer meth                                                              
      integer i, ib, nq, nqm1, nqp1                                             
      double precision elco, tesco                                              
      double precision agamq, fnq, fnqm1, pc, pint, ragq,                       
     1   rqfac, rq1fac, tsign, xpin                                             
      dimension elco(13,12), tesco(3,12)                                        
c-----------------------------------------------------------------------        
c cfode is called by the integrator routine to set coefficients                 
c needed there.  the coefficients for the current method, as                    
c given by the value of meth, are set for all orders and saved.                 
c the maximum order assumed here is 12 if meth = 1 and 5 if meth = 2.           
c (a smaller value of the maximum order is also allowed.)                       
c cfode is called once at the beginning of the problem,                         
c and is not called again unless and until meth is changed.                     
c                                                                               
c the elco array contains the basic method coefficients.                        
c the coefficients el(i), 1 .le. i .le. nq+1, for the method of                 
c order nq are stored in elco(i,nq).  they are given by a genetrating           
c polynomial, i.e.,                                                             
c     l(x) = el(1) + el(2)*x + ... + el(nq+1)*x**nq.                            
c for the implicit adams methods, l(x) is given by                              
c     dl/dx = (x+1)*(x+2)*...*(x+nq-1)/factorial(nq-1),    l(-1) = 0.           
c for the bdf methods, l(x) is given by                                         
c     l(x) = (x+1)*(x+2)* ... *(x+nq)/k,                                        
c where         k = factorial(nq)*(1 + 1/2 + ... + 1/nq).                       
c                                                                               
c the tesco array contains test constants used for the                          
c local error test and the selection of step size and/or order.                 
c at order nq, tesco(k,nq) is used for the selection of step                    
c size at order nq - 1 if k = 1, at order nq if k = 2, and at order             
c nq + 1 if k = 3.                                                              
c-----------------------------------------------------------------------        
      dimension pc(12)                                                          
c                                                                               
      go to (100, 200), meth                                                    
c                                                                               
 100  elco(1,1) = 1.0d0                                                         
      elco(2,1) = 1.0d0                                                         
      tesco(1,1) = 0.0d0                                                        
      tesco(2,1) = 2.0d0                                                        
      tesco(1,2) = 1.0d0                                                        
      tesco(3,12) = 0.0d0                                                       
      pc(1) = 1.0d0                                                             
      rqfac = 1.0d0                                                             
      do 140 nq = 2,12                                                          
c-----------------------------------------------------------------------        
c the pc array will contain the coefficients of the polynomial                  
c     p(x) = (x+1)*(x+2)*...*(x+nq-1).                                          
c initially, p(x) = 1.                                                          
c-----------------------------------------------------------------------        
        rq1fac = rqfac                                                          
        rqfac = rqfac/dfloat(nq)                                                
        nqm1 = nq - 1                                                           
        fnqm1 = dfloat(nqm1)                                                    
        nqp1 = nq + 1                                                           
c form coefficients of p(x)*(x+nq-1). ----------------------------------        
        pc(nq) = 0.0d0                                                          
        do 110 ib = 1,nqm1                                                      
          i = nqp1 - ib                                                         
 110      pc(i) = pc(i-1) + fnqm1*pc(i)                                         
        pc(1) = fnqm1*pc(1)                                                     
c compute integral, -1 to 0, of p(x) and x*p(x). -----------------------        
        pint = pc(1)                                                            
        xpin = pc(1)/2.0d0                                                      
        tsign = 1.0d0                                                           
        do 120 i = 2,nq                                                         
          tsign = -tsign                                                        
          pint = pint + tsign*pc(i)/dfloat(i)                                   
 120      xpin = xpin + tsign*pc(i)/dfloat(i+1)                                 
c store coefficients in elco and tesco. --------------------------------        
        elco(1,nq) = pint*rq1fac                                                
        elco(2,nq) = 1.0d0                                                      
        do 130 i = 2,nq                                                         
 130      elco(i+1,nq) = rq1fac*pc(i)/dfloat(i)                                 
        agamq = rqfac*xpin                                                      
        ragq = 1.0d0/agamq                                                      
        tesco(2,nq) = ragq                                                      
        if (nq .lt. 12) tesco(1,nqp1) = ragq*rqfac/dfloat(nqp1)                 
        tesco(3,nqm1) = ragq                                                    
 140    continue                                                                
      return                                                                    
c                                                                               
 200  pc(1) = 1.0d0                                                             
      rq1fac = 1.0d0                                                            
      do 230 nq = 1,5                                                           
c-----------------------------------------------------------------------        
c the pc array will contain the coefficients of the polynomial                  
c     p(x) = (x+1)*(x+2)*...*(x+nq).                                            
c initially, p(x) = 1.                                                          
c-----------------------------------------------------------------------        
        fnq = dfloat(nq)                                                        
        nqp1 = nq + 1                                                           
c form coefficients of p(x)*(x+nq). ------------------------------------        
        pc(nqp1) = 0.0d0                                                        
        do 210 ib = 1,nq                                                        
          i = nq + 2 - ib                                                       
 210      pc(i) = pc(i-1) + fnq*pc(i)                                           
        pc(1) = fnq*pc(1)                                                       
c store coefficients in elco and tesco. --------------------------------        
        do 220 i = 1,nqp1                                                       
 220      elco(i,nq) = pc(i)/pc(2)                                              
        elco(2,nq) = 1.0d0                                                      
        tesco(1,nq) = rq1fac                                                    
        tesco(2,nq) = dfloat(nqp1)/elco(1,nq)                                   
        tesco(3,nq) = dfloat(nq+2)/elco(1,nq)                                   
        rq1fac = rq1fac/fnq                                                     
 230    continue                                                                
      return                                                                    
c----------------------- end of subroutine cfode -----------------------        
      end                                                                       
      subroutine intdy (t, k, yh, nyh, dky, iflag)                              
clll. optimize                                                                  
      integer k, nyh, iflag                                                     
      integer iownd, iowns,                                                     
     1   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,                
     2   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu                 
      integer i, ic, j, jb, jb2, jj, jj1, jp1                                   
      double precision t, yh, dky                                               
      double precision rowns,                                                   
     1   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround                          
      double precision c, r, s, tp                                              
      dimension yh(nyh,1), dky(1)                                               
      common /ls0001/ rowns(209),                                               
     2   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround,                         
     3   iownd(14), iowns(6),                                                   
     4   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,                
     5   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu                 
c-----------------------------------------------------------------------        
c intdy computes interpolated values of the k-th derivative of the              
c dependent variable vector y, and stores it in dky.  this routine              
c is called within the package with k = 0 and t = tout, but may                 
c also be called by the user for any k up to the current order.                 
c (see detailed instructions in the usage documentation.)                       
c-----------------------------------------------------------------------        
c the computed values in dky are gotten by interpolation using the              
c nordsieck history array yh.  this array corresponds uniquely to a             
c vector-valued polynomial of degree nqcur or less, and dky is set              
c to the k-th derivative of this polynomial at t.                               
c the formula for dky is..                                                      
c              q                                                                
c  dky(i)  =  sum  c(j,k) * (t - tn)**(j-k) * h**(-j) * yh(i,j+1)               
c             j=k                                                               
c where  c(j,k) = j*(j-1)*...*(j-k+1), q = nqcur, tn = tcur, h = hcur.          
c the quantities  nq = nqcur, l = nq+1, n = neq, tn, and h are                  
c communicated by common.  the above sum is done in reverse order.              
c iflag is returned negative if either k or t is out of bounds.                 
c-----------------------------------------------------------------------        
      iflag = 0                                                                 
      if (k .lt. 0 .or. k .gt. nq) go to 80                                     
      tp = tn - hu -  100.0d0*uround*(tn + hu)                                  
      if ((t-tp)*(t-tn) .gt. 0.0d0) go to 90                                    
c                                                                               
      s = (t - tn)/h                                                            
      ic = 1                                                                    
      if (k .eq. 0) go to 15                                                    
      jj1 = l - k                                                               
      do 10 jj = jj1,nq                                                         
 10     ic = ic*jj                                                              
 15   c = dfloat(ic)                                                            
      do 20 i = 1,n                                                             
 20     dky(i) = c*yh(i,l)                                                      
      if (k .eq. nq) go to 55                                                   
      jb2 = nq - k                                                              
      do 50 jb = 1,jb2                                                          
        j = nq - jb                                                             
        jp1 = j + 1                                                             
        ic = 1                                                                  
        if (k .eq. 0) go to 35                                                  
        jj1 = jp1 - k                                                           
        do 30 jj = jj1,j                                                        
 30       ic = ic*jj                                                            
 35     c = dfloat(ic)                                                          
        do 40 i = 1,n                                                           
 40       dky(i) = c*yh(i,jp1) + s*dky(i)                                       
 50     continue                                                                
      if (k .eq. 0) return                                                      
 55   r = h**(-k)                                                               
      do 60 i = 1,n                                                             
 60     dky(i) = r*dky(i)                                                       
      return                                                                    
c                                                                               
 80   call xerrwv(30hintdy--  k (=i1) illegal      ,                            
     1   30, 51, 0, 1, k, 0, 0, 0.0d0, 0.0d0)                                   
      iflag = -1                                                                
      return                                                                    
 90   call xerrwv(30hintdy--  t (=r1) illegal      ,                            
     1   30, 52, 0, 0, 0, 0, 1, t, 0.0d0)                                       
      call xerrwv(                                                              
     1  60h      t not in interval tcur - hu (= r1) to tcur (=r2)      ,        
     1   60, 52, 0, 0, 0, 0, 2, tp, tn)                                         
      iflag = -2                                                                
      return                                                                    
c----------------------- end of subroutine intdy -----------------------        
      end                                                                       
      subroutine xerrwv (msg, nmes, nerr, level, ni, i1, i2, nr, r1, r2)        
      integer msg, nmes, nerr, level, ni, i1, i2, nr,                           
     1   i, lun, lunit, mesflg, ncpw, nch, nwds                                 
      double precision r1, r2                                                   
      dimension msg(nmes)                                                       
c-----------------------------------------------------------------------        
c subroutines xerrwv, xsetf, and xsetun, as given here, constitute              
c a simplified version of the slatec error handling package.                    
c written by a. c. hindmarsh at llnl.  version of march 30, 1987.               
c this version is in double precision.                                          
c                                                                               
c all arguments are input arguments.                                            
c                                                                               
c msg    = the message (hollerith literal or integer array).                    
c nmes   = the length of msg (number of characters).                            
c nerr   = the error number (not used).                                         
c level  = the error level..                                                    
c          0 or 1 means recoverable (control returns to caller).                
c          2 means fatal (run is aborted--see note below).                      
c ni     = number of integers (0, 1, or 2) to be printed with message.          
c i1,i2  = integers to be printed, depending on ni.                             
c nr     = number of reals (0, 1, or 2) to be printed with message.             
c r1,r2  = reals to be printed, depending on nr.                                
c                                                                               
c note..  this routine is machine-dependent and specialized for use             
c in limited context, in the following ways..                                   
c 1. the number of hollerith characters stored per word, denoted                
c    by ncpw below, is a data-loaded constant.                                  
c 2. the value of nmes is assumed to be at most 60.                             
c    (multi-line messages are generated by repeated calls.)                     
c 3. if level = 2, control passes to the statement   stop                       
c    to abort the run.  this statement may be machine-dependent.                
c 4. r1 and r2 are assumed to be in double precision and are printed            
c    in d21.13 format.                                                          
c 5. the common block /eh0001/ below is data-loaded (a machine-                 
c    dependent feature) with default values.                                    
c    this block is needed for proper retention of parameters used by            
c    this routine which the user can reset by calling xsetf or xsetun.          
c    the variables in this block are as follows..                               
c       mesflg = print control flag..                                           
c                1 means print all messages (the default).                      
c                0 means no printing.                                           
c       lunit  = logical unit number for messages.                              
c                the default is 6 (machine-dependent).                          
c-----------------------------------------------------------------------        
c the following are instructions for installing this routine                    
c in different machine environments.                                            
c                                                                               
c to change the default output unit, change the data statement                  
c in the block data subprogram below.                                           
c                                                                               
c for a different number of characters per word, change the                     
c data statement setting ncpw below, and format 10.  alternatives for           
c various computers are shown in comment cards.                                 
c                                                                               
c for a different run-abort command, change the statement following             
c statement 100 at the end.                                                     
c-----------------------------------------------------------------------        
      common /eh0001/ mesflg, lunit                                             
c-----------------------------------------------------------------------        
c the following data-loaded value of ncpw is valid for the cdc-6600             
c and cdc-7600 computers.                                                       
c     data ncpw/10/                                                             
c the following is valid for the cray-1 computer.                               
c     data ncpw/8/                                                              
c the following is valid for the burroughs 6700 and 7800 computers.             
c     data ncpw/6/                                                              
c the following is valid for the pdp-10 computer.                               
c     data ncpw/5/                                                              
c the following is valid for the vax computer with 4 bytes per integer,         
c and for the ibm-360, ibm-370, ibm-303x, and ibm-43xx computers.               
      data ncpw/4/                                                              
c the following is valid for the pdp-11, or vax with 2-byte integers.           
c     data ncpw/2/                                                              
c-----------------------------------------------------------------------        
      if (mesflg .eq. 0) go to 100                                              
c get logical unit number. ---------------------------------------------        
      lun = lunit                                                               
c get number of words in message. --------------------------------------        
      nch = min0(nmes,60)                                                       
      nwds = nch/ncpw                                                           
      if (nch .ne. nwds*ncpw) nwds = nwds + 1                                   
c write the message. ---------------------------------------------------        
      write (lun, 10) (msg(i),i=1,nwds)                                         
c-----------------------------------------------------------------------        
c the following format statement is to have the form                            
c 10  format(1x,mmann)                                                          
c where nn = ncpw and mm is the smallest integer .ge. 60/ncpw.                  
c the following is valid for ncpw = 10.                                         
c 10  format(1x,6a10)                                                           
c the following is valid for ncpw = 8.                                          
c 10  format(1x,8a8)                                                            
c the following is valid for ncpw = 6.                                          
c 10  format(1x,10a6)                                                           
c the following is valid for ncpw = 5.                                          
c 10  format(1x,12a5)                                                           
c the following is valid for ncpw = 4.                                          
  10  format(1x,15a4)                                                           
c the following is valid for ncpw = 2.                                          
c 10  format(1x,30a2)                                                           
c-----------------------------------------------------------------------        
      if (ni .eq. 1) write (lun, 20) i1                                         
 20   format(6x,23hin above message,  i1 =,i10)                                 
      if (ni .eq. 2) write (lun, 30) i1,i2                                      
 30   format(6x,23hin above message,  i1 =,i10,3x,4hi2 =,i10)                   
      if (nr .eq. 1) write (lun, 40) r1                                         
 40   format(6x,23hin above message,  r1 =,d21.13)                              
      if (nr .eq. 2) write (lun, 50) r1,r2                                      
 50   format(6x,15hin above,  r1 =,d21.13,3x,4hr2 =,d21.13)                     
c abort the run if level = 2. ------------------------------------------        
 100  if (level .ne. 2) return                                                  
      stop                                                                      
c----------------------- end of subroutine xerrwv ----------------------        
      end                                                                       
      block data                                                                
c-----------------------------------------------------------------------        
c this data subprogram loads variables into the internal common                 
c blocks used by the odepack solvers.  the variables are                        
c defined as follows..                                                          
c   illin  = counter for the number of consecutive times the package            
c            was called with illegal input.  the run is stopped when            
c            illin reaches 5.                                                   
c   ntrep  = counter for the number of consecutive times the package            
c            was called with istate = 1 and tout = t.  the run is               
c            stopped when ntrep reaches 5.                                      
c   mesflg = flag to control printing of error messages.  1 means print,        
c            0 means no printing.                                               
c   lunit  = default value of logical unit number for printing of error         
c            messages.                                                          
c-----------------------------------------------------------------------        
      integer illin, iduma, ntrep, idumb, iowns, icomm, mesflg, lunit           
      double precision rowns, rcomm                                             
      common /ls0001/ rowns(209), rcomm(9),                                     
     1   illin, iduma(10), ntrep, idumb(2), iowns(6), icomm(19)                 
      common /eh0001/ mesflg, lunit                                             
      data illin/0/, ntrep/0/                                                   
      data mesflg/1/, lunit/6/                                                  
c                                                                               
c----------------------- end of block data -----------------------------        
      end                                                                       
      subroutine ewset (n, itol, rtol, atol, ycur, ewt)                         
clll. optimize                                                                  
c-----------------------------------------------------------------------        
c this subroutine sets the error weight vector ewt according to                 
c     ewt(i) = rtol(i)*abs(ycur(i)) + atol(i),  i = 1,...,n,                    
c with the subscript on rtol and/or atol possibly replaced by 1 above,          
c depending on the value of itol.                                               
c-----------------------------------------------------------------------        
      integer n, itol                                                           
      integer i                                                                 
      double precision rtol, atol, ycur, ewt                                    
      dimension rtol(1), atol(1), ycur(n), ewt(n)                               
c                                                                               
      go to (10, 20, 30, 40), itol                                              
 10   continue                                                                  
      do 15 i = 1,n                                                             
 15     ewt(i) = rtol(1)*dabs(ycur(i)) + atol(1)                                
      return                                                                    
 20   continue                                                                  
      do 25 i = 1,n                                                             
 25     ewt(i) = rtol(1)*dabs(ycur(i)) + atol(i)                                
      return                                                                    
 30   continue                                                                  
      do 35 i = 1,n                                                             
 35     ewt(i) = rtol(i)*dabs(ycur(i)) + atol(1)                                
      return                                                                    
 40   continue                                                                  
      do 45 i = 1,n                                                             
 45     ewt(i) = rtol(i)*dabs(ycur(i)) + atol(i)
      return
c----------------------- end of subroutine ewset -----------------------
      end
c
      double precision function dzasum(n,zx,incx)
c
c     takes the sum of the absolute values.
c     jack dongarra, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double complex zx(*)
      double precision stemp,dcabs1
      integer i,incx,ix,n
c
      dzasum = 0.0d0
      stemp = 0.0d0
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      do 10 i = 1,n
        stemp = stemp + dcabs1(zx(ix))
        ix = ix + incx
   10 continue
      dzasum = stemp
      return
c
c        code for increment equal to 1
c
   20 do 30 i = 1,n
        stemp = stemp + dcabs1(zx(i))
   30 continue
      dzasum = stemp
      return
      end
      integer function izamax(n,zx,incx)
c
c     finds the index of element having max. absolute value.
c     jack dongarra, 1/15/85.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double complex zx(*)
      double precision smax
      integer i,incx,ix,n
      double precision dcabs1
c
      izamax = 0
      if( n.lt.1 .or. incx.le.0 )return
      izamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      smax = dcabs1(zx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(dcabs1(zx(ix)).le.smax) go to 5
         izamax = i
         smax = dcabs1(zx(ix))
    5    ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 smax = dcabs1(zx(1))
      do 30 i = 2,n
         if(dcabs1(zx(i)).le.smax) go to 30
         izamax = i
         smax = dcabs1(zx(i))
   30 continue
      return
      end
      subroutine zaxpy(n,za,zx,incx,zy,incy)
c
c     constant times a vector plus a vector.
c     jack dongarra, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double complex zx(*),zy(*),za
      integer i,incx,incy,ix,iy,n
      double precision dcabs1
      if(n.le.0)return
      if (dcabs1(za) .eq. 0.0d0) return
      if (incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        zy(iy) = zy(iy) + za*zx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
   20 do 30 i = 1,n
        zy(i) = zy(i) + za*zx(i)
   30 continue
      return
      end
      double complex function zdotc(n,zx,incx,zy,incy)
c
c     forms the dot product of a vector.
c     jack dongarra, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double complex zx(*),zy(*),ztemp
      integer i,incx,incy,ix,iy,n
      ztemp = (0.0d0,0.0d0)
      zdotc = (0.0d0,0.0d0)
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        ztemp = ztemp + dconjg(zx(ix))*zy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      zdotc = ztemp
      return
c
c        code for both increments equal to 1
c
   20 do 30 i = 1,n
        ztemp = ztemp + dconjg(zx(i))*zy(i)
   30 continue
      zdotc = ztemp
      return
      end
      subroutine  zdscal(n,da,zx,incx)
c
c     scales a vector by a constant.
c     jack dongarra, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double complex zx(*)
      double precision da
      integer i,incx,ix,n
c
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      do 10 i = 1,n
        zx(ix) = dcmplx(da,0.0d0)*zx(ix)
        ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 do 30 i = 1,n
        zx(i) = dcmplx(da,0.0d0)*zx(i)
   30 continue
      return
      end
      subroutine  zscal(n,za,zx,incx)
c
c     scales a vector by a constant.
c     jack dongarra, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double complex za,zx(*)
      integer i,incx,ix,n
c
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      do 10 i = 1,n
        zx(ix) = za*zx(ix)
        ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 do 30 i = 1,n
        zx(i) = za*zx(i)
   30 continue
      return
      end
      double precision function dcabs1(z)
      double complex z,zz
      double precision t(2)
      equivalence (zz,t(1))
      zz = z
      dcabs1 = dabs(t(1)) + dabs(t(2))
      return
      end
      subroutine zgeco(a,lda,n,ipvt,rcond,z)
      integer lda,n,ipvt(1)
      complex*16 a(lda,1),z(1)
      double precision rcond
c
c     zgeco factors a complex*16 matrix by gaussian elimination
c     and estimates the condition of the matrix.
c
c     if  rcond  is not needed, zgefa is slightly faster.
c     to solve  a*x = b , follow zgeco by zgesl.
c     to compute  inverse(a)*c , follow zgeco by zgesl.
c     to compute  determinant(a) , follow zgeco by zgedi.
c     to compute  inverse(a) , follow zgeco by zgedi.
c
c     on entry
c
c        a       complex*16(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        rcond   double precision
c                an estimate of the reciprocal condition of  a .
c                for the system  a*x = b , relative perturbations
c                in  a  and  b  of size  epsilon  may cause
c                relative perturbations in  x  of size  epsilon/rcond .
c                if  rcond  is so small that the logical expression
c                           1.0 + rcond .eq. 1.0
c                is true, then  a  may be singular to working
c                precision.  in particular,  rcond  is zero  if
c                exact singularity is detected or the estimate
c                underflows.
c
c        z       complex*16(n)
c                a work vector whose contents are usually unimportant.
c                if  a  is close to a singular matrix, then  z  is
c                an approximate null vector in the sense that
c                norm(a*z) = rcond*norm(a)*norm(z) .
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     linpack zgefa
c     blas zaxpy,zdotc,zdscal,dzasum
c     fortran dabs,dmax1,dcmplx,dconjg
c
c     internal variables
c
      complex*16 zdotc,ek,t,wk,wkm
      double precision anorm,s,dzasum,sm,ynorm
      integer info,j,k,kb,kp1,l
c
      complex*16 zdum,zdum1,zdum2,csign1
      double precision cabs1
      double precision dreal,dimag
      complex*16 zdumr,zdumi
      dreal(zdumr) = zdumr
      dimag(zdumi) = (0.0d0,-1.0d0)*zdumi
      cabs1(zdum) = dabs(dreal(zdum)) + dabs(dimag(zdum))
      csign1(zdum1,zdum2) = cabs1(zdum1)*(zdum2/cabs1(zdum2))
c
c     compute 1-norm of a
c
      anorm = 0.0d0
      do 10 j = 1, n
         anorm = dmax1(anorm,dzasum(n,a(1,j),1))
   10 continue
c
c     factor
c
      call zgefa(a,lda,n,ipvt,info)
c
c     rcond = 1/(norm(a)*(estimate of norm(inverse(a)))) .
c     estimate = norm(z)/norm(y) where  a*z = y  and  ctrans(a)*y = e .
c     ctrans(a)  is the conjugate transpose of a .
c     the components of  e  are chosen to cause maximum local
c     growth in the elements of w  where  ctrans(u)*w = e .
c     the vectors are frequently rescaled to avoid overflow.
c
c     solve ctrans(u)*w = e
c
      ek = (1.0d0,0.0d0)
      do 20 j = 1, n
         z(j) = (0.0d0,0.0d0)
   20 continue
      do 100 k = 1, n
         if (cabs1(z(k)) .ne. 0.0d0) ek = csign1(ek,-z(k))
         if (cabs1(ek-z(k)) .le. cabs1(a(k,k))) go to 30
            s = cabs1(a(k,k))/cabs1(ek-z(k))
            call zdscal(n,s,z,1)
            ek = dcmplx(s,0.0d0)*ek
   30    continue
         wk = ek - z(k)
         wkm = -ek - z(k)
         s = cabs1(wk)
         sm = cabs1(wkm)
         if (cabs1(a(k,k)) .eq. 0.0d0) go to 40
            wk = wk/dconjg(a(k,k))
            wkm = wkm/dconjg(a(k,k))
         go to 50
   40    continue
            wk = (1.0d0,0.0d0)
            wkm = (1.0d0,0.0d0)
   50    continue
         kp1 = k + 1
         if (kp1 .gt. n) go to 90
            do 60 j = kp1, n
               sm = sm + cabs1(z(j)+wkm*dconjg(a(k,j)))
               z(j) = z(j) + wk*dconjg(a(k,j))
               s = s + cabs1(z(j))
   60       continue
            if (s .ge. sm) go to 80
               t = wkm - wk
               wk = wkm
               do 70 j = kp1, n
                  z(j) = z(j) + t*dconjg(a(k,j))
   70          continue
   80       continue
   90    continue
         z(k) = wk
  100 continue
      s = 1.0d0/dzasum(n,z,1)
      call zdscal(n,s,z,1)
c
c     solve ctrans(l)*y = w
c
      do 120 kb = 1, n
         k = n + 1 - kb
         if (k .lt. n) z(k) = z(k) + zdotc(n-k,a(k+1,k),1,z(k+1),1)
         if (cabs1(z(k)) .le. 1.0d0) go to 110
            s = 1.0d0/cabs1(z(k))
            call zdscal(n,s,z,1)
  110    continue
         l = ipvt(k)
         t = z(l)
         z(l) = z(k)
         z(k) = t
  120 continue
      s = 1.0d0/dzasum(n,z,1)
      call zdscal(n,s,z,1)
c
      ynorm = 1.0d0
c
c     solve l*v = y
c
      do 140 k = 1, n
         l = ipvt(k)
         t = z(l)
         z(l) = z(k)
         z(k) = t
         if (k .lt. n) call zaxpy(n-k,t,a(k+1,k),1,z(k+1),1)
         if (cabs1(z(k)) .le. 1.0d0) go to 130
            s = 1.0d0/cabs1(z(k))
            call zdscal(n,s,z,1)
            ynorm = s*ynorm
  130    continue
  140 continue
      s = 1.0d0/dzasum(n,z,1)
      call zdscal(n,s,z,1)
      ynorm = s*ynorm
c
c     solve  u*z = v
c
      do 160 kb = 1, n
         k = n + 1 - kb
         if (cabs1(z(k)) .le. cabs1(a(k,k))) go to 150
            s = cabs1(a(k,k))/cabs1(z(k))
            call zdscal(n,s,z,1)
            ynorm = s*ynorm
  150    continue
         if (cabs1(a(k,k)) .ne. 0.0d0) z(k) = z(k)/a(k,k)
         if (cabs1(a(k,k)) .eq. 0.0d0) z(k) = (1.0d0,0.0d0)
         t = -z(k)
         call zaxpy(k-1,t,a(1,k),1,z(1),1)
  160 continue
c     make znorm = 1.0
      s = 1.0d0/dzasum(n,z,1)
      call zdscal(n,s,z,1)
      ynorm = s*ynorm
c
      if (anorm .ne. 0.0d0) rcond = ynorm/anorm
      if (anorm .eq. 0.0d0) rcond = 0.0d0
      return
      end
      subroutine zgefa(a,lda,n,ipvt,info)
      integer lda,n,ipvt(1),info
      complex*16 a(lda,1)
c
c     zgefa factors a complex*16 matrix by gaussian elimination.
c
c     zgefa is usually called by zgeco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for zgeco) = (1 + 9/n)*(time for zgefa) .
c
c     on entry
c
c        a       complex*16(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that zgesl or zgedi will divide by zero
c                     if called.  use  rcond  in zgeco for a reliable
c                     indication of singularity.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas zaxpy,zscal,izamax
c     fortran dabs
c
c     internal variables
c
      complex*16 t
      integer izamax,j,k,kp1,l,nm1
c
      complex*16 zdum
      double precision cabs1
      double precision dreal,dimag
      complex*16 zdumr,zdumi
      dreal(zdumr) = zdumr
      dimag(zdumi) = (0.0d0,-1.0d0)*zdumi
      cabs1(zdum) = dabs(dreal(zdum)) + dabs(dimag(zdum))
c
c     gaussian elimination with partial pivoting
c
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
c
c        find l = pivot index
c
         l = izamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
c
c        zero pivot implies this column already triangularized
c
         if (cabs1(a(l,k)) .eq. 0.0d0) go to 40
c
c           interchange if necessary
c
            if (l .eq. k) go to 10
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
   10       continue
c
c           compute multipliers
c
            t = -(1.0d0,0.0d0)/a(k,k)
            call zscal(n-k,t,a(k+1,k),1)
c
c           row elimination with column indexing
c
            do 30 j = kp1, n
               t = a(l,j)
               if (l .eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
               call zaxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (cabs1(a(n,n)) .eq. 0.0d0) info = n
      return
      end
c
      subroutine zgesl(a,lda,n,ipvt,b,job)
      integer lda,n,ipvt(1),job
      complex*16 a(lda,1),b(1)
c
c     zgesl solves the complex*16 system
c     a * x = b  or  ctrans(a) * x = b
c     using the factors computed by zgeco or zgefa.
c
c     on entry
c
c        a       complex*16(lda, n)
c                the output from zgeco or zgefa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    integer(n)
c                the pivot vector from zgeco or zgefa.
c
c        b       complex*16(n)
c                the right hand side vector.
c
c        job     integer
c                = 0         to solve  a*x = b ,
c                = nonzero   to solve  ctrans(a)*x = b  where
c                            ctrans(a)  is the conjugate transpose.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains a
c        zero on the diagonal.  technically this indicates singularity
c        but it is often caused by improper arguments or improper
c        setting of lda .  it will not occur if the subroutines are
c        called correctly and if zgeco has set rcond .gt. 0.0
c        or zgefa has set info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call zgeco(a,lda,n,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           do 10 j = 1, p
c              call zgesl(a,lda,n,ipvt,c(1,j),0)
c        10 continue
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas zaxpy,zdotc
c     fortran dconjg
c
c     internal variables
c
      complex*16 zdotc,t
      integer k,kb,l,nm1
      double precision dreal,dimag
      complex*16 zdumr,zdumi
      dreal(zdumr) = zdumr
      dimag(zdumi) = (0.0d0,-1.0d0)*zdumi
c
      nm1 = n - 1
      if (job .ne. 0) go to 50
c
c        job = 0 , solve  a * x = b
c        first solve  l*y = b
c
         if (nm1 .lt. 1) go to 30
         do 20 k = 1, nm1
            l = ipvt(k)
            t = b(l)
            if (l .eq. k) go to 10
               b(l) = b(k)
               b(k) = t
   10       continue
            call zaxpy(n-k,t,a(k+1,k),1,b(k+1),1)
   20    continue
   30    continue
c
c        now solve  u*x = y
c
         do 40 kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/a(k,k)
            t = -b(k)
            call zaxpy(k-1,t,a(1,k),1,b(1),1)
   40    continue
      go to 100
   50 continue
c
c        job = nonzero, solve  ctrans(a) * x = b
c        first solve  ctrans(u)*y = b
c
         do 60 k = 1, n
            t = zdotc(k-1,a(1,k),1,b(1),1)
            b(k) = (b(k) - t)/dconjg(a(k,k))
   60    continue
c
c        now solve ctrans(l)*x = y
c
         if (nm1 .lt. 1) go to 90
         do 80 kb = 1, nm1
            k = n - kb
            b(k) = b(k) + zdotc(n-k,a(k+1,k),1,b(k+1),1)
            l = ipvt(k)
            if (l .eq. k) go to 70
               t = b(l)
               b(l) = b(k)
               b(k) = t
   70       continue
   80    continue
   90    continue
  100 continue
      return
      end
c
      subroutine cbabk2(nm,n,low,igh,scale,m,zr,zi)
c
      integer i,j,k,m,n,ii,nm,igh,low
      double precision scale(n),zr(nm,m),zi(nm,m)
      double precision s
c
c     this subroutine is a translation of the algol procedure
c     cbabk2, which is a complex version of balbak,
c     num. math. 13, 293-304(1969) by parlett and reinsch.
c     handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).
c
c     this subroutine forms the eigenvectors of a complex general
c     matrix by back transforming those of the corresponding
c     balanced matrix determined by  cbal.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        low and igh are integers determined by  cbal.
c
c        scale contains information determining the permutations
c          and scaling factors used by  cbal.
c
c        m is the number of eigenvectors to be back transformed.
c
c        zr and zi contain the real and imaginary parts,
c          respectively, of the eigenvectors to be
c          back transformed in their first m columns.
c
c     on output
c
c        zr and zi contain the real and imaginary parts,
c          respectively, of the transformed eigenvectors
c          in their first m columns.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      if (m .eq. 0) go to 200
      if (igh .eq. low) go to 120
c
      do 110 i = low, igh
         s = scale(i)
c     .......... left hand eigenvectors are back transformed
c                if the foregoing statement is replaced by
c                s=1.0d0/scale(i). ..........
         do 100 j = 1, m
            zr(i,j) = zr(i,j) * s
            zi(i,j) = zi(i,j) * s
  100    continue
c
  110 continue
c     .......... for i=low-1 step -1 until 1,
c                igh+1 step 1 until n do -- ..........
  120 do 140 ii = 1, n
         i = ii
         if (i .ge. low .and. i .le. igh) go to 140
         if (i .lt. low) i = low - ii
         k = scale(i)
         if (k .eq. i) go to 140
c
         do 130 j = 1, m
            s = zr(i,j)
            zr(i,j) = zr(k,j)
            zr(k,j) = s
            s = zi(i,j)
            zi(i,j) = zi(k,j)
            zi(k,j) = s
  130    continue
c
  140 continue
c
  200 return
      end
c
      subroutine cbal(nm,n,ar,ai,low,igh,scale)
c
      integer i,j,k,l,m,n,jj,nm,igh,low,iexc
      double precision ar(nm,n),ai(nm,n),scale(n)
      double precision c,f,g,r,s,b2,radix
      logical noconv
c
c     this subroutine is a translation of the algol procedure
c     cbalance, which is a complex version of balance,
c     num. math. 13, 293-304(1969) by parlett and reinsch.
c     handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).
c
c     this subroutine balances a complex matrix and isolates
c     eigenvalues whenever possible.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        ar and ai contain the real and imaginary parts,
c          respectively, of the complex matrix to be balanced.
c
c     on output
c
c        ar and ai contain the real and imaginary parts,
c          respectively, of the balanced matrix.
c
c        low and igh are two integers such that ar(i,j) and ai(i,j)
c          are equal to zero if
c           (1) i is greater than j and
c           (2) j=1,...,low-1 or i=igh+1,...,n.
c
c        scale contains information determining the
c           permutations and scaling factors used.
c
c     suppose that the principal submatrix in rows low through igh
c     has been balanced, that p(j) denotes the index interchanged
c     with j during the permutation step, and that the elements
c     of the diagonal matrix used are denoted by d(i,j).  then
c        scale(j) = p(j),    for j = 1,...,low-1
c                 = d(j,j)       j = low,...,igh
c                 = p(j)         j = igh+1,...,n.
c     the order in which the interchanges are made is n to igh+1,
c     then 1 to low-1.
c
c     note that 1 is returned for igh if igh is zero formally.
c
c     the algol procedure exc contained in cbalance appears in
c     cbal  in line.  (note that the algol roles of identifiers
c     k,l have been reversed.)
c
c     arithmetic is real throughout.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      radix = 16.0d0
c
      b2 = radix * radix
      k = 1
      l = n
      go to 100
c     .......... in-line procedure for row and
c                column exchange ..........
   20 scale(m) = j
      if (j .eq. m) go to 50
c
      do 30 i = 1, l
         f = ar(i,j)
         ar(i,j) = ar(i,m)
         ar(i,m) = f
         f = ai(i,j)
         ai(i,j) = ai(i,m)
         ai(i,m) = f
   30 continue
c
      do 40 i = k, n
         f = ar(j,i)
         ar(j,i) = ar(m,i)
         ar(m,i) = f
         f = ai(j,i)
         ai(j,i) = ai(m,i)
         ai(m,i) = f
   40 continue
c
   50 go to (80,130), iexc
c     .......... search for rows isolating an eigenvalue
c                and push them down ..........
   80 if (l .eq. 1) go to 280
      l = l - 1
c     .......... for j=l step -1 until 1 do -- ..........
  100 do 120 jj = 1, l
         j = l + 1 - jj
c
         do 110 i = 1, l
            if (i .eq. j) go to 110
            if (ar(j,i) .ne. 0.0d0 .or. ai(j,i) .ne. 0.0d0) go to 120
  110    continue
c
         m = l
         iexc = 1
         go to 20
  120 continue
c
      go to 140
c     .......... search for columns isolating an eigenvalue
c                and push them left ..........
  130 k = k + 1
c
  140 do 170 j = k, l
c
         do 150 i = k, l
            if (i .eq. j) go to 150
            if (ar(i,j) .ne. 0.0d0 .or. ai(i,j) .ne. 0.0d0) go to 170
  150    continue
c
         m = k
         iexc = 2
         go to 20
  170 continue
c     .......... now balance the submatrix in rows k to l ..........
      do 180 i = k, l
  180 scale(i) = 1.0d0
c     .......... iterative loop for norm reduction ..........
  190 noconv = .false.
c
      do 270 i = k, l
         c = 0.0d0
         r = 0.0d0
c
         do 200 j = k, l
            if (j .eq. i) go to 200
            c = c + dabs(ar(j,i)) + dabs(ai(j,i))
            r = r + dabs(ar(i,j)) + dabs(ai(i,j))
  200    continue
c     .......... guard against zero c or r due to underflow ..........
         if (c .eq. 0.0d0 .or. r .eq. 0.0d0) go to 270
         g = r / radix
         f = 1.0d0
         s = c + r
  210    if (c .ge. g) go to 220
         f = f * radix
         c = c * b2
         go to 210
  220    g = r * radix
  230    if (c .lt. g) go to 240
         f = f / radix
         c = c / b2
         go to 230
c     .......... now balance ..........
  240    if ((c + r) / f .ge. 0.95d0 * s) go to 270
         g = 1.0d0 / f
         scale(i) = scale(i) * f
         noconv = .true.
c
         do 250 j = k, n
            ar(i,j) = ar(i,j) * g
            ai(i,j) = ai(i,j) * g
  250    continue
c
         do 260 j = 1, l
            ar(j,i) = ar(j,i) * f
            ai(j,i) = ai(j,i) * f
  260    continue
c
  270 continue
c
      if (noconv) go to 190
c
  280 low = k
      igh = l
      return
      end
c
      subroutine cdiv(ar,ai,br,bi,cr,ci)
      double precision ar,ai,br,bi,cr,ci
c
c     complex division, (cr,ci) = (ar,ai)/(br,bi)
c
      double precision s,ars,ais,brs,bis
      s = dabs(br) + dabs(bi)
      ars = ar/s
      ais = ai/s
      brs = br/s
      bis = bi/s
      s = brs**2 + bis**2
      cr = (ars*brs + ais*bis)/s
      ci = (ais*brs - ars*bis)/s
      return
      end
      subroutine cg(nm,n,ar,ai,wr,wi,matz,zr,zi,fv1,fv2,fv3,ierr)
c
      integer n,nm,is1,is2,ierr,matz
      double precision ar(nm,n),ai(nm,n),wr(n),wi(n),zr(nm,n),zi(nm,n),
     x       fv1(n),fv2(n),fv3(n)
c
c     this subroutine calls the recommended sequence of
c     subroutines from the eigensystem subroutine package (eispack)
c     to find the eigenvalues and eigenvectors (if desired)
c     of a complex general matrix.
c
c     on input
c
c        nm  must be set to the row dimension of the two-dimensional
c        array parameters as declared in the calling program
c        dimension statement.
c
c        n  is the order of the matrix  a=(ar,ai).
c
c        ar  and  ai  contain the real and imaginary parts,
c        respectively, of the complex general matrix.
c
c        matz  is an integer variable set equal to zero if
c        only eigenvalues are desired.  otherwise it is set to
c        any non-zero integer for both eigenvalues and eigenvectors.
c
c     on output
c
c        wr  and  wi  contain the real and imaginary parts,
c        respectively, of the eigenvalues.
c
c        zr  and  zi  contain the real and imaginary parts,
c        respectively, of the eigenvectors if matz is not zero.
c
c        ierr  is an integer output variable set equal to an error
c           completion code described in the documentation for comqr
c           and comqr2.  the normal completion code is zero.
c
c        fv1, fv2, and  fv3  are temporary storage arrays.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      if (n .le. nm) go to 10
      ierr = 10 * n
      go to 50
c
   10 call  cbal(nm,n,ar,ai,is1,is2,fv1)
      call  corth(nm,n,is1,is2,ar,ai,fv2,fv3)
      if (matz .ne. 0) go to 20
c     .......... find eigenvalues only ..........
      call  comqr(nm,n,is1,is2,ar,ai,wr,wi,ierr)
      go to 50
c     .......... find both eigenvalues and eigenvectors ..........
   20 call  comqr2(nm,n,is1,is2,fv2,fv3,ar,ai,wr,wi,zr,zi,ierr)
      if (ierr .ne. 0) go to 50
      call  cbabk2(nm,n,is1,is2,fv1,n,zr,zi)
   50 return
      end
      subroutine comqr(nm,n,low,igh,hr,hi,wr,wi,ierr)
c
      integer i,j,l,n,en,ll,nm,igh,itn,its,low,lp1,enm1,ierr
      double precision hr(nm,n),hi(nm,n),wr(n),wi(n)
      double precision si,sr,ti,tr,xi,xr,yi,yr,zzi,zzr,norm,tst1,tst2,
     x       pythag
c
c     this subroutine is a translation of a unitary analogue of the
c     algol procedure  comlr, num. math. 12, 369-376(1968) by martin
c     and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 396-403(1971).
c     the unitary analogue substitutes the qr algorithm of francis
c     (comp. jour. 4, 332-345(1962)) for the lr algorithm.
c
c     this subroutine finds the eigenvalues of a complex
c     upper hessenberg matrix by the qr method.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        low and igh are integers determined by the balancing
c          subroutine  cbal.  if  cbal  has not been used,
c          set low=1, igh=n.
c
c        hr and hi contain the real and imaginary parts,
c          respectively, of the complex upper hessenberg matrix.
c          their lower triangles below the subdiagonal contain
c          information about the unitary transformations used in
c          the reduction by  corth, if performed.
c
c     on output
c
c        the upper hessenberg portions of hr and hi have been
c          destroyed.  therefore, they must be saved before
c          calling  comqr  if subsequent calculation of
c          eigenvectors is to be performed.
c
c        wr and wi contain the real and imaginary parts,
c          respectively, of the eigenvalues.  if an error
c          exit is made, the eigenvalues should be correct
c          for indices ierr+1,...,n.
c
c        ierr is set to
c          zero       for normal return,
c          j          if the limit of 30*n iterations is exhausted
c                     while the j-th eigenvalue is being sought.
c
c     calls cdiv for complex division.
c     calls csroot for complex square root.
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      ierr = 0
      if (low .eq. igh) go to 180
c     .......... create real subdiagonal elements ..........
      l = low + 1
c
      do 170 i = l, igh
         ll = min0(i+1,igh)
         if (hi(i,i-1) .eq. 0.0d0) go to 170
         norm = pythag(hr(i,i-1),hi(i,i-1))
         yr = hr(i,i-1) / norm
         yi = hi(i,i-1) / norm
         hr(i,i-1) = norm
         hi(i,i-1) = 0.0d0
c
         do 155 j = i, igh
            si = yr * hi(i,j) - yi * hr(i,j)
            hr(i,j) = yr * hr(i,j) + yi * hi(i,j)
            hi(i,j) = si
  155    continue
c
         do 160 j = low, ll
            si = yr * hi(j,i) + yi * hr(j,i)
            hr(j,i) = yr * hr(j,i) - yi * hi(j,i)
            hi(j,i) = si
  160    continue
c
  170 continue
c     .......... store roots isolated by cbal ..........
  180 do 200 i = 1, n
         if (i .ge. low .and. i .le. igh) go to 200
         wr(i) = hr(i,i)
         wi(i) = hi(i,i)
  200 continue
c
      en = igh
      tr = 0.0d0
      ti = 0.0d0
      itn = 30*n
c     .......... search for next eigenvalue ..........
  220 if (en .lt. low) go to 1001
      its = 0
      enm1 = en - 1
c     .......... look for single small sub-diagonal element
c                for l=en step -1 until low d0 -- ..........
  240 do 260 ll = low, en
         l = en + low - ll
         if (l .eq. low) go to 300
         tst1 = dabs(hr(l-1,l-1)) + dabs(hi(l-1,l-1))
     x            + dabs(hr(l,l)) + dabs(hi(l,l))
         tst2 = tst1 + dabs(hr(l,l-1))
         if (tst2 .eq. tst1) go to 300
  260 continue
c     .......... form shift ..........
  300 if (l .eq. en) go to 660
      if (itn .eq. 0) go to 1000
      if (its .eq. 10 .or. its .eq. 20) go to 320
      sr = hr(en,en)
      si = hi(en,en)
      xr = hr(enm1,en) * hr(en,enm1)
      xi = hi(enm1,en) * hr(en,enm1)
      if (xr .eq. 0.0d0 .and. xi .eq. 0.0d0) go to 340
      yr = (hr(enm1,enm1) - sr) / 2.0d0
      yi = (hi(enm1,enm1) - si) / 2.0d0
      call csroot(yr**2-yi**2+xr,2.0d0*yr*yi+xi,zzr,zzi)
      if (yr * zzr + yi * zzi .ge. 0.0d0) go to 310
      zzr = -zzr
      zzi = -zzi
  310 call cdiv(xr,xi,yr+zzr,yi+zzi,xr,xi)
      sr = sr - xr
      si = si - xi
      go to 340
c     .......... form exceptional shift ..........
  320 sr = dabs(hr(en,enm1)) + dabs(hr(enm1,en-2))
      si = 0.0d0
c
  340 do 360 i = low, en
         hr(i,i) = hr(i,i) - sr
         hi(i,i) = hi(i,i) - si
  360 continue
c
      tr = tr + sr
      ti = ti + si
      its = its + 1
      itn = itn - 1
c     .......... reduce to triangle (rows) ..........
      lp1 = l + 1
c
      do 500 i = lp1, en
         sr = hr(i,i-1)
         hr(i,i-1) = 0.0d0
         norm = pythag(pythag(hr(i-1,i-1),hi(i-1,i-1)),sr)
         xr = hr(i-1,i-1) / norm
         wr(i-1) = xr
         xi = hi(i-1,i-1) / norm
         wi(i-1) = xi
         hr(i-1,i-1) = norm
         hi(i-1,i-1) = 0.0d0
         hi(i,i-1) = sr / norm
c
         do 490 j = i, en
            yr = hr(i-1,j)
            yi = hi(i-1,j)
            zzr = hr(i,j)
            zzi = hi(i,j)
            hr(i-1,j) = xr * yr + xi * yi + hi(i,i-1) * zzr
            hi(i-1,j) = xr * yi - xi * yr + hi(i,i-1) * zzi
            hr(i,j) = xr * zzr - xi * zzi - hi(i,i-1) * yr
            hi(i,j) = xr * zzi + xi * zzr - hi(i,i-1) * yi
  490    continue
c
  500 continue
c
      si = hi(en,en)
      if (si .eq. 0.0d0) go to 540
      norm = pythag(hr(en,en),si)
      sr = hr(en,en) / norm
      si = si / norm
      hr(en,en) = norm
      hi(en,en) = 0.0d0
c     .......... inverse operation (columns) ..........
  540 do 600 j = lp1, en
         xr = wr(j-1)
         xi = wi(j-1)
c
         do 580 i = l, j
            yr = hr(i,j-1)
            yi = 0.0d0
            zzr = hr(i,j)
            zzi = hi(i,j)
            if (i .eq. j) go to 560
            yi = hi(i,j-1)
            hi(i,j-1) = xr * yi + xi * yr + hi(j,j-1) * zzi
  560       hr(i,j-1) = xr * yr - xi * yi + hi(j,j-1) * zzr
            hr(i,j) = xr * zzr + xi * zzi - hi(j,j-1) * yr
            hi(i,j) = xr * zzi - xi * zzr - hi(j,j-1) * yi
  580    continue
c
  600 continue
c
      if (si .eq. 0.0d0) go to 240
c
      do 630 i = l, en
         yr = hr(i,en)
         yi = hi(i,en)
         hr(i,en) = sr * yr - si * yi
         hi(i,en) = sr * yi + si * yr
  630 continue
c
      go to 240
c     .......... a root found ..........
  660 wr(en) = hr(en,en) + tr
      wi(en) = hi(en,en) + ti
      en = enm1
      go to 220
c     .......... set error -- all eigenvalues have not
c                converged after 30*n iterations ..........
 1000 ierr = en
 1001 return
      end
c
      subroutine comqr2(nm,n,low,igh,ortr,orti,hr,hi,wr,wi,zr,zi,ierr)
C  MESHED overflow control WITH vectors of isolated roots (10/19/89 BSG)
C  MESHED overflow control WITH triangular multiply (10/30/89 BSG)
c
      integer i,j,k,l,m,n,en,ii,jj,ll,nm,nn,igh,ip1,
     x        itn,its,low,lp1,enm1,iend,ierr
      double precision hr(nm,n),hi(nm,n),wr(n),wi(n),zr(nm,n),zi(nm,n),
     x       ortr(igh),orti(igh)
      double precision si,sr,ti,tr,xi,xr,yi,yr,zzi,zzr,norm,tst1,tst2,
     x       pythag
c
c     this subroutine is a translation of a unitary analogue of the
c     algol procedure  comlr2, num. math. 16, 181-204(1970) by peters
c     and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).
c     the unitary analogue substitutes the qr algorithm of francis
c     (comp. jour. 4, 332-345(1962)) for the lr algorithm.
c
c     this subroutine finds the eigenvalues and eigenvectors
c     of a complex upper hessenberg matrix by the qr
c     method.  the eigenvectors of a complex general matrix
c     can also be found if  corth  has been used to reduce
c     this general matrix to hessenberg form.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        low and igh are integers determined by the balancing
c          subroutine  cbal.  if  cbal  has not been used,
c          set low=1, igh=n.
c
c        ortr and orti contain information about the unitary trans-
c          formations used in the reduction by  corth, if performed.
c          only elements low through igh are used.  if the eigenvectors
c          of the hessenberg matrix are desired, set ortr(j) and
c          orti(j) to 0.0d0 for these elements.
c
c        hr and hi contain the real and imaginary parts,
c          respectively, of the complex upper hessenberg matrix.
c          their lower triangles below the subdiagonal contain further
c          information about the transformations which were used in the
c          reduction by  corth, if performed.  if the eigenvectors of
c          the hessenberg matrix are desired, these elements may be
c          arbitrary.
c
c     on output
c
c        ortr, orti, and the upper hessenberg portions of hr and hi
c          have been destroyed.
c
c        wr and wi contain the real and imaginary parts,
c          respectively, of the eigenvalues.  if an error
c          exit is made, the eigenvalues should be correct
c          for indices ierr+1,...,n.
c
c        zr and zi contain the real and imaginary parts,
c          respectively, of the eigenvectors.  the eigenvectors
c          are unnormalized.  if an error exit is made, none of
c          the eigenvectors has been found.
c
c        ierr is set to
c          zero       for normal return,
c          j          if the limit of 30*n iterations is exhausted
c                     while the j-th eigenvalue is being sought.
c
c     calls cdiv for complex division.
c     calls csroot for complex square root.
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated october 1989.
c
c     ------------------------------------------------------------------
c
      ierr = 0
c     .......... initialize eigenvector matrix ..........
      do 101 j = 1, n
c
         do 100 i = 1, n
            zr(i,j) = 0.0d0
            zi(i,j) = 0.0d0
  100    continue
         zr(j,j) = 1.0d0
  101 continue
c     .......... form the matrix of accumulated transformations
c                from the information left by corth ..........
      iend = igh - low - 1
      if (iend) 180, 150, 105
c     .......... for i=igh-1 step -1 until low+1 do -- ..........
  105 do 140 ii = 1, iend
         i = igh - ii
         if (ortr(i) .eq. 0.0d0 .and. orti(i) .eq. 0.0d0) go to 140
         if (hr(i,i-1) .eq. 0.0d0 .and. hi(i,i-1) .eq. 0.0d0) go to 140
c     .......... norm below is negative of h formed in corth ..........
         norm = hr(i,i-1) * ortr(i) + hi(i,i-1) * orti(i)
         ip1 = i + 1
c
         do 110 k = ip1, igh
            ortr(k) = hr(k,i-1)
            orti(k) = hi(k,i-1)
  110    continue
c
         do 130 j = i, igh
            sr = 0.0d0
            si = 0.0d0
c
            do 115 k = i, igh
               sr = sr + ortr(k) * zr(k,j) + orti(k) * zi(k,j)
               si = si + ortr(k) * zi(k,j) - orti(k) * zr(k,j)
  115       continue
c
            sr = sr / norm
            si = si / norm
c
            do 120 k = i, igh
               zr(k,j) = zr(k,j) + sr * ortr(k) - si * orti(k)
               zi(k,j) = zi(k,j) + sr * orti(k) + si * ortr(k)
  120       continue
c
  130    continue
c
  140 continue
c     .......... create real subdiagonal elements ..........
  150 l = low + 1
c
      do 170 i = l, igh
         ll = min0(i+1,igh)
         if (hi(i,i-1) .eq. 0.0d0) go to 170
         norm = pythag(hr(i,i-1),hi(i,i-1))
         yr = hr(i,i-1) / norm
         yi = hi(i,i-1) / norm
         hr(i,i-1) = norm
         hi(i,i-1) = 0.0d0
c
         do 155 j = i, n
            si = yr * hi(i,j) - yi * hr(i,j)
            hr(i,j) = yr * hr(i,j) + yi * hi(i,j)
            hi(i,j) = si
  155    continue
c
         do 160 j = 1, ll
            si = yr * hi(j,i) + yi * hr(j,i)
            hr(j,i) = yr * hr(j,i) - yi * hi(j,i)
            hi(j,i) = si
  160    continue
c
         do 165 j = low, igh
            si = yr * zi(j,i) + yi * zr(j,i)
            zr(j,i) = yr * zr(j,i) - yi * zi(j,i)
            zi(j,i) = si
  165    continue
c
  170 continue
c     .......... store roots isolated by cbal ..........
  180 do 200 i = 1, n
         if (i .ge. low .and. i .le. igh) go to 200
         wr(i) = hr(i,i)
         wi(i) = hi(i,i)
  200 continue
c
      en = igh
      tr = 0.0d0
      ti = 0.0d0
      itn = 30*n
c     .......... search for next eigenvalue ..........
  220 if (en .lt. low) go to 680
      its = 0
      enm1 = en - 1
c     .......... look for single small sub-diagonal element
c                for l=en step -1 until low do -- ..........
  240 do 260 ll = low, en
         l = en + low - ll
         if (l .eq. low) go to 300
         tst1 = dabs(hr(l-1,l-1)) + dabs(hi(l-1,l-1))
     x            + dabs(hr(l,l)) + dabs(hi(l,l))
         tst2 = tst1 + dabs(hr(l,l-1))
         if (tst2 .eq. tst1) go to 300
  260 continue
c     .......... form shift ..........
  300 if (l .eq. en) go to 660
      if (itn .eq. 0) go to 1000
      if (its .eq. 10 .or. its .eq. 20) go to 320
      sr = hr(en,en)
      si = hi(en,en)
      xr = hr(enm1,en) * hr(en,enm1)
      xi = hi(enm1,en) * hr(en,enm1)
      if (xr .eq. 0.0d0 .and. xi .eq. 0.0d0) go to 340
      yr = (hr(enm1,enm1) - sr) / 2.0d0
      yi = (hi(enm1,enm1) - si) / 2.0d0
      call csroot(yr**2-yi**2+xr,2.0d0*yr*yi+xi,zzr,zzi)
      if (yr * zzr + yi * zzi .ge. 0.0d0) go to 310
      zzr = -zzr
      zzi = -zzi
  310 call cdiv(xr,xi,yr+zzr,yi+zzi,xr,xi)
      sr = sr - xr
      si = si - xi
      go to 340
c     .......... form exceptional shift ..........
  320 sr = dabs(hr(en,enm1)) + dabs(hr(enm1,en-2))
      si = 0.0d0
c
  340 do 360 i = low, en
         hr(i,i) = hr(i,i) - sr
         hi(i,i) = hi(i,i) - si
  360 continue
c
      tr = tr + sr
      ti = ti + si
      its = its + 1
      itn = itn - 1
c     .......... reduce to triangle (rows) ..........
      lp1 = l + 1
c
      do 500 i = lp1, en
         sr = hr(i,i-1)
         hr(i,i-1) = 0.0d0
         norm = pythag(pythag(hr(i-1,i-1),hi(i-1,i-1)),sr)
         xr = hr(i-1,i-1) / norm
         wr(i-1) = xr
         xi = hi(i-1,i-1) / norm
         wi(i-1) = xi
         hr(i-1,i-1) = norm
         hi(i-1,i-1) = 0.0d0
         hi(i,i-1) = sr / norm
c
         do 490 j = i, n
            yr = hr(i-1,j)
            yi = hi(i-1,j)
            zzr = hr(i,j)
            zzi = hi(i,j)
            hr(i-1,j) = xr * yr + xi * yi + hi(i,i-1) * zzr
            hi(i-1,j) = xr * yi - xi * yr + hi(i,i-1) * zzi
            hr(i,j) = xr * zzr - xi * zzi - hi(i,i-1) * yr
            hi(i,j) = xr * zzi + xi * zzr - hi(i,i-1) * yi
  490    continue
c
  500 continue
c
      si = hi(en,en)
      if (si .eq. 0.0d0) go to 540
      norm = pythag(hr(en,en),si)
      sr = hr(en,en) / norm
      si = si / norm
      hr(en,en) = norm
      hi(en,en) = 0.0d0
      if (en .eq. n) go to 540
      ip1 = en + 1
c
      do 520 j = ip1, n
         yr = hr(en,j)
         yi = hi(en,j)
         hr(en,j) = sr * yr + si * yi
         hi(en,j) = sr * yi - si * yr
  520 continue
c     .......... inverse operation (columns) ..........
  540 do 600 j = lp1, en
         xr = wr(j-1)
         xi = wi(j-1)
c
         do 580 i = 1, j
            yr = hr(i,j-1)
            yi = 0.0d0
            zzr = hr(i,j)
            zzi = hi(i,j)
            if (i .eq. j) go to 560
            yi = hi(i,j-1)
            hi(i,j-1) = xr * yi + xi * yr + hi(j,j-1) * zzi
  560       hr(i,j-1) = xr * yr - xi * yi + hi(j,j-1) * zzr
            hr(i,j) = xr * zzr + xi * zzi - hi(j,j-1) * yr
            hi(i,j) = xr * zzi - xi * zzr - hi(j,j-1) * yi
  580    continue
c
         do 590 i = low, igh
            yr = zr(i,j-1)
            yi = zi(i,j-1)
            zzr = zr(i,j)
            zzi = zi(i,j)
            zr(i,j-1) = xr * yr - xi * yi + hi(j,j-1) * zzr
            zi(i,j-1) = xr * yi + xi * yr + hi(j,j-1) * zzi
            zr(i,j) = xr * zzr + xi * zzi - hi(j,j-1) * yr
            zi(i,j) = xr * zzi - xi * zzr - hi(j,j-1) * yi
  590    continue
c
  600 continue
c
      if (si .eq. 0.0d0) go to 240
c
      do 630 i = 1, en
         yr = hr(i,en)
         yi = hi(i,en)
         hr(i,en) = sr * yr - si * yi
         hi(i,en) = sr * yi + si * yr
  630 continue
c
      do 640 i = low, igh
         yr = zr(i,en)
         yi = zi(i,en)
         zr(i,en) = sr * yr - si * yi
         zi(i,en) = sr * yi + si * yr
  640 continue
c
      go to 240
c     .......... a root found ..........
  660 hr(en,en) = hr(en,en) + tr
      wr(en) = hr(en,en)
      hi(en,en) = hi(en,en) + ti
      wi(en) = hi(en,en)
      en = enm1
      go to 220
c     .......... all roots found.  backsubstitute to find
c                vectors of upper triangular form ..........
  680 norm = 0.0d0
c
      do 720 i = 1, n
c
         do 720 j = i, n
            tr = dabs(hr(i,j)) + dabs(hi(i,j))
            if (tr .gt. norm) norm = tr
  720 continue
c
      if (n .eq. 1 .or. norm .eq. 0.0d0) go to 1001
c     .......... for en=n step -1 until 2 do -- ..........
      do 800 nn = 2, n
         en = n + 2 - nn
         xr = wr(en)
         xi = wi(en)
         hr(en,en) = 1.0d0
         hi(en,en) = 0.0d0
         enm1 = en - 1
c     .......... for i=en-1 step -1 until 1 do -- ..........
         do 780 ii = 1, enm1
            i = en - ii
            zzr = 0.0d0
            zzi = 0.0d0
            ip1 = i + 1
c
            do 740 j = ip1, en
               zzr = zzr + hr(i,j) * hr(j,en) - hi(i,j) * hi(j,en)
               zzi = zzi + hr(i,j) * hi(j,en) + hi(i,j) * hr(j,en)
  740       continue
c
            yr = xr - wr(i)
            yi = xi - wi(i)
            if (yr .ne. 0.0d0 .or. yi .ne. 0.0d0) go to 765
               tst1 = norm
               yr = tst1
  760          yr = 0.01d0 * yr
               tst2 = norm + yr
               if (tst2 .gt. tst1) go to 760
  765       continue
            call cdiv(zzr,zzi,yr,yi,hr(i,en),hi(i,en))
c     .......... overflow control ..........
            tr = dabs(hr(i,en)) + dabs(hi(i,en))
            if (tr .eq. 0.0d0) go to 780
            tst1 = tr
            tst2 = tst1 + 1.0d0/tst1
            if (tst2 .gt. tst1) go to 780
            do 770 j = i, en
               hr(j,en) = hr(j,en)/tr
               hi(j,en) = hi(j,en)/tr
  770       continue
c
  780    continue
c
  800 continue
c     .......... end backsubstitution ..........
c     .......... vectors of isolated roots ..........
      do  840 i = 1, N
         if (i .ge. low .and. i .le. igh) go to 840
c
         do 820 j = I, n
            zr(i,j) = hr(i,j)
            zi(i,j) = hi(i,j)
  820    continue
c
  840 continue
c     .......... multiply by transformation matrix to give
c                vectors of original full matrix.
c                for j=n step -1 until low do -- ..........
      do 880 jj = low, N
         j = n + low - jj
         m = min0(j,igh)
c
         do 880 i = low, igh
            zzr = 0.0d0
            zzi = 0.0d0
c
            do 860 k = low, m
               zzr = zzr + zr(i,k) * hr(k,j) - zi(i,k) * hi(k,j)
               zzi = zzi + zr(i,k) * hi(k,j) + zi(i,k) * hr(k,j)
  860       continue
c
            zr(i,j) = zzr
            zi(i,j) = zzi
  880 continue
c
      go to 1001
c     .......... set error -- all eigenvalues have not
c                converged after 30*n iterations ..........
 1000 ierr = en
 1001 return
      end
c
      subroutine corth(nm,n,low,igh,ar,ai,ortr,orti)
c
      integer i,j,m,n,ii,jj,la,mp,nm,igh,kp1,low
      double precision ar(nm,n),ai(nm,n),ortr(igh),orti(igh)
      double precision f,g,h,fi,fr,scale,pythag
c
c     this subroutine is a translation of a complex analogue of
c     the algol procedure orthes, num. math. 12, 349-368(1968)
c     by martin and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).
c
c     given a complex general matrix, this subroutine
c     reduces a submatrix situated in rows and columns
c     low through igh to upper hessenberg form by
c     unitary similarity transformations.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        low and igh are integers determined by the balancing
c          subroutine  cbal.  if  cbal  has not been used,
c          set low=1, igh=n.
c
c        ar and ai contain the real and imaginary parts,
c          respectively, of the complex input matrix.
c
c     on output
c
c        ar and ai contain the real and imaginary parts,
c          respectively, of the hessenberg matrix.  information
c          about the unitary transformations used in the reduction
c          is stored in the remaining triangles under the
c          hessenberg matrix.
c
c        ortr and orti contain further information about the
c          transformations.  only elements low through igh are used.
c
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      la = igh - 1
      kp1 = low + 1
      if (la .lt. kp1) go to 200
c
      do 180 m = kp1, la
         h = 0.0d0
         ortr(m) = 0.0d0
         orti(m) = 0.0d0
         scale = 0.0d0
c     .......... scale column (algol tol then not needed) ..........
         do 90 i = m, igh
   90    scale = scale + dabs(ar(i,m-1)) + dabs(ai(i,m-1))
c
         if (scale .eq. 0.0d0) go to 180
         mp = m + igh
c     .......... for i=igh step -1 until m do -- ..........
         do 100 ii = m, igh
            i = mp - ii
            ortr(i) = ar(i,m-1) / scale
            orti(i) = ai(i,m-1) / scale
            h = h + ortr(i) * ortr(i) + orti(i) * orti(i)
  100    continue
c
         g = dsqrt(h)
         f = pythag(ortr(m),orti(m))
         if (f .eq. 0.0d0) go to 103
         h = h + f * g
         g = g / f
         ortr(m) = (1.0d0 + g) * ortr(m)
         orti(m) = (1.0d0 + g) * orti(m)
         go to 105
c
  103    ortr(m) = g
         ar(m,m-1) = scale
c     .......... form (i-(u*ut)/h) * a ..........
  105    do 130 j = m, n
            fr = 0.0d0
            fi = 0.0d0
c     .......... for i=igh step -1 until m do -- ..........
            do 110 ii = m, igh
               i = mp - ii
               fr = fr + ortr(i) * ar(i,j) + orti(i) * ai(i,j)
               fi = fi + ortr(i) * ai(i,j) - orti(i) * ar(i,j)
  110       continue
c
            fr = fr / h
            fi = fi / h
c
            do 120 i = m, igh
               ar(i,j) = ar(i,j) - fr * ortr(i) + fi * orti(i)
               ai(i,j) = ai(i,j) - fr * orti(i) - fi * ortr(i)
  120       continue
c
  130    continue
c     .......... form (i-(u*ut)/h)*a*(i-(u*ut)/h) ..........
         do 160 i = 1, igh
            fr = 0.0d0
            fi = 0.0d0
c     .......... for j=igh step -1 until m do -- ..........
            do 140 jj = m, igh
               j = mp - jj
               fr = fr + ortr(j) * ar(i,j) - orti(j) * ai(i,j)
               fi = fi + ortr(j) * ai(i,j) + orti(j) * ar(i,j)
  140       continue
c
            fr = fr / h
            fi = fi / h
c
            do 150 j = m, igh
               ar(i,j) = ar(i,j) - fr * ortr(j) - fi * orti(j)
               ai(i,j) = ai(i,j) + fr * orti(j) - fi * ortr(j)
  150       continue
c
  160    continue
c
         ortr(m) = scale * ortr(m)
         orti(m) = scale * orti(m)
         ar(m,m-1) = -g * ar(m,m-1)
         ai(m,m-1) = -g * ai(m,m-1)
  180 continue
c
  200 return
      end
c
      subroutine csroot(xr,xi,yr,yi)
      double precision xr,xi,yr,yi
c
c     (yr,yi) = complex dsqrt(xr,xi) 
c     branch chosen so that yr .ge. 0.0 and sign(yi) .eq. sign(xi)
c
      double precision s,tr,ti,pythag
      tr = xr
      ti = xi
      s = dsqrt(0.5d0*(pythag(tr,ti) + dabs(tr)))
      if (tr .ge. 0.0d0) yr = s
      if (ti .lt. 0.0d0) s = -s
      if (tr .le. 0.0d0) yi = s
      if (tr .lt. 0.0d0) yr = 0.5d0*(ti/yi)
      if (tr .gt. 0.0d0) yi = 0.5d0*(ti/yr)
      return
      end
c
      double precision function pythag(a,b)
      double precision a,b
c
c     finds dsqrt(a**2+b**2) without overflow or destructive underflow
c
      double precision p,r,s,t,u
      p = dmax1(dabs(a),dabs(b))
      if (p .eq. 0.0d0) go to 20
      r = (dmin1(dabs(a),dabs(b))/p)**2
   10 continue
         t = 4.0d0 + r
         if (t .eq. 4.0d0) go to 20
         s = r/t
         u = 1.0d0 + 2.0d0*s
         p = u*p
         r = (s/u)**2 * r
      go to 10
   20 pythag = p
      return
      end
c
C
C
      subroutine dgesl(a,lda,n,ipvt,b,job)
      integer lda,n,ipvt(1),job
      double precision a(lda,1),b(1)
c
c     dgesl solves the double precision system
c     a * x = b  or  trans(a) * x = b
c     using the factors computed by dgeco or dgefa.
c
c     on entry
c
c        a       double precision(lda, n)
c                the output from dgeco or dgefa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    integer(n)
c                the pivot vector from dgeco or dgefa.
c
c        b       double precision(n)
c                the right hand side vector.
c
c        job     integer
c                = 0         to solve  a*x = b ,
c                = nonzero   to solve  trans(a)*x = b  where
c                            trans(a)  is the transpose.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains a
c        zero on the diagonal.  technically this indicates singularity
c        but it is often caused by improper arguments or improper
c        setting of lda .  it will not occur if the subroutines are
c        called correctly and if dgeco has set rcond .gt. 0.0
c        or dgefa has set info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call dgeco(a,lda,n,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           do 10 j = 1, p
c              call dgesl(a,lda,n,ipvt,c(1,j),0)
c        10 continue
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,ddot
c
c     internal variables
c
      double precision ddot,t
      integer k,kb,l,nm1
c
      nm1 = n - 1
      if (job .ne. 0) go to 50
c
c        job = 0 , solve  a * x = b
c        first solve  l*y = b
c
         if (nm1 .lt. 1) go to 30
         do 20 k = 1, nm1
            l = ipvt(k)
            t = b(l)
            if (l .eq. k) go to 10
               b(l) = b(k)
               b(k) = t
   10       continue
            call daxpy(n-k,t,a(k+1,k),1,b(k+1),1)
   20    continue
   30    continue
c
c        now solve  u*x = y
c
         do 40 kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/a(k,k)
            t = -b(k)
            call daxpy(k-1,t,a(1,k),1,b(1),1)
   40    continue
      go to 100
   50 continue
c
c        job = nonzero, solve  trans(a) * x = b
c        first solve  trans(u)*y = b
c
         do 60 k = 1, n
            t = ddot(k-1,a(1,k),1,b(1),1)
            b(k) = (b(k) - t)/a(k,k)
   60    continue
c
c        now solve trans(l)*x = y
c
         if (nm1 .lt. 1) go to 90
         do 80 kb = 1, nm1
            k = n - kb
            b(k) = b(k) + ddot(n-k,a(k+1,k),1,b(k+1),1)
            l = ipvt(k)
            if (l .eq. k) go to 70
               t = b(l)
               b(l) = b(k)
               b(k) = t
   70       continue
   80    continue
   90    continue
  100 continue
      return
      end
      double precision function ddot(n,dx,incx,dy,incy)
c
c     forms the dot product of two vectors.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(1),dy(1),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      ddot = 0.0d0
      dtemp = 0.0d0
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dtemp + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      ddot = dtemp
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dtemp + dx(i)*dy(i)
   30 continue
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +
     *   dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
   50 continue
   60 ddot = dtemp
      return
      end
      subroutine daxpy(n,da,dx,incx,dy,incy)
c
c     constant times a vector plus a vector.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(1),dy(1),da
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if (da .eq. 0.0d0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dy(i) + da*dx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        dy(i) = dy(i) + da*dx(i)
        dy(i + 1) = dy(i + 1) + da*dx(i + 1)
        dy(i + 2) = dy(i + 2) + da*dx(i + 2)
        dy(i + 3) = dy(i + 3) + da*dx(i + 3)
   50 continue
      return
      end
      subroutine dgeco(a,lda,n,ipvt,rcond,z)
      integer lda,n,ipvt(1)
      double precision a(lda,1),z(1)
      double precision rcond
c
c     dgeco factors a double precision matrix by gaussian elimination
c     and estimates the condition of the matrix.
c
c     if  rcond  is not needed, dgefa is slightly faster.
c     to solve  a*x = b , follow dgeco by dgesl.
c     to compute  inverse(a)*c , follow dgeco by dgesl.
c     to compute  determinant(a) , follow dgeco by dgedi.
c     to compute  inverse(a) , follow dgeco by dgedi.
c
c     on entry
c
c        a       double precision(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        rcond   double precision
c                an estimate of the reciprocal condition of  a .
c                for the system  a*x = b , relative perturbations
c                in  a  and  b  of size  epsilon  may cause
c                relative perturbations in  x  of size  epsilon/rcond .
c                if  rcond  is so small that the logical expression
c                           1.0 + rcond .eq. 1.0
c                is true, then  a  may be singular to working
c                precision.  in particular,  rcond  is zero  if
c                exact singularity is detected or the estimate
c                underflows.
c
c        z       double precision(n)
c                a work vector whose contents are usually unimportant.
c                if  a  is close to a singular matrix, then  z  is
c                an approximate null vector in the sense that
c                norm(a*z) = rcond*norm(a)*norm(z) .
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     linpack dgefa
c     blas daxpy,ddot,dscal,dasum
c     fortran dabs,dmax1,dsign
c
c     internal variables
c
      double precision ddot,ek,t,wk,wkm
      double precision anorm,s,dasum,sm,ynorm
      integer info,j,k,kb,kp1,l
c
c
c     compute 1-norm of a
c
      anorm = 0.0d0
      do 10 j = 1, n
         anorm = dmax1(anorm,dasum(n,a(1,j),1))
   10 continue
c
c     factor
c
      call dgefa(a,lda,n,ipvt,info)
c
c     rcond = 1/(norm(a)*(estimate of norm(inverse(a)))) .
c     estimate = norm(z)/norm(y) where  a*z = y  and  trans(a)*y = e .
c     trans(a)  is the transpose of a .  the components of  e  are
c     chosen to cause maximum local growth in the elements of w  where
c     trans(u)*w = e .  the vectors are frequently rescaled to avoid
c     overflow.
c
c     solve trans(u)*w = e
c
      ek = 1.0d0
      do 20 j = 1, n
         z(j) = 0.0d0
   20 continue
      do 100 k = 1, n
         if (z(k) .ne. 0.0d0) ek = dsign(ek,-z(k))
         if (dabs(ek-z(k)) .le. dabs(a(k,k))) go to 30
            s = dabs(a(k,k))/dabs(ek-z(k))
            call dscal(n,s,z,1)
            ek = s*ek
   30    continue
         wk = ek - z(k)
         wkm = -ek - z(k)
         s = dabs(wk)
         sm = dabs(wkm)
         if (a(k,k) .eq. 0.0d0) go to 40
            wk = wk/a(k,k)
            wkm = wkm/a(k,k)
         go to 50
   40    continue
            wk = 1.0d0
            wkm = 1.0d0
   50    continue
         kp1 = k + 1
         if (kp1 .gt. n) go to 90
            do 60 j = kp1, n
               sm = sm + dabs(z(j)+wkm*a(k,j))
               z(j) = z(j) + wk*a(k,j)
               s = s + dabs(z(j))
   60       continue
            if (s .ge. sm) go to 80
               t = wkm - wk
               wk = wkm
               do 70 j = kp1, n
                  z(j) = z(j) + t*a(k,j)
   70          continue
   80       continue
   90    continue
         z(k) = wk
  100 continue
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
c
c     solve trans(l)*y = w
c
      do 120 kb = 1, n
         k = n + 1 - kb
         if (k .lt. n) z(k) = z(k) + ddot(n-k,a(k+1,k),1,z(k+1),1)
         if (dabs(z(k)) .le. 1.0d0) go to 110
            s = 1.0d0/dabs(z(k))
            call dscal(n,s,z,1)
  110    continue
         l = ipvt(k)
         t = z(l)
         z(l) = z(k)
         z(k) = t
  120 continue
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
c
      ynorm = 1.0d0
c
c     solve l*v = y
c
      do 140 k = 1, n
         l = ipvt(k)
         t = z(l)
         z(l) = z(k)
         z(k) = t
         if (k .lt. n) call daxpy(n-k,t,a(k+1,k),1,z(k+1),1)
         if (dabs(z(k)) .le. 1.0d0) go to 130
            s = 1.0d0/dabs(z(k))
            call dscal(n,s,z,1)
            ynorm = s*ynorm
  130    continue
  140 continue
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
      ynorm = s*ynorm
c
c     solve  u*z = v
c
      do 160 kb = 1, n
         k = n + 1 - kb
         if (dabs(z(k)) .le. dabs(a(k,k))) go to 150
            s = dabs(a(k,k))/dabs(z(k))
            call dscal(n,s,z,1)
            ynorm = s*ynorm
  150    continue
         if (a(k,k) .ne. 0.0d0) z(k) = z(k)/a(k,k)
         if (a(k,k) .eq. 0.0d0) z(k) = 1.0d0
         t = -z(k)
         call daxpy(k-1,t,a(1,k),1,z(1),1)
  160 continue
c     make znorm = 1.0
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
      ynorm = s*ynorm
c
      if (anorm .ne. 0.0d0) rcond = ynorm/anorm
      if (anorm .eq. 0.0d0) rcond = 0.0d0
      return
      end
      subroutine dgefa(a,lda,n,ipvt,info)
      integer lda,n,ipvt(1),info
      double precision a(lda,1)
c
c     dgefa factors a double precision matrix by gaussian elimination.
c
c     dgefa is usually called by dgeco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for dgeco) = (1 + 9/n)*(time for dgefa) .
c
c     on entry
c
c        a       double precision(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that dgesl or dgedi will divide by zero
c                     if called.  use  rcond  in dgeco for a reliable
c                     indication of singularity.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,dscal,idamax
c
c     internal variables
c
      double precision t
      integer idamax,j,k,kp1,l,nm1
c
c
c     gaussian elimination with partial pivoting
c
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
c
c        find l = pivot index
c
         l = idamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
c
c        zero pivot implies this column already triangularized
c
         if (a(l,k) .eq. 0.0d0) go to 40
c
c           interchange if necessary
c
            if (l .eq. k) go to 10
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
   10       continue
c
c           compute multipliers
c
            t = -1.0d0/a(k,k)
            call dscal(n-k,t,a(k+1,k),1)
c
c           row elimination with column indexing
c
            do 30 j = kp1, n
               t = a(l,j)
               if (l .eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
               call daxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (a(n,n) .eq. 0.0d0) info = n
      return
      end
      subroutine  dscal(n,da,dx,incx)
c
c     scales a vector by a constant.
c     uses unrolled loops for increment equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified to correct problem with negative increment, 8/21/90.
c
      double precision da,dx(1)
      integer i,incx,ix,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      do 10 i = 1,n
        dx(ix) = da*dx(ix)
        ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dx(i) = da*dx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
   50 continue
      return
      end
      integer function idamax(n,dx,incx)
c
c     finds the index of element having max. absolute value.
c     jack dongarra, linpack, 3/11/78.
c     modified to correct problem with negative increment, 8/21/90.
c
      double precision dx(1),dmax
      integer i,incx,ix,n
c
      idamax = 0
      if( n .lt. 1 ) return
      idamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      dmax = dabs(dx(ix))
      ix = ix + incx
      do 10 i = 2,n
         if(dabs(dx(ix)).le.dmax) go to 5
         idamax = i
         dmax = dabs(dx(ix))
    5    ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 dmax = dabs(dx(1))
      do 30 i = 2,n
         if(dabs(dx(i)).le.dmax) go to 30
         idamax = i
         dmax = dabs(dx(i))
   30 continue
      return
      end
      double precision function dasum(n,dx,incx)
c
c     takes the sum of the absolute values.
c     uses unrolled loops for increment equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified to correct problem with negative increment, 8/21/90.
c
      double precision dx(1),dtemp
      integer i,incx,ix,m,mp1,n
c
      dasum = 0.0d0
      dtemp = 0.0d0
      if(n.le.0)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      do 10 i = 1,n
        dtemp = dtemp + dabs(dx(ix))
        ix = ix + incx
   10 continue
      dasum = dtemp
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,6)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dtemp + dabs(dx(i))
   30 continue
      if( n .lt. 6 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,6
        dtemp = dtemp + dabs(dx(i)) + dabs(dx(i + 1)) + dabs(dx(i + 2))
     *  + dabs(dx(i + 3)) + dabs(dx(i + 4)) + dabs(dx(i + 5))
   50 continue
   60 dasum = dtemp
      return
      end
C
      SUBROUTINE ZERO(RM,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION RM(*)
      DO 1 I=1,N
         RM(I)=0.D0
    1 CONTINUE
      RETURN
      END
C
C *** INVERSION DE LA MATRICE [A] ***
C
      SUBROUTINE INVERS(A,B,N,IR,IC)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(N,N),B(N,N)
      DIMENSION W1(N),W2(N)
      DIMENSION IR(N),IC(N)
c
      CALL LUPIV(A,N,N,IR,IC,IPER,W2)
c
      DO 1 I=1,N
         DO 2 J=1,N
            W2(J)=0.
   2     CONTINUE
         W2(I)=1.
         CALL SYSTEM1(N,N,A,W1,W2,IR,IC)
         DO 3 J=1,N
C        WRITE(6,*) ' W1',W1(J)
            B(J,I)=W1(J)
   3     CONTINUE
   1  CONTINUE
c
      RETURN
      END
C
C *** INVERSION DE LA MATRICE [A] ***
C
      SUBROUTINE PRODMV(A,B,C,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(N,N),B(N),C(N)
c
      DO 1 I=1,N
         som=0.D0
         DO 2 J=1,N
            som=som+A(I,J)*B(J)
   2     CONTINUE
         C(I)=som
   1  CONTINUE
c
      RETURN
      END
