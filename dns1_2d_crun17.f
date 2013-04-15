! --------------------------------------------------------- MARS 2001 .
! ------------------------------------------------------ M.MARQUILLIE .
!
! =====================================================================
! =====================================================================
! =  SIMULATION SPATIALE D UN ECOULEMENT DE TYPE COUCHE LIMITE        =
! =                   EN GEOMETRIE COMPLEXE                           =
! =====================================================================
! =====================================================================
!
      program dns2d
      implicit double precision (a-h,o-z)
!
! ---------------------------------------------------------------------
!     NX   : NOMBRE DE POINTS EN DIFFERENCES FINIES SUIVANT X
!     NY : NOMBRE DE POINTS DE COLLOCATIONS DE CHEBYSCHEV SUIVANT Y
!     YSL  : COEFFICIENT DE REPARTITION DES POINTS DE COLLOCATIONS 
!            DE CHEBYSCHEV SUIVANT Y
!     YMAX : BORNE SUPERIEURE DE LA DISCRETISATION SUIVANT Y
!     PASX : VALEUR DE CHAQUE SUBDIVISION SUIVANT X
!     PAST : VALEUR DE L INCREMENTATION EN TEMPS
!     XA   : ABCISSE X A L ENTREE DU DOMAINE
!     XB   : ABCISSE X A LA SORTIE DU DOMAINE
!     REY  : NOMBRE DE REYNOLDS
! ---------------------------------------------------------------------
!
      include 'par.f'
      parameter (delta=1.7207678d0)
      parameter (nm=2)
!
!     NM=1 : FFT
!     NM=2 : MULTIPLICATION MATRICIELLE
!
! --------------------------------------------------------------------
!     DX   : MATRICE DERIVEE PREMIERE EN DIFFERENCES FINIES SUIVANT X
!     DDX  : MATRICE DERIVEE SECONDE EN DIFFERENCES FINIES SUIVANT X
!     DDDDX: MATRICE DERIVEE QUATRIEME EN DIFFERENCES FINIES SUIVANT X
!     Y    : ORDONNEE DES POINTS DE COLLOCATION SUIVANT Y
!     DY   : MATRICE DE COLLOCATION DE LA DERIVEE PREMIERE SUIVANT Y
!     DDY  : MATRICE DE COLLOCATION DE LA DERIVEE SECONDE SUIVANT Y
! --------------------------------------------------------------------
!
      dimension dx(7,8),ddx(0:4,5),ddddx(10,10)
      dimension y(ny),dy(ny,ny),ddy(ny,ny) 
      dimension ddxn(7,5),ddyn(ny,ny),dyn(ny,ny)
      dimension dy1(ny,ny),ddy1(ny,ny)
!
! --------------------------------------------------------------------
!     F     : MATRICE REPRESENTANT LE SECOND MEMBRE 
!             DU PROBLEME DE TYPE HELMHOLTZ
!     U,V,P : MATRICES REPRESENTANT LES SOLUTIONS
!             DU PROBLEME DE TYPE HELMHOLTZ
! --------------------------------------------------------------------
!
      dimension wsave1(3*ny+15),tr(ny,2)
      dimension xfc(20),xmesh(nx,ny),ymesh(nx,ny)
      dimension ps(nx,ny),wz(nx,ny)
      dimension umean(nx,ny),vmean(nx,ny)
      dimension uabs(nx,ny),vabs(nx,ny)
      dimension tauu(nx,ny),tauv(nx,ny)
      dimension pmean(nx,ny),tmpp(nx,ny)
      dimension tmpum(nx,ny),tmpvm(nx,ny)
      dimension tmpua(nx,ny),tmpva(nx,ny)
      dimension tmptu(nx,ny),tmptv(nx,ny)
      dimension tmpmin(nx)
!
      dimension etha(nx),bosse(nwal)
      dimension dxetha(nx),ddxetha(nx)
!
      dimension u(nx,ny,2)
      dimension dxu(nx,ny,2),dyu(nx,ny,2) 
      dimension dxdyu(nx,ny,2),ddyu(nx,ny,2)
!
      dimension v(nx,ny,2)
      dimension dxv(nx,ny,2),dyv(nx,ny,2) 
      dimension dxdyv(nx,ny,2),ddyv(nx,ny,2)
!
      dimension p(nx,ny,2)
      dimension dxp(nx,ny,2),dyp(nx,ny,2) 
      dimension dxdyp(nx,ny,2),ddyp(nx,ny,2)
!
      dimension pl(nx,ny)
      dimension dxpl(nx,ny),dypl(nx,ny) 
! 
      dimension phi(nx,ny)
      dimension dxphi(nx,ny),dyphi(nx,ny)
      dimension dxdyphi(nx,ny),ddyphi(nx,ny)
!
      dimension ubl(ny),vbl(ny),yubl(ny)
      dimension clpxen(ny),clpxsn(ny),clpypn(nx)
      dimension clpyin(nx) 
      dimension cl1(nx,4),cl2(ny,4) ! paroi,infini,entree,sortie
!
      dimension clubxe(ny)
! 
      dimension f(ny-2,nx-2),fp(ny-2,nx-2)
      dimension fu(ny-2,nx-2),fv(ny-2,nx-2)
      dimension phel1(nx-2,ny-2),phel2(nx-2,ny-2)
      dimension phel3(nx-2,ny-2),phel4(nx-2,ny-2)
      dimension phel5(nx-2,ny-2),pterm(ny-2)
      dimension uhel1(nx-2,ny-2),uhel2(nx-2,ny-2)
      dimension uhel3(nx-2,ny-2),uhel4(nx-2,ny-2)
      dimension uhel5(nx-2,ny-2),uterm(ny-2)
!
      dimension q(ny-2,ny-2),b(ny-2,ny-2)
      dimension xmat(ny-2,ny-2)
      dimension xwr(ny-2),xwi(ny-2),xmodw(ny-2)
      dimension irq(ny-2),icq(ny-2),workq(ny-2) 
!
      dimension qn(ny-2,ny-2),bn(ny-2,ny-2)
      dimension xmatn(ny-2,ny-2) 
      dimension xwrn(ny-2),xwin(ny-2),xmodwn(ny-2) 
      dimension irqn(ny-2),icqn(ny-2),workqn(ny-2)
!
      dimension q1(ny-2,ny-2),q2(ny-2,ny-2)
      dimension dxu1(nx,ny),dxv1(nx,ny)
!
      external blas
!
! =====================================================================
! =                   PARAMETRES DE LA SIMULATION                     =
! =====================================================================
!
      pi=4.d0*datan(1.d0)
!
      write(6,*) 'PRE-CALCULS :'
!
!-> lecture des parametres dans le fichier param.dat et 
!   lecture du graphe de la bosse dans le fichier bossen
!
      open(42,file='param.dat',status='old',form='formatted')
      read(42,*) rey,ysl,past,ymax
      close(42)
!
      open(11,file='bossen',status='old',form='formatted')
      read(11,*)pasx,bosse 
      close(11)
!
      pasx=0.2d0
!      do i=1,nwal
!      bosse(i)=0.d0 
!      enddo
!
!-> calcul de la longueur du domaine, du Reynolds a la sortie
!   et du Reynolds au sommet de la bosse 
!
      xa=0. 
      xb=xa+(nx-1)*pasx 
      reyb=dsqrt(((xb*(delta**2)*rey)+(rey**2))) 
      xbo=xa+(nin-1)*pasx
      reybo=dsqrt(((xbo*(delta**2)*rey)+(rey**2)))
!
!-> ecriture des parametres dans le fichier fort.60
!
      write(60,*)
      write(60,*)'donnees :'
      write(60,*)'    ysl          = ', ysl
      write(60,*)'    ymax         = ', ymax
      write(60,*)'    nx, ny       = ',nx,ny
      write(60,*)'    xa et xb     = ',xa,xb
      write(60,*)'    reyb         = ',reyb
      write(60,*)'    reybosse     = ',reybo
      write(60,*)'    pasx et past = ',pasx,past
      write(60,*)'    rey inflow   = ',rey
      write(60,*)'    nin et nout  = ',nin,nout
!
! =====================================================================
! =           CALCUL DES MATRICES DERIVEES PAR CHEBYSHEV              =
! =           ET CALCUL DES POINTS DE MAILLAGE DU DOMAINE             =
! =====================================================================
!
      write(6,*)'    -> maillage, matrices derivees.'
!
!-> calcul du maillage et des matrices derivees 
!
      call dermat (ysl,pasx,y,dx,ddx,ddddx,dy,ddy,nx,ny,ymax,tr)
!
      do j=1,ny
         write(60,'(i3,e17.8)') j,y(j)
      enddo
!
!-> coefficients de derivation en differences finies pour 
!   les conditions de Neuman 
! 
      do l=2,5 
      ddxn(1,l)=ddx(0,l)-ddx(0,1)*dx(1,l)/dx(1,1) 
      ddxn(2,l)=ddx(1,l)-ddx(1,1)*dx(1,l)/dx(1,1) 
      ddxn(3,l)=ddx(2,l)-ddx(2,1)*dx(1,l)/dx(1,1) 
      enddo 
      do l=1,5 
      ddxn(4,l)=ddx(2,l) 
      enddo 
      do l=1,4 
      ddxn(5,l)=ddx(2,l)-ddx(2,5)*dx(7,l)/dx(7,5) 
      ddxn(6,l)=ddx(3,l)-ddx(3,5)*dx(7,l)/dx(7,5) 
      ddxn(7,l)=ddx(4,l)-ddx(4,5)*dx(7,l)/dx(7,5) 
      enddo  
!
!-> matrices derivees de Chebyshev pour les conditions de Neuman
! 
      do j=1,ny 
      do l=1,ny 
      ddyn(j,l)=0.d0 
      dyn(j,l)=0.d0 
      enddo 
      enddo 
! 
      do j=1,ny 
      do l=2,ny-1 
      aux1=dy(1,1)*dy(ny,ny)-dy(1,ny)*dy(ny,1) 
      aux2=dy(1,ny)*dy(ny,l)-dy(ny,ny)*dy(1,l) 
      aux3=dy(ny,1)*dy(1,l)-dy(1,1)*dy(ny,l) 
      ddyn(j,l)=ddy(j,l)+ddy(j,1)*aux2/aux1+ddy(j,ny)*aux3/aux1 
      dyn(j,l)=dy(j,l)+dy(j,1)*aux2/aux1+dy(j,ny)*aux3/aux1 
      enddo 
      enddo  
!
!-> transposition des matrices derivees Chebyshev pour l'utilisation
!   de la subroutine MATMUL
!  
      do j=1,ny
      do l=1,ny
      dy1(l,j)=dy(j,l)
      ddy1(l,j)=ddy(j,l)
      enddo
      enddo 
!
! =====================================================================
! =         PARAMETRES DU CHANGEMENT DE VARIABLE : ETHA               =
! =====================================================================
!
      write(6,*)'    -> changement de variable.'
      etha=0.d0
      dxetha=0.d0
      ddxetha=0.d0
!
!-> affectation du graphe de la bosse dans etha
! 
      do i=1,nin 
      etha(i)=0.d0 
      enddo
! 
      do i=nin+1,nout-2 
      etha(i)=bosse(i-nin)
      enddo
!
      do i=nout-1,nx 
      etha(i)=0.d0
      enddo
      
!      do i=1,nx
!         xi=(i-1)*pasx
!         etha(i)=10.d0*(-0.216d0/cosh(0.28d0*(xi-2.5d0*10.d0)))
!      enddo
! 
!-> calcul des derivees de etha et localisation du sommet 
!   de la bosse
!
      ethamax=0.d0 
      do i=2,nx-1 
         call derdxcl(i,dx,etha,nx,dxe) 
!         call derddxcl(i,ddx,etha,nx,ddxe) 
         dxetha(i)=dxe 
!         ddxetha(i)=ddxe 
         if (etha(i).ge.ethamax) then 
         nvit=i 
         ethamax=etha(i) 
         endif 
      enddo
      do i=2,nx-1 
         call derdxcl(i,dx,dxetha,nx,ddxe) 
         ddxetha(i)=ddxe 
      enddo
!
      do i=1,nx  
!      dxetha(i)=0.d0
!      ddxetha(i)=0.d0 
      enddo
      write(60,*)'Sommet de la bosse :',nvit
      close(60) 
      do i=1,nx
      write(65,'(4es17.8)')(i-1)*pasx,etha(i),dxetha(i),ddxetha(i)
      enddo
      close(65)
!      stop
!
!-> ecriture du maillage dans le fichier mesh2d
!
!      goto 1201
      open(600,file='mesh2dbump1')
      open(601,file='mesh2dbump2')
      do j=1,ny,1
         do i=1,nx,1
            xmesh(i,j)=xa+(i-1)*pasx
            ymesh(i,j)=y(j)+etha(i)
            write(600,'(2es17.8)')xmesh(i,j),ymesh(i,j)
         enddo
         write(600,*)
      enddo
      do i=1,nx,1
         do j=1,ny,1
            xmesh(i,j)=xa+(i-1)*pasx
            ymesh(i,j)=y(j)+etha(i)
            write(601,'(2es17.8)')xmesh(i,j),ymesh(i,j)
         enddo
         write(601,*)
      enddo
      close(600)
      close(601)
!      stop
! 1201 continue
!
! =====================================================================
! =         INITIALISATION DES VARIABLES ET CALCUL DU PROFIL          =
! =                       IMPOSE A L ENTREE                           =
! =====================================================================
!
      write(6,*)'    -> initialisation des variables.'
!
      do j=1,ny 
         do i=1,nx 
            u(i,j,1)=0.d0 
            u(i,j,2)=0.d0 
            v(i,j,1)=0.d0 
            v(i,j,2)=0.d0 
            p(i,j,1)=0.d0 
            p(i,j,2)=0.d0
            phi(i,j)=0.d0
!========= AJOUT CARO
            umean(i,j)=0.d0
            vmean(i,j)=0.d0
            uabs(i,j)=0.d0
            vabs(i,j)=0.d0
            tauu(i,j)=0.d0
            tauv(i,j)=0.d0
            pmean(i,j)=0.d0
            tmpum(i,j)=0.d0
            tmpvm(i,j)=0.d0
            tmpua(i,j)=0.d0
            tmpva(i,j)=0.d0
            tmptu(i,j)=0.d0
            tmptv(i,j)=0.d0
            tmpp(i,j)=0.d0

!===================
         enddo
      enddo
! 
      call blas1(y,ubl,vbl,ny)
      do j=1,ny
      clubxe(j)=ubl(j)
!      write(33,'(2es17.8)')y(j),clubxe(j)
      enddo 
 
!      do i=1,nx
!      write(34,'(2es17.8)')(i-1)*pasx,etha(i)
!      enddo

!
! =====================================================================
! =        DETERMINATION DES VALEURS PROPRES ET VECTEURS PROPRES      =
! =                 DE LA MATRICE A DIAGONALISER                      =
! =====================================================================
!
      write(6,*)'    -> initialisation du probleme de Helmholtz.'
!
      do j=1,ny-2
         do i=1,ny-2
            xmat(i,j)=0.d0
            b(i,j)=0.d0 
            xmatn(i,j)=0.d0 
            bn(i,j)=0.d0
         enddo
      enddo
!
      do j=2,ny-1
         b(j-1,j-1)=1.0d0 
         bn(j-1,j-1)=1.0d0
         do i=2,ny-1
            xmat(j-1,i-1)=ddy(j,i) 
            xmatn(j-1,i-1)=ddyn(j,i)
         enddo
      enddo
!
      matz=1
      call rgg(ny-2,ny-2,xmat,b,xwr,xwi,xmodw,matz,q,ierr) 
      call rgg(ny-2,ny-2,xmatn,bn,xwrn,xwin,xmodwn,matz,qn,ierr)
!
! --------------------------------------------------------------------
! xwr et xwi contiennent respectivement les parties relles et 
! imaginaires des valeurs propres.                 
! xmodw representent les modules des valeurs propres.         
! matz = 0 donnent les valeurs propres                        
! sinon donnent les valeurs et les vecteurs propres. 
! Q contient les vecteurs propres en colonnes.                
! ierr different de zero indique que le calcul de la j-eme    
! valeur propre a depasse les 30*(NY) iterations.  
! --------------------------------------------------------------------    
!
      if (ierr.ne.0) then
         write(6,*) 'ierr = ',ierr
         stop
      endif
!
      do j=1,ny-2
         do i=1,ny-2
            b(i,j)=q(i,j) 
            bn(i,j)=qn(i,j)
         enddo
      enddo
!
      call lupiv(q,ny-2,ny-2,irq,icq,iper,workq) 
      call lupiv(qn,ny-2,ny-2,irqn,icqn,ipern,workqn)
!
      ntest=0
      test=abs(xwrn(1)) 
      do l=1,ny-2 
         if (abs(xwrn(l)).le.test) then 
            test=abs(xwrn(l)) 
            nprp=l 
         endif     
         write(61,*)xwrn(l),xwin(l) 
         if (xwin(l).ne.0.d0) then 
            write(*,*) 'probleme sur les valeurs propres' 
            ntest=1
         endif 
      enddo 
      write(61,*)nprp,test 
      xwrn(nprp)=0.d0
      close(61)
      if (ntest.eq.1) stop
!
!-> calcul de l inverse de q et qn
!
      do j=1,ny-2
         do k=1,ny-2
            q1(j,k)=0.d0
         enddo
      enddo
      do k=1,ny-2
         q1(k,k)=1.d0
      enddo
      call resolq(ny-2,ny-2,q,q2,q1,irq,icq)
      do j=1,ny-2
         do k=1,ny-2
            q(j,k)=q2(j,k)
         enddo
      enddo
      call resolq(ny-2,ny-2,qn,q2,q1,irqn,icqn)
      do j=1,ny-2
         do k=1,ny-2
            qn(j,k)=q2(j,k)
         enddo
      enddo
!
! =====================================================================
! =            CALCUL DE LA MATRICE BANDE DE HELMHOLTZ                =
! =====================================================================
!
      write(6,*)'    -> matrice bande de Helmholtz.'
!
!-> phel* : matrice bande pour la pression p
!
      sigma=0.D0
!
      do ml=1,ny-2
!
         phel3(1,ml)=ddxn(2,2)+(xwrn(ml)/xmodwn(ml))-sigma
         phel4(1,ml)=ddxn(2,3)
         phel5(1,ml)=ddxn(2,4)
         palpha     =ddxn(2,5)
!
         phel2(2,ml)=ddxn(3,2)
         phel3(2,ml)=ddxn(3,3)+(xwrn(ml)/xmodwn(ml))-sigma
         phel4(2,ml)=ddxn(3,4)
         phel5(2,ml)=ddxn(3,5)
!
         do l=3,nx-4
            phel1(l,ml)=ddxn(4,1)
            phel2(l,ml)=ddxn(4,2)
            phel3(l,ml)=ddxn(4,3)+(xwrn(ml)/xmodwn(ml))-sigma
            phel4(l,ml)=ddxn(4,4)
            phel5(l,ml)=ddxn(4,5)
         enddo
!
         phel1(nx-3,ml)=ddxn(5,1)
         phel2(nx-3,ml)=ddxn(5,2)
         phel3(nx-3,ml)=ddxn(5,3)+(xwrn(ml)/xmodwn(ml))-sigma
         phel4(nx-3,ml)=ddxn(5,4)
!
         pterm(ml)     =ddxn(6,1)
         phel1(nx-2,ml)=ddxn(6,2)
         phel2(nx-2,ml)=ddxn(6,3)
         phel3(nx-2,ml)=ddxn(6,4)+(xwrn(ml)/xmodwn(ml))-sigma 
! 
	if (ml.eq.nprp) then  
         phel1(nx-3,ml)=ddx(2,1) 
         phel2(nx-3,ml)=ddx(2,2) 
         phel3(nx-3,ml)=ddx(2,3)+(xwrn(ml)/xmodwn(ml))-sigma 
         phel4(nx-3,ml)=ddx(2,4) 
! 
         pterm(ml)     =ddx(3,1) 
         phel1(nx-2,ml)=ddx(3,2) 
         phel2(nx-2,ml)=ddx(3,3) 
         phel3(nx-2,ml)=ddx(3,4)+(xwrn(ml)/xmodwn(ml))-sigma
	endif
!
         call pentalu(phel1(1,ml),phel2(1,ml),phel3(1,ml),phel4(1,ml),
     &phel5(1,ml),palpha,pterm(ml))
!
      enddo
!
!-> uhel* : matrice bande pour les vitesses u et v
!
      sigma=1.5d0*rey/past 
!
      do ml=1,ny-2
!
         uhel3(1,ml)=ddx(1,2)+(xwr(ml)/xmodw(ml))-sigma
         uhel4(1,ml)=ddx(1,3)
         uhel5(1,ml)=ddx(1,4)
         ualpha     =ddx(1,5)
!
         uhel2(2,ml)=ddx(2,2)
         uhel3(2,ml)=ddx(2,3)+(xwr(ml)/xmodw(ml))-sigma
         uhel4(2,ml)=ddx(2,4)
         uhel5(2,ml)=ddx(2,5)
!
         do l=3,nx-4
            uhel1(l,ml)=ddx(2,1)
            uhel2(l,ml)=ddx(2,2)
            uhel3(l,ml)=ddx(2,3)+(xwr(ml)/xmodw(ml))-sigma
            uhel4(l,ml)=ddx(2,4)
            uhel5(l,ml)=ddx(2,5)
         enddo
!
         uhel1(nx-3,ml)=ddx(2,1)
         uhel2(nx-3,ml)=ddx(2,2)
         uhel3(nx-3,ml)=ddx(2,3)+(xwr(ml)/xmodw(ml))-sigma
         uhel4(nx-3,ml)=ddx(2,4)
!
         uterm(ml)     =ddx(3,1)
         uhel1(nx-2,ml)=ddx(3,2)
         uhel2(nx-2,ml)=ddx(3,3)
         uhel3(nx-2,ml)=ddx(3,4)+(xwr(ml)/xmodw(ml))-sigma
!
         call pentalu(uhel1(1,ml),uhel2(1,ml),uhel3(1,ml),uhel4(1,ml),
     &uhel5(1,ml),ualpha,uterm(ml))
!
      enddo
!
! =====================================================================
! =              INITIALISATION : FFT, iact, ianc.                    =
! =              LECTURE ET CALCUL DES DERIVEES POUR                  =
! =              LE REDEMARAGE D UN CALCUL.                           =
! =              PERTURBATION INITIALE.                               =
! =====================================================================
!
      write(6,*)'    -> initialisations du calcul.'
! 
!-> initialisation pour les fft 
!
      call vcosti(ny,wsave1) 
!
!-> initialisation de iact et ianc
!
      iact=1 
      ianc=2
!
!-> lecture du fichier de redemarage
!
      do k=1,1
         read(801)iact,ianc,u,v,p,cl1,cl2,phi
      enddo
!      write(921)iact,ianc,u,v,p,cl1,cl2,phi
!      stop
!
!      nx1=1000
!      do l=1,3
!      read(801)iact,ianc,u(1:nx1,1:ny,1:2),v(1:nx1,1:ny,1:2),
!     &p(1:nx1,1:ny,1:2),cl1(1:nx1,1:4),cl2(1:ny,1:4),
!     &phi(1:nx1,1:ny)
!      enddo
!      do i=nx1+1,nx
!         u(i,1:ny,1:2)=u(nx1,1:ny,1:2)
!         v(i,1:ny,1:2)=v(nx1,1:ny,1:2)
!         p(i,1:ny,1:2)=p(nx1,1:ny,1:2)
!         cl1(i,1:4)=cl1(nx1,1:4)
!         phi(i,1:ny)=phi(nx1,1:ny)
!      enddo
!      print*,maxval(u),maxval(v)
!
!-> perturbation initiale
!
!      sigx=0.2d0
!      sigy=0.1d0
!      am=5.d-1
!      xo=60.d0
!      yo=1.5d0
!      do j=2,ny-1
!      do i=2,nx-1
!      xi=float(i-1)*pasx
!      yi=y(j)
!      a1=-am*(yi-yo)*dexp(-(xi-xo)*(xi-xo)/(2.d0*sigx*sigx)
!     &-(yi-yo)*(yi-yo)/(2.d0*sigy*sigy))
!      u(i,j,ianc)=u(i,j,iact)+a1
!      a2=am*sigy*sigy/(sigx*sigx)*(xi-xo)
!     &*dexp(-(xi-xo)*(xi-xo)/(2.d0*sigx*sigx)
!     &-(yi-yo)*(yi-yo)/(2.d0*sigy*sigy))
!      v(i,j,ianc)=v(i,j,iact)+a2
!      enddo
!      enddo
!      stop
!
!-> calcul des derivees de la solution initiale
!
      call derx_1(dx,p(1,1,iact),ny-2,nx-2,dxp(1,1,iact)) 
      call dery(tr,dy1,ddy1,p(1,1,iact),nx,ny,wsave1,
     &dyp(1,1,iact),ddyp(1,1,iact),2,nm) 
      call derx_1(dx,dyp(1,1,iact),ny-2,nx-2,dxdyp(1,1,iact)) 
! 
      call derx_1(dx,u(1,1,iact),ny-2,nx-2,dxu(1,1,iact))
      call dery(tr,dy1,ddy1,u(1,1,iact),nx,ny,wsave1,
     &dyu(1,1,iact),ddyu(1,1,iact),2,nm)  
      call derx_1(dx,dyu(1,1,iact),ny-2,nx-2,dxdyu(1,1,iact)) 
! 
      call derx_1(dx,v(1,1,iact),ny-2,nx-2,dxv(1,1,iact)) 
      call dery(tr,dy1,ddy1,v(1,1,iact),nx,ny,wsave1,
     &dyv(1,1,iact),ddyv(1,1,iact),2,nm)  
      call derx_1(dx,dyv(1,1,iact),ny-2,nx-2,dxdyv(1,1,iact)) 
! 
      call derx_1(dx,p(1,1,ianc),ny-2,nx-2,dxp(1,1,ianc))
      call dery(tr,dy1,ddy1,p(1,1,ianc),nx,ny,wsave1,
     &dyp(1,1,ianc),ddyp(1,1,ianc),2,nm)  
      call derx_1(dx,dyp(1,1,ianc),ny-2,nx-2,dxdyp(1,1,ianc))  
!  
      call derx_1(dx,u(1,1,ianc),ny-2,nx-2,dxu(1,1,ianc)) 
      call dery(tr,dy1,ddy1,u(1,1,ianc),nx,ny,wsave1,
     &dyu(1,1,ianc),ddyu(1,1,ianc),2,nm)  
      call derx_1(dx,dyu(1,1,ianc),ny-2,nx-2,dxdyu(1,1,ianc)) 
! 
      call derx_1(dx,v(1,1,ianc),ny-2,nx-2,dxv(1,1,ianc))
      call dery(tr,dy1,ddy1,v(1,1,ianc),nx,ny,wsave1,
     &dyv(1,1,ianc),ddyv(1,1,ianc),2,nm)  
      call derx_1(dx,dyv(1,1,ianc),ny-2,nx-2,dxdyv(1,1,ianc)) 
 

!      do j=1,ny
!         do i=1,nx
!            edx=dxetha(i)
!            dxu1(i,j)=dxu(i,j,ianc)-edx*dyu(i,j,ianc)
!            dxv1(i,j)=dxv(i,j,ianc)-edx*dyv(i,j,ianc)
!         enddo
!      enddo
!      open(980,file='500cavity1.dat',access='direct',
!     &recl=8+((nx)*(ny))*8*8)
!      write(980,rec=1)nx,ny,xmesh,ymesh,
!     &u(1:nx,1:ny,ianc:ianc),v(1:nx,1:ny,ianc:ianc),
!     &dxu1(1:nx,1:ny),dxv1(1:nx,1:ny),
!     &dyu(1:nx,1:ny,ianc:ianc),dyv(1:nx,1:ny,ianc:ianc)
!      close(980)
!      stop
!
! =====================================================================
! =                        CALCUL GENERAL                             =
! =====================================================================
!

      write(6,*)
      write(6,*) 'CALCUL GENERAL :'
!
!-> debut de la boucle en temps
!   ite : nombres d iterations en temps
!

      nitei=0
      nitef=20000
      do 10000 ite=nitei,nitef 
! 
!-> conditions limites a l entree et a l infini
! 
      do j=1,ny 
         u(1,j,iact)=clubxe(j) 
      enddo 
      do i=1,nx 
         u(i,ny,iact)=1.d0 
      enddo
!####################################################################
!-> CONTROL: 
!
      xite=float(ite)*past
!     blowing-suction sinus function:
!      xl=2.d0

!      do i=1,10
!         ampl=-0.19d0
!         xfc(iact)=dcos(2.d0*pi*float(iact)*pasx/xl)
!         xfc(i)=dcos(2.d0*pi*float(i)*pasx/xl) 
!         xfc(i)=dcos(2.d0*pi*float(ite)*past)
!         u(nvit-10+i,1,iact)=ampl*xfc(i)
!         v(nvit-10+i,1,iact)=ampl*xfc(i)
!      enddo
!         write(210,*) xite,u(nvit-10+1,1,iact),v(nvit-10+1,1,iact)

!     Step function
      icontb=200
      icontc=1500
c      write(*,*) ite
      if ((ite.ge.icontb).and.(ite.le.icontc)) then
c      write(*,*) 'ITE=',ite
         ampl=0.1d0
         do i=1,10
            xfc(i)=1
            u(nvit+10-i,1,iact)=ampl*xfc(i) !*xfc(i)/sqrt(2.d0)
            v(nvit+10-i,1,iact)=ampl*xfc(i) !*xfc(i)/sqrt(2.d0)
         enddo
         else
         do i=1,10
            u(nvit+10-i,1,iact)=0.d0
            v(nvit+10-i,1,iact)=0.d0
         enddo
      end if
      write(210,*) xite,u(nvit+10-1,1,iact),v(nvit+10-1,1,iact)

!#######################################################################
! 
!-> condition d advection a la sortie
!
      som=0.d0 
      do j=1,ny-1
      som=som+((u(nx,j+1,ianc)+u(nx,j,ianc))/2.d0)* 
     &(y(j+1)-y(j)) 
      yc=y(j+1) 
      if (u(nx,j,ianc).gt.(0.5d0)) goto 855 
      enddo 
 855  continue 
      uc=som/yc 
!      uc=0.1
      do j=2,ny-1 
         a1=4.d0*u(nx,j,ianc)-cl2(j,2)
         a2=past*2.d0 
         a3=2.d0*dxu(nx,j,ianc)*uc 
         a4=dxu(nx,j,iact)*uc 
         u(nx,j,iact)=(a1-a2*(a3-a4))/3.d0 
!
         a1=4.d0*v(nx,j,ianc)-cl2(j,4) 
         a2=past*2.d0 
         a3=2.d0*dxv(nx,j,ianc)*uc 
         a4=dxv(nx,j,iact)*uc 
         v(nx,j,iact)=(a1-a2*(a3-a4))/3.d0
      enddo

      u(nx,ny,iact)=1.d0 
      v(nx,ny,iact)=0.d0 
!
! =====================================================================
! =             CALCUL DES SECONDS MEMBRES EXPLICITES                 =
! =====================================================================
!
      do j=2,ny-1
         do i=2,nx-1
!
            p2dy=dyp(i,j,ianc)
            p2ddy=ddyp(i,j,ianc)
            p2dxdy=dxdyp(i,j,ianc)
            p1dy=dyp(i,j,iact)
            p1ddy=ddyp(i,j,iact)
            p1dxdy=dxdyp(i,j,iact)
!
            u2dx=dxu(i,j,ianc)
            u2dy=dyu(i,j,ianc)
            u2ddy=ddyu(i,j,ianc)
            u2dxdy=dxdyu(i,j,ianc)
            u1dx=dxu(i,j,iact)
            u1dy=dyu(i,j,iact)
            u1ddy=ddyu(i,j,iact)
            u1dxdy=dxdyu(i,j,iact)
!
            v2dx=dxv(i,j,ianc)
            v2dy=dyv(i,j,ianc)
            v2ddy=ddyv(i,j,ianc)
            v2dxdy=dxdyv(i,j,ianc)
            v1dx=dxv(i,j,iact)
            v1dy=dyv(i,j,iact)
            v1ddy=ddyv(i,j,iact)
            v1dxdy=dxdyv(i,j,iact)
!
            dtu=0.5d0*(4.d0*u(i,j,ianc)-u(i,j,iact))/past
            dtv=0.5d0*(4.d0*v(i,j,ianc)-v(i,j,iact))/past
! 
            edx=dxetha(i) 
            eddx=ddxetha(i)
!
      poids1a=4.d0*(u2dx*v2dy-u2dy*v2dx)-2.d0*(u1dx*v1dy-u1dy*v1dx)
      poids1b= 
     &2.d0*(eddx*p2dy+2.d0*edx*p2dxdy-(edx*edx)*p2ddy)
     &-1.d0*(eddx*p1dy+2.d0*edx*p1dxdy-(edx*edx)*p1ddy)
!
      poids2a=-2.d0*(u(i,j,ianc)*u2dx+v(i,j,ianc)*u2dy- 
     &u(i,j,ianc)*edx*u2dy)+ 
     &(u(i,j,iact)*u1dx+v(i,j,iact)*u1dy- 
     &u(i,j,iact)*edx*u1dy)
!  
      poids2b=-2.d0*((eddx*u2dy+2.d0*edx*u2dxdy-(edx*edx)*u2ddy)/rey)+ 
     &((eddx*u1dy+2.d0*edx*u1dxdy-(edx*edx)*u1ddy)/rey) 
!
      poids3a=-2.d0*(u(i,j,ianc)*v2dx+v(i,j,ianc)*v2dy- 
     &u(i,j,ianc)*edx*v2dy)+ 
     &(u(i,j,iact)*v1dx+v(i,j,iact)*v1dy- 
     &u(i,j,iact)*edx*v1dy) 
! 
      poids3b=-2.d0*((eddx*v2dy+2.d0*edx*v2dxdy-(edx*edx)*v2ddy)/rey)+ 
     &((eddx*v1dy+2.d0*edx*v1dxdy-(edx*edx)*v1ddy)/rey) 
!
!-> calcul du second membre explicite de p
!
      fp(j-1,i-1)=(poids1b+poids1a)
!
!-> calcul du second membre explicite de u
!
      fu(j-1,i-1)=-rey*(dtu+poids2b+poids2a)
!
!-> calcul du second membre explicite de v
!
      fv(j-1,i-1)=-rey*(dtv+poids3b+poids3a)
!
      enddo
      enddo
!
! =====================================================================
! =       CONDITIONS DE NEUMAN POUR LA PAROI ET L'INFINI              =
! =====================================================================
!
      do i=2,nx-1
      edx=dxetha(i)
      eddx=ddxetha(i)
      t_pres=edx*(2.d0*(dxp(i,ny,ianc)-edx*dyp(i,ny,ianc))-
     &(dxp(i,ny,iact)-edx*dyp(i,ny,iact)))
      t_nl=2.d0*u(i,ny,ianc)*edx*(dyv(i,ny,ianc) 
     &-edx*dyu(i,ny,ianc))-
     &cl1(i,2)*edx*(dyv(i,ny,iact)-edx*dyu(i,ny,iact))
      t_dif=-2.d0*(eddx*dyv(i,ny,ianc)+dxdyu(i,ny,ianc)+edx*
     &dxdyv(i,ny,ianc))/rey+
     &(eddx*dyv(i,ny,iact)+dxdyu(i,ny,iact) 
     &+edx*dxdyv(i,ny,iact))/rey
      t_dt=0.5d0*(-edx*(3.d0*u(i,ny,iact) 
     &-4.d0*u(i,ny,ianc)+cl1(i,2))+
     &(3.d0*v(i,ny,iact)-4.d0*v(i,ny,ianc)+cl1(i,4)))/past
      clpyin(i)=t_pres+t_nl+t_dif-t_dt
!
      t_pres=edx*(2.d0*(dxp(i,1,ianc)-edx*dyp(i,1,ianc))-
     &(dxp(i,1,iact)-edx*dyp(i,1,iact)))
      t_dif=-2.d0*(eddx*dyv(i,1,ianc)+dxdyu(i,1,ianc)+edx*
     &dxdyv(i,1,ianc))/rey+
     &(eddx*dyv(i,1,iact)+dxdyu(i,1,iact)+edx*dxdyv(i,1,iact))/rey
      t_dt=0.5d0*(-edx*(3.d0*u(i,1,iact)-4.d0*u(i,1,ianc)+cl1(i,1))+
     &(3.d0*v(i,1,iact)-4.d0*v(i,1,ianc)+cl1(i,3)))/past
      clpypn(i)=t_pres+t_dif-t_dt
      enddo
!
! =====================================================================
! =       CONDITIONS DE NEUMAN POUR L ENTREE ET LA SORTIE             =
! =====================================================================
!
      do j=2,ny-1
      clpxen(j)=2.d0*(-dxdyv(1,j,ianc)/REY+
     &ddyu(1,j,ianc)/REY-u(1,j,ianc)*dxu(1,j,ianc))-
     &(-dxdyv(1,j,iact)/REY+ddyu(1,j,iact)/REY 
     &-cl2(j,1)*dxu(1,j,iact))
     &-0.5d0*(3.d0*u(1,j,iact)-4.d0*u(1,j,ianc)+cl2(j,1))/past
!
      clpxsn(j)=2.d0*(-dxdyv(nx,j,ianc)/REY
     &+ddyu(nx,j,ianc)/REY
     &-u(nx,j,ianc)*dxu(nx,j,ianc)-
     &v(nx,j,ianc)*dyu(nx,j,ianc))-
     &(-dxdyv(nx,j,iact)/REY+ddyu(nx,j,iact)/REY 
     &-cl2(j,2)*dxu(nx,j,iact)-
     &cl2(j,4)*dyu(nx,j,iact))
     &-0.5d0*(3.d0*u(nx,j,iact)-
     &4.d0*u(nx,j,ianc)+cl2(j,2))/past
      enddo
!
! =====================================================================
! =               RESOLUTION DU PROBLEME DE HELMHOLTZ                 =
! =                      LAPLACIEN ( P* ) = F                         =
! =====================================================================
!
!-> calcul du second membre du probleme de helmholtz
!
      do j=2,ny-1
         do i=2,nx-1
            f(j-1,i-1)=fp(j-1,i-1)
         enddo
      enddo
! 
      aux1=dy(1,1)*dy(ny,ny)-dy(1,ny)*dy(ny,1) 
      do j=2,ny-1  
         f(j-1,1)=f(j-1,1)-ddx(1,1)*clpxen(j)/dx(1,1)
         f(j-1,2)=f(j-1,2)-ddx(2,1)*clpxen(j)/dx(1,1) 
         f(j-1,nx-3)=f(j-1,nx-3)-ddx(2,5)*clpxsn(j)/dx(7,5)
         f(j-1,nx-2)=f(j-1,nx-2)-ddx(3,5)*clpxsn(j)/dx(7,5)
         do i=2,nx-1
         aux2=ddy(j,1)*dy(ny,ny)-
     &ddy(j,ny)*dy(ny,1)
         aux3=ddy(j,ny)*dy(1,1)-ddy(j,1)*dy(1,ny)
         f(j-1,i-1)=f(j-1,i-1)-clpypn(i)*aux2/aux1-clpyin(i)*aux3/aux1
         enddo
      enddo
!
!-> resolution du probleme de helmholtz
! 
      call helmholtz(phel1,phel2,phel3,phel4,phel5,palpha,pterm, 
     &bn,qn,irqn,icqn,f,pl)
! 
      call neumx(nx,ny,pl,dx,clpxen,clpxsn) 
      call neumy(nx,ny,pl,dy,clpypn,clpyin)
!
      call derx_1(dx,pl,ny-2,nx-2,dxpl) 
      call dery(tr,dy1,ddy1,pl,nx,ny,wsave1,
     &dypl,ddypl,1,nm)
!
! =====================================================================
! =               RESOLUTION DU PROBLEME DE HELMHOLTZ                 =
! =                [ LAPLACIEN - SIGMA ] ( U* ) = F                   =
! =====================================================================
!
!-> calcul du second membre du probleme de helmholtz
!
      do j=2,ny-1
         do i=2,nx-1
            f(j-1,i-1)=rey*(dxpl(i,j)-dxetha(i)*dypl(i,j))
     &                 +fu(j-1,i-1) 
         enddo
      enddo
!
      do j=2,ny-1
         f(j-1,1)=f(j-1,1)-ddx(1,1)*u(1,j,iact)
         f(j-1,2)=f(j-1,2)-ddx(2,1)*u(1,j,iact)
         f(j-1,nx-3)=f(j-1,nx-3)-ddx(2,5)*u(nx,j,iact)
         f(j-1,nx-2)=f(j-1,nx-2)-ddx(3,5)*u(nx,j,iact)
         do i=2,nx-1
         f(j-1,i-1)=f(j-1,i-1)-ddy(j,1)*u(i,1,iact) 
     &   -ddy(j,ny)*u(i,ny,iact)
         enddo
      enddo
!
!-> resolution du probleme de helmholtz
!
      call helmholtz(uhel1,uhel2,uhel3,uhel4,uhel5,ualpha,uterm,
     &b,q,irq,icq,f,u(1,1,iact))
!
! =====================================================================
! =               RESOLUTION DU PROBLEME DE HELMHOLTZ                 =
! =                [ LAPLACIEN - SIGMA ] ( V* ) = F                   =
! =====================================================================
!
!-> calcul du second membre du probleme de helmholtz
!
      do j=2,ny-1
         do i=2,nx-1
            f(j-1,i-1)=rey*dypl(i,j)+fv(j-1,i-1)
         enddo
      enddo
!
      do j=2,ny-1
         f(j-1,1)=f(j-1,1)-ddx(1,1)*v(1,j,iact)
         f(j-1,2)=f(j-1,2)-ddx(2,1)*v(1,j,iact)
         f(j-1,nx-3)=f(j-1,nx-3)-ddx(2,5)*v(nx,j,iact)
         f(j-1,nx-2)=f(j-1,nx-2)-ddx(3,5)*v(nx,j,iact)
         do i=2,nx-1
         f(j-1,i-1)=f(j-1,i-1)-ddy(j,1)*v(i,1,iact) 
     &   -ddy(j,ny)*v(i,ny,iact)
         enddo
      enddo
!
!-> resolution du probleme de helmholtz
!
      call helmholtz(uhel1,uhel2,uhel3,uhel4,uhel5,ualpha,uterm,
     &b,q,irq,icq,f,v(1,1,iact))
!
! =====================================================================
! =               RESOLUTION DU PROBLEME DE HELMHOLTZ                 =
! =                      LAPLACIEN ( PHI ) =  F                       =
! =====================================================================
!
      call derx_1(dx,u(1,1,iact),ny-2,nx-2,dxu(1,1,iact)) 
      call dery(tr,dy1,ddy1,v(1,1,iact),nx,ny,wsave1,
     &dyv(1,1,iact),ddyv(1,1,iact),1,nm)  
      call dery(tr,dy1,ddy1,u(1,1,iact),nx,ny,wsave1,
     &dyu(1,1,iact),ddyu(1,1,iact),1,nm)
!
!-> debut du procede iteratif
!
      ites=0
      do k=1,50
!
!-> calcul du second membre du probleme de helmholtz
!
      call derx_1(dx,phi,ny-2,nx-2,dxphi)
      call dery(tr,dy1,ddy1,phi,nx,ny,wsave1,
     &dyphi,ddyphi,2,nm)
      call derx_1(dx,dyphi,ny-2,nx-2,dxdyphi) 
!
      do j=2,ny-1
         do i=2,nx-1 
         edx=dxetha(i)
         eddx=ddxetha(i)
         phidy=dyphi(i,j)
         phiddy=ddyphi(i,j)
         phidxdy=dxdyphi(i,j)
         f(j-1,i-1)=1.5d0*(dxu(i,j,iact)+dyv(i,j,iact)-
     &   dxetha(i)*dyu(i,j,iact))/past
     &   +(eddx*phidy+2.d0*edx*phidxdy-(edx*edx)*phiddy)
         enddo
      enddo
!
      do i=2,nx-1
         edx=dxetha(i)
         clpyin(i)=edx*dxphi(i,ny)/(1.d0+edx*edx)
         clpypn(i)=edx*dxphi(i,1)/(1.d0+edx*edx)
      enddo
!
      do j=1,ny
         clpxen(j)=0.d0
         clpxsn(j)=0.d0
      enddo
!
      aux1=dy(1,1)*dy(ny,ny)-dy(1,ny)*dy(ny,1)  
      do j=2,ny-1 
         f(j-1,1)=f(j-1,1)-ddx(1,1)*clpxen(j)/dx(1,1)
         f(j-1,2)=f(j-1,2)-ddx(2,1)*clpxen(j)/dx(1,1) 
         f(j-1,nx-3)=f(j-1,nx-3)-ddx(2,5)*clpxsn(j)/dx(7,5)
         f(j-1,nx-2)=f(j-1,nx-2)-ddx(3,5)*clpxsn(j)/dx(7,5)
         do i=2,nx-1
         aux2=ddy(j,1)*dy(ny,ny)-
     &ddy(j,ny)*dy(ny,1)
         aux3=ddy(j,ny)*dy(1,1)-ddy(j,1)*dy(1,ny)
         f(j-1,i-1)=f(j-1,i-1)-clpypn(i)*aux2/aux1-clpyin(i)*aux3/aux1
         enddo
      enddo
!
!-> resolution du probleme de helmholtz
!
      call helmholtz(phel1,phel2,phel3,phel4,phel5,palpha,pterm, 
     &bn,qn,irqn,icqn,f,p(1,1,iact)) 
!
      call neumx(nx,ny,p(1,1,iact),dx,clpxen,clpxsn) 
      call neumy(nx,ny,p(1,1,iact),dy,clpypn,clpyin)
!
!-> test de convergence
!
      pmax=0.d0
      do j=1,ny
         do i=1,nx
         if (dabs(p(i,j,iact)).ge.pmax) pmax=dabs(p(i,j,iact))
         enddo
      enddo
!
      testmax=0.d0
      do j=1,ny
         do i=1,nx
         test=dabs((p(i,j,iact)-phi(i,j))/pmax)
         if (test.ge.testmax) testmax=test
         enddo
      enddo
!
      do j=1,ny
         do i=1,nx
         phi(i,j)=p(i,j,iact)
         enddo
      enddo
!
      prec=past*past
      ites=ites+1
!      write(*,*)ites,pmax,testmax
      if (testmax.le.prec) exit
!
!-> fin du procede iteratif
!
      enddo
!      write(55,*)ite,ites,testmax
!
! =====================================================================
! =                CALCUL DE U, V, P ET DE LEUR DERIVEES              =
! =                           AU TEMPS N+1                            =
! =====================================================================
!
!-> calcul de u, v et p
!
      call derx_1(dx,p(1,1,iact),ny-2,nx-2,dxp(1,1,iact)) 
      call dery(tr,dy1,ddy1,p(1,1,iact),nx,ny,wsave1,
     &dyp(1,1,iact),ddyp(1,1,iact),1,nm) 
!
      do j=2,ny-1  
      do i=2,nx-1
      u(i,j,iact)=u(i,j,iact)-2.d0*past*(dxp(i,j,iact)-dxetha(i)*
     &dyp(i,j,iact))/3.d0
      v(i,j,iact)=v(i,j,iact)-2.d0*past*dyp(i,j,iact)/3.d0
      enddo 
      enddo 
      do j=1,ny
      do i=1,nx
      p(i,j,iact)=pl(i,j)+p(i,j,iact)
      enddo
      enddo 
!
!-> calcul des derivees de u,v et p
!
      call derx_1(dx,p(1,1,iact),ny-2,nx-2,dxp(1,1,iact))
      call dery(tr,dy1,ddy1,p(1,1,iact),nx,ny,wsave1,
     &dyp(1,1,iact),ddyp(1,1,iact),2,nm)
      call derx_1(dx,dyp(1,1,iact),ny-2,nx-2,dxdyp(1,1,iact)) 
!
      call derx_1(dx,u(1,1,iact),ny-2,nx-2,dxu(1,1,iact))
      call dery(tr,dy1,ddy1,u(1,1,iact),nx,ny,wsave1,
     &dyu(1,1,iact),ddyu(1,1,iact),2,nm) 
      call derx_1(dx,dyu(1,1,iact),ny-2,nx-2,dxdyu(1,1,iact)) 
!  
      call derx_1(dx,v(1,1,iact),ny-2,nx-2,dxv(1,1,iact)) 
      call dery(tr,dy1,ddy1,v(1,1,iact),nx,ny,wsave1,
     &dyv(1,1,iact),ddyv(1,1,iact),2,nm)  
      call derx_1(dx,dyv(1,1,iact),ny-2,nx-2,dxdyv(1,1,iact))
!
! =====================================================================
! =                  ENREGISTREMENTS DES RESULTATS                    =
! =====================================================================
!
!-> SENSORS
c      xite=float(ite)*past
c      write(200,*) xite,dyu(int(30.D0/pasx),1,iact)

!
!-> 
!
!      nu?=nvit+5 !(x,y)=(25.8,1.957)
      nu1=nvit+10 !(x,y)=(26.8,1.83) (nvit=125)
      nu2=nvit+16 !(x,y)=(28,1.60)
      nu3=nvit+21 !(x,y)=(29,1.356)
      nu4=nvit+26 !(x,y)=(30,1.089)
      nu5=nvit+31 !(x,y)=(31,0.82)
      nu6=nvit+36 !(x,y)=(32,0.57)
      nu7=nvit+41 !(x,y)=(32.8,0.396)
      nu8=nvit+46 !(x,y)=(34,0.18)
      nu9=nvit+56 !(x,y)=(36,0.01)
      nu10=198 !(x,y)=(39.4,0.00088)
      nu11=241 !(x,y)=(48,0.000357)
      xite=float(ite)*past
      write(301,108)xite,u(nu1,2,iact),v(nu1,2,iact),
     &dyu(nu1,2,iact),dxv(nu1,2,iact),p(nu1,2,iact)
      write(302,108)xite,u(nu2,2,iact),v(nu2,2,iact),
     &dyu(nu2,2,iact),dxv(nu2,2,iact),p(nu2,2,iact)
      write(303,108)xite,u(nu3,2,iact),v(nu3,2,iact),
     &dyu(nu3,2,iact),dxv(nu3,2,iact),p(nu3,2,iact)
      write(304,108)xite,u(nu4,2,iact),v(nu4,2,iact),
     &dyu(nu4,2,iact),dxv(nu4,2,iact),p(nu4,2,iact)
      write(305,108)xite,u(nu5,2,iact),v(nu5,2,iact),
     &dyu(nu5,2,iact),dxv(nu5,2,iact),p(nu5,2,iact)
      write(306,108)xite,u(nu6,2,iact),v(nu6,2,iact),
     &dyu(nu6,2,iact),dxv(nu6,2,iact),p(nu6,2,iact)
      write(307,108)xite,u(nu7,2,iact),v(nu7,2,iact),
     &dyu(nu7,2,iact),dxv(nu7,2,iact),p(nu7,2,iact)
      write(308,108)xite,u(nu8,2,iact),v(nu8,2,iact),
     &dyu(nu8,2,iact),dxv(nu8,2,iact),p(nu8,2,iact)
      write(309,108)xite,u(nu9,2,iact),v(nu9,2,iact),
     &dyu(nu9,2,iact),dxv(nu9,2,iact),p(nu9,2,iact)
      write(310,108)xite,u(nu10,2,iact),v(nu10,2,iact),
     &dyu(nu10,2,iact),dxv(nu10,2,iact),p(nu10,2,iact)
      write(311,108)xite,u(nu11,2,iact),v(nu11,2,iact),
     &dyu(nu11,2,iact),dxv(nu11,2,iact),p(nu11,2,iact)

 108  format(6e23.14)

!
!-> evolution temporelle de la vitesse u en differentes locations
!
c      nu1=nx/6           
c      nu2=2*nx/6
c      nu3=3*nx/6
c      nu4=4*nx/6
c      nu5=5*nx/6
c      xite=float(ite)*past
c      write(92,108)xite,u(nu1,2,iact),u(nu1,5,iact),
c     &u(nu1,17,iact),u(nu1,45,iact),u(nu1,60,iact)
c      write(93,108)xite,u(nu2,2,iact),u(nu2,5,iact),
c     &u(nu2,17,iact),u(nu2,45,iact),u(nu2,60,iact)
c      write(94,108)xite,u(nu3,2,iact),u(nu3,5,iact),
c     &u(nu3,17,iact),u(nu3,45,iact),u(nu3,60,iact)
c      write(95,108)xite,u(nu4,2,iact),u(nu4,5,iact),
c     &u(nu4,17,iact),u(nu4,45,iact),u(nu4,60,iact)
c      write(96,108)xite,u(nu5,2,iact),u(nu5,5,iact),
c     &u(nu5,17,iact),u(nu5,45,iact),u(nu5,60,iact)


!
!-> test de convergence pour u et v
!
      umax=0.d0
      do j=2,ny-1 
         do i=2,nx-1 
            if (dabs(u(i,j,iact)).ge.umax) umax=dabs(u(i,j,iact)) 
         enddo
      enddo
      conv1=0.d0 
      conv2=0.d0
      do j=2,ny-1 
         do i=2,nx-1
            convu=dabs((u(i,j,iact)-u(i,j,ianc))/(umax*past)) 
            convv=dabs((v(i,j,iact)-v(i,j,ianc))/(umax*past)) 
            if (convu.ge.conv1) conv1=convu 
            if (convv.ge.conv2) conv2=convv 
         enddo
      enddo
!
!-> calcul de la divergence de u et v
!
      xmax1=0.d0
      xmax2=0.d0
      do i=2,nx-1
         xdiv=0.d0
         do j=2,ny-1
         test1=dabs(dxu(i,j,iact)-dxetha(i)* 
     &   dyu(i,j,iact)+dyv(i,j,iact))/umax
         if (test1.gt.xmax1) xmax1=test1
         enddo
!         if (ite-10000*(ite/10000).eq.0) then
!            xi=xa+(i-1)*pasx
!            write(3350,*) xi,xmax1
!         endif
      enddo
!
      write(64,67)ite,conv1,conv2,xmax1 
      write(*,'(2i6,3es21.14)')ite,ites,conv1,conv2,xmax1 
!      if (ite==nitef) 
!     &write(*,'(2i6,3es21.14)')ite,ites,conv1,conv2,xmax1
 67   format(i6,3es17.8)
      if (conv1.eq.(0.d0)) stop
!
! =====================================================================
! =                      PASSAGE : N+1 -> N                           =
! =                                N -> N-1                           =
! =====================================================================
! 
	do i=1,nx 
           cl1(i,1)=u(i,1,ianc) 
           cl1(i,2)=u(i,ny,ianc) 
           cl1(i,3)=v(i,1,ianc) 
           cl1(i,4)=v(i,ny,ianc) 
	enddo   
	do j=1,ny 
           cl2(j,1)=u(1,j,ianc) 
           cl2(j,2)=u(nx,j,ianc) 
           cl2(j,3)=v(1,j,ianc) 
           cl2(j,4)=v(nx,j,ianc) 
	enddo  
	iaux1=iact 
	iaux2=ianc 
	iact=iaux2 
	ianc=iaux1
!
! =====================================================================
! =            ENREGISTREMENT DU FICHIER DE REDEMARAGE                =
! =====================================================================
! 
        ipas1=50000
        if (ite-ipas1*(ite/ipas1).eq.0.and.ite>0) then
           if (ite.eq.ipas1) ifich=0
           if (ite>ipas1) ifich=ifich+1
           write(902+ifich)iact,ianc,u,v,p,cl1,cl2,phi
           close(902+ifich)
        endif
!        ipas2=5000
        ipas2=500
        if (ite-ipas2*(ite/ipas2).eq.0.and.ite>0) then
!
           do i=1,nx
              do j=1,ny
                 som=0.d0
                 do k=1,j-1
                    som=som+(u(i,k,ianc)+u(i,k+1,ianc))
     &                   *(y(k+1)-y(k))*0.5d0
                 enddo
                 ps(i,j)=som
              enddo
           enddo

           do i=1,nx
              do j=1,ny
                 wz(i,j)=dxv(i,j,ianc)-dxetha(i)*dyv(i,j,ianc)
     &                -dyu(i,j,ianc)
!                 umean(i,j)=umean(i,j)+u(i,j,ianc)
              enddo
           enddo
           call stockage(ps,xmesh,ymesh,nx,ny,ite,past,'Ps')
           call stockage(wz,xmesh,ymesh,nx,ny,ite,past,'Wz')
        endif
!================================================
!     CONVERGENCE MEAN ET VALEUR ABSOLUE DE U ET V
!================================================
           do i=1,nx
              do j=1,ny
                 tmpum(i,j)=tmpum(i,j)+u(i,j,ianc)
                 tmpvm(i,j)=tmpvm(i,j)+v(i,j,ianc)
                 tmpua(i,j)=tmpua(i,j)+abs(u(i,j,ianc))
                 tmpva(i,j)=tmpva(i,j)+abs(v(i,j,ianc))
                 tmptu(i,j)=tmptu(i,j)+abs(dyu(i,j,ianc))
                 tmptv(i,j)=tmptv(i,j)+abs(dxv(i,j,ianc))
                 tmpp(i,j)=tmpp(i,j)+p(i,j,ianc)

              enddo
           enddo
           
           if (ite-ipas2*(ite/ipas2).eq.0.and.ite>0) then
              do i=1,nx
                 do j=1,ny

                    umean(i,j)=tmpum(i,j)/float(ite)
                    vmean(i,j)=tmpvm(i,j)/float(ite)
                    uabs(i,j)=tmpua(i,j)/float(ite)
                    vabs(i,j)=tmpva(i,j)/float(ite)
                    tauu(i,j)=tmptu(i,j)/float(ite)
                    tauv(i,j)=tmptv(i,j)/float(ite)
                    pmean(i,j)=tmpp(i,j)/float(ite)
         write(400,'(es17.8,i3,2es17.8)') i*pasx,j,umean(i,j),uabs(i,j)
c     write(60,'(i3,e17.8)') j,y(j)
c                    write(400,'(i4,i3,4es17.8)') i,j,umean(i,j),vmean(i,
c     &                   j),uabs(i,j),vabs(i,j)
                 enddo
                 write(400,*)
              enddo
              
              call stockage(umean,xmesh,ymesh,nx,ny,ite,past,'Um') !range (-0.2:1.2)
              call stockage(vmean,xmesh,ymesh,nx,ny,ite,past,'Vm') !range (-0.06:0.01)
              call stockage(tauu,xmesh,ymesh,nx,ny,ite,past,'tu') !range (0:0.7)
              call stockage(tauv,xmesh,ymesh,nx,ny,ite,past,'tv')  !range (0:0.004)
              call stockage(pmean,xmesh,ymesh,nx,ny,ite,past,'Pm')!range (-0.06:0.003)
              call stockage(uabs,xmesh,ymesh,nx,ny,ite,past,'Ua') !range (0:1.2)
              call stockage(vabs,xmesh,ymesh,nx,ny,ite,past,'Va') !range (0:0.01)
c           endif
!=================================================================


!================================================
!     Calcul lf
!================================================
!!!!!! USING changement signe of Umean
           do i=nvit-20,nx-10
              if (umean(i,2).le.0.d0) then
                 nx1=i                 
                 exit
              endif              
           enddo           
           do i=nx1+1,nx-10             
              if (umean(i,2).ge.0.d0) then
                 nx2=i                 
                 exit                 
              endif            
           enddo
           xlf=(nx2-nx1)*pasx
        write(100,'(2es17.8,i4,i4)') xite,xlf,nx1,nx2

!!!!!! USING tauu (research of first and second minimum)
        do i=nvit+12,nx-10
           if (tauu(i-1,2).ge.tauu(i,2).and.tauu(i+1,2).ge.tauu(i,2)) 
     &          then
              nxmin1=i
              exit
           endif
        enddo             
        do i=nxmin1+1,nx-10
           if (tauu(i-1,2).ge.tauu(i,2).and.tauu(i+1,2).ge.tauu(i,2)) 
     &          then
              nxmin2=i                                  
              exit                 
           endif            
        enddo
        xlf=(nxmin2-(nvit+11))*pasx
!     write(110,'(2es17.8)') float(i),tauu(i,2)

        write(110,'(2es17.8,i4,i4)') xite,xlf,nxmin1,nxmin2
      endif

!
!-> fin de la boucle en temps
!

10000 continue 
!============ STOCKAGE MEAN |Tauu| =========
        do i=nvit+12,nx-10
           write(120,'(2es17.8)') float(i),tauu(i,2)
        enddo


      do j=1,ny
         do i=1,nx
            edx=dxetha(i)
            dxu1(i,j)=dxu(i,j,ianc)-edx*dyu(i,j,ianc)
            dxv1(i,j)=dxv(i,j,ianc)-edx*dyv(i,j,ianc)
         enddo
      enddo

!      open(980,file='450ds.dat',access='direct',
!     &recl=8+((nx)*(ny))*8*8)
!      write(980,rec=1)nx,ny,xmesh,ymesh,
!     &u(1:nx,1:ny,ianc:ianc),v(1:nx,1:ny,ianc:ianc),
!     &dxu1(1:nx,1:ny),dxv1(1:nx,1:ny),
!     &dyu(1:nx,1:ny,ianc:ianc),dyv(1:nx,1:ny,ianc:ianc)
!      close(980)
      stop

      end

!
! =====================================================================
! =                             STOCKAGE                              =
! =====================================================================
!
      subroutine stockage(u,x,y,nx,ny,ite,past,abrev)
      implicit double precision(a-h,o-z)
      dimension u(nx,ny),x(nx,ny),y(nx,ny)
      character(6) num
      character(8) nom
      character(2) abrev

      write(num,'(i6)')100000+int(ite*past)
      write(nom,'(a2,a6)')abrev,num(1:6)
      nlenr=2*4+nx*ny*8*3
      open(unit=904,file=nom,access='direct',recl=nlenr)
      write(904,rec=1)nx,ny,x,y,u
      close(904)


      end
