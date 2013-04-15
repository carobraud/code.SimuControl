! =====================================================================
! =           RESOLUTION D UN PROBLEME DE TYPE HELMHOLTZ              =
! =                 [ DD - sigma ] [ U ] = [ F ]                      =
! =       DD represente l operateur Laplacien en bidimensionnel       =
! =====================================================================
!
      subroutine helmholtz(diag1,diag2,diag3,diag4,diag5,dalpha,dterm,
     &b,q,irq,icq,f,u)
      implicit double precision (a-h,o-z)
      include 'par.f'
!
! --------------------------------------------------------------------
!     NX : NOMBRE DE POINTS EN DIFFERENCES FINIES SUIVANT X
!     NY : NOMBRE DE POINTS DE COLLOCATIONS DE CHEBYSCHEV SUIVANT Y
!     B    : MATRICE SERVANT A  R G G (MATRICE IDENTITE) ET 
!            AU CALCUL DE LA MATRICE U
!     F    : MATRICE REPRESENTANT LE SECOND MEMBRE 
!            DU PROBLEME DE TYPE HELMHOLTZ
!     Q    : MATRICE CONTENANT LES VECTEURS PROPRES CALCULES PAR RGG
!     R    : MATRICE RESULTANT DU SYSTEME [Q][R]=[F]
!     THETA: MATRICE RESULTANT DU SYSTEME [A][THETA]=[R]
!     U    : MATRICE REPRESENTANT LA SOLUTION
!            DU PROBLEME DE TYPE HELMHOLTZ
!     IRA, ICA, IRQ, ICQ, W1, W2, WORK1 : VECTEURS SERVANT A LA
!            RESOLUTION D UN SYSTEME DE TYPE [A] (X) = (B) PAR LA
!            METHODE DE DECOMPOSITION L U .
! --------------------------------------------------------------------
!
      dimension b(ny-2,ny-2),work1(nx-2)
      dimension q(ny-2,ny-2),r(ny-2,nx-2)
      dimension theta(ny-2,nx-2),u1(ny-2,nx-2)
      dimension u(nx,ny),f(ny-2,nx-2)
      dimension irq(ny-2),icq(ny-2)
      dimension diag1(nx-2,ny-2),diag2(nx-2,ny-2)
      dimension diag3(nx-2,ny-2),diag4(nx-2,ny-2)
      dimension diag5(nx-2,ny-2),dterm(ny-2)
!
! =====================================================================
! = 1.) DETERMINATION DE LA MATRICE [R] :  [R]=[Q][F]                 =
! =====================================================================
!
!      r=matmul(q,f)
      call dgemm('n','n',ny-2,nx-2,ny-2,1.d0,q,ny-2,f,ny-2,0.d0,r,ny-2)
!
! =====================================================================
! = 2.) DETERMINATION DE LA MATRICE [THETA] : [A][THETA]=[R]          =
! =====================================================================
! 
      do k=1,ny-2 
         theta(k,1)=r(k,1) 
      enddo 
      do k=1,ny-2 
         theta(k,2)=r(k,2)-theta(k,1)*diag2(2,k) 
      enddo 
      do i=3,nx-3 
         do k=1,ny-2 
            theta(k,i)=r(k,i)-diag1(i,k)*theta(k,i-2) 
     &-diag2(i,k)*theta(k,i-1) 
         enddo 
      enddo 
      do k=1,ny-2 
         theta(k,nx-2)=r(k,nx-2)-diag1(nx-2,k)*
     &theta(k,nx-4) 
     &-diag2(nx-2,k)*theta(k,nx-3)-dterm(k)*theta(k,nx-5) 
      enddo 
! 
      do k=1,ny-2 
         theta(k,nx-2)=theta(k,nx-2)/diag3(nx-2,k) 
      enddo 
      do k=1,ny-2 
         theta(k,nx-3)=(theta(k,nx-3) 
     &-diag4(nx-3,k)*theta(k,nx-2))/diag3(nx-3,k) 
      enddo 
      do i=nx-4,2,-1 
         do k=1,ny-2 
            theta(k,i)=(theta(k,i)-diag5(i,k)*theta(k,i+2) 
     &-diag4(i,k)*theta(k,i+1))/diag3(i,k) 
         enddo 
      enddo 
      do k=1,ny-2 
         theta(k,1)=(theta(k,1)-diag5(1,k)*theta(k,3) 
     &-diag4(1,k)*theta(k,2)-dalpha*theta(k,4))/diag3(1,k) 
      enddo
!
! =====================================================================
! = 3.)  DETERMINATION DE LA MATRICE [U] :  [U]=[Q][THETA]            =
! =====================================================================
!
!      u1=matmul(b,theta)
      call dgemm('n','n',ny-2,nx-2,ny-2,1.d0,b,ny-2,theta,ny-2,0.d0,
     &u1,ny-2)
!     
      do j=2,ny-1
         do i=2,nx-1
            u(i,j)=u1(j-1,i-1)
         enddo
      enddo
!
      return
      end
!
! =====================================================================
! =                 RESOLUTION DU PROBLEME [A](X)=(B)                 =
! =                  OU [A] EST QUASI-PENTADIAGONALE                  =
! =====================================================================
!
      subroutine pentalu(p,f,c,d,e,alpha,beta)
      implicit double precision (a-h,o-z)
      include 'par.f'
      dimension p(nx-2),f(nx-2),c(nx-2)
      dimension d(nx-2),e(nx-2)
!
! -> etape 1
!
      f(2)=f(2)/c(1)
      c(2)=c(2)-d(1)*f(2)
      d(2)=d(2)-e(1)*f(2)
      e(2)=e(2)-alpha*f(2)
!
! -> etape 2
!
      i=3
      p(i)=p(i)/c(i-2)
      f(i)=(f(i)-p(i)*d(i-2))/c(i-1)
      c(i)=c(i)-d(i-1)*f(i)-e(i-2)*p(i)
      d(i)=d(i)-e(i-1)*f(i)-alpha*p(i)
!
      do 50 i=4,nx-3
         p(i)=p(i)/c(i-2)
         f(i)=(f(i)-p(i)*d(i-2))/c(i-1)
         c(i)=c(i)-e(i-2)*p(i)-d(i-1)*f(i)
         d(i)=d(i)-e(i-1)*f(i)
 50   continue
!
! -> etape 3
!
      beta=beta/c(nx-5)
      p(nx-2)=p(nx-2)-beta*d(nx-5)
      f(nx-2)=f(nx-2)-beta*e(nx-5)
      p(nx-2)=p(nx-2)/c(nx-4)
      f(nx-2)=(f(nx-2)-p(nx-2)*d(nx-4))/c(nx-3)
      c(nx-2)=c(nx-2)-e(nx-4)*p(nx-2)-
     &d(nx-3)*f(nx-2)
!
      return
      end
!
! =====================================================================
!
      subroutine pentares(p,f,c,d,e,alpha,beta,x,b)
      implicit double precision (a-h,o-z)
      include 'par.f'
      dimension p(nx-2),f(nx-2),c(nx-2),d(nx-2)
      dimension e(nx-2),x(nx-2),b(nx-2)
!
! -> etape 4
!
      x(1)=b(1)
      x(2)=b(2)-x(1)*f(2)
!
! -> etape 5
!
      do 100 i=3,nx-3
         x(i)=b(i)-p(i)*x(i-2)-f(i)*x(i-1)
 100  continue
!
      i=nx-2
      x(i)=b(i)-p(i)*x(i-2)-f(i)*x(i-1)-beta*x(i-3)
!
! -> etape 6
!
      x(nx-2)=x(nx-2)/c(nx-2)
      x(nx-3)=(x(nx-3)-d(nx-3)*x(nx-2))/c(nx-3)
!
! -> etape 7
!
      do 150 i=nx-4,2,-1
         x(i)=(x(i)-e(i)*x(i+2)-d(i)*x(i+1))/c(i)
 150  continue
!
      i=1
      x(i)=(x(i)-e(i)*x(i+2)-d(i)*x(i+1)-alpha*x(i+3))/c(i)
!
      return
      end
!
! =====================================================================
! =  CALCUL DES DERIVEES COLLOCATION-CHEBYSHEV : CHOIX DE LA METHODE  =
! =         NM=1 : FFT                                                =
! =         NM=2 : MULTIPLICATION MATRICIELLE                         =
! =====================================================================
!
      subroutine dery(tz,dy,ddy,a,nx,ny,wsave,dya,ddya,nc1,nm)
      implicit double precision (a-h,o-z)
      dimension dy(ny,ny),ddy(ny,ny)
      dimension a(nx,ny),dya(nx,ny)
      dimension ddya(nx,ny)
      dimension tz(ny,2)
      dimension wsave(3*(ny)+15),xt(nx,ny-1)
      if (nm.eq.1) then
      call deryfft(tz,a,nx,ny,wsave,dya,ddya,nc1)
      elseif (nm.eq.2) then
         if (nc1.eq.1) then
            call dery_1(dy,a,ny,nx,dya)
         elseif (nc1.eq.2) then
            call dery_1(dy,a,ny,nx,dya)
            call dery_2(ddy,a,ny,nx,ddya)
         endif 
      endif
      end
!
! =====================================================================
! =  CALCUL DES DERIVEES COLLOCATION-CHEBYSHEV : FFT                  =
! =====================================================================
!
      subroutine deryfft(tz,u,nx,ny,wsave,dya,ddya,nchoix)
      implicit double precision(a-h,o-z)
      dimension tz(0:ny-1,2),u(0:nx-1,0:ny-1) 
      dimension a(0:nx-1,0:ny-1),dya(0:nx-1,0:ny-1) 
      dimension ddya(0:nx-1,0:ny-1)
      dimension wsave(3*ny+15),xt(nx,ny-1)
!->calcul des coefficients de chebyshev  
      do j=0,ny-1
         do i=0,nx-1
            a(i,j)=u(i,j)
         enddo
      enddo
      call vcost(nx,ny,a,xt,nx,wsave)
!->
      do j=1,ny-2
         do i=0,nx-1
            a(i,j)=a(i,j)*sqrt(2.d0/real((ny-1)))
         enddo
      enddo
      do i=0,nx-1
         a(i,0)=a(i,0)*0.5d0*sqrt(2.d0/real(ny-1)) 
         a(i,ny-1)=a(i,ny-1)*0.5d0*sqrt(2.d0/real(ny-1))
      enddo 
!
!->calcul de la derivee premiere 
! 
      do i=0,nx-1
         dya(i,ny-1)=0.d0
         dya(i,ny-2)=2.d0*real(ny-1)*a(i,ny-1)
      enddo
      do j=ny-3,0,-1
         do i=0,nx-1
            dya(i,j)=dya(i,j+2)+2.d0*(real(j)+1.d0)*a(i,j+1)
         enddo
      enddo
      call vcost(nx,ny,dya,xt,nx,wsave)
      do j=0,ny-1
         do i=0,nx-1
            dya(i,j)=-dya(i,j)*0.5d0*tz(j,1)*sqrt(2.d0*real(ny-1))
         enddo
      enddo 
! 
!->calcul de la derivee seconde 
!  
      if (nchoix.eq.2) then 
         do i=0,(nx-1) 
            ddya(i,ny-1)=0.d0 
            ddya(i,ny-2)=0.D0 
            ddya(i,ny-3)=4.d0*real((ny-1)*(ny-2))*a(i,ny-1) 
            ddya(i,ny-4)=4.d0*real((ny-2)*(ny-3))*
     &a(i,ny-2) 
         enddo 
         do j=ny-5,0,-1 
            do i=0,nx-1 
               ddya(i,j)=(-real(j+1)*ddya(i,j+4) 
     &+2.d0*real(j+2)*ddya(i,j+2))/real(j+3) 
     &+4.d0*real((j+2)*(j+1))*a(i,j+2) 
            enddo 
         enddo 
         call vcost(nx,ny,ddya,xt,nx,wsave) 
         do j=0,ny-1 
            do i=0,nx-1 
               ddya(i,j)=ddya(i,j)*0.5d0*sqrt(2.d0*real(ny-1)) 
            enddo 
         enddo 
         do j=0,ny-1
            do i=0,nx-1 
               ddya(i,j)=ddya(i,j)*tz(j,1)**2-tz(j,2)*dya(i,j) 
            enddo 
         enddo 
      endif 
      end
!
! =====================================================================
! =  CALCUL DE LA DERIVEE PREMIERE COLLOCATION-CHEBYSHEV :            =
! =               MULTIPLICATION MATRICIELLE                          =
! =====================================================================
!
      subroutine dery_1(dy,a,ny,nx,dya)
      implicit double precision (a-h,o-z)
      dimension dy(ny,ny)
      dimension a(nx,ny),dya(nx,ny)
!      dya=matmul(a,dy)
      call dgemm('n','n',nx,ny,ny,1.d0,a,nx,dy,ny,0.d0,dya,nx)
      return
      end
!
! =====================================================================
! =  CALCUL DE LA DERIVEE SECONDE COLLOCATION-CHEBYSHEV :             =
! =               MULTIPLICATION MATRICIELLE                          =
! =====================================================================
!
      subroutine dery_2(ddy,a,ny,nx,ddya)
      implicit double precision (a-h,o-z)
      dimension ddy(ny,ny)
      dimension a(nx,ny),ddya(nx,ny)
!      ddya=matmul(a,ddy)
      call dgemm('n','n',nx,ny,ny,1.d0,a,nx,ddy,ny,0.d0,ddya,nx)
      return
      end
!
! =====================================================================
! = CALCUL DE LA DERIVEE PREMIERE POUR LES DIFFERENCES FINIES         =
! =====================================================================
!
      subroutine derx_1(dx,a,nj,ni,dxa)
      implicit double precision (a-h,o-z)
      dimension dx(7,8)
      dimension a(ni+2,nj+2),dxa(ni+2,nj+2)
!
      do j=1,nj+2
!
      dxa(1,j)=dx(1,1)*a(1,j)+dx(1,2)*a(2,j)+
     &dx(1,3)*a(3,j)+dx(1,4)*a(4,j)+dx(1,5)*a(5,j)
!
      dxa(2,j)=dx(2,1)*a(1,j)+dx(2,2)*a(2,j)+
     &dx(2,3)*a(3,j)+dx(2,4)*a(4,j)+dx(2,5)*a(5,j)
!
      dxa(3,j)=dx(3,1)*a(1,j)+dx(3,2)*a(2,j)+
     &dx(3,3)*a(4,j)+dx(3,4)*a(5,j)
!
      dxa(4,j)=dx(4,1)*a(1,j)+dx(4,2)*a(2,j)+
     &dx(4,3)*a(3,j)+dx(4,4)*a(5,j)+dx(4,5)*a(6,j)+
     &dx(4,6)*a(7,j)
!
      i=ni-1
      dxa(i,j)=dx(4,1)*a(i-3,j)+dx(4,2)*a(i-2,j)+
     &dx(4,3)*a(i-1,j)+dx(4,4)*a(i+1,j)+dx(4,5)*a(i+2,j)+
     &dx(4,6)*a(i+3,j)
!
      i=ni
      dxa(i,j)=dx(3,1)*a(i-2,j)+dx(3,2)*a(i-1,j)+
     &dx(3,3)*a(i+1,j)+dx(3,4)*a(i+2,j)
!
      i=ni+1
      dxa(i,j)=dx(6,1)*a(i-3,j)+dx(6,2)*a(i-2,j)+
     &dx(6,3)*a(i-1,j)+dx(6,4)*a(i,j)+dx(6,5)*a(i+1,j)
!
      i=ni+2
      dxa(i,j)=dx(7,1)*a(i-4,j)+dx(7,2)*a(i-3,j)+
     &dx(7,3)*a(i-2,j)+dx(7,4)*a(i-1,j)+dx(7,5)*a(i,j)
!
      enddo
!
      do j=1,nj+2
      do i=5,ni-2
      dxa(i,j)=dx(5,1)*a(i-4,j)+
     &dx(5,2)*a(i-3,j)+dx(5,3)*a(i-2,j)+dx(5,4)*a(i-1,j)+
     &dx(5,5)*a(i+1,j)+dx(5,6)*a(i+2,j)+dx(5,7)*a(i+3,j)+
     &dx(5,8)*a(i+4,j)
      enddo
      enddo
!
      return
      end
!
! =====================================================================
! = CALCUL DE LA DERIVEE PREMIERE D'UN ELEMENT D'UN VECTEUR i         =
! = POUR LES DIFFERENCES FINIES                                       =
! =====================================================================
!
      subroutine derdxcl(i,dx,a,ni,res)
      implicit double precision (a-h,o-z)
      dimension dx(7,8),ddx(0:4,5),a(ni)
!
      if (i.eq.1) then
         res=dx(1,1)*a(i)+dx(1,2)*a(i+1)+
     &dx(1,3)*a(i+2)+dx(1,4)*a(i+3)+dx(1,5)*a(i+4)
         return
         endif
!
      if (i.eq.2) then
         res=dx(2,1)*a(i-1)+dx(2,2)*a(i)+
     &dx(2,3)*a(i+1)+dx(2,4)*a(i+2)+dx(2,5)*a(i+3)
         return
         endif
!
      if (i.eq.3) then
         res=dx(3,1)*a(i-2)+dx(3,2)*a(i-1)
     &+dx(3,3)*a(i+1)+dx(3,4)*a(i+2)
         return
         endif
!
      if (i.eq.4) then
         res=dx(4,1)*a(i-3)+dx(4,2)*a(i-2)
     &+dx(4,3)*a(i-1)+dx(4,4)*a(i+1)+dx(4,5)*a(i+2)
     &+dx(4,6)*a(i+3)
         return
         endif
!
      if ((i.ge.5).and.(i.le.(ni-4))) then
         res=dx(5,1)*a(i-4)
     &+dx(5,2)*a(i-3)+dx(5,3)*a(i-2)+dx(5,4)*a(i-1)
     &+dx(5,5)*a(i+1)+dx(5,6)*a(i+2)+dx(5,7)*a(i+3)
     &+dx(5,8)*a(i+4)
         return
         endif
!
      if (i.eq.(ni-3)) then
         res=dx(4,1)*a(i-3)+dx(4,2)*a(i-2)
     &+dx(4,3)*a(i-1)+dx(4,4)*a(i+1)+dx(4,5)*a(i+2)
     &+dx(4,6)*a(i+3)
         return
         endif
!
      if (i.eq.(ni-2)) then
         res=dx(3,1)*a(i-2)+dx(3,2)*a(i-1)
     &+dx(3,3)*a(i+1)+dx(3,4)*a(i+2)
         return
         endif
!
      if (i.eq.(ni-1)) then
         res=dx(6,1)*a(i-3)+dx(6,2)*a(i-2)
     &+dx(6,3)*a(i-1)+dx(6,4)*a(i)+dx(6,5)*a(i+1)
         return
         endif
!
      if (i.eq.ni) then
         res=dx(7,1)*a(i-4)+dx(7,2)*a(i-3)
     &+dx(7,3)*a(i-2)+dx(7,4)*a(i-1)+dx(7,5)*a(i)
         return
         endif
!
!
      return
      end
!
! =====================================================================
! = CALCUL DE LA DERIVEE SECONDE D'UN ELEMENT D'UN VECTEUR i          =
! = POUR LES DIFFERENCES FINIES                                       =
! =====================================================================
!
      subroutine derddxcl(i,ddx,a,ni,res)
      implicit double precision (a-h,o-z)
      dimension ddx(0:4,5),a(ni)
!
      if (i.eq.1) then
         res=ddx(0,1)*a(i)+ddx(0,2)*a(i+1)+
     &ddx(0,3)*a(i+2)+ddx(0,4)*a(i+3)+ddx(0,5)*a(i+4)
         return
         endif
!
      if (i.eq.2) then
         res=ddx(1,1)*a(i-1)+ddx(1,2)*a(i)+
     &ddx(1,3)*a(i+1)+ddx(1,4)*a(i+2)+ddx(1,5)*a(i+3)
         return
         endif
!
      if ((i.ge.3).and.(i.le.(ni-2))) then
         res=ddx(2,1)*a(i-2)
     &+ddx(2,2)*a(i-1)+ddx(2,3)*a(i)+ddx(2,4)*a(i+1)
     &+ddx(2,5)*a(i+2)
         return
         endif
!
      if (i.eq.(ni-1)) then
         res=ddx(3,1)*a(i-3)+ddx(3,2)*a(i-2)
     &+ddx(3,3)*a(i-1)+ddx(3,4)*a(i)+ddx(3,5)*a(i+1)
         return
         endif
!
      if (i.eq.ni) then
         res=ddx(4,1)*a(i-4)+ddx(4,2)*a(i-3)
     &+ddx(4,3)*a(i-2)+ddx(4,4)*a(i-1)+ddx(4,5)*a(i)
         return
         endif
!
!
      return
      end
!
! =====================================================================
! = CALCUL DES CONDITIONS LIMITES POUR LA PRESSION EN X               =
! =====================================================================
! 
      subroutine neumx(nx,ny,p,dx,clxen,clxsn) 
      implicit double precision (a-h,o-z) 
      dimension p(nx,ny),dx(7,8)
      dimension clxen(ny),clxsn(ny) 
! 
      do j=2,ny-1
      p(1,j)=-(dx(1,2)*p(2,j)+dx(1,3)*p(3,j)+dx(1,4)*p(4,j)+ 
     &dx(1,5)*p(5,j))/dx(1,1)+clxen(j)/dx(1,1)
      p(nx,j)=-(dx(7,1)*p(nx-4,j)+dx(7,2)*p(nx-3,j)+ 
     &dx(7,3)*p(nx-2,j)+dx(7,4)*p(nx-1,j))
     &/dx(7,5)+clxsn(j)/dx(7,5)
      enddo
! 
      end 
!  
! =====================================================================
! = CALCUL DES CONDITIONS LIMITES POUR LA PRESSION EN Y               =
! =====================================================================
! 
      subroutine neumy(nx,ny,p,dy,clypn,clyin) 
      implicit double precision (a-h,o-z) 
      dimension p(nx,ny),dy(ny,ny)
      dimension clypn(nx),clyin(nx) 
!  
      aux1=dy(1,1)*dy(ny,ny)-dy(1,ny)*dy(ny,1)
      do i=2,nx-1
      p(i,1)=(dy(ny,ny)*clypn(i)-
     &dy(1,ny)*clyin(i))/aux1
      p(i,ny)=(dy(1,1)*clyin(i)-
     &dy(ny,1)*clypn(i))/aux1
      do l=2,ny-1
      p(i,1)=p(i,1)+((dy(1,ny)*dy(ny,l)
     &-dy(ny,ny)*dy(1,l))*
     &p(i,l))/aux1
      p(i,ny)=p(i,ny)+((dy(ny,1)*
     &dy(1,l)-dy(1,1)*dy(ny,l))*
     &p(i,l))/aux1
      enddo
      enddo 
      end 

