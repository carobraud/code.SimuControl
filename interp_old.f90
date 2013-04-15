subroutine interp(nx,ny,xmesh,ymesh,u,v,num)
  implicit none
  integer,intent(in) :: nx,ny
  real(8),intent(in) :: xmesh(nx,ny),ymesh(nx,ny) 
  real(8),intent(inout) :: u(nx,ny),v(nx,ny) 
  character(*) :: num
! interpolation controle
  integer :: nxi,nyi
  real(4) :: xia,xib,yia,yib
  real(4),allocatable :: xi(:,:),yi(:,:)
  real(4),allocatable :: ui(:,:),vi(:,:)
  integer :: i,j,ii,ji,idy
  real(4),allocatable :: xii(:,:),yii(:,:),uii(:,:),vii(:,:)
  integer,allocatable :: idx(:)
  character(100) :: pathresults,fileout
! control
  integer :: nxc,ica,icb
  real(4) :: xc(nx),uc(nx),vc(nx),xca,xcb

! parametres interpolation
  nxi=10
  nyi=10
  xia=20.
  xib=120.
  yia=0.
  yib=15.
  allocate(xi(nxi,nyi),yi(nxi,nyi))
  allocate(ui(nxi,nyi),vi(nxi,nyi))

! calcul point maillage interpolation
  do j=1,nyi
     do i=1,nxi
        xi(i,j)=xia+real(i-1)*(xib-xia)/(real(nxi)-1.)
        yi(i,j)=yia+real(j-1)*(yib-yia)/(real(nyi)-1.)
     enddo
  enddo
  
! interpolation in x direction
  allocate(xii(nxi,ny),yii(nxi,ny),uii(nxi,ny),vii(nxi,ny)) 
  allocate(idx(nxi))
  do i=1,nxi
     do j=1,nx-1
        if (xi(i,1)>=xmesh(j,1).and.xi(i,1)<=xmesh(j+1,1)) then
           idx(i)=j
           exit
        endif
     enddo
  enddo
  do j=1,ny
     do i=1,nxi
        uii(i,j)=u(idx(i),j)+(xi(i,1)-xmesh(idx(i),1))*&
             (u(idx(i)+1,j)-u(idx(i),j))*&
             (xmesh(idx(i)+1,j)-xmesh(idx(i),j))
        vii(i,j)=v(idx(i),j)+(xi(i,1)-xmesh(idx(i),1))*&
             (v(idx(i)+1,j)-v(idx(i),j))*&
             (xmesh(idx(i)+1,j)-xmesh(idx(i),j))
        yii(i,j)=ymesh(idx(i),j)+(xi(i,1)-xmesh(idx(i),1))*&
             (ymesh(idx(i)+1,j)-ymesh(idx(i),j))*&
             (xmesh(idx(i)+1,j)-xmesh(idx(i),j))
     enddo
  enddo

  ! interpolation in y direction
  do ji=1,nyi
     do ii=1,nxi
        ! Are we inside the original mesh ?
        if (yi(ii,ji)>=yii(ii,1)) then
           ! inside mesh -> interp
           do j=1,ny
              if (yi(ii,ji)>=yii(ii,j).and.yi(ii,ji)<=yii(ii,j+1)) then
                 idy=j
                 exit
              endif
           enddo
           ui(ii,ji)=uii(ii,idy)+(yi(ii,ji)-yii(ii,idy))*&
                (ui(ii,idy+1)-ui(ii,idy))*&
                (yii(ii,idy+1)-yii(ii,idy))
           vi(ii,ji)=vii(ii,idy)+(yi(ii,ji)-yii(ii,idy))*&
                (vi(ii,idy+1)-vi(ii,idy))*&
                (yii(ii,idy+1)-yii(ii,idy))
        else
           ! not inside mesh -> zero
           ui(ii,ji)=0.
           vi(ii,ji)=0.
        endif
     enddo
  enddo
  pathresults='./output/'        
  fileout='Uinst'//num//'.nc'

  !-> write netcdf 
  call Write_Ncdf(pathresults,fileout,nxi,nyi,xi,yi(1,:),ui,vi)

  !-> controle
  call controle(nxi,nyi,xi,yi(1,:),ui,vi,nxc,xc,uc,vc,xca,xcb)
  print*,xca,xcb,xc(1:nxc),uc(1:nxc),vc(1:nxc)

  !-> interpolation wall control 
  !-> search ica and icb
  do i=1,nx
     if (xca>=xmesh(i,1).and.xca<=xmesh(i+1,1)) ica=i
     if (xcb>=xmesh(i,1).and.xcb<=xmesh(i+1,1)) icb=i
  enddo

  do i=ica,icb
     do ii=1,nxc
        if (xmesh(i,1)>=xc(ii).and.xmesh(i,1)<=xc(ii+1)) then
           idx(1)=ii
        endif
     enddo
     u(i,1)=uc(idx(1))+(xmesh(i,1)-xc(idx(1)))*&
          (uc(idx(1)+1)-uc(idx(1)))*&
          (xc(idx(1)+1)-xc(idx(1)))
     v(i,1)=vc(idx(1))+(xmesh(i,1)-xc(idx(1)))*&
          (vc(idx(1)+1)-vc(idx(1)))*&
          (xc(idx(1)+1)-xc(idx(1)))
     print*,xmesh(i,1),xc(idx(1)),xc(idx(1)+1),vc(idx(1)),v(i,1),vc(idx(1)+1)
  enddo

  deallocate(xi,yi,ui,vi)
  deallocate(xii,yii,uii,vii) 
  deallocate(idx)


end subroutine interp
