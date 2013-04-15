#!/usr/bin/python
import matplotlib
#matplotlib.use("WXAgg")
from pylab import*
import numpy
import scipy.interpolate
import matplotlib.colors
import sys
from scipy.io.numpyio import fwrite, fread

filename=sys.argv[1]
#slide_number=int(sys.argv[2])
#print dir,slide_number

#rcParams['font.colors']='w'


fd=open(filename, 'rb')
datatype = 'i'
size = 2
shape=(2)
size=fread(fd, size, datatype)
nx=size[0]
ny=size[1]
data = zeros((nx,ny,3))

datatype='d'
size=nx*ny*3
shape=(3,ny,nx)
data=fread(fd, size, datatype)
data=data.reshape(shape)
    
print "Dimensions",nx,ny

nx1=0
nx2=800
ny1=0
ny2=60
 
# read the grids
grid_x = data[0,ny1:ny2,nx1:nx2]
grid_y = data[1,ny1:ny2,nx1:nx2]

# read the data 
u = data[2,ny1:ny2,nx1:nx2]
    
n1=nx2-nx1
n2=ny2-ny1
x1=grid_x
x2=grid_y
f=u[:,:]

print f.shape

n1i=16*n1
n2i=16*n2
x1i=numpy.linspace(x1.min(),x1.max(),n1i)
x2i=numpy.linspace(x2.min(),x2.max(),n2i)
x2i_i=numpy.array(zeros((n2,n1i),numpy.float32))
u2di=numpy.array(zeros((n2,n1i),numpy.float32))
u2df=numpy.array(zeros((n2i,n1i),numpy.float32))

for j in range(n2):
    uint=scipy.interpolate.interp1d(x1[j,:],f[j,:])
    u2di[j,:]=uint(x1i[:])
    yint=scipy.interpolate.interp1d(x1[j,:],x2[j,:])
    x2i_i[j,:]=yint(x1i[:])

for i in range(n1i):
    uint=scipy.interpolate.interp1d(x2i_i[:,i],u2di[:,i],kind='linear',bounds_error=False)
    u2df[:,i]=uint(x2i[:])

ratio=(x1i.max()-x1i.min())/(x2i.max()-x2i.min())*0.9

# choose size of plot
nw=12
nh=nw/ratio
fig=figure(1,(nw,nh),edgecolor='w')

ax = axes((0.05,0.,1.04,1.))
#ax.set_aspect(1) 

# Choose palette and empty color
#pal= cm.gist_stern
#pal= cm.hot
#pal= cm.jet
#pal= cm.Greens
#pal= cm.Blues_r
pal= cm.RdBu_r

# change number of color in palette
cdict=pal._segmentdata.copy()
for key in ('red','green','blue'):
    D = array(cdict[key])
#    print D
my_cmap = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,1024)

# mask empty values (for not plotting them)
u2df=numpy.ma.array(u2df,mask=numpy.isnan(u2df))

#hold(True)
my_cmap.set_bad('black',1.)

ax.set_xlabel( r'$x$',fontsize=20)
ax.set_ylabel( r'$y$',fontsize=20)
nc=(linspace(0,u2df.max(),u2df.max()+1))
print nc

# Creating image
#contourf(u2df,nc,cmap=pal,linewidths=1.,extent=[x1i.min(),x1i.max(),x2i.min(),x2i.max()])
#pcolor(u2df,shading='flat') 
imshow(u2df,cmap=my_cmap, interpolation='bicubic',origin='lower', extent=[x1i.min(),x1i.max(),x2i.min(),x2i.max()],vmin=-1,vmax=1.)
#imshow(u2df,cmap=my_cmap,origin='lower', extent=[x1i.min(),x1i.max(),x2i.min(),x2i.max()],vmin=u2df.min(),vmax=0)


colorbar(shrink=0.6)
show()
ext = ".png"
fig.savefig('image'+ext,dpi=100,facecolor='black')

clf()
#plot(x1i[:],u2df[n2i/3-1,:], 'kx-')
#plot(x1i[:],u2df[n2i/3,:], 'k+-')
#plot(x1i[:],u2df[n2i/3+1,:], 'ko-')
#print u2df[n2i/3-1,:]
#print u2df[n2i/3,:]
#print u2df[n2i/3+1,:]
#show()

