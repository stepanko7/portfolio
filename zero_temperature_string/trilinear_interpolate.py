#!/usr/bin/env python

"""This is fully numpy-vectorized implementation of the trilinear interpolation.
   https://en.wikipedia.org/wiki/trilinear_interpolation"""

__author__ = "Stepan Hlushak"
__copyright__ = "2019, Stepan Hlushak"
__license__ = "GPL"
__version__ = "0.0.1"
__maintainer__ = "Stepan Hlushak"
__email__ = "stepan.hlushak@gmail.com"
__status__ = "Production"

import numpy as np

# This is fully numpy-vectorized implementation of the trilinear interpolation.
# https://en.wikipedia.org/wiki/trilinear_interpolation



# xp,yp,zp are numpy 1D-arrays giving a path long which we have to interpolate
# the tabulated function f.
# f is presented as 3D-array.
# x,y and z represent grid on which f is tabulated

#reverse = False is compatible with meshgrid x,y,z order (axis0=z, axis1=y, axis2=x)
#
#reverse = True is compatible with readgnuplot4d(reverse=True) x,y,z order (axis0=x, axis1=y, axis2=z)
def trilin_int(xp,yp,zp,x,y,z,f,reverse=False):

  xpm = np.min(x)
  ypm = np.min(y)
  zpm = np.min(z)

  xpx = np.max(x)
  ypx = np.max(y)
  zpx = np.max(z)


  dx = x[1]-x[0]
  dy = y[1]-y[0]
  dz = z[1]-z[0]

# here we calculate arrays of indices of the f 2D-array
# which correspond to the lower left vertex of the cube
# in which we interpolate
  xpi = np.array(np.floor((xp-xpm)/dx),int)
  ypi = np.array(np.floor((yp-ypm)/dy),int)
  zpi = np.array(np.floor((zp-zpm)/dz),int)



  # f function values on the vertices of the cube
  if reverse:
    f000 = f[xpi,ypi,zpi]      #reverse
    f100 = f[xpi+1,ypi,zpi]
    f010 = f[xpi,ypi+1,zpi]
    f001 = f[xpi,ypi,zpi+1]
    f101 = f[xpi+1,ypi,zpi+1]
    f110 = f[xpi+1,ypi+1,zpi]
    f011 = f[xpi,ypi+1,zpi+1]
    f111 = f[xpi+1,ypi+1,zpi+1]
  else:
    f000 = f[zpi,ypi,xpi]
    f100 = f[zpi+1,ypi,xpi]
    f010 = f[zpi,ypi+1,xpi]
    f001 = f[zpi,ypi,xpi+1]
    f101 = f[zpi+1,ypi,xpi+1]
    f110 = f[zpi+1,ypi+1,xpi]
    f011 = f[zpi,ypi+1,xpi+1]
    f111 = f[zpi+1,ypi+1,xpi+1]


  xf = (xp-xpm)/dx - xpi # relative postion of the path points in the cube
  yf = (yp-ypm)/dy - ypi
  zf = (zp-zpm)/dz - zpi

# finally interpolating f along the path
# Note that all opertaions are vectorized! We don't have any loops!
  if reverse:
    fp = f000*(1-xf)*(1-yf)*(1-zf) \
          + f100*xf*(1-yf)*(1-zf)  \
          + f010*(1-xf)*yf*(1-zf)  \
          + f001*(1-xf)*(1-yf)*zf  \
          + f110*xf*yf*(1-zf)      \
          + f101*xf*(1-yf)*zf      \
          + f011*(1-xf)*yf*zf      \
          + f111*xf*yf*zf
  else:
    fp = f000*(1-zf)*(1-yf)*(1-xf) \
          + f100*zf*(1-yf)*(1-xf)  \
          + f010*(1-zf)*yf*(1-xf)  \
          + f001*(1-zf)*(1-yf)*xf  \
          + f110*zf*yf*(1-xf)      \
          + f101*zf*(1-yf)*xf      \
          + f011*(1-zf)*yf*xf      \
          + f111*zf*yf*xf

  return fp


if(__name__=='__main__'):
  x = np.array([1.,2.,3.,4.])
  x = x.reshape((4,1,1))
  y = x.reshape((1,4,1))
  z = x.reshape((1,1,4))
  f = x*y*z
  x = x.reshape((4))
  y = x.reshape((4))
  z = x.reshape((4))
  path = [(1.5,1.5,1.5), (2.5,2.5,2.5), (3.5,3.5,3.5),]
  print trilin_int(path[0],path[1],path[2],x,y,z,f)

