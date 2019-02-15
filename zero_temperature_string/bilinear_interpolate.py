#!/usr/bin/env python

"""This is fully numpy-vectorized implementation of the bilinear interpolation.
   https://en.wikipedia.org/wiki/Bilinear_interpolation"""

__author__ = "Stepan Hlushak"
__copyright__ = "2019, Stepan Hlushak"
__license__ = "GPL"
__version__ = "0.0.1"
__maintainer__ = "Stepan Hlushak"
__email__ = "stepan.hlushak@gmail.com"
__status__ = "Production"


import numpy as np


# This is fully numpy-vectorized implementation of the bilinear interpolation.
# https://en.wikipedia.org/wiki/Bilinear_interpolation

# xp and yp are numpy 1D-arrays giving a path long which we have to interpolate
# the tabulated function f
# f is presented as 2D-array.
# x and y represent grid on which f is tabulated

def bilin_int(xp,yp,x,y,f,reverse=False):
#reverse = False is compatible with meshgrid x,y,z order
#reverse = True is compatible with readgnuplot4d x,y,z order

  xpm = np.min(x)
  ypm = np.min(y)

  xpx = np.max(x)
  ypx = np.max(y)


  dx = x[1]-x[0]
  dy = y[1]-y[0]


# here we calculate arrays of indices of the f 2D-array
# which correspond to the lower left vertex of the square
# in which we interpolate
  xpi = np.array(np.floor((xp-xpm)/dx),int)
  ypi = np.array(np.floor((yp-ypm)/dy),int)


  # f values on the vertices of the square
  if reverse:
    f00 = f[xpi,ypi]
    f10 = f[xpi+1,ypi]
    f01 = f[xpi,ypi+1]
    f11 = f[xpi+1,ypi+1]
  else:
    f00 = f[ypi,xpi]
    f10 = f[ypi+1,xpi]
    f01 = f[ypi,xpi+1]
    f11 = f[ypi+1,xpi+1]

  xf = (xp-xpm)/dx - xpi # relative postion of the path points
  yf = (yp-ypm)/dy - ypi # in the squre

# finally interpolating f along the path
# Note that all opertaions are vectorized! We don't have any loops!
  if reverse:
    fp = f00*(1-xf)*(1-yf) \
       + f10*xf*(1-yf)  \
       + f01*(1-xf)*yf  \
       + f11*xf*yf
  else:
    fp = f00*(1-yf)*(1-xf) \
       + f10*yf*(1-xf)  \
       + f01*(1-yf)*xf  \
       + f11*yf*xf

  return fp



if(__name__=='__main__'):
  x = np.array([1.,2.,3.,4.])
  y = np.array([0.,.25,.5,.75,1.])
  x = x.reshape((1,4))
  y = y.reshape((5,1))
  f = x*y
  x = x.reshape((4))
  y = y.reshape((5))
  path = [(2,2.5,3.5), (0.1,0.25,0.5)]
  print bilin_int(path[0],path[1],x,y,f)

