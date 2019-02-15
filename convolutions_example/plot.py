import sys
import numpy as np
#import scipy as sp
import matplotlib.pyplot as plt
# import glob

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cmx
import matplotlib.colors as clrsm
from pylab import *


def read_conv2d(file_in):
  if type(file_in) == type(sys.stdin):
    readfile = file_in
  else:
    readfile = file(file_in,'r')
  lines = readfile.readlines()

  words = [s for s in lines[0].split() if len(s)>0]
  nx,ny = (int(words[0]),int(words[1]))

  nlns = nx*ny

  data = []
  for l in lines[1:]:
    words = [float(s) for s in l.split() if len(s)>0]
    data.append(words)
  data = np.array(data)
  rho_2d = np.zeros((nx,ny))
  w3_2d = np.zeros((nx,ny))
  conv_2d = np.zeros((nx,ny))

  for i in range(nx):
    for j in range(ny):
      ig = i + j*nx
      rho_2d[i,j] = data[ig,0]
      w3_2d[i,j] = data[ig,1]
      conv_2d[i,j] = data[ig,2]

  return nx, ny, rho_2d, w3_2d, conv_2d


if(__name__=='__main__'):

    files = ["convolution_1d.dat"]
    for f in files:
      t = np.loadtxt(f)
      n,m = t.shape
      indx = t[:,0]
      rho = t[:,1]
      w3 = t[:,2]
      conv = t[:,3]
      plt.plot(indx, rho, "k", label="sqaure 1")
      plt.plot(indx, w3, "g", label="square 2" )
      plt.plot(indx, conv, "r", label="triangle")

    lgnd = plt.legend(loc=0, labelspacing=-0.001)
    for t in lgnd.texts:
      t.set_fontsize(10)

    v = plt.axis()
    xmin, xmax, ymin, ymax = v[0],v[1],v[2],v[3]
    plt.axis([xmin, xmax, ymin, ymax])
    plt.savefig('convolution_1.png')
    plt.close()




    nx, ny, rho_2d, w3_2d, conv_2d = read_conv2d("convolution_2d.dat")

    x = np.array(range(100))*0.1
    y = np.array(range(100))*0.1
    xx,yy = np.meshgrid(x,y)

    xsize = 5
    ysize = xsize*0.7/0.9
    fig = plt.figure(1,(xsize, ysize),dpi=300)
    plt.subplots_adjust(bottom=0.130,left=0.120,right=0.97,top=0.915, hspace=0.001)
    pcolor(xx, yy, rho_2d, vmax=1)
    colorbar()
    savefig('rho_2d.png', dpi=300)
    close()

    xsize = 5
    ysize = xsize*0.7/0.9
    fig = plt.figure(1,(xsize, ysize),dpi=300)
    plt.subplots_adjust(bottom=0.130,left=0.120,right=0.97,top=0.915, hspace=0.001)
    pcolor(xx, yy, w3_2d, vmax=1)
    colorbar()
    savefig('w3_2d.png', dpi=300)
    close()

    xsize = 5
    ysize = xsize*0.7/0.9
    fig = plt.figure(1,(xsize, ysize),dpi=300)
    plt.subplots_adjust(bottom=0.130,left=0.120,right=0.97,top=0.915, hspace=0.001)
    pcolor(xx, yy, conv_2d, vmax=1)
    colorbar()
    savefig('conv_2d.png', dpi=300)
    close()

