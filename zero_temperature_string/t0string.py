#Zero-temperature string method (ZTS) code 
#applied to the Mueller potential.

import os
import numpy as np
import scipy as sp
from scipy.interpolate import interp1d
import scipy.integrate as spint
import bilinear_interpolate as bi
import trilinear_interpolate as ti
# import matplotlib.pyplot as plt






def get_t0string_3Dpath(x,y,z, xg,yg,zg, dVdx,dVdy,dVdz, nstepmax=20000, tol1=1e-7, h=1e-4, reverse=False, endfix=True):

  # dx,dy,dz -- x,y,z steps of the grid
  # dVdx  -- partial derivative of V on x
  # dVdy  -- partial derivative of V on y
  # dVdz  -- partial derivative of V on z
  # xg,yg,zg -- x,y,z 1D coordinates of the grid from whch 2D coordintes are formed with np.meshgrid
  # xx,yy,zz -- x,y,z 2D coordinates of the grid
  # V     -- FEP given on a grid

  n1 = len(x)
  g1 = np.linspace(0,1,n1);


  dx = x-np.roll(x,1)
  dy = y-np.roll(y,1)
  dz = z-np.roll(z,1)
  dx[0] = 0
  dy[0] = 0
  dz[0] = 0
  
  dl = np.sqrt(dx**2+dy**2+dz**2)
  lxyz = np.cumsum(dl)
  lxyz = lxyz/lxyz[-1]


  # intepolate
  x = np.interp(g1,lxyz,x)  
  y = np.interp(g1,lxyz,y)
  z = np.interp(g1,lxyz,z)


  # tangent vectors
  dx = x-np.roll(x,1)
  dy = y-np.roll(y,1)
  dz = z-np.roll(z,1)
  dx[0] = 0
  dy[0] = 0
  dz[0] = 0
  dx = np.roll(dx,-1) #now zero is dx[-1]
  dy = np.roll(dy,-1) #now zero is dy[-1]
  dz = np.roll(dz,-1) #now zero is dz[-1]
  
  dl = np.sqrt(dx**2+dy**2+dz**2)
  dl[-1] = 1.0 #one instead of 0


  vx = dx/dl
  vy = dy/dl
  vz = dz/dl
  vx[-1], vy[-1], vz[-1] = vx[-2], vy[-2], vz[-2]
  vx2 = np.roll(vx,1)
  vy2 = np.roll(vy,1)
  vz2 = np.roll(vz,1)
  vx2[0], vy2[0], vz2[0] = vx2[1], vy2[1], vz2[1]
  vx = (vx + vx2)/2
  vy = (vy + vy2)/2
  vz = (vz + vz2)/2
  vl = np.sqrt(vx**2 + vy**2 + vz**2)
  vx = vx/vl
  vy = vy/vl
  vz = vz/vl



  xi = x
  yi = y
  zi = z

  dVx = 0 #force components
  dVy = 0
  dVz = 0

  # Main loop
  for nstep in range(1,nstepmax+1):
#   calculation of the x,y and z-components of the force, dVx, dVy, and dVz respectively
      dVx = ti.trilin_int(x,y,z,xg,yg,zg,dVdx,reverse)
      dVy = ti.trilin_int(x,y,z,xg,yg,zg,dVdy,reverse)
      dVz = ti.trilin_int(x,y,z,xg,yg,zg,dVdz,reverse)
#       print "dVdx.shape: ",dVdx.shape
#       print "dVx.shape: ",dVx.shape
     
      x0 = x;
      y0 = y;
      z0 = z;

  # fix endpoints
      if endfix:
        dVx[0],dVy[0],dVz[0] = 0.0,0.0,0.0
        dVx[-1],dVy[-1],dVz[-1] = 0.0,0.0,0.0
      else:
        vx[0], vy[0], vz[0] = 0.0, 0.0, 0.0
        vx[-1], vy[-1], vz[-1] = 0.0, 0.0, 0.0

  # string steps:

  # 1. evolve
      dVv = vx*dVx + vy*dVy + vz*dVz #projection at tangent vector (vx,vy,vz)

      x = x - h*(dVx - vx*dVv);
      y = y - h*(dVy - vy*dVv);
      z = z - h*(dVz - vz*dVv);


  # 2. reparametrize  
      dx = x-np.roll(x,1)
      dy = y-np.roll(y,1)
      dz = z-np.roll(z,1)
      dx[0] = 0
      dy[0] = 0
      dz[0] = 0

      dl = np.sqrt(dx**2+dy**2+dz**2)
      lxyz = np.cumsum(dl)
      lxyz = lxyz/lxyz[-1]

      # reintepolate
#       print "g1.shape: ",g1.shape
#       print "lxyz.shape: ",lxyz.shape
#       print "x.shape: ",x.shape
      x = np.interp(g1,lxyz,x)  
      y = np.interp(g1,lxyz,y)
      z = np.interp(g1,lxyz,z)

      # new tangent vectors
      dx = x-np.roll(x,1)
      dy = y-np.roll(y,1)
      dz = z-np.roll(z,1)
      dx[0] = 0
      dy[0] = 0
      dz[0] = 0
      dx = np.roll(dx,-1) #now zero is dx[-1]
      dy = np.roll(dy,-1) #now zero is dy[-1]
      dz = np.roll(dz,-1) #now zero is dz[-1]
      
      dl = np.sqrt(dx**2+dy**2+dz**2)
      dl[-1] = 1.0 #one instead of 0
    
      vx = dx/dl
      vy = dy/dl
      vz = dz/dl
      vx[-1], vy[-1], vz[-1] = vx[-2], vy[-2], vz[-2]
      vx2 = np.roll(vx,1)
      vy2 = np.roll(vy,1)
      vz2 = np.roll(vz,1)
      vx2[0], vy2[0], vz2[0] = vx2[1], vy2[1], vz2[1]
      vx = (vx + vx2)/2
      vy = (vy + vy2)/2
      vz = (vz + vz2)/2
      vl = np.sqrt(vx**2 + vy**2 + vz**2)
      vx = vx/vl
      vy = vy/vl
      vz = vz/vl

#      print 'x=',   x
#      print 'y=',   y
#      print 'z=',   z
#      print 'dl=', dl

      tol = (np.linalg.norm(x-x0) + np.linalg.norm(y-y0) + np.linalg.norm(z-z0))/n1;
      print 'nstep, tol = ', nstep, tol  
      if tol <= tol1:
          break

  return x, y, z, tol, nstep, dVx, dVy, dVz









def get_t0string_2Dpath(x,y, xg,yg, dVdx,dVdy, nstepmax=2000, tol1=1e-7, h=1e-4, endfix=True):

  # dx,dy -- x and y steps of the grid
  # dVdx  -- partial derivative of V on x
  # dVdy  -- partial derivative of V on y
  # xg,yg -- x and y 1D coordinates of the grid from whch 2D coordintes are formed with np.meshgrid
  # xx,yy -- x and y 2D coordinates of the grid
  # V     -- FEP given on a grid

  n1 = len(x)
  g1 = np.linspace(0,1,n1);

  dx = x-np.roll(x,1)
  dy = y-np.roll(y,1)
  dx[0] = 0
  dy[0] = 0
  
  
  dl = np.sqrt(dx**2+dy**2)
  lxy = np.cumsum(dl)
  lxy = lxy/lxy[-1]

  #interpolate
  x = np.interp(g1,lxy,x)  
  y = np.interp(g1,lxy,y)


  # tangent vectors
  dx = x-np.roll(x,1)
  dy = y-np.roll(y,1)
  dx[0] = 0
  dy[0] = 0
  dx = np.roll(dx,-1) #now zero is dx[-1]
  dy = np.roll(dy,-1) #now zero is dy[-1]
  
  dl = np.sqrt(dx**2+dy**2)
  dl[-1] = 1.0 #one instead of 0


  vx = dx/dl
  vy = dy/dl
  vx[-1], vy[-1] = vx[-2], vy[-2]
  vx2 = np.roll(vx,1)
  vy2 = np.roll(vy,1)
  vx2[0], vy2[0] = vx2[1], vy2[1]
  vx = (vx + vx2)/2
  vy = (vy + vy2)/2
  vl = np.sqrt(vx**2 + vy**2)
  vx = vx/vl
  vy = vy/vl



  xi = x
  yi = y

  dVx = 0 #force components
  dVy = 0

  # Main loop
  for nstep in range(1,nstepmax+1):
#   calculation of the x and y-components of the force, dVx and dVy respectively
      dVx = bi.bilin_int(x,y,xg,yg,dVdx)
      dVy = bi.bilin_int(x,y,xg,yg,dVdy)
     
      x0 = x;
      y0 = y;


  # fix endpoints
      if endfix:
        dVx[0], dVy[0] = 0.0, 0.0
        dVx[-1], dVy[-1] = 0.0, 0.0
      else:
        vx[0], vy[0] = 0.0, 0.0
        vx[-1], vy[-1] = 0.0, 0.0

  # string steps:

  # 1. evolve
#      x = x - h*dVx;
#      y = y - h*dVy;

      dVv = vx*dVx + vy*dVy #projection at tangent vector (vx,vy,vz)

      x = x - h*(dVx - vx*dVv);
      y = y - h*(dVy - vy*dVv);

  # 2. reparametrize  
      dx = x-np.roll(x,1) 
      dy = y-np.roll(y,1) 
      dx[0] = 0
      dy[0] = 0

      dl = np.sqrt(dx**2+dy**2)
      lxy = np.cumsum(dl)
      lxy = lxy/lxy[-1]

      # reintepolate
      x = np.interp(g1,lxy,x)  
      y = np.interp(g1,lxy,y)

      # new tangent vectors
      dx = x-np.roll(x,1)
      dy = y-np.roll(y,1)
      dx[0] = 0
      dy[0] = 0
      dx = np.roll(dx,-1) #now zero is dx[-1]
      dy = np.roll(dy,-1) #now zero is dy[-1]
      
      dl = np.sqrt(dx**2+dy**2)
      dl[-1] = 1.0 #one instead of 0
    
      vx = dx/dl
      vy = dy/dl
      vx[-1], vy[-1] = vx[-2], vy[-2]
      vx2 = np.roll(vx,1)
      vy2 = np.roll(vy,1)
      vx2[0], vy2[0] = vx2[1], vy2[1]
      vx = (vx + vx2)/2
      vy = (vy + vy2)/2
      vl = np.sqrt(vx**2 + vy**2)
      vx = vx/vl
      vy = vy/vl


      tol = (np.linalg.norm(x-x0) + np.linalg.norm(y-y0))/n1;
      print 'nstep, tol = ', nstep, tol  
      
      if tol <= tol1:
          break

  return x, y, tol, nstep, dVx, dVy









if(__name__=="__main__"):
  import matplotlib.pyplot as plt
  from mpl_toolkits.mplot3d import Axes3D
  import matplotlib.cm as cmx
  import matplotlib.colors as clrsm
  from pylab import *



  if(False):  #test 3D case

    # max number of iterations
    nstepmax = 2000;
    # frequency of plotting
    nstepplot = 1e1;

    # plot string every nstepplot if flag1 = 1 
    flag1 = 1;

    # parameter used as stopping criterion
    tol1 = 1e-6;


    # number of images along the string (try from  n1 = 3 up to n1 = 1e4)
    n1 = 50;

    # time-step (limited by the ODE step on line 83 & 84 but independent of n1)
    h = 1e-4;

    # end points of the initial string
    # notice that they do NOT have to be at minima of V -- the method finds
    # those automatically
    xa = -1;
    ya = 1;
    za = 0.4;

    xb = 0.7;
    yb = 0.5;
    zb = 0.75;

    # initialization
    g1 = np.linspace(0,1,n1);
    x = (xb-xa)*g1+xa;
    y = (x-xa)*(yb-ya)/(xb-xa)+ya;
    z = (x-xa)*(zb-za)/(xb-xa)+ya;


    dx = x-np.roll(x,1) 
    dy = y-np.roll(y,1) 
    dz = z-np.roll(z,1) 
    dx[0] = 0
    dy[0] = 0
    dz[0] = 0



    lxyz = np.cumsum(np.sqrt(dx**2+dy**2+dz**2))
    lxyz = lxyz/lxyz[-1]

    #x = interp1(lxy,x,g1)   #vq = interp1(x,v,xq) -- from matlab man pages
    #y = interp1(lxy,y,g1)
    x = np.interp(g1,lxyz,x)  
    y = np.interp(g1,lxyz,y)
    z = np.interp(g1,lxyz,z)

    xi = x
    yi = y
    zi = z

    # parameters in Mueller potential
    aa = [-1, -1, -6.5, 0.7]
    bb = [0, 0, 11, 0.6]
    cc = [-10, -10, -6.5, 0.7]
    dd = [-100, -10, -6.5, 0.7]
    AA = [-200, -100, -170, 15]

    XX = [1, 0, -0.5, -1]
    YY = [0, 0.5, 1.5, 1]
    ZZ = [0.5, 0.5, 0.5, 0.5]



    #xx,yy = meshgrid(-1.5:0.01:1.2,-0.2:0.01:2)
    xg = np.arange(-1.5, 1.2, 0.01)
    yg = np.arange(-0.2,   2, 0.01)
    zg = np.arange(0,   1, 0.01)
#    xx,yy,zz = np.meshgrid(xg,yg,zg) #todo: fix here
    zz,yy,xx = np.mgrid[0:1:0.01,-0.2:2:0.01,-1.5:1.2:0.01] 

    V1 = AA[1]*np.exp( aa[0]*(xx-XX[0])**2 + bb[0]*(xx-XX[0])*(yy-YY[0]) +  cc[0]*(yy-YY[0])**2 + dd[0]*(zz-ZZ[0])**2 );
    for j in range(1,4):
        V1 = V1 + AA[j]*np.exp( aa[j]*(xx-XX[j])**2 + bb[j]*(xx-XX[j])*(yy-YY[j]) + cc[j]*(yy-YY[j])**2 );




    ee = AA[0]*np.exp(aa[0]*(xx-XX[0])**2 + bb[0]*(xx-XX[0])*(yy-YY[0]) + cc[0]*(yy-YY[0])**2 + dd[0]*(zz-ZZ[0])**2)
    dVdx = (2*aa[0]*(xx-XX[0]) + bb[0]*(yy-YY[0]))*ee
    dVdy = (bb[0]*(xx-XX[0]) + 2*cc[0]*(yy-YY[0]))*ee
    dVdz = (2*dd[0]*(zz-ZZ[0]))*ee
    for j in range(1,4):
        ee = AA[j]*np.exp(aa[j]*(xx-XX[j])**2+bb[j]*(xx-XX[j])*(yy-YY[j])+cc[j]*(yy-YY[j])**2)
        dVdx = dVdx + (2*aa[j]*(xx-XX[j])+bb[j]*(yy-YY[j]))*ee
        dVdy = dVdy + (bb[j]*(xx-XX[j])+2*cc[j]*(yy-YY[j]))*ee




    # initialization
    g1 = np.linspace(0,1,n1);
    x = (xb-xa)*g1+xa;
    y = (yb-ya)*g1+ya;
    z = (zb-za)*g1+za;

    x,y,z,tol,nstep,dVx,dVy,dVz = get_t0string_3Dpath(x,y,z, xg,yg,zg, dVdx,dVdy,dVdz, 
                                                      nstepmax=20000, tol1=1e-6, h=1e-4, 
                                                      reverse=False, endfix=True)



    # Output
    print '\n' 
    print '\n' 
    print 'ZTS calculation with %d images\n'%n1 
    if tol > tol1:
        print 'The calculation failed to converge after %d iterations\n'%nstep 
    else:
        print 'The calculation terminated after %d iterations\n'%nstep 

    
    print 'x=',x
    print 'y=',y
    print 'z=',z

    xsize = 5
    ysize = xsize*0.7/0.9
    fig = plt.figure(1,(xsize, ysize),dpi=300)
    plt.subplots_adjust(bottom=0.130,left=0.120,right=0.97,top=0.915, hspace=0.001)
    pcolor(xx[50,:,:], yy[50,:,:], V1[50,:,:], vmax=200)
    plot(x,y,'w-')
    plot(x,y,'wo', ms=2)
    colorbar()
    savefig('final3D.png', dpi=300)
    close()




    # Energy along MEP
    tx = np.roll(x,-1) - np.roll(x,1)
    ty = np.roll(y,-1) - np.roll(y,1)
    tz = np.roll(z,-1) - np.roll(z,1)

    # potential computed as integral of projection of gradV on string tangent
    Vz = spint.cumtrapz(tx*dVx + ty*dVy + tz*dVz)
    Vz = 0.5*Vz;
    Vz = Vz-np.min(Vz);
    Vz = np.hstack((0,Vz))

    print len(tx), len(dVx), len(Vz)


    ntxy = np.sqrt(tx*tx+ty*ty);
    tx = tx/ntxy;
    ty = ty/ntxy;


    # err is an estimate of the error between disrectized MEP and actual one 
    # err scale as 1/n1
    err = np.trapz(1-(tx*dVx+ty*dVy)**2/(dVx*dVx+dVy*dVy))/(n1-1);
    print 'Estimate of difference between discretized MEP and actual one: %f\n'%err

    #
    # Reinterpolate string with lots of points to get accurate energy along it
    g0 = np.linspace(0,1,1000)
    x0 = np.interp(g0,lxyz,x)
    y0 = np.interp(g0,lxyz,y)
    z0 = np.interp(g0,lxyz,z)
    ee = AA[0]*np.exp( aa[0]*(x0-XX[0])**2 + bb[0]*(x0-XX[0])*(y0-YY[0]) + cc[0]*(y0-YY[0])**2 + dd[0]*(z0-ZZ[0])**2 )
    V0=ee
    dVx = (2*aa[0]*(x0-XX[0]) + bb[0]*(y0-YY[0]))*ee
    dVy = (bb[0]*(x0-XX[0]) + 2*cc[0]*(y0-YY[0]))*ee
    dVz = (2*dd[0]*(z0-ZZ[0]))*ee
    for j in range(1,4):
        ee = AA[j]*np.exp(aa[j]*(x0-XX[j])**2 + bb[j]*(x0-XX[j])*(y0-YY[j]) + cc[j]*(y0-YY[j])**2)
        V0 = V0 + ee
        dVx = dVx + (2*aa[j]*(x0-XX[j]) + bb[j]*(y0-YY[j])  )*ee
        dVy = dVy + (  bb[j]*(x0-XX[j]) + 2*cc[j]*(y0-YY[j]))*ee

    V0=V0-np.min(V0);


#    xsize = 5
#    ysize = xsize*0.7/0.9
#    fig = plt.figure(1,(xsize, ysize),dpi=300)
#    plt.subplots_adjust(bottom=0.130,left=0.120,right=0.97,top=0.915, hspace=0.001)
#    plot(g1,Vz,'w-')
#    plot(g0,V0,'r-')
#    savefig('fep1d.png', dpi=300)
#    close()





  if(True):  #test 2D case

    # max number of iterations
    nstepmax = 2000;
    # frequency of plotting
    nstepplot = 1e1;

    # plot string every nstepplot if flag1 = 1 
    flag1 = 1;

    # parameter used as stopping criterion
    tol1 = 1e-7;


    # number of images along the string (try from  n1 = 3 up to n1 = 1e4)
    n1 = 50;

    # time-step (limited by the ODE step on line 83 & 84 but independent of n1)
    h = 1e-4;

    # end points of the initial string
    # notice that they do NOT have to be at minima of V -- the method finds
    # those automatically
    xa = -1;
    ya = 0.5;

    xb = 0.7;
    yb = 0.5;

    # initialization
    g1 = np.linspace(0,1,n1);
    x = (xb-xa)*g1+xa;
    y = (x-xa)*(yb-ya)/(xb-xa)+ya;


    #dx = x-circshift(x,[0 1]); #replace this with python equivalent
    #dy = y-circshift(y,[0 1]); 
    #dx(1) = 0;
    #dy(1) = 0;
    #
    dx = x-np.roll(x,1) 
    dy = y-np.roll(y,1) 
    dx[0] = 0
    dy[0] = 0



    lxy = np.cumsum(np.sqrt(dx**2+dy**2))
    lxy = lxy/lxy[-1]

    #x = interp1(lxy,x,g1)   #vq = interp1(x,v,xq) -- from matlab man pages
    #y = interp1(lxy,y,g1)
    x = np.interp(g1,lxy,x)  
    y = np.interp(g1,lxy,y)

    xi = x
    yi = y

    # parameters in Mueller potential
    aa = [-1, -1, -6.5, 0.7]
    bb = [0, 0, 11, 0.6]
    cc = [-10, -10, -6.5, 0.7]
    AA = [-200, -100, -170, 15]

    XX = [1, 0, -0.5, -1]
    YY = [0, 0.5, 1.5, 1]



    #xx,yy = meshgrid(-1.5:0.01:1.2,-0.2:0.01:2)
    xg = np.arange(-1.5, 1.2, 0.01)
    yg = np.arange(-0.2,   2, 0.01)
    xx,yy = np.meshgrid(xg,yg)

    V1 = AA[1]*np.exp(aa[0]*(xx-XX[0])**2+bb[0]*(xx-XX[0])*(yy-YY[0])+cc[0]*(yy-YY[0])**2);
    for j in range(1,4):
        V1 = V1 + AA[j]*np.exp(aa[j]*(xx-XX[j])**2+bb[j]*(xx-XX[j])*(yy-YY[j])+cc[j]*(yy-YY[j])**2);




    ee = AA[0]*np.exp(aa[0]*(xx-XX[0])**2 + bb[0]*(xx-XX[0])*(yy-YY[0]) + cc[0]*(yy-YY[0])**2)
    dVdx = (2*aa[0]*(xx-XX[0]) + bb[0]*(yy-YY[0]))*ee
    dVdy = (bb[0]*(xx-XX[0]) + 2*cc[0]*(yy-YY[0]))*ee
    for j in range(1,4):
        ee = AA[j]*np.exp(aa[j]*(xx-XX[j])**2+bb[j]*(xx-XX[j])*(yy-YY[j])+cc[j]*(yy-YY[j])**2)
        dVdx = dVdx + (2*aa[j]*(xx-XX[j])+bb[j]*(yy-YY[j]))*ee
        dVdy = dVdy + (bb[j]*(xx-XX[j])+2*cc[j]*(yy-YY[j]))*ee



    # initialization
    g1 = np.linspace(0,1,n1);
    x = (xb-xa)*g1+xa;
    y = (x-xa)*(yb-ya)/(xb-xa)+ya;

    x,y,tol,nstep,dVx,dVy = get_t0string_2Dpath(x,y, xg,yg, dVdx,dVdy, nstepmax=20000, h=1e-5, tol1=1e-30, endfix=False)



    # Output
    print '\n' 
    print '\n' 
    print 'ZTS calculation with %d images\n'%n1 
    if tol > tol1:
        print 'The calculation failed to converge after %d iterations\n'%nstep 
    else:
        print 'The calculation terminated after %d iterations\n'%nstep 



    xsize = 5
    ysize = xsize*0.7/0.9
    fig = plt.figure(1,(xsize, ysize),dpi=300)
    plt.subplots_adjust(bottom=0.130,left=0.120,right=0.97,top=0.915, hspace=0.001)
    pcolor(xx, yy, V1, vmax=50)
    plot(x,y,'w-')
    plot(x,y,'wo', ms=2)
    colorbar()
    savefig('final.png', dpi=300)
    close()




    # Energy along MEP
    #tx = circshift(x,[0 -1])-circshift(x,[0 1]);
    #ty = circshift(y,[0 -1])-circshift(y,[0 1]);
    #
    tx = np.roll(x,-1) - np.roll(x,1)
    ty = np.roll(y,-1) - np.roll(y,1)

    # potential computed as integral of projection of gradV on string tangent
    #Vz = cumtrapz(tx*dVx + ty*dVy);
    Vz = spint.cumtrapz(tx*dVx + ty*dVy)
    Vz = 0.5*Vz;
    Vz = Vz-np.min(Vz);
    Vz = np.hstack((0,Vz))

    print len(tx), len(dVx), len(Vz)


    ntxy = np.sqrt(tx*tx+ty*ty);
    tx = tx/ntxy;
    ty = ty/ntxy;


    # err is an estimate of the error between disrectized MEP and actual one 
    # err scale as 1/n1
    err = np.trapz(1-(tx*dVx+ty*dVy)**2/(dVx*dVx+dVy*dVy))/(n1-1);
    print 'Estimate of difference between discretized MEP and actual one: %f\n'%err

    #
    # Reinterpolate string with lots of points to get accurate energy along it
    g0 = np.linspace(0,1,1000)
    x0 = np.interp(g0,lxy,x)
    y0 = np.interp(g0,lxy,y)
    ee = AA[0]*np.exp(aa[0]*(x0-XX[0])**2+bb[0]*(x0-XX[0])*(y0-YY[0])+cc[0]*(y0-YY[0])**2)
    V0=ee
    dVx = (2*aa[0]*(x0-XX[0])+bb[0]*(y0-YY[0]))*ee
    dVy = (bb[0]*(x0-XX[0])+2*cc[0]*(y0-YY[0]))*ee
    for j in range(1,4):
        ee = AA[j]*np.exp(aa[j]*(x0-XX[j])**2 + bb[j]*(x0-XX[j])*(y0-YY[j]) + cc[j]*(y0-YY[j])**2)
        V0 = V0 + ee
        dVx = dVx + (2*aa[j]*(x0-XX[j]) + bb[j]*(y0-YY[j])  )*ee
        dVy = dVy + (  bb[j]*(x0-XX[j]) + 2*cc[j]*(y0-YY[j]))*ee

    V0=V0-np.min(V0);



    xsize = 5
    ysize = xsize*0.7/0.9
    fig = plt.figure(1,(xsize, ysize),dpi=300)
    plt.subplots_adjust(bottom=0.130,left=0.120,right=0.97,top=0.915, hspace=0.001)
    #print len(g1), len(Vz)
    plot(g1,Vz,'w-')
    plot(g0,V0,'r-')
    savefig('fep1d.png', dpi=300)
    close()

