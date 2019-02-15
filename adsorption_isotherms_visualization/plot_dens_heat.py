import os
import glob
from cStringIO import StringIO
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
# import matplotlib.cm as cmx
# import matplotlib.colors as clrsm

markers = ['','--','-.',':','-','--']
clrs = ['k','r','g','b','y','c','m']

#########################################################################  
##  Physical Constants
#########################################################################  

kb = 1.3806488E-23 #J/K
Na = 6.0221413E+23 #mol-1
epsff = 1.0 #94.45 #K   # N2
Patm = 10.0**5 #Pa

sgmffA = 3.81 # in A for CH4
sgmsfA = 3.60 # in A for CH4
sgmssA = 3.38 

# % rhoSTP:  0.0441128644794  mol/l
# % 1/rhoSTP:  22.6691241161  l/mol
rhoSTP = 0.0441128644794 #mol/l
vSTP = 22.6691241161  #l/mol




def loadcv(f):
  linestring = open(f, 'r').read()
  strlist = linestring.split('\n')
  s = "\n".join(strlist)
  fs = StringIO(s)
  t = np.loadtxt(fs)
  return t


if(__name__=='__main__'):

#loading isotherms for different pores into a ddt list
  files = glob.glob('data/density*.dat') 
  ddt = []
  for f in files:
    print "file:",f
    t = loadcv(f)
    ddt.append((t,f))


  def cmpW(x,y):
    Hx = (x[0])[0,19]
    Hy = (y[0])[0,19]
    if(Hx>Hy):
      return 1
    elif(Hx<Hy):
      return -1
    else:
      return 0

  def cmpWM(x,y):
    return -cmpW(x,y)

  ddt.sort(cmp=cmpW) #sorting the isotherms according the width of the pores




#########################################################################  
##  Main data processing loop
#########################################################################  

  r_a = []
  lgP_a = []
  dens_a = []
  lgPsel = [-2, -1, 0, 1, 2, 3, 0, 1.477121, 1.544068, 1.812913]
  lgPsel2 = [0.0, 1.477121] + [i*0.025-2 for i in range(201)] 
  Rasel_a = []
  Rasel2_a = []
  fi=1
  for t,f in ddt:
    factor = kb*Na*epsff/1000
    #all the values extraced below are arrays, so, the operations are vectorized
    temp = t[:,18]  
    rhob = t[:,17] #gas density in the pore isotherm data
    rhoav = t[:,2] 
    etaav = t[:,3]
    bPb = t[:,16]
    bUb = t[:,15]
    Nfa = t[:,14]
    bUa = t[:,13]
    R = t[:,19]
    H = 2*R
    hb = temp*(bUb+bPb)/rhob

    print 'r', R[0]

    P = temp*kb*bPb/(1e-10)**3  # Convert Pressure to Pa   # 1A = 1e-10m


    H = H - sgmssA  #width of the pore
    R = R - sgmssA/2


    qmax = 30.0
    qmin = 10.0

    lgP = np.log10(P/Patm) # calculating logarithm of the pressure


    lgPmin = -2
    lgPmax = 3.00
    lgP_ = np.linspace(lgPmin, lgPmax, 100) # intepolating the logarithm of pressure into more (100) data points
    P_ = 10**(lgP_) * Patm
    bPb_ = P_*(1e-10)**3 / (temp[0]*kb)


    dens = Nfa/(np.pi*R**2) # calculating density of the gasmolecules in the cylindrical pore
    dens = dens*sgmffA**3 # converting to dimensionles value

    densf = interp1d(lgP, dens) # creating scipy function for interpolation to new pressure values
    dens_ = densf(lgP_) # interpolating density into 100 data points
    Rasel = densf(lgPsel) # interpolating density to get values at several selected pressures
    Rasel2 = densf(lgPsel2) # interpolating density to get values at 200 selected pressures


    Rf = interp1d(lgP, R) 
    R_ = Rf(lgP_)

    r_a.append(R_)
    lgP_a.append(lgP_)  # we will store the interpolated arrays for later
    dens_a.append(dens_)# and later we will vertically stack them with tomake 2D-arrays
    Rasel_a.append(Rasel)
    Rasel2_a.append(Rasel2)




  #########################################################################  
  ##  Optimal pore size plot
  #########################################################################  
  r = np.array([r[0] for r in r_a])
  Rasel2m = np.vstack(tuple(Rasel2_a)) # vertically stacking the densities of the isotherms to get 2D-array
  Rasel2m = Rasel2m*10**27/(Na*sgmffA**3) # converting to mol/l

  # calculating optimal pore size by finding maximal delivered density
  lgP0, Rasel0 = lgPsel2[0], Rasel2m[:,0]
  optsizes = []
  i = 0
  for lgP in lgPsel2:
    if(lgP>0):
      val = Rasel2m[:,i]-Rasel0
      sel = (r>2.4)*(r<20)
      val = val[sel]
      rsel = r[sel]
      optsize=rsel[np.argmax(val)]
      optsizes.append((lgP,10**lgP,optsize))
    i+=1
  np.savetxt('Optsizes.dat',optsizes)
  
  optsizesa = np.array(optsizes)

  #plotting the calculated optimal size of the pore
  xsize = 5
  ysize = xsize*0.7/0.9
  fig = plt.figure(1,(xsize, ysize),dpi=300)
  plt.subplots_adjust(bottom=0.130,left=0.120,right=0.97,top=0.915, hspace=0.001)
  plt.semilogx(optsizesa[3:-1,1], optsizesa[3:-1,2], 'k-', ms=4)

  plt.xlabel('$P$ $[bar]$')
  plt.ylabel('Optimal Pore Radius, $R$ $[\\AA]$')
  v = plt.axis()
  xmin, xmax, ymin, ymax = v[0],v[1],v[2],v[3]
  plt.axis([-0, 3, ymin, ymax]) #R=3
  plt.axis([1, 1000, 1, 10]) #R=3

  savefig('optsize.png', dpi=300)
  close()





  #########################################################################  
  ##  Colour-contour plot of the delivered methane density
  #########################################################################  


  r = np.array([r[0] for r in r_a])
  lgP = lgP_a[0]
  densm = np.vstack(tuple(dens_a)) # vertically stacking separate isotherms into a 2D-array
  
  densm = densm*(10.0**27/(Na*sgmffA**3))  #mol/cm^3


  lgPm,rm = meshgrid(lgP, r)
  sel = lgP >= -0.1
  rm = rm[:,sel]
  lgPm = lgPm[:,sel]
  densm = densm[:,sel]
  densm = densm[:,:]-densm[:,0:1] #subtracting density at pressure 1 atm to get delivered density
 

  #plotting thedelivered density as a 2D color-contour plot
  xsize = 5
  ysize = xsize*0.7/0.9
  fig = plt.figure(1,(xsize, ysize),dpi=300)
  plt.subplots_adjust(bottom=0.130,left=0.120,right=0.97,top=0.915, hspace=0.001)

  pcolor(lgPm, rm, densm, vmin=0, vmax=25, cmap='gist_earth') #cmap='bwr')
  colorbar()

  # setting contour isovalue lines
  CS = plt.contour(lgPm, rm, densm, [1, 2, 5, 10, 15, 20, 25 ],
        colors=['w', 'w', 'w', 'w', 'w', 'w', 'w'],
                     )
  plt.clabel(CS, inline=1, fontsize=8, fmt="%2.0f")

  lgPo = optsizesa[3:-1,0]
  Hopt = optsizesa[3:-1,2]
  sel = (Hopt >= 3) * (Hopt < 20)
  plt.plot(lgPo[sel], Hopt[sel], 'w-', lw=2) # plotting optimal size on top of color map plot

  plt.ylabel('Pore Radius, $R$ $[\\AA]$')
  plt.xlabel('Bulk Pressure, $P_b$ $[bar]$')
  plt.title('Delivered Dens. $\\rho_s-\\rho_d$ $[mol/dm^3]$')
  v = plt.axis()
  plt.axis([0, 3, 1.975, 11])

  plt.axes().xaxis.set_major_locator(MaxNLocator(3))
  locs, labels = plt.xticks()
  new_labels = ['10$^{%i}$'%int(i) for i in locs]
  xticks(locs, new_labels)

  plt.text( -0.5 ,11.25, 'b)', fontsize=16)
  savefig('deliver_dens_heat.png', dpi=300)
  close()








  #########################################################################  
  ##  Plotting several intersections of the above plot at fixed pressures
  #########################################################################  

  xsize = 4
  fig = plt.figure(1,(xsize, xsize*0.56/0.9))
  plt.subplots_adjust(bottom=0.100,left=0.25,right=0.84,top=0.95, hspace=0.001)

  r = np.array([r[0] for r in r_a])
  Raselm = np.vstack(tuple(Rasel_a))
  Raselm = Raselm*10**27/(Na*sgmffA**3) 

  new_labels = ['$P_b=10^{%i}$ bar'%(i) for i in lgPsel]

  plt.plot(r, Raselm[:,3]-Raselm[:,6], 'k', label=new_labels[3], ms=4)
  plt.plot(r, Raselm[:,4]-Raselm[:,6], 'r', label=new_labels[4], ms=4)
  plt.plot(r, Raselm[:,7]-Raselm[:,6], 'b', label=new_labels[6], ms=4)
  plt.plot(r, Raselm[:,8]-Raselm[:,6], 'c', label=new_labels[8], ms=4)
  plt.plot(r, Raselm[:,9]-Raselm[:,6], 'm', label=new_labels[9], ms=4)

  val = Raselm[:,7]-Raselm[:,6]
  sel = (r>3)*(r<10)
  val = val[sel]
  rsel = r[sel]
  print 'Optimum R=', rsel[np.argmax(val)]

  plt.ylabel('Deliv. Dens. $\\rho_s-\\rho_d$ $[mol/dm^3]$', fontsize=11)
  v = plt.axis()
  xmin, xmax, ymin, ymax = v[0],v[1],v[2],v[3]
  xmin, xmax, ymin, ymax = 2, 11, ymin, 20
  plt.axis([xmin, xmax, ymin, ymax]) #R=3

  plt.text( -0.75,19.5, 'f)', fontsize=16)
  plt.plot([6.3457143,6.3457143], [-5000,5000], 'k:')

  ax1 = plt.gca()
  ax2 = ax1.twinx()
  ax2.set_ylabel('$v_{STP}/v$', color='k')
  ax2.tick_params('y', colors='k')
  plt.axis([xmin, xmax, ymin, ymax]) #R=3
  ax2.set_ylim(0, vSTP*ymax)

  savefig('DRplot.png', dpi=300)
  close()



