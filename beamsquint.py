# beamsquint.py
# this is to help me understand Chu and Turrin 1972, IEEE AP-21, 339.

import matplotlib
matplotlib.use('GTKAgg')
import matplotlib.pyplot as pyplot
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Circle
from matplotlib.patches import Rectangle
from mpl_toolkits.mplot3d import Axes3D


import numpy
import math
# lambda is alternate way of defining simple function
sin = lambda arg: numpy.sin(arg)
cos = lambda arg: numpy.cos(arg)

      
# plot polarization vectors
def plotPolars( X, Y, Exx, Exy, Eyx, Eyy, radius=1. ) :
    pyplot.ioff()
    pp = PdfPages("polars.pdf")
    fig = pyplot.figure()
    ax = pyplot.subplot(1,1,1)
    pyplot.quiver( X, Y, Exx, Exy, pivot='mid', width=.004, edgecolor='none', headaxislength=0, headlength=0 )
    pyplot.quiver( X, Y, Eyx, Eyy, pivot='mid', width=.004, edgecolor='none', headaxislength=0, headlength=0, color='r' )
    left,right,bot,top = pyplot.axis()
    dx = right-left
    dy = top - bot
    pyplot.axis( [left-0.05*dx, right+0.05*dx, bot-0.05*dy, top+0.05*dy] )
    ax.set_aspect('equal')

  # i and j are coords of middle pixelm (if len(i) is odd), which will be center of circle
    i = len(X)/2
    j = len(X[i])/2
    circ = Circle( (X[i][j],Y[i][j]), radius, color='r', fill=False )
    ax.add_patch(circ)

    pyplot.savefig( pp, format='pdf')
    pyplot.show()
    pp.close()

# theta0 = angle (radians) between axis of feed and axis of parabalic reflector
# thetaC = half-angle subtended by reflector
# f = focal length of parabaloid

def computePolars( theta0=1., thetaC=0.5, f=5. ) :
  # begin with uniform grid across the beam from the feed (at some distance rho0 from focus)
    x = numpy.array( numpy.arange(-.9*thetaC,thetaC,.3*thetaC) )
    xprime,yprime = numpy.meshgrid( x, x)
      # produce 2D arrays containing xprime and yprime coordinates

  # convert to thetaPrime,phiPrime polar coordinates (rho=0 at focus of feed)
    thetaPrime = numpy.sqrt( pow(xprime,2) + pow(yprime,2) )
    phiPrime = numpy.arctan2( yprime,xprime)

  # compute Ex and Ey in reflected signal
    Exx = numpy.ones( [len(x),len(x)] )
    Exy = numpy.ones( [len(x),len(x)] )
    Xproj = numpy.ones( [len(x),len(x)] )
    Yproj = numpy.ones( [len(x),len(x)] )

    for i in range(0,len(thetaPrime)) :
      for j in range(0,len(thetaPrime)) :
        thetaP = thetaPrime[i][j]
        phiP = phiPrime[i][j]
        rho = 2*f / ( 1. + math.cos(thetaP)*math.cos(theta0) - \
          math.sin(thetaP)*math.sin(theta0)*math.cos(phiP) )
            # eqn (12)

      # X-polarized feed
        EthetaP = math.cos(phiP)/rho
        EphiP = -math.sin(phiP)/rho

        t = 1 + math.cos(thetaP)*math.cos(theta0) - math.sin(thetaP)*math.sin(theta0)*math.cos(phiP)
          # eqn (4)

      # projected locations of points in aperture plane, eqns 12,10,11
        Xproj[i][j] = rho * (math.cos(theta0)*math.sin(thetaP)*math.cos(phiP) \
	      + math.sin(theta0)*math.cos(thetaP) )
        Yproj[i][j] = rho * math.sin(thetaP) * math.sin(phiP) 

        Exx[i][j] = ((sin(thetaP)*sin(theta0) - cos(phiP)*(1.+cos(thetaP)*cos(theta0)) ) * EthetaP \
          + sin(phiP)*(cos(thetaP) + cos(theta0))*EphiP) / t
        Exy[i][j] = ( -1.*sin(phiP)*(cos(thetaP) + cos(theta0)) * EthetaP  \
          + (sin(thetaP)*sin(theta0) - cos(phiP)*(1.+cos(thetaP)*cos(theta0)) ) * EphiP) / t

    print Exx
    print Exy
    Eyx = numpy.ones( [len(x),len(x)] )
    Eyy = numpy.ones( [len(x),len(x)] )

    for i in range(0,len(thetaPrime)) :
      for j in range(0,len(thetaPrime)) :
        thetaP = thetaPrime[i][j]
        phiP = phiPrime[i][j]
        rho = 2*f / ( 1. + math.cos(thetaP)*math.cos(theta0) - \
          math.sin(thetaP)*math.sin(theta0)*math.cos(phiP) )
            # eqn (12)

      # Y-polarized feed
        EthetaP = math.sin(phiP)/rho
        EphiP = math.cos(phiP)/rho

        t = 1 + math.cos(thetaP)*math.cos(theta0) - math.sin(thetaP)*math.sin(theta0)*math.cos(phiP)
          # eqn (4)

        Eyx[i][j] = ((sin(thetaP)*sin(theta0) - cos(phiP)*(1.+cos(thetaP)*cos(theta0)) ) * EthetaP \
          + sin(phiP)*(cos(thetaP) + cos(theta0))*EphiP) / t
        Eyy[i][j] = ( -1.*sin(phiP)*(cos(thetaP) + cos(theta0)) * EthetaP  \
          + (sin(thetaP)*sin(theta0) - cos(phiP)*(1.+cos(thetaP)*cos(theta0)) ) * EphiP) / t


    rProj = 2*f*math.sin(thetaC)/( math.cos(thetaC) + math.cos(theta0) )
    print "rProj = %.3f" % rProj
    plotPolars( Xproj, Yproj, Exx, Exy, Eyx, Eyy, radius=rProj )

def Er( thetaP, phiP, theta0=0., f=1. ) :
    t = 1. + cos(thetaP)*cos(theta0) - sin(thetaP)*sin(theta0)*cos(phiP)
    rho = 2*f / ( 1. + cos(thetaP)*cos(theta0) - sin(thetaP)*sin(theta0)*cos(phiP) )
    EthetaP = cos(phiP)/rho
    EphiP = -1.*sin(phiP)/rho
    ErX = ((sin(thetaP)*sin(theta0) - cos(phiP)*(1.+cos(thetaP)*cos(theta0)) ) * EthetaP \
      + sin(phiP)*(cos(thetaP) + cos(theta0))*EphiP) / t
    ErY = ( -1.*sin(phiP)*(cos(thetaP) + cos(theta0)) * EthetaP  \
      + (sin(thetaP)*sin(theta0) - cos(phiP)*(1.+cos(thetaP)*cos(theta0)) ) * EphiP) / t
    print "rho = %.3f, Er = %.5f X + %.5f Y" % (rho,ErX,ErY)

    
def p3D( ) :
   npnts = 256
   x1 = numpy.array( numpy.arange(-1.,1.01,.01) )
   x,y = numpy.meshgrid( x1, x1 )     
   z = x * x + y * y
   fig = pyplot.figure()
   ax = fig.add_subplot(111, projection='3d')
   ax.plot_wireframe( x, y, z, rstride=10, cstride=10 )
   ax.set_aspect('equal')
   ax.view_init(elev=0., azim=0.)
   pyplot.show()

