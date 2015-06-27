# optics.py
# used to visualize Gaussian beams on front-fed or side-fed Dragone reflectors
# the reflector parameters are calculated using formulae from Chang and Prata 2004, 
#    IEEE Trans Antennas Prop., 52, 12
# ---------------------------------------------------------------------------------------------------------- #

#from Numeric import *
import numpy 
import math
import cmath
import random
import pol				# def Jrot
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import MultipleLocator
from mpl_toolkits.mplot3d import Axes3D

# -----  here are the parameters defining the setup ----
D=8.
thetaEdeg=30.
thetaCdeg=-160.  # -170.
ell=8.
alphadeg=-60.    # -60.
signFF=1. 
fGHz = 90.
F2F = 15.

# -----------------------------------------------------

def parab( ):
  plt.ion()
  fig = plt.figure( figsize=(20,10), facecolor='white', tight_layout=True )
  ax = fig.add_subplot( 1, 1, 1, projection='3d', title="parab" )
  #ax.set_axis_off()
  #ax.view_init(elev=viewEl,azim=viewAz)
  x = numpy.arange(-5,5,.05)
  y = numpy.arange(-5,5,.05)
  xx,yy = numpy.meshgrid(x, y)
  rsq = (xx**2 + yy**2)
  print rsq
  z = 0.02 * rsq 
  surf = ax.plot_surface( xx, yy, z, rstride=10, cstride=10, alpha=0.5 )
  ax.set_aspect('equal')
  ax.set_zlim(0,10.)
  plt.show()
  #plt.savefig( pp, format='pdf', bbox_inches='tight' )

# compute Dragone parameters from steps 1-5 in Chang and Prata 2004
# D = diameter of primary mirror (4*w for Gaussian beam, to avoid spillover)
# thetaE = halfwidth angle of feed pattern (2*w for Gaussian beam)
# thetaC = aiming angle of feed axis relative to final axis of beam
# ell = distance from center of subreflector to center of primary
# alpha = angle from feed axis to (hyperboloidal) subreflector focus
# signFF = 1. for front fed, -1. for side fed config
# 
# distances can be given in any units (e.g., wavelengths or inches)
# angles are in degrees (converted to radians inside this subroutine

#def dragone( D=100., thetaEdeg=30, thetaCdeg=-170., ell=100, alphadeg=-60., signFF=1.) :
def dragone( D=D, thetaEdeg=thetaEdeg, thetaCdeg=thetaCdeg, ell=ell, alphadeg=alphadeg, signFF=signFF) :
  dtor = math.pi/180.
  thetaE = thetaEdeg * dtor
  thetaC = thetaCdeg * dtor
  alpha = alphadeg * dtor

  beta = thetaC - alpha															  # eqn 45
  theta0 = beta - 2* math.atan( pow( math.tan(alpha/2.), 2.)/ math.tan( beta/2.) ) 	  # eqn 46
  e = math.sin( (beta+alpha)/2. ) / math.sin( (beta-alpha)/2. )					  # eqn (47)
  Feq = D/ (4. * math.tan( thetaE/2. ) )											# eqn (48)
  F = D * math.sin(beta) / (4. * math.sin(alpha) * math.tan(thetaE/2.) )					# eqn (49)
  c = 0.5 * ( 2.*F/(1.+math.cos(theta0)) - ell) \
	    * math.sin( thetaC - theta0 ) / math.sin( alpha )							# eqn (51)
  d0 = 2.*F * math.tan( abs(theta0)/2. )
  df = 2.*c*math.sin( abs(beta) )															# eqn (37)
  #dcf = abs(d0) - D/2. - df															# eqn (36,52)
  factor = 2.*F * math.tan( abs(theta0)/2. ) - 2 * c * math.sin( abs(beta) )      
  dcf = signFF * factor - D/2.														# eqn (54)
  dcs = factor - D/2 \
        + ( pow(e,2.) - 1.) * c/e * math.sin( abs(thetaC) + signFF*thetaE ) \
        / (e * math.cos( abs(alpha) + signFF*thetaE ) + 1.)							# eqn (55)
  
  print "beta = %.3f deg" % (beta/dtor)
  print "theta0 = %.3f deg" % (theta0/dtor)
  print "e = %.5f" % e
  print "Feq = %.5f" % Feq
  print "F = %.5f" % F
  print "2c = %.5f" % (2.*c)
  print "d0 = %.5f" % d0
  print "df = %.5f" % df
  print "dcf = %.5f" % dcf
  print "dcs = %.5f" % dcs
  return [D,thetaE,thetaC,ell,alpha,beta,theta0,e,Feq,F,c,d0,df,dcf]

# make data file of cross section of gaussian beam to plot with wip, put in talk
def gCrossSection( fGHz=70., w0in=0.3 ) :
  #fout = open("gCrossSection.dat", "w" )
  #fout.write("# output of optics.gCrosssection for fGHz = %.3f, w0in = %.3f\n" % (fGHz,w0in) )
  lambdain = (29.98/fGHz)/2.54
  print "lambda = %.3f inches" % lambdain
  for z in numpy.arange(-10.,10.05,.05) :
    w = w0in * math.sqrt( 1. + pow( lambdain*z/(math.pi * pow(w0in,2.)),2. ) )
    #fout.write("%8.3f  %8.3f  %8.3f\n" % (z, w, -1.*w) )    
  fout.close()

def parBeam( ax, w0in, d0, surf1, F2F, fGHz=fGHz ) :
  lambdain = (29.98/fGHz)/2.54
  gx = numpy.array( numpy.arange(0.,25.,.01) )
  gy = w0in * numpy.sqrt( 1. + pow( lambdain*gx/(math.pi * pow(w0in,2.)),2. ) )
  bmx = F2F/2.-gx
  bmy = d0 + gy
  nmin,imin = trunc( numpy.array( [bmx,bmy] ), surf1 )
  ax.plot( bmx[0:nmin], bmy[0:nmin], color="red", linestyle='--', linewidth=3 )
  bmy = d0 - gy
  nmin,imin = trunc( numpy.array( [bmx,bmy] ), surf1 )
  ax.plot( bmx[0:nmin], bmy[0:nmin], color="red", linestyle='--', linewidth=3 )

  gy = 2.*w0in * numpy.sqrt( 1. + pow( lambdain*gx/(math.pi * pow(w0in,2.)),2. ) )
  bmx = F2F/2.-gx
  bmy = d0 + gy
  nmin,imin = trunc( numpy.array( [bmx,bmy] ), surf1 )
  ax.plot( bmx[0:nmin], bmy[0:nmin], color="blue", linestyle='--', linewidth=3 )
  bmy = d0 - gy
  nmin,imin = trunc( numpy.array( [bmx,bmy] ), surf1 )
  ax.plot( bmx[0:nmin], bmy[0:nmin], color="blue", linestyle='--', linewidth=3 )
  
  
def gBeam( ax, x, y, theta, surf1, surf2, w0in=0.23, fGHz=fGHz ) :
  lambdain = (29.98/fGHz)/2.54
  gx = numpy.array( numpy.arange(0.,25.,.01) )

  gy = w0in * numpy.sqrt( 1. + pow( lambdain*gx/(math.pi * pow(w0in,2.)),2. ) )
  tiltbeam = pol.Jrot( numpy.array( [gx, gy] ), -1.*theta )
    # rotate the beam by angle theta
  beam = numpy.array( [tiltbeam[0]+x, tiltbeam[1]+y] )
    # offset the beam by x,y
  ybm1 = plotRay( ax, beam, surf1, surf2, color="red" )
  print "ybm1 = ", ybm1

  gy = -gy
  tiltbeam = pol.Jrot( numpy.array( [gx, gy] ), -1.*theta )
  beam = numpy.array( [tiltbeam[0]+x, tiltbeam[1]+y] )
  ybm2 = plotRay( ax, beam, surf1, surf2, color="red" )
  print "ybm2 = ", ybm2

  gy = -2.*gy
  tiltbeam = pol.Jrot( numpy.array( [gx, gy] ), -1.*theta )
  beam = numpy.array( [tiltbeam[0]+x, tiltbeam[1]+y] )
  plotRay( ax, beam, surf1, surf2, color="blue" )

  gy = -gy
  tiltbeam = pol.Jrot( numpy.array( [gx, gy] ), -1.*theta )
  beam = numpy.array( [tiltbeam[0]+x, tiltbeam[1]+y] )
  plotRay( ax, beam, surf1, surf2, color="blue" )
 
  return (ybm1-ybm2)/2.

# plot ray reflecting off 2 surfaces (or truncated at 2nd surface)
def plotRay( ax, ray, surf1, surf2, color='black', truncate=True ) :
  nmin,imin = trunc( ray, surf1 )
    # figure out where ray intersects surface1; nmin is ray index, imin is surface index
  ax.plot( ray[0][0:nmin], ray[1][0:nmin], color=color, linestyle='--', linewidth=3 )
    # draw the ray up to the first surface
  if (len(ray[0])-nmin) < 2 :
    print "ray too short"
    return
  normal = norm( None, surf1, imin )
    # norm is angle normal to the surface at point of intersection
  ray = reflected( ax, ray, nmin, normal ) 
    # new ray is the reflected ray from surface1
  nmin,imin = trunc( ray, surf2 )
    # figure out where ray intersects surface2
  ysave = ray[1][nmin]
    # may need this to get new beamwaist
  ax.plot( ray[0][0:nmin], ray[1][0:nmin], color=color, linestyle='--', linewidth=3 )
    # plot ray from surface1 to surface2
  if not(truncate) and (len(ray[0])-nmin) > 1 :
    normal = norm( None, surf2, imin )
      # norm is angle normal to surf2 at point of intersection with ray
    ray = reflected( ax, ray, nmin, normal ) 
      # new ray is the reflected ray from surface2
    ax.plot( ray[0][0:nmin], ray[1][0:nmin], color=color, linestyle='--', linewidth=3 )
  return ysave
  
def reflected( ax, beam, nmin, normal ) :
  bm = numpy.array( [beam[0][nmin:]-beam[0][nmin], beam[1][nmin:]-beam[1][nmin]] )
    # this is the reflected part of the ray, with point of reflection as the origin
  bmrot = pol.Jrot( bm, 180./math.pi*normal-90. )
    # rotate into the frame of the normal to the surface
  bm2 = pol.Jrot( numpy.array( [bmrot[0],-bmrot[1]] ), -180./math.pi*normal+90. ) 
  return numpy.array( [ bm2[0]+beam[0][nmin], bm2[1]+beam[1][nmin]] ) 
    # rotate back into plot frame, move origin to point of reflection on surface

def trunc( beam, surf ) :
  # for each point in beam, compute dsq to every point on surface
  nmin = 0
  imin = 0
  dsqmin = 1000.
  for n in range(0, len(beam[0])) :
    dsq = pow( surf[0] - beam[0][n], 2 ) + pow( surf[1] - beam[1][n], 2 )
      # distance from each point on surface to this point on beam
    i = numpy.argmin(dsq)
      # return index of minimum; this is the point on surface that is closest
    if (dsq[i] < dsqmin ) :
      dsqmin = dsq[i]
      nmin = n 
      imin = i
      #print nmin, dsqmin
  #print "returning nmin = ", nmin
  return nmin,imin
  
 
# draw optics for test setup in 2D
# F2F is spacing between foci of transmitting and receiving paraboloids
def layout2D( F2F=F2F ) :
  # pp = PdfPages( 'ComplexLeaks.pdf' )
  # pyplot.ioff()
  D,thetaE,thetaC,ell,alpha,beta,theta0,e,Feq,F,c,d0,df,dcf = dragone()
  plt.ion()
  fig = plt.figure( figsize=(20,10), facecolor='white', tight_layout=True )
  ax = fig.add_subplot( 1, 1, 1 )
  #minorLocator = MultipleLocator( 1. )
  #ax.yaxis.set_minor_locator( minorLocator )
  #ax.xaxis.set_minor_locator( minorLocator )
  ax.set_aspect('equal')
  #theta1 = -1.*theta0 - .28
  #theta2 = -1.*theta0 + .28
  theta1 = 2. * math.atan( (d0 - D/2.)/(2.*F) ) 
  theta2 = 2. * math.atan( (d0 + D/2.)/(2.*F) ) 
  print "theta1,2 = ",theta1,theta2
  pcurve = parab2D( ax, F, theta1, theta2, F2F ) 
  hcurve = hyperb2D( ax, alpha, thetaE, e, c, beta, F2F )
  focus( ax, c, beta, F2F )
  hornOutline( ax, 2.*c*math.cos(math.pi-beta), 2.*c*math.sin(math.pi+beta), 180.+thetaCdeg, F2F )
  newBW = gBeam( ax, 2.*c*math.cos(math.pi-beta), 2.*c*math.sin(math.pi+beta), thetaCdeg, hcurve, pcurve )
  parBeam( ax, newBW, d0, pcurve, F2F, fGHz=fGHz )
  #test( ax, c, e )
  #plt.grid( which='both' )
  #plt.xlim( -7., 22.)
  plt.grid( True )
  plt.show()
   
# draw dot at focus
def focus( ax, c, beta, F2F ) :
  x = 2.*c*math.cos( math.pi + beta )    # plus sign because beta is negative
  y = 2.*c*math.sin( math.pi + beta )
  ax.plot( x, y, 'rx' )
  x = F2F - x
  ax.plot( x, y, 'rx' )

def parab2D( ax, F, theta1, theta2, F2F ) :
  dTheta=0.005   # radians 
  theta = numpy.arange(theta1,theta2+dTheta, dTheta) ;
  rho = 2 * F / (1 + numpy.cos(theta) )
  #print "theta = ",theta
  #print "rho = ",rho
  x = rho * numpy.cos( math.pi - theta)
  y = rho * numpy.sin( math.pi - theta) 
  surf = numpy.array( [x,y] )
  nmidpoint = len(x)/2
     # note: points evenly spaced in thetaS, so this may not be true midpoint
  dx = 0.7 * math.cos( math.pi + norm( None, surf, nmidpoint ) )
  dy = 0.7 * math.sin( math.pi + norm( None, surf, nmidpoint ) )
  mbx = numpy.array( [x[-1]+dx, x[0]+dx] )
  mby = numpy.array( [y[-1]+dy, y[0]+dy] )
  mirrorx = numpy.concatenate( (x, mbx) )
  mirrory = numpy.concatenate( (y, mby) )
  #print "mbx = ",mirrorx
  #print "mby = ",mirrory  
  ax.fill( mirrorx, mirrory, alpha=0.5, color="blue" )
  ax.plot( x, y, color='black' )
  pcurve = numpy.array( [x,y] )
  mirrorx = F2F - mirrorx
  x = F2F - x
  surf = numpy.array( [x,y] )
  ax.plot( x, y, color='black' )
  ax.fill( mirrorx, mirrory, alpha=0.5, color="blue" )
  return pcurve
  
def hyperb2D( ax, alpha, thetaE, e, c, beta, F2F ) :
  dTheta=0.005   # radians 
  thetaBeta =  numpy.arange( alpha-thetaE, alpha+thetaE+dTheta, dTheta ) 
 	 # from Fig 3, this spans the complete subreflector
  #thetaS = 2. * numpy.arctan( (1-e)/(1+e) * numpy.tan(thetaBeta/2.) ) 
     # this is Chang and Prata version of eqn 30; it doesn't work
  thetaS =  math.pi - 2. * numpy.arctan( (1-e)/(1+e) * numpy.tan(thetaBeta/2.) )
     # this is my version of eqn 30
  rhoS = ((1 - pow(e,2.))*c/e) / (e * numpy.cos(thetaS)+ 1. )  
  #print "thetaBeta = ", (180./math.pi * thetaBeta) 
  #print "argumment = ", (1-e)/(1+e) * numpy.tan(thetaBeta/2.)  
  #print "thetaS = ",(180/math.pi * thetaS)
  #print "rhoS = ",rhoS
  angFrame = -1.*(math.pi + beta - thetaS)
	 # angle relative to paraboloid axis = x-axis of my plot
  x = -1. * rhoS * numpy.cos( angFrame )
  y = rhoS * numpy.sin( angFrame ) 
     # coordinates of surface in x-y frame
  ax.plot( x, y, color='black' )
  surf = numpy.array( [x,y] )
  nmidpoint = len(x)/2
     # note: points evenly spaced in thetaS, so this is not the true midpoint
  dx = 0.7 * math.cos( math.pi + norm( None, surf, nmidpoint ) )
  dy = 0.7 * math.sin( math.pi + norm( None, surf, nmidpoint ) )
  mbx = numpy.array( [x[-1]+dx, x[0]+dx] )
  mby = numpy.array( [y[-1]+dy, y[0]+dy] )
  mirrorx = numpy.concatenate( (x, mbx) )
  mirrory = numpy.concatenate( (y, mby) )
  ax.fill( mirrorx, mirrory, alpha=0.5, color="blue" )
  mirrorx = F2F - mirrorx
  ax.fill( mirrorx, mirrory, alpha=0.5, color="blue" )
    # draw the receiving mirror
  return surf

# finds normal to curve[n] in radians, pointing toward concave side
def norm( ax, curve, n ) :
  x = curve[0]
  y = curve[1]
  slope1 = math.atan( (y[n]-y[n-1])/(x[n]-x[n-1]) )
  slope2 = math.atan( (y[n+1]-y[n])/(x[n+1]-x[n]) )
  # sign of 2nd derivative determines concave direction
  if (slope2 -slope1)/(x[n+1]-x[n]) > 0. :
    normal = (slope1+slope2)/2. + math.pi/2.
  else :
    normal = (slope1+slope2)/2. - math.pi/2.
  #print " "
  #print "x = ",x[n-1],x[n],x[n+1]
  #print "y = ",y[n-1],y[n],y[n+1]
  #print "slope1 = ", 180./math.pi * slope1
  #print "slope2 = ", 180./math.pi * slope2
  #print "normal = ", 180./math.pi * normal
  xp = [ x[n], x[n] + math.cos(normal) ] 
  yp = [ y[n], y[n] + math.sin(normal) ] 
  if ax :
    ax.plot( xp, yp, color='red' )
  print "normal = ", 180./math.pi * normal
  return normal


# draw outline of horn and waveguide flange
# x,y are center of horn aperture
# theta = angle of horn axis relative to center of aperture
# note: theta=0 -> horn opens to the left, along the -x axis
def hornOutline( ax, x, y, theta, F2F ) :
  hx = numpy.array( [0.00, 0.00, 1.25, 2.70, 2.70, 2.80, 2.80] )
  hy = numpy.array( [0.00, 0.45, 0.17, 0.17, 0.38, 0.38, 0.00] )
  hornx = numpy.concatenate( (hx, hx) )
  horny = numpy.concatenate( (hy, -hy) )
  horn = pol.Jrot( numpy.array( [hornx, horny] ), -1.*theta )
  ax.fill( horn[0]+x, horn[1]+y, alpha=0.5, color="blue" )
  ax.fill( F2F-horn[0]-x, horn[1]+y, alpha=0.5, color="blue" )
 
def test( ax, c, e ) :
  theta = numpy.arange( 2.8, 3.5 , .05)
  rho = c*(e*e -1.)/(1.-e*numpy.cos(theta))
  x = rho * numpy.cos( theta)
  y = rho * numpy.sin( theta) 
  ax.plot( x, y, color='black' )
