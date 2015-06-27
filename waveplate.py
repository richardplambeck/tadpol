# ---------------------------------------------------------------------------------------------------------- #
# waveplate.py

#from Numeric import *
import numpy 
import math
import cmath
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# ---------------------------------------------------------------------------------------------------------- #
# define basis vectors in reference coordinate system 
# each vector consists of 2 complex numbers [ (Re_x,Im_x), (Re_y,Im_y) ]
# R and L propagate in +Z direction according to right hand rule
# note: x.conjugate() gives complex conjugate of a number
# ---------------------------------------------------------------------------------------------------------- #
X = numpy.array( [1,0], dtype=complex )
Y = numpy.array( [0,1], dtype=complex )
R = numpy.array( [1./math.sqrt(2.), -1j*(1./math.sqrt(2.))] ) 
L = numpy.array( [1./math.sqrt(2.), +1j*(1./math.sqrt(2.))] ) 
clight = 29.9792458	# speed of light, cm/nanosec

# ---------------------------------------------------------------------------------------------------------- #
# returns new basis vector in a coordinate system rotated by thetaDegrees
# ---------------------------------------------------------------------------------------------------------- #
def Jrot ( vec, thetaDegrees ) :
  rad = math.pi * thetaDegrees/180.
  rotmat = numpy.array( [[math.cos(rad),math.sin(rad)],[-math.sin(rad),math.cos(rad)]] )
  return numpy.dot(rotmat,vec)

# ---------------------------------------------------------------------------------------------------------- #
# returns basis vector after passing through birefringent dielectric
# component 2 (Y-axis) is advanced by delayDegrees relative to component 1 (X-axis)
# ---------------------------------------------------------------------------------------------------------- #
def Jdelay ( vec, delayDegrees ) :
  rad = math.pi * delayDegrees/180.
  rotmat = numpy.array( [ [1, 0], [0, cmath.exp(-1.j * rad)] ] )
  return numpy.dot( rotmat,vec )

# --- delays both components --- #
def J2delay ( vec, delayDegrees ) :
  rad = math.pi * delayDegrees/180.
  rotmat = numpy.array( [ [cmath.exp(-1.j * rad), 0], [0, cmath.exp(-1.j * rad)] ] )
  return numpy.dot( rotmat,vec )

# ---------------------------------------------------------------------------------------------------------- #
# compute phase delay in degrees of Y-pol vs X-pol signals traveling through thickness tmm of birefringent 
#   dielectric; dimensions in mm 
# ---------------------------------------------------------------------------------------------------------- #
def dphi( tcm, fGHz=150., no = 3.123, ne=3.389 ) :
	# default refractive indices taken from Savini2006
  dphi = 360. * tcm * fGHz / clight * (no - ne)
  return dphi

# 5plates: 7, 36, 102, 36, 7 degrees
 
  
def nsapphire( fGHz ) :
  no = 3.053 + (4.7e-4 * fGHz) + (2.2e-10 * pow(fGHz,2.)) + (1.1e-12 * pow(fGHz,3.))
  ne = 3.387 + 1.3e-5 * fGHz
  return [no,ne]

def doit2() :
  stack2 = [ [22.5, .37], [22.5, .37] ]
  stack3 = [ [15., .37], [75., .37], [15.,.37]]
  #v2 = propagate( 150., X, stack1 )
  makepath( 150., X, stack1 )
  #print "Stokes for X: ", stokes(X)
  #print "Stokes for Y: ", stokes(Y)
  #print "Stokes for R: ", stokes(R)
  #print "Stokes for L: ", stokes(L)
  #print "Stokes for rotated R: ", stokes( Jrot( R, 45 ))

# propagate polarization state v1 through waveplate 'plate'
# v1 is a complex matrix in the X and Y basis
# plate is a list of angles (relative to the X axis) and thicknesses
def propagate( fGHz, vin, stack ) :
  v1 = vin
  for plate in stack :
    print "propagate through plate with angle %.2f degrees, thickness %.3f cm" % (plate[0],plate[1])
    v2 = Jrot( v1, plate[0] )
    v3 = Jdelay( v2, dphi( plate[1], fGHz ) )
    v1 = Jrot( v3, -plate[0])
  return v1
  

# Stokes Q,U,V are x,y,z of the poincare sphere!
def stokes( vin ) :
  Ex = vin[0]
  Ey = vin[1]
  Q = (Ex * Ex.conjugate() - Ey * Ey.conjugate())
  U = (Ex * Ey.conjugate() + Ex.conjugate() * Ey)
  V = (cmath.exp(-1j * math.pi/2.) * ( Ex * Ey.conjugate() - Ex.conjugate() * Ey ))
  return numpy.real( [Q,U,V] )

# poincare plot

# generates array [x,y,z] = [Q,U,V] showing path from start to final 
def makepath( fGHz, vin, stack ) :
  path = []
  v1 = vin
  for plate in stack :
    print v1
    v2 = Jrot( v1, plate[0] )                   # v2 is rotated starting point for this segment
    for t in numpy.arange( 0., plate[1], .005 ) :
      v3 = Jdelay( v2, dphi(t, fGHz) )      
      v4 = Jrot( v3, -plate[0] )	            # v4 is endpoint for thickness t
      path.append( stokes( v4 ) )
    v3 = Jdelay( v2, dphi( plate[1], fGHz ) )  	
    v1 = Jrot( v3, -plate[0] )                  # v1 is endpoint for thickness tplate[1]
  path.append( stokes( v1 ) )				    # endpoint of entire stack
  return numpy.array(path)

def makeendpoints( fGHz, vin, stack ) :
  path = []
  v1 = vin
  for plate in stack :
    v2 = Jrot( v1, plate[0] )                   # v2 is rotated starting point for this segment
    v3 = Jdelay( v2, dphi( plate[1], fGHz ) )  	
    v1 = Jrot( v3, -plate[0] )                  # v1 is endpoint for thickness tplate[1]
    path.append( stokes( v1 ) )				    # endpoint of entire stack
  return numpy.array(path)

def makeequator( ax ) :
  npts = 256
  x = numpy.zeros( npts )
  y = numpy.zeros( npts )
  z = numpy.zeros( npts )
  n = 0
  for phi in numpy.arange( 0., 2*math.pi, 2*math.pi/npts ) :
    x[n] = math.cos(phi)
    y[n] = math.sin(phi)
    n = n+1
  q = Axes3D.plot(ax, x, y, zs=z, c='gray' )
  q = Axes3D.plot(ax, x, z, zs=y, c='gray' )
  q = Axes3D.plot(ax, z, y, zs=x, c='gray' )

def makeaxes( ax ) :
  x = numpy.array( [-1, 1] )
  y = numpy.array( [0,0] )
  q = Axes3D.plot(ax, x, y, zs=y, c='black' )
  q = Axes3D.plot(ax, y, x, zs=y, c='black' )
  q = Axes3D.plot(ax, y, y, zs=x, c='black' )
  
def addcurve( ax, path, endpoints, color ) :
  px = path[ :,0 ]
  py = path[ :,1 ]
  pz = path[ :,2 ]
  n = len(px) - 1
  q = Axes3D.plot(ax, px, py, zs=pz, c=color, linewidth=2 )
  px = endpoints[ :,0 ]
  py = endpoints[ :,1 ]
  pz = endpoints[ :,2 ]
  print px, py, pz
  q = Axes3D.scatter(ax, px,py, zs=pz, c=color, marker='o', s=60)
     
def pplot( p1, p2, p3, e1, e2, e3 ) :
  fig = plt.figure()
  ax = Axes3D(fig)
  ax.set_aspect('equal')
  ax.set_axis_off()
  makeequator( ax )
  makeaxes( ax )
  addcurve( ax, p1, e1, 'red' )
  addcurve( ax, p2, e2, 'green' )
  addcurve( ax, p3, e3, 'blue' )
  ax.view_init(elev=15,azim=160)
  plt.show()

def doit( thetadeg ) :
  theta = thetadeg * math.pi/180.
  vin = numpy.array( [math.cos(theta),math.sin(theta)], dtype=complex )
  tt = .37
  stack1 = [ [45.,tt] ]
  stack2 = [ [22.5,tt], [67.5,tt] ]
  stack3 = [ [15.,tt],  [75.,tt],   [15.,tt] ]
  stack5 = [ [7.,tt],   [36.,tt],  [102.,tt], [36.,tt], [7.,tt] ]
  which = stack3
  p1 = makepath( 130., vin, which )
  e1 = makeendpoints( 130., vin, which )
  p2 = makepath( 150., vin, which )
  e2 = makeendpoints( 150., vin, which )
  p3 = makepath( 170., vin, which )
  e3 = makeendpoints( 170., vin, which )
  pplot( p1, p2, p3, e1, e2, e3 )
  
def drawcylinder( ax ) :
  npts = 2048 
  xa = numpy.array( [-1, 1] )
  ya = numpy.array( [0.,0.] )
  za = numpy.array( [40.,40.] )
  q = Axes3D.plot(ax, xa, za, zs=ya, c='gray' )
  q = Axes3D.plot(ax, ya, za, zs=xa, c='gray' )
  za = numpy.array( [-40.,-40.] )
  q = Axes3D.plot(ax, xa, za, zs=ya, c='gray' )
  q = Axes3D.plot(ax, ya, za, zs=xa, c='gray' )
  xa = numpy.array( [0.,0.] )
  za = numpy.array( [0.,0.] )
  ya = numpy.array( [-40.,40.] )
  q = Axes3D.plot(ax, xa, ya, zs=za, c='gray', linestyle='--' )

  x = numpy.zeros( npts )
  y = -30.*numpy.ones( npts )
  z = numpy.zeros( npts )
  xx = numpy.zeros( npts )
  yy = numpy.zeros( npts )
  zz = numpy.zeros( npts )

  n = 0
  for phi in numpy.arange( 0., 2*math.pi, 2*math.pi/npts ) :
    x[n] = math.cos(phi)
    z[n] = math.sin(phi)
    n = n+1
  #q = Axes3D.plot(ax, x, y, zs=z, c='gray' )
  y = 30.* numpy.ones( npts )
  #q = Axes3D.plot(ax, x, y, zs=z, c='gray' )
  n = 0
  for phi in numpy.arange( -30., 30.0, (60./npts) ) :
    yy[n] = phi
    amp = math.cos(phi)
    zz[n] = amp * math.cos( .007 * (30.-phi) )
    xx[n] = -1.* amp * math.sin( .007 * (30.-phi) )
    #zz[n] = math.cos(phi)
    n = n+1
  q = Axes3D.plot(ax, xx, yy, zs=zz, c='black' )
  xxx = numpy.zeros( 2 )
  yyy = 40.*numpy.ones( 2 )
  zzz = numpy.zeros( 2 )
  zzz[0] = -1.
  zzz[1] = 1. 
  q = Axes3D.plot(ax, xxx, yyy, zs=zzz, c='black' )
  xxxx = numpy.zeros( 2 )
  yyyy = -40.*numpy.ones( 2 )
  zzzz = numpy.zeros( 2 )
  xxxx[0] = -1.* math.sin( .007 * (60.) )
  zzzz[0] = math.cos( .007 * (60.) )
  xxxx[1] = 1.* math.sin( .007 * (60.) )
  zzzz[1] = -1.* math.cos( .007 * (60.) )
  q = Axes3D.plot(ax, xxxx, yyyy, zs=zzzz, c='black' )
  

def cylinder() :
  fig = plt.figure()
  ax = Axes3D(fig)
  #ax.set_aspect('equal')
  #ax.set_axis_off()
  drawcylinder( ax )
  ax.auto_scale_xyz( [-5,5], [-40,40], [-5,5] )
  ax.set_aspect('equal')
  ax.view_init(25,-75)
  #ax.set_axis_off()
  #q = fig.gca(projection='3d')
  #q._axis3don=False
  plt.show()
