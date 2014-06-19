# marsPA.py 

import math
import time
import numpy
import subprocess
import shlex
import string
import sys
from matplotlib import pyplot
from matplotlib.patches import Circle 
from matplotlib.patches import Ellipse 

# figure out deviation from radial of mars PA vectors


def readArray( imageFile, arcsecBox ) :
  # dump selected region of image to a logfile
  p = subprocess.Popen( ( shlex.split('imtab in=%s region=arcsec,box(%d) log=imtablog' % ( imageFile, arcsecBox ) ) ), \
     stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.STDOUT) 
  result = p.communicate()[0]
  print result
  # open the logfile, fill out lists
  x = []
  y = []
  z = []
  fin = open("imtablog", "r")
  for line in fin :
    a = line.split()
    x.append( -1.*float(a[0]) )   # reverse RA to get arcsec on sky
    y.append( float(a[1]) )
    z.append( float(a[2]) )
  fin.close()
  return [ numpy.array(x), numpy.array(y), numpy.array(z) ]

# read in as an n^2 x 1 array, print as n x n
def printArray( array ) :
  n = int( math.sqrt(len(array)) + .0001 )
  if pow(n,2) != len(array) :
    print "error - not a square array"
  else :
    asq = numpy.reshape(array, (n,n) )
    for i in range( n-1, -1, -1 ) :
      print " " + numpy.array_str( asq[i], precision=2, suppress_small=True, max_line_width=200 )
       
def plotArray( p, p2, array, arcsecBox, vrange=0. ) :
  nsq = int( math.sqrt(len(array)) + .0001 )
  p.tick_params( labelsize=10 )
  if vrange > 0. :
    imgplot = p.imshow(numpy.reshape(array, (nsq,nsq) ), origin='lower', \
        extent=[-arcsecBox,arcsecBox,-arcsecBox,arcsecBox], vmin=-1.*vrange,vmax=vrange )
  else : 
    imgplot = p.imshow(numpy.reshape(array, (nsq,nsq) ), origin='lower', \
        extent=[-arcsecBox,arcsecBox,-arcsecBox,arcsecBox] )
  if not p2 == None :
    pyplot.colorbar( imgplot, cax=p2  )
 
def getBeam( mapName ) :
  '''reads image file, returns bmaj,bmin,bpa'''
  imageFile = mapName+".U.cm"
  p = subprocess.Popen( ( shlex.split('imlist in=%s' % imageFile ) ), \
     stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.STDOUT) 
  result = p.communicate()[0]
  bmaj = bmin = bpa = 0.
  n = result.find( "bmaj" ) 
  if n > 0 :
    a = result[n:].split()
    b = a[2].split(":")
    print "bmaj :", b[2]
    bmaj = float( b[2] )
  n = result.find( "bmin" )
  if n > 0 :
    a = result[n:].split()
    b = a[2].split(":")
    print "bmin:", b[2]
    bmin = float(b[2])
  n = result.find( "bpa" )
  if n > 0 :
    a = result[n:].split()
    print "bpa :", a[2]
    bpa = float(a[2])
  return [ bmaj, bmin, bpa ]

def plotHisto( p, array, wgt ) :
  p.tick_params( labelsize=10 )
  n, bins, patches = p.hist( array, bins=40, range=(-20.,20.), normed=True, weights=wgt, facecolor='g')
  #pyplot.axis( [-10.,10., 0, .1] )
  #pyplot.show()


# minPoli = minimum polarized intensity in Jy/beam that will be used in the calculation
# pldiam = mean planet diameter in arcsec

def calc( mapName, arcsecBox, minPoli, pldiam=7.76 ) :
  [x,y,U] = readArray( mapName+".U.cm", arcsecBox )
  [x,y,Q] = readArray( mapName+".Q.cm", arcsecBox )
  #[x,y,pa] = readArray( "Mars.pa.cm", arcsecBox )

  #print "\nx:"
  #printArray( x )
  #print "\ny:"
  #printArray( y )

  PAradial = numpy.arctan(y/x) - math.pi/2.   # x,y are -RA,DEC coordinates on sky
  for n in range( 0, len(PAradial) ) :
    if PAradial[n] < -math.pi/2. :
      PAradial[n] = PAradial[n] + math.pi
  #print "\npa-radial:"
  #printArray( PAradial * 180./math.pi )

  #print "\npa-miriad:"
  #printArray( pa )

  #print "\nU:"
  #printArray( U )

  #print "\nQ:"
  #printArray( Q )

  PAmeas = 0.5 * numpy.arctan2(U,Q)
  #print "\nPAmeas"
  #printArray( PAmeas * 180./math.pi )

  wgt = numpy.sqrt( pow(U,2) + pow(Q,2) )
  # or should it be wgt =( pow(U,2) + pow(Q,2) )
  #print "\nwgt:"
  #printArray( wgt )

  wgtMasked = numpy.ma.masked_less(wgt, minPoli)
    # ... mask off values below threshold
  
  devDeg = ( PAmeas - PAradial ) * 180./math.pi
  for n in range( 0, len(devDeg) ) :
    if devDeg[n] > 90. :
      devDeg[n] = devDeg[n] - 180. 
    if devDeg[n] < -90. :
      devDeg[n] = devDeg[n] + 180.
    if numpy.isnan(devDeg[n]) :
      devDeg[n] = 0.
      wgt[n] = 0.
  #print "\ndevDeg:"
  #printArray( devDeg )

  devDegMasked = numpy.ma.masked_array( devDeg, mask=wgtMasked.mask )
  avgdev = numpy.ma.average(devDegMasked, weights=wgt)
  print "\nmean deviation = %.2f deg" % avgdev

  Smax = numpy.max( [ U.max(), -1.*U.min(), Q.max(), -1.*Q.min() ] ) 

  fig = pyplot.figure()
  #p = pyplot.subplot(2, 2, 1)
  p = fig.add_axes( [.08,.55,.35,.35])
  p2 = fig.add_axes([.41,.55,.015,.35])
  p2.tick_params( labelsize=10 )
  plotArray( p, p2, Q, arcsecBox, vrange=Smax )
  p.set_title( 'Stokes Q (Jy/beam)', fontsize=12)
  circ = Circle( (0.,0.), pldiam/2., color='r', fill=False)
  p.add_patch(circ)
  [ bmaj, bmin, bpa ] = getBeam( mapName )
  ell = Ellipse( (-.75*arcsecBox,-.75*arcsecBox), bmaj, bmin, angle=bpa+90., \
    edgecolor="m", label="beam", hatch="/", fill=False)
  p.add_patch(ell)
  p.set_xlabel( 'arcsec', fontsize=10 )
  p.set_ylabel( 'arcsec', fontsize=10 )

  #p = pyplot.subplot(2, 2, 2)
  p = fig.add_axes( [.55, .55,.35,.35])
  p2 = fig.add_axes( [.88, .55, .015, .35] ) 
  p2.tick_params( labelsize=10 )
  plotArray( p, p2, U, arcsecBox, vrange=Smax )
  p.set_title( 'Stokes U (Jy/beam)',fontsize=12)
  circ = Circle( (0.,0.), pldiam/2., color='r', fill=False)
  p.add_patch(circ)
  p.set_xlabel( 'arcsec', fontsize=10 )
  p.set_ylabel( 'arcsec', fontsize=10 )

  #p = pyplot.subplot(2, 2, 3)
  p = fig.add_axes( [.08,.1,.35,.35])
  p2 = fig.add_axes([.41,.1,.015,.35])
  p2.tick_params( labelsize=10 )
  plotArray( p, p2, wgt*devDegMasked, arcsecBox, vrange=10. )
  p.set_title( 'PA dev from radial (deg)', fontsize=12 )
  circ = Circle( (0.,0.), pldiam/2., color='r', fill=False)
  p.add_patch(circ)
  p.set_xlabel( 'arcsec', fontsize=10 )
  p.set_ylabel( 'arcsec', fontsize=10 )

  #p = pyplot.subplot(2, 2, 4)
  p = fig.add_axes( [.58,.1,.33,.35])
  plotHisto( p, devDegMasked, wgt )
  p.grid( True )
  #p.set_aspect( 'equal' )
  p.text(.2, .9,  "mean", horizontalalignment='center', transform=p.transAxes)
  p.text(.2, .8,  "%.2f deg" % avgdev, horizontalalignment='center', transform=p.transAxes)
  p.set_title( 'PA dev from radial', fontsize=12 )
  p.set_xlabel( 'dev from radial, degrees', fontsize=10 )
  p.set_ylabel( 'probability', fontsize=10 )

  pyplot.suptitle(mapName)
  pyplot.subplots_adjust( wspace=.4, hspace=.4 )
  pyplot.show()
