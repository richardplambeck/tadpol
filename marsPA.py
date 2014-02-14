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
       
def plotArray( array, arcsecBox ) :
  nsq = int( math.sqrt(len(array)) + .0001 )
  imgplot = pyplot.imshow(numpy.reshape(array, (nsq,nsq) ), origin='lower', \
        extent=[-arcsecBox,arcsecBox,-arcsecBox,arcsecBox] )
  #circ = Circle( (0.,0.), 4., linestyle='-', color='r', fill=False)
  #imgplot.add_patch(circ)
  pyplot.show()
 
def plotHisto( array, wgt ) :
  n, bins, patches = pyplot.hist( array, bins=100, normed=True, weights=wgt, facecolor='g')
  #pyplot.axis( [-10.,10., 0, .1] )
  pyplot.show()


def calc( arcsecBox ) :
  [x,y,U] = readArray( "Mars.U.cm", arcsecBox )
  [x,y,Q] = readArray( "Mars.Q.cm", arcsecBox )
  #[x,y,pa] = readArray( "Mars.pa.cm", arcsecBox )

  print "\nx:"
  printArray( x )
  print "\ny:"
  printArray( y )

  PAradial = numpy.arctan(y/x) - math.pi/2.   # x,y are -RA,DEC coordinates on sky
  for n in range( 0, len(PAradial) ) :
    if PAradial[n] < -math.pi/2. :
      PAradial[n] = PAradial[n] + math.pi
  print "\npa-radial:"
  printArray( PAradial * 180./math.pi )

  #print "\npa-miriad:"
  #printArray( pa )

  print "\nU:"
  printArray( U )

  print "\nQ:"
  printArray( Q )

  PAmeas = 0.5 * numpy.arctan2(U,Q)
  print "\nPAmeas"
  printArray( PAmeas * 180./math.pi )

  wgt = numpy.sqrt( pow(U,2) + pow(Q,2) )
  print "\nwgt:"
  printArray( wgt )

  devDeg = ( PAmeas - PAradial ) * 180./math.pi
  for n in range( 0, len(devDeg) ) :
    if devDeg[n] > 90. :
      devDeg[n] = devDeg[n] - 180. 
    if devDeg[n] < -90. :
      devDeg[n] = devDeg[n] + 180.
    if numpy.isnan(devDeg[n]) :
      devDeg[n] = 0.
      wgt[n] = 0.
  print "\ndevDeg:"
  printArray( devDeg )

  #plotArray( x, arcsecBox )
  #plotArray( y, arcsecBox )
  plotArray( wgt, arcsecBox )
  #plotArray( PAradial, arcsecBox )
  #plotArray( devDeg, arcsecBox )
  plotArray( wgt*devDeg, arcsecBox )
  plotHisto( devDeg, wgt )

  
  print "\nmean deviation = ", numpy.ma.average(devDeg, weights=wgt)
