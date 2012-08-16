import numpy
import pylab

xymax = 0.18
fmin = 210.
fmax = 235.


def readdata( infile, xcol, ycol ) :
  x = []
  y = []
  fin = open( infile, "r" )
  for line in fin :
    a = line.split()
    fmean = 0.5 * (float(a[0]) + float(a[1]))
    if (fmean > fmin) and (fmean < fmax) :
      x.append( a[xcol-1] )
      y.append( a[ycol-1] )
  fin.close()
  return [ numpy.array(x), numpy.array(y) ]
    
for n in range(1,16,1) :
  #pylab.figure(n)
  pylab.subplot(3,5,n)
  pylab.axis( "scaled" )
  pylab.axis( [-1.*xymax, xymax, -1.*xymax, xymax], fontsize=5 )
  infile = "/fringe2/plambeck/pol/leakage/12aug2012/lk%d" % n
  [x,y] = readdata( infile, 7, 8 )  
  pylab.plot( x, y )

  infile = "/fringe2/plambeck/c0913/28jul2012/lk%d" % n
  [x,y] = readdata( infile, 7, 8 )  
  pylab.plot( x, y, color="r" )

  #pylab.axhline( linewidth=1 )
  #pylab.axvline( linewidth=1 )
  pylab.grid()
  pylab.title( ("C%d" % n) )

pylab.show()
