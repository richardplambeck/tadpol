# sim.py

import math
import numpy
import matplotlib.pyplot as pyplot

# pick random point on surface of a sphere, compute PA of projection along the line of sight
# coordinate system: +Z points at observer; +X, +Y as in polarization measurements; PA = atan(y/x);
#    ... note: atan, not atan2, because we want all PAs within 180 degree span
# method (suggested by Wolfram Mathworld): generate 3 gaussian random variables x,y,z; then
#    compute 1/sqrt(x^2 + y^2 +z^2) * [x y z]

def paRandom( ) :
  x = numpy.random.normal()
  y = numpy.random.normal()
  z = numpy.random.normal()  
  length = math.sqrt(x*x + y*y + z*z)
  pa = 180./math.pi * math.atan(y/x)	# note: atan, not atan2
  return [ pa, x/length, y/length, z/length ]

def paDifRandom( outfile ) :
  fout = open( outfile, "w" )
  for n in range(0,10000) :
    [pa1,x,y,z] = paRandom()
    [pa2,x,y,z] = paRandom()
    paDif = abs(pa1 - pa2)
    if paDif > 90. :
      paDif = 180. - paDif
    fout.write("%6d  %10.3f  %10.3f  %10.3f\n" % (n, pa1, pa2, paDif) )
  fout.close() 
  
# choose 2 points randomly on surface of sphere; compute PA projected on the plane of
#  the sky for each of them; then compute angular difference of the two directions;
#  reorder and compute cumulative distribution function

def paDifUniform( outfile ) :
  ntrials = 100000
  nsamples = 1000
  pdif = numpy.zeros( ntrials, dtype=float ) 
  for n in range( 0, ntrials ) :
    [pa1,x,y,z] = paRandom()
    [pa2,x,y,z] = paRandom()
    pdif[n] = abs(pa1 - pa2)
    if pdif[n] > 90. :
      pdif[n] = 180. - pdif[n]
  pdifSorted = numpy.msort(pdif)
  nskip = ntrials/nsamples
  fout = open( outfile, "w" )
  for n in range( 0, ntrials, nskip ) :
    fout.write("%10.3f  %10.3f\n" % (pdifSorted[n], float(n)/float(ntrials) ) )
  fout.write("%10.3f  %10.3f\n" % (pdifSorted[ntrials-1], 1.0 ) )

# K-S statistic: 
# N = number of data points (> 20 desirable)
# D = max deviation from theoretical distribution
# probKS is the probability that D > observed value if the observed data
#   were drawn from the theoretical distribution; small values of probKS
#   mean that the data differ significantly from the model
# see Numerical Recipes p 474

def probKS( N, D ) :
  x = math.sqrt(N) * D
  sum = 0.
  for n in range(1, 100) :
    sum = sum + 2.* pow(-1,(n-1)) * math.exp(-2.*n*n*x*x)
    print n, sum
  
# choose first point randomly on surface of sphere; keep choosing second point
#   randomly until it lies between thetaMin and thetaMax of first point;
#   compute position angles  pa1 and pa2 of these 2 vectors projected on the 
#   plane of the sky; compute angular difference between them; then reorder
#   and generate the cumulative distribution function
# thetaMin is minimum angle between vector1 and vector2
# thetaMax is maximum angle between vector1 and vector2
# ntrials is the number of trials
# nbins is the number of bins for the cumulative distribution function
#   (it could be as large as ntrials but then wip couldn't plot it - wip
#   will handle at most 20000 points)

def paDifCone( thetaMin, thetaMax, outfile, ntrials=100000, nbins=1000 ) :
  nsamps = 0
  pdif = numpy.zeros( ntrials, dtype=float ) 
  costhetaMin = math.cos( math.pi*thetaMin/180. )
  costhetaMax = math.cos( math.pi*thetaMax/180. )
  for n in range( 0, ntrials ) :
    [pa1,x1,y1,z1] = paRandom()
    cosphi = -1.01	# impossible value
    while (cosphi > costhetaMin) or (cosphi < costhetaMax) :
      [pa2,x2,y2,z2] = paRandom()
      cosphi = x1*x2 + y1*y2 + z1*z2	# yeah, this should be dot(a,b)
      nsamps = nsamps + 1
    pdif[n] = abs(pa1 - pa2)
    if pdif[n] > 90. :
      pdif[n] = 180. - pdif[n]          # complement of angle to keep angle acute
  print "ntrials = %d, nsamps = %d" % (ntrials,nsamps)
  pdifSorted = numpy.msort(pdif)
  nskip = ntrials/nbins
  if nskip == 0 :
    nskip = ntrials
  fout = open( outfile, "w" )
  for n in range( 0, ntrials, nskip ) :
    fout.write("%10.3f  %10.3f\n" % (pdifSorted[n], float(n)/float(ntrials) ) )
  fout.write("%10.3f  %10.3f\n" % (pdifSorted[ntrials-1], 1.0 ) )
  return pdifSorted

# bimodal distribution: fraction f1 lie between thetaMin1 and thetaMax2, franction (1-f1) 
#   between thetaMin2 and thetaMax2

def bimodal( f1, thetaMin1, thetaMax1, thetaMin2, thetaMax2, outfile, ntrials=100000, nbins=1000) :
  ntrials1 = f1 * ntrials
  ntrials2 = ntrials - ntrials1
  pdif1 = paDifCone( thetaMin1, thetaMax1, 'dummy1', ntrials=ntrials1, nbins=nbins ) 
  pdif2 = paDifCone( thetaMin2, thetaMax2, 'dummy2', ntrials=ntrials2, nbins=nbins ) 
  pdifSorted = numpy.union1d(pdif1,pdif2)
  nskip = ntrials/nbins
  if nskip == 0 :
    nskip = ntrials
  fout = open( outfile, "w" )
  for n in range( 0, ntrials, nskip ) :
    fout.write("%10.3f  %10.3f\n" % (pdifSorted[n], float(n)/float(ntrials) ) )
  fout.write("%10.3f  %10.3f\n" % (pdifSorted[ntrials-1], 1.0 ) )
   
  

#  paDifCone( 180, "paDif180.dat" )		# outflows randomly oriented w/ respect to field
#  paDifCone( 20, "paDif20.dat" )		# outflows within 20 degrees of field direction
#  paDifCone2( 70, 110, "paDif90.dat" )   # outflows perpendicular to field
#  paDifCone2( 0, 20, "paDif20B.dat" )   

def examples( ) :
  bimodal( 0.5, 0., 15., 75., 105., "bimodal15", ntrials=100000, nbins=1000)


# this is a quick test of Karto's intensity mapping result where average of 3 noisy
#   channels seems to give a mean that is much more consistent than I would expect
def karto() :
  pyplot.clf()
  size = 45
  y = numpy.random.normal( size=size )     # array of 45 Gaussian numbers
  x = numpy.arange(0, size, 1)             # index ranges from 0 to size-1
  print x,y
  fig = pyplot.subplot(1,1,1)
  fig.axis( [-2,size+2,-3.,3.] )
  fig.grid( True, linewidth=0.1, color="0.15" ) 
  fig.plot(x,y,'o')
  yy = numpy.mean(numpy.reshape( y, (-1,3) ), axis=1)
  print yy 
  x = numpy.arange(1, len(y), 3)		   # new index 1,3,7... for yy
  print x
  yerr = numpy.std(yy)
  fig.errorbar(x, yy, yerr=yerr, xerr=1., color='red', fmt="+" )
  pyplot.show()

