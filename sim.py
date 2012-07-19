# sim.py

import math
import numpy

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
#   randomly until it lies within theta of first point; compute position angles
#   pa1 and pa2 of these 2 vectors projected on the plane of the sky; compute
#   angular difference between them; then reorder and generate the cumulative
#   distribution function
# thetaDegrees is half angle of cone around first point
# ntrials is the number of trials
# nbins is the number of bins for the cumulative distribution function
#   (it could be as large as ntrials but then wip couldn't plot it - wip
#   will handle at most 20000 points)

def paDifCone( thetaDegrees, outfile, ntrials=100000, nbins=1000 ) :
  nsamps = 0
  pdif = numpy.zeros( ntrials, dtype=float ) 
  costheta = math.cos( math.pi*thetaDegrees/180. )
  for n in range( 0, ntrials ) :
    [pa1,x1,y1,z1] = paRandom()
    cosphi = -1.01	# impossible value
    while cosphi < costheta :
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

def examples( ) :
  paDifCone( 180, "paDif180.dat" )		# outflows randomly oriented w/ respect to field
  paDifCone( 20, "paDif20.dat" )		# outflows within 20 degrees of field direction
