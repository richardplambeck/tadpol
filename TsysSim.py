# find mean and rms of noise power 
# assume independent samples of noise waveform have gaussian voltage distribution
# then what is the mean and variance of v*v?
# expect: mean = sigma^2; variance propto 2*sigma^2


import random
import numpy
import math


v = numpy.zeros( 100000, dtype=float )
for i in range(0, len(v)) :
  v[i] = random.gauss( 0., 2. )     # mean = 0, sigma = 1.
powerMean = numpy.average( v*v )
powerVariance = numpy.var( v*v )
print powerMean, powerVariance, math.sqrt(powerVariance)/powerMean
