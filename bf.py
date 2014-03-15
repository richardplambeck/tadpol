# bf.py
# is it possible to phase LSB and USB,then add the 2 voltages together, without losing s/n ?

import math
import numpy

snr = 0.1  # snr for 1 and 2; 
snr3 = 1.   # snr for 3

NSAMP = 10000

# there are 4 nyquist intervals per walsh sequence
# voltages from ants 1,2,3 in each of 4 nyquist intervals
v1 = numpy.zeros( 4, dtype=complex)   
v2 = numpy.zeros( 4, dtype=complex)
v3L = numpy.zeros( 4, dtype=complex)
v3U = numpy.zeros( 4, dtype=complex)

xL = numpy.zeros( 4, dtype=complex)
xU = numpy.zeros( 4, dtype=complex)
xLmeas = numpy.zeros( 4, dtype=complex)
xUmeas = numpy.zeros( 4, dtype=complex)
xL2meas = numpy.zeros( 4, dtype=complex)
xU2meas = numpy.zeros( 4, dtype=complex)

# Walsh LO1 phase switch sequence for each of these 4 nyquist intervals
w1 = numpy.array( [ 1, 0+1j, 1,    0+1j ] )
w2 = numpy.array( [ 1, 1   , 0+1j, 0+1j ] )
w1 = numpy.array( [1,1,1,1], dtype=complex )
w2 = numpy.array( [1,1,1,1], dtype=complex )

for ns in range(0, NSAMP) :

  if (ns % 1000) == 0 :
    print "ns = ",ns,"/",NSAMP

# LSB and USB signals in each nyquist interval (note that signal is random noise)
  sL = numpy.zeros( 4, dtype=complex)   
  sU = numpy.zeros( 4, dtype=complex)   
  for n in range(0, 4) :
    sL[n] = numpy.random.normal() + 1j * numpy.random.normal()
    sU[n] = numpy.random.normal() + 1j * numpy.random.normal()

# signal voltages from the receivers
  v1 = sL * w1 + sU * numpy.conj(w1)
  v2 = sL * w2 + sU * numpy.conj(w2)

# phase wrap for image sideband, v3
  phi = (ns/100.)*math.pi

# add random noise
  for n in range(0,4) :
    v1[n] = v1[n] + (numpy.random.normal() + 1j * numpy.random.normal())/snr
    v2[n] = v2[n] + (numpy.random.normal() + 1j * numpy.random.normal())/snr
    v3L[n] = sL[n] + sU[n] * numpy.exp(1j*phi) + (numpy.random.normal() + 1j * numpy.random.normal())/snr3
    v3U[n] = sL[n] * numpy.exp(-1j*phi) + sU[n] + (numpy.random.normal() + 1j * numpy.random.normal())/snr3

# LSB and USB beamformed sums
  vL = v1 * numpy.conj(w1) + v2 * numpy.conj(w2)
  vU = v1 * w1 + v2 * w2

# phase wrap for image sideband, v3
  phi = (ns/100.)*math.pi
  
# find accumulated cross-correlated power in each category
  xL = xL + sL * numpy.conj(sL)     # true LSB signal power
  xU = xU + sU * numpy.conj(sU)     # true USB signal power
  xLmeas = xLmeas + vL * numpy.conj(v3L)
  xUmeas = xUmeas + vU * numpy.conj(v3U)
  xL2meas = xL2meas + (vL+vU) * numpy.conj(v3L)
  xU2meas = xU2meas + (vL+vU) * numpy.conj(v3U)
   
# compute the means
xL = xL/NSAMP
xU = xU/NSAMP
xLmeas = xLmeas/NSAMP
xUmeas = xUmeas/NSAMP
xL2meas = xL2meas/NSAMP
xU2meas = xU2meas/NSAMP
 
print "xL      ", xL
print "xU      ", xU
print "xLmeas  ", xLmeas
print "xUmeas  ", xUmeas
print "xL2meas ", xL2meas
print "xU2meas ", xU2meas

print "xL      %.3f" % abs(numpy.mean(xL))
print "xU      %.3f" % abs(numpy.mean(xU))
print "xLmeas  %.3f" % abs(numpy.mean(xLmeas))
print "xUmeas  %.3f" % abs(numpy.mean(xUmeas))
print "xL2meas %.3f" % abs(numpy.mean(xL2meas))
print "xU2meas %.3f" % abs(numpy.mean(xU2meas))
