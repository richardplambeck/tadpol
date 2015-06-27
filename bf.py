# bf.py
# is it possible to phase LSB and USB,then add the 2 voltages together, without losing s/n ?

import math
import numpy

snr = 10.     # snr for beamformed antennas
snrref = 10.   # snr for refant

NSAMP = 10
phref = 2.*math.pi*numpy.random.random()
phref = 0.

# I tried to model just 2 antennas, but ran into problems simulating the 
# case where vL and vU are summed before going to the beamformer: for one
# piece of the Walsh sequence, vL+vU = 0! no power sent to recorder! so
# simulate 8 antennas with typical Walsh sequence 

# LSB and USB signal voltages in each Walsh interval
sL = numpy.zeros( 8, dtype=complex)   
sU = numpy.zeros( 8, dtype=complex)   

# beamformed sums for each Walsh interval
vL = numpy.zeros( 8, dtype=complex)   
vU = numpy.zeros( 8, dtype=complex)   

# refernece antenna voltage for each Walsh interval
vrefL = numpy.zeros( 8, dtype=complex)   
vrefU = numpy.zeros( 8, dtype=complex)   

# unnormalized amplitude of signal
xL = numpy.zeros( 8, dtype=complex)
xU = numpy.zeros( 8, dtype=complex)

# normalized measured cross correaltions between beamformer and refant
xLmeas = numpy.zeros( 8, dtype=complex)
xUmeas = numpy.zeros( 8, dtype=complex)

# normalized measured cross correaltions between (vL+vU) and [vrefL, vrefU]
xSLmeas = numpy.zeros( 8, dtype=complex)
xSUmeas = numpy.zeros( 8, dtype=complex)

# Walsh LO1 phase switch sequence for each of these 4 nyquist intervals
#w1 = numpy.array( [ 1, 0+1j, 1,    0+1j ] )
#w2 = numpy.array( [ 1, 1   , 0+1j, 0+1j ] )
#w1 = numpy.array( [1,1,1,1], dtype=complex )
#w2 = numpy.array( [1,1,1,1], dtype=complex )

# normalized cross correaltions between each dewalshed antenna and refant
xLant = numpy.zeros( (8,8), dtype=complex )
xUant = numpy.zeros( (8,8), dtype=complex )

# antenna voltges including noise in ech walsh interval
v = numpy.zeros( (8,8), dtype=complex )

w = numpy.zeros( (8,8), dtype=complex )
w[0] = numpy.array( [  1,  1,  1,  1,  1,  1,  1,  1 ] )   # ref antenna
w[1] = numpy.array( [  1,  1,  1,  1, 1j, 1j, 1j, 1j ] )
w[2] = numpy.array( [  1,  1, 1j, 1j, 1j, 1j,  1,  1 ] )
w[3] = numpy.array( [ 1j, 1j,  1,  1, 1j, 1j,  1,  1 ] )
w[4] = numpy.array( [ 1j,  1,  1, 1j, 1j,  1,  1, 1j ] )
w[5] = numpy.array( [  1, 1j, 1j,  1, 1j,  1,  1, 1j ] )
w[6] = numpy.array( [  1, 1j,  1, 1j, 1j,  1, 1j,  1 ] )
w[7] = numpy.array( [ 1j,  1, 1j,  1, 1j,  1, 1j,  1 ] )

for ns in range(0, NSAMP) :

  if (ns % 1000) == 0 :
    print "ns = ",ns,"/",NSAMP

# LSB and USB signals in each nyquist interval (note that signal is random noise)
  for n in range(0, 8) :
    sL[n] = numpy.random.normal() + 1j * numpy.random.normal()
    sU[n] = numpy.random.normal() + 1j * numpy.random.normal()
    sL[n] = 1+0j
    sU[n] = 0.001+0j

# signal voltages from the receivers
  for ant in range(1,8) :
    v[ant] = sL[n] * w[ant] + sU * numpy.conj(w[ant])

# phase wrap for image sideband, vref; fixed for each Walsh interval
  phi = (ns/100.)*math.pi

# add random noise
  for n in range(0,8) :
    for ant in range(1,8) :
      v[ant][n] = v[ant][n] #+ (numpy.random.normal() + 1j * numpy.random.normal())/snr
    vrefL[n] = sL[n] #*numpy.exp(1j*phref) + 0.*sU[n]*numpy.exp(1j*phi) #+ (numpy.random.normal() + 1j * numpy.random.normal())/snrref
    vrefU[n] = sU[n] #0.*sL[n]*numpy.exp(-1j*phi) + sU[n]*numpy.exp(-1j*phref) #+ (numpy.random.normal() + 1j * numpy.random.normal())/snrref

# LSB and USB beamformed sums
  vL = numpy.zeros( 8, dtype=complex)   
  vU = numpy.zeros( 8, dtype=complex)   
  for ant in range(1,8) :
    vL = vL + v[ant] * numpy.conj(w[ant])  
    vU = vU + v[ant] * w[ant]
  vS = vL + vU
  tpL = abs(vL)
  tpU = abs(vU)
  tpS = abs(vS)

# find accumulated cross-correlated power in each category
  for ant in range(1,8) :
    xLant[ant] = xLant[ant] + v[ant] * numpy.conj( w[ant] ) * numpy.conj(vrefL) / (abs(v[ant])*abs(vrefL)) 
    xUant[ant] = xUant[ant] + v[ant] * w[ant] * numpy.conj(vrefU) / (abs(v[ant])*abs(vrefU) )
  xL = xL + sL * numpy.conj(sL)     # true LSB signal power
  xU = xU + sU * numpy.conj(sU)     # true USB signal power
  xLmeas = xLmeas + vL * numpy.conj(vrefL) / (tpL * abs(vrefL))
  xUmeas = xUmeas + vU * numpy.conj(vrefU) / (tpU * abs(vrefU))
  print "xLmeas: ", xLmeas
  print "xUmeas: ", xUmeas
  for n in range(0,8):
    if tpS[n] == 0. :
      tpS[n] = 1.
  xSLmeas = xSLmeas + (vL+vU) * numpy.conj(vrefL) / (tpS * abs(vrefL))
  xSUmeas = xSUmeas + (vL+vU) * numpy.conj(vrefU) / (tpS * abs(vrefU))
   
print "vL :", vL
print "vU :", vU
# compute the means
xLantx = numpy.delete( xLant, 0, 0 )  # first row is all zeros (there is no zero antenna)
xUantx = numpy.delete( xUant, 0, 0 )
xLantx = xLantx/NSAMP
xUantx = xUantx/NSAMP
xL = xL/NSAMP
xU = xU/NSAMP
xLmeas = xLmeas/NSAMP
xUmeas = xUmeas/NSAMP
xSUmeas = xSUmeas/NSAMP
xSLmeas = xSLmeas/NSAMP
 
#print "xL      ", xL
#print "xU      ", xU
print "xLant   ", xLantx
print "xUant   ", xUantx
print "xLmeas  ", xLmeas
print "xUmeas  ", xUmeas
print "xSLmeas ", xSLmeas
print "xSUmeas ", xSUmeas

#print "xL      %.3f" % abs(numpy.mean(xL))
#print "xU      %.3f" % abs(numpy.mean(xU))
print "xLant   %.3f" % abs(numpy.mean(xLantx))
print "xUant   %.3f" % abs(numpy.mean(xUantx))
print "xLmeas  %.3f" % abs(numpy.mean(xLmeas))
print "xUmeas  %.3f" % abs(numpy.mean(xUmeas))
print "xSLmeas %.3f" % abs(numpy.mean(xSLmeas))
print "xSUmeas %.3f" % abs(numpy.mean(xSUmeas))
