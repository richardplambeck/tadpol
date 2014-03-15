# bf.py
# is it possible to phase LSB and USB,then add the 2 voltages together, without losing s/n ?

import math
import numpy

snr12 = 1000.  # snr for 1 and 2; 
snr3 = 1000   # snr for 3

NSAMP = 100000

n1 = numpy.zeros( NSAMP, dtype=complex )
v1 = numpy.zeros( NSAMP, dtype=complex)
n2 = numpy.zeros( NSAMP, dtype=complex)
v2 = numpy.zeros( NSAMP, dtype=complex)
vL = numpy.zeros( 4, dtype=complex)
vU = numpy.zeros( 4, dtype=complex)

# Walsh sequences for antennas 1 and 2
w1 = numpy.array( [1., 0+1j, 1., 0+1j] )
w2 = numpy.array( [1., 1., 0+1j, 0+1j] )
p1 = numpy.array( [1., 0-1j, 1., 0-1j] )
p2 = numpy.array( [1., 1., 0-1j, 0-1j] )

# pick random signal voltage in LSB and USB
sL = numpy.random.normal() + 1j * numpy.random.normal()
sU = numpy.random.normal() + 1j * numpy.random.normal()

# accumulators
vsl  =  0 + 0j 
vsl2 =  0 + 0j
vsu  =  0 + 0j
vsu2 =  0 + 0j
 
# what signals come from the antennas in each interval
for ns in range(0,NSAMP):
  n = ns % 4.
  if (ns % 10000) == 0 :
    print "ns = ",ns," / ",NSAMP
  v1 =(numpy.random.normal() + 1j * numpy.random.normal())/snr12 + w1[n]*sL + numpy.conj(w1[n])*sU
  v2 =(numpy.random.normal() + 1j * numpy.random.normal())/snr12 + w2[n]*sL + numpy.conj(w2[n])*sU
  vL = v1 * numpy.conj(w1[n]) + v2 * numpy.conj(w2[n])   # LSB phased sum
  vU = v1 * w1[n] + v2 * w2[n]                           # USB phased sum
 
  n3 = (numpy.random.normal() + 1j * numpy.random.normal())/snr3
  v3L = sL + numpy.exp(1j*(ns%200)/100.*math.pi)*sU + n3   # reference antenna LSB corr (USB phase wraps)
  v3U = sU + numpy.exp(1j*(ns%200)/100.*math.pi)*sL + n3   # reference antenna USB corr (LSB phase wraps)

  vsl = vsl + vL * numpy.conj(v3L)
  vsl2 = vsl2 + (vL+vU) * numpy.conj(v3L)
  vsu = vsu + vU * numpy.conj(v3U)
  vsu2 = vsu2 + (vL+vU) * numpy.conj(v3U)

vsu = vsu/(2.*NSAMP) 
vsl = vsl/(2.*NSAMP) 
vsu2 = vsu2/(2.*NSAMP) 
vsl2 = vsl2/(2.*NSAMP) 

print sL*numpy.conj(sL), vsl, vsl2, abs(sL*numpy.conj(sL) - vsl), abs(sL*numpy.conj(sL) - vsl2)
print sU*numpy.conj(sU), vsu, vsu2, abs(sU*numpy.conj(sU) - vsu), abs(sU*numpy.conj(sU) - vsu2)

