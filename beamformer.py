# beamformer.py
# compute s/n ratio for various weightings

import math
import numpy

wgt = numpy.ones( 8, dtype=float )
Tsys = 500.*numpy.ones( 8, dtype=float)
JperKnom = [ 65., 65., 65., 65., 65., 145., 145., 145. ]
g = [ 1.1, 1.1, 1.1, 1.1, 1.1, 1., 1., 1. ]
JperK = numpy.array(JperK) * pow( numpy.array(g), 2. )

def VSNR( wgt=wgt, JperK=JperK, Tsys=Tsys ) :
  vsig = 0.
  pnoise = 0.
  for n in range(0,8) :
    vsig = vsig + wgt[n] / math.sqrt(JperK[n])
    pnoise = pnoise + pow( wgt[n], 2.) * Tsys[n]
  vsnr = vsig / math.sqrt(pnoise)
  vsnr0 = 1. / math.sqrt(JperK[0] * Tsys[0] )
  print "wgt   :",wgt
  print "Tsys  :", Tsys
  print "JperK :", JperK
  print "vsnr = %.3e, Psnr/Psnr0 = %.3f" % (vsnr, pow(vsnr/vsnr0,2.))

print "\noptimum weighting: "
wgtopt = 1./ (numpy.sqrt(JperK) * Tsys)
VSNR(wgt=wgtopt )

print "\ntpwr weighting: "
wgtpwr = 1./ numpy.sqrt(Tsys)
VSNR(wgt=wgtpwr)
