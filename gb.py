# gb.py - gaussian beam related calculations for PolarBear
# work out formulas 14a and 14b from Goldsmith

import numpy
import math
import sys

f = 3    # focal length of 3mm dewar window in inches
freq = 90.
lmbda = 300./freq   # wavelength in mm

# compute output waist diameter and distance as a function
#  of distance of lens from feed horn aperture
# calculation is all in inches

def lens( freq, d1, w01=0.2 ) :
    lmbda = (29.98/freq)/2.54    # wavelength in inches
    factor = pow(d1/f - 1., 2.) + pow(math.pi*w01*w01/(lmbda*f),2.) 

  # compute location and size of next beamwaist
    w02 = w01 / math.sqrt( factor )
    d2 = f * (1. + (d1/f - 1.)/factor)
    print "d2 = %.3f, w02 = %.3f" % (d2,w02)

  # now dump out beamwaist as a function of distance along the wavefront
    fout = open("beamprofile.dat", "w")
    fout.write("# flens = %.3f in, freq = %.2f GHz\n  d1 = %.3f" % (f,freq,d1) )
    for z in numpy.arange(0.,d1, 0.1) :
      w = w01 * math.sqrt(1. + pow( lmbda*z/(math.pi*w01*w01), 2.) )
        # Goldsmith eq 5
      fout.write("%6.3f %6.3f %5.3f\n" % ( z-d1, w, -w) )
      print z,w
    for z in numpy.arange(0., 2.*d2 + 0.1, 0.1) :
      w = w02 * math.sqrt(1. + pow( lmbda*(z-d2)/(math.pi*w02*w02), 2.) )
      fout.write("%6.3f %6.3f %5.3f\n" % ( z, w, -w) )
      print z,w
    fout.close()
    print "closing file"

# compute z(r) for 3mm lens, based on coefficients in bima memo
# surf = 0 for horn lens, 1 or 2 for dewar window lens
# R is effective radius of curvature - see PolarBear notebook, p. 21
def curv( nsurf, r ) :
  print "computing z and R for surface %d" % nsurf
  c = [1.21386, .43786, .39337]
  s = [-.29075, .16395, -3.13973]
  a4 = [.14593, 0., 0.]
  a6 = [.00483, .00004, -.00103]
  a8 = [.00819, .00000, .00041]
  a10 = [.0, .0, -.00007]
  z = c[nsurf]*pow(r,2)/(1. + math.sqrt(1. - s[nsurf]*pow(c[nsurf]*r,2)))
  z = z + a4[nsurf] * pow(r,4)
  z = z + a6[nsurf] * pow(r,6)
  z = z + a8[nsurf] * pow(r,8)
  z = z + a10[nsurf] * pow(r,10)
  R = (z*z + r*r)/(2.*z)
  print "z = %.3f, R = %.3f, f = %.3f" % (z,R,R/.442)

def plotbeam( infile, outfile ) :
  fin = open( infile, "r" )
  fout = open( outfile, "w" )
  x = 0.
  for line in fin :
    a = line.split()
    if a[2] == "space" :
      x = x + float(a[3])
    elif a[0] == "theta=" :
      y = float( a[5][0:-1] )
      print x, y, -y
      fout.write("%8.2f %10.7f %10.7f\n" % (x, y, -y) )
      x = 0.
  fin.close()
  fout.close()
