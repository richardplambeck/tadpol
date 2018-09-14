# corrHorn.py
# find effective focus and FWHM of radiation pattern for a corrugated feed horn

import math

clight = 2.998e10   # cm/sec
freq = 90      # GHz

# standard WR10 CARMA horn; aperture .6986", flare angle 12 deg
a = .6986/2.        # radius of aperture, inches
flareAngle = 12     # half-angle of flare, degrees

# smaller corrugated horn
a = .5872/2.
flareAngle = 15

Rcap = a/math.sin(math.radians(flareAngle))    # inches
hornLength = Rcap * math.cos( math.radians(flareAngle) )

# Wylde eqn 23
k = 2.*math.pi/(clight/(freq*1.e9))    # cm^-1
M = k * pow(2.54*a,2) / (2.*2.54*Rcap)

T = 1/(1. + 6.769*pow(M/(2.*math.pi),2) )   # (phase center - apex distance)/(horn length)
  # Wylde eqn 45

w0 = math.sqrt( pow(.6435*a,2)/(1 + pow(0.6435*0.6435*M,2) ) )
  # Wylde eqn 27

print "hornLength = %.3f inches" % hornLength 
print "aperture diameter = %.3f inches" % (2.*a)
print "flare angle = %.1f degrees" % flareAngle
print "focus at %.0f GHz is %.3f inches behind aperture" % (freq,(1-T)*hornLength)
print "beam waist at focus is %.3f inches" % w0
print "far field FWHM = %.1f degrees" % (math.degrees( 1.18 * clight/(freq*1.e9) / (math.pi * 2.54 * w0) ) )


