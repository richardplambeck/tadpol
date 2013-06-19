import scipy.optimize
import numpy
import math
import random
import pylab

PA0in = math.pi * 70./180.
P0in = 3.5*.05
RMin = 5.e5
QUnoise = 0.02

freqlist = numpy.array([210., 211., 212, 213., 227., 228., 229., 230.])

# Create the data
x = (.30*.30)/(freqlist*freqlist) - (.30*.30)/(220.*220.) 
  # lambdasq array, units m^2
PAin = PA0in + RMin *  x 
   # noiseless PA array
#print 180./math.pi * PAin

if QUnoise > 1.e-7 :
  qr = numpy.random.normal(loc=0.,scale=QUnoise,size=len(freqlist))
  ur = numpy.random.normal(loc=0.,scale=QUnoise,size=len(freqlist))
else :
  qr = numpy.zeros( len(freqlist) )
  ur = numpy.zeros( len(freqlist) )

Q = P0in * numpy.cos(2.*PAin) + qr
U = P0in * numpy.sin(2.*PAin) + ur
Pol = Q + 1j*U

ydata = numpy.concatenate( (Q,U) )
xdata = numpy.concatenate( (x,x) )

def func( xdata, p0, pa0, rm ) :
  xloc = xdata[0: len(xdata)/2.]
  qmodel = p0 * numpy.cos(2 * (pa0 + rm*xloc))
  umodel = p0 * numpy.sin(2 * (pa0 + rm*xloc))
  return numpy.concatenate( (qmodel,umodel) )

print ""
print "model data:"
for freq,pol1 in zip( freqlist, Pol ) :
  print freq, pol1, numpy.abs(pol1), 0.5*numpy.angle( pol1, deg=True) 

# model complex polarization array P as a function of 1/lambda^2
# ... x = c^2/freq^2 - c^2/freq0^2, in m^2; an array
# ... p0 is polarized intensity in Jy; a real number
# ... pa0 is the position angle at the reference wavelength x=0, in radians; real number
# ... rm is the rotation measure in radians/m^2; real number


popt,pcov = scipy.optimize.curve_fit( func, xdata, ydata, p0=[0.,0.,0.], sigma=None )
if popt[0] < 0. :
  popt[0] = -1. * popt[0]
  popt[1] = popt[1] + math.pi/2.
print "p0 = %.3f (%.3f)" % (popt[0], math.sqrt(pcov[0][0]))
print "pa0 = %.2f (%.2f)" % (popt[1] * 180./math.pi, math.sqrt(pcov[1][1])*180./math.pi)
print "RM = %.3e (%.3e)" % (popt[2], math.sqrt(pcov[2][2]))


def plotfit( freqlist, freq0, Q, U, p0, pa0, rm ) :
  pylab.ion()
  fig = pylab.subplot(1,1,1)
  frange = freqlist.max() - freqlist.min()
  fmin = freqlist.min() - .05 * frange
  fmax = freqlist.max() + .05 * frange
  Ymax = abs(1.05 * numpy.concatenate( (Q,U,-1.*Q,-1*U) ).max())
  fig.axis( [fmin, fmax, -1.*Ymax, Ymax] )
  fig.grid( True )
  fig.plot( freqlist, Q, 'ro')
  fig.plot( freqlist, U, 'bo')
  fx = numpy.arange(fmin,fmax,.1)
  xp = (.30*.30)/(fx*fx) - (.30*.30)/(freq0*freq0) 
  data = func( numpy.concatenate( (xp,xp) ), p0, pa0, rm)
  Qfit = data[0:len(data)/2]
  Ufit = data[len(data)/2 :] 
  fig.plot( fx, Qfit, 'r-' )
  fig.plot( fx, Ufit, 'b-' )
  pylab.show()

plotfit( freqlist, 220., Q, U, P0in, PA0in, RMin)
