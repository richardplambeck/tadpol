# RM.py
#
# bits from paPlot.py, etc

import math
import time
import numpy
import subprocess
import shlex
import string
import leakSolve
import paPlot
import scipy.optimize
import pylab

     
# Extract Stokes (I,Q,U, or V) amplitude and rms from one line of uvflux output 
def parseLine( line ) :
  print line
  a = line.split()
  if (len(a) == 9) :
    src = a[0]
    stokes = float(a[3])
    ncorrs = int(a[8])
    rms = float(a[5])/math.sqrt(ncorrs)
  else :
    src = ' '
    stokes = float(a[2])
    ncorrs = int(a[7])
    rms = float(a[4])/math.sqrt(ncorrs)
  return [ src, stokes, rms, ncorrs ]

# Run uvflux on selected data; return I, Q, U, V, and rms errors
def getStokes2( infile, selectString, lineString ) :
  p= subprocess.Popen( ( shlex.split('uvflux vis=%s select=%s line=%s stokes=I,Q,U,V' \
     % (infile, selectString, lineString) ) ), \
     stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.STDOUT) 
  result = p.communicate()[0]
  lines = result.split("\n")
  I = rmsI = pa = sigmaPa = poli = sigmaPoli = V = rmsV = 0.	# in case we get error
  Q = rmsQ = U = rmsU = V = rmsV = 0.
  nlines = len(lines)
  if nlines >= 10 :
    [src, I, rmsI, ncorrs] = parseLine( lines[nlines-6] )
    [dummy, Q, rmsQ, ncorrs] = parseLine( lines[nlines-5] )
    [dummy, U, rmsU, ncorrs] = parseLine( lines[nlines-4] )
    [dummy, V, rmsV, ncorrs] = parseLine( lines[nlines-3] )
  print "\n%7.4f %5.4f  %6.4f %5.4f  %6.4f %5.4f  %6.4f %5.4f" % (I, rmsI, Q, rmsQ, U, rmsU, V, rmsV) 
  return [ I, rmsI, Q, rmsQ, U, rmsU, V, rmsV ]

# makes list of lineStrings in LkFile
def makeStrList( LkFile ) :
  fin = open( LkFile, "r" )
  strList = []
  for line in fin :
    if not line.startswith("#") and len(line) > 1 :
      a = line.split()
      if a[0] == "C01" :
        strList.append( a[9] )
  fin.close()
  return strList


# Generate table of I,Q,U,V vs channel range in file vis["fileName"], using leakages from file LkFile
def makeIQUspectrum( vis, outfile, LkFile ) :
  frqlist = []
  Stokes = []
  strList = makeStrList( LkFile )				# make list of available channel ranges in LkFile
  leakSolve.makeFtable( vis )					# fill in the chfreq and chwidth items in vis dictionary

  fout = open( outfile, "a" )
  fout.write("# Lkfile: %s   VisFile: %s,  selectStr: %s\n" % (LkFile, vis["fileName"], vis["selectStr"]) )
  fout.write("#  f1       f2           I        rmsI           Q        rmsQ           U        rmsU           V        rmsV\n" )
  fout.close()

  for lineStr in strList :
    [fstartLeak,fstopLeak] = leakSolve.restoreLk( LkFile, lineStr, vis["fileName"] )
    a = lineStr.split(",")   # "chan,1,49,8"
    ch1 = int(a[2])
    ch2 = ch1 + int(a[3]) - 1
    fstart = vis["chfreq"][ch1-1] - 0.5*vis["chwidth"][ch1-1]
    fstop = vis["chfreq"][ch2-1] + 0.5*vis["chwidth"][ch2-1]
    print "  Leak %7.3f - %7.3f" % (fstartLeak, fstopLeak)
    print "  Data %7.3f - %7.3f" % (fstart, fstop)
    [ I, rmsI, Q, rmsQ, U, rmsU, V, rmsV ] = getStokes2( vis["fileName"], vis["selectStr"], lineStr )
    fout = open( outfile, "a" )
    fout.write("%8.3f %8.3f   %10.3e %10.3e   %10.3e %10.3e   %10.3e %10.3e   %10.3e %10.3e\n" \
         % (fstart,fstop,I,rmsI,Q,rmsQ,U,rmsU,V,rmsV) )
    fout.close()
    frqlist.append( 0.5*(fstart + fstop) )
    Stokes.append( [I, rmsI, Q, rmsQ, U, rmsU] )

  frq = numpy.array(frqlist)
  freq0 = 226.
  Ia = numpy.array(Stokes).transpose()[0]
  Iavg = numpy.average( Ia )      # could be weighted
  Qa = numpy.array(Stokes).transpose()[2]
  rmsQa = numpy.array(Stokes).transpose()[3]
  Ua = numpy.array(Stokes).transpose()[4]
  rmsUa = numpy.array(Stokes).transpose()[5]
  [p0, p0rms, pa0, pa0rms, rm, rmrms] = fitPARM( frq, freq0, Qa, Ua, rmsQa, rmsUa )
  fout = open( outfile, "a" )
  fout.write("# results: poli = %.3f (%.3f), PA = %.2f (%.2f),  RM = %.3e (%.3e)\n" \
       % (p0, p0rms, pa0, pa0rms, rm, rmrms) )
  fout.close()
  plotfit( frq, freq0, Qa, Ua, p0, math.pi*pa0/180., rm )
  
# model Stokes Q, U as a function of 1/lambda^2
# ... x = c^2/freq^2 - c^2/freq0^2, in m^2; each x given twice, once for Q, once for U
# ... p0 is polarized intensity in Jy; a real number
# ... pa0 is the position angle at the reference wavelength x=0, in radians; real number
# ... rm is the rotation measure in radians/m^2; real number

def func( x, p0, pa0, rm ) :
  xhalf = x[0 : len(x)/2]
  qmodel = p0 * numpy.cos(2 * (pa0 + rm*xhalf))
  umodel = p0 * numpy.sin(2 * (pa0 + rm*xhalf))
  return numpy.concatenate( (qmodel,umodel) )
  

# least squares fit of Q,U vs freq to get PA and RM
def fitPARM( freq, freq0, Q, U, rmsQ, rmsU ) :
  xhalf = (.30*.30)/(freq*freq) - (.30*.30)/(freq0*freq0)   # lambdasq - lambdasq0, m^2
  x = numpy.concatenate( (xhalf,xhalf) )
  y = numpy.concatenate( (Q,U) )
  sigma = numpy.concatenate( (rmsQ,rmsU) )
  popt,pcov = scipy.optimize.curve_fit( func, x, y, p0=[0.,0.,0.], sigma=sigma  )
  if popt[0] < 0. :                 # if fit amplitude is negative, add 180 to 2 * PA
    popt[0] = -1. * popt[0]
    popt[1] = popt[1] + math.pi/2.
  p0,pa0,rm = popt
  pa0 = pa0 * 180./math.pi
  try :
    p0rms = math.sqrt(pcov[0][0])
    pa0rms = math.sqrt(pcov[1][1]) * 180./math.pi
    rmrms = math.sqrt(pcov[2][2]) 
  except :
    p0rms = pa0rms = rmrms = 0.
  print "p0 = %.3f (%.3f)" % (p0, p0rms)
  print "pa0 = %.2f (%.2f)" % (pa0, pa0rms)
  print "RM = %.3e (%.3e)" % (rm, rmrms)
  return [p0, p0rms, pa0, pa0rms, rm, rmrms]

def plotfit( freqlist, freq0, Q, U, p0, pa0, rm ) :
  pylab.clf()
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

  
def RMsearch( infile, outfile ) :
  RMmin = -5.e9
  RMmax = 5.e9
  RMstep = 5.e5

# read in the data
  plist = []
  flist = []
  fin = open( infile, "r" )
  for line in fin :
    if not line.startswith("#") :
      a = line.split()
      flist.append( (float(a[0]) + float(a[1]))/2. )
      Q = float(a[4])
      U = float(a[6])
      plist.append( Q + 1j*U )
  fin.close()  

  fGHz = numpy.array(flist)
  lambdasq = (.2998*.2998)/(fGHz*fGHz)
  pol = numpy.array(plist)
  norm = 1./len(pol)

  fout = open( outfile, "w" ) 
  for RM in range( RMmin, RMmax+RMstep, RMstep ) :
    rot = numpy.exp( 1j * RM * lambdasq )
    #print numpy.angle(rot, deg=True)
    P = norm * numpy.abs(numpy.sum( pol * rot )) 
    print "RM = %.3e, P = %5.3f" % (RM,P)
    fout.write( "%.3e  %11.4e\n" % ( RM, P ) )
  fout.close()
   
def test() :
  frq =  numpy.array([ 220.57, 219.07, 223.82, 223.27, 232.08, 233.57, 228.82, 229.369] )
  freq0 = 226.
  Qa = numpy.array( [.0953, .0935, .106, .0938, .0869, .08834, .07027, .09538] )
  Ua = numpy.array( [-.1177, -.1342, -.1321, -.1308,-.1305,-.1489,-.1456,-.1379] )
  rmsQa = rmsUa = numpy.array( [0.,0.,0.,0.,0.,0.,0.,0.] )
  [p0, p0rms, pa0, pa0rms, rm, rmrms] = fitPARM( frq, freq0, Qa, Ua, rmsQa, rmsUa )
  plotfit( frq, freq0, Qa, Ua, p0, math.pi*pa0/180., rm )

