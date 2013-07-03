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
import matplotlib.pyplot as pyplot
from matplotlib.backends.backend_pdf import PdfPages

     
# Extract Stokes (I,Q,U, or V) amplitude and rms from one line of uvflux output 
def parseLine( line ) :
  #print line
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
  #print "\n%7.4f %5.4f  %6.4f %5.4f  %6.4f %5.4f  %6.4f %5.4f" % (I, rmsI, Q, rmsQ, U, rmsU, V, rmsV) 
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

# time average selected data from visFile into 'sstmp', with options=nopol
# this saves time in creating a multichannel SS object because input file
#    is re-read for every LkFile channel interval
def avgtoSStmp( visFile, selectStr ) :
  p= subprocess.Popen( ( shlex.split('rm -rf sstmp' ) ), \
     stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.STDOUT) 
  time.sleep(1)
  p= subprocess.Popen( ( shlex.split('uvaver vis=%s select=%s out=sstmp options=nopol' \
     % (visFile, selectStr) ) ), \
     stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.STDOUT) 
  result = p.communicate()[0]
  print result

# an SS ("StokeSpectrum") object is a collection of I,Q,U,V measurements and uncertainties vs freq
#   for a particular source over a particular time range
# the SS object may include a fit to the polarized flux, PA (at some ref freq) and rotation measure

class SS :
  
  def __init__(self) :
    self.visFile = None
    self.strList = []
    self.f1 = numpy.zeros( 1 )
    self.f2 = numpy.zeros( 1 )
    self.I = numpy.zeros( 1 )
    self.rmsI = numpy.zeros( 1 )
    self.Q = numpy.zeros( 1 )
    self.rmsQ = numpy.zeros( 1 )
    self.U = numpy.zeros( 1 )
    self.rmsU = numpy.zeros( 1 )
    self.V = numpy.zeros( 1 )
    self.rmsV = numpy.zeros( 1 )
    self.p0 = 0 
    self.p0rms = 0.
    self.pa0 = 0.
    self.pa0rms = 0.
    self.RM = 0. 
    self.RMrms = 0.
    self.freq0 = 0.

  def make( self, visFile, selectStr, LkFile) :
    'generate multichannel SS object from visibility data using leakages from LkFile'
    self.visFile = visFile        
    self.selectStr = selectStr
    self.LkFile = LkFile
    print "creating temporary file sstmp"
    avgtoSStmp( visFile, selectStr ) 
    vis = { "fileName" : "sstmp", "selectStr" : selectStr }
    frq1 = []
    frq2 = []
    Stokes = []									# effectively, a temporary work area
    self.strList = makeStrList( LkFile )	    # make list of available channel ranges in LkFile
    leakSolve.makeFtable( vis )					# fill in the chfreq and chwidth items in vis dictionary
    for lineStr in self.strList :
      [fstartLeak,fstopLeak] = leakSolve.restoreLk( LkFile, lineStr, "sstmp" )
      a = lineStr.split(",")   # "chan,1,49,8"
      ch1 = int(a[2])
      ch2 = ch1 + int(a[3]) - 1
      fstart = vis["chfreq"][ch1-1] - 0.5*vis["chwidth"][ch1-1]
      fstop = vis["chfreq"][ch2-1] + 0.5*vis["chwidth"][ch2-1]
      if (abs(fstartLeak - fstart) > .0005) or (abs(fstopLeak - fstop) > .0005) :
        print "    Leak %10.3f - %10.3f" % (fstartLeak, fstopLeak)
        print "    Data %10.3f - %10.3f" % (fstart, fstop)
      result = getStokes2( vis["fileName"], vis["selectStr"], lineStr )
      if numpy.isnan(result[0]) or numpy.isnan(result[2]) or numpy.isnan(result[4]) :
        print "skipping data with nan"
      else :
        frq1.append( fstart )
        frq2.append( fstop )
        Stokes.append( result )

    self.visFile = visFile               # change visfile name back to originating file
    self.f1 = numpy.array(frq1)
    self.f2 = numpy.array(frq2)
    self.I = numpy.array(Stokes).transpose()[0]
    self.rmsI = numpy.array(Stokes).transpose()[1]
    self.Q = numpy.array(Stokes).transpose()[2]
    self.rmsQ = numpy.array(Stokes).transpose()[3]
    self.U = numpy.array(Stokes).transpose()[4]
    self.rmsU = numpy.array(Stokes).transpose()[5]
    self.V = numpy.array(Stokes).transpose()[6]
    self.rmsV = numpy.array(Stokes).transpose()[7]
    self.UT, self.parang, self.HA = paPlot.getUTPAHA( self.visFile, self.selectStr )

  # writes to file; if outfile=None, prints to screen
  def dump( self, outfile ) :
    outStr = "# visFile     : %s\n" % self.visFile 
    outStr = outStr + "# LkFile      : %s\n" % self.LkFile 
    outStr = outStr + "# selectStr   : %s\n" % self.selectStr 
    outStr = outStr + "# avg UT      : %.3f hrs\n" % self.UT 
    outStr = outStr + "# avg HA      : %.3f hrs\n" % self.HA 
    outStr = outStr + "# avg parang  : %.3f deg\n" % self.parang 
    outStr = outStr + "# p0          : %.4f  ( %.4f ) Jy\n" % (self.p0, self.p0rms) 
    outStr = outStr + "# PAfit       : %.2f  ( %.2f ) deg\n" % (self.pa0,self.pa0rms)     
    outStr = outStr + "# RMfit       : %.4f  ( %.4f ) x 1.e5 rad/m^2\n" % (self.RM/1.e5, self.RMrms/1.e5) 
    outStr = outStr + "# PAfreq0     : %.3f GHz\n" % self.freq0 
    outStr = outStr + \
      "#  f1       f2           I        rmsI           Q        rmsQ           U        rmsU           V        rmsV\n"
    for f1,f2,I,rmsI,Q,rmsQ,U,rmsU,V,rmsV,selectStr in zip( self.f1, self.f2, self.I, self.rmsI, self.Q, \
        self.rmsQ, self.U, self.rmsU, self.V, self.rmsV, self.strList ) : 
      outStr = outStr + "%8.3f %8.3f   %10.3e %10.3e   %10.3e %10.3e   %10.3e %10.3e   %10.3e %10.3e   %s\n" \
         % (f1,f2,I,rmsI,Q,rmsQ,U,rmsU,V,rmsV,selectStr) 
    if outfile == None :
      print outStr
    else :
      fout = open( outfile, "a" )
      fout.write( "%s" % outStr )
      fout.write( "# ----------------------------\n" )       # this delimiter needed by readAll
      fout.close()

  # model Stokes Q, U as a function of 1/lambda^2
  # ... x = c^2/freq^2 - c^2/freq0^2, in m^2; each x given twice, once for Q, once for U
  # ... p0 is polarized intensity in Jy; a real number
  # ... pa0 is the position angle at the reference wavelength x=0, in radians; real number
  # ... rm is the rotation measure in radians/m^2; real number
  def func( self, x, p0, pa0, rm ) :
    xhalf = x[0 : len(x)/2]
    qmodel = p0 * numpy.cos(2 * (pa0 + rm*xhalf))
    umodel = p0 * numpy.sin(2 * (pa0 + rm*xhalf))
    return numpy.concatenate( (qmodel,umodel) )

  def fitPARM( self, freq0 ) :
    freq = 0.5 * (self.f1 + self.f2)
    lambdasq = (.30*.30)/(freq*freq) - (.30*.30)/(freq0*freq0)
    x = numpy.concatenate( (lambdasq,lambdasq) )
    y = numpy.concatenate( (self.Q,self.U) )
    sigma = numpy.concatenate( (self.rmsQ,self.rmsU) )
    Qavg = numpy.average( self.Q )
    Uavg = numpy.average( self.U )
    p0 = math.sqrt( Qavg*Qavg + Uavg*Uavg) 
    pa0 = 0.5 * math.atan(Uavg/Qavg)
    print "start with p0 = %.3f, pa0 = %.2f" % (p0, pa0 * 180./math.pi)
    popt,pcov = scipy.optimize.curve_fit( self.func, x, y, p0=[p0,pa0,0.], sigma=sigma  )
    if popt[0] < 0. :                 # if fit amplitude is negative, add 180 to 2 * PA
      popt[0] = -1. * popt[0]
      popt[1] = popt[1] + math.pi/2.
    p0,pa0,rm = popt
    pa0 = pa0 * 180./math.pi
    # compute chisq and scale covariances to get 'true' error estimates - see scipy thread
    #     by Christoph Deil; makes errors unreasonably small, so I am dropping it
    #  chi = y - self.func( x, p0, pa0, rm)/sigma
    #  chisq = (chi**2).sum()
    #  print "chisq = %.2e" % chisq
    #  dof = len(x) - len(popt)
    #  pcov = pcov/(chisq/dof)
    try :
      p0rms = math.sqrt(pcov[0][0])
      pa0rms = math.sqrt(pcov[1][1]) * 180./math.pi
      rmrms = math.sqrt(pcov[2][2]) 
    except :
      p0rms = pa0rms = rmrms = 0.
    print "p0 = %.3f (%.3f)" % (p0, p0rms)
    print "pa0 = %.2f (%.2f)" % (pa0, pa0rms)
    print "RM = %.3e (%.3e)" % (rm, rmrms)
    self.p0 = p0
    self.p0rms = p0rms
    self.pa0 = pa0
    self.pa0rms = pa0rms
    self.RM = rm
    self.RMrms = rmrms
    self.freq0 = freq0

  def polm( self ) :
    'compute fractional polarization and uncertainty'
    #Iavg = numpy.average(ss.I, weights=ss.rmsI )
    Iavg = numpy.average( self.I )
    Istd = numpy.std( self.I, ddof=1 )/math.sqrt(len(self.I))
    polm = self.p0/Iavg
    polmstd = polm * math.sqrt( pow(Istd/Iavg,2.) + pow(self.p0rms/self.p0,2.) )
    return [polm, polmstd]

  def plot( self ) :
    pyplot.clf()
    pyplot.ion()
    Ymax = abs(1.2 * numpy.concatenate( (self.Q, self.U, -1.*self.Q, -1*self.U) ).max())
    fig = pyplot.subplot(1,1,1)
    self.plot2( fig, Ymax, labelsize=12 ) 
    pyplot.draw()                 # plots and continues...

  def plot2( self, fig, Ymax, labelsize ) :
    fmin = numpy.concatenate( (self.f1,self.f2) ).min()
    fmax = numpy.concatenate( (self.f1,self.f2) ).max()
    fmin = fmin - 0.1 * (fmax - fmin)
    fmax = fmax + 0.1 * (fmax - fmin)
    if Ymax == 0. :
      Ymax = abs(1.3 * numpy.concatenate( (self.Q, self.U, -1.*self.Q, -1*self.U) ).max())
    freq = 0.5 * (self.f1 + self.f2)
    dfreq = 0.5 * (self.f1 - self.f2)
    fig.axis( [fmin, fmax, -1.*Ymax, Ymax], size=3 )
    fig.tick_params(axis='both', which='major', labelsize=labelsize )
    fig.grid( True )
    fig.errorbar( freq, self.Q, xerr=dfreq, yerr=self.rmsQ, fmt='ro', markersize=3 )
    fig.errorbar( freq, self.U, xerr=dfreq, yerr=self.rmsU, fmt='bo', markersize=3 )
    fx = numpy.arange(fmin,fmax,.1)
    xp = (.30*.30)/(fx*fx) - (.30*.30)/(self.freq0 * self.freq0)
    data = self.func( numpy.concatenate( (xp,xp) ), self.p0, math.pi*self.pa0/180., self.RM)
    Qfit = data[0:len(data)/2]
    Ufit = data[len(data)/2 :]
    fig.plot( fx, Qfit, 'r-' )
    fig.plot( fx, Ufit, 'b-' )
    data = self.func( numpy.concatenate( (xp,xp) ), self.p0, \
        math.pi*(self.pa0)/180., (self.RM - self.RMrms) )
    Qfit = data[0:len(data)/2]
    Ufit = data[len(data)/2 :]
    fig.plot( fx, Qfit, 'r--' )
    fig.plot( fx, Ufit, 'b--' )
    data = self.func( numpy.concatenate( (xp,xp) ), self.p0, \
        math.pi*(self.pa0)/180., (self.RM + self.RMrms) )
    Qfit = data[0:len(data)/2]
    Ufit = data[len(data)/2 :]
    fig.plot( fx, Qfit, 'r--' )
    fig.plot( fx, Ufit, 'b--' )
    polm, polmrms = self.polm()
    textStr = "pol = %.2f (%.2f)%%\nPA = %.1f (%.1f)\nRM = %.1e (%.1e)" % \
      ( 100.*polm, 100*polmrms, self.pa0 ,self.pa0rms, self.RM, self.RMrms)
    props = dict(boxstyle='round', facecolor='wheat', alpha=1 )	# patch.Patch properties
    fig.text( .5, .05, textStr, transform=fig.transAxes, horizontalalignment='center', 
      verticalalignment='bottom', fontsize=labelsize, bbox=props )
    pyplot.title( self.selectStr, size=labelsize, y=.9+.05*labelsize/5. )  

# generate table of SS objects
def makeTable( visFile, LkFile, srcName, extra, nint, maxgap, outfile, plot=False ) :
  selectList = paPlot.makeSelectList( visFile, srcName, nint, maxgap )
  nlist = len(selectList)
  for n, selectString in enumerate(selectList) :
    if len( extra ) > 0 :
      selectString = selectString + "," + extra
    print " "
    print "  %d/%d  %s" % (n+1,nlist,selectString)
    ss = SS()
    ss.make( visFile, selectString, LkFile )
    ss.fitPARM( 226. )
    if plot : ss.plot()
    ss.dump( outfile )
                                                                                                                  
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
   

# read in series of SSs from input file
def readAll(  infile, ssList ) :
  fin = open( infile, "r" )
  ss = SS()
  f1, f2, I, rmsI, Q, rmsQ, U, rmsU, V, rmsV = ( [] for i in range(10) )

  for line in fin :
    a = line.split()
    if line == "# ----------------------------\n" :     # this delimits end of each SS list
      ss.f1 = numpy.array(f1)
      ss.f2 = numpy.array(f2)
      ss.I = numpy.array(I)
      ss.rmsI = numpy.array(rmsI)
      ss.Q = numpy.array(Q)
      ss.rmsQ = numpy.array(rmsQ)
      ss.U = numpy.array(U)
      ss.rmsU = numpy.array(rmsU)
      ss.V = numpy.array(V)
      ss.rmsV = numpy.array(rmsV)
      ssList.append(ss)       # add completed SS object to the list
      ss = SS()               # initialize a new SS object
      f1, f2, I, rmsI, Q, rmsQ, U, rmsU, V, rmsV = ( [] for i in range(10) )
    elif line.startswith("# visFile") : ss.visFile = a[3]
    elif line.startswith("# LkFile") : ss.LkFile = a[3]
    elif line.startswith("# selectStr") : ss.selectStr = a[3]
    elif line.startswith("# avg UT") : ss.UT = float( a[4] )
    elif line.startswith("# avg HA") : ss.HA = float( a[4] )
    elif line.startswith("# avg parang") : ss.parang = float( a[4] )
    elif line.startswith("# p0") : 
      ss.p0 = float( a[3] )
      ss.p0rms = float( a[5] )
    elif line.startswith("# PAfit") : 
      ss.pa0 = float( a[3] )
      ss.pa0rms = float( a[5] )
    elif line.startswith("# RMfit") : 
      ss.RM = float( a[3] ) * 1.e5
      ss.RMrms = float( a[5] ) * 1.e5
    elif line.startswith("# PAfreq0") : ss.freq0 = float( a[3] )
    elif not line.startswith("#") :
      f1.append( float( a[0] ) )
      f2.append( float( a[1] ) )
      I.append( float(a[2]) )
      rmsI.append( float( a[3] ) )
      Q.append( float( a[4] ) )
      rmsQ.append( float( a[5] ) )
      U.append( float( a[6] ) )
      rmsU.append( float( a[7] ) )
      V.append( float( a[8] ) )
      rmsV.append( float( a[9] ) )
      ss.strList.append( a[10] ) 

# read in one or more PARM files, write summary file for wip
def plotPA( paList, outfile ) :
  srcList = [ "0234+285", "3c84", "3c111", "0721+713", "0854+201", "OJ287", "3c273", "M87", \
	"3c279", "1337-129", "1633+382", "3c345", "NRAO530", "SgrA", "1749+096", "1921-293", \
    "bllac", "titan" ]
  colorList = [ 3,2,14,4,3,3,4,8,1,9,10,11,12,13,5,15,1,2]
  nstart = 1
  fin = open( paList, "r" )
  for line in fin :
    a = line.split()
    if not line.startswith("#") and len(a) > 0 :
      ssList = []
      infile = a[0]
      for src,col in zip (srcList,colorList) :
        if (a[0].find(src) >= 0) or (a[0].find(src.capitalize()) >= 0) :
          srcName = src
          color = col
      fout = open( outfile, "a" )
      fout.write("#\n")
      fout.write("# %s\n" % infile)
      fout.write("#  dechr   parang     HA        S    sigma     poli  sigma     PA  sigma      RM   sigma  col   selectString\n")
      fout.close()
      readAll( infile, ssList )
      for ss in ssList : 
        #Iavg = numpy.average(ss.I, weights=ss.rmsI )
        Iavg = numpy.average(ss.I )
        Istd = numpy.std(ss.I, ddof=1)
        fout = open(outfile, 'a')
        fout.write("%8.3f %8.2f %8.3f %8.3f %6.3f %8.3f %6.3f %7.1f %5.1f %8.3f %6.3f  %2d   %s\n" % \
          (ss.UT, ss.parang, ss.HA, Iavg, Istd, ss.p0, ss.p0rms, ss.pa0, ss.pa0rms, ss.RM/1.e5, ss.RMrms/1.e5, color, ss.selectStr) )
        fout.close()
      nstop = nstart + len(ssList) + 2
      print "1src %d %d %d %s" % (nstart, nstop, color, srcName )
      nstart = nstop + 1
  fin.close()
        
def replot( paList, Ymax=0., nrows=2, ncols=1 ) :
  pyplot.ioff()
  pyplot.clf()
  pp = PdfPages( 'multipage.pdf' )
  fin = open( paList, "r" )
  nplot = 1
  for line in fin :
    if not line.startswith("#") :
      ssList = []
      a = line.split()
      paFile = a[0]
      readAll( paFile, ssList )
      pyplot.subplots_adjust( hspace=0.25 )
      for ss in ssList :
        print "plotting file %s, string %s" % (paFile, ss.selectStr)
        if nplot > nrows*ncols : 
           pyplot.savefig( pp, format='pdf' )
           nplot = 1
           pyplot.clf()
        fig = pyplot.subplot( nrows, ncols, nplot )
        ss.plot2( fig, Ymax, labelsize=10./ncols )
        nplot = nplot + 1
  pyplot.savefig( pp, format='pdf' )
  pp.close()
  
def test() :
  ssList = []
  readAll("3c84.c3.PA", ssList)
  for ss in ssList :
    polm,polmstd = ss.polm()
    print polm,polmstd 
