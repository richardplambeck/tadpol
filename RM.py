# RM.py
#
# bits from paPlot.py, etc

import math
import sys
import time
import numpy
import subprocess
import shlex
import string
import leakSolve
import paPlot
import pickle
import copy
import scipy.optimize
import matplotlib.pyplot as pyplot
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.mplot3d import Axes3D


     
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
  print result
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
#def avgtoSStmp( visFile, selectStr ) :
#  'time-average from visFile to sstmp over interval specified by selectStr'
#  p= subprocess.Popen( ( shlex.split('rm -rf sstmp' ) ), \
#     stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.STDOUT) 
#  time.sleep(1)
#  p= subprocess.Popen( ( shlex.split('uvaver vis=%s select=%s out=sstmp options=nopol interval=10000' \
#     % (visFile, selectStr) ) ), \
#     stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.STDOUT) 
#  result = p.communicate()[0]
#  print result

# an SS ("StokeSpectrum") object is a collection of I,Q,U,V measurements and their uncertainties vs freq
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
    self.p0 = 0. 
    self.p0rms = 0.
    self.frac = 0. 
    self.fracrms = 0.
    self.pa0 = 0.
    self.pa0rms = 0.
    self.RM = 0. 
    self.RMrms = 0.
    self.freq0 = 0.

  def make( self, visFile, selectStr, LkFile, preavg="copy", f0=226. ) :
    'generate multichannel SS object from visibility data using leakages from LkFile'
    self.visFile = visFile        
    self.selectStr = selectStr
    self.LkFile = LkFile
    if preavg == "avg" :
      print "Preaveraging to create temporary file sstmp" 
      avgtoSStmp( visFile, selectStr ) 
        # preaveraging is not recommended because it leaves only 1 parallactic angle!
    elif preavg == "copy" :
      print "Copying (but not averaging) selected data to temporary file sstmp"
      leakSolve.copytoSStmp( visFile, selectStr ) 
        # copying is preferred

  # to take advantage of various leakSolve routines, create a (temporary) visibility dictionary
  # note that this is NOT the same as the SS object!
    vis = { "fileName" : "sstmp", "selectStr" : selectStr }
       # selectStr is superfluous since only data specified by selectStr were copied to sstmp
    leakSolve.makeFtable( vis )					# fill in the chfreq and chwidth items in vis dictionary
    if not vis["chfreq"] :						# this checks for empty list
      print "ABORTING THIS INTERVAL: sstmp is an empty file - perhaps all data are flagged"
      return False

  # process data in sstmp one channel range at a time
    self.strList = makeStrList( LkFile )	    # make list of available channel ranges in LkFile
    frq1 = []
    frq2 = []
    Stokes = []									# temporary work area
    for lineStr in self.strList :
      [fstartLeak,fstopLeak] = leakSolve.restoreLk( LkFile, lineStr, "sstmp" )
         # ... copy leakage for one lineStr from LkFile to sstmp

    # check that freq range for chans 'lineStr' is same in LkFile and in sstmp
      a = lineStr.split(",")   # "chan,1,49,8"
      ch1 = int(a[2])
      ch2 = ch1 + int(a[3]) - 1
      [fstart,fstop] = leakSolve.fminfmax( vis, ch1, ch2 )
      if (abs(fstartLeak - fstart) > .002) or (abs(fstopLeak - fstop) > .002) :
        print "    warning: freq mismatch:  Leak %.3f - %.3f,  Data %.3f - %.3f"  % (fstartLeak, fstopLeak, fstart, fstop) 

    # compute the stokes parameters for this 'lineStr'
      result = getStokes2( vis["fileName"], vis["selectStr"], lineStr )
      if numpy.isnan(result[0]) or numpy.isnan(result[2]) or numpy.isnan(result[4]) :
        print "skipping data with nan"
      else :
        frq1.append( fstart )
        frq2.append( fstop )
        Stokes.append( result )

    #print "Stokes = ",numpy.array(Stokes)
    #print "Stokes.transpose = ",numpy.array(Stokes).transpose()

  # copy data from temporary work area to SS object
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
    self.fitPARM( f0 )                 # fits p0, PA, RM, frac
    return True

  # -----------------------------------------------------------------------------------------------------------------------#
  # writes to file; if outfile=None, prints to screen
  def dump( self, outfile ) :
    outStr = "# visFile     : %s\n" % self.visFile 
    outStr = outStr + "# LkFile      : %s\n" % self.LkFile 
    outStr = outStr + "# selectStr   : %s\n" % self.selectStr 
    outStr = outStr + "# avg UT      : %.3f hrs\n" % self.UT 
    outStr = outStr + "# avg HA      : %.3f hrs\n" % self.HA 
    outStr = outStr + "# avg parang  : %.3f deg\n" % self.parang 
    outStr = outStr + "# p0          : %.4f  ( %.4f ) Jy\n" % (self.p0, self.p0rms) 
    outStr = outStr + "# frac        : %.3f  ( %.3f ) \n" % (self.frac, self.fracrms) 
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

  # -----------------------------------------------------------------------------------------------------------------------#
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
    print "entering fitPARM"
    freq = 0.5 * (self.f1 + self.f2)
    lambdasq = (.30*.30)/(freq*freq) - (.30*.30)/(freq0*freq0)
    x = numpy.concatenate( (lambdasq,lambdasq) )
    y = numpy.concatenate( (self.Q,self.U) )
    sigma = numpy.concatenate( (self.rmsQ,self.rmsU) )
    Qavg = numpy.average( self.Q )
    Uavg = numpy.average( self.U )
    #print "# Qavg,Uvg = ", Qavg, Uavg
    p0 = math.sqrt( Qavg*Qavg + Uavg*Uavg) 
    pa0 = 0.5 * math.atan2(Uavg,Qavg)
    print "# initial guess: p0 = %.3f, pa0 = %.2f, RM = 0." % (p0, pa0 * 180./math.pi)
    
    # --- fit data if there are measurements at at least 2 freqs --- #
    if len(self.f1) > 1 :
      popt,pcov = scipy.optimize.curve_fit( self.func, x, y, p0=[p0,pa0,0.], sigma=sigma  )
      print "popt = ", popt
      if popt[0] < 0. :                 # if fit amplitude is negative, add 180 to 2 * PA
        popt[0] = -1. * popt[0]
        popt[1] = popt[1] + math.pi/2.
      p0,pa0,rm = popt
      # compute chisq and scale covariances to get 'true' error estimates - see scipy thread
      #     by Christoph Deil; makes errors unreasonably small, so I am dropping it
      chi = (y - self.func( x, p0, pa0, rm))/sigma
      print "sigma = ",sigma
      print "chi = ",chi
      chisq = (chi**2).sum()
      print "chisq = %.2e" % chisq
      dof = len(x) - len(popt)
      #  pcov = pcov/(chisq/dof)
      try :
        p0rms = math.sqrt(pcov[0][0])
        pa0rms = math.sqrt(pcov[1][1])
        rmrms = math.sqrt(pcov[2][2]) 
      except :
        p0rms = pa0rms = rmrms = 0.
      #self.ran( x, y, p0, p0rms, pa0, pa0rms, rm, rmrms, sigma )  

    # --- for DSB data, one freq only, do not fit RM, estimate errors from Qrms and Urms ---
    else :
      p0rms = self.rmsI
      pa0 = pa0 
      pa0rms = 0.5 * math.sqrt( pow(self.Q*self.rmsU, 2.) + pow(self.U*self.rmsQ, 2.)) \
        / ( pow(self.Q,2.) + pow(self.U,2.) )
      pa0rms = pa0rms
      rm = 0.
      rmrms = 0.
      freq0 = freq
    
    self.p0 = p0
    self.p0rms = p0rms
    self.pa0 = pa0 * 180./math.pi
    self.pa0rms = pa0rms * 180./math.pi
    self.RM = rm
    self.RMrms = rmrms
    self.freq0 = freq0
    self.calcFrac()                            # fills in fractional pol
    print "# p0 = %.3f (%.3f)" % (self.p0, self.p0rms)
    print "# frac = %.3f (%.3f)" % (self.frac, self.fracrms)
    print "# pa0 = %.2f (%.2f)" % (self.pa0, self.pa0rms)
    print "# RM = %.3e (%.3e)" % (self.RM, self.RMrms)
    print ""

# returns changes that double chisq
# x and y are the independent and dependent variables
# p0, pa0, rm are the best fit values
  def ran( self, x, y, p0, p0rms, pa0, pa0rms, rm, rmrms, sigma ) : 
    for p0test in numpy.arange( p0-3.*p0rms, p0+4*p0rms, p0rms) :
      for pa0test in numpy.arange( pa0-3.*pa0rms, pa0+4.*pa0rms, pa0rms ) :
        for rmtest in numpy.arange( rm-3.*rmrms, rm+4.*rmrms, rmrms ) :
          chi = (y - self.func( x, p0test, pa0test, rmtest))/sigma
          chisq = (chi**2).sum()
          print " %5.3f %5.1f %7.2f %9.1f" % (p0test,pa0test*180./math.pi,rmtest/1.e5,chisq)

  def calcFrac( self ) :
    'compute fractional polarization and uncertainty'
    #Iavg = numpy.average(ss.I, weights=ss.rmsI )
    Iavg = numpy.average( self.I )
    Istd = numpy.std( self.I, ddof=1 )/math.sqrt(len(self.I))
    self.frac = self.p0/Iavg
    self.fracrms = self.frac * math.sqrt( pow(Istd/Iavg,2.) + pow(self.p0rms/self.p0,2.) )
    #print "Iavg,Istd,frac,fracrms = %.4f, %.4f, %.4f, %.4f" % (Iavg, Istd, self.frac,self.fracrms)

  def plot( self ) :
    pyplot.clf()
    pyplot.ion()
    Ymax = abs(1.5 * numpy.concatenate( (self.Q, self.U, -1.*self.Q, -1*self.U) ).max())
    fig = pyplot.subplot(1,1,1)
    self.plot2( fig, Ymax, labelsize=12 ) 
    pyplot.draw()                 # plots and continues...

  def plot2( self, fig, Ymax, labelsize ) :
    fmin = numpy.concatenate( (self.f1,self.f2) ).min()
    fmax = numpy.concatenate( (self.f1,self.f2) ).max()
    fmin = fmin - 0.1 * (fmax - fmin)
    fmax = fmax + 0.1 * (fmax - fmin)
    if Ymax == 0. :
      Ymax = abs(1.5 * numpy.concatenate( (self.Q, self.U, -1.*self.Q, -1*self.U) ).max())
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
    textStr = "pol = %.2f (%.2f)%%\nPA = %.1f (%.1f)\nRM = %.1e (%.1e)" % \
      ( 100.*self.frac, 100*self.fracrms, self.pa0 ,self.pa0rms, self.RM, self.RMrms)
    props = dict(boxstyle='round', facecolor='wheat', alpha=1 )	# patch.Patch properties
    fig.text( .5, .05, textStr, transform=fig.transAxes, horizontalalignment='center', 
      verticalalignment='bottom', fontsize=0.8*labelsize, bbox=props )
    pyplot.title( self.selectStr, size=labelsize, y=.9+.05*labelsize/5. )  

# -----------------------------------------------------------------------------------------------------------------------#
# generate series of SS objects (usually for one source, with different Lks, sometimes for different time intervals)
# pickle the results and append to outfile
# renamed - was makeTable
# change: LkListFile can be either the name of a single LkFile, or a list of LkFile names

def appendToSSfile( visFile, LkListFile, srcName, extra, nint, maxgap, outfile, schedFile=None, plot=False ) :
  # Generate list of leakages
  LkList = []
  first = True
  fin = open( LkListFile, "r" )
  for line in fin :
    if first and line.startswith("# input data") :   # then LkListFile is itself a single Lk file
      LkList.append( LkListFile )
      break
    else :                                           # then LkListFile is a list of Lk file names
      first = False
      if not line.startswith("#") :
        a = line.split()
        if len(a) > 0 :
          LkList.append( a[0] )
  fin.close()
  print "LkList = ", LkList
  # Generate list of selectStrings (time ranges, usually)
  if schedFile :
    print "making list according to vlbi schedFile %s" % schedFile 
    selectList = paPlot.makeSelectList2( schedFile, srcName )
  else :
    print "include up to %d integrations; no gaps > %.1f minutes" % (nint, maxgap)
    selectList = paPlot.makeSelectList( visFile, srcName, nint, maxgap )
  print "selectList = ",selectList
  # Generate SS object for each time range, each leakage; append to outfile 
  nlist = len(selectList)
  for n, selectString in enumerate(selectList) :
    if len( extra ) > 0 :
      selectString = selectString + "," + extra
    print " "
    print "%d/%d  %s" % (n+1,nlist,selectString)
    print "Copying (but not averaging) selected data to temporary file sstmp"
    leakSolve.copytoSStmp( visFile, selectString ) 
    for Lk in LkList :
      ss = SS()
      if ss.make( "sstmp", selectString, Lk, preavg=None ) :		# returns False if it fails
        fout = open(outfile,"ab")
        pickle.dump( ss, fout )
        fout.close()
        if plot : 
          ss.plot()
                                                                                                                   
# read in series of SSs from input pickle file, append to ssList (was readAll2)
def readSSfile( infile, ssList ) :
  fin = open( infile, "rb" )
  while True :
    try :
      ss = pickle.load(fin)
      ss.dump(None)
      ssList.append(ss)
    except (EOFError, pickle.UnpicklingError ) :
      break
  fin.close()

def findUniqueStrings( ssList ) :
  uniqueStrList = []
  for ss in ssList :
    if not ss.selectStr in uniqueStrList :
      uniqueStrList.append(ss.selectStr)
  return uniqueStrList

# read in many SS spectra - usually, these were fit with different leakages
# for each unique selectStr, concatenate all the I,Q,U,V data into one SS, fit p, PA, frac, RM 
# the hope was that this global fit would be more accurate than individual fits
def fitAll( infile, outfile, plot=False ) :
  ssList = []
  readSSfile( infile, ssList )
  for ss in ssList :
    print ss.selectStr
  # Make list of unique selectStr; note this will usually be different time range
  uniqueList = findUniqueStrings( ssList )
  for uniqueStr in uniqueList :
    first = True
    for ss in ssList :
      if ss.selectStr == uniqueStr :
        if first :
          ssWork = copy.deepcopy(ss)
          first = False
        else :
          print "appending values from %s" % ss.LkFile
          tmp = numpy.append( ssWork.f1, ss.f1 )
          ssWork.f1 = tmp
          tmp = numpy.append( ssWork.f2, ss.f2 )
          ssWork.f2 = tmp
          tmp = numpy.append( ssWork.I, ss.I )
          ssWork.I = tmp
          tmp = numpy.append( ssWork.rmsI, ss.rmsI )
          ssWork.rmsI = tmp
          tmp = numpy.append( ssWork.Q, ss.Q )
          ssWork.Q = tmp
          tmp = numpy.append( ssWork.rmsQ, ss.rmsQ )
          ssWork.rmsQ = tmp
          tmp = numpy.append( ssWork.U, ss.U )
          ssWork.U = tmp
          tmp = numpy.append( ssWork.rmsU, ss.rmsU )
          ssWork.rmsU = tmp
          tmp = numpy.append( ssWork.V, ss.V )
          ssWork.V = tmp
          tmp = numpy.append( ssWork.rmsV, ss.rmsV )
          ssWork.rmsV = tmp
    f0 = 226.
    ssWork.fitPARM( f0 )                 # fits p0, PA, RM, frac
    ssWork.LkFile = "various"
    fout = open(outfile,"ab")
    pickle.dump( ssWork, fout )
    fout.close()
    if plot : 
      ssWork.plot()

# weighted mean and its error (Bevington+Robinson, pp 58-59)
# .. weighted mean = sum(x_i/sigma_i^2) / sum(1./sigma_i^2)
# .. error in the weighted mean = sqrt(1./sum(1/sigma_i^2))
def weightedMeanAndError( data, rms ) :
  wgts = numpy.array( 1./(rms*rms) )
  avg = numpy.average( data, weights=wgts )
  avgRms = math.sqrt(1./numpy.sum(wgts))
  return (avg, avgRms)

# read in file of pickled SS spectra
# for each unique selectStr, tabulate the fit results, compute means and rms scatter
def summarizeFits( infile, outfile, plot=False ) :
  ssList = []
  readSSfile( infile, ssList )
  for ss in ssList :
    print ss.selectStr
  fout = open( outfile, "w")
  # Make list of unique selectStr - may be different sources, time ranges, etc.
  uniqueList = findUniqueStrings( ssList )
  for uniqueStr in uniqueList :
    I = []
    V = []
    p0 = []
    pa = []
    rm = []
    frac = []
    first = True
    for ss in ssList :
      if ss.selectStr == uniqueStr :
        if first :
          fout.write("\n#> %s   avgUT = %.3f\n" % (uniqueStr,ss.UT))
          fout.write("#     S    sigma   poli  sigma     PA  sigma      RM  sigma   frac  sigma   V sigma  LkFile\n")
          first = False
        [Iavg,Istd] = weightedMeanAndError( ss.I, ss.rmsI )
        [Vavg,Vstd] = weightedMeanAndError( ss.V, ss.rmsV )
        Iavg = numpy.mean(ss.I)
        Istd = numpy.std(ss.I)
        Vavg = numpy.mean(ss.V)
        Vstd = numpy.std(ss.V)
        I.append(Iavg)
        V.append(Vavg)
        p0.append(ss.p0)
        pa.append(ss.pa0)
        rm.append(ss.RM)
        frac.append(ss.frac)
        fout.write(" %8.3f %6.3f %7.3f %5.3f %7.1f %5.1f %8.2f %5.2f %7.3f %5.3f %7.3f %5.3f  %s \n" % \
          ( Iavg, Istd, ss.p0, ss.p0rms, ss.pa0, ss.pa0rms, ss.RM/1.e5, ss.RMrms/1.e5, \
          ss.frac, ss.fracrms, Vavg, Vstd, ss.LkFile) )
  # Print averages and rms; here the rms is from the dispersion in the means!
    fout.write("#* %6.3f %6.3f %7.3f %5.3f %7.1f %5.1f %8.2f %5.2f %7.3f %5.3f %7.3f %5.3f  AVG\n" % \
      ( numpy.mean(I), numpy.std(I), numpy.mean(p0), numpy.std(p0), numpy.mean(pa), numpy.std(pa), \
      numpy.mean(rm)/1.e5, numpy.std(rm)/1.e5, numpy.mean(frac), numpy.std(frac), numpy.mean(V), numpy.std(V) ) )
  # 3D plot of fits
    #print "begin figure"
    #fig = pyplot.figure()
    #ax = Axes3D(fig)
    #q = Axes3D.scatter(ax,numpy.array(p0),numpy.array(pa),numpy.array(rm) )
    #print "show figure"
    #pyplot.show()
  fout.close()

# read in series of SSs from input file (old style, pre-pickle)
#def readAll(  infile, ssList ) :
#  fin = open( infile, "r" )
#  ss = SS()
#  f1, f2, I, rmsI, Q, rmsQ, U, rmsU, V, rmsV = ( [] for i in range(10) )
#
#  for line in fin :
#    a = line.split()
#    if line == "# ----------------------------\n" :     # this delimits end of each SS list
#      ss.f1 = numpy.array(f1)
#      ss.f2 = numpy.array(f2)
#      ss.I = numpy.array(I)
#      ss.rmsI = numpy.array(rmsI)
#      ss.Q = numpy.array(Q)
#      ss.rmsQ = numpy.array(rmsQ)
#      ss.U = numpy.array(U)
#      ss.rmsU = numpy.array(rmsU)
#      ss.V = numpy.array(V)
#      ss.rmsV = numpy.array(rmsV)
#      ssList.append(ss)       # add completed SS object to the list
#      ss = SS()               # initialize a new SS object
#      f1, f2, I, rmsI, Q, rmsQ, U, rmsU, V, rmsV = ( [] for i in range(10) )
#    elif line.startswith("# visFile") : ss.visFile = a[3]
#    elif line.startswith("# LkFile") : ss.LkFile = a[3]
#    elif line.startswith("# selectStr") : ss.selectStr = a[3]
#    elif line.startswith("# avg UT") : ss.UT = float( a[4] )
#    elif line.startswith("# avg HA") : ss.HA = float( a[4] )
#    elif line.startswith("# avg parang") : ss.parang = float( a[4] )
#    elif line.startswith("# p0") : 
#      ss.p0 = float( a[3] )
#      ss.p0rms = float( a[5] )
#    elif line.startswith("# frac") : 
#      ss.frac = float( a[3] )
#      ss.fracrms = float( a[5] )
#    elif line.startswith("# PAfit") : 
#      ss.pa0 = float( a[3] )
#      ss.pa0rms = float( a[5] )
#    elif line.startswith("# RMfit") : 
#      ss.RM = float( a[3] ) * 1.e5
#      ss.RMrms = float( a[5] ) * 1.e5
#    elif line.startswith("# PAfreq0") : ss.freq0 = float( a[3] )
#    elif not line.startswith("#") :
#      f1.append( float( a[0] ) )
#      f2.append( float( a[1] ) )
#      I.append( float(a[2]) )
#      rmsI.append( float( a[3] ) )
#      Q.append( float( a[4] ) )
#      rmsQ.append( float( a[5] ) )
#      U.append( float( a[6] ) )
#      rmsU.append( float( a[7] ) )
#      V.append( float( a[8] ) )
#      rmsV.append( float( a[9] ) )
#      ss.strList.append( a[10] ) 

# read in one or more PA files, write summary file for wip (this is geared toward vlbi)
def summary( paList, outfile ) :
  srcList = [ "0234+285", "3c84", "3c111", "0721+713", "0854+201", "OJ287", "3c273", "M87", \
	"3c279", "1337-129", "1633+382", "3c345", "NRAO530", "SgrA", "1749+096", "1921-293", \
    "bllac", "titan", "1743-038", "3c286" ]
  colorList = [ 3,2,14,4,3,3,4,8,1,9,10,11,12,13,5,15,1,2,3,4]
  nstart = 1
  fin = open( paList, "r" )
  for line in fin :
    a = line.split()
    if not line.startswith("#") and len(a) > 0 :
      ssList = []
      infile = a[0]
      srcName = a[0]
      color = 1
      for src,col in zip (srcList,colorList) :
        if (a[0].find(src) >= 0) or (a[0].find(src.upper()) >= 0) :
          srcName = src
          color = col
      fout = open( outfile, "a" )
      fout.write("#\n")
      fout.write("# %s\n" % infile)
      fout.write("#  dechr  parang    HA        S    sigma   poli  sigma     PA  sigma     RM  sigma   frac  sigma  col   selectString\n")
      fout.close()

      fout = open( "QUVsummary", "a" )
      fout.write("#\n")
      fout.write("# %s\n" % infile)
      fout.write("#  dechr  parang    HA        S    sigma     Q   sigma      U  sigma      V  sigma    col   selectString\n")
      fout.close()

      readSSfile( infile, ssList )
      for ss in ssList : 
        #Iavg = numpy.average(ss.I, weights=ss.rmsI )
        Iavg = numpy.average(ss.I )
        Istd = numpy.std(ss.I, ddof=1)
        fout = open(outfile, 'a')
        fout.write("%8.3f %7.2f %7.3f %8.3f %6.3f %7.3f %5.3f %7.1f %5.1f %8.2f %5.2f %7.3f %5.3f %2d   %s\n" % \
          (ss.UT, ss.parang, ss.HA, Iavg, Istd, ss.p0, ss.p0rms, ss.pa0, ss.pa0rms, ss.RM/1.e5, ss.RMrms/1.e5, \
          ss.frac, ss.fracrms, color, ss.selectStr) )
        fout.close()

        fout = open("QUVsummary", 'a')
        #fout.write("%8.3f %7.2f %7.3f %8.3f %6.3f %7.3f %5.3f %7.3f %5.3f %7.3f %5.3f  %2d   %s\n" % \
        #  (ss.UT, ss.parang, ss.HA, Iavg, Istd, ss.Q, ss.rmsQ, ss.U, ss.rmsU, ss.V, ss.rmsV, color, ss.selectStr) )
        fout.close()

      nstop = nstart + len(ssList) + 2
      print "1src %d %d %d %s" % (nstart, nstop, color, srcName )
      nstart = nstop + 1
  fin.close()
        
# read in one or more PA files, write summary file with Q,U for wip
def DSBsummary( paList, outfile ) :
  srcList = [ "0234+285", "3c84", "3c111", "0721+713", "0854+201", "OJ287", "3c273", "M87", \
	"3c279", "1337-129", "1633+382", "3c345", "NRAO530", "SgrA", "1749+096", "1921-293", \
    "bllac", "titan", "1743-038", "3c286" ]
  colorList = [ 3,2,14,4,3,3,4,8,1,9,10,11,12,13,5,15,1,2,3,4]
  nstart = 1
  fin = open( paList, "r" )
  for line in fin :
    a = line.split()
    if not line.startswith("#") and len(a) > 0 :
      ssList = []
      infile = a[0]
      srcName = a[0]
      color = 1
      for src,col in zip (srcList,colorList) :
        if (a[0].find(src) >= 0) or (a[0].find(src.upper()) >= 0) :
          srcName = src
          color = col
      fout = open( outfile, "a" )
      fout.write("#\n")
      fout.write("# %s\n" % infile)
      fout.write("#  dechr  parang    HA        S    sigma   poli  sigma     PA  sigma     RM  sigma   frac  sigma  col   selectString\n")
      fout.close()
      readAll( infile, ssList )
      for ss in ssList : 
        Iavg = numpy.average(ss.I )
        Istd = numpy.average(ss.rmsI)
        Q = numpy.average(ss.Q)
        rmsQ = numpy.average(ss.rmsQ)
        U = numpy.average(ss.U)
        rmsU = numpy.average(ss.rmsU)
        V = numpy.average(ss.V)
        rmsV = numpy.average(ss.rmsV)
        fout = open(outfile, 'a')
        outStr = "%8.3f %7.2f %7.3f %8.3f %6.3f %7.3f %5.3f %7.1f %5.1f %7.3f %5.3f %2d" % \
          (ss.UT, ss.parang, ss.HA, Iavg, Istd, ss.p0, ss.p0rms, ss.pa0, ss.pa0rms, \
          ss.frac, ss.fracrms, color )
        outStr = outStr + "    %10.3e %10.3e   %10.3e %10.3e   %10.3e %10.3e   %s\n" \
         % (Q,rmsQ,U,rmsU,V,rmsV,ss.selectStr) 
        fout.write( "%s" % outStr )
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
    a = line.split()
    if (len(a) > 0) and (not line.startswith("#")) :
      ssList = []
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
  
def replot2( infile, Ymax=0., nrows=2, ncols=1 ) :
  ssList = []
  readSSfile( infile, ssList )
  pyplot.ioff()
  pp = PdfPages( 'multipage.pdf' )
  pyplot.subplots_adjust( hspace=0.25 )
  nplot = 0
  for ss in ssList :
    if nplot > nrows*ncols : 
       pyplot.savefig( pp, format='pdf' )
       nplot = 1
       pyplot.clf()
    fig = pyplot.subplot( nrows, ncols, nplot )
    ss.plot2( fig, Ymax, labelsize=10./ncols )
    nplot = nplot + 1
  pyplot.savefig( pp, format='pdf' )
  pp.close()
  
# plot Q,U with error bars over the course of a night
# color code the points (and arrows between them) by time?
# paFile contains a set of stokes objects

def QUtrajectory ( paFile ) :
  pyplot.ion()
  pyplot.clf()
  ssList = []
  readAll( paFile, ssList )
  UT = []
  Q = []
  Qerr = [] 
  U = []
  Uerr = [] 
  for ss in ssList : 
     UT.append(ss.UT)
     Q.append(numpy.average(ss.Q))
     #Qerr.append(numpy.average(ss.Q))
     Qerr.append(0.02)
     U.append(numpy.average(ss.U))      
     #Uerr.append(numpy.average(ss.U))      
     Uerr.append(0.03)      
  ymax = 1.1*numpy.amax(numpy.absolute([Q,U]))

  fig = pyplot.subplot(1,1,1)
  fig.axis( [-1.*ymax, ymax, -1.*ymax, ymax], size=3 )
  fig.tick_params(axis='both', which='major', labelsize=10 )
  fig.grid( True )
  fig.errorbar(Q,U,xerr=Qerr,yerr=Uerr,color="r",alpha=1, fmt=None,capsize=0, ecolor="0.2",linewidth=1,mew=0, zorder=1)
  fig.plot(Q,U,color="blue",alpha=0.3, zorder=2)
  fig.scatter(Q,U, s=80, c=UT, marker="o", linewidths=1, alpha=1, zorder=3, edgecolor="0.6" )
  #pyplot.colorbar(ax=fig,ticks=[min(UT),max(UT)] )
  #for n in range(0,len(Q)-1) :
    #pyplot.arrow(Q[n],U[n],Q[n+1]-Q[n],U[n+1]-U[n],length_includes_head=True, head_width=0.02*ymax)
    #pyplot.annotate( "", [Q[n],U[n]], xytext=[Q[n+1],U[n+1]], arrowprops=dict(linewidth=1) )
  pyplot.show()

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
   
