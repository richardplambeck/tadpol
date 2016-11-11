# ob.py
# these are routines for dealing with the Optics Bench scans of various kinds

import numpy
import math
import sys
import string
import matplotlib
matplotlib.use('GTKAgg')
import matplotlib.pyplot as pyplot
from matplotlib.backends.backend_pdf import PdfPages
import pol

clight = 29.9792458 # speed of light, cm/nanosec


def findminima( x, p, fGHz ):
    delta = 0.45*300./fGHz
      # anticipated spacing between minima
    xminlist = []
    pminlist = []
    while numpy.ma.count( x ) > 0 :
      nmin = numpy.argmin( p )
        # nmin is index of point with  minimum power
      xmin = x[nmin]
      xminlist.append( xmin )
      pminlist.append( p[nmin] )

    # now mask all points that are within +/-delta of the min
      x2 = numpy.ma.masked_inside( x, xmin-delta, xmin+delta )
      p2 = numpy.ma.masked_where( numpy.ma.getmask(x2), p )

    # copy the arrays to allow easy iteration
      x = x2
      p = p2

    print "finished - all masked"
    indexsorted = numpy.argsort(xminlist)
    xminsorted = []
    pminsorted = []
    for n in indexsorted :
      xminsorted.append( xminlist[n] )
      pminsorted.append( pminlist[n] )
    print numpy.array_str( numpy.array(xminsorted), precision=4 ) 
    gap = []
    for n in range(1,len(xminsorted)) :
      gap.append( xminsorted[n] - xminsorted[n-1] )
    print numpy.array_str( numpy.array(gap), precision=4 ) 
    return [xminsorted, pminsorted]

def processFile( infile ) :

    [x,p,fGHz] = readscan( infile )
    print "fGHz = %.2f; lambda/2 = %.4f mm" % (fGHz,  0.5*299.8/fGHz)

    x1 = 0.
    while (x1 >= 0.) :
      xminsorted,pminsorted = findminima( x, p, fGHz )
      print "enter -1 to get out of loop"
      startstop = raw_input( "startpos stoppos (mm): ") 
      a = startstop.split()
      try :
        x1 = float(a[0])
      except:
        x1 = 0.
      try :
        x2 = float(a[1])
      except :
        x2 = 20.
      print "averaging values from %.3f to %.3f mm" % (x1,x2)
      navg = 0
      avg = 0.
      for i in range(0,len(x)) :
        if (x[i] >= x1) and (x[i] <= x2) :
            avg = avg + p[i]
            navg = navg + 1
      print "%s: n = %d, avg = %.5f" % (infile, navg, avg/navg )
      print " "

    # append values to csv file
      fout = open( "summary2.csv", "a" )
      fout.write("%s,%.4f,%.4f,%d,%.5f\n" % (infile, x1, x2, navg, avg/navg))
      fout.close()

def readscan( infile ) :
    xin = []
    pin = []
    try :
      n = infile.find("_")
      fGHz = float( infile[0:n] )
    except :
      fGHz = float( raw_input( "fGHz: ") )
    fin = open( infile, "r" )
    for line in fin :
      if line.startswith("#") :
        print line
      else :
        a = line.split()
        xin.append( float(a[1]) )
        pin.append( float(a[2]) )
    return numpy.array(xin), numpy.array(pin), fGHz
    
def process( scanList=["100_out","100_posA","100_posB"], plotName="100" ) :
    pyplot.ioff()
    pp = PdfPages("%s.pdf" % plotName )
    npanels = 3
    fig = pyplot.figure( figsize=(8,11) )

    npanel = 0
    for scan in scanList: 
      npanel = npanel + 1
      [x,p,fGHz] = readscan( scan )
      s = pyplot.subplot(npanels,1,npanel)
      s.grid( True, linewidth=0.05, color="0.01" )   # color=0.1 is a light gray

      if npanel == 1 :
        xmin = x.min() - 0.04*(x.max()-x.min()) 
        xmax = x.max() + 0.04*(x.max()-x.min()) 
      pmin = p.min() - 0.04*(p.max()-p.min()) 
      pmax = p.max() + 0.06*(p.max()-p.min()) 

      s.axis( [xmin, xmax, pmin, pmax] )
      s.plot( x, p, color="r", linewidth=0.5 )
      xminsorted,pminsorted = findminima( x, p, fGHz )
      s.plot( xminsorted, pminsorted, "bo" )
      s.text(.03, .92,  "%s" % scan, horizontalalignment='left', transform=s.transAxes, fontsize=12 )

    pyplot.savefig( pp, format='pdf' )
    pp.close()
    pyplot.show()

def calib( infile, outfile ) :
  fin = open( infile, "r" )
  fout = open( outfile, "w" )
  fout.write("# calibration from %s\n" % infile)
  first = True
  for line in fin :
    a = line.split()
    if a[0] == "#" :
      dB = float( a[1] )
      if first :
        dBref = dB  
      v = []
    else :
      v.append( float( a[0] ))
      if len(v) >= 2 :
        vavg = (v[0] + v[1])/2.
        trueTrans = pow(10., -(dB-dBref)/10.)
        if first :
          vref = vavg 
          first = False 
        measRatio = vavg/vref
        print "%.2f   %.5f   %.5f" % (dB,trueTrans,measRatio)
        fout.write( "%.2f   %.5f   %.5f\n" % (dB,trueTrans,measRatio) )
  fout.close()

def doit() :
  process( ["70_out1", "70_posA", "70_posB"], "70a" )
  process( ["70_out2", "70_posA", "70_posB"], "70b" )
  process( ["80_out", "80_posA", "80_posB"], "80" )
  process( ["90_out", "90_posA", "90_posB"], "90" )
  process( ["100_out", "100_posA", "100_posB"], "100" )
  process( ["110_out", "110_posA", "110_posB"], "110a" )
  process( ["110_outrepeat", "110_posA", "110_posB"], "110b" )

def doit2() :
  calib( "cal_70", "calib_70" )
  calib( "cal_80", "calib_80" )
  calib( "cal_90", "calib_90" )
  calib( "cal_100", "calib_100" )
  calib( "cal_110", "calib_110" )

def doit3() :
  processFile( "70_out1" )
  processFile( "70_out2" )
  processFile( "70_posA" )
  processFile( "70_posB" )
  processFile( "80_out" )
  processFile( "80_posA" )
  processFile( "80_posB" )
  processFile( "90_out" )
  processFile( "90_posA" )
  processFile( "90_posB" )
  processFile( "100_out" )
  processFile( "100_posA" )
  processFile( "100_posB" )
  processFile( "110_outrepeat" )
  processFile( "110_out" )
  processFile( "110_posA" )
  processFile( "110_posB" )

# finds offset of phase center from scan chord
def fitPhase( infile ) :
  fin = open( infile, "r" )
  xin = []
  cin = []
  sin = []
  for line in fin :
    if not line.startswith( "#" ) :
      a = line.split()
      x.append( float(a[0]) )
      c.append( float(a[1]) )
      s.append( float(a[2]) )
  fin.close()
  
  x = numpy.array(xin)
  c = numpy.array(cin)
  s = numpy.array(sin)

  coffset = c.max() + c.min()
  crange = c.max() - c.min()
  soffset = s.max() + s.min()
  srange = s.max() - s.min()
  print "cos max,min = %.6f,%.6f" % (c.max(),c.min())
  print "sin max,min = %.6f,%.6f" % (s.max(),s.min())
 
  phase = numpy.unwrap(numpy.atan2( c-coffset, crange/srange*(s-soffset)))
  

#  tpuckcm = .2525 * 2.54
#  nrefrac = 3.105
#      dphs_est = 360.*fGHz/29.98 * tpuckcm * (nrefrac - 1.)
#      phs1_est = (phs0 + dphs_est) % 360.
#      delta1 = phs1 - phs1_est
#      delta = phs1 - phs0
#      while delta1 < -180. : 
#        delta1 = delta1 + 360.
#      while delta1 > 180. :
#        delta1 = delta1 - 360.
#      while delta < -180. : 
#        delta = delta + 360.
#      while delta > 180. :
#        delta = delta - 360.
#      print fGHz, delta
#  fout = open("/o/plambeck/PolarBear/OpticsBench/20oct2016/alumina.dat", "w")

# reads summary file and returns arrays of freq, ampratio, deltaPhase
def processSummary_20oct( infile="/o/plambeck/PolarBear/OpticsBench/20oct2016/summ" ) :
  fin = open(infile, "r")
  fGHz = []
  trans = []
  dphi = []
  nseq = 0
  for line in fin :
    nseq = nseq + 1
    a = line.split()
    if nseq == 1 :
      if not line.startswith("freq:") :
        print "WARNING - file is out of sync"
      fGHz.append( float(a[1])/1000. )
    elif nseq == 2 : 
      phs0 = float( a[5] )  
      amp0 = float( a[4] )
    elif nseq == 3 :
      phs1 = float( a[5] )  
      amp1 = float( a[4] )
      nseq = 0
      trans.append( amp1/amp0 )
      delta = phs1 - phs0
      while delta < -180. : 
        delta = delta + 360.
      while delta > 180. :
        delta = delta - 360.
      dphi.append( delta )
  return numpy.array(fGHz), numpy.array(trans), numpy.array(dphi)
  fin.close()

def processSummary( infile="/o/plambeck/PolarBear/OpticsBench/09nov2016/transmission.dat" ) :
  fin = open(infile, "r")
  fGHz = []
  trans = []
  dphi = []
  for line in fin :
    a = line.split()
    fGHz.append( float( a[0] ) )
    trans.append( float( a[1] ) )
    dphi.append( float( a[2] ) )
  return numpy.array(fGHz), numpy.array(trans), numpy.array(dphi)
  fin.close()
    
#numpy.var( numpy.unwrap( deltaPhase ) )

# returns expected transmission, deltaPhase for 1pass through plate (no FP reflections)
# fGHz = freq in GHz
# nrefrac = estimated refractive index
# tcm = puck thickness in cm
# angIdeg = angle of incidence, degrees
# tanDelta = estimated loss tangent
# npar = True if incident radiation is polarized parallel to plane of incidence, False if perp

def fit1pass( fGHz, nrefrac, tcm, angIdeg, tanDelta, npar=True ) :
    thetaI = numpy.radians(angIdeg)
    thetaT = math.asin((math.sin(numpy.radians(angIdeg)))/nrefrac)
    tpar1, tperp1, rpar, rperp = pol.fresnel( 1., thetaI, nrefrac, thetaT )
    tpar2, tperp2, rpar, rperp = pol.fresnel( nrefrac, thetaT, 1., thetaI )
    dphs_est = 360.*fGHz/clight * tcm/math.cos(thetaT) * (nrefrac - 1.)
    return dphs_est
 
def fitFP( fGHz, nrefrac, tcm, angIdeg, tanDelta, npar=True ) :
    phs = [] 
    trans = []
    for freq in fGHz :
      ampTpar,ampTperp,ampRpar,ampRperp = pol.plateTrans( freq, nrefrac, nrefrac, angIdeg=angIdeg, tcm=tcm, tanDelta=tanDelta)
      trans.append( abs( ampTpar ) )
      phs.append( numpy.angle( ampTpar, deg=True ) - 360.*freq/clight*tcm )
    return unwrap(numpy.array(phs) ), numpy.array(trans)

def unwrap( phi ) :
    for n in range(0,len(phi)) :
      while phi[n] < -180. : 
        phi[n] = phi[n] + 360.
      while phi[n] > 180. :
        phi[n] = phi[n] - 360.
    return phi
     
# fits refractive index of dielectric slab from phase differences (slab in - slab out)
# retrieves measurements using processSummary() - edit that routine as necessary
# models amp and phase of transmitted signal assuming 1 pass through the slab,
#    OR multiple passes using Fabry Perot etalon formulae
# resolves 2pi ambiguities by keeping data-model phase differences between -180 and 180
# inputs: tcm = thickness of slab in cm; angIdeg = angle of incidence in degrees;
#    tanDelta = estimate of loss tangent; guess = estimate of refractive index,
#    plotType = 0 for none, 1 for phases and residuals vs freq; 2 for resids vs fit, 3 for raw data
#
def nrefracFit( tcm=.6406, angIdeg=0., nguess=3.11, tanDelta=4.e-4, plotType=0 ) :
  pyplot.ioff()
  pp = PdfPages("puck.pdf")
  fGHz,trans,dphs_meas = processSummary()
     # read in the data; edit as necessary as data format changes
  if plotType == 3 :
       pyplot.figure( figsize=(11,8) )
       ax = pyplot.subplot(2,1,1)
       ax.axis( [72.,114.,-195.,195.] )
       ax.plot( fGHz, dphs_meas, "o" )
       #ax.plot( fGHz, dphs_est, "r-" )
       ax.set_ylabel("$\Delta\phi$ (deg)", fontsize=14)
       pyplot.title( "PB2bc alumina sample, .2522 inch thick", fontsize=14 ) 
       pyplot.grid(True)
       ax = pyplot.subplot(2,1,2)
       ax.axis( [72.,114.,0.,1.5] )
       ax.plot( fGHz, trans, "o" )
       ax.plot( [10.,1000.],[1.,1.], "--" )
       ax.set_xlabel("freq (GHz)", fontsize=14)
       ax.set_ylabel("transmission", fontsize=14)
       pyplot.grid(True)
       pyplot.savefig( pp, format="pdf" )
       pyplot.show()
  nr = []
  avg = []
  std = []
  imin = 0
  for nrefrac in numpy.arange( nguess-.05, nguess+.06, .05) :
     nr.append( nrefrac )
     # dphs_est = unwrap( fit1pass( fGHz, nrefrac, tcm, angIdeg, tanDelta ) )
     dphs_est,trans_est = fitFP( fGHz, nrefrac, tcm, angIdeg, tanDelta )
     dphi = unwrap( dphs_meas-dphs_est )
       # unwrap keeps phase in range -180 to 180 in all cases
     avg.append( numpy.average(dphi) )
     if abs(avg[-1]) < abs(avg[imin]) :
       imin = len(avg)-1
     std.append( numpy.std(dphi) / math.sqrt(len(fGHz)) )
     print nrefrac, avg[-1], std[-1]
     if plotType == 1 or plotType == 4 :
       pyplot.figure( figsize=(11,8) )
       ax = pyplot.subplot(2,1,1)
       ax.axis( [72.,114.,-195.,195.] )
       ax.plot( fGHz, dphs_meas, "o" )
       ax.plot( fGHz, dphs_est, "r-" )
       ax.set_ylabel("phase (deg)", fontsize=14)
       pyplot.title( "PB2bc alumina sample 0.2522 in thick", fontsize=14 ) 
       pyplot.grid(True)
       if plotType == 1 :
         ax = pyplot.subplot(2,1,2)
         ax.axis( [72.,114.,-50.,50.] )
         ax.plot( fGHz, dphi, "o" )
         ax.plot( [75.,112.],[avg[-1],avg[-1]],"r-" )
         ax.set_xlabel("freq (GHz)", fontsize=14)
         ax.set_ylabel("residual (deg)", fontsize=14)
         ax.text( 0.05, .82, "n = %.3f" % nrefrac, transform=ax.transAxes, \
           horizontalalignment='left', fontsize=16, color="black", rotation='horizontal' )
         pyplot.grid(True)
       else :
         ax = pyplot.subplot(2,1,2)
         ax.axis( [72.,114.,0.,1.5] )
         ax.plot( fGHz, trans, "o" )
         ax.plot( fGHz, trans_est, "r-" )
         ax.set_xlabel("freq (GHz)", fontsize=14)
         ax.set_ylabel("transmission", fontsize=14)
         ax.text( 0.05, .82, "n = %.3f" % nrefrac, transform=ax.transAxes, \
           horizontalalignment='left', fontsize=16, color="black", rotation='horizontal' )
         pyplot.grid(True)

       pyplot.savefig( pp, format="pdf" )
       pyplot.show()
  slope = (avg[imin+1] - avg[imin-1])/(nr[imin+1] - nr[imin-1])
  nrfit = nr[imin] - avg[imin]/slope  
  fitunc = abs(2.*std[imin]/slope)
  print "best fit is nrefrac = %.3f (+/- %.3f)" % (nrfit,fitunc)
  
  if plotType == 2 :
    pyplot.figure( figsize=(11,8) )
    ax = pyplot.subplot(1,1,1)
    ax.axis( [3.05,3.17,-50.,50.], linewidth=2 )
    ax.errorbar( nr, avg, yerr=std, linewidth=2 )
    ax.errorbar( nrfit, 0., xerr=fitunc, color='red', elinewidth=6 )
    #ax.set_xticklabels( fontsize=16 )
    ax.set_xlabel("refractive index", fontsize=18)
    ax.set_ylabel("mean phase residual (deg)", fontsize=18)
    ax.text( 0.95, .92, "best fit = %.3f (+/- %.3f)" % (nrfit,fitunc), transform=ax.transAxes, \
          horizontalalignment='right', fontsize=18, color="black", rotation='horizontal' )
    pyplot.title( "PB2bc alumina sample 0.2522 in thick", fontsize=18 ) 
    pyplot.grid(True)
    pyplot.savefig( pp, format="pdf" )
    pyplot.show()
  if plotType > 0 :
    pp.close()


# this is just used to make drawing
def sineWaves() :
  x = numpy.arange(0,50.,.01)
  y = numpy.sin(x)
  pyplot.ioff()
  pp = PdfPages("sinewave.pdf")
  pyplot.figure( figsize=(11,8) )
  ax = pyplot.subplot(3,1,2)
  ax.axis( [-5.,55.,-1.5,1.5] )
  ax.plot( x, y, "-", linewidth=2 )
  pyplot.savefig( pp, format="pdf" )
  pyplot.show()
  pp.close()
  
