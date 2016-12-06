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
#import hou

# lambda is alternate way of defining simple function
sq = lambda arg: pow(arg, 2)
sin = lambda arg: numpy.sin(arg)
cos = lambda arg: numpy.cos(arg)

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


# reads *.trans file and returns arrays of freq, ampratio, deltaPhase, title, tcm, angIdeg
def readTransFile( infile ) :
  fin = open(infile, "r")
  fGHz = []
  trans = []
  dphi = []
  tcm = 0.
  title = ""
  angIdeg = 0.
  for line in fin :
    if not line.startswith( "#" ) :
      a = line.split()
      fGHz.append( float( a[0] ) )
      trans.append( float( a[1] ) )
      dphi.append( float( a[2] ) )
    elif "tcm" in line :
      tcm = float(line[line.find("=")+1:])
    elif "angIdeg" in line :
      angIdeg = float( line[line.find("=")+1:] )
    elif "title" in line :
      title = line[line.find("=")+1:]
      print len(title)
      title = title[2:-2]
      print len(title)
  return numpy.array(fGHz), numpy.array(trans), numpy.array(dphi), title, tcm, angIdeg
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
 
def fitFP( fGHz, nrefrac, tcm, angIdeg, tanDelta, npar=False ) :
    # print "tanDelta = %.2e" % tanDelta
    theta1 = numpy.radians(angIdeg)
    theta2 = math.asin( theta1/nrefrac )
    airEquivPath = tcm * math.cos(theta1-theta2)/math.cos(theta2)
       # air equivalent path is LESS than thickness of dielectric if dielectric slab is tilted !
    phs = [] 
    trans = []
    for freq in fGHz :
      ampTpar,ampTperp,ampRpar,ampRperp = pol.plateTrans( freq, nrefrac, nrefrac, angIdeg=angIdeg, tcm=tcm, tanDelta=tanDelta)
      trans.append( abs( ampTperp ) )
      phs.append( numpy.angle( ampTperp, deg=True ) - 360.*freq/clight*airEquivPath )
        # ampTperp gives phase through plate; must subtract phase through airEquivPath 
    return unwrap(numpy.array(phs) ), numpy.array(trans)

# === under construction === #
# use Born and Wolf [13.4] eqn (7) and (4) to compute transmission into lossy dielectric
# this computes complex Fresnel coefficients 
#def complexFresnel( n1, kappa1, ):
#    factor1 = sq(n2)*(1-sq(kappa2)) - sq(n1)*sq(sin(theta1))
#    factor2 = math.sqrt( sq(sq(n2)*(1-sq(kappa2) - sq(n1)*sq(sin(theta1)))) + 4.*pow(n2,4)*sq(kappa2) )
#    u = 0.5 * math.sqrt( factor1 + factor2 )
#    v = 0.5 * math.sqrt( -1.*factor1 + factor2 )
#    tpar =  
# === under construction === #

def unwrap( phi ) :
    for n in range(0,len(phi)) :
      while phi[n] < -180. : 
        phi[n] = phi[n] + 360.
      while phi[n] > 180. :
        phi[n] = phi[n] - 360.
    return phi
     
# fits refractive index of dielectric slab from phase differences (slab in - slab out)
# retrieves measurements using readTransFile() - edit that routine as necessary
# models amp and phase of transmitted signal assuming 1 pass through the slab,
#    OR multiple passes using Fabry Perot etalon formulae
# resolves 2pi ambiguities by keeping data-model phase differences between -180 and 180
# inputs: tcm = thickness of slab in cm; angIdeg = angle of incidence in degrees;
#    tanDelta = estimate of loss tangent; guess = estimate of refractive index,
#    plotType = 0 for none, 1 for phases and residuals vs freq; 2 for resids vs fit, 3 for raw data
#
def nrefracFit( infile, nguess=1.765, tanDelta=0., plotType=0 ) :
    pyplot.ioff()
    pp = PdfPages("refrac.pdf")
    fGHz,trans,dphs_meas,title,tcm,angIdeg = readTransFile( infile )
       # read in the data; edit as necessary as data format changes
    print "\n%s" % title
    print "tcm = %.5f" % tcm
    print "angIdeg = %.2f\n" % angIdeg
    nr = []
    avg = []
    std = []
    imin = 0

  # if nguess=0, find best fit refractive index in range [1,3.5]
    if abs(nguess) < 1. :
      for nrefrac in numpy.arange( 1., 3.55, .05) :
        nr.append( nrefrac )
        dphs_est,trans_est = fitFP2( fGHz, nrefrac, tcm, angIdeg, tanDelta )
        #dphs_est,trans_est = fitFP2( fGHz, nrefrac, tcm, angIdeg, tanDelta )
        dphi = unwrap( dphs_meas-dphs_est )
          # unwrap keeps phase in range -180 to 180 in all cases
        avg.append( numpy.average(dphi) )
        std.append( numpy.std(dphi) / math.sqrt(len(fGHz)) )
        if abs(std[-1]) < abs(std[imin]) :
          imin = len(avg)-1
        print ".. n = %6.3f  avg = %7.3f  std = %7.3f" % (nrefrac, avg[-1], std[-1])
      if imin == 0:
        imin = 1
      slope = (avg[imin] - avg[imin-1])/(nr[imin] - nr[imin-1])
      nrfit = nr[imin] - avg[imin]/slope  
      fitunc = abs(2.*std[imin]/slope)
      print "best fit is nrefrac = %.3f (+/- %.3f)" % (nrfit,fitunc)
    else : 
      nrfit = nguess

  # plot mean phase residual over the range [nrfit-.1, nrfit+.1]
    nr = []
    avg = []
    std = []
    imin = 0
    for nrefrac in numpy.arange( nrfit-.1,nrfit+.11,.02 ) :
      nr.append( nrefrac )
      dphs_est,trans_est = fitFP2( fGHz, nrefrac, tcm, angIdeg, tanDelta )
      dphi = unwrap( dphs_meas-dphs_est )
        # unwrap keeps phase in range -180 to 180 in all cases
      avg.append( numpy.average(dphi) )
      std.append( numpy.std(dphi) / math.sqrt(len(fGHz)) )
      if abs(std[-1]) < abs(std[imin]) :
        imin = len(avg)-1
      print ".. n = %6.3f  avg = %7.3f  std = %7.3f" % (nrefrac, avg[-1], std[-1])
    if imin == 0:
      imin = 1
    slope = (avg[imin] - avg[imin-1])/(nr[imin] - nr[imin-1])
    nrfit = nr[imin] - avg[imin]/slope  
    fitunc = abs(2.*std[imin]/slope)
    print "best fit is nrefrac = %.3f (+/- %.3f)" % (nrfit,fitunc)

  # begin the plot
    fig =  pyplot.figure( figsize=(11,8) )
    pyplot.suptitle( "%s (%s)" % (title,infile), fontsize=14 )

  # lower left plot will show mean phase resid vs refractive index
    ax = fig.add_axes( [.1,.1,.38,.2] )
    xmin,xmax = minmax( numpy.array(nr) )
    ymin,ymax = minmax( numpy.array(avg) )
    ax.axis( [xmin, xmax, ymin, ymax] )
    ax.tick_params( labelsize=10 )
    ax.errorbar( nr, avg, yerr=std, linewidth=2 )
    ax.errorbar( nrfit, 0., xerr=fitunc, color='red', elinewidth=6 )
    ax.set_xlabel("refractive index", fontsize=10)
    ax.set_ylabel("mean phase residual (deg)", fontsize=10)
    ax.text( 0.92, .82, "best fit = %.3f (+/- %.3f)" % (nrfit,fitunc), transform=ax.transAxes, \
      horizontalalignment='right', fontsize=10, color="black", rotation='horizontal' )
    pyplot.grid(True)

  # compute amp and phase residuals for best fit
    dphs_est,trans_est = fitFP2( fGHz, nrfit, tcm, angIdeg, tanDelta )
    dphi = unwrap( dphs_meas-dphs_est )
 
  # top left panel is phase vs freq, measured and fit
    ax = fig.add_axes( [.1,.6,.38,.32] )
    xmin,xmax = minmax( fGHz )
    ymin,ymax = minmax( dphs_meas )
    ax.axis( [xmin, xmax, ymin, ymax] )
    ax.tick_params( labelsize=10 )
    ax.plot( fGHz, dphs_meas, "o" )
    ax.plot( fGHz, dphs_est, "r-" )
    ax.set_ylabel("phase (deg)", fontsize=10)
    pyplot.grid(True)
  
  # middle left panel shows phase residual vs freq; avg shown as horiz dashed line
    ax = fig.add_axes( [.1,.37,.38,.2]) 
    ymin,ymax = minmax( dphi, margin=.2 )
    if abs(ymin) < abs(ymax) :
      ymin = math.copysign(ymax,ymin)
    if abs(ymax) < abs(ymin) :
      ymax = math.copysign(ymin,ymax)
    ax.axis( [xmin, xmax, ymin, ymax] )
    ax.plot( fGHz, dphi, "o-" )
    ax.tick_params( labelsize=10 )
    dphiavg = numpy.average(dphi)
    ax.plot( [fGHz[0],fGHz[-1]],[dphiavg,dphiavg],"r--" )
    ax.set_ylabel("residual (deg)", fontsize=10)
    ax.set_xlabel("freq (GHz)", fontsize=10)
    pyplot.grid(True)

  # top right panel is transmission vs freq, measured and theoretical
    ax = fig.add_axes( [.57,.6,.38,.32] )
    ymin,ymax = minmax( trans )
    ax.axis( [xmin, xmax, ymin, ymax] )
    ax.plot( fGHz, trans, "o" )
    ax.plot( fGHz, trans_est, "r-" )
    ax.tick_params( labelsize=10 )
    ax.set_ylabel("transmission", fontsize=10)
    pyplot.grid(True)

  # middle right panel is measured-theoretical trans
    ax = fig.add_axes( [.57,.37,.38,.2])
    ax.tick_params( labelsize=10 )
    ymin,ymax = minmax( trans-trans_est, margin=0.2 )
    ax.axis( [xmin, xmax, ymin, ymax] )
    ax.plot( fGHz, trans-trans_est, "o-" )
    ax.set_xlabel("freq (GHz)", fontsize=10)
    ax.set_ylabel("residual", fontsize=10)
    pyplot.grid(True)

  # fictitious lower right panel gives room for text
    ax = fig.add_axes( [.57,.1,.38,.2])
    ax.text( 0.03, .84, "n = %.3f (+/- %.4f)" % (nrfit,fitunc), transform=ax.transAxes, \
      horizontalalignment='left', fontsize=14, color="black", rotation='horizontal' )
    ax.text( 0.03, .7, "tcm = %.4f" % tcm, transform=ax.transAxes, \
      horizontalalignment='left', fontsize=14, color="black", rotation='horizontal' )
    ax.text( 0.03, .56, "angle = %.1f deg" % angIdeg, transform=ax.transAxes, \
      horizontalalignment='left', fontsize=14, color="black", rotation='horizontal' )
    ax.text( 0.03, .42, "tanDelta = %.1e" % tanDelta, transform=ax.transAxes, \
      horizontalalignment='left', fontsize=14, color="black", rotation='horizontal' )
    pyplot.axis('off')
    pyplot.savefig( pp, format="pdf" )
    pyplot.show()
    pp.close()

def minmax( varray, margin=.05 ) :
    '''find ymin and ymax for a plot'''
    vmin = varray.min()
    vmax = varray.max()
    diff = vmax - vmin
    return [ vmin - margin*diff, vmax + margin*diff ]

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
  
def compareMethods() :
  nuArr = numpy.array(numpy.arange(76.e9,90.e9,1.e9))
  print nuArr
  noArr = numpy.array([1.585]) #Actually measured
  neArr = numpy.array([1.585]) #Actually measured
  oltArr = numpy.array([1.5e-3]) #Estimated from Dick's and my measurement
  eltArr = numpy.array([1.5e-3]) #Estimated from Dick's and my measurement
  dArr = numpy.array([6.3652e-3]) #[m], Actually measured
  plateAngles = numpy.array([0.]) #[deg]
  stackRotAngle = 0. #[deg], arbitrary (but maybe one of the angles that Dick and I measured?)
  polAngle = 90. # pol parallel to plane of incidence
  stackTilt = 16.2 #[deg]

  #P trans, S trans, P refl, S refl
  pT, sT, pR, sR = hou.anisotropicTransmissionPol(nuArr, noArr, neArr, oltArr, eltArr, dArr, plateAngles, stackRotAngle, stackTilt, polAngle)
  print "pT : ", pT
  print " "
  print "sT : ", sT 

# transmission from characteristic matrix for 1 isotropic, lossless layer, E in plane of incidence
# based on Hecht eq 9.91
def fitFP2( fGHz, nrefrac, tcm, angIdeg, tanDelta ) :
    theta1 = numpy.radians(angIdeg)
    theta2 = math.asin( theta1/nrefrac )
    h = nrefrac * tcm * cos(theta2)
    Y0 = 1./cos(theta1)
    #Y0 = cos(theta1)
    Y1 = nrefrac/cos(theta2)
    phs = [] 
    trans = []
    for freq in fGHz :
      k0 = -2.*math.pi*freq/clight
      M = numpy.matrix( [[cos(k0*h), 1j*sin(k0*h)/Y1],[1j*Y1*sin(k0*h), cos(k0*h)]] )
      ampTperp = 2.*Y0/(Y0*M.item(0,0) + Y0*Y0*M.item(0,1) + M.item(1,0) + Y0*M.item(1,1)) \
               / numpy.exp(-1j * k0 * tcm * cos(theta2) )
      trans.append( abs( ampTperp ) )
      phs.append( numpy.angle( ampTperp, deg=True )) # - 360.* freq * tcm * cos(theta2) / clight )
    return unwrap(numpy.array(phs) ), numpy.array(trans)
      

