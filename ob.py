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

# readList is a helper routine for readTransFile
# it interprets an input string (whatever is after the "=" sign on an input line),
#    converting it to a list
# if input string is a const, result is [ x ]
# if input string is a triplet x,y,z, result is [ numpy.arange(x,y,z) ]
def readList( inString ) :
  try :
    a = inString.split(",")
    if len(a) == 1 :
      return numpy.array( [ float(a[0]) ] )
    elif len(a) == 2 :
      return numpy.array( [float(a[0]), float(a[1]) ] )
    elif len(a) == 3 :
      return numpy.array(numpy.arange( float(a[0]), float(a[1]), float(a[2]) ))
    else :
      print "readList: error interpreting input string '%s'" % inString
      return [0.]
  except :
      print "readList: error interpreting input string '%s'" % inString
      return [0.]

# reads *.trans file and returns arrays of freq, ampratio, deltaPhase, title, tcm, angIdeg
def readTransFile( infile ) :
  fin = open(infile, "r")
  fGHz = []
  trans = []
  dphi = []
  totcmRange = [0.,1000.]
  tcmList = []
  tanDeltaList = []
  nrList = []
  title = ""
  angIdeg = 0.
  for line in fin :
    if not line.startswith( "#" ) :
      a = line.split()
      fGHz.append( float( a[0] ) )
      trans.append( float( a[1] ) )
      dphi.append( float( a[2] ) )
    elif "totcm" in line :
      totcmRange = readList( line[line.find("=")+1:] )
    elif "totin" in line :
      totcmRange = 2.54*readList( line[line.find("=")+1:] )
    elif "angIdeg" in line :
      angIdeg = float( line[line.find("=")+1:] )
    elif "tcm" in line :
      tcmList.append( readList(line[line.find("=")+1:]) )
    elif "tin" in line :
      tcmList.append( readList(line[line.find("=")+1:]) )
      if tcmList[-1][-1] > 0. :
        tcmList[-1] = 2.54*tcmList[-1]			# to avoid screwing up mirroring
    elif "tanDelta" in line :
      tanDeltaList.append( readList(line[line.find("=")+1:]) )
    elif "nr" in line :
      nrList.append( readList(line[line.find("=")+1:]) )
    elif "title" in line :
      title = line[line.find("=")+1:]
      title = title[2:-2]
  fin.close()
  print "%d layers found" % len(nrList)
  for i in range(0,len(nrList)) :
    print "  layer %d: tcmList = " % i, tcmList[i]
    print "            nrList =", nrList[i] 
    print "      tanDeltaList = ", tanDeltaList[i]
    print " "
  return numpy.array(fGHz), numpy.array(trans), numpy.array(dphi), title, angIdeg, totcmRange, tcmList, nrList, tanDeltaList

# tcmList is a list of lists: tcmList[i] is a list of possible thicknesses for layer i
# not all combinations of thicknesses are possible, however, since total thickness 
#   of the sample is constrained to lie in the range totcm[0,1]
# when computing transmission, skip all combinations which fail this test

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
 
# in 09nov setup, E field is vertical, perpendicular to plane of incidence
def fitFP( fGHz, nrefrac, tcm, angIdeg, tanDelta, npar=False ) :
    # print "tanDelta = %.2e" % tanDelta
    theta1 = math.radians( angIdeg )
    theta2 = math.asin( math.sin(theta1)/nrefrac )
    # print "theta1,2 = ",theta1,theta2
    airEquivPath = tcm / math.cos(theta2) * math.cos(theta1-theta2)
       # air equivalent path is LESS than thickness of dielectric if dielectric slab is tilted !
    phs = [] 
    trans = []
    for freq in fGHz :
      ampTpar,ampTperp,ampRpar,ampRperp = pol.plateTrans( freq, nrefrac, nrefrac, angIdeg=angIdeg, tcm=tcm, tanDelta=tanDelta)
      # trans.append( abs( ampTpar ) )
      # phs.append( numpy.angle( ampTpar, deg=True ) - 360.*freq/clight*airEquivPath )
      trans.append( abs( ampTperp ) )
      #print "theta1,2 = ",theta1,theta2
      #print "airEquivPath = ",airEquivPath," cm"
      #print "air equiv phase = ", 360.*freq/clight * airEquivPath
      phs.append( numpy.angle( ampTperp, deg=True ) - 360.*freq/clight*airEquivPath )
        # ampTperp gives phase through plate; must subtract phase through airEquivPath 
    return unwrap(numpy.array(phs) ), numpy.array(trans)

# use Born and Wolf [13.4] eqn (7) and (4) to compute n2*cos(theta2) for lossy dielectic
def complexNC( n1, theta1Radians, n2, tanDelta ):
    if tanDelta == 0. :
      return n2 * math.cos( math.asin(n1*math.sin(theta1Radians)/n2) )
    else :
      kappa2 = 0.5*tanDelta
      factor1 = sq(n2)*(1-sq(kappa2)) - sq(n1)*sq(sin(theta1Radians))
      factor2 = numpy.sqrt( sq(sq(n2)*(1-sq(kappa2)) - sq(n1)*sq(sin(theta1Radians))) + 4.*pow(n2,4)*sq(kappa2) )
      u2 = numpy.sqrt( 0.5 * (factor1 + factor2) )
      v2 = numpy.sqrt( 0.5 * (-1.*factor1 + factor2) )
      return u2 + 1j*v2

def unwrap( phi ) :
    for n in range(0,len(phi)) :
      while phi[n] < -180. : 
        phi[n] = phi[n] + 360.
      while phi[n] > 180. :
        phi[n] = phi[n] - 360.
    return phi
     
# fits refractive index of 1 layer in stack of dielectric slabs from phase differences (slab in - slab out)
# retrieves measurements and stack parameters using readTransFile() - edit that routine as necessary
# models amp and phase of transmitted signal 
# resolves 2pi ambiguities by keeping data-model phase differences between -180 and 180
# inputs: tcm = thickness of slab in cm; angIdeg = angle of incidence in degrees;
#    tanDelta = estimate of loss tangent; guess = estimate of refractive index,
#
# input lists: tcmList = thicknesses of layers in cm
#              nrList = refractive indices of layers; 1st layer may have nr = 0, which triggers 
#                         search for best fit
#              tanDeltaList = tan deltas of layers (not used yet)
#

def nrefracFit( infile, nrSrchList=numpy.arange(1.3,2.5,.05), tcmSrchList=numpy.arange(0.03696,.06,1.), tanDelta=0., angIdegOver=0., tcmOver=0 ) :
    pyplot.ioff()
    pp = PdfPages("refrac.pdf")
    print " "
    fGHz,trans,dphs_meas,title,angIdeg,tcmTotal,tcmList,nrList,tanDeltaList = readTransFile( infile )
    # if angIdegOver != 0. :
    #   angIdeg = angIdegOver
    # if tcmOver != 0. :
    #   tcm = tcmOver
    print "\n%s" % title
    print "angIdeg = %.2f" % angIdeg
    print "total thickness cm\n" % tcmTotal

    nr = []
    tcm = []
    avg = []
    std = []
    imin = 0
    nrsave = 0.
    fitunc = 0.

  # if nr=0 for 1st layer, search for best fit refractive index and layer thickness in ranges nrSrchList and tcmSrchList
    if abs(nrList[0]) < 1. :
      for nrefrac in nrSrchList :
        for tcm0 in tcmSrchList :
          nr.append( nrefrac )
          tcm.append( tcm0 )

        # enter these values into the parameters for layer 0
          nrList[0] = nrefrac
          tcmList[0] = tcm0

        # adjust layer 1 thickness to keep total thickness constant
          tcmList[1] = tcmTotal - tcmList[0]

        # model this stack, compute phase and amp differences
          dphs_est,trans_est = solveStack( fGHz, angIdeg, tcmList, nrList, tanDeltaList ) 
          dphi = unwrap( dphs_meas-dphs_est )
            # unwrap keeps phase in range -180 to 180 in all cases
          avg.append( numpy.average(dphi) )
          std.append( numpy.std(dphi) / math.sqrt(len(fGHz)) )
          avgT = numpy.average( trans - trans_est )
          stdT = numpy.std( trans - trans_est )
        
          if abs(avg[-1]) < abs(avg[imin]) :
            imin = len(avg)-1
          print ".. n = %6.3f  tcm0 = %.4f avg = %7.3f  std = %6.3f  avgT = %6.3f  stdT = %6.3f" % (nrefrac, tcm0, avg[-1], std[-1], avgT, stdT)

      if imin == 0:
        imin = 1
      slope = (avg[imin] - avg[imin-1])/(nr[imin] - nr[imin-1])
      nrfit = nr[imin] - avg[imin]/slope  
      fitunc = abs(2.*std[imin]/slope)
      print "best fit is nrefrac = %.3f (+/- %.3f)" % (nrfit,fitunc)

      #nrfit = nr[imin]
      tcmfit = tcm[imin]
      #print "best fit is nr0 = %.3f, tcm0 = %.4f" % (nrfit,tcmfit)
    
    # at this point, fix tcm values for both layers
      tcmList[0] = tcm0
      tcmList[1] = tcmTotal - tcmList[0]

  # if nrList[0] was non-zero, set nrfit to this value; leave input thicknesses unchanged
    else : 
      nrfit = nrList[0]
      nrsave = nrList[0]

  # now search narrower range for best nr, keeping thicknesses constant
  # always plot mean phase residual over the range [nrfit-.1, nrfit+.1]

    nr = []
    avg = []
    std = []
    imin = 0
    for nrefrac in numpy.arange( nrfit-.1,nrfit+.11,.02 ) :
      nr.append( nrefrac )
      nrList[0] = nrefrac
      #dphs_est,trans_est = fitFP2( fGHz, nrefrac, tcm, angIdeg, tanDelta )
      # dphs_est,trans_est = solveStack( fGHz, angIdeg, [tcm], [nrefrac], [tanDelta] ) 
      dphs_est,trans_est = solveStack( fGHz, angIdeg, tcmList, nrList, tanDeltaList ) 
      dphi = unwrap( dphs_meas-dphs_est )
        # unwrap keeps phase in range -180 to 180 in all cases
      avg.append( numpy.average(dphi) )
      std.append( numpy.std(dphi) / math.sqrt(len(fGHz)) )
      avgT = numpy.average( trans - trans_est )
      stdT = numpy.std( trans - trans_est )
      if abs(avg[-1]) < abs(avg[imin]) :
        imin = len(avg)-1
      print ".. n = %6.3f  avg = %7.3f  std = %7.3f  avgT = %7.3f  stdT = %7.3f" % (nrefrac, avg[-1], std[-1], avgT,stdT)
    if imin == 0:
      imin = 1
    slope = (avg[imin] - avg[imin-1])/(nr[imin] - nr[imin-1])
    nrfit = nr[imin] - std[imin]/slope  
    fitunc = abs(2.*std[imin]/slope)
    #nrfit = nr[imin]
    #print "best fit is nrefrac = %.3f (+/- %.3f)" % (nrfit,fitunc)

  # in case nrList[0] was given, restore it so that final plots show this value (even if it was not the best fit)
    nrList[0] = nrfit
  #  if nrsave > 0. : nrList[0] = nrsave

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
    # dphs_est,trans_est = fitFP2( fGHz, nrfit, tcm, angIdeg, tanDelta )
    # dphs_est,trans_est = solveStack( fGHz, angIdeg, [tcm], [nrfit], [tanDelta] ) 
    dphs_est,trans_est = solveStack( fGHz, angIdeg, tcmList, nrList, tanDeltaList ) 
    dphi = unwrap( dphs_meas-dphs_est )
 
  # top left panel is phase vs freq, measured and fit
    ax = fig.add_axes( [.1,.6,.38,.32] )
    xmin,xmax = minmax( fGHz )
    ymin,ymax = minmax( numpy.concatenate((dphs_meas,dphs_est)))
    ax.axis( [xmin, xmax, ymin, ymax] )
    ax.tick_params( labelsize=10 )
    ax.plot( fGHz, dphs_meas, "o" )
    ax.plot( fGHz, dphs_est, "r-" )
    ax.set_ylabel("phase (deg)", fontsize=10)
    pyplot.grid(True)
  
  # middle left panel shows phase residual vs freq; avg shown as horiz dashed line
    ax = fig.add_axes( [.1,.37,.38,.2]) 
    ymin,ymax = minmax( dphi, margin=.2 )
    #if abs(ymin) < abs(ymax) :
    #  ymin = math.copysign(ymax,ymin)
    #elif abs(ymax) < abs(ymin) :
    #  ymax = math.copysign(ymin,ymax)
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
    ymin,ymax = minmax( numpy.concatenate((trans,trans_est)) )
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
    if fitunc == 0. :
      ax.text( 0.03, .84, "n = %.3f" % nrsave, transform=ax.transAxes, \
        horizontalalignment='left', fontsize=14, color="black", rotation='horizontal' )
    else :
      ax.text( 0.03, .84, "n = %.3f (+/- %.4f)" % (nrfit,fitunc), transform=ax.transAxes, \
        horizontalalignment='left', fontsize=14, color="black", rotation='horizontal' )
    ax.text( 0.03, .7, "tcm = %.4f" % tcmList[0], transform=ax.transAxes, \
      horizontalalignment='left', fontsize=14, color="black", rotation='horizontal' )
    ax.text( 0.03, .56, "angle = %.1f deg" % angIdeg, transform=ax.transAxes, \
      horizontalalignment='left', fontsize=14, color="black", rotation='horizontal' )
    ax.text( 0.03, .42, "tanDelta = %.1e" % tanDeltaList[0], transform=ax.transAxes, \
      horizontalalignment='left', fontsize=14, color="black", rotation='horizontal' )
    pyplot.axis('off')
    pyplot.savefig( pp, format="pdf" )
    pyplot.show()
    pp.close()


# expandTree uses meshgrid to expand the thickness and refractive index arrays
#
# input: 
#   tcmRange = [tcmMin,tcmMax]  - specifies allowed range of total thicknesses
#   tcmListIn = [ [tcm1list], [tcm2list], ...] - lists of thicknesses to model for each layer
#   nrListIn = [ [nr1list], [nr2list], ... ] - lists of refractive indices to model for each layer
# output :
#   tcmListOut = [ [tcm1,tcm2,...], [tcm1,tcm2,...], ... ] - expanded lists of thicknesses
#   nrListOut = [ [nr1,nr2,...], [nr1,nr2,...], ... ] - expanded lists of refractive indices
#
# thus: tcmlistOut[i] and nrListOut[i] specify one set of thicknesses and refractive indices
#   that should be modeled; note that same list of thicknesses may appear many times, one
#   for each possible set of refractive indices; same for refractive indices
# total thickness of tcmListOut is guaranteed to be in the range [tcmMin,tcmMax]
 
def expandTree( tcmRange, tcmListIn, nrListIn) :

    nlayers = len(tcmListIn)
    #print "tcmListIn = ",tcmListIn
    #print "nlayers = ", nlayers
    #print "tcmRange = ", tcmRange

  # meshgrid arrays of thicknesses AND refractive indices
  # this is lame, but I can't figure out how to do this elegantly, without if structure
    if nlayers == 1 :
      a = numpy.meshgrid( tcmListIn[0], nrListIn[0] )
    elif nlayers == 2 :
      a = numpy.meshgrid( tcmListIn[0],tcmListIn[1],nrListIn[0],nrListIn[1] ) 
    elif nlayers == 3 :
      a = numpy.meshgrid( tcmListIn[0],tcmListIn[1],tcmListIn[2], \
            nrListIn[0],nrListIn[1],nrListIn[2])
    elif nlayers == 4 :
      a = numpy.meshgrid( tcmListIn[0],tcmListIn[1],tcmListIn[2],tcmListIn[3],\
            nrListIn[0],nrListIn[1],nrListIn[2],nrListIn[3])
    elif nlayers == 5 :
      a = numpy.meshgrid( tcmListIn[0],tcmListIn[1],tcmListIn[2],tcmListIn[3],tcmListIn[4], \
            nrListIn[0],nrListIn[1],nrListIn[2],nrListIn[3],nrListIn[4] )
    else :
      print "nlayers = %d not allowed by expandTree" % nlayers
      return

  # flatten each array contained in a 
    b = []
    for i in range(0,len(a)) :
      b.append( a[i].flatten() )

  # constrained layers are indicated by -L, where L is the corresponding layer
    for n in range(0,nlayers) :
      if b[n][0] < 0. :
        m = int( abs(b[n][0]) + 0.5 )
        try :
          b[n] = b[m-1]
          print "mirrored row %d into row %d" % (m-1,n)
        except :
          print "mirroring failed; tried to copy row %d into row %d" % (m-1,n)
          return
      if b[n+nlayers][0] < 0. :
        m = int( abs(b[n+nlayers][0]) + 0.5 ) + nlayers
        try :
          b[n+nlayers] = b[m-1]
          print "mirrored row %d into row %d" % (m-1,n+nlayers)
        except :
          print "mirroring failed; tried to copy row %d into row %d" % (m-1,n+nlayers)
          return

  # compute total thickness for each possible stack
    totcm = numpy.zeros( len(b[0]) )
    for i in range(0,nlayers) :
      totcm = totcm + b[i]
    #print "totcm = ",totcm
    #print " ========== "
  
  # now form transpose, retain only rows where thickness is within allowed range
    c = numpy.transpose(b)
    d = []
    for n in range(0,len(c)) :
      if (totcm[n] >= tcmRange[0]) and (totcm[n] <= tcmRange[1]) :
        d.append( c[n] )
        #print d[-1]    

  # return separate tcm and nr arrays
    return numpy.hsplit( numpy.array(d),2)


def nrefracFit2( infile ) :
    pyplot.ioff()
    pp = PdfPages("refrac.pdf")
    print " "
    fGHz,trans,dphs_meas,title,angIdeg,tcmRange,tcmListIn,nrListIn,tanDeltaList = readTransFile( infile )
    tcmList,nrList = expandTree( tcmRange, tcmListIn, nrListIn )
    print "\n%s" % title
    print "angIdeg = %.2f" % angIdeg
    print "total thickness in range ", tcmRange

  # model all stacks in list, save avg and std phase and amp residuals
    nlayers = len( tcmList[0] )
    tanDelta = []
    for n in range(0,nlayers) :
      tanDelta.append( 0. )
    ilen = len( tcmList )
    print "modeling %d layers, %d possible stacks" % (nlayers, ilen)
    avgPhResid = 1.e6*numpy.ones( ilen )
    stdPhResid = 1.e6*numpy.ones( ilen )
    avgAmpResid = 1.e6*numpy.ones( ilen )
    stdAmpResid = 1.e6*numpy.ones( ilen )
    for i in range(0,ilen) :
      if i%10 == 0 :
        print ".. finished %d/%d stacks" % (i,ilen)
      #print "i, tcmList, nrList: %d" % i, tcmList[i], nrList[i]
      dphs_est,trans_est = solveStack( fGHz, angIdeg, tcmList[i], nrList[i], tanDelta ) 
      dphi = unwrap( dphs_meas-dphs_est )
        # unwrap keeps phase in range -180 to 180 in all cases
      avgPhResid[i] = numpy.average(dphi) 
      stdPhResid[i] = numpy.std(dphi) / math.sqrt(len(fGHz)) 
      avgAmpResid[i] = numpy.average( trans - trans_est )
      stdAmpResid[i] = numpy.std( trans - trans_est )
      print tcmList[i], nrList[i], ".. %7.3f, %7.3f" \
             % (avgPhResid[i], stdPhResid[i])

    #merit = stdPhResid + stdAmpResid
    merit = avgPhResid
    nbest = numpy.argmin(merit)
    print "\nbest fit:"
    for nl in range(0,nlayers) :
      print " layer %d  tin = %.5f, nr = %.3f" % (nl, tcmList[nbest][nl]/2.54, nrList[nbest][nl] ) 
      
  # begin the plot
    fig =  pyplot.figure( figsize=(11,8) )
    pyplot.suptitle( "%s (%s)" % (title,infile), fontsize=14 )

  # lower left plot will show mean phase resid vs refractive index
  #  ax = fig.add_axes( [.1,.1,.38,.2] )
  #  xmin,xmax = minmax( numpy.array(nr) )
  #  ymin,ymax = minmax( numpy.array(avg) )
  #  ax.axis( [xmin, xmax, ymin, ymax] )
  #  ax.tick_params( labelsize=10 )
  #  ax.errorbar( nr, avg, yerr=std, linewidth=2 )
  #  ax.errorbar( nrfit, 0., xerr=fitunc, color='red', elinewidth=6 )
  #  ax.set_xlabel("refractive index", fontsize=10)
  #  ax.set_ylabel("mean phase residual (deg)", fontsize=10)
  #  ax.text( 0.92, .82, "best fit = %.3f (+/- %.3f)" % (nrfit,fitunc), transform=ax.transAxes, \
  #    horizontalalignment='right', fontsize=10, color="black", rotation='horizontal' )
  #  pyplot.grid(True)

  # compute amp and phase residuals for best fit
    dphs_est,trans_est = solveStack( fGHz, angIdeg, tcmList[nbest], nrList[nbest], tanDelta ) 
    dphi = unwrap( dphs_meas-dphs_est )
 
  # top left panel is phase vs freq, measured and fit
    ax = fig.add_axes( [.1,.6,.38,.32] )
    xmin,xmax = minmax( fGHz )
    ymin,ymax = minmax( numpy.concatenate((dphs_meas,dphs_est)))
    ax.axis( [xmin, xmax, ymin, ymax] )
    ax.tick_params( labelsize=10 )
    ax.plot( fGHz, dphs_meas, "o" )
    ax.plot( fGHz, dphs_est, "r-" )
    ax.set_ylabel("phase (deg)", fontsize=10)
    pyplot.grid(True)
  
  # middle left panel shows phase residual vs freq; avg shown as horiz dashed line
    ax = fig.add_axes( [.1,.37,.38,.2]) 
    ymin,ymax = minmax( dphi, margin=.2 )
    #if abs(ymin) < abs(ymax) :
    #  ymin = math.copysign(ymax,ymin)
    #elif abs(ymax) < abs(ymin) :
    #  ymax = math.copysign(ymin,ymax)
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
    ymin,ymax = minmax( numpy.concatenate((trans,trans_est)) )
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
    ylab = 0.84
    for nl in range(0,nlayers) :
      ax.text( 0.03, ylab, "layer %d:  t = %.5f in  nr = %.3f" % \
        ( nl+1, tcmList[nbest][nl]/2.54, nrList[nbest][nl]), \
        transform=ax.transAxes, horizontalalignment='left', fontsize=14, \
        color="black", rotation='horizontal' )
      ylab = ylab - 0.14
    #ax.text( 0.03, .56, "angle = %.1f deg" % angIdeg, transform=ax.transAxes, \
    #  horizontalalignment='left', fontsize=14, color="black", rotation='horizontal' )
    #ax.text( 0.03, .42, "tanDelta = %.1e" % tanDelta[0], transform=ax.transAxes, \
    #  horizontalalignment='left', fontsize=14, color="black", rotation='horizontal' )
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

# transmission from characteristic matrix for 1 isotropic, lossless layer, E perpendicular to plane of incidence
# based on Hecht eq 9.91
def fitFP2( fGHz, nrefrac, tcm, angIdeg, tanDelta ) :
    theta1 = math.radians(angIdeg)
    nCosTheta2 = complexNC( 1., theta1, nrefrac, tanDelta )
    #theta2 = math.asin( math.sin(theta1)/nrefrac )
    #h = nrefrac * tcm * cos(theta2)
    h = tcm * nCosTheta2
    Y0 = cos(theta1)
    #Y1 = nrefrac * cos(theta2)  # for E perpendicular to plane of incidence
    Y1 = nCosTheta2
    phs = [] 
    trans = []
    for freq in fGHz :
      k0 = -2.*math.pi*freq/clight
      M = numpy.matrix( [[cos(k0*h), 1j*sin(k0*h)/Y1],[1j*Y1*sin(k0*h), cos(k0*h)]] )
      ampTperp = 2.*Y0/(Y0*M.item(0,0) + Y0*Y0*M.item(0,1) + M.item(1,0) + Y0*M.item(1,1)) 
              #  / numpy.exp(-1j * k0 * tcm * cos(theta1) )   alternate way of correcting for air path
      trans.append( abs( ampTperp ) )
      phs.append( numpy.angle( ampTperp, deg=True ) - 360.* freq * tcm * cos(theta1) / clight )
          # 2nd term is the phase through air 
    return unwrap(numpy.array(phs) ), numpy.array(trans)
      
# computes theoretical transmission for arbitrary stack of ISOTROPIC layers
# assumes E is perpendicular to plane of incidence ("s wave")
# fGHz is array of frequencies
# angIdeg is angle of incidence to stack, in degrees
# tcm is list of thicknesses
# nrefrac is list of refractive indices
# tanDelta is list of loss tangents

def solveStack( fGHzArr, angIdeg, tcmArr, nrArr, tanDeltaArr ) :
    theta1 = math.radians(angIdeg)
    Y0 = cos(theta1)
    phs = [] 
    trans = []
    for freq in fGHzArr :
      M = numpy.identity(2)
      k0 = -2.*math.pi*freq/clight
    # compute transfer matrix for the stack; keep track of total thickness
      tcmtotal = 0.
      for tcm,nrefrac,tanDelta in zip( tcmArr,nrArr,tanDeltaArr ) :
        # theta2 = math.asin( math.sin(theta1)/nrefrac )
        # h = nrefrac * tcm * cos(theta2)
        # Y1 = nrefrac * cos(theta2)  # for E perpendicular to plane of incidence
        nCosTheta2 = complexNC( 1., theta1, nrefrac, tanDelta )
        h = tcm * nCosTheta2
        Y1 = nCosTheta2
        M = numpy.dot( numpy.matrix( [[cos(k0*h), 1j*sin(k0*h)/Y1],[1j*Y1*sin(k0*h), cos(k0*h)]] ), M)
        tcmtotal = tcmtotal + tcm
    # solve for ampTperp (complex); then compute transmission and phase change relative to air 
      ampTperp = 2.*Y0/(Y0*M.item(0,0) + Y0*Y0*M.item(0,1) + M.item(1,0) + Y0*M.item(1,1)) 
              #  / numpy.exp(-1j * k0 * tcm * cos(theta1) )   alternate way of correcting for air path
      trans.append( abs( ampTperp ) )
      phs.append( numpy.angle( ampTperp, deg=True ) - 360.* freq * tcmtotal * cos(theta1) / clight )
    return unwrap(numpy.array(phs) ), numpy.array(trans)
      
#angIdegList = numpy.arange( 15.,17.,.2)
#tcmList = [ [1.], 
#def fitStack( fGHzArr, angIdegList, tcmList, nrList, tanDeltaMat ):
#  nlayers = len(tcmList)
  # check that all lists have same number of layers
  # marginalized probability... 

#def fit2Stack( fGHzArr, angIdeg, tcmList, nrList, tanDelta) :
#  for tcm in tcmList :
#    for nr in nrList :
#      fitStack (blah)
#      compute rms
#      find best (if more than 1)
#      produce plot

def compFP( nrefrac=3.11, tcm=0.6, angIdeg=10., tanDelta=0. ) :
    fGHz = numpy.arange(80.,115.01,.1)
    #dphs1,trans1 = fitFP( fGHz, nrefrac, tcm, angIdeg, tanDelta )
    dphs1,trans1 = fitStack( fGHz, angIdeg, [.33,tcm],[1.8,nrefrac],[0.,0.] )
    dphs2,trans2 = fitStack( fGHz, angIdeg, [.33,tcm],[1.85,nrefrac],[0.,0.] )
    #dphs3,trans3 = fitFP2( fGHz, 1., tcm, angIdeg, tanDelta )
    fig =  pyplot.figure()
    ax = pyplot.subplot(2,1,1)
    ax.plot( fGHz, trans1, "b-" )
    ax.plot( fGHz, trans2, "r-" )
    ax = pyplot.subplot(2,1,2)
    ax.plot( fGHz, dphs1, "b-" )
    ax.plot( fGHz, dphs2, "r-" )
    #ax.plot( fGHz, dphs2-dphs3, "r-" )
    pyplot.show() 
    
