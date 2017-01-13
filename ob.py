# ob.py
# these are routines for dealing with the Optics Bench scans of various kinds

import numpy
import math
import sys
import string
import matplotlib
import pickle
matplotlib.use('GTKAgg')
import matplotlib.pyplot as pyplot
from matplotlib.backends.backend_pdf import PdfPages
import pol
#import scipy
from scipy.interpolate import splrep, splev
#import hou

# lambda is alternate way of defining simple function
sin = lambda arg: numpy.sin(arg)
cos = lambda arg: numpy.cos(arg)
inv = lambda Mat: numpy.linalg.inv(Mat)
dr = lambda theta: math.pi*theta/180.
pow = lambda arg, index: numpy.power(arg, index)
sq = lambda arg: pow(arg, 2)


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

# readTransFile inputs data from transmission file created by OpticsBench/09nov2016/unpack.py
# lines not beginning with "#" should contain fGHz, transmission, and phase change when
#    dielectric is inserted into beam
# lines beginning with "#" may contain information related to fitting, such as the angle
#    of incidence, the total thickness limits of the dielectric stack, title for a plot,
#    and range of thickness, refractive index, and tanDelta to search for each layer;
#    all these variables are indicated by keywords followed by "=" 
# lines beginning with "##" are ignored
# tcmList is a list of lists: tcmList[i] is a list of possible thicknesses for layer i
# not all combinations of thicknesses are possible, however, since total thickness 
#   of the sample is constrained to lie in the range totcm[0]-totcm[1]
# when computing transmission, skip all combinations which fail this test
#
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
    if line.startswith( "##" ) :     # ignore lines beginning with double hash mark
      continue
    elif not line.startswith( "#" ) :    # lines without hash marks are data
      a = line.split()
      fGHz.append( float( a[0] ) )
      trans.append( float( a[1] ) )
      dphi.append( float( a[2] ) )
    else :                               # lines with single hash marks are search parameters
      if "totcm" in line :
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
  print "%d layers found\n" % len(nrList)
  for i in range(0,len(nrList)) :
    print "  layer %d: tcmList = " % i, tcmList[i]
    print "            nrList = ", nrList[i] 
    print "          tanDList = ", tanDeltaList[i]
    print " "
  return numpy.array(fGHz), numpy.array(trans), numpy.array(dphi), title, angIdeg, totcmRange, tcmList, nrList, tanDeltaList


# returns expected transmission, deltaPhase for 1pass through plate (no FP reflections)
# fGHz = freq in GHz
# nrefrac = estimated refractive index
# tcm = puck thickness in cm
# angIdeg = angle of incidence, degrees
# tanDelta = estimated loss tangent
# npar = True if incident radiation is polarized parallel to plane of incidence, False if perp
#
def fit1pass( fGHz, nrefrac, tcm, angIdeg, tanDelta, npar=True ) :
    thetaI = numpy.radians(angIdeg)
    thetaT = math.asin((math.sin(numpy.radians(angIdeg)))/nrefrac)
    tpar1, tperp1, rpar, rperp = pol.fresnel( 1., thetaI, nrefrac, thetaT )
    tpar2, tperp2, rpar, rperp = pol.fresnel( nrefrac, thetaT, 1., thetaI )
    dphs_est = 360.*fGHz/clight * tcm/math.cos(thetaT) * (nrefrac - 1.)
    return dphs_est
 
# fitFP uses pol.plateTrans to compute transmission through a dielectric slab; this is the method
#   where I follow the wave rattling back and forth through the slab; it is useful only for
#   a single layer
# in 09nov setup, E field is vertical, perpendicular to plane of incidence
def fitFP( fGHz, nrefrac, tcm, angIdeg, tanDelta ) :
    theta1 = math.radians( angIdeg )
    theta2 = math.asin( math.sin(theta1)/nrefrac )
    airEquivPath = tcm / math.cos(theta2) * math.cos(theta1-theta2)
       # air equivalent path is LESS than thickness of dielectric if dielectric slab is tilted !
    phs = [] 
    trans = []
    for freq in fGHz :
      ampTpar,ampTperp,ampRpar,ampRperp = pol.plateTrans( freq, nrefrac, nrefrac, angIdeg=angIdeg, tcm=tcm, tanDelta=tanDelta)
      trans.append( abs( ampTperp ) )
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
# 
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
          print "mirrored thicknesses for layer %d into layer row %d" % (m,n+1)
        except :
          print "mirroring failed; tried to copy row %d into row %d" % (m-1,n)
          return
      if b[n+nlayers][0] < 0. :
        m = int( abs(b[n+nlayers][0]) + 0.5 )    # e.g., -1 -> 1
        try :
          b[n+nlayers] = b[m+nlayers-1]
          print "mirrored refractive indices for layer %d into layer %d" % (m,n+1)
        except :
          print "mirroring failed; tried to copy row %d into row %d" % (m+nlayers-1,n+nlayers)
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

# prints parameters for ndump stacks with smallest absolute deviations
def printBest( merit, tcmList, nrList, ndump=10 ) :
    ntrials = len( tcmList )
    nlayers = len( tcmList[0] )
    if ndump > ntrials :
      ndump = ntrials
    merit2 = numpy.abs(merit)
    nbest = merit2.argsort()[:ndump]    # finds indexes of ndump smallest values
    print nbest
    for nb in nbest:
      print tcmList[nb]/2.54, nrList[nb], merit[nb]
    return nbest[0]

def nrefracFit2( infile ) :
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
    print "modeling %d layers, %d trials" % (nlayers, ilen)
    avgPhResid = 1.e6*numpy.ones( ilen )
    stdPhResid = 1.e6*numpy.ones( ilen )
    avgAmpResid = 1.e6*numpy.ones( ilen )
    stdAmpResid = 1.e6*numpy.ones( ilen )
    for i in range(0,ilen) :
      if i%100 == 0 :
        print ".. finished %d/%d trials" % (i,ilen)
      #print "i, tcmList, nrList: %d" % i, tcmList[i], nrList[i]
      dphs_est,trans_est = solveStack( fGHz, angIdeg, tcmList[i], nrList[i], tanDelta ) 
      dphi = unwrap( dphs_meas-dphs_est )
        # unwrap keeps phase in range -180 to 180 in all cases
      avgPhResid[i] = numpy.average(dphi) 
      stdPhResid[i] = numpy.std(dphi) / math.sqrt(len(fGHz)) 
      avgAmpResid[i] = numpy.average( trans - trans_est )
      stdAmpResid[i] = numpy.std( trans - trans_est )
      #print tcmList[i], nrList[i], ".. %7.3f, %7.3f" \
      #       % (avgPhResid[i], stdPhResid[i])

    #merit = stdPhResid + stdAmpResid
   
  # save the results for further analysis
    fout = open( infile+"Fits.pickle", "ab" )
    pickle.dump( [tcmList,nrList,avgPhResid,stdPhResid,avgAmpResid,stdAmpResid], fout ) 
    fout.close 

  # save best fits in file
    # fout = open( infile+".fit" )
    print "\n10 best fits minimizing std residual amps:"
    # fout.write("\n10 best fits minimizing std residual amps:\n")
    printBest( stdAmpResid, tcmList, nrList ) 

    print "\n10 best fits minimizing avg residual phases:"
    # fout.write("\n10 best fits minimizing avg residual phases:\n")
    printBest( avgPhResid, tcmList, nrList ) 

    print "\n10 best fits minimizing std residual phases:"
    # fout.write("\n10 best fits minimizing std residual phases:\n")
    nbest = printBest( stdPhResid, tcmList, nrList ) 
    # fout.close()

  # begin the plot
    pyplot.ioff()
    pp = PdfPages("refrac.pdf")
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
# note: all data in OpticsBench/09nov2016 were taken with E field perpendicular to plane of incidence (s-polarized)
def fitFP2( fGHz, nrefrac, tcm, angIdeg, tanDelta ) :
    theta1 = math.radians(angIdeg)
    Y0 = cos(theta1)
    #theta2 = math.asin( math.sin(theta1)/nrefrac )
    #h = nrefrac * tcm * cos(theta2)
    #Y1 = nrefrac * cos(theta2)  # for E perpendicular to plane of incidence
    nCosTheta2 = complexNC( 1., theta1, nrefrac, tanDelta )
    h = tcm * nCosTheta2
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
      
# compare transmission and phase predicted by different subroutines
# note that all cases are for S-polarized radiation (E field perpendicular to plane of incidence)
def compFP( nrefrac=3.11, tcm=0.6, angIdeg=10., tanDelta=1.e-3 ) :
    theta1 = numpy.radians(angIdeg) 
    fGHz = numpy.arange(80.,115.01,0.5)
    dphs1,trans1 = fitFP( fGHz, nrefrac, tcm, angIdeg, tanDelta )
    dphs2,trans2 = fitFP2( fGHz, nrefrac, tcm, angIdeg, tanDelta ) 
    dphs3,trans3 = solveStack( fGHz, angIdeg, [tcm], [nrefrac], [tanDelta] )
    # pT,sT,pR,sR = anisotropicTransmissionPol(fGHz, [nrefrac], [nrefrac], [tanDelta], [tanDelta], \
    sT = anisotropicTransmissionPol(fGHz, [nrefrac], [nrefrac], [tanDelta], [tanDelta], \
       [tcm], [0.], rho=0., incAngle=angIdeg, polAngle=90.)
    dphs4 = unwrap(-numpy.angle( sT, deg=True ) - 360.* fGHz * tcm * cos(theta1) / clight )
    trans4 = numpy.abs(sT)
    fig =  pyplot.figure()
    ax = pyplot.subplot(2,1,1)
    ax.plot( fGHz, trans1, "b-" )
    ax.plot( fGHz, trans2, "ro" )
    ax.plot( fGHz, trans4, "c-" )
    ax = pyplot.subplot(2,1,2)
    ax.plot( fGHz, dphs1, "b-" )
    ax.plot( fGHz, dphs2, "ro" )
    ax.plot( fGHz, dphs4, "c-" )
    #ax.plot( fGHz, dphs2-dphs3, "r-" )
    pyplot.show() 
    

# ********** Transmission through a stack of birefringent slabs (Tom's Code) **********

#Transfer function for stratified medium
#This function should only be run within the context of the function "anisotropicTransmissionPol"
# changes to Tom's code:
#   - input angles are degrees, but internally everything is radians; lambda functions cos,sin
#        were defined accordingly at the top
#   - nu is now in GHz, not Hz
#   - t is now in cm, not meters
#
def TransMat(nu, noArr, neArr, oLTArr, eLTArr, tArr, chiArr, rho, incAngle, incN):
  # First medium is air
    theta1 = numpy.radians(incAngle)    # internal to this routine, all angles are radians
    n1 = incN
    
  # Create transfer matrix
    Tr = numpy.identity(4)

  # Loop over the layers
    for i in range(len(chiArr)):
        no = noArr[i]
        ne = neArr[i]
        oLT = oLTArr[i]
        eLT = eLTArr[i]
        t = tArr[i]
        chi = numpy.radians( chiArr[i]+rho )   #Plate orientation w.r.t x-axis

      # Vacuum wavevector
        k0 = (2*math.pi)/(clight/nu)
        
        nP = no
        nPP = ne*numpy.sqrt(1+(pow(ne,-2)-pow(no,-2))*sq(n1)*sq(sin(theta1))*sq(cos(chi)))
        
        R = lambda chi: numpy.matrix([[cos(chi),-sin(chi),0],[sin(chi),cos(chi),0],[0,0,1]])
        epsilonP = R(chi)*numpy.matrix([[sq(ne),0,0],[0,sq(no),0],[0,0,sq(no)]])*R(-chi)
        
        LTp = oLT
        nPT = nP*numpy.sqrt(1-1j*LTp)
     
      # I added this to handle isotropic dielectrics
        if ne != no :
          LTpp = oLT + ((eLT - oLT)/(ne - no))*(nPP - no)
        else :
          LTpp = oLT

        nPPT = nPP*numpy.sqrt(1-1j*LTpp)
        
        thetaP = numpy.arcsin(n1*sin(theta1)/nP)
        thetaPP = numpy.arcsin(n1*sin(theta1)/nPP)
        
        deltaP = nPT*t*cos(thetaP)
        deltaPP = nPPT*t*cos(thetaPP)
        
        DPt1 = pow(sq(cos(thetaP))+sq(sin(thetaP))*sq(sin(chi)),-0.5) * \
           numpy.array([-sin(chi)*cos(thetaP),cos(chi)*cos(thetaP),sin(chi)*sin(thetaP)])
        DPPt1 = pow(sq(cos(chi))*sq(cos(thetaP))+sq(sin(chi))*sq(cos(thetaP-thetaPP)),-0.5) * \
           numpy.array([cos(chi)*cos(thetaP)*cos(thetaPP),sin(chi)*(sin(thetaP)*sin(thetaPP)+cos(thetaP)*cos(thetaPP)), \
                  -cos(chi)*cos(thetaP)*sin(thetaPP)])
        HPt1 = pow(sq(cos(thetaP))*sq(cos(chi))+sq(sin(chi)),-0.5) * \
           numpy.array([-sq(cos(thetaP))*cos(chi),-sin(chi),cos(thetaP)*sin(thetaP)*cos(chi)])
        HPPt1 = pow(sq(cos(thetaP-thetaPP))*sq(sin(chi))+sq(cos(thetaP))*sq(cos(chi)),-0.5) * \
           numpy.array([-cos(thetaP-thetaPP)*cos(thetaPP)*sin(chi),cos(thetaP)*cos(chi),cos(thetaP-thetaPP)*sin(thetaPP)*sin(chi)])
        Phi1 = numpy.matrix([[DPt1.item(0),DPPt1.item(0),DPt1.item(0),DPPt1.item(0)], \
          [(1./nP)*HPt1.item(1),(1./nPP)*HPPt1.item(1),(-1./nP)*HPt1.item(1),(-1./nPP)*HPPt1.item(1)], \
          [DPt1.item(1),DPPt1.item(1),DPt1.item(1),DPPt1.item(1)], \
          [(-1./nP)*HPt1.item(0),(-1./nPP)*HPPt1.item(0),(1./nP)*HPt1.item(0),(1./nPP)*HPPt1.item(0)]])
        Psi1 = numpy.matrix([[inv(epsilonP).item(0,0),0,inv(epsilonP).item(0,1),0],[0,1.,0,0], \
          [inv(epsilonP).item(1,0),0,inv(epsilonP).item(1,1),0],[0,0,0,1.]])
        DeltaP = 1j*k0*deltaP
        DeltaPP = 1j*k0*deltaPP
        P = numpy.matrix([[numpy.exp(-DeltaP),0,0,0],[0,numpy.exp(-DeltaPP),0,0], \
          [0,0,numpy.exp(DeltaP),0],[0,0,0,numpy.exp(DeltaPP)]])
        
      # Store the transfer matrix fragment into an array
        Tr = Psi1*Phi1*inv(P)*inv(Phi1)*inv(Psi1)*Tr
        
  # Return the transfer matrix
    return Tr

# my version: calculate frequency-independent part of transfer matrix just once
def TransMatrix( noArr, neArr, oLTArr, eLTArr, tArr, chiArr, rho, incAngle, incN):
  # First medium is air
    theta1 = numpy.radians(incAngle)    # internal to this routine, all angles are radians
    n1 = incN
    
  # Create transfer matrix
    Tr = numpy.identity(4)

  # Loop over the layers
    for i in range(len(chiArr)):
        no = noArr[i]
        ne = neArr[i]
        oLT = oLTArr[i]
        eLT = eLTArr[i]
        t = tArr[i]
        chi = numpy.radians( chiArr[i]+rho )   #Plate orientation w.r.t x-axis

        nP = no
        nPP = ne*numpy.sqrt(1+(pow(ne,-2)-pow(no,-2))*sq(n1)*sq(sin(theta1))*sq(cos(chi)))
        
        R = lambda chi: numpy.matrix([[cos(chi),-sin(chi),0],[sin(chi),cos(chi),0],[0,0,1]])
        epsilonP = R(chi)*numpy.matrix([[sq(ne),0,0],[0,sq(no),0],[0,0,sq(no)]])*R(-chi)
        
        LTp = oLT
        nPT = nP*numpy.sqrt(1-1j*LTp)
     
      # I added this to handle isotropic dielectrics
        if ne != no :
          LTpp = oLT + ((eLT - oLT)/(ne - no))*(nPP - no)
        else :
          LTpp = oLT

        nPPT = nPP*numpy.sqrt(1-1j*LTpp)
        
        thetaP = numpy.arcsin(n1*sin(theta1)/nP)
        thetaPP = numpy.arcsin(n1*sin(theta1)/nPP)
        
        deltaP = nPT*t*cos(thetaP)
        deltaPP = nPPT*t*cos(thetaPP)
        
        DPt1 = pow(sq(cos(thetaP))+sq(sin(thetaP))*sq(sin(chi)),-0.5) * \
           numpy.array([-sin(chi)*cos(thetaP),cos(chi)*cos(thetaP),sin(chi)*sin(thetaP)])
        DPPt1 = pow(sq(cos(chi))*sq(cos(thetaP))+sq(sin(chi))*sq(cos(thetaP-thetaPP)),-0.5) * \
           numpy.array([cos(chi)*cos(thetaP)*cos(thetaPP),sin(chi)*(sin(thetaP)*sin(thetaPP)+cos(thetaP)*cos(thetaPP)), \
                  -cos(chi)*cos(thetaP)*sin(thetaPP)])
        HPt1 = pow(sq(cos(thetaP))*sq(cos(chi))+sq(sin(chi)),-0.5) * \
           numpy.array([-sq(cos(thetaP))*cos(chi),-sin(chi),cos(thetaP)*sin(thetaP)*cos(chi)])
        HPPt1 = pow(sq(cos(thetaP-thetaPP))*sq(sin(chi))+sq(cos(thetaP))*sq(cos(chi)),-0.5) * \
           numpy.array([-cos(thetaP-thetaPP)*cos(thetaPP)*sin(chi),cos(thetaP)*cos(chi),cos(thetaP-thetaPP)*sin(thetaPP)*sin(chi)])
        Phi1 = numpy.matrix([[DPt1.item(0),DPPt1.item(0),DPt1.item(0),DPPt1.item(0)], \
          [(1./nP)*HPt1.item(1),(1./nPP)*HPPt1.item(1),(-1./nP)*HPt1.item(1),(-1./nPP)*HPPt1.item(1)], \
          [DPt1.item(1),DPPt1.item(1),DPt1.item(1),DPPt1.item(1)], \
          [(-1./nP)*HPt1.item(0),(-1./nPP)*HPPt1.item(0),(1./nP)*HPt1.item(0),(1./nPP)*HPPt1.item(0)]])
        Psi1 = numpy.matrix([[inv(epsilonP).item(0,0),0,inv(epsilonP).item(0,1),0],[0,1.,0,0], \
          [inv(epsilonP).item(1,0),0,inv(epsilonP).item(1,1),0],[0,0,0,1.]])

    return deltaP, deltaPP, Psi1, Phi1, inv(Phi1), inv(Psi1), Tr

# this is the frequency-dependent part
def TransMat1( nu, deltaP0, deltaPP0, Psi1, Phi1, invPhi1, invPsi1, Tr0 ):
    k0 = (2*math.pi)/(clight/nu)
    DeltaP = 1j*k0*deltaP0
    DeltaPP = 1j*k0*deltaPP0
    P = numpy.matrix([[numpy.exp(-DeltaP),0,0,0],[0,numpy.exp(-DeltaPP),0,0], \
          [0,0,numpy.exp(DeltaP),0],[0,0,0,numpy.exp(DeltaPP)]])
    Tr = Psi1*Phi1*inv(P)*invPhi1*invPsi1*Tr0
    return Tr
        

#Arguments (all dielectric arrays should NOT include the air layers)
#  - nuArr: array of frequencies to be analyzed [GHz *changed from Tom's code]
#  - noArr: array of ordinary dielectric constants
#  - neArr: array of extraordinary dielectric constants
#  - oLTArr: array of ordinary loss tangents
#  - eLTArr: array of extraordinary loss tangents
#  - tArr: array of the thicknesses of the dielectric layers [m]
#  - chiArr: array of the angles of the extraordinary axes w.r.t. the stack angle (see "rho" below) [degrees]
#  - rho: angle of the stack (i.e. the "stack angle") w.r.t. the plane of incidence, as measured in the plane of rotation of the stack [degrees]
#  - incAngle: angle of incidence [degrees]
#  - polAngle: polarization angle as measured from the plane of incidence, as measured in the plane of the propogating wave [degrees]
#Returns
#  - p-polarized transmission
#  - s-polarized transmission
#  - p-polarized reflection
#  - s-polarized reflection
#
def anisotropicTransmissionPol(nuArr, noArr, neArr, oLTArr, eLTArr, tArr, chiArr, rho=0., incAngle=0., polAngle=0.):

  # Check array integrity
    if not (len(chiArr) == len(noArr) == len(neArr) == len(tArr) == len(oLTArr) == len(eLTArr)):
        raise NameError("Length of 'anisotropicTransmission' arrays need to be the same")
    if not len(nuArr):
        raise NameError("Length 'nuArr' needs to be non-zero")

  # First barrier is an air barrier
    n1 = 1.
    theta1 = numpy.radians(incAngle)

  # Create transfer matrix for this stack
    # T = lambda nu: TransMat(nu, noArr, neArr, oLTArr, eLTArr, tArr, chiArr, rho, incAngle, n1)

  # Calculate frequency-independent components of transfer matrix once
    deltaP0, deltaPP0, Psi1, Phi1, invPhi1, invPsi1, Tr0 \
       = TransMatrix( noArr, neArr, oLTArr, eLTArr, tArr, chiArr, rho, incAngle, n1)

  # Full transfer matrix depends on frequency
    T = lambda nu: TransMat1( nu, deltaP0, deltaPP0, Psi1, Phi1, invPhi1, invPsi1, Tr0 )
        
  # Calculate interface to the final air barrier
    n3 = 1.
    theta3 = numpy.arcsin((n1/n3)*sin(theta1))

  # Transmission coefficients
    alphaT = lambda nu: (T(nu).item(0,0)*cos(theta3)+T(nu).item(0,1)*n3)/cos(theta1)
    betaT = lambda nu: (T(nu).item(0,2)+T(nu).item(0,3)*n3*cos(theta3))/cos(theta1)
    gammaT = lambda nu: (T(nu).item(1,0)*cos(theta3)+T(nu).item(1,1)*n3)/n1
    deltaT = lambda nu: (T(nu).item(1,2)+T(nu).item(1,3)*n3*cos(theta3))/n1
    etaT = lambda nu: (T(nu).item(2,0)*cos(theta3)+T(nu).item(2,1)*n3)
    kappaT = lambda nu: (T(nu).item(2,2)+T(nu).item(2,3)*n3*cos(theta3))
    rhoT = lambda nu: (T(nu).item(3,0)*cos(theta3)+T(nu).item(3,1)*n3)/(n1*cos(theta1))
    sigmaT = lambda nu: (T(nu).item(3,2)+T(nu).item(3,3)*n3*cos(theta3))/(n1*cos(theta1))
    GammaT = lambda nu: pow((alphaT(nu)+gammaT(nu))*(kappaT(nu)+sigmaT(nu))-(betaT(nu)+deltaT(nu))*(etaT(nu)+rhoT(nu)),-1.)
    
  # Jones Transmission and Reflection Matrices
    JT = lambda nu: 2.*GammaT(nu)*numpy.matrix(
        [[(kappaT(nu)+sigmaT(nu)),
          (-betaT(nu)-deltaT(nu))],
         [(-etaT(nu)-rhoT(nu)),
          (alphaT(nu)+gammaT(nu))]]
        )
    JR = lambda nu: GammaT(nu)*numpy.matrix(
        [[(gammaT(nu)-alphaT(nu))*(kappaT(nu)+sigmaT(nu))-(deltaT(nu)-betaT(nu))*(etaT(nu)+rhoT(nu)),
          2.*(alphaT(nu)*deltaT(nu)-gammaT(nu)*betaT(nu))],
         [2.*(etaT(nu)*sigmaT(nu)-rhoT(nu)*kappaT(nu)),
          (alphaT(nu)+gammaT(nu))*(kappaT(nu)-sigmaT(nu))-(betaT(nu)+deltaT(nu))*(etaT(nu)-rhoT(nu))]]
        )

  # Incident Electric Field Amplitude
    Ei = numpy.array([abs(cos(numpy.radians(polAngle))), abs(sin(numpy.radians(polAngle)))])
    
  # Transmitted Electric Field Amplitude
    EPT = lambda nu: JT(nu).item(0,0)*Ei[0] + JT(nu).item(0,1)*Ei[1]
    EST = lambda nu: JT(nu).item(1,0)*Ei[0] + JT(nu).item(1,1)*Ei[1]

  # Reflected Electric Field Amplitude
    EPR = lambda nu: JR(nu).item(0,0)*Ei[0] + JR(nu).item(0,1)*Ei[1]
    ESR = lambda nu: JR(nu).item(1,0)*Ei[0] + JR(nu).item(1,1)*Ei[1]

  # Loop over the frequencies and calculate transmission for each
    pT = []
    sT = []
    pR = []
    sR = []
    for nu in nuArr:
        #print "Analyzing %.1f GHz" % (nu)
        pT.append(EPT(nu))
        sT.append(EST(nu))
        pR.append(EPR(nu))
        sR.append(ESR(nu))
        
  # Return the amplitude arrays
    return numpy.array(pT), numpy.array(sT), numpy.array(pR), numpy.array(sR)

# =====================================================================================
# stackDesc describes the stack of dielectric layers
# to maintain compatibility with Tom's code, stackDesc is a list of arrays
# each of these arrays contains the relevant parameter for each layer in the stack --
#   so, for example, each array will have 2 elements if there are 2 layers in the stack
# data files 100_30deg.dat and similar are from June; in these, 3mm receiver was horizontally
#   polarized (polAngle=90) and 1mm receiver was LCP; plate thickness was 1.009 cm
# data files 220_pos10.dat and simlar are from Sep; in these, both bands were vertically 
#   polarized (polAngle=0); plate thickness was .376 cm
# parameters for jul 2015 data
#
noArr = numpy.array([3.06]) #Actually measured
neArr = numpy.array([3.39]) #Actually measured
oltArr = numpy.array([0.5e-4]) #Estimated from Dick's and my measurement
eltArr = numpy.array([1.5e-4]) #Estimated from Dick's and my measurement
dArr = numpy.array([1.009]) #[cm], Actually measured
chiArr = numpy.array([0.]) #[deg]
stackDesc = [ noArr, neArr, oltArr, eltArr, dArr, chiArr ]
#
# 3 other parameters are required:
#   stackRotAngle is the azimuthal rotation of the entire stack
#   stackTilt is the angle of incidence of radiation to the stack
#   polAngle describes the angle of linear polarization relative to the plane of incidence;
#     e.g., polAngle = 90 describes "S" polarization, perpendicular to plane of incidence
#
stackRotAngle = 40. #[deg], arbitrary (but maybe one of the angles that Dick and I measured?)
# polAngle = 90. # 3mm polarization (horiz) is perp to the plane of incidence 
stackTilt = 42.5 #[deg]
#
# parameters for sep 2015 data:
# dArr = numpy.array([3.76e-1]) #[cm], Actually measured
# polAngle = 0. #Polarization is aligned with the plane of incidence (maximize transmission)
# =====================================================================================

# saphModel is designed to model thermal emission data from summer 2015, where LSB and USB
#   are folded together
def saphModel( LOGHz, IFfrqArray, stackDesc, stackRotAngle, stackTilt, polAngle, outFile  ) :
   noArr = stackDesc[0]
   neArr = stackDesc[1]
   oltArr = stackDesc[2]
   eltArr = stackDesc[3]
   dArr = stackDesc[4]
   chiArr = stackDesc[5]
   if outFile :
     fout = open( outFile, "w" )
     fout.write("# LOGHz = %.3f\n" % LOGHz)
     for n in range(0,len(noArr)) :
       fout.write("# layer %d: no=%.3f  ne=%.3f  olt=%.2e  elt=%.2e  dcm=%.3f  chi=%.1f\n" % \
          (n,noArr[n],neArr[n],oltArr[n],eltArr[n],dArr[n],chiArr[n]) )
     fout.write("# stackRotAngle = %.2f\n" % stackRotAngle)
     fout.write("# stackTilt = %.2f\n" % stackTilt)
     fout.write("# polAngle = %.2f\n" % polAngle) 
     fout.write("# frq  P3  P4  P5\n")
   nuLSBArr = LOGHz - IFfrqArray
   nuUSBArr = LOGHz + IFfrqArray 
   pTLSB, sTLSB, pRLSB, sRLSB = anisotropicTransmissionPol(nuLSBArr, noArr, neArr, oltArr, eltArr, dArr, \
      chiArr, stackRotAngle, stackTilt, polAngle)
   pTUSB, sTUSB, pRUSB, sRUSB = anisotropicTransmissionPol(nuUSBArr, noArr, neArr, oltArr, eltArr, dArr, \
      chiArr, stackRotAngle, stackTilt, polAngle)
   T = 0.5 * ( (pTLSB * numpy.conj(pTLSB) + sTLSB * numpy.conj(sTLSB) ).astype(float) \
             + (pTUSB * numpy.conj(pTUSB) + sTUSB * numpy.conj(sTUSB) ).astype(float) )
   R = 0.5 * ( (pRLSB * numpy.conj(pRLSB) + sRLSB * numpy.conj(sRLSB) ).astype(float) \
             + (pRUSB * numpy.conj(pRUSB) + sRUSB * numpy.conj(sRUSB) ).astype(float) )
   A = numpy.ones( len(T), dtype=float) - T - R 
   P3 = 295.*T + 77.*R + 295.*A
   P4 = 77.*T + 295.*R + 295.*A
   P5 = 77.*T + 77.*R + 295.*A
   if outFile :
     for n in range(0,len(T)) :
       fout.write("%8.3f  %8.2f  %8.2f  %8.2f\n" % (IFfrqArray[n],P3[n],P4[n],P5[n]))
     fout.close()
   return P3,P4,P5

# read P3,P4,P5 temperatures from input file, vs freq
def saphReadData( infile ) :
    fin = open( infile, "r" )
    LOGHz = 0.
    phiOffset = 0.
    polAngle = 0.
    LOGHzMissing = phiOffsetMissing = polAngleMissing = True
    fGHz = []
    P3 = []
    P4 = []
    P5 = []
    P6 = []
    for line in fin :
      if line.startswith( "##" ) :     # ignore lines beginning with double hash mark
        continue
      elif not line.startswith( "#" ) :    # lines without hash marks are data
        a = line.split()
        fGHz.append( float( a[0] ) )
        P3.append( float(a[1]) )
        P4.append( float(a[2]) )
        P5.append( float(a[3]) )
        P6.append( float(a[4]) )
      else :                               # lines with single hash marks are search parameters
        if "LOGHz" in line :
          LOGHz = float( line[line.find("=")+1:] )
          LOGHzMissing = False
        elif "phiOffset" in line :
          phiOffset = float( line[line.find("=")+1:] )
          phiOffsetMissing = False
        elif "polAngle" in line :
          polAngle = float( line[line.find("=")+1:] )
          polAngleMissing = False
    fin.close()
    if LOGHzMissing or phiOffsetMissing or polAngleMissing :
      print "... saphReadData WARNING - missing parameter" 
    print "... read data from %s -- npts = %d, LOGHz = %.1f, phiOffset = %.1f deg, polAngle = %.1f deg" % \
       (infile, len(P3), LOGHz, phiOffset, polAngle)
    return numpy.array(fGHz), numpy.array(P3), numpy.array(P4), numpy.array(P5), numpy.array(P6), LOGHz, phiOffset, polAngle

# use spline fits to find smoothed values at frequencies in frqArray (units: GHz)
def saphSmooth( fGHzIn, Tin, frqArray ) :
    #tck = scipy.interpolate.splrep( fGHzIn, Tin, k=3, task=-1, t=frqArray )
    #Tsm = scipy.interpolate.splev( frqArray, tck )
    tck = splrep( fGHzIn, Tin, k=3, task=-1, t=frqArray )
    Tsm = splev( frqArray, tck )
    return Tsm

def saphPlot( infile,theta, fGHz,P3,P4,P5, frqArray,P3m,P4m,P5m ) :
    #pyplot.figure(0)
    #pyplot.plot( fGHz, P3, color='b'  )
    #pyplot.plot( frqArray, P3m, color='b', linestyle='--' )
    #pyplot.plot( fGHz, P4, color='g'  )
    #pyplot.plot( frqArray, P4m, color='g', linestyle='--' )
    #pyplot.plot( fGHz, P5, 'r'  )
    #pyplot.plot( frqArray, P5m, 'r-' )
    #pyplot.show()
    pyplot.ioff()
    pp = PdfPages("saphFit.pdf")
    pyplot.figure( figsize=(11,8) )
    ax = pyplot.subplot(1,1,1)
    ax.plot( fGHz, P3, color='b'  )
    ax.plot( frqArray, P3m, color='b', linestyle='--' )
    ax.plot( fGHz, P4, color='g'  )
    ax.plot( frqArray, P4m, color='g', linestyle='--' )
    ax.plot( fGHz, P5, 'r'  )
    ax.plot( frqArray, P5m, 'r-' )
    pyplot.suptitle( "%s   theta=%.1f deg" % (infile,theta), fontsize=14 )
    pyplot.savefig( pp, format="pdf" )
    pp.close()
    pyplot.show()

# this is my first attempt to model jun 2015 data; vary only theta
def saphFit( infile, thetaList=numpy.arange(0.,90.1,10.) ) :
    fGHz,P3,P4,P5,LOGHz,phiOffset,polAngle = saphReadData( infile )
    frqArray = numpy.array( numpy.arange(.3,9.91,.3) )
    P3sm = saphSmooth( fGHz, P3, frqArray )
    P4sm = saphSmooth( fGHz, P4, frqArray )
    P5sm = saphSmooth( fGHz, P5, frqArray )
    if len(thetaList) > 1 :
      bestResid = 1.e6
      thetaBest = -1000.
      for theta in thetaList :
        # P3m, P4m, P5m =  saphModel( LOGHz, frqArray, None, stackRotAngle=theta ) 
        P3m, P4m, P5m = saphModel( LOGHz, frqArray, stackDesc, theta+phiOffset, stackTilt, polAngle, None  ) 
        resid3 = numpy.std( P3sm - P3m )
        resid4 = numpy.std( P4sm - P4m )
        print "... theta = %.1f  resid3 = %.3f  resid4 = %.3f" % (theta,resid3,resid4)
        saphPlot( infile, theta, fGHz,P3,P4,P5, frqArray,P3m,P4m,P5m ) 
        if (resid3 + resid4) < bestResid :
          bestResid = resid3 + resid4
          thetaBest = theta
      print "theta = %.1f gives best fit" % thetaBest
    else :
      thetaBest = thetaList[0]
    P3m, P4m, P5m = saphModel( LOGHz, frqArray, stackDesc, thetaBest+phiOffset, stackTilt, polAngle, None  ) 
    resid3 = numpy.std( P3sm - P3m )
    resid4 = numpy.std( P4sm - P4m )
    print "... theta = %.1f  resid3 = %.3f  resid4 = %.3f" % (thetaBest,resid3,resid4)
    saphPlot( infile, thetaBest, fGHz,P3,P4,P5, frqArray,P3m,P4m,P5m ) 
      
    
# saphFit2 is adapted from pol.doit6()
# it is designed to fit data taken at multiple angles; looks for lowest variance overall
# I imagine doing this one freq at a time, so plots 4 panels with final results
# ... dataFileList = list of data file names, e.g., { 'file1.dat', 'file2.dat' ]
# ... srchList is a dictionary of search ranges; if values are not filled out in it, srchDefault value is used
# ... T5fit to find min variance in T5 (normally fits T3 and T4)
# ... T6override to manually set temperature seen by reflection off the plate; this is important only for T5fit=True
# unlike doit6, spline fit is used to obtain smoothed, measured values for freqArray
# srcDefault ranges are below

stackDesc = [ noArr, neArr, oltArr, eltArr, dArr, chiArr ]

tcmSrch = [1.009]            # plate thickness in cm
noSrch = [3.050]             # ordinary index
neSrch = [3.400]             # extraordinary index
oltSrch = [0.5e-4]           # ordinary loss tangent
eltSrch = [1.5e-4]           # extraordinary loss tangent
angISrch = [43.,44.,45.]     # angle of incidence ('stackTilt')
thetaSrch = [87.]            # stack rotation angle
srchList = [ noSrch, neSrch, oltSrch, eltSrch, tcmSrch, angISrch, thetaSrch ]

def saphFit2( dataFileList, srchList, T5fit=False, T6override=None,  ) :
    fout = open("saphFit2.log", "a")
    fout.write("#\n#\n# %s\n" % (dataFileList) )
    fout.close()

  # read in and smooth the data files; each file corresponds to a particular LOfreq, phiOffset, polAngle
    frqArray = numpy.array( numpy.arange(.3,9.91,.3) )
    P3sm = []
    P4sm = []
    P5sm = []
    P6avg = []
    LOfrq = []
    phiOff = []
    polAng = []
    for dataFile in dataFileList :  
      fGHz,P3,P4,P5,P6,LOGHz,phiOffset,polAngle = saphReadData( dataFile )
      P3sm.append( saphSmooth( fGHz, P3, frqArray ) )
      P4sm.append( saphSmooth( fGHz, P4, frqArray ) )
      P5sm.append( saphSmooth( fGHz, P5, frqArray ) )
      LOfrq.append( LOGHz )
      phiOff.append( phiOffset )
      polAng.append( polAngle )
    if T6override :
      P6avg.append = T6override
    else :
      P6avg.append( numpy.mean( P6 ) )
    if T5fit :
      print "... looking for min variance in T5; using T6avg = %.2f" % T6avg
    
  # fit all, searching for min variance
    ntot = len(srchList[0]) * len(srchList[1]) * len(srchList[2]) * len(srchList[3]) \
             * len(srchList[4]) * len(srchList[5]) * len(srchList[6])
    print "... modeling %d possibilities" % ntot
    nfinished = 0
    varsave = 1.e6
    for no in srchList[0] :
      for ne in srchList[1] :
        for olt in srchList[2] :
          for elt in srchList[3] :
            for tcm in srchList[4] :
              for angI in srchList[5] :
                for theta in srchList[6] :

                # form the stackDesc arrays (one value per layer, in this case) needed by saphModel
                  stackDesc = [ no, ne, olt, elt, tcm, 0. ]    

                # model each data file separately (could have different LOfrq, phiOff, or polAng)
                  var = 0.
                  var5 = 0.
                  for n in range(0,len(LOGHz)) :
                    P3m,P4m,P5m = saphModel( LOGHz, frqArray, stackDesc, theta+phiOff[n], angI, polAng[n], None  ) 
                    var = var + numpy.mean( pow( T3m-T3sm[n], 2. )) + numpy.mean( pow( T4m-T4sm[n], 2. ))
                    var5 =  var5 + numpy.mean( pow( T5m-T5sm[n], 2. ))
                  print "tcm = %.3f, no = %.3f, ne = %.3f, angI = %.1f, phi = %.1f, alphao,e = %.2e, %.2e, var34 = %.1f, var5 = %.2f" % \
                      (tcm, nO, nE, angI, phi, alphaO, alphaE, variance, var5)
                  fout = open("doit6.log", "a")
                  fout.write( "tcm: %.3f  nO: %.3f  nE: %.3f  angI: %.1f  phi: %.1f  alphaO,E: %.3f %.3f  variance: %.1f  var5: %.2f\n" % \
                    (tcm, nO, nE, angI, phi, alphaO, alphaE, variance, var5))
                  fout.close()
                  if T5fit and (var5 < varsave) :
                    tcmbest = tcm
                    nObest = nO
                    nEbest = nE
                    angIbest = angI
                    phibest = phi 
                    alphaObest = alphaO
                    alphaEbest = alphaE
                    varsave = var5
                   elif (variance < varsave) :
                  tcmbest = tcm
                  nObest = nO
                  nEbest = nE
                  angIbest = angI
                  phibest = phi 
                  alphaObest = alphaO
                  alphaEbest = alphaE
                  varsave = variance

  print "# BEST: tcm = %.3f, nO = %.3f, nE = %.3f, angI = %.1f, phi = %.1f, alphaO,E = %.3f, %.3f, var34 = %.1f" % \
             (tcmbest, nObest, nEbest, angIbest, phibest, alphaObest, alphaEbest, varsave )
  fout = open("doit6.log", "a")
  fout.write("# BEST tcm: %.3f  nO: %.3f  nE: %.3f  angI: %.1f  phi: %.1f  alphaO,E: %.3f %.3f  var34: %.1f" % \
             (tcmbest, nObest, nEbest, angIbest, phibest, alphaObest, alphaEbest, varsave) )
  fout.close()

# plot data for all files; also dump results to "fit.dat"
  plt.ioff()
  nplot = 1
  nmax = len( dataFileList)
  nm = 1
  fsize = 12
  if nmax > 1 :
    nm = 2 
    fsize = 6

  fout = open("fit.dat", "w")
  fout.write( '# created by pol.doit6()\n' )
  fout.write( '# dataFileList:  ' + ', '.join(dataFileList) + '\n' )
  fout.write( '# phiOffsetList: ' + ', '.join( [str(q) for q in phiOffsetList] ) + '\n#\n' )
  fout.write( "# tcm    : [ %.3f, %.3f, %.3f ] best %.3f\n" % \
     ( srch['tcm'][0], srch['tcm'][1], srch['tcm'][2], tcmbest ) )
  fout.write( "# nO     : [ %.3f, %.3f, %.3f ] best %.3f\n" % \
     ( srch['nO'][0], srch['nO'][1], srch['nO'][2], nObest ) )
  fout.write( "# nE     : [ %.3f, %.3f, %.3f ] best %.3f\n" % \
     ( srch['nE'][0], srch['nE'][1], srch['nE'][2], nEbest ) )
  fout.write( "# angI   : [ %.1f, %.1f, %.1f ] best %.1f\n" % \
     ( srch['angI'][0], srch['angI'][1], srch['angI'][2], angIbest ) )
  fout.write( "# phi    : [ %.1f, %.1f, %.1f ] best %.1f\n" % \
     ( srch['phi'][0], srch['phi'][1], srch['phi'][2], phibest ) )
  fout.write( "# alphaO : [ %.3f, %.3f, %.3f ] best %.3f\n" % \
     ( srch['alphaO'][0], srch['alphaO'][1], srch['alphaO'][2], alphaObest ) )
  fout.write( "# alphaE is constrained to equal alphaO\n" )
  if (T5fit) :
    fout.write( "# minimized (T5-T5dat)^2\n")
  else :
    fout.write( "# minimized (T3-T3dat)^2 + (T4-T4dat)^2\n")

    
  for dataFile, phiOff in zip( dataFileList, phiOffsetList) :  
    fdat, T3dat, T4dat, T5dat, T6dat = readDataFile( dataFile )
    if nSmooth > 0 :
      fdat = smArray( fdat, nSmooth )
      T3dat = smArray( T3dat, nSmooth )
      T4dat = smArray( T4dat, nSmooth ) 
      T5dat = smArray( T5dat, nSmooth ) 
      T6dat = smArray( T6dat, nSmooth ) 
    phOff = phiOff * numpy.ones( len(fdat), dtype=float )
    T3,T4,T5 = Tfit2( LOGHz, fdat, phOff, angIbest, phibest, alphaObest, alphaEbest, nObest, nEbest, tcmbest, T5fit, T6avg )  

    fout.write( '#\n# %s\n' % dataFile )
    for ifq in range(0,len(fdat)) :
      fout.write("%8.5f  %6.2f %6.2f  %6.2f %6.2f  %6.2f %6.2f  %6.2f %6.2f\n" % \
        ( fdat[ifq], T3dat[ifq], T3[ifq], T4dat[ifq], T4[ifq], T5dat[ifq], T5[ifq], T6dat[ifq], T6avg ) )
    fout.close()

    fig = plt.subplot(nm,nm,nplot)
    fig.plot( fdat, T3, color='red' )
    fig.plot( fdat, T3dat, color='red' )
    fig.plot( fdat, T4, color='blue' )
    fig.plot( fdat, T4dat, color='blue' )
    fig.plot( fdat, T5, color='purple' )
    fig.plot( fdat, T5dat, color='purple' )
    fig.grid( True )
    fig.set_title( dataFile, fontsize=12 )
    fig.set_ylim(0,300)
    fig.text( .05, .1, "tcm = %.3f, nO = %.3f, nE = %.3f, variance = %.1f" % \
       ( tcmbest, nObest, nEbest, varsave), \
       transform=fig.transAxes, horizontalalignment="left", verticalalignment="center", fontsize=fsize )
    fig.text( .05, .05, "phi = %.1f, angI = %.1f, alphaO = %.3f, alphaE = %.3f" % (phibest,angIbest,alphaObest,alphaEbest), \
       transform=fig.transAxes, horizontalalignment="left", verticalalignment="center", fontsize=fsize )
    nplot = nplot + 1
  plt.show( )

# converts alpha to tanDelta, using average index of 3.25

def tanD (fGHz, alpha) :
  print "tanDelta = %.2e" %  (alpha * clight/ (2.*math.pi*3.25*fGHz) )
"""
