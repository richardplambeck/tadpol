# ob.py
# these are routines for dealing with the Optics Bench scans of various kinds

import numpy
import math
import sys
import time
import string
import pickle
import matplotlib
import datetime
#matplotlib.use('GTK')
matplotlib.use('GTKAgg')
import matplotlib.pyplot as pyplot
from matplotlib.backends.backend_pdf import PdfPages
#import pol
import scipy
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
    elif nlayers == 6 :
      a = numpy.meshgrid( tcmListIn[0],tcmListIn[1],tcmListIn[2],tcmListIn[3],tcmListIn[4],tcmListIn[5], \
            nrListIn[0],nrListIn[1],nrListIn[2],nrListIn[3],nrListIn[4],nrListIn[5] )
    elif nlayers == 7 :
      a = numpy.meshgrid( tcmListIn[0],tcmListIn[1],tcmListIn[2],tcmListIn[3],tcmListIn[4],tcmListIn[5],tcmListIn[6], \
            nrListIn[0],nrListIn[1],nrListIn[2],nrListIn[3],nrListIn[4],nrListIn[5],nrListIn[6] )
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
    Resid = 1.e6*numpy.ones( ilen )
    secs0 = time.time()
    for i in range(0,ilen) :
      dphs_est,trans_est = solveStack( fGHz, angIdeg, tcmList[i], nrList[i], tanDelta ) 
      dphi = unwrap( dphs_meas-dphs_est )
        # unwrap keeps phase in range -180 to 180 in all cases
      avgPhResid[i] = numpy.average(dphi) 
      stdPhResid[i] = numpy.std(dphi) 
      avgAmpResid[i] = numpy.average( trans - trans_est )
      stdAmpResid[i] = numpy.std( trans - trans_est )
      # Resid[i] = avgPhResid[i] + stdPhResid[i] + stdAmpResid[i]
    # phases assumed uncertain by 2 degree, amps by .02; so upweight amps by factor of 100
      Resid[i] =  math.sqrt( pow(avgPhResid[i],2) + pow(stdPhResid[i],2) + 1.e4*pow(stdAmpResid[i],2) )
      if i > 0 and i%200 == 0 :
        secselapsed = time.time() - secs0
        secsremaining = (ilen-i)*(secselapsed/float(i))
        if secsremaining/60. > 60. :
          print "... %d/%d finished; %.2f hours remaining" % (i,ilen,secsremaining/3600.)
        else :
          print "... %d/%d finished; %.2f minutes remaining" % (i,ilen,secsremaining/60.)
      #print tcmList[i], nrList[i], ".. %7.3f, %7.3f" \
      #       % (avgPhResid[i], stdPhResid[i])

    #merit = stdPhResid + stdAmpResid
   
  # save the results for further analysis
    fout = open( infile+"Fits.pickle", "wb" )     # write over previous file
    pickle.dump( [tcmList,nrList,avgPhResid,stdPhResid,stdAmpResid,Resid], fout ) 
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
    printBest( stdPhResid, tcmList, nrList ) 
    # fout.close()

    print "\n10 best fits minimizing overall residuals:"
    nbest = printBest( Resid, tcmList, nrList ) 

  # begin the plot
    pyplot.ioff()
    pp = PdfPages("refrac.pdf")
    fig =  pyplot.figure( figsize=(11,8) )
    pyplot.suptitle( "%s (%s)" % (title,infile), fontsize=14 )

  # lower left plot will show mean phase resid vs refractive index
    #ax = fig.add_axes( [.1,.1,.38,.2] )
    #xmin,xmax = minmax( numpy.array(nrList[0]) )
    #ymin,ymax = minmax( numpy.array(stdPhResid) )
    #ax.axis( [xmin, xmax, ymin, ymax] )
    #ax.tick_params( labelsize=10 )
    # ax.errorbar( nrList[0], stdPhResid, yerr=std, linewidth=2 )
    # ax.errorbar( nrfit, 0., xerr=fitunc, color='red', elinewidth=6 )
    #ax.set_xlabel("refractive index", fontsize=10)
    #ax.set_ylabel("mean phase residual (deg)", fontsize=10)
    #ax.text( 0.92, .82, "best fit = %.3f (+/- %.3f)" % (nrfit,fitunc), transform=ax.transAxes, \
    #  horizontalalignment='right', fontsize=10, color="black", rotation='horizontal' )
    #pyplot.grid(True)

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
      ax.text( 0.01, ylab, "layer %d:  t=%.4f\"  nr=%.3f  eps=%.3f" % \
        ( nl+1, tcmList[nbest][nl]/2.54, nrList[nbest][nl], pow(nrList[nbest][nl],2) ),  \
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
       [tcm], [0.], rho=0., incAngle=angIdeg, Ei=numpy.array([0.,1.]) )
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
#  - replaced by Ei = numpy.array( Epar, Eperp )
#Returns
#  - p-polarized transmission
#  - s-polarized transmission
#  - p-polarized reflection
#  - s-polarized reflection
#
def anisotropicTransmissionPol(nuArr, noArr, neArr, oLTArr, eLTArr, tArr, chiArr, rho=0., incAngle=0., Ei=numpy.array([0.,1.]) ):

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
    # Ei = numpy.array([abs(cos(numpy.radians(polAngle))), abs(sin(numpy.radians(polAngle)))])
    
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

# I am searching for some kind of conceptual problem in dealing with circular pol
def doit7() :
   nuArr = numpy.array( [220.] )
   noArr = numpy.array( [3.066] )
   neArr = numpy.array( [3.404] )
   oltArr = numpy.array( [2.e-4] )
   eltArr = numpy.array( [2.e-4] )
   dArr = numpy.array( [1.008] )
   chiArr = numpy.array( [0.] )
   stackRotAngle = 0. 
   stackTilt = 44.
   for Ei in [ numpy.array([0.,1.]), numpy.array([1.,0.]), numpy.array([.707,1j*.707]), numpy.array([.707,-1j*.707]) ] :
     pT, sT, pR, sR = anisotropicTransmissionPol(nuArr, noArr, neArr, oltArr, eltArr, dArr, \
        chiArr, stackRotAngle, stackTilt, Ei)
     print "\nEi: ", Ei
     print "pT: %.4f" % (pT*numpy.conj(pT)).astype(float), pT
     print "sT: %.4f" % (sT*numpy.conj(sT)).astype(float), sT
     print "pR: %.4f" % (pR*numpy.conj(pR)).astype(float), pR
     print "sR: %.4f" % (sR*numpy.conj(sR)).astype(float), sR 

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
#noArr = numpy.array([3.06]) #Actually measured
#neArr = numpy.array([3.39]) #Actually measured
#oltArr = numpy.array([0.5e-4]) #Estimated from Dick's and my measurement
#eltArr = numpy.array([1.5e-4]) #Estimated from Dick's and my measurement
#dArr = numpy.array([1.009]) #[cm], Actually measured
#chiArr = numpy.array([0.]) #[deg]
#stackDesc = [ noArr, neArr, oltArr, eltArr, dArr, chiArr ]
#
# 3 other parameters are required:
#   stackRotAngle is the azimuthal rotation of the entire stack
#   stackTilt is the angle of incidence of radiation to the stack
#   polAngle describes the angle of linear polarization relative to the plane of incidence;
#     e.g., polAngle = 90 describes "S" polarization, perpendicular to plane of incidence
#
#stackRotAngle = 40. #[deg], arbitrary (but maybe one of the angles that Dick and I measured?)
# polAngle = 90. # 3mm polarization (horiz) is perp to the plane of incidence 
#stackTilt = 42.5 #[deg]
#
# parameters for sep 2015 data:
# dArr = numpy.array([3.76e-1]) #[cm], Actually measured
# polAngle = 0. #Polarization is aligned with the plane of incidence (maximize transmission)
# =====================================================================================

# computes P3,P4,P5 for one sideband
def saphModel1( nuArr, stackDesc, stackRotAngle, stackTilt, Ei ) :
   noArr = stackDesc[0]
   neArr = stackDesc[1]
   oltArr = stackDesc[2]
   eltArr = stackDesc[3]
   dArr = stackDesc[4]
   chiArr = stackDesc[5]
   pT, sT, pR, sR = anisotropicTransmissionPol(nuArr, noArr, neArr, oltArr, eltArr, dArr, \
      chiArr, stackRotAngle, stackTilt, Ei)
   T = (pT * numpy.conj(pT) + sT * numpy.conj(sT) ).astype(float) 
   R = (pR * numpy.conj(pR) + sR * numpy.conj(sR) ).astype(float)

   # verified that the following expressions give exactly same answer as those above,
   #   but take slightly longer (7.3 minutes vs 6.8 mins, for example)
   # T = numpy.power( numpy.absolute( pT), 2.) + numpy.power( numpy.absolute(sT), 2.)
   # R = numpy.power( numpy.absolute( pR), 2.) + numpy.power( numpy.absolute(sR), 2.)

   A = numpy.ones( len(T), dtype=float) - T - R
   P3 = 291.*T + 78.*R + 294.*A
   P4 = 80.*T + 294.*R + 294.*A
   P5 = 80.*T + 78.*R + 294.*A
   return numpy.array( [P3,P4,P5] )

# compute model DSB temperatures, and variance from smoothed measured data,
#    given model LSB and USB temperatures
# TmL, TmU, and Tsm are arrays each containing [T3,T4,T5]
# Lfrac is the LSB fraction, which should be close to 0.5
# if Lfrac is negative, saphDSB searches for the best fit value between 0.4 and 0.6 
#    that gives lowest variance
#
def saphDSB( TmL, TmU, Tsm, Lfrac=0.5 ):
    if Lfrac < 0. :
      varmin = 1.e6
      for Lf in numpy.arange( 0.4, 0.61, .02 ) :
        Tm = Lf*TmL + (1.-Lf)*TmU
        var34 = numpy.mean( numpy.power( Tm[0]-Tsm[0], 2.) ) + numpy.mean( numpy.power( Tm[1]-Tsm[1], 2. ))
        var5 = numpy.mean( numpy.power( Tm[2]-Tsm[2], 2.) )
        print "Lf = %.2f, var34 = %.1f" % (Lf,var34)
        if var34 < varmin :
          varmin = var34
          Lfbest = Lf
    else :
      Lfbest = Lfrac
    Tm = Lfbest*TmL + (1.-Lfbest)*TmU
    var34 = numpy.mean( numpy.power( Tm[0]-Tsm[0], 2.) ) + numpy.mean( numpy.power( Tm[1]-Tsm[1], 2. ))
    var5 = numpy.mean( numpy.power( Tm[2]-Tsm[2], 2.) )
    return Tm, Lfbest, var34, var5

# read P3,P4,P5 temperatures from input file, vs freq
def saphReadData( infile, logfile="saphFit2.log" ) :
    fin = open( infile, "r" )
    LOGHz = 0.
    rhoOffset = 0.
    Lfrac = 0.5
    LOGHzMissing = rhoOffsetMissing = rcvrPolMissing = LfracMissing = True
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
        elif "Lfrac" in line :
          Lfrac = float( line[line.find("=")+1:] )
          LfracMissing = False
        elif "rhoOffset" in line :
          rhoOffset = float( line[line.find("=")+1:] )
          rhoOffsetMissing = False
      # receiver polarization (H,V,L,R) - will be converted to Ei = [Epar, Eperp]
      # plane of incidence is V (vertical reflection into bucket); this is polAngle=0 in Tom's code
      # 2/20/2017: changed def of L and R to make Y pol have same sign as H
        elif "rcvrPol" in line :
          if 'H' in line[line.find("=")+1:].upper() :
            Ei = numpy.array( [0.,1.], dtype=complex )
          elif 'V' in line[line.find("=")+1:].upper() :
            Ei = numpy.array( [1.,0.], dtype=complex )
          elif 'L' in line[line.find("=")+1:].upper() :
            #Ei = numpy.array( [1/math.sqrt(2.),1j/math.sqrt(2)], dtype=complex )
            Ei = numpy.array( [-1j/math.sqrt(2.),1./math.sqrt(2)], dtype=complex )
          elif 'R' in line[line.find("=")+1:].upper() :
            #Ei = numpy.array( [1/math.sqrt(2.),-1j/math.sqrt(2)], dtype=complex )
            Ei = numpy.array( [1j/math.sqrt(2.),1/math.sqrt(2)], dtype=complex )
          rcvrPolMissing = False
    fin.close()
    if LOGHzMissing or rhoOffsetMissing or rcvrPolMissing :
      print "... saphReadData WARNING - missing parameter" 
    fout = open( logfile, "a")
    fout.write( "# ... read data from %s -- npts = %d, LOGHz = %.1f, rhoOffset = %.1f deg, Lfrac = %.2f, Ei = [%.2f%+.2fj, %.2f%+.2fj]\n" % \
       (infile, len(P3), LOGHz, rhoOffset, Lfrac, Ei[0].real, Ei[0].imag, Ei[1].real, Ei[1].imag))
    fout.close()
    print " ... read data from %s -- npts = %d, LOGHz = %.1f, rhoOffset = %.1f deg, Lfrac = %.2f, Ei = [%.2f%+.2fj, %.2f%+.2fj]" % \
       (infile, len(P3), LOGHz, rhoOffset, Lfrac, Ei[0].real, Ei[0].imag, Ei[1].real, Ei[1].imag)
    return numpy.array(fGHz), numpy.array(P3), numpy.array(P4), numpy.array(P5), numpy.array(P6), LOGHz, rhoOffset, Lfrac, Ei 

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

    
# saphFit2 searches for best fit for thermal emission data (jun and sep 2015),
#    which are input temperature vs IF frequency measured with DSB receivers
# it is designed to fit multiple data sets simultaneously (e.g, taken at different
#    LO frequencies, different stack rotation angles, and different polarizations)
#    these are given by dataFileList
# the LO frequency, the stack rotation offset, and the input polarization for each
#    data set within dataFileList are read from comment lines in the data file
# input data are smoothed using spline fit so that smaller number of freqs can be 
#    modeled
# search ranges are specified with srchList, which contains a list of values to be
#    checked for each of the parameters - note that this is just a list, not the
#    inputs to numpy.arange
# other variables:
# ... T5fit to find min variance in T5 (normally fits T3 and T4)
# ... T6override to manually set temperature seen by reflection off the plate; this is important only for T5fit=True

# tcmSrch = numpy.array( [ 1.007, 1.009, 1.011, 1.013])                       
# tcmSrch = numpy.array( [1.009] )
# tcmSrch = numpy.array( [.3750,.3755,.3760] )
# tcmSrch = numpy.array( [.3750] )
tcmSrch = numpy.arange( 1.0085,1.0101,.0005)

# noSrch = [3.020, 3.025, 3.030, 3.035, 3.040]             # ordinary index
noSrch = numpy.arange(3.063,3.0671,.001)             # ordinary index
# noSrch = numpy.array( [3.065] )

# neSrch = [3.380, 3.385, 3.390, 3.395, 3.400]             # extraordinary index
neSrch = numpy.arange(3.402,3.4061,.001)             # extraordinary index
# neSrch = numpy.array( [3.405] )

# oltSrch = numpy.arange(1.5e-4, 5.1e-4, 0.5e-4)
oltSrch = numpy.array( [3.0e-4 ] )          # ordinary loss tangent
# eltSrch = numpy.arange(1.5e-4, 5.1e-4, 0.5e-4)
eltSrch = numpy.array( [3.0e-4 ] )          # extraordinary loss tangent 

# angISrch = [ 40., 40.5, 41.,41.5, 42., 42.5 ]     # angle of incidence ('stackTilt')
# angISrch = numpy.arange( 42.,48.1,1.)
angISrch = numpy.array( [ 43.,44.,45.,46.,47. ] )     # angle of incidence ('stackTilt')
# angISrch = numpy.array( [ 44.] )     # angle of incidence ('stackTilt')

rhoSrch = numpy.arange(-90.,88.1,4.)              # stack rotation angle relative to plane of incidence ('stackAngle')
#rhoSrch = numpy.array( [38.] )

srchList = [ noSrch, neSrch, oltSrch, eltSrch, tcmSrch, angISrch, rhoSrch ]

#dataFileList = ["100_0deg.dat", "100_30deg.dat", "100_60deg.dat", "100_90deg.dat", \
#                "110_0deg.dat", "110_30deg.dat", "110_60deg.dat", "110_90deg.dat" ]
#                "220_0deg.dat", "220_30deg.dat", "220_60deg.dat", "220_90deg.dat", \
#                "250_0deg.dat", "250_30deg.dat", "250_60deg.dat", "250_90deg.dat" ]

dataFileList = [ "250_0deg.dat", "250_30deg.dat", "250_60deg.dat", "250_90deg.dat" ]
#dataFileList = [ "220_0deg.dat", "220_30deg.dat", "220_60deg.dat", "220_90deg.dat" ]
#dataFileList = [ "110_0deg.dat", "110_30deg.dat", "110_60deg.dat", "110_90deg.dat" ]
#dataFileList = [ "100_0deg.dat", "100_30deg.dat", "100_60deg.dat", "100_90deg.dat" ]
#dataFileList = [ "220_pos10.dat" ]
#dataFileList = [ "229.5_pos12.dat" ]

# doit6 makes it quick to run saphFit2 without remembering the input variables
def doit6() :
    saphFit2( dataFileList, srchList)

# does forward modeling of data given in dataFileList, using parameters given in srchList;
#   computes variances for each set of parameters, writes results to outFile
# when finished, writes bestfit parameters for lowest var34 and var5 to outFile,
#   and plots results for lowest var34 on screen
#
def saphFit2( dataFileList, srchList=srchList, T6override=None, outFile="saphFit2.log" ) :
    fout = open( outFile, "a")
    fout.write("#\n# %s\n" % datetime.datetime.now().strftime("%m/%d/%Y %H:%M") )
    fout.write("# dataFileList: %s\n" % (dataFileList) )
    fout.write("# tcmSrch:  %s\n" % numpy.array_str(tcmSrch, precision=4, max_line_width=200 ) )
    fout.write("# noSrch:   %s\n" % numpy.array_str(noSrch, precision=3, max_line_width=200 ) )
    fout.write("# neSrch:   %s\n" % numpy.array_str(neSrch, precision=3, max_line_width=200 ) )
    fout.write("# oltSrch:  %s\n" % numpy.array_str(oltSrch, precision=5, max_line_width=200 ) )
    fout.write("# eltSrch:  %s\n" % numpy.array_str(eltSrch, precision=5, max_line_width=200 ) )
    fout.write("# angISrch: %s\n" % numpy.array_str(angISrch, precision=2, max_line_width=200 ) )
    fout.write("# rhoSrch:  %s\n" % numpy.array_str(rhoSrch, precision=2, max_line_width=200 ) )
    fout.write("#\n")
    fout.close()

    frqArray = numpy.array( numpy.arange(.3,9.91,.4) )

    ntot = len(srchList[0]) * len(srchList[1]) * len(srchList[2]) * len(srchList[3]) \
             * len(srchList[4]) * len(srchList[5]) * len(srchList[6])

  # search for fit that gives min variance if ntot > 1
    if ntot > 1 :
      print "search for best fit -- modeling %d possibilities" % ntot
      print "read and smooth the data"
      T3sm = []       # list of arrays
      T4sm = []       # list of arrays
      T5sm = []       # list of arrays
      T6avg = []      # list of numbers
      LOfrq = []      # list of numbers
      rhoOff = []     # list of numbers
      Ei = []         # list of vectors
      Lfrac = []	    # list of numbers (lower sideband power fraction; 0.5 should be normal)
      for dataFile in dataFileList :  
        fGHz,T3,T4,T5,T6,LOGHz,rhoOffset,Lfraction,Eivec = saphReadData( dataFile )
        T3sm.append( saphSmooth( fGHz, T3, frqArray ) )
        T4sm.append( saphSmooth( fGHz, T4, frqArray ) )
        T5sm.append( saphSmooth( fGHz, T5, frqArray ) )
        LOfrq.append( LOGHz )
        rhoOff.append( rhoOffset )
        Ei.append( Eivec )
        Lfrac.append( Lfraction )
      if T6override :
        T6avg.append = T6override
      else :
        T6avg.append( numpy.mean( T6 ) )
      Tsm = numpy.array( [T3sm,T4sm,T5sm] )
    
      secs0 = time.time()
      nfinished = 0
      var5min = var34min = 1.e6
      for no in srchList[0] :
        for ne in srchList[1] :
          for olt in srchList[2] :
            for elt in srchList[3] :
              for tcm in srchList[4] :
                for angI in srchList[5] :
                  for rho in srchList[6] :

                  # form the stackDesc arrays (one value per layer, in this case) needed by saphModel
                    stackDesc = [ [no], [ne], [olt], [elt], [tcm], [0.] ]    

                  # model each dataFile separately (same stack, but LOfrq, rhoOff, or polAng could be different)
                    var34 = []
                    var5 = []
                    for n in range(0,len(dataFileList)) :
                      TmL = saphModel1( LOfrq[n]-frqArray, stackDesc, rho-rhoOff[n], angI, Ei[n] ) 
                      TmU = saphModel1( LOfrq[n]+frqArray, stackDesc, rho-rhoOff[n], angI, Ei[n] ) 
                      Tm,Lfbest,var34DSB,var5DSB = saphDSB( TmL, TmU, Tsm, Lfrac[n] )
                      var34.append( var34DSB ) 
                      var5.append( var5DSB )
                      print "... %20s   var34: %.2f   var5: %.2f   Lfrac: %.2f " % (dataFileList[n], var34[-1], var5[-1], Lfbest)
                      fout = open(outFile, "a")
                      fout.write("...  %20s   var34: %.2f   var5: %.2f  Lfbest: %.2f\n" % (dataFileList[n], var34[-1], var5[-1], Lfbest))
                      fout.close()
                    nfinished = nfinished + 1
                    secselapsed = time.time() - secs0
                    secsremaining = (ntot-nfinished)*(secselapsed/nfinished)
                    if secsremaining/60. > 60. :
                      print "%d/%d finished; %.2f hours remaining" % (nfinished,ntot,secsremaining/3600.)
                    else :
                      print "%d/%d finished; %.2f minutes remaining" % (nfinished,ntot,secsremaining/60.)

                  # var34sum and var5sum are summed over the independent data files
                    var34sum = numpy.sum(var34)
                    var5sum = numpy.sum(var5)
                    print "tcm = %.4f, no = %.3f, ne = %.3f, angI = %.1f, rho = %.1f, olt,elt = %.2e, %.2e, var34 = %.1f, var5 = %.2f\n" % \
                        (tcm, no, ne, angI, rho, olt, elt, var34sum, var5sum )
                    fout = open(outFile, "a")
                    fout.write( "tcm: %.4f  nO: %.3f  nE: %.3f  angI: %.1f  rho: %.1f  olt,elt: %.2e %.2e  var34: %.1f  var5: %.2f\n" % \
                       (tcm, no, ne, angI, rho, olt, elt, var34sum, var5sum ))
                    fout.close()
                    if (var34sum < var34min) :
                      stackDesc34 = stackDesc
                      rho34 = rho
                      angI34 = angI
                      var34min = var34sum
                    if (var5sum < var5min) :
                      stackDesc5 = stackDesc
                      rho5 = rho
                      angI5 = angI
                      var5min = var5sum
  
    # when all possibilities have been modeled, print out best parameters
      print "# lowest var34: tcm = %.3f, no = %.3f, ne = %.3f, angI = %.1f, rho = %.1f, olt,elt = %.2e, %.2e, var34 = %.1f" % \
         ( stackDesc34[4][0], stackDesc34[0][0], stackDesc34[1][0], angI34, rho34, stackDesc34[2][0], stackDesc34[3][0], var34min ) 
      print "# lowest var5 : tcm = %.3f, no = %.3f, ne = %.3f, angI = %.1f, rho = %.1f, olt,elt = %.2e, %.2e, var5 = %.1f" % \
         ( stackDesc5[4][0], stackDesc5[0][0], stackDesc5[1][0], angI5, rho5, stackDesc5[2][0], stackDesc5[3][0], var5min ) 
      fout = open( outFile, "a")
      fout.write("# lowest var34: tcm: %.4f  nO: %.3f  nE: %.3f  angI: %.1f  rho: %.1f  olt,elt: %.2e %.2e  var34: %.1f\n" % \
         ( stackDesc34[4][0], stackDesc34[0][0], stackDesc34[1][0], angI34, rho34, stackDesc34[2][0], stackDesc34[3][0], var34min ))
      fout.write("# lowest var5:  tcm: %.4f  nO: %.3f  nE: %.3f  angI: %.1f  rho: %.1f  olt,elt: %.2e %.2e  var5: %.1f\n" % \
         ( stackDesc5[4][0], stackDesc5[0][0], stackDesc5[1][0], angI5, rho5, stackDesc5[2][0], stackDesc5[3][0], var5min ) )
      fout.close()

  # if there is only one possibility, skip optimization, go directly to plot
    else :
      stackDesc34 = [ srchList[0], srchList[1], srchList[2], srchList[3], srchList[4], [0.] ]
      angI34 = srchList[5][0]
      rho34 = srchList[6][0]

  # plot best var34 fit 
    plotFit( dataFileList, stackDesc34, rho34, angI34  ) 


# converts alpha to tanDelta, using average index of 3.25
def tanD (fGHz, alpha) :
  print "tanDelta = %.2e" %  (alpha * clight/ (2.*math.pi*3.25*fGHz) )

# plot saphFit2 model for data files in dataFileList
# model is specified by stackDesc, rho, angI
# if showSB, plot USB and LSB values (each normalized to 1) with faint dotted lines
#
def plotFit( dataFileList, stackDesc, rho, angI, showSB=False ) :
    frqArray = numpy.array( numpy.arange(.3,9.91,.2) )
    pyplot.ioff()
    pp = PdfPages("BestFit.pdf")
    pyplot.figure( figsize=(11,8) )
    nplot = 1
    nmax = len( dataFileList)
    nm = 1
    fsize = 12
    if nmax > 1 :
      nm = 2 
      fsize = 6
    for n in range(0,len(dataFileList)) :
      fGHz,T3,T4,T5,T6,LOGHz,rhoOffset,Lfraction,Eivec = saphReadData( dataFileList[n] )
         # reread each data file so we can plot UNSMOOTHED data; smoothed data already in T3sm, T4sm, T5sm arrays
      T3sm = saphSmooth( fGHz, T3, frqArray )
      T4sm = saphSmooth( fGHz, T4, frqArray )
      T5sm = saphSmooth( fGHz, T5, frqArray )
      Tsm = numpy.array( [T3sm,T4sm,T5sm] )
      TmL = saphModel1( LOGHz-frqArray, stackDesc, rho-rhoOffset, angI, Eivec ) 
      TmU = saphModel1( LOGHz+frqArray, stackDesc, rho-rhoOffset, angI, Eivec ) 
      Tm,Lfbest,var34DSB,var5DSB = saphDSB( TmL, TmU, Tsm, Lfraction )
      
      fig = pyplot.subplot(nm,nm,nplot)
      fig.plot( fGHz, T3, color='red' )
      fig.plot( frqArray, Tm[0], color='red' )
      if showSB :
        fig.plot( frqArray, TmU[0], color='red', linestyle='dashed', lw=0.5 )
        fig.plot( frqArray, TmL[0], color='red', linestyle='dotted', lw=0.5 )
      fig.plot( fGHz, T4, color='blue' )
      fig.plot( frqArray, Tm[1], color='blue' )
      if showSB :
        fig.plot( frqArray, TmU[1], color='blue', linestyle='dashed', lw=0.5 )
        fig.plot( frqArray, TmL[1], color='blue', linestyle='dotted', lw=0.5 )
      fig.plot( fGHz, T5, color='purple' )
      fig.plot( frqArray, Tm[2], color='purple' )
      fig.grid( True )
      fig.set_title( dataFileList[n], fontsize=12 )
      fig.set_ylim(0,300)
      fig.text( .05, .1, "tcm = %.4f, no = %.3f, ne = %.3f, var34 = %.1f, var5 = %.1f" % \
          ( stackDesc[4][0], stackDesc[0][0], stackDesc[1][0], var34DSB, var5DSB ), \
          transform=fig.transAxes, horizontalalignment="left", verticalalignment="center", fontsize=fsize )
      fig.text( .05, .05, "rho = %.1f, angI = %.1f, olt = %.2e, elt = %.2e, Lfrac = %.2f" % \
          ( rho, angI, stackDesc[2][0], stackDesc[3][0], Lfbest), \
          transform=fig.transAxes, horizontalalignment="left", verticalalignment="center", fontsize=fsize )
      nplot = nplot + 1
      if nplot == 5 :
        pyplot.savefig( pp, format="pdf" )
        print "plotting data on screen; delete plot window to continue"
        pyplot.show()
        pyplot.figure( figsize=(11,8) )
        nplot = 1
    if (nplot > 1) :
      pyplot.savefig( pp, format="pdf" )
      print "plotting data on screen; delete plot window to finish"
      pyplot.show()
    pp.close()
    print "BestFit.pdf written to disk"


# reads fit parameters from logFile (e.g., saphFit2.log) and plots model curves
#   vs data for files in dataFileList
# uses fitReadLog to retrieve the fit parameters; these will be taken from last
#   line beginning with "*tcm:" if such a line exists, otherwise from line with
#   smallest var34 
# note 1: last parameter returned by fitReadLog is an ascii string with the
#   dataFileList used by the fitting routine; this parameter is ignored by
#   plotFinal, which plots the data in the list of files given by dataFileList;
#   this latter is a true list,not an ascii string; also, it makes it possible
#   to test fit parameters from one plot to fit other data files
# note 2: fitReadLog returns lists of variables no, ne, etc; unfortunately
#   these are NOT the same as the lists required by anisotropicTransmissionPol;
#   in the first case, no[i] is the ordinary index used for a particular trial
#   fit; in the second case, no[i] is the ordinary index of the ith layer
#
def plotFinal( dataFileList=dataFileList, resultsLog="saphFit2.log", plot5=False, showSB=False ) :
    tcm,no,ne,angI,rho,olt,elt,var34,var5,iplot,fitList = fitReadLog( resultsLog, plot5 )
    stackDesc = [ [no[iplot]], [ne[iplot]], [olt[iplot]], [elt[iplot]], [tcm[iplot]], [0.] ]    
    print stackDesc, rho[iplot], angI[iplot]
    plotFit( dataFileList, stackDesc, rho[iplot], angI[iplot], showSB=showSB ) 
        

# read in results from saphFit2.log
# returns vectors of fit values; also iplot and dataFileList
# iplot is the index of the value that should be plotted:
#   - if plot5 == True, plot line corresponding to min(var5)
#   - else if file contains lines beginning with "*tcm:", iplot is index of LAST of these
#   - default is to plot line corresponding to min(var34)
#   - 
def fitReadLog( infile, plot5 ) :
  fin = open( infile, "r" )
  tcm = []
  no = []
  ne = []
  angI = []
  rho = []
  olt = []
  elt = []
  var34 = []
  var5 = []
  iplot = -1
  dataFileList = ""
  for line in fin :
    if line.startswith("tcm:") or line.startswith("*tcm:") :
      a = line.split()
      tcm.append( float(a[1]) )
      no.append( float(a[3]) )
      ne.append( float(a[5]) )
      angI.append( float(a[7]) )
      rho.append( float(a[9]) )
      olt.append( float(a[11]) )
      elt.append( float(a[12]) )
      var34.append( float(a[14]) )
      var5.append( float(a[16]) )
    if line.startswith("*tcm:") :
      iplot = len(tcm) - 1         # line number of vector 
    if line.startswith("# dataFileList:") :    
      dataFileList = line[ line.find("[")+1 : line.find("]")-1 ]
  fin.close()
  if plot5 :
    iplot=var5.index( min(var5) )
  elif (iplot == -1) :
    iplot = var34.index( min(var34) )
  print "read in %d lines from file %s" % (len(tcm),infile)
  print var34[iplot]
  return tcm,no,ne,angI,rho,olt,elt,var34,var5,iplot,dataFileList
    # first 9 items are lists, last 2 are values

# show correlations between params
def covarPlot( infile='saphFit2.log', plot5=False, vmin=0., vmax=0. ) :
    label = ['tcm','no','ne','angI','rho','olt','elt','var34','var5']
    vec = fitReadLog( infile, plot5 )
    ibest = vec[9]
       # ibest is the index of the point through which all planes will pass
    if plot5 :
      jz = 8      # var5
    else :
      jz = 7      # var34
    pyplot.ioff()
    pp = PdfPages("Covar.pdf")
    pyplot.figure( figsize=(11,8) )
    nplot = 1
    for jx in range(0,7) :
      for jy in range( jx+1,7 ) :
        x,y,z = stripout( vec, ibest, jx, jy, jz )
        jbest = numpy.argmin(z)
        if (len(numpy.unique(x)) > 1) and (len(numpy.unique(y)) > 1) :
          fig = pyplot.subplot(4,3,nplot)
          xmin,xmax = minmax(x,margin=0.01)
          ymin,ymax = minmax(y,margin=0.01)
          if (vmin == 0.) and (vmax == 0.) :
            vmin,vmax = minmax(z,margin=0.0)
        # section below copied from web example
          xi,yi = numpy.linspace( xmin, xmax, 50), numpy.linspace(ymin, ymax, 50)
          xi,yi = numpy.meshgrid(xi,yi)
          zi = scipy.interpolate.griddata( (x,y) ,z, (xi,yi), method='cubic')
          xibest,yibest = numpy.unravel_index(numpy.nanargmin(zi),zi.shape)
              # return position of interpolated minimum
              # zi[xibest,yibest] occurs at x[xibest,yibest],y[xibest,yibest]
          # print "xibest,yibest,zibest =",xibest,yibest,zi[xibest,yibest]
          fig.imshow( zi, aspect='auto', origin='lower', vmin=vmin, vmax=vmax, \
              extent=[xmin,xmax,ymin,ymax], cmap='terrain' )
          fig.scatter(x,y,s=10, marker='x', c='black')
          fig.scatter(x[jbest],y[jbest],s=20, marker='x', c='white')
          fig.scatter(xi[xibest,yibest],yi[xibest,yibest],s=40, marker='s', c='red')
          fig.ticklabel_format(useOffset=False)
          fig.set_xlabel( label[jx], fontsize=6 )
          fig.set_ylabel( label[jy], fontsize=6 )
          #fig.scatter(x,y,c=z,s=200,edgecolors='none',marker='s',vmin=0.,vmax=.3*z.max())
          #fig.scatter(x,y,c=z,s=200,edgecolors='none',marker='s',vmin=50.,vmax=200.)
          fig.axis( [xmin, xmax, ymin, ymax],fontsize=6 )
          pyplot.locator_params( axis='x', nbins=5)
          pyplot.locator_params( axis='y', nbins=5)
          fig.tick_params( axis='both', which='major', labelsize=6 )
          pyplot.tight_layout( pad=3, h_pad=1)

        # increment the subplot; begin new page if too many panels
          nplot = nplot + 1
          if nplot > 12 :
            pyplot.savefig( pp, format="pdf" )
            pyplot.show()
            pyplot.figure( figsize=(11,8) )
            nplot = 1

  # final panel gives point through which planes pass
    fig = pyplot.subplot(4,3,nplot)
    fig.text( 0.03, .92, "dataFileList: %s" % vec[10], \
        transform=fig.transAxes, horizontalalignment='left', fontsize=8, \
        color="black", rotation='horizontal' )
    ylab = .82 
    fig.text( 0.03, ylab, "planes pass through the following point: (record #%d)" % ibest, \
      transform=fig.transAxes, horizontalalignment='left', fontsize=8, \
      color="black", rotation='horizontal' )
    for jx in range(0,8) :
      ylab = ylab - .1
      fig.text( 0.04, ylab, "%s = %.4f" % (label[jx],vec[jx][ibest]), \
        transform=fig.transAxes, horizontalalignment='left', fontsize=8, \
        color="black", rotation='horizontal' )
      pyplot.axis('off')
    if (nplot > 1) :
      pyplot.savefig( pp, format="pdf" )
      pyplot.show()
    pp.close()
    print "Covar.pdf written to disk"

# strip out plane corresponding to variables jx,jy from a multidimensional data cube, passing through point with lowest variance
# that is, identify rows in vec where all variables other than jx,jy are held at their optimum value
def stripout( vec, ibest, jx, jy, jz ) :
    print "\nplot plane passing through %.5f, %.5f (var = %.5f)" % (vec[jx][ibest],vec[jy][ibest],vec[jz][ibest])
    x = []
    y = []
    z = []
    for i in range(0,len(vec[0])) :  
      usePt = True         
      for j in range(0,4) :
      # all variables except jx and jy must be equal to best values
        if (j != jx) and (j != jy) and (vec[j][i] != vec[j][ibest] ) :
          usePt = False
      if usePt :
        print vec[0][i], vec[1][i], vec[2][i], vec[3][i], vec[4][i], vec[5][i], vec[6][i], vec[jz][i]
        x1 = vec[jx][i]
        y1 = vec[jy][i]
        z1 = vec[jz][i]
        x.append(x1)
        y.append(y1)
        z.append(z1)
    return numpy.array(x),numpy.array(y),numpy.array(z)

# quickly compute expected relationship between n, tcm, and angI
def sim() :
   n0 = 3.06
   tcm0 = 1.009
   lambda0 = 30./250.
   angI = numpy.arange(40.,50.,1.) 
   theta2 =  numpy.arcsin(numpy.sin( numpy.radians(angI) )/n0)  
   print angI
   print numpy.degrees(theta2)
   dphi = n0 * tcm0 / (lambda0 * numpy.cos(theta2))
   tcmBest = lambda0 * dphi[5] * numpy.cos(theta2) / n0
   print tcmBest
   
# compare amp and phase for different indices and thicknesses; can I measure the difference?
def checkSensitivity() :
    fGHzArr = numpy.arange(70.,230.5,.5)
    #fGHzArr = numpy.arange(210.,230.5,0.1)
    angIdeg = 8.
    tcmArr1 = 2.54 * numpy.array( [.008,.009,.250,.009,.008] )
    tcmArr2 = 2.54 * numpy.array( [.009,.008,.250,.008,.009] )
    nrArr1 = numpy.array( [1.45, 2.42, 3.114, 2.42, 1.45] )
    nrArr2 = nrArr1 
    tanDeltaArr = numpy.array( [0.] ) 
    phs1,trans1 = solveStack(fGHzArr, angIdeg, tcmArr1, nrArr1, tanDeltaArr ) 
    phs2,trans2 = solveStack(fGHzArr, angIdeg, tcmArr2, nrArr2, tanDeltaArr ) 
    print phs1
    print phs2
    print trans1
    print trans2

  # begin the plot
    pyplot.ioff()
    pp = PdfPages("testSens.pdf")
    fig =  pyplot.figure( figsize=(11,8) )
    #pyplot.suptitle( "%s (%s)" % (title,infile), fontsize=14 )
    dphi = unwrap( phs1-phs2 )

  # top left panel is phase vs freq, measured and fit
    ax = fig.add_axes( [.1,.6,.38,.32] )
    xmin,xmax = minmax( fGHzArr )
    ymin,ymax = minmax( numpy.concatenate((phs1,phs2)))
    ax.axis( [xmin, xmax, ymin, ymax] )
    ax.tick_params( labelsize=10 )
    ax.plot( fGHzArr, phs1, "b-" )
    ax.plot( fGHzArr, phs2, "r-" )
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
    ax.plot( fGHzArr, dphi, "b-" )
    pyplot.grid(True)
                                                       
  # top right panel is transmission vs freq, measured and theoretical
    ax = fig.add_axes( [.57,.6,.38,.32] )
    ymin,ymax = minmax( numpy.concatenate((trans1,trans2)) )
    ax.axis( [xmin, xmax, ymin, ymax] )
    ax.plot( fGHzArr, trans1, "-", color='blue' )
    ax.plot( fGHzArr, trans2, "r-" )
    ax.tick_params( labelsize=10 )
    ax.set_ylabel("transmission", fontsize=10)
    pyplot.grid(True)

  # middle right panel is measured-theoretical trans
    ax = fig.add_axes( [.57,.37,.38,.2])
    ax.tick_params( labelsize=10 )
    ymin,ymax = minmax( trans1-trans2, margin=0.2 )
    ax.axis( [xmin, xmax, ymin, ymax] )
    ax.plot( fGHzArr, trans1-trans2, "-", color='blue' )
    ax.set_xlabel("freq (GHz)", fontsize=10)
    ax.set_ylabel("residual", fontsize=10)
    pyplot.grid(True)

  # fictitious lower right panel gives room for text
  #  ax = fig.add_axes( [.57,.1,.38,.2])
  #  ylab = 0.84
  #  for nl in range(0,nlayers) :
  #    ax.text( 0.03, ylab, "layer %d:  t = %.5f in  nr = %.3f" % \
  #      ( nl+1, tcmList[nbest][nl]/2.54, nrList[nbest][nl]), \
  #      transform=ax.transAxes, horizontalalignment='left', fontsize=14, \
  #      color="black", rotation='horizontal' )
  #    ylab = ylab - 0.14
  #  #ax.text( 0.03, .56, "angle = %.1f deg" % angIdeg, transform=ax.transAxes, \
  #  #  horizontalalignment='left', fontsize=14, color="black", rotation='horizontal' )
  #  #ax.text( 0.03, .42, "tanDelta = %.1e" % tanDelta[0], transform=ax.transAxes, \
  #  #  horizontalalignment='left', fontsize=14, color="black", rotation='horizontal' )
  #  pyplot.axis('off')
    pyplot.savefig( pp, format="pdf" )
    pyplot.show()
    pp.close()

def checkSensitivity2() :
    # fGHzArr = numpy.append( numpy.arange(75.,115.5,.1), numpy.arange(210.,230.5,0.1) )
    fGHzArr = numpy.arange(70.,115.5,.5)
    fGHzArr = numpy.arange(220.,230.5,.2)
    angIdeg = 8.
    tcmArr1 = 2.54 * numpy.array( [.009,.2467] )
    tcmArr2 = 2.54 * numpy.array( [.0088,.2467] )
    nrArr1 = numpy.array( [2.49, 3.114] )
    nrArr2 = numpy.array( [2.49, 3.114] )
    tanDeltaArr = numpy.array( [0.,0.] ) 
    phs1,trans1 = solveStack(fGHzArr, angIdeg, tcmArr1, nrArr1, tanDeltaArr ) 
    phs2,trans2 = solveStack(fGHzArr, angIdeg, tcmArr2, nrArr2, tanDeltaArr ) 
    print phs1
    print phs2
    print trans1
    print trans2

  # begin the plot
    pyplot.ioff()
    pp = PdfPages("testSens.pdf")
    fig =  pyplot.figure( figsize=(11,8) )
    #pyplot.suptitle( "%s (%s)" % (title,infile), fontsize=14 )
    dphi = unwrap( phs1-phs2 )

  # top left panel is phase vs freq, measured and fit
    ax = fig.add_axes( [.1,.6,.38,.32] )
    xmin,xmax = minmax( fGHzArr )
    ymin,ymax = minmax( numpy.concatenate((phs1,phs2)))
    ax.axis( [xmin, xmax, ymin, ymax] )
    ax.tick_params( labelsize=10 )
    ax.plot( fGHzArr, phs1, "b-" )
    ax.plot( fGHzArr, phs2, "r-" )

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
    ax.plot( fGHzArr, dphi, "b-" )
    pyplot.grid(True)
                                                       
  # top right panel is transmission vs freq, measured and theoretical
    ax = fig.add_axes( [.57,.6,.38,.32] )
    ymin,ymax = minmax( numpy.concatenate((trans1,trans2)) )
    ax.axis( [xmin, xmax, ymin, ymax] )
    ax.plot( fGHzArr, trans1, "-", color='blue' )
    ax.plot( fGHzArr, trans2, "r-" )
    ax.tick_params( labelsize=10 )
    ax.set_ylabel("transmission", fontsize=10)
    pyplot.grid(True)

  # middle right panel is measured-theoretical trans
    ax = fig.add_axes( [.57,.37,.38,.2])
    ax.tick_params( labelsize=10 )
    ymin,ymax = minmax( trans1-trans2, margin=0.2 )
    ax.axis( [xmin, xmax, ymin, ymax] )
    ax.plot( fGHzArr, trans1-trans2, "-", color='blue' )
    ax.set_xlabel("freq (GHz)", fontsize=10)
    ax.set_ylabel("residual", fontsize=10)
    pyplot.grid(True)

  # fictitious lower right panel gives room for text
  #  ax = fig.add_axes( [.57,.1,.38,.2])
  #  ylab = 0.84
  #  for nl in range(0,nlayers) :
  #    ax.text( 0.03, ylab, "layer %d:  t = %.5f in  nr = %.3f" % \
  #      ( nl+1, tcmList[nbest][nl]/2.54, nrList[nbest][nl]), \
  #      transform=ax.transAxes, horizontalalignment='left', fontsize=14, \
  #      color="black", rotation='horizontal' )
  #    ylab = ylab - 0.14
  #  #ax.text( 0.03, .56, "angle = %.1f deg" % angIdeg, transform=ax.transAxes, \
  #  #  horizontalalignment='left', fontsize=14, color="black", rotation='horizontal' )
  #  #ax.text( 0.03, .42, "tanDelta = %.1e" % tanDelta[0], transform=ax.transAxes, \
  #  #  horizontalalignment='left', fontsize=14, color="black", rotation='horizontal' )
  #  pyplot.axis('off')
    pyplot.savefig( pp, format="pdf" )
    pyplot.show()
    pp.close()

# figure out xband freqs, if any, that will allow me to phaselock one gunn osc
#   at approx 105-115 GHz (for 2nd harmonic mixer), and the other at 70-77 GHz
#   (for tripler used as transmitter); for 1mm optics bench measurements
def quickFind() :
    for xGHz in numpy.arange(8.0,12.5,.01) :
      for n in range(0,20) :
        for m in range(0,20) :
          f1 = n*xGHz + .08
          f2 = m*xGHz + .04
          if (f1 > 100.) and (f1 < 115.) and (f2 > 67.) and (f2 < 77.) and \
             ( abs(3.*f2 - 2.*f1) < 1.) :
            print xGHz, n, f1, m, f2, 1000.*(2.*f1-3.*f2)

# show correlations between params for refractive index fits - specialized to just 2 parameters
def covarPlot2( infile, jopt=7 ) :
    fin = open( infile, "r" )
    [tcmList,nrList,avgPhResid,stdPhResid,stdAmpResid,Resid] = pickle.load( fin )

  # create array of results with 7 columns (t1,t2,nr1,nr2,avgPh,stdPh,stdAmp,Resid)
    vec = numpy.concatenate( ( [tcmList[:,0]/2.54], [tcmList[:,1]/2.54], [nrList[:,0]], [nrList[:,1]], \
      [avgPhResid], [stdPhResid], [stdAmpResid], [Resid] ), axis=0 )
    label = ["t1", "t2", "nr1", "nr2"]
  
  # find index of point with minimum avg or variance (whatever we are using)
    ibest = numpy.argmin( numpy.abs(vec[jopt]) )
    print len(vec[0]), ibest, vec[:,ibest]

    pyplot.ioff()
    pp = PdfPages("Covar.pdf")
    pyplot.figure( figsize=(11,8) )
    nplot = 1
    for jx in range(0,4) :
      for jy in range( jx+1,4 ) :
        x,y,z = stripout( vec, ibest, jx, jy, jopt )
        print "\n"
        jbest = numpy.argmin(z)
        if (len(numpy.unique(x)) > 1) and (len(numpy.unique(y)) > 1) :
          fig = pyplot.subplot(2,3,nplot)
          xmin,xmax = minmax(x,margin=0.01)
          ymin,ymax = minmax(y,margin=0.01)
          vmin,vmax = minmax(z,margin=0.0)
        # section below copied from web example
          xi,yi = numpy.linspace( xmin, xmax, 50), numpy.linspace(ymin, ymax, 50)
          xi,yi = numpy.meshgrid(xi,yi)
          zi = scipy.interpolate.griddata( (x,y) ,z, (xi,yi), method='cubic')
          xibest,yibest = numpy.unravel_index(numpy.nanargmin(zi),zi.shape)
              # return position of interpolated minimum
              # zi[xibest,yibest] occurs at x[xibest,yibest],y[xibest,yibest]
          # print "xibest,yibest,zibest =",xibest,yibest,zi[xibest,yibest]
          fig.imshow( zi, aspect='auto', origin='lower', vmin=vmin, vmax=vmax, \
              extent=[xmin,xmax,ymin,ymax], cmap='terrain' )
          fig.scatter(x,y,s=10, marker='x', c='black')
          fig.scatter(x[jbest],y[jbest],s=20, marker='x', c='white')
          fig.scatter(xi[xibest,yibest],yi[xibest,yibest],s=40, marker='s', c='red')
          fig.ticklabel_format(useOffset=False)
          fig.set_xlabel( label[jx], fontsize=8 )
          fig.set_ylabel( label[jy], fontsize=8 )
          #fig.scatter(x,y,c=z,s=200,edgecolors='none',marker='s',vmin=0.,vmax=.3*z.max())
          #fig.scatter(x,y,c=z,s=200,edgecolors='none',marker='s',vmin=50.,vmax=200.)
          fig.axis( [xmin, xmax, ymin, ymax],fontsize=6 )
          pyplot.locator_params( axis='x', nbins=5)
          pyplot.locator_params( axis='y', nbins=5)
          fig.tick_params( axis='both', which='major', labelsize=8 )
          pyplot.tight_layout( pad=3, h_pad=1)

        # increment the subplot; begin new page if too many panels
          nplot = nplot + 1
          if nplot > 6 :
            pyplot.savefig( pp, format="pdf" )
            pyplot.show()
            pyplot.figure( figsize=(11,8) )
            nplot = 1
    pyplot.savefig( pp, format="pdf" )
    pyplot.show()
    pp.close()

def junk() :
  # final panel gives point through which planes pass
    fig = pyplot.subplot(4,3,nplot)
    fig.text( 0.03, .92, "dataFileList: %s" % vec[10], \
        transform=fig.transAxes, horizontalalignment='left', fontsize=8, \
        color="black", rotation='horizontal' )
    ylab = .82 
    fig.text( 0.03, ylab, "planes pass through the following point: (record #%d)" % ibest, \
      transform=fig.transAxes, horizontalalignment='left', fontsize=8, \
      color="black", rotation='horizontal' )
    for jx in range(0,8) :
      ylab = ylab - .1
      fig.text( 0.04, ylab, "%s = %.4f" % (label[jx],vec[jx][ibest]), \
        transform=fig.transAxes, horizontalalignment='left', fontsize=8, \
        color="black", rotation='horizontal' )
      pyplot.axis('off')
    if (nplot > 1) :
      pyplot.savefig( pp, format="pdf" )
      pyplot.show()
    pp.close()
    print "Covar.pdf written to disk"

# generate an array of perfect fake data
# will input this to nrefracFit2, use it to figure out how to weight amps and phases
def fakeData( outFile ):
    fout = open( outFile, "w" )
    fGHzArr = numpy.arange(70.,115.1,.1)
    angIdeg = 8.
    tcmArr1 = 2.54 * numpy.array( [.008,.250] )
    nrArr1 = numpy.array( [2.0, 3.114] )
    tanDeltaArr = numpy.array( [0.,0.] ) 
    fout.write("# file generated by ob.fakeData()\n")
    fout.write("# -----------------------------\n")
    fout.write("# title = %s\n" % outFile)
    totin = (tcmArr1[0]+tcmArr1[1])/2.54
    fout.write("# totin = %.4f,%.4f\n" % (totin-.001,totin+.001))
    fout.write("# layer 1\n")
    fout.write("## true tin = %.5f\n" % (tcmArr1[0]/2.54) )
    fout.write("## true nr = %.3f\n" % nrArr1[0] )
    fout.write("# tin = %.4f,%.4f,%.4f\n" % (tcmArr1[0]/2.54-.001, tcmArr1[0]/2.54+.001,.0005))
    fout.write("# tanDelta = 0.\n")
    fout.write("# layer 2\n")
    fout.write("## true tin = %.5f\n" % (tcmArr1[1]/2.54) )
    fout.write("## true nr = %.3f\n" % nrArr1[1] )
    fout.write("# tin = %.4f,%.4f,%.4f\n" % (tcmArr1[1]/2.54-.001, tcmArr1[0]/2.54+.001,.0005))
    fout.write("# tanDelta = 0.\n")
    fout.write("# -----------------------------\n")
    phs1,trans1 = solveStack(fGHzArr, angIdeg, tcmArr1, nrArr1, tanDeltaArr ) 
    for n in range(0,len(fGHzArr)) :
      fout.write("%8.2f  %8.6f  %8.3f  0.005\n" % ( fGHzArr[n], trans1[n], phs1[n] ) )
    fout.close()
