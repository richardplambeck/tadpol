# ob2.py 
# these are routines for dealing with the Optics Bench scans of AR coatings
# it was cloned from ob.py on 28 nov 2017; I then deleted all sapphire-related routines
# all dielectrics are assumed to be isotropic!
# electric field is presumed to be perp to plane of incidence
# 

import numpy
import math
import sys
import time
import string
import cPickle as pickle
import matplotlib
import datetime
#matplotlib.use('GTK')
matplotlib.use('GTKAgg')
import matplotlib.pyplot as pyplot
from matplotlib.backends.backend_pdf import PdfPages
#import pol
import scipy
from scipy.interpolate import splrep, splev
import scipy.fftpack
import scipy.stats
import os
import bisect

# scipy.stats.chi2.isf(.046,2.)    computes chisq for 2 degrees of freedom, .046 probability of exceeding


# lambda is alternate way of defining simple function
sin = lambda arg: numpy.sin(arg)
cos = lambda arg: numpy.cos(arg)
inv = lambda Mat: numpy.linalg.inv(Mat)
dr = lambda theta: math.pi*theta/180.
pow = lambda arg, index: numpy.power(arg, index)
sq = lambda arg: pow(arg, 2)


clight = 29.9792458 # speed of light, cm/nanosec

# makeList is a helper routine for readTransFile
# it interprets an input string (whatever is after the "=" sign on an input line),
#    converting it to a list
# if input string is a const, result is [ x ]
# if input string is a triplet x,y,z, result is [ numpy.arange(x,y,z) ]
def makeList( inString ) :
  try :
    a = inString.split(",")
    if len(a) == 1 :
      return numpy.array( [ float(a[0]) ] )
    elif len(a) == 2 :
      return numpy.array( [float(a[0]), float(a[1]) ] )
    elif len(a) == 3 :
      return numpy.array(numpy.arange( float(a[0]), float(a[1]), float(a[2]) ))
    else :
      print "makeList: error interpreting input string '%s'" % inString
      return [0.]
  except :
      print "makeList: error interpreting input string '%s'" % inString
      return [0.]

# readTransFile inputs data from transmission file created by OpticsBench/09nov2016/unpack.py;
#   or from files averaged with binit, which have rms uncertainty in column 4; can also be used
#   to read Oliver's FTS data if it is reformatted with phase=0. in column 3
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
# 9jan2018 - modify so that this file will read rms uncertainty from averaged data files;
#   use readUnc=True to do this; readUnc=False for older data files
#
def readTransFile( infile, readUnc=True ) :
  fin = open(infile, "r")
  fGHz = []
  trans = []
  dphi = []
  unc = []
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
      if readUnc :                    # averaged data file; data must be preceded by "# readUnc" line
        unc.append( float(a[3]) )        
      else :                          # otherwise set all uncertainties to 1.
        unc.append( 1. )              
    else :                               # lines with single hash marks are search parameters
      if "totcm" in line :
        totcmRange = makeList( line[line.find("=")+1:] )
      elif "totin" in line :
        totcmRange = 2.54*makeList( line[line.find("=")+1:] )
      elif "angIdeg" in line :
        angIdeg = float( line[line.find("=")+1:] )
      elif "tcm" in line :
        tcmList.append( makeList(line[line.find("=")+1:]) )
      elif "tin" in line :
        tcmList.append( makeList(line[line.find("=")+1:]) )
        if tcmList[-1][-1] > 0. :
          tcmList[-1] = 2.54*tcmList[-1]			# apply factor of 2.54 only for unmirrored layers!
      elif "tanDelta" in line :
        tanDeltaList.append( makeList(line[line.find("=")+1:]) )
    # this is an alternative way of specifying mirrored layer
      elif "mirror" in line :
        tcmList.append( makeList(line[line.find("=")+1:]) )
        nrList.append( makeList(line[line.find("=")+1:]) )
        tanDeltaList.append( makeList(line[line.find("=")+1:]) )
      elif "nr" in line :
        nrList.append( makeList(line[line.find("=")+1:]) )
      elif "title" in line :
        title = line[line.find("=")+1:]
        title = title[2:-2]
      elif "readUnc" in line :
        readUnc = True
  # print "%d layers found\n" % len(nrList)
  # for i in range(0,len(nrList)) :
  #   print "  layer %d: tcmList = " % i, tcmList[i]
  #   print "            nrList = ", nrList[i]
  #   print "          tanDList = ", tanDeltaList[i]
  #   print " "
  fitConstraints = { "title" : title, \
                     "angIdeg" : angIdeg, \
                     "tcmRange" : totcmRange, \
                     "tcmList" : tcmList, \
                     "nrList" : nrList, \
                     "tanDeltaList" : tanDeltaList }
  return numpy.array(fGHz), numpy.array(trans), numpy.array(dphi), numpy.array(unc), fitConstraints


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

# returns min and max in save array
def tLimits( save, jlayer ) :
    vmin = 100.
    vmax = -100.
    if jlayer > len(save[0]) :
      print "error in findLimits"
    for n in range(0,len(save)) :
      if save[n][jlayer] < vmin :
        vmin = save[n][jlayer]
      if save[n][jlayer] > vmax :
        vmax = save[n][jlayer]
    return [vmin,vmax]
    
def nrefracFit2( infile ) :

    fullName = infile + ".dat"
    fGHz,trans,dphs_meas,unc,fitConstraints = readTransFile( fullName )
    tcmList,nrList = expandTree( fitConstraints["tcmRange"], fitConstraints["tcmList"], fitConstraints["nrList"] )
      #   tcmList = [ [tcm1,tcm2,...], [tcm1,tcm2,...], ... ] - expanded lists of thicknesses
      #   nrList = [ [nr1,nr2,...], [nr1,nr2,...], ... ] - expanded lists of refractive indices
      #   thus: tcmlist[i] and nrList[i] specify one set of thicknesses and refractive indices to be modeled

  # model all stacks in list, save avg and std phase and amp residuals
    nlayers = len( tcmList[0] )
    ntrials = len( tcmList )
    print "modeling %d layers, %d trials" % (nlayers, ntrials)

    tanDelta = []
    for n in range(0,nlayers) :
      tanDelta.append( fitConstraints["tanDeltaList"][n][0] )
    avgPhResid = 1.e6*numpy.ones( ntrials )
    stdPhResid = 1.e6*numpy.ones( ntrials )
    avgAmpResid = 1.e6*numpy.ones( ntrials )
    stdAmpResid = 1.e6*numpy.ones( ntrials )
    Resid = 1.e6*numpy.ones( ntrials )
    secs0 = time.time()
    for i in range(0,ntrials) :
      dphs_est,trans_est = solveStack( fGHz, fitConstraints["angIdeg"], tcmList[i], nrList[i], tanDelta ) 
      dphi = unwrap( dphs_meas-dphs_est )
        # unwrap keeps phase in range -180 to 180 in all cases
      avgPhResid[i] = numpy.average(dphi) 
      stdPhResid[i] = numpy.std(dphi) 
      avgAmpResid[i] = numpy.average( trans - trans_est )
      stdAmpResid[i] = numpy.std( trans - trans_est )

    # now attempt to make Resid ~ chisq/Npts; assume sigma = .02 for both real and imaginary differences
      measReal = trans * numpy.cos(dphs_meas*math.pi/180.)
      measImag = trans * numpy.sin(dphs_meas*math.pi/180.)
      modelReal = trans_est * numpy.cos(dphs_est*math.pi/180.)
      modelImag = trans_est * numpy.sin(dphs_est*math.pi/180.)
      Resid[i] = numpy.average( pow(measReal-modelReal,2)/0.0006 + pow(measImag-modelImag,2)/.0006 )

      if i > 0 and i%200 == 0 :
        secselapsed = time.time() - secs0
        secsremaining = (ntrials-i)*(secselapsed/float(i))
        if secsremaining/60. > 60. :
          print "... %d/%d finished; %.2f hours remaining" % (i,ntrials,secsremaining/3600.)
        else :
          print "... %d/%d finished; %.2f minutes remaining" % (i,ntrials,secsremaining/60.)

  # save the results for further analysis
    fout = open( infile+"fit.pickle", "wb" )     # write over previous file
    pickle.dump( [fitConstraints,tcmList,nrList,avgPhResid,stdPhResid,stdAmpResid,Resid], fout ) 
    fout.close() 

  # print out the best fits
    print "\n10 best fits minimizing avg residual amps:"
    printBest( avgAmpResid, tcmList, nrList ) 

    print "\n10 best fits minimizing std residual amps:"
    printBest( stdAmpResid, tcmList, nrList ) 

    print "\n10 best fits minimizing avg residual phases:"
    printBest( avgPhResid, tcmList, nrList ) 

    print "\n10 best fits minimizing std residual phases:"
    printBest( stdPhResid, tcmList, nrList ) 

    print "\n10 best fits minimizing overall residuals:"
    nbest = printBest( Resid, tcmList, nrList ) 

  # plot the fit
    time.sleep(1)
    plotFit( infile )

def minmax( varray, margin=.05 ) :
    '''find ymin and ymax for a plot'''
    vmin = varray.min()
    vmax = varray.max()
    diff = vmax - vmin
    return [ vmin - margin*diff, vmax + margin*diff ]

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
      

# converts alpha to tanDelta, using average index of 3.25
def tanD (fGHz, alpha) :
  print "tanDelta = %.2e" %  (alpha * clight/ (2.*math.pi*3.25*fGHz) )


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
        #print vec[0][i], vec[1][i], vec[2][i], vec[3][i], vec[4][i], vec[5][i], vec[6][i], vec[jz][i]
        x1 = vec[jx][i]
        y1 = vec[jy][i]
        z1 = vec[jz][i]
        x.append(x1)
        y.append(y1)
        z.append(z1)
    return numpy.array(x),numpy.array(y),numpy.array(z)

# compare amp and phase for different indices and thicknesses; can I measure the difference?
def checkSensitivity() :
    fGHzArr = numpy.arange(70.,115.01,.05)
    #fGHzArr = numpy.arange(210.,230.5,0.1)
    angIdeg = 5.2
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
    ax.set_ylim( [-9.,9.] )
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
    ax.set_ylim( [-9.,9.] )
    pyplot.grid(True)

    pyplot.savefig( pp, format="pdf" )
    pyplot.show()
    pp.close()

def checkSensitivity2() :
    # fGHzArr = numpy.append( numpy.arange(75.,115.5,.1), numpy.arange(210.,230.5,0.1) )
    fGHzArr = numpy.arange(70.,115.01,.04)
    angIdeg = 8.
    tcmArr1 = 2.54 * numpy.array( [.0093,.2461] )
    tcmArr2 = 2.54 * numpy.array( [.0093,.2451] )
    nrArr1 = numpy.array( [2.21, 3.113] )
    nrArr2 = numpy.array( [2.35, 3.115] )
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

def quickFind2( fGHz ):
    NoSolution = True
    for m in range(4,15) :
      f2 = fGHz/3.
      xGHz = (f2 - .04)/m
      if xGHz > 8.1 and xGHz < 12.4 :
        for n in range(4,15) :
          f1 = n*xGHz + .08
          if abs(2.*f1 - fGHz) < .05 :
            print fGHz, m, n, xGHz, f1, f2, 2*f1-fGHz
            NoSolution = False
    if NoSolution :
       print "no solution"
         

# show correlations between params for refractive index fits - specialized to just 2 parameters
def covarPlot2( infile, jopt=7, FTS=False ) :
    fin = open( infile, "r" )
    if FTS :
      [tcmList,nrList,Resid] = pickle.load( fin )
    # create array of results with 5 columns (t1,t2,nr1,nr2,Resid)
      vec = numpy.concatenate( ( [tcmList[:,0]/2.54], [tcmList[:,1]/2.54], [nrList[:,0]], [nrList[:,1]], \
        [Resid] ), axis=0 )
      jopt = 4
    else :
      [fitConstraints,tcmList,nrList,avgPhResid,stdPhResid,stdAmpResid,Resid] = pickle.load( fin )
      print fitConstraints
    # create array of results with 7 columns (t1,t2,nr1,nr2,avgPh,stdPh,stdAmp,Resid)
      vec = numpy.concatenate( ( [tcmList[:,0]/2.54], [tcmList[:,1]/2.54], [nrList[:,0]], [nrList[:,1]], \
        [avgPhResid], [stdPhResid], [stdAmpResid], [Resid] ), axis=0 )

    label = ["t1", "t2", "nr1", "nr2"]
  
  # find index of point with minimum avg or variance (whatever we are using)
    ibest = numpy.argmin( numpy.abs(vec[jopt]) )
    print len(vec[0]), ibest, vec[:,ibest]

  # create array of vectors for which vec[:,jopt] < 2.*vec[:,ibest]  
  # in other words, within 2 sigma contour
    print "all within 2 sigma"
    save = []
    for n in range(0,len(vec[0])) :
      if vec[jopt,n] < 2.*vec[jopt,ibest] :
        save.append( vec[:,n] )
    print save 
    print len(save)
    for j in range(0,4) :
      vmin = 100.
      vmax = -100.
      for n in range(0,len(save)) :
        if save[n][j] < vmin :
          vmin = save[n][j]
        elif save[n][j] > vmax :
          vmax = save[n][j]
      if j < 2 :
        print vmin,vmax
      else :
        print vmin, vmax, vmin*vmin, vmax*vmax

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
          zibest = numpy.nanmin(zi)
              # return position of interpolated minimum
              # zi[xibest,yibest] occurs at x[xibest,yibest],y[xibest,yibest]
          # print "xibest,yibest,zibest =",xibest,yibest,zi[xibest,yibest]
          # fig.imshow( zi, aspect='auto', origin='lower', vmin=vmin, vmax=vmax, \
          #    extent=[xmin,xmax,ymin,ymax], cmap='terrain' )
          clevels = [ 0.999*zibest,2.*zibest,3.*zibest ]
          colors = ["blue","green"]
          contour = pyplot.contourf(xi,yi,zi,clevels,colors=colors)
          #pyplot.colorbar(contour)
          fig.scatter(x,y,s=10, marker='x', c='black')
          fig.scatter(x[jbest],y[jbest],s=20, marker='x', c='white')
          fig.scatter(xi[xibest,yibest],yi[xibest,yibest],s=40, marker='s', c='red')
          fig.ticklabel_format(useOffset=False)
          fig.set_xlabel( label[jx], fontsize=12 )
          fig.set_ylabel( label[jy], fontsize=12 )
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
    if nplot > 1 :
      pyplot.savefig( pp, format="pdf" )
      pyplot.show()
    pp.close()

# generate an array of perfect fake data, for single layer or 2 layers
# can input this to nrefracFit2, use it to figure out how to weight amps and phases
def fakeData( outName, t1=.25, t2=0., n1=3.147, n2=1., sigma=0.0, tanDelta=6.e-4 ):
    fout = open( outName+".dat", "w" )
    fGHz1 = numpy.arange(74.,115.1,.04)
    fGHz2 = numpy.arange(220.,230.1,.1)
    #fGHzArr = numpy.concatenate( (fGHz1,fGHz2) )
    fGHzArr = fGHz1
    angIdeg = 5.16
    #tcmArr1 = 2.54 * numpy.array( [.008,.250] )
    #tcmArr1 = 2.54 * numpy.array( [.0093,.2461] )
    tcmArr1 = 2.54 * numpy.array( [t1] )
    nrArr1 = numpy.array( [n1] )
    tanDeltaArr = numpy.array( [tanDelta] ) 
    fout.write("# file generated by ob.fakeData() with sigma = %.4f\n" % sigma )
    fout.write("# angIdeg = %.3f\n" % angIdeg )
    fout.write("# -----------------------------\n")
    fout.write("# title = \"%s\"\n" % outName)
    #totin = (tcmArr1[0]+tcmArr1[1])/2.54
    totin = tcmArr1[0]/2.54
    fout.write("# totin = %.4f,%.4f\n" % (totin-.001,totin+.001))
    fout.write("# layer 1\n")
    fout.write("## true tin = %.5f\n" % (tcmArr1[0]/2.54) )
    fout.write("## true nr = %.3f\n" % nrArr1[0] )
    fout.write("#   tin = %.4f,%.4f,%.4f\n" % (tcmArr1[0]/2.54-.001, tcmArr1[0]/2.54+.0011,.0005))
    fout.write("#   nr = %.4f,%.4f,%.4f\n" % (nrArr1[0]-.5, nrArr1[0]+.505,.02))
    fout.write("#   tanDelta = %.2e\n" % tanDelta)
    #fout.write("# layer 2\n")
    #fout.write("## true tin = %.5f\n" % (tcmArr1[1]/2.54) )
    #fout.write("## true nr = %.3f\n" % nrArr1[1] )
    #fout.write("#   tin = %.4f,%.4f,%.4f\n" % (tcmArr1[1]/2.54-.001, tcmArr1[1]/2.54+.0011,.0005))
    #fout.write("#   nr = %.4f,%.4f,%.4f\n" % (nrArr1[1]-.002, nrArr1[1]+.0021,.002))
    #fout.write("#   tanDelta = 0.\n")
    #fout.write("# -----------------------------\n")
    phs1,trans1 = solveStack(fGHzArr, angIdeg, tcmArr1, nrArr1, tanDeltaArr ) 
    if sigma > .00001 :
      measReal = trans1 * numpy.cos(phs1*math.pi/180.) + numpy.random.normal( 0., sigma, len(trans1) )
      measImag = trans1 * numpy.sin(phs1*math.pi/180.) + numpy.random.normal( 0., sigma, len(phs1) )
      measAmp = numpy.absolute( (measReal + 1j * measImag ) )
      measPhs = numpy.angle( (measReal + 1j * measImag), deg=True ) 
    else :
      measAmp = trans1
      measPhs = phs1
    
    for n in range(0,len(fGHzArr)) :
      fout.write("%8.2f  %8.6f  %8.3f  0.005\n" % ( fGHzArr[n], measAmp[n], measPhs[n] ) )
    fout.close()

# produce Fourier transform of measurements; convert amp,phase to complex
# if FTS=False, 1st 3 columns of data file are interpreted as fGHz, trans_amp, trans_phase(deg)
# if FTS=True, 1st 3 columns of data file are interpreted as fGHz, trans_pwr, uncertainty 
# important note: fft returns zero freq term (sum of signal) as element 0, positive freq terms
#   as elements 1-N/2, neg freq terms as N/2+1:
#
def FFT( infile, FTS=False ) :
    fullName = infile + ".dat"
    fGHz,trans,phs,title,angIdeg,tcmRange,tcmListIn,nrListIn,tanDeltaList = readTransFile( fullName )
    if FTS :
      yf = numpy.fft.fft(trans) 
    else :
      measReal = trans * numpy.cos(phs*math.pi/180.)
      measImag = trans * numpy.sin(phs*math.pi/180.)
      yf = numpy.abs( numpy.fft.fft(measReal + 1j*measImag) )
      yf2 = numpy.abs( numpy.fft.fft( trans*trans) )
    delta = fGHz[1]-fGHz[0]
    tnsec = numpy.linspace( 0., 1./(2.*delta), len(yf)//2 )
    print "delta = ",delta
    print "tnsec = ",tnsec
    fig = pyplot.subplot(1,1,1)
    fig.plot( tnsec, yf[0:len(yf)//2], "-",color="r", linewidth=0.5 )
    if not FTS :
      fig.plot( tnsec, yf2[:len(yf)//2], "-",color="b", linewidth=0.5 )
    else :
      fig.plot( tnsec, numpy.abs(yf)[0:len(yf)//2], "-", color="b", linewidth=0.5)
    pyplot.show()
    

def FT2() :
    N=600
    T=1.0/800.
    x = numpy.linspace(0.,N*T,N)
    print x
    y = numpy.sin(50.*2.*numpy.pi*x)
    yf = scipy.fftpack.fft(y)
    xf = numpy.linspace(0.,1.0/(2.*T),N/2)
    fig = pyplot.subplot(1,1,1)
    fig.plot( xf, numpy.abs(yf[0:N//2]), color="r", linewidth=0.5 )
    pyplot.show()
    
# find rms of transmission amplitudes and phases in cols 2 and 3 of input file; use to
#   compute reproducibility of measurements using 2nd open slot, or of last-first reference measurement
def rms( infile ) :
    amp = []
    phs = []
    fin = open( infile, "r" )
    for line in fin:
      if not line.startswith("#") :
        a = line.split()
        amp.append( float(a[1]) )
        phs.append( float(a[2]) )
    fin.close()
    print numpy.std( numpy.array(amp) ), numpy.std( numpy.array(phs) )
    

# generate text string inidicating search range
def makeString( oneList, inches=False ) :
    if inches:
      if len(oneList) == 1 :
        return "%.4f\"" % oneList[0]
      elif len(oneList) == 2 :
        return "%.4f, %.4f\"" % (oneList[0],oneList[1])
      elif len(oneList) == 3 :
        return "%.4f\", %.4f, %.4f\"" % (oneList[0],oneList[1],oneList[2])
      else :
        return "%.4f, %.4f,... %.4f\"" % (oneList[0],oneList[1],oneList[-1])
    else :
      if len(oneList) == 1 :
        return "%.3f" % oneList[0]
      elif len(oneList) == 2 :
        return "%.3f, %.3f" % (oneList[0],oneList[1])
      elif len(oneList) == 3 :
        return "%.3f, %.3f, %.3f" % (oneList[0],oneList[1],oneList[2])
      else :
        return "%.3f, %.3f,... %.3f" % (oneList[0],oneList[1],oneList[-1])
1
# model Oliver's data
def nrefracFit3( infile ) :
    print " "
    fullName = infile + ".dat"
    fGHz,pwr,dummyphs,unc,fitConstraints = readTransFile( fullName )
    tcmList,nrList = expandTree( fitConstraints["tcmRange"], fitConstraints["tcmList"], fitConstraints["nrList"] )
    angIdeg = 0.
    print "\n%s" % fitConstraints["title"]
    print "total thickness in range ", fitConstraints["tcmRange"]

  # model all stacks in list, save avg and std phase and amp residuals
    nlayers = len( tcmList[0] )

  # caution: for now I am assuming that there is only 1 tanDelta per layer
    tanDelta = []
    for n in range(0,nlayers) :
      tanDelta.append( fitConstraints["tanDeltaList"][n][0] )

    ilen = len( tcmList )
    print "modeling %d layers, %d trials" % (nlayers, ilen)
    stdAmpResid = 1.e6*numpy.ones( ilen )

    secs0 = time.time()
    for i in range(0,ilen) :
      dphs_est,trans_est = solveStack( fGHz, angIdeg, tcmList[i], nrList[i], tanDelta ) 
         # assuming only one tanDelta per layer
      pwr_est = numpy.power(trans_est, 2)
      stdAmpResid[i] = numpy.std( (pwr-pwr_est)/unc )
      if i > 0 and i%200 == 0 :
        secselapsed = time.time() - secs0
        secsremaining = (ilen-i)*(secselapsed/float(i))
        if secsremaining/60. > 60. :
          print "... %d/%d finished; %.2f hours remaining" % (i,ilen,secsremaining/3600.)
        else :
          print "... %d/%d finished; %.2f minutes remaining" % (i,ilen,secsremaining/60.)

  # save the results for further analysis
    fout = open( infile+"fit.pickle", "wb" )     # write over previous file
    pickle.dump( [fitConstraints,tcmList,nrList,stdAmpResid], fout ) 
    fout.close() 

  # print out the best fits
    print "\n10 best fits minimizing std residual amps:"
    nbest = printBest( stdAmpResid, tcmList, nrList ) 

  # recalculate the model amps and phases for the best fit parameters
    dphs_est,trans_est = solveStack( fGHz, angIdeg, tcmList[nbest], nrList[nbest], tanDelta ) 
    pwr_est = numpy.power(trans_est, 2)

  # check numpy.std calculation of chisq
    var = 0.
    for i in range(0,len(pwr)) :
      var = var + pow( (pwr[i]-pwr_est[i])/unc[i], 2 )
    var = var/len(pwr)  
    print "my std = %.3f; numpy.std = %.3f" % (math.sqrt(var), stdAmpResid[nbest])
    plotFit( infile, FTS=True ) 
   

# dumps out file of differences vs freq for duplicates in a single dataset
def differ( inFile, outFile ):
    lastFreq = 0.
    lastAmp = 0.
    lastPhs = 0.
    varAmp = 0.
    varPhs = 0.
    npts = 0
    fin = open( inFile, "r" )
    fout = open( outFile, "w" )
    for line in fin:
      if not line.startswith("#") :
        a = line.split()
        if float(a[0]) == lastFreq :
          amp = float(a[1])
          phs = float(a[2])
          varAmp = varAmp + pow( amp-lastAmp, 2)
          varPhs = varPhs + pow( phs-lastPhs, 2)
          npts = npts+1
          fout.write("%8.2f  %8.6f %8.3f   %8.6f %8.3f   %8.6f %8.3f\n" % (lastFreq, lastAmp, lastPhs, amp, phs, amp-lastAmp, phs-lastPhs)  )
        else :
          lastFreq = float(a[0])
          lastAmp = float(a[1])
          lastPhs = float(a[2])
    print math.sqrt(varAmp/npts), math.sqrt(varPhs/npts)
    fin.close()
    fout.close()
  
# plot the fit on raw data, list parameters
# read from pickleFile
# if pdfOnly=True, do not plot to computer screen
#
def plotFit( infile, FTS=False, pdfOnly=False ) :

  # retrieve the raw data; ignore constraints given in data file because they might have changed since fit was done
    fullName = infile + ".dat"
    print "loading raw data from file %s" % fullName
    if FTS :
      fGHz,pwr,dummyphs,unc,dummyConstraints = readTransFile( fullName )
    else :
      fGHz,trans,dphs_meas,unc,dummyConstraints = readTransFile( fullName )
      print unc

  # retrieve the fit results
    fin = open( infile+"fit.pickle", "rb" )
    print "loading fit results from file %s" % (infile+"fit.pickle")
    if FTS :
      [fitConstraints,tcmList,nrList,Resid] = pickle.load(fin)
    else :
      [fitConstraints,tcmList,nrList,avgPhResid,stdPhResid,stdAmpResid,Resid] = pickle.load(fin)
    fin.close()
    nlayers = len( tcmList[0] )

    print "\n10 best fits minimizing overall residuals:"
    nbest = printBest( Resid, tcmList, nrList ) 

  # create savet and saven lists, containing all configurations with Resid < 2.* minResid
    savet = []
    savenr = []
    for n in range(0,len(Resid)) :
      if Resid[n] < 2.*Resid[nbest]:
        savet.append( tcmList[n] )
        savenr.append( nrList[n] )

  # begin the plot
    pyplot.ioff()
    pp = PdfPages( infile + ".pdf")
    fig =  pyplot.figure( figsize=(11,8) )
    pyplot.suptitle( "%s (%s)" % (fitConstraints["title"],fullName), fontsize=14 )

  # recalculate the model amps and phases for the best fit parameters; tanDelta is not currently optimized, so use input array
    tanDelta = []
    for n in range(0,nlayers) :
      tanDelta.append( fitConstraints["tanDeltaList"][n][0] )
    print tanDelta
    #dphs_est,trans_est = solveStack( fGHz, fitConstraints["angIdeg"], tcmList[nbest], nrList[nbest], fitConstraints["tanDeltaList"] ) 
    dphs_est,trans_est = solveStack( fGHz, fitConstraints["angIdeg"], tcmList[nbest], nrList[nbest], tanDelta ) 

    if FTS :
      pwr_est = numpy.power(trans_est, 2)

    # top panel is transmission vs freq, measured and theoretical
      ax = fig.add_axes( [.1,.6,.85,.32] )
      xmin,xmax = minmax( fGHz )
      ymin,ymax = minmax( numpy.concatenate((pwr,pwr_est)) )
      ax.axis( [xmin, xmax, ymin, ymax] )
      ax.errorbar( fGHz, pwr, yerr=unc, elinewidth=2, capsize=0. )
      ax.plot( fGHz, pwr_est, "r-" )
      ax.tick_params( labelsize=10 )
      ax.set_ylabel("pwr transmission", fontsize=10)
      pyplot.grid(True)

    # middle panel is measured-theoretical trans
      ax = fig.add_axes( [.1,.37,.85,.2])
      ax.tick_params( labelsize=10 )
      dpwr = pwr-pwr_est
      ymin,ymax = minmax( dpwr, margin=0.2 )
      ax.axis( [xmin, xmax, ymin, ymax] )
      ax.plot( fGHz, dpwr, "o-", ms=1. )
      ax.set_xlabel("freq (GHz)", fontsize=10)
      ax.set_ylabel("residual", fontsize=10)
      pyplot.grid(True)

    else :
      dphi = unwrap( dphs_meas-dphs_est )

    # top left panel is phase vs freq, measured and fit
      ax = fig.add_axes( [.1,.6,.38,.32] )
      xmin,xmax = minmax( fGHz )
      ymin,ymax = minmax( numpy.concatenate((dphs_meas,dphs_est)))
      ax.axis( [xmin, xmax, ymin, ymax] )
      ax.tick_params( labelsize=10 )
      ax.plot( fGHz, dphs_meas, "o", ms=2. )
      ax.plot( fGHz, dphs_est, "r-" )
      ax.set_ylabel("phase (deg)", fontsize=10)
      ax.set_ylim( [-195.,195.] )
      pyplot.grid(True)
  
    # middle left panel shows phase residual vs freq; avg shown as horiz dashed line
      ax = fig.add_axes( [.1,.37,.38,.2]) 
      ymin,ymax = minmax( dphi, margin=.2 )
      ax.axis( [xmin, xmax, ymin, ymax] )
      ax.plot( fGHz, dphi, "o-", ms=1. )
      ax.tick_params( labelsize=10 )
      dphiavg = numpy.average(dphi)
      ax.plot( [fGHz[0],fGHz[-1]],[dphiavg,dphiavg],"r--" )
      ax.set_ylabel("residual (deg)", fontsize=10)
      ax.set_xlabel("freq (GHz)", fontsize=10)
      ax.set_ylim( [-9.,9.] )
      pyplot.grid(True)

    # top right panel is transmission vs freq, measured and theoretical
      ax = fig.add_axes( [.57,.6,.38,.32] )
      ymin,ymax = minmax( numpy.concatenate((trans,trans_est)) )
      ax.axis( [xmin, xmax, ymin, ymax] )
      if numpy.mean(unc) > .5 :
        ax.plot( fGHz, trans, "o", ms=2. )
      else :
        ax.errorbar( fGHz, trans, yerr=unc, elinewidth=0.5, capsize=0 )
      ax.plot( fGHz, trans_est, "r-" )
      ax.tick_params( labelsize=10 )
      ax.set_ylabel("transmission", fontsize=10)
      ax.set_ylim( [.5,1.05] )
      pyplot.grid(True)

    # middle right panel is measured-theoretical trans
      ax = fig.add_axes( [.57,.37,.38,.2])
      ax.tick_params( labelsize=10 )
      dtrans = trans-trans_est
      ymin,ymax = minmax( dtrans, margin=0.2 )
      ax.axis( [xmin, xmax, ymin, ymax] )
      ax.plot( fGHz, dtrans, "o-", ms=1. )
      ax.set_xlabel("freq (GHz)", fontsize=10)
      ax.set_ylabel("residual", fontsize=10)
      ax.set_ylim( [-.09,.09] )
      pyplot.grid(True)

  # fictitious lower left panel gives room for fit ranges
    fontSize = 10
    ystep = .12
    ax = fig.add_axes( [.1,.1,.38,.2])
    ylab = 1 - ystep
    ax.text( 0.00, ylab, "fit constraints:", fontsize=fontSize)
    ylab = ylab - ystep
    ax.text( 0.03, ylab, "angle = %.2f deg" % fitConstraints["angIdeg"], transform=ax.transAxes, \
      horizontalalignment='left', fontsize=fontSize, color="black", rotation='horizontal' )
    ylab = ylab - ystep 
    ax.text( 0.03, ylab, "allowed thickness = [%s]" % makeString( fitConstraints["tcmRange"]/2.54, inches=True), transform=ax.transAxes, \
      horizontalalignment='left', fontsize=fontSize, color="black", rotation='horizontal' )
    tcmListIn = fitConstraints["tcmList"]
    nrListIn = fitConstraints["nrList"]
    ylab = ylab - ystep 
    ax.text( 0.00, ylab, "search ranges [t] [n] [tanDelta]:", fontsize=fontSize)
    for nl in range(0,nlayers) : 
      ylab = ylab - ystep 
      if nrListIn[nl][0] < 0. :
        ax.text( 0.03, ylab, "%d: mirror layer %d" % (nl+1, -1*nrListIn[nl][0]), \
        transform=ax.transAxes, horizontalalignment='left', fontsize=fontSize, \
        color="black", rotation='horizontal' )
      else :
        tstring = makeString( tcmListIn[nl]/2.54, inches=True )
        nstring = makeString( nrListIn[nl] )
        tanDstring = "%.4f" % tanDelta[nl]
        ax.text( 0.03, ylab, "%d: [%s]  [%s]  [%s]" % (nl+1,tstring,nstring,tanDstring), \
          transform=ax.transAxes, horizontalalignment='left', fontsize=fontSize, \
          color="black", rotation='horizontal' )
    pyplot.axis('off')

  # fictitious lower right panel gives room for text
    ax = fig.add_axes( [.57,.1,.38,.2])
    ylab = 1.-ystep
    ax.text( 0.00, ylab, "fit results:", fontsize=fontSize)
    ylab = ylab - ystep
    ax.text( 0.03, ylab, "reduced chisq = %.3f " % Resid[nbest], transform=ax.transAxes, \
      horizontalalignment='left', fontsize=fontSize, color="black", rotation='horizontal' )
    ylab = ylab - ystep 
    for nl in range(0,nlayers) :
      if nrListIn[nl][0] > 0. :
        tmin,tmax = tLimits( savet, nl)
        ax.text( 0.03, ylab, "layer %d:  t=%.4f\"   [%.4f, %.4f\"]" % \
          ( nl+1, tcmList[nbest][nl]/2.54, tmin/2.54, tmax/2.54) ,  \
          transform=ax.transAxes, horizontalalignment='left', fontsize=fontSize, \
          color="black", rotation='horizontal' )
        ylab = ylab - ystep
        nrmin,nrmax = tLimits( savenr, nl )
        ax.text( 0.03, ylab, "              n=%.4f  [%.4f, %.4f]" % \
          ( nrList[nbest][nl], nrmin, nrmax ), \
          transform=ax.transAxes, horizontalalignment='left', fontsize=fontSize, \
          color="black", rotation='horizontal' )
        ylab = ylab - ystep 
        ax.text( 0.03, ylab, "              e=%.3f  [%.3f, %.3f]" % \
          ( pow(nrList[nbest][nl],2), pow(nrmin,2), pow(nrmax,2) ),  \
          transform=ax.transAxes, horizontalalignment='left', fontsize=fontSize, \
          color="black", rotation='horizontal' )
        ylab = ylab - ystep
    pyplot.axis('off')
    #ax.text( 0.03, .42, "tanDelta = %.1e" % tanDelta[0], transform=ax.transAxes, \

    pyplot.savefig( pp, format="pdf" )
    if not pdfOnly :
      pyplot.show()
    pp.close()

# 9jan2018 - new version of nrefracFit which treats data as vectors, and which computes reduced
#   chisq from the measured uncertainties (these may be 1.0 if data were not averaged)

def nrefracFit4( infile, pdfOnly=False ) :
    fullName = infile + ".dat"
    fGHz,trans,dphs_meas,unc,fitConstraints = readTransFile( fullName )
    vecmeas = trans * numpy.exp(1j*numpy.radians(dphs_meas))   
    tcmList,nrList = expandTree( fitConstraints["tcmRange"], fitConstraints["tcmList"], fitConstraints["nrList"] )
      #   tcmList = [ [tcm1,tcm2,...], [tcm1,tcm2,...], ... ] - expanded lists of thicknesses
      #   nrList = [ [nr1,nr2,...], [nr1,nr2,...], ... ] - expanded lists of refractive indices
      #   thus: tcmlist[i] and nrList[i] specify one set of thicknesses and refractive indices to be modeled

  # model all stacks in list, save avg and std phase and amp residuals
    nlayers = len( tcmList[0] )
    ntrials = len( tcmList )
    print "modeling %d layers, %d trials" % (nlayers, ntrials)

  # for now, I am modeling only 1 value of tanDelta
    tanDelta = []
    for n in range(0,nlayers) :
      tanDelta.append( fitConstraints["tanDeltaList"][n][0] )

  # produce vecmodel for each set of parameters, calc reduced chisq
    chisq = 1.e6*numpy.ones( ntrials )
    secs0 = time.time()
    for i in range(0,ntrials) :
      phs,amp = solveStack( fGHz, fitConstraints["angIdeg"], tcmList[i], nrList[i], tanDelta ) 
      vecmodel = amp * numpy.exp(1j*numpy.radians(phs))
      chisq[i] = numpy.var( (vecmeas-vecmodel)/unc )
      print i, chisq[i]

      if i > 0 and i%200 == 0 :
        secselapsed = time.time() - secs0
        secsremaining = (ntrials-i)*(secselapsed/float(i))
        if secsremaining/60. > 60. :
          print "... %d/%d finished; %.2f hours remaining" % (i,ntrials,secsremaining/3600.)
        else :
          print "... %d/%d finished; %.2f minutes remaining" % (i,ntrials,secsremaining/60.)

  # save the results for further analysis
    fout = open( infile+"fit.pickle", "wb" )     # write over previous file
    pickle.dump( [fitConstraints,tcmList,nrList,chisq,chisq,chisq,chisq], fout ) 
    fout.close() 

  # plot the fit
    time.sleep(1)
    plotFit( infile, pdfOnly=pdfOnly )

# readFitFile is helper routine for nrefracFit5
# it reads file "sample.fitFile", where "sample" is the name of the sample
# this file should contain the following information: 
#   title = puck 1824
#   data = ../02jan2018/1822_D.avg.dat
#   data = ../10jan2018/1822_F.avg.dat
#   totin = ...
#   tin = ...
#   nr = ...
#   tanDelta = ...
#      
# it returns dictionary fitParams containing names of data files and parameter search ranges 
#
def readFitFile( sampleName ) :
  title = ""
  datafileList = []
  totcmRange = [0.,1000.]      
  tcmList = []
  nrList = []
  tanDeltaList = []
  fin = open( sampleName + ".fitFile", "r" )
  for line in fin :
    if line.find("#") > -1:
      line = line[0:line.find("#")]      # strip off comments, if any
    if len(line) > 0 :                   # ignore blank lines (or those beginning with #)
      if "title" in line :
        title = line[line.find("=")+1:].strip()
      elif "data" in line :
        datafileList.append( line[line.find("=")+1:].strip() )
      elif "totcm" in line :
        totcmRange = makeList( line[line.find("=")+1:] )
      elif "totin" in line :
        totcmRange = 2.54*makeList( line[line.find("=")+1:] )
      elif "tcm" in line :
        tcmList.append( makeList(line[line.find("=")+1:]) )
      elif "tin" in line :
        tcmList.append( makeList(line[line.find("=")+1:]) )
        if tcmList[-1][-1] > 0. :
          tcmList[-1] = 2.54*tcmList[-1]			# apply factor of 2.54 only for unmirrored layers!
      elif "nr" in line :
        nrList.append( makeList(line[line.find("=")+1:]) )
      elif "tanDelta" in line :
        tanDeltaList.append( makeList(line[line.find("=")+1:]) )
    # this is an alternative way of specifying mirrored layer
      elif "mirror" in line :
        tcmList.append( makeList(line[line.find("=")+1:]) )
        nrList.append( makeList(line[line.find("=")+1:]) )
        tanDeltaList.append( makeList(line[line.find("=")+1:]) )
  fin.close()
  fitParams = { "datafileList": datafileList, \
                "title" : title, \
                "tcmRange" : totcmRange, \
                "tcmList" : tcmList, \
                "nrList" : nrList, \
                "tanDeltaList" : tanDeltaList }
  return fitParams

# selects data points for pretest; these are points near the min and max transmission
# input: full data set
# output: sparse data set
# ideally, one would find period of the transmission curve using Fourier transform,
#   then fit to maxima and minima with polynomial fits to produce sparse data
# for now, I am identifying peaks by finding maxima, then masking out 0.8 of the
#   estimated period surrounding that peak; after all points are masked, I find
#   typical spacing between peaks, take midpoints to identify minima
#
def pretest( fGHz, trans, phase, unc ) :
    period = 7.5					    # estimated period of transmission data, GHz
    tma = numpy.ma.asarray( trans )     # create masked array of transmission data

  # create list of maxima, one per period
    nlist = []
    while ( tma.count() > 0) :	        # loop until all elements are masked
      nmax = numpy.ma.argmax(tma)       # maximum remaining element
      print "max at n = %d" % nmax
      nlist.append( nmax )              
      for i in range(0, len(tma) ) :
        if abs(fGHz[i] - fGHz[nmax]) < 0.75 * period :
          tma[i] = numpy.ma.masked
    nmaxlist = numpy.sort(nlist) 
    print "nmaxlist =",nmaxlist

  # assume midpoints between maxima are close to minima
    if nmaxlist[-1] == len(trans) - 1 :
      deltan = int( (nmaxlist[-2]-nmaxlist[-3])/2. + 0.5)
    else :
      deltan = int( (nmaxlist[-1]-nmaxlist[-2])/2. + 0.5)
    print "avg delta = %d" % deltan
    if (nmaxlist[0]-deltan) > 0 :
      nlist.append( nmaxlist[0] - deltan )
    for n in range(0,len(nmaxlist)-1) :
      nlist.append(nmaxlist[n] + deltan)
    if (nmaxlist[-1] + deltan) < len(tma) :
      nlist.append(nmaxlist[-1] + deltan)      
    
  # return sparse arrays
    nsorted = numpy.sort(nlist)
    print "nsorted =",nsorted
    fsparse = []
    tsparse = []
    psparse = []
    usparse = []
    for n in range(0,len(nsorted)) :
      fsparse.append( fGHz[nsorted[n]] )
      tsparse.append( trans[nsorted[n]] )
      psparse.append( phase[nsorted[n]] )
      usparse.append( unc[nsorted[n]] )
      print fsparse[-1],tsparse[-1],psparse[-1],usparse[-1]
    return numpy.array(fsparse), numpy.array(tsparse), numpy.array(psparse), numpy.array(usparse)      
    
def printElapsed( ntrials, i, secs0 ) :
    if i > 0 and i%200 == 0 :
      secselapsed = time.time() - secs0
      secsremaining = (ntrials-i)*(secselapsed/float(i))
      if secsremaining/60. > 60. :
        print "... %d/%d finished; %.2f hours remaining" % (i,ntrials,secsremaining/3600.)
      else :
        print "... %d/%d finished; %.2f minutes remaining" % (i,ntrials,secsremaining/60.)

# nrefracFit5 handles joint fits to multiple data sets, each with its own angIdeg and freq coverage
#   (but currently it will not handle both heterodyne and FTS data)
# it uses a pretest: calculates "pchisq" for a very sparse dataset that samples just max and min values;
#   then does the calculation for the full dataset only for those parameters where 
#   pchisq < pCutoff * pchisqBest
# fitFile lists input files, parameter ranges to search etc; but angIdeg comes from individual files
#
def nrefracFit5( fitFile, pCutoff=10., pdfOnly=False ) :
    fitParams = readFitFile( fitFile )
    print fitParams

  # read the data file(s); store as LISTS of arrays; these are not all appended into one array
  #   in case angIdeg differs for them
    fGHz = []
    vecmeas = []
    unc = []
    angIdeg = []
    for i, datafile in enumerate( fitParams["datafileList"] ) :
      f,trans,dphs_meas,u,fitConstraints = readTransFile( datafile )     # read single data file

    # for pretest, create sparse data set for each full dataset, probing peaks and valleys
      pf, ptrans, pphs, punc =  pretest( f, trans, dphs_meas, u ) 
      fGHz.append( pf )
      vecmeas.append( ptrans * numpy.exp(1j*numpy.radians(pphs)) )
      unc.append( punc )
      angIdeg.append( fitConstraints["angIdeg"] )  # this is the only fitConstraint from the data file that we use

    # now append full data set to lists
      fGHz.append(f)
      vecmeas.append( trans * numpy.exp(1j*numpy.radians(dphs_meas)) )  
      unc.append(u)
      angIdeg.append( fitConstraints["angIdeg"] )  # this is the only fitConstraint from the data file that we use

    print "read in %d data files" % len(fitParams["datafileList"])
    print "created %d datasets" % len(fGHz)

  # now that the data are in hand, create tree of parameters to search
    tcmList,nrList = expandTree( fitParams["tcmRange"], fitParams["tcmList"], fitParams["nrList"] )
      #   tcmList = [ [tcm1,tcm2,...], [tcm1,tcm2,...], ... ] - expanded lists of thicknesses
      #   nrList = [ [nr1,nr2,...], [nr1,nr2,...], ... ] - expanded lists of refractive indices
      #   thus: tcmlist[i] and nrList[i] specify one set of thicknesses and refractive indices to be modeled
    nlayers = len( tcmList[0] )
    ntrials = len( tcmList )
    print "modeling %d layers, %d trials" % (nlayers, ntrials)

  # for now, I am allowing only 1 value of tanDelta per layer
    tanDelta = []
    for n in range(0,nlayers) :
      tanDelta.append( fitParams["tanDeltaList"][n][0] )

  # for Pretest, model odd-numbered layers only, calculate reduced chisq for every possible trial
    pchisqCutoff = 1.e6
    secs0 = time.time()       
    pchisq = 1.e6*numpy.ones( ntrials )
    print "\nbeginning Pretest"
    for i in range(0,ntrials) :
      errs = []
      for j in range(0,len(fGHz),2) :
        phs,amp = solveStack( fGHz[j], angIdeg[j], tcmList[i], nrList[i], tanDelta ) 
        vecmodel = amp * numpy.exp(1j*numpy.radians(phs))
        errs.append( (vecmeas[j]-vecmodel)/unc[j] )
      pchisq[i] = numpy.var( numpy.concatenate(errs) )
      printElapsed( ntrials, i, secs0 )
    nbest = printBest( pchisq, tcmList, nrList ) 
    pchisqCutoff = pCutoff*pchisq[nbest]

  # figure out how many trials will be needed for full test
    ntrials2 = bisect.bisect_left( numpy.sort(pchisq), pchisqCutoff )
    print "full data set will be modeled for %d out of %d trials" % (ntrials2,ntrials)
        
  # now model full dataset, only for trials where pchisq < pchisqCutoff
    secs0 = time.time()       
    chisq = 1.e6*numpy.ones( ntrials )
    imodel = 0
    print "\nbeginning full test"
    for i in range(0,ntrials) :
      printElapsed( ntrials2, imodel, secs0 )
      if pchisq[i] > pchisqCutoff :
        continue
      else :
        imodel = imodel+1
        errs = []
        for j in range(1, len(fGHz), 2) :
          phs,amp = solveStack( fGHz[j], angIdeg[j], tcmList[i], nrList[i], tanDelta ) 
          vecmodel = amp * numpy.exp(1j*numpy.radians(phs))
          errs.append( (vecmeas[j]-vecmodel)/unc[j] )
      chisq[i] = numpy.var( numpy.concatenate(errs) )
    nbest = printBest( chisq, tcmList, nrList ) 
          
  # save results for further analysis in a unique pickle file (don't overwrite old file)
    n = 1
    while os.path.isfile( fitFile + ".pickle;%d" % n) :
      n = n+1
    fout = open( fitFile + ".pickle;%d" % n, "wb" )     # create new file
    pickle.dump( [fitParams,tcmList,nrList,pchisq,chisq], fout ) 
    fout.close() 

  # plot the fit
  #  time.sleep(1)
  #  plotFit5( infileList, pdfOnly=pdfOnly )
