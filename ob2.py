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
from mpl_toolkits.mplot3d import Axes3D
import scipy
from scipy.interpolate import splrep, splev, interp1d
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

# readTransFile inputs data
# default - assume created by unpack3.py: col1 = fGHz, col2 = amp, col3 = phase, col4 = vector uncertainty
# following line "# FTS", interpret col2 = pwr, col3 = pwr unc
# default - assume angle of incidence = 0; otherwise expect "# angIdeg = xx.xx"
#
def readTransFile( infile, path="/o/plambeck/PolarBear/OpticsBench/" ) :
  fin = open( path + infile, "r" )
  fGHz = []
  trans = []
  dphi = []
  unc = []
  angIdeg = -1000.
  FTS = False
  for line in fin :
    if not line.startswith( "#" ) :    # lines without hash marks are data
      a = line.split()
      fGHz.append( float( a[0] ) )
      trans.append( float( a[1] ) )
      if FTS :
        dphi.append( 0. )
        unc.append( float( a[2] ) )
      else :
        dphi.append( float( a[2] ) )
        unc.append( float(a[3]) )        
    else :                               # lines with single hash marks are search parameters
      if "angIdeg" in line :
        angIdeg = float( line[line.find("=")+1:] )
      elif "FTS" in line :
        FTS = True
  if angIdeg == -1000. :
    sys.exit( "FATAL ERROR: angIdeg is not given in data file %s" % (path+infile) )
  else :
    return numpy.array(fGHz), numpy.array(trans), numpy.array(dphi), numpy.array(unc), angIdeg, FTS


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
     
# expandTree2 uses meshgrid to expand the thickness and refractive index arrays
#
# input: 
#   tcmRange = [tcmMin,tcmMax]  - specifies allowed range of total thicknesses
#   tcmListIn = [ [tcm1list], [tcm2list], ...] - lists of thicknesses to model for each layer
#   nrListIn = [ [nr1list], [nr2list], ... ] - lists of refractive indices to model for each layer
#   tanDListIN - list of refractive indices to model for each layer
# output :
#   tcmListOut = [ [tcm1,tcm2,...], [tcm1,tcm2,...], ... ] - expanded lists of thicknesses
#   nrListOut = [ [nr1,nr2,...], [nr1,nr2,...], ... ] - expanded lists of refractive indices
#
# thus: tcmlistOut[i], nrListOut[i], tanDListOut[i] specify one set of thicknesses, refractive indices,
#   and loss tangents that should be modeled; 
#    ote that same list of thicknesses may appear many times, one
#   for each possible set of refractive indices; same for refractive indices
# total thickness of tcmListOut is guaranteed to be in the range [tcmMin,tcmMax]
# 
# the deficiency of this approach is that maybe only a small set of thicknesses is possible; it's
#   dumb to create full very large expanded tree, then prune out all layers with unacceptable thicknesses
#
def expandTree2( tcmRange, tcmListIn, nrListIn, tanDListIn ) :
    # print "tcmListIn = ",tcmListIn
    # print "nrListIn = ",tcmListIn
    # print "tanDListIn = ", tanDListIn
    nlayers = len(tcmListIn)
    # print "nlayers = ", nlayers
    # print "tcmRange = ", tcmRange

  # before calling meshgrid, calculate size of arrays that it will produce
  # also, create 1st list of parameters to model, just in case there is is for modeling
    ntrials = 1
    tcmList = []
    nrList = []
    tanDList = []
    for n in range(0,nlayers) :
      ntrials = ntrials * len(tcmListIn[n]) * len(nrListIn[n]) * len(tanDListIn[n])
      tcmList.append( tcmListIn[n][0] )
      nrList.append( nrListIn[n][0] )
      tanDList.append( tanDListIn[n][0] )
   
  # in the special case where ntrials=1 (as in modeling a stack), we are finished!
    if ntrials == 1 :
      print "modeling only 1 trial"
      return numpy.array([tcmList]),numpy.array([nrList]),numpy.array([tanDList])
    elif ntrials > 1.e12 :
      print "exiting! too many possible trials, estimated at %.2e!" % ntrials
      return

  # meshgrid arrays of thicknesses, refractive indices, and tanDs
  # this is lame, but I can't figure out how to do this elegantly, without if structure
    if nlayers == 1 :
      a = numpy.meshgrid( tcmListIn[0], nrListIn[0], tanDListIn[0] )
    elif nlayers == 2 :
      a = numpy.meshgrid( tcmListIn[0],tcmListIn[1],nrListIn[0],nrListIn[1],tanDListIn[0],tanDListIn[1] ) 
    elif nlayers == 3 :
      a = numpy.meshgrid( tcmListIn[0],tcmListIn[1],tcmListIn[2], \
            nrListIn[0],nrListIn[1],nrListIn[2], \
            tanDListIn[0],tanDListIn[1],tanDListIn[2] )
    elif nlayers == 4 :
      a = numpy.meshgrid( tcmListIn[0],tcmListIn[1],tcmListIn[2],tcmListIn[3],\
            nrListIn[0],nrListIn[1],nrListIn[2],nrListIn[3], \
            tanDListIn[0],tanDListIn[1],tanDListIn[2],tanDListIn[3] )
    elif nlayers == 5 :
      a = numpy.meshgrid( tcmListIn[0],tcmListIn[1],tcmListIn[2],tcmListIn[3],tcmListIn[4], \
            nrListIn[0],nrListIn[1],nrListIn[2],nrListIn[3],nrListIn[4], \
            tanDListIn[0],tanDListIn[1],tanDListIn[2],tanDListIn[3],tanDListIn[4] )
    elif nlayers == 6 :
      a = numpy.meshgrid( tcmListIn[0],tcmListIn[1],tcmListIn[2],tcmListIn[3],tcmListIn[4],tcmListIn[5], \
            nrListIn[0],nrListIn[1],nrListIn[2],nrListIn[3],nrListIn[4],nrListIn[5], \
            tanDListIn[0],tanDListIn[1],tanDListIn[2],tanDListIn[3],tanDListIn[4],tanDListIn[5] )
    elif nlayers == 7 :
      a = numpy.meshgrid( tcmListIn[0],tcmListIn[1],tcmListIn[2],tcmListIn[3],tcmListIn[4],tcmListIn[5],tcmListIn[6], \
            nrListIn[0],nrListIn[1],nrListIn[2],nrListIn[3],nrListIn[4],nrListIn[5],nrListIn[6], \
            tanDListIn[0],tanDListIn[1],tanDListIn[2],tanDListIn[3],tanDListIn[4],tanDListIn[5],tanDListIn[6] )
    else :
      print "nlayers = %d not allowed by expandTree" % nlayers
      return

  # flatten each array contained in a; this should create an array with length 3*nlayers 
    b = []
    for i in range(0,len(a)) :
      b.append( a[i].flatten() )

  # constrained layers are indicated by -L, where L is the corresponding layer
    for n in range(0,nlayers) :
      if b[n][0] < 0. :
        m = int( abs(b[n][0]) + 0.5 )
        try :
          b[n] = b[m-1]
          # print "mirrored thicknesses for layer %d into layer row %d" % (m,n+1)
        except :
          print "mirroring failed; tried to copy row %d into row %d" % (m-1,n)
          return
      if b[n+nlayers][0] < 0. :
        m = int( abs(b[n+nlayers][0]) + 0.5 )    # e.g., -1 -> 1
        try :
          b[n+nlayers] = b[m+nlayers-1]
          # print "mirrored refractive indices for layer %d into layer %d" % (m,n+1)
        except :
          print "mirroring failed; tried to copy row %d into row %d" % (m+nlayers-1,n+nlayers)
          return
      if b[n+2*nlayers][0] < 0. :
        m = int( abs(b[n+2*nlayers][0]) + 0.5 )    # e.g., -1 -> 1
        try :
          b[n+2*nlayers] = b[m+2*nlayers-1]
          # print "mirrored tanDeltas for layer %d into layer %d" % (m,n+1)
        except :
          print "mirroring failed; tried to copy row %d into row %d" % (m+2*nlayers-1,n+2*nlayers)
          return

  # compute total thickness for each possible stack
    totcm = numpy.zeros( len(b[0]) )
    for i in range(0,nlayers) :
      totcm = totcm + b[i]
  
  # now form transpose, retain only rows where thickness is within allowed range
    c = numpy.transpose(b)
    d = []
    for n in range(0,len(c)) :
      if (totcm[n] >= tcmRange[0]) and (totcm[n] <= tcmRange[1]) :
        d.append( c[n] )

  # send error message and quit if no rows meet the thickness criterion
    if len(d) == 0 :
      sys.exit("\n** exiting: none of the proposed stacks meet thickness criterion **\n\n")

  # return separate tcm, nr, and tanD arrays
    # print "final ntrials = %d" % len(numpy.hsplit(numpy.array(d),3)[0])
    return numpy.hsplit( numpy.array(d),3)


# prints parameters for ndump stacks with smallest absolute deviations; returns index of best fit
def printBest( merit, tcmList, nrList, tanDeltaList, ndump=10 ) :
    ntrials = len( tcmList )
    nlayers = len( tcmList[0] )
    if ndump > ntrials :
      ndump = ntrials
    merit2 = numpy.abs(merit)
    nbest = merit2.argsort()[:ndump]    # finds indexes of ndump smallest values
    print nbest
    for nb in nbest:
      print tcmList[nb]/2.54, nrList[nb], tanDeltaList[nb], merit[nb]
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
    
# return best fit and uncertainties for each parameter
# returns tcmFit = [ [layer1_min,layer1_best,layer1_max], [layer2_min,layer2_best,layer2_max] ...], etc
def summarize( merit, tcmList, nrList, tanDeltaList, cutoff=2. ) :
    nlayers = len(tcmList[0])

  # fill out fit arrays with -1., just to spot trouble in case it occurs
    tcmFit = -1.*numpy.ones( (nlayers,3) )
    nrFit = -1.*numpy.ones( (nlayers,3) )
    tanDFit = -1.*numpy.ones( (nlayers,3) )

  # sort indexes in order of merit; best fit parameters have index nsorted[0]
    nsorted = merit.argsort()   
    tcmsave = []
    nrsave = []
    tanDsave = []

  # save parameter arrays for which  merit <= cutoff * merit_min 
    tcmsave = []
    nrsave = []
    tanDsave = []
    for n in range(0,len(nsorted)) :
      if merit[n] <= cutoff*merit[nsorted[0]] :
        tcmsave.append( tcmList[n] )
        nrsave.append( nrList[n] )
        tanDsave.append( tanDeltaList[n] )
    
  # for each layer, find min and max for each parameter
  # must convert lists to numpy arrays order to index them as [:,nl]
    for nl in range(0, nlayers ) :
      tcmFit[nl,0] = numpy.array(tcmsave)[:,nl].min()  
      tcmFit[nl,1] = tcmList[nsorted[0],nl]
      tcmFit[nl,2] = numpy.array(tcmsave)[:,nl].max()  
      nrFit[nl,0] = numpy.array(nrsave)[:,nl].min()  
      nrFit[nl,1] = nrList[nsorted[0],nl]
      nrFit[nl,2] = numpy.array(nrsave)[:,nl].max()  
      tanDFit[nl,0] = numpy.array(tanDsave)[:,nl].min()  
      tanDFit[nl,1] = tanDeltaList[nsorted[0],nl]
      tanDFit[nl,2] = numpy.array(tanDsave)[:,nl].max()  

    #for nb in nbest:
    #  print tcmList[nb]/2.54, nrList[nb], tanDeltaList[nb], merit[nb]
    return merit[nsorted[0]], tcmFit, nrFit, tanDFit

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
        # print vec[0][i], vec[1][i], vec[2][i], vec[3][i], vec[4][i]
        x.append( vec[jx][i] )
        y.append( vec[jy][i] )
        z.append( vec[jz][i] )
    return numpy.array(x),numpy.array(y),numpy.array(z)

# compare predicted transmission for 2 pucks
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

  # top left panel compares phase vs freq
    ax = fig.add_axes( [.1,.6,.38,.32] )
    xmin,xmax = minmax( fGHzArr )
    ymin,ymax = minmax( numpy.concatenate((phs1,phs2)))
    ax.axis( [xmin, xmax, ymin, ymax] )
    ax.tick_params( labelsize=10 )
    ax.plot( fGHzArr, phs1, "b-" )
    ax.plot( fGHzArr, phs2, "r-" )

    ax.set_ylabel("phase (deg)", fontsize=10)
    pyplot.grid(True)

  # middle left panel shows phase difference
    ax = fig.add_axes( [.1,.37,.38,.2])
    ymin,ymax = minmax( dphi, margin=.2 )
    #if abs(ymin) < abs(ymax) :
    #  ymin = math.copysign(ymax,ymin)
    #elif abs(ymax) < abs(ymin) :
    #  ymax = math.copysign(ymin,ymax)
    ax.axis( [xmin, xmax, ymin, ymax] )
    ax.plot( fGHzArr, dphi, "b-" )
    pyplot.grid(True)
                                                       
  # top right panel compares transmission vs freq
    ax = fig.add_axes( [.57,.6,.38,.32] )
    ymin,ymax = minmax( numpy.concatenate((trans1,trans2)) )
    ax.axis( [xmin, xmax, ymin, ymax] )
    ax.plot( fGHzArr, trans1, "-", color='blue' )
    ax.plot( fGHzArr, trans2, "r-" )
    ax.tick_params( labelsize=10 )
    ax.set_ylabel("transmission", fontsize=10)
    pyplot.grid(True)

  # middle right panel is transmission difference
    ax = fig.add_axes( [.57,.37,.38,.2])
    ax.tick_params( labelsize=10 )
    ymin,ymax = minmax( trans1-trans2, margin=0.2 )
    ax.axis( [xmin, xmax, ymin, ymax] )
    ax.plot( fGHzArr, trans1-trans2, "-", color='blue' )
    ax.set_xlabel("freq (GHz)", fontsize=10)
    ax.set_ylabel("residual", fontsize=10)
    pyplot.grid(True)

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
         

# show correlations between params for refractive index fits - specialized to just 2 layers
def covarPlot2( pickleFile, chisqIndex=-1, worstRatio=10., FTS=False ) :
    fin = open( pickleFile, "r" )

  # create array of results with 5 columns (t1,t2,nr1,nr2,chisq)
    if FTS :
      [tcmList,nrList,Resid] = pickle.load( fin )
      vec = numpy.concatenate( ( [tcmList[:,0]/2.54], [tcmList[:,1]/2.54], [nrList[:,0]], [nrList[:,1]], \
        [Resid] ), axis=0 )
    else :
      [fitParams,tcmList,nrList,tanDeltaList,chisq] = pickle.load( fin )
      vec = numpy.concatenate( ( [tcmList[:,0]/2.54], [tcmList[:,1]/2.54], [nrList[:,0]], [nrList[:,1]], \
        [chisq[chisqIndex]]), axis=0 )

  # this routine is specialized to just 2 layers
    nlayers = len( tcmList[0] )
    if nlayers != 2 :
      print "nlayers = %d; this routine is specialized to exactly 2 layers; exiting" % nlayers
      return
  
    label = ["t1", "t2", "nr1", "nr2"]
  
  # find index of point with minimum chisq
    ibest = printBest( chisq[chisqIndex], tcmList, nrList, tanDeltaList ) 
    chisqMin = chisq[chisqIndex][ibest]

    pyplot.ioff()
    pp = PdfPages("Covar.pdf")
    pyplot.figure( figsize=(11,8) )
    nplot = 0
    for jx in range(0,4) :
      for jy in range( jx+1,4 ) :
        x,y,z = stripout( vec, ibest, jx, jy, 4 )
        print "\n"
        jbest = numpy.argmin(z)
        zbest = numpy.nanmin(z)
        print "zbest = ",zbest

        if (len(numpy.unique(x)) > 1) and (len(numpy.unique(y)) > 1) :
          nplot = nplot + 1
          print "nplot = ",nplot
          if nplot > 6 :
            pyplot.savefig( pp, format="pdf" )
            pyplot.show()
            pyplot.figure( figsize=(11,8) )
            nplot = 1
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
          print "zibest = ",zibest
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
          fig.scatter(x[jbest],y[jbest],s=20, marker='x', c='red')
          fig.scatter(xi[xibest,yibest],yi[xibest,yibest],s=40, marker='s', c='red')
          fig.ticklabel_format(useOffset=False)
          fig.set_xlabel( label[jx], fontsize=12 )
          fig.set_ylabel( label[jy], fontsize=12 )
          fig.scatter(x,y,c=z,s=200,edgecolors='none',marker='s',vmin=0.,vmax=.3*z.max())
          #fig.scatter(x,y,c=z,s=200,edgecolors='none',marker='s',vmin=50.,vmax=200.)
          fig.axis( [xmin, xmax, ymin, ymax],fontsize=6 )
          pyplot.locator_params( axis='x', nbins=5)
          pyplot.locator_params( axis='y', nbins=5)
          fig.tick_params( axis='both', which='major', labelsize=8 )
          pyplot.tight_layout( pad=3, h_pad=1)

    if nplot > 1 :
      pyplot.savefig( pp, format="pdf" )
      pyplot.show()
    pp.close()

# plot chisq for each combination of variables, "marginalized" over all the others (I think)
# that is, for each point on the var1,var2 plane, plot the minimum chisq, for any combination
#    of the other variables
# this differs from my original implementation, where I plotted the chisq for the values
#    of the other variables that matched the best fit
#
# bug: tries to plot mirrored parameters against each other
#
def covarPlot3( pickleFile, chisqIndex=-1, worstRatio=10. ) :

    fin = open( pickleFile, "r" )
    [fitParams,tcmList,nrList,tanDeltaList,chisq] = pickle.load( fin )
    nlayers = len(tcmList[0])
    label = []
    for nl in range(0,nlayers) :
      label.append("t%d" % (nl+1))
    for nl in range(0,nlayers) :
      label.append("nr%d" % (nl+1))
    for nl in range(0,nlayers) :
      label.append("tanD%d" % (nl+1))
    var = numpy.concatenate( (tcmList/2.54,nrList,tanDeltaList), axis=1 )
      # this makes a flat list of variables for each trial
    print "labels: ", label 
    print "shape(tcmList) = ", numpy.shape(tcmList)
    print "shape(nrList) = ", numpy.shape(nrList)
    print "shape(tanDeltaList) = ", numpy.shape(tanDeltaList)
    print "shape(chisq) = ", numpy.shape(chisq)

  # sort by chisq
    isort = numpy.argsort( chisq[chisqIndex] )
    print isort[0:10]
    chisqMin = chisq[chisqIndex][isort[0]]
    print "chisqMin = %.3f" % chisqMin
    nbest = printBest( chisq[chisqIndex], tcmList, nrList, tanDeltaList ) 

  # now plot one panel for each pair of non-trivial parameters
    pyplot.ioff()
    pp = PdfPages( pickleFile+".covar.pdf")
    pyplot.figure( figsize=(11,8) )
    nplot = 0
    for jx in range(0,len(label)) :
      for jy in range( jx+1,len(label)) :
        if (len(numpy.unique(var[:,jx])) > 1) and (len(numpy.unique(var[:,jy])) > 1) :     # both parameters have >1 value
          npts = 0
          xl = []
          yl = []
          zl = []
          for i in isort :
            npts = npts + 1
            if chisq[chisqIndex][i] > worstRatio * chisqMin :
              # print "%d / %d points meet chisq criterion" % (npts, len(isort) )
              break
            else :
              newPoint = True
              for j in range(0,len(xl)) :
                if ( var[i][jx] == xl[j] ) and ( var[i][jy] == yl[j] ) :
                  newPoint = False   # lower chisq already plotted for this x,y
              if newPoint :
                xl.append( var[i][jx] ) 
                yl.append( var[i][jy] ) 
                zl.append( chisq[chisqIndex][i])    
                # print i, xl[-1], yl[-1], zl[-1]
          x = numpy.array(xl)
          y = numpy.array(yl)
          z = numpy.array(zl)
    
          nplot = nplot + 1
          if nplot > 6 :
            pyplot.suptitle( pickleFile, fontsize=14 )
            pyplot.savefig( pp, format="pdf" )
            pyplot.show()
            pyplot.figure( figsize=(11,8) )
            nplot = 1
          fig = pyplot.subplot(2,3,nplot)
          xmin,xmax = minmax(x,margin=0.01)
          ymin,ymax = minmax(y,margin=0.01)

        # section below copied from web example
          xi,yi = numpy.linspace( xmin, xmax, 50), numpy.linspace(ymin, ymax, 50)
          xi,yi = numpy.meshgrid(xi,yi)
          try :
            zi = scipy.interpolate.griddata( (x,y) ,z, (xi,yi), method='linear')
          except :
            zi = scipy.interpolate.griddata( (x,y) ,z, (xi,yi), method='nearest')
          xibest,yibest = numpy.unravel_index(numpy.nanargmin(zi),zi.shape)
          zibest = numpy.nanmin(zi)
          print label[jx], label[jy], " plot %d points  zibest = %.2f" % ( len(xl), zibest )
              # return position of interpolated minimum
              # zi[xibest,yibest] occurs at x[xibest,yibest],y[xibest,yibest]
          # print "xibest,yibest,zibest =",xibest,yibest,zi[xibest,yibest]
          # fig.imshow( zi, aspect='auto', origin='lower', vmin=vmin, vmax=vmax, \
          #    extent=[xmin,xmax,ymin,ymax], cmap='terrain' )
          clevels = [ 0.999*z[0],2.*z[0],3.*z[0] ]
          colors = ["blue","green"]
          contour = pyplot.contourf(xi,yi,zi,clevels,colors=colors)
          #pyplot.colorbar(contour)
          fig.scatter(x,y,s=10, marker='x', c='black')
          fig.scatter(x[0],y[0],s=40, marker='s', c='red')
          fig.scatter(xi[xibest,yibest],yi[xibest,yibest],s=40, marker='x', c='red')
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

    if nplot > 1 :
      pyplot.suptitle( pickleFile, fontsize=14 )
      pyplot.savefig( pp, format="pdf" )
      pyplot.show()
    pp.close()

# produce Fourier transform of measurements; convert amp,phase to complex
# if FTS=False, 1st 3 columns of data file are interpreted as fGHz, trans_amp, trans_phase(deg)
# if FTS=True, 1st 3 columns of data file are interpreted as fGHz, trans_pwr, uncertainty 
# important note: fft returns zero freq term (sum of signal) as element 0, positive freq terms
#   as elements 1-N/2, neg freq terms as N/2+1:
#
def FFT( infile, FTS=False ) :
    fGHz,trans,phs,u,angI,FTS = readTransFile( infile )     # read single data file
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
    #fig.plot( tnsec, yf[0:len(yf)//2], "-",color="r", linewidth=0.5 )
    if not FTS :
      fig.plot( tnsec[1:], yf2[1:len(yf)//2], "-",color="b", linewidth=0.5 )
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
def makeString( oneList, fmt=4 ) :
    if fmt==4 :
      if len(oneList) == 1 :
        return "%.4f" % oneList[0]
      elif len(oneList) == 2 :
        return "%.4f, %.4f" % (oneList[0],oneList[1])
      elif len(oneList) == 3 :
        return "%.4f\", %.4f, %.4f" % (oneList[0],oneList[1],oneList[2])
      else :
        return "%.4f, %.4f,... %.4f" % (oneList[0],oneList[1],oneList[-1])
    else :
      if len(oneList) == 1 :
        return "%.3f" % oneList[0]
      elif len(oneList) == 2 :
        return "%.3f, %.3f" % (oneList[0],oneList[1])
      elif len(oneList) == 3 :
        return "%.3f, %.3f, %.3f" % (oneList[0],oneList[1],oneList[2])
      else :
        return "%.3f, %.3f,... %.3f" % (oneList[0],oneList[1],oneList[-1])

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
  
# readFitFile is helper routine for nrefracFit, fuzzy1, fuzzy2
# for nrefracFit, it specifies the search ranges to fit
# for fuzzy1 and fuzzy2 it specifies the parameters to model
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
  deltaTinList = []
  fin = open( sampleName + ".fitFile", "r" )
  for line in fin :
    if line.find("#") > -1:
      line = line[0:line.find("#")]      # strip off comments, if any
    if len(line) > 0 :                   # ignore blank lines (or those beginning with #)
      if "title" in line :
        title = line[line.find("=")+1:].strip("\n")
        print "title = %s" % title
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
      elif "deltaTin" in line :
        deltaTinList.append( float(line[line.find("=")+1:] ) )    # allow only 1 value!
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
                "tanDeltaList" : tanDeltaList,
                "deltaTinList" : deltaTinList }
  return fitParams

def findPeriod( fGHz, trans ) :
    #fGHz,trans,phs,u,angI,FTS = readTransFile( infile )     # read single data file
    yf = numpy.abs(numpy.fft.fft(trans))
    delta = fGHz[1]-fGHz[0]
    tnsec = numpy.linspace( 0., 1./(2.*delta), len(yf)//2 )
    index = numpy.argmax(yf[1:len(yf)/2]) + 1
    print yf[:5]
    print tnsec[:5]
    period = 1./tnsec[index]
    print "period = %.2f GHz" % period
    return period

# selects data points for prefit; these are points near the min and max transmission
# input: full data set
# output: sparse data set
# ideally, one would find period of the transmission curve using Fourier transform,
#   then perhaps use polynomial fits to produce sparse data
# for now, I am identifying peaks by finding maxima, then masking out 0.8 of the
#   estimated period surrounding that peak; after all points are masked, I find
#   typical spacing between peaks, take midpoints as minima
#
def selectPrefitData( fGHz, trans, phase, unc ) :
    period = findPeriod( fGHz, trans)					    # estimated period of transmission data, GHz
    tma = numpy.ma.asarray( trans )     # create masked array of transmission data

  # create list of maxima, one per period
    nlist = []
    while ( tma.count() > 0) :	        # loop until all elements are masked
      nmax = numpy.ma.argmax(tma)       # maximum remaining element
      nlist.append( nmax )              
      for i in range(0, len(tma) ) :
        if abs(fGHz[i] - fGHz[nmax]) < 0.75 * period :
          tma[i] = numpy.ma.masked
    nmaxlist = numpy.sort(nlist) 
    # print "selectPrefitData: nmaxlist =",nmaxlist

  # assume midpoints between maxima are close to minima
    if nmaxlist[-1] == len(trans) - 1 :
      deltan = int( (nmaxlist[-2]-nmaxlist[-3])/2. + 0.5)
    else :
      deltan = int( (nmaxlist[-1]-nmaxlist[-2])/2. + 0.5)
    # print "selectPrefitData: avg delta = %d" % deltan
    if (nmaxlist[0]-deltan) > 0 :
      nlist.append( nmaxlist[0] - deltan )
    for n in range(0,len(nmaxlist)-1) :
      nlist.append(nmaxlist[n] + deltan)
    if (nmaxlist[-1] + deltan) < len(tma) :
      nlist.append(nmaxlist[-1] + deltan)      
    
  # return sparse arrays
    nsorted = numpy.sort(nlist)
    print "selectPrefitData: nsorted =",nsorted
    fsparse = []
    tsparse = []
    psparse = []
    usparse = []
    for n in range(0,len(nsorted)) :
      fsparse.append( fGHz[nsorted[n]] )
      tsparse.append( trans[nsorted[n]] )
      psparse.append( phase[nsorted[n]] )
      usparse.append( unc[nsorted[n]] )
      # print fsparse[-1],tsparse[-1],psparse[-1],usparse[-1]
    return numpy.array(fsparse), numpy.array(tsparse), numpy.array(psparse), numpy.array(usparse)      

# used to create 5x finer grid around grid point with index i
# step size is 5x finer than calculated from fitParams
# returns new lists tcmList, nrList, tanDlist just around that point
def expandPoint( fitParams, index, tcmList, nrList, tanDList ) :
    nlayers = len(tcmList[index])
    tcmListNew = []
    nrListNew = []
    tanDListNew = []
    for nl in range(0,nlayers) :

      if len(fitParams["tcmList"][nl]) > 1 :
        step = fitParams["tcmList"][nl][1] - fitParams["tcmList"][nl][0]    # original step size
        tcmListNew.append( [ tcmList[index][nl]-0.4*step, tcmList[index][nl]-0.2*step, \
          tcmList[index][nl]+0.2*step, tcmList[index][nl]+0.4*step] )
      else :
        tcmListNew.append( fitParams["tcmList"][nl] )  # fixed value if positive, mirror if negative

      if len(fitParams["nrList"][nl]) > 1 :
        step = fitParams["nrList"][nl][1] - fitParams["nrList"][nl][0]
        nrListNew.append( [ nrList[index][nl]-0.4*step, nrList[index][nl]-0.2*step, \
          nrList[index][nl]+0.2*step, nrList[index][nl]+0.4*step] )
      else :
        nrListNew.append( fitParams["nrList"][nl] )   # fixed value if positive, mirror if negative

      if len(fitParams["tanDeltaList"][nl]) > 1 :
        step = fitParams["tanDeltaList"][nl][1] - fitParams["tanDeltaList"][nl][0]
        tanDListNew.append( [ tanDList[index][nl]-0.4*step, tanDList[index][nl]-0.2*step, \
          tanDList[index][nl]+0.2*step, tanDList[index][nl]+0.4*step] )
      else :
        tanDListNew.append( fitParams["tanDeltaList"][nl] )   # fixed value if positive, mirror if negative

    # print tcmListNew
    # print nrListNew
    # print tanDListNew
    return tcmListNew, nrListNew, tanDListNew 
    
def printElapsed( ntrials, i, nblock, secs0 ) :
    if i > 0 and i%nblock == 0 :
      secselapsed = time.time() - secs0
      secsremaining = (ntrials-i)*(secselapsed/float(i))
      if secsremaining/60. > 60. :
        print "... %d/%d finished; %.2f hours remaining" % (i,ntrials,secsremaining/3600.)
      else :
        print "... %d/%d finished; %.2f minutes remaining" % (i,ntrials,secsremaining/60.)

# nrefracFit handles joint fits to multiple data sets, each with its own angIdeg and freq coverage
#   (but currently it will not handle both heterodyne and FTS data)
# it uses a prefit: calculates chisq[0] for a very sparse dataset that samples just max and min values;
#   then does the calculation for the full dataset only for those parameters where 
#   chisq[0] < pCutoff * pchisqBest
# fitFile lists input files, parameter ranges to search etc; but angIdeg comes from individual files
#
def nrefracFit( fitFile, pCutoff=5., pdfOnly=False, Iterate=False, path="/o/plambeck/PolarBear/OpticsBench/" ) :

  # fitParams dictionary contains list of data files, parameter ranges to search, etc
    fitParams = readFitFile( fitFile )
    nDataSets = len( fitParams["datafileList"] )

  # read the data files; store as LISTS of arrays; do not append into one array in case angIdeg differs
    fGHz = []
    vecmeas = []
    unc = []
    angIdeg = []
    FTSdata = []    # True if FTS (power), False if heterodyne (amplitude)

    for datafile in fitParams["datafileList"] :
      print "\nreading data file %s" % datafile
      f,amp,phs,u,angI,FTS = readTransFile( datafile, path=path )     # read single data file

    # create sparse data set for prefit, probing peaks and valleys; save as an odd-numbered data set
      pf, pamp, pphs, pu =  selectPrefitData( f, amp, phs, u ) 
      fGHz.append( pf )
      vecmeas.append( pamp * numpy.exp(1j*numpy.radians(pphs)) )
         # note: although FTS data are scalars (power transmission) they are saved as vectors anyway
      unc.append( pu )
      angIdeg.append( angI )  
      FTSdata.append( FTS )

    # now append full data set to lists; these will be even-numbered
      fGHz.append(f)
      vecmeas.append( amp * numpy.exp(1j*numpy.radians(phs)) )  
      unc.append(u)
      angIdeg.append( angI )  
      FTSdata.append( FTS )

  # create tree of parameters to search
    tcmList,nrList,tanDeltaList = expandTree2( fitParams["tcmRange"], fitParams["tcmList"], fitParams["nrList"], fitParams["tanDeltaList"] )
      # tcmList = [ [tcm1,tcm2,...], [tcm1,tcm2,...], ... ] - expanded lists of thicknesses
      # nrList = [ [nr1,nr2,...], [nr1,nr2,...], ... ] - expanded lists of refractive indices
      # thus: { tcmlist[i], nrList[i], tanDeltaList[i] } specify one set of thicknesses, refractive indices, tanDs to be modeled
    ntrials = len( tcmList )
    nlayers = len( tcmList[0] )

  # one chisq array for prefit, one for fit to each dataset, one for joint fit
    chisq = 1.e6*numpy.ones( (nDataSets+2, ntrials) )

  # for prefit, model odd-numbered data sets only, calculate reduced chisq for every possible trial
    secs0 = time.time()       
    print "\nbeginning prefit"
    print "modeling %d layers, %d trials" % (nlayers, ntrials)
    for i in range(0,ntrials) :
      printElapsed( ntrials, i, 1000, secs0 )
      errs = []
      for j in range(0,2*nDataSets,2) :
        phs,amp = solveStack( fGHz[j], angIdeg[j], tcmList[i], nrList[i], tanDeltaList[i] ) 
        if FTSdata[j] :
          vecmodel = amp*amp
        else :
          vecmodel = amp * numpy.exp(1j*numpy.radians(phs))
        errs.append( (vecmeas[j]-vecmodel)/unc[j] )
      chisq[0,i] = numpy.var( numpy.concatenate(errs) )
    nbest = printBest( chisq[0], tcmList, nrList, tanDeltaList )      # nbest is INDEX of parameters with smallest chisq
    pchisqCutoff = pCutoff*chisq[0,nbest]

  # if we are iterating, define new trials around each trial with chisq < pCutoff
    if (Iterate) :
      for i in range(0,ntrials) :
        if chisq[0,i] < pchisqCutoff :
          tx,nx,tDx = expandPoint( fitParams, i, tcmList, nrList, tanDeltaList ) 
          tcmExtra,nrExtra,tanDextra = expandTree2( fitParams["tcmRange"], tx, nx, tDx )
        # append to existing list of trials (should check for duplicates here!)
          tcmList = numpy.concatenate( (tcmList, tcmExtra) )
          nrList = numpy.concatenate( (nrList, nrExtra) )
          tanDeltaList = numpy.concatenate( (tanDeltaList, tanDextra) )

    # compute number of new tests  
      ntrials2 = len(tcmList) - ntrials 
      print "\nfitting %d additional trials" % ntrials2
      chisqExtra = 1.e6*numpy.ones( (nDataSets+2, ntrials2) )
      chisq = numpy.concatenate( (chisq, chisqExtra), axis=1 )
      secs0 = time.time()       
      for i in range(ntrials+1, len(tcmList) ) :
        printElapsed( ntrials2, i-ntrials, 1000, secs0 )
        errs = []
        for j in range(0,2*nDataSets,2) :
          phs,amp = solveStack( fGHz[j], angIdeg[j], tcmList[i], nrList[i], tanDeltaList[i] ) 
          if FTSdata[j] :
            vecmodel = amp*amp
          else :
            vecmodel = amp * numpy.exp(1j*numpy.radians(phs))
          errs.append( (vecmeas[j]-vecmodel)/unc[j] )
        chisq[0,i] = numpy.var( numpy.concatenate(errs) )

    # find new best fits
      nbest = printBest( chisq[0], tcmList, nrList, tanDeltaList )      # nbest is INDEX of parameters with smallest chisq
      pchisqCutoff = pCutoff*chisq[0,nbest]

  # figure out how many trials will be needed for full test
    ntrials3 = bisect.bisect_left( numpy.sort(chisq[0]), pchisqCutoff )
        
  # now model full datasets, only for trials where pchisq < pchisqCutoff
    ntrials = len( tcmList )
    secs0 = time.time()       
    imodel = 0
    imodelLast = 0
    print "\nbeginning full test"
    print "full data set will be modeled for %d out of %d trials" % (ntrials3,ntrials)
    for i in range(0,ntrials) :
      if imodel > imodelLast :
        printElapsed( ntrials3, imodel, 200, secs0 )
        imodelLast = imodel
      if chisq[0,i] > pchisqCutoff :
        continue
      else :
        imodel = imodel+1
        errs = []
        for j in range(1, len(fGHz), 2) :
          phs,amp = solveStack( fGHz[j], angIdeg[j], tcmList[i], nrList[i], tanDeltaList[i] ) 
          if FTSdata[j] :
            vecmodel = amp*amp
          else :
            vecmodel = amp * numpy.exp(1j*numpy.radians(phs))
          errs.append( (vecmeas[j]-vecmodel)/unc[j] )
          chisq[(j+1)/2, i] = numpy.var(errs[-1])      # this is variance for fit to single data set only!
          # print "chisq[%d, %d] = %.3f" % ( (j+1)/2, i, numpy.var(errs[-1]))
      chisq[-1,i] = numpy.var( numpy.concatenate(errs) )   # this is variance for joint fit to all data
    nbest = printBest( chisq[-1], tcmList, nrList, tanDeltaList ) 

  # save the date and time (added 31-oct-2018)
    now = datetime.datetime.now()
    dateString = now.strftime("%Y-%m-%d   %H:%M")
        
  # save results for further analysis in a unique pickle file (don't overwrite old file)
    nseq = 1
    while os.path.isfile( fitFile + ".pickle.%d" % nseq) :
      nseq = nseq + 1
    pickleFile = fitFile + ".pickle.%d" % nseq
    fout = open( pickleFile, "wb" )
    pickle.dump( [fitParams,tcmList,nrList,tanDeltaList,chisq,dateString], fout ) 
      # note: angIdeg was used in fit, but it is stored in each data file
    fout.close() 

  # plot the fit
    time.sleep(1)
    plotFit( pickleFile, pdfOnly=pdfOnly, path=path )

def plotFit( pickleFile, FTS=False, pdfOnly=False, path="/o/plambeck/PolarBear/OpticsBench/" ) :

  # read fit results from pickleFile
    print "loading fit results from file %s" % pickleFile

  # dateString was added on 31-oct-2018; must be able to handle earlier files without it
    try :
      fin = open( pickleFile, "rb" )
      [fitParams,tcmList,nrList,tanDeltaList,chisq,dateString] = pickle.load( fin )
    except :
      fin = open( pickleFile, "rb" )
      [fitParams,tcmList,nrList,tanDeltaList,chisq] = pickle.load( fin )
      dateString = " "
    fin.close() 
    nlayers = len( tcmList[0] )

    print "\n10 best fits minimizing overall residuals:"
    nbest = printBest( chisq[-1], tcmList, nrList, tanDeltaList ) 

  # create savet,saven,savetanD lists, containing all configurations with chisq< 2.* minchisq
    savet = []
    savenr = []
    savetanD = []
    for n in range(0,len(chisq[-1])) :
      if chisq[-1,n] < 2.*chisq[-1,nbest]:
        savet.append( tcmList[n] )
        savenr.append( nrList[n] )
        savetanD.append( tanDeltaList[n] )

  # begin the plot
    pyplot.ioff()
    #pyplot.clf()
    pp = PdfPages( pickleFile + ".pdf")

  # one page per dataset
    for datafile in fitParams["datafileList"] :
      fig =  pyplot.figure( figsize=(11,8) )
      pyplot.suptitle( "%s   %s   %s" % (fitParams["title"],pickleFile,datafile), fontsize=12 )

    # retrieve the raw data; also, identify points used for prefit
      fGHz,trans,dphs_meas,unc,angIdeg,FTS = readTransFile( datafile, path=path )
      pf, pamp, pphs, pu =  selectPrefitData( fGHz, trans, dphs_meas, unc ) 

    # recalculate the model amps and phases for the best fit parameters
      dphs_est,trans_est = solveStack( fGHz, angIdeg, tcmList[nbest], nrList[nbest], tanDeltaList[nbest] ) 

    # plot phase unless FTS
      xmin,xmax = minmax( fGHz )
      if not FTS :
        dphi = unwrap( dphs_meas-dphs_est )

      # top left panel is phase vs freq, measured and fit
        ax = fig.add_axes( [.1,.6,.38,.32] )
        ymin,ymax = minmax( numpy.concatenate((dphs_meas,dphs_est)))
        ax.axis( [xmin, xmax, ymin, ymax] )
        ax.tick_params( labelsize=10 )
        ax.plot( fGHz, dphs_meas, "bo", ms=2. )
        ax.plot( pf, pphs, "kx", ms=5. )
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
      if FTS :
        trans_est = trans_est*trans_est
        ax = fig.add_axes( [.1,.6,.85,.32] )
        ax.set_ylabel("transmitted pwr", fontsize=10)
      else :
        ax = fig.add_axes( [.57,.6,.38,.32] )
        ax.set_ylabel("transmitted amplitude", fontsize=10)
      ymin,ymax = minmax( numpy.concatenate((trans,trans_est)) )
      ax.axis( [xmin, xmax, ymin, ymax] )
      if numpy.mean(unc) > .5 :
        ax.plot( fGHz, trans, "o", ms=2. )
      else :
        ax.errorbar( fGHz, trans, yerr=unc, elinewidth=0.5, capsize=0 )
      ax.plot( pf, pamp, "kx", ms=5. )
      ax.plot( fGHz, trans_est, "r-" )
      ax.tick_params( labelsize=10 )
      if FTS :
        ax.set_ylim( [.3,1.05] )
      else :
        ax.set_ylim( [.5,1.05] )
      pyplot.grid(True)

    # middle right panel is measured-theoretical trans
      if FTS :
        ax = fig.add_axes( [.1,.37,.85,.2] )
      else :
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
      ax.text( 0.03, ylab, "angle = %.2f deg" % angIdeg, transform=ax.transAxes, \
        horizontalalignment='left', fontsize=fontSize, color="black", rotation='horizontal' )
      ylab = ylab - ystep 
      ax.text( 0.03, ylab, "allowed thickness = [%s]" % makeString( fitParams["tcmRange"]/2.54 ), transform=ax.transAxes, \
        horizontalalignment='left', fontsize=fontSize, color="black", rotation='horizontal' )
      tcmListIn = fitParams["tcmList"]
      nrListIn = fitParams["nrList"]
      tanDeltaListIn = fitParams["tanDeltaList"]
      ylab = ylab - ystep 
      ax.text( 0.00, ylab, "search ranges [tinches] [n] [tanDelta]:", fontsize=fontSize)
      for nl in range(0,nlayers) : 
        ylab = ylab - ystep 
        if nrListIn[nl][0] < 0. :
          ax.text( 0.03, ylab, "%d: mirror layer %d" % (nl+1, -1*nrListIn[nl][0]), \
            transform=ax.transAxes, horizontalalignment='left', fontsize=fontSize, \
            color="black", rotation='horizontal' )
        else :
          tstring = makeString( tcmListIn[nl]/2.54 )
          nstring = makeString( nrListIn[nl], fmt=4 )
          tanDstring = makeString( tanDeltaListIn[nl] )
          ax.text( 0.03, ylab, "%d: [%s]" % (nl+1,tstring), \
            transform=ax.transAxes, horizontalalignment='left', fontsize=fontSize, \
            color="black", rotation='horizontal' )
          ylab = ylab - ystep 
          ax.text( 0.03, ylab, "    [%s]" % (nstring), \
            transform=ax.transAxes, horizontalalignment='left', fontsize=fontSize, \
            color="black", rotation='horizontal' )
          ylab = ylab - ystep 
          ax.text( 0.03, ylab, "    [%s]" % (tanDstring), \
            transform=ax.transAxes, horizontalalignment='left', fontsize=fontSize, \
            color="black", rotation='horizontal' )
      pyplot.axis('off')

    # fictitious lower right panel gives room for text
      ax = fig.add_axes( [.57,.1,.38,.2])
      ylab = 1.-ystep
      ax.text( 0.00, ylab, "fit results:", fontsize=fontSize)
      ylab = ylab - ystep
      ax.text( 0.03, ylab, "reduced chisq = %.3f " % chisq[-1,nbest], transform=ax.transAxes, \
        horizontalalignment='left', fontsize=fontSize, color="black", rotation='horizontal' )
      ylab = ylab - ystep 
      for nl in range(0,nlayers) :
        if nrListIn[nl][0] > 0. :
          tmin,tmax = tLimits( savet, nl)
          ax.text( 0.03, ylab, "layer %d:  t = %.4f\"   [%.4f, %.4f\"]" % \
            ( nl+1, tcmList[nbest][nl]/2.54, tmin/2.54, tmax/2.54) ,  \
            transform=ax.transAxes, horizontalalignment='left', fontsize=fontSize, \
            color="black", rotation='horizontal' )
          ylab = ylab - ystep
          nrmin,nrmax = tLimits( savenr, nl )
          ax.text( 0.03, ylab, "              n = %.4f  [%.4f, %.4f]" % \
            ( nrList[nbest][nl], nrmin, nrmax ), \
            transform=ax.transAxes, horizontalalignment='left', fontsize=fontSize, \
            color="black", rotation='horizontal' )
          ylab = ylab - ystep 
          ax.text( 0.03, ylab, "              e = %.3f  [%.3f, %.3f]" % \
            ( pow(nrList[nbest][nl],2), pow(nrmin,2), pow(nrmax,2) ),  \
            transform=ax.transAxes, horizontalalignment='left', fontsize=fontSize, \
            color="black", rotation='horizontal' )
          ylab = ylab - ystep 
          tanDmin,tanDmax = tLimits( savetanD, nl )
          ax.text( 0.03, ylab, "              tanD = %.4f  [%.4f, %.4f]" % \
            ( tanDeltaList[nbest][nl], tanDmin, tanDmax ), \
            transform=ax.transAxes, horizontalalignment='left', fontsize=fontSize, \
            color="black", rotation='horizontal' )
          ylab = ylab - ystep
      pyplot.axis('off')
    
    # print date and time in lower right corner
      ax = fig.add_axes( [.1,.05,.85,.1])
      ax.text(.5, 0., "%s" % dateString, \
            transform=ax.transAxes, horizontalalignment='center', fontsize=fontSize, \
            color="black", rotation='horizontal' )
      pyplot.axis('off')

      pyplot.savefig( pp, format="pdf" )
      if not pdfOnly :
        pyplot.show()
      else :
        print "plot saved to file %s" % (pickleFile + ".pdf")
    pp.close()

def makeString4( srchmin, low, best, high, srchmax) :
    lowstring = "   ----   "
    if low > srchmin :
       lowstring = "%6.4f" % low
    histring = "    ----    "
    if high < srchmax :
       histring = "%6.4f" % high
    return "%s   %6.4f   %s" % (lowstring,best,histring)
    
  
# new version of plot that combines search ranges and results for each layer onto 1 line
def plotFit2( pickleFile, FTS=False, pdfOnly=False, path="/o/plambeck/PolarBear/OpticsBench/" ) :

  # read fit results from pickleFile
    print "loading fit results from file %s" % pickleFile

  # dateString was added on 31-oct-2018; must be able to handle earlier files without it
    try :
      fin = open( pickleFile, "rb" )
      [fitParams,tcmList,nrList,tanDeltaList,chisq,dateString] = pickle.load( fin )
    except :
      fin = open( pickleFile, "rb" )
      [fitParams,tcmList,nrList,tanDeltaList,chisq] = pickle.load( fin )
      dateString = " "
    fin.close() 
    nlayers = len( tcmList[0] )

    print "\n10 best fits minimizing overall residuals:"
    nbest = printBest( chisq[-1], tcmList, nrList, tanDeltaList ) 

  # create savet,saven,savetanD lists, containing all configurations with chisq< 2.* minchisq
    savet = []
    savenr = []
    savetanD = []
    for n in range(0,len(chisq[-1])) :
      if chisq[-1,n] < 2.*chisq[-1,nbest]:
        savet.append( tcmList[n] )
        savenr.append( nrList[n] )
        savetanD.append( tanDeltaList[n] )

  # begin the plot
    pyplot.ioff()
    #pyplot.clf()
    pp = PdfPages( pickleFile + ".pdf")

  # one page per dataset
    for datafile in fitParams["datafileList"] :
      fig =  pyplot.figure( figsize=(11,8) )
      pyplot.suptitle( "%s   %s   %s" % (fitParams["title"],pickleFile,datafile), fontsize=12 )

    # retrieve the raw data; also, identify points used for prefit
      fGHz,trans,dphs_meas,unc,angIdeg,FTS = readTransFile( datafile, path=path )
      pf, pamp, pphs, pu =  selectPrefitData( fGHz, trans, dphs_meas, unc ) 

    # recalculate the model amps and phases for the best fit parameters
      dphs_est,trans_est = solveStack( fGHz, angIdeg, tcmList[nbest], nrList[nbest], tanDeltaList[nbest] ) 

    # plot phase unless FTS
      xmin,xmax = minmax( fGHz )
      if not FTS :
        dphi = unwrap( dphs_meas-dphs_est )

      # top left panel is phase vs freq, measured and fit
        ax = fig.add_axes( [.1,.6,.38,.32] )
        ymin,ymax = minmax( numpy.concatenate((dphs_meas,dphs_est)))
        ax.axis( [xmin, xmax, ymin, ymax] )
        ax.tick_params( labelsize=10 )
        ax.plot( fGHz, dphs_meas, "bo", ms=2. )
        ax.plot( pf, pphs, "kx", ms=5. )
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
      if FTS :
        trans_est = trans_est*trans_est
        ax = fig.add_axes( [.1,.6,.85,.32] )
        ax.set_ylabel("transmitted pwr", fontsize=10)
      else :
        ax = fig.add_axes( [.57,.6,.38,.32] )
        ax.set_ylabel("transmitted amplitude", fontsize=10)
      ymin,ymax = minmax( numpy.concatenate((trans,trans_est)) )
      ax.axis( [xmin, xmax, ymin, ymax] )
      if numpy.mean(unc) > .5 :
        ax.plot( fGHz, trans, "o", ms=2. )
      else :
        ax.errorbar( fGHz, trans, yerr=unc, elinewidth=0.5, capsize=0 )
      ax.plot( pf, pamp, "kx", ms=5. )
      ax.plot( fGHz, trans_est, "r-" )
      ax.tick_params( labelsize=10 )
      if FTS :
        ax.set_ylim( [.3,1.05] )
      else :
        ax.set_ylim( [.5,1.05] )
      pyplot.grid(True)

    # middle right panel is measured-theoretical trans
      if FTS :
        ax = fig.add_axes( [.1,.37,.85,.2] )
      else :
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

    # (fictitious) bottom panel lists search ranges, limits, best values
      fontSize = 12
      ystep = .12
      ax = fig.add_axes( [.1,.1,.85,.2])
      ylab = 1 - ystep
      ax.text( 0.0, ylab, "angle = %.2f deg; allowed thickness = [%s]" % \
        ( angIdeg, makeString(fitParams["tcmRange"]/2.54) ), transform=ax.transAxes, \
        horizontalalignment='left', fontsize=fontSize, color="black", rotation='horizontal' )
      ax.text( 1.0, ylab, "reduced chisq = %.3f " % chisq[-1,nbest], transform=ax.transAxes, \
        horizontalalignment='right', fontsize=fontSize, color="black", rotation='horizontal' )

      tcmListIn = fitParams["tcmList"]
      nrListIn = fitParams["nrList"]
      tanDeltaListIn = fitParams["tanDeltaList"]

      ylab = ylab - ystep - .05
      ax.text( 0.03, ylab, "layer", style='oblique', horizontalalignment='center')
      ax.text( 0.18, ylab, "thickness", style='oblique', horizontalalignment='center')
      ax.text( 0.42, ylab, "index", style='oblique', horizontalalignment='center')
      ax.text( 0.65, ylab, "epsilon", style='oblique', horizontalalignment='center')
      ax.text( 0.87, ylab, "tanDelta", style='oblique', horizontalalignment='center')
      ylab = ylab - 0.02
     
      for nl in range(0,nlayers) : 
        ylab = ylab - ystep 
        if nrListIn[nl][0] < 0. :
          ax.text( 0.0, ylab, "%d: mirror layer %d" % (nl+1, -1*nrListIn[nl][0]), \
            transform=ax.transAxes, horizontalalignment='left', fontsize=fontSize, \
            color="black", rotation='horizontal' )
        else :
          ax.text( 0.03, ylab, "%d" % (nl+1), horizontalalignment='center')

          tmin,tmax = tLimits( savet, nl)
          lowString = "  ----   "
          if tcmListIn[nl][0] < tmin :
            lowString = "%.4f" % (tmin/2.54)
          highString = "  ----   "
          if tcmListIn[nl][-1] > tmax :
            highString = "%.4f" % (tmax/2.54)
          ax.text( 0.105, ylab, "%s" % lowString, horizontalalignment='center') 
          ax.text( 0.18, ylab, "%.4f" % (tcmList[nbest][nl]/2.54), color='blue', horizontalalignment='center') 
          ax.text( 0.255, ylab, "%s" % highString, horizontalalignment='center') 

          nrmin,nrmax = tLimits( savenr, nl )
          lowString = "  ----  "
          epslowString = "  ----  "
          if nrListIn[nl][0] < nrmin :
            lowString = "%.3f" % nrmin
            epslowString = "%.3f" % (nrmin*nrmin)
          highString = "  ----  "
          epshighString = "  ----  "
          if nrListIn[nl][-1] > nrmax :
            highString = "%.3f" % nrmax
            epshighString = "%.3f" % (nrmax*nrmax)

          ax.text( 0.355, ylab, "%s" % lowString, horizontalalignment='center') 
          ax.text( 0.42, ylab, "%.3f" % (nrList[nbest][nl]), color='blue', horizontalalignment='center') 
          ax.text( 0.485, ylab, "%s" % highString, horizontalalignment='center') 
 
          ax.text( 0.585, ylab, "%s" % epslowString, horizontalalignment='center') 
          ax.text( 0.65, ylab, "%.3f" % (pow(nrList[nbest][nl],2)), color='blue', horizontalalignment='center') 
          ax.text( 0.715, ylab, "%s" % epshighString, horizontalalignment='center') 


          tanDmin,tanDmax = tLimits( savetanD, nl )
          lowString = "  ----   "
          if tanDeltaListIn[nl][0] < tanDmin :
            lowString = "%.4f" % tanDmin
          highString = "  ----   "
          if tanDeltaListIn[nl][-1] > tmax :
            highString = "%.4f" % tanDmax
          ax.text( 0.8, ylab, "%s" % lowString, horizontalalignment='center') 
          ax.text( 0.87, ylab, "%.4f" % (tanDeltaList[nbest][nl]), color='blue', horizontalalignment='center') 
          ax.text( 0.94, ylab, "%s" % highString, horizontalalignment='center') 


          
          #ax.text( 0.01, ylab, "%d:     %s      %s      %s"  % \
          #  (nl+1, makeString4(tcmListIn[nl][0]/2.54, tmin/2.54, tcmList[nbest][nl]/2.54, tmax/2.54, tcmListIn[nl][-1]/2.54), \
          #  makeString4(nrListIn[nl][0], nrmin, nrList[nbest][nl], nrmax, nrListIn[nl][-1]), \
          #  makeString4(tanDeltaListIn[nl][0], tanDmin, tanDeltaList[nbest][nl], tanDmax, tanDeltaListIn[nl][-1]) ), \
          #  transform=ax.transAxes, horizontalalignment='left', fontsize=fontSize, \
          #  color="black", rotation='horizontal' )

         

          #ax.text( 0.2, ylab, "%s" % \
          #  makeString4(tcmListIn[nl][0]/2.54, tmin/2.54, tcmList[nbest][nl]/2.54, tmax/2.54, tcmListIn[nl][-1]/2.54), \
          #  transform=ax.transAxes, horizontalalignment='center', fontsize=fontSize, \
          #  color="black", rotation='horizontal' )


      pyplot.axis('off')

    # print date and time at the bottom
      ax = fig.add_axes( [.1,.05,.85,.1])
      ax.text(.5, 0., "%s" % dateString, \
            transform=ax.transAxes, horizontalalignment='center', fontsize=fontSize, \
            color="black", rotation='horizontal' )
      pyplot.axis('off')

      pyplot.savefig( pp, format="pdf" )
      if not pdfOnly :
        pyplot.show()
      else :
        print "plot saved to file %s" % (pickleFile + ".pdf")
    pp.close()

# dump fit results to csv file to make it easier to fill out google sheet
# chisq array includes fits to individual data sets and to joint fit;
#   use this to tabulate best fit for each dataset and for the joint fit
#
def tabulateResults( pickleFile, outFile="summary.csv" ) :

  # read fit results from pickleFile
    print "loading fit results from file %s" % pickleFile
    try :
      fin = open( pickleFile, "rb" )
      [fitParams,tcmList,nrList,tanDeltaList,chisq,dateString] = pickle.load( fin )
    except :
      fin = open( pickleFile, "rb" )
      [fitParams,tcmList,nrList,tanDeltaList,chisq] = pickle.load( fin )
      dateString = " "
    fin.close() 
    nlayers = len( tcmList[0] )
 
  # print fit results, uncertainties, and chisq for each dataset individually 
    fout = open( outFile, "a" )
    for j in range(1,len(chisq)):

    # for fits to a single data set, don't write out joint fit results, as they are the same as individual fit
      if (len(chisq) == 3) and (j == 2) :
        break

      name = "joint fit"
      if j < len(chisq)-1 :
        name = fitParams["datafileList"][j-1].strip(".avg.dat")

      print "\n====="
      print "\n10 best fits for %s" % name
      nbest = printBest( chisq[j], tcmList, nrList, tanDeltaList ) 

      best,tcmFit,nrFit,tanDfit = summarize( chisq[j], tcmList, nrList, tanDeltaList )
      fout.write("\n%s,%s,%s,%.3f," % (pickleFile,fitParams["title"], name, best ) )   # col 1 is puck name, col 2 is dataFile name
      print best
      print tcmFit/2.54
      print nrFit
      print tanDfit
      for nl in range(0,nlayers) :
        print fitParams["nrList"][nl][0]
        print nrFit[nl,0],nrFit[nl,1],nrFit[nl,2]
        print fitParams["nrList"][0][-1]
        fout.write(", %.4f,%.4f,%.4f,%.4f,%.4f," % (fitParams["tcmList"][nl][0]/2.54, tcmFit[nl,0]/2.54, \
          tcmFit[nl,1]/2.54, tcmFit[nl,2]/2.54, fitParams["tcmList"][nl][-1]/2.54 ) )
        fout.write(", %.4f,%.4f,%.4f,%.4f,%.4f," % (fitParams["nrList"][nl][0], nrFit[nl,0], \
          nrFit[nl,1], nrFit[nl,2], fitParams["nrList"][nl][-1]) )
        fout.write(", %.4f,%.4f,%.4f,%.4f,%.4f," % (fitParams["tanDeltaList"][nl][0], tanDfit[nl,0], \
          tanDfit[nl,1], tanDfit[nl,2], fitParams["tanDeltaList"][nl][-1]) )
    fout.write("\n")
    fout.close()

# vary thickness of each layer over a range, predict transmission curve
# fitFile contains total thickness range of each layer; nr and tanDelta should be single-valued
# fuzzy1 does (weighted) vector average of ALL allowed thicknesses to get final transmission amp and phase;
def fuzzy1( fitFile, fGHz=numpy.arange(70.,115.1,.1), angIdeg=5.16 ) :
    fitParams = readFitFile( fitFile )
    tcmList,nrList,tanDeltaList = expandTree2( fitParams["tcmRange"], fitParams["tcmList"], fitParams["nrList"], fitParams["tanDeltaList"] )
    ntrials = len(tcmList)
    print "averaging results of %d trials" % ntrials
    vecmodel = numpy.zeros( len(fGHz), dtype=complex )
    for n in range(0,ntrials) :
      phs,amp = solveStack( fGHz, angIdeg, tcmList[n], nrList[n], tanDeltaList[n] ) 
      vecmodel = vecmodel + amp * numpy.exp(1j*numpy.radians(phs))
    vecmodel = vecmodel/float(ntrials)
    amplitude = numpy.absolute( vecmodel )
    phase = numpy.angle( vecmodel, deg=True )
    pwr = numpy.power( amplitude, 2 )
    print "average power transmission = %.3f" % (numpy.average(pwr))

  # copy fitFile parameters to output file
    fin = open( fitFile+".fitFile", "r")
    fout = open( fitFile+".model", "w" )
    fout.write("# %s.model\n" % fitFile)
    fout.write("# ---------------------------------------------------- \n")
    for line in fin:
      fout.write("# %s" % line)
    fout.write("# ---------------------------------------------------- \n")
    fout.write("# averaged results of %d trials\n" % ntrials)
    fout.write("# average power transmission = %.3f\n#\n" % (numpy.average(pwr)) )

  # now write out the model data, in format compatible with readTransData
    unc = .02
    for i in range(0,len(fGHz)) :
      fout.write("%8.2f  %8.6f  %8.3f  %6.4f   1  %7.3f\n" % (fGHz[i], amplitude[i], phase[i], unc, pwr[i] ) )
    fout.close()

# fuzzy2 creates model transmission data for the case where the surfaces are rough on
#   scales less than a wavelength, so that there is a gradation in the refractive index
#   and loss tangent from one layer to the next
# normally the input fitFile should specify exactly one model, with fixed thicknesses, 
#   refractive indices, and tanDeltas; the thicknesses are interpreted to be the nominal
#   thicknesses of the layers, without accounting for the transition regions; that is,
#   the total thicknesses of the pure layers should add up to be the total puck thickness
# deltaTinList specifies the thickness of each transition (e.g., 0 for an abrupt transition,
#   0.002 for transition with .002" p-p variation, etc); note that for n layers there must
#   be (n+1) transitions specified
#
def fuzzy2( fitFile, angIdeg=5.16 ) :
    #fGHz = numpy.concatenate( (numpy.arange(76.,115.,.2),numpy.arange(205.,230.,.2)) )
    #fGHz = numpy.arange(210.,219.6,.1)
    #fGHz = numpy.arange(70.,230.1,.2)
    fGHz = numpy.arange(110.,240.,1.)
    fitParams = readFitFile( fitFile )
    tcmList,nrList,tanDeltaList = expandTree2( fitParams["tcmRange"], fitParams["tcmList"], fitParams["nrList"], fitParams["tanDeltaList"] )
    ntrials = len(tcmList)
    if len(fitParams["deltaTinList"]) == len(tcmList[0])+1 :
      print "using deltaTinList =", fitParams["deltaTinList"]
      deltaTcmList = 2.54*numpy.array( fitParams["deltaTinList"] )
    else :
      deltaTcmList = numpy.zeros( len(tcmList[0])+1 )
      print "using deltaTinList =", deltaTcmList/2.54
    print "averaging results of %d trials" % ntrials
    vecmodel = numpy.zeros( len(fGHz), dtype=complex )
    for n in range(0,ntrials) :
      [tcmListSm, nrListSm, tanDeltaListSm ] = smoothTransitions( tcmList[n], nrList[n], tanDeltaList[n], deltaTcmList, fitFile )
      phs,amp = solveStack( fGHz, angIdeg, tcmListSm, nrListSm, tanDeltaListSm ) 
      vecmodel = vecmodel + amp * numpy.exp(1j*numpy.radians(phs))
    vecmodel = vecmodel/float(ntrials)
    amplitude = numpy.absolute( vecmodel )
    phase = numpy.angle( vecmodel, deg=True )
    pwr = numpy.power( amplitude, 2 )
    print "average power transmission = %.3f" % (numpy.average(pwr))

  # copy fitFile parameters to output file
    fin = open( fitFile+".fitFile", "r")
    fout = open( fitFile+".model", "w" )
    fout.write("# %s.model\n" % fitFile)
    fout.write("# ---------------------------------------------------- \n")
    for line in fin:
      fout.write("# %s" % line)
    fout.write("# angIdeg = %.2f\n" % angIdeg)
    fout.write("# ---------------------------------------------------- \n")
    fout.write("# average power transmission = %.3f\n#\n" % (numpy.average(pwr)) )
    for i in range(0,len(fGHz)) :
      fout.write("%9.3f  %8.6f  %8.3f  .01   %8.3f\n" % (fGHz[i], amplitude[i], phase[i], pwr[i] ) )
    fout.close()
  # compute power transmission averaged over band

# smoothTransitions expands the number of layers in each stack to model gradual changes
#  in the index of refraction - due, for example, to surface roughness on a scale much
#  less than a wavelength; it interpolates the refractive index and tanDelta through
#  these transition layers
# deltaTcmList is a list of the transition layer thicknesses, including the interfaces
#  to air on either side of the puck; so its length should be len(tcmList) + 1;
# tcmList[m] is the NOMINAL thickness of layer m; smoothTransitions reduces this
#  thickness appropriately to account for the transition layers
# to aid in visualizing the variation in refractive index, this routine writes out a
#  "profile" of refractive index vs distance through the puck
#
def smoothTransitions( tcmList, nrList, tanDeltaList, deltaTcmList, fitFile, stepsizeIn=.0005 ) :
    stepsize = stepsizeIn*2.54
    tcmListSm = []  
    nrListSm = []
    tanDeltaListSm = []
    nrPrev = 1.
    tdPrev = 0.
    xcm = 0.     # begin profile at xcm
    fout = open(fitFile+".profile","w")
    fout.write("%9.6f %9.6f %9.6f\n" % (xcm/2.54,nrPrev,0.))

  # reduce thicknesses of "pure" layers to allow for transition regions
    for n in range(0,len(tcmList)) :
      tcmList[n] = tcmList[n] - deltaTcmList[n]/2. 
      tcmList[n] = tcmList[n] - deltaTcmList[n+1]/2. 
    tcmList[0] = tcmList[0] - deltaTcmList[0]/2.
    tcmList[-1] = tcmList[-1] - deltaTcmList[-1]/2.

  # now work through the stack, to create tcmListsm, etc.
    for n in range(0,len(tcmList)) :

    # deal with the transition at the beginning of each layer 
      if deltaTcmList[n] > 0. :
        nsteps = int( round(deltaTcmList[n]/stepsize) )
        deltastep = deltaTcmList[n]/nsteps       # deltastep is thickness of each transition layer
        nrstep = (nrList[n]-nrPrev)/(nsteps+1)
        tdstep = (tanDeltaList[n] - tdPrev)/(nsteps+1)
        for ns in range(1, nsteps+1) :
          tcmListSm.append( deltastep )
          nrListSm.append( nrPrev + ns*nrstep )
          tanDeltaListSm.append( tdPrev + ns*tdstep )
          fout.write("%9.6f %9.6f %9.6f\n" % (xcm/2.54,nrListSm[-1],tanDeltaListSm[-1]))
             # this is the beginning of the step
          xcm = xcm + tcmListSm[-1]
          fout.write("%9.6f %9.6f %9.6f\n" % (xcm/2.54,nrListSm[-1],tanDeltaListSm[-1]))
             # this is the end of the step

    # put in the "pure" layer
      tcmListSm.append( tcmList[n] )
      nrListSm.append( nrList[n] )
      tanDeltaListSm.append( tanDeltaList[n] )
      fout.write("%9.6f %9.6f %9.6f\n" % (xcm/2.54,nrListSm[-1],tanDeltaListSm[-1]))
         # this is the beginning of the step
      xcm = xcm + tcmListSm[-1]
      fout.write("%9.6f %9.6f %9.6f\n" % (xcm/2.54,nrListSm[-1],tanDeltaListSm[-1]))
         # this is the end of the step
      nrPrev = nrList[n]
      tdPrev = tanDeltaList[n]

  # final transition back to air
    if deltaTcmList[-1] > 0. :
      nsteps = int( round(deltaTcmList[-1]/stepsize) )
      deltastep = deltaTcmList[-1]/nsteps
      nrstep = (1.-nrPrev)/(nsteps+1)
      tdstep = -1.*tdPrev/(nsteps+1)
      for ns in range(1, nsteps+1) :
        tcmListSm.append( deltastep )
        nrListSm.append( nrPrev + ns*nrstep )
        tanDeltaListSm.append( tdPrev + ns*tdstep )
        fout.write("%9.6f %9.6f %9.6f\n" % (xcm/2.54,nrListSm[-1],tanDeltaListSm[-1]))
        xcm = xcm + tcmListSm[-1]
        fout.write("%9.6f %9.6f %9.6f\n" % (xcm/2.54,nrListSm[-1],tanDeltaListSm[-1]))
    fout.write("%9.6f %9.6f %9.6f\n" % (xcm/2.54,1.,0.))
    fout.close()
    return [tcmListSm, nrListSm, tanDeltaListSm ]

# for 2 layers only; make 3Dplot of chisq for nr1,nr2 above the (tin1,tin2) plane
def my3Dplot( pickleFile, chisqIndex=-1, worstRatio=10. ) :

  # read fit results from pickleFile
    print "loading fit results from file %s" % pickleFile
    fin = open( pickleFile, "rb" )
    [fitParams,tcmList,nrList,tanDeltaList,chisq] = pickle.load( fin )
    fin.close() 
    ntrials = len( tcmList )
    nlayers = len( tcmList[0] )

  # this routine is specialized to just 2 layers
    if nlayers != 2 :
      print "nlayers = %d; this routine is specialized to exactly 2 layers; exiting" % nlayers
      return

    print "\n10 best fits minimizing overall residuals:"
    nbest = printBest( chisq[-1], tcmList, nrList, tanDeltaList ) 

    xplot = []
    yplot = []
    zplot = []
    c = []
    for n in range(0,ntrials) :
      if chisq[1][n] < worstRatio*chisq[chisqIndex][nbest] :
        xplot.append( tcmList[n][0]/2.54 )
        yplot.append( tcmList[n][1]/2.54 )
        zplot.append( nrList[n][0] )
        c.append( chisq[chisqIndex][n] )
    fig = pyplot.figure()
    ax = Axes3D(fig)
    ax.scatter3D( xplot, yplot, zplot, c=c, marker="o", alpha=0.6, edgecolors="none")
    #for n in range(0,ntrials) :
    #  if chisq[1][n] < 10.*chisq[0][nbest] :
    #    xplot.append( tcmList[n][0]/2.54 )
    #    yplot.append( tcmList[n][1]/2.54 )
    #    zplot.append( nrList[n][1] )
    #    c.append( chisq[chisqIndex][n] )
    #ax.scatter3D( xplot, yplot, zplot, c=c, marker="o", alpha=0.6, edgecolors="none")
    ax.view_init(elev=15,azim=160)
    pyplot.show()

# for 2 layers only; make 3Dplot of for nr1,nr2 above the (tin1,tin2) plane, for the smallest chisq for each tin1,tin2
def my3Dplot2( pickleFile, chisqIndex=-1, worstRatio=10. ) :

  # read fit results from pickleFile
    print "loading fit results from file %s" % pickleFile
    fin = open( pickleFile, "rb" )
    [fitParams,tcmList,nrList,tanDeltaList,chisq] = pickle.load( fin )
    fin.close() 
    ntrials = len( tcmList )
    nlayers = len( tcmList[0] )

  # this routine is specialized to just 2 layers
    if nlayers != 2 :
      print "nlayers = %d; this routine is specialized to exactly 2 layers; exiting" % nlayers
      return
  
  # sort by chisq
    isort = numpy.argsort( chisq[chisqIndex] )
    print isort[0:10]
    chisqMin = chisq[chisqIndex][isort[0]]
    nbest = printBest( chisq[chisqIndex], tcmList, nrList, tanDeltaList ) 
    print "chisqMin = %.3f" % chisqMin

    tcm0 = []
    tcm1 = []
    nr0 = []
    nr1 = []
    chisqNorm = []
    npts = 0
    for i in isort :
      npts = npts + 1
      if chisq[chisqIndex][i] > worstRatio * chisqMin :
        print "%d / %d points meet chisq criterion" % (npts, len(isort) )
        break
      else :
        newPoint = True
        for j in range(0,len(tcm0)) :
          if ( tcmList[i][0] == tcm0[j] ) and ( tcmList[i][1] == tcm1[j] ) :
            newPoint = False
        if newPoint :
          tcm0.append( tcmList[i][0] ) 
          tcm1.append( tcmList[i][1] ) 
          nr0.append( nrList[i][0] )
          nr1.append( nrList[i][1] )
          #chisqNorm.append( chisq[chisqIndex][i]/ chisqMin )
          chisqNorm.append( chisq[chisqIndex][i])
          print i, tcm0[-1]/2.54, tcm1[-1]/2.54, chisqNorm[-1]

    print "will be plotting %d points" % ( len(tcm0) )
    
    tin0 = numpy.array(tcm0)/2.54
    tin1 = numpy.array(tcm1)/2.54

    fig = pyplot.figure()
    ax = Axes3D(fig)
    ax.scatter3D( tin0, tin1, nr0, c=chisqNorm, marker="o", alpha=0.6, edgecolors="none")
    #ax.plot_trisurf( tin0, tin1, nr0, alpha=0.1  )
    #ax.scatter3D( tin0, tin1, nr1, c=chisqNorm, marker="o", alpha=0.6, edgecolors="none")
    ax.view_init(elev=15,azim=160)
    pyplot.show()

def testExpand( fitFile, index) :
    fitParams = readFitFile( fitFile )
    tcmList,nrList,tanDeltaList = expandTree2( fitParams["tcmRange"], fitParams["tcmList"], fitParams["nrList"], fitParams["tanDeltaList"] )
    print "\nstart with:", tcmList[index], nrList[index], tanDeltaList[index]
    print "\nexpand to:"
    expandPoint( fitParams, index, tcmList, nrList, tanDeltaList ) 

# designed to find most likely freq offset between heterodyne and FTS data
# resamples FTS data to heterodyne resolution using interp1d, covering a slightly wider freq range than heterodyne data;
#   then find what offset (in channels) gives highest crosscorrelation, convert this to frequency
#
def frqOffset( FTSfile, hetfile, fstart=205.2, fstop=229.8, fstep=0.05 ) :
    fGHz1,trans1,dphs_meas,unc,angIdeg,FTS = readTransFile( hetfile+".avg.dat", path="/o/plambeck/PolarBear/OpticsBench/" )
    fGHz2,trans2,dphs_meas,unc,angIdeg,FTS = readTransFile( FTSfile+"_transmission.txt", path="/o/plambeck/PolarBear/ARcoatings/FTS/" )

  # resample heterodyne data
    fout = open("interp.dat","w")
    fGHzHet = numpy.arange(fstart,fstop+fstep,fstep)
    hetInterp = interp1d( fGHz1, numpy.power(trans1,2.))
    pwrHet = hetInterp( fGHzHet )
    fout.write("# --- het data --- #\n")
    for f,p in zip(fGHzHet,pwrHet) :
      fout.write( "%10.4f %10.4f\n" % (f,p) )

  # resample FTS data over slightly wider freq range
    fGHzFTS = numpy.arange( fstart-30.*fstep,fstop+31.*fstep, fstep)
    FTSinterp = interp1d( fGHz2, trans2 )
    pwrFTS = FTSinterp( fGHzFTS )
    fout.write("# --- FTS data --- #\n")
    for f,p in zip(fGHzFTS, pwrFTS) :
      fout.write( "%10.4f, %10.4f\n" % (f,p) )
    fout.close()

    ccor = numpy.correlate( pwrFTS-numpy.average(pwrFTS), pwrHet-numpy.average(pwrHet), mode='valid')
    foffset = (numpy.argmax(ccor) - len(ccor)/2 )*fstep
       # so if signals are correlated with zero delay, argmax = 5 (counting from 0, not 1),
       # and len(ccor)/2 = 5 using integer arithmetic, so result should be 0
    # print foffset, ccor
    return foffset

# take out linear slope from FTS data to match heterodyne; also shift freq scale
def shiftFTS( FTSname, het3mmfile, het1mmfile, offsetGHz=-.4) :
    fGHz1,trans1,dphs_meas,unc,angIdeg,FTS = readTransFile( het3mmfile+".avg.dat", path="/o/plambeck/PolarBear/OpticsBench/" )
    fGHz2,trans2,dphs_meas,unc,angIdeg,FTS = readTransFile( het1mmfile+".avg.dat", path="/o/plambeck/PolarBear/OpticsBench/" )
    fGHz3,trans3,dphs_meas,unc,angIdeg,FTS = readTransFile( FTSname+"_transmission.txt", path="/o/plambeck/PolarBear/ARcoatings/FTS/" )
    f1,h1 = calcTavg( 90., 115.,fGHz1,numpy.power(trans1,2.))    
    f2,h2 = calcTavg( 205., 230.,fGHz2,numpy.power(trans2,2.))    
    f3,F3 = calcTavg( 90., 115., fGHz3, trans3 )
    f4,F4 = calcTavg( 205., 230., fGHz3, trans3 )
    fout = open( FTSname+"_shifted.txt", "w" )
    fout.write("# output from ob2.shiftFTS\n")
    fout.write("# %s with tilt removed, shifted by %.2f GHz\n" % (FTSname+"_transmission.txt", offsetGHz))
    fout.write("# FTS\n")
    fout.write("# angIdeg = 0.\n")
    for n in range(0, len(fGHz3)) :
      if fGHz3[n] >= 110. and fGHz3[n] <= 240. :
        Fav = F3 + (fGHz3[n]-f3)*(F4-F3)/(f4-f3)
        hav = h1 + (fGHz3[n]-f1)*(h2-h1)/(f2-f1)
        fout.write("%10.3f  %10.4f  %10.4f\n" % (fGHz3[n]+offsetGHz, trans3[n]-Fav+hav, unc[n]) )
    fout.close()

# computes mean freq and power for data between fstart and fstop
def calcTavg( fstart, fstop, fGHzArr, pwrArr ) : 
   npts = 0
   favg = 0.
   pavg = 0.
   for f,p in zip(fGHzArr,pwrArr) :
     if f >= fstart and f<=fstop :
       npts = npts + 1
       favg = favg + f
       pavg = pavg + p
   return (favg/npts),(pavg/npts)
   
# make fitFile to model a coating with an increasing density profile,
#   to mimic a Fillite coating; then use fuzzy2 to model
# want nr = (1 - exp(-x/xscale)) * nr / (1-exp(-t/xscale))
# t = total thickness of coating in inches
# nr is refractive index at the boundary with underlying alumina puck
# tscale is the 1/e scale factor for the refractive index
# dt = step size in inches
# tanD = tanDelta (assumed equal for all steps)
#
def profile( t, nr, tscale, dt=.0005, tanD=.01 ) :
    x = numpy.arange(dt,t+.0001,dt)
    ncoeff = (nr-1.)/(1.-math.exp(-t/tscale))
    print ncoeff
    nr = 1.+ncoeff * (1. - numpy.exp( -1.*(x-dt/2.)/tscale))
       # this is nr at the MIDPOINT of the layer 
    fout1 = open("profile.dat","w")
    fout1.write(" %8.5f %8.5f %8.3f\n" % (-.01, 1, 0.))
    fout1.write(" %8.5f %8.5f %8.3f\n" % (0., 1, 0.))
    fout2 = open("profile.fitFile","w")
    fout2.write("title = profile\n")
    fout2.write("totin = %.5f,%.5f\n" % (.25+t-.001,.25+t+.001))
    for n in range(0,len(x)) :
      fout1.write(" %8.5f %8.5f %8.3f\n" % (x[n],nr[n],tanD))
      fout2.write("\n  tin = %.5f\n" % dt)
      fout2.write("  nr = %.4f\n" % nr[n] )
      fout2.write("  tanDelta = %.3f\n" % tanD )

  # now include the alumina puck
    fout2.write("\n  tin = .25\n")
    fout2.write("  nr = 3.116\n")
    fout2.write("  tanDelta = 0.0004\n")
    fout1.write(" %8.5f %8.5f %8.3f\n" % (x[-1],3.116,.0004))
    fout1.write(" %8.5f %8.5f %8.3f\n" % (x[-1]+.25,3.116,.0004))
    fout1.write(" %8.5f %8.5f %8.3f\n" % (x[-1]+.25,1,0))
    fout1.write(" %8.5f %8.5f %8.3f\n" % (x[-1]+.26,1,0))
    fout1.close()
    fout2.close()

# ==== UNFINISHED!! ====
'''
# csvOut reads pickle file written by nrefracFit, returns csv lines giving:
#   line 1:  puckLabel,picklefile,dataFile(s),chisq,{layer1 params}
#   line 2:  ,,,,{layer2 params}
#   line N:  ,,,,{layerN params} 
#   where layer params are {tsrchmin,tmin,tbest,tmax,tsrchmax,nrsrchmin,nrmin,nrbest,nrmax,nrsrchmax,tanDsrchmin,etc}
def csvOut{ pickleFile }
    try :
      fin = open( pickleFile, "rb" )
      [fitParams,tcmList,nrList,tanDeltaList,chisq,dateString] = pickle.load( fin )
    except :
      fin = open( pickleFile, "rb" )
      [fitParams,tcmList,nrList,tanDeltaList,chisq] = pickle.load( fin )
      dateString = " "
    fin.close() 
    nlayers = len( tcmList[0] )

    print "\n10 best fits minimizing overall residuals:"
    nbest = printBest( chisq[-1], tcmList, nrList, tanDeltaList ) 

  # create savet,saven,savetanD lists, containing all configurations with chisq< 2.* minchisq
    savet = []
    savenr = []
    savetanD = []
    for n in range(0,len(chisq[-1])) :
      if chisq[-1,n] < 2.*chisq[-1,nbest]:
        savet.append( tcmList[n] )
        savenr.append( nrList[n] )
        savetanD.append( tanDeltaList[n] )

  fitParams = { "datafileList": datafileList, \
                "title" : title, \
                "tcmRange" : totcmRange, \
                "tcmList" : tcmList, \
                "nrList" : nrList, \
                "tanDeltaList" : tanDeltaList,
                "deltaTinList" : deltaTinList }
  # write the output string
    outString1 = "%s,%s,%s,%s" % (pickleFile[:pickleFile.index('.')],pickleFile,fitParams['datafileList'],chisq[-1])


# summaryPlot will make bar graph plot showing results for a 2 layer puck
def summaryPlot( infile ) :
    pyplot.ioff()
    pp = PdfPages( "summary.pdf")
    fig =  pyplot.figure( figsize=(11,8) )

    ymin = .1
    yrange = .8
    ax1 = fig.add_axes( [.05,ymin,.05,yrange] )   # puck labels
    ax2 = fig.add_axes( [.10,ymin,.20,yrange] )   # eps1
    ax3 = fig.add_axes( [.30,ymin,.20,yrange] )   # t1
    ax4 = fig.add_axes( [.50,ymin,.20,yrange] )   # eps2
    ax5 = fig.add_axes( [.70,ymin,.20,yrange] )   # t2
    
'''
